/* Phase 1: device-side BWA FM-index (Occ) primitives.
 *
 * Faithful CUDA mirrors of bwt_occ4 / bwt_2occ4 / bwt_occ / bwt_match_exact
 * from bwt.c, operating directly on the existing on-disk BWT layout:
 *   - bwt->bwt is a uint32 array of 64-byte buckets (one cache line each):
 *       [ 4 x uint64 Occ checkpoints ][ 8 x uint32 packed 2-bit BWT (128 bases) ]
 *   - OCC_INTERVAL = 128 (OCC_INTV_SHIFT = 7), so this code assumes that layout
 *     exactly as bwt.h does.
 * The suffix array is intentionally NOT needed here: bwa aln emits SA *intervals*,
 * and SA->coordinate happens later in samse on the host.
 *
 * cnt_table[256] and L2[5] live in constant memory (see fm_device.cu definitions).
 */
#ifndef FM_DEVICE_CUH
#define FM_DEVICE_CUH

#include <stdint.h>

#define D_OCC_INTV_MASK 0x7fULL   /* OCC_INTERVAL(128) - 1 */

/* constant-memory tables; defined once in the translation unit that sets
 * FM_DEVICE_DEFINE_CONST before including this header. */
#ifdef FM_DEVICE_DEFINE_CONST
__constant__ uint32_t c_cnt_table[256];
__constant__ uint64_t c_L2[5];
#else
extern __constant__ uint32_t c_cnt_table[256];
extern __constant__ uint64_t c_L2[5];
#endif

/* lightweight handle passed to kernels (scalars + device BWT pointer) */
struct fmidx_dev {
	const uint32_t *bwt;   /* device copy of bwt->bwt */
	uint64_t primary;      /* S^{-1}(0) */
	uint64_t seq_len;      /* number of symbols in the BWT */
};

/* __occ_aux4: count A/C/G/T in a 32-bit word (16 bases) via the packed cnt_table.
 * Returns 4 byte-packed counts (A in bits 0-7, C in 8-15, G in 16-23, T in 24-31). */
__device__ __forceinline__ uint32_t d_occ_aux4(uint32_t b)
{
	return c_cnt_table[b & 0xff] + c_cnt_table[(b >> 8) & 0xff]
	     + c_cnt_table[(b >> 16) & 0xff] + c_cnt_table[b >> 24];
}

/* mirror of bwt_occ4(): cumulative A/C/G/T counts in BWT[0..k]. */
__device__ __forceinline__ void d_bwt_occ4(const uint32_t *__restrict__ bwt,
                                            uint64_t primary, uint64_t k,
                                            uint64_t cnt[4])
{
	if (k == (uint64_t)-1) { cnt[0] = cnt[1] = cnt[2] = cnt[3] = 0; return; }
	k -= (k >= primary);                         /* $ is not in bwt */
	const uint32_t *p = bwt + ((k >> 7) << 4);   /* bucket base (bwt_occ_intv) */
	/* 4 x uint64 checkpoint counts at the bucket head */
	const uint64_t *pc = (const uint64_t *)p;
	cnt[0] = __ldg(pc + 0); cnt[1] = __ldg(pc + 1);
	cnt[2] = __ldg(pc + 2); cnt[3] = __ldg(pc + 3);
	p += 8;                                       /* skip the 4 uint64 (= 8 uint32) */
	uint64_t x = 0;
	const uint32_t *end = p + ((k >> 4) - ((k & ~D_OCC_INTV_MASK) >> 4));
	for (; p < end; ++p) x += d_occ_aux4(__ldg(p));
	uint32_t tmp = __ldg(p) & ~((1U << ((~k & 15) << 1)) - 1);
	x += (uint64_t)d_occ_aux4(tmp) - (~k & 15);
	cnt[0] += x & 0xff; cnt[1] += (x >> 8) & 0xff;
	cnt[2] += (x >> 16) & 0xff; cnt[3] += x >> 24;
}

/* mirror of bwt_2occ4(): Occ4 at both ends of an SA interval, with the shared-bucket
 * fast path (one bucket load + a single scan) when k and l land in the same 128-base
 * interval -- the common case as SA intervals shrink during search. */
__device__ __forceinline__ void d_bwt_2occ4(const uint32_t *__restrict__ bwt,
                                             uint64_t primary, uint64_t k, uint64_t l,
                                             uint64_t cntk[4], uint64_t cntl[4])
{
	uint64_t _k = k - (k >= primary), _l = l - (l >= primary);
	if ((_k >> 7) != (_l >> 7) || k == (uint64_t)-1 || l == (uint64_t)-1) {
		d_bwt_occ4(bwt, primary, k, cntk);
		d_bwt_occ4(bwt, primary, l, cntl);
		return;
	}
	k -= (k >= primary); l -= (l >= primary);
	const uint32_t *p = bwt + ((k >> 7) << 4);
	const uint64_t *pc = (const uint64_t *)p;
	cntk[0] = __ldg(pc + 0); cntk[1] = __ldg(pc + 1);
	cntk[2] = __ldg(pc + 2); cntk[3] = __ldg(pc + 3);
	p += 8;
	const uint32_t *endk = p + ((k >> 4) - ((k & ~D_OCC_INTV_MASK) >> 4));
	const uint32_t *endl = p + ((l >> 4) - ((l & ~D_OCC_INTV_MASK) >> 4));
	uint64_t x = 0, y;
	for (; p < endk; ++p) x += d_occ_aux4(__ldg(p));
	y = x;
	{ uint32_t tmp = __ldg(p) & ~((1U << ((~k & 15) << 1)) - 1); x += (uint64_t)d_occ_aux4(tmp) - (~k & 15); }
	for (; p < endl; ++p) y += d_occ_aux4(__ldg(p));
	{ uint32_t tmp = __ldg(p) & ~((1U << ((~l & 15) << 1)) - 1); y += (uint64_t)d_occ_aux4(tmp) - (~l & 15); }
	cntl[0] = cntk[0]; cntl[1] = cntk[1]; cntl[2] = cntk[2]; cntl[3] = cntk[3];
	cntk[0] += x & 0xff; cntk[1] += (x >> 8) & 0xff; cntk[2] += (x >> 16) & 0xff; cntk[3] += x >> 24;
	cntl[0] += y & 0xff; cntl[1] += (y >> 8) & 0xff; cntl[2] += (y >> 16) & 0xff; cntl[3] += y >> 24;
}

/* mirror of bwt_occ() for a single character c. occ4(seq_len) already yields the
 * total count of c (== L2[c+1]-L2[c]), so no special seq_len branch is needed. */
__device__ __forceinline__ uint64_t d_bwt_occ1(const uint32_t *__restrict__ bwt,
                                                uint64_t primary, uint64_t k, int c)
{
	if (k == (uint64_t)-1) return 0;
	uint64_t cnt[4];
	d_bwt_occ4(bwt, primary, k, cnt);
	return cnt[c];
}

/* mirror of bwt_match_exact(): backward search of str[0..len-1].
 * Returns SA interval size (l-k+1) and the interval in *sa_b/*sa_e, 0 on no match. */
__device__ __forceinline__ int d_bwt_match_exact(const fmidx_dev fm,
                                                  const uint8_t *str, int len,
                                                  uint64_t *sa_b, uint64_t *sa_e)
{
	uint64_t k = 0, l = fm.seq_len, ok, ol;
	for (int i = len - 1; i >= 0; --i) {
		uint8_t c = str[i];
		if (c > 3) return 0;                 /* N -> no match */
		ok = d_bwt_occ1(fm.bwt, fm.primary, k - 1, c);
		ol = d_bwt_occ1(fm.bwt, fm.primary, l,     c);
		k = c_L2[c] + ok + 1;
		l = c_L2[c] + ol;
		if (k > l) break;
	}
	if (k > l) return 0;
	if (sa_b) *sa_b = k;
	if (sa_e) *sa_e = l;
	return (int)(l - k + 1);
}

/* mirror of bwt_match_exact_alt(): backward search of str[0..len-1] starting from a
 * given SA interval (*k0,*l0). Returns interval size or 0; updates *k0,*l0 on match. */
__device__ __forceinline__ int d_bwt_match_exact_alt(const fmidx_dev fm,
                                                      const uint8_t *str, int len,
                                                      uint64_t *k0, uint64_t *l0)
{
	uint64_t k = *k0, l = *l0;
	for (int i = len - 1; i >= 0; --i) {
		uint8_t c = str[i];
		if (c > 3) return 0;
		uint64_t ck[4], cl[4];
		d_bwt_2occ4(fm.bwt, fm.primary, k - 1, l, ck, cl); /* shared-bucket fast path */
		k = c_L2[c] + ck[c] + 1;
		l = c_L2[c] + cl[c];
		if (k > l) return 0;
	}
	*k0 = k; *l0 = l;
	return (int)(l - k + 1);
}

/* mirror of int_log2() in bwtgap.c */
__device__ __forceinline__ int d_int_log2(uint32_t v)
{
	int c = 0;
	if (v & 0xffff0000u) { v >>= 16; c |= 16; }
	if (v & 0xff00) { v >>= 8; c |= 8; }
	if (v & 0xf0)   { v >>= 4; c |= 4; }
	if (v & 0xc)    { v >>= 2; c |= 2; }
	if (v & 0x2) c |= 1;
	return c;
}

#endif /* FM_DEVICE_CUH */
