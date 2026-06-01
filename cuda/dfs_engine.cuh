/* Shared warp-cooperative two-level-stack DFS engine (the production engine).
 * Include AFTER fm_device.cuh. Used by both cuda/dfstest.cu (validation) and
 * cuda/aln_gpu.cu (the streaming bwa-aln tool). See cuda/PROGRESS.md.
 */
#ifndef DFS_ENGINE_CUH
#define DFS_ENGINE_CUH

#include <stdint.h>

#ifndef DFS_STATE_DEFS
#define DFS_STATE_DEFS
#define STATE_M 0
#define STATE_I 1
#define STATE_D 2
#endif

/* per-read parameters that vary (other opts are kernel scalars, identical for all reads) */
struct ReadParam { uint32_t seq_off, w_off; int len, max_diff; };

/* DFS node packed into a u32: i(0-8) | n_mm(9-14) | n_gapo(15-18) | n_gape(19-22) | state(23-24) */
__device__ __forceinline__ uint32_t pack_node(int i, int mm, int go, int ge, int st)
{ return (uint32_t)i | (mm<<9) | (go<<15) | (ge<<19) | (st<<23); }

/* Two-level stack: per-warp shared top-window (CAP_SM) + per-warp GLOBAL backing.
 * Detects has_hit (existence) for a read; bushy reads spill their deep frontier to global and
 * stay on the GPU; only the 2M budget (or frontier > CAP_GL) flags to the CPU. has_hit-only ->
 * traversal order/shape is irrelevant -> bit-exact when paired with exact CPU reconcile. */
__device__ inline int d_dfs_has_hit_warp2(const fmidx_dev fm, const uint8_t *seq, int len,
                                   const uint64_t *w_w, const int *w_bid, int max_diff,
                                   int max_gapo, int max_gape, int mode, int indel_end_skip,
                                   int max_del_occ, uint64_t *sk, uint64_t *sl, uint32_t *sn, int CAP_SM,
                                   uint64_t *gk, uint64_t *gl, uint32_t *gn, int CAP_GL,
                                   unsigned long long budget, int *flagged, unsigned long long *nn_out)
{
	const unsigned FULL = 0xffffffffu; const int CHUNK = 128;
	int lane = threadIdx.x & 31;
	int nN = 0;
	for (int j = 0; j < len; ++j) if (seq[j] > 3) ++nN;
	if (nN > max_diff) { if (lane==0) *nn_out = 0; return 0; }

	if (lane == 0) { sk[0]=0; sl[0]=fm.seq_len; sn[0]=pack_node(len,0,0,0,STATE_M); }
	int sp = 1, gsp = 0;
	unsigned long long nn = 0;
	__syncwarp();

	for (;;) {
		if (sp == 0) {
			if (gsp == 0) { if (lane==0) *nn_out = nn; return 0; }
			int c = gsp < CHUNK ? gsp : CHUNK;
			for (int j = lane; j < c; j += 32) { sk[j]=gk[gsp-c+j]; sl[j]=gl[gsp-c+j]; sn[j]=gn[gsp-c+j]; }
			gsp -= c; sp = c; __syncwarp();
		}
		int room = CAP_SM - sp;
		if (room < 9) {
			if (gsp + CHUNK > CAP_GL) { if (lane==0){ *flagged=1; *nn_out=nn; } return 1; }
			for (int j = lane; j < CHUNK; j += 32) { gk[gsp+j]=sk[j]; gl[gsp+j]=sl[j]; gn[gsp+j]=sn[j]; }
			__syncwarp();
			for (int j = lane; j < sp - CHUNK; j += 32) { sk[j]=sk[j+CHUNK]; sl[j]=sl[j+CHUNK]; sn[j]=sn[j+CHUNK]; }
			gsp += CHUNK; sp -= CHUNK; __syncwarp();
			room = CAP_SM - sp;
		}
		int n_active = sp < 32 ? sp : 32;
		int r9 = room / 9; if (n_active > r9) n_active = r9; if (n_active < 1) n_active = 1;
		bool active = lane < n_active;
		uint64_t k=0, l=0; int i=0, e_mm=0, e_go=0, e_ge=0, e_st=0;
		if (active) { int idx = sp-1-lane; k=sk[idx]; l=sl[idx]; uint32_t nd=sn[idx];
			i=nd&0x1ff; e_mm=(nd>>9)&0x3f; e_go=(nd>>15)&0xf; e_ge=(nd>>19)&0xf; e_st=(nd>>23)&0x3; }
		sp -= n_active; nn += n_active;
		if (nn > budget) { if (lane==0){ *flagged=1; *nn_out=nn; } return 1; }

		uint64_t cck[9], ccl[9]; uint32_t ccn[9]; int nc = 0; bool hit = false;
		if (active) {
			int m = max_diff - (e_mm + e_go);
			if (mode & BWA_MODE_GAPE) m -= e_ge;
			bool prune = (m < 0) || (i > 0 && m < w_bid[i-1]);
			if (!prune) {
				if (i == 0) hit = true;
				else if (m == 0 && (e_st==STATE_M || (mode&BWA_MODE_GAPE) || e_ge==max_gape)) {
					uint64_t kk=k, ll=l;
					if (d_bwt_match_exact_alt(fm, seq, i, &kk, &ll)) hit = true; else prune = true;
				}
				if (!hit && !prune) {
					int xi = i - 1;
					uint64_t cntk[4], cntl[4];
					d_bwt_2occ4(fm.bwt, fm.primary, k - 1, l, cntk, cntl);
					uint64_t occ = l - k + 1;
					int allow_diff = 1, allow_M = 1;
					if (xi > 0) {
						if (w_bid[xi-1] > m-1) allow_diff = 0;
						else if (w_bid[xi-1]==m-1 && w_bid[xi]==m-1 && w_w[xi-1]==w_w[xi]) allow_M = 0;
					}
					int tmp = (mode & BWA_MODE_LOGGAP) ? d_int_log2(e_ge+e_go)/2+1 : e_go+e_ge;
					#define ADD2(_i,_k,_l,_mm,_go,_ge,_st) do { \
						cck[nc]=(_k); ccl[nc]=(_l); ccn[nc]=pack_node((_i),(_mm),(_go),(_ge),(_st)); ++nc; } while(0)
					if (allow_diff && xi >= indel_end_skip+tmp && len-xi >= indel_end_skip+tmp) {
						if (e_st == STATE_M) {
							if (e_go < max_gapo) {
								ADD2(xi, k, l, e_mm, e_go+1, e_ge, STATE_I);
								for (int j=0;j<4;++j){ uint64_t dk=c_L2[j]+cntk[j]+1, dl=c_L2[j]+cntl[j];
									if (dk<=dl) ADD2(xi+1, dk, dl, e_mm, e_go+1, e_ge, STATE_D); }
							}
						} else if (e_st == STATE_I) {
							if (e_ge < max_gape) ADD2(xi, k, l, e_mm, e_go, e_ge+1, STATE_I);
						} else {
							if (e_ge < max_gape && (e_ge+e_go < max_diff || occ < (uint64_t)max_del_occ))
								for (int j=0;j<4;++j){ uint64_t dk=c_L2[j]+cntk[j]+1, dl=c_L2[j]+cntl[j];
									if (dk<=dl) ADD2(xi+1, dk, dl, e_mm, e_go, e_ge+1, STATE_D); }
						}
					}
					if (allow_diff && allow_M) {
						for (int j=1;j<=4;++j){ int c=(seq[xi]+j)&3; int is_mm=(j!=4||seq[xi]>3);
							uint64_t mk=c_L2[c]+cntk[c]+1, ml=c_L2[c]+cntl[c];
							if (mk<=ml) ADD2(xi, mk, ml, e_mm+is_mm, e_go, e_ge, STATE_M); }
					} else if (seq[xi] < 4) {
						int c=seq[xi]&3; uint64_t mk=c_L2[c]+cntk[c]+1, ml=c_L2[c]+cntl[c];
						if (mk<=ml) ADD2(xi, mk, ml, e_mm, e_go, e_ge, STATE_M);
					}
					#undef ADD2
				}
			}
		}
		if (__any_sync(FULL, hit)) { if (lane==0) *nn_out = nn; return 1; }

		int incl = nc;
		for (int d=1; d<32; d<<=1) { int y = __shfl_up_sync(FULL, incl, d); if (lane >= d) incl += y; }
		int total = __shfl_sync(FULL, incl, 31);
		int myoff = incl - nc;
		for (int j=0;j<nc;++j){ int d = sp + myoff + j; sk[d]=cck[j]; sl[d]=ccl[j]; sn[d]=ccn[j]; }
		sp += total;
		__syncwarp();
	}
}

/* persistent work-pool: warps pull reads via an atomic counter until the queue drains */
__global__ void k_dfs_warp2(fmidx_dev fm, const uint8_t *seq, const uint64_t *w_w, const int *w_bid,
                            const ReadParam *rp, int nreads, int max_gapo, int max_gape, int mode,
                            int indel_end_skip, int max_del_occ, int CAP_SM, int CAP_GL,
                            uint64_t *Gk, uint64_t *Gl, uint32_t *Gn, uint8_t *has_hit,
                            int *workctr, unsigned long long *npop, unsigned long long budget,
                            int *nflag, int wpb)
{
	extern __shared__ unsigned char smem[];
	int warp_in_block = threadIdx.x >> 5, lane = threadIdx.x & 31;
	uint64_t *sk = (uint64_t*)smem + (size_t)warp_in_block * CAP_SM;
	uint64_t *sl = (uint64_t*)smem + (size_t)wpb * CAP_SM + (size_t)warp_in_block * CAP_SM;
	uint32_t *sn = (uint32_t*)((uint64_t*)smem + (size_t)wpb * CAP_SM * 2) + (size_t)warp_in_block * CAP_SM;
	int gslot = blockIdx.x * wpb + warp_in_block;
	uint64_t *gk = Gk + (size_t)gslot * CAP_GL;
	uint64_t *gl = Gl + (size_t)gslot * CAP_GL;
	uint32_t *gn = Gn + (size_t)gslot * CAP_GL;
	for (;;) {
		int r;
		if (lane == 0) r = atomicAdd(workctr, 1);
		r = __shfl_sync(0xffffffffu, r, 0);
		if (r >= nreads) break;
		ReadParam p = rp[r];
		int flagged = 0; unsigned long long nn = 0;
		int hh = d_dfs_has_hit_warp2(fm, seq + p.seq_off, p.len, w_w + p.w_off, w_bid + p.w_off,
			p.max_diff, max_gapo, max_gape, mode, indel_end_skip, max_del_occ, sk, sl, sn, CAP_SM,
			gk, gl, gn, CAP_GL, budget, &flagged, &nn);
		if (lane == 0) { has_hit[r] = (uint8_t)hh; npop[r] = nn; if (flagged) atomicAdd(nflag, 1); }
	}
}

#endif /* DFS_ENGINE_CUH */
