/* Phase 2 step 1: bit-exact device DFS (hit-detection) for bwa aln.
 *
 * The device DFS reproduces bwt_match_gap's bounded search tree and reports, per read,
 * whether ANY hit exists within the initial diff bound (has_hit). Rationale (see
 * ../CUDA_PORT_PLAN.md sec 6a): for reads that never reach a hit -- and when the
 * max_entries cap is never hit -- best-first and DFS enumerate the identical node set,
 * so existence-of-hit on the GPU == (n_aln>0) on the CPU. The ~0.5% reads with hits are
 * reconciled by the exact CPU bwt_match_gap to yield a bit-exact .sai.
 *
 * Host side reuses bwa's real read I/O + preprocessing (bwt_cal_width, complement, per-read
 * max_diff) so the GPU and CPU see identical (seq, width, opt).
 *
 * Build: make dfstest
 * Run:   ./dfstest <ref.fa> <reads.fq> [golden.sai]
 *        replicates `bwa aln -l 1024 -n 0.01 -o 2 -t 16`.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cuda_runtime.h>

extern "C" {
#include "bwt.h"
#include "bwtaln.h"
#include "bwtgap.h"
int bwt_cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width); /* not in a header */
}

#define FM_DEVICE_DEFINE_CONST
#include "fm_device.cuh"

#define CK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { \
	fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e_)); \
	exit(1); } } while (0)

#define STATE_M 0
#define STATE_I 1
#define STATE_D 2

/* per-read parameters that vary (the rest are kernel scalars, identical for all reads) */
struct ReadParam { uint32_t seq_off, w_off; int len, max_diff; };

/* DFS node packed into a u32: i(0-8) | n_mm(9-14) | n_gapo(15-18) | n_gape(19-22) | state(23-24) */
__device__ __forceinline__ uint32_t pack_node(int i, int mm, int go, int ge, int st)
{ return (uint32_t)i | (mm<<9) | (go<<15) | (ge<<19) | (st<<23); }

/* ---- device DFS over a COALESCED SoA stack: returns 1 if any hit exists within bound.
 * Stack entry d for this thread lives at [(size_t)d*P + slot] in each of sk/sl/sn, so a
 * warp's 32 threads touch 32 contiguous slots per push/pop (coalesced). ---- */
__device__ int d_dfs_has_hit(const fmidx_dev fm, const uint8_t *seq, int len,
                             const uint64_t *w_w, const int *w_bid, int max_diff,
                             int max_gapo, int max_gape, int mode, int indel_end_skip,
                             int max_del_occ, uint64_t *sk, uint64_t *sl, uint32_t *sn,
                             int P, int slot, int cap, int *overflow,
                             unsigned long long budget, int *nflag, unsigned long long &nn)
{
	nn = 0;
	int nN = 0;
	for (int j = 0; j < len; ++j) if (seq[j] > 3) ++nN;
	if (nN > max_diff) return 0;

	#define SIDX(d) ((size_t)(d) * P + slot)
	int sp = 0;
	/* current node held in registers: only SIBLINGS are pushed to the global stack, and
	 * the stack is read only on backtrack. This removes the push->immediately-pop global
	 * round-trip that dominated runtime (stack working set was blowing L2). */
	uint64_t k = 0, l = fm.seq_len;
	int i = len, e_mm = 0, e_go = 0, e_ge = 0, e_st = STATE_M;
	bool have_cur = true;

	for (;;) {
		if (!have_cur) {                       /* backtrack: pop a sibling */
			if (sp == 0) break;
			--sp;
			k = sk[SIDX(sp)]; l = sl[SIDX(sp)];
			uint32_t nd = sn[SIDX(sp)];
			i = nd & 0x1ff; e_mm = (nd>>9)&0x3f; e_go = (nd>>15)&0xf; e_ge = (nd>>19)&0xf; e_st = (nd>>23)&0x3;
		}
		have_cur = false;                      /* consume current; re-set only if we descend */
		++nn;
		if (nn > budget) { atomicAdd(nflag, 1); return 1; }

		int m = max_diff - (e_mm + e_go);
		if (mode & BWA_MODE_GAPE) m -= e_ge;
		if (m < 0) continue;
		if (i > 0 && m < w_bid[i-1]) continue;

		int hit_found = 0;
		if (i == 0) hit_found = 1;
		else if (m == 0 && (e_st == STATE_M || (mode & BWA_MODE_GAPE) || e_ge == max_gape)) {
			uint64_t kk = k, ll = l;
			if (d_bwt_match_exact_alt(fm, seq, i, &kk, &ll)) hit_found = 1;
			else continue;
		}
		if (hit_found) return 1; /* existence is all we need */

		int xi = i - 1;
		uint64_t cntk[4], cntl[4];
		d_bwt_2occ4(fm.bwt, fm.primary, k - 1, l, cntk, cntl);
		uint64_t occ = l - k + 1;
		int allow_diff = 1, allow_M = 1;
		if (xi > 0) {
			if (w_bid[xi-1] > m - 1) allow_diff = 0;
			else if (w_bid[xi-1] == m-1 && w_bid[xi] == m-1 && w_w[xi-1] == w_w[xi]) allow_M = 0;
		}
		int tmp = (mode & BWA_MODE_LOGGAP) ? d_int_log2(e_ge + e_go)/2 + 1 : e_go + e_ge;

		/* one child stays in registers (the last EMITted, descended into); EMIT pushes
		 * the previously-staged child so exactly nchild-1 hit the global stack. */
		bool have_next = false;
		uint64_t nk_=0, nl_=0; int ni_=0, nmm_=0, ngo_=0, nge_=0, nst_=0;
		#define EMIT(_i,_k,_l,_mm,_go,_ge,_st) do { \
			if (have_next) { if (sp >= cap) { *overflow = 1; return 1; } \
				sk[SIDX(sp)]=nk_; sl[SIDX(sp)]=nl_; sn[SIDX(sp)]=pack_node(ni_,nmm_,ngo_,nge_,nst_); ++sp; } \
			nk_=(_k); nl_=(_l); ni_=(_i); nmm_=(_mm); ngo_=(_go); nge_=(_ge); nst_=(_st); have_next=true; \
		} while (0)

		if (allow_diff && xi >= indel_end_skip + tmp && len - xi >= indel_end_skip + tmp) {
			if (e_st == STATE_M) {
				if (e_go < max_gapo) {
					EMIT(xi, k, l, e_mm, e_go + 1, e_ge, STATE_I);          /* insertion */
					for (int j = 0; j < 4; ++j) {                            /* deletion */
						uint64_t dk = c_L2[j] + cntk[j] + 1, dl = c_L2[j] + cntl[j];
						if (dk <= dl) EMIT(xi + 1, dk, dl, e_mm, e_go + 1, e_ge, STATE_D);
					}
				}
			} else if (e_st == STATE_I) {
				if (e_ge < max_gape) EMIT(xi, k, l, e_mm, e_go, e_ge + 1, STATE_I);
			} else { /* STATE_D */
				if (e_ge < max_gape && (e_ge + e_go < max_diff || occ < (uint64_t)max_del_occ)) {
					for (int j = 0; j < 4; ++j) {
						uint64_t dk = c_L2[j] + cntk[j] + 1, dl = c_L2[j] + cntl[j];
						if (dk <= dl) EMIT(xi + 1, dk, dl, e_mm, e_go, e_ge + 1, STATE_D);
					}
				}
			}
		}
		if (allow_diff && allow_M) {
			for (int j = 1; j <= 4; ++j) {
				int c = (seq[xi] + j) & 3;
				int is_mm = (j != 4 || seq[xi] > 3);
				uint64_t mk = c_L2[c] + cntk[c] + 1, ml = c_L2[c] + cntl[c];
				if (mk <= ml) EMIT(xi, mk, ml, e_mm + is_mm, e_go, e_ge, STATE_M);
			}
		} else if (seq[xi] < 4) {
			int c = seq[xi] & 3;
			uint64_t mk = c_L2[c] + cntk[c] + 1, ml = c_L2[c] + cntl[c];
			if (mk <= ml) EMIT(xi, mk, ml, e_mm, e_go, e_ge, STATE_M);
		}
		#undef EMIT

		if (have_next) {                       /* descend into last-staged child (registers) */
			k = nk_; l = nl_; i = ni_; e_mm = nmm_; e_go = ngo_; e_ge = nge_; e_st = nst_;
			have_cur = true;
		}
	}
	#undef SIDX
	return 0;
}

/* persistent work-pool kernel: each thread atomically pulls the next read until the queue
 * drains, so a thread finishing a tiny tree immediately grabs more work (load balancing). */
__global__ void k_dfs(fmidx_dev fm, const uint8_t *seq, const uint64_t *w_w, const int *w_bid,
                      const ReadParam *rp, int nreads, int max_gapo, int max_gape, int mode,
                      int indel_end_skip, int max_del_occ, uint64_t *sk, uint64_t *sl, uint32_t *sn,
                      int cap, uint8_t *has_hit, int *overflow, int *workctr, unsigned long long *npop,
                      unsigned long long budget, int *nflag)
{
	int slot = blockIdx.x * blockDim.x + threadIdx.x;
	int P = gridDim.x * blockDim.x;
	for (;;) {
		int r = atomicAdd(workctr, 1);
		if (r >= nreads) break;
		ReadParam p = rp[r];
		unsigned long long nn = 0;
		has_hit[r] = (uint8_t)d_dfs_has_hit(fm, seq + p.seq_off, p.len,
			w_w + p.w_off, w_bid + p.w_off, p.max_diff,
			max_gapo, max_gape, mode, indel_end_skip, max_del_occ,
			sk, sl, sn, P, slot, cap, overflow, budget, nflag, nn);
		npop[r] = nn;
	}
}

/* ---- write a .sai exactly like bwa_aln_core ---- */
static void write_sai(const char *fn, const gap_opt_t *opt, int n, int *n_aln, bwt_aln1_t **aln)
{
	FILE *f = fopen(fn, "wb");
	fwrite(SAI_MAGIC, 1, 4, f);
	fwrite(opt, sizeof(gap_opt_t), 1, f);
	for (int i = 0; i < n; ++i) {
		fwrite(&n_aln[i], 4, 1, f);
		if (n_aln[i]) fwrite(aln[i], sizeof(bwt_aln1_t), n_aln[i], f);
	}
	fclose(f);
}

static int md5cmp(const char *a, const char *b)
{
	char cmd[2048]; snprintf(cmd, sizeof cmd, "md5sum '%s' '%s'", a, b);
	fprintf(stderr, "  "); fflush(stderr); return system(cmd);
}

int main(int argc, char **argv)
{
	if (argc < 3) { fprintf(stderr, "usage: %s <ref.fa> <reads.fq> [golden.sai]\n", argv[0]); return 1; }
	const char *prefix = argv[1], *reads_fn = argv[2], *golden = argc > 3 ? argv[3] : NULL;

	/* opt for: bwa aln -l 1024 -n 0.01 -o 2 -t 16 */
	gap_opt_t *opt = gap_init_opt();
	opt->seed_len = 1024; opt->fnr = 0.01; opt->max_diff = -1; opt->max_gapo = 2; opt->n_threads = 16;

	char bwt_fn[4096]; snprintf(bwt_fn, sizeof bwt_fn, "%s.bwt", prefix);
	fprintf(stderr, "[load] %s\n", bwt_fn);
	bwt_t *bwt = bwt_restore_bwt(bwt_fn);
	if (!bwt) { fprintf(stderr, "failed to load bwt\n"); return 1; }

	/* upload index + constants */
	uint32_t *d_bwt = NULL;
	CK(cudaMalloc(&d_bwt, bwt->bwt_size * sizeof(uint32_t)));
	CK(cudaMemcpy(d_bwt, bwt->bwt, bwt->bwt_size * sizeof(uint32_t), cudaMemcpyHostToDevice));
	CK(cudaMemcpyToSymbol(c_cnt_table, bwt->cnt_table, sizeof(uint32_t) * 256));
	CK(cudaMemcpyToSymbol(c_L2, bwt->L2, sizeof(uint64_t) * 5));
	fmidx_dev fm{ d_bwt, bwt->primary, bwt->seq_len };

	/* read ALL reads in one batch (sub10k < 0x40000), exactly as the driver would */
	bwa_seqio_t *ks = bwa_seq_open(reads_fn);
	int n_seqs = 0;
	bwa_seq_t *seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual);
	bwa_seq_close(ks);
	fprintf(stderr, "[reads] %d\n", n_seqs);

	/* replicate bwa_cal_sa_reg_gap preprocessing + CPU reference */
	gap_opt_t local_opt = *opt;
	int max_len = 0; for (int i = 0; i < n_seqs; ++i) if (seqs[i].len > max_len) max_len = seqs[i].len;
	if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
	if (local_opt.max_diff < local_opt.max_gapo) local_opt.max_gapo = local_opt.max_diff;
	gap_stack_t *stack = gap_init_stack(local_opt.max_diff, local_opt.max_gapo, local_opt.max_gape, &local_opt);
	gap_opt_t base = local_opt; /* stable max_gapo/mode/...; max_diff & seed_len overridden per read */
	/* full per-read CPU reference (checks 1&3) is slow; default on only for small sets */
	bool full_cpu = getenv("DFS_FULLCPU") || n_seqs <= 20000;

	std::vector<uint8_t> seq_flat; std::vector<uint64_t> w_flat; std::vector<int> bid_flat;
	std::vector<ReadParam> rp(n_seqs);
	std::vector<int> n_aln_cpu(n_seqs);
	std::vector<bwt_aln1_t*> aln_cpu(n_seqs, NULL);
	bwt_width_t *w = NULL; int max_w = 0;

	for (int i = 0; i < n_seqs; ++i) {
		bwa_seq_t *p = &seqs[i];
		if (max_w < p->len) { max_w = p->len; w = (bwt_width_t*)realloc(w, (max_w + 1) * sizeof(bwt_width_t)); }
		memset(w, 0, (p->len + 1) * sizeof(bwt_width_t));
		bwt_cal_width(bwt, p->len, p->seq, w);
		if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
		local_opt.seed_len = opt->seed_len < p->len ? opt->seed_len : 0x7fffffff;
		for (int j = 0; j < p->len; ++j) p->seq[j] = p->seq[j] > 3 ? 4 : 3 - p->seq[j]; /* complement */

		/* record inputs for the GPU (identical seq + width + per-read max_diff) */
		rp[i].seq_off = seq_flat.size(); rp[i].w_off = w_flat.size();
		rp[i].len = p->len; rp[i].max_diff = local_opt.max_diff;
		for (int j = 0; j < p->len; ++j) seq_flat.push_back(p->seq[j]);
		for (int j = 0; j <= p->len; ++j) { w_flat.push_back(w[j].w); bid_flat.push_back(w[j].bid); }

		/* optional full CPU reference for checks 1 & 3 (len <= seed_len -> seed_w NULL) */
		if (full_cpu) {
			int na = 0;
			aln_cpu[i] = bwt_match_gap(bwt, p->len, p->seq, w, (bwt_width_t*)0, &local_opt, &na, stack);
			n_aln_cpu[i] = na;
		}
	}
	free(w);

	if (full_cpu) { /* check 1: CPU-only .sai == golden (validates the harness replicates the driver) */
		write_sai("/tmp/dfstest.cpu.sai", opt, n_seqs, n_aln_cpu.data(), aln_cpu.data());
		long n_hit_cpu = 0; for (int i = 0; i < n_seqs; ++i) if (n_aln_cpu[i] > 0) ++n_hit_cpu;
		fprintf(stderr, "[cpu] reads with hits: %ld (%.2f%%)\n", n_hit_cpu, 100.0 * n_hit_cpu / n_seqs);
		if (golden) { fprintf(stderr, "[check1] CPU-only .sai vs golden:\n"); md5cmp("/tmp/dfstest.cpu.sai", golden); }
	} else fprintf(stderr, "[cpu] full reference skipped (set DFS_FULLCPU=1 to enable checks 1&3)\n");

	/* ---- GPU DFS ---- */
	uint8_t *d_seq; uint64_t *d_ww; int *d_wbid; ReadParam *d_rp; uint8_t *d_hit; int *d_ovf;
	CK(cudaMalloc(&d_seq, seq_flat.size()));
	CK(cudaMalloc(&d_ww, w_flat.size() * 8));
	CK(cudaMalloc(&d_wbid, bid_flat.size() * 4));
	CK(cudaMalloc(&d_rp, n_seqs * sizeof(ReadParam)));
	CK(cudaMalloc(&d_hit, n_seqs));
	CK(cudaMalloc(&d_ovf, 4)); CK(cudaMemset(d_ovf, 0, 4));
	CK(cudaMemcpy(d_seq, seq_flat.data(), seq_flat.size(), cudaMemcpyHostToDevice));
	CK(cudaMemcpy(d_ww, w_flat.data(), w_flat.size() * 8, cudaMemcpyHostToDevice));
	CK(cudaMemcpy(d_wbid, bid_flat.data(), bid_flat.size() * 4, cudaMemcpyHostToDevice));
	CK(cudaMemcpy(d_rp, rp.data(), n_seqs * sizeof(ReadParam), cudaMemcpyHostToDevice));

	int blockDim = 128, cap = 1024;
	int numSM = 0, maxBlocks = 0;
	CK(cudaDeviceGetAttribute(&numSM, cudaDevAttrMultiProcessorCount, 0));
	CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&maxBlocks, k_dfs, blockDim, 0));
	int nblocks = maxBlocks * numSM, P = nblocks * blockDim;
	size_t poolsz = (size_t)P * cap * (8 + 8 + 4);
	fprintf(stderr, "[gpu] %d SM x %d blk/SM x %d thr = %d slots; stack cap %d -> coalesced pool %.0f MB\n",
		numSM, maxBlocks, blockDim, P, cap, poolsz / 1e6);
	uint64_t *d_sk, *d_sl; uint32_t *d_sn; int *d_wc;
	CK(cudaMalloc(&d_sk, (size_t)P * cap * 8));
	CK(cudaMalloc(&d_sl, (size_t)P * cap * 8));
	CK(cudaMalloc(&d_sn, (size_t)P * cap * 4));
	CK(cudaMalloc(&d_wc, 4)); CK(cudaMemset(d_wc, 0, 4));
	unsigned long long *d_npop; CK(cudaMalloc(&d_npop, (size_t)n_seqs * 8)); CK(cudaMemset(d_npop, 0, (size_t)n_seqs * 8));
	int *d_nflag; CK(cudaMalloc(&d_nflag, 4)); CK(cudaMemset(d_nflag, 0, 4));
	unsigned long long budget = getenv("DFS_BUDGET") ? strtoull(getenv("DFS_BUDGET"), NULL, 10) : 2000000ULL;
	fprintf(stderr, "[gpu] per-read work budget = %llu pops (overflow -> CPU reconcile)\n", budget);

	cudaEvent_t t0, t1; CK(cudaEventCreate(&t0)); CK(cudaEventCreate(&t1));
	CK(cudaEventRecord(t0));
	k_dfs<<<nblocks, blockDim>>>(fm, d_seq, d_ww, d_wbid, d_rp, n_seqs,
		local_opt.max_gapo, local_opt.max_gape, local_opt.mode, local_opt.indel_end_skip,
		local_opt.max_del_occ, d_sk, d_sl, d_sn, cap, d_hit, d_ovf, d_wc, d_npop, budget, d_nflag);
	CK(cudaEventRecord(t1)); CK(cudaEventSynchronize(t1)); CK(cudaGetLastError());
	float ms = 0; CK(cudaEventElapsedTime(&ms, t0, t1));
	{ int nflag = 0; CK(cudaMemcpy(&nflag, d_nflag, 4, cudaMemcpyDeviceToHost));
	  fprintf(stderr, "[gpu] budget-flagged reads (sent to CPU) = %d (%.3f%%)\n", nflag, 100.0*nflag/n_seqs); }

	{ std::vector<unsigned long long> np(n_seqs); CK(cudaMemcpy(np.data(), d_npop, (size_t)n_seqs*8, cudaMemcpyDeviceToHost));
	  unsigned long long tot=0,mx=0; int amx=0; for(int i=0;i<n_seqs;i++){ tot+=np[i]; if(np[i]>mx){mx=np[i];amx=i;} }
	  fprintf(stderr, "[gpu] node-pops total=%llu mean=%.0f max=%llu (read %d, len %d)\n",
	          tot, (double)tot/n_seqs, mx, amx, rp[amx].len);
	  /* top-5 heaviest reads */
	  std::vector<int> idx(n_seqs); for(int i=0;i<n_seqs;i++) idx[i]=i;
	  std::sort(idx.begin(), idx.end(), [&](int a,int b){return np[a]>np[b];});
	  fprintf(stderr, "[gpu] heaviest reads (pops): "); for(int t=0;t<5&&t<n_seqs;t++) fprintf(stderr,"%llu ", np[idx[t]]); fprintf(stderr,"\n"); }

	std::vector<uint8_t> has_hit(n_seqs); int ovf = 0;
	CK(cudaMemcpy(has_hit.data(), d_hit, n_seqs, cudaMemcpyDeviceToHost));
	CK(cudaMemcpy(&ovf, d_ovf, 4, cudaMemcpyDeviceToHost));
	fprintf(stderr, "[gpu] DFS time %.1f ms  (%.0f reads/s)  stack_overflow=%d\n",
		ms, n_seqs / (ms / 1e3), ovf);

	long gpu_hit = 0; for (int i = 0; i < n_seqs; ++i) if (has_hit[i]) ++gpu_hit;
	long fn = 0;
	if (full_cpu) { /* check 3: has_hit vs CPU n_aln>0 */
		long fp = 0;
		for (int i = 0; i < n_seqs; ++i) {
			int cpu = n_aln_cpu[i] > 0;
			if (cpu && !has_hit[i]) { if (fn < 10) fprintf(stderr, "  FALSE-NEG read %d (cpu n_aln=%d)\n", i, n_aln_cpu[i]); ++fn; }
			if (!cpu && has_hit[i]) ++fp;
		}
		fprintf(stderr, "[check3] gpu has_hit=%ld  false_neg=%ld (MUST be 0)  false_pos=%ld (harmless)\n",
			gpu_hit, fn, fp);
	}

	/* the real production hybrid: reconcile the GPU-flagged reads via the exact CPU search,
	 * reconstructing each read's pristine width from the stored flat arrays. */
	std::vector<int> n_aln_hyb(n_seqs, 0);
	std::vector<bwt_aln1_t*> aln_hyb(n_seqs, NULL);
	std::vector<bwt_width_t> wtmp(max_len + 1);
	clock_t rc0 = clock();
	for (int i = 0; i < n_seqs; ++i) {
		if (!has_hit[i]) continue;
		int len = rp[i].len;
		for (int j = 0; j <= len; ++j) { wtmp[j].w = w_flat[rp[i].w_off + j]; wtmp[j].bid = bid_flat[rp[i].w_off + j]; }
		gap_opt_t lo = base; lo.max_diff = rp[i].max_diff;
		lo.seed_len = opt->seed_len < len ? opt->seed_len : 0x7fffffff;
		int na = 0;
		aln_hyb[i] = bwt_match_gap(bwt, len, seq_flat.data() + rp[i].seq_off, wtmp.data(), (bwt_width_t*)0, &lo, &na, stack);
		n_aln_hyb[i] = na;
	}
	double rcs = (double)(clock() - rc0) / CLOCKS_PER_SEC;
	fprintf(stderr, "[hybrid] CPU reconcile of %ld flagged reads: %.2f s (1 thread; parallelizable)\n", gpu_hit, rcs);

	write_sai("/tmp/dfstest.hybrid.sai", opt, n_seqs, n_aln_hyb.data(), aln_hyb.data());
	if (golden) { fprintf(stderr, "[check2] hybrid .sai vs golden:\n"); md5cmp("/tmp/dfstest.hybrid.sai", golden); }

	fprintf(stderr, "[result] GPU %.0f reads/s; hybrid %s\n",
		n_seqs / (ms / 1e3), (golden ? "compared above (md5 must match)" : "(provide golden.sai to verify)"));
	return (full_cpu && fn) ? 2 : 0;
}
