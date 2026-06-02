/* GPU bwa-aln / fused alnse (the `bwa gpualn` subcommand) -- GPU BWA-backtrack
   port for ancient DNA.
   Copyright (C) 2026  teepean  <https://github.com/teepean>
   Derived from bwa (Heng Li; Broad Institute / Dana-Farber / Genome Research Ltd.).

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>. */

/* bwa-aln-gpu: full-file GPU bwa-aln (BWA-backtrack) producing a bit-exact .sai.
 *
 * Streams the FASTQ in bwa's native 0x40000-read chunks (so chunking matches the CPU driver
 * exactly). Per chunk: MT host preprocessing (bwt_cal_width + complement + per-read max_diff,
 * identical to bwa_cal_sa_reg_gap) -> warp2 GPU has_hit -> MT CPU reconcile of the ~0.015%
 * flagged/hit reads via the exact bwt_match_gap -> write the chunk's records in read order.
 * Output is byte-identical to `bwa aln -l 1024 -n 0.01 -o 2`.
 *
 * Build: make bwa-aln-gpu
 * Run:   ./bwa-aln-gpu [-l N -n F -o N -t N -f out.sai] <ref.fa> <in.fq>   (default opts as above)
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <thread>
#include <chrono>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <unistd.h>
#include <cuda_runtime.h>

extern "C" {
#include "bwt.h"
#include "bwtaln.h"
#include "bwtgap.h"
#include "bntseq.h"
#include "bwase.h"
#include "bwa.h"
#include "utils.h"
int bwt_cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width);
/* in bwase.c but not in bwase.h */
void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi);
void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2);
extern char *bwa_pg;   /* @PG line printed by bwa_print_sam_hdr if set */
}

#define FM_DEVICE_DEFINE_CONST
#include "fm_device.cuh"
#include "dfs_engine.cuh"

#define CK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { \
	fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e_)); exit(1); } } while (0)

static double now_s() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

/* one chunk's host data, handed from the GPU producer thread to the CPU finisher thread */
struct Chunk {
	bwa_seq_t *seqs; int n_seqs;
	std::vector<ReadParam> rp;
	std::vector<uint8_t> seq_flat; std::vector<uint64_t> w_flat; std::vector<int> bid_flat;
	std::vector<uint8_t> has_hit;
	gap_opt_t base; int stack_maxdiff, max_len;
};

/* C-callable entry: usable as a bwa subcommand (`bwa gpualn ...`) or standalone (ALN_GPU_MAIN).
 * argv[0] is the program/subcommand name; options are parsed from argv[1..] (bwa convention). */
extern "C" int bwa_alnse_gpu(int argc, char **argv)
{
	gap_opt_t *opt = gap_init_opt();
	opt->seed_len = 1024; opt->fnr = 0.01; opt->max_diff = -1; opt->max_gapo = 2; opt->n_threads = 16;
	const char *out_fn = NULL; int wpb = 4, CAP_SM = 512, CAP_GL = 16384;
	int sam_mode = 0, n_occ = 3; char *rg_line = NULL;   /* -S: fused alnse (SAM out); -r: RG */
	int c;
	while ((c = getopt(argc, argv, "l:n:o:t:f:Sr:")) >= 0) {
		if (c=='l') opt->seed_len = atoi(optarg);
		else if (c=='n') { if (strstr(optarg,".")) { opt->fnr=atof(optarg); opt->max_diff=-1; } else { opt->max_diff=atoi(optarg); opt->fnr=-1; } }
		else if (c=='o') opt->max_gapo = atoi(optarg);
		else if (c=='t') opt->n_threads = atoi(optarg);
		else if (c=='f') out_fn = optarg;
		else if (c=='S') sam_mode = 1;
		else if (c=='r') { if ((rg_line = bwa_set_rg(optarg)) == 0) return 1; }
	}
	if (optind + 2 > argc) { fprintf(stderr, "usage: %s [-l N -n F -o N -t N -f out.sai] <ref.fa> <in.fq>\n", argv[0]); return 1; }
	const char *prefix = argv[optind], *fq = argv[optind+1];
	if (getenv("DFS_WARP_CAP")) CAP_SM = atoi(getenv("DFS_WARP_CAP"));
	if (getenv("DFS_WARP_GCAP")) CAP_GL = atoi(getenv("DFS_WARP_GCAP"));
	unsigned long long budget = getenv("DFS_BUDGET") ? strtoull(getenv("DFS_BUDGET"),NULL,10) : 2000000ULL;
	int nT = opt->n_threads > 0 ? opt->n_threads : 1;

	char bwt_fn[4096]; snprintf(bwt_fn, sizeof bwt_fn, "%s.bwt", prefix);
	fprintf(stderr, "[aln-gpu] loading %s\n", bwt_fn);
	bwt_t *bwt = bwt_restore_bwt(bwt_fn);
	if (!bwt) { fprintf(stderr, "failed to load bwt\n"); return 1; }

	uint32_t *d_bwt = NULL;
	CK(cudaMalloc(&d_bwt, bwt->bwt_size * sizeof(uint32_t)));
	CK(cudaMemcpy(d_bwt, bwt->bwt, bwt->bwt_size * sizeof(uint32_t), cudaMemcpyHostToDevice));
	CK(cudaMemcpyToSymbol(c_cnt_table, bwt->cnt_table, sizeof(uint32_t)*256));
	CK(cudaMemcpyToSymbol(c_L2, bwt->L2, sizeof(uint64_t)*5));
	fmidx_dev fm{ d_bwt, bwt->primary, bwt->seq_len };

	/* occupancy + one-time global backing allocation */
	int numSM = 0; CK(cudaDeviceGetAttribute(&numSM, cudaDevAttrMultiProcessorCount, 0));
	int bdim = wpb * 32; size_t shbytes = (size_t)wpb * CAP_SM * 20;
	CK(cudaFuncSetAttribute(k_dfs_warp2, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)shbytes));
	int mb = 0; CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&mb, k_dfs_warp2, bdim, shbytes));
	int nblocks = mb > 0 ? mb * numSM : numSM; size_t nwarps = (size_t)nblocks * wpb;
	uint64_t *Gk, *Gl; uint32_t *Gn;
	CK(cudaMalloc(&Gk, nwarps*CAP_GL*8)); CK(cudaMalloc(&Gl, nwarps*CAP_GL*8)); CK(cudaMalloc(&Gn, nwarps*CAP_GL*4));
	fprintf(stderr, "[aln-gpu] %d blk/SM x %d warps, CAP_SM=%d CAP_GL=%d, backing %.0f MB; %d CPU threads\n",
		mb, wpb, CAP_SM, CAP_GL, nwarps*CAP_GL*20.0/1e6, nT);

	/* output: .sai (default) or fused alnse SAM (-S) */
	FILE *out = NULL; bntseq_t *bns = NULL; ubyte_t *pacseq = NULL;
	if (sam_mode) {
		bwase_initialize();
		bns = bns_restore(prefix);
		srand48(bns->seed);                       /* exact samse RNG seeding (repeat-hit selection) */
		char sa_fn[4096]; snprintf(sa_fn, sizeof sa_fn, "%s.sa", prefix);
		bwt_restore_sa(sa_fn, bwt);                /* SA for bwa_sa2pos; bwt kept (not destroyed) */
		pacseq = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
		err_rewind(bns->fp_pac); err_fread_noeof(pacseq, 1, bns->l_pac/4+1, bns->fp_pac);
		if (!bwa_pg) { /* @PG provenance (standalone only; as a bwa subcommand, main.c sets bwa_pg) */
			char pg[8192]; int o = snprintf(pg, sizeof pg, "@PG\tID:bwa-aln-gpu\tPN:bwa-aln-gpu\tVN:gpu\tCL:");
			for (int i=0;i<argc && o<(int)sizeof pg-2;i++) o += snprintf(pg+o, sizeof pg-o, "%s%s", i?" ":"", argv[i]);
			bwa_pg = strdup(pg);
		}
		bwa_print_sam_hdr(bns, rg_line);           /* @HD/@SQ/@RG/@PG header to stdout */
		fprintf(stderr, "[aln-gpu] fused alnse (SAM) mode; SA+pac loaded\n");
	} else {
		out = out_fn ? fopen(out_fn, "wb") : stdout;
		if (!out) { fprintf(stderr, "cannot open %s\n", out_fn); return 1; }
		err_fwrite(SAI_MAGIC, 1, 4, out);
		err_fwrite(opt, sizeof(gap_opt_t), 1, out);
	}

	bwa_seqio_t *ks = bwa_seq_open(fq);

	/* reusable device buffers, grown on demand */
	uint8_t *d_seq=NULL; uint64_t *d_ww=NULL; int *d_wbid=NULL; ReadParam *d_rp=NULL;
	uint8_t *d_hit=NULL; unsigned long long *d_npop=NULL; int *d_wc=NULL, *d_nflag=NULL;
	size_t cap_seq=0, cap_w=0, cap_n=0;

	long long tot=0, tot_flag=0; double t0=now_s();

	/* ---- CPU/GPU overlap (#5): the main thread owns the GPU (read -> preprocess -> upload ->
	 * kernel -> download has_hit, serial on one stream so the shared global backing is never
	 * raced). Each finished chunk is handed to ONE in-order finisher thread that does the CPU
	 * reconcile + output, running concurrently with the next chunk's GPU work. Single ordered
	 * consumer => drand48/output order preserved => bit-exact. (For multi-GPU later: replicate the
	 * GPU producer stage per device and round-robin chunks; the finisher stays a single ordered
	 * consumer to keep output/RNG order.) ---- */
	std::mutex qmu; std::condition_variable q_ready, q_free; std::queue<Chunk*> Q; bool done=false;
	const size_t QCAP = 2;
	auto finisher = [&](){
		for (;;) {
			Chunk *c;
			{ std::unique_lock<std::mutex> lk(qmu); q_ready.wait(lk, [&]{ return !Q.empty() || done; });
			  if (Q.empty()) break; c = Q.front(); Q.pop(); }
			q_free.notify_one();
			int nseq = c->n_seqs;
			std::vector<int> idx; for (int i=0;i<nseq;i++) if (c->has_hit[i]) idx.push_back(i);
			tot_flag += idx.size();
			std::vector<int> n_aln(nseq,0); std::vector<bwt_aln1_t*> aln(nseq,NULL);
			std::vector<std::thread> ths;
			for (int t=0;t<nT;t++) ths.emplace_back([&,t](){
				gap_stack_t *st = gap_init_stack(c->stack_maxdiff, c->base.max_gapo, c->base.max_gape, &c->base);
				std::vector<bwt_width_t> w(c->max_len+1);
				for (size_t x=t;x<idx.size();x+=nT){
					int i=idx[x], len=c->rp[i].len;
					for (int j=0;j<=len;j++){ w[j].w=c->w_flat[c->rp[i].w_off+j]; w[j].bid=c->bid_flat[c->rp[i].w_off+j]; }
					gap_opt_t lo=c->base; lo.max_diff=c->rp[i].max_diff; lo.seed_len = opt->seed_len<len?opt->seed_len:0x7fffffff;
					int na=0; aln[i]=bwt_match_gap(bwt, len, c->seq_flat.data()+c->rp[i].seq_off, w.data(), (bwt_width_t*)0, &lo, &na, st);
					n_aln[i]=na;
				}
				gap_destroy_stack(st);
			});
			for (auto&th:ths) th.join();
			if (!sam_mode) {
				for (int i=0;i<nseq;i++){ err_fwrite(&n_aln[i],4,1,out); if (n_aln[i]) err_fwrite(aln[i],sizeof(bwt_aln1_t),n_aln[i],out); free(aln[i]); }
			} else {
				for (int i=0;i<nseq;i++){ bwa_aln2seq_core(n_aln[i], aln[i], &c->seqs[i], 1, n_occ); free(aln[i]); }
				ths.clear();
				for (int t=0;t<nT;t++) ths.emplace_back([&,t](){
					for (int i=t;i<nseq;i+=nT){ bwa_seq_t *p=&c->seqs[i];
						bwa_cal_pac_pos_core(bns, bwt, p, opt->max_diff, opt->fnr);
						int strand, nm2=0;
						for (int j=0;j<p->n_multi;j++){ bwt_multi1_t *q=p->multi+j;
							q->pos = bwa_sa2pos(bns, bwt, q->pos, p->len+q->ref_shift, &strand); q->strand=strand;
							if (q->pos != p->pos && q->pos != (bwtint_t)-1) p->multi[nm2++]=*q; }
						p->n_multi=nm2; }
				});
				for (auto&th:ths) th.join();
				bwa_refine_gapped(bns, nseq, c->seqs, pacseq);
				for (int i=0;i<nseq;i++) bwa_print_sam1(bns, &c->seqs[i], 0, opt->mode, opt->max_top2);
			}
			tot += nseq;
			bwa_free_read_seq(nseq, c->seqs);
			delete c;
			fprintf(stderr, "\r[aln-gpu] %lld reads done (%.0f reads/s)   ", tot, tot/(now_s()-t0));
		}
	};
	std::thread cons(finisher);

	int n_seqs; bwa_seq_t *seqs;
	while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual)) != 0) {
		Chunk *c = new Chunk();
		c->seqs = seqs; c->n_seqs = n_seqs; c->base = *opt; c->max_len = 0;
		for (int i=0;i<n_seqs;i++) if (seqs[i].len > c->max_len) c->max_len = seqs[i].len;
		if (opt->fnr > 0.0) c->base.max_diff = bwa_cal_maxdiff(c->max_len, BWA_AVG_ERR, opt->fnr);
		if (c->base.max_diff < c->base.max_gapo) c->base.max_gapo = c->base.max_diff;
		c->stack_maxdiff = c->base.max_diff;
		c->rp.resize(n_seqs);
		size_t so=0, wo=0;
		for (int i=0;i<n_seqs;i++){ c->rp[i].seq_off=so; c->rp[i].w_off=wo; c->rp[i].len=seqs[i].len; so+=seqs[i].len; wo+=seqs[i].len+1; }
		c->seq_flat.resize(so); c->w_flat.resize(wo); c->bid_flat.resize(wo); c->has_hit.resize(n_seqs);

		std::vector<std::thread> ths;   /* MT preprocess (width pre-complement + complement into flat) */
		for (int t=0;t<nT;t++) ths.emplace_back([&,t](){
			std::vector<bwt_width_t> w(c->max_len+1);
			for (int i=t;i<n_seqs;i+=nT){ bwa_seq_t *p=&seqs[i];
				memset(w.data(), 0, (p->len+1)*sizeof(bwt_width_t));
				bwt_cal_width(bwt, p->len, p->seq, w.data());
				c->rp[i].max_diff = (opt->fnr>0.0)? bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr) : c->base.max_diff;
				for (int j=0;j<p->len;j++) c->seq_flat[c->rp[i].seq_off+j] = p->seq[j]>3?4:3-p->seq[j];
				for (int j=0;j<=p->len;j++){ c->w_flat[c->rp[i].w_off+j]=w[j].w; c->bid_flat[c->rp[i].w_off+j]=w[j].bid; }
			}
		});
		for (auto&th:ths) th.join();

		/* GPU stage (serial on main thread) */
		if (so>cap_seq){ if(d_seq)cudaFree(d_seq); CK(cudaMalloc(&d_seq, so)); cap_seq=so; }
		if (wo>cap_w){ if(d_ww)cudaFree(d_ww); if(d_wbid)cudaFree(d_wbid); CK(cudaMalloc(&d_ww,wo*8)); CK(cudaMalloc(&d_wbid,wo*4)); cap_w=wo; }
		if ((size_t)n_seqs>cap_n){ if(d_rp)cudaFree(d_rp); if(d_hit)cudaFree(d_hit); if(d_npop)cudaFree(d_npop);
			CK(cudaMalloc(&d_rp,n_seqs*sizeof(ReadParam))); CK(cudaMalloc(&d_hit,n_seqs)); CK(cudaMalloc(&d_npop,n_seqs*8)); cap_n=n_seqs; }
		if (!d_wc){ CK(cudaMalloc(&d_wc,4)); CK(cudaMalloc(&d_nflag,4)); }
		CK(cudaMemcpy(d_seq, c->seq_flat.data(), so, cudaMemcpyHostToDevice));
		CK(cudaMemcpy(d_ww, c->w_flat.data(), wo*8, cudaMemcpyHostToDevice));
		CK(cudaMemcpy(d_wbid, c->bid_flat.data(), wo*4, cudaMemcpyHostToDevice));
		CK(cudaMemcpy(d_rp, c->rp.data(), n_seqs*sizeof(ReadParam), cudaMemcpyHostToDevice));
		CK(cudaMemset(d_wc,0,4)); CK(cudaMemset(d_nflag,0,4));
		k_dfs_warp2<<<nblocks, bdim, shbytes>>>(fm, d_seq, d_ww, d_wbid, d_rp, n_seqs,
			c->base.max_gapo, c->base.max_gape, c->base.mode, c->base.indel_end_skip, c->base.max_del_occ,
			CAP_SM, CAP_GL, Gk, Gl, Gn, d_hit, d_wc, d_npop, budget, d_nflag, wpb);
		CK(cudaDeviceSynchronize()); CK(cudaGetLastError());
		CK(cudaMemcpy(c->has_hit.data(), d_hit, n_seqs, cudaMemcpyDeviceToHost));

		{ std::unique_lock<std::mutex> lk(qmu); q_free.wait(lk, [&]{ return Q.size() < QCAP; }); Q.push(c); }
		q_ready.notify_one();
	}
	{ std::lock_guard<std::mutex> lk(qmu); done = true; } q_ready.notify_one();
	cons.join();
	double total = now_s()-t0;
	fprintf(stderr, "\n[aln-gpu] DONE: %lld reads in %.1f s = %.0f reads/s; flagged->CPU %lld (%.3f%%)\n",
		tot, total, tot/total, tot_flag, 100.0*tot_flag/tot);

	if (out_fn) fclose(out);
	bwa_seq_close(ks);
	bwt_destroy(bwt);
	return 0;
}

#ifdef ALN_GPU_MAIN
int main(int argc, char **argv) { return bwa_alnse_gpu(argc, argv); }
#endif
