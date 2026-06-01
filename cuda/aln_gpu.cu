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

int main(int argc, char **argv)
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
		{ /* @PG provenance line */
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

	long long tot=0, tot_flag=0; double t_pre=0,t_gpu=0,t_rec=0,t_io=0; double t0=now_s();
	int n_seqs;
	bwa_seq_t *seqs;
	while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual)) != 0) {
		double a = now_s();
		/* chunk-level opt (matches bwa_cal_sa_reg_gap stack sizing) */
		gap_opt_t base = *opt; int max_len = 0;
		for (int i=0;i<n_seqs;i++) if (seqs[i].len > max_len) max_len = seqs[i].len;
		if (opt->fnr > 0.0) base.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
		if (base.max_diff < base.max_gapo) base.max_gapo = base.max_diff;
		int stack_maxdiff = base.max_diff;

		/* offsets (serial, cheap) */
		std::vector<ReadParam> rp(n_seqs);
		size_t so=0, wo=0;
		for (int i=0;i<n_seqs;i++){ rp[i].seq_off=so; rp[i].w_off=wo; rp[i].len=seqs[i].len; so+=seqs[i].len; wo+=seqs[i].len+1; }
		std::vector<uint8_t> seq_flat(so); std::vector<uint64_t> w_flat(wo); std::vector<int> bid_flat(wo);

		/* MT preprocess: width (pre-complement) + per-read max_diff + complement */
		std::vector<std::thread> ths;
		for (int t=0;t<nT;t++) ths.emplace_back([&,t](){
			std::vector<bwt_width_t> w(max_len+1);
			for (int i=t;i<n_seqs;i+=nT){
				bwa_seq_t *p=&seqs[i];
				memset(w.data(), 0, (p->len+1)*sizeof(bwt_width_t));
				bwt_cal_width(bwt, p->len, p->seq, w.data());
				rp[i].max_diff = (opt->fnr>0.0)? bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr) : base.max_diff;
				/* complement into the flat GPU buffer ONLY; leave p->seq in read-seq orientation
				 * (samse needs it; reconcile uses seq_flat) */
				for (int j=0;j<p->len;j++) seq_flat[rp[i].seq_off+j] = p->seq[j]>3?4:3-p->seq[j];
				for (int j=0;j<=p->len;j++){ w_flat[rp[i].w_off+j]=w[j].w; bid_flat[rp[i].w_off+j]=w[j].bid; }
			}
		});
		for (auto&th:ths) th.join();
		t_pre += now_s()-a; a=now_s();

		/* (re)alloc device buffers */
		if (so>cap_seq){ if(d_seq)cudaFree(d_seq); CK(cudaMalloc(&d_seq, so)); cap_seq=so; }
		if (wo>cap_w){ if(d_ww)cudaFree(d_ww); if(d_wbid)cudaFree(d_wbid); CK(cudaMalloc(&d_ww,wo*8)); CK(cudaMalloc(&d_wbid,wo*4)); cap_w=wo; }
		if ((size_t)n_seqs>cap_n){ if(d_rp)cudaFree(d_rp); if(d_hit)cudaFree(d_hit); if(d_npop)cudaFree(d_npop);
			CK(cudaMalloc(&d_rp,n_seqs*sizeof(ReadParam))); CK(cudaMalloc(&d_hit,n_seqs)); CK(cudaMalloc(&d_npop,n_seqs*8)); cap_n=n_seqs; }
		if (!d_wc){ CK(cudaMalloc(&d_wc,4)); CK(cudaMalloc(&d_nflag,4)); }
		CK(cudaMemcpy(d_seq, seq_flat.data(), so, cudaMemcpyHostToDevice));
		CK(cudaMemcpy(d_ww, w_flat.data(), wo*8, cudaMemcpyHostToDevice));
		CK(cudaMemcpy(d_wbid, bid_flat.data(), wo*4, cudaMemcpyHostToDevice));
		CK(cudaMemcpy(d_rp, rp.data(), n_seqs*sizeof(ReadParam), cudaMemcpyHostToDevice));
		CK(cudaMemset(d_wc,0,4)); CK(cudaMemset(d_nflag,0,4));

		k_dfs_warp2<<<nblocks, bdim, shbytes>>>(fm, d_seq, d_ww, d_wbid, d_rp, n_seqs,
			base.max_gapo, base.max_gape, base.mode, base.indel_end_skip, base.max_del_occ,
			CAP_SM, CAP_GL, Gk, Gl, Gn, d_hit, d_wc, d_npop, budget, d_nflag, wpb);
		CK(cudaDeviceSynchronize()); CK(cudaGetLastError());
		std::vector<uint8_t> has_hit(n_seqs);
		CK(cudaMemcpy(has_hit.data(), d_hit, n_seqs, cudaMemcpyDeviceToHost));
		t_gpu += now_s()-a; a=now_s();

		/* MT reconcile flagged/hit reads via exact CPU search */
		std::vector<int> idx; for (int i=0;i<n_seqs;i++) if (has_hit[i]) idx.push_back(i);
		tot_flag += idx.size();
		std::vector<int> n_aln(n_seqs,0); std::vector<bwt_aln1_t*> aln(n_seqs,NULL);
		ths.clear();
		for (int t=0;t<nT;t++) ths.emplace_back([&,t](){
			gap_stack_t *st = gap_init_stack(stack_maxdiff, base.max_gapo, base.max_gape, &base);
			std::vector<bwt_width_t> w(max_len+1);
			for (size_t x=t;x<idx.size();x+=nT){
				int i=idx[x], len=rp[i].len;
				for (int j=0;j<=len;j++){ w[j].w=w_flat[rp[i].w_off+j]; w[j].bid=bid_flat[rp[i].w_off+j]; }
				gap_opt_t lo=base; lo.max_diff=rp[i].max_diff; lo.seed_len = opt->seed_len<len?opt->seed_len:0x7fffffff;
				int na=0; aln[i]=bwt_match_gap(bwt, len, seq_flat.data()+rp[i].seq_off, w.data(), (bwt_width_t*)0, &lo, &na, st);
				n_aln[i]=na;
			}
			gap_destroy_stack(st);
		});
		for (auto&th:ths) th.join();
		t_rec += now_s()-a; a=now_s();

		if (!sam_mode) {
			/* write .sai records in read order (matches bwa_aln_core) */
			for (int i=0;i<n_seqs;i++){ err_fwrite(&n_aln[i],4,1,out); if (n_aln[i]) err_fwrite(aln[i],sizeof(bwt_aln1_t),n_aln[i],out); free(aln[i]); }
		} else {
			/* fused samse: serial aln2seq (drand48 order) -> MT SA-lookup (RNG-free) ->
			 * serial refine_gapped + print (bit-exact, ordered output). */
			for (int i=0;i<n_seqs;i++) { bwa_aln2seq_core(n_aln[i], aln[i], &seqs[i], 1, n_occ); free(aln[i]); }
			ths.clear();
			for (int t=0;t<nT;t++) ths.emplace_back([&,t](){
				for (int i=t;i<n_seqs;i+=nT){
					bwa_seq_t *p=&seqs[i];
					bwa_cal_pac_pos_core(bns, bwt, p, opt->max_diff, opt->fnr);
					int strand, nm2=0;
					for (int j=0;j<p->n_multi;j++){ bwt_multi1_t *q=p->multi+j;
						q->pos = bwa_sa2pos(bns, bwt, q->pos, p->len+q->ref_shift, &strand); q->strand=strand;
						if (q->pos != p->pos && q->pos != (bwtint_t)-1) p->multi[nm2++]=*q; }
					p->n_multi=nm2;
				}
			});
			for (auto&th:ths) th.join();
			bwa_refine_gapped(bns, n_seqs, seqs, pacseq);
			for (int i=0;i<n_seqs;i++) bwa_print_sam1(bns, &seqs[i], 0, opt->mode, opt->max_top2);
		}
		t_io += now_s()-a;
		tot += n_seqs;
		bwa_free_read_seq(n_seqs, seqs);
		fprintf(stderr, "\r[aln-gpu] %lld reads processed (%.0f reads/s overall)   ", tot, tot/(now_s()-t0));
	}
	double total = now_s()-t0;
	fprintf(stderr, "\n[aln-gpu] DONE: %lld reads in %.1f s = %.0f reads/s; flagged->CPU %lld (%.3f%%)\n",
		tot, total, tot/total, tot_flag, 100.0*tot_flag/tot);
	fprintf(stderr, "[aln-gpu] time: preprocess %.1fs  gpu %.1fs  reconcile %.1fs  io %.1fs\n", t_pre, t_gpu, t_rec, t_io);

	if (out_fn) fclose(out);
	bwa_seq_close(ks);
	bwt_destroy(bwt);
	return 0;
}
