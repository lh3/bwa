#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include "bwtaln.h"
#include "kvec.h"
#include "bntseq.h"
#include "utils.h"
#include "stdaln.h"
#include "bwatpx.h"

#include "khash.h"
KHASH_MAP_INIT_INT64(64, poslist_t)

#include "ksort.h"
KSORT_INIT_GENERIC(uint64_t)

#define MIN_HASH_WIDTH 1000

kh_64_t *g_hash[MAX_CPUS];

uint64_t n_tot[2], n_mapped[2];

extern int g_log_n[]; // in bwase.c

extern char bwaversionstr[];
extern char bwablddatestr[];

void bwase_initialize();
void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi);
void bwa_aln2seq(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s);
void bwa_refine_gapped(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, ubyte_t *_pacseq, bntseq_t *ntbns);
int  bwa_approx_mapQ(const bwa_seq_t *p, int mm);
void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2);
bntseq_t *bwa_open_nt(const char *prefix);
void bwa_print_sam_SQ(const bntseq_t *bns);
void bwa_print_sam_PG(void);
extern void bwa_read_seq2_tpx(bwa_seqio_t *ks1, bwa_seqio_t *ks2, int n_needed, int *n,
				int mode1, int mode2, int trim_qual1, int trim_qual2,
				bwa_seq_t **seq1, bwa_seq_t **seq2, pe_data_t *d[MAX_CPUS], 
				aln_buf_t **buf1, aln_buf_t **buf2, int nexti, FILE *fp_sa[2]);
extern void bwa_read_seq2_wait_tpx(void);
extern void bwa_se_tpx(int iidx, bwt_t *bwt[2], int n_seqs1, int n_seqs2, bwa_seq_t *seqs[2], const gap_opt_t *gopt);
extern int  bwa_pe_tpx(int iidx, bwt_t *bwt[2], int n_seqs1, int n_seqs2, bwa_seq_t *seqs[2], isize_info_t *ii,
			const pe_opt_t *opt, const gap_opt_t *gopt, pe_data_t *d[MAX_CPUS], aln_buf_t *buf[2]);
extern void bwa_sw_tpx(int iidx, const bntseq_t *bns, const ubyte_t *pacseq, int n_seqs1, int n_seqs2, bwa_seq_t *seqs[2], 
			const pe_opt_t *popt, const isize_info_t *ii);
extern void bwa_print2_tpx(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs[2], const gap_opt_t opt);
extern void bwa_print2_wait_tpx(void);

int num_sampe_threads = 1;
THR_BWA_SE_TPX thr_bwa_se_info[MAX_CPUS];
THR_BWA_PE_TPX thr_bwa_pe_info[MAX_CPUS];
THR_BWA_SW_TPX thr_bwa_sw_info[MAX_CPUS];

int adj_n_needed = 1;
int async_read_seq = 1;
int async_print_res = 1;
int use_soap2_format = 0;
int soap2_qual_adj = 33;

static clock_t read_aln_clocks = 0;

// -------------------

pe_opt_t *bwa_init_pe_opt()
{
	pe_opt_t *po;
	po = (pe_opt_t*)calloc(1, sizeof(pe_opt_t));
	po->max_isize = 500;
	po->force_isize = 0;
	po->max_occ = 100000;
	po->n_multi = 3;
	po->N_multi = 10;
	po->type = BWA_PET_STD;
	po->is_sw = 1;
	po->ap_prior = 1e-5;
	return po;
}

static inline uint64_t hash_64(uint64_t key)
{
	key += ~(key << 32);
	key ^= (key >> 22);
	key += ~(key << 13);
	key ^= (key >> 8);
	key += (key << 3);
	key ^= (key >> 15);
	key += ~(key << 27);
	key ^= (key >> 31);
	return key;
}
/*
static double ierfc(double x) // inverse erfc(); iphi(x) = M_SQRT2 *ierfc(2 * x);
{
	const double a = 0.140012;
	double b, c;
	b = log(x * (2 - x));
	c = 2./M_PI/a + b / 2.;
	return sqrt(sqrt(c * c - b / a) - c);
}
*/

// for normal distribution, this is about 3std
#define OUTLIER_BOUND 2.0

static int infer_isize(int n_seqs, bwa_seq_t *seqs[2], isize_info_t *ii, double ap_prior, int64_t L)
{
	uint64_t x, *isizes, n_ap = 0;
	int n, i, tot, p25, p75, p50, max_len = 1, tmp;
	double skewness = 0.0, kurtosis = 0.0, y;

	ii->avg = ii->std = -1.0;
	ii->low = ii->high = ii->high_bayesian = 0;
	isizes = (uint64_t*)calloc(n_seqs, 8);
	for (i = 0, tot = 0; i != n_seqs; ++i) {
		bwa_seq_t *p[2];
		p[0] = seqs[0] + i; p[1] = seqs[1] + i;
		if (p[0]->mapQ >= 20 && p[1]->mapQ >= 20) {
			x = (p[0]->pos < p[1]->pos)? p[1]->pos + p[1]->len - p[0]->pos : p[0]->pos + p[0]->len - p[1]->pos;
			if (x < 100000) isizes[tot++] = x;
		}
		if (p[0]->len > max_len) max_len = p[0]->len;
		if (p[1]->len > max_len) max_len = p[1]->len;
	}
	if (tot < 20) {
		fprintf(stderr, "[infer_isize] fail to infer insert size: too few good pairs\n");
		free(isizes);
		return -1;
	}
	ks_introsort(uint64_t, tot, isizes);
	p25 = isizes[(int)(tot*0.25 + 0.5)];
	p50 = isizes[(int)(tot*0.50 + 0.5)];
	p75 = isizes[(int)(tot*0.75 + 0.5)];
	tmp  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
	ii->low = tmp > max_len? tmp : max_len; // ii->low is unsigned
	ii->high = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);
	for (i = 0, x = n = 0; i < tot; ++i)
		if (isizes[i] >= ii->low && isizes[i] <= ii->high)
			++n, x += isizes[i];
	ii->avg = (double)x / n;
	for (i = 0; i < tot; ++i) {
		if (isizes[i] >= ii->low && isizes[i] <= ii->high) {
			double tmp = (isizes[i] - ii->avg) * (isizes[i] - ii->avg);
			ii->std += tmp;
			skewness += tmp * (isizes[i] - ii->avg);
			kurtosis += tmp * tmp;
		}
	}
	kurtosis = kurtosis/n / (ii->std / n * ii->std / n) - 3;
	ii->std = sqrt(ii->std / n); // it would be better as n-1, but n is usually very large
	skewness = skewness / n / (ii->std * ii->std * ii->std);
	for (y = 1.0; y < 10.0; y += 0.01)
		if (.5 * erfc(y / M_SQRT2) < ap_prior / L * (y * ii->std + ii->avg)) break;
	ii->high_bayesian = (bwtint_t)(y * ii->std + ii->avg + .499);
	for (i = 0; i < tot; ++i)
		if (isizes[i] > ii->high_bayesian) ++n_ap;
	ii->ap_prior = .01 * (n_ap + .01) / tot;
	if (ii->ap_prior < ap_prior) ii->ap_prior = ap_prior;
	free(isizes);
	fprintf(stderr, "[infer_isize] (25, 50, 75) percentile: (%d, %d, %d)\n", p25, p50, p75);
	if (isnan(ii->std) || p75 > 100000) {
		ii->low = ii->high = ii->high_bayesian = 0; ii->avg = ii->std = -1.0;
		fprintf(stderr, "[infer_isize] fail to infer insert size: weird pairing\n");
		return -1;
	}
	for (y = 1.0; y < 10.0; y += 0.01)
		if (.5 * erfc(y / M_SQRT2) < ap_prior / L * (y * ii->std + ii->avg)) break;
	ii->high_bayesian = (bwtint_t)(y * ii->std + ii->avg + .499);
	fprintf(stderr, "[infer_isize] low and high boundaries: %d and %d for estimating avg and std\n", ii->low, ii->high);
	fprintf(stderr, "[infer_isize] inferred external isize from %d pairs: %.3lf +/- %.3lf\n", n, ii->avg, ii->std);
	fprintf(stderr, "[infer_isize] skewness: %.3lf; kurtosis: %.3lf; ap_prior: %.2e\n", skewness, kurtosis, ii->ap_prior);
	fprintf(stderr, "[infer_isize] inferred maximum insert size: %d (%.2lf sigma)\n", ii->high_bayesian, y);
	return 0;
}

// -------------------

void thr_bwa_se_tpx(long idx)
{
  int iidx = (int)idx;

  bwa_se_tpx(iidx,
             thr_bwa_se_info[iidx].bwt,
             thr_bwa_se_info[iidx].start,
             thr_bwa_se_info[iidx].end,
             thr_bwa_se_info[iidx].seqs,
             thr_bwa_se_info[iidx].gopt);

  return;
}

// -------------------

void thr_bwa_pe_tpx(long idx)
{
  int iidx = (int)idx;

  thr_bwa_pe_info[iidx].cnt_chg = bwa_pe_tpx(iidx,
                                             thr_bwa_pe_info[iidx].bwt,
                                             thr_bwa_pe_info[iidx].start,
                                             thr_bwa_pe_info[iidx].end,
                                             thr_bwa_pe_info[iidx].seqs,
                                             thr_bwa_pe_info[iidx].ii,
                                             thr_bwa_pe_info[iidx].opt,
                                             thr_bwa_pe_info[iidx].gopt,
                                             thr_bwa_pe_info[iidx].d,
                                             thr_bwa_pe_info[iidx].buf);

  return;
}

// -------------------

int bwa_cal_pac_pos_pe(const char *prefix, bwt_t *const _bwt[2], int n_seqs, bwa_seq_t *seqs[2], FILE *fp_sa[2], isize_info_t *ii,
			const pe_opt_t *opt, const gap_opt_t *gopt, const isize_info_t *last_ii, pe_data_t *d[MAX_CPUS], aln_buf_t *buf[2])
{
	int i, j, cnt_chg = 0;
	char str[1024];
	bwt_t *bwt[2];
	// clock_t t;
#ifdef HAVE_PTHREAD
	int srtn = 0;
	long delta = 0L;
	pthread_t tid;
#endif // HAVE_PTHREAD

#if 0
	// ---------------

	fprintf(stderr, "[bwa_sai2sam_pe_core] read alignments... ");

	t = clock();

	// SE
	for (i = 0; i < n_seqs; ++i) {
		bwa_seq_t *p[2];
		for (j = 0; j < 2; ++j) {
			int n_aln;
			p[j] = seqs[j] + i;
			p[j]->n_multi = 0;
			p[j]->extra_flag |= SAM_FPD | (j == 0? SAM_FR1 : SAM_FR2);

			fread(&n_aln, 4, 1, fp_sa[j]);

			if (n_aln > kv_max(d[0]->aln[j]))
				kv_resize(bwt_aln1_t, d[0]->aln[j], n_aln);

			d[0]->aln[j].n = n_aln;

			fread(d[0]->aln[j].a, sizeof(bwt_aln1_t), n_aln, fp_sa[j]);

			kv_copy(bwt_aln1_t, buf[j][i].aln, d[0]->aln[j]); // backup d[0]->aln[j]

			// generate SE alignment and mapping quality
			bwa_aln2seq(n_aln, d[0]->aln[j].a, p[j]);
		}
	}

	read_aln_clocks = clock() - t;

	fprintf(stderr, "%.2f sec\n", (float)(read_aln_clocks) / CLOCKS_PER_SEC); t = clock();
#endif
	// ---------------

	fprintf(stderr, "[bwa_sai2sam_pe_core] convert to sequence coordinate... \n");

	if (_bwt[0] == 0) { // load forward SA
		strcpy(str, prefix); strcat(str, ".bwt");  bwt[0] = bwt_restore_bwt(str);
		strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt[0]);
		strcpy(str, prefix); strcat(str, ".rbwt"); bwt[1] = bwt_restore_bwt(str);
		strcpy(str, prefix); strcat(str, ".rsa"); bwt_restore_sa(str, bwt[1]);
	} else bwt[0] = _bwt[0], bwt[1] = _bwt[1];

	// tpx

        // -----------------

#ifdef HAVE_PTHREAD
        if(num_sampe_threads > 1){

          delta = n_seqs / num_sampe_threads;
       
          for(i=0;i<num_sampe_threads;i++){
            thr_bwa_se_info[i].end = delta * (i+1);
            thr_bwa_se_info[i].bwt[0] = bwt[0];
            thr_bwa_se_info[i].bwt[1] = bwt[1];
            thr_bwa_se_info[i].seqs[0] = seqs[0];
            thr_bwa_se_info[i].seqs[1] = seqs[1];
            thr_bwa_se_info[i].gopt = gopt;
          }
       
          thr_bwa_se_info[num_sampe_threads-1].end = n_seqs;
       
          thr_bwa_se_info[0].start = 0;
       
          for(i=1;i<num_sampe_threads;i++){
            thr_bwa_se_info[i].start = thr_bwa_se_info[i-1].end;
          }
       
          for(i=0;i<num_sampe_threads;i++){
            srtn = pthread_create(&tid,NULL,(void *(*)(void *))thr_bwa_se_tpx,(void *)(long)i);
            if(srtn != 0){
              fprintf(stderr,"[%s] pthread_create thr_bwa_se_tpx error %d\n", __func__, srtn);
              exit(1);
            }
            thr_bwa_se_info[i].tid = tid;
          }
       
          for(i=0;i<num_sampe_threads;i++){
            pthread_join(thr_bwa_se_info[i].tid,NULL);
          }

        }else{

	  bwa_se_tpx(0, bwt, 0, n_seqs, seqs, gopt);

	}
#else // HAVE_PTHREAD
	bwa_se_tpx(0, bwt, 0, n_seqs, seqs, gopt);
#endif // HAVE_PTHREAD

        // -----------------

	// tpx

	// infer isize
	infer_isize(n_seqs, seqs, ii, opt->ap_prior, bwt[0]->seq_len);

	if (ii->avg < 0.0 && last_ii->avg > 0.0) *ii = *last_ii;

	if (opt->force_isize) {
		fprintf(stderr, "[%s] discard insert size estimate as user's request.\n", __func__);
		ii->low = ii->high = 0; ii->avg = ii->std = -1.0;
	}

	// tpx

        // -----------------

	// PE
#ifdef HAVE_PTHREAD
        if(num_sampe_threads > 1){

	  for(i=1; i<num_sampe_threads; i++){
            for(j=0; j<2; j++){
              d[i]->aln[j].n = d[0]->aln[j].n;
              kv_copy(bwt_aln1_t, d[i]->aln[j], d[0]->aln[j]);
            }
	  }

          delta = n_seqs / num_sampe_threads;
       
          for(i=0;i<num_sampe_threads;i++){
            thr_bwa_pe_info[i].end = delta * (i+1);
            thr_bwa_pe_info[i].cnt_chg = 0;
            thr_bwa_pe_info[i].bwt[0] = bwt[0];
            thr_bwa_pe_info[i].bwt[1] = bwt[1];
            thr_bwa_pe_info[i].seqs[0] = seqs[0];
            thr_bwa_pe_info[i].seqs[1] = seqs[1];
            thr_bwa_pe_info[i].ii = ii;
            thr_bwa_pe_info[i].opt = opt;
            thr_bwa_pe_info[i].gopt = gopt;
            for(j=0;j<num_sampe_threads;j++){
              thr_bwa_pe_info[i].d[j] = d[j];
            }
            thr_bwa_pe_info[i].buf[0] = buf[0];
            thr_bwa_pe_info[i].buf[1] = buf[1];
          }
       
          thr_bwa_pe_info[num_sampe_threads-1].end = n_seqs;
       
          thr_bwa_pe_info[0].start = 0;
       
          for(i=1;i<num_sampe_threads;i++){
            thr_bwa_pe_info[i].start = thr_bwa_pe_info[i-1].end;
          }
       
          for(i=0;i<num_sampe_threads;i++){
            srtn = pthread_create(&tid,NULL,(void *(*)(void *))thr_bwa_pe_tpx,(void *)(long)i);
            if(srtn != 0){
              fprintf(stderr,"[%s] pthread_create thr_bwa_pe_tpx error %d\n", __func__, srtn);
              exit(1);
            }
            thr_bwa_pe_info[i].tid = tid;
          }
       
          for(i=0;i<num_sampe_threads;i++){
            pthread_join(thr_bwa_pe_info[i].tid,NULL);
            cnt_chg += thr_bwa_pe_info[i].cnt_chg;
          }

        }else{

	  cnt_chg = bwa_pe_tpx(0, bwt, 0, n_seqs, seqs, ii, opt, gopt, d, buf);

        }
#else // HAVE_PTHREAD
	cnt_chg = bwa_pe_tpx(0, bwt, 0, n_seqs, seqs, ii, opt, gopt, d, buf);
#endif // HAVE_PTHREAD

        // -----------------

	// tpx

	// free
	for (i = 0; i < n_seqs; ++i) {
		kv_destroy(buf[0][i].aln);
		kv_destroy(buf[1][i].aln);
	}

	if (_bwt[0] == 0) {
		bwt_destroy(bwt[0]); bwt_destroy(bwt[1]);
	}

        for(i=0;i<num_sampe_threads;i++){
		kv_destroy(d[i]->arr);
		kv_destroy(d[i]->pos[0]); kv_destroy(d[i]->pos[1]);
		kv_destroy(d[i]->aln[0]); kv_destroy(d[i]->aln[1]);
	}

	return cnt_chg;
}

// tpx #define SW_MIN_MATCH_LEN 20
// tpx #define SW_MIN_MAPQ 17

// cnt = n_mm<<16 | n_gapo<<8 | n_gape
bwa_cigar_t *bwa_sw_core(bwtint_t l_pac, const ubyte_t *pacseq, int len, const ubyte_t *seq, int64_t *beg, int reglen,
					  int *n_cigar, uint32_t *_cnt)
{
	bwa_cigar_t *cigar = 0;
	ubyte_t *ref_seq;
	bwtint_t k, x, y, l;
	int path_len, ret;
	AlnParam ap = aln_param_bwa;
	path_t *path, *p;

	// check whether there are too many N's
	if (reglen < SW_MIN_MATCH_LEN || (int64_t)l_pac - *beg < len) return 0;
	for (k = 0, x = 0; k < len; ++k)
		if (seq[k] >= 4) ++x;
	if ((float)x/len >= 0.25 || len - x < SW_MIN_MATCH_LEN) return 0;

	// get reference subsequence
	ref_seq = (ubyte_t*)calloc(reglen, 1);
	for (k = *beg, l = 0; l < reglen && k < l_pac; ++k)
		ref_seq[l++] = pacseq[k>>2] >> ((~k&3)<<1) & 3;
	path = (path_t*)calloc(l+len, sizeof(path_t));

	// do alignment
	ret = aln_local_core(ref_seq, l, (ubyte_t*)seq, len, &ap, path, &path_len, 1, 0);
	if (ret < 0) {
		free(path); free(cigar); free(ref_seq); *n_cigar = 0;
		return 0;
	}
	cigar = bwa_aln_path2cigar(path, path_len, n_cigar);

	// check whether the alignment is good enough
	for (k = 0, x = y = 0; k < *n_cigar; ++k) {
		bwa_cigar_t c = cigar[k];
		if (__cigar_op(c) == FROM_M) x += __cigar_len(c), y += __cigar_len(c);
		else if (__cigar_op(c) == FROM_D) x += __cigar_len(c);
		else y += __cigar_len(c);
	}
	if (x < SW_MIN_MATCH_LEN || y < SW_MIN_MATCH_LEN) { // not good enough
		free(path); free(cigar); free(ref_seq);
		*n_cigar = 0;
		return 0;
	}

	{ // update cigar and coordinate;
		int start, end;
		p = path + path_len - 1;
		*beg += (p->i? p->i : 1) - 1;
		start = (p->j? p->j : 1) - 1;
		end = path->j;
		cigar = (bwa_cigar_t*)realloc(cigar, sizeof(bwa_cigar_t) * (*n_cigar + 2));
		if (start) {
			memmove(cigar + 1, cigar, sizeof(bwa_cigar_t) * (*n_cigar));
			cigar[0] = __cigar_create(3, start);
			++(*n_cigar);
		}
		if (end < len) {
			/*cigar[*n_cigar] = 3<<14 | (len - end);*/
			cigar[*n_cigar] = __cigar_create(3, (len - end));
			++(*n_cigar);
		}
	}

	{ // set *cnt
		int n_mm, n_gapo, n_gape;
		n_mm = n_gapo = n_gape = 0;
		p = path + path_len - 1;
		x = p->i? p->i - 1 : 0; y = p->j? p->j - 1 : 0;
		for (k = 0; k < *n_cigar; ++k) {
			bwa_cigar_t c = cigar[k];
			if (__cigar_op(c) == FROM_M) {
				for (l = 0; l < (__cigar_len(c)); ++l)
					if (ref_seq[x+l] < 4 && seq[y+l] < 4 && ref_seq[x+l] != seq[y+l]) ++n_mm;
				x += __cigar_len(c), y += __cigar_len(c);
			} else if (__cigar_op(c) == FROM_D) {
				x += __cigar_len(c), ++n_gapo, n_gape += (__cigar_len(c)) - 1;
			} else if (__cigar_op(c) == FROM_I) {
				y += __cigar_len(c), ++n_gapo, n_gape += (__cigar_len(c)) - 1;
			}
		}
		*_cnt = (uint32_t)n_mm<<16 | n_gapo<<8 | n_gape;
	}
	
	free(ref_seq); free(path);
	return cigar;
}

// -------------------

void thr_bwa_sw_tpx(long idx)
{
  int iidx = (int)idx;

  bwa_sw_tpx(iidx,
             thr_bwa_sw_info[iidx].bns,
             thr_bwa_sw_info[iidx].pacseq,
             thr_bwa_sw_info[iidx].start,
             thr_bwa_sw_info[iidx].end,
             thr_bwa_sw_info[iidx].seqs,
             thr_bwa_sw_info[iidx].popt,
             thr_bwa_sw_info[iidx].ii);

  return;
}

// -------------------

ubyte_t *bwa_paired_sw(const bntseq_t *bns, const ubyte_t *_pacseq, int n_seqs, bwa_seq_t *seqs[2], 
                       const pe_opt_t *popt, const isize_info_t *ii)
{
	ubyte_t *pacseq;
	// tpx int i;
#ifdef HAVE_PTHREAD
	int i = 0;
	int srtn = 0;
	long delta = 0L;
	pthread_t tid;
#endif // HAVE_PTHREAD

	// load reference sequence
	if (_pacseq == 0) {
		pacseq = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
		rewind(bns->fp_pac);
		fread(pacseq, 1, bns->l_pac/4+1, bns->fp_pac);
	} else pacseq = (ubyte_t*)_pacseq;
	if (!popt->is_sw || ii->avg < 0.0) return pacseq;

	// perform mate alignment
	n_tot[0] = n_tot[1] = n_mapped[0] = n_mapped[1] = 0;

	// tpx

        // -----------------

#ifdef HAVE_PTHREAD
        if(num_sampe_threads > 1){

          delta = n_seqs / num_sampe_threads;
       
          for(i=0;i<num_sampe_threads;i++){
            thr_bwa_sw_info[i].end = delta * (i+1);
            thr_bwa_sw_info[i].bns = bns;
            thr_bwa_sw_info[i].pacseq = pacseq;
            thr_bwa_sw_info[i].seqs[0] = seqs[0];
            thr_bwa_sw_info[i].seqs[1] = seqs[1];
            thr_bwa_sw_info[i].popt = popt;
            thr_bwa_sw_info[i].ii = ii;
          }
       
          thr_bwa_sw_info[num_sampe_threads-1].end = n_seqs;
       
          thr_bwa_sw_info[0].start = 0;
       
          for(i=1;i<num_sampe_threads;i++){
            thr_bwa_sw_info[i].start = thr_bwa_sw_info[i-1].end;
          }
       
          for(i=0;i<num_sampe_threads;i++){
            srtn = pthread_create(&tid,NULL,(void *(*)(void *))thr_bwa_sw_tpx,(void *)(long)i);
            if(srtn != 0){
              fprintf(stderr,"[%s] pthread_create thr_bwa_sw_tpx error %d\n", __func__, srtn);
              exit(1);
            }
            thr_bwa_sw_info[i].tid = tid;
          }
       
          for(i=0;i<num_sampe_threads;i++){
            pthread_join(thr_bwa_sw_info[i].tid,NULL);
          }

        }else{

          bwa_sw_tpx(0, bns, pacseq, 0, n_seqs, seqs, popt, ii);

        }
#else // HAVE_PTHREAD
        bwa_sw_tpx(0, bns, pacseq, 0, n_seqs, seqs, popt, ii);
#endif // HAVE_PTHREAD

        // -----------------

	// tpx

	fprintf(stderr, "[bwa_paired_sw] %lld out of %lld Q%d singletons are mated.\n",
			(long long)n_mapped[1], (long long)n_tot[1], SW_MIN_MAPQ);
	fprintf(stderr, "[bwa_paired_sw] %lld out of %lld Q%d discordant pairs are fixed.\n",
			(long long)n_mapped[0], (long long)n_tot[0], SW_MIN_MAPQ);
	return pacseq;
}

void bwa_sai2sam_pe_core(const char *prefix, char *const fn_sa[2], char *const fn_fa[2], pe_opt_t *popt)
{
	int i, j, n_seqs[3], tot_seqs = 0;
	bwa_seq_t *seqs[3][2];
	bwa_seqio_t *ks[2];
	clock_t t;
	clock_t t2;
	bntseq_t *bns, *ntbns = 0;
	FILE *fp_sa[2];
	gap_opt_t opt, opt0;
	isize_info_t last_ii; // this is for the last batch of reads
	char str[1024];
	bwt_t *bwt[2];
	uint8_t *pac;
	int n_needed;
	int nexti1;
	int nexti2;
	int seqsd[3];
	int max_threads = 1;
	clock_t tio;
	int first = 1;
	pe_data_t *d[3][MAX_CPUS];
	aln_buf_t *alnbuf[3][2];

	extern bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa);
      
#ifdef _SC_NPROCESSORS_ONLN
	max_threads = sysconf(_SC_NPROCESSORS_ONLN);
#else
	max_threads = MAX_CPUS;
#endif
 
	if(max_threads > MAX_CPUS)
		max_threads = MAX_CPUS;
       
	if(max_threads < 1)
		max_threads = 1;
       
	if(num_sampe_threads > max_threads)
		num_sampe_threads = max_threads;

	if(num_sampe_threads < 1)
		num_sampe_threads = max_threads;

	n_needed = 262144;
	if(adj_n_needed)
		n_needed = 1048576;

	fprintf(stderr, "[bwa_sai2sam_pe_core] version: %s (%s)\n",
			bwaversionstr, bwablddatestr);
	fprintf(stderr, "[bwa_sai2sam_pe_core] num threads: %d (max: %d)\n", 
			num_sampe_threads, max_threads);

	// initialization
	bwase_initialize(); // initialize g_log_n[] in bwase.c
	pac = 0; bwt[0] = bwt[1] = 0;
	for (i = 1; i != 256; ++i) g_log_n[i] = (int)(4.343 * log(i) + 0.5);
	bns = bns_restore(prefix);
	srand48(bns->seed);
	fp_sa[0] = xopen(fn_sa[0], "r");
	fp_sa[1] = xopen(fn_sa[1], "r");
#ifndef _USE_LOCAL_GHASH
	for(i=0; i<num_sampe_threads; i++){
		g_hash[i] = kh_init(64);
	}
#endif // ! _USE_LOCAL_GHASH
	last_ii.avg = -1.0;

	fread(&opt, sizeof(gap_opt_t), 1, fp_sa[0]);
	ks[0] = bwa_open_reads(opt.mode, fn_fa[0]);
	opt0 = opt;
	fread(&opt, sizeof(gap_opt_t), 1, fp_sa[1]); // overwritten!
	ks[1] = bwa_open_reads(opt.mode, fn_fa[1]);
	if (!(opt.mode & BWA_MODE_COMPREAD)) {
		popt->type = BWA_PET_SOLID;
		ntbns = bwa_open_nt(prefix);
	} else { // for Illumina alignment only
		if (popt->is_preload) {
			strcpy(str, prefix); strcat(str, ".bwt");  bwt[0] = bwt_restore_bwt(str);
			strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt[0]);
			strcpy(str, prefix); strcat(str, ".rbwt"); bwt[1] = bwt_restore_bwt(str);
			strcpy(str, prefix); strcat(str, ".rsa"); bwt_restore_sa(str, bwt[1]);
			pac = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
			rewind(bns->fp_pac);
			fread(pac, 1, bns->l_pac/4+1, bns->fp_pac);
		}
	}

	if( (opt0.mode & BWA_MODE_IL13) && (opt.mode & BWA_MODE_IL13) ){
		// fprintf(stderr,"[bwa_sai2sam_pe_core] aln files in Illumina 1.3+ format\n");
		soap2_qual_adj = 64;
	}

	// core loop
	bwa_print_sam_SQ(bns);
	bwa_print_sam_PG();

	first = 1;
	seqsd[0] = 0;
	seqsd[1] = 0;
	seqsd[2] = 0;
	nexti1 = 0;
	nexti2 = 0;
	n_seqs[0] = 0;
	n_seqs[1] = 0;
	n_seqs[2] = 0;
	seqs[0][0] = NULL;
	seqs[0][1] = NULL;
	seqs[1][0] = NULL;
	seqs[1][1] = NULL;
	seqs[2][0] = NULL;
	seqs[2][1] = NULL;
	alnbuf[0][0] = NULL;
	alnbuf[0][1] = NULL;
	alnbuf[1][0] = NULL;
	alnbuf[1][1] = NULL;
	alnbuf[2][0] = NULL;
	alnbuf[2][1] = NULL;

	for(i=0; i<num_sampe_threads; i++){
		d[0][i] = (pe_data_t*)calloc(1, sizeof(pe_data_t));
		d[1][i] = (pe_data_t*)calloc(1, sizeof(pe_data_t));
		d[2][i] = (pe_data_t*)calloc(1, sizeof(pe_data_t));
	}

	t = tio = clock();

	bwa_read_seq2_tpx(ks[0], ks[1], n_needed, &n_seqs[nexti1], opt0.mode, opt.mode, 
				opt0.trim_qual, opt.trim_qual, &seqs[nexti1][0], &seqs[nexti1][1], 
				d[nexti1], &alnbuf[nexti1][0], &alnbuf[nexti1][1], nexti1, fp_sa);

        while(1){

		int cnt_chg;
		isize_info_t ii;
		ubyte_t *pacseq;

                // ---------------

		if( ( (async_read_seq) && (num_sampe_threads > 1) ) || (!first) ){
			tio = clock();
		}

		bwa_read_seq2_wait_tpx();

		if(seqs[nexti1][0] == NULL){
			break;
		}

		tot_seqs += n_seqs[nexti1];

		seqsd[nexti1] = 1;

		nexti2 = nexti1 + 1;
		if(nexti2 > 2){
			nexti2 = 0;
		}

		bwa_read_seq2_tpx(ks[0], ks[1], n_needed, &n_seqs[nexti2], opt0.mode, opt.mode, 
					opt0.trim_qual, opt.trim_qual, &seqs[nexti2][0], &seqs[nexti2][1], 
					d[nexti2], &alnbuf[nexti2][0], &alnbuf[nexti2][1], nexti2, fp_sa);

		fprintf(stderr, "[bwa_sai2sam_pe_core] bwa_read_seq2... %.2f sec", (float)(clock() - tio) / CLOCKS_PER_SEC);
		if( (async_read_seq) && (num_sampe_threads > 1) && (!first) ){
			fprintf(stderr," (async, bsize=%dk)\n", n_needed / 1024);
		}else{
			fprintf(stderr," (bsize=%dk)\n", n_needed / 1024);
		}

		first = 0;

                // ---------------

		t = clock();

		cnt_chg = bwa_cal_pac_pos_pe(prefix, bwt, n_seqs[nexti1], seqs[nexti1], fp_sa, &ii, popt, 
						&opt, &last_ii, d[nexti1], alnbuf[nexti1]);
		t2 = clock() - read_aln_clocks;
		fprintf(stderr, "[bwa_sai2sam_pe_core] time elapses: %.2f sec\n", (float)(t2 - t) / CLOCKS_PER_SEC); t = clock();

                // ---------------

		fprintf(stderr, "[bwa_sai2sam_pe_core] changing coordinates of %d alignments.\n", cnt_chg);

		fprintf(stderr, "[bwa_sai2sam_pe_core] align unmapped mate...\n");
		pacseq = bwa_paired_sw(bns, pac, n_seqs[nexti1], seqs[nexti1], popt, &ii);
		fprintf(stderr, "[bwa_sai2sam_pe_core] time elapses: %.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

                // ---------------

		fprintf(stderr, "[bwa_sai2sam_pe_core] refine gapped alignments... ");
		for (j = 0; j < 2; ++j)
			bwa_refine_gapped(bns, n_seqs[nexti1], seqs[nexti1][j], pacseq, ntbns);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		if (pac == 0) free(pacseq);

		// ---------------

		fprintf(stderr, "[bwa_sai2sam_pe_core] print alignments... ");
		bwa_print2_tpx(bns, n_seqs[nexti1], seqs[nexti1], opt);
		fprintf(stderr, "%.2f sec", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();
		if( (async_print_res) && (num_sampe_threads > 1) ){
			fprintf(stderr," (async)\n");
		}else{
			fprintf(stderr,"\n");
		}

		// ---------------

		fprintf(stderr, "[bwa_sai2sam_pe_core] %d sequences have been processed.\n", tot_seqs);

		last_ii = ii;

		// ---------------

		nexti1 = nexti2;

		nexti2 = nexti1 + 1;
		if(nexti2 > 2){
			nexti2 = 0;
		}

		if(seqsd[nexti2]){
			for (j = 0; j < 2; ++j){
				bwa_free_read_seq(n_seqs[nexti2], seqs[nexti2][j]);
			}
			seqsd[nexti2] = 0;
		}

	}

	// ---------------

	if( (async_print_res) && (num_sampe_threads > 1) ){
		t = clock();
		fprintf(stderr, "[bwa_sai2sam_pe_core] wait for final print alignments... ");
	}

	bwa_print2_wait_tpx();

	if( (async_print_res) && (num_sampe_threads > 1) ){
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();
	}

	for(nexti1=0; nexti1<3; nexti1++){
		if(seqsd[nexti1]){
			for (j = 0; j < 2; ++j){
				bwa_free_read_seq(n_seqs[nexti1], seqs[nexti1][j]);
				if(alnbuf[nexti1][j] != NULL){
					free(alnbuf[nexti1][j]);
					alnbuf[nexti1][j] = NULL;
				}
			}
			seqsd[nexti1] = 0;
		}
	}

	for(i=0; i<num_sampe_threads; i++){
		free(d[0][i]);
		free(d[1][i]);
		free(d[2][i]);
	}

	// ---------------

	// destroy
	bns_destroy(bns);
	if (ntbns) bns_destroy(ntbns);
	for (i = 0; i < 2; ++i) {
		bwa_seq_close(ks[i]);
		fclose(fp_sa[i]);
	}

#ifndef _USE_LOCAL_GHASH
	for(i=0; i<num_sampe_threads; i++){
		khint_t iter;
		for (iter = kh_begin(g_hash[i]); iter != kh_end(g_hash[i]); ++iter) {
			if (kh_exist(g_hash[i], iter)) 
				free(kh_val(g_hash[i], iter).a);
		}
		kh_destroy(64, g_hash[i]);
	}
#endif // ! _USE_LOCAL_GHASH

	if (pac) {
		free(pac); bwt_destroy(bwt[0]); bwt_destroy(bwt[1]);
	}

	return;
}

int bwa_sai2sam_pe(int argc, char *argv[])
{
	extern char *bwa_rg_line, *bwa_rg_id;
	extern int bwa_set_rg(const char *s);
	int c;
	pe_opt_t *popt;
	struct timeval st;
	uint64_t s1, e1;
	double total_time = 0.0;

	gettimeofday(&st, NULL);
	s1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;

	popt = bwa_init_pe_opt();

	while ((c = getopt(argc, argv, "a:t:o:ISTXYsPn:N:c:f:Ar:")) >= 0) {
		switch (c) {
		case 'r':
			if (bwa_set_rg(optarg) < 0) {
				fprintf(stderr, "[%s] malformated @RG line\n", __func__);
				return 1;
			}
			break;
		case 't': num_sampe_threads = atoi(optarg); break;
		case 'a': popt->max_isize = atoi(optarg); break;
		case 'o': popt->max_occ = atoi(optarg); break;
		case 's': popt->is_sw = 0; break;
		case 'P': popt->is_preload = 1; break;
		case 'n': popt->n_multi = atoi(optarg); break;
		case 'N': popt->N_multi = atoi(optarg); break;
		case 'c': popt->ap_prior = atof(optarg); break;
		case 'f': xreopen(optarg, "w", stdout); break;
		case 'A': popt->force_isize = 1; break;
		case 'T': adj_n_needed = 0; break;
		case 'X': async_read_seq = 0; break;
		case 'Y': async_print_res = 0; break;
		case 'S': use_soap2_format = 1; break;
		case 'I': soap2_qual_adj = 64; break;
		default: return 1;
		}
	}

	if (optind + 5 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa sampe [options] <prefix> <in1.sai> <in2.sai> <in1.fq> <in2.fq>\n\n");
		fprintf(stderr, "Options: -a INT   maximum insert size [%d]\n", popt->max_isize);
		fprintf(stderr, "         -o INT   maximum occurrences for one end [%d]\n", popt->max_occ);
		fprintf(stderr, "         -n INT   maximum hits to output for paired reads [%d]\n", popt->n_multi);
		fprintf(stderr, "         -N INT   maximum hits to output for discordant pairs [%d]\n", popt->N_multi);
		fprintf(stderr, "         -c FLOAT prior of chimeric rate (lower bound) [%.1le]\n", popt->ap_prior);
        	fprintf(stderr, "         -f FILE  sam file to output results to [stdout]\n");
		fprintf(stderr, "         -r STR   read group header line such as `@RG\\tID:foo\\tSM:bar' [null]\n");
		fprintf(stderr, "         -P       preload index into memory (for base-space reads only)\n");
		fprintf(stderr, "         -s       disable Smith-Waterman for the unmapped mate\n");
		fprintf(stderr, "         -A       disable insert size estimate (force -s)\n");
		fprintf(stderr, "         -t INT   number of threads [%d] (use <=0 for all)\n", num_sampe_threads);
		fprintf(stderr, "         -I       adjust quality to Illumina v1.3+ format\n");
		fprintf(stderr, "         -S       output SOAP2 format\n");
		fprintf(stderr, "         -T       use original read buffer size\n");
		fprintf(stderr, "         -X       disable async read seq/aln method\n");
		fprintf(stderr, "         -Y       disable async print results method\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "Notes: 1. For SOLiD reads, <in1.fq> corresponds R3 reads and <in2.fq> to F3.\n");
		fprintf(stderr, "       2. For reads shorter than 30bp, applying a smaller -o is recommended to\n");
		fprintf(stderr, "          to get a sensible speed at the cost of pairing accuracy.\n");
		fprintf(stderr, "\n");
		return 1;
	}

	bwa_sai2sam_pe_core(argv[optind], argv + optind + 1, argv + optind+3, popt);

	free(bwa_rg_line); free(bwa_rg_id);
	free(popt);

	// cant use getrusage for ru.maxrss until kernel 2.6.36 ...
	int srtn = 0;
	long maxrsskb = 0L;
	srtn = getmaxrss(&maxrsskb);
	if(srtn == 0){
		fprintf(stderr,"[bwa_sai2sam_pe] mem rss max = %ld (mb)\n",maxrsskb / 1024L);
	}

	gettimeofday(&st, NULL);
	e1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;
	total_time = (double)((double)e1 - (double)s1) / 1000000.0;

	fprintf(stderr,"[bwa_sai2sam_pe] total time = %lf (sec)\n",total_time);

	return 0;
}
