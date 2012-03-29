#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <sys/time.h>
#include <time.h>
#include "stdaln.h"
#include "bwase.h"
#include "bwtaln.h"
#include "bntseq.h"
#include "utils.h"
#include "kstring.h"
#include "bwatpx.h"

int g_log_n[256];
char *bwa_rg_line, *bwa_rg_id;

extern char bwaversionstr[];
extern char bwablddatestr[];

void bwa_print_sam_PG(void);
extern void bwa_rg_tpx(int iidx, const bntseq_t *bns, int n_seqs1, int n_seqs2,
                       bwa_seq_t *seqs, ubyte_t *pacseq, bntseq_t *ntbns);
extern void bwa_read_seq1_tpx(bwa_seqio_t *ks1, int n_needed, int *n,
                              int mode1, int trim_qual1,
                              bwa_seq_t **seq1, bwt_aln1_t **aln, int n_occ, FILE *fp_sa, int m_aln);
extern void bwa_read_seq1_wait_tpx(void);
extern void bwa_print1_tpx(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, const gap_opt_t opt);
extern void bwa_print1_wait_tpx(void);
extern void bwa_cal_pac_posf2_tpx(const bwt_t *bwt,int n_seqs1,int n_seqs2,
                                 bwa_seq_t *seqs,const int max_mm, const float fnr);
extern void bwa_cal_pac_posr2_tpx(const bwt_t *bwt,int n_seqs1,int n_seqs2,
                                 bwa_seq_t *seqs,const int max_mm, const float fnr);

extern int num_sampe_threads;
THR_BWA_RG_TPX thr_bwa_rg_info[MAX_CPUS];
THR_BWA_SE_PAC_TPX thr_bwa_se_pac_info[MAX_CPUS];

extern int adj_n_needed;
extern int async_read_seq;
extern int async_print_res;
extern int use_soap2_format;
extern int soap2_qual_adj;

static char *samf1 = NULL;
static int samf1_len = 0;
static char *samf3 = NULL;
static int samf3_len = 0;
static char *samf6 = NULL;
static int samf6_len = 0;
static char *samf7 = NULL;
static int samf7_len = 0;
static char *samf10 = NULL;
static int samf10_len = 0;
static char *samf11 = NULL;
static int samf11_len = 0;
static char *sam_cigar = NULL;
static int samc_len = 0;
static char *sam_mdzstr = NULL;
static int sammdz_len = 0;
static char *sam_mutstr = NULL;
static int sammut_len = 0;

// -------------------

void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi)
{
	int i, cnt, best;
	if (n_aln == 0) {
		s->type = BWA_TYPE_NO_MATCH;
		s->c1 = s->c2 = 0;
		return;
	}

	if (set_main) {
		best = aln[0].score;
		for (i = cnt = 0; i < n_aln; ++i) {
			const bwt_aln1_t *p = aln + i;
			if (p->score > best) break;
			if (drand48() * (p->l - p->k + 1 + cnt) > (double)cnt) {
				s->n_mm = p->n_mm; s->n_gapo = p->n_gapo; s->n_gape = p->n_gape; s->strand = p->a;
				s->score = p->score;
				s->sa = p->k + (bwtint_t)((p->l - p->k + 1) * drand48());
			}
			cnt += p->l - p->k + 1;
		}
		s->c1 = cnt;
		for (; i < n_aln; ++i) cnt += aln[i].l - aln[i].k + 1;
		s->c2 = cnt - s->c1;
		s->type = s->c1 > 1? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE;
	}

	if (n_multi) {
		int k, rest, n_occ, z = 0;
		for (k = n_occ = 0; k < n_aln; ++k) {
			const bwt_aln1_t *q = aln + k;
			n_occ += q->l - q->k + 1;
		}
		if (s->multi) free(s->multi);
		if (n_occ > n_multi + 1) { // if there are too many hits, generate none of them
			s->multi = 0; s->n_multi = 0;
			return;
		}
		/* The following code is more flexible than what is required
		 * here. In principle, due to the requirement above, we can
		 * simply output all hits, but the following samples "rest"
		 * number of random hits. */
		rest = n_occ > n_multi + 1? n_multi + 1 : n_occ; // find one additional for ->sa
		s->multi = calloc(rest, sizeof(bwt_multi1_t));
		for (k = 0; k < n_aln; ++k) {
			const bwt_aln1_t *q = aln + k;
			if (q->l - q->k + 1 <= rest) {
				bwtint_t l;
				for (l = q->k; l <= q->l; ++l) {
					s->multi[z].pos = l;
					s->multi[z].gap = q->n_gapo + q->n_gape;
					s->multi[z].mm = q->n_mm;
					s->multi[z++].strand = q->a;
				}
				rest -= q->l - q->k + 1;
			} else { // Random sampling (http://code.activestate.com/recipes/272884/). In fact, we never come here. 
				int j, i, k;
				for (j = rest, i = q->l - q->k + 1, k = 0; j > 0; --j) {
					double p = 1.0, x = drand48();
					while (x < p) p -= p * j / (i--);
					s->multi[z].pos = q->l - i;
					s->multi[z].gap = q->n_gapo + q->n_gape;
					s->multi[z].mm = q->n_mm;
					s->multi[z++].strand = q->a;
				}
				rest = 0;
				break;
			}
		}
		s->n_multi = z;
		for (k = z = 0; k < s->n_multi; ++k)
			if (s->multi[k].pos != s->sa)
				s->multi[z++] = s->multi[k];
		s->n_multi = z < n_multi? z : n_multi;
	}
}

void bwa_aln2seq(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s)
{
	bwa_aln2seq_core(n_aln, aln, s, 1, 0);
}

int bwa_approx_mapQ(const bwa_seq_t *p, int mm)
{
	int n;
	if (p->c1 == 0) return 23;
	if (p->c1 > 1) return 0;
	if (p->n_mm == mm) return 25;
	if (p->c2 == 0) return 37;
	n = (p->c2 >= 255)? 255 : p->c2;
	return (23 < g_log_n[n])? 0 : 23 - g_log_n[n];
}

// -------------------

void thr_bwa_se_pacf_tpx(long idx)
{
  int iidx = (int)idx;

  bwa_cal_pac_posf2_tpx(thr_bwa_se_pac_info[iidx].bwt,
             thr_bwa_se_pac_info[iidx].start,
             thr_bwa_se_pac_info[iidx].end,
             thr_bwa_se_pac_info[iidx].seqs,
             thr_bwa_se_pac_info[iidx].max_mm,
             thr_bwa_se_pac_info[iidx].fnr);

  return;
}

// -------------------

void thr_bwa_se_pacr_tpx(long idx)
{
  int iidx = (int)idx;

  bwa_cal_pac_posr2_tpx(thr_bwa_se_pac_info[iidx].bwt,
             thr_bwa_se_pac_info[iidx].start,
             thr_bwa_se_pac_info[iidx].end,
             thr_bwa_se_pac_info[iidx].seqs,
             thr_bwa_se_pac_info[iidx].max_mm,
             thr_bwa_se_pac_info[iidx].fnr);

  return;
}

// -------------------

void bwa_cal_pac_posf_tpx(const bwt_t *forward_bwt,int n_seqs, bwa_seq_t *seqs,const int max_mm, const float fnr)
{
#ifdef HAVE_PTHREAD
	int i;
        int srtn = 0;
        long delta = 0L;
        pthread_t tid;
#endif // HAVE_PTHREAD

#ifdef HAVE_PTHREAD
        if(num_sampe_threads > 1){

          delta = n_seqs / num_sampe_threads;

          for(i=0;i<num_sampe_threads;i++){
            thr_bwa_se_pac_info[i].end = delta * (i+1);
            thr_bwa_se_pac_info[i].bwt = forward_bwt;
            thr_bwa_se_pac_info[i].seqs = seqs;
            thr_bwa_se_pac_info[i].max_mm = max_mm;
            thr_bwa_se_pac_info[i].fnr = fnr;
          }

          thr_bwa_se_pac_info[num_sampe_threads-1].end = n_seqs;

          thr_bwa_se_pac_info[0].start = 0;

          for(i=1;i<num_sampe_threads;i++){
            thr_bwa_se_pac_info[i].start = thr_bwa_se_pac_info[i-1].end;
          }

          for(i=0;i<num_sampe_threads;i++){
            srtn = pthread_create(&tid,NULL,(void *(*)(void *))thr_bwa_se_pacf_tpx,(void *)(long)i);
            if(srtn != 0){
              fprintf(stderr,"[%s] pthread_create thr_bwa_se_pacf_tpx error %d\n", __func__, srtn);
              exit(1);
            }
            thr_bwa_se_pac_info[i].tid = tid;
          }

          for(i=0;i<num_sampe_threads;i++){
            pthread_join(thr_bwa_se_pac_info[i].tid,NULL);
          }

        }else{

	  bwa_cal_pac_posf2_tpx(forward_bwt,0,n_seqs,seqs,max_mm,fnr);

        }
#else // HAVE_PTHREAD
	bwa_cal_pac_posf2_tpx(forward_bwt,0,n_seqs,seqs,max_mm,fnr);
#endif // HAVE_PTHREAD

	return;
}

// -------------------

void bwa_cal_pac_posr_tpx(const bwt_t *reverse_bwt,int n_seqs, bwa_seq_t *seqs,const int max_mm, const float fnr)
{
#ifdef HAVE_PTHREAD
	int i;
        int srtn = 0;
        long delta = 0L;
        pthread_t tid;
#endif // HAVE_PTHREAD

#ifdef HAVE_PTHREAD
        if(num_sampe_threads > 1){

          delta = n_seqs / num_sampe_threads;

          for(i=0;i<num_sampe_threads;i++){
            thr_bwa_se_pac_info[i].end = delta * (i+1);
            thr_bwa_se_pac_info[i].bwt = reverse_bwt;
            thr_bwa_se_pac_info[i].seqs = seqs;
            thr_bwa_se_pac_info[i].max_mm = max_mm;
            thr_bwa_se_pac_info[i].fnr = fnr;
          }

          thr_bwa_se_pac_info[num_sampe_threads-1].end = n_seqs;

          thr_bwa_se_pac_info[0].start = 0;

          for(i=1;i<num_sampe_threads;i++){
            thr_bwa_se_pac_info[i].start = thr_bwa_se_pac_info[i-1].end;
          }

          for(i=0;i<num_sampe_threads;i++){
            srtn = pthread_create(&tid,NULL,(void *(*)(void *))thr_bwa_se_pacr_tpx,(void *)(long)i);
            if(srtn != 0){
              fprintf(stderr,"[%s] pthread_create thr_bwa_se_pacr_tpx error %d\n", __func__, srtn);
              exit(1);
            }
            thr_bwa_se_pac_info[i].tid = tid;
          }

          for(i=0;i<num_sampe_threads;i++){
            pthread_join(thr_bwa_se_pac_info[i].tid,NULL);
          }

        }else{

	  bwa_cal_pac_posr2_tpx(reverse_bwt,0,n_seqs,seqs,max_mm,fnr);

        }
#else // HAVE_PTHREAD
	bwa_cal_pac_posr2_tpx(reverse_bwt,0,n_seqs,seqs,max_mm,fnr);
#endif // HAVE_PTHREAD

	return;
}

/**
 * Derive the actual position in the read from the given suffix array
 * coordinates. Note that the position will be approximate based on
 * whether indels appear in the read and whether calculations are
 * performed from the start or end of the read.
 */
void bwa_cal_pac_pos_core(const bwt_t *forward_bwt, const bwt_t *reverse_bwt, bwa_seq_t *seq, const int max_mm, const float fnr)
{
	int max_diff;
	if (seq->type != BWA_TYPE_UNIQUE && seq->type != BWA_TYPE_REPEAT) return;
	max_diff = fnr > 0.0? bwa_cal_maxdiff(seq->len, BWA_AVG_ERR, fnr) : max_mm;
	if (seq->strand) { // reverse strand only
		seq->pos = bwt_sa(forward_bwt, seq->sa);
		seq->seQ = seq->mapQ = bwa_approx_mapQ(seq, max_diff);
	} else { // forward strand only
		/* NB: For gapped alignment, p->pos may not be correct, which
		 *     will be fixed in refine_gapped_core(). This line also
		 *     determines the way "x" is calculated in
		 *     refine_gapped_core() when (ext < 0 && is_end == 0). */
		seq->pos = reverse_bwt->seq_len - (bwt_sa(reverse_bwt, seq->sa) + seq->len);
		seq->seQ = seq->mapQ = bwa_approx_mapQ(seq, max_diff);
	}
}

void bwa_cal_pac_pos(const char *prefix, int n_seqs, bwa_seq_t *seqs, int max_mm, float fnr)
{
	char str[1024];
	bwt_t *bwt;

	// load forward SA
	strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
	strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt);

	bwa_cal_pac_posf_tpx(bwt, n_seqs, seqs, max_mm, fnr);

	bwt_destroy(bwt);

	// load reverse BWT and SA
	strcpy(str, prefix); strcat(str, ".rbwt"); bwt = bwt_restore_bwt(str);
	strcpy(str, prefix); strcat(str, ".rsa"); bwt_restore_sa(str, bwt);

	bwa_cal_pac_posr_tpx(bwt, n_seqs, seqs, max_mm, fnr);

	bwt_destroy(bwt);

	return;
}

/* is_end_correct == 1 if (*pos+len) gives the correct coordinate on
 * forward strand. This happens when p->pos is calculated by
 * bwa_cal_pac_pos(). is_end_correct==0 if (*pos) gives the correct
 * coordinate. This happens only for color-converted alignment. */
bwa_cigar_t *refine_gapped_core(bwtint_t l_pac, const ubyte_t *pacseq, int len, const ubyte_t *seq, bwtint_t *_pos,
									int ext, int *n_cigar, int is_end_correct)
{
	bwa_cigar_t *cigar = 0;
	ubyte_t *ref_seq;
	int l = 0, path_len, ref_len;
	AlnParam ap = aln_param_bwa;
	path_t *path;
	int64_t k, __pos = *_pos > l_pac? (int64_t)((int32_t)*_pos) : *_pos;

	ref_len = len + abs(ext);
	if (ext > 0) {
		ref_seq = (ubyte_t*)calloc(ref_len, 1);
		for (k = __pos; k < __pos + ref_len && k < l_pac; ++k)
			ref_seq[l++] = pacseq[k>>2] >> ((~k&3)<<1) & 3;
	} else {
		int64_t x = __pos + (is_end_correct? len : ref_len);
		ref_seq = (ubyte_t*)calloc(ref_len, 1);
		for (l = 0, k = x - ref_len > 0? x - ref_len : 0; k < x && k < l_pac; ++k)
			ref_seq[l++] = pacseq[k>>2] >> ((~k&3)<<1) & 3;
	}
	path = (path_t*)calloc(l+len, sizeof(path_t));

	aln_global_core(ref_seq, l, (ubyte_t*)seq, len, &ap, path, &path_len);
	cigar = bwa_aln_path2cigar(path, path_len, n_cigar);
	
	if (ext < 0 && is_end_correct) { // fix coordinate for reads mapped on the forward strand
		for (l = k = 0; k < *n_cigar; ++k) {
			if (__cigar_op(cigar[k]) == FROM_D) l -= __cigar_len(cigar[k]);
			else if (__cigar_op(cigar[k]) == FROM_I) l += __cigar_len(cigar[k]);
		}
		__pos += l;
	}

	if (__cigar_op(cigar[0]) == FROM_D) { // deletion at the 5'-end
		__pos += __cigar_len(cigar[0]);
		for (k = 0; k < *n_cigar - 1; ++k) cigar[k] = cigar[k+1];
		--(*n_cigar);
	}
	if (__cigar_op(cigar[*n_cigar-1]) == FROM_D) --(*n_cigar); // deletion at the 3'-end

	// change "I" at either end of the read to S. just in case. This should rarely happen...
	if (__cigar_op(cigar[*n_cigar-1]) == FROM_I) cigar[*n_cigar-1] = __cigar_create(3, (__cigar_len(cigar[*n_cigar-1])));
	if (__cigar_op(cigar[0]) == FROM_I) cigar[0] = __cigar_create(3, (__cigar_len(cigar[0])));

	*_pos = (bwtint_t)__pos;
	free(ref_seq); free(path);
	return cigar;
}

char *bwa_cal_md1(int n_cigar, bwa_cigar_t *cigar, int len, bwtint_t pos, ubyte_t *seq,
				  bwtint_t l_pac, ubyte_t *pacseq, kstring_t *str, int *_nm)
{
	bwtint_t x, y;
	int z, u, c, nm = 0;
	str->l = 0; // reset
	x = pos; y = 0;
	if (cigar) {
		int k, l;
		for (k = u = 0; k < n_cigar; ++k) {
			l = __cigar_len(cigar[k]);
			if (__cigar_op(cigar[k]) == FROM_M) {
				for (z = 0; z < l && x+z < l_pac; ++z) {
					c = pacseq[(x+z)>>2] >> ((~(x+z)&3)<<1) & 3;
					if (c > 3 || seq[y+z] > 3 || c != seq[y+z]) {
						ksprintf(str, "%d", u);
						kputc("ACGTN"[c], str);
						++nm;
						u = 0;
					} else ++u;
				}
				x += l; y += l;
/*		        } else if (cigar[k]>>14 == FROM_I || cigar[k]>>14 == 3) { */
                        } else if (__cigar_op(cigar[k]) == FROM_I || __cigar_op(cigar[k]) == FROM_S) {
				y += l;
				if (__cigar_op(cigar[k]) == FROM_I) nm += l;
			} else if (__cigar_op(cigar[k]) == FROM_D) {
				ksprintf(str, "%d", u);
				kputc('^', str);
				for (z = 0; z < l && x+z < l_pac; ++z)
					kputc("ACGT"[pacseq[(x+z)>>2] >> ((~(x+z)&3)<<1) & 3], str);
				u = 0;
				x += l; nm += l;
			}
		}
	} else { // no gaps
		for (z = u = 0; z < (bwtint_t)len && x+z < l_pac; ++z) {
			c = pacseq[(x+z)>>2] >> ((~(x+z)&3)<<1) & 3;
			if (c > 3 || seq[y+z] > 3 || c != seq[y+z]) {
				ksprintf(str, "%d", u);
				kputc("ACGTN"[c], str);
				++nm;
				u = 0;
			} else ++u;
		}
	}
	ksprintf(str, "%d", u);
	*_nm = nm;
	return strdup(str->s);
}

void bwa_correct_trimmed(bwa_seq_t *s)
{
	if (s->len == s->full_len) return;
	if (s->strand == 0) { // forward
		if (s->cigar && __cigar_op(s->cigar[s->n_cigar-1]) == FROM_S) { // the last is S
			s->cigar[s->n_cigar-1] += s->full_len - s->len;
		} else {
			if (s->cigar == 0) {
				s->n_cigar = 2;
				s->cigar = calloc(s->n_cigar, sizeof(bwa_cigar_t));
				s->cigar[0] = __cigar_create(0, s->len);
			} else {
				++s->n_cigar;
				s->cigar = realloc(s->cigar, s->n_cigar * sizeof(bwa_cigar_t));
			}
			s->cigar[s->n_cigar-1] = __cigar_create(3, (s->full_len - s->len));
		}
	} else { // reverse
		if (s->cigar && __cigar_op(s->cigar[0]) == FROM_S) { // the first is S
			s->cigar[0] += s->full_len - s->len;
		} else {
			if (s->cigar == 0) {
				s->n_cigar = 2;
				s->cigar = calloc(s->n_cigar, sizeof(bwa_cigar_t));
				s->cigar[1] = __cigar_create(0, s->len);
			} else {
				++s->n_cigar;
				s->cigar = realloc(s->cigar, s->n_cigar * sizeof(bwa_cigar_t));
				memmove(s->cigar + 1, s->cigar, (s->n_cigar-1) * sizeof(bwa_cigar_t));
			}
			s->cigar[0] = __cigar_create(3, (s->full_len - s->len));
		}
	}
	s->len = s->full_len;
}

// -------------------

void thr_bwa_rg_tpx(long idx)
{
  int iidx = (int)idx;

  bwa_rg_tpx(iidx,
             thr_bwa_rg_info[iidx].bns,
             thr_bwa_rg_info[iidx].start,
             thr_bwa_rg_info[iidx].end,
             thr_bwa_rg_info[iidx].seqs,
             thr_bwa_rg_info[iidx].pacseq,
             thr_bwa_rg_info[iidx].ntbns);

  return;
}

// -------------------

void bwa_refine_gapped(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, ubyte_t *_pacseq, bntseq_t *ntbns)
{
	ubyte_t *pacseq, *ntpac = 0;
#ifdef HAVE_PTHREAD
	int i;
	int srtn = 0;
	pthread_t tid;
	long delta = 0L;
#endif // HAVE_PTHREAD

	if (ntbns) { // in color space
		ntpac = (ubyte_t*)calloc(ntbns->l_pac/4+1, 1);
		rewind(ntbns->fp_pac);
		fread(ntpac, 1, ntbns->l_pac/4 + 1, ntbns->fp_pac);
	}

	if (!_pacseq) {
		pacseq = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
		rewind(bns->fp_pac);
		fread(pacseq, 1, bns->l_pac/4+1, bns->fp_pac);
	} else pacseq = _pacseq;

	// tpx

	// -----------------

#ifdef HAVE_PTHREAD
        if(num_sampe_threads > 1){

          delta = n_seqs / num_sampe_threads;
       
          for(i=0;i<num_sampe_threads;i++){
            thr_bwa_rg_info[i].end = delta * (i+1);
            thr_bwa_rg_info[i].bns = bns;
            thr_bwa_rg_info[i].seqs = seqs;
            thr_bwa_rg_info[i].pacseq = pacseq;
            thr_bwa_rg_info[i].ntbns = ntbns;
          }
       
          thr_bwa_rg_info[num_sampe_threads-1].end = n_seqs;
       
          thr_bwa_rg_info[0].start = 0;
       
          for(i=1;i<num_sampe_threads;i++){
            thr_bwa_rg_info[i].start = thr_bwa_rg_info[i-1].end;
          }
       
          for(i=0;i<num_sampe_threads;i++){
            srtn = pthread_create(&tid,NULL,(void *(*)(void *))thr_bwa_rg_tpx,(void *)(long)i);
            if(srtn != 0){
              fprintf(stderr,"[%s] pthread_create thr_bwa_rg_tpx error %d\n", __func__, srtn);
              exit(1);
            }
            thr_bwa_rg_info[i].tid = tid;
          }
       
          for(i=0;i<num_sampe_threads;i++){
            pthread_join(thr_bwa_rg_info[i].tid,NULL);
          }

        }else{

	  bwa_rg_tpx(0, bns, 0, n_seqs, seqs, pacseq, ntbns);

	}
#else // HAVE_PTHREAD
	bwa_rg_tpx(0, bns, 0, n_seqs, seqs, pacseq, ntbns);
#endif // HAVE_PTHREAD

	// -----------------

	// tpx

	if (!_pacseq) free(pacseq);

	free(ntpac);

	return;
}

int64_t pos_end(const bwa_seq_t *p)
{
	if (p->cigar) {
		int j;
		int64_t x = p->pos;
		for (j = 0; j != p->n_cigar; ++j) {
			int op = __cigar_op(p->cigar[j]);
			if (op == 0 || op == 2) x += __cigar_len(p->cigar[j]);
		}
		return x;
	} else return p->pos + p->len;
}

int64_t pos_end_multi(const bwt_multi1_t *p, int len) // analogy to pos_end()
{
	if (p->cigar) {
		int j;
		int64_t x = p->pos;
		for (j = 0; j != p->n_cigar; ++j) {
			int op = __cigar_op(p->cigar[j]);
			if (op == 0 || op == 2) x += __cigar_len(p->cigar[j]);
		}
		return x;
	} else return p->pos + len;
}

static int64_t pos_5(const bwa_seq_t *p)
{
	if (p->type != BWA_TYPE_NO_MATCH)
		return p->strand? pos_end(p) : p->pos;
	return -1;
}

void bwa_print_sam1a(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2)
{
	int j;

	if (p->type != BWA_TYPE_NO_MATCH || (mate && mate->type != BWA_TYPE_NO_MATCH)) {

		int seqid, nn, am = 0, flag = p->extra_flag;
		char XT;

		if (p->type == BWA_TYPE_NO_MATCH) {
			p->pos = mate->pos;
			p->strand = mate->strand;
			flag |= SAM_FSU;
			j = 1;
		} else j = pos_end(p) - p->pos; // j is the length of the reference in the alignment

		// get seqid
		nn = bns_coor_pac2real(bns, p->pos, j, &seqid);
		if (p->type != BWA_TYPE_NO_MATCH && p->pos + j - bns->anns[seqid].offset > bns->anns[seqid].len)
			flag |= SAM_FSU; // flag UNMAP as this alignment bridges two adjacent reference sequences

		// update flag and print it
		if (p->strand) flag |= SAM_FSR;
		if (mate) {
			if (mate->type != BWA_TYPE_NO_MATCH) {
				if (mate->strand) flag |= SAM_FMR;
			} else flag |= SAM_FMU;
		}

		printf("%s\t%d\t%s\t", p->name, flag, bns->anns[seqid].name);
		printf("%d\t%d\t", (int)(p->pos - bns->anns[seqid].offset + 1), p->mapQ);

		// print CIGAR
		if (p->cigar) {
			for (j = 0; j != p->n_cigar; ++j) {
				printf("%d%c", __cigar_len(p->cigar[j]), "MIDS"[__cigar_op(p->cigar[j])]);
			}
		} else if (p->type == BWA_TYPE_NO_MATCH) {
			printf("*");
		} else {
			printf("%dM", p->len);
		}

		// print mate coordinate
		if (mate && mate->type != BWA_TYPE_NO_MATCH) {
			int m_seqid, m_is_N;
			long long isize;
			am = mate->seQ < p->seQ? mate->seQ : p->seQ; // smaller single-end mapping quality
			// redundant calculation here, but should not matter too much
			m_is_N = bns_coor_pac2real(bns, mate->pos, mate->len, &m_seqid);
			printf("\t%s\t", (seqid == m_seqid)? "=" : bns->anns[m_seqid].name);
			isize = (seqid == m_seqid)? pos_5(mate) - pos_5(p) : 0;
			if (p->type == BWA_TYPE_NO_MATCH) isize = 0;
			printf("%d\t%lld\t", (int)(mate->pos - bns->anns[m_seqid].offset + 1), isize);
		} else if (mate) {
			printf("\t=\t%d\t0\t", (int)(p->pos - bns->anns[seqid].offset + 1));
		} else {
			printf("\t*\t0\t0\t");
		}

		// print sequence and quality
		if (p->strand == 0) {
			for (j = 0; j != p->full_len; ++j) {
				putchar("ACGTN"[(int)p->seq[j]]);
			}
		} else {
			for (j = 0; j != p->full_len; ++j) {
				putchar("TGCAN"[p->seq[p->full_len - 1 - j]]);
			}
		}

		putchar('\t');

		if (p->qual) {
			if (p->strand) {
				seq_reverse(p->len, p->qual, 0); // reverse quality
			}
			printf("%s", p->qual);
		} else {
			printf("*");
		}

		if (bwa_rg_id) {
			printf("\tRG:Z:%s", bwa_rg_id);
		}
		if (p->bc[0]) {
			printf("\tBC:Z:%s", p->bc);
		}
		if (p->clip_len < p->full_len) {
			printf("\tXC:i:%d", p->clip_len);
		}
		if (p->type != BWA_TYPE_NO_MATCH) {
			int i;
			// calculate XT tag
			XT = "NURM"[p->type];
			if (nn > 10) XT = 'N';
			// print tags
			printf("\tXT:A:%c\t%s:i:%d", XT, (mode & BWA_MODE_COMPREAD)? "NM" : "CM", p->nm);
			if (nn) printf("\tXN:i:%d", nn);
			if (mate) printf("\tSM:i:%d\tAM:i:%d", p->seQ, am);
			if (p->type != BWA_TYPE_MATESW) { // X0 and X1 are not available for this type of alignment
				printf("\tX0:i:%d", p->c1);
				if (p->c1 <= max_top2) printf("\tX1:i:%d", p->c2);
			}
			printf("\tXM:i:%d\tXO:i:%d\tXG:i:%d", p->n_mm, p->n_gapo, p->n_gapo+p->n_gape);
			if (p->md) printf("\tMD:Z:%s", p->md);
			// print multiple hits
			if (p->n_multi) {
				printf("\tXA:Z:");
				for (i = 0; i < p->n_multi; ++i) {
					bwt_multi1_t *q = p->multi + i;
					int k;
					j = pos_end_multi(q, p->len) - q->pos;
					nn = bns_coor_pac2real(bns, q->pos, j, &seqid);
					printf("%s,%c%d,", bns->anns[seqid].name, q->strand? '-' : '+',
						   (int)(q->pos - bns->anns[seqid].offset + 1));
					if (q->cigar) {
						for (k = 0; k < q->n_cigar; ++k)
							printf("%d%c", __cigar_len(q->cigar[k]), "MIDS"[__cigar_op(q->cigar[k])]);
					} else printf("%dM", p->len);
					printf(",%d;", q->gap + q->mm);
				}
			}
		}

		putchar('\n');

	} else { // this read has no match

		ubyte_t *s = p->strand? p->rseq : p->seq;
		int flag = p->extra_flag | SAM_FSU;
		if (mate && mate->type == BWA_TYPE_NO_MATCH) flag |= SAM_FMU;

		printf("%s\t%d\t*\t0\t0\t*\t*\t0\t0\t", p->name, flag);

		for (j = 0; j != p->len; ++j) {
			putchar("ACGTN"[(int)s[j]]);
		}

		putchar('\t');

		if (p->qual) {
			if (p->strand) {
				seq_reverse(p->len, p->qual, 0); // reverse quality
			}
			printf("%s", p->qual);
		} else {
			printf("*");
		}

		if (bwa_rg_id) {
			printf("\tRG:Z:%s", bwa_rg_id);
		}
		if (p->bc[0]) {
			printf("\tBC:Z:%s", p->bc);
		}
		if (p->clip_len < p->full_len) {
			printf("\tXC:i:%d", p->clip_len);
		}

		putchar('\n');
	}

	return;
}

// ------------------

void bwa_print_soap1a(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2)
{
	int i, j;
	int slen;
	int indx;
	int samf2, samf4, samf5, samf8, samf9;
	int samf2a, samf5a;
	char samf2b, samf2c;
	char tmpstr[2000];
	char tmpstr2[1000][500];
	// char tmp_cigar[2000];
	int num_mut;
	int inside_del;
	static int first = 1;

	if(first){
		first = 0;
		samc_len = 10000;
		sam_cigar = (char *)malloc(samc_len);
		if(sam_cigar == NULL){
			fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
			exit(1);
		}
		sammdz_len = 10000;
		sam_mdzstr = (char *)malloc(sammdz_len);
		if(sam_mdzstr == NULL){
			fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
			exit(1);
		}
		sammut_len = 10000;
		sam_mutstr = (char *)malloc(sammut_len);
		if(sam_mutstr == NULL){
			fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
			exit(1);
		}
	}

	num_mut = 0;
	// tmp_cigar[0] = '\0';
	sam_cigar[0] = '\0';
	sam_mdzstr[0] = '\0';
	sam_mutstr[0] = '\0';

	if (p->type != BWA_TYPE_NO_MATCH || (mate && mate->type != BWA_TYPE_NO_MATCH)) {

		int seqid, nn, am = 0, flag = p->extra_flag;

		if (p->type == BWA_TYPE_NO_MATCH) {
			p->pos = mate->pos;
			p->strand = mate->strand;
			flag |= SAM_FSU;
			j = 1;
		} else j = pos_end(p) - p->pos; // j is the length of the reference in the alignment

		// get seqid
		nn = bns_coor_pac2real(bns, p->pos, j, &seqid);
		if (p->type != BWA_TYPE_NO_MATCH && p->pos + j - bns->anns[seqid].offset > bns->anns[seqid].len)
			flag |= SAM_FSU; // flag UNMAP as this alignment bridges two adjacent reference sequences

		// update flag and print it
		if (p->strand) flag |= SAM_FSR;
		if (mate) {
			if (mate->type != BWA_TYPE_NO_MATCH) {
				if (mate->strand) flag |= SAM_FMR;
			} else flag |= SAM_FMU;
		}

		slen = (int)strlen(p->name);
		if(slen >= samf1_len){
			samf1_len = slen + 10000;
			samf1 = (char *)realloc(samf1, samf1_len);
			if(samf1 == NULL){
				fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
				exit(1);
			}
		}
		strcpy(samf1, p->name);

		samf2 = flag;

		slen = (int)strlen(bns->anns[seqid].name);
		if(slen >= samf3_len){
			samf3_len = slen + 10000;
			samf3 = (char *)realloc(samf3, samf3_len);
			if(samf3 == NULL){
				fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
				exit(1);
			}
		}
		strcpy(samf3, bns->anns[seqid].name);

		samf4 = (int)(p->pos - bns->anns[seqid].offset + 1);
		samf5 = p->mapQ;

#if 0
		// print CIGAR
		if (p->cigar) {
			for (j = 0; j != p->n_cigar; ++j) {
				printf("%d%c", __cigar_len(p->cigar[j]), "MIDS"[__cigar_op(p->cigar[j])]);
			}
		} else if (p->type == BWA_TYPE_NO_MATCH) {
			printf("*");
		} else {
			printf("%dM", p->len);
		}
#endif

		if (p->cigar) {
			indx = 0;
			for (j = 0; j != p->n_cigar; ++j) {
				// printf("%d%c", __cigar_len(p->cigar[j]), "MIDS"[__cigar_op(p->cigar[j])]);
				// sprintf(tmpstr,"%d%c", __cigar_len(p->cigar[j]), "MIDS"[__cigar_op(p->cigar[j])]);
				// strcat(tmp_cigar,tmpstr);
				indx += __cigar_len(p->cigar[j]);
				if ( ( __cigar_op(p->cigar[j]) == 1 ) || ( __cigar_op(p->cigar[j]) == 2 ) ) {
					if ( __cigar_op(p->cigar[j]) == 1 ) {        // I
						sprintf(tmpstr,"\t101\t%d", indx-1);
					} else if ( __cigar_op(p->cigar[j]) == 2 ) { // D
						sprintf(tmpstr,"\t201\t%d", indx-1);
					}
					slen += (int)strlen(tmpstr);
					if(slen >= samc_len){
						samc_len = slen + 10000;
						sam_cigar = (char *)realloc(sam_cigar, samc_len);
						if(sam_cigar == NULL){
							fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
							exit(1);
						}
					}
					strcat(sam_cigar, tmpstr);
				}
			}
		}else{
			// sprintf(tmp_cigar,"%d",p->len);
		}

		// print mate coordinate
		if (mate && mate->type != BWA_TYPE_NO_MATCH) {

			int m_seqid, m_is_N;
			long long isize;
			am = mate->seQ < p->seQ? mate->seQ : p->seQ; // smaller single-end mapping quality
			// redundant calculation here, but should not matter too much
			m_is_N = bns_coor_pac2real(bns, mate->pos, mate->len, &m_seqid);

			if(seqid == m_seqid){
				slen = 1;
			}else{
				slen = (int)strlen(bns->anns[m_seqid].name);
			}
			if(slen >= samf7_len){
				samf7_len = slen + 10000;
				samf7 = (char *)realloc(samf7, samf7_len);
				if(samf7 == NULL){
					fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
					exit(1);
				}
			}
			if(seqid == m_seqid){
				strcpy(samf7, "=");
			}else{
				strcpy(samf7, bns->anns[m_seqid].name);
			}

			isize = (seqid == m_seqid)? pos_5(mate) - pos_5(p) : 0;
			if (p->type == BWA_TYPE_NO_MATCH) isize = 0;

			samf8 = (int)(mate->pos - bns->anns[m_seqid].offset + 1);
			samf9 = isize;

		} else if (mate) {

			slen = 1;
			if(slen >= samf7_len){
				samf7_len = slen + 10000;
				samf7 = (char *)realloc(samf7, samf7_len);
				if(samf7 == NULL){
					fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
					exit(1);
				}
			}
			strcpy(samf7, "=");

			samf8 = (int)(p->pos - bns->anns[seqid].offset + 1);
			samf9 = 0;

		} else {

			slen = 1;
			if(slen >= samf7_len){
				samf7_len = slen + 10000;
				samf7 = (char *)realloc(samf7, samf7_len);
				if(samf7 == NULL){
					fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
					exit(1);
				}
			}
			strcpy(samf7, "=");

			samf8 = 0;
			samf9 = 0;

		}

		// print sequence and quality
		slen = p->full_len;
		if(slen >= samf10_len){
			samf10_len = slen + 10000;
			samf10 = (char *)realloc(samf10, samf10_len);
			if(samf10 == NULL){
				fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
				exit(1);
			}
		}
		if (p->strand == 0) {
			for (j = 0; j != p->full_len; ++j) {
				samf10[j] = "ACGTN"[(int)p->seq[j]];
			}
		} else {
			for (j = 0; j != p->full_len; ++j) {
				samf10[j] = "TGCAN"[p->seq[p->full_len - 1 - j]];
			}
		}
		samf10[j] = '\0';

		if (p->qual) {
			if (p->strand) {
				seq_reverse(p->len, p->qual, 0); // reverse quality
			}
			slen = p->len;
			if(slen >= samf11_len){
				samf11_len = slen + 10000;
				samf11 = (char *)realloc(samf11, samf11_len);
				if(samf11 == NULL){
					fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
					exit(1);
				}
			}
			strcpy(samf11, (char *)p->qual);
		} else {
			slen = 1;
			if(slen >= samf11_len){
				samf11_len = slen + 10000;
				samf11 = (char *)realloc(samf11, samf11_len);
				if(samf11 == NULL){
					fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
					exit(1);
				}
			}
			strcpy(samf11, "*");
		}

#if 0
		if (bwa_rg_id) {
			printf("\tRG:Z:%s", bwa_rg_id);
		}
		if (p->bc[0]) {
			printf("\tBC:Z:%s", p->bc);
		}
		if (p->clip_len < p->full_len) {
			printf("\tXC:i:%d", p->clip_len);
		}
#endif

		if (p->type != BWA_TYPE_NO_MATCH) {

#if 0
			// calculate XT tag
			char XT;
			XT = "NURM"[p->type];
			if (nn > 10) XT = 'N';
			// print tags
			printf("\tXT:A:%c\t%s:i:%d", XT, (mode & BWA_MODE_COMPREAD)? "NM" : "CM", p->nm);
			if (nn) printf("\tXN:i:%d", nn);
			if (mate) printf("\tSM:i:%d\tAM:i:%d", p->seQ, am);
			if (p->type != BWA_TYPE_MATESW) { // X0 and X1 are not available for this type of alignment
				printf("\tX0:i:%d", p->c1);
				if (p->c1 <= max_top2) printf("\tX1:i:%d", p->c2);
			}
			printf("\tXM:i:%d\tXO:i:%d\tXG:i:%d", p->n_mm, p->n_gapo, p->n_gapo+p->n_gape);
#endif

			if (p->md) {

				// printf("\tMD:Z:%s", p->md);

				slen = (int)strlen(p->md);
				if(slen >= sammdz_len){
					sammdz_len = slen + 10000;
					sam_mdzstr = (char *)realloc(sam_mdzstr, sammdz_len);
					if(sam_mdzstr == NULL){
						fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
						exit(1);
					}
				}
				strcpy(sam_mdzstr, p->md);

				j = 0;
				indx = 0;
				num_mut = 0;
				inside_del = 0;
				for(i=0; i<slen; i++){
					if(isdigit(p->md[i])){
						inside_del = 0;
						tmpstr[j] = p->md[i];
						j++;
						tmpstr[j] = '\0';
					}else{
						if(p->md[i] == '^') {
							inside_del = 1;
							indx += atoi(tmpstr);
							j = 0;
						}else if(!inside_del){
							indx += atoi(tmpstr);
							j = 0;
							int qual_val = 0;
							if (p->qual) {
								qual_val = (int)samf11[indx] - soap2_qual_adj;
							}
							sprintf(tmpstr2[num_mut],"\t%c->%d%c%d",p->md[i],indx,samf10[indx],qual_val);
							num_mut++;
							if(num_mut >= 1000){
								fprintf(stderr,"Error, too many mutations [1000]\n");
								exit(1);
							}
							indx++;
						}
					}
				}

				if(num_mut > 0){
					sprintf(tmpstr, "%d", num_mut);
					slen = (int)strlen(tmpstr);
					if(slen >= sammut_len){
						sammut_len = slen + 10000;
						sam_mutstr = (char *)realloc(sam_mutstr, sammut_len);
						if(sam_mutstr == NULL){
							fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
							exit(1);
						}
					}
					strcpy(sam_mutstr, tmpstr);

					for(i=0; i<num_mut; i++){
						slen += (int)strlen(tmpstr2[num_mut]);
						if(slen >= sammut_len){
							sammut_len = slen + 10000;
							sam_mutstr = (char *)realloc(sam_mutstr, sammut_len);
							if(sam_mutstr == NULL){
								fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
								exit(1);
							}
						}
						strcat(sam_mutstr, tmpstr2[i]);
					}
				}

				if(sam_cigar[0] != '\0'){
					slen = (int)strlen(sam_cigar) + 8;
					if(slen >= sammut_len){
						sammut_len = slen + 10000;
						sam_mutstr = (char *)realloc(sam_mutstr, sammut_len);
						if(sam_mutstr == NULL){
							fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
							exit(1);
						}
					}
					strcat(sam_mutstr, sam_cigar);
				}

			} else {

				sprintf(tmpstr, "%d", p->full_len);
				slen = (int)strlen(tmpstr);
				if(slen >= sammdz_len){
					sammdz_len = slen + 10000;
					sam_mdzstr = (char *)realloc(sam_mdzstr, sammdz_len);
					if(sam_mdzstr == NULL){
						fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
						exit(1);
					}
				}
				strcpy(sam_mdzstr, tmpstr);

			}

#if 0
			// print multiple hits
			if (p->n_multi) {
				printf("\tXA:Z:");
				for (i = 0; i < p->n_multi; ++i) {
					bwt_multi1_t *q = p->multi + i;
					int k;
					j = pos_end_multi(q, p->len) - q->pos;
					nn = bns_coor_pac2real(bns, q->pos, j, &seqid);
					printf("%s,%c%d,", bns->anns[seqid].name, q->strand? '-' : '+',
						   (int)(q->pos - bns->anns[seqid].offset + 1));
					if (q->cigar) {
						for (k = 0; k < q->n_cigar; ++k)
							printf("%d%c", __cigar_len(q->cigar[k]), "MIDS"[__cigar_op(q->cigar[k])]);
					} else printf("%dM", p->len);
					printf(",%d;", q->gap + q->mm);
				}
			}
#endif

		}

		// write out soap fields from sam fields here ...

		samf2a = 0;
		samf2b = 'c';
		samf2c = '+';
		if(samf2 & 0x10){
			samf2c = '-';
		}
		if(samf2 & 0x40){
			samf2a = 1;
			samf2b = 'a';
		}else if(samf2 & 0x80){
			samf2a = 2;
			samf2b = 'b';
		}
		samf5a = 0;
		if(samf5 >= 60){
			samf5a = 1;
		}

		if(sam_mdzstr[0] == '\0'){
			sprintf(tmpstr, "%d", p->full_len);
			slen = (int)strlen(tmpstr);
			if(slen >= sammdz_len){
				sammdz_len = slen + 10000;
				sam_mdzstr = (char *)realloc(sam_mdzstr, sammdz_len);
				if(sam_mdzstr == NULL){
					fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
					exit(1);
				}
			}
			strcpy(sam_mdzstr, tmpstr);
		}

		if(num_mut > 0){
			printf("%s/%d\t%s\t%s\t%d\t%c\t%d\t%c\t%s\t%d\t%s\t%dM\t%s\n",
				samf1,samf2a,samf10,samf11,samf5a,samf2b,p->full_len,samf2c,samf3,samf4,sam_mutstr,p->full_len,sam_mdzstr);
		}else{
			printf("%s/%d\t%s\t%s\t%d\t%c\t%d\t%c\t%s\t%d\t0\t%dM\t%s\n",
				samf1,samf2a,samf10,samf11,samf5a,samf2b,p->full_len,samf2c,samf3,samf4,p->full_len,sam_mdzstr);
		}

	} else { // this read has no match

		ubyte_t *s = p->strand? p->rseq : p->seq;
		int flag = p->extra_flag | SAM_FSU;
		if (mate && mate->type == BWA_TYPE_NO_MATCH) flag |= SAM_FMU;

		slen = (int)strlen(p->name);
		if(slen >= samf1_len){
			samf1_len = slen + 10000;
			samf1 = (char *)realloc(samf1, samf1_len);
			if(samf1 == NULL){
				fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
				exit(1);
			}
		}
		strcpy(samf1, p->name);

		samf2 = flag;

		slen = 1;
		if(slen >= samf3_len){
			samf3_len = slen + 10000;
			samf3 = (char *)realloc(samf3, samf3_len);
			if(samf3 == NULL){
				fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
				exit(1);
			}
		}
		strcpy(samf3, "*");

		samf4 = 0;
		samf5 = 0;

		// sprintf(tmp_cigar,"%d",p->len);

		slen = 1;
		if(slen >= samf6_len){
			samf6_len = slen + 10000;
			samf6 = (char *)realloc(samf6, samf6_len);
			if(samf6 == NULL){
				fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
				exit(1);
			}
		}
		strcpy(samf6, "*");

		slen = 1;
		if(slen >= samf7_len){
			samf7_len = slen + 10000;
			samf7 = (char *)realloc(samf7, samf7_len);
			if(samf7 == NULL){
				fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
				exit(1);
			}
		}
		strcpy(samf7, "*");

		samf8 = 0;
		samf9 = 0;

		slen = p->len;
		if(slen >= samf10_len){
			samf10_len = slen + 10000;
			samf10 = (char *)realloc(samf10, samf10_len);
			if(samf10 == NULL){
				fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
				exit(1);
			}
		}
		for (j = 0; j != p->len; ++j) {
			samf10[j] = "ACGTN"[(int)s[j]];
		}
		samf10[j] = '\0';

		if (p->qual) {
			if (p->strand) {
				seq_reverse(p->len, p->qual, 0); // reverse quality
			}
			slen = p->len;
			if(slen >= samf11_len){
				samf11_len = slen + 10000;
				samf11 = (char *)realloc(samf11, samf11_len);
				if(samf11 == NULL){
					fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
					exit(1);
				}
			}
			strcpy(samf11, (char *)p->qual);
		} else {
			slen = 1;
			if(slen >= samf11_len){
				samf11_len = slen + 10000;
				samf11 = (char *)realloc(samf11, samf11_len);
				if(samf11 == NULL){
					fprintf(stderr,"Error, unable to alloc mem for expanding output string\n");
					exit(1);
				}
			}
			strcpy(samf11, "*");
		}

#if 0
		if (bwa_rg_id) {
			printf("\tRG:Z:%s", bwa_rg_id);
		}
		if (p->bc[0]) {
			printf("\tBC:Z:%s", p->bc);
		}
		if (p->clip_len < p->full_len) {
			printf("\tXC:i:%d", p->clip_len);
		}
#endif

		// write out soap fields from sam fields here ...

		samf2a = 0;
		samf2b = 'c';
		samf2c = '+';
		if(samf2 & 0x10){
			samf2c = '-';
		}
		if(samf2 & 0x40){
			samf2a = 1;
			samf2b = 'a';
		}else if(samf2 & 0x80){
			samf2a = 2;
			samf2b = 'b';
		}
		samf5a = 0;
		if(samf5 >= 60){
			samf5a = 1;
		}
		
		printf("%s/%d\t%s\t%s\t%d\t%c\t%d\t%c\t%s\t%d\t0\t%dM\t%d\n",
			samf1,samf2a,samf10,samf11,samf5a,samf2b,p->full_len,samf2c,samf3,samf4,p->full_len,p->len);
	}

	return;
}

void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2)
{
	if(use_soap2_format){
		 bwa_print_soap1a(bns, p, mate, mode, max_top2);
	}else{
		 bwa_print_sam1a(bns, p, mate, mode, max_top2);
	}

	return;
}

// ------------------

bntseq_t *bwa_open_nt(const char *prefix)
{
	bntseq_t *ntbns;
	char *str;
	str = (char*)calloc(strlen(prefix) + 10, 1);
	strcat(strcpy(str, prefix), ".nt");
	ntbns = bns_restore(str);
	free(str);
	return ntbns;
}

void bwa_print_sam_SQ(const bntseq_t *bns)
{
	int i;

	if(use_soap2_format){
		return;
	}

	for (i = 0; i < bns->n_seqs; ++i)
		printf("@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
	if (bwa_rg_line) printf("%s\n", bwa_rg_line);
}

void bwase_initialize() 
{
	int i;
	for (i = 1; i != 256; ++i) g_log_n[i] = (int)(4.343 * log(i) + 0.5);
}

char *bwa_escape(char *s)
{
	char *p, *q;
	for (p = q = s; *p; ++p) {
		if (*p == '\\') {
			++p;
			if (*p == 't') *q++ = '\t';
			else if (*p == 'n') *q++ = '\n';
			else if (*p == 'r') *q++ = '\r';
			else if (*p == '\\') *q++ = '\\';
		} else *q++ = *p;
	}
	*q = '\0';
	return s;
}

int bwa_set_rg(const char *s)
{
	char *p, *q, *r;
	if (strstr(s, "@RG") != s) return -1;
	if (bwa_rg_line) free(bwa_rg_line);
	if (bwa_rg_id) free(bwa_rg_id);
	bwa_rg_line = strdup(s);
	bwa_rg_id = 0;
	bwa_escape(bwa_rg_line);
	p = strstr(bwa_rg_line, "\tID:");
	if (p == 0) return -1;
	p += 4;
	for (q = p; *q && *q != '\t' && *q != '\n'; ++q);
	bwa_rg_id = calloc(q - p + 1, 1);
	for (q = p, r = bwa_rg_id; *q && *q != '\t' && *q != '\n'; ++q)
		*r++ = *q;
	return 0;
}

void bwa_sai2sam_se_core(const char *prefix, const char *fn_sa, const char *fn_fa, int n_occ)
{
	int m_aln[3];
	int tot_seqs = 0;
	int n_seqs[3];
	bwt_aln1_t *aln[3];
	bwa_seq_t *seqs[3];
	bwa_seqio_t *ks;
	clock_t t;
	bntseq_t *bns, *ntbns = 0;
	FILE *fp_sa;
	gap_opt_t opt;
	int n_needed;
	int nexti1;
	int nexti2;
	int seqsd[3];
	int max_threads = 1;
	clock_t tio;
	int first = 1;

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

	fprintf(stderr, "[bwa_sai2sam_se_core] version: %s (%s)\n",
		bwaversionstr, bwablddatestr);
	fprintf(stderr, "[bwa_sai2sam_se_core] num threads: %d (max: %d)\n",
		num_sampe_threads, max_threads);

	// initialization
	bwase_initialize();
	bns = bns_restore(prefix);
	srand48(bns->seed);
	fp_sa = xopen(fn_sa, "r");

	fread(&opt, sizeof(gap_opt_t), 1, fp_sa);
	if (!(opt.mode & BWA_MODE_COMPREAD)) // in color space; initialize ntpac
		ntbns = bwa_open_nt(prefix);

	if(opt.mode & BWA_MODE_IL13){
		// fprintf(stderr,"[bwa_sai2sam_se_core] aln file in Illumina 1.3+ format\n")
		soap2_qual_adj = 64;
	}

	bwa_print_sam_SQ(bns);
	bwa_print_sam_PG();

	// set ks
	ks = bwa_open_reads(opt.mode, fn_fa);

	first = 1;
	seqsd[0] = 0;
	seqsd[1] = 0;
	seqsd[2] = 0;
	nexti1 = 0;
	nexti2 = 0;
	n_seqs[0] = 0;
	n_seqs[1] = 0;
	n_seqs[2] = 0;
	seqs[0] = NULL;
	seqs[1] = NULL;
	seqs[2] = NULL;
	aln[0] = NULL;
	aln[1] = NULL;
	aln[2] = NULL;
	m_aln[0] = 0;
	m_aln[1] = 0;
	m_aln[2] = 0;

	t = tio = clock();

	bwa_read_seq1_tpx(ks, n_needed, &n_seqs[nexti1], opt.mode, opt.trim_qual, &seqs[nexti1], &aln[nexti1], n_occ, fp_sa, m_aln[nexti1]);

	// core loop
	while(1){

		if( ( (async_read_seq) && (num_sampe_threads > 1) ) || (!first) ){
			tio = clock();
		}

		bwa_read_seq1_wait_tpx();

		if(seqs[nexti1] == NULL){
			break;
		}

		tot_seqs += n_seqs[nexti1];

		seqsd[nexti1] = 1;

		nexti2 = nexti1 + 1;
		if(nexti2 > 2){
			nexti2 = 0;
		}

		bwa_read_seq1_tpx(ks, n_needed, &n_seqs[nexti2], opt.mode, opt.trim_qual, &seqs[nexti2], &aln[nexti2], n_occ, fp_sa, m_aln[nexti2]);

		fprintf(stderr, "[bwa_sai2sam_se_core] bwa_read_seq1... %.2f sec", (float)(clock() - tio) / CLOCKS_PER_SEC);
		if( (async_read_seq) && (num_sampe_threads > 1) && (!first) ){
			fprintf(stderr," (async, bsize=%dk)\n", n_needed / 1024);
		}else{
			fprintf(stderr," (bsize=%dk)\n", n_needed / 1024);
		}

		first = 0;

#if 0
		// ---------------

		fprintf(stderr, "[bwa_sai2sam_se_core] read alignments... ");

		t = clock();

		// read alignment
		for (i = 0; i < n_seqs[nexti1]; ++i) {
			bwa_seq_t *p = seqs[nexti1] + i;
			int n_aln;

			fread(&n_aln, 4, 1, fp_sa);

			if (n_aln > m_aln) {
				m_aln = n_aln;
				aln = (bwt_aln1_t*)realloc(aln, sizeof(bwt_aln1_t) * m_aln);
			}

			fread(aln, sizeof(bwt_aln1_t), n_aln, fp_sa);

			bwa_aln2seq_core(n_aln, aln, p, 1, n_occ);
		}

		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); 
#endif
		// ---------------

		t = clock();

		fprintf(stderr, "[bwa_sai2sam_se_core] convert to sequence coordinate... ");
		bwa_cal_pac_pos(prefix, n_seqs[nexti1], seqs[nexti1], opt.max_diff, opt.fnr); // forward bwt will be destroyed here
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		// ---------------

		fprintf(stderr, "[bwa_sai2sam_se_core] refine gapped alignments... ");
		bwa_refine_gapped(bns, n_seqs[nexti1], seqs[nexti1], 0, ntbns);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		// ---------------

		fprintf(stderr, "[bwa_sai2sam_se_core] print alignments... ");
		bwa_print1_tpx(bns, n_seqs[nexti1], seqs[nexti1], opt);
		fprintf(stderr, "%.2f sec", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();
		if( (async_print_res) && (num_sampe_threads > 1) ){
			fprintf(stderr," (async)\n");
		}else{
			fprintf(stderr,"\n");
		}

		// ---------------

		fprintf(stderr, "[bwa_sai2sam_se_core] %d sequences have been processed.\n", tot_seqs);

		// ---------------

		nexti1 = nexti2;

		nexti2 = nexti1 + 1;
		if(nexti2 > 2){
			nexti2 = 0;
		}

		if(seqsd[nexti2]){
			bwa_free_read_seq(n_seqs[nexti2], seqs[nexti2]);
			seqsd[nexti2] = 0;
		}

	}

	// ---------------

	if( (async_print_res) && (num_sampe_threads > 1) ){
		t = clock();
		fprintf(stderr, "[bwa_sai2sam_se_core] wait for final print alignments... ");
	}

	bwa_print1_wait_tpx();

	if( (async_print_res) && (num_sampe_threads > 1) ){
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();
	}

	for(nexti1=0; nexti1<3; nexti1++){
		if(seqsd[nexti1]){
			bwa_free_read_seq(n_seqs[nexti1], seqs[nexti1]);
			seqsd[nexti1] = 0;
		}
		if(aln[nexti1]){
			free(aln[nexti1]);
			aln[nexti1] = 0;
		}
	}

	// ---------------

	// destroy
	bwa_seq_close(ks);
	if (ntbns) bns_destroy(ntbns);
	bns_destroy(bns);
	fclose(fp_sa);

	return;
}

int bwa_sai2sam_se(int argc, char *argv[])
{
	int c, n_occ = 3;
	struct timeval st;
	uint64_t s1, e1;
	double total_time = 0.0;

	gettimeofday(&st, NULL);
	s1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;

	while ((c = getopt(argc, argv, "ISTXYhn:f:r:t:")) >= 0) {
		switch (c) {
		case 't': num_sampe_threads = atoi(optarg); break;
		case 'h': break;
		case 'r':
			if (bwa_set_rg(optarg) < 0) {
				fprintf(stderr, "[%s] malformated @RG line\n", __func__);
				return 1;
			}
			break;
		case 'n': n_occ = atoi(optarg); break;
		case 'f': xreopen(optarg, "w", stdout); break;
		case 'T': adj_n_needed = 0; break;
		case 'X': async_read_seq = 0; break;
		case 'Y': async_print_res = 0; break;
		case 'S': use_soap2_format = 1; break;
		case 'I': soap2_qual_adj = 64; break;
		default: return 1;
		}
	}

	if (optind + 3 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: bwa samse [options] <prefix> <in.sai> <in.fq>\n\n");
		fprintf(stderr, "         -n       max_occ\n");
		fprintf(stderr, "         -r       RG_line\n");
		fprintf(stderr, "         -f FILE  sam file to output results to [stdout]\n");
		fprintf(stderr, "         -t INT   number of threads [%d] (use <=0 for all)\n", num_sampe_threads);
		fprintf(stderr, "         -I       adjust quality to Illumina v1.3+ format\n");
		fprintf(stderr, "         -S       output SOAP2 format\n");
		fprintf(stderr, "         -T       use original read buffer size\n");
		fprintf(stderr, "         -X       disable async read seq/aln method\n");
		fprintf(stderr, "         -Y       disable async print results method\n");
		fprintf(stderr, "\n");
		return 1;
	}

	bwa_sai2sam_se_core(argv[optind], argv[optind+1], argv[optind+2], n_occ);

	free(bwa_rg_line); free(bwa_rg_id);

	// cant use getrusage for ru.maxrss until kernel 2.6.36 ...
	int srtn = 0;
	long maxrsskb = 0L;
	srtn = getmaxrss(&maxrsskb);
	if(srtn == 0){
		fprintf(stderr,"[bwa_sai2sam_se] mem rss max = %ld (mb)\n",maxrsskb / 1024L);
	}

	gettimeofday(&st, NULL);
	e1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;
	total_time = (double)((double)e1 - (double)s1) / 1000000.0;

	fprintf(stderr,"[bwa_sai2sam_se] total time = %lf (sec)\n",total_time);

	return 0;
}
