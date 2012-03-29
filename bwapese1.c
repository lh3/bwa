#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#include "bwatpx.h"

extern bwa_cigar_t *refine_gapped_core(bwtint_t l_pac, const ubyte_t *pacseq, int len, const ubyte_t *seq, bwtint_t *_pos,
                                       int ext, int *n_cigar, int is_end_correct);
extern char *bwa_cal_md1(int n_cigar, bwa_cigar_t *cigar, int len, bwtint_t pos, ubyte_t *seq,
                         bwtint_t l_pac, ubyte_t *pacseq, kstring_t *str, int *_nm);
extern void bwa_correct_trimmed(bwa_seq_t *s);

#ifdef HAVE_PTHREAD
extern pthread_mutex_t pe_lock;
#endif // HAVE_PTHREAD

// -------------------

void bwa_rg_tpx(int iidx, const bntseq_t *bns, int n_seqs1, int n_seqs2, 
                bwa_seq_t *seqs, ubyte_t *pacseq, bntseq_t *ntbns)
{
	ubyte_t *ntpac = 0;
	int i, j;
	kstring_t *str;
#ifdef _TIMING
	struct timeval st;
	uint64_t s1, e1;
	double pos1_time = 0.0;
#endif

#ifdef _TIMING
	gettimeofday(&st, NULL);
	s1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;
#endif

	for (i = n_seqs1; i < n_seqs2; ++i) {
		bwa_seq_t *s = seqs + i;
		seq_reverse(s->len, s->seq, 0); // IMPORTANT: s->seq is reversed here!!!
		for (j = 0; j < s->n_multi; ++j) {
			bwt_multi1_t *q = s->multi + j;
			int n_cigar;
			if (q->gap == 0) continue;
			q->cigar = refine_gapped_core(bns->l_pac, pacseq, s->len, q->strand? s->rseq : s->seq, &q->pos,
										  (q->strand? 1 : -1) * q->gap, &n_cigar, 1);
			q->n_cigar = n_cigar;
		}
		if (s->type == BWA_TYPE_NO_MATCH || s->type == BWA_TYPE_MATESW || s->n_gapo == 0) continue;
		s->cigar = refine_gapped_core(bns->l_pac, pacseq, s->len, s->strand? s->rseq : s->seq, &s->pos,
									  (s->strand? 1 : -1) * (s->n_gapo + s->n_gape), &s->n_cigar, 1);
	}

	if (ntbns) { // in color space
		for (i = n_seqs1; i < n_seqs2; ++i) {
			bwa_seq_t *s = seqs + i;
			bwa_cs2nt_core(s, bns->l_pac, ntpac);
			for (j = 0; j < s->n_multi; ++j) {
				bwt_multi1_t *q = s->multi + j;
				int n_cigar;
				if (q->gap == 0) continue;
				free(q->cigar);
				q->cigar = refine_gapped_core(bns->l_pac, ntpac, s->len, q->strand? s->rseq : s->seq, &q->pos,
											  (q->strand? 1 : -1) * q->gap, &n_cigar, 0);
				q->n_cigar = n_cigar;
			}
			if (s->type != BWA_TYPE_NO_MATCH && s->cigar) { // update cigar again
				free(s->cigar);
				s->cigar = refine_gapped_core(bns->l_pac, ntpac, s->len, s->strand? s->rseq : s->seq, &s->pos,
											  (s->strand? 1 : -1) * (s->n_gapo + s->n_gape), &s->n_cigar, 0);
			}
		}
	}

	// generate MD tag
	str = (kstring_t*)calloc(1, sizeof(kstring_t));

	for (i = n_seqs1; i < n_seqs2; ++i) {
		bwa_seq_t *s = seqs + i;
		if (s->type != BWA_TYPE_NO_MATCH) {
			int nm;
			s->md = bwa_cal_md1(s->n_cigar, s->cigar, s->len, s->pos, s->strand? s->rseq : s->seq,
								bns->l_pac, ntbns? ntpac : pacseq, str, &nm);
			s->nm = nm;
		}
	}

        free(str->s); free(str);

	// correct for trimmed reads
	if (!ntbns) // trimming is only enabled for Illumina reads
		for (i = n_seqs1; i < n_seqs2; ++i) bwa_correct_trimmed(seqs + i);
       
#ifdef _TIMING
	gettimeofday(&st, NULL);
	e1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;
	pos1_time = (double)((double)e1 - (double)s1) / 1000000.0;

# ifdef HAVE_PTHREAD
	pthread_mutex_lock(&pe_lock);
# endif // HAVE_PTHREAD
	fprintf(stderr,"bwapese1 time = %lf (sec)\n",pos1_time);
# ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&pe_lock);
# endif // HAVE_PTHREAD
#endif
 
	return;
}

