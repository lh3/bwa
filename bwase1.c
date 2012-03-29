#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <sys/time.h>
#include <time.h>
#include "bwatpx.h"

extern int bwa_approx_mapQ(const bwa_seq_t *p, int mm);
extern bwtint_t bwa_sa2pos(const bntseq_t *bns, const bwt_t *bwt, bwtint_t sapos, int len, int *strand);

// -------------------

void bwa_cal_pac_pos2_tpx(const bntseq_t *bns, const bwt_t *bwt, int n_seqs1, int n_seqs2,
                          bwa_seq_t *seqs, const int max_mm, const float fnr)
{
	int i, j;
	int strand;
	int max_diff;
	int n_multi;

	for (i = n_seqs1; i < n_seqs2; ++i) {

		if (seqs[i].type == BWA_TYPE_UNIQUE || seqs[i].type == BWA_TYPE_REPEAT) {
			max_diff = fnr > 0.0 ? bwa_cal_maxdiff(seqs[i].len, BWA_AVG_ERR, fnr) : max_mm;
			seqs[i].pos = bwa_sa2pos(bns, bwt, seqs[i].sa, seqs[i].len, &strand);
			seqs[i].strand = strand;
			seqs[i].seQ = seqs[i].mapQ = bwa_approx_mapQ(&seqs[i], max_diff);
		}

		n_multi = 0;
                for (j = 0; j < seqs[i].n_multi; ++j) {
                        bwt_multi1_t *q = seqs[i].multi + j;
                        q->pos = bwa_sa2pos(bns, bwt, q->pos, seqs[i].len, &strand);
                        q->strand = strand;
			if (q->pos != seqs[i].pos)
				seqs[i].multi[n_multi++] = *q;
                }
		seqs[i].n_multi = n_multi;

        }

	return;
}

