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

void bwa_cal_pac_posf2_tpx(const bwt_t *forward_bwt,int n_seqs1,int n_seqs2,bwa_seq_t *seqs,const int max_mm, const float fnr)
{
	int i, j;
	int max_diff;

	for (i = n_seqs1; i < n_seqs2; ++i) {

                if (seqs[i].strand) {
			if (seqs[i].type == BWA_TYPE_UNIQUE || seqs[i].type == BWA_TYPE_REPEAT) {
				max_diff = fnr > 0.0? bwa_cal_maxdiff(seqs[i].len, BWA_AVG_ERR, fnr) : max_mm;
				seqs[i].pos = bwt_sa(forward_bwt, seqs[i].sa);
				seqs[i].seQ = seqs[i].mapQ = bwa_approx_mapQ(&seqs[i], max_diff);
			}
		}

                for (j = 0; j < seqs[i].n_multi; ++j) {
                        bwt_multi1_t *p = seqs[i].multi + j;
                        if (p->strand) p->pos = bwt_sa(forward_bwt, p->pos);
                }

        }

	return;
}

void bwa_cal_pac_posr2_tpx(const bwt_t *reverse_bwt,int n_seqs1,int n_seqs2,bwa_seq_t *seqs,const int max_mm, const float fnr)
{
	int i, j;
	int max_diff;

	for (i = n_seqs1; i < n_seqs2; ++i) {

                if (!seqs[i].strand) {
			if (seqs[i].type == BWA_TYPE_UNIQUE || seqs[i].type == BWA_TYPE_REPEAT) {
				max_diff = fnr > 0.0? bwa_cal_maxdiff(seqs[i].len, BWA_AVG_ERR, fnr) : max_mm;
				seqs[i].pos = reverse_bwt->seq_len - (bwt_sa(reverse_bwt, seqs[i].sa) + seqs[i].len);
				seqs[i].seQ = seqs[i].mapQ = bwa_approx_mapQ(&seqs[i], max_diff);
			}
		}

                for (j = 0; j < seqs[i].n_multi; ++j) {
                        bwt_multi1_t *p = seqs[i].multi + j;
                        if (!p->strand) p->pos = reverse_bwt->seq_len - (bwt_sa(reverse_bwt, p->pos) + seqs[i].len);
                }

        }

	return;
}

