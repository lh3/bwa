#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#include "bwatpx.h"

extern int bwa_approx_mapQ(const bwa_seq_t *p, int mm);

#ifdef HAVE_PTHREAD
extern pthread_mutex_t pe_lock;
#endif // HAVE_PTHREAD

// -------------------

void bwa_se_tpx(int iidx, bwt_t *bwt[2], int n_seqs1, int n_seqs2, bwa_seq_t *seqs[2], const gap_opt_t *gopt)
{
	int i, j;
#ifdef _TIMING
	struct timeval st;
	uint64_t s1, e1;
	double pos1_time = 0.0;
#endif

#ifdef _DEBUG
# ifdef HAVE_PTHREAD
        pthread_mutex_lock(&pe_lock);
# endif // HAVE_PTHREAD
        fprintf(stderr,"bwape3: indx=%d start=%d end=%d\n",iidx,n_seqs1,n_seqs2);
# ifdef HAVE_PTHREAD
        pthread_mutex_unlock(&pe_lock);
# endif // HAVE_PTHREAD
#endif

#ifdef _TIMING
	gettimeofday(&st, NULL);
	s1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;
#endif

	for (i = n_seqs1; i < n_seqs2; ++i) {

		bwa_seq_t *p[2];

		for (j = 0; j < 2; ++j) {

			p[j] = seqs[j] + i;

			if (p[j]->type == BWA_TYPE_UNIQUE || p[j]->type == BWA_TYPE_REPEAT) {

				int max_diff = gopt->fnr > 0.0? bwa_cal_maxdiff(p[j]->len, BWA_AVG_ERR, gopt->fnr) : gopt->max_diff;

				// tpx - all the SE time is here ...

				p[j]->pos = p[j]->strand? bwt_sa(bwt[0], p[j]->sa)
						: bwt[1]->seq_len - (bwt_sa(bwt[1], p[j]->sa) + p[j]->len);

				// tpx ------

				p[j]->seQ = p[j]->mapQ = bwa_approx_mapQ(p[j], max_diff);

			}

		}

	}

#ifdef _TIMING
	gettimeofday(&st, NULL);
	e1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;
        pos1_time = (double)((double)e1 - (double)s1) / 1000000.0;

# ifdef HAVE_PTHREAD
        pthread_mutex_lock(&pe_lock);
# endif // HAVE_PTHREAD
        fprintf(stderr,"bwape3 time = %lf (sec)\n",pos1_time);
# ifdef HAVE_PTHREAD
        pthread_mutex_unlock(&pe_lock);
# endif // HAVE_PTHREAD
#endif

	return;
}

