#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#include "bwatpx.h"

extern uint64_t n_tot[], n_mapped[];

extern bwa_cigar_t *bwa_sw_core(bwtint_t l_pac, const ubyte_t *pacseq, int len, const ubyte_t *seq, 
                                int64_t *beg, int reglen, int *n_cigar, uint32_t *_cnt);

#ifdef HAVE_PTHREAD
extern pthread_mutex_t pe_lock;
#endif // HAVE_PTHREAD

// -------------------

void bwa_sw_tpx(int iidx, const bntseq_t *bns, const ubyte_t *pacseq, int n_seqs1, int n_seqs2, 
                bwa_seq_t *seqs[2], const pe_opt_t *popt, const isize_info_t *ii)
{
	int i;
#ifdef _TIMING
	struct timeval st;
	uint64_t s1, e1;
	double pos1_time = 0.0;
#endif

	// perform mate alignment
	n_tot[0] = n_tot[1] = n_mapped[0] = n_mapped[1] = 0;

#ifdef _DEBUG
# ifdef HAVE_PTHREAD
        pthread_mutex_lock(&pe_lock);
# endif // HAVE_PTHREAD
        fprintf(stderr,"bwape2: indx=%d start=%d end=%d\n",iidx,n_seqs1,n_seqs2);
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
		p[0] = seqs[0] + i; p[1] = seqs[1] + i;
		if ((p[0]->mapQ >= SW_MIN_MAPQ || p[1]->mapQ >= SW_MIN_MAPQ) && (p[0]->extra_flag&SAM_FPP) == 0) { // unpaired and one read has high mapQ
			int k, n_cigar[2], is_singleton, mapQ = 0, mq_adjust[2];
			int64_t beg[2], end[2];
			bwa_cigar_t *cigar[2];
			uint32_t cnt[2];

			/* In the following, _pref points to the reference read
			 * which must be aligned; _pmate points to its mate which is
			 * considered to be modified. */

#define __set_rght_coor(_a, _b, _pref, _pmate) do {						\
				(_a) = (int64_t)_pref->pos + ii->avg - 3 * ii->std - _pmate->len * 1.5; \
				(_b) = (_a) + 6 * ii->std + 2 * _pmate->len;			\
				if ((_a) < (int64_t)_pref->pos + _pref->len) (_a) = _pref->pos + _pref->len; \
				if ((_b) > bns->l_pac) (_b) = bns->l_pac;				\
			} while (0)

#define __set_left_coor(_a, _b, _pref, _pmate) do {						\
				(_a) = (int64_t)_pref->pos + _pref->len - ii->avg - 3 * ii->std - _pmate->len * 0.5; \
				(_b) = (_a) + 6 * ii->std + 2 * _pmate->len;			\
				if ((_a) < 0) (_a) = 0;									\
				if ((_b) > _pref->pos) (_b) = _pref->pos;				\
			} while (0)
			
#define __set_fixed(_pref, _pmate, _beg, _cnt) do {						\
				_pmate->type = BWA_TYPE_MATESW;							\
				_pmate->pos = _beg;										\
				_pmate->seQ = _pref->seQ;								\
				_pmate->strand = (popt->type == BWA_PET_STD)? 1 - _pref->strand : _pref->strand; \
				_pmate->n_mm = _cnt>>16; _pmate->n_gapo = _cnt>>8&0xff; _pmate->n_gape = _cnt&0xff; \
				_pmate->extra_flag |= SAM_FPP;							\
				_pref->extra_flag |= SAM_FPP;							\
			} while (0)

			mq_adjust[0] = mq_adjust[1] = 255; // not effective
			is_singleton = (p[0]->type == BWA_TYPE_NO_MATCH || p[1]->type == BWA_TYPE_NO_MATCH)? 1 : 0;

			++n_tot[is_singleton];
			cigar[0] = cigar[1] = 0;
			n_cigar[0] = n_cigar[1] = 0;
			if (popt->type != BWA_PET_STD && popt->type != BWA_PET_SOLID) continue; // other types of pairing is not considered
			for (k = 0; k < 2; ++k) { // p[1-k] is the reference read and p[k] is the read considered to be modified
				ubyte_t *seq;
				if (p[1-k]->type == BWA_TYPE_NO_MATCH) continue; // if p[1-k] is unmapped, skip
				if (popt->type == BWA_PET_STD) {
					if (p[1-k]->strand == 0) { // then the mate is on the reverse strand and has larger coordinate
						__set_rght_coor(beg[k], end[k], p[1-k], p[k]);
						seq = p[k]->rseq;
					} else { // then the mate is on forward stand and has smaller coordinate
						__set_left_coor(beg[k], end[k], p[1-k], p[k]);
						seq = p[k]->seq;
						seq_reverse(p[k]->len, seq, 0); // because ->seq is reversed; this will reversed back shortly
					}
				} else { // BWA_PET_SOLID
					if (p[1-k]->strand == 0) { // R3-F3 pairing
						if (k == 0) __set_left_coor(beg[k], end[k], p[1-k], p[k]); // p[k] is R3
						else __set_rght_coor(beg[k], end[k], p[1-k], p[k]); // p[k] is F3
						seq = p[k]->rseq;
						seq_reverse(p[k]->len, seq, 0); // because ->seq is reversed
					} else { // F3-R3 pairing
						if (k == 0) __set_rght_coor(beg[k], end[k], p[1-k], p[k]); // p[k] is R3
						else __set_left_coor(beg[k], end[k], p[1-k], p[k]); // p[k] is F3
						seq = p[k]->seq;
					}
				}
				// perform SW alignment
				cigar[k] = bwa_sw_core(bns->l_pac, pacseq, p[k]->len, seq, &beg[k], end[k] - beg[k], &n_cigar[k], &cnt[k]);
				if (cigar[k] && p[k]->type != BWA_TYPE_NO_MATCH) { // re-evaluate cigar[k]
					int s_old, clip = 0, s_new;
					if (__cigar_op(cigar[k][0]) == 3) clip += __cigar_len(cigar[k][0]);
					if (__cigar_op(cigar[k][n_cigar[k]-1]) == 3) clip += __cigar_len(cigar[k][n_cigar[k]-1]);
					s_old = (int)((p[k]->n_mm * 9 + p[k]->n_gapo * 13 + p[k]->n_gape * 2) / 3. * 8. + .499);
					s_new = (int)(((cnt[k]>>16) * 9 + (cnt[k]>>8&0xff) * 13 + (cnt[k]&0xff) * 2 + clip * 3) / 3. * 8. + .499);
					s_old += -4.343 * log(ii->ap_prior / bns->l_pac);
					s_new += (int)(-4.343 * log(.5 * erfc(M_SQRT1_2 * 1.5) + .499)); // assume the mapped isize is 1.5\sigma
					if (s_old < s_new) { // reject SW alignment
						mq_adjust[k] = s_new - s_old;
						free(cigar[k]); cigar[k] = 0; n_cigar[k] = 0;
					} else mq_adjust[k] = s_old - s_new;
				}
				// now revserse sequence back such that p[*]->seq looks untouched
				if (popt->type == BWA_PET_STD) {
					if (p[1-k]->strand == 1) seq_reverse(p[k]->len, seq, 0);
				} else {
					if (p[1-k]->strand == 0) seq_reverse(p[k]->len, seq, 0);
				}
			}
			k = -1; // no read to be changed
			if (cigar[0] && cigar[1]) {
				k = p[0]->mapQ < p[1]->mapQ? 0 : 1; // p[k] to be fixed
				mapQ = abs(p[1]->mapQ - p[0]->mapQ);
			} else if (cigar[0]) k = 0, mapQ = p[1]->mapQ;
			else if (cigar[1]) k = 1, mapQ = p[0]->mapQ;
			if (k >= 0 && p[k]->pos != beg[k]) {
				++n_mapped[is_singleton];
				{ // recalculate mapping quality
					int tmp = (int)p[1-k]->mapQ - p[k]->mapQ/2 - 8;
					if (tmp <= 0) tmp = 1;
					if (mapQ > tmp) mapQ = tmp;
					p[k]->mapQ = p[1-k]->mapQ = mapQ;
					p[k]->seQ = p[1-k]->seQ = p[1-k]->seQ < mapQ? p[1-k]->seQ : mapQ;
					if (p[k]->mapQ > mq_adjust[k]) p[k]->mapQ = mq_adjust[k];
					if (p[k]->seQ > mq_adjust[k]) p[k]->seQ = mq_adjust[k];
				}
				// update CIGAR
				free(p[k]->cigar); p[k]->cigar = cigar[k]; cigar[k] = 0;
				p[k]->n_cigar = n_cigar[k];
				// update the rest of information
				__set_fixed(p[1-k], p[k], beg[k], cnt[k]);
			}
			free(cigar[0]); free(cigar[1]);
		}
	}

#ifdef _TIMING
	gettimeofday(&st, NULL);
	e1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;
        pos1_time = (double)((double)e1 - (double)s1) / 1000000.0;

# ifdef HAVE_PTHREAD
        pthread_mutex_lock(&pe_lock);
# endif // HAVE_PTHREAD
        fprintf(stderr,"bwape2 time = %lf (sec)\n",pos1_time);
# ifdef HAVE_PTHREAD
        pthread_mutex_unlock(&pe_lock);
# endif // HAVE_PTHREAD
#endif

	return;
}

