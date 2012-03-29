#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#include "bwatpx.h"

#include "khash.h"
KHASH_MAP_INIT_INT64(64, poslist_t)

#include "ksort.h"
KSORT_INIT_GENERIC(uint64_t)

extern void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi);

extern kh_64_t *g_hash[]; // in bwape.c

extern int g_log_n[]; // in bwase.c

#ifdef HAVE_PTHREAD
pthread_mutex_t pe_lock = PTHREAD_MUTEX_INITIALIZER;
#endif // HAVE_PTHREAD

// -------------------

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

// -------------------

static int pairing(int iidx, bwa_seq_t *p[2], pe_data_t *d[MAX_CPUS], const pe_opt_t *opt, int s_mm, const isize_info_t *ii)
{
	int i, j, o_n, subo_n, cnt_chg = 0, low_bound = ii->low, max_len;
	uint64_t last_pos[2][2], o_pos[2], subo_score, o_score;
	max_len = p[0]->full_len;
	if (max_len < p[1]->full_len) max_len = p[1]->full_len;
	if (low_bound < max_len) low_bound = max_len;

	// here v>=u. When ii is set, we check insert size with ii; otherwise with opt->max_isize
#define __pairing_aux(u,v) do {											\
		bwtint_t l = ((v)>>32) + p[(v)&1]->len - ((u)>>32);				\
		if ((u) != (uint64_t)-1 && (v)>>32 > (u)>>32 && l >= max_len	\
			&& ((ii->high && l <= ii->high_bayesian) || (ii->high == 0 && l <= opt->max_isize))) \
		{																\
			uint64_t s = d[iidx]->aln[(v)&1].a[(uint32_t)(v)>>1].score + d[iidx]->aln[(u)&1].a[(uint32_t)(u)>>1].score; \
			s *= 10;													\
			if (ii->high) s += (int)(-4.343 * log(.5 * erfc(M_SQRT1_2 * fabs(l - ii->avg) / ii->std)) + .499); \
			s = s<<32 | (uint32_t)hash_64((u)>>32<<32 | (v)>>32);		\
			if (s>>32 == o_score>>32) ++o_n;							\
			else if (s>>32 < o_score>>32) { subo_n += o_n; o_n = 1; }	\
			else ++subo_n;												\
			if (s < o_score) subo_score = o_score, o_score = s, o_pos[(u)&1] = (u), o_pos[(v)&1] = (v); \
			else if (s < subo_score) subo_score = s;					\
		}																\
	} while (0)

#define __pairing_aux2(q, w) do {										\
		const bwt_aln1_t *r = d[iidx]->aln[(w)&1].a + ((uint32_t)(w)>>1);		\
		(q)->extra_flag |= SAM_FPP;										\
		if ((q)->pos != (w)>>32 || (q)->strand != r->a) {				\
			(q)->n_mm = r->n_mm; (q)->n_gapo = r->n_gapo; (q)->n_gape = r->n_gape; (q)->strand = r->a; \
			(q)->score = r->score;										\
			(q)->pos = (w)>>32;											\
			if ((q)->mapQ > 0) ++cnt_chg;								\
		}																\
	} while (0)

	o_score = subo_score = (uint64_t)-1;
	o_n = subo_n = 0;
	ks_introsort(uint64_t, d[iidx]->arr.n, d[iidx]->arr.a);
	for (j = 0; j < 2; ++j) last_pos[j][0] = last_pos[j][1] = (uint64_t)-1;
	if (opt->type == BWA_PET_STD) {
		for (i = 0; i < d[iidx]->arr.n; ++i) {
			uint64_t x = d[iidx]->arr.a[i];
			int strand = d[iidx]->aln[x&1].a[(uint32_t)x>>1].a;
			if (strand == 1) { // reverse strand, then check
				int y = 1 - (x&1);
				__pairing_aux(last_pos[y][1], x);
				__pairing_aux(last_pos[y][0], x);
			} else { // forward strand, then push
				last_pos[x&1][0] = last_pos[x&1][1];
				last_pos[x&1][1] = x;
			}
		}
	} else if (opt->type == BWA_PET_SOLID) {
		for (i = 0; i < d[iidx]->arr.n; ++i) {
			uint64_t x = d[iidx]->arr.a[i];
			int strand = d[iidx]->aln[x&1].a[(uint32_t)x>>1].a;
			if ((strand^x)&1) { // push
				int y = 1 - (x&1);
				__pairing_aux(last_pos[y][1], x);
				__pairing_aux(last_pos[y][0], x);
			} else { // check
				last_pos[x&1][0] = last_pos[x&1][1];
				last_pos[x&1][1] = x;
			}
		}
	} else {
		fprintf(stderr, "[paring] not implemented yet!\n");
		exit(1);
	}

	// set pairing
	//fprintf(stderr, "[%d, %d, %d, %d]\n", d[iidx]->arr.n, (int)(o_score>>32), (int)(subo_score>>32), o_n);

	if (o_score != (uint64_t)-1) {
		int mapQ_p = 0; // this is the maximum mapping quality when one end is moved
		int rr[2];
		//fprintf(stderr, "%d, %d\n", o_n, subo_n);
		if (o_n == 1) {
			if (subo_score == (uint64_t)-1) mapQ_p = 29; // no sub-optimal pair
			else if ((subo_score>>32) - (o_score>>32) > s_mm * 10) mapQ_p = 23; // poor sub-optimal pair
			else {
				int n = subo_n > 255? 255 : subo_n;
				mapQ_p = ((subo_score>>32) - (o_score>>32)) / 2 - g_log_n[n];
				if (mapQ_p < 0) mapQ_p = 0;
			}
		}
		rr[0] = d[iidx]->aln[o_pos[0]&1].a[(uint32_t)o_pos[0]>>1].a;
		rr[1] = d[iidx]->aln[o_pos[1]&1].a[(uint32_t)o_pos[1]>>1].a;
		if ((p[0]->pos == o_pos[0]>>32 && p[0]->strand == rr[0]) && (p[1]->pos == o_pos[1]>>32 && p[1]->strand == rr[1])) { // both ends not moved
			if (p[0]->mapQ > 0 && p[1]->mapQ > 0) {
				int mapQ = p[0]->mapQ + p[1]->mapQ;
				if (mapQ > 60) mapQ = 60;
				p[0]->mapQ = p[1]->mapQ = mapQ;
			} else {
				if (p[0]->mapQ == 0) p[0]->mapQ = (mapQ_p + 7 < p[1]->mapQ)? mapQ_p + 7 : p[1]->mapQ;
				if (p[1]->mapQ == 0) p[1]->mapQ = (mapQ_p + 7 < p[0]->mapQ)? mapQ_p + 7 : p[0]->mapQ;
			}
		} else if (p[0]->pos == o_pos[0]>>32 && p[0]->strand == rr[0]) { // [1] moved
			p[1]->seQ = 0; p[1]->mapQ = p[0]->mapQ;
			if (p[1]->mapQ > mapQ_p) p[1]->mapQ = mapQ_p;
		} else if (p[1]->pos == o_pos[1]>>32 && p[1]->strand == rr[1]) { // [0] moved
			p[0]->seQ = 0; p[0]->mapQ = p[1]->mapQ;
			if (p[0]->mapQ > mapQ_p) p[0]->mapQ = mapQ_p;
		} else { // both ends moved
			p[0]->seQ = p[1]->seQ = 0;
			mapQ_p -= 20;
			if (mapQ_p < 0) mapQ_p = 0;
			p[0]->mapQ = p[1]->mapQ = mapQ_p;
		}
		__pairing_aux2(p[0], o_pos[0]);
		__pairing_aux2(p[1], o_pos[1]);
	}
	return cnt_chg;
}

// -------------------

int bwa_pe_tpx(int iidx, const bwt_t *bwt[2], int n_seqs1, int n_seqs2, bwa_seq_t *seqs[2], isize_info_t *ii,
               const pe_opt_t *opt, const gap_opt_t *gopt, pe_data_t *d[MAX_CPUS], aln_buf_t *buf[2])
{
	int i, j, cnt_chg = 0;
#ifdef _TIMING
	struct timeval st;
	uint64_t s1, e1;
	double pos1_time = 0.0;
#endif

#ifdef _DEBUG
# ifdef HAVE_PTHREAD
	pthread_mutex_lock(&pe_lock);
# endif // HAVE_PTHREAD
	fprintf(stderr,"bwape1: indx=%d start=%d end=%d\n",iidx,n_seqs1,n_seqs2);
# ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&pe_lock);
# endif // HAVE_PTHREAD
#endif

#ifdef _TIMING
	gettimeofday(&st, NULL);
	s1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;
#endif

#ifdef _USE_LOCAL_GHASH
	g_hash[iidx] = kh_init(64);
#endif // _USE_LOCAL_GHASH

	// PE
	for (i = n_seqs1; i < n_seqs2; ++i) {

		bwa_seq_t *p[2];

		for (j = 0; j < 2; ++j) {
			p[j] = seqs[j] + i;
			kv_copy(bwt_aln1_t, d[iidx]->aln[j], buf[j][i].aln);
		}

		if ((p[0]->type == BWA_TYPE_UNIQUE || p[0]->type == BWA_TYPE_REPEAT)
			&& (p[1]->type == BWA_TYPE_UNIQUE || p[1]->type == BWA_TYPE_REPEAT))
		{ // only when both ends mapped

			uint64_t x;
			int k;
			long n_occ[2];

			for (j = 0; j < 2; ++j) {
				n_occ[j] = 0;
				for (k = 0; k < d[iidx]->aln[j].n; ++k)
					n_occ[j] += d[iidx]->aln[j].a[k].l - d[iidx]->aln[j].a[k].k + 1;
			}

			if (n_occ[0] > opt->max_occ || n_occ[1] > opt->max_occ) continue;

			d[iidx]->arr.n = 0;
			for (j = 0; j < 2; ++j) {
				for (k = 0; k < d[iidx]->aln[j].n; ++k) {

					bwt_aln1_t *r = d[iidx]->aln[j].a + k;
					bwtint_t l;

					if (r->l - r->k + 1 >= MIN_HASH_WIDTH) { // then check hash table
						uint64_t key = (uint64_t)r->k<<32 | r->l;
						int ret;
						khint_t iter = kh_put(64, g_hash[iidx], key, &ret);
						if (ret) { // not in the hash table; ret must equal 1 as we never remove elements
							poslist_t *z = &kh_val(g_hash[iidx], iter);
							z->n = r->l - r->k + 1;
							z->a = (bwtint_t*)malloc(sizeof(bwtint_t) * z->n);
							for (l = r->k; l <= r->l; ++l)
								z->a[l - r->k] = r->a? bwt_sa(bwt[0], l) : bwt[1]->seq_len - (bwt_sa(bwt[1], l) + p[j]->len);
						}
						for (l = 0; l < kh_val(g_hash[iidx], iter).n; ++l) {
							x = kh_val(g_hash[iidx], iter).a[l];
							x = x<<32 | k<<1 | j;
							kv_push(uint64_t, d[iidx]->arr, x);
						}
					} else { // then calculate on the fly
						for (l = r->k; l <= r->l; ++l) {
							x = r->a? bwt_sa(bwt[0], l) : bwt[1]->seq_len - (bwt_sa(bwt[1], l) + p[j]->len);
							x = x<<32 | k<<1 | j;
							kv_push(uint64_t, d[iidx]->arr, x);
						}
					}
				}
			}

			cnt_chg += pairing(iidx, p, d, opt, gopt->s_mm, ii);

		}

		if (opt->N_multi || opt->n_multi) {
			for (j = 0; j < 2; ++j) {
				if (p[j]->type != BWA_TYPE_NO_MATCH) {

					int k;

					if (!(p[j]->extra_flag&SAM_FPP) && p[1-j]->type != BWA_TYPE_NO_MATCH) {
						bwa_aln2seq_core(d[iidx]->aln[j].n, d[iidx]->aln[j].a, p[j], 0, 
									(p[j]->c1+p[j]->c2-1 > opt->N_multi? opt->n_multi : opt->N_multi));
					} else {
						bwa_aln2seq_core(d[iidx]->aln[j].n, d[iidx]->aln[j].a, p[j], 0, opt->n_multi);
					}

					for (k = 0; k < p[j]->n_multi; ++k) {
						bwt_multi1_t *q = p[j]->multi + k;
						q->pos = q->strand? bwt_sa(bwt[0], q->pos) : bwt[1]->seq_len - (bwt_sa(bwt[1], q->pos) + p[j]->len);
					}
				}
			}
		}

	}

#ifdef _USE_LOCAL_GHASH
	khint_t iter;
	for (iter = kh_begin(g_hash[iidx]); iter != kh_end(g_hash[iidx]); ++iter){
		if (kh_exist(g_hash[iidx], iter)) free(kh_val(g_hash[iidx], iter).a);
	}

	kh_destroy(64, g_hash[iidx]);
#endif // _USE_LOCAL_GHASH

#ifdef _TIMING
	gettimeofday(&st, NULL);
	e1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;
	pos1_time = (double)((double)e1 - (double)s1) / 1000000.0;

# ifdef HAVE_PTHREAD
	pthread_mutex_lock(&pe_lock);
# endif // HAVE_PTHREAD
	fprintf(stderr,"bwape1 time = %lf (sec)\n",pos1_time);
# ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&pe_lock);
# endif // HAVE_PTHREAD
#endif

	return cnt_chg;
}

