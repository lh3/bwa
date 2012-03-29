#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include "bwatpx.h"

#include "khash.h"
KHASH_INIT(b128, b128_t, poslist_t, 1, b128_hash, b128_eq)

#include "ksort.h"
KSORT_INIT(b128, b128_t, b128_lt)

extern void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi);
extern bwtint_t bwa_sa2pos(const bntseq_t *bns, const bwt_t *bwt, bwtint_t sapos, int len, int *strand);

extern kh_b128_t *g_hash[]; // in bwape.c

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

static int pairing(bwa_seq_t *p[2], pe_data_t *d, const pe_opt_t *opt, int s_mm, const isize_info_t *ii)
{
	int i, j, o_n, subo_n, cnt_chg = 0, low_bound = ii->low, max_len;
	uint64_t o_score, subo_score;
	b128_t last_pos[2][2], o_pos[2];
	max_len = p[0]->full_len;
	if (max_len < p[1]->full_len) max_len = p[1]->full_len;
	if (low_bound < max_len) low_bound = max_len;

	// here v>=u. When ii is set, we check insert size with ii; otherwise with opt->max_isize
#define __pairing_aux(u,v) do { \
		bwtint_t l = (v).x + p[(v).y&1]->len - ((u).x); \
		if ((u).x != (uint64_t)-1 && (v).x > (u).x && l >= max_len \
			&& ((ii->high && l <= ii->high_bayesian) || (ii->high == 0 && l <= opt->max_isize))) \
		{ \
			uint64_t s = d->aln[(v).y&1].a[(v).y>>2].score + d->aln[(u).y&1].a[(u).y>>2].score; \
			s *= 10; \
			if (ii->high) s += (int)(-4.343 * log(.5 * erfc(M_SQRT1_2 * fabs(l - ii->avg) / ii->std)) + .499); \
			s = s<<32 | (uint32_t)hash_64((u).x<<32 | (v).x); \
			if (s>>32 == o_score>>32) ++o_n; \
			else if (s>>32 < o_score>>32) { subo_n += o_n; o_n = 1; } \
			else ++subo_n; \
			if (s < o_score) subo_score = o_score, o_score = s, o_pos[(u).y&1] = (u), o_pos[(v).y&1] = (v); \
			else if (s < subo_score) subo_score = s; \
		} \
	} while (0)

#define __pairing_aux2(q, w) do { \
		const bwt_aln1_t *r = d->aln[(w).y&1].a + ((w).y>>2); \
		(q)->extra_flag |= SAM_FPP; \
		if ((q)->pos != (w).x || (q)->strand != ((w).y>>1&1)) { \
			(q)->n_mm = r->n_mm; (q)->n_gapo = r->n_gapo; (q)->n_gape = r->n_gape; (q)->strand = (w).y>>1&1; \
			(q)->score = r->score; \
			(q)->pos = (w).x; \
			if ((q)->mapQ > 0) ++cnt_chg; \
		} \
	} while (0)

	o_score = subo_score = (uint64_t)-1;
	o_n = subo_n = 0;
	ks_introsort(b128, d->arr.n, d->arr.a);
	for (j = 0; j < 2; ++j) last_pos[j][0].x = last_pos[j][0].y = last_pos[j][1].x = last_pos[j][1].y = (uint64_t)-1;
	if (opt->type == BWA_PET_STD) {
		for (i = 0; i < d->arr.n; ++i) {
			b128_t x = d->arr.a[i];
			int strand = x.y>>1&1;
			if (strand == 1) { // reverse strand, then check
				int y = 1 - (x.y&1);
				__pairing_aux(last_pos[y][1], x);
				__pairing_aux(last_pos[y][0], x);
			} else { // forward strand, then push
				last_pos[x.y&1][0] = last_pos[x.y&1][1];
				last_pos[x.y&1][1] = x;
			}
		}
	} else if (opt->type == BWA_PET_SOLID) {
		for (i = 0; i < d->arr.n; ++i) {
			b128_t x = d->arr.a[i];
			int strand = x.y>>1&1;
			if ((strand^x.y)&1) { // push
				int y = 1 - (x.y&1);
				__pairing_aux(last_pos[y][1], x);
				__pairing_aux(last_pos[y][0], x);
			} else { // check
				last_pos[x.y&1][0] = last_pos[x.y&1][1];
				last_pos[x.y&1][1] = x;
			}
		}
	} else {
		fprintf(stderr, "[paring] not implemented yet!\n");
		exit(1);
	}
	// set pairing
	//fprintf(stderr, "[%ld, %d, %d, %d]\n", d->arr.n, (int)(o_score>>32), (int)(subo_score>>32), o_n);
	if (o_score != (uint64_t)-1) {
		int mapQ_p = 0; // this is the maximum mapping quality when one end is moved
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
		if ((p[0]->pos == o_pos[0].x && p[0]->strand == (o_pos[0].y>>1&1)) && (p[1]->pos == o_pos[1].x && p[1]->strand == (o_pos[1].y>>1&1))) { // both ends not moved
			if (p[0]->mapQ > 0 && p[1]->mapQ > 0) {
				int mapQ = p[0]->mapQ + p[1]->mapQ;
				if (mapQ > 60) mapQ = 60;
				p[0]->mapQ = p[1]->mapQ = mapQ;
			} else {
				if (p[0]->mapQ == 0) p[0]->mapQ = (mapQ_p + 7 < p[1]->mapQ)? mapQ_p + 7 : p[1]->mapQ;
				if (p[1]->mapQ == 0) p[1]->mapQ = (mapQ_p + 7 < p[0]->mapQ)? mapQ_p + 7 : p[0]->mapQ;
			}
		} else if (p[0]->pos == o_pos[0].x && p[0]->strand == (o_pos[0].y>>1&1)) { // [1] moved
			p[1]->seQ = 0; p[1]->mapQ = p[0]->mapQ;
			if (p[1]->mapQ > mapQ_p) p[1]->mapQ = mapQ_p;
		} else if (p[1]->pos == o_pos[1].x && p[1]->strand == (o_pos[1].y>>1&1)) { // [0] moved
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

int bwa_pe_tpx(int iidx, const bntseq_t *bns, const bwt_t *bwt, int n_seqs1, int n_seqs2, bwa_seq_t *seqs[2], isize_info_t *ii,
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
	g_hash[iidx] = kh_init(b128);
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

			b128_t x;
			int j, k;
			long long n_occ[2];

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
					if (0 && r->l - r->k + 1 >= MIN_HASH_WIDTH) { // then check hash table
						b128_t key;
						int ret;
						key.x = r->k; key.y = r->l;
						khint_t iter = kh_put(b128, g_hash[iidx], key, &ret);
						if (ret) { // not in the hash table; ret must equal 1 as we never remove elements
							poslist_t *z = &kh_val(g_hash[iidx], iter);
							z->n = r->l - r->k + 1;
							z->a = (bwtint_t*)malloc(sizeof(bwtint_t) * z->n);
							for (l = r->k; l <= r->l; ++l) {
								int strand;
								z->a[l - r->k] = bwa_sa2pos(bns, bwt, l, p[j]->len, &strand)<<1;
								z->a[l - r->k] |= strand;
							}
						}
						for (l = 0; l < kh_val(g_hash[iidx], iter).n; ++l) {
							x.x = kh_val(g_hash[iidx], iter).a[l]>>1;
							x.y = k<<2 | (kh_val(g_hash[iidx], iter).a[l]&1)<<1 | j;
							kv_push(b128_t, d[iidx]->arr, x);
						}
					} else { // then calculate on the fly
						for (l = r->k; l <= r->l; ++l) {
							int strand;
							x.x = bwa_sa2pos(bns, bwt, l, p[j]->len, &strand);
							x.y = k<<2 | strand<<1 | j;
							kv_push(b128_t, d[iidx]->arr, x);
						}
					}
				}
			}

			cnt_chg += pairing(p, d[iidx], opt, gopt->s_mm, ii);
		}

                if (opt->N_multi || opt->n_multi) {
                        for (j = 0; j < 2; ++j) {
                                if (p[j]->type != BWA_TYPE_NO_MATCH) {

                                        int k, n_multi;

                                        if (!(p[j]->extra_flag&SAM_FPP) && p[1-j]->type != BWA_TYPE_NO_MATCH) {
                                                bwa_aln2seq_core(d[iidx]->aln[j].n, d[iidx]->aln[j].a, p[j], 0, p[j]->c1+p[j]->c2-1 > opt->N_multi? opt->n_multi : opt->N_multi);
                                        } else {
						bwa_aln2seq_core(d[iidx]->aln[j].n, d[iidx]->aln[j].a, p[j], 0, opt->n_multi);
					}

                                        for (k = 0, n_multi = 0; k < p[j]->n_multi; ++k) {
                                                int strand;
                                                bwt_multi1_t *q = p[j]->multi + k;
                                                q->pos = bwa_sa2pos(bns, bwt, q->pos, p[j]->len, &strand);
                                                q->strand = strand;
                                                if (q->pos != p[j]->pos)
                                                        p[j]->multi[n_multi++] = *q;
                                        }
                                        p[j]->n_multi = n_multi;
                                }
                        }
                }

	}

#ifdef _USE_LOCAL_GHASH
	khint_t iter;
	for (iter = kh_begin(g_hash[iidx]); iter != kh_end(g_hash[iidx]); ++iter){
		if (kh_exist(g_hash[iidx], iter)) free(kh_val(g_hash[iidx], iter).a);
	}

	kh_destroy(b128, g_hash[iidx]);
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

