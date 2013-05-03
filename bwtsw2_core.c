#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/resource.h>
#include <assert.h>
#include "bwt_lite.h"
#include "bwtsw2.h"
#include "bwt.h"
#include "kvec.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

typedef struct {
	bwtint_t k, l;
} qintv_t;

#define qintv_eq(a, b) ((a).k == (b).k && (a).l == (b).l)
#define qintv_hash(a) ((a).k>>7^(a).l<<17)

#include "khash.h"
KHASH_INIT(qintv, qintv_t, uint64_t, 1, qintv_hash, qintv_eq)
KHASH_MAP_INIT_INT64(64, uint64_t)

#define MINUS_INF -0x3fffffff
#define MASK_LEVEL 0.90f

struct __mempool_t;
static void mp_destroy(struct __mempool_t*);
typedef struct {
	bwtint_t qk, ql;
	int I, D, G;
	uint32_t pj:2, qlen:30;
	int tlen;
	int ppos, upos;
	int cpos[4];
} bsw2cell_t;

#include "ksort.h"
KSORT_INIT_GENERIC(int)
#define __hitG_lt(a, b) (((a).G + ((int)(a).n_seeds<<2)) > (b).G + ((int)(b).n_seeds<<2))
KSORT_INIT(hitG, bsw2hit_t, __hitG_lt)

static const bsw2cell_t g_default_cell = { 0, 0, MINUS_INF, MINUS_INF, MINUS_INF, 0, 0, 0, -1, -1, {-1, -1, -1, -1} };

typedef struct {
	int n, max;
	uint32_t tk, tl; // this is fine
	bsw2cell_t *array;
} bsw2entry_t, *bsw2entry_p;

/* --- BEGIN: Stack operations --- */
typedef struct {
	int n_pending;
	kvec_t(bsw2entry_p) stack0, pending;
	struct __mempool_t *pool;
} bsw2stack_t;

#define stack_isempty(s) (kv_size(s->stack0) == 0 && s->n_pending == 0)
static void stack_destroy(bsw2stack_t *s) { mp_destroy(s->pool); kv_destroy(s->stack0); kv_destroy(s->pending); free(s); }
inline static void stack_push0(bsw2stack_t *s, bsw2entry_p e) { kv_push(bsw2entry_p, s->stack0, e); }
inline static bsw2entry_p stack_pop(bsw2stack_t *s)
{
	assert(!(kv_size(s->stack0) == 0 && s->n_pending != 0));
	return kv_pop(s->stack0);
}
/* --- END: Stack operations --- */

/* --- BEGIN: memory pool --- */
typedef struct __mempool_t {
	int cnt; // if cnt!=0, then there must be memory leak
	kvec_t(bsw2entry_p) pool;
} mempool_t;
inline static bsw2entry_p mp_alloc(mempool_t *mp)
{
	++mp->cnt;
	if (kv_size(mp->pool) == 0) return (bsw2entry_t*)calloc(1, sizeof(bsw2entry_t));
	else return kv_pop(mp->pool);
}
inline static void mp_free(mempool_t *mp, bsw2entry_p e)
{
	--mp->cnt; e->n = 0;
	kv_push(bsw2entry_p, mp->pool, e);
}
static void mp_destroy(struct __mempool_t *mp)
{
	int i;
	for (i = 0; i != kv_size(mp->pool); ++i) {
		free(kv_A(mp->pool, i)->array);
		free(kv_A(mp->pool, i));
	}
	kv_destroy(mp->pool);
	free(mp);
}
/* --- END: memory pool --- */

/* --- BEGIN: utilities --- */
static khash_t(64) *bsw2_connectivity(const bwtl_t *b)
{
	khash_t(64) *h;
	uint32_t k, l, cntk[4], cntl[4]; // this is fine
	uint64_t x;
	khiter_t iter;
	int j, ret;
	kvec_t(uint64_t) stack;

	kv_init(stack);
	h = kh_init(64);
	kh_resize(64, h, b->seq_len * 4);
	x = b->seq_len;
	kv_push(uint64_t, stack, x);
	while (kv_size(stack)) {
		x = kv_pop(stack);
		k = x>>32; l = (uint32_t)x;
		bwtl_2occ4(b, k-1, l, cntk, cntl);
		for (j = 0; j != 4; ++j) {
			k = b->L2[j] + cntk[j] + 1;
			l = b->L2[j] + cntl[j];
			if (k > l) continue;
			x = (uint64_t)k << 32 | l;
			iter = kh_put(64, h, x, &ret);
			if (ret) { // if not present
				kh_value(h, iter) = 1;
				kv_push(uint64_t, stack, x);
			} else ++kh_value(h, iter);
		}
	}
	kv_destroy(stack);
	//fprintf(stderr, "[bsw2_connectivity] %u nodes in the DAG\n", kh_size(h));
	return h;
}
// pick up top T matches at a node
static void cut_tail(bsw2entry_t *u, int T, bsw2entry_t *aux)
{
	int i, *a, n, x;
	if (u->n <= T) return;
	if (aux->max < u->n) {
		aux->max = u->n;
		aux->array = (bsw2cell_t*)realloc(aux->array, aux->max * sizeof(bsw2cell_t));
	}
	a = (int*)aux->array;
	for (i = n = 0; i != u->n; ++i)
		if (u->array[i].ql && u->array[i].G > 0)
			a[n++] = -u->array[i].G;
	if (n <= T) return;
	x = -ks_ksmall(int, n, a, T);
	n = 0;
	for (i = 0; i < u->n; ++i) {
		bsw2cell_t *p = u->array + i;
		if (p->G == x) ++n;
		if (p->G < x || (p->G == x && n >= T)) {
			p->qk = p->ql = 0; p->G = 0;
			if (p->ppos >= 0) u->array[p->ppos].cpos[p->pj] = -1;
		}
	}
}
// remove duplicated cells
static inline void remove_duplicate(bsw2entry_t *u, khash_t(qintv) *hash)
{
	int i, ret, j;
	khiter_t k;
	qintv_t key;
	kh_clear(qintv, hash);
	for (i = 0; i != u->n; ++i) {
		bsw2cell_t *p = u->array + i;
		if (p->ql == 0) continue;
		key.k = p->qk; key.l = p->ql;
		k = kh_put(qintv, hash, key, &ret);
		j = -1;
		if (ret == 0) {
			if ((uint32_t)kh_value(hash, k) >= p->G) j = i;
			else {
				j = kh_value(hash, k)>>32;
				kh_value(hash, k) = (uint64_t)i<<32 | p->G;
			}
		} else kh_value(hash, k) = (uint64_t)i<<32 | p->G;
		if (j >= 0) {
			p = u->array + j;
			p->qk = p->ql = 0; p->G = 0;
			if (p->ppos >= 0) u->array[p->ppos].cpos[p->pj] = -3;
		}
	}
}
// merge two entries
static void merge_entry(const bsw2opt_t * __restrict opt, bsw2entry_t *u, bsw2entry_t *v, bwtsw2_t *b)
{
	int i;
	if (u->n + v->n >= u->max) {
		u->max = u->n + v->n;
		u->array = (bsw2cell_t*)realloc(u->array, u->max * sizeof(bsw2cell_t));
	}
	for (i = 0; i != v->n; ++i) {
		bsw2cell_t *p = v->array + i;
		if (p->ppos >= 0) p->ppos += u->n;
		if (p->cpos[0] >= 0) p->cpos[0] += u->n;
		if (p->cpos[1] >= 0) p->cpos[1] += u->n;
		if (p->cpos[2] >= 0) p->cpos[2] += u->n;
		if (p->cpos[3] >= 0) p->cpos[3] += u->n;
	}
	memcpy(u->array + u->n, v->array, v->n * sizeof(bsw2cell_t));
	u->n += v->n;
}

static inline bsw2cell_t *push_array_p(bsw2entry_t *e)
{
	if (e->n == e->max) {
		e->max = e->max? e->max<<1 : 256;
		e->array = (bsw2cell_t*)realloc(e->array, sizeof(bsw2cell_t) * e->max);
	}
	return e->array + e->n;
}

static inline double time_elapse(const struct rusage *curr, const struct rusage *last)
{
	long t1 = (curr->ru_utime.tv_sec - last->ru_utime.tv_sec) + (curr->ru_stime.tv_sec - last->ru_stime.tv_sec);
	long t2 = (curr->ru_utime.tv_usec - last->ru_utime.tv_usec) + (curr->ru_stime.tv_usec - last->ru_stime.tv_usec);
	return (double)t1 + t2 * 1e-6;
}
/* --- END: utilities --- */

/* --- BEGIN: processing partial hits --- */
static void save_hits(const bwtl_t *bwt, int thres, bsw2hit_t *hits, bsw2entry_t *u)
{
	int i;
	uint32_t k; // this is fine
	for (i = 0; i < u->n; ++i) {
		bsw2cell_t *p = u->array + i;
		if (p->G < thres) continue;
		for (k = u->tk; k <= u->tl; ++k) {
			int beg, end;
			bsw2hit_t *q = 0;
			beg = bwt->sa[k]; end = beg + p->tlen;
			if (p->G > hits[beg*2].G) {
				hits[beg*2+1] = hits[beg*2];
				q = hits + beg * 2;
			} else if (p->G > hits[beg*2+1].G) q = hits + beg * 2 + 1;
			if (q) {
				q->k = p->qk; q->l = p->ql; q->len = p->qlen; q->G = p->G;
				q->beg = beg; q->end = end; q->G2 = q->k == q->l? 0 : q->G;
				q->flag = q->n_seeds = 0;
			}
		}
	}
}
/* "narrow hits" are node-to-node hits that have a high score and
 * are not so repetitive (|SA interval|<=IS). */
static void save_narrow_hits(const bwtl_t *bwtl, bsw2entry_t *u, bwtsw2_t *b1, int t, int IS)
{
	int i;
	for (i = 0; i < u->n; ++i) {
		bsw2hit_t *q;
		bsw2cell_t *p = u->array + i;
		if (p->G >= t && p->ql - p->qk + 1 <= IS) { // good narrow hit
			if (b1->max == b1->n) {
				b1->max = b1->max? b1->max<<1 : 4;
				b1->hits = realloc(b1->hits, b1->max * sizeof(bsw2hit_t));
			}
			q = &b1->hits[b1->n++];
			q->k = p->qk; q->l = p->ql;
			q->len = p->qlen;
			q->G = p->G; q->G2 = 0;
			q->beg = bwtl->sa[u->tk]; q->end = q->beg + p->tlen;
			q->flag = 0;
			// delete p
			p->qk = p->ql = 0; p->G = 0;
			if (p->ppos >= 0) u->array[p->ppos].cpos[p->pj] = -3;
		}
	}
}
/* after this, "narrow SA hits" will be expanded and the coordinates
 * will be obtained and stored in b->hits[*].k. */
int bsw2_resolve_duphits(const bntseq_t *bns, const bwt_t *bwt, bwtsw2_t *b, int IS)
{
	int i, j, n, is_rev;
	if (b->n == 0) return 0;
	if (bwt && bns) { // convert to chromosomal coordinates if requested
		int old_n = b->n;
		bsw2hit_t *old_hits = b->hits;
		for (i = n = 0; i < b->n; ++i) { // compute the memory to allocated
			bsw2hit_t *p = old_hits + i;
			if (p->l - p->k + 1 <= IS) n += p->l - p->k + 1;
			else if (p->G > 0) ++n;
		}
		b->n = b->max = n;
		b->hits = calloc(b->max, sizeof(bsw2hit_t));
		for (i = j = 0; i < old_n; ++i) {
			bsw2hit_t *p = old_hits + i;
			if (p->l - p->k + 1 <= IS) { // the hit is no so repetitive
				bwtint_t k;
				if (p->G == 0 && p->k == 0 && p->l == 0 && p->len == 0) continue;
				for (k = p->k; k <= p->l; ++k) {
					b->hits[j] = *p;
					b->hits[j].k = bns_depos(bns, bwt_sa(bwt, k), &is_rev);
					b->hits[j].l = 0;
					b->hits[j].is_rev = is_rev;
					if (is_rev) b->hits[j].k -= p->len - 1;
					++j;
				}
			} else if (p->G > 0) {
				b->hits[j] = *p;
				b->hits[j].k = bns_depos(bns, bwt_sa(bwt, p->k), &is_rev);
				b->hits[j].l = 0;
				b->hits[j].flag |= 1;
				b->hits[j].is_rev = is_rev;
				if (is_rev) b->hits[j].k -= p->len - 1;
				++j;
			}
		}
		free(old_hits);
	}
	for (i = j = 0; i < b->n; ++i) // squeeze out empty elements
		if (b->hits[i].G) b->hits[j++] = b->hits[i];
	b->n = j;
	ks_introsort(hitG, b->n, b->hits);
	for (i = 1; i < b->n; ++i) {
		bsw2hit_t *p = b->hits + i;
		for (j = 0; j < i; ++j) {
			bsw2hit_t *q = b->hits + j;
			int compatible = 1;
			if (p->is_rev != q->is_rev) continue; // hits from opposite strands are not duplicates
			if (p->l == 0 && q->l == 0) {
				int qol = (p->end < q->end? p->end : q->end) - (p->beg > q->beg? p->beg : q->beg); // length of query overlap
				if (qol < 0) qol = 0;
				if ((float)qol / (p->end - p->beg) > MASK_LEVEL || (float)qol / (q->end - q->beg) > MASK_LEVEL) {
					int64_t tol = (int64_t)(p->k + p->len < q->k + q->len? p->k + p->len : q->k + q->len)
						- (int64_t)(p->k > q->k? p->k : q->k); // length of target overlap
					if ((double)tol / p->len > MASK_LEVEL || (double)tol / q->len > MASK_LEVEL)
						compatible = 0;
				}
			}
			if (!compatible) {
				p->G = 0;
				if (q->G2 < p->G2) q->G2 = p->G2;
				break;
			}
		}
	}
	n = i;
	for (i = j = 0; i < n; ++i) {
		if (b->hits[i].G == 0) continue;
		if (i != j) b->hits[j++] = b->hits[i];
		else ++j;
	}
	b->n = j;
	return b->n;
}

int bsw2_resolve_query_overlaps(bwtsw2_t *b, float mask_level)
{
	int i, j, n;
	if (b->n == 0) return 0;
	ks_introsort(hitG, b->n, b->hits);
	{ // choose a random one
		int G0 = b->hits[0].G;
		for (i = 1; i < b->n; ++i)
			if (b->hits[i].G != G0) break;
		j = (int)(i * drand48());
		if (j) {
			bsw2hit_t tmp;
			tmp = b->hits[0]; b->hits[0] = b->hits[j]; b->hits[j] = tmp;
		}
	}
	for (i = 1; i < b->n; ++i) {
		bsw2hit_t *p = b->hits + i;
		int all_compatible = 1;
		if (p->G == 0) break;
		for (j = 0; j < i; ++j) {
			bsw2hit_t *q = b->hits + j;
			int64_t tol = 0;
			int qol, compatible = 0;
			float fol;
			if (q->G == 0) continue;
			qol = (p->end < q->end? p->end : q->end) - (p->beg > q->beg? p->beg : q->beg);
			if (qol < 0) qol = 0;
			if (p->l == 0 && q->l == 0) {
				tol = (int64_t)(p->k + p->len < q->k + q->len? p->k + p->len : q->k + q->len)
					- (p->k > q->k? p->k : q->k);
				if (tol < 0) tol = 0;
			}
			fol = (float)qol / (p->end - p->beg < q->end - q->beg? p->end - p->beg : q->end - q->beg);
			if (fol < mask_level || (tol > 0 && qol < p->end - p->beg && qol < q->end - q->beg)) compatible = 1;
			if (!compatible) {
				if (q->G2 < p->G) q->G2 = p->G;
				all_compatible = 0;
			}
		}
		if (!all_compatible) p->G = 0;
	}
	n = i;
	for (i = j = 0; i < n; ++i) {
		if (b->hits[i].G == 0) continue;
		if (i != j) b->hits[j++] = b->hits[i];
		else ++j;
	}
	b->n = j;
	return j;
}
/* --- END: processing partial hits --- */

/* --- BEGIN: global mem pool --- */
bsw2global_t *bsw2_global_init()
{
	bsw2global_t *pool;
	bsw2stack_t *stack;
	pool = calloc(1, sizeof(bsw2global_t));
	stack = calloc(1, sizeof(bsw2stack_t));
	stack->pool = (mempool_t*)calloc(1, sizeof(mempool_t));
	pool->stack = (void*)stack;
	return pool;
}

void bsw2_global_destroy(bsw2global_t *pool)
{
	stack_destroy((bsw2stack_t*)pool->stack);
	free(pool->aln_mem);
	free(pool);
}
/* --- END: global mem pool --- */

static inline int fill_cell(const bsw2opt_t *o, int match_score, bsw2cell_t *c[4])
{
	int G = c[3]? c[3]->G + match_score : MINUS_INF;
	if (c[1]) {
		c[0]->I = c[1]->I > c[1]->G - o->q? c[1]->I - o->r : c[1]->G - o->qr;
		if (c[0]->I > G) G = c[0]->I;
	} else c[0]->I = MINUS_INF;
	if (c[2]) {
		c[0]->D = c[2]->D > c[2]->G - o->q? c[2]->D - o->r : c[2]->G - o->qr;
		if (c[0]->D > G) G = c[0]->D;
	} else c[0]->D = MINUS_INF;
	return(c[0]->G = G);
}

static void init_bwtsw2(const bwtl_t *target, const bwt_t *query, bsw2stack_t *s)
{
	bsw2entry_t *u;
	bsw2cell_t *x;

	u = mp_alloc(s->pool);
	u->tk = 0; u->tl = target->seq_len;
	x = push_array_p(u);
	*x = g_default_cell;
	x->G = 0; x->qk = 0; x->ql = query->seq_len;
	u->n++;
	stack_push0(s, u);
}
/* On return, ret[1] keeps not-so-repetitive hits (narrow SA hits); ret[0] keeps all hits (right?) */
bwtsw2_t **bsw2_core(const bntseq_t *bns, const bsw2opt_t *opt, const bwtl_t *target, const bwt_t *query, bsw2global_t *pool)
{
	bsw2stack_t *stack = (bsw2stack_t*)pool->stack;
	bwtsw2_t *b, *b1, **b_ret;
	int i, j, score_mat[16], *heap, heap_size, n_tot = 0;
	struct rusage curr, last;
	khash_t(qintv) *rhash;
	khash_t(64) *chash;

	// initialize connectivity hash (chash)
	chash = bsw2_connectivity(target);
	// calculate score matrix
	for (i = 0; i != 4; ++i)
		for (j = 0; j != 4; ++j)
			score_mat[i<<2|j] = (i == j)? opt->a : -opt->b;
	// initialize other variables
	rhash = kh_init(qintv);
	init_bwtsw2(target, query, stack);
	heap_size = opt->z;
	heap = calloc(heap_size, sizeof(int));
	// initialize the return struct
	b = (bwtsw2_t*)calloc(1, sizeof(bwtsw2_t));
	b->n = b->max = target->seq_len * 2;
	b->hits = calloc(b->max, sizeof(bsw2hit_t));
	b1 = (bwtsw2_t*)calloc(1, sizeof(bwtsw2_t));
	b_ret = calloc(2, sizeof(void*));
	b_ret[0] = b; b_ret[1] = b1;
	// initialize timer
	getrusage(0, &last);
	// the main loop: traversal of the DAG
	while (!stack_isempty(stack)) {
		int old_n, tj;
		bsw2entry_t *v;
		uint32_t tcntk[4], tcntl[4];
		bwtint_t k, l;

		v = stack_pop(stack); old_n = v->n;
		n_tot += v->n;

		for (i = 0; i < v->n; ++i) { // test max depth and band width
			bsw2cell_t *p = v->array + i;
			if (p->ql == 0) continue;
			if (p->tlen - (int)p->qlen > opt->bw || (int)p->qlen - p->tlen > opt->bw) {
				p->qk = p->ql = 0;
				if (p->ppos >= 0) v->array[p->ppos].cpos[p->pj] = -5;
			}
		}

		// get Occ for the DAG
		bwtl_2occ4(target, v->tk - 1, v->tl, tcntk, tcntl);
		for (tj = 0; tj != 4; ++tj) { // descend to the children
			bwtint_t qcntk[4], qcntl[4];
			int qj, *curr_score_mat = score_mat + tj * 4;
			khiter_t iter;
			bsw2entry_t *u;

			k = target->L2[tj] + tcntk[tj] + 1;
			l = target->L2[tj] + tcntl[tj];
			if (k > l) continue;
			// update counter
			iter = kh_get(64, chash, (uint64_t)k<<32 | l);
			--kh_value(chash, iter);
			// initialization
			u = mp_alloc(stack->pool);
			u->tk = k; u->tl = l;
			memset(heap, 0, sizeof(int) * opt->z);
			// loop through all the nodes in v
		    for (i = 0; i < v->n; ++i) {
				bsw2cell_t *p = v->array + i, *x, *c[4]; // c[0]=>current, c[1]=>I, c[2]=>D, c[3]=>G
				int is_added = 0;
				if (p->ql == 0) continue; // deleted node
				c[0] = x = push_array_p(u);
				x->G = MINUS_INF;
				p->upos = x->upos = -1;
				if (p->ppos >= 0) { // parent has been visited
					c[1] = (v->array[p->ppos].upos >= 0)? u->array + v->array[p->ppos].upos : 0;
					c[3] = v->array + p->ppos; c[2] = p;
					if (fill_cell(opt, curr_score_mat[p->pj], c) > 0) { // then update topology at p and x
						x->ppos = v->array[p->ppos].upos; // the parent pos in u
						p->upos = u->n++; // the current pos in u
						if (x->ppos >= 0) u->array[x->ppos].cpos[p->pj] = p->upos; // the child pos of its parent in u
						is_added = 1;
					}
				} else {
					x->D = p->D > p->G - opt->q? p->D - opt->r : p->G - opt->qr;
					if (x->D > 0) {
						x->G = x->D;
						x->I = MINUS_INF; x->ppos = -1;
						p->upos = u->n++;
						is_added = 1;
					}
				}
				if (is_added) { // x has been added to u->array. fill the remaining variables
					x->cpos[0] = x->cpos[1] = x->cpos[2] = x->cpos[3] = -1;
					x->pj = p->pj; x->qk = p->qk; x->ql = p->ql; x->qlen = p->qlen; x->tlen = p->tlen + 1;
					if (x->G > -heap[0]) {
						heap[0] = -x->G;
						ks_heapadjust(int, 0, heap_size, heap);
					}
				}
				if ((x->G > opt->qr && x->G >= -heap[0]) || i < old_n) { // good node in u, or in v
					if (p->cpos[0] == -1 || p->cpos[1] == -1 || p->cpos[2] == -1 || p->cpos[3] == -1) {
						bwt_2occ4(query, p->qk - 1, p->ql, qcntk, qcntl);
						for (qj = 0; qj != 4; ++qj) { // descend to the prefix trie
							if (p->cpos[qj] != -1) continue; // this node will be visited later
							k = query->L2[qj] + qcntk[qj] + 1;
							l = query->L2[qj] + qcntl[qj];
							if (k > l) { p->cpos[qj] = -2; continue; }
							x = push_array_p(v);
							p = v->array + i; // p may not point to the correct position after realloc
							x->G = x->I = x->D = MINUS_INF;
							x->qk = k; x->ql = l; x->pj = qj; x->qlen = p->qlen + 1; x->ppos = i; x->tlen = p->tlen;
							x->cpos[0] = x->cpos[1] = x->cpos[2] = x->cpos[3] = -1;
							p->cpos[qj] = v->n++;
						} // ~for(qj)
					} // ~if(p->cpos[])
				} // ~if
			} // ~for(i)
			if (u->n) save_hits(target, opt->t, b->hits, u);
			{ // push u to the stack (or to the pending array)
				uint32_t cnt, pos;
				cnt = (uint32_t)kh_value(chash, iter);
				pos = kh_value(chash, iter)>>32;
				if (pos) { // something in the pending array, then merge
					bsw2entry_t *w = kv_A(stack->pending, pos-1);
					if (u->n) {
						if (w->n < u->n) { // swap
							w = u; u = kv_A(stack->pending, pos-1); kv_A(stack->pending, pos-1) = w;
						}
						merge_entry(opt, w, u, b);
					}
					if (cnt == 0) { // move from pending to stack0
						remove_duplicate(w, rhash);
						save_narrow_hits(target, w, b1, opt->t, opt->is);
						cut_tail(w, opt->z, u);
						stack_push0(stack, w);
						kv_A(stack->pending, pos-1) = 0;
						--stack->n_pending;
					}
					mp_free(stack->pool, u);
				} else if (cnt) { // the first time
					if (u->n) { // push to the pending queue
						++stack->n_pending;
						kv_push(bsw2entry_p, stack->pending, u);
						kh_value(chash, iter) = (uint64_t)kv_size(stack->pending)<<32 | cnt;
					} else mp_free(stack->pool, u);
				} else { // cnt == 0, then push to the stack
					bsw2entry_t *w = mp_alloc(stack->pool);
					save_narrow_hits(target, u, b1, opt->t, opt->is);
					cut_tail(u, opt->z, w);
					mp_free(stack->pool, w);
					stack_push0(stack, u);
				}
			}
		} // ~for(tj)
		mp_free(stack->pool, v);
	} // while(top)
	getrusage(0, &curr);
	for (i = 0; i < 2; ++i)
		for (j = 0; j < b_ret[i]->n; ++j)
			b_ret[i]->hits[j].n_seeds = 0;
	bsw2_resolve_duphits(bns, query, b, opt->is);
	bsw2_resolve_duphits(bns, query, b1, opt->is);
	//fprintf(stderr, "stats: %.3lf sec; %d elems\n", time_elapse(&curr, &last), n_tot);
	// free
	free(heap);
	kh_destroy(qintv, rhash);
	kh_destroy(64, chash);
	stack->pending.n = stack->stack0.n = 0;
	return b_ret;
}
