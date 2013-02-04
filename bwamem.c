#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "bwamem.h"
#include "kvec.h"
#include "bntseq.h"
#include "ksw.h"

void mem_fill_scmat(int a, int b, int8_t mat[25])
{
	int i, j, k;
	for (i = k = 0; i < 5; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? a : -b;
		mat[k++] = 0; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = 0;
}

mem_opt_t *mem_opt_init()
{
	mem_opt_t *o;
	o = calloc(1, sizeof(mem_opt_t));
	o->a = 1; o->b = 5; o->q = 8; o->r = 1; o->w = 100;
	o->min_seed_len = 17;
	o->max_occ = 10;
	o->max_chain_gap = 10000;
	mem_fill_scmat(o->a, o->b, o->mat);
	return o;
}

/***************************
 * SMEM iterator interface *
 ***************************/

struct __smem_i {
	const bwt_t *bwt;
	const uint8_t *query;
	int start, len;
	bwtintv_v *matches; // matches; to be returned by smem_next()
	bwtintv_v *sub;     // sub-matches inside the longest match; temporary
	bwtintv_v *tmpvec[2]; // temporary arrays
};

smem_i *smem_itr_init(const bwt_t *bwt)
{
	smem_i *itr;
	itr = calloc(1, sizeof(smem_i));
	itr->bwt = bwt;
	itr->tmpvec[0] = calloc(1, sizeof(bwtintv_v));
	itr->tmpvec[1] = calloc(1, sizeof(bwtintv_v));
	itr->matches   = calloc(1, sizeof(bwtintv_v));
	itr->sub       = calloc(1, sizeof(bwtintv_v));
	return itr;
}

void smem_itr_destroy(smem_i *itr)
{
	free(itr->tmpvec[0]->a); free(itr->tmpvec[0]);
	free(itr->tmpvec[1]->a); free(itr->tmpvec[1]);
	free(itr->matches->a);   free(itr->matches);
	free(itr->sub->a);       free(itr->sub);
	free(itr);
}

void smem_set_query(smem_i *itr, int len, const uint8_t *query)
{
	itr->query = query;
	itr->start = 0;
	itr->len = len;
}

const bwtintv_v *smem_next(smem_i *itr, int split_len)
{
	int i, max, max_i;
	itr->tmpvec[0]->n = itr->tmpvec[1]->n = itr->matches->n = itr->sub->n = 0;
	if (itr->start >= itr->len || itr->start < 0) return 0;
	while (itr->start < itr->len && itr->query[itr->start] > 3) ++itr->start; // skip ambiguous bases
	if (itr->start == itr->len) return 0;
	itr->start = bwt_smem1(itr->bwt, itr->len, itr->query, itr->start, 1, itr->matches, itr->tmpvec); // search for SMEM
	if (itr->matches->n == 0) return itr->matches; // well, in theory, we should never come here
	for (i = max = 0, max_i = 0; i < itr->matches->n; ++i) { // look for the longest match
		bwtintv_t *p = &itr->matches->a[i];
		int len = (uint32_t)p->info - (p->info>>32);
		if (max < len) max = len, max_i = i;
	}
	if (split_len > 0 && max >= split_len && itr->matches->a[max_i].x[2] == 1) { // if the longest SMEM is unique and long
		int j;
		bwtintv_v *a = itr->tmpvec[0]; // reuse tmpvec[0] for merging
		bwtintv_t *p = &itr->matches->a[max_i];
		bwt_smem1(itr->bwt, itr->len, itr->query, ((uint32_t)p->info + (p->info>>32))>>1, 2, itr->sub, itr->tmpvec); // starting from the middle of the longest MEM
		i = j = 0; a->n = 0;
		while (i < itr->matches->n && j < itr->sub->n) { // ordered merge
			if (itr->matches->a[i].info < itr->sub->a[j].info) {
				kv_push(bwtintv_t, *a, itr->matches->a[i]);
				++i;
			} else {
				kv_push(bwtintv_t, *a, itr->sub->a[j]);
				++j;
			}
		}
		for (; i < itr->matches->n; ++i) kv_push(bwtintv_t, *a, itr->matches->a[i]);
		for (; j < itr->sub->n; ++j)     kv_push(bwtintv_t, *a, itr->sub->a[j]);
		kv_copy(bwtintv_t, *itr->matches, *a);
	}
	return itr->matches;
}

#include "kbtree.h"

#define chain_cmp(a, b) ((a).pos - (b).pos)
KBTREE_INIT(chn, mem_chain1_t, chain_cmp)

static int test_and_merge(const mem_opt_t *opt, mem_chain1_t *c, const mem_seed_t *p)
{
	int64_t qend, rend, x, y;
	const mem_seed_t *last = &c->seeds[c->n-1];
	qend = last->qbeg + last->len;
	rend = last->rbeg + last->len;
	if (p->qbeg >= c->seeds[0].qbeg && p->qbeg + p->len <= qend && p->rbeg >= c->seeds[0].rbeg && p->rbeg + p->len <= rend)
		return 1; // contained seed; do nothing
	x = p->qbeg - last->qbeg; // always positive
	y = p->rbeg - last->rbeg;
	if (y > 0 && x - y <= opt->w && y - x <= opt->w && x - last->len < opt->max_chain_gap && y - last->len < opt->max_chain_gap) { // grow the chain
		if (c->n == c->m) {
			c->m <<= 1;
			c->seeds = realloc(c->seeds, c->m * sizeof(mem_seed_t));
		}
		c->seeds[c->n++] = *p;
		return 1;
	}
	return 0; // request to add a new chain
}

static void mem_insert_seed(const mem_opt_t *opt, kbtree_t(chn) *tree, smem_i *itr)
{
	const bwtintv_v *a;
	while ((a = smem_next(itr, opt->min_seed_len<<1)) != 0) { // to find all SMEM and some internal MEM
		int i;
		for (i = 0; i < a->n; ++i) { // go through each SMEM/MEM up to itr->start
			bwtintv_t *p = &a->a[i];
			int slen = (uint32_t)p->info - (p->info>>32); // seed length
			int64_t k;
			if (slen < opt->min_seed_len || p->x[2] > opt->max_occ) continue; // ignore if too short or too repetitive
			for (k = 0; k < p->x[2]; ++k) {
				mem_chain1_t tmp, *lower, *upper;
				mem_seed_t s;
				int to_add = 0;
				s.rbeg = tmp.pos = bwt_sa(itr->bwt, p->x[0] + k); // this is the base coordinate in the forward-reverse reference
				s.qbeg = p->info>>32;
				s.len  = slen;
				if (kb_size(tree)) {
					kb_intervalp(chn, tree, &tmp, &lower, &upper); // find the closest chain
					if (!lower || !test_and_merge(opt, lower, &s)) to_add = 1;
				} else to_add = 1;
				if (to_add) { // add the seed as a new chain
					tmp.n = 1; tmp.m = 4;
					tmp.seeds = calloc(tmp.m, sizeof(mem_seed_t));
					tmp.seeds[0] = s;
					kb_putp(chn, tree, &tmp);
				}
			}
		}
	}
}

mem_chain_t mem_chain(const mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq)
{
	mem_chain_t chain;
	smem_i *itr;
	kbtree_t(chn) *tree;

	memset(&chain, 0, sizeof(mem_chain_t));
	if (len < opt->min_seed_len) return chain; // if the query is shorter than the seed length, no match
	tree = kb_init(chn, KB_DEFAULT_SIZE);
	itr = smem_itr_init(bwt);
	smem_set_query(itr, len, seq);
	mem_insert_seed(opt, tree, itr);

	chain.m = kb_size(tree); chain.n = 0;
	chain.chains = malloc(chain.m * sizeof(mem_chain1_t));

	#define traverse_func(p_) (chain.chains[chain.n++] = *(p_))
	__kb_traverse(mem_chain1_t, tree, traverse_func);
	#undef traverse_func

	smem_itr_destroy(itr);
	kb_destroy(chn, tree);
	return chain;
}

static inline int cal_max_gap(const mem_opt_t *opt, int qlen)
{
	int l = (int)((double)(qlen * opt->a - opt->q) / opt->r + 1.);
	return l > 1? l : 1;
}

mem_aln_t mem_chain2aln(const mem_opt_t *opt, int64_t l_pac, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain1_t *c)
{
	mem_aln_t a;
	int i, j, qbeg, qend;
	int64_t rlen, rbeg, rend, rmax[2], tmp;
	mem_seed_t *s;
	uint8_t *rseq = 0;

	memset(&a, 0, sizeof(mem_aln_t));
	// get the start and end of the seeded region
	rbeg = c->seeds[0].rbeg; qbeg = c->seeds[0].qbeg;
	s = &c->seeds[c->n-1];
	rend = s->rbeg + s->len; qend = s->qbeg + s->len;
	// get the max possible span
	rmax[0] = rbeg - (qbeg + cal_max_gap(opt, qbeg));
	rmax[1] = rend + ((l_query - qend) + cal_max_gap(opt, l_query - qend));
	if (rmax[0] < 0) rmax[0] = 0;
	if (rmax[1] > l_pac<<1) rmax[1] = l_pac<<1;
	// retrieve the reference sequence
	rseq = bns_get_seq(l_pac, pac, rmax[0], rmax[1], &rlen);

	if (qbeg) { // left extension of the first seed
		uint8_t *rs, *qs;
		int qle, tle;
		qs = malloc(qbeg);
		for (i = 0; i < qbeg; ++i) qs[i] = query[qbeg - 1 - i];
		tmp = rbeg - rmax[0];
		rs = malloc(tmp);
		for (i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i];
		a.score = ksw_extend(qbeg, qs, tmp, rs, 5, opt->mat, opt->q, opt->r, opt->w, c->seeds[0].len * opt->a, 0, &qle, &tle);
		a.qb = qbeg - qle; a.rb = rbeg - tle;
		free(qs); free(rs);
	} else a.score = c->seeds[0].len * opt->a, a.qb = 0, a.rb = rbeg;

	s = &c->seeds[0];
	if (s->qbeg + s->len != l_query) { // right extension of the first seed
		int qle, tle, qe, re;
		qe = s->qbeg + s->len; re = s->rbeg + s->len - rmax[0];
		for (j = 0; j < l_query - qe; ++j) putchar("ACGTN"[(int)query[j+qe]]); putchar('\n');
		for (j = 0; j < rmax[1] - rmax[0] - re; ++j) putchar("ACGTN"[(int)rseq[j+re]]); putchar('\n');
		a.score = ksw_extend(l_query - qe, query + qe, rmax[1] - rmax[0] - re, rseq + re, 5, opt->mat, opt->q, opt->r, opt->w, a.score, 0, &qle, &tle);
		a.qe = qe + qle; a.re = rmax[0] + re + tle;
	} else a.qe = l_query, a.re = s->rbeg + s->len;

	a.is_all = 1;
	if (c->n > 1) { // check if all the seeds have been included
		s = &c->seeds[c->n - 1];
		if (s->qbeg + s->len > a.qe) a.is_all = 0;
	}

	printf("[%d] score=%d\t[%d,%d) <=> [%lld,%lld)\tis_all=%d\n", c->n, a.score, a.qb, a.qe, a.rb, a.re, a.is_all);

	free(rseq);
	return a;
}
