#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "bwamem.h"

memopt_t *mem_opt_init()
{
	memopt_t *o;
	o = calloc(1, sizeof(memopt_t));
	o->a = 1; o->b = 9; o->q = 16; o->r = 1; o->w = 100;
	o->min_seed_len = 17;
	o->max_occ = 10;
	o->max_chain_gap = 10000;
	return o;
}

/***************************
 * SMEM iterator interface *
 ***************************/

smem_i *smem_itr_init(const bwt_t *bwt)
{
	smem_i *itr;
	itr = calloc(1, sizeof(smem_i));
	itr->bwt = bwt;
	itr->tmpvec[0] = calloc(1, sizeof(bwtintv_v));
	itr->tmpvec[1] = calloc(1, sizeof(bwtintv_v));
	itr->matches   = calloc(1, sizeof(bwtintv_v));
	return itr;
}

void smem_itr_destroy(smem_i *itr)
{
	free(itr->tmpvec[0]->a);
	free(itr->tmpvec[1]->a);
	free(itr->matches->a);
	free(itr);
}

void smem_set_query(smem_i *itr, int min_intv, int len, const uint8_t *query)
{
	itr->query = query;
	itr->start = 0;
	itr->len = len;
	itr->min_intv = min_intv;
}

int smem_next(smem_i *itr)
{
	itr->tmpvec[0]->n = itr->tmpvec[1]->n = itr->matches->n = 0;
	if (itr->start >= itr->len || itr->start < 0) return -1;
	while (itr->start < itr->len && itr->query[itr->start] > 3) ++itr->start; // skip ambiguous bases
	if (itr->start == itr->len) return -1;
	itr->start = bwt_smem1(itr->bwt, itr->len, itr->query, itr->start, itr->min_intv, itr->matches, itr->tmpvec);
	return itr->start;
}

#include "kbtree.h"

#define chain_lt(a, b) ((a).pos < (b).pos)
KBTREE_INIT(chn, memchain1_t, chain_lt)

static int test_and_merge(const memopt_t *opt, memchain1_t *c, const memseed_t *p)
{
	int64_t qend, rend, x, y;
	const memseed_t *last = &c->seeds[c->n-1];
	qend = last->qbeg + last->len;
	rend = last->rbeg + last->len;
	if (p->qbeg > c->seeds[0].qbeg && p->qbeg + p->len < qend && p->rbeg > c->seeds[0].rbeg && p->rbeg + p->len < rend)
		return 1; // contained seed; do nothing
	x = p->qbeg - last->qbeg; // always positive
	y = p->rbeg - last->rbeg;
	if (y > 0 && x - y <= opt->w && y - x <= opt->w && x - last->len < opt->max_chain_gap && y - last->len < opt->max_chain_gap) {
		if (c->n == c->m) {
			c->m <<= 1;
			c->seeds = realloc(c->seeds, c->m * sizeof(memseed_t));
		}
		c->seeds[c->n++] = *p;
		return 1;
	}
	return 0;
}

void mem_insert_seed(const memopt_t *opt, kbtree_t(chn) *tree, smem_i *itr)
{
	while (smem_next(itr) > 0) {
		int i;
		for (i = 0; i < itr->matches->n; ++i) {
			bwtintv_t *p = &itr->matches->a[i];
			int slen = (uint32_t)p->info - (p->info>>32); // seed length
			int64_t k;
			if (slen >= opt->min_seed_len || p->x[2] > opt->max_occ) continue;
			for (k = 0; k < p->x[2]; ++k) {
				memchain1_t tmp, *lower, *upper;
				memseed_t c1;
				int to_add = 0;
				c1.rbeg = tmp.pos = bwt_sa(itr->bwt, p->x[0] + k);
				c1.qbeg = p->info>>32;
				c1.len  = slen;
				if (kb_size(tree)) {
					kb_intervalp(chn, tree, &tmp, &lower, &upper);
					if (!test_and_merge(opt, lower, &c1)) to_add = 1;
				} to_add = 1;
				if (to_add) {
					tmp.n = 1; tmp.m = 4;
					tmp.seeds = calloc(tmp.m, sizeof(memseed_t));
					kb_putp(chn, tree, &tmp);
				}
			}
		}
	}
}

memchain_t mem_collect_seed(const memopt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq)
{
	memchain_t chain;
	smem_i *itr;
	kbtree_t(chn) *tree;

	memset(&chain, 0, sizeof(memchain_t));
	if (len < opt->min_seed_len) return chain; // if the query is shorter than the seed length, no match
	tree = kb_init(chn, KB_DEFAULT_SIZE);
	itr = smem_itr_init(bwt);
	smem_set_query(itr, 1, len, seq);

	smem_itr_destroy(itr);
	kb_destroy(chn, tree);
	return chain;
}
