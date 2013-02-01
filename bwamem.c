#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "bwamem.h"
#include "kvec.h"
#include "bntseq.h"

mem_opt_t *mem_opt_init()
{
	mem_opt_t *o;
	o = calloc(1, sizeof(mem_opt_t));
	o->a = 1; o->b = 9; o->q = 16; o->r = 1; o->w = 100;
	o->min_seed_len = 17;
	o->max_occ = 10;
	o->max_chain_gap = 10000;
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

mem_aln_t mem_chain2aln(int64_t l_pac, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain1_t *c)
{
	mem_aln_t a;
	int i, j;
	int64_t len;
	for (i = 0; i < c->n; ++i) {
		mem_seed_t *s = &c->seeds[i];
		uint8_t *seq = bns_get_seq(l_pac, pac, s->rbeg, s->rbeg + s->len, &len);
		for (j = 0; j < len; ++j) putchar("ACGTN"[seq[j]]); putchar('\n');
		for (j = 0; j < s->len; ++j) putchar("ACGTN"[query[j+s->qbeg]]); putchar('\n');
		free(seq);
	}
	return a;
}
