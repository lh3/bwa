#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif
#include "kstring.h"
#include "bwamem.h"
#include "bntseq.h"
#include "ksw.h"
#include "ksort.h"

#define MAPQ_COEF 40.

int mem_debug = 0;

void mem_fill_scmat(int a, int b, int8_t mat[25])
{
	int i, j, k;
	for (i = k = 0; i < 4; ++i) {
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
	o->min_seed_len = 19;
	o->max_occ = 50;
	o->max_chain_gap = 10000;
	o->mask_level = 0.50;
	o->chain_drop_ratio = 0.50;
	o->chunk_size = 10000000;
	o->n_threads = 1;
	o->pe_dir = 0<<1|1;
	o->is_pe = 0;
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
			int64_t xi = itr->matches->a[i].info>>32<<32 | (itr->len - (uint32_t)itr->matches->a[i].info);
			int64_t xj = itr->sub->a[j].info>>32<<32 | (itr->len - (uint32_t)itr->sub->a[j].info);
			if (xi < xj) {
				kv_push(bwtintv_t, *a, itr->matches->a[i]);
				++i;
			} else if ((uint32_t)itr->sub->a[j].info - (itr->sub->a[j].info>>32) >= max>>1) {
				kv_push(bwtintv_t, *a, itr->sub->a[j]);
				++j;
			} else ++j;
		}
		for (; i < itr->matches->n; ++i) kv_push(bwtintv_t, *a, itr->matches->a[i]);
		for (; j < itr->sub->n; ++j)
			if ((uint32_t)itr->sub->a[j].info - (itr->sub->a[j].info>>32) >= max>>1)
				kv_push(bwtintv_t, *a, itr->sub->a[j]);
		kv_copy(bwtintv_t, *itr->matches, *a);
	}
	return itr->matches;
}

/********************************
 * Chaining while finding SMEMs *
 ********************************/

#include "kbtree.h"

#define chain_cmp(a, b) ((a).pos - (b).pos)
KBTREE_INIT(chn, mem_chain_t, chain_cmp)

static int test_and_merge(const mem_opt_t *opt, mem_chain_t *c, const mem_seed_t *p)
{
	int64_t qend, rend, x, y;
	const mem_seed_t *last = &c->seeds[c->n-1];
	qend = last->qbeg + last->len;
	rend = last->rbeg + last->len;
	if (p->qbeg >= c->seeds[0].qbeg && p->qbeg + p->len <= qend && p->rbeg >= c->seeds[0].rbeg && p->rbeg + p->len <= rend)
		return 1; // contained seed; do nothing
	x = p->qbeg - last->qbeg; // always non-negtive
	y = p->rbeg - last->rbeg;
	if (y >= 0 && x - y <= opt->w && y - x <= opt->w && x - last->len < opt->max_chain_gap && y - last->len < opt->max_chain_gap) { // grow the chain
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
				mem_chain_t tmp, *lower, *upper;
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

void mem_print_chain(const bntseq_t *bns, mem_chain_v *chn)
{
	int i, j;
	for (i = 0; i < chn->n; ++i) {
		mem_chain_t *p = &chn->a[i];
		printf("%d", p->n);
		for (j = 0; j < p->n; ++j) {
			bwtint_t pos;
			int is_rev, ref_id;
			pos = bns_depos(bns, p->seeds[j].rbeg, &is_rev);
			if (is_rev) pos -= p->seeds[j].len - 1;
			bns_cnt_ambi(bns, pos, p->seeds[j].len, &ref_id);
			printf("\t%d,%d,%ld(%s:%c%ld)", p->seeds[j].len, p->seeds[j].qbeg, (long)p->seeds[j].rbeg, bns->anns[ref_id].name, "+-"[is_rev], (long)(pos - bns->anns[ref_id].offset) + 1);
		}
		putchar('\n');
	}
}

mem_chain_v mem_chain(const mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq)
{
	mem_chain_v chain;
	smem_i *itr;
	kbtree_t(chn) *tree;

	kv_init(chain);
	if (len < opt->min_seed_len) return chain; // if the query is shorter than the seed length, no match
	tree = kb_init(chn, KB_DEFAULT_SIZE);
	itr = smem_itr_init(bwt);
	smem_set_query(itr, len, seq);
	mem_insert_seed(opt, tree, itr);

	kv_resize(mem_chain_t, chain, kb_size(tree));

	#define traverse_func(p_) (chain.a[chain.n++] = *(p_))
	__kb_traverse(mem_chain_t, tree, traverse_func);
	#undef traverse_func

	smem_itr_destroy(itr);
	kb_destroy(chn, tree);
	return chain;
}

/********************
 * Filtering chains *
 ********************/

typedef struct {
	int beg, end, w;
	void *p, *p2;
} flt_aux_t;

#define flt_lt(a, b) ((a).w > (b).w)
KSORT_INIT(mem_flt, flt_aux_t, flt_lt)

int mem_chain_flt(const mem_opt_t *opt, int n_chn, mem_chain_t *chains)
{
	flt_aux_t *a;
	int i, j, n;
	if (n_chn <= 1) return n_chn; // no need to filter
	a = malloc(sizeof(flt_aux_t) * n_chn);
	for (i = 0; i < n_chn; ++i) {
		mem_chain_t *c = &chains[i];
		int w = 0;
		for (j = 0; j < c->n; ++j) w += c->seeds[j].len; // FIXME: take care of seed overlaps
		a[i].beg = c->seeds[0].qbeg;
		a[i].end = c->seeds[c->n-1].qbeg + c->seeds[c->n-1].len;
		a[i].w = w; a[i].p = c; a[i].p2 = 0;
	}
	ks_introsort(mem_flt, n_chn, a);
	{ // reorder chains such that the best chain appears first
		mem_chain_t *swap;
		swap = malloc(sizeof(mem_chain_t) * n_chn);
		for (i = 0; i < n_chn; ++i) {
			swap[i] = *((mem_chain_t*)a[i].p);
			a[i].p = &chains[i]; // as we will memcpy() below, a[i].p is changed
		}
		memcpy(chains, swap, sizeof(mem_chain_t) * n_chn);
		free(swap);
	}
	for (i = 1, n = 1; i < n_chn; ++i) {
		for (j = 0; j < n; ++j) {
			int b_max = a[j].beg > a[i].beg? a[j].beg : a[i].beg;
			int e_min = a[j].end < a[i].end? a[j].end : a[i].end;
			if (e_min > b_max) { // have overlap
				int min_l = a[i].end - a[i].beg < a[j].end - a[j].beg? a[i].end - a[i].beg : a[j].end - a[j].beg;
				if (e_min - b_max >= min_l * opt->mask_level) { // significant overlap
					if (a[j].p2 == 0) a[j].p2 = a[i].p;
					if (a[i].w < a[j].w * opt->chain_drop_ratio)
						break;
				}
			}
		}
		if (j == n) a[n++] = a[i]; // if have no significant overlap with better chains, keep it.
	}
	for (i = 0; i < n; ++i) { // mark chains to be kept
		mem_chain_t *c = (mem_chain_t*)a[i].p;
		if (c->n > 0) c->n = -c->n;
		c = (mem_chain_t*)a[i].p2;
		if (c && c->n > 0) c->n = -c->n;
	}
	free(a);
	for (i = 0; i < n_chn; ++i) { // free discarded chains
		mem_chain_t *c = &chains[i];
		if (c->n >= 0) {
			free(c->seeds);
			c->n = c->m = 0;
		} else c->n = -c->n;
	}
	for (i = n = 0; i < n_chn; ++i) { // squeeze out discarded chains
		if (chains[i].n > 0) {
			if (n != i) chains[n++] = chains[i];
			else ++n;
		}
	}
	return n;
}

#define alnreg_lt(a, b) ((a).score > (b).score || ((a).score == (b).score && ((a).rb < (b).rb || ((a).rb == (b).rb && (a).qb < (b).qb))))
KSORT_INIT(mem_ar, mem_alnreg_t, alnreg_lt)

int mem_choose_alnreg_se(const mem_opt_t *opt, int n, mem_alnreg_t *a)
{ // similar to the loop in mem_chain_flt()
	int i, j, m;
	if (n <= 1) return n;
	ks_introsort(mem_ar, n, a);
	for (i = 1; i < n; ++i) { // mark identical hits
		if (a[i].score == a[i-1].score && a[i].rb == a[i-1].rb && a[i].qb == a[i-1].qb)
			a[i].score = 0;
	}
	for (i = 1, m = 1; i < n; ++i) // exclude identical hits
		if (a[i].score > 0) a[m++] = a[i];
	n = m;
	for (i = 0; i < n; ++i) a[i].sub = 0;
	for (i = 1, m = 1; i < n; ++i) {
		for (j = 0; j < m; ++j) {
			int b_max = a[j].qb > a[i].qb? a[j].qb : a[i].qb;
			int e_min = a[j].qe < a[i].qe? a[j].qe : a[i].qe;
			if (e_min > b_max) { // have overlap
				int min_l = a[i].qe - a[i].qb < a[j].qe - a[j].qb? a[i].qe - a[i].qb : a[j].qe - a[j].qb;
				if (e_min - b_max >= min_l * opt->mask_level) { // significant overlap
					if (a[j].sub == 0) a[j].sub = a[i].score;
					break;
				}
			}
		}
		if (j == m) a[m++] = a[i];
	}
	return m;
}

/****************************************
 * Construct the alignment from a chain *
 ****************************************/

static inline int cal_max_gap(const mem_opt_t *opt, int qlen)
{
	int l = (int)((double)(qlen * opt->a - opt->q) / opt->r + 1.);
	return l > 1? l : 1;
}

void mem_chain2aln(const mem_opt_t *opt, int64_t l_pac, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_t *a)
{ // FIXME: in general, we SHOULD check funny seed patterns such as contained seeds. When that happens, we should use a SW or extend more seeds
	int i, k;
	int64_t rlen, rmax[2], tmp, max = 0, max_i = 0;
	const mem_seed_t *s;
	uint8_t *rseq = 0;
	mem_alnreg_t best;

	memset(&best, 0, sizeof(mem_alnreg_t));
	// get the max possible span
	rmax[0] = l_pac<<1; rmax[1] = 0;
	for (i = 0; i < c->n; ++i) {
		int64_t b, e;
		const mem_seed_t *t = &c->seeds[i];
		b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
		e = t->rbeg + t->len + ((l_query - t->qbeg - t->len) + cal_max_gap(opt, l_query - t->qbeg - t->len));
		rmax[0] = rmax[0] < b? rmax[0] : b;
		rmax[1] = rmax[1] > e? rmax[1] : e;
		if (t->len > max) max = t->len, max_i = i;
	}
	// retrieve the reference sequence
	rseq = bns_get_seq(l_pac, pac, rmax[0], rmax[1], &rlen);

	for (k = 0; k < c->n;) {
		s = &c->seeds[k];
		memset(a, 0, sizeof(mem_alnreg_t));
		if (s->qbeg) { // left extension
			uint8_t *rs, *qs;
			int qle, tle;
			qs = malloc(s->qbeg);
			for (i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];
			tmp = s->rbeg - rmax[0];
			rs = malloc(tmp);
			for (i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i];
			a->score = ksw_extend(s->qbeg, qs, tmp, rs, 5, opt->mat, opt->q, opt->r, opt->w, s->len * opt->a, &qle, &tle);
			a->qb = s->qbeg - qle; a->rb = s->rbeg - tle;
			free(qs); free(rs);
		} else a->score = s->len * opt->a, a->qb = 0, a->rb = s->rbeg;

		if (s->qbeg + s->len != l_query) { // right extension of the first seed
			int qle, tle, qe, re;
			qe = s->qbeg + s->len;
			re = s->rbeg + s->len - rmax[0];
			a->score = ksw_extend(l_query - qe, query + qe, rmax[1] - rmax[0] - re, rseq + re, 5, opt->mat, opt->q, opt->r, opt->w, a->score, &qle, &tle);
			a->qe = qe + qle; a->re = rmax[0] + re + tle;
		} else a->qe = l_query, a->re = s->rbeg + s->len;
		if (mem_debug >= 2) printf("[%d] score=%d\t[%d,%d) <=> [%ld,%ld)\n", k, a->score, a->qb, a->qe, (long)a->rb, (long)a->re);
		// check how many seeds have been covered
		for (i = k + 1; i < c->n; ++i) {
			const mem_seed_t *t = &c->seeds[i];
			if (t->rbeg + t->len > a->re || t->qbeg + t->len > a->qe)
				break;
		}
		if (i >= c->n) break; // all seeds are included; no need to proceed
		if (a->score > best.score) best = *a;
		k = i;
	}
	if (a->score < best.score) *a = best;
	free(rseq);

	// compute seedcov
	if (c->n > 1) {
		for (i = 0, a->seedcov = 0; i < c->n; ++i) {
			s = &c->seeds[i];
			if (s->qbeg >= a->qb && s->qbeg + s->len <= a->qe && s->rbeg >= a->rb && s->rbeg + s->len <= a->re) // seed fully contained
				a->seedcov += s->len; // this is not very accurate, but for approx. mapQ, this is good enough
		}
	} else a->seedcov = c->seeds[0].len;
}

uint32_t *mem_gen_cigar(const mem_opt_t *opt, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar)
{
	uint32_t *cigar = 0;
	uint8_t tmp, *rseq;
	int i, w;
	int64_t rlen;
	*n_cigar = 0;
	if (l_query <= 0 || rb >= re || (rb < l_pac && re > l_pac)) return 0; // reject if negative length or bridging the forward and reverse strand
	rseq = bns_get_seq(l_pac, pac, rb, re, &rlen);
	if (re - rb != rlen) goto ret_gen_cigar; // possible if out of range
	if (rb >= l_pac) { // then reverse both query and rseq; this is to ensure indels to be placed at the leftmost position
		for (i = 0; i < l_query>>1; ++i)
			tmp = query[i], query[i] = query[l_query - 1 - i], query[l_query - 1 - i] = tmp;
		for (i = 0; i < rlen>>1; ++i)
			tmp = rseq[i], rseq[i] = rseq[rlen - 1 - i], rseq[rlen - 1 - i] = tmp;
	}
	//printf("[Q] "); for (i = 0; i < l_query; ++i) putchar("ACGTN"[(int)query[i]]); putchar('\n');
	//printf("[R] "); for (i = 0; i < re - rb; ++i) putchar("ACGTN"[(int)rseq[i]]); putchar('\n');
	// set the band-width
	w = (int)((double)(l_query * opt->a - opt->q) / opt->r + 1.);
	w = w < 1? w : 1;
	w = w < opt->w? w : opt->w;
	w += abs(rlen - l_query);
	// NW alignment
	*score = ksw_global(l_query, query, rlen, rseq, 5, opt->mat, opt->q, opt->r, w, n_cigar, &cigar);
	if (rb >= l_pac) // reverse back query
		for (i = 0; i < l_query>>1; ++i)
			tmp = query[i], query[i] = query[l_query - 1 - i], query[l_query - 1 - i] = tmp;

ret_gen_cigar:
	free(rseq);
	return cigar;
}

/************************
 * Integrated interface *
 ************************/

static inline int approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a)
{
	int mapq, l, sub = a->sub? a->sub : opt->min_seed_len * opt->a;
	double identity;
	l = a->qe - a->qb > a->re - a->rb? a->qe - a->qb : a->re - a->rb;
	mapq = a->score? (int)(MAPQ_COEF * (1. - (double)sub / a->score) * log(a->seedcov) + .499) : 0;
	identity = 1. - (double)(l * opt->a - a->score) / (opt->a + opt->b) / l;
	mapq = identity < 0.95? (int)(mapq * identity * identity + .499) : mapq;
	if (mapq > 60) mapq = 60;
	return mapq;
}

void mem_sam_se(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a)
{
	int i, k, m;
	kstring_t str;
	char *seq;

	str.l = str.m = 0; str.s = 0;
	m = mem_choose_alnreg_se(opt, a->n, a->a);
	seq = malloc(s->l_seq);
	if (m == 0) { // no seeds found
		for (i = 0; i < s->l_seq; ++i) seq[i] = "ACGTN"[(int)s->seq[i]];
		kputs(s->name, &str); kputs("\t8\t*\t0\t0\t*\t*\t0\t0\t", &str);
		kputsn(seq, s->l_seq, &str);
		if (s->qual) kputsn(s->qual, s->l_seq, &str);
		else kputc('*', &str);
		kputc('\n', &str);
		goto ret_sam_se;
	}
	for (k = 0; k < m; ++k) {
		uint32_t *cigar = 0;
		int score, is_rev, nn, rid, flag = 0, n_cigar = 0, mapq = 0;
		int64_t pos;
		mem_alnreg_t *p = &a->a[k];
		cigar = mem_gen_cigar(opt, bns->l_pac, pac, p->qe - p->qb, (uint8_t*)&s->seq[p->qb], p->rb, p->re, &score, &n_cigar);
		pos = bns_depos(bns, p->rb < bns->l_pac? p->rb : p->re - 1, &is_rev);
		nn = bns_cnt_ambi(bns, pos, p->re - p->rb, &rid);
		flag |= is_rev? 16 : 0;
		if (n_cigar == 0) flag |= 8;
		kputs(s->name, &str); kputc('\t', &str); kputw(flag, &str); kputc('\t', &str);
		kputs(bns->anns[rid].name, &str); kputc('\t', &str); kputuw(pos - bns->anns[rid].offset + 1, &str); kputc('\t', &str);
		mapq = approx_mapq_se(opt, p);
		kputw(mapq, &str); kputc('\t', &str);
		if (n_cigar) {
			int clip5, clip3;
			clip5 = is_rev? s->l_seq - p->qe : p->qb;
			clip3 = is_rev? p->qb : s->l_seq - p->qe;
			if (clip5) { kputw(clip5, &str); kputc('S', &str); }
			for (i = 0; i < n_cigar; ++i) {
				kputw(cigar[i]>>4, &str); kputc("MIDSH"[cigar[i]&0xf], &str);
			}
			if (clip3) { kputw(clip3, &str); kputc('S', &str); }
		} else kputc('*', &str);
		kputsn("\t*\t0\t0\t", 7, &str);
		if (is_rev) for (i = s->l_seq - 1; i >= 0; --i) seq[i] = "TGCAN"[(int)s->seq[i]];
		else for (i = 0; i < s->l_seq; ++i) seq[i] = "ACGTN"[(int)s->seq[i]];
		kputsn(seq, s->l_seq, &str); kputc('\t', &str);
		if (s->qual) kputsn(s->qual, s->l_seq, &str);
		else kputc('*', &str);
		kputsn("\tAS:i:", 6, &str); kputw(p->score, &str);
		kputsn("\tss:i:", 6, &str); kputw(p->sub, &str);
		kputsn("\tnw:i:", 6, &str); kputw(score, &str);
		kputc('\n', &str);
		free(cigar);
	}

ret_sam_se:
	free(seq);
	s->sam = str.s;
}

static mem_alnreg_v find_alnreg(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s)
{
	int i;
	mem_chain_v chn;
	mem_alnreg_v regs;
	for (i = 0; i < s->l_seq; ++i)
		s->seq[i] = nst_nt4_table[(int)s->seq[i]];
	chn = mem_chain(opt, bwt, s->l_seq, (uint8_t*)s->seq);
	chn.n = mem_chain_flt(opt, chn.n, chn.a);
	if (mem_debug >= 1) mem_print_chain(bns, &chn);
	regs.n = regs.m = chn.n;
	regs.a = malloc(regs.n * sizeof(mem_alnreg_t));
	for (i = 0; i < chn.n; ++i) {
		mem_chain2aln(opt, bns->l_pac, pac, s->l_seq, (uint8_t*)s->seq, &chn.a[i], &regs.a[i]);
		free(chn.a[i].seeds);
	}
	free(chn.a);
	return regs;
}

typedef struct {
	int start, step, n;
	const mem_opt_t *opt;
	const bwt_t *bwt;
	const bntseq_t *bns;
	const uint8_t *pac;
	bseq1_t *seqs;
	mem_alnreg_v *regs;
} worker_t;

static void *worker1(void *data)
{
	worker_t *w = (worker_t*)data;
	int i;
	for (i = w->start; i < w->n; i += w->step)
		w->regs[i] = find_alnreg(w->opt, w->bwt, w->bns, w->pac, &w->seqs[i]);
	return 0;
}

static void *worker2(void *data)
{
	worker_t *w = (worker_t*)data;
	int i;
	if (!w->opt->is_pe) {
		for (i = 0; i < w->n; i += w->step) {
			mem_sam_se(w->opt, w->bns, w->pac, &w->seqs[i], &w->regs[i]);
			free(w->regs[i].a);
		}
	} else {
		for (i = 0; i < w->n>>1; i += w->step) { // not implemented yet
			free(w->regs[i<<1|0].a); free(w->regs[i<<1|1].a);
		}
	}
	return 0;
}

int mem_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int n, bseq1_t *seqs)
{
	int i;
	worker_t *w;
	mem_alnreg_v *regs;
	w = calloc(opt->n_threads, sizeof(worker_t));
	regs = malloc(n * sizeof(mem_alnreg_v));
	for (i = 0; i < opt->n_threads; ++i) {
		worker_t *p = &w[i];
		p->start = i; p->step = opt->n_threads; p->n = n;
		p->opt = opt; p->bwt = bwt; p->bns = bns; p->pac = pac;
		p->seqs = seqs; p->regs = regs;
	}
#ifdef HAVE_PTHREAD
	if (opt->n_threads == 1) {
		worker1(w); worker2(w);
	} else {
		pthread_t *tid;
		tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
		for (i = 0; i < opt->n_threads; ++i) pthread_create(&tid[i], 0, worker1, &w[i]);
		for (i = 0; i < opt->n_threads; ++i) pthread_join(tid[i], 0);
		for (i = 0; i < opt->n_threads; ++i) pthread_create(&tid[i], 0, worker2, &w[i]);
		for (i = 0; i < opt->n_threads; ++i) pthread_join(tid[i], 0);
		free(tid);
	}
#else
	worker1(w); worker2(w);
#endif
	for (i = 0; i < n; ++i) {
		fputs(seqs[i].sam, stdout);
		free(seqs[i].name); free(seqs[i].comment); free(seqs[i].seq); free(seqs[i].qual); free(seqs[i].sam);
	}
	free(regs); free(w);
	return 0;
}
