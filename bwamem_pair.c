/* The MIT License

   Copyright (c) 2018-     Dana-Farber Cancer Institute
                 2009-2018 Broad Institute, Inc.
                 2008-2009 Genome Research Ltd. (GRL)

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "kstring.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
#include "ksw.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif


#define MIN_RATIO     0.8
#define MIN_DIR_CNT   10
#define MIN_DIR_RATIO 0.05
#define OUTLIER_BOUND 2.0
#define MAPPING_BOUND 3.0
#define MAX_STDDEV    4.0

static inline int mem_infer_dir(int64_t l_pac, int64_t b1, int64_t b2, int64_t *dist)
{
	int64_t p2;
	int r1 = (b1 >= l_pac), r2 = (b2 >= l_pac);
	p2 = r1 == r2? b2 : (l_pac<<1) - 1 - b2; // p2 is the coordinate of read 2 on the read 1 strand
	*dist = p2 > b1? p2 - b1 : b1 - p2;
	return (r1 == r2? 0 : 1) ^ (p2 > b1? 0 : 3);
}

static int cal_sub(const mem_opt_t *opt, mem_alnreg_v *r)
{
	int j;
	for (j = 1; j < r->n; ++j) { // choose unique alignment
		int b_max = r->a[j].qb > r->a[0].qb? r->a[j].qb : r->a[0].qb;
		int e_min = r->a[j].qe < r->a[0].qe? r->a[j].qe : r->a[0].qe;
		if (e_min > b_max) { // have overlap
			int min_l = r->a[j].qe - r->a[j].qb < r->a[0].qe - r->a[0].qb? r->a[j].qe - r->a[j].qb : r->a[0].qe - r->a[0].qb;
			if (e_min - b_max >= min_l * opt->mask_level) break; // significant overlap
		}
	}
	return j < r->n? r->a[j].score : opt->min_seed_len * opt->a;
}

void mem_pestat(const mem_opt_t *opt, int64_t l_pac, int n, const mem_alnreg_v *regs, mem_pestat_t pes[4])
{
	int i, d, max;
	uint64_v isize[4];
	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	memset(isize, 0, sizeof(kvec_t(int)) * 4);
	for (i = 0; i < n>>1; ++i) {
		int dir;
		int64_t is;
		mem_alnreg_v *r[2];
		r[0] = (mem_alnreg_v*)&regs[i<<1|0];
		r[1] = (mem_alnreg_v*)&regs[i<<1|1];
		if (r[0]->n == 0 || r[1]->n == 0) continue;
		if (cal_sub(opt, r[0]) > MIN_RATIO * r[0]->a[0].score) continue;
		if (cal_sub(opt, r[1]) > MIN_RATIO * r[1]->a[0].score) continue;
		if (r[0]->a[0].rid != r[1]->a[0].rid) continue; // not on the same chr
		dir = mem_infer_dir(l_pac, r[0]->a[0].rb, r[1]->a[0].rb, &is);
		if (is && is <= opt->max_ins) kv_push(uint64_t, isize[dir], is);
	}
	if (bwa_verbose >= 3) fprintf(stderr, "[M::%s] # candidate unique pairs for (FF, FR, RF, RR): (%ld, %ld, %ld, %ld)\n", __func__, isize[0].n, isize[1].n, isize[2].n, isize[3].n);
	for (d = 0; d < 4; ++d) { // TODO: this block is nearly identical to the one in bwtsw2_pair.c. It would be better to merge these two.
		mem_pestat_t *r = &pes[d];
		uint64_v *q = &isize[d];
		int p25, p50, p75, x;
		if (q->n < MIN_DIR_CNT) {
			fprintf(stderr, "[M::%s] skip orientation %c%c as there are not enough pairs\n", __func__, "FR"[d>>1&1], "FR"[d&1]);
			r->failed = 1;
			free(q->a);
			continue;
		} else fprintf(stderr, "[M::%s] analyzing insert size distribution for orientation %c%c...\n", __func__, "FR"[d>>1&1], "FR"[d&1]);
		ks_introsort_64(q->n, q->a);
		p25 = q->a[(int)(.25 * q->n + .499)];
		p50 = q->a[(int)(.50 * q->n + .499)];
		p75 = q->a[(int)(.75 * q->n + .499)];
		r->low  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
		if (r->low < 1) r->low = 1;
		r->high = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);
		fprintf(stderr, "[M::%s] (25, 50, 75) percentile: (%d, %d, %d)\n", __func__, p25, p50, p75);
		fprintf(stderr, "[M::%s] low and high boundaries for computing mean and std.dev: (%d, %d)\n", __func__, r->low, r->high);
		for (i = x = 0, r->avg = 0; i < q->n; ++i)
			if (q->a[i] >= r->low && q->a[i] <= r->high)
				r->avg += q->a[i], ++x;
		r->avg /= x;
		for (i = 0, r->std = 0; i < q->n; ++i)
			if (q->a[i] >= r->low && q->a[i] <= r->high)
				r->std += (q->a[i] - r->avg) * (q->a[i] - r->avg);
		r->std = sqrt(r->std / x);
		fprintf(stderr, "[M::%s] mean and std.dev: (%.2f, %.2f)\n", __func__, r->avg, r->std);
		r->low  = (int)(p25 - MAPPING_BOUND * (p75 - p25) + .499);
		r->high = (int)(p75 + MAPPING_BOUND * (p75 - p25) + .499);
		if (r->low  > r->avg - MAX_STDDEV * r->std) r->low  = (int)(r->avg - MAX_STDDEV * r->std + .499);
		if (r->high < r->avg + MAX_STDDEV * r->std) r->high = (int)(r->avg + MAX_STDDEV * r->std + .499);
		if (r->low < 1) r->low = 1;
		fprintf(stderr, "[M::%s] low and high boundaries for proper pairs: (%d, %d)\n", __func__, r->low, r->high);
		free(q->a);
	}
	for (d = 0, max = 0; d < 4; ++d)
		max = max > isize[d].n? max : isize[d].n;
	for (d = 0; d < 4; ++d)
		if (pes[d].failed == 0 && isize[d].n < max * MIN_DIR_RATIO) {
			pes[d].failed = 1;
			fprintf(stderr, "[M::%s] skip orientation %c%c\n", __func__, "FR"[d>>1&1], "FR"[d&1]);
		}
}

int mem_matesw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], const mem_alnreg_t *a, int l_ms, const uint8_t *ms, mem_alnreg_v *ma)
{
	extern int mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, int n, mem_alnreg_t *a);
	int64_t l_pac = bns->l_pac;
	int i, r, skip[4], n = 0, rid;
	for (r = 0; r < 4; ++r)
		skip[r] = pes[r].failed? 1 : 0;
	for (i = 0; i < ma->n; ++i) { // check which orinentation has been found
		int64_t dist;
		r = mem_infer_dir(l_pac, a->rb, ma->a[i].rb, &dist);
		if (dist >= pes[r].low && dist <= pes[r].high)
			skip[r] = 1;
	}
	if (skip[0] + skip[1] + skip[2] + skip[3] == 4) return 0; // consistent pair exist; no need to perform SW
	for (r = 0; r < 4; ++r) {
		int is_rev, is_larger;
		uint8_t *seq, *rev = 0, *ref = 0;
		int64_t rb, re;
		if (skip[r]) continue;
		is_rev = (r>>1 != (r&1)); // whether to reverse complement the mate
		is_larger = !(r>>1); // whether the mate has larger coordinate
		if (is_rev) {
			rev = malloc(l_ms); // this is the reverse complement of $ms
			for (i = 0; i < l_ms; ++i) rev[l_ms - 1 - i] = ms[i] < 4? 3 - ms[i] : 4;
			seq = rev;
		} else seq = (uint8_t*)ms;
		if (!is_rev) {
			rb = is_larger? a->rb + pes[r].low : a->rb - pes[r].high;
			re = (is_larger? a->rb + pes[r].high: a->rb - pes[r].low) + l_ms; // if on the same strand, end position should be larger to make room for the seq length
		} else {
			rb = (is_larger? a->rb + pes[r].low : a->rb - pes[r].high) - l_ms; // similarly on opposite strands
			re = is_larger? a->rb + pes[r].high: a->rb - pes[r].low;
		}
		if (rb < 0) rb = 0;
		if (re > l_pac<<1) re = l_pac<<1;
		if (rb < re) ref = bns_fetch_seq(bns, pac, &rb, (rb+re)>>1, &re, &rid);
		if (a->rid == rid && re - rb >= opt->min_seed_len) { // no funny things happening
			kswr_t aln;
			mem_alnreg_t b;
			int tmp, xtra = KSW_XSUBO | KSW_XSTART | (l_ms * opt->a < 250? KSW_XBYTE : 0) | (opt->min_seed_len * opt->a);
			aln = ksw_align2(l_ms, seq, re - rb, ref, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, xtra, 0);
			memset(&b, 0, sizeof(mem_alnreg_t));
			if (aln.score >= opt->min_seed_len && aln.qb >= 0) { // something goes wrong if aln.qb < 0
				b.rid = a->rid;
				b.is_alt = a->is_alt;
				b.qb = is_rev? l_ms - (aln.qe + 1) : aln.qb;                                                                                                                                                                              
				b.qe = is_rev? l_ms - aln.qb : aln.qe + 1; 
				b.rb = is_rev? (l_pac<<1) - (rb + aln.te + 1) : rb + aln.tb;
				b.re = is_rev? (l_pac<<1) - (rb + aln.tb) : rb + aln.te + 1;
				b.score = aln.score;
				b.csub = aln.score2;
				b.secondary = -1;
				b.seedcov = (b.re - b.rb < b.qe - b.qb? b.re - b.rb : b.qe - b.qb) >> 1;
//				printf("*** %d, [%lld,%lld], %d:%d, (%lld,%lld), (%lld,%lld) == (%lld,%lld)\n", aln.score, rb, re, is_rev, is_larger, a->rb, a->re, ma->a[0].rb, ma->a[0].re, b.rb, b.re);
				kv_push(mem_alnreg_t, *ma, b); // make room for a new element
				// move b s.t. ma is sorted
				for (i = 0; i < ma->n - 1; ++i) // find the insertion point
					if (ma->a[i].score < b.score) break;
				tmp = i;
				for (i = ma->n - 1; i > tmp; --i) ma->a[i] = ma->a[i-1];
				ma->a[i] = b;
			}
			++n;
		}
		if (n) ma->n = mem_sort_dedup_patch(opt, 0, 0, 0, ma->n, ma->a);
		if (rev) free(rev);
		free(ref);
	}
	return n;
}

int mem_pair(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], bseq1_t s[2], mem_alnreg_v a[2], int id, int *sub, int *n_sub, int z[2], int n_pri[2])
{
	pair64_v v, u;
	int r, i, k, y[4], ret; // y[] keeps the last hit
	int64_t l_pac = bns->l_pac;
	kv_init(v); kv_init(u);
	for (r = 0; r < 2; ++r) { // loop through read number
		for (i = 0; i < n_pri[r]; ++i) {
			pair64_t key;
			mem_alnreg_t *e = &a[r].a[i];
			key.x = e->rb < l_pac? e->rb : (l_pac<<1) - 1 - e->rb; // forward position
			key.x = (uint64_t)e->rid<<32 | (key.x - bns->anns[e->rid].offset);
			key.y = (uint64_t)e->score << 32 | i << 2 | (e->rb >= l_pac)<<1 | r;
			kv_push(pair64_t, v, key);
		}
	}
	ks_introsort_128(v.n, v.a);
	y[0] = y[1] = y[2] = y[3] = -1;
	//for (i = 0; i < v.n; ++i) printf("[%d]\t%d\t%c%ld\n", i, (int)(v.a[i].y&1)+1, "+-"[v.a[i].y>>1&1], (long)v.a[i].x);
	for (i = 0; i < v.n; ++i) {
		for (r = 0; r < 2; ++r) { // loop through direction
			int dir = r<<1 | (v.a[i].y>>1&1), which;
			if (pes[dir].failed) continue; // invalid orientation
			which = r<<1 | ((v.a[i].y&1)^1);
			if (y[which] < 0) continue; // no previous hits
			for (k = y[which]; k >= 0; --k) { // TODO: this is a O(n^2) solution in the worst case; remember to check if this loop takes a lot of time (I doubt)
				int64_t dist;
				int q;
				double ns;
				pair64_t *p;
				if ((v.a[k].y&3) != which) continue;
				dist = (int64_t)v.a[i].x - v.a[k].x;
				//printf("%d: %lld\n", k, dist);
				if (dist > pes[dir].high) break;
				if (dist < pes[dir].low)  continue;
				ns = (dist - pes[dir].avg) / pes[dir].std;
				q = (int)((v.a[i].y>>32) + (v.a[k].y>>32) + .721 * log(2. * erfc(fabs(ns) * M_SQRT1_2)) * opt->a + .499); // .721 = 1/log(4)
				if (q < 0) q = 0;
				p = kv_pushp(pair64_t, u);
				p->y = (uint64_t)k<<32 | i;
				p->x = (uint64_t)q<<32 | (hash_64(p->y ^ id<<8) & 0xffffffffU);
				//printf("[%lld,%lld]\t%d\tdist=%ld\n", v.a[k].x, v.a[i].x, q, (long)dist);
			}
		}
		y[v.a[i].y&3] = i;
	}
	if (u.n) { // found at least one proper pair
		int tmp = opt->a + opt->b;
		tmp = tmp > opt->o_del + opt->e_del? tmp : opt->o_del + opt->e_del;
		tmp = tmp > opt->o_ins + opt->e_ins? tmp : opt->o_ins + opt->e_ins;
		ks_introsort_128(u.n, u.a);
		i = u.a[u.n-1].y >> 32; k = u.a[u.n-1].y << 32 >> 32;
		z[v.a[i].y&1] = v.a[i].y<<32>>34; // index of the best pair
		z[v.a[k].y&1] = v.a[k].y<<32>>34;
		ret = u.a[u.n-1].x >> 32;
		*sub = u.n > 1? u.a[u.n-2].x>>32 : 0;
		for (i = (long)u.n - 2, *n_sub = 0; i >= 0; --i)
			if (*sub - (int)(u.a[i].x>>32) <= tmp) ++*n_sub;
	} else ret = 0, *sub = 0, *n_sub = 0;
	free(u.a); free(v.a);
	return ret;
}

void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m);
void mem_reorder_primary5(int T, mem_alnreg_v *a);

#define raw_mapq(diff, a) ((int)(6.02 * (diff) / (a) + .499))

int mem_sam_pe(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_pestat_t pes[4], uint64_t id, bseq1_t s[2], mem_alnreg_v a[2])
{
	extern int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id);
	extern int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a);
	extern void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m);
	extern char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_alnreg_v *a, int l_query, const char *query);

	int n = 0, i, j, z[2], o, subo, n_sub, extra_flag = 1, n_pri[2], n_aa[2];
	kstring_t str;
	mem_aln_t h[2], g[2], aa[2][2];

	str.l = str.m = 0; str.s = 0;
	memset(h, 0, sizeof(mem_aln_t) * 2);
	memset(g, 0, sizeof(mem_aln_t) * 2);
	n_aa[0] = n_aa[1] = 0;
	if (!(opt->flag & MEM_F_NO_RESCUE)) { // then perform SW for the best alignment
		mem_alnreg_v b[2];
		kv_init(b[0]); kv_init(b[1]);
		for (i = 0; i < 2; ++i)
			for (j = 0; j < a[i].n; ++j)
				if (a[i].a[j].score >= a[i].a[0].score  - opt->pen_unpaired)
					kv_push(mem_alnreg_t, b[i], a[i].a[j]);
		for (i = 0; i < 2; ++i)
			for (j = 0; j < b[i].n && j < opt->max_matesw; ++j)
				n += mem_matesw(opt, bns, pac, pes, &b[i].a[j], s[!i].l_seq, (uint8_t*)s[!i].seq, &a[!i]);
		free(b[0].a); free(b[1].a);
	}
	n_pri[0] = mem_mark_primary_se(opt, a[0].n, a[0].a, id<<1|0);
	n_pri[1] = mem_mark_primary_se(opt, a[1].n, a[1].a, id<<1|1);
	if (opt->flag & MEM_F_PRIMARY5) {
		mem_reorder_primary5(opt->T, &a[0]);
		mem_reorder_primary5(opt->T, &a[1]);
	}
	if (opt->flag&MEM_F_NOPAIRING) goto no_pairing;
	// pairing single-end hits
	if (n_pri[0] && n_pri[1] && (o = mem_pair(opt, bns, pac, pes, s, a, id, &subo, &n_sub, z, n_pri)) > 0) {
		int is_multi[2], q_pe, score_un, q_se[2];
		char **XA[2];
		// check if an end has multiple hits even after mate-SW
		for (i = 0; i < 2; ++i) {
			for (j = 1; j < n_pri[i]; ++j)
				if (a[i].a[j].secondary < 0 && a[i].a[j].score >= opt->T) break;
			is_multi[i] = j < n_pri[i]? 1 : 0;
		}
		if (is_multi[0] || is_multi[1]) goto no_pairing; // TODO: in rare cases, the true hit may be long but with low score
		// compute mapQ for the best SE hit
		score_un = a[0].a[0].score + a[1].a[0].score - opt->pen_unpaired;
		//q_pe = o && subo < o? (int)(MEM_MAPQ_COEF * (1. - (double)subo / o) * log(a[0].a[z[0]].seedcov + a[1].a[z[1]].seedcov) + .499) : 0;
		subo = subo > score_un? subo : score_un;
		q_pe = raw_mapq(o - subo, opt->a);
		if (n_sub > 0) q_pe -= (int)(4.343 * log(n_sub+1) + .499);
		if (q_pe < 0) q_pe = 0;
		if (q_pe > 60) q_pe = 60;
		q_pe = (int)(q_pe * (1. - .5 * (a[0].a[0].frac_rep + a[1].a[0].frac_rep)) + .499);
		// the following assumes no split hits
		if (o > score_un) { // paired alignment is preferred
			mem_alnreg_t *c[2];
			c[0] = &a[0].a[z[0]]; c[1] = &a[1].a[z[1]];
			for (i = 0; i < 2; ++i) {
				if (c[i]->secondary >= 0)
					c[i]->sub = a[i].a[c[i]->secondary].score, c[i]->secondary = -2;
				q_se[i] = mem_approx_mapq_se(opt, c[i]);
			}
			q_se[0] = q_se[0] > q_pe? q_se[0] : q_pe < q_se[0] + 40? q_pe : q_se[0] + 40;
			q_se[1] = q_se[1] > q_pe? q_se[1] : q_pe < q_se[1] + 40? q_pe : q_se[1] + 40;
			extra_flag |= 2;
			// cap at the tandem repeat score
			q_se[0] = q_se[0] < raw_mapq(c[0]->score - c[0]->csub, opt->a)? q_se[0] : raw_mapq(c[0]->score - c[0]->csub, opt->a);
			q_se[1] = q_se[1] < raw_mapq(c[1]->score - c[1]->csub, opt->a)? q_se[1] : raw_mapq(c[1]->score - c[1]->csub, opt->a);
		} else { // the unpaired alignment is preferred
			z[0] = z[1] = 0;
			q_se[0] = mem_approx_mapq_se(opt, &a[0].a[0]);
			q_se[1] = mem_approx_mapq_se(opt, &a[1].a[0]);
		}
		for (i = 0; i < 2; ++i) {
			int k = a[i].a[z[i]].secondary_all;
			if (k >= 0 && k < n_pri[i]) { // switch secondary and primary if both of them are non-ALT
				assert(a[i].a[k].secondary_all < 0);
				for (j = 0; j < a[i].n; ++j)
					if (a[i].a[j].secondary_all == k || j == k)
						a[i].a[j].secondary_all = z[i];
				a[i].a[z[i]].secondary_all = -1;
			}
		}
		if (!(opt->flag & MEM_F_ALL)) {
			for (i = 0; i < 2; ++i)
				XA[i] = mem_gen_alt(opt, bns, pac, &a[i], s[i].l_seq, s[i].seq);
		} else XA[0] = XA[1] = 0;
		// write SAM
		for (i = 0; i < 2; ++i) {
			h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, &a[i].a[z[i]]);
			h[i].mapq = q_se[i];
			h[i].flag |= 0x40<<i | extra_flag;
			h[i].XA = XA[i]? XA[i][z[i]] : 0;
			aa[i][n_aa[i]++] = h[i];
			if (n_pri[i] < a[i].n) { // the read has ALT hits
				mem_alnreg_t *p = &a[i].a[n_pri[i]];
				if (p->score < opt->T || p->secondary >= 0 || !p->is_alt) continue;
				g[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, p);
				g[i].flag |= 0x800 | 0x40<<i | extra_flag;
				g[i].XA = XA[i]? XA[i][n_pri[i]] : 0;
				aa[i][n_aa[i]++] = g[i];
			}
		}
		for (i = 0; i < n_aa[0]; ++i)
			mem_aln2sam(opt, bns, &str, &s[0], n_aa[0], aa[0], i, &h[1]); // write read1 hits
		s[0].sam = strdup(str.s); str.l = 0;
		for (i = 0; i < n_aa[1]; ++i)
			mem_aln2sam(opt, bns, &str, &s[1], n_aa[1], aa[1], i, &h[0]); // write read2 hits
		s[1].sam = str.s;
		if (strcmp(s[0].name, s[1].name) != 0) err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n", s[0].name, s[1].name);
		// free
		for (i = 0; i < 2; ++i) {
			free(h[i].cigar); free(g[i].cigar);
			if (XA[i] == 0) continue;
			for (j = 0; j < a[i].n; ++j) free(XA[i][j]);
			free(XA[i]);
		}
	} else goto no_pairing;
	return n;

no_pairing:
	for (i = 0; i < 2; ++i) {
		int which = -1;
		if (a[i].n) {
			if (a[i].a[0].score >= opt->T) which = 0;
			else if (n_pri[i] < a[i].n && a[i].a[n_pri[i]].score >= opt->T)
				which = n_pri[i];
		}
		if (which >= 0) h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, &a[i].a[which]);
		else h[i] = mem_reg2aln(opt, bns, pac, s[i].l_seq, s[i].seq, 0);
	}
	if (!(opt->flag & MEM_F_NOPAIRING) && h[0].rid == h[1].rid && h[0].rid >= 0) { // if the top hits from the two ends constitute a proper pair, flag it.
		int64_t dist;
		int d;
		d = mem_infer_dir(bns->l_pac, a[0].a[0].rb, a[1].a[0].rb, &dist);
		if (!pes[d].failed && dist >= pes[d].low && dist <= pes[d].high) extra_flag |= 2;
	}
	mem_reg2sam(opt, bns, pac, &s[0], &a[0], 0x41|extra_flag, &h[1]);
	mem_reg2sam(opt, bns, pac, &s[1], &a[1], 0x81|extra_flag, &h[0]);
	if (strcmp(s[0].name, s[1].name) != 0) err_fatal(__func__, "paired reads have different names: \"%s\", \"%s\"\n", s[0].name, s[1].name);
	free(h[0].cigar); free(h[1].cigar);
	return n;
}
