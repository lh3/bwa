#include <stdio.h>
#include "bwtsw2.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

typedef struct {
	uint32_t tbeg, tend;
	int qbeg, qend;
	uint32_t flag:1, idx:31;
	int chain; // also reuse as a counter
} hsaip_t;

#define _hsaip_lt(a, b) ((a).qbeg < (b).qbeg)

#include "ksort.h"
KSORT_INIT(hsaip, hsaip_t, _hsaip_lt)

static int chaining(const bsw2opt_t *opt, int shift, int n, hsaip_t *z, hsaip_t *chain)
{
	int j, k, m = 0;
	ks_introsort(hsaip, n, z);
	for (j = 0; j < n; ++j) {
		hsaip_t *p = z + j;
		for (k = m - 1; k >= 0; --k) {
			hsaip_t *q = chain + k;
			int x = p->qbeg - q->qbeg; // always positive
			int y = p->tbeg - q->tbeg;
			if (y > 0 && x < opt->max_chain_gap && y < opt->max_chain_gap && x - y <= opt->bw && y - x <= opt->bw) { // chained
				if (p->qend > q->qend) q->qend = p->qend;
				if (p->tend > q->tend) q->tend = p->tend;
				++q->chain;
				p->chain = shift + k;
				break;
			} else if (q->chain > opt->t_seeds * 2) k = 0; // if the chain is strong enough, do not check the previous chains
		}
		if (k < 0) { // not added to any previous chains
			chain[m] = *p;
			chain[m].chain = 1;
			chain[m].idx = p->chain = shift + m;
			++m;
		}
	}
	return m;
}

void bsw2_chain_filter(const bsw2opt_t *opt, int len, bwtsw2_t *b[2])
{
	hsaip_t *z[2], *chain[2];
	int i, j, k, n[2], m[2], thres = opt->t_seeds * 2;
	char *flag;
	// initialization
	n[0] = b[0]->n; n[1] = b[1]->n;
	z[0] = calloc(n[0] + n[1], sizeof(hsaip_t));
	z[1] = z[0] + n[0];
	chain[0] = calloc(n[0] + n[1], sizeof(hsaip_t));
	for (k = j = 0; k < 2; ++k) {
		for (i = 0; i < b[k]->n; ++i) {
			bsw2hit_t *p = b[k]->hits + i;
			hsaip_t *q = z[k] + i;
			q->flag = k; q->idx = i;
			q->tbeg = p->k; q->tend = p->k + p->len;
			q->chain = -1;
			q->qbeg = p->beg; q->qend = p->end;
		}
	}
	// chaining
	m[0] = chaining(opt, 0,    n[0], z[0], chain[0]);
	chain[1] = chain[0] + m[0];
	m[1] = chaining(opt, m[0], n[1], z[1], chain[1]);	
	// change query coordinate on the reverse strand
	for (k = 0; k < m[1]; ++k) {
		hsaip_t *p = chain[1] + k;
		int tmp = p->qbeg;
		p->qbeg = len - p->qend; p->qend = len - tmp;
	}
	//for (k = 0; k < m[0]; ++k) printf("%d, [%d,%d), [%d,%d)\n", chain[0][k].chain, chain[0][k].tbeg, chain[0][k].tend, chain[0][k].qbeg, chain[0][k].qend);
	// filtering
	flag = calloc(m[0] + m[1], 1);
	ks_introsort(hsaip, m[0] + m[1], chain[0]);
	for (k = 1; k < m[0] + m[1]; ++k) {
		hsaip_t *p = chain[0] + k;
		for (j = 0; j < k; ++j) {
			hsaip_t *q = chain[0] + j;
			if (flag[q->idx]) continue;
			if (q->qend >= p->qend && q->chain > p->chain * thres && p->chain < thres) {
				flag[p->idx] = 1;
				break;
			}
		}
	}
	for (k = 0; k < n[0] + n[1]; ++k) {
		hsaip_t *p = z[0] + k;
		if (flag[p->chain])
			b[p->flag]->hits[p->idx].G = 0;
	}
	free(flag);
	// squeeze out filtered elements in b[2]
	for (k = 0; k < 2; ++k) {
		for (j = i = 0; j < n[k]; ++j) {
			bsw2hit_t *p = b[k]->hits + j;
			if (p->G) {
				if (i != j) b[k]->hits[i++] = *p;
				else ++i;
			}
		}
		b[k]->n = i;
	}
	// free
	free(z[0]); free(chain[0]);
}
