#include <stdlib.h>
#include "kstring.h"
#include "bwamem.h"
#include "kvec.h"

#define MIN_RATIO     0.8

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
	int i;
	kvec_t(int) isize[4];
	memset(isize, 0, sizeof(kvec_t(int)) * 4);
	for (i = 0; i < n>>1; i += 2) {
		int dir;
		int64_t is, pos[2];
		mem_alnreg_v *r[2];
		r[0] = (mem_alnreg_v*)&regs[i<<1|0];
		r[1] = (mem_alnreg_v*)&regs[i<<1|1];
		if (r[0]->n == 0 || r[1]->n == 0) continue;
		if (cal_sub(opt, r[0]) > MIN_RATIO * r[0]->a[0].score) continue;
		if (cal_sub(opt, r[1]) > MIN_RATIO * r[1]->a[0].score) continue;
		pos[0] = r[0]->a[0].rb < l_pac? r[0]->a[0].rb : (l_pac<<1) - 1 - r[0]->a[0].rb; // forward coordinate
		pos[1] = r[1]->a[0].rb < l_pac? r[1]->a[0].rb : (l_pac<<1) - 1 - r[1]->a[0].rb;
		if (pos[0] < pos[1]) dir = (r[0]->a[0].rb >= l_pac)<<1 | (r[1]->a[0].rb >= l_pac);
		else dir = (r[1]->a[0].rb >= l_pac)<<1 | (r[0]->a[0].rb >= l_pac);
		is = abs(pos[0] - pos[1]);
		if (is <= opt->max_ins) kv_push(int, isize[dir], is);
	}
	if (mem_verbose >= 3) fprintf(stderr, "[M::%s] # candidates unique pairs for (FF, FR, RF, RR): (%ld, %ld, %ld, %ld)\n", __func__, isize[0].n, isize[1].n, isize[2].n, isize[3].n);
}
