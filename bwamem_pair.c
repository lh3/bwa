#include <stdlib.h>
#include <math.h>
#include "kstring.h"
#include "bwamem.h"
#include "kvec.h"

#define MIN_RATIO     0.8
#define MIN_DIR_CNT   10
#define MIN_DIR_RATIO 0.05
#define OUTLIER_BOUND 2.0
#define MAPPING_BOUND 3.0
#define MAX_STDDEV    4.0
#define EXT_STDDEV    4.0

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

typedef kvec_t(uint64_t) vec64_t;

void mem_pestat(const mem_opt_t *opt, int64_t l_pac, int n, const mem_alnreg_v *regs, mem_pestat_t pes[4])
{
	extern void ks_introsort_uint64_t(size_t n, uint64_t *a);
	int i, d, max;
	vec64_t isize[4];
	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	memset(isize, 0, sizeof(kvec_t(int)) * 4);
	for (i = 0; i < n>>1; ++i) {
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
		if (is <= opt->max_ins) kv_push(uint64_t, isize[dir], is);
	}
	if (mem_verbose >= 3) fprintf(stderr, "[M::%s] # candidate unique pairs for (FF, FR, RF, RR): (%ld, %ld, %ld, %ld)\n", __func__, isize[0].n, isize[1].n, isize[2].n, isize[3].n);
	for (d = 0; d < 4; ++d) { // TODO: this block is nearly identical to the one in bwtsw2_pair.c. It would be better to merge these two.
		mem_pestat_t *r = &pes[d];
		vec64_t *q = &isize[d];
		int p25, p50, p75, x;
		if (q->n < MIN_DIR_CNT) {
			fprintf(stderr, "[M::%s] skip orientation %c%c as there are not enough pairs\n", __func__, "FR"[d>>1&1], "FR"[d&1]);
			r->failed = 1;
			continue;
		} else fprintf(stderr, "[M::%s] analyzing insert size distribution for orientation %c%c...\n", __func__, "FR"[d>>1&1], "FR"[d&1]);
		ks_introsort_uint64_t(q->n, q->a);
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
		if (r->high < r->avg - MAX_STDDEV * r->std) r->high = (int)(r->avg + MAX_STDDEV * r->std + .499);
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
