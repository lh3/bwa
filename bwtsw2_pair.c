#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "bwt.h"
#include "bntseq.h"
#include "bwtsw2.h"
#include "ksw.h"

#define MAX_INS       20000
#define MIN_RATIO     0.8
#define OUTLIER_BOUND 2.0
#define MAX_STDDEV    4.0

typedef struct {
	int low, high;
	double avg, std;
} bsw2pestat_t;

bsw2pestat_t bwtsw2_stat(int n, bwtsw2_t **buf)
{
	extern void ks_introsort_uint64_t(size_t n, uint64_t *a);
	int i, k, x, p25, p50, p75, tmp, max_len = 0;
	uint64_t *isize;
	bsw2pestat_t r;

	isize = calloc(n, 8);
	for (i = k = 0; i < n; i += 2) {
		bsw2hit_t *t[2];
		int l;
		if (buf[i]->n != 1 || buf[i+1]->n != 1) continue; // more than 1 hits
		t[0] = &buf[i]->hits[0]; t[1] = &buf[i+1]->hits[0];
		if (t[0]->G2 > 0.8 * t[0]->G) continue; // the best hit is not good enough
		if (t[1]->G2 > 0.8 * t[1]->G) continue; // the best hit is not good enough
		l = t[0]->k > t[1]->k? t[0]->k - t[1]->k + (t[1]->end - t[1]->beg) : t[1]->k - t[0]->k + (t[0]->end - t[0]->beg);
		max_len = max_len > t[0]->end - t[0]->beg? max_len : t[0]->end - t[0]->beg;
		max_len = max_len > t[1]->end - t[1]->beg? max_len : t[1]->end - t[1]->beg;
		isize[k++] = l;
	}
	ks_introsort_uint64_t(k, isize);
	p25 = isize[(int)(.25 * k + .499)];
	p50 = isize[(int)(.50 * k + .499)];
	p75 = isize[(int)(.75 * k + .499)];
	tmp    = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
	r.low  = tmp > max_len? tmp : max_len;
	r.high = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);
	fprintf(stderr, "[%s] (25, 50, 75) percentile: (%d, %d, %d)\n", __func__, p25, p50, p75);
	fprintf(stderr, "[%s] low and high boundaries for computing mean and std.dev: (%d, %d)\n", __func__, r.low, r.high);
	for (i = x = 0, r.avg = 0; i < k; ++i)
		if (isize[i] >= r.low && isize[i] <= r.high)
			r.avg += isize[i], ++x;
	r.avg /= x;
	for (i = 0, r.std = 0; i < k; ++i)
		if (isize[i] >= r.low && isize[i] <= r.high)
			r.std += (isize[i] - r.avg) * (isize[i] - r.avg);
	r.std = sqrt(r.std / x);
	fprintf(stderr, "[%s] mean and std.dev: (%.2f, %.2f)\n", __func__, r.avg, r.std);
	tmp  = (int)(p25 - 3. * (p75 - p25) + .499);
	r.low  = tmp > max_len? tmp : max_len;
	r.high = (int)(p75 + 3. * (p75 - p25) + .499);
	if (r.low > r.avg - MAX_STDDEV * 4.) r.low = (int)(r.avg - MAX_STDDEV * 4. + .499);
	r.low = tmp > max_len? tmp : max_len;
	if (r.high < r.avg - MAX_STDDEV * 4.) r.high = (int)(r.avg + MAX_STDDEV * 4. + .499);
	fprintf(stderr, "[%s] low and high boundaries for proper pairs: (%d, %d)\n", __func__, r.low, r.high);
	return r;
}

void bwtsw2_pair(const uint8_t *pac, int n, bsw2seq1_t *seq, bwtsw2_t **hits)
{
	bsw2pestat_t pes;
	pes = bwtsw2_stat(n, hits);
}
