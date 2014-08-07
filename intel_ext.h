#ifndef INTEL_EXT_H
#define INTEL_EXT_H

#include <stdint.h>

#define INTEL_MIN_LEN  16
#define INTEL_MAX_LEN  126

#ifdef __cplusplus
extern "C" {
#endif

	void intel_init();
	void intel_filter(uint8_t *refSeq, int refLen, uint8_t *querySeq, int queryLen, int initScore, int endBonus,
					  int *alignedQLen, int *alignedRLen, int *score, float *confidence);
	int intel_extend(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int gapo, int gape,
					 int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off);

#ifdef __cplusplus
}
#endif

#endif
