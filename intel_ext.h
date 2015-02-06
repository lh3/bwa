#ifndef INTEL_EXT_H
#define INTEL_EXT_H

#include <stdint.h>
#include <stdbool.h>
#define INTEL_MIN_LEN  16
#define INTEL_MAX_LEN  255


#ifdef __cplusplus
extern "C" {
#endif

	void intel_init();
        void intel_destroy();
	void intel_filter(uint8_t *refSeq, int refLen, uint8_t *querySeq, int queryLen, int initScore, int endBonus);
        int (*intel_extend)(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off);

        int intel_filter_and_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off);
#ifdef __cplusplus
}
#endif


#endif
