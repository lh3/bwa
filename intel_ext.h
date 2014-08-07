#ifndef INTEL_EXT_H
#define INTEL_EXT_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

	void intel_init();
	void intel_filter(uint8_t *refSeq, int refLen, uint8_t *querySeq, int queryLen, int initScore, int endBonus,
					  int *alignedQLen, int *alignedRLen, int *score, float *confidence);

#ifdef __cplusplus
}
#endif

#endif
