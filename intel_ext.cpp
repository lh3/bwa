#include <stdio.h>
#include "intel_ext.h"

//Check if it is __x86_64 HW and supports SSE 4.2
#ifdef __x86_64
inline bool is_cpuid_ecx_bit_set(int eax, int bitidx) 
{ 
  int ecx = 0, edx = 0, ebx = 0; 
  __asm__ ("cpuid" 
	   :"=b" (ebx), 
	    "=c" (ecx), 
	    "=d" (edx) 
	   :"a" (eax) 
	   ); 
  return (((ecx >> bitidx)&1) == 1); 
}

inline bool is_sse42_supported() 
{ 
  return is_cpuid_ecx_bit_set(1, 20); 
}
#endif

extern "C" {
int ksw_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off);
}

// Needs to be called only once per program to set up static arrays
void init_ed_dist() ;

// The filtering function:
// Inputs: Reference and query sequences,	
//         Initial score (i.e. h0),
//         endBonus (i.e. the extra score if the whole query is aligned)
// Outputs:
//         Alignment length in query
//         Alignment length in reference
//         Alignment score
//         Confidence: For now, it's either 0.0 or 1.0, corresponding to no/full confidence in outputs
// Usage:
//         If confidence == 0.0: Partial alignment - need to rerun ksw_extend(...)
//         If confidence == 1.0
//               if alignedQLen == queryLen: Full alignment
//               if alignedQLen == 0: No alignment
// Notes: 
//        For now, the costMatrix and gap penalties are hardcoded. 
//
void extend_with_edit_dist(uint8_t* refSeq, int refLen, uint8_t* querySeq, int queryLen,
     			   int initScore, int endBonus, 
                           int& alignedQLen, int& alignedRLen, int& score, float& confidence) ;


// Filter-and-extend function:
// Inputs: Reference and query sequences,	
//         Initial score (i.e. h0),
//         endBonus (i.e. the extra score if the whole query is aligned)
//         zdrop value passed to ksw_extend
// Outputs:
//         Alignment length in query
//         Alignment length in reference
//         Alignment score
// Behavior:
//         The filtering function will be called internally first.
//         If there is an obvious result, it will be returned.
//         If not, ksw_extend() will be called with feedback from filtering function,
//         and its result will be returned.
// Notes: 
//        For now, the costMatrix and gap penalties are hardcoded. 
//        It is assumed that ksw_extend(...) function is linked. In this version,
//        it is defined in bwa_extend.cpp file.
void filter_and_extend(const uint8_t* refSeq, int refLen, const uint8_t* querySeq, int queryLen,
		       int initScore, int endBonus, int zdrop,
		       int& alignedQLen, int& alignedRLen, int& score) ;


void intel_init()
{
#ifdef __x86_64
    if (is_sse42_supported()) {
	fprintf(stderr,"Using Intel's filter_and_extend\n");
	intel_extend = intel_filter_and_extend;
	init_ed_dist();
    } else
	intel_extend = ksw_extend;
#else	
    intel_extend = ksw_extend;
#endif

}

void intel_filter(uint8_t *refSeq, int refLen, uint8_t *querySeq, int queryLen, int initScore, int endBonus,
				  int *alignedQLen, int *alignedRLen, int *score, float *confidence)
{
	if (queryLen < INTEL_MIN_LEN || queryLen > INTEL_MAX_LEN) {
		*confidence = 0.0;
		return;
	}
	extend_with_edit_dist(refSeq, refLen, querySeq, queryLen, initScore, endBonus, *alignedQLen, *alignedRLen, *score, *confidence);
}

int intel_filter_and_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off)
{
	if (qlen >= INTEL_MIN_LEN && qlen <= INTEL_MAX_LEN) {
		int score, a = mat[0];
		filter_and_extend(target, tlen, query, qlen, h0/a, end_bonus/a, zdrop/a, *qle, *tle, score);
		*max_off = 0; *gtle = *tle; *gscore = score;
		score *= a;
		/*
		int score2 = score, qle2 = *qle, tle2 = *tle;
		score = ksw_extend(qlen, query, tlen, target, m, mat, gapo, gape, w, end_bonus, zdrop, h0, qle, tle, gtle, gscore, max_off);
		fprintf(stderr, "[%d,%d,%d] %d:%d; %d:%d; %d:%d\n", qlen, tlen, h0, score, score2, *qle, qle2, *tle, tle2);
		*/
		return score;
	} else return ksw_extend(qlen, query, tlen, target, m, mat, gapo, gape, w, end_bonus, zdrop, h0, qle, tle, gtle, gscore, max_off);
}
