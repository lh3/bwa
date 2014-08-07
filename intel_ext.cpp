#include "intel_ext.h"

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
void filter_and_extend(uint8_t* refSeq, int refLen, uint8_t* querySeq, int queryLen,
		       int initScore, int endBonus, int zdrop,
		       int& alignedQLen, int& alignedRLen, int& score) ;


void intel_init()
{	
	init_ed_dist();
}

void intel_filter(uint8_t *refSeq, int refLen, uint8_t *querySeq, int queryLen, int initScore, int endBonus,
				  int *alignedQLen, int *alignedRLen, int *score, float *confidence)
{
	extend_with_edit_dist(refSeq, refLen, querySeq, queryLen, initScore, endBonus, *alignedQLen, *alignedRLen, *score, *confidence);
}
