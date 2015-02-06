#include <stdio.h>
#include "intel_ext.h"

inline bool is_sse42_supported() 
{ 
    if(__builtin_cpu_supports("sse4.2"))
	return true;
    else
	return false;
}

inline bool is_avx2_supported()
{
    if(__builtin_cpu_supports("avx2"))
	return true;
    else
	return false;
}


extern "C" {
int ksw_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off);

}
extern const int defaultCostMatrixRowCnt ;
extern const int8_t defaultCostMatrix[];
void fast_filter_and_extend(const uint8_t* refSeq, int refLen, const uint8_t* querySeq, int queryLen,
			    int initScore, int endBonus, int zdrop,
			    int& alignedQLen, int& alignedRLen, int& score,
			    int costMatrixRowCnt=defaultCostMatrixRowCnt, 
			    const int8_t* costMatrix=defaultCostMatrix,
			    int gapWt=1, int gapOpenWt=6) ;
void fast_filter(const uint8_t* refSeq, int refLen, const uint8_t* querySeq, int queryLen,
                 int initScore, int endBonus,
                 int& alignedQLen, int& alignedRLen, int& score, float& confidence,
                 int mismatchWt=4, int gapWt=1, int gapOpenWt=6, int ambigWt=1) ;

void init_fast_extend();

void intel_init()
{
#ifdef __x86_64
    if (is_avx2_supported()) {
	intel_extend = intel_filter_and_extend;
	init_fast_extend();
	fprintf(stderr,"Initializing Intel's filter_and_extend with AVX2\n");
    } else if (is_sse42_supported()) {
	intel_extend = intel_filter_and_extend;
	init_fast_extend();
	fprintf(stderr,"Initializing Intel's filter_and_extend with SSE4.2\n");
    } else {
	intel_extend = ksw_extend;
    }
#else	
    intel_extend = ksw_extend;
#endif

}
void intel_destroy()
{
#ifdef __x86_64

#endif
}
void intel_filter(uint8_t* refSeq, int refLen, uint8_t* querySeq, int queryLen,
		  int initScore, int endBonus)
{
    int alignedQLen, alignedRLen, score;
    float confidence;
	if (queryLen < INTEL_MIN_LEN || queryLen > INTEL_MAX_LEN) {
		confidence = 0.0;
		return;
	}

	int mismatchWt = 4 ;
	int gapExtendWt = 1 ;
	int gapOpenWt = 6 ;
	int ambigWt = 1 ;
	fast_filter(refSeq, refLen, querySeq, queryLen, initScore, endBonus, alignedQLen, alignedRLen, score, confidence, mismatchWt, gapExtendWt, gapOpenWt, ambigWt);
}

int intel_filter_and_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off)
{
	if (qlen >= INTEL_MIN_LEN && qlen <= INTEL_MAX_LEN) {
		int score, a = mat[0];
		int gapExtendWt = 1 ;
		int gapOpenWt = 6 ;
		int ii = 0;
#ifdef DEBUG0
		for (ii=0;ii<tlen;ii++) fprintf(stderr,"%d", target[ii]);
		fprintf(stderr, " %d ", tlen);
		for (ii=0;ii<qlen;ii++) fprintf(stderr, "%d", query[ii]); 
		fprintf(stderr, " %d ", qlen);
		fprintf(stderr, "%d ", h0/a);
		fprintf(stderr, "%d ", end_bonus/a);
		fprintf(stderr, "%d \n", zdrop/a);
#endif

		fast_filter_and_extend(target, tlen, query, qlen, h0/a, end_bonus/a, zdrop/a, *qle, *tle, score, defaultCostMatrixRowCnt, defaultCostMatrix, gapExtendWt, gapOpenWt);
		// Convert it to the form ksw_extend uses:
		
		if (*qle == qlen) { // global score is chosen over local score
		    *gscore = score; 
		    *gtle = *tle;
		} else {
		    *gscore = 0;
		    *gtle = 0;
		}

		*max_off = 0; //*gtle = *tle; *gscore = score;
		score *= a;
		/*
		int score2 = score, qle2 = *qle, tle2 = *tle;
		score = ksw_extend(qlen, query, tlen, target, m, mat, gapo, gape, w, end_bonus, zdrop, h0, qle, tle, gtle, gscore, max_off);
		fprintf(stderr, "[%d,%d,%d] %d:%d; %d:%d; %d:%d\n", qlen, tlen, h0, score, score2, *qle, qle2, *tle, tle2);
		*/
		return score;
	} else {
	    return ksw_extend(qlen, query, tlen, target, m, mat, gapo, gape, w, end_bonus, zdrop, h0, qle, tle, gtle, gscore, max_off);}
}
