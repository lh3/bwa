/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <iostream>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#ifdef USE_AVX2
#include "fast_extend_vec256.h"
#include "fast_extend_bitv256.h"
#endif

#include "fast_extend_vec128x2.h"
#include "fast_extend_vec128.h"
#include "fast_extend_vec64.h"
#include "fast_extend_vec32.h"
#include "fast_extend_bitv64.h"
#include "fast_extend_bitv64x2.h"
#include "fast_extend_bitv128.h"
#include "fast_extend_bitv128x2.h"
#include "fast_extend_fine.h"
#include "fast_extend_dist.h"
#include "fast_extend.h"

// Perform fine-grain alignment because it improves quality of cases with very short extension
#define FINE_ALIGN

//#define DEBUG0
//#define DEBUG1
//#define DEBUG2

//#ifdef SW_FILTER_AND_EXTEND
//#define SW_FEEDBACK 
//#endif
#define SW_FILTER_AND_EXTEND
#define SW_FEEDBACK

extern "C" {
int ksw_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off);
}

const int defaultCostMatrixRowCnt = 5 ;

const int8_t defaultCostMatrix[] = {
  1, -4, -4, -4, -1,
  -4,  1, -4, -4, -1,
  -4, -4,  1, -4, -1,
  -4, -4, -4,  1, -1,
  -1, -1, -1, -1, -1 
} ;

struct SWFeedback {
  int maxQLen ;
  int maxRLen ;
  int maxBand ;
  int minQLen ;
  int minRLen ;
  int minScore ;

  SWFeedback():maxQLen(0), maxRLen(0), maxBand(0), minQLen(0), minRLen(0), minScore(0) {} ;

  SWFeedback(int qlen, int rlen, int band, int origQLen, int origRLen)
    : maxQLen(qlen), maxRLen(rlen), maxBand(band), 
      minQLen(origQLen), minRLen(origRLen), minScore(0) {}

  void updateMax(int qlen, int rlen, int band) {
    maxQLen = qlen ;
    maxRLen = rlen ;
    maxBand = band ;
  }

  void updateMin(int currQLen, int currRLen, int currScore) {
    minQLen = currQLen ;
    minRLen = currRLen ;
    minScore = currScore ;
  }
} ;

inline ostream& operator<<(ostream& os, const SWFeedback& swfb) {
  os << "{"  << swfb.minQLen << ":" << swfb.maxQLen << ", " << swfb.minRLen << ":"
     << swfb.maxRLen << " " << swfb.maxBand << " " << swfb.minScore << "}" << endl ;
  return os ;
}

struct PartialAlignment {
  int partialQLen, partialRLen ;
  int scoreLow, scoreHigh ;
#ifdef SW_FEEDBACK
  int estimatedBand ;
#endif
  
  PartialAlignment() {}
  PartialAlignment(int qlen, int rlen, int scoreL, int scoreH)
    : partialQLen(qlen), partialRLen(rlen), scoreLow(scoreL), scoreHigh(scoreH) {
#ifdef SW_FEEDBACK
    estimatedBand = 0 ;
#endif
  }

} ;

inline ostream& operator<<(ostream& os, const PartialAlignment& align) {
  cerr << "Partial alignment (" << align.partialQLen << ", " << align.partialRLen 
       << "): score = [" << align.scoreLow << ", " << align.scoreHigh << "]"
#ifdef SW_FEEDBACK
       << ", estimatedBand = " << align.estimatedBand
#endif
    ;
  return os ;
}


// Calls the original ksw_extend function
int run_ksw_extend(const uint8_t* refSeq, int refLen, const uint8_t* querySeq, int queryLen,
		   int bandW, int initScore, int endBonus, int zdrop, int costMatrixRowCnt,
		   const int8_t* costMatrix, int gapo, int gape,
		   int& alignedQLen, int& alignedRLen) {

  int qle, tle, gtle, gscore, max_off ;

  int pscore = 0 ;

  if (bandW < 8) 
    pscore = ksw_extend_sse_8(queryLen, querySeq, refLen, refSeq, costMatrixRowCnt, costMatrix,
			    gapo, gape, bandW, endBonus, zdrop, initScore,
			    &qle, &tle, &gtle, &gscore, &max_off) ;
    
  else if (bandW < 16) 
    pscore = ksw_extend_16(queryLen, querySeq, refLen, refSeq, costMatrixRowCnt, costMatrix,
			     gapo, gape, bandW, endBonus, zdrop, initScore,
			     &qle, &tle, &gtle, &gscore, &max_off) ;
  
  else if (bandW < 32) {
    pscore = ksw_extend_32(queryLen, querySeq, refLen, refSeq, costMatrixRowCnt, costMatrix,
				gapo, gape, bandW, endBonus, zdrop, initScore,
				&qle, &tle, &gtle, &gscore, &max_off) ;

  }else { // if no vectorized solution exists, run the serial function 
    pscore = ksw_extend(queryLen, querySeq, refLen, refSeq, costMatrixRowCnt, costMatrix,
			gapo, gape, bandW, endBonus, zdrop, initScore,
			&qle, &tle, &gtle, &gscore, &max_off) ;
  }


  if (gscore <= 0 || gscore <= pscore - endBonus) {
    alignedQLen = qle ;
    alignedRLen = tle ;
    return pscore ;
  }
  else {
    alignedQLen = queryLen ;
    alignedRLen = gtle ;
    return gscore ;
  }
}

// Function used to pack query bits into uint64_t
inline void orBitAtPosition(uint64_t& w, bool newBit, int index) {
  w |= (((uint64_t)newBit) << index) ;
}

template<class BitVec>
inline BitVec computeMatchVec(const BitVec& Q0, const BitVec& Q1, const BitVec& Q2,
			      const BitVec& R0, const BitVec& R1, const BitVec& R2) {
  return ~Q2 & ~R2 & ((Q1 ^ R1) | (Q0 ^ R0)) ;
}

// This can be made faster if needed
int computeHammingDistance(const uint8_t* refSeq, int refLen, const uint8_t* querySeq, int queryLen, int upperBound, int& lastConsecutiveMatching) {
  
  int cnt = min(refLen, queryLen) ;
  int dist = 0 ;
  lastConsecutiveMatching = 0 ;
  for (int ci=0; ci < cnt; ++ci) {

    ++lastConsecutiveMatching ;
    if (refSeq[ci] != querySeq[ci]) {
      ++dist ;
      lastConsecutiveMatching = 0 ;
      if (dist >= upperBound)
	return dist ;
    }
  }

  return dist ;
}

inline bool isHammingDistanceSmallEnough(const uint8_t* refSeq, int refLen, const uint8_t* querySeq, int queryLen,
					 int& alignedQLenLow, int& alignedQLenHigh,
					 int& alignedRLenLow, int& alignedRLenHigh,
					 int& maxScoreLow, int& maxScoreHigh,
					 int mismatchWt, int gapOpenWt) {

	int minSeqLen = std::min(refLen, queryLen) ;

  // This heuristic is applicable to only queries that are long enough (e.g. >= 20 bases)
  if (minSeqLen < 20)
    return false ;

  // If there is a mismatch towards the end of the query, it is possible that partial
  // alignment is better than full alignment. So, the last consecutive match is returned in the
  // function below.
  int lastConsecutiveMatching=0 ;
  
  // If Hamming dist is greater than upper bound, no need to process the whole query
  // The value below can be tuned.
  int upperBound = minSeqLen / gapOpenWt ; 
  
  int hdist = computeHammingDistance(refSeq, refLen, querySeq, queryLen, upperBound, lastConsecutiveMatching) ;
    
  // If lastConsecutiveMatching is <= gapOpenWt, we don't have high confidence in
  // full alignment. 
  if (hdist < upperBound && lastConsecutiveMatching > gapOpenWt) {
#ifdef DEBUG0
    cerr << "The Hamming distance is " << hdist  << " with " << lastConsecutiveMatching
	 << " bases matches at the end of query. " << endl ;
    cerr << "Returning full alignment directly." << endl ;
#endif
    alignedQLenLow = alignedQLenHigh = minSeqLen ;
    alignedRLenLow = alignedRLenHigh = minSeqLen ;
    maxScoreLow = maxScoreHigh = (minSeqLen - hdist) - hdist * mismatchWt ;
    
    return true ;
  }

  return false ;
}
				  

template<class BitVec, class EDVec>
void fast_extend(const uint8_t* refSeq, int refLen, const uint8_t* querySeq, int queryLen,
		 int initScore, int endBonus, 
		 int& alignedQLenLow, int& alignedQLenHigh,
		 int& alignedRLenLow, int& alignedRLenHigh,
		 int& maxScoreLow, int& maxScoreHigh, SWFeedback& swfb,
		 int mismatchWt=4, int gapWt=1, int gapOpenWt=6, int ambigWt = 1) {


  assert (BitVec::isQueryLengthOk(queryLen)) ;

#ifdef DEBUG0
  cerr << "Query and ref lens for SW problem: " << queryLen << " " << refLen << endl ;
#endif

  // If full alignment is obvious, return immediately
  if (isHammingDistanceSmallEnough(refSeq, refLen, querySeq, queryLen,
																	 alignedQLenLow, alignedQLenHigh, 
																	 alignedRLenLow, alignedRLenHigh,
																	 maxScoreLow, maxScoreHigh,
																	 mismatchWt, gapOpenWt)) {
		maxScoreLow += initScore;
		maxScoreHigh += initScore;
    return ;
  }

  // Define the bit vectors that will be used in fast edit distance computation
  BitVec L0, L1, U0, U1, UP0, UP1 ;
  BitVec Q0, Q1, Q2 ;
  BitVec D[5] ; // profile vecs

  bool initL0 = 1 ; 
  bool initL1 = 0 ;

  bool initU0 = 1 ;
  bool initU1 = 0 ;

#ifdef FINE_ALIGN
  const int FINE_ALIGN_CNT = 16 ; // align the first 16 bases in a fine-grain way
  assert (FINE_ALIGN_CNT <= 16) ; // the rest of the code assumes this.
  int ambigBaseCntsFineArr[FINE_ALIGN_CNT] ;
  for (int i=0; i < FINE_ALIGN_CNT; ++i)
    ambigBaseCntsFineArr[i] = 0 ;
#endif

  // Keep track of the ambig bases corresponding to each score probe
  EDVec ambigBaseCnts(0) ;
  int ambigBaseCntsArr[EDVec::WORD_CNT()] ;
  for (int i=0; i < EDVec::WORD_CNT(); ++i)
    ambigBaseCntsArr[i] = 0 ;
  int ambigWordIndex = 0 ;

  // Store the query in bit vectors
  for (int wIndex=0; wIndex * 64 < queryLen; ++wIndex) {
    uint64_t bit0Partial = 0x0 ;
    uint64_t bit1Partial = 0x0 ;
    uint64_t bit2Partial = 0x0 ;

    for (int bIndex=0; bIndex < 64; ++bIndex) {
      int qi = wIndex * 64 + bIndex ;
      if (qi >= queryLen)
	break ;
      uint8_t qbase = querySeq[qi] ;
      if (qbase == 4) {

	while (qi >= EDVec::getDistAt(ambigWordIndex, queryLen))
	  ++ambigWordIndex ;
	++ambigBaseCntsArr[ambigWordIndex] ;

#ifdef FINE_ALIGN
	if (qi < FINE_ALIGN_CNT)
	  ++ambigBaseCntsFineArr[qi] ;
#endif
      }

      orBitAtPosition(bit0Partial, qbase & 0x1, bIndex) ;
      orBitAtPosition(bit1Partial, (qbase & 0x2) >> 1, bIndex) ;
      orBitAtPosition(bit2Partial, (qbase & 0x4) >> 2, bIndex) ;
    }
    Q0.setWord(wIndex, bit0Partial) ;
    Q1.setWord(wIndex, bit1Partial) ;
    Q2.setWord(wIndex, bit2Partial) ;
  }
  
  
  // Compute running sum of ambig base cnts
  for (int wi=1; wi < EDVec::WORD_CNT(); ++wi) {
    ambigBaseCntsArr[wi] += ambigBaseCntsArr[wi-1] ;
  }
  ambigBaseCnts.setWords(ambigBaseCntsArr) ;
  

#ifdef FINE_ALIGN
  for (int wi=1; wi < FINE_ALIGN_CNT; ++wi) {
    ambigBaseCntsFineArr[wi] += ambigBaseCntsFineArr[wi-1] ;
  }
#endif

  // pre-compute the profile vectors
  BitVec B0, B1 ;
  B0.setAllZeroes() ;
  B1.setAllOnes() ;
  
  D[0] = computeMatchVec(Q0, Q1, Q2, B0, B0, B0) ;
  D[1] = computeMatchVec(Q0, Q1, Q2, B1, B0, B0) ;
  D[2] = computeMatchVec(Q0, Q1, Q2, B0, B1, B0) ;
  D[3] = computeMatchVec(Q0, Q1, Q2, B1, B1, B0) ;
  D[4] = computeMatchVec(Q0, Q1, Q2, B0, B0, B1) ;
  

  // Initialize boundary conditions
  BitVec initL0Vec, initL1Vec ;
  initL0Vec.setLSBClearRest(initL0) ;
  initL1Vec.setLSBClearRest(initL1) ;

  U0.setAllBits(initU0) ;
  U1.setAllBits(initU1) ;

  // Init the edit distance vec.
  DistVec<BitVec,EDVec> accumDist(queryLen) ;
  accumDist.initDist(initU0, initU1) ;

#ifdef SW_FEEDBACK
  EDVec bestAccumDist = accumDist.getVec() ;
#endif

  int INF_SCORE = EDVec::getMaxWordVal() ;

  EDVec currColVec(0) ;
  EDVec bestScore(0), bestColVec(0) ;

  EDVec queryLenVec ;
  queryLenVec.setWordsAsDist(queryLen, 1) ;

#ifdef FINE_ALIGN
  FineAlignment16 fineAlign(queryLen, endBonus, mismatchWt, ambigWt) ;
#endif

  EDVec mismatchWtP1Vec(mismatchWt+1), // P1: plus 1
    gapWtVec(gapWt), gapWtP1Vec(gapWt+1), gapOpenWtVec(gapOpenWt), ambigWtP1Vec(ambigWt+1), 
    zeroVec(0) ;
  
  EDVec endBonusVec(0), deltaScoreVec(initScore), badScoreVec(0) ;
  endBonusVec.setWordsAsEndBonus(queryLen, endBonus) ;
  deltaScoreVec = deltaScoreVec + queryLenVec - ambigBaseCnts * ambigWtP1Vec + endBonusVec ;
  badScoreVec.setWordsAsBadScore(queryLen, 0, INF_SCORE) ;

  // Every time the edit distance is increased by 1, we add -(mismatchWt+1) to score.
  // For an INSERT (downwards move) in the lower-triangle part of the matrix:
  //       * Increasing column index means reducing the insert count. This should increase score
  //         by (gapWt+1), i.e. replacing one insertion by one match.
  //       * However, if this change is not reflected in the edit distance, this means that
  //         the insertion was replaced by a mismatch.
  //       * So, the delta cost for INSERT is: (gapWt+1) - (mismatchWt+1) = gapWt - mismatchWt
  // 
  // For a DELETE (right move) in the upper-triangle part of the matrix:
  //       * Increasing column index means increasing the delete count. This should decrease score
  //         by gapWt.
  //       * However, we also pay for mismatch cost of (mismatchWt+1), assuming that edit dist has
  //         increased by 1. We should not pay for this mismatch cost.
  //       * So, the delta cost for DELETE is: -gapWt + (mismatchWt+1)

  EDVec deltaCostInsert((gapWt-mismatchWt)) ; // insert costs decrease as ci increases

  EDVec deltaCostDelete(0) ; // delete costs increase as ci increases
  EDVec deltaCostDistWt(-(mismatchWt+1)) ;
  EDVec deltaColVec(1) ;
  
  EDVec insertVecModifier(0), deleteVecModifier(0), gapOpenDelta(0) ;
  insertVecModifier.setFirstWord(-(gapWt-mismatchWt)) ; // will cancel one entry when added
  deleteVecModifier.setFirstWord(-gapWt + mismatchWt + 1) ; // will add deleteWt to one entry
  gapOpenDelta.setFirstWord(gapOpenWt) ; // will remove/add gapOpenWt from/to one diag entry

  EDVec firstColCost = deltaScoreVec - queryLenVec * gapWtP1Vec - gapOpenWtVec ;

  EDVec currScore = firstColCost ;
  EDVec deltaMismatch(0) ;

#ifdef DEBUG1
  cerr << " The init score is: " << initScore << endl ;
  cerr << " The cost of the 0th column: " << firstColCost << endl ;
#endif

  // Score is computed as follows:
  // (Note that ambigCnt and gapCnt are NOT included in editDist.)
  // (If the gap type is deletion and there are ambig bases in query, the score estimated
  //  is lower than the actual score. Some of the ambig bases may be skipped in deletions, but
  //  we don't capture this. This can be fixed later by subtracting delete-cnt from ambigCnt
  //  at every column less than qlen.)
  // score =  (qlen - editDist - ambigCnt - insertCnt) // score corresponding to matches
  //          - ambigCnt * ambigWt // score corresponding to ambig bases
  //          - (deleteCnt+insertCnt) * gapWt // score corresponding to gap extends
  //          - editDist * mismatchWt // score corresponding to mismatches
  // Reordering the constant terms, we have:
  // score =  - insertCnt * (gapWt+1)     // variable term
  //          - deleteCnt * gapWt         // variable term
  //          - editDist * (mismatchWt+1) // variable term
  //          + qlen - ambigCnt * (ambigWt+1) // constant terms


#define PRINT_UPDATE_COST_DEBUG(__ci) {					\
    cerr << "In iteration " << __ci << ", the cost values are updated to:" << endl ; \
    cerr << "L0: " << std::hex << L0 << std::dec << endl ;		\
    cerr << "L1: " << std::hex << L1 << std::dec << endl ;		\
    cerr << "currColVec: " << currColVec << endl ;			\
    cerr << "queryLenVec: " << queryLenVec << endl ;			\
    cerr << "deltaCostInsert: " << deltaCostInsert << endl ;	        \
    cerr << "deltaCostDelete: " << deltaCostDelete << endl ;      	\
    cerr << "deltaMismatch: " << deltaMismatch << endl ;		\
    cerr << "deltaCostDistWt: " << deltaCostDistWt << endl ;		\
    cerr << "currScore: " << currScore << endl ;			\
    cerr << "bestCol: " << bestColVec << endl ;				\
    cerr << "bestScore: " << bestScore << endl ;			\
  }

#define PRINT_DELTA_MODIFICATIONS_DEBUG(__ci) {			       \
    cerr << "The delta costs have been updated after processing column " << __ci << endl ; \
    cerr << "New delta cost insert: " << deltaCostInsert << endl ;	\
    cerr << "New delta cost delete: " << deltaCostDelete << endl ;	\
    cerr << "using the modifiers: " << endl ;				\
    cerr << "Insert vec modifier: " << insertVecModifier << endl ;	\
    cerr << "Delete vec modifier: " << deleteVecModifier << endl ;	\
  }


#define COMPUTE_LU(__ci) {			\
    uint8_t refChar = refSeq[ci] ;		\
    const BitVec& currD = D[refChar] ;		\
    BitVec SL1 = U0.andnot(currD) ;		\
    SL1.shiftLeftAndInsert(initL1) ;		\
						\
    BitVec INJ = U0 & SL1 ;			\
    BitVec SUM = INJ + U0 ;			\
    BitVec CIN = SUM ^ INJ ^ U0 ;		\
						\
    L1 = SL1 | CIN ;				\
    L0 = U1 | currD.andnot(L1).andnot(U0) ;	\
    L0.shiftLeftAndInsert(initL0) ;		\
						\
    BitVec TU = currD.andnot(U1) ;		\
    U1 = L0.andnot(TU) ;			\
    U0 = L1 | TU.andnot(L0) ;			\
						\
    /* Add delta costs corresponding to the current mismatch delta (-1, 0, or 1) */ \
    deltaMismatch = accumDist.getDeltaDistVec(L0, L1) ;			\
    accumDist.addDist(deltaMismatch) ;					\
    									\
  }

#define UPDATE_COSTS(__ci) {						\
    currColVec += deltaColVec ;						\
    EDVec gapDelta = deltaCostInsert + deltaCostDelete ;		\
    /*gapDelta.addThirdIfFirstGTSecond(gapDelta, zeroVec, gapOpenWtVec, zeroVec) ; */ \
    currScore += gapDelta ;						\
									\
    currScore += (deltaMismatch * deltaCostDistWt) ;			\
    bestScore.setMax(currScore, bestColVec, currColVec) ;		\
  }

#define UPDATE_COSTS_AND_CHECK(__ci) {					\
    currColVec += deltaColVec ;						\
    EDVec gapDelta = deltaCostInsert + deltaCostDelete ;		\
    /*gapDelta.addThirdIfFirstGTSecond(gapDelta, zeroVec, gapOpenWtVec, zeroVec) ; */ \
    currScore += gapDelta ;						\
									\
    /* Add delta costs corresponding to the current mismatch delta (-1, 0, or 1) */  \
    deltaMismatch = accumDist.getDeltaDistVec(L0, L1) ;			\
									\
    currScore += (deltaMismatch * deltaCostDistWt) ;			\
    updateFlag |= bestScore.setMaxAndReturnFlag(currScore, bestColVec, currColVec) ; \
  }

#define UPDATE_COSTS_AND_STORE_BEST_ACCUM_DIST(__ci) {			\
    currColVec += deltaColVec ;						\
    EDVec gapDelta = deltaCostInsert + deltaCostDelete ;		\
    /*gapDelta.addThirdIfFirstGTSecond(gapDelta, zeroVec, gapOpenWtVec, zeroVec) ; */ \
    currScore += gapDelta ;						\
									\
    /* Add delta costs corresponding to the current mismatch delta (-1, 0, or 1) */  \
    deltaMismatch = accumDist.getDeltaDistVec(L0, L1) ;			\
									\
    currScore += (deltaMismatch * deltaCostDistWt) ;			\
    bestScore.setMax(currScore, bestColVec, currColVec, bestAccumDist, accumDist.getVec()) ; \
  }

#define UPDATE_COSTS_AND_CHECK_AND_STORE_BEST_ACCUM_DIST(__ci) {	\
    currColVec += deltaColVec ;						\
    EDVec gapDelta = deltaCostInsert + deltaCostDelete ;		\
    /*gapDelta.addThirdIfFirstGTSecond(gapDelta, zeroVec, gapOpenWtVec, zeroVec) ; */ \
    currScore += gapDelta ;						\
									\
    /* Add delta costs corresponding to the current mismatch delta (-1, 0, or 1) */  \
    deltaMismatch = accumDist.getDeltaDistVec(L0, L1) ;		\
									\
    currScore += (deltaMismatch * deltaCostDistWt) ;			\
    updateFlag |= bestScore.setMaxAndReturnFlag(currScore, bestColVec, currColVec, bestAccumDist, accumDist.getVec()) ; \
  }


#define REMOVE_GAP_OPEN_PENALTY_FOR_DIAG(__ci) {			\
    currScore += gapOpenDelta ;						\
  }

#define ADD_GAP_OPEN_PENALTY_TO_DIAG(__ci) {				\
    currScore -= gapOpenDelta ;						\
    									\
    /* shift the modifier so that they modify the next lane in the next update */ \
    gapOpenDelta.shiftWordsLeftByOne() ;				\
  }

#define PERFORM_DELTA_MODIFICATIONS(__ci) {				\
    deltaCostInsert += insertVecModifier ;				\
    deltaCostDelete += deleteVecModifier ;				\
    									\
    /* shift the modifiers so that they modify the next lane in the next update */ \
    insertVecModifier.shiftWordsLeftByOne() ;				\
    deleteVecModifier.shiftWordsLeftByOne() ;				\
    									\
  }

#define UPDATE_FINE_ALIGN(__ci) {		\
    if (__ci < FINE_ALIGN_CNT)						\
      fineAlign.update(U0.getLow16Bits(), U1.getLow16Bits(), __ci+1,	\
		       ambigBaseCntsFineArr[__ci]) ;			\
    									\
  }

#define CHECK_TERMINATION(__ci) {				\
    if (currScore.allLessThanOrEqualTo(badScoreVec)) {		\
      goto loopEnd ;						\
    }								\
  }

#define CHECK_COSTS(__ci) {						\
    EDVec deleteCntVec = currColVec.subSat(queryLenVec) ;		\
    EDVec insertCntVec = queryLenVec.subSat(currColVec) ;		\
    EDVec gapCntVec = insertCntVec + deleteCntVec ;			\
									\
    EDVec gapPenaltyVec = gapCntVec * gapWtVec + insertCntVec ;		\
    gapPenaltyVec.addThirdIfFirstGTSecond(gapCntVec, zeroVec, gapOpenWtVec, zeroVec) ; \
									\
    EDVec currScoreVecTest =						\
      deltaScoreVec - gapPenaltyVec - accumDist.getVec() * mismatchWtP1Vec ; \
    EDVec currScoreVecTest =						\
      deltaScoreVec - gapPenaltyVec - accumDist.getVec() * mismatchWtP1Vec ;	\
									\
    if (!(currScoreVecTest == currScore) || !(currScoreVecTest == currScore)) { \
      cerr << "Error: Mismatch in cost values computed using two methods at column " << ci << endl ; \
      cerr << "deleteCntVec: " << deleteCntVec << endl ;		\
      cerr << "insertCntVec: " << insertCntVec << endl ;		\
      cerr << "gapPenaltyVec: " << gapPenaltyVec << endl ;		\
      cerr << "score (test): " << currScoreVecTest << endl ;		\
      exit(0) ;								\
    }									\
  }

  // If there are 3 words in EDVec, flagMask will be 0....0111111 in binary.
  // The flagMask has 2 bits per word, because of SSE limitations.
  uint32_t flagMask = (uint32_t) ((((uint64_t)1) << (2*(EDVec::getLastWordIndexFor(queryLen)+1))) - 1) ;

  int firstProbeOffset = min(refLen-1, accumDist.getProbeColumn(0)) ;
  int numFullIters ;
  
  int ci = 0 ;
  for (; ci < firstProbeOffset; ++ci) {
    COMPUTE_LU(ci) ;

#ifdef SW_FEEDBACK
    UPDATE_COSTS_AND_STORE_BEST_ACCUM_DIST(ci) ;
#else
    UPDATE_COSTS(ci) ;
#endif

#ifdef FINE_ALIGN
    UPDATE_FINE_ALIGN(ci) ;
#endif
#ifdef DEBUG1
    PRINT_UPDATE_COST_DEBUG(ci) ;
#endif
    //CHECK_COSTS(ci) ;
  }
  
  REMOVE_GAP_OPEN_PENALTY_FOR_DIAG(ci) ;
  
  if (firstProbeOffset >= 0) {
    COMPUTE_LU(ci) ;

#ifdef SW_FEEDBACK
    UPDATE_COSTS_AND_STORE_BEST_ACCUM_DIST(ci) ;
#else
    UPDATE_COSTS(ci) ;
#endif

#ifdef FINE_ALIGN
    UPDATE_FINE_ALIGN(ci) ;
#endif
#ifdef DEBUG1
    PRINT_UPDATE_COST_DEBUG(ci) ;
#endif
    //CHECK_COSTS(ci) ;
    ++ci ; 
  }
  ADD_GAP_OPEN_PENALTY_TO_DIAG(ci-1) ;

  PERFORM_DELTA_MODIFICATIONS(ci-1) ;
#ifdef DEBUG1
  PRINT_DELTA_MODIFICATIONS_DEBUG(ci-1) ;
#endif
  CHECK_TERMINATION(ci-1) ;


	numFullIters = (min(queryLen, refLen)-ci) /EDVec::PERIOD() ;
  for (int iterIndex = 0; iterIndex < numFullIters; ++iterIndex) {

    for (int subIterIndex=0; subIterIndex < EDVec::PERIOD()-1; ++subIterIndex, ++ci) {
      COMPUTE_LU(ci) ;
#ifdef SW_FEEDBACK
    UPDATE_COSTS_AND_STORE_BEST_ACCUM_DIST(ci) ;
#else
    UPDATE_COSTS(ci) ;
#endif

#ifdef FINE_ALIGN
      UPDATE_FINE_ALIGN(ci) ;
#endif
#ifdef DEBUG1
      PRINT_UPDATE_COST_DEBUG(ci) ;
#endif
      //CHECK_COSTS(ci) ;
    }
    REMOVE_GAP_OPEN_PENALTY_FOR_DIAG(ci) ;
    COMPUTE_LU(ci) ;

#ifdef SW_FEEDBACK
    UPDATE_COSTS_AND_STORE_BEST_ACCUM_DIST(ci) ;
#else
    UPDATE_COSTS(ci) ;
#endif

#ifdef FINE_ALIGN
    UPDATE_FINE_ALIGN(ci) ;
#endif
#ifdef DEBUG1
    PRINT_UPDATE_COST_DEBUG(ci) ;
#endif
    //CHECK_COSTS(ci) ;
    ADD_GAP_OPEN_PENALTY_TO_DIAG(ci) ;

    PERFORM_DELTA_MODIFICATIONS(ci) ;
#ifdef DEBUG1
    PRINT_DELTA_MODIFICATIONS_DEBUG(ci) ;
#endif

    CHECK_TERMINATION(ci) ;
    ++ci ;
  }


  // If there are 3 words in EDVec, flagMask will be 0....0111111 in binary
  // flagMask has 2 bits per word, because of SSE limitations

  // Next, process the columns with indices greater than or equal to query length
  numFullIters = (refLen-ci) / EDVec::PERIOD() ;
  for (int iterIndex = 0; iterIndex < numFullIters; ++iterIndex) {
    uint32_t updateFlag = 0x0 ; // 2 bits per word, because of SSE limitations

    for (int subIterIndex=0; subIterIndex < EDVec::PERIOD(); ++subIterIndex, ++ci) {
      COMPUTE_LU(ci) ;

#ifdef SW_FEEDBACK
      UPDATE_COSTS_AND_CHECK_AND_STORE_BEST_ACCUM_DIST(ci) ;
#else
      UPDATE_COSTS_AND_CHECK(ci) ;
#endif


#ifdef DEBUG1
      PRINT_UPDATE_COST_DEBUG(ci) ;
#endif
      //CHECK_COSTS(ci) ; // will fail for the padded entries
    }

    updateFlag &= flagMask ;

    if (!updateFlag) // if the best score hasn't been updated in any iteration
      goto loopEnd ;

    CHECK_TERMINATION(ci) ;
  }


  // Next, process the remaining columns
  for (; ci < refLen; ++ci) {
    COMPUTE_LU(ci) ;

#ifdef SW_FEEDBACK
    UPDATE_COSTS_AND_STORE_BEST_ACCUM_DIST(ci) ;
#else
    UPDATE_COSTS(ci) ;
#endif

#ifdef DEBUG1
    PRINT_UPDATE_COST_DEBUG(ci) ;
#endif
    // CHECK_COSTS(ci) ; // will fail for the padded entries
  }

 loopEnd:
#ifdef DEBUG0
  if (ci < refLen)
    cerr << "Loop terminated after iteration " << ci << endl ;
#endif

  PartialAlignment partials[EDVec::WORD_CNT()+2] ;

  // no alignment
  partials[0] = PartialAlignment(0, 0, initScore, initScore) ;
  maxScoreHigh = initScore ;
  maxScoreLow = initScore ;

#ifdef SW_FEEDBACK
  partials[0].estimatedBand = 0 ;
#endif

  int partialCnt = 1 ;

#ifdef FINE_ALIGN
  // result of fine_align
  int alignedLenFine, scoreFine ;
  fineAlign.getBest(alignedLenFine, scoreFine) ;
  scoreFine += initScore ;
  partials[1] = PartialAlignment(alignedLenFine, alignedLenFine, scoreFine, scoreFine) ;
#ifdef SW_FEEDBACK
  partials[1].estimatedBand = 1 ;
#endif

  if (scoreFine > maxScoreHigh) {
    maxScoreHigh = scoreFine ;
  }
  if (scoreFine > maxScoreLow) {
    maxScoreLow = scoreFine ;
  }

  ++partialCnt ;
#endif // #ifdef FINE_ALIGN

  for (int di=0; di < accumDist.getValidDistCnt(); ++di) {
    int currQLen = queryLenVec.getWord(di) ;
    int currRLen = bestColVec.getWord(di) ;
    int currScoreHigh = bestScore.getWord(di) ;
    
    int gapCnt = abs(currQLen - currRLen) ;

    int currScoreLow = (gapCnt == 0) ? currScoreHigh :
      (currScoreHigh - gapCnt * (mismatchWt - gapWt) + gapOpenWt) ;

    if (currScoreHigh > maxScoreHigh) {
      maxScoreHigh = currScoreHigh ;
    }
    if (currScoreLow > maxScoreLow) {
      maxScoreLow = currScoreLow ;
    }

    partials[partialCnt] = PartialAlignment(currQLen, currRLen, currScoreLow, currScoreHigh) ;

#ifdef SW_FEEDBACK
    int currEditDist = bestAccumDist.getWord(di) ;
    partials[partialCnt].estimatedBand = gapCnt + max(2, (currEditDist - gapCnt) / 4) ;
#endif

    ++partialCnt ;
  }
  
  int scoreToleranceForMin = EDVec::PERIOD()-1 ;
  int scoreToleranceForMax = EDVec::PERIOD()/2 ;

#ifdef SW_FEEDBACK
    int lastValidBand = 0 ;
    int lastValidScore = initScore ;
#endif

  alignedQLenLow = queryLen+1 ;
  alignedQLenHigh = 0 ;
  alignedRLenLow = 0 ;
  alignedRLenHigh = 0 ;
  for (int ai=0; ai < partialCnt; ++ai) {

#ifdef DEBUG0
    cerr << partials[ai] << endl ;
#endif
  
    int currToleranceForMin = scoreToleranceForMin ;
    int currToleranceForMax = scoreToleranceForMax ;

#ifdef FINE_ALIGN
    if (ai != 1 && partials[ai].partialQLen < FINE_ALIGN_CNT && partials[ai].scoreHigh < initScore)
      currToleranceForMax = 0 ;
#endif // #ifdef FINE_ALIGN
    
#ifdef SW_FEEDBACK
    if (partials[ai].estimatedBand >= lastValidBand + 8 && 
	partials[ai].scoreHigh <= lastValidScore) {

      continue ;
    }
    lastValidBand = partials[ai].estimatedBand ;
    lastValidScore = partials[ai].scoreHigh ;
#endif

    if ((partials[ai].scoreLow >= maxScoreLow - currToleranceForMin) ||
       	(partials[ai].scoreHigh >= maxScoreHigh - currToleranceForMin)) {

      if (partials[ai].partialQLen < alignedQLenLow) {

	alignedQLenLow = partials[ai].partialQLen ;
	alignedRLenLow = partials[ai].partialRLen ;
	
#ifdef SW_FEEDBACK
	swfb.updateMin(alignedQLenLow, alignedRLenLow, 
		       (partials[ai].scoreLow + partials[ai].scoreHigh)/2) ;
#endif
      }
    }

    if ((partials[ai].scoreLow >= maxScoreLow - currToleranceForMax) ||
	(partials[ai].scoreHigh >= maxScoreHigh - currToleranceForMax)) {

      if (partials[ai].partialQLen > alignedQLenHigh) {
	alignedQLenHigh = partials[ai].partialQLen ;
	alignedRLenHigh = partials[ai].partialRLen ;
#ifdef SW_FEEDBACK
	swfb.updateMax(alignedQLenHigh, alignedRLenHigh, partials[ai].estimatedBand) ;
#endif      
      }
    }
  }

  
#ifdef SW_FEEDBACK
  // Do adjustments on the feedback values
  if (swfb.maxQLen > 0)
    swfb.maxQLen = min(queryLen, swfb.maxQLen + EDVec::PERIOD()-1) ;

  swfb.maxRLen = max(swfb.maxRLen, min(refLen, swfb.maxQLen + swfb.maxBand)) ;

  if (swfb.maxQLen < swfb.minQLen)
    std::swap(swfb.maxQLen, swfb.minQLen) ;
  
  if (swfb.maxRLen < swfb.minRLen)
    std::swap(swfb.maxRLen, swfb.minRLen) ;


  if (abs(swfb.minQLen-swfb.minRLen) > 8) { 
    // We don't have high confidence for the long gaps. To be conservative, start ksw_extend from 
    // the beginning, i.e. set the min query and ref indices to 0:
    swfb.updateMin(0, 0, initScore) ;
  }

#endif

  assert(alignedQLenLow <= queryLen) ;
  assert(alignedQLenHigh >= 0) ;

#ifdef DEBUG0
  cerr << "Max scores computed: [" << maxScoreLow << ", " << maxScoreHigh << "]" << endl ;
  cerr << "Query alignment range: [" << alignedQLenLow << ", " << alignedQLenHigh << "]" << endl ;
  cerr << "Ref alignment range: [" << alignedRLenLow << ", " << alignedRLenHigh << "]" << endl ;
#ifdef SW_FEEDBACK
  cerr << "The feedback computed: " << swfb << endl ;
#endif // #ifdef SW_FEEDBACK
#endif // #ifdef DEBUG0


  // Adjust the alignedRLen value for full alignment for obvious cases
  if (alignedQLenHigh == alignedQLenLow && alignedQLenHigh == queryLen && 
      alignedRLenHigh == alignedRLenLow && alignedQLenHigh >= 4 && alignedRLenHigh >= 5) {

    int origMatchCnt = 
      (refSeq[alignedRLenHigh-1] == querySeq[alignedQLenHigh-1])
      + (refSeq[alignedRLenHigh-2] == querySeq[alignedQLenHigh-2])
      + (refSeq[alignedRLenHigh-3] == querySeq[alignedQLenHigh-3])
      + (refSeq[alignedRLenHigh-4] == querySeq[alignedQLenHigh-4]) ;

    int back1MatchCnt =
      (refSeq[alignedRLenHigh-2] == querySeq[alignedQLenHigh-1])
      + (refSeq[alignedRLenHigh-3] == querySeq[alignedQLenHigh-2])
      + (refSeq[alignedRLenHigh-4] == querySeq[alignedQLenHigh-3]) 
      + (refSeq[alignedRLenHigh-5] == querySeq[alignedQLenHigh-4]) ;
    
    int forw1MatchCnt = 0 ;
    if (alignedRLenHigh < refLen-1) {
      forw1MatchCnt = 
	(refSeq[alignedRLenHigh] == querySeq[alignedQLenHigh-1])
	+ (refSeq[alignedRLenHigh-1] == querySeq[alignedQLenHigh-2])
	+ (refSeq[alignedRLenHigh-2] == querySeq[alignedQLenHigh-3])
	+ (refSeq[alignedRLenHigh-3] == querySeq[alignedQLenHigh-4]) ;
    }

    if (back1MatchCnt > origMatchCnt && back1MatchCnt >= forw1MatchCnt) {
      --alignedRLenHigh ;
      --alignedRLenLow ;
    }
    else if (forw1MatchCnt > origMatchCnt) {
      ++alignedRLenHigh ;
      ++alignedRLenLow ;
    }
    
#ifdef DEBUG0
    cerr << "Orig/Back-by-1/Forward-by-1 match counts: " << origMatchCnt << " "
	 << back1MatchCnt << " " << forw1MatchCnt << endl ;
    cerr << "The aligned ref len is now: " << alignedRLenHigh << endl ;
#endif

  }

  // In BWA, endBonus is considered during comparisons, but it is not included in the final
  // score value returned. To be consistent with BWA, we do the following:
  if (alignedQLenHigh == queryLen)
    maxScoreHigh -= endBonus ;
  if (alignedQLenLow == queryLen)
    maxScoreLow -= endBonus ;

#ifdef DEBUG0
  cerr << "Init score = " << initScore << endl ;
  cerr << "Max score-high = " << maxScoreHigh << " for alignedQLen = " 
       << alignedQLenHigh << ", alignedRLen = " << alignedRLenHigh << endl ;
  cerr << "Max score-low = " << maxScoreLow << " for alignedQLen = " 
       << alignedQLenLow << ", alignedRLen = " << alignedRLenLow << endl ;
#endif

}




void fast_extend_u32(const uint8_t* refSeq, int refLen, const uint8_t* querySeq, int queryLen, 
		     int initScore, int endBonus, 
		     int& alignedQLenLow, int& alignedQLenHigh, 
		     int& alignedRLenLow, int& alignedRLenHigh, 
		     int& scoreLow, int& scoreHigh, SWFeedback& swFeedback,
		     int mismatchWt=4, int gapWt=1, int gapOpenWt=6, int ambigWt = 1) {
  
  fast_extend<BitVec64, EDVec32Every4>(refSeq, refLen, querySeq, queryLen, 
				       initScore, endBonus, 
				       alignedQLenLow, alignedQLenHigh, 
				       alignedRLenLow, alignedRLenHigh, 
				       scoreLow, scoreHigh, swFeedback,
				       mismatchWt, gapWt, gapOpenWt, ambigWt) ;
}

void fast_extend_u64(const uint8_t* refSeq, int refLen, const uint8_t* querySeq, int queryLen, 
		     int initScore, int endBonus, 
		     int& alignedQLenLow, int& alignedQLenHigh, 
		     int& alignedRLenLow, int& alignedRLenHigh, 
		     int& scoreLow, int& scoreHigh, SWFeedback& swFeedback,
		     int mismatchWt=4, int gapWt=1, int gapOpenWt=6, int ambigWt = 1) {
  
  fast_extend<BitVec64, EDVec64Every8>(refSeq, refLen, querySeq, queryLen, 
				       initScore, endBonus, 
				       alignedQLenLow, alignedQLenHigh, 
				       alignedRLenLow, alignedRLenHigh, 
				       scoreLow, scoreHigh, swFeedback,
				       mismatchWt, gapWt, gapOpenWt, ambigWt) ;
}

void fast_extend_u64x2(const uint8_t* refSeq, int refLen, const uint8_t* querySeq, int queryLen, 
		       int initScore, int endBonus, 
		       int& alignedQLenLow, int& alignedQLenHigh, 
		       int& alignedRLenLow, int& alignedRLenHigh, 
		       int& scoreLow, int& scoreHigh, SWFeedback& swFeedback,
		       int mismatchWt=4, int gapWt=1, int gapOpenWt=6, int ambigWt = 1) {

  fast_extend<BitVec64x2, EDVec128Every16>(refSeq, refLen, querySeq, queryLen, 
					   initScore, endBonus, 
					   alignedQLenLow, alignedQLenHigh, 
					   alignedRLenLow, alignedRLenHigh, 
					   scoreLow, scoreHigh, swFeedback,
					   mismatchWt, gapWt, gapOpenWt, ambigWt) ;
}

void fast_extend_u128(const uint8_t* refSeq, int refLen, const uint8_t* querySeq, int queryLen, 
		      int initScore, int endBonus, 
		      int& alignedQLenLow, int& alignedQLenHigh, 
		      int& alignedRLenLow, int& alignedRLenHigh, 
		      int& scoreLow, int& scoreHigh, SWFeedback& swFeedback,
		      int mismatchWt=4, int gapWt=1, int gapOpenWt=6, int ambigWt = 1) {

  fast_extend<BitVec128, EDVec128Every16>(refSeq, refLen, querySeq, queryLen, 
					  initScore, endBonus, 
					  alignedQLenLow, alignedQLenHigh, 
					  alignedRLenLow, alignedRLenHigh, 
					  scoreLow, scoreHigh, swFeedback,
					  mismatchWt, gapWt, gapOpenWt, ambigWt) ;
}

void fast_extend_u256(const uint8_t* refSeq, int refLen, const uint8_t* querySeq, int queryLen, 
		      int initScore, int endBonus, 
		      int& alignedQLenLow, int& alignedQLenHigh, 
		      int& alignedRLenLow, int& alignedRLenHigh, 
		      int& scoreLow, int& scoreHigh, SWFeedback& swFeedback,
		      int mismatchWt=4, int gapWt=1, int gapOpenWt=6, int ambigWt = 1) {

#ifdef USE_AVX2
    if (useavx2) {
	fast_extend<BitVec256, EDVec256Every16> 
    (refSeq, refLen, querySeq, queryLen, initScore, endBonus, alignedQLenLow, 
     alignedQLenHigh, alignedRLenLow, alignedRLenHigh, scoreLow, scoreHigh, swFeedback,
     mismatchWt, gapWt, gapOpenWt, ambigWt) ;
    } else {
#endif
    fast_extend<BitVec128x2, EDVec128x2Every16>
    (refSeq, refLen, querySeq, queryLen, initScore, endBonus, alignedQLenLow, 
     alignedQLenHigh, alignedRLenLow, alignedRLenHigh, scoreLow, scoreHigh, swFeedback,
     mismatchWt, gapWt, gapOpenWt, ambigWt) ;
#ifdef USE_AVX2
	}
#endif

}


void fast_extend(const uint8_t* refSeq, int refLen, const uint8_t* querySeq, int queryLen,
		 int initScore, int endBonus, 
		 int& alignedQLenLow, int& alignedQLenHigh, 
		 int& alignedRLenLow, int& alignedRLenHigh, 
		 int& scoreLow, int& scoreHigh, SWFeedback& swFeedback,
		 int mismatchWt=4, int gapWt=1, int gapOpenWt=6, int ambigWt = 1) {
  
  if (queryLen < 15) {
    cerr << "Warning: It is not recommended to use SW filtering API for queries shorter than 15 bases!" << endl ;
  }

  if (queryLen <= 32)
    fast_extend_u32(refSeq, refLen, querySeq, queryLen, initScore, endBonus,
		    alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh,
		    scoreLow, scoreHigh, swFeedback, mismatchWt, gapWt, gapOpenWt, ambigWt) ;

  else if (queryLen <= 63)
    fast_extend_u64(refSeq, refLen, querySeq, queryLen, initScore, endBonus,
		    alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh,
		    scoreLow, scoreHigh, swFeedback, mismatchWt, gapWt, gapOpenWt, ambigWt) ;


  else if (queryLen <= 127)
    fast_extend_u128(refSeq, refLen, querySeq, queryLen, initScore, endBonus,
		     alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh,
		     scoreLow, scoreHigh, swFeedback, mismatchWt, gapWt, gapOpenWt, ambigWt) ;

  else if (queryLen <= 255)
    fast_extend_u256(refSeq, refLen, querySeq, queryLen, initScore, endBonus,
		     alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh,
		     scoreLow, scoreHigh, swFeedback, mismatchWt, gapWt, gapOpenWt, ambigWt) ;

  else {
    // The return values indicate an ambigious match
    alignedQLenLow = alignedRLenLow = 0 ;
    alignedQLenHigh = queryLen ;
    alignedRLenHigh = refLen ;
    scoreLow = 0 ;
    scoreHigh = queryLen ;
    cerr << "Warning: Skipping edit distance computation because reads longer than 126 bases not supported yet" << endl ;
  }
}


void fast_filter(const uint8_t* refSeq, int refLen, const uint8_t* querySeq, int queryLen,
		 int initScore, int endBonus, 
		 int& alignedQLen, int& alignedRLen, int& score, float& confidence,
		 int mismatchWt, int gapWt, int gapOpenWt, int ambigWt) {
  
  int alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh, scoreLow, scoreHigh ;
  SWFeedback swFeedback ;

  fast_extend(refSeq, refLen, querySeq, queryLen, initScore, endBonus,
	      alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh,
	      scoreLow, scoreHigh, swFeedback, 
	      mismatchWt, gapWt, gapOpenWt, ambigWt) ;
  
  alignedQLen = (alignedQLenLow + alignedQLenHigh)/2 ;
  alignedRLen = (alignedRLenLow + alignedRLenHigh)/2 ;
  score = (scoreLow + scoreHigh)/2 ;
  confidence = 
    (alignedQLenLow == alignedQLenHigh && alignedRLenLow == alignedRLenHigh) ? 1.0 : 0.0 ;

}

void fast_filter_and_extend(const uint8_t* refSeq, int refLen, const uint8_t* querySeq, int queryLen,
			    int initScore, int endBonus, int zdrop,
			    int& alignedQLen, int& alignedRLen, int& score,
			    int costMatrixRowCnt, const int8_t* costMatrix, 
			    int gapWt, int gapOpenWt) {

  int mismatchWt = -1 * costMatrix[1] ;
  int ambigWt = -1 * costMatrix[costMatrixRowCnt*costMatrixRowCnt-1] ;

  int alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh, scoreLow, scoreHigh ;
  SWFeedback swFeedback ;

  fast_extend(refSeq, refLen, querySeq, queryLen, initScore, endBonus,
	      alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh,
	      scoreLow, scoreHigh, swFeedback,
	      mismatchWt, gapWt, gapOpenWt, ambigWt) ;
  
  score = (scoreLow + scoreHigh) / 2 ;
  alignedQLen = (alignedQLenLow + alignedQLenHigh)/2 ;
  alignedRLen = (alignedRLenLow + alignedRLenHigh)/2 ;  

  bool partialQ = (alignedQLenLow != alignedQLenHigh) || 
    (alignedQLenLow > 0 && alignedQLenLow < queryLen && alignedRLenLow < refLen) ;
  bool partialR = (alignedRLenLow != alignedRLenHigh) ;

  if (partialQ || partialR) {

#ifdef DEBUG0
    cerr << "Partial alignment: Running ksw_extend with feedback: " << swFeedback << endl ;
#endif

    if (swFeedback.maxQLen == swFeedback.minQLen) {
      alignedQLen = swFeedback.minQLen ;
      alignedRLen = swFeedback.minRLen ;
      score = swFeedback.minScore ;
#ifdef DEBUG0
      cerr << "Skipped ksw_extend because of zero query length" << endl ;
#endif
    }
    else {
      int endBonusUpdated = (swFeedback.maxQLen == queryLen) ? endBonus : 0 ;
      score = run_ksw_extend(refSeq+swFeedback.minRLen, swFeedback.maxRLen-swFeedback.minRLen, 
			     querySeq+swFeedback.minQLen, swFeedback.maxQLen-swFeedback.minQLen,
			     swFeedback.maxBand, swFeedback.minScore, endBonusUpdated, zdrop, 
			     costMatrixRowCnt, costMatrix, gapOpenWt, gapWt,
			     alignedQLen, alignedRLen) ;

      alignedQLen += swFeedback.minQLen ;
      alignedRLen += swFeedback.minRLen ;
    }
    
#ifdef DEBUG1
    cerr << "The final alignment: " << alignedQLen << " " << alignedRLen 
	 << " with score " << score << endl ;

    int alignedQLenOrig=-1, alignedRLenOrig=-1 ;
    int scoreOrig = run_ksw_extend(refSeq, refLen, querySeq, queryLen, queryLen, initScore, endBonus,
				   zdrop, costMatrixRowCnt, costMatrix, gapOpenWt, gapWt,
				   alignedQLenOrig, alignedRLenOrig) ;
    cerr << "The original alignment was: " << alignedQLenOrig << " " << alignedRLenOrig 
	 << " with score " << scoreOrig << endl ;

    if (alignedQLenOrig != alignedQLen) 
      cerr << "Mismatch in query lengths: " << alignedQLen << " vs. " << alignedQLenOrig 
	   << " with scores " << score << " vs. " << scoreOrig << endl ;
#endif    
  }
#ifdef DEBUG1
  else {
    cerr << ((alignedQLenHigh == 0) ? "No alignment" : "Full alignment") << endl ;

    int alignedQLenOrig=-1, alignedRLenOrig=-1 ;

    int scoreOrig = run_ksw_extend(refSeq, refLen, querySeq, queryLen, queryLen, initScore, endBonus,
				   zdrop, costMatrixRowCnt, costMatrix, gapOpenWt, gapWt,
				   alignedQLenOrig, alignedRLenOrig) ;

    cerr << "The original alignment was: " << alignedQLenOrig << " " << alignedRLenOrig 
	 << " with score " << scoreOrig << endl ;

  }
#endif

  assert(alignedQLen >= 0 && alignedQLen <= queryLen) ;
  assert(alignedRLen >= 0 && alignedRLen <= refLen) ;
}



void init_fast_extend(bool avx2present) {
#ifdef FINE_ALIGN
  BitCount8::init_lut() ;
#endif
  //Set ksw_extend functions
#ifdef USE_AVX2
  if (avx2present) {
      useavx2 = avx2present;
      ksw_extend_16 = ksw_extend_avx2_16;
      ksw_extend_32 = ksw_extend_avx2_32;
  } else {
#endif
      ksw_extend_16 = ksw_extend_sse_16;
      ksw_extend_32 = ksw_extend;
#ifdef USE_AVX2
  }
#endif
}


uint8_t BitCount8::lut_[] ;
void BitCount8::init_lut() {

  for (int bitV=0; bitV < 256; ++bitV) {

    // The following is from Brian Kernighan
    // See: "Bit Twiddling Hacks"
    int cnt ;
    int val = bitV ;
    for (cnt=0; val ; ++cnt) {
      val &= val - 1 ; // clear the least significant bit set
    }
    lut_[bitV] = cnt ;
  }

}

