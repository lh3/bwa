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

#include "ed_intrav64.h"
#include "ed_intrav64x2.h"

#include "ed_intrav.h"
#include "ed_fine.h"
#include "filter.h"


//using namespace std ;

#define DIST_TYPE uint16_t

//#define MAIN
//#define DEBUG0
//#define DEBUG1
//#define DEBUG2

#define FINE_ALIGN

#define FILTER_USING_HAMMING_DISTANCE


#ifdef SW_FILTER_AND_EXTEND
#define SW_FEEDBACK 
#endif


int run_ksw_extend(uint8_t* refSeq, int refLen, uint8_t* querySeq, int queryLen,
		   int bandW, int initScore, int endBonus, int zdrop, int costMatrixRowCnt,
		   const int8_t* costMatrix, int gapo, int gape,
		   int& alignedQLen, int& alignedRLen) {

  int qle, tle, gtle, gscore, max_off ;

  int pscore = ksw_extend(queryLen, querySeq, refLen, refSeq, costMatrixRowCnt, costMatrix,
			  gapo, gape, bandW, endBonus, zdrop, initScore,
			  &qle, &tle, &gtle, &gscore, &max_off) ;

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


#ifdef SW_FEEDBACK

template<class EDVec>
void compute_sw_feedback(const EDVec& bestLowScoreVec, const EDVec& bestColVecLow,
			 const EDVec& bestHighScoreVec, const EDVec& bestColVecHigh,
			 const EDVec& queryLenVec, const EDVec& gapWtVec, 
			 const EDVec& gapOpenWtVec, const EDVec& ambigBaseCntsVec,
			 const EDVec& ambigWtP1Vec, 
			 int mismatchWtLow, int mismatchWtHigh, int validDistCnt,
			 SWFeedback& feedback) {

  EDVec deleteCntLow = bestColVecLow.subSat(queryLenVec) ;
  EDVec insertCntLow = queryLenVec.subSat(bestColVecLow) ;
  EDVec deleteCntHigh = bestColVecHigh.subSat(queryLenVec) ;
  EDVec insertCntHigh = queryLenVec.subSat(bestColVecHigh) ;

  // Recompute the real edit distance.
  // score =  - insertCnt * (gapWt+1)     
  //          - deleteCnt * gapWt         
  //          - (gapCnt > 0) ? gapOpenWt  
  //          - editDist * (mismatchWt+1) 
  //          + qlen - ambigCnt * (ambigWt+1) 
  //
  // edistDist * (mismatchWt+1) = qlen - ambigCnt * (ambigWt+1)
  //                              - deleteCnt * (gapWt+1)     
  //                              - insertCnt * gapWt    
  //                              - (gapCnt > 0) ? gapOpenWt 
  //                              - score 

  EDVec gapCntLow = insertCntLow + deleteCntLow ;
  EDVec gapCntHigh = insertCntHigh + deleteCntHigh ;
  EDVec zeroVec(0) ;

  EDVec gapPenaltyVecLow = gapCntLow * gapWtVec + insertCntLow ;
  gapPenaltyVecLow.addThirdIfFirstGTSecond(gapCntLow, zeroVec, gapOpenWtVec, zeroVec) ;

  EDVec gapPenaltyVecHigh = gapCntHigh * gapWtVec + insertCntHigh ;
  gapPenaltyVecHigh.addThirdIfFirstGTSecond(gapCntHigh, zeroVec, gapOpenWtVec, zeroVec) ;

  // Best scores were computed using the following formula:
  // score = deltaScore - gapPenalty - accumDist * mismatchWtP1 ;
  // Then:
  // accumDist * mismatchWtP1 = deltaScore - gapPenalty - score 
  //

  EDVec ambigPenaltyVec = ambigBaseCntsVec * ambigWtP1Vec ;

  EDVec weightedEDLow = queryLenVec - ambigPenaltyVec - gapPenaltyVecHigh - bestLowScoreVec ;
  EDVec weightedEDHigh = queryLenVec - ambigPenaltyVec - gapPenaltyVecLow - bestHighScoreVec ;

  //cout << queryLenVec << endl << ambigPenaltyVec << endl << gapPenaltyVecHigh << endl
  //     << bestLowScoreVec << endl << weightedEDLow << endl ;

  feedback.maxQLen = 0 ;
  feedback.maxRLen = 0 ;
  feedback.maxBand = 0 ;

  for (int iter=0; iter < 2; ++iter) {
    const EDVec& bestScoreVec = (iter == 0) ? bestLowScoreVec : bestHighScoreVec ;
    const EDVec& weightedEDVec = (iter == 0) ? weightedEDLow : weightedEDHigh ;
    const EDVec& bestColVec = (iter == 0) ? bestColVecLow : bestColVecHigh ;
    int mismatchWt = ((iter == 0) ? mismatchWtHigh : mismatchWtLow) + 1 ;

    for (int di=0; di < validDistCnt; ++di) {

      // We will consider only the entries with a non-negative score
      if (bestScoreVec.getWord(di) > 0) {

#ifdef DEBUG1
	cout << "Edit distance of " << ((iter==0) ? "low" : "high") << " score for qlen = "
	     << queryLenVec.getWord(di) << " is " << weightedEDVec.getWord(di) / mismatchWt
	     << ", where mismatchWt = " << mismatchWt << endl ;
#endif
	
	feedback.maxBand = std::max(feedback.maxBand, 
				    weightedEDVec.getWord(di) / mismatchWt) ;

	feedback.maxQLen = std::max(feedback.maxQLen, queryLenVec.getWord(di)) ;
	feedback.maxRLen = std::max(feedback.maxRLen, bestColVec.getWord(di)) ;
	
      }
    }
  }

  feedback.maxBand += 1 ; // band should be at least 1 more than the edit distance
}

#endif // #ifdef SW_FEEDBACK

template<class BitVec>
int getValueForDebug(const BitVec& B0, const BitVec& B1, int index) {
  int bit0 = (int) B0.getBit(index) ;
  int bit1 = (int) B1.getBit(index) ;
  //assert (bit0 == 0 || bit1 == 0) ;
  //assert (bit0 <= 1 && bit1 <= 1) ;

  return (bit0 == 1) ? 1 : ((bit1 == 1) ? -1 : 0) ;
  //return (bit1 << 1) | bit0 ;
}

template<class BitVec>
int getValueForDebug(const BitVec& B0, const BitVec& B1, const BitVec& B2, int index) {
  int bit0 = (int) B0.getBit(index) ;
  int bit1 = (int) B1.getBit(index) ;
  int bit2 = (int) B2.getBit(index) ;
  //assert (bit0 == 0 || bit1 == 0) ;
  //assert (bit0 <= 1 && bit1 <= 1) ;

  //return (bit0 == 1) ? 1 : ((bit1 == 1) ? -1 : 0) ;
  return (bit2 << 2) | (bit1 << 1) | bit0 ;
}


template<class BitVec>
void printVecForDebug(const BitVec& B0, const BitVec& B1, int cnt) {
  for (int i=0; i < cnt; ++i) {
    cout << getValueForDebug(B0, B1, i) ;
    if (i % 8 == 7)
      cout << "|" ;
    if (i % 16 == 15)
      cout << "|" ;
    if (i % 32 == 31)
      cout << "|" ;
    if (i % 64 == 63)
      cout << "|" ;
    cout << " " ;
  }
}

template<class BitVec>
void printVecForDebug(const char* msg, const BitVec& B0, int cnt) {

  cout << msg ;
  for (int i=0; i < cnt; ++i) {
    cout << B0.getBit(i) << " " ;
  }
  cout << endl ;
}




template<class BitVec>
void printVecForDebug(const BitVec& B0, const BitVec& B1, const BitVec& B2, int cnt) {
  for (int i=0; i < cnt; ++i) {
    cout << getValueForDebug(B0, B1, i) << " " ;
  }
}


template<class BitVec>
void printVecsForDebug(const BitVec& L0, const BitVec& L1,
		       const BitVec& U0, const BitVec& U1, int cnt) {

  cout << "L: " ;
  printVecForDebug(L0, L1, cnt) ;
  cout << endl ;

  cout << "U: " ;
  printVecForDebug(U0, U1, cnt) ;
  cout << endl ;
}

inline void orBitAtPosition(uint64_t& w, bool newBit, int index) {
  w |= (((uint64_t)newBit) << index) ;
}

template<class BitVec>
inline BitVec computeMatchVec(const BitVec& Q0, const BitVec& Q1, const BitVec& Q2,
			      const BitVec& R0, const BitVec& R1, const BitVec& R2) {

#ifdef MAIN
  return Q2 | R2 | (Q1 ^ R1) | (Q0 ^ R0) ;
#else
  return ~Q2 & ~R2 & ((Q1 ^ R1) | (Q0 ^ R0)) ;
#endif
}

/*
template<class BitVec, class EDVec>
class ScoreProbe {
  DistVec<BitVec,EDVec> accumDist ;
  EDVec ambigBaseCnts ;
  
  EDVec bestLowScoreVec, bestHighScoreVec, bestColVecLow, bestColVecHigh, currNegLowScoreVec,
    currNegHighScoreVec, currColVec ;
  EDVec queryLenVec ;
} ;
*/

int computeHammingDistance(const uint8_t* refSeq, int refLen, uint8_t* querySeq, int queryLen, int upperBound, int& lastConsecutiveMatching) {
  
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


template<class BitVec, class EDVec>
DIST_TYPE edit_dist(uint8_t* refSeq, int refLen, uint8_t* querySeq, int queryLen,
		    int initScore, int endBonus, 
		    int& alignedQLenLow, int& alignedQLenHigh,
		    int& alignedRLenLow, int& alignedRLenHigh,
		    int& maxScoreLow, int& maxScoreHigh, SWFeedback& swfb,
		    bool printDebug = false) {

  assert (BitVec::isQueryLengthOk(queryLen)) ;

#ifdef DEBUG0
  cout << "Query and ref lens for SW problem: " << queryLen << " " << refLen << endl ;
#endif

#ifdef FILTER_USING_HAMMING_DISTANCE
  
  if (queryLen >= 20) {
    int lastConsecutiveMatching=0 ;
    int upperBound = queryLen / 6 ;

    int hdist = computeHammingDistance(refSeq, refLen, querySeq, queryLen, upperBound, lastConsecutiveMatching) ;
    
    if (hdist < upperBound && lastConsecutiveMatching > 6) {
#ifdef DEBUG0
      cout << "The Hamming distance is " << hdist  << " with "
	   << lastConsecutiveMatching
	   << " bases matches at the end of query. " << endl ;
      cout << "Returning full alignment directly." << endl ;
#endif
      alignedQLenLow = alignedQLenHigh = queryLen ;
      alignedRLenLow = alignedRLenHigh = queryLen ;
      maxScoreLow = queryLen - hdist - hdist * 6 ;
      maxScoreHigh = queryLen - hdist - hdist * 4 ;

      return hdist ;
    }
  }
  
  


#endif // #ifdef FILTER_USING_HAMMING_DISTANCE


  BitVec L0, L1, U0, U1, UP0, UP1 ;
  BitVec Q0, Q1, Q2 ;
  BitVec D[5] ; // profile vecs

#ifdef MAIN
  bool initL0 = 1 ;
  bool initL1 = 0 ;

  bool initU0 = 1 ;
  bool initU1 = 0 ;
#else
  bool initL0 = 0 ;
  bool initL1 = 0 ;

  bool initU0 = 0 ;
  bool initU1 = 0 ;
#endif


#ifdef FINE_ALIGN
  const int FINE_ALIGN_CNT = 16 ; // align the first 16 bases in a fine-grain way
  assert (FINE_ALIGN_CNT <= 16) ; // the rest of the code assumes this.
  int ambigBaseCntsFineArr[FINE_ALIGN_CNT] ;
  for (int i=0; i < FINE_ALIGN_CNT; ++i)
    ambigBaseCntsFineArr[i] = 0 ;
#endif

  EDVec ambigBaseCnts(0) ;
  int ambigBaseCntsArr[EDVec::WORD_CNT()] ;
  for (int i=0; i < EDVec::WORD_CNT(); ++i)
    ambigBaseCntsArr[i] = 0 ;
  
  int ambigWordIndex = 0 ;
  //int totalAmbigBaseCnt = 0 ;
  //int lastWordUpdated = 0 ;
  // Store the query in bit vectors
  for (int wIndex=0; wIndex * 63 < queryLen; ++wIndex) {
    uint64_t bit0Partial = 0x0 ;
    uint64_t bit1Partial = 0x0 ;
    uint64_t bit2Partial = 0x0 ;

    for (int bIndex=0; bIndex < 63; ++bIndex) {
      int qi = wIndex * 63 + bIndex ;
      if (qi >= queryLen)
	break ;
      uint8_t qbase = querySeq[qi] ;
      if (qbase == 4) {
	//lastWordUpdated = EDVec::getWordIndexFor(qi) ;
	//ambigBaseCntsArr[lastWordUpdated] = ++totalAmbigBaseCnt ;
	//++ambigBaseCntsArr[EDVec::getWordIndexFor(qi)] ;

	while (qi >= EDVec::getDistAt(ambigWordIndex, queryLen))
	  ++ambigWordIndex ;
	++ambigBaseCntsArr[ambigWordIndex] ;

	//cout << "Ambig base at " << wIndex << " " << bIndex << " " << qi << endl ;
	//cout << "Word index for ambig base: " << ambigWordIndex << endl ;

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
  EDVec bestAccumDistLow = accumDist.getVec() ;
  EDVec bestAccumDistHigh = accumDist.getVec() ;

  int INF_SCORE = EDVec::getMaxWordVal() ;

  EDVec currColVec(0) ;
  EDVec bestLowScore(0), bestHighScore(0), bestColVecLow(0), bestColVecHigh(0) ;

  EDVec queryLenVec ;
  queryLenVec.setWordsAsDist(queryLen, 1) ;

  int mismatchWtLow = 3; // 3 ;
  int mismatchWtHigh = 6 ; // 6 ;
  int gapWt = 1 ;
  int gapOpenWt = 6 ; // 2
  int ambigWt = 1 ;

#ifdef FINE_ALIGN
  FineAlignment16 fineAlign(queryLen, endBonus, mismatchWtLow, mismatchWtHigh, ambigWt) ;
#endif

  EDVec mismatchWtLowVec(mismatchWtLow), mismatchWtHighVec(mismatchWtHigh),
    mismatchWtLowP1Vec(mismatchWtLow+1), mismatchWtHighP1Vec(mismatchWtHigh+1), // P1: plus 1
    gapWtVec(gapWt), gapWtP1Vec(gapWt+1), gapOpenWtVec(gapOpenWt), ambigWtP1Vec(ambigWt+1), 
    zeroVec(0) ;
  
  initScore -= (6-gapOpenWt) ;
  initScore += endBonus ;

  EDVec endBonusVec(0), deltaScoreVec(initScore), badScoreVec(0) ;
  endBonusVec.setWordsAsEndBonus(queryLen, endBonus) ;
  deltaScoreVec = deltaScoreVec + queryLenVec - ambigBaseCnts * ambigWtP1Vec + endBonusVec ;
  badScoreVec.setWordsAsBadScore(queryLen, 0, INF_SCORE) ;

  EDVec deltaCostInsert((gapWt+1)) ; // insert costs decrease as ci increases
  EDVec deltaCostDelete(0) ; // delete costs increase as ci increases
  EDVec deltaCostDistWtLow(-(mismatchWtLow+1)), deltaCostDistWtHigh(-(mismatchWtHigh+1)) ;
  EDVec deltaColVec(1) ;
  
  EDVec insertVecModifier(0), deleteVecModifier(0), gapOpenDelta(0) ;
  insertVecModifier.setFirstWord(-(gapWt+1)) ; // will cancel one entry when added
  deleteVecModifier.setFirstWord(-gapWt) ; // will add deleteWt to one entry
  gapOpenDelta.setFirstWord(gapOpenWt) ; // will remove/add gapOpenWt from/to one diag entry

  EDVec firstColCost = deltaScoreVec - queryLenVec * gapWtP1Vec - gapOpenWtVec ;

  EDVec currLowScore = firstColCost ;
  EDVec currHighScore = firstColCost ;
  EDVec deltaMismatch(0) ;

#ifdef DEBUG1
  cout << " The init score is: " << initScore << endl ;
  cout << " The cost of the 0th column: " << firstColCost << endl ;
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
    cout << "In iteration " << __ci << ", the cost values are updated to:" << endl ; \
    cout << "L0: " << std::hex << L0 << std::dec << endl ;		\
    cout << "L1: " << std::hex << L1 << std::dec << endl ;		\
    cout << "currColVec: " << currColVec << endl ;			\
    cout << "queryLenVec: " << queryLenVec << endl ;			\
    cout << "deltaCostInsert: " << deltaCostInsert << endl ;	        \
    cout << "deltaCostDelete: " << deltaCostDelete << endl ;      	\
    cout << "deltaMismatch: " << deltaMismatch << endl ;		\
    cout << "deltaCostDistWtHigh: " << deltaCostDistWtHigh << endl ;	\
    cout << "deltaCostDistWtLow: " << deltaCostDistWtLow << endl ;	\
    cout << "currLowScore: " << currLowScore << endl ;			\
    cout << "currHighScore: " << currHighScore << endl ;		\
    cout << "bestColLow: " << bestColVecLow << endl ;			\
    cout << "bestColHigh: " << bestColVecHigh << endl ;			\
    cout << "bestLowScore: " << bestLowScore << endl ;			\
    cout << "bestHighScore: " << bestHighScore << endl ;		\
  }

#define PRINT_DELTA_MODIFICATIONS_DEBUG(__ci) {			       \
    cout << "The delta costs have been updated after processing column " << __ci << endl ; \
    cout << "New delta cost insert: " << deltaCostInsert << endl ;	\
    cout << "New delta cost delete: " << deltaCostDelete << endl ;	\
    cout << "using the modifiers: " << endl ;				\
    cout << "Insert vec modifier: " << insertVecModifier << endl ;	\
    cout << "Delete vec modifier: " << deleteVecModifier << endl ;	\
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
    currLowScore += gapDelta ;						\
    currHighScore += gapDelta ;						\
									\
    currLowScore += (deltaMismatch * deltaCostDistWtHigh) ;		\
    currHighScore += (deltaMismatch * deltaCostDistWtLow) ;		\
    bestLowScore.setMax(currLowScore, bestColVecLow, currColVec) ;	\
    bestHighScore.setMax(currHighScore, bestColVecHigh, currColVec) ;	\
  }

#define UPDATE_COSTS_AND_CHECK(__ci) {					\
    currColVec += deltaColVec ;						\
    EDVec gapDelta = deltaCostInsert + deltaCostDelete ;		\
    /*gapDelta.addThirdIfFirstGTSecond(gapDelta, zeroVec, gapOpenWtVec, zeroVec) ; */ \
    currLowScore += gapDelta ;						\
    currHighScore += gapDelta ;						\
									\
    /* Add delta costs corresponding to the current mismatch delta (-1, 0, or 1) */  \
    deltaMismatch = accumDist.getDeltaDistVec(L0, L1) ;		\
									\
    currLowScore += (deltaMismatch * deltaCostDistWtHigh) ;		\
    currHighScore += (deltaMismatch * deltaCostDistWtLow) ;		\
    updateFlag |= bestLowScore.setMaxAndReturnFlag(currLowScore, bestColVecLow, currColVec) ; \
    updateFlag |= bestHighScore.setMaxAndReturnFlag(currHighScore, bestColVecHigh, currColVec) ; \
  }

#define UPDATE_COSTS_AND_STORE_BEST_ACCUM_DIST(__ci) {			\
    currColVec += deltaColVec ;						\
    EDVec gapDelta = deltaCostInsert + deltaCostDelete ;		\
    /*gapDelta.addThirdIfFirstGTSecond(gapDelta, zeroVec, gapOpenWtVec, zeroVec) ; */ \
    currLowScore += gapDelta ;						\
    currHighScore += gapDelta ;						\
									\
    /* Add delta costs corresponding to the current mismatch delta (-1, 0, or 1) */  \
    deltaMismatch = accumDist.getDeltaDistVec(L0, L1) ;		\
									\
    currLowScore += (deltaMismatch * deltaCostDistWtHigh) ;		\
    currHighScore += (deltaMismatch * deltaCostDistWtLow) ;		\
    bestLowScore.setMax(currLowScore, bestColVecLow, currColVec, bestAccumDistLow, accumDist.getVec()) ; \
    bestHighScore.setMax(currHighScore, bestColVecHigh, currColVec, bestAccumDistHigh, accumDist.getVec()) ; \
  }

#define UPDATE_COSTS_AND_CHECK_AND_STORE_BEST_ACCUM_DIST(__ci) {	\
    currColVec += deltaColVec ;						\
    EDVec gapDelta = deltaCostInsert + deltaCostDelete ;		\
    /*gapDelta.addThirdIfFirstGTSecond(gapDelta, zeroVec, gapOpenWtVec, zeroVec) ; */ \
    currLowScore += gapDelta ;						\
    currHighScore += gapDelta ;						\
									\
    /* Add delta costs corresponding to the current mismatch delta (-1, 0, or 1) */  \
    deltaMismatch = accumDist.getDeltaDistVec(L0, L1) ;		\
									\
    currLowScore += (deltaMismatch * deltaCostDistWtHigh) ;		\
    currHighScore += (deltaMismatch * deltaCostDistWtLow) ;		\
    updateFlag |= bestLowScore.setMaxAndReturnFlag(currLowScore, bestColVecLow, currColVec, bestAccumDistLow, accumDist.getVec()) ; \
    updateFlag |= bestHighScore.setMaxAndReturnFlag(currHighScore, bestColVecHigh, currColVec, bestAccumDistHigh, accumDist.getVec()) ; \
  }


#define REMOVE_GAP_OPEN_PENALTY_FOR_DIAG(__ci) {			\
    currLowScore += gapOpenDelta ;					\
    currHighScore += gapOpenDelta ;					\
  }

#define ADD_GAP_OPEN_PENALTY_TO_DIAG(__ci) {				\
    currLowScore -= gapOpenDelta ;					\
    currHighScore -= gapOpenDelta ;					\
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
    if (currHighScore.allLessThanOrEqualTo(badScoreVec)) {	\
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
    EDVec currLowScoreVecTest =						\
      deltaScoreVec - gapPenaltyVec - accumDist.getVec() * mismatchWtHighP1Vec ; \
    EDVec currHighScoreVecTest =							\
      deltaScoreVec - gapPenaltyVec - accumDist.getVec() * mismatchWtLowP1Vec ;	\
									\
    if (!(currLowScoreVecTest == currLowScore) || !(currHighScoreVecTest == currHighScore)) { \
      cout << "Error: Mismatch in cost values computed using two methods at column " << ci << endl ; \
      cout << "deleteCntVec: " << deleteCntVec << endl ;		\
      cout << "insertCntVec: " << insertCntVec << endl ;		\
      cout << "gapPenaltyVec: " << gapPenaltyVec << endl ;		\
      cout << "low score (test): " << currLowScoreVecTest << endl ;	\
      cout << "high score (test): " << currHighScoreVecTest << endl ;	\
      exit(0) ;								\
    }									\
  }

  uint16_t flagMask = (1 << (EDVec::getLastWordIndexFor(queryLen)+1)) - 1 ;
  int firstProbeOffset = accumDist.getProbeColumn(0) ;
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

  numFullIters = (min(queryLen, refLen)-ci+1) /EDVec::PERIOD() ;
  for (int iterIndex = 0; iterIndex < numFullIters; ++iterIndex) {

    //cout << "OUTER ITERATION: " << iterIndex << " of " << numFullIters << endl ;

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


  // If there are 3 words in EDVec, flagMask will be 0....0111 in binary
  
  // Next, process the columns with indices greater than or equal to query length
  numFullIters = (refLen-ci) / EDVec::PERIOD() ;
  for (int iterIndex = 0; iterIndex < numFullIters; ++iterIndex) {
    uint16_t updateFlag = 0x0 ;

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

    //cout << ci << " " << queryLen << " " << refLen << " " 
    //	 << std::hex << updateFlag << std::dec << endl ;
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
    cout << "Loop terminated after iteration " << ci << endl ;
#endif
  
  //maxScoreHigh = initScore ;
  //maxScoreLow = initScore ;

  //const int SCORE_INF = 10000 ;
  //maxScoreHigh = -SCORE_INF ;
  //maxScoreLow = -SCORE_INF ;

  maxScoreHigh = initScore ;
  maxScoreLow = initScore ;

  alignedQLenLow = 0 ;
  alignedQLenHigh = 0 ;
  alignedRLenLow = 0 ;
  alignedRLenHigh = 0 ;

#ifdef FINE_ALIGN
  int alignedLenLowFine, alignedLenHighFine, scoreLowFine, scoreHighFine ;
  fineAlign.getBestLow(alignedLenLowFine, scoreLowFine) ;
  fineAlign.getBestHigh(alignedLenHighFine, scoreHighFine) ;

  scoreLowFine += initScore ;
  scoreHighFine += initScore ;

  if (scoreLowFine > maxScoreLow) {
    maxScoreLow = scoreLowFine ;
    alignedQLenLow = alignedLenLowFine ;
    alignedRLenLow = alignedLenLowFine ;
  }

  if (scoreHighFine > maxScoreHigh) {
    maxScoreHigh = scoreHighFine ;
    alignedQLenHigh = alignedLenHighFine ;
    alignedRLenHigh = alignedLenHighFine ;
  }

#ifdef DEBUG0 
  cout << "The fine-grain alignment: " << fineAlign << endl ;
#endif

#endif // #ifdef FINE_ALIGN


#ifdef SW_FEEDBACK
  swfb = SWFeedback(alignedQLenHigh, alignedRLenHigh, 1) ;
  const int scoreTolerance = EDVec::PERIOD() ;
#endif

  for (int di=0; di < accumDist.getValidDistCnt(); ++di) {
  //for (int di=0; di < 1; ++di) {

    int currQLen = queryLenVec.getWord(di) ;
    int currRLenLow = bestColVecLow.getWord(di) ;
    int currRLenHigh = bestColVecHigh.getWord(di) ;

    // if (currQLen < 5) {
    //  continue ; // skip the partial alignments that are too short
    //}
    
    int currLowScore = bestLowScore.getWord(di) ; 
    //currLowScore -= (currQLen > 16) ? 5 : 0 ;

    if (currLowScore > maxScoreLow) {
      maxScoreLow = currLowScore ;
      alignedQLenLow = currQLen ;
      alignedRLenLow = currRLenLow ;
    }

    int currHighScore = bestHighScore.getWord(di) ;//+ ((currQLen == queryLen) ? endBonus : 0) ;
    //currHighScore += (currQLen > 16) ? 5 : 0 ;

    if (currHighScore > maxScoreHigh) {
      maxScoreHigh = currHighScore ;
      alignedQLenHigh = currQLen ; 
      alignedRLenHigh = currRLenHigh ; 
    }

#ifdef SW_FEEDBACK

    int currDistLow = bestAccumDistLow.getWord(di) ;
    int currDistHigh = bestAccumDistHigh.getWord(di) ;
    int currGapLow = (currQLen > currRLenLow) ? (currQLen - currRLenLow) : (currRLenLow - currQLen) ;
    int currGapHigh = (currQLen > currRLenHigh) ? (currQLen - currRLenHigh) : (currRLenHigh - currQLen) ;

    int currQLenUpdated = min(queryLen, currQLen+EDVec::PERIOD()-1) ;

    int currRLenLowUpdated = min(refLen, max(currRLenLow, currQLenUpdated)+currDistLow) ;
    int currRLenHighUpdated = min(refLen, max(currRLenHigh, currQLenUpdated)+currDistHigh) ;



    if (currLowScore > maxScoreLow - scoreTolerance) {
      swfb.updateMax(currQLenUpdated, currRLenLowUpdated, currDistLow + currGapLow + 1) ;
    }
    if (currHighScore > maxScoreHigh - scoreTolerance) {
      swfb.updateMax(currQLenUpdated, currRLenHighUpdated, currDistHigh + currGapHigh + 1) ;
    }
#endif

#ifdef DEBUG0
    cout << "For partialQLen " << queryLenVec.getWord(di) << ": "
	 << "score = [" << currLowScore << ", " << currHighScore << "], "
 	 << "partial ref len = [" << bestColVecLow.getWord(di) << ", "
	 << bestColVecHigh.getWord(di) << "], "
	 << "qlen in vector = " << queryLenVec.getWord(di) << " "
	 << "ambigBaseCnt = " << ambigBaseCnts.getWord(di) << ", "
#ifdef SW_FEEDBACK
	 << "editDist = [" << currDistLow << ", " << currDistHigh << "]"
#endif
	 << endl ;

#ifdef SW_FEEDBACK
    cout << "The feedback to be sent: " << swfb << endl ;
#endif

#endif
  }

  
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
    cout << "Orig/Back-by-1/Forward-by-1 match counts: " << origMatchCnt << " "
	 << back1MatchCnt << " " << forw1MatchCnt << endl ;
    cout << "The aligned ref len is now: " << alignedRLenHigh << endl ;
#endif

  }


#ifdef DEBUG0
  cout << "Init score = " << initScore << endl ;
  cout << "Max score-high = " << maxScoreHigh << " for alignedQLen = " 
       << alignedQLenHigh << ", alignedRLen = " << alignedRLenHigh << endl ;
  cout << "Max score-low = " << maxScoreLow << " for alignedQLen = " 
       << alignedQLenLow << ", alignedRLen = " << alignedRLenLow << endl ;
#endif

#if 0
#ifdef SW_FEEDBACK
  if (!((alignedQLenHigh == queryLen && alignedQLenLow == queryLen) ||
	(alignedQLenHigh == 0 && alignedQLenLow == 0))) {

    // TODO: Also consider FINE_ALIGN
    compute_sw_feedback<EDVec>(bestLowScore, bestColVecLow, bestHighScore, bestColVecHigh,
			       queryLenVec, gapWtVec, gapOpenWtVec, ambigBaseCnts,
			       ambigWtP1Vec, mismatchWtLow, mismatchWtHigh, 
			       accumDist.getValidDistCnt(), swFeedback) ;

#ifdef DEBUG0
    cout << "The feedback sent to SW: " << swFeedback.maxQLen << " " << swFeedback.maxRLen
	 << " " << swFeedback.maxBand << endl ;
#endif

  }
#endif
#endif

  return accumDist.getTotalDist() ;
}


bool unit_test1() {

  uint8_t ref[] = {4, 4, 1, 3, 2, 2} ;
  uint8_t query[] = {4, 3, 0, 2, 2} ;

  //uint8_t ref[] = {1, 3, 3, 2, 2} ;
  //uint8_t query[] = {3, 3, 2, 2} ;

  //uint8_t ref[]   = {3, 3, 2, 2} ;
  //uint8_t query[] = {3, 3, 2, 2} ;

  //uint8_t ref[] = {2} ;
  //uint8_t query[] = {2} ;

  int alignedQLenLow=0, alignedQLenHigh=0, alignedRLenLow=0, alignedRLenHigh=0 ;
  int scoreLow, scoreHigh ;
  SWFeedback swFeedback ;
  int dist = edit_dist<BitVec64, EDVec32Every4>(ref, 6, query, 5, 0, 0, alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh, scoreLow, scoreHigh, swFeedback, false) ;

  cout << "Unit test1 result: " << dist << endl ;

  return dist == 4 ;
}

bool unit_test2() {

  //uint8_t query[] = {3,1,3,0,2,0,2,2,0,2,2,0,0} ;
  //uint8_t ref[]   = {0,1,3,0,2,0,2,2,0,2,2,0,0} ;

  uint8_t query[] = {3,1,3} ; 
  uint8_t ref[]   = {0,1,3} ; 

  int alignedQLenLow=0, alignedQLenHigh=0, alignedRLenLow=0, alignedRLenHigh=0 ;
  int scoreLow, scoreHigh ;
  SWFeedback swFeedback ;
  int dist = edit_dist<BitVec64, EDVec32Every4>(ref, 3, query, 3, 0, 0, alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh, scoreLow, scoreHigh, swFeedback, false) ;
  cout << "Unit test2 result: " << dist << endl ;

  return dist == 1 ;

}

bool unit_test3() {

  uint8_t query[] = {2,0,3,3,3,1,3,3,1,0,2,3,0,2,2,0,1,1,2,0,2,3,1} ;
  uint8_t ref[] =   {3,0,3,3,3,1,3,3,1,0,2,3,0,0,2,0,1,0,2,0,2,1,1} ;
  // qlen = rlen = 23

  //uint8_t query[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0} ;
  //uint8_t ref[] =   {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0} ;
  // qlen = rlen = 17

  int alignedQLenLow=0, alignedQLenHigh=0, alignedRLenLow=0, alignedRLenHigh=0 ;
  int scoreLow, scoreHigh ;
  SWFeedback swFeedback ;
  int dist = edit_dist<BitVec64x2, EDVec128Every16>(ref, 23, query, 23, 21, 5, alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh, scoreLow, scoreHigh, swFeedback, false) ; 

  cout << "Unit test3 result: " << dist << endl ;

  return true ;
}


DIST_TYPE edit_dist_u32(uint8_t* refSeq, int refLen, uint8_t* querySeq, int queryLen, 
			int initScore, int endBonus, 
			int& alignedQLenLow, int& alignedQLenHigh, 
			int& alignedRLenLow, int& alignedRLenHigh, 
			int& scoreLow, int& scoreHigh, SWFeedback& swFeedback,
			bool printDebug=false) {

  return edit_dist<BitVec64, EDVec32Every4>(refSeq, refLen, querySeq, queryLen, 
  //return edit_dist<BitVec64, EDVec128Every16>(refSeq, refLen, querySeq, queryLen, 
					    initScore, endBonus, 
					    alignedQLenLow, alignedQLenHigh, 
					    alignedRLenLow, alignedRLenHigh, 
					    scoreLow, scoreHigh, 
					    swFeedback, printDebug) ;
}


DIST_TYPE edit_dist_u64(uint8_t* refSeq, int refLen, uint8_t* querySeq, int queryLen, 
			int initScore, int endBonus, 
			int& alignedQLenLow, int& alignedQLenHigh, 
			int& alignedRLenLow, int& alignedRLenHigh, 
			int& scoreLow, int& scoreHigh, SWFeedback& swFeedback,
			bool printDebug=false) {

  return edit_dist<BitVec64, EDVec64Every8>(refSeq, refLen, querySeq, queryLen, 
  //return edit_dist<BitVec64, EDVec128Every16>(refSeq, refLen, querySeq, queryLen, 
					    initScore, endBonus, 
					    alignedQLenLow, alignedQLenHigh, 
					    alignedRLenLow, alignedRLenHigh, 
					    scoreLow, scoreHigh, 
					    swFeedback, printDebug) ;
}


DIST_TYPE edit_dist_u64x2(uint8_t* refSeq, int refLen, uint8_t* querySeq, int queryLen,
			  int initScore, int endBonus, 
			  int& alignedQLenLow, int& alignedQLenHigh, 
			  int& alignedRLenLow, int& alignedRLenHigh, 
			  int& scoreLow, int& scoreHigh, SWFeedback& swFeedback,
			  bool printDebug=false) {

  return edit_dist<BitVec64x2, EDVec128Every16>(refSeq, refLen, querySeq, queryLen, 
  //return edit_dist<BitVec64x2, EDVec64x2>(refSeq, refLen, querySeq, queryLen, 
						initScore, endBonus, 
						alignedQLenLow, alignedQLenHigh, 
						alignedRLenLow, alignedRLenHigh, 
						scoreLow, scoreHigh, 
						swFeedback, printDebug) ;
}

void edit_dist(uint8_t* refSeq, int refLen, uint8_t* querySeq, int queryLen,
		  int initScore, int endBonus, 
		  int& alignedQLenLow, int& alignedQLenHigh, 
		  int& alignedRLenLow, int& alignedRLenHigh, 
		  int& scoreLow, int& scoreHigh, SWFeedback& swFeedback,
		  bool printDebug=false) {
  
  if (queryLen < 15) {
    cout << "Warning: It is not recommended to use SW filtering API for queries shorter than 15 bases!" << endl ;
  }

  if (queryLen <= 32)
    edit_dist_u32(refSeq, refLen, querySeq, queryLen, initScore, endBonus,
		  alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh,
		  scoreLow, scoreHigh, swFeedback, printDebug) ;
  else if (queryLen <= 63)
    edit_dist_u64(refSeq, refLen, querySeq, queryLen, initScore, endBonus,
		  alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh,
		  scoreLow, scoreHigh, swFeedback, printDebug) ;
  else if (queryLen <= 126)
    edit_dist_u64x2(refSeq, refLen, querySeq, queryLen, initScore, endBonus,
		    alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh,
		    scoreLow, scoreHigh, swFeedback, printDebug) ;
  else {
    alignedQLenLow = alignedRLenLow = 0 ;
    alignedQLenHigh = queryLen ;
    alignedRLenHigh = refLen ;
    scoreLow = 0 ;
    scoreHigh = queryLen ;
    cout << "Warning: Skipping edit distance computation because reads longer than 126 bases not supported yet" << endl ;
  }
  
}


void extend_with_edit_dist(uint8_t* refSeq, int refLen, uint8_t* querySeq, int queryLen,
			   int initScore, int endBonus, 
			   int& alignedQLen, int& alignedRLen, int& score,
			   float& confidence) {
  
  int alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh, scoreLow, scoreHigh ;
  SWFeedback swFeedback ;

  edit_dist(refSeq, refLen, querySeq, queryLen, initScore, endBonus,
	    alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh,
	    scoreLow, scoreHigh, swFeedback) ;
  
  alignedQLen = (alignedQLenLow + alignedQLenHigh)/2 ;
  alignedRLen = (alignedRLenLow + alignedRLenHigh)/2 ;
  score = (scoreLow + scoreHigh)/2 ;
  confidence = 
    (alignedQLenLow == alignedQLenHigh && alignedRLenLow == alignedRLenHigh) ? 1.0 : 0.0 ;

}

void filter_and_extend(uint8_t* refSeq, int refLen, uint8_t* querySeq, int queryLen,
		       int initScore, int endBonus, int zdrop,
		       int& alignedQLen, int& alignedRLen, int& score) {

  int alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh, scoreLow, scoreHigh ;
  SWFeedback swFeedback ;

  edit_dist(refSeq, refLen, querySeq, queryLen, initScore, endBonus,
	    alignedQLenLow, alignedQLenHigh, alignedRLenLow, alignedRLenHigh,
	    scoreLow, scoreHigh, swFeedback) ;
  
  score = (scoreLow + scoreHigh) / 2 ;
  alignedQLen = (alignedQLenLow + alignedQLenHigh)/2 ;
  alignedRLen = (alignedRLenLow + alignedRLenHigh)/2 ;  
  
  bool partialQ = (alignedQLenLow != alignedQLenHigh) || 
    (alignedQLenLow > 0 && alignedQLenLow < queryLen) ;
  bool partialR = (alignedRLenLow != alignedRLenHigh) ;

  if (partialQ || partialR) {
    int endBonusUpdated = (swFeedback.maxQLen == queryLen) ? endBonus : 0 ;
    score = run_ksw_extend(refSeq, swFeedback.maxRLen, querySeq, swFeedback.maxQLen,
			   swFeedback.maxBand, initScore, endBonusUpdated, zdrop, 
			   costMatrixRowCnt, costMatrix, gapo, gape,
			   alignedQLen, alignedRLen) ;
  }

}



void init_ed_dist() {

#ifdef FINE_ALIGN
  BitCount8::init_lut() ;
#endif
  
}


#ifdef MAIN
int main() {
  
  init_ed_dist() ;



  if (!unit_test1()) {
    cout << "Failed unit test 1" << endl ;
    exit(0) ;
  }

  if (!unit_test2()) {
    cout << "Failed unit test 2" << endl ;
    exit(0) ;
  }
  
  unit_test3() ;

  return 0 ;
}
#endif

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

