/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <iostream>
#include <iomanip>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include "extend_vec128.h"
#include "extend_vec128x2.h"

#ifdef USE_AVX2
#include "extend_vec256.h"
#include "extend_vec256x2.h"
#endif

//#define EXT_AVX2_MAIN
//#define DEBUG

using std::cout ;
using std::endl ;

const int AMBIG_BASE = 4 ;
const int INVALID_BASE_Q = 5 ;
const int INVALID_BASE_R = 6 ;

class ExtendProblem {

  const uint8_t* query_ ;
  int queryLen_ ;
  const uint8_t* ref_ ;
  int refLen_ ;
  int width_ ;
  int initScore_ ;

public:

  ExtendProblem(const uint8_t* query, int queryLen, const uint8_t* ref, int refLen, 
		int width, int initScore) :
    query_(query), queryLen_(queryLen), ref_(ref), refLen_(refLen),
    width_(width), initScore_(initScore) {}
  
  inline uint8_t getQ(int index) {
    return (index >= queryLen_) ? INVALID_BASE_Q : query_[index] ;
  }

  inline uint8_t getR(int index) {
    return (index >= refLen_) ? INVALID_BASE_R : ref_[index] ;
  }

} ;



// To prepare the vecs that hold the next query/ref chars  
// nextQV: ..................... [qi+2] [qi+1] [qi]
// nextRV: [ri] [ri+1] [ri+2] .....................
#define POPULATE_NEXTQV() {					\
    for (int vi=0; vi < EC::MAX_VEC_LEN; ++vi) {		\
      arr[vi] = prob.getQ(qi++) ;				\
    }								\
    nextQV = EC::vec_set(arr) ;					\
  }

#define POPULATE_NEXTRV() {					\
    for (int vi = EC::MAX_VEC_LEN-1; vi >= 0; --vi) {		\
      arr[vi] = prob.getR(ri++) ;				\
    }								\
    nextRV = EC::vec_set(arr) ;					\
  }								\

#define POPULATE_NEXTQV_VAR(__cnt) {				\
    for (int vi=0; vi < __cnt; ++vi) {				\
      arr[vi] = prob.getQ(qi++) ;				\
    }								\
    for (int vi=__cnt; vi < EC::MAX_VEC_LEN; ++vi) {		\
      arr[vi] = INVALID_BASE_Q ;				\
    }								\
    nextQV = EC::vec_set(arr) ;					\
  }

#define POPULATE_NEXTRV_VAR(__cnt) {					\
    for (int vi = EC::MAX_VEC_LEN-1; vi >= (EC::MAX_VEC_LEN-__cnt) ; --vi) { \
      arr[vi] = prob.getR(ri++) ;					\
    }									\
    for (int vi = EC::MAX_VEC_LEN-__cnt-1; vi >= 0; --vi) {		\
      arr[vi] = INVALID_BASE_R ;					\
    }									\
    nextRV = EC::vec_set(arr) ;						\
  }									\


#define COMPUTE_SCOREV() {						\
    SCOREV = EC::vec_blend(MISMATCHSV, MATCHSV, EC::vec_compare_eq(RV, QV)) ; \
    SCOREV = EC::vec_blend(SCOREV, AMBIGSV,															\
													 EC::vec_or(EC::vec_compare_eq(RV, AMBIGBV),	\
																			EC::vec_compare_eq(QV, AMBIGBV))) ;	\
    SCOREV = EC::vec_blend(SCOREV, INVALIDSV,														\
													 EC::vec_or(EC::vec_compare_eq(RV, INVALIDRBV),	\
																			EC::vec_compare_eq(QV, INVALIDQBV))) ; \
		}

#define COMPUTE_B() {							\
    BH = EC::vec_max(EC::vec_add(AH, DESV), EC::vec_add(AD, DOSV)) ;	\
    BV = EC::vec_max(EC::vec_add(EC::vec_shift_left_and_insert(AV, BADV), \
				 IESV),					\
		     EC::vec_add(EC::vec_shift_left_and_insert(AD, BADV), \
				 IOSV)) ;				\
    BD = EC::vec_max(EC::vec_add(BD, SCOREV), EC::vec_max(BH, BV)) ;	\
  }

#define COMPUTE_A() {							\
    AV = EC::vec_max(EC::vec_add(BV, IESV), EC::vec_add(BD, IOSV)) ;	\
    AH = EC::vec_max(EC::vec_add(EC::vec_shift_right_and_insert(BH, BADV), \
				 DESV),					\
		     EC::vec_add(EC::vec_shift_right_and_insert(BD, BADV), \
				 DOSV)) ;				\
    AD = EC::vec_max(EC::vec_add(AD, SCOREV), EC::vec_max(AH, AV)) ;	\
  }

#define UPDATE_BEST_B() {						\
    typename EC::Vec __MASK = EC::vec_compare_gt(BD, MAX_B) ;		\
    MAX_B = EC::vec_blend(MAX_B, BD, __MASK) ;				\
    BEST_B_DIAGV = EC::vec_blend(BEST_B_DIAGV, DIAGV, __MASK) ;		\
  }

#define UPDATE_BEST_A() {						\
    typename EC::Vec __MASK = EC::vec_compare_gt(AD, MAX_A) ;		\
    typename EC::Vec __MASK2 = EC::vec_compare_gt(ZEROV, AD) ;		\
    terminate = EC::vec_test_all_ones(__MASK2, ALL_ONESV) ;		\
    MAX_A = EC::vec_blend(MAX_A, AD, __MASK) ;				\
    BEST_A_DIAGV = EC::vec_blend(BEST_A_DIAGV, DIAGV, __MASK) ;		\
  }

#define UPDATE_BEST_B_WITH_END_BONUS() {				\
    typename EC::Vec __BD_EB = EC::vec_add(BD, endBonusV) ;		\
    typename EC::Vec __MASK = EC::vec_compare_gt(__BD_EB, MAX_B) ;	\
    MAX_B = EC::vec_blend(MAX_B, __BD_EB, __MASK) ;			\
    BEST_B_DIAGV = EC::vec_blend(BEST_B_DIAGV, DIAGV, __MASK) ;		\
  }

#define UPDATE_BEST_A_WITH_END_BONUS() {				\
    typename EC::Vec __AD_EB = EC::vec_add(AD, endBonusV) ;		\
    typename EC::Vec __MASK = EC::vec_compare_gt(__AD_EB, MAX_A) ;	\
    typename EC::Vec __MASK2 = EC::vec_compare_gt(ZEROV, __AD_EB) ;	\
    terminate = EC::vec_test_all_ones(__MASK2, ALL_ONESV) ;		\
    MAX_A = EC::vec_blend(MAX_A, __AD_EB, __MASK) ;			\
    BEST_A_DIAGV = EC::vec_blend(BEST_A_DIAGV, DIAGV, __MASK) ;		\
  }


#define PROCESS_B() {							\
  									\
    /* Shift the ref vec before computing B vectors*/			\
    RV = EC::vec_shift_left_and_insert(RV, nextRV) ;			\
    nextRV = EC::vec_shift_left(nextRV) ;				\
									\
    /* Compute base score vectors (based on match/mismatch/ambig) */	\
    COMPUTE_SCOREV() ;							\
									\
    /* Compute B vecs */						\
    COMPUTE_B() ;							\
									\
    /* Update max scores and best indices corresponding to B vector */	\
    UPDATE_BEST_B() ;							\
  }
  
  // Repeat the corresponding operations for A
#define PROCESS_A() {					\
    QV = EC::vec_shift_right_and_insert(QV, nextQV) ;	\
    nextQV = EC::vec_shift_right(nextQV) ;		\
    COMPUTE_SCOREV() ;					\
    COMPUTE_A() ;					\
    UPDATE_BEST_A() ;					\
  }

#define PERFORM_SINGLE_ITER() {			\
    PROCESS_B() ;				\
    PROCESS_A() ;				\
    /* Increment diag vector */			\
    DIAGV = EC::vec_add(DIAGV, ONEV) ;		\
  }


#define PROCESS_B_WITH_END_BONUS() {			\
    RV = EC::vec_shift_left_and_insert(RV, nextRV) ;	\
    nextRV = EC::vec_shift_left(nextRV) ;		\
    COMPUTE_SCOREV() ;					\
    COMPUTE_B() ;					\
    UPDATE_BEST_B_WITH_END_BONUS() ;			\
  }

#define PROCESS_A_WITH_END_BONUS() {			\
    QV = EC::vec_shift_right_and_insert(QV, nextQV) ;	\
    nextQV = EC::vec_shift_right(nextQV) ;		\
    COMPUTE_SCOREV() ;					\
    COMPUTE_A() ;					\
    UPDATE_BEST_A_WITH_END_BONUS() ;			\
  }

#define PERFORM_SINGLE_ITER_WITH_END_BONUS() {		\
    PROCESS_B_WITH_END_BONUS() ;			\
    endBonusV = EC::vec_shift_right(endBonusV) ;	\
    PROCESS_A_WITH_END_BONUS() ;			\
    /* Increment diag vector */				\
    DIAGV = EC::vec_add(DIAGV, ONEV) ;			\
  }

#define PRINT_VECS_FOR_DEBUG() {		\
    cout << "DiagV: " << endl ;			\
    EC::vec_print(DIAGV) ;			\
    cout << "RV: " << endl ;			\
    EC::vec_print(RV) ;				\
    cout << "QV: " << endl ;			\
    EC::vec_print(QV) ;				\
    cout << "SCOREV: " << endl ;		\
    EC::vec_print(SCOREV) ;			\
    cout << "EndBonus: " << endl ;		\
    EC::vec_print(endBonusV) ;			\
    cout << "AD: " << endl ;			\
    EC::vec_print(AD) ;				\
    cout << "AH: " << endl ;			\
    EC::vec_print(AH) ;				\
    cout << "AV: " << endl ;			\
    EC::vec_print(AV) ;				\
    cout << "BD: " << endl ;			\
    EC::vec_print(BD) ;				\
    cout << "BH: " << endl ;			\
    EC::vec_print(BH) ;				\
    cout << "BV: " << endl ;			\
    EC::vec_print(BV) ;				\
    cout << "MAX_A: " << endl ;			\
    EC::vec_print(MAX_A) ;			\
    cout << "BEST_A: " << endl ;		\
    EC::vec_print(BEST_A_DIAGV) ;		\
    cout << "MAX_B: " << endl ;			\
    EC::vec_print(MAX_B) ;			\
    cout << "BEST_B: " << endl ;		\
    EC::vec_print(BEST_B_DIAGV) ;		\
  }


template<class EC>
int extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore, int *_max_off)
{

  int diagLen = w + 1 ;
  if (diagLen > EC::MAX_VEC_LEN) {
    cout << "Error: This implementation cannot handle w = " << w << endl ;
    cout << "The max w allowed = " << EC::MAX_VEC_LEN - 1 << endl ;
    exit(0) ;
  }

  assert(m == 5) ; // 5 bases, including ambig base
  int matchScore = mat[0] ;
  int mismatchScore = mat[1] ;
  int ambigScore = mat[AMBIG_BASE] ;

#ifdef DEBUG
  // Double check that the matrix is consistent
  for (int i=0; i <= AMBIG_BASE; ++i) 
    for (int j=0; j <= AMBIG_BASE; ++j) {
      if (i == AMBIG_BASE || j == AMBIG_BASE)
	assert(mat[i*m+j] == ambigScore) ;
      else if (i == j)
	assert(mat[i*m+j] == matchScore) ;
      else
	assert(mat[i*m+j] == mismatchScore) ;
    }
#endif

  int midElt = w / 2 ; // rounded down

  ExtendProblem prob(query, qlen, target, tlen, w, h0) ;

  typename EC::Vec AH, AV, AD, BH, BV, BD ; // the diagonals
  typename EC::Vec QV, RV ; // query and ref vectors
  typename EC::Vec nextQV=EC::vec_set1(0), nextRV=EC::vec_set1(0) ; // the vectors holding the next query/ref chars
  typename EC::Vec MAX_A ; // to keep track of the max score at each lane
  typename EC::Vec MAX_B ; // to keep track of the max score at each lane
  typename EC::Vec BEST_A_DIAGV ; // to keep track of the best diag index at each lane
  typename EC::Vec BEST_B_DIAGV ; // to keep track of the best diag index at each lane
  typename EC::Vec DIAGV ; // each lane has the same diagIndex value
  typename EC::Vec endBonusV = EC::vec_set1(0) ;

  // constant vectors
  static const typename EC::Vec ZEROV = EC::vec_set1(0) ;
  static const typename EC::Vec ALL_ONESV = EC::vec_set1(-1) ;
  static const typename EC::Vec ONEV = EC::vec_set1(1) ;
  static const typename EC::Vec BADV = EC::vec_set1(EC::BAD_SCORE) ;
  static const typename EC::Vec DESV = EC::vec_set1(-e_del) ; // delete-extend
  static const typename EC::Vec DOSV = EC::vec_set1(-o_del-e_del) ; // delete-open+extend
  static const typename EC::Vec IESV = EC::vec_set1(-e_ins) ; // insert-extend
  static const typename EC::Vec IOSV = EC::vec_set1(-o_ins-e_ins) ; // insert-open+extend
  static const typename EC::Vec MATCHSV = EC::vec_set1(matchScore) ;
  static const typename EC::Vec MISMATCHSV = EC::vec_set1(mismatchScore) ;
	static const typename EC::Vec INVALIDSV = EC::vec_set1(-o_ins-o_del-e_ins-e_del) ;
  static const typename EC::Vec AMBIGSV = EC::vec_set1(ambigScore) ;
  static const typename EC::Vec AMBIGBV = EC::vec_set1(AMBIG_BASE) ;
  static const typename EC::Vec INVALIDRBV = EC::vec_set1(INVALID_BASE_R) ;
  static const typename EC::Vec INVALIDQBV = EC::vec_set1(INVALID_BASE_Q) ;

  typename EC::Word arr[EC::MAX_VEC_LEN] ;

  // initialize the Q vector
  for (int vi=0; vi <= midElt; ++vi) {
    arr[vi] = INVALID_BASE_Q ;
  }

  int qi = 0 ;
  for (int vi=midElt+1; vi < EC::MAX_VEC_LEN; ++vi) {
    arr[vi] = prob.getQ(qi++) ;
  }
  QV = EC::vec_set(arr) ;

  // initialize the RV vector
  for (int vi=EC::MAX_VEC_LEN-1; vi >= midElt; --vi) {
    arr[vi] = INVALID_BASE_R ;
  }
  int ri = 0 ;
  for (int vi=midElt-1; vi >=0; --vi) {
    arr[vi] = prob.getR(ri++) ;
  }
  RV = EC::vec_set(arr) ;

  // Initialize the diagonal vectors
  AH = BADV ;
  AV = BADV ;
  BH = BADV ;
  BV = BADV ;
  BD = BADV ;

  bool terminate = false ;

  for (int vi=0; vi < EC::MAX_VEC_LEN; ++vi) {
    arr[vi] = EC::BAD_SCORE ;
  }
  arr[midElt] = h0 ; // initial score
  AD = EC::vec_set(arr) ;

  // Initialize the max scores and best indices
  MAX_A = AD ;
  BEST_A_DIAGV = ZEROV ;
  MAX_B = BD ;
  BEST_B_DIAGV = ZEROV ;

  DIAGV = ONEV ;

  typename EC::Vec SCOREV ; /*temp vec to compute the match score*/  


  // We define 3 phases as follows:
  // Phase 1: When the last elt of vec A is above the bottom row
  //          # of iters = qlen - (diagLen-midElt-1) - 1 (may be negative)
  // Phase 2: When the last elt of vec A is at the bottom row
  //          # of iters = 1 (typical) or 0 (in case the diagLen is too large compared to qlen)
  //                     = 1 if (# of iters in phase1) >= 0 ;
  //                       0 otherwise
  // Phase 3: When the last elt of vec A is below the bottom row
  //          # of iters = min(diagLen, floor((tlen-qlen)/2)+diagLen-midElt 
  //                              The first term is selected when we have a long enough reference.
  //                              The second term is selected when the query and ref lengths are close.
  //                       + min(0, # of iters in phase1 +1),
  //                              If the # of iters in phase1 is < -1, we adjust the iter cnt.
  //                              (In a typical case, we assume phase 3 starts when the last elt
  //                               of vec A is 1 row below the last row. If the # of iters in phase1
  //                               is < -1, we have to adjust this accordingly.)
  // 
  // Another upper bound for # of iterations is for cases where tlen < qlen:
  //       # iters <= tlen + (diagLen - midElt)
  //       Note: The equality above happens when the last elt of vec A is to the right of last column
  // 
  //
  // While running these 3 phases, we need to fill in nextQV and nextRV vectors with EC::MAX_VEC_LEN
  // elements. So, we need to divide these phases further to refill nextQV and nextRV occasionally.

  int numItersUpperBound = tlen + diagLen - midElt ; // upper bound due to target length

  // Phase 1: 
  int numItersPhase1 = std::min(numItersUpperBound, qlen - (diagLen - midElt)) ;
  int numItersPhase3 = std::min(diagLen,(tlen-qlen)/2 + diagLen - midElt) + std::min(0, numItersPhase1+1) ;
  if (numItersPhase1+1 >= numItersUpperBound)
    numItersPhase3 = 0 ;
  

  int numRemainingEltsInNextQR = 0 ;

  if (numItersPhase1 > 0) {
    int numFullIters = numItersPhase1 / EC::MAX_VEC_LEN ;
    int numPartialIters = numItersPhase1 % EC::MAX_VEC_LEN ;
    

    for (int fullIterIndex=0; fullIterIndex < numFullIters; ++fullIterIndex) {
      POPULATE_NEXTQV() ;
      POPULATE_NEXTRV() ;

      // bvi: base vector index (indexing into nextRV and nextQV)
      for (int bvi=0; bvi < EC::MAX_VEC_LEN; ++bvi) {
	PERFORM_SINGLE_ITER() ;
#ifdef DEBUG
	PRINT_VECS_FOR_DEBUG() ;
#endif
	if (terminate) goto LoopEnd ;
      }
    }
    POPULATE_NEXTQV() ;
    POPULATE_NEXTRV() ; 

    for (int bvi=0; bvi < numPartialIters; ++bvi) {
      PERFORM_SINGLE_ITER() ;
#ifdef DEBUG
      PRINT_VECS_FOR_DEBUG() ;
#endif
      if (terminate) goto LoopEnd ;
    }

    numRemainingEltsInNextQR = EC::MAX_VEC_LEN - numPartialIters ;
  }

  // In phase 2 and 3, we should add the end_bonus to each score at the last row.
  // In phase 2, the last elt of vec A is at the last row.
  // Initialize the end bonus vec:
  for (int vi=0; vi < EC::MAX_VEC_LEN; ++vi)
    arr[vi] = 0 ;
  arr[diagLen-1] = end_bonus ;
  endBonusV = EC::vec_set(arr) ;


  // Phase 2: 
  if (numItersPhase1 >= 0 && numItersPhase1 < numItersUpperBound) {
    assert(numRemainingEltsInNextQR >= 0) ;
    if (numRemainingEltsInNextQR == 0) {
      POPULATE_NEXTQV() ;
      POPULATE_NEXTRV() ;
      numRemainingEltsInNextQR = EC::MAX_VEC_LEN ;
    }

    PROCESS_B() ;
    PROCESS_A_WITH_END_BONUS() ;    
    DIAGV = EC::vec_add(DIAGV, ONEV) ;

#ifdef DEBUG
    PRINT_VECS_FOR_DEBUG() ;
#endif
    --numRemainingEltsInNextQR ;

    if (terminate) goto LoopEnd ;
  }

  // Phase 3:
  assert(numItersPhase3 <= EC::MAX_VEC_LEN) ;
  
  if (numItersPhase3 > 0) {
    
    if (numRemainingEltsInNextQR == 0) {
      POPULATE_NEXTQV_VAR(numItersPhase3) ;
      POPULATE_NEXTRV_VAR(numItersPhase3) ;
      numRemainingEltsInNextQR = numItersPhase3 ;
    }

    // We may need to refill nextQV and nextRV vecs in the middle. So, divide phase3 into 2 stages.
    int numItersPhase3a = std::min(numRemainingEltsInNextQR, numItersPhase3) ;
    int numItersPhase3b = std::max(0, numItersPhase3 - numRemainingEltsInNextQR) ;

    for (int iter=0; iter<numItersPhase3a; ++iter) {
      PERFORM_SINGLE_ITER_WITH_END_BONUS() ;
#ifdef DEBUG
      PRINT_VECS_FOR_DEBUG() ;
#endif
      if (terminate) goto LoopEnd ;
    }

    POPULATE_NEXTQV_VAR(numItersPhase3b) ;
    POPULATE_NEXTRV_VAR(numItersPhase3b) ;

    for (int iter=0; iter<numItersPhase3b; ++iter) {
      PERFORM_SINGLE_ITER_WITH_END_BONUS() ;
#ifdef DEBUG
      PRINT_VECS_FOR_DEBUG() ;
#endif
      if (terminate) goto LoopEnd ;
    }
  }
  
 LoopEnd:

  // Compute the best score
  typename EC::Word* arrMaxB = ((typename EC::Word*) &MAX_B) ;
  typename EC::Word* arrMaxA = ((typename EC::Word*) &MAX_A) ;
  typename EC::Word* arrBestBDiag = ((typename EC::Word*) &BEST_B_DIAGV) ;
  typename EC::Word* arrBestADiag = ((typename EC::Word*) &BEST_A_DIAGV) ;

  typename EC::Word maxLocalScore = h0 ;
  typename EC::Word bestRIndex = -1 ;
  typename EC::Word bestQIndex = -1 ;

  for (int lane=0; lane < diagLen; ++lane) {

    if (arrMaxA[lane] >= maxLocalScore) {
#ifdef DEBUG
      cout << lane << ": " << maxLocalScore << " < " << arrMaxA[lane] 
      	   << ", d = " << arrBestADiag[lane] << endl ;
#endif
      int d = arrBestADiag[lane] ; // the diag index for this score
      int currRIndex = d - 1 + midElt - lane ;
      int currQIndex = d - 1 - midElt + lane ;

      // If arrMaxA[lane] == maxLocalScore, then we favor smaller query
      // indices to be consistent with the BWA implementation.
      if (arrMaxA[lane] > maxLocalScore || currQIndex < bestQIndex) {
	bestRIndex = currRIndex ;
	bestQIndex = currQIndex ;
	maxLocalScore = arrMaxA[lane] ;
      }

#ifdef DEBUG
      cout << "best indices: " << bestRIndex << " " << bestQIndex << endl ;
#endif
    }

    if (arrMaxB[lane] >= maxLocalScore) {

      int d = arrBestBDiag[lane] ; // the diag index for this score
      int currRIndex = d - 1 + midElt - lane ;
      int currQIndex = d - 2 - midElt + lane ;

      // If arrMaxB[lane] == maxLocalScore, then we favor smaller query
      // indices to be consistent with the BWA implementation.
      if (arrMaxB[lane] > maxLocalScore || currQIndex < bestQIndex) {
	bestRIndex = currRIndex ;
	bestQIndex = currQIndex ;
	maxLocalScore = arrMaxB[lane] ;
      }
    }

  }

  // Convert it to the form ksw_extend uses:
  int gscore = 0, lscore = 0, tle = 0, gtle = 0 ;
  if (bestQIndex+1 == qlen) { // global score is chosen over local score
    gscore = maxLocalScore - end_bonus ; // end_bonus was added to maxLocalScore before
    gtle = std::min(tlen, bestRIndex + 1) ;
    // Because of band limitations, it's possible that max score ends up beyond ref index.
    // This can happen when tlen < queryLen and the band that overlaps with the last query
    // row is beyond column tlen.
  }
  else {
    lscore = maxLocalScore ;
    tle = std::min(tlen, bestRIndex + 1) ;
    // Because of band limitations, it's possible that max score ends up beyond ref index.
    // The min operation probably is not needed in case of local alignment, but keeping it
    // here for safety.
  }

#ifdef DEBUG
  cout << "Results: " << gscore << " " <<  bestQIndex+1 << " " << gtle << " " << lscore << " " << tle << endl ;
#endif

  if (_qle) *_qle = bestQIndex+1 ;
  if (_tle) *_tle = tle ;
  if (_gtle) *_gtle = gtle ; 
  if (_gscore) *_gscore = gscore ; 

	assert(*_qle <= qlen) ;
	assert(*_tle <= tlen) ;

  return lscore ;
}


// wrappers

int ksw_extend_sse_8(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off) {

  return extend<EC16BitSSE>(qlen, query, tlen, target, m, mat, gapo, gape, gapo,
			    gape, w, end_bonus, zdrop, h0, qle, tle, gtle,
			    gscore, max_off) ;
}

int ksw_extend_sse_16(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off) {

  return extend<EC16BitSSE_16>(qlen, query, tlen, target, m, mat, gapo, gape, gapo,
			       gape, w, end_bonus, zdrop, h0, qle, tle, gtle,
			       gscore, max_off) ;
}


#ifdef USE_AVX2
int ksw_extend_avx2_16(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off) {

  return extend<EC16BitAvx2>(qlen, query, tlen, target, m, mat, gapo, gape, gapo,
			     gape, w, end_bonus, zdrop, h0, qle, tle, gtle,
			     gscore, max_off) ;
}


int ksw_extend_avx2_32(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int end_bonus, int zdrop, int h0, int *qle, int *tle, int *gtle, int *gscore, int *max_off) {

  return extend<EC16BitAvx2_32>(qlen, query, tlen, target, m, mat, gapo, gape, gapo,
				 gape, w, end_bonus, zdrop, h0, qle, tle, gtle,
				 gscore, max_off) ;
}

#endif // #ifdef USE_AVX2




#ifdef EXT_AVX2_MAIN

const int costMatrixRowCnt = 5 ;
const int8_t costMatrix[] = {
  1, -4, -4, -4, -1,
  -4, 1, -4, -4, -1,
  -4, -4, 1, -4, -1,
  -4, -4, -4, 1, -1,
  -1, -1, -1, -1, -1 
} ;

const int gapo = 6 ;
const int gape = 1 ; 

bool test1() {

  const uint8_t query[] = {2, 0, 3, 1, 2} ;
  const uint8_t ref[] = {2, 0, 3, 1, 2} ;
  

  int qlen = 5 ;
  int rlen = 5 ;
  int zdrop = 0 ;
  int h0 = 0 ;
  int w = 1 ;
  int endBonus = 0 ;
  

  int qle, tle, gtle, gscore ;

  int localScore = 
    extend<EC16BitAvx2> (qlen, query, rlen, ref, costMatrixRowCnt, costMatrix,
			 gapo, gape, gapo, gape, w, endBonus, zdrop, h0,
			 &qle, &tle, &gtle, &gscore, NULL) ;

  cout << "Local score computed = " << localScore << endl ;
  cout << "Query and ref alignment= " << qle << " and " << tle << endl ;

  if (qle != qlen || tle != qlen) {
    cout << "Error: Test 1 failed" << endl ;
  }

}

bool test2() {

  const uint8_t query[] = {2, 0, 3, 1, 2, 3, 1, 2, 0} ;
  const uint8_t ref[] = {2, 4, 2, 0, 3, 1, 2, 3, 1, 2, 0} ;
  

  int qlen = 9 ;
  int rlen = 11 ;
  int zdrop = 0 ;
  int h0 = 0 ;
  int w = 15 ;
  int endBonus = 0 ;

  int qle, tle, gtle, gscore ;

  int localScore = 
    extend<EC16BitAvx2> (qlen, query, rlen, ref, costMatrixRowCnt, costMatrix,
			 gapo, gape, gapo, gape, w, endBonus, zdrop, h0,
			 &qle, &tle, &gtle, &gscore, NULL) ;

  cout << "Local score computed = " << localScore << endl ;
  cout << "Query and ref alignment= " << qle << " and " << tle << endl ;

  if (qle != qlen || tle != rlen) {
    cout << "Error: Test 1 failed" << endl ;
  }

}



int main() {

  test2() ;

  return 0 ;
}

#endif // #ifdef EXT_AVX2_MAIN
