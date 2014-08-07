/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/



#ifndef _ED_INTRAVED_H
#define _ED_INTRAVED_H


#include <immintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>

#ifdef USE_SSE4
#include <tmmintrin.h>
#include <smmintrin.h>

#define BLENDV(__v1, __v2, __pred)		\
  _mm_blendv_epi8(__v1, __v2, __pred)

#else // #ifdef USE_SSE4

#define BLENDV(__v1, __v2, __pred)					\
  _mm_or_si128(_mm_andnot_si128(__pred, __v1), _mm_and_si128(__v2, __pred))

#endif // #ifdef USE_SSE4

#include "ed_intrav64.h"
#include "ed_intrav64x2.h"

		    /*
inline void printVec8(const __m128i& v) {
  
  cout << "{" ;
  cout << " " << _mm_extract_epi8(v, 0) ;
  cout << " " <<_mm_extract_epi8(v, 1) ;
  cout << " " <<_mm_extract_epi8(v, 2) ;
  cout << " " << _mm_extract_epi8(v, 3) ;
  cout << " " <<_mm_extract_epi8(v, 4) ;
  cout << " " <<_mm_extract_epi8(v, 5) ;
  cout << " " <<_mm_extract_epi8(v, 6) ;
  cout << " " <<_mm_extract_epi8(v, 7) ;
  cout << " " <<_mm_extract_epi8(v, 8) ;
  cout << " " <<_mm_extract_epi8(v, 9) ;
  cout << " " <<_mm_extract_epi8(v, 10) ;
  cout << " " <<_mm_extract_epi8(v, 11) ;
  cout << " " <<_mm_extract_epi8(v, 12) ;
  cout << " " <<_mm_extract_epi8(v, 13) ;
  cout << " " <<_mm_extract_epi8(v, 14) ;
  cout << " " <<_mm_extract_epi8(v, 15) ;

  cout << "}" << endl ;
}
		    */
inline void printVec16(const __m128i& v) {
  
  cout << "{" ;
  cout << " " << _mm_extract_epi16(v, 0) ;
  cout << " " << _mm_extract_epi16(v, 1) ;
  cout << " " << _mm_extract_epi16(v, 2) ;
  cout << " " << _mm_extract_epi16(v, 3) ;
  cout << " " <<_mm_extract_epi16(v, 4) ;
  cout << " " <<_mm_extract_epi16(v, 5) ;
  cout << " " <<_mm_extract_epi16(v, 6) ;
  cout << " " <<_mm_extract_epi16(v, 7) ;

  cout << "}" << endl ;
}


class EDVec128Every16 {
  __m128i vec_ ;
  static const int WORD_CNT_ = 8 ;
  static const int PERIOD_ = 16 ;

  // Update equality operator if more members are added

 public:
  EDVec128Every16() {}

  explicit EDVec128Every16(__m128i v):vec_(v) {}
  explicit EDVec128Every16(int16_t val):vec_(_mm_set1_epi16(val)) {}
  explicit EDVec128Every16(const BitVec64& bv): vec_(_mm_set_epi64x(0, bv.bitV)) {}
  explicit EDVec128Every16(const BitVec64x2& bv): vec_(_mm_set_epi64x(bv.bitV[1], bv.bitV[0])) {}

  inline void setAll(int16_t val) {
    vec_ = _mm_set1_epi16(val) ;    
  }

  static int16_t getMaxWordVal() { return numeric_limits<int16_t>::max() ;}

  inline int getWord(int index) const {
    switch(index) {
    case 0: return _mm_extract_epi16(vec_, 0) ;
    case 1: return  _mm_extract_epi16(vec_, 1) ;
    case 2: return  _mm_extract_epi16(vec_, 2) ;
    case 3: return  _mm_extract_epi16(vec_, 3) ;
    case 4: return  _mm_extract_epi16(vec_, 4) ;
    case 5: return  _mm_extract_epi16(vec_, 5) ;
    case 6: return  _mm_extract_epi16(vec_, 6) ;
    case 7: return  _mm_extract_epi16(vec_, 7) ;
    default: assert(false); return 0 ;
    }
  }

  bool allLessThanOrEqualTo(const EDVec128Every16& other) const {
    // each elt in pred will be set to 0xFFFF if the element in this vector 
    // is greater than the corresponding element in other
    __m128i pred = _mm_cmpgt_epi16(vec_, other.vec_) ; 
    int flags = _mm_movemask_epi8(pred) ; // creates a mask from MSB of each 8-bit elt in pred

    return flags == 0 ; // if there's at least one elt greater than "other", flags != 0
  }

  inline void setWords (int* v) {
    vec_ = _mm_set_epi16(v[7], v[6], v[5], v[4], v[3], v[2], v[1], v[0]) ;
  }

  inline void setFirstWord(int val) {
    vec_ = _mm_insert_epi16(vec_, val, 0) ;
  }

  inline void setWordsAsMask () {
    vec_ = _mm_set1_epi16(0x1) ;
  }

  static int getLastWordIndexFor (int queryLen) {
    int baseIndex = queryLen - 1 ;
    return (baseIndex >= 63) ? 
      (((baseIndex - 63) / PERIOD_) + WORD_CNT_/2) :
      (baseIndex/PERIOD_) ;
  }

  // If the query length is greater than 63, we need two 64-words.
  // In that case, we will have:
  // bitOffset(lowerWord probes) = bitOffset(upper word probes)-1
  // Consider the following example:
  // query len = 79, probeOffset = (79-64) % 16 = 15. The bit offset of the 
  // probes in the upper 64-bit word will be 15, while the bit offsets of the probes 
  // in the lower 64-bit word will be 14 in this case. So, the probes will be at
  // indices: 78 (in upper word), 62, 46, 30, 14 (in lower word). Note that if the
  // bit offset were also 15 in the lower word, one probe would correspond to position
  // 63, which has no corresponding base (because bit 63 stores the carryout bit of the
  // lower word). We can generalize this notion for all query lengths. 
  // As a potential corner case, consider query length of 64. In this case, the 
  // probes in the upper 64-bit word will have bitOffset = 0, and the lower
  // ones will have bitOffset = -1 (conceptually). The probe locations will be at:
  // 63 (upper word), 47, 31, 15, -1 (lower word). The corresponding lengths of
  // the subsequences will be: 64, 48, 32, 16, 0. Hence the last probe will be ignored.
  //
  // In the following function, the probe offset corresponding to the upper word will be
  // returned. Internally, we'll use probe offset one less for the lower word.
  // If the query length is less than 63, we will have only the lower word, and the
  // last probe location should match the MSB. So, the return value of the following 
  // function should be one more than the real probe bit offset.
  // Consider an example where query length is 60. The probe locations should be at 
  // bit offsets: 59, 43, 27, 11. The return value of the following function in this 
  // case will be 12, but internally, bit offset 11 will be used downstream.
  static int getProbeOffsetFor (int queryLen) {
    return queryLen % PERIOD_ ;
    // It is not (queryLen-1) % PERIOD, because:
    // For queries longer than 63, one bit in the bit vector is unused (the one at position 63).
    // For queries shorter than or equal to 63, the value returned is one more than the real
    // offset (see the comments above).
  }

  void setWordsAsEndBonus(int queryLen, int endBonus) {
    int v[WORD_CNT_] = {0, 0, 0, 0, 0, 0, 0, 0} ;
    
    int lastWordIndex = getLastWordIndexFor(queryLen) ;    
    v[lastWordIndex] = endBonus ;

    setWords(v) ;
  }

  void setWordsAsBadScore(int queryLen, int thr, int infScore) {
    int v[WORD_CNT_] = {thr, thr, thr, thr, thr, thr, thr, thr} ;
    
    int lastWordIndex = getLastWordIndexFor(queryLen) ;    
    for (int i=lastWordIndex+1; i < WORD_CNT_; ++i) {
      v[lastWordIndex] = infScore ;
    }

    setWords(v) ;
  }

  static int getDistAt(int wordIndex, int queryLen) {
    return wordIndex * PERIOD_ + getProbeOffsetFor(queryLen) ;
  }

  void setWordsAsDist (int queryLen, int distWt) {
    //int16_t v[WORD_CNT_] = {0, 0, 0, 0, 0, 0, 0, 0} ;
    int v[WORD_CNT_] ;
    int bitOffset = getProbeOffsetFor(queryLen) ;

    for (int wi=0; wi < WORD_CNT_; ++wi) {
      v[wi] = distWt * (wi * PERIOD_ + bitOffset) ;
    }
    
    /*
    for (int wi=0; wi < WORD_CNT_/2 - 1; ++wi) {
      v[wi] = distWt * (wi * PERIOD_ + bitOffset + 1) ;
    }

    for (int wi=WORD_CNT_/2-1; wi < WORD_CNT_; ++wi) {
      v[wi] = distWt * (wi * PERIOD_ + bitOffset) ; // ... + PERIOD_-1 + bitOffset+1
    }
    */

    setWords(v) ;
  }

  void setMin(const EDVec128Every16& other) {
    vec_ = _mm_min_epi16(vec_, other.vec_) ;
  }

  void setMin(const EDVec128Every16& other, EDVec128Every16& bestIndices, const EDVec128Every16& otherIndices) {
    __m128i pred = _mm_cmplt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other < this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
  }

  void setMax(const EDVec128Every16& other, EDVec128Every16& bestIndices, const EDVec128Every16& otherIndices) {
    __m128i pred = _mm_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
  }

  void setMax(const EDVec128Every16& other, EDVec128Every16& bestIndices, const EDVec128Every16& otherIndices, EDVec128Every16& bestAccumDist, const EDVec128Every16& otherAccumDist) {
    __m128i pred = _mm_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
    bestAccumDist.vec_ = BLENDV(bestAccumDist.vec_, otherAccumDist.vec_, pred) ;
  }

  // The return value will be nonzero iff at least one of the entries is updated
  uint16_t setMaxAndReturnFlag(const EDVec128Every16& other, EDVec128Every16& bestIndices, const EDVec128Every16& otherIndices) {
    __m128i pred = _mm_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
    
    return (uint16_t) _mm_movemask_epi8(pred) ;
  }

  // The return value will be nonzero iff at least one of the entries is updated
  uint16_t setMaxAndReturnFlag(const EDVec128Every16& other, EDVec128Every16& bestIndices, const EDVec128Every16& otherIndices, EDVec128Every16& bestAccumDist, const EDVec128Every16& otherAccumDist) {
    __m128i pred = _mm_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
    bestAccumDist.vec_ = BLENDV(bestAccumDist.vec_, otherAccumDist.vec_, pred) ;
    
    return (uint16_t) _mm_movemask_epi8(pred) ;
  }

  void addThirdIfFirstGTSecond(const EDVec128Every16& first, const EDVec128Every16& second, 
			       const EDVec128Every16& third, const EDVec128Every16& zero) {

    __m128i pred = _mm_cmpgt_epi16(first.vec_, second.vec_) ; 
    // pred is 0xFFFF if first > second

    // If pred is 0, then blendv will add zero; o/w it will add third
    vec_ = _mm_add_epi16(vec_, BLENDV(zero.vec_, third.vec_, pred)) ;
  }

  /*
  EDVec128Every16 abs(const EDVec128Every16& other) const {
    return EDVec128Every16(_mm_abs_epi16(_mm_sub_epi16(this->vec_, other.vec_))) ;
  }

  EDVec128Every16 abs() const {
    return EDVec128Every16(_mm_abs_epi16(this->vec_)) ;
  }
  */

  EDVec128Every16 subSat(const EDVec128Every16& other) const {
    return EDVec128Every16(_mm_subs_epu16(this->vec_, other.vec_)) ;
  }

  inline EDVec128Every16 shiftBitsRightWithinWords(int shiftVal) {
    return EDVec128Every16(_mm_srli_epi16(vec_, shiftVal)) ;
  }

  inline void shiftWordsLeftByOne() {
    vec_ = _mm_slli_si128(vec_, 2) ; // 16-bit shift = 2 * 8-bit shift
  }

  inline bool operator == (const EDVec128Every16& other) {
    return _mm_movemask_epi8(_mm_cmpeq_epi16(vec_, other.vec_)) == 0xFFFF ;
  }

  inline EDVec128Every16& operator += (const EDVec128Every16& other) {
    vec_ = _mm_add_epi16(vec_, other.vec_) ;
    return *this ;
  }

  inline EDVec128Every16& operator -= (const EDVec128Every16& other) {
    vec_ = _mm_sub_epi16(vec_, other.vec_) ;
    return *this ;
  }

  inline EDVec128Every16 operator+ (const EDVec128Every16& other) const {
    return EDVec128Every16(_mm_add_epi16(vec_, other.vec_)) ;
  }

  inline EDVec128Every16 operator- (const EDVec128Every16& other) const {
    return EDVec128Every16(_mm_sub_epi16(vec_, other.vec_)) ;
  }

  inline EDVec128Every16 operator* (const EDVec128Every16& other) const {
    return EDVec128Every16(_mm_mullo_epi16(vec_, other.vec_)) ;
  }

  inline EDVec128Every16 operator& (const EDVec128Every16& other) const {
    return EDVec128Every16(_mm_and_si128(vec_, other.vec_)) ;
  }

  static int WORD_CNT() {
    return WORD_CNT_ ;
  }

  static int PERIOD() {
    return PERIOD_ ;
  }

  friend std::ostream& operator<< (std::ostream& os, const EDVec128Every16& v) {
    os << "{" ;
    for (int i=0; i < WORD_CNT_; ++i) {
      os << v.getWord(i) << " " ;
    }
    os << "}" ;
    return os ;
  }

} ;


//
// Edit distance vector for a 64-bit query.
// 16-bit words will be stored in an SSE vector in interleaved order as follows:
//           w7 w5 w3 w1 w6 w4 w2 w0
//
class EDVec64Every8 {
  __m128i vec_ ;
  __m128i shiftLeftIndices_ ; // can be made static later
  static const int WORD_CNT_ = 8 ;
  static const int PERIOD_ = 8 ;

  void init_() {
    shiftLeftIndices_ = _mm_set_epi8(7, 6, 5, 4, 3, 2, 1, 0, 13, 12, 11, 10, 9, 8, 0xFF, 0xFF) ;
  }


 public:
  EDVec64Every8() {}

  explicit EDVec64Every8(__m128i v):vec_(v) {init_() ;}
  explicit EDVec64Every8(int16_t val):vec_(_mm_set1_epi16(val)) {init_() ;}
  explicit EDVec64Every8(const BitVec64& bv): vec_(_mm_set_epi64x(bv.bitV>>8, bv.bitV)) {init_() ;}


  inline void setAll(int16_t val) {
    vec_ = _mm_set1_epi16(val) ;    
  }

  static int16_t getMaxWordVal() { return numeric_limits<int16_t>::max() ;}

  inline int getWord(int index) const {
    switch(index) {
    case 0: return _mm_extract_epi16(vec_, 0) ;
    case 1: return  _mm_extract_epi16(vec_, 4) ;
    case 2: return  _mm_extract_epi16(vec_, 1) ;
    case 3: return  _mm_extract_epi16(vec_, 5) ;
    case 4: return  _mm_extract_epi16(vec_, 2) ;
    case 5: return  _mm_extract_epi16(vec_, 6) ;
    case 6: return  _mm_extract_epi16(vec_, 3) ;
    case 7: return  _mm_extract_epi16(vec_, 7) ;
    default: assert(false); return 0 ;
    }
  }

  bool allLessThanOrEqualTo(const EDVec64Every8& other) const {
    // each elt in pred will be set to 0xFFFF if the element in this vector 
    // is greater than the corresponding element in other
    __m128i pred = _mm_cmpgt_epi16(vec_, other.vec_) ; 
    int flags = _mm_movemask_epi8(pred) ; // creates a mask from MSB of each 8-bit elt in pred

    return flags == 0 ; // if there's at least one elt greater than "other", flags != 0
  }


  inline void setWords (int* v) {
    vec_ = _mm_set_epi16(v[7], v[5], v[3], v[1], v[6], v[4], v[2], v[0]) ;
  }

  inline void setFirstWord(int val) {
    vec_ = _mm_insert_epi16(vec_, val, 0) ;
  }

  inline void setWordsAsMask () {
    vec_ = _mm_set1_epi16(0x1) ;
  }

  static int getLastWordIndexFor (int queryLen) {
    return (queryLen-1) / PERIOD_ ; // logical word index - no interleaving at this level
  }

  static int getProbeOffsetFor (int queryLen) {
    return (queryLen-1) % PERIOD_ ;
  }

  void setWordsAsEndBonus(int queryLen, int endBonus) {
    int v[WORD_CNT_] = {0, 0, 0, 0, 0, 0, 0, 0} ;
    
    int lastWordIndex = getLastWordIndexFor(queryLen) ;    
    v[lastWordIndex] = endBonus ;

    setWords(v) ;
  }

  void setWordsAsBadScore(int queryLen, int thr, int infScore) {
    int v[WORD_CNT_] = {thr, thr, thr, thr, thr, thr, thr, thr} ;
    
    int lastWordIndex = getLastWordIndexFor(queryLen) ;    
    for (int i=lastWordIndex+1; i < WORD_CNT_; ++i) {
      v[i] = infScore ;
    }

    setWords(v) ;
  }

  static int getDistAt(int wordIndex, int queryLen) {
    return wordIndex * PERIOD_ + getProbeOffsetFor(queryLen) + 1 ;
  }

  void setWordsAsDist (int queryLen, int distWt) {
    int v[WORD_CNT_] ;
    int bitOffset = getProbeOffsetFor(queryLen) ;
    
    for (int wi=0; wi < WORD_CNT_; ++wi) {
      v[wi] = distWt * (wi * PERIOD_ + bitOffset + 1) ;
    }

    setWords(v) ;
  }

  void setMin(const EDVec64Every8& other) {
    vec_ = _mm_min_epi16(vec_, other.vec_) ;
  }

  void setMin(const EDVec64Every8& other, EDVec64Every8& bestIndices, const EDVec64Every8& otherIndices) {
    __m128i pred = _mm_cmplt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other < this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
  }

  void setMax(const EDVec64Every8& other, EDVec64Every8& bestIndices, const EDVec64Every8& otherIndices) {
    __m128i pred = _mm_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
  }

  void setMax(const EDVec64Every8& other, EDVec64Every8& bestIndices, const EDVec64Every8& otherIndices, EDVec64Every8& bestAccumDist, const EDVec64Every8& otherAccumDist) {
    __m128i pred = _mm_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
    bestAccumDist.vec_ = BLENDV(bestAccumDist.vec_, otherAccumDist.vec_, pred) ;
  }


  // The return value will be nonzero iff at least one of the entries is updated
  uint16_t setMaxAndReturnFlag(const EDVec64Every8& other, EDVec64Every8& bestIndices, const EDVec64Every8& otherIndices) {
    __m128i pred = _mm_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
    return (uint16_t) _mm_movemask_epi8(pred) ;
  }

  // The return value will be nonzero iff at least one of the entries is updated
  uint16_t setMaxAndReturnFlag(const EDVec64Every8& other, EDVec64Every8& bestIndices, const EDVec64Every8& otherIndices, EDVec64Every8& bestAccumDist, const EDVec64Every8& otherAccumDist) {
    __m128i pred = _mm_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
    bestAccumDist.vec_ = BLENDV(bestAccumDist.vec_, otherAccumDist.vec_, pred) ;
 
   return (uint16_t) _mm_movemask_epi8(pred) ;
  }



  void addThirdIfFirstGTSecond(const EDVec64Every8& first, const EDVec64Every8& second, 
			       const EDVec64Every8& third, const EDVec64Every8& zero) {

    __m128i pred = _mm_cmpgt_epi16(first.vec_, second.vec_) ; 
    // pred is 0xFFFF if first > second

    // If pred is 0, then blendv will add zero; o/w it will add third
    vec_ = _mm_add_epi16(vec_, BLENDV(zero.vec_, third.vec_, pred)) ;
  }


  /*
  EDVec64Every8 abs(const EDVec64Every8& other) const {
    return EDVec64Every8(_mm_abs_epi16(_mm_sub_epi16(this->vec_, other.vec_))) ;
  }

  EDVec64Every8 abs() const {
    return EDVec64Every8(_mm_abs_epi16(this->vec_)) ;
  }
  */

  EDVec64Every8 subSat(const EDVec64Every8& other) const {
    return EDVec64Every8(_mm_subs_epu16(this->vec_, other.vec_)) ;
  }

  inline EDVec64Every8 shiftBitsRightWithinWords(int shiftVal) {
    return EDVec64Every8(_mm_srli_epi16(vec_, shiftVal)) ;
  }

  inline void shiftWordsLeftByOne() {
    //vec_ = _mm_slli_si128(vec_, 2) ; // 16-bit shift = 2 * 8-bit shift

#ifdef USE_SSE4
    vec_ = _mm_shuffle_epi8(vec_, shiftLeftIndices_) ;
#else

    // The layout of the logical 16-bit words before shift:
    // 7 5 3 1 6 4 2 0
    // Consider the physical 16-bit words before shift: 
    // 7 6 5 4 3 2 1 0
    // These words should be shuffled into the following for logical shift:
    // 3 2 1 0 6 5 4 X
    
    // First, do the following: 7 6 5 4 3 2 1 0 --> 6 5 4 7 3 2 1 0
    // The binary imm value: 10 01 00 11 (i.e. 2 1 0 3)
    vec_ = _mm_shufflehi_epi16(vec_, 0x93) ;

    // Then, do the following: 6 5 4 7 3 2 1 0 --> 3 2 1 0 6 5 4 7
    // The binary imm value: 01 00 11 10 (i.e. 1 0 3 2)
    vec_ = _mm_shuffle_epi32(vec_, 0x4E) ;

#endif

  }

  inline bool operator == (const EDVec64Every8& other) {
    return _mm_movemask_epi8(_mm_cmpeq_epi16(vec_, other.vec_)) == 0xFFFF ;
  }

  inline EDVec64Every8& operator += (const EDVec64Every8& other) {
    vec_ = _mm_add_epi16(vec_, other.vec_) ;
    return *this ;
  }

  inline EDVec64Every8& operator -= (const EDVec64Every8& other) {
    vec_ = _mm_sub_epi16(vec_, other.vec_) ;
    return *this ;
  }

  inline EDVec64Every8 operator+ (const EDVec64Every8& other) const {
    return EDVec64Every8(_mm_add_epi16(vec_, other.vec_)) ;
  }

  inline EDVec64Every8 operator- (const EDVec64Every8& other) const {
    return EDVec64Every8(_mm_sub_epi16(vec_, other.vec_)) ;
  }

  inline EDVec64Every8 operator* (const EDVec64Every8& other) const {
    return EDVec64Every8(_mm_mullo_epi16(vec_, other.vec_)) ;
  }

  inline EDVec64Every8 operator& (const EDVec64Every8& other) const {
    return EDVec64Every8(_mm_and_si128(vec_, other.vec_)) ;
  }

  static int WORD_CNT() {
    return WORD_CNT_ ;
  }

  static int PERIOD() {
    return PERIOD_ ;
  }

  friend std::ostream& operator<< (std::ostream& os, const EDVec64Every8& v) {
    os << "{" ;
    for (int i=0; i < WORD_CNT_; ++i) {
      os << v.getWord(i) << " " ;
    }
    os << "}" ;
    return os ;
  }

} ;


//
// Edit distance vector for a 32-bit query.
// 16-bit words will be stored in an SSE vector in interleaved order as follows:
//           w7 w3 w6 w2 w5 w1 w4 w0
//
class EDVec32Every4 {
  __m128i vec_ ;
  __m128i shiftLeftIndices_ ; // can be made static later
  static const int WORD_CNT_ = 8 ;
  static const int PERIOD_ = 4 ;

  void init_() {
    shiftLeftIndices_ = _mm_set_epi8(11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 13, 12, 0xFF, 0xFF) ;
  }

 public:
  EDVec32Every4() {}

  explicit EDVec32Every4(__m128i v):vec_(v) {init_();}
  explicit EDVec32Every4(int16_t val):vec_(_mm_set1_epi16(val)) {init_();}
  explicit EDVec32Every4(const BitVec64& bv)
    : vec_(_mm_set_epi32(bv.bitV>>12, bv.bitV>>8, bv.bitV>>4, bv.bitV)) {init_();}

  inline void setAll(int16_t val) {
    vec_ = _mm_set1_epi16(val) ;    
  }

  static int16_t getMaxWordVal() { return numeric_limits<int16_t>::max() ;}

  inline int getWord(int index) const {
    switch(index) {
    case 0: return _mm_extract_epi16(vec_, 0) ;
    case 1: return  _mm_extract_epi16(vec_, 2) ;
    case 2: return  _mm_extract_epi16(vec_, 4) ;
    case 3: return  _mm_extract_epi16(vec_, 6) ;
    case 4: return  _mm_extract_epi16(vec_, 1) ;
    case 5: return  _mm_extract_epi16(vec_, 3) ;
    case 6: return  _mm_extract_epi16(vec_, 5) ;
    case 7: return  _mm_extract_epi16(vec_, 7) ;
    default: assert(false); return 0 ;
    }
  }

  bool allLessThanOrEqualTo(const EDVec32Every4& other) const {
    // each elt in pred will be set to 0xFFFF if the element in this vector 
    // is greater than the corresponding element in other
    __m128i pred = _mm_cmpgt_epi16(vec_, other.vec_) ; 
    int flags = _mm_movemask_epi8(pred) ; // creates a mask from MSB of each 8-bit elt in pred

    return flags == 0 ; // if there's at least one elt greater than "other", flags != 0
  }


  inline void setWords (int* v) {
    vec_ = _mm_set_epi16(v[7], v[3], v[6], v[2], v[5], v[1], v[4], v[0]) ;
  }

  inline void setFirstWord(int val) {
    vec_ = _mm_insert_epi16(vec_, val, 0) ;
  }

  inline void setWordsAsMask () {
    vec_ = _mm_set1_epi16(0x1) ;
  }

  static int getLastWordIndexFor (int queryLen) {
    return (queryLen-1) / PERIOD_ ; // logical word index - no interleaving at this level
  }

  static int getProbeOffsetFor (int queryLen) {
    return (queryLen-1) % PERIOD_ ;
  }

  void setWordsAsEndBonus(int queryLen, int endBonus) {
    int v[WORD_CNT_] = {0, 0, 0, 0, 0, 0, 0, 0} ;
    
    int lastWordIndex = getLastWordIndexFor(queryLen) ;    
    v[lastWordIndex] = endBonus ;

    setWords(v) ;
  }

  void setWordsAsBadScore(int queryLen, int thr, int infScore) {
    int v[WORD_CNT_] = {thr, thr, thr, thr, thr, thr, thr, thr} ;
    
    int lastWordIndex = getLastWordIndexFor(queryLen) ;    
    for (int i=lastWordIndex+1; i < WORD_CNT_; ++i)
      v[lastWordIndex] = infScore ;

    setWords(v) ;
  }

  static int getDistAt(int wordIndex, int queryLen) {
    return wordIndex * PERIOD_ + getProbeOffsetFor(queryLen) + 1 ;
  }

  void setWordsAsDist (int queryLen, int distWt) {
    int v[WORD_CNT_] ;
    int bitOffset = getProbeOffsetFor(queryLen) ;
    
    for (int wi=0; wi < WORD_CNT_; ++wi) {
      v[wi] = distWt * (wi * PERIOD_ + bitOffset + 1) ;
    }

    setWords(v) ;
  }

  void setMin(const EDVec32Every4& other) {
    vec_ = _mm_min_epi16(vec_, other.vec_) ;
  }

  void setMin(const EDVec32Every4& other, EDVec32Every4& bestIndices, const EDVec32Every4& otherIndices) {
    __m128i pred = _mm_cmplt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other < this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
  }

  void setMax(const EDVec32Every4& other, EDVec32Every4& bestIndices, const EDVec32Every4& otherIndices) {
    __m128i pred = _mm_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
  }

  void setMax(const EDVec32Every4& other, EDVec32Every4& bestIndices, const EDVec32Every4& otherIndices, EDVec32Every4& bestAccumDist, const EDVec32Every4& otherAccumDist) {
    __m128i pred = _mm_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
    bestAccumDist.vec_ = BLENDV(bestAccumDist.vec_, otherAccumDist.vec_, pred) ;
  }

  // The return value will be nonzero iff at least one of the entries is updated
  uint16_t setMaxAndReturnFlag(const EDVec32Every4& other, EDVec32Every4& bestIndices, const EDVec32Every4& otherIndices) {
    __m128i pred = _mm_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
    return (uint16_t) _mm_movemask_epi8(pred) ;
  }

  // The return value will be nonzero iff at least one of the entries is updated
  uint16_t setMaxAndReturnFlag(const EDVec32Every4& other, EDVec32Every4& bestIndices, const EDVec32Every4& otherIndices, EDVec32Every4& bestAccumDist, const EDVec32Every4& otherAccumDist) {
    __m128i pred = _mm_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
    bestAccumDist.vec_ = BLENDV(bestAccumDist.vec_, otherAccumDist.vec_, pred) ; 

   return (uint16_t) _mm_movemask_epi8(pred) ;
  }


  void addThirdIfFirstGTSecond(const EDVec32Every4& first, const EDVec32Every4& second, 
			       const EDVec32Every4& third, const EDVec32Every4& zero) {

    __m128i pred = _mm_cmpgt_epi16(first.vec_, second.vec_) ; 
    // pred is 0xFFFF if first > second

    // If pred is 0, then blendv will add zero; o/w it will add third
    vec_ = _mm_add_epi16(vec_, BLENDV(zero.vec_, third.vec_, pred)) ;
  }


  /*
  EDVec32Every4 abs(const EDVec32Every4& other) const {
    return EDVec32Every4(_mm_abs_epi16(_mm_sub_epi16(this->vec_, other.vec_))) ;
  }

  EDVec32Every4 abs() const {
    return EDVec32Every4(_mm_abs_epi16(this->vec_)) ;
  }
  */

  EDVec32Every4 subSat(const EDVec32Every4& other) const {
    return EDVec32Every4(_mm_subs_epu16(this->vec_, other.vec_)) ;
  }

  inline EDVec32Every4 shiftBitsRightWithinWords(int shiftVal) {
    return EDVec32Every4(_mm_srli_epi16(vec_, shiftVal)) ;
  }

  inline void shiftWordsLeftByOne() {
    //vec_ = _mm_slli_si128(vec_, 2) ; // 16-bit shift = 2 * 8-bit shift

#ifdef USE_SSE4
    vec_ = _mm_shuffle_epi8(vec_, shiftLeftIndices_) ;

#else

    // The layout of the logical 16-bit words before shift:
    // 7 3 6 2 5 1 4 0
    // Consider the physical 16-bit words before shift: 
    // 7 6 5 4 3 2 1 0
    // These words should be shuffled into the following for logical shift:
    // 5 4 3 2 1 0 6 X
    
    // First, do the following: 7 6 5 4 3 2 1 0 --> 5 4 3 2 1 0 7 6
    // The binary imm value: 10 01 00 11 (i.e. 2 1 0 3)
    vec_ = _mm_shuffle_epi32(vec_, 0x93) ;

    // Then, do the following: 5 4 3 2 1 0 7 6 --> 5 4 3 2 1 0 6 7
    // The binary imm value: 11 10 00 01 (i.e. 3 2 0 1)
    vec_ = _mm_shufflelo_epi16(vec_, 0xE1) ;
#endif

  }

  inline bool operator == (const EDVec32Every4& other) {
    //cout << "Comparing: " << endl ;
    //cout << *this << endl ;
    //cout << other << endl ;
    //cout << EDVec32Every4(_mm_cmpeq_epi16(vec_, other.vec_)) << endl ;
    //cout << std::hex << _mm_movemask_epi8(_mm_cmpeq_epi16(vec_, other.vec_))
    //     << std::dec << endl ;

    return _mm_movemask_epi8(_mm_cmpeq_epi16(vec_, other.vec_)) == 0xFFFF ;
  }

  inline EDVec32Every4& operator += (const EDVec32Every4& other) {
    vec_ = _mm_add_epi16(vec_, other.vec_) ;
    return *this ;
  }

  inline EDVec32Every4& operator -= (const EDVec32Every4& other) {
    vec_ = _mm_sub_epi16(vec_, other.vec_) ;
    return *this ;
  }

  inline EDVec32Every4 operator+ (const EDVec32Every4& other) const {
    return EDVec32Every4(_mm_add_epi16(vec_, other.vec_)) ;
  }

  inline EDVec32Every4 operator- (const EDVec32Every4& other) const {
    return EDVec32Every4(_mm_sub_epi16(vec_, other.vec_)) ;
  }

  inline EDVec32Every4 operator* (const EDVec32Every4& other) const {
    return EDVec32Every4(_mm_mullo_epi16(vec_, other.vec_)) ;
  }

  inline EDVec32Every4 operator& (const EDVec32Every4& other) const {
    return EDVec32Every4(_mm_and_si128(vec_, other.vec_)) ;
  }

  static int WORD_CNT() {
    return WORD_CNT_ ;
  }

  static int PERIOD() {
    return PERIOD_ ;
  }

  friend std::ostream& operator<< (std::ostream& os, const EDVec32Every4& v) {
    os << "{" ;
    for (int i=0; i < WORD_CNT_; ++i) {
      os << v.getWord(i) << " " ;
    }
    os << "}" ;
    return os ;
  }

} ;

#if 0
//
// Edit distance vector for a 16-bit query.
// 8-bit words will be stored in an SSE vector in interleaved order as follows:
//           w7 w3 w6 w2 w5 w1 w4 w0
//
class EDVec16Every1 {
  __m128i vec_ ;
  static const int WORD_CNT_ = 16 ;
  static const int PERIOD_ = 1 ;

  static uint64_t expand8To64LUT_[256] ; // to expand 8 bits to 64 bits

 public:
  EDVec16Every1() {}

  explicit EDVec16Every1(__m128i v):vec_(v) {}
  explicit EDVec16Every1(int16_t val):vec_(_mm_set1_epi8((int8_t) val)) {}

  explicit EDVec16Every1(const BitVec64& bv)
    :vec_(_mm_set_epi64(expand8To64LUT_[bv.bitV >> 8], expand8To64LUT_[bv.bitV])) {}

  /*
  explicit EDVec16Every1(const BitVec64& bv):vec_(_mm_set_epi64(bv.bitV >> 8, bv.bitV)) {
    // The initialization above separates the 16-bit vector into upper and lower 64-bits.
    // The 8-bit words of 128-bit vec look like this:
    // _ _ _ _ _ _ _ b[15:8] _ _ _ _ _ _ _ b[7:0]

    // Further separate them into 4-bit words:
    vec_ = _mm_or_si128(vec_, _mm_slli_si128(_mm_srli_epi16(vec_, 4), 4)) ;
    // Now, we have:
    // _ _ _ b[15:12] _ _ _ b[11:8] _ _ _ b[7:4] _ _ _ b[3:0]
    
    // Further separate them into 2-bit words:
    vec_ = _mm_or_si128(vec_, _mm_slli_si128(_mm_srli_epi16(vec_, 2), 2)) ;
    // Now, we have:
    // _ b[15:14] _ b[13:12] _ b[11:10] _ b[9:8] _ b[7:6] _ b[5:4] _ b[3:2] _ b[1:0]

    // Further separate them into 1-bit words:
    vec_ = _mm_or_si128(vec_, _mm_slli_si128(_mm_srli_epi16(vec_, 1), 1)) ;

  }
  */

  static void initLUT() ;

  inline void setAll(int16_t val) {
    vec_ = _mm_set1_epi8((int8_t)val) ;    
  }

  static int16_t getMaxWordVal() { return numeric_limits<int8_t>::max() ;}

  inline int getWord(int index) const {
    switch(index) {
    case 0: return _mm_extract_epi8(vec_, 0) ;
    case 1: return  _mm_extract_epi8(vec_, 1) ;
    case 2: return  _mm_extract_epi8(vec_, 2) ;
    case 3: return  _mm_extract_epi8(vec_, 3) ;
    case 4: return  _mm_extract_epi8(vec_, 4) ;
    case 5: return  _mm_extract_epi8(vec_, 5) ;
    case 6: return  _mm_extract_epi8(vec_, 6) ;
    case 7: return  _mm_extract_epi8(vec_, 7) ;
    case 8: return  _mm_extract_epi8(vec_, 8) ;
    case 9: return  _mm_extract_epi8(vec_, 9) ;
    case 10: return  _mm_extract_epi8(vec_, 10) ;
    case 11: return  _mm_extract_epi8(vec_, 11) ;
    case 12: return  _mm_extract_epi8(vec_, 12) ;
    case 13: return  _mm_extract_epi8(vec_, 13) ;
    case 14: return  _mm_extract_epi8(vec_, 14) ;
    case 15: return  _mm_extract_epi8(vec_, 15) ;

    default: assert(false); return 0 ;
    }
  }

  inline void setWords (int* v) {
    vec_ = _mm_set_epi8(v[15], v[14], v[13], v[12], v[11], v[10], v[9], v[8],
			v[7], v[6], v[5], v[4], v[3], v[2], v[1], v[0]) ;
  }

  inline void setFirstWord(int val) {
    vec_ = _mm_insert_epi16(vec_, val, 0) ;
  }

  inline void setWordsAsMask () {
    vec_ = _mm_set1_epi8(0x1) ;
  }

  static int getLastWordIndexFor (int queryLen) {
    return (queryLen-1) / PERIOD_ ; // logical word index - no interleaving at this level
  }

  static int getProbeOffsetFor (int queryLen) {
    return (queryLen-1) % PERIOD_ ;
  }

  void setWordsAsEndBonus(int queryLen, int endBonus) {
    int v[WORD_CNT_] = {0, 0, 0, 0, 0, 0, 0, 0} ;
    
    int lastWordIndex = getLastWordIndexFor(queryLen) ;    
    v[lastWordIndex] = endBonus ;

    setWords(v) ;
  }

  void setWordsAsBadScore(int queryLen, int thr, int infScore) {
    int v[WORD_CNT_] = {thr, thr, thr, thr, thr, thr, thr, thr} ;
    
    int lastWordIndex = getLastWordIndexFor(queryLen) ;    
    for (int i=lastWordIndex+1; i < WORD_CNT_; ++i)
      v[lastWordIndex] = infScore ;

    setWords(v) ;
  }

  static int getDistAt(int wordIndex, int queryLen) {
    return wordIndex * PERIOD_ + getProbeOffsetFor(queryLen) + 1 ;
  }

  void setWordsAsDist (int queryLen, int distWt) {
    int v[WORD_CNT_] ;
    int bitOffset = getProbeOffsetFor(queryLen) ;
    
    for (int wi=0; wi < WORD_CNT_; ++wi) {
      v[wi] = distWt * (wi * PERIOD_ + bitOffset + 1) ;
    }

    setWords(v) ;
  }

  void setMin(const EDVec16Every1& other) {
    vec_ = _mm_min_epi8(vec_, other.vec_) ;
  }

  void setMin(const EDVec16Every1& other, EDVec16Every1& bestIndices, const EDVec16Every1& otherIndices) {
    __m128i pred = _mm_cmplt_epi8(other.vec_, vec_) ; // pred is 0xFFFF if other < this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
  }

  void setMax(const EDVec16Every1& other, EDVec16Every1& bestIndices, const EDVec16Every1& otherIndices) {
    __m128i pred = _mm_cmpgt_epi8(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
  }

  void setMax(const EDVec16Every1& other, EDVec16Every1& bestIndices, const EDVec16Every1& otherIndices, EDVecEvery1& bestAccumDist, const EDVecEvery1& otherAccumDist) {
    __m128i pred = _mm_cmpgt_epi8(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
    bestAccumDist.vec_ = BLENDV(bestAccumDist.vec_, otherAccumDist.vec_, pred) ;    
  }


  // The return value will be nonzero iff at least one of the entries is updated
  uint16_t setMaxAndReturnFlag(const EDVec16Every1& other, EDVec16Every1& bestIndices, const EDVec16Every1& otherIndices) {
    __m128i pred = _mm_cmpgt_epi8(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
    return (uint16_t) _mm_movemask_epi8(pred) ;
  }

  EDVec16Every1 subSat(const EDVec16Every1& other) const {
    return EDVec16Every1(_mm_subs_epu8(this->vec_, other.vec_)) ;
  }

  inline bool operator == (const EDVec16Every1& other) {
    return _mm_movemask_epi8(_mm_cmpeq_epi16(vec_, other.vec_)) == 0xFFFFFFFF ;
  }

  inline EDVec16Every1& operator += (const EDVec16Every1& other) {
    vec_ = _mm_add_epi8(vec_, other.vec_) ;
    return *this ;
  }

  inline EDVec16Every1& operator -= (const EDVec16Every1& other) {
    vec_ = _mm_sub_epi8(vec_, other.vec_) ;
    return *this ;
  }

  inline EDVec16Every1 operator+ (const EDVec16Every1& other) const {
    return EDVec16Every1(_mm_add_epi8(vec_, other.vec_)) ;
  }

  inline EDVec16Every1 operator- (const EDVec16Every1& other) const {
    return EDVec16Every1(_mm_sub_epi8(vec_, other.vec_)) ;
  }

  inline EDVec16Every1 operator* (const EDVec16Every1& other) const {
    return EDVec16Every1(_mm_mullo_epi8(vec_, other.vec_)) ;
  }

  inline EDVec16Every1 operator& (const EDVec16Every1& other) const {
    return EDVec16Every1(_mm_and_si128(vec_, other.vec_)) ;
  }

  static int WORD_CNT() {
    return WORD_CNT_ ;
  }

  static int PERIOD() {
    return PERIOD_ ;
  }

  friend std::ostream& operator<< (std::ostream& os, const EDVec16Every1& v) {
    os << "{" ;
    for (int i=0; i < WORD_CNT_; ++i) {
      os << v.getWord(i) << " " ;
    }
    os << "}" ;
    return os ;
  }

} ;

#endif // #if 0



#endif // #ifndef _ED_INTRAVED_H
