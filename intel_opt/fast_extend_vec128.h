/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _FAST_EXTEND_VEC128_H
#define _FAST_EXTEND_VEC128_H

#include "fast_extend_vec.h"
#include "fast_extend_bitv64.h"
#include "fast_extend_bitv64x2.h"
#include "fast_extend_bitv128.h"

//
// Edit distance vector for a 128-bit query.
// The probes are at every 16 bases. Each probe is 16 bits.
// There are 8 probes in total.
//
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
  explicit EDVec128Every16(const BitVec128& bv): vec_(bv.bitV) {}


  inline void setAll(int16_t val) {
    vec_ = _mm_set1_epi16(val) ;    
  }

  static int16_t getMaxWordVal() { return std::numeric_limits<int16_t>::max() ;}

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
    // Each elt in pred will be set to 0xFFFF if the element in this vector 
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

  inline void setAllWords (uint16_t val) {
    vec_ = _mm_set1_epi16(val) ;
  }

  inline void setWordsAsMask () {
    vec_ = _mm_set1_epi16(0x1) ;
  }

  static int getLastWordIndexFor (int queryLen) {
    int baseIndex = queryLen - 1 ;
    return baseIndex/PERIOD_ ;
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
    return 1+ wordIndex * PERIOD_ + getProbeOffsetFor(queryLen) ;
  }

  void setWordsAsDist (int queryLen, int distWt) {
    int v[WORD_CNT_] ;
    int bitOffset = getProbeOffsetFor(queryLen) ;

    for (int wi=0; wi < WORD_CNT_; ++wi) {
      v[wi] = distWt * (1 + wi * PERIOD_ + bitOffset) ;
    }
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
  // 2 bits per word is returned because of SSE limitations.
  uint32_t setMaxAndReturnFlag(const EDVec128Every16& other, EDVec128Every16& bestIndices, const EDVec128Every16& otherIndices) {
    __m128i pred = _mm_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
    
    return (uint32_t) _mm_movemask_epi8(pred) ;
  }

  // The return value will be nonzero iff at least one of the entries is updated
  // 2 bits per word is returned because of SSE limitations.
  uint32_t setMaxAndReturnFlag(const EDVec128Every16& other, EDVec128Every16& bestIndices, const EDVec128Every16& otherIndices, EDVec128Every16& bestAccumDist, const EDVec128Every16& otherAccumDist) {
    __m128i pred = _mm_cmpgt_epi16(other.vec_, vec_) ; // pred is 0xFFFF if other > this
    vec_ = BLENDV(vec_, other.vec_, pred) ; // if pred is 0, choose this; o/w choose other
    bestIndices.vec_ = BLENDV(bestIndices.vec_, otherIndices.vec_, pred) ;
    bestAccumDist.vec_ = BLENDV(bestAccumDist.vec_, otherAccumDist.vec_, pred) ;
    
    return (uint32_t) _mm_movemask_epi8(pred) ;
  }

  void addThirdIfFirstGTSecond(const EDVec128Every16& first, const EDVec128Every16& second, 
			       const EDVec128Every16& third, const EDVec128Every16& zero) {

    __m128i pred = _mm_cmpgt_epi16(first.vec_, second.vec_) ; 
    // pred is 0xFFFF if first > second

    // If pred is 0, then blendv will add zero; o/w it will add third
    vec_ = _mm_add_epi16(vec_, BLENDV(zero.vec_, third.vec_, pred)) ;
  }
 
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



#endif // #ifndef _FAST_EXTEND_VEC128_H
