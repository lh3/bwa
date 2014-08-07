/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#ifndef _ED_INTRAV64x2_H
#define _ED_INTRAV64x2_H


#include <stdint.h>
#include <stdlib.h>
#include <iostream>

using namespace std ;

class BitVec64x2 {

  uint64_t bitV[2] ;


public:
  
  BitVec64x2() {
    bitV[0] = 0x0 ;
    bitV[1] = 0x0 ;
  }

  inline bool getLSB() const {
    return bitV[0] & 0x1 ;
  }

  inline uint16_t getLow16Bits() const {
    return (uint16_t) bitV[0] ;
  }

  inline void setAllOnes () {
    bitV[0] = 0xFFFFFFFFFFFFFFFF ;
    bitV[1] = 0xFFFFFFFFFFFFFFFF ;
  }

  inline void setAllZeroes () {
    bitV[0] = 0x0000000000000000 ;
    bitV[1] = 0x0000000000000000 ;
  }

  inline void setAllBits (bool bit) {
    if (bit) setAllOnes() ;
    else setAllZeroes() ;
  }

  inline void setLSBClearRest (bool bit) {
    bitV[0] = bit ;
    bitV[1] = 0x0 ;
  }

  inline void shiftLeftAndInsert (bool newBit) {
    bitV[0] <<= 1 ;
    bitV[1] <<= 1 ;
    bitV[1] |= (bitV[0] >> 63) ;
    bitV[0] |= ((uint8_t) newBit) ;
  }

  inline void shiftLeft () {
    bitV[0] <<= 1 ;
    bitV[1] <<= 1 ;
    bitV[1] |= (bitV[0] >> 63) ;
  }

  /*
  inline BitVec64x2 shiftBitsRight(int offset) const {
    BitVec64x2 r ;
    r.bitV[0] = this->bitV[0] >> offset ;
    r.bitV[1] = this->bitV[1] >> offset ;
    return r ;
  }
  */

  inline BitVec64x2 shiftProbesRight(int probeOffset) const {
    BitVec64x2 r ;
    
    // Copy LSB of bitV[1] to MSB of bitV[0] in case probeOffset is 15
    // L' of base 62 is equal to L of base 63, which is stored at the LSB of bitV[1]
    r.bitV[0] = (this->bitV[0] & 0x7FFFFFFFFFFFFFFF) | (this->bitV[1] << 63) ;

    // The carryout bit is at location probeOffset+1
    r.bitV[0] = r.bitV[0] >> probeOffset ; // probes in the lower word have offset-1
    r.bitV[1] = this->bitV[1] >> (probeOffset+1) ; 
    return r ;
  }


  void setWord (int index, uint64_t w) {
    bitV[index] = w ;
  }

  inline BitVec64x2 operator~ () const {
    BitVec64x2 r ;
    r.bitV[0] = ~this->bitV[0] ;
    r.bitV[1] = ~this->bitV[1] ;
    return r ;
  }

  inline BitVec64x2 operator& (const BitVec64x2& other) const {
    BitVec64x2 r ;
    r.bitV[0] = this->bitV[0] & other.bitV[0] ;
    r.bitV[1] = this->bitV[1] & other.bitV[1] ;
    return r ;
  }

  inline BitVec64x2 andnot (const BitVec64x2& other) const {
    return *this & (~other) ;
  }

  inline BitVec64x2 operator| (const BitVec64x2& other) const {
    BitVec64x2 r ;
    r.bitV[0] = this->bitV[0] | other.bitV[0] ;
    r.bitV[1] = this->bitV[1] | other.bitV[1] ;
    return r ;
  }

  inline BitVec64x2 operator^ (const BitVec64x2& other) const {
    BitVec64x2 r ;
    r.bitV[0] = this->bitV[0] ^ other.bitV[0] ;
    r.bitV[1] = this->bitV[1] ^ other.bitV[1] ;
    return r ;
  }

  inline BitVec64x2 operator+ (const BitVec64x2& other) const {
    BitVec64x2 r ;
    r.bitV[0] = (this->bitV[0] & 0x7FFFFFFFFFFFFFFF) 
      +         (other.bitV[0] & 0x7FFFFFFFFFFFFFFF) ;

    r.bitV[1] = this->bitV[1] + other.bitV[1] ;

    r.bitV[1] += ((r.bitV[0] >> 63) & 0x1) ;

    /*
    cout << "Add operation: " << endl ;
    cout << std::hex << this->bitV[0] << " + " << other.bitV[0] << endl ;
    cout << this->bitV[1] << " + " << other.bitV[1] << endl ;
    cout << "Output: " << endl ;
    cout << r.bitV[0] << endl ;
    cout << r.bitV[1] << std::dec << endl ;
    */

    return r ;
  }


  static bool isQueryLengthOk(int qlen) {
    return qlen <= 126 ;
  }

  bool getBit (int bitIndex) const {
    const uint64_t& bitVchosen = (bitIndex < 63) ? bitV[0] : bitV[1] ;
    int bitIndexMod = (bitIndex < 63) ? bitIndex : (bitIndex - 63) ;

    return (bitVchosen >> bitIndexMod) & 0x1 ;
  }

  friend std::ostream& operator<< (std::ostream& os, const BitVec64x2& b) {
    os << "{" << b.bitV[0] << ", " << b.bitV[1] << "}" ;
    return os ;
  }

  friend class EDVec64x2 ;
  friend class EDVec128Every16 ;
} ;



class EDVec64x2 {
  int64_t vals_[2] ;

 public:
  EDVec64x2() {} 

  explicit EDVec64x2(int val) {
    vals_[0] = val ;
    vals_[1] = val ;
  }
  
  explicit EDVec64x2(const BitVec64x2& bv) {
    vals_[0] = bv.bitV[0] ;
    vals_[1] = bv.bitV[1] ;
  }

  inline void set(int val, int index) {
    vals_[index] = val ;
  }

  static int16_t getMaxWordVal() { return numeric_limits<int16_t>::max() ;}

  inline int getWord(int index) const {
    return (int) vals_[index] ;
  }

  bool allLessThanOrEqualTo(const EDVec64x2& other) const {
    return vals_[0] <= other.vals_[0] && vals_[1] <= other.vals_[1] ;
  }

  inline void setWords (int* vals) {
    vals_[0] = vals[0] ;
    vals_[1] = vals[1] ;
  }

  inline void setFirstWord(int val) {
    vals_[0] = val ;
  }

  inline void setWordsAsMask () {
    vals_[0] = 0x1 ;
    vals_[1] = 0x1 ;
  }

  static int getLastWordIndexFor(int queryLen) {
    int baseIndex = queryLen - 1 ;
    return (baseIndex < 63) ? 0 : 1 ;
  }

  static int getProbeOffsetFor (int queryLen) {
    return (queryLen-1) % 63 ;
  }
  
  static int getBitOffsetFor(int baseIndex) {
    return (baseIndex < 63) ? baseIndex : (baseIndex - 63) ;
  }

  void setWordsAsEndBonus(int queryLen, int endBonus) {
    vals_[0] = vals_[1] = 0 ;
    if (queryLen <= 63)
      vals_[0] = endBonus ;
    else
      vals_[1] = endBonus ;
  }


  inline void setWordsAsDist(int queryLen, int distWt) {
    int bitOffset = getBitOffsetFor(queryLen-1) ;
    vals_[0] = (bitOffset + 1) * distWt ;
    vals_[1] = (63 + bitOffset + 1) * distWt ;
  }

  void setWordsAsBadScore(int queryLen, int thr, int infScore) {
    vals_[0] = thr ;
    vals_[1] = thr ;
    vals_[getLastWordIndexFor(queryLen-1)] = infScore ;
  }

  static int getDistAt(int wordIndex, int queryLen) {
    return wordIndex * PERIOD() + getProbeOffsetFor(queryLen) ;
  }

  void setMin(const EDVec64x2& other) {
    this->vals_[0] = min(this->vals_[0], other.vals_[0]) ;
    this->vals_[1] = min(this->vals_[1], other.vals_[1]) ;
  }
  
  void setMin(const EDVec64x2& other, EDVec64x2& bestIndices, const EDVec64x2& otherIndices) {
    for (int i=0; i < 2; ++i) {
      if (other.vals_[i] < this->vals_[i]) {
	this->vals_[i] = other.vals_[i] ;
	bestIndices.vals_[i] = otherIndices.vals_[i] ;
      }
    }
  }

  void setMax(const EDVec64x2& other, EDVec64x2& bestIndices, const EDVec64x2& otherIndices) {
    for (int i=0; i < 2; ++i) {
      if (other.vals_[i] > this->vals_[i]) {
	this->vals_[i] = other.vals_[i] ;
	bestIndices.vals_[i] = otherIndices.vals_[i] ;
      }
    }
  }

  uint16_t setMaxAndReturnFlag(const EDVec64x2& other, EDVec64x2& bestIndices, const EDVec64x2& otherIndices) {

    uint16_t flag = 0x0 ;
    
    for (int i=0; i < 2; ++i) {
      if (other.vals_[i] > this->vals_[i]) {
	this->vals_[i] = other.vals_[i] ;
	bestIndices.vals_[i] = otherIndices.vals_[i] ;
	flag |= (0x1 << i) ;
      }
    }

    return flag ;
  }


  void addThirdIfFirstGTSecond(const EDVec64x2& first, const EDVec64x2& second, 
			       const EDVec64x2& third, const EDVec64x2& zero) {
    
    for (int i=0; i < 2; ++i) {
      if (first.vals_[i] > second.vals_[i])
	this->vals_[i] += third.vals_[i] ;
    }
  }


  /*
  EDVec64x2 abs(const EDVec64x2& other) const {
    EDVec64x2 r ;
    r.vals_[0] = labs(this->vals_[0] - other.vals_[0]) ;
    r.vals_[1] = labs(this->vals_[1] - other.vals_[1]) ;

    return r ;
  }

  EDVec64x2 abs() const {
    EDVec64x2 r ;
    r.vals_[0] = labs(this->vals_[0]) ;
    r.vals_[1] = labs(this->vals_[1]) ;
    return r ;
  }
  */

  EDVec64x2 subSat(const EDVec64x2& other) const {
    EDVec64x2 r ;
    int64_t diff0 = this->vals_[0] - other.vals_[0] ;
    int64_t diff1 = this->vals_[1] - other.vals_[1] ;
    
    r.vals_[0] = (diff0 < 0 ? 0 : diff0) ;
    r.vals_[1] = (diff1 < 0 ? 0 : diff1) ;
    return r ;
  }

  inline EDVec64x2 shiftBitsRightWithinWords(int shiftVal) {
    EDVec64x2 r ;
    r.vals_[0] = this->vals_[0] >> shiftVal ;
    r.vals_[1] = this->vals_[1] >> shiftVal ;
    return r ;
  }

  inline void shiftWordsLeftByOne() {
    vals_[1] = vals_[0] ;
    vals_[0] = 0 ;
  }

  inline bool operator == (const EDVec64x2& other) const {
    return this->vals_[0] == other.vals_[0] && this->vals_[1] == other.vals_[1] ;
  }

  inline EDVec64x2& operator += (const EDVec64x2& other) {
    this->vals_[0] += other.vals_[0] ;
    this->vals_[1] += other.vals_[1] ;
    return *this ;
  }

  inline EDVec64x2& operator -= (const EDVec64x2& other) {
    this->vals_[0] -= other.vals_[0] ;
    this->vals_[1] -= other.vals_[1] ;
    return *this ;
  }


  inline EDVec64x2 operator+ (const EDVec64x2& other) const {
    EDVec64x2 r ;
    r.vals_[0] = this->vals_[0] + other.vals_[0] ;
    r.vals_[1] = this->vals_[1] + other.vals_[1] ;
    return r ;
  }

  inline EDVec64x2 operator- (const EDVec64x2& other) const {
    EDVec64x2 r ;
    r.vals_[0] = this->vals_[0] - other.vals_[0] ;
    r.vals_[1] = this->vals_[1] - other.vals_[1] ;
    return r ;
  }


  inline EDVec64x2 operator* (const EDVec64x2& other) const {
    EDVec64x2 r ;
    r.vals_[0] = this->vals_[0] * other.vals_[0] ;
    r.vals_[1] = this->vals_[1] * other.vals_[1] ;
    return r ;
  }

  inline EDVec64x2 operator& (const EDVec64x2& other) const {
    EDVec64x2 r ;
    r.vals_[0] = this->vals_[0] & other.vals_[0] ;
    r.vals_[1] = this->vals_[1] & other.vals_[1] ;
    return r ;
  }

  static int PERIOD() {
    return 64 ;
  }

  static int WORD_CNT() {
    return 2 ;
  }

  friend std::ostream& operator<< (std::ostream& os, const EDVec64x2& v) {
    os << "{" << v.vals_[0] << ", " << v.vals_[1] << "}" ;
    return os ;
  }

} ;

#endif // #ifndef _ED_INTRAV64x2_H

