/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/


#ifndef _ED_INTRAV64_H
#define _ED_INTRAV64_H


#include <stdint.h>
#include <stdlib.h>
#include <limits>

class BitVec64 {

  uint64_t bitV ;

  explicit BitVec64(uint64_t val): bitV(val) {}
  //  cout << "val = " << val << ", bitV = " << bitV << endl ;

public:
  
  BitVec64(): bitV(0x0) {}

  inline bool getLSB() const {
    return bitV & 0x1 ;
  }
  
  inline uint16_t getLow16Bits() const {
    return (uint16_t) bitV ;
  }

  inline void setAllOnes () {
    bitV = 0xFFFFFFFFFFFFFFFF ;
  }

  inline void setAllZeroes () {
    bitV = 0x0000000000000000 ;
  }

  inline void setAllBits (bool bit) {
    bitV = bit ? 0xFFFFFFFFFFFFFFFF : 0x0 ;
  }

  inline void setLSBClearRest (bool bit) {
    bitV = bit ;
  }

  inline void setBit(int bitIndex, bool bit) {
    uint64_t mask = ((uint64_t)0x1) << bitIndex ;
    bitV = bit ? (bitV | mask) : (bitV & ~mask) ;
  }

  inline void shiftLeftAndInsert (bool newBit) {
    bitV = (bitV << 1) | ((uint8_t) newBit) ;
  }

  inline void shiftLeft () {
    bitV <<= 1 ;
  }

  /*
  inline BitVec64 shiftBitsRight(int offset) const {
    return BitVec64(bitV >> offset) ;
  }
  */

  inline BitVec64 shiftProbesRight(int probeOffset) const {
    // The carryout bit is at location probeOffset+1
    return BitVec64(bitV >> (probeOffset+1)) ;
  }


  void setWord (int index, uint64_t w) {
    bitV = w ;
  }

  inline void orBitAtPosition (bool newBit, int index) {
    bitV |= (((uint64_t)newBit) << index) ;
  }
  
  inline bool getMSB() const {
    return bitV & 0x8000000000000000 ;
  }

  inline BitVec64 operator~ () const {
    return BitVec64(~this->bitV) ;
  }

  inline BitVec64 operator& (const BitVec64& other) const {
    return BitVec64(this->bitV & other.bitV) ;
  }

  inline BitVec64 andnot (const BitVec64& other) const {
    return *this & (~other) ;
  }

  inline BitVec64 operator| (const BitVec64& other) const {
    return BitVec64(this->bitV | other.bitV) ;
  }

  inline BitVec64 operator^ (const BitVec64& other) const {
    return BitVec64(this->bitV ^ other.bitV) ;
  }

  inline BitVec64 operator+ (const BitVec64& other) const {
    return BitVec64(this->bitV + other.bitV) ;
  }

  //inline BitVec64 operator>>(int shiftVal) const {
  //  return BitVec64(this->bitV >> shiftVal) ;
  //}

  static bool isQueryLengthOk(int qlen) {
    return qlen < 64 ;
  }

  bool getBit (int index) const {
    return (bitV >> index) & 0x1 ;
  }

  friend std::ostream& operator<< (std::ostream& os, const BitVec64& b) {
    os << b.bitV ;
    return os ;
  }

  friend class EDVec64 ;
  friend class EDVec128Every16 ;
  friend class EDVec64Every8 ;
  friend class EDVec32Every4 ;
  friend class EDVec16Every1 ;

} ;


class EDVec64 {
  int64_t val_ ;
  
 public:
  EDVec64() {}
  explicit EDVec64(int val): val_(val) {}
  explicit EDVec64(const BitVec64& bv): val_(bv.bitV) {}

  inline void setAll(int val) {
    val_ = val ;
  }

  static int16_t getMaxWordVal() { return std::numeric_limits<int16_t>::max() ;}

  inline int getWord(int index) const {
    return (int) val_ ;
  }

  bool allLessThanOrEqualTo(const EDVec64& other) const {
    return val_ <= other.val_ ;
  }

  inline void setWords (int* vals) {
    val_ = vals[0] ;
  }

  inline void setFirstWord(int val) {
    val_ = val ;
  }

  inline void setWordsAsMask () {
    val_ = 0x1 ;
  }

  static int getLastWordIndexFor(int queryLen) {
    return 0 ;
  }
  
  static int getProbeOffsetFor (int queryLen) {
    return queryLen-1 ;
  }

  static int getBitOffsetFor(int baseIndex) {
    return baseIndex ;
  }

  void setWordsAsEndBonus(int queryLen, int endBonus) {
    val_ = endBonus ;
  }
 

  void setWordsAsDist (int queryLen, int distWt) {
    val_ = queryLen * distWt ;
  }

  void setWordsAsBadScore(int queryLen, int thr, int infScore) {
    val_ = infScore ;
  }

  static int getDistAt(int wordIndex, int queryLen) {
    return queryLen ;
  }

  void setMin(const EDVec64& other) {
    this->val_ = std::min(this->val_, other.val_) ;
  }

  void setMin(const EDVec64& other, EDVec64& bestIndices, const EDVec64& otherIndices) {
    if (other.val_ < this->val_) {
      this->val_ = other.val_ ;
      bestIndices.val_ = otherIndices.val_ ;
    }
  }

  void setMax(const EDVec64& other, EDVec64& bestIndices, const EDVec64& otherIndices) {
    if (other.val_ > this->val_) {
      this->val_ = other.val_ ;
      bestIndices.val_ = otherIndices.val_ ;
    }
  }

  uint16_t setMaxAndReturnFlag(const EDVec64& other, EDVec64& bestIndices, const EDVec64& otherIndices) {
    if (other.val_ > this->val_) {
      this->val_ = other.val_ ;
      bestIndices.val_ = otherIndices.val_ ;
      return 0x1 ;
    }
    
    return 0x0 ;
  }
  

  void addThirdIfFirstGTSecond(const EDVec64& first, const EDVec64& second, 
			       const EDVec64& third, const EDVec64& zero) {
    if (first.val_ > second.val_)
      this->val_ += third.val_ ;
  }

  /*
  EDVec64 abs(const EDVec64& other) const {
    return EDVec64(labs((uint64_t)this->val_ - other.val_)) ;
  }

  EDVec64 abs() const {
    return EDVec64(labs(this->val_)) ;
  }
  */

  EDVec64 subSat(const EDVec64& other) const {
    int64_t diff = this->val_ - other.val_ ;

    return EDVec64(diff < 0 ? 0 : diff) ;
  }

  inline EDVec64 shiftBitsRightWithinWords(int shiftVal) {
    return EDVec64(this->val_ >> shiftVal) ;
  }

  inline void shiftWordsLeftByOne() {
    val_ = 0 ;
  }

  inline bool operator == (const EDVec64& other) const {
    return this->val_ == other.val_ ;
  }

  inline EDVec64& operator += (const EDVec64& other) {
    this->val_ += other.val_ ;
    return *this ;
  }

  inline EDVec64& operator -= (const EDVec64& other) {
    this->val_ -= other.val_ ;
    return *this ;
  }

  inline EDVec64 operator+ (const EDVec64& other) const {
    return EDVec64(this->val_ + other.val_) ;
  }

  inline EDVec64 operator- (const EDVec64& other) const {
    return EDVec64(this->val_ - other.val_) ;
  }

  inline EDVec64 operator* (const EDVec64& other) const {
    return EDVec64(this->val_ * other.val_) ;
  }

  inline EDVec64 operator& (const EDVec64& other) const {
    return EDVec64(this->val_ & other.val_) ;
  }

  static int PERIOD() {
    return 64 ;
  }

  static int WORD_CNT() {
    return 1 ;
  }

  friend std::ostream& operator<< (std::ostream& os, const EDVec64& v) {
    os << v.val_ ;
    return os ;
  }

} ;


#endif // #ifndef _ED_INTRAV64_H
