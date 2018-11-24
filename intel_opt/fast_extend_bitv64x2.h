/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _FAST_EXTEND_BITV64x2_H
#define _FAST_EXTEND_BITV64x2_H

#include <stdint.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

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

  inline void setAll64BitWords(uint64_t val) {
    bitV[0] = bitV[1] = val ;
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
    bitV[1] <<= 1 ;
    bitV[1] |= (bitV[0] >> 63) ;
    bitV[0] <<= 1 ;
    bitV[0] |= ((uint8_t) newBit) ;
  }

  inline void shiftLeft () {
    bitV[1] <<= 1 ;
    bitV[1] |= (bitV[0] >> 63) ;
    bitV[0] <<= 1 ;
  }

  
  inline BitVec64x2 shiftProbesRight(int probeOffset, const BitVec64x2& shiftCnt) const {
    BitVec64x2 r ;

    r.bitV[0] = bitV[0] >> (probeOffset + 1) ;

    if (probeOffset == 15) 
      r.bitV[0] |= ((bitV[1] & 0x1) << 48) ; 

    r.bitV[1] = bitV[1] >> (probeOffset + 1) ;
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
    r.bitV[0] = this->bitV[0] + other.bitV[0] ;
    r.bitV[1] = this->bitV[1] + other.bitV[1] + (r.bitV[0] < this->bitV[0]) ;
    return r ;
  }


  static bool isQueryLengthOk(int qlen) {
    return qlen <= 127 ;
  }

  friend std::ostream& operator<< (std::ostream& os, const BitVec64x2& b) {
    os << "{" << b.bitV[0] << ", " << b.bitV[1] << "}" ;
    return os ;
  }

  friend class EDVec64x2 ;
  friend class EDVec128Every16 ;
} ;



#endif // #ifndef _FAST_EXTEND_BITV64x2_H
