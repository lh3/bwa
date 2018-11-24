/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include <stdint.h>
#include <iostream>

// Alternatively, can use popcnt instruction
class BitCount8 {
  static uint8_t lut_[256] ;

 public:

  static void init_lut() ;
  
  static int getCnt(uint8_t bitVec) {
    return lut_[bitVec] ;
  }

} ;

// Alternatively, can use popcnt instruction
class BitCount16 {
  static uint8_t lut_[256] ;

 public:

  static void init_lut() ;
  
  static int getCnt(uint16_t bitVec) {
    return BitCount8::getCnt((uint8_t) (bitVec >> 8)) + 
      BitCount8::getCnt((uint8_t) bitVec) ;
  }
  
} ;


// Fine-grain alignment for up to the first 16 bases.
// Note that this returns the correct edit distance at cells (i,j) for i=j
class FineAlignment16 {

  int queryLen_, endBonus_ ;
  
  int bestAlignedLen_ ;
  int maxScore_ ;

  int mismatchWtP1_ ; // mismatch weight plus one
  int ambigWtP1_ ; // ambig weight plus one
  

 public:

 FineAlignment16(int queryLen, int endBonus, int mismatchWt, int ambigWt)
   : queryLen_(queryLen), endBonus_(endBonus), bestAlignedLen_(0), maxScore_(0),
    mismatchWtP1_(mismatchWt+1), ambigWtP1_(ambigWt+1) {}
  
  void update(uint16_t UP0, uint16_t UP1, int partialLen, int ambigCnt) {
    uint16_t mask = ((uint32_t)(1) << partialLen) - 1 ; // clear higher bits beyond partialLen
    UP0 &= mask ;
    UP1 &= mask ;
    
    int editDist = partialLen + 
      ((partialLen <= 8) ? 
       (BitCount8::getCnt((uint8_t)UP0) - BitCount8::getCnt((uint8_t)UP1)) :
       (BitCount16::getCnt(UP0) - BitCount16::getCnt(UP1))) ;
    
    int bonus = (partialLen == queryLen_) ? endBonus_ : 0 ;
    
    int score = partialLen - editDist * mismatchWtP1_ - ambigCnt * ambigWtP1_ + bonus ;

    if (score > maxScore_) {
      maxScore_ = score ;
      bestAlignedLen_ = partialLen ;
    }
  }
  
  void getBest(int& alignedLen, int& score) {
    alignedLen = bestAlignedLen_ ;
    score = maxScore_ ;
  }

  friend std::ostream& operator<< (std::ostream& os, const FineAlignment16& a) {
    os << "{maxScore = " << a.maxScore_ << ", best col: " << a.bestAlignedLen_ << "}" << endl ;
    return os ;
  }

} ;
