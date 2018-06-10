/* 
The MIT License (MIT)
Copyright (c) 2014 Intel Corp.
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _FAST_EXTEND_DIST_H
#define _FAST_EXTEND_DIST_H

//
// This class handles the distance computations by reading the probes at BitVec, and
// storing the computed distance values in EDVec.
//
template<class BitVec, class EDVec>
class DistVec {

  EDVec dist_, mask_ ;
  BitVec probeShiftCnt_ ;
  int queryLen_, msWordIndex_, probeOffset_ ;
  
 public: 
  
 DistVec(int queryLen):
  queryLen_(queryLen), msWordIndex_(EDVec::getLastWordIndexFor(queryLen)),
    probeOffset_(EDVec::getProbeOffsetFor(queryLen)) {

    probeShiftCnt_.setAll64BitWords(probeOffset_+1) ;
    mask_.setWordsAsMask() ;
  }

  const EDVec& getVec() const {
    return dist_ ;
  }

  void initDist(int initU0, int initU1) {
    dist_.setWordsAsDist(queryLen_, initU0 - initU1) ;
  }

  void addDist (const BitVec& LP0, const BitVec& LP1) {
 
    dist_ += (EDVec(LP0.shiftProbesRight(probeOffset_, probeShiftCnt_)) & mask_) ;
    dist_ -= (EDVec(LP1.shiftProbesRight(probeOffset_, probeShiftCnt_)) & mask_) ;
  }

  EDVec getDeltaDistVec(const BitVec& LP0, const BitVec& LP1) {
    EDVec r = EDVec(LP0.shiftProbesRight(probeOffset_, probeShiftCnt_)) & mask_ ;
    r -= EDVec(LP1.shiftProbesRight(probeOffset_,probeShiftCnt_)) & mask_ ;

    return r ;
  }

  void addDist (const EDVec& deltaDist) {
    dist_ += deltaDist ;
  }

  int getValidDistCnt() const {
    return msWordIndex_ + 1 ;
  }

  int getProbeColumn(int index) const {
    return EDVec::getDistAt(index, queryLen_) - 1 ;
  }

  int getDist(int index) const {
    return (int) dist_.getWord(index) ;
  }

  int getTotalDist() const {
    return (int) dist_.getWord(msWordIndex_) ;
  }

  void setMin (const DistVec& other) {
    this->dist_.setMin(other.dist_) ;
  }

  friend std::ostream& operator<< (std::ostream& os, const DistVec& d) {
    os << d.dist_ ;
    return os ;
  }

} ;


#endif // #ifndef _FAST_EXTEND_DIST_H
