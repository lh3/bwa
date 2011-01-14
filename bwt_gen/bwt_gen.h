/*

   BWTConstruct.h		BWT-Index Construction

   This module constructs BWT and auxiliary data structures.

   Copyright (C) 2004, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#ifndef BWT_GEN_H
#define BWT_GEN_H

#define ALPHABET_SIZE				4
#define BIT_PER_CHAR				2
#define CHAR_PER_WORD				16
#define CHAR_PER_BYTE				4

#define BITS_IN_WORD 32
#define BITS_IN_BYTE 8
#define BYTES_IN_WORD 4

#define ALL_ONE_MASK 0xFFFFFFFF
#define DNA_OCC_CNT_TABLE_SIZE_IN_WORD	65536

#define BITS_PER_OCC_VALUE			16
#define OCC_VALUE_PER_WORD			2
#define OCC_INTERVAL				256
#define OCC_INTERVAL_MAJOR			65536

#define TRUE    1
#define FALSE   0

#define BWTINC_INSERT_SORT_NUM_ITEM 7

#define average(value1, value2)					( ((value1) & (value2)) + ((value1) ^ (value2)) / 2 )
#define min(value1, value2)						( ((value1) < (value2)) ? (value1) : (value2) )
#define max(value1, value2)						( ((value1) > (value2)) ? (value1) : (value2) )
#define med3(a, b, c)							( a<b ? (b<c ? b : a<c ? c : a) : (b>c ? b : a>c ? c : a))
#define swap(a, b, t);							t = a; a = b; b = t;
#define truncateLeft(value, offset)				( (value) << (offset) >> (offset) )
#define truncateRight(value, offset)			( (value) >> (offset) << (offset) )
#define DNA_OCC_SUM_EXCEPTION(sum)			((sum & 0xfefefeff) == 0)

typedef struct SaIndexRange {
	unsigned int startSaIndex;
	unsigned int endSaIndex;
} SaIndexRange;

typedef struct BWT {
	unsigned int textLength;			// length of the text
	unsigned int saInterval;			// interval between two SA values stored explicitly
	unsigned int inverseSaInterval;		// interval between two inverse SA stored explicitly
	unsigned int inverseSa0;			// SA-1[0]
	unsigned int *cumulativeFreq;		// cumulative frequency
	unsigned int *bwtCode;				// BWT code
	unsigned int *occValue;				// Occurrence values stored explicitly
	unsigned int *occValueMajor;		// Occurrence values stored explicitly
	unsigned int *saValue;				// SA values stored explicitly
	unsigned int *inverseSa;			// Inverse SA stored explicitly
	SaIndexRange *saIndexRange;			// SA index range
	int saIndexRangeNumOfChar;			// Number of characters indexed in SA index range
	unsigned int *saValueOnBoundary;	// Pre-calculated frequently referred data
	unsigned int *decodeTable;			// For decoding BWT by table lookup
	unsigned int decodeTableGenerated;	// == TRUE if decode table is generated on load and will be freed
	unsigned int bwtSizeInWord;			// Temporary variable to hold the memory allocated
	unsigned int occSizeInWord;			// Temporary variable to hold the memory allocated
	unsigned int occMajorSizeInWord;	// Temporary variable to hold the memory allocated
	unsigned int saValueSize;			// Temporary variable to hold the memory allocated
	unsigned int inverseSaSize;			// Temporary variable to hold the memory allocated
	unsigned int saIndexRangeSize;		// Temporary variable to hold the memory allocated
} BWT;

typedef struct BWTInc {
	BWT *bwt;
	unsigned int numberOfIterationDone;
	unsigned int *cumulativeCountInCurrentBuild;
	unsigned int availableWord;
	unsigned int targetTextLength;
	float targetNBit;
	unsigned int buildSize;
	unsigned int initialMaxBuildSize;
	unsigned int incMaxBuildSize;
	unsigned int firstCharInLastIteration;
	unsigned int *workingMemory;
	unsigned int *packedText;
	unsigned char *textBuffer;
	unsigned int *packedShift;
} BWTInc;

#endif
