/* The MIT License

   Copyright (c) 2003-2006, 2008, 2009, by Heng Li <lh3lh3@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "stdaln.h"

/* char -> 17 (=16+1) nucleotides */
unsigned char aln_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,16 /*'-'*/,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15, 1,14, 4, 11,15,15, 2, 13,15,15,10, 15, 5,15,15,
	15,15, 3, 6,  8,15, 7, 9,  0,12,15,15, 15,15,15,15,
	15, 1,14, 4, 11,15,15, 2, 13,15,15,10, 15, 5,15,15,
	15,15, 3, 6,  8,15, 7, 9,  0,12,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};
char *aln_nt16_rev_table = "XAGRCMSVTWKDYHBN-";

/* char -> 5 (=4+1) nucleotides */
unsigned char aln_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 2,  4, 4, 4, 1,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 2,  4, 4, 4, 1,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};
char *aln_nt4_rev_table = "AGCTN-";

/* char -> 22 (=20+1+1) amino acids */
unsigned char aln_aa_table[256] = {
	21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
	21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
	21,21,21,21, 21,21,21,21, 21,21,20,21, 21,22 /*'-'*/,21,21,
	21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
	21, 0,21, 4,  3, 6,13, 7,  8, 9,21,11, 10,12, 2,21,
	14, 5, 1,15, 16,21,19,17, 21,18,21,21, 21,21,21,21,
	21, 0,21, 4,  3, 6,13, 7,  8, 9,21,11, 10,12, 2,21,
	14, 5, 1,15, 16,21,19,17, 21,18,21,21, 21,21,21,21,
	21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
	21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
	21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
	21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
	21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
	21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
	21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21,
	21,21,21,21, 21,21,21,21, 21,21,21,21, 21,21,21,21
};
char *aln_aa_rev_table = "ARNDCQEGHILKMFPSTWYV*X-";
                       /* 01234567890123456789012 */

/* translation table. They are useless in stdaln.c, but when you realize you need it, you need not write the table again. */
unsigned char aln_trans_table_eu[66] = {
	11,11, 2, 2,  1, 1,15,15, 16,16,16,16,  9,12, 9, 9,
	 6, 6, 3, 3,  7, 7, 7, 7,  0, 0, 0, 0, 19,19,19,19,
	 5, 5, 8, 8,  1, 1, 1, 1, 14,14,14,14, 10,10,10,10,
	20,20,18,18, 20,17, 4, 4, 15,15,15,15, 10,10,13,13, 21, 22
};
char *aln_trans_table_eu_char = "KKNNRRSSTTTTIMIIEEDDGGGGAAAAVVVVQQHHRRRRPPPPLLLL**YY*WCCSSSSLLFFX";
                              /* 01234567890123456789012345678901234567890123456789012345678901234 */
int aln_sm_blosum62[] = {
/*	 A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  *  X */
	 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-4, 0,
	-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-4,-1,
	-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,-4,-1,
	-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,-4,-1,
	 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-4,-2,
	-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,-4,-1,
	-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,-4,-1,
	 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-4,-1,
	-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,-4,-1,
	-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-4,-1,
	-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-1,
	-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,-4,-1,
	-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-4,-1,
	-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-4,-1,
	-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-4,-2,
	 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2,-4, 0,
	 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-4, 0,
	-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-2,
	-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-4,-1,
	 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-4,-1,
	-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1,-4,
	 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-4,-1
};

int aln_sm_blosum45[] = {
/*	 A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  *  X */
	 5,-2,-1,-2,-1,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-2,-2, 0,-5, 0,
	-2, 7, 0,-1,-3, 1, 0,-2, 0,-3,-2, 3,-1,-2,-2,-1,-1,-2,-1,-2,-5,-1,
	-1, 0, 6, 2,-2, 0, 0, 0, 1,-2,-3, 0,-2,-2,-2, 1, 0,-4,-2,-3,-5,-1,
	-2,-1, 2, 7,-3, 0, 2,-1, 0,-4,-3, 0,-3,-4,-1, 0,-1,-4,-2,-3,-5,-1,
	-1,-3,-2,-3,12,-3,-3,-3,-3,-3,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1,-5,-2,
	-1, 1, 0, 0,-3, 6, 2,-2, 1,-2,-2, 1, 0,-4,-1, 0,-1,-2,-1,-3,-5,-1,
	-1, 0, 0, 2,-3, 2, 6,-2, 0,-3,-2, 1,-2,-3, 0, 0,-1,-3,-2,-3,-5,-1,
	 0,-2, 0,-1,-3,-2,-2, 7,-2,-4,-3,-2,-2,-3,-2, 0,-2,-2,-3,-3,-5,-1,
	-2, 0, 1, 0,-3, 1, 0,-2,10,-3,-2,-1, 0,-2,-2,-1,-2,-3, 2,-3,-5,-1,
	-1,-3,-2,-4,-3,-2,-3,-4,-3, 5, 2,-3, 2, 0,-2,-2,-1,-2, 0, 3,-5,-1,
	-1,-2,-3,-3,-2,-2,-2,-3,-2, 2, 5,-3, 2, 1,-3,-3,-1,-2, 0, 1,-5,-1,
	-1, 3, 0, 0,-3, 1, 1,-2,-1,-3,-3, 5,-1,-3,-1,-1,-1,-2,-1,-2,-5,-1,
	-1,-1,-2,-3,-2, 0,-2,-2, 0, 2, 2,-1, 6, 0,-2,-2,-1,-2, 0, 1,-5,-1,
	-2,-2,-2,-4,-2,-4,-3,-3,-2, 0, 1,-3, 0, 8,-3,-2,-1, 1, 3, 0,-5,-1,
	-1,-2,-2,-1,-4,-1, 0,-2,-2,-2,-3,-1,-2,-3, 9,-1,-1,-3,-3,-3,-5,-1,
	 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-3,-1,-2,-2,-1, 4, 2,-4,-2,-1,-5, 0,
	 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-1,-1, 2, 5,-3,-1, 0,-5, 0,
	-2,-2,-4,-4,-5,-2,-3,-2,-3,-2,-2,-2,-2, 1,-3,-4,-3,15, 3,-3,-5,-2,
	-2,-1,-2,-2,-3,-1,-2,-3, 2, 0, 0,-1, 0, 3,-3,-2,-1, 3, 8,-1,-5,-1,
	 0,-2,-3,-3,-1,-3,-3,-3,-3, 3, 1,-2, 1, 0,-3,-1, 0,-3,-1, 5,-5,-1,
	-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5, 1,-5,
	 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,-2,-1,-1,-5,-1
};

int aln_sm_nt[] = {
/*	 X  A  G  R  C  M  S  V  T  W  K  D  Y  H  B  N */
	-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
	-2, 2,-1, 1,-2, 1,-2, 0,-2, 1,-2, 0,-2, 0,-2, 0,
	-2,-1, 2, 1,-2,-2, 1, 0,-2,-2, 1, 0,-2,-2, 0, 0,
	-2, 1, 1, 1,-2,-1,-1, 0,-2,-1,-1, 0,-2, 0, 0, 0,
	-2,-2,-2,-2, 2, 1, 1, 0,-1,-2,-2,-2, 1, 0, 0, 0,
	-2, 1,-2,-1, 1, 1,-1, 0,-2,-1,-2, 0,-1, 0, 0, 0,
	-2,-2, 1,-1, 1,-1, 1, 0,-2,-2,-1, 0,-1, 0, 0, 0,
	-2, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0,
	-2,-2,-2,-2,-1,-2,-2,-2, 2, 1, 1, 0, 1, 0, 0, 0,
	-2, 1,-2,-1,-2,-1,-2, 0, 1, 1,-1, 0,-1, 0, 0, 0,
	-2,-2, 1,-1,-2,-2,-1, 0, 1,-1, 1, 0,-1, 0, 0, 0,
	-2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	-2,-2,-2,-2, 1,-1,-1, 0, 1,-1,-1, 0, 1, 0, 0, 0,
	-2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

int aln_sm_read[] = {
/*	  X   A   G   R   C   M   S   V   T   W   K   D   Y   H   B   N  */
	-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,-17,
	-17,  2,-17,  1,-17,  1,-17,  0,-17,  1,-17,  0,-17,  0,-17,  0,
	-17,-17,  2,  1,-17,-17,  1,  0,-17,-17,  1,  0,-17,-17,  0,  0,
	-17,  1,  1,  1,-17,-17,-17,  0,-17,-17,-17,  0,-17,  0,  0,  0,
	-17,-17,-17,-17,  2,  1,  1,  0,-17,-17,-17,-17,  1,  0,  0,  0,
	-17,  1,-17,-17,  1,  1,-17,  0,-17,-17,-17,  0,-17,  0,  0,  0,
	-17,-17,  1,-17,  1,-17,  1,  0,-17,-17,-17,  0,-17,  0,  0,  0,
	-17,  0,  0,  0,  0,  0,  0,  0,-17,  0,  0,  0,  0,  0,  0,  0,
	-17,-17,-17,-17,-17,-17,-17,-17,  2,  1,  1,  0,  1,  0,  0,  0,
	-17,  1,-17,-17,-17,-17,-17,  0,  1,  1,-17,  0,-17,  0,  0,  0,
	-17,-17,  1,-17,-17,-17,-17,  0,  1,-17,  1,  0,-17,  0,  0,  0,
	-17,  0,  0,  0,-17,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	-17,-17,-17,-17,  1,-17,-17,  0,  1,-17,-17,  0,  1,  0,  0,  0,
	-17,  0,-17,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	-17,-17,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
	-17,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
};

int aln_sm_hs[] = {
/*     A    G    C    T    N */
	  91, -31,-114,-123, -44,
	 -31, 100,-125,-114, -42,
	-123,-125, 100, -31, -42,
	-114,-114, -31,  91, -42,
	 -44, -42, -42, -42, -43
};

int aln_sm_maq[] = {
	11, -19, -19, -19, -13,
	-19, 11, -19, -19, -13,
	-19, -19, 11, -19, -13,
	-19, -19, -19, 11, -13,
	-13, -13, -13, -13, -13
};

int aln_sm_blast[] = {
	1, -3, -3, -3, -2,
	-3, 1, -3, -3, -2,
	-3, -3, 1, -3, -2,
	-3, -3, -3, 1, -2,
	-2, -2, -2, -2, -2
};

/********************/
/* START OF align.c */
/********************/

AlnParam aln_param_blast   = {  5,  2,  2, aln_sm_blast, 5, 50 };
AlnParam aln_param_bwa     = { 26,  9,  5, aln_sm_maq, 5, 50 };
AlnParam aln_param_nt2nt   = {  8,  2,  2, aln_sm_nt, 16, 75 };
AlnParam aln_param_rd2rd   = {  1, 19, 19, aln_sm_read, 16, 75 };
AlnParam aln_param_aa2aa   = { 10,  2,  2, aln_sm_blosum62, 22, 50 };

AlnAln *aln_init_AlnAln()
{
	AlnAln *aa;
	aa = (AlnAln*)malloc(sizeof(AlnAln));
	aa->path = 0;
	aa->out1 = aa->out2 = aa->outm = 0;
	aa->path_len = 0;
	return aa;
}
void aln_free_AlnAln(AlnAln *aa)
{
	free(aa->path); free(aa->cigar32);
	free(aa->out1); free(aa->out2); free(aa->outm);
	free(aa);
}

/***************************/
/* START OF common_align.c */
/***************************/

#define LOCAL_OVERFLOW_THRESHOLD 32000
#define LOCAL_OVERFLOW_REDUCE 16000
#define NT_LOCAL_SCORE int
#define NT_LOCAL_SHIFT 16
#define NT_LOCAL_MASK 0xffff

#define SET_INF(s) (s).M = (s).I = (s).D = MINOR_INF;

#define set_M(MM, cur, p, sc)							\
{														\
	if ((p)->M >= (p)->I) {								\
		if ((p)->M >= (p)->D) {							\
			(MM) = (p)->M + (sc); (cur)->Mt = FROM_M;	\
		} else {										\
			(MM) = (p)->D + (sc); (cur)->Mt = FROM_D;	\
		}												\
	} else {											\
		if ((p)->I > (p)->D) {							\
			(MM) = (p)->I + (sc); (cur)->Mt = FROM_I;	\
		} else {										\
			(MM) = (p)->D + (sc); (cur)->Mt = FROM_D;	\
		}												\
	}													\
}
#define set_I(II, cur, p)								\
{														\
	if ((p)->M - gap_open > (p)->I) {					\
		(cur)->It = FROM_M;								\
		(II) = (p)->M - gap_open - gap_ext;				\
	} else {											\
		(cur)->It = FROM_I;								\
		(II) = (p)->I - gap_ext;						\
	}													\
}
#define set_end_I(II, cur, p)							\
{														\
	if (gap_end >= 0) {									\
		if ((p)->M - gap_open > (p)->I) {				\
			(cur)->It = FROM_M;							\
			(II) = (p)->M - gap_open - gap_end;			\
		} else {										\
			(cur)->It = FROM_I;							\
			(II) = (p)->I - gap_end;					\
		}												\
	} else set_I(II, cur, p);							\
}
#define set_D(DD, cur, p)								\
{														\
	if ((p)->M - gap_open > (p)->D) {					\
		(cur)->Dt = FROM_M;								\
		(DD) = (p)->M - gap_open - gap_ext;				\
	} else {											\
		(cur)->Dt = FROM_D;								\
		(DD) = (p)->D - gap_ext;						\
	}													\
}
#define set_end_D(DD, cur, p)							\
{														\
	if (gap_end >= 0) {									\
		if ((p)->M - gap_open > (p)->D) {				\
			(cur)->Dt = FROM_M;							\
			(DD) = (p)->M - gap_open - gap_end;			\
		} else {										\
			(cur)->Dt = FROM_D;							\
			(DD) = (p)->D - gap_end;					\
		}												\
	} else set_D(DD, cur, p);							\
}

typedef struct
{
	unsigned char Mt:3, It:2, Dt:2;
} dpcell_t;

typedef struct
{
	int M, I, D;
} dpscore_t;

/* build score profile for accelerating alignment, in theory */
void aln_init_score_array(unsigned char *seq, int len, int row, int *score_matrix, int **s_array)
{
	int *tmp, *tmp2, i, k;
	for (i = 0; i != row; ++i) {
		tmp = score_matrix + i * row;
		tmp2 = s_array[i];
		for (k = 0; k != len; ++k)
			tmp2[k] = tmp[seq[k]];
	}
}
/***************************
 * banded global alignment *
 ***************************/
int aln_global_core(unsigned char *seq1, int len1, unsigned char *seq2, int len2, const AlnParam *ap,
					path_t *path, int *path_len)
{
	register int i, j;
	dpcell_t **dpcell, *q;
	dpscore_t *curr, *last, *s;
	path_t *p;
	int b1, b2, tmp_end;
	int *mat, end, max;
	unsigned char type, ctype;

	int gap_open, gap_ext, gap_end, b;
	int *score_matrix, N_MATRIX_ROW;

	/* initialize some align-related parameters. just for compatibility */
	gap_open = ap->gap_open;
	gap_ext = ap->gap_ext;
	gap_end = ap->gap_end;
	b = ap->band_width;
	score_matrix = ap->matrix;
	N_MATRIX_ROW = ap->row;
	
	if (len1 == 0 || len2 == 0) {
		*path_len = 0;
		return 0;
	}
	/* calculate b1 and b2 */
	if (len1 > len2) {
		b1 = len1 - len2 + b;
		b2 = b;
	} else {
		b1 = b;
		b2 = len2 - len1 + b;
	}
	if (b1 > len1) b1 = len1;
	if (b2 > len2) b2 = len2;
	--seq1; --seq2;

	/* allocate memory */
	end = (b1 + b2 <= len1)? (b1 + b2 + 1) : (len1 + 1);
	dpcell = (dpcell_t**)malloc(sizeof(dpcell_t*) * (len2 + 1));
	for (j = 0; j <= len2; ++j)
		dpcell[j] = (dpcell_t*)malloc(sizeof(dpcell_t) * end);
	for (j = b2 + 1; j <= len2; ++j)
		dpcell[j] -= j - b2;
	curr = (dpscore_t*)malloc(sizeof(dpscore_t) * (len1 + 1));
	last = (dpscore_t*)malloc(sizeof(dpscore_t) * (len1 + 1));
	
	/* set first row */
	SET_INF(*curr); curr->M = 0;
	for (i = 1, s = curr + 1; i < b1; ++i, ++s) {
		SET_INF(*s);
		set_end_D(s->D, dpcell[0] + i, s - 1);
	}
	s = curr; curr = last; last = s;

	/* core dynamic programming, part 1 */
	tmp_end = (b2 < len2)? b2 : len2 - 1;
	for (j = 1; j <= tmp_end; ++j) {
		q = dpcell[j]; s = curr; SET_INF(*s);
		set_end_I(s->I, q, last);
		end = (j + b1 <= len1 + 1)? (j + b1 - 1) : len1;
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		++s; ++q;
		for (i = 1; i != end; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]); /* this will change s->M ! */
			set_I(s->I, q, last + i);
			set_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + i - 1, mat[seq1[i]]);
		set_D(s->D, q, s - 1);
		if (j + b1 - 1 > len1) { /* bug fixed, 040227 */
			set_end_I(s->I, q, last + i);
		} else s->I = MINOR_INF;
		s = curr; curr = last; last = s;
	}
	/* last row for part 1, use set_end_D() instead of set_D() */
	if (j == len2 && b2 != len2 - 1) {
		q = dpcell[j]; s = curr; SET_INF(*s);
		set_end_I(s->I, q, last);
		end = (j + b1 <= len1 + 1)? (j + b1 - 1) : len1;
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		++s; ++q;
		for (i = 1; i != end; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]); /* this will change s->M ! */
			set_I(s->I, q, last + i);
			set_end_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + i - 1, mat[seq1[i]]);
		set_end_D(s->D, q, s - 1);
		if (j + b1 - 1 > len1) { /* bug fixed, 040227 */
			set_end_I(s->I, q, last + i);
		} else s->I = MINOR_INF;
		s = curr; curr = last; last = s;
		++j;
	}

	/* core dynamic programming, part 2 */
	for (; j <= len2 - b2 + 1; ++j) {
		SET_INF(curr[j - b2]);
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		end = j + b1 - 1;
		for (i = j - b2 + 1, q = dpcell[j] + i, s = curr + i; i != end; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]);
			set_I(s->I, q, last + i);
			set_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + i - 1, mat[seq1[i]]);
		set_D(s->D, q, s - 1);
		s->I = MINOR_INF;
		s = curr; curr = last; last = s;
	}

	/* core dynamic programming, part 3 */
	for (; j < len2; ++j) {
		SET_INF(curr[j - b2]);
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		for (i = j - b2 + 1, q = dpcell[j] + i, s = curr + i; i < len1; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]);
			set_I(s->I, q, last + i);
			set_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + len1 - 1, mat[seq1[i]]);
		set_end_I(s->I, q, last + i);
		set_D(s->D, q, s - 1);
		s = curr; curr = last; last = s;
	}
	/* last row */
	if (j == len2) {
		SET_INF(curr[j - b2]);
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		for (i = j - b2 + 1, q = dpcell[j] + i, s = curr + i; i < len1; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]);
			set_I(s->I, q, last + i);
			set_end_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + len1 - 1, mat[seq1[i]]);
		set_end_I(s->I, q, last + i);
		set_end_D(s->D, q, s - 1);
		s = curr; curr = last; last = s;
	}

	/* backtrace */
	i = len1; j = len2;
	q = dpcell[j] + i;
	s = last + len1;
	max = s->M; type = q->Mt; ctype = FROM_M;
	if (s->I > max) { max = s->I; type = q->It; ctype = FROM_I; }
	if (s->D > max) { max = s->D; type = q->Dt; ctype = FROM_D; }

	p = path;
	p->ctype = ctype; p->i = i; p->j = j; /* bug fixed 040408 */
	++p;
	do {
		switch (ctype) {
			case FROM_M: --i; --j; break;
			case FROM_I: --j; break;
			case FROM_D: --i; break;
		}
		q = dpcell[j] + i;
		ctype = type;
		switch (type) {
			case FROM_M: type = q->Mt; break;
			case FROM_I: type = q->It; break;
			case FROM_D: type = q->Dt; break;
		}
		p->ctype = ctype; p->i = i; p->j = j;
		++p;
	} while (i || j);
	*path_len = p - path - 1;

	/* free memory */
	for (j = b2 + 1; j <= len2; ++j)
		dpcell[j] += j - b2;
	for (j = 0; j <= len2; ++j)
		free(dpcell[j]);
	free(dpcell);
	free(curr); free(last);
	
	return max;
}
/*************************************************
 * local alignment combined with banded strategy *
 *************************************************/
int aln_local_core(unsigned char *seq1, int len1, unsigned char *seq2, int len2, const AlnParam *ap,
				   path_t *path, int *path_len, int _thres, int *_subo)
{
	register NT_LOCAL_SCORE *s;
	register int i;
	int q, r, qr, tmp_len, qr_shift;
	int **s_array, *score_array;
	int e, f;
	int is_overflow, of_base;
	NT_LOCAL_SCORE *eh, curr_h, last_h, curr_last_h;
	int j, start_i, start_j, end_i, end_j;
	path_t *p;
	int score_f, score_r, score_g;
	int start, end, max_score;
	int thres, *suba, *ss;

	int gap_open, gap_ext, b;
	int *score_matrix, N_MATRIX_ROW;

	/* initialize some align-related parameters. just for compatibility */
	gap_open = ap->gap_open;
	gap_ext = ap->gap_ext;
	b = ap->band_width;
	score_matrix = ap->matrix;
	N_MATRIX_ROW = ap->row;
	thres = _thres > 0? _thres : -_thres;

	if (len1 == 0 || len2 == 0) return -1;

	/* allocate memory */
	suba = (int*)malloc(sizeof(int) * (len2 + 1));
	eh = (NT_LOCAL_SCORE*)malloc(sizeof(NT_LOCAL_SCORE) * (len1 + 1));
	s_array = (int**)malloc(sizeof(int*) * N_MATRIX_ROW);
	for (i = 0; i != N_MATRIX_ROW; ++i)
		s_array[i] = (int*)malloc(sizeof(int) * len1);
	/* initialization */
	aln_init_score_array(seq1, len1, N_MATRIX_ROW, score_matrix, s_array);
	q = gap_open;
	r = gap_ext;
	qr = q + r;
	qr_shift = (qr+1) << NT_LOCAL_SHIFT;
	tmp_len = len1 + 1;
	start_i = start_j = end_i = end_j = 0;
	for (i = 0, max_score = 0; i != N_MATRIX_ROW * N_MATRIX_ROW; ++i)
		if (max_score < score_matrix[i]) max_score = score_matrix[i];
	/* convert the coordinate */
	--seq1; --seq2;
	for (i = 0; i != N_MATRIX_ROW; ++i) --s_array[i];

	/* forward dynamic programming */
	for (i = 0, s = eh; i != tmp_len; ++i, ++s) *s = 0;
	score_f = 0;
	is_overflow = of_base = 0;
	suba[0] = 0;
	for (j = 1, ss = suba + 1; j <= len2; ++j, ++ss) {
		int subo = 0;
		last_h = f = 0;
		score_array = s_array[seq2[j]];
		if (is_overflow) { /* adjust eh[] array if overflow occurs. */
			/* If LOCAL_OVERFLOW_REDUCE is too small, optimal alignment might be missed.
			 * If it is too large, this block will be excuted frequently and therefore
			 * slow down the whole program.
			 * Acually, smaller LOCAL_OVERFLOW_REDUCE might also help to reduce the
			 * number of assignments because it sets some cells to zero when overflow
			 * happens. */
			int tmp, tmp2;
			score_f -= LOCAL_OVERFLOW_REDUCE;
			of_base += LOCAL_OVERFLOW_REDUCE;
			is_overflow = 0;
			for (i = 1, s = eh; i <= tmp_len; ++i, ++s) {
				tmp = *s >> NT_LOCAL_SHIFT; tmp2 = *s & NT_LOCAL_MASK;
				if (tmp2 < LOCAL_OVERFLOW_REDUCE) tmp2 = 0;
				else tmp2 -= LOCAL_OVERFLOW_REDUCE;
				if (tmp < LOCAL_OVERFLOW_REDUCE) tmp = 0;
				else tmp -= LOCAL_OVERFLOW_REDUCE;
				*s = (tmp << NT_LOCAL_SHIFT) | tmp2;
			}
		}
		for (i = 1, s = eh; i != tmp_len; ++i, ++s) {
			/* prepare for calculate current h */
			curr_h = (*s >> NT_LOCAL_SHIFT) + score_array[i];
			if (curr_h < 0) curr_h = 0;
			if (last_h > 0) { /* initialize f */
				f = (f > last_h - q)? f - r : last_h - qr;
				if (curr_h < f) curr_h = f;
			}
			if (*(s+1) >= qr_shift) { /* initialize e */
				curr_last_h = *(s+1) >> NT_LOCAL_SHIFT;
				e = ((*s & NT_LOCAL_MASK) > curr_last_h - q)? (*s & NT_LOCAL_MASK) - r : curr_last_h - qr;
				if (curr_h < e) curr_h = e;
				*s = (last_h << NT_LOCAL_SHIFT) | e;
			} else *s = last_h << NT_LOCAL_SHIFT; /* e = 0 */
			last_h = curr_h;
			if (subo < curr_h) subo = curr_h;
			if (score_f < curr_h) {
				score_f = curr_h; end_i = i; end_j = j;
				if (score_f > LOCAL_OVERFLOW_THRESHOLD) is_overflow = 1;
			}
		}
		*s = last_h << NT_LOCAL_SHIFT;
		*ss = subo + of_base;
	}
	score_f += of_base;

	if (score_f < thres) { /* no matching residue at all, 090218 */
		if (path_len) *path_len = 0;
		goto end_func;
	}
	if (path == 0) goto end_func; /* skip path-filling */

	/* reverse dynamic programming */
	for (i = end_i, s = eh + end_i; i >= 0; --i, --s) *s = 0;
	if (end_i == 0 || end_j == 0) goto end_func; /* no local match */
	score_r = score_matrix[seq1[end_i] * N_MATRIX_ROW + seq2[end_j]];
	is_overflow = of_base = 0;
	start_i = end_i; start_j = end_j;
	eh[end_i] = ((NT_LOCAL_SCORE)(qr + score_r)) << NT_LOCAL_SHIFT; /* in order to initialize f and e, 040408 */
	start = end_i - 1;
	end = end_i - 3;
	if (end <= 0) end = 0;

	/* second pass DP can be done in a band, speed will thus be enhanced */
	for (j = end_j - 1; j != 0; --j) {
		last_h = f = 0;
		score_array = s_array[seq2[j]];
		if (is_overflow) { /* adjust eh[] array if overflow occurs. */
			int tmp, tmp2;
			score_r -= LOCAL_OVERFLOW_REDUCE;
			of_base += LOCAL_OVERFLOW_REDUCE;
			is_overflow = 0;
			for (i = start, s = eh + start + 1; i >= end; --i, --s) {
				tmp = *s >> NT_LOCAL_SHIFT; tmp2 = *s & NT_LOCAL_MASK;
				if (tmp2 < LOCAL_OVERFLOW_REDUCE) tmp2 = 0;
				else tmp2 -= LOCAL_OVERFLOW_REDUCE;
				if (tmp < LOCAL_OVERFLOW_REDUCE) tmp = 0;
				else tmp -= LOCAL_OVERFLOW_REDUCE;
				*s = (tmp << NT_LOCAL_SHIFT) | tmp2;
			}
		}
		for (i = start, s = eh + start + 1; i != end; --i, --s) {
			/* prepare for calculate current h */
			curr_h = (*s >> NT_LOCAL_SHIFT) + score_array[i];
			if (curr_h < 0) curr_h = 0;
			if (last_h > 0) { /* initialize f */
				f = (f > last_h - q)? f - r : last_h - qr;
				if (curr_h < f) curr_h = f;
			}
			curr_last_h = *(s-1) >> NT_LOCAL_SHIFT;
			e = ((*s & NT_LOCAL_MASK) > curr_last_h - q)? (*s & NT_LOCAL_MASK) - r : curr_last_h - qr;
			if (e < 0) e = 0;
			if (curr_h < e) curr_h = e;
			*s = (last_h << NT_LOCAL_SHIFT) | e;
			last_h = curr_h;
			if (score_r < curr_h) {
				score_r = curr_h; start_i = i; start_j = j;
				if (score_r + of_base - qr == score_f) {
					j = 1; break;
				}
				if (score_r > LOCAL_OVERFLOW_THRESHOLD) is_overflow = 1;
			}
		}
		*s = last_h << NT_LOCAL_SHIFT;
		/* recalculate start and end, the boundaries of the band */
		if ((eh[start] >> NT_LOCAL_SHIFT) <= qr) --start;
		if (start <= 0) start = 0;
		end = start_i - (start_j - j) - (score_r + of_base + (start_j - j) * max_score) / r - 1;
		if (end <= 0) end = 0;
	}

	if (_subo) {
		int tmp2 = 0, tmp = (int)(start_j - .33 * (end_j - start_j) + .499);
		for (j = 1; j <= tmp; ++j)
			if (tmp2 < suba[j]) tmp2 = suba[j];
		tmp = (int)(end_j + .33 * (end_j - start_j) + .499);
		for (j = tmp; j <= len2; ++j)
			if (tmp2 < suba[j]) tmp2 = suba[j];
		*_subo = tmp2;
	}

	if (path_len == 0) {
		path[0].i = start_i; path[0].j = start_j;
		path[1].i = end_i; path[1].j = end_j;
		goto end_func;
	}

	score_r += of_base;
	score_r -= qr;

#ifdef DEBUG
	/* this seems not a bug */
	if (score_f != score_r)
		fprintf(stderr, "[aln_local_core] unknown flaw occurs: score_f(%d) != score_r(%d)\n", score_f, score_r);
#endif

	if (_thres > 0) { /* call global alignment to fill the path */
		score_g = 0;
		j = (end_i - start_i > end_j - start_j)? end_i - start_i : end_j - start_j;
		++j; /* j is the maximum band_width */
		for (i = ap->band_width;; i <<= 1) {
			AlnParam ap_real = *ap;
			ap_real.gap_end = -1;
			ap_real.band_width = i;
			score_g = aln_global_core(seq1 + start_i, end_i - start_i + 1, seq2 + start_j,
									  end_j - start_j + 1, &ap_real, path, path_len);
			if (score_g == score_r || score_f == score_g) break;
			if (i > j) break;
		}
		if (score_r > score_g && score_f > score_g) {
			fprintf(stderr, "[aln_local_core] Potential bug: (%d,%d) > %d\n", score_f, score_r, score_g);
			score_f = score_r = -1;
		} else score_f = score_g;

		/* convert coordinate */
		for (p = path + *path_len - 1; p >= path; --p) {
			p->i += start_i - 1;
			p->j += start_j - 1;
		}
	} else { /* just store the start and end */
		*path_len = 2;
		path[1].i = start_i; path[1].j = start_j;
		path->i = end_i; path->j = end_j;
	}

end_func:
	/* free */
	free(eh); free(suba);
	for (i = 0; i != N_MATRIX_ROW; ++i) {
		++s_array[i];
		free(s_array[i]);
	}
	free(s_array);
	return score_f;
}
AlnAln *aln_stdaln_aux(const char *seq1, const char *seq2, const AlnParam *ap,
					   int type, int thres, int len1, int len2)
{
	unsigned char *seq11, *seq22;
	int score;
	int i, j, l;
	path_t *p;
	char *out1, *out2, *outm;
	AlnAln *aa;

	if (len1 < 0) len1 = strlen(seq1);
	if (len2 < 0) len2 = strlen(seq2);

	aa = aln_init_AlnAln();
	seq11 = (unsigned char*)malloc(sizeof(unsigned char) * len1);
	seq22 = (unsigned char*)malloc(sizeof(unsigned char) * len2);
	aa->path = (path_t*)malloc(sizeof(path_t) * (len1 + len2 + 1));

	if (ap->row < 10) { /* 4-nucleotide alignment */
		for (i = 0; i < len1; ++i)
			seq11[i] = aln_nt4_table[(int)seq1[i]];
		for (j = 0; j < len2; ++j)
			seq22[j] = aln_nt4_table[(int)seq2[j]];
	} else if (ap->row < 20) { /* 16-nucleotide alignment */
		for (i = 0; i < len1; ++i)
			seq11[i] = aln_nt16_table[(int)seq1[i]];
		for (j = 0; j < len2; ++j)
			seq22[j] = aln_nt16_table[(int)seq2[j]];
	} else { /* amino acids */
		for (i = 0; i < len1; ++i)
			seq11[i] = aln_aa_table[(int)seq1[i]];
		for (j = 0; j < len2; ++j)
			seq22[j] = aln_aa_table[(int)seq2[j]];
	}
	
	if (type == ALN_TYPE_GLOBAL) score = aln_global_core(seq11, len1, seq22, len2, ap, aa->path, &aa->path_len);
	else if (type == ALN_TYPE_LOCAL) score = aln_local_core(seq11, len1, seq22, len2, ap, aa->path, &aa->path_len, thres, &aa->subo);
	else if (type == ALN_TYPE_EXTEND)  score = aln_extend_core(seq11, len1, seq22, len2, ap, aa->path, &aa->path_len, 1, 0);
	else {
		free(seq11); free(seq22); free(aa->path);
		aln_free_AlnAln(aa);
		return 0;
	}
	aa->score = score;

	if (thres > 0) {
		out1 = aa->out1 = (char*)malloc(sizeof(char) * (aa->path_len + 1));
		out2 = aa->out2 = (char*)malloc(sizeof(char) * (aa->path_len + 1));
		outm = aa->outm = (char*)malloc(sizeof(char) * (aa->path_len + 1));

		--seq1; --seq2;
		--seq11; --seq22;

		p = aa->path + aa->path_len - 1;

		for (l = 0; p >= aa->path; --p, ++l) {
			switch (p->ctype) {
			case FROM_M: out1[l] = seq1[p->i]; out2[l] = seq2[p->j];
				outm[l] = (seq11[p->i] == seq22[p->j] && seq11[p->i] != ap->row)? '|' : ' ';
				break;
			case FROM_I: out1[l] = '-'; out2[l] = seq2[p->j]; outm[l] = ' '; break;
			case FROM_D: out1[l] = seq1[p->i]; out2[l] = '-'; outm[l] = ' '; break;
			}
		}
		out1[l] = out2[l] = outm[l] = '\0';
		++seq11; ++seq22;
	}

	free(seq11);
	free(seq22);

	p = aa->path + aa->path_len - 1;
	aa->start1 = p->i? p->i : 1;
	aa->end1 = aa->path->i;
	aa->start2 = p->j? p->j : 1;
	aa->end2 = aa->path->j;
	aa->cigar32 = aln_path2cigar32(aa->path, aa->path_len, &aa->n_cigar);

	return aa;
}
AlnAln *aln_stdaln(const char *seq1, const char *seq2, const AlnParam *ap, int type, int thres)
{
	return aln_stdaln_aux(seq1, seq2, ap, type, thres, -1, -1);
}

/* for backward compatibility */
uint16_t *aln_path2cigar(const path_t *path, int path_len, int *n_cigar)
{
	uint32_t *cigar32;
	uint16_t *cigar;
	int i;
	cigar32 = aln_path2cigar32(path, path_len, n_cigar);
	cigar = (uint16_t*)cigar32;
	for (i = 0; i < *n_cigar; ++i)
		cigar[i] = (cigar32[i]&0xf)<<14 | (cigar32[i]>>4&0x3fff);
	return cigar;
}

/* newly added functions (2009-07-21) */

int aln_extend_core(unsigned char *seq1, int len1, unsigned char *seq2, int len2, const AlnParam *ap,
					path_t *path, int *path_len, int G0, uint8_t *_mem)
{
	int q, r, qr, tmp_len;
	int32_t **s_array, *score_array;
	int is_overflow, of_base;
	uint32_t *eh;
	int i, j, end_i, end_j;
	int score, start, end;
	int *score_matrix, N_MATRIX_ROW;
	uint8_t *mem, *_p;

	/* initialize some align-related parameters. just for compatibility */
	q = ap->gap_open;
	r = ap->gap_ext;
	qr = q + r;
	score_matrix = ap->matrix;
	N_MATRIX_ROW = ap->row;

	if (len1 == 0 || len2 == 0) return -1;

	/* allocate memory */
	mem = _mem? _mem : calloc((len1 + 2) * (N_MATRIX_ROW + 1), 4);
	_p = mem;
	eh = (uint32_t*)_p, _p += 4 * (len1 + 2);
	s_array = calloc(N_MATRIX_ROW, sizeof(void*));
	for (i = 0; i != N_MATRIX_ROW; ++i)
		s_array[i] = (int32_t*)_p, _p += 4 * len1;
	/* initialization */
	aln_init_score_array(seq1, len1, N_MATRIX_ROW, score_matrix, s_array);
	tmp_len = len1 + 1;
	start = 1; end = 2;
	end_i = end_j = 0;
	score = 0;
	is_overflow = of_base = 0;
	/* convert the coordinate */
	--seq1; --seq2;
	for (i = 0; i != N_MATRIX_ROW; ++i) --s_array[i];

	/* dynamic programming */
	memset(eh, 0, 4 * (len1 + 2));
	eh[1] = (uint32_t)G0<<16;
	for (j = 1; j <= len2; ++j) {
		int _start, _end;
		int h1 = 0, f = 0;
		score_array = s_array[seq2[j]];
		/* set start and end */
		_start = j - ap->band_width;
		if (_start < 1) _start = 1;
		if (_start > start) start = _start;
		_end = j + ap->band_width;
		if (_end > len1 + 1) _end = len1 + 1;
		if (_end < end) end = _end;
		if (start == end) break;
		/* adjust eh[] array if overflow occurs. */
		if (is_overflow) {
			int tmp, tmp2;
			score -= LOCAL_OVERFLOW_REDUCE;
			of_base += LOCAL_OVERFLOW_REDUCE;
			is_overflow = 0;
			for (i = start; i <= end; ++i) {
				uint32_t *s = &eh[i];
				tmp = *s >> 16; tmp2 = *s & 0xffff;
				if (tmp2 < LOCAL_OVERFLOW_REDUCE) tmp2 = 0;
				else tmp2 -= LOCAL_OVERFLOW_REDUCE;
				if (tmp < LOCAL_OVERFLOW_REDUCE) tmp = 0;
				else tmp -= LOCAL_OVERFLOW_REDUCE;
				*s = (tmp << 16) | tmp2;
			}
		}
		_start = _end = 0;
		/* the inner loop */
		for (i = start; i < end; ++i) {
			/* At the beginning of each cycle:
			     eh[i] -> h[j-1,i-1]<<16 | e[j,i]
				 f     -> f[j,i]
				 h1    -> h[j,i-1]
			*/
			uint32_t *s = &eh[i];
			int h = (int)(*s >> 16);
			int e = *s & 0xffff; /* this is e[j,i] */
			*s = (uint32_t)h1 << 16; /* eh[i] now stores h[j,i-1]<<16 */
			h += h? score_array[i] : 0; /* this is left_core() specific */
			/* calculate h[j,i]; don't need to test 0, as {e,f}>=0 */
			h = h > e? h : e;
			h = h > f? h : f; /* h now is h[j,i] */
			h1 = h;
			if (h > 0) {
				if (_start == 0) _start = i;
				_end = i;
				if (score < h) {
					score = h; end_i = i; end_j = j;
					if (score > LOCAL_OVERFLOW_THRESHOLD) is_overflow = 1;
				}
			}
			/* calculate e[j+1,i] and f[j,i+1] */
			h -= qr;
			h = h > 0? h : 0;
			e -= r;
			e = e > h? e : h;
			f -= r;
			f = f > h? f : h;
			*s |= e;
		}			
		eh[end] = h1 << 16;
		/* recalculate start and end, the boundaries of the band */
		if (_end <= 0) break; /* no cell in this row has a positive score */
		start = _start;
		end = _end + 3;
	}

	score += of_base - 1;
	if (score <= 0) {
		if (path_len) *path_len = 0;
		goto end_left_func;
	}

	if (path == 0) goto end_left_func;

	if (path_len == 0) {
		path[0].i = end_i; path[0].j = end_j;
		goto end_left_func;
	}

	{ /* call global alignment to fill the path */
		int score_g = 0;
		j = (end_i - 1 > end_j - 1)? end_i - 1 : end_j - 1;
		++j; /* j is the maximum band_width */
		for (i = ap->band_width;; i <<= 1) {
			AlnParam ap_real = *ap;
			ap_real.gap_end = -1;
			ap_real.band_width = i;
			score_g = aln_global_core(seq1 + 1, end_i, seq2 + 1, end_j, &ap_real, path, path_len);
			if (score == score_g) break;
			if (i > j) break;
		}
		if (score > score_g)
			fprintf(stderr, "[aln_left_core] no suitable bandwidth: %d < %d\n", score_g, score);
		score = score_g;
	}

end_left_func:
	/* free */
	free(s_array);
	if (!_mem) free(mem);
	return score;
}

uint32_t *aln_path2cigar32(const path_t *path, int path_len, int *n_cigar)
{
	int i, n;
	uint32_t *cigar;
	unsigned char last_type;

	if (path_len == 0 || path == 0) {
		*n_cigar = 0;
		return 0;
	}

	last_type = path->ctype;
	for (i = n = 1; i < path_len; ++i) {
		if (last_type != path[i].ctype) ++n;
		last_type = path[i].ctype;
	}
	*n_cigar = n;
	cigar = (uint32_t*)malloc(*n_cigar * 4);

	cigar[0] = 1u << 4 | path[path_len-1].ctype;
	last_type = path[path_len-1].ctype;
	for (i = path_len - 2, n = 0; i >= 0; --i) {
		if (path[i].ctype == last_type) cigar[n] += 1u << 4;
		else {
			cigar[++n] = 1u << 4 | path[i].ctype;
			last_type = path[i].ctype;
		}
	}

	return cigar;
}

#ifdef STDALN_MAIN
int main()
{
	AlnAln *aln_local, *aln_global, *aln_left;
	int i;

	aln_local  = aln_stdaln("CGTGCGATGCactgCATACGGCTCGCCTAGATCA", "AAGGGATGCTCTGCATCgCTCGGCTAGCTGT", &aln_param_blast, 0, 1);
	aln_global = aln_stdaln("CGTGCGATGCactgCATACGGCTCGCCTAGATCA", "AAGGGATGCTCTGCATCGgCTCGGCTAGCTGT", &aln_param_blast, 1, 1);
//	aln_left   = aln_stdaln(     "GATGCACTGCATACGGCTCGCCTAGATCA",     "GATGCTCTGCATCGgCTCGGCTAGCTGT", &aln_param_blast, 2, 1);
	aln_left   = aln_stdaln("CACCTTCGACTCACGTCTCATTCTCGGAGTCGAGTGGACGGTCCCTCATACACGAACAGGTTC",
							"CACCTTCGACTTTCACCTCTCATTCTCGGACTCGAGTGGACGGTCCCTCATCCAAGAACAGGGTCTGTGAAA", &aln_param_blast, 2, 1);

	printf(">%d,%d\t%d,%d\n", aln_local->start1, aln_local->end1, aln_local->start2, aln_local->end2);
	printf("%s\n%s\n%s\n", aln_local->out1, aln_local->outm, aln_local->out2);

	printf(">%d,%d\t%d,%d\t", aln_global->start1, aln_global->end1, aln_global->start2, aln_global->end2);
	for (i = 0; i != aln_global->n_cigar; ++i)
		printf("%d%c", aln_global->cigar32[i]>>4, "MID"[aln_global->cigar32[i]&0xf]);
	printf("\n%s\n%s\n%s\n", aln_global->out1, aln_global->outm, aln_global->out2);

	printf(">%d\t%d,%d\t%d,%d\t", aln_left->score, aln_left->start1, aln_left->end1, aln_left->start2, aln_left->end2);
	for (i = 0; i != aln_left->n_cigar; ++i)
		printf("%d%c", aln_left->cigar32[i]>>4, "MID"[aln_left->cigar32[i]&0xf]);
	printf("\n%s\n%s\n%s\n", aln_left->out1, aln_left->outm, aln_left->out2);

	aln_free_AlnAln(aln_local);
	aln_free_AlnAln(aln_global);
	aln_free_AlnAln(aln_left);
	return 0;
}
#endif
