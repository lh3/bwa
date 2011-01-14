/* The MIT License

   Copyright (c) 2003-2006, 2008, by Heng Li <lh3lh3@gmail.com>

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

/*
  2009-07-23, 0.10.0

  - Use 32-bit to store CIGAR

  - Report suboptimal aligments

  - Implemented half-fixed-half-open DP

  2009-04-26, 0.9.10

  - Allow to set a threshold for local alignment

  2009-02-18, 0.9.9

  - Fixed a bug when no residue matches

  2008-08-04, 0.9.8

  - Fixed the wrong declaration of aln_stdaln_aux()

  - Avoid 0 coordinate for global alignment

  2008-08-01, 0.9.7

  - Change gap_end penalty to 5 in aln_param_bwa

  - Add function to convert path_t to the CIGAR format

  2008-08-01, 0.9.6

  - The first gap now costs (gap_open+gap_ext), instead of
    gap_open. Scoring systems are modified accordingly.

  - Gap end is now correctly handled. Previously it is not correct.

  - Change license to MIT.

 */

#ifndef LH3_STDALN_H_
#define LH3_STDALN_H_


#define STDALN_VERSION 0.11.0

#include <stdint.h>

#define FROM_M 0
#define FROM_I 1
#define FROM_D 2
#define FROM_S 3

#define ALN_TYPE_LOCAL  0
#define ALN_TYPE_GLOBAL 1
#define ALN_TYPE_EXTEND 2

/* This is the smallest integer. It might be CPU-dependent in very RARE cases. */
#define MINOR_INF -1073741823

typedef struct
{
	int gap_open;
	int gap_ext;
	int gap_end;

	int *matrix;
	int row;
	int band_width;
} AlnParam;

typedef struct
{
	int i, j;
	unsigned char ctype;
} path_t;

typedef struct
{
	path_t *path; /* for advanced users... :-) */
	int path_len; /* for advanced users... :-) */
	int start1, end1; /* start and end of the first sequence, coordinations are 1-based */
	int start2, end2; /* start and end of the second sequence, coordinations are 1-based */
	int score, subo; /* score */

	char *out1, *out2; /* print them, and then you will know */
	char *outm;

	int n_cigar;
	uint32_t *cigar32;
} AlnAln;

#ifdef __cplusplus
extern "C" {
#endif

	AlnAln *aln_stdaln_aux(const char *seq1, const char *seq2, const AlnParam *ap,
						   int type, int do_align, int len1, int len2);
	AlnAln *aln_stdaln(const char *seq1, const char *seq2, const AlnParam *ap, int type, int do_align);
	void aln_free_AlnAln(AlnAln *aa);

	int aln_global_core(unsigned char *seq1, int len1, unsigned char *seq2, int len2, const AlnParam *ap,
						path_t *path, int *path_len);
	int aln_local_core(unsigned char *seq1, int len1, unsigned char *seq2, int len2, const AlnParam *ap,
					   path_t *path, int *path_len, int _thres, int *_subo);
	int aln_extend_core(unsigned char *seq1, int len1, unsigned char *seq2, int len2, const AlnParam *ap,
						path_t *path, int *path_len, int G0, uint8_t *_mem);
	uint16_t *aln_path2cigar(const path_t *path, int path_len, int *n_cigar);
	uint32_t *aln_path2cigar32(const path_t *path, int path_len, int *n_cigar);

#ifdef __cplusplus
}
#endif

/********************
 * global variables *
 ********************/

extern AlnParam aln_param_bwa;   /* = { 37,  9,  0, aln_sm_maq, 5, 50 }; */
extern AlnParam aln_param_blast; /* = {  5,  2,  2, aln_sm_blast, 5, 50 }; */
extern AlnParam aln_param_nt2nt; /* = { 10,  2,  2, aln_sm_nt, 16, 75 }; */
extern AlnParam aln_param_aa2aa; /* = { 20, 19, 19, aln_sm_read, 16, 75 }; */
extern AlnParam aln_param_rd2rd; /* = { 12,  2,  2, aln_sm_blosum62, 22, 50 }; */

/* common nucleotide score matrix for 16 bases */
extern int           aln_sm_nt[], aln_sm_bwa[];

/* BLOSUM62 and BLOSUM45 */
extern int           aln_sm_blosum62[], aln_sm_blosum45[];

/* common read for 16 bases. note that read alignment is quite different from common nucleotide alignment */
extern int           aln_sm_read[];

/* human-mouse score matrix for 4 bases */
extern int           aln_sm_hs[];

#endif
