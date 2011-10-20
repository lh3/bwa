/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

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

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "bntseq.h"
#include "utils.h"
#include "main.h"
#include "bwt.h"

#ifdef _DIVBWT
#include "divsufsort.h"
#endif

int is_bwt(ubyte_t *T, int n);

int64_t bwa_seq_len(const char *fn_pac)
{
	FILE *fp;
	int64_t pac_len;
	ubyte_t c;
	fp = xopen(fn_pac, "rb");
	fseek(fp, -1, SEEK_END);
	pac_len = ftell(fp);
	fread(&c, 1, 1, fp);
	fclose(fp);
	return (pac_len - 1) * 4 + (int)c;
}

bwt_t *bwt_pac2bwt(const char *fn_pac, int use_is)
{
	bwt_t *bwt;
	ubyte_t *buf, *buf2;
	int i, pac_size;
	FILE *fp;

	// initialization
	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	bwt->seq_len = bwa_seq_len(fn_pac);
	bwt->bwt_size = (bwt->seq_len + 15) >> 4;
	fp = xopen(fn_pac, "rb");

	// prepare sequence
	pac_size = (bwt->seq_len>>2) + ((bwt->seq_len&3) == 0? 0 : 1);
	buf2 = (ubyte_t*)calloc(pac_size, 1);
	fread(buf2, 1, pac_size, fp);
	fclose(fp);
	memset(bwt->L2, 0, 5 * 4);
	buf = (ubyte_t*)calloc(bwt->seq_len + 1, 1);
	for (i = 0; i < bwt->seq_len; ++i) {
		buf[i] = buf2[i>>2] >> ((3 - (i&3)) << 1) & 3;
		++bwt->L2[1+buf[i]];
	}
	for (i = 2; i <= 4; ++i) bwt->L2[i] += bwt->L2[i-1];
	free(buf2);

	// Burrows-Wheeler Transform
	if (use_is) {
		bwt->primary = is_bwt(buf, bwt->seq_len);
	} else {
#ifdef _DIVBWT
		bwt->primary = divbwt(buf, buf, 0, bwt->seq_len);
#else
		err_fatal_simple("libdivsufsort is not compiled in.");
#endif
	}
	bwt->bwt = (u_int32_t*)calloc(bwt->bwt_size, 4);
	for (i = 0; i < bwt->seq_len; ++i)
		bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);
	free(buf);
	return bwt;
}

int bwa_pac2bwt(int argc, char *argv[])
{
	bwt_t *bwt;
	int c, use_is = 1;
	while ((c = getopt(argc, argv, "d")) >= 0) {
		switch (c) {
		case 'd': use_is = 0; break;
		default: return 1;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: bwa pac2bwt [-d] <in.pac> <out.bwt>\n");
		return 1;
	}
	bwt = bwt_pac2bwt(argv[optind], use_is);
	bwt_dump_bwt(argv[optind+1], bwt);
	bwt_destroy(bwt);
	return 0;
}

#define bwt_B00(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

void bwt_bwtupdate_core(bwt_t *bwt)
{
	bwtint_t i, k, c[4], n_occ;
	uint32_t *buf;

	n_occ = (bwt->seq_len + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
	bwt->bwt_size += n_occ * sizeof(bwtint_t); // the new size
	buf = (uint32_t*)calloc(bwt->bwt_size, 4); // will be the new bwt
	c[0] = c[1] = c[2] = c[3] = 0;
	for (i = k = 0; i < bwt->seq_len; ++i) {
		if (i % OCC_INTERVAL == 0) {
			memcpy(buf + k, c, sizeof(bwtint_t) * 4);
			k += sizeof(bwtint_t); // in fact: sizeof(bwtint_t)=4*(sizeof(bwtint_t)/4)
		}
		if (i % 16 == 0) buf[k++] = bwt->bwt[i/16]; // 16 == sizeof(uint32_t)/2
		++c[bwt_B00(bwt, i)];
	}
	// the last element
	memcpy(buf + k, c, sizeof(bwtint_t) * 4);
	xassert(k + sizeof(bwtint_t) == bwt->bwt_size, "inconsistent bwt_size");
	// update bwt
	free(bwt->bwt); bwt->bwt = buf;
}

int bwa_bwtupdate(int argc, char *argv[])
{
	bwt_t *bwt;
	if (argc < 2) {
		fprintf(stderr, "Usage: bwa bwtupdate <the.bwt>\n");
		return 1;
	}
	bwt = bwt_restore_bwt(argv[1]);
	bwt_bwtupdate_core(bwt);
	bwt_dump_bwt(argv[1], bwt);
	bwt_destroy(bwt);
	return 0;
}

const int nst_color_space_table[] = { 4, 0, 0, 1, 0, 2, 3, 4, 0, 3, 2, 4, 1, 4, 4, 4};

/* this function is not memory efficient, but this will make life easier
   Ideally we should also change .amb files as one 'N' in the nucleotide
   sequence leads to two ambiguous colors. I may do this later... */
uint8_t *bwa_pac2cspac_core(const bntseq_t *bns)
{
	uint8_t *pac, *cspac;
	bwtint_t i;
	int c1, c2;
	pac = (uint8_t*)calloc(bns->l_pac/4 + 1, 1);
	cspac = (uint8_t*)calloc(bns->l_pac/4 + 1, 1);
	fread(pac, 1, bns->l_pac/4+1, bns->fp_pac);
	rewind(bns->fp_pac);
	c1 = pac[0]>>6; cspac[0] = c1<<6;
	for (i = 1; i < bns->l_pac; ++i) {
		c2 = pac[i>>2] >> (~i&3)*2 & 3;
		cspac[i>>2] |= nst_color_space_table[(1<<c1)|(1<<c2)] << (~i&3)*2;
		c1 = c2;
	}
	free(pac);
	return cspac;
}

int bwa_pac2cspac(int argc, char *argv[])
{
	bntseq_t *bns;
	uint8_t *cspac, ct;
	char *str;
	FILE *fp;

	if (argc < 3) {
		fprintf(stderr, "Usage: bwa pac2cspac <in.nt.prefix> <out.cs.prefix>\n");
		return 1;
	}
	bns = bns_restore(argv[1]);
	cspac = bwa_pac2cspac_core(bns);
	bns_dump(bns, argv[2]);
	// now write cspac
	str = (char*)calloc(strlen(argv[2]) + 5, 1);
	strcat(strcpy(str, argv[2]), ".pac");
	fp = xopen(str, "wb");
	fwrite(cspac, 1, bns->l_pac/4 + 1, fp);
	ct = bns->l_pac % 4;
	fwrite(&ct, 1, 1, fp);	
	fclose(fp);
	bns_destroy(bns);
	free(cspac);
	return 0;
}

int bwa_bwt2sa(int argc, char *argv[])
{
	bwt_t *bwt;
	int c, sa_intv = 32;
	while ((c = getopt(argc, argv, "i:")) >= 0) {
		switch (c) {
		case 'i': sa_intv = atoi(optarg); break;
		default: return 1;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: bwa bwt2sa [-i %d] <in.bwt> <out.sa>\n", sa_intv);
		return 1;
	}
	bwt = bwt_restore_bwt(argv[optind]);
	bwt_cal_sa(bwt, sa_intv);
	bwt_dump_sa(argv[optind+1], bwt);
	bwt_destroy(bwt);
	return 0;
}
