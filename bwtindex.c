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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>
#include "bntseq.h"
#include "bwt.h"
#include "main.h"
#include "utils.h"

bwt_t *bwt_pac2bwt(const char *fn_pac, int use_is);
void bwa_pac_rev_core(const char *fn, const char *fn_rev);

int bwa_index(int argc, char *argv[])
{
	char *prefix = 0, *str, *str2, *str3;
	int c, algo_type = 3, is_color = 0;
	clock_t t;

	while ((c = getopt(argc, argv, "ca:p:")) >= 0) {
		switch (c) {
		case 'a':
			if (strcmp(optarg, "div") == 0) algo_type = 1;
			else if (strcmp(optarg, "bwtsw") == 0) algo_type = 2;
			else if (strcmp(optarg, "is") == 0) algo_type = 3;
			else err_fatal(__func__, "unknown algorithm: '%s'.", optarg);
			break;
		case 'p': prefix = strdup(optarg); break;
		case 'c': is_color = 1; break;
		default: return 1;
		}
	}

	if (optind + 1 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa index [-a bwtsw|div|is] [-c] <in.fasta>\n\n");
		fprintf(stderr, "Options: -a STR    BWT construction algorithm: bwtsw or is [is]\n");
		fprintf(stderr, "         -p STR    prefix of the index [same as fasta name]\n");
		fprintf(stderr, "         -c        build color-space index\n\n");
		fprintf(stderr,	"Warning: `-a bwtsw' does not work for short genomes, while `-a is' and\n");
		fprintf(stderr, "         `-a div' do not work not for long genomes. Please choose `-a'\n");
		fprintf(stderr, "         according to the length of the genome.\n\n");
		return 1;
	}
	if (prefix == 0) prefix = strdup(argv[optind]);
	str  = (char*)calloc(strlen(prefix) + 10, 1);
	str2 = (char*)calloc(strlen(prefix) + 10, 1);
	str3 = (char*)calloc(strlen(prefix) + 10, 1);

	if (is_color == 0) { // nucleotide indexing
		gzFile fp = xzopen(argv[optind], "r");
		t = clock();
		fprintf(stderr, "[bwa_index] Pack FASTA... ");
		bns_fasta2bntseq(fp, prefix);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		gzclose(fp);
	} else { // color indexing
		gzFile fp = xzopen(argv[optind], "r");
		strcat(strcpy(str, prefix), ".nt");
		t = clock();
		fprintf(stderr, "[bwa_index] Pack nucleotide FASTA... ");
		bns_fasta2bntseq(fp, str);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		gzclose(fp);
		{
			char *tmp_argv[3];
			tmp_argv[0] = argv[0]; tmp_argv[1] = str; tmp_argv[2] = prefix;
			t = clock();
			fprintf(stderr, "[bwa_index] Convert nucleotide PAC to color PAC... ");
			bwa_pac2cspac(3, tmp_argv);
			fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		}
	}
	{
		strcpy(str, prefix); strcat(str, ".pac");
		strcpy(str2, prefix); strcat(str2, ".rpac");
		t = clock();
		fprintf(stderr, "[bwa_index] Reverse the packed sequence... ");
		bwa_pac_rev_core(str, str2);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	{
		strcpy(str, prefix); strcat(str, ".pac");
		strcpy(str2, prefix); strcat(str2, ".bwt");
		t = clock();
		fprintf(stderr, "[bwa_index] Construct BWT for the packed sequence...\n");
		if (algo_type == 2) bwt_bwtgen(str, str2);
		else if (algo_type == 1 || algo_type == 3) {
			bwt_t *bwt;
			bwt = bwt_pac2bwt(str, algo_type == 3);
			bwt_dump_bwt(str2, bwt);
			bwt_destroy(bwt);
		}
		fprintf(stderr, "[bwa_index] %.2f seconds elapse.\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	{
		strcpy(str, prefix); strcat(str, ".rpac");
		strcpy(str2, prefix); strcat(str2, ".rbwt");
		t = clock();
		fprintf(stderr, "[bwa_index] Construct BWT for the reverse packed sequence...\n");
		if (algo_type == 2) bwt_bwtgen(str, str2);
		else if (algo_type == 1 || algo_type == 3) {
			bwt_t *bwt;
			bwt = bwt_pac2bwt(str, algo_type == 3);
			bwt_dump_bwt(str2, bwt);
			bwt_destroy(bwt);
		}
		fprintf(stderr, "[bwa_index] %.2f seconds elapse.\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	{
		bwt_t *bwt;
		strcpy(str, prefix); strcat(str, ".bwt");
		t = clock();
		fprintf(stderr, "[bwa_index] Update BWT... ");
		bwt = bwt_restore_bwt(str);
		bwt_bwtupdate_core(bwt);
		bwt_dump_bwt(str, bwt);
		bwt_destroy(bwt);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	{
		bwt_t *bwt;
		strcpy(str, prefix); strcat(str, ".rbwt");
		t = clock();
		fprintf(stderr, "[bwa_index] Update reverse BWT... ");
		bwt = bwt_restore_bwt(str);
		bwt_bwtupdate_core(bwt);
		bwt_dump_bwt(str, bwt);
		bwt_destroy(bwt);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	{
		bwt_t *bwt;
		strcpy(str, prefix); strcat(str, ".bwt");
		strcpy(str3, prefix); strcat(str3, ".sa");
		t = clock();
		fprintf(stderr, "[bwa_index] Construct SA from BWT and Occ... ");
		bwt = bwt_restore_bwt(str);
		bwt_cal_sa(bwt, 32);
		bwt_dump_sa(str3, bwt);
		bwt_destroy(bwt);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	{
		bwt_t *bwt;
		strcpy(str, prefix); strcat(str, ".rbwt");
		strcpy(str3, prefix); strcat(str3, ".rsa");
		t = clock();
		fprintf(stderr, "[bwa_index] Construct SA from reverse BWT and Occ... ");
		bwt = bwt_restore_bwt(str);
		bwt_cal_sa(bwt, 32);
		bwt_dump_sa(str3, bwt);
		bwt_destroy(bwt);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	free(str3); free(str2); free(str); free(prefix);
	return 0;
}
