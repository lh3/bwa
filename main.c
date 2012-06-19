#include <stdio.h>
#include <string.h>
#include "main.h"
#include "utils.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.6.2-r126"
#endif

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: bwa (alignment via Burrows-Wheeler transformation)\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Heng Li <lh3@sanger.ac.uk>\n\n");
	fprintf(stderr, "Usage:   bwa <command> [options]\n\n");
	fprintf(stderr, "Command: index         index sequences in the FASTA format\n");
	fprintf(stderr, "         aln           gapped/ungapped alignment\n");
	fprintf(stderr, "         samse         generate alignment (single ended)\n");
	fprintf(stderr, "         sampe         generate alignment (paired ended)\n");
	fprintf(stderr, "         bwasw         BWA-SW for long queries\n");
	fprintf(stderr, "         fastmap       identify super-maximal exact matches\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "         fa2pac        convert FASTA to PAC format\n");
	fprintf(stderr, "         pac2bwt       generate BWT from PAC\n");
	fprintf(stderr, "         pac2bwtgen    alternative algorithm for generating BWT\n");
	fprintf(stderr, "         bwtupdate     update .bwt to the new format\n");
	fprintf(stderr, "         bwt2sa        generate SA from BWT and Occ\n");
	fprintf(stderr, "         pac2cspac     convert PAC to color-space PAC\n");
	fprintf(stderr, "         stdsw         standard SW/NW alignment\n");
	fprintf(stderr, "\n");
	return 1;
}

void bwa_print_sam_PG()
{
	printf("@PG\tID:bwa\tPN:bwa\tVN:%s\n", PACKAGE_VERSION);
}

int main(int argc, char *argv[])
{
	int i, ret;
	double t_real;
	t_real = realtime();
	if (argc < 2) return usage();
	if (strcmp(argv[1], "fa2pac") == 0) ret = bwa_fa2pac(argc-1, argv+1);
	else if (strcmp(argv[1], "pac2bwt") == 0) ret = bwa_pac2bwt(argc-1, argv+1);
	else if (strcmp(argv[1], "pac2bwtgen") == 0) ret = bwt_bwtgen_main(argc-1, argv+1);
	else if (strcmp(argv[1], "bwtupdate") == 0) ret = bwa_bwtupdate(argc-1, argv+1);
	else if (strcmp(argv[1], "bwt2sa") == 0) ret = bwa_bwt2sa(argc-1, argv+1);
	else if (strcmp(argv[1], "index") == 0) ret = bwa_index(argc-1, argv+1);
	else if (strcmp(argv[1], "aln") == 0) ret = bwa_aln(argc-1, argv+1);
	else if (strcmp(argv[1], "sw") == 0) ret = bwa_stdsw(argc-1, argv+1);
	else if (strcmp(argv[1], "samse") == 0) ret = bwa_sai2sam_se(argc-1, argv+1);
	else if (strcmp(argv[1], "sampe") == 0) ret = bwa_sai2sam_pe(argc-1, argv+1);
	else if (strcmp(argv[1], "pac2cspac") == 0) ret = bwa_pac2cspac(argc-1, argv+1);
	else if (strcmp(argv[1], "stdsw") == 0) ret = bwa_stdsw(argc-1, argv+1);
	else if (strcmp(argv[1], "bwtsw2") == 0) ret = bwa_bwtsw2(argc-1, argv+1);
	else if (strcmp(argv[1], "dbwtsw") == 0) ret = bwa_bwtsw2(argc-1, argv+1);
	else if (strcmp(argv[1], "bwasw") == 0) ret = bwa_bwtsw2(argc-1, argv+1);
	else if (strcmp(argv[1], "fastmap") == 0) ret = main_fastmap(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	err_fflush(stdout);
	err_fclose(stdout);
	if (ret == 0) {
		fprintf(stderr, "[%s] Version: %s\n", __func__, PACKAGE_VERSION);
		fprintf(stderr, "[%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fprintf(stderr, "\n[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	}
	return 0;
}
