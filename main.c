#include <stdio.h>
#include <string.h>
#include "main.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.5.9-r16"
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
	fprintf(stderr, "\n");
	fprintf(stderr, "         fa2pac        convert FASTA to PAC format\n");
	fprintf(stderr, "         pac2bwt       generate BWT from PAC\n");
	fprintf(stderr, "         pac2bwtgen    alternative algorithm for generating BWT\n");
	fprintf(stderr, "         bwtupdate     update .bwt to the new format\n");
	fprintf(stderr, "         pac_rev       generate reverse PAC\n");
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
	if (argc < 2) return usage();
	if (strcmp(argv[1], "fa2pac") == 0) return bwa_fa2pac(argc-1, argv+1);
	else if (strcmp(argv[1], "pac2bwt") == 0) return bwa_pac2bwt(argc-1, argv+1);
	else if (strcmp(argv[1], "pac2bwtgen") == 0) return bwt_bwtgen_main(argc-1, argv+1);
	else if (strcmp(argv[1], "bwtupdate") == 0) return bwa_bwtupdate(argc-1, argv+1);
	else if (strcmp(argv[1], "pac_rev") == 0) return bwa_pac_rev(argc-1, argv+1);
	else if (strcmp(argv[1], "bwt2sa") == 0) return bwa_bwt2sa(argc-1, argv+1);
	else if (strcmp(argv[1], "index") == 0) return bwa_index(argc-1, argv+1);
	else if (strcmp(argv[1], "aln") == 0) return bwa_aln(argc-1, argv+1);
	else if (strcmp(argv[1], "sw") == 0) return bwa_stdsw(argc-1, argv+1);
	else if (strcmp(argv[1], "samse") == 0) return bwa_sai2sam_se(argc-1, argv+1);
	else if (strcmp(argv[1], "sampe") == 0) return bwa_sai2sam_pe(argc-1, argv+1);
	else if (strcmp(argv[1], "pac2cspac") == 0) return bwa_pac2cspac(argc-1, argv+1);
	else if (strcmp(argv[1], "stdsw") == 0) return bwa_stdsw(argc-1, argv+1);
	else if (strcmp(argv[1], "bwtsw2") == 0) return bwa_bwtsw2(argc-1, argv+1);
	else if (strcmp(argv[1], "dbwtsw") == 0) return bwa_bwtsw2(argc-1, argv+1);
	else if (strcmp(argv[1], "bwasw") == 0) return bwa_bwtsw2(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;
}
