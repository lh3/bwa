#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "bwt.h"
#include "bwtsw2.h"
#include "utils.h"
#include "bwa.h"

int bwa_bwtsw2(int argc, char *argv[])
{
	bsw2opt_t *opt;
	bwaidx_t *idx;
	int c;

	opt = bsw2_init_opt();
	srand48(11);
	while ((c = getopt(argc, argv, "q:r:a:b:t:T:w:d:z:m:s:c:N:Hf:MI:SG:C")) >= 0) {
		switch (c) {
		case 'q': opt->q = atoi(optarg); break;
		case 'r': opt->r = atoi(optarg); break;
		case 'a': opt->a = atoi(optarg); break;
		case 'b': opt->b = atoi(optarg); break;
		case 'w': opt->bw = atoi(optarg); break;
		case 'T': opt->t = atoi(optarg); break;
		case 't': opt->n_threads = atoi(optarg); break;
		case 'z': opt->z = atoi(optarg); break;
		case 's': opt->is = atoi(optarg); break;
		case 'm': opt->mask_level = atof(optarg); break;
		case 'c': opt->coef = atof(optarg); break;
		case 'N': opt->t_seeds = atoi(optarg); break;
		case 'M': opt->multi_2nd = 1; break;
		case 'H': opt->hard_clip = 1; break;
		case 'f': xreopen(optarg, "w", stdout); break;
		case 'I': opt->max_ins = atoi(optarg); break;
		case 'S': opt->skip_sw = 1; break;
		case 'C': opt->cpy_cmt = 1; break;
		case 'G': opt->max_chain_gap = atoi(optarg); break;
		default: return 1;
		}
	}
	opt->qr = opt->q + opt->r;

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa bwasw [options] <target.prefix> <query.fa> [query2.fa]\n\n");
		fprintf(stderr, "Options: -a INT   score for a match [%d]\n", opt->a);
		fprintf(stderr, "         -b INT   mismatch penalty [%d]\n", opt->b);
		fprintf(stderr, "         -q INT   gap open penalty [%d]\n", opt->q);
		fprintf(stderr, "         -r INT   gap extension penalty [%d]\n", opt->r);
		fprintf(stderr, "         -w INT   band width [%d]\n", opt->bw);
		fprintf(stderr, "         -m FLOAT mask level [%.2f]\n", opt->mask_level);
		fprintf(stderr, "\n");
		fprintf(stderr, "         -t INT   number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "         -f FILE  file to output results to instead of stdout\n");
		fprintf(stderr, "         -H       in SAM output, use hard clipping instead of soft clipping\n");
		fprintf(stderr, "         -C       copy FASTA/Q comment to SAM output\n");
		fprintf(stderr, "         -M       mark multi-part alignments as secondary\n");
		fprintf(stderr, "         -S       skip Smith-Waterman read pairing\n");
		fprintf(stderr, "         -I INT   ignore pairs with insert >=INT for inferring the size distr [%d]\n", opt->max_ins);
		fprintf(stderr, "\n");
		fprintf(stderr, "         -T INT   score threshold divided by a [%d]\n", opt->t);
		fprintf(stderr, "         -c FLOAT coefficient of length-threshold adjustment [%.1f]\n", opt->coef);
		fprintf(stderr, "         -z INT   Z-best [%d]\n", opt->z);
		fprintf(stderr, "         -s INT   maximum seeding interval size [%d]\n", opt->is);
		fprintf(stderr, "         -N INT   # seeds to trigger rev aln; 2*INT is also the chaining threshold [%d]\n", opt->t_seeds);
		fprintf(stderr, "         -G INT   maximum gap size during chaining [%d]\n", opt->max_chain_gap);
		fprintf(stderr, "\n");
		fprintf(stderr, "Note: For long Illumina, 454 and Sanger reads, assembly contigs, fosmids and\n");
		fprintf(stderr, "      BACs, the default setting usually works well. For the current PacBio\n");
		fprintf(stderr, "      reads (end of 2010), '-b5 -q2 -r1 -z10' is recommended. One may also\n");
		fprintf(stderr, "      increase '-z' for better sensitivity.\n");
		fprintf(stderr, "\n");

		return 1;
	}

	// adjust opt for opt->a
	opt->t *= opt->a;
	opt->coef *= opt->a;

	if ((idx = bwa_idx_load(argv[optind], BWA_IDX_BWT|BWA_IDX_BNS)) == 0) return 1;
	bsw2_aln(opt, idx->bns, idx->bwt, argv[optind+1], optind+2 < argc? argv[optind+2] : 0);
	bwa_idx_destroy(idx);
	free(opt);
	
	return 0;
}
