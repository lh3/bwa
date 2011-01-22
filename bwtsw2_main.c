#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "bwt.h"
#include "bwtsw2.h"
#include "utils.h"

int bwa_bwtsw2(int argc, char *argv[])
{
	bsw2opt_t *opt;
	bwt_t *target[2];
	char buf[1024];
	bntseq_t *bns;
	int c;

	opt = bsw2_init_opt();
	srand48(11);
	while ((c = getopt(argc, argv, "q:r:a:b:t:T:w:d:z:m:y:s:c:N:Hf:")) >= 0) {
		switch (c) {
		case 'q': opt->q = atoi(optarg); break;
		case 'r': opt->r = atoi(optarg); break;
		case 'a': opt->a = atoi(optarg); break;
		case 'b': opt->b = atoi(optarg); break;
		case 'w': opt->bw = atoi(optarg); break;
		case 'T': opt->t = atoi(optarg); break;
		case 't': opt->n_threads = atoi(optarg); break;
		case 'z': opt->z = atoi(optarg); break;
		case 'y': opt->yita = atof(optarg); break;
		case 's': opt->is = atoi(optarg); break;
		case 'm': opt->mask_level = atof(optarg); break;
		case 'c': opt->coef = atof(optarg); break;
		case 'N': opt->t_seeds = atoi(optarg); break;
		case 'H': opt->hard_clip = 1; break;
		case 'f': xreopen(optarg, "w", stdout); break;
		}
	}
	opt->qr = opt->q + opt->r;

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa bwasw [options] <target.prefix> <query.fa>\n\n");
		fprintf(stderr, "Options: -a INT   score for a match [%d]\n", opt->a);
		fprintf(stderr, "         -b INT   mismatch penalty [%d]\n", opt->b);
		fprintf(stderr, "         -q INT   gap open penalty [%d]\n", opt->q);
		fprintf(stderr, "         -r INT   gap extension penalty [%d]\n", opt->r);
//		fprintf(stderr, "         -y FLOAT error recurrence coef. (4..16) [%.1f]\n", opt->yita);
		fprintf(stderr, "\n");
		fprintf(stderr, "         -t INT   number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "         -s INT   size of a chunk of reads [%d]\n", opt->chunk_size);
		fprintf(stderr, "\n");
		fprintf(stderr, "         -w INT   band width [%d]\n", opt->bw);
		fprintf(stderr, "         -m FLOAT mask level [%.2f]\n", opt->mask_level);
		fprintf(stderr, "\n");
		fprintf(stderr, "         -T INT   score threshold divided by a [%d]\n", opt->t);
		fprintf(stderr, "         -s INT   maximum seeding interval size [%d]\n", opt->is);
		fprintf(stderr, "         -z INT   Z-best [%d]\n", opt->z);
		fprintf(stderr, "         -N INT   # seeds to trigger reverse alignment [%d]\n", opt->t_seeds);
		fprintf(stderr, "         -c FLOAT coefficient of length-threshold adjustment [%.1f]\n", opt->coef);
		fprintf(stderr, "         -H       in SAM output, use hard clipping rather than soft\n");
        fprintf(stderr, "         -f FILE  file to output results to instead of stdout\n\n");
		fprintf(stderr, "Note: For long Illumina, 454 and Sanger reads, assembly contigs, fosmids and\n");
		fprintf(stderr, "      BACs, the default setting usually works well. For the current PacBio\n");
		fprintf(stderr, "      reads (end of 2010), '-b5 -q2 -r1 -z10' is recommended. One may also\n");
		fprintf(stderr, "      increase '-z' for better sensitivity.\n");
		fprintf(stderr, "\n");

		if (0) {
			double c, theta, eps, delta;
			c = opt->a / log(opt->yita);
			theta = exp(-opt->b / c) / opt->yita;
			eps = exp(-opt->q / c);
			delta = exp(-opt->r / c);
			fprintf(stderr, "mismatch: %lf, gap_open: %lf, gap_ext: %lf\n\n",
					theta, eps, delta);
		}
		return 1;
	}

	// adjust opt for opt->a
	opt->t *= opt->a;
	opt->coef *= opt->a;

	strcpy(buf, argv[optind]); target[0] = bwt_restore_bwt(strcat(buf, ".bwt"));
	strcpy(buf, argv[optind]); bwt_restore_sa(strcat(buf, ".sa"), target[0]);
	strcpy(buf, argv[optind]); target[1] = bwt_restore_bwt(strcat(buf, ".rbwt"));
	strcpy(buf, argv[optind]); bwt_restore_sa(strcat(buf, ".rsa"), target[1]);
	bns = bns_restore(argv[optind]);

	bsw2_aln(opt, bns, target, argv[optind+1]);

	bns_destroy(bns);
	bwt_destroy(target[0]); bwt_destroy(target[1]);
	free(opt);
	
	return 0;
}
