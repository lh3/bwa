#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include "main.h"
#include "utils.h"

// -------------------

#ifndef PACKAGE_VERS
# define PACKAGE_VERS 0.6.1-r104-tpx
#endif

#define mkstr(s) #s
#define mkxstr(s) mkstr(s)
#define PACKAGE_VERSION mkxstr(PACKAGE_VERS)
#ifndef BLDDATE
# define BLDDATE unknown
#endif
#ifndef SVNURL
# define SVNURL unknown
#endif
#ifndef SVNREV
# define SVNREV unknown
#endif
char __attribute__((used)) svnid[] = mkxstr(@(#)$Id: bwa PACKAGE_VERS build-date: BLDDATE svn-url: SVNURL svn-rev: SVNREV $);

time_t _prog_start = 1;
char bwaversionstr[200] = { "" };
char bwablddatestr[200] = { "" };

// -------------------

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
	struct timeval st;
	int srtn = 0;
	int64_t maxrsskb = 0L;

	t_real = realtime();

	if (argc < 2) return usage();

	// ---------------

        gettimeofday(&st, NULL);
        _prog_start = st.tv_sec * 1000000L + (time_t)st.tv_usec;

        sprintf(bwaversionstr,"%s-%s",mkxstr(PACKAGE_VERS),mkxstr(SVNREV));
        sprintf(bwablddatestr,"%s",mkxstr(BLDDATE));

        for(i=1;i<argc;i++){
          if(strncmp(argv[i],"-ver",4) == 0){
            fprintf(stdout,"BWA program (%s)\n", bwaversionstr);
            return 0;
          }
        }

	// ---------------

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

		fprintf(stderr, "\n");

		// ---------------

		srtn = getmaxrss(&maxrsskb);
		if(srtn == 0)
			fprintf(stderr,"[%s] Mem RSS max: %ld mb\n",__func__, (long)maxrsskb / 1024L);

		// ---------------

		fprintf(stderr, "[%s] Real time: %.3f sec; CPU: %.3f sec\n", __func__, realtime() - t_real, cputime());
	}

	return 0;
}
