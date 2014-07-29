#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"

#include "utils.h"

#include "tbb_mem.h"

#include "tbb/task_scheduler_init.h"
#include "tbb/pipeline.h"


extern "C" int main_mem_tbb(int argc, char *argv[])
{
	mem_opt_t *opt, opt0;
	int fd, fd2, i, c, n, copy_comment = 0;
	gzFile fp, fp2 = 0;
	kseq_t *ks, *ks2 = 0;
	bseq1_t *seqs;
	bwaidx_t *idx;
	char *p, *rg_line = 0;
	void *ko = 0, *ko2 = 0;
	int64_t n_processed = 0;
	mem_pestat_t pes[4], *pes0 = 0;

	memset(pes, 0, 4 * sizeof(mem_pestat_t));
	for (i = 0; i < 4; ++i) pes[i].failed = 1;

	opt = mem_opt_init();
	opt0.a = opt0.b = opt0.o_del = opt0.e_del = opt0.o_ins = opt0.e_ins = opt0.pen_unpaired = -1;
	opt0.pen_clip5 = opt0.pen_clip3 = opt0.zdrop = opt0.T = -1;
	while ((c = getopt(argc, argv, "epaMCSPHk:c:v:s:r:t:R:A:B:O:E:U:w:L:d:T:Q:D:m:I:")) >= 0) {
		if (c == 'k') opt->min_seed_len = atoi(optarg);
		else if (c == 'w') opt->w = atoi(optarg);
		else if (c == 'A') opt->a = atoi(optarg), opt0.a = 1;
		else if (c == 'B') opt->b = atoi(optarg), opt0.b = 1;
		else if (c == 'T') opt->T = atoi(optarg), opt0.T = 1;
		else if (c == 'U') opt->pen_unpaired = atoi(optarg), opt0.pen_unpaired = 1;
		else if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		else if (c == 'p') opt->flag |= MEM_F_PE;
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
		else if (c == 'S') opt->flag |= MEM_F_NO_RESCUE;
		else if (c == 'e') opt->flag |= MEM_F_NO_EXACT;
		else if (c == 'c') opt->max_occ = atoi(optarg);
		else if (c == 'd') opt->zdrop = atoi(optarg), opt0.zdrop = 1;
		else if (c == 'v') bwa_verbose = atoi(optarg);
		else if (c == 'r') opt->split_factor = atof(optarg);
		else if (c == 'D') opt->chain_drop_ratio = atof(optarg);
		else if (c == 'm') opt->max_matesw = atoi(optarg);
		else if (c == 's') opt->split_width = atoi(optarg);
		else if (c == 'C') copy_comment = 1;
		else if (c == 'Q') {
			opt->mapQ_coef_len = atoi(optarg);
			opt->mapQ_coef_fac = opt->mapQ_coef_len > 0? log(opt->mapQ_coef_len) : 0;
		} else if (c == 'O') {
			opt0.o_del = opt0.o_ins = 1;
			opt->o_del = opt->o_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->o_ins = strtol(p+1, &p, 10);
		} else if (c == 'E') {
			opt0.e_del = opt0.e_ins = 1;
			opt->e_del = opt->e_ins = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->e_ins = strtol(p+1, &p, 10);
		} else if (c == 'L') {
			opt0.pen_clip5 = opt0.pen_clip3 = 1;
			opt->pen_clip5 = opt->pen_clip3 = strtol(optarg, &p, 10);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				opt->pen_clip3 = strtol(p+1, &p, 10);
		} else if (c == 'R') {
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1; // FIXME: memory leak
		} else if (c == 'I') { // specify the insert size distribution
			pes0 = pes;
			pes[1].failed = 0;
			pes[1].avg = strtod(optarg, &p);
			pes[1].std = pes[1].avg * .1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].std = strtod(p+1, &p);
			pes[1].high = (int)(pes[1].avg + 4. * pes[1].std + .499);
			pes[1].low  = (int)(pes[1].avg - 4. * pes[1].std + .499);
			if (pes[1].low < 1) pes[1].low = 1;
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].high = (int)(strtod(p+1, &p) + .499);
			if (*p != 0 && ispunct(*p) && isdigit(p[1]))
				pes[1].low  = (int)(strtod(p+1, &p) + .499);
			if (bwa_verbose >= 3)
				fprintf(stderr, "[M::%s] mean insert size: %.3f, stddev: %.3f, max: %d, min: %d\n",
						__func__, pes[1].avg, pes[1].std, pes[1].high, pes[1].low);
		}
		else return 1;
	}
	if (opt->n_threads < 1) opt->n_threads = 1;
	if (optind + 1 >= argc || optind + 3 < argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT        number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -k INT        minimum seed length [%d]\n", opt->min_seed_len);
		fprintf(stderr, "       -w INT        band width for banded alignment [%d]\n", opt->w);
		fprintf(stderr, "       -d INT        off-diagonal X-dropoff [%d]\n", opt->zdrop);
		fprintf(stderr, "       -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
//		fprintf(stderr, "       -s INT        look for internal seeds inside a seed with less than INT occ [%d]\n", opt->split_width);
		fprintf(stderr, "       -c INT        skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
		fprintf(stderr, "       -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [%.2f]\n", opt->chain_drop_ratio);
		fprintf(stderr, "       -m INT        perform at most INT rounds of mate rescues for each read [%d]\n", opt->max_matesw);
		fprintf(stderr, "       -S            skip mate rescue\n");
		fprintf(stderr, "       -P            skip pairing; mate rescue performed unless -S also in use\n");
		fprintf(stderr, "       -e            discard full-length exact matches\n");
		fprintf(stderr, "       -A INT        score for a sequence match, which scales [-TdBOELU] unless overridden [%d]\n", opt->a);
		fprintf(stderr, "       -B INT        penalty for a mismatch [%d]\n", opt->b);
		fprintf(stderr, "       -O INT[,INT]  gap open penalties for deletions and insertions [%d,%d]\n", opt->o_del, opt->o_ins);
		fprintf(stderr, "       -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [%d,%d]\n", opt->e_del, opt->e_ins);
		fprintf(stderr, "       -L INT[,INT]  penalty for 5'- and 3'-end clipping [%d,%d]\n", opt->pen_clip5, opt->pen_clip3);
		fprintf(stderr, "       -U INT        penalty for an unpaired read pair [%d]\n", opt->pen_unpaired);
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p            first query file consists of interleaved paired-end sequences\n");
		fprintf(stderr, "       -R STR        read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "       -T INT        minimum score to output [%d]\n", opt->T);
		fprintf(stderr, "       -a            output all alignments for SE or unpaired PE\n");
		fprintf(stderr, "       -C            append FASTA/FASTQ comment to SAM output\n");
		fprintf(stderr, "       -M            mark shorter split hits as secondary\n\n");
		fprintf(stderr, "       -I FLOAT[,FLOAT[,INT[,INT]]]\n");
		fprintf(stderr, "                     specify the mean, standard deviation (10%% of the mean if absent), max\n");
		fprintf(stderr, "                     (4 sigma from the mean if absent) and min of the insert size distribution.\n");
		fprintf(stderr, "                     FR orientation only. [inferred]\n");
		fprintf(stderr, "\nNote: Please read the man page for detailed description of the command line and options.\n");
		fprintf(stderr, "\n");
		free(opt);
		return 1;
	}

	if (opt0.a == 1) { // matching score is changed
		if (opt0.b != 1) opt->b *= opt->a;
		if (opt0.T != 1) opt->T *= opt->a;
		if (opt0.o_del != 1) opt->o_del *= opt->a;
		if (opt0.e_del != 1) opt->e_del *= opt->a;
		if (opt0.o_ins != 1) opt->o_ins *= opt->a;
		if (opt0.e_ins != 1) opt->e_ins *= opt->a;
		if (opt0.zdrop != 1) opt->zdrop *= opt->a;
		if (opt0.pen_clip5 != 1) opt->pen_clip5 *= opt->a;
		if (opt0.pen_clip3 != 1) opt->pen_clip3 *= opt->a;
		if (opt0.pen_unpaired != 1) opt->pen_unpaired *= opt->a;
	}
	bwa_fill_scmat(opt->a, opt->b, opt->mat);
	if ((idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1; // FIXME: memory leak

	ko = kopen(argv[optind + 1], &fd);
	if (ko == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 1]);
		return 1;
	}
	fp = gzdopen(fd, "r");
	ks = kseq_init(fp);
	if (optind + 2 < argc) {
		if (opt->flag&MEM_F_PE) {
			if (bwa_verbose >= 2)
				fprintf(stderr, "[W::%s] when '-p' is in use, the second query file will be ignored.\n", __func__);
		} else {
			ko2 = kopen(argv[optind + 2], &fd2);
			if (ko2 == 0) {
				if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to open file `%s'.\n", __func__, argv[optind + 2]);
				return 1;
			}
			fp2 = gzdopen(fd2, "r");
			ks2 = kseq_init(fp2);
			opt->flag |= MEM_F_PE;
		}
	}

	std::string rg_line_;
	if(rg_line) rg_line_=rg_line;
	int n_threads=1;
	if(opt->n_threads>1) n_threads=opt->n_threads;
	tbb::task_scheduler_init tbb_init(n_threads);
	std::vector<kseq_t*> kss(2);
	//kss
	bwa::LogFilter<bwa::ReadBSeqFilter,tbb::flow_control&,bwa::BSeq*> read_bseq("ReadBSeq",bwa::ReadBSeqFilter({ks,ks2},*opt,copy_comment));
	bwa::LogFilter<bwa::MemParallelFilter,bwa::BSeq*,bwa::BSeq*> mem("MemParallel",bwa::MemParallelFilter(*idx->bwt,*idx->bns,idx->pac,*opt,pes0));
	bwa::LogFilter<bwa::WriteBSeqFilter,bwa::BSeq*,void> write_bseq("WriteBSeq",bwa::WriteBSeqFilter(std::cout,*idx->bns,rg_line_));
	write_bseq->write_hdr();

	tbb::parallel_pipeline(n_threads*2,
		tbb::make_filter<void,bwa::BSeq*>(tbb::filter::mode::serial_in_order,read_bseq)&
		tbb::make_filter<bwa::BSeq*,bwa::BSeq*>(tbb::filter::mode::parallel,mem)&
		tbb::make_filter<bwa::BSeq*,void>(tbb::filter::mode::serial_in_order,write_bseq));

	free(opt);
	bwa_idx_destroy(idx);
	kseq_destroy(ks);
	err_gzclose(fp); kclose(ko);
	if (ks2) {
		kseq_destroy(ks2);
		err_gzclose(fp2); kclose(ko2);
	}
	return 0;
}
