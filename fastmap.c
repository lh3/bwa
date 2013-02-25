#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

extern unsigned char nst_nt4_table[256];

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

int main_mem(int argc, char *argv[])
{
	mem_opt_t *opt;
	int fd, fd2, i, c, n, copy_comment = 0;
	gzFile fp, fp2 = 0;
	kseq_t *ks, *ks2 = 0;
	bseq1_t *seqs;
	bwaidx_t *idx;
	char *rg_line = 0;
	void *ko = 0, *ko2 = 0;

	opt = mem_opt_init();
	while ((c = getopt(argc, argv, "paMCPHk:c:v:s:r:t:R:")) >= 0) {
		if (c == 'k') opt->min_seed_len = atoi(optarg);
		else if (c == 't') opt->n_threads = atoi(optarg), opt->n_threads = opt->n_threads > 1? opt->n_threads : 1;
		else if (c == 'P') opt->flag |= MEM_F_NOPAIRING;
		else if (c == 'H') opt->flag |= MEM_F_HARDCLIP;
		else if (c == 'a') opt->flag |= MEM_F_ALL;
		else if (c == 'p') opt->flag |= MEM_F_PE;
		else if (c == 'M') opt->flag |= MEM_F_NO_MULTI;
		else if (c == 'c') opt->max_occ = atoi(optarg);
		else if (c == 'v') bwa_verbose = atoi(optarg);
		else if (c == 'r') opt->split_factor = atof(optarg);
		else if (c == 'C') copy_comment = 1;
		else if (c == 'R') {
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1; // FIXME: memory leak
		} else if (c == 's') opt->split_width = atoi(optarg);
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]\n\n");
		fprintf(stderr, "Algorithm options:\n\n");
		fprintf(stderr, "       -t INT     number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "       -k INT     minimum seed length [%d]\n", opt->min_seed_len);
		fprintf(stderr, "       -r FLOAT   look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
		fprintf(stderr, "       -s INT     look for internal seeds inside a seed with less than INT occ [%d]\n", opt->split_width);
		fprintf(stderr, "       -c INT     skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
		fprintf(stderr, "       -P         skip pairing; perform mate SW only\n");
		fprintf(stderr, "\nInput/output options:\n\n");
		fprintf(stderr, "       -p         first query file consists of interleaved paired-end sequences\n");
		fprintf(stderr, "       -R STR     read group header line such as '@RG\\tID:foo\\tSM:bar' [null]\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "       -v INT     verbose level: 1=error, 2=warning, 3=message, 4+=debugging [%d]\n", bwa_verbose);
		fprintf(stderr, "       -a         output all alignments for SE or unpaired PE\n");
		fprintf(stderr, "       -C         append FASTA/FASTQ comment to SAM output\n");
		fprintf(stderr, "       -H         hard clipping\n");
		fprintf(stderr, "       -M         mark shorter split hits as secondary (for Picard/GATK compatibility)\n");
		fprintf(stderr, "\n");
		free(opt);
		return 1;
	}

	mem_fill_scmat(opt->a, opt->b, opt->mat);
	if ((idx = bwa_idx_load(argv[optind], BWA_IDX_ALL)) == 0) return 1; // FIXME: memory leak
	bwa_print_sam_hdr(idx->bns, rg_line);

	ko = kopen(argv[optind + 1], &fd);
	fp = gzdopen(fd, "r");
	ks = kseq_init(fp);
	if (optind + 2 < argc) {
		if (opt->flag&MEM_F_PE) {
			if (bwa_verbose >= 2)
				fprintf(stderr, "[W::%s] when '-p' is in use, the second query file will be ignored.\n", __func__);
		} else {
			ko2 = kopen(argv[optind + 2], &fd2);
			fp2 = gzdopen(fd2, "r");
			ks2 = kseq_init(fp2);
			opt->flag |= MEM_F_PE;
		}
	}
	while ((seqs = bseq_read(opt->chunk_size * (ko2? 2 : 1), &n, ks, ks2)) != 0) {
		int64_t size = 0;
		if (!copy_comment)
			for (i = 0; i < n; ++i) {
				free(seqs[i].comment); seqs[i].comment = 0;
			}
		for (i = 0; i < n; ++i) size += seqs[i].l_seq;
		if (bwa_verbose >= 3)
			fprintf(stderr, "[M::%s] read %d sequences (%ld bp)...\n", __func__, n, (long)size);
		mem_process_seqs(opt, idx->bwt, idx->bns, idx->pac, n, seqs);
		free(seqs);
	}

	free(opt);
	bwa_idx_destroy(idx);
	kseq_destroy(ks);
	gzclose(fp); kclose(ko);
	if (ks2) {
		kseq_destroy(ks2);
		gzclose(fp2); kclose(ko2);
	}
	return 0;
}

int main_fastmap(int argc, char *argv[])
{
	int c, i, min_iwidth = 20, min_len = 17, print_seq = 0, split_width = 0;
	kseq_t *seq;
	bwtint_t k;
	gzFile fp;
	smem_i *itr;
	const bwtintv_v *a;
	bwaidx_t *idx;

	while ((c = getopt(argc, argv, "w:l:ps:")) >= 0) {
		switch (c) {
			case 's': split_width = atoi(optarg); break;
			case 'p': print_seq = 1; break;
			case 'w': min_iwidth = atoi(optarg); break;
			case 'l': min_len = atoi(optarg); break;
		}
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "Usage: bwa fastmap [-p] [-s splitWidth=%d] [-l minLen=%d] [-w maxSaSize=%d] <idxbase> <in.fq>\n", split_width, min_len, min_iwidth);
		return 1;
	}

	fp = gzopen(argv[optind + 1], "r");
	seq = kseq_init(fp);
	idx = bwa_idx_load(argv[optind], BWA_IDX_BWT|BWA_IDX_BNS);
	itr = smem_itr_init(idx->bwt);
	while (kseq_read(seq) >= 0) {
		printf("SQ\t%s\t%ld", seq->name.s, seq->seq.l);
		if (print_seq) {
			putchar('\t');
			puts(seq->seq.s);
		} else putchar('\n');
		for (i = 0; i < seq->seq.l; ++i)
			seq->seq.s[i] = nst_nt4_table[(int)seq->seq.s[i]];
		smem_set_query(itr, seq->seq.l, (uint8_t*)seq->seq.s);
		while ((a = smem_next(itr, min_len<<1, split_width)) != 0) {
			for (i = 0; i < a->n; ++i) {
				bwtintv_t *p = &a->a[i];
				if ((uint32_t)p->info - (p->info>>32) < min_len) continue;
				printf("EM\t%d\t%d\t%ld", (uint32_t)(p->info>>32), (uint32_t)p->info, (long)p->x[2]);
				if (p->x[2] <= min_iwidth) {
					for (k = 0; k < p->x[2]; ++k) {
						bwtint_t pos;
						int len, is_rev, ref_id;
						len  = (uint32_t)p->info - (p->info>>32);
						pos = bns_depos(idx->bns, bwt_sa(idx->bwt, p->x[0] + k), &is_rev);
						if (is_rev) pos -= len - 1;
						bns_cnt_ambi(idx->bns, pos, len, &ref_id);
						printf("\t%s:%c%ld", idx->bns->anns[ref_id].name, "+-"[is_rev], (long)(pos - idx->bns->anns[ref_id].offset) + 1);
					}
				} else fputs("\t*", stdout);
				putchar('\n');
			}
		}
		puts("//");
	}

	smem_itr_destroy(itr);
	bwa_idx_destroy(idx);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
