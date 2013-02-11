#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "bntseq.h"
#include "bwt.h"
#include "bwamem.h"
#include "kvec.h"
#include "bseq.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

extern unsigned char nst_nt4_table[256];

int main_mem(int argc, char *argv[])
{
	mem_opt_t *opt;
	bwt_t *bwt;
	bntseq_t *bns;
	int c, n;
	gzFile fp, fp2 = 0;
	kseq_t *ks, *ks2 = 0;
	uint8_t *pac = 0;
	bseq1_t *seqs;

	opt = mem_opt_init();
	while ((c = getopt(argc, argv, "k:c:v:s:")) >= 0) {
		if (c == 'k') opt->min_seed_len = atoi(optarg);
		else if (c == 'c') opt->max_occ = atoi(optarg);
		else if (c == 'v') mem_verbose = atoi(optarg);
		else if (c == 's') opt->split_width = atoi(optarg);
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa mem [options] <idxbase> <in.fq>\n\n");
		fprintf(stderr, "Options: -k INT     minimum seed length [%d]\n", opt->min_seed_len);
		fprintf(stderr, "         -c INT     skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
		fprintf(stderr, "         -s INT     look for internal seeds inside a seed with less than INT occ [%d]\n", opt->split_width);
		fprintf(stderr, "         -v INT     verbose level [%d]\n", mem_verbose);
		fprintf(stderr, "\n");
		free(opt);
		return 1;
	}
	mem_fill_scmat(opt->a, opt->b, opt->mat);
	{ // load the packed sequences, BWT and SA
		char *tmp = calloc(strlen(argv[optind]) + 5, 1);
		strcat(strcpy(tmp, argv[optind]), ".bwt");
		bwt = bwt_restore_bwt(tmp);
		strcat(strcpy(tmp, argv[optind]), ".sa");
		bwt_restore_sa(tmp, bwt);
		free(tmp);
		bns = bns_restore(argv[optind]);
		pac = calloc(bns->l_pac/4+1, 1);
		fread(pac, 1, bns->l_pac/4+1, bns->fp_pac);
	}

	fp = strcmp(argv[optind + 1], "-")? gzopen(argv[optind + 1], "r") : gzdopen(fileno(stdin), "r");
	ks = kseq_init(fp);
	if (optind + 2 < argc) {
		fp2 = gzopen(argv[optind + 2], "r");
		ks2 = kseq_init(fp2);
		opt->is_pe = 1;
	}
	while ((seqs = bseq_read(opt->n_threads * opt->chunk_size, &n, ks, ks2)) != 0) {
		mem_process_seqs(opt, bwt, bns, pac, n, seqs);
		free(seqs);
	}

	free(opt); free(pac);
	bns_destroy(bns);
	bwt_destroy(bwt);
	kseq_destroy(ks);
	gzclose(fp);
	if (ks2) {
		kseq_destroy(ks2);
		gzclose(fp2);
	}
	return 0;
}

int main_fastmap(int argc, char *argv[])
{
	int c, i, min_iwidth = 20, min_len = 17, print_seq = 0, split_width = 0;
	kseq_t *seq;
	bwtint_t k;
	gzFile fp;
	bwt_t *bwt;
	bntseq_t *bns;
	smem_i *itr;
	const bwtintv_v *a;

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
	{ // load the packed sequences, BWT and SA
		char *tmp = calloc(strlen(argv[optind]) + 5, 1);
		strcat(strcpy(tmp, argv[optind]), ".bwt");
		bwt = bwt_restore_bwt(tmp);
		strcat(strcpy(tmp, argv[optind]), ".sa");
		bwt_restore_sa(tmp, bwt);
		free(tmp);
		bns = bns_restore(argv[optind]);
	}
	itr = smem_itr_init(bwt);
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
						pos = bns_depos(bns, bwt_sa(bwt, p->x[0] + k), &is_rev);
						if (is_rev) pos -= len - 1;
						bns_cnt_ambi(bns, pos, len, &ref_id);
						printf("\t%s:%c%ld", bns->anns[ref_id].name, "+-"[is_rev], (long)(pos - bns->anns[ref_id].offset) + 1);
					}
				} else fputs("\t*", stdout);
				putchar('\n');
			}
		}
		puts("//");
	}

	smem_itr_destroy(itr);
	bns_destroy(bns);
	bwt_destroy(bwt);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
