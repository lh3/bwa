#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "bntseq.h"
#include "bwt.h"
#include "bwamem.h"
#include "kvec.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

extern unsigned char nst_nt4_table[256];

int main_mem(int argc, char *argv[])
{
	mem_opt_t *opt;
	bwt_t *bwt;
	bntseq_t *bns;
	int i, j, c;
	gzFile fp;
	kseq_t *seq;
	uint8_t *pac = 0;

	opt = mem_opt_init();
	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]\n");
		fprintf(stderr, "\n");
		free(opt);
		return 1;
	}
	mem_fill_scmat(opt->a, opt->b, opt->mat);
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
		pac = calloc(bns->l_pac/4+1, 1);
		fread(pac, 1, bns->l_pac/4+1, bns->fp_pac);
	}
	while (kseq_read(seq) >= 0) {
		mem_chain_t chain;
		printf(">%s\n", seq->name.s);
		for (i = 0; i < seq->seq.l; ++i)
			seq->seq.s[i] = nst_nt4_table[(int)seq->seq.s[i]];
		chain = mem_chain(opt, bwt, seq->seq.l, (uint8_t*)seq->seq.s);
		for (i = 0; i < chain.n; ++i) {
			mem_chain1_t *p = &chain.chains[i];
			mem_aln_t a;
			mem_chain2aln(opt, bns->l_pac, pac, seq->seq.l, (uint8_t*)seq->seq.s, p, &a);
			printf("%d\t%d", i, p->n);
			for (j = 0; j < p->n; ++j) {
				bwtint_t pos;
				int is_rev, ref_id;
				pos = bns_depos(bns, p->seeds[j].rbeg, &is_rev);
				if (is_rev) pos -= p->seeds[j].len - 1;
				bns_cnt_ambi(bns, pos, p->seeds[j].len, &ref_id);
				printf("\t%d,%d,%s:%c%ld", p->seeds[j].len, p->seeds[j].qbeg, bns->anns[ref_id].name, "+-"[is_rev], (long)(pos - bns->anns[ref_id].offset) + 1);
			}
			putchar('\n');
		}
		puts("//");
		for (i = 0; i < chain.n; ++i) free(chain.chains[i].seeds);
		free(chain.chains);
	}

	free(pac); free(opt);
	bns_destroy(bns);
	bwt_destroy(bwt);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

int main_fastmap(int argc, char *argv[])
{
	int c, i, min_iwidth = 20, min_len = 17, print_seq = 0, split_long = 0;
	kseq_t *seq;
	bwtint_t k;
	gzFile fp;
	bwt_t *bwt;
	bntseq_t *bns;
	smem_i *itr;
	const bwtintv_v *a;

	while ((c = getopt(argc, argv, "w:l:ps")) >= 0) {
		switch (c) {
			case 's': split_long = 1; break;
			case 'p': print_seq = 1; break;
			case 'w': min_iwidth = atoi(optarg); break;
			case 'l': min_len = atoi(optarg); break;
		}
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "Usage: bwa fastmap [-ps] [-l minLen=%d] [-w maxSaSize=%d] <idxbase> <in.fq>\n", min_len, min_iwidth);
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
		while ((a = smem_next(itr, split_long? min_len<<1 : 0)) != 0) {
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
