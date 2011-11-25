#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "bntseq.h"
#include "bwt.h"
#include "kvec.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

extern unsigned char nst_nt4_table[256];

typedef struct {
	const bwt_t *bwt;
	const uint8_t *query;
	int start, len;
	bwtintv_v *tmpvec[2], *matches;
} smem_i;

smem_i *smem_iter_init(const bwt_t *bwt)
{
	smem_i *iter;
	iter = calloc(1, sizeof(smem_i));
	iter->bwt = bwt;
	iter->tmpvec[0] = calloc(1, sizeof(bwtintv_v));
	iter->tmpvec[1] = calloc(1, sizeof(bwtintv_v));
	iter->matches   = calloc(1, sizeof(bwtintv_v));
	return iter;
}

void smem_iter_destroy(smem_i *iter)
{
	free(iter->tmpvec[0]->a);
	free(iter->tmpvec[1]->a);
	free(iter->matches->a);
	free(iter);
}

void smem_set_query(smem_i *iter, int len, const uint8_t *query)
{
	iter->query = query;
	iter->start = 0;
	iter->len = len;
}

int smem_next(smem_i *iter)
{
	iter->tmpvec[0]->n = iter->tmpvec[1]->n = iter->matches->n = 0;
	if (iter->start >= iter->len || iter->start < 0) return -1;
	while (iter->start < iter->len && iter->query[iter->start] > 3) ++iter->start; // skip ambiguous bases
	if (iter->start == iter->len) return -1;
	iter->start = bwt_smem1(iter->bwt, iter->len, iter->query, iter->start, iter->matches, iter->tmpvec);
	return iter->start;
}

int main_fastmap(int argc, char *argv[])
{
	int c, i, min_iwidth = 20, min_len = 17, print_seq = 0;
	kseq_t *seq;
	bwtint_t k;
	gzFile fp;
	bwt_t *bwt;
	bntseq_t *bns;
	smem_i *iter;

	while ((c = getopt(argc, argv, "w:l:s")) >= 0) {
		switch (c) {
			case 's': print_seq = 1; break;
			case 'w': min_iwidth = atoi(optarg); break;
			case 'l': min_len = atoi(optarg); break;
		}
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "Usage: bwa fastmap [-s] [-l minLen=%d] [-w maxSaSize=%d] <idxbase> <in.fq>\n", min_len, min_iwidth);
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
	iter = smem_iter_init(bwt);
	while (kseq_read(seq) >= 0) {
		printf("SQ\t%s\t%ld", seq->name.s, seq->seq.l);
		if (print_seq) {
			putchar('\t');
			puts(seq->seq.s);
		} else putchar('\n');
		for (i = 0; i < seq->seq.l; ++i)
			seq->seq.s[i] = nst_nt4_table[(int)seq->seq.s[i]];
		smem_set_query(iter, seq->seq.l, (uint8_t*)seq->seq.s);
		while (smem_next(iter) > 0) {
			for (i = 0; i < iter->matches->n; ++i) {
				bwtintv_t *p = &iter->matches->a[i];
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

	smem_iter_destroy(iter);
	bns_destroy(bns);
	bwt_destroy(bwt);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
