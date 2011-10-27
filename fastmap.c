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

int main_fastmap(int argc, char *argv[])
{
	int c, i, min_iwidth = 20, min_len = 17;
	kseq_t *seq;
	bwtint_t k;
	gzFile fp;
	bwt_t *bwt;
	bntseq_t *bns;
	bwtintv_v a[3], mem, *tvec[3];

	while ((c = getopt(argc, argv, "w:l:")) >= 0) {
		switch (c) {
			case 'w': min_iwidth = atoi(optarg); break;
			case 'l': min_len = atoi(optarg); break;
		}
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "Usage: bwa fastmap [-l minLen=%d] [-w maxSaSize=%d] <idxbase> <in.fq>\n", min_len, min_iwidth);
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
	for (i = 0; i < 3; ++i) { // initiate the temporary array
		kv_init(a[i]);
		tvec[i] = &a[i];
	}
	kv_init(mem);
	while (kseq_read(seq) >= 0) {
		for (i = 0; i < seq->seq.l; ++i)
			seq->seq.s[i] = nst_nt4_table[(int)seq->seq.s[i]];
		bwt_smem(bwt, seq->seq.l, (uint8_t*)seq->seq.s, &mem, tvec);
		printf("SQ\t%s\t%ld\n", seq->name.s, seq->seq.l);
		for (i = 0; i < mem.n; ++i) {
			bwtintv_t *p = &mem.a[i];
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
			}
			putchar('\n');
		}
		puts("//");
	}

	free(mem.a);
	for (i = 0; i < 3; ++i) free(a[i].a);
	bns_destroy(bns);
	bwt_destroy(bwt);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
