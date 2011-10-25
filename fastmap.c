#include <zlib.h>
#include <unistd.h>
#include <stdio.h>
#include "bwt.h"
#include "kvec.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

extern unsigned char nst_nt4_table[256];

int main_fastmap(int argc, char *argv[])
{
	int c, i;
	kseq_t *seq;
	gzFile fp;
	bwt_t *bwt;
	bwtintv_v a[3], mem, *tvec[3];
	while ((c = getopt(argc, argv, "")) >= 0) {
		switch (c) {
		}
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "bwa fastmap <idxbase> <in.fq>\n");
		return 1;
	}
	fp = gzopen(argv[optind + 1], "r");
	seq = kseq_init(fp);
	{ // load the BWT
		char *tmp = calloc(strlen(argv[optind]) + 5, 1);
		strcat(strcpy(tmp, argv[optind]), ".bwt");
		bwt = bwt_restore_bwt(tmp);
		free(tmp);
	}
	for (i = 0; i < 3; ++i) {
		kv_init(a[i]);
		tvec[i] = &a[i];
	}
	kv_init(mem);
	while (kseq_read(seq) >= 0) {
		for (i = 0; i < seq->seq.l; ++i)
			seq->seq.s[i] = nst_nt4_table[(int)seq->seq.s[i]];
		bwt_smem(bwt, seq->seq.l, (uint8_t*)seq->seq.s, &mem, tvec);
		printf(">%s\t%ld\n", seq->name.s, mem.n);
		for (i = 0; i < mem.n; ++i) {
			bwtintv_t *p = &mem.a[i];
			printf("%d\t%d\t%ld\n", (uint32_t)(p->info>>32), (uint32_t)p->info, (long)p->x[2]);
		}
		puts("//");
	}
	free(mem.a);
	for (i = 0; i < 3; ++i) free(a[i].a);
	bwt_destroy(bwt);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
