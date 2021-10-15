#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <unistd.h>
#include "bwa.h"
#include "bwamem.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

int main_maxk(int argc, char *argv[])
{
	int i, c, self = 0, max_len = 0;
	uint8_t *cnt = 0;
	uint64_t hist[256];
	bwt_t *bwt;
	kseq_t *ks;
	smem_i *itr;
	gzFile fp;

	while ((c = getopt(argc, argv, "s")) >= 0) {
		if (c == 's') self = 1;
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: bwa maxk [-s] <index.prefix> <seq.fa>\n");
		return 1;
	}
	fp = strcmp(argv[optind+1], "-")? gzopen(argv[optind+1], "rb") : gzdopen(fileno(stdin), "rb");
	ks = kseq_init(fp);
	bwt = bwt_restore_bwt(argv[optind]);
	itr = smem_itr_init(bwt);
	if (self) smem_config(itr, 2, INT_MAX, 0);
	memset(hist, 0, 8 * 256);

	while (kseq_read(ks) >= 0) {
		const bwtintv_v *a;
		if (ks->seq.l > max_len) {
			max_len = ks->seq.l;
			kroundup32(max_len);
			cnt = realloc(cnt, max_len);
		}
		memset(cnt, 0, ks->seq.l);
		for (i = 0; i < ks->seq.l; ++i)
			ks->seq.s[i] = nst_nt4_table[(int)ks->seq.s[i]];
		smem_set_query(itr, ks->seq.l, (uint8_t*)ks->seq.s);
		while ((a = smem_next(itr)) != 0) {
			for (i = 0; i < a->n; ++i) {
				bwtintv_t *p = &a->a[i];
				int j, l, start = p->info>>32, end = (uint32_t)p->info;
				l = end - start < 255? end - start : 255;
				for (j = start; j < end; ++j)
					cnt[j] = cnt[j] > l? cnt[j] : l;
			}
		}
		for (i = 0; i < ks->seq.l; ++i) ++hist[cnt[i]];
	}
	for (i = 0; i < 256; ++i)
		printf("%d\t%lld\n", i, (long long)hist[i]);
	free(cnt);

	smem_itr_destroy(itr);
	bwt_destroy(bwt);
	kseq_destroy(ks);
	gzclose(fp);
	return 0;
}
