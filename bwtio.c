#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "bwt.h"
#include "utils.h"

void bwt_dump_bwt(const char *fn, const bwt_t *bwt)
{
	FILE *fp = NULL;
	fp = xopen(fn, "wb");
	err_fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fwrite(bwt->L2+1, sizeof(bwtint_t), 4, fp);
	err_fwrite(bwt->bwt, 4, bwt->bwt_size, fp);
	err_fclose(fp);
}

void bwt_dump_sa(const char *fn, const bwt_t *bwt)
{
	FILE *fp;
	fp = xopen(fn, "wb");
	err_fwrite(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fwrite(bwt->L2+1, sizeof(bwtint_t), 4, fp);
	err_fwrite(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	err_fwrite(&bwt->seq_len, sizeof(bwtint_t), 1, fp);
	err_fwrite(bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp);
	err_fclose(fp);
}

void bwt_restore_sa(const char *fn, bwt_t *bwt)
{
	char skipped[256];
	FILE *fp;
	bwtint_t primary;

	fp = xopen(fn, "rb");
	err_fread_noeof(&primary, sizeof(bwtint_t), 1, fp);
	xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
	err_fread_noeof(skipped, sizeof(bwtint_t), 4, fp); // skip
	err_fread_noeof(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	err_fread_noeof(&primary, sizeof(bwtint_t), 1, fp);
	xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");

	bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = (bwtint_t*)xcalloc(bwt->n_sa, sizeof(bwtint_t));
	bwt->sa[0] = -1;

	err_fread_noeof(bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp);
	err_fclose(fp);
}

bwt_t *bwt_restore_bwt(const char *fn)
{
	bwt_t *bwt;
	FILE *fp;

	bwt = (bwt_t*)xcalloc(1, sizeof(bwt_t));
	fp = xopen(fn, "rb");
	err_fseek(fp, 0, SEEK_END);
	bwt->bwt_size = (err_ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
	bwt->bwt = (uint32_t*)xcalloc(bwt->bwt_size, 4);
	err_fseek(fp, 0, SEEK_SET);
	err_fread_noeof(&bwt->primary, sizeof(bwtint_t), 1, fp);
	err_fread_noeof(bwt->L2+1, sizeof(bwtint_t), 4, fp);
	err_fread_noeof(bwt->bwt, 4, bwt->bwt_size, fp);
	bwt->seq_len = bwt->L2[4];
	err_fclose(fp);
	bwt_gen_cnt_table(bwt);

	return bwt;
}

void bwt_destroy(bwt_t *bwt)
{
	if (bwt == 0) return;
	free(bwt->sa); free(bwt->bwt);
	free(bwt);
}
