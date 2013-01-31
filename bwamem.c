#include <stdlib.h>
#include "bwamem.h"

memopt_t *mem_opt_init()
{
	memopt_t *o;
	o = calloc(1, sizeof(memopt_t));
	o->a = 1; o->b = 9; o->q = 16; o->r = 1; o->w = 100;
	o->min_seed_len = 17;
	o->max_occ = 10;
	return o;
}

/***************************
 * SMEM iterator interface *
 ***************************/

smem_i *smem_itr_init(const bwt_t *bwt)
{
	smem_i *itr;
	itr = calloc(1, sizeof(smem_i));
	itr->bwt = bwt;
	itr->tmpvec[0] = calloc(1, sizeof(bwtintv_v));
	itr->tmpvec[1] = calloc(1, sizeof(bwtintv_v));
	itr->matches   = calloc(1, sizeof(bwtintv_v));
	return itr;
}

void smem_itr_destroy(smem_i *itr)
{
	free(itr->tmpvec[0]->a);
	free(itr->tmpvec[1]->a);
	free(itr->matches->a);
	free(itr);
}

void smem_set_query(smem_i *itr, int min_intv, int len, const uint8_t *query)
{
	itr->query = query;
	itr->start = 0;
	itr->len = len;
	itr->min_intv = min_intv;
}

int smem_next(smem_i *itr)
{
	itr->tmpvec[0]->n = itr->tmpvec[1]->n = itr->matches->n = 0;
	if (itr->start >= itr->len || itr->start < 0) return -1;
	while (itr->start < itr->len && itr->query[itr->start] > 3) ++itr->start; // skip ambiguous bases
	if (itr->start == itr->len) return -1;
	itr->start = bwt_smem1(itr->bwt, itr->len, itr->query, itr->start, itr->min_intv, itr->matches, itr->tmpvec);
	return itr->start;
}

