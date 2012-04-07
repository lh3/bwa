#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "bwa.h"
#include "bwt.h"
#include "bwtgap.h"
#include "bntseq.h"

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

extern unsigned char nst_nt4_table[256];

bwa_opt_t bwa_def_opt = { 11, 4, -1, 1, 6, 32, 2, 0.04 };

struct bwa_idx_t {
	bwt_t *bwt;
	bntseq_t *bns;
	uint8_t *pac;
};

struct bwa_aux_t {
	int max_buf;
	gap_stack_t *stack;
	gap_opt_t *opt;
	int *diff_tab;
	uint8_t *buf;
};

bwa_idx_t *bwa_idx_load(const char *prefix)
{
	bwa_idx_t *p;
	int l;
	char *str;
	l = strlen(prefix);
	p = calloc(1, sizeof(bwa_idx_t));
	str = malloc(l + 10);
	strcpy(str, prefix);
	p->bns = bns_restore(str);
	strcpy(str + l, ".bwt");
	p->bwt = bwt_restore_bwt(str);
	str[l] = 0;
	strcpy(str + l, ".sa");
	bwt_restore_sa(str, p->bwt);
	free(str);
	p->pac = calloc(p->bns->l_pac/4+1, 1);
	fread(p->pac, 1, p->bns->l_pac/4+1, p->bns->fp_pac);
	fclose(p->bns->fp_pac);
	p->bns->fp_pac = 0;
	return p;
}

void bwa_idx_destroy(bwa_idx_t *p)
{
	bns_destroy(p->bns);
	bwt_destroy(p->bwt);
	free(p->pac);
	free(p);
}

bwa_aux_t *bwa_aux_init(const bwa_opt_t *opt, int max_score)
{
	extern gap_opt_t *gap_init_opt(void);
	extern int bwa_cal_maxdiff(int l, double err, double thres);
	int i;
	bwa_aux_t *p;
	p = malloc(sizeof(bwa_aux_t));
	p->stack = gap_init_stack2(max_score);
	p->opt = gap_init_opt();
	p->opt->s_gapo = opt->s_gapo;
	p->opt->s_gape = opt->s_gape;
	p->opt->max_diff = opt->max_diff;
	p->opt->max_gapo = opt->max_gapo;
	p->opt->max_gape = opt->max_gape;
	p->opt->seed_len = opt->seed_len;
	p->opt->max_seed_diff = opt->max_seed_diff;
	p->opt->fnr = opt->fnr;
	p->diff_tab = calloc(BWA_MAX_QUERY_LEN, sizeof(int));
	for (i = 1; i < BWA_MAX_QUERY_LEN; ++i)
		p->diff_tab[i] = bwa_cal_maxdiff(i, BWA_AVG_ERR, opt->fnr);
	return p;
}

void bwa_aux_destroy(bwa_aux_t *p)
{
	gap_destroy_stack(p->stack);
	free(p->diff_tab); free(p->opt);
	free(p);
}

bwa_alnpre_t *bwa_aln_pre(const bwa_idx_t *idx, bwa_aux_t *aux, const char *seq, int *n_aln)
{
	extern int bwt_cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width);
	extern void seq_reverse(int len, uint8_t *seq, int is_comp);
	int i, seq_len, buf_len;
	bwt_width_t *w, *seed_w;
	uint8_t *s;
	gap_opt_t opt2 = *aux->opt;

	seq_len = strlen(seq);
	// estimate the buffer length
	buf_len = (aux->opt->seed_len + seq_len + 1) * sizeof(bwt_width_t) + seq_len;
	if (buf_len > aux->max_buf) {
		aux->max_buf = buf_len;
		kroundup32(aux->max_buf);
		aux->buf = realloc(aux->buf, aux->max_buf);
	}
	memset(aux->buf, 0, buf_len);
	seed_w = (bwt_width_t*)aux->buf;
	w = seed_w + aux->opt->seed_len;
	s = (uint8_t*)(w + seq_len + 1);
	if (opt2.fnr > 0.) opt2.max_diff = aux->diff_tab[seq_len];
	// copy the sequence
	for (i = 0; i < seq_len; ++i)
		s[i] = nst_nt4_table[(int)seq[i]];
	seq_reverse(seq_len, s, 0);
	// mapping
	bwt_cal_width(idx->bwt, seq_len, s, w);
	if (opt2.seed_len >= seq_len) opt2.seed_len = 0x7fffffff;
	if (seq_len > aux->opt->seed_len)
		bwt_cal_width(idx->bwt, aux->opt->seed_len, s + (seq_len - aux->opt->seed_len), seed_w);
	for (i = 0; i < seq_len; ++i) // complement; I forgot why...
		s[i] = s[i] > 3? 4 : 3 - s[i];
	return (bwa_alnpre_t*)bwt_match_gap(idx->bwt, seq_len, s, w, seq_len <= aux->opt->seed_len? 0 : seed_w, &opt2, n_aln, aux->stack);
}
