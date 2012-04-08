#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "bwa.h"
#include "bwt.h"
#include "bwtgap.h"
#include "bntseq.h"

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

extern unsigned char nst_nt4_table[256];
extern void seq_reverse(int len, uint8_t *seq, int is_comp);

bwa_opt_t bwa_def_opt = { 11, 4, -1, 1, 6, 32, 2, 0.04 };

struct bwa_idx_t {
	bwt_t *bwt;
	bntseq_t *bns;
	uint8_t *pac;
};

struct bwa_buf_t {
	int max_buf;
	bwa_pestat_t pes;
	gap_stack_t *stack;
	gap_opt_t *opt;
	int *diff_tab;
	uint8_t *buf;
	int *logn;
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

bwa_buf_t *bwa_buf_init(const bwa_opt_t *opt, int max_score)
{
	extern gap_opt_t *gap_init_opt(void);
	extern int bwa_cal_maxdiff(int l, double err, double thres);
	int i;
	bwa_buf_t *p;
	p = malloc(sizeof(bwa_buf_t));
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
	p->logn = calloc(256, sizeof(int));
	for (i = 1; i != 256; ++i)
		p->logn[i] = (int)(4.343 * log(i) + 0.499);
	return p;
}

void bwa_buf_destroy(bwa_buf_t *p)
{
	gap_destroy_stack(p->stack);
	free(p->diff_tab); free(p->logn); free(p->opt);
	free(p);
}

bwa_sai_t bwa_sai(const bwa_idx_t *idx, bwa_buf_t *buf, const char *seq)
{
	extern int bwt_cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width);
	int i, seq_len, buf_len;
	bwt_width_t *w, *seed_w;
	uint8_t *s;
	gap_opt_t opt2 = *buf->opt;
	bwa_sai_t sai;

	seq_len = strlen(seq);
	// estimate the buffer length
	buf_len = (buf->opt->seed_len + seq_len + 1) * sizeof(bwt_width_t) + seq_len;
	if (buf_len > buf->max_buf) {
		buf->max_buf = buf_len;
		kroundup32(buf->max_buf);
		buf->buf = realloc(buf->buf, buf->max_buf);
	}
	memset(buf->buf, 0, buf_len);
	seed_w = (bwt_width_t*)buf->buf;
	w = seed_w + buf->opt->seed_len;
	s = (uint8_t*)(w + seq_len + 1);
	if (opt2.fnr > 0.) opt2.max_diff = buf->diff_tab[seq_len];
	// copy the sequence
	for (i = 0; i < seq_len; ++i)
		s[i] = nst_nt4_table[(int)seq[i]];
	seq_reverse(seq_len, s, 0);
	// mapping
	bwt_cal_width(idx->bwt, seq_len, s, w);
	if (opt2.seed_len >= seq_len) opt2.seed_len = 0x7fffffff;
	if (seq_len > buf->opt->seed_len)
		bwt_cal_width(idx->bwt, buf->opt->seed_len, s + (seq_len - buf->opt->seed_len), seed_w);
	for (i = 0; i < seq_len; ++i) // complement; I forgot why...
		s[i] = s[i] > 3? 4 : 3 - s[i];
	sai.sai = (bwa_sai1_t*)bwt_match_gap(idx->bwt, seq_len, s, w, seq_len <= buf->opt->seed_len? 0 : seed_w, &opt2, &sai.n, buf->stack);
	return sai;
}

static void compute_NM(const uint8_t *pac, uint64_t l_pac, uint8_t *seq, int64_t pos, int n_cigar, uint32_t *cigar, int *n_mm, int *n_gaps)
{
	uint64_t x = pos, z;
	int k, y = 0;
	*n_mm = *n_gaps = 0;
	for (k = 0; k < n_cigar; ++k) {
		int l = cigar[k]>>4;
		int op = cigar[k]&0xf;
		if (op == 0) { // match/mismatch
			for (z = 0; z < l && x + z < l_pac; ++z) {
				int c = pac[(x+z)>>2] >> ((~(x+z)&3)<<1) & 3;
				if (c > 3 || seq[y+z] > 3 || c != seq[y+z]) ++(*n_mm);
			}
		}
		if (op == 1 || op == 2) (*n_gaps) += l;
		if (op == 0 || op == 2) x += l;
		if (op == 0 || op == 1 || op == 4) y += l;
	}
}

void bwa_sa2aln(const bwa_idx_t *idx, bwa_buf_t *buf, const char *seq, uint64_t sa, int n_gaps, bwa_aln_t *aln)
{
	extern bwtint_t bwa_sa2pos(const bntseq_t *bns, const bwt_t *bwt, bwtint_t sapos, int len, int *strand);
	extern bwa_cigar_t *bwa_refine_gapped_core(bwtint_t l_pac, const ubyte_t *pacseq, int len, const uint8_t *seq, bwtint_t *_pos, int ext, int *n_cigar, int is_end_correct);
	int strand, seq_len, i, n_gap, n_mm;
	uint64_t pos3, pac_pos;
	uint8_t *s[2];

	memset(aln, 0, sizeof(bwa_aln_t));
	seq_len = strlen(seq);
	if (seq_len<<1 > buf->max_buf) {
		buf->max_buf = seq_len<<1;
		kroundup32(buf->max_buf);
		buf->buf = realloc(buf->buf, buf->max_buf);
	}
	s[0] = buf->buf;
	s[1] = s[0] + seq_len;
	for (i = 0; i < seq_len; ++i)
		s[0][i] = s[1][i] = nst_nt4_table[(int)seq[i]];
	seq_reverse(seq_len, s[1], 1);
	pac_pos = bwa_sa2pos(idx->bns, idx->bwt, sa, seq_len, &strand);
	if (strand) aln->flag |= 16;
	if (n_gaps) { // only for gapped alignment
		int n_cigar;
		bwa_cigar_t *cigar16;
		cigar16 = bwa_refine_gapped_core(idx->bns->l_pac, idx->pac, seq_len, s[strand], &pac_pos, strand? n_gaps : -n_gaps, &n_cigar, 1);
		aln->n_cigar = n_cigar;
		aln->cigar = malloc(n_cigar * 4);
		for (i = 0, pos3 = pac_pos; i < n_cigar; ++i) {
			int op = cigar16[i]>>14;
			int len = cigar16[i]&0x3fff;
			if (op == 3) op = 4; // the 16-bit CIGAR is different from the 32-bit CIGAR
			aln->cigar[i] = len<<4 | op;
			if (op == 0 || op == 2) pos3 += len;
		}
		free(cigar16);
	} else { // ungapped
		aln->n_cigar = 1;
		aln->cigar = malloc(4);
		aln->cigar[0] = seq_len<<4 | 0;
		pos3 = pac_pos + seq_len;
	}
	aln->n_n = bns_cnt_ambi(idx->bns, pac_pos, pos3 - pac_pos, &aln->ref_id);
	aln->offset = pac_pos - idx->bns->anns[aln->ref_id].offset;
	if (pos3 - idx->bns->anns[aln->ref_id].offset > idx->bns->anns[aln->ref_id].len) // read mapped beyond the end of a sequence
		aln->flag |= 4; // read unmapped
	compute_NM(idx->pac, idx->bns->l_pac, s[strand], pac_pos, aln->n_cigar, aln->cigar, &n_mm, &n_gap);
	aln->n_mm = n_mm;
	aln->n_gap = n_gap;
}

/************************
 * Single-end alignment *
 ************************/

bwa_one_t *bwa_se(const bwa_idx_t *idx, bwa_buf_t *buf, const char *seq, int gen_cigar)
{
	bwa_one_t *one;
	int best, cnt, i, seq_len;

	seq_len = strlen(seq);
	one = calloc(1, sizeof(bwa_one_t));
	one->sai = bwa_sai(idx, buf, seq);
	if (one->sai.n == 0) return one;
	// count number of hits; randomly select one alignment
	best = one->sai.sai[0].score;
	for (i = cnt = 0; i < one->sai.n; ++i) {
		bwa_sai1_t *p = &one->sai.sai[i];
		if (p->score > best) break;
		if (drand48() * (p->l - p->k + 1 + cnt) > (double)cnt) {
			one->which = p;
			one->sa = p->k + (bwtint_t)((p->l - p->k + 1) * drand48());
		}
		cnt += p->l - p->k + 1;
	}
	one->c1 = cnt;
	for (; i < one->sai.n; ++i)
		cnt += one->sai.sai[i].l - one->sai.sai[i].k + 1;
	one->c2 = cnt - one->c1;
	// estimate single-end mapping quality
	one->mapQs = -1;
	if (one->c1 == 0) one->mapQs = 23; // FIXME: is it possible?
	else if (one->c1 > 1) one->mapQs = 0;
	else {
		int diff = one->which->n_mm + one->which->n_gapo + one->which->n_gape;
		if (diff >= buf->diff_tab[seq_len]) one->mapQs = 25;
		else if (one->c2 == 0) one->mapQs = 37;
	}
	if (one->mapQs < 0) {
		cnt = (one->c2 >= 255)? 255 : one->c2;
		one->mapQs = 23 < buf->logn[cnt]? 0 : 23 - buf->logn[cnt];
	}
	one->mapQ = one->mapQs;
	// compute CIGAR on request
	one->one.ref_id = -1;
	if (gen_cigar) bwa_sa2aln(idx, buf, seq, one->sa, one->which->n_gapo + one->which->n_gape, &one->one);
	return one;
}

void bwa_one_destroy(bwa_one_t *one)
{
	free(one->sai.sai);
	free(one->one.cigar);
	free(one);
}

/************************
 * Paired-end alignment *
 ************************/

void bwa_pestat(bwa_buf_t *buf, int n, bwa_one_t **o[2])
{
}

void bwa_pe(const bwa_idx_t *idx, bwa_buf_t *buf, const char *seq[2], bwa_one_t *o[2])
{
}
