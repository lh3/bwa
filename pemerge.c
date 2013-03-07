#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <zlib.h>
#include "ksw.h"
#include "kseq.h"

#ifdef _PEM_MAIN
KSEQ_INIT(gzFile, gzread)

unsigned char nst_nt4_table[128] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void bwa_fill_scmat(int a, int b, int8_t mat[25])
{
	int i, j, k;
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? a : -b;
		mat[k++] = 0; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = 0;
}
#else
#include "bwa.h"
KSEQ_DECLARE(gzFile)
#endif

#define MAX_SCORE_RATIO 0.9f
#define MAX_ERR 8

static const char *err_msg[MAX_ERR+1] = {
	"successful merges",
	"low-scoring pairs",
	"pairs where the best SW alignment is not an overlap (long left end)",
	"pairs where the best SW alignment is not an overlap (long right end)",
	"pairs with large 2nd best SW score",
	"pairs with gapped overlap",
	"pairs where the end-to-end ungapped alignment score is not high enough",
	"pairs potentially with tandem overlaps",
	"pairs with high sum of errors"
};

typedef struct {
	kstring_t n, s, q;
} mseq_t;

typedef struct {
	int a, b, q, r, w;
	int q_def, q_thres;
	int T;
	int8_t mat[25];
} pem_opt_t;

pem_opt_t *pem_opt_init()
{
	pem_opt_t *opt;
	opt = calloc(1, sizeof(pem_opt_t));
	opt->a = 5; opt->b = 4; opt->q = 2, opt->r = 17; opt->w = 20;
	opt->T = opt->a * 10;
	opt->q_def = 20;
	opt->q_thres = 70;
	bwa_fill_scmat(opt->a, opt->b, opt->mat);
	return opt;
}

int bwa_pemerge(const pem_opt_t *opt, mseq_t x[2], uint8_t **seq_, uint8_t **qual_)
{
	uint8_t *s[2], *q[2], *seq, *qual;
	int i, xtra, l, l_seq, sum_q;
	kswr_t r;

	*seq_ = *qual_ = 0;
	s[0] = malloc(x[0].s.l); q[0] = malloc(x[0].s.l);
	s[1] = malloc(x[1].s.l); q[1] = malloc(x[1].s.l);
	for (i = 0; i < x[0].s.l; ++i) {
		int c = x[0].s.s[i];
		s[0][i] = c < 0 || c > 127? 4 : c <= 4? c : nst_nt4_table[c];
		q[0][i] = x[0].q.l? x[0].q.s[i] - 33 : opt->q_def;
	}
	for (i = 0; i < x[1].s.l; ++i) {
		int c = x[1].s.s[x[1].s.l - 1 - i];
		c = c < 0 || c > 127? 4 : c < 4? c : nst_nt4_table[c];
		s[1][i] = c < 4? 3 - c : 4;
		q[1][i] = x[1].q.l? x[1].q.s[x[1].q.l - 1 - i] - 33 : opt->q_def;
	}

	xtra = KSW_XSTART | KSW_XSUBO;
	r = ksw_align(x[1].s.l, s[1], x[0].s.l, s[0], 5, opt->mat, opt->q, opt->r, xtra, 0);
	++r.qe; ++r.te; // change to the half-close-half-open coordinates

	if (r.score < opt->T) return -1; // poor alignment
	if (r.tb < r.qb) return -2; // no enough space for the left end
	if (x[0].s.l - r.te > x[1].s.l - r.qe) return -3; // no enough space for the right end
	if ((double)r.score2 / r.score >= MAX_SCORE_RATIO) return -4; // the second best score is too large
	if (r.qe - r.qb != r.te - r.tb) return -5; // we do not allow gaps

	{ // test tandem match; O(n^2)
		int max_m, max_m2, min_l, max_l, max_l2;
		max_m = max_m2 = 0; max_l = max_l2 = 0;
		min_l = x[0].s.l < x[1].s.l? x[0].s.l : x[1].s.l;
		for (l = 0; l < min_l; ++l) {
			int m = 0, o = x[0].s.l - l;
			int a = opt->a, b = -opt->b;
			for (i = 0; i < l; ++i)
				m += (s[0][o + i] == s[1][i])? a : b;
			if (m > max_m) max_m2 = max_m, max_m = m, max_l2 = max_l, max_l = l;
			else if (m > max_m2) max_m2 = m, max_l2 = l;
		}
		if (max_m < opt->T) return -6;
		if (max_l2 < max_l && max_m2 >= opt->T && (double)(max_m2 + (max_l - max_l2) * opt->a) / max_m >= MAX_SCORE_RATIO)
			return -7;
		//printf("*** %d,%d; %d,%d\n", max_m, max_m2, max_l, max_l2);
	}

	l = x[0].s.l - (r.tb - r.qb); // length to merge
	l_seq = x[0].s.l + x[1].s.l - l;
	seq = malloc(l_seq + 1);
	qual = malloc(l_seq + 1);
	memcpy(seq,  s[0], x[0].s.l); memcpy(seq  + x[0].s.l, &s[1][l], x[1].s.l - l);
	memcpy(qual, q[0], x[0].s.l); memcpy(qual + x[0].s.l, &q[1][l], x[1].s.l - l);
	for (i = 0, sum_q = 0; i < l; ++i) {
		int k = x[0].s.l - l + i;
		if (s[0][k] == 4) { // ambiguous
			seq[k]  = s[1][i];
			qual[k] = q[1][i];
		} else if (s[1][i] == 4) { // do nothing
		} else if (s[0][k] == s[1][i]) {
			qual[k] = qual[k] > q[1][i]? qual[k] : q[1][i];
		} else { // s[0][k] != s[1][i] and neither is N
			int qq = q[0][k] < q[1][i]? q[0][k] : q[1][i];
			sum_q += qq >= 3? qq<<1 : 1;
			seq[k]  = q[0][k] > q[1][i]? s[0][k] : s[1][i];
			qual[k] = abs((int)q[0][k] - (int)q[1][i]);
		}
	}
	if (sum_q>>1 > opt->q_thres) { // too many mismatches
		free(seq); free(qual);
		return -8;
	}

	for (i = 0; i < l_seq; ++i) seq[i] = "ACGTN"[(int)seq[i]], qual[i] += 33;
	seq[l_seq] = qual[l_seq] = 0;
	*seq_ = seq, *qual_ = qual;
	return l_seq;
}

static inline void kstrcpy(kstring_t *dst, const kstring_t *src)
{
	dst->l = 0;
	if (src->l == 0) return;
	if (dst->m < src->l + 1) {
		dst->m = src->l + 2;
		kroundup32(dst->m);
		dst->s = realloc(dst->s, dst->m);
	}
	dst->l = src->l;
	memcpy(dst->s, src->s, src->l + 1);
}

static inline void kseq2mseq(mseq_t *ms, const kseq_t *ks)
{
	kstrcpy(&ms->n, &ks->name);
	kstrcpy(&ms->s, &ks->seq);
	kstrcpy(&ms->q, &ks->qual);
}

static inline void print_seq(const char *n, const char *s, const char *q)
{
	putchar(q? '@' : '>');
	puts(n); puts(s);
	if (q) {
		puts("+"); puts(q);
	}
}

#ifdef _PEM_MAIN
int main(int argc, char *argv[])
#else
int main_pemerge(int argc, char *argv[])
#endif
{
	int c, flag = 0, i;
	int64_t cnt[MAX_ERR+1];
	gzFile fp, fp2 = 0;
	kseq_t *ks, *ks2 = 0;
	mseq_t r[2];
	pem_opt_t *opt;

	opt = pem_opt_init();
	while ((c = getopt(argc, argv, "muQ:")) >= 0) {
		if (c == 'm') flag |= 1;
		else if (c == 'u') flag |= 2;
		else if (c == 'Q') opt->q_thres = atoi(optarg);
	}
	if (flag == 0) flag = 3;

	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa pemerge [-mu] <read1.fq> [read2.fq]\n\n");
		fprintf(stderr, "Options: -m       output merged reads only\n");
		fprintf(stderr, "         -u       output unmerged reads only\n");
		fprintf(stderr, "         -Q INT   max sum of errors [%d]\n", opt->q_thres);
		fprintf(stderr, "\n");
		free(opt);
		return 1;
	}

	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	ks = kseq_init(fp);
	if (optind + 1 < argc) {
		fp2 = strcmp(argv[optind+1], "-")? gzopen(argv[optind+1], "r") : gzdopen(fileno(stdin), "r");
		ks2 = kseq_init(fp2);
	}

	memset(r, 0, sizeof(mseq_t)<<1);
	memset(cnt, 0, 8 * (MAX_ERR+1));
	while (kseq_read(ks) >= 0) {
		uint8_t *seq, *qual;
		int l_seq;
		kseq2mseq(&r[0], ks);
		if (ks2) {
			if (kseq_read(ks2) < 0) break;
			kseq2mseq(&r[1], ks2);
		} else {
			if (kseq_read(ks) < 0) break;
			kseq2mseq(&r[1], ks);
		}
		l_seq = bwa_pemerge(opt, r, &seq, &qual);
		if (l_seq > 0) {
			++cnt[0];
			if (flag & 1) {
				if (r[0].n.l > 2 && (r[0].n.s[r[0].n.l-1] == '1' || r[0].n.s[r[0].n.l-1] == '2') && r[0].n.s[r[0].n.l-2] == '/')
					r[0].n.s[r[0].n.l-2] = 0, r[0].n.l -= 2;
				print_seq(r[0].n.s, (char*)seq, (char*)qual);
			}
		} else {
			++cnt[-l_seq];
			if (flag & 2) {
				printf("*** %d\n", l_seq);
				print_seq(r[0].n.s, r[0].s.s, r[0].q.l? r[0].q.s : 0);
				print_seq(r[1].n.s, r[1].s.s, r[1].q.l? r[1].q.s : 0);
			}
		}
	}

	fprintf(stderr, "%12ld %s\n", (long)cnt[0], err_msg[0]);
	for (i = 1; i <= MAX_ERR; ++i)
		fprintf(stderr, "%12ld %s\n", (long)cnt[i], err_msg[i]);
	kseq_destroy(ks);
	gzclose(fp);
	if (ks2) {
		kseq_destroy(ks2);
		gzclose(fp2);
	}
	free(opt);
	return 0;
}
