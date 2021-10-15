#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <zlib.h>
#include <pthread.h>
#include <errno.h>
#include "ksw.h"
#include "kseq.h"
#include "kstring.h"
#include "bwa.h"
#include "utils.h"
KSEQ_DECLARE(gzFile)

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
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
	"pairs where the end-to-end alignment is inconsistent with SW",
	"pairs potentially with tandem overlaps",
	"pairs with high sum of errors"
};

typedef struct {
	int a, b, q, r, w;
	int q_def, q_thres;
	int T;
	int chunk_size;
	int n_threads;
	int flag; // bit 1: print merged; 2: print unmerged
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
	opt->chunk_size = 10000000;
	opt->n_threads = 1;
	opt->flag = 3;
	bwa_fill_scmat(opt->a, opt->b, opt->mat);
	return opt;
}

int bwa_pemerge(const pem_opt_t *opt, bseq1_t x[2])
{
	uint8_t *s[2], *q[2], *seq, *qual;
	int i, xtra, l, l_seq, sum_q, ret = 0;
	kswr_t r;

	s[0] = malloc(x[0].l_seq); q[0] = malloc(x[0].l_seq);
	s[1] = malloc(x[1].l_seq); q[1] = malloc(x[1].l_seq);
	for (i = 0; i < x[0].l_seq; ++i) {
		int c = x[0].seq[i];
		s[0][i] = c < 0 || c > 127? 4 : c <= 4? c : nst_nt4_table[c];
		q[0][i] = x[0].qual? x[0].qual[i] - 33 : opt->q_def;
	}
	for (i = 0; i < x[1].l_seq; ++i) {
		int c = x[1].seq[x[1].l_seq - 1 - i];
		c = c < 0 || c > 127? 4 : c < 4? c : nst_nt4_table[c];
		s[1][i] = c < 4? 3 - c : 4;
		q[1][i] = x[1].qual? x[1].qual[x[1].l_seq - 1 - i] - 33 : opt->q_def;
	}

	xtra = KSW_XSTART | KSW_XSUBO;
	r = ksw_align(x[1].l_seq, s[1], x[0].l_seq, s[0], 5, opt->mat, opt->q, opt->r, xtra, 0);
	++r.qe; ++r.te; // change to the half-close-half-open coordinates

	if (r.score < opt->T) { ret = -1; goto pem_ret; } // poor alignment
	if (r.tb < r.qb) { ret = -2; goto pem_ret; } // no enough space for the left end
	if (x[0].l_seq - r.te > x[1].l_seq - r.qe) { ret = -3; goto pem_ret; } // no enough space for the right end
	if ((double)r.score2 / r.score >= MAX_SCORE_RATIO) { ret = -4; goto pem_ret; } // the second best score is too large
	if (r.qe - r.qb != r.te - r.tb) { ret = -5; goto pem_ret; } // we do not allow gaps

	{ // test tandem match; O(n^2)
		int max_m, max_m2, min_l, max_l, max_l2;
		max_m = max_m2 = 0; max_l = max_l2 = 0;
		min_l = x[0].l_seq < x[1].l_seq? x[0].l_seq : x[1].l_seq;
		for (l = 1; l < min_l; ++l) {
			int m = 0, o = x[0].l_seq - l;
			uint8_t *s0o = &s[0][o], *s1 = s[1];
			for (i = 0; i < l; ++i) // TODO: in principle, this can be done with SSE2. It is the bottleneck!
				m += opt->mat[(s1[i]<<2) + s1[i] + s0o[i]]; // equivalent to s[1][i]*5 + s[0][o+i]
			if (m > max_m) max_m2 = max_m, max_m = m, max_l2 = max_l, max_l = l;
			else if (m > max_m2) max_m2 = m, max_l2 = l;
		}
		if (max_m < opt->T || max_l != x[0].l_seq - (r.tb - r.qb)) { ret = -6; goto pem_ret; }
		if (max_l2 < max_l && max_m2 >= opt->T && (double)(max_m2 + (max_l - max_l2) * opt->a) / max_m >= MAX_SCORE_RATIO) {
			ret = -7; goto pem_ret;
		}
		if (max_l2 > max_l && (double)max_m2 / max_m >= MAX_SCORE_RATIO) { ret = -7; goto pem_ret; }
	}

	l = x[0].l_seq - (r.tb - r.qb); // length to merge
	l_seq = x[0].l_seq + x[1].l_seq - l;
	seq = malloc(l_seq + 1);
	qual = malloc(l_seq + 1);
	memcpy(seq,  s[0], x[0].l_seq); memcpy(seq  + x[0].l_seq, &s[1][l], x[1].l_seq - l);
	memcpy(qual, q[0], x[0].l_seq); memcpy(qual + x[0].l_seq, &q[1][l], x[1].l_seq - l);
	for (i = 0, sum_q = 0; i < l; ++i) {
		int k = x[0].l_seq - l + i;
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
		ret = -8; goto pem_ret;
	}

	for (i = 0; i < l_seq; ++i) seq[i] = "ACGTN"[(int)seq[i]], qual[i] += 33;
	seq[l_seq] = qual[l_seq] = 0;

	free(x[1].name); free(x[1].seq); free(x[1].qual); free(x[1].comment);
	memset(&x[1], 0, sizeof(bseq1_t));
	free(x[0].seq); free(x[0].qual);
	x[0].l_seq = l_seq; x[0].seq = (char*)seq; x[0].qual = (char*)qual;

pem_ret:
	free(s[0]); free(s[1]); free(q[0]); free(q[1]);
	return ret;
}

static inline void print_bseq(const bseq1_t *s, int rn)
{
	err_putchar(s->qual? '@' : '>');
	err_fputs(s->name, stdout);
	if (rn == 1 || rn == 2) {
		err_putchar('/'); err_putchar('0' + rn); err_putchar('\n');
	} else err_puts(" merged");
	err_puts(s->seq);
	if (s->qual) {
		err_puts("+"); err_puts(s->qual);
	}
}

typedef struct {
	int n, start;
	bseq1_t *seqs;
	int64_t cnt[MAX_ERR+1];
	const pem_opt_t *opt;
} worker_t;

void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	int i;
	for (i = w->start; i < w->n>>1; i += w->opt->n_threads)
		++w->cnt[-bwa_pemerge(w->opt, &w->seqs[i<<1])];
	return 0;
}

static void process_seqs(const pem_opt_t *opt, int n_, bseq1_t *seqs, int64_t cnt[MAX_ERR+1])
{
	int i, j, n = n_>>1<<1;
	worker_t *w;

	w = calloc(opt->n_threads, sizeof(worker_t));
	for (i = 0; i < opt->n_threads; ++i) {
		worker_t *p = &w[i];
		p->start = i; p->n = n;
		p->opt = opt;
		p->seqs = seqs;
	}
	if (opt->n_threads == 1) {
		worker(w);
	} else {
		pthread_t *tid;
		tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
		for (i = 0; i < opt->n_threads; ++i) pthread_create(&tid[i], 0, worker, &w[i]);
		for (i = 0; i < opt->n_threads; ++i) pthread_join(tid[i], 0);
		free(tid);
	}
	for (i = 0; i < opt->n_threads; ++i) {
		worker_t *p = &w[i];
		for (j = 0; j <= MAX_ERR; ++j) cnt[j] += p->cnt[j];
	}
	free(w);
	for (i = 0; i < n>>1; ++i) {
		if (seqs[i<<1|1].l_seq != 0) {
			if (opt->flag&2) {
				print_bseq(&seqs[i<<1|0], 1);
				print_bseq(&seqs[i<<1|1], 2);
			}
		} else if (opt->flag&1)
			print_bseq(&seqs[i<<1|0], 0);
	}
	for (i = 0; i < n; ++i) {
		bseq1_t *s = &seqs[i];
		free(s->name); free(s->seq); free(s->qual); free(s->comment);
	}
}

int main_pemerge(int argc, char *argv[])
{
	int c, flag = 0, i, n, min_ovlp = 10;
	int64_t cnt[MAX_ERR+1];
	bseq1_t *bseq;
	gzFile fp, fp2 = 0;
	kseq_t *ks, *ks2 = 0;
	pem_opt_t *opt;

	opt = pem_opt_init();
	while ((c = getopt(argc, argv, "muQ:t:T:")) >= 0) {
		if (c == 'm') flag |= 1;
		else if (c == 'u') flag |= 2;
		else if (c == 'Q') opt->q_thres = atoi(optarg);
		else if (c == 't') opt->n_threads = atoi(optarg);
		else if (c == 'T') min_ovlp = atoi(optarg);
		else return 1;
	}
	if (flag == 0) flag = 3;
	opt->flag = flag;
	opt->T = opt->a * min_ovlp;

	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bwa pemerge [-mu] <read1.fq> [read2.fq]\n\n");
		fprintf(stderr, "Options: -m       output merged reads only\n");
		fprintf(stderr, "         -u       output unmerged reads only\n");
		fprintf(stderr, "         -t INT   number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, "         -T INT   minimum end overlap [%d]\n", min_ovlp);
		fprintf(stderr, "         -Q INT   max sum of errors [%d]\n", opt->q_thres);
		fprintf(stderr, "\n");
		free(opt);
		return 1;
	}

	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	if (NULL == fp) {
		fprintf(stderr, "Couldn't open %s : %s\n",
				strcmp(argv[optind], "-") ? argv[optind] : "stdin",
				errno ? strerror(errno) : "Out of memory");
		exit(EXIT_FAILURE);
	}
	ks = kseq_init(fp);
	if (optind + 1 < argc) {
		fp2 = strcmp(argv[optind+1], "-")? gzopen(argv[optind+1], "r") : gzdopen(fileno(stdin), "r");
		if (NULL == fp) {
			fprintf(stderr, "Couldn't open %s : %s\n",
					strcmp(argv[optind+1], "-") ? argv[optind+1] : "stdin",
					errno ? strerror(errno) : "Out of memory");
			exit(EXIT_FAILURE);
		}
		ks2 = kseq_init(fp2);
	}

	memset(cnt, 0, 8 * (MAX_ERR+1));
	while ((bseq = bseq_read(opt->n_threads * opt->chunk_size, &n, ks, ks2)) != 0) {
		process_seqs(opt, n, bseq, cnt);
		free(bseq);
	}

	fprintf(stderr, "%12ld %s\n", (long)cnt[0], err_msg[0]);
	for (i = 1; i <= MAX_ERR; ++i)
		fprintf(stderr, "%12ld %s\n", (long)cnt[i], err_msg[i]);
	kseq_destroy(ks);
	err_gzclose(fp);
	if (ks2) {
		kseq_destroy(ks2);
		err_gzclose(fp2);
	}
	free(opt);

	err_fflush(stdout);

	return 0;
}
