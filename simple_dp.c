#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>
#include "stdaln.h"
#include "utils.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	int l;
	unsigned char *s;
	char *n;
} seq1_t;

typedef struct {
	int n_seqs, m_seqs;
	seq1_t *seqs;
} seqs_t;

unsigned char aln_rev_table[256] = {
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','T','V','G', 'H','N','N','C', 'D','N','N','M', 'N','K','N','N',
	'N','N','Y','S', 'A','N','B','W', 'X','R','N','N', 'N','N','N','N',
	'N','t','v','g', 'h','N','N','c', 'd','N','N','m', 'N','k','N','N',
	'N','N','y','s', 'a','N','b','w', 'x','r','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N'
};

static int g_is_global = 0, g_thres = 1, g_strand = 0, g_aa = 0;
static AlnParam g_aln_param;

static void revseq(int len, uint8_t *seq)
{
	int i;
	for (i = 0; i < len>>1; ++i) {
		uint8_t tmp = aln_rev_table[seq[len-1-i]];
		seq[len-1-i] = aln_rev_table[seq[i]];
		seq[i] = tmp;
	}
	if (len&1) seq[i] = aln_rev_table[seq[i]];
}

static seqs_t *load_seqs(const char *fn)
{
	seqs_t *s;
	seq1_t *p;
	gzFile fp;
	int l;
	kseq_t *seq;

	fp = xzopen(fn, "r");
	seq = kseq_init(fp);
	s = (seqs_t*)calloc(1, sizeof(seqs_t));
	s->m_seqs = 256;
	s->seqs = (seq1_t*)calloc(s->m_seqs, sizeof(seq1_t));
	while ((l = kseq_read(seq)) >= 0) {
		if (s->n_seqs == s->m_seqs) {
			s->m_seqs <<= 1;
			s->seqs = (seq1_t*)realloc(s->seqs, s->m_seqs * sizeof(seq1_t));
		}
		p = s->seqs + (s->n_seqs++);
		p->l = seq->seq.l;
		p->s = (unsigned char*)malloc(p->l + 1);
		memcpy(p->s, seq->seq.s, p->l);
		p->s[p->l] = 0;
		p->n = strdup((const char*)seq->name.s);
	}
	kseq_destroy(seq);
	gzclose(fp);
	fprintf(stderr, "[load_seqs] %d sequences are loaded.\n", s->n_seqs);
	return s;
}

static void aln_1seq(const seqs_t *ss, const char *name, int l, const char *s, char strand)
{
	int i;
	for (i = 0; i < ss->n_seqs; ++i) {
		AlnAln *aa;
		seq1_t *p = ss->seqs + i;
		g_aln_param.band_width = l + p->l;
		aa = aln_stdaln_aux(s, (const char*)p->s, &g_aln_param, g_is_global, g_thres, l, p->l);
		if (aa->score >= g_thres || g_is_global) {
			printf(">%s\t%d\t%d\t%s\t%c\t%d\t%d\t%d\t%d\t", p->n, aa->start1? aa->start1 : 1, aa->end1, name, strand,
				   aa->start2? aa->start2 : 1, aa->end2, aa->score, aa->subo);
			// NB: I put the short sequence as the first sequence in SW, an insertion to
			// the reference becomes a deletion from the short sequence. Therefore, I use
			// "MDI" here rather than "MID", and print ->out2 first rather than ->out1.
			for (i = 0; i != aa->n_cigar; ++i)
				printf("%d%c", aa->cigar32[i]>>4, "MDI"[aa->cigar32[i]&0xf]);
			printf("\n%s\n%s\n%s\n", aa->out2, aa->outm, aa->out1);
		}
		aln_free_AlnAln(aa);
	}
}

static void aln_seqs(const seqs_t *ss, const char *fn)
{
	gzFile fp;
	kseq_t *seq;
	int l;

	fp = xzopen(fn, "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		if (g_strand&1) aln_1seq(ss, (char*)seq->name.s, l, seq->seq.s, '+');
		if (g_strand&2) {
			revseq(l, (uint8_t*)seq->seq.s);
			aln_1seq(ss, (char*)seq->name.s, l, seq->seq.s, '-');
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
}

int bwa_stdsw(int argc, char *argv[])
{
	int c;
	seqs_t *ss;

	while ((c = getopt(argc, argv, "gT:frp")) >= 0) {
		switch (c) {
		case 'g': g_is_global = 1; break;
		case 'T': g_thres = atoi(optarg); break;
		case 'f': g_strand |= 1; break;
		case 'r': g_strand |= 2; break;
		case 'p': g_aa = 1; break;
		}
	}
	if (g_strand == 0) g_strand = 3;
	if (g_aa) g_strand = 1;
	if (optind + 1 >= argc) {
		fprintf(stderr, "\nUsage:   bwa stdsw [options] <seq1.long.fa> <seq2.short.fa>\n\n");
		fprintf(stderr, "Options: -T INT    minimum score [%d]\n", g_thres);
		fprintf(stderr, "         -p        protein alignment (suppressing -r)\n");
		fprintf(stderr, "         -f        forward strand only\n");
		fprintf(stderr, "         -r        reverse strand only\n");
		fprintf(stderr, "         -g        global alignment\n\n");
		fprintf(stderr, "Note: This program is specifically designed for alignment between multiple short\n");
		fprintf(stderr, "      sequences and ONE long sequence. It outputs the suboptimal score on the long\n");
		fprintf(stderr, "      sequence.\n\n");
		return 1;
	}
	g_aln_param = g_aa? aln_param_aa2aa : aln_param_blast;
	g_aln_param.gap_end = 0;
	ss = load_seqs(argv[optind]);
	aln_seqs(ss, argv[optind+1]);
	return 0;
}
