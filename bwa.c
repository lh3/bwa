#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <assert.h>
#include "bntseq.h"
#include "bwa.h"
#include "ksw.h"
#include "utils.h"
#include "kstring.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

int bwa_verbose = 3;
char bwa_rg_id[256];

/************************
 * Batch FASTA/Q reader *
 ************************/

#include "kseq.h"
KSEQ_DECLARE(gzFile)

static inline void trim_readno(kstring_t *s)
{
	if (s->l > 2 && s->s[s->l-2] == '/' && isdigit(s->s[s->l-1]))
		s->l -= 2, s->s[s->l] = 0;
}

static inline void kseq2bseq1(const kseq_t *ks, bseq1_t *s)
{ // TODO: it would be better to allocate one chunk of memory, but probably it does not matter in practice
	s->name = strdup(ks->name.s);
	s->comment = ks->comment.l? strdup(ks->comment.s) : 0;
	s->seq = strdup(ks->seq.s);
	s->qual = ks->qual.l? strdup(ks->qual.s) : 0;
	s->l_seq = strlen(s->seq);
}

bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_)
{
	kseq_t *ks = (kseq_t*)ks1_, *ks2 = (kseq_t*)ks2_;
	int size = 0, m, n;
	bseq1_t *seqs;
	m = n = 0; seqs = 0;
	while (kseq_read(ks) >= 0) {
		if (ks2 && kseq_read(ks2) < 0) { // the 2nd file has fewer reads
			fprintf(stderr, "[W::%s] the 2nd file has fewer sequences.\n", __func__);
			break;
		}
		if (n >= m) {
			m = m? m<<1 : 256;
			seqs = realloc(seqs, m * sizeof(bseq1_t));
		}
		trim_readno(&ks->name);
		kseq2bseq1(ks, &seqs[n]);
		size += seqs[n++].l_seq;
		if (ks2) {
			trim_readno(&ks2->name);
			kseq2bseq1(ks2, &seqs[n]);
			size += seqs[n++].l_seq;
		}
		if (size >= chunk_size && (n&1) == 0) break;
	}
	if (size == 0) { // test if the 2nd file is finished
		if (ks2 && kseq_read(ks2) >= 0)
			fprintf(stderr, "[W::%s] the 1st file has fewer sequences.\n", __func__);
	}
	*n_ = n;
	return seqs;
}

/*****************
 * CIGAR related *
 *****************/

void bwa_fill_scmat(int a, int b, int8_t mat[25])
{
	int i, j, k;
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? a : -b;
		mat[k++] = -1; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = -1;
}

// Generate CIGAR when the alignment end points are known
uint32_t *bwa_gen_cigar(const int8_t mat[25], int q, int r, int w_, int64_t l_pac, const uint8_t *pac, int l_query, uint8_t *query, int64_t rb, int64_t re, int *score, int *n_cigar, int *NM)
{
	uint32_t *cigar = 0;
	uint8_t tmp, *rseq;
	int i;
	int64_t rlen;
	kstring_t str;
	const char *int2base;

	*n_cigar = 0; *NM = -1;
	if (l_query <= 0 || rb >= re || (rb < l_pac && re > l_pac)) return 0; // reject if negative length or bridging the forward and reverse strand
	rseq = bns_get_seq(l_pac, pac, rb, re, &rlen);
	if (re - rb != rlen) goto ret_gen_cigar; // possible if out of range
	if (rb >= l_pac) { // then reverse both query and rseq; this is to ensure indels to be placed at the leftmost position
		for (i = 0; i < l_query>>1; ++i)
			tmp = query[i], query[i] = query[l_query - 1 - i], query[l_query - 1 - i] = tmp;
		for (i = 0; i < rlen>>1; ++i)
			tmp = rseq[i], rseq[i] = rseq[rlen - 1 - i], rseq[rlen - 1 - i] = tmp;
	}
	if (l_query == re - rb && w_ == 0) { // no gap; no need to do DP
		cigar = malloc(4);
		cigar[0] = l_query<<4 | 0;
		*n_cigar = 1;
		for (i = 0, *score = 0; i < l_query; ++i)
			*score += mat[rseq[i]*5 + query[i]];
	} else {
		int w, max_gap, min_w;
		//printf("[Q] "); for (i = 0; i < l_query; ++i) putchar("ACGTN"[(int)query[i]]); putchar('\n');
		//printf("[R] "); for (i = 0; i < re - rb; ++i) putchar("ACGTN"[(int)rseq[i]]); putchar('\n');
		// set the band-width
		max_gap = (int)((double)(((l_query+1)>>1) * mat[0] - q) / r + 1.);
		max_gap = max_gap > 1? max_gap : 1;
		w = (max_gap + abs(rlen - l_query) + 1) >> 1;
		w = w < w_? w : w_;
		min_w = abs(rlen - l_query) + 3;
		w = w > min_w? w : min_w;
		// NW alignment
		*score = ksw_global(l_query, query, rlen, rseq, 5, mat, q, r, w, n_cigar, &cigar);
	}
	{// compute NM
		int k, x, y, u, n_mm = 0, n_gap = 0;
		str.l = str.m = *n_cigar * 4; str.s = (char*)cigar; // append MD to CIGAR
		int2base = rb < l_pac? "ACGTN" : "TGCAN";
		for (k = 0, x = y = u = 0; k < *n_cigar; ++k) {
			int op, len;
			cigar = (uint32_t*)str.s;
			op  = cigar[k]&0xf, len = cigar[k]>>4;
			if (op == 0) { // match
				for (i = 0; i < len; ++i) {
					if (query[x + i] != rseq[y + i]) {
						kputw(u, &str);
						kputc(int2base[rseq[y+i]], &str);
						++n_mm; u = 0;
					} else ++u;
				}
				x += len; y += len;
			} else if (op == 2) { // deletion
				kputw(u, &str); kputc('^', &str);
				for (i = 0; i < len; ++i)
					kputc(int2base[rseq[y+i]], &str);
				u = 0;
				y += len, n_gap += len;
			} else if (op == 1) x += len, n_gap += len; // insertion
		}
		kputw(u, &str); kputc(0, &str);
		*NM = n_mm + n_gap;
		cigar = (uint32_t*)str.s;
	}
	if (rb >= l_pac) // reverse back query
		for (i = 0; i < l_query>>1; ++i)
			tmp = query[i], query[i] = query[l_query - 1 - i], query[l_query - 1 - i] = tmp;

ret_gen_cigar:
	free(rseq);
	return cigar;
}

int bwa_fix_xref(const int8_t mat[25], int q, int r, int w, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, int *qb, int *qe, int64_t *rb, int64_t *re)
{
	int is_rev;
	int64_t cb, ce, fm;
	bntann1_t *ra;
	if (*rb < bns->l_pac && *re > bns->l_pac) { // cross the for-rev boundary; actually with BWA-MEM, we should never come to here
		*qb = *qe = *rb = *re = -1;
		return -1; // unable to fix
	}
	fm = bns_depos(bns, (*rb + *re) >> 1, &is_rev); // coordinate of the middle point on the forward strand
	ra = &bns->anns[bns_pos2rid(bns, fm)]; // annotation of chr corresponding to the middle point
	cb = is_rev? (bns->l_pac<<1) - (ra->offset + ra->len) : ra->offset; // chr start on the mapping strand
	ce = cb + ra->len; // chr end
	if (cb > *rb || ce < *re) { // fix is needed
		int i, score, n_cigar, y, NM;
		uint32_t *cigar;
		int64_t x;
		cb = cb > *rb? cb : *rb;
		ce = ce < *re? ce : *re;
		cigar = bwa_gen_cigar(mat, q, r, w, bns->l_pac, pac, *qe - *qb, query + *qb, *rb, *re, &score, &n_cigar, &NM);
		for (i = 0, x = *rb, y = *qb; i < n_cigar; ++i) {
			int op = cigar[i]&0xf, len = cigar[i]>>4;
			if (op == 0) {
				if (x <= cb && cb < x + len)
					*qb = y + (cb - x), *rb = cb;
				if (x < ce && ce <= x + len) {
					*qe = y + (ce - x), *re = ce;
					break;
				} else x += len, y += len;
			} else if (op == 1) {
				y += len;
			} else if (op == 2) {
				if (x <= cb && cb < x + len)
					*qb = y, *rb = x + len;
				if (x < ce && ce <= x + len) {
					*qe = y, *re = x;
					break;
				} else x += len;
			} else abort(); // should not be here
		}
		free(cigar);
	}
	return (*qb == *qe || *rb == *re)? -2 : 0;
}

/*********************
 * Full index reader *
 *********************/

char *bwa_idx_infer_prefix(const char *hint)
{
	char *prefix;
	int l_hint;
	FILE *fp;
	l_hint = strlen(hint);
	prefix = malloc(l_hint + 3 + 4 + 1);
	strcpy(prefix, hint);
	strcpy(prefix + l_hint, ".64.bwt");
	if ((fp = fopen(prefix, "rb")) != 0) {
		fclose(fp);
		prefix[l_hint + 3] = 0;
		return prefix;
	} else {
		strcpy(prefix + l_hint, ".bwt");
		if ((fp = fopen(prefix, "rb")) == 0) {
			free(prefix);
			return 0;
		} else {
			fclose(fp);
			prefix[l_hint] = 0;
			return prefix;
		}
	}
}

bwt_t *bwa_idx_load_bwt(const char *hint)
{
	char *tmp, *prefix;
	bwt_t *bwt;
	prefix = bwa_idx_infer_prefix(hint);
	if (prefix == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to locate the index files\n", __func__);
		return 0;
	}
	tmp = calloc(strlen(prefix) + 5, 1);
	strcat(strcpy(tmp, prefix), ".bwt"); // FM-index
	bwt = bwt_restore_bwt(tmp);
	strcat(strcpy(tmp, prefix), ".sa");  // partial suffix array (SA)
	bwt_restore_sa(tmp, bwt);
	free(tmp); free(prefix);
	return bwt;
}

bwaidx_t *bwa_idx_load(const char *hint, int which)
{
	bwaidx_t *idx;
	char *prefix;
	prefix = bwa_idx_infer_prefix(hint);
	if (prefix == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] fail to locate the index files\n", __func__);
		return 0;
	}
	idx = calloc(1, sizeof(bwaidx_t));
	if (which & BWA_IDX_BWT) idx->bwt = bwa_idx_load_bwt(hint);
	if (which & BWA_IDX_BNS) {
		idx->bns = bns_restore(prefix);
		if (which & BWA_IDX_PAC) {
			idx->pac = calloc(idx->bns->l_pac/4+1, 1);
			err_fread_noeof(idx->pac, 1, idx->bns->l_pac/4+1, idx->bns->fp_pac); // concatenated 2-bit encoded sequence
			err_fclose(idx->bns->fp_pac);
			idx->bns->fp_pac = 0;
		}
	}
	free(prefix);
	return idx;
}

void bwa_idx_destroy(bwaidx_t *idx)
{
	if (idx == 0) return;
	if (idx->bwt) bwt_destroy(idx->bwt);
	if (idx->bns) bns_destroy(idx->bns);
	if (idx->pac) free(idx->pac);
	free(idx);
}

/***********************
 * SAM header routines *
 ***********************/

void bwa_print_sam_hdr(const bntseq_t *bns, const char *rg_line)
{
	int i;
	extern char *bwa_pg;
	for (i = 0; i < bns->n_seqs; ++i)
		err_printf("@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
	if (rg_line) err_printf("%s\n", rg_line);
	err_printf("%s\n", bwa_pg);
}

static char *bwa_escape(char *s)
{
	char *p, *q;
	for (p = q = s; *p; ++p) {
		if (*p == '\\') {
			++p;
			if (*p == 't') *q++ = '\t';
			else if (*p == 'n') *q++ = '\n';
			else if (*p == 'r') *q++ = '\r';
			else if (*p == '\\') *q++ = '\\';
		} else *q++ = *p;
	}
	*q = '\0';
	return s;
}

char *bwa_set_rg(const char *s)
{
	char *p, *q, *r, *rg_line = 0;
	memset(bwa_rg_id, 0, 256);
	if (strstr(s, "@RG") != s) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] the read group line is not started with @RG\n", __func__);
		goto err_set_rg;
	}
	rg_line = strdup(s);
	bwa_escape(rg_line);
	if ((p = strstr(rg_line, "\tID:")) == 0) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] no ID at the read group line\n", __func__);
		goto err_set_rg;
	}
	p += 4;
	for (q = p; *q && *q != '\t' && *q != '\n'; ++q);
	if (q - p + 1 > 256) {
		if (bwa_verbose >= 1) fprintf(stderr, "[E::%s] @RG:ID is longer than 255 characters\n", __func__);
		goto err_set_rg;
	}
	for (q = p, r = bwa_rg_id; *q && *q != '\t' && *q != '\n'; ++q)
		*r++ = *q;
	return rg_line;

err_set_rg:
	free(rg_line);
	return 0;
}

