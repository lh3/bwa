#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "bwase.h"
#include "bwtaln.h"
#include "bntseq.h"
#include "utils.h"
#include "kstring.h"
#include "bwa.h"
#include "ksw.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

int g_log_n[256];

void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi)
{
	int i, cnt, best;
	if (n_aln == 0) {
		s->type = BWA_TYPE_NO_MATCH;
		s->c1 = s->c2 = 0;
		return;
	}

	if (set_main) {
		best = aln[0].score;
		for (i = cnt = 0; i < n_aln; ++i) {
			const bwt_aln1_t *p = aln + i;
			if (p->score > best) break;
			if (drand48() * (p->l - p->k + 1 + cnt) > (double)cnt) {
				s->n_mm = p->n_mm; s->n_gapo = p->n_gapo; s->n_gape = p->n_gape;
				s->ref_shift = (int)p->n_del - (int)p->n_ins;
				s->score = p->score;
				s->sa = p->k + (bwtint_t)((p->l - p->k + 1) * drand48());
			}
			cnt += p->l - p->k + 1;
		}
		s->c1 = cnt;
		for (; i < n_aln; ++i) cnt += aln[i].l - aln[i].k + 1;
		s->c2 = cnt - s->c1;
		s->type = s->c1 > 1? BWA_TYPE_REPEAT : BWA_TYPE_UNIQUE;
	}

	if (n_multi) {
		int k, rest, n_occ, z = 0;
		for (k = n_occ = 0; k < n_aln; ++k) {
			const bwt_aln1_t *q = aln + k;
			n_occ += q->l - q->k + 1;
		}
		if (s->multi) free(s->multi);
		if (n_occ > n_multi + 1) { // if there are too many hits, generate none of them
			s->multi = 0; s->n_multi = 0;
			return;
		}
		/* The following code is more flexible than what is required
		 * here. In principle, due to the requirement above, we can
		 * simply output all hits, but the following samples "rest"
		 * number of random hits. */
		rest = n_occ > n_multi + 1? n_multi + 1 : n_occ; // find one additional for ->sa
		s->multi = calloc(rest, sizeof(bwt_multi1_t));
		for (k = 0; k < n_aln; ++k) {
			const bwt_aln1_t *q = aln + k;
			if (q->l - q->k + 1 <= rest) {
				bwtint_t l;
				for (l = q->k; l <= q->l; ++l) {
					s->multi[z].pos = l;
					s->multi[z].gap = q->n_gapo + q->n_gape;
					s->multi[z].ref_shift = (int)q->n_del - (int)q->n_ins;
					s->multi[z++].mm = q->n_mm;
				}
				rest -= q->l - q->k + 1;
			} else { // Random sampling (http://code.activestate.com/recipes/272884/). In fact, we never come here. 
				int j, i;
				for (j = rest, i = q->l - q->k + 1; j > 0; --j) {
					double p = 1.0, x = drand48();
					while (x < p) p -= p * j / (i--);
					s->multi[z].pos = q->l - i;
					s->multi[z].gap = q->n_gapo + q->n_gape;
					s->multi[z].ref_shift = (int)q->n_del - (int)q->n_ins;
					s->multi[z++].mm = q->n_mm;
				}
				rest = 0;
				break;
			}
		}
		s->n_multi = z;
	}
}

void bwa_aln2seq(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s)
{
	bwa_aln2seq_core(n_aln, aln, s, 1, 0);
}

int bwa_approx_mapQ(const bwa_seq_t *p, int mm)
{
	int n;
	if (p->c1 == 0) return 23;
	if (p->c1 > 1) return 0;
	if (p->n_mm == mm) return 25;
	if (p->c2 == 0) return 37;
	n = (p->c2 >= 255)? 255 : p->c2;
	return (23 < g_log_n[n])? 0 : 23 - g_log_n[n];
}

bwtint_t bwa_sa2pos(const bntseq_t *bns, const bwt_t *bwt, bwtint_t sapos, int ref_len, int *strand)
{
	bwtint_t pos_f;
	int is_rev;
	*strand = 0; // initialise strand to 0 otherwise we could return without setting it
	pos_f = bwt_sa(bwt, sapos); // position on the forward-reverse coordinate
	if (pos_f < bns->l_pac && bns->l_pac < pos_f + ref_len) return (bwtint_t)-1;
	pos_f = bns_depos(bns, pos_f, &is_rev); // position on the forward strand; this may be the first base or the last base
	*strand = !is_rev;
	if (is_rev) pos_f = pos_f + 1 < ref_len? 0 : pos_f - ref_len + 1; // position of the first base
	return pos_f; // FIXME: it is possible that pos_f < bns->anns[ref_id].offset
}

/**
 * Derive the actual position in the read from the given suffix array
 * coordinates. Note that the position will be approximate based on
 * whether indels appear in the read and whether calculations are
 * performed from the start or end of the read.
 */
void bwa_cal_pac_pos_core(const bntseq_t *bns, const bwt_t *bwt, bwa_seq_t *seq, const int max_mm, const float fnr)
{
	int max_diff, strand;
	if (seq->type != BWA_TYPE_UNIQUE && seq->type != BWA_TYPE_REPEAT) return;
	max_diff = fnr > 0.0? bwa_cal_maxdiff(seq->len, BWA_AVG_ERR, fnr) : max_mm;
	seq->seQ = seq->mapQ = bwa_approx_mapQ(seq, max_diff);
	//fprintf(stderr, "%d\n", seq->ref_shift);
	seq->pos = bwa_sa2pos(bns, bwt, seq->sa, seq->len + seq->ref_shift, &strand);
	seq->strand = strand;
	seq->seQ = seq->mapQ = bwa_approx_mapQ(seq, max_diff);
	if (seq->pos == (bwtint_t)-1) seq->type = BWA_TYPE_NO_MATCH;
}

void bwa_cal_pac_pos(const bntseq_t *bns, const char *prefix, int n_seqs, bwa_seq_t *seqs, int max_mm, float fnr)
{
	int i, j, strand, n_multi;
	char str[1024];
	bwt_t *bwt;
	// load forward SA
	strcpy(str, prefix); strcat(str, ".bwt");  bwt = bwt_restore_bwt(str);
	strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt);
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p = &seqs[i];
		bwa_cal_pac_pos_core(bns, bwt, p, max_mm, fnr);
		for (j = n_multi = 0; j < p->n_multi; ++j) {
			bwt_multi1_t *q = p->multi + j;
			q->pos = bwa_sa2pos(bns, bwt, q->pos, p->len + q->ref_shift, &strand);
			q->strand = strand;
			if (q->pos != p->pos && q->pos != (bwtint_t)-1)
				p->multi[n_multi++] = *q;
		}
		p->n_multi = n_multi;
	}
	bwt_destroy(bwt);
}

#define SW_BW 50

bwa_cigar_t *bwa_refine_gapped_core(bwtint_t l_pac, const ubyte_t *pacseq, int len, ubyte_t *seq, int ref_shift, bwtint_t *_rb, int *n_cigar)
{
	bwa_cigar_t *cigar = 0;
	uint32_t *cigar32 = 0;
	ubyte_t *rseq;
	int64_t k, rb, re, rlen;
	int8_t mat[25];
	int w;

	bwa_fill_scmat(1, 3, mat);
	rb = *_rb; re = rb + len + ref_shift;
	assert(re <= l_pac);
	rseq = bns_get_seq(l_pac, pacseq, rb, re, &rlen);
	assert(re - rb == rlen);
	w = abs((int)rlen - len) * 1.5;
	ksw_global(len, seq, rlen, rseq, 5, mat, 5, 1, SW_BW > w? SW_BW : w, n_cigar, &cigar32);
	assert(*n_cigar > 0);
	if ((cigar32[*n_cigar - 1]&0xf) == 1) cigar32[*n_cigar - 1] = (cigar32[*n_cigar - 1]>>4<<4) | 3; // change endding ins to soft clipping
	if ((cigar32[0]&0xf) == 1) cigar32[0] = (cigar32[0]>>4<<4) | 3; // change beginning ins to soft clipping
	if ((cigar32[*n_cigar - 1]&0xf) == 2) --*n_cigar; // delete endding del
	if ((cigar32[0]&0xf) == 2) { // delete beginning del
		*_rb += cigar32[0]>>4;
		--*n_cigar;
		memmove(cigar32, cigar32+1, (*n_cigar) * 4);
	}
	cigar = (bwa_cigar_t*)cigar32;
	for (k = 0; k < *n_cigar; ++k)
		cigar[k] = __cigar_create((cigar32[k]&0xf), (cigar32[k]>>4));
	free(rseq);
	return cigar;
}

char *bwa_cal_md1(int n_cigar, bwa_cigar_t *cigar, int len, bwtint_t pos, ubyte_t *seq,
				  bwtint_t l_pac, ubyte_t *pacseq, kstring_t *str, int *_nm)
{
	bwtint_t x, y;
	int z, u, c, nm = 0;
	str->l = 0; // reset
	x = pos; y = 0;
	if (cigar) {
		int k, l;
		for (k = u = 0; k < n_cigar; ++k) {
			l = __cigar_len(cigar[k]);
			if (__cigar_op(cigar[k]) == FROM_M) {
				for (z = 0; z < l && x+z < l_pac; ++z) {
					c = pacseq[(x+z)>>2] >> ((~(x+z)&3)<<1) & 3;
					if (c > 3 || seq[y+z] > 3 || c != seq[y+z]) {
						ksprintf(str, "%d", u);
						kputc("ACGTN"[c], str);
						++nm;
						u = 0;
					} else ++u;
				}
				x += l; y += l;
			} else if (__cigar_op(cigar[k]) == FROM_I || __cigar_op(cigar[k]) == FROM_S) {
				y += l;
				if (__cigar_op(cigar[k]) == FROM_I) nm += l;
			} else if (__cigar_op(cigar[k]) == FROM_D) {
				ksprintf(str, "%d", u);
				kputc('^', str);
				for (z = 0; z < l && x+z < l_pac; ++z)
					kputc("ACGT"[pacseq[(x+z)>>2] >> ((~(x+z)&3)<<1) & 3], str);
				u = 0;
				x += l; nm += l;
			}
		}
	} else { // no gaps
		for (z = u = 0; z < (bwtint_t)len && x+z < l_pac; ++z) {
			c = pacseq[(x+z)>>2] >> ((~(x+z)&3)<<1) & 3;
			if (c > 3 || seq[y+z] > 3 || c != seq[y+z]) {
				ksprintf(str, "%d", u);
				kputc("ACGTN"[c], str);
				++nm;
				u = 0;
			} else ++u;
		}
	}
	ksprintf(str, "%d", u);
	*_nm = nm;
	return strdup(str->s);
}

void bwa_correct_trimmed(bwa_seq_t *s)
{
	if (s->len == s->full_len) return;
	if (s->strand == 0) { // forward
		if (s->cigar && __cigar_op(s->cigar[s->n_cigar-1]) == FROM_S) { // the last is S
			s->cigar[s->n_cigar-1] += s->full_len - s->len;
		} else {
			if (s->cigar == 0) {
				s->n_cigar = 2;
				s->cigar = calloc(s->n_cigar, sizeof(bwa_cigar_t));
				s->cigar[0] = __cigar_create(0, s->len);
			} else {
				++s->n_cigar;
				s->cigar = realloc(s->cigar, s->n_cigar * sizeof(bwa_cigar_t));
			}
			s->cigar[s->n_cigar-1] = __cigar_create(3, (s->full_len - s->len));
		}
	} else { // reverse
		if (s->cigar && __cigar_op(s->cigar[0]) == FROM_S) { // the first is S
			s->cigar[0] += s->full_len - s->len;
		} else {
			if (s->cigar == 0) {
				s->n_cigar = 2;
				s->cigar = calloc(s->n_cigar, sizeof(bwa_cigar_t));
				s->cigar[1] = __cigar_create(0, s->len);
			} else {
				++s->n_cigar;
				s->cigar = realloc(s->cigar, s->n_cigar * sizeof(bwa_cigar_t));
				memmove(s->cigar + 1, s->cigar, (s->n_cigar-1) * sizeof(bwa_cigar_t));
			}
			s->cigar[0] = __cigar_create(3, (s->full_len - s->len));
		}
	}
	s->len = s->full_len;
}

void bwa_refine_gapped(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, ubyte_t *_pacseq)
{
	ubyte_t *pacseq;
	int i, j, k;
	kstring_t *str;

	if (!_pacseq) {
		pacseq = (ubyte_t*)calloc(bns->l_pac/4+1, 1);
		err_rewind(bns->fp_pac);
		err_fread_noeof(pacseq, 1, bns->l_pac/4+1, bns->fp_pac);
	} else pacseq = _pacseq;
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *s = seqs + i;
		seq_reverse(s->len, s->seq, 0); // IMPORTANT: s->seq is reversed here!!!
		for (j = k = 0; j < s->n_multi; ++j) {
			bwt_multi1_t *q = s->multi + j;
			int n_cigar;
			if (q->gap) { // gapped alignment
				q->cigar = bwa_refine_gapped_core(bns->l_pac, pacseq, s->len, q->strand? s->rseq : s->seq, q->ref_shift, &q->pos, &n_cigar);
				q->n_cigar = n_cigar;
				if (q->cigar) s->multi[k++] = *q;
			} else s->multi[k++] = *q;
		}
		s->n_multi = k; // this squeezes out gapped alignments which failed the CIGAR generation
		if (s->type == BWA_TYPE_NO_MATCH || s->type == BWA_TYPE_MATESW || s->n_gapo == 0) continue;
		s->cigar = bwa_refine_gapped_core(bns->l_pac, pacseq, s->len, s->strand? s->rseq : s->seq, s->ref_shift, &s->pos, &s->n_cigar);
		if (s->cigar == 0) s->type = BWA_TYPE_NO_MATCH;
	}
	// generate MD tag
	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *s = seqs + i;
		if (s->type != BWA_TYPE_NO_MATCH) {
			int nm;
			s->md = bwa_cal_md1(s->n_cigar, s->cigar, s->len, s->pos, s->strand? s->rseq : s->seq, bns->l_pac, pacseq, str, &nm);
			s->nm = nm;
		}
	}
	free(str->s); free(str);

	// correct for trimmed reads
	for (i = 0; i < n_seqs; ++i) bwa_correct_trimmed(seqs + i);

	if (!_pacseq) free(pacseq);
}

int64_t pos_end(const bwa_seq_t *p)
{
	if (p->cigar) {
		int j;
		int64_t x = p->pos;
		for (j = 0; j != p->n_cigar; ++j) {
			int op = __cigar_op(p->cigar[j]);
			if (op == 0 || op == 2) x += __cigar_len(p->cigar[j]);
		}
		return x;
	} else return p->pos + p->len;
}

int64_t pos_end_multi(const bwt_multi1_t *p, int len) // analogy to pos_end()
{
	if (p->cigar) {
		int j;
		int64_t x = p->pos;
		for (j = 0; j != p->n_cigar; ++j) {
			int op = __cigar_op(p->cigar[j]);
			if (op == 0 || op == 2) x += __cigar_len(p->cigar[j]);
		}
		return x;
	} else return p->pos + len;
}

static int64_t pos_5(const bwa_seq_t *p)
{
	if (p->type != BWA_TYPE_NO_MATCH)
		return p->strand? pos_end(p) : p->pos;
	return -1;
}

void bwa_print_seq(FILE *stream, bwa_seq_t *seq) {
	char buffer[4096];
	const int bsz = sizeof(buffer);
	int i, j, l;
	
	if (seq->strand == 0) {
		for (i = 0; i < seq->full_len; i += bsz) {
			l = seq->full_len - i > bsz ? bsz : seq->full_len - i;
			for (j = 0; j < l; j++) buffer[j] = "ACGTN"[seq->seq[i + j]];
			err_fwrite(buffer, 1, l, stream);
		}
	} else {
		for (i = seq->full_len - 1; i >= 0; i -= bsz) {
			l = i + 1 > bsz ? bsz : i + 1;
			for (j = 0; j < l; j++) buffer[j] = "TGCAN"[seq->seq[i - j]];
			err_fwrite(buffer, 1, l, stream);
		}
	}
}

void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2)
{
	int j;
	if (p->type != BWA_TYPE_NO_MATCH || (mate && mate->type != BWA_TYPE_NO_MATCH)) {
		int seqid, nn, am = 0, flag = p->extra_flag;
		char XT;

		if (p->type == BWA_TYPE_NO_MATCH) {
			p->pos = mate->pos;
			p->strand = mate->strand;
			flag |= SAM_FSU;
			j = 1;
		} else j = pos_end(p) - p->pos; // j is the length of the reference in the alignment

		// get seqid
		nn = bns_cnt_ambi(bns, p->pos, j, &seqid);
		if (p->type != BWA_TYPE_NO_MATCH && p->pos + j - bns->anns[seqid].offset > bns->anns[seqid].len)
			flag |= SAM_FSU; // flag UNMAP as this alignment bridges two adjacent reference sequences

		// update flag and print it
		if (p->strand) flag |= SAM_FSR;
		if (mate) {
			if (mate->type != BWA_TYPE_NO_MATCH) {
				if (mate->strand) flag |= SAM_FMR;
			} else flag |= SAM_FMU;
		}
		err_printf("%s\t%d\t%s\t", p->name, flag, bns->anns[seqid].name);
		err_printf("%d\t%d\t", (int)(p->pos - bns->anns[seqid].offset + 1), p->mapQ);

		// print CIGAR
		if (p->cigar) {
			for (j = 0; j != p->n_cigar; ++j)
				err_printf("%d%c", __cigar_len(p->cigar[j]), "MIDS"[__cigar_op(p->cigar[j])]);
		} else if (p->type == BWA_TYPE_NO_MATCH) err_printf("*");
		else err_printf("%dM", p->len);

		// print mate coordinate
		if (mate && mate->type != BWA_TYPE_NO_MATCH) {
			int m_seqid;
			long long isize;
			am = mate->seQ < p->seQ? mate->seQ : p->seQ; // smaller single-end mapping quality
			// redundant calculation here, but should not matter too much
			bns_cnt_ambi(bns, mate->pos, mate->len, &m_seqid);
			err_printf("\t%s\t", (seqid == m_seqid)? "=" : bns->anns[m_seqid].name);
			isize = (seqid == m_seqid)? pos_5(mate) - pos_5(p) : 0;
			if (p->type == BWA_TYPE_NO_MATCH) isize = 0;
			err_printf("%d\t%lld\t", (int)(mate->pos - bns->anns[m_seqid].offset + 1), isize);
		} else if (mate) err_printf("\t=\t%d\t0\t", (int)(p->pos - bns->anns[seqid].offset + 1));
		else err_printf("\t*\t0\t0\t");

		// print sequence and quality
		bwa_print_seq(stdout, p);
		err_putchar('\t');
		if (p->qual) {
			if (p->strand) seq_reverse(p->len, p->qual, 0); // reverse quality
			err_printf("%s", p->qual);
		} else err_printf("*");

		if (bwa_rg_id[0]) err_printf("\tRG:Z:%s", bwa_rg_id);
		if (p->bc[0]) err_printf("\tBC:Z:%s", p->bc);
		if (p->clip_len < p->full_len) err_printf("\tXC:i:%d", p->clip_len);
		if (p->type != BWA_TYPE_NO_MATCH) {
			int i;
			// calculate XT tag
			XT = "NURM"[p->type];
			if (nn > 10) XT = 'N';
			// print tags
			err_printf("\tXT:A:%c\t%s:i:%d", XT, (mode & BWA_MODE_COMPREAD)? "NM" : "CM", p->nm);
			if (nn) err_printf("\tXN:i:%d", nn);
			if (mate) err_printf("\tSM:i:%d\tAM:i:%d", p->seQ, am);
			if (p->type != BWA_TYPE_MATESW) { // X0 and X1 are not available for this type of alignment
				err_printf("\tX0:i:%d", p->c1);
				if (p->c1 <= max_top2) err_printf("\tX1:i:%d", p->c2);
			}
			err_printf("\tXM:i:%d\tXO:i:%d\tXG:i:%d", p->n_mm, p->n_gapo, p->n_gapo+p->n_gape);
			if (p->md) err_printf("\tMD:Z:%s", p->md);
			// print multiple hits
			if (p->n_multi) {
				err_printf("\tXA:Z:");
				for (i = 0; i < p->n_multi; ++i) {
					bwt_multi1_t *q = p->multi + i;
					int k;
					j = pos_end_multi(q, p->len) - q->pos;
					nn = bns_cnt_ambi(bns, q->pos, j, &seqid);
					err_printf("%s,%c%d,", bns->anns[seqid].name, q->strand? '-' : '+',
						   (int)(q->pos - bns->anns[seqid].offset + 1));
					if (q->cigar) {
						for (k = 0; k < q->n_cigar; ++k)
							err_printf("%d%c", __cigar_len(q->cigar[k]), "MIDS"[__cigar_op(q->cigar[k])]);
					} else err_printf("%dM", p->len);
					err_printf(",%d;", q->gap + q->mm);
				}
			}
		}
		err_putchar('\n');
	} else { // this read has no match
		//ubyte_t *s = p->strand? p->rseq : p->seq;
		int flag = p->extra_flag | SAM_FSU;
		if (mate && mate->type == BWA_TYPE_NO_MATCH) flag |= SAM_FMU;
		err_printf("%s\t%d\t*\t0\t0\t*\t*\t0\t0\t", p->name, flag);
		//Why did this work differently to the version above??
		//for (j = 0; j != p->len; ++j) putchar("ACGTN"[(int)s[j]]);
		bwa_print_seq(stdout, p);
		err_putchar('\t');
		if (p->qual) {
			if (p->strand) seq_reverse(p->len, p->qual, 0); // reverse quality
			err_printf("%s", p->qual);
		} else err_printf("*");
		if (bwa_rg_id[0]) err_printf("\tRG:Z:%s", bwa_rg_id);
		if (p->bc[0]) err_printf("\tBC:Z:%s", p->bc);
		if (p->clip_len < p->full_len) err_printf("\tXC:i:%d", p->clip_len);
		err_putchar('\n');
	}
}

void bwase_initialize() 
{
	int i;
	for (i = 1; i != 256; ++i) g_log_n[i] = (int)(4.343 * log(i) + 0.5);
}

void bwa_sai2sam_se_core(const char *prefix, const char *fn_sa, const char *fn_fa, int n_occ, const char *rg_line)
{
	extern bwa_seqio_t *bwa_open_reads(int mode, const char *fn_fa);
	int i, n_seqs, m_aln;
	long long tot_seqs = 0;
	bwt_aln1_t *aln = 0;
	bwa_seq_t *seqs;
	bwa_seqio_t *ks;
	clock_t t;
	bntseq_t *bns;
	FILE *fp_sa;
	gap_opt_t opt;
	char magic[4];

	// initialization
	bwase_initialize();
	bns = bns_restore(prefix);
	srand48(bns->seed);
	fp_sa = xopen(fn_sa, "r");

	m_aln = 0;
	err_fread_noeof(magic, 1, 4, fp_sa);
	if (strncmp(magic, SAI_MAGIC, 4) != 0) {
		fprintf(stderr, "[E::%s] Unmatched SAI magic. Please re-run `aln' with the same version of bwa.\n", __func__);
		exit(1);
	}
	err_fread_noeof(&opt, sizeof(gap_opt_t), 1, fp_sa);
	bwa_print_sam_hdr(bns, rg_line);
	// set ks
	ks = bwa_open_reads(opt.mode, fn_fa);
	// core loop
	while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt.mode, opt.trim_qual)) != 0) {
		tot_seqs += n_seqs;
		t = clock();

		// read alignment
		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p = seqs + i;
			int n_aln;
			err_fread_noeof(&n_aln, 4, 1, fp_sa);
			if (n_aln > m_aln) {
				m_aln = n_aln;
				aln = (bwt_aln1_t*)realloc(aln, sizeof(bwt_aln1_t) * m_aln);
			}
			err_fread_noeof(aln, sizeof(bwt_aln1_t), n_aln, fp_sa);
			bwa_aln2seq_core(n_aln, aln, p, 1, n_occ);
		}

		fprintf(stderr, "[bwa_aln_core] convert to sequence coordinate... ");
		bwa_cal_pac_pos(bns, prefix, n_seqs, seqs, opt.max_diff, opt.fnr); // forward bwt will be destroyed here
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		fprintf(stderr, "[bwa_aln_core] refine gapped alignments... ");
		bwa_refine_gapped(bns, n_seqs, seqs, 0);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC); t = clock();

		fprintf(stderr, "[bwa_aln_core] print alignments... ");
		for (i = 0; i < n_seqs; ++i)
			bwa_print_sam1(bns, seqs + i, 0, opt.mode, opt.max_top2);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);

		bwa_free_read_seq(n_seqs, seqs);
		fprintf(stderr, "[bwa_aln_core] %lld sequences have been processed.\n", tot_seqs);
	}

	// destroy
	bwa_seq_close(ks);
	bns_destroy(bns);
	err_fclose(fp_sa);
	free(aln);
}

int bwa_sai2sam_se(int argc, char *argv[])
{
	int c, n_occ = 3;
	char *prefix, *rg_line = 0;
	while ((c = getopt(argc, argv, "hn:f:r:")) >= 0) {
		switch (c) {
		case 'h': break;
		case 'r':
			if ((rg_line = bwa_set_rg(optarg)) == 0) return 1;
			break;
		case 'n': n_occ = atoi(optarg); break;
		case 'f': xreopen(optarg, "w", stdout); break;
		default: return 1;
		}
	}

	if (optind + 3 > argc) {
		fprintf(stderr, "Usage: bwa samse [-n max_occ] [-f out.sam] [-r RG_line] <prefix> <in.sai> <in.fq>\n");
		return 1;
	}
	if ((prefix = bwa_idx_infer_prefix(argv[optind])) == 0) {
		fprintf(stderr, "[%s] fail to locate the index\n", __func__);
		return 1;
	}
	bwa_sai2sam_se_core(prefix, argv[optind+1], argv[optind+2], n_occ, rg_line);
	free(prefix);
	return 0;
}
