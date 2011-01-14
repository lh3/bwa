/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "bntseq.h"
#include "main.h"
#include "utils.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

unsigned char nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void bns_dump(const bntseq_t *bns, const char *prefix)
{
	char str[1024];
	FILE *fp;
	int i;
	{ // dump .ann
		strcpy(str, prefix); strcat(str, ".ann");
		fp = xopen(str, "w");
		fprintf(fp, "%lld %d %u\n", (long long)bns->l_pac, bns->n_seqs, bns->seed);
		for (i = 0; i != bns->n_seqs; ++i) {
			bntann1_t *p = bns->anns + i;
			fprintf(fp, "%d %s", p->gi, p->name);
			if (p->anno[0]) fprintf(fp, " %s\n", p->anno);
			else fprintf(fp, "\n");
			fprintf(fp, "%lld %d %d\n", (long long)p->offset, p->len, p->n_ambs);
		}
		fclose(fp);
	}
	{ // dump .amb
		strcpy(str, prefix); strcat(str, ".amb");
		fp = xopen(str, "w");
		fprintf(fp, "%lld %d %u\n", (long long)bns->l_pac, bns->n_seqs, bns->n_holes);
		for (i = 0; i != bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			fprintf(fp, "%lld %d %c\n", (long long)p->offset, p->len, p->amb);
		}
		fclose(fp);
	}
}

bntseq_t *bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename)
{
	char str[1024];
	FILE *fp;
	bntseq_t *bns;
	long long xx;
	int i;
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	{ // read .ann
		fp = xopen(ann_filename, "r");
		fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
		bns->l_pac = xx;
		bns->anns = (bntann1_t*)calloc(bns->n_seqs, sizeof(bntann1_t));
		for (i = 0; i < bns->n_seqs; ++i) {
			bntann1_t *p = bns->anns + i;
			char *q = str;
			int c;
			// read gi and sequence name
			fscanf(fp, "%u%s", &p->gi, str);
			p->name = strdup(str);
			// read fasta comments 
			while ((c = fgetc(fp)) != '\n' && c != EOF) *q++ = c;
			*q = 0;
			if (q - str > 1) p->anno = strdup(str + 1); // skip leading space
			else p->anno = strdup("");
			// read the rest
			fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
			p->offset = xx;
		}
		fclose(fp);
	}
	{ // read .amb
		int64_t l_pac;
		int32_t n_seqs;
		fp = xopen(amb_filename, "r");
		fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
		l_pac = xx;
		xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
		bns->ambs = (bntamb1_t*)calloc(bns->n_holes, sizeof(bntamb1_t));
		for (i = 0; i < bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			fscanf(fp, "%lld%d%s", &xx, &p->len, str);
			p->offset = xx;
			p->amb = str[0];
		}
		fclose(fp);
	}
	{ // open .pac
		bns->fp_pac = xopen(pac_filename, "rb");
	}
	return bns;
}

bntseq_t *bns_restore(const char *prefix)
{  
	char ann_filename[1024], amb_filename[1024], pac_filename[1024];
	strcat(strcpy(ann_filename, prefix), ".ann");
	strcat(strcpy(amb_filename, prefix), ".amb");
	strcat(strcpy(pac_filename, prefix), ".pac");
	return bns_restore_core(ann_filename, amb_filename, pac_filename);
}

void bns_destroy(bntseq_t *bns)
{
	if (bns == 0) return;
	else {
		int i;
		if (bns->fp_pac) fclose(bns->fp_pac);
		free(bns->ambs);
		for (i = 0; i < bns->n_seqs; ++i) {
			free(bns->anns[i].name);
			free(bns->anns[i].anno);
		}
		free(bns->anns);
		free(bns);
	}
}

void bns_fasta2bntseq(gzFile fp_fa, const char *prefix)
{
	kseq_t *seq;
	char name[1024];
	bntseq_t *bns;
	bntamb1_t *q;
	int l_buf;
	unsigned char buf[0x10000];
	int32_t m_seqs, m_holes, l, i;
	FILE *fp;

	// initialization
	seq = kseq_init(fp_fa);
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	bns->seed = 11; // fixed seed for random generator
	srand48(bns->seed);
	m_seqs = m_holes = 8;
	bns->anns = (bntann1_t*)calloc(m_seqs, sizeof(bntann1_t));
	bns->ambs = (bntamb1_t*)calloc(m_holes, sizeof(bntamb1_t));
	q = bns->ambs;
	l_buf = 0;
	strcpy(name, prefix); strcat(name, ".pac");
	fp = xopen(name, "wb");
	memset(buf, 0, 0x10000);
	// read sequences
	while ((l = kseq_read(seq)) >= 0) {
		bntann1_t *p;
		int lasts;
		if (bns->n_seqs == m_seqs) {
			m_seqs <<= 1;
			bns->anns = (bntann1_t*)realloc(bns->anns, m_seqs * sizeof(bntann1_t));
		}
		p = bns->anns + bns->n_seqs;
		p->name = strdup((char*)seq->name.s);
		p->anno = seq->comment.s? strdup((char*)seq->comment.s) : strdup("(null)");
		p->gi = 0; p->len = l;
		p->offset = (bns->n_seqs == 0)? 0 : (p-1)->offset + (p-1)->len;
		p->n_ambs = 0;
		for (i = 0, lasts = 0; i < l; ++i) {
			int c = nst_nt4_table[(int)seq->seq.s[i]];
			if (c >= 4) { // N
				if (lasts == seq->seq.s[i]) { // contiguous N
					++q->len;
				} else {
					if (bns->n_holes == m_holes) {
						m_holes <<= 1;
						bns->ambs = (bntamb1_t*)realloc(bns->ambs, m_holes * sizeof(bntamb1_t));
					}
					q = bns->ambs + bns->n_holes;
					q->len = 1;
					q->offset = p->offset + i;
					q->amb = seq->seq.s[i];
					++p->n_ambs;
					++bns->n_holes;
				}
			}
			lasts = seq->seq.s[i];
			{ // fill buffer
				if (c >= 4) c = lrand48()&0x3;
				if (l_buf == 0x40000) {
					fwrite(buf, 1, 0x10000, fp);
					memset(buf, 0, 0x10000);
					l_buf = 0;
				}
				buf[l_buf>>2] |= c << ((3 - (l_buf&3)) << 1);
				++l_buf;
			}
		}
		++bns->n_seqs;
		bns->l_pac += seq->seq.l;
	}
	xassert(bns->l_pac, "zero length sequence.");
	{ // finalize .pac file
		ubyte_t ct;
		fwrite(buf, 1, (l_buf>>2) + ((l_buf&3) == 0? 0 : 1), fp);
		// the following codes make the pac file size always (l_pac/4+1+1)
		if (bns->l_pac % 4 == 0) {
			ct = 0;
			fwrite(&ct, 1, 1, fp);
		}
		ct = bns->l_pac % 4;
		fwrite(&ct, 1, 1, fp);
		// close .pac file
		fclose(fp);
	}
	bns_dump(bns, prefix);
	bns_destroy(bns);
	kseq_destroy(seq);
}

int bwa_fa2pac(int argc, char *argv[])
{
	gzFile fp;
	if (argc < 2) {
		fprintf(stderr, "Usage: bwa fa2pac <in.fasta> [<out.prefix>]\n");
		return 1;
	}
	fp = xzopen(argv[1], "r");
	bns_fasta2bntseq(fp, (argc < 3)? argv[1] : argv[2]);
	gzclose(fp);
	return 0;
}

int bns_coor_pac2real(const bntseq_t *bns, int64_t pac_coor, int len, int32_t *real_seq)
{
	int left, mid, right, nn;
	if (pac_coor >= bns->l_pac)
		err_fatal("bns_coor_pac2real", "bug! Coordinate is longer than sequence (%lld>=%lld).", pac_coor, bns->l_pac);
	// binary search for the sequence ID. Note that this is a bit different from the following one...
	left = 0; mid = 0; right = bns->n_seqs;
	while (left < right) {
		mid = (left + right) >> 1;
		if (pac_coor >= bns->anns[mid].offset) {
			if (mid == bns->n_seqs - 1) break;
			if (pac_coor < bns->anns[mid+1].offset) break;
			left = mid + 1;
		} else right = mid;
	}
	*real_seq = mid;
	// binary search for holes
	left = 0; right = bns->n_holes; nn = 0;
	while (left < right) {
		int64_t mid = (left + right) >> 1;
		if (pac_coor >= bns->ambs[mid].offset + bns->ambs[mid].len) left = mid + 1;
		else if (pac_coor + len <= bns->ambs[mid].offset) right = mid;
		else { // overlap
			if (pac_coor >= bns->ambs[mid].offset) {
				nn += bns->ambs[mid].offset + bns->ambs[mid].len < pac_coor + len?
					bns->ambs[mid].offset + bns->ambs[mid].len - pac_coor : len;
			} else {
				nn += bns->ambs[mid].offset + bns->ambs[mid].len < pac_coor + len?
					bns->ambs[mid].len : len - (bns->ambs[mid].offset - pac_coor);
			}
			break;
		}
	}
	return nn;
}
