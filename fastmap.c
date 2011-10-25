#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include "bntseq.h"
#include "bwt.h"
#include "kvec.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	int qbeg, len;
	int64_t tbeg;
} fmhit_t;

#define hit_qlt(a, b) ((a).qbeg < (b).qbeg)
#define hit_tlt(a, b) ((a).tbeg < (b).tbeg)

#include "ksort.h"
KSORT_INIT(hitq, fmhit_t, hit_qlt)
KSORT_INIT(hitt, fmhit_t, hit_tlt)

extern unsigned char nst_nt4_table[256];

typedef struct { size_t n, m; fmhit_t *a; } fmhit_v;

static uint64_t cluster_qhits(fmhit_v *hits, int beg, int end, int max_dist, int *_score)
{
	int qend, score, max = 0, cbeg = beg;
	size_t i;
	uint64_t max_cluster = 0;
	ks_introsort(hitq, end - beg, hits->a + beg);
	qend = hits->a[beg].qbeg + hits->a[beg].len;
	score = hits->a[beg].len;
	for (i = beg + 1; i < end; ++i) {
		fmhit_t *p = &hits->a[i];
		if (p->qbeg - qend > max_dist) {
			if (score > max) max = score, max_cluster = (uint64_t)cbeg<<32 | i;
			score = 0; cbeg = i;
		}
		score += p->len;
		qend = p->qbeg + p->len;
	}
	if (score > max) max = score, max_cluster = (uint64_t)cbeg<<32 | i;
	*_score = score;
	return max_cluster;
}

static int cluster_hits(fmhit_v *hits, int max_dist)
{
	size_t i;
	int64_t tend;
	uint64_t cluster, max_cluster = 0;
	int j, score, max = 0, max2 = 0, cbeg = 0, cend;
	ks_introsort(hitt, hits->n, hits->a);
	tend = hits->a[0].tbeg + hits->a[0].len;
	for (i = 1; i < hits->n; ++i) {
		fmhit_t *p = &hits->a[i];
		if (p->tbeg - tend > max_dist) {
			cluster = cluster_qhits(hits, cbeg, i, max_dist, &score);
			if (score > max) max2 = max, max = score, max_cluster = cluster;
			else if (score > max2) max2 = score;
			cbeg = i;
		}
		tend = p->tbeg + p->len;
	}
	cluster = cluster_qhits(hits, cbeg, i, max_dist, &score);
	if (score > max) max2 = max, max = score, max_cluster = cluster;
	else if (score > max2) max2 = score;
	cbeg = max_cluster>>32; cend = (uint32_t)max_cluster;
	for (i = 0, j = cbeg; j < cend; ++j) hits->a[i++] = hits->a[j];
	hits->n = i;
	return (int)(200.0 * (max - max2) / max + .499);
}

static int fake_cigar(const bntseq_t *bns, fmhit_v *hits, int beg, int end, int len, uint32_t *cigar, int64_t *pos, int *is_rev)
{
	size_t i;
	int qbeg, qend, n_cigar = 0;
	int64_t tbeg, tend, tmp;
	qbeg = len; qend = 0;
	tbeg = bns->l_pac<<1; tend = 0;
	for (i = beg; i < end; ++i) {
		fmhit_t *p = &hits->a[i];
		if (p->qbeg < qbeg) qbeg = p->qbeg;
		if (p->qbeg + p->len > qend) qend = p->qbeg + p->len;
		if (p->tbeg < tbeg) tbeg = p->tbeg;
		if (p->tbeg + p->len > qend) tend = p->tbeg + p->len;
	}
	if (tbeg >= bns->l_pac) {
		tmp = tend; tend = bns->l_pac*2 - tbeg; tbeg = bns->l_pac*2 - tmp;
		tmp = qend; qend = len - qbeg; qbeg = len - tmp;
		*is_rev = 1;
	} else *is_rev = 0;
	*pos = tbeg;
	if (qbeg) cigar[n_cigar++] = qbeg<<4|4;
	if (tend - tbeg < qend - qbeg) { // reference is shorter
		cigar[n_cigar++] = (uint32_t)(tend - tbeg)<<4 | 0;
		cigar[n_cigar++] = (uint32_t)((qend - qbeg) - (tend - tbeg))<<4 | 2;
	} else if (tend - tbeg > qend - qbeg) { // query is shorter
		cigar[n_cigar++] = (uint32_t)(qend - qbeg)<<4 | 0;
		cigar[n_cigar++] = (uint32_t)((tend - tbeg) - (qend - qbeg))<<4 | 1;
	} else cigar[n_cigar++] = (uint32_t)(qend - qbeg)<<4 | 0;
	if (len > qend) cigar[n_cigar++] = (uint32_t)(len - qend)<<4|4;
	return n_cigar;
}

int main_fastmap(int argc, char *argv[])
{
	int c, i, min_iwidth = 3, min_len = 17, max_dist = 100, mem_only = 0;
	kseq_t *seq;
	bwtint_t k;
	gzFile fp;
	bwt_t *bwt;
	bntseq_t *bns;
	bwtintv_v a[3], mem, *tvec[3];
	fmhit_v hits;
	uint32_t cigar[1024];

	while ((c = getopt(argc, argv, "w:l:d:p")) >= 0) {
		switch (c) {
			case 'w': min_iwidth = atoi(optarg); break;
			case 'l': min_len = atoi(optarg); break;
			case 'd': max_dist = atoi(optarg); break;
			case 'p': mem_only = 1; break;
		}
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "bwa fastmap <idxbase> <in.fq>\n");
		return 1;
	}

	fp = gzopen(argv[optind + 1], "r");
	seq = kseq_init(fp);
	{ // load the packed sequences, BWT and SA
		char *tmp = calloc(strlen(argv[optind]) + 5, 1);
		strcat(strcpy(tmp, argv[optind]), ".bwt");
		bwt = bwt_restore_bwt(tmp);
		strcat(strcpy(tmp, argv[optind]), ".sa");
		bwt_restore_sa(tmp, bwt);
		free(tmp);
		bns = bns_restore(argv[optind]);
	}
	for (i = 0; i < 3; ++i) { // initiate the temporary array
		kv_init(a[i]);
		tvec[i] = &a[i];
	}
	kv_init(mem); kv_init(hits);
	while (kseq_read(seq) >= 0) {
		for (i = 0; i < seq->seq.l; ++i)
			seq->seq.s[i] = nst_nt4_table[(int)seq->seq.s[i]];
		bwt_smem(bwt, seq->seq.l, (uint8_t*)seq->seq.s, &mem, tvec);
		if (mem_only) printf("SQ\t%s\t%ld\n", seq->name.s, seq->seq.l);
		for (i = 0, hits.n = 0; i < mem.n; ++i) {
			bwtintv_t *p = &mem.a[i];
			if ((uint32_t)p->info - (p->info>>32) < min_len) continue;
			if (mem_only) printf("EM\t%d\t%d\t%ld", (uint32_t)(p->info>>32), (uint32_t)p->info, (long)p->x[2]);
			if (p->x[2] <= min_iwidth) {
				for (k = 0; k < p->x[2]; ++k) {
					fmhit_t z;
					z.tbeg = bwt_sa(bwt, p->x[0] + k);
					z.qbeg = p->info>>32;
					z.len  = (uint32_t)p->info - z.qbeg;
					kv_push(fmhit_t, hits, z);
					if (mem_only) {
						int is_rev, ref_id;
						int64_t pos = bns_depos(bns, z.tbeg, &is_rev);
						if (is_rev) pos -= z.len - 1;
						bns_cnt_ambi(bns, pos, z.len, &ref_id);
						printf("\t%s:%c%ld", bns->anns[ref_id].name, "+-"[is_rev], (long)(pos - bns->anns[ref_id].offset) + 1);
					}
				}
			}
			if (mem_only) putchar('\n');
		}
		if (!mem_only) {
			int64_t pos;
			int n_cigar, is_rev, ref_id, mapq;
			mapq = cluster_hits(&hits, max_dist);
			n_cigar = fake_cigar(bns, &hits, 0, hits.n, seq->seq.l, cigar, &pos, &is_rev);
			bns_cnt_ambi(bns, pos, 1, &ref_id);
			printf("%s\t%d\t%s\t%ld\t%d\t", seq->name.s, is_rev?16:0,  bns->anns[ref_id].name, (long)(pos - bns->anns[ref_id].offset) + 1, mapq);
			for (i = 0; i < n_cigar; ++i)
				printf("%d%c", cigar[i]>>4, "MIDNSHP"[cigar[i]&0xf]);
			printf("\t*\t0\t0\t*\t*\n");
		} else puts("//");
	}

	free(mem.a);
	for (i = 0; i < 3; ++i) free(a[i].a);
	bns_destroy(bns);
	bwt_destroy(bwt);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}
