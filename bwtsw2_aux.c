#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif
#include "bntseq.h"
#include "bwt_lite.h"
#include "utils.h"
#include "bwtsw2.h"
#include "stdaln.h"
#include "kstring.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "ksort.h"
#define __left_lt(a, b) ((a).end > (b).end)
KSORT_INIT(hit, bsw2hit_t, __left_lt)

extern unsigned char nst_nt4_table[256];

unsigned char nt_comp_table[256] = {
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','T','V','G', 'H','N','N','C', 'D','N','N','M', 'N','K','N','N',
	'N','N','Y','S', 'A','N','B','W', 'X','R','N','N', 'N','N','N','N',
	'n','t','v','g', 'h','n','n','c', 'd','n','n','m', 'n','k','n','n',
	'n','n','y','s', 'a','n','b','w', 'x','r','n','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N',
	'N','N','N','N', 'N','N','N','N', 'N','N','N','N', 'N','N','N','N'
};

extern int bsw2_resolve_duphits(const bwt_t *bwt, bwtsw2_t *b, int IS);
extern int bsw2_resolve_query_overlaps(bwtsw2_t *b, float mask_level);

bsw2opt_t *bsw2_init_opt()
{
	bsw2opt_t *o = (bsw2opt_t*)calloc(1, sizeof(bsw2opt_t));
	o->a = 1; o->b = 3; o->q = 5; o->r = 2; o->t = 30;
	o->bw = 50;
	o->z = 1; o->is = 3; o->t_seeds = 5; o->hard_clip = 0;
	o->mask_level = 0.50f; o->yita = 5.5f; o->coef = 5.5f;
	o->qr = o->q + o->r; o->n_threads = 1; o->chunk_size = 10000000;
	return o;
}

void bsw2_destroy(bwtsw2_t *b)
{
	int i;
	if (b == 0) return;
	if (b->cigar)
		for (i = 0; i < b->n; ++i) free(b->cigar[i]);
	free(b->cigar); free(b->n_cigar); free(b->hits);
	free(b);
}

#define __gen_ap(par, opt) do {									\
		int i;													\
		for (i = 0; i < 25; ++i) (par).matrix[i] = -(opt)->b;	\
		for (i = 0; i < 4; ++i) (par).matrix[i*5+i] = (opt)->a; \
		(par).gap_open = (opt)->q; (par).gap_ext = (opt)->r;	\
		(par).gap_end = (opt)->r;								\
		(par).row = 5; (par).band_width = opt->bw;				\
	} while (0)

#define __rpac(pac, l, i) (pac[(l-i-1)>>2] >> (~(l-i-1)&3)*2 & 0x3)

void bsw2_extend_left(const bsw2opt_t *opt, bwtsw2_t *b, uint8_t *_query, int lq, uint8_t *pac, uint32_t l_pac, int is_rev, uint8_t *_mem)
{
	int i, matrix[25];
	bwtint_t k;
	uint8_t *target = 0, *query;
	AlnParam par;

	par.matrix = matrix;
	__gen_ap(par, opt);
	query = calloc(lq, 1);
	// sort according to the descending order of query end
	ks_introsort(hit, b->n, b->hits);
	target = calloc(((lq + 1) / 2 * opt->a + opt->r) / opt->r + lq, 1);
	// reverse _query
	for (i = 0; i < lq; ++i) query[lq - i - 1] = _query[i];
	// core loop
	for (i = 0; i < b->n; ++i) {
		bsw2hit_t *p = b->hits + i;
		int lt = ((p->beg + 1) / 2 * opt->a + opt->r) / opt->r + lq;
		int score, j;
		path_t path;
		p->n_seeds = 1;
		if (p->l || p->k == 0) continue;
		for (j = score = 0; j < i; ++j) {
			bsw2hit_t *q = b->hits + j;
			if (q->beg <= p->beg && q->k <= p->k && q->k + q->len >= p->k + p->len) {
				if (q->n_seeds < (1<<14) - 2) ++q->n_seeds;
				++score;
			}
		}
		if (score) continue;
		if (lt > p->k) lt = p->k;
		if (is_rev) {
			for (k = p->k - 1, j = 0; k > 0 && j < lt; --k) // FIXME: k=0 not considered!
				target[j++] = __rpac(pac, l_pac, k);
		} else {
			for (k = p->k - 1, j = 0; k > 0 && j < lt; --k) // FIXME: k=0 not considered!
				target[j++] = pac[k>>2] >> (~k&3)*2 & 0x3;
		}
		lt = j;
		score = aln_extend_core(target, lt, query + lq - p->beg, p->beg, &par, &path, 0, p->G, _mem);
		if (score > p->G) { // extensible
			p->G = score;
			p->len += path.i;
			p->beg -= path.j;
			p->k -= path.i;
		}
	}
	free(query); free(target);
}

void bsw2_extend_rght(const bsw2opt_t *opt, bwtsw2_t *b, uint8_t *query, int lq, uint8_t *pac, uint32_t l_pac, int is_rev, uint8_t *_mem)
{
	int i, matrix[25];
	uint32_t k;
	uint8_t *target;
	AlnParam par;
	
	par.matrix = matrix;
	__gen_ap(par, opt);
	target = calloc(((lq + 1) / 2 * opt->a + opt->r) / opt->r + lq, 1);
	for (i = 0; i < b->n; ++i) {
		bsw2hit_t *p = b->hits + i;
		int lt = ((lq - p->beg + 1) / 2 * opt->a + opt->r) / opt->r + lq;
		int j, score;
		path_t path;
		if (p->l) continue;
		if (is_rev) {
			for (k = p->k, j = 0; k < p->k + lt && k < l_pac; ++k)
				target[j++] = __rpac(pac, l_pac, k);
		} else {
			for (k = p->k, j = 0; k < p->k + lt && k < l_pac; ++k)
				target[j++] = pac[k>>2] >> (~k&3)*2 & 0x3;
		}
		lt = j;
		score = aln_extend_core(target, lt, query + p->beg, lq - p->beg, &par, &path, 0, 1, _mem);
//		if (score < p->G) fprintf(stderr, "[bsw2_extend_hits] %d < %d\n", score, p->G);
		if (score >= p->G) {
			p->G = score;
			p->len = path.i;
			p->end = path.j + p->beg;
		}
	}
	free(target);
}

/* generate CIGAR array(s) in b->cigar[] */
static void gen_cigar(const bsw2opt_t *opt, int lq, uint8_t *seq[2], uint8_t *pac, bwtsw2_t *b)
{
	uint8_t *target;
	int i, matrix[25];
	AlnParam par;
	path_t *path;

	par.matrix = matrix;
	__gen_ap(par, opt);
	i = ((lq + 1) / 2 * opt->a + opt->r) / opt->r + lq; // maximum possible target length
	target = calloc(i, 1);
	path = calloc(i + lq, sizeof(path_t));
	// memory clean up for b
	if (b->n < b->max) {
		b->max = b->n;
		b->hits = realloc(b->hits, b->n * sizeof(bsw2hit_t));
	}
	if (b->cigar) free(b->cigar);
	if (b->n_cigar) free(b->n_cigar);
	b->cigar = (uint32_t**)calloc(b->max, sizeof(void*));
	b->n_cigar = (int*)calloc(b->max, sizeof(int));
	// generate CIGAR
	for (i = 0; i < b->n; ++i) {
		bsw2hit_t *p = b->hits + i;
		uint8_t *query;
		uint32_t k;
		int score, path_len, beg, end;
		if (p->l) continue;
		beg = (p->flag & 0x10)? lq - p->end : p->beg;
		end = (p->flag & 0x10)? lq - p->beg : p->end;
		query = seq[(p->flag & 0x10)? 1 : 0] + beg;
		for (k = p->k; k < p->k + p->len; ++k) // in principle, no out-of-boundary here
			target[k - p->k] = pac[k>>2] >> (~k&3)*2 & 0x3;
		score = aln_global_core(target, p->len, query, end - beg, &par, path, &path_len);
		b->cigar[i] = aln_path2cigar32(path, path_len, &b->n_cigar[i]);
		if (beg != 0 || end < lq) { // write soft clipping
			b->cigar[i] = realloc(b->cigar[i], 4 * (b->n_cigar[i] + 2));
			if (beg != 0) {
				memmove(b->cigar[i] + 1, b->cigar[i], b->n_cigar[i] * 4);
				b->cigar[i][0] = beg<<4 | 4;
				++b->n_cigar[i];
			}
			if (end < lq) {
				b->cigar[i][b->n_cigar[i]] = (lq - end)<<4 | 4;
				++b->n_cigar[i];
			}
		}
	}
	free(target); free(path);
}

/* this is for the debugging purpose only */
void bsw2_debug_hits(const bwtsw2_t *b)
{
	int i;
	printf("# raw hits: %d\n", b->n);
	for (i = 0; i < b->n; ++i) {
		bsw2hit_t *p = b->hits + i;
		if (p->l == 0)
			printf("%d, %d, %d, %u, %u\n", p->G, p->beg, p->end, p->k, p->l);
	}
}

static void merge_hits(bwtsw2_t *b[2], int l, int is_reverse)
{
	int i;
	if (b[0]->n + b[1]->n > b[0]->max) {
		b[0]->max = b[0]->n + b[1]->n;
		b[0]->hits = realloc(b[0]->hits, b[0]->max * sizeof(bsw2hit_t));
	}
	for (i = 0; i < b[1]->n; ++i) {
		bsw2hit_t *p = b[0]->hits + b[0]->n + i;
		*p = b[1]->hits[i];
		if (is_reverse) {
			int x = p->beg;
			p->beg = l - p->end;
			p->end = l - x;
			p->flag |= 0x10;
		}
	}
	b[0]->n += b[1]->n;
	bsw2_destroy(b[1]);
	b[1] = 0;
}
/* seq[0] is the forward sequence and seq[1] is the reverse complement. */
static bwtsw2_t *bsw2_aln1_core(const bsw2opt_t *opt, const bntseq_t *bns, uint8_t *pac, const bwt_t *target,
								int l, uint8_t *seq[2], int is_rev, bsw2global_t *pool)
{
	extern void bsw2_chain_filter(const bsw2opt_t *opt, int len, bwtsw2_t *b[2]);
	bwtsw2_t *b[2], **bb[2];
	int k;
	for (k = 0; k < 2; ++k) {
		bwtl_t *query = bwtl_seq2bwtl(l, seq[k]);
		bb[k] = bsw2_core(opt, query, target, pool);
		bwtl_destroy(query);
	}
	b[0] = bb[0][1]; b[1] = bb[1][1]; // bb[*][1] are "narrow SA hits"
	bsw2_chain_filter(opt, l, b);
	for (k = 0; k < 2; ++k) {
		bsw2_extend_left(opt, bb[k][1], seq[k], l, pac, bns->l_pac, is_rev, pool->aln_mem);
		merge_hits(bb[k], l, 0); // bb[k][1] is merged to bb[k][0] here
		bsw2_resolve_duphits(0, bb[k][0], 0);
		bsw2_extend_rght(opt, bb[k][0], seq[k], l, pac, bns->l_pac, is_rev, pool->aln_mem);
		b[k] = bb[k][0];
		free(bb[k]);		
	}
	merge_hits(b, l, 1); // again, b[1] is merged to b[0]
	bsw2_resolve_query_overlaps(b[0], opt->mask_level);
	return b[0];
}

/* set ->flag to records the origin of the hit (to forward bwt or reverse bwt) */
static void flag_fr(bwtsw2_t *b[2])
{
	int i, j;
	for (i = 0; i < b[0]->n; ++i) {
		bsw2hit_t *p = b[0]->hits + i;
		p->flag |= 0x10000;
	}
	for (i = 0; i < b[1]->n; ++i) {
		bsw2hit_t *p = b[1]->hits + i;
		p->flag |= 0x20000;
	}
	for (i = 0; i < b[0]->n; ++i) {
		bsw2hit_t *p = b[0]->hits + i;
		for (j = 0; j < b[1]->n; ++j) {
			bsw2hit_t *q = b[1]->hits + j;
			if (q->beg == p->beg && q->end == p->end && q->k == p->k && q->len == p->len && q->G == p->G) {
				q->flag |= 0x30000; p->flag |= 0x30000;
				break;
			}
		}
	}
}

typedef struct {
	int l, tid;
	char *name, *seq, *qual, *sam;
} bsw2seq1_t;

typedef struct {
	int n, max;
	bsw2seq1_t *seq;
} bsw2seq_t;

#ifdef HAVE_PTHREAD
static pthread_mutex_t g_dbwtsw_lock = PTHREAD_MUTEX_INITIALIZER;
#endif

static int fix_cigar(const char *qname, const bntseq_t *bns, bsw2hit_t *p, int n_cigar, uint32_t *cigar)
{
	// FIXME: this routine does not work if the query bridge three reference sequences
	int32_t coor, refl, lq;
	int x, y, i, seqid;
	bns_coor_pac2real(bns, p->k, p->len, &seqid);
	coor = p->k - bns->anns[seqid].offset;
	refl = bns->anns[seqid].len;
	x = coor; y = 0;
	// test if the alignment goes beyond the boundary
	for (i = 0; i < n_cigar; ++i) {
		int op = cigar[i]&0xf, ln = cigar[i]>>4;
		if (op == 1 || op == 4 || op == 5) y += ln;
		else if (op == 2) x += ln;
		else x += ln, y += ln;
	}
	lq = y; // length of the query sequence
	if (x > refl) { // then fix it
		int j, nc, mq[2], nlen[2];
		uint32_t *cn, kk = 0;
		nc = mq[0] = mq[1] = nlen[0] = nlen[1] = 0;
		cn = calloc(n_cigar + 3, 4);
		x = coor; y = 0;
		for (i = j = 0; i < n_cigar; ++i) {
			int op = cigar[i]&0xf, ln = cigar[i]>>4;
			if (op == 4 || op == 5 || op == 1) { // ins or clipping
				y += ln;
				cn[j++] = cigar[i];
			} else if (op == 2) { // del
				if (x + ln >= refl && nc == 0) {
					cn[j++] = (uint32_t)(lq - y)<<4 | 4;
					nc = j;
					cn[j++] = (uint32_t)y<<4 | 4;
					kk = p->k + (x + ln - refl);
					nlen[0] = x - coor;
					nlen[1] = p->len - nlen[0] - ln;
				} else cn[j++] = cigar[i];
				x += ln;
			} else if (op == 0) { // match
				if (x + ln >= refl && nc == 0) {
					// FIXME: not consider a special case where a split right between M and I
					cn[j++] = (uint32_t)(refl - x)<<4 | 0; // write M
					cn[j++] = (uint32_t)(lq - y - (refl - x))<<4 | 4; // write S
					nc = j;
					mq[0] += refl - x;
					cn[j++] = (uint32_t)(y + (refl - x))<<4 | 4;
					if (x + ln - refl) cn[j++] = (uint32_t)(x + ln - refl)<<4 | 0;
					mq[1] += x + ln - refl;
					kk = bns->anns[seqid].offset + refl;
					nlen[0] = refl - coor;
					nlen[1] = p->len - nlen[0];
				} else {
					cn[j++] = cigar[i];
					mq[nc?1:0] += ln;
				}
				x += ln; y += ln;
			}
		}
		if (mq[0] > mq[1]) { // then take the first alignment
			n_cigar = nc;
			memcpy(cigar, cn, 4 * nc);
			p->len = nlen[0];
		} else {
			p->k = kk; p->len = nlen[1];
			n_cigar = j - nc;
			memcpy(cigar, cn + nc, 4 * (j - nc));
		}
		free(cn);
	}
	return n_cigar;
}

/* generate SAM lines for a sequence in ks with alignment stored in
 * b. ks->name and ks->seq will be freed and set to NULL in the end. */
static void print_hits(const bntseq_t *bns, const bsw2opt_t *opt, bsw2seq1_t *ks, bwtsw2_t *b)
{
	int i, k;
	kstring_t str;
	memset(&str, 0, sizeof(kstring_t));
	if (b == 0 || b->n == 0) { // no hits
		ksprintf(&str, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t", ks->name);
		for (i = 0; i < ks->l; ++i) kputc(ks->seq[i], &str);
		if (ks->qual) {
			kputc('\t', &str);
			for (i = 0; i < ks->l; ++i) kputc(ks->qual[i], &str);
		} else kputs("\t*", &str);
		kputc('\n', &str);
	}
	for (i = 0; b && i < b->n; ++i) {
		bsw2hit_t *p = b->hits + i;
		int32_t seqid = -1, coor = -1;
		int j, qual, nn = 0;
		int beg, end;
		if (p->l == 0) {
			b->n_cigar[i] = fix_cigar(ks->name, bns, p, b->n_cigar[i], b->cigar[i]);
			nn = bns_coor_pac2real(bns, p->k, p->len, &seqid);
			coor = p->k - bns->anns[seqid].offset;
		}
		ksprintf(&str, "%s\t%d", ks->name, p->flag&0x10);
		ksprintf(&str, "\t%s\t%d", seqid>=0? bns->anns[seqid].name : "*", coor + 1);
		if (p->l == 0) {
			{ // estimate mapping quality
				float c = 1.0;	
				int subo = p->G2 > opt->t? p->G2 : opt->t;
				if (p->flag>>16 == 1 || p->flag>>16 == 2) c *= .5;
				if (p->n_seeds < 2) c *= .2;
				qual = (int)(c * (p->G - subo) * (250.0 / p->G + 0.03 / opt->a) + .499);
				if (qual > 250) qual = 250;
				if (p->flag&1) qual = 0;
			}
			ksprintf(&str, "\t%d\t", qual);
			for (k = 0; k < b->n_cigar[i]; ++k)
				ksprintf(&str, "%d%c", b->cigar[i][k]>>4, (opt->hard_clip? "MIDNHHP" : "MIDNSHP")[b->cigar[i][k]&0xf]);
		} else ksprintf(&str, "\t0\t*");
		ksprintf(&str, "\t*\t0\t0\t");
		beg = 0; end = ks->l;
		if (opt->hard_clip) {
			if ((b->cigar[i][0]&0xf) == 4) beg += b->cigar[i][0]>>4;
			if ((b->cigar[i][b->n_cigar[i]-1]&0xf) == 4) end -= b->cigar[i][b->n_cigar[i]-1]>>4;
		}
		for (j = beg; j < end; ++j) {
			if (p->flag&0x10) kputc(nt_comp_table[(int)ks->seq[ks->l - 1 - j]], &str);
			else kputc(ks->seq[j], &str);
		}
		if (ks->qual) {
			kputc('\t', &str);
			for (j = beg; j < end; ++j) {
				if (p->flag&0x10) kputc(ks->qual[ks->l - 1 - j], &str);
				else kputc(ks->qual[j], &str);
			}
		} else ksprintf(&str, "\t*");
		ksprintf(&str, "\tAS:i:%d\tXS:i:%d\tXF:i:%d\tXE:i:%d\tXN:i:%d", p->G, p->G2, p->flag>>16, p->n_seeds, nn);
		if (p->l) ksprintf(&str, "\tXI:i:%d", p->l - p->k + 1);
		kputc('\n', &str);
	}
	ks->sam = str.s;
	free(ks->seq); ks->seq = 0;
	free(ks->qual); ks->qual = 0;
	free(ks->name); ks->name = 0;
}

/* Core routine to align reads in _seq. It is separated from
 * process_seqs() to realize multi-threading */ 
static void bsw2_aln_core(int tid, bsw2seq_t *_seq, const bsw2opt_t *_opt, const bntseq_t *bns, uint8_t *pac, bwt_t * const target[2])
{
	int x;
	bsw2opt_t opt = *_opt;
	bsw2global_t *pool = bsw2_global_init();
	for (x = 0; x < _seq->n; ++x) {
		bsw2seq1_t *p = _seq->seq + x;
		uint8_t *seq[2], *rseq[2];
		int i, l, k;
		bwtsw2_t *b[2];
		l = p->l;

#ifdef HAVE_PTHREAD
		if (_opt->n_threads > 1) {
			pthread_mutex_lock(&g_dbwtsw_lock);
			if (p->tid < 0) p->tid = tid;
			else if (p->tid != tid) {
				pthread_mutex_unlock(&g_dbwtsw_lock);
				continue;
			} // in pinciple else should not happen
			pthread_mutex_unlock(&g_dbwtsw_lock);
		}
#endif

		// set opt->t
		opt.t = _opt->t;
		if (opt.t < log(l) * opt.coef) opt.t = (int)(log(l) * opt.coef + .499);
		if (pool->max_l < l) { // then enlarge working space for aln_extend_core()
			int tmp = ((l + 1) / 2 * opt.a + opt.r) / opt.r + l;
			pool->max_l = l;
			pool->aln_mem = realloc(pool->aln_mem, (tmp + 2) * 24);
		}
		// set opt->bw
		opt.bw = _opt->bw;
		k = (l * opt.a - 2 * opt.q) / (2 * opt.r + opt.a);
		i = (l * opt.a - opt.a - opt.t) / opt.r;
		if (k > i) k = i;
		if (k < 1) k = 1; // I do not know if k==0 causes troubles
		opt.bw = _opt->bw < k? _opt->bw : k;
		// set seq[2] and rseq[2]
		seq[0] = calloc(l * 4, 1);
		seq[1] = seq[0] + l;
		rseq[0] = seq[1] + l; rseq[1] = rseq[0] + l;
		// convert sequences to 2-bit representation
		for (i = k = 0; i < l; ++i) {
			int c = nst_nt4_table[(int)p->seq[i]];
			if (c >= 4) { c = (int)(drand48() * 4); ++k; } // FIXME: ambiguous bases are not properly handled
			seq[0][i] = c;
			seq[1][l-1-i] = 3 - c;
			rseq[0][l-1-i] = c;
			rseq[1][i] = 3 - c;
		}
		if (l - k < opt.t) { // too few unambiguous bases
			print_hits(bns, &opt, p, 0);
			free(seq[0]); continue;
		}
		// alignment
		b[0] = bsw2_aln1_core(&opt, bns, pac, target[0], l, seq, 0, pool);
		for (k = 0; k < b[0]->n; ++k)
			if (b[0]->hits[k].n_seeds < opt.t_seeds) break;
		if (k < b[0]->n) {
			b[1] = bsw2_aln1_core(&opt, bns, pac, target[1], l, rseq, 1, pool);
			for (i = 0; i < b[1]->n; ++i) {
				bsw2hit_t *p = b[1]->hits + i;
				int x = p->beg;
				p->beg = l - p->end;
				p->end = l - x;
				if (p->l == 0) p->k = bns->l_pac - (p->k + p->len);
			}
			flag_fr(b);
			merge_hits(b, l, 0);
			bsw2_resolve_duphits(0, b[0], 0);
			bsw2_resolve_query_overlaps(b[0], opt.mask_level);
		} else b[1] = 0;
		// generate CIGAR and print SAM
		gen_cigar(&opt, l, seq, pac, b[0]);
		print_hits(bns, &opt, p, b[0]);
		// free
		free(seq[0]);
		bsw2_destroy(b[0]);
	}
	bsw2_global_destroy(pool);
}

#ifdef HAVE_PTHREAD
typedef struct {
	int tid;
	bsw2seq_t *_seq;
	const bsw2opt_t *_opt;
	const bntseq_t *bns;
	uint8_t *pac;
	bwt_t *target[2];
} thread_aux_t;

/* another interface to bsw2_aln_core() to facilitate pthread_create() */
static void *worker(void *data)
{
	thread_aux_t *p = (thread_aux_t*)data;
	bsw2_aln_core(p->tid, p->_seq, p->_opt, p->bns, p->pac, p->target);
	return 0;
}
#endif

/* process sequences stored in _seq, generate SAM lines for these
 * sequences and reset _seq afterwards. */
static void process_seqs(bsw2seq_t *_seq, const bsw2opt_t *opt, const bntseq_t *bns, uint8_t *pac, bwt_t * const target[2])
{
	int i;

#ifdef HAVE_PTHREAD
	if (opt->n_threads <= 1) {
		bsw2_aln_core(0, _seq, opt, bns, pac, target);
	} else {
		pthread_t *tid;
		pthread_attr_t attr;
		thread_aux_t *data;
		int j;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
		data = (thread_aux_t*)calloc(opt->n_threads, sizeof(thread_aux_t));
		tid = (pthread_t*)calloc(opt->n_threads, sizeof(pthread_t));
		for (j = 0; j < opt->n_threads; ++j) {
			thread_aux_t *p = data + j;
			p->tid = j; p->_seq = _seq; p->_opt = opt; p->bns = bns;
			p->pac = pac; p->target[0] = target[0]; p->target[1] = target[1];
			pthread_create(&tid[j], &attr, worker, p);
		}
		for (j = 0; j < opt->n_threads; ++j) pthread_join(tid[j], 0);
		free(data); free(tid);
	}
#else
	bsw2_aln_core(0, _seq, opt, bns, pac, target);
#endif

	// print and reset
	for (i = 0; i < _seq->n; ++i) {
		bsw2seq1_t *p = _seq->seq + i;
		if (p->sam) printf("%s", p->sam);
		free(p->name); free(p->seq); free(p->qual); free(p->sam);
		p->tid = -1; p->l = 0;
		p->name = p->seq = p->qual = p->sam = 0;
	}
	fflush(stdout);
	_seq->n = 0;
}

void bsw2_aln(const bsw2opt_t *opt, const bntseq_t *bns, bwt_t * const target[2], const char *fn)
{
	gzFile fp;
	kseq_t *ks;
	int l, size = 0;
	uint8_t *pac;
	bsw2seq_t *_seq;

	pac = calloc(bns->l_pac/4+1, 1);
	if (pac == 0) {
		fprintf(stderr, "[bsw2_aln] insufficient memory!\n");
		return;
	}
	for (l = 0; l < bns->n_seqs; ++l)
		printf("@SQ\tSN:%s\tLN:%d\n", bns->anns[l].name, bns->anns[l].len);
	fread(pac, 1, bns->l_pac/4+1, bns->fp_pac);
	fp = xzopen(fn, "r");
	ks = kseq_init(fp);
	_seq = calloc(1, sizeof(bsw2seq_t));
	while ((l = kseq_read(ks)) >= 0) {
		bsw2seq1_t *p;
		if (_seq->n == _seq->max) {
			_seq->max = _seq->max? _seq->max<<1 : 1024;
			_seq->seq = realloc(_seq->seq, _seq->max * sizeof(bsw2seq1_t));
		}
		p = &_seq->seq[_seq->n++];
		p->tid = -1;
		p->l = l;
		p->name = strdup(ks->name.s);
		p->seq = strdup(ks->seq.s);
		p->qual = ks->qual.l? strdup(ks->qual.s) : 0;
		p->sam = 0;
		size += l;
		if (size > opt->chunk_size) {
			fprintf(stderr, "[bsw2_aln] read %d sequences (%d bp)...\n", _seq->n, size);
			process_seqs(_seq, opt, bns, pac, target);
			size = 0;
		}
	}
	fprintf(stderr, "[bsw2_aln] read %d sequences (%d bp)...\n", _seq->n, size);
	process_seqs(_seq, opt, bns, pac, target);
	free(_seq->seq); free(_seq);
	kseq_destroy(ks);
	gzclose(fp);
	free(pac);
}
