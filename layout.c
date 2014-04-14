#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include "ksort.h"
#include "kvec.h"
#include "khash.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

static int lo_verbose = 3;

/***********************
 * Core data structure *
 ***********************/

#define LO_T_C1 8
#define LO_T_C2 12
#define LO_T_I  16

typedef struct {
	int min_ext, fuzzy_dist;
	int min_n_ovlp;
	float min_aln_ratio;
} lo_opt_t;

typedef struct {
	uint32_t dist:31, reduced:1;
	uint32_t id:31, ori:1;
} lo_nei_t;

typedef kvec_t(lo_nei_t) lo_nei_v;

#define nei_lt(a, b) ((a).dist < (b).dist)
KSORT_INIT(nei, lo_nei_t, nei_lt)

#define LO_VF_CONTAINED 0x1
#define LO_VF_LACK_OVLP 0x2

typedef struct {
	int id;
	int flag:16, state:16;
	char *name;
	lo_nei_v *nei[2];
} vertex_t;

typedef kvec_t(vertex_t) vertex_v;

#define lo_skipped(v) ((v)->flag & (LO_VF_CONTAINED|LO_VF_LACK_OVLP))

typedef struct {
	int type:16, reduced:16;
	int l[2], s[2], e[2], d[2]; // length, start and end
	float usc;
} edgeinfo_t;

KHASH_MAP_INIT_INT64(edge, edgeinfo_t)
typedef khash_t(edge) ehash_t;

KHASH_MAP_INIT_STR(name, int)
typedef khash_t(name) nhash_t;

typedef struct {
	nhash_t *n;
	ehash_t *e;
	vertex_v v;
} ograph_t;

void lo_opt_init(lo_opt_t *opt)
{
	opt->min_ext = 50;
	opt->min_aln_ratio = 0.9;
	opt->min_n_ovlp = 1;
	opt->fuzzy_dist = 100;
}

const char *lo_edge_label[] = {
	">>", "><", "<>", "<<",
	"??", "??", "??", "??",
	"C1", "??", "??", "??",
	"C2", "??", "??", "??", "IN" };

/**********
 * Parser *
 **********/

ograph_t *lo_graph_init()
{
	ograph_t *g;
	g = calloc(1, sizeof(ograph_t));
	g->n = kh_init(name);
	g->e = kh_init(edge);
	return g;
}

void lo_print_edge(const ograph_t *g)
{
	khint_t k;
	for (k = 0; k != kh_end(g->e); ++k) {
		if (kh_exist(g->e, k) && !kh_val(g->e, k).reduced) {
			int id[2];
			edgeinfo_t *e = &kh_val(g->e, k);
			id[0] = kh_key(g->e, k)>>32;
			id[1] = (uint32_t)kh_key(g->e, k);
			if (e->s[1] < e->e[1])
				printf("%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%.3f\n", lo_edge_label[e->type],
					   g->v.a[id[0]].name, e->l[0], e->s[0], e->e[0],
					   g->v.a[id[1]].name, e->l[1], e->s[1], e->e[1], e->usc);
			else
				printf("%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%.3f\n", lo_edge_label[e->type],
					   g->v.a[id[0]].name, e->l[0], e->e[0], e->s[0],
					   g->v.a[id[1]].name, e->l[1], e->e[1], e->s[1], e->usc);
		}
	}
}

int lo_infer_edge_type(const lo_opt_t *opt, int l[2], int s[2], int e[2], int d[2])
{
	int el, x[2], a[2], r[2]; // x: eXtended length, a: Aligned length; r: Remaining length
	int t[2][2], type;

	t[0][1] = s[1], t[1][1] = l[1] - e[1];
	if (s[0] < e[0]) t[0][0] = s[0], t[1][0] = l[0] - e[0];
	else t[0][0] = l[0] - s[0], t[1][0] = e[0];

	x[0] = a[0] = abs(e[0] - s[0]);
	x[1] = a[1] = e[1] - s[1];
	r[0] = t[0][0] - t[0][1];
	r[1] = t[1][1] - t[1][0];
	el  = r[0] < 0? t[0][0] : t[0][1];
	el += r[1] < 0? t[1][1] : t[1][0];
	x[0] += el, x[1] += el;
	d[0] = d[1] = -1;
	if ((double)a[0] / x[0] >= opt->min_aln_ratio && (double)a[1] / x[1] >= opt->min_aln_ratio) {
		if ((r[0] >= opt->min_ext && r[1] >= opt->min_ext) || (r[0] <= -opt->min_ext && r[1] <= -opt->min_ext)) { // suffix-prefix match
			type = s[0] < e[0]? 0 : 2;
			if (r[0] < 0) type ^= 3, d[0] = -r[1], d[1] = -r[0]; // reverse the direction
			else d[0] = r[0], d[1] = r[1];
		} else type = (double)x[0] / l[0] > (double)x[1] / l[1]? LO_T_C1 : LO_T_C2;
	} else type = LO_T_I; // internal local match; not a suffix-prefix match
	return type;
}

ograph_t *lo_graph_parse(const lo_opt_t *opt, kstream_t *ks)
{
	ograph_t *g;
	kstring_t str = {0,0,0};
	char *p, *q;
	khint_t k;
	int dret, absent;
	g = lo_graph_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		int i, id[2];
		edgeinfo_t e;
		for (p = q = str.s, i = 0;; ++q) {
			if (*q != '\t' && *q != 0) continue;
			if (i == 0 || i == 4) {
				int c = *q;
				*q = 0;
				k = kh_get(name, g->n, p);
				if (k == kh_end(g->n)) { // a new entry
					vertex_t *z;
					z = kv_pushp(vertex_t, g->v);
					z->id = kh_size(g->n);
					z->name = strdup(p);
					z->flag = z->state = 0;
					z->nei[0] = z->nei[1] = 0; // don't initialize the neighbor list right now
					k = kh_put(name, g->n, z->name, &absent);
					assert(absent);
					kh_val(g->n, k) = z->id;
				}
				id[(i==4)] = kh_val(g->n, k);
				*q = c;
			}
			else if (i == 1) e.l[0] = strtol(p, &p, 10);
			else if (i == 2) e.s[0] = strtol(p, &p, 10);
			else if (i == 3) e.e[0] = strtol(p, &p, 10);
			else if (i == 5) e.l[1] = strtol(p, &p, 10);
			else if (i == 6) e.s[1] = strtol(p, &p, 10);
			else if (i == 7) e.e[1] = strtol(p, &p, 10);
			else if (i == 8) e.usc  = strtod(p, &p);
			++i;
			p = q + 1;
			if (*q == 0) break;
		}
		if (i < 9) continue; // not enough fields
		e.type = lo_infer_edge_type(opt, e.l, e.s, e.e, e.d);
		if (e.type == LO_T_C1) {
			g->v.a[id[0]].flag |= LO_VF_CONTAINED;
		} else if (e.type == LO_T_C2) {
			g->v.a[id[1]].flag |= LO_VF_CONTAINED;
		} else if (e.type < 4) { // a suffix-prefix overlap
			uint64_t x = (uint64_t)id[0]<<32 | id[1];
			int sc_new, sc_old;
			edgeinfo_t *f;
			k = kh_put(edge, g->e, x, &absent);
			f = &kh_val(g->e, k);
			sc_old = f->usc * (abs(f->s[0] - f->e[0]) > abs(f->s[1] - f->e[1])? abs(f->s[0] - f->e[0]) : abs(f->s[1] - f->e[1]));
			sc_new = e.usc * (abs(e.s[0] - e.e[0]) > abs(e.s[1] - e.e[1])? abs(e.s[0] - e.e[0]) : abs(e.s[1] - e.e[1]));
			if (absent || sc_old < sc_new) // TODO: compare the total score, not unit score!
				kh_val(g->e, k) = e;
		}
	}
	free(str.s);
	if (lo_verbose >= 3)
		fprintf(stderr, "[M::%s] read %d edges\n", __func__, kh_size(g->e));
	return g;
}

void lo_graph_destroy(ograph_t *g)
{
	int i;
	for (i = 0; i < g->v.n; ++i) {
		if (g->v.a[i].nei[0]) free(g->v.a[i].nei[0]->a);
		if (g->v.a[i].nei[1]) free(g->v.a[i].nei[1]->a);
		free(g->v.a[i].name);
	}
	free(g->v.a);
	kh_destroy(edge, g->e);
	kh_destroy(name, g->n);
	free(g);
}

/******************
 * Graph routines *
 ******************/

#define lo_swap(tmp, a, b) ((tmp) = (a), (a) = (b), (b) = (tmp))

static inline void lo_flip_edge(edgeinfo_t *e)
{
	int tmp;
	lo_swap(tmp, e->l[0], e->l[1]);
	lo_swap(tmp, e->s[0], e->s[1]);
	lo_swap(tmp, e->e[0], e->e[1]);
	lo_swap(tmp, e->d[0], e->d[1]);
	e->type = ((e->type&1)<<1 | (e->type&2)>>1) ^ 3;
}

void lo_add_missing(ograph_t *g)
{
	khint_t k, l;
	int absent;
	ehash_t *added;
	added = kh_init(edge);
	for (k = 0; k != kh_end(g->e); ++k) {
		int id[2];
		if (!kh_exist(g->e, k)) continue;
		id[0] = kh_key(g->e, k)>>32;
		id[1] = (uint32_t)kh_key(g->e, k);
		if (!lo_skipped(&g->v.a[id[0]]) && !lo_skipped(&g->v.a[id[1]])) {
			uint64_t key2 = (uint64_t)id[1]<<32|id[0];
			l = kh_get(edge, g->e, key2);
			if (l == kh_end(g->e)) {
				l = kh_put(edge, added, key2, &absent);
				kh_val(added, l) = kh_val(g->e, k);
				lo_flip_edge(&kh_val(added, l));
			}
		}
	}
	for (k = 0; k != kh_end(added); ++k) {
		if (!kh_exist(added, k)) continue;
		l = kh_put(edge, g->e, kh_key(added, k), &absent);
		assert(absent);
		kh_val(g->e, l) = kh_val(added, k);
	}
	if (lo_verbose >= 3)
		fprintf(stderr, "[M::%s] added %d missing edges. %d edges remain in total.\n", __func__, kh_size(added), kh_size(g->e));
	kh_destroy(edge, added);
}

void lo_rm_skipped(ograph_t *g)
{
	khint_t k;
	int n_del = 0;
	for (k = 0; k != kh_end(g->e); ++k) {
		int id[2];
		if (!kh_exist(g->e, k)) continue;
		id[0] = kh_key(g->e, k)>>32;
		id[1] = (uint32_t)kh_key(g->e, k);
		if (lo_skipped(&g->v.a[id[0]]) || lo_skipped(&g->v.a[id[1]])) {
			++n_del;
			kh_del(edge, g->e, k); // kh_del() will not trigger rehash
		}
	}
	if (lo_verbose >= 3)
		fprintf(stderr, "[M::%s] removed %d edges; %d remain\n", __func__, n_del, kh_size(g->e));
}

void lo_rm_conflict(ograph_t *g)
{
}

void lo_mark_lack_ovlp(ograph_t *g, int min_n_ovlp) // can only be called before lo_populate_nei()
{
	int *count, i, n_marked = 0;
	khint_t k;
	count = calloc(g->v.n<<1, sizeof(int));
	for (k = 0; k != kh_end(g->e); ++k) {
		int id[2];
		edgeinfo_t *e;
		if (!kh_exist(g->e, k) || kh_val(g->e, k).type >= 4) continue;
		e = &kh_val(g->e, k);
		id[0] = kh_key(g->e, k)>>32;
		id[1] = (uint32_t)kh_key(g->e, k);
		++count[id[0]<<1|(e->type>>1^1)];
		++count[id[1]|(e->type&1)];
	}
	for (i = 0; i < g->v.n; ++i)
		if (!lo_skipped(&g->v.a[i]) && (count[i<<1|0] < min_n_ovlp || count[i<<1|1] < min_n_ovlp))
			g->v.a[i].flag |= LO_VF_LACK_OVLP, ++n_marked;
	free(count);
	if (lo_verbose >= 3) fprintf(stderr, "[M::%s] %d vertices to be dropped due to lack of overlaps\n", __func__, n_marked);
}

void lo_print_nei(ograph_t *g)
{
	int i;
	for (i = 0; i < g->v.n; ++i) {
		vertex_t *p = &g->v.a[i];
		if (!lo_skipped(p) && p->nei[0]) {
			int j, k;
			printf("%s\t%ld,%ld", p->name, p->nei[0]->n, p->nei[1]->n);
			for (j = 0; j < 2; ++j) {
				if (p->nei[j]->n) {
					putchar('\t');
					for (k = 0; k < p->nei[j]->n; ++k) {
						lo_nei_t *q = &p->nei[j]->a[k];
						if (k) putchar(',');
						printf("%s%c%d:%c", g->v.a[q->id].name, "<>"[q->ori], q->dist, "+-"[q->reduced]);
					}
				} else printf("\t*");
			}
			putchar('\n');
		}
	}
}

void lo_populate_nei(ograph_t *g)
{
	int i;
	khint_t k;
	for (i = 0; i < g->v.n; ++i) {
		if (lo_skipped(&g->v.a[i])) continue;
		g->v.a[i].nei[0] = calloc(1, sizeof(lo_nei_v));
		g->v.a[i].nei[1] = calloc(1, sizeof(lo_nei_v));
	}
	for (k = 0; k != kh_end(g->e); ++k) {
		int id[2];
		edgeinfo_t *e;
		lo_nei_t *p;
		if (!kh_exist(g->e, k)) continue;
		id[0] = kh_key(g->e, k)>>32;
		id[1] = (uint32_t)kh_key(g->e, k);
		if (id[0] > id[1]) continue;
		e = &kh_val(g->e, k);
		p = kv_pushp(lo_nei_t, *g->v.a[id[0]].nei[e->type>>1^1]);
		p->dist = e->d[0], p->id = id[1], p->ori = e->type&1, p->reduced = 0;
		p = kv_pushp(lo_nei_t, *g->v.a[id[1]].nei[e->type&1]);
		p->dist = e->d[1], p->id = id[0], p->ori = e->type>>1^1, p->reduced = 0;
	}
	for (i = 0; i < g->v.n; ++i) {
		vertex_t *p = &g->v.a[i];
		if (p->nei[0]) ks_introsort(nei, p->nei[0]->n, p->nei[0]->a);
		if (p->nei[1]) ks_introsort(nei, p->nei[1]->n, p->nei[1]->a);
	}
}

static inline edgeinfo_t *lo_get_edge(ehash_t *e, int id0, int id1)
{
	khint_t k;
	k = kh_get(edge, e, (uint64_t)id0<<32 | id1);
	return k == kh_end(e)? 0 : &kh_val(e, k);
}

void lo_trans_reduce(ograph_t *g, int fd) // fd: fuzzy distance
{
	int i, j, k, l;
	for (i = 0; i < g->v.n; ++i)
		g->v.a[i].state = 0;
	for (i = 0; i < g->v.n; ++i) {
		vertex_t *pi = &g->v.a[i];
		if (pi->nei[0] == 0) continue;
		for (j = 0; j < 2; ++j) {
			int max;
			lo_nei_v *q = pi->nei[j];
			if (q->n == 0) continue;
			for (k = 0; k < q->n; ++k)
				g->v.a[q->a[k].id].state = 1;
			max = q->a[q->n - 1].dist + fd;
			// loop between line 9--14
			for (k = 0; k < q->n; ++k) {
				vertex_t *pk = &g->v.a[q->a[k].id];
				if (pk->state == 1) {
					lo_nei_v *r = pk->nei[q->a[k].ori^1];
					for (l = 0; l < r->n; ++l)
						if (r->a[l].dist + q->a[k].dist < max && g->v.a[r->a[l].id].state == 1)
							g->v.a[r->a[l].id].state = 2;
				}
			}
			// loop between line 20--23
			for (k = 0; k < q->n; ++k) {
				if (g->v.a[q->a[k].id].state == 2) {
					edgeinfo_t *e;
					e = lo_get_edge(g->e, i, q->a[k].id);
					e->reduced = 1;
					e = lo_get_edge(g->e, q->a[k].id, i);
					e->reduced = 1;
				}
				g->v.a[q->a[k].id].state = 0;
			}
		}
	}
	for (i = 0; i < g->v.n; ++i) {
		vertex_t *pi = &g->v.a[i];
		if (pi->nei[0] == 0) continue;
		for (j = 0; j < 2; ++j) {
			lo_nei_v *q = pi->nei[j];
			for (k = 0; k < q->n; ++k)
				q->a[k].reduced = lo_get_edge(g->e, i, q->a[k].id)->reduced;
		}
	}
	if (lo_verbose == 4) lo_print_edge(g);
	if (lo_verbose == 5) lo_print_nei(g);
}

/*****************
 * Main function *
 *****************/

int main_layout(int argc, char *argv[])
{
	gzFile fp;
	kstream_t *ks;
	lo_opt_t opt;
	ograph_t *g;
	int c;

	lo_opt_init(&opt);
	while ((c = getopt(argc, argv, "v:d:o:")) >= 0) {
		if (c == 'v') lo_verbose = atoi(optarg);
		else if (c == 'd') opt.fuzzy_dist = atoi(optarg);
		else if (c == 'o') opt.min_n_ovlp = atoi(optarg);
	}
	if (argc == optind && isatty(fileno(stdin))) {
		fprintf(stderr, "Usage: bwa layout <in.ovlp>\n");
		return 1;
	}
	fp = (optind == argc && !isatty(fileno(stdin))) || strcmp(argv[optind], "-") == 0? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	ks = ks_init(fp);
	g = lo_graph_parse(&opt, ks);
	lo_rm_skipped(g);
	lo_add_missing(g);
	lo_rm_conflict(g);
	lo_mark_lack_ovlp(g, opt.min_n_ovlp);
	lo_rm_skipped(g);
	lo_populate_nei(g);
	if (lo_verbose == 6) lo_print_edge(g);
	lo_trans_reduce(g, opt.fuzzy_dist);
	lo_graph_destroy(g);
	ks_destroy(ks);
	gzclose(fp);
	return 0;
}
