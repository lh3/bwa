#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
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
	int min_ext;
	float min_aln_ratio;
} lo_opt_t;

KHASH_SET_INIT_INT(int)
typedef khash_t(int) inthash_t;

typedef struct {
	int id;
	int contained;
	char *name;
	inthash_t *nei;
} vertex_t;

typedef kvec_t(vertex_t) vertex_v;

typedef struct {
	int type, score, l[2], s[2], e[2]; // length, start and end
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
}

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

int lo_infer_edge_type(const lo_opt_t *opt, int l[2], int s[2], int e[2])
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
	if ((float)a[0] / x[0] >= opt->min_aln_ratio && (float)a[1] / x[1] >= opt->min_aln_ratio) {
		if ((r[0] >= opt->min_ext && r[1] >= opt->min_ext) || (r[0] <= -opt->min_ext && r[1] <= -opt->min_ext)) { // suffix-prefix match
			type = s[0] < e[0]? 0 : 2;
			if (r[0] < 0) type ^= 3; // reverse the direction
		} else type = x[0] / l[0] > x[1] / l[1]? LO_T_C1 : LO_T_C2;
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
					z->contained = 0;
					z->nei = 0; // don't initialize the hash table right now
					k = kh_put(name, g->n, z->name, &absent);
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
			else if (i == 8) e.score= strtol(p, &p, 10);
			++i;
			p = q + 1;
			if (*q == 0) break;
		}
		if (i < 9) continue; // not enough fields
		e.type = lo_infer_edge_type(opt, e.l, e.s, e.e);
		if (e.type == LO_T_C1) {
			g->v.a[id[0]].contained = 1;
		} else if (e.type == LO_T_C2) {
			g->v.a[id[1]].contained = 1;
		} else if (e.type < 4) { // a suffix-prefix overlap
			uint64_t x = (uint64_t)id[0]<<32 | id[1];
			k = kh_put(edge, g->e, x, &absent);
			if (absent || kh_val(g->e, k).score < e.score) 
				kh_val(g->e, k) = e;
		}
//		printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", id[0], e.l[0], e.s[0], e.e[0], id[1], e.l[1], e.s[1], e.e[1]);
	}
	free(str.s);
	if (lo_verbose >= 3) fprintf(stderr, "[M::%s] read %d edges\n", __func__, kh_size(g->e));
	return g;
}

void lo_graph_destroy(ograph_t *g)
{
	int i;
	for (i = 0; i < g->v.n; ++i) {
		if (g->v.a[i].nei) kh_destroy(int, g->v.a[i].nei);
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

void lo_rm_contained(ograph_t *g)
{
	khint_t k;
	for (k = 0; k != kh_end(g->e); ++k) {
		int id[2];
		if (!kh_exist(g->e, k)) continue;
		id[0] = kh_key(g->e, k)>>32;
		id[1] = (uint32_t)kh_key(g->e, k);
		if (g->v.a[id[0]].contained || g->v.a[id[1]].contained)
			kh_del(edge, g->e, k);
	}
	if (lo_verbose >= 3) fprintf(stderr, "[M::%s] %d edges remain after removing contained vertices\n", __func__, kh_size(g->e));
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
	while ((c = getopt(argc, argv, "")) >= 0) {
	}
	if (argc == optind && isatty(fileno(stdin))) {
		fprintf(stderr, "Usage: bwa layout <in.ovlp>\n");
		return 1;
	}
	fp = (optind == argc && !isatty(fileno(stdin))) || strcmp(argv[optind], "-") == 0? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	ks = ks_init(fp);
	g = lo_graph_parse(&opt, ks);
	lo_rm_contained(g);
	lo_graph_destroy(g);
	ks_destroy(ks);
	gzclose(fp);
	return 0;
}
