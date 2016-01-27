#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "rle.h"
#include "rope.h"

/*******************
 *** Memory Pool ***
 *******************/

#define MP_CHUNK_SIZE 0x100000 // 1MB per chunk

typedef struct { // memory pool for fast and compact memory allocation (no free)
	int size, i, n_elems;
	int64_t top, max;
	uint8_t **mem;
} mempool_t;

static mempool_t *mp_init(int size)
{
	mempool_t *mp;
	mp = calloc(1, sizeof(mempool_t));
	mp->size = size;
	mp->i = mp->n_elems = MP_CHUNK_SIZE / size;
	mp->top = -1;
	return mp;
}

static void mp_destroy(mempool_t *mp)
{
	int64_t i;
	for (i = 0; i <= mp->top; ++i) free(mp->mem[i]);
	free(mp->mem); free(mp);
}

static inline void *mp_alloc(mempool_t *mp)
{
	if (mp->i == mp->n_elems) {
		if (++mp->top == mp->max) {
			mp->max = mp->max? mp->max<<1 : 1;
			mp->mem = realloc(mp->mem, sizeof(void*) * mp->max);
		}
		mp->mem[mp->top] = calloc(mp->n_elems, mp->size);
		mp->i = 0;
	}
	return mp->mem[mp->top] + (mp->i++) * mp->size;
}

/***************
 *** B+ rope ***
 ***************/

rope_t *rope_init(int max_nodes, int block_len)
{
	rope_t *rope;
	rope = calloc(1, sizeof(rope_t));
	if (block_len < 32) block_len = 32;
	rope->max_nodes = (max_nodes+ 1)>>1<<1;
	rope->block_len = (block_len + 7) >> 3 << 3;
	rope->node = mp_init(sizeof(rpnode_t) * rope->max_nodes);
	rope->leaf = mp_init(rope->block_len);
	rope->root = mp_alloc(rope->node);
	rope->root->n = 1;
	rope->root->is_bottom = 1;
	rope->root->p = mp_alloc(rope->leaf);
	return rope;
}

void rope_destroy(rope_t *rope)
{
	mp_destroy(rope->node);
	mp_destroy(rope->leaf);
	free(rope);
}

static inline rpnode_t *split_node(rope_t *rope, rpnode_t *u, rpnode_t *v)
{ // split $v's child. $u is the first node in the bucket. $v and $u are in the same bucket. IMPORTANT: there is always enough room in $u
	int j, i = v - u;
	rpnode_t *w; // $w is the sibling of $v
	if (u == 0) { // only happens at the root; add a new root
		u = v = mp_alloc(rope->node);
		v->n = 1; v->p = rope->root; // the new root has the old root as the only child
		memcpy(v->c, rope->c, 48);
		for (j = 0; j < 6; ++j) v->l += v->c[j];
		rope->root = v;
	}
	if (i != u->n - 1) // then make room for a new node
		memmove(v + 2, v + 1, sizeof(rpnode_t) * (u->n - i - 1));
	++u->n; w = v + 1;
	memset(w, 0, sizeof(rpnode_t));
	w->p = mp_alloc(u->is_bottom? rope->leaf : rope->node);
	if (u->is_bottom) { // we are at the bottom level; $v->p is a string instead of a node
		uint8_t *p = (uint8_t*)v->p, *q = (uint8_t*)w->p;
		rle_split(p, q);
		rle_count(q, w->c);
	} else { // $v->p is a node, not a string
		rpnode_t *p = v->p, *q = w->p; // $v and $w are siblings and thus $p and $q are cousins
		p->n -= rope->max_nodes>>1;
		memcpy(q, p + p->n, sizeof(rpnode_t) * (rope->max_nodes>>1));
		q->n = rope->max_nodes>>1; // NB: this line must below memcpy() as $q->n and $q->is_bottom are modified by memcpy()
		q->is_bottom = p->is_bottom;
		for (i = 0; i < q->n; ++i)
			for (j = 0; j < 6; ++j)
				w->c[j] += q[i].c[j];
	}
	for (j = 0; j < 6; ++j) // compute $w->l and update $v->c
		w->l += w->c[j], v->c[j] -= w->c[j];
	v->l -= w->l; // update $v->c
	return v;
}

int64_t rope_insert_run(rope_t *rope, int64_t x, int a, int64_t rl, rpcache_t *cache)
{ // insert $a after $x symbols in $rope and the returns rank(a, x)
	rpnode_t *u = 0, *v = 0, *p = rope->root; // $v is the parent of $p; $u and $v are at the same level and $u is the first node in the bucket
	int64_t y = 0, z = 0, cnt[6];
	int n_runs;
	do { // top-down update. Searching and node splitting are done together in one pass.
		if (p->n == rope->max_nodes) { // node is full; split
			v = split_node(rope, u, v); // $v points to the parent of $p; when a new root is added, $v points to the root
			if (y + v->l < x) // if $v is not long enough after the split, we need to move both $p and its parent $v
				y += v->l, z += v->c[a], ++v, p = v->p;
		}
		u = p;
		if (v && x - y > v->l>>1) { // then search backwardly for the right node to descend
			p += p->n - 1; y += v->l; z += v->c[a];
			for (; y >= x; --p) y -= p->l, z -= p->c[a];
			++p;
		} else for (; y + p->l < x; ++p) y += p->l, z += p->c[a]; // then search forwardly
		assert(p - u < u->n);
		if (v) v->c[a] += rl, v->l += rl; // we should not change p->c[a] because this may cause troubles when p's child is split
		v = p; p = p->p; // descend
	} while (!u->is_bottom);
	rope->c[a] += rl; // $rope->c should be updated after the loop as adding a new root needs the old $rope->c counts
	if (cache) {
		if (cache->p != (uint8_t*)p) memset(cache, 0, sizeof(rpcache_t));
		n_runs = rle_insert_cached((uint8_t*)p, x - y, a, rl, cnt, v->c, &cache->beg, cache->bc);
		cache->p = (uint8_t*)p;
	} else n_runs = rle_insert((uint8_t*)p, x - y, a, rl, cnt, v->c);
	z += cnt[a];
	v->c[a] += rl; v->l += rl; // this should be after rle_insert(); otherwise rle_insert() won't work
	if (n_runs + RLE_MIN_SPACE > rope->block_len) {
		split_node(rope, u, v);
		if (cache) memset(cache, 0, sizeof(rpcache_t));
	}
	return z;
}

static rpnode_t *rope_count_to_leaf(const rope_t *rope, int64_t x, int64_t cx[6], int64_t *rest)
{
	rpnode_t *u, *v = 0, *p = rope->root;
	int64_t y = 0;
	int a;

	memset(cx, 0, 48);
	do {
		u = p;
		if (v && x - y > v->l>>1) {
			p += p->n - 1; y += v->l;
			for (a = 0; a != 6; ++a) cx[a] += v->c[a];
			for (; y >= x; --p) {
				y -= p->l;
				for (a = 0; a != 6; ++a) cx[a] -= p->c[a];
			}
			++p;
		} else {
			for (; y + p->l < x; ++p) {
				y += p->l;
				for (a = 0; a != 6; ++a) cx[a] += p->c[a];
			}
		}
		v = p; p = p->p;
	} while (!u->is_bottom);
	*rest = x - y;
	return v;
}

void rope_rank2a(const rope_t *rope, int64_t x, int64_t y, int64_t *cx, int64_t *cy)
{
	rpnode_t *v;
	int64_t rest;
	v = rope_count_to_leaf(rope, x, cx, &rest);
	if (y < x || cy == 0) {
		rle_rank1a((const uint8_t*)v->p, rest, cx, v->c);
	} else if (rest + (y - x) <= v->l) {
		memcpy(cy, cx, 48);
		rle_rank2a((const uint8_t*)v->p, rest, rest + (y - x), cx, cy, v->c);
	} else {
		rle_rank1a((const uint8_t*)v->p, rest, cx, v->c);
		v = rope_count_to_leaf(rope, y, cy, &rest);
		rle_rank1a((const uint8_t*)v->p, rest, cy, v->c);
	}
}

/*********************
 *** Rope iterator ***
 *********************/

void rope_itr_first(const rope_t *rope, rpitr_t *i)
{
	memset(i, 0, sizeof(rpitr_t));
	i->rope = rope;
	for (i->pa[i->d] = rope->root; !i->pa[i->d]->is_bottom;) // descend to the leftmost leaf
		++i->d, i->pa[i->d] = i->pa[i->d - 1]->p;
}

const uint8_t *rope_itr_next_block(rpitr_t *i)
{
	const uint8_t *ret;
	assert(i->d < ROPE_MAX_DEPTH); // a B+ tree should not be that tall
	if (i->d < 0) return 0;
	ret = (uint8_t*)i->pa[i->d][i->ia[i->d]].p;
	while (i->d >= 0 && ++i->ia[i->d] == i->pa[i->d]->n) i->ia[i->d--] = 0; // backtracking
	if (i->d >= 0)
		while (!i->pa[i->d]->is_bottom) // descend to the leftmost leaf
			++i->d, i->pa[i->d] = i->pa[i->d - 1][i->ia[i->d - 1]].p;
	return ret;
}

/***********
 *** I/O ***
 ***********/

void rope_print_node(const rpnode_t *p)
{
	if (p->is_bottom) {
		int i;
		putchar('(');
		for (i = 0; i < p->n; ++i) {
			uint8_t *block = (uint8_t*)p[i].p;
			const uint8_t *q = block + 2, *end = block + 2 + *rle_nptr(block);
			if (i) putchar(',');
			while (q < end) {
				int c = 0;
				int64_t j, l;
				rle_dec1(q, c, l);
				for (j = 0; j < l; ++j) putchar("$ACGTN"[c]);
			}
		}
		putchar(')');
	} else {
		int i;
		putchar('(');
		for (i = 0; i < p->n; ++i) {
			if (i) putchar(',');
			rope_print_node(p[i].p);
		}
		putchar(')');
	}
}

void rope_dump_node(const rpnode_t *p, FILE *fp)
{
	int16_t i, n = p->n;
	uint8_t is_bottom = p->is_bottom;
	fwrite(&is_bottom, 1, 1, fp);
	fwrite(&n, 2, 1, fp);
	if (is_bottom) {
		for (i = 0; i < n; ++i) {
			fwrite(p[i].c, 8, 6, fp);
			fwrite(p[i].p, 1, *rle_nptr(p[i].p) + 2, fp);
		}
	} else {
		for (i = 0; i < p->n; ++i)
			rope_dump_node(p[i].p, fp);
	}
}

void rope_dump(const rope_t *r, FILE *fp)
{
	fwrite(&r->max_nodes, 4, 1, fp);
	fwrite(&r->block_len, 4, 1, fp);
	rope_dump_node(r->root, fp);
}

rpnode_t *rope_restore_node(const rope_t *r, FILE *fp, int64_t c[6])
{
	uint8_t is_bottom, a;
	int16_t i, n;
	rpnode_t *p;
	fread(&is_bottom, 1, 1, fp);
	fread(&n, 2, 1, fp);
	p = mp_alloc(r->node);
	p->is_bottom = is_bottom, p->n = n;
	if (is_bottom) {
		for (i = 0; i < n; ++i) {
			uint16_t *q;
			p[i].p = mp_alloc(r->leaf);
			q = rle_nptr(p[i].p);
			fread(p[i].c, 8, 6, fp);
			fread(q, 2, 1, fp);
			fread(q + 1, 1, *q, fp);
		}
	} else {
		for (i = 0; i < n; ++i)
			p[i].p = rope_restore_node(r, fp, p[i].c);
	}
	memset(c, 0, 48);
	for (i = 0; i < n; ++i) {
		p[i].l = 0;
		for (a = 0; a < 6; ++a)
			c[a] += p[i].c[a], p[i].l += p[i].c[a];
	}
	return p;
}

rope_t *rope_restore(FILE *fp)
{
	rope_t *r;
	r = calloc(1, sizeof(rope_t));
	fread(&r->max_nodes, 4, 1, fp);
	fread(&r->block_len, 4, 1, fp);
	r->node = mp_init(sizeof(rpnode_t) * r->max_nodes);
	r->leaf = mp_init(r->block_len);
	r->root = rope_restore_node(r, fp, r->c);
	return r;
}
