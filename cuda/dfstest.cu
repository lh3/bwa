/* Device-DFS bit-exactness validation harness -- part of the bwa `gpualn` GPU
   BWA-backtrack port for ancient DNA.
   Copyright (C) 2026  teepean  <https://github.com/teepean>
   Derived from bwa (Heng Li; Broad Institute / Dana-Farber / Genome Research Ltd.).

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>. */

/* Phase 2 step 1: bit-exact device DFS (hit-detection) for bwa aln.
 *
 * The device DFS reproduces bwt_match_gap's bounded search tree and reports, per read,
 * whether ANY hit exists within the initial diff bound (has_hit). Rationale (see
 * ../CUDA_PORT_PLAN.md sec 6a): for reads that never reach a hit -- and when the
 * max_entries cap is never hit -- best-first and DFS enumerate the identical node set,
 * so existence-of-hit on the GPU == (n_aln>0) on the CPU. The ~0.5% reads with hits are
 * reconciled by the exact CPU bwt_match_gap to yield a bit-exact .sai.
 *
 * Host side reuses bwa's real read I/O + preprocessing (bwt_cal_width, complement, per-read
 * max_diff) so the GPU and CPU see identical (seq, width, opt).
 *
 * Build: make dfstest
 * Run:   ./dfstest <ref.fa> <reads.fq> [golden.sai]
 *        replicates `bwa aln -l 1024 -n 0.01 -o 2 -t 16`.
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <ctime>
#include <thread>
#include <chrono>
#include <cuda_runtime.h>

extern "C" {
#include "bwt.h"
#include "bwtaln.h"
#include "bwtgap.h"
int bwt_cal_width(const bwt_t *bwt, int len, const ubyte_t *str, bwt_width_t *width); /* not in a header */
}

#define FM_DEVICE_DEFINE_CONST
#include "fm_device.cuh"

#define CK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { \
	fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e_)); \
	exit(1); } } while (0)

#define STATE_M 0
#define STATE_I 1
#define STATE_D 2

/* per-read parameters that vary (the rest are kernel scalars, identical for all reads) */
struct ReadParam { uint32_t seq_off, w_off; int len, max_diff; };

/* DFS node packed into a u32: i(0-8) | n_mm(9-14) | n_gapo(15-18) | n_gape(19-22) | state(23-24) */
__device__ __forceinline__ uint32_t pack_node(int i, int mm, int go, int ge, int st)
{ return (uint32_t)i | (mm<<9) | (go<<15) | (ge<<19) | (st<<23); }

/* ---- device DFS over a COALESCED SoA stack: returns 1 if any hit exists within bound.
 * Stack entry d for this thread lives at [(size_t)d*P + slot] in each of sk/sl/sn, so a
 * warp's 32 threads touch 32 contiguous slots per push/pop (coalesced). ---- */
__device__ int d_dfs_has_hit(const fmidx_dev fm, const uint8_t *seq, int len,
                             const uint64_t *w_w, const int *w_bid, int max_diff,
                             int max_gapo, int max_gape, int mode, int indel_end_skip,
                             int max_del_occ, uint64_t *sk, uint64_t *sl, uint32_t *sn,
                             int P, int slot, int cap, int *overflow,
                             unsigned long long budget, int *nflag, unsigned long long &nn, int &maxsp)
{
	nn = 0; maxsp = 0;
	int nN = 0;
	for (int j = 0; j < len; ++j) if (seq[j] > 3) ++nN;
	if (nN > max_diff) return 0;

	#define SIDX(d) ((size_t)(d) * P + slot)
	int sp = 0;
	/* current node held in registers: only SIBLINGS are pushed to the global stack, and
	 * the stack is read only on backtrack. This removes the push->immediately-pop global
	 * round-trip that dominated runtime (stack working set was blowing L2). */
	uint64_t k = 0, l = fm.seq_len;
	int i = len, e_mm = 0, e_go = 0, e_ge = 0, e_st = STATE_M;
	bool have_cur = true;

	for (;;) {
		if (!have_cur) {                       /* backtrack: pop a sibling */
			if (sp == 0) break;
			--sp;
			k = sk[SIDX(sp)]; l = sl[SIDX(sp)];
			uint32_t nd = sn[SIDX(sp)];
			i = nd & 0x1ff; e_mm = (nd>>9)&0x3f; e_go = (nd>>15)&0xf; e_ge = (nd>>19)&0xf; e_st = (nd>>23)&0x3;
		}
		have_cur = false;                      /* consume current; re-set only if we descend */
		++nn;
		if (nn > budget) { atomicAdd(nflag, 1); return 1; }

		int m = max_diff - (e_mm + e_go);
		if (mode & BWA_MODE_GAPE) m -= e_ge;
		if (m < 0) continue;
		if (i > 0 && m < w_bid[i-1]) continue;

		int hit_found = 0;
		if (i == 0) hit_found = 1;
		else if (m == 0 && (e_st == STATE_M || (mode & BWA_MODE_GAPE) || e_ge == max_gape)) {
			uint64_t kk = k, ll = l;
			if (d_bwt_match_exact_alt(fm, seq, i, &kk, &ll)) hit_found = 1;
			else continue;
		}
		if (hit_found) return 1; /* existence is all we need */

		int xi = i - 1;
		uint64_t cntk[4], cntl[4];
		d_bwt_2occ4(fm.bwt, fm.primary, k - 1, l, cntk, cntl);
		uint64_t occ = l - k + 1;
		int allow_diff = 1, allow_M = 1;
		if (xi > 0) {
			if (w_bid[xi-1] > m - 1) allow_diff = 0;
			else if (w_bid[xi-1] == m-1 && w_bid[xi] == m-1 && w_w[xi-1] == w_w[xi]) allow_M = 0;
		}
		int tmp = (mode & BWA_MODE_LOGGAP) ? d_int_log2(e_ge + e_go)/2 + 1 : e_go + e_ge;

		/* one child stays in registers (the last EMITted, descended into); EMIT pushes
		 * the previously-staged child so exactly nchild-1 hit the global stack. */
		bool have_next = false;
		uint64_t nk_=0, nl_=0; int ni_=0, nmm_=0, ngo_=0, nge_=0, nst_=0;
		#define EMIT(_i,_k,_l,_mm,_go,_ge,_st) do { \
			if (have_next) { if (sp >= cap) { *overflow = 1; return 1; } \
				sk[SIDX(sp)]=nk_; sl[SIDX(sp)]=nl_; sn[SIDX(sp)]=pack_node(ni_,nmm_,ngo_,nge_,nst_); ++sp; \
				if (sp > maxsp) maxsp = sp; } \
			nk_=(_k); nl_=(_l); ni_=(_i); nmm_=(_mm); ngo_=(_go); nge_=(_ge); nst_=(_st); have_next=true; \
		} while (0)

		if (allow_diff && xi >= indel_end_skip + tmp && len - xi >= indel_end_skip + tmp) {
			if (e_st == STATE_M) {
				if (e_go < max_gapo) {
					EMIT(xi, k, l, e_mm, e_go + 1, e_ge, STATE_I);          /* insertion */
					for (int j = 0; j < 4; ++j) {                            /* deletion */
						uint64_t dk = c_L2[j] + cntk[j] + 1, dl = c_L2[j] + cntl[j];
						if (dk <= dl) EMIT(xi + 1, dk, dl, e_mm, e_go + 1, e_ge, STATE_D);
					}
				}
			} else if (e_st == STATE_I) {
				if (e_ge < max_gape) EMIT(xi, k, l, e_mm, e_go, e_ge + 1, STATE_I);
			} else { /* STATE_D */
				if (e_ge < max_gape && (e_ge + e_go < max_diff || occ < (uint64_t)max_del_occ)) {
					for (int j = 0; j < 4; ++j) {
						uint64_t dk = c_L2[j] + cntk[j] + 1, dl = c_L2[j] + cntl[j];
						if (dk <= dl) EMIT(xi + 1, dk, dl, e_mm, e_go, e_ge + 1, STATE_D);
					}
				}
			}
		}
		if (allow_diff && allow_M) {
			for (int j = 1; j <= 4; ++j) {
				int c = (seq[xi] + j) & 3;
				int is_mm = (j != 4 || seq[xi] > 3);
				uint64_t mk = c_L2[c] + cntk[c] + 1, ml = c_L2[c] + cntl[c];
				if (mk <= ml) EMIT(xi, mk, ml, e_mm + is_mm, e_go, e_ge, STATE_M);
			}
		} else if (seq[xi] < 4) {
			int c = seq[xi] & 3;
			uint64_t mk = c_L2[c] + cntk[c] + 1, ml = c_L2[c] + cntl[c];
			if (mk <= ml) EMIT(xi, mk, ml, e_mm, e_go, e_ge, STATE_M);
		}
		#undef EMIT

		if (have_next) {                       /* descend into last-staged child (registers) */
			k = nk_; l = nl_; i = ni_; e_mm = nmm_; e_go = ngo_; e_ge = nge_; e_st = nst_;
			have_cur = true;
		}
	}
	#undef SIDX
	return 0;
}

/* ===================== Warp-cooperative engine (Idea #1 + #2) =====================
 * One read per WARP; 32 lanes co-explore the search tree as a frontier; the per-warp stack
 * lives in SHARED memory (SoA: k/l/packed), removing global stack traffic entirely.
 * Each iteration: pop up to 32 frontier nodes (top of stack, one per lane) -> EXPAND each in
 * parallel -> warp-scan COMPACT all children -> push back to the shared stack (SPILL). Returns
 * 1 on first hit (warp ballot) or on budget/overflow (-> CPU reconcile); 0 if exhausted.
 * All 32 lanes stay converged across the collective ops; sp/nn are warp-uniform. */
__device__ int d_dfs_has_hit_warp(const fmidx_dev fm, const uint8_t *seq, int len,
                                  const uint64_t *w_w, const int *w_bid, int max_diff,
                                  int max_gapo, int max_gape, int mode, int indel_end_skip,
                                  int max_del_occ, uint64_t *ks, uint64_t *ls, uint32_t *ns,
                                  int CAP, unsigned long long budget, int *flagged,
                                  unsigned long long *nn_out)
{
	const unsigned FULL = 0xffffffffu;
	int lane = threadIdx.x & 31;

	int nN = 0;
	for (int j = 0; j < len; ++j) if (seq[j] > 3) ++nN;
	if (nN > max_diff) { if (lane==0) *nn_out = 0; return 0; }

	if (lane == 0) { ks[0]=0; ls[0]=fm.seq_len; ns[0]=pack_node(len,0,0,0,STATE_M); }
	int sp = 1;                       /* warp-uniform stack pointer */
	unsigned long long nn = 0;
	int have_carry = 0;               /* lane-0 register-continue node (drain mode) */
	uint64_t carry_k = 0, carry_l = 0; uint32_t carry_nd = 0;
	__syncwarp();

	for (;;) {
		if (!have_carry && sp == 0) { if (lane==0) *nn_out = nn; return 0; }
		/* WAVE MODE: pop a wide batch (up to 32) while there is room for every popped node's
		 * <=9 children (room-bounded -> wave never overflows). DRAIN MODE: near capacity, pop
		 * 1 and keep one child in registers (carry) so single-child chains never touch the
		 * stack, burning frontier depth without growth. Correctness-neutral (has_hit only). */
		int n_active = have_carry ? 1 : (sp < 32 ? sp : 32);
		if (!have_carry) {
			int room = (int)(((long)CAP - sp) / 9);
			if (n_active > room) n_active = room;
			if (n_active < 1) n_active = 1;   /* progress; near-full handled by carry/flag below */
		}
		bool active = lane < n_active;
		uint64_t k=0, l=0; int i=0, e_mm=0, e_go=0, e_ge=0, e_st=0;
		if (have_carry) {
			if (lane == 0) { k=carry_k; l=carry_l; uint32_t nd=carry_nd;
				i=nd&0x1ff; e_mm=(nd>>9)&0x3f; e_go=(nd>>15)&0xf; e_ge=(nd>>19)&0xf; e_st=(nd>>23)&0x3; }
			have_carry = 0;
		} else {
			int idx = sp - 1 - lane;       /* this lane's popped node (top of stack) */
			if (active) { k = ks[idx]; l = ls[idx]; uint32_t nd = ns[idx];
				i = nd & 0x1ff; e_mm=(nd>>9)&0x3f; e_go=(nd>>15)&0xf; e_ge=(nd>>19)&0xf; e_st=(nd>>23)&0x3; }
			sp -= n_active;
		}
		nn += n_active;
		if (nn > budget) { if (lane==0){ *flagged=1; *nn_out=nn; } return 1; }

		/* expand this lane's node into local child list (mirrors the serial expansion) */
		uint64_t cck[9], ccl[9]; uint32_t ccn[9]; int nc = 0; bool hit = false;
		if (active) {
			int m = max_diff - (e_mm + e_go);
			if (mode & BWA_MODE_GAPE) m -= e_ge;
			bool prune = (m < 0) || (i > 0 && m < w_bid[i-1]);
			if (!prune) {
				if (i == 0) hit = true;
				else if (m == 0 && (e_st==STATE_M || (mode&BWA_MODE_GAPE) || e_ge==max_gape)) {
					uint64_t kk=k, ll=l;
					if (d_bwt_match_exact_alt(fm, seq, i, &kk, &ll)) hit = true;
					else prune = true;
				}
				if (!hit && !prune) {
					int xi = i - 1;
					uint64_t cntk[4], cntl[4];
					d_bwt_2occ4(fm.bwt, fm.primary, k - 1, l, cntk, cntl);
					uint64_t occ = l - k + 1;
					int allow_diff = 1, allow_M = 1;
					if (xi > 0) {
						if (w_bid[xi-1] > m-1) allow_diff = 0;
						else if (w_bid[xi-1]==m-1 && w_bid[xi]==m-1 && w_w[xi-1]==w_w[xi]) allow_M = 0;
					}
					int tmp = (mode & BWA_MODE_LOGGAP) ? d_int_log2(e_ge+e_go)/2+1 : e_go+e_ge;
					#define ADD(_i,_k,_l,_mm,_go,_ge,_st) do { \
						cck[nc]=(_k); ccl[nc]=(_l); ccn[nc]=pack_node((_i),(_mm),(_go),(_ge),(_st)); ++nc; } while(0)
					if (allow_diff && xi >= indel_end_skip+tmp && len-xi >= indel_end_skip+tmp) {
						if (e_st == STATE_M) {
							if (e_go < max_gapo) {
								ADD(xi, k, l, e_mm, e_go+1, e_ge, STATE_I);
								for (int j=0;j<4;++j){ uint64_t dk=c_L2[j]+cntk[j]+1, dl=c_L2[j]+cntl[j];
									if (dk<=dl) ADD(xi+1, dk, dl, e_mm, e_go+1, e_ge, STATE_D); }
							}
						} else if (e_st == STATE_I) {
							if (e_ge < max_gape) ADD(xi, k, l, e_mm, e_go, e_ge+1, STATE_I);
						} else {
							if (e_ge < max_gape && (e_ge+e_go < max_diff || occ < (uint64_t)max_del_occ))
								for (int j=0;j<4;++j){ uint64_t dk=c_L2[j]+cntk[j]+1, dl=c_L2[j]+cntl[j];
									if (dk<=dl) ADD(xi+1, dk, dl, e_mm, e_go, e_ge+1, STATE_D); }
						}
					}
					if (allow_diff && allow_M) {
						for (int j=1;j<=4;++j){ int c=(seq[xi]+j)&3; int is_mm=(j!=4||seq[xi]>3);
							uint64_t mk=c_L2[c]+cntk[c]+1, ml=c_L2[c]+cntl[c];
							if (mk<=ml) ADD(xi, mk, ml, e_mm+is_mm, e_go, e_ge, STATE_M); }
					} else if (seq[xi] < 4) {
						int c=seq[xi]&3; uint64_t mk=c_L2[c]+cntk[c]+1, ml=c_L2[c]+cntl[c];
						if (mk<=ml) ADD(xi, mk, ml, e_mm, e_go, e_ge, STATE_M);
					}
					#undef ADD
				}
			}
		}

		if (__any_sync(FULL, hit)) { if (lane==0){ *nn_out = nn; } return 1; }

		/* warp-exclusive prefix sum of per-lane child counts -> contiguous push offsets */
		int incl = nc;
		for (int d=1; d<32; d<<=1) { int y = __shfl_up_sync(FULL, incl, d); if (lane >= d) incl += y; }
		int total = __shfl_sync(FULL, incl, 31);
		int myoff = incl - nc;
		/* DRAIN register-continue: when near-full and processing a single node (only lane 0 has
		 * children), keep its last child in registers (carry) instead of pushing+re-popping. */
		bool drain_keep = (n_active == 1) && (sp >= (CAP >> 1));
		if (drain_keep) {
			int push_n = total > 0 ? total - 1 : 0;          /* keep 1 child in carry */
			if (sp + push_n > CAP) { if (lane==0){ *flagged=1; *nn_out=nn; } return 1; }
			if (lane == 0 && nc > 0) {
				for (int j=0;j<nc-1;++j){ int d=sp+j; ks[d]=cck[j]; ls[d]=ccl[j]; ns[d]=ccn[j]; }
				carry_k=cck[nc-1]; carry_l=ccl[nc-1]; carry_nd=ccn[nc-1];
			}
			sp += push_n;
			have_carry = (total > 0) ? 1 : 0;                /* uniform: total is warp-broadcast */
		} else {
			if (sp + total > CAP) { if (lane==0){ *flagged=1; *nn_out=nn; } return 1; } /* overflow -> CPU */
			for (int j=0;j<nc;++j){ int d = sp + myoff + j; ks[d]=cck[j]; ls[d]=ccl[j]; ns[d]=ccn[j]; }
			sp += total;
		}
		__syncwarp();
	}
}

/* ===== Two-level stack engine (Idea: small shared top-window + per-warp GLOBAL backing) =====
 * Shared holds the hot TOP of the stack (CAP_SM entries); the cold BOTTOM spills to a
 * pre-allocated per-warp global region via warp-parallel COALESCED copies (128-entry chunks).
 * Invariant: global[0..gsp) = oldest entries (bottom), shared[0..sp) = newest (top); stack
 * top = shared[sp-1]. 95%+ of reads never spill (stay in shared at high occupancy); bushy
 * reads spill their deep frontier to global but stay ON THE GPU (no CPU flag). Only the 2M
 * budget (or an enormous frontier > CAP_GL) flags to CPU. has_hit-only -> bit-exact. */
__device__ int d_dfs_has_hit_warp2(const fmidx_dev fm, const uint8_t *seq, int len,
                                   const uint64_t *w_w, const int *w_bid, int max_diff,
                                   int max_gapo, int max_gape, int mode, int indel_end_skip,
                                   int max_del_occ, uint64_t *sk, uint64_t *sl, uint32_t *sn, int CAP_SM,
                                   uint64_t *gk, uint64_t *gl, uint32_t *gn, int CAP_GL,
                                   unsigned long long budget, int *flagged, unsigned long long *nn_out)
{
	const unsigned FULL = 0xffffffffu; const int CHUNK = 128;
	int lane = threadIdx.x & 31;
	int nN = 0;
	for (int j = 0; j < len; ++j) if (seq[j] > 3) ++nN;
	if (nN > max_diff) { if (lane==0) *nn_out = 0; return 0; }

	if (lane == 0) { sk[0]=0; sl[0]=fm.seq_len; sn[0]=pack_node(len,0,0,0,STATE_M); }
	int sp = 1, gsp = 0;              /* shared sp, global sp (both warp-uniform) */
	unsigned long long nn = 0;
	__syncwarp();

	for (;;) {
		if (sp == 0) {
			if (gsp == 0) { if (lane==0) *nn_out = nn; return 0; }
			int c = gsp < CHUNK ? gsp : CHUNK;               /* UNSPILL: global top -> shared bottom */
			for (int j = lane; j < c; j += 32) { sk[j]=gk[gsp-c+j]; sl[j]=gl[gsp-c+j]; sn[j]=gn[gsp-c+j]; }
			gsp -= c; sp = c; __syncwarp();
		}
		int room = CAP_SM - sp;
		if (room < 9) {                                      /* SPILL bottom CHUNK -> global */
			if (gsp + CHUNK > CAP_GL) { if (lane==0){ *flagged=1; *nn_out=nn; } return 1; }
			for (int j = lane; j < CHUNK; j += 32) { gk[gsp+j]=sk[j]; gl[gsp+j]=sl[j]; gn[gsp+j]=sn[j]; }
			__syncwarp();
			for (int j = lane; j < sp - CHUNK; j += 32) { sk[j]=sk[j+CHUNK]; sl[j]=sl[j+CHUNK]; sn[j]=sn[j+CHUNK]; }
			gsp += CHUNK; sp -= CHUNK; __syncwarp();
			room = CAP_SM - sp;
		}
		int n_active = sp < 32 ? sp : 32;
		int r9 = room / 9; if (n_active > r9) n_active = r9; if (n_active < 1) n_active = 1;
		bool active = lane < n_active;
		uint64_t k=0, l=0; int i=0, e_mm=0, e_go=0, e_ge=0, e_st=0;
		if (active) { int idx = sp-1-lane; k=sk[idx]; l=sl[idx]; uint32_t nd=sn[idx];
			i=nd&0x1ff; e_mm=(nd>>9)&0x3f; e_go=(nd>>15)&0xf; e_ge=(nd>>19)&0xf; e_st=(nd>>23)&0x3; }
		sp -= n_active; nn += n_active;
		if (nn > budget) { if (lane==0){ *flagged=1; *nn_out=nn; } return 1; }

		uint64_t cck[9], ccl[9]; uint32_t ccn[9]; int nc = 0; bool hit = false;
		if (active) {
			int m = max_diff - (e_mm + e_go);
			if (mode & BWA_MODE_GAPE) m -= e_ge;
			bool prune = (m < 0) || (i > 0 && m < w_bid[i-1]);
			if (!prune) {
				if (i == 0) hit = true;
				else if (m == 0 && (e_st==STATE_M || (mode&BWA_MODE_GAPE) || e_ge==max_gape)) {
					uint64_t kk=k, ll=l;
					if (d_bwt_match_exact_alt(fm, seq, i, &kk, &ll)) hit = true; else prune = true;
				}
				if (!hit && !prune) {
					int xi = i - 1;
					uint64_t cntk[4], cntl[4];
					d_bwt_2occ4(fm.bwt, fm.primary, k - 1, l, cntk, cntl);
					uint64_t occ = l - k + 1;
					int allow_diff = 1, allow_M = 1;
					if (xi > 0) {
						if (w_bid[xi-1] > m-1) allow_diff = 0;
						else if (w_bid[xi-1]==m-1 && w_bid[xi]==m-1 && w_w[xi-1]==w_w[xi]) allow_M = 0;
					}
					int tmp = (mode & BWA_MODE_LOGGAP) ? d_int_log2(e_ge+e_go)/2+1 : e_go+e_ge;
					#define ADD2(_i,_k,_l,_mm,_go,_ge,_st) do { \
						cck[nc]=(_k); ccl[nc]=(_l); ccn[nc]=pack_node((_i),(_mm),(_go),(_ge),(_st)); ++nc; } while(0)
					if (allow_diff && xi >= indel_end_skip+tmp && len-xi >= indel_end_skip+tmp) {
						if (e_st == STATE_M) {
							if (e_go < max_gapo) {
								ADD2(xi, k, l, e_mm, e_go+1, e_ge, STATE_I);
								for (int j=0;j<4;++j){ uint64_t dk=c_L2[j]+cntk[j]+1, dl=c_L2[j]+cntl[j];
									if (dk<=dl) ADD2(xi+1, dk, dl, e_mm, e_go+1, e_ge, STATE_D); }
							}
						} else if (e_st == STATE_I) {
							if (e_ge < max_gape) ADD2(xi, k, l, e_mm, e_go, e_ge+1, STATE_I);
						} else {
							if (e_ge < max_gape && (e_ge+e_go < max_diff || occ < (uint64_t)max_del_occ))
								for (int j=0;j<4;++j){ uint64_t dk=c_L2[j]+cntk[j]+1, dl=c_L2[j]+cntl[j];
									if (dk<=dl) ADD2(xi+1, dk, dl, e_mm, e_go, e_ge+1, STATE_D); }
						}
					}
					if (allow_diff && allow_M) {
						for (int j=1;j<=4;++j){ int c=(seq[xi]+j)&3; int is_mm=(j!=4||seq[xi]>3);
							uint64_t mk=c_L2[c]+cntk[c]+1, ml=c_L2[c]+cntl[c];
							if (mk<=ml) ADD2(xi, mk, ml, e_mm+is_mm, e_go, e_ge, STATE_M); }
					} else if (seq[xi] < 4) {
						int c=seq[xi]&3; uint64_t mk=c_L2[c]+cntk[c]+1, ml=c_L2[c]+cntl[c];
						if (mk<=ml) ADD2(xi, mk, ml, e_mm, e_go, e_ge, STATE_M);
					}
					#undef ADD2
				}
			}
		}
		if (__any_sync(FULL, hit)) { if (lane==0) *nn_out = nn; return 1; }

		int incl = nc;
		for (int d=1; d<32; d<<=1) { int y = __shfl_up_sync(FULL, incl, d); if (lane >= d) incl += y; }
		int total = __shfl_sync(FULL, incl, 31);
		int myoff = incl - nc;
		/* room-bounded: total <= room, so this push always fits the shared window */
		for (int j=0;j<nc;++j){ int d = sp + myoff + j; sk[d]=cck[j]; sl[d]=ccl[j]; sn[d]=ccn[j]; }
		sp += total;
		__syncwarp();
	}
}

__global__ void k_dfs_warp2(fmidx_dev fm, const uint8_t *seq, const uint64_t *w_w, const int *w_bid,
                            const ReadParam *rp, int nreads, int max_gapo, int max_gape, int mode,
                            int indel_end_skip, int max_del_occ, int CAP_SM, int CAP_GL,
                            uint64_t *Gk, uint64_t *Gl, uint32_t *Gn, uint8_t *has_hit,
                            int *workctr, unsigned long long *npop, unsigned long long budget,
                            int *nflag, int wpb)
{
	extern __shared__ unsigned char smem[];
	int warp_in_block = threadIdx.x >> 5, lane = threadIdx.x & 31;
	uint64_t *sk = (uint64_t*)smem + (size_t)warp_in_block * CAP_SM;
	uint64_t *sl = (uint64_t*)smem + (size_t)wpb * CAP_SM + (size_t)warp_in_block * CAP_SM;
	uint32_t *sn = (uint32_t*)((uint64_t*)smem + (size_t)wpb * CAP_SM * 2) + (size_t)warp_in_block * CAP_SM;
	int gslot = blockIdx.x * wpb + warp_in_block;        /* this warp's global backing slot */
	uint64_t *gk = Gk + (size_t)gslot * CAP_GL;
	uint64_t *gl = Gl + (size_t)gslot * CAP_GL;
	uint32_t *gn = Gn + (size_t)gslot * CAP_GL;
	for (;;) {
		int r;
		if (lane == 0) r = atomicAdd(workctr, 1);
		r = __shfl_sync(0xffffffffu, r, 0);
		if (r >= nreads) break;
		ReadParam p = rp[r];
		int flagged = 0; unsigned long long nn = 0;
		int hh = d_dfs_has_hit_warp2(fm, seq + p.seq_off, p.len, w_w + p.w_off, w_bid + p.w_off,
			p.max_diff, max_gapo, max_gape, mode, indel_end_skip, max_del_occ, sk, sl, sn, CAP_SM,
			gk, gl, gn, CAP_GL, budget, &flagged, &nn);
		if (lane == 0) { has_hit[r] = (uint8_t)hh; npop[r] = nn; if (flagged) atomicAdd(nflag, 1); }
	}
}

__global__ void k_dfs_warp(fmidx_dev fm, const uint8_t *seq, const uint64_t *w_w, const int *w_bid,
                           const ReadParam *rp, int nreads, int max_gapo, int max_gape, int mode,
                           int indel_end_skip, int max_del_occ, int CAP, uint8_t *has_hit,
                           int *workctr, unsigned long long *npop, unsigned long long budget,
                           int *nflag, int wpb)
{
	extern __shared__ unsigned char smem[];
	int warp_in_block = threadIdx.x >> 5, lane = threadIdx.x & 31;
	uint64_t *ks = (uint64_t*)smem + (size_t)warp_in_block * CAP;
	uint64_t *ls = (uint64_t*)smem + (size_t)wpb * CAP + (size_t)warp_in_block * CAP;
	uint32_t *ns = (uint32_t*)((uint64_t*)smem + (size_t)wpb * CAP * 2) + (size_t)warp_in_block * CAP;
	for (;;) {
		int r;
		if (lane == 0) r = atomicAdd(workctr, 1);
		r = __shfl_sync(0xffffffffu, r, 0);
		if (r >= nreads) break;
		ReadParam p = rp[r];
		int flagged = 0; unsigned long long nn = 0;
		int hh = d_dfs_has_hit_warp(fm, seq + p.seq_off, p.len, w_w + p.w_off, w_bid + p.w_off,
			p.max_diff, max_gapo, max_gape, mode, indel_end_skip, max_del_occ, ks, ls, ns,
			CAP, budget, &flagged, &nn);
		if (lane == 0) { has_hit[r] = (uint8_t)hh; npop[r] = nn; if (flagged) atomicAdd(nflag, 1); }
	}
}

/* persistent work-pool kernel: each thread atomically pulls the next read until the queue
 * drains, so a thread finishing a tiny tree immediately grabs more work (load balancing). */
__global__ void k_dfs(fmidx_dev fm, const uint8_t *seq, const uint64_t *w_w, const int *w_bid,
                      const ReadParam *rp, int nreads, int max_gapo, int max_gape, int mode,
                      int indel_end_skip, int max_del_occ, uint64_t *sk, uint64_t *sl, uint32_t *sn,
                      int cap, uint8_t *has_hit, int *overflow, int *workctr, unsigned long long *npop,
                      unsigned long long budget, int *nflag, int *maxsp_out)
{
	/* block-local stack: this block's slice is contiguous; within it, depth levels are
	 * blockDim apart (a warp's 32 lanes stay coalesced AND a thread's adjacent depths stay
	 * close -> L1-friendly), instead of the global stride-P layout that thrashed L2. */
	int tib = threadIdx.x, stride = blockDim.x;
	size_t bb = (size_t)blockIdx.x * (size_t)blockDim.x * cap;
	uint64_t *bsk = sk + bb, *bsl = sl + bb; uint32_t *bsn = sn + bb;
	for (;;) {
		int r = atomicAdd(workctr, 1);
		if (r >= nreads) break;
		ReadParam p = rp[r];
		unsigned long long nn = 0; int maxsp = 0;
		has_hit[r] = (uint8_t)d_dfs_has_hit(fm, seq + p.seq_off, p.len,
			w_w + p.w_off, w_bid + p.w_off, p.max_diff,
			max_gapo, max_gape, mode, indel_end_skip, max_del_occ,
			bsk, bsl, bsn, stride, tib, cap, overflow, budget, nflag, nn, maxsp);
		npop[r] = nn; maxsp_out[r] = maxsp;
	}
}

/* ---- write a .sai exactly like bwa_aln_core ---- */
static void write_sai(const char *fn, const gap_opt_t *opt, int n, int *n_aln, bwt_aln1_t **aln)
{
	FILE *f = fopen(fn, "wb");
	fwrite(SAI_MAGIC, 1, 4, f);
	fwrite(opt, sizeof(gap_opt_t), 1, f);
	for (int i = 0; i < n; ++i) {
		fwrite(&n_aln[i], 4, 1, f);
		if (n_aln[i]) fwrite(aln[i], sizeof(bwt_aln1_t), n_aln[i], f);
	}
	fclose(f);
}

static int md5cmp(const char *a, const char *b)
{
	char cmd[2048]; snprintf(cmd, sizeof cmd, "md5sum '%s' '%s'", a, b);
	fprintf(stderr, "  "); fflush(stderr); return system(cmd);
}

int main(int argc, char **argv)
{
	if (argc < 3) { fprintf(stderr, "usage: %s <ref.fa> <reads.fq> [golden.sai]\n", argv[0]); return 1; }
	const char *prefix = argv[1], *reads_fn = argv[2], *golden = argc > 3 ? argv[3] : NULL;

	/* opt for: bwa aln -l 1024 -n 0.01 -o 2 -t 16 */
	gap_opt_t *opt = gap_init_opt();
	opt->seed_len = 1024; opt->fnr = 0.01; opt->max_diff = -1; opt->max_gapo = 2; opt->n_threads = 16;

	char bwt_fn[4096]; snprintf(bwt_fn, sizeof bwt_fn, "%s.bwt", prefix);
	fprintf(stderr, "[load] %s\n", bwt_fn);
	bwt_t *bwt = bwt_restore_bwt(bwt_fn);
	if (!bwt) { fprintf(stderr, "failed to load bwt\n"); return 1; }

	/* upload index + constants */
	uint32_t *d_bwt = NULL;
	CK(cudaMalloc(&d_bwt, bwt->bwt_size * sizeof(uint32_t)));
	CK(cudaMemcpy(d_bwt, bwt->bwt, bwt->bwt_size * sizeof(uint32_t), cudaMemcpyHostToDevice));
	CK(cudaMemcpyToSymbol(c_cnt_table, bwt->cnt_table, sizeof(uint32_t) * 256));
	CK(cudaMemcpyToSymbol(c_L2, bwt->L2, sizeof(uint64_t) * 5));
	fmidx_dev fm{ d_bwt, bwt->primary, bwt->seq_len };

	/* read ALL reads in one batch (sub10k < 0x40000), exactly as the driver would */
	bwa_seqio_t *ks = bwa_seq_open(reads_fn);
	int n_seqs = 0;
	bwa_seq_t *seqs = bwa_read_seq(ks, 0x40000, &n_seqs, opt->mode, opt->trim_qual);
	bwa_seq_close(ks);
	fprintf(stderr, "[reads] %d\n", n_seqs);

	/* replicate bwa_cal_sa_reg_gap preprocessing + CPU reference */
	gap_opt_t local_opt = *opt;
	int max_len = 0; for (int i = 0; i < n_seqs; ++i) if (seqs[i].len > max_len) max_len = seqs[i].len;
	if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
	if (local_opt.max_diff < local_opt.max_gapo) local_opt.max_gapo = local_opt.max_diff;
	gap_stack_t *stack = gap_init_stack(local_opt.max_diff, local_opt.max_gapo, local_opt.max_gape, &local_opt);
	gap_opt_t base = local_opt; /* stable max_gapo/mode/...; max_diff & seed_len overridden per read */
	/* full per-read CPU reference (checks 1&3) is slow; default on only for small sets */
	bool full_cpu = getenv("DFS_FULLCPU") || n_seqs <= 20000;

	std::vector<uint8_t> seq_flat; std::vector<uint64_t> w_flat; std::vector<int> bid_flat;
	std::vector<ReadParam> rp(n_seqs);
	std::vector<int> n_aln_cpu(n_seqs);
	std::vector<bwt_aln1_t*> aln_cpu(n_seqs, NULL);
	bwt_width_t *w = NULL; int max_w = 0;

	for (int i = 0; i < n_seqs; ++i) {
		bwa_seq_t *p = &seqs[i];
		if (max_w < p->len) { max_w = p->len; w = (bwt_width_t*)realloc(w, (max_w + 1) * sizeof(bwt_width_t)); }
		memset(w, 0, (p->len + 1) * sizeof(bwt_width_t));
		bwt_cal_width(bwt, p->len, p->seq, w);
		if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
		local_opt.seed_len = opt->seed_len < p->len ? opt->seed_len : 0x7fffffff;
		for (int j = 0; j < p->len; ++j) p->seq[j] = p->seq[j] > 3 ? 4 : 3 - p->seq[j]; /* complement */

		/* record inputs for the GPU (identical seq + width + per-read max_diff) */
		rp[i].seq_off = seq_flat.size(); rp[i].w_off = w_flat.size();
		rp[i].len = p->len; rp[i].max_diff = local_opt.max_diff;
		for (int j = 0; j < p->len; ++j) seq_flat.push_back(p->seq[j]);
		for (int j = 0; j <= p->len; ++j) { w_flat.push_back(w[j].w); bid_flat.push_back(w[j].bid); }

		/* optional full CPU reference for checks 1 & 3 (len <= seed_len -> seed_w NULL) */
		if (full_cpu) {
			int na = 0;
			aln_cpu[i] = bwt_match_gap(bwt, p->len, p->seq, w, (bwt_width_t*)0, &local_opt, &na, stack);
			n_aln_cpu[i] = na;
		}
	}
	free(w);

	if (full_cpu) { /* check 1: CPU-only .sai == golden (validates the harness replicates the driver) */
		write_sai("/tmp/dfstest.cpu.sai", opt, n_seqs, n_aln_cpu.data(), aln_cpu.data());
		long n_hit_cpu = 0; for (int i = 0; i < n_seqs; ++i) if (n_aln_cpu[i] > 0) ++n_hit_cpu;
		fprintf(stderr, "[cpu] reads with hits: %ld (%.2f%%)\n", n_hit_cpu, 100.0 * n_hit_cpu / n_seqs);
		if (golden) { fprintf(stderr, "[check1] CPU-only .sai vs golden:\n"); md5cmp("/tmp/dfstest.cpu.sai", golden); }
	} else fprintf(stderr, "[cpu] full reference skipped (set DFS_FULLCPU=1 to enable checks 1&3)\n");

	/* ---- GPU DFS ---- */
	uint8_t *d_seq; uint64_t *d_ww; int *d_wbid; ReadParam *d_rp; uint8_t *d_hit; int *d_ovf;
	CK(cudaMalloc(&d_seq, seq_flat.size()));
	CK(cudaMalloc(&d_ww, w_flat.size() * 8));
	CK(cudaMalloc(&d_wbid, bid_flat.size() * 4));
	CK(cudaMalloc(&d_rp, n_seqs * sizeof(ReadParam)));
	CK(cudaMalloc(&d_hit, n_seqs));
	CK(cudaMalloc(&d_ovf, 4)); CK(cudaMemset(d_ovf, 0, 4));
	CK(cudaMemcpy(d_seq, seq_flat.data(), seq_flat.size(), cudaMemcpyHostToDevice));
	CK(cudaMemcpy(d_ww, w_flat.data(), w_flat.size() * 8, cudaMemcpyHostToDevice));
	CK(cudaMemcpy(d_wbid, bid_flat.data(), bid_flat.size() * 4, cudaMemcpyHostToDevice));
	CK(cudaMemcpy(d_rp, rp.data(), n_seqs * sizeof(ReadParam), cudaMemcpyHostToDevice));

	int numSM = 0;
	CK(cudaDeviceGetAttribute(&numSM, cudaDevAttrMultiProcessorCount, 0));
	int *d_wc; CK(cudaMalloc(&d_wc, 4)); CK(cudaMemset(d_wc, 0, 4));
	unsigned long long *d_npop; CK(cudaMalloc(&d_npop, (size_t)n_seqs * 8)); CK(cudaMemset(d_npop, 0, (size_t)n_seqs * 8));
	int *d_nflag; CK(cudaMalloc(&d_nflag, 4)); CK(cudaMemset(d_nflag, 0, 4));
	int *d_maxsp = NULL;
	unsigned long long budget = getenv("DFS_BUDGET") ? strtoull(getenv("DFS_BUDGET"), NULL, 10) : 2000000ULL;
	const char *engine = getenv("DFS_ENGINE") ? getenv("DFS_ENGINE") : "thread";
	fprintf(stderr, "[gpu] engine=%s  per-read work budget = %llu pops (overflow -> CPU reconcile)\n", engine, budget);

	cudaEvent_t t0, t1; CK(cudaEventCreate(&t0)); CK(cudaEventCreate(&t1));
	CK(cudaEventRecord(t0));
	if (strcmp(engine, "warp") == 0) {
		int wpb = getenv("DFS_WARP_WPB") ? atoi(getenv("DFS_WARP_WPB")) : 4;
		int CAPW = getenv("DFS_WARP_CAP") ? atoi(getenv("DFS_WARP_CAP")) : 768;
		int bdim = wpb * 32;
		size_t shbytes = (size_t)wpb * CAPW * (8 + 8 + 4);
		CK(cudaFuncSetAttribute(k_dfs_warp, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)shbytes));
		int mb = 0; CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&mb, k_dfs_warp, bdim, shbytes));
		int nb = mb > 0 ? mb * numSM : numSM;
		fprintf(stderr, "[gpu] WARP: %d warps/blk CAP=%d shared=%.1f KB/blk -> %d blk/SM (%d warps/SM, %d reads in flight)\n",
			wpb, CAPW, shbytes/1024.0, mb, mb*wpb, nb*wpb);
		k_dfs_warp<<<nb, bdim, shbytes>>>(fm, d_seq, d_ww, d_wbid, d_rp, n_seqs,
			local_opt.max_gapo, local_opt.max_gape, local_opt.mode, local_opt.indel_end_skip,
			local_opt.max_del_occ, CAPW, d_hit, d_wc, d_npop, budget, d_nflag, wpb);
	} else if (strcmp(engine, "warp2") == 0) {
		int wpb = getenv("DFS_WARP_WPB") ? atoi(getenv("DFS_WARP_WPB")) : 4;
		int CAP_SM = getenv("DFS_WARP_CAP") ? atoi(getenv("DFS_WARP_CAP")) : 512;
		int CAP_GL = getenv("DFS_WARP_GCAP") ? atoi(getenv("DFS_WARP_GCAP")) : 16384;
		int bdim = wpb * 32;
		size_t shbytes = (size_t)wpb * CAP_SM * (8 + 8 + 4);
		CK(cudaFuncSetAttribute(k_dfs_warp2, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)shbytes));
		int mb = 0; CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&mb, k_dfs_warp2, bdim, shbytes));
		int nb = mb > 0 ? mb * numSM : numSM;
		size_t nwarps = (size_t)nb * wpb;
		uint64_t *Gk, *Gl; uint32_t *Gn;
		CK(cudaMalloc(&Gk, nwarps * CAP_GL * 8));
		CK(cudaMalloc(&Gl, nwarps * CAP_GL * 8));
		CK(cudaMalloc(&Gn, nwarps * CAP_GL * 4));
		fprintf(stderr, "[gpu] WARP2: %d w/blk CAP_SM=%d shared=%.1fKB/blk -> %d blk/SM (%d w/SM, %zu reads in flight); "
			"global backing %d/warp -> %.0f MB\n", wpb, CAP_SM, shbytes/1024.0, mb, mb*wpb, nwarps,
			CAP_GL, nwarps*CAP_GL*20.0/1e6);
		k_dfs_warp2<<<nb, bdim, shbytes>>>(fm, d_seq, d_ww, d_wbid, d_rp, n_seqs,
			local_opt.max_gapo, local_opt.max_gape, local_opt.mode, local_opt.indel_end_skip,
			local_opt.max_del_occ, CAP_SM, CAP_GL, Gk, Gl, Gn, d_hit, d_wc, d_npop, budget, d_nflag, wpb);
	} else {
		int blockDim = 128, cap = 512; /* DFS depth max ~386; overflow -> CPU reconcile (safe) */
		int maxBlocks = 0;
		CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&maxBlocks, k_dfs, blockDim, 0));
		int nblocks = maxBlocks * numSM, P = nblocks * blockDim;
		fprintf(stderr, "[gpu] THREAD: %d SM x %d blk/SM x %d thr = %d slots; cap %d -> pool %.0f MB\n",
			numSM, maxBlocks, blockDim, P, cap, (size_t)P*cap*20/1e6);
		uint64_t *d_sk, *d_sl; uint32_t *d_sn;
		CK(cudaMalloc(&d_sk, (size_t)P * cap * 8));
		CK(cudaMalloc(&d_sl, (size_t)P * cap * 8));
		CK(cudaMalloc(&d_sn, (size_t)P * cap * 4));
		CK(cudaMalloc(&d_maxsp, (size_t)n_seqs * 4)); CK(cudaMemset(d_maxsp, 0, (size_t)n_seqs * 4));
		k_dfs<<<nblocks, blockDim>>>(fm, d_seq, d_ww, d_wbid, d_rp, n_seqs,
			local_opt.max_gapo, local_opt.max_gape, local_opt.mode, local_opt.indel_end_skip,
			local_opt.max_del_occ, d_sk, d_sl, d_sn, cap, d_hit, d_ovf, d_wc, d_npop, budget, d_nflag, d_maxsp);
	}
	CK(cudaEventRecord(t1)); CK(cudaEventSynchronize(t1)); CK(cudaGetLastError());
	float ms = 0; CK(cudaEventElapsedTime(&ms, t0, t1));
	{ int nflag = 0; CK(cudaMemcpy(&nflag, d_nflag, 4, cudaMemcpyDeviceToHost));
	  fprintf(stderr, "[gpu] budget-flagged reads (sent to CPU) = %d (%.3f%%)\n", nflag, 100.0*nflag/n_seqs); }

	{ std::vector<unsigned long long> np(n_seqs); CK(cudaMemcpy(np.data(), d_npop, (size_t)n_seqs*8, cudaMemcpyDeviceToHost));
	  unsigned long long tot=0,mx=0; int amx=0; for(int i=0;i<n_seqs;i++){ tot+=np[i]; if(np[i]>mx){mx=np[i];amx=i;} }
	  fprintf(stderr, "[gpu] node-pops total=%llu mean=%.0f max=%llu (read %d, len %d)\n",
	          tot, (double)tot/n_seqs, mx, amx, rp[amx].len);
	  /* top-5 heaviest reads */
	  std::vector<int> idx(n_seqs); for(int i=0;i<n_seqs;i++) idx[i]=i;
	  std::sort(idx.begin(), idx.end(), [&](int a,int b){return np[a]>np[b];});
	  fprintf(stderr, "[gpu] heaviest reads (pops): "); for(int t=0;t<5&&t<n_seqs;t++) fprintf(stderr,"%llu ", np[idx[t]]); fprintf(stderr,"\n"); }

	if (d_maxsp) { std::vector<int> ms(n_seqs); CK(cudaMemcpy(ms.data(), d_maxsp, (size_t)n_seqs*4, cudaMemcpyDeviceToHost));
	  int mx=0; double mean=0; std::vector<int> b(13,0); /* log2 bins */
	  for(int i=0;i<n_seqs;i++){ int v=ms[i]; if(v>mx)mx=v; mean+=v; int bb=0; while((1<<bb)<v && bb<12)++bb; b[bb]++; }
	  fprintf(stderr, "[gpu] DFS stack DEPTH: mean=%.1f max=%d  (sizes per-warp shared stack)\n", mean/n_seqs, mx);
	  fprintf(stderr, "[gpu] depth histogram <=2^b: "); for(int bb=0;bb<13;bb++) if(b[bb]) fprintf(stderr,"[%d]=%d ", 1<<bb, b[bb]); fprintf(stderr,"\n"); }

	std::vector<uint8_t> has_hit(n_seqs); int ovf = 0;
	CK(cudaMemcpy(has_hit.data(), d_hit, n_seqs, cudaMemcpyDeviceToHost));
	CK(cudaMemcpy(&ovf, d_ovf, 4, cudaMemcpyDeviceToHost));
	fprintf(stderr, "[gpu] DFS time %.1f ms  (%.0f reads/s)  stack_overflow=%d\n",
		ms, n_seqs / (ms / 1e3), ovf);

	long gpu_hit = 0; for (int i = 0; i < n_seqs; ++i) if (has_hit[i]) ++gpu_hit;
	long fn = 0;
	if (full_cpu) { /* check 3: has_hit vs CPU n_aln>0 */
		long fp = 0;
		for (int i = 0; i < n_seqs; ++i) {
			int cpu = n_aln_cpu[i] > 0;
			if (cpu && !has_hit[i]) { if (fn < 10) fprintf(stderr, "  FALSE-NEG read %d (cpu n_aln=%d)\n", i, n_aln_cpu[i]); ++fn; }
			if (!cpu && has_hit[i]) ++fp;
		}
		fprintf(stderr, "[check3] gpu has_hit=%ld  false_neg=%ld (MUST be 0)  false_pos=%ld (harmless)\n",
			gpu_hit, fn, fp);
	}

	if (getenv("DFS_NORECON")) { fprintf(stderr, "[hybrid] reconcile skipped (DFS_NORECON); gpu_hit=%ld\n", gpu_hit); return 0; }

	/* the real production hybrid: reconcile the GPU-flagged reads via the exact CPU search,
	 * reconstructing each read's pristine width from the stored flat arrays. */
	std::vector<int> n_aln_hyb(n_seqs, 0);
	std::vector<bwt_aln1_t*> aln_hyb(n_seqs, NULL);
	std::vector<int> flagged_idx;
	for (int i = 0; i < n_seqs; ++i) if (has_hit[i]) flagged_idx.push_back(i);

	/* reconcile is embarrassingly parallel (independent bwt_match_gap calls, each thread its own
	 * stack); in the streaming engine this also overlaps the next GPU batch. */
	int nT = opt->n_threads > 0 ? opt->n_threads : 1;
	int stack_maxdiff = (opt->fnr > 0.0) ? bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr) : base.max_diff;
	if (stack_maxdiff < base.max_gapo) stack_maxdiff = base.max_gapo; /* match driver stack sizing */
	auto t_rc = std::chrono::steady_clock::now();
	std::vector<std::thread> ths;
	for (int t = 0; t < nT; ++t) ths.emplace_back([&, t]() {
		gap_stack_t *st = gap_init_stack(stack_maxdiff, base.max_gapo, base.max_gape, &base);
		std::vector<bwt_width_t> wt(max_len + 1);
		for (size_t x = t; x < flagged_idx.size(); x += nT) {
			int i = flagged_idx[x], len = rp[i].len;
			for (int j = 0; j <= len; ++j) { wt[j].w = w_flat[rp[i].w_off + j]; wt[j].bid = bid_flat[rp[i].w_off + j]; }
			gap_opt_t lo = base; lo.max_diff = rp[i].max_diff;
			lo.seed_len = opt->seed_len < len ? opt->seed_len : 0x7fffffff;
			int na = 0;
			aln_hyb[i] = bwt_match_gap(bwt, len, seq_flat.data() + rp[i].seq_off, wt.data(), (bwt_width_t*)0, &lo, &na, st);
			n_aln_hyb[i] = na;
		}
		gap_destroy_stack(st);
	});
	for (auto &th : ths) th.join();
	double rcs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_rc).count();
	fprintf(stderr, "[hybrid] CPU reconcile of %zu flagged reads: %.2f s wall (%d threads)\n",
		flagged_idx.size(), rcs, nT);

	write_sai("/tmp/dfstest.hybrid.sai", opt, n_seqs, n_aln_hyb.data(), aln_hyb.data());
	if (golden) { fprintf(stderr, "[check2] hybrid .sai vs golden:\n"); md5cmp("/tmp/dfstest.hybrid.sai", golden); }

	fprintf(stderr, "[result] GPU %.0f reads/s; hybrid %s\n",
		n_seqs / (ms / 1e3), (golden ? "compared above (md5 must match)" : "(provide golden.sai to verify)"));
	return (full_cpu && fn) ? 2 : 0;
}
