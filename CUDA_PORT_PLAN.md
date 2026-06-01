# Porting `bwa aln` (BWA-backtrack) to CUDA — Research & Plan

Status: research / design draft (no code yet).
Scope target: **single-end only** (ancient-DNA reads are merged before alignment, so no
paired-end / `sampe` logic is needed). Goal is a GPU `bwa aln` that produces a `.sai`
identical (or equivalent) to CPU BWA, fast even with **seeding disabled (`-l 1024`)**,
which is the case that matters for aDNA and the case where BarraCUDA collapsed.

---

## 1. What the algorithm actually is (read the code, not the legend)

The core is `bwt_match_gap()` in `bwtgap.c`, driven by `bwa_cal_sa_reg_gap()` in `bwtaln.c`,
on top of the FM-index in `bwt.c`. Three facts dominate every design decision:

### 1a. It is a *best-first branch-and-bound* search, not a plain DFS
`gap_stack_t` is **not** a LIFO stack. It is an array of buckets indexed by alignment score:

```
score = n_mm*s_mm + n_gapo*s_gapo + n_gape*s_gape      (s_mm=3, s_gapo=11, s_gape=4)
```

`gap_pop()` always pops from the **lowest-score non-empty bucket** (`stack->best`). So the
search expands partial alignments in increasing-cost order — a Dijkstra/A\*-style priority
queue over the space of partial alignments. Combined with the *top2* early termination
(`best_score`, `best_cnt`, `max_top2`, and the `max_diff = best_diff+1` tightening once the
first hit is found), this guarantees BWA's output semantics: the best hit and a bounded count
of equally/next-best hits.

> **This is the single biggest thing BarraCUDA changed and got wrong for our use case.**
> To save memory it replaced best-first with a depth-first search keeping only "the branch
> with the local best hit score". DFS visits nodes in the wrong order, so it re-expands /
> revisits, and it does not get the priority-queue pruning for free. That is tolerable when a
> short exact seed prunes the tree to almost nothing; it is catastrophic when seeding is off.

### 1b. Each search node is a partial alignment over an FM-index SA interval
A node (`gap_entry_t`, ~32 bytes) holds: SA interval `[k,l]`, read position `i`, counts
`n_mm/n_gapo/n_gape/n_ins/n_del`, `state ∈ {M,I,D}`, and `last_diff_pos`. Expanding a node:

1. one random-access FM-index probe: `bwt_2occ4(bwt, k-1, l, cnt_k, cnt_l)` — the four base
   counts at both interval ends;
2. push up to **9 children**: 1 insertion + 4 deletion + 4 (mis)match. Hence the
   O(9^n) worst-case search space the BarraCUDA paper cites for gapped alignment (O(4^n)
   ungapped).

Pruning comes from the precomputed **`width[]`** array (`bwt_cal_width`), a lower bound on the
number of differences still required for the remaining prefix (`m < width[i-1].bid` ⇒ prune),
plus the optional **`seed_width[]`** bound. The tail of an alignment, once no diffs remain, is
finished by the exact-match routine `bwt_match_exact_alt()`.

### 1c. `-l 1024` (aDNA) deliberately removes the strongest pruning
With `seed_len ≥ read_len`, `bwa_cal_sa_reg_gap` passes `seed_w = NULL`, so **only the global
`width[]` bound applies** and `max_diff` is driven by `-n` (e.g. `-n 0.03/0.04` ⇒ several diffs
on a 35–60 bp read). The tree is large and irregular. `max_entries` (default **2,000,000**)
is the safety valve: if a read's queue exceeds it, the read is abandoned. On a CPU this is one
reused, resizable priority queue. **You cannot give each of tens of thousands of GPU threads a
2M-entry queue** — this is the crux of the whole port.

### 1d. The `aln` path needs the FM-index but NOT the suffix array
`bwt_match_gap` outputs SA **intervals** `(k,l)` in `bwt_aln1_t`; the suffix array is only
consulted later in `bwa samse` to turn intervals into coordinates. Confirmed: nothing in
`bwtaln.c`/`bwtgap.c` calls `bwt_sa`. **The GPU kernel only needs `bwt->bwt` (BWT+Occ
checkpoints), `L2[5]`, `primary`, and `cnt_table[256]`.** For a human genome that is ~0.8 GB
BWT + Occ overhead (a few GB total) — fits comfortably in ≥8 GB device memory, and the multi-GB
SA never has to touch the GPU.

### 1e. The on-disk BWT is already a GPU-friendly FM-index
`bwt.h` packs the index in 128-base buckets: each bucket = 4×`uint64` Occ checkpoints + 128×2-bit
symbols = **64 bytes = one cache line** (`bwt_occ_intv`, `__occ_aux4` with the `cnt_table`
popcount trick). One Occ query = one 64-byte load + popcounts. This is exactly the layout the
"FM-index on GPU" / "Boosting the FM-index on the GPU" papers converge on. **We can upload the
existing array essentially as-is**; the random-access latency is hidden by running many warps
concurrently (latency hiding is what GPUs are *good* at), not by avoiding the randomness.

---

## 2. Why BarraCUDA was slow with seeding off (and what to learn)

From the BarraCUDA paper (Klus et al. 2012) and follow-ups:

- **One read per thread, DFS.** Tens of thousands of threads, each its own search. Fine when
  a seed kills the tree; with `-l 1024` each thread runs a huge, deeply divergent search.
- **Whole-warp stall on the worst read.** 32 reads share a warp; the warp runs until the
  *slowest* read's tree is exhausted. aDNA read trees vary by orders of magnitude ⇒ ~31 lanes
  idle most of the time. This is the dominant loss with seeding off.
- **Per-thread stack in local memory.** A per-thread DFS stack spills to (uncoalesced) local/
  global memory; register pressure caps occupancy.
- **Fragmenting reads into 32 bp chunks + multiple kernel launches + host-side re-ranking** to
  emulate best-first. This round-trips to the host and only papers over the divergence.
- **Texture memory** for the BWT — reasonable in 2012, superseded by the L1/L2 + `__ldg`/
  read-only path today.

Takeaway: keep BWA's **best-first** semantics, and attack **load imbalance / divergence** and
the **per-thread stack**, not the FM-index layout (which is already good).

### 2a. Two reference points that frame the whole design
- **mapAD** (mpieva, Rust) — the *algorithmic* reference. It is "pure backtracking on top of the
  bidirectional FMD-index," explicitly inspired by BWA-backtrack, and is a state-of-the-art aDNA
  mapper. Two things to take from it: (i) it confirms backtracking-over-FM-index (not seed-extend)
  is the right approach for short damaged reads; (ii) it replaces BWA's flat penalties with a
  **damage-aware `SequenceDifferenceModel`** — position-dependent C→T deamination at fragment
  ends, single/double-stranded library models, and PHRED base-quality integration — and reports
  MAPQ 3 (not 0) for equally-scoring multi-mappers. Our GPU scoring should be written so this
  generalized, position/base/quality-dependent difference model can drop in where `aln_score`
  and the mismatch test currently sit (a flat-penalty special case of it). Known mapAD limits:
  high memory, no paired-end — neither matters for our single-end aDNA scope.
- **NVBIO / nvBowtie** (and the RbowtieCuda wrapper) — the *engine* reference. Its **algorithm**
  is seed-and-extend (FM-index seed ranking → batched SA-range locating → banded Myers
  verification), which is exactly the approach that *fails* for `-l 1024` aDNA, so we do **not**
  copy the algorithm. But its **GPU plumbing is the architecture we want** and proves it scales:
  it works on ~1M-read batches as a pipeline of simple, deeply parallel stages communicating
  through **ping-pong work queues** (reads still "alive" are pushed to an output queue that
  becomes the next input queue); parallelism is spread *below* the read level (many candidate
  hits per read at once); and it load-balances by **spreading variable-size SA ranges across
  threads so each thread handles one unit of work**. That is precisely the frontier/work-pool
  engine in §4 — we borrow NVBIO's plumbing and run BWA's best-first backtracking through it.

## 3. The modern CUDA toolbox (latest Programming Guide, CUDA 12/13)

- **Dynamic Parallelism** (the HAL-Inria backtracking paper's approach): a kernel can launch
  child grids for newly discovered subtrees. But device-side `cudaDeviceSynchronize()` was
  **removed in CUDA 12**; per-launch overhead and pending-launch tracking memory are
  non-trivial. Literature consensus (and NVIDIA's own guidance): for irregular tree/backtracking
  workloads, **persistent threads + a global work queue usually beat dynamic parallelism** by
  up to an order of magnitude; DP shines only when subtree granularity is large and launches are
  few. → We treat DP as a fallback, not the primary design.
- **Persistent threads + global work pool / "Broker Queue".** A fixed grid of resident warps
  that pull work items from a global queue, expand them, and push children back. This is the
  textbook fix for *exactly* BarraCUDA's failure mode: a warp that drains a small tree
  immediately grabs more work instead of idling. This is the load-balancer.
- **Warp-level primitives** (`__ballot_sync`, `__shfl_sync`, `__match_any_sync`, `__activemask`,
  `__syncwarp`) + Volta+ **Independent Thread Scheduling**: let a warp cooperate on one read's
  frontier, do warp-aggregated queue push/pop, and parallelize the 4-symbol Occ probe across
  lanes — without unsafe implicit warp-synchronous assumptions (compute masks explicitly).
- **Cooperative Groups / grid sync**, **L2 persistence** (pin hot Occ checkpoints / `cnt_table`
  in L2), **CUDA streams** to overlap host read-parsing + `width[]` computation with device
  search, and **constant memory** for `L2[5]`, `primary`, `cnt_table`, and `gap_opt_t`.

## 4. Design options considered

| Option | Idea | Pro | Con |
|---|---|---|---|
| A. Read/thread, best-first | BarraCUDA but keep priority queue | simple | per-thread 2M queue impossible; divergence |
| B. Read/**warp** | 32 lanes cooperate on one read; queue in shared mem, spill to global; parallel Occ | kills intra-warp divergence; bigger queue; parallel probes | shared mem tiny vs 2M entries; inter-warp imbalance remains |
| C. **Persistent warps + global work pool** | resident warps pull (read or node) work items, expand, push back | best load balancing for irregular trees; no host round-trips | needs a fast concurrent global queue; best-first ordering is global/approximate |
| D. **Frontier / bulk-synchronous best-first** | expand the whole batch's frontier level-by-level in SoA global arrays, bucketed by score | no per-thread stacks; huge coalesced parallelism; latency hidden by occupancy | frontier memory can blow up; per-read top2 termination needs care |

**Recommended primary architecture: D as the engine, with C's pull-pool for load balancing,
and B's warp-cooperation as an intra-node optimization.** Concretely:

- Search state lives in **global SoA frontier arrays** (one "slot" = a `gap_entry_t`-equivalent
  plus its owning read id), never in per-thread stacks. This removes the register/local-memory
  spill problem outright.
- Maintain **score buckets globally** (a small number of priority levels, like the CPU's score
  buckets). Each engine step expands the current best non-empty bucket(s) across *all* reads at
  once → preserves best-first order globally and gives millions of independent nodes to chew on.
- Each read carries its **own bound state** (`best_score`, `best_cnt`, `max_diff`, `n_aln`,
  abandoned-flag) in a small per-read struct, so top2 termination and the `max_entries` cap are
  enforced per read exactly as on CPU.
- Node expansion (the `bwt_2occ4` probe + up to 9 children) is the hot kernel: embarrassingly
  parallel, one warp can co-probe the 4 symbols, children are appended via warp-aggregated atomics
  into the next bucket.
- A **persistent-thread variant** of the same expansion can replace bulk-synchronous level steps
  if profiling shows the level barriers waste time on long-tail reads.

This keeps BWA semantics (best-first + top2 + width pruning + exact-match tail) while turning the
"irregular per-read tree" into "one big, well-load-balanced pool of nodes."

## 5. Open design decisions (need data + a couple of user calls)

1. **Bit-exact vs. semantically-equivalent output.** BWA's tie-break order among equal-score
   hits affects `c1/c2`/mapQ. Do we need byte-identical `.sai`, or "same alignments, mapQ within
   tolerance"? (Bit-exact is a much stronger constraint on bucket/expansion ordering.)
2. **Frontier memory budget & spill policy.** What global cap replaces the per-read 2M `max_entries`?
   Need the measured node-count distribution on real aDNA at `-l 1024 -n 0.03/0.04` to size it.
3. **Target GPU / compute capability** (sets shared-mem size, ITS assumptions, occupancy tuning).
4. Engine choice **D vs C** finalized by Phase-0 profiling (tree-size distribution & long-tail).

## 6. Phased plan

- **Phase 0 — Reference & profiling harness.**
  Build a standalone driver that runs the *unmodified* `bwt_match_gap` over a real merged aDNA
  batch with `-l 1024`. Dump golden inputs/outputs (BWT handle, `width[]`, reads, `gap_opt_t`,
  resulting `bwt_aln1_t` lists) as test vectors. Instrument `stack->n_entries` (peak) and total
  nodes per read to get the **node-count distribution** — this sizes the frontier and decides D vs C.

- **Phase 1 — GPU FM-index.**
  Upload `bwt->bwt`, `L2`, `primary`, `cnt_table` to device. Implement `__device__` `bwt_2occ4`,
  `bwt_occ4`, `bwt_2occ`, `bwt_match_exact_alt` against the existing 64-byte bucket layout
  (read-only/`__ldg`, `cnt_table` in constant mem). Validate Occ values bit-for-bit vs CPU.

- **Phase 2 — Single-read device search (correctness first, speed later).**
  Port `bwt_match_gap` verbatim to a `__device__` function, one read per thread, per-read priority
  queue in a global arena. Goal: **bit-exact match to Phase-0 golden vectors**. This de-risks the
  fiddly semantics (top2, `width`/`seed_width` pruning, `gap_shadow`, exact-match tail, tandem-repeat
  dedup) *before* restructuring for parallelism.

- **Phase 3 — Scalable engine.**
  Implement the §4 frontier/best-first engine (SoA frontier, global score buckets, per-read bound
  state, warp-aggregated expansion). Add the persistent-thread / work-pool variant if the long tail
  warrants it. Tune occupancy, L2 persistence, bucket granularity.

- **Phase 4 — Integration & pipelining.**
  Wire as the `bwa aln` device path emitting identical `.sai`. Keep read parsing + `width[]` on the
  host, overlapped with device search via streams/double-buffering across the existing 0x40000-read
  batches. Single-end only.

- **Phase 5 — Validation on real aDNA.**
  Compare `.sai`/downstream `.sam` vs CPU BWA on damaged short-read sets at `-l 1024 -n 0.03/0.04`;
  benchmark vs multi-core CPU BWA and (if available) BarraCUDA, especially with seeding off.

## 6a. Phase-0 findings (measured on real data — 2026-06-01)

Sample `AVA1B.combined.fq.gz` (merged aDNA, reads 30–91 bp, mean 46 bp), ref `hs37d5.fa`,
exact production command `bwa aln -l 1024 -n 0.01 -o 2`. Instrumented `bwt_match_gap`
(`#ifdef ALN_PROFILE` in `bwtgap.c`); produced a **bit-exact** `.sai` vs the clean build
(md5 `54410c75…`), so the instrumentation is verified non-perturbing. Index sizes: `.bwt`
3.0 GB (goes to GPU), `.sa` 1.5 GB (stays on host). Toolchain: CUDA 13.2, RTX 3090 (sm_86),
24 GB; clangd 22.1.6 available.

Representative 20k strided sample (matches head-10k):
- **Mapping rate 0.52%** (from the gold BAM: 20,701 / 3,948,528). This is a very-low-endogenous
  ancient sample → **99.3% of reads have zero hits**, and those unmapped reads dominate runtime.
- **No read hit the `max_entries` (2,000,000) cap.** Worst peak queue = 1.08 M entries; worst
  case = 1.27 M node expansions. Mean ≈ 11.5 k peak entries / 40 k expansions per read.
- Tree size spans **six orders of magnitude** (2 … 1.27 M expansions). Mode ≈ 32 k expansions.

**Consequences for the design (this supersedes the §4 "global best-first frontier" lean):**
1. **DFS, not a priority queue, for the bulk.** For a read with no hit within the bound — and
   with the cap never reached — best-first and DFS enumerate the *identical* node set; ordering
   only changes the post-first-hit pruning (`max_diff` tightening + top2). 99.3% of reads never
   reach that, so the priority-queue's value is ~nil for them. A DFS needs only an
   **O(read_len + gaps)-depth stack (~4 KB/thread)** — it lives in (L1-cached) local memory, not
   a multi-MB global frontier. **BarraCUDA's central memory problem disappears.** (BarraCUDA also
   chose DFS; its losses were warp divergence, 32 bp fragmentation, and host round-trips — not DFS
   itself.)
2. **Bit-exactness via a cheap hybrid.** Run the GPU DFS to exhaustion of the bounded tree
   (trivially bit-exact for the n_aln=0 majority). Flag the ~0.5–0.7% of reads that find ≥1 hit
   and re-run *only those* through the unmodified CPU `bwt_match_gap` (or a dedicated best-first
   kernel) to reproduce the exact aln list + c1/c2 + ordering. Re-processing ~0.5% of reads on
   the CPU costs seconds.
3. **Load balancing is the #1 problem, by far.** 6-orders-of-magnitude tree-size variance ⇒
   one-read-per-thread (BarraCUDA) wastes ~31/32 lanes on the slowest read in a warp. Use
   **persistent threads pulling reads from a global work queue** so a lane that finishes a tiny
   tree immediately grabs another read. This is the single most important fix.
4. **Memory budget is comfortable.** 3.0 GB index + tiny per-thread DFS stacks + read/out buffers.
   No 6–34 MB per-read frontier arenas. Easily within 24 GB even with tens of thousands of reads
   in flight; can batch the full 0.26 M-read I/O chunks.

**Revised target architecture:** persistent-thread **work-pool** over reads; each worker runs an
**iterative DFS** of the bounded search tree with a small local-memory stack, reusing the existing
64-byte-bucket FM-index `Occ` probes; reads that find any hit are flagged and reconciled by the
exact CPU/best-first path for bit-exact `.sai`. Keep best-first only where it actually pays off.

**Resolved decisions:** bit-exact `.sai` is the goal (fallback to mapAD model only if unachievable);
target sm_86 / 24 GB confirmed; engine = DFS work-pool (not global best-first frontier).
Open: exact stack layout (recompute vs store `cnt_k/cnt_l` on backtrack), batch size tuning,
and the precise "find a hit ⇒ reprocess" trigger (any hit vs score within window).

## 7. Key references

- BWA-backtrack source: `bwtgap.c`, `bwtaln.c`, `bwt.c`, `bwt.h` (this repo).
- **mapAD** — backtracking aDNA mapper on the bidirectional FMD-index (algorithm + damage model reference): https://github.com/mpieva/mapAD
- **NVBIO / nvBowtie** — GPU FM-index pipeline plumbing (engine reference): https://nvlabs.github.io/nvbio/ , https://nvlabs.github.io/nvbio/fmmap_page.html ; wrapper: https://github.com/FranckRICHARD01/RbowtieCuda
- Klus et al., *BarraCUDA — a fast short read sequence aligner using GPUs*, BMC Res Notes 2012 — PMC3278344.
- Heng Li, *Why is bwa-aln still used (for aDNA)?* — https://lh3.github.io/2024/09/28/why-is-bwa-aln-still-used
- *GPU-accelerated backtracking using CUDA Dynamic Parallelism* — https://inria.hal.science/hal-01919514/
- "FM-index on GPU: cooperative scheme…" and "Boosting the FM-index on the GPU" (PMID 26451818) — memory-access mitigation.
- NVIDIA CUDA C Programming Guide (v13.x): Dynamic Parallelism (§4.18), Cooperative Groups (§4.4),
  warp-level primitives & ITS, L2 persistence (§4.13) — https://docs.nvidia.com/cuda/cuda-programming-guide/
- NVIDIA, *Using CUDA Warp-Level Primitives* — https://developer.nvidia.com/blog/using-cuda-warp-level-primitives/
- *A Study of Persistent Threads Style GPU Programming* (GTC) — persistent threads vs dynamic parallelism for irregular work.
