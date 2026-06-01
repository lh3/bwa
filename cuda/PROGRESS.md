# CUDA port — progress log

Goal: massively speed up `bwa aln` (BWA-backtrack) for ancient DNA, producing a **bit-exact**
`.sai` vs CPU bwa. Single-end only. Target GPU: RTX 3090 (sm_86, 24 GB). Gold-standard command:
`bwa aln -l 1024 -n 0.01 -o 2`. See `../CUDA_PORT_PLAN.md` for the full research + architecture.

## Test assets (not committed; under `../test_data/`, regenerate as needed)
- Reads: `/home/dnastorage/aDNApipeline/AVA1B/AVA1B.combined.fq.gz` (merged aDNA, 30–91 bp, mean 46).
- Index: `/home/dnastorage/aDNApipeline/hs37d5.fa` (.bwt 3.14 GB; .sa NOT needed on GPU).
- Subsets: `sub10k.fq`, `sub100k.fq` (deterministic heads), `strided20k.fq` (representative).
- Golden CPU `.sai` for sub10k: md5 `54410c755509b8389ae992953dd3476c`.

## Phase 0 — baseline & profiling — DONE
Instrumented `bwt_match_gap` (`#ifdef ALN_PROFILE` in `bwtgap.c`; build `bwa-prof`). Verified
bit-exact `.sai` vs clean build. Findings (representative 20k sample):
- Mapping rate **0.52%** → **99.3% of reads have zero hits** and dominate runtime.
- **No read hit the 2M `max_entries` cap** (worst peak queue 1.08 M, worst 1.27 M expansions; mean ~40 k).
- Tree size spans **6 orders of magnitude**.
- ⇒ For zero-hit reads (cap never reached), DFS and best-first visit the IDENTICAL node set; the
  priority queue is unnecessary for the bulk. Architecture: persistent-thread **work-pool** + per-read
  iterative DFS (small local stack); flag the ~0.5% reads that find a hit and reconcile via the exact
  CPU best-first path for bit-exactness. Load balancing is the #1 problem.

## Phase 1 — device FM-index (Occ) — DONE
`cuda/fm_device.cuh`: device mirrors of `bwt_occ4` / `bwt_2occ4` / `bwt_occ` / `bwt_match_exact`,
operating on the existing 64-byte-bucket layout; `cnt_table`/`L2` in constant memory; `.bwt` uploaded
to global memory (3.14 GB); reads via `__ldg`. `cuda/fmtest.cu` validates against CPU.
Build: `make fmtest`  ·  Run: `./fmtest <ref.fa> [reads.fq] [n_random_occ]`.

Result on hs37d5 + sub10k:
- Test A (random Occ4, 10 M probes): **PASS, 0 mismatches, 2.28 G-Occ4/s**.
- Test B (full backward exact-search, 10 k real reads): **PASS, 0 mismatches** (42 full-length hits).

## Phase 2 step 1 — bit-exact device DFS (hit detection) — DONE (correctness), perf baseline only
`cuda/dfstest.cu` + `cuda/fm_device.cuh` (added `d_bwt_match_exact_alt`, `d_int_log2`). Host reuses
bwa's real I/O + `bwt_cal_width` + complement + per-read `max_diff`, then the device DFS
(`d_dfs_has_hit`) reproduces `bwt_match_gap`'s bounded node set and reports has_hit; the ~0.5% hit
reads are reconciled by the exact CPU `bwt_match_gap` for a bit-exact `.sai`.
Build `make dfstest`; run `./dfstest <ref.fa> <reads.fq> [golden.sai]` (replicates `-l 1024 -n 0.01 -o 2 -t 16`).

Result on sub10k vs golden (md5 `54410c75…`):
- check1 CPU-only `.sai` == golden ✓ (harness == driver)
- check3 GPU has_hit=52, **false_neg=0**, false_pos=0 ✓ (DFS detection complete & exact)
- check2 **hybrid `.sai` == golden** ✓ — **bit-exact achieved**
- stack_overflow=0 (cap 8192 entries sufficient).

**Performance baseline (naive, one-read-per-thread grid-stride, P=8192, global-memory stack):**
**100 reads/s** — vs ~3,300 reads/s on the 16-core CPU → **~33x SLOWER**. This reproduces the known
BarraCUDA-style failure mode. Diagnosed causes for step 2:
1. **Warp divergence**: 32 reads per warp with 6-orders-of-magnitude tree-size variance → warp waits
   for its slowest read; one 1.27M-node read stalls 31 idle lanes.
2. **Low occupancy / poor latency hiding**: only 8192 threads; FM-index probes are ~500 ns random
   global reads that need far more in-flight warps to hide.
3. **Pathological stack memory pattern**: per-thread stacks at `slot*cap` stride (196 KB apart) →
   every push/pop is 32 uncoalesced cache lines; ~5-9 stack writes per node dominate traffic.
4. `d_bwt_2occ4` = 2 full occ4 (no shared-bucket fast path).

Reference ceiling: ~1e9 occ4 for the whole 10k set; at the measured 2.28 G-Occ4/s that's ~0.44 s if
perfectly parallel → the naive run is ~227x off the FM-index ceiling, i.e. essentially all loss is
scheduling/memory, not the index. Headroom is enormous.

## Phase 2 step 2 — performance — IN PROGRESS (bit-exact maintained throughout)

Changes to `cuda/dfstest.cu`: coalesced SoA stack (`[sp*P+slot]`), persistent-thread work-pool
(atomic work counter), occupancy-sized launch (82 SM x 7 blk x 128 = 73472 slots), per-read **work
budget** (default 2M pops via `DFS_BUDGET`) that flags pathological reads to the exact CPU path,
and an opt-in full CPU reference (`DFS_FULLCPU=1`, auto-on for <=20k reads) — large runs reconcile
only the flagged/hit reads from stored flat data (the real production hybrid).

Findings (the key one): the naive 100s was **dominated by ONE read** doing 33.1M node-pops on a
single thread. The CPU's `max_entries=2M` frontier cap makes `bwt_match_gap` bail on such reads
(returning n_aln=0), so CPU stays fast; the GPU had no equivalent. Coalesced stack + occupancy +
work-pool ALONE changed nothing (still ~100s) because time was that one serial read, not
scheduling. The **work budget** fixed it — and it's free for bit-exactness because flagged reads
are reconciled by the exact CPU search anyway.

Measured (RTX 3090, hs37d5, `-l 1024 -n 0.01 -o 2`, budget 2M):
| set | naive | +budget | flagged | bit-exact |
|-----|-------|---------|---------|-----------|
| sub2k  | 101.8 s (20 r/s) | 4.5 s (440 r/s) | 2 (0.1%) | yes |
| sub100k | — | **12.08 s (8276 r/s)** + 2.37 s reconcile (520 reads, 1 thread) | 11 budget + ~509 hits | **yes (md5 eecf35c1)** |

CPU baseline on sub100k, same box: **`bwa aln -t 16` = 23.47 s (4261 r/s)**.
=> GPU kernel is **~1.9x the full 16-core CPU** already (reconcile parallelizes/overlaps away);
**~83x** over the naive GPU baseline. Throughput 338M pops/s ~= 38% of the FM-index occ4 ceiling,
so meaningful headroom remains.

Next levers: `bwt_2occ4` shared-bucket fast path (~halve probes), cut register pressure to raise
occupancy beyond 7 blk/SM, cheaper exact-match tail, then batch streaming over the full file +
multi-thread the reconcile. Goal: push well past 2x toward the occ4 ceiling.
