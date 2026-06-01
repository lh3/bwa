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

### step-2 optimization log (sub100k, RTX 3090, bit-exact throughout)
| change | reads/s | note |
|--------|---------|------|
| work-budget engine (above) | 8276 | baseline for this log |
| + `bwt_2occ4` shared-bucket fast path | 8224 | **no change** -> not occ-load bound |
| maxrregcount 48 (40 warps) | 6445 | **slower** |
| maxrregcount 40 (48 warps) | 5243 | **slower** -> NOT occupancy/latency bound |
| + in-register-continue DFS (stack only for siblings, pop only on backtrack) | **9934** | +20%; 69 regs |

Diagnosis (consulted CUDA Ampere tuning guide + NVIDIA forums + Volkov latency-hiding + the 2025
arXiv "N-Queens GPU iterative DFS" paper): the kernel is **bound by the global-memory DFS stack**,
not occ loads and not occupancy. Evidence: the 2occ4 fast path didn't help, and *raising* occupancy
*hurt* (more concurrent threads -> the per-thread live stack working set exceeds the 6 MB L2 ->
DRAM). The literature is explicit: one-thread-per-task DFS with a global stack scales poorly;
assign a subtree to a group of threads sharing fast memory, and keep the stack in shared/registers.
The in-register-continue change applies the register part (descend in registers, push only siblings,
read stack only on backtrack) for +20%.

Current standing: **9934 reads/s = ~2.3x the 16-core CPU (4261 r/s), ~27x one core (363 r/s),
99x over the naive GPU baseline (100 r/s); bit-exact.** Throughput ~408M pops/s ~= 45% of the
occ4 ceiling.

Next big levers (toward "massive"): (1) **one-read-per-warp / subtree-per-warp cooperation** with the
hot stack in shared memory (the literature's recommended structure; should cut both divergence and
the L2-thrashing stack); (2) shared-memory hot-stack window backing the global spill; (3) full-file
streaming with overlapped, multithreaded CPU reconcile. References saved in cuda/REFERENCES.md.

### step-2 cont. (data that pins down the next move)
- **DFS stack DEPTH (sub100k): mean 111, max 386.** Histogram concentrated at <=128 (90,876 reads),
  <=256 (7,612), <=512 (223). So a **512-entry stack covers every read** -> a per-warp shared-memory
  stack is feasible (~512*20 B = 10 KB/warp). Reduced cap 1024->512 (0 overflow).
- **Block-local stack re-striding ([blockIdx*blockDim*cap + d*blockDim + tib]): NO win** (9165 vs
  9934 r/s, bit-exact). Negative result -> the global stack is bound by traffic VOLUME/bandwidth,
  not layout/locality. Re-striding can't reduce volume; only removing global stack traffic can.
- A per-*thread* full shared stack can't fit (128 thr * 386 * 20 B ~= 1 MB/block). So the shared
  stack **requires one read per warp** (one stack/warp) -> which forces lane cooperation = Idea #1.
  Decision: implement the warp-cooperative engine with a per-warp shared-memory stack.

## Phase 2 step 2 — WARP-COOPERATIVE engine (Idea #1+#2) — DONE (kernel), bit-exact
`DFS_ENGINE=warp` in `cuda/dfstest.cu`: one read per warp; 32 lanes co-explore the tree as a
frontier; per-warp stack in SHARED memory (SoA). Wave = pop up to 32 top nodes -> expand in
parallel -> warp prefix-sum compact children -> push to shared stack; hit via `__any_sync`;
budget/overflow -> CPU reconcile. Removes global stack traffic.

**sub100k GPU kernel: 29,708 reads/s** (CAP=640, 4 warps/blk) = **3x the thread engine (9934),
~7x the 16-core CPU (4261 r/s)**; bit-exact (md5 eecf35c1). sub2k: 6304 vs 547 r/s (11.5x) —
monster reads now spread across 32 lanes instead of serializing on one thread.

### CURRENT BOTTLENECK (documented) — CPU reconcile of overflow-flagged reads
The warp wave pops 32 nodes and pushes ALL their children, so the live frontier is **BFS-wide**,
not bounded by the serial DFS depth (386). With the shared stack CAP=640 it **overflows for 4.6%
of reads**, which are flagged to the CPU. Reconciling those 4,835 bushy reads on **one CPU thread
takes 84 s** — it now *dwarfs* the 3.4 s GPU kernel (end-to-end 3.4 + 84 = 87 s, worse than CPU's
23.5 s). So the bottleneck has MOVED from "GPU global-stack traffic" to "CPU reconcile of
overflow-flagged reads"; the GPU itself is no longer the limiter.

Levers (both needed):
1. **Cut the overflow flag rate** toward the real-hit floor (~0.5%): raise CAP (bounded by shared
   mem -> occupancy), and/or bound the warp frontier so it stays near the serial DFS depth (e.g.
   pop fewer per wave / drain depth-first when near capacity) so CAP=512 suffices.
2. **Multithread + GPU-overlap the reconcile** (Idea #5): the flagged-read CPU work is embarrassingly
   parallel; 16 threads -> ~5x, and it can overlap the next GPU batch via streams. Even at 4.6%,
   84 s/16 ~= 5.3 s; with flag rate driven to ~0.5% it becomes negligible and the engine is
   GPU-bound at ~29.7k r/s (~7x CPU).
Added `DFS_NORECON` to skip reconcile during perf sweeps.

### Resolution of the bottleneck (CAP tuning + multithreaded reconcile)
CAP sweep (sub100k, warp, 4 warps/blk): CAP=640 -> 4.6% flagged; CAP=1024 (80 KB) -> 0.28%;
CAP=1216 (95 KB) -> 0.065%. Larger CAP = GPU does the bushy reads itself (fewer flags) but is a
bit slower (more GPU work); the CAP=640 "34k r/s" was illusory (it punted bushy reads to a 84 s
CPU reconcile). Default CAP set to 1024.
**Multithreaded the reconcile** (std::thread, per-thread gap_stack; bwt_match_gap is thread-safe):
751 flagged reads reconciled in **1.57 s wall on 16 threads** (was 12.1 s single-thread).

**Engine standing (sub100k, RTX 3090, bit-exact md5 eecf35c1):**
- Warp GPU kernel: **21,432 reads/s** (CAP=1024).
- End-to-end GPU + 16-thread reconcile (sequential): 4.67 + 1.57 = 6.24 s = **16,026 reads/s**.
- With streaming overlap (reconcile hidden behind next GPU batch): GPU-bound ~**21,432 r/s**.
- CPU `bwa aln -t 16` = 4,261 r/s on the same box -> **~3.8x (sequential) to ~5x (overlapped)**.
- vs the naive GPU baseline (100 r/s): ~210x.

Remaining levers: (a) shared(small)+global-backing two-level stack -> high occupancy (currently
only 4 warps/SM due to the 80 KB shared stack) without flagging, to push the GPU beyond 21k;
(b) full-file streaming engine with double-buffered batches + reconcile overlapping the next
kernel (to realize the GPU-bound rate and process the whole 3.95M-read file end-to-end).

### Drain-down + carry register-continue (Lever 1) + the end-to-end ceiling
Implemented dynamic pop-count: WAVE mode pops up to 32 nodes but only while room=(CAP-sp)/9 allows
every node's <=9 children (so wave mode NEVER overflows); DRAIN mode (near full) pops 1 and keeps a
child in registers (carry) so single-chains never touch the stack. Gated carry to sp>=CAP/2 to keep
the wave ramp-up for normal reads. MAX_CHILDREN confirmed = 9 (state M: 1 ins + 4 del + 4 mm).
All bit-exact (sub2k false-pos 90->1).

CAP sweep (carry-drain, sub100k, GPU-only): CAP=256 -> 27% flag, 59k r/s (16 warps/SM);
384 -> 5.6%, 35.7k (12 w/SM); 512 -> 2.2%, 21k (8 w/SM); 768 -> 0.34%, 22.2k (4 w/SM).
With the overlapped 16-thread reconcile (~444 flagged-reads/s), end-to-end is reconcile-bound for
small CAP and GPU-bound (~22k) for CAP>=768. **KEY: the high-occupancy small-CAP GPU speed is
illusory — it punts the bushy reads (2-5%) to the CPU. Those reads cost similarly on GPU (low
occupancy) or CPU (reconcile), pinning end-to-end at ~22k r/s (~5x the 16-core CPU) regardless of
CAP.** Default CAP set to 768 (0.34% flag, GPU-bound 22.2k).

To break the ~22k ceiling, the bushy reads must run ON THE GPU AT HIGH OCCUPANCY -> the two-level
**shared(small, e.g. 256 -> 16 warps/SM) + per-warp global-backing** stack: 95%+ of reads stay in
shared (59k base rate); the 2-5% bushy reads spill their deep frontier to global (handled on GPU,
NOT flagged). Only true 2M-budget reads -> CPU. This is the next build.
