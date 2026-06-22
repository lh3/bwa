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

## Phase 2 step 2 — TWO-LEVEL stack engine (`DFS_ENGINE=warp2`) — DONE, bit-exact
Per-warp shared top-window (CAP_SM) + pre-allocated per-warp GLOBAL backing; warp-parallel
coalesced spill/unspill of 128-entry chunks (invariant: global=oldest bottom, shared=newest top).
Bushy reads spill to global and STAY ON THE GPU instead of flagging to CPU.

Result (sub100k, bit-exact md5 eecf35c1): **flag rate driven to ~0.015%** (vs 0.34% single-level)
-> CPU reconcile negligible (520 reads, 0.61 s). CAP_SM sweep: 256->15.9k r/s (16 w/SM),
384->18.8k (12), 512->21.5k (8), 768->22.1k (4).

**KEY FINDING — we are at the FM-index occ4 ceiling, NOT occupancy-bound.** Occupancy from 4 to 16
warps/SM does NOT raise throughput; all engines/configs plateau at **~22k reads/s**. This work is
~8e9 occ4 (4.1e9 node-pops x ~2 probes); at the measured 2.28 G-occ4/s that is a ~3.6 s floor vs
~4.5 s actual (~80% of ceiling, ~5x the 16-core CPU). The 40-50k target is not reachable on this
GPU because the kernel is bound by FM-index probe throughput, not latency hiding. Raw speedup
beyond this needs FEWER occ4 per node (algorithmic; constrained by bit-exactness) or a faster
FM-index (k-step / occ caching) -- diminishing returns. Default CAP_SM=512.

**warp2 is the engine to ship**: same ceiling speed as warp1 but ~0 CPU reconcile tail, which keeps
the streaming engine clean. Next: full-file streaming (#2) -- the deliverable, not a speed lever.

## Phase 2 step 2 — STREAMING full-file tool `bwa-aln-gpu` — DONE, bit-exact end-to-end
`cuda/aln_gpu.cu` + `cuda/dfs_engine.cuh` (extracted warp2). Streams the FASTQ in bwa's native
0x40000-read chunks; per chunk: MT preprocess (bwt_cal_width + complement + per-read max_diff) ->
warp2 GPU has_hit -> MT CPU reconcile of flagged/hit reads -> write records in read order.
`make bwa-aln-gpu`; `./bwa-aln-gpu [-l -n -o -t -f] <ref.fa> <in.fq>`.

**FULL FILE (AVA1B.combined.fq.gz, 3,948,528 reads, hs37d5, -l 1024 -n 0.01 -o 2, RTX 3090):**
- bwa-aln-gpu: **175.9 s = 22,445 reads/s** (preprocess 1.4s, gpu 156.5s, reconcile 16.7s/16thr, io 0.1s);
  0.531% reconciled on CPU. Output .sai = 16.85 MB.
- CPU `bwa aln -t 16` golden on the same file: **875 s**.
- => **4.97x end-to-end speedup, BIT-EXACT**: both .sai = 16,852,676 bytes, md5
  `a09a26dd894690031da42a00862f7a3d` (GPU == CPU). Full 3.95M-read file verified byte-identical.

## Phase 3 — FUSED `alnse` (aln+samse, multithreaded) — `bwa-aln-gpu -S`
Added `-S` (SAM out) + `-r RG` to `cuda/aln_gpu.cu`: fuses aln and samse in one command
(GeoGenetics-style alnse). No `.sai` on disk; the FASTQ is read ONCE (samse reuses the loaded
`seqs[]`); `bns`/SA/`pac` loaded once. Per chunk after GPU aln + reconcile:
- serial `bwa_aln2seq_core` (preserves `drand48` repeat-hit-selection order; the only samse RNG --
  `lrand48` at bntseq.c:266 is index-build only),
- **MT** `bwa_cal_pac_pos_core` SA-lookup (RNG-free) across `n_threads`,
- serial `bwa_refine_gapped` + `bwa_print_sam1` (ordered, bit-exact output).
Header via `bwa_print_sam_hdr` (@HD/@SQ/@RG) + a `@PG` for bwa-aln-gpu.

Validated: fused SAM **alignment records are byte-identical** to `bwa samse` on the bit-exact `.sai`
(sub100k md5 `ee2dd6ef…`; header differs only in the expected @PG CL). `... -S | samtools sort` -> BAM.

**FULL FILE fused alnse (3,948,528 reads, RTX 3090): 174.2 s = 22,672 reads/s**, SAM records
BIT-EXACT vs reference samse (md5 `e2a6e1c9…`). preprocess 1.4s, gpu 153.9s, reconcile 16.8s, io 0.8s.
The fused run (174 s) is FASTER than aln-only(176s)+separate samse(20s)=196s -- samse folds in for
free (single FASTQ read, MT SA-lookup, samse compute hidden). vs CPU `bwa aln -t16`+samse = 895 s
=> **~5.1x end-to-end, one command, no intermediate .sai, byte-identical alignments.**
Usage: `bwa-aln-gpu -S -r '@RG\t...' ref.fa reads.fq.gz | samtools sort -O bam -o out.bam -`

## Phase 3b — CPU/GPU overlap (#5) — DONE, bit-exact
`cuda/aln_gpu.cu` restructured into a producer/consumer pipeline: the main thread owns the GPU
(read -> MT preprocess -> upload -> kernel -> download has_hit, serial on one stream so the shared
global backing is never raced); each chunk is handed to ONE in-order finisher thread that does the
CPU reconcile + output (.sai or samse), running concurrently with the next chunk's GPU work. Single
ordered consumer preserves drand48/output order -> bit-exact. (Multi-GPU-ready: replicate the GPU
producer stage per device + round-robin chunks; keep one ordered finisher.)

Full file (3.95M reads, fused alnse SAM): **174.2 s -> 160.9 s = 24,541 reads/s** (overlap hid ~13s
of the ~18s CPU tail), records BIT-EXACT (md5 e2a6e1c9). Both modes re-verified bit-exact under the
pipeline (sub100k .sai eecf35c1, SAM records identical). vs CPU `bwa aln -t16`+samse = 895 s ->
**~5.56x end-to-end.**

## Phase 3c — native `bwa gpualn` subcommand — DONE
`cuda/aln_gpu.cu`'s entry is now `extern "C" int bwa_alnse_gpu(int,char**)`; `main.c` dispatches
`bwa gpualn` under `#ifdef HAVE_CUDA`. Builds:
- `make`           -> CPU-only `bwa` (NO CUDA dependency; `gpualn` absent) -- fork still builds anywhere.
- `make bwa-gpu`   -> CUDA `bwa` with the `gpualn` subcommand (main.c -DHAVE_CUDA + nvcc aln_gpu.o,
                      linked via nvcc; reuses the existing AOBJS for samse/aln CPU functions).
- `make bwa-aln-gpu` -> standalone tool (unchanged; -DALN_GPU_MAIN).
Usage:
  bwa gpualn [-l 1024 -n 0.01 -o 2 -t 16] ref.fa reads.fq.gz > out.sai           # .sai (like bwa aln)
  bwa gpualn -S -r '@RG\t...' ref.fa reads.fq.gz | samtools sort -O bam -o out.bam -   # fused alnse -> BAM
Verified: `bwa gpualn` .sai md5 == golden (eecf35c1); `bwa gpualn -S` SAM records == `bwa samse`;
@PG is proper bwa provenance (ID:bwa). CPU-only `bwa` regression-checked (no gpualn, no CUDA link).

## Validation vs the PRODUCTION BAM (/home/dnastorage/aDNApipeline/AVA1B_aln/)
Reference: `AVA1B_aln.short.bam` = production `bwa aln -l 1024 -n 0.01 -o 2` + `samse` + sort, built
with bwa **0.7.18-r1243-dirty** on the same 3,948,528 reads (20,701 mapped, 0.52%). Compared to
`bwa gpualn -S` output (this tree): extracted mapped records, dropped the RG:Z tag (production used
a different @RG ID), sorted by QNAME:
- mapped count identical (20,701 vs 20,701); mapped read SET identical (same QNAMEs) -> no read
  changed mapped/unmapped status.
- **ALL 20,701 mapped records + ALL tags (FLAG/RNAME/POS/MAPQ/CIGAR, XT/NM/X0/X1/XM/XO/XG/MD)
  byte-identical.** Holds across the 0.7.18-dirty -> 0.7.19 gap (aln/bwtgap core unchanged).
Only the @RG ID and the SAM @PG line differ -- nothing in the alignments. End-to-end equivalence to
the production pipeline confirmed.

## STATUS: GOAL ACHIEVED
GPU `bwa aln` (BWA-backtrack) for ancient DNA at `-l 1024 -n 0.01 -o 2`, single-end:
**~5x the 16-core CPU on a full real file, byte-identical .sai, GPU-bound at the FM-index ceiling.**
Engine: warp-cooperative two-level-stack DFS (`cuda/dfs_engine.cuh`) + MT CPU reconcile of ~0.5%;
tool: `bwa-aln-gpu` (`cuda/aln_gpu.cu`). The earlier BarraCUDA failure mode (slow with seeding off)
is resolved: this is fastest precisely in the seeding-off aDNA regime.

Possible future work (optional, diminishing/oncosting returns): double-buffered stream overlap of
preprocess/reconcile with the kernel (~10% — preprocess+reconcile is only ~18s of 176s); faster
FM-index (k-step / L2-pinned occ) to lift the ~22k occ4 ceiling; multi-GPU; wrap as a real
`bwa aln` subcommand/flag in main.c; libdeflate FASTQ decode (not currently a bottleneck, ~1.2s).
