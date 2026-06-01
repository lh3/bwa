# Optimization roadmap (user ideas + critical assessment)

Current standing: `cuda/dfstest.cu`, 9934 reads/s on sub100k = ~2.3x the 16-core CPU, bit-exact.
Bottleneck = global-memory DFS stack (raising occupancy hurt; 2occ4 fast path didn't help).

## 1. Warp-level parallel-branch DFS + shared-memory stack  — ADOPTED (next move)
One warp per read; the per-warp stack lives in **shared memory** (spill to global only on
overflow); the warp expands multiple frontier nodes per iteration (32 lanes) instead of one
thread serializing one path. Directly attacks the diagnosed bottleneck: stack moves off global
memory, and stack ops become warp-aggregated shared-memory ops.
- Occupancy note: one-read-per-warp gives ~2.3k reads in flight (vs ~73k now) but the SAME
  ~28 warps/SM, so latency hiding is preserved while each warp keeps ~32 FM-index probes in
  flight. Lane utilization is fine here because the 99.3% zero-hit reads still have ~40k nodes
  to spread across lanes (they are deep, not shallow).
- Must stay bit-exact: has_hit detection is order-independent, so any traversal covering the same
  bounded node set (and early-exiting on any hit) is valid. Keep the work-budget -> CPU reconcile.
- Combines user ideas #1 (warp branching) and #2 (register + SMem stack window).

## 2. Register + shared-memory stack window — ADOPTED (part of #1)
Top of stack in registers (already done via in-register-continue); next slice in per-warp shared
memory; global only on SMem overflow. 16 entries x 32 lanes x 20 B ~= 10 KB/warp — fits.

## 3. "Zero-hit" exact-seed fast path — REJECTED (breaks correctness for -l 1024)
Proposal: if the first l (e.g. 32) bases don't match exactly, mark the read unmappable and skip.
**This is incorrect for our use case.** With `-l 1024` seeding is DISABLED precisely because
ancient-DNA reads carry mismatches/deamination anywhere, including the first 32 bp. A read can map
with several mismatches in that region; an exact-seed filter would wrongly skip it -> false
negatives -> non-bit-exact `.sai`. (This is the whole reason aDNA pipelines use -l 1024.)
- A *conservative* prefilter that only skips reads PROVABLY unmappable within max_diff (e.g. a
  q-gram/minimizer count lower bound) could be bit-exact-safe and is worth exploring later, but the
  exact-seed version is not safe. Deferred, low priority.

## 4. Vectorized Occ loads — PARTIAL / minor
The BWT bucket is already a 64 B cache line; `__occ_aux4` already unpacks 4 nt per 32-bit word via
cnt_table (so the "unpack 4 per load" is done). Possible minor win: fetch the 4xuint64 checkpoint
as two `int4`/`uint4` `__ldg` instead of 4 scalar loads. Low priority; do after #1.

## 5. Pipelining: overlap CPU reconcile with GPU via streams — ADOPTED (for the full-file engine)
Currently the harness is synchronous. The production full-file path must: stream read batches,
run the DFS kernel on stream A while a host thread pool runs `bwt_match_gap` on the previous
batch's flagged reads, double-buffering. The 1-thread 2.37 s reconcile becomes ~0.15 s on 16
threads and overlaps the kernel. Do when building the end-to-end streaming engine.

## 6. #pragma unroll / template on read length — LOW priority
The DFS is a tree, not a length-bounded loop, so templating read length gives little. The bounded
inner loops (within-bucket occ, 0..8) are already unrollable. Revisit only if profiling shows it.

## Order of execution
1) #1+#2 warp-cooperative DFS with shared-memory stack (biggest lever).
2) #5 streaming + multithreaded overlapped reconcile (end-to-end wall-clock).
3) #4 vectorized checkpoint load; then re-profile.
Re-check CUDA docs + NVIDIA forums at each step (cuda/REFERENCES.md).
