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

Next: Phase 2 — device iterative DFS reproducing `bwt_match_gap`'s node set, validated bit-exact on
sub10k, then wrap in the persistent-thread work-pool (the actual speedup lever).
