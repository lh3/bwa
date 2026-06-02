# Development notes (LLM-assisted) & references

`bwa gpualn` — the CUDA port of BWA-backtrack in this fork — was built as a
**human-directed, AI-assisted** effort. This note records how it was made and
the literature it drew on, for transparency and reproducibility.

## How it was built

A three-way collaboration:

- **Human lead — teepean (project owner).** Set the goal (a GPU `bwa aln` for
  ancient DNA that is fast with seeding disabled, `-l 1024`), provided the
  hardware (RTX 3090), the reference/index (hs37d5) and real aDNA samples, the
  production pipeline, and the references below; steered the architecture; and
  validated results against real data.

- **Implementation, profiling & measurement — Anthropic's Claude**
  (Claude Code, Opus 4.8). Wrote the device FM-index (`Occ`/`2occ4`), the
  warp-cooperative two-level-stack DFS engine, the hybrid GPU-detect + CPU-exact
  reconcile, the streaming CPU/GPU overlap, and the `bwa gpualn` subcommand; ran
  the profiling and benchmarks; and pushed back with data when a proposal didn't
  hold — e.g. empirically **disproving** the "40k reads/s" target by showing the
  kernel is at the GPU's random-access FM-index bandwidth ceiling (not
  occupancy- or scheduling-bound), and **rejecting** an exact-seed prefilter
  because it would break bit-exactness under `-l 1024`.

- **Architectural sparring & optimization ideas — z.ai's GLM-5.1**
  (relayed through teepean). Proposed the design hypotheses that drove several
  key steps: the warp-level parallel-branch traversal, the drain-down heuristic
  with warp-level in-register-continue, the two-level shared+global stack with
  warp-parallel coalesced spill, and the CPU/GPU overlap and multi-GPU framing.

The working loop was: **GLM-5.1 (via teepean) proposed an optimization →
Claude implemented, measured, and reported it (often with a caveat or a
disproof) → next idea.** The decisive wins — the per-read work-budget cap, the
two-level stack that eliminated CPU-flagging, and the ~5×, bit-exact result —
came out of that back-and-forth. Nothing here is asserted without measurement:
output is **bit-exact** to CPU `bwa aln`/`samse`, verified up to a 47-million-read
sample (944,922 mapped records byte-identical).

The full design/profiling log is in [PROGRESS.md](PROGRESS.md); the
optimization decisions (adopted and rejected, with reasons) are in
[OPTIMIZATION_IDEAS.md](OPTIMIZATION_IDEAS.md); CUDA-tuning specifics are in
[REFERENCES.md](REFERENCES.md).

## References

### Provided by teepean (at project start)
- NVIDIA CUDA C++ Programming Guide — https://docs.nvidia.com/cuda/cuda-programming-guide/index.html
- GPU-accelerated backtracking using CUDA dynamic parallelism (INRIA / HAL) — https://inria.hal.science/hal-01919514/
- "Trouble porting backtracking algorithm to CUDA kernel" (StackOverflow) — https://stackoverflow.com/questions/64084225/trouble-porting-backtracking-algorithm-to-cuda-kernel
- "Non-recursive backtracking in CUDA" (NVIDIA Developer Forums) — https://forums.developer.nvidia.com/t/non-recursive-backtracking-in-cuda/15351
- BackTrack 4 CUDA Guide (Offensive Security) — https://www.offsec.com/documentation/backtrack-4-cuda-guide.pdf
- A GPU-Based Backtracking Algorithm for Permutation Combinatorial Problems (ResearchGate, 310810766)
- parallel-sudoku-solver (GitHub) — https://github.com/vduan/parallel-sudoku-solver
- "GPU strategies for fine-grained combinatorial" (figure, academia.edu/figures/14331374)

### Provided by teepean (during the work)
- **mapAD** — backtracking aDNA mapper on the bidirectional FMD-index (algorithm/damage-model reference) — https://github.com/mpieva/mapAD
- **RbowtieCuda / NVBIO** — GPU FM-index pipeline (engine reference) — https://github.com/FranckRICHARD01/RbowtieCuda , https://nvlabs.github.io/nvbio/
- **GeoGenetics/bwa** — `alnse` (fused aln+samse) — https://github.com/GeoGenetics/bwa

### Found during the work
- Klus et al., *BarraCUDA — a fast short read sequence aligner using GPUs*, BMC Res Notes 2012 — PMC3278344 (the prior attempt; the cautionary baseline)
- Heng Li, *Why is bwa-aln still used (for ancient DNA)?* — https://lh3.github.io/2024/09/28/why-is-bwa-aln-still-used
- *FM-index on GPU: a cooperative scheme…* and *Boosting the FM-index on the GPU* (PMID 26451818) — memory-access mitigation
- V. Volkov, *Understanding Latency Hiding on GPUs* — UCB EECS-2016-143 (ILP vs occupancy)
- *High-Performance N-Queens Solver on GPU: Iterative DFS with Zero Bank Conflicts* — arXiv 2511.12009 (shared-memory DFS stack)
- *A Study of Persistent Threads Style GPU Programming for GPGPU Workloads* (GTC) — persistent threads vs dynamic parallelism
- NVIDIA CUDA Ampere GPU Architecture Tuning Guide — https://docs.nvidia.com/cuda/ampere-tuning-guide/
- NVIDIA, *Using CUDA Warp-Level Primitives* — https://developer.nvidia.com/blog/using-cuda-warp-level-primitives/

## License
This GPU code is **GPLv3** (same as bwa). The `gpualn` port is a derivative work
of bwa's GPL-licensed core — it `#include`s and links `bwtgap.c` (`bwt_match_gap`),
`bwtaln.c` (`bwt_cal_width` etc.), `bwase.c` (the `samse` functions) and
`bwaseqio.c`, which carry no separate MIT header and so fall under bwa's GPLv3
(`COPYING`). The new CUDA source files (`cuda/fm_device.cuh`, `cuda/dfs_engine.cuh`,
`cuda/aln_gpu.cu`, `cuda/dfstest.cu`, `cuda/fmtest.cu`) carry the standard GPLv3
header; the full license text is in `COPYING`. (Some bwa files are MIT-licensed,
but a work that derives from and links the GPL parts is GPLv3 as a whole.)
This is an engineering reading of the file headers, not legal advice.

## Reproducibility / disclosure
LLM-generated and LLM-assisted code and analysis were reviewed and empirically
validated against stock `bwa` before acceptance (bit-exact `.sai`/SAM checks at
every step). Commits authored during this work carry a `Co-Authored-By:` trailer
for the AI assistant. Treat the design rationale here as a record of *why*
choices were made; the authoritative behavior is "matches CPU `bwa aln`."
