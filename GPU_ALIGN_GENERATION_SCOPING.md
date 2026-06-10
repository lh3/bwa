# Scoping: GPU-side alignment generation (eliminating the reconcile floor)

## Why this is on the table

On the 64–69 bp (d=5) band the GPU *existence* kernel runs at ~24,500 reads/s — **7.2× the
32-thread CPU** (3,387 reads/s). But end-to-end gpualn only *ties* the CPU there, because
every read with a hit (~37% of the band) has its alignment **regenerated on the CPU** by
`bwt_match_gap` — the ~120 s reconcile floor that no budget setting can move. Cashing in
the 7× headroom means producing the alignment **on the GPU** so there is no reconcile.

## What "the alignment" actually is (the bit-exact target)

`bwt_match_gap` (bwtgap.c) returns, per read, an array `aln[]` of `bwt_aln1_t`:

```c
struct bwt_aln1_t { n_mm:8, n_gapo:8, n_gape:8, score:20, n_ins:10, n_del:10; bwtint_t k, l; };
```

i.e. a list of SA intervals `[k,l]` each tagged with its mismatch/gap counts and score.
`.sai` stores `n_aln` then these entries **in discovery order**. To stay byte-identical to
`bwa aln`, the GPU must produce the *same entries, same order, same count*.

## The roadblock: the search is best-first, and the output is order-sensitive

The existence engine is a **DFS** — correct precisely because existence is *order-free*
(any traversal yields the same yes/no). `bwt_match_gap` is the opposite: a **best-first
priority queue** (`gap_stack_t`: bins keyed by `score = n_mm·3 + n_gapo·11 + n_gape·4`, pop
the lowest-score bin, LIFO within a bin). Five behaviors make its output depend on the exact
pop order:

1. **Best-diff lowering (top-2).** On the *first* hit, `max_diff` drops to `best_diff+1`
   (bwtgap.c:227). This prunes the rest of the search — and is only correct because
   best-first guarantees the first hit found is the best-scoring one.
2. **Early termination.** `if (e.info>>21 > best_score + s_mm) break;` (line 198) and
   `else if (best_cnt > max_top2) break;` (line 230) — both stop the search based on
   best-first pop order and a running best-score count.
3. **`gap_shadow` mutation.** Each accepted hit *mutates the `width[]` bound array*
   (line 238) — changing pruning for every subsequent pop. A sequential data dependency
   *between* hit events.
4. **Discovery-order output.** `aln[]` is appended in pop order, so the array order itself
   is order-sensitive.
5. **`max_entries` cap (2,000,000 queue entries).** On overflow the search aborts and
   returns *partial* results — so even the failure mode must match exactly.

None of these exist in the existence kernel. Reproducing them is the whole job.

## Two implementation routes

### Route 1 — one read per *thread*, faithful best-first (bit-exact by construction)
Each GPU thread runs the exact `bwt_match_gap` loop with its own priority queue in global
memory. Same algorithm → same order → bit-exact, trivially.
- **Pro:** correctness is free; it *is* the reference algorithm, just many reads at once.
- **Con A — divergence returns.** One-read-per-thread is exactly the SIMT heavy-tail stall
  the warp engine was built to avoid: a 2M-entry-queue read freezes its warp while 31 lanes
  idle. Mitigate with the existing budget→flag→CPU-reconcile fallback for the heaviest few %.
- **Con B — memory.** A per-thread priority queue (bins × growable, peak up to `max_entries`)
  is large. Needs the same CAP + spill discipline the existence engine uses, per thread
  rather than per warp.
- **Con C — throughput is unproven.** ~10⁴ slow, divergent GPU threads vs 32 fast CPU cores
  at branchy pointer-chasing PQ ops. Plausibly wins on cheap-dominated mixes, **uncertain on
  the d=5 band** — which is the only band that matters here. *This is the open question.*

### Route 2 — warp-cooperative best-first (research risk, likely not bit-exact)
Keep 32 lanes per read but reproduce best-first order. Behaviors 1–3 above are *sequential
dependencies* (mutate-then-prune); parallelizing the pop while preserving bwa's exact order
appears to require serializing exactly the parts that carry the parallelism. A
"collect-all-then-replay-sequentially" hybrid has to first enumerate the *full* max_diff tree
(more work than bwa does, and on d=5 that's the explosive tree we're trying to avoid). I do
not have a credible bit-exact warp-parallel best-first design; I flag it as unproven.

## Recommendation

Route 1 is the only path I'm confident is *correct*; its *throughput* on the d=5 band is the
real unknown (Con C). Before committing engine work, answer that question in isolation — does
a per-thread best-first collector on the GPU beat a multi-threaded CPU on a heavy-tailed
workload, while staying bit-exact with the sequential reference? That is exactly what the
companion test case (`bestfirst_selftest.cu`) measures, on synthetic data, with no dependence
on this codebase. If it wins there, Route 1 is worth building; if it only ties (as full-engine
gpualn does today), the reconcile floor is fundamental and the cutoff-64 routing stands.
