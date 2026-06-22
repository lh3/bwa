# Porting an unbounded-budget backtracking search to CUDA

*A self-contained engineering report. Domain-neutral: this is purely about running a
bounded-error backtracking string search on an FM-index, and what breaks when you put
that search on a GPU.*

---

## 1. The computational problem

We have a large immutable text `T` (hundreds of millions to billions of symbols over a
4-letter alphabet `Σ = {0,1,2,3}`) preprocessed into an **FM-index** (a BWT plus rank
checkpoints). We are given a stream of short **query strings** `Q` (tens of symbols
each, millions of them) and must, for each query, find the set of positions in `T` that
match it **within an error budget** `d` (substitutions + a small number of short gaps).

### 1.1 The primitive: backward extension

The FM-index represents any matched pattern as a **suffix-array interval** `[k, l]`
(half-open count `l - k + 1` = number of occurrences). The one operation it supports is
*backward extension*: given the interval for pattern `P`, the interval for `cP` (prepend
symbol `c`) is computed in O(1) from two **rank queries** (`Occ`):

```
k' = C[c] + Occ(c, k-1) + 1
l' = C[c] + Occ(c, l)
```

`C[c]` is a length-Σ prefix-sum table; `Occ(c, i)` = number of `c`'s in `BWT[0..i]`,
answered from the nearest checkpoint plus a short popcount-style scan. **Exact** search
of a length-`m` query is therefore `m` rank queries — trivial.

### 1.2 Approximate search = a DFS over an implicit tree

Allowing up to `d` errors turns the single backward path into a **branch-and-bound DFS**
over an implicit tree. At each query position, from the current interval you may:

- **match** the query symbol (1 child, no cost),
- **substitute** (3 children — the other symbols — each costs 1 mismatch),
- **open / extend an insertion or deletion** (gap children, cost 1 each, capped).

Every edge is one O(1) rank query that narrows `[k, l]`. A node is a tuple
`(i, k, l, n_mm, n_gapo, n_gape, state)` where `i` is the remaining-prefix length.

The search is made tractable by a **bound function** computed once per query: an array
`bid[i]` = the *minimum* number of errors any match of the prefix `Q[0..i]` must incur.
Before expanding a node you check `remaining_budget >= bid[i]`; if not, prune the whole
subtree. This is the admissible heuristic that keeps the tree from being `4^m`.

The number of surviving nodes is approximately

```
nodes(m, d)  ≈  C(m, d) · (Σ-1)^d   =   C(m, d) · 3^d
```

— **exponential in the error budget `d`**.

---

## 2. What `-l 1024` does (the seeding switch)

Normally the search is **seeded**: only a short prefix of each query (the *seed*, e.g.
the first `s = 32` symbols) is searched with the error-tolerant DFS, and only a tiny
budget is allowed inside the seed; the rest of the query is required to extend with no
extra freedom. Seeding is the safety valve — it caps `d` *inside the explosive region*
to a small constant regardless of how long the query is.

`-l 1024` sets the seed length to **1024**, larger than any query in the stream. This
**disables seeding entirely**: the *whole* query is now searched with the full error
budget, and the budget itself is a function of query length,

```
d = maxdiff(m)        // grows with m; a step function in practice
```

So two things happen at once as queries get longer: `m` grows **and** `d` grows. Plugging
into `nodes ≈ C(m, d)·3^d`, the per-query cost doesn't grow smoothly — it has **cliffs**
at each `d` increment. Measured `maxdiff` thresholds (for our error model) land at:

| query length `m` | error budget `d` |
|---|---|
| ≤ 21  | 2 |
| 22–41 | 3 |
| 42–63 | 4 |
| 64–89 | 5 |
| ≥ 90  | 6 |

The jump from `d=4` to `d=5` at **`m = 64`** is the killer: `C(64,5)·3^5` is ~2 orders of
magnitude above `C(63,4)·3^4`. Below the cliff the average query expands tens to hundreds
of nodes; above it, a non-trivial fraction expand **millions**.

**The defining property of the workload is therefore an extremely heavy-tailed,
*unbounded* and *unpredictable* per-query cost.** ~99.98% of queries are cheap; a
fraction of a percent are catastrophic; and you cannot tell which is which without
running the search.

---

## 3. Why this is hard on a GPU specifically

A GPU is a throughput machine built for **regular, balanced, statically-schedulable**
work. Every property of §2 is the opposite.

1. **SIMT lockstep → warp divergence.** 32 lanes of a warp share one program counter.
   The obvious mapping — *one query per thread, 32 queries per warp* — runs each warp
   until its **slowest lane** is done. With a heavy-tailed cost distribution, one
   exploding query freezes 31 idle lanes for millions of steps. Utilization collapses to
   ~1/32 exactly on the queries that dominate runtime. This is the central problem.

2. **No recursion, tiny per-thread memory.** DFS wants a call stack. GPU threads have no
   growable stack and only a few KB of fast state. A per-thread explicit stack sized for
   the worst-case frontier (millions of nodes) × thousands of resident threads cannot
   fit anywhere.

3. **Irregular, uncoalesced memory.** Every tree edge is a rank query: a near-random
   gather into an index far larger than L2. Lanes in a warp touch unrelated checkpoints →
   uncoalesced loads and cache thrash. There is no natural spatial locality to exploit
   across lanes.

4. **No static partition possible.** Because cost is unknown until you run it, you cannot
   bucket queries into balanced batches ahead of time. Any fixed assignment of queries to
   threads is wrong for the tail.

5. **Bit-exactness constraint.** The output must be byte-identical to a trusted scalar
   reference. So we may not approximate, may not cap quality, and may not reorder anything
   in a way that changes the result — which rules out the usual GPU trick of "just bound
   the work and accept slightly worse answers."

---

## 4. The solution

Six ideas compose into a hybrid GPU+CPU engine that is fast **and** exact.

### 4.1 One query per *warp*, not per *thread* (kills divergence)

The whole warp cooperates on a **single** query's tree. The frontier (stack of
unexpanded nodes) is shared by the warp. Each iteration, up to 32 lanes each **pop one
node and expand it in parallel**, producing 0–9 children each; a warp-scan compacts the
variable child counts contiguously back onto the shared stack.

Intra-tree parallelism becomes lane parallelism. There is no longer any divergence from
*unequal trees* — a warp only ever works one tree at a time — and the expensive query
gets 32 hands instead of 1.

```c
// node packed into a u32: i(0-8) | n_mm(9-14) | n_gapo(15-18) | n_gape(19-22) | state(23-24)
__device__ __forceinline__ uint32_t pack_node(int i,int mm,int go,int ge,int st)
{ return (uint32_t)i | (mm<<9) | (go<<15) | (ge<<19) | (st<<23); }

// --- core warp loop (abridged) ---
int sp = 1;                                  // shared-stack depth (the frontier)
for (;;) {
    // 1. pop up to 32 nodes, one per lane, bounded by stack room for children
    int n_active = min(sp, 32);
    int r9 = room / 9; if (n_active > r9) n_active = r9;   // 9 = max children per node
    bool active = lane < n_active;
    if (active) { int idx = sp-1-lane; load node idx into registers; }
    sp -= n_active;
    nn += n_active;                          // node-expansion counter (budget, §4.3)

    // 2. each active lane expands its node -> up to 9 children in registers (cck/ccl/ccn)
    int nc = 0; bool hit = false;
    if (active) expand(node, /*out*/ children, &nc, &hit);

    // 3. existence short-circuit: ANY lane found an acceptable match -> done (§4.4)
    if (__any_sync(FULL, hit)) return 1;

    // 4. warp-scan the per-lane child counts, then scatter children contiguously
    int incl = nc;
    for (int d=1; d<32; d<<=1) { int y=__shfl_up_sync(FULL,incl,d); if (lane>=d) incl+=y; }
    int total = __shfl_sync(FULL, incl, 31);
    int myoff = incl - nc;
    for (int j=0;j<nc;++j){ int d=sp+myoff+j; push child j at slot d; }
    sp += total;
    __syncwarp();
    if (sp == 0) break;                      // (plus the global-backing refill, §4.2)
}
```

`expand()` is the branch-and-bound body: it applies the `bid[]` prune, emits the match /
3 substitution / gap children via rank queries, and sets `hit` when `i` reaches 0 within
budget.

### 4.2 Two-level stack: shared "top window" + global backing (keeps bushy queries resident)

The hot top of the frontier lives in **shared memory** (`CAP_SM` nodes per warp — fast,
but small). When a bushy query overflows it, a fixed `CHUNK` of the *bottom* of the
window is **spilled** to a per-warp slab in **global memory** (`CAP_GL` nodes); when the
shared window later drains, a chunk is **pulled back**. Shared-memory use stays a small
constant per warp, yet a query with a huge frontier stays resident on the GPU instead of
failing.

```c
int room = CAP_SM - sp;
if (room < 9) {                              // not enough room for one node's children
    if (gsp + CHUNK > CAP_GL) { flag_to_cpu(); return 1; }      // even the slab is full
    cooperatively copy sk/sl/sn[0..CHUNK) -> global slab[gsp..]; // spill the window bottom
    shift the remaining shared window down by CHUNK;
    gsp += CHUNK; sp -= CHUNK;
}
// ... and on underflow:
if (sp == 0) {
    if (gsp == 0) return 0;                   // truly empty -> no hit
    int c = min(gsp, CHUNK);
    cooperatively copy global slab top chunk -> shared window;   // refill
    gsp -= c; sp = c;
}
```

### 4.3 Hard budget cap + exact CPU reconcile (bounds the tail)

Every query carries a hard **node-expansion budget** (`DFS_BUDGET`, default 2,000,000).
If a query exceeds it — *or* its frontier would overflow the global slab — the warp
**aborts it and sets a `flagged` bit** instead of finishing it:

```c
sp -= n_active; nn += n_active;
if (nn > budget) { *flagged = 1; return 1; }   // give up on the GPU, hand to CPU
```

The aborted set is tiny and is finished by a **multithreaded CPU pass running the exact
reference search** on just those queries. Measured flag rate: **~0.015%**. The GPU never
pays the catastrophic tail; the CPU absorbs it on a few thousand queries out of tens of
millions. Throughput stays high *and* the answer stays exact.

### 4.4 Existence-only detection → order-independence → exactness (the keystone)

The crucial realization: **the GPU does not need to reproduce the reference's traversal
order or its full result.** It only needs to decide, per query, a single bit:

> *"Does this query have at least one acceptable match, reachable within budget?"*

Existence is **order-independent** — any traversal of the tree yields the same yes/no, so
the warp engine is free to expand nodes in whatever order is convenient for the hardware.
The actual *exact, order-sensitive* result (best hits, tie-breaking, secondary records)
is regenerated **on the CPU, in input order**, only for the queries that the GPU marked
as hit-or-flagged. That clean split — fast order-free filter on the GPU, exact ordered
finish on the CPU — is what makes the hybrid **bit-identical** to the pure-CPU reference.

```c
// the kernel's per-query output is just one byte:
has_hit[r] = hh;          // 1 = has a match (or was flagged); 0 = provably no match
// downstream: only reads with has_hit[r]==1 are re-searched exactly on the CPU,
//             in order, by a single ordered consumer (preserves any RNG/output state).
```

### 4.5 Persistent work-pool kernel (dynamic load balance, no host round-trips)

We launch a **fixed** grid sized to fully occupy the device (from the occupancy
calculator), then let every warp **pull the next query itself** via a single global
`atomicAdd` counter, looping until the queue drains. A warp that finishes a cheap query
in microseconds immediately grabs another; a warp stuck on a bushy query just holds its
slot. This is *self-balancing* and needs no host involvement — directly answering the
"can't statically partition" problem (§3.4).

```c
__global__ void k_dfs_warp2(/* index, queries, params, budget, outputs */) {
    // each warp owns one shared top-window + one global slab (by blockIdx/warp id)
    for (;;) {
        int r;
        if (lane == 0) r = atomicAdd(workctr, 1);     // claim the next query
        r = __shfl_sync(FULL, r, 0);
        if (r >= nreads) break;
        int flagged = 0; unsigned long long nn = 0;
        int hh = dfs_has_hit_warp2(/* query r */, budget, &flagged, &nn);
        if (lane == 0) { has_hit[r] = hh; if (flagged) atomicAdd(nflag, 1); }
    }
}
```

### 4.6 CPU/GPU pipeline overlap (hide the reconcile + I/O)

The stream is processed in fixed-size chunks. The **main thread owns the GPU** end to end
for a chunk — preprocess → upload → kernel → download the `has_hit` bytes — kept serial
on **one** stream so the shared global slabs are never raced between concurrent kernels.
Each finished chunk is handed to a **single in-order finisher thread** that runs the CPU
reconcile + output *concurrently with the next chunk's GPU work*. Because there is exactly
one ordered consumer, output order (and any RNG state threaded through it) is preserved —
again, exactness for free.

```
main thread:   [preprocess|upload|kernel|download]   [preprocess|upload|kernel|download] ...
finisher:                        [reconcile+emit chunk N-1] [reconcile+emit chunk N]   ...
                                  ^ single ordered consumer  => deterministic, bit-exact
```

### 4.7 (Upstream) cost-model partitioning — don't send the cliff to the GPU at all

Because `nodes ≈ C(m,d)·3^d` with `d = maxdiff(m)` is *predictable from length alone*,
the input can be cheaply pre-split at the `d=4→5` cliff (`m = 64`): queries below the
cliff go to this GPU engine; the genuinely explosive long tail is routed to a different
algorithm entirely. The split point is auto-chosen by walking the length histogram and
accumulating `Σ count(m)·cost(m)` until a runtime budget is hit, hard-capped at the
cliff. This keeps the GPU fed with work it is good at.

```python
def maxdiff(m, fnr, err=0.02):                  # smallest d with P(>d errors) < fnr
    el = math.exp(-m*err); s = el; y = 1.0; x = 1
    for k in range(1, 1000):
        y *= m*err; x *= k; s += el*y/x
        if 1.0 - s < fnr: return k
    return 2

def cost(m):                                    # tree-cost proxy per query of length m
    d = maxdiff(m, fnr); return comb(m, d) * 3**d if m >= d else 1

cliff = next(m for m in range(1, maxlen+2) if maxdiff(m, fnr) > 4)   # d=4->5 boundary
cum = 0; best = 1
for L in range(1, cliff+1):                     # largest cutoff under the runtime budget
    cum += cost(L-1) * hist.get(L-1, 0)
    best = L if cum <= budget else best
    if cum > budget: break
```

---

## 5. Results

On a single 24 GB Ampere-class GPU with a 32-thread host:

- **~17,000 queries/s sustained**, ≈ **5×** a fully multithreaded CPU baseline.
- **~0.015%** of queries flagged to the CPU reconcile path.
- **Byte-identical** output to the scalar reference (verified by direct comparison).
- Shared-memory footprint bounded to a small constant per warp regardless of query
  bushiness; the global slabs are a one-time fixed allocation sized by occupancy.

## 6. Transferable lessons

1. **Match the parallelism granularity to the work-item variance.** When per-item cost is
   heavy-tailed, *one item per thread* is a trap on SIMT hardware. Drop to **one item per
   warp** and parallelize *inside* the item so a warp's lanes always share fate.
2. **Bound the unbounded, then repair exactly elsewhere.** A hard budget on the GPU plus a
   tiny exact CPU reconcile turns an unbounded worst case into a bounded common case
   without sacrificing correctness.
3. **Weaken the GPU's contract to the part that is order-independent.** Computing an
   *existence bit* (order-free) on the GPU and regenerating the *ordered* result on the
   CPU is what makes a parallel reordered engine bit-exact with a scalar reference.
4. **Two-level stacks** (fast shared top window + global backing, with chunked
   spill/refill) let irregular DFS stay resident on-device without sizing shared memory
   for the worst case.
5. **Persistent kernels + an atomic work counter** beat static partitioning whenever
   per-item cost is unknowable in advance.
6. **A closed-form cost model is leverage.** If you can predict cost from a cheap feature
   (here, length → `C(m,d)·3^d`), you can route work to the right engine *before* paying
   for it.
