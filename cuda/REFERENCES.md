# CUDA references consulted (keep checking these as tuning proceeds)

Per user request, re-check CUDA docs and NVIDIA developer forums periodically during optimization.

## Primary docs
- CUDA C++ Programming Guide (13.x): https://docs.nvidia.com/cuda/cuda-c-programming-guide/
- CUDA Ampere GPU Architecture Tuning Guide (sm_86): https://docs.nvidia.com/cuda/ampere-tuning-guide/
  - sm_86: **48 max resident warps/SM**, 64K 32-bit regs/SM, ≤255 regs/thread, 100 KB shared/SM.
  - full 48 warps needs ≤ ~42 regs/thread (we use 69 -> 7 blocks = 28 warps).
- Best Practices Guide / CUDA warp-level primitives: https://developer.nvidia.com/blog/using-cuda-warp-level-primitives/
- L2 persistence, Cooperative Groups, Dynamic Parallelism: in the Programming Guide special-topics.

## Directly applicable findings
- **Volkov, "Understanding Latency Hiding on GPUs"** (UCB EECS-2016-143): latency is hidden by
  ILP *and* occupancy; a dependent-chain kernel (ILP≈1) leans on occupancy. But our measurements
  show raising occupancy *hurts* -> we are stack-memory-capacity bound, not latency bound.
- **"High-Performance N-Queens Solver on GPU: Iterative DFS with Zero Bank Conflicts"**, arXiv
  2511.12009 (2025): iterative DFS with the **stack mapped to shared memory**, bank-conflict-free.
  Target structure for our next big lever.
- NVIDIA/GTC fundamental-optimization decks + forum consensus: **one-thread-per-task tree search
  scales poorly (load imbalance + warp divergence)**; prefer **assigning a subtree to a group of
  threads sharing fast memory**; keep DFS stack in shared/registers (global ~290 cyc vs L1 ~33 cyc).
- Persistent threads >> dynamic parallelism for irregular work (GTC persistent-threads study);
  device-side cudaDeviceSynchronize removed in CUDA 12.

## Open questions to take to the forums next
- Best shared-memory stack window size vs occupancy trade for sm_86 with ~46 bp reads.
- Warp-cooperative FM-index backtracking: split the 4-symbol Occ probe across lanes vs
  one-read-per-lane-but-grouped; any published BWA-on-GPU warp schemes beyond BarraCUDA/NVBIO.
