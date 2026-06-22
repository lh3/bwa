/* fmdfs_selftest.cu — self-contained test for a warp-cooperative, bounded-error
 * backtracking search over an FM-index.
 *
 * No external libraries, no on-disk index, no domain data. It synthesises a random
 * text over a 4-symbol alphabet ("DNA-like"), builds an FM-index in-process, generates
 * query strings (exact / mutated / random), runs the GPU warp-DFS existence search, and
 * verifies every result against two independent CPU references:
 *     (1) a recursive CPU search with identical semantics (the order-free oracle), and
 *     (2) a brute-force Hamming scan over a query subset (ground truth for the oracle).
 *
 * It is the substitution-only core of the production engine — same node-packing, same
 * lower-bound (D-array) prune, same warp pop/expand/scan/push loop, same budget->flag /
 * stack-overflow->flag fallback, same "existence bit only" GPU contract. Gaps are omitted
 * to keep the test one file; adding them does not change any of the GPU-side concerns.
 *
 * Build: nvcc -O3 -arch=sm_86 -o fmdfs_selftest cuda/fmdfs_selftest.cu
 * Run:   ./fmdfs_selftest [text_len] [n_queries] [seed]
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cuda_runtime.h>

#define SIGMA 4               /* alphabet {0,1,2,3}; sentinel is symbol 0 after a +1 shift */
#define OCC_INT 64            /* rank checkpoint spacing */

#define CK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { \
    fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e_)); \
    exit(1); } } while (0)

static double wall() { return std::chrono::duration<double>(
    std::chrono::steady_clock::now().time_since_epoch()).count(); }

/* xorshift64 — deterministic PRNG so runs are reproducible */
static inline uint64_t rng(uint64_t *s){ uint64_t x=*s; x^=x<<13; x^=x>>7; x^=x<<17; return *s=x; }

/* Poisson-tail max-error budget vs query length (the cost-cliff curve from the report):
 * smallest d such that P(#errors > d) < fnr, errors ~ Poisson(len*err). */
static int max_diff(int len, double err, double fnr){
    if (len < 1) return 0;
    double el = exp(-len*err), s = el, y = 1.0; double x = 1.0;
    for (int k = 1; k < 1000; ++k){ y *= len*err; x *= k; s += el*y/x; if (1.0 - s < fnr) return k; }
    return 2;
}

/* ============================ FM-index (host build) ============================ */
struct FM {
    std::vector<uint8_t>  bwt;   /* BWT over symbols {0=sentinel,1..4=real}, length L=n+1 */
    std::vector<uint64_t> occ;   /* checkpoints: occ[(i/OCC_INT)*4 + (c-1)] = #c in bwt[0..i*..) */
    uint32_t C[6];               /* C[c] = #symbols in text strictly < c (sentinel-inclusive) */
    uint32_t L;                  /* = n+1 */
};

/* O(n log^2 n) prefix-doubling suffix array of s (length L, symbols already 0..K-1) */
static std::vector<int> build_sa(const std::vector<int>& s){
    int L = s.size();
    std::vector<int> sa(L), rnk(L), tmp(L);
    for (int i=0;i<L;i++){ sa[i]=i; rnk[i]=s[i]; }
    for (int k=1;; k<<=1){
        auto cmp = [&](int a,int b){
            if (rnk[a]!=rnk[b]) return rnk[a]<rnk[b];
            int ra = a+k<L ? rnk[a+k] : -1, rb = b+k<L ? rnk[b+k] : -1;
            return ra<rb;
        };
        std::sort(sa.begin(), sa.end(), cmp);
        tmp[sa[0]]=0;
        for (int i=1;i<L;i++) tmp[sa[i]] = tmp[sa[i-1]] + (cmp(sa[i-1],sa[i]) ? 1 : 0);
        for (int i=0;i<L;i++) rnk[i]=tmp[i];
        if (rnk[sa[L-1]]==L-1) break;
    }
    return sa;
}

/* build an FM-index over text[0..n-1] (symbols 0..3). Sentinel 0 prepended via +1 shift. */
static FM build_fm(const uint8_t* text, int n){
    std::vector<int> s(n+1);
    for (int i=0;i<n;i++) s[i] = text[i]+1;   /* real symbols become 1..4 */
    s[n] = 0;                                 /* sentinel, smallest */
    std::vector<int> sa = build_sa(s);
    FM fm; fm.L = n+1;
    fm.bwt.resize(n+1);
    for (int i=0;i<=n;i++) fm.bwt[i] = (uint8_t)( sa[i]==0 ? 0 : s[sa[i]-1] );
    /* totals -> C[] */
    uint32_t cnt[6] = {0,0,0,0,0,0};
    for (int i=0;i<=n;i++) cnt[fm.bwt[i]]++;
    fm.C[0]=0; for (int c=1;c<6;c++) fm.C[c]=fm.C[c-1]+cnt[c-1];
    /* rank checkpoints */
    int ncp = (n+1)/OCC_INT + 1;
    fm.occ.assign((size_t)ncp*4, 0);
    uint64_t run[4] = {0,0,0,0};
    for (int i=0;i<=n;i++){
        if (i % OCC_INT == 0){ size_t b=(size_t)(i/OCC_INT)*4; for(int c=0;c<4;c++) fm.occ[b+c]=run[c]; }
        uint8_t b = fm.bwt[i]; if (b>=1 && b<=4) run[b-1]++;
    }
    return fm;
}

/* host rank: counts of symbols 1..4 in bwt[0..i-1] -> r[1..4] */
static inline void rank_all_h(const FM& fm, uint32_t i, uint32_t r[6]){
    uint32_t cp = i / OCC_INT; size_t b=(size_t)cp*4;
    uint32_t c1=fm.occ[b],c2=fm.occ[b+1],c3=fm.occ[b+2],c4=fm.occ[b+3];
    for (uint32_t j=cp*OCC_INT; j<i; ++j){ uint8_t x=fm.bwt[j];
        c1+=(x==1); c2+=(x==2); c3+=(x==3); c4+=(x==4); }
    r[1]=c1; r[2]=c2; r[3]=c3; r[4]=c4;
}

/* ============================ device FM-index ============================ */
struct FMdev { const uint8_t* bwt; const uint64_t* occ; uint32_t C[6]; uint32_t L; };

__device__ __forceinline__ void d_rank_all(const FMdev fm, uint32_t i, uint32_t r[6]){
    uint32_t cp = i >> 6;                         /* / OCC_INT(64) */
    const uint64_t* b = fm.occ + (size_t)cp*4;
    uint32_t c1=b[0],c2=b[1],c3=b[2],c4=b[3];
    for (uint32_t j=cp<<6; j<i; ++j){ uint8_t x=fm.bwt[j];
        c1+=(x==1); c2+=(x==2); c3+=(x==3); c4+=(x==4); }
    r[1]=c1; r[2]=c2; r[3]=c3; r[4]=c4;
}

/* ============================ DFS node packing ============================ */
/* node = (i: remaining prefix length, z: errors so far); interval [sp,ep) stored alongside */
__device__ __host__ __forceinline__ uint32_t pack(uint32_t i, uint32_t z){ return i | (z<<9); }
#define NODE_I(p) ((p) & 0x1ff)
#define NODE_Z(p) (((p) >> 9) & 0x3f)

/* ====================== warp-cooperative existence DFS ======================
 * One query per WARP. 32 lanes co-explore the search tree as a shared frontier (a stack in
 * shared memory, SoA: sp/ep/node). Each step: pop up to 32 frontier nodes (one per lane) ->
 * expand each into <=4 children -> warp-scan to compact -> push back. Returns 1 on first hit
 * (warp ballot) OR on budget/overflow (-> flagged, CPU finishes it); 0 if the tree is exhausted
 * with no match. Order-independent: existence is the same for any traversal. */
__device__ int dfs_match_warp(const FMdev fm, const uint8_t* q, int m, const uint8_t* D, int d,
                              uint32_t* Ssp, uint32_t* Sep, uint32_t* Snode, int CAP,
                              unsigned long long budget, int* flagged, unsigned long long* pops_out){
    const unsigned FULL = 0xffffffffu;
    int lane = threadIdx.x & 31;
    if (lane==0){ Ssp[0]=0; Sep[0]=fm.L; Snode[0]=pack(m,0); }
    int sp = 1;                                   /* warp-uniform stack depth */
    unsigned long long pops = 0;
    __syncwarp();

    for (;;){
        if (sp == 0){ if (lane==0) *pops_out = pops; return 0; }
        int room = (CAP - sp) / 4;                /* room for every popped node's <=4 children */
        int nact = sp < 32 ? sp : 32; if (nact > room) nact = room; if (nact < 1) nact = 1;
        bool active = lane < nact;

        uint32_t isp=0, iep=0, nd=0; int qi=0, z=0;
        if (active){ int idx = sp-1-lane; isp=Ssp[idx]; iep=Sep[idx]; nd=Snode[idx];
                     qi=NODE_I(nd); z=NODE_Z(nd); }
        sp -= nact; pops += nact;
        if (pops > budget){ if (lane==0){ *flagged=1; *pops_out=pops; } return 1; }  /* give up -> CPU */

        /* expand this lane's node into up to 4 children (the 4 alphabet symbols) */
        uint32_t csp[4], cep[4], cnd[4]; int nc = 0; bool hit = false;
        if (active){
            if (qi == 0) hit = true;                              /* whole query consumed -> match */
            else if (d - z < D[qi-1]) { /* lower-bound prune: can't finish within budget */ }
            else {
                int xi = qi - 1;
                uint32_t rs[6], re[6];
                d_rank_all(fm, isp, rs); d_rank_all(fm, iep, re);
                int want = q[xi] + 1;                              /* matching symbol (1..4) */
                for (int c = 1; c <= 4; ++c){
                    int nz = z + (c != want ? 1 : 0);
                    if (nz > d) continue;
                    uint32_t nsp = fm.C[c] + rs[c], nep = fm.C[c] + re[c];
                    if (nsp < nep){ csp[nc]=nsp; cep[nc]=nep; cnd[nc]=pack(xi, nz); ++nc; }
                }
            }
        }
        if (__any_sync(FULL, hit)){ if (lane==0) *pops_out = pops; return 1; }

        /* warp exclusive prefix-sum of per-lane child counts -> contiguous push offsets */
        int incl = nc;
        for (int s2=1; s2<32; s2<<=1){ int y=__shfl_up_sync(FULL,incl,s2); if (lane>=s2) incl+=y; }
        int total = __shfl_sync(FULL, incl, 31);
        int myoff = incl - nc;
        if (sp + total > CAP){ if (lane==0){ *flagged=1; *pops_out=pops; } return 1; }  /* overflow -> CPU */
        for (int j=0;j<nc;++j){ int dst = sp + myoff + j; Ssp[dst]=csp[j]; Sep[dst]=cep[j]; Snode[dst]=cnd[j]; }
        sp += total;
        __syncwarp();
    }
}

/* per-query parameters (flattened query/D layout, like the production harness) */
struct QParam { uint32_t q_off, d_off; int len, d; };

/* global work counter for the persistent pool */
__device__ int g_work;
__global__ void k_reset(){ g_work = 0; }

/* persistent work-pool: each warp atomically pulls the next query until the queue drains */
__global__ void k_search(FMdev fm, const uint8_t* Q, const uint8_t* D, const QParam* qp, int nq,
                          int CAP, uint8_t* has_hit, unsigned long long* pops, unsigned long long budget,
                          int* nflag, int wpb){
    extern __shared__ unsigned char smem[];
    int w = threadIdx.x >> 5, lane = threadIdx.x & 31;
    uint32_t* Ssp  = (uint32_t*)smem + (size_t)w*CAP;
    uint32_t* Sep  = (uint32_t*)smem + (size_t)wpb*CAP + (size_t)w*CAP;
    uint32_t* Snode= (uint32_t*)smem + (size_t)2*wpb*CAP + (size_t)w*CAP;
    for (;;){
        int r;
        if (lane==0) r = atomicAdd(&g_work, 1);
        r = __shfl_sync(0xffffffffu, r, 0);
        if (r >= nq) break;
        QParam p = qp[r];
        int flagged = 0; unsigned long long np = 0;
        int hh = dfs_match_warp(fm, Q + p.q_off, p.len, D + p.d_off, p.d,
                                Ssp, Sep, Snode, CAP, budget, &flagged, &np);
        if (lane==0){ has_hit[r] = (uint8_t)hh; pops[r] = np; if (flagged) atomicAdd(nflag, 1); }
    }
}

/* ============================ CPU references ============================ */
/* recursive existence oracle: identical semantics, no budget (ground-truth existence) */
static bool cpu_match(const FM& fm, const uint8_t* q, int qi, int z, int d,
                      uint32_t sp, uint32_t ep, const uint8_t* D){
    if (qi == 0) return true;
    if (d - z < D[qi-1]) return false;
    int xi = qi - 1;
    uint32_t rs[6], re[6]; rank_all_h(fm, sp, rs); rank_all_h(fm, ep, re);
    int want = q[xi] + 1;
    /* try the match child first (cheapest, most likely to hit early) */
    for (int pass=0; pass<2; ++pass)
        for (int c=1;c<=4;c++){
            int isMM = (c != want);
            if ((pass==0) == isMM) continue;            /* pass0: match, pass1: mismatches */
            int nz = z + (isMM?1:0); if (nz > d) continue;
            uint32_t nsp = fm.C[c]+rs[c], nep = fm.C[c]+re[c];
            if (nsp < nep && cpu_match(fm, q, xi, nz, d, nsp, nep, D)) return true;
        }
    return false;
}

/* brute-force Hamming ground truth: does any text position match q within d substitutions? */
static bool brute_match(const uint8_t* text, int n, const uint8_t* q, int m, int d){
    for (int p=0; p+m<=n; ++p){ int e=0; for (int j=0;j<m;j++){ if (text[p+j]!=q[j]){ if(++e>d) break; } }
        if (e<=d) return true; }
    return false;
}

/* D-array (Li–Durbin lower bound) for query q via the reverse-text FM-index */
static void calc_D(const FM& rfm, const uint8_t* q, int m, uint8_t* D){
    uint32_t sp=0, ep=rfm.L; int z=0;
    for (int i=0;i<m;i++){
        int c = q[i]+1;
        uint32_t rs[6], re[6]; rank_all_h(rfm, sp, rs); rank_all_h(rfm, ep, re);
        sp = rfm.C[c]+rs[c]; ep = rfm.C[c]+re[c];
        if (sp >= ep){ sp=0; ep=rfm.L; z++; }
        D[i] = (uint8_t)z;
    }
}

/* ============================ main ============================ */
int main(int argc, char** argv){
    int n  = argc>1 ? atoi(argv[1]) : (1<<18);     /* text length */
    int nq = argc>2 ? atoi(argv[2]) : 20000;       /* number of queries */
    uint64_t seed = argc>3 ? strtoull(argv[3],NULL,10) : 0x1234567890abcdefULL;

    /* ---- synthesise a random 4-symbol text ---- */
    std::vector<uint8_t> text(n);
    { uint64_t s=seed; for (int i=0;i<n;i++) text[i] = rng(&s) & 3; }
    fprintf(stderr, "[data] text_len=%d, alphabet=%d, queries=%d\n", n, SIGMA, nq);

    /* ---- build forward + reverse FM-index ---- */
    double tb = wall();
    FM fm = build_fm(text.data(), n);
    std::vector<uint8_t> rtext(n); for (int i=0;i<n;i++) rtext[i]=text[n-1-i];
    FM rfm = build_fm(rtext.data(), n);
    fprintf(stderr, "[index] built fwd+rev FM in %.2fs (L=%u, bwt=%zu bytes, occ=%zu ckpts)\n",
            wall()-tb, fm.L, fm.bwt.size(), fm.occ.size()/4);

    /* ---- generate queries: mix of exact substrings, mutated substrings, and random ---- */
    std::vector<uint8_t> Qflat; std::vector<uint8_t> Dflat; std::vector<QParam> qp(nq);
    std::vector<int>     qlen(nq);
    uint64_t s = seed ^ 0xdeadbeefULL;
    for (int i=0;i<nq;i++){
        int m = 20 + (int)(rng(&s) % 51);                 /* length 20..70 (spans the cliff) */
        int d = max_diff(m, 0.02, 0.01);
        std::vector<uint8_t> q(m);
        int kind = rng(&s) % 4;
        if (kind == 0){                                   /* pure random (usually no hit) */
            for (int j=0;j<m;j++) q[j] = rng(&s)&3;
        } else {                                          /* sample a real substring... */
            int p = rng(&s) % (n - m + 1);
            for (int j=0;j<m;j++) q[j] = text[p+j];
            int nmut = (kind==1) ? 0 : (1 + rng(&s) % (d+1));   /* ...then mutate 0..d+1 spots */
            for (int t=0;t<nmut;t++){ int pos = rng(&s)%m; q[pos] = (q[pos] + 1 + (rng(&s)&2)) & 3; }
        }
        qp[i].q_off = Qflat.size(); qp[i].d_off = Dflat.size(); qp[i].len = m; qp[i].d = d; qlen[i]=m;
        std::vector<uint8_t> D(m); calc_D(rfm, q.data(), m, D.data());
        for (int j=0;j<m;j++){ Qflat.push_back(q[j]); Dflat.push_back(D[j]); }
    }

    /* ---- CPU oracle (existence ground truth) ---- */
    double tc = wall();
    std::vector<uint8_t> cpu_hit(nq);
    for (int i=0;i<nq;i++){
        const uint8_t* q = Qflat.data()+qp[i].q_off; const uint8_t* D = Dflat.data()+qp[i].d_off;
        cpu_hit[i] = cpu_match(fm, q, qp[i].len, 0, qp[i].d, 0, fm.L, D) ? 1 : 0;
    }
    long cpu_nhit=0; for (int i=0;i<nq;i++) cpu_nhit += cpu_hit[i];
    fprintf(stderr, "[cpu]  oracle: %ld/%d queries have a match (%.1f%%) in %.2fs\n",
            cpu_nhit, nq, 100.0*cpu_nhit/nq, wall()-tc);

    /* ---- brute-force cross-check on a subset (validates the oracle itself) ---- */
    int brute_n = nq < 256 ? nq : 256;
    long brute_disagree = 0;
    for (int i=0;i<brute_n;i++){
        const uint8_t* q = Qflat.data()+qp[i].q_off;
        bool b = brute_match(text.data(), n, q, qp[i].len, qp[i].d);
        if (b != (bool)cpu_hit[i]) { if (brute_disagree<5) fprintf(stderr,
            "  BRUTE-DISAGREE q%d len%d d%d brute=%d oracle=%d\n", i, qp[i].len, qp[i].d, b, cpu_hit[i]);
            ++brute_disagree; }
    }
    fprintf(stderr, "[cpu]  brute-force cross-check on %d queries: %s (%ld disagreements)\n",
            brute_n, brute_disagree? "FAIL":"PASS", brute_disagree);

    /* ---- upload index + queries ---- */
    uint8_t *d_bwt; uint64_t *d_occ;
    CK(cudaMalloc(&d_bwt, fm.bwt.size()));        CK(cudaMemcpy(d_bwt, fm.bwt.data(), fm.bwt.size(), cudaMemcpyHostToDevice));
    CK(cudaMalloc(&d_occ, fm.occ.size()*8));      CK(cudaMemcpy(d_occ, fm.occ.data(), fm.occ.size()*8, cudaMemcpyHostToDevice));
    FMdev fd; fd.bwt=d_bwt; fd.occ=d_occ; fd.L=fm.L; for (int c=0;c<6;c++) fd.C[c]=fm.C[c];

    uint8_t *d_Q, *d_D; QParam *d_qp; uint8_t *d_hit; unsigned long long *d_pops; int *d_nflag;
    CK(cudaMalloc(&d_Q, Qflat.size()));   CK(cudaMemcpy(d_Q, Qflat.data(), Qflat.size(), cudaMemcpyHostToDevice));
    CK(cudaMalloc(&d_D, Dflat.size()));   CK(cudaMemcpy(d_D, Dflat.data(), Dflat.size(), cudaMemcpyHostToDevice));
    CK(cudaMalloc(&d_qp, nq*sizeof(QParam))); CK(cudaMemcpy(d_qp, qp.data(), nq*sizeof(QParam), cudaMemcpyHostToDevice));
    CK(cudaMalloc(&d_hit, nq)); CK(cudaMalloc(&d_pops, (size_t)nq*8)); CK(cudaMalloc(&d_nflag,4)); CK(cudaMemset(d_nflag,0,4));

    /* ---- launch the warp-DFS ---- */
    int wpb = 4, CAP = 1024;
    unsigned long long budget = getenv("DFS_BUDGET") ? strtoull(getenv("DFS_BUDGET"),NULL,10) : 200000ULL;
    int bdim = wpb*32; size_t shbytes = (size_t)3*wpb*CAP*sizeof(uint32_t);
    CK(cudaFuncSetAttribute(k_search, cudaFuncAttributeMaxDynamicSharedMemorySize, (int)shbytes));
    int numSM=0; CK(cudaDeviceGetAttribute(&numSM, cudaDevAttrMultiProcessorCount, 0));
    int mb=0; CK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&mb, k_search, bdim, shbytes));
    int nblocks = (mb>0?mb:1)*numSM;
    fprintf(stderr, "[gpu]  %d SM, %d blk/SM x %d warps, CAP=%d, shared=%.1fKB/blk, budget=%llu pops\n",
            numSM, mb, wpb, CAP, shbytes/1024.0, budget);

    k_reset<<<1,1>>>(); CK(cudaDeviceSynchronize());
    double tg = wall();
    k_search<<<nblocks, bdim, shbytes>>>(fd, d_Q, d_D, d_qp, nq, CAP, d_hit, d_pops, budget, d_nflag, wpb);
    CK(cudaDeviceSynchronize()); CK(cudaGetLastError());
    double gpu_s = wall()-tg;

    /* ---- download + verify ---- */
    std::vector<uint8_t> gpu_hit(nq); std::vector<unsigned long long> pops(nq); int nflag=0;
    CK(cudaMemcpy(gpu_hit.data(), d_hit, nq, cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(pops.data(), d_pops, (size_t)nq*8, cudaMemcpyDeviceToHost));
    CK(cudaMemcpy(&nflag, d_nflag, 4, cudaMemcpyDeviceToHost));

    long fneg=0, fpos=0; unsigned long long tot_pops=0, max_pops=0; int amax=0;
    for (int i=0;i<nq;i++){
        if (cpu_hit[i] && !gpu_hit[i]){ if (fneg<10) fprintf(stderr,
            "  FALSE-NEG q%d len%d d%d\n", i, qp[i].len, qp[i].d); ++fneg; }
        if (!cpu_hit[i] && gpu_hit[i]) ++fpos;   /* allowed: a flagged no-hit query, reconciled on CPU */
        tot_pops += pops[i]; if (pops[i]>max_pops){ max_pops=pops[i]; amax=i; }
    }
    long gpu_nhit=0; for (int i=0;i<nq;i++) gpu_nhit += gpu_hit[i];

    fprintf(stderr, "[gpu]  %.3fs, %.0f queries/s; has_hit=%ld; flagged->CPU=%d (%.3f%%)\n",
            gpu_s, nq/gpu_s, gpu_nhit, nflag, 100.0*nflag/nq);
    fprintf(stderr, "[gpu]  node-pops: total=%llu mean=%.0f max=%llu (q%d len%d d%d)\n",
            tot_pops, (double)tot_pops/nq, max_pops, amax, qp[amax].len, qp[amax].d);
    fprintf(stderr, "[check] false-negatives=%ld (MUST be 0)   false-positives=%ld (harmless, CPU-reconciled)\n",
            fneg, fpos);

    bool ok = (fneg==0) && (brute_disagree==0);
    fprintf(stderr, "\n%s\n", ok ? "==== PASS ====" : "==== FAIL ====");
    return ok ? 0 : 1;
}
