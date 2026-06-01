/* Phase 1 validation harness for the device FM-index.
 *
 * Loads a real bwa .bwt index, uploads bwt->bwt + cnt_table + L2 to the GPU, then:
 *   Test A: Occ4 at millions of random positions, compared bit-for-bit to CPU bwt_occ4.
 *   Test B: full backward exact-search of real reads, compared to CPU bwt_match_exact.
 * Reports mismatches (must be zero) and a first Occ-throughput datapoint.
 *
 * Build:  make fmtest      Run:  ./fmtest <ref.fa> [reads.fq] [n_random_occ]
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <cuda_runtime.h>

extern "C" {
#include "bwt.h"
}

#define FM_DEVICE_DEFINE_CONST
#include "fm_device.cuh"

#define CK(call) do { cudaError_t e_ = (call); if (e_ != cudaSuccess) { \
	fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(e_)); \
	exit(1); } } while (0)

/* A/a->0 C/c->1 G/g->2 T/t->3 else 4 */
static unsigned char nt4[256];
static void init_nt4(void){ for(int i=0;i<256;i++) nt4[i]=4;
	nt4['A']=nt4['a']=0; nt4['C']=nt4['c']=1; nt4['G']=nt4['g']=2; nt4['T']=nt4['t']=3; }

static inline uint64_t xorshift64(uint64_t *s){ uint64_t x=*s; x^=x<<13; x^=x>>7; x^=x<<17; return *s=x; }

/* ---- kernels ---- */
__global__ void k_occ4(const uint32_t *bwt, uint64_t primary,
                       const uint64_t *ks, uint64_t n, uint64_t *out /* n*4 */)
{
	uint64_t i = blockIdx.x * (uint64_t)blockDim.x + threadIdx.x;
	if (i >= n) return;
	uint64_t cnt[4];
	d_bwt_occ4(bwt, primary, ks[i], cnt);
	out[i*4+0]=cnt[0]; out[i*4+1]=cnt[1]; out[i*4+2]=cnt[2]; out[i*4+3]=cnt[3];
}

__global__ void k_match(fmidx_dev fm, const uint8_t *seqs, const uint64_t *off,
                        const int *len, int nreads, uint64_t *outk, uint64_t *outl, int *outc)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= nreads) return;
	uint64_t k=0,l=0; int c;
	c = d_bwt_match_exact(fm, seqs + off[i], len[i], &k, &l);
	outk[i]=k; outl[i]=l; outc[i]=c;
}

int main(int argc, char **argv)
{
	if (argc < 2) { fprintf(stderr,"usage: %s <ref.fa> [reads.fq] [n_random_occ]\n", argv[0]); return 1; }
	init_nt4();
	const char *prefix = argv[1];
	const char *reads_fn = argc > 2 ? argv[2] : NULL;
	uint64_t n_occ = argc > 3 ? strtoull(argv[3],NULL,10) : 5000000ULL;

	char bwt_fn[4096]; snprintf(bwt_fn, sizeof bwt_fn, "%s.bwt", prefix);
	fprintf(stderr, "[load] %s\n", bwt_fn);
	bwt_t *bwt = bwt_restore_bwt(bwt_fn);
	if (!bwt) { fprintf(stderr,"failed to load bwt\n"); return 1; }
	fprintf(stderr, "[load] seq_len=%llu primary=%llu bwt_size=%llu (%.2f GB)\n",
		(unsigned long long)bwt->seq_len, (unsigned long long)bwt->primary,
		(unsigned long long)bwt->bwt_size, bwt->bwt_size*4.0/1e9);

	/* upload index */
	uint32_t *d_bwt=NULL;
	CK(cudaMalloc(&d_bwt, bwt->bwt_size * sizeof(uint32_t)));
	CK(cudaMemcpy(d_bwt, bwt->bwt, bwt->bwt_size * sizeof(uint32_t), cudaMemcpyHostToDevice));
	CK(cudaMemcpyToSymbol(c_cnt_table, bwt->cnt_table, sizeof(uint32_t)*256));
	CK(cudaMemcpyToSymbol(c_L2, bwt->L2, sizeof(uint64_t)*5));
	fmidx_dev fm{ d_bwt, bwt->primary, bwt->seq_len };

	/* ---------- Test A: random Occ4 ---------- */
	fprintf(stderr, "[A] random Occ4: %llu probes\n", (unsigned long long)n_occ);
	std::vector<uint64_t> hk(n_occ);
	uint64_t seed = 0x9e3779b97f4a7c15ULL;
	/* FMTEST_KRANGE caps the probed k-range to test L1/L2/HBM working-set effects (locality ceiling) */
	uint64_t krange = getenv("FMTEST_KRANGE") ? strtoull(getenv("FMTEST_KRANGE"),NULL,10) : (bwt->seq_len + 1);
	if (krange > bwt->seq_len + 1) krange = bwt->seq_len + 1;
	fprintf(stderr, "[A] k-range = %llu (%.1f MB of BWT touched)\n", (unsigned long long)krange, krange/4.0/1e6);
	for (uint64_t i=0;i<n_occ;i++) hk[i] = xorshift64(&seed) % krange;
	uint64_t *d_ks=NULL, *d_out=NULL;
	CK(cudaMalloc(&d_ks, n_occ*sizeof(uint64_t)));
	CK(cudaMalloc(&d_out, n_occ*4*sizeof(uint64_t)));
	CK(cudaMemcpy(d_ks, hk.data(), n_occ*sizeof(uint64_t), cudaMemcpyHostToDevice));
	int tpb=128; uint64_t blocks=(n_occ+tpb-1)/tpb;
	cudaEvent_t t0,t1; CK(cudaEventCreate(&t0)); CK(cudaEventCreate(&t1));
	CK(cudaEventRecord(t0));
	k_occ4<<<blocks,tpb>>>(d_bwt, bwt->primary, d_ks, n_occ, d_out);
	CK(cudaEventRecord(t1)); CK(cudaEventSynchronize(t1));
	CK(cudaGetLastError());
	float ms=0; CK(cudaEventElapsedTime(&ms,t0,t1));
	std::vector<uint64_t> hout(n_occ*4);
	CK(cudaMemcpy(hout.data(), d_out, n_occ*4*sizeof(uint64_t), cudaMemcpyDeviceToHost));
	uint64_t mism=0; for (uint64_t i=0;i<n_occ;i++){ bwtint_t cnt[4]; bwt_occ4(bwt, hk[i], cnt);
		for(int j=0;j<4;j++) if ((uint64_t)cnt[j]!=hout[i*4+j]){ if(mism<10) fprintf(stderr,
			"  MISMATCH k=%llu c=%d cpu=%llu gpu=%llu\n",(unsigned long long)hk[i],j,
			(unsigned long long)cnt[j],(unsigned long long)hout[i*4+j]); mism++; } }
	fprintf(stderr, "[A] %s  mismatches=%llu  time=%.1f ms  throughput=%.1f M-Occ4/s\n",
		mism? "FAIL":"PASS", (unsigned long long)mism, ms, n_occ/1e6/(ms/1e3));

	/* ---------- Test B: exact backward search of real reads ---------- */
	if (reads_fn) {
		FILE *fp=fopen(reads_fn,"r"); if(!fp){ fprintf(stderr,"cannot open %s\n",reads_fn); return 1; }
		std::vector<uint8_t> flat; std::vector<uint64_t> off; std::vector<int> len;
		char *line=NULL; size_t cap=0; ssize_t nl; long ln=0;
		while ((nl=getline(&line,&cap,fp))>=0){
			if ((ln & 3)==1){ /* sequence line */
				int L=0; while (line[L] && line[L]!='\n' && line[L]!='\r') ++L;
				off.push_back(flat.size());
				for (int i=0;i<L;i++) flat.push_back(nt4[(unsigned char)line[i]]);
				len.push_back(L);
			}
			++ln;
		}
		free(line); fclose(fp);
		int nreads=len.size();
		fprintf(stderr, "[B] exact backward-search on %d reads from %s\n", nreads, reads_fn);
		uint8_t *d_seq=NULL; uint64_t *d_off=NULL; int *d_len=NULL;
		uint64_t *d_k=NULL,*d_l=NULL; int *d_c=NULL;
		CK(cudaMalloc(&d_seq, flat.size())); CK(cudaMemcpy(d_seq, flat.data(), flat.size(), cudaMemcpyHostToDevice));
		CK(cudaMalloc(&d_off, off.size()*8)); CK(cudaMemcpy(d_off, off.data(), off.size()*8, cudaMemcpyHostToDevice));
		CK(cudaMalloc(&d_len, len.size()*4)); CK(cudaMemcpy(d_len, len.data(), len.size()*4, cudaMemcpyHostToDevice));
		CK(cudaMalloc(&d_k,nreads*8)); CK(cudaMalloc(&d_l,nreads*8)); CK(cudaMalloc(&d_c,nreads*4));
		int bb=(nreads+tpb-1)/tpb;
		k_match<<<bb,tpb>>>(fm,d_seq,d_off,d_len,nreads,d_k,d_l,d_c);
		CK(cudaDeviceSynchronize()); CK(cudaGetLastError());
		std::vector<uint64_t> gk(nreads),gl(nreads); std::vector<int> gc(nreads);
		CK(cudaMemcpy(gk.data(),d_k,nreads*8,cudaMemcpyDeviceToHost));
		CK(cudaMemcpy(gl.data(),d_l,nreads*8,cudaMemcpyDeviceToHost));
		CK(cudaMemcpy(gc.data(),d_c,nreads*4,cudaMemcpyDeviceToHost));
		uint64_t bmis=0, nhit=0;
		for (int i=0;i<nreads;i++){ bwtint_t k,l; int c=bwt_match_exact(bwt, len[i], flat.data()+off[i], &k,&l);
			if (c>0) nhit++;
			uint64_t ek = c? (uint64_t)k:gk[i], el = c? (uint64_t)l:gl[i]; /* k,l undefined on no-match */
			if (c!=gc[i] || (c>0 && (ek!=gk[i]||el!=gl[i]))){ if(bmis<10) fprintf(stderr,
				"  MISMATCH read %d: cpu(c=%d k=%llu l=%llu) gpu(c=%d k=%llu l=%llu)\n",
				i,c,(unsigned long long)k,(unsigned long long)l,gc[i],
				(unsigned long long)gk[i],(unsigned long long)gl[i]); bmis++; } }
		fprintf(stderr, "[B] %s  reads=%d  full-length-exact-hits=%llu  mismatches=%llu\n",
			bmis?"FAIL":"PASS", nreads, (unsigned long long)nhit, (unsigned long long)bmis);
	}

	cudaFree(d_bwt); cudaFree(d_ks); cudaFree(d_out);
	bwt_destroy(bwt);
	return 0;
}
