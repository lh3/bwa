#ifndef BWATPX_H
#define BWATPX_H

#include "bwtaln.h"
#include "bntseq.h"
#include "kvec.h"
#include "kstring.h"

#ifdef HAVE_PTHREAD
# include <pthread.h>
 typedef pthread_t my_pthread_tid_t;
#else // HAVE_PTHREAD
 typedef long      my_pthread_tid_t;
#endif // HAVE_PTHREAD

#define MAX_CPUS 24

#define MIN_HASH_WIDTH 1000

#define SW_MIN_MATCH_LEN 20
#define SW_MIN_MAPQ      17

// mck - Aug 2011 - tpx - threaded sampe

typedef struct {
        int n;
        bwtint_t *a;
} poslist_t;

typedef struct {
        double avg, std, ap_prior;
        bwtint_t low, high, high_bayesian;
} isize_info_t;

typedef struct {
        kvec_t(uint64_t) arr;
        kvec_t(uint64_t) pos[2];
        kvec_t(bwt_aln1_t) aln[2];
} pe_data_t;

typedef struct {
        kvec_t(bwt_aln1_t) aln;
} aln_buf_t;

typedef struct {
        long start;
        long end;
        int cnt_chg;
        bwt_t *bwt[2];
        bwa_seq_t *seqs[2];
        isize_info_t *ii;
        const pe_opt_t *opt;
        const gap_opt_t *gopt;
        pe_data_t *d[MAX_CPUS];
        aln_buf_t *buf[2];
        my_pthread_tid_t tid;
} THR_BWA_PE_TPX;

typedef struct {
        long start;
        long end;
        const bntseq_t *bns;
        ubyte_t *pacseq;
        bwa_seq_t *seqs[2];
        const pe_opt_t *popt;
        const isize_info_t *ii;
        my_pthread_tid_t tid;
} THR_BWA_SW_TPX;

typedef struct {
       long start;
       long end;
       bwt_t *bwt[2];
       bwa_seq_t *seqs[2];
       const gap_opt_t *gopt;
       my_pthread_tid_t tid;
} THR_BWA_SE_TPX;

typedef struct {
       long start;
       long end;
       const bntseq_t *bns;
       bwa_seq_t *seqs;
       ubyte_t *pacseq;
       bntseq_t *ntbns;
       my_pthread_tid_t tid;
} THR_BWA_RG_TPX;

typedef struct {
       long start;
       long end;
       const bwt_t *bwt;
       bwa_seq_t *seqs;
       int max_mm;
       float fnr;
       my_pthread_tid_t tid;
} THR_BWA_SE_PAC_TPX;

#endif

