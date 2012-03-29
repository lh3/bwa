#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#include "bwatpx.h"

extern int num_sampe_threads;
extern int async_read_seq;

#ifdef HAVE_PTHREAD
extern pthread_mutex_t io_lock;
static pthread_t readse_tid = 0;
#endif // HAVE_PTHREAD

static bwa_seqio_t *ks1_copy = NULL;
static int *n_read_addr = 0;
static int mode1_copy = 0;
static int trim1_copy = 0;
static bwa_seq_t **seq1_addr = NULL;
static bwt_aln1_t **aln_addr = NULL;
static int n_occ_copy = 0;
static FILE *fp_sa_addr = NULL;
static int m_aln_copy = 0;

extern void bwa_aln2seq_core(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s, int set_main, int n_multi);

// -------------------

static void thr_bwa_read_seq1_tpx(long n_needed)
{
	int i;
	int n1;

	*seq1_addr = bwa_read_seq(ks1_copy, (int)n_needed, &n1, mode1_copy, trim1_copy);

	if(*seq1_addr == NULL){
		return;
	}

	*n_read_addr = n1;

	// read alignment
	for (i = 0; i < n1; ++i) {
		bwa_seq_t *p = *seq1_addr + i;
		int n_aln;

		fread(&n_aln, 4, 1, fp_sa_addr);

		if (n_aln > m_aln_copy) {
			m_aln_copy = n_aln;
			*aln_addr = (bwt_aln1_t*)realloc(*aln_addr, sizeof(bwt_aln1_t) * m_aln_copy);
		}

		fread(*aln_addr, sizeof(bwt_aln1_t), n_aln, fp_sa_addr);

		bwa_aln2seq_core(n_aln, *aln_addr, p, 1, n_occ_copy);
	}

	return;
}

// -------------------

void bwa_read_seq1_wait_tpx(void)
{
#ifdef HAVE_PTHREAD
        if( (async_read_seq) && (readse_tid != 0) ){
                pthread_join(readse_tid, NULL);
                readse_tid = 0;
        }
#endif // HAVE_PTHREAD

        return;
}

// -------------------

void bwa_read_seq1_tpx(bwa_seqio_t *ks1, int n_needed, int *n,
                       int mode1, int trim_qual1,
                       bwa_seq_t **seq1, bwt_aln1_t **aln, int n_occ, FILE *fp_sa, int m_aln)
{
#ifdef _TIMING
	struct timeval st;
	uint64_t s1, e1;
	double pos1_time = 0.0;
#endif

#ifdef _DEBUG
	pthread_mutex_lock(&io_lock);
	fprintf(stderr,"bwaseio1: n_needed = %d\n",n_needed);
	pthread_mutex_unlock(&io_lock);
#endif

#ifdef _TIMING
	gettimeofday(&st, NULL);
	s1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;
#endif

#ifdef HAVE_PTHREAD
	if( (async_read_seq) && (num_sampe_threads > 1) ){

		if(readse_tid != 0){
			pthread_join(readse_tid, NULL);
			readse_tid = 0;
		}

		ks1_copy = ks1;
		n_read_addr = n;
		mode1_copy = mode1;
		trim1_copy = trim_qual1;
		seq1_addr = seq1;
		aln_addr = aln;
		n_occ_copy = n_occ;
		fp_sa_addr = fp_sa;
		m_aln_copy = m_aln;

		int srtn;
		srtn = pthread_create(&readse_tid,NULL,(void *(*)(void *))thr_bwa_read_seq1_tpx,(void *)(long)n_needed);
		if(srtn != 0){
			fprintf(stderr,"[bwa_read_seq1_tpx] pthread_create thr_bwa_read_seq1_tpx error %d\n",srtn);
			exit(1);
		}

	}else{

		ks1_copy = ks1;
		n_read_addr = n;
		mode1_copy = mode1;
		trim1_copy = trim_qual1;
		seq1_addr = seq1;
		aln_addr = aln;
		n_occ_copy = n_occ;
		fp_sa_addr = fp_sa;
		m_aln_copy = m_aln;

		thr_bwa_read_seq1_tpx(n_needed);

	}
#else // HAVE_PTHREAD
	ks1_copy = ks1;
	n_read_addr = n;
	mode1_copy = mode1;
	trim1_copy = trim_qual1;
	seq1_addr = seq1;
	aln_addr = aln;
	n_occ_copy = n_occ;
	fp_sa_addr = fp_sa;
	m_aln_copy = m_aln;

	thr_bwa_read_seq1_tpx(n_needed);
#endif // HAVE_PTHREAD

#ifdef _TIMING
	gettimeofday(&st, NULL);
	e1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;
	pos1_time = (double)((double)e1 - (double)s1) / 1000000.0;

# ifdef HAVE_PTHREAD
	pthread_mutex_lock(&io_lock);
# endif // HAVE_PTHREAD
	fprintf(stderr,"bwaseio1 time = %lf (sec)\n",pos1_time);
# ifdef HAVE_PTHREAD
	pthread_mutex_unlock(&io_lock);
# endif // HAVE_PTHREAD
#endif

	return;
}
