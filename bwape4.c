#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#include "bwatpx.h"

extern void bwa_print_sam1(const bntseq_t *bns, bwa_seq_t *p, const bwa_seq_t *mate, int mode, int max_top2);

extern int num_sampe_threads;
extern int async_print_res;

#ifdef HAVE_PTHREAD
extern pthread_mutex_t pe_lock;
static pthread_t printpe_tid = 0;
#endif // HAVE_PTHREAD

static const bntseq_t *bns_copy = NULL;
static bwa_seq_t *seqs_copy[2] = { NULL, NULL };
static int opt_mode = 0;
static int opt_max_top2 = 0;

// -------------------

void thr_bwa_print2_tpx(long n_seqs)
{
	int i = 0;

	for (i = 0; i < n_seqs; ++i) {

		bwa_seq_t *p[2];

		p[0] = seqs_copy[0] + i; 
		p[1] = seqs_copy[1] + i;

		if (p[0]->bc[0] || p[1]->bc[0]) {
			strcat(p[0]->bc, p[1]->bc);
			strcpy(p[1]->bc, p[0]->bc);
		}

		bwa_print_sam1(bns_copy, p[0], p[1], opt_mode, opt_max_top2);
		bwa_print_sam1(bns_copy, p[1], p[0], opt_mode, opt_max_top2);

	}

	return;
}

// -------------------

void bwa_print2_wait_tpx(void)
{
#ifdef HAVE_PTHREAD
	if( (async_print_res) && (printpe_tid != 0) ){
		pthread_join(printpe_tid, NULL);
		printpe_tid = 0;
	}
#endif // HAVE_PTHREAD

	return;
}

// -------------------

void bwa_print2_tpx(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs[2], const gap_opt_t opt)
{
#ifdef _TIMING
	struct timeval st;
	uint64_t s1, e1;
	double pos1_time = 0.0;
#endif

#ifdef _DEBUG
# ifdef HAVE_PTHREAD
        pthread_mutex_lock(&pe_lock);
# endif // HAVE_PTHREAD
        fprintf(stderr,"bwape4:\n");
# ifdef HAVE_PTHREAD
        pthread_mutex_unlock(&pe_lock);
# endif // HAVE_PTHREAD
#endif

#ifdef _TIMING
	gettimeofday(&st, NULL);
	s1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;
#endif

# ifdef HAVE_PTHREAD
	if( (async_print_res) && (num_sampe_threads > 1) ){

		if(printpe_tid != 0){	
			pthread_join(printpe_tid, NULL);
			printpe_tid = 0;
		}

		seqs_copy[0] = seqs[0];	
		seqs_copy[1] = seqs[1];	
		bns_copy = bns;
		opt_mode = opt.mode;
		opt_max_top2 = opt.max_top2;

		int srtn;	
		srtn = pthread_create(&printpe_tid,NULL,(void *(*)(void *))thr_bwa_print2_tpx,(void *)(long)n_seqs);
		if(srtn != 0){
			fprintf(stderr,"[bwa_print2_tpx] pthread_create thr_bwa_print2 error %d\n",srtn);
			exit(1);
		}

	}else{

		seqs_copy[0] = seqs[0];	
		seqs_copy[1] = seqs[1];	
		bns_copy = bns;
		opt_mode = opt.mode;
		opt_max_top2 = opt.max_top2;

		thr_bwa_print2_tpx(n_seqs);

	}
# else // HAVE_PTHREAD
	seqs_copy[0] = seqs[0];	
	seqs_copy[1] = seqs[1];	
	bns_copy = bns;
	opt_mode = opt.mode;
	opt_max_top2 = opt.max_top2;

	thr_bwa_print2_tpx(n_seqs);
# endif // HAVE_PTHREAD

#ifdef _TIMING
	gettimeofday(&st, NULL);
	e1 = st.tv_sec * 1000000L + (time_t)st.tv_usec;
        pos1_time = (double)((double)e1 - (double)s1) / 1000000.0;

# ifdef HAVE_PTHREAD
        pthread_mutex_lock(&pe_lock);
# endif // HAVE_PTHREAD
        fprintf(stderr,"bwape4 time = %lf (sec)\n",pos1_time);
# ifdef HAVE_PTHREAD
        pthread_mutex_unlock(&pe_lock);
# endif // HAVE_PTHREAD
#endif

	return;
}

