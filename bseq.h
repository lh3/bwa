#ifndef BATCHSEQ_H_
#define BATCHSEQ_H_

typedef struct {
	int l_seq;
	char *name, *comment, *seq, *qual;
} bseq1_t;

bseq1_t *bseq_read(int chunk_size, int *n_, void *ks1_, void *ks2_);

#endif
