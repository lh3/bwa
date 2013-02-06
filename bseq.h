#ifndef BATCHSEQ_H_
#define BATCHSEQ_H_

typedef struct {
	char *name, *comment, *seq, *qual;
} bseq1_t;

typedef struct {
	int n, m;
	bseq1_t *seqs;
} bseq_t;

int bseq_read(int chunk_size, bseq_t *bs);

#endif
