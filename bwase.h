#ifndef BWASE_H
#define BWASE_H

#include "bntseq.h"
#include "bwt.h"
#include "bwtaln.h"

#ifdef __cplusplus
extern "C" {
#endif

	// Initialize mapping tables in the bwa single-end mapper.
	void bwase_initialize();
	// Calculate the approximate position of the sequence from the specified bwt with loaded suffix array.
	void bwa_cal_pac_pos_core(const bwt_t* forward_bwt, const bwt_t* reverse_bwt, bwa_seq_t* seq, const int max_mm, const float fnr);
	// Refine the approximate position of the sequence to an actual placement for the sequence.
	void bwa_refine_gapped(const bntseq_t *bns, int n_seqs, bwa_seq_t *seqs, ubyte_t *_pacseq, bntseq_t *ntbns);
	// Backfill certain alignment properties mainly centering around number of matches.
	void bwa_aln2seq(int n_aln, const bwt_aln1_t *aln, bwa_seq_t *s);
	// Calculate the end position of a read given a certain sequence.
	int64_t pos_end(const bwa_seq_t *p);

#ifdef __cplusplus
}
#endif

#endif // BWASE_H
