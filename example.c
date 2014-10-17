#include <stdio.h>
#include <zlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include "bwamem.h"
#include "kseq.h" // for the FASTA/Q parser
KSEQ_DECLARE(gzFile)

int main(int argc, char *argv[])
{
	bwaidx_t *idx;
	gzFile fp;
	kseq_t *ks;
	mem_opt_t *opt;

	if (argc < 3) {
		fprintf(stderr, "Usage: bwamem-lite <idx.base> <reads.fq>\n");
		return 1;
	}

	idx = bwa_idx_load(argv[1], BWA_IDX_ALL); // load the BWA index
	if (NULL == idx) {
		fprintf(stderr, "Index load failed.\n");
		exit(EXIT_FAILURE);
	}
	fp = strcmp(argv[2], "-")? gzopen(argv[2], "r") : gzdopen(fileno(stdin), "r");
	if (NULL == fp) {
		fprintf(stderr, "Couldn't open %s : %s\n",
				strcmp(argv[2], "-") ? argv[2] : "stdin",
				errno ? strerror(errno) : "Out of memory");
		exit(EXIT_FAILURE);
	}
	ks = kseq_init(fp); // initialize the FASTA/Q parser
	opt = mem_opt_init(); // initialize the BWA-MEM parameters to the default values

	while (kseq_read(ks) >= 0) { // read one sequence
		mem_alnreg_v ar;
		int i, k;
		ar = mem_align1(opt, idx->bwt, idx->bns, idx->pac, ks->seq.l, ks->seq.s); // get all the hits
		for (i = 0; i < ar.n; ++i) { // traverse each hit
			mem_aln_t a;
			if (ar.a[i].secondary >= 0) continue; // skip secondary alignments
			a = mem_reg2aln(opt, idx->bns, idx->pac, ks->seq.l, ks->seq.s, &ar.a[i]); // get forward-strand position and CIGAR
			// print alignment
			printf("%s\t%c\t%s\t%ld\t%d\t", ks->name.s, "+-"[a.is_rev], idx->bns->anns[a.rid].name, (long)a.pos, a.mapq);
			for (k = 0; k < a.n_cigar; ++k) // print CIGAR
				printf("%d%c", a.cigar[k]>>4, "MIDSH"[a.cigar[k]&0xf]);
			printf("\t%d\n", a.NM); // print edit distance
			free(a.cigar); // don't forget to deallocate CIGAR
		}
		free(ar.a); // and deallocate the hit list
	}

	free(opt);
	kseq_destroy(ks);
	gzclose(fp);
	bwa_idx_destroy(idx);
	return 0;
}
