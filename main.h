#ifndef BWA_MAIN_H
#define BWA_MAIN_H

#ifdef __cplusplus
extern "C" {
#endif

	int bwa_fa2pac(int argc, char *argv[]);
	int bwa_pac2cspac(int argc, char *argv[]);
	int bwa_pac2bwt(int argc, char *argv[]);
	int bwa_bwtupdate(int argc, char *argv[]);
	int bwa_bwt2sa(int argc, char *argv[]);
	int bwa_index(int argc, char *argv[]);
	int bwa_aln(int argc, char *argv[]);
	int bwt_bwtgen_main(int argc, char *argv[]);

	int bwa_sai2sam_se(int argc, char *argv[]);
	int bwa_sai2sam_pe(int argc, char *argv[]);

	int bwa_stdsw(int argc, char *argv[]);

	int bwa_bwtsw2(int argc, char *argv[]);

	int main_fastmap(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif
