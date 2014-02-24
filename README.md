###Information about the BWA Cilk version

This repository contains an adapted version of BWA where the main working loop in BWA aln is parallelised using the Cilk Plus language extension for C, instead of the original pthread version. Depending on work loads and machine configuration, this can yield a performance improvement of up to a factor two. See http://www.exascience.com/wp-content/uploads/2013/12/Herzeel-BWAReport.pdf for a technical report that summarises our findings. Note: The improvements only apply to BWA aln, not BWA mem!

You need a recent Intel C compiler with built-in support for Cilk Plus to compile this code. See http://software.intel.com/en-us/intel-sdp-home for information about Intel compilers. The code is also likely to work with the experimental open source implementation of Cilk Plus for gcc available at https://www.cilkplus.org/download#gcc-development-branch - but we haven’t tested this ourselves.

After cloning this repository, you have to check out the cilk branch (git checkout cilk). You can edit the Makefile to switch between Cilk and pthread execution. For Cilk Plus, the number of worker threads does not need to be specified at the command line, because the Cilk runtime will figure out the number of available cores itself. For scaling experiments, you can use the CILK_NWORKERS environment variable to specify the number of Cilk worker threads. Please check the Cilk Plus documentation at https://www.cilkplus.org/cilk-documentation-full for more details.

###Getting started

	git clone https://github.com/exascience/bwa.git
	cd bwa; git checkout cilk; make
	cd bwa; make
	./bwa aln ...

###Introduction

BWA is a software package for mapping low-divergent sequences against a large
reference genome, such as the human genome. It consists of three algorithms:
BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina
sequence reads up to 100bp, while the rest two for longer sequences ranged from
70bp to 1Mbp. BWA-MEM and BWA-SW share similar features such as the support of
long reads and chimeric alignment, but BWA-MEM, which is the latest, is
generally recommended for high-quality queries as it is faster and more
accurate. BWA-MEM also has better performance than BWA-backtrack for 70-100bp
Illumina reads.

For all the algorithms, BWA first needs to construct the FM-index for the
reference genome (the **index** command). Alignment algorithms are invoked with
different sub-commands: **aln/samse/sampe** for BWA-backtrack,
**bwasw** for BWA-SW and **mem** for the BWA-MEM algorithm.

###Availability

BWA is released under [GPLv3][1]. The latest souce code is [freely
available][2] at github. Released packages can [be downloaded ][3] at
SourceForge. After you acquire the source code, simply use `make` to compile
and copy the single executable `bwa` to the destination you want. The only
dependency of BWA is [zlib][14].

###Seeking helps

The detailed usage is described in the man page available together with the
source code. You can use `man ./bwa.1` to view the man page in a terminal. The
[HTML version][4] of the man page can be found at the [BWA website][5]. If you
have questions about BWA, you may [sign up the mailing list][6] and then send
the questions to [bio-bwa-help@sourceforge.net][7]. You may also ask questions
in forums such as [BioStar][8] and [SEQanswers][9].

###Citing BWA

* Li H. and Durbin R. (2009) Fast and accurate short read alignment with
 Burrows-Wheeler transform. *Bioinformatics*, **25**, 1754-1760. [PMID:
 [19451168][10]]. (if you use the BWA-backtrack algorithm)

* Li H. and Durbin R. (2010) Fast and accurate long-read alignment with
 Burrows-Wheeler transform. *Bioinformatics*, **26**, 589-595. [PMID:
 [20080505][11]]. (if you use the BWA-SW algorithm)

* Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs
 with BWA-MEM. [arXiv:1303.3997v2][12] [q-bio.GN]. (if you use the BWA-MEM
 algorithm or the **fastmap** command, or want to cite the whole BWA package)

Please note that the last reference is a preprint hosted at [arXiv.org][13]. I
do not have plan to submit it to a peer-reviewed journal in the near future.



[1]: http://en.wikipedia.org/wiki/GNU_General_Public_License
[2]: https://github.com/lh3/bwa
[3]: http://sourceforge.net/projects/bio-bwa/files/
[4]: http://bio-bwa.sourceforge.net/bwa.shtml
[5]: http://bio-bwa.sourceforge.net/
[6]: https://lists.sourceforge.net/lists/listinfo/bio-bwa-help
[7]: mailto:bio-bwa-help@sourceforge.net
[8]: http://biostars.org
[9]: http://seqanswers.com/
[10]: http://www.ncbi.nlm.nih.gov/pubmed/19451168
[11]: http://www.ncbi.nlm.nih.gov/pubmed/20080505
[12]: http://arxiv.org/abs/1303.3997
[13]: http://arxiv.org/
[14]: http://zlib.net/
