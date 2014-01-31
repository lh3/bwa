###Getting started

	git clone https://github.com/lh3/bwa.git
	cd bwa; make
	./bwa index ref.fa
	./bwa mem ref.fa read-se.fq.gz | gzip -3 > aln-se.sam.gz
	./bwa mem ref.fa read1.fq read2.fq | gzip -3 > aln-pe.sam.gz

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
