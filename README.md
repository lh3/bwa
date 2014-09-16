##Getting started

	git clone https://github.com/lh3/bwa.git
	cd bwa; make
	./bwa index ref.fa
	./bwa mem ref.fa read-se.fq.gz | gzip -3 > aln-se.sam.gz
	./bwa mem ref.fa read1.fq read2.fq | gzip -3 > aln-pe.sam.gz

##Introduction

BWA is a software package for mapping DNA sequences against a large reference
genome, such as the human genome. It consists of three algorithms:
BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina
sequence reads up to 100bp, while the rest two for longer sequences ranged from
70bp to a few megabases. BWA-MEM and BWA-SW share similar features such as the
support of long reads and chimeric alignment, but BWA-MEM, which is the latest,
is generally recommended as it is faster and more accurate. BWA-MEM also has
better performance than BWA-backtrack for 70-100bp Illumina reads.

For all the algorithms, BWA first needs to construct the FM-index for the
reference genome (the **index** command). Alignment algorithms are invoked with
different sub-commands: **aln/samse/sampe** for BWA-backtrack,
**bwasw** for BWA-SW and **mem** for the BWA-MEM algorithm.

##Availability

BWA is released under [GPLv3][1]. The latest source code is [freely
available at github][2]. Released packages can [be downloaded][3] at
SourceForge. After you acquire the source code, simply use `make` to compile
and copy the single executable `bwa` to the destination you want. The only
dependency required to build BWA is [zlib][14].

##Seeking helps

The detailed usage is described in the man page available together with the
source code. You can use `man ./bwa.1` to view the man page in a terminal. The
[HTML version][4] of the man page can be found at the [BWA website][5]. If you
have questions about BWA, you may [sign up the mailing list][6] and then send
the questions to [bio-bwa-help@sourceforge.net][7]. You may also ask questions
in forums such as [BioStar][8] and [SEQanswers][9].

##Citing BWA

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

##Frequently asked questions (FAQs)

1. [What types of data does BWA work with?](#type)
2. [Why does a read appear multiple times in the output SAM?](#multihit)
3. [Does BWA work on reference sequences longer than 4GB in total?](#4gb)
4. [Why can one read in a pair has high mapping quality but the other has zero?](#pe0)
5. [How can a BWA-backtrack alignment stands out of the end of a chromosome?](#endref)
6. [How to map sequences to GRCh38 with ALT contigs?](#h38)

####<a name="type"></a>1. What types of data does BWA work with?

BWA works with a variety types of DNA sequence data, though the optimal
algorithm and setting may vary. The following list gives the recommended
settings:

* Illumina/454/IonTorrent single-end reads longer than ~70bp or assembly
  contigs up to a few megabases mapped to a close related reference genome:

		bwa mem ref.fa reads.fq > aln.sam

* Illumina single-end reads no longer than ~70bp:

		bwa aln ref.fa reads.fq > reads.sai; bwa samse ref.fa reads.sai reads.fq > aln-se.sam

* Illumina/454/IonTorrent paired-end reads longer than ~70bp:

		bwa mem ref.fa read1.fq read2.fq > aln-pe.sam

* Illumina paired-end reads no longer than ~70bp:

		bwa aln ref.fa read1.fq > read1.sai; bwa aln ref.fa read2.fq > read2.sai
		bwa sampe ref.fa read1.sai read2.sai read1.fq read2.fq > aln-pe.sam

* PacBio subreads to a reference genome:

		bwa mem -x pacbio ref.fa reads.fq > aln.sam

BWA-MEM is recommended for query sequences longer than ~70bp for a variety of
error rates (or sequence divergence). Generally, BWA-MEM is more tolerant with
errors given longer query sequences as the chance of missing all seeds is small.
As is shown above, with non-default settings, BWA-MEM works with PacBio subreads
with a sequencing error rate as high as ~15%.

####<a name="multihit"></a>2. Why does a read appear multiple times in the output SAM?

BWA-SW and BWA-MEM perform local alignments. If there is a translocation, a gene
fusion or a long deletion, a read bridging the break point may have two hits,
occupying two lines in the SAM output. With the default setting of BWA-MEM, one
and only one line is primary and is soft clipped; other lines are tagged with
0x800 SAM flag (supplementary alignment) and are hard clipped.

####<a name="4gb"></a>3. Does BWA work on reference sequences longer than 4GB in total?

Yes. Since 0.6.x, all BWA algorithms work with a genome with total length over
4GB. However, individual chromosome should not be longer than 2GB.

####<a name="pe0"></a>4. Why can one read in a pair has high mapping quality but the other has zero?

This is correct. Mapping quality is assigned for individual read, not for a read
pair. It is possible that one read can be mapped unambiguously, but its mate
falls in a tandem repeat and thus its accurate position cannot be determined.

####<a name="endref"></a>5. How can a BWA-backtrack alignment stands out of the end of a chromosome?

Internally BWA concatenates all reference sequences into one long sequence. A
read may be mapped to the junction of two adjacent reference sequences. In this
case, BWA-backtrack will flag the read as unmapped (0x4), but you will see
position, CIGAR and all the tags. A similar issue may occur to BWA-SW alignment
as well. BWA-MEM does not have this problem.



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
[15]: https://github.com/lh3/bwa/tree/mem
[16]: ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/
