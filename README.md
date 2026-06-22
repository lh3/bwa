[![Build Status](https://github.com/lh3/bwa/actions/workflows/ci.yaml/badge.svg)](https://github.com/lh3/bwa/actions)
[![SourceForge Downloads](https://img.shields.io/sourceforge/dt/bio-bwa.svg?label=SF%20downloads)](https://sourceforge.net/projects/bio-bwa/files/?source=navbar)
[![GitHub Downloads](https://img.shields.io/github/downloads/lh3/bwa/total.svg?style=flat&label=GitHub%20downloads)](https://github.com/lh3/bwa/releases)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/bwa.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/bwa)

> [!Note]
> BWA-MEM has been replaced by [minimap2][minimap2] for long reads
> and will be replaced by [minibwa][minibwa] for short reads.
> If you still prefer the exact legacy output of BWA-MEM, give [BWA-MEM2][bwa-mem2] a try.
> The original bwa-aln algorithm is still unique due to its sensitivity to very short reads.

[minimap2]: https://github.com/lh3/minimap2
[bwa-mem2]: https://github.com/bwa-mem2/bwa-mem2
[minibwa]: https://github.com/lh3/minibwa

## Getting started

	git clone https://github.com/lh3/bwa.git
	cd bwa; make
	./bwa index ref.fa
	./bwa mem ref.fa read-se.fq.gz | gzip -3 > aln-se.sam.gz
	./bwa mem ref.fa read1.fq read2.fq | gzip -3 > aln-pe.sam.gz

## GPU acceleration: `bwa gpualn` (this fork)

This fork adds **`gpualn`**, a CUDA implementation of BWA-backtrack (`bwa aln`)
for single-end reads, targeting **ancient DNA (aDNA)** where seeding is disabled
(`-l 1024`) — the regime where seed-and-extend aligners (and the older BarraCUDA)
struggle. It produces a **bit-exact** `.sai`/SAM (byte-identical to CPU
`bwa aln`/`samse`) while running **~5x faster than 16 CPU threads** on a single
NVIDIA RTX 3090 (full 3.95M-read aDNA file: 161 s vs 895 s; `.sai` md5 identical).

Build (needs the CUDA toolkit + an NVIDIA GPU; the default `make` is unchanged
and has **no** CUDA dependency):

	make bwa-gpu          # CUDA-enabled `bwa` with the `gpualn` subcommand
	# CUDA_ARCH defaults to sm_86 (RTX 30xx); e.g. `make bwa-gpu CUDA_ARCH=sm_80`

Usage — drop-in for `bwa aln` (writes a `.sai`):

	./bwa-gpu gpualn -l 1024 -n 0.01 -o 2 -t 16 ref.fa reads.fq.gz > out.sai
	./bwa-gpu samse -r '@RG\tID:S\tSM:S' ref.fa out.sai reads.fq.gz | samtools sort -o out.bam -

Or the **fused alnse** mode (`-S`): aln + samse in one command, no intermediate
`.sai`, the FASTQ read once, multithreaded `samse`, SAM to stdout:

	./bwa-gpu gpualn -S -r '@RG\tID:S\tSM:S\tPL:illumina' ref.fa reads.fq.gz \
	    | samtools sort -@16 -O bam -o out.bam -

Options mirror `bwa aln`: `-l` seed length (`-l 1024` disables seeding for aDNA),
`-n` max diff / missing-prob, `-o` max gap opens, `-t` CPU threads (preprocessing
+ the small CPU reconcile), `-f` output file; plus `-S` (SAM/alnse) and `-r`
(read-group line).

**Performance (short-read alignment, `gpualn` vs CPU `bwa aln`+`samse`, single
RTX 3090 vs 16 CPU threads, `-l 1024 -n 0.01 -o 2`, hs37d5):**

| sample (single-end aDNA) | reads | `gpualn` | CPU `bwa aln`+`samse` | speedup |
|--------------------------|-------|----------|-----------------------|---------|
| low-endogenous (~0.5% mapped)    | 3,948,528  | 174 s (≈22.7k reads/s) | ~895 s            | **~5.1×** |
| higher-endogenous (~2% mapped)   | 40,571,012 | 1,947 s (≈20.8k reads/s) | ~9,790 s (wall) | **~5.0×** |

Throughput is ~20,000–28,000 reads/s depending on how much backtracking each
read needs — damaged/non-matching reads (the bulk of aDNA, and seeding is off)
explore larger search trees, so a sample's effective rate tracks its tree sizes,
not just read count. The result is **bit-exact**: on the 40.6 M-read sample all
944,922 mapped records were byte-identical to CPU `bwa aln`. The kernel sits at
~80% of the GPU's random-access FM-index `Occ` bandwidth ceiling, i.e. it is
memory-bound on the index probes rather than compute- or scheduling-bound.

Why it is bit-exact: the GPU runs the BWA-backtrack search as a warp-cooperative
depth-first traversal of the FM-index (one read per warp, a two-level
shared+global stack) and detects which reads have a hit; the ~0.5% that map are
reconciled on the CPU with the **unmodified** `bwt_match_gap`, so output matches
stock `bwa`. The standard `bwa index` is used as-is (the suffix array is only
needed by `samse`, not on the GPU); ~3.1 GB of GPU memory for a human-genome
index (fits an 8 GB card). Design, profiling and validation:
[cuda/PROGRESS.md](cuda/PROGRESS.md). How it was built (LLM-assisted) and the
references it drew on: [cuda/LLM_AND_REFERENCES.md](cuda/LLM_AND_REFERENCES.md).

### Pipeline integration and the short/long split cutoff

`gpualn` is a bit-exact drop-in for the short-read aligner in the aDNA mapping
pipeline ([cuda/adna_aligner.sh](cuda/adna_aligner.sh), `-g`): reads **below** a
length cutoff go to `gpualn` (BWA-backtrack), reads **at/above** it to
`bwa mem`/`mem3`. The cutoff **defaults to a fixed 64 bp**; `-s auto` instead
picks it from the AdapterRemoval3 length report (clamped to a **floor of 60**),
and `-s INT` forces a value.

Why 64 with a floor of 60: BWA-backtrack (`aln`) is most precise on short reads,
and `bwa mem -k19 -r2.5` *loses* precision / inflates reference bias for reads
**< 70 bp**, worst at **30–60 bp** — the bulk of aDNA (Oliva et al. 2021). A
cutoff of 64 keeps that precision-critical band on `aln`/`gpualn` while sitting at
backtrack's cost cliff (`max_diff` reaches 5 at ~64 bp for `-n 0.01`). A real-data
check (VILA01, 11.06 M reads, hs37d5) confirmed that *raising* the cutoff is
counter-productive here: 64 → 70 → 75 monotonically **lost mappings**
(389,722 → 385,453 → 376,414) and **added wall time** (240 → 326 → 470 s) — the
64+ band is slow on `gpualn`'s d=5 search and `mem` recovers more of it (at lower
precision). The floor of 60 stops `auto` from ever exiling the precision-critical
band to `mem`. Use `-s 70` only for a study that is specifically reference-bias-
sensitive and will pay the ~+36% runtime.

### GPU alignment *generation* (investigation, not shipped)

Production `gpualn` is GPU hit-*detection* + a small CPU reconcile. For the
expensive d=5 band (64–70 bp) the CPU reconcile of every mapped read is a floor
that pins `gpualn` at *parity* with CPU `bwa aln` on that band. We investigated
generating the alignments on the GPU to remove that floor:

- `bwt_match_gap`'s output is **order-sensitive** (best-first + the `gap_shadow`
  bound mutation), and order-sensitivity turns out to be **equivalent to
  multi-mapping** (`c1>1`, the `XT:A:R` reads) — validated *with gaps* against a
  brute-force edit-distance oracle at ~99.97–100 % recall. So unique-best reads
  (86 % overall, 90 % of the 64–69 band on real data) are order-free and *can* be
  generated bit-exactly on the GPU.
- A **warp-cooperative collect-and-sort prototype is bit-exact but
  throughput-counter-productive**: unique reads have small search trees that waste
  31/32 lanes, while the big-tree reads warp-cooperation would help are exactly
  the multi-mappers routed to the CPU anyway. Per-thread best-first is the right
  mapping and is only a modest ~1.2–1.4× over CPU on a cache-resident benchmark.

Conclusion: bit-exact GPU *generation* of the 64–70 band is achievable but the
realistic speedup is modest — the **cutoff above is the practical lever**, not a
new engine. Full scoping:
[GPU_ALIGN_GENERATION_SCOPING.md](GPU_ALIGN_GENERATION_SCOPING.md); domain-neutral
engine writeup: [GPU_BOUNDED_DFS_REPORT.md](GPU_BOUNDED_DFS_REPORT.md). Opt-in
per-length-band instrumentation: run `gpualn` with `GPUALN_HISTO=1` (prints a
node-pop/flag histogram by read length + GPU-kernel-vs-CPU-reconcile timing).
Self-contained CUDA self-tests of the engines (build with `nvcc`, no other deps):
[cuda/fmdfs_selftest.cu](cuda/fmdfs_selftest.cu) (existence search, bit-exact vs a
CPU oracle + brute force) and
[cuda/bestfirst_selftest.cu](cuda/bestfirst_selftest.cu) (order-sensitive
best-first + the warp-cooperative prototype + the `c1>1` validation).

## Introduction

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

## Availability

BWA is released under [GPLv3][1]. The latest source code is [freely
available at github][2]. Released packages can [be downloaded][3] at
SourceForge. After you acquire the source code, simply use `make` to compile
and copy the single executable `bwa` to the destination you want. The only
dependency required to build BWA is [zlib][14].

Since 0.7.11, precompiled binary for x86\_64-linux is available in [bwakit][17].
In addition to BWA, this self-consistent package also comes with bwa-associated
and 3rd-party tools for proper BAM-to-FASTQ conversion, mapping to ALT contigs,
adapter triming, duplicate marking, HLA typing and associated data files.

## Seeking help

The detailed usage is described in the man page available together with the
source code. You can use `man ./bwa.1` to view the man page in a terminal. The
[HTML version][4] of the man page can be found at the [BWA website][5]. If you
have questions about BWA, you may [sign up the mailing list][6] and then send
the questions to [bio-bwa-help@sourceforge.net][7]. You may also ask questions
in forums such as [BioStar][8] and [SEQanswers][9].

## Citing BWA

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

## Frequently asked questions (FAQs)

1. [What types of data does BWA work with?](#type)
2. [Why does a read appear multiple times in the output SAM?](#multihit)
3. [Does BWA work on reference sequences longer than 4GB in total?](#4gb)
4. [Why can one read in a pair has high mapping quality but the other has zero?](#pe0)
5. [How can a BWA-backtrack alignment stands out of the end of a chromosome?](#endref)
6. [Does BWA work with ALT contigs in the GRCh38 release?](#altctg)
7. [Can I just run BWA-MEM against GRCh38+ALT without post-processing?](#postalt)
8. [Why does BWA use a lot of memory?](#largemem)

#### <a name="type"></a>1. What types of data does BWA work with?

BWA works with a variety types of DNA sequence data, though the optimal
algorithm and setting may vary. The following list gives the recommended
settings:

* Illumina/454/IonTorrent single-end reads longer than ~70bp or assembly
  contigs up to a few megabases mapped to a closely related reference genome:

		bwa mem ref.fa reads.fq > aln.sam

* Illumina single-end reads shorter than ~70bp:

		bwa aln ref.fa reads.fq > reads.sai; bwa samse ref.fa reads.sai reads.fq > aln-se.sam

* Illumina/454/IonTorrent paired-end reads longer than ~70bp:

		bwa mem ref.fa read1.fq read2.fq > aln-pe.sam

* Illumina paired-end reads shorter than ~70bp:

		bwa aln ref.fa read1.fq > read1.sai; bwa aln ref.fa read2.fq > read2.sai
		bwa sampe ref.fa read1.sai read2.sai read1.fq read2.fq > aln-pe.sam

* PacBio subreads or Oxford Nanopore reads to a reference genome:

		bwa mem -x pacbio ref.fa reads.fq > aln.sam
		bwa mem -x ont2d ref.fa reads.fq > aln.sam

BWA-MEM is recommended for query sequences longer than ~70bp for a variety of
error rates (or sequence divergence). Generally, BWA-MEM is more tolerant with
errors given longer query sequences as the chance of missing all seeds is small.
As is shown above, with non-default settings, BWA-MEM works with Oxford Nanopore
reads with a sequencing error rate over 20%.

#### <a name="multihit"></a>2. Why does a read appear multiple times in the output SAM?

BWA-SW and BWA-MEM perform local alignments. If there is a translocation, a gene
fusion or a long deletion, a read bridging the break point may have two hits,
occupying two lines in the SAM output. With the default setting of BWA-MEM, one
and only one line is primary and is soft clipped; other lines are tagged with
0x800 SAM flag (supplementary alignment) and are hard clipped.

#### <a name="4gb"></a>3. Does BWA work on reference sequences longer than 4GB in total?

Yes. Since 0.6.x, all BWA algorithms work with a genome with total length over
4GB. However, individual chromosome should not be longer than 2GB.

#### <a name="pe0"></a>4. Why can one read in a pair have a high mapping quality but the other has zero?

This is correct. Mapping quality is assigned for individual read, not for a read
pair. It is possible that one read can be mapped unambiguously, but its mate
falls in a tandem repeat and thus its accurate position cannot be determined.

#### <a name="endref"></a>5. How can a BWA-backtrack alignment stand out of the end of a chromosome?

Internally BWA concatenates all reference sequences into one long sequence. A
read may be mapped to the junction of two adjacent reference sequences. In this
case, BWA-backtrack will flag the read as unmapped (0x4), but you will see
position, CIGAR and all the tags. A similar issue may occur to BWA-SW alignment
as well. BWA-MEM does not have this problem.

#### <a name="altctg"></a>6. Does BWA work with ALT contigs in the GRCh38 release?

Yes, since 0.7.11, BWA-MEM officially supports mapping to GRCh38+ALT.
BWA-backtrack and BWA-SW don't properly support ALT mapping as of now. Please
see [README-alt.md][18] for details. Briefly, it is recommended to use
[bwakit][17], the binary release of BWA, for generating the reference genome
and for mapping.

#### <a name="postalt"></a>7. Can I just run BWA-MEM against GRCh38+ALT without post-processing?

If you are not interested in hits to ALT contigs, it is okay to run BWA-MEM
without post-processing. The alignments produced this way are very close to
alignments against GRCh38 without ALT contigs. Nonetheless, applying
post-processing helps to reduce false mappings caused by reads from the
diverged part of ALT contigs and also enables HLA typing. It is recommended to
run the post-processing script.

### <a name="largemem"></a>8. Why does BWA use a lot of memory?

This is typically caused by FASTQ generated from a coordinate-sorted BAM.
BWA uses a lot more memory for centromeric reads than for unique reads.
In a FASTQ file generated from a sequencing run, centromeric reads are rare in each batch and rarely cause troubles.
However, in a coordinate-sorted FASTQ file, a whole batch could consist of centromeric reads.
Such a batch will take a lot more memory and time to map; the insert size estimate will be distorted as well.
General rule: ***NEVER*** use Picard SamToFastq on coordiate-sorted BAM;
use samtools [collate+fastq][remap] instead.

[remap]: https://lh3.github.io/2021/07/06/remapping-an-aligned-bam
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
[17]: http://sourceforge.net/projects/bio-bwa/files/bwakit/
[18]: https://github.com/lh3/bwa/blob/master/README-alt.md
