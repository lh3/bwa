## For the Impatient

```sh
# Download bwakit (or from <http://sourceforge.net/projects/bio-bwa/files/bwakit/> manually)
wget -O- http://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.12_x64-linux.tar.bz2/download \
  | gzip -dc | tar xf -
# Generate the GRCh38+ALT+decoy+HLA and create the BWA index
bwa.kit/run-gen-ref hs38DH   # download GRCh38 and write hs38DH.fa
bwa.kit/bwa index hs38DH.fa  # create BWA index
# mapping
bwa.kit/run-bwamem -o out -H hs38DH.fa read1.fq read2.fq | sh  # skip "|sh" to show command lines
```

This generates `out.aln.bam` as the final alignment, `out.hla.top` for best HLA
genotypes on each gene and `out.hla.all` for other possible HLA genotypes.
Please check out [bwa/bwakit/README.md][kithelp] for details.

## Background

GRCh38 consists of several components: chromosomal assembly, unlocalized contigs
(chromosome known but location unknown), unplaced contigs (chromosome unknown)
and ALT contigs (long clustered variations). The combination of the first three
components is called the *primary assembly*. It is recommended to use the
complete primary assembly for all analyses. Using ALT contigs in read mapping is
tricky.

GRCh38 ALT contigs are totaled 109Mb in length, spanning 60Mbp of the primary
assembly. However, sequences that are highly diverged from the primary assembly
only contribute a few million bp. Most subsequences of ALT contigs are nearly
identical to the primary assembly. If we align sequence reads to GRCh38+ALT
blindly, we will get many additional reads with zero mapping quality and miss
variants on them. It is crucial to make mappers aware of ALTs.

BWA-MEM is ALT-aware. It essentially computes mapping quality across the
non-redundant content of the primary assembly plus the ALT contigs and is free
of the problem above.

## Methods

### Sequence alignment

As of now, ALT mapping is done in two separate steps: BWA-MEM mapping and
postprocessing. The `bwa.kit/run-bwamem` script performs the two steps when ALT
contigs are present. The following picture shows an example about how BWA-MEM
infers mapping quality and reports alignment after step 2:

![](http://lh3lh3.users.sourceforge.net/images/alt-demo.png)

#### Step 1: BWA-MEM mapping

At this step, BWA-MEM reads the ALT contig names from "*idxbase*.alt", ignoring
the ALT-to-ref alignment, and labels a potential hit as *ALT* or *non-ALT*,
depending on whether the hit lands on an ALT contig or not. BWA-MEM then reports
alignments and assigns mapQ following these two rules:

1. The mapQ of a non-ALT hit is computed across non-ALT hits only. The mapQ of
   an ALT hit is computed across all hits.

2. If there are no non-ALT hits, the best ALT hit is outputted as the primary
   alignment. If there are both ALT and non-ALT hits, non-ALT hits will be
   primary and ALT hits be supplementary (SAM flag 0x800).

In theory, non-ALT alignments from step 1 should be identical to alignments
against the reference genome with ALT contigs. In practice, the two types of
alignments may differ in rare cases due to seeding heuristics. When an ALT hit
is significantly better than non-ALT hits, BWA-MEM may miss seeds on the
non-ALT hits.

If we don't care about ALT hits, we may skip postprocessing (step 2).
Nonetheless, postprocessing is recommended as it improves mapQ and gives more
information about ALT hits.

#### Step 2: Postprocessing

Postprocessing is done with a separate script `bwa-postalt.js`. It reads all
potential hits reported in the XA tag, lifts ALT hits to the chromosomal
positions using the ALT-to-ref alignment, groups them based on overlaps between
their lifted positions, and then re-estimates mapQ across the best scoring hit
in each group. Being aware of the ALT-to-ref alignment, this script can greatly
improve mapQ of ALT hits and occasionally improve mapQ of non-ALT hits. It also
writes each hit overlapping the reported hit into a separate SAM line. This
enables variant calling on each ALT contig independent of others.

### On the completeness of GRCh38+ALT

While GRCh38 is much more complete than GRCh37, it is still missing some true
human sequences. To make sure every piece of sequence in the reference assembly
is correct, the [Genome Reference Consortium][grc] (GRC) require each ALT contig
to have enough support from multiple sources before considering to add it to the
reference assembly. This careful and sophisticated procedure has left out some
sequences, one of which is [this example][novel], a 10kb contig assembled from
CHM1 short reads and present also in NA12878. You can try [BLAT][blat] or
[BLAST][blast] to see where it maps.

For a more complete reference genome, we compiled a new set of decoy sequences
from GenBank clones and the de novo assembly of 254 public [SGDP][sgdp] samples.
The sequences are included in `hs38DH-extra.fa` from the [BWA binary
package][res].

In addition to decoy, we also put multiple alleles of HLA genes in
`hs38DH-extra.fa`. These genomic sequences were acquired from [IMGT/HLA][hladb],
version 3.18.0 and are used to collect reads sequenced from these genes.

### HLA typing

HLA genes are known to be associated with many autoimmune diseases, infectious
diseases and drug responses. They are among the most important genes but are
rarely studied by WGS projects due to the high sequence divergence between
HLA genes and the reference genome in these regions.

By including the HLA gene regions in the reference assembly as ALT contigs, we
are able to effectively identify reads coming from these genes. We also provide
a pipeline, which is included in the [BWA binary package][res], to type the
several classic HLA genes. The pipeline is conceptually simple. It de novo
assembles sequence reads mapped to each gene, aligns exon sequences of each
allele to the assembled contigs and then finds the pairs of alleles that best
explain the contigs. In practice, however, the completeness of IMGT/HLA and
copy-number changes related to these genes are not so straightforward to
resolve. HLA typing may not always be successful. Users may also consider to use
other programs for typing such as [Warren et al (2012)][hla4], [Liu et al
(2013)][hla2], [Bai et al (2014)][hla3] and [Dilthey et al (2014)][hla1], though
most of them are distributed under restrictive licenses.

## Preliminary Evaluation

To check whether GRCh38 is better than GRCh37, we mapped the CHM1 and NA12878
unitigs to GRCh37 primary (hs37), GRCh38 primary (hs38) and GRCh38+ALT+decoy
(hs38DH), and called small variants from the alignment. CHM1 is haploid.
Ideally, heterozygous calls are false positives (FP). NA12878 is diploid. The
true positive (TP) heterozygous calls from NA12878 are approximately equal
to the difference between NA12878 and CHM1 heterozygous calls. A better assembly
should yield higher TP and lower FP. The following table shows the numbers for
these assemblies:

|Assembly|hs37   |hs38   |hs38DH|CHM1_1.1|  huref|
|:------:|------:|------:|------:|------:|------:|
|FP      | 255706| 168068| 142516|307172 | 575634|
|TP      |2142260|2163113|2150844|2167235|2137053|

With this measurement, hs38 is clearly better than hs37. Genome hs38DH reduces
FP by ~25k but also reduces TP by ~12k. We manually inspected variants called
from hs38 only and found the majority of them are associated with excessive read
depth, clustered variants or weak alignment. We believe most hs38-only calls are
problematic. In addition, if we compare two NA12878 replicates from HiSeq X10
with nearly identical library construction, the difference is ~140k, an order
of magnitude higher than the difference between hs38 and hs38DH. ALT contigs,
decoy and HLA genes in hs38DH improve variant calling and enable the analyses of
ALT contigs and HLA typing at little cost.

## Problems and Future Development

There are some uncertainties about ALT mappings - we are not sure whether they
help biological discovery and don't know the best way to analyze them. Without
clear demand from downstream analyses, it is very difficult to design the
optimal mapping strategy. The current BWA-MEM method is just a start. If it
turns out to be useful in research, we will probably rewrite bwa-postalt.js in C
for performance; if not, we may make changes. It is also possible that we might
make breakthrough on the representation of multiple genomes, in which case, we
can even get rid of ALT contigs for good.



[res]: https://sourceforge.net/projects/bio-bwa/files/bwakit
[sb]: https://github.com/GregoryFaust/samblaster
[grc]: http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/
[novel]: https://gist.github.com/lh3/9935148b71f04ba1a8cc
[blat]: https://genome.ucsc.edu/cgi-bin/hgBlat
[blast]: http://blast.st-va.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
[sgdp]: http://www.simonsfoundation.org/life-sciences/simons-genome-diversity-project/
[hladb]: http://www.ebi.ac.uk/ipd/imgt/hla/
[grcdef]: http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/info/definitions.shtml
[hla1]: http://biorxiv.org/content/early/2014/07/08/006973
[hlalink]: http://www.hladiseaseassociations.com
[hlatools]: https://www.biostars.org/p/93245/
[hla2]: http://nar.oxfordjournals.org/content/41/14/e142.full.pdf+html
[hla3]: http://www.biomedcentral.com/1471-2164/15/325
[hla4]: http://genomemedicine.com/content/4/12/95
[kithelp]: https://github.com/lh3/bwa/tree/master/bwakit
