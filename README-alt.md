## Getting Started

Since version 0.7.11, BWA-MEM supports read mapping against a reference genome
with long alternative haplotypes present in separate ALT contigs. To use the
ALT-aware mode, users need to provide pairwise ALT-to-reference alignment in the
SAM format and rename the file to "*idxbase*.alt". For GRCh38, this alignment
is available from the [binary package of BWA][res].

```sh
# Generate the GRCh38+ALT+decoy+HLA and create the BWA index
wget -O- http://sourceforge.net/projects/bio-bwa/files/bwakit-0.7.11_x64-linux.tar.bz2/download \
  | gzip -dc | tar xf -
bwa.kit/run-gen-hs38d6       # download GRCh38 and write hs38d6.fa
bwa.kit/bwa index hs38d6.fa  # create BWA index
# mapping
bwa.kit/run-bwamem hs38d6.fa read1.fq read2.fq
```

In the final alignment, a read may be placed on the [primary assembly][grcdef]
and multiple overlapping ALT contigs at the same time (on multiple SAM lines).
Mapping quality (mapQ) is properly adjusted by the postprocessing script
`bwa-postalt.js` using the ALT-to-reference alignment `hs38a.fa.alt`. For
details, see the [Methods section](#methods).

## Background

GRCh38 consists of several components: chromosomal assembly, unlocalized contigs
(chromosome known but location unknown), unplaced contigs (chromosome unknown)
and ALT contigs (long clustered variations). The combination of the first three
components is called the *primary assembly*. You can find the more exact
definitions from the [GRC website][grcdef].

GRCh38 ALT contigs are totaled 109Mb in length, spanning 60Mbp genomic regions.
However, sequences that are highly diverged from the primary assembly only
contribute a few million bp. Most subsequences of ALT contigs are nearly
identical to the primary assembly. If we align sequence reads to GRCh38+ALT
treating ALT equal to the primary assembly, we will get many reads with zero
mapping quality and lose variants on them. It is crucial to make the mapper
aware of ALTs.

BWA-MEM is designed to minimize the interference of ALT contigs such that on the
primary assembly, the ALT-aware alignment is highly similar to the alignment
without using ALT contigs in the index. This design choice makes it almost
always safe to map reads to GRCh38+ALT. Although we don't know yet how much
variations on ALT contigs contribute to phenotypes, we would not get the answer
without mapping large cohorts to these extra sequences. We hope our current
implementation encourages researchers to use ALT contigs soon and often.

## Methods

### Sequence alignment

As of now, ALT mapping is done in two separate steps: BWA-MEM mapping and
postprocessing.

#### Step 1: BWA-MEM mapping

At this step, BWA-MEM reads the ALT contig names from "*idxbase*.alt", ignoring
the ALT-to-ref alignment, and labels a potential hit as *ALT* or *non-ALT*,
depending on whether the hit lands on an ALT contig or not. BWA-MEM then reports
alignments and assigns mapQ following these two rules:

1. The mapQ of a non-ALT hit is computed across non-ALT hits only. The mapQ of
   an ALT hit is computed across all hits.

2. If there are no non-ALT hits, the best ALT hit is outputted as the primary
   alignment. If there are both ALT and non-ALT hits, non-ALT hits will be
   primary. ALT hits are reported as supplementary alignments (flag 0x800) only
   if they are better than all overlapping non-ALT hits.

In theory, non-ALT alignments from step 1 should be identical to alignments
against a reference genome with ALT contigs. In practice, the two types of
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
The sequences are included in `hs38d4-extra.fa` from the [BWA resource bundle
for GRCh38][res].

In addition to decoy, we also put multiple alleles of HLA genes in
`hs38d4-extra.fa`. These genomic sequences were acquired from [IMGT/HLA][hladb],
version 3.18.0. Script `bwa-postalt.js` also helps to genotype HLA genes, though
not to high resolution for now.

### More on HLA typing

It is [well known][hlalink] that HLA genes are associated with many autoimmunity
infectious diseases and drug responses. However, many HLA alleles are highly
diverged from the reference genome. If we map whole-genome shotgun (WGS) reads
to the reference only, many allele-informative will get lost. As a result, the
vast majority of WGS projects have ignored these important genes.

We recommend to include the genomic regions of classical HLA genes in the BWA
index. This way we will be able to get a more complete collection of reads
mapped to HLA. We can then isolate these reads with little computational cost
and type HLA genes with another program, such as [Warren et al (2012)][hla4],
[Liu et al (2013)][hla2], [Bai et al (2014)][hla3], [Dilthey et al (2014)][hla1]
or others from [this list][hlatools].

### Evaluating ALT Mapping

(Coming soon...)

## Problems and Future Development

There are some uncertainties about ALT mappings - we are not sure whether they
help biological discovery and don't know the best way to analyze them. Without
clear demand from downstream analyses, it is very difficult to design the
optimal mapping strategy. The current BWA-MEM method is just a start. If it
turns out to be useful in research, we will probably rewrite bwa-postalt.js in C
for performance; if not, we may make changes. It is also possible that we might
make breakthrough on the representation of multiple genomes, in which case, we
can even get rid of ALT contigs for good.



[res]: https://sourceforge.net/projects/bio-bwa/files/
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
