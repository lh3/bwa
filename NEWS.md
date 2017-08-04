Release 0.7.16 (30 July 2017)
-----------------------------

This release added a couple of minor features and incorporated multiple pull
requests, including:

 * Added option -5, which is useful to some Hi-C pipelines.

 * Fixed an error with samtools sorting (#129). Updated download link for
   GRCh38 (#123). Fixed README MarkDown formatting (#70). Addressed multiple
   issues via a collected pull request #139 by @jmarshall. Avoid malformatted
   SAM header when -R is used with TAB (#84). Output mate CIGAR (#138).

(0.7.16: 30 July 2017, r1180)



Release 0.7.15 (31 May 2016)
----------------------------

Fixed a long existing bug which potentially leads to underestimated insert size
upper bound. This bug should have little effect in practice.

(0.7.15: 31 May 2016, r1140)



Release 0.7.14 (4 May 2016)
---------------------------

In the ALT mapping mode, this release adds the "AH:*" header tag to SQ lines
corresponding to alternate haplotypes.

(0.7.14: 4 May 2016, r1136)



Release 0.7.13 (23 Feburary 2016)
---------------------------------

This release fixes a few minor bugs in the previous version and adds a few
minor features. All BWA algorithms should produce identical output to 0.7.12
when there are no ALT contigs.

Detailed changes:

 * Fixed a bug in "bwa-postalt.js". The old version may produce 0.5% of wrong
   bases for reads mapped to the ALT contigs.

 * Fixed a potential bug in the multithreading mode. It may occur when mapping
   is much faster than file reading, which should almost never happen in
   practice.

 * Changed the download URL of GRCh38.

 * Removed the read overlap mode. It is not working well.

 * Added the ropebwt2 algorithm as an alternative to index large genomes.
   Ropebwt2 is slower than the "bwtsw" algorithm, but it has a permissive
   license. This allows us to create an Apache2-licensed BWA (in the "Apache2"
   branch) for commercial users who are concerned with GPL.

(0.7.13: 23 Feburary 2016, r1126)



Release 0.7.12 (28 December 2014)
---------------------------------

This release fixed a bug in the pair-end mode when ALT contigs are present. It
leads to undercalling in regions overlapping ALT contigs.

(0.7.12: 28 December 2014, r1039)



Release 0.7.11 (23 December, 2014)
----------------------------------

A major change to BWA-MEM is the support of mapping to ALT contigs in addition
to the primary assembly. Part of the ALT mapping strategy is implemented in
BWA-MEM and the rest in a postprocessing script for now. Due to the extra
layer of complexity on generating the reference genome and on the two-step
mapping, we start to provide a wrapper script and precompiled binaries since
this release. The package may be more convenient to some specific use cases.
For general uses, the single BWA binary still works like the old way.

Another major addition to BWA-MEM is HLA typing, which made possible with the
new ALT mapping strategy. Necessary data and programs are included in the
binary release. The wrapper script also optionally performs HLA typing when HLA
genes are included in the reference genome as additional ALT contigs.

Other notable changes to BWA-MEM:

 * Added option `-b` to `bwa index`. This option tunes the batch size used in
   the construction of BWT. It is advised to use large `-b` for huge reference
   sequences such as the BLAST *nt* database.

 * Optimized for PacBio data. This includes a change to scoring based on a
   study done by Aaron Quinlan and a heuristic speedup. Further speedup is
   possible, but needs more careful investigation.

 * Dropped PacBio read-to-read alignment for now. BWA-MEM is good for finding
   the best hit, but is not very sensitive to suboptimal hits. Option `-x pbread`
   is still available, but hidden on the command line. This may be removed in
   future releases.

 * Added a new pre-setting for Oxford Nanopore 2D reads. LAST is still a little
   more sensitive on older bacterial data, but bwa-mem is as good on more
   recent data and is times faster for mapping against mammalian genomes.

 * Added LAST-like seeding. This improves the accuracy for longer reads.

 * Added option `-H` to insert arbitrary header lines.

 * Smarter option `-p`. Given an interleaved FASTQ stream, old bwa-mem identifies
   the 2i-th and (2i+1)-th reads as a read pair. The new verion identifies
   adjacent reads with the same read name as a read pair. It is possible to mix
   single-end and paired-end reads in one FASTQ.

 * Improved parallelization. Old bwa-mem waits for I/O. The new version puts
   I/O on a separate thread. It performs mapping while reading FASTQ and
   writing SAM. This saves significant wall-clock time when reading from
   or writing to a slow Unix pipe.

With the new release, the recommended way to map Illumina reads to GRCh38 is to
use the bwakit binary package:

    bwa.kit/run-gen-ref hs38DH
    bwa.kit/bwa index hs38DH.fa
    bwa.kit/run-bwamem -t8 -H -o out-prefix hs38DH.fa read1.fq.gz read2.fq.gz | sh

Please check bwa.kit/README.md for details and command line options.

(0.7.11: 23 December 2014, r1034)



Release 0.7.10 (13 July, 2014)
------------------------------

Notable changes to BWA-MEM:

 * Fixed a segmentation fault due to an alignment bridging the forward-reverse
   boundary. This is a bug.

 * Use the PacBio heuristic to map contigs to the reference genome. The old
   heuristic evaluates the necessity of full extension for each chain. This may
   not work in long low-complexity regions. The PacBio heuristic performs
   SSE2-SW around each short seed. It works better. Note that the heuristic is
   only applied to long query sequences. For Illumina reads, the output is
   identical to the previous version.

(0.7.10: 13 July 2014, r789)



Release 0.7.9 (19 May, 2014)
----------------------------

This release brings several major changes to BWA-MEM. Notably, BWA-MEM now
formally supports PacBio read-to-reference alignment and experimentally supports
PacBio read-to-read alignment. BWA-MEM also runs faster at a minor cost of
accuracy. The speedup is more significant when GRCh38 is in use. More
specifically:

 * Support PacBio subread-to-reference alignment. Although older BWA-MEM works
   with PacBio data in principle, the resultant alignments are frequently
   fragmented. In this release, we fine tuned existing methods and introduced
   new heuristics to improve PacBio alignment. These changes are not used by
   default. Users need to add option "-x pacbio" to enable the feature.

 * Support PacBio subread-to-subread alignment (EXPERIMENTAL). This feature is
   enabled with option "-x pbread". In this mode, the output only gives the
   overlapping region between a pair of reads without detailed alignment.

 * Output alternative hits in the XA tag if there are not so many of them. This
   is a BWA-backtrack feature.

 * Support mapping to ALT contigs in GRCh38 (EXPERIMENTAL). We provide a script
   to postprocess hits in the XA tag to adjust the mapping quality and generate
   new primary alignments to all overlapping ALT contigs. We would *NOT*
   recommend this feature for production uses.

 * Improved alignments to many short reference sequences. Older BWA-MEM may
   generate an alignment bridging two or more adjacent reference sequences.
   Such alignments are split at a later step as postprocessing. This approach
   is complex and does not always work. This release forbids these alignments
   from the very beginning. BWA-MEM should not produce an alignment bridging
   two or more reference sequences any more.

 * Reduced the maximum seed occurrence from 10000 to 500. Reduced the maximum
   rounds of Smith-Waterman mate rescue from 100 to 50. Added a heuristic to
   lower the mapping quality if a read contains seeds with excessive
   occurrences. These changes make BWA-MEM faster at a minor cost of accuracy
   in highly repetitive regions.

 * Added an option "-Y" to use soft clipping for supplementary alignments.

 * Bugfix: incomplete alignment extension in corner cases.

 * Bugfix: integer overflow when aligning long query sequences.

 * Bugfix: chain score is not computed correctly (almost no practical effect)

 * General code cleanup

 * Added FAQs to README

Changes in BWA-backtrack:

 * Bugfix: a segmentation fault when an alignment stands out of the end of the
   last chromosome.

(0.7.9: 19 May 2014, r783)



Release 0.7.8 (31 March, 2014)
------------------------------

Changes in BWA-MEM:

 * Bugfix: off-diagonal X-dropoff (option -d) not working as intended.
   Short-read alignment is not affected.

 * Bugfix: unnecessarily large bandwidth used during global alignment,
   which reduces the mapping speed by -5% for short reads. Results are not
   affected.

 * Bugfix: when the matching score is not one, paired-end mapping quality is
   inaccurate.

 * When the matching score (option -A) is changed, scale all score-related
   options accordingly unless overridden by users.

 * Allow to specify different gap open (or extension) penalties for deletions
   and insertions separately.

 * Allow to specify the insert size distribution.

 * Better and more detailed debugging information.

With the default setting, 0.7.8 and 0.7.7 gave identical output on one million
100bp read pairs.

(0.7.8: 31 March 2014, r455)



Release 0.7.7 (25 Feburary, 2014)
---------------------------------

This release fixes incorrect MD tags in the BWA-MEM output.

A note about short-read mapping to GRCh38. The new human reference genome
GRCh38 contains 60Mbp program generated alpha repeat arrays, some of which are
hard masked as they cannot be localized. These highly repetitive arrays make
BWA-MEM -50% slower. If you are concerned with the performance of BWA-MEM, you
may consider to use option "-c2000 -m50". On simulated data, this setting helps
the performance at a very minor cost on accuracy. I may consider to change the
default in future releases.

(0.7.7: 25 Feburary 2014, r441)



Release 0.7.6 (31 Januaray, 2014)
---------------------------------

Changes in BWA-MEM:

 * Changed the way mapping quality is estimated. The new method tends to give
   the same alignment a higher mapping quality. On paired-end reads, the change
   is minor as with pairing, the mapping quality is usually high. For short
   single-end reads, the difference is considerable.

 * Improved load balance when many threads are spawned. However, bwa-mem is
   still not very thread efficient, probably due to the frequent heap memory
   allocation. Further improvement is a little difficult and may affect the
   code stability.

 * Allow to use different clipping penalties for 5'- and 3'-ends. This helps
   when we do not want to clip one end.

 * Print the @PG line, including the command line options.

 * Improved the band width estimate: a) fixed a bug causing the band
   width extimated from extension not used in the final global alignment; b)
   try doubled band width if the global alignment score is smaller.
   Insufficient band width leads to wrong CIGAR and spurious mismatches/indels.

 * Added a new option -D to fine tune a heuristic on dropping suboptimal hits.
   Reducing -D increases accuracy but decreases the mapping speed. If unsure,
   leave it to the default.

 * Bugfix: for a repetitive single-end read, the reported hit is not randomly
   distributed among equally best hits.

 * Bugfix: missing paired-end hits due to unsorted list of SE hits.

 * Bugfix: incorrect CIGAR caused by a defect in the global alignment.

 * Bugfix: incorrect CIGAR caused by failed SW rescue.

 * Bugfix: alignments largely mapped to the same position are regarded to be
   distinct from each other, which leads to underestimated mapping quality.

 * Added the MD tag.

There are no changes to BWA-backtrack in this release. However, it has a few
known issues yet to be fixed. If you prefer BWA-track, It is still advised to
use bwa-0.6.x.

While I developed BWA-MEM, I also found a few issues with BWA-SW. It is now
possible to improve BWA-SW with the lessons learned from BWA-MEM. However, as
BWA-MEM is usually better, I will not improve BWA-SW until I find applications
where BWA-SW may excel.

(0.7.6: 31 January 2014, r432)



Release 0.7.5a (30 May, 2013)
-----------------------------

Fixed a bug in BWA-backtrack which leads to off-by-one mapping errors in rare
cases.

(0.7.5a: 30 May 2013, r405)



Release 0.7.5 (29 May, 2013)
----------------------------

Changes in all components:

 * Improved error checking on memory allocation and file I/O. Patches provided
   by Rob Davies.

 * Updated README.

 * Bugfix: return code is zero upon errors.

Changes in BWA-MEM:

 * Changed the way a chimeric alignment is reported (conforming to the upcoming
   SAM spec v1.5). With 0.7.5, if the read has a chimeric alignment, the paired
   or the top hit uses soft clipping and is marked with neither 0x800 nor 0x100
   bits. All the other hits part of the chimeric alignment will use hard
   clipping and be marked with 0x800 if option "-M" is not in use, or marked
   with 0x100 otherwise.

 * Other hits part of a chimeric alignment are now reported in the SA tag,
   conforming to the SAM spec v1.5.

 * Better method for resolving an alignment bridging two or more short
   reference sequences. The current strategy maps the query to the reference
   sequence that covers the middle point of the alignment. For most
   applications, this change has no effects.

Changes in BWA-backtrack:

 * Added a magic number to .sai files. This prevents samse/sampe from reading
   corrupted .sai (e.g. a .sai file containing LSF log) or incompatible .sai
   generated by a different version of bwa.

 * Bugfix: alignments in the XA:Z: tag were wrong.

 * Keep track of #ins and #del during backtracking. This simplifies the code
   and reduces errors in rare corner cases. I should have done this in the
   early days of bwa.

In addition, if you use BWA-MEM or the fastmap command of BWA, please cite:

 - Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs
   with BWA-MEM. arXiv:1303.3997v2 [q-bio.GN].

Thank you.

(0.7.5: 29 May 2013, r404)



Release 0.7.4 (23 April, 2013)
------------------------------

This is a bugfix release. Most of bugs are considered to be minor which only
occur very rarely.

 * Bugfix: wrong CIGAR when a query sequence bridges three or more target
   sequences. This only happens when aligning reads to short assembly contigs.

 * Bugfix: leading "D" operator in CIGAR.

 * Extend more seeds for better alignment around tandem repeats. This is also
   a cause of the leading "D" operator in CIGAR.

 * Bugfix: SSE2-SSW may occasionally find incorrect query starting position
   around tandem repeat. This will lead to a suboptimal CIGAR in BWA-MEM and
   a wrong CIGAR in BWA.

 * Bugfix: clipping penalty does not work as is intended when there is a gap
   towards the end of a read.

 * Fixed an issue caused by a bug in the libc from Mac/Darwin. In Darwin,
   fread() is unable to read a data block longer than 2GB due to an integer
   overflow bug in its implementation.

Since version 0.7.4, BWA-MEM is considered to reach similar stability to
BWA-backtrack for short-read mapping.

(0.7.4: 23 April, r385)



Release 0.7.3a (15 March, 2013)
-------------------------------

In 0.7.3, the wrong CIGAR bug was only fixed in one scenario, but not fixed
in another corner case.

(0.7.3a: 15 March 2013, r367)



Release 0.7.3 (15 March, 2013)
------------------------------

Changes to BWA-MEM:

 * Bugfix: pairing score is inaccurate when option -A does not take the default
   value. This is a very minor issue even if it happens.

 * Bugfix: occasionally wrong CIGAR. This happens when in the alignment there
   is a 1bp deletion and a 1bp insertion which are close to the end of the
   reads, and there are no other substitutions or indels. BWA-MEM would not do
   a gapped alignment due to the bug.

 * New feature: output other non-overlapping alignments in the XP tag such that
   we can see the entire picture of alignment from one SAM line. XP gives the
   position, CIGAR, NM and mapQ of each aligned subsequence of the query.

BWA-MEM has been used to align -300Gbp 100-700bp SE/PE reads. SNP/indel calling
has also been evaluated on part of these data. BWA-MEM generally gives better
pre-filtered SNP calls than BWA. No significant issues have been observed since
0.7.2, though minor improvements or bugs (e.g. the bug fixed in this release)
are still possible. If you find potential issues, please send bug reports to
<bio-bwa-help@lists.sourceforge.net> (free registration required).

In addition, more detailed description of the BWA-MEM algorithm can be found at
<https://github.com/lh3/mem-paper>.

(0.7.3: 15 March 2013, r366)



Release 0.7.2 (9 March, 2013)
-----------------------------

Emergent bug fix: 0.7.0 and 0.7.1 give a wrong sign to TLEN. In addition,
flagging 'properly paired' also gets improved a little.

(0.7.2: 9 March 2013, r351)



Release 0.7.1 (8 March, 2013)
-----------------------------

Changes to BWA-MEM:

 * Bugfix: rare segmentation fault caused by a partial hit to the end of the
   last sequence.

 * Bugfix: occasional mis-pairing given an interleaved fastq.

 * Bugfix: wrong mate information when the mate is unmapped. SAM generated by
   BWA-MEM can now be validated with Picard.

 * Improved the performance and accuracy for ultra-long query sequences.
   Short-read alignment is not affected.

Changes to other components:

 * In BWA-backtrack and BWA-SW, replaced the code for global alignment,
   Smith-Waterman and SW extension. The performance and accuracy of the two
   algorithms stay the same.

 * Added an experimental subcommand to merge overlapping paired ends. The
   algorithm is very conservative: it may miss true overlaps but rarely makes
   mistakes.

An important note is that like BWA-SW, BWA-MEM may output multiple primary
alignments for a read, which may cause problems to some tools. For aligning
sequence reads, it is advised to use '-M' to flag extra hits as secondary. This
option is not the default because multiple primary alignments are theoretically
possible in sequence alignment.

(0.7.1: 8 March 2013, r347)



Beta Release 0.7.0 (28 Feburary, 2013)
--------------------------------------

This release comes with a new alignment algorithm, BWA-MEM, for 70bp-1Mbp query
sequences. BWA-MEM essentially seeds alignments with a variant of the fastmap
algorithm and extends seeds with banded affine-gap-penalty dynamic programming
(i.e. the Smith-Waterman-Gotoh algorithm). For typical Illumina 100bp reads or
longer low-divergence query sequences, BWA-MEM is about twice as fast as BWA
and BWA-SW and is more accurate. It also supports split alignments like BWA-SW
and may optionally output multiple hits like BWA. BWA-MEM does not guarantee
to find hits within a certain edit distance, but BWA is not efficient for such
task given longer reads anyway, and the edit-distance criterion is arguably
not as important in long-read alignment.

In addition to the algorithmic improvements, BWA-MEM also implements a few
handy features in practical aspects:

 1. BWA-MEM automatically switches between local and glocal (global wrt reads;
    local wrt reference) alignment. It reports the end-to-end glocal alignment
    if the glocal alignment is not much worse than the optimal local alignment.
    Glocal alignment reduces reference bias.

 2. BWA-MEM automatically infers pair orientation from a batch of single-end
    alignments. It allows more than one orientations if there are sufficient
    supporting reads. This feature has not been tested on reads from Illumina
    jumping library yet. (EXPERIMENTAL)

 3. BWA-MEM optionally takes one interleaved fastq for paired-end mapping. It
    is possible to convert a name-sorted BAM to an interleaved fastq on the fly
    and feed the data stream to BWA-MEM for mapping.

 4. BWA-MEM optionally copies FASTA/Q comments to the final SAM output, which
    helps to transfer individual read annotations to the output.

 5. BWA-MEM supports more advanced piping. Users can now run:
    (bwa mem ref.fa '<bzcat r1.fq.bz2' '<bzcat r2.fq.bz2') to map bzip'd read
    files without replying on bash features.

 6. BWA-MEM provides a few basic APIs for single-end mapping. The 'example.c'
    program in the source code directory implements a full single-end mapper in
    50 lines of code.

The BWA-MEM algorithm is in the beta phase. It is not advised to use BWA-MEM
for production use yet. However, when the implementation becomes stable after a
few release cycles, existing BWA users are recommended to migrate to BWA-MEM
for 76bp or longer Illumina reads and long query sequences. The original BWA
short-read algorithm will not deliver satisfactory results for 150bp+ Illumina
reads. Change of mappers will be necessary sooner or later.

(0.7.0 beta: 28 Feburary 2013, r313)



Release 0.6.2 (19 June, 2012)
-----------------------------

This is largely a bug-fix release. Notable changes in BWA-short and BWA-SW:

 * Bugfix: BWA-SW may give bad alignments due to incorrect band width.

 * Bugfix: A segmentation fault due to an out-of-boundary error. The fix is a
   temporary solution. The real cause has not been identified.

 * Attempt to read index from prefix.64.bwt, such that the 32-bit and 64-bit
   index can coexist.

 * Added options '-I' and '-S' to control BWA-SW pairing.

(0.6.2: 19 June 2012, r126)



Release 0.6.1 (28 November, 2011)
---------------------------------

Notable changes to BWA-short:

 * Bugfix: duplicated alternative hits in the XA tag.

 * Bugfix: when trimming enabled, bwa-aln trims 1bp less.

 * Disabled the color-space alignment. 0.6.x is not working with SOLiD reads at
   present.

Notable changes to BWA-SW:

 * Bugfix: segfault due to excessive ambiguous bases.

 * Bugfix: incorrect mate position in the SE mode.

 * Bugfix: rare segfault in the PE mode

 * When macro _NO_SSE2 is in use, fall back to the standard Smith-Waterman
   instead of SSE2-SW.

 * Optionally mark split hits with lower alignment scores as secondary.

Changes to fastmap:

 * Bugfix: infinite loop caused by ambiguous bases.

 * Optionally output the query sequence.

(0.6.1: 28 November 2011, r104)



Release 0.5.10 and 0.6.0 (12 November, 2011)
--------------------------------------------

The 0.6.0 release comes with two major changes. Firstly, the index data
structure has been changed to support genomes longer than 4GB. The forward and
reverse backward genome is now integrated in one index. This change speeds up
BWA-short by about 20% and BWA-SW by 90% with the mapping acccuracy largely
unchanged. A tradeoff is BWA requires more memory, but this is the price almost
all mappers that index the genome have to pay.

Secondly, BWA-SW in 0.6.0 now works with paired-end data. It is more accurate
for highly unique reads and more robust to long indels and structural
variations. However, BWA-short still has edges for reads with many suboptimal
hits. It is yet to know which algorithm is the best for variant calling.

0.5.10 is a bugfix release only and is likely to be the last release in the 0.5
branch unless I find critical bugs in future.

Other notable changes:

 * Added the 'fastmap' command that finds super-maximal exact matches. It does
   not give the final alignment, but runs much faster. It can be a building
   block for other alignment algorithms. [0.6.0 only]

 * Output the timing information before BWA exits. This also tells users that
   the task has been finished instead of being killed or aborted. [0.6.0 only]

 * Sped up multi-threading when using many (>20) CPU cores.

 * Check I/O error.

 * Increased the maximum barcode length to 63bp.

 * Automatically choose the indexing algorithm.

 * Bugfix: very rare segfault due to an uninitialized variable. The bug also
   affects the placement of suboptimal alignments. The effect is very minor.

This release involves quite a lot of tricky changes. Although it has been
tested on a few data sets, subtle bugs may be still hidden. It is *NOT*
recommended to use this release in a production pipeline. In future, however,
BWA-SW may be better when reads continue to go longer. I would encourage users
to try the 0.6 release. I would also like to hear the users' experience. Thank
you.

(0.6.0: 12 November 2011, r85)



Beta Release 0.5.9 (24 January, 2011)
-------------------------------------

Notable changes:

 * Feature: barcode support via the '-B' option.

 * Feature: Illumina 1.3+ read format support via the '-I' option.

 * Bugfix: RG tags are not attached to unmapped reads.

 * Bugfix: very rare bwasw mismappings

 * Recommend options for PacBio reads in bwasw help message.


Also, since January 13, the BWA master repository has been moved to github:

  https://github.com/lh3/bwa

The revision number has been reset. All recent changes will be first
committed to this repository.

(0.5.9: 24 January 2011, r16)



Beta Release Candidate 0.5.9rc1 (10 December, 2010)
---------------------------------------------------

Notable changes in bwasw:

 * Output unmapped reads.

 * For a repetitive read, choose a random hit instead of a fixed
   one. This is not well tested.

Notable changes in bwa-short:

 * Fixed a bug in the SW scoring system, which may lead to unexpected
   gaps towards the end of a read.

 * Fixed a bug which invalidates the randomness of repetitive reads.

 * Fixed a rare memory leak.

 * Allowed to specify the read group at the command line.

 * Take name-grouped BAM files as input.

Changes to this release are usually safe in that they do not interfere
with the key functionality. However, the release has only been tested on
small samples instead of on large-scale real data. If anything weird
happens, please report the bugs to the bio-bwa-help mailing list.

(0.5.9rc1: 10 December 2010, r1561)



Beta Release 0.5.8 (8 June, 2010)
---------------------------------

Notable changes in bwasw:

 * Fixed an issue of missing alignments. This should happen rarely and
   only when the contig/read alignment is multi-part. Very rarely, bwasw
   may still miss a segment in a multi-part alignment. This is difficult
   to fix, although possible.

Notable changes in bwa-short:

 * Discard the SW alignment when the best single-end alignment is much
   better. Such a SW alignment may caused by structural variations and
   forcing it to be aligned leads to false alignment. This fix has not
   been tested thoroughly. It would be great to receive more users
   feedbacks on this issue.

 * Fixed a typo/bug in sampe which leads to unnecessarily large memory
   usage in some cases.

 * Further reduced the chance of reporting 'weird pairing'.

(0.5.8: 8 June 2010, r1442)



Beta Release 0.5.7 (1 March, 2010)
----------------------------------

This release only has an effect on paired-end data with fat insert-size
distribution. Users are still recommended to update as the new release
improves the robustness to poor data.

 * The fix for 'weird pairing' was not working in version 0.5.6, pointed
   out by Carol Scott. It should work now.

 * Optionally output to a normal file rather than to stdout (by Tim
   Fennel).

(0.5.7: 1 March 2010, r1310)



Beta Release 0.5.6 (10 Feburary, 2010)
--------------------------------------

Notable changes in bwa-short:

 * Report multiple hits in the SAM format at a new tag XA encoded as:
   (chr,pos,CIGAR,NM;)*. By default, if a paired or single-end read has
   4 or fewer hits, they will all be reported; if a read in a anomalous
   pair has 11 or fewer hits, all of them will be reported.

 * Perform Smith-Waterman alignment also for anomalous read pairs when
   both ends have quality higher than 17. This reduces false positives
   for some SV discovery algorithms.

 * Do not report "weird pairing" when the insert size distribution is
   too fat or has a mean close to zero.

 * If a read is bridging two adjacent chromsomes, flag it as unmapped.

 * Fixed a small but long existing memory leak in paired-end mapping.

 * Multiple bug fixes in SOLiD mapping: a) quality "-1" can be correctly
   parsed by solid2fastq.pl; b) truncated quality string is resolved; c)
   SOLiD read mapped to the reverse strand is complemented.

 * Bwa now calculates skewness and kurtosis of the insert size
   distribution.

 * Deploy a Bayesian method to estimate the maximum distance for a read
   pair considered to be paired properly. The method is proposed by
   Gerton Lunter, but bwa only implements a simplified version.

 * Export more functions for Java bindings, by Matt Hanna (See:
   http://www.broadinstitute.org/gsa/wiki/index.php/Sting_BWA/C_bindings)

 * Abstract bwa CIGAR for further extension, by Rodrigo Goya.

(0.5.6: 10 Feburary 2010, r1303)



Beta Release 0.5.5 (10 November, 2009)
--------------------------------------

This is a bug fix release:

 * Fixed a serious bug/typo in aln which does not occur given short
   reads, but will lead to segfault for >500bp reads. Of course, the aln
   command is not recommended for reads longer than 200bp, but this is a
   bug anyway.

 * Fixed a minor bug/typo which leads to incorrect single-end mapping
   quality when one end is moved to meet the mate-pair requirement.

 * Fixed a bug in samse for mapping in the color space. This bug is
   caused by quality filtration added since 0.5.1.

(0.5.5: 10 November 2009, r1273)



Beta Release 0.5.4 (9 October, 2009)
------------------------------------

Since this version, the default seed length used in the "aln" command is
changed to 32.

Notable changes in bwa-short:

 * Added a new tag "XC:i" which gives the length of clipped reads.

 * In sampe, skip alignments in case of a bug in the Smith-Waterman
   alignment module.

 * In sampe, fixed a bug in pairing when the read sequence is identical
   to its reverse complement.

 * In sampe, optionally preload the entire FM-index into memory to
   reduce disk operations.

Notable changes in dBWT-SW/BWA-SW:

 * Changed name dBWT-SW to BWA-SW.

 * Optionally use "hard clipping" in the SAM output.

(0.5.4: 9 October 2009, r1245)



Beta Release 0.5.3 (15 September, 2009)
---------------------------------------

Fixed a critical bug in bwa-short: reads mapped to the reverse strand
are not complemented.

(0.5.3: 15 September 2009, r1225)



Beta Release 0.5.2 (13 September, 2009)
---------------------------------------

Notable changes in bwa-short:

 * Optionally trim reads before alignment. See the manual page on 'aln
   -q' for detailed description.

 * Fixed a bug in calculating the NM tag for a gapped alignment.

 * Fixed a bug given a mixture of reads with some longer than the seed
   length and some shorter.

 * Print SAM header.

Notable changes in dBWT-SW:

 * Changed the default value of -T to 30. As a result, the accuracy is a
   little higher for short reads at the cost of speed.

(0.5.2: 13 September 2009, r1223)



Beta Release 0.5.1 (2 September, 2009)
--------------------------------------

Notable changes in the short read alignment component:

 * Fixed a bug in samse: do not write mate coordinates.

Notable changes in dBWT-SW:

 * Randomly choose one alignment if the read is a repetitive.

 * Fixed a flaw when a read is mapped across two adjacent reference
   sequences. However, wrong alignment reports may still occur rarely in
   this case.

 * Changed the default band width to 50. The speed is slower due to this
   change.

 * Improved the mapping quality a little given long query sequences.

(0.5.1: 2 September 2009, r1209)



Beta Release 0.5.0 (20 August, 2009)
------------------------------------

This release implements a novel algorithm, dBWT-SW, specifically
designed for long reads. It is 10-50 times faster than SSAHA2, depending
on the characteristics of the input data, and achieves comparable
alignment accuracy while allowing chimera detection. In comparison to
BLAT, dBWT-SW is several times faster and much more accurate especially
when the error rate is high. Please read the manual page for more
information.

The dBWT-SW algorithm is kind of developed for future sequencing
technologies which produce much longer reads with a little higher error
rate. It is still at its early development stage. Some features are
missing and it may be buggy although I have evaluated on several
simulated and real data sets. But following the "release early"
paradigm, I would like the users to try it first.

Other notable changes in BWA are:

 * Fixed a rare bug in the Smith-Waterman alignment module.

 * Fixed a rare bug about the wrong alignment coordinate when a read is
   poorly aligned.

 * Fixed a bug in generating the "mate-unmap" SAM tag when both ends in
   a pair are unmapped.

(0.5.0: 20 August 2009, r1200)



Beta Release 0.4.9 (19 May, 2009)
---------------------------------

Interestingly, the integer overflow bug claimed to be fixed in 0.4.7 has
not in fact. Now I have fixed the bug. Sorry for this and thank Quan
Long for pointing out the bug (again).

(0.4.9: 19 May 2009, r1075)



Beta Release 0.4.8 (18 May, 2009)
---------------------------------

One change to "aln -R". Now by default, if there are no more than '-R'
equally best hits, bwa will search for suboptimal hits. This change
affects the ability in finding SNPs in segmental duplications.

I have not tested this option thoroughly, but this simple change is less
likely to cause new bugs. Hope I am right.

(0.4.8: 18 May 2009, r1073)



Beta Release 0.4.7 (12 May, 2009)
---------------------------------

Notable changes:

 * Output SM (single-end mapping quality) and AM (smaller mapping
   quality among the two ends) tag from sam output.

 * Improved the functionality of stdsw.

 * Made the XN tag more accurate.

 * Fixed a very rare segfault caused by integer overflow.

 * Improve the insert size estimation.

 * Fixed compiling errors for some Linux systems.

(0.4.7: 12 May 2009, r1066)



Beta Release 0.4.6 (9 March, 2009)
----------------------------------

This release improves the SOLiD support. First, a script for converting
SOLiD raw data is provided. This script is adapted from solid2fastq.pl
in the MAQ package. Second, a nucleotide reference file can be directly
used with 'bwa index'. Third, SOLiD paired-end support is
completed. Fourth, color-space reads will be converted to nucleotides
when SAM output is generated. Color errors are corrected in this
process. Please note that like MAQ, BWA cannot make use of the primer
base and the first color.

In addition, the calculation of mapping quality is also improved a
little bit, although end-users may barely observe the difference.

(0.4.6: 9 March 2009, r915)



Beta Release 0.4.5 (18 Feburary, 2009)
--------------------------------------

Not much happened, but I think it would be good to let the users use the
latest version.

Notable changes (Thank Bob Handsaker for catching the two bugs):

 * Improved bounary check. Previous version may still give incorrect
   alignment coordinates in rare cases.

 * Fixed a bug in SW alignment when no residue matches. This only
   affects the 'sampe' command.

 * Robustly estimate insert size without setting the maximum on the
   command line. Since this release 'sampe -a' only has an effect if
   there are not enough good pairs to infer the insert size
   distribution.

 * Reduced false PE alignments a little bit by using the inferred insert
   size distribution. This fix may be more important for long insert
   size libraries.

(0.4.5: 18 Feburary 2009, r829)



Beta Release 0.4.4 (15 Feburary, 2009)
--------------------------------------

This is mainly a bug fix release. Notable changes are:

 * Imposed boundary check for extracting subsequence from the
   genome. Previously this causes memory problem in rare cases.

 * Fixed a bug in failing to find whether an alignment overlapping with
   N on the genome.

 * Changed MD tag to meet the latest SAM specification.

(0.4.4: 15 Feburary 2009, r815)



Beta Release 0.4.3 (22 January, 2009)
------------------------------------

Notable changes:

 * Treat an ambiguous base N as a mismatch. Previous versions will not
   map reads containing any N.

 * Automatically choose the maximum allowed number of differences. This
   is important when reads of different lengths are mixed together.

 * Print mate coordinate if only one end is unmapped.

 * Generate MD tag. This tag encodes the mismatching positions and the
   reference bases at these positions. Deletions from the reference will
   also be printed.

 * Optionally dump multiple hits from samse, in another concise format
   rather than SAM.

 * Optionally disable iterative search. This is VERY SLOOOOW, though.

 * Fixed a bug in generate SAM.

(0.4.3: 22 January 2009, r787)



Beta Release 0.4.2 (9 January, 2009)
------------------------------------

Aaron Quinlan found a bug in the indexer: the bwa indexer segfaults if
there are no comment texts in the FASTA header. This is a critical
bug. Nothing else was changed.

(0.4.2: 9 January 2009, r769)



Beta Release 0.4.1 (7 January, 2009)
------------------------------------

I am sorry for the quick updates these days. I like to set a milestone
for BWA and this release seems to be. For paired end reads, BWA also
does Smith-Waterman alignment for an unmapped read whose mate can be
mapped confidently. With this strategy BWA achieves similar accuracy to
maq. Benchmark is also updated accordingly.

(0.4.1: 7 January 2009, r760)



Beta Release 0.4.0 (6 January, 2009)
------------------------------------

In comparison to the release two days ago, this release is mainly tuned
for performance with some tricks I learnt from Bowtie. However, as the
indexing format has also been changed, I have to increase the version
number to 0.4.0 to emphasize that *DATABASE MUST BE RE-INDEXED* with
'bwa index'.

 * Improved the speed by about 20%.

 * Added multi-threading to 'bwa aln'.

(0.4.0: 6 January 2009, r756)



Beta Release 0.3.0 (4 January, 2009)
------------------------------------

 * Added paired-end support by separating SA calculation and alignment
   output.

 * Added SAM output.

 * Added evaluation to the documentation.

(0.3.0: 4 January 2009, r741)



Beta Release 0.2.0 (15 Augusst, 2008)
-------------------------------------

 * Take the subsequence at the 5'-end as seed. Seeding strategy greatly
   improves the speed for long reads, at the cost of missing a few true
   hits that contain many differences in the seed. Seeding also increase
   the memory by 800MB.

 * Fixed a bug which may miss some gapped alignments. Fixing the bug
   also slows the speed a little.

(0.2.0: 15 August 2008, r428)



Beta Release 0.1.6 (08 Augusst, 2008)
-------------------------------------

 * Give accurate CIGAR string.

 * Add a simple interface to SW/NW alignment

(0.1.6: 08 August 2008, r414)



Beta Release 0.1.5 (27 July, 2008)
----------------------------------

 * Improve the speed. This version is expected to give the same results.

(0.1.5: 27 July 2008, r400)



Beta Release 0.1.4 (22 July, 2008)
----------------------------------

 * Fixed a bug which may cause missing gapped alignments.

 * More clearly define what alignments can be found by BWA (See
   manual). Now BWA runs a little slower because it will visit more
   potential gapped alignments.

 * A bit code clean up.

(0.1.4: 22 July 2008, r387)



Beta Release 0.1.3 (21 July, 2008)
----------------------------------

Improve the speed with some tricks on retrieving occurences. The results
should be exactly the same as that of 0.1.2.

(0.1.3: 21 July 2008, r382)



Beta Release 0.1.2 (17 July, 2008)
----------------------------------

Support gapped alignment. Codes for ungapped alignment has been removed.

(0.1.2: 17 July 2008, r371)



Beta Release 0.1.1 (03 June, 2008)
-----------------------------------

This is the first release of BWA, Burrows-Wheeler Alignment tool. Please
read man page for more information about this software.

(0.1.1: 03 June 2008, r349)
