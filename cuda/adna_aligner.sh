#!/bin/bash
#
#            *** adna_aligner.sh -- aDNA mapping with aln/mem read-length split ***
#
# Combines the DNAharvester (NBISweden) "bwa-aln-mem" idea with the local
# aDNApipeline scripts:
#
#   1. AdapterRemoval3  (merge overlapping pairs, --preserve5p, --min-length 30)
#   2. Split the processed reads at a length threshold (default 70 bp):
#        - SHORT reads  (< threshold)  -> bwa aln (-l 1024 -n 0.01 -o 2) + samse,
#                                         or GPU bwa gpualn with -g (bit-exact, ~5x faster)
#        - LONG  reads  (>= threshold) -> one of:  bwa mem | bwa-mem3 mem | parabricks fq2bam
#   3. samtools merge the two BAMs back into one
#   4. (optional) sambamba markdup
#
# Short reads always use bwa aln (BWA-backtrack), the correct aligner for
# damaged/short ancient inserts. Only the LONG-read aligner is selectable.
#
# Single-end and paired-end inputs are both supported. After adapter removal
# everything is treated as single-end (merged/collapsed reads), exactly as in
# paired_mem3.sh -- so no -p / interleaving is used.
#
# ---------------------------------------------------------------------------
# Usage:
#   ./adna_aligner.sh -1 R1.fq.gz [-2 R2.fq.gz] -i SAMPLE [options]
#
#   -1 FILE   R1 fastq (or the single-end fastq)            [required]
#   -2 FILE   R2 fastq -> enables paired-end mode           [optional]
#   -i NAME   individual/sample name (= output directory)   [required]
#   -p NAME   population name (for the log only)            [default: unknown]
#   -t INT    threads                                       [default: 8]
#   -a TOOL   long-read aligner: mem | mem3 | fq2bam        [default: mem]
#   -s VAL    short/long split length in bp, or 'auto'      [default: 64]
#   -r FILE   reference fasta                               [default: hs37d5.fa]
#   -O STR    long-read bwa mem / mem3 options              [default: "-k 19 -r 2.5 -L 15"]
#   -g        map SHORT reads on the GPU (bwa gpualn) instead of CPU bwa aln
#   -A        skip the aln step: map ALL reads with the -a aligner (quick first pass)
#   -d        run duplicate marking (sambamba markdup) on the merged BAM
#   -h        show this help
#
# Examples:
#   ./adna_aligner.sh -1 ind.fq.gz            -i IND1 -t 16
#   ./adna_aligner.sh -1 R1.fq.gz -2 R2.fq.gz -i IND2 -t 16 -a mem3 -d
#   ./adna_aligner.sh -1 R1.fq.gz -2 R2.fq.gz -i IND3 -a fq2bam -s 80
# ---------------------------------------------------------------------------

set -o pipefail

# --- Tunable external tools (override via environment if needed) ------------
BWA_MEM3="${BWA_MEM3:-/home/teemu/sorsa/bwa-mem3-main/bwa-mem3}"
PARABRICKS_IMAGE="${PARABRICKS_IMAGE:-nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1}"
BWA_GPU="${BWA_GPU:-bwa-gpu}"     # CUDA-enabled bwa with `gpualn` (used when -g is given)

# --- Defaults ---------------------------------------------------------------
FASTQ1=""
FASTQ2=""
INDNAME=""
POPNAME="unknown"
THREADS=8
ALIGNER="mem"
SPLIT=64                          # fixed default cutoff (bp); or 'auto' to pick from the AdapterRemoval report (floored at 60)
REF="hs37d5.fa"
MEM_OPTS="-k 19 -r 2.5 -L 15"     # ancient bwa mem / mem3 params (DNAharvester defaults)
DEDUP=0
GPU_ALN=0                         # -g : map SHORT reads on the GPU (bwa gpualn) instead of CPU bwa aln
SKIP_ALN=0                        # -A : skip the short-read aln step; map ALL reads with the -a aligner

usage() {
    cat <<'EOF'
adna_aligner.sh -- aDNA mapping with a bwa aln / mem read-length split

Usage:
  ./adna_aligner.sh -1 R1.fq.gz [-2 R2.fq.gz] -i SAMPLE [options]

  -1 FILE   R1 fastq (or the single-end fastq)            [required]
  -2 FILE   R2 fastq -> enables paired-end mode           [optional]
  -i NAME   individual/sample name (= output directory)   [required]
  -p NAME   population name (for the log only)            [default: unknown]
  -t INT    threads                                       [default: 8]
  -a TOOL   long-read aligner: mem | mem3 | fq2bam        [default: mem]
  -s VAL    short/long split length: an INT in bp, or 'auto' [default: 64]
  -r FILE   reference fasta                               [default: hs37d5.fa]
  -O STR    long-read bwa mem / mem3 options              [default: "-k 19 -r 2.5 -L 15"]
  -g        map SHORT reads on the GPU (bwa gpualn) instead of CPU bwa aln [CPU]
  -A        skip the short-read aln step entirely: map ALL reads with the -a
            aligner (mem/mem3/fq2bam). Fast first pass for a new sample with no
            prior data; ignores -s and -g (no reads go to aln/gpualn).
  -d        run duplicate marking (sambamba markdup) on the merged BAM
  -h        show this help

Short reads (< split) map with BWA-backtrack (-l 1024 -n 0.01 -o 2): on the CPU
(bwa aln + samse) by default, or on the GPU with -g (bwa gpualn, bit-exact, ~5x
faster on an RTX 3090; needs the CUDA-enabled bwa-gpu, override path via BWA_GPU).

Split length (-s): default is a FIXED 64 bp. This keeps the precision-critical 30-60 bp
band (the bulk of aDNA, where aln is most precise and mem -k19 -r2.5 loses precision and
inflates reference bias -- Oliva et al.) on backtracking 'aln', and sits at aln's cost
cliff (reads whose max_diff would reach 5, ~64 bp at -n 0.01). Pass 'auto' to instead pick
the cutoff from AdapterRemoval3's post-trim length report ({sample}.json) via a tree-cost
budget (override PICK_SPLIT_BUDGET, default 2e14); auto is CLAMPED to a floor of 60 so the
runtime guard can never push the precision-critical band onto mem. Or pass any integer
(e.g. -s 70) to force a fixed cutoff (70 leans toward precision at higher cost).
The long-read aligner is selected with -a. Note: parabricks fq2bam only
accepts a subset of bwa mem flags (-M -Y -C -T -B -U -L -I -K), so the
ancient '-k 19 -r 2.5' seeding params are dropped for fq2bam.

Examples:
  ./adna_aligner.sh -1 ind.fq.gz            -i IND1 -t 16
  ./adna_aligner.sh -1 R1.fq.gz -2 R2.fq.gz -i IND2 -t 16 -a mem3 -d
  ./adna_aligner.sh -1 R1.fq.gz -2 R2.fq.gz -i IND3 -a fq2bam -s 80
EOF
    exit "${1:-0}"
}

# --- Parse options ----------------------------------------------------------
while getopts ":1:2:i:p:t:a:s:r:O:gAdh" opt; do
    case "$opt" in
        1) FASTQ1="$OPTARG" ;;
        2) FASTQ2="$OPTARG" ;;
        i) INDNAME="$OPTARG" ;;
        p) POPNAME="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        a) ALIGNER="$OPTARG" ;;
        s) SPLIT="$OPTARG" ;;
        r) REF="$OPTARG" ;;
        O) MEM_OPTS="$OPTARG" ;;
        g) GPU_ALN=1 ;;
        A) SKIP_ALN=1 ;;
        d) DEDUP=1 ;;
        h) usage 0 ;;
        :) echo "ERROR: -$OPTARG requires an argument" >&2; usage 1 ;;
        \?) echo "ERROR: unknown option -$OPTARG" >&2; usage 1 ;;
    esac
done

# --- Validate ---------------------------------------------------------------
[ -z "$FASTQ1" ] && { echo "ERROR: -1 <fastq> is required" >&2; usage 1; }
[ -z "$INDNAME" ] && { echo "ERROR: -i <sample name> is required" >&2; usage 1; }
case "$ALIGNER" in
    mem|mem3|fq2bam) ;;
    *) echo "ERROR: -a must be one of: mem, mem3, fq2bam (got '$ALIGNER')" >&2; exit 1 ;;
esac
case "$SPLIT" in
    auto) ;;
    ''|*[!0-9]*) echo "ERROR: -s must be a positive integer or 'auto' (got '$SPLIT')" >&2; exit 1 ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ "$GPU_ALN" -eq 1 ] && ! command -v "$BWA_GPU" >/dev/null 2>&1; then
    echo "ERROR: -g given but '$BWA_GPU' not found on PATH (build with 'make bwa-gpu' or set \$BWA_GPU)" >&2
    exit 1
fi

if [ -n "$FASTQ2" ]; then MODE="paired"; else MODE="single"; fi

SHORT_ALN_DESC=$([ "$GPU_ALN" -eq 1 ] && echo "GPU bwa gpualn" || echo "CPU bwa aln")

REF_PREFIX="${REF%.fa}"; REF_PREFIX="${REF_PREFIX%.fasta}"; REF_PREFIX="${REF_PREFIX%.fna}"

echo ""
echo "            *** adna_aligner.sh ***"
echo "  Sample      : $INDNAME   (population: $POPNAME)"
echo "  Mode        : $MODE-end"
echo "  Reference   : $REF"
if [ "$SKIP_ALN" -eq 1 ]; then
    echo "  Split length: (none -- aln skipped with -A; ALL reads -> $ALIGNER)"
elif [ "$SPLIT" = "auto" ]; then
    echo "  Split length: auto   (chosen from the AdapterRemoval length report; short -> $SHORT_ALN_DESC, long -> $ALIGNER)"
else
    echo "  Split length: $SPLIT bp   (short -> $SHORT_ALN_DESC, long -> $ALIGNER)"
fi
echo "  Threads     : $THREADS"
echo ""

# --- Reference + indexes ----------------------------------------------------
if [ ! -f "$REF" ]; then
    echo "You need to have $REF in the current directory."
    echo "Do you want to download hs37d5? (y/n)"
    read -p "" choice
    if [ "$choice" = "y" ] || [ "$choice" = "Y" ]; then
        curl --insecure -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
        bgzip -d -@ "$THREADS" hs37d5.fa.gz
    else
        exit 1
    fi
fi

if [ ! -f "${REF}.fai" ]; then
    echo "Indexing reference (samtools faidx)"
    samtools faidx "$REF" || exit 1
fi

# Short reads always use bwa aln -> we always need the classic bwa index.
if [ ! -f "${REF}.bwt" ]; then
    echo "No bwa index (${REF}.bwt) found. Build it now? This can take hours. (y/n)"
    read -p "" choice
    if [ "$choice" = "y" ] || [ "$choice" = "Y" ]; then
        bwa index "$REF" || exit 1
    else
        exit 1
    fi
fi

# bwa-mem3 needs its own index.
if [ "$ALIGNER" = "mem3" ] && [ ! -f "${REF}.bwt.2bit.64" ]; then
    echo "No bwa-mem3 index (${REF}.bwt.2bit.64) found. Build it now? (y/n)"
    read -p "" choice
    if [ "$choice" = "y" ] || [ "$choice" = "Y" ]; then
        "$BWA_MEM3" index --max-memory 80G "$REF" || exit 1
    else
        exit 1
    fi
fi

# --- Output directory -------------------------------------------------------
if [ -d "$INDNAME" ]; then
    echo "Destination directory '$INDNAME' exists. Overwrite (delete it)? (y/n)"
    read -p "" choice
    if [ "$choice" = "y" ] || [ "$choice" = "Y" ]; then
        rm -rf "$INDNAME"; mkdir "$INDNAME"
    else
        echo "Quitting"; exit 1
    fi
else
    mkdir "$INDNAME"
fi

WORKDIR="$(pwd)"

# Resolve input paths (accept relative or absolute).
if [ -f "$WORKDIR/$FASTQ1" ]; then FASTQONE="$WORKDIR/$FASTQ1"; else FASTQONE="$FASTQ1"; fi
if [ "$MODE" = "paired" ]; then
    if [ -f "$WORKDIR/$FASTQ2" ]; then FASTQTWO="$WORKDIR/$FASTQ2"; else FASTQTWO="$FASTQ2"; fi
fi

RG="@RG\tID:ILLUMINA-${INDNAME}\tSM:${INDNAME}\tPL:illumina\tPU:ILLUMINA-${INDNAME}-${MODE^^}"
# ${MODE^^} -> SINGLE/PAIRED; trim to SE/PE for the PU tag.
[ "$MODE" = "single" ] && RG="${RG/-SINGLE/-SE}" || RG="${RG/-PAIRED/-PE}"

# ===========================================================================
# 1. Adapter removal
# ===========================================================================
echo ""; echo ">>> [1/5] Adapter removal (AdapterRemoval3)"; echo ""

PROC_FQ="${INDNAME}/${INDNAME}.proc.fq.gz"

if [ "$MODE" = "paired" ]; then
    if [[ "$FASTQONE" == *.bz2 ]] && [[ "$FASTQTWO" == *.bz2 ]]; then
        IN1=(--in-file1 <(lbzcat "$FASTQONE")); IN2=(--in-file2 <(lbzcat "$FASTQTWO"))
    else
        IN1=(--in-file1 "$FASTQONE"); IN2=(--in-file2 "$FASTQTWO")
    fi
    # Merge overlapping read pairs; all outputs collapse into one fastq, as in
    # paired_mem3.sh -- collapsed inserts + singletons become single-end reads.
    adapterremoval3 \
        "${IN1[@]}" "${IN2[@]}" \
        --out-file1   "$PROC_FQ" \
        --out-file2   "$PROC_FQ" \
        --out-merged  "$PROC_FQ" \
        --out-singleton "$PROC_FQ" \
        --out-html "${INDNAME}/${INDNAME}.html" \
        --out-json "${INDNAME}/${INDNAME}.json" \
        --threads "$THREADS" \
        --merge-quality-max 41 \
        --merge \
        --preserve5p \
        --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA \
        --min-length 30 || { echo "ERROR: AdapterRemoval failed" >&2; exit 1; }
else
    if [[ "$FASTQONE" == *.bz2 ]]; then
        IN1=(--in-file1 <(lbzcat "$FASTQONE"))
    else
        IN1=(--in-file1 "$FASTQONE")
    fi
    adapterremoval3 \
        "${IN1[@]}" \
        --out-file1 "$PROC_FQ" \
        --out-html "${INDNAME}/${INDNAME}.html" \
        --out-json "${INDNAME}/${INDNAME}.json" \
        --threads "$THREADS" \
        --min-length 30 \
        --trim-qualities \
        --preserve5p \
        --adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
        --quality-format sam || { echo "ERROR: AdapterRemoval failed" >&2; exit 1; }
fi

[ -s "$PROC_FQ" ] || { echo "ERROR: processed fastq $PROC_FQ is empty" >&2; exit 1; }

# Resolve an automatic split length from the AdapterRemoval length report.
# (-n 0.01 below must match the fnr used by the short-read aln step.)
if [ "$SKIP_ALN" -ne 1 ] && [ "$SPLIT" = "auto" ]; then
    AR_JSON="${INDNAME}/${INDNAME}.json"
    SPLIT="$(python3 "$SCRIPT_DIR/pick_split.py" "$AR_JSON" 0.01 ${PICK_SPLIT_BUDGET:+$PICK_SPLIT_BUDGET})"
    echo ">>> auto split length: ${SPLIT} bp  (from $AR_JSON; reads <${SPLIT} -> $SHORT_ALN_DESC, >=${SPLIT} -> $ALIGNER)"
fi

# ===========================================================================
# 2. Split short vs long reads at $SPLIT bp
# ===========================================================================
SHORT_FQ="${INDNAME}/${INDNAME}.short.fq.gz"
LONG_FQ="${INDNAME}/${INDNAME}.long.fq.gz"
has_reads() { [ "$(zcat -f "$1" 2>/dev/null | head -n 1 | wc -l)" -gt 0 ]; }
SHORT_HAS=0; LONG_HAS=0

if [ "$SKIP_ALN" -eq 1 ]; then
    # -A: no split -- every read goes to the long-read aligner (quick first pass, no aln).
    echo ""; echo ">>> [2/5] Skipping the aln split (-A) -- mapping ALL reads with: $ALIGNER"; echo ""
    LONG_FQ="$PROC_FQ"
    has_reads "$LONG_FQ" && LONG_HAS=1
    [ $LONG_HAS = 0 ] && { echo "ERROR: no reads to map" >&2; exit 1; }
else
    echo ""; echo ">>> [2/5] Splitting reads at ${SPLIT} bp"; echo ""

    # Compress the (transient, deleted-after-alignment) split fastqs with parallel BGZF
    # instead of single-threaded gzip -- ~6x faster on a many-core box, and BGZF is valid
    # gzip so bwa / mem3 / fq2bam read it natively. Two bgzip pipes run at once, so give
    # each half the threads.
    BGZIP_T=$(( THREADS > 1 ? THREADS / 2 : 1 ))
    : | bgzip -@ "$BGZIP_T" > "$SHORT_FQ"
    : | bgzip -@ "$BGZIP_T" > "$LONG_FQ"

    zcat -f "$PROC_FQ" | awk -v minlen="$SPLIT" \
        -v short_out="bgzip -@ $BGZIP_T -l 2 > $SHORT_FQ" -v long_out="bgzip -@ $BGZIP_T -l 2 > $LONG_FQ" '
        {
            if      (NR%4==1) header=$0;
            else if (NR%4==2) seq=$0;
            else if (NR%4==3) plus=$0;
            else if (NR%4==0) {
                qual=$0;
                rec = header "\n" seq "\n" plus "\n" qual;
                if (length(seq) < minlen) print rec | short_out;
                else                      print rec | long_out;
            }
        }'

    has_reads "$SHORT_FQ" && SHORT_HAS=1
    has_reads "$LONG_FQ"  && LONG_HAS=1
    echo "  short reads present: $([ $SHORT_HAS = 1 ] && echo yes || echo no)"
    echo "  long  reads present: $([ $LONG_HAS  = 1 ] && echo yes || echo no)"
    [ $SHORT_HAS = 0 ] && [ $LONG_HAS = 0 ] && { echo "ERROR: no reads after split" >&2; exit 1; }
fi

BAMS=()

# ===========================================================================
# 3. SHORT reads -> bwa aln (BWA-backtrack)
# ===========================================================================
if [ $SHORT_HAS = 1 ]; then
    echo ""; echo ">>> [3/5] Mapping SHORT reads with $SHORT_ALN_DESC (-l 1024 -n 0.01 -o 2)"; echo ""
    SHORT_BAM="${INDNAME}/${INDNAME}.short.bam"
    if [ "$GPU_ALN" -eq 1 ]; then
        # GPU BWA-backtrack: fused aln+samse (-S) in one pass (no .sai, single FASTQ read,
        # multithreaded samse). Output is bit-exact to the CPU bwa aln+samse path.
        "$BWA_GPU" gpualn -S -l 1024 -n 0.01 -o 2 -t "$THREADS" -r "$RG" "$REF" "$SHORT_FQ" \
            | samtools sort --no-PG -@ "$THREADS" -O bam -o "$SHORT_BAM" - \
            || { echo "ERROR: bwa gpualn failed" >&2; exit 1; }
    else
        bwa aln -l 1024 -n 0.01 -o 2 -t "$THREADS" "$REF" "$SHORT_FQ" > "${INDNAME}/${INDNAME}.short.sai" \
            || { echo "ERROR: bwa aln failed" >&2; exit 1; }
        bwa samse -r "$RG" "$REF" "${INDNAME}/${INDNAME}.short.sai" "$SHORT_FQ" \
            | samtools sort --no-PG -@ "$THREADS" -O bam -o "$SHORT_BAM" - \
            || { echo "ERROR: bwa samse/sort failed" >&2; exit 1; }
    fi
    BAMS+=("$SHORT_BAM")
fi

# ===========================================================================
# 4. LONG reads -> selected aligner
# ===========================================================================
if [ $LONG_HAS = 1 ]; then
    echo ""; echo ">>> [4/5] Mapping LONG reads with: $ALIGNER"; echo ""
    LONG_BAM="${INDNAME}/${INDNAME}.long.bam"

    case "$ALIGNER" in
        mem)
            # Classic bwa mem. No -p: merged reads are single-end records.
            bwa mem -t "$THREADS" $MEM_OPTS -R "$RG" "$REF" "$LONG_FQ" \
                | samtools sort --no-PG -@ "$THREADS" -m2G -O bam -o "$LONG_BAM" - \
                || { echo "ERROR: bwa mem failed" >&2; exit 1; }
            ;;
        mem3)
            # bwa-mem3 with the local hard-cap fix, as in paired_mem3.sh.
            "$BWA_MEM3" mem -t "$THREADS" $MEM_OPTS --supp-rep-hard-cap 5 -R "$RG" \
                "$REF" "$LONG_FQ" \
                | samtools sort --no-PG -@ "$THREADS" -m2G -O bam -o "$LONG_BAM" - \
                || { echo "ERROR: bwa-mem3 mem failed" >&2; exit 1; }
            ;;
        fq2bam)
            # NVIDIA Parabricks GPU fq2bam.
            #
            # NOTE: fq2bam's --bwa-options only accepts a subset of bwa mem
            # flags: -M -Y -C -T -B -U -L -I -K. The ancient seeding params
            # '-k 19 -r 2.5' are NOT supported and are therefore dropped here;
            # only -L (clip penalty) carries over. This is acceptable because
            # fq2bam handles the LONG-read fraction, which does not need the
            # aggressive short-read seeding.
            #
            # --no-markdups (NOT --align-only): skips fq2bam's duplicate marking
            # so dedup is deferred to the single sambamba pass in step 5, BUT
            # still coordinate-sorts the output (proper @HD SO:coordinate) and
            # leaves no orphan "PG:Z:MarkDuplicates" read tags -- so the merge
            # is clean with no header/sort fix-ups needed. (--align-only would
            # skip sorting and stamp every read with an orphan PG tag.)
            FQ2BAM_BWA_OPTS="-L 15"
            echo "  (fq2bam: passing --bwa-options=\"$FQ2BAM_BWA_OPTS\"; -k/-r are unsupported and omitted)"
            docker run --rm --runtime=nvidia \
                --volume "${WORKDIR}:/workdir" \
                --volume "${WORKDIR}/${INDNAME}:/outputdir" \
                --workdir /workdir \
                "$PARABRICKS_IMAGE" \
                pbrun fq2bam --no-markdups \
                --ref "/workdir/${REF}" \
                --in-se-fq "/outputdir/$(basename "$LONG_FQ")" \
                --out-bam "/outputdir/$(basename "$LONG_BAM")" \
                --read-group-sm "${INDNAME}" \
                --read-group-lb "ILLUMINA-${INDNAME}" \
                --read-group-pl "illumina" \
                --bwa-options="$FQ2BAM_BWA_OPTS" \
                --low-memory || { echo "ERROR: parabricks fq2bam failed" >&2; exit 1; }
            ;;
    esac
    BAMS+=("$LONG_BAM")
fi

# ===========================================================================
# 5. Merge + index (+ optional dedup)
# ===========================================================================
echo ""; echo ">>> [5/5] Merging BAMs and indexing"; echo ""

MAPPED_BAM="${INDNAME}/${INDNAME}_mapped.bam"

if [ "${#BAMS[@]}" -eq 1 ]; then
    # Only one fraction had reads -- pass it through sort to normalise the
    # header/ordering (cheap; input is already coordinate-sorted).
    samtools sort --no-PG -@ "$THREADS" -m2G -O bam -o "$MAPPED_BAM" "${BAMS[0]}" \
        || { echo "ERROR: samtools sort failed" >&2; exit 1; }
else
    # Merge the short (bwa aln) and long BAMs, then re-sort the stream. All
    # inputs are individually coordinate-sorted (bwa aln/mem via samtools sort,
    # fq2bam via --no-markdups), so a plain `samtools merge` would also work;
    # re-sorting is a cheap, defensive normalisation that guarantees a clean
    # coordinate order and a single consistent header regardless of source.
    samtools merge -f -u -@ "$THREADS" - "${BAMS[@]}" \
        | samtools sort --no-PG -@ "$THREADS" -m2G -O bam -o "$MAPPED_BAM" - \
        || { echo "ERROR: samtools merge/sort failed" >&2; exit 1; }
fi
samtools index -@ "$THREADS" "$MAPPED_BAM"

FINAL_BAM="$MAPPED_BAM"
if [ "$DEDUP" -eq 1 ]; then
    echo ""; echo ">>> Duplicate marking (sambamba markdup)"; echo ""
    DEDUP_BAM="${INDNAME}/${INDNAME}_markdup.bam"
    sambamba markdup -p -t "$THREADS" "$MAPPED_BAM" "$DEDUP_BAM" \
        || { echo "ERROR: sambamba markdup failed" >&2; exit 1; }
    samtools index -@ "$THREADS" "$DEDUP_BAM"
    FINAL_BAM="$DEDUP_BAM"
fi

# --- Cleanup intermediates --------------------------------------------------
rm -f "${INDNAME}/${INDNAME}.short.sai" "$SHORT_FQ" "$LONG_FQ" "$PROC_FQ"

# --- Log --------------------------------------------------------------------
{
    echo "adna_aligner.sh run on $(date)"
    echo "Mode: $MODE-end"
    echo "Reference: $REF"
    if [ "$SKIP_ALN" -eq 1 ]; then
        echo "Split length: none (aln skipped with -A; ALL reads->$ALIGNER)"
    else
        echo "Split length: $SPLIT bp (short->$SHORT_ALN_DESC, long->$ALIGNER)"
    fi
    echo "Long-read options: $MEM_OPTS"
    echo "Dedup: $([ $DEDUP -eq 1 ] && echo yes || echo no)"
    echo "Command: $0 $*"
    echo "System: $(uname -a)"
} > "${INDNAME}/pipeline.log"

echo ""
echo "Done. Final BAM: $FINAL_BAM"
echo ""
exit 0
