#!/bin/bash

if [ $# -lt 2 ]; then
	echo "Usage: $0 <prefix> <gene>"
	exit 1
fi

root=`dirname $0`
pre=$1.$2
len=`$root/seqtk comp $pre.fq | awk '{++x;y+=$2}END{printf("%.0f\n", y/x)}'`

# de novo assembly
$root/fermi2.pl unitig -t2 -l$len -p $pre.tmp $pre.fq > $pre.tmp.mak
make -f $pre.tmp.mak

# get contigs overlapping HLA exons
$root/bwa mem -B1 -O1 -E1 $root/HLA-ALT.fa $pre.tmp.mag.gz | grep -v ^@ | sort -k3,3 -k4,4n | bgzip > $pre.tmp.alt.gz
$root/tabix -Bpsam $pre.tmp.alt.gz $root/HLA-approx.anno | cut -f1 | sort | uniq | $root/seqtk subseq $pre.tmp.mag.gz - | gzip -1 > $pre.tmp.fq.gz

# map HLA exons to de novo contigs
$root/bwa index -p $pre.tmp $pre.tmp.fq.gz
$root/bwa mem -aD.1 -t2 $pre.tmp $root/hla.fa | egrep "^(@|$2)" > $pre.tmp.sam

# type
$root/k8 $root/bwa-typeHLA.js $pre.tmp.sam > $pre.gt

# delete temporary files
rm -f $pre.tmp.*
