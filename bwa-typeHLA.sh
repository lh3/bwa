#!/bin/bash

if [ $# -lt 2 ]; then
	echo "Usage: $0 <prefix> <gene>"
	exit 1
fi

root=`dirname $0`
pre=$1.$2
len=`$root/seqtk comp $pre.fq | awk '{++x;y+=$2}END{printf("%.0f\n", y/x)}'`

$root/fermi2.pl unitig -t2 -l$len -p $pre.utg $pre.fq > $pre.utg.mak
make -f $pre.utg.mak
$root/bwa mem -B1 -O1 -E1 $root/HLA-ALT.fa $pre.utg.mag.gz | grep -v ^@ | sort -k3,3 -k4,4n | bgzip > $pre.utg.alt.gz
$root/tabix -Bpsam $pre.utg.alt.gz $root/HLA-approx.anno | cut -f1 | sort | uniq | $root/seqtk subseq $pre.utg.mag.gz - | gzip -1 > $pre.utg.fq.gz
rm $pre.utg.{flt,ec,raw,pre}.* $pre.utg.mag* $pre.utg.mak $pre.utg.alt.gz
$root/bwa index -p $pre.utg $pre.utg.fq.gz
$root/bwa mem -aD.1 -t2 $pre.utg $root/hla.fa | egrep "^(@|$2)" > $pre.utg.sam
rm $pre.utg.{ann,amb,bwt,sa,pac}
