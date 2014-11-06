#!/bin/bash

is_ctg=0

if [ $# -gt 1 ] && [ $1 == '-A' ]; then
	is_ctg=1
	shift
fi

if [ $# -lt 2 ]; then
	echo "Usage: $0 [-A] <prefix> <gene>"
	exit 1
fi

root=`dirname $0`
pre=$1.$2

# de novo assembly
if [ $is_ctg -eq 0 ]; then
	len=`$root/seqtk comp $pre.fq | awk '{++x;y+=$2}END{printf("%.0f\n", y/x)}'`
	$root/fermi2.pl unitig -t2 -l$len -p $pre.tmp $pre.fq > $pre.tmp.mak
	make -f $pre.tmp.mak EXE_FERMI2=$root/fermi2 EXE_ROPEBWT2=$root/ropebwt2
else
	ln -sf $pre.fq $pre.tmp.mag.gz
fi

# get contigs overlapping HLA exons
(ls $root/HLA-idx/*.fa | xargs -i $root/bwa mem -t2 -B1 -O1 -E1 {} $pre.tmp.mag.gz 2>/dev/null) | grep -v ^@ | sort -k3,3 -k4,4n | bgzip > $pre.tmp.sam.gz
$root/tabix -Bpsam $pre.tmp.sam.gz $root/HLA-approx.anno | cut -f1 | sort | uniq > $pre.tmp.kept
gzip -dc $pre.tmp.sam.gz | perl -ane 'print "$F[0]\n" if /AS:i:(\d+).*XS:i:(\d+)/&&$1==$2' | cut -f1 | sort | uniq > $pre.tmp.dropped
awk -v f=$pre.tmp.dropped 'BEGIN{while((getline<f)>0)l[$1]=1}!l[$1]' $pre.tmp.kept | $root/seqtk subseq $pre.tmp.mag.gz - | gzip -1 > $pre.tmp.fq.gz

# map HLA exons to de novo contigs
$root/bwa index -p $pre.tmp $pre.tmp.fq.gz 2> /dev/null
$root/bwa mem -aD.1 -t2 $pre.tmp $root/hla.fa | egrep "^(@|$2)" > $pre.sam

# type HLA
$root/k8 $root/bwa-typeHLA.js $pre.sam > $pre.gt

# delete temporary files
rm -f $pre.tmp.*
