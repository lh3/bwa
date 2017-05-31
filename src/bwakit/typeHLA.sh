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

preres="resource-human-HLA"
root=`dirname $0`
pre=$1.$2
touch $pre.gt

if [ ! -s $pre.fq ]; then
	echo '** Empty input file. Abort!' >&2
	exit 0
fi

if [ $is_ctg -eq 0 ]; then
	echo "** De novo assembling..." >&2
	len=`$root/seqtk comp $pre.fq | awk '{++x;y+=$2}END{printf("%.0f\n", y/x)}'`
	$root/fermi2.pl unitig -f $root/fermi2 -r $root/ropebwt2 -t2 -l$len -p $pre.tmp $pre.fq > $pre.tmp.mak
	make -f $pre.tmp.mak >&2
	cp $pre.tmp.mag.gz $pre.mag.gz
else
	rm -f $pre.tmp.mag.gz
	ln -s $pre.fq $pre.tmp.mag.gz
fi

echo "** Selecting contigs overlapping target exons..." >&2
(ls $root/$preres/HLA-ALT-idx/*.fa.bwt | sed s,.bwt,, | xargs -i $root/bwa mem -t2 -B1 -O1 -E1 {} $pre.tmp.mag.gz 2>/dev/null) | grep -v ^@ | sort -k3,3 -k4,4n | gzip > $pre.tmp.ALT.sam.gz
$root/k8 $root/typeHLA-selctg.js $2 $root/$preres/HLA-ALT-exons.bed $pre.tmp.ALT.sam.gz | $root/seqtk subseq $pre.tmp.mag.gz - | gzip -1 > $pre.tmp.fq.gz

echo "** Mapping exons to de novo contigs..." >&2
$root/bwa index -p $pre.tmp $pre.tmp.fq.gz 2>/dev/null
$root/seqtk comp $root/$preres/HLA-CDS.fa | cut -f1 | grep ^$2 | $root/seqtk subseq $root/$preres/HLA-CDS.fa - | $root/bwa mem -aD.1 -t2 $pre.tmp - 2>/dev/null | gzip -1 > $pre.sam.gz

echo "** Typing..." >&2
$root/k8 $root/typeHLA.js $pre.sam.gz > $pre.gt

# delete temporary files
rm -f $pre.tmp.*
[ $is_ctg -eq 1 ] && rm -f $pre.mag.gz
