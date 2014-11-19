#!/bin/bash

ctg_opt=""
if [ $# -gt 1 ] && [ $1 == '-A' ]; then
	ctg_opt="-A"
	shift
fi

if [ $# -eq 0 ]; then
	echo "Usage: $0 <prefix>"
	exit 1
fi

for f in $1.HLA-*.fq; do
	gene=`echo $f | perl -pe 's/^.*(HLA-[A-Z]+[0-9]*).*fq$/$1/'`
	echo -e "\n*** Processing gene $gene...\n" >&2
	`dirname $0`/typeHLA.sh $ctg_opt $1 $gene
done

ls $1.HLA-*.gt | xargs -i echo grep ^GT {} \| head -1 | sh | sed "s,^GT,$1,"
