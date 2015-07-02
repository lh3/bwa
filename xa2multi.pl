#!/usr/bin/env perl

use strict;
use warnings;

while (<>) {
	if (/\tXA:Z:(\S+)/) {
		my $l = $1;
		print;
		my @t = split("\t");
		while ($l =~ /([^,;]+),([-+]\d+),([^,]+),(\d+);/g) {
			my $mchr = ($t[6] eq $1)? '=' : $t[6]; # FIXME: TLEN/ISIZE is not calculated!
			my $seq = $t[9];
			my $phred = $t[10];
			# if alternative alignment has other orientation than primary, 
			# then print the reverse (complement) of sequence and phred string
			if ((($t[1]&0x10)>0) xor ($2<0)) {
				$seq = reverse $seq;
				$seq =~ tr/ACGTacgt/TGCAtgca/;
				$phred = reverse $phred;
			}
			print(join("\t", $t[0], 0x100|($t[1]&0x6e9)|($2<0?0x10:0), $1, abs($2), 0, $3, @t[6..7], 0, $seq, $phred, "NM:i:$4"), "\n");
		}
	} else { print; }
}
