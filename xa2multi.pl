#!/usr/bin/perl -w

use strict;
use warnings;

while (<>) {
	if (/\tXA:Z:(\S+)/) {
		my $l = $1;
		print;
		my @t = split("\t");
		while ($l =~ /([^,;]+),([-+]\d+),([^,]+),(\d+);/g) {
			my $mchr = ($t[6] eq $1)? '=' : $t[6]; # FIXME: TLEN/ISIZE is not calculated!
			print(join("\t", $t[0], 0x100|($t[1]&0x6e9)|($2<0?0x10:0), $1, abs($2), 0, $3, @t[6..7], 0, @t[9..10], "NM:i:$4"), "\n");
		}
	} else { print; }
}
