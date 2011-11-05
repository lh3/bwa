#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;

my %opts = (a=>700, b=>100);
getopts('a:b:', \%opts);
die("Usage: pairsam.pl [-a $opts{a}] <read1.sam> <read2.sam>\n") if (@ARGV < 2);

my ($fh0, $fh1, $l0, $l1);

open($fh0, $ARGV[0] =~ /\.gz$/? "gzip -dc $ARGV[0] |" : $ARGV[0]) || die;
open($fh1, $ARGV[1] =~ /\.gz$/? "gzip -dc $ARGV[1] |" : $ARGV[1]) || die;

while ($l0 = <$fh0>) { last if $l0 !~ /^@/; print $l0 }
while ($l1 = <$fh1>) { last if $l1 !~ /^@/ }

while (defined($l0) && defined($l1)) {
	my ($r0, $r1) = &pair_line(\$l0, \$l1, $opts{a}, $opts{b});
	while ($l0 = <$fh0>) { last if $l0 !~ /^$r0/; print $l0; }
	while ($l1 = <$fh1>) { last if $l1 !~ /^$r1/; print $l1; }
}

close($fh0); close($fh1);

sub pair_line {
	my ($l0, $l1, $max_ins, $min_ins) = @_;
	my @t0 = split("\t", $$l0);
	my @t1 = split("\t", $$l1);
	my ($n0, $n1) = ($t0[0], $t1[0]);
	my ($cigar, $a0, $a1, $p0, $p1);
	# length in alignment
	$cigar = $t0[5]; $a0 = 0; $cigar =~ s/(\d+)[MI]/$a0 += $1/eg;
	$cigar = $t1[5]; $a1 = 0; $cigar =~ s/(\d+)[MI]/$a1 += $1/eg;
	# 5'-end alignment position on the read
	$p0 = $t0[1] == 16? $t0[3] + $a0 : $t0[3];
	$p1 = $t1[1] == 16? $t1[3] + $a1 : $t1[3];
	# adjust mapping quality
	if ($t0[2] eq $t1[2] && $t0[1]+$t1[1] == 16) { # on the same chr and forward-reverse
		if (abs($p0 - $p1) <= $max_ins && abs($p0 - $p1) >= $min_ins) { # within the right insert size distribution
			$t0[1] |= 2; $t1[1] |= 2; # flag as paired
			if ($t0[4] < $t1[4]) { # increase mapQ
				$t0[4] = $t0[4] + 10 < $t1[4]? $t0[4] + 10 : $t1[4];
			} else {
				$t1[4] = $t1[4] + 10 < $t0[4]? $t1[4] + 10 : $t0[4];
			}
		}
	}
	unless ($t0[1]&2) { # decrease mapQ if unpaired
		$t0[4] = $t0[4] > 10? $t0[4] - 10 : 0;
		$t1[4] = $t1[4] > 10? $t1[4] - 10 : 0;
	}
	# strip off /[12]
	$t0[0] =~ s/\/[12]$//; $t1[0] =~ s/\/[12]$//;
	# update FLAG
	$t0[1] |= 0x41 | (($t1[1]&16)? 0x20 : 0) | (($t1[1]&4)? 0x8 : 0);
	$t1[1] |= 0x81 | (($t0[1]&16)? 0x20 : 0) | (($t0[1]&4)? 0x8 : 0);
	# update mate positions
	if ($t0[2] eq $t1[2]) {
		$t0[6] = $t1[6] = '=';
		$t0[8] = $p1 - $p0; $t1[8] = $p0 - $p1;
	} else { $t0[6] = $t1[2]; $t1[6] = $t0[2]; }
	$t0[7] = $t1[3]; $t1[7] = $t0[3];
	# print out
	print join("\t", @t0);
	print join("\t", @t1);
	return ($n0, $n1);
}
