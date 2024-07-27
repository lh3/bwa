#!/usr/bin/env perl

use strict;
use warnings;

die("Usage: qualfa2fq.pl <in.fasta> <in.qual>\n") if (@ARGV != 2);

my ($fhs, $fhq, $q);
open($fhs, ($ARGV[0] =~ /\.gz$/)? "gzip -dc $ARGV[0] |" : $ARGV[0]) || die;
open($fhq, ($ARGV[1] =~ /\.gz$/)? "gzip -dc $ARGV[1] |" : $ARGV[1]) || die;

$/ = ">"; <$fhs>; <$fhq>; $/ = "\n";
while (<$fhs>) {
  $q = <$fhq>;
  print "\@$_";
  $/ = ">";
  $_ = <$fhs>; $q = <$fhq>;
  chomp; chomp($q);
  $q =~ s/\s*(\d+)\s*/chr($1+33)/eg;
  print $_, "+\n";
  for (my $i = 0; $i < length($q); $i += 60) {
	print substr($q, $i, 60), "\n";
  }
  $/ = "\n";
}

close($fhs); close($fhq);
