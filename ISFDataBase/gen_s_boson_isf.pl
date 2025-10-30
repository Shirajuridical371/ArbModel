#!/usr/bin/perl -w

use strict;

my $n1=1001; # start with this seniority
my $n2=10000; # end with this seniority
my ($i,$filename);

for ($i=$n1;$i<=$n2;$i++) {
  $filename=sprintf("disf_un01_n%02d.dat",$i);
  print "$filename\n";
  system "echo 1 > $filename";
}
