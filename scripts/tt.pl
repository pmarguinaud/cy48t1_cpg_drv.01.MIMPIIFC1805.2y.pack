#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my $F90 = shift;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);

#my @expr = &F ('//named-E[string(.)="SUM(ZTESTSAVE(KIDIA:KFDIA))"]', $d);
 my @expr = &F ('//named-E[translate(string(.)," ","")="SUM(ZTESTSAVE(KIDIA:KFDIA))"]', $d);

print &Dumper (\@expr);

