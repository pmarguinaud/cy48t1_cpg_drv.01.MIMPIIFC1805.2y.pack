#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use List::MoreUtils qw (uniq);
use lib $Bin;
use Fxtran;

my $d = &Fxtran::fxtran (location => $ARGV[0], fopts => [qw (-line-length 400)]);

my $s = 'YDMF_PHYS_TMP';

my @expr = &F ('.//named-E[string(N)="?"]', $s, $d);

my %ct;

for my $expr (@expr)
  {
    my @ct = &F ('./R-LT/component-R/ct', $expr, 1);
    $ct{join ('%', @ct)}++;
  }

my @ct = sort keys (%ct);

print &Dumper (\%ct);
