#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use List::MoreUtils qw (uniq all);
use lib $Bin;
use Fxtran;
use Associate;


my $F90 = shift;

my $doc = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);
&Associate::resolveAssociates ($doc);

my ($pu0) = &F ('.//program-unit[not(ancestor::program-unit)]', $doc);

my @N0 = (&F ('.//EN-N', $pu0, 1), &F ('./use-stmt//use-N', $pu0, 1));
my %N0 = map { ($_, 1) } @N0;

my @pu = &F ('.//program-unit[ancestor::program-unit]', $doc);

for my $pu (@pu)
  {
    my @N = (&F ('.//EN-N', $pu, 1), &F ('./use-stmt//use-N', $pu, 1));
    my %N = map { ($_, 1) } @N;

    my @n = &F ('.//named-E/N/n/text()', $pu, 1);

    my %seen;

    for my $n (@n)
      {
        next if ($N{$n});
        next if ($seen{$n}++);
        print $n, "\n";
      }

    print "\n";

  }



