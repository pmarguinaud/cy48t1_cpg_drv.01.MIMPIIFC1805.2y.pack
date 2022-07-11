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

my $t0 = $d->textContent;

my @assoc = &F ('.//associate-construct', $d);

for my $assoc (@assoc)
  {
    my $stmt = $assoc->firstChild;
    &Fxtran::expand ($stmt);

    my @N = &F ('./associate-LT/associate/associate-N', $stmt);

    for my $N (@N)
      {
        my $n = $N->textContent;
        my @expr = &F ('.//named-E[string(N)="?"]', $n, $d);
        next if (@expr);
        my $associate = $N->parentNode;
        &Fxtran::removeListElement ($associate);
      }

    &Fxtran::fold ($stmt);

  }

my $t1 = $d->textContent;

'FileHandle'->new (">$F90.new")->print ($t1) if ($t0 ne $t1);
