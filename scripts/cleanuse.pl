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

my @use = &F ('.//use-stmt', $d);

for my $use (@use)
  {
    my @n = &F ('.//use-N/N/n', $use);
    for my $n (@n)
      {
        my @expr = &F ('.//named-E[string(N)="?"]', $n->textContent, $d);
        my @type = &F ('.//T-N[string(.)="?"]', $n->textContent, $d);
        next if (@expr || @type);
        my ($rename) = &F ('ancestor::rename', $n);
        $rename->unbindNode ();
      }
  }

my $t1 = $d->textContent;

'FileHandle'->new (">$F90.new")->print ($t1) if ($t0 ne $t1);
