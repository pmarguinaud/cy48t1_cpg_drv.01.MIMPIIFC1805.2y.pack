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

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 800)]);

my @da = &F ('//dummy-arg-LT/arg-N', $d, 1);
my %da = map { ($_, 1) } @da;

my @en_decl = &F ('.//ANY-stmt//EN-decl', $d);

for my $en_decl (@en_decl)
  {
    my ($n) = &F ('./EN-N', $en_decl, 1);
    next if ($da{$n});
    my @expr = &F ('.//named-E[string(N)="?"]', $n, $d);
    unless (@expr)
      {
        $en_decl->unbindNode ();
      }
  }



'FileHandle'->new (">$F90.new")->print ($d->textContent);
