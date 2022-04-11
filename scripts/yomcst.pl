#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my ($F90) = @ARGV;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 800)]);

my @cst = &F ('.//use-stmt[string(module-N)="YOMCST"]//use-N', $d, 1);

for my $cst (@cst)
  {
    my @expr = &F ('.//named-E[string(N)="?"]/N/n/text()', $cst, $d);
    for my $expr (@expr)
      {
        $expr->setData ("YDCST%$cst");
      }
  }

'FileHandle'->new (">$F90.new")->print ($d->textContent);

