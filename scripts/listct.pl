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

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 400)]);

my @expr = &F ('.//named-E[./R-LT/component-R]', $d);

for my $expr (@expr)
  {
    my @r = &F ('./R-LT/ANY-R', $expr);
    if (@r && (($r[-1]->nodeName eq 'parens-R') || ($r[-1]->nodeName eq 'array-R')))
      {
        $r[-1]->unbindNode ();
      }
  }

for my $expr (@expr)
  {
    print $expr->textContent, "\n";
  }


