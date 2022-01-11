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

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 500)]);

my @use = &F ('.//use-stmt', $d);

for my $use (@use)
  {
    my ($M) = &F ('./module-N', $use, 1);
    my @n = &F ('.//use-N/N/n', $use);
    for my $n (@n)
      {
        my @ren = &F ('.//use-stmt[string(module-N)="?"]//rename[.//use-N[string(N)="?"]]', 
                      $M, $n->textContent, $d);
        for my $ren (@ren)
          {
            my ($nr) = &F ('./use-N/N/n', $ren);
            next if ($nr->isSameNode ($n));
            $ren->nextSibling && $ren->nextSibling->unbindNode;
            $ren->unbindNode;
          }
      }
  }


'FileHandle'->new (">$F90.new")->print ($d->textContent);
