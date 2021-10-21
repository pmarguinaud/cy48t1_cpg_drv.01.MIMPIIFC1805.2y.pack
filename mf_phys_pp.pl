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

my @ap = &F ('//pointer-a-stmt', $d);

for my $ap (@ap)
  {
    my ($E1) = &F ('./E-1/ANY-E', $ap);
    my ($E2) = &F ('./E-2/ANY-E', $ap);
    next if ($E2->textContent =~ m/\(/o);


    my ($N1) = &F ('./N', $E1, 1);

    my @expr = &F ('//named-E/N[string(.)="?"]', $N1, $d);

    for my $expr (@expr)
      {
        my $stmt = &Fxtran::stmt ($expr);
        next if ($stmt->isSameNode ($ap));
        if ($stmt->nodeName eq 'call-stmt')
          {
            my ($p) = &F ('./procedure-designator', $stmt, 1);
#           next if ($p =~ m/^(?:APLPAR|APL_AROME)$/o);
          }
        $expr->replaceNode (&t ($E2->textContent));
      }

  }

'FileHandle'->new (">$F90.new")->print ($d->textContent);

