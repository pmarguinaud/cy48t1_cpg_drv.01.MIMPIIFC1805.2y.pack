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

my @ac = &F ('.//associate-construct', $d);

for my $ac (@ac)
  {
    my $stmt = $ac->firstChild;
    my @as = &F ('./associate-LT/associate', $stmt);

    for my $as (@as)
      {
        my ($N) = &F ('./associate-N', $as, 1);
        my @f = &F ('.//associate-construct/associate-stmt/associate-LT/associate[string(associate-N)="?"]', $N, $ac);
        for my $f (@f)
          {
            $f->unbindNode ();
          }
      }

  }

'FileHandle'->new (">$F90.new")->print ($d->textContent);
