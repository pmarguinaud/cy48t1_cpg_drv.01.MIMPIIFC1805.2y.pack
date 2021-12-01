#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my $d = &Fxtran::fxtran (location => 'src/local/arpifs/phys_dmn/aplpar_init.F90', fopts => [qw (-line-length 800)]);

my @da = &F ('//dummy-arg-LT/arg-N', $d, 1);

my %del;

for my $da (@da)
  {
    my @expr = &F ('.//named-E[string(N)="?"]', $da, $d);
    unless (@expr)
      {
        $del{$da} = 1;
      }
  }

for (qw (KIDIA KFDIA))
  {
    delete $del{$_};
  }

for my $n (keys (%del))
  {
    my ($decl) = &F ('.//ANY-stmt[.//EN-decl[string(EN-N)="?"]]', $n, $d);
    $decl->unbindNode ();
    my ($da) = &F ('//dummy-arg-LT/arg-N[string(.)="?"]', $n, $d);
    if ($da->nextSibling)
      {
        $da->nextSibling->unbindNode ();
      }
    $da->unbindNode ();
  }


'FileHandle'->new (">src/local/arpifs/phys_dmn/aplpar_init.F90.new")->print ($d->textContent);
