#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my $d1 = &Fxtran::fxtran (location => 'src/local/arpifs/phys_dmn/mf_phys.F90',    fopts => [qw (-line-length 800)]);
my $d2 = &Fxtran::fxtran (location => 'src/local/arpifs/phys_dmn/initaplpar.F90', fopts => [qw (-line-length 800)]);

my @aa = &F ('//call-stmt[string(procedure-designator)="INITAPLPAR"]/arg-spec/arg/*', $d1, 1);

my @da = &F ('//dummy-arg-LT/arg-N', $d2, 1);

for my $i (0 .. $#da)
  {
    if ($aa[$i] =~ m/%TMP%APLPAR%/o)
      {
        my @stmt = &F ('.//a-stmt[E-1/named-E[string(N)="?"]]', $da[$i], $d2);
        for my $stmt (@stmt)
          {
            $stmt->parentNode->insertBefore (&t ("IF (LAPLPAR) "), $stmt);
          }
      }
  }

'FileHandle'->new (">src/local/arpifs/phys_dmn/initaplpar.F90.new")->print ($d2->textContent);
