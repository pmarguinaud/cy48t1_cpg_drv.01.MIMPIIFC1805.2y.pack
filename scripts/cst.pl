#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FileHandle;
use FindBin qw ($Bin);
use lib $Bin;
use Fxtran;
use File::Basename;
use Data::Dumper;
use List::MoreUtils qw (uniq);

my @proc = qw (
GNHD3
GNHEE_REFINE_PREH
GNHEE_TNDLAGADIAB_GW
GNHEE_TNDLAGADIAB_SVD
GNHQE_TNDLAGADIAB_GW
GNHQE_TNDLAGADIAB_SVD
GNH_TNDLAGADIAB_GW
GNH_TNDLAGADIAB_SVD
GPCTY_FORC
);


my %proc = map { ($_, 1) } @proc;

for my $f (qw (src/local/arpifs/adiab/cpg_gp.F90 
               src/local/arpifs/adiab/cpg_gp_hyd.F90
               src/local/arpifs/adiab/cpg_gp_nhee.F90
               src/local/arpifs/adiab/cpg_gp_nhqe.F90))
  {
    my $d = &Fxtran::fxtran (location => $f, fopts => [qw (-line-length 300)]);
    my @call = &F ('.//call-stmt', $d);

    for my $call (@call)
      {
        my ($proc) = &F ('./procedure-designator', $call, 1);
        next unless ($proc{$proc});
        my ($argspec) = &F ('./arg-spec', $call);
        $argspec->insertBefore (&t (', '), $argspec->firstChild);
        $argspec->insertBefore (&n ('<arg><named-E><N><n>YDMODEL%YRCST</n></N></named-E></arg>'), $argspec->firstChild);
      }

    'FileHandle'->new (">$f.new")->print ($d->textContent);
  }


