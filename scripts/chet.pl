#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my ($F90_1, $F90_2) = ('src/local/arpifs/phys_dmn/writephysio.F90', 'src/local/arpifs/phys_dmn/profilechet.F90');

my $d1 = &Fxtran::fxtran (location => $F90_1, fopts => [qw (-line-length 800)]);
my $d2 = &Fxtran::fxtran (location => $F90_2, fopts => [qw (-line-length 800)]);

my @aa1 = &F ('//call-stmt[string(procedure-designator)="WRITEPROFILE"]/arg-spec/arg/*', $d1, 1);
my @aa2 = &F ('//call-stmt[string(procedure-designator)="WRITEPROFILE"]/arg-spec/arg/*', $d2, 1);

for my $i (0 .. $#aa1)
  {
    printf ("%8d | %-30s %-30s\n", $i, $aa1[$i], $aa2[$i]) if ($aa1[$i] ne $aa2[$i]);
  }

