#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my ($F90_1, $F90_2) = @ARGV;

my $d1 = &Fxtran::fxtran (location => $F90_1, fopts => [qw (-line-length 300)]);
my $d2 = &Fxtran::fxtran (location => $F90_2, fopts => [qw (-line-length 300)]);


my %k = qw (UVH 0 XYB 1 RCP 1 CTY 0);

my @dummy = &F ('//dummy-arg-LT/arg-N', $d2, 1);
my ($name) = &F ('//subroutine-N', $d2, 1);

for my $call (&F ('.//call-stmt[string(procedure-designator)="?"]', $name, $d1))
  {
    my @actual = &F ('./arg-spec/arg/*', $call);

    print scalar (@dummy), " " , scalar (@actual), "\n";

    for my $i (0 .. $#actual)
      {
        my $actual = $actual[$i]->textContent;
        print $dummy[$i], " = ", $actual, "\n";
      }

    print "\n" x 3;


  }





