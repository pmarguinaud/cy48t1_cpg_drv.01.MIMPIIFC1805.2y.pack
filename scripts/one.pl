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

my @proc = &F ('.//call-stmt/procedure-designator', $d, 1);

for my $proc (@proc)
  {
    my @arg = &F ('.//call-stmt[string(procedure-designator)="?"]/arg-spec/arg/named-E/N', $proc, $d, 1);

    print "==> $proc <==\n";

    for my $arg (@arg)
      {
        my ($decl) = &F ('.//T-decl-stmt[.//EN-decl[string(EN-N)="?"]/array-spec/shape-spec-LT[string(shape-spec)="YDCPG_DIM%KLON"]]'
                         , $arg, $d);
        next unless ($decl);

        my @expr = &F ('.//named-E[string(N)="?"]', $arg, $d);
        next if (@expr > 1);
        print $arg, "\n";
      }

  }

