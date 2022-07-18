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
    
my @expr = &F ('.//E-1/named-E[string(N)="YRDYNA"]', $d);

if (@expr)
{
  print $F90, "\n" ;
  for my $expr (@expr)
    {
      my $stmt = &Fxtran::stmt ($expr);
      print $stmt->textContent, "\n";
    }

}

