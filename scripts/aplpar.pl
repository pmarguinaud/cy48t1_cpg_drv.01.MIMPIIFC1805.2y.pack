#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my $F90 = "src/local/arpifs/phys_dmn/mf_phys.F90";

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);


for my $r (qw (CPTENDSM))
  {
    my @call = &F ('//call-stmt[string (procedure-designator)="?"]/arg-spec', $r, $d);
    
    for my $i (0 .. $#call)
      {
        my @args = map { $_->textContent } &F ('./arg', $call[$i]);
        'FileHandle'->new (">$r.$i.txt")->print (join ("\n", @args));
      }
  }
    
    
    
    
    
    
    
    
