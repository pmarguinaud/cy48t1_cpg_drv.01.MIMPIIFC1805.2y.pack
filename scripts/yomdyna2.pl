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
    

my ($ydmodel) = &F ('.//arg-N[string(.)="YDMODEL"]', $d);
my ($ydml_dyn) = &F ('.//arg-N[string(.)="YDML_DYN"]', $d);

exit (0) unless ($ydmodel || $ydml_dyn);

my ($use) = &F ('.//use-stmt[.//use-N[string(.)="YRDYNA"]]', $d);

exit (0) unless ($use);

my @expr = &F ('.//named-E[string(N)="YRDYNA"]', $d);

for my $expr (@expr)
  {
    my ($N) = &F ('./N/n/text()', $expr);
    if ($ydmodel)
      {
        $N->setData ("YDMODEL%YRML_DYN%YRDYNA");
      }
    elsif ($ydml_dyn)
      {
        $N->setData ("YDML_DYN%YRDYNA");
      }
  }


$use->unbindNode;

'FileHandle'->new (">$F90.new")->print ($d->textContent);
