#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my $F90 = "src/local/arpifs/phys_dmn/aplpar.F90";

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);
    
my ($call) = &F ('//call-stmt[string(procedure-designator)="RECMWF_DUM"]', $d);

my @cnt = &F ('.//cnt', $call);

for (@cnt)
  {
    $_->unbindNode ();
  }

for my $n (&F ('.//text()', $call))
  {
    my $t = $n->getData;
    $t =~ s/(?:^\s*|\s*$)//go;
    $n->setData ($t);
  }

for my $arg (&F ('./arg-spec/arg', $call))
  {
    $arg->parentNode->insertBefore (&t ("\n"), $arg);
  }

'FileHandle'->new (">$F90")->print ($d->textContent);


