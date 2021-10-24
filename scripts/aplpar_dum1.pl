#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my $F90 = "src/local/arpifs/phys_dmn/aplpar_dum.F90";

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 200)]);
    
my ($sub) = &F ('.//subroutine-stmt[string(subroutine-N)="APLPAR_DUM"]', $d);

my @cnt = &F ('.//cnt', $sub);

for (@cnt)
  {
    $_->unbindNode ();
  }

for my $n (&F ('.//text()', $sub))
  {
    my $t = $n->getData;
    $t =~ s/(?:^\s*|\s*$)//go;
    $n->setData ($t);
  }

for my $arg (&F ('./dummy-arg-LT/arg-N', $sub))
  {
    $arg->parentNode->insertBefore (&t ("\n"), $arg);
  }

#print $sub->textContent, "\n";

'FileHandle'->new (">$F90")->print ($d->textContent);


