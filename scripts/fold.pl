#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my ($F90, $xpath) = @ARGV;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 500)]);


sub normalizeCommas
{
  my $x = shift;
  my @c = &F ('.//text()[normalize-space(.)=","]', $x);
  for my $c (@c)
    {
      my $t = $c->getData;
#     print ">$t<\n";
      $t =~ s/\s*,\s*/, /go;
#     print ">$t<\n";
      $c->setData ($t);
    }
}

for my $x (&F ($xpath, $d))
  {
    &Fxtran::expand ($x);
    &normalizeCommas ($x);
    print $x->textContent, "\n";
    &Fxtran::fold ($x);
  }

'FileHandle'->new (">$F90.new")->print ($d->textContent);



