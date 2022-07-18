#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;
use Associate;

my $F90 = shift;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 500)]);

my ($geom) = &F ('.//EN-decl[string(EN-N)="YDGEOMETRY"]', $d);

exit (0) unless ($geom);


my @expr = &F ('.//named-E[string(N)="YDGEOMETRY"]', $d);

my $c = 0;

for my $expr (@expr)
  {
    my $t = $expr->textContent;
    next if ($t =~ m/NFLEVG|NPROMA/o);
    $c++;
  }

print "$F90\n" unless ($c);

