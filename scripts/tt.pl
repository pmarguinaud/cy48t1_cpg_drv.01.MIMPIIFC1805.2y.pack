#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use List::MoreUtils qw (uniq);
use lib $Bin;
use Fxtran;

for my $F90 (@ARGV)
  {

    my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 400)]);

    my @x = &F ('.//named-E[string(N)="PGMV" or string(N)="PGMVT1"]/R-LT/parens-R/element-LT/element[last()]/named-E', $d, 1);

    @x = sort &uniq (@x);

    print &Dumper (\@x);
    
#   'FileHandle'->new (">$F90.new")->print ($d->textContent);
    

  }
