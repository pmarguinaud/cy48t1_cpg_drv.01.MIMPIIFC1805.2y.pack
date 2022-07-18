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
    
my @use = &F ('.//use-stmt[string(module-N)="YOMDYNA"]', $d);

for my $use (@use)
  {
    my ($pu) = &F ('ancestor::program-unit', $use);
    my @N = &F ('.//use-N', $use, 1);

    for my $N (@N)
      {
        my @expr = &F ('.//named-E/N/n/text()[string(.)="?"]', $N, $pu);
        for my $expr (@expr)
          {
            $expr->setData ("YRDYNA%$N");
          }
      }

    my ($lt) = &F ('./rename-LT', $use);
    $lt->replaceNode (&t ("YRDYNA"));
  }


'FileHandle'->new (">$F90.new")->print ($d->textContent);


