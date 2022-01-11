#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my ($F90_1, $F90_2) = @ARGV;

my $d1 = &Fxtran::fxtran (location => $F90_1, fopts => [qw (-line-length 800)]);
my $d2 = &Fxtran::fxtran (location => $F90_2, fopts => [qw (-line-length 800)]);


my @dummy = &F ('//dummy-arg-LT/arg-N', $d2);

my ($name) = &F ('//subroutine-N', $d2, 1);

for my $call (&F ('.//call-stmt[string(procedure-designator)="?"]', $name, $d1))
  {
    my @actual = &F ('./arg-spec/arg/*', $call);


    my %d2a;
    for my $i (0 .. $#actual)
      {
        my $actual = $actual[$i];
        $d2a{$dummy[$i]->textContent} = $actual->textContent;
        die $actual->textContent if (&F ('//EN-N[string(N)="?"]', $actual->textContent, $d2) && ($dummy[$i]->textContent ne $actual->textContent ));
      }

    while (my ($k, $v) = each (%d2a))
      {
        my @expr = &F ('.//named-E[string(N)="?"]', $k, $d2);
        for my $expr (@expr)
          {
            my ($N) = &F ('./N', $expr); my $NN = $N->textContent;
            $N->replaceNode (&t ($v));
          }
        my @en_decl = &F ('.//EN-N[string(N)="?"]', $k, $d2);
        for my $en_decl (@en_decl)
          {
            $en_decl->replaceNode (&t ($v));
          }
        my @dum = &F ('//dummy-arg-LT/arg-N[string(.)="?"]', $k, $d2);
        for my $dum (@dum)
          {
            $dum->replaceNode (&t ($v));
          }
      }

    last;
  }

'FileHandle'->new (">$F90_2.new")->print ($d2->textContent);

