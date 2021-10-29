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

my $d1 = &Fxtran::fxtran (location => $F90_1, fopts => [qw (-line-length 300)]);
my $d2 = &Fxtran::fxtran (location => $F90_2, fopts => [qw (-line-length 300)]);


my %k = qw (UVH 0 XYB 1 RCP 1 CTY 0);

my @dummy = &F ('//dummy-arg-LT/arg-N', $d2, 1);
my ($name) = &F ('//subroutine-N', $d2, 1);

for my $call (&F ('.//call-stmt[string(procedure-designator)="?"]', $name, $d1))
  {
    my @actual = &F ('./arg-spec/arg/*', $call);
    my %d2a;
    for my $i (0 .. $#actual)
      {
        my $actual = $actual[$i]->textContent;
        next if ($actual =~ m/\(/o);
#       next unless ($actual =~ m/^(?:YDCPG_DYN0|YDCPG_DYN9|YDCPG_PHY0|YDCPG_PHY9|YDCPG_TND|YDCPG_MISC)%/o);
        next unless ($actual =~ m/^(?:YLAPLPAR%)/o);
        $d2a{$dummy[$i]} = $actual;
      }

    while (my ($k, $v) = each (%d2a))
      {
        my @expr = &F ('.//named-E[string(N)="?"]', $k, $d2);
        for my $expr (@expr)
          {
            my ($N) = &F ('./N', $expr); my $NN = $N->textContent;
            my ($r) = &F ('./R-LT/ANY-R', $expr);
            if ($r && grep { ("P${_}0" eq $NN) || ("P${_}9" eq $NN) } keys (%k))
              {
                my ($n) = ($NN =~ m/^P(\w+)\d$/o);
                if ($r->nodeName eq 'parens-R')
                  {
                    my @e = &F ('./element-LT/element/ANY-E', $r);
                    if (($e[0]->textContent eq '1') && ($e[1]->textContent eq '1'))
                      {
                        $e[0]->replaceNode (&t (':'));
                        if ($k{$n} eq 1)
                          {
                            $e[1]->replaceNode (&t (':'));
                          }
                        else
                          {
                            $e[1]->replaceNode (&t ('1:'));
                          }
                      }
                    elsif (($e[0]->textContent eq '1') && ($e[1]->textContent eq '0'))
                      {
                        $e[0]->replaceNode (&t (':'));
                        if ($k{$n} eq 0)
                          {
                            $e[1]->replaceNode (&t (':'));
                          }
                        else
                          {
                            die $expr->textContent;
                          }
                      }
                  }
                else
                  {
                    die $r->textContent;
                  }
              }
            $N->replaceNode (&t ($v));
          }
      }

  }




'FileHandle'->new (">$F90_2.new")->print ($d2->textContent);

