#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;
use Storable;
use Associate;

my $arg = shift;

my $sub = shift;

my @F90 = ("src/local/arpifs/control/gp_model.F90",
           "src/local/arpifs/adiab/cpg_drv.F90",
           "src/local/arpifs/adiab/cpg.F90",
           "src/local/arpifs/phys_dmn/mf_phys.F90",
           $sub ? ("src/local/arpifs/phys_dmn/$sub") : ());

@F90 = reverse (@F90);
          

my @d = map { &Fxtran::fxtran (location => $_, fopts => [qw (-line-length 300)]) } @F90;

for my $d (@d)
  {
    &Associate::resolveAssociates ($d);
  }

my @n = map { &F ('//program-unit[not (ancestor::program-unit)]/subroutine-stmt/subroutine-N', $_, 1) } @d;

my @trace;

for my $i (0 .. $#d)
  {
    my $F90 = $F90[$i];
    my $d = $d[$i];
    my $n = $n[$i];


    my @a = &F ('//dummy-arg-LT/arg-N', $d, 1);

    if ($i == 0)
      {
        my ($rank) = grep { $a[$_] eq $arg } (0 .. $#a);
        my ($decl) = &F ('//T-decl-stmt[.//EN-decl[string (EN-N)="?"]', $arg, $d);
        unshift @trace, [[$F90, $n, $rank, $decl->textContent]];
      }
    else
      {
        for my $t (@{ $trace[0] })
          {
            my ($nd, $rankd) = @{$t}[1,2];
            next unless (defined ($rankd));

            my @call = &F ('//call-stmt[string (procedure-designator)="?"]', $nd, $d);
            my @t;
            for my $call (@call)
              {
                my @arga = &F ('.//arg-spec/arg/node()', $call);
                my $arga = $arga[$rankd];

                die unless ($arga->nodeName eq 'named-E');

                my ($na) = &F ('./N', $arga, 1);

                my ($ranka) = grep { $a[$_] eq $na } (0 .. $#a);

                my ($decl) = &F ('//T-decl-stmt[.//EN-decl[string (EN-N)="?"]', $na, $d);

                push @t, [$F90, $n, $ranka, $decl && $decl->textContent, $arga->textContent];
       
              }
 
            my @tt = map { &Storable::nfreeze ($_) } @t;
            my %seen;
            @t = map { $t[$_] } grep { ! ($seen{$tt[$_]}++) } (0 .. $#t);

            unshift (@trace, \@t);
          }
      }

  }


print &Dumper (\@trace);

