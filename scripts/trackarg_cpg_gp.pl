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

my @arg = @ARGV;

@arg = qw (ZPRE9L ZPRE9M ZPRE9M ZGDW0 ZGDW9 ZGWHT0 ZGWHT9 ZOROGLL ZOROGMM ZOROGLM);


my %arg;

my @F90 = ("src/local/arpifs/adiab/cpg_gp.F90");

@F90 = reverse (@F90);

my @d = map { &Fxtran::fxtran (location => $_, fopts => [qw (-line-length 300)]) } @F90;

my $D = $d[0]->cloneNode (1);

for my $d (@d)
  {
    &Associate::resolveAssociates ($d);
  }

my @n = map { &F ('//program-unit[not (ancestor::program-unit)]/subroutine-stmt/subroutine-N', $_, 1) } @d;

my %d;

for my $i (0 .. $#d)
  {
    $d{$n[$i]} = $d[$i];
  }

MAIN : for my $arg (@arg)
  {

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
            next MAIN unless ($decl);
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
    
    
    my $f = shift (@trace);
    
    die unless (@$f == 1);
    
    $f = $f->[0];
    
    my ($n) = $f->[4];
    my $d = $d{$f->[1]};
    
    my ($e) = &F ('//pointer-a-stmt[E-1/named-E[string(N)="?"]]/E-2/named-E', $n, $d);
    next unless ($e);
    
    print "$arg ", $e->textContent, "\n";
    
    $d = $D;
    
    my @expr = &F ('//named-E[string(N)="?"]', $arg, $d);
    
    for my $expr (@expr)
      {
        my ($N) = &F ('./N', $expr);
        $N->replaceNode (&t ($e->textContent));
      }
    
    $arg{$arg}++;
  }

'FileHandle'->new (">$F90[0].new")->print ($D->textContent);
    
@arg = sort { length ($b) <=> length ($a) } grep { $arg{$_} } @arg;

print join ("\\|", map { "\\<$_\\>" } @arg), "\n";

@arg = map { s/^P/Z/o; $_ } @arg;
    
print join ("\\|", map { "\\<$_\\>" } @arg), "\n";
    
    
