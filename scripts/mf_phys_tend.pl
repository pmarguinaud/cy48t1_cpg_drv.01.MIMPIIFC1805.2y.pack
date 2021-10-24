#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my @d = qw (
             src/local/arpifs/phys_dmn
             .vimpack/src=/arpifs/adiab
             .vimpack/src=/arpifs/phys_dmn
           );

my $d = &Fxtran::fxtran (location => "src/local/arpifs/phys_dmn/mf_phys.F90", 
                         fopts => [qw (-line-length 300)]);



my @v = grep { m/TEND/o } &F ('//T-decl-stmt//EN-decl[array-spec/shape-spec-LT]/EN-N', $d, 1);

my %gfl;

for my $v (@v)
  {
    next if ($v eq 'ZTENDGFL');
    my @sc2prg = &F ('//call-stmt[string (procedure-designator)="SC2PRG"][./arg-spec/arg[string (.)="?"]]', $v, $d);
    if (@sc2prg)
      {
        $gfl{$v} = 1;
      }
  }

my @call = &F ('//call-stmt[.//named-E['
             . join (' or ', map { "string (N)=\"$_\"" } @v)
             . ']]', $d);

my @sub = map { &F ('./procedure-designator', $_, 1) } @call;

my %skip = map { ($_, 1) } qw (DR_HOOK SC2PRG);
my @ind = grep { ! $skip{$sub[$_]} } (0 .. $#sub);
@sub  = @sub [@ind];
@call = @call[@ind];


my %intent;

for my $sub (@sub)
  {
    my @f = map { "$_/" . lc ($sub) . '.F90' } @d;
    my ($f) = grep { -f } @f;

    my $d = &Fxtran::fxtran (location => $f, fopts => [qw (-line-length 300)]);
    my @arg = &F ('.//program-unit/subroutine-stmt/dummy-arg-LT/arg-N', $d, 1);
    my @int;
    for my $arg (@arg)
      {
        my ($intent) = &F ('.//program-unit/T-decl-stmt[.//EN-N[string (.)="?"]]/attribute/intent-spec', $arg, $d, 1);
        push @int, $intent;
      }
    $intent{$sub} = \@int;
  }



for my $i (0 .. $#sub)
  {
    printf (" %3d | %s\n", $i, $sub[$i]);
  }

my @cc;

for my $v (@v)
  {
    my $cc = '';
    for my $i (0 .. $#call)
      {
        my $call = $call[$i];
        my $sub  = $sub[$i];

        my @arg = &F ('./arg-spec/arg/ANY-E', $call);

        my @ind = grep { my $e = $arg[$_]; ($e->nodeName eq 'named-E') && ((&F ('./N', $e, 1))[0] eq $v) } (0 .. $#arg);

        my $int = (@ind == 0) ? ' ' : (@ind == 1) ? $intent{$sub}[$ind[0]] : '?';

        my %in = (IN => 'I ', OUT => ' O', INOUT => 'IO');
        $int = $in{$int} || '';

        $cc .= sprintf (" %2s |", $int);
      }
    $cc .= sprintf (" | %-20s", $v . ($gfl{$v} ? '  (ZTENDGFL)' : ''));
    $cc .= sprintf ("\n");
    push @cc, $cc;
  }

@cc = sort @cc;

my $count = '';

for my $i (0 .. $#sub)
  {
    $count .= sprintf ("%4d ", $i);
  }

$count .= "\n";

for my $i (0 .. $#cc)
  {
    if (($i % 20) == 0)
      {
        print $count;
      }
    print $cc[$i];
  }
