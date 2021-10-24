#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my $d = &Fxtran::fxtran (location => "src/local/arpifs/phys_dmn/mf_phys.F90", 
                         fopts => [qw (-line-length 300)]);

my @v = grep { m/TEND/o } &F ('//T-decl-stmt//EN-decl[array-spec/shape-spec-LT]/EN-N', $d, 1);

my @sub = &F ('//call-stmt[.//named-E['
            . join (' or ', map { "string (N)=\"$_\"" } @v)
            . ']]/procedure-designator', $d, 1);



my %skip = map { ($_, 1) } qw (DR_HOOK SC2PRG);
@sub = grep { ! $skip{$_} } @sub;

for my $i (0 .. $#sub)
  {
    printf (" %3d | %s\n", $i, $sub[$i]);
  }

my @cc;

for my $v (@v)
  {
    my @call = &F ('//call-stmt[.//named-E[string (N)="?"]]/procedure-designator', $v, $d, 1);
    my $cc = '';
    for my $sub (@sub)
      {
        $cc .= sprintf (" %s |", (grep ({ $_ eq $sub } @call) ? 'X' : ' '));
      }
    $cc .= sprintf (" | %-20s", $v);
    $cc .= sprintf ("\n");
    push @cc, $cc;
  }

@cc = sort @cc;

my $count = '';

for my $i (0 .. $#sub)
  {
    $count .= sprintf ("%3d ", $i);
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
