#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

# EN-decl><EN-N><N><n>ZFGELS</n></N></EN-N><array-spec>(<shape-spec-LT><shape-spec><upper-bound
# T-decl-stmt><_T-spec_><derived-T-spec>TYPE(<T-N><N><n>GEOMETRY</n></N></T-N>)</derived-T-spec></_T-spec_>    ,<attribute><attribute-N>INTENT

my $d = &Fxtran::fxtran (location => "src/local/arpifs/adiab/cpg.F90", 
                         fopts => [qw (-line-length 300)]);

my @v = &F ('//T-decl-stmt[not (.//attribute-N[string (.)="INTENT"])]//EN-decl[array-spec/shape-spec-LT[starts-with(string(shape-spec),"YDGEOMETRY%YRDIM%NPROM")]]/EN-N', $d, 1);


my @sub = &F ('//call-stmt/procedure-designator', $d, 1);
my %skip = map { ($_, 1) } qw (DR_HOOK CP_PTRSLB1);
@sub = grep { ! $skip{$_} } @sub;

for my $i (0 .. $#sub)
  {
    printf (" %3d | %s\n", $i, $sub[$i]);
  }


my $count = '';

for my $i (0 .. $#sub)
  {
    $count .= sprintf ("%3d ", $i);
  }

$count .= "\n";


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

for my $i (0 .. $#cc)
  {
    if (($i % 20) == 0)
      {
        print $count;
      }
    print $cc[$i];
  }
