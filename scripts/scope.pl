#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use List::Util;
use List::MoreUtils qw (uniq);
use lib $Bin;
use Fxtran;

my $F90 = shift;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 400)], xopts => [qw (line_numbers 1)]);


my @skip = qw (SC2PRG DR_HOOK YLMF_PHYS_BASE_STATE%INIT YLMF_PHYS_NEXT_STATE%INIT);
my %skip = map { ($_, 1) } @skip;

my %size =
(
  NGFL_EXT => 0,
  NUMFLDS  => 18,
  KSW      => 0,
  NLEVS    => 1,
  KVCLIS     => 0,
  KDTPREC    => 0,
  KDTPREC2   => 0,
  KTSSG      => 0,
  KGRADIENTS => 0,
  KMAXDRAFT  => 0,
  KFLEVG     => 105,
  '0:KFLEVG' => 106,
  'KTSSG+1'  => 1,
  N2D        => 1,
  '0:1'      => 2,
);

my @vars;

for my $en_decl (&F ('.//ANY-stmt[not(./attribute[string(attribute-N)="INTENT"])]'
                   . '//EN-decl'
                   . '[./array-spec/shape-spec-LT'
                   . '[string(shape-spec)="YDCPG_DIM%KLON"]'
                   . ']'
                   , $d))
  {
    my ($N) = &F ('./EN-N', $en_decl, 1);

    push @vars, $N;

    my @ss = &F ('./array-spec/shape-spec-LT/shape-spec', $en_decl);
    shift (@ss);
  
    @ss = map { $_->textContent } @ss;    
    for (@ss)
      {
        1 while (s/\w+%//o);
        s/\b1://go;
      }

   
    my $size = 1;
    for my $ss (@ss)
      {
        die $ss unless (exists $size{$ss} || ($ss =~ m/^\d+$/o));
        $size *= exists $size{$ss} ? $size{$ss} : $ss;
      }

    $size{$N} = $size;
  }


my %vscope;

for my $var (@vars)
  {
    my @expr = &F ('.//named-E[string(N)="?"]', $var, $d);
    my ($expr1, $expr2) = ($expr[0], $expr[-1]);
    my $line1 = $expr1->line_number ();
    my $line2 = $expr2->line_number ();
    $vscope{$var} = [$line1, $line2];
  }


my %cscope;

my @call;

for my $call (&F ('.//call-stmt', $d))
  {
    my ($N) = &F ('./procedure-designator', $call, 1);
    next if ($skip{$N});
    my $line1 = $call->line_number ();
    my $lines = ($call->textContent =~ m/(\n)/goms);
    my $line2 = $line1 + $lines;
    $cscope{$N} = [$line1, $line2];
    push @call, $N;
  }


my @lines = sort { $a <=> $b } &uniq (map { @$_ } (values (%cscope), values (%vscope)));

my ($line1, $line2) = ($lines[0], $lines[-1]);

my $Size = 0;
for my $var (@vars)
  {
    $Size += $size{$var};
  }


my $fhd = 'FileHandle'->new (">aplpar.dat");
my $fhl = 'FileHandle'->new (">aplpar.lab");

my $sizeMax = 0;

my $Call = '';
my $count = 0;

for my $line ($line1 .. $line2)
  {
    my $size = 0;
    for my $var (@vars)
      {
        if (($vscope{$var}[0] <= $line) && ($line <= $vscope{$var}[1]))
          {
            $size += $size{$var};
          }
      }

    $sizeMax = $size if ($size > $sizeMax);

    $fhd->printf ("%8d %10d %10d\n", $line, $size, $Size);

    for my $call (@call)
      {
        if (($cscope{$call}[0] <= $line) && ($line <= $cscope{$call}[1]))
          {
            if ($call ne $Call)
              {
                (my $str = $call) =~ s/_/\\\\_/go;
                if ($count)
                  {
                    $fhl->printf ("set label \"%s\" at %d,%d rotate by -90\n", $str, $line, $Size);
                  }
                else
                  {
                    $fhl->printf ("set label \"%s\" at %d,%d rotate by +90\n", $str, $line, 1000);
                  }
                $Call = $call;
                $count = ($count + 1) % 2;
              }
          }
      }

  }

$fhd->close ();

