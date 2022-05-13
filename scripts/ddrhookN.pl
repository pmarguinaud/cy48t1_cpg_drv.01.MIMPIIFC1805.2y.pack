#!/usr/bin/perl -w

use strict;
use FileHandle;
use Data::Dumper;


sub slurp
{
  my $f = shift;
  my $fh = 'FileHandle'->new ("<$f");
  my @line = <$fh>;
  chomp for (@line);
  return @line;
}

sub parse
{
  my $f = shift;
  my @line = &slurp ($f);

  pop (@line);

  while (@line)
    {
      last if (shift (@line) =~ m/Avg-/o);
    }
  
  my %t;
  for (@line)
    {
      s/(?:^\s*|\s*$)//go;
      my @x = split (m/\s+/o);

      (my $k = $x[-1]) =~ s/_PARALLEL\b//o;
      $t{$k} = $x[1];
    }

  
  return \%t;
}

my @t = map { &parse ($_) } @ARGV;

my @k = keys (%{ {map {%$_} @t} });

@k = sort { (exists ($t[-1]->{$b}) ? $t[-1]->{$b} : -9e99) <=> (exists ($t[-1]->{$a}) ? $t[-1]->{$a} : -9e99) } @k;

my ($t, $tm, $tp) = (0, 0, 0);

for my $k (@k)
  {
    my ($T0, $T1) = ($t[0]->{$k} || 0, $t[-1]->{$k} || 0);
    my @T = map { $_->{$k} || 0.0 } @t;

    my $r = $T0 ? $T1 / $T0 : ($T1 ? 9e99 : 1.);

    if (abs ($r - 1) > 0.1)
      {
        printf("%12.2f %-60s ", $t, $k);
        for (@T)
          {
            printf (" %12.5f", $_);
          }
        print "\n";
      }

    my $dt = $T1 - $T0;
    $t += $dt;
    $tp += $dt if ($dt > 0);
    $tm += $dt if ($dt < 0);
  }

print "\n" x 2;
printf (" %+12.5f %+12.5f %+12.5f\n", $t, $tp, $tm);
print "\n" x 2;


