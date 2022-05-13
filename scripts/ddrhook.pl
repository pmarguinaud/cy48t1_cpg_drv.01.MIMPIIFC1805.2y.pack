#!/usr/bin/perl -w

use strict;
use FileHandle;
use Data::Dumper;

my $f2f = do ('./f2f.pl');


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


my $t0 = &parse ($ARGV[0]);
my $t1 = &parse ($ARGV[1]);

#print &Dumper ($t0);
#print &Dumper ($t1);

my @k = keys (%{ {%$t0, %$t1} });

@k = sort { (exists ($t1->{$b}) ? $t1->{$b} : -9e99) <=> (exists ($t1->{$a}) ? $t1->{$a} : -9e99) } @k;

my ($t, $tm, $tp) = (0, 0, 0);

for my $k (@k)
  {
    my ($T0, $T1) = ($t0->{$k} || 0, $t1->{$k} || 0);

    my $r = $T0 ? $T1 / $T0 : ($T1 ? 9e99 : 1.);

    printf("%12.2f | %-60s |  %12.5f -> %12.5f | %10s | %+12.5f\n", 
           $t, $k, $T0, $T1, ($T0 ? sprintf ("%10.2f", $T1 / $T0) : ($T1 ? 'Infinity' : sprintf ("%10.2f", 1.0))), $T1 - $T0) if (abs ($r - 1) > 0.1);
    my $dt = $T1 - $T0;
    $t += $dt;
    $tp += $dt if ($dt > 0);
    $tm += $dt if ($dt < 0);
  }

print "\n" x 2;
printf (" %+12.5f %+12.5f %+12.5f\n", $t, $tp, $tm);
print "\n" x 2;


