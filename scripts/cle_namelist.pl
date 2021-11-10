#!/usr/bin/perl -w
#
use strict;
use FileHandle;
use File::Find;
use Data::Dumper;
use File::Basename;

my @nam;
my %k2n;

sub wanted
{
  return unless (m/\.nam\.h$/o);
  return if (m/\.include/o);
  my $nam = $File::Find::name;
  push @nam, $nam;
  my @text = do { my $fh = 'FileHandle'->new ("<$nam"); <$fh> };
  for (@text)
    {
      s/!.*//o;
    }
  @text = grep { !/::/o } @text;

  (my $namelist = &basename ($nam)) =~ s/\.nam.h$//o;
  $namelist = uc ($namelist);

  my @w = (join ('', @text) =~ m/(\w+)/goms);

  for (@w)
    {
      $k2n{$_} = $namelist;
    }
}

&find ({wanted => \&wanted, no_chdir => 1}, 'src/main/');

my %c;

my $cle = shift;
my $fh = 'FileHandle'->new ("<$cle");

while (my $line = <$fh>)
  { 
    next if ($line =~ m/^\s*$/o);
    next if ($line =~ m/^\s*#/o);
    my ($k, $v) = split (m/\s*=\s*/o, $line);
    for ($k, $v)
      {
        s/^\s*//o;
        s/\s*$//o;
      }
    if (($v =~ m/^(?:TRUE|FALSE)$/o) || ($v =~ m/^\'/o))
      {
        $c{$k} = ".$v.";
        next;
      }
    print STDERR "$k = $v\n";
  }

$fh->close ();

my %n2k;

for my $c (keys (%c))
  {
    unless ($k2n{$c})
      {
        print STDERR "$c not found\n";
        next;
      }
    $n2k{$k2n{$c}}{$c} = 1;
  }


for my $n (sort keys (%n2k))
  {
    print "&$n\n";
    for my $k (sort keys (%{ $n2k{$n} }))
      {
        print "  $k=$c{$k},\n";
      }
    print "/\n";
  }






