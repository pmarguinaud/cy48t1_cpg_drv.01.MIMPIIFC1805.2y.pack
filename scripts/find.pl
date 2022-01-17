#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use File::Find;
use lib $Bin;

sub slurp
{
  my $file = shift;
  my $fh = 'FileHandle'->new ("<$file");
  local $/ = undef;
  my $data = <$fh>;
  return $data;
}

my %f2f;
my %seen;

sub wanted
{
  return unless (m/\.F90$/o);
  my $f = $File::Find::name;
  return if ($seen{&basename ($f)}++);
  my $code = &slurp ($f);
  my ($s) = ($code =~ m/.*(?:MODULE|FUNCTION|PROGRAM|SUBROUTINE)\s+(\w+)/goms);
  return unless ($s);
  $s = uc ($s);
  $f2f{$s} = $f;
}

my @view = do { my $fh = 'FileHandle'->new ("<.gmkview"); my @v = <$fh>; chomp for (@v); @v };

for my $view (@view)
  {
    &find ({wanted => \&wanted, no_chdir => 1}, "src/$view/");
  }

'FileHandle'->new (">f2f.pl")->print (&Dumper (\%f2f));


