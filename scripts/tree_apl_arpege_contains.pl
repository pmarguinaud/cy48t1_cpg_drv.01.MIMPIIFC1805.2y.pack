#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;
use GraphViz2;

my $f2f = do './f2f.pl';


my @q = ('APL_ARPEGE');

sub add
{
  my $name = shift;
  my $file = $f2f->{$name};
  unless ($file)
    {
      return;
    }
}

for my $q (@q)
  {
    &add ($q);
  }


my %seen;

while (my $name = shift (@q))
  {
#   next if ($D{$name} > 2);

    my $file = $f2f->{$name};

    next unless ($file);

    my $d = &Fxtran::fxtran (location => $file, fopts => [qw (-line-length 300)], dir => '/tmp');
    
    my @contains = &F ('.//contains-stmt', $d);

    if (@contains)
      {
        print $file, "\n";
      }

    my @call = 
       grep { ! m/GSTAT/o } grep { ! m/^MPL_/o } grep { ! m/^JFH_/o } grep { ! /DR_HOOK/o }
       grep { ! m/ARO_GROUND_/o } grep { $_ ne 'VEXP' } grep { $_ ne 'ABOR1' }
       &F ('//call-stmt/procedure-designator', $d, 1);

    for my $call (@call)
      {
        &add ($call);
        push @q, $call unless ($seen{$call}++);
      }
    
    
  }



