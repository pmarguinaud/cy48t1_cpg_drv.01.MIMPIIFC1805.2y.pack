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
my $mod2var = do './mod2var.pl';

my (%g, %L, %D, %U, %O, %V);

my @q = ('APL_ARPEGE');

sub add
{
  my $name = shift;
  return if (exists $g{$name});
  $g{$name} ||= [];
  my $file = $f2f->{$name};
  unless ($file)
    {
      $L{$name} = '';
      return;
    }
  my @code = do { my $fh = 'FileHandle'->new ("<$file"); <$fh> };
  $L{$name} = scalar (@code);
}

for my $q (@q)
  {
    &add ($q);
    $D{$q} = 0;
  }


my %seen;

while (my $name = shift (@q))
  {
#   next if ($D{$name} > 2);

    my $file = $f2f->{$name};

    next unless ($file);

    my $d = &Fxtran::fxtran (location => $file, fopts => [qw (-line-length 300)], dir => '/tmp');
    
    my @call = 
       grep { ! m/GSTAT/o } grep { ! m/^MPL_/o } grep { ! m/^JFH_/o } grep { ! /DR_HOOK/o }
       grep { ! m/ARO_GROUND_/o } grep { $_ ne 'VEXP' } grep { $_ ne 'ABOR1' }
       &F ('//call-stmt/procedure-designator', $d, 1);

    my @use = &F ('.//use-stmt', $d);
    my %use;

    my %skip = map { ($_, 1) } qw (LHOOK);
    for my $use (@use)
      {
        my ($mod) = &F ('./module-N', $use, 1);
        my @ren = &F ('.//rename', $use);
        $O{$name}{$mod}++ unless (@ren);
        for my $ren (@ren)
          {
            my ($N) = &F ('./N', $ren, 1);
            ($N) = &F ('./use-N', $ren, 1) unless ($N);
            $use{$mod}{$N}++;

            if ($mod2var->{$mod}{$N} && (! $skip{$N}))
              {
                my @expr = &F ('.//named-E[string(N)="?"]', $N, $d);
                $V{$file}{$mod}{$N} = scalar (@expr);
              }
          }
      }

    $U{$name} = \%use;

    for my $call (@call)
      {
        &add ($call);
        push @{ $g{$name} }, $call;
        $D{$call} = $D{$name} + 1;

        push @q, $call unless ($seen{$call}++);
      }
    
  }


die if (%O);

# {local $Data::Dumper::Terse = 1, print &Dumper (\%V); }
#

for my $file (sort keys (%V))
  {
    printf ("%-40s ", &basename ($file));
    for my $mod (sort keys (%{ $V{$file} }))
      {
        printf (" | %s : %s", $mod, join (', ', sort keys (%{ $V{$file}{$mod} })));
      }
    print "\n";
  }



