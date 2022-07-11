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

my (%g, %L, %D);

my @q = @ARGV;

my %skip1 = map { ($_, 1) } qw (DR_HOOK ABOR1_SFX ABOR1 WRSCMR SC2PRG VERDISINT VEXP MPL_ALLREDUCE NEW_ADD_FIELD_3D ADD_FIELD_3D FMLOOK_LL
                                FMWRIT LES_MEAN_SUBGRID SECOND_MNH ABORT_SURF SURF_INQ SHIFT ABORT LFAECRR FLUSH LFAFER LFAECRI LFAOUV LFAECRC
                                LFAPRECR NEW_ADD_FIELD_2D ADD_FIELD_2D WRITEPROFILE WRITEMUSC WRITEPHYSIO CONVECT_SATMIXRATIO COMPUTE_FRAC_ICE
                                TRIDIA LFAPRECI PPP2DUST GET_LUOUT ALLOCATE DEALLOCATE PUT GET SAVE_INPUTS GET_FRAC_N DGEMM SGEMM ARO_GROUND_PARAM ARO_GROUND_DIAG);

#%skip1 = (%skip1, map { ($_, 1) } qw (GPHLWI GPHLUV GPHPRE SITNU SIGAM SISEVE CP_PTRSLB1 ETENC));

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
       grep { ! m/GSTAT/o } grep { ! m/^MPL_/o } grep { ! m/^JFH_/o }
       grep { ! $skip1{$_} } &F ('//call-stmt/procedure-designator', $d, 1);

    for my $call (@call)
      {
        &add ($call);
        push @{ $g{$name} }, $call;
        $D{$call} = $D{$name} + 1;

        push @q, $call unless ($seen{$call}++);
      }
    
  }


#print &Dumper (\%ARP);

my @root;

while (my ($k, $v) = each (%g))
  {
    my %seen;
    @$v = grep { ! ($seen{$_}++) } grep { $g{$_} } @$v;
  }

@root = keys (%g);

while (my ($k, $v) = each (%g))
  {
    my %v = map { ($_, 1) } @$v;
    @root = grep { ! $v{$_} } @root;
  }

sub color
{
  my $name = shift;

  return ();
}

my $root = join ('-', sort @root);

my $g = 'GraphViz2'->new (graph => {rankdir => 'LR', ordering => 'out'}, global => {rank => 'source'});
#my $g = 'GraphViz2'->new (graph => {rankdir => 'TB', ordering => 'out'}, global => {rank => 'source'});
#

my %k;

while (my ($k, $v) = each (%g))
  {
     next if ($k =~ m/%/o);
    $g->add_node (name => $k, label => "$k\n$L{$k}", shape => 'box', &color ($k));
#   $g->add_node (name => $k, label => "$k", shape => 'box', &color ($k));
    for my $v (@$v)
      {   
        next if ($v =~ m/%/o);
        $g->add_edge (from => $k, to => $v);
#print "$k -- $v\n";
        $k{$k}++;
        $k{$v}++;
      }   
  }

for (sort keys (%k))
  {
    print $f2f->{$_} || $_, "\n";
  }

$g->run (format => 'svg', output_file => "$root.svg");


