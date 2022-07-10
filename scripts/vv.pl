#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FileHandle;
use FindBin qw ($Bin);
use lib $Bin;
use Fxtran;
use File::Basename;
use Data::Dumper;
use List::MoreUtils qw (uniq);

my @view = do { my $fh = 'FileHandle'->new ('<.gmkview'); <$fh> };
chomp for (@view);

my @dir = qw (arpifs/adiab aladin/coupling);

sub resolve
{
  my $f = shift;
  for my $view (@view)
    {
      for my $dir (@dir)
        {
          my $g = "src/$view/$dir/$f";
          return $g if (-f $g);
        }
    }
  die $f;
}

my @skip = qw (DR_HOOK CPG_GP_NHEE_T9 CPG_GP_NHQE_T9 ABOR1 CPG_GP_HYD_T9 SC2PRG CP_FORCING UPDATE_GWDIAG GPPOPER UPDSST ETENC GPINOZST);
my %skip = map { ($_, 1) } @skip;
my @type = qw (MODEL_GENERAL_CONF_TYPE CPG_BNDS_TYPE CPG_DYN_TYPE CPG_OPTS_TYPE CPG_TND_TYPE DR_HOOK FIELD_VARIABLES GEOMETRY 
               JPIM JPRB JPRD LHOOK MODEL TCST SETDEFAULTL TDIMV TDPHY TDYN TEPHY TGMV TPHY2 TVAB TVERTICAL_GEOM TVETA TVFE TOZO TPTRSLB2 
               TSIMPHL TSPNG TYPE_GFLD);
my %type = map { ($_, 1) } @type;


my @call;

my @f = (qw (src/local/arpifs/adiab/cpg_gp.F90 
             src/local/arpifs/adiab/cpg_gp_hyd.F90
             src/local/arpifs/adiab/cpg_gp_nhee.F90
             src/local/arpifs/adiab/cpg_gp_nhqe.F90));
for my $f (@f)
  {
    my $d = &Fxtran::fxtran (location => $f, fopts => [qw (-line-length 300)]);
    push @call, &F ('.//call-stmt/procedure-designator', $d, 1);
  }

@call = &uniq (sort grep { ! $skip{$_} } @call);

my $fh = 'FileHandle'->new (">list.txt");

my %use;

for my $call (@call)
  {
    $call = lc ($call);
    my $f = &resolve ("$call.F90");

    next if (grep { $_ eq $f } @f);

    my $k = 0;

    my $d = &Fxtran::fxtran (location => $f, fopts => [qw (-line-length 300)]);
    
    my @use = &F ('.//use-stmt', $d);

    for my $use (@use)
      {
        my ($N) = &F ('./module-N', $use, 1);
        next if (($N eq 'YOMLSFORC') or ($N eq 'YOMGWDIAG') or ($N eq 'YOMCT3'));
        my @N = &F ('.//use-N', $use, 1);
        for (@N)
          {
            if ($use{$_}) 
              {
                die unless ($use{$_} eq $N);
              }
            $use{$_} = $N;
            $k++ unless ($type{$_});
          }
      }

    $fh->print ("$f\n") if ($k);

  }

$fh->close ();

my %mod;

for my $N (sort keys (%use))
  {
    next if ($type{$N});
    printf ("%-20s %-20s\n", $N, $use{$N});
    push @{ $mod{$use{$N}} }, $N;
  }

__END__

for my $mod (sort keys (%mod))
  {
    my @v = sort @{ $mod{$mod} };
    print "USE $mod, ONLY : " . join (', ', @v)  . "\n";
  }








