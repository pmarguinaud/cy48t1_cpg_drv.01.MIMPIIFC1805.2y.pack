#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

sub slurp
{
  my $f = shift;
  my $data = do { my $fh = 'FileHandle'->new ("<$f"); local $/ = undef; <$fh> };
  return $data;
}

sub eqf
{
  &slurp ($_[0]) eq &slurp ($_[1]);
}

my ($F90) = @ARGV;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 800)]);

my @v = qw (
LAPRXPK              
LGWADV               
LNHDYN               
LNHX                 
LNHXDER              
LPC_FULL             
LRPLANE              
LRUBC                
LSLAG                
LTWOTL               
LVERTFE              
LVFE_ECMWF           
LVFE_GW              
LVFE_LAPL_BC         
L_RDRY_VD            
NDLNPR               
NPDVAR               
NVDVAR               
RHYDR0               
TOPPRES              
);

sub dummy
{
  my $n = shift;

  if ($n =~ s/^L/LD/o)
    {
      return $n;
    }

  if ($n =~ s/^N/K/o)
    {
      return $n;
    }
  
  return $n;
}

sub type
{
  my $n = shift;

  if ($n =~ m/^L/o)
    {
      return 'LOGICAL';
    }

  if ($n =~ m/^[NK]/o)
    {
      return 'INTEGER (KIND=JPIM)';
    }
  
  return 'REAL (KIND=JPRB)';
}

my @args;

for my $v (reverse (@v))
  {
    my ($un) = &F ('.//use-stmt//use-N[string(.)="?"]', $v, $d);

    next unless ($un);

    $un->unbindNode ();

    my ($dlt) = &F ('.//dummy-arg-LT', $d);
  
    my ($arg0) = &F ('.//arg-N', $dlt, 1);
  
    unshift (@args, $v);

    my $dum = &dummy ($v);
    my $type = &type ($v);

    $dlt->insertBefore (&t (', '), $dlt->firstChild);
    $dlt->insertBefore (&n ("<arg-N><N><n>$dum</n></N></arg-N>"), $dlt->firstChild);
  
    my ($stmt) = &F ('.//T-decl-stmt[.//EN-N[string(.)="?"]]', $arg0, $d);

    my $stmt1 = &Fxtran::fxtran (statement => "$type, INTENT (IN) :: $dum");
    $stmt->parentNode->insertBefore ($stmt1, $stmt);
    $stmt->parentNode->insertBefore (&t ("\n"), $stmt);
  
  
    my @expr = &F ('.//named-E[string(N)="?"]/N/n/text()', $v, $d);
    for my $expr (@expr)
      {
        $expr->setData ($dum);
      }
  }    

my @use = &F ('.//use-stmt[not(.//use-N)]', $d);

for my $use (@use)
  {
    $use->unbindNode ();
  }

    
'FileHandle'->new (">$F90.new")->print ($d->textContent);
    
my @f = (qw (src/local/arpifs/adiab/cpg_gp.F90 
             src/local/arpifs/adiab/cpg_gp_hyd.F90
             src/local/arpifs/adiab/cpg_gp_nhee.F90
             src/local/arpifs/adiab/cpg_gp_nhqe.F90
             src/local/arpifs/adiab/cpg_gp_hyd_t9.F90
             src/local/arpifs/adiab/cpg_gp_nhee_t9.F90
             src/local/arpifs/adiab/cpg_gp_nhqe_t9.F90
             src/local/arpifs/adiab/gnhee_refine_grp.F90
             src/local/arpifs/adiab/gpcty.F90
             src/local/arpifs/adiab/gpgrgeo.F90
             src/local/arpifs/adiab/gpgrp.F90
             src/local/arpifs/adiab/gpgrxyb.F90
             src/local/arpifs/adiab/gphpre.F90
             src/local/arpifs/adiab/gpmpfc_expl.F90
             src/local/arpifs/adiab/gpmpfc_pgfl.F90
             src/local/arpifs/pp_obs/pos.F90
            ));

(my $sub = &basename ($F90)) =~ s/\.F90$//o;
$sub = uc ($sub);

print "args=@args\n";

for my $f (@f)
  {
    my $d = &Fxtran::fxtran (location => $f, fopts => [qw (-line-length 800)]);
    my @call = &F ('.//call-stmt[string(procedure-designator)="?"]', $sub, $d);

    my $k = 0;
    for my $call (@call)
      {
        my ($argspec) = &F ('./arg-spec', $call);
        for my $arg (reverse (@args))
          {
            $argspec->insertBefore (&t (', '), $argspec->firstChild);
            $argspec->insertBefore (&n ("<arg><named-E><N><n>$arg</n></N></named-E></arg>"), $argspec->firstChild);
          }
        $k++;
      }

    'FileHandle'->new (">$f.new")->print ($d->textContent) if ($k);

  }
    
