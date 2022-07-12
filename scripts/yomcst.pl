#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my ($F90) = @ARGV;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 800)]);

my @s = qw (
  GP_KAPPAT
  LATTES
  LATTE_NL
  LAVENT
  SIDD
  SIGAM
  SILKAP
  SINHQE_SEVE
  SINHQE_VDERI
  SIPTP
  SISKAP
  SITNU
  ESPCSI
  ESPNHEESI
  ESPNHQESI
  ESPNHSI
  SILKAPI
  SPCSI
  SPNHEESI
  SPNHQESI
  SPNHSI
  SUNHEEBMAT
  SUNHQEBMAT
  SUNHQESI
  ESPCSI
  ESPNHEESI
  ESPNHQESI
  ESPNHSI
  SI_CCCOR
  SILKAPI
  SPCSI
  SPNHEESI
  SPNHQESI
  SPNHSI
  SUBMAT
  SUDYN_STABILITY
  SUNHEEBMAT
  SUNHQEBMAT
  SUNHQESI
  SUNHSI_TESTCONV
  SUPRECOV

SUNHEESI
ESPCSI
SUSI
SUNHBMAT
SUNHEEBMAT
SUNHSI

LATTE_KAPPA
CPEULDYN


);



for my $ms (['YOMCST', 'CST'], ['YOETHF', 'THF'])
  {
    my ($mod, $str) = @$ms;
    my ($rlt) = &F ('.//use-stmt[string(module-N)="?"]/rename-LT', $mod, $d);

    next unless ($rlt);

    for my $func (qw (fcttre fctdoi fcttrm))
      {
        my @inc = &F ('.//include/filename[string(.)="?"]/text()', "$func.func.h", $d);
        
        for my $inc (@inc)
          {
            $inc->setData ("$func.ycst.h");
          }
      }
    
    my ($dlt) = &F ('.//dummy-arg-LT', $d);
  
    my ($arg0) = &F ('.//arg-N', $dlt, 1);
  
    $dlt->insertBefore (&t (', '), $dlt->firstChild);
    $dlt->insertBefore (&n ("<arg-N><N><n>YD$str</n></N></arg-N>"), $dlt->firstChild);
  
    my ($stmt) = &F ('.//T-decl-stmt[.//EN-N[string(.)="?"]]', $arg0, $d);

    my $stmt1 = &Fxtran::fxtran (statement => "TYPE (T$str), INTENT (IN) :: YD$str");
    $stmt->parentNode->insertBefore ($stmt1, $stmt);
    $stmt->parentNode->insertBefore (&t ("\n"), $stmt);
  
  
    my @cst = &F ('.//use-stmt[string(module-N)="?"]//use-N', $mod, $d, 1);
    
    for my $cst (@cst)
      {
        my @expr = &F ('.//named-E[string(N)="?"]/N/n/text()', $cst, $d);
        for my $expr (@expr)
          {
            $expr->setData ("YD$str%$cst");
          }
      }
    
    
    for ($rlt->childNodes)
      {
        $_->unbindNode ();
      }
    
    $rlt->appendChild (&n ("<rename><use-N><N><n>T$str</n></N></use-N></rename>"));
    
  }    


for my $s (@s)
  {
    my @call = &F ('.//call-stmt[string(procedure-designator)="?"]', $s, $d);
    for my $call (@call)
      {
        my ($argspec) = &F ('./arg-spec', $call);
        my ($arg) = &F ('./arg', $argspec);
        next if ($arg->textContent eq 'YDCST');
        $argspec->insertBefore (&t (', '), $argspec->firstChild);
        $argspec->insertBefore (&n ('<arg><named-E><N><n>YDCST</n></N></named-E></arg>'), $argspec->firstChild);
      }
  }


    
'FileHandle'->new (">$F90.new")->print ($d->textContent);
    
    
