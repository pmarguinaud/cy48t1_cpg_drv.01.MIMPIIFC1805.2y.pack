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

my ($model) = &F ('.//dummy-arg-LT//arg-N[string(.)="YDMODEL"]', $d);


for my $func (qw (fctdoi fcttrm))
  {
    my @inc = &F ('.//include/filename[string(.)="?"]/text()', "$func.func.h", $d);
    
    for my $inc (@inc)
      {
        $inc->setData ("$func.ycst.h");
      }
  }

if ($model)
  {

    my @cst = &F ('.//use-stmt[string(module-N)="YOMCST"]//use-N', $d, 1);

    for my $cst (@cst)
      {
        my @expr = &F ('.//named-E[string(N)="?"]/N/n/text()', $cst, $d);
        for my $expr (@expr)
          {
            $expr->setData ("YDCST%$cst");
          }
      }
    
    my @use = &F ('.//use-stmt[string(module-N)="YOMCST"]', $d);
    for (@use)
      {
        $_->unbindNode ();
      }
 
    my ($hook1, $hook2) = &F ('.//if-stmt[.//call-stmt[string(procedure-designator)="DR_HOOK"]]', $d);

    $hook1->parentNode->insertAfter (&t ("\n\nASSOCIATE (YDCST => YDMODEL%YRCST)"), $hook1);
    $hook2->parentNode->insertBefore (&t ("END ASSOCIATE\n\n"), $hook2);

  }
else
  {
    my ($dlt) = &F ('.//dummy-arg-LT', $d);

    my ($arg0) = &F ('.//arg-N', $dlt, 1);

    $dlt->insertBefore (&t (', '), $dlt->firstChild);
    $dlt->insertBefore (&n ("<arg-N><N><n>YDCST</n></N></arg-N>"), $dlt->firstChild);

    my ($stmt) = &F ('.//T-decl-stmt[.//EN-N[string(.)="?"]]', $arg0, $d);
    $stmt->parentNode->insertBefore (&t ("TYPE (TCST), INTENT (IN) :: YDCST\n"), $stmt);


    my @cst = &F ('.//use-stmt[string(module-N)="YOMCST"]//use-N', $d, 1);
    
    for my $cst (@cst)
      {
        my @expr = &F ('.//named-E[string(N)="?"]/N/n/text()', $cst, $d);
        for my $expr (@expr)
          {
            $expr->setData ("YDCST%$cst");
          }
      }
    
    my ($rlt) = &F ('.//use-stmt[string(module-N)="YOMCST"]/rename-LT', $d);
    
    for ($rlt->childNodes)
      {
        $_->unbindNode ();
      }
    
    $rlt->appendChild (&n ("<rename><use-N><N><n>TCST</n></N></use-N></rename>"));

  }


'FileHandle'->new (">$F90.new")->print ($d->textContent);


