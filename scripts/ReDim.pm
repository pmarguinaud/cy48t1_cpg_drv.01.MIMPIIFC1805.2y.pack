package ReDim;
#
use strict;
use FileHandle;
use Data::Dumper;

use FindBin qw ($Bin);
use lib $Bin;
use Fxtran;

sub reDim
{
  my $d = shift;
  
  my @en_decl = &F ('.//EN-decl[./array-spec/shape-spec-LT[string(shape-spec)="?"]]', 'KLON', $d);
  
  EN_DECL : for my $en_decl (@en_decl)
    {
      my ($N) = &F ('./EN-N', $en_decl, 1);
      my ($stmt) = &Fxtran::stmt ($en_decl);

      next if (&F ('.//attribute-N[string(.)="INTENT"]', $stmt));
      next if (&F ('.//call-stmt[.//named-E[string(N)="?"]', $N, $d));
  
# The following is used sometimes (eg acbl89.F90)
#
# ZTESTM=SUM(ZTESTSAVE(KIDIA:KFDIA))
#
# DO JJLEV=JLEV,KTDIAN+1,-1
#   IF (ZTESTM > 0.0_JPRB) THEN
#     DO JLON=KIDIA,KFDIA
#       ZDLUP1   = ZGDZF(JLON,JJLEV)
#       ZZTHVL   =(ZTHETA (JLON,JJLEV)+ZTHETA (JLON,JJLEV-1))/2.0_JPRB
#       ...
#       ZEN   (JLON)=ZEN  (JLON)-ZINCR*ZTEST0
#     ENDDO
#     ZTESTM=SUM(ZTESTSAVE(KIDIA:KFDIA))
#   ENDIF
# ENDDO

      my @sum = &F ('//E-2/named-E[translate(string(.)," ","")="?"]', "SUM($N(KIDIA:KFDIA))", $d);

# So we need to replace the SUM by a scalar assignment

      my $N_JLON = &e ("$N(JLON)");
      for my $sum (@sum)
        {
          $sum->replaceNode ($N_JLON->cloneNode (1));
        }


      my ($as) = &F ('./array-spec', $en_decl);

      my @ss = &F ('./shape-spec-LT/shape-spec', $as);

      for my $ss (@ss[1..$#ss])
        {
          my ($lb) = &F ('./lower-bound/*', $ss);
          my ($ub) = &F ('./upper-bound/*', $ss);
          for ($lb, $ub)
            {
              next unless ($_);
              next EN_DECL unless ($_->nodeName eq 'literal-E');
            }
        }

      if (@ss > 1)
        {
          $ss[0]->nextSibling->unbindNode ();
          $ss[0]->unbindNode ();  
        }
      else
        {
          $as->unbindNode ();
        }

      my @rlt = &F ('.//named-E[string(N)="?"]/R-LT', $N, $d);
  
      for my $rlt (@rlt)
        {
          my ($r) = &F ('./ANY-R[1]', $rlt);
          unless ($r->nodeName eq 'parens-R')
            {
              my $expr = $rlt->parentNode;
              die $expr->textContent;
            }
          my @e = &F ('./element-LT/element', $r);
          if (scalar (@e) > 1)
            {
              $e[0]->nextSibling->unbindNode ();
              $e[0]->unbindNode ();
            }
          else
            {
              $rlt->unbindNode ();
            }
        }
  
    }
}

1;
