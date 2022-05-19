package Vector;

use Fxtran;
use strict;

sub hoistJlonLoops
{
  my $doc = shift;
  
  my @pu = &f ('./f:object/f:file/f:program-unit', $doc);
  
  my ($JLON, $KIDIA, $KFDIA) = ('JLON', 'KIDIA', 'KFDIA');

  for my $pu (@pu)
    {
      # Find DO loops with JLON variable
  
      my @do = &f ('.//f:do-construct[./f:do-stmt/f:do-V/f:named-E/f:N/f:n/text ()="?"]', $JLON, $pu);

      my %lh;
  
      # For each JLON DO loop, find the outermost loop enclosing this JLON DO loop when possible
  
      for my $do (@do)
        {
          my @doo = &f ('ancestor-or-self::f:do-construct[./f:do-stmt/f:do-V/f:named-E/f:N/f:n/text ()!="?"]', $JLON, $do);


          my $doo;

          for my $dooo (@doo)
            {
              my ($var) = &f ('./f:do-stmt/f:do-V', $dooo, 1);
              $doo = $dooo if ($var && ($var ne $JLON) && ($var ne 'JITER'));
              last if ($var && ($var ne 'JITER'));
            }


          if ($doo)
            {
              my ($write) = &f ('.//f:write-stmt', $doo);
              $lh{$doo->unique_key} = $doo unless ($write);
            }
          
        }

      my @lh = values (%lh);

      for my $lh (@lh)
        {
          # Find current level of indentation
  
          my $sp = $lh->previousSibling;
          $sp = $sp->textContent;
          $sp =~ s/^\s*\n//o;
  
          # Create a JLON loop nest
  
          my @dob = &n (<< "EOF");
<do-construct><do-stmt>DO <do-V><named-E><N><n>$JLON</n></N></named-E></do-V> = <lower-bound><named-E><N><n>$KIDIA</n></N></named-E></lower-bound>, <upper-bound><named-E><N><n>$KFDIA</n></N></named-E></upper-bound></do-stmt>
$sp<C/>
$sp<end-do-stmt>ENDDO</end-do-stmt></do-construct>
EOF
  
          # Inject the nest before the outermost loop
  
          for my $dob (@dob)
            {
              $lh->parentNode->insertBefore ($dob, $lh);
            }
  
          # Re-nest outermost loop
          my ($dob)= &f ('.//f:C', $dob[0]);
  
          $dob->replaceNode ($lh);
  
          # Remove innermost JLON DO loops
  
          my @do = &f ('descendant-or-self::f:do-construct[./f:do-stmt/f:do-V/f:named-E/f:N/f:n/text ()="?"]', $JLON, $lh);
  
          for my $do  (@do)
            {
              $do->firstChild->unbindNode ();
              $do->lastChild->unbindNode ();
              # Use a pseudo target to remove the loop construct
              my $C = &n ('<C/>');
              $do->replaceNode ($C);
              for my $c ($do->childNodes ())
                {
                  $C->parentNode->insertBefore ($c, $C);
                }
              $C->unbindNode ();
            } 
        }

    }
  
}

sub addDirectives
{
  my $doc = shift;
  
  my @pu = &f ('./f:object/f:file/f:program-unit', $doc);
  
  my ($JLON, $KIDIA, $KFDIA) = ('JLON', 'KIDIA', 'KFDIA');

  for my $pu (@pu)
    {
      my @do = &f ('.//f:do-construct[./f:do-stmt/f:do-V/f:named-E/f:N/f:n/text ()="?"]', $JLON, $pu);
  
      for my $do (@do)
        {
          my $sp = &Fxtran::getIndent ($do);
          $do->parentNode->insertBefore (&n ('<C>!$acc loop vector</C>'), $do);
          $do->parentNode->insertBefore (&t ("\n" . (' ' x $sp)), $do);
        }
      
      my ($subroutine) = &f ('./f:subroutine-stmt', $pu);

      my $sp = &Fxtran::getIndent ($subroutine);
      my ($name) = &f ('./f:subroutine-N/f:N/f:n/text ()', $subroutine, 1);
      $subroutine->parentNode->insertBefore (&n ("<C>!\$acc routine ($name) vector</C>"), $subroutine);
      $subroutine->parentNode->insertBefore (&t ("\n" . (' ' x $sp)), $subroutine);

    }
  
}



1;
