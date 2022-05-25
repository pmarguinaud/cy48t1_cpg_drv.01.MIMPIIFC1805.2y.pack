package SingleDirective;
#
use strict;
use FileHandle;
use Data::Dumper;

use Fxtran;


sub hoistJlonLoops
{
  my $par = shift;
  
  my ($do_jblk) = &F ('./do-construct[./do-stmt[string(do-V)="JBLK"]]', $par);
  die $par unless ($do_jblk);

  # Find DO loops with JLON variable & remove them

  my @do = &F ('.//do-construct[./do-stmt/do-V/named-E/N/n/text()="JLON"]', $par);

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

  $par->normalize ();

  my ($update) = &F ('./call-stmt[string(procedure-designator)="YLCPG_BNDS%UPDATE"]', $do_jblk);
  my @e = &F ('following-sibling::node()', $update);
  pop (@e) for (1 .. 2); # ENDDO JBLK

  my $indent = &Fxtran::getIndent ($update);

  my ($dob) = &Fxtran::fxtran (fragment => << "EOF");
DO JLON = YLCPG_BNDS%KIDIA, YLCPG_BNDS%KFDIA
  ! coucou
ENDDO
EOF
  
 
  &Fxtran::reIndent ($dob, $indent);
  $do_jblk->insertAfter ($dob, $update);
  $do_jblk->insertAfter (&t ("\n" . (' ' x $indent)), $update);


  my ($C)= &F ('.//C', $dob);
  for my $e (reverse (@e))
    {
      &Fxtran::reIndent ($e, $indent);
      $dob->insertAfter ($e, $C);
      $dob->insertAfter (&t (' ' x ($indent)), $C);
    }

  $C->unbindNode ();

  $par->normalize ();

}

sub addParallelLoopDirectives
{
  my $d = shift;

  my @pu = &f ('./f:object/f:file/f:program-unit', $d);
  
  # Insert OpenACC parallel directives
  
  for my $pu (@pu)
    {
      my @do = &f ('.//f:do-construct[./f:do-stmt/f:do-V/f:named-E/f:N/f:n/text()="JLON"]', $d);
  
      for my $do (@do)
        {
          my %p;
  
          # Loop variables & temporary scalars; these are meant to be private
          my @s = &f ('.//f:E-1/f:named-E[not (./f:R-LT)]/f:N/f:n/text ()', $do);
          my @v = &f ('descendant-or-self::f:do-construct/f:do-stmt/f:do-V/f:named-E/f:N/f:n/text()', $do);
       
          for (@s, @v)
            {
              $p{$_->textContent}++;
            }
  
          my @p = sort keys (%p);
  
          my $sp = $do->previousSibling;
          ($sp = $sp->textContent) =~ s/^\s*\n//o;
          $do->parentNode->insertBefore (&n ('<C>!$acc parallel loop gang vector '
#                                          . (@p ? 'private (' . join (', ', @p) . ') ' : '')
#                                          . 'default (none)'
                                           . '</C>'), $do);
          $do->parentNode->insertBefore (&t ("\n"), $do);
          $do->parentNode->insertBefore (&t ($sp), $do);
        }
  
    }

}

1;
