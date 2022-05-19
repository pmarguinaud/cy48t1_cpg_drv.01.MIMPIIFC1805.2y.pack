package SingleBlock;
#
use strict;
use FileHandle;
use Data::Dumper;

use Fxtran;


sub hoistJlonLoops
{
  my $doc = shift;
  
  my @pu = &f ('./f:object/f:file/f:program-unit', $doc);
  
  for my $pu (@pu)
    {
      # Find DO loops with JLON variable
  
      my @do = &f ('.//f:do-construct[./f:do-stmt/f:do-V/f:named-E/f:N/f:n/text ()="JLON"]', $pu);

      my %lh;
  
      # For each JLON DO loop, find the outermost loop enclosing this JLON DO loop when possible
  
      for my $do (@do)
        {
          my ($doo) = &f ('ancestor-or-self::f:do-construct[./f:do-stmt/f:do-V/f:named-E/f:N/f:n/text ()!="JLON"]', $do);
          $lh{$doo->unique_key} = $doo if ($doo);
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
<do-construct><do-stmt>DO <do-V><named-E><N><n>JLON</n></N></named-E></do-V> = <lower-bound><named-E><N><n>KIDIA</n></N></named-E></lower-bound>, <upper-bound><named-E><N><n>KFDIA</n></N></named-E></upper-bound></do-stmt>
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
  
          my @do = &f ('descendant-or-self::f:do-construct[./f:do-stmt/f:do-V/f:named-E/f:N/f:n/text ()="JLON"]', $lh);
  
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
