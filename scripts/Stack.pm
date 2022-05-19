package Stack;

use Fxtran;
use strict;
use Data::Dumper;

sub removeListElement
{
  my $x = shift;

  my $nn = $x->nodeName;

  my ($p) = $x->parentNode;
  
  my @cf = &F ('following-sibling::text()[contains(.,",")]', $x);   
  my @cp = &F ('preceding-sibling::text()[contains(.,",")]', $x);   
  
  if (@cf)
    {   
      $cf[+0]->unbindNode (); 
    }   
  elsif (@cp)
    {   
      $cp[-1]->unbindNode (); 
    }   
  
  $x->parentNode->appendChild (&t (' '));
  my $l = $x->parentNode->lastChild;
  
  $x->unbindNode (); 
  
  while ($l)
    {   
      last if (($l->nodeName ne '#text') && ($l->nodeName ne 'cnt'));
      $l = $l->previousSibling;
      last unless ($l);
      $l->nextSibling->unbindNode;
    }   

  return &F ("./$nn", $p) ? 0 : 1;
}



sub addStack
{
  my $d = shift;

  my @call = &F ('.//call-stmt[string(procedure-designator)!="ABOR1" and string(procedure-designator)!="REDUCE"]', $d);

  for my $call (@call)
    {
      my ($argspec) = &F ('./arg-spec', $call);
      $argspec->appendChild (&t (', '));
      $argspec->appendChild (&n ("<named-E><N><n>YLSTACK</n></N></named-E>"));
    }

  my ($dummy_arg_lt) = &F ('.//subroutine-stmt/dummy-arg-LT', $d);

  my ($last) = &F ('./arg-N[last()]', $dummy_arg_lt, 1);

  $dummy_arg_lt->appendChild (&t (', '));
  $dummy_arg_lt->appendChild (&n ("<arg-N><N><n>YDSTACK</n></N></arg-N>"));

  my ($use) = &F ('.//use-stmt[last()]', $d);
  $use->parentNode->insertAfter (&n ("<include>#include &quot;<filename>stack.h</filename>&quot;</include>"), $use);
  $use->parentNode->insertAfter (&t ("\n"), $use);
  $use->parentNode->insertAfter (&n ("<use-stmt>USE <module-N><N><n>STACK_MOD</n></N></module-N></use-stmt>"), $use);
  $use->parentNode->insertAfter (&t ("\n"), $use);


  my ($decl) = &F ('.//T-decl-stmt[.//EN-N[string(.)="?"]]', $last, $d);
  $decl->parentNode->insertAfter (&n (
'<T-decl-stmt><_T-spec_><derived-T-spec>TYPE(<T-N><N><n>STACK</n></N></T-N>)</derived-T-spec></_T-spec_> ' .
':: <EN-decl-LT><EN-decl><EN-N><N><n>YDSTACK</n></N></EN-N></EN-decl>' .
', <EN-decl><EN-N><N><n>YLSTACK</n></N></EN-N></EN-decl></EN-decl-LT></T-decl-stmt>'), $decl);
  $decl->parentNode->insertAfter (&t ("\n"), $decl);

  
  my ($noexec) = do 
  {
    my ($exec) = grep { &Fxtran::stmt_is_executable ($_) } &F ('.//ANY-stmt', $d);
    my @prev = &F ('preceding::*', $exec);

    my $prev;
    for my $p (reverse (@prev))
      {
        next if ($p->nodeName eq '#text');
        next if ($p->nodeName eq 'C');
        $prev = $p;
        last;
      }

    my @anc = &F ('ancestor::*', $prev);

    for my $anc (reverse (@anc))
      {
        if (($anc->nodeName =~ m/-(?:construct|stmt)$/o) || ($anc->nodeName eq 'include'))
          {
            $prev = $anc;
          }
      }

    $prev
  };

  my $C = &n ("<C/>");

  $noexec->parentNode->insertAfter (&t ("\n"), $noexec);
  $noexec->parentNode->insertAfter ($C, $noexec);

  $C->parentNode->insertBefore (&t ("\n"), $C);
  $C->parentNode->insertBefore (&t ("\n"), $C);
  $C->parentNode->insertBefore (&n (
'<a-stmt><E-1><named-E><N><n>YLSTACK</n></N></named-E></E-1><a>=</a>' .
'<E-2><named-E><N><n>YDSTACK</n></N></named-E></E-2></a-stmt>'), $C);
  $C->parentNode->insertBefore (&t ("\n"), $C);
  $C->parentNode->insertBefore (&t ("\n"), $C);

  my @KLON = qw (KLON);

  for my $KLON (@KLON)
    {
      my @en_decl = &F ('.//T-decl-stmt[not(string(.//attribute-N)="INTENT")]'
                      . '//EN-decl[./array-spec/shape-spec-LT/shape-spec[string(./upper-bound)="?"]]', 
                      $KLON, $d);
      
      for my $en_decl (@en_decl)
        {
          my ($n) = &F ('./EN-N', $en_decl, 1);

          my $stmt = &Fxtran::stmt ($en_decl);
      
          my ($t) = &F ('./_T-spec_',   $stmt);     &Fxtran::expand ($t); $t = $t->textContent;
          my ($s) = &F ('./array-spec', $en_decl);  &Fxtran::expand ($s); $s = $s->textContent;
      
          $stmt->parentNode->insertBefore (my $temp = &t ("temp ($t, $n, $s)"), $stmt);
      
          if (&removeListElement ($en_decl))
            {
              $stmt->unbindNode ();
            }
          else
            {
              $temp->parentNode->insertAfter (&t ("\n"), $temp);
            }
      
          $C->parentNode->insertBefore (&t ("alloc ($n)\n"), $C);

        }

    }

  $C->unbindNode ();

}

1;
