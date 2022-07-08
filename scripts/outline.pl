#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

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


sub resolveAssociates
{
  my $node = shift;
  

AGAIN:
  my @N = &F ('.//named-E/N', $node);

  my $count = 0;
    
  for my $N (@N)
    {
      my ($expr) = $N->parentNode;
      my @s = &F ('ancestor::associate-construct/associate-stmt/associate-LT/associate[associate-N[string(.)="?"]]/selector/*', $N->textContent, $node);
      die if (@s > 1);
      next unless (my $s = $s[0]);
      $s = $s->cloneNode (1);
      my @r = &F ('./R-LT/*', $expr);
      my ($rlt) = &F ('./R-LT', $s);
      for my $r (@r)
        {
          $rlt->appendChild ($r->cloneNode (1));
        }
      $expr->replaceNode ($s);
      $count++;
    }

  goto AGAIN if ($count);


}

sub sortArgs
{
  my ($d, @args) = @_;

  my @dum = &F ('.//dummy-arg-LT/arg-N', $d, 1);

  my %dum = map { ($dum[$_], $_) } (0 .. $#dum);

  my @d = grep {   defined ($dum{$_}) } @args;
  my @l = grep { ! defined ($dum{$_}) } @args;

  @d = sort { $dum{$a} <=> $dum{$b} } @d;
  @l = sort { $a cmp $b } @l;

  return (@d, @l);
}

my @INTRINSIC = qw (SIGN MAX MIN MOD JPRB REAL EXP ASIN FOLH SQRT PRESENT);
my %INTRINSIC = map { ($_, 1) } @INTRINSIC;


my $F90 = shift;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 400)]);

my @C = &F ('//C[starts-with(string (.),"!* outline")]', $d);

for my $C (@C)
  {
    (my $SUB = $C->textContent) =~ s/^!\* outline\s+//o;
    my $sub = lc ($SUB);

    my $o = &Fxtran::fxtran (string => << "EOF");
SUBROUTINE $SUB ()

USE PARKIND1, ONLY : JPIM, JPRB

!

IMPLICIT NONE

!

END SUBROUTINE
EOF

    my ($include) = &F ('.//include[last()]', $d);
    if ($include)
      {
        $include->parentNode->insertAfter (&n ("<include>#include &quot;<filename>$sub.intfb.h</filename>&quot;</include>"), $include);
        $include->parentNode->insertAfter (&t ("\n"), $include);
      }


    my @node;
    for (my $node = $C; ; $node = $node->nextSibling)
      {
        $node or die;
        push @node, $node;
        if (($node->nodeName eq 'C') && (index ($node->textContent, '!* end outline') == 0))
          {
            last;
          }
      }

    for my $node (@node)
      {
        &resolveAssociates ($node);
      }
    
    my (%N, %do, %call, %loc);
    for my $node (@node)
      {
        my @expr = &F ('.//named-E', $node);
        for my $expr (@expr)
          {
            my ($n) = &F ('./N/n', $expr, 1);
            next if ($INTRINSIC{$n});
            if ($expr->parentNode->nodeName eq 'do-V')
              {
                $do{$n} = 1;
              }
            elsif ($expr->parentNode->nodeName eq 'procedure-designator')
              {
                $call{$n} = 1;
              }
            $N{$n} = 1;
          }
      }

    my @N = &sortArgs ($d, grep { (! $do{$_}) && (! $call{$_}) } keys (%N));

    for my $N (@N)
      {
        my @expr_d = &F ('.//named-E[string(N)="?"]', $N, $d); 
        my @expr_n;
        for my $node (@node)
          {
            push @expr_n, &F ('.//named-E[string(N)="?"]', $N, $node);
          }
        $loc{$N} = @expr_d > @expr_n ? 0 : 1;
        if (my ($decl) = &F ('.//T-decl-stmt[.//EN-decl[string(EN-N)="?"]]', $N, $d))
          {
             my ($intent) = &F ('./attribute[string(attribute-N)="INTENT"]', $decl);
             if ($intent)
               {
                 $loc{$N} = 0;
               }
          }
        else
          {
            die $N;
          }
      }

    my %rank = map { ($N[$_], $_) } (0 .. $#N);

    @N = sort { ($loc{$a} <=> $loc{$b}) or ($rank{$a} <=> $rank{$b}) } @N;


    my %N2M;
    
    my %S = qw (Z P LL LD I K J K N K L LD YL YD R P Y Y);
    my @S = qw (Z LL I J N L YL R Y);

    my %M;

    for my $N (@N)
      {
        for my $v (values (%S))
          {
            if (index ($N, $v) == 0)
              {
                $N2M{$N} = $N;
                $M{$N} = 1;
                last;
              }
          }
      }
    
    for my $N (@N)
      {
        next if ($N2M{$N});
        for my $k (@S)
          {
            if (index ($N, $k) == 0)
              {
                my $v = $S{$k};
                my $M = $N;
                $M =~ s/^$k/$v/ unless ($loc{$N});
                while ($M{$M})
                  {
                    $M .= '_';
                  }
                $N2M{$N} = $M;

                last;
              }
          }
      }


    $C->replaceNode (&n ("<call-stmt>CALL <procedure-designator><named-E><N><n>$SUB</n></N></named-E></procedure-designator> (<arg-spec>"  
                       . join (', ', map { "<arg><named-E><N><n>$_</n></N></named-E></arg>" } grep { ! $loc{$_} } @N) . '</arg-spec>)</call-stmt>'));

    my ($C1, $C2) = &F ('.//C', $o);
    
    $node[-1]->unbindNode ();
    shift (@node); pop (@node);
   
    for my $node (reverse (@node))
      {
        $C2->parentNode->insertAfter ($node, $C2);
      }

    my ($dummy_arg_LT) = &F ('.//dummy-arg-LT', $o);
    for my $N (@N)
      {

        $N2M{$N} or die $N;

        if (! $loc{$N})
          {
            $dummy_arg_LT->appendChild (&n ("<arg-N><N><n>$N2M{$N}</n></N></arg-N>"));
            $dummy_arg_LT->appendChild (&t (", ")) if ($N ne $N[-1]);
          }

        my @n = &F ('.//named-E[string(N)="?"]/N/n/text()', $N, $o);
        for my $n (@n)
          {
            $n->replaceNode (&t ($N2M{$N}));
          }

        if (my ($decl) = &F ('.//T-decl-stmt[.//EN-decl[string(EN-N)="?"]]', $N, $d))
          {
            $decl = $decl->cloneNode (1);
            my ($n) = &F ('.//EN-decl[string(EN-N)="?"]/EN-N/N/n/text()', $N, $decl);
            my @ts = &F ('.//derived-T-spec/T-N', $decl, 1);

            my @en_decl = &F ('.//EN-decl', $decl);

            for my $en_decl (@en_decl)
              {
                my ($EN_N) = &F ('./EN-N', $en_decl, 1);
                &removeListElement ($en_decl) if ($EN_N ne $N);
              }


            for my $ts (@ts)
              {
                my ($use) = &F ('.//use-stmt[.//use-N[string(.)="?"]]', $ts, $d);
                $use = $use->cloneNode (1);

                $C1->parentNode->insertBefore ($use->cloneNode (1), $C1);
                $C1->parentNode->insertBefore (&t ("\n"), $C1);
              }

            $n->replaceNode (&t ($N2M{$N}));
            $C2->parentNode->insertBefore ($decl, $C2);
          }
        else
          {
            if ($N =~ m/^[KIJMN]/o)
              {
                $C2->parentNode->insertBefore (&Fxtran::fxtran (statement => "INTEGER (KIND=JPIM), INTENT (IN) :: $N2M{$N}"), $C2);
              }
            elsif ($N =~ m/^L/o)
              {
                $C2->parentNode->insertBefore (&Fxtran::fxtran (statement => "LOGICAL, INTENT (IN) :: $N2M{$N}"), $C2);
              }
            elsif ($N =~ m/^R/o)
              {
                $C2->parentNode->insertBefore (&Fxtran::fxtran (statement => "REAL (KIND=JPRB), INTENT (IN) :: $N2M{$N}"), $C2);
              }
            else
              {
                die $N;
              }
          }
        $C2->parentNode->insertBefore (&t ("\n"), $C2);

      }

    $C2->parentNode->insertBefore (&t ("\n"), $C2);

    for my $call (sort keys (%call))
      {
        $C2->parentNode->insertBefore (&Fxtran::fxtran (statement => '#include "' . lc ($call) . '.intfb.h"'), $C2);
        $C2->parentNode->insertBefore (&t ("\n"), $C2);
      }

    $C2->parentNode->insertBefore (&t ("\n"), $C2);

    for my $v (sort keys (%do))
      {
        $C2->parentNode->insertBefore (&Fxtran::fxtran (statement => "INTEGER (KIND=JPIM) :: $v"), $C2);
        $C2->parentNode->insertBefore (&t ("\n"), $C2);
      }

    $C1->unbindNode ();
    $C2->unbindNode ();

    for my $N (@N)
      {
        next unless ($loc{$N});
        my ($en_decl) = &F ('.//EN-decl[string(EN-N)="?"]', $N, $d);
        die $N unless ($en_decl);
        my $stmt = &Fxtran::stmt ($en_decl);
        my ($intent) = &F ('./attribute[string(attribute-N)="INTENT"]', $stmt);
        unless ($intent)
          {
            $stmt->unbindNode () if (&removeListElement ($en_decl));
          }
      }

    my $dir = &dirname ($F90);
    'FileHandle'->new (">$dir/$sub.F90")->print ($o->textContent);
  }


'FileHandle'->new (">$F90.new")->print ($d->textContent);





