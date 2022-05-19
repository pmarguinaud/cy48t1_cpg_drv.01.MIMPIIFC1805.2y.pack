package Loop;
#
use strict;
use FileHandle;
use Data::Dumper;

use FindBin qw ($Bin);
use lib $Bin;
use Fxtran;

my $conf0 = 
{
  ind => ['JLON', 'JLEV'],
  dim2ind => {'YDCPG_OPTS%KLON' => 'JLON', 'YDCPG_OPTS%KFLEVG' => 'JLEV'},
  ind2bnd => {'JLON' => ['YDCPG_BNDS%KIDIA', 'YDCPG_BNDS%KFDIA'], 'JLEV' => ['1', 'YDCPG_OPTS%KFLEVG']},
  block => 1,
};

sub arraySyntaxLoop
{
  my ($d, $t, $conf) = @_;

  $conf ||= $conf0;

  # Transform simple array syntax in loops
  
  my @stmt = &F ('.//a-stmt[./E-1/named-E[not(R-LT)]]|.//a-stmt[./E-1/named-E/R-LT/array-R]', $d);

  for my $stmt (@stmt)
    {
      my ($N1) = &F ('./E-1/named-E/N', $stmt, 1);
      my $s1 = $t->{$N1};

      next unless ($s1->{nproma});

      my %ind;
      my @ind = @{ $conf->{ind} };
      my %dim2ind = %{ $conf->{dim2ind} };
      my %ind2bnd = %{ $conf->{ind2bnd} };

      # For all expressions in this statement, involving arrays whose first dimension if NIT, replace subscript sections by indices

      my @bb;

      my @expr = &F ('.//named-E', $stmt);
      for my $expr (@expr)
        {
          my ($N) = &F ('./N', $expr, 1);
          my $s = $t->{$N};

          next unless ($s->{nproma});

          my $as = $s->{as};

          my @ss1 = &F ('./shape-spec-LT/shape-spec', $as);

          my ($r) = &F ('./R-LT/array-R', $expr);

          # The array reference does not exist; add it
          
          unless ($r)
            {
              my $rlt = &n ('<R-LT/>');
              $expr->appendChild ($rlt);
              $r = &n ('<array-R/>');
              $rlt->appendChild ($r);
              my $sslt = &n ('<section-subscript-LT/>');
              $r->appendChild (&t ('('));
              $r->appendChild ($sslt);
              $r->appendChild (&t (')'));
              for my $i (0 .. $#ss1)
                {
                  $sslt->appendChild (&n ('<section-subscript>:</section-subscript>'));
                  $sslt->appendChild (&t (', ')) if ($i != $#ss1);
                }
            }

          $r->setNodeName ('parens-R');
          my ($sslt) = &F ('./section-subscript-LT', $r);
          $sslt->setNodeName ('element-LT');

          my @ss2 = &F ('./section-subscript', $sslt);
 
          if ($conf->{block})
            {
              die unless ('JBLK' ne pop (@ss2));
            }

          die &Dumper ([$expr->textContent, [map { $_->textContent } @ss1], [map { $_->textContent } @ss2]]) . ' ' unless (scalar (@ss1) == scalar (@ss2));

          for my $i (0 .. $#ss1)
            {
              my ($lb1) = &F ('./lower-bound', $ss1[$i], 1); $lb1 = 1 unless (defined ($lb1));
              my ($ub1) = &F ('./upper-bound', $ss1[$i], 1);

              $ss2[$i]->setNodeName ('element');
              if ($ss2[$i]->textContent eq ':') 
                {
                  my $ind = $dim2ind{$ub1};
                  die $expr->textContent . ' ' . $ss1[$i]->textContent unless ($ind);
                  $ind{$ind}++;
                  $ss2[$i]->firstChild->replaceNode (&n ("<named-E><N><n>$ind</n></N></named-E>"));

                  if ($i > 0)
                    {
                      $ind2bnd{$ind}[0] = $lb1;
                      $ind2bnd{$ind}[1] = $ub1;
                    }

                }
              elsif ($ss2[$i]->textContent =~ m/:/o)
                {
                  die &Dumper ([$stmt->textContent, $sslt->textContent, $i, $ss2[$i]->textContent]) . ' ';
                }
            }
          

        }

      @ind = grep { $ind{$_} } @ind;

      my $indent = &Fxtran::getIndent ($stmt);

      # Generate the loop and inserts the assignment statement 

      my $loop = join ('', 
                        map ({ my $i = $_; my $ind = $ind[$#ind-$i]; ('  ' x $i) . "DO $ind = $ind2bnd{$ind}[0], $ind2bnd{$ind}[1]\n" } (0 .. $#ind)),
                        ('  ' x scalar (@ind)) . "!\n",
                        reverse (map ({ my $i = $_; ('  ' x $i) . "ENDDO\n" } (0 .. $#ind)))
                      );

      ($loop) = &Fxtran::fxtran (fragment => $loop);

      &Fxtran::reIndent ($loop, $indent);

      my ($C) = &F ('.//C', $loop);
      $C->replaceNode ($stmt->cloneNode (1));

      $stmt->replaceNode ($loop);



    }

}

sub removeJlonLoops
{
  my $d = shift;
  
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
  
  
  $noexec->parentNode->insertAfter (&Fxtran::fxtran (statement => "JLON = KIDIA\n"), $noexec);
  $noexec->parentNode->insertAfter (&t ("\n"), $noexec);
  $noexec->parentNode->insertAfter (&t ("\n"), $noexec);
  
  
  my @do = &F ('.//do-construct[./do-stmt[string(do-V)="JLON"]]', $d);
  
  for my $do (@do)
    {
      $do->firstChild->unbindNode;
      $do->lastChild->unbindNode;
      my @nodes = &F ('./node()', $do);
      for (@nodes)
        {
          $do->parentNode->insertBefore ($_, $do);
        }
      $do->unbindNode ();
    }
  
}

1;
