#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my ($F90_1, $F90_2) = @ARGV;

my $d1 = &Fxtran::fxtran (location => $F90_1, fopts => [qw (-line-length 800)]);
my $d2 = &Fxtran::fxtran (location => $F90_2, fopts => [qw (-line-length 800)]);


my %k = qw (UVH 0 XYB 1 RCP 1 CTY 0);

my @dummy = &F ('//dummy-arg-LT/arg-N', $d2);
my @remove = map { 0 } @dummy;

my ($name) = &F ('//subroutine-N', $d2, 1);
print "name=$name\n";

for my $call (&F ('.//call-stmt[string(procedure-designator)="?"]', $name, $d1))
  {
    my @actual = &F ('./arg-spec/arg/*', $call);

    my %d2a;
    for my $i (0 .. $#actual)
      {
#       next if ($dummy[$i]->textContent =~ m/^Y/o);
        my $actual = $actual[$i]->textContent;
        next if ($actual =~ m/\(/o);
#       next unless ($actual =~ m/^(?:YDCPG_DYN0|YDCPG_DYN9|YDCPG_PHY0|YDCPG_PHY9|YDCPG_TND|YDCPG_MISC)%/o);
#       next unless ($actual =~ m/^(?:YLAPLPAR%)/o);
#       next unless ($actual =~ m/^(?:YDMF_PHYS%TMP%RDT%)/o);
#       next unless ($actual =~ m/^(?:YDCPG_DYN|YDCPG_PHY|YDMF_PHYS_SURF|YDVARS|YDMF_PHYS|YDCPG_MISC|YLMF_PHYS_STATE)/o);
#       next unless ($actual =~ m/^(?:YDCPG_MISC|YDCPG_PHY0|YDCPG_PHY9|YDMF_PHYS|YDCPG_DYN0|YDCPG_DYN9|YDMF_PHYS_SURF|YDVARS|YDMF_PHYS_TMP|YDAPLPAR_TMP)/o);
#       next unless ($actual =~ m/^(?:YDCFU|YDCSGEOM|YDGSGEOM|YDOROG)%/goms);
#       next unless ($actual =~ m/^(?:YLMF_PHYS_BASE_STATE)/goms);
#       next unless ($actual =~ m/^(?:YDCPG_DYN|YDCPG_PHY|YDVARS)/goms);
#       next unless ($actual =~ m/^(?:YDCPG_DIM|YDMF_PHYS|YDCPG_MISC)/goms);
#       next unless ($actual =~ m/^(?:YDCPG_OPTS)/goms);
        next unless ($actual =~ m/^(?:YDMODEL)/goms);
#       next unless ($actual =~ m/^(?:YDVARS)/goms);
#       next unless ($actual =~ m/^(?:YDCPG_BNDS|YDCPG_OPTS|YDMF_PHYS_BASE_STATE|YDCPG_MISC|YDMF_PHYS_SURF|YDMF_PHYS)/goms);
        $actual =~ s/^YL/YD/o;
        $d2a{$dummy[$i]->textContent} = $actual;
        if ($actual[$i]->parentNode->nextSibling)
          {
            $actual[$i]->parentNode->nextSibling->unbindNode;
          }
        elsif ($actual[$i]->parentNode->previousSibling)
          {
            $actual[$i]->parentNode->previousSibling->unbindNode;
          }
        $actual[$i]->parentNode->unbindNode;
        $remove[$i] = 1;
      }

print &Dumper (\%d2a);

    while (my ($k, $v) = each (%d2a))
      {
        my @expr = &F ('.//named-E[string(N)="?"]', $k, $d2);
        for my $expr (@expr)
          {
            my ($N) = &F ('./N', $expr); my $NN = $N->textContent;
            my ($r) = &F ('./R-LT/ANY-R', $expr);
            if ($r && grep { ("P${_}0" eq $NN) || ("P${_}9" eq $NN) } keys (%k))
              {
                my ($n) = ($NN =~ m/^P(\w+)\d$/o);
                if ($r->nodeName eq 'parens-R')
                  {
                    my @e = &F ('./element-LT/element/ANY-E', $r);
                    if (($e[0]->textContent eq '1') && ($e[1]->textContent eq '1'))
                      {
                        $e[0]->replaceNode (&t (':'));
                        if ($k{$n} eq 1)
                          {
                            $e[1]->replaceNode (&t (':'));
                          }
                        else
                          {
                            $e[1]->replaceNode (&t ('1:'));
                          }
                      }
                    elsif (($e[0]->textContent eq '1') && ($e[1]->textContent eq '0'))
                      {
                        $e[0]->replaceNode (&t (':'));
                        if ($k{$n} eq 0)
                          {
                            $e[1]->replaceNode (&t (':'));
                          }
                        else
                          {
                            die $expr->textContent;
                          }
                      }
                  }
                else
                  {
                    die $r->textContent;
                  }
              }
            $N->replaceNode (&t ($v));
          }
      }

  }

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


for my $i (0 .. $#dummy)
  {
    next unless ($remove[$i]);
    next if ($dummy[$i]->textContent =~ m/^Y/o);

    &removeListElement ($dummy[$i]);

    my $N = $dummy[$i]->textContent;
    my ($en_decl) = &F ('.//EN-decl[string(EN-N)="?"]', $N, $d2);
    unless ($en_decl)
      {
#       'FileHandle'->new (">$F90_2.new")->print ($d2->textContent);
        die "N=$N\n";
      }

    my $stmt = &Fxtran::stmt ($en_decl);

    if (&removeListElement ($en_decl))
      {
        $stmt->unbindNode;
      }
  }


'FileHandle'->new (">$F90_1.new")->print ($d1->textContent);
'FileHandle'->new (">$F90_2.new")->print ($d2->textContent);

