#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use List::MoreUtils qw (uniq);
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

sub ct2a
{
  my $ct = shift;
  (my $a1 = $ct) =~ s/%/_/go; $a1 = "P$a1";
  return $a1;
}


my $d0 = &Fxtran::fxtran (location => $ARGV[0], fopts => [qw (-line-length 400)]);
my $d1 = &Fxtran::fxtran (location => $ARGV[1], fopts => [qw (-line-length 400)]);

my $s = 'YDAPLPAR_TMP';

my %decl = %{ do ("./h.pl") };

# Work on second document

my ($arg1) = &F ('.//subroutine-stmt/dummy-arg-LT/arg-N[string(.)="?"]', $s, $d1);
my ($dec1) = &F ('.//T-decl-stmt[.//EN-decl[string(EN-N)="?"]]', $s, $d1);
my ($typ1) = &F ('./_T-spec_/derived-T-spec/T-N', $dec1, 1);


my ($n1) = &F ('.//subroutine-N', $d1, 1);

my @expr = &F ('.//named-E[string(N)="?"]', $s, $d1);

my %ct;

for my $expr (@expr)
  {
    my @ct = &F ('./R-LT/component-R/ct', $expr);
    my $ct = join ('%', map { $_->textContent } @ct);
    $ct{$ct} = 1;
    for (@ct)
      {
        $_->parentNode->unbindNode ();
      }
    my ($n) = &F ('./N/n/text()', $expr);
    my $a1 = &ct2a ($ct);
    $n->replaceNode (&t ($a1));
  }

my @ct = sort keys (%ct);

for my $ct (reverse (@ct))
  {
    my $a1 = &ct2a ($ct);
    $arg1->parentNode->insertAfter (&n ("<arg-N><N><n>$a1</n></N></arg-N>"), $arg1);
    $arg1->parentNode->insertAfter (&t (", "), $arg1);
   
    (my $decl = $decl{"${typ1}%${ct}"}) or die ("Cannot find ${typ1}%${ct}");
    $decl = &Fxtran::fxtran (statement => $decl);

    my ($n) = &F ('.//EN-N/N/n/text()', $decl);
    $n->replaceNode (&t ("$a1"));
 
    $dec1->parentNode->insertAfter ($decl, $dec1);
    $dec1->parentNode->insertAfter (&t ("\n"), $dec1);
    
  }

&removeListElement ($arg1);
$dec1->unbindNode ();

# Work on first document

my ($call) = &F ('.//call-stmt[string(procedure-designator)="?"]', $n1, $d0);

my ($arg0) = &F ('./arg-spec/arg[string(.)="?"]', $s, $call);


for my $ct (reverse (@ct))
  {
    my $a0 = $arg0->cloneNode (1);
    my $e0 = $a0->firstChild;
    $e0->appendChild (my $rlt = &n ('<R-LT/>'));
    for my $c (split (m/%/o, $ct))
      {
        $rlt->appendChild (&n ("<component-R>%<ct>$c</ct></component-R>"));
      }
    $arg0->parentNode->insertAfter ($a0, $arg0);
    $arg0->parentNode->insertAfter (&t (', '), $arg0);
  }

&removeListElement ($arg0);

'FileHandle'->new (">$ARGV[0].new")->print ($d0->textContent);
'FileHandle'->new (">$ARGV[1].new")->print ($d1->textContent);

