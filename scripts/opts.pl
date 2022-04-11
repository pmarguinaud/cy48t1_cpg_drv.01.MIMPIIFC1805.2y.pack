#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;
use Associate;

my $F90 = shift;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 400)]);

my @da = &F ('.//dummy-arg-LT/arg-N[string(N)="YDCPG_DIM"]', $d);

exit (0) unless (@da);

for my $da (@da)
  {
    my $lt = $da->parentNode;
    $lt->insertAfter (&n ("<arg-N><N><n>YDCPG_OPTS</n></N></arg-N>"), $da);
    $lt->insertAfter (&t (', '), $da);
  }


my @use = &F ('.//use-stmt[string(module-N)="CPG_DIM_TYPE_MOD"][./rename-LT/rename[string(use-N)="CPG_DIM_TYPE"]]', $d);

for my $use (@use)
  {
    my ($rlt) = &F ('./rename-LT', $use);
    $rlt->appendChild (&t (', '));
    $rlt->appendChild (&n ("<rename><use-N><N><n>CPG_OPTS_TYPE</n></N></use-N></rename>"));
  }

my @aa = &F ('.//call-stmt/arg-spec/arg[string(.)="YDCPG_DIM"]', $d);

for my $aa (@aa)
  {
    $aa->parentNode->insertAfter (&n ("<arg><named-E><N><n>YDCPG_OPTS</n></N></named-E></arg>"), $aa);
    $aa->parentNode->insertAfter (&t (', '), $aa);
  }

my @en_decl = &F ('.//EN-decl[string(EN-N)="YDCPG_DIM"]', $d);

for my $en_decl (@en_decl)
  {
    my ($stmt) = &Fxtran::stmt ($en_decl);
    my $decl = $stmt->cloneNode (1);
    $stmt->parentNode->insertAfter ($decl, $stmt);
    $stmt->parentNode->insertAfter (&t ("\n"), $stmt);
    my ($N) = &F ('.//EN-N/N/n/text()', $decl);
    $N->replaceNode (&t ("YDCPG_OPTS"));
    my ($T) = &F ('./_T-spec_/derived-T-spec/T-N/N/n/text()', $decl);
    $T->replaceNode (&t ("CPG_OPTS_TYPE"));
    my ($ts) = &F ('./_T-spec_', $decl); $ts = $ts->nextSibling;
    my $tt = $ts->textContent;
    if ($tt =~ s/^,\s/,/o)
      {
        $ts->setData ($tt);
      }
    else
      {
        my ($lt) = &F ('./EN-decl-LT', $decl);
        my $ss = $lt->previousSibling;
        $tt = $ss->textContent;
        $tt =~ s/\s::/::/o;
        $ss->setData ($tt);
      }
    
  }

'FileHandle'->new (">$F90.new")->print ($d->textContent);

