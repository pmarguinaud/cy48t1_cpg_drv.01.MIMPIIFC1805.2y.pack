#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w

use strict;
use FindBin qw ($Bin);
use lib $Bin;
use FileHandle;
use File::Copy;
use File::Basename;
use File::stat;
use File::Path;
use Getopt::Long;

my $suffix = '_openacc';

sub newer
{
  return 1;
  my ($f1, $f2)  = @_;
  die unless (-f $f1);
  return 1 unless (-f $f2);
  return stat ($f1)->mtime > stat ($f2)->mtime;
}

sub copyIfNewer
{
  my ($f1, $f2) = @_;

  if (&newer ($f1, $f2))
    {
      print "Copy $f1 to $f2\n"; 
      &copy ($f1, $f2); 
    }
}

sub saveToFile
{
  my ($x, $f) = @_;

  unless (-d (my $d = &dirname ($f)))
    {
      &mkpath ($d);
    }

  'FileHandle'->new (">$f")->print ($x->textContent ());
  'FileHandle'->new (">$f.xml")->print ($x->toString ());
}

sub addSeqDirective
{
  my $d = shift;
  my ($pu) = &F ('./object/file/program-unit', $d);
  my ($N) = &F ('./subroutine-stmt/subroutine-N', $pu, 1);
  $pu->insertAfter (&n ("<C>!\$acc routine ($N) seq</C>"), $pu->firstChild);
  $pu->insertAfter (&t ("\n"), $pu->firstChild);
}

sub addSuffix
{
  my ($d, $suffix) = @_;

  $suffix = uc ($suffix);

  my @pde = &F ('.//call-stmt/procedure-designator[not(string(.)="DR_HOOK")]/N/n/text()', $d);
  my @sub = &F ('.//subroutine-N/N/n/text()', $d);

  for (@pde, @sub)
    {
      $_->setData ($_->data . $suffix);
    }

  my @inc = &F ('.//include/filename/text()', $d);

  $suffix = lc ($suffix);

  for (@inc)
    {
      my $f = $_->data;
      $f =~ s/\.intfb\.h$/$suffix.intfb.h/goms;
      $_->setData ($f);
    }

}

sub prune
{
  use Apply;

  my $d = shift;

  &Apply::apply 
  (
    $d, 
    '//named-E[string(N)="LMUSCLFA"]' => &n ('<literal-E>.FALSE.</literal-E>'),
    '//named-E[string(N)="YDLDDH"][./R-LT/component-R[string(ct)="LFLEXDIA"]]', &n ('<literal-E>.FALSE.</literal-E>'),
    '//named-E[string(N)="YDLDDH"][./R-LT/component-R[string(ct)="LFLEXDIA"]]', &n ('<literal-E>.FALSE.</literal-E>'),
  );
  

}

sub preProcessIfNewer
{
  use Inline;
  use Associate;
  use Fxtran;
  use Blocks;
  use SingleBlock;
  use Vector;
  use Stack;
  use Loop;
  use ReDim;
  use DrHook;
  use SymbolTable;
  use Construct;

  my ($f1, $f2) = @_;

  my $conf = 
  {
    ind => ['JLON', 'JLEV'],
    dim2ind => {'KLON' => 'JLON', 'KLEV' => 'JLEV'},
    ind2bnd => {'JLON' => ['KIDIA', 'KFDIA'], 'JLEV' => ['1', 'KLEV']},
  };

  if (&newer ($f1, $f2))
    {
      print "Preprocess $f1\n";

      my $d = &Fxtran::fxtran (location => $f1);
      &saveToFile ($d, "tmp/$f2");

      my @hdr = &F ('//include/filename', $d);
      my @pde = &F ('//procedure-designator', $d, 1);

      my $t = &SymbolTable::getSymbolTable ($d, {NPROMA => 'KLON'});

      &Inline::inlineContainedSubroutines ($d);
      &saveToFile ($d, "tmp/inlineContainedSubroutines/$f2");

      &Associate::resolveAssociates ($d);
      &saveToFile ($d, "tmp/resolveAssociates/$f2");

      &Construct::changeIfStatementsInIfConstructs ($d);
      &saveToFile ($d, "tmp/changeIfStatementsInIfConstructs/$f2");

      &prune ($d);
      &saveToFile ($d, "tmp/prune/$f2");

      &Loop::arraySyntaxLoop ($d, $t, $conf);
      &saveToFile ($d, "tmp/arraySyntaxLoop/$f2");

      &Loop::removeJlonLoops ($d);
      &saveToFile ($d, "tmp/removeJlonLoops/$f2");

      &ReDim::reDim ($d);
      &saveToFile ($d, "tmp/reDim/$f2");

      &addSeqDirective ($d);

      &Stack::addStack ($d);
      &saveToFile ($d, "tmp/addStack/$f2");

      for my $hdr (@hdr)
        {
          (my $f = $hdr->textContent) =~ s/\.intfb\.h//o;
          $f = uc ($f);
          next unless (grep { $_ eq $f } @pde);
          my $inc = $hdr->parentNode;
          $inc->unbindNode ();
        }

      &addSuffix ($d, $suffix);
      &saveToFile ($d, "tmp/addSuffix/$f2");

      &DrHook::remove ($d);
      &saveToFile ($d, "tmp/DrHook/$f2");

      'FileHandle'->new (">$f2")->print ($d->textContent ());

    }
}


my ($f, $g) = @ARGV;

unless ($g)
  {
    $g = 'src/local/openacc/' . &basename ($f);
    $g =~ s/\.F90$/$suffix.F90/;
  }

&preProcessIfNewer ($f, $g);

&Fxtran::intfb ($g, 'src/local/openacc');


