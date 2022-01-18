#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

sub parseDirectives
{
  my $d = shift;

  my @C = &F ('//C[starts-with(string (.),"!=")]', $d);
  
  while (my $C  = shift (@C))
    {
      (my $bdir = $C->textContent) =~ s/^!=\s*//o;
      $bdir = lc ($bdir);
      my ($tag) = ($bdir =~ m/^(\w+)/o);
  
      my @node;
      for (my $node = $C->nextSibling; ; $node = $node->nextSibling)
        {
          $node or die;
          if (($node->nodeName eq 'C') && (index ($node->textContent, '!=') == 0))
            {
              my $C = shift (@C);
              (my $edir = $C->textContent) =~ s/^!=\s*//o;
              $edir = lc ($edir);

              die unless ($edir =~ s/^end\s+//o);
              die unless ($edir eq $tag);

              $C->unbindNode ();
              
              last;
            }

          push @node, $node;

        }

      my $e = &n ("<$tag-directive/>");
 
      for my $node (@node)
        {
          $e->appendChild ($node);
        }

      $C->replaceNode ($e);

    }
}

sub parallel
{
  shift;
  my ($par, $doc) = @_;

  my $indent = ' ' x &Fxtran::getIndent ($par);

  my ($loop) = &Fxtran::fxtran (fragment => << "EOF");
DO JBLK = 1, YDCPG_DIM%KGPBLKS
${indent}ENDDO
EOF

  my ($enddo) = &F ('.//end-do-stmt', $loop);
  my $p = $enddo->parentNode;

  for my $node ($par->childNodes ())
    {
      $p->insertBefore (&t (' ' x ($indent + 2)), $enddo);
      &Fxtran::reIndent ($node, $indent + 2);
      $p->insertBefore ($node, $enddo);
    }
  
  $par->replaceNode ($loop);

}

sub skip
{
  shift;
  $_[0]->unbindNode ();
}

sub processDirectives
{
  my $doc = shift;

  my @dir = &F ('.//ANY-directive', $doc);

  for my $dir (@dir)
    {
      (my $name = $dir->nodeName) =~ s/-directive$//o;
      __PACKAGE__->can ($name) or die ("Unknown directive $name");
      __PACKAGE__->$name ($dir, $doc);
    }

}


my $F90 = shift;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);

my @name = &F ('.//subroutine-N/N/n/text()', $d);

for (@name)
  {
    $_->setData ($_->textContent . '_OPENMP');
  }

&parseDirectives ($d);

&processDirectives ($d);

print $d->textContent;
    
    
