#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my $F90 = shift;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);

my @r = &F ('.//named-E/R-LT/component-R[string(ct)="ZVIEW"]', $d);

for my $r (@r)
  {
    next unless (my $ra = $r->nextSibling);
    if ($ra->nodeName eq 'array-R')
      {
        my ($x) = &F ('./section-subscript-LT/section-subscript[last()]/lower-bound/ANY-E', $ra);
        if ($x->textContent =~ m/(\d)(?:_PHY)?%M_(\w+)$/o)
          {
            my ($t, $f) = ($1, $2);
            $t = "T$t";
            my $expr = $r->parentNode->parentNode;
            $x->parentNode->parentNode->previousSibling->unbindNode;
            $x->parentNode->parentNode->unbindNode;
            my ($ct) = &F ('./ct', $r);
            $ct->replaceNode (&t ("$f"));
            my @x = &F ('./section-subscript-LT/section-subscript/node()', $ra);
            unless (grep { $_->textContent ne ':' } @x)
              {
                $ra->unbindNode ();
              }
          }
        else
          {
            my $expr = $r->parentNode->parentNode;
            die $expr->textContent;
          }
      }
    elsif ($ra->nodeName eq 'parens-R')
      {
        my ($x) = &F ('./element-LT/element[last()]/ANY-E', $ra);
        if ($x->textContent =~ m/(\d)(?:_PHY)?%M_(\w+)$/o)
          {
            my ($t, $f) = ($1, $2);
            $t = "T$t";
            my $expr = $r->parentNode->parentNode;
            $x->parentNode->previousSibling->unbindNode;
            $x->parentNode->unbindNode;
            my ($ct) = &F ('./ct', $r);
            $ct->replaceNode (&t ("$f"));
          }
        else
          {
            my $expr = $r->parentNode->parentNode;
            die $expr->textContent;
          }
      }
    else
      {
        die $ra->textContent;
      }
  }


'FileHandle'->new (">$F90.new")->print ($d->textContent);

