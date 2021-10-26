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


for my $en_decl (&F ('.//EN-decl', $d))
  {

    my ($N) = &F ('./EN-N', $en_decl, 1);

    print "$N\n";

    my @expr = &F ('.//named-E[string(N)="?"]', $N, $d);

    next if (@expr);
   
    if ($en_decl->nextSibling)
      {
        $en_decl->nextSibling->unbindNode;
        $en_decl->unbindNode;
      }
    elsif ($en_decl->previousSibling)
      {
        $en_decl->nextSibling->unbindNode;
        $en_decl->unbindNode;
      }
    else
      {
        my $stmt = &Fxtran::stmt ($en_decl);
        $stmt->unbindNode;
      }

    my ($arg) = &F ('.//dummy-arg-LT/arg-N[string(.)="?"]', $N, $d);
    
    if (!$arg)
      {
      }
    elsif ($arg->nextSibling)
      {
        $arg->nextSibling->unbindNode;
        $arg->unbindNode;
      }
    elsif ($arg->previousSibling)
      {
        $arg->previousSibling->unbindNode;
        $arg->unbindNode;
      }



  }

'FileHandle'->new (">$F90.new")->print ($d->textContent);

