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

my @assoc = &F ('.//associate-construct', $d);

for my $assoc (@assoc)
  {
    my $stmt = $assoc->firstChild;

    my @N = &F ('./associate-LT/associate/associate-N', $stmt);

    for my $N (@N)
      {
        my $n = $N->textContent;
        my @expr = &F ('.//named-E[string(N)="?"]', $n, $d);
        next if (@expr);
        my $associate = $N->parentNode;
        if ($associate->nextSibling)
          {
            $associate->nextSibling->unbindNode;
            $associate->unbindNode;
          }
        elsif ($associate->previousSibling)
          {
            $associate->previousSibling->unbindNode;
            $associate->unbindNode;
          }
      }

  }

for my $cnt (&F ('.//cnt', $d))
  {
    next unless ($cnt->previousSibling);
    next unless ($cnt->nextSibling);
    next unless ($cnt->previousSibling->nodeName eq '#text');
    next unless ($cnt->previousSibling->data =~ m/^\s*\n\s*$/o);
    next unless ($cnt->nextSibling->nodeName eq '#text');
    next unless ($cnt->nextSibling->data =~ m/^\s*$/o);
    next unless ($cnt->nextSibling->nextSibling->nodeName eq 'cnt');
    my $sp = $cnt->previousSibling->data;
    $sp =~ s/\n\s*$//o;
    $cnt->previousSibling->setData ($sp);
    $cnt->nextSibling->nextSibling->unbindNode;
    $cnt->nextSibling->unbindNode;
    $cnt->unbindNode;
  }


'FileHandle'->new (">$F90.new")->print ($d->textContent);
