#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FileHandle;
use FindBin qw ($Bin);
use lib $Bin;
use Fxtran;
use File::Basename;
use Data::Dumper;
use List::MoreUtils qw (uniq);
use List::Util qw (min max);

my $F90 = shift;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);

my $t0 = $d->textContent;

my @use = &F ('.//use-stmt', $d);

for my $use (@use)
  {
    my $count = 0;
    my @n = &F ('.//use-N/N/n', $use);
    for my $n (@n)
      {
        my @expr = &F ('.//named-E[string(N)="?"]', $n->textContent, $d);
        my @type = &F ('.//T-N[string(.)="?"]', $n->textContent, $d);
        if (@expr || @type)
          {
            $count++;
            next;
          }
        my ($rename) = &F ('ancestor::rename', $n);
        $rename->unbindNode ();
      }
    $use->unbindNode unless ($coount);
  }

@use = &F ('.//use-stmt', $d);


my %use;

for my $use (@use)
  {
    my ($N) = &F ('./module-N', $use, 1);
    my @U = &F ('.//use-N/N/n', $use, 1);
    for my $U (@U)
      {
        $use{$N}{$U}++;
      }
  }

my ($len) = sort { $b <=> $a } map { length ($_) } keys (%use);



my $t1 = $d->textContent;

'FileHandle'->new (">$F90.new")->print ($t1) if ($t0 ne $t1);

    
