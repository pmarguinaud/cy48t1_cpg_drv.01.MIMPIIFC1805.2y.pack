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

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 200)]);

my @C = &F ('//C[string (.)="!* outline"]', $d);

for my $C (@C)
  {
    my @node;
    for (my $node = $C; ; $node = $node->nextSibling)
      {
        $node or die;
        push @node, $node;
        if (($node->nodeName eq 'C') && ($node->textContent eq '!* end outline'))
          {
            last;
          }
      }
    my %N;
    for my $node (@node)
      {
        my @N = &F ('.//named-E/N/n/text()', $node);
        for (@N)
          {
            $N{$_->textContent} = 1;
          }
      }

    my @N = sort keys (%N);

    my %N2M;
    
    my %S = qw (Z P LL LD I K J K N K L LD YL YD);
    my @S = qw (Z LL I J N L YL);

    my %M;

    for my $N (@N)
      {
        for my $v (values (%S))
          {
            if (index ($N, $v) == 0)
              {
                $N2M{$N} = $N;
                $M{$N} = 1;
                last;
              }
          }
      }
    
    for my $N (@N)
      {
        next if ($N2M{$N});
        for my $k (@S)
          {
            if (index ($N, $k) == 0)
              {
                my $v = $S{$k};
                (my $M = $N) =~ s/^$k/$v/;
                while ($M{$M})
                  {
                    $M .= '_';
                  }
                $N2M{$N} = $M;
                last;
              }
          }
      }

    print &Dumper (\%N2M);
    
  }


'FileHandle'->new (">$F90.new")->print ($d->textContent);






