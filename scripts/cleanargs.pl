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

while (1)
  {
    my $count = 0;

    my @arg = &F ('.//dummy-arg-LT/arg-N', $d);

    for my $arg (@arg)
      {
        my $N = $arg->textContent;
        my @expr = &F ('.//named-E[string(N)="?"]', $N, $d);
        next if (@expr);
        &Fxtran::removeListElement ($arg);

        my ($en_decl) = &F ('.//T-decl-stmt//EN-decl[string(EN-N)="?"]', $N, $d);
        my $decl = &Fxtran::stmt ($en_decl);
        if (&Fxtran::removeListElement ($en_decl))
          {
            $decl->unbindNode ();
          }
        $count++;
      }
  
    last unless ($count);
  }



'FileHandle'->new (">$F90.new")->print ($d->textContent);
