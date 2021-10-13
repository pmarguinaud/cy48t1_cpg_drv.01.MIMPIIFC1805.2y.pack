#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my $F90 = 'src/local/arpifs/adiab/cpg_gp.F90';

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);

my @ap = &F ('//pointer-a-stmt', $d);

for my $ap (@ap)
  {
    my ($E1) = &F ('./E-1/ANY-E', $ap);
    my ($E2) = &F ('./E-2/ANY-E', $ap);
    next if ($E2->textContent =~ m/\(/o);


    my ($N1) = &F ('./N', $E1, 1);

    my @expr = &F ('//named-E/N[string(.)="?"]', $N1, $d);

    next unless (scalar (@expr) == 1);

    print $N1, "\n";

    my ($en_decl) = &F ('//EN-decl[string(EN-N)="?"]', $N1, $d);

    if ($en_decl->nextSibling)
      {
        $en_decl->nextSibling->unbindNode ();
        $en_decl->unbindNode ();
      }
    elsif ($en_decl->previousSibling)
      {
        $en_decl->previousSibling->unbindNode ();
        $en_decl->unbindNode ();
      }
    else
      {
        &Fxtran::stmt ($en_decl)->unbindNode ();
      }

    $ap->unbindNode ();
  }

'FileHandle'->new (">$F90.new")->print ($d->textContent);

