#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my @F90 = qw (src/local/arpifs/adiab/cpg_gp_hyd.F90
              src/local/arpifs/adiab/cpg_gp_nhqe.F90
              src/local/arpifs/adiab/cpg_gp_nhee.F90);

for my $F90 (@F90)
  {

    my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);
 
    for my $s (qw (PCTY0))
      {    
    
        my ($k, $t) = ($s =~ m/P(\w+)(\d)$/o);

        for my $expr (&F ('//named-E[string (N)="?"]', $s, $d))
          {
        #   print $expr->textContent, "\n";
            my ($N) = &F ('./N', $expr);
            my ($r) = &F ('./R-LT/ANY-R', $expr);
        
            next unless ($r);

            my $p;
            if ($r->nodeName eq 'parens-R')
              {
                ($p) = &F ('./element-LT/element[3]', $r);
              }
            elsif ($r->nodeName eq 'array-R')
              {
                ($p) = &F ('./section-subscript-LT/section-subscript[3]', $r);
              }
        
            $p->previousSibling->unbindNode;
            $p->unbindNode;
            $p = $p->textContent;
            $p =~ s/^YY\w+%M_//o;    
        
            $N->replaceNode (&t ("YDCPG_DYN$t%$k%$p"));
        
            if ($r->textContent eq '(1,1)')
              {
                $r->replaceNode (&t ('(:,1:)'));
              }
        
            if ($r->textContent eq '(1,0)')
              {
                $r->unbindNode ();
              }
        
            if ($r->textContent eq '(:,:)')
              {
                $r->unbindNode ();
              }
        
          }
      }
    
    
    'FileHandle'->new (">$F90.new")->print ($d->textContent);

  }
