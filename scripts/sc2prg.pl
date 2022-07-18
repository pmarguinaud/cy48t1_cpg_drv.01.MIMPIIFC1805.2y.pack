#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

for my $F90 (@ARGV)
  {
    my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 200)]);
    
    my @call = &F ('//call-stmt[string (procedure-designator)="SC2PRG"]', $d);
    
    for my $call (@call)
      {
        my ($argspec) = &F ('./arg-spec', $call);
        my @arg = &F ('./arg/ANY-E', $argspec);
        next unless (scalar (@arg) == 3);
    
       my ($p) = &F ('./N', $arg[2], 1);
        my @ss = &F ('.//T-decl-stmt//EN-decl[string(EN-N)="?"]/array-spec/shape-spec-LT/shape-spec', $p, $d);
        my $nd = scalar (@ss);

        if ($arg[1]->textContent =~ m/^(PGMV\w*)/o)
          {
            $call->unbindNode ();
            my $pgmv = $1;
            my @expr = &F ('.//named-E[string(N)="?"]', $p, $d);

            for my $expr (@expr)
              {
                my ($n) = &F ('./N/n/text()', $expr);
                $n->setData ($pgmv);

                my ($rlt) = &F ('./R-LT', $expr);
                unless ($rlt)
                  {
                    $rlt = &n ('<R-LT><array-R>(<section-subscript-LT>' 
                               . join (',', ('<section-subscript>:</section-subscript>') x $nd)
                               . '</section-subscript-LT>)</array-R></R-LT>');
                    $expr->appendChild ($rlt);
                  }
                my ($r) = &F ('./ANY-R', $rlt);
                my ($sslt) = &F ('./ANY-LT', $r);
                $sslt->appendChild (&t (', '));
                if ($r->nodeName eq 'array-R')
                  {
                    $sslt->appendChild (&n ("<section-subscript><lower-bound>$arg[0]</lower-bound></section-subscript>"));
                  }
                elsif ($r->nodeName eq 'parens-R')
                  {
                    $sslt->appendChild (&n ("<element>$arg[0]</element>"));
                  }
              }
            
          }
        elsif ($arg[1]->textContent =~ m/^PB1$/o)
          {
            $call->unbindNode ();
            my @expr = &F ('.//named-E[string(N)="?"]', $p, $d);

            for my $expr (@expr)
              {
                print $expr->textContent, "\n";
                my ($n) = &F ('./N/n/text()', $expr);
                $n->setData ('PB1');

                my ($rlt) = &F ('./R-LT', $expr);
                unless ($rlt)
                  {
                    $rlt = &n ('<R-LT><array-R>(<section-subscript-LT>' 
                               . join (',', ('<section-subscript>:</section-subscript>') x (1 + $nd))
                               . '</section-subscript-LT>)</array-R></R-LT>');
                    $expr->appendChild ($rlt);
                  }
                my ($r) = &F ('./ANY-R', $rlt);
                my ($sslt) = &F ('./ANY-LT', $r);
                if ($r->nodeName eq 'array-R')
                  {
                    $expr->replaceNode (&e ("PB1 (:," . $arg[0]->textContent .  ")"));
                  }
                elsif ($r->nodeName eq 'parens-R')
                  {
                    die;
                  }
                
              }
          }
        elsif ($arg[1]->textContent =~ m/^PB2$/o)
          {
            my $cc = 0;
            my @expr = &F ('.//named-E[string(N)="?"]', $p, $d);

            for my $expr (@expr)
              {
                print $expr->textContent, "\n";
                my ($n) = &F ('./N/n/text()', $expr);
                $n->setData ('PB2');

                my ($rlt) = &F ('./R-LT', $expr);
                unless ($rlt)
                  {
                    $rlt = &n ('<R-LT><array-R>(<section-subscript-LT>' 
                               . join (',', ('<section-subscript>:</section-subscript>') x (1 + $nd))
                               . '</section-subscript-LT>)</array-R></R-LT>');
                    $expr->appendChild ($rlt);
                  }
                my ($r) = &F ('./ANY-R', $rlt);
                my ($sslt) = &F ('./ANY-LT', $r);
                if ($r->nodeName eq 'array-R')
                  {
                    $expr->replaceNode (&e ("PB2 (:," . $arg[0]->textContent .  ")"));
                  }
                elsif ($r->nodeName eq 'parens-R')
                  {
                    $cc++;
                    die;
                  }
                
              }

            $call->unbindNode () unless ($cc);
          }
    
      }
        
    
    'FileHandle'->new (">$F90.new")->print ($d->textContent);
    

  }
