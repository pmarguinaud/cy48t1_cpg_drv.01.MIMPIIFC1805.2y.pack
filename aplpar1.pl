#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my $mf_phys = &Fxtran::fxtran (location => "src/local/arpifs/phys_dmn/mf_phys.F90", 
                               fopts => [qw (-line-length 300)]);

my $r = 'APL_AROME';


 print "\n" x 3, $r, "\n" x 3;

 my $func  = &Fxtran::fxtran (location => "src/local/arpifs/phys_dmn/" . lc ($r) . ".F90", 
                                fopts => [qw (-line-length 300)]);
 
 my @argd = &F ('//dummy-arg-LT/arg-N', $func, 1);
 
 my @call = &F ('//call-stmt[string (procedure-designator)="?"]/arg-spec', $r, $mf_phys);
 
 my %pointer;
 
 for my $pointer (&F ('//pointer-a-stmt', $mf_phys))
   {
     my ($E1) = &F ('./E-1/ANY-E', $pointer);
     my ($E2) = &F ('./E-2/ANY-E', $pointer);
     $pointer{$E1->textContent} = $E2;
   }
 
 my %args;
 
 for my $i (0 .. $#call)
   {
     my @arga = map { $_->textContent } &F ('./arg', $call[$i]);
     for my $j (0 .. $#argd)
       {
         if (exists $args{$argd[$j]})
           {
             if (defined ($args{$argd[$j]}) && ($args{$argd[$j]} ne $arga[$j]))
               {
                 $args{$argd[$j]} = undef;
               }
           }
         else
           {
             $args{$argd[$j]} = $arga[$j];
           }
       }
   }
 
     
 my %ok;
 my $count = 80;
 for my $argd (@argd)
   {
     next unless (my $arg = $args{$argd});
     next unless (my $E2 = $pointer{$arg});
     $ok{$argd} = 1;
     last unless ($count--);
   }

 $ok{PEZDIAG} = 0;

 for my $argd (@argd)
   {
     next unless (my $arg = $args{$argd});
     next unless (my $E2 = $pointer{$arg});
     unless ($ok{$argd})
       {
         print "Skip $argd, $arg, ", $E2->textContent, "\n";
         next;
       }

     print $argd, " ", $arg, " ", $E2->textContent, "\n";
     

     my @expr = &F ('//named-E[string (N)="?"]', $argd, $func);
     for my $expr (@expr)
       {
         my @r = &F ('./R-LT/ANY-R', $expr);
         my $e2 = $E2->cloneNode (1);
         my ($rlt2) = &F ('./R-LT', $e2);
         for my $r (@r)
           {
             $rlt2->appendChild ($r);
           }
         $expr->replaceNode ($e2);
       }

     my ($darg) = &F ('//dummy-arg-LT/arg-N[string (.)="?"]', $argd, $func);
     $darg->nextSibling && $darg->nextSibling->unbindNode ();
     $darg->unbindNode ();

     my ($EN) = &F ('//EN-decl[string (EN-N)="?"]', $argd, $func);
     my $decl = &Fxtran::stmt ($EN);
     $decl->unbindNode ();


      for my $call (@call)
        {
          my ($arga) = &F ('./arg[string (.)="?"]', $arg, $call);
          $arga->nextSibling && $arga->nextSibling->unbindNode ();
          $arga->unbindNode ();
        }

   }



 
'FileHandle'->new (">src/local/arpifs/phys_dmn/" . lc ($r) . ".F90.new")->print ($func->textContent);
'FileHandle'->new (">src/local/arpifs/phys_dmn/mf_phys.F90.new")->print ($mf_phys->textContent);

