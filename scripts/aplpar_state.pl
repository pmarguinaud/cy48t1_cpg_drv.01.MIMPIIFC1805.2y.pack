#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my $F90 = "src/local/arpifs/phys_dmn/mf_phys.F90";

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);


my ($call1, $call2) = &F ('//call-stmt[string (procedure-designator)="?"]/arg-spec', 'APL_AROME', $d);

my @arg1 = &F ('./arg', $call1);
my @arg2 = &F ('./arg', $call2);

my $count = 0;

for my $i (0 .. $#arg1)
  {
    my ($arg1, $arg2) = ($arg1[$i], $arg2[$i]);
    my ($a1, $a2) = ($arg1->textContent, $arg2->textContent);
    if ($a1 ne $a2)
      {
        my $f = $a1;
        if ($a1 =~ m/^YDVARS%/o)
          {
            for ($f)
              {
                s/^YDVARS%//o;
                s/%T0$//o;
              }
          }
        elsif ($a1 =~ m/^YDCPG_/o)
          {
            for ($f)
              {
                s/^YD/Y/o;
                s/[09]%/%/o;
              }
          }
        elsif ($a1 =~ m/^ZP/o)
          {
            for ($f)
              {
                s/^Z//o;
                s/[09]$//o;
              }
          }
        elsif ($a1 =~ m/^YDMF_PHYS_SURF/o)
          {
            for ($f)
              {
                s/^YDMF_PHYS_SURF%/Y/o;
                s/_T[09]$//o;
                s/%P/%/o;
              }
          }
        else
          {
            print "$a1 $a2\n";
            next;
          }
        $arg1->replaceNode (&t ("YLAPLPAR%$f"));
        $arg2->replaceNode (&t ("YLAPLPAR%$f"));
#       last if ($count++ > 7);
      }
  }

# 10 KO
#  7 KO
#  6 OK
#  5 OK


'FileHandle'->new (">$F90.new")->print ($d->textContent);
    
    
    
    
    
    
    
    
