#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

# EN-decl><EN-N><N><n>ZFGELS</n></N></EN-N><array-spec>(<shape-spec-LT><shape-spec><upper-bound
# T-decl-stmt><_T-spec_><derived-T-spec>TYPE(<T-N><N><n>GEOMETRY</n></N></T-N>)</derived-T-spec></_T-spec_>    ,<attribute><attribute-N>INTENT


my @F90 = qw (
src/local/arpifs/phys_dmn/aplpar.F90
src/local/arpifs/phys_dmn/apl_arome.F90
src/local/arpifs/phys_dmn/mf_phys.F90
src/local/arpifs/adiab/cpg.F90
src/local/arpifs/adiab/cpg_dyn.F90
src/local/arpifs/adiab/cpg_end.F90
src/local/arpifs/adiab/cpg_dia.F90
src/local/arpifs/adiab/cpg_gp.F90
);

for my $F90 (@F90)
  {
    print "$F90\n";

    my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);
    
    # KLON, DGEOMETRY%YRDIM%NPROMM, YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIM%NPROMNH
    
    my @v = &F ('//T-decl-stmt[not (.//attribute-N[string (.)="INTENT"])]'
              . '//EN-decl/array-spec/shape-spec-LT'
              . '[starts-with(string(shape-spec),"YDGEOMETRY%YRDIM%NPROM") or '
              . 'starts-with(string(shape-spec),"KLON") or '
              . 'starts-with(string(shape-spec),"KFDIA")]',
              $d);
    
    
    my %sz;
    
    for my $v (@v)
      {
        my @ss = &F ('./shape-spec', $v);
        shift (@ss);
        my $sz = join (' | ', map { $_->textContent } @ss);
        $sz{$sz}++;
      }
    
    print &Dumper (\%sz);
  }
