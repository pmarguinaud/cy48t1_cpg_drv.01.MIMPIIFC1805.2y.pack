#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FileHandle;
use FindBin qw ($Bin);
use lib $Bin;
use Fxtran;
use File::Basename;
use Data::Dumper;

my $NPROMA = 'YDGEOMETRY%YRDIM%NPROMA';
my $NFLEVG = 'YDGEOMETRY%YRDIMV%NFLEVG';

my $AS_void                 = "$NPROMA";
my $AS_0_nflevg             = "$NPROMA,0:$NFLEVG";
my $AS_1_nflevg             = "$NPROMA,$NFLEVG";
my $AS_1_nflevg_2           = "$NPROMA,$NFLEVG,2";


my %dims =
(
  "GDW"            =>  $AS_1_nflevg,        
  "GWHT"           =>  $AS_0_nflevg,        
  "OROGLL"         =>  $AS_void,            
  "OROGMM"         =>  $AS_void,            
  "OROGLM"         =>  $AS_void,            
  "DELNHPRE"       =>  $AS_1_nflevg,        
  "DPHYCTY"        =>  $AS_1_nflevg,        
  "DVER"           =>  $AS_1_nflevg,        
  "DVERL"          =>  $AS_1_nflevg,        
  "DVERM"          =>  $AS_1_nflevg,        
  "DVERW"          =>  $AS_1_nflevg,        
  "EIQCHAF"        =>  $AS_1_nflevg,        
  "EQCHAF"         =>  $AS_1_nflevg,        
  "EQCHAH"         =>  $AS_0_nflevg,        
  "GPHL"           =>  $AS_0_nflevg,        
  "GPHM"           =>  $AS_0_nflevg,        
  "GWHL"           =>  $AS_0_nflevg,        
  "GWHM"           =>  $AS_0_nflevg,        
  "KAPH"           =>  $AS_0_nflevg,        
  "LNNHPREFL"      =>  $AS_1_nflevg,        
  "LNNHPREFM"      =>  $AS_1_nflevg,        
  "NHPPI"          =>  $AS_1_nflevg,        
  "NHPREL"         =>  $AS_1_nflevg,        
  "NHPREM"         =>  $AS_1_nflevg,        
  "NHXD_T"         =>  $AS_1_nflevg,        
  "NHXS"           =>  $AS_1_nflevg,        
  "NHXS_T"         =>  $AS_1_nflevg,        
  "NHXT"           =>  $AS_1_nflevg,        
  "PDEP"           =>  $AS_1_nflevg,        
  "PDEPS"          =>  $AS_void,            
  "PSGRTL"         =>  $AS_1_nflevg,        
  "PSGRTM"         =>  $AS_1_nflevg,        
  "R0T"            =>  $AS_1_nflevg,        
  "R9T"            =>  $AS_1_nflevg,        
  "RDT"            =>  $AS_1_nflevg,        
  "RDTL"           =>  $AS_1_nflevg,        
  "RDTM"           =>  $AS_1_nflevg,        
  "RNHPPI"         =>  $AS_1_nflevg,        
  "RPREF"          =>  $AS_1_nflevg,        
  "RRED"           =>  $AS_1_nflevg,        
  "RT"             =>  $AS_1_nflevg,        
  "RTR"            =>  $AS_1_nflevg,        
  "SGRTL"          =>  $AS_1_nflevg,        
  "SGRTM"          =>  $AS_1_nflevg,        
  "SGRTSL"         =>  $AS_void,            
  "SGRTSM"         =>  $AS_void,            
  "SVDINCR13"      =>  $AS_1_nflevg,        
  "TNDGWF_LAP"     =>  $AS_1_nflevg,        
  "TNDGWF_OTH"     =>  $AS_1_nflevg,        
  "TNDGWH_LAP"     =>  $AS_0_nflevg,        
  "TNDGWH_OTH"     =>  $AS_0_nflevg,        
  "TNDUS"          =>  $AS_void,            
  "TNDVS"          =>  $AS_void,            
  "US"             =>  $AS_void,            
  "US_L"           =>  $AS_void,            
  "US_M"           =>  $AS_void,            
  "VS"             =>  $AS_void,            
  "VS_L"           =>  $AS_void,            
  "VS_M"           =>  $AS_void,            
  "Z3DIVG"         =>  $AS_1_nflevg,        
  "WH2F"           =>  $AS_1_nflevg_2,      
  "ZERO"           =>  $AS_void,            
  "ONE"            =>  $AS_1_nflevg,        
  "UVH_UH"         =>  $AS_0_nflevg,  
  "UVH_VH"         =>  $AS_0_nflevg,  
  "UVH_WWI"        =>  $AS_0_nflevg,  
  "XYBDER_LNPRL"   =>  $AS_1_nflevg,  
  "XYBDER_LNPRM"   =>  $AS_1_nflevg,  
  "XYBDER_ALPHL"   =>  $AS_1_nflevg,  
  "XYBDER_ALPHM"   =>  $AS_1_nflevg,  
  "XYBDER_ALPHPLL" =>  $AS_1_nflevg,  
  "XYBDER_ALPHPLM" =>  $AS_1_nflevg,  
  "XYBDER_COEFD"   =>  $AS_1_nflevg,  
  "XYBDER_COEFA"   =>  $AS_1_nflevg,  
  "XYBDER_COEFAPL" =>  $AS_1_nflevg,  
);


my %v;
my %f;

for my $f (qw (src/local/arpifs/adiab/cpg_gp.F90 
               src/local/arpifs/adiab/cpg_gp_hyd.F90
               src/local/arpifs/adiab/cpg_gp_nhee.F90
               src/local/arpifs/adiab/cpg_gp_nhqe.F90))
  {
    my $g = &basename ($f);

    my $d = &Fxtran::fxtran (location => $f, fopts => [qw (-line-length 300)]);
    
    my @expr = &F ('.//named-E[string(N)="YDTMP"]', $d);
    
    my %vv;

    for my $expr (@expr)
      {
        my ($n) = &F ('./N/n/text()', $expr);
        my @r = &F ('./R-LT/component-R', $expr);
        next unless (@r);
    
        my @n;
        for my $r (@r)
          {
            push @n, $r->textContent;
            $r->unbindNode;
          }
    
        for (@n)
          {
            s/^%//o;
          }
    
        push @n, shift (@n);
    
        my $N = join ('_', 'Z', @n);
    
        $n->setData ($N);

        $v{$g}{$N}++;
        $f{$N}{$g}++;

        $vv{$N}++;
      }
    
    my ($drhook) = &F ('.//if-stmt[./action-stmt/call-stmt[string(procedure-designator)="DR_HOOK"]]', $d);
    for my $v (sort keys (%vv))
      {
        (my $x = $v) =~ s/(?:^Z_|_T\d$)//go;
        my $dim = $dims{$x};
        die $x unless ($dim);
        my $stmt = &Fxtran::fxtran (statement => "REAL (KIND=JPRB) :: $v ($dim)");
        $drhook->parentNode->insertBefore ($stmt, $drhook);
        $drhook->parentNode->insertBefore (&t ("\n"), $drhook);
      }

    $drhook->parentNode->insertBefore (&t ("\n"), $drhook) for (1 .. 2);

    'FileHandle'->new (">$g")->print ($d->textContent);
  }

print &Dumper (\%v);

for my $N (keys (%{ $v{'cpg_gp.F90'} }))
  {
    print "$N:";
    for my $f (qw (cpg_gp_hyd.F90 cpg_gp_nhee.F90 cpg_gp_nhqe.F90))
      {
        print " $f" if ($v{$f}{$N});
      }
    print "\n";
  }


__END__
