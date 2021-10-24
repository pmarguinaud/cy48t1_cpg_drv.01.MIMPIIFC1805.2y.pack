#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my %F90 = 
(
  MF_PHYS   => 'src/local/arpifs/phys_dmn/mf_phys.F90',
  CPG_DIA   => 'src/local/arpifs/adiab/cpg_dia.F90',
  APLPAR    => 'src/local/arpifs/phys_dmn/aplpar.F90',
  APL_AROME => 'src/local/arpifs/phys_dmn/apl_arome.F90',
  CPG       => 'src/local/arpifs/adiab/cpg.F90',
);
my @r = sort keys (%F90);

my %d = map { my $r = $_; ($r, &Fxtran::fxtran (location => $F90{$r}, fopts => [qw (-line-length 300)])) } @r;

my %i;
my %y;

for my $r (@r)
  {
    my $d = $d{$r};
    my @args = &F ('./object/file/program-unit/subroutine-stmt/dummy-arg-LT/arg-N', $d, 1);
    $i{$r} = \@args;
  }



my @x = qw (
  ZDIFCQ ZDIFCQL ZDIFCQN ZDIFCS ZDIFTQ ZDIFTQL ZDIFTQN ZDIFTS ZFCCQL
  ZFCCQN ZFCSQL ZFCSQN ZFCQLNG ZFCQNNG ZFCQNG ZFPLCL ZFPLCN ZFPLCG
  ZFPLCH ZFPLSL ZFPLSN ZFPLSG ZFPLSH ZFPFPSL ZFPFPSN ZFPFPSG ZFPFPCL
  ZFPFPCN ZFPEVPSL ZFPEVPSN ZFPEVPSG ZFPEVPCL ZFPEVPCN ZFPEVPCG ZFRSO ZFRTH
  ZSTRCU ZSTRCV ZSTRDU ZSTRDV ZSTRTU ZSTRTV ZSTRMU ZSTRMV ZDIFCQLC
  ZDIFCQIC ZFIMCC  ZFEDQLC ZFEDQIC ZFEDQRC ZFEDQSC ZFCNEGQLC ZFCNEGQIC ZFCNEGQRC
  ZFCNEGQSC ZFRMH ZFCHOZ ZFDIS  ZFHPSL ZFHPSN ZFHPSG ZFHPCL ZFHPCN
  ZFHSCL ZFHSCN ZFHPCG ZFHSSL ZFHSSN ZFEPFP ZFCMPCQ ZFCMPSN ZFHSSG
  ZFCMPSL ZFRSOC ZFRTHC ZFCHSP ZFCLL ZFCLN ZFCS ZFEVL ZFEVN
  ZFEVV ZFLWSP ZFTR ZFRSODS ZFRSOPS ZFRSDNI ZFRSGNI ZFRTHDS ZFONTE
  ZFGEL ZFGELS ZALB ZGZ0 ZGZ0H ZRUISL ZRUISP ZRUISS ZFRSOPT
  ZFRSOLU ZQCLS ZTCLS ZUCLS ZVCLS ZNUCLS ZNVCLS ZRHCLS ZMRT
  ZCLCH ZCLCM ZCLCL ZCLCC ZCLPH ZVEIN ZDRNSHF ZCAPE ZCTOP
  ZMOCON ZUGST ZVGST ZCT ZTENDU ZTENDV ZFCQRNG ZFCQSNG ZFCQGNG
  ZDIAGH ZVISICLD ZFLASH ZVISIHYD ZMXCLWC ZTPWCLS ZCUCONVCA ZNLCONVCA
);

my %X = (CPG => {map { ($_, $_) } @x});

for my $p (qw (CPG MF_PHYS))
  {
    my @call = &F ('//call-stmt', $d{$p});

    for my $call (@call)
      {
        my ($r) = &F ('./procedure-designator', $call, 1);
        next unless ($F90{$r});
    
        my @args = &F ('./arg-spec/arg', $call, 1);
    
        for my $i (0 .. $#args)
          {
            $y{$r}[$i] = $X{$p}{$args[$i]};
            $X{$r}{$i{$r}[$i]} = $X{$p}{$args[$i]};
          }
        
      }

  }

for my $r (@r)
  { 
    my $d = $d{$r};
    while (my ($x, $y) = each (%{ $X{$r} }))
      {
        next unless ($y);
        (my $z = $y) =~ s/^Z//o;
        my @expr = &F ('//named-E[string(N)="?"]/N', $x, $d);
        for my $expr (@expr)
          {
            my $stmt = &Fxtran::stmt ($expr);
            if ($stmt->nodeName eq 'call-stmt')
              {
                my ($p) = &F ('./procedure-designator', $stmt, 1);
                next if ($d{$r});
              }
            $expr->replaceNode (&t ("YDMF_PHYS%$z"));
          }
      }
    'FileHandle'->new (">$F90{$r}")->print ($d->textContent);
  }





