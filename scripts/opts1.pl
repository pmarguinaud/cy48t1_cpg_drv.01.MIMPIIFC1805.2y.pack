#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;
use Associate;

my $F90 = shift;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 500)]);


my @n = qw (LFLASH LXCLP LXTGST LXXGST ND4SYS NDLNPR NFNUDG NINDAT NPDVAR NVDVAR L3DTURB  LAPRXPK LAROME 
            LCALLSFX LCORWAT LELAM LFLASH LGRADSP LGWADV LNHDYN LNHEE LNHEE_REFINE_GRP LNHEE_REFINE_PREH_BBC
            LNHEE_SVDLAPL_FIRST LNHQE LNHX LNHXDER LNUDG LPC_CHEAP LPC_FULL LRDBBC L_RDRY_VD LRPLANE LRUBC 
            LSFORC LSFORCS  LSLHD LSOMEGA_FRC LSPRT LSPS_FRC LSW_FRC LTWOTL LVEREGINT LVERTFE NINDAT RPLDARE 
            RPLRG LVFE_ECMWF LVFE_GW LVFE_LAPL_BC LXCLP LXTGST LXXGST RHYDR0 RPLDARE RPLRG TOPPRES XPNUDG);
@n = qw (NSTEP LSLAG);

my %n = map { ($_, 1) } @n;

for my $expr (&F ('.//named-E', $d))
  {
    my ($t) = &F ('./N/n/text()', $expr);
    my $n = $t->textContent;
    next unless ($n{$n});
    $t->setData ("YDCPG_OPTS%$n");
  }

'FileHandle'->new (">$F90.new")->print ($d->textContent);


