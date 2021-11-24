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

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 200)]);

die unless ($d);

my @w = qw (
ZKAP0H
ZEQCHA0F
ZEQCHA0H
ZEIQCHA0F
ZRPRE0F
ZRT0
ZR0T9
ZR9T9
ZDVER0
ZDVER0L
ZDVER0M
ZDVER9
ZUS0
ZVS0
ZUS0_L
ZVS0_L
ZUS0_M
ZVS0_M
ZUS9
ZVS9
ZSGRTSL
ZSGRTSM
ZPSGRTL
ZPSGRTM
ZSGRTL
ZSGRTM
ZTNDUS
ZTNDVS
ZGPHL
ZGPHM
ZGWHL
ZGWHM
ZNHXS_T0
ZNHXD_T0
ZNHXS_T9
ZNHXD_T9
ZONE
);



for my $w (@w)
  {
    (my $n = $w) =~ s/^Z//o;

    my $x = reverse ($n);
    my ($t) = ($x =~ m/([09])/o); $t = defined ($t) ? $t : '0';
    $x =~ s/([09])//o;

    $n = reverse ($x);

    my @tt = &F ('.//named-E[string(N)="?"]/N/n/text()', $w, $d);

    for my $tt (@tt)
      {
        $tt->setData ("YDTMP%T$t%$n");
      }

  }

'FileHandle'->new (">$F90.new")->print ($d->textContent);
