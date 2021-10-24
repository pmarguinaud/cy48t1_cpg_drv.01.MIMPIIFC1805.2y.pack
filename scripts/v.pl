#!/usr/bin/perl -w
#
use strict;

use FileHandle;


my %x = (
  'PGMVS(JL,YT0%MSP)'      , 'Z0MSP(JL)',
  'PGMVS(JL,YT0%MSPL)'     , 'Z0MSPL(JL)',
  'PGMVS(JL,YT0%MSPM)'     , 'Z0MSPM(JL)',
  'PGMVS(JL,YT9%MSP)'      , 'Z9MSP(JL)',
  'PGMVS(JL,YT9%MSPL)'     , 'Z9MSPL(JL)',
  'PGMVS(JL,YT9%MSPM)'     , 'Z9MSPM(JL)',
  'PGMV(JL,JK,YT0%MDIV)'   , 'Z0MDIV(JL,JK)',
  'PGMV(JL,JK,YT0%MNHX)'   , 'Z0MNHX(JL,JK)',
  'PGMV(JL,JK,YT0%MSPD)'   , 'Z0MSPD(JL,JK)',
  'PGMV(JL,JK,YT0%MSPDL)'  , 'Z0MSPDL(JL,JK)',
  'PGMV(JL,JK,YT0%MSPDM)'  , 'Z0MSPDM(JL,JK)',
  'PGMV(JL,JK,YT0%MSVD)'   , 'Z0MSVD(JL,JK)',
  'PGMV(JL,JK,YT0%MSVDL)'  , 'Z0MSVDL(JL,JK)',
  'PGMV(JL,JK,YT0%MSVDM)'  , 'Z0MSVDM(JL,JK)',
  'PGMV(JL,JK,YT0%MT)'     , 'Z0MT(JL,JK)',
  'PGMV(JL,JK,YT0%MTL)'    , 'Z0MTL(JL,JK)',
  'PGMV(JL,JK,YT0%MTM)'    , 'Z0MTM(JL,JK)',
  'PGMV(JL,JK,YT0%MU)'     , 'Z0MU(JL,JK)',
  'PGMV(JL,JK,YT0%MV)'     , 'Z0MV(JL,JK)',
  'PGMV(JL,JK,YT9%MDIV)'   , 'Z9MDIV(JL,JK)',
  'PGMV(JL,JK,YT9%MNHX)'   , 'Z9MNHX(JL,JK)',
  'PGMV(JL,JK,YT9%MSPD)'   , 'Z9MSPD(JL,JK)',
  'PGMV(JL,JK,YT9%MSPDL)'  , 'Z9MSPDL(JL,JK)',
  'PGMV(JL,JK,YT9%MSPDM)'  , 'Z9MSPDM(JL,JK)',
  'PGMV(JL,JK,YT9%MSVD)'   , 'Z9MSVD(JL,JK)',
  'PGMV(JL,JK,YT9%MSVDL)'  , 'Z9MSVDL(JL,JK)',
  'PGMV(JL,JK,YT9%MSVDM)'  , 'Z9MSVDM(JL,JK)',
  'PGMV(JL,JK,YT9%MT)'     , 'Z9MT(JL,JK)',
  'PGMV(JL,JK,YT9%MTL)'    , 'Z9MTL(JL,JK)',
  'PGMV(JL,JK,YT9%MTM)'    , 'Z9MTM(JL,JK)',
  'PGMV(JL,JK,YT9%MU)'     , 'Z9MU(JL,JK)',
  'PGMV(JL,JK,YT9%MV)'     , 'Z9MV(JL,JK)',
);


my $f = shift;

my $code = do { local $/ = undef; my $fh = 'FileHandle'->new ("<$f"); <$fh> };

while (my ($k, $v) = each (%x))
  {
    $k = quotemeta ($k);
    $code =~ s/$k/$v/gms;
  }

'FileHandle'->new (">$f.new")->print ($code);
