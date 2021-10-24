#!/usr/bin/perl -w

use strict;
use FileHandle;

my $f = shift;

my $code = do { local $/ = undef; my $fh = 'FileHandle'->new ("<$f"); <$fh> };

my %x = (
   'PGMV(JROF,JLEV,YT0%MU)'     , 'Z0U(JROF,JLEV)'              ,
   'PGMV(JROF,JLEV,YT0%MV)'     , 'Z0V(JROF,JLEV)'              ,
   'PGMV(JROF,JLEV,YT0%MDIV)'   , 'Z0DIV(JROF,JLEV)'            ,
   'PGMV(JROF,JLEV,YT0%MTL)'    , 'Z0TL(JROF,JLEV)'             ,
   'PGMV(JROF,JLEV,YT0%MTM)'    , 'Z0TM(JROF,JLEV)'             ,
   'PGMV(JROF,JLEV,YT9%MU)'     , 'Z9U(JROF,JLEV)'              ,
   'PGMV(JROF,JLEV,YT9%MV)'     , 'Z9V(JROF,JLEV)'              ,
   'PGMV(JROF,JLEV,YT0%MUL)'    , 'Z0UL(JROF,JLEV)'             ,
   'PGMV(JROF,JLEV,YT0%MVL)'    , 'Z0VL(JROF,JLEV)'             ,
   'PGMV(JROF,JLEV,YT0%MVOR)'   , 'Z0VOR(JROF,JLEV)'            ,
   'PGMV(JROF,JLEV,YT0%MSPDL)'  , 'Z0SPDL(JROF,JLEV)'           ,
   'PGMV(JROF,JLEV,YT0%MSPDM)'  , 'Z0SPDM(JROF,JLEV)'           ,
   'PGMV(JROF,JLEV,YT0%MSVDL)'  , 'Z0SVDL(JROF,JLEV)'           ,
   'PGMV(JROF,JLEV,YT0%MSVDM)'  , 'Z0SVDM(JROF,JLEV)'           ,
   'PGMV(JROF,JLEV,YT0%MNHXL)'  , 'Z0NHXL(JROF,JLEV)'           ,
   'PGMV(JROF,JLEV,YT0%MNHXM)'  , 'Z0NHXM(JROF,JLEV)'           ,
   'PGMV(JROF,JLEV,YT9%MDIV)'   , 'Z9DIV(JROF,JLEV)'            ,
   'PGMV(JROF,JLEV,YT9%MTL)'    , 'Z9TL(JROF,JLEV)'             ,
   'PGMV(JROF,JLEV,YT9%MTM)'    , 'Z9TM(JROF,JLEV)'             ,
   'PGMV(JROF,JLEV,YT9%MSPDL)'  , 'Z9SPDL(JROF,JLEV)'           ,
   'PGMV(JROF,JLEV,YT9%MSPDM)'  , 'Z9SPDM(JROF,JLEV)'           ,
   'PGMV(JROF,JLEV,YT9%MSVDL)'  , 'Z9SVDL(JROF,JLEV)'           ,
   'PGMV(JROF,JLEV,YT9%MSVDM)'  , 'Z9SVDM(JROF,JLEV)'           ,
   'PGMVS(JROF,YT0%MSPL)'       , 'Z0SPL(JROF)'                 ,
   'PGMVS(JROF,YT0%MSPM)'       , 'Z0SPM(JROF)'                 ,
   'PGMVS(JROF,YT9%MSPL)'       , 'Z9SPL(JROF)'                 ,
   'PGMVS(JROF,YT9%MSPM)'       , 'Z9SPM(JROF)'                 ,
);

while (my ($k, $v) = each (%x))
  {
    $k = quotemeta ($k);
    $code =~ s/$k/$v/gms;
  }

'FileHandle'->new (">$f.new")->print ($code);
