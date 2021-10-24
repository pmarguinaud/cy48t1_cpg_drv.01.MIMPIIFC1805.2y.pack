#!/usr/bin/perl -w
#
use strict;
use FileHandle;

my $f = shift;

my $code = do { local $/ = undef; my $fh = 'FileHandle'->new ("<$f"); <$fh> };

for ($code)
  {
    s/PXYB\(([^,]+,[^,]+),YYTXYB%M_(\w+)\)/P$2($1)/goms;
  }

'FileHandle'->new (">$f.new")->print ($code);
