#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my $F90 = 'src/local/arpifs/phys_dmn/mf_phys.F90';

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);

my %aplpar    = map { ($_, 1) } &F ('//call-stmt[string(procedure-designator)="APLPAR"]/arg-spec/arg/named-E/N', $d, 1);
my %apl_arome = map { ($_, 1) } &F ('//call-stmt[string(procedure-designator)="APL_AROME"]/arg-spec/arg/named-E/N', $d, 1);

my @a = sort grep { m/^Z/o } grep { $aplpar{$_} } keys (%apl_arome);
my @o = sort grep { m/^Z/o } keys (%{ {%aplpar, %apl_arome} });

print &Dumper (\@a);
print &Dumper (\@o);


my @decl = &F ('//EN-decl[starts-with(string (EN-N),"Z")][array-spec]/EN-N', $d, 1);

#print &Dumper (\@decl);
