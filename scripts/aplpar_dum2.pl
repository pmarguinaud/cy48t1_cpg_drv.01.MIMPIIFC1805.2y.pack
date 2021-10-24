#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my $F90 = "src/local/arpifs/phys_dmn/aplpar_dum.F90";

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 200)]);
    
my ($sub) = &F ('.//subroutine-stmt[string(subroutine-N)="APLPAR_DUM"]', $d);

my @arg = &F ('./dummy-arg-LT/arg-N', $sub);
my %ind = map { ($arg[$_]->textContent, $_) } (0 .. $#arg);


for (&F ('//C', $d))
  {
    $_->unbindNode ();
  }

my @en = &F ('.//EN-N', $d);

my %en;

for my $en (@en)
  {
    my $stmt = &Fxtran::stmt ($en);
    $en{$en->textContent} = $stmt;
  }

@en = map { $_->textContent } @en;

for my $en (@en)
  {
    die ($en . "\n") unless (defined $ind{$en});
  }

@en = sort { $ind{$a} <=> $ind{$b} } @en;

for (values (%en))
  {
    $_->unbindNode ();
  }

my @stmt = map { $en{$_} } @en;

my ($end) = &F ('//end-subroutine-stmt', $d);


for (@stmt)
  {
    $end->parentNode->insertBefore ($_, $end);
    $end->parentNode->insertBefore (&t ("\n"), $end);
  }


print $d->textContent;


'FileHandle'->new (">$F90")->print ($d->textContent);


