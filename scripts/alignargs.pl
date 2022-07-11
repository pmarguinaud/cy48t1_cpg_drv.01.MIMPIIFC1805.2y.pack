#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FileHandle;
use FindBin qw ($Bin);
use lib $Bin;
use Fxtran;
use File::Basename;
use Data::Dumper;
use List::MoreUtils qw (uniq);
use List::Util qw (min max);

my $F90 = shift;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);

my $t0 = $d->textContent;

my @decl = &F ('.//T-decl-stmt[.//attribute-N[string(.)="INTENT"]', $d);

my %len;
my %att;

for my $decl (@decl)
  {
    my ($tspec) = &F ('./_T-spec_', $decl, 1); $tspec =~ s/\s+//go;
    
    $len{type} = &max ($len{type} || 0, length ($tspec));

    my @attr = &F ('.//attribute', $decl);
    for my $attr (@attr)
      {
        my ($N) = &F ('./attribute-N', $attr, 1);
        $attr = $attr->textContent; $attr =~ s/\s+//go;
        $att{$N} = 1;
        $len{$N} = &max (($len{$N} || 0), length ($attr));
      }

  }


my @att = sort keys (%att);

for (values (%len))
  {
    $_++;
  }


for my $decl (@decl)
  {
    my ($tspec) = &F ('./_T-spec_', $decl, 1); $tspec =~ s/\s+//go;
    my ($endlt) = &F ('./EN-decl-LT', $decl, 1);


    my @attr = &F ('.//attribute', $decl);
    my %attr = map { my ($N) = &F ('./attribute-N', $_, 1); ($N, $_->textContent) } @attr;

    for (values (%attr))
      {
        s/\s+//go;
      }
  
    my $code = sprintf ("%-$len{type}s", $tspec);

    for my $att (@att)
      {
        if ($attr{$att})
          {
            $code .= sprintf (",%-$len{$att}s", $attr{$att});
          }
        else
          {
            $code .= ' ' x ($len{$att} + 1);
          }
      }

    $code .= ' :: ' . $endlt;

    my $stmt = &Fxtran::fxtran (statement => $code, fopts => [qw (-line-length 300)]);

    $decl->replaceNode ($stmt);

  }

my $t1 = $d->textContent;

'FileHandle'->new (">$F90.new")->print ($t1) if ($t0 ne $t1);

    
