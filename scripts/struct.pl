#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use List::MoreUtils qw (uniq);
use lib $Bin;
use Fxtran;
use GraphViz2;


my $name = shift;

my @F90 = @ARGV;

my @view = do { my $fh = 'FileHandle'->new ("<.gmkview"); my @v = <$fh>; chomp for (@v); @v };

for my $F90 (@F90)
  {
    $F90 =~ s,^src\/\w+/,,o;
    for my $view (@view)
      {
        if (-f "src/$view/$F90")
          {
            $F90 = "src/$view/$F90";
            last;
          }
      }
  }

@F90 = &uniq (@F90);

my $code = '';

for my $F90 (@F90)
  {
    $code .= do { my $fh = 'FileHandle'->new ("<$F90"); local $/ = undef; <$fh> };
  }

my %type;

my $d = &Fxtran::fxtran (string => $code, fopts => [qw (-line-length 400)]);

my @type = &F ('.//T-construct', $d);

for my $type (@type)
  {
    my ($N) = &F ('./T-stmt/T-N', $type, 1);

    my @c = &F ('./component-decl-stmt', $type);
    my @a = &F ('./component-decl-stmt[string(attribute)="POINTER" or string(attribute)="ALLOCATABLE"]', $type);
    my @s = &uniq (&F ('./component-decl-stmt/_T-spec_/derived-T-spec/T-N', $type, 1));

    $type{$N} = [$N, scalar (@c), scalar (@a), \@s];

  }

my %seen;

sub walk
{
  my ($t, $g) = @_;

  return if ($seen{$t->[0]}++);

  $g->add_node (name => $t->[0], label => "$t->[0]\n$t->[1]", shape => 'box', 
                fillcolor => ($t->[2] ? 'red' : 'white'), style => 'filled');

  for my $n (@{ $t->[3] })
    {
      my $t = $type{$n} || [$n, 0, 0, []];
      &walk ($t, $g);
    }

  for my $tt (@{ $t->[3] })
    {
      $g->add_edge (from => $t->[0], to => $tt);
    }

}


my $g = 'GraphViz2'->new (graph => {rankdir => 'TB', ordering => 'out'}, global => {rank => 'source'});

&walk ($type{$name}, $g);

$g->run (format => 'svg', output_file => "$name.svg");


