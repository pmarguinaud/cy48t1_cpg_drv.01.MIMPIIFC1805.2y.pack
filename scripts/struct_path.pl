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
my $cle = shift;

my %cle;

{
  my $fh = 'FileHandle'->new ("<$cle");
  while (my $line = <$fh>)
    {
      next if ($line =~ m/^\s*#/o);
      for ($line)
        {
          chomp;
          s/(?:^\s*|\s*$)//go;
        }
      my ($k, $v) = split (m/\s*=\s*/o, $line);
      $cle{$k} = $v;
    }
  $fh->close ();
}

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

    my @n = &F ('./component-decl-stmt//EN-N', $type);

    my %n;

    for my $n (@n)
      {
        my ($stmt) = &Fxtran::stmt ($n);
        my ($t) = &F ('./_T-spec_/derived-T-spec/T-N', $stmt, 1);
        $n{uc ($n->textContent)} = $t ? uc ($t) : $t;
      }

    $type{$N} = \%n;

  }


for my $v (values (%type))
  {
    for my $t (values (%$v))
      {
        $t = $t ? $type{$t} : $t;
      }
  }

my %seen;

sub walk
{
  my ($type, $path, $cle) = @_;

  for my $k (sort keys (%{ $type }))
    {
next if ($k eq 'PREVIOUS');
      if ($type->{$k})
        {
          &walk ($type->{$k}, $path . '%' . $k, $cle);
        }
      else
        {
          if (exists $cle->{$k})
            {
              if ($cle->{$k} =~ m/^(?:TRUE|FALSE)$/o)
                {
                  printf ("IF (%-40s %-8s %-20s) CALL ABOR1('SUAPL_ARPEGE: $k MISMATCH')\n", "$path%$k", '.NEQV.', '.' . $cle->{$k} . '.');
                }
              else
                {
                  printf ("IF (%-40s %-8s %-20s) CALL ABOR1('SUAPL_ARPEGE: $k MISMATCH')\n", "$path%$k", '/=', $cle->{$k});
                }
              $seen{$k}++;
            }
        }
    }
}

&walk ($type{$name}, 'YDMODEL', \%cle);

print "\n";

for my $k (sort keys (%cle))
  {
    next if ($seen{$k});
    if ($cle{$k} =~ m/^(?:TRUE|FALSE)$/o)
      {
        printf ("IF (%-40s %-8s %-20s) CALL ABOR1('SUAPL_ARPEGE: $k MISMATCH')\n", $k, '.NEQV.', '.' . $cle{$k} . '.');
      }
    else
      {
        printf ("IF (%-40s %-8s %-20s) CALL ABOR1('SUAPL_ARPEGE: $k MISMATCH')\n", $k, '/=', $cle{$k});
      }
  }  



