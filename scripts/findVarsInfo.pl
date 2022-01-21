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

my $F90 = shift;

my $code = do { my $fh = 'FileHandle'->new ("<$F90"); local $/ = undef; <$fh> };

my $d = &Fxtran::fxtran (string => $code, fopts => [qw (-line-length 400)]);

# Remove pointers

for (&F ('.//T-construct//component-decl-stmt[string(attribute)="POINTER"]', $d))
  {
    $_->unbindNode ();
  }

# Remove scalars (not derived types)

for (&F ('.//T-construct//component-decl-stmt[not(_T-spec_/derived-T-spec)][not(.//array-spec)]', $d))
  {
    $_->unbindNode ();
  }

my %decl;
my %Decl;

for my $T (&F ('.//T-construct', $d))
  {
    my ($N) = &F ('./T-stmt/T-N', $T, 1);
    for my $cd (&F ('./component-decl-stmt', $T))
      {
        my ($n) = &F ('.//EN-N', $cd, 1);
        my ($type) = &F ('./_T-spec_/derived-T-spec/T-N', $cd, 1);
        if ($type)
          {
            if ($type =~ m/^VARIABLE_(\d)D$/o)
              {
                $decl{$N}{$n} = $1;
                $Decl{$N}{$n} = $1;
              }
            else
              {
                $decl{$type} ||= {};
                $decl{$N}{$n} = $decl{$type};
                $Decl{$type} ||= {};
                $Decl{$N}{$n} = $Decl{$type};
              }
          }
      }
  }

sub walk
{
  my ($path, $r, $h) = @_;

  if (ref ($r))
    {
      while (my ($k, $v) = each (%$r))
        {
          &walk ($path . '%' . $k, $v, $h);
        }
    }
  else
    {
      my ($n) = ($path =~ m/%(\w+)$/o);

      for my $suf (qw (T0 T1 T9 PH9 DL DM DL9 DM9))
        {
          if ($r == 2)
            {
              $h->{"$path%$suf"} = "REAL (KIND=JPRB) :: $n (YDCPG_DIM%KLON)";
            }
          elsif ($r == 3)
            {
              $h->{"$path%$suf"} = "REAL (KIND=JPRB) :: $n (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)";
            }
        }

    }
  
}

my %h;

while (my ($k, $v) = each (%decl))
  {
    &walk ($k, $v, \%h);
  }


#%h = (%{ -f 'h.pl' ? do ('./h.pl') : {} }, %h);

local $Data::Dumper::Terse = 1;

'FileHandle'->new ('>h.pl')->print (&Dumper (\%h));

my @D = 
(
  undef,
  undef,
  &Fxtran::fxtran (statement => 'REAL (KIND=JPRB) :: X (YDCPG_DIM%KLON)'),
  &Fxtran::fxtran (statement => 'REAL (KIND=JPRB) :: X (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)'),
);

sub Walk
{
  my ($path, $r, $H) = @_;

  if (ref ($r) eq 'HASH')
    {
      while (my ($k, $v) = each (%$r))
        {
          &Walk ($path . '%' . $k, $v, $H);
        }
    }
  else
    {
      my ($n) = ($path =~ m/%(\w+)$/o);

      for my $suf (qw (T0 T1 T9 PH9 DL DM DL9 DM9))
        {
          if ($r < scalar (@D))
            {
              my $decl = $D[$r]->cloneNode (1);
              my ($N) = &F ('.//EN-N/N/n/text()', $decl);
              $N->replaceNode (&t ($suf));
              $H->{"$path%$suf"} = $decl;
            }
        }

    }
  
}


my %H;

while (my ($k, $v) = each (%Decl))
  {
    &Walk ($k, $v, \%H);
  }

my ($doc, $list);

my $indent = 1;

if (-f 'h.xml')
  {
    $doc = 'XML::LibXML'->load_xml (location => 'h.xml');
    $list = $doc->documentElement ();

    for my $key (sort keys (%H))
      {
        my ($decl) = &F ('./decl[@key="?"]', $key, $list);
        $decl && $decl->unbindNode ();
      }

  }
else
  {
    $doc = 'XML::LibXML::Document'->new ();
    $doc->setDocumentElement ($list = &n ('<list/>'));
    $list->appendChild (&t ("\n" . (' ' x 4))) if ($indent);
  }


for my $key (sort keys (%H))
  {
    my $val = $H{$key};
    my $decl = &n ('<decl/>');

    $list->appendChild ($decl);
    $list->appendChild (&t ("\n" . (' ' x 4))) if ($indent);
    $decl->appendChild ($val->cloneNode (1));
    $decl->setAttribute (key => $key);
  }

$doc->toFile ('h.xml');


