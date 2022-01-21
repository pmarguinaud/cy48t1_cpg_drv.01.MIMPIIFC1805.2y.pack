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


my $code = '';

for my $F90 (@ARGV)
  {
    $code .= do { my $fh = 'FileHandle'->new ("<$F90"); local $/ = undef; <$fh> };
  }

$code =~ s/^!>/  /goms;


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
        $cd->setNodeName ('T-decl-stmt');

        my ($n) = &F ('.//EN-N', $cd, 1);
        my ($ts) = &F ('./_T-spec_', $cd);
        $ts->nextSibling->replaceNode (&t (' :: '));

        my ($type) = &F ('./_T-spec_/derived-T-spec/T-N', $cd, 1);
        if ($type)
          {
            $decl{$type} ||= {};
            $decl{$N}{$n} = $decl{$type};
            $Decl{$type} ||= {};
            $Decl{$N}{$n} = $Decl{$type};
          }
        else
          {
            $Decl{$N}{$n} = $cd;
            for ($decl{$N}{$n} = $cd->textContent)
              {
                s/\bNLEV\b/YDCPG_DIM%KFLEVG/goms;
                s/\bNPROMA\b/YDCPG_DIM%KLON/goms;
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
      %$h = (%$h, $path, $r);
    }
  
}

my %h;

while (my ($k, $v) = each (%decl))
  {
    &walk ($k, $v, \%h);
  }

%h = (%{ -f 'h.pl' ? do ('./h.pl') : {} }, %h);
local $Data::Dumper::Terse = 1;
'FileHandle'->new ('>h.pl')->print (&Dumper (\%h));


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
      %$H = (%$H, $path, $r);
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

