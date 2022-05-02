package Object;

use strict;
use Fxtran;
use Data::Dumper;

{

my $h = do ('./h.pl');
my %decl;

sub getObjectDecl
{
  my $key = shift;

  unless ($decl{$key}) 
    {
      $h->{$key} or die $key;
      ($decl{$key}) = &Fxtran::fxtran (statement => $h->{$key});
    }

  return $decl{$key};
}

my %type;

sub getObjectType
{
  my ($s, $obj) = @_;

  unless ($type{$obj})
    {
      ($type{$obj}) = &F ('./T-N', $s->{ts}, 1);
    }

  return $type{$obj};
}

}

sub getObjectAS
{
  my ($key) = @_;

  my $decl = &getObjectDecl ($key);

  return &asFromDecl ($decl);
}

sub asFromDecl
{
  my $decl = shift;

  my ($as) = &F ('.//EN-decl/array-spec', $decl);

  return $as;
}

sub getFieldFromObjectComponents
{
  my ($obj, @ctl) = @_;

  if ($ctl[-1] =~ m/^(?:T[019]|(?:DM|DL)[019]?)$/o)
    {
      $ctl[-1] = 'F' . $ctl[-1]; 
    }
  elsif ($obj eq 'YDMF_PHYS_SURF') 
    {
      if ($ctl[-1] =~ m/^P(\w+)_T[019]$/o)
        {
          $ctl[-1] =~ s/^P/F_/o;
        }
      else
        {
          $ctl[-1] =~ s/^P/F_/o;
        }
    }
  else
    {
      $ctl[-1] = 'F_' . $ctl[-1]; 
    }

  return &n ("<named-E><N><n>$obj</n></N><R-LT>" 
           . join ('', map { "<component-R>%<ct>$_</ct></component-R>" } @ctl)
           . "</R-LT></named-E>");
}

1;
