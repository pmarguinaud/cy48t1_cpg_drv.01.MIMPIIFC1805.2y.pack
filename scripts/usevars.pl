#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;
use File::Find;

my %file = %{ do ("./file.pl") };

sub wanted
{
  return unless ((my $f = $File::Find::name) =~ m/\.(?:F|F90)$/o);
  (my $g = $f) =~ s,^src/(\w+)/,,o;
  $file{$g} = $f;
}


if (0)
  {
    my @view = do { my $fh = 'FileHandle'->new ("<.gmkview"); <$fh> };
    chomp for (@view);
    
    for my $view (@view)
      {
        &find ({wanted => \&wanted, no_chdir => 1}, "src/$view/");
      }
    
    {local $Data::Dumper::Terse = 1; 'FileHandle'->new (">file.pl")->print (&Dumper (\%file));}
  }


my %mod2var;


my %skip = map { ($_, 1) } qw (append_l2c.F90 mode_gltools_interp.F90);

for my $file (values (%file))
  {
    next if ($skip{&basename ($file)});
    next unless ($file =~ m/\.F90$/o);
    next unless (my $d = &Fxtran::fxtran (location => $file, fopts => [qw (-line-length 300)], dir => '/tmp'));

    my @pu = &F ('.//program-unit[./module-stmt]', $d);

    for my $pu (@pu)
      {
        my ($mod) = &F ('./module-stmt/module-N', $pu, 1);
        my @var = &F ('./T-decl-stmt[not(./attribute[string(./attribute-N)="PARAMETER"])]//EN-N', $pu, 1);
        for my $var (@var)
          {
            $mod2var{$mod}{$var}++;
          }
      }
  }

{local $Data::Dumper::Terse = 1; 'FileHandle'->new (">mod2var.pl")->print (&Dumper (\%mod2var)); }


