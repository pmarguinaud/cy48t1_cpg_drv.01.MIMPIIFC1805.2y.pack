#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use lib $Bin;
use Fxtran;
use List::MoreUtils qw (uniq);
use File::Path;


my @view = do { my $fh = 'FileHandle'->new ('<.gmkview'); <$fh> };
chomp for (@view);

my @F90 = <src/*/arpifs/module/*.F90>;


for (@F90)
  {
    s,^src/\w+/,,o;
  }

@F90 = &uniq (sort @F90);

@F90 = map 
         { 
           my $F90 = $_; 
           for my $view (@view) 
             { 
               if (-f "src/$view/$F90") { $F90 = "src/$view/$F90"; last; }
             }
           $F90
         } @F90;


sub name2list
{
  my ($file, $name2list, $name2file) = @_;

  my $doc = &Fxtran::fxtran (location => $file, dir => 'tmp', fopts => [qw (-line-length 800)]);
  die $file unless ($doc);

  my ($pu) = &F ('.//program-unit[./module-stmt]', $doc);

  my @tconst = &F ('./T-construct', $pu);

  for my $tconst (@tconst)
    {
      my ($name) = &F ('.//T-stmt/T-N/N/n/text()', $tconst, 1);
      my @list = &F ('./component-decl-stmt/_T-spec_/derived-T-spec/T-N/N/n/text()', $tconst, 1);
      $name2list->{$name} = \@list;
      $name2file->{$name} = $file;
    }

}

sub walk
{
  my ($name, $level, $name2list, $name2file, $seen) = @_;

  print '  ' x $level, $name, ' ', ($name2file->{$name} || '?'), "\n";
  return if ($seen->{$name}++);

  for (@{ $name2list->{$name} || [] })
    {
      &walk ($_, $level+1, $name2list, $name2file, $seen);
    }

}

my $name2list = {};
my $name2file = {};

for my $F90 (@F90)
  {
    &name2list ($F90, $name2list, $name2file);
  }

my $seen = {};

&walk ('MODEL', 0, $name2list, $name2file, $seen);

