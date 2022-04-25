#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my ($F90) = @ARGV;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 800)]);

for my $ms (['YOMCLI', 'CLI'])
  {
    my ($mod, $str) = @$ms;
    my ($rlt) = &F ('.//use-stmt[string(module-N)="?"]/rename-LT', $mod, $d);

    next unless ($rlt);

    my @cst = &F ('.//use-stmt[string(module-N)="?"]//use-N', $mod, $d, 1);
    
    for my $cst (@cst)
      {
        my @expr = &F ('.//named-E[string(N)="?"]/N/n/text()', $cst, $d);
        for my $expr (@expr)
          {
            $expr->setData ("YR$str%$cst");
          }
      }
    
    
    for ($rlt->childNodes)
      {
        $_->unbindNode ();
      }
    
    $rlt->appendChild (&n ("<rename><use-N><N><n>YR$str</n></N></use-N></rename>"));
    
  }    

    
'FileHandle'->new (">$F90.new")->print ($d->textContent);
    
    
