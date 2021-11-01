#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

my %c;

while (<DATA>)
  { 
    my ($k, $v) = split (m/\s*=\s*/o);
    for ($k, $v)
      {
        s/^\s*//o;
        s/\s*$//o;
      }
    next unless ($v =~ m/^(?:TRUE|FALSE)$/o);
    $c{$k} = $v;
  }

#my $F90 = "src/local/arpifs/phys_dmn/aplpar.F90";
my $F90 = shift;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);

while (my ($k, $v) = each (%c))
  {
    my @expr = &F ('//named-E[string(N)="?"]', $k, $d);
    for my $expr (@expr)
      {
        $expr->replaceNode (&n ("<literal-E>.$v.</literal-E>"));
      } 
  }

sub simplify
{
  my $e1 = shift;
  my $e = $e1->parentNode;
  return unless ($e->nodeName =~ m/-E$/o);
  
  my $nn;

  if (my ($op) = &F ('./op', $e, 1))
    {
      my ($e2) = grep { $_->unique_key != $e1->unique_key } &F ('./ANY-E', $e);
      if ($e1->textContent eq '.TRUE.')
        {
          if ($op eq '.AND.')
            {
              $nn = $e2;
            }
          elsif ($op eq '.OR.')
            {
              $nn = $e1;
            }
          elsif ($op eq '.NOT.')
            {
              $nn = &n ('<literal-E>.FALSE.</literal-E>');
            }
        }
      if ($e1->textContent eq '.FALSE.')
        {
          if ($op eq '.AND.')
            {
              $nn = $e1;
            }
          elsif ($op eq '.OR.')
            {
              $nn = $e2;
            }
          elsif ($op eq '.NOT.')
            {
              $nn = &n ('<literal-E>.TRUE.</literal-E>');
            }
        }
    }
  elsif ($e->nodeName eq 'parens-E')
    {
      $nn = $e1;
    }

  unless ($nn)
     {
#      print $e->textContent, "\n";
       return;
     }

  $e->replaceNode ($nn); 

  if ($nn->textContent =~ m/^\.(?:TRUE|FALSE)\.$/o)
    {
      &simplify ($nn);
    }

}


my @expr = &F ('//literal-E[string(.)=".FALSE." or string(.)=".TRUE."]', $d);

for my $expr (@expr)
  {
    next unless (&Fxtran::stmt ($expr));
    &simplify ($expr);
  }

'FileHandle'->new (">$F90.1.new")->print ($d->textContent);


for my $ce (&F ('.//if-stmt/condition-E[string(.)=".TRUE." or string(.)=".FALSE."]', $d))
  {
    my $stmt = $ce->parentNode;
    if ($stmt->nodeName eq 'if-stmt')
      {
        if ($ce->textContent eq '.TRUE.')
          {
            my ($stmt1) = &F ('./action-stmt/ANY-stmt', $stmt);
            $stmt->replaceNode ($stmt1);
          }
        else
          {
            $stmt->unbindNode ();
          }
      }
  }

my @if_construct = &F ('.//if-construct[./if-block/ANY-stmt/condition-E[string(.)=".TRUE." or string(.)=".FALSE."]]', $d);

for my $i (0 .. $#if_construct)
  {
    my $if_construct = $if_construct[$i];


    my ($file) = &F ('ancestor::file', $if_construct);
    next unless ($file);

    my $str = $if_construct->textContent;

    if (my ($if_block) = &F ('./if-block[./ANY-stmt/condition-E[string(.)=".TRUE."]]', $if_construct))
      {
        if (my @if_block = &F ('following-sibling::if-block', $if_block))
          {
            for (@if_block)
              {
                $_->unbindNode ();
              }
            $if_block->appendChild (&n ('<end-if-stmt>ENDIF</end-if-stmt>'));
          }
      }

    if (my @if_block = &F ('./if-block[./ANY-stmt/condition-E[string(.)=".FALSE."]]', $if_construct))
      {
        for my $if_block (@if_block)
          {
            if (my $prev = $if_block->previousSibling)
              {
                if ($if_block->nextSibling)
                  {
                    $prev->appendChild (&n ('<end-if-stmt>ENDIF</end-if-stmt>'));
                  }
                $if_block->unbindNode ();
              }
            elsif (my $next = $if_block->nextSibling)
              {
                my ($stmt) = &F ('./ANY-stmt', $next);
                if ($stmt->nodeName eq 'else-stmt')
                  {
                    $if_block->unbindNode ();
                    $stmt->replaceNode (&n ("<if-then-stmt>IF (<condition-E><literal-E>.TRUE.</literal-E></condition-E>) THEN</if-then-stmt>"));
                  }
                elsif ($stmt->nodeName eq 'else-if-stmt')
                  {
                    $if_block->unbindNode ();
                    $stmt->setNodeName ('if-then-stmt');
                    my $tt = $stmt->firstChild;
                    $tt->replaceNode (&t ("IF ("));
                  }
              }
          }
      }

    my @if_block = &F ('./if-block', $if_construct);

    if (scalar (@if_block) == 1)
      {
        my $if_block = $if_block[0];
        if (&F ('.//ANY-stmt/condition-E[string(.)=".TRUE."]', $if_block))
          {
            $if_block->firstChild->unbindNode ();
            $if_block->lastChild->unbindNode ();
          }
        elsif (&F ('.//ANY-stmt/condition-E[string(.)=".FALSE."]', $if_block))
          {
            $if_construct->replaceNode (&t (''));
          }
      }


    if ($str ne $if_construct->textContent)
      {
        my $iii = sprintf ('%3.3d', $i);
        'FileHandle'->new (">if_construct.$iii.F90.old")->print ($str);
        'FileHandle'->new (">if_construct.$iii.F90.new")->print ($if_construct->textContent);
      }


  }


my @include = &F ('//include', $d);

for my $include (@include)
  {
    my ($sub) = &F ('./filename', $include, 2);
    for ($sub)
      {
       s/\.intfb\.h$// ; s/\.h$//o; $_ = uc ($_) 
      }
    next unless ($sub =~ m/^\w+$/o);
    next if (my @call = &F ('//call-stmt[string(procedure-designator)="?"]', $sub, $d));
    $include->unbindNode ();
  }





'FileHandle'->new (">$F90.2.new")->print ($d->textContent);

__END__
RDECRD = 20000.
LMPHYS = TRUE
LOZONE = FALSE
LO3FL = TRUE   
NOZOCL = 1   
LRAYFM = TRUE
NRADFR = 15
LNEWSTAT = TRUE
LRPROX = FALSE
LMDUST = FALSE    
NGFL_EXT = 0
LCONDWT = TRUE
L3MT = FALSE    
LSTRAPRO = FALSE  
LPROCLD = TRUE
LGRAPRO = FALSE
LGPCMT = FALSE
LCVPRO = FALSE
LTHERMO = TRUE
LSFORCS = FALSE
LMSE = TRUE
LSOLV = FALSE
NSWB_MNH = 6  
NSW = 6  
LRAYFM15 = FALSE
LRAY = FALSE
LCOEFKTKE = FALSE 
LCOEFK_TOMS = FALSE
LHMTO = TRUE
LVDIF = TRUE  
LGWD = TRUE
LCOEFKSURF = FALSE
LPTKE = FALSE
LCOEFK_PTTE = FALSE
LRAFTUR = TRUE
LNEBR = FALSE  
LECT = TRUE
CGMIXLEN = 'AY'
LFLUSO = TRUE
LCVPPKF = TRUE
LEDKF = FALSE
LEDMFI = FALSE  
LECSHAL = TRUE
LECDEEP = TRUE
LSNV = FALSE
LVGSN = FALSE
LZ0HSREL = TRUE  
LALBMERCLIM = FALSE
NAER = 1 
LAEROSUL = FALSE
LAEROVOL = FALSE
LRAYSP = FALSE  
LRSTAER = TRUE
LCVPGY = FALSE
LNCVPGY = FALSE
LNEBECT = FALSE
LCVCSD = FALSE
LMUSCLFA = FALSE
LRCOEF = FALSE
NCALLRAD = 2
LRAYLU = TRUE
NRAY = 1
LXMRT = FALSE
LACDIFUS = FALSE
LELAM = FALSE
LDIFCONS = TRUE   
LNOIDFQC = TRUE
L3DTURB = FALSE
LCVGQD = FALSE
LVDIFSPNL = FALSE 
LSTRA = FALSE
LSTARS = FALSE
LCVRA = FALSE  
LCVRAV3 = FALSE
LCOMOD = TRUE
LCVTDK = TRUE
LCAPE = FALSE
LCAMOD = FALSE
LCVPRO = FALSE
LRCVOTT = FALSE
LCDDPRO = FALSE
LNSDO = FALSE
LADJCLD = TRUE
LNEBCO = FALSE
LGWDC = FALSE
LNORGWD = FALSE
LFPCOR = TRUE
LSFHYD = TRUE
LRRMES = TRUE
LPHSPSH = FALSE
LXVISI = TRUE 
LXVISI2 = TRUE
LDPRECIPS = TRUE  
LDPRECIPS2 = TRUE
NDTPRECCUR = 1
