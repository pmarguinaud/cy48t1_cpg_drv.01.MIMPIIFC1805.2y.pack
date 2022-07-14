#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

for my $F90 (@ARGV)
  {

    my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 400)]);
    
    my $f = &Fxtran::fxtran (location => 'fypp/field_variables_mod.F90', fopts => [qw (-line-length 200)]);
    
    my ($ft) = &F ('//T-construct[.//T-stmt[string (T-N)="FIELD_VARIABLES"]]', $f);
    my @fn = &F ('./component-decl-stmt//EN-N//text()', $ft, 1);
    
    my @expr = &F ('//named-E[string(N)="PGMV" or string(N)="PGFL" or string(N)="PGMVS" or string(N)="PGMVT1"]', $d);

    for my $expr (@expr)
      {
#       print $expr->textContent, "\n";
        my ($r) = &F ('./R-LT/ANY-R', $expr);
        goto NEXT unless ($r);

        my $ind;
        my $indp;
        if ($r->nodeName eq 'array-R')
          {
            ($ind) = &F ('./section-subscript-LT/section-subscript[last()]/lower-bound/ANY-E', $r);
            $indp = $ind->parentNode->parentNode;

#           print $ind->textContent, "\n";
          }
        elsif ($r->nodeName eq 'parens-R')
          {
            ($ind) = &F ('./element-LT/element[last()]/ANY-E', $r);
            $indp = $ind->parentNode;

#           print $ind->textContent, "\n";
          }

        my ($N) = &F ('./N', $expr);

        if ($N->textContent =~ m/^PGMV/o)
          {
            my ($t, $y) = split (m/%/o, $ind->textContent);
            $t =~ s/^Y//o;
            $y =~ s/^M//o; 
            if ($y =~ s/(L|M)$//o)
              {
                my $d = $1;
                next unless (($t eq 'T0') || ($t eq 'T9'));
                if ($t eq 'T0')
                  {
                    $t = "D$d";
                  }
                else
                  {
                    $t = "D${d}9";
                  }
              }
    

            if (($y eq 'SGRT') && ($t =~ m/^(?:DL|DM)9?$/o))
              {
                my ($d, $T) = ($t =~ m/^D(M|L)(9)?$/o);
                $T ||= '0';
                $y = "$y$d";
                $t = "T$T";
              }

            $N->replaceNode (&t ("YDVARS%$y%$t"));
          }
        elsif ($N->textContent eq 'PGFL')
          {
            my ($y, $t) = split (m/%/o, $ind->textContent);
            $t =~ s/^MP//o or die "t=$t";
            $t ||= '0';
            $y =~ s/^Y//o;
    
            next if ($y =~ m/\(/o);
            die "$y" unless (grep { $_ eq $y } @fn);
    
            if ($t =~ m/^\d$/o)
              {
                $t = "T$t";
              }
            elsif (($t eq 'L') or ($t eq 'M'))
              {
                $t = "D$t";
              }
            else
              {
                die;
              }
    
            $N->replaceNode (&t ("YDVARS%$y%$t"));
    

          }

        $indp->previousSibling->unbindNode;
        $indp->unbindNode;

        $r->unbindNode if ($r->textContent eq '(:,:)');
        $r->unbindNode if ($r->textContent eq '(:)');

NEXT:
#       print "\n";
      }
 

    @expr = &F ('.//named-E[string(N)="PB1"][R-LT]', $d);

    for my $expr (@expr) 
      {
        my ($r) = &F ('./R-LT/parens-R', $expr);
        my ($rlt) = $r->parentNode;

        die unless ($r);
        my @ind = &F ('./element-LT/element/ANY-E', $r);
        die $r unless (scalar (@ind) == 2);
        die $expr unless (($ind[0]->textContent eq 'JROF') && ($ind[1]->nodeName eq 'op-E'));
        my @e = &F ('.//named-E', $ind[1]);
        die $ind[1] unless (scalar (@e) == 3);
        die unless (($e[1]->textContent eq 'JLEV') && ($e[2]->textContent eq 'NFLSA'));
        my $M = $e[0]->textContent;
        die unless ($M =~ s/^MSLB1//o);

        my ($N) = &F ('./N/n/text()', $expr);

        $N->setData ('YDCPG_SL1');
        $rlt->insertBefore (&n ("<component-R>%<ct>$M</ct></component-R>"), $rlt->firstChild);
        $ind[1]->replaceNode (&n ('<named-E><N><n>JLEV</n></N></named-E>'));
      }

    
    @expr = &F ('.//named-E[string(N)="PB2"][R-LT]', $d);

    for my $expr (@expr) 
      {
        my ($r) = &F ('./R-LT/parens-R', $expr);
        my ($rlt) = $r->parentNode;

        die unless ($r);
        my @ind = &F ('./element-LT/element/ANY-E', $r);
        die $r unless (scalar (@ind) == 2);
        if ($ind[1]->nodeName eq 'named-E')
          {
            die &Dumper ([$expr->textContent, [map { $_->textContent } @ind]]) 
              unless (($ind[0]->textContent eq 'JROF') && ($ind[1]->nodeName eq 'named-E'));

            my $M = $ind[1]->textContent;
            die unless ($M =~ s/^MSLB2//o);

            my ($N) = &F ('./N/n/text()', $expr);

            $N->setData ('YDCPG_SL2');
            $rlt->insertBefore (&n ("<component-R>%<ct>$M</ct></component-R>"), $rlt->firstChild);
            $ind[1]->parentNode->previousSibling->unbindNode;
            $ind[1]->parentNode->unbindNode;
          }
        else
          {
            die &Dumper ([$expr->textContent, [map { $_->textContent } @ind]]) 
              unless (($ind[0]->textContent eq 'JROF') && ($ind[1]->nodeName eq 'op-E'));
            my @e = &F ('.//named-E', $ind[1]);
            die $ind[1] unless (scalar (@e) == 2);
            die unless ($e[1]->textContent eq 'JLEV');
            my $M = $e[0]->textContent;
            die unless ($M =~ s/^MSLB2//o);

            my ($N) = &F ('./N/n/text()', $expr);

            $N->setData ('YDCPG_SL2');
            $rlt->insertBefore (&n ("<component-R>%<ct>$M</ct></component-R>"), $rlt->firstChild);
            $ind[1]->replaceNode (&n ('<named-E><N><n>JLEV</n></N></named-E>'));
          }
      }

    
#   print $d->textContent;
    
    'FileHandle'->new (">$F90.new")->print ($d->textContent);
    

  }
