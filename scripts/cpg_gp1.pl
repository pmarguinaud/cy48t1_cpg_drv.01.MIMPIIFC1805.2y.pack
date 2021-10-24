#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

for my $F90 ("src/local/arpifs/adiab/cpg_gp.F90")
  {

    my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 400)]);
    
    my $f = &Fxtran::fxtran (location => 'fypp/field_variables_mod.F90', fopts => [qw (-line-length 200)]);
    
    my ($ft) = &F ('//T-construct[.//T-stmt[string (T-N)="FIELD_VARIABLES"]]', $f);
    my @fn = &F ('./component-decl-stmt//EN-N//text()', $ft, 1);
    
    my @expr = &F ('//named-E[string(N)="PGMV" or string(N)="PGFL" or string(N)="PGMVS"]', $d);

    for my $expr (@expr)
      {
        print $expr->textContent, "\n";
        my ($r) = &F ('./R-LT/ANY-R', $expr);
        goto NEXT unless ($r);

        my $ind;
        my $indp;
        if ($r->nodeName eq 'array-R')
          {
            ($ind) = &F ('./section-subscript-LT/section-subscript[last()]/lower-bound/ANY-E', $r);
            $indp = $ind->parentNode->parentNode;

            print $ind->textContent, "\n";
          }
        elsif ($r->nodeName eq 'parens-R')
          {
            ($ind) = &F ('./element-LT/element[last()]/ANY-E', $r);
            $indp = $ind->parentNode;

            print $ind->textContent, "\n";
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
        print "\n";
      }
 

=pod

    for my $call (@call)
      {
        my ($argspec) = &F ('./arg-spec', $call);
        my @arg = &F ('./arg/ANY-E', $argspec);
        next unless (scalar (@arg) == 3);
    
        my $stmt;
    
        if ($arg[1]->textContent =~ m/^(?:PSP_(?:RR|SB|SG)|PSD_(?:VF|VA|VD|VH|VK|VP|VV))$/o)
          {
            my ($p, $x, $v) = map ({ $_->textContent } @arg);
            $x =~ s/^PS(P|D)_//o; my $z = $1;
            $p =~ s/^.*%Y//o;
    
            my ($f, $t) = split (m/%/o, $p);
            $t =~ s/^MP/T/o;
            $t = "T0" if ($t eq 'T');

            $t = '' if ($z eq 'D');
    
            $t = $t ? "_$t" : "";
            
            my @ss = &F ('//EN-decl[string (EN-N)="?"]//shape-spec', $v, $d);
#           print &Dumper ([$f, $t, $x, $v, scalar (@ss)]);
    
    
            $stmt = &n ("<pointer-a-stmt><E-1><named-E><N><n>$v</n></N></named-E></E-1> "
                      . "<a>=&gt;</a> <E-2>"
                      . "<named-E><N><n>YDMF_PHYS_SURF</n></N><R-LT><component-R>%<ct>GS${z}_$x</ct></component-R>"
                      . "<component-R>%<ct>P${f}$t</ct></component-R></R-LT></named-E>"
                      . "</E-2></pointer-a-stmt>");
    
    
          }
        elsif ($arg[1]->textContent =~ m/^PGMV/o)
          {
            (my $p, undef, my $v) = map ({ $_->textContent } @arg);


            my ($t, $y) = split (m/%/o, $p);
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

            $stmt = &n ("<pointer-a-stmt><E-1><named-E><N><n>$v</n></N></named-E></E-1> <a>=&gt;</a> "
                      . "<E-2><named-E><N><n>YDVARS</n></N><R-LT><component-R>%<ct>$y</ct></component-R>"
                      . "<component-R>%<ct>$t</ct></component-R></R-LT></named-E></E-2></pointer-a-stmt>");
    
          }
        elsif ($arg[1]->textContent =~ m/^PGFL/o)
          {
            (my $p, undef, my $v) = map ({ $_->textContent } @arg);
            my ($y, $t) = split (m/%/o, $p);
            $t =~ s/^MP//o or die "t=$t";
            $t ||= '0';
            $y =~ s/^Y//o;
    
    
            die unless (grep { $_ eq $y } @fn);
    
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
    
            $stmt = &n ("<pointer-a-stmt><E-1><named-E><N><n>$v</n></N></named-E></E-1> <a>=&gt;</a> "
                      . "<E-2><named-E><N><n>YDVARS</n></N><R-LT><component-R>%<ct>$y</ct></component-R>"
                      . "<component-R>%<ct>$t</ct></component-R></R-LT></named-E></E-2></pointer-a-stmt>");
    
          }
    
        if ($stmt)
          {
            my ($E1) = &F ('./E-1', $stmt);
            my $sp = $E1->nextSibling;
            $sp->setData (' ' x (12 - length ($E1->textContent)));
            $call->replaceNode ($stmt);
          }


      }
        
=cut
    
    'FileHandle'->new (">$F90.new")->print ($d->textContent);
    

  }
