#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use List::MoreUtils qw (uniq all);
use lib $Bin;
use Fxtran;

my @obj = qw (YDMF_PHYS_BASE_STATE YDMF_PHYS_NEXT_STATE YDCPG_MISC YDCPG_PHY9
              YDCPG_PHY0 YDMF_PHYS YDCPG_DYN9 YDCPG_DYN0 YDMF_PHYS_SURF YDVARS);

sub parseDirectives
{
  my $d = shift;

  my @C = &F ('//C[starts-with(string (.),"!=")]', $d);
  
  while (my $C  = shift (@C))
    {
      (my $bdir = $C->textContent) =~ s/^!=\s*//o;
      $bdir = lc ($bdir);
      my ($tag) = ($bdir =~ m/^(\w+)/o);
  
      my @node;
      for (my $node = $C->nextSibling; ; $node = $node->nextSibling)
        {
          $node or die $C->textContent;
          if (($node->nodeName eq 'C') && (index ($node->textContent, '!=') == 0))
            {
              my $C = shift (@C);
              (my $edir = $C->textContent) =~ s/^!=\s*//o;
              $edir = lc ($edir);

              die unless ($edir =~ s/^end\s+//o);
              die unless ($edir eq $tag);

              $C->unbindNode ();
              
              last;
            }

          push @node, $node;

        }

      my $e = &n ("<$tag-directive/>");
 
      for my $node (@node)
        {
          $e->appendChild ($node);
        }

      $C->replaceNode ($e);

    }
}

sub serial
{
  shift;
}

sub parallel
{
  shift;
  my ($par, $doc) = @_;

  my $indent = &Fxtran::getIndent ($par);

  my $str = ' ' x $indent;

  my ($loop) = &Fxtran::fxtran (fragment => << "EOF");
DO JBLK = 1, YDCPG_DIM%KGPBLKS

${str}  YLCPG_DIM = YDCPG_DIM
${str}  CALL YLCPG_DIM%UPDATE (JBLK)
${str}ENDDO
EOF

  my ($enddo) = &F ('.//end-do-stmt', $loop);
  my $p = $enddo->parentNode;

  for my $node ($par->childNodes ())
    {
      $p->insertBefore (&t (' ' x ($indent + 2)), $enddo);
      &Fxtran::reIndent ($node, $indent + 2);
      $p->insertBefore ($node, $enddo);
    }
  
  $par->appendChild ($loop);

  my @expr = &F ('.//named-E/N/n[string(.)="YDCPG_DIM"]/text()', $par);

  shift (@expr) for (1 .. 2);

  for my $expr (@expr)
    {
      $expr->setData ('YLCPG_DIM');
    }

}

sub skip
{
  shift;
  $_[0]->unbindNode ();
}

sub processDirectives
{
  my $doc = shift;

  my @dir = &F ('.//ANY-directive', $doc);

  for my $dir (@dir)
    {
      (my $name = $dir->nodeName) =~ s/-directive$//o;
      __PACKAGE__->can ($name) or die ("Unknown directive $name");
      __PACKAGE__->$name ($dir, $doc);
    }

}

sub addBlkDimensionToLocals
{
  my $doc = shift;

  my @skip = qw (PGFL PGFLT1 PGMVT1);
  my %skip = map { ($_, 1) } @skip;

  # Prepare (:,:,JBLK) references

  my @ref;
  
  for my $i (1 .. 3)
    {
      my ($v) = &Fxtran::fxtran (expr => 'X(' . join (',', (':') x ($i-1), 'JBLK') . ')');
      my ($r) = &F ('./R-LT/ANY-R', $v);
      $ref[$i] = $r;
    }
  
  for my $en_decl (&F ('.//EN-decl[./array-spec/shape-spec-LT/shape-spec[string(.)="YDCPG_DIM%KLON"]]', $doc))
    {
      my $decl = &Fxtran::stmt ($en_decl);

      my ($name) = &F ('./EN-N', $en_decl, 1);
     
      next if ($skip{$name});

      my @par = &F (
                    './/parallel-directive[.//named-E[string(N)="?"]]|'
                  . './/parallel-call-directive[.//named-E[string(N)="?"]]', 
                    $name, $name, $doc);

      my ($intent) = &F ('.//attribute-N[string(.)="INTENT"]', $decl);

      next if ((! $intent) && scalar (@par) < 2);

      my ($ssl) = &F ('./array-spec/shape-spec-LT', $en_decl);
      my ($ngpblks) = &Fxtran::fxtran (expr => 'YDCPG_DIM%KGPBLKS');
      $ssl->appendChild (&t (','));
      $ssl->appendChild (&n ('<shape-spec><upper-bound>' . $ngpblks . '</upper-bound></shape-spec>'));
  
      my @ss = &F ('./shape-spec', $ssl);
  
      my @expr = &F ('.//named-E[string(N)="?"]', $name, $doc);
      for my $expr (@expr)
        {
          my ($par) = &F ('ancestor::parallel-directive', $expr);
          next unless ($par);
          my ($r) = &F ('./R-LT/ANY-R', $expr);
          if (! $r)
            {
              # No reference was found, add one
              $expr->appendChild (my $rlt = &n ('<R-LT/>'));
              my $r = $ref[scalar (@ss)]->cloneNode (1);
              $rlt->appendChild ($r);
            }
          elsif ($r->nodeName eq 'parens-R')
            {
              my ($elt) = &F ('./element-LT', $r);
              $elt->appendChild (&t (','));
              $elt->appendChild (&n ('<named-E><N><n>JBLK</n></N></named-E>'));
            }
          else
            {
              # An array reference was found; add JBLK index
              die $expr;
            }
        }
    }

}

sub addBlkDimensionToObjects
{
  my $doc = shift;

  my @par = &F ('.//parallel-directive/do-construct', $doc);

  my $h = do ('./h.pl');

  my $H = 'XML::LibXML'->load_xml (location => 'h.xml');

  my %decl;

  my %ptr;

  for my $par (@par)
    {
      my ($do) = $par->firstChild;

      my $indent = &Fxtran::getIndent ($do);

      my @o;

      for my $obj (@obj)
        {

          my @expr = &F ('.//named-E[string(N)="?"]', $obj, $par);
          next unless (@expr);
         
          my ($typ) = &F ('.//T-decl-stmt[.//EN-decl[string(EN-N)="?"]]/_T-spec_/derived-T-spec/T-N', $obj, $doc, 1);

          my %p;

          for my $expr (@expr)
            {
              my ($name) = &F ('./N/n/text()', $expr); 
              my @ct = &F ('./R-LT/component-R/ct', $expr, 1);
              my $e = $expr->cloneNode (1);
              my @r = &F ('./R-LT/component-R', $expr);
              $_->unbindNode for (@r);
              my $ptr = join ('_', 'Z', $obj, @ct);
              $name->setData ($ptr);
              $p{$ptr} = {ctl => \@ct};
            }

          for my $ptr (sort keys (%p))
            {
              $ptr{$ptr} = $p{$ptr};

              my $ctl = $ptr{$ptr}{ctl};
              my $key = join ('%', $typ, @$ctl);
              $ptr{$ptr}{key} = $key;

              my $decl;
              unless ($decl = $decl{$key}) 
                {
                  ($decl) = &F ('./list/decl[@key="?"]/T-decl-stmt', $key, $H); 
                  $decl{$key} = $decl;
                }

              my @ss = &F ('.//shape-spec', $decl);
              my @lb;
              for my $ss (@ss)
                {
                  my ($lb) = &F ('./lower-bound', $ss);
                  push @lb, $lb ? $lb->textContent : '1';
                }

              
              my $stmt;

              my @ctl = @$ctl;

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

              my $p = join ('%', $obj, @ctl);
              my $d = (all { $_ eq '1' } @lb) ? '' : "(" . join (',', map { "$_:" } @lb) . ")";
              $stmt = &Fxtran::fxtran (statement => "IF (ASSOCIATED ($p)) $ptr $d => $p%GET_VIEW (JBLK)", fopts => [qw (-line-length 300)]);

              $par->insertAfter ($stmt, $do);
              $par->insertAfter (&t ("\n" . (' ' x ($indent +2))), $do);
            }

          

       }
    }

  my %dim;

  for my $ptr (sort keys (%ptr))
    {
      my $key = $ptr{$ptr}{key};
      my $decl = $decl{$key}->cloneNode (1);
      my ($ts) = &F ('./_T-spec_', $decl);

      my @ss = &F ('.//shape-spec', $decl);
      my ($name) = &F ('.//EN-N/N/n/text()', $decl);
 
      $name->setData ($ptr);
 
      for my $ss (@ss)
        {
          $ss->replaceNode (&n ('<shape-spec>:</shape-spec>'));
        }

      for my $attr (qw (POINTER CONTIGUOUS))
        {
          $decl->insertAfter (&n ("<attribute><attribute-N>$attr</attribute-N></attribute>"), $ts);
          $decl->insertAfter (&t (', '), $ts);
        }
      
      &addVariable ($doc, $decl);

      &addInit ($doc, "$ptr => ZDUM" . scalar (@ss));

      $dim{scalar (@ss)} = 1;

    }

  for my $dim (sort keys (%dim))
    {
      &addVariable ($doc, "REAL (KIND=JPRB), SAVE, TARGET :: ZDUM$dim (" . join (',', ('1') x $dim) . ")");
    }

}

sub renameSubroutine
{
  my ($doc) = @_;

  my $suf = '_OPENMP';

  my @name = &F ('.//subroutine-N/N/n/text()', $doc);
  my $name = $name[0]->textContent;

  for (@name)
    {
      $_->setData ($_->textContent . $suf);
    }

  my @drhook = &F ('.//call-stmt[string(procedure-designator)="DR_HOOK"]', $doc);

  for my $drhook (@drhook)
    {
      next unless (my ($S) = &F ('./arg-spec/arg/string-E/S/text()', $drhook));
      my $str = $S->textContent;
      $str =~ s/$name/$name$suf/o;
      $S->setData ($str);
    }
  
}

sub callOpenMPRoutines
{
  my $doc = shift;
  my @call = &F ('.//call-stmt[not(ancestor::parallel-directive)]', $doc);

  my $suf = '_OPENMP';

  my %proc;

  for my $call (@call)
    {
      next if (&F ('.//procedure-designator/named-E/R-LT', $call)); # Skip objects calling methods
      my ($proc) = &F ('./procedure-designator/named-E/N/n/text()', $call);

      my $par = &n ('<parallel-call-directive/>');
      $call->replaceNode ($par);
      $par->appendChild ($call);
      

      my @arg = &F ('./arg-spec/arg/named-E/N/n/text()', $call, 1);

      for my $arg (@arg)
        {
          goto FOUND if (grep { $_ eq $arg } @obj);
          my @ss = &F ('.//EN-decl[string(EN-N)="?"]/array-spec//shape-spec', $arg, $doc, 1);
          goto FOUND if (@ss && ('YDCPG_DIM%KLON' eq $ss[0]));
        }

      next;
FOUND:
      

      $proc{$proc->textContent} = 1;
      $proc->setData ($proc->textContent . $suf);
    }

  my @proc = sort keys (%proc);

  for my $proc (@proc)
    {
      next unless (my ($include) = &F ('.//include/filename[string(.)="?"]/text()', lc ($proc) . '.intfb.h', $doc));
      (my $str = $include->textContent) =~ s{^(\w+)\.intfb.h$}{$1 . lc ($suf) . '.intfb.h'}exo;
      $include->setData ($str);
    }

}

sub cleanInterfaces
{
  my $doc = shift;
  my @include = &F ('.//include', $doc);
  for my $include (@include)
    {
      my ($filename) = &F ('./filename', $include, 2);
      next unless ($filename =~ m/^(\w+)\.intfb.h$/o);
      my $proc = uc ($1);
      my @call = &F ('.//call-stmt[string(procedure-designator)="?"]', $proc, $doc);
      next if (@call);
      $include->unbindNode ();
    }
}

sub cleanParallelDirectives
{
  my $doc = shift;
  my @par = &F ('.//parallel-directive|.//parallel-call-directive|.//serial-directive', $doc);
  for my $par (@par)
    {
      for my $n ($par->childNodes ())
        {
          $par->parentNode->insertBefore ($n, $par);
        }
      $par->unbindNode ();
    }
}

sub addInit
{
  my ($doc, $decl) = @_;
  
  my ($drhook_call) = &F ('.//call-stmt[string(.//procedure-designator)="DR_HOOK"]', $doc);

  my $stmt = ref ($decl) ? $decl : &Fxtran::fxtran (statement => $decl);
  
  $drhook_call->parentNode->insertAfter ($stmt, $drhook_call);
  $drhook_call->parentNode->insertAfter (&t ("\n"), $drhook_call);
}

sub addVariable
{
  my ($doc, $decl, $cr) = @_;
  my ($zhook_decl) = &F ('.//T-decl-stmt[string(.//EN-N)="ZHOOK_HANDLE"]', $doc);

  my $stmt = ref ($decl) ? $decl : &Fxtran::fxtran (statement => $decl);

  $zhook_decl->parentNode->insertBefore ($stmt, $zhook_decl);
  $zhook_decl->parentNode->insertBefore (&t ("\n"), $zhook_decl);
  $zhook_decl->parentNode->insertBefore (&t ("\n"), $zhook_decl) if ($cr);
}

sub addVariables
{
  my $doc = shift;
  &addVariable ($doc, 'INTEGER(KIND=JPIM) :: JBLK', 1);
  &addVariable ($doc, 'TYPE(CPG_DIM_TYPE) :: YLCPG_DIM', 1);
}

sub reduceVariableScope
{
  my $doc = shift;
  my @en_decl = &F ('.//T-decl-stmt[not(string(.//attribute-N)="INTENT")]//EN-decl'
                  . '[.//shape-spec[string(.)="YDCPG_DIM%KGPBLKS"]]'
                  , $doc);
  
  for my $en_decl (@en_decl)
    {
      my ($name) = &F ('./EN-N', $en_decl, 1);

      my ($stmt) = &Fxtran::stmt ($en_decl);
      my @ss = &F ('.//shape-spec', $en_decl);
      my ($ts) = &F ('./_T-spec_', $stmt);
      $stmt->insertAfter (&n ('<attribute><attribute-N>ALLOCATABLE</attribute-N></attribute>'), $ts);
      $stmt->insertAfter (&t (', '), $ts);

      my @SS = map { $_->textContent } @ss;

      for my $ss (@ss)
        {
          $ss->replaceNode (&t (':'));
        }

      my @par = &F ('.//parallel-directive[.//named-E[string(N)="?"]]|'
                  . './/parallel-call-directive[.//named-E[string(N)="?"]]', $name, $name, $doc);

      @par or die $name;

      my $par1 = $par[0]; my $par2 = $par[-1];     


      ($par2) = &F ('ancestor-or-self::*[not(.//if-construct)]', $par2);

      my $ind1 = &Fxtran::getIndent ($par1);
      my $ind2 = &Fxtran::getIndent ($par2);

      $par1->insertBefore (&t ("\n" . (' ' x $ind1)), $par1->firstChild) unless ($par1->firstChild->nodeName eq 'allocate-stmt');
      $par1->insertBefore (&t ("\n" . (' ' x $ind1)), $par1->firstChild);
      $par1->insertBefore (&Fxtran::fxtran (statement => "IF (.NOT. ALLOCATED ($name)) ALLOCATE ($name (" . join (',', @SS)  . "))"), $par1->firstChild);

      $par2->insertAfter (&t ("\n" . (' ' x $ind2)), $par2->lastChild) unless ($par2->lastChild->nodeName eq 'deallocate-stmt');
      $par2->insertAfter (&t ("\n" . (' ' x $ind2)), $par2->lastChild);
      $par2->insertAfter (&Fxtran::fxtran (statement => "IF (ALLOCATED ($name)) DEALLOCATE ($name)"), $par2->lastChild);
 
    }

}

sub addOpenMPDirectives
{
  my $doc = shift;
  my @do = &F ('.//parallel-directive/do-construct', $doc);
  
  my %obj = map { ($_, 1) } @obj;

  for my $do (@do)
    {
      my @var = &uniq (&F ('.//named-E/N', $do, 1));

      my %prv;

      for my $var (@var)
        {
          my ($decl) = &F ('.//T-decl-stmt[.//EN-decl[string(EN-N)="?"]]', $var, $doc);
          next unless ($decl);

          my ($arg) = &F ('.//attribute-N[string(.)="INTENT"]', $decl);
          my ($dim) = &F ('.//EN-decl[string(EN-N)="?"]/array-spec', $var, $decl);
          my ($blk) = $dim && &F ('.//shape-spec[string(.)="YDCPG_DIM%KGPBLKS"]', $dim);
          my ($lon) = $dim && &F ('.//shape-spec[string(.)="YDCPG_DIM%KLON"]', $dim);
          next if ($arg || $blk);

          my ($stmt) = &F ('.//a-stmt[./E-1/named-E[string(N)="?"]]', $var, $do);
          my ($loop) = &F ('.//do-V/named-E[string(N)="?"]', $var, $do);

          next unless ($stmt || $loop || $lon);

          $prv{$var}++;
        }


      my @prv = sort keys (%prv);


#     print join ('|', @prv), "\n";
#     print &Dumper ([\@prv]);
#     print $do->textContent, "\n" x 2;
#     print "\n\n", '-' x 80, "\n\n";
    }
}

sub removeUnusedArrays
{
  my $doc = shift;
  my @en_decl = &F ('.//T-decl-stmt[not(.//attribute[string(.)="INTENT"])]//EN-decl[.//shape-spec[string(.)="YDCPG_DIM%KLON"]]', $doc);

return;

  for my $en_decl (@en_decl)
    {
      my ($var) = &F ('./EN-N', $en_decl, 1);
      my @expr = &F ('.//named-E[string(N)="?"]', $var, $doc);
      next if (@expr);
      my $stmt = &Fxtran::stmt ($en_decl);
      $stmt->unbindNode ();
    }
}

my $F90 = shift;

my $doc = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);

&renameSubroutine ($doc);

&parseDirectives ($doc);

&processDirectives ($doc);

&callOpenMPRoutines ($doc);

&addBlkDimensionToLocals ($doc);

&addBlkDimensionToObjects ($doc);

&cleanInterfaces ($doc);

&addVariables ($doc);

&addOpenMPDirectives ($doc);

&reduceVariableScope ($doc);

&cleanParallelDirectives ($doc);

&removeUnusedArrays ($doc);


$F90 =~ s/.F90$/_openmp.F90/o;

'FileHandle'->new (">$F90")->print ($doc->textContent);

