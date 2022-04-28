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
my @skip = qw (PGFL PGFLT1 PGMVT1 PGPSDT2D);
my %skip = map { ($_, 1) } @skip;


sub removeListElement
{

# Remove element from list, take care of removing comma before or after the element

  my $x = shift;

  my $nn = $x->nodeName;

  my ($p) = $x->parentNode;
  
  my @cf = &F ('following-sibling::text()[contains(.,",")]', $x);   
  my @cp = &F ('preceding-sibling::text()[contains(.,",")]', $x);   
  
  if (@cf)
    {   
      $cf[+0]->unbindNode (); 
    }   
  elsif (@cp)
    {   
      $cp[-1]->unbindNode (); 
    }   
  
  $x->parentNode->appendChild (&t (' '));
  my $l = $x->parentNode->lastChild;
  
  $x->unbindNode (); 
  
  while ($l)
    {   
      last if (($l->nodeName ne '#text') && ($l->nodeName ne 'cnt'));
      $l = $l->previousSibling;
      last unless ($l);
      $l->nextSibling->unbindNode;
    }   

  return &F ("./$nn", $p) ? 0 : 1;
}

sub singleDecl
{

# Single declaration statement per entity

  my $doc = shift;

# Select all entity lists with several entities

  my @en_decl_lst = &F ('.//EN-decl-LT[count(./EN-decl)>1]', $doc);

  for my $en_decl_lst (@en_decl_lst)
    {
      my $stmt = &Fxtran::stmt ($en_decl_lst);
      my $indent = &Fxtran::getIndent ($stmt);
      my @en_decl = &F ('./EN-decl', $en_decl_lst);
      for my $en_decl (@en_decl)
        {
          my $s = $stmt->cloneNode (1);
          my ($l) = &F ('.//EN-decl-LT', $s);
          for ($l->childNodes ())
            {
              $_->unbindNode ();
            }
          $l->appendChild ($en_decl->cloneNode (1));
          $stmt->parentNode->insertAfter ($s, $stmt);
          $stmt->parentNode->insertAfter (&t ("\n" . (' ' x $indent)), $stmt);
        }
      $stmt->unbindNode ();
    }

}


sub parseDirectives
{

# Add tags for each section

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

# Do nothing

  shift;
}

sub parallel
{

# Add a loop on blocks

  shift;
  my ($par, $doc) = @_;

  my ($stmt) = &F ('.//ANY-stmt', $par);

  my $indent = &Fxtran::getIndent ($stmt);


  my $str = ' ' x $indent;

  my ($loop) = &Fxtran::fxtran (fragment => << "EOF");
DO JBLK = 1, YDCPG_OPTS%KGPBLKS
${str}  YLCPG_BNDS = YDCPG_BNDS
${str}  CALL YLCPG_BNDS%UPDATE (JBLK)
${str}ENDDO
EOF


  my ($enddo) = &F ('.//end-do-stmt', $loop);
  my $p = $enddo->parentNode;

  for my $node ($par->childNodes ())
    {
      $p->insertBefore (&t (' ' x (2)), $enddo);
      &Fxtran::reIndent ($node, 2);
      $p->insertBefore ($node, $enddo);
    }
  $p->insertBefore (&t (' ' x $indent), $enddo);
  
  $par->appendChild ($loop);

  my @expr = &F ('.//named-E/N/n[string(.)="YDCPG_BNDS"]/text()', $par);

  shift (@expr);

  for my $expr (@expr)
    {
      $expr->setData ('YLCPG_BNDS');
    }


}

sub skip
{

# Remove node

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

sub cleanParallelDirectives
{

# Remove parallel tags

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
  my ($doc, $obj) = @_;

  unless ($type{$obj})
    {
      ($type{$obj}) = &F ('.//T-decl-stmt[.//EN-decl[string(EN-N)="?"]]/_T-spec_/derived-T-spec/T-N', $obj, $doc, 1);
    }

  return $type{$obj};
}

sub getObjectDims
{
  my ($doc, $key) = @_;

  my $decl = &getObjectDecl ($key);

  return &dimsFromDecl ($decl);
}

}

sub dimsFromDecl
{
  my $decl = shift;

  my @ss = &F ('.//EN-decl/array-spec/shape-spec-LT/shape-spec', $decl);

  my @dims;

  for my $ss (@ss)
    {
      my @b = &F ('./ANY-bound/ANY-E', $ss);
      push @dims, [map { $_->cloneNode (1) } @b];
    }
  push @dims, [&n ('<named-E><N><n>YDCPG_OPTS</n></N><R-LT><component-R>%<ct>NBLKS</ct></component-R></R-LT></named-E>')];

  return @dims;
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

  return join ('%', $obj, @ctl);
}

sub addPointerAttr
{
  my $stmt = shift;
  my $ts = $stmt->firstChild;

  for my $attr (qw (CONTIGUOUS POINTER))
    {
      $stmt->insertAfter (&n ("<attribute><attribute-N>$attr</attribute-N></attribute>"), $ts);
      $stmt->insertAfter (&t (', '), $ts);
    }
}

sub getSubroutineInterface
{
  my $proc = shift;
  $proc = lc ($proc);

  my $dir = 'src/local/.intfb/arpifs';

  my $code = do { local $/ = undef; my $fh = 'FileHandle'->new ("<$dir/$proc.intfb.h"); <$fh> };

  my ($intf) = &Fxtran::fxtran (fragment => $code);

  return $intf;
}

sub getArgumentIntent
{
  my ($call, $expr) = @_;
  my @args = &F ('./arg-spec/arg/ANY-E', $call);

  my $rank;
  for my $i (0 .. $#args)
    {
      if ($expr->isEqual ($args[$i]))
        {
          $rank = $i;
          last;
        }
    }

  return unless ($rank);

  my ($proc) = &F ('./procedure-designator', $call, 1);

  my $intf = &getSubroutineInterface ($proc);

  my ($unit) = &F ('.//program-unit[./subroutine-stmt[string(subroutine-N)="?"]]', $proc, $intf);

  my ($stmt) = &F ('./subroutine-stmt', $unit);

  @args = &F ('./dummy-arg-LT/arg-N', $stmt, 1);

  die unless ($rank < @args);

  my ($intent) = &F ('.//T-decl-stmt[.//EN-decl[string(EN-N)="?"]]//intent-spec', $args[$rank], $unit, 1);

  return $intent;
}

sub makeParallel
{
  my $doc = shift;

  # Argument list
  
  my @args = &F ('.//subroutine-stmt/dummy-arg-LT/arg-N/N/n/text()', $doc);
  my %args = map { ($_->textContent, 1) } @args;

  # Objects passed as arguments
  
  my %obj = map { ($_, 1) } @obj;
  
  my %dims;     # Array dimensions
  my %arr;      # True if this is a NPROMA array
  my %arr2fld;  # Field object associated to this array

  # Local arrays: turn them into CONTIGUOUS POINTERs with an extra dimension
  # For each local array ZARR, add a FIELD object (YL_ZARR)
  # For each argument array PARR, add a FIELD object (YD_PARR)

  my @init_yl;

  my @en_decl = &F ('.//EN-decl[./array-spec/shape-spec-LT/shape-spec[string(.)="YDCPG_OPTS%KLON"]]', $doc);
  for my $en_decl (@en_decl)
    {
      my $stmt = &Fxtran::stmt ($en_decl);

      my ($N) = &F ('./EN-N', $en_decl, 1);

      next if ($skip{$N});

      $arr{$N} = 1;

      my ($sslt) = &F ('./array-spec/shape-spec-LT', $en_decl);

# TODO : Create a symbol table and use it

      my @dims = &dimsFromDecl ($stmt);
      $dims{$N} = \@dims;

# Turn array declaration into a POINTER declaration with an extra dimension

      my @ss = &F ('./shape-spec', $sslt);

      for ($sslt->childNodes ())
        {
          $_->unbindNode ();
        }
      for my $i (0 .. $#ss+1)
        {
          $sslt->appendChild (&n ('<shape-spec>:</shape-spec>'));
          $sslt->appendChild (&t (',')) if ($i <= $#ss);
        }

      &addPointerAttr ($stmt);

      my $nd = scalar (@dims);
      my $decl;

      if ($args{$N})
        {
          my ($intent) = &F ('.//intent-spec', $stmt);
          ($decl) = &Fxtran::fxtran (statement => "TYPE (FIELD_${nd}D), INTENT(" . $intent->textContent . ") :: YD_$N");

          # Remove INTENT from POINTER declaration
          ($intent) = &F ('ancestor::attribute', $intent);
          $intent->previousSibling->unbindNode ();
          $intent->unbindNode ();

          my ($arg) = grep { $_->textContent eq $N } @args;
          $arg->setData ("YD_$N");
          $arr2fld{$N} = "YD_$N";
        }
      else
        {
          ($decl) = &Fxtran::fxtran (statement => "TYPE (FIELD_${nd}D), POINTER :: YL_$N");
          $arr2fld{$N} = "YL_$N";
        }


      $stmt->parentNode->insertBefore ($decl, $stmt);
      $stmt->parentNode->insertBefore (&t ("\n"), $stmt);

      push @init_yl, &Fxtran::fxtran (statement => "$arr2fld{$N} => NULL ()") unless ($args{$N});
    }

  &addInit ($doc, @init_yl);

  my %ptr2decl; # For each pointer, the type declaration statement, with type, kind and dimensions

  my @par = &F ('.//parallel-directive', $doc);
  for my $par (@par)
    {
      my %ptr2dims; # Dimensions for each pointer
      my %ptr2vars; # Field for each pointer
      my %ptr2cond; # Field associated or not for each pointer
      my %ptr2r;    # Pointer read
      my %ptr2w;    # Pointer written

      # For each named expression in the parallel block
      my @expr = &F ('.//named-E', $par);
      for my $expr (@expr)
        {
          my ($N) = &F ('./N', $expr, 1);
          next if ($skip{$N});

          my $nd;

          # This comes from a NPROMA array
          if ($arr{$N})
            {
              $nd = scalar (@{ $dims{$N} }) - 1;
              $ptr2vars{$N} = $arr2fld{$N};                                # Field for pointer
              $ptr2dims{$N} = $dims{$N};
              $ptr2cond{$N} = 0;
            }
          # Here we deal with members of objects, such a YDVARS
          elsif ($obj{$N})
            {
              # Change expression using the POINTER instead of the object and members

              my ($name) = &F ('./N/n/text()', $expr); 
              my @ctl = &F ('./R-LT/component-R/ct', $expr, 1);
              my $e = $expr->cloneNode (1);
              my @r = &F ('./R-LT/component-R', $expr);
              $_->unbindNode for (@r);
              my $ptr = join ('_', 'Z', $N, @ctl);
              $name->setData ($ptr);

              my $key = join ('%', &getObjectType ($doc, $N), @ctl);

              my @dims = &getObjectDims ($doc, $key);

              $ptr2decl{$ptr} = &getObjectDecl ($key);
              $ptr2dims{$ptr} = \@dims;
              $ptr2vars{$ptr} = &getFieldFromObjectComponents ($N, @ctl);  # Field for pointer
              $ptr2cond{$ptr} = 1;

              $nd = scalar (@dims) - 1;

              $N = $ptr;
            }
          else
            {
              next;
            }

          # Add reference list if needed
          
          my ($rlt) = &F ('./R-LT', $expr);
          unless ($rlt)
            {
              $expr->appendChild ($rlt = &n ('<R-LT/>'));
            }

          # Add array reference if needed

          my $r = $rlt->lastChild;
          unless ($r)
            {
              $rlt->appendChild ($r = &n ('<array-R>(<section-subscript-LT>' . join (',', ('<section-subscript>:</section-subscript>') x $nd) 
                                        . '</section-subscript-LT>)</array-R>'));
            }

          # Add extra block dimension

          if ($r->nodeName eq 'array-R')
            {
              my ($sslt) = &F ('./section-subscript-LT', $r);
              $sslt->appendChild (&t (','));
              $sslt->appendChild (&n ('<section-subscript><lower-bound><named-E><N><n>JBLK</n></N></named-E></lower-bound></section-subscript>'));
            }
          elsif ($r->nodeName eq 'parens-R')
            {
              my ($elt) = &F ('./element-LT', $r);
              $elt->appendChild (&t (','));
              $elt->appendChild (&n ('<named-E><N><n>JBLK</n></N></named-E>'));
            }
          else
            {
              die $r;
            }

 
          if ($expr->parentNode->nodeName eq 'E-1')
            {
              $ptr2w{$N} = 1;
            }
          else
            {
              my $stmt = &Fxtran::stmt ($expr);
              if ($stmt->nodeName eq 'call-stmt')
                {
                  my $intent = &getArgumentIntent ($stmt, $expr) || 'INOUT';
                  if ($intent =~ m/IN/o)
                    {
                      $ptr2r{$N} = 1;
                    }
                  if ($intent =~ m/OUT/o)
                    {
                      $ptr2w{$N} = 1;
                    }
                }
              else
                {
                  $ptr2r{$N} = 1;
                }
            }
         

        }
      

      my ($do) = &F ('./do-construct', $par);

      $par->normalize ();

      my $enddo = $do->lastChild;
      my $indent = &Fxtran::getIndent ($enddo, 1);

      for my $ptr (sort keys (%ptr2vars))
        {
          my $var = $ptr2vars{$ptr};

          $ptr2r{$ptr} or $ptr2w{$ptr} or die "$ptr";

          my $access = $ptr2r{$ptr} && $ptr2w{$ptr} ? 'RDWR' : $ptr2r{$ptr} ? 'RDONLY' : 'WRONLY';

          my $stmt = "$ptr => GET_HOST_DATA_$access ($var)";
          $stmt = "IF (ASSOCIATED ($var)) $stmt" if ($ptr2cond{$ptr});
          $stmt = &Fxtran::fxtran (statement => $stmt);
          $par->insertBefore (&t ("\n" . (' ' x $indent)), $do);
          $par->insertBefore ($stmt, $do);
        }
      $par->insertBefore (&t ("\n" . (' ' x $indent)), $do);
 

    }



  # Add pointers associated with objects components

  my @decl;
  my @init;
  my %zdum;

  for my $ptr (sort keys (%ptr2decl))
    {
      my $decl = $ptr2decl{$ptr}->cloneNode (1);
      my ($N) = &F ('.//EN-N/N/n/text()', $decl);
      $N->setData ($ptr);
      my ($sslt) = &F ('.//EN-decl/array-spec/shape-spec-LT', $decl);
      my @ss = &F ('./shape-spec', $sslt);
      for my $ss (@ss)
        {
          $ss->replaceNode (&n ('<shape-spec>:</shape-spec>'));
        }
      $sslt->appendChild (&t (','));
      $sslt->appendChild (&n ('<shape-spec>:</shape-spec>'));
      &addPointerAttr ($decl);
      push @decl, $decl;

      push @init, &Fxtran::fxtran (statement => "$ptr => ZDUM" . (scalar (@ss) + 1));
      $zdum{scalar (@ss)+1} = 1;
    }


  &addDecl ($doc, 0, @decl);
  &addInit ($doc,    @init);

  &addDecl ($doc, 0, map { &Fxtran::fxtran (statement => "REAL (KIND=JPRB), SAVE, TARGET :: ZDUM$_(" . join (',', (1)x$_) . ")") } sort keys (%zdum));

  &manageFields ($doc, \%dims, \%arr2fld);

}

sub callOpenMPRoutines
{

# Process CALL statements outside PARALLEL sections; add a parallel-call-directive node
# - suffix procedure designator with a `_OPENMP'
# - suffix procedure designator in interface include with a `_openmp'
# - replace NPROMA array arguments by field descriptor arguments; no array section allowed

  my $doc = shift;
 
  # Argument list

  my @args = &F ('.//subroutine-stmt/dummy-arg-LT/arg-N/N/n/text()', $doc);
  my %args = map { ($_->textContent, 1) } @args;

  my @call = &F ('.//call-stmt[not(ancestor::parallel-directive)]', $doc);

  my $suf = '_OPENMP';

  my %proc;

  for my $call (@call)
    {
      next if (&F ('.//procedure-designator/named-E/R-LT', $call)); # Skip objects calling methods
      my ($proc) = &F ('./procedure-designator/named-E/N/n/text()', $call);

      next if ($proc eq 'DR_HOOK');

      my $text = $call->textContent;

      my $par = &n ('<parallel-call-directive/>');
      $call->replaceNode ($par);
      $par->appendChild ($call);

      my @arg = &F ('./arg-spec/arg/named-E/N/n/text()', $call);

      my $found = 0;
      for my $arg (@arg)
        {

          $found++ if (grep { $_ eq $arg } @obj);

# TODO : build a symbol table and use it to check whether $arg is a NPROMA array

          my @ss = &F ('.//EN-decl[string(EN-N)="?"]/array-spec//shape-spec', $arg->textContent, $doc, 1); 

          # Is the actual argument a dummy argument of the current routine ?
          my $isArg = $args{$arg->textContent};

          if (@ss && ('YDCPG_OPTS%KLON' eq $ss[0]))
            {
              my ($expr) = &Fxtran::expr ($arg);
              die ("No array reference allowed in CALL statement:\n$text\n") if (&F ('./R-LT', $expr));
              $arg->replaceNode (&t (($isArg ? 'YD_' : 'YL_') . $arg->textContent));
              $found++;
            }
        }

      next unless ($found);

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

sub addInit
{
  my ($doc, @decl) = @_;
  
  my ($drhook_call) = &F ('.//call-stmt[string(.//procedure-designator)="DR_HOOK"]', $doc);

  for my $decl (@decl)
    {
      my $stmt = ref ($decl) ? $decl : &Fxtran::fxtran (statement => $decl);
      $drhook_call->parentNode->insertAfter ($stmt, $drhook_call);
      $drhook_call->parentNode->insertAfter (&t ("\n"), $drhook_call);
    }
}

sub addDecl
{
  my ($doc, $cr, @decl) = @_;
  my ($zhook_decl) = &F ('.//T-decl-stmt[string(.//EN-N)="ZHOOK_HANDLE"]', $doc);

  for my $decl (@decl)
    {
      my $stmt = ref ($decl) ? $decl : &Fxtran::fxtran (statement => $decl);
      $zhook_decl->parentNode->insertBefore ($stmt, $zhook_decl);
      $zhook_decl->parentNode->insertBefore (&t ("\n"), $zhook_decl);
      $zhook_decl->parentNode->insertBefore (&t ("\n"), $zhook_decl) if ($cr);
     }
}

sub manageFields
{
  my ($doc, $dims, $arr2fld) = @_;

  my @args = &F ('.//subroutine-stmt/dummy-arg-LT/arg-N', $doc, 1);
  my %args = map { ($_, 1) } @args;

  for my $name (sort keys (%$dims))
    {
      next if ($skip{$name});
      next if ($args{$arr2fld->{$name}});

      my $dim = $dims->{$name};

      my @par = &F ('.//parallel-directive[.//named-E[string(N)="?"]]|'
                  . './/parallel-call-directive[.//named-E[string(N)="?"]]', "YL_$name", "YL_$name", $doc);

      next unless (@par);

      my $par1 = $par[0]; my $par2 = $par[-1];     

      my @par2 = &F ('ancestor-or-self::*[not(ancestor-or-self::if-construct)]', $par2);

      $par2 = $par2[-1];

      my ($stmt1) = &F ('.//ANY-stmt', $par1);
      my ($stmt2) = &F ('.//ANY-stmt', $par2);

      my $ind1 = &Fxtran::getIndent ($stmt1);
      my $ind2 = &Fxtran::getIndent ($stmt2);

      if ($par1->firstChild->nodeName =~ m/-stmt$/o)
        {
          $par1->insertBefore (&t ("\n" . (' ' x $ind1)), $par1->firstChild) for (1 .. 2);
        }

      my @dd;


      for my $d (@{ $dim })
        {
          my @b = map { $_->textContent } @$d;
          if (@b == 2)
            {
              if ($b[0] eq '1')
                {
                  push @dd, $b[1];
                }
              elsif ($b[0] eq '0')
                {
                  push @dd, $b[1] . '+1';
                }
              else
                {
                  push @dd, $b[1] . '-' . $b[0] . '+1';
                }
            }
          else
            {
              push @dd, @b;
            }
        }

      my $args = '';

      $args .= ", NLEV=$dd[1]"  if (@dd > 2);
      $args .= ", NDIM=$dd[2]"  if (@dd > 3);
      $args .= ", NDIM2=$dd[3]" if (@dd > 4);


      $par1->insertBefore (&Fxtran::fxtran (statement => "IF (.NOT. ASSOCIATED (YL_$name)) YL_$name => CREATE_TEMPORARY (YDGEOMETRY, PERSISTENT=.TRUE.$args)"), $par1->firstChild);
      $par1->insertBefore (&t ("\n" . (' ' x $ind1)), $par1->firstChild);

      $par2->insertAfter (&t ("\n" . (' ' x $ind2)), $par2->lastChild) unless ($par2->lastChild->nodeName eq 'deallocate-stmt');

      $par2->insertAfter (&t ("\n" . (' ' x $ind2)), $par2->lastChild);
      $par2->insertAfter (&Fxtran::fxtran (statement => "IF (ASSOCIATED (YL_$name)) CALL DELETE_TEMPORARY (YL_$name)"), $par2->lastChild);
      $par2->insertAfter (&t ("\n" . (' ' x $ind2)), $par2->lastChild);
      $par2->insertAfter (&Fxtran::fxtran (statement => "YL_$name => NULL ()"), $par2->lastChild);
 
    }

}

sub cleanInterfaces
{

# Remove unused included interfaces (.intfb.h)

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

my $F90 = shift;

my $doc = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);


&singleDecl ($doc);

&renameSubroutine ($doc);

&parseDirectives ($doc);

&processDirectives ($doc);

&callOpenMPRoutines ($doc);

&makeParallel ($doc);

&cleanInterfaces ($doc);

&cleanParallelDirectives ($doc);

&addDecl ($doc, 1, 
          'INTEGER(KIND=JPIM) :: JBLK',
          'TYPE(CPG_BNDS_TYPE) :: YLCPG_BNDS');

my ($parkind1) = &F ('.//use-stmt[string(module-N)="PARKIND1"]', $doc);

for my $mod (qw (FIELD_MODULE FIELD_REGISTRY_MOD))
  {
    next if (&F ('.//use-stmt[string(module-N)="?"]', $mod, $doc));
    $parkind1->parentNode->insertAfter (&n ("<use-stmt>USE <module-N><N><n>$mod</n></N></module-N></use-stmt>"), $parkind1);
    $parkind1->parentNode->insertAfter (&t ("\n"), $parkind1);
  }

$F90 =~ s/.F90$/_openmp.F90/o;

'FileHandle'->new (">$F90")->print ($doc->textContent);

