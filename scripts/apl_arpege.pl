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
use Decl;
use Object;
use Loop;
use SymbolTable;
use Associate;

sub updateFile
{
  my ($F90, $code) = @_;

  my $c = do { local $/ = undef; my $fh = 'FileHandle'->new ("<$F90"); $fh ? <$fh> : undef };
  
  if ((! defined ($c)) || ($c ne $code))
    {
      'FileHandle'->new (">$F90")->print ($code);
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

      ($bdir, my $opts) = ($bdir =~ m/^(\w+)\s*(?:\(\s*(\S.*\S)\s*\)\s*)?$/goms);
      my %opts = $opts ? split (m/\s*[=,]\s*/o, $opts) : ();

      $opts{TARGET} ||= 'OPENMP';

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

      my $e = &n ("<$tag-section " . join (' ', map { sprintf ('%s="%s"', lc ($_), $opts{$_}) } keys (%opts))  . "/>");
 
      for my $node (@node)
        {
          $e->appendChild ($node);
        }

      $C->replaceNode ($e);

    }
}

sub fieldifyDecl
{
  my ($doc, $t) = @_;

# First step : process all NPROMA arrays declarations
# - local arrays are added an extra dimension for blocks
# - NPROMA arrays are complemented by a field object
# - NPROMA arguments arrays are replaced by field objects

  for my $N (sort keys (%$t))
    {
      my $s = $t->{$N};
      next if ($s->{skip});
  
      next unless ($s->{nproma});
  
      my $en_decl = delete $s->{en_decl};
  
      my $stmt = &Fxtran::stmt ($en_decl);
      my ($sslt) = &F ('./array-spec/shape-spec-LT', $en_decl);
      
      # Use implicit shape, with an extra dimension for blocks
  
      for ($sslt->childNodes ())
        {
          $_->unbindNode ();
        }
      for my $i (0 .. $s->{nd})
        {
          $sslt->appendChild (&n ('<shape-spec>:</shape-spec>'));
          $sslt->appendChild (&t (',')) if ($i < $s->{nd});
        }
  
      &SymbolTable::addAttributes ($stmt, qw (POINTER));
  
      my $type_fld = &SymbolTable::getFieldType ($s->{nd}, $s->{ts});
      $type_fld or die "Unknown type : " . $s->{ts}->textContent;
  
      my $decl_fld;
  
      if ($s->{arg})
        {
          &SymbolTable::removeAttributes ($stmt, 'INTENT');
          $s->{arg}->setData ("YD_$N");
          ($decl_fld) = &Fxtran::fxtran (statement => "TYPE ($type_fld), POINTER :: YD_$N");
          $s->{field} = &n ("<named-E><N><n>YD_$N</n></N></named-E>");
        }
      else
        {
          ($decl_fld) = &Fxtran::fxtran (statement => "TYPE ($type_fld), POINTER :: YL_$N");
          $s->{field} = &n ("<named-E><N><n>YL_$N</n></N></named-E>");
        }
  
      $stmt->parentNode->insertBefore ($decl_fld, $stmt);
      $stmt->parentNode->insertBefore (&t ("\n"), $stmt);
  
  
    }

}

sub makeParallel
{
  my ($par, $t, $hook_suffix) = @_;

  # Add a loop nest on blocks

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

  my %intent;

  # Process each expression (only NPROMA local arrays and FIELD API backed data) in the parallel section

  for my $expr (&F ('.//named-E', $par))
    {
      my ($N) = &F ('./N', $expr, 1);
      my $s = $t->{$N};
      next if ($s->{skip});

      # Object wrapping fields : replace by a pointer to data with all blocks
      if ($s->{object})
        {
          my ($name) = &F ('./N/n/text()', $expr); 
          my @ctl = &F ('./R-LT/component-R/ct', $expr, 1);
          my $e = $expr->cloneNode (1);
          my @r = &F ('./R-LT/component-R', $expr);
          $_->unbindNode for (@r);
          my $ptr = join ('_', 'Z', $N, @ctl);
          $name->setData ($ptr);

          # Create new entry in symbol table 
          # we record the pointer wich will be used to access the object component
          unless ($t->{$ptr})
            {
              my $key = join ('%', &Object::getObjectType ($s, $N), @ctl);
              my $decl;
              eval
                {
                  $decl = &Object::getObjectDecl ($key);
                };
              if (my $c = $@)
                {
                  my $stmt = &Fxtran::stmt ($expr); 
                  die $c . $stmt->textContent;
                }
              my ($as) = &F ('.//array-spec', $decl);
              my ($ts) = &F ('./_T-spec_/*', $decl);
              my @ss = &F ('./shape-spec-LT/shape-spec', $as);
              my $nd = scalar (@ss);
              $t->{$ptr} = {
                             object => 0,
                             skip => 0,
                             nproma => 1,
                             arg => 0,
                             ts => $ts,
                             as => $as,
                             nd => $nd,
                             field => &Object::getFieldFromObjectComponents ($N, @ctl),
                             object_based => 1, # postpone pointer declaration
                           };
            }
          $N = $ptr;
          $s = $t->{$ptr};
        }
      # Local NPROMA array
      elsif ($s->{nproma})
        {
        }
      # Other: skip
      else
        {
          next;
        }

      &addExtraIndex ($expr, &n ("<named-E><N><n>JBLK</n></N></named-E>"), $s);

      &SymbolTable::grokIntent ($expr, \$intent{$N});
    }

  my %intent2access = qw (IN RDONLY INOUT RDWR OUT WRONLY);

  $par->insertBefore (&t ("\n" . (' ' x $indent)), $loop);

  $par->insertBefore (&s ("IF (LHOOK) CALL DR_HOOK ('GET_DATA$hook_suffix',0,ZHOOK_HANDLE_FIELD_API)"), $loop);
  $par->insertBefore (&t ("\n" . (' ' x $indent)), $loop);

  $par->insertAfter (&s ("IF (LHOOK) CALL DR_HOOK ('NULLIFY$hook_suffix',1,ZHOOK_HANDLE_FIELD_API)"), $loop);
  $par->insertAfter (&t ("\n" . (' ' x $indent)), $loop);

  for my $ptr (reverse (sort keys (%intent)))
    {
      my $s = $t->{$ptr};
      my $access = $intent2access{$intent{$ptr}};
      my $var = $s->{field}->textContent;
      my $stmt = &Fxtran::fxtran (statement => "$ptr => GET_HOST_DATA_$access ($var)");
      $par->insertBefore ($stmt, $loop);
      $par->insertBefore (&t ("\n" . (' ' x $indent)), $loop);

      $par->insertAfter (&s ("$ptr => NULL ()"), $loop);
      $par->insertAfter (&t ("\n" . (' ' x $indent)), $loop);
    }
  $par->insertBefore (&s ("IF (LHOOK) CALL DR_HOOK ('GET_DATA$hook_suffix',1,ZHOOK_HANDLE_FIELD_API)"), $loop);
  $par->insertBefore (&t ("\n" . (' ' x $indent)), $loop);

  $par->insertAfter (&s ("IF (LHOOK) CALL DR_HOOK ('NULLIFY$hook_suffix',0,ZHOOK_HANDLE_FIELD_API)"), $loop);
  $par->insertAfter (&t ("\n" . (' ' x $indent)), $loop);

  $par->insertBefore (&t ("\n" . (' ' x $indent)), $loop);
  $par->insertAfter (&t ("\n" . (' ' x $indent)), $loop);


}

sub addExtraIndex
{
  my ($expr, $ind, $s) = @_;

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
      $rlt->appendChild ($r = &n ('<array-R>(<section-subscript-LT>' . join (',', ('<section-subscript>:</section-subscript>') x $s->{nd}) 
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


}

sub callParallelRoutine
{
# Process CALL statement outside PARALLEL sections
# Replace NPROMA array arguments by field descriptor arguments; no array section allowed

  my ($call, $t) = @_;

  my $text = $call->textContent;

  my @arg = &F ('./arg-spec/arg/named-E/N/n/text()', $call);

  my $found = 0;
  for my $arg (@arg)
    {
      my $s = $t->{$arg};
      next if ($s->{skip});
      $found++ if ($s->{object});
   

      # Is the actual argument a dummy argument of the current routine ?
      my $isArg = $s->{arg};


      if ($s->{nproma})
        {
          my ($expr) = &Fxtran::expr ($arg);
          die ("No array reference allowed in CALL statement:\n$text\n") if (&F ('./R-LT', $expr));
          $s->{field} or die &Dumper ([$arg->textContent, $s, $text]);
          $expr->replaceNode ($s->{field}->cloneNode (1));
          $found++;
        }
      elsif ($s->{object})
        {
          my ($expr) = &Fxtran::expr ($arg);
          my @ctl = &F ('./R-LT/component-R/ct', $expr, 1);
          if (@ctl)
            {
              my $e = &Object::getFieldFromObjectComponents ($arg->textContent, @ctl);
              $expr->replaceNode ($e);
            }
        }
    }

  return $found;
}

sub removeUnusedIncludes
{
  my $doc = shift;
  for my $include (&F ('.//include', $doc))
    {
      my ($filename) = &F ('./filename', $include, 2);
      (my $name = $filename) =~ s/\.intfb.h$//o;
      $name = uc ($name);
      next if (&F ('.//call-stmt[string(procedure-designator)="?"]', $name, $doc));
      my $next = $include->nextSibling;
      $next->unbindNode () if ($next->textContent eq "\n");
      $include->unbindNode ();
    }
}

sub setupLocalFields
{
# Use CREATE_TEMPORARY at the beginning of the routine (after the call to DR_HOOK) to create
# FIELD API objects backing local NPROMA arrays
# Use DELETE_TEMPORARY to delete these FIELD API objects

  my ($doc, $t, $hook_suffix) = @_;

  my @drhook = &F ('.//if-stmt[.//call-stmt[string(.//procedure-designator)="DR_HOOK"]]', $doc);
  @drhook = @drhook[0,-1];

  my ($drhook1, $drhook2) = @drhook;
  my ($ind1, $ind2) = map { &Fxtran::getIndent ($_) } ($drhook1, $drhook2);
  my ($p1  , $p2  ) = map { $_->parentNode          } ($drhook1, $drhook2);

  $p1->insertAfter (&s ("IF (LHOOK) CALL DR_HOOK ('CREATE_TEMPORARIES$hook_suffix',1,ZHOOK_HANDLE_FIELD_API)"), $drhook1);
  $p1->insertAfter (&t ("\n" . (' ' x $ind1)), $drhook1);

  $p2->insertBefore (&s ("IF (LHOOK) CALL DR_HOOK ('DELETE_TEMPORARIES$hook_suffix',0,ZHOOK_HANDLE_FIELD_API)"), $drhook2);
  $p2->insertBefore (&t ("\n" . (' ' x $ind2)), $drhook2);

  for my $n (sort keys (%$t))
    {
      my $s = $t->{$n};

      next unless ($s->{nproma});
      next if ($s->{object_based} || $s->{arg});
      my @ss = &F ('./shape-spec-LT/shape-spec', $s->{as});

      my (@lb, @ub);
    
      for my $i (0 .. $#ss)
        {
          my $ss = $ss[$i];
          my @b = map { $_->textContent } &F ('./ANY-bound', $ss);
          unshift (@b, '1') if (@b == 1);
          push @lb, $b[0];
          push @ub, $b[1];
        }

      push @lb, '1';
      push @ub, 'YDCPG_OPTS%KGPBLKS';
      
      my $f = $s->{field}->textContent;

      my $ubounds = 'UBOUNDS=[' . join (', ', @ub) . '], ';
      my $lbounds = grep ({ $_ ne '1' } @lb) 
                  ? 'LBOUNDS=[' . join (', ', @lb) . '], '
                  : '';

      $p1->insertAfter (&s ("CALL CREATE_TEMPORARY_LU ($f, ${ubounds}${lbounds}PERSISTENT=.TRUE.)"), $drhook1);
      $p1->insertAfter (&t ("\n" . (' ' x $ind1)), $drhook1);

      $p2->insertBefore (&s ("IF (ASSOCIATED ($f)) CALL DELETE_TEMPORARY ($f)"), $drhook2);
      $p2->insertBefore (&t ("\n" . (' ' x $ind2)), $drhook2);
      

    }

  $p1->insertAfter (&s ("IF (LHOOK) CALL DR_HOOK ('CREATE_TEMPORARIES$hook_suffix',0,ZHOOK_HANDLE_FIELD_API)"), $drhook1);
  $p1->insertAfter (&t ("\n" . (' ' x $ind1)), $drhook1);

  $p2->insertBefore (&s ("IF (LHOOK) CALL DR_HOOK ('DELETE_TEMPORARIES$hook_suffix',1,ZHOOK_HANDLE_FIELD_API)"), $drhook2);
  $p2->insertBefore (&t ("\n" . (' ' x $ind2)), $drhook2);


}

sub setOpenMPDirective
{
  my ($par, $t) = @_;

  my @N = &F ('.//named-E/N', $par);

  my %priv;

  for my $N (@N)
    {
      my $n = $N->textContent;
      my $s = $t->{$n};
      next if ($s->{nproma});
      my $expr = $N->parentNode;
      my $p = $expr->parentNode;
      
      $priv{$n}++ if (($p->nodeName eq 'E-1') || ($p->nodeName eq 'do-V'));
    }

  my ($do) = &F ('.//do-construct[./do-stmt[string(do-V)="JBLK"]]', $par);

  my $indent = &Fxtran::getIndent ($do);

  my $C = &n ('<C>!$OMP PARALLEL DO PRIVATE (' . join (', ', sort keys (%priv))  . ')</C>');
  
  $do->parentNode->insertBefore ($C, $do);
  $do->parentNode->insertBefore (&t ("\n" . (' ' x $indent)), $do);

}


my $suffix = '_parallel';

my $F90 = shift;

my $doc = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);

# Prepare the code

&Associate::resolveAssociates ($doc);
&Decl::forceSingleDecl ($doc);

&parseDirectives ($doc);

# Add modules

&SymbolTable::useModule ($doc, qw (FIELD_MODULE FIELD_REGISTRY_MOD FIELD_HELPER_MODULE));

# Add local variables

&SymbolTable::addDecl ($doc, 1, 
          'INTEGER(KIND=JPIM) :: JBLK',
          'TYPE(CPG_BNDS_TYPE) :: YLCPG_BNDS', 
          'REAL(KIND=JPRB) :: ZHOOK_HANDLE_FIELD_API');

my $t = &SymbolTable::getSymbolTable ($doc);

for my $v (qw (JLON JLEV))
  {
    &SymbolTable::addDecl ($doc, 1, "INTEGER(KIND=JPIM) :: $v") unless ($t->{$v});
  }

# Remove SKIP sections

for (&F ('.//skip-section', $doc))
  {
    $_->unbindNode ();
  }

# Transform NPROMA fields into a pair of (FIELD API object, Fortran pointer)

&fieldifyDecl ($doc, $t);

# Process parallel sections

my @par = &F ('.//parallel-section', $doc);


for my $par (@par)
  {
    &makeParallel ($par, $t, '');
  }

# Process call to parallel routines

my @call = &F ('.//call-stmt[not(ancestor::parallel-section)]' # Skip calls in parallel sections
            . '[not(string(procedure-designator)="DR_HOOK")]'  # Skip DR_HOOK calls
            . '[not(procedure-designator/named-E/R-LT)]'       # Skip objects calling methods
            . '[not(ancestor::serial-section)]', $doc);        # Skip calls in serial sections

my %seen;

for my $call (@call)
  {
    if (&callParallelRoutine ($call, $t))
      {
        # Add include for the parallel CALL
        my ($name) = &F ('./procedure-designator/named-E/N/n/text()', $call);
        unless ($seen{$name->textContent}++)
          {
            my ($include) = &F ('.//include[./filename[string(.)="?"]]', lc ($name) . '.intfb.h', $doc);
            $include->parentNode->insertAfter (&n ('<include>#include "<filename>' . lc ($name) . '_parallel.intfb.h</filename>"</include>'), $include);
            $include->parentNode->insertAfter (&t ("\n"), $include);
          }
        $name->setData ($name->data . uc ($suffix));
      }
  }

# Declare pointers required for parallel sections

my @decl;

for my $n (sort keys (%$t))
  {
    my $s = $t->{$n};
    next unless ($s->{object_based});
    my $decl = &Fxtran::fxtran (statement => $s->{ts}->textContent . ", POINTER :: " . $n . "(" . join (',', (':') x ($s->{nd} + 1)) . ")");
    push @decl, $decl;
  }


&SymbolTable::addDecl ($doc, 0, @decl);

# Create/delete fields for local arrays

&setupLocalFields ($doc, $t, '');

&removeUnusedIncludes ($doc);


for my $par (@par)
  {
    my $target = $par->getAttribute ('target');
    if ($target eq 'OPENACC')
      {
        use SingleDirective;
        &Loop::arraySyntaxLoop ($par, $t);
        &SingleDirective::hoistJlonLoops ($par);
      }
    elsif ($target eq 'OPENMP')
      {
        &setOpenMPDirective ($par, $t);
      }
  }


&SymbolTable::renameSubroutine ($doc, sub { return $_[0] . uc ($suffix) });

$F90 =~ s/.F90$/$suffix.F90/o;

&updateFile ($F90, $doc->textContent);



