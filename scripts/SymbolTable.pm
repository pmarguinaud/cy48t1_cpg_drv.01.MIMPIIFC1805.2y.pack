package SymbolTable;

use strict;
use Fxtran;

my @object = qw (YDMF_PHYS_BASE_STATE YDMF_PHYS_NEXT_STATE YDCPG_MISC YDCPG_PHY9
                 YDCPG_PHY0 YDMF_PHYS YDCPG_DYN9 YDCPG_DYN0 YDMF_PHYS_SURF YDVARS);
my %object = map { ($_, 1) } @object;
my @skip = qw (PGFL PGFLT1 PGMVT1 PGPSDT2D);
my %skip = map { ($_, 1) } @skip;

sub getSymbolTable
{
  my $doc = shift;

  my @args = &F ('.//subroutine-stmt/dummy-arg-LT/arg-N/N/n/text()', $doc);
  my %args = map { ($_->textContent, $_) } @args;

  my @en_decl = &F ('.//EN-decl', $doc);

  my %t;

  for my $en_decl (@en_decl)
    {
      my ($N) = &F ('.//EN-N', $en_decl, 1);
      my ($stmt) = &Fxtran::stmt ($en_decl);
      my ($ts) = &F ('./_T-spec_/*', $stmt);
      my ($as) = &F ('./array-spec', $en_decl);
      my @ss = $as ? &F ('./shape-spec-LT/shape-spec', $as) : ();
      my $nd = scalar (@ss);
      $t{$N} = {
                 object => $object{$N},
                 skip => $skip{$N},
                 nproma => $as && $ss[0]->textContent eq 'YDCPG_OPTS%KLON',
                 arg => $args{$N} || 0, 
                 ts => $ts->cloneNode (1), 
                 as => $as ? $as->cloneNode (1) : undef, 
                 nd => $nd,
                 en_decl => $en_decl,
               };
    }

  return \%t;
}

sub getAttributes
{
  my ($stmt, @attr) = @_;

  my @v;

  for my $attr (@attr)
    {
      my ($v) = &F ('.//attribute-N[string(.)="?"]', $attr, $stmt);
      push @v, $v;
    }

  return @v;
}

sub addAttributes
{
  my ($stmt, @attr) = @_;
  my $ts = $stmt->firstChild;

  for my $attr (@attr)
    {
      $stmt->insertAfter (&n ("<attribute><attribute-N>$attr</attribute-N></attribute>"), $ts);
      $stmt->insertAfter (&t (', '), $ts);
    }
}

sub removeAttributes
{
  my ($stmt, @attr) = @_;

  my @v;

  for my $attr (@attr)
    {
      next unless (my ($x) = &F ('.//attribute-N[string(.)="?"]', $attr, $stmt));
      push @v, $x->parentNode;
      $x->parentNode->previousSibling->unbindNode ();
      $x->parentNode->unbindNode ();
    }

  return @v;
}

sub getFieldType
{
  my ($nd, $ts) = @_;

  $nd++;

  ($ts = $ts->textContent) =~ s/\s+//go;

  my %ts = ('INTEGER(KIND=JPIM)' => 'INT', 'REAL(KIND=JPRB)' => '');

  return unless (defined ($ts{$ts}));

  return "FIELD_$ts{$ts}${nd}D";
}

sub useModule
{
  my ($doc, @mod) = @_;

  my ($implicit_none) = &F ('.//implicit-none-stmt', $doc);

  for my $mod (@mod)
    {
      next if (&F ('.//use-stmt[string(module-N)="?"]', $mod, $doc));
      $implicit_none->parentNode->insertBefore (&n ("<use-stmt>USE <module-N><N><n>$mod</n></N></module-N></use-stmt>"), $implicit_none);
      $implicit_none->parentNode->insertBefore (&t ("\n"), $implicit_none);
    }
}

sub renameSubroutine
{
  my ($doc, $sub) = @_;

  my @name = &F ('.//subroutine-N/N/n/text()', $doc);
  my $name = $name[0]->textContent;

  my $name1 = $sub->($name);

  for (@name)
    {
      $_->setData ($name1);
    }

  my @drhook = &F ('.//call-stmt[string(procedure-designator)="DR_HOOK"]', $doc);

  for my $drhook (@drhook)
    {
      next unless (my ($S) = &F ('./arg-spec/arg/string-E/S/text()', $drhook));
      my $str = $S->textContent;
      $str =~ s/$name/$name1/;
      $S->setData ($str);
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


1;
