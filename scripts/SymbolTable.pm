package SymbolTable;

use strict;
use Fxtran;

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

  for my $attr (@attr)
    {
      next unless (my ($x) = &F ('.//attribute-N[string(.)="?"]', $attr, $stmt));
      $x->parentNode->previousSibling->unbindNode ();
      $x->parentNode->unbindNode ();
    }

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

1;
