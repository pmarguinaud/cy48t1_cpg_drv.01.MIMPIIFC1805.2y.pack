#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

#my ($var, $typ, $mod) = qw (YDVARS FIELD_VARIABLES FIELD_VARIABLES_MOD);
my ($var, $typ, $mod) = qw (YDSURFVARS SURFACE_VARIABLES SURFACE_VARIABLES_MOD);

my @r = qw (APLPAR APL_AROME MF_PHYS CPG CPG_GP CPG_DYN CPG_DIA);

my $F90 = shift;

my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 200)]);

# Dummy argument

{ 
  my ($ydgeometry_dummy) = &F ('//dummy-arg-LT/arg-N[string (.)="YDGEOMETRY"]', $d);
  my $comma = $ydgeometry_dummy->nextSibling ();
  my $lst = $ydgeometry_dummy->parentNode ();

  $lst->insertAfter (&n ("<arg-N><N><n>$var</n></N></arg-N>"), $ydgeometry_dummy);
  $lst->insertAfter ($comma->cloneNode (1), $ydgeometry_dummy);
}

# CALL statements

for my $r (@r)
  {
    my @call = &F ('//call-stmt[string (procedure-designator)="?"]', $r, $d);

    for my $call (@call)
      {
        my ($ydgeometry_actual) = &F ('.//arg[string (.)="YDGEOMETRY"]', $call);
        my $comma = $ydgeometry_actual->nextSibling ();
        my $lst = $ydgeometry_actual->parentNode ();

        $lst->insertAfter (&n ("<arg><named-E><N><n>$var</n></N></named-E></arg>"), $ydgeometry_actual);
        $lst->insertAfter ($comma->cloneNode (1), $ydgeometry_actual);
      }
  }

# Declaration 

{
  my ($ydgeometry_decl) = &F ('//T-decl-stmt[.//EN-N[string (.)="YDGEOMETRY"]', $d);
  my ($cr) = &f ('following::text ()[contains (.,"' . "\n" . '")]', $ydgeometry_decl);
  
  my $decl = $ydgeometry_decl->cloneNode (1);
  
  my ($N) = &F ('.//EN-N/N/n/text()', $decl);
  $N->setData ("$var");
  
  my ($intent) = &F ('.//attribute[string (attribute-N)="INTENT"]/intent-spec/text()', $decl);
  $intent->setData ('INOUT');
  
  my $attribute = $intent->parentNode->parentNode;
  my $space1 = $attribute->nextSibling;
  
  if ($space1->nodeName eq '#text')
    {
      my $text = $space1->data;
      $text =~ s/^ //o for (1 .. 3);
      $space1->setData ($text);
    }
  
  my ($spec) = &F ('./_T-spec_', $decl);
  my ($type) = &F ('./derived-T-spec/T-N/N/n/text()', $spec);
  
  $type->setData ($typ);
  
  my $space2 = $spec->nextSibling;
  
  if ($space2->nodeName eq '#text')
    {
      my $text = $space2->data;
      for (1 .. length ($typ) - length ("GEOMETRY"))
        {
          ($text =~ s/ $//o) or ($text =~ s/^ //o);
        }
      $space2->setData ($text);
    }
  
  
  $cr->parentNode->insertAfter (&t ("\n"), $cr);
  $cr->parentNode->insertAfter ($decl, $cr);
}

# Use statement

{
  my ($use_geometry) = &F ('.//use-stmt[string (module-N)="GEOMETRY_MOD"]', $d);
  my $use_field_variables = $use_geometry->cloneNode (1);

  my ($moduleN) = &F ('.//module-N', $use_field_variables);
  my ($N) = &F ('./N/n/text()', $moduleN);
  $N->setData ($mod);
  my $space = $moduleN->nextSibling;
  my $text = $space->data;
  for (1 .. length ($typ) - length ("GEOMETRY"))
    {
      ($text =~ s/^ //o) or ($text =~ s/ ONLY/ONLY/o);
    }

  $space->setData ($text);

  my ($renameLT) = &F ('.//rename-LT', $use_field_variables);
  for ($renameLT->childNodes)
    {
      $_->unbindNode ();
    }
  $renameLT->appendChild (&n ("<rename><use-N><N><n>$typ</n></N></use-N></rename>"));

  my ($cr) = &f ('following::text ()[contains (.,"' . "\n" . '")]', $use_geometry);
  
  $cr->parentNode->insertAfter ($use_field_variables, $use_geometry);
  $cr->parentNode->insertAfter (&t ("\n"), $use_geometry);

}




'FileHandle'->new (">$F90.new")->print ($d->textContent);






