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
use SymbolTable;
use Associate;

my @obj = qw (YDMF_PHYS_BASE_STATE YDMF_PHYS_NEXT_STATE YDCPG_MISC YDCPG_PHY9
              YDCPG_PHY0 YDMF_PHYS YDCPG_DYN9 YDCPG_DYN0 YDMF_PHYS_SURF YDVARS);
my @skip = qw (PGFL PGFLT1 PGMVT1 PGPSDT2D);
my %skip = map { ($_, 1) } @skip;



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

      my $e = &n ("<$tag-section/>");
 
      for my $node (@node)
        {
          $e->appendChild ($node);
        }

      $C->replaceNode ($e);

    }
}

sub fielFifyDecl
{
  my ($doc, $t) = @_;

# First step : process all NPROMA arrays declarations
# - local arrays are added an extra dimension for blocks
# - NPROMA arrays are complemented by a field object
# - NPROMA arguments arrays are replaced by field objects

  for my $N (keys (%$t))
    {
      next if ($skip{$N});
      my $s = $t->{$N};
  
      next unless (my $nd = $s->{nd});
  
      my @ss = &F ('./shape-spec-LT/shape-spec',, $s->{as});
  
      next unless ($ss[0]->textContent eq 'YDCPG_OPTS%KLON');
  
      my $en_decl = delete $s->{en_decl};
  
      my $stmt = &Fxtran::stmt ($en_decl);
      my ($sslt) = &F ('./array-spec/shape-spec-LT', $en_decl);
      
      # Use implicit shape, with an extra dimension for blocks
  
      for ($sslt->childNodes ())
        {
          $_->unbindNode ();
        }
      for my $i (0 .. $#ss+1)
        {
          $sslt->appendChild (&n ('<shape-spec>:</shape-spec>'));
          $sslt->appendChild (&t (',')) if ($i <= $#ss);
        }
  
      &SymbolTable::addAttributes ($stmt, qw (CONTIGUOUS POINTER));
  
      my $type_fld = &SymbolTable::getFieldType ($nd, $s->{ts});
      $type_fld or die "Unknown type : " . $s->{ts}->textContent;
  
      my $decl_fld;
  
      if ($s->{arg})
        {
          my ($intent) = &SymbolTable::removeAttributes ('INTENT');
          &SymbolTable::removeAttributes ($stmt, 'INTENT');
          $s->{arg}->setData ("YD_$N");
          ($decl_fld) = &Fxtran::fxtran (statement => "TYPE ($type_fld), INTENT(" . $intent->textContent . ") :: YD_$N");
          $s->{field} = "YD_$N";
        }
      else
        {
          ($decl_fld) = &Fxtran::fxtran (statement => "TYPE ($type_fld), POINTER :: YL_$N");
          $s->{field} = "YL_$N";
        }
  
      $stmt->parentNode->insertBefore ($decl_fld, $stmt);
      $stmt->parentNode->insertBefore (&t ("\n"), $stmt);
  
  
    }

}


my $F90 = shift;

my $doc = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 300)]);

# Prepare the code

&Associate::resolveAssociates ($doc);
&Decl::forceSingleDecl ($doc);

&parseDirectives ($doc);

# Add modules

&SymbolTable::useModule ($doc, qw (FIELD_MODULE FIELD_REGISTRY_MOD));

my $t = &SymbolTable::getSymbolTable ($doc);

&fielFifyDecl ($doc, $t);

my @par = &F ('.//parallel-section', $doc);



for (&F ('.//parallel-section', $doc), &F ('.//call-stmt[not(string(procedure-designator)="DR_HOOK")]', $doc), &F ('.//skip-section', $doc), &F ('.//include', $doc))
  {
    $_->unbindNode ();
  }

&SymbolTable::renameSubroutine ($doc, sub { return $_[0] . '_OPENMP' });

$F90 =~ s/.F90$/_openmp.F90/o;

'FileHandle'->new (">$F90")->print ($doc->textContent);


