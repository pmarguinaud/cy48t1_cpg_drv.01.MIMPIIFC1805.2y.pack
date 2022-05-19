package Construct;

use strict;
use Fxtran;

sub changeIfStatementsInIfConstructs
{
  my $d = shift;

  # Change if statements into if constructs

  my @if_stmt = &F ('.//if-stmt', $d);

  for my $if_stmt (@if_stmt)
    {
      my $indent = &Fxtran::getIndent ($if_stmt);

      my ($action) = &F ('./action-stmt', $if_stmt);
      my ($stmt) = &F ('./ANY-stmt', $action);

      next if ($stmt->nodeName eq 'call-stmt' && $stmt->textContent =~ m/DR_HOOK/o);

      $action->unbindNode ();

      $stmt = $stmt->textContent;
      my $if = $if_stmt->textContent;

      my ($if_construct) = &Fxtran::fxtran (fragment => << "EOF");
$if THEN
  $stmt
ENDIF
EOF

     &Fxtran::reIndent ($if_construct, $indent);

     $if_stmt->replaceNode ($if_construct);
      
    }

}

1;
