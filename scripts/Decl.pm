package Decl;

use strict;
use Fxtran;

sub forceSingleDecl
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

1;
