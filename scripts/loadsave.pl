#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use lib $Bin;
use Fxtran;
use Getopt::Long;
use File::Path;


sub process_decl
{
  my ($opts, $en_decl, $sname, $prefix, $BODY_SAVE, $BODY_LOAD, $BODY_COPY, $BODY_SIZE, $U, $J, $L, $B) = @_;

  my (@BODY_SAVE, @BODY_LOAD, @BODY_COPY, @BODY_SIZE);
  my (%U, %J, %L, %B);

  my $stmt = &Fxtran::stmt ($en_decl);

  return if ($stmt->nodeName eq 'final-stmt');

  my %attr = map { ($_, 1) } &f ('.//f:attribute/f:attribute-N/text ()', $stmt);

  return if ($attr{PARAMETER});

  my ($name) = &f ('.//f:EN-N/f:N/f:n/text ()', $en_decl, 1);
  
  my $skip = grep { "$sname$name" eq $_ } @{ $opts->{skip} };

  if ($skip)
    {
      if ($attr{POINTER})
        {
          push @BODY_LOAD, "NULLIFY ($prefix$name)";
        }
      goto RETURN;
    }
  

  my ($tspec) = &f ('./f:_T-spec_', $stmt);

  die $stmt->toString unless ($tspec);

  return if (uc ($tspec->textContent) eq 'PROCEDURE');


  my ($intrinsic) = &f ('./f:intrinsic-T-spec', $tspec);
  my ($tname) = &f ('./f:derived-T-spec/f:T-N/f:N/f:n/text ()', $tspec);
  
  $tname && ($U{$tname} = 1);
  
  
  my @ss = &f ('./f:array-spec/f:shape-spec-LT/f:shape-spec', $en_decl, 1);

  unless (@ss)
    {
      my $stmt = &Fxtran::stmt ($en_decl);
      @ss = &F ('.//attribute/array-spec/shape-spec-LT/shape-spec', $stmt, 1);
    }
  
  
  if ($attr{POINTER} || $attr{ALLOCATABLE})
    {
      my $func = $attr{POINTER} ? 'ASSOCIATED' : 'ALLOCATED';
      push @BODY_SAVE, "L$name = $func ($prefix$name)\n";
      push @BODY_COPY, "L$name = $func ($prefix$name)\n";
      push @BODY_SIZE, "L$name = $func ($prefix$name)\n";
      push @BODY_SAVE, "WRITE (KLUN) L$name\n";
      push @BODY_LOAD, "READ (KLUN) L$name\n";
      $L{$name} = 1;
      push @BODY_SAVE, "IF (L$name) THEN\n";
      push @BODY_LOAD, "IF (L$name) THEN\n";
      push @BODY_COPY, "IF (L$name) THEN\n";
      push @BODY_SIZE, "IF (L$name) THEN\n";
      if (@ss)
        {
          push @BODY_SAVE, "WRITE (KLUN) LBOUND ($prefix$name)\n";
          push @BODY_SAVE, "WRITE (KLUN) UBOUND ($prefix$name)\n";
          $B{scalar (@ss)} = 1;
          my $r = scalar (@ss);
          push @BODY_LOAD, "READ (KLUN) IL" . scalar (@ss) . "\n";
          push @BODY_LOAD, "READ (KLUN) IU" . scalar (@ss) . "\n";
          push @BODY_LOAD, "ALLOCATE ($prefix$name (" . join (', ', map { "IL$r($_):IU$r($_)" } (1 .. $#ss+1) ) . "))\n";
        }
      else
        {
          push @BODY_LOAD, "ALLOCATE ($prefix$name)\n";
        }
      push @BODY_COPY, "!\$acc enter data create ($prefix$name)\n";
      push @BODY_COPY, "!\$acc update device ($prefix$name)\n",
    }
  

  if ($intrinsic)
    {
      my ($tn) = &F ('./T-N', $intrinsic, 1);
      push @BODY_SAVE, "WRITE (KLUN) $prefix$name\n";
      push @BODY_LOAD, "READ (KLUN) $prefix$name\n";
      my $size = "ISIZE = KIND ($prefix$name)"; 
      $size .= " * SIZE ($prefix$name)" if (@ss);
      $size .= " * LEN ($prefix$name)" if ($tn eq 'CHARACTER');
      push @BODY_SIZE, $size . "\n", 
                       "IF (LDPRINT) THEN\n", 
                       "WRITE (*, '(I10,\" \")', ADVANCE='NO') ISIZE\n", 
                       "WRITE (*, *) TRIM (CDPATH)//'%$name'\n", 
                       "ENDIF\n", 
                       "KSIZE = KSIZE + ISIZE\n";
    }
  else 
    {
      push @BODY_SIZE, "JSIZE = 0\n";
      for (my $i = $#ss+1; $i >= 1; $i--)
        {
          $J{"J$i"} = 1;
          my $do = "DO J$i = LBOUND ($prefix$name, $i), UBOUND ($prefix$name, $i)\n";
          push @BODY_SAVE, $do;
          push @BODY_LOAD, $do;
          push @BODY_COPY, $do;
          push @BODY_SIZE, $do;
        }
      my @J = map { "J$_"  } (1 .. $#ss+1);
      my $J = @ss ? " (" . join (', ', @J) . ")" : '';
#     my $LDPRINT = @J ? ".FALSE." : "LDPRINT";
      my $LDPRINT = '.FALSE.';
      push @BODY_SAVE, ('  ' x scalar (@ss)) 
                   . "CALL SAVE (KLUN, $prefix$name" . $J . ")\n";
      push @BODY_LOAD, ('  ' x scalar (@ss)) 
                   . "CALL LOAD (KLUN, $prefix$name" . $J . ")\n";
      push @BODY_COPY, ('  ' x scalar (@ss)) 
                   . "CALL COPY ($prefix$name" . $J . ")\n";
      push @BODY_SIZE, ('  ' x scalar (@ss))
                   . "ISIZE = SIZE ($prefix$name" . $J . ", CDPATH//'%$name', $LDPRINT)\n", 
                     "JSIZE = JSIZE + ISIZE\n",
                     "KSIZE = KSIZE + ISIZE\n";
      for (my $i = $#ss; $i >= 0; $i--)
        {
          push @BODY_SAVE, "ENDDO\n";
          push @BODY_LOAD, "ENDDO\n";
          push @BODY_COPY, "ENDDO\n";
          push @BODY_SIZE, "ENDDO\n";
        }
      push @BODY_SIZE, 
                     "IF (LDPRINT) THEN\n", 
                     "WRITE (*, '(I10,\" \")', ADVANCE='NO') JSIZE\n",
                     "WRITE (*, *) TRIM (CDPATH)//'%$name'\n",
                     "ENDIF\n", 
    }
  
  if ($attr{POINTER} || $attr{ALLOCATABLE})
    {
      push @BODY_SAVE, "ENDIF\n";
      if ($attr{POINTER})
        {
          push @BODY_LOAD, "ELSE\n", "NULLIFY ($prefix$name)\n";
        }
      push @BODY_LOAD, "ENDIF\n";
      push @BODY_COPY, "!\$acc enter data attach ($prefix$name)\n";
      push @BODY_COPY, "ENDIF\n";
      push @BODY_SIZE, "ENDIF\n";
    }
  push @BODY_COPY, "\n";

RETURN:
  
  push @$BODY_SAVE, @BODY_SAVE;
  push @$BODY_LOAD, @BODY_LOAD;
  push @$BODY_COPY, @BODY_COPY;
  push @$BODY_SIZE, @BODY_SIZE;

  %$U = (%$U, %U); %$J = (%$J, %J); 
  %$L = (%$L, %L); %$B = (%$B, %B); 

}

sub indent
{
  my $n = 0;
  for (@_)
    {
      chomp;
      s/^\s*//o;
      $n-- if (m/^\s*(?:ELSE|ENDIF|ENDDO)\b/o);
      $_ = ('  ' x $n) . $_;
      $n++ if (m/^\s*(?:ELSE|IF|DO)\b/o);
    }
}

sub r
{
  my $f = shift;
  return '' unless (-f $f);
  return do { local $/ = undef; my $fh = 'FileHandle'->new ("<$f"); <$fh> };
}

sub w
{
  my $f = shift;
  my $t = &r ($f);
  return if ($t eq $_[0]);
  'FileHandle'->new (">$f")->print ($_[0]);
}

sub process_types
{
  my ($doc, $opts) = @_;

  my ($mod) = &f ('.//f:module-stmt/f:module-N/f:N/f:n/text ()', $doc);
  
  my @tconst = &f ('.//f:T-construct', $doc);
  
  for my $tconst (@tconst)
    {
      my ($INTERFACE_SAVE, $CONTAINS_SAVE) = ('', '');
      my ($INTERFACE_LOAD, $CONTAINS_LOAD) = ('', '');
      my ($INTERFACE_COPY, $CONTAINS_COPY) = ('', '');
      my ($INTERFACE_SIZE, $CONTAINS_SIZE) = ('', '');
  
      my ($name) = &f ('.//f:T-stmt/f:T-N/f:N/f:n/text ()', $tconst, 1);
      my $tname = $name;
  
  
      $INTERFACE_SAVE .= "MODULE PROCEDURE SAVE_$name\n";
      $INTERFACE_LOAD .= "MODULE PROCEDURE LOAD_$name\n";
      $INTERFACE_COPY .= "MODULE PROCEDURE COPY_$name\n";
      $INTERFACE_SIZE .= "MODULE PROCEDURE SIZE_$name\n";
  
      my (@BODY_SAVE, @BODY_LOAD, @BODY_COPY, @BODY_SIZE);

      push @BODY_SIZE, 'KSIZE = 0';
    
      my (%U, %J, %L, %B);
  
      my @en_decl = &f ('.//f:EN-decl', $tconst);
      for my $en_decl (@en_decl)
        {
          &process_decl ($opts, $en_decl, "$tname%", 'YD%', \@BODY_SAVE, \@BODY_LOAD, \@BODY_COPY, \@BODY_SIZE, \%U, \%J, \%L, \%B);
        }
  
      my $DECL_SAVE = '';
      my $DECL_LOAD = '';
      my $DECL_COPY = '';
      my $DECL_SIZE = "INTEGER*8 :: ISIZE, JSIZE\n";
  
      if (%J)
        {
          $DECL_SAVE .= "INTEGER :: " . join (', ', sort keys (%J)) . "\n";
          $DECL_LOAD .= "INTEGER :: " . join (', ', sort keys (%J)) . "\n";
          $DECL_COPY .= "INTEGER :: " . join (', ', sort keys (%J)) . "\n";
          $DECL_SIZE .= "INTEGER :: " . join (', ', sort keys (%J)) . "\n";
        }
      if (%B)
        {
          $DECL_LOAD .= "INTEGER :: " . join (', ', map  { ("IL$_($_)", "IU$_($_)") } sort keys (%B)) . "\n";
        }
      if (%L)
        {
          $DECL_SAVE .= "LOGICAL :: " . join (', ', map { "L$_" } sort keys (%L)) . "\n";
          $DECL_LOAD .= "LOGICAL :: " . join (', ', map { "L$_" } sort keys (%L)) . "\n";
          $DECL_COPY .= "LOGICAL :: " . join (', ', map { "L$_" } sort keys (%L)) . "\n";
          $DECL_SIZE .= "LOGICAL :: " . join (', ', map { "L$_" } sort keys (%L)) . "\n";
        }
  
      my $USE_SAVE = join ('', map { "USE UTIL_${_}_MOD\n" } grep { $_ ne $name } sort keys (%U));
      my $USE_LOAD = join ('', map { "USE UTIL_${_}_MOD\n" } grep { $_ ne $name } sort keys (%U));
      my $USE_COPY = join ('', map { "USE UTIL_${_}_MOD\n" } grep { $_ ne $name } sort keys (%U));
      my $USE_SIZE = join ('', map { "USE UTIL_${_}_MOD\n" } grep { $_ ne $name } sort keys (%U));
  
      for ($USE_SAVE, $USE_SAVE, $USE_COPY, $USE_SIZE, $DECL_SAVE, $DECL_LOAD)
        {
          chomp ($_);
        }
  
      $CONTAINS_SAVE .= << "EOF";
SUBROUTINE SAVE_$name (KLUN, YD)
$USE_SAVE
IMPLICIT NONE
INTEGER, INTENT (IN) :: KLUN
TYPE ($name), INTENT (IN) :: YD
EOF

      $CONTAINS_LOAD .= << "EOF";
SUBROUTINE LOAD_$name (KLUN, YD)
$USE_LOAD
IMPLICIT NONE
INTEGER, INTENT (IN) :: KLUN
TYPE ($name), INTENT (OUT) :: YD
EOF

      $CONTAINS_COPY .= << "EOF";
SUBROUTINE COPY_$name (YD)
$USE_COPY
IMPLICIT NONE
TYPE ($name), INTENT (IN) :: YD
EOF

      $CONTAINS_SIZE .= << "EOF";
INTEGER*8 FUNCTION SIZE_$name (YD, CDPATH, LDPRINT) RESULT (KSIZE)
$USE_SIZE
IMPLICIT NONE
TYPE ($name),     INTENT (IN) :: YD
CHARACTER(LEN=*), INTENT (IN) :: CDPATH
LOGICAL,          INTENT (IN) :: LDPRINT
EOF

      &indent (@BODY_SAVE);
      &indent (@BODY_LOAD);
      &indent (@BODY_COPY);
      &indent (@BODY_SIZE);


      $CONTAINS_SAVE .= $DECL_SAVE . "\n" . join ("\n", @BODY_SAVE, '') . "END SUBROUTINE\n";
      $CONTAINS_LOAD .= $DECL_LOAD . "\n" . join ("\n", @BODY_LOAD, '') . "END SUBROUTINE\n";
      $CONTAINS_COPY .= $DECL_COPY . "\n" . join ("\n", @BODY_COPY, '') . "END SUBROUTINE\n";
      $CONTAINS_SIZE .= $DECL_SIZE . "\n" . join ("\n", @BODY_SIZE, '') . "END FUNCTION\n";

      for ($CONTAINS_SAVE, $CONTAINS_SAVE, $CONTAINS_COPY, $CONTAINS_SIZE, $INTERFACE_SAVE, $INTERFACE_LOAD, $INTERFACE_COPY, $INTERFACE_SIZE)
        {
          chomp ($_);
        }
  
      my $n = lc ($name);

      $CONTAINS_SAVE = '' unless ($opts->{save});
      $CONTAINS_LOAD = '' unless ($opts->{load});
      $CONTAINS_COPY = '' unless ($opts->{copy});
      $CONTAINS_SIZE = '' unless ($opts->{size});

      $INTERFACE_SAVE = "INTERFACE SAVE\n$INTERFACE_SAVE\nEND INTERFACE\n";
      $INTERFACE_LOAD = "INTERFACE LOAD\n$INTERFACE_LOAD\nEND INTERFACE\n";
      $INTERFACE_COPY = "INTERFACE COPY\n$INTERFACE_COPY\nEND INTERFACE\n";
      $INTERFACE_SIZE = "INTERFACE SIZE\n$INTERFACE_SIZE\nEND INTERFACE\n";
  
      $INTERFACE_SAVE = '' unless ($opts->{save});
      $INTERFACE_LOAD = '' unless ($opts->{load});
      $INTERFACE_COPY = '' unless ($opts->{copy});
      $INTERFACE_SIZE = '' unless ($opts->{size});

      &w ("$opts->{dir}/util_${n}_mod.F90", << "EOF");
MODULE UTIL_${name}_MOD

USE $mod, ONLY : $name

$INTERFACE_SAVE
$INTERFACE_LOAD
$INTERFACE_COPY
$INTERFACE_SIZE

CONTAINS

$CONTAINS_SAVE
$CONTAINS_LOAD
$CONTAINS_COPY
$CONTAINS_SIZE

END MODULE
EOF

    }
}

my %opts;

&GetOptions
(
  'skip=s@' => \$opts{skip}, 'dir=s' => \$opts{dir},
  size => \$opts{size}, save => \$opts{save}, load => \$opts{load}, copy => \$opts{copy},
);

( -d $opts{dir}) or &mkpath ($opts{dir});

my $F90 = shift;

my $doc = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 800)]);

&process_types ($doc, \%opts);

