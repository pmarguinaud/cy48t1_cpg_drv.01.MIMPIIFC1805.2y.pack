#!/home/gmap/mrpm/marguina/install/perl-5.32.0/bin/perl -w
#
use strict;
use FindBin qw ($Bin);
use Data::Dumper;
use FileHandle;
use File::Basename;
use lib $Bin;
use Fxtran;

for my $F90 ("src/local/arpifs/phys_dmn/mf_phys.F90", "src/local/arpifs/adiab/cpg_gp.F90")
  {

    my $d = &Fxtran::fxtran (location => $F90, fopts => [qw (-line-length 200)]);
    
    my $f = &Fxtran::fxtran (location => '.vimpack/src=/.fypp/arpifs/module/field_variables_mod.F90', fopts => [qw (-line-length 200)]);
    
    my ($ft) = &F ('//T-construct[.//T-stmt[string (T-N)="FIELD_VARIABLES"]]', $f);
    my @fn = &F ('./component-decl-stmt//EN-N//text()', $ft, 1);
    
    my @call = &F ('//call-stmt[string (procedure-designator)="SC2PRG"]', $d);
    
    for my $call (@call)
      {
        my ($argspec) = &F ('./arg-spec', $call);
        my @arg = &F ('./arg/ANY-E', $argspec);
        next unless (scalar (@arg) == 3);
    
        my $stmt;
    
        if ($arg[1]->textContent =~ m/^(?:PSP_(?:RR|SB|SG)|PSD_(?:VF|VA|VD|VH|VK|VP|VV))$/o)
          {
            my ($p, $x, $v) = map ({ $_->textContent } @arg);
            $x =~ s/^PS(P|D)_//o; my $z = $1;
            $p =~ s/^.*%Y//o;
    
            my ($f, $t) = split (m/%/o, $p);
            $t =~ s/^MP/T/o;
            $t = "T0" if ($t eq 'T');
    
            
            my @ss = &F ('//EN-decl[string (EN-N)="?"]//shape-spec', $v, $d);
#           print &Dumper ([$f, $t, $x, $v, scalar (@ss)]);
    
    
            $stmt = &n ("<pointer-a-stmt><E-1><named-E><N><n>$v</n></N></named-E></E-1> "
                      . "<a>=&gt;</a> <E-2><named-E><N><n>YDSURFVARS</n></N><R-LT><component-R>%<ct>GS${z}_$x</ct></component-R>"
                      . "<component-R>%<ct>V$f</ct></component-R><component-R>%<ct>F$t</ct></component-R>"
                      . "<component-R>%<ct>DATA</ct></component-R> <array-R>(<section-subscript-LT>"
                      . "<section-subscript><lower-bound><named-E><N><n>KBL</n></N></named-E></lower-bound>"
                      . "</section-subscript></section-subscript-LT>)</array-R></R-LT></named-E></E-2></pointer-a-stmt>");
    
            my ($lt) = &F ('.//section-subscript-LT', $stmt);
            
            for (@ss)
              {
                $lt->insertBefore (&t (', '), $lt->firstChild);
                $lt->insertBefore (&n ("<section-subscript>:</section-subscript>"), $lt->firstChild);
              }
    
          }
        elsif ($arg[1]->textContent =~ m/^PGMV/o)
          {
            (my $p, undef, my $v) = map ({ $_->textContent } @arg);
            my ($t, $y) = split (m/%/o, $p);
            $t =~ s/^Y//o;
            $y =~ s/^M//o; 
            if ($y =~ s/(L|M)$//o)
              {
                my $d = $1;
                next unless (($t eq 'T0') || ($t eq 'T9'));
                $t = "D$d";
              }
    
            $stmt = &n ("<pointer-a-stmt><E-1><named-E><N><n>$v</n></N></named-E></E-1> <a>=&gt;</a> "
                      . "<E-2><named-E><N><n>YDVARS</n></N><R-LT><component-R>%<ct>$y</ct></component-R>"
                      . "<component-R>%<ct>$t</ct></component-R></R-LT></named-E></E-2></pointer-a-stmt>");
    
          }
        elsif ($arg[1]->textContent =~ m/^PGFL/o)
          {
            (my $p, undef, my $v) = map ({ $_->textContent } @arg);
            my ($y, $t) = split (m/%/o, $p);
            $t =~ s/^MP//o or die;
            $t ||= '0';
            $y =~ s/^Y//o;
    
    
=pod
    
      REAL(KIND=JPRB), POINTER :: P(:)  => NULL()  ! Basic field at t
      REAL(KIND=JPRB), POINTER :: T0(:) => NULL()  ! Basic field at t (alias of P)
      REAL(KIND=JPRB), POINTER :: T1(:) => NULL()  ! Basic field at t+dt
      REAL(KIND=JPRB), POINTER :: T9(:) => NULL()  ! Basic field at t-dt
      REAL(KIND=JPRB), POINTER :: PH9(:)=> NULL()  ! Basic field for physics
      REAL(KIND=JPRB), POINTER :: DL(:) => NULL()  ! Zonal derivative field
      REAL(KIND=JPRB), POINTER :: DM(:) => NULL()  ! Meridional derivative field
    
    
=cut
    
            die unless (grep { $_ eq $y } @fn);
    
            if ($t =~ m/^\d$/o)
              {
                $t = "T$t";
              }
            elsif (($t eq 'L') or ($t eq 'M'))
              {
                $t = "D$t";
              }
            else
              {
                die;
              }
    
            $stmt = &n ("<pointer-a-stmt><E-1><named-E><N><n>$v</n></N></named-E></E-1> <a>=&gt;</a> "
                      . "<E-2><named-E><N><n>YDVARS</n></N><R-LT><component-R>%<ct>$y</ct></component-R>"
                      . "<component-R>%<ct>$t</ct></component-R></R-LT></named-E></E-2></pointer-a-stmt>");
    
          }
    
        if ($stmt)
          {
            my ($E1) = &F ('./E-1', $stmt);
            my $sp = $E1->nextSibling;
            $sp->setData (' ' x (12 - length ($E1->textContent)));
            $call->replaceNode ($stmt);
          }
      }
        
    
    'FileHandle'->new (">$F90.new")->print ($d->textContent);
    

  }
