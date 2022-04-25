SUBROUTINE INCLI9(YDGEOMETRY)

!**** *INCLI9*

!     PURPOSE.
!     --------

!      Interpolation of aerosol fields

!**   INTERFACE.
!     ----------

!     CALL INCLI9

!     METHOD.
!     -------

!      The data are read on unit 11, if available, then interpolated.
!      Finally they are written on ARPEGE clim file (unit 3).

!     EXTERNALS.
!     ----------

!     GEO923
!     INTER2
!     FA-LFI package (FAITOU,FACILE,FAIENC,LFILAF,FAIRME)
!     CHIEN

!     AUTHORS.
!     --------
!      D. GIARD  SEP 04 (from INCLIR modified by F. BOUYSSEL)

!     MODIFICATIONS.
!     --------------
!      Modified 03-10-01 by M.Hamrud : CY28 Cleaning
!      Modified 04-06-15 by F. Bouyssel : change of format and scaling for input
!        climatological data, AEROS.URBAN is renamed by AEROS.SOOT accordingly
!      Modified 04-09-15 by D. Giard : new version renamed INCLI9 and phased
!                                      INCLIR kept for older versions
!      K. Yessad (Jan 2010): externalisation of group EGGX in XRD/IFSAUX
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMCT0   , ONLY : NQUAD
USE YOMVERT  , ONLY : VP00
USE YOMCLI   , ONLY : YRCLI

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN)   :: YDGEOMETRY
INTEGER(KIND=JPIM) :: JPBX
INTEGER(KIND=JPIM) :: JPBY
INTEGER(KIND=JPIM) :: JPNCH

PARAMETER (JPBX=2,JPBY=2,JPNCH=4)

!     JPBY : Extra-latitudes 
!          2 when INTER2 or INTER3 are used (even NDATY)
!          1 when INTER4 or INTER5 are used (odd  NDATY)
!     JPNCH: Number of fields to be processed (caution if not 4)

REAL(KIND=JPRB) :: ZS(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON,JPNCH)
REAL(KIND=JPRB) :: ZBO(YDGEOMETRY%YRDIM%NDGLG+1),ZLONG(YDGEOMETRY%YRDIM%NDGLG),&
 & ZMU(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON),ZSLA(YDGEOMETRY%YRDIM%NDGLG*YRCLI%NPINT)&
 & ,ZSLO(YDGEOMETRY%YRDIM%NDLON*YRCLI%NPINT,YDGEOMETRY%YRDIM%NDGLG),ZCLO(YDGEOMETRY%YRDIM%NDLON*YRCLI%NPINT,YDGEOMETRY%YRDIM%NDGLG)  
REAL(KIND=JPRB),ALLOCATABLE :: ZFLD(:,:), ZZLA(:), ZZLO(:)
REAL(KIND=JPRB) :: ZEPS, ZSCALE
REAL(KIND=JPRB) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: IARI, IARP, ICO, IDATX, IDATY, IFLD, IINF, ILENE, ILENT,&
 & IMES, INIQ, INIV, INJQ, INUM, IOS, IREP, JF, J, JJ  

CHARACTER :: CLNOMC*16, CLNOMF*10, CLFORM*12
CHARACTER :: CLPRE(JPNCH)*8, CLSUF(JPNCH)*12

LOGICAL :: LLCOSP, LLIMST, LLPOLE, LLDATA, LLOPEN

!     ------------------------------------------------------------------

#include "chien.h"

#include "abor1.intfb.h"
#include "geo923.intfb.h"
#include "inter2.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INCLI9',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  &
 & YDCSGLEG=>YDGEOMETRY%YRCSGLEG, YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE,  &
 & YDSTA=>YDGEOMETRY%YRSTA, YDLAP=>YDGEOMETRY%YRLAP, YDVSPLIP=>YDGEOMETRY%YRVSPLIP,  &
 & YDVSLETA=>YDGEOMETRY%YRVSLETA, YDHSLMER=>YDGEOMETRY%YRHSLMER, YDCSGEOM=>YDGEOMETRY%YRCSGEOM, &
  & YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB, YDGSGEOM=>YDGEOMETRY%YRGSGEOM, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB,  &
  & YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NDGENG=>YDDIM%NDGENG, NDGLG=>YDDIM%NDGLG, NDGSAG=>YDDIM%NDGSAG, &
 & NDLON=>YDDIM%NDLON, NSMAX=>YDDIM%NSMAX, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NHTYP=>YDGEM%NHTYP, NLOENG=>YDGEM%NLOENG, NMENG=>YDGEM%NMENG, &
 & NSTTYP=>YDGEM%NSTTYP, RLOCEN=>YDGEM%RLOCEN, RMUCEN=>YDGEM%RMUCEN, &
 & RSTRET=>YDGEM%RSTRET)
!     ------------------------------------------------------------------
!*
!     1. SET INITIAL VALUES.
!        -------------------

!     Constants

ZEPS=1.E-10_JPRB
ZSCALE=1._JPRB

!     Grids

ICO=YRCLI%NPINT*YRCLI%NPINT
INJQ=NDGLG*YRCLI%NPINT
INIQ=NDLON*YRCLI%NPINT
ILENT=NDGLG*NDLON
IDATX=YRCLI%NDATX+2*JPBX
IDATY=YRCLI%NDATY+2*JPBY

ILENE=0
CALL GEO923(YDGEOMETRY,YRCLI%NPINT,ILENE,ZBO,ZLONG,ZMU,ZSLA,ZSLO,ZCLO)

!     ARPEGE fields

INIV=0
LLCOSP=.FALSE.
DO JF=1,JPNCH
  CLPRE(JF)='SURF    '
ENDDO
CLSUF(1)='AEROS.SEA   '
CLSUF(2)='AEROS.LAND  '
CLSUF(3)='AEROS.SOOT  '
CLSUF(4)='AEROS.DESERT'
ZS(:,:)=0.0_JPRB

!     ARPEGE file identifiers

INUM=3
IREP=0
IARI=0
IARP=1
CLNOMF='Const.Clim'
CLNOMC='Const.Clim.Surfa'
LLIMST=.TRUE.
IMES=1

IINF=0
LLPOLE=.TRUE.

!     Data type

IF (YRCLI%LIEEE) THEN
  CLFORM='UNFORMATTED'
  CALL ABOR1('INCLI9 : STANDARD INPUT FILES NOT YET READY !')
ELSE
  CLFORM='FORMATTED'
ENDIF

!     ------------------------------------------------------------------
!*
!     2. OPEN FILES.
!        -----------

!     ARPEGE file

CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'OLD',.TRUE.,LLIMST,&
 & IMES,IARP,IARI,CLNOMC)  
CALL CHIEN(CLNOMC,NSTTYP,RMUCEN,RLOCEN,RSTRET,NSMAX,&
 & NDGLG,NDLON,NLOENG,NMENG,NHTYP,NFLEVG,VP00,YDVAB%VALH,YDVAB%VBH,&
 & NQUAD,IINF,NDGSAG,NDGENG,ZEPS,LLPOLE,NULOUT)  
IF (LLPOLE) CALL ABOR1(' CLIM. FILES MUST NOT HAVE POLES !')

!     External input data

LLDATA=.FALSE.
LLOPEN=.FALSE.
IOS= 0
INQUIRE(FILE='aero_GL',IOSTAT=IOS,EXIST=LLDATA,OPENED=LLOPEN)
LLDATA= LLDATA .AND. (IOS == 0) .AND. .NOT.LLOPEN

IF (LLDATA) THEN
  OPEN(UNIT=11,FILE='aero_GL',FORM=CLFORM)
ELSE
  CALL ABOR1('INCLI9 : INPUT DATA NOT AVAILABLE !')
ENDIF

!     ------------------------------------------------------------------
!*
!     3. READ AND INTERPOLATE DATA.
!        --------------------------

!     Allocate temporary space

ALLOCATE ( ZFLD(IDATX,IDATY) )
ALLOCATE ( ZZLA(YRCLI%NDATY) )
ALLOCATE ( ZZLO(YRCLI%NDATX) )

!     Read file header

READ(11,9901) ZZLA
READ(11,9901) ZZLO
9901 FORMAT(5F15.4)

!     Loop on fields

DO JF=1,JPNCH

!     Read and scale data
  ZFLD(:,:)=0.0_JPRB
  DO JJ=JPBY+1,JPBY+YRCLI%NDATY
    READ (11,9902) (ZFLD(J,JJ),J=JPBX+1,JPBX+YRCLI%NDATX)
  ENDDO

!     Horizontal interpolation (12 points)
  IFLD=1
  CALL INTER2(NDGLG,NDLON,IDATY,IDATX,IFLD,YRCLI%NPINT,ILENE,ICO,INJQ,INIQ,&
   & NLOENG(1:),ILENT,ZS(1,JF),ZFLD,ZSLA,ZSLO,ZCLO)  

!    Control and scale interpolated data
  DO J=1,ILENE
    ZS(J,JF)=MAX(0.0_JPRB,ZS(J,JF))/ZSCALE
  ENDDO

ENDDO
9902 FORMAT(12F10.6)

!     Deallocate temporary space

DEALLOCATE ( ZFLD )
DEALLOCATE ( ZZLA )
DEALLOCATE ( ZZLO )

!     Close input file

CLOSE(11)

!     ------------------------------------------------------------------
!*
!     4. WRITE FIELDS.
!        -------------

!     Write fields

DO JF=1,JPNCH
  CALL FAIENC (IREP,INUM,CLPRE(JF),INIV,CLSUF(JF),ZS(1,JF),LLCOSP)
ENDDO

!     Close ARPEGE file

CALL LFILAF (IREP,INUM,.TRUE.)
CALL FAIRME (IREP,INUM,'KEEP')

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INCLI9',1,ZHOOK_HANDLE)
END SUBROUTINE INCLI9
