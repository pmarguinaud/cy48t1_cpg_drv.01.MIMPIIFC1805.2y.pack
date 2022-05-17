SUBROUTINE INCLI8(YDGEOMETRY)

!**** *INCLI8*

!     PURPOSE.
!     --------

!     Interpolation of the 3 coefficients of the profile of ozone
!     This routine reads a,b,c on a 2.5 climatological grid

!**   INTERFACE.
!     ----------

!     CALL INCLI8

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     GEO923
!     INTER2
!     FA-LFI package (FAITOU,FACILE,FAIENC,LFILAF,FAIRME)
!     CHIEN

!     AUTHORS.
!     --------

!          Y.Bouteloup      29/03/2002

!     MODIFICATIONS.
!     --------------
!        M.Hamrud    01-Oct-2003  CY28 Cleaning
!        Y.Bouteloup 01-Jul-2004  Modification of interpolation 
!      & F. Bouyssel              algorithm : INTER1 -> INTER2
!        D. Giard    15-Sep-2004  Phasing to 28T3  
!        K. Yessad (Jan 2010): externalisation of group EGGX in XRD/IFSAUX
!        G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCT0   , ONLY : NQUAD
USE YOMLUN   , ONLY : NULOUT
USE YOMVERT  , ONLY : VP00
USE YOMCLI   , ONLY : YRCLI   

IMPLICIT NONE

! JPBX/JPBY : number of extra longitudes/latitudes on each side of the data
!          grid required by interpolation (4 or 12 points -> 2/2)

TYPE(GEOMETRY), INTENT(IN)   :: YDGEOMETRY
INTEGER(KIND=JPIM) :: JPBX
INTEGER(KIND=JPIM) :: JPBY

PARAMETER(JPBX=2,JPBY=2)

!  Initial and interpolated datasets
REAL(KIND=JPRB),ALLOCATABLE :: ZA(:,:), ZB(:,:), ZC(:,:)
!  Geometry
REAL(KIND=JPRB) :: ZBO(YDGEOMETRY%YRDIM%NDGLG+1),ZLONG(YDGEOMETRY%YRDIM%NDGLG),&
 & ZMU(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON),ZSLA(YDGEOMETRY%YRDIM%NDGLG*YRCLI%NPINT)&
 & ,ZSLO(YDGEOMETRY%YRDIM%NDLON*YRCLI%NPINT,YDGEOMETRY%YRDIM%NDGLG),ZCLO(YDGEOMETRY%YRDIM%NDLON*YRCLI%NPINT,YDGEOMETRY%YRDIM%NDGLG)&
 & ,ZS(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON,3)  
!   Local variables
REAL(KIND=JPRB) :: ZEPS, ZLAT, ZLON
REAL(KIND=JPRB) :: ZHOOK_HANDLE

CHARACTER :: CLNOMF*16,CLNOMC*16,CLFORM*12

INTEGER(KIND=JPIM) :: IARI, IARP, ICO, IDATX, IDATY, IFLD, IINF,&
 & ILENE, ILENT, IMES, INIQ, INIV, INJQ, INUM, IOS, IREP, JX, JY

LOGICAL :: LLCOSP, LLIMST, LLOP, LLOZ, LLPOLE

!     ------------------------------------------------------------------

#include "chien.h"

#include "abor1.intfb.h"
#include "geo923.intfb.h"
#include "inter2.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INCLI8',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDVAB=>YDGEOMETRY%YRVAB&
& )
ASSOCIATE(NDGENG=>YDDIM%NDGENG, NDGLG=>YDDIM%NDGLG, NDGSAG=>YDDIM%NDGSAG,   NDLON=>YDDIM%NDLON, NSMAX=>YDDIM%NSMAX,   &
& NFLEVG=>YDDIMV%NFLEVG,   NHTYP=>YDGEM%NHTYP, NLOENG=>YDGEM%NLOENG, NMENG=>YDGEM%NMENG,   NSTTYP=>YDGEM%NSTTYP,      &
& RLOCEN=>YDGEM%RLOCEN, RMUCEN=>YDGEM%RMUCEN,   RSTRET=>YDGEM%RSTRET)
!     ------------------------------------------------------------------

!     1. SET INITIAL VALUES.
!        -------------------

!     1.1 Constants

ZEPS=1.E-6_JPRB

!     1.3 ARPEGE and data files

INUM=3
IMES=1
IARP=31
IREP=0
IARI=0
CLNOMF='Const.Clim'
CLNOMC='Const.Clim.Surfa'
LLIMST=.TRUE.

IINF=-1
LLPOLE=.TRUE.
LLCOSP=.FALSE.

IF (YRCLI%LIEEE) THEN
  CLFORM='UNFORMATTED'
  CALL ABOR1(' LIEEE=.TRUE. : Input files not yet ready !')
ELSE
  CLFORM='FORMATTED'
ENDIF

!     1.4 Geometry

ICO=YRCLI%NPINT*YRCLI%NPINT
INJQ=NDGLG*YRCLI%NPINT
INIQ=NDLON*YRCLI%NPINT
ILENT=NDGLG*NDLON
IDATX=YRCLI%NDATX+2*JPBX
IDATY=YRCLI%NDATY+2*JPBY
IFLD=1
ILENE=0
CALL GEO923(YDGEOMETRY,YRCLI%NPINT,ILENE,ZBO,ZLONG,ZMU,ZSLA,ZSLO,ZCLO)

!     ------------------------------------------------------------------

!     2. CHECK REQUIRED FILES.
!        ----------------------

!     2.1 Check datasets

!  input file
LLOZ=.FALSE.
LLOP=.FALSE.
IOS= 0
INQUIRE(FILE='abc_coef',IOSTAT=IOS,EXIST=LLOZ,OPENED=LLOP)
LLOZ= LLOZ .AND. (IOS == 0) .AND. .NOT.LLOP

!     2.2 Control and print

IF (.NOT.LLOZ) THEN
  CALL ABOR1(' INCLI8 : No input file !')
ENDIF

!     2.3 Open ARPEGE file

CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'OLD',.TRUE.,LLIMST,&
 & IMES,IARP,IARI,CLNOMC)  
CALL CHIEN(CLNOMC,NSTTYP,RMUCEN,RLOCEN,RSTRET,NSMAX,&
 & NDGLG,NDLON,NLOENG,NMENG,NHTYP,NFLEVG,VP00,YDVAB%VALH,YDVAB%VBH,&
 & NQUAD,IINF,NDGSAG,NDGENG,ZEPS,LLPOLE,NULOUT)  
IF (LLPOLE) CALL ABOR1(' CLIM. FILES MUST NOT HAVE POLES !')

!     ------------------------------------------------------------------

!     3. READ AND INTERPOLATE NEW FIELDS.
!        --------------------------------

!     3.1 Allocate and initialize arrays for input data

ALLOCATE ( ZA(IDATX,IDATY) )
ALLOCATE ( ZB(IDATX,IDATY) )
ALLOCATE ( ZC(IDATX,IDATY) )
ZA(:,:)=0.0_JPRB
ZB(:,:)=0.0_JPRB
ZC(:,:)=0.0_JPRB

!     3.2 Read input data

OPEN(UNIT=11,FILE='abc_coef',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  CALL ABOR1(' LIEEE=.TRUE. : Not yet !')
ELSE
  DO JY=JPBY+1,JPBY+YRCLI%NDATY
    DO JX=JPBX+1,JPBX+YRCLI%NDATX
       READ(11,*) ZLON,ZLAT,ZA(JX,JY),ZB(JX,JY),ZC(JX,JY)
    ENDDO
  ENDDO  
ENDIF
CLOSE(11)

!     3.3 Interpolate a,b and c

CALL INTER2(NDGLG,NDLON,IDATY,IDATX,IFLD,YRCLI%NPINT,ILENE,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,1),ZA,ZSLA,ZSLO,ZCLO)

CALL INTER2(NDGLG,NDLON,IDATY,IDATX,IFLD,YRCLI%NPINT,ILENE,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,2),ZB,ZSLA,ZSLO,ZCLO)

CALL INTER2(NDGLG,NDLON,IDATY,IDATX,IFLD,YRCLI%NPINT,ILENE,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,3),ZC,ZSLA,ZSLO,ZCLO)

!     3.4 Deallocate local arrays

DEALLOCATE ( ZA )
DEALLOCATE ( ZB )
DEALLOCATE ( ZC )

!     ------------------------------------------------------------------

!     5. WRITE.
!        ------

!     5.1 Write ozone coefficients : a, b, c

INIV=0
CALL FAIENC(IREP,INUM,'SURF',INIV,'A.OF.OZONE  ',ZS(1,1),LLCOSP)
CALL FAIENC(IREP,INUM,'SURF',INIV,'B.OF.OZONE  ',ZS(1,2),LLCOSP)
CALL FAIENC(IREP,INUM,'SURF',INIV,'C.OF.OZONE  ',ZS(1,3),LLCOSP)

!     5.2 Close ARPEGE file

CALL LFILAF(IREP,INUM,.TRUE.)
CALL FAIRME(IREP,INUM,'KEEP')

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INCLI8',1,ZHOOK_HANDLE)
END SUBROUTINE INCLI8
