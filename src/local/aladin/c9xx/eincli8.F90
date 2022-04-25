SUBROUTINE EINCLI8(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF)

!**** *EINCLI8*

!     PURPOSE.
!     --------

!     Interpolation of the 3 coef of the profil of ozone
!     This routine read the tree coef a,b,c on a 2.5 climatological grid

!**   INTERFACE.
!     ----------

!     CALL EINCLI8

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     FA-LFI package (FAITOU,FACILE,FAIENC,LFILAF,FAIRME)
!     CCHIEN, EGEO3
!     EINTER1 for interpolations

!     AUTHORS.
!     --------

!          Y.Bouteloup      12/07/2004 (from incli8)

!     MODIFICATIONS.
!     --------------
!        D. Giard 04-09-15 cleaning
!S.Ivatek-Sahdan : 04-11-30 pre29 dbg IFLD not initialised
!        D. Giard 05-04-07 new call to EBICLI
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE YOEPHY       , ONLY : TEPHY
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMLUN   , ONLY : NULOUT
USE YOMCLI   , ONLY : YRCLI
USE YOM_YGFL , ONLY : TYPE_GFLD

!     ------------------------------------------------------------------

IMPLICIT NONE

!     JPBY : Extra-latitudes
!          2 when INTER2 or INTER3 are used (even NDATY)
!          1 when INTER4 or INTER5 are used (odd  NDATY)
TYPE(GEOMETRY), INTENT(INOUT)   :: YDGEOMETRY
TYPE(TYPE_GFLD)    ,INTENT(INOUT):: YDGFL
TYPE(TEPHY)    ,INTENT(INOUT)   :: YDEPHY
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
INTEGER(KIND=JPIM) :: JPBX
INTEGER(KIND=JPIM) :: JPBY
PARAMETER (JPBX=2,JPBY=2)

!  Initial and interpolated datasets
REAL(KIND=JPRB),ALLOCATABLE :: ZA(:,:), ZB(:,:), ZC(:,:), ZRES(:)
!  Final data
REAL(KIND=JPRB) :: ZS(YDGEOMETRY%YRDIM%NDLON*YDGEOMETRY%YRDIM%NDGLG,3)
REAL(KIND=JPRB) :: ZLAT, ZLON
REAL(KIND=JPRB) :: ZHOOK_HANDLE

CHARACTER :: CLNOMF*16,CLNOMC*16,CLFORM*12
CHARACTER :: CLPREF(3)*8,CLSUFF(3)*12

INTEGER(KIND=JPIM) :: INIVL(3)
INTEGER(KIND=JPIM) :: IARI, IARP, ICO, IDATX, IDATY, IFLD, IINF, IMES,&
 & INUM, IOS, IREP, ITFING, IXFING, IYFING, J, JX, JY

LOGICAL :: LLBIP(3),LLWRI(3),LLPAC(3)
LOGICAL :: LLIMST, LLOP, LLOZ

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "cchien.intfb.h"
#include "ebicli.intfb.h"
#include "egeo923.intfb.h"
#include "einter2.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EINCLI8',0,ZHOOK_HANDLE)

ASSOCIATE( &
 & NDGUNG=>YDGEOMETRY%YRDIM%NDGUNG, NDGUXG=>YDGEOMETRY%YRDIM%NDGUXG, NDLUNG=>YDGEOMETRY%YRDIM%NDLUNG, &
 & NDLUXG=>YDGEOMETRY%YRDIM%NDLUXG, NDGLG=>YDGEOMETRY%YRDIM%NDGLG, NDLON=>YDGEOMETRY%YRDIM%NDLON, &
 & EDELX=>YDGEOMETRY%YREGEO%EDELX, EDELY=>YDGEOMETRY%YREGEO%EDELY)
!     ------------------------------------------------------------------

!     1. SET INITIAL VALUES.
!        -------------------

!     1.1 ARPEGE and data files

INUM=3
IMES=1
IARP=31
IREP=0
IARI=0
CLNOMF='Const.Clim'
CLNOMC='Const.Clim.Surfa'
LLIMST=.TRUE.

IINF=-1

IF (YRCLI%LIEEE) THEN
  CLFORM='UNFORMATTED'
ELSE
  CLFORM='FORMATTED'
ENDIF

!     1.2 Dimensions for data interpolation

! Original grid
IDATX=YRCLI%NDATX+2*JPBX
IDATY=YRCLI%NDATY+2*JPBY

! Final grid:
IXFING=NDLUXG-NDLUNG+1
IYFING=NDGUXG-NDGUNG+1
ITFING=IXFING*IYFING

! NPINT: size of the interpolation boxes, overwritten
YRCLI%NPINT=1
ICO=YRCLI%NPINT*YRCLI%NPINT

CALL EGEO923(YDGEOMETRY%YRGEM)

!     ------------------------------------------------------------------

!     2. CHECK REQUIRED FILES.
!        ---------------------

!     2.1 Check datasets

LLOZ=.FALSE.
LLOP=.FALSE.
IOS= 0
INQUIRE(FILE='abc_coef',IOSTAT=IOS,EXIST=LLOZ,OPENED=LLOP)
LLOZ= LLOZ .AND. (IOS == 0) .AND. .NOT.LLOP

!     2.2 Control and print

IF (.NOT.LLOZ) THEN
  CALL ABOR1(' EINCLI8 : No input file !')
ENDIF

!     2.3 Open ARPEGE fields

CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'OLD',.TRUE.,LLIMST,&
 & IMES,IARP,IARI,CLNOMC)
CALL CCHIEN(YDGEOMETRY,CLNOMC,INUM,IINF)           

!     ------------------------------------------------------------------

!     3. READ AND INTERPOLATE NEW FIELDS.
!        --------------------------------

!     3.1 Allocate and initialize arrays

ALLOCATE ( ZA(IDATX,IDATY) )
ALLOCATE ( ZB(IDATX,IDATY) )
ALLOCATE ( ZC(IDATX,IDATY) )
ALLOCATE ( ZRES(ITFING) )

ZA(:,:)=0._JPRB
ZB(:,:)=0._JPRB
ZC(:,:)=0._JPRB
ZS(:,:)=0.0_JPRB

!     3.2 Read input data

OPEN(UNIT=11,FILE='abc_coef',FORM=CLFORM)

IF (YRCLI%LIEEE) THEN
  CALL ABOR1(' LIEEE=.TRUE. : Not yet !')
ELSE
  DO JY=JPBY+YRCLI%NDATY,JPBY+1,-1
     DO JX=JPBX+1,JPBX+YRCLI%NDATX
        READ(11,*) ZLON,ZLAT,ZA(JX,JY),ZB(JX,JY),ZC(JX,JY)
     ENDDO
  ENDDO      
ENDIF

CLOSE(11)


!     3.3 Interpolate (a, b, c) coefficients

! a
IFLD=1
ZRES(:)=0.0_JPRB
CALL EINTER2(YDGEOMETRY%YREGEO,ZA,IDATY,IDATX,IFLD,ZRES,&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)
DO J=1,ITFING
  ZS(J,1)=ZRES(J)
ENDDO

! b
IFLD=1
ZRES(:)=0.0_JPRB
CALL EINTER2(YDGEOMETRY%YREGEO,ZB,IDATY,IDATX,IFLD,ZRES,&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)
DO J=1,ITFING
  ZS(J,2)=ZRES(J)
ENDDO

! c
IFLD=1
ZRES(:)=0.0_JPRB
CALL EINTER2(YDGEOMETRY%YREGEO,ZC,IDATY,IDATX,IFLD,ZRES,&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)
DO J=1,ITFING
  ZS(J,3)=ZRES(J)
ENDDO

!     3.4 Deallocate temporary space

DEALLOCATE ( ZA )
DEALLOCATE ( ZB )
DEALLOCATE ( ZC )
DEALLOCATE ( ZRES )

!     ------------------------------------------------------------------

!     4. BIPERIODIZATION AND WRITING.
!        ----------------------------

!     4.1 Fields description

CLPREF(:)='SURF'
CLSUFF(1)='A.OF.OZONE  '
CLSUFF(2)='B.OF.OZONE  '
CLSUFF(3)='C.OF.OZONE  '
LLBIP(:)=.TRUE.
LLWRI(:)=.TRUE.
LLPAC(:)=.TRUE.
INIVL(:)= 0

!     4.2 Write

IFLD=3
CALL EBICLI(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF,IFLD,INIVL,CLPREF,CLSUFF,INUM,ZS,LLBIP,LLWRI,LLPAC)

!     4.3 Close ALADIN file

CALL LFILAF(IREP,INUM,.TRUE.)
CALL FAIRME(IREP,INUM,'KEEP')

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EINCLI8',1,ZHOOK_HANDLE)
END SUBROUTINE EINCLI8
