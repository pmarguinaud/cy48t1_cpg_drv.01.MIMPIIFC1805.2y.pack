SUBROUTINE EINTER2(YDEGEO,PZ0,KLAT0,KLON0,KNCH,PZI,KLAT,KLON,KNPT,&
 & PDELX,PDELY,KPP,KNCO,KULOUT)  

!     PURPOSE.
!     --------

!       Interpolation of a field from a finer grid.

!**   INTERFACE.
!     ----------
!     CALL EINTER2(...)

!        PZ0   = initial fields
!        KLAT0 = number of latitudes  of the initial grid
!                  2  latitudes are added beyond each pole
!        KLON0 = number of longitudes of the initial grid
!                  4 longitudes are duplicated : 2 on each side
!        KNCH  = number of fields to interpolate
!        PZI   = interpolated fields                           (output)
!        KLAT  = number of latitudes  of the final grid
!        KLON  = number of longitudes of the final grid
!        KNPT  = number of points     of the final grid (KLAT*KLON)
!        PDELX,PDELY=size of final grid (m or rad.)
!        KPP   = size of the box
!        KNCO  = KPP*KPP
!        KULOUT= unit of output file

!     METHOD.
!     -------

!       The field is first calculated on a subgrid of the final grid.
!       Each value is interpolated from the 12 nearest points of the initial
!       grid.
!       Then the values of the box (subgrid) are averaged

!     EXTERNALS.
!     ----------

!       EGGPACK

!     AUTHORS.
!     --------
!      J.F. GELEYN  M. DEQUE  1 FEB 91.

!     MODIFICATION
!     ------------
!     M. Janousek 01-04-16: Replace EGGRVS by EGGPACK functions
!     JD. Gril 18-Nov-2004 Modifs for Mercator Rotated Tilted
!      ----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YEMCLI   , ONLY : NILCLI1  ,NJLCLI1
USE YOMCLI   , ONLY : YRCLI
USE YOMCST   , ONLY : RPI      ,RA
USE YEMGEO   , ONLY : TEGEO
USE YOMCT0   , ONLY : LRPLANE    
USE EGGANGLES , ONLY : ANGLE_DOMAIN 
USE EGGPACK  , ONLY : LOLA,PARAM_PROJ,XY,XY_TO_LATLON,REF_DATAS, &
 & XY_NEW_TO_STD_ORIGIN

!      ----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TEGEO)       ,INTENT(IN)    :: YDEGEO
INTEGER(KIND=JPIM),INTENT(IN)    :: KLAT0 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON0 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNCH 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNPT 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0(KLON0,KLAT0,KNCH) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZI(KNPT,KNCH) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLAT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELY 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNCO 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 

!      ----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZWGT12(KNCO,12)
REAL(KIND=JPRB) :: ZXCO  (KNCO*KLAT*KLON),ZYCO  (KNCO*KLAT*KLON)
INTEGER(KIND=JPIM) :: IXX(KNCO),IYY(KNCO)

INTEGER(KIND=JPIM) :: ICO, II1, IJ1, IM, INDICE, IPAUT, IXM, IYM,&
 & IZ, J, JCH, JFL, JI, JJ, JK, JL, JM  

REAL(KIND=JPRB) :: ZCONO, ZCONO2, ZCX, ZCY, ZGX, ZGY, ZX, ZXX, ZY, ZYY, ZEPS

TYPE (LOLA)                :: YL_TLKRES, YL_TLREFC, YL_NWORIG
TYPE (LOLA) , ALLOCATABLE  :: YL_TLGRID_LOLA(:)
TYPE (PARAM_PROJ)          :: YL_TLMODDOM
TYPE (XY)   , ALLOCATABLE  :: YL_TLGRID_XY(:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------------

#include "abor1.intfb.h"

!      ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EINTER2',0,ZHOOK_HANDLE)

ASSOCIATE( &
 & ELON0=>YDEGEO%ELON0, ELAT0=>YDEGEO%ELAT0, ELAT1=>YDEGEO%ELAT1, &
 & ELON1=>YDEGEO%ELON1, ELONC=>YDEGEO%ELONC, ELATC=>YDEGEO%ELATC, &
 & LMRT=>YDEGEO%LMRT )
!      ----------------------------------------------------------------------

!     1. Setting initial values to zero and extending the initial grid

ZEPS=EPSILON(1.0_JPRB)*100.0_JPRB
DO JCH=1,KNCH
  DO J=1,KNPT
    PZI(J,JCH)=0.0_JPRB
  ENDDO
ENDDO

IF (YRCLI%LGLOBE) THEN
  II1=1
  IJ1=1
  ZGX=REAL(KLON0-4,JPRB)
  ZGY=REAL(KLAT0-4,JPRB)
ELSE
  II1=NILCLI1
  IJ1=NJLCLI1
  ZGX=REAL(YRCLI%NGLOBX,JPRB)
  ZGY=REAL(YRCLI%NGLOBY,JPRB)
ENDIF

IF ( YRCLI%LGLOBE .AND. ((KLON0-4-YRCLI%NGLOBX)*(KLAT0-4-YRCLI%NGLOBY) /= 0) )&
 & WRITE(KULOUT,'('' EINTER2 : CAUTION !'',/,&
 & '' EITHER YOU ARE USING OLD E923, OR'',&
 & '' THERE IS A BUG IN THE SIZE OF DATASET (ARGUMENT 2 OR 3)'')')  

IF (YRCLI%LGLOBE) THEN
! global dataset
  IPAUT=(KLON0-4)/2
  DO JFL=1,KNCH
    DO JI=3,IPAUT+2
      PZ0(JI,1,JFL)=PZ0(JI+IPAUT,4,JFL)
      PZ0(JI+IPAUT,1,JFL)=PZ0(JI,4,JFL)
      PZ0(JI,2,JFL)=PZ0(JI+IPAUT,3,JFL)
      PZ0(JI+IPAUT,2,JFL)=PZ0(JI,3,JFL)
      PZ0(JI,KLAT0,JFL)=PZ0(JI+IPAUT,KLAT0-3,JFL)
      PZ0(JI+IPAUT,KLAT0,JFL)=PZ0(JI,KLAT0-3,JFL)
      PZ0(JI,KLAT0-1,JFL)=PZ0(JI+IPAUT,KLAT0-2,JFL)
      PZ0(JI+IPAUT,KLAT0-1,JFL)=PZ0(JI,KLAT0-2,JFL)
    ENDDO

    DO JJ=1,KLAT0
      PZ0(1,JJ,JFL)=PZ0(KLON0-3,JJ,JFL)
      PZ0(2,JJ,JFL)=PZ0(KLON0-2,JJ,JFL)
      PZ0(KLON0-1,JJ,JFL)=PZ0(3,JJ,JFL)
      PZ0(KLON0,JJ,JFL)=PZ0(4,JJ,JFL)
    ENDDO
  ENDDO
ELSE
! non-global input grid: just duplicate border values into margins
  DO JFL=1,KNCH
    DO JI=3,KLON0-3
      PZ0(JI,1,JFL)=PZ0(JI,3,JFL)
      PZ0(JI,2,JFL)=PZ0(JI,3,JFL)
      PZ0(JI,KLAT0-1,JFL)=PZ0(JI,KLAT0-2,JFL)
      PZ0(JI,KLAT0,JFL)=PZ0(JI,KLAT0-2,JFL)
    ENDDO
    DO JJ=1,KLAT0
      PZ0(1,JJ,JFL)=PZ0(3,JJ,JFL)
      PZ0(2,JJ,JFL)=PZ0(3,JJ,JFL)
      PZ0(KLON0-1,JJ,JFL)=PZ0(KLON0-2,JJ,JFL)
      PZ0(KLON0,JJ,JFL)=PZ0(KLON0-2,JJ,JFL)
    ENDDO
  ENDDO
ENDIF

!     2. Bulk summation.

ICO=KNCO*KLAT*KLON
ALLOCATE(YL_TLGRID_XY(ICO))
ALLOCATE(YL_TLGRID_LOLA(ICO))
INDICE=0
ZCONO=1.0_JPRB/REAL(KPP,JPRB)
DO JJ=1,KLAT
  DO JI=1,KLON

!     Calculation of the coordinates of the points in the box
    DO JL=1,KPP
      DO JK=1,KPP
        IZ=(JL-1)*KPP + JK + ((JI-1+(JJ-1)*KLON))*KPP*KPP
        YL_TLGRID_XY(IZ)%X=(REAL(JI,JPRB)-1.5_JPRB)*PDELX &
         & +(REAL(JK,JPRB)-0.5_JPRB)*PDELX*ZCONO  
        YL_TLGRID_XY(IZ)%Y=(REAL(JJ,JPRB)-1.5_JPRB)*PDELY &
         & +(REAL(JL,JPRB)-0.5_JPRB)*PDELY*ZCONO  
      ENDDO
    ENDDO

  ENDDO
ENDDO

!     Coordinate transform.
IF (LRPLANE) THEN
  IF (LMRT .AND. (ABS(ELAT0) >= ZEPS)) THEN
    WRITE(KULOUT,*) 'EINTER2: ELAT0=',ELAT0,&
    & ' MUST BE EQUAL ZERO IF LMRT IS TRUE!'
    CALL ABOR1('EINTER2: LMRT & ELAT0 INCONSISTENT')
  ENDIF 
  YL_TLKRES%LON=ELON0*180._JPRB/RPI
  YL_TLKRES%LAT=ELAT0*180._JPRB/RPI
  YL_TLKRES=ANGLE_DOMAIN(YL_TLKRES,RPI,'+-','D')
  YL_NWORIG%LON=ELON1*180._JPRB/RPI
  YL_NWORIG%LAT=ELAT1*180._JPRB/RPI
  YL_NWORIG=ANGLE_DOMAIN(YL_NWORIG,RPI,'+-','D')
  YL_TLREFC%LON=ELONC*180._JPRB/RPI
  YL_TLREFC%LAT=ELATC*180._JPRB/RPI
  YL_TLREFC=ANGLE_DOMAIN(YL_TLREFC,RPI,'+-','D')
  YL_TLMODDOM=REF_DATAS(YL_TLKRES,RA,YL_TLREFC,LMRT)
  YL_TLGRID_LOLA=&
   & XY_TO_LATLON(XY_NEW_TO_STD_ORIGIN(YL_NWORIG,YL_TLGRID_XY,YL_TLMODDOM,RPI),&
   & YL_TLMODDOM,RPI)  
ELSE
  YL_TLGRID_LOLA(:)%LON=ELON1+YL_TLGRID_XY(:)%X
  YL_TLGRID_LOLA(:)%LAT=ELAT1+YL_TLGRID_XY(:)%Y
ENDIF
ZXCO(:)=YL_TLGRID_LOLA(:)%LON
ZYCO(:)=YL_TLGRID_LOLA(:)%LAT
DEALLOCATE(YL_TLGRID_XY)
DEALLOCATE(YL_TLGRID_LOLA)

DO JJ=1,KLAT
  DO JI=1,KLON

!     Loop in the box
    DO JM=1,KNCO
! the interpolation values are localised.
      IM=JM+(JI-1+(JJ-1)*KLON)*KNCO
      ZX=MOD( ZGX/RPI/2.0_JPRB*(ZXCO(IM)+2.0_JPRB*RPI) , ZGX )
      IXX(JM)=INT( 3.5_JPRB + ZX - II1 )
      ZXX=3.5_JPRB + ZX - II1 - IXX(JM)
      ZCX=1.0_JPRB-ZXX
      ZY=ZGY/RPI*(RPI/2.0_JPRB+ZYCO(IM))
      IYY(JM)=INT( 3.5_JPRB + ZY - IJ1 )
      ZYY=3.5_JPRB + ZY - IJ1 - IYY(JM)
      ZCY=1.0_JPRB-ZYY
! the linear interpolation weights are set (multiplied by 24 to save time).
      ZWGT12(JM,1)=3._JPRB*(8._JPRB*ZXX*ZYY+ZXX*ZCX*(ZXX-ZCX+2.0_JPRB*ZYY)&
       & +ZYY*ZCY*(ZYY-ZCY+2.0_JPRB*ZXX))  
      ZWGT12(JM,2)=3._JPRB*(8._JPRB*ZCX*ZYY+ZCX*ZXX*(ZCX-ZXX+2.0_JPRB*ZYY)&
       & +ZYY*ZCY*(ZYY-ZCY+2.0_JPRB*ZCX))  
      ZWGT12(JM,3)=3._JPRB*(8._JPRB*ZCX*ZCY+ZCX*ZXX*(ZCX-ZXX+2.0_JPRB*ZCY)&
       & +ZCY*ZYY*(ZCY-ZYY+2.0_JPRB*ZCX))  
      ZWGT12(JM,4)=3._JPRB*(8._JPRB*ZXX*ZCY+ZXX*ZCX*(ZXX-ZCX+2.0_JPRB*ZCY)&
       & +ZCY*ZYY*(ZCY-ZYY+2.0_JPRB*ZXX))  
      ZWGT12(JM, 5)=ZXX*ZCX*(ZCX-ZXX-6._JPRB*ZYY)
      ZWGT12(JM, 6)=ZYY*ZCY*(ZCY-ZYY-6._JPRB*ZXX)
      ZWGT12(JM, 7)=ZYY*ZCY*(ZCY-ZYY-6._JPRB*ZCX)
      ZWGT12(JM, 8)=ZCX*ZXX*(ZXX-ZCX-6._JPRB*ZYY)
      ZWGT12(JM, 9)=ZCX*ZXX*(ZXX-ZCX-6._JPRB*ZCY)
      ZWGT12(JM,10)=ZCY*ZYY*(ZYY-ZCY-6._JPRB*ZCX)
      ZWGT12(JM,11)=ZCY*ZYY*(ZYY-ZCY-6._JPRB*ZXX)
      ZWGT12(JM,12)=ZXX*ZCX*(ZCX-ZXX-6._JPRB*ZCY)
    ENDDO

    INDICE=INDICE+1
    DO JM=1,KNCO
      IXM=IXX(JM)
      IYM=IYY(JM)
      DO JFL=1,KNCH
        PZI(INDICE,JFL)=PZI(INDICE,JFL)+&
         & (ZWGT12(JM, 1)*PZ0(IXM+1,IYM+1,JFL)&
         & +ZWGT12(JM, 2)*PZ0(IXM  ,IYM+1,JFL)&
         & +ZWGT12(JM, 3)*PZ0(IXM  ,IYM  ,JFL)&
         & +ZWGT12(JM, 4)*PZ0(IXM+1,IYM  ,JFL)&
         & +ZWGT12(JM, 5)*PZ0(IXM+2,IYM+1,JFL)&
         & +ZWGT12(JM, 6)*PZ0(IXM+1,IYM+2,JFL)&
         & +ZWGT12(JM, 7)*PZ0(IXM  ,IYM+2,JFL)&
         & +ZWGT12(JM, 8)*PZ0(IXM-1,IYM+1,JFL)&
         & +ZWGT12(JM, 9)*PZ0(IXM-1,IYM  ,JFL)&
         & +ZWGT12(JM,10)*PZ0(IXM  ,IYM-1,JFL)&
         & +ZWGT12(JM,11)*PZ0(IXM+1,IYM-1,JFL)&
         & +ZWGT12(JM,12)*PZ0(IXM+2,IYM  ,JFL))  
      ENDDO
    ENDDO

  ENDDO
ENDDO

!     3. Normalisation.

ZCONO2=1.0_JPRB/REAL(24*KNCO,JPRB)
DO JFL=1,KNCH
  DO J=1,KNPT
    PZI(J,JFL)=PZI(J,JFL)*ZCONO2
  ENDDO
ENDDO

!      ----------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EINTER2',1,ZHOOK_HANDLE)
END SUBROUTINE EINTER2
