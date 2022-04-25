SUBROUTINE EINTER1(YDEGEO,PZ0,KLAT0,KLON0,KNPT0,KNCH,PZI,KLAT,KLON,KNPT,&
 & PDELX,PDELY,KPP,KNCO,KULOUT)  

!     PURPOSE.
!     --------

!       Interpolation of a field from a finer grid.

!**   INTERFACE.
!     ----------

!     CALL EINTER1 (...)

!        PZ0   = initial fields
!        KLAT0 = number of latitudes  of the initial grid
!        KLON0 = number of longitudes of the initial grid
!        KNPT0 = number of points     of the initial grid (KLAT0*KLON0)
!        KNCH  = number of fields to interpolate
!        PZI   = interpolated fields (output)
!        KLAT  = number of latitudes  of the  final  grid
!        KLON  = number of longitudes of the  final  grid
!        KNPT  = number of points     of the  final  grid (KLAT*KLON)
!        PDELX,PDELY=size of final grid (m or rad.)
!        KPP   = size of the box
!        KNCO  = KPP*KPP
!        KULOUT= unit of output file

!     METHOD.
!     -------

!       The field is first calculated on a subgrid of the final grid.
!       Each value is taken at the nearest point of the initial grid.
!       Then the values of the box (subgrid) are averaged

!     EXTERNAL
!     ---------

!     EGGPACK

!     AUTHORS.
!     --------
!      J.F. GELEYN  M. DEQUE  1 FEB 91.

!     MODIFICATION
!     ------------
!     M. Janousek 01-04-16: Replace EGGRVS by EGGPACK functions
!     JD. Gril 18-Nov-2004 Modifs for Mercator Rotated Tilted
!     F. Taillefer May 2009 : correction in the longitude calculation for local data grid
!     ------------------------------------------------------------------

USE PARKIND1   ,ONLY : JPIM     ,JPRB
USE YOMHOOK    ,ONLY : LHOOK,   DR_HOOK

USE YEMCLI    , ONLY : NILCLI1  ,NILCLI2  ,NJLCLI1
USE YOMCLI    , ONLY : YRCLI
USE YOMCST    , ONLY : RPI      ,RA
USE YEMGEO    , ONLY : TEGEO
USE YOMCT0    , ONLY : LRPLANE
USE EGGANGLES , ONLY : ANGLE_DOMAIN
USE EGGPACK  , ONLY : LOLA,PARAM_PROJ,XY,XY_TO_LATLON,REF_DATAS, &
 & XY_NEW_TO_STD_ORIGIN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TEGEO),       INTENT(IN)    :: YDEGEO
INTEGER(KIND=JPIM),INTENT(IN)    :: KNPT0 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNCH 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNPT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0(KNPT0,KNCH) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLAT0 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON0 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZI(KNPT,KNCH) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLAT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELY 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNCO 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) , ALLOCATABLE :: ZXCO(:), ZYCO(:)
INTEGER(KIND=JPIM) , ALLOCATABLE ::  INBLAT(:), IICO(:)
INTEGER(KIND=JPIM) :: IXY(KNCO)

INTEGER(KIND=JPIM) :: II1, II2, IJ1, ILONG, IM, INDEX, INDICE,&
 & ISUMLAT, IX, IY, IZ, J, JCH, JI, JJ, JK, JL, JM, JS  

REAL(KIND=JPRB) :: ZCONO, ZCONO1, ZGX, ZGY, ZX, ZY, ZEPS

TYPE (LOLA)                :: YL_TLKRES, YL_TLREFC, YL_NWORIG
TYPE (LOLA) , ALLOCATABLE  :: YL_TLGRID_LOLA(:)
TYPE (PARAM_PROJ)          :: YL_TLMODDOM
TYPE (XY)   , ALLOCATABLE  :: YL_TLGRID_XY(:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EINTER1',0,ZHOOK_HANDLE)

ASSOCIATE( &
 & ELON0=>YDEGEO%ELON0, ELAT0=>YDEGEO%ELAT0, ELAT1=>YDEGEO%ELAT1, &
 & ELON1=>YDEGEO%ELON1, ELONC=>YDEGEO%ELONC, ELATC=>YDEGEO%ELATC, &
 & LMRT=>YDEGEO%LMRT )
!     ------------------------------------------------------------------

ZEPS=EPSILON(1.0_JPRB)*100.0_JPRB
ILONG=INT(KLAT/YRCLI%NSLICE)
ALLOCATE (INBLAT(0:YRCLI%NSLICE))
ALLOCATE (IICO(1:YRCLI%NSLICE))
INBLAT(0)=0

!     1.0 Check kind of grid, set lower left corner global grid coordinate

IF (YRCLI%LGLOBE) THEN
  II1=1
  IJ1=1
  ZGX=REAL(KLON0,JPRB)
  ZGY=REAL(KLAT0,JPRB)
ELSE
  II1=NILCLI1
  II2=NILCLI2
  IJ1=NJLCLI1
  ZGX=REAL(YRCLI%NGLOBX,JPRB)
  ZGY=REAL(YRCLI%NGLOBY,JPRB)
ENDIF

IF ( YRCLI%LGLOBE .AND. ((KLON0-YRCLI%NGLOBX)*(KLAT0-YRCLI%NGLOBY) /= 0) )&
 & WRITE(KULOUT,'('' EINTER1 : CAUTION !'',/,&
 & '' EITHER YOU ARE USING OLD E923, OR'',&
 & '' THERE IS A BUG IN THE SIZE OF DATASET (ARGUMENT 2 OR 3)'')')  

!     1.1 Set initial values to zero

PZI(:,:)=0.0_JPRB

!     1.2 Bulk summation.

!     Initialize the number of latitudes and the number of
!     points for each packet.

ISUMLAT=0
DO JS=1,YRCLI%NSLICE-1
  INBLAT(JS)=ILONG
  ISUMLAT=ISUMLAT+INBLAT(JS)
  IICO(JS)=KNCO*KLON*INBLAT(JS)
ENDDO
INBLAT(YRCLI%NSLICE)=KLAT-ISUMLAT
IICO(YRCLI%NSLICE)=KNCO*KLON*INBLAT(YRCLI%NSLICE)

!     Calculation of the coordinates of the points in the box

INDICE=0
INDEX=1
ZCONO=1.0_JPRB/REAL(KPP,JPRB)

DO JS=1,YRCLI%NSLICE

  ALLOCATE(ZXCO(IICO(JS)))
  ALLOCATE(ZYCO(IICO(JS)))
  ALLOCATE(YL_TLGRID_XY(IICO(JS)))
  ALLOCATE(YL_TLGRID_LOLA(IICO(JS)))

  INDEX=INDEX+INBLAT(JS-1)

  DO JJ=INDEX,INDEX-1+INBLAT(JS)
    DO JI=1,KLON

!      Calculation of the coordinates of the points in the box.

      DO JL=1,KPP
        DO JK=1,KPP
          IZ=(JL-1)*KPP + JK + ((JI-1+(JJ-INDEX)*KLON))*KPP*KPP
          YL_TLGRID_XY(IZ)%X=(REAL(JI,JPRB)-1.5_JPRB)*PDELX &
           & +(REAL(JK,JPRB)-0.5_JPRB)*PDELX*ZCONO  
          YL_TLGRID_XY(IZ)%Y=(REAL(JJ,JPRB)-1.5_JPRB)*PDELY &
           & +(REAL(JL,JPRB)-0.5_JPRB)*PDELY*ZCONO  
        ENDDO
      ENDDO

    ENDDO
  ENDDO
!      Coordinate transform

  IF (LRPLANE) THEN
    IF (LMRT .AND. (ABS(ELAT0) >= ZEPS)) THEN
      WRITE(KULOUT,*) 'EINTER1: ELAT0=',ELAT0,&
       & ' MUST BE EQUAL ZERO IF LMRT IS TRUE!'
      CALL ABOR1('EINTER1: LMRT & ELAT0 INCONSISTENT')
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

!     Loop in the box : the nearest value is taken.
  IF (YRCLI%LGLOBE) THEN
    DO JJ=INDEX,INDEX-1+INBLAT(JS)
      DO JI=1,KLON

        DO JM=1,KNCO
! ZXCO(1:KNCO), ZYCO(1:KNCO)= geog coord (rad) of the pts in the box
! find the nearest pt in initial grid: (IX,IY)
! first compute zx,zy distances from lower left corner of global grid
! the point (1,1) would have coord ZX=ii1-0.5, zy=ij1-0.5
! NGLOBX grid points cover 2*RPI radians longitude
! ZX=II1-0.5 is at grid position ixx=1
! NGLOBY points cover PI radians latitude
! Origin of global grid is at Latitude -PI/2
          IM=JM+(JI-1+(JJ-INDEX)*KLON)*KNCO
          ZXCO(IM)=MOD(ZXCO(IM),2.0_JPRB*RPI)
          ZX=MOD(ZGX/RPI/2.0_JPRB*(ZXCO(IM)+2.0_JPRB*RPI),ZGX)
          IX=NINT(1.5_JPRB + ZX - II1)
          ZY=ZGY/RPI*(RPI/2.0_JPRB+ZYCO(IM))
          IY=NINT(1.5_JPRB + ZY - IJ1)
          IXY(JM)=(IY-1)*KLON0+IX
        ENDDO

        INDICE=INDICE+1
        DO JCH=1,KNCH
          DO JM=1,KNCO
            PZI(INDICE,JCH)=PZI(INDICE,JCH)+PZ0(IXY(JM),JCH)
          ENDDO
        ENDDO

      ENDDO
    ENDDO
  ELSE
    DO JJ=INDEX,INDEX-1+INBLAT(JS)
      DO JI=1,KLON

        DO JM=1,KNCO
          IM=JM+(JI-1+(JJ-INDEX)*KLON)*KNCO
          ZXCO(IM)=MOD(ZXCO(IM),2.0_JPRB*RPI)
          ZX=MOD(ZGX/RPI/2.0_JPRB*(ZXCO(IM)+2.0_JPRB*RPI),ZGX)
          IF(II2 > II1) THEN
            IX=NINT(1.5_JPRB + ZX - II1)
          ELSE
            IX=NINT(1.5_JPRB + ZX - II1)+(1+SIGN(1,II2-NINT(0.5_JPRB + ZX)))&
             & /2*YRCLI%NGLOBX  
          ENDIF
          ZY=ZGY/RPI*(RPI/2.0_JPRB+ZYCO(IM))
          IY=NINT(1.5_JPRB + ZY - IJ1)
          IXY(JM)=(IY-1)*KLON0+IX
        ENDDO

        INDICE=INDICE+1
        DO JCH=1,KNCH
          DO JM=1,KNCO
            PZI(INDICE,JCH)=PZI(INDICE,JCH)+PZ0(IXY(JM),JCH)
          ENDDO
        ENDDO

      ENDDO
    ENDDO
  ENDIF

  DEALLOCATE(ZXCO)
  DEALLOCATE(ZYCO)

ENDDO

!     1.3 Normalisation.

ZCONO1=1.0_JPRB/REAL(KNCO,JPRB)
DO JCH=1,KNCH
  DO J=1,KNPT
    PZI(J,JCH)=PZI(J,JCH)*ZCONO1
  ENDDO
ENDDO

DEALLOCATE(IICO)
DEALLOCATE(INBLAT)

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EINTER1',1,ZHOOK_HANDLE)
END SUBROUTINE EINTER1
