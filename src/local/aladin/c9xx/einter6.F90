SUBROUTINE EINTER6(YDEGEO,PZ0,KLAT0,KLON0,KNCH,PZI,KLAT,KLON,KNPT,&
 & PDELX,PDELY,KPP,KNCO,KULOUT,PSEA0,PSEA,PSPCM)  

!     PURPOSE.
!     --------

!       Interpolation of a field from a finer grid, with sea/land contrast

!**   INTERFACE.
!     ----------
!     CALL EINTER6(...)

!        PZ0   = initial field
!        KLAT0 = number of latitudes  of the initial grid
!                  2  latitudes are added beyond each pole
!        KLON0 = number of longitudes of the initial grid
!                  4 longitudes are duplicated : 2 on each side
!        KNCH  = number of fields to interpolate
!        PZI   = interpolated field                             (output)
!        KLAT  = number of latitudes  of the final grid
!        KLON  = number of longitudes of the final grid
!        KNPT  = number of points     of the final grid (KLAT*KLON)
!        PDEX,PDELY = size of final grid (m or rad.)
!        KPP   = size of the box
!        KNCO  = number of points in the box (KPP*KPP)
!        KULOUT= unit of output file
!        PSEA0 = fraction of sea on the initial grid
!        PSEA  = fraction of sea on the  final  grid (reference value)
!        PSPCM = threshold to define continents from the above reference

!     METHOD.
!     -------

!       The field is first calculated on a subgrid of the final grid.
!       Each value is interpolated from the 12 or 4 nearest points of the
!       initial grid.
!       Then the values of the box (subgrid) are averaged (only for the points
!       which belong to the same sea/land category as the final point)

!     EXTERNALS.
!     ----------

!       EGGPACK

!     AUTHORS.
!     --------
!      J.F. GELEYN  M. DEQUE  1 FEB 91.

!     MODIFICATIONS.
!     --------------
!     M. Janousek 01-04-16: Replace EGGRVS by EGGPACK functions
!     M.Hamrud      01-Oct-2003 CY28 Cleaning
!     JD. Gril 19-Nov-2004 Modifs for Mercator Rotated Tilted
!         ---------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : RPI      ,RA
USE YEMCLI   , ONLY : NILCLI1  ,NJLCLI1
USE YOMCLI   , ONLY : YRCLI
USE YEMGEO   , ONLY : TEGEO
USE YOMCT0   , ONLY : LRPLANE
USE EGGANGLES , ONLY : ANGLE_DOMAIN 
USE EGGPACK  , ONLY : LOLA,PARAM_PROJ,XY,XY_TO_LATLON,REF_DATAS, &
 & XY_NEW_TO_STD_ORIGIN

!         ---------------------------------------------------------------------

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
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSEA0(KLON0,KLAT0) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSEA(KNPT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPCM 

!         ---------------------------------------------------------------------

REAL(KIND=JPRB) :: ZWGT12(KNCO,12),ZXCO(KNCO*KLAT*KLON),ZYCO(KNCO*KLAT*KLON),&
 & ZXX(KNCO),ZYY(KNCO)  
INTEGER(KIND=JPIM) :: IXX(KNCO),IYY(KNCO)

INTEGER(KIND=JPIM) :: ICO, II1, IJ1, IM, INDICE, IPAUT, IXM, IYM,&
 & IZ, J, JCH, JFL, JI, JJ, JK, JL, JM  

REAL(KIND=JPRB) :: ZCONO, ZCONO2, ZCX, ZCY, ZEPS12, ZEPS4, ZGX,&
 & ZGY, ZWCO1, ZWCO2, ZWGT, ZX, ZXM, ZY, ZYM, ZEPS 

TYPE (LOLA)                :: YL_TLKRES, YL_TLREFC, YL_NWORIG
TYPE (LOLA) , ALLOCATABLE  :: YL_TLGRID_LOLA(:)
TYPE (PARAM_PROJ)          :: YL_TLMODDOM
TYPE (XY)   , ALLOCATABLE  :: YL_TLGRID_XY(:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!         ---------------------------------------------------------------------

#include "abor1.intfb.h"

!         ---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EINTER6',0,ZHOOK_HANDLE)

ASSOCIATE( &
 & ELON0=>YDEGEO%ELON0, ELAT0=>YDEGEO%ELAT0, ELAT1=>YDEGEO%ELAT1, &
 & ELON1=>YDEGEO%ELON1, ELONC=>YDEGEO%ELONC, ELATC=>YDEGEO%ELATC, LMRT=>YDEGEO%LMRT )
!         ---------------------------------------------------------------------

!     1. PRE-PROCESSING.
!        ---------------

!     1.1 Setting initial values.

ZEPS12=.25_JPRB
ZEPS4=1.E-8_JPRB
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
 & WRITE(KULOUT,'('' EINTER6 : CAUTION !'',/,&
 & '' EITHER YOU ARE USING OLD E923, OR'',&
 & '' THERE IS A BUG IN THE SIZE OF DATASET (ARGUMENT 2 OR 3)'')')  

!     1.2 Extending the initial grid according to LGLOBE.

IF (YRCLI%LGLOBE) THEN
! global dataset
  IPAUT=(KLON0-4)/2
  DO JFL=1,KNCH
    DO JI=3,IPAUT+2
      PZ0(JI      ,1      ,JFL)=PZ0(JI+IPAUT,4      ,JFL)
      PZ0(JI+IPAUT,1      ,JFL)=PZ0(JI      ,4      ,JFL)
      PZ0(JI      ,2      ,JFL)=PZ0(JI+IPAUT,3      ,JFL)
      PZ0(JI+IPAUT,2      ,JFL)=PZ0(JI      ,3      ,JFL)
      PZ0(JI      ,KLAT0  ,JFL)=PZ0(JI+IPAUT,KLAT0-3,JFL)
      PZ0(JI+IPAUT,KLAT0  ,JFL)=PZ0(JI      ,KLAT0-3,JFL)
      PZ0(JI      ,KLAT0-1,JFL)=PZ0(JI+IPAUT,KLAT0-2,JFL)
      PZ0(JI+IPAUT,KLAT0-1,JFL)=PZ0(JI      ,KLAT0-2,JFL)
    ENDDO
    DO JJ=1,KLAT0
      PZ0(1      ,JJ,JFL)=PZ0(KLON0-3,JJ,JFL)
      PZ0(2      ,JJ,JFL)=PZ0(KLON0-2,JJ,JFL)
      PZ0(KLON0-1,JJ,JFL)=PZ0(3      ,JJ,JFL)
      PZ0(KLON0  ,JJ,JFL)=PZ0(4      ,JJ,JFL)
    ENDDO
  ENDDO

  DO JI=3,IPAUT+2
    PSEA0(JI      ,1      )=PSEA0(JI+IPAUT,4      )
    PSEA0(JI+IPAUT,1      )=PSEA0(JI      ,4      )
    PSEA0(JI      ,2      )=PSEA0(JI+IPAUT,3      )
    PSEA0(JI+IPAUT,2      )=PSEA0(JI      ,3      )
    PSEA0(JI      ,KLAT0  )=PSEA0(JI+IPAUT,KLAT0-3)
    PSEA0(JI+IPAUT,KLAT0  )=PSEA0(JI      ,KLAT0-3)
    PSEA0(JI      ,KLAT0-1)=PSEA0(JI+IPAUT,KLAT0-2)
    PSEA0(JI+IPAUT,KLAT0-1)=PSEA0(JI      ,KLAT0-2)
  ENDDO
  DO JJ=1,KLAT0
    PSEA0(1      ,JJ)=PSEA0(KLON0-3,JJ)
    PSEA0(2      ,JJ)=PSEA0(KLON0-2,JJ)
    PSEA0(KLON0-1,JJ)=PSEA0(3      ,JJ)
    PSEA0(KLON0  ,JJ)=PSEA0(4      ,JJ)
  ENDDO

ELSE
! local dataset covering the whole ALADIN domain
! -> just duplicate border values into margins

  DO JFL=1,KNCH
    DO JI=3,KLON0-3
      PZ0(JI,1      ,JFL)=PZ0(JI,3      ,JFL)
      PZ0(JI,2      ,JFL)=PZ0(JI,3      ,JFL)
      PZ0(JI,KLAT0-1,JFL)=PZ0(JI,KLAT0-2,JFL)
      PZ0(JI,KLAT0  ,JFL)=PZ0(JI,KLAT0-2,JFL)
    ENDDO
    DO JJ=1,KLAT0
      PZ0(1      ,JJ,JFL)=PZ0(3      ,JJ,JFL)
      PZ0(2      ,JJ,JFL)=PZ0(3      ,JJ,JFL)
      PZ0(KLON0-1,JJ,JFL)=PZ0(KLON0-2,JJ,JFL)
      PZ0(KLON0  ,JJ,JFL)=PZ0(KLON0-2,JJ,JFL)
    ENDDO
  ENDDO

  DO JI=3,KLON0-3
    PSEA0(JI,1      )=PSEA0(JI,3      )
    PSEA0(JI,2      )=PSEA0(JI,3      )
    PSEA0(JI,KLAT0-1)=PSEA0(JI,KLAT0-2)
    PSEA0(JI,KLAT0  )=PSEA0(JI,KLAT0-2)
  ENDDO
  DO JJ=1,KLAT0
    PSEA0(1      ,JJ)=PSEA0(3      ,JJ)
    PSEA0(2      ,JJ)=PSEA0(3      ,JJ)
    PSEA0(KLON0-1,JJ)=PSEA0(KLON0-2,JJ)
    PSEA0(KLON0  ,JJ)=PSEA0(KLON0-2,JJ)
  ENDDO

ENDIF

!     ------------------------------------------------------------------

!     2. LOOP ON FINAL GRID.
!        -------------------

ZCONO=1.0_JPRB/REAL(KPP,JPRB)

ICO=KNCO*KLAT*KLON
ALLOCATE(YL_TLGRID_XY(ICO))
ALLOCATE(YL_TLGRID_LOLA(ICO))
INDICE=0
DO JJ=1,KLAT
  DO JI=1,KLON

!     2.1 Calculation of the coordinates of the points in the box.

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

!     2.2 Coordinate transform.

IF (LRPLANE) THEN
  IF (LMRT .AND. (ABS(ELAT0) >= ZEPS)) THEN
    WRITE(KULOUT,*) 'EINTER6: ELAT0=',ELAT0,&
    & ' MUST BE EQUAL ZERO IF LMRT IS TRUE!'
    CALL ABOR1('EINTER6: LMRT & ELAT0 INCONSISTENT')
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

!     2.3 Loop in the box : the interpolation values are localised.

DO JJ=1,KLAT
  DO JI=1,KLON

    DO JM=1,KNCO
      IM=JM+(JI-1+(JJ-1)*KLON)*KNCO
      ZX=MOD( ZGX/RPI/2.0_JPRB*(ZXCO(IM)+2.0_JPRB*RPI) , ZGX )
      ZX=3.5_JPRB+ZX-II1
      IXX(JM)=INT(ZX)
      ZXX(JM)=ZX-IXX(JM)
      ZY=ZGY/RPI*(RPI/2.0_JPRB+ZYCO(IM))
      ZY=3.5_JPRB+ZY-IJ1
      IYY(JM)=INT(ZY)
      ZYY(JM)=ZY-IYY(JM)
    ENDDO

!     2.4 Checking sea or land.

    INDICE=INDICE+1

    IF (PSEA(INDICE) <= PSPCM) THEN
! land
      ZWCO1=1.0_JPRB
      ZWCO2=-1.0_JPRB
    ELSE
! sea
      ZWCO1=0.0_JPRB
      ZWCO2=1.0_JPRB
    ENDIF

!     2.5 Loop in the box : the linear 12 pts interpolation weights are set.
!         (multiplied by 24 to save time)

    DO JM=1,KNCO
      ZXM=ZXX(JM)
      ZYM=ZYY(JM)
      IXM=IXX(JM)
      IYM=IYY(JM)
      ZCX=1.0_JPRB-ZXM
      ZCY=1.0_JPRB-ZYM
      ZWGT12(JM, 1)=3._JPRB*(8._JPRB*ZXM*ZYM+ZXM*ZCX*(ZXM-ZCX+2.0_JPRB*ZYM)&
       & +ZYM*ZCY*(ZYM-ZCY+2.0_JPRB*ZXM))*(ZWCO1+ZWCO2*PSEA0(IXM+1,IYM+1))  
      ZWGT12(JM, 2)=3._JPRB*(8._JPRB*ZCX*ZYM+ZCX*ZXM*(ZCX-ZXM+2.0_JPRB*ZYM)&
       & +ZYM*ZCY*(ZYM-ZCY+2.0_JPRB*ZCX))*(ZWCO1+ZWCO2*PSEA0(IXM  ,IYM+1))  
      ZWGT12(JM, 3)=3._JPRB*(8._JPRB*ZCX*ZCY+ZCX*ZXM*(ZCX-ZXM+2.0_JPRB*ZCY)&
       & +ZCY*ZYM*(ZCY-ZYM+2.0_JPRB*ZCX))*(ZWCO1+ZWCO2*PSEA0(IXM  ,IYM  ))  
      ZWGT12(JM, 4)=3._JPRB*(8._JPRB*ZXM*ZCY+ZXM*ZCX*(ZXM-ZCX+2.0_JPRB*ZCY)&
       & +ZCY*ZYM*(ZCY-ZYM+2.0_JPRB*ZXM))*(ZWCO1+ZWCO2*PSEA0(IXM+1,IYM  ))  
      ZWGT12(JM, 5)=ZXM*ZCX*(ZCX-ZXM-6._JPRB*ZYM)&
       & *(ZWCO1+ZWCO2*PSEA0(IXM+2,IYM+1))  
      ZWGT12(JM, 6)=ZYM*ZCY*(ZCY-ZYM-6._JPRB*ZXM)&
       & *(ZWCO1+ZWCO2*PSEA0(IXM+1,IYM+2))  
      ZWGT12(JM, 7)=ZYM*ZCY*(ZCY-ZYM-6._JPRB*ZCX)&
       & *(ZWCO1+ZWCO2*PSEA0(IXM  ,IYM+2))  
      ZWGT12(JM, 8)=ZCX*ZXM*(ZXM-ZCX-6._JPRB*ZYM)&
       & *(ZWCO1+ZWCO2*PSEA0(IXM-1,IYM+1))  
      ZWGT12(JM, 9)=ZCX*ZXM*(ZXM-ZCX-6._JPRB*ZCY)&
       & *(ZWCO1+ZWCO2*PSEA0(IXM-1,IYM  ))  
      ZWGT12(JM,10)=ZCY*ZYM*(ZYM-ZCY-6._JPRB*ZCX)&
       & *(ZWCO1+ZWCO2*PSEA0(IXM  ,IYM-1))  
      ZWGT12(JM,11)=ZCY*ZYM*(ZYM-ZCY-6._JPRB*ZXM)&
       & *(ZWCO1+ZWCO2*PSEA0(IXM+1,IYM-1))  
      ZWGT12(JM,12)=ZXM*ZCX*(ZCX-ZXM-6._JPRB*ZCY)&
       & *(ZWCO1+ZWCO2*PSEA0(IXM+2,IYM  ))  
    ENDDO

    ZWGT=0.0_JPRB
    DO JK=1,12
      DO JM=1,KNCO
        ZWGT=ZWGT+ZWGT12(JM,JK)
      ENDDO
    ENDDO
    ZWGT=ZWGT/REAL(24*KNCO,JPRB)
    IF (ZWGT >= ZEPS12) THEN
      ZCONO2=1.0_JPRB/REAL(24*KNCO,JPRB)/ZWGT
    ELSE

!     2.6 Loop in the box : the linear 4 pts interpolation weights are set.

      ZWGT=0.0_JPRB
      DO JM=1,KNCO
        ZXM=ZXX(JM)
        ZYM=ZYY(JM)
        IXM=IXX(JM)
        IYM=IYY(JM)
        ZCX=1.0_JPRB-ZXM
        ZCY=1.0_JPRB-ZYM
        ZWGT12(JM, 1)=ZXM*ZYM*(ZWCO1+ZWCO2*PSEA0(IXM+1,IYM+1))
        ZWGT12(JM, 2)=ZCX*ZYM*(ZWCO1+ZWCO2*PSEA0(IXM  ,IYM+1))
        ZWGT12(JM, 3)=ZCX*ZCY*(ZWCO1+ZWCO2*PSEA0(IXM  ,IYM  ))
        ZWGT12(JM, 4)=ZXM*ZCY*(ZWCO1+ZWCO2*PSEA0(IXM+1,IYM  ))
        ZWGT12(JM, 5)=0.0_JPRB
        ZWGT12(JM, 6)=0.0_JPRB
        ZWGT12(JM, 7)=0.0_JPRB
        ZWGT12(JM, 8)=0.0_JPRB
        ZWGT12(JM, 9)=0.0_JPRB
        ZWGT12(JM,10)=0.0_JPRB
        ZWGT12(JM,11)=0.0_JPRB
        ZWGT12(JM,12)=0.0_JPRB
        ZWGT=ZWGT+ZWGT12(JM,1)+ZWGT12(JM,2)+ZWGT12(JM,3)+ZWGT12(JM,4)
      ENDDO
      ZWGT=ZWGT/REAL(KNCO,JPRB)

      IF (ZWGT < ZEPS4) THEN
! the 4 pts interpolation weights are recalculated by inverting the mask.
        DO JM=1,KNCO
          ZXM=ZXX(JM)
          ZYM=ZYY(JM)
          ZCX=1.0_JPRB-ZXM
          ZCY=1.0_JPRB-ZYM
          ZWGT12(JM,1)=ZXM*ZYM-ZWGT12(JM,1)
          ZWGT12(JM,2)=ZCX*ZYM-ZWGT12(JM,2)
          ZWGT12(JM,3)=ZCX*ZCY-ZWGT12(JM,3)
          ZWGT12(JM,4)=ZXM*ZCY-ZWGT12(JM,4)
        ENDDO
        ZWGT=1.0_JPRB-ZWGT
      ENDIF
      ZCONO2=1.0_JPRB/REAL(KNCO,JPRB)/ZWGT

    ENDIF

!     2.7 Interpolation.

    DO JM=1,KNCO
      DO JFL=1,KNCH
        IXM=IXX(JM)
        IYM=IYY(JM)
        PZI(INDICE,JFL)=PZI(INDICE,JFL)&
         & +ZWGT12(JM, 1)*PZ0(IXM+1,IYM+1,JFL)&
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
         & +ZWGT12(JM,12)*PZ0(IXM+2,IYM  ,JFL)  
      ENDDO
    ENDDO
    DO JFL=1,KNCH
      PZI(INDICE,JFL)=PZI(INDICE,JFL)*ZCONO2
    ENDDO

  ENDDO
ENDDO

!         ---------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EINTER6',1,ZHOOK_HANDLE)
END SUBROUTINE EINTER6
