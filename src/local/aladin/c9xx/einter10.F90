SUBROUTINE EINTER10(YDEGEO,PZ0,KLAT0,KLON0,KNCH,PZI,PZS,KLAT,KLON,KNPT,&
 & PDELX,PDELY,KPP,KNCO,KULOUT,PSEA0,PSEA1,PSEA2,PSEA3,PSPCM)  

!     PURPOSE.
!     --------

!       Interpolation of a field from a finer grid :
!        fields with sea/land contrast and major/minor averages
!       As EINTER7 but with 4 points interpolation instead of 12

!**   INTERFACE.
!     ----------
!     CALL EINTER10(...)

!        PZ0   = initial field
!        KLAT0 = number of latitudes  of the initial grid
!                  2  latitudes are added beyond each pole
!        KLON0 = number of longitudes of the initial grid
!                  4 longitudes are duplicated : 2 on each side
!        KNCH  = number of fields to interpolate
!        PZI   = interpolated fields                           (output)
!        PZS   = secondary interpolated fields                 (output)
!        KLAT  = number of latitudes  of the final grid
!        KLON  = number of longitudes of the final grid
!        KNPT  = number of points     of the final grid (KLAT*KLON)
!        PDEX,PDELY=size of final grid (m or rad.)
!        KPP   = size of the box
!        KNCO  = KPP*KPP
!        KULOUT= unit of output file
!        PSEA0 = fraction of sea on initial grid
!        PSEA1 = fraction of sea on  final  grid (result of 4 pt. inter.)
!        PSEA2 = fraction of sea on  final  grid (reference value)
!        PSEA3 = set equal to PSEA1                            (output)
!        PSPCM = threshold to define continents from the above reference

!     METHOD.
!     -------

!     The field is first calculated on a subgrid of the final grid.
!     Each value is interpolated from the 4 nearest points of the initial grid.
!     Then the values of the box (subgrid) are averaged (only for the points
!     which belong to the same "sea/land" category as the final point)
!       If there is not enough points in the neighbourhood to calculate
!       the secondary average (area < 10%), an alternate (back and forth)
!       research is done to find out a value in this category among the further
!       neighbours in the target grid.

!     EXTERNALS.
!     ----------

!       EGGPACK

!     AUTHORS.
!     --------
!     A. Dziedzic (27 OCT 95.)
!                 from INTRER10 (M. DEQUE 19 JULY 95.)
!                 and  EINTRER7 (J.F. GELEYN, M. DEQUE 1 FEB 91.)

!     MODIFICATION
!     ------------
!     M. Janousek 01-04-16: Replace EGGRVS by EGGPACK functions
!     JD. Gril 19-Nov-2004 Modifs for Mercator Rotated Tilted
!     --------------------------------------------------------------------

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

!     --------------------------------------------------------------------

IMPLICIT NONE

TYPE(TEGEO)       ,INTENT(IN)    :: YDEGEO
INTEGER(KIND=JPIM),INTENT(IN)    :: KLAT0 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON0 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNCH 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNPT 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZ0(KLON0,KLAT0,KNCH) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZI(KNPT,KNCH) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZS(KNPT,KNCH) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLAT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELY 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNCO 
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSEA0(KLON0,KLAT0) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSEA1(KNPT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSEA2(KNPT) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSEA3(KNPT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPCM 

!     --------------------------------------------------------------------

REAL(KIND=JPRB) :: ZWGT4(KNCO,4),ZXCO(KNCO*KLAT*KLON),ZYCO(KNCO*KLAT*KLON),&
 & ZZQ(KNPT,KNCH)  
INTEGER(KIND=JPIM) :: IXX(KNCO),IYY(KNCO)

INTEGER(KIND=JPIM) :: ICO, IDIF, II1, IJ1, IM, INDICE, IPAUT, ISMAS,&
 & IXM, IYM, IZ, J, JCH, JFL, JI, JJ, JK, JL, JM, JX, JY  

REAL(KIND=JPRB) :: ZCONO, ZCONO2, ZCX, ZCY, ZGX, ZGY, ZX, ZXX, ZY, ZYY, ZEPS

TYPE (LOLA)                :: YL_TLKRES, YL_TLREFC, YL_NWORIG
TYPE (LOLA) , ALLOCATABLE  :: YL_TLGRID_LOLA(:)
TYPE (PARAM_PROJ)          :: YL_TLMODDOM
TYPE (XY)   , ALLOCATABLE  :: YL_TLGRID_XY(:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     --------------------------------------------------------------------

#include "abor1.intfb.h"

!     --------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EINTER10',0,ZHOOK_HANDLE)

ASSOCIATE( &
 & ELON0=>YDEGEO%ELON0, ELAT0=>YDEGEO%ELAT0, ELAT1=>YDEGEO%ELAT1, &
 & ELON1=>YDEGEO%ELON1, ELONC=>YDEGEO%ELONC, ELATC=>YDEGEO%ELATC, LMRT=>YDEGEO%LMRT)
!     --------------------------------------------------------------------

!     1. Setting initial values to zero and extending the initial grid.

ZEPS=EPSILON(1.0_JPRB)*100.0_JPRB
DO JCH=1,KNCH
  DO J=1,KNPT
    PZI(J,JCH)=0.0_JPRB
    PZS(J,JCH)=0.0_JPRB
    ZZQ(J,JCH)=0.0_JPRB
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
 & WRITE(KULOUT,'('' EINTER10 : CAUTION !'',/,&
 & '' EITHER YOU ARE USING OLD E923, OR'',&
 & '' THERE IS A BUG IN THE SIZE OF DATASET (ARGUMENT 2 OR 3)'')')  

IF (YRCLI%LGLOBE) THEN
! global dataset : set extra latitudes and longitudes
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
ELSE
! non-global input grid: just duplicate border values into margins
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
ENDIF

!     2. Bulk summation.

!  Calculation of the coordinates of the points in the box

ICO=KNCO*KLAT*KLON
ALLOCATE(YL_TLGRID_XY(ICO))
ALLOCATE(YL_TLGRID_LOLA(ICO))
INDICE=0
ZCONO=1.0_JPRB/REAL(KPP,JPRB)
DO JJ=1,KLAT
  DO JI=1,KLON

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

!  Coordinate transform.

IF (LRPLANE) THEN
  IF (LMRT .AND. (ABS(ELAT0) >= ZEPS)) THEN
    WRITE(KULOUT,*) 'EINTER10: ELAT0=',ELAT0,&
    & ' MUST BE EQUAL ZERO IF LMRT IS TRUE!'
    CALL ABOR1('EINTER10: LMRT & ELAT0 INCONSISTENT')
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

!  Loop in the box
    DO JM=1,KNCO
! the interpolation values are localised.
      IM=JM+(JI-1+(JJ-1)*KLON)*KNCO
      ZX=MOD( ZGX/RPI/2.0_JPRB*(ZXCO(IM)+2.0_JPRB*RPI) , ZGX )
      ZX=ZX - II1 + 3.5_JPRB
      IXX(JM)=INT(ZX)
      ZXX=ZX-IXX(JM)
      ZCX=1.0_JPRB-ZXX
      ZY=ZGY/RPI*(RPI/2.0_JPRB+ZYCO(IM))
      ZY=ZY - IJ1 + 3.5_JPRB
      IYY(JM)=INT(ZY)
      ZYY=ZY-IYY(JM)
      ZCY=1.0_JPRB-ZYY

! the linear interpolation weights are set
      ZWGT4(JM,1)=ZXX*ZYY
      ZWGT4(JM,2)=ZCX*ZYY
      ZWGT4(JM,3)=ZXX*ZCY
      ZWGT4(JM,4)=ZCX*ZCY
    ENDDO

    INDICE=INDICE+1
    DO JM=1,KNCO
      IXM=IXX(JM)
      IYM=IYY(JM)
      DO JFL=1,KNCH
        PZI(INDICE,JFL)=PZI(INDICE,JFL)+&
         & (ZWGT4(JM, 1)*PZ0(IXM+1,IYM+1,JFL)&
         & +ZWGT4(JM, 2)*PZ0(IXM  ,IYM+1,JFL)&
         & +ZWGT4(JM, 3)*PZ0(IXM+1,IYM  ,JFL)&
         & +ZWGT4(JM, 4)*PZ0(IXM  ,IYM  ,JFL))  

        ZZQ(INDICE,JFL)=ZZQ(INDICE,JFL)+&
         & (ZWGT4(JM, 1)*PZ0(IXM+1,IYM+1,JFL)*PSEA0(IXM+1,IYM+1)&
         & +ZWGT4(JM, 2)*PZ0(IXM  ,IYM+1,JFL)*PSEA0(IXM  ,IYM+1)&
         & +ZWGT4(JM, 3)*PZ0(IXM+1,IYM  ,JFL)*PSEA0(IXM+1,IYM  )&
         & +ZWGT4(JM, 4)*PZ0(IXM  ,IYM  ,JFL)*PSEA0(IXM  ,IYM  ))  
      ENDDO
    ENDDO

  ENDDO
ENDDO

!     3. Normalisation.

ZCONO2=1.0_JPRB/REAL(KNCO,JPRB)

!  Securities
ISMAS=0
IDIF=0
DO J=1,KNPT
  IF ( (PSEA2(J) <= PSPCM.AND.PSEA1(J) > 0.9_JPRB) .OR.&
     & (PSEA2(J) > PSPCM.AND.PSEA1(J) < 0.1_JPRB) ) THEN  
    IDIF=IDIF+1
  ENDIF
  IF (PSEA2(J) <= PSPCM) ISMAS=ISMAS+1
ENDDO
WRITE(6,'('' NUMBER OF CONFLICT POINT BETWEEN '',&
 & ''PSEA1 (Int. Mask) AND PSEA2 (Final Mask): I=''&
 & ,I7)')IDIF  
WRITE(6,'('' ISMAS='',I6,'' KNPT='',I6,'' PSPCM='',F5.2)')ISMAS,KNPT,PSPCM

IF (ISMAS == 0.OR.ISMAS == KNPT) THEN
  WRITE(6,'('' SORRY : NEVER SECONDARY TYPE '')')
  WRITE(0,'('' EINTER10 : NEVER SECONDARY TYPE '')')
  WRITE(0,'(''     ---> CRASH AVOIDED BUT LOWER SECURITY LEVEL'')')
  DO J=1,KNPT
    WRITE(6,'('' PSEA2('',I6,'')='',F9.4)') J,PSEA2(J)
  ENDDO
! avoid crash --> comment ABORT -FT 02/2013-
!  CALL ABOR1('EINTER10 : SECURITIES, NEVER SECONDARY TYPE')
ENDIF

DO JFL=1,KNCH

  DO J=1,KNPT
    PZS(J,JFL)=-9999._JPRB
    IF (PSEA2(J) <= PSPCM) THEN
!           Land on the final grid
      IF (PSEA1(J) <= 0.9_JPRB) THEN
        PZI(J,JFL)=(PZI(J,JFL)-ZZQ(J,JFL))*ZCONO2/(1.0_JPRB-PSEA1(J))
      ELSE
        PZI(J,JFL)=-9999._JPRB
      ENDIF
      IF (PSEA1(J) >= 0.1_JPRB) THEN
        PZS(J,JFL)=ZZQ(J,JFL)*ZCONO2/PSEA1(J)
      ELSE
        PZS(J,JFL)=-9999._JPRB
      ENDIF
    ELSE
!           Sea on the final grid
      IF (PSEA1(J) <= 0.9_JPRB) THEN
        PZS(J,JFL)=(PZI(J,JFL)-ZZQ(J,JFL))*ZCONO2/(1.0_JPRB-PSEA1(J))
      ELSE
        PZS(J,JFL)=-9999._JPRB
      ENDIF
      IF (PSEA1(J) >= 0.1_JPRB) THEN
        PZI(J,JFL)=ZZQ(J,JFL)*ZCONO2/PSEA1(J)
      ELSE
        PZI(J,JFL)=-9999._JPRB
      ENDIF
    ENDIF
  ENDDO

  DO J=1,KNPT
    IF ((PZI(J,JFL) <= -9998).AND.(PSEA2(J) <= PSPCM)) THEN
!            Land on the final grid
      JY=J+1
      JY=MAX(1,MIN(KNPT,JY))
      IF (PSEA2(JY) <= PSPCM) THEN
        PZI(J,JFL)=PZI(JY,JFL)
      ELSE
        JY=J-1
        JY=MAX(1,MIN(KNPT,JY))
        IF (PSEA2(JY) <= PSPCM) THEN
          PZI(J,JFL)=PZI(JY,JFL)
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  DO J=KNPT,1,-1
    IF ((PZI(J,JFL) <= -9998).AND.(PSEA2(J) <= PSPCM)) THEN
!            Land on the final grid
      JY=J-1
      JY=MAX(1,MIN(KNPT,JY))
      IF ((PSEA2(JY) <= PSPCM).AND.(PZI(JY,JFL) > -9998)) THEN
        PZI(J,JFL)=PZI(JY,JFL)
      ELSE
        JY=J+1
        JY=MAX(1,MIN(KNPT,JY))
        IF ((PSEA2(JY) <= PSPCM).AND.(PZI(JY,JFL) > -9998)) THEN
          PZI(J,JFL)=PZI(JY,JFL)
        ENDIF
      ENDIF
    ENDIF
  ENDDO

  DO J=1,KNPT
    IF (PZS(J,JFL) <= -9998) THEN
!             Alternate research of a neighbour for the secondary average
      JX=1
      JY=J+1
      JK=MAX(1,MIN(KNPT,JY))
!               beginning of the test
      135 CONTINUE
      IF(((PSEA2(J) <= PSPCM.AND.PSEA2(JK) > PSPCM).OR.&
       & (PSEA2(J) > PSPCM.AND.PSEA2(JK) <= PSPCM)).AND.&
       & (PZI(JK,JFL) >= -9998._JPRB)) GO TO 137  

      IF (PZS(JK,JFL) >= -9998) GO TO 136
!               next point to test
      JX=SIGN(IABS(JX)+1,-JX)
      JY=JY+JX
      JK=MAX(1,MIN(KNPT,JY))
      GOTO 135
!                The point is of the same category and PZS <> -9999
      136 CONTINUE
      IF (PZS(J,JFL) <= -9998.) PZS(J,JFL)=PZS(JK,JFL)
      IF (PZI(J,JFL) <= -9998._JPRB) PZI(J,JFL)=PZI(JK,JFL)
      GO TO 138
!                The point is of the opposite category
      137 CONTINUE
      PZS(J,JFL)=PZI(JK,JFL)
      138 CONTINUE
    ENDIF
  ENDDO

ENDDO

DO J=1,KNPT
  PSEA3(J)=PSEA1(J)
ENDDO

WRITE(6,'(''----- fin de EINTER10 - OK!'')')

!     --------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EINTER10',1,ZHOOK_HANDLE)
END SUBROUTINE EINTER10
