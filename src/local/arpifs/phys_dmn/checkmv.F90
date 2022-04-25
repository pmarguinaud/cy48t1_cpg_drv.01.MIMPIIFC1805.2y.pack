!OPTIONS XOPT(NOEVAL)
SUBROUTINE CHECKMV(KINDAT, YDCST, YDRIP,YDPHY0,YDPHY2,KIDIA,KFDIA,KLON,KLEV,KSTEP,PAPHI,PAPHIF &
& ,PAPRS,PAPRSF,PGELAM,PGEMU,PMU0,PLSM,PT,PQ,PTS)

!**** *CHECKMV * - Check magnitude of model variables, creates abort if out of bounds.

!**   Interface.
!     ----------
!        *CALL* *CHECKMV*

!-----------------------------------------------------------------------
! -   ARGUMENTS D'ENTREE.
! -   INPUT ARGUMENTS.
!     -------------------

! - DIMENSIONS.

! KIDIA, KFDIA : START/END OF HORIZONTAL LOOP  (IST,IEND IN *CPG*).
! KLON : HORIZONTAL DIMENSION                  (NPROMA IN *CPG*).
! KLEV : END OF VERTICAL LOOP AND VERTICAL DIMENSION(NFLEVG IN *CPG*).

! - 2D (0:KLEV) .

! PAPHI      : GEOPOTENTIAL ON HALF-LEVELS.
! PAPRS      : PRESSURE ON HALF-LEVELS.

! - 2D (1:KLEV) .

! PAPHIF     : GEOPOTENTIAL ON FULL LEVELS.
! PAPRSF     : PRESSURE ON FULL LEVELS.
! PT         : TEMPERATURE.
! PQ         : SPECIFIC HUMIDITY OF WATER VAPOUR.

! - 1D .
! PLSM       : LAND/SEA MASK.
! PTS        : SURFACE LAYER TEMPERATURE.
!-----------------------------------------------------------------------

! -   IMPLICIT ARGUMENTS.
!     ---------------------

!-----------------------------------------------------------------------

!     Externals.
!     ----------

!     Method.
!     -------

!     Author.
!     -------
!     2016-01-11, J.M. Piriou. 

!     Modifications.
!     --------------

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK
USE YOMCST   , ONLY : TCST
USE YOMPHY0  , ONLY : TPHY0
USE YOMPHY2  , ONLY : TPHY2
USE YOMLUN_IFSAUX, ONLY : NULERR
USE YOMRIP   , ONLY : TRIP

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KINDAT
TYPE(TPHY0)       ,INTENT(IN)    :: YDPHY0
TYPE(TPHY2)       ,INTENT(IN)    :: YDPHY2
TYPE (TCST), INTENT (IN) :: YDCST
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHI(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHIF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PTS(KLON)

REAL(KIND=JPRB) :: ZTOACI,ZDT,ZPROD,ZLOCST,ZCONRD,ZLON,ZLAT
INTEGER(KIND=JPIM) :: ILEV,ILON,JLEV,JLON,IVAR,IBIN

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CHECKMV',0,ZHOOK_HANDLE)
ASSOCIATE(RSTATI=>YDRIP%RSTATI,RHGMT=>YDRIP%RHGMT)
!     ------------------------------------------------------------------
ILEV=0
ILON=0
IVAR=0
DO JLEV=1,KLEV
  DO JLON=KIDIA,KFDIA
    !
    !-------------------------------------------------
    ! If T minus T_OACI is outside [ Tmin, Tmax ], ILEV will receive the corresponding level.
    !-------------------------------------------------
    !
    ZTOACI=YDCST%RTT+15._JPRB-0.0065*PAPHIF(JLON,JLEV)/YDCST%RG
    ZDT=PT(JLON,JLEV)-ZTOACI
    ZPROD=(ZDT-YDPHY0%GCHECKMV_TMIN)*(ZDT-YDPHY0%GCHECKMV_TMAX)
    IBIN=NINT(MAX(0._JPRB,SIGN(1._JPRB,ZPROD)))&
      & *NINT(MAX(0._JPRB,SIGN(1._JPRB,100000._JPRB-PAPHIF(JLON,JLEV))))
    ILEV=IBIN*JLEV+(1-IBIN)*ILEV
    ILON=IBIN*JLON+(1-IBIN)*ILON
    IVAR=IBIN*   1+(1-IBIN)*IVAR
    !
    !-------------------------------------------------
    ! If qv is outside [ qv_min, qv_max ], ILEV will receive the corresponding level.
    !-------------------------------------------------
    !
    ZPROD=(PQ(JLON,JLEV)-YDPHY0%GCHECKMV_QVMIN)*(PQ(JLON,JLEV)-YDPHY0%GCHECKMV_QVMAX)
    IBIN=NINT(MAX(0._JPRB,SIGN(1._JPRB,ZPROD)))&
      & *NINT(MAX(0._JPRB,SIGN(1._JPRB,100000._JPRB-PAPHIF(JLON,JLEV))))
    ILEV=IBIN*JLEV+(1-IBIN)*ILEV
    ILON=IBIN*JLON+(1-IBIN)*ILON
    IVAR=IBIN*   2+(1-IBIN)*IVAR
  ENDDO
ENDDO

IF(ILEV > 0) THEN
  !
  !-------------------------------------------------
  ! An out of bounds has occurred.
  !-------------------------------------------------
  !
  WRITE(NULERR,*) ' '
  IF(IVAR == 1) THEN
    WRITE(NULERR,*) 'CHECKMV/ERROR : T out of physical bounds ! '
    WRITE(NULERR,*) 'T(',ILON,',',ILEV,')=',PT(ILON,ILEV),' !!'
  ELSEIF(IVAR == 2) THEN
    WRITE(NULERR,*) 'CHECKMV/ERROR : qv out of physical bounds ! '
    WRITE(NULERR,*) 'qv(',ILON,',',ILEV,')=',PQ(ILON,ILEV),' !!'
  ENDIF
  WRITE(NULERR,*) ' '
  WRITE(NULERR,*) 'Time of error: '
  WRITE(NULERR,*) '  NINDAT=',KINDAT
  WRITE(NULERR,*) '  KSTEP=',KSTEP
  WRITE(NULERR,*) '  TSPHY (physical time step)=',YDPHY2%TSPHY
  WRITE(NULERR,*) '  RSTATI (time in seconds since model integration start)=',RSTATI
  WRITE(NULERR,*) '  RHGMT (GMT time in seconds between 0 and 86400)=',RHGMT
  WRITE(NULERR,*) '  GMT time (in hours)=',RHGMT/3600._JPRB
  !
  !-------------------------------------------------
  ! ZLOCST : local solar time in hours.
  ! The constant 0.261799387 rad is 15 degrees (the Earth rotates 15 degrees per hour).
  !-------------------------------------------------
  !
  ZLOCST=MODULO(RHGMT+PGELAM(ILON)/0.261799387_JPRB*3600._JPRB,86400._JPRB)/3600._JPRB
  WRITE(NULERR,*) '  Local solar time (in hours)=',ZLOCST
  !
  !-------------------------------------------------
  ! Location of error.
  !-------------------------------------------------
  !
  WRITE(NULERR,*) ' '
  WRITE(NULERR,*) 'Vertical location of error: '
  WRITE(NULERR,*) '  ILEV=',ILEV
  WRITE(NULERR,*) '  KLEV=',KLEV
  WRITE(NULERR,*) '  PAPHIF(',ILON,',',ILEV,')=',PAPHIF(ILON,ILEV)
  WRITE(NULERR,*) '  Geopotential at surface : PAPHI(',ILON,',KLEV)=',PAPHIF(ILON,KLEV)
  WRITE(NULERR,*) '  Land-sea mask PLSM=',PLSM(ILON)
  WRITE(NULERR,*) ' '
  WRITE(NULERR,*) 'Horizontal location of error: '
  WRITE(NULERR,*) '  ILON=',ILON
  ZCONRD=45._JPRB/ATAN(1._JPRB)
  ZLON=PGELAM(ILON)*ZCONRD
  IF(ZLON > 180._JPRB) ZLON=ZLON-360._JPRB
  ZLAT=ZCONRD*ASIN(PGEMU(ILON))
  WRITE(NULERR,*) '  Latitude, Longitude (in degrees) : ',ZLAT,' , ',ZLON
  WRITE(NULERR,*) '  PMU0 (cosine of Sun zenithal angle)=',PMU0(ILON)
  !
  !-------------------------------------------------
  ! Vertical profile.
  !-------------------------------------------------
  !
  WRITE(NULERR,*) ' '
  WRITE(NULERR,*) 'Vertical profile of T (K), qv (kg/kg), z (in m), p (hPa), level : '
  DO JLEV=1,KLEV
    WRITE(NULERR,FMT='(F8.2,F12.8,F10.1,F9.2,I5)') PT(ILON,JLEV),PQ(ILON,JLEV)&
      & ,PAPHIF(ILON,JLEV)/YDCST%RG,PAPRSF(ILON,JLEV)/100._JPRB,JLEV
  ENDDO
  !
  !-------------------------------------------------
  ! Surface temperature.
  !-------------------------------------------------
  !
  WRITE(NULERR,*) ' '
  WRITE(NULERR,*) 'Surface temperature : ',PTS(ILON)
  !
  !-------------------------------------------------
  ! Generates abort.
  !-------------------------------------------------
  !
  CALL ABOR1('CHECKMV : T or qv out of physical bounds !')
ENDIF


!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CHECKMV',1,ZHOOK_HANDLE)
END ASSOCIATE
END SUBROUTINE CHECKMV
