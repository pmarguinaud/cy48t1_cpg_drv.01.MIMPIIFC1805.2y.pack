SUBROUTINE GNHY(&
 ! ----- INPUT ---------------------------------------------------------------
 & YDGEOMETRY,KPROMA,KSTART,KEND,POROGL,POROGM,&
 & PSP,PUF,PVF,&
 ! ----- OUTPUT --------------------------------------------------------------
 & PNHY)

! -----
! REMARKS:
!  -
! -----

! GNHY - Diagnose NHY-term

! Purpose
! -------
!   Diagnose NHY-term

! Interface
! ---------
! * INPUT:
!   YDGEOMETRY   : structure containing all geometry.
!   KPROMA       : length of work
!   KSTART       : start of work
!   KEND         : end of work
!   POROGL       : zonal gradient of surface geopotential
!   POROGM       : meridional gradient of surface geopotential
!   PSP          : ln(prehyds)
!   PUF          : U-wind at full levels
!   PVF          : V-wind at full levels

! * OUTPUT:
!   PNHY         : NHY-term

! Externals
! ---------

! Method
! ------

! Reference
! ---------

! Author
! ------
!   28-Nov-2019 F. Voitus (CNRM/GMAP/ALGO)

!---------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDIMV      , ONLY : TDIMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCT0       , ONLY : LNHEE
USE YOMDYNA      , ONLY : YRDYNA
USE YOMCVER      , ONLY : LVERTFE, LVFE_GW
USE INTDYN_MOD   , ONLY : YYTHW, YYTXYB

! -----------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY),    INTENT(IN)       :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)       :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)       :: KSTART
INTEGER(KIND=JPIM),INTENT(IN)       :: KEND
REAL(KIND=JPRB)   ,INTENT(IN)       :: POROGL (KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)       :: POROGM (KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)       :: PSP    (KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)       :: PUF    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)       :: PVF    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)      :: PNHY   (KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF

REAL(KIND=JPRB) :: ZXYB  (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)  ! contains "delta", "alpha"
REAL(KIND=JPRB) :: ZPREH  (KPROMA, 0:YDGEOMETRY%YRDIMV%NFLEVG)           ! prehyd at half levels
REAL(KIND=JPRB) :: ZUVH(KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTHW%NDIM)    ! hor wind at half levels

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

#include "gphluv.intfb.h"
#include "gphlwi.intfb.h"
#include "gphpre.intfb.h"
#include "gpyy.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHY',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,YDDIMV=>YDGEOMETRY%YRDIMV)

! -----------------------------------------------------------------------------

! -----
! computation of alpha, delta
! -----

! half level pressure
DO JROF=KSTART,KEND
  ZPREH(JROF,NFLEVG)=EXP(PSP(JROF))
ENDDO
CALL GPHPRE(KPROMA,NFLEVG,KSTART,KEND,YDVAB,ZPREH,PXYB=ZXYB,LDELP=.FALSE.,LRTGR=.FALSE.,&
 & LRPP=.FALSE.)

! -----
! computation of half level wind if relevant
! -----

! compute interpolation weights
CALL GPHLWI(YDDIMV,KPROMA,KSTART,KEND,ZXYB(1,1,YYTXYB%M_LNPR),&
          & ZXYB(1,1,YYTXYB%M_ALPH),ZUVH(1,1,YYTHW%M_WWI),LDVERINT=YRDYNA%LVEREGINT)
! interpolate wind into half levels
CALL GPHLUV(YDDIMV,KPROMA,KSTART,KEND,PUF,PVF,ZUVH)

! -----
! computation of NHY-term
! -----
CALL GPYY(LVERTFE, LVFE_GW, YDGEOMETRY,NFLEVG,KPROMA,KSTART,KEND,POROGL,POROGM,&
 & ZUVH(1,0,YYTHW%M_UH),ZUVH(1,0,YYTHW%M_VH),PNHY)

! -----------------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHY',1,ZHOOK_HANDLE)

END SUBROUTINE GNHY

