SUBROUTINE GNHEE_SVDINCR13_CALLER(&
 ! ----- INPUT ---------------------------------------------------------------
 & YDGEOMETRY,KPROMA,KSTART,KEND,POROGL,POROGM,&
 & PSP,PSPL,PSPM,PT,PTL,PTM,&
 & PQ,PQL,PQM,PL,PLL,PLM,PI,PIL,PIM,PR,PS,PG,PUF,PVF,PSPD,&
 ! ----- OUTPUT --------------------------------------------------------------
 & PSVDINCR13)

! -----
! REMARKS:
!  - Variables PLL, PLM, PIL, PIM are not used currently, but they
!    might be needed in corrected version of GPRT.
!  - The structure of this routine must remain consistent with the one of CPG_GP_NHEE.
! -----

! GNHEE_SVDINCR13_CALLER - caller for GNHEE_SVDINCR13
!                          compute d13 - dver in NHEE model
!                          This routine is called for NVDVAR=13 or 14

! Purpose
! -------
!  See comments in GNHEE_SVDINCR13.

! Remarks:
! - structure of GNHEE_SVDINCR13_CALLER follows GNHX one.
! - to avoid inconsistencies in parts doing conversions between d13 and [gW],
!   the increment (d13 - dver) is defined with the same value of R as 'dver',
!   according to the value of L_RDRY_VD.

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
!   PSPL         : zonal gradient of ln(prehyds)
!   PSPM         : meridional gradient of ln(prehyds)
!   PT           : temperature T
!   PTL          : zonal gradient of temperature
!   PTM          : meridional gradient of temperature
!   PQ           : specific humidity q
!   PQL          : zonal gradient of q
!   PQM          : meridional gradient of q
!   PL           : specific mass of liquid ql
!   PLL          : zonal gradient of ql
!   PLM          : meridional gradient of ql
!   PI           : specific mass of ice qi
!   PIL          : zonal gradient of qi
!   PIM          : meridional gradient of qi
!   PR           : specific mass of rain qr
!   PS           : specific mass of snow qs
!   PG           : specific mass of graupel qg
!   PUF          : U-wind at full levels
!   PVF          : V-wind at full levels
!   PSPD         : NH pressure departure variable

! * OUTPUT:
!   PSVDINCR13   : (d13 - dver) at full levels.

! Externals
! ---------

! Method
! ------

! Reference
! ---------

! Author
! ------
!   K. Yessad (August 2018), after GNHX.

! Modifications
! -------------
! End Modifications
!---------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDIMV      , ONLY : TDIMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : RD, RV
USE YOMCT0       , ONLY : LSPRT
USE YOMCVER      , ONLY : LVERTFE, LVFE_GW
USE YOMDYNA      , ONLY : L_RDRY_VD, LVEREGINT
USE YOMDYNA      , ONLY : NVDVAR, NPDVAR
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
REAL(KIND=JPRB)   ,INTENT(IN)       :: PSPL   (KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PSPM   (KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PT     (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PTL    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PTM    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PQ     (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PQL    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PQM    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PL     (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PLL    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)       :: PLM    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)       :: PI     (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PIL    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)       :: PIM    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)       :: PR     (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PS     (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PG     (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PUF    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PVF    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)      :: PSVDINCR13(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN) :: PSPD   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF, JLEV

REAL(KIND=JPRB) :: ZHOOK_HANDLE

REAL(KIND=JPRB) :: ZPIS  (KPROMA)          ! surface pressure prehyds
REAL(KIND=JPRB) :: ZPISL (KPROMA)          ! zonal derivative of prehyds 
REAL(KIND=JPRB) :: ZPISM (KPROMA)          ! meridional derivative of prehyds

REAL(KIND=JPRB) :: ZXYB  (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)  ! contains "delta", "alpha"
REAL(KIND=JPRB) :: ZPIH  (KPROMA, 0:YDGEOMETRY%YRDIMV%NFLEVG)! prehyd at half levels
REAL(KIND=JPRB) :: ZPIF  (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! prehyd at full levels
REAL(KIND=JPRB) :: ZR    (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! gas "constatnt" of air R
REAL(KIND=JPRB) :: ZRT   (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! RT
REAL(KIND=JPRB) :: ZRTL  (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! zonal gradient of RT
REAL(KIND=JPRB) :: ZRTM  (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! meridional gradient of RT

REAL(KIND=JPRB) :: ZUVH(KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTHW%NDIM) ! hor wind at half levels

REAL(KIND=JPRB) :: ZNHPPI (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! "pre/prehyd" at full levels.

! -----------------------------------------------------------------------------

#include "abor1.intfb.h"
#include "gnhpre.intfb.h"
#include "gphluv.intfb.h"
#include "gphlwi.intfb.h"
#include "gphpre.intfb.h"
#include "gprcp_qlirsg.intfb.h"
#include "gprt.intfb.h"
#include "gnhee_svdincr13.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHEE_SVDINCR13_CALLER',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,YDDIMV=>YDGEOMETRY%YRDIMV)

! -----------------------------------------------------------------------------

! * convert ln(prehyds) to prehyds (same for derivatives)
DO JROF=KSTART,KEND
  ZPIS (JROF)=EXP(PSP(JROF))
  ZPISL(JROF)=ZPIS(JROF)*PSPL(JROF)
  ZPISM(JROF)=ZPIS(JROF)*PSPM(JROF)
ENDDO

! * half level pressure and ZXYB
DO JROF=KSTART,KEND
  ZPIH(JROF,NFLEVG)=ZPIS(JROF)
ENDDO
CALL GPHPRE(KPROMA,NFLEVG,KSTART,KEND,YDVAB,ZPIH,PXYB=ZXYB,PRESF=ZPIF)

! * computation of RT (R is the one used in definition of 'dver').
IF (L_RDRY_VD) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KSTART,KEND
      ZRT(JROF,JLEV)=RD*PT(JROF,JLEV)
    ENDDO
  ENDDO
ELSE
  CALL GPRCP_QLIRSG(KPROMA,KSTART,KEND,NFLEVG,PQ=PQ,PQI=PI,PQL=PL,PQR=PR,PQS=PS,PQG=PG,PR=ZR)  
  CALL GPRT(LSPRT,KPROMA,KSTART,KEND,NFLEVG,RD,RV,ZR,PT,PTL,PTM,PQL,PQM,ZRT,ZRTL,ZRTM)  
ENDIF

! * computation of pre/prehyd.
CALL GNHPRE(NPDVAR,YDGEOMETRY,KPROMA,NFLEVG,KSTART,KEND,PSPD,ZPIF,PNHPPI=ZNHPPI)

! * computation of half level wind if relevant
IF (LVERTFE .AND. LVFE_GW) THEN
  ! case not yet coded.
  CALL ABOR1(' GNHEE_SVDINCR13_CALLER: LVFE_GW not yet coded! ')
  ! ZUVH is useless, set to 0 for safety.
  ZUVH(KSTART:KEND,0:NFLEVG,1:YYTHW%NDIM)=0.0_JPRB
ELSE
  ! compute interpolation weights
  CALL GPHLWI(YDDIMV,KPROMA,KSTART,KEND,ZXYB(1,1,YYTXYB%M_LNPR),ZXYB(1,1,YYTXYB%M_ALPH),&
            & ZUVH(1,1,YYTHW%M_WWI),LDVERINT=LVEREGINT)
  ! interpolate wind into half levels
  CALL GPHLUV(YDDIMV,KPROMA,KSTART,KEND,PUF,PVF,ZUVH)
ENDIF

! * computation of (d13 - dver).
CALL GNHEE_SVDINCR13(LVERTFE,LVFE_GW,YDGEOMETRY,NFLEVG,KPROMA,KSTART,KEND,POROGL,POROGM,&
 & ZUVH(1,0,YYTHW%M_UH),ZUVH(1,0,YYTHW%M_VH),&
 & PLNPR=ZXYB(1,1,YYTXYB%M_LNPR),PRT=ZRT,PNHPPI=ZNHPPI,&
 & PSVDINCR13=PSVDINCR13)

! -----------------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHEE_SVDINCR13_CALLER',1,ZHOOK_HANDLE)

END SUBROUTINE GNHEE_SVDINCR13_CALLER

