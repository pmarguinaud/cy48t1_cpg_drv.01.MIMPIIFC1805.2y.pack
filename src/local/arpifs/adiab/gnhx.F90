SUBROUTINE GNHX(&
 ! ----- INPUT ---------------------------------------------------------------
 & YDGEOMETRY,KPROMA,KSTART,KEND,POROG,POROGL,POROGM,&
 & PSP,PSPL,PSPM,PT,PTL,PTM,&
 & PQ,PQL,PQM,PL,PLL,PLM,PI,PIL,PIM,PR,PS,PG,PUF,PVF,&
 ! ----- OUTPUT --------------------------------------------------------------
 & PNHX,&
 ! ----- OPTIONAL INPUT ------------------------------------------------------
 & PSPD,PSPDL,PSPDM,PPIH)

! -----
! REMARKS:
!  - Variables PLL, PLM, PIL, PIM are not used currently, but they
!    might be needed in corrected version of GPRT.
!  - The structure of this routine must remain consistent with the one of CPG_GP.
!  - In the NHQE model, prognostic T is a modified temperature.
! -----

! GNHX - Diagnose NHX-term

! Purpose
! -------
!   Diagnose NHX-term

! Interface
! ---------
! * INPUT:
!   YDGEOMETRY   : structure containing all geometry.
!   KPROMA       : length of work
!   KSTART       : start of work
!   KEND         : end of work
!   POROG        : surface geopotential
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

! * OUTPUT:
!   PNHX         : NHX-term

! * OPTIONAL INPUT:
!   PSPD         : NH pressure departure (required for NHEE only)
!   PSPDL        : zonal gradient of NH pressure departure (required for NHEE only)
!   PSPDM        : meridional gradient of NH pressure departure (required for NHEE only)

! Externals
! ---------

! Method
! ------

! Reference
! ---------

! Author
! ------
!   08-Mar-2002 J. Masek (SHMI)

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (March 2009): correct false comments for LRWSDLG=T
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   K. Yessad (Dec 2011): Use GPHPRE.
!   K. Yessad (June 2017): Introduce NHQE model.
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   H Petithomme (Dec 2020): gphpre optimization, use of pointers
!
! End Modifications
!---------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDIMV      , ONLY : TDIMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : RD, RV
USE YOMCT0       , ONLY : LSPRT, LNHEE
USE YOMCVER      , ONLY : LVERTFE
USE YOMDYNA      , ONLY : LVEREGINT
USE YOMDYNA      , ONLY : NVDVAR, NPDVAR
USE INTDYN_MOD   , ONLY : YYTHW, YYTXYBDER, YYTXYB

! -----------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY),    INTENT(IN)       :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)       :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)       :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)       :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)       :: POROG  (KPROMA) 
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
REAL(KIND=JPRB)   ,INTENT(OUT)      :: PNHX   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN) :: PSPD   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN) :: PSPDL  (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN) :: PSPDM  (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB),OPTIONAL,TARGET,INTENT(IN) :: PPIH(KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF

REAL(KIND=JPRB) :: ZHOOK_HANDLE

REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZPIS(:),ZPIH(:,:) ! pressure prehyd and prehyds
REAL(KIND=JPRB) :: ZPISL (KPROMA)          ! zonal derivative of prehyds 
REAL(KIND=JPRB) :: ZPISM (KPROMA)          ! meridional derivative of prehyds

REAL(KIND=JPRB) :: ZXYBDER(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYBDER%NDIM) ! cf. PXYBDER in GPGRXYB
REAL(KIND=JPRB) :: ZXYB  (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)  ! contains "delta", "alpha"
REAL(KIND=JPRB),TARGET :: ZPIH0(KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) ! prehyd at half levels
REAL(KIND=JPRB) :: ZPIF  (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! prehyd at full levels
REAL(KIND=JPRB) :: ZCP   (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! specific heat of air cp
REAL(KIND=JPRB) :: ZR    (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! gas "constatnt" of air R
REAL(KIND=JPRB) :: ZKAPPA(KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! kappa = R/cp
REAL(KIND=JPRB) :: ZRT   (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! RT
REAL(KIND=JPRB) :: ZRTL  (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! zonal gradient of RT
REAL(KIND=JPRB) :: ZRTM  (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! meridional gradient of RT

REAL(KIND=JPRB) :: ZPHIHL(KPROMA, 0:YDGEOMETRY%YRDIMV%NFLEVG)  ! zonal comp of grad(gz) (half lay)
REAL(KIND=JPRB) :: ZPHIHM(KPROMA, 0:YDGEOMETRY%YRDIMV%NFLEVG)  ! merid comp of grad(gz) (half lay)
REAL(KIND=JPRB) :: ZPHIFL(KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)    ! zonal comp of grad(gz) (full lay)
REAL(KIND=JPRB) :: ZPHIFM(KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)    ! merid comp of grad(gz) (full lay)
REAL(KIND=JPRB) :: ZUVH(KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTHW%NDIM) ! hor wind at half levels

REAL(KIND=JPRB) :: ZNHPREF(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! "pre" at full levels.
REAL(KIND=JPRB) :: ZNHPPI (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! "pre/prehyd" at full levels.
REAL(KIND=JPRB) :: ZRNHPPI(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! "prehyd/pre" at full levels.
REAL(KIND=JPRB) :: ZNHPREL(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! zon comp of "grad pre" (full lev)
REAL(KIND=JPRB) :: ZNHPREM(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! mer comp of "grad pre" (full lev)
REAL(KIND=JPRB) :: ZLNNHPREFL(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) ! zon comp grad(ZNHPREF)/ZNHPREF
REAL(KIND=JPRB) :: ZLNNHPREFM(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) ! mer comp grad(ZNHPREF)/ZNHPREF
REAL(KIND=JPRB) :: ZQCHAL(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)     ! zonal comp grad(log(pre/prehyd))
REAL(KIND=JPRB) :: ZQCHAM(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)     ! merid comp grad(log(pre/prehyd))

! -----------------------------------------------------------------------------

#include "abor1.intfb.h"
#include "gnhgrpre.intfb.h"
#include "gnhpre.intfb.h"
#include "gpgrgeo.intfb.h"
#include "gpgrxyb.intfb.h"
#include "gphluv.intfb.h"
#include "gphlwi.intfb.h"
#include "gphpre.intfb.h"
#include "gprcp_qlirsg.intfb.h"
#include "gprt.intfb.h"
#include "gpxx.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHX',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,YDDIMV=>YDGEOMETRY%YRDIMV)

! -----------------------------------------------------------------------------

! -----
! computation of pressures
! -----

! half level pressure
IF (PRESENT(PPIH)) THEN
  ! prehyd and prehyds already computed, only compute derivatives
  ZPIH => PPIH(:,:)
  ZPIS => PPIH(:,NFLEVG)

  DO JROF=KSTART,KEND
    ZPISL(JROF)=ZPIS(JROF)*PSPL(JROF)
    ZPISM(JROF)=ZPIS(JROF)*PSPM(JROF)
  ENDDO
ELSE
  ! convert ln(prehyds) to prehyds (same for derivatives)
  ZPIH => ZPIH0(:,:)
  ZPIS => ZPIH0(:,NFLEVG)

  DO JROF=KSTART,KEND
    ZPIS(JROF)=EXP(PSP(JROF))
    ZPISL(JROF)=ZPIS(JROF)*PSPL(JROF)
    ZPISM(JROF)=ZPIS(JROF)*PSPM(JROF)
  ENDDO
ENDIF

CALL GPHPRE(KPROMA,NFLEVG,KSTART,KEND,YDVAB,ZPIH,PXYB=ZXYB,PRESF=ZPIF,LDELP=LVERTFE,&
 & LHSET=PRESENT(PPIH),LALPHA=.NOT.LVERTFE)

! additional auxiliary quantities (grad(delta) and grad(alpha) for ex.)
CALL GPGRXYB(KPROMA,KSTART,KEND,NFLEVG,.FALSE.,YDVAB,ZPISL,ZPISM,ZXYB,ZXYBDER)

! -----
! computation of R, cp, kappa
! -----

CALL GPRCP_QLIRSG(KPROMA,KSTART,KEND,NFLEVG,PQ,PI,PL,PR,PS,PG,ZCP,ZR,ZKAPPA)  

! -----
! computation of RT, grad(RT)
! -----

CALL GPRT(LSPRT,KPROMA,KSTART,KEND,NFLEVG,RD,RV,ZR,PT,PTL,PTM,PQL,PQM,ZRT,ZRTL,ZRTM)  

! -----
! computation of half level wind if relevant
! -----

IF (.NOT.LVERTFE) THEN
  ! compute interpolation weights
  CALL GPHLWI(YDDIMV,KPROMA,KSTART,KEND,ZXYB(1,1,YYTXYB%M_LNPR), &
   & ZXYB(1,1,YYTXYB%M_ALPH),ZUVH(1,1,YYTHW%M_WWI),LDVERINT=LVEREGINT)
  ! interpolate wind into half levels
  CALL GPHLUV(YDDIMV,KPROMA,KSTART,KEND,PUF,PVF,ZUVH)
ENDIF

! -----
! computation of pre/prehyd and some other "pressure departure" quantities for NHEE model
! computation of grad(gz), then NHX
! -----

IF (LNHEE) THEN
  IF (.NOT.(PRESENT(PSPD).AND.PRESENT(PSPDL).AND.PRESENT(PSPDM)) ) THEN
    CALL ABOR1(' GNHX: missing input PSPD, PSPDL, PSPDM !!!')
  ENDIF

  CALL GNHPRE(NPDVAR,YDGEOMETRY,KPROMA,NFLEVG,KSTART,KEND,PSPD,ZPIF,PNHPREF=ZNHPREF,&
   & PNHPPI=ZNHPPI,PRNHPPI=ZRNHPPI)

  CALL GNHGRPRE(NPDVAR,YDGEOMETRY,KPROMA,NFLEVG,KSTART,KEND,ZXYB(1,1,YYTXYB%M_RTGR),ZPISL,ZPISM,&
   & ZNHPREF,PSPDL,PSPDM,&
   & ZNHPREL,ZNHPREM,ZLNNHPREFL,ZLNNHPREFM,ZQCHAL,ZQCHAM)

  CALL GPGRGEO(YDGEOMETRY,KPROMA,KSTART,KEND,NFLEVG,&
   & ZRT,ZRTL,ZRTM,&
   & ZXYB(1,1,YYTXYB%M_LNPR),ZXYB(1,1,YYTXYB%M_ALPH),ZXYBDER,&
   & POROGL,POROGM,&
   & ZPHIFL,ZPHIFM,ZPHIHL,ZPHIHM,&
   & LDNHEE=LNHEE,PRNHPPI=ZRNHPPI,PQCHAL=ZQCHAL,PQCHAM=ZQCHAM)

  CALL GPXX(LVERTFE, YDGEOMETRY,NFLEVG,KPROMA,KSTART,KEND,ZPHIHL,ZPHIHM,ZPHIFL,ZPHIFM,ZXYB(1,1,YYTXYB%M_LNPR),&
   & ZRT,PUF,PVF,ZUVH(1,0,YYTHW%M_UH),ZUVH(1,0,YYTHW%M_VH),PNHX,PNHPPI=ZNHPPI)
ELSE
  ! NHQE, also valid for HYD.
  CALL GPGRGEO(YDGEOMETRY,KPROMA,KSTART,KEND,NFLEVG,&
   & ZRT,ZRTL,ZRTM,&
   & ZXYB(1,1,YYTXYB%M_LNPR),ZXYB(1,1,YYTXYB%M_ALPH),ZXYBDER,&
   & POROGL,POROGM,&
   & ZPHIFL,ZPHIFM,ZPHIHL,ZPHIHM)

  CALL GPXX(LVERTFE,YDGEOMETRY,NFLEVG,KPROMA,KSTART,KEND,ZPHIHL,ZPHIHM,ZPHIFL,ZPHIFM,ZXYB(1,1,YYTXYB%M_LNPR),&
   & ZRT,PUF,PVF,ZUVH(1,0,YYTHW%M_UH),ZUVH(1,0,YYTHW%M_VH),PNHX)
ENDIF

! -----------------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHX',1,ZHOOK_HANDLE)

END SUBROUTINE GNHX

