SUBROUTINE GNHQE_NHX(&
 ! ----- INPUT ---------------------------------------------------------------
 & YDGEOMETRY,KPROMA,KSTART,KEND,POROG,POROGL,POROGM,&
 & PSP,PSPL,PSPM,PT,PTL,PTM,&
 & PQ,PQL,PQM,PL,PLL,PLM,PI,PIL,PIM,PR,PS,PG,PUF,PVF,PDIV,&
 ! ----- OUTPUT --------------------------------------------------------------
 & PNHX,PPIH)

! -----
! REMARKS:
!  - Variables PLL, PLM, PIL, PIM are not used currently, but they
!    might be needed in corrected version of GPRT.
!  - The structure of this routine must remain consistent with the one of CPG_GP_NHQE.
! -----

! GNHQE_NHX - Diagnose NHX-term in the NHQE model.

! Purpose
! -------
!   Diagnose NHX-term in the NHQE model.
!   NHX writes NHX=NHX_s+NHX_d
!   NHX_s has an expression close to the NHEE-model total NHX-term, and is linked to the
!    interaction between wind shear and orography gradient.
!   NHX_d is an additional term involving the horizontal divergence.
!   In the NHQE model, d4 = dver + NHX, but D3 = D + d4 - NHX_d = D + dver + NHX_s

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
!   PT           : Tt (modified temperature T)
!   PTL          : zonal gradient of Tt
!   PTM          : meridional gradient of Tt
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
!   PDIV         : horizontal divergence

! * OUTPUT:
!   PNHX         : NHX-term

! Externals
! ---------

! Method
! ------

! Reference
! ---------

! Author
! ------
!   K. Yessad and F. Voitus, after GNHX (July 2017).

! Modifications
! -------------
!   H Petithomme (Dec 2020): gphpre optimization, use of pointers
! End Modifications
!---------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDIMV      , ONLY : TDIMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : RD, RV
USE YOMCT0       , ONLY : LSPRT
USE YOMCVER      , ONLY : LVERTFE
USE INTDYN_MOD   , ONLY : YYTHW, YYTXYBDER, YYTXYB
USE YOMDYNA      , ONLY : LRUBC

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
REAL(KIND=JPRB)   ,INTENT(IN)       :: PDIV   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)      :: PNHX   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB),OPTIONAL,TARGET,INTENT(IN) :: PPIH(KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF, JLEV

REAL(KIND=JPRB) :: ZHOOK_HANDLE

REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZPIS(:),ZPIH(:,:) ! pressure prehyd and prehyds
REAL(KIND=JPRB) :: ZPISL (KPROMA)          ! zonal derivative of prehyds 
REAL(KIND=JPRB) :: ZPISM (KPROMA)          ! meridional derivative of prehyds

REAL(KIND=JPRB) :: ZXYBDER(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYBDER%NDIM) ! cf. PXYBDER in GPGRXYB
REAL(KIND=JPRB) :: ZXYB  (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)  ! contains "delta", "alpha"
REAL(KIND=JPRB),TARGET :: ZPIH0(KPROMA, 0:YDGEOMETRY%YRDIMV%NFLEVG) ! prehyd at half levels
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

REAL(KIND=JPRB) :: ZNHXS(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZNHXD(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 

! -----------------------------------------------------------------------------

#include "gpgrgeo.intfb.h"
#include "gpgrxyb.intfb.h"
#include "gphluv.intfb.h"
#include "gphlwi.intfb.h"
#include "gphpre.intfb.h"
#include "gprcp_qlirsg.intfb.h"
#include "gprt.intfb.h"
#include "gpxx.intfb.h"
#include "gnhqe_xxd.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHQE_NHX',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,YDDIMV=>YDGEOMETRY%YRDIMV)

! -----------------------------------------------------------------------------

! -----
! computation of pressures
! -----

! half level pressure
IF (PRESENT(PPIH)) THEN
  ! prehyds already computed, only compute derivatives
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

CALL GPHPRE(KPROMA,NFLEVG,KSTART,KEND,YDVAB,ZPIH,PXYB=ZXYB,PRESF=ZPIF,&
 & LHSET=PRESENT(PPIH),LALPHA=.NOT.LVERTFE)

! additional auxiliary quantities (grad(delta) and grad(alpha) for ex.)
CALL GPGRXYB(KPROMA,KSTART,KEND,NFLEVG,.FALSE.,YDVAB,ZPISL,ZPISM,ZXYB,ZXYBDER)

! -----
! computation of R, cp, kappa
! -----

CALL GPRCP_QLIRSG(KPROMA,KSTART,KEND,NFLEVG,PQ=PQ,PQI=PI,PQL=PL,PQR=PR,PQS=PS,PQG=PG,&
 & PCP=ZCP,PR=ZR,PKAP=ZKAPPA)  

! -----
! computation of RT, grad(RT)
! -----

CALL GPRT(LSPRT,KPROMA,KSTART,KEND,NFLEVG,RD,RV,ZR,PT,PTL,&
 & PTM,PQL,PQM,ZRT,ZRTL,ZRTM)  

! -----
! computation of half level wind if relevant
! -----

IF (.NOT.LVERTFE) THEN
  ! compute interpolation weights
  CALL GPHLWI(YDDIMV,KPROMA,KSTART,KEND,ZXYB(1,1,YYTXYB%M_LNPR),ZXYB(1,1,YYTXYB%M_ALPH),ZUVH(1,1,YYTHW%M_WWI))
  ! interpolate wind into half levels
  CALL GPHLUV(YDDIMV,KPROMA,KSTART,KEND,PUF,PVF,ZUVH)
ENDIF

! -----
! computation of grad(gz)
! -----
CALL GPGRGEO(YDGEOMETRY,KPROMA,KSTART,KEND,NFLEVG,&
 & ZRT,ZRTL,ZRTM,&
 & ZXYB(1,1,YYTXYB%M_LNPR),ZXYB(1,1,YYTXYB%M_ALPH),ZXYBDER,&
 & POROGL,POROGM,&
 & ZPHIFL,ZPHIFM,ZPHIHL,ZPHIHM)

! -----
! computation of NHX
! -----

CALL GPXX(YDGEOMETRY,NFLEVG,KPROMA,KSTART,KEND,ZPHIHL,ZPHIHM,ZPHIFL,ZPHIFM,ZXYB(1,1,YYTXYB%M_LNPR),&
  & ZRT,PUF,PVF,ZUVH(1,0,YYTHW%M_UH),ZUVH(1,0,YYTHW%M_VH),ZNHXS)

CALL GNHQE_XXD(LRUBC,LVERTFE,YDGEOMETRY,NFLEVG,KPROMA,KSTART,KEND,PUF,PVF,&
    & ZXYB(:,:,YYTXYB%M_LNPR), ZXYB(:,:,YYTXYB%M_RDELP), ZXYB(:,:,YYTXYB%M_ALPH), &
    & ZXYB(:,:,YYTXYB%M_RTGR), ZPISL,ZPISM,ZKAPPA,ZNHXD)

DO JLEV=1,NFLEVG
  DO JROF=KSTART,KEND
    PNHX(JROF,JLEV)=ZNHXS(JROF,JLEV)+ZNHXD(JROF,JLEV)
  ENDDO
ENDDO

! -----------------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHQE_NHX',1,ZHOOK_HANDLE)

END SUBROUTINE GNHQE_NHX

