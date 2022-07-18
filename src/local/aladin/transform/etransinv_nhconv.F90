!OCL NOEVAL
SUBROUTINE ETRANSINV_NHCONV(YDGEOMETRY,YGFL,YDSP,&
 ! ----- OUTPUT --------------------------------------------------------------
 & PU,PV,PDIV,PT,PTL,PTM,PSPD,PSPDL,PSPDM,PSVD,PSP,PSPL,PSPM,&
 & PQ,PQL,PQM,PL,PLL,PLM,PI,PIL,PIM,PR,PS,PG &
 &)

!* ETRANSINV_NHCONV - Conversion of NH variables (model to file and vice-versa)
!                     Inverse transforms for ALADIN

! Purpose
! -------

! Interface
! ---------
!  * OUTPUT:
!     PU     : U-wind
!     PV     : V-wind
!     PDIV   : horizontal divergence
!     PT     : temperature
!     PTL    : zonal comp of grad(temperature)
!     PTM    : merid comp of grad(temperature)
!     PSPD   : NH pressure departure variable "spd"
!     PSPDL  : zonal comp of grad(spd)
!     PSPDM  : merid comp of grad(spd)
!     PSVD   : NH vertical divergence variable "svd"
!     PSP    : log(prehyd)
!     PSPL   : zonal comp of grad(log(prehyd))
!     PSPM   : merid comp of grad(log(prehyd))
!     PQ     : moisture
!     PQL    : zonal comp of grad(moisture)
!     PQM    : merid comp of grad(moisture)
!     PL     : liquid water
!     PLL    : zonal comp of grad(liquid water)
!     PLM    : merid comp of grad(liquid water)
!     PI     : ice
!     PIL    : zonal comp of grad(ice)
!     PIM    : merid comp of grad(ice)
!     PR     : rain
!     PS     : snow
!     PG     : graupels

! Externals
! ---------

! Method
! ------
!   See documentation

! Reference
! ---------

! Author
! ------
!   27-Jan-2005 K. Yessad (after old part 2 of GNHPDVDCONV)

! Modifications
! -------------
!   K. Yessad (Dec 2009): bug correction.
!   K. Yessad (Jan 2011): merge ESPEREE_DER and ESPEREE.
!   O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Introduce NHQE model.
!   K. Yessad (Feb 2018): remove deep-layer formulations.
! End Modifications
!------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK

USE YOMCT0   , ONLY : LNHQE
USE YOMDYNA  , ONLY : YRDYNA
USE YOM_YGFL , ONLY : TYPE_GFLD
USE SPECTRAL_FIELDS_DATA, ONLY: SPECTRAL_FIELD

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN)  :: YDGEOMETRY
TYPE(TYPE_GFLD),INTENT(IN)  :: YGFL
TYPE(SPECTRAL_FIELD), INTENT(IN) :: YDSP
REAL(KIND=JPRB),INTENT(OUT) :: PU   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PV   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PDIV (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PT   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PTL  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PTM  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PSPD (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PSPDL(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PSPDM(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PSVD (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PSP  (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB),INTENT(OUT) :: PSPL (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB),INTENT(OUT) :: PSPM (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB),INTENT(OUT) :: PQ   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PQL  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PQM  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PL   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PLL  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PLM  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PI   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PIL  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PIM  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PR   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PS   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT) :: PG   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
!------------------------------------------------------------------------------

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------

#include "esperee.intfb.h"
#include "espeuv.intfb.h"

!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ETRANSINV_NHCONV',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE( &
 & NGPTOT=>YDGEM%NGPTOT, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLSUR=>YDDIMV%NFLSUR, &
 & YQ=>YGFL%YQ, YL=>YGFL%YL, YI=>YGFL%YI)
!------------------------------------------------------------------------------

!* 1. TRANSFORM SPECTRAL FIELDS TO GRID POINT SPACE
!     ---------------------------------------------

!  1.1 GMV variables

! * vorticity, divergence, wind
CALL ESPEUV(YDGEOMETRY,YDSP%VOR,YDSP%DIV,YDSP%UB,YDSP%VB,PU,PV,NFLSUR,NFLEVG,1)
IF (LNHQE) CALL ESPEREE(YDGEOMETRY,NFLSUR,NFLEVG,YDSP%DIV,PDIV)

! * temperature
IF ( YRDYNA%NVDVAR == 4 .OR. YRDYNA%NVDVAR == 5 ) THEN
  CALL ESPEREE(YDGEOMETRY,NFLSUR,NFLEVG,YDSP%T,PT,PREELL=PTL,PREELM=PTM)
ELSEIF ( YRDYNA%NVDVAR == 3 ) THEN
  CALL ESPEREE(YDGEOMETRY,NFLSUR,NFLEVG,YDSP%T,PT)
ENDIF

! * NH pressure departure variable
IF ( YRDYNA%NVDVAR == 4 .OR. YRDYNA%NVDVAR == 5 ) THEN
  CALL ESPEREE(YDGEOMETRY,NFLSUR,NFLEVG,YDSP%SPD,PSPD,PREELL=PSPDL,PREELM=PSPDM)
ELSEIF ( YRDYNA%NVDVAR == 3 ) THEN
  CALL ESPEREE(YDGEOMETRY,NFLSUR,NFLEVG,YDSP%SPD,PSPD)
ENDIF

! * NH vertical divergence variable
CALL ESPEREE(YDGEOMETRY,NFLSUR,NFLEVG,YDSP%SVD,PSVD)


!  1.2 GMVS variables

! * surface pressure
IF ( YRDYNA%NPDVAR == 2 ) THEN
  IF ( YRDYNA%NVDVAR == 4 .OR. YRDYNA%NVDVAR == 5 ) THEN
    CALL ESPEREE(YDGEOMETRY,1,1,YDSP%SP,PSP,PREELL=PSPL,PREELM=PSPM)
  ELSEIF ( YRDYNA%NVDVAR == 3 ) THEN
    CALL ESPEREE(YDGEOMETRY,1,1,YDSP%SP,PSP)
  ENDIF
ENDIF


!  1.3 GFL variables

! remark: this part has to be rewritten in the future with
!  "global treatment of GFL".

IF ( YRDYNA%NVDVAR == 4 .OR. YRDYNA%NVDVAR == 5 ) THEN

  ! * specific humidity
  IF (YQ%LSP) THEN
    CALL ESPEREE(YDGEOMETRY,NFLSUR,NFLEVG,YDSP%Q,PQ,PREELL=PQL,PREELM=PQM)
  ELSE
    PQ (1:NGPTOT,1:NFLEVG)=0.0_JPRB
    PQM(1:NGPTOT,1:NFLEVG)=0.0_JPRB
    PQL(1:NGPTOT,1:NFLEVG)=0.0_JPRB
  ENDIF

  ! * specific mass of liquid
  IF (YL%LSP) THEN
    CALL ESPEREE(YDGEOMETRY,NFLSUR,NFLEVG,YDSP%L,PL,PREELL=PLL,PREELM=PLM)
  ELSE
    PL (1:NGPTOT,1:NFLEVG)=0.0_JPRB
    PLM(1:NGPTOT,1:NFLEVG)=0.0_JPRB
    PLL(1:NGPTOT,1:NFLEVG)=0.0_JPRB
  ENDIF

  ! * specific mass of ice
  IF (YI%LSP) THEN
    CALL ESPEREE(YDGEOMETRY,NFLSUR,NFLEVG,YDSP%I,PI,PREELL=PIL,PREELM=PIM)
  ELSE
    PI (1:NGPTOT,1:NFLEVG)=0.0_JPRB
    PIM(1:NGPTOT,1:NFLEVG)=0.0_JPRB
    PIL(1:NGPTOT,1:NFLEVG)=0.0_JPRB
  ENDIF

  ! * Rain, graupel and snow for AROME
  PR(1:NGPTOT,1:NFLEVG) =0.0_JPRB
  PS(1:NGPTOT,1:NFLEVG) =0.0_JPRB
  PG(1:NGPTOT,1:NFLEVG) =0.0_JPRB

ENDIF ! NVDVAR

!------------------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ETRANSINV_NHCONV',1,ZHOOK_HANDLE)

END SUBROUTINE ETRANSINV_NHCONV

