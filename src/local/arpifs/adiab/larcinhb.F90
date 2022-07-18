SUBROUTINE LARCINHB(YDGEOMETRY,YDML_DYN,KST,KPROF,KASLB1,KL0H, &
 & PLSCAWH,PRSCAWH,PB1,PGWF,PGWF_NL)

! LARCINHB - semi-LAgrangian scheme (Interpolation)

! Purpose
! -------
!   LARCINHB - semi-Lagrangian scheme: origin point interpolations
!   for quantities situated at interface levels, i.e. vertical wind.
!   This subroutine is called only if LNHDYN=LGWADV=.TRUE. and LVFE_GW=F.

!   Abbreviation "vwv" stands for "vertical wind variable".

! Interface
! ---------
!   INPUT:
!     KST         - first elt of arrays where computations are performed.
!     KPROF       - depth of work.
!     KASLB1      - horizontal dimension of SL fields.
!     KL0H        - index of the four western points of the 16 points interpolation grid.
!     PLSCAWH     - linear weights (distances) for interpolations.
!     PRSCAWH     - non-linear weights for interpolations.
!     PB1         - "SLBUF1" buffer for interpolations.

!   OUTPUT:
!     PGWF        - Interpolated 3D quantity at O in the "gw" eqn.
!     PGWF_NL     - Interpolated 3D quantity at O in the "gw" eqn.

! Externals
! ---------

! Method
! ------

! Reference
! ---------

! Author
! ------
!      Aug-2002 C. SMITH, based on subroutine LARCINB.

! Modifications
! -------------
!   K. Yessad 07-03-2007: Remove useless (gw)_surf interpolations in NH+LGWADV.
!   28-Aug-2007 F. Vana    removing splines
!   30-Jun-2008 J. Masek   Dataflow for new SLHD interpolators.
!   K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!   F. Vana  15-Oct-2009: option NSPLTHOI
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   F. Vana  21-Feb-2011: horiz. turbulence, N[x]LAG=4
!   G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!   S. Malardel (Nov 2013): COMAD weights for SL interpolations
!   F. Vana  13-Feb-2014  Distinguish between heat and momentum SLHD.
!   K. Yessad (March 2017): simplify level numbering in interpolator.
!   K. Yessad (June 2017): alternate "cheap" FULL PC.
!   J. Vivoda (July 2018): LSETTLS with LPC_CHEAP.
! End Modifications
!------------------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK


!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KASLB1
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0H(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSCAWH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAWH%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSCAWH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAWH%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB1(KASLB1,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGWF(YDGEOMETRY%YRDIM%NPROMA,1:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGWF_NL(YDGEOMETRY%YRDIM%NPROMA,1:YDGEOMETRY%YRDIMV%NFLEVG)

!------------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IWK, IWKM
INTEGER(KIND=JPIM) :: INOWENOH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZGWF(YDGEOMETRY%YRDIM%NPROMA,1:YDGEOMETRY%YRDIMV%NFLEVG)
LOGICAL :: LLHOI2, LLMOM, LL_LIN
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------

#include "laitre_gmv.intfb.h"
#include "laitli.intfb.h"

!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LARCINHB',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDDYN=>YDML_DYN%YRDYN,YDPTRSLB1=>YDML_DYN%YRPTRSLB1)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & LQMHVD=>YDDYN%LQMHVD, LQMVD=>YDDYN%LQMVD, NSLDIMK=>YDDYN%NSLDIMK, &
 & LVWENO_SVD=>YDDYN%LVWENO_SVD, WENO_ALPHA_SVD=>YDDYN%WENO_ALPHA_SVD, &
 & NSPLTHOI=>YDDYN%NSPLTHOI, NSVDLAG=>YDDYN%NSVDLAG, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, NSITER=>YDDYN%NSITER, &
 & MSLB1VD0=>YDPTRSLB1%MSLB1VD0, MSLB1VD9=>YDPTRSLB1%MSLB1VD9, &
 & MSLB1VDF9=>YDPTRSLB1%MSLB1VDF9, NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & MSLB1VD9_NL=>YDPTRSLB1%MSLB1VD9_NL)
!------------------------------------------------------------------------------

! 0. PRELIMINARY INITIALISATIONS.

! Pointers for interpolation weights controlling the horizontal part of 3D turbulence
!  IWKM - momentum variables
IF (YDML_DYN%YRDYNA%L3DTURB) THEN
  IWKM=2
ELSE
  IWKM=1
ENDIF

LLMOM=.TRUE.

IF (YDML_DYN%YRDYNA%LPC_FULL.AND.YDML_DYN%YRDYNA%LPC_CHEAP2.AND.(NCURRENT_ITER<NSITER)) THEN
  LL_LIN=.TRUE.
ELSE
  LL_LIN=.FALSE.
ENDIF

! 1. ORIGIN POINT INTERPOLATIONS.

! * Vertical divergence equation.
CALL LAITRE_GMV(YDML_DYN,KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,IWKM,&
 & LL_LIN,YDML_DYN%YRDYNA%LCOMAD_SVD,YDML_DYN%YRDYNA%LSLHD_SVD,LVWENO_SVD,WENO_ALPHA_SVD,LLMOM,LQMVD,LQMHVD,INOWENOH,KL0H(1,1,0),PLSCAWH,PRSCAWH,&
 & PB1(1,MSLB1VD9),PGWF)

IF ((NSPLTHOI /= 0).AND.(NSVDLAG<4)) THEN
  ! --- part of the RHS requiring a "high order" interpolation
  !     not affected by a diffusion
  IF (NSPLTHOI == 1) THEN
    LLHOI2=.FALSE.
    IWK=NSLDIMK
  ELSE
    LLHOI2=YDML_DYN%YRDYNA%LSLHD_SVD
    IWK=IWKM
  ENDIF
  CALL LAITRE_GMV(YDML_DYN,KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,IWK,&
   & LL_LIN,YDML_DYN%YRDYNA%LCOMAD_SVD,LLHOI2,LVWENO_SVD,WENO_ALPHA_SVD,LLMOM,LQMVD,LQMHVD,INOWENOH,KL0H(1,1,0),PLSCAWH,PRSCAWH,&
   & PB1(1,MSLB1VDF9),ZGWF)
  PGWF(KST:KPROF,1:NFLEVG)=PGWF(KST:KPROF,1:NFLEVG)+ZGWF(KST:KPROF,1:NFLEVG)
ENDIF

IF (NSVDLAG >= 3) THEN
  IF (YDML_DYN%YRDYNA%LCOMAD_SVD) THEN
    CALL LAITLI(KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
     & PLSCAWH(1,1,YDML_DYN%YYTLSCAWH%M_WDLAMAD),PLSCAWH(1,1,YDML_DYN%YYTLSCAWH%M_WDLOMAD+1),&
     & KL0H(1,1,1),PLSCAWH(1,1,YDML_DYN%YYTLSCAWH%M_WDVERMAD),&
     & PB1(1,MSLB1VD0),ZGWF)
  ELSE
    CALL LAITLI(KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
     & PLSCAWH(1,1,YDML_DYN%YYTLSCAWH%M_WDLAT),PLSCAWH(1,1,YDML_DYN%YYTLSCAWH%M_WDLO+1),&
     & KL0H(1,1,1),PLSCAWH(1,1,YDML_DYN%YYTLSCAWH%M_WDVER),&
     & PB1(1,MSLB1VD0),ZGWF)
  ENDIF
  PGWF(KST:KPROF,1:NFLEVG)=PGWF(KST:KPROF,1:NFLEVG)+ZGWF(KST:KPROF,1:NFLEVG)
ENDIF

IF (YDML_DYN%YRDYNA%LSETTLS .AND. YDML_DYN%YRDYNA%LPC_CHEAP .AND. NSVDLAG == 3) THEN
  IF (YDML_DYN%YRDYNA%LCOMAD_SVD) THEN
    CALL LAITLI(KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
     & PLSCAWH(1,1,YDML_DYN%YYTLSCAWH%M_WDLAMAD),PLSCAWH(1,1,YDML_DYN%YYTLSCAWH%M_WDLOMAD+1),&
     & KL0H(1,1,1),PLSCAWH(1,1,YDML_DYN%YYTLSCAWH%M_WDVERMAD),&
     & PB1(1,MSLB1VD9_NL),ZGWF)
  ELSE
    CALL LAITLI(KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
     & PLSCAWH(1,1,YDML_DYN%YYTLSCAWH%M_WDLAT),PLSCAWH(1,1,YDML_DYN%YYTLSCAWH%M_WDLO+1),&
     & KL0H(1,1,1),PLSCAWH(1,1,YDML_DYN%YYTLSCAWH%M_WDVER),&
     & PB1(1,MSLB1VD9_NL),ZGWF)
  ENDIF
  PGWF_NL(KST:KPROF,1:NFLEVG)=ZGWF(KST:KPROF,1:NFLEVG)
ENDIF

!------------------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LARCINHB',1,ZHOOK_HANDLE)
END SUBROUTINE LARCINHB
