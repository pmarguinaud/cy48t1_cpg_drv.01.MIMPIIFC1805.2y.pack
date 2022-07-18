SUBROUTINE LARCINBTL(YDGEOMETRY,YDGMV,YGFL,YDML_DYN,KST,KPROF,KASLB1,LD2TLFF1,KL0,KLH0,&
 & PLSCAW,PRSCAW,PB1,PB15,KNOWENO,&
 & PGMVF,PDP,PGFLF,PCF,PUFZ,PVFZ,&
 & PUF5,PVF5,PUFZ5,PVFZ5,PLSCAW5,PRSCAW5)

!**** *LARCINBTL - semi-LAgrangian scheme:(Interpolation)
!                                            (tangent-linear version)

!     Purpose.
!     --------
!       Does the interpolations necessary in the semi-Lagrangian scheme.

!**   Interface.
!     ----------
!        *CALL* *LARCINBTL(......)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST     - first element where computations are performed.
!          KPROF   - depth of work.
!          KASLB1  - horizontal dimension of SL fields
!          LD2TLFF1- .T./.F.: Refined treatement of (2*Omega Vec r) at
!                    the origin point when there is t-dt (or t in SL2TL)
!                    physics / Other cases.
!          KL0     - index of the four western points
!                    of the 16 points interpolation grid.
!          KLH0    - second value of index of the four western points
!                    of the 16 points interpolation grid if needed.
!          PLSCAW  - linear weights (distances) for interpolations.
!          PRSCAW  - non-linear weights for interpolations.
!          PB1     - "SLBUF1" buffer for interpolations.
!          PB15    - "SLBUF1" buffer for interpolations (trajectory).
!          KNOWENO - boundary condition treatment for WENO

!        OUTPUT:
!          PGMVF   - Interpolated quantity at O for GMV.
!          PDP     - Interpolated 2D term at O for continuity equation.
!          PGFLF   - Interpolated quantity at O for GFL
!          PCF     - Interpolated 3D quantity at O in the continuity equation.
!          PUFZ and PVFZ for SL2TL scheme:
!          PUFZ    - Interpolated U-wind at O for L2TLFF=.T. option.
!          PVFZ    - Interpolated V-wind at O for L2TLFF=.T. option.
! ---------------------------- trajectory variables --------------------------
!          PUF5    - Trajectory for PGMVF(.,.,YGP%MU).
!          PVF5    - Trajectory for PGMVF(.,.,YGP%MV).
!          PUFZ5   - Trajectory for PUFZ.
!          PVFZ5   - Trajectory for PVFZ.
! ----------------------------------------------------------------------------

!        TRAJECTORY INPUT:
!          PLSCAW5 - linear weights (distances) for interpolations (trajectory).
!          PRSCAW5 - non-linear weights for interpolations (trajectory).

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!        Called by LAPINEBTL.
!        Calls LAITRE_GMV_TL,LAITLITL,LAIDDITL

!     Reference.
!     ----------

!     Author.
!     -------
!      C. Temperton (ECMWF)
!      Original: 99/07/12

!     Modifications.
!     --------------
!      Modified 01-08-30 by K. YESSAD: pruning and some other cleanings.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      D.Salmond     22-Dec-2003 Karim's Larcin cleaning
!      Modified 05-01-19 by C. Temperton: cleaning.
!      Modified 19-Mar-2007 K. Yessad: reorder according to LARCINB.
!      E.Holm      07/06/27: Option to advect without wind increments (LADV5)
!      F.Vana   30-Jun-2008: Weights driven interpolation
!      K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!      K. Yessad Nov 2009: rationalisation of dummy argument interfaces: PGMVF
!      K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!      G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      F. Vana  13-Feb-2014  Distinguish between heat and mometum SLHD
!      M. Diamantakis (Feb 2016): Introduce call to LAITRE_GFL_TL for GFL vars 
!      F. Vana 25-Feb-2019: Vertical quintic interpolation
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE YOMGMV             , ONLY : TGMV
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE YOM_YGFL           , ONLY : TYPE_GFLD


!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KASLB1
LOGICAL           ,INTENT(IN)    :: LD2TLFF1 
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLH0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB1(KASLB1,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB15(KASLB1,YDML_DYN%YRPTRSLB15%NFLDSLB15)
INTEGER(KIND=JPIM),INTENT(IN)    :: KNOWENO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGMVF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%YGP%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGFLF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YGFL%NUMFLDS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUFZ(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVFZ(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUFZ5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVFZ5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM)
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZGMVF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG), &
 & ZGMVF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDUM5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: JLEV, JROF,JGFL

LOGICAL         :: LLADV5, LLSLHD_P, LLMOM

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "laidditl.intfb.h"
#include "laitlitl.intfb.h"
#include "laitre_gmv_tl.intfb.h"
#include "laitre_gfl_tl.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LARCINBTL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDMP=>YDGEOMETRY%YRMP, YDDYN=>YDML_DYN%YRDYN,YDPTRSLB1=>YDML_DYN%YRPTRSLB1,YDPTRSLB15=>YDML_DYN%YRPTRSLB15, &
 & YDTLSCAW=>YDML_DYN%YYTLSCAW,YDTRSCAW=>YDML_DYN%YYTRSCAW)
ASSOCIATE(NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & LSLHDHEAT=>YDDYN%LSLHDHEAT, NTLAG=>YDDYN%NTLAG, NVLAG=>YDDYN%NVLAG, &
 & LVWENO_W=>YDDYN%LVWENO_W,LVWENO_T=>YDDYN%LVWENO_T,LVWENO_SP=>YDDYN%LVWENO_SP, &
 & WENO_ALPHA_W=>YDDYN%WENO_ALPHA_W,WENO_ALPHA_T=>YDDYN%WENO_ALPHA_T, &
 & WENO_ALPHA_SP=>YDDYN%WENO_ALPHA_SP, &
 & NWLAG=>YDDYN%NWLAG, &
 & YGP=>YDGMV%YGP, &
 & MSLB1C9=>YDPTRSLB1%MSLB1C9, MSLB1GFL9=>YDPTRSLB1%MSLB1GFL9, &
 & MSLB1SP9=>YDPTRSLB1%MSLB1SP9, MSLB1T0=>YDPTRSLB1%MSLB1T0, &
 & MSLB1T9=>YDPTRSLB1%MSLB1T9, MSLB1U0=>YDPTRSLB1%MSLB1U0, &
 & MSLB1U9=>YDPTRSLB1%MSLB1U9, MSLB1UR9=>YDPTRSLB1%MSLB1UR9, &
 & MSLB1V0=>YDPTRSLB1%MSLB1V0, MSLB1V9=>YDPTRSLB1%MSLB1V9, &
 & MSLB1VR9=>YDPTRSLB1%MSLB1VR9, NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & MSLB1C95=>YDPTRSLB15%MSLB1C95, MSLB1GFL95=>YDPTRSLB15%MSLB1GFL95, &
 & MSLB1SP95=>YDPTRSLB15%MSLB1SP95, MSLB1T05=>YDPTRSLB15%MSLB1T05, &
 & MSLB1T95=>YDPTRSLB15%MSLB1T95, MSLB1U05=>YDPTRSLB15%MSLB1U05, &
 & MSLB1U95=>YDPTRSLB15%MSLB1U95, MSLB1UR95=>YDPTRSLB15%MSLB1UR95, &
 & MSLB1V05=>YDPTRSLB15%MSLB1V05, MSLB1V95=>YDPTRSLB15%MSLB1V95, &
 & MSLB1VR95=>YDPTRSLB15%MSLB1VR95, NFLDSLB15=>YDPTRSLB15%NFLDSLB15)
!     ------------------------------------------------------------------

!*       2.    ORIGIN POINT INTERPOLATIONS FOR GMV VARIABLES.
!              ----------------------------------------------

! * U-Momentum equation:

LLMOM=.TRUE.
! --- part of the RHS requiring a "high order" interpolation:
LLADV5=.FALSE.
CALL LAITRE_GMV_TL(YDML_DYN,KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
 & YDML_DYN%YRDYNA%LSLHD_W,LVWENO_W,WENO_ALPHA_W,LLMOM,LLADV5,KNOWENO,KL0,PLSCAW,PRSCAW,PLSCAW5,PRSCAW5,&
 & PB1(1,MSLB1U9),PB15(1,MSLB1U95),PGMVF(1,1,YGP%MU),PUF5)

IF(NWLAG == 3) THEN
  ! --- save U in some cases:
  IF (LD2TLFF1) THEN
    ! Interpolate U at the origin point and use PUFZ
    ! to store the interpolated quantities.
    LLADV5=.FALSE.
    CALL LAITRE_GMV_TL(YDML_DYN,KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
     & YDML_DYN%YRDYNA%LSLHD_W,LVWENO_W,WENO_ALPHA_W,LLMOM,LLADV5,KNOWENO,KL0,PLSCAW,PRSCAW,PLSCAW5,PRSCAW5,&
     & PB1(1,MSLB1UR9),PB15(1,MSLB1UR95),PUFZ,PUFZ5)
  ELSE
    ! Save PGMVF(.,.,YGP%MU) in case of refined treatment of
    ! Coriolis term and lagged physics in 2TL scheme:
    PUFZ(KST:KPROF,1:NFLEVG)=PGMVF(KST:KPROF,1:NFLEVG,YGP%MU)
    PUFZ5(KST:KPROF,1:NFLEVG)=PUF5(KST:KPROF,1:NFLEVG)
  ENDIF

  ! --- part of the RHS requiring a trilinear interpolation:
  CALL LAITLITL(KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PLSCAW(1,1,YDTLSCAW%M_WDLO+1),KL0(1,1,1),PLSCAW(1,1,YDTLSCAW%M_WDVER),&
   & PLSCAW5(1,1,YDTLSCAW%M_WDLAT),PLSCAW5(1,1,YDTLSCAW%M_WDLO+1),PLSCAW5(1,1,YDTLSCAW%M_WDVER),&
   & PB1(1,MSLB1U0),PB15(1,MSLB1U05),ZGMVF,ZGMVF5)  
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      PGMVF(JROF,JLEV,YGP%MU)=PGMVF(JROF,JLEV,YGP%MU)+ZGMVF(JROF,JLEV)
      PUF5(JROF,JLEV)=PUF5(JROF,JLEV)+ZGMVF5(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

! * V-Momentum equation:

! --- part of the RHS requiring a "high order" interpolation:
LLADV5=.FALSE.
CALL LAITRE_GMV_TL(YDML_DYN,KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
 & YDML_DYN%YRDYNA%LSLHD_W,LVWENO_W,WENO_ALPHA_W,LLMOM,LLADV5,KNOWENO,KL0,PLSCAW,PRSCAW,PLSCAW5,PRSCAW5,&
 & PB1(1,MSLB1V9),PB15(1,MSLB1V95),PGMVF(1,1,YGP%MV),PVF5)

IF(NWLAG == 3) THEN
  ! --- save V in some cases:
  IF (LD2TLFF1) THEN
    ! Interpolate V at the origin point and use PVFZ
    ! to store the interpolated quantities.
    LLADV5=.FALSE.
    CALL LAITRE_GMV_TL(YDML_DYN,KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
     & YDML_DYN%YRDYNA%LSLHD_W,LVWENO_W,WENO_ALPHA_W,LLMOM,LLADV5,KNOWENO,KL0,PLSCAW,PRSCAW,PLSCAW5,PRSCAW5,&
     & PB1(1,MSLB1VR9),PB15(1,MSLB1VR95),PVFZ,PVFZ5)
  ELSE
    ! Save PGMVF(.,.,YGP%MV) in case of refined treatment of
    ! Coriolis term and lagged physics in 2TL scheme:
    PVFZ(KST:KPROF,1:NFLEVG)=PGMVF(KST:KPROF,1:NFLEVG,YGP%MV)
    PVFZ5(KST:KPROF,1:NFLEVG)=PVF5(KST:KPROF,1:NFLEVG)
  ENDIF

  ! --- part of the RHS requiring a trilinear interpolation:
  CALL LAITLITL(KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PLSCAW(1,1,YDTLSCAW%M_WDLO+1),KL0(1,1,1),PLSCAW(1,1,YDTLSCAW%M_WDVER),&
   & PLSCAW5(1,1,YDTLSCAW%M_WDLAT),PLSCAW5(1,1,YDTLSCAW%M_WDLO+1),PLSCAW5(1,1,YDTLSCAW%M_WDVER),&
   & PB1(1,MSLB1V0),PB15(1,MSLB1V05),ZGMVF,ZGMVF5)  
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      PGMVF(JROF,JLEV,YGP%MV)=PGMVF(JROF,JLEV,YGP%MV)+ZGMVF(JROF,JLEV)
      PVF5(JROF,JLEV)=PVF5(JROF,JLEV)+ZGMVF5(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

! * Temperature equation:

LLMOM=.FALSE. .OR. (.NOT.LSLHDHEAT)
! --- part of the RHS requiring a "high order" interpolation:
LLADV5=.FALSE.
CALL LAITRE_GMV_TL(YDML_DYN,KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
 & YDML_DYN%YRDYNA%LSLHD_T,LVWENO_T,WENO_ALPHA_T,LLMOM,LLADV5,KNOWENO,KL0,PLSCAW,PRSCAW,PLSCAW5,PRSCAW5,&
 & PB1(1,MSLB1T9),PB15(1,MSLB1T95),PGMVF(1,1,YGP%MT),ZDUM5)

IF(NTLAG == 3) THEN
  ! --- part of the RHS requiring a trilinear interpolation:
  CALL LAITLITL(KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PLSCAW(1,1,YDTLSCAW%M_WDLO+1),KL0(1,1,1),PLSCAW(1,1,YDTLSCAW%M_WDVER),&
   & PLSCAW5(1,1,YDTLSCAW%M_WDLAT),PLSCAW5(1,1,YDTLSCAW%M_WDLO+1),PLSCAW5(1,1,YDTLSCAW%M_WDVER),&
   & PB1(1,MSLB1T0),PB15(1,MSLB1T05),ZGMVF,ZDUM5)  
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      PGMVF(JROF,JLEV,YGP%MT)=PGMVF(JROF,JLEV,YGP%MT)+ZGMVF(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*       3.    ORIGIN POINT INTERPOLATIONS FOR GFL VARIABLES.
!              ----------------------------------------------

! * The interpolations are tri-dimensional ones.

DO JGFL=1,NUMFLDS
  IF(YCOMP(JGFL)%LADV) THEN
    IF (TRIM(YCOMP(JGFL)%CNAME) == 'TKE' ) THEN
      LLMOM=.TRUE.  ! the only momentum GFL variable so far
    ELSE
      LLMOM=.FALSE. .OR. (.NOT.LSLHDHEAT)
    ENDIF
    LLADV5=YCOMP(JGFL)%LADV5
    CALL LAITRE_GFL_TL(YDML_DYN,KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
     & YCOMP(JGFL)%CSLINT,YCOMP(JGFL)%WENO_ALPHA,LLMOM,LLADV5,KNOWENO,KL0,PLSCAW,PRSCAW,PLSCAW5,PRSCAW5,&
     & PB1(1,MSLB1GFL9+(YCOMP(JGFL)%MP_SL1-1)*(NFLEN-NFLSA+1)),&
     & PB15(1,MSLB1GFL95+(YCOMP(JGFL)%MP_SL1-1)*(NFLEN-NFLSA+1)),&
     & PGFLF(1,1,JGFL),ZDUM5)
  ENDIF
ENDDO

!     ------------------------------------------------------------------

!*       4.    ORIGIN POINT INTERPOLATIONS FOR GMVS VARIABLES.
!              -----------------------------------------------

LLMOM=.FALSE. .OR. (.NOT.LSLHDHEAT)
! The interpolations are tri-dimensional and 2D ones.

! * Continuity equation:

! --- 3D part of the RHS requiring a "high order" interpolation:
IF(NVLAG == 2) THEN
  LLADV5=.FALSE.
  LLSLHD_P=.FALSE.
  CALL LAITRE_GMV_TL(YDML_DYN,KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
   & LLSLHD_P,LVWENO_SP,WENO_ALPHA_SP,LLMOM,LLADV5,KNOWENO,KL0,PLSCAW,PRSCAW,PLSCAW5,PRSCAW5,&
   & PB1(1,MSLB1C9),PB15(1,MSLB1C95),PCF,ZDUM5)
ENDIF

! --- 3D part of the RHS requiring a trilinear interpolation:
IF(NVLAG == 3) THEN
  CALL LAITLITL(KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
   & PLSCAW(1,1,YDTLSCAW%M_WDLAT),PLSCAW(1,1,YDTLSCAW%M_WDLO+1),KL0(1,1,1),PLSCAW(1,1,YDTLSCAW%M_WDVER),&
   & PLSCAW5(1,1,YDTLSCAW%M_WDLAT),PLSCAW5(1,1,YDTLSCAW%M_WDLO+1),PLSCAW5(1,1,YDTLSCAW%M_WDVER),&
   & PB1(1,MSLB1C9),PB15(1,MSLB1C95),PCF,ZDUM5)  
ENDIF

! --- 2D part of the RHS requiring a "high order" interpolation:
CALL LAIDDITL(KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
 & PRSCAW(1,1,YDTRSCAW%M_WCLA(1)),PLSCAW(1,1,YDTLSCAW%M_WDLO),PRSCAW(1,1,YDTRSCAW%M_WCLO(1)),&
 & PRSCAW5(1,1,YDTRSCAW%M_WCLA(1)),PLSCAW5(1,1,YDTLSCAW%M_WDLO),PRSCAW5(1,1,YDTRSCAW%M_WCLO(1)),KLH0,&
 & PB1(1,MSLB1SP9),PB15(1,MSLB1SP95),PDP,ZDUM5)  

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LARCINBTL',1,ZHOOK_HANDLE)
END SUBROUTINE LARCINBTL

