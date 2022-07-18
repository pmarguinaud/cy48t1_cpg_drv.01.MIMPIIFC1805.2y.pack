SUBROUTINE LARCINB5(&
 ! --- INPUT ---------------------------------------------------
 & YDGEOMETRY,YDML_DYN,KST,KPROF,KASLB1,KNOWENO,KL0,PLSCAW5,PRSCAW5,PB15,&
 ! --- OUTPUT --------------------------------------------------
 & PUF,PVF,PUFZ,PVFZ)

!**** *LARCINB5    -  semi-LAgrangian scheme:(Interpolation)
!      (Stripped-down version of LARCINB, used by SL adjoint)

!     Purpose.
!     --------
!       Does the interpolations necessary in the semi-Lagrangian scheme.

!**   Interface.
!     ----------
!        *CALL* *LARCINB5(......)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST     - first element where computations are performed.
!          KPROF   - depth of work.
!          KASLB1  - horizontal dimension of SL fields
!          KNOWENO - specifies stencil for quintic vertical interpolation
!          KL0     - index of the four western points
!                    of the 16 points interpolation grid.
!          PLSCAW5 - linear weights (distances) for interpolations.
!          PRSCAW5 - non-linear weights for interpolations.
!          PB15    - "SLBUF1" buffer for interpolations (trajectory).

!        OUTPUT:
!          PUF     - Interpolated quantity at O in the U-wind eqn.
!          PVF     - Interpolated quantity at O in the V-wind eqn.
!          PUFZ    - Interpolated U-wind at O for L2TLFF=.T. option.
!          PVFZ    - Interpolated V-wind at O for L2TLFF=.T. option.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!      D. SALMOND, from LARCINB.
!      Original : DEC 2003.

!     Modifications.
!     --------------
!     11-Jan-2005 C. Temperton  Streamlining.
!     30-Jun-2008 J. Masek      New dataflow for modified LAITRI.
!     14-Aug-2008 F. Vana       Weights driven interpolation
!     K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!     K. Yessad Nov 2008: interpolation routines: merge QM with not-QM version.
!     K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!     G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!     F. Vana   March 2019:  make quintic interpolation also available
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK


!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KASLB1
INTEGER(KIND=JPIM),INTENT(IN)    :: KNOWENO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB15(KASLB1,YDML_DYN%YRPTRSLB15%NFLDSLB15)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUFZ(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVFZ(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZUF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZVF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: JLEV, JROF, IWKM

LOGICAL :: LL_LIN, LLCOMAD_W, LLMOM, LLQMHW

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "laitli.intfb.h"
#include "laitre_gmv.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LARCINB5',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
 & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDPTRSLB15=>YDML_DYN%YRPTRSLB15, &
 & YDTLSCAW=>YDML_DYN%YYTLSCAW,YDTRSCAW=>YDML_DYN%YYTRSCAW,YDDYN=>YDML_DYN%YRDYN)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & LQMW=>YDDYN%LQMW,LQMHW=>YDDYN%LQMHW,WENO_ALPHA_W=>YDDYN%WENO_ALPHA_W,&
 & LVWENO_W=>YDDYN%LVWENO_W,&
 & MSLB1U05=>YDPTRSLB15%MSLB1U05, MSLB1U95=>YDPTRSLB15%MSLB1U95, &
 & MSLB1V05=>YDPTRSLB15%MSLB1V05, MSLB1V95=>YDPTRSLB15%MSLB1V95, &
 & NFLDSLB15=>YDPTRSLB15%NFLDSLB15)
!     ------------------------------------------------------------------

!*       6.    ORIGIN POINT INTERPOLATIONS FOR U AND V.
!              ----------------------------------------

!  (U and V only; used in adjoint of semi-Lagrangian. Assumes NWLAG=3,
!  weights driven cubic/quintic and linear interpolations)

LL_LIN=.false.     ! Assume high order interpolation
LLCOMAD_W=.false.  ! Not coded so far in TL/AD
LLMOM=.true.       ! Momentum variables only
IWKM=1             ! No 3D turbulemnce available
LLQMHW=.false.

! High order interpolation
CALL LAITRE_GMV(YDML_DYN,KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,IWKM,&
 & LL_LIN,LLCOMAD_W,YDML_DYN%YRDYNA%LSLHD_W,LVWENO_W,WENO_ALPHA_W,LLMOM,LQMW,LLQMHW,KNOWENO,KL0,PLSCAW5,PRSCAW5,&
 & PB15(1,MSLB1U95),PUF)
CALL LAITRE_GMV(YDML_DYN,KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,IWKM,&
 & LL_LIN,LLCOMAD_W,YDML_DYN%YRDYNA%LSLHD_W,LVWENO_W,WENO_ALPHA_W,LLMOM,LQMW,LLQMHW,KNOWENO,KL0,PLSCAW5,PRSCAW5,&
 & PB15(1,MSLB1V95),PVF)
! * Save PUF,PVF in case of refined treatment of
!   Coriolis term and lagged physics in 2TL scheme:
PUFZ(KST:KPROF,1:NFLEVG)=PUF(KST:KPROF,1:NFLEVG)
PVFZ(KST:KPROF,1:NFLEVG)=PVF(KST:KPROF,1:NFLEVG)
! Linear interpolation
CALL LAITLI(KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
 & PLSCAW5(1,1,YDTLSCAW%M_WDLAT),PLSCAW5(1,1,YDTLSCAW%M_WDLO+1),KL0(1,1,1),PLSCAW5(1,1,YDTLSCAW%M_WDVER),&
 & PB15(1,MSLB1U05),ZUF)  
CALL LAITLI(KASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,NFLEN,&
 & PLSCAW5(1,1,YDTLSCAW%M_WDLAT),PLSCAW5(1,1,YDTLSCAW%M_WDLO+1),KL0(1,1,1),PLSCAW5(1,1,YDTLSCAW%M_WDVER),&
 & PB15(1,MSLB1V05),ZVF)  
DO JLEV=1,NFLEVG
  DO JROF=KST,KPROF
    PUF(JROF,JLEV)=PUF(JROF,JLEV) +ZUF(JROF,JLEV)
    PVF(JROF,JLEV)=PVF(JROF,JLEV) +ZVF(JROF,JLEV)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LARCINB5',1,ZHOOK_HANDLE)
END SUBROUTINE LARCINB5
