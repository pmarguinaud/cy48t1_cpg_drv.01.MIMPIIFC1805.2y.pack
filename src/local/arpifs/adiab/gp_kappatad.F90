SUBROUTINE GP_KAPPATAD(&
 & YDGEOMETRY, YDDYN,KPROMA,KSTART,KPROF,KFLEV,&
 & PDTS2,PTT05,PTT0L5,PTT0M5,PSPT05,PSPT0L5,PSPT0M5,&
 & PTT0,PTT0L,PTT0M,PSPT0,PSPT0L,PSPT0M,PSLHDA,PSLHDD0,&
 & PKAPPA)

!**** *GP_KAPPATAD* ADJOINT MODEL
!                Computation of the kappa function ("coefficient
!                of diffusion") relevant to the heat variables.
!                Used by SLHD scheme.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *GP_KAPPATAD(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KPROMA      - horizontal dimension.
!          KSTART      - first element of work.
!          KPROF       - depth of work.
!          KFLEV       - number of layers.
!          PDTS2       - 0.5*"pdt", where "pdt" =
!                        time step for the first time-integration step of
!                        a leap-frog scheme or all time-integration steps of
!                        a two-time level scheme; 2*time step for the following
!                        time-integration steps of a leap-frog scheme.
!          PTT0        - temperature at time t.
!          PTT0L       - zonal derivative of temperature at time t.
!          PTT0M       - meridional derivative of temperature at time t.
!          PSPT0       - "ln(prehyds)" at t.
!          PSPT0L      - zonal derivative of "ln(prehyds)" at t.
!          PSPT0M      - merid derivative of "ln(prehyds)" at t.
!          PTT05       - trajectory of temperature at time t.
!          PTT0L5      - trajectory of zonal derivative of temperature at time t.
!          PTT0M5      - trajectory of meridional derivative of temperature at time t.
!          PSPT05      - trajectory of "ln(prehyds)" at t.
!          PSPT0L5     - trajectory of zonal derivative of "ln(prehyds)" at t.
!          PSPT0M5     - trajectory of merid derivative of "ln(prehyds)" at t.
!          PSLHDA      - Scaling factor of the deformation in f(d) function
!                        (including the model resolution correction)
!          PSLHDD0     - Treshold for deformation tensor enhancement
!          PKAPPA      - Kappa heat function ("coefficient of SLHD")
!                        based on the rescaled horizontal grad of
!                        T-alpha(ln Pi_s) evaluated at instant "t".

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           none
!           Called by LATTE_KAPPAAD.

!     Reference.
!     ----------
!             IFS documentation about semi-Lagrangian scheme.

!  Author.
!  -------
!    Original F. VANA : AUGUST 2013.

!  Modifications.
!  --------------
!   21-May-2019 F. Vana     Computation restricted to levels where required.

! End Modifications
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM     ,JPRB
USE YOMHOOK      , ONLY : LHOOK    ,DR_HOOK

USE YOMDYN       , ONLY : TDYN
USE YOMCST       , ONLY : RKAPPA
USE INTDYN_MOD   , ONLY : YYTXYBT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT05(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT0L5(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT0M5(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT05(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0L5(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0M5(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTT0(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTT0L(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTT0M(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPT0(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPT0L(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPT0M(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDA(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDD0(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PKAPPA(KPROMA,KFLEV) 

!     -------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPRB) :: ZABSDTS2, ZTX, ZTY, ZT, ZRATA, ZRRATD
REAL(KIND=JPRB) :: ZTX5(KPROMA,KFLEV), ZTY5(KPROMA,KFLEV), ZT5(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZMAX(KPROMA,KFLEV), ZMAXE(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZKAPPA51(KPROMA,KFLEV), ZKAPPA52(KPROMA,KFLEV),&
 &  ZKAPPA53(KPROMA,KFLEV)

REAL(KIND=JPRB) :: ZPREF5(KPROMA,KFLEV),ZPREH5 (KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZEXN5  (KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZEXN51(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZPREF(KPROMA,KFLEV),ZPREH (KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZEXN  (KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZRTGR(KPROMA,KFLEV),ZRPRE(KPROMA,KFLEV),ZRPP(KPROMA,KFLEV)
REAL(KIND=JPRB), TARGET :: ZXYB5(KPROMA,KFLEV,YYTXYBT%NDIM)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     -------------------------------------------------------
#include "gphpre.intfb.h"
#include "gphpread.intfb.h"

!     -------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GP_KAPPATAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDSTA=>YDGEOMETRY%YRSTA)
ASSOCIATE(SLHDA0=>YDDYN%SLHDA0, SLHDA0T=>YDDYN%SLHDA0T, SLHDBT=>YDDYN%SLHDBT,   SLHDD00=>YDDYN%SLHDD00, &
& SLHDD00T=>YDDYN%SLHDD00T,   NLEV_SPONGE=>YDDYN%NLEV_SPONGE, STPREH=>YDSTA%STPREH)
!    --------------------------------------------------------

!  Positive time step (required in DFI)
ZABSDTS2=ABS(PDTS2)

! Scaling ratio for A and D specific to T
ZRATA=SLHDA0T/SLHDA0
ZRRATD=SLHDD00/SLHDD00T  ! Note this one is reversed

! Trajectory

!*       1.   Compute pressure related quantities

ZPREH5(KSTART:KPROF,KFLEV) = EXP(PSPT05(KSTART:KPROF))
CALL GPHPRE(KPROMA,KFLEV,KSTART,KPROF,YDVAB,ZPREH5,PXYB=ZXYB5,PRESF=ZPREF5)
ZEXN5(KSTART:KPROF,1:NLEV_SPONGE)=(STPREH(KFLEV)/ZPREF5(KSTART:KPROF,1:NLEV_SPONGE))**RKAPPA

!*       2.   Kappa function computation

DO JLEV=1,NLEV_SPONGE
  DO JROF=KSTART,KPROF
    ZT5(JROF,JLEV) =-RKAPPA*(YDVAB%VBH(JLEV-1)+YDVAB%VBH(JLEV))*0.5_JPRB&
       & *PTT05(JROF,JLEV)*ZPREH5(JROF,KFLEV)/ZPREF5(JROF,JLEV)
    ZEXN51(JROF,JLEV)=ZEXN5(JROF,JLEV)*ZT5(JROF,JLEV)  ! just to save some CPU
    ZTX5(JROF,JLEV)=ZEXN5(JROF,JLEV)*(PTT0L5(JROF,JLEV) + ZT5(JROF,JLEV) * PSPT0L5(JROF))
    ZTY5(JROF,JLEV)=ZEXN5(JROF,JLEV)*(PTT0M5(JROF,JLEV) + ZT5(JROF,JLEV) * PSPT0M5(JROF))
    !ZTX5(JROF,JLEV)= PTT0L5(JROF,JLEV) - RCORDIF(JLEV)*PSPT0L5(JROF)
    !ZTY5(JROF,JLEV)= PTT0M5(JROF,JLEV) - RCORDIF(JLEV)*PSPT0M5(JROF)
    ZKAPPA51(JROF,JLEV) = SQRT(ZTX5(JROF,JLEV)*ZTX5(JROF,JLEV) + ZTY5(JROF,JLEV)*ZTY5(JROF,JLEV))
    ZMAX(JROF,JLEV)=MAX(1.0_JPRB,ZKAPPA51(JROF,JLEV)*ZRRATD/PSLHDD0(JROF))
    ZMAXE(JROF,JLEV)=ZMAX(JROF,JLEV)**SLHDBT
    ZKAPPA52(JROF,JLEV)=2.0_JPRB*ZABSDTS2*ZRATA*PSLHDA(JROF)&
     & *ZKAPPA51(JROF,JLEV)*ZMAXE(JROF,JLEV)
    ZKAPPA53(JROF,JLEV)=ZKAPPA52(JROF,JLEV)/(ZKAPPA52(JROF,JLEV)+1.0_JPRB)
  ENDDO
ENDDO


! Adjoint
ZPREF(:,:) = 0.0_JPRB
ZPREH(:,:) = 0.0_JPRB
ZRTGR(:,:) = 0.0_JPRB
ZRPRE(:,:) = 0.0_JPRB
ZRPP(:,:)  = 0.0_JPRB

!*         Kappa function computation
DO JLEV=1,NLEV_SPONGE
!DEC$ IVDEP
  DO JROF=KSTART,KPROF
    PKAPPA(JROF,JLEV)= PKAPPA(JROF,JLEV)&
     & /(ZKAPPA52(JROF,JLEV)+1.0_JPRB)**2._JPRB

    PKAPPA(JROF,JLEV)=2.0_JPRB*ZABSDTS2*ZRATA*PSLHDA(JROF)* (&
     & PKAPPA(JROF,JLEV)*ZMAXE(JROF,JLEV) + ZKAPPA51(JROF,JLEV)*SLHDBT*(ZMAXE(JROF,JLEV)/ZMAX(JROF,JLEV))&
     & *(0.5_JPRB-SIGN(0.5_JPRB,1._JPRB-ZKAPPA51(JROF,JLEV)*ZRRATD/PSLHDD0(JROF)))&
     & * PKAPPA(JROF,JLEV)*ZRRATD/PSLHDD0(JROF)              )

    ZTX = ZTX5(JROF,JLEV)/ZKAPPA51(JROF,JLEV) * PKAPPA(JROF,JLEV)
    ZTY = ZTY5(JROF,JLEV)/ZKAPPA51(JROF,JLEV) * PKAPPA(JROF,JLEV)
    ZEXN(JROF,JLEV)=(PTT0L5(JROF,JLEV) + ZT5(JROF,JLEV) * PSPT0L5(JROF))*ZTX&
     &            + (PTT0M5(JROF,JLEV) + ZT5(JROF,JLEV) * PSPT0M5(JROF))*ZTY
    PTT0L(JROF,JLEV)=PTT0L(JROF,JLEV) + ZEXN5(JROF,JLEV)*ZTX
    PTT0M(JROF,JLEV)=PTT0M(JROF,JLEV) + ZEXN5(JROF,JLEV)*ZTY
    ZT  = ZEXN5(JROF,JLEV)*PSPT0L5(JROF)*ZTX + ZEXN5(JROF,JLEV)*PSPT0M5(JROF)*ZTY
    PSPT0L(JROF) = PSPT0L(JROF) + ZEXN51(JROF,JLEV)*ZTX
    PSPT0M(JROF) = PSPT0M(JROF) + ZEXN51(JROF,JLEV)*ZTY
    !PTT0L(JROF,JLEV) = PTT0L(JROF,JLEV) + ZTX
    !PTT0M(JROF,JLEV) = PTT0M(JROF,JLEV) + ZTY
    !PSPT0L(JROF)     = PSPT0L(JROF) - RCORDIF(JLEV)*ZTX
    !PSPT0M(JROF)     = PSPT0M(JROF) - RCORDIF(JLEV)*ZTY
    PTT0 (JROF,JLEV)  = PTT0 (JROF,JLEV)  + ZT5(JROF,JLEV)/PTT05(JROF,JLEV)*ZT
    ZPREH (JROF,KFLEV)= ZPREH(JROF,KFLEV) + ZT5(JROF,JLEV)/ZPREH5(JROF,KFLEV)*ZT
    ZPREF (JROF,JLEV) = ZPREF(JROF,JLEV)  - ZT5(JROF,JLEV)/ZPREF5(JROF,JLEV)*ZT
  ENDDO
ENDDO

!*         Compute pressure related quantities
ZPREF(KSTART:KPROF,1:NLEV_SPONGE) = ZPREF(KSTART:KPROF,1:NLEV_SPONGE)&
 & -RKAPPA*ZEXN5(KSTART:KPROF,1:NLEV_SPONGE)/ZPREF5(KSTART:KPROF,1:NLEV_SPONGE)*ZEXN(KSTART:KPROF,1:NLEV_SPONGE)
CALL GPHPREAD(KPROMA,KFLEV,KSTART,KPROF,YDVAB,ZPREH,ZPREH5,PXYB5=ZXYB5,PRESF=ZPREF)
PSPT0(KSTART:KPROF)=PSPT0(KSTART:KPROF)+EXP(PSPT05(KSTART:KPROF))*ZPREH (KSTART:KPROF,KFLEV)

!In any case the output for PKAPPA should be 0.

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GP_KAPPATAD',1,ZHOOK_HANDLE)
END SUBROUTINE GP_KAPPATAD

