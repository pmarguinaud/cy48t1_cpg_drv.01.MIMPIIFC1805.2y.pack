SUBROUTINE GP_KAPPAT(&
 ! --- INPUT --------------------------------------------------
 & YDCST, YDGEOMETRY, YDDYN,KPROMA,KSTART,KPROF,KFLEV,&
 & PDTS2,PTT0,PTT0L,PTT0M,PSPT0,PSPT0L,PSPT0M,PSLHDA,PSLHDD0,&
 ! --- OUTPUT -------------------------------------------------
 & PKAPPA)

!**** *GP_KAPPAT* input for SLHD diffusion.
!                Computation of the kappa function ("coefficient
!                of diffusion") relevant to the heat variables.
!                Used by SLHD scheme.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *GP_KAPPAT(..)

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
!          PSLHDA      - Scaling factor of the deformation in f(d) function
!                        (including the model resolution correction)
!          PSLHDD0     - Treshold for deformation tensor enhancement

!        OUTPUT:
!          PKAPPA     - Kappa heat function ("coefficient of SLHD")
!                       based on the rescaled horizontal grad of
!                       T-alpha(ln Pi_s) evaluated at instant "t".

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           none
!           Called by LATTE_KAPPA.

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

USE YOMDYNA      , ONLY : LSLHD_STATIC
USE YOMDYN       , ONLY : TDYN
USE YOMCST       , ONLY : TCST
USE INTDYN_MOD   , ONLY : YYTXYB


!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT0(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT0L(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT0M(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0L(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0M(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDA(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDD0(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKAPPA(KPROMA,KFLEV) 

!     -------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPRB) :: ZABSDTS2, ZTX, ZTY,&
 &  ZRATA, ZRRATD, ZC, ZAUX, ZDT
REAL(KIND=JPRB) :: ZEXN (KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZPREF(KPROMA,KFLEV),ZPREH (KPROMA,0:KFLEV)
REAL(KIND=JPRB), TARGET :: ZXYB(KPROMA,KFLEV,YYTXYB%NDIM)
!REAL(KIND=JPRB) :: ZDELP(KPROMA,KFLEV),ZRDELP(KPROMA,KFLEV)
!REAL(KIND=JPRB) :: ZLNPR(KPROMA,KFLEV),ZALPH (KPROMA,KFLEV)
!REAL(KIND=JPRB) :: ZDUM (KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZTHETAH (KPROMA,0:KFLEV)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     -------------------------------------------------------
#include "gphpre.intfb.h"

!     -------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GP_KAPPAT',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDSTA=>YDGEOMETRY%YRSTA)
ASSOCIATE(SLHDA0=>YDDYN%SLHDA0, SLHDA0T=>YDDYN%SLHDA0T, SLHDBT=>YDDYN%SLHDBT,   SLHDD00=>YDDYN%SLHDD00,  &
& SLHDD00T=>YDDYN%SLHDD00T, SLHDHOR=>YDDYN%SLHDHOR,   NLEV_SPONGE=>YDDYN%NLEV_SPONGE,STPREH=>YDSTA%STPREH&
& )
!    --------------------------------------------------------

!*       0.   Basic settings

!  Positive time step (required in DFI)
ZABSDTS2=ABS(PDTS2)

! Scaling ratio for A and D specific to T
ZRATA=SLHDA0T/SLHDA0
ZRRATD=SLHDD00/SLHDD00T  ! Note this one is reversed


!*       1.   Compute pressure related quantities
ZPREH (KSTART:KPROF,KFLEV) = EXP(PSPT0(KSTART:KPROF))

CALL GPHPRE(KPROMA,KFLEV,KSTART,KPROF,YDVAB,ZPREH,PXYB=ZXYB,PRESF=ZPREF)

ZEXN(KSTART:KPROF,1:NLEV_SPONGE)=(STPREH(KFLEV)/ZPREF(KSTART:KPROF,1:NLEV_SPONGE))**YDCST%RKAPPA

!*     Potential temperature on half levels
DO JLEV=1,KFLEV-1
  DO JROF=KSTART,KPROF
    ZTHETAH(JROF,JLEV)=0.5_JPRB*(PTT0(JROF,JLEV)+PTT0(JROF,JLEV+1))*&
      &  (STPREH(KFLEV)/ZPREH(JROF,JLEV))**YDCST%RKAPPA
  ENDDO
ENDDO
 !  Boundaries, extrapolated only
DO JROF=KSTART,KPROF
  ZTHETAH(JROF,0)=(1.5_JPRB*PTT0(JROF,1)-0.5_JPRB*PTT0(JROF,2))*&
   & (STPREH(KFLEV)/MAX(1._JPRB,ZPREH(JROF,0)))**YDCST%RKAPPA
  ZTHETAH(JROF,KFLEV)=(1.5_JPRB*PTT0(JROF,KFLEV)-0.5_JPRB*PTT0(JROF,KFLEV-1))*&
   & (STPREH(KFLEV)/ZPREH(JROF,KFLEV))**YDCST%RKAPPA
ENDDO



!*       2.   Kappa function computation

IF (LSLHD_STATIC) THEN
  ! Trivial solution
  PKAPPA(KSTART:KPROF,1:KFLEV)=1.0_JPRB
ELSE

  DO JLEV=1,NLEV_SPONGE
    ZC=0.5_JPRB*SLHDHOR*(YDVAB%VBH(JLEV)+YDVAB%VBH(JLEV-1))
    DO JROF=KSTART,KPROF

      ! fields used for theta derivatives along pressure surfaces
      ZAUX=ZC*ZXYB(JROF,JLEV,YYTXYB%M_RDELP)*ZPREH(JROF,KFLEV)
      ZDT =ZTHETAH(JROF,JLEV)-ZTHETAH(JROF,JLEV-1)

      ! compute kappa function
      ZTX=PTT0L(JROF,JLEV)-YDCST%RKAPPA*(YDVAB%VBH(JLEV-1)+YDVAB%VBH(JLEV))*0.5_JPRB&
       & *PTT0(JROF,JLEV)*ZPREH(JROF,KFLEV)/ZPREF(JROF,JLEV)*PSPT0L(JROF)
      ZTX=ZTX*ZEXN(JROF,JLEV) - ZAUX*PSPT0L(JROF)*ZDT
      ZTY=PTT0M(JROF,JLEV)-YDCST%RKAPPA*(YDVAB%VBH(JLEV-1)+YDVAB%VBH(JLEV))*0.5_JPRB&
       & *PTT0(JROF,JLEV)*ZPREH(JROF,KFLEV)/ZPREF(JROF,JLEV)*PSPT0M(JROF)
      ZTY=ZTY*ZEXN(JROF,JLEV) - ZAUX*PSPT0M(JROF)*ZDT
      !ZTX= PTT0L(JROF,JLEV) - RCORDIF(JLEV)*PSPT0L(JROF)
      !ZTY= PTT0M(JROF,JLEV) - RCORDIF(JLEV)*PSPT0M(JROF)
      PKAPPA(JROF,JLEV)=  SQRT(ZTX*ZTX + ZTY*ZTY)

      ! rescale and normalize it
      PKAPPA(JROF,JLEV)=2.0_JPRB*ZABSDTS2*&
        & ZRATA*PSLHDA(JROF)*PKAPPA(JROF,JLEV)*&
        & (MAX(1.0_JPRB,PKAPPA(JROF,JLEV)*ZRRATD/PSLHDD0(JROF)))**SLHDBT
      PKAPPA(JROF,JLEV)= PKAPPA(JROF,JLEV)/(PKAPPA(JROF,JLEV)+1.0_JPRB)

    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GP_KAPPAT',1,ZHOOK_HANDLE)
END SUBROUTINE GP_KAPPAT
