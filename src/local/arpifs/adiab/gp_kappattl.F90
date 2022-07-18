SUBROUTINE GP_KAPPATTL(&
 ! --- INPUT --------------------------------------------------
 & YDGEOMETRY, YDDYN,KPROMA,KSTART,KPROF,KFLEV,&
 & PDTS2,PTT05,PTT0L5,PTT0M5,PSPT05,PSPT0L5,PSPT0M5,&
 & PTT0,PTT0L,PTT0M,PSPT0,PSPT0L,PSPT0M,PSLHDA,PSLHDD0,&
 ! --- OUTPUT -------------------------------------------------
 & PKAPPA)

!**** *GP_KAPPATTL* input for SLHD diffusion.
!                Computation of the kappa function ("coefficient
!                of diffusion") relevant to the heat variables.
!                Used by SLHD scheme.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *GP_KAPPATL(..)

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
!           Called by LATTE_KAPPAT.

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

USE YOMDYNA      , ONLY : YRDYNA
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
REAL(KIND=JPRB) :: ZABSDTS2, ZTX, ZTY, ZTX5, ZTY5, ZT, ZT5
REAL(KIND=JPRB) :: ZMAX, ZMAXE, ZRATA, ZRRATD
REAL(KIND=JPRB) :: ZKAPPA5(KPROMA,KFLEV)

REAL(KIND=JPRB) :: ZEXN (KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZPREF(KPROMA,KFLEV),ZPREH (KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZEXN5(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZPREF5(KPROMA,KFLEV),ZPREH5 (KPROMA,0:KFLEV)
REAL(KIND=JPRB), TARGET :: ZXYB5(KPROMA,KFLEV,YYTXYBT%NDIM)



REAL(KIND=JPRB) :: ZHOOK_HANDLE
!     -------------------------------------------------------

#include "gphpre.intfb.h"
#include "gphpretl.intfb.h"

!     -------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GP_KAPPATTL',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDSTA=>YDGEOMETRY%YRSTA)
ASSOCIATE(SLHDA0=>YDDYN%SLHDA0, SLHDA0T=>YDDYN%SLHDA0T, SLHDBT=>YDDYN%SLHDBT,   SLHDD00=>YDDYN%SLHDD00, &
& SLHDD00T=>YDDYN%SLHDD00T,   NLEV_SPONGE=>YDDYN%NLEV_SPONGE, STPREH=>YDSTA%STPREH)
!    --------------------------------------------------------

!*       0.   Basic settings

!  Positive time step (required in DFI)
ZABSDTS2=ABS(PDTS2)

! Scaling ratio for A and D specific to T
ZRATA=SLHDA0T/SLHDA0
ZRRATD=SLHDD00/SLHDD00T  ! Note this one is reversed

!*       1.   Compute pressure related quantities
ZPREH (KSTART:KPROF,KFLEV) = EXP(PSPT05(KSTART:KPROF))*PSPT0(KSTART:KPROF)
ZPREH5(KSTART:KPROF,KFLEV) = EXP(PSPT05(KSTART:KPROF))

CALL GPHPRE(KPROMA,KFLEV,KSTART,KPROF,YDVAB,ZPREH5,PXYB=ZXYB5,PRESF=ZPREF5)
CALL GPHPRETL(KPROMA,KFLEV,KSTART,KPROF,YDVAB,ZPREH,ZPREH5,PXYB5=ZXYB5,PRESF=ZPREF)

ZEXN5(KSTART:KPROF,1:NLEV_SPONGE)=(STPREH(KFLEV)/ZPREF5(KSTART:KPROF,1:NLEV_SPONGE))**RKAPPA
ZEXN(KSTART:KPROF,1:NLEV_SPONGE)=-RKAPPA*ZEXN5(KSTART:KPROF,1:NLEV_SPONGE)/ZPREF5(KSTART:KPROF,1:NLEV_SPONGE)&
 &  *  ZPREF(KSTART:KPROF,1:NLEV_SPONGE)


!*       2.   Kappa function computation
IF (YRDYNA%LSLHD_STATIC) THEN
  PKAPPA(KSTART:KPROF,1:KFLEV)=0.0_JPRB
ELSE

  DO JLEV=1,NLEV_SPONGE
    DO JROF=KSTART,KPROF

      ! compute kappa function
      ZT5 =-RKAPPA*(YDVAB%VBH(JLEV-1)+YDVAB%VBH(JLEV))*0.5_JPRB&
       & *PTT05(JROF,JLEV)*ZPREH5(JROF,KFLEV)/ZPREF5(JROF,JLEV)
      ZT  = ZT5/PTT05(JROF,JLEV)   *PTT0 (JROF,JLEV)&
       & +  ZT5/ZPREH5(JROF,KFLEV) *ZPREH (JROF,KFLEV)&
       & -  ZT5/ZPREF5(JROF,JLEV)  *ZPREF (JROF,JLEV) 

      ZTX5=ZEXN5(JROF,JLEV)*(PTT0L5(JROF,JLEV) + ZT5 * PSPT0L5(JROF))
      ZTX = ZEXN(JROF,JLEV)*(PTT0L5(JROF,JLEV) + ZT5 * PSPT0L5(JROF))&
       & + ZEXN5(JROF,JLEV)*(PTT0L(JROF,JLEV)  + ZT  * PSPT0L5(JROF)&
       &                                       + ZT5 * PSPT0L(JROF) )

      ZTY5=ZEXN5(JROF,JLEV)*(PTT0M5(JROF,JLEV) + ZT5 * PSPT0M5(JROF))
      ZTY = ZEXN(JROF,JLEV)*(PTT0M5(JROF,JLEV) + ZT5 * PSPT0M5(JROF))&
       & + ZEXN5(JROF,JLEV)*(PTT0M(JROF,JLEV)  + ZT  * PSPT0M5(JROF)&
       &                                       + ZT5 * PSPT0M(JROF) )


      !ZTX5= PTT0L5(JROF,JLEV) - RCORDIF(JLEV)*PSPT0L5(JROF)
      !ZTY5= PTT0M5(JROF,JLEV) - RCORDIF(JLEV)*PSPT0M5(JROF)
      ZKAPPA5(JROF,JLEV)=  SQRT(ZTX5*ZTX5 + ZTY5*ZTY5)
      !ZTX = PTT0L(JROF,JLEV)  - RCORDIF(JLEV)*PSPT0L(JROF)
      !ZTY = PTT0M(JROF,JLEV)  - RCORDIF(JLEV)*PSPT0M(JROF)
      PKAPPA(JROF,JLEV) =   (ZTX5*ZTX + ZTY5*ZTY)/ZKAPPA5(JROF,JLEV)

      ! rescale and normalize it
      ZMAX=MAX(1.0_JPRB,ZKAPPA5(JROF,JLEV)*ZRRATD/PSLHDD0(JROF))
      ZMAXE=ZMAX**SLHDBT
      PKAPPA(JROF,JLEV)=2.0_JPRB*ZABSDTS2*ZRATA*PSLHDA(JROF)* (&
       & PKAPPA(JROF,JLEV)*ZMAXE + ZKAPPA5(JROF,JLEV)*SLHDBT*(ZMAXE/ZMAX)&
       & *(0.5_JPRB-SIGN(0.5_JPRB,1._JPRB-ZKAPPA5(JROF,JLEV)*ZRRATD/PSLHDD0(JROF)))&
       & * PKAPPA(JROF,JLEV)*ZRRATD/PSLHDD0(JROF)              )
      ZKAPPA5(JROF,JLEV)=2.0_JPRB*ZABSDTS2*ZRATA*PSLHDA(JROF)&
       & *ZKAPPA5(JROF,JLEV)*ZMAXE

      PKAPPA(JROF,JLEV)= PKAPPA(JROF,JLEV)&
       & /(ZKAPPA5(JROF,JLEV)+1.0_JPRB)**2._JPRB
      ZKAPPA5(JROF,JLEV)=  ZKAPPA5(JROF,JLEV)/(ZKAPPA5(JROF,JLEV)+1.0_JPRB)

    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GP_KAPPATTL',1,ZHOOK_HANDLE)
END SUBROUTINE GP_KAPPATTL
