SUBROUTINE GP_KAPPATL( &
 ! --- INPUT --------------------------------------------------
 & YDDYN,KPROMA,KSTART,KPROF,KFLEV,&
 & PDTS2, &
 & PUT0L5,PVT0L5,PVORT05,PDIVT05,&
 & PUT0L,PVT0L,PVORT0,PDIVT0,&
 & PSLHDA,PSLHDD0,&
 ! --- OUTPUT -------------------------------------------------
 & PKAPPA)

!**** *GP_KAPPATL*      tangent-linear version 
!                Computation of the kappa function ("coefficient
!                of diffusion") based on rescaled horizontal
!                deformation of the flow.
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
!   --------------------- trajectory variables -----------------------
!          PUT0L5      - zonal derivative of U-wind at time t.
!          PVT0L5      - zonal derivative of V-wind at time t.
!          PVORT05     - vorticity at time t.
!          PDIVT05     - divergence at t.
!   ------------------------ end of traj. variables -------------------
!          PUT0L       - zonal derivative of U-wind at time t.
!          PVT0L       - zonal derivative of V-wind at time t.
!          PVORT0      - vorticity at time t.
!          PDIVT0      - divergence at t.
!          PSLHDA      - Scaling factor of the deformation in f(d) function
!                        (including the model resolution correction)
!          PSLHDD0     - Treshold for deformation tensor enhancement

!        OUTPUT:
!          PKAPPA     - Kappa function ("coefficient of SLHD")
!                       based on the rescaled horizontal deformation
!                       of the flow evaluated at instant "t".

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           none
!           Called by LATTE_KAPPATL.

!     Reference.
!     ----------
!             Arpege documentation about semi-Lagrangian scheme.

!  Author.
!  -------
!    Original F. VANA : DECEMBER 2008.

!  Modifications.
!  --------------
!   21-May-2019 F. Vana     Computation restricted to levels where required.

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK

USE YOMDYNA  , ONLY : LSLHD_STATIC
USE YOMDYN   , ONLY : TDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT0L5(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVT0L5(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVORT05(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIVT05(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT0L(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVT0L(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVORT0(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIVT0(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDA(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDD0(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKAPPA(KPROMA,KFLEV) 

!     -------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV,JROF
REAL(KIND=JPRB) :: ZABSDTS2
REAL(KIND=JPRB) :: ZKAP5, ZMAX, ZMAXE
REAL(KIND=JPRB) :: ZDT, ZDS, ZDT5, ZDS5
REAL(KIND=JPRB) :: ZKAPPA5(KPROMA,KFLEV) 
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     -------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GP_KAPPATL',0,ZHOOK_HANDLE)
ASSOCIATE(SLHDB=>YDDYN%SLHDB, SLHDDIV=>YDDYN%SLHDDIV, &
 & SLHDRATDDIV=>YDDYN%SLHDRATDDIV,NLEV_SPONGE=>YDDYN%NLEV_SPONGE)
!    --------------------------------------------------------

!*       0.   Basic settings

!  Positive time step (required in DFI)
ZABSDTS2=ABS(PDTS2)


!*       1.   Kappa function computation

IF (LSLHD_STATIC) THEN
  PKAPPA(KSTART:KPROF,1:KFLEV) =0.0_JPRB
ELSE

  DO JLEV=1,NLEV_SPONGE
!DEC$ IVDEP
    DO JROF=KSTART,KPROF

      ! diagnose wind derivatives 
      ! (here rather shearing and tensions terms of deformation)
      ZDT= 2._JPRB*PUT0L(JROF,JLEV) - PDIVT0(JROF,JLEV)
      ZDS= 2._JPRB*PVT0L(JROF,JLEV) - PVORT0(JROF,JLEV)

      ZDT5= 2._JPRB*PUT0L5(JROF,JLEV) - PDIVT05(JROF,JLEV)
      ZDS5= 2._JPRB*PVT0L5(JROF,JLEV) - PVORT05(JROF,JLEV)

      ! compute kappa function
      ZKAP5 = SQRT(ZDT5*ZDT5+ZDS5*ZDS5)
      PKAPPA(JROF,JLEV)=(1.0_JPRB-SLHDDIV)*(ZDT5*ZDT + ZDS5*ZDS)/ZKAP5&
       & + SLHDDIV*SLHDRATDDIV*SIGN(1.0_JPRB,PDIVT05(JROF,JLEV))&
       & * PDIVT0(JROF,JLEV)
      ZKAPPA5(JROF,JLEV)=  (1.0_JPRB-SLHDDIV)* ZKAP5&
       & + SLHDDIV*SLHDRATDDIV*ABS(PDIVT05(JROF,JLEV))

      ZMAX=MAX(1.0_JPRB,ZKAPPA5(JROF,JLEV)/PSLHDD0(JROF))
      ZMAXE=ZMAX**SLHDB
      PKAPPA(JROF,JLEV)=2.0_JPRB*ZABSDTS2*PSLHDA(JROF)* (&
       & PKAPPA(JROF,JLEV)*ZMAXE + ZKAPPA5(JROF,JLEV)*SLHDB*(ZMAXE/ZMAX)&
       & *(0.5_JPRB-SIGN(0.5_JPRB,1._JPRB-ZKAPPA5(JROF,JLEV)/PSLHDD0(JROF)))&
       & * PKAPPA(JROF,JLEV)/PSLHDD0(JROF)              )
      ZKAPPA5(JROF,JLEV)=2.0_JPRB*ZABSDTS2*PSLHDA(JROF)&
       & *ZKAPPA5(JROF,JLEV)*ZMAXE

      PKAPPA(JROF,JLEV)= PKAPPA(JROF,JLEV)&
       & /(ZKAPPA5(JROF,JLEV)+1.0_JPRB)**2._JPRB
      !ZKAPPA5(JROF,JLEV)=ZKAPPA5(JROF,JLEV)/(ZKAPPA5(JROF,JLEV)+1.0_JPRB)

    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GP_KAPPATL',1,ZHOOK_HANDLE)
END SUBROUTINE GP_KAPPATL
