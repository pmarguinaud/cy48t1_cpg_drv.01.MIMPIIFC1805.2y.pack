SUBROUTINE LATTESTL(YDGEOMETRY,YDGMV,YDGMV5,YDRIP,YDML_DYN,KST,KPROF,PDTS2,PBDT,PESGP,PESGM,&
 & POROGL,POROGM,PSDIV0,PSDVBC,PRES0,PGMVS,&
 & PGMV,PGMVT1S,PB1,PB2,&
 & PSDVBC5,PRES5,PGMV5,PGMV5S)

!**** *LATTESTL* Semi-Lagrangian scheme.  (tangent-linear version)
!                Computation of the t and t-dt useful quantities
!                 at grid-points. Equations for bi-dimensional
!                 variables (continuity equation).

!     Purpose.
!     --------
!        * This subroutine computes the equation quantities to be
!          interpolated at each grid-point of the colocation grid.
!          The terms considered here are the explicit terms and
!          the explicit part of the semi-implicit scheme (terms
!          previously computed in LASSIE or LANHSI).
!          Equations considered here are equations for bi-dimensional
!          variables: continuity equation.
!        * Remark 1: when an alternate averaging is used for linear terms
!          in the 2TL SL scheme, the first timestep is treated differently
!          (first order uncentering), no uncentering is applied to the
!          total term ("cursive A") and the term saved in P[X]NLT9
!          is [ (Delta t/2) ("cursive A" - (1 + xidt) beta "cursive B") ]
!          instead of [ (Delta t/2) ("cursive A" - beta "cursive B") ].
!        * Remark 2: for lsettls=true, uncentering is applied to
!          the 'stable' extrapolation if vesl > 0 to avoid instability
!          in the momentum equation.
!        * Remark 3: for PC schemes:
!          - this routine is called for nsiter=0.
!          - this routine is called for the predictor of lpc_full.
!          - this routine is called for the corrector of lpc_full.

!**   Interface.
!     ----------
!        *CALL* *LATTESTL(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST         - first element of work.
!          KPROF       - depth of work.
!          PDTS2       - 0.5*"pdt", where "pdt" =
!                        time step for the first time-integration step of
!                        a leap-frog scheme or all time-integration steps of
!                        a two-time level scheme; 2*time step for the following
!                        time-integration steps of a leap-frog scheme.
!          PBDT        - zbt if semi-implicit scheme with unreduced
!                        divergence, zbt*(c**2/GM**2) if semi-implicit
!                        scheme with reduced divergence, where zbt
!                        is equal to PDTS2*BETADT (BETADT is in YOMDYN) or zero
!                        according to the value of configuration.
!          PESGP       - (1 + uncentering factor).
!          PESGM       - (1 - uncentering factor).
!          POROGL      - Zonal derivative of orography.
!          POROGM      - Meridian derivative of orography.
!          PSDIV0      - SI term at time t for continuity equation (Nu*D).
!          PSDVBC      - vertical integral of divergence computed in "gpcty",
!                        including the "lrubc" and "delta m=1" contributions
!                        of (etadot d prehyd/d eta).
!          PRES0       - hydrostatic pressure at half levels at t.
!          PGMVS       - GMVS variables at t-dt and t.

!        INPUT/OUTPUT:
!          PGMV        - GMV variables at t-dt and t.
!          PGMVT1S     - GMVS variables at t+dt.
!          PB1         - "SLBUF1" buffer for interpolations.
!          PB2         - "SLBUF2" buffer.

!        INPUT-TRAJECTORY:
!          PSDVBC5     - trajectory for PSDVBC.
!          PRES5       - trajectory for PRES0.
!          PGMV5,PGMV5S: trajectory for PGMV,PGMVS.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           none
!           Called by LACDYNTL.

!     Reference.
!     ----------
!             Arpege documentation about semi-Lagrangian scheme.

!     Author.
!     -------
!      C. Temperton (ECMWF)
!      Original : 99/08/09.

!     Modifications.
!     --------------
!      K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!      K. Yessad May 2009: rewrite to match the current direct code.
!      K.Yessad (Aug 2009): remove NTRSLTYPE/=2 cases
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (Oct 2013): allow NESC without ICI-PC scheme.
!      K. Yessad (July 2014): Move some variables, rename some variables.
!      K. Yessad (Dec 2016): Prune obsolete options.
!      K. Yessad (Feb 2018): remove deep-layer formulations.
! End Modifications
!     ------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE YOMGMV             , ONLY : TGMV
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE YOMCST             , ONLY : RD
USE YOMSTA             , ONLY : RTSUR
USE YOMCT0             , ONLY : LTWOTL
USE YOMDYNA            , ONLY : YRDYNA
USE YOMCT3             , ONLY : NSTEP
USE YOMCVER            , ONLY : LVERTFE
USE YOMRIP             , ONLY : TRIP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)           ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)               ,INTENT(INOUT) :: YDGMV
TYPE(TGMV)               ,INTENT(INOUT) :: YDGMV5
TYPE(TRIP)               ,INTENT(IN)    :: YDRIP
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN)    :: YDML_DYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBDT(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESGP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESGM 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSDIV0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSDVBC(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRES0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1S(YDGEOMETRY%YRDIM%NPROMA,YDGMV%YT1%NDIMS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSDVBC5(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRES5(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMV5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV5%YT5%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMV5S(YDGEOMETRY%YRDIM%NPROMA,YDGMV5%YT5%NDIMS)
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZMOY1SP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSPNLT0(YDGEOMETRY%YRDIM%NPROMA)
!nyc-pc REAL(KIND=JPRB) :: ZSPNLT1(NPROMA)
REAL(KIND=JPRB) :: ZSPNLT_FE(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSPTB(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZSPT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZRPRES5

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZCMSLP, ZXIGP
REAL(KIND=JPRB) :: ZXIDT0, ZXIDT9

LOGICAL :: LLCT, LLCTC
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "verint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LATTESTL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA,   &
& YDVFE=>YDGEOMETRY%YRVFE,  YDDYN=>YDML_DYN%YRDYN,YDPTRSLB1=>YDML_DYN%YRPTRSLB1,    YDPTRSLB2=>YDML_DYN%YRPTRSLB2 &
& )

ASSOCIATE(NPROMA=>YDDIM%NPROMA,   NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA,   NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
& NVLAG=>YDDYN%NVLAG,   RCMSLP0=>YDDYN%RCMSLP0, XIDT=>YDDYN%XIDT,   NFOST=>YDRIP%NFOST,   YT0=>YDGMV%YT0,           &
& YT1=>YDGMV%YT1, YT9=>YDGMV%YT9,   YT5=>YDGMV5%YT5,   MSLB1C9=>YDPTRSLB1%MSLB1C9,   MSLB1SP9=>YDPTRSLB1%MSLB1SP9,  &
& MSLB2SPSI=>YDPTRSLB2%MSLB2SPSI)
!     ------------------------------------------------------------------

!############################################
! 1. AUXILIARITIES
!############################################

LLCT = YRDYNA%LPC_FULL .AND. NCURRENT_ITER > 0 ! corrector for LPC_FULL
LLCTC = YRDYNA%LPC_CHEAP .AND. NCURRENT_ITER > 0

ZXIDT0=1.0_JPRB+XIDT
ZXIDT9=1.0_JPRB+XIDT
ZXIGP=1.0_JPRB+XIDT

ZCMSLP=RCMSLP0/(RD*RTSUR)

!     ------------------------------------------------------------------

!############################################
! 2. NONLINEAR MODEL
!############################################

DO JROF=KST,KPROF
  ZRPRES5=1.0_JPRB/PRES5(JROF,NFLEVG)
  ZSPTB(JROF)=ZRPRES5*PSDVBC(JROF,NFLEVG)&
   & - ZRPRES5*ZRPRES5*PSDVBC5(JROF,NFLEVG)*PRES0(JROF,NFLEVG)
ENDDO

! linear model at time t

DO JROF=KST,KPROF
  PB2(JROF,MSLB2SPSI)=PBDT(JROF)*PSDIV0(JROF)
  ZSPTB(JROF)=-PDTS2*ZSPTB(JROF)
ENDDO

! nonlinear model at time t

DO JLEV=1,NFLEVG
  DO JROF=KST,KPROF
    ZMOY1SP(JROF,JLEV)=PDTS2*(&
     & PGMV5(JROF,JLEV,YT5%MU)*PGMVS(JROF,YT0%MSPL)&
     & +PGMV5S(JROF,YT5%MSPL)*PGMV(JROF,JLEV,YT0%MU)&
     & +ZCMSLP*POROGL(JROF)*PGMV(JROF,JLEV,YT0%MU)&
     & +PGMV5(JROF,JLEV,YT5%MV)*PGMVS(JROF,YT0%MSPM)&
     & +PGMV5S(JROF,YT5%MSPM)*PGMV(JROF,JLEV,YT0%MV)&
     & +ZCMSLP*POROGM(JROF)*PGMV(JROF,JLEV,YT0%MV) )
  ENDDO
  DO JROF=KST,KPROF
    ZMOY1SP(JROF,JLEV)=ZMOY1SP(JROF,JLEV)+ZSPTB(JROF)
  ENDDO
ENDDO

!*       *   Continuity equation.

IF (LTWOTL) THEN

  IF(.NOT.LLCT)THEN

    !############################################
    ! Predictor for LPC_FULL
    ! or case nsiter=0.
    !############################################

    IF (YRDYNA%LPC_CHEAP.AND.(.NOT.YRDYNA%LNESC)) THEN
      ! LPC_CHEAP=T is currently coded for LNESC only
      !  (for LNESC=F interpolated quantities are not the same ones
      !  at the predictor and corrector steps, and the LPC_CHEAP code
      !  currently assumes that they are identical).
      CALL ABOR1('LATTESTL: If LPC_CHEAP=T, LNESC should be T')
    ENDIF

    !---------------------------------------------------------------
    ! 3D buffers
    !---------------------------------------------------------------
    DO JLEV=1,NFLEVG

      ! nonlinear residual time t
      DO JROF=KST,KPROF
        ZSPNLT0(JROF)=ZMOY1SP(JROF,JLEV)+ZXIDT0*PB2(JROF,MSLB2SPSI)
      ENDDO

      ! Fill PB1(.,MSLB1C9),PGMVT1S(.,YT1%MSP),ZSPNLT_FE:
      IF (NVLAG == 2 .OR. NVLAG == 3) THEN

        IF ( NSTEP <= NFOST .OR. YRDYNA%LNESC ) THEN

          DO JROF=KST,KPROF
            PB1(JROF,MSLB1C9+JLEV-NFLSA)=PB1(JROF,MSLB1C9+JLEV-NFLSA)&
             & +PESGM*ZMOY1SP(JROF,JLEV)
          ENDDO
          IF (LVERTFE) THEN
            DO JROF=KST,KPROF
              ZSPNLT_FE(JROF,JLEV)=ZSPNLT0(JROF)&
               & *YDVAB%VDELB(JLEV)*YDVETA%VFE_RDETAH(JLEV)*PESGP
            ENDDO
          ELSE
            DO JROF=KST,KPROF
              PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)&
               & +PESGP*YDVAB%VDELB(JLEV)*ZSPNLT0(JROF)
            ENDDO
          ENDIF

        ELSEIF (YRDYNA%LSETTLS) THEN

          DO JROF=KST,KPROF
            PB1(JROF,MSLB1C9+JLEV-NFLSA)=PB1(JROF,MSLB1C9+JLEV-NFLSA)&
             & +PESGM*ZMOY1SP(JROF,JLEV)&
             & +(ZSPNLT0(JROF)-PGMV(JROF,JLEV,YT9%MSPNL))
          ENDDO
          IF( LVERTFE )THEN
            DO JROF=KST,KPROF
              ZSPNLT_FE(JROF,JLEV)=&
               & ZSPNLT0(JROF)*PESGP*YDVAB%VDELB(JLEV)*YDVETA%VFE_RDETAH(JLEV)
            ENDDO
          ELSE
            DO JROF=KST,KPROF
              PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)&
               & +PESGP*YDVAB%VDELB(JLEV)*ZSPNLT0(JROF)
            ENDDO
          ENDIF

        ELSE

          ! remaining case: ldsettls=false, lnesc=false.
          DO JROF=KST,KPROF
            PB1(JROF,MSLB1C9+JLEV-NFLSA)=PB1(JROF,MSLB1C9+JLEV-NFLSA)&
             & +PESGM*ZMOY1SP(JROF,JLEV)&
             & +0.5_JPRB*PESGM*(ZSPNLT0(JROF)-PGMV(JROF,JLEV,YT9%MSPNL))
          ENDDO
          IF (LVERTFE) THEN
            DO JROF=KST,KPROF
              ZSPNLT_FE(JROF,JLEV)=&
               & (1.5_JPRB*ZSPNLT0(JROF)-0.5_JPRB*PGMV(JROF,JLEV,YT9%MSPNL))&
               & *YDVAB%VDELB(JLEV)*YDVETA%VFE_RDETAH(JLEV)*PESGP
            ENDDO
          ELSE
            DO JROF=KST,KPROF
              PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)&
               & +PESGP*YDVAB%VDELB(JLEV)&
               & *(1.5_JPRB*ZSPNLT0(JROF)-0.5_JPRB*PGMV(JROF,JLEV,YT9%MSPNL))
            ENDDO
          ENDIF

        ENDIF

      ENDIF

      ! Save quantities for corrector step
      !nyc-pc  IF( LPC_FULL )THEN
      !nyc-pc    ! save nonlinear model at time t
      !nyc-pc    PGMV(KST:KPROF,JLEV,YT9%MCSPNL)=ZMOY1SP(KST:KPROF,JLEV)
      !nyc-pc  ENDIF

      IF( .NOT.YRDYNA%LNESC )THEN
        ! save of nonlinear residual at time t
        ! to be used as nonlinear residual at time t-dt next time step
        DO JROF=KST,KPROF
          PGMV(JROF,JLEV,YT9%MSPNL)=ZMOY1SP(JROF,JLEV)+ZXIDT9*PB2(JROF,MSLB2SPSI)
        ENDDO
      ENDIF

    ENDDO

    !------------------------------------------------------------------------
    ! 3D additional actions
    !------------------------------------------------------------------------
    IF (LVERTFE) THEN
      IF(NVLAG == 2 .OR. NVLAG == 3)THEN
        CALL VERINT(NPROMA,KST,KPROF,NFLEVG,NFLEVG+1,YDVFE%RINTE,ZSPNLT_FE,ZSPT1,0)
        DO JROF=KST,KPROF
          PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)+ZSPT1(JROF,NFLEVG+1)
        ENDDO
      ENDIF
    ENDIF

    !------------------------------------------------------------------------
    ! 2D buffers
    !------------------------------------------------------------------------
    IF(NVLAG == 2 .OR. NVLAG == 3)THEN
      DO JROF=KST,KPROF
        PB1(JROF,MSLB1SP9)=PB1(JROF,MSLB1SP9)+PGMVS(JROF,YT0%MSP)
      ENDDO
    ENDIF

  ELSE

    !############################################
    ! Corrector for LPC_FULL
    !############################################

    ! ky: piece of code not yet validated.
    !     Commented for the time being.

    !---------------------------------------------------------------
    ! 3D buffers
    !---------------------------------------------------------------
    !nyc-pc DO JLEV=1,NFLEVG

      !nyc-pc ! nonlinear residual time t

      !nyc-pc DO JROF=KST,KPROF
      !nyc-pc   ZSPNLT1(JROF)=ZMOY1SP(JROF,JLEV)+ZXIDT9*PB2(JROF,MSLB2SPSI)
      !nyc-pc ENDDO

      !nyc-pc IF (.NOT.LLCTC) THEN
      !nyc-pc   DO JROF=KST,KPROF
      !nyc-pc     PB1(JROF,MSLB1C9+JLEV-NFLSA)=PB1(JROF,MSLB1C9+JLEV-NFLSA) &
      !nyc-pc      & +PESGM*PGMV(JROF,JLEV,YT9%MCSPNL)
      !nyc-pc   ENDDO
      !nyc-pc ENDIF
      !nyc-pc IF (LVERTFE) THEN
      !nyc-pc   DO JROF=KST,KPROF
      !nyc-pc     ZSPNLT_FE(JROF,JLEV)=ZSPNLT1(JROF) &
      !nyc-pc      & *YDVAB%VDELB(JLEV)*YDVETA%VFE_RDETAH(JLEV)*PESGP
      !nyc-pc   ENDDO
      !nyc-pc ELSE
      !nyc-pc   DO JROF=KST,KPROF
      !nyc-pc     PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP) &
      !nyc-pc      & +YDVAB%VDELB(JLEV)*ZSPNLT1(JROF)*PESGP
      !nyc-pc   ENDDO
      !nyc-pc ENDIF

    !nyc-pc ENDDO

    !------------------------------------------------------------------------
    ! 3D additional actions
    !------------------------------------------------------------------------
    !nyc-pc IF (LVERTFE) THEN
    !nyc-pc   CALL VERINT(NPROMA,KST,KPROF,NFLEVG,NFLEVG+1,YDVFE%RINTE,ZSPNLT_FE,ZSPT1,0)
    !nyc-pc   DO JROF=KST,KPROF
    !nyc-pc     PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)+ZSPT1(JROF,NFLEVG+1)
    !nyc-pc   ENDDO
    !nyc-pc ENDIF

    !------------------------------------------------------------------------
    ! 2D buffers
    !------------------------------------------------------------------------
    !nyc-pc IF (.NOT.LLCTC) THEN
    !nyc-pc   DO JROF=KST,KPROF
    !nyc-pc     PB1(JROF,MSLB1SP9)=PB1(JROF,MSLB1SP9)+PGMVS(JROF,YT9%MSP)
    !nyc-pc   ENDDO
    !nyc-pc ENDIF

  ENDIF

  ! add preliminary quantity for LAPGHY physics

  IF( XIDT > 0.0_JPRB )THEN
    DO JROF=KST,KPROF
      PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)-ZXIGP*PB2(JROF,MSLB2SPSI)
      PB2(JROF,MSLB2SPSI)=ZXIGP*PB2(JROF,MSLB2SPSI)
    ENDDO
  ELSE
    DO JROF=KST,KPROF
      PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP)-PESGP*PB2(JROF,MSLB2SPSI)
      PB2(JROF,MSLB2SPSI)=PESGP*PB2(JROF,MSLB2SPSI)
    ENDDO
  ENDIF

ELSE

  ! TL of SL3TL currently not coded.
  CALL ABOR1(' LATTESTL: SL3TL not coded')

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LATTESTL',1,ZHOOK_HANDLE)
END SUBROUTINE LATTESTL

