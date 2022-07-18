SUBROUTINE LATTESAD(YDGEOMETRY,YDGMV,YDGMV5,YDRIP,YDML_DYN,KST,KPROF,PDTS2,PBDT,PESGP,PESGM,&
 & POROGL,POROGM,PSDIV0,PSDVBC,PRES0,PGMVS,&
 & PGMV,PGMVT1S,PB1,PB2,&
 & PSDVBC5,PRES5,PGMV5,PGMV5S)

!**** *LATTESAD* Semi-Lagrangian scheme.  (adjoint version)
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
!        *CALL* *LATTESAD(..)

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
!                        according to the value of CDCONF.
!          PESGP       - (1 + uncentering factor).
!          PESGM       - (1 - uncentering factor).
!          POROGL      - Zonal derivative of orography.
!          POROGM      - Meridian derivative of orography.

!        INPUT/OUTPUT:
!          PSDIV0      - SI term at time t for continuity equation (Nu*D).
!          PSDVBC      - vertical integral of divergence computed in "gpcty",
!                        including the "lrubc" and "delta m=1" contributions
!                        of (etadot d prehyd/d eta).
!          PRES0       - hydrostatic pressure at half levels at t.
!          PGMVS       - GMVS variables at t-dt and t.
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
!           Called by LACDYNAD.

!     Reference.
!     ----------
!             Arpege documentation about semi-Lagrangian scheme.

!     Author.
!     -------
!      C. Temperton (ECMWF)
!      Original : 99/10/12.

!     Modifications.
!     --------------
!      K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!      K. Yessad May 2009: rewrite to match the current direct code.
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
USE YOMCT0             , ONLY : LTWOTL

USE YOMCT3             , ONLY : NSTEP
USE YOMCVER            , ONLY : LVERTFE
USE YOMRIP             , ONLY : TRIP
USE YOMSTA             , ONLY : RTSUR

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
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSDIV0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSDVBC(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRES0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS)
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
#include "verintad.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LATTESAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA,   &
& YDDYN=>YDML_DYN%YRDYN,YDPTRSLB1=>YDML_DYN%YRPTRSLB1,    YDPTRSLB2=>YDML_DYN%YRPTRSLB2)

ASSOCIATE(NPROMA=>YDDIM%NPROMA,   NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA,   NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
& NVLAG=>YDDYN%NVLAG,   RCMSLP0=>YDDYN%RCMSLP0, XIDT=>YDDYN%XIDT,   NFOST=>YDRIP%NFOST,   YT0=>YDGMV%YT0,           &
& YT1=>YDGMV%YT1, YT9=>YDGMV%YT9,   YT5=>YDGMV5%YT5,   MSLB1C9=>YDPTRSLB1%MSLB1C9,   MSLB1SP9=>YDPTRSLB1%MSLB1SP9,  &
& MSLB2SPSI=>YDPTRSLB2%MSLB2SPSI)
!     ------------------------------------------------------------------

!############################################
! 1. AUXILIARITIES
!############################################

LLCT = YDML_DYN%YRDYNA%LPC_FULL .AND. NCURRENT_ITER > 0 ! corrector for LPC_FULL
LLCTC = YDML_DYN%YRDYNA%LPC_CHEAP .AND. NCURRENT_ITER > 0

ZXIDT0=1.0_JPRB+XIDT
ZXIDT9=1.0_JPRB+XIDT
ZXIGP=1.0_JPRB+XIDT

ZCMSLP=RCMSLP0/(RD*RTSUR)

! Initial zeroing of local quantities.
ZMOY1SP=0.0_JPRB
ZSPTB=0.0_JPRB
ZSPNLT_FE=0.0_JPRB
ZSPT1=0.0_JPRB

!     ------------------------------------------------------------------

!############################################
! 1b. TRAJECTORY
!############################################

!     ------------------------------------------------------------------

!############################################
! 2. NONLINEAR MODEL
!############################################

!*       *   Continuity equation.
    
IF (LTWOTL) THEN

  ! add preliminary quantity for LAPGHY physics

  IF( XIDT > 0.0_JPRB )THEN
    DO JROF=KST,KPROF
      PB2(JROF,MSLB2SPSI)=ZXIGP*PB2(JROF,MSLB2SPSI)
      PB2(JROF,MSLB2SPSI)=PB2(JROF,MSLB2SPSI)-ZXIGP*PGMVT1S(JROF,YT1%MSP)
    ENDDO
  ELSE 
    DO JROF=KST,KPROF
      PB2(JROF,MSLB2SPSI)=PESGP*PB2(JROF,MSLB2SPSI)
      PB2(JROF,MSLB2SPSI)=PB2(JROF,MSLB2SPSI)-PESGP*PGMVT1S(JROF,YT1%MSP)
    ENDDO
  ENDIF

  IF(.NOT.LLCT)THEN

    !############################################
    ! Predictor for LPC_FULL
    ! or case nsiter=0.
    !############################################

    IF (YDML_DYN%YRDYNA%LPC_CHEAP.AND.(.NOT.YDML_DYN%YRDYNA%LNESC)) THEN
      ! LPC_CHEAP=T is currently coded for LNESC only
      !  (for LNESC=F interpolated quantities are not the same ones
      !  at the predictor and corrector steps, and the LPC_CHEAP code
      !  currently assumes that they are identical).
      CALL ABOR1('LATTESAD: If LPC_CHEAP=T, LNESC should be T')
    ENDIF

    !------------------------------------------------------------------------
    ! 2D buffers
    !------------------------------------------------------------------------
    IF(NVLAG == 2 .OR. NVLAG == 3)THEN
      DO JROF=KST,KPROF
        PGMVS(JROF,YT0%MSP)=PGMVS(JROF,YT0%MSP)+PB1(JROF,MSLB1SP9)
      ENDDO
    ENDIF

    !------------------------------------------------------------------------
    ! 3D additional actions
    !------------------------------------------------------------------------
    IF (LVERTFE) THEN
      IF(NVLAG == 2 .OR. NVLAG == 3)THEN
        DO JROF=KST,KPROF
          ZSPT1(JROF,NFLEVG+1)=ZSPT1(JROF,NFLEVG+1)+PGMVT1S(JROF,YT1%MSP)
        ENDDO
        CALL VERINTAD(NPROMA,KST,KPROF,NFLEVG,ZSPNLT_FE,ZSPT1,0,YDGEOMETRY%YRVFE)
      ENDIF
    ENDIF

    !---------------------------------------------------------------
    ! 3D buffers
    !---------------------------------------------------------------
    DO JLEV=1,NFLEVG

      ZSPNLT0=0.0_JPRB

      IF( .NOT.YDML_DYN%YRDYNA%LNESC )THEN
        ! save of nonlinear residual at time t
        ! to be used as nonlinear residual at time t-dt next time step
        DO JROF=KST,KPROF
          ZMOY1SP(JROF,JLEV)=ZMOY1SP(JROF,JLEV)+PGMV(JROF,JLEV,YT9%MSPNL)
          PB2(JROF,MSLB2SPSI)=PB2(JROF,MSLB2SPSI)&
           & +ZXIDT9*PGMV(JROF,JLEV,YT9%MSPNL)
          PGMV(JROF,JLEV,YT9%MSPNL)=0.0_JPRB
        ENDDO
      ENDIF

      ! Save quantities for corrector step
      !nyc-pc  IF( LPC_FULL )THEN
      !nyc-pc    ! save nonlinear model at time t
      !nyc-pc    DO JROF=KST,KPROF
      !nyc-pc      ZMOY1SP(JROF,JLEV)=ZMOY1SP(JROF,JLEV)+PGMV(JROF,JLEV,YT9%MCSPNL)
      !nyc-pc    ENDDO
      !nyc-pc  ENDIF

      ! Fill PB1(.,MSLB1C9),PGMVT1S(.,YT1%MSP),ZSPNLT_FE:
      IF (NVLAG == 2 .OR. NVLAG == 3) THEN

        IF ( NSTEP <= NFOST .OR. YDML_DYN%YRDYNA%LNESC ) THEN

          IF (LVERTFE) THEN
            DO JROF=KST,KPROF
              ZSPNLT0(JROF)=ZSPNLT0(JROF)&
               & +YDVAB%VDELB(JLEV)*YDVETA%VFE_RDETAH(JLEV)*PESGP*ZSPNLT_FE(JROF,JLEV)
            ENDDO
          ELSE
            DO JROF=KST,KPROF
              ZSPNLT0(JROF)=ZSPNLT0(JROF)&
               & +PESGP*YDVAB%VDELB(JLEV)*PGMVT1S(JROF,YT1%MSP)
            ENDDO
          ENDIF
          DO JROF=KST,KPROF
            ZMOY1SP(JROF,JLEV)=ZMOY1SP(JROF,JLEV)&
             & +PESGM*PB1(JROF,MSLB1C9+JLEV-NFLSA)
          ENDDO

        ELSEIF (YDML_DYN%YRDYNA%LSETTLS) THEN

          IF( LVERTFE )THEN
            DO JROF=KST,KPROF
              ZSPNLT0(JROF)=ZSPNLT0(JROF)&
               & +PESGP*YDVAB%VDELB(JLEV)*YDVETA%VFE_RDETAH(JLEV)*ZSPNLT_FE(JROF,JLEV)
            ENDDO
          ELSE
            DO JROF=KST,KPROF
              ZSPNLT0(JROF)=ZSPNLT0(JROF)&
               & +PESGP*YDVAB%VDELB(JLEV)*PGMVT1S(JROF,YT1%MSP)
            ENDDO
          ENDIF
          DO JROF=KST,KPROF
            ZMOY1SP(JROF,JLEV)=ZMOY1SP(JROF,JLEV)&
             & +PESGM*PB1(JROF,MSLB1C9+JLEV-NFLSA)
            ZSPNLT0(JROF)=ZSPNLT0(JROF)+PB1(JROF,MSLB1C9+JLEV-NFLSA)
            PGMV(JROF,JLEV,YT9%MSPNL)=PGMV(JROF,JLEV,YT9%MSPNL)&
             & -PB1(JROF,MSLB1C9+JLEV-NFLSA)
          ENDDO

        ELSE

          ! remaining case: ldsettls=false, lnesc=false.
          IF (LVERTFE) THEN
            DO JROF=KST,KPROF
              ZSPNLT0(JROF)=ZSPNLT0(JROF)+1.5_JPRB*&
               & YDVAB%VDELB(JLEV)*YDVETA%VFE_RDETAH(JLEV)*PESGP*ZSPNLT_FE(JROF,JLEV)
              PGMV(JROF,JLEV,YT9%MSPNL)=PGMV(JROF,JLEV,YT9%MSPNL)-0.5_JPRB*&
               & YDVAB%VDELB(JLEV)*YDVETA%VFE_RDETAH(JLEV)*PESGP*ZSPNLT_FE(JROF,JLEV)
            ENDDO
          ELSE 
            DO JROF=KST,KPROF
              ZSPNLT0(JROF)=ZSPNLT0(JROF)&
               & +1.5_JPRB*PESGP*YDVAB%VDELB(JLEV)*PGMVT1S(JROF,YT1%MSP)
              PGMV(JROF,JLEV,YT9%MSPNL)=PGMV(JROF,JLEV,YT9%MSPNL)&
               & -0.5_JPRB*PESGP*YDVAB%VDELB(JLEV)*PGMVT1S(JROF,YT1%MSP)
            ENDDO
          ENDIF

          DO JROF=KST,KPROF
            ZMOY1SP(JROF,JLEV)=ZMOY1SP(JROF,JLEV)&
             & +PESGM*PB1(JROF,MSLB1C9+JLEV-NFLSA)
            ZSPNLT0(JROF)=ZSPNLT0(JROF)&
             & +0.5_JPRB*PESGM*PB1(JROF,MSLB1C9+JLEV-NFLSA)
            PGMV(JROF,JLEV,YT9%MSPNL)=PGMV(JROF,JLEV,YT9%MSPNL)&
             & -0.5_JPRB*PESGM*PB1(JROF,MSLB1C9+JLEV-NFLSA)
          ENDDO

        ENDIF

      ENDIF

      ! nonlinear residual time t
      DO JROF=KST,KPROF
        ZMOY1SP(JROF,JLEV)=ZMOY1SP(JROF,JLEV)+ZSPNLT0(JROF)
        PB2(JROF,MSLB2SPSI)=PB2(JROF,MSLB2SPSI)+ZXIDT0*ZSPNLT0(JROF)
      ENDDO

    ENDDO

  ELSE

    !############################################
    ! Corrector for LPC_FULL
    !############################################

    ! ky: piece of code not yet validated.
    !     Commented for the time being.

    !------------------------------------------------------------------------
    ! 2D buffers
    !------------------------------------------------------------------------
    !nyc-pc IF (.NOT.LLCTC) THEN
    !nyc-pc   DO JROF=KST,KPROF
    !nyc-pc     PGMVS(JROF,YT9%MSP)=PGMVS(JROF,YT9%MSP)+PB1(JROF,MSLB1SP9)
    !nyc-pc   ENDDO
    !nyc-pc ENDIF

    !------------------------------------------------------------------------
    ! 3D additional actions
    !------------------------------------------------------------------------
    !nyc-pc IF (LVERTFE) THEN
    !nyc-pc   DO JROF=KST,KPROF
    !nyc-pc     ZSPT1(JROF,NFLEVG+1)=ZSPT1(JROF,NFLEVG+1)+PGMVT1S(JROF,YT1%MSP)
    !nyc-pc   ENDDO
    !nyc-pc   CALL VERINTAD(NPROMA,KST,KPROF,NFLEVG,ZSPNLT_FE,ZSPT1,0)
    !nyc-pc ENDIF

    !---------------------------------------------------------------
    ! 3D buffers
    !---------------------------------------------------------------
    !nyc-pc DO JLEV=1,NFLEVG

      !nyc-pc ZSPNLT1=0.0_JPRB

      !nyc-pc ! nonlinear residual time t

      !nyc-pc IF (LVERTFE) THEN
      !nyc-pc   DO JROF=KST,KPROF
      !nyc-pc     ZSPNLT1(JROF)=ZSPNLT1(JROF) &
      !nyc-pc      & +YDVAB%VDELB(JLEV)*YDVETA%VFE_RDETAH(JLEV)*PESGP*ZSPNLT_FE(JROF,JLEV)
      !nyc-pc   ENDDO
      !nyc-pc ELSE
      !nyc-pc   DO JROF=KST,KPROF
      !nyc-pc     PGMVT1S(JROF,YT1%MSP)=PGMVT1S(JROF,YT1%MSP) &
      !nyc-pc      & +YDVAB%VDELB(JLEV)*PESGP*ZSPNLT1(JROF)
      !nyc-pc   ENDDO
      !nyc-pc ENDIF
      !nyc-pc IF (.NOT.LLCTC) THEN
      !nyc-pc   DO JROF=KST,KPROF
      !nyc-pc     PGMV(JROF,JLEV,YT9%MCSPNL)=PGMV(JROF,JLEV,YT9%MCSPNL) &
      !nyc-pc      & +PESGM*PB1(JROF,MSLB1C9+JLEV-NFLSA)
      !nyc-pc   ENDDO
      !nyc-pc ENDIF

      !nyc-pc DO JROF=KST,KPROF
      !nyc-pc ZMOY1SP(KST:KPROF,JLEV)=ZMOY1SP(KST:KPROF,JLEV) &
      !nyc-pc  & +ZSPNLT1(KST:KPROF)
      !nyc-pc PB2(KST:KPROF,MSLB2SPSI)=PB2(KST:KPROF,MSLB2SPSI) &
      !nyc-pc  & +ZXIDT9*ZSPNLT1(KST:KPROF)
      !nyc-pc ENDDO

    !nyc-pc ENDDO

  ENDIF

ELSE
  
  ! AD of SL3TL currently not coded.
  CALL ABOR1(' LATTESAD: SL3TL not coded')

ENDIF

! nonlinear model at time t

DO JLEV=1,NFLEVG
  DO JROF=KST,KPROF
    ZSPTB(JROF)=ZSPTB(JROF)+ZMOY1SP(JROF,JLEV)
  ENDDO
  DO JROF=KST,KPROF
    PGMVS(JROF,YT0%MSPL)=PGMVS(JROF,YT0%MSPL)&
     & +PDTS2*PGMV5(JROF,JLEV,YT5%MU)*ZMOY1SP(JROF,JLEV)
    PGMV(JROF,JLEV,YT0%MU)=PGMV(JROF,JLEV,YT0%MU)&
     & +PDTS2*(PGMV5S(JROF,YT5%MSPL)+ZCMSLP*POROGL(JROF))*ZMOY1SP(JROF,JLEV)
    PGMVS(JROF,YT0%MSPM)=PGMVS(JROF,YT0%MSPM)&
     & +PDTS2*PGMV5(JROF,JLEV,YT5%MV)*ZMOY1SP(JROF,JLEV)
    PGMV(JROF,JLEV,YT0%MV)=PGMV(JROF,JLEV,YT0%MV)&
     & +PDTS2*(PGMV5S(JROF,YT5%MSPM)+ZCMSLP*POROGM(JROF))*ZMOY1SP(JROF,JLEV)
  ENDDO
ENDDO

! linear model at time t
DO JROF=KST,KPROF
  ZSPTB(JROF)=-PDTS2*ZSPTB(JROF)
  PSDIV0(JROF)=PSDIV0(JROF)+PBDT(JROF)*PB2(JROF,MSLB2SPSI)
ENDDO
PB2(KST:KPROF,MSLB2SPSI)=0.0_JPRB

DO JROF=KST,KPROF
  ZRPRES5=1.0_JPRB/PRES5(JROF,NFLEVG)
  PSDVBC(JROF,NFLEVG)=PSDVBC(JROF,NFLEVG)+ZRPRES5*ZSPTB(JROF)
  PRES0(JROF,NFLEVG)=PRES0(JROF,NFLEVG)&
   & -ZRPRES5*ZRPRES5*PSDVBC5(JROF,NFLEVG)*ZSPTB(JROF)
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LATTESAD',1,ZHOOK_HANDLE)
END SUBROUTINE LATTESAD

