SUBROUTINE LATTEX_DNT_AD(YDGEOMETRY,YDRIP,YDDYN,KST,KPROF,LDSETTLS,KXLAG,PESGP,PESGM,PXT0, &
 !nyc-pc & PXT9, &
 & PMOY1X,PXSI,PXNLT9,PXT1,PXL0,PXL9,&
 !nyc-pc & PCXNLT9, &
 & LDNESC)

!------------------------------------------------------------------------------
! LATTEX_DNT_AD - Semi-Lagrangian scheme.
!               Computation of the t and t-dt useful quantities
!               at grid-points. Equations for tri-dimensional
!               variables for 2TL scheme.
!               Adjoint code.

! Purpose
! -------
!       Adjoint code of LATTEX_DNT: see LATTEX_DNT for more information.
!       Some options present in LATTEX_DNT are not coded or not completely
!        validated in LATTEX_DNT_AD (for ex: LPC_FULL, KXLAG=2).

! Interface
! ---------
!   CALL LATTEX_DNT_AD(..)

! Explicit arguments :
! --------------------

!          KST     - first element of work.
!          KPROF   - depth of work.
!          LDSETTLS- .T./.F.: Stable/Conventional extrapolations for SL2TL.
!          KXLAG   - type of SL discretisation 
!          PESGP   - (1 + uncentering factor).
!          PESGM   - (1 - uncentering factor).
!          PXT0    - prognostic variable time t (predictor or SI),
!                    preliminary t+dt (corrector)
!          PXT9    - prognostic variable time t (corrector),
!                    not used for predictor or SI
!          PMOY1X  - full nonlinear model at time t [(Delta t/2) "cursive" A]
!          PXSI    - semi-implicit linear model at time t
!                    [- (Delta t/2) "cursive" B]
!          PXNLT9  - buffer used during predictor resp. SI time step
!          PXT1    - t+dt term and other final point terms 
!          PXL0    - second SL quantity to be interpolated (linearly) at O
!                    if KXLAG=3
!          PXL9    - SL quantity to be interpolated at O 
!          PCXNLT9 - buffer used during corrector
!          LDNESC  - non-extrapolating option for SL2TL.

! Externals
! ---------
!   none
!   Called by LATTEXAD.

! Method
! ------

! Reference
! ---------
!   Arpege documentation about semi-Lagrangian scheme.

! Author
! ------
!      K. Yessad (March 2007) after LATTEXAD.

! Modifications
! -------------
!  K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!  K. Yessad May 2009: complete recoding according to LATTEX_DNT-CY36.
!  K. Yessad (Oct 2013): allow NESC without ICI-PC scheme.
!  K. Yessad (July 2014): Move some variables, rename some variables.
! End Modifications
!------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMDYNA      , ONLY : LNESC, LPC_FULL, LPC_CHEAP
USE YOMCT3       , ONLY : NSTEP
USE YOMRIP       , ONLY : TRIP
USE YOMDYN       , ONLY : TDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)            ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)                ,INTENT(IN)    :: YDDYN
TYPE(TRIP)                ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM)        ,INTENT(IN)    :: KST 
INTEGER(KIND=JPIM)        ,INTENT(IN)    :: KPROF 
LOGICAL                   ,INTENT(IN)    :: LDSETTLS 
INTEGER(KIND=JPIM)        ,INTENT(IN)    :: KXLAG 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PESGP 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PESGM 
REAL(KIND=JPRB)           ,INTENT(INOUT) :: PXT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
!nyc-pc REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXT9(NPROMA,NFLEVG)
REAL(KIND=JPRB)           ,INTENT(INOUT) :: PMOY1X(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(INOUT) :: PXSI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(INOUT) :: PXNLT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PXT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PXL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN) 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PXL9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN) 
!nyc-pc REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCXNLT9(NPROMA,NFLEVG)
LOGICAL,OPTIONAL  ,INTENT(IN)    :: LDNESC

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV,JROF

! * ZXNLT0 (resp. ZXNLT1) the non linear term at t (resp. t+dt).
!nyc-pc REAL(KIND=JPRB) :: ZXNLT1(NPROMA)
REAL(KIND=JPRB) :: ZXNLT0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZXIDT0,ZXIDT9,ZXIGP

LOGICAL :: LLCT,LLNESC

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LATTEX_DNT_AD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & NFOST=>YDRIP%NFOST, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, XIDT=>YDDYN%XIDT)
!     ------------------------------------------------------------------

!*      1. AUXILIARY VARIABLES.
!       -----------------------

IF (PRESENT(LDNESC)) THEN
  LLNESC=LDNESC
ELSE
  LLNESC=LNESC
ENDIF

LLCT = (LPC_FULL .AND. NCURRENT_ITER > 0)

ZXIDT0=1.0_JPRB+XIDT
ZXIDT9=1.0_JPRB+XIDT
ZXIGP=1.0_JPRB+XIDT

!     ------------------------------------------------------------------

!*      2. MAIN CALCULATIONS.
!       ---------------------

!########################################################
! 2.3 Addition of preliminary quantity for LAGPHY physics
!########################################################

DO JLEV=1,NFLEVG
  IF(XIDT > 0.0_JPRB)THEN
    DO JROF=KST,KPROF
      PXSI(JROF,JLEV)=ZXIGP*PXSI(JROF,JLEV)
      PXSI(JROF,JLEV)=PXSI(JROF,JLEV)-ZXIGP*PXT1(JROF,JLEV)
    ENDDO
  ELSE
    DO JROF=KST,KPROF
      PXSI(JROF,JLEV)=PESGP*PXSI(JROF,JLEV)
      PXSI(JROF,JLEV)=PXSI(JROF,JLEV)-PESGP*PXT1(JROF,JLEV)
    ENDDO
  ENDIF
ENDDO

IF( .NOT.LLCT )THEN

  !############################################
  ! 2.1 Predictor for LPC_FULL,  
  !     or case NSITER=0.
  !############################################

  IF (LPC_CHEAP.AND.(.NOT.LLNESC)) THEN
    ! LPC_CHEAP=T is currently coded for LNESC only
    !  (for LNESC=F interpolated quantities are not the same ones
    !  at the predictor and corrector steps, and the LPC_CHEAP code
    !  currently assumes that they are identical). 
    CALL ABOR1('LATTEX_DNT_AD: If LPC_CHEAP=T, LNESC should be T')
  ENDIF

  DO JLEV=1,NFLEVG

    ! Preliminary fill with zeros.
    ZXNLT0(KST:KPROF)=0.0_JPRB
    PMOY1X(KST:KPROF,JLEV)=0.0_JPRB

    IF( .NOT.LLNESC )THEN
      ! save of nonlinear residual at time t
      ! to be used as nonlinear residual at time t-dt next time step
      DO JROF=KST,KPROF
        PMOY1X(JROF,JLEV)=PMOY1X(JROF,JLEV)+PXNLT9(JROF,JLEV)
        PXSI(JROF,JLEV)=PXSI(JROF,JLEV)+ZXIDT9*PXNLT9(JROF,JLEV)
      ENDDO
    ENDIF
    PXNLT9(KST:KPROF,JLEV)=0.0_JPRB

    ! save quantities for corrector step
    !nyc-pc IF( LPC_FULL )THEN
    !nyc-pc   ! save nonlinear model at time t
    !nyc-pc   PMOY1X(KST:KPROF,JLEV)=PMOY1X(KST:KPROF,JLEV)+PCXNLT9(KST:KPROF,JLEV)
    !nyc-pc ENDIF

    ! Fill PXL9,PXL0,PXT1.
    IF (KXLAG == 2) THEN     !! AD NOT YET DONE !!

    ELSEIF (KXLAG == 3) THEN
      DO JROF=KST,KPROF
        PXT0(JROF,JLEV)=PXT0(JROF,JLEV)+PXL9(JROF,JLEV)
      ENDDO
      IF (NSTEP <= NFOST .OR. LLNESC) THEN
        DO JROF=KST,KPROF
          ZXNLT0(JROF)=ZXNLT0(JROF)+PESGP*PXT1(JROF,JLEV)
          PMOY1X(JROF,JLEV)=PMOY1X(JROF,JLEV)+PESGM*PXL0(JROF,JLEV)
        ENDDO
      ELSEIF (LDSETTLS) THEN
        DO JROF=KST,KPROF
          ZXNLT0(JROF)=ZXNLT0(JROF)+PESGP*PXT1(JROF,JLEV)
          PMOY1X(JROF,JLEV)=PMOY1X(JROF,JLEV)+PESGM*PXL0(JROF,JLEV)
          ZXNLT0(JROF)=ZXNLT0(JROF)+PXL0(JROF,JLEV)
          PXNLT9(JROF,JLEV)=PXNLT9(JROF,JLEV)-PXL0(JROF,JLEV)
        ENDDO
      ELSE
        ! * remaining case: ldsettls=false, lnesc=false.
        DO JROF=KST,KPROF
          ZXNLT0(JROF)=ZXNLT0(JROF)+1.5_JPRB*PESGP*PXT1(JROF,JLEV)
          PXNLT9(JROF,JLEV)=PXNLT9(JROF,JLEV)-0.5_JPRB*PESGP*PXT1(JROF,JLEV)
          PMOY1X(JROF,JLEV)=PMOY1X(JROF,JLEV)+PESGM*PXL0(JROF,JLEV)
          ZXNLT0(JROF)=ZXNLT0(JROF)+0.5_JPRB*PESGM*PXL0(JROF,JLEV)
          PXNLT9(JROF,JLEV)=PXNLT9(JROF,JLEV)-0.5_JPRB*PESGM*PXL0(JROF,JLEV)
        ENDDO
      ENDIF
    ENDIF

    ! nonlinear residual time t
    DO JROF=KST,KPROF
      PMOY1X(JROF,JLEV)=PMOY1X(JROF,JLEV)+ZXNLT0(JROF)
      PXSI(JROF,JLEV)=PXSI(JROF,JLEV)+ZXIDT0*ZXNLT0(JROF)
    ENDDO

  ENDDO

ELSE

  !############################################
  ! 2.2 Corrector for LPC_FULL
  !############################################

  ! ky: piece of code not yet validated.
  !     Commented for the time being.

  !nyc-pc DO JLEV=1,NFLEVG

  !nyc-pc   ! Preliminary fill with zeros.
  !nyc-pc   ZXNLT1(KST:KPROF)=0.0_JPRB

  !nyc-pc   IF (KXLAG == 2) THEN     !! AD NOT YET DONE !!

  !nyc-pc   ELSEIF (KXLAG == 3) THEN
  !nyc-pc     DO JROF=KST,KPROF
  !nyc-pc       ZXNLT1(JROF)=ZXNLT1(JROF)+PESGP*PXT1(JROF,JLEV)
  !nyc-pc     ENDDO
  !nyc-pc     IF (.NOT. LPC_CHEAP) THEN
  !nyc-pc       DO JROF=KST,KPROF
  !nyc-pc         PXT9(JROF,JLEV)=PXT9(JROF,JLEV)+PXL9(JROF,JLEV)
  !nyc-pc         PCXNLT9(JROF,JLEV)=PCXNLT9(JROF,JLEV)+PESGM*PXL0(JROF,JLEV)
  !nyc-pc       ENDDO
  !nyc-pc     ENDIF
  !nyc-pc   ENDIF

  !nyc-pc   ! nonlinear residual time t+dt (preliminary state from predictor)
  !nyc-pc   DO JROF=KST,KPROF
  !nyc-pc     PMOY1X(JROF,JLEV)=PMOY1X(JROF,JLEV)+ZXNLT1(JROF)
  !nyc-pc     PXSI(JROF,JLEV)=PXSI(JROF,JLEV)+ZXIDT9*ZXNLT1(JROF) 
  !nyc-pc   ENDDO

  !nyc-pc ENDDO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LATTEX_DNT_AD',1,ZHOOK_HANDLE)
END SUBROUTINE LATTEX_DNT_AD

