#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE LATTE_NL(KSTEP, LDRPLANE, YDCST, YDGEOMETRY,YDML_DYN,YDRIP,&
 ! --- INPUT --------------------------------------------------
 & KST,KPROF,PDTS2,PBDT,&
 & PU,PV,PEVEL,PRDELP,&
 ! --- OUTPUT -------------------------------------------------
 & PMIXNL,KSETTLOFF)

!**** *LATTE_NL*   Semi-Lagrangian scheme.
!
!        Determine the way NL residuals are extrapolated in 2TL scheme
!        - stable points - SETTLS
!        - unstable points - NESC

!**   Interface.
!     ----------
!        *CALL* *LATTE_NL(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST         - first element of work.
!          KPROF       - depth of work.
!          PDTS2       - 0.5*time step for the first time-integration step of
!                        a leap-frog scheme or all time-integration steps of
!                        a two-time level scheme; time step for the following
!                        time-integration steps of a leap-frog scheme.
!          PBDT        - PDTS2*BETADT if semi-implicit scheme with unreduced
!                        divergence, PDTS2*BETADT*(c**2/GM**2) if semi-implicit
!                        scheme with reduced divergence.
!          PU          - U-wind
!          PV          - V-wind
!          PEVEL       - "etadot (d prehyd/d eta)" at half levels at t.
!          PRDELP      - 1/[Delta prehyd] at full levels.

!        INPUT/OUTPUT:
!          PMIXNL      - quantity used to control way we extrapolate quantities
!                        PMIXNL = 1 - SETTLS, PMIXNL = 0 - NESC
!          KSETTLOFF   - counter to accumulate informations about PMIXNL

!        OUTPUT:

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        none
!        Called by LACDYN.

!     Reference.
!     ----------
!        Arpege documentation about semi-Lagrangian scheme.

!     Author.
!     -------
!        J. Vivoda (SHMU/LACE)
!        Original : JUNE 2018

! Modifications
! -------------
! End Modifications
!------------------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE PARKIND1           , ONLY : JPIM, JPRB


USE YOMCST             , ONLY : TCST

USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE YOMRIP             , ONLY : TRIP

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER (KIND=JPIM), INTENT (IN) :: KSTEP
LOGICAL, INTENT (IN) :: LDRPLANE
TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY), INTENT(IN)       :: YDGEOMETRY
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN) :: YDML_DYN
TYPE(TRIP) ,    INTENT(IN)       :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBDT(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVEL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMIXNL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KSETTLOFF(YDGEOMETRY%YRDIMV%NFLEVG)
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JLEV, JROF, IOFF

REAL(KIND=JPRB) :: ZCFL, ZU, ZV, ZCFL_SETTLS, ZCFL_NESC, ZARG
REAL(KIND=JPRB) :: ZCFLV, ZETADOT, ZDETA
REAL(KIND=JPRB) :: ZWRL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     -------------------------------------------------------

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LATTE_NL',0,ZHOOK_HANDLE)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,   YDEGEO=>YDGEOMETRY%YREGEO,YDVETA=>YDGEOMETRY%YRVETA,&
& YDDYN=>YDML_DYN%YRDYN)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG,NCURRENT_ITER=>YDDYN%NCURRENT_ITER,TSTEP=>YDRIP%TSTEP,   EDELX=>YDEGEO%EDELX,&
& EDELY=>YDEGEO%EDELY,VETAH=>YDVETA%VETAH)

!     ------------------------------------------------------------------

!*      1.  PRELIMINARY INITIALISATIONS.
!       --------------------------------

! etadot at full levels
DO JLEV=1,NFLEVG
  DO JROF=KST,KPROF
    ZWRL0(JROF,JLEV)=&
     & 0.5_JPRB*(PEVEL(JROF,JLEV)+PEVEL(JROF,JLEV-1))&
     & *(VETAH(JLEV)-VETAH(JLEV-1))*PRDELP(JROF,JLEV)
  ENDDO
ENDDO


!*      2.  WEIGHTS CALCULATION.
!       ------------------------
IF (YDML_DYN%YRDYNA%LMIXETTLS .AND. NCURRENT_ITER == 0) THEN
  IF (LDRPLANE) THEN

    DO JLEV=1,NFLEVG
      IOFF = JLEV
      KSETTLOFF(IOFF) = 0
      ZDETA = VETAH(JLEV) - VETAH(JLEV-1)
      DO JROF=KST,KPROF
        ZU                      = ABS(PU(JROF, JLEV))
        ZV                      = ABS(PV(JROF, JLEV))
        ZETADOT                 = ABS(ZWRL0(JROF, JLEV))
        ZCFL                    = MAX(ZU * TSTEP / EDELX , ZV * TSTEP / EDELY)
        ZCFLV                   = ZETADOT * TSTEP / ZDETA
        ZCFL                    = MAX(ZCFLV, ZCFL)
        ZCFL_SETTLS             = 1.0_JPRB   ! under this value we use SETTLS
        ZCFL_NESC               = 1.01_JPRB  ! above this value we use NESC
        ZARG                    = (2.0_JPRB * ZCFL - ZCFL_SETTLS - ZCFL_NESC) /&
         &                        (ZCFL_SETTLS-ZCFL_NESC)
        PMIXNL(JROF, JLEV)      = 0.5_JPRB * (1.0_JPRB + TANH( YDCST%RPI * ZARG))
        PMIXNL(JROF, JLEV)      = MIN(1.0_JPRB, MAX(0.0_JPRB, PMIXNL(JROF, JLEV)))
        IF (PMIXNL(JROF, JLEV) < YDML_DYN%YRDYNA%RMIXNL_TRH) THEN
          KSETTLOFF(IOFF)    = KSETTLOFF(IOFF) + 1         
        ENDIF
      ENDDO
    ENDDO

  ELSE
    CALL ABOR1('LATTE_NL: GLOBAL GEOMETRY NOT YET CODED FOR MIXED TTL SCHEME.')
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LATTE_NL',1,ZHOOK_HANDLE)
END SUBROUTINE LATTE_NL

