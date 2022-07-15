SUBROUTINE LANHSIB(YDGEOMETRY,YDVARS,YDCPG_SL1,YDCPG_SL2,&
 ! --- INPUT ------------------------------------------------------------------
 & YDML_DYN,KST,KPROF,PBT,PBDT,PESGP,&
 & PGAGT0L,PGAGT0M,PTOD0,PSPDS0,PLPD0,PSDIV0,PNHXT0,PMIXNL)

!**** *LANHSIB*   Semi-Lagrangian scheme.
!        Storage and memory transfers for linear terms.
!        NH models with LGWADV=T.

!     Purpose.
!     --------
!        Storage and memory transfers for linear terms.
!        NH models with LGWADV=T.

!**   Interface.
!     ----------
!        *CALL* *LANHSIB(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST      - first element of work.
!          KPROF    - depth of work.
!          PBT      - PDTS2*BETADT (BETADT is in YOMDYN).
!          PBDT     - PBT if semi-implicit scheme with unreduced
!                        divergence, PBT*(c**2/GM**2) if semi-implicit
!                        scheme with reduced divergence.
!          PESGP    - first  order decentering factor
!          PGAGT0L  - semi-implicit term at time t for U-wind equation.
!          PGAGT0M  - semi-implicit term at time t for V-wind equation.
!          PTOD0    - semi-implicit term at time t for temperature equation.
!          PSPDS0   - semi-implicit term at time t for "pressure departure variable" equation (NHEE).
!          PLPD0    - semi-implicit term at time t for "vertical divergence variable" equation.
!          PSDIV0   - semi-implicit term at time t for "ln(pi_s)" equation.
!          PNHXT0   - "Xterm" on layers at t.
!          PMIXNL   - quantity used to control way we extrapolate quantities.

!        INPUT-OUTPUT:
!          PGMV     - GMV variables.
!          PB1      - SLBUF1 buffer for interpolations.

!        OUTPUT:
!          PB2      - SLBUF2 buffer.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation.

!     Externals.
!     ----------
!           Called by LACDYN.

!     Reference.
!     ----------
!             Arpege documentation about semi-Lagrangian scheme.

!     Author.
!     -------
!        J. Vivoda (SHMI/LACE) 
!        Original : Feb 2005

!     Modifications.
!     --------------
!   N. Wedi and K. Yessad (Nov 2007): ND4SYS to improve NH model stability
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (July 2014): Move some variables.
!   K. Yessad (June 2017): Introduce NHQE model.
!   J. Vivoda (July 2018): mixed NESC/SETTLS scheme.
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_SL1_TYPE, CPG_SL2_TYPE
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE YOMGMV             , ONLY : TGMV
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE YOMCT0             , ONLY : LNHEE
USE YOMCT3             , ONLY : NSTEP
USE YOMDYNA            , ONLY : NVDVAR, ND4SYS, LSLINL, LSLINLC2
USE FIELD_VARIABLES_MOD,ONLY : FIELD_VARIABLES

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(CPG_SL1_TYPE)    ,INTENT(INOUT)  :: YDCPG_SL1
TYPE(CPG_SL2_TYPE)    ,INTENT(INOUT)  :: YDCPG_SL2
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBDT(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESGP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAGT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAGT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTOD0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDS0(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLPD0(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSDIV0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHXT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMIXNL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPRB)    :: ZBT,ZBDT(YDGEOMETRY%YRDIM%NPROMA),ZMIXNL
REAL(KIND=JPRB)    :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LANHSIB',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDDYN=>YDML_DYN%YRDYN,YDPTRSLB1=>YDML_DYN%YRPTRSLB1,&
& YDPTRSLB2=>YDML_DYN%YRPTRSLB2)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, NCURRENT_ITER=>YDDYN%NCURRENT_ITER, XIDT=>YDDYN%XIDT)
!     ------------------------------------------------------------------

!     1. SEMI-IMPLICIT LINEAR MODEL CORRECTED WITH TIME STEP AND
!       CENTERING FACTORS
!     ----------------------------------------------------------

IF (XIDT > 0.0_JPRB)THEN
  ZBT=(1.0_JPRB+XIDT)*PBT
  ZBDT(KST:KPROF)=(1.0_JPRB+XIDT)*PBDT(KST:KPROF)
ELSE
  ZBT=PBT
  ZBDT(KST:KPROF)=PBDT(KST:KPROF)
ENDIF

! * multiplication by time step
DO JLEV = 1,NFLEVG
  DO JROF = KST,KPROF
    YDCPG_SL2%USI(JROF,JLEV)=ZBT*PGAGT0L(JROF,JLEV)
    YDCPG_SL2%VSI(JROF,JLEV)=ZBT*PGAGT0M(JROF,JLEV)
    YDCPG_SL2%TSI(JROF,JLEV)=ZBDT(JROF)*PTOD0(JROF,JLEV)
    YDCPG_SL2%VDSI(JROF,JLEV)=ZBT*PLPD0(JROF,JLEV)
  ENDDO
  IF (LNHEE) THEN
    DO JROF = KST,KPROF
      YDCPG_SL2%PDSI(JROF,JLEV)=ZBDT(JROF)*PSPDS0(JROF,JLEV)
    ENDDO
  ENDIF
ENDDO
DO JROF = KST,KPROF
  YDCPG_SL2%SPSI(JROF)=ZBDT(JROF)*PSDIV0(JROF)
ENDDO

! * first-order decentering factor
IF(PESGP > 1.0_JPRB)THEN
  DO JROF = KST,KPROF
    YDCPG_SL2%SPSI(JROF)=PESGP*YDCPG_SL2%SPSI(JROF)
  ENDDO
  DO JLEV=1,NFLEVG
    DO JROF = KST,KPROF
      YDCPG_SL2%USI(JROF,JLEV)=PESGP*YDCPG_SL2%USI(JROF,JLEV)
      YDCPG_SL2%VSI(JROF,JLEV)=PESGP*YDCPG_SL2%VSI(JROF,JLEV)
      YDCPG_SL2%TSI(JROF,JLEV)=PESGP*YDCPG_SL2%TSI(JROF,JLEV)
      YDCPG_SL2%VDSI(JROF,JLEV)=PESGP*YDCPG_SL2%VDSI(JROF,JLEV)
    ENDDO
    IF (LNHEE) THEN
      DO JROF = KST,KPROF
        YDCPG_SL2%PDSI(JROF,JLEV)=PESGP*YDCPG_SL2%PDSI(JROF,JLEV)
      ENDDO
    ENDIF
  ENDDO
ENDIF

! * add NHX-term to the equation for 
!   (K.Y.: I have serious doubts about the multiplication of NHX by PESGP).
IF( NVDVAR == 4 .OR. NVDVAR == 5 )THEN
  IF (ND4SYS==1 .OR. (ND4SYS==2 .AND. NCURRENT_ITER==0)) THEN
    DO JLEV=1,NFLEVG
      DO JROF = KST,KPROF
        YDCPG_SL2%VDSI(JROF,JLEV)=YDCPG_SL2%VDSI(JROF,JLEV)&
         & +PESGP*PNHXT0(JROF,JLEV)
      ENDDO
    ENDDO
  ELSEIF (ND4SYS==2 .AND. NCURRENT_ITER>0) THEN
    DO JLEV=1,NFLEVG
      DO JROF = KST,KPROF
        YDCPG_SL2%VDSI(JROF,JLEV)=YDCPG_SL2%VDSI(JROF,JLEV)&
         & +PESGP*YDVARS%NHX%T9(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!     2. FILL PB1 WHEN REQUIRED 
!     -------------------------

! If NSTEP <= 0, values at t-dt are equal to values a 't', what is added in PB1
! is simply 0 in this case.
IF (NSTEP > 0 .AND. NCURRENT_ITER==0) THEN
  ! * Put 0.5 Dt betadt (Lin(t)-Lin(t-dt)) if XIDT=0 in PB1
  !   Put 0.5 Dt (1+XIDT) betadt (Lin(t)-Lin(t-dt)) in PB1 if XIDT>0
  IF (LSLINLC2) THEN
    DO JLEV = 1,NFLEVG
      DO JROF = KST,KPROF
        ZMIXNL = PMIXNL(JROF, JLEV)
        YDCPG_SL1%C9_SI(JROF,JLEV)=YDCPG_SL1%C9_SI(JROF,JLEV)&
         & +ZMIXNL*(ZBDT(JROF)*PSDIV0(JROF)-YDVARS%SPNL_SI%T9(JROF,JLEV))
      ENDDO
    ENDDO
  ENDIF
  IF (LSLINL) THEN
    DO JLEV = 1,NFLEVG
      DO JROF = KST,KPROF
        ZMIXNL = PMIXNL(JROF, JLEV)
        YDCPG_SL1%U9_SI(JROF,JLEV)=YDCPG_SL1%U9_SI(JROF,JLEV)&
         & + ZMIXNL*(ZBT*PGAGT0L(JROF,JLEV)-YDVARS%UNL_SI%T9(JROF,JLEV))
        YDCPG_SL1%V9_SI(JROF,JLEV)=YDCPG_SL1%V9_SI(JROF,JLEV)&
         & + ZMIXNL*(ZBT*PGAGT0M(JROF,JLEV)-YDVARS%VNL_SI%T9(JROF,JLEV))
        YDCPG_SL1%T9_SI(JROF,JLEV)=YDCPG_SL1%T9_SI(JROF,JLEV)&
         & + ZMIXNL*(ZBDT(JROF)*PTOD0(JROF,JLEV)-YDVARS%TNL_SI%T9(JROF,JLEV))
        YDCPG_SL1%VD9_SI(JROF,JLEV)=YDCPG_SL1%VD9_SI(JROF,JLEV)&
         & + ZMIXNL*(ZBT*PLPD0(JROF,JLEV)-YDVARS%SVDNL_SI%T9(JROF,JLEV))
      ENDDO
      IF (LNHEE) THEN
        DO JROF = KST,KPROF
          ZMIXNL = PMIXNL(JROF, JLEV)
          YDCPG_SL1%PD9_SI(JROF,JLEV)=YDCPG_SL1%PD9_SI(JROF,JLEV)&
           & + ZMIXNL*(ZBDT(JROF)*PSPDS0(JROF,JLEV)-YDVARS%SPDNL_SI%T9(JROF,JLEV))
        ENDDO
      ENDIF
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!     3. SAVE 't' VALUES WHEN REQUIRED 
!     --------------------------------

! * Save 't' values of linear terms (with 'beta Dt' but without uncentering
!   factor) when required:
IF (LSLINLC2 .AND. NCURRENT_ITER==0) THEN
  DO JLEV = 1,NFLEVG
    DO JROF = KST,KPROF
      YDVARS%SPNL_SI%T9(JROF,JLEV)=ZBDT(JROF)*PSDIV0(JROF)
    ENDDO
  ENDDO
ENDIF
IF (LSLINL .AND. NCURRENT_ITER==0) THEN
  DO JLEV = 1,NFLEVG
    DO JROF = KST,KPROF
      YDVARS%UNL_SI%T9(JROF,JLEV)=ZBT*PGAGT0L(JROF,JLEV)
      YDVARS%VNL_SI%T9(JROF,JLEV)=ZBT*PGAGT0M(JROF,JLEV)
      YDVARS%TNL_SI%T9(JROF,JLEV)=ZBDT(JROF)*PTOD0(JROF,JLEV)
      YDVARS%SVDNL_SI%T9(JROF,JLEV)=ZBT*PLPD0(JROF,JLEV)
    ENDDO
    IF (LNHEE) THEN
      DO JROF = KST,KPROF
        YDVARS%SPDNL_SI%T9(JROF,JLEV)=ZBDT(JROF)*PSPDS0(JROF,JLEV)
      ENDDO
    ENDIF
  ENDDO
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LANHSIB',1,ZHOOK_HANDLE)
END SUBROUTINE LANHSIB

