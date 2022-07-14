SUBROUTINE LANHSIB(YDGEOMETRY,YDGMV, &
 ! --- INPUT ------------------------------------------------------------------
 & YDML_DYN,KST,KPROF,PBT,PBDT,PESGP,&
 & PGAGT0L,PGAGT0M,PTOD0,PSPDS0,PLPD0,PSDIV0,PNHXT0,PMIXNL,&
 ! --- INPUT-OUTPUT -----------------------------------------------------------
 & PGMV,PB1,&
 ! --- OUTPUT -----------------------------------------------------------------
 & PB2)

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
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE YOMGMV             , ONLY : TGMV
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE YOMCT0             , ONLY : LNHEE
USE YOMCT3             , ONLY : NSTEP
USE YOMDYNA            , ONLY : NVDVAR, ND4SYS, LSLINL, LSLINLC2

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
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
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPRB)    :: ZBT,ZBDT(YDGEOMETRY%YRDIM%NPROMA),ZMIXNL
REAL(KIND=JPRB)    :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LANHSIB',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
  & YDMP=>YDGEOMETRY%YRMP, YDDYN=>YDML_DYN%YRDYN,YDPTRSLB1=>YDML_DYN%YRPTRSLB1,YDPTRSLB2=>YDML_DYN%YRPTRSLB2)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, NPROMNH=>YDDIM%NPROMNH, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, XIDT=>YDDYN%XIDT, &
 & NDIMGMV=>YDGMV%NDIMGMV, YT9=>YDGMV%YT9, &
 & MSLB1C9_SI=>YDPTRSLB1%MSLB1C9_SI, MSLB1PD9_SI=>YDPTRSLB1%MSLB1PD9_SI, &
 & MSLB1T9_SI=>YDPTRSLB1%MSLB1T9_SI, MSLB1U9_SI=>YDPTRSLB1%MSLB1U9_SI, &
 & MSLB1V9_SI=>YDPTRSLB1%MSLB1V9_SI, MSLB1VD9_SI=>YDPTRSLB1%MSLB1VD9_SI, &
 & NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & MSLB2PDSI=>YDPTRSLB2%MSLB2PDSI, MSLB2SPSI=>YDPTRSLB2%MSLB2SPSI, &
 & MSLB2TSI=>YDPTRSLB2%MSLB2TSI, MSLB2USI=>YDPTRSLB2%MSLB2USI, &
 & MSLB2VDSI=>YDPTRSLB2%MSLB2VDSI, MSLB2VSI=>YDPTRSLB2%MSLB2VSI, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2)
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
    PB2(JROF,MSLB2USI-1+JLEV)=ZBT*PGAGT0L(JROF,JLEV)
    PB2(JROF,MSLB2VSI-1+JLEV)=ZBT*PGAGT0M(JROF,JLEV)
    PB2(JROF,MSLB2TSI-1+JLEV)=ZBDT(JROF)*PTOD0(JROF,JLEV)
    PB2(JROF,MSLB2VDSI-1+JLEV)=ZBT*PLPD0(JROF,JLEV)
  ENDDO
  IF (LNHEE) THEN
    DO JROF = KST,KPROF
      PB2(JROF,MSLB2PDSI-1+JLEV)=ZBDT(JROF)*PSPDS0(JROF,JLEV)
    ENDDO
  ENDIF
ENDDO
DO JROF = KST,KPROF
  PB2(JROF,MSLB2SPSI)=ZBDT(JROF)*PSDIV0(JROF)
ENDDO

! * first-order decentering factor
IF(PESGP > 1.0_JPRB)THEN
  DO JROF = KST,KPROF
    PB2(JROF,MSLB2SPSI)=PESGP*PB2(JROF,MSLB2SPSI)
  ENDDO
  DO JLEV=1,NFLEVG
    DO JROF = KST,KPROF
      PB2(JROF,MSLB2USI-1+JLEV)=PESGP*PB2(JROF,MSLB2USI-1+JLEV)
      PB2(JROF,MSLB2VSI-1+JLEV)=PESGP*PB2(JROF,MSLB2VSI-1+JLEV)
      PB2(JROF,MSLB2TSI-1+JLEV)=PESGP*PB2(JROF,MSLB2TSI-1+JLEV)
      PB2(JROF,MSLB2VDSI-1+JLEV)=PESGP*PB2(JROF,MSLB2VDSI-1+JLEV)
    ENDDO
    IF (LNHEE) THEN
      DO JROF = KST,KPROF
        PB2(JROF,MSLB2PDSI-1+JLEV)=PESGP*PB2(JROF,MSLB2PDSI-1+JLEV)
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
        PB2(JROF,MSLB2VDSI-1+JLEV)=PB2(JROF,MSLB2VDSI-1+JLEV)&
         & +PESGP*PNHXT0(JROF,JLEV)
      ENDDO
    ENDDO
  ELSEIF (ND4SYS==2 .AND. NCURRENT_ITER>0) THEN
    DO JLEV=1,NFLEVG
      DO JROF = KST,KPROF
        PB2(JROF,MSLB2VDSI-1+JLEV)=PB2(JROF,MSLB2VDSI-1+JLEV)&
         & +PESGP*PGMV(JROF,JLEV,YT9%MNHX)
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
        PB1(JROF,MSLB1C9_SI+JLEV-NFLSA)=PB1(JROF,MSLB1C9_SI+JLEV-NFLSA)&
         & +ZMIXNL*(ZBDT(JROF)*PSDIV0(JROF)-PGMV(JROF,JLEV,YT9%MSPNL_SI))
      ENDDO
    ENDDO
  ENDIF
  IF (LSLINL) THEN
    DO JLEV = 1,NFLEVG
      DO JROF = KST,KPROF
        ZMIXNL = PMIXNL(JROF, JLEV)
        PB1(JROF,MSLB1U9_SI+JLEV-NFLSA)=PB1(JROF,MSLB1U9_SI+JLEV-NFLSA)&
         & + ZMIXNL*(ZBT*PGAGT0L(JROF,JLEV)-PGMV(JROF,JLEV,YT9%MUNL_SI))
        PB1(JROF,MSLB1V9_SI+JLEV-NFLSA)=PB1(JROF,MSLB1V9_SI+JLEV-NFLSA)&
         & + ZMIXNL*(ZBT*PGAGT0M(JROF,JLEV)-PGMV(JROF,JLEV,YT9%MVNL_SI))
        PB1(JROF,MSLB1T9_SI+JLEV-NFLSA)=PB1(JROF,MSLB1T9_SI+JLEV-NFLSA)&
         & + ZMIXNL*(ZBDT(JROF)*PTOD0(JROF,JLEV)-PGMV(JROF,JLEV,YT9%MTNL_SI))
        PB1(JROF,MSLB1VD9_SI+JLEV-NFLSA)=PB1(JROF,MSLB1VD9_SI+JLEV-NFLSA)&
         & + ZMIXNL*(ZBT*PLPD0(JROF,JLEV)-PGMV(JROF,JLEV,YT9%MSVDNL_SI))
      ENDDO
      IF (LNHEE) THEN
        DO JROF = KST,KPROF
          ZMIXNL = PMIXNL(JROF, JLEV)
          PB1(JROF,MSLB1PD9_SI+JLEV-NFLSA)=PB1(JROF,MSLB1PD9_SI+JLEV-NFLSA)&
           & + ZMIXNL*(ZBDT(JROF)*PSPDS0(JROF,JLEV)-PGMV(JROF,JLEV,YT9%MSPDNL_SI))
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
      PGMV(JROF,JLEV,YT9%MSPNL_SI)=ZBDT(JROF)*PSDIV0(JROF)
    ENDDO
  ENDDO
ENDIF
IF (LSLINL .AND. NCURRENT_ITER==0) THEN
  DO JLEV = 1,NFLEVG
    DO JROF = KST,KPROF
      PGMV(JROF,JLEV,YT9%MUNL_SI)=ZBT*PGAGT0L(JROF,JLEV)
      PGMV(JROF,JLEV,YT9%MVNL_SI)=ZBT*PGAGT0M(JROF,JLEV)
      PGMV(JROF,JLEV,YT9%MTNL_SI)=ZBDT(JROF)*PTOD0(JROF,JLEV)
      PGMV(JROF,JLEV,YT9%MSVDNL_SI)=ZBT*PLPD0(JROF,JLEV)
    ENDDO
    IF (LNHEE) THEN
      DO JROF = KST,KPROF
        PGMV(JROF,JLEV,YT9%MSPDNL_SI)=ZBDT(JROF)*PSPDS0(JROF,JLEV)
      ENDDO
    ENDIF
  ENDDO
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LANHSIB',1,ZHOOK_HANDLE)
END SUBROUTINE LANHSIB

