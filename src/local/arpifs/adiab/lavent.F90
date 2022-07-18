!OCL  NOEVAL
SUBROUTINE LAVENT(YDCST, YDGEOMETRY,YDVARS,YDCPG_SL1, YDCPG_SL2, YDRIP,YDML_DYN,YDPARAR,KST,KPROF,KSETTLOFF,LD2TLFF1, &
 & PDTS2,PRDELP,PEVEL,PATND,PWRL9)

!**** *LAVENT*   Semi-Lagrangian scheme.
!                Computation of wind components necessary 
!                to find the SL trajectory.

!     Purpose.
!     --------
!          This subroutine computes the wind components at t
!          at each grid-point of the colocation grid (Gauss
!          grid and poles at all levels). These quantities will
!          be used in routine (E)LARMES to find the semi-Lagrangian trajectory.

!**   Interface.
!     ----------
!        *CALL* *LAVENT(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST       - first element of work.
!          KPROF     - depth of work.
!          KIBL      - index into YDGSGEOM instance in YDGEOMETRY (NPROMA chunks counter)
!          LD2TLFF1  - .T./.F.: Refined treatement of (2*Omega Vec r) at
!                      the origin point when there is t-dt (or t in SL2TL)
!                      physics / Other cases.
!          PDTS2     - 0.5*(time step) for the first time-integration step of
!                      a leap-frog scheme or all time-integration steps of
!                      a two-time level scheme; time step for the following
!                      time-integration steps of a leap-frog scheme.
!          PRDELP    - "1/(pressure depth of layers)" at t.
!          PEVEL     - "etadot (d prehyd/d eta)" at t.
!          PATND     - adiabatic Lagrangian tendencies.

!        INPUT/OUTPUT:
!          PGMV      - GMV variables at t-dt and t.
!          PB1       - "SLBUF1" buffer for interpolations.

!        OUTPUT:
!          PB2       - "SLBUF2" buffer.
!          KSETTLOFF - counter for SETTLSTF=T (# of points new scheme activated at each lev)
!          PWRL9     - vertical velocity at t - dt. (to be eventually stored in trajectory)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           none
!           Called by LACDYN.

!     Reference.
!     ----------
!             Arpege documentation about semi-Lagrangian scheme.

!     Author.
!     -------
!      K. YESSAD (METEO FRANCE/CNRM/GMAP) after old loops 3202 and 3212 of LACDYN.
!      Original : AUGUST 1995.

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad Aug 2008: rationalisation of dummy argument interfaces + optimis.
!   F. Vana  15-Oct-2009: NSPLTHOI option
!   K. Yessad Nov 2009: rm case (lsettlst,lpc_nesct)=(f,f) for sl2tl
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   F. Vana  22-Feb-2011: NWLAG=4
!   K. Yessad (Nov 2011): use PATND for RHS of wind equation.
!   K. Yessad (Oct 2013): allow NESCT or SETTLST for predictor of LPC_CHEAP.
!   K. Yessad (Oct 2013): use NESCV or SETTLSV for vertical displacement.
!   K. Yessad (July 2014): Rename some variables, move some variables.
!   M. Diamantakis (Feb 2014): code for LSETTLSVF=T (SETTLST filter)
!   F. Vana and M. Diamantakis (Aug 2016): regularization of LSETTLSVF=T
!   F. Vana and J. Masek   21-Nov-2017: Option LSLDP_CURV
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   F. Vana July 2018: RK4 scheme for trajectory research.
!   H Petithomme (Dec 2020): turn ZUCOMP to scalar
! End Modifications
!------------------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE CPG_TYPE_MOD       , ONLY : CPG_SL1_TYPE, CPG_SL2_TYPE
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE YOMCST             , ONLY : TCST
USE YOMCT0             , ONLY : LAROME, LTWOTL
USE YOMCT3             , ONLY : NSTEP
USE YOMCVER            , ONLY : LVERTFE

USE YOMRIP             , ONLY : TRIP
USE YOMPARAR           , ONLY : TPARAR
USE INTDYN_MOD         , ONLY : YYTTND

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(CPG_SL1_TYPE)    ,INTENT(INOUT)  :: YDCPG_SL1
TYPE(CPG_SL2_TYPE)    ,INTENT(INOUT)  :: YDCPG_SL2
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
TYPE(TPARAR)      ,INTENT(IN)    :: YDPARAR
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KSETTLOFF(YDGEOMETRY%YRDIMV%NFLEVG)
LOGICAL           ,INTENT(IN)    :: LD2TLFF1 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVEL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PATND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTTND%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWRL9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) ::JROF, JLEV
REAL(KIND=JPRB) :: ZWRL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZUSA, ZALPHA,ZUCOMP,ZZU,ZZV
LOGICAL :: LLCT, LLSETTLST, LLNESCT, LLSETTLSV, LLNESCV, LLCOND

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAVENT',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDVETA=>YDGEOMETRY%YRVETA,   YDDYN=>YDML_DYN%YRDYN,&
& YDPTRSLB1=>YDML_DYN%YRPTRSLB1,YDPTRSLB2=>YDML_DYN%YRPTRSLB2)

ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, LSETFSTAT=>YDDYN%LSETFSTAT, LSETTLSVF=>YDDYN%LSETTLSVF, LSVTSM=>YDDYN%LSVTSM,   &
& RSCALE=>YDDYN%RSCALE,RSCALEOFF=>YDDYN%RSCALEOFF, LSLDP_CURV=>YDDYN%LSLDP_CURV,   LSLDP_RK=>YDDYN%LSLDP_RK,     &
& NCURRENT_ITER=>YDDYN%NCURRENT_ITER, NFLEVSF=>YDDYN%NFLEVSF,   NSPLTHOI=>YDDYN%NSPLTHOI, NWLAG=>YDDYN%NWLAG,    &
& NFOST=>YDRIP%NFOST,   LSQUALL=>YDPARAR%LSQUALL, VSQUALL=>YDPARAR%VSQUALL)
!     ------------------------------------------------------------------

ZUSA=1.0_JPRB/YDCST%RA

LLCT = (YDML_DYN%YRDYNA%LPC_FULL .AND. NCURRENT_ITER > 0)

LLSETTLST=YDML_DYN%YRDYNA%LSETTLST.AND.(.NOT.YDML_DYN%YRDYNA%LELTRA)
LLNESCT=YDML_DYN%YRDYNA%LNESCT
LLSETTLSV=YDML_DYN%YRDYNA%LSETTLSV.AND.(.NOT.YDML_DYN%YRDYNA%LELTRA)
LLNESCV=YDML_DYN%YRDYNA%LNESCV
LLCOND = (NSPLTHOI /= 0).AND.(NWLAG<4)

!     * Computation of "etadot" at full levels.

IF(LVERTFE) THEN
  ! PEVEL is at full levels.
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      ZWRL0(JROF,JLEV)=PEVEL(JROF,JLEV)&
       & *(YDVETA%VETAH(JLEV)-YDVETA%VETAH(JLEV-1))*PRDELP(JROF,JLEV)
    ENDDO
  ENDDO
ELSE
  ! PEVEL is at half levels.
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      ZWRL0(JROF,JLEV)=&
       & 0.5_JPRB*(PEVEL(JROF,JLEV)+PEVEL(JROF,JLEV-1))&
       & *(YDVETA%VETAH(JLEV)-YDVETA%VETAH(JLEV-1))*PRDELP(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

!     * Computation of wind components necessary to find the SL trajectory.

IF(LTWOTL.AND.LSVTSM) THEN
  DO JLEV=1,NFLEVG
    YDCPG_SL1%WRA(KST:KPROF,JLEV)=ZWRL0(KST:KPROF,JLEV)
  ENDDO
ENDIF


IF (LTWOTL.AND.YDML_DYN%YRDYNA%LELTRA) THEN

  ! Use the RHS of momentum equation.
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      IF (LLCOND) THEN
        ZZU = YDCPG_SL1%UF9(JROF,JLEV)
        ZZV = YDCPG_SL1%VF9(JROF,JLEV)
      ELSE
        ZZU = YDCPG_SL1%U9(JROF,JLEV)
        ZZV = YDCPG_SL1%V9(JROF,JLEV)
      ENDIF
      YDCPG_SL1%UR0(JROF,JLEV)=&
       & YDVARS%U%T0(JROF,JLEV)+0.5_JPRB*ZZU&
       & +PDTS2*PATND(JROF,JLEV,YYTTND%M_TNDU)
      YDCPG_SL1%VR0(JROF,JLEV)=&
       & YDVARS%V%T0(JROF,JLEV)+0.5_JPRB*ZZV&
       & +PDTS2*PATND(JROF,JLEV,YYTTND%M_TNDV)
    ENDDO
  ENDDO

  DO JLEV=1,NFLEVG
    YDCPG_SL1%WR0(KST:KPROF,JLEV)=ZWRL0(KST:KPROF,JLEV)
    YDCPG_SL2%URL(KST:KPROF,JLEV)=YDCPG_SL1%UR0(KST:KPROF,JLEV)
    YDCPG_SL2%VRL(KST:KPROF,JLEV)=YDCPG_SL1%VR0(KST:KPROF,JLEV)
    YDCPG_SL2%WRL(KST:KPROF,JLEV)=YDCPG_SL1%WR0(KST:KPROF,JLEV)
  ENDDO

  IF (LD2TLFF1) THEN
    DO JLEV=1,NFLEVG
      YDCPG_SL1%UR9(KST:KPROF,JLEV)=YDVARS%U%T0(KST:KPROF,JLEV)
      YDCPG_SL1%VR9(KST:KPROF,JLEV)=YDVARS%V%T0(KST:KPROF,JLEV)
    ENDDO
  ENDIF

ELSEIF (LTWOTL.AND.(.NOT.YDML_DYN%YRDYNA%LELTRA)) THEN

  IF (.NOT.LLCT) THEN

    ! * Cases:
    !   NSITER=0
    !   Predictor step of LPC_FULL

    ! -- Horizontal displacement:
    IF ( NSTEP <= NFOST .OR. LLNESCT ) THEN
      ! * non-extrapolating scheme X(t+dt) = X(t) 
      DO JLEV=1,NFLEVG
        YDCPG_SL1%UR0(KST:KPROF,JLEV)=YDVARS%U%T0(KST:KPROF,JLEV)
        YDCPG_SL1%VR0(KST:KPROF,JLEV)=YDVARS%V%T0(KST:KPROF,JLEV)
      ENDDO
      DO JLEV=1,NFLEVG
        YDCPG_SL2%URL(KST:KPROF,JLEV)=YDCPG_SL1%UR0(KST:KPROF,JLEV)
        YDCPG_SL2%VRL(KST:KPROF,JLEV)=YDCPG_SL1%VR0(KST:KPROF,JLEV)
      ENDDO
    ELSEIF (LLSETTLST) THEN
      ! * SETTLS-extrapolating schemes X(t+dt) = f( X(t),X(t-dt) )
      DO JLEV=1,NFLEVG
        DO JROF=KST,KPROF
          YDCPG_SL1%UR0(JROF,JLEV)=2.0_JPRB*YDVARS%U%T0(JROF,JLEV)&
           & -YDVARS%U%T9(JROF,JLEV)
          YDCPG_SL1%VR0(JROF,JLEV)=2.0_JPRB*YDVARS%V%T0(JROF,JLEV)&
           & -YDVARS%V%T9(JROF,JLEV)
        ENDDO
      ENDDO
      DO JLEV=1,NFLEVG
        YDCPG_SL2%URL(KST:KPROF,JLEV)=YDVARS%U%T0(KST:KPROF,JLEV)
        YDCPG_SL2%VRL(KST:KPROF,JLEV)=YDVARS%V%T0(KST:KPROF,JLEV)
      ENDDO
    ENDIF

    ! -- Vertical displacement:
    IF ( NSTEP <= NFOST .OR. LLNESCV ) THEN
      ! * non-extrapolating scheme X(t+dt) = X(t) 
      DO JLEV=1,NFLEVG
        YDCPG_SL1%WR0(KST:KPROF,JLEV)=ZWRL0(KST:KPROF,JLEV)
        YDCPG_SL2%WRL(KST:KPROF,JLEV)=YDCPG_SL1%WR0(KST:KPROF,JLEV)
      ENDDO
      IF (NSTEP == 0 .AND. LSETTLSVF) PWRL9(KST:KPROF,1:NFLEVG)=0.0_JPRB
    ELSEIF (LLSETTLSV) THEN
      ! * SETTLS-extrapolating schemes X(t+dt) = f( X(t),X(t-dt) )
      IF (LSETTLSVF) THEN
        ! detect noisy gridpoints which should not be extrapolated
        DO JLEV=1,NFLEVSF
          KSETTLOFF(JLEV)=0
          DO JROF=KST,KPROF
          ! this is equal to 1-2-1 average for tn-1, tn, tn+1 etadot values where tn+1
          ! is obtained with 2nd order t-extrapolation 
            ZALPHA=0.5_JPRB*(1._JPRB - &
              & TANH(-RSCALE*(ZWRL0(JROF,JLEV)*YDVARS%EDOT%T9(JROF,JLEV)+RSCALEOFF)))

            YDCPG_SL1%WR0(JROF,JLEV)=ZWRL0(JROF,JLEV) + &
              & ZALPHA*(ZWRL0(JROF,JLEV)-YDVARS%EDOT%T9(JROF,JLEV))
            IF (LSETFSTAT) KSETTLOFF(JLEV)=KSETTLOFF(JLEV)+NINT(1.0_JPRB-ZALPHA)
          ENDDO
        ENDDO
        DO JLEV=NFLEVSF+1,NFLEVG
          DO JROF=KST,KPROF
            YDCPG_SL1%WR0(JROF,JLEV)=2.0_JPRB*ZWRL0(JROF,JLEV)&
             & -YDVARS%EDOT%T9(JROF,JLEV)
          ENDDO
        ENDDO
        PWRL9(KST:KPROF,1:NFLEVG)=YDVARS%EDOT%T9(KST:KPROF,1:NFLEVG)
      ELSE
        DO JLEV=1,NFLEVG
          DO JROF=KST,KPROF
            YDCPG_SL1%WR0(JROF,JLEV)=2.0_JPRB*ZWRL0(JROF,JLEV)&
             & -YDVARS%EDOT%T9(JROF,JLEV)
          ENDDO
        ENDDO
#ifdef WITH_MGRIDS
        IF( YDML_DYN%YR_MGRIDS_ADVECTION%ACTIVE() ) THEN
          PWRL9(KST:KPROF,1:NFLEVG)=PGMV(KST:KPROF,1:NFLEVG,YT9%MEDOT)
        ENDIF
#endif
      ENDIF
      DO JLEV = 1, NFLEVG
        YDCPG_SL2%WRL(KST:KPROF,JLEV)=ZWRL0(KST:KPROF,JLEV)
      ENDDO
    ENDIF
    YDVARS%EDOT%T9(KST:KPROF,1:NFLEVG)=ZWRL0(KST:KPROF,1:NFLEVG)

  ELSE 

    ! * Corrector step of LPC_FULL:

    DO JLEV=1,NFLEVG
      YDCPG_SL1%UR0(KST:KPROF,JLEV)=YDVARS%U%T9(KST:KPROF,JLEV)
      YDCPG_SL1%VR0(KST:KPROF,JLEV)=YDVARS%V%T9(KST:KPROF,JLEV)
    ENDDO
    DO JLEV=1,NFLEVG
      YDCPG_SL1%WR0(KST:KPROF,JLEV)=YDVARS%EDOT%T9(KST:KPROF,JLEV)
    ENDDO

    DO JLEV=1,NFLEVG
      YDCPG_SL2%URL(KST:KPROF,JLEV)=YDVARS%U%T0(KST:KPROF,JLEV)
      YDCPG_SL2%VRL(KST:KPROF,JLEV)=YDVARS%V%T0(KST:KPROF,JLEV)
    ENDDO
    DO JLEV=1,NFLEVG
      YDCPG_SL2%WRL(KST:KPROF,JLEV)=ZWRL0(KST:KPROF,JLEV)
    ENDDO

  ENDIF

  IF (LD2TLFF1) THEN
    DO JLEV=1,NFLEVG
      YDCPG_SL1%UR9(KST:KPROF,JLEV)=YDVARS%U%T0(KST:KPROF,JLEV)
      YDCPG_SL1%VR9(KST:KPROF,JLEV)=YDVARS%V%T0(KST:KPROF,JLEV)
    ENDDO
  ENDIF

ELSE

  ! SL3TL:

  DO JLEV=1,NFLEVG
    YDCPG_SL1%UR0(KST:KPROF,JLEV)=YDVARS%U%T0(KST:KPROF,JLEV)
    YDCPG_SL1%VR0(KST:KPROF,JLEV)=YDVARS%V%T0(KST:KPROF,JLEV)
  ENDDO
  DO JLEV=1,NFLEVG
    YDCPG_SL1%WR0(KST:KPROF,JLEV)=ZWRL0(KST:KPROF,JLEV)
    YDCPG_SL2%URL(KST:KPROF,JLEV)=YDCPG_SL1%UR0(KST:KPROF,JLEV)
    YDCPG_SL2%VRL(KST:KPROF,JLEV)=YDCPG_SL1%VR0(KST:KPROF,JLEV)
    YDCPG_SL2%WRL(KST:KPROF,JLEV)=YDCPG_SL1%WR0(KST:KPROF,JLEV)
  ENDDO

  IF (LD2TLFF1) THEN
    DO JLEV=1,NFLEVG
      YDCPG_SL1%UR9(KST:KPROF,JLEV)=YDVARS%U%T9(KST:KPROF,JLEV)
      YDCPG_SL1%VR9(KST:KPROF,JLEV)=YDVARS%V%T9(KST:KPROF,JLEV)
    ENDDO
  ENDIF

ENDIF

IF (LAROME.AND.LSQUALL) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      YDCPG_SL1%VR0(JROF,JLEV)=YDCPG_SL1%VR0(JROF,JLEV)+VSQUALL
      YDCPG_SL2%VRL(JROF,JLEV)=YDCPG_SL2%VRL(JROF,JLEV)+VSQUALL
    ENDDO
  ENDDO
ENDIF

! Filling extra buffers for wind interpolation with quantity valid at time t
! Makes sense for the SETTLS option. In the other case it duplicates the MSLB1[x]R0 content.
IF (LSLDP_RK) THEN
  DO JLEV=1,NFLEVG
    YDCPG_SL1%UR00(KST:KPROF,JLEV)=YDCPG_SL2%URL(KST:KPROF,JLEV)
    YDCPG_SL1%VR00(KST:KPROF,JLEV)=YDCPG_SL2%VRL(KST:KPROF,JLEV)
    YDCPG_SL1%WR00(KST:KPROF,JLEV)=YDCPG_SL2%WRL(KST:KPROF,JLEV)
  ENDDO
ENDIF

! Transforming horizontal wind components from lat,lon to 3D cartesian x,y,z coordinates
!    to mitigate the pole singularity.
! NOTE: This transformation has to be of exact match with the reverse transformation in LARCINA
IF (LSLDP_CURV) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      ZUCOMP=YDCPG_SL1%UR0(JROF,JLEV)
      YDCPG_SL1%ZR0(JROF,JLEV)=YDCPG_SL1%VR0(JROF,JLEV)*YDVARS%GEOMETRY%GSQM2%T0(JROF)
      YDCPG_SL1%UR0(JROF,JLEV)=-YDCPG_SL1%UR0(JROF,JLEV)*YDVARS%GEOMETRY%GESLO%T0(JROF)&
       & -YDCPG_SL1%VR0(JROF,JLEV)*YDVARS%GEOMETRY%GECLO%T0(JROF)*YDVARS%GEOMETRY%GEMU%T0(JROF)
      YDCPG_SL1%VR0(JROF,JLEV)= ZUCOMP*YDVARS%GEOMETRY%GECLO%T0(JROF) &
       & -YDCPG_SL1%VR0(JROF,JLEV)*YDVARS%GEOMETRY%GESLO%T0(JROF)*YDVARS%GEOMETRY%GEMU%T0(JROF)
    ENDDO
  ENDDO
  IF (LSLDP_RK) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KPROF
        ZUCOMP=YDCPG_SL1%UR00(JROF,JLEV)
        YDCPG_SL1%ZR00(JROF,JLEV)= YDCPG_SL1%VR00(JROF,JLEV)*YDVARS%GEOMETRY%GSQM2%T0(JROF)
        YDCPG_SL1%UR00(JROF,JLEV)=-YDCPG_SL1%UR00(JROF,JLEV)*YDVARS%GEOMETRY%GESLO%T0(JROF)&
         & -YDCPG_SL1%VR00(JROF,JLEV)*YDVARS%GEOMETRY%GECLO%T0(JROF)*YDVARS%GEOMETRY%GEMU%T0(JROF)
        YDCPG_SL1%VR00(JROF,JLEV)= ZUCOMP*YDVARS%GEOMETRY%GECLO%T0(JROF) &
         & -YDCPG_SL1%VR00(JROF,JLEV)*YDVARS%GEOMETRY%GESLO%T0(JROF)*YDVARS%GEOMETRY%GEMU%T0(JROF)
      ENDDO
    ENDDO
  ENDIF
ENDIF
!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LAVENT',1,ZHOOK_HANDLE)
END SUBROUTINE LAVENT

