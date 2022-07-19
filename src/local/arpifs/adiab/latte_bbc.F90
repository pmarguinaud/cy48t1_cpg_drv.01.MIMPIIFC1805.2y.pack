SUBROUTINE LATTE_BBC(YDGEOMETRY,YDVARS,YDCPG_SL1,YDCPG_SL2, YDCPG_BNDS, YDCPG_OPTS,&
 ! --- INPUT --------------------------------------------------
 & YDML_DYN,PDTS2,PESGP,PESGM,PDBBC,PRDPHI,PMIXNL,&
 ! --- INPUT/OUTPUT -------------------------------------------
 & PGWS)

!------------------------------------------------------------------------------
!**** *LATTE_BBC*   Semi-Lagrangian scheme for NH model.
!                   Computation of the t and t-dt useful quantities at grid-points.
!                   Additional variables required in the NH
!                   vertical divergence variable when (LRDBBC,LGWADV)=(T,F).

!                   Remark: this routine must not be called for LGWADV=T
!                   because the additional quantities to be interpolated
!                   computed in this routine are useless if LGWADV=T.
!                   When LGWADV=T, [Gw]_surf(t+dt) is computed by a diagnostic
!                   relationship in routine LAPINEB, using some provisional
!                   "(U,V)(t+dt)" data available after the interpolations.

! Purpose
! -------

! Interface
! ---------
!   INPUT:
!     KST         - first element of work.
!     KPROF       - depth of work.
!     PDTS2       - 0.5*time step for the first time-integration step of
!                   a leap-frog scheme or all time-integration steps of
!                   a two-time level scheme; time step for the following
!                   time-integration steps of a leap-frog scheme.
!     PESGP       - (1 + uncentering factor).
!     PESGM       - (1 - uncentering factor).
!     PDBBC       - [D (Gw)_surf / Dt]_adiab.
!     PRDPHI      - NHEE: contains pre/(R T prehyd [Delta log(prehyd)]) at t.
!                   NHQE: contains 1/(R Tt [Delta log(prehyd)]) at t.
!                   "R" is the version of R (may be Rdry or Rmoist) used in the definition of vertical divergence "dver".
!     PMIXNL      - extrapolation control variable for mixed NESC/SETTLS scheme.

!   INPUT/OUTPUT:
!     PGWS        - [Gw]_surf at t (LRDBBC only).
!     PGMVS       - GMVS variables at t-dt and t.
!     PGMV        - GMV variables at t-dt and t.

!   OUTPUT:
!     PB1         - "SLBUF1" buffer for interpolations.
!     PB2         - "SLBUF2" buffer.

! Externals
! ---------
!   none
!   Called by LACDYN.

! Method
! ------

! Reference
! ---------
!   Arpege documentation about semi-Lagrangian scheme.

! Author
! ------
!      Aug-2003 P. SMOLIKOVA, K. YESSAD

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Aug 2008): simplify XIDT treatment with PC + cleanings
!   K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!   K. Yessad (Dec 2009): LRWSDLR=T in NH model for LGWADV=F.
!   K. Yessad (Oct 2013): allow NESC without ICI-PC scheme.
!   K. Yessad (July 2014): Move some variables, rename some variables.
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   J. Vivoda (July 2018): mixed NESC/SETTLS scheme.
! End Modifications
!------------------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE FIELD_VARIABLES_MOD,ONLY : FIELD_VARIABLES
USE CPG_TYPE_MOD,ONLY : CPG_SL1_TYPE, CPG_SL2_TYPE
USE CPG_OPTS_TYPE_MOD, ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE


!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(FIELD_VARIABLES) ,INTENT(INOUT)  :: YDVARS
TYPE(CPG_SL1_TYPE)    ,INTENT(INOUT)  :: YDCPG_SL1
TYPE(CPG_SL2_TYPE)    ,INTENT(INOUT)  :: YDCPG_SL2
TYPE(CPG_BNDS_TYPE)   ,INTENT(IN)     :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE)   ,INTENT(IN)     :: YDCPG_OPTS
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
 
 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESGP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESGM 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDBBC(YDGEOMETRY%YRDIM%NPROMNH) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDPHI(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMIXNL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGWS(YDGEOMETRY%YRDIM%NPROMNH) 
!     -------------------------------------------------------
INTEGER(KIND=JPIM) :: JLEV,JROF

LOGICAL :: LLCT,LLCTC
REAL(KIND=JPRB) :: ZRDPHI(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZMIXNL, ZSETTLS, ZNESC

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LATTE_BBC',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDDYN=>YDML_DYN%YRDYN,YDPTRSLB1=>YDML_DYN%YRPTRSLB1,&
& YDPTRSLB2=>YDML_DYN%YRPTRSLB2)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, NCURRENT_ITER=>YDDYN%NCURRENT_ITER)
!     ------------------------------------------------------------------

!*       1.   Preliminary initializations.

LLCT = YDML_DYN%YRDYNA%LPC_FULL .AND. NCURRENT_ITER > 0  ! corrector step
LLCTC= YDML_DYN%YRDYNA%LPC_CHEAP .AND. NCURRENT_ITER > 0

! compute "pre / (R T prehyd delta)" in ZRDPHI.
ZRDPHI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)=PRDPHI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)

!     ------------------------------------------------------------------

!*       2.   2TLSL scheme.

IF (YDCPG_OPTS%LTWOTL) THEN

  ! * SL2TL:

  IF( .NOT.LLCT )THEN

    !############################################
    ! 2.1 Predictor for LPC_FULL,
    !     or case NSITER=0.
    !############################################

    ! * Save PGWS in PGMVS(.,YT9%MGWS):
    IF (YDML_DYN%YRDYNA%LPC_FULL) YDVARS%GWS%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=PGWS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)

    ! * Calculation of "DPHIT1" and "DBBCT1":
    IF ((YDCPG_OPTS%NSTEP <= 0).OR.YDML_DYN%YRDYNA%LNESC.OR.YDML_DYN%YRDYNA%LSETTLS) THEN
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        YDCPG_SL2%DPHI1(JROF)=PESGP*ZRDPHI(JROF,NFLEVG)
        YDCPG_SL2%DBBC1(JROF)=PESGP*(PDTS2*PDBBC(JROF)*ZRDPHI(JROF,NFLEVG))
      ENDDO
    ELSE
      ! * remaining case: lsettls=false, lnesc=false.
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        YDCPG_SL2%DPHI1(JROF)=PESGP*&
         & (1.5_JPRB*ZRDPHI(JROF,NFLEVG)-0.5_JPRB*YDVARS%DPHI%T9(JROF,NFLEVG))
        YDCPG_SL2%DBBC1(JROF)=PESGP*&
         & (1.5_JPRB*(PDTS2*PDBBC(JROF)*ZRDPHI(JROF,NFLEVG))&
         & -0.5_JPRB*PDTS2*YDVARS%DBBC%T9(JROF)*YDVARS%DPHI%T9(JROF,NFLEVG))
      ENDDO
    ENDIF

    ! * Calculation of PB1(.,MSLB1DPHI9) and PB1(.,MSLB1DBBC9):
    IF ((YDCPG_OPTS%NSTEP <= 0).OR.YDML_DYN%YRDYNA%LNESC) THEN
      DO JLEV=1,NFLEVG
        DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          YDCPG_SL1%DPHI9(JROF,JLEV)=PESGM*ZRDPHI(JROF,JLEV)
          YDCPG_SL1%DBBC9(JROF,JLEV)=&
           & PESGM*PDTS2*PDBBC(JROF)*ZRDPHI(JROF,JLEV)
        ENDDO
      ENDDO
    ELSEIF (YDML_DYN%YRDYNA%LSETTLS) THEN
      DO JLEV=1,NFLEVG
        DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          ZMIXNL  = PMIXNL(JROF,JLEV)
          ZNESC   = PESGM*ZRDPHI(JROF,JLEV)
          ZSETTLS = (1.0_JPRB+PESGM)*ZRDPHI(JROF,JLEV)-YDVARS%DPHI%T9(JROF,JLEV)
          YDCPG_SL1%DPHI9(JROF,JLEV)=ZMIXNL*ZSETTLS+(1.0_JPRB-ZMIXNL)*ZNESC

          ZNESC   = PESGM*PDTS2*PDBBC(JROF)*ZRDPHI(JROF,JLEV)
          ZSETTLS = (1.0_JPRB+PESGM)*PDTS2*PDBBC(JROF)*ZRDPHI(JROF,JLEV)&
                   & -PDTS2*YDVARS%DBBC%T9(JROF)*YDVARS%DPHI%T9(JROF,JLEV) 
          YDCPG_SL1%DBBC9(JROF,JLEV)=ZMIXNL*ZSETTLS+(1.0_JPRB-ZMIXNL)*ZNESC 
        ENDDO
      ENDDO
    ELSE
      ! * remaining case: lsettls=false, lnesc=false.
      DO JLEV=1,NFLEVG
        DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          YDCPG_SL1%DPHI9(JROF,JLEV)=PESGM*&
           & (1.5_JPRB*ZRDPHI(JROF,JLEV)-0.5_JPRB*YDVARS%DPHI%T9(JROF,JLEV))
          YDCPG_SL1%DBBC9(JROF,JLEV)=PESGM*&
           & (1.5_JPRB*(PDTS2*PDBBC(JROF)*ZRDPHI(JROF,JLEV))&
           & -0.5_JPRB*PDTS2*YDVARS%DBBC%T9(JROF)*YDVARS%DPHI%T9(JROF,JLEV))
        ENDDO
      ENDDO
    ENDIF

    ! * Calculation of PB1(.,MSLB1GWS9):
    YDCPG_SL1%GWS9=PGWS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)

    ! * Update t-dt arrays for next time step and corrector
    !   (PGMVS(.,YT9%MDBBC) and PGMV(.,.,YT9%MDPHI)).
    YDVARS%DBBC%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=PDBBC(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
    YDVARS%DPHI%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)=ZRDPHI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)

  ELSE

    !############################################
    ! 2.2 Corrector for LPC_FULL
    !############################################

    ! * Save PGMVS(.,YT9%MGWS) in PGWS:
    PGWS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDVARS%GWS%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)

    ! * Calculation of "DPHIT1" and "DBBCT1":
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDCPG_SL2%DPHI1(JROF)=PESGP*ZRDPHI(JROF,NFLEVG)
      YDCPG_SL2%DBBC1(JROF)=PESGP*(PDTS2*PDBBC(JROF)*ZRDPHI(JROF,NFLEVG))
    ENDDO

    ! * Calculation of PB1(.,MSLB1DPHI9) and PB1(.,MSLB1DBBC9):
    IF (.NOT.LLCTC) THEN
      DO JLEV=1,NFLEVG
        DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          YDCPG_SL1%DPHI9(JROF,JLEV)=PESGM*YDVARS%DPHI%T9(JROF,JLEV)
          YDCPG_SL1%DBBC9(JROF,JLEV)=&
           & PESGM*PDTS2*YDVARS%DBBC%T9(JROF)*YDVARS%DPHI%T9(JROF,JLEV)
        ENDDO
      ENDDO
    ENDIF

    ! * Calculation of PB1(.,MSLB1GWS9):
    IF (.NOT.LLCTC) YDCPG_SL1%GWS9=YDVARS%GWS%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)

  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       3.   3TLSL scheme.

IF (.NOT.YDCPG_OPTS%LTWOTL) THEN

  ! * SL3TL:
  ! * Remark: the case LLCT is not yet coded.

  ! * Calculation of "DPHIT1" and "DBBCT1":
  DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    YDCPG_SL2%DPHI1(JROF)=PESGP*ZRDPHI(JROF,NFLEVG)
    YDCPG_SL2%DBBC1(JROF)=PESGP*(PDTS2*PDBBC(JROF)*ZRDPHI(JROF,NFLEVG))
  ENDDO

  ! * Calculation of PB1(.,MSLB1DPHI9) and PB1(.,MSLB1DBBC9):
  IF (.NOT.LLCT) THEN
    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        YDCPG_SL1%DPHI9(JROF,JLEV)=PESGM*ZRDPHI(JROF,JLEV)
        YDCPG_SL1%DBBC9(JROF,JLEV)=&
         & PESGM*(PDTS2*PDBBC(JROF)*ZRDPHI(JROF,JLEV))
      ENDDO
    ENDDO
  ENDIF

  ! * Calculation of PB1(.,MSLB1GWS9):
  IF (.NOT.LLCT) YDCPG_SL1%GWS9=PGWS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)

ENDIF

!     ------------------------------------------------------------------

!*       4.   Final calculations.

! * Fill the jlev=kflev+1 and jlev=0 parts of PB1(.,MSLB1DPHI9)
!   and PB1(.,MSLB1DBBC9) with zeros
!   (as it is done for the other variables of "PB1" in routine LAVABO).
IF (.NOT.LLCTC) THEN
  YDCPG_SL1%DPHI9=0.0_JPRB
  YDCPG_SL1%DPHI9(:,NFLEVG+1)=0.0_JPRB
  YDCPG_SL1%DBBC9=0.0_JPRB
  YDCPG_SL1%DBBC9(:,NFLEVG+1)=0.0_JPRB
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LATTE_BBC',1,ZHOOK_HANDLE)
END SUBROUTINE LATTE_BBC

