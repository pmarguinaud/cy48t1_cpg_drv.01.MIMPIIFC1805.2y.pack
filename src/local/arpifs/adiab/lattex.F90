#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE LATTEX(YDCST,YDGEOMETRY,YDVARS,YDCPG_SL1,YDCPG_SL2,YDCPG_BNDS,YDCPG_OPTS,&
 ! --- INPUT --------------------------------------------------
 & YDLDDH,YDMDDH,YDML_GCONF,YDML_DYN,PDTS2,PBT,PBDT,PESGP,PESGM,&
 & POROGL,POROGM,&
 & PGAGT0L,PGAGT0M,PTOD0,PSPDS0,PLPD0,&
 & PGAGT9L,PGAGT9M,PTOD9,PSPDS9,PLPD9,&
 & PRDELP,PEVEL,PATND,&
 & PNHXT0,PGWT0,PNHYT0,PNHXT9,PGWT9,PNHYT9,PGFL,PMIXNL,&
 ! --- INPUT/OUTPUT -------------------------------------------
 & PGMVTNDSI,PB1)

!**** *LATTEX*   Semi-Lagrangian scheme.
!                Computation of the t and t-dt useful quantities at grid-points.
!                 Equations for tri-dimensional variables.
!                Abbreviation "vwv" stands for "vertical wind variable".

!     Purpose.
!     --------
!        * This subroutine computes the equation quantities to be
!          interpolated at each grid-point of the colocation grid
!          (Gauss grid at all levels).
!          The terms considered here are the explicit terms and
!          the explicit part of the semi-implicit scheme (terms
!          previously computed in LASSIE or LANHSI).
!          Equations considered here are equations for tri-dimensional
!          variables: momentum, temperature, NH variables, GFL.
!        * Remark 1: when an alternate averaging is used for linear terms
!          in the 2TL SL scheme, the first timestep is treated differently
!          (first order uncentering), no uncentering is applied to the
!          total term ("cursive A") and the term saved in PGMV(1,1,YT9%M[X]NL)
!          is [ (Delta t/2) ("cursive A" - (1 + xidt) beta "cursive B") ]
!          instead of [ (Delta t/2) ("cursive A" - beta "cursive B") ].
!        * Remark 2: for lsettls=true, uncentering is applied to
!          the 'stable' extrapolation if vesl > 0 to avoid instability
!          in the momentum equation.
!        * Remark 3: for lgwadv=true in the NH models, variable
!          "gw" is advected instead of vertical divergence; it is advected at
!          half levels if lvfe_gw=f, full levels if lvfe_gw=t.
!          That means that PLPD0,PLPD9,PATND(.,.,YYTTND%M_TNDGW),PGWT0,PGWT9 contain:
!          - half level values 0 to nflevg-1 if LVFE_GW=F.
!          - full level values 1 to nflevg if LVFE_GW=T.
!        * Remark 4: for PC schemes:
!          - this routine is called for nsiter=0.
!          - this routine is called for the predictor of lpc_full.
!          - this routine is called for the corrector of lpc_full.

!**   Interface.
!     ----------
!        *CALL* *LATTEX(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST         - first element of work.
!          KPROF       - depth of work.
!          PDTS2       - 0.5*time step for the first time-integration step of
!                        a leap-frog scheme or all time-integration steps of
!                        a two-time level scheme; time step for the following
!                        time-integration steps of a leap-frog scheme.
!          PBT         - PDTS2*BETADT (BETADT is in YOMDYN).
!          PBDT        - PBT if semi-implicit scheme with unreduced
!                        divergence, PBT*(c**2/GM**2) if semi-implicit
!                        scheme with reduced divergence.
!          PESGP       - (1 + uncentering factor).
!          PESGM       - (1 - uncentering factor).
!          KIBL        - index into YDGSGEOM instance in YDGEOMETRY
!          POROGL      - zonal component of "grad(surf orography)"
!          POROGM      - meridian component of "grad(surf orography)"
!          PGAGT0L     - semi-implicit term at time t for U-wind equation.
!          PGAGT0M     - semi-implicit term at time t for V-wind equation.
!          PTOD0       - semi-implicit term at time t for temperature equation.
!          PSPDS0      - semi-implicit term at time t for press dep equation.
!          PLPD0       - semi-implicit term at time t for vert div equation.
!          PGAGT9L     - semi-implicit term at time t-dt for U-wind equation.
!          PGAGT9M     - semi-implicit term at time t-dt for V-wind equation.
!          PTOD9       - semi-implicit term at time t-dt for temperature eqn.
!          PSPDS9      - semi-implicit term at time t-dt for press dep equation.
!          PLPD9       - semi-implicit term at time t-dt for vert div equation.
!          PRDELP      - "1/(pressure depth of layers)" at t.
!          PEVEL       - "etadot (d prehyd/d eta)" at half levels at t.
!          PATND       - adiabatic Lagrangian tendencies.
!          PNHXT0      - "X" at full levels at t, diagnosed in CPG_GP.
!          PGWT0       - "gw" at t (LGWADV=T only; layers: see in CPG_GP).
!          PNHYT0      - "Y" at half levels at t, diagnosed in CPG_GP
!          PNHXT9      - "X" at full levels at t-dt, diagnosed in CPG_GP.
!          PGWT9       - "gw" at t-dt (LGWADV=T only; layers: see in CPG_GP).
!          PNHYT9      - "Y" at half levels at t-dt, diagnosed in CPG_GP
!          PGFL        - unified_treatment grid-point fields
!          PMIXNL      - extrapolation control variable for mixed NESC/SETTLS scheme

!        INPUT/OUTPUT:
!          PGMV        - GMV variables at t-dt and t.
!          PGMVT1      - GMV variables at t+dt.
!          PGMVTNDSI   - GMV: tendency due to linear terms (for DDH).
!          PB1         - "SLBUF1" buffer for interpolations.

!        OUTPUT:
!          PB2         - "SLBUF2" buffer.

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
!      K. YESSAD (METEO FRANCE/CNRM/GMAP) after old part 3.2 of LACDYN. 
!      Original : AUGUST 1995.

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!   P. Bechtold+A. Untch 26-10-2008: add LEGWWMS switch for non-orogr. GWD
!   F. Vana  15-Oct-2009: option NSPLTHOI
!   K. Yessad (Nov 2009): cleanings, DT/Dt now pre-computed in CPG_GP.
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   F. Vana  22-Feb-2011: diff of phys. tendencies and LTDIABLIN attribute
!   K. Yessad (Dec 2011): various contributions.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (Oct 2013): allow NESC without ICI-PC scheme.
!   K. Yessad (July 2014): Rename some variables, move some variables.
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Introduce NHQE model.
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   J. Vivoda (July 2018): mixed NESC/SETTLS scheme.
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
! End Modifications
!------------------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD     , ONLY : MODEL_DYNAMICS_TYPE
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK
USE YOMCST                 , ONLY : TCST
USE YOMSTA                 , ONLY : RTSUR
USE YOMCT0                 , ONLY : LNHDYN, LTWOTL, LNHEE
USE YOMCT3                 , ONLY : NSTEP
USE YOMCVER                , ONLY : LVERTFE

USE INTDYN_MOD             , ONLY : YYTTND
USE YOMMDDH                , ONLY : TMDDH
USE YOMLDDH                , ONLY : TLDDH
USE FIELD_VARIABLES_MOD,ONLY : FIELD_VARIABLES
USE CPG_OPTS_TYPE_MOD, ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE CPG_TYPE_MOD, ONLY : CPG_SL1_TYPE, CPG_SL2_TYPE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCST),        INTENT(IN)    :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(CPG_SL1_TYPE)    ,INTENT(INOUT)  :: YDCPG_SL1
TYPE(CPG_SL2_TYPE)    ,INTENT(INOUT)  :: YDCPG_SL2
TYPE(CPG_BNDS_TYPE)   ,INTENT(IN)   :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE)   ,INTENT(IN)   :: YDCPG_OPTS
TYPE(TLDDH)       ,INTENT(IN)    :: YDLDDH
TYPE(TMDDH)       ,INTENT(IN)    :: YDMDDH
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
 
 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBDT(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESGP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESGM 

REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAGT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAGT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTOD0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDS0(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLPD0(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAGT9L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAGT9M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTOD9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDS9(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLPD9(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVEL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PATND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTTND%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHXT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWT0(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHYT0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHXT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWT9(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHYT9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMIXNL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVTNDSI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMDDH%NDIMSIGMV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB1%NFLDSLB1)
!     ------------------------------------------------------------------
REAL(KIND=JPRB)               :: ZWT0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) , ALLOCATABLE :: ZUSI9(:,:)
REAL(KIND=JPRB) , ALLOCATABLE :: ZVSI9(:,:)
REAL(KIND=JPRB) , ALLOCATABLE :: ZTSI9(:,:)
REAL(KIND=JPRB) , ALLOCATABLE :: ZSPDSI9(:,:)
REAL(KIND=JPRB) , ALLOCATABLE :: ZVWVSI9(:,:)
REAL(KIND=JPRB)               :: ZMOY1U(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)               :: ZMOY1V(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)               :: ZMOY1T(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)               :: ZMOY1SPD(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)               :: ZMOY1VWV(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)               :: ZUSELESS(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
REAL(KIND=JPRB)               :: ZT0NHY(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IPX
INTEGER(KIND=JPIM) :: JLEV, JGFL, JROF

LOGICAL :: LLSETTLSW, LLCT, LLCTC
LOGICAL :: LLTDIABLIN, LLVD5

!     -------------------------------------------------------

REAL(KIND=JPRB) :: ZCMSLP

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "lattex_dnt.intfb.h"
#include "lattex_tnt.intfb.h"
#include "lattex_gfl_2tl.intfb.h"
#include "lattex_gfl_3tl.intfb.h"
#include "lattex_gfl_vspltrans.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LATTEX',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGSGEOM=>YDGEOMETRY%YRGSGEOM(YDCPG_BNDS%KBL),                     &
& YDOROG=>YDGEOMETRY%YROROG(YDCPG_BNDS%KBL), YDDYN=>YDML_DYN%YRDYN,   YDPTRSLB1=>YDML_DYN%YRPTRSLB1,YDRIP=>YDML_GCONF%YRRIP,   &
& YDPTRSLB2=>YDML_DYN%YRPTRSLB2,YGFL=>YDML_GCONF%YGFL)

ASSOCIATE(NPROMA=>YDDIM%NPROMA,   NPROMNH=>YDDIM%NPROMNH,   NFLEVG=>YDDIMV%NFLEVG, LADVF=>YDDYN%LADVF,                    &
& NCURRENT_ITER=>YDDYN%NCURRENT_ITER, NSPDLAG=>YDDYN%NSPDLAG,   NSVDLAG=>YDDYN%NSVDLAG, NTLAG=>YDDYN%NTLAG,               &
& NWLAG=>YDDYN%NWLAG, RCMSLP0=>YDDYN%RCMSLP0, RCORDIF=>YDDYN%RCORDIF,   RCORDIH=>YDDYN%RCORDIH,   NFOST=>YDRIP%NFOST,     &
& LRSIDDH=>YDLDDH%LRSIDDH,   MSIDDH_PD0=>YDMDDH%MSIDDH_PD0, MSIDDH_PD9=>YDMDDH%MSIDDH_PD9,   MSIDDH_T0=>YDMDDH%MSIDDH_T0, &
& MSIDDH_T9=>YDMDDH%MSIDDH_T9,   MSIDDH_U0=>YDMDDH%MSIDDH_U0, MSIDDH_U9=>YDMDDH%MSIDDH_U9,   MSIDDH_V0=>YDMDDH%MSIDDH_V0, &
& MSIDDH_V9=>YDMDDH%MSIDDH_V9,   MSIDDH_VD0=>YDMDDH%MSIDDH_VD0, MSIDDH_VD9=>YDMDDH%MSIDDH_VD9)

!     ------------------------------------------------------------------

!*      1.  PRELIMINARY INITIALISATIONS.
!       --------------------------------

!       1.1  Allocate arrays for 3TL scheme:

IF(.NOT.LTWOTL)THEN
  ALLOCATE( ZUSI9(NPROMA,NFLEVG))
  ALLOCATE( ZVSI9(NPROMA,NFLEVG))
  ALLOCATE( ZTSI9(NPROMA,NFLEVG))
  IF( LNHEE ) ALLOCATE( ZSPDSI9(NPROMNH,NFLEVG))
  IF( LNHDYN ) ALLOCATE( ZVWVSI9(NPROMNH,NFLEVG))
ENDIF

!       1.2  Scalar initialisations:

LLCT =YDML_DYN%YRDYNA%LPC_FULL .AND. NCURRENT_ITER > 0  ! corrector step
LLCTC=YDML_DYN%YRDYNA%LPC_CHEAP .AND. NCURRENT_ITER > 0

LLVD5=(YDML_DYN%YRDYNA%NVDVAR==5)

ZCMSLP=RCMSLP0/(YDCST%RD*RTSUR)

!       1.3  reset to zero ddh arrays and pointers

IF (LRSIDDH) THEN
  PGMVTNDSI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG, MSIDDH_U0)=0.0_JPRB
  PGMVTNDSI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG, MSIDDH_V0)=0.0_JPRB
  PGMVTNDSI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG, MSIDDH_T0)=0.0_JPRB
  DO JLEV=1,NFLEVG
    YDCPG_SL1%U9_SI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0.0_JPRB
    YDCPG_SL1%V9_SI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0.0_JPRB
    YDCPG_SL1%T9_SI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0.0_JPRB
  ENDDO
  IF (LNHDYN) THEN !! ?????? LNHDYN or LNHEE?
    PGMVTNDSI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG, MSIDDH_PD0)=0.0_JPRB
  ENDIF
  IF (LNHDYN) THEN
    PGMVTNDSI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG, MSIDDH_VD0)=0.0_JPRB
  ENDIF
  IF (LNHEE) THEN
    DO JLEV=1,NFLEVG
      YDCPG_SL1%PD9_SI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0.0_JPRB
    ENDDO
  ENDIF
  IF (LNHDYN) THEN
    DO JLEV=1,NFLEVG
      YDCPG_SL1%VD9_SI(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0.0_JPRB
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*      2.  TREATMENT OF GMV VARIABLES.
!       -------------------------------

!*       2.1   Momentum equation.

! * LSETTLS is replaced by LLSETTLSW=.FALSE. for wind-eqn if VESL>0 because
!   stable extrapolation deteriorates scores without improving stability.

IF (PESGP > PESGM) THEN
  LLSETTLSW=.FALSE.
ELSE
  LLSETTLSW=YDML_DYN%YRDYNA%LSETTLS
ENDIF

DO JLEV=1,NFLEVG

  DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    YDCPG_SL2%USI(JROF,JLEV)=PBT*PGAGT0L(JROF,JLEV)
    YDCPG_SL2%VSI(JROF,JLEV)=PBT*PGAGT0M(JROF,JLEV)
  ENDDO

  ! * Add pressure gradient term + Rayleigh friction in the wind equation.
  DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZMOY1U(JROF,JLEV)=PDTS2*PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)
    ZMOY1V(JROF,JLEV)=PDTS2*PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)
  ENDDO

  ! * Add "- 2 Omega vec V" explicit contribution when required.
  IF (.NOT.LADVF) THEN
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZMOY1U(JROF,JLEV)=ZMOY1U(JROF,JLEV)&
       & +PDTS2*YDGSGEOM%RCORI(JROF)*YDVARS%V%T0(JROF,JLEV)
      ZMOY1V(JROF,JLEV)=ZMOY1V(JROF,JLEV)&
       & -PDTS2*YDGSGEOM%RCORI(JROF)*YDVARS%U%T0(JROF,JLEV)
    ENDDO
  ENDIF

ENDDO

IF (LTWOTL) THEN
  
  CALL LATTEX_DNT(YDGEOMETRY,YDLDDH,YDRIP,YDDYN,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,LLSETTLSW,NWLAG,PESGP,PESGM,&
   & YDVARS%U%T0,YDVARS%U%T9,ZMOY1U,PMIXNL,&
   & YDCPG_SL2%USI,YDVARS%UNL%T9,YDVARS%U%T1,&
   & YDCPG_SL1%U0,YDCPG_SL1%U9,YDCPG_SL1%UF9,YDVARS%CUNL%T9,&
   & PGMVTNDSI(:,:, MSIDDH_U0),PGMVTNDSI(:,:, MSIDDH_U9),YDCPG_SL1%U9_SI,&
   & YDCPG_SL1%U9_NL)
  
  CALL LATTEX_DNT(YDGEOMETRY,YDLDDH,YDRIP,YDDYN,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,LLSETTLSW,NWLAG,PESGP,PESGM,&
   & YDVARS%V%T0,YDVARS%V%T9,ZMOY1V,PMIXNL,&
   & YDCPG_SL2%VSI,YDVARS%VNL%T9,YDVARS%V%T1,&
   & YDCPG_SL1%V0,YDCPG_SL1%V9,YDCPG_SL1%VF9,YDVARS%CVNL%T9,&
   & PGMVTNDSI(:,:, MSIDDH_V0),PGMVTNDSI(:,:, MSIDDH_V9),YDCPG_SL1%V9_SI,&
   & YDCPG_SL1%V9_NL)
  
ELSE
  
  DO JLEV=1,NFLEVG
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZUSI9(JROF,JLEV)=PBT*PGAGT9L(JROF,JLEV)
      ZVSI9(JROF,JLEV)=PBT*PGAGT9M(JROF,JLEV)
    ENDDO
  ENDDO
  
  CALL LATTEX_TNT(YDGEOMETRY,YDLDDH,YDDYN,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,NWLAG,PESGP,PESGM,YDVARS%U%T9,ZMOY1U,ZUSI9,&
   & YDCPG_SL2%USI,YDVARS%U%T1,YDCPG_SL1%U0,YDCPG_SL1%U9,&
   & YDCPG_SL1%UF9,PGMVTNDSI(:,:, MSIDDH_U0),YDCPG_SL1%U9_SI)
  CALL LATTEX_TNT(YDGEOMETRY,YDLDDH,YDDYN,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,NWLAG,PESGP,PESGM,YDVARS%V%T9,ZMOY1V,ZVSI9,&
   & YDCPG_SL2%VSI,YDVARS%V%T1,YDCPG_SL1%V0,YDCPG_SL1%V9,&
   & YDCPG_SL1%VF9,PGMVTNDSI(:,:, MSIDDH_V0),YDCPG_SL1%V9_SI)
  
ENDIF
  
IF(LADVF) THEN
  DO JLEV=1,NFLEVG
    !dir$ ivdep
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDVARS%U%T1(JROF,JLEV)=YDVARS%U%T1(JROF,JLEV)-YDGSGEOM%GOMVRL(JROF)
      YDVARS%V%T1(JROF,JLEV)=YDVARS%V%T1(JROF,JLEV)-YDGSGEOM%GOMVRM(JROF)
    ENDDO
  ENDDO
ENDIF

!*       2.2   Temperature equation.

DO JLEV=1,NFLEVG

  DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    YDCPG_SL2%TSI(JROF,JLEV)=PBDT(JROF)*PTOD0(JROF,JLEV)
  ENDDO

  ! * compute ZMOY1T
  IF(LVERTFE) THEN
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZWT0(JROF)=PEVEL(JROF,JLEV)*PRDELP(JROF,JLEV)
    ENDDO
  ELSE
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZWT0(JROF)=0.5_JPRB*(PEVEL(JROF,JLEV)&
       & +PEVEL(JROF,JLEV-1))*PRDELP(JROF,JLEV)
    ENDDO
  ENDIF
  DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZMOY1T(JROF,JLEV)=&
     & PDTS2*ZCMSLP*RCORDIF(JLEV)*YDVARS%U%T0(JROF,JLEV)*POROGL(JROF)&
     & +PDTS2*ZCMSLP*RCORDIF(JLEV)*YDVARS%V%T0(JROF,JLEV)*POROGM(JROF)
  ENDDO
  DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZMOY1T(JROF,JLEV)=ZMOY1T(JROF,JLEV)&
     & +PDTS2*PATND(JROF,JLEV,YYTTND%M_TNDT)&
     & +PDTS2*ZCMSLP*(RCORDIH(JLEV)-RCORDIH(JLEV-1))*ZWT0(JROF)*YDOROG%OROG(JROF)
  ENDDO

ENDDO

IF (LTWOTL) THEN
  
  CALL LATTEX_DNT(YDGEOMETRY,YDLDDH,YDRIP,YDDYN,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,YDML_DYN%YRDYNA%LSETTLS,NTLAG,PESGP,PESGM,&
   & YDVARS%T%T0,YDVARS%T%T9,ZMOY1T,PMIXNL,&
   & YDCPG_SL2%TSI,YDVARS%TNL%T9,YDVARS%T%T1,&
   & YDCPG_SL1%T0,YDCPG_SL1%T9,YDCPG_SL1%TF9,YDVARS%CTNL%T9,&
   & PGMVTNDSI(:,:, MSIDDH_T0),PGMVTNDSI(:,:, MSIDDH_T9),YDCPG_SL1%T9_SI,&
   & YDCPG_SL1%T9_NL)
  
ELSE
  
  DO JLEV=1,NFLEVG
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZTSI9(JROF,JLEV)=PBDT(JROF)*PTOD9(JROF,JLEV)
    ENDDO
  ENDDO
  
  CALL LATTEX_TNT(YDGEOMETRY,YDLDDH,YDDYN,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,NTLAG,PESGP,PESGM,YDVARS%T%T9,ZMOY1T,ZTSI9,&
   & YDCPG_SL2%TSI,YDVARS%T%T1,YDCPG_SL1%T0,YDCPG_SL1%T9,&
   & YDCPG_SL1%TF9,PGMVTNDSI(:,:, MSIDDH_T0),YDCPG_SL1%T9_SI)
  
ENDIF

DO JLEV=1,NFLEVG
  !dir$ ivdep
  DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    YDVARS%T%T1(JROF,JLEV)=YDVARS%T%T1(JROF,JLEV)&
     & -RCORDIF(JLEV)*ZCMSLP*YDOROG%OROG(JROF)
  ENDDO
  IF (.NOT.LLCTC) THEN
    !dir$ ivdep
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDCPG_SL1%T9(JROF,JLEV)=YDCPG_SL1%T9(JROF,JLEV)&
       & +RCORDIF(JLEV)*ZCMSLP*YDOROG%OROG(JROF)
    ENDDO
  ENDIF
ENDDO

!*       2.3   Anhydrostatic variables equations: "pressure departure".

IF (LNHEE) THEN

  DO JLEV=1,NFLEVG

    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDCPG_SL2%PDSI(JROF,JLEV)=PBDT(JROF)*PSPDS0(JROF,JLEV)
      ZMOY1SPD(JROF,JLEV)=PDTS2*PATND(JROF,JLEV,YYTTND%M_TNDPD)
    ENDDO

  ENDDO

  IF (LTWOTL) THEN
  
    CALL LATTEX_DNT(YDGEOMETRY,YDLDDH,YDRIP,YDDYN,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,YDML_DYN%YRDYNA%LSETTLS,NSPDLAG,PESGP,PESGM,&
     & YDVARS%SPD%T0,YDVARS%SPD%T9,ZMOY1SPD,PMIXNL,&
     & YDCPG_SL2%PDSI,YDVARS%SPDNL%T9,YDVARS%SPD%T1,&
     & YDCPG_SL1%PD0,YDCPG_SL1%PD9,ZUSELESS,YDVARS%CSPDNL%T9,&
     & PGMVTNDSI(:,:, MSIDDH_PD0),PGMVTNDSI(:,:, MSIDDH_PD9),YDCPG_SL1%PD9_SI,&
     & YDCPG_SL1%PD9_NL)
  ELSE
  
    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZSPDSI9(JROF,JLEV)=PBDT(JROF)*PSPDS9(JROF,JLEV)
      ENDDO
    ENDDO
  
    CALL LATTEX_TNT(YDGEOMETRY,YDLDDH,YDDYN,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,NSPDLAG,PESGP,PESGM,YDVARS%SPD%T9,&
     & ZMOY1SPD,ZSPDSI9,&
     & YDCPG_SL2%PDSI,YDVARS%SPD%T1,YDCPG_SL1%PD0,YDCPG_SL1%PD9,&
     & ZUSELESS,PGMVTNDSI(:,:, MSIDDH_PD0),YDCPG_SL1%PD9_SI)
  ENDIF

ENDIF

!*       2.4a  Anhydrostatic variables equations: "vertical divergence".

IF (LNHDYN.AND.(.NOT.YDML_DYN%YRDYNA%LGWADV)) THEN

  DO JLEV=1,NFLEVG
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDCPG_SL2%VDSI(JROF,JLEV)=PBT*PLPD0(JROF,JLEV)
      ZMOY1VWV(JROF,JLEV)=PDTS2*PATND(JROF,JLEV,YYTTND%M_TNDVD)
    ENDDO
  ENDDO

  IF (LTWOTL) THEN

    CALL LATTEX_DNT(YDGEOMETRY,YDLDDH,YDRIP,YDDYN,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,YDML_DYN%YRDYNA%LSETTLS,NSVDLAG,PESGP,PESGM,&
     & YDVARS%SVD%T0,YDVARS%SVD%T9,ZMOY1VWV,PMIXNL,&
     & YDCPG_SL2%VDSI,YDVARS%VWVNL%T9,YDVARS%SVD%T1,&
     & YDCPG_SL1%VD0,YDCPG_SL1%VD9,YDCPG_SL1%VDF9,YDVARS%CVWVNL%T9,&
     & PGMVTNDSI(:,:, MSIDDH_VD0),PGMVTNDSI(:,:, MSIDDH_VD9),YDCPG_SL1%VD9_SI,&
     & YDCPG_SL1%VD9_NL)

    IF ((YDML_DYN%YRDYNA%NVDVAR==4 .OR. YDML_DYN%YRDYNA%NVDVAR==5) .AND. YDML_DYN%YRDYNA%ND4SYS==1) THEN

      IF (.NOT.LLCT) THEN
        ! * predictor step for LPC_FULL or normal SI step.
        IF(YDML_DYN%YRDYNA%LNESC.OR.(NSTEP <= NFOST)) THEN
          ! * LNESC=T predictor or NSTEP<=0.
          DO JLEV=1,NFLEVG
            !dir$ ivdep
            DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
              YDVARS%SVD%T1(JROF,JLEV)=YDVARS%SVD%T1(JROF,JLEV)&
               & +PNHXT0(JROF,JLEV)
              YDVARS%NHX%T1(JROF,JLEV)=PNHXT0(JROF,JLEV)
              YDCPG_SL1%VD9(JROF,JLEV)=YDCPG_SL1%VD9(JROF,JLEV)&
               & -PNHXT0(JROF,JLEV)
              YDCPG_SL1%NHX9(JROF,JLEV)=0.0_JPRB
            ENDDO
          ENDDO
        ELSEIF (.NOT.YDML_DYN%YRDYNA%LNESC.AND.(NSTEP > NFOST).AND.YDML_DYN%YRDYNA%LSETTLS) THEN
          ! * LSETTLS=T predictor if NSTEP>0.
          DO JLEV=1,NFLEVG
            !dir$ ivdep
            DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
              YDVARS%SVD%T1(JROF,JLEV)=YDVARS%SVD%T1(JROF,JLEV)&
               & +PNHXT0(JROF,JLEV)
              YDVARS%NHX%T1(JROF,JLEV)=PNHXT0(JROF,JLEV)
              YDCPG_SL1%VD9(JROF,JLEV)=YDCPG_SL1%VD9(JROF,JLEV)&
               & -YDVARS%NHX%T9(JROF,JLEV)
              YDCPG_SL1%NHX9(JROF,JLEV)=PNHXT0(JROF,JLEV)&
               & -YDVARS%NHX%T9(JROF,JLEV)
            ENDDO
          ENDDO
        ELSE
          ! * predictor LNESC=F, NSTEP>0, LSETTLS=F.
          DO JLEV=1,NFLEVG
            !dir$ ivdep
            DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
              YDVARS%SVD%T1(JROF,JLEV)=YDVARS%SVD%T1(JROF,JLEV)&
               & +1.5_JPRB*PNHXT0(JROF,JLEV)-0.5_JPRB*YDVARS%NHX%T9(JROF,JLEV)
              YDVARS%NHX%T1(JROF,JLEV)=&
               & 1.5_JPRB*PNHXT0(JROF,JLEV)-0.5_JPRB*YDVARS%NHX%T9(JROF,JLEV)
              YDCPG_SL1%VD9(JROF,JLEV)=YDCPG_SL1%VD9(JROF,JLEV)&
               & -0.5_JPRB*PNHXT0(JROF,JLEV)-0.5_JPRB*YDVARS%NHX%T9(JROF,JLEV)
              YDCPG_SL1%NHX9(JROF,JLEV)=&
               & 0.5_JPRB*(PNHXT0(JROF,JLEV)-YDVARS%NHX%T9(JROF,JLEV))
            ENDDO
          ENDDO
        ENDIF

        !dir$ ivdep
        YDVARS%NHX%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)=PNHXT0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)
      ELSE
        ! * corrector step for LPC_FULL.
        DO JLEV=1,NFLEVG
          !dir$ ivdep
          DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            YDVARS%SVD%T1(JROF,JLEV)=YDVARS%SVD%T1(JROF,JLEV)&
             & +PNHXT0(JROF,JLEV)
            YDVARS%NHX%T1(JROF,JLEV)=PNHXT0(JROF,JLEV)
          ENDDO
          IF (.NOT.LLCTC) THEN
            DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
              YDCPG_SL1%VD9(JROF,JLEV)=YDCPG_SL1%VD9(JROF,JLEV)&
               & -YDVARS%NHX%T9(JROF,JLEV)
              YDCPG_SL1%NHX9(JROF,JLEV)=0.0_JPRB
            ENDDO
          ENDIF
        ENDDO
      ENDIF

    ELSEIF ((YDML_DYN%YRDYNA%NVDVAR==4 .OR. YDML_DYN%YRDYNA%NVDVAR==5) .AND. YDML_DYN%YRDYNA%ND4SYS==2) THEN

      IF (.NOT.LLCT) THEN
        DO JLEV=1,NFLEVG
          !dir$ ivdep
          DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            YDVARS%SVD%T1(JROF,JLEV)=YDVARS%SVD%T1(JROF,JLEV)&
             & +PNHXT0(JROF,JLEV)
            YDVARS%NHX%T1(JROF,JLEV)=PNHXT0(JROF,JLEV)
            YDCPG_SL1%VD9(JROF,JLEV)=YDCPG_SL1%VD9(JROF,JLEV)&
             & -PNHXT0(JROF,JLEV)
            YDCPG_SL1%NHX9(JROF,JLEV)=0.0_JPRB
          ENDDO
        ENDDO
        YDVARS%NHX%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)=PNHXT0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)
      ELSE
        ! * corrector step for LPC_FULL (in this case PNHXT9=PGMV(.,.,YT9%MNHX))
        DO JLEV=1,NFLEVG
          !dir$ ivdep
          DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            YDVARS%SVD%T1(JROF,JLEV)=YDVARS%SVD%T1(JROF,JLEV)&
             & +YDVARS%NHX%T9(JROF,JLEV)
            YDVARS%NHX%T1(JROF,JLEV)=PNHXT0(JROF,JLEV)
          ENDDO
          IF (.NOT.LLCTC) THEN
            DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
              YDCPG_SL1%VD9(JROF,JLEV)=YDCPG_SL1%VD9(JROF,JLEV)&
               & -YDVARS%NHX%T9(JROF,JLEV)
              YDCPG_SL1%NHX9(JROF,JLEV)=0.0_JPRB
            ENDDO
          ENDIF
        ENDDO
      ENDIF

    ENDIF

  ELSE
   
    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZVWVSI9(JROF,JLEV)=PBT*PLPD9(JROF,JLEV)
      ENDDO
    ENDDO
  
    CALL LATTEX_TNT(YDGEOMETRY,YDLDDH,YDDYN,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,NSVDLAG,PESGP,PESGM,YDVARS%SVD%T9,&
     & ZMOY1VWV,ZVWVSI9,&
     & YDCPG_SL2%VDSI,YDVARS%SVD%T1,YDCPG_SL1%VD0,YDCPG_SL1%VD9,&
     & YDCPG_SL1%VDF9,PGMVTNDSI(:,:, MSIDDH_VD0),YDCPG_SL1%VD9_SI)

    ! LPC_FULL not yet coded, only the predictor is provided.
    IF ((YDML_DYN%YRDYNA%NVDVAR==4 .OR. YDML_DYN%YRDYNA%NVDVAR==5) .AND. YDML_DYN%YRDYNA%ND4SYS==1) THEN
      DO JLEV=1,NFLEVG
        !dir$ ivdep
        DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          YDVARS%SVD%T1(JROF,JLEV)=YDVARS%SVD%T1(JROF,JLEV)&
           & +PNHXT0(JROF,JLEV)
          YDVARS%NHX%T1(JROF,JLEV)=PNHXT0(JROF,JLEV)
          YDCPG_SL1%VD9(JROF,JLEV)=YDCPG_SL1%VD9(JROF,JLEV)&
           & +PNHXT0(JROF,JLEV)-2.0_JPRB*PNHXT9(JROF,JLEV)
          YDCPG_SL1%NHX9(JROF,JLEV)=0.0_JPRB 
        ENDDO
      ENDDO
    ELSEIF ((YDML_DYN%YRDYNA%NVDVAR==4 .OR. YDML_DYN%YRDYNA%NVDVAR==5) .AND. YDML_DYN%YRDYNA%ND4SYS==2) THEN
      DO JLEV=1,NFLEVG
        !dir$ ivdep
        DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          YDVARS%SVD%T1(JROF,JLEV)=YDVARS%SVD%T1(JROF,JLEV)&
           & +PNHXT9(JROF,JLEV)
          YDVARS%NHX%T1(JROF,JLEV)=PNHXT0(JROF,JLEV)
          YDCPG_SL1%VD9(JROF,JLEV)=YDCPG_SL1%VD9(JROF,JLEV)&
           & -PNHXT9(JROF,JLEV)
          YDCPG_SL1%NHX9(JROF,JLEV)=0.0_JPRB 
        ENDDO
      ENDDO
    ENDIF

  ENDIF

ENDIF

!*       2.4b  Anhydrostatic variables equations: "gw".

IF (LNHDYN.AND.YDML_DYN%YRDYNA%LGWADV) THEN

  ! To avoid the addition of an extra global 3D array, PB2(.,MSLB2VDSI)
  ! initially holds the SI term for the gw equation. At the end of grid-point
  ! calculations it will be converted to the SI term for the "dver" equation.
  ! The SI stage then follows exactly as usual for the "dver" equation.

  DO JLEV=1,NFLEVG
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDCPG_SL2%VDSI(JROF,JLEV)=PBT*PLPD0(JROF,JLEV)
      ZMOY1VWV(JROF,JLEV)=PDTS2*PATND(JROF,JLEV,YYTTND%M_TNDGW)
    ENDDO
  ENDDO

  IF (LTWOTL) THEN

    ! store "gw(t)" into buffer for corrector step
    IF( YDML_DYN%YRDYNA%LPC_FULL.AND.NCURRENT_ITER==0)THEN
      YDVARS%GW%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)  = PGWT0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)

      IF (LLVD5) YDVARS%NHY%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)=PNHYT0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,0:NFLEVG-1)
    ENDIF

    IF (YDML_DYN%YRDYNA%NVDVAR == 5) THEN
       ZT0NHY(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)=PNHYT0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,0:NFLEVG-1)

       CALL LATTEX_DNT(YDGEOMETRY,YDLDDH,YDRIP,YDDYN,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,YDML_DYN%YRDYNA%LSETTLS,NSVDLAG,PESGP,PESGM,&
        & PGWT0,YDVARS%GW%T9,ZMOY1VWV,PMIXNL,&
        & YDCPG_SL2%VDSI,YDVARS%VWVNL%T9,YDVARS%SVD%T1,&
        & YDCPG_SL1%VD0,YDCPG_SL1%VD9,YDCPG_SL1%VDF9,YDVARS%CVWVNL%T9,&
        & PGMVTNDSI(1,1,MSIDDH_VD0),PGMVTNDSI(1,1,MSIDDH_VD9),YDCPG_SL1%VD9_SI,&
        & YDCPG_SL1%VD9_NL,LDVD5=LLVD5,PDYT0=ZT0NHY,PDYT9=YDVARS%NHY%T9)
    ELSE
       CALL LATTEX_DNT(YDGEOMETRY,YDLDDH,YDRIP,YDDYN,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,YDML_DYN%YRDYNA%LSETTLS,NSVDLAG,PESGP,PESGM,&
        & PGWT0,YDVARS%GW%T9,ZMOY1VWV,PMIXNL,&
        & YDCPG_SL2%VDSI,YDVARS%VWVNL%T9,YDVARS%SVD%T1,&
        & YDCPG_SL1%VD0,YDCPG_SL1%VD9,YDCPG_SL1%VDF9,YDVARS%CVWVNL%T9,&
        & PGMVTNDSI(1,1,MSIDDH_VD0),PGMVTNDSI(1,1,MSIDDH_VD9),YDCPG_SL1%VD9_SI,&
        & YDCPG_SL1%VD9_NL)
    ENDIF
   
    IF ((YDML_DYN%YRDYNA%NVDVAR==4 .OR. YDML_DYN%YRDYNA%NVDVAR==5) .AND. YDML_DYN%YRDYNA%ND4SYS==2) THEN
      IF (.NOT.LLCT) THEN
        DO JLEV=1,NFLEVG
          YDCPG_SL1%NHX9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=PNHXT0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
        ENDDO

        ! ky: this saving does not seem to be necessary for ND4SYS=1,
        !     but this must be checked.
        YDVARS%NHX%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)=PNHXT0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)
      ELSE
        DO JLEV=1,NFLEVG
          YDCPG_SL1%NHX9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=YDVARS%NHX%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
        ENDDO
      ENDIF
    ENDIF

  ELSE

    ! This piece of code has not been validated
    !  LGWADV with SL3TL => ABOR1 in the setup, missing interpolations.

    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZVWVSI9(JROF,JLEV)=PBT*PLPD9(JROF,JLEV)
      ENDDO
    ENDDO

    CALL LATTEX_TNT(YDGEOMETRY,YDLDDH,YDDYN,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,NSVDLAG,PESGP,PESGM,PGWT9,ZMOY1VWV,ZVWVSI9,&
     & YDCPG_SL2%VDSI,YDVARS%SVD%T1,YDCPG_SL1%VD0,YDCPG_SL1%VD9,&
     & YDCPG_SL1%VDF9,PGMVTNDSI(:,:, MSIDDH_VD0),YDCPG_SL1%VD9_SI)

    ! LPC_FULL not yet coded, only the predictor is provided.
    IF ((YDML_DYN%YRDYNA%NVDVAR==4 .OR. YDML_DYN%YRDYNA%NVDVAR==5) .AND. YDML_DYN%YRDYNA%ND4SYS==2) THEN
      DO JLEV=1,NFLEVG
        YDCPG_SL1%NHX9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=PNHXT9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
      ENDDO
    ENDIF

  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*      3.  TREATMENT OF GFL VARIABLES.
!       --------------------------------

IF (LTWOTL) THEN
  CALL LATTEX_GFL_2TL (YDGEOMETRY, YDCPG_BNDS, YDML_GCONF, YDML_DYN, PGFL, PB1, LLCT, LLCTC, LLTDIABLIN)
ELSE
  CALL LATTEX_GFL_3TL (YDGEOMETRY, YDCPG_BNDS, YDML_GCONF, YDML_DYN, PGFL, PB1, LLTDIABLIN)
ENDIF

! * Transform the fields to be interpolated with cubic spline in the
!   vertical to "B-spline space".

IF (.NOT.LLCTC) THEN
  CALL LATTEX_GFL_VSPLTRANS (YDGEOMETRY, YDCPG_BNDS, YDML_GCONF, YDML_DYN, PB1, LLTDIABLIN)
ENDIF

!     ------------------------------------------------------------------

!*      4.  DEALLOCATIONS.
!       ------------------

IF(.NOT.LTWOTL)THEN
  IF( ALLOCATED( ZUSI9 ) ) DEALLOCATE( ZUSI9 )
  IF( ALLOCATED( ZVSI9 ) ) DEALLOCATE( ZVSI9 )
  IF( ALLOCATED( ZTSI9 ) ) DEALLOCATE( ZTSI9 )
  IF( ALLOCATED( ZVWVSI9 ) ) DEALLOCATE( ZVWVSI9 )
  IF( ALLOCATED( ZSPDSI9 ) ) DEALLOCATE( ZSPDSI9 )
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LATTEX',1,ZHOOK_HANDLE)
END SUBROUTINE LATTEX
