#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE CPG_DYN(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDCPG_TND, YDCPG_MISC, YDCPG_DYN0, YDCPG_DYN9, &
& YDVARS, YDGMV, YDMODEL, PSLHDA, PSLHDD0, PGFL, KSETTLOFF, PGFLPC, PGMV, PGMVS, PB1,   &
& PB2, PGMVT1, PGMVT1S, PGFLT1, PGMVTNDSI, PTRAJ_SLAG)

! -------------------------------------------------------------
!**** *CPG_DYN* - Grid point calculations: non lagged dynamics.

!     Purpose.
!     --------
!           Grid point calculations non lagged dynamics.
!           + simple physics.

!           Abbreviation "vwv" stands for "vertical wind variable".

!**   Interface.
!     ----------
!        *CALL* *CPG_DYN(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        KST       : first element of work.
!        KEND      : last element of work.
!        PBETADT   : BETADT or 0 according to configuration.
!        PDT       : For a leap-frog scheme (three time level scheme):
!                     'dt' at the first time-step, '2 dt' otherwise.
!                    For a 2TL SL scheme: timestep 'dt'.
!        PSLHDA    : Scaling factor of the deformation in f(d) function
!                    (including the model resolution correction)
!        PSLHDD0   : Threshold for deformation tensor enhancement
!        KIBL      : index into YRCSGEOM/YRCSGEOM instances in YDGEOMETRY
!        POROGL    : zonal component of "grad(surf orography)"
!        POROGM    : meridian component of "grad(surf orography)"
!        PRE0      : hydrostatic pressure "prehyd" at half levels at time t.
!        PRDELP0   : 1/(hydrostatic pressure depth of layers) at t.
!        PUVH0     : horizontal wind at time t at half levels.
!        PHI0      : geopotential height "gz" at t at half levels.
!        PHIF0     : geopotential height "gz" at t at full levels.
!        PRCP0     : contains "cp", "R" and "Kap=R/Cp" at t.
!        PCTY0     : contains vertical velocities, vertical integral of divergence at t.
!        PKENE0    : kinetic energy at t.
!        PLSM      : land-sea mask.
!        PTSOL     : surface temperature.
!        PQSOL     : surface specific humidity.
!        PATND     : adiabatic Lagrangian tendencies.
!        PDBBC     : [D (Gw)_surf / Dt]_adiab .
!        PRDPHI    : HYD: not used.
!                    NHEE: contains pre/(R T prehyd [Delta log(prehyd)]) at t.
!                    NHQE: contains 1/(R Tt [Delta log(prehyd)]) at t.
!                    "R" is the version of R (may be Rdry or Rmoist) used in the definition of vertical divergence "dver".
!        PGWT0     : [Gw] at t (LGWADV=T only; levels: see in CPG_GP).
!        PGWT9     : [Gw] at t-dt (LGWADV=T only; levels: see in CPG_GP).
!        PNHXT0    : "X" at full levels at t, diagnosed in CPG_GP.
!        PNHXT9    : "X" at full levels at t-dt, diagnosed in CPG_GP.
!        PNHYT0    : "Y" at half levels at t, diagnosed in CPG_GP.
!        PNHYT9    : "Y" at half levels at t-dt, diagnosed in CPG_GP.
!        PGFL      : unified_treatment grid-point fields at t.

!     INPUT/OUTPUT:
!     -------------
!        PGFLPC    : unified_treatment grid-point fields at t (3TL PC only)
!        PGWS      : [Gw]_surf at t (LRDBBC only).
!        PGMV      : GMV variables at t-dt and t.
!        PGMVS     : GMVS variables at t-dt and t.
!        PB1       : "SLBUF1" buffer.
!        PB2       : "SLBUF2" buffer.
!        PGMVT1    : GMV variables at t+dt.
!        PGMVT1S   : GMVS variables at t+dt.
!        PGFLT1    : t+dt term for the unified_treatment of grid-point fields
!                    (as input contains the non lagged physics).
!        PGMVTNDSI : GMV: tendency due to linear terms (for DDH).
!        PTRAJ_SLAG: Stored SL quantities for TL/AD models

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE documentation about Eulerian and semi-Lagrangian schemes.

! Author
! ------
!   09-Aug-2001 K. YESSAD after part 7 of CPG. 

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   09-Sep-2008 J. Masek  Dataflow for flow deformation along pressure levels
!   K. Yessad (Dec 2008): remove dummy CDLOCK + cleanings
!   F. Vana 13-Jan-2009: Storing KAPPA for TL/AD of SLHD
!   K. Yessad (Nov 2009): cleanings, DT/Dt now pre-computed in CPG_GP.
!   K. Yessad (Nov 2009): prune lpc_old.
!   K. Yessad (Dec 2009): deep-layer NH dynamics for LGWADV=F.
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   E. Bazile  19-01-2011 : Modified CP_FORCING input.
!   K. Yessad (Nov 2011): various contributions, call CP_FORCING moved in CPG_GP.
!   F. Vana    19-03-2012 : supression of storing not available values
!   K. Yessad (Nov 2012): simplify testings.
!   F. Vana  28-Nov-2013 : Redesigned trajectory handling
!   F. Vana    13-Feb-2014 Specific buffer for KAPPAT.
!   M. Diamantakis (Dec 2013): code for LSETTLSVF=true
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (July 2014): Move some variables.
!   F. Vana    21-Nov-2017: Option LSLDP_CURV
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   F. Vana July 2018: RK4 scheme for trajectory research.
! End Modifications
!------------------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE CPG_TYPE_MOD,ONLY : CPG_DYN_TYPE, CPG_MISC_TYPE, CPG_TND_TYPE
USE CPG_OPTS_TYPE_MOD, ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD,ONLY : MF_PHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD,ONLY : FIELD_VARIABLES
USE YOMGMV   , ONLY : TGMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCT0   , ONLY : LSLAG
USE YOMDYNA  , ONLY : LSLHD, LSLHD_STATIC
USE YOMTRAJ  , ONLY : TRAJ_SLAG_TYPE,  LPRTTRAJ, LTRAJSAVE, LTRAJSLAG

!    ------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)       ,INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_BNDS_TYPE)   ,INTENT(IN)   :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE)   ,INTENT(IN)   :: YDCPG_OPTS
TYPE(CPG_TND_TYPE)   ,INTENT(INOUT) :: YDCPG_TND
TYPE(CPG_MISC_TYPE)  ,INTENT(INOUT) :: YDCPG_MISC
TYPE(CPG_DYN_TYPE)   ,INTENT(INOUT) :: YDCPG_DYN0
TYPE(CPG_DYN_TYPE)   ,INTENT(INOUT) :: YDCPG_DYN9
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(TGMV)           ,INTENT(INOUT) :: YDGMV
TYPE(MODEL)          ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM)   ,INTENT(OUT)   :: KSETTLOFF(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PSLHDA(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)      ,INTENT(IN)    :: PSLHDD0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGFLPC(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIMPC)
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV) 
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS) 
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PB1(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1) 
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2) 
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGMVT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%YT1%NDIM) 
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGMVT1S(YDGEOMETRY%YRDIM%NPROMA,YDGMV%YT1%NDIMS) 
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGFLT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1) 
REAL(KIND=JPRB)      ,INTENT(INOUT) :: PGMVTNDSI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
 & YDMODEL%YRML_DIAG%YRMDDH%NDIMSIGMV)
TYPE (TRAJ_SLAG_TYPE),INTENT(INOUT) :: PTRAJ_SLAG
!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "cpg_dyn_slg.intfb.h"
#include "cpg_dyn_eul.intfb.h"
#include "cpg_dyn_tra.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPG_DYN', 0, ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDVAB=>YDGEOMETRY%YRVAB, YDDYN=>YDMODEL%YRML_DYN%YRDYN, &
& YDPTRSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1, YDPTRSLB2=>YDMODEL%YRML_DYN%YRPTRSLB2, YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH,     &
& YDRIP=>YDMODEL%YRML_GCONF%YRRIP, YGFL=>YDMODEL%YRML_GCONF%YGFL, YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,                 &
& YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, YDPHLC=>YDMODEL%YRML_PHY_SLIN%YRPHLC, YDVDF=>YDMODEL%YRML_PHY_G%YRVDF             &
& )

ASSOCIATE(  NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, NPROMA=>YDDIM%NPROMA, NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG,         &
& NFLSA=>YDDIMV%NFLSA, L2TLFF=>YDDYN%L2TLFF, LADVF=>YDDYN%LADVF, LSETTLSVF=>YDDYN%LSETTLSVF, LSLHDHEAT=>YDDYN%LSLHDHEAT,        &
& LSLDP_CURV=>YDDYN%LSLDP_CURV, LSLDP_RK=>YDDYN%LSLDP_RK, LEPHYS=>YDEPHY%LEPHYS, LSPHLC=>YDPHLC%LSPHLC,                         &
& MSLB1C9=>YDPTRSLB1%MSLB1C9, MSLB1GFL9=>YDPTRSLB1%MSLB1GFL9, MSLB1SP9=>YDPTRSLB1%MSLB1SP9, MSLB1T0=>YDPTRSLB1%MSLB1T0,         &
& MSLB1T9=>YDPTRSLB1%MSLB1T9, MSLB1U0=>YDPTRSLB1%MSLB1U0, MSLB1U9=>YDPTRSLB1%MSLB1U9, MSLB1UR0=>YDPTRSLB1%MSLB1UR0,             &
& MSLB1V0=>YDPTRSLB1%MSLB1V0, MSLB1V9=>YDPTRSLB1%MSLB1V9, MSLB1VR0=>YDPTRSLB1%MSLB1VR0, MSLB1WR0=>YDPTRSLB1%MSLB1WR0,           &
& MSLB1ZR0=>YDPTRSLB1%MSLB1ZR0, MSLB1UR00=>YDPTRSLB1%MSLB1UR00, MSLB1VR00=>YDPTRSLB1%MSLB1VR00, MSLB1ZR00=>YDPTRSLB1%MSLB1ZR00, &
& MSLB1WR00=>YDPTRSLB1%MSLB1WR00, MSLB2KAPPA=>YDPTRSLB2%MSLB2KAPPA, MSLB2KAPPAT=>YDPTRSLB2%MSLB2KAPPAT,                         &
& MSLB2URL=>YDPTRSLB2%MSLB2URL, MSLB2VRL=>YDPTRSLB2%MSLB2VRL, MSLB2WRL=>YDPTRSLB2%MSLB2WRL, TSTEP=>YDRIP%TSTEP,                 &
& LMPHYS=>YDPHY%LMPHYS)
!     ------------------------------------------------------------------

!*       2.    Semi-lagrangian scheme
!              ----------------------

!* Cases nsiter=0, predictor and corrector for lpc_full.

IF (LSLAG) THEN
  CALL CPG_DYN_SLG (YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDCPG_TND, YDCPG_MISC, YDCPG_DYN0, YDCPG_DYN9, &
                  & YDVARS, YDGMV, YDMODEL, PSLHDA, PSLHDD0, PGFL, KSETTLOFF, PGMV, PGMVS, PB1, PB2, &
                  & PGMVT1, PGMVT1S, PGMVTNDSI)
ENDIF

!     ------------------------------------------------------------------

!*       3.    Eulerian scheme.
!              ----------------

!* Cases nsiter=0, predictor and corrector for lpc_full.

IF (.NOT.LSLAG) THEN
  CALL CPG_DYN_EUL (YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDCPG_TND, YDCPG_MISC, YDCPG_DYN0, YDCPG_DYN9, &
                  & YDVARS, YDGMV, YDMODEL, PGFL, PGFLPC, PGMV, PGMVS, PB2, PGMVT1, PGMVT1S, PGFLT1)
ENDIF


!     ------------------------------------------------------------------

!*       4.   Stores some variables for AD/TL.
!             --------------------------------

IF (LTRAJSAVE .AND. LTRAJSLAG) THEN
  CALL CPG_DYN_TRA (YDGEOMETRY, YDCPG_BNDS, YDCPG_DYN9, YDVARS, YDMODEL, PB1, PB2, PTRAJ_SLAG)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG_DYN', 1, ZHOOK_HANDLE)
END SUBROUTINE CPG_DYN
