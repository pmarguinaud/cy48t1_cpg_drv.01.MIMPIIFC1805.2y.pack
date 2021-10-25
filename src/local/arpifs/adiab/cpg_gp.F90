#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE CPG_GP(YDGEOMETRY,YDCPG_DYN0,YDCPG_DYN9,YDMF_PHYS_SURF,YDVARS,YDMODEL,YDFIELDS,&
 !---------------------------------------------------------------------
 ! - INPUT .
 & LD_DFISTEP,KST,KEND,KSTC,KENDC,KBL,KSTGLO,LDLFSTEP,LDLDIAB,PDT,PTE,&
 & KIBL,&
 !---------------------------------------------------------------------
 ! - INPUT/OUTPUT .
 & PGFL,&
 & PGMVTNDSL,PGFLTNDSL,&
 !---------------------------------------------------------------------
 ! - OUTPUT .
 & PGPAR,PB2,&
 & PGFLT1,&
 & PKOZO,PQS,PQICE,PQLI,PQRAIN,PQSNOW,&
 & PATND,&
 & PEXTRA,YDDDH)

!**** *CPG_GP* - Grid point calculations:
!                initial part of not lagged grid-point calculations.

!     Purpose.
!     --------
!           Grid point calculations:
!           initial part of not lagged grid-point calculations.
!           - get data in buffers.
!           - multiply p-order horizontal derivatives by M**p.
!           - second part of the temporal filter.
!           - grid-point calculations for nudging.
!           - calls some GP... routines to initialise some auxiliary variables.
!           - sets-up and PB2.

!           Abbreviation "vwv" stands for "vertical wind variable".

!**   Interface.
!     ----------
!        *CALL* *CPG_GP(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        LD_DFISTEP   : 'D' -> DFI computations
!        KST       : first element of work.
!        KEND      : last element of work.
!        KSTC,KENDC: the same as KST,KEND but including zone C for ALADIN.
!        KBL       : block number.
!        KSTGLO    : global offset.
!        LDLFSTEP  : .T.: first time-step?
!        LDLDIAB   : .T. if complete physics is activated and predictor step.
!        PDT       : For a leap-frog scheme (three time level scheme):
!                     'dt' at the first time-step, '2 dt' otherwise.
!                    For a 2TL SL scheme: timestep 'dt'.
!        PTE       : 1. or 0. according to different configurations.
!        KIBL      : index into YRGSGEOM/YRCSGEOM types in YDGEOMETRY

!     INPUT/OUTPUT:
!     -------------
!        PGFL      : unified_treatment grid-point fields at t
!        PGMV      : upper air GMV variables at time t and t-dt.
!        PGMVS     : surface GMV variables at time t and t-dt.
!        PGMVTNDSL : GMV(t+dt,F)-GMV(t or t-dt,O) for DDH
!        PGFLTNDSL : GFL(t+dt,F)-GFL(t or t-dt,O) for DDH

!     OUTPUT:
!     -------
!        POROGL    : zonal component of "grad(surf orography)"
!        POROGM    : meridian component of "grad(surf orography)"
!        PRE0      : hydrostatic pressure "prehyd" at half levels at time t.
!        PRE0L     : zonal component of "grad prehyds" at t.
!        PRE0M     : meridian component of "grad prehyds" at t.
!        PRE0F     : hydrostatic pressure "prehyd" at full levels at time t.
!        PNHPRE0F  : "pre" at full levels (time t).
!        PNHPRE0H  : "pre" at half levels (time t).
!        PXYB0     : contains pressure depth, "delta", "alpha" at t.
!        PUVH0     : horizontal wind at time t at half levels.
!        PHI0      : geopotential height "gz" at t at half levels.
!        PHIF0     : geopotential height "gz" at t at full levels.
!        PHI0FL    : zonal component of "grad (gz)" at t at full levels.
!        PHI0FM    : meridian component of "grad (gz)" at t at full levels.
!        PRCP0     : contains "cp", "R" and "Kap=R/Cp" at t.
!        PCTY0     : contains vertical velocities, vertical integral of divergence at t.
!        PGWFT0    : [Gw] at full layers at t.
!        PKENE0    : kinetic energy at t.
!        PRT0L     : zonal component of "grad RT" at full levels.
!        PRT0M     : meridian component of "grad RT" at full levels.
!        PRE9      : hydrostatic pressure "prehyd" at half levels at t-dt.
!        PRE9F     : hydrostatic pressure "prehyd" at full levels at time t-dt.
!        PNHPRE9F  : "pre" at full levels (time t-dt).
!        PNHPRE9H  : "pre" at half levels (time t-dt).
!        PXYB9     : contains pressure depth, "delta", "alpha" at t-dt.
!        PHI9      : geopotential height "gz" at t-dt at half levels.
!        PHIF9     : geopotential height "gz" at t-dt at full levels.
!        PRCP9     : contains "cp", "R" and "Kap=R/Cp" at t-dt.
!        PGWFT9    : [Gw] at full levels at t-dt.
!        PGPAR     : surface fields for AROME.
!        PB2       : "SLB2" buffer.
!        PGMVT1    : upper air GMV variables at t+dt.
!        PGMVT1S   : surface GMV variables at t+dt.
!        PGFLT1    : GFL variables at t+dt.
!        PKOZO     : fields for photochemistery of ozon.
!        PQS       : specific humidity at surface level.
!        PQICE     : specific humidity of solid water for radiation.
!        PQLI      : specific humidity of liquid water for radiation.
!        PQRAIN    : specific humidity of rain for radiation.
!        PQSNOW    : specific humidity of snow for radiation.
!        PATND     : adiabatic Lagrangian tendencies.
!        PDBBC     : [D (Gw)_surf / Dt]_adiab.
!        PRDPHI    : HYD: not used.
!                    NHEE: contains pre/(R T prehyd [Delta log(prehyd)]) at t.
!                    NHQE: contains 1/(R Tt [Delta log(prehyd)]) at t.
!                    "R" is the version of R (may be Rdry or Rmoist) used in the definition of vertical divergence "dver".
!        PGWT0     : [Gw] at t (LGWADV=T only).
!        PGWT9     : [Gw] at t-dt (LGWADV=T only).
!                    for PGWT0 and PGWT9:
!                    * half level 0 to nflevg-1 values if LVFE_GW=F.
!                    * full level 1 to nflevg values if LVFE_GW=T.
!        PGWS      : [Gw]_surf at t (LRDBBC only).
!        PGWFL     : zonal comp grad(Gw) at full level at t.
!        PGWFM     : meridian comp grad(Gw) at full level at t.
!        PNHXT0    : term 'NHX' at t.
!        PNHXT9    : term 'NHX' at t-dt.
!        PNHYT0    : term 'NHY' at t.    (NHY = gW - gw)
!        PNHYT9    : term 'NHY' at t-dt. (NHY = gW - gw)
!        PQCHA0L   : zonal comp grad(log(pre/prehyd)).
!        PQCHA0M   : merid comp grad(log(pre/prehyd)).
!        PEXTRA    : additional quantity for diagnostics.
!        YDDDH     : diagnostic superstructure

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE documentation vol 2 ch 1 and vol 3 ch 6

!     Author.
!     -------
!        K. YESSAD, after parts 1, 2 and 3 of old CPG.
!        Original : 16-08-2001

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Sep 2008): prune enhanced diffusion (lfrein).
!   F. Bouyssel (Nov 2008): Removal of LPROCLD protection
!   K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (March 2009): correct false comments for LRWSDLG=T
!   K. Yessad (Aug 2009): remove LPC_OLD in TL and AD codes.
!   K. Yessad (Nov 2009): cleanings, DT/Dt now pre-computed in CPG_GP.
!   K. Yessad (Nov 2009): prune lpc_old.
!   K. Yessad (Dec 2009): LRWSDLW,LRWSDLR,LRWSDLG=T,T,T in NH model for LGWADV=F.
!   A. Alias  (Mar 2011) No nudging of SST when LSME=.T.
!   K. Yessad (Nov 2011): various contributions.
!   N. Wedi   (Nov 2011): add LGRADSP
!   M. Ahlgrimm  31-Oct-2011 add rain, snow and PEXTRA to DDH output
!   K. Yessad (Nov 2012): simplify testings.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!   K. Yessad (July 2014): Move some variables.
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Introduce NHQE model.
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   R. Roehrig (Sept 2018): add argument to CP_FORCING + add possibility to
!                           impose (time-evolving) surface pressure in MUSC
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
!   F. Voitus (Dec 2019): NVDVAR=5.
!   P. Smolikova (Sep 2020): Remove obsolete calculations in hyd.model.
! End Modifications
!     ------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE CPG_TYPE_MOD   , ONLY : CPG_DYN_TYPE
USE MFPHYS_SURFACE_TYPE_MOD,ONLY : MFPHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE SURFACE_FIELDS_MIX , ONLY : GPPOPER
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE SC2PRG_MOD,ONLY : SC2PRG
USE YOMCT0   , ONLY : LSLAG, LTWOTL, LNHDYN, LSFORC, LNHEE, LNHQE
USE YOMCT3   , ONLY : NSTEP
USE YOMDYNA  , ONLY : NVDVAR, ND4SYS, LPC_FULL, LNHX, LSLHD
USE YOMGPPB  , ONLY : GPARBUF
USE YOMNUD   , ONLY : NFNUDG, LNUDG
USE YOMLSFORC, ONLY : LSPS_FRC
USE YOMSNU   , ONLY : XPNUDG
USE INTDYN_MOD,ONLY : YYTTND, YYTHW0, YYTHW9,&
 & YYTCTY0, YYTXYBDER0, YYTRCP0, YYTRCP9, YYTXYB0, YYTXYB9
USE YEMLBC_INIT,ONLY : LTENC
USE DDH_MIX, ONLY   : TYP_DDH

! Added ----
USE FIELDS_MOD   ,ONLY : FIELDS
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_DYN_TYPE),INTENT(INOUT) :: YDCPG_DYN0
TYPE(CPG_DYN_TYPE),INTENT(INOUT) :: YDCPG_DYN9
TYPE(MFPHYS_SURF_TYPE),INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
!---
TYPE(FIELDS)      ,INTENT(INOUT) :: YDFIELDS
!---
LOGICAL           ,INTENT(IN)    :: LD_DFISTEP
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTC
INTEGER(KIND=JPIM),INTENT(IN)    :: KENDC
INTEGER(KIND=JPIM),INTENT(IN)    :: KBL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTGLO
LOGICAL           ,INTENT(IN)    :: LDLFSTEP
LOGICAL           ,INTENT(IN)    :: LDLDIAB
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTE
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVTNDSL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,2+YDMODEL%YRML_GCONF%YRDIMF%NFTHER)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLTNDSL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGPAR(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_PHY_MF%YRPARAR%NGPAR+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGFLT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKOZO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG*YDMODEL%YRML_PHY_G%YRDPHY%NVCLIS+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQICE(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQLI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQRAIN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQSNOW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PATND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTTND%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEXTRA(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_PHY_G%YRDPHY%NVEXTRDYN)
TYPE(TYP_DDH)     ,INTENT(INOUT) :: YDDDH

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,IFLAG, JGFL
!!! INTEGER(KIND=JPIM) :: INHFIELDS
LOGICAL :: LLSTR, LLGPXX, LLUVH
REAL(KIND=JPRB) :: ZEPS

REAL(KIND=JPRB), POINTER :: ZP1FORC(:,:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

REAL(KIND=JPRB) :: ZATND_Q(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! 1D model adiab Lagr tendency for "q"
REAL(KIND=JPRB) :: ZGM2

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "cp_forcing.intfb.h"
#include "cp_forcing_ps.intfb.h"
#include "etenc.intfb.h"
#include "gpinislb.intfb.h"
#include "gpinozst.intfb.h"
#include "gpmpfc.intfb.h"
#include "gpnspng.intfb.h"
#include "gp_spv.intfb.h"
#include "gptf2.intfb.h"
#include "updsst.intfb.h"
#include "cpg_gp_hyd.intfb.h"
#include "cpg_gp_nhee.intfb.h"
#include "cpg_gp_nhqe.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPG_GP',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDOROG=>YDGEOMETRY%YROROG(KIBL), &
 & YDPTRSLB2=>YDMODEL%YRML_DYN%YRPTRSLB2,YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR,YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH, &
 & YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY,YDDIMF=>YDMODEL%YRML_GCONF%YRDIMF,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
 & YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH,YGFL=>YDMODEL%YRML_GCONF%YGFL, &
 & YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY, &
 & YDDYN=>YDMODEL%YRML_DYN%YRDYN,YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL, &
 & YDPHLC=>YDMODEL%YRML_PHY_SLIN%YRPHLC)
ASSOCIATE(NDIM=>YGFL%NDIM, NDIM1=>YGFL%NDIM1, NUMFLDS=>YGFL%NUMFLDS, &
 & YCOMP=>YGFL%YCOMP, YFORC=>YGFL%YFORC, YQ=>YGFL%YQ, &
 & NPROMA=>YDDIM%NPROMA, &
 & NPROMNH=>YDDIM%NPROMNH, &
 & NFTHER=>YDDIMF%NFTHER, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NTOZ1D=>YDDPHY%NTOZ1D, NTOZ2D=>YDDPHY%NTOZ2D, NTOZ3D=>YDDPHY%NTOZ3D, &
 & NVCLIS=>YDDPHY%NVCLIS, NVEXTRDYN=>YDDPHY%NVEXTRDYN, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
 & LEO3CH=>YDEPHY%LEO3CH, &
 & RSTRET=>YDGEM%RSTRET, &
 & NDIMGMV=>YDFIELDS%YRGMV%NDIMGMV, NDIMGMVS=>YDFIELDS%YRGMV%NDIMGMVS, YT0=>YDFIELDS%YRGMV%YT0, &
 & YT1=>YDFIELDS%YRGMV%YT1, YT9=>YDFIELDS%YRGMV%YT9, &
 & LRSLDDH=>YDLDDH%LRSLDDH, &
 & MSLDDH_NHX=>YDMDDH%MSLDDH_NHX, MSLDDH_PD=>YDMDDH%MSLDDH_PD, &
 & MSLDDH_T=>YDMDDH%MSLDDH_T, MSLDDH_U=>YDMDDH%MSLDDH_U, &
 & MSLDDH_V=>YDMDDH%MSLDDH_V, MSLDDH_VD=>YDMDDH%MSLDDH_VD, &
 & LSPHLC=>YDPHLC%LSPHLC, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2, &
 & YSD_VF=>YDFIELDS%YRSURF%YSD_VF, YSD_VFD=>YDFIELDS%YRSURF%YSD_VFD, YSD_VV=>YDFIELDS%YRSURF%YSD_VV, &
 & YSD_VVD=>YDFIELDS%YRSURF%YSD_VVD, YSP_RR=>YDFIELDS%YRSURF%YSP_RR, YSP_RRD=>YDFIELDS%YRSURF%YSP_RRD, &
 & YSP_SB=>YDFIELDS%YRSURF%YSP_SB, YSP_SBD=>YDFIELDS%YRSURF%YSP_SBD, YSP_SG=>YDFIELDS%YRSURF%YSP_SG, &
 & YSP_SGD=>YDFIELDS%YRSURF%YSP_SGD, &
 & LSIMPH=>YDSIMPHL%LSIMPH, LMSE=>YDARPHY%LMSE, LMPA=>YDARPHY%LMPA, LMPHYS=>YDPHY%LMPHYS, &
 & NGPAR=>YDPARAR%NGPAR)

ZP1FORC     => YDVARS%FORC(1)%T0

!     ------------------------------------------------------------------

!*       0.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

ZEPS=100.0_JPRB*TINY(1.0_JPRB)
LLSTR=(ABS(RSTRET-1.0_JPRB)>ZEPS)

!     ------------------------------------------------------------------

!*       1.    INTERFACE TO GLOBAL ARRAYS/WORK FILES.
!              --------------------------------------


!*     1.1  TENDENCY COUPLING FOR SURFACE PRESSURE.

IF (LTENC.AND..NOT.LTWOTL) THEN
  CALL ABOR1('CPG_GP: ABOR1 CALLED: LTENC.AND..NOT.LTWOTL')
ENDIF
IF (NSTEP > 0.AND.LTENC) THEN
  CALL ETENC(YDGEOMETRY,YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_LBC,YDFIELDS%YRELBC_FIELDS,&
   & LD_DFISTEP,KENDC,KSTGLO,YDVARS%SP%T0,&
   & YDVARS%SP%DL,YDVARS%SP%DM,&
   & YDVARS%SP%T9,YDVARS%SP%DL,YDVARS%SP%DM)
ENDIF

!*     1.3   SPONGE AT THE TOP OF THE MODEL (ACADEMIC SIMULATIONS).

! new sponge for 2D and 3D models (grid-point GFL only).
IF (YDMODEL%YRML_DYN%YRSPNG%LNSPONGE) CALL GPNSPNG(YGFL,YDMODEL%YRML_DYN%YRSPNG,NPROMA,NFLEVG,KST,KEND,PGFL)

!*     1.4   PART 2 OF TIME FILTER.

IF(NCURRENT_ITER == 0) THEN
  CALL GPTF2(YDGEOMETRY,YDFIELDS%YRGMV,&
   ! --- INPUT ---------------------------------------------------------
   & YDMODEL%YRML_GCONF,YDDYN,KST,KEND,LDLFSTEP,&
   ! --- INPUT for P.T0.; INPUT/OUTPUT for P.T9. -----------------------
   & PGFL=PGFL, &
   & P0SP=YDVARS%SP%T0,    P0SPL=YDVARS%SP%DL,    P0SPM=YDVARS%SP%DM,    P9SP=YDVARS%SP%T9,     &
   & P9SPL=YDVARS%SP%DL9,  P9SPM=YDVARS%SP%DM9,   P0DIV=YDVARS%DIV%T0,   P0NHX=YDVARS%NHX%T0,   &
   & P0SPD=YDVARS%SPD%T0,  P0SPDL=YDVARS%SPD%DL,  P0SPDM=YDVARS%SPD%DM,  P0SVD=YDVARS%SVD%T0,   &
   & P0SVDL=YDVARS%SVD%DL, P0SVDM=YDVARS%SVD%DM,  P0T=YDVARS%T%T0,       P0TL=YDVARS%T%DL,      &
   & P0TM=YDVARS%T%DM,     P0U=YDVARS%U%T0,       P0V=YDVARS%V%T0,       P9DIV=YDVARS%DIV%T9,   &
   & P9NHX=YDVARS%NHX%T9,  P9SPD=YDVARS%SPD%T9,   P9SPDL=YDVARS%SPD%DL9, P9SPDM=YDVARS%SPD%DM9, &
   & P9SVD=YDVARS%SVD%T9,  P9SVDL=YDVARS%SVD%DL9, P9SVDM=YDVARS%SVD%DM9, P9T=YDVARS%T%T9,       &
   & P9TL=YDVARS%T%DL9,    P9TM=YDVARS%T%DM9,     P9U=YDVARS%U%T9,       P9V=YDVARS%V%T9)
ENDIF

!*     1.5   MAP FACTOR AND SURFACE PRESSURE VARIABLES.

YDCPG_DYN0%OROGL(KST:KEND)=YDOROG%OROGL(KST:KEND)*YDGSGEOM%GM(KST:KEND)
YDCPG_DYN0%OROGM(KST:KEND)=YDOROG%OROGM(KST:KEND)*YDGSGEOM%GM(KST:KEND)
IF(LNHDYN) THEN
  DO JROF=KST,KEND
    ZGM2 = YDGSGEOM%GM(JROF)**2
    YDCPG_DYN0%OROGLL(JROF)=YDOROG%OROGLL(JROF)*ZGM2
    YDCPG_DYN0%OROGMM(JROF)=YDOROG%OROGMM(JROF)*ZGM2
    YDCPG_DYN0%OROGLM(JROF)=YDOROG%OROGLM(JROF)*ZGM2
  ENDDO
ENDIF

IF(LLSTR) THEN
  IFLAG=0
  CALL GPMPFC(YDFIELDS%YRGMV,YDMODEL%YRML_GCONF,YDDYN,NPROMA,NFLEVG,KST,KEND,IFLAG,YDGSGEOM%GM,PGFL=PGFL,&
            & P0U=YDVARS%U%T0, P0V=YDVARS%V%T0, P0DIV=YDVARS%DIV%T0, P0TL=YDVARS%T%DL, &
            & P0TM=YDVARS%T%DM, P9U=YDVARS%U%T9, P9V=YDVARS%V%T9, P0UL=YDVARS%U%DL, &
            & P0VL=YDVARS%V%DL, P0VOR=YDVARS%VOR%T0, &
            & P0SPDL=YDVARS%SPD%DL, P0SPDM=YDVARS%SPD%DM, &
            & P0SVDL=YDVARS%SVD%DL, P0SVDM=YDVARS%SVD%DM, P0NHXL=YDVARS%NHX%DL, P0NHXM=YDVARS%NHX%DM, &
            & P9DIV=YDVARS%DIV%T9, P9TL=YDVARS%T%DL9, P9TM=YDVARS%T%DM9, P9SPDL=YDVARS%SPD%DL9, &
            & P9SPDM=YDVARS%SPD%DM9, P9SVDL=YDVARS%SVD%DL9, P9SVDM=YDVARS%SVD%DM9, &
            & P0SPL=YDVARS%SP%DL, P0SPM=YDVARS%SP%DM, &
            & P9SPL=YDVARS%SP%DL9, P9SPM=YDVARS%SP%DM9)
ENDIF

CALL GP_SPV(YDGEOMETRY,YDDYN,YDSIMPHL,.FALSE.,KST,KEND,&
 & YDVARS%SP%T0,YDVARS%SP%DL,YDVARS%SP%DM,&
 & YDVARS%SP%T9,YDVARS%SP%DL,YDVARS%SP%DM,&
 & YDCPG_DYN0%PRE,YDCPG_DYN0%PREL,YDCPG_DYN0%PREM,&
 & YDCPG_DYN9%PRE,YDCPG_DYN9%PREL,YDCPG_DYN9%PREM)

IF (LSFORC .AND. LSPS_FRC) THEN
  CALL CP_FORCING_PS(YDGEOMETRY,KST,KEND,YDCPG_DYN0%PRE)
ENDIF

IF (NCURRENT_ITER == 0) THEN

!*     1.6   SURFACE VARIABLES.

  IF (LDLDIAB.OR.LSPHLC.OR.LSIMPH) THEN
    PQS(KST:KEND)=0.0_JPRB
    IF(LDLFSTEP .AND. .NOT. LTWOTL) THEN
      CALL GPPOPER(YDDYN,'SET9TO0',YDFIELDS%YRSURF,KBL=KBL)
    ENDIF
  ENDIF

  ! for AROME
  IF (LMSE.AND.NGPAR/=0) THEN
    PGPAR(1:NPROMA,1:NGPAR) = GPARBUF (:, :, 1+(KSTGLO-1)/YDGEOMETRY%YRDIM%NPROMA)
  ENDIF

!*     1.8   OZON PHYSICO-CHEMICAL PROPERTIES AND SUBGRID SURFACE TEMPERATURE.

  IF (LDLDIAB) THEN
    IF(NTOZ3D > 0.OR.NTOZ2D > 0.OR.NTOZ1D > 0.AND.(.NOT.LEO3CH)) THEN
      CALL GPINOZST(YDGEOMETRY,YDMODEL%YRML_CHEM%YROZO,YDDPHY,KST,KEND,KSTGLO,PKOZO)
    ENDIF
  ENDIF

!*     1.9  NUDGING.

 IF(LNUDG.AND.(.NOT.LMSE)) THEN
    CALL UPDSST(YDMODEL%YRML_AOC%YRMCC,YDMODEL%YRML_PHY_MF%YRPHY1,NPROMA,KST,KEND,NFNUDG,YSP_SBD%NLEVS,KBL,&
     & XPNUDG,&
     & YDMF_PHYS_SURF%GSP_RR%PT_T0,&
     & YDMF_PHYS_SURF%GSP_SB%PT_T0,YDMF_PHYS_SURF%GSP_RR%PW_T0,&
     & YDMF_PHYS_SURF%GSP_SB%PQ_T0,YDMF_PHYS_SURF%GSP_SG%PF_T0,&
     & YDMF_PHYS_SURF%GSP_SG%PA_T0,YDMF_PHYS_SURF%GSP_SG%PR_T0,&
     & YDMF_PHYS_SURF%GSP_RR%PIC_T0,YDMF_PHYS_SURF%GSP_SB%PTL_T0,&
     & YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS_SURF%GSD_VF%PALBF,&
     & YDMF_PHYS_SURF%GSD_VF%PEMISF,YDMF_PHYS_SURF%GSD_VF%PZ0F,&
     & YDMF_PHYS_SURF%GSD_VV%PZ0H,YDMF_PHYS_SURF%GSD_VV%PD2,&
     & YDMF_PHYS_SURF%GSD_VV%PIVEG)
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       2.    INTERFACE TO GLOBAL ARRAYS/WORK FILES.
!              --------------------------------------

!     ------------------------------------------------------------------

!*       3.    INITIALISE AUXILIARY VARIABLES.
!              -------------------------------

IF (LNHQE) THEN
  ! * NHQE model:
  LLGPXX=.NOT.(LSLAG.AND.NVDVAR == 4.AND.ND4SYS==2)
  CALL CPG_GP_NHQE(YDGEOMETRY,YDVARS,YDCPG_DYN0,YDCPG_DYN9,YDMODEL,YDFIELDS%YRGMV,&
   !---------------------------------------------------------------------
   ! - INPUT .
   & KST,KEND,KSTGLO,LLGPXX,LDLDIAB,LMPA,&
   & KIBL,YDOROG%OROG,&
   !---------------------------------------------------------------------
   ! - INPUT/OUTPUT .
   & PGFL,&
   !---------------------------------------------------------------------
   ! - OUTPUT .
   & PATND)
ELSEIF (LNHEE) THEN
  ! * NHEE model:
  LLGPXX=.NOT.(LSLAG.AND.(NVDVAR == 4 .OR. NVDVAR == 5).AND.ND4SYS==2)
  CALL CPG_GP_NHEE(YDGEOMETRY,YDVARS,YDCPG_DYN0,YDCPG_DYN9,YDMODEL,YDFIELDS%YRGMV,&
   !---------------------------------------------------------------------
   ! - INPUT .
   & KST,KEND,KSTGLO,LLGPXX,LDLDIAB,LMPA,&
   & KIBL,YDOROG%OROG,&
   !---------------------------------------------------------------------
   ! - INPUT/OUTPUT .
   & PGFL,&
   !---------------------------------------------------------------------
   ! - OUTPUT .
   & PATND)
ELSE
  ! * Hydrostatic model:
  YDCPG_DYN0%NHX(:,:)=0.0_JPRB
  YDCPG_DYN9%NHX(:,:)=0.0_JPRB
  YDCPG_DYN0%GWFT(:,:)=0.0_JPRB
  YDCPG_DYN9%GWFT(:,:)=0.0_JPRB
  YDCPG_DYN0%GWHT(:,:)=0.0_JPRB
  YDCPG_DYN9%GWHT(:,:)=0.0_JPRB
  LLGPXX=(LMPHYS.OR.LSIMPH) .AND. NCURRENT_ITER == 0
  LLUVH=LLGPXX.OR.LSLHD
  CALL CPG_GP_HYD(YDGEOMETRY,YDVARS,YDCPG_DYN0,YDCPG_DYN9,YDMODEL,YDFIELDS%YRGMV,&
   !---------------------------------------------------------------------
   ! - INPUT .
   & KST,KEND,KSTGLO,LLUVH,LLGPXX,LDLDIAB,LMPA,&
   & KIBL,YDOROG%OROG,&
   !---------------------------------------------------------------------
   ! - INPUT/OUTPUT .
   & PGFL,&
   !---------------------------------------------------------------------
   ! - OUTPUT .
   & PATND)
ENDIF

!     ------------------------------------------------------------------

!*       4.    TIME DIMENSION.
!              ---------------

CALL GPINISLB(&
 & YDGEOMETRY,YGFL,YDDYN,YDPTRSLB2,KST,KEND,PTE,.FALSE.,&
 & PGFL,&
 & YDVARS%U%T9,YDVARS%V%T9,YDVARS%T%T9,&
 & YDVARS%SPD%T9,YDVARS%SVD%T9,YDCPG_DYN9%NHX,&
 & YDVARS%SP%T9,&
 & YDCPG_DYN0%CTY%VVEL,YDCPG_DYN0%PREF,&
 & YDVARS%U%T1,YDVARS%V%T1,YDVARS%T%T1,PGFLT1,&
 & YDVARS%SPD%T1,YDVARS%SVD%T1,YDVARS%NHX%T1,&
 & YDVARS%SP%T1,&
 & PB2,PQICE,PQLI,PQRAIN,PQSNOW,&
 & PGWFT0=YDCPG_DYN0%GWFT,PGDW0=YDCPG_DYN0%GDW,PGWS0=YDCPG_DYN0%GWHT(:,NFLEVG)) 

IF (NCURRENT_ITER == 0) THEN
  IF (LDLDIAB.OR.LSIMPH) THEN
    CALL GPPOPER(YDDYN,'SET1TO9',YDFIELDS%YRSURF,KBL=KBL)
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       5.    TENDENCIES FOR 1D MODEL.
!              ------------------------

! Add large scale forcing for SCUM (Single Column Unified Model)
! - incrementation of PATND for GMV variables.
! - fills ZATND_Q for humidity.

IF (LSFORC) THEN
  IF (NCURRENT_ITER > 0 .AND. LPC_FULL) THEN
    ! this call to CP_FORCING has not been well checked for the corrector step,
    !  but the corrector step cannot work for GFL variables.
    CALL ABOR1('CPG_GP: part 5')
  ENDIF

  IF (YQ%MP /= YQ%MP1) THEN
    CALL ABOR1('CPG_GP: CP_FORCING NEEDS DOUBLE CHECK ON Q POINTER. REK')
  ENDIF
  CALL CP_FORCING(YDGEOMETRY,YDMODEL,KST,KEND,PDT,YDCPG_DYN0%PRE,YDCPG_DYN0%PREF,ZP1FORC,&
   & YDVARS%U%T0,YDVARS%V%T0,&
   & YDVARS%T%T0,YDVARS%Q%T0,&
   ! most likely not correct : & PATND,ZATND_Q,PGFLT1(1,1,YQ%MP),YDDDH)
   & PATND,ZATND_Q,PGFLT1(:,:,YQ%MP1),YDDDH)
ENDIF

!     ------------------------------------------------------------------

!*       6.    STORE T9 ARRAYS FOR NH TRAJECTORY.
!              ----------------------------------

!!! Obsolete way to manage NH trajectory, to recode properly if developpement of some TL/AD NH code.
!!! IF (LNHDYN.AND.LTRAJNH) THEN
!!!   INHFIELDS=NG3NH95*NFLEVG+1
!!!   CALL WRNHTRAJM(YDGEOMETRY,YDMODEL%YRML_DYN%YRTNH,NSTEP,KST,KEND,KSTGLO,INHFIELDS,&
!!!    & ZT0SPD,ZT9SVD,&
!!!    & ZT9DIV,ZT9U,ZT9V,ZT9T,&
!!!    & ZT9SP)
!!! ENDIF

!     ------------------------------------------------------------------

!*       7.    INITIALZE SL DDH TENDENCY ARRAYS.
!              ----------------------------------

IF (LRSLDDH.AND.LSLAG) THEN

  IF (LTWOTL) THEN
    PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_U)=-YDVARS%U%T0(KST:KEND,1:NFLEVG)
    PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_V)=-YDVARS%V%T0(KST:KEND,1:NFLEVG)
    IF (NFTHER > 0) PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_T)=-YDVARS%T%T0(KST:KEND,1:NFLEVG)
    IF (LNHDYN)THEN !!! ?????? ky: LNHDYN or LNHEE there?
      PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_PD)=-YDVARS%SPD%T0(KST:KEND,1:NFLEVG)
    ENDIF
    IF (LNHDYN)THEN
      PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_VD)=-YDVARS%SVD%T0(KST:KEND,1:NFLEVG)
      IF (LNHX) PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_NHX)=-YDCPG_DYN0%NHX(KST:KEND,1:NFLEVG)
    ENDIF
    DO JGFL=1,NUMFLDS
      IF (YCOMP(JGFL)%LADV) THEN
        PGFLTNDSL(KST:KEND,1:NFLEVG,JGFL)=-PGFL(KST:KEND,1:NFLEVG,YCOMP(JGFL)%MP)
      ENDIF
    ENDDO
  ELSE
    PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_U)=-YDVARS%U%T9(KST:KEND,1:NFLEVG)
    PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_V)=-YDVARS%V%T9(KST:KEND,1:NFLEVG)
    IF (NFTHER > 0) PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_T)=-YDVARS%T%T9(KST:KEND,1:NFLEVG)
    IF (LNHDYN)THEN !!! ?????? ky: LNHDYN or LNHEE there?
      PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_PD)=-YDVARS%SPD%T9(KST:KEND,1:NFLEVG)
    ENDIF
    IF (LNHDYN)THEN
      PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_VD)=-YDVARS%SVD%T9(KST:KEND,1:NFLEVG)
      IF (LNHX) PGMVTNDSL(KST:KEND,1:NFLEVG,MSLDDH_NHX)=-YDCPG_DYN9%NHX(KST:KEND,1:NFLEVG)
    ENDIF
    DO JGFL=1,NUMFLDS
      IF (YCOMP(JGFL)%LADV) THEN
        PGFLTNDSL(KST:KEND,1:NFLEVG,JGFL)=-PGFL(KST:KEND,1:NFLEVG,YCOMP(JGFL)%MP9)
      ENDIF
    ENDDO
  ENDIF

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG_GP',1,ZHOOK_HANDLE)
END SUBROUTINE CPG_GP
