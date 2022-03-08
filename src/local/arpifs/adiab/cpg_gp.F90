#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE CPG_GP(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDTMP, YDCPG_TND, YDCPG_MISC, YDCPG_GPAR,     &
& YDCPG_DYN0, YDCPG_DYN9, YDMF_PHYS_SURF, YDVARS, YDMODEL, YDFIELDS, LD_DFISTEP, LDLFSTEP, LDLDIAB, &
& PDT, PTE, PGFL, PGMVTNDSL, PGFLTNDSL, PB2, PGFLT1, PKOZO, YDDDH)

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
USE CPG_TYPE_MOD   , ONLY : CPG_DYN_TYPE, CPG_MISC_TYPE, &
                          & CPG_TND_TYPE, CPG_GP_TMP_TYPE, &
                          & CPG_GPAR_TYPE
USE CPG_OPTS_TYPE_MOD   , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD,ONLY : MF_PHYS_SURF_TYPE
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
TYPE(CPG_BNDS_TYPE),INTENT(IN)   :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE),INTENT(IN)   :: YDCPG_OPTS
TYPE(CPG_GP_TMP_TYPE),INTENT(INOUT):: YDTMP
TYPE(CPG_TND_TYPE),INTENT(INOUT) :: YDCPG_TND
TYPE(CPG_MISC_TYPE),INTENT(INOUT):: YDCPG_MISC
TYPE(CPG_GPAR_TYPE),INTENT(INOUT):: YDCPG_GPAR
TYPE(CPG_DYN_TYPE),INTENT(INOUT) :: YDCPG_DYN0
TYPE(CPG_DYN_TYPE),INTENT(INOUT) :: YDCPG_DYN9
TYPE(MF_PHYS_SURF_TYPE),INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES),INTENT(INOUT)  :: YDVARS
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
!---
TYPE(FIELDS)      ,INTENT(INOUT) :: YDFIELDS
!---
LOGICAL           ,INTENT(IN)    :: LD_DFISTEP
LOGICAL           ,INTENT(IN)    :: LDLFSTEP
LOGICAL           ,INTENT(IN)    :: LDLDIAB
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTE

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVTNDSL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,2+YDMODEL%YRML_GCONF%YRDIMF%NFTHER)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLTNDSL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NUMFLDS)

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGFLT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKOZO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_PHY_G%YRDPHY%NVCLIS)
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

IF (LHOOK) CALL DR_HOOK('CPG_GP', 0, ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, YDGSGEOM=>YDGEOMETRY%YRGSGEOM(YDCPG_BNDS%KBL), &
& YDOROG=>YDGEOMETRY%YROROG(YDCPG_BNDS%KBL), YDPTRSLB2=>YDMODEL%YRML_DYN%YRPTRSLB2, YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR,             &
& YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH, YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY, YDDIMF=>YDMODEL%YRML_GCONF%YRDIMF,                             &
& YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, YGFL=>YDMODEL%YRML_GCONF%YGFL,                                &
& YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY, YDDYN=>YDMODEL%YRML_DYN%YRDYN,                              &
& YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL, YDPHLC=>YDMODEL%YRML_PHY_SLIN%YRPHLC)
ASSOCIATE(NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, YQ=>YGFL%YQ, NPROMA=>YDDIM%NPROMA, NFTHER=>YDDIMF%NFTHER,                     &
& NFLEVG=>YDDIMV%NFLEVG, NTOZ1D=>YDDPHY%NTOZ1D, NTOZ2D=>YDDPHY%NTOZ2D, NTOZ3D=>YDDPHY%NTOZ3D, NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
& LEO3CH=>YDEPHY%LEO3CH, RSTRET=>YDGEM%RSTRET, LRSLDDH=>YDLDDH%LRSLDDH, MSLDDH_NHX=>YDMDDH%MSLDDH_NHX,                            &
& MSLDDH_PD=>YDMDDH%MSLDDH_PD, MSLDDH_T=>YDMDDH%MSLDDH_T, MSLDDH_U=>YDMDDH%MSLDDH_U, MSLDDH_V=>YDMDDH%MSLDDH_V,                   &
& MSLDDH_VD=>YDMDDH%MSLDDH_VD, LSPHLC=>YDPHLC%LSPHLC, YSP_SBD=>YDFIELDS%YRSURF%YSP_SBD, LSIMPH=>YDSIMPHL%LSIMPH,                  &
& LMSE=>YDARPHY%LMSE, LMPA=>YDARPHY%LMPA, LMPHYS=>YDPHY%LMPHYS, NGPAR=>YDPARAR%NGPAR)

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
  CALL ETENC(YDGEOMETRY, YDMODEL%YRML_GCONF%YRRIP, YDMODEL%YRML_LBC, YDFIELDS%YRELBC_FIELDS,                 &
  & LD_DFISTEP, YDCPG_BNDS%KFDIE, YDCPG_BNDS%KSTGLO, YDVARS%SP%T0, YDVARS%SP%DL, YDVARS%SP%DM, YDVARS%SP%T9, &
  & YDVARS%SP%DL, YDVARS%SP%DM)
ENDIF

!*     1.3   SPONGE AT THE TOP OF THE MODEL (ACADEMIC SIMULATIONS).

! new sponge for 2D and 3D models (grid-point GFL only).
IF (YDMODEL%YRML_DYN%YRSPNG%LNSPONGE) CALL GPNSPNG(YGFL, YDMODEL%YRML_DYN%YRSPNG, NPROMA, NFLEVG, YDCPG_BNDS%KIDIA, &
                                      & YDCPG_BNDS%KFDIA, PGFL)

!*     1.4   PART 2 OF TIME FILTER.

IF(NCURRENT_ITER == 0) THEN
  CALL GPTF2(YDGEOMETRY, YDFIELDS%YRGMV, YDMODEL%YRML_GCONF, YDDYN, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,          &
  & LDLFSTEP, PGFL=PGFL, P0SP=YDVARS%SP%T0, P0SPL=YDVARS%SP%DL, P0SPM=YDVARS%SP%DM, P9SP=YDVARS%SP%T9,           &
  & P9SPL=YDVARS%SP%DL9, P9SPM=YDVARS%SP%DM9, P0DIV=YDVARS%DIV%T0, P0NHX=YDVARS%NHX%T0, P0SPD=YDVARS%SPD%T0,     &
  & P0SPDL=YDVARS%SPD%DL, P0SPDM=YDVARS%SPD%DM, P0SVD=YDVARS%SVD%T0, P0SVDL=YDVARS%SVD%DL, P0SVDM=YDVARS%SVD%DM, &
  & P0T=YDVARS%T%T0, P0TL=YDVARS%T%DL, P0TM=YDVARS%T%DM, P0U=YDVARS%U%T0, P0V=YDVARS%V%T0, P9DIV=YDVARS%DIV%T9,  &
  & P9NHX=YDVARS%NHX%T9, P9SPD=YDVARS%SPD%T9, P9SPDL=YDVARS%SPD%DL9, P9SPDM=YDVARS%SPD%DM9, P9SVD=YDVARS%SVD%T9, &
  & P9SVDL=YDVARS%SVD%DL9, P9SVDM=YDVARS%SVD%DM9, P9T=YDVARS%T%T9, P9TL=YDVARS%T%DL9, P9TM=YDVARS%T%DM9,         &
  & P9U=YDVARS%U%T9, P9V=YDVARS%V%T9)
ENDIF

!*     1.5   MAP FACTOR AND SURFACE PRESSURE VARIABLES.

YDCPG_DYN0%OROGL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDOROG%OROGL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)*YDGSGEOM%GM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
YDCPG_DYN0%OROGM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDOROG%OROGM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)*YDGSGEOM%GM(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)
IF(LNHDYN) THEN
  DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZGM2 = YDGSGEOM%GM(JROF)**2
    YDTMP%T0%OROGLL(JROF)=YDOROG%OROGLL(JROF)*ZGM2
    YDTMP%T0%OROGMM(JROF)=YDOROG%OROGMM(JROF)*ZGM2
    YDTMP%T0%OROGLM(JROF)=YDOROG%OROGLM(JROF)*ZGM2
  ENDDO
ENDIF

IF(LLSTR) THEN
  IFLAG=0
  CALL GPMPFC(YDFIELDS%YRGMV, YDMODEL%YRML_GCONF, YDDYN, NPROMA, NFLEVG, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,      &
  & IFLAG, YDGSGEOM%GM, PGFL=PGFL, P0U=YDVARS%U%T0, P0V=YDVARS%V%T0, P0DIV=YDVARS%DIV%T0, P0TL=YDVARS%T%DL,       &
  & P0TM=YDVARS%T%DM, P9U=YDVARS%U%T9, P9V=YDVARS%V%T9, P0UL=YDVARS%U%DL, P0VL=YDVARS%V%DL, P0VOR=YDVARS%VOR%T0,  &
  & P0SPDL=YDVARS%SPD%DL, P0SPDM=YDVARS%SPD%DM, P0SVDL=YDVARS%SVD%DL, P0SVDM=YDVARS%SVD%DM, P0NHXL=YDVARS%NHX%DL, &
  & P0NHXM=YDVARS%NHX%DM, P9DIV=YDVARS%DIV%T9, P9TL=YDVARS%T%DL9, P9TM=YDVARS%T%DM9, P9SPDL=YDVARS%SPD%DL9,       &
  & P9SPDM=YDVARS%SPD%DM9, P9SVDL=YDVARS%SVD%DL9, P9SVDM=YDVARS%SVD%DM9, P0SPL=YDVARS%SP%DL, P0SPM=YDVARS%SP%DM,  &
  & P9SPL=YDVARS%SP%DL9, P9SPM=YDVARS%SP%DM9)
ENDIF

CALL GP_SPV(YDGEOMETRY, YDDYN, YDSIMPHL, .FALSE., YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDVARS%SP%T0, YDVARS%SP%DL, &
& YDVARS%SP%DM, YDVARS%SP%T9, YDVARS%SP%DL, YDVARS%SP%DM, YDCPG_DYN0%PRE, YDCPG_DYN0%PREL, YDCPG_DYN0%PREM,       &
& YDCPG_DYN9%PRE, YDCPG_DYN9%PREL, YDCPG_DYN9%PREM)

IF (LSFORC .AND. LSPS_FRC) THEN
  CALL CP_FORCING_PS(YDGEOMETRY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_DYN0%PRE)
ENDIF

IF (NCURRENT_ITER == 0) THEN

!*     1.6   SURFACE VARIABLES.

  IF (LDLDIAB.OR.LSPHLC.OR.LSIMPH) THEN
    YDCPG_MISC%QS(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0.0_JPRB
    IF(LDLFSTEP .AND. .NOT. LTWOTL) THEN
      CALL GPPOPER(YDDYN, 'SET9TO0', YDFIELDS%YRSURF, KBL=YDCPG_BNDS%KBL)
    ENDIF
  ENDIF

  ! for AROME
  IF (LMSE.AND.NGPAR/=0) THEN
    YDCPG_GPAR%ZVIEW(1:NPROMA,1:NGPAR) = GPARBUF (:, :, 1+(YDCPG_BNDS%KSTGLO-1)/YDGEOMETRY%YRDIM%NPROMA)
  ENDIF

!*     1.8   OZON PHYSICO-CHEMICAL PROPERTIES AND SUBGRID SURFACE TEMPERATURE.

  IF (LDLDIAB) THEN
    IF(NTOZ3D > 0.OR.NTOZ2D > 0.OR.NTOZ1D > 0.AND.(.NOT.LEO3CH)) THEN
      CALL GPINOZST(YDGEOMETRY, YDMODEL%YRML_CHEM%YROZO, YDDPHY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, &
      & YDCPG_BNDS%KSTGLO, PKOZO)
    ENDIF
  ENDIF

!*     1.9  NUDGING.

 IF(LNUDG.AND.(.NOT.LMSE)) THEN
    CALL UPDSST(YDMODEL%YRML_AOC%YRMCC, YDMODEL%YRML_PHY_MF%YRPHY1, NPROMA, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,            &
    & NFNUDG, YSP_SBD%NLEVS, YDCPG_BNDS%KBL, XPNUDG, YDMF_PHYS_SURF%GSP_RR%PT_T0, YDMF_PHYS_SURF%GSP_SB%PT_T0,             &
    & YDMF_PHYS_SURF%GSP_RR%PW_T0, YDMF_PHYS_SURF%GSP_SB%PQ_T0, YDMF_PHYS_SURF%GSP_SG%PF_T0, YDMF_PHYS_SURF%GSP_SG%PA_T0,  &
    & YDMF_PHYS_SURF%GSP_SG%PR_T0, YDMF_PHYS_SURF%GSP_RR%PIC_T0, YDMF_PHYS_SURF%GSP_SB%PTL_T0, YDMF_PHYS_SURF%GSD_VF%PLSM, &
    & YDMF_PHYS_SURF%GSD_VF%PALBF, YDMF_PHYS_SURF%GSD_VF%PEMISF, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H,   &
    & YDMF_PHYS_SURF%GSD_VV%PD2, YDMF_PHYS_SURF%GSD_VV%PIVEG)
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
  CALL CPG_GP_NHQE(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDTMP, YDVARS, YDCPG_DYN0, YDCPG_DYN9, &
  & YDMODEL, YDFIELDS%YRGMV, LLGPXX, LDLDIAB, LMPA, YDOROG%OROG, PGFL, YDCPG_TND%ZVIEW)
ELSEIF (LNHEE) THEN
  ! * NHEE model:
  LLGPXX=.NOT.(LSLAG.AND.(NVDVAR == 4 .OR. NVDVAR == 5).AND.ND4SYS==2)
  CALL CPG_GP_NHEE(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDTMP, YDVARS, YDCPG_DYN0, YDCPG_DYN9, &
  & YDMODEL, YDFIELDS%YRGMV, LLGPXX, LDLDIAB, LMPA, YDOROG%OROG, PGFL, YDCPG_TND%ZVIEW)
ELSE
  ! * Hydrostatic model:
  YDCPG_DYN0%NHX(:,:)=0.0_JPRB
  YDCPG_DYN0%GWFT(:,:)=0.0_JPRB
  YDCPG_DYN9%NHX(:,:)=0.0_JPRB
  YDCPG_DYN9%GWFT(:,:)=0.0_JPRB
  LLGPXX=(LMPHYS.OR.LSIMPH) .AND. NCURRENT_ITER == 0
  LLUVH=LLGPXX.OR.LSLHD
  CALL CPG_GP_HYD(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDTMP, YDVARS, YDCPG_DYN0, YDCPG_DYN9, YDMODEL, &
  & YDFIELDS%YRGMV, LLUVH, LLGPXX, LDLDIAB, LMPA, YDOROG%OROG, PGFL, YDCPG_TND%ZVIEW)
ENDIF

!     ------------------------------------------------------------------

!*       4.    TIME DIMENSION.
!              ---------------

CALL GPINISLB(  YDGEOMETRY, YGFL, YDDYN, YDPTRSLB2, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, PTE, .FALSE.,      &
& PGFL, YDVARS%U%T9, YDVARS%V%T9, YDVARS%T%T9, YDVARS%SPD%T9, YDVARS%SVD%T9, YDCPG_DYN9%NHX, YDVARS%SP%T9, &
& YDCPG_DYN0%CTY%VVEL, YDCPG_DYN0%PREF, YDVARS%U%T1, YDVARS%V%T1, YDVARS%T%T1, PGFLT1, YDVARS%SPD%T1,      &
& YDVARS%SVD%T1, YDVARS%NHX%T1, YDVARS%SP%T1, PB2, YDCPG_MISC%QICE, YDCPG_MISC%QLI, YDCPG_MISC%QRAIN,      &
& YDCPG_MISC%QSNOW, PGWFT0=YDCPG_DYN0%GWFT, PGDW0=YDTMP%T0%GDW, PGWS0=YDTMP%T0%GWHT(:, NFLEVG)) 

IF (NCURRENT_ITER == 0) THEN
  IF (LDLDIAB.OR.LSIMPH) THEN
    CALL GPPOPER(YDDYN, 'SET1TO9', YDFIELDS%YRSURF, KBL=YDCPG_BNDS%KBL)
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
  CALL CP_FORCING(YDGEOMETRY, YDMODEL, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, PDT, YDCPG_DYN0%PRE, YDCPG_DYN0%PREF, &
  & ZP1FORC, YDVARS%U%T0, YDVARS%V%T0, YDVARS%T%T0, YDVARS%Q%T0, YDCPG_TND%ZVIEW, ZATND_Q, PGFLT1(:, :, YQ%MP1), &
  & YDDDH)
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
    PGMVTNDSL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,MSLDDH_U)=-YDVARS%U%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)
    PGMVTNDSL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,MSLDDH_V)=-YDVARS%V%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)
    IF (NFTHER > 0) PGMVTNDSL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,MSLDDH_T)=-YDVARS%T%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)
    IF (LNHDYN)THEN !!! ?????? ky: LNHDYN or LNHEE there?
      PGMVTNDSL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,MSLDDH_PD)=-YDVARS%SPD%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)
    ENDIF
    IF (LNHDYN)THEN
      PGMVTNDSL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,MSLDDH_VD)=-YDVARS%SVD%T0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)
      IF (LNHX) PGMVTNDSL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,MSLDDH_NHX)=-YDCPG_DYN0%NHX(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)
    ENDIF
    DO JGFL=1,NUMFLDS
      IF (YCOMP(JGFL)%LADV) THEN
        PGFLTNDSL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,JGFL)=-PGFL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,YCOMP(JGFL)%MP)
      ENDIF
    ENDDO
  ELSE
    PGMVTNDSL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,MSLDDH_U)=-YDVARS%U%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)
    PGMVTNDSL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,MSLDDH_V)=-YDVARS%V%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)
    IF (NFTHER > 0) PGMVTNDSL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,MSLDDH_T)=-YDVARS%T%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)
    IF (LNHDYN)THEN !!! ?????? ky: LNHDYN or LNHEE there?
      PGMVTNDSL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,MSLDDH_PD)=-YDVARS%SPD%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)
    ENDIF
    IF (LNHDYN)THEN
      PGMVTNDSL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,MSLDDH_VD)=-YDVARS%SVD%T9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)
      IF (LNHX) PGMVTNDSL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,MSLDDH_NHX)=-YDCPG_DYN9%NHX(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG)
    ENDIF
    DO JGFL=1,NUMFLDS
      IF (YCOMP(JGFL)%LADV) THEN
        PGFLTNDSL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,JGFL)=-PGFL(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:NFLEVG,YCOMP(JGFL)%MP9)
      ENDIF
    ENDDO
  ENDIF

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG_GP', 1, ZHOOK_HANDLE)
END SUBROUTINE CPG_GP
