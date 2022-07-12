#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE CPG_DIA(YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDCPG_TND, YDCPG_MISC, YDCPG_GPAR, YDMF_PHYS, &
& YDCPG_DYN0, YDMF_PHYS_SURF, YDVARS, YDGMV, YDSURF, YDCFU, YDXFU, YDMODEL, KDDHI, &
& PGMV, PGFL, POMEGA, PEXT, PCOVPTOT, PGMVTNDSI, PGMVTNDHD, PGFLTNDHD, PDHCV,  &
& YDDDH, PFTCNS)

!**** *CPG_DIA* - Grid point calculations: diagnostics.

!     Purpose.
!     --------
!           Grid point calculations: diagnostics using non lagged dynamics
!           and non lagged physics (DDH, CFU, XFU).

!**   Interface.
!     ----------
!        *CALL* *CPG_DIA(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        KST       : first element of work.
!        KEND      : last element of work.
!        KSTC,KENDC: the same as KST,KEND but including zone C for ALADIN.
!        KSTGLO    : global offset.
!        KDDHI     : internal domain to which each point belongs (DDH).
!        LDCONFX   : (see in CPG)
!        LDLDIAB   : .T. if complete physics is activated.
!        LDLFSTEP  : .T. if first step.
!        LD_DFISTEP: 'D' -> DFI computations
!        KIBL      : index into YRGSGEOM/YRCSGEOM types in YDGEOMETRY
!        POROGL    : zonal component of "grad(surf orography)"
!        POROGM    : meridian component of "grad(surf orography)"
!        PGMV      : upper air GMV variables at time t and t-dt.
!        PGFL      : unified_treatment grid-point fields at t.
!        PRE0      : hydrostatic pressure "prehyd" at half levels at time t.
!        PRE0F     : hydrostatic pressure "prehyd" at full levels at time t.
!        PXYB0     : contains pressure depth, "delta", "alpha" at t.
!        PHI0      : geopotential height "gz" at t at half levels.
!        PHIF0     : geopotential height "gz" at t at full levels.
!        PHI0FL    : zonal component of "grad(gz)" at full levels.
!        PHI0FM    : meridian component of "grad(gz)" at full levels.
!        PRCP0     : contains "cp", "R" and "Kap=R/Cp" at t.
!        PCTY0     : contains vertical velocities, vertical integral of divergence at t.
!        POMEGA    : "omega" at full levels at t, including the "lrubc" and "delta m=1" effects.
!        PKENE0    : kinetic energy at t.
!        PRT0L     : zonal component of "grad(RT)" at full levels.
!        PRT0M     : meridian component of "grad(RT)" at full levels.
!        PRE0L     : zonal component of "grad(prehyds)" at full levels.
!        PRE0M     : meridian component of "grad(prehyds)" at full levels.
!        PNHPRE0F  : "pre" at full levels (time t).
!        PNHPRE0H  : "pre" at half levels (time t).
!        PQCHA0L   : zonal comp grad(log(pre/prehyd)).
!        PQCHA0M   : merid comp grad(log(pre/prehyd)).
!        PDHSF     : distribution of horizontal mean weights used for simplified radiation scheme.
!        ---------------------- output of aplpar ------------------------------
!        PALB      : model surface shortwave albedo.
!        PCAPE     : CAPE.
!        PCTOP     : pressure of top of convective cloud (diagnostic).
!        PCLCC     : convective cloud cover (diagnostic).
!        PCLCH     : high cloud cover (diagnostic).
!        PCLCL     : low cloud cover (diagnostic).
!        PCLCM     : medium cloud cover (diagnostic).
!        PCLCT     : total cloud cover (diagnostic).
!        PCLPH     : height (in meters) of the PBL.
!        PVEIN     : ventilation index in the PBL.
!        PCT       : thermical coefficient of soil-vegetation middle.
!        PDIFCQ    : convective flux of specific humidity (not rain/snow).
!        PDIFCQN   : convective flux of solid water (not rain/snow).
!        PDIFCQL   : convective flux of liquid water (not rain/snow).
!        PDIFCS    : convective flux of enthalpy (not rain/snow).
!        PDIFTQ    : turbulent flux (inc. "q" negative) of specific humidity.
!        PDIFTQN   : turbulent flux (inc. "q" negative) of solid water.
!        PDIFTQL   : turbulent flux (inc. "q" negative) of liquid water.
!        PDIFTS    : turbulent flux of enthalpy (or dry static energy).
!        PFCCQL    : convective condensation flux for liquid water.
!        PFCCQN    : convective condensation flux for ice.
!        PFCHSP    : heat flux from surface to deep soil.
!        PFCLL     : latent heat flux over liquid water (or wet soil).
!        PFCLN     : latent heat flux over snow (or ice).
!        PFCQNNG   : pseudo-flux of ice to correct for "qi"<0.
!        PFCQLNG   : pseudo-flux of liquid water to correct for "ql"<0.
!        PFCQNG    : pseudo-flux of water to correct for Q<0.
!        PFCS      : sensible heat flux at surface level.
!        PFCSQL    : stratiform condensation flux for liquid water.
!        PFCSQN    : stratiform condensation flux for ice.
!        PFEVL     : water vapour flux over liquid water (or wet soil).
!        PFEVN     : water vapour flux over snow (or ice) and frozen soil.
!        PFEVV     : evapotranspiration flux.
!        PFGEL     : freezing flux of soil water.
!        PFGELS    : freezing flux of soil water at surface level.
!        PFLWSP    : water flux from surface to deep soil.
!        PFONTE    : water flux corresponding to surface snow melt.
!        PFPLCL    : convective precipitation as rain.
!        PFPLCN    : convective precipitation as snow.
!        PFPLCG    : convective precipitation as graupel.
!        PFPLCH    : convective precipitation as hail.
!        PFPLSL    : stratiform precipitation as rain.
!        PFPLSN    : stratiform precipitation as snow.
!        PFPLSG    : stratiform precipitation as graupel.
!        PFPLSH    : stratiform precipitation as hail.
!        PMRT      : mean radiant temperature.
!        PFRMH     : mesospheric enthalpy flux.
!        PFRSO     : shortwave radiative flux.
!        PFRSOC    : shortwave clear sky radiative flux.
!        PFRSODS   : surface downwards solar flux.
!        PFRSOLU   : downward lunar flux at surface.
!        PFRSGNI   : surface global normal irradiance.
!        PFRSDNI   : surface direct normal irradiance.
!        PFRSOPS   : surface parallel solar flux.
!        PFRSOPT   : top parallel solar flux.
!        PFRTH     : longwave radiative flux.
!        PFRTHC    : longwave clear sky radiative flux.
!        PFRTHDS   : surface downwards IR flux.
!        PFTR      : transpiration flux.
!        PGZ0      : g*roughness lenght (current).
!        PGZ0H     : current g*thermal roughness lenght (if KVCLIV >= 8).
!        PNEB      : fractional cloudiness for radiation.
!        PQICE     : specific humidity of solid water for radiation.
!        PQLI      : specific humidity of liquid water for radiation.
!        PQRAIN    : specific humidity of rain for radiation.
!        PQSNOW    : specific humidity of snow for radiation.
!        PQS       : specific humidity at surface level.
!        PEXT      : extra diagnostics from physics.
!        PRH       : relative humidity.
!        PRHCLS    : relative humidity at 2 meters (diagnostic).
!        PRUISL    : run-off flux out the interception water-tank.
!        PRUISP    : run-off flux in soil.
!        PRUISS    : run-off flux at surface level.
!        PSTRCU    : convective flux of momentum "U".
!        PSTRCV    : convective flux of momentum "V".
!        PSTRDU    : gravity wave drag flux "U".
!        PSTRDV    : gravity wave drag flux "V".
!        PSTRTU    : turbulent flux of momentum "U".
!        PSTRTV    : turbulent flux of momentum "V".
!        PUCLS     : U-component of wind at 10 meters (diagnostic).
!        PVCLS     : V-component of wind at 10 meters (diagnostic).
!        PNUCLS    : U-component of neutral wind at 10 meters (diagnostic).
!        PNVCLS    : V-component of neutral wind at 10 meters (diagnostic).
!        PTCLS     : temperature at 2 meters (diagnostic).
!        PQCLS     : specific humidity at 2 meters (diagnostic).
!        PTPWCLS   : wet bulb temperature at 2 meters (diagnostic).
!        PUGST     : U-component of gusts (diagnostic).
!        PVGST     : V-component of gusts (diagnostic).
!        PVISICLD  : visibility due to cloud (ice and liquid) water
!        PVISIHYRO : visibility due to precipitation
!        ---------------------- end of output of aplpar -----------------------
!        PMOCON    : moisture convergence.
!        PFDIS     : enthalpy flux due to dissipation of kinetic energy. 
!        PFHPCL    : liquid water convective condensation enthalpy flux.
!        PFHPCN    : snow convective condensation enthalpy flux.
!        PFHPSL    : liquid water stratiform condensation enthalpy flux.
!        PFHPSN    : snow stratiform condensation enthalpy flux.
!        PFHSCL    : sensible heat flux due to liquid convective precipitations
!        PFHSCN    : sensible heat flux due to snow convective precipitations.
!        PFHSSL    : sensible heat flux due to liquid stratiform precipitations
!        PFHSSN    : sensible heat flux due to snow stratiform precipitations.
!        PFPFPSL   : flux of liquid resolved precipitation: the generation term.
!        PFPFPSN   : flux of solid resolved precipitation: the generation term.
!        PFPFPCL   : flux of liquid conv. precipitation: the generation term.
!        PFPFPCN   : flux of solid conv. precipitation: the generation term.
!        PFPEVPSL  : resolved precipitation flux due to evaporation.
!        PFPEVPSN  : resolved precipitation flux due to sublimation.
!        PFPEVPCL  : convective precipitation flux due to evaporation.
!        PFPEVPCN  : convective precipitation flux due to sublimation.
!        PGMVTNDSI : tendencies of semi-implicit scheme.
!        PGMVTNDHD : tendencies of horizontal diffusion scheme for GMV.
!        PGFLTNDHD : tendencies of horizontal diffusion for spectrally treated GFL.
!        PATND     : adiabatic Lagrangian tendencies.

!        PCOVPTOT  : precip fraction 

!     INPUT/OUTPUT:
!     -------------
!        PDHCV     : result array containing flux divergences and fluxes.
!                    Entering this routine this array is assumed to be
!                    pre-initialised with zeros.
!        PDRNSHF   : Derivative of non solar surface fluxes
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
!        ARPEGE documentation about Eulerian and semi-Lagrangian schemes.

! Author
! ------
!   13-Aug-2001 K. YESSAD after part 8 of CPG. 

! Modifications
! -------------
!   04-Mar-2008 Y. Seity  Cleaning MTS pictures (IR and WV) replaced by Fullpos
!   01-Oct-2008 O.riviere Introduction of new data flow for DDH
!   K. Yessad (Dec 2008): remove dummy CDLOCK + cleanings
!   18-May-2009 S. Riette : mean cls wind added
!   14-Oct-2009 A. Alias   ZUCLS, ZVCLS added in CUMCPLDM (E.Maisonnave)
!   04-Jan-2011 Y. Seity Add PDIAGH (AROME hail daignostic)
!   K. Yessad (Nov 2009): cleanings, DT/Dt now pre-computed in CPG_GP.
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   2011-10-04 J.M. Piriou: do not call CPPHDDH if flexible DDH.
!   K. Yessad (Nov 2011): various modifications.
!   M. Ahlgrimm  31-Oct-2011 add rain, snow and PEXTRA to DDH output
!   M. Ahlgrimm Apr 2014: Add lake variables and precip fraction to DDH output
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (July 2014): Move some variables, rename some variables.
!   R. El Khatib 07-Mar-2016 Pruning of ISP
!   Y. Seity (4-Sept-2017) add LDLRESET_GST* to replace KXXGSTTS
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   R. El Khatib 05-Jun-2018 computation of periods moved from cnt4 (OOPS refactoring)
!   R. Brozkova (Jul 2018): Dataflow for global normal irradiance and mean
!                           radiant temperature.
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
!   I. Etchevers (Sept 2018) : Add visisbilites
!   I. Etchevers (Janv 2019) : Add precipitation type
!   M. Hrastinski (Sep-2019): Dataflow for TKE and TTE terms in ALARO DDH (PFTCNS)
! End Modifications
!     ------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE MF_PHYS_TYPE_MOD   , ONLY : MF_PHYS_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_DYN_TYPE, CPG_MISC_TYPE, CPG_TND_TYPE, CPG_GPAR_TYPE
USE CPG_OPTS_TYPE_MOD   , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD,ONLY : MF_PHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMCFU             , ONLY : TCFU
USE YOMXFU             , ONLY : TXFU
USE YOMGMV             , ONLY : TGMV
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE DDH_MIX            , ONLY : RESET_DDHFLEX, TYP_DDH

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_BNDS_TYPE),INTENT(IN)   :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE),INTENT(IN)   :: YDCPG_OPTS
TYPE(CPG_TND_TYPE),INTENT(INOUT) :: YDCPG_TND
TYPE(CPG_MISC_TYPE),INTENT(INOUT):: YDCPG_MISC
TYPE(CPG_GPAR_TYPE),INTENT(INOUT):: YDCPG_GPAR
TYPE(MF_PHYS_TYPE),INTENT(INOUT) :: YDMF_PHYS
TYPE(CPG_DYN_TYPE),INTENT(INOUT) :: YDCPG_DYN0
TYPE(MF_PHYS_SURF_TYPE),INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
TYPE(TCFU)        ,INTENT(INOUT) :: YDCFU
TYPE(TXFU)        ,INTENT(INOUT) :: YDXFU
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KDDHI(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POMEGA(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_PHY_G%YRDPHY%NVEXTR)

REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMVTNDSI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,2+YDMODEL%YRML_GCONF%YRDIMF%NFTHER)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMVTNDHD(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,2+YDMODEL%YRML_GCONF%YRDIMF%NFTHER)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFLTNDHD(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDHCV(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_DIAG%YRMDDH%NDHCVSUN) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCOVPTOT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,1)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PFTCNS(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,6)
TYPE(TYP_DDH)     ,INTENT(INOUT) :: YDDDH 
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JLEV,  JROF, JCV
INTEGER (KIND=JPIM) :: JEXT



REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!#include "posinsub.intfb.h"

!     ------------------------------------------------------------------

#include "cpg_dia_ddh.intfb.h"
#include "cpg_dia_flu.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPG_DIA', 0, ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDVAB=>YDGEOMETRY%YRVAB, YDGSGEOM=>YDGEOMETRY%YRGSGEOM(YDCPG_BNDS%KBL), &
& YDCSGEOM=>YDGEOMETRY%YRCSGEOM(YDCPG_BNDS%KBL), YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY,              &
& YDSDDH=>YDMODEL%YRML_DIAG%YRSDDH, YDMCC=>YDMODEL%YRML_AOC%YRMCC, YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH,                                  &
& YDRIP=>YDMODEL%YRML_GCONF%YRRIP, YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY, YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL,                         &
& YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, YGFL=>YDMODEL%YRML_GCONF%YGFL, YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,                                &
& YDDIMF=>YDMODEL%YRML_GCONF%YRDIMF, YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR)

ASSOCIATE(NDIM=>YGFL%NDIM, YA=>YGFL%YA, YI=>YGFL%YI, YL=>YGFL%YL, YO3=>YGFL%YO3, YQ=>YGFL%YQ, YR=>YGFL%YR,                 &
& YS=>YGFL%YS, LCUMFU=>YDCFU%LCUMFU, NPROMA=>YDDIM%NPROMA, NPROMM=>YDDIM%NPROMM, NPROMNH=>YDDIM%NPROMNH,                   &
& NFTHER=>YDDIMF%NFTHER, NFLEVG=>YDDIMV%NFLEVG, NTSSG=>YDDPHY%NTSSG, NVEXTR=>YDDPHY%NVEXTR, LAGPHY=>YDEPHY%LAGPHY,         &
& LEPHYS=>YDEPHY%LEPHYS, NDIMGMV=>YDGMV%NDIMGMV, YT0=>YDGMV%YT0, LFLEXDIA=>YDLDDH%LFLEXDIA, LSDDH=>YDLDDH%LSDDH,           &
& LDDH_OMP=>YDLDDH%LDDH_OMP, NFRCPL=>YDMCC%NFRCPL, NDHCSSU=>YDMDDH%NDHCSSU, NDHCVSUN=>YDMDDH%NDHCVSUN,                     &
& NSTOP=>YDRIP%NSTOP, LHDPAS=>YDSDDH%LHDPAS, NDTPREC=>YDPHY%YRDPRECIPS%NDTPREC, NDTPREC2=>YDPHY%YRDPRECIPS%NDTPREC2,       &
& YSD_VF=>YDSURF%YSD_VF, YSD_VFD=>YDSURF%YSD_VFD, YSP_RR=>YDSURF%YSP_RR, YSP_RRD=>YDSURF%YSP_RRD, YSP_SB=>YDSURF%YSP_SB,   &
& YSP_SBD=>YDSURF%YSP_SBD, YSP_SG=>YDSURF%YSP_SG, YSP_SGD=>YDSURF%YSP_SGD, YSD_XP=>YDSURF%YSD_XP, YSD_XP2=>YDSURF%YSD_XP2, &
& LDPRECIPS=>YDPHY%LDPRECIPS, LDPRECIPS2=>YDPHY%LDPRECIPS2, LXFU=>YDXFU%LXFU, LDIRCLSMOD=>YDDPHY%LDIRCLSMOD,               &
& LTRAJPS=>YDSIMPHL%LTRAJPS, LSIMPH=>YDSIMPHL%LSIMPH, LMPHYS=>YDPHY%LMPHYS, LMPA=>YDARPHY%LMPA, MVQS=>YDPARAR%MVQS,        &
& MVTS=>YDPARAR%MVTS)
      
!     ------------------------------------------------------------------

!*       1.    DDH:
!              ----


CALL CPG_DIA_DDH (YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDCPG_TND, YDCPG_MISC, YDCPG_GPAR, YDMF_PHYS, YDCPG_DYN0, &
                & YDMF_PHYS_SURF, YDVARS, YDGMV, YDMODEL, KDDHI, PGMV, PGFL, POMEGA, PEXT, PCOVPTOT, PGMVTNDSI, &
                & PGMVTNDHD, PGFLTNDHD, PDHCV, YDDDH, PFTCNS)

CALL CPG_DIA_FLU (YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDCPG_MISC, YDMF_PHYS, YDCPG_DYN0, YDMF_PHYS_SURF, YDVARS, &
                & YDCFU, YDXFU, YDMODEL)

!     ------------------------------------------------------------------

!*       3.    RESETS DDH FLEXIBLE STRUCTURES:
!              -------------------------------

! resets ddh flexible structures if call only to phys

IF (YDCPG_OPTS%LCONFX.AND.LFLEXDIA.AND.(.NOT.LDDH_OMP)) CALL RESET_DDHFLEX

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG_DIA', 1, ZHOOK_HANDLE)
END SUBROUTINE CPG_DIA
