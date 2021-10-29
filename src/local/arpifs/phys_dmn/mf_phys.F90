#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE MF_PHYS(YDGEOMETRY,YDCPG_MISC,YDCPG_PHY0,YDCPG_PHY9,YDMF_PHYS,YDCPG_DYN0,YDCPG_DYN9,&
 & YDMF_PHYS_SURF,YDVARS,YDGMV,YDSURF,YDCFU,YDXFU,YDMODEL,&
 !---------------------------------------------------------------------
 ! - INPUT and INOUT.
 & KBL,KGPCOMP,KST,KEND,KGL1,KGL2,KSTGLO,&
 & LDCONFX,PDTPHY,&
 & KIBL,&
 & PGFL,&
 & PKOZO,PGP2DSDT,&
 & PB1,PB2,PGMVT1,PGFLT1,&
 & PGPAR,&
 !---------------------------------------------------------------------
 ! - OUTPUT .
 & PTRAJ_PHYS,YDDDH,PFTCNS)

!**** *MF_PHYS* METEO-FRANCE PHYSICS.

!     Purpose.
!     --------
!         Call METEO-FRANCE physics and physical tendencies.

!**   Interface.
!     ----------
!        *CALL* *MF_PHYS(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        KBL       : NPROMA-packets number
!        KGPCOMP   : total number of grid points in the domain
!        KST       : first element of work.
!        KEND      : last element of work.
!        KGL1,KGL2 : first and last latitude of computations.
!                    - bounds NGPTOT-packets for DM calculations.
!        KSTGLO    : global offset.
!        LDCONFX   : (see in CPG)
!        PDTPHY    : timestep used in the physics.
!        KIBL      : index into YRCSGEOM/YRGSGEOM types in YDGEOMETRY
!        POROGL,POROGM: components of grad(orography).
!        PCUCONVCA : CA array for interaction with the physics
!        PNLCONVCA : CA array for interaction with the physics
!        PGMV      : GMV at time t and t-dt.
!        PGMVS     : GMVS at time t and t-dt.
!        PGFL      : GFL at time t and t-dt.
!        PWT0      : w-wind time t.
!        PWT0L     : zonal derivative of w-wind at time t.
!        PWT0M     : merid derivative of w-wind at time t.
!        PRCP0     : contains "cp", "R" and "Kap=R/Cp" at t.
!        PHI0      : geopotential height at half levels at time t.
!        PHIF0     : geopotential height at full levels at time t.
!        PRE0      : hydrostatic pressure "prehyd" at half levels at time t.
!        PRE0F     : hydrostatic pressure "prehyd" at full levels at time t.
!        PREPHY0   : input pressure "pre" for AROME at half levels at time t.
!        PREPHY0F  : input pressure "pre" for AROME at full levels at time t.
!        PXYB0     : contains pressure depth, "delta", "alpha" at time t.
!        PRKQVH    : Rasch-Kristjansson scheme - water vapour tendency
!        PRKQCH    : Rasch-Kristjansson scheme - condensates tendency
!        PWT9      : Vertical wind time t-dt.
!        PRCP9     : contains "cp", "R" and "Kap=R/Cp" at t-dt.
!        PHI9      : geopotential height at half levels at time t-dt.
!        PHIF9     : geopotential height at full levels at time t-dt.
!        PRE9      : hydrostatic pressure "prehyd" at half levels at time t-dt.
!        PRE9F     : hydrostatic pressure "prehyd" at full levels at time t-dt.
!        PREPHY9   : input pressure "pre" for AROME at half levels at time t-dt.
!        PREPHY9F  : input pressure "pre" for AROME at full levels at time t-dt.
!        PXYB9     : contains pressure depth, "delta", "alpha" at time t-dt.
!        PKOZO     : fields for photochemistery of ozon.
!        PGP2DSDT  : stochastic physics random pattern.
!        PGRADH_PHY: horizontal gradients for physics

!     INPUT/OUTPUT:
!     -------------
!        PCTY0     : contains vertical velocities, vertical integral of divergence at t.
!        PB1       : "SLB1"-buffer, used for interpolations in the SL scheme.
!        PB2       : "SLB2"-buffer.
!        PGFLT1    : GFL t+dt
!        PGPAR     : surface fields for AROME.
!        PGDEOSI   : DESCENDING INCREMENTAL OPTICAL DEPTHS, SOLAR
!        PGUEOSI   : ASCENDING  INCREMENTAL OPTICAL DEPTHS, SOLAR
!        PGMU0     : COSINE OF SOLAR ZENITH ANGLE, APPROXIMATE ACTUAL VALUE
!        PGMU0_MIN : COSINE OF SOLAR ZENITH ANGLE, MIN VALUE
!        PGMU0_MAX : COSINE OF SOLAR ZENITH ANGLE, MAX VALUE
!        PGDEOTI   : descending incremental optical depths, dB/dT(T0) weights
!        PGDEOTI2  : descending incremental optical depths, B weights with
!                    linear T_e correction
!        PGUEOTI   : ascending incremental optical depths, dB/dT(T0) weights
!        PGUEOTI2  : ascending incremental optical depths, B weights with
!                    linear T_e correction
!        PGEOLT    : local optical depths, dB/dT(T0) weights
!        PGEOXT    : maximum optical depths for EBL-EAL, dB/dT(T0) weights
!        PGRPROX   : correction term for adjacent exchanges
!        PGMIXP    : non-statistical weights for bracketing
!        PGFLUXC   : out of bracket part of clearsky EBL, resp. EBL-EAL flux
!        PGRSURF   : corrective ratio for surface cts contribution

!     OUTPUT:
!     -------
!        PDHSF     : distribution of horizontal mean weights used for
!                    simplified radiation scheme.
!        ---------------------- output of aplpar ------------------------------
!        PALBDG      : modele surface shortwave albedo (diagnostic).
!        PCAPE     : CAPE.
!        PCTOP     : top of convective nebulosity (diagnostic).
!        PCLCC     : convective cloud cover (diagnostic).
!        PCLCH     : high cloud cover (diagnostic).
!        PCLCL     : low cloud cover (diagnostic).
!        PCLCM     : medium cloud cover (diagnostic).
!        PCLCT     : total cloud cover (diagnostic).
!        PCLPH     : height (in meters) of the PBL.
!        PVEIN     : ventilation index in the PBL.
!        PCT       : thermical coefficient of soil-vegetation middle.
!        PDIFCQ    : convective flux of specific humidity (not rain/snow).
!        PDIFCQI   : convective flux of solid water (not rain/snow).
!        PDIFCQL   : convective flux of liquid water (not rain/snow).
!        PDIFCS    : convective flux of enthalpy (not rain/snow).
!        PDIFTQ    : turbulent flux (inc. "q" negative) of specific humidity.
!        PDIFTQI   : turbulent flux (inc. "q" negative) of solid water.
!        PDIFTQL   : turbulent flux (inc. "q" negative) of liquid water.
!        PDIFTS    : turbulent flux of enthalpy (or dry static energy).
!        PFCCQL    : convective condensation flux for liquid water.
!        PFCCQN    : convective condensation flux for ice.
!        PFCHOZ    : ozon photo-chemical flux.
!        PFPFPSL   : flux of liquid resol. precipitation: the generation term. 
!        PFPFPSN   : flux of solid resolved precipitation: the generation term. 
!        PFPFPCL   : flux of liquid conv. precipitation: the generation term. 
!        PFPFPCN   : flux of solid conv. precipitation: the generation term. 
!        PFPEVPSL  : resolved precipitation flux due to evaporation.
!        PFPEVPSN  : resolved precipitation flux due to sublimation.
!        PFPEVPCL  : convective precipitation flux due to evaporation.
!        PFPEVPCN  : convective precipitation flux due to sublimation.
!        PFCHSP    : heat flux from surface to deep soil.
!        PFCLL     : latent heat flux over liquid water (or wet soil).
!        PFCLN     : latent heat flux over snow (or ice).
!        PFCQING   : pseudo-flux of ice to correct for "qi"<0.
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
!        PFPLCHL   : convective precipitation as hail.
!        PFPLSL    : stratiform precipitation as rain.
!        PFPLSN    : stratiform precipitation as snow.
!        PFPLSG    : stratiform precipitation as graupel.
!        PFPLSHL   : stratiform precipitation as hail.
!        PMRT      : mean radiant temperature.
!        PFRMH     : mesospheric enthalpy flux.
!        PFRSO     : shortwave radiative flux.
!        PFRSOC    : shortwave clear sky radiative flux.
!        PFRSODS   : surface downwards solar flux.
!        PFRSOLU   : downward lunar flux at surface.
!        PFRSGNI   : Global normal irradiance
!        PFRSDNI   : Direct normal irradiance
!        PFRSOPS   : surface parallel solar flux.
!        PFRSOPT   : top parallel solar flux.
!        PFRTH     : longwave radiative flux.
!        PFRTHC    : longwave clear sky radiative flux.
!        PFRTHDS   : surface downwards IR flux.
!        PFTR      : transpiration flux.
!        PGZ0      : g*roughness length (current).
!        PGZ0H     : current g*thermal roughness length (if KVCLIV >=8).
!        PNEB      : fractional cloudiness for radiation.
!        PQCLS     : specific humidity at 2 meters (diagnostic).
!        PQICE     : specific humidity of solid water for radiation.
!        PQLI      : specific humidity of liquid water for radiation.
!        PQS       : specific humidity at surface level.
!        PRH       : relative humidity.
!        PRHCLS    : relative humidity at 2 meters (diagnostic).
!        PRUISL    : run-off flux out the interception water-tank.
!        PRUISP    : run-off flux in soil.
!        PRUISS    : run-off flux at surface level.
!        PSTRCU    : convective flux of momentum "U".
!        PSTRCV    : convective flux of momentum "V".
!        PSTRDU    : gravity wave drag flux "U".
!        PSTRDV    : gravity wave drag flux "V".
!        PSTRMU    : mesospheric flux for "U"-momentum.
!        PSTRMV    : mesospheric flux for "V"-momentum.
!        PSTRTU    : turbulent flux of momentum "U".
!        PSTRTV    : turbulent flux of momentum "V".
!        PDIFCQLC to PFCNEGQSC:
!        PUCLS     : U-component of wind at 10 meters (diagnostic).
!        PVCLS     : V-component of wind at 10 meters (diagnostic).
!        PNUCLS    : U-component of neutral wind at 10 meters (diagnostic).
!        PNVCLS    : V-component of neutral wind at 10 meters (diagnostic).
!        PTCLS     : temperature at 2 meters (diagnostic).
!        PTPWCLS   : wet-bulb temperature at 2 meters (diagnostic)
!        PUGST     : U-component of gusts (diagnostic).
!        PVGST     : V-component of gusts (diagnostic).
!        PDERNSHF  : derivative of the non solar surface with respect to Tsurf
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
!        PTENDU    : "U"-wind tendency due to physics.
!        PTENDV    : "V"-wind tendency due to physics.
!        PQSOL     : surface specific humidity used in case "delta m=1".
!        PFCQRNG   : pseudo-flux of rain to correct for Q<0
!        PFCQSNG   : pseudo-flux of snow to correct for Q<0
!        PDIAGH    : Add Hail diagnostic PDIAGH (AROME)
!        PFLASH    : Add lightening density (fl/ km2 /s )
!        PVISICLD  : Visibility due to ice and/or water cloud
!        PVISIHYDRO : Vsibility due to precipitations(rain, graupel, snow)
!        PMXCLWC   : Cloud Water Liquid Content at HVISI meters

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!       2000-12-04: F. Bouyssel & J.M. Piriou

!     Modifications.
!     --------------
!       04-Mar-2009 A.Alias : call CPTEND/INITAPLPAR modified to add
!                         Humidity Mesopheric flux (ZFRMQ).
!                     and IVCLIA removed and call to CPNUDG modified as
!                         Nuding mask is now in SD_VF group
!                         call HL_APLPAR modified to add PFCQNG for acdifus
!                         call APL_AROME modified to add Sulfate/Volcano aerosols for radaer
!       K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!       2009-10-15 Y. Bouteloup : Store radiative cloud water and ice in GFL (YIRAD and YLRAD)
!       F. Vana   15-Oct-2009 : NSPLTHOI option
!       K.Yessad (Feb 2010): use YM_RADTC and RFORADTC
!       2010-03-26 Y. Bouteloup : Store radiative cloud water and ice in GFL (AROME case)
!       2010-04-26 Y. Bouteloup : Only one call to cputqy, cputqys and cputqy_arome
!            This need the use of ZTENDGFL as argument of cptend, cptend_new and apl_arome.
!       2010-05-11 F. Bouyssel : Use of PINDX, PINDY
!       2010-05-28 C. Geijo    : Fix error in IPTR array element referencing 
!       2010-06-21 O.Riviere/F. Bouyssel : Fix to have Ts evolving in Fa files with Surfex
!       Dec 2010 A.Alias   : ZMU0N added to call CPPHINP/APLPAR/APL_AROME/HL_APLPAR
!                            CALL to CPNUDG with or with LMSE (A.Voldoire)
!       K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!       L. Bengtsson-Sedlar & F. Vana 18-Feb-2011 : CA scheme for convection
!       F. Vana   22-Feb-2011 : 3D turbulence
!       2011-02-01 M. Mokhtari: Add LMDUST and PEXTT9 and PEXTT0 IN APLPAR
!                             (treatment of the desert aerosols) 
!       2011-03 A.Alias  : new argument to  for sunshine hours YSD_VD%YSUND 
!                      CPNUDG if LMSE=.T. or LMSE=.F. (bugfix)
!                      debug ozone GFL (IPO3) (D. St-Martin)
!                      Humidity Mesopheric flux (ZFRMQ) added in CPTEND_NEW
!       F.Bouyssel (26-03-2011): Fix to have Snow in hist file with surfex
!       2011-06: M. Jerczynski - some cleaning to meet norms
!       E. Bazile 2011-08-26 : Output for MUSC 1D with LFA files with WRITEPHYSIO
!         used previously for extracting profiles from 3D (now also available for AROME).
!       K. Yessad (Dec 2011): use YDOROG, YDGSGEOM and YDCSGEOM.
!       2011-11-21 JF Gueremy : dry convective adjustment (LAJUCV)
!       F. Vana  26-Jan-2012 : historic Qs for TOM's BBC.
!       F.Bouttier Jul 2012: stochastic physics for AROME
!       Z. SASSI  : 07-Mar-2013   INITIALIZING THE WEIGHT VECTORS PDHSF(NPROMA)
!       [DISTRIBUTION OF HORIZONTAL MEANS WEIGHTS]
!       F. Vana  28-Nov-2013 : Redesigned trajectory handling
!       T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!       2013-11, D. Degrauwe: Flexible interface CPTEND_FLEX.
!       2013-11, J. Masek: Passing intermittency arrays for ACRANEB2.
!       K. Yessad (July 2014): Move some variables.
!       2016-04, J. Masek: Passing sunshine duration to APL_AROME.
!       2016-09, M. Mokhtari & A. Ambar: replacement of ZEXT and ZEZDIAG by PGFL
!                                        in aplpar.F90 argument.
!       2016-10, P. Marguinaud : Port to single precision
!       K. Yessad (Dec 2016): Prune obsolete options.
!       K. Yessad (June 2017): Introduce NHQE model.
!       2017-09, J. Masek: Shifted dimensioning of PGMU0.
!       K. Yessad (Feb 2018): remove deep-layer formulations.
!       K. Yessad (Apr 2018): introduce key L_RDRY_VD (ensure consistent definition of "dver" everywhere).
!       2018-09, F. Duruisseau: add rconv and sconv in gfl for bayrad
!       2018-09, R. Brozkova: Passing of diagnostic hail, global normal
!         irradiance and mean radiant temperature from APLPAR.
!       2018-09, D. St-Martin : add NOGWD inputs in aplpar
!       2018-09, M. Michou : add ARPEGE-Climat chemistry call in aplpar  
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
!       2019-05, I. Etchevers : add visibilities and precipitation type
!   R. El Khatib 27-02-2019 memory bandwidth savings.
!   R. El Khatib 30-Oct-2018 IMAXDRAFT
!       2019-09, M. Hrastinski: Dataflow for TKE and TTE terms in ALARO DDH (PFTCNS).
!       2019-09, J. Masek: Modified call to APL_AROME (added argument NFRRC).
!       2019-12, Y. Bouteloup: Introduction of ZTENDU and ZTENDV for computation of ZDEC in cputqy
!                diferent from PTENDU and PTENDV in the case of the use of Tiedtke scheme to avoid double counting
!       2020-12, U. Andrae : Introduce SPP for HARMONIE-AROME
!       2021-01, R. Brozkova: ALARO graupel fix.
! End Modifications
!-------------------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE MF_PHYS_TYPE_MOD   , ONLY : MF_PHYS_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_DYN_TYPE, CPG_PHY_TYPE, &
                              & CPG_MISC_TYPE
USE MFPHYS_SURFACE_TYPE_MOD,ONLY : MFPHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMGMV             , ONLY : TGMV
USE YOMCFU             , ONLY : TCFU
USE YOMXFU             , ONLY : TXFU
USE TYPE_MODEL         , ONLY : MODEL
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK
USE SC2PRG_MOD         , ONLY : SC2PRG

USE YOMCT0             , ONLY : LSLAG, LTWOTL, LNHDYN, LAROME, LSFORCS, LNHQE,LCORWAT
USE YOMCT3             , ONLY : NSTEP 
USE YOMCVER            , ONLY : LVERTFE  ,LVFE_GWMPA 
USE YOMDYNA            , ONLY : LGWADV, L3DTURB, L_RDRY_VD
USE YOMNUD             , ONLY : NFNUDG   ,LNUDG 
USE YOMSNU             , ONLY : XPNUDG
USE MODULE_RADTC_MIX   , ONLY : YM_RADTC
USE YOMSCM             , ONLY : LGSCM
USE YOMCST             , ONLY : RG, RD
USE YOMCHET            , ONLY : GCHETN
USE INTDYN_MOD         , ONLY : YYTCTY0  ,YYTRCP0  ,YYTRCP9  ,YYTXYB0_PHY,YYTXYB9_PHY
USE YOMLSFORC          , ONLY : LMUSCLFA
USE YOMSPSDT           , ONLY : YSPPT
USE YOMTRAJ            , ONLY : TRAJ_PHYS_TYPE, LPRTTRAJ
USE YOMLUN             , ONLY : NULOUT

USE DDH_MIX            , ONLY : TYP_DDH
USE INTFLEX_MOD        , ONLY : LINTFLEX, TYPE_INTPROCSET, NEWINTPROCSET, CLEANINTPROCSET
USE SPP_MOD , ONLY : YSPP_CONFIG,YSPP
USE APLPAR_STATE_TYPE_MOD, ONLY : APLPAR_STATE_TYPE
!     -------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_MISC_TYPE),INTENT(INOUT):: YDCPG_MISC
TYPE(CPG_PHY_TYPE),INTENT(INOUT), TARGET :: YDCPG_PHY0
TYPE(CPG_PHY_TYPE),INTENT(INOUT), TARGET :: YDCPG_PHY9
TYPE(MF_PHYS_TYPE),INTENT(INOUT) :: YDMF_PHYS
TYPE(CPG_DYN_TYPE),INTENT(INOUT), TARGET :: YDCPG_DYN0
TYPE(CPG_DYN_TYPE),INTENT(INOUT), TARGET :: YDCPG_DYN9
TYPE(MFPHYS_SURF_TYPE),INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
TYPE(TCFU)        ,INTENT(INOUT) :: YDCFU
TYPE(TXFU)        ,INTENT(INOUT) :: YDXFU
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KBL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPCOMP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGL1 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGL2 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTGLO 
LOGICAL           ,INTENT(IN)    :: LDCONFX
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTPHY 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL


REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKOZO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG*YDMODEL%YRML_PHY_G%YRDPHY%NVCLIS+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP2DSDT(YDGEOMETRY%YRDIM%NPROMA,YSPPT%YGPSDT(1)%NG2D)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%YT1%NDIM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGPAR(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_PHY_MF%YRPARAR%NGPAR+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFTCNS(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,6)
TYPE (TRAJ_PHYS_TYPE), INTENT(INOUT) :: PTRAJ_PHYS
TYPE(TYP_DDH)     ,INTENT(INOUT) :: YDDDH

!     ------------------------------------------------------------------
LOGICAL :: LLDIAB
LOGICAL :: LL_SAVE_PHSURF
LOGICAL :: LLXFUMSE

INTEGER(KIND=JPIM) :: IFIELDSS
INTEGER(KIND=JPIM) :: IBLK
INTEGER(KIND=JPIM) :: IPQ,IPO3,ITDIA,IPTREXT,IPTR_CONT,IEFB1,IEFB2,IEFB3
INTEGER(KIND=JPIM) :: IPGFL(YDMODEL%YRML_GCONF%YGFL%NUMFLDS),IPTR(YDMODEL%YRML_GCONF%YGFL%NUMFLDS)

INTEGER(KIND=JPIM) :: INSTEP_DEB,INSTEP_FIN
INTEGER(KIND=JPIM) :: JLEV, JGFL
INTEGER(KIND=JPIM) :: JROF
INTEGER(KIND=JPIM) :: ISLB1U9  ,ISLB1V9  ,ISLB1T9  ,ISLB1GFL9, ISLB1VD9

INTEGER(KIND=JPIM) :: ICLPH(YDGEOMETRY%YRDIM%NPROMA)          ! cf. KCLPH in APLPAR.

!     --- UPPER AIR PHYSICAL TENDENCIES.
REAL(KIND=JPRB) :: ZTENDH(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)     ! Enthalpy tendency.
REAL(KIND=JPRB) :: ZTENDQ(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)     ! Moisture tendency.
REAL(KIND=JPRB) :: ZTENDPTKE(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)  ! Pseudo progn. TKE

! GFL tendencies for APL_AROME (assumes YDMODEL%YRML_GCONF%YGFL%NUMFLDS>=YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
! for now, use Jovi's trick :
REAL(KIND=JPRB) :: ZTENDGFLR(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG,0:YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
REAL(KIND=JPRB) :: ZTENDGFL(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NUMFLDS)   ! GFL tendencies

REAL(KIND=JPRB) :: ZUDOM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZUDAL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDDOM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZDDAL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZUNEBH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZENTCH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     --- UPPER AIR PHYSICAL TENDENCIES FOR AROME.
!       (the previous one are not used in AROME)
REAL(KIND=JPRB) :: ZTENDT (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)        ! temperature tendency
REAL(KIND=JPRB) :: ZTENDW (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)        ! W  tendency
REAL(KIND=JPRB) :: ZTENDD (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)        ! d  tendency
REAL(KIND=JPRB) :: ZTENDEXT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)      ! GFL EXTRA tendency
REAL(KIND=JPRB) :: ZTENDEXT_DEP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)  ! GFL EXTRA tendency
REAL(KIND=JPRB) :: ZDIFEXT(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)     ! Extra-GFL fluxes.

REAL(KIND=JPRB) :: ZTENDU (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! U tendency without deep convection contribution
REAL(KIND=JPRB) :: ZTENDV (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! V tendency without deep convection contribution
!     --- SURFACE AND DEEP RESERVOIR PHYSICAL TENDENCIES.
REAL(KIND=JPRB) :: ZTDTS(YDGEOMETRY%YRDIM%NPROMM)             ! Surface temperature tendency.
REAL(KIND=JPRB) :: ZTDTP(YDGEOMETRY%YRDIM%NPROMM,YDSURF%YSP_SBD%NLEVS)        ! Deep temperature tendency.
REAL(KIND=JPRB) :: ZTDWS(YDGEOMETRY%YRDIM%NPROMM)             ! Surface water res. tendency.
REAL(KIND=JPRB) :: ZTDWP(YDGEOMETRY%YRDIM%NPROMM)             ! Deep water res. tendency.
REAL(KIND=JPRB) :: ZTDWL(YDGEOMETRY%YRDIM%NPROMM)             ! Interception water res. tendency
REAL(KIND=JPRB) :: ZTDSNS(YDGEOMETRY%YRDIM%NPROMM)            ! Snow res. tendency. 
REAL(KIND=JPRB) :: ZTDWPI(YDGEOMETRY%YRDIM%NPROMM)            ! Deep ice res. tendency.
REAL(KIND=JPRB) :: ZTDWSI(YDGEOMETRY%YRDIM%NPROMM)            ! Surface ice res. tendency.
REAL(KIND=JPRB) :: ZTDALBNS(YDGEOMETRY%YRDIM%NPROMA)          ! Snow albedo tendency.
REAL(KIND=JPRB) :: ZTDRHONS(YDGEOMETRY%YRDIM%NPROMA)          ! Snow density tendency.

!     --- FLUXES FROM PARAMETERISATIONS AND TENDENCIES COMP.
REAL(KIND=JPRB) :: ZFP(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)      ! Total rainfall flux.
REAL(KIND=JPRB) :: ZFPLCH(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)   ! convective precipitation flux,
                                             ! local
REAL(KIND=JPRB) :: ZFPLSH(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)   ! stratiform precipitation flux,
                                             ! local
REAL(KIND=JPRB) :: ZFPLSN(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)   ! all solid stratiform precipitation flux,
                                             ! local
REAL(KIND=JPRB) :: ZFTKE(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)    ! TKE flux.
REAL(KIND=JPRB) :: ZFTKEI(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)    ! TKE flux.
REAL(KIND=JPRB) :: ZFEFB1(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)    ! EFB1 flux.
REAL(KIND=JPRB) :: ZFEFB2(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)    ! EFB2 flux.
REAL(KIND=JPRB) :: ZFEFB3(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)    ! EFB3 flux.

!     --- Fields of potential use for HIRLAM .
REAL(KIND=JPRB) :: ZCVGQL(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)     ! local array for cloud water tendency 
REAL(KIND=JPRB) :: ZCVGQI(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)     ! local array for cloud ice tendency
REAL(KIND=JPRB) :: ZCVGT(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)      ! local array for temperature tendency

!     --- MISCELLANEOUS PARAMETERISATIONS, 2D ARRAYS ---
REAL(KIND=JPRB) :: ZCVGQ(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)      ! convergence of humidity
                                             ! ("Kuo" condition).
REAL(KIND=JPRB) :: ZFHP(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)     ! Total enthalpy flux
                                             ! + sensible heat flux.
REAL(KIND=JPRB) :: ZLCVQ(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)      ! limited physical contribution to
                                             ! the diagnostic of moisture
                                             ! convergence (remove
                                             ! large scale precip. contrib.).
REAL(KIND=JPRB) :: ZLSCPE(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)     ! cf. PLSCPE in APLPAR.
REAL(KIND=JPRB) :: ZLH(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)        ! cf. PLH in APLPAR.
REAL(KIND=JPRB) :: ZQW(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)        ! cf. PQW in APLPAR.
REAL(KIND=JPRB) :: ZTW(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)        ! cf. PTW in APLPAR.

REAL(KIND=JPRB) :: ZFRMQ(YDGEOMETRY%YRDIM%NPROMM,0:YDGEOMETRY%YRDIMV%NFLEVG)    ! cf. MESOSPHERIC humidity flux in APLPAR.

!     --- 1D DIAGNOSTIC FIELDS, SURFACE FLUXES
REAL(KIND=JPRB) :: ZCD(YDGEOMETRY%YRDIM%NPROMM)               ! cf. PCD in APLPAR.
REAL(KIND=JPRB) :: ZCDN(YDGEOMETRY%YRDIM%NPROMM)              ! cf. PCDN in APLPAR.
REAL(KIND=JPRB) :: ZCH(YDGEOMETRY%YRDIM%NPROMM)               ! cf. PCH in APLPAR.
REAL(KIND=JPRB) :: ZEMIS(YDGEOMETRY%YRDIM%NPROMM)             ! cf. PEMIS in APLPAR.
REAL(KIND=JPRB) :: ZFEVI(YDGEOMETRY%YRDIM%NPROMM,YDMODEL%YRML_PHY_G%YRDPHY%NTSSG+1)     ! cf. PFEVI in APLPAR.
REAL(KIND=JPRB) :: ZNEIJ(YDGEOMETRY%YRDIM%NPROMM)             ! cf. PNEIJ in APLPAR.
REAL(KIND=JPRB) :: ZVEG(YDGEOMETRY%YRDIM%NPROMM)              ! cf. PVEG in APLPAR.
REAL(KIND=JPRB) :: ZQSAT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)      ! specific humidity at saturation.
REAL(KIND=JPRB) :: ZQSATS(YDGEOMETRY%YRDIM%NPROMA)            ! cf. PQSATS in APLPAR.
REAL(KIND=JPRB) :: ZQS1(YDGEOMETRY%YRDIM%NPROMM)              ! cf. PQS1 in APLPARS.

!     --- DIAGNOSTIC FIELDS STATE OF SURFACE AIR ---
REAL(KIND=JPRB) :: ZC1(YDGEOMETRY%YRDIM%NPROMM)               ! cf. PC1 in APLPAR.
REAL(KIND=JPRB) :: ZC2(YDGEOMETRY%YRDIM%NPROMM)               ! cf. PC2 in APLPAR.
REAL(KIND=JPRB) :: ZCPS(YDGEOMETRY%YRDIM%NPROMM)              ! cf. PCPS in APLPAR.
REAL(KIND=JPRB) :: ZLHS(YDGEOMETRY%YRDIM%NPROMM)              ! cf. PLHS in APLPAR.
REAL(KIND=JPRB) :: ZRS(YDGEOMETRY%YRDIM%NPROMM)               ! cf. PRS in APLPAR.

!     --- RADIATION COEFFICIENTS FOR SIMPLIFIED PHYSICS IN GRID-POINT ---
REAL(KIND=JPRB) :: ZAC(YDGEOMETRY%YRDIM%NPROMM,(YDGEOMETRY%YRDIMV%NFLEVG+1)*(YDGEOMETRY%YRDIMV%NFLEVG+1))   ! Curtis matrix.
REAL(KIND=JPRB) :: ZAC_HC(YDGEOMETRY%YRDIMV%NFLEVG+1,YDGEOMETRY%YRDIMV%NFLEVG+1)           ! horizontally-constant field for ZAC.
REAL(KIND=JPRB) :: ZRADTC(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG,YM_RADTC%NDIM) ! others than Curtis matrix.

!     --- GEOMETRY FOR RADIATION ---
REAL(KIND=JPRB) :: ZMMU0(YDGEOMETRY%YRDIM%NPROMA)  ! mean solar angle (for given day for
                                  ! simpl. rad. scheme)
REAL(KIND=JPRB) :: ZMU0(YDGEOMETRY%YRDIM%NPROMA)   ! local cosine of instantaneous solar zenith
                                  ! angle.
REAL(KIND=JPRB) :: ZMU0LU(YDGEOMETRY%YRDIM%NPROMM) ! local cosine of instantaneous lunar zenith
                                  ! angle.
REAL(KIND=JPRB) :: ZMU0M(YDGEOMETRY%YRDIM%NPROMA)  ! local cosine of averaged solar zenith angle
REAL(KIND=JPRB) :: ZMU0N(YDGEOMETRY%YRDIM%NPROMA)  ! same as ZMU0 for next time step (used for YDMODEL%YRML_PHY_MF%YRARPHY%LMSE)    

!     ---FOR AROME PHYSICS  ---
REAL(KIND=JPRB) :: ZDT   !pour cputqy_arome, a changer peut etre plus tard...
REAL(KIND=JPRB) :: ZGWT1(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)     ! vertical velocity calculated by 
                                              ! cputqy_arome before convertion in d
REAL(KIND=JPRB) :: ZTT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)       ! Temperature at t1

! ZRTT1: appropriate version of R*T at t1 for gnhgw2svd
!  Version of R must be consistent with definition of vertical divergence.
REAL(KIND=JPRB) :: ZRTT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     --- Buffers to save the initial value of   ---
!     --- some pseudo-historical surface buffers ---
REAL(KIND=JPRB) :: ZHV(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZGZ0F(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZGZ0HF(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZPBLH(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZFHPS(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZQSH (YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB) :: ZUDGRO(YDGEOMETRY%YRDIM%NPROMM)

!     --- FOR BAYRAD ALLSKY FRAMEWORK ---
REAL(KIND=JPRB) :: ZQRCONV(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZQSCONV(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)

!  Horizontal exchange coefficients for 3D turbulence
REAL(KIND=JPRB) :: ZKUROV_H(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZKTROV_H(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!  Empty arrays for 3TL scheme (to be coded later)
!!later REAL(KIND=JPRB) :: ZDIVT9(NPROMA,NFLEVG)
!!later REAL(KIND=JPRB) :: ZUT9L(NPROMA,NFLEVG)
!!later REAL(KIND=JPRB) :: ZVT9L(NPROMA,NFLEVG)
REAL(KIND=JPRB) :: ZWT9L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZWT9M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDPRECIPS(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS%NDTPREC)
REAL(KIND=JPRB) :: ZDPRECIPS2(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS%NDTPREC2)

REAL(KIND=JPRB) :: ZTAUX(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 &  ZDTAJU(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZEDR(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

! To save Tt for NHQE model
REAL(KIND=JPRB) :: ZTT0_SAVE(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTT0L_SAVE(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTT0M_SAVE(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTT9_SAVE(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)

! required for INTFLEX
TYPE(TYPE_INTPROCSET) :: YLPROCSET

! SPP
REAL(KIND=JPRB) :: ZGP2DSPP(YDGEOMETRY%YRDIM%NPROMA,YSPP%N2D)

INTEGER(KIND=JPIM) :: IMAXDRAFT

INTEGER(KIND=JPIM) :: IPTRLIMA
INTEGER(KIND=JPIM) :: IRR ! pointer of 1st hydrometeors in ZTENDGFLR
INTEGER(KIND=JPIM) :: IPTRTKE ! pointer of TKE in ZTENDGFLR

REAL(KIND=JPRB), POINTER :: ZMCOR(:,:), ZMRAB3C(:,:), ZMRAB3N(:,:), ZMRAB4C(:,:)
REAL(KIND=JPRB), POINTER :: ZMRAB4N(:,:), ZMRAB6C(:,:), ZMRAB6N(:,:), ZMRAT1C(:,:)
REAL(KIND=JPRB), POINTER :: ZMRAT1N(:,:), ZMRAT2C(:,:), ZMRAT2N(:,:), ZMRAT3C(:,:)
REAL(KIND=JPRB), POINTER :: ZMRAT3N(:,:), ZMRAT4C(:,:), ZMRAT4N(:,:), ZMRAT5C(:,:)
REAL(KIND=JPRB), POINTER :: ZMRAT5N(:,:)
REAL(KIND=JPRB), POINTER :: ZP1CHEM0(:,:), ZP1CHEM9(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDEFB11(:,:), ZPTENDEFB21(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDEFB31(:,:)
REAL(KIND=JPRB), POINTER :: ZP1EXT0(:,:,:), ZP1EXT9(:,:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDG1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDICONV1(:,:), ZPTENDI1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDLCONV1(:,:)
REAL(KIND=JPRB), POINTER :: ZP1LIMA0(:,:), ZP1LIMA9(:,:), ZPTENDL1(:,:)
REAL(KIND=JPRB), POINTER :: ZP1EZDIAG(:,:)
REAL(KIND=JPRB), POINTER :: ZP1NOGW0(:,:), ZP1NOGW9(:,:), ZP2NOGW0(:,:), ZP2NOGW9(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDQ1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDRCONV1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDR1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDSCONV1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDS1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDTKE1(:,:)

TYPE (APLPAR_STATE_TYPE) :: YLAPLPAR

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "apl_arome.intfb.h"
#include "aplpar.intfb.h"
#include "aplpar2intflex.intfb.h"
#include "aplpars.intfb.h"
#include "aplpassh.intfb.h"
#include "cpchet.intfb.h"
#include "cpmvvps.intfb.h"
#include "cpnudg.intfb.h"
#include "cpozo.intfb.h"
#include "cpphinp.intfb.h"
#include "cppsolan.intfb.h"
#include "cpqsol.intfb.h"
#include "cptend_new.intfb.h"
#include "cptend_flex.intfb.h"
#include "cptends.intfb.h"
#include "cptendsm.intfb.h"
#include "cputqy_arome.intfb.h"
#include "cputqy.intfb.h"
#include "cputqys.intfb.h"
#include "cpwts.intfb.h"
#include "cp_ptrslb1.intfb.h"
#include "gnhgw2svdarome.intfb.h"
#include "initaplpar.intfb.h"
#include "profilechet.intfb.h"
#include "rdradcoef.intfb.h"
#include "writephysio.intfb.h"
#include "wrphtrajm.intfb.h"
#include "wrradcoef.intfb.h"
#include "acajucv.intfb.h"
#include "gnhqe_conv_tempe.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MF_PHYS',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDVAB=>YDGEOMETRY%YRVAB, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL), YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL),  YDOROG=>YDGEOMETRY%YROROG(KIBL), &
 & YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,YDPTRSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1,YDPTRSLB2=>YDMODEL%YRML_DYN%YRPTRSLB2, &
 & YDTOPH=>YDMODEL%YRML_PHY_MF%YRTOPH,YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL,YDRIP=>YDMODEL%YRML_GCONF%YRRIP, &
 & YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH,YDRCOEF=>YDMODEL%YRML_PHY_RAD%YRRCOEF, &
 & YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY,YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY,YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, &
 & YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2,YGFL=>YDMODEL%YRML_GCONF%YGFL,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,  &
 & YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR, YDPRECIPS=>YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS, &
 & YDSTOPH=>YDMODEL%YRML_PHY_STOCH%YRSTOPH)

ASSOCIATE(MVTS=>YDPARAR%MVTS, CMF_UPDRAFT=>YDPARAR%CMF_UPDRAFT,&
 & TSPHY=>YDPHY2%TSPHY, &
 & NPROMA=>YDDIM%NPROMA, &
 & NTSSG=>YDDPHY%NTSSG, NVCLIS=>YDDPHY%NVCLIS, &
 & LMDUST=>YDARPHY%LMDUST, LMPA=>YDARPHY%LMPA, LMSE=>YDARPHY%LMSE, LMFSHAL=>YDARPHY%LMFSHAL,&
 & YI=>YGFL%YI, YH=>YGFL%YH, YEZDIAG=>YGFL%YEZDIAG, YL=>YGFL%YL, &
 & YEXT=>YGFL%YEXT, &
 & YUEN=>YGFL%YUEN, YG=>YGFL%YG, YMXL=>YGFL%YMXL, &
 & YQ=>YGFL%YQ, YR=>YGFL%YR, YSCONV=>YGFL%YSCONV, YS=>YGFL%YS, &
 & YEFB3=>YGFL%YEFB3, YEFB2=>YGFL%YEFB2, YEFB1=>YGFL%YEFB1, YRKTH=>YGFL%YRKTH, &
 & YFQTUR=>YGFL%YFQTUR, YFSTUR=>YGFL%YFSTUR, &
 & YUNEBH=>YGFL%YUNEBH, YO3=>YGFL%YO3, YNOGW=>YGFL%YNOGW, &
 & LCHEM_ARPCLIM=>YDMODEL%YRML_CHEM%YRCHEM%LCHEM_ARPCLIM,&
 & YCHEM=>YGFL%YCHEM, NGFL_EXT=>YGFL%NGFL_EXT, YRKTQC=>YGFL%YRKTQC, YTKE=>YGFL%YTKE, &
 & YDAL=>YGFL%YDAL, YCOMP=>YGFL%YCOMP, &
 & YRKTQV=>YGFL%YRKTQV, YUAL=>YGFL%YUAL, &
 & YTTE=>YGFL%YTTE, YLCONV=>YGFL%YLCONV, YSHTUR=>YGFL%YSHTUR, &
 & YRCONV=>YGFL%YRCONV, YICONV=>YGFL%YICONV, &
 & YSP_SBD=>YDSURF%YSP_SBD, &
 & YSD_VF=>YDSURF%YSD_VF, YSD_VH=>YDSURF%YSD_VH, &
 & YSD_VK=>YDSURF%YSD_VK, &
 & YSD_VV=>YDSURF%YSD_VV, &
 & YSD_VVD=>YDSURF%YSD_VVD, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, NFLSUL=>YDDIMV%NFLSUL, &
 & LTRAJPS=>YDSIMPHL%LTRAJPS, LSIMPH=>YDSIMPHL%LSIMPH, LRAYSP=>YDSIMPHL%LRAYSP, &
 & LNEBR=>YDPHY%LNEBR, LCDDPRO=>YDPHY%LCDDPRO, LNEBN=>YDPHY%LNEBN, &
 & LMPHYS=>YDPHY%LMPHYS, &
 & LCVPRO=>YDPHY%LCVPRO, &
 & LSTRAPRO=>YDPHY%LSTRAPRO, LPTKE=>YDPHY%LPTKE, &
 & NDPSFI=>YDPHY%NDPSFI, LOZONE=>YDPHY%LOZONE, L3MT=>YDPHY%L3MT, &
 & LGPCMT=>YDPHY%LGPCMT, LAJUCV=>YDPHY%LAJUCV, LCVPGY=>YDPHY%LCVPGY, &
 & LRRGUST=>YDPHY%LRRGUST, LEDR=>YDPHY%LEDR, &
 & NTAJUC=>YDTOPH%NTAJUC, NTPLUI=>YDTOPH%NTPLUI, &
 & NUMFLDS=>YGFL%NUMFLDS, LAGPHY=>YDEPHY%LAGPHY, &
 & LEPHYS=>YDEPHY%LEPHYS, &
 & LDPRECIPS=>YDPHY%LDPRECIPS, LDPRECIPS2=>YDPHY%LDPRECIPS2, &
 & NDTPRECCUR=>YDPRECIPS%NDTPRECCUR, NDTPRECCUR2=>YDPRECIPS%NDTPRECCUR2,&
 & LSDDH=>YDLDDH%LSDDH, &
 & HDSF=>YDMDDH%HDSF, &
 & MSLB1SP9=>YDPTRSLB1%MSLB1SP9, &
 & MSLB2KAPPAH=>YDPTRSLB2%MSLB2KAPPAH, &
 & MSLB2KAPPAM=>YDPTRSLB2%MSLB2KAPPAM, &
 & LRCOEF=>YDRCOEF%LRCOEF, &
 & LTLADDIA=>YDRCOEF%LTLADDIA, &
 & NG3SR=>YDRCOEF%NG3SR, &
 & NLIMA=>YGFL%NLIMA, YLIMA=>YGFL%YLIMA)

CALL SC2PRG(1,YCHEM(:)%MP,PGFL  ,ZP1CHEM0)    !  YDVARS%CHEM(1)%T0
CALL SC2PRG(1,YCHEM(:)%MP9,PGFL ,ZP1CHEM9)    !  YDVARS%CHEM(1)%T9
CALL SC2PRG(YEFB1%MP1 ,ZTENDGFL ,ZPTENDEFB11)
CALL SC2PRG(YEFB2%MP1 ,ZTENDGFL ,ZPTENDEFB21)
CALL SC2PRG(YEFB3%MP1 ,ZTENDGFL ,ZPTENDEFB31)
!ALL SC2PRG(1,YEXT(:)%MP,PGFL   ,ZP1EXT0)     !  YDVARS%EXT(1)%T0
!ALL SC2PRG(1,YEXT(:)%MP9,PGFL  ,ZP1EXT9)     !  YDVARS%EXT(1)%T9
CALL SC2PRG(1,YEXT(:)%MP ,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT,PGFL,ZP1EXT0) !  YDVARS%EXT(1)%T0
CALL SC2PRG(1,YEXT(:)%MP9,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT,PGFL,ZP1EXT9) !  YDVARS%EXT(1)%T9
CALL SC2PRG(1,YEZDIAG(:)%MP,PGFL,ZP1EZDIAG)
CALL SC2PRG(YG%MP1    ,ZTENDGFL ,ZPTENDG1)
CALL SC2PRG(YICONV%MP1,ZTENDGFL ,ZPTENDICONV1)
CALL SC2PRG(YI%MP1    ,ZTENDGFL ,ZPTENDI1)
CALL SC2PRG(YLCONV%MP1,ZTENDGFL ,ZPTENDLCONV1)
CALL SC2PRG(1,YLIMA(:)%MP,PGFL  ,ZP1LIMA0)    !  YDVARS%LIMA(1)%T0
CALL SC2PRG(1,YLIMA(:)%MP9,PGFL ,ZP1LIMA9)    !  YDVARS%LIMA(1)%T9
CALL SC2PRG(YL%MP1    ,ZTENDGFL ,ZPTENDL1)
CALL SC2PRG(1,YNOGW(:)%MP,PGFL  ,ZP1NOGW0)    !  YDVARS%NOGW(1)%T0
CALL SC2PRG(1,YNOGW(:)%MP9,PGFL ,ZP1NOGW9)    !  YDVARS%NOGW(1)%T9 
CALL SC2PRG(2,YNOGW(:)%MP,PGFL  ,ZP2NOGW0)    !  YDVARS%NOGW(2)%T0
CALL SC2PRG(2,YNOGW(:)%MP9,PGFL ,ZP2NOGW9)    !  YDVARS%NOGW(2)%T9 
CALL SC2PRG(YQ%MP1    ,ZTENDGFL ,ZPTENDQ1)
CALL SC2PRG(YRCONV%MP1,ZTENDGFL ,ZPTENDRCONV1)
CALL SC2PRG(YR%MP1    ,ZTENDGFL ,ZPTENDR1)
CALL SC2PRG(YSCONV%MP1,ZTENDGFL ,ZPTENDSCONV1)
CALL SC2PRG(YS%MP1    ,ZTENDGFL ,ZPTENDS1)
CALL SC2PRG(YTKE%MP1  ,ZTENDGFL ,ZPTENDTKE1)
CALL SC2PRG(YM_RADTC%MCOR,ZRADTC   ,ZMCOR)
CALL SC2PRG(YM_RADTC%MRAB3C,ZRADTC ,ZMRAB3C)
CALL SC2PRG(YM_RADTC%MRAB3N,ZRADTC ,ZMRAB3N)
CALL SC2PRG(YM_RADTC%MRAB4C,ZRADTC ,ZMRAB4C)
CALL SC2PRG(YM_RADTC%MRAB4N,ZRADTC ,ZMRAB4N)
CALL SC2PRG(YM_RADTC%MRAB6C,ZRADTC ,ZMRAB6C)
CALL SC2PRG(YM_RADTC%MRAB6N,ZRADTC ,ZMRAB6N)
CALL SC2PRG(YM_RADTC%MRAT1C,ZRADTC ,ZMRAT1C)
CALL SC2PRG(YM_RADTC%MRAT1N,ZRADTC ,ZMRAT1N)
CALL SC2PRG(YM_RADTC%MRAT2C,ZRADTC ,ZMRAT2C)
CALL SC2PRG(YM_RADTC%MRAT2N,ZRADTC ,ZMRAT2N)
CALL SC2PRG(YM_RADTC%MRAT3C,ZRADTC ,ZMRAT3C)
CALL SC2PRG(YM_RADTC%MRAT3N,ZRADTC ,ZMRAT3N)
CALL SC2PRG(YM_RADTC%MRAT4C,ZRADTC ,ZMRAT4C)
CALL SC2PRG(YM_RADTC%MRAT4N,ZRADTC ,ZMRAT4N)
CALL SC2PRG(YM_RADTC%MRAT5C,ZRADTC ,ZMRAT5C)
CALL SC2PRG(YM_RADTC%MRAT5N,ZRADTC ,ZMRAT5N)

CALL YLAPLPAR%INIT (LTWOTL, YDCPG_DYN0, YDCPG_DYN9, YDCPG_PHY0, YDCPG_PHY9, YDVARS, YDMF_PHYS_SURF, &
                  & ZP1EXT0, ZP1EXT9, ZP1CHEM0, ZP1CHEM9, ZP1NOGW0, ZP1NOGW9, ZP2NOGW0, ZP2NOGW9)

!     ------------------------------------------------------------------

!        0.    constructor for procset
IF (LINTFLEX) YLPROCSET=NEWINTPROCSET()

!        1.    Preliminary calculations necessary
!              for all types of physics.
!              ------------------------------------

IF (LSLAG) CALL CP_PTRSLB1(YDMODEL%YRML_DYN%YRDYN,YDPTRSLB1,ISLB1U9,ISLB1V9,ISLB1T9,ISLB1VD9,ISLB1GFL9)

IBLK=(KSTGLO-1)/NPROMA+1
INSTEP_DEB=1
INSTEP_FIN=1

! initialisation for surfex if XFU
LLXFUMSE=.FALSE.
IF (LDCONFX) THEN
  LLXFUMSE=.TRUE.
ENDIF

! SPP 
IF ( YSPP_CONFIG%LSPP ) THEN
 DO JROF=1,YSPP%N2D
   ZGP2DSPP(:,JROF) = YSPP%GP_ARP(JROF)%GP2D(:,1,KIBL)
 ENDDO
ENDIF

! Complete physics is called.
LLDIAB=(LMPHYS.OR.LEPHYS).AND.(.NOT.LAGPHY)

! In the NHQE model, MF_PHYS enters with Tt and grad(Tt), where Tt = T * exp(-(R/cp) log(pre/prehyd)).
! But calculations of MF_PHYS must use T and grad(T).
! So we do a conversion Tt -> T.
IF (LNHQE) THEN
  ! Valid for NPDVAR=2 only.
  ! At instant t (with the derivatives):
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZTT0_SAVE(JROF,JLEV)=YDVARS%T%T0(JROF,JLEV)
      ZTT0L_SAVE(JROF,JLEV)=YDVARS%T%DL(JROF,JLEV)
      ZTT0M_SAVE(JROF,JLEV)=YDVARS%T%DM(JROF,JLEV)
    ENDDO
  ENDDO
  CALL GNHQE_CONV_TEMPE(YDGEOMETRY,.TRUE.,YDMODEL%YRML_GCONF%YGFL%NDIM,KST,KEND,&
   & YDVARS%SPD%T0,YDVARS%SP%T0,YDVARS%T%T0,&
   & KGFLTYP=0,PGFL=PGFL,KDDER=2,PQCHAL=YDVARS%SPD%DL,PQCHAM=YDVARS%SPD%DM,&
   & PTL=YDVARS%T%DL,PTM=YDVARS%T%DM)
  ! At instant t-dt for leap-frog advections (without the derivatives):
  IF (.NOT.LTWOTL) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        ZTT9_SAVE(JROF,JLEV)=YDVARS%T%T9(JROF,JLEV)
      ENDDO
    ENDDO
    CALL GNHQE_CONV_TEMPE(YDGEOMETRY,.TRUE.,YDMODEL%YRML_GCONF%YGFL%NDIM,KST,KEND,&
     & YDVARS%SPD%T9,YDVARS%SP%T9,YDVARS%T%T9,&
     & KGFLTYP=9,PGFL=PGFL)
  ENDIF
ENDIF

IF (LRAYSP.AND.(NSTEP >= INSTEP_DEB .AND. NSTEP <= INSTEP_FIN)) THEN  
  CALL CPPSOLAN(YDGEOMETRY%YRDIM,KST,KEND,YDGSGEOM%GEMU,YDGSGEOM%GELAM,ZMMU0)
  IF (.NOT.LSDDH) THEN
!-----------INITIALIZING THE WEIGHT VECTORS-------------
DO JROF=KST,KEND
  YDCPG_MISC%DHSF(JROF)=HDSF(JROF+KSTGLO-1)
ENDDO
!-------------------------------------------------------
  ENDIF
ENDIF

IF (LLDIAB.OR.LSIMPH) THEN
  CALL CPPHINP(YDGEOMETRY,YDMODEL,KST,KEND,&
   & YDGSGEOM%GEMU,YDGSGEOM%GELAM,&
   & YDVARS%U%T0,YDVARS%V%T0,&
   & YDVARS%Q%T0,YDVARS%Q%DL,YDVARS%Q%DM,YDVARS%CVGQ%DL,YDVARS%CVGQ%DM,&
   & YDCPG_PHY0%XYB%RDELP,YDCPG_DYN0%CTY%EVEL,YDVARS%CVGQ%T0,&
   & ZMU0,ZMU0LU,ZMU0M,ZMU0N,ZCVGQ)
  ZLCVQ(KST:KEND,1:NFLEVG)=ZCVGQ(KST:KEND,1:NFLEVG)
ENDIF

DO JLEV=1,NFLEVG
  DO JROF=KST,KEND
    ZCVGQL(JROF,JLEV)=0._JPRB
    ZCVGQI(JROF,JLEV)=0._JPRB
    ZCVGT (JROF,JLEV)=0._JPRB
  ENDDO
ENDDO

DO JROF=KST,KEND
  ZQSATS(JROF)=0.0_JPRB
  ZFPLCH(JROF,0)=0.0_JPRB
  ZFPLSH(JROF,0)=0.0_JPRB
ENDDO

IF(LNEBN.OR.LNEBR.OR.LRRGUST) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZFPLCH(JROF,JLEV)=YDVARS%CPF%T0(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF
IF(LRRGUST) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZFPLSH(JROF,JLEV)=YDVARS%SPF%T0(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

! * In some cases, some pseudo-historic surface buffers (like z0) should
!   not be modified between the entrance and the output of MF_PHYS
!   (this is the case for example if LDCONFX=T).
!   For the time being, we must save:
!   - HV (group VV) : resistance to evapotranspiration
!   - Z0F (group VD): gravity * surface roughness length
!   - Z0H (group VV): gravity * roughness length for heat
!   - PBLH (group VH): PBL height
!   - SPSH (group VH):
!   - QSH (group VH):
LL_SAVE_PHSURF=LLDIAB.AND.LDCONFX
IF (LL_SAVE_PHSURF) THEN
  IF(YSD_VV%YHV%LSET) ZHV(1:NPROMA)=YDMF_PHYS_SURF%GSD_VV%PHV(1:NPROMA)
  IF(YSD_VF%YZ0F%LSET) ZGZ0F(1:NPROMA)=YDMF_PHYS_SURF%GSD_VF%PZ0F(1:NPROMA)
  IF(YSD_VV%YZ0H%LSET) ZGZ0HF(1:NPROMA)=YDMF_PHYS_SURF%GSD_VV%PZ0H(1:NPROMA)
  IF(YSD_VH%YPBLH%LSET) ZPBLH(1:NPROMA)=YDMF_PHYS_SURF%GSD_VH%PPBLH(1:NPROMA)
  IF(YSD_VH%YSPSH%LSET) ZFHPS(1:NPROMA)=YDMF_PHYS_SURF%GSD_VH%PSPSH(1:NPROMA)
  IF(YSD_VH%YQSH%LSET) ZQSH(1:NPROMA)=YDMF_PHYS_SURF%GSD_VH%PQSH(1:NPROMA)
  IF(YSD_VK%YUDGRO%LSET) ZUDGRO(1:NPROMA)=YDMF_PHYS_SURF%GSD_VK%PUDGRO(1:NPROMA)
  IF(LCVPRO.OR.LGPCMT) THEN
    ZUDAL(:,:)=YDVARS%UAL%T0(:,:)
    ZUDOM(:,:)=YDVARS%UOM%T0(:,:)
    IF(LCDDPRO) THEN
      ZDDAL(:,:)=YDVARS%DAL%T0(:,:)
      ZDDOM(:,:)=YDVARS%DOM%T0(:,:)
    ENDIF
  ENDIF
  IF(YUNEBH%LACTIVE) ZUNEBH(:,:)=YDVARS%UNEBH%T0(:,:)
  IF(YUEN%LACTIVE)   ZENTCH(:,:)=YDVARS%UEN%T0(:,:)
ENDIF

CALL INITAPLPAR ( YGFL, YDARPHY, KST, KEND, NPROMA, NFLEVG, NTSSG, YSP_SBD%NLEVS,&
  & YDMF_PHYS_SURF%GSD_VF%PVEG ,&
  & YDMF_PHYS%OUT%DIFCQ , YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS , ZDIFEXT, YDMF_PHYS%OUT%DIFTQ , YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL,&
  & YDMF_PHYS%OUT%DIFTS , YDMF_PHYS%OUT%FCCQL , YDMF_PHYS%OUT%FCCQN , YDMF_PHYS%OUT%FCSQL , YDMF_PHYS%OUT%FCSQN , YDMF_PHYS%OUT%FCQNG , YDMF_PHYS%OUT%FCQNNG,&
  & YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG,YDMF_PHYS%OUT%FCQGNG,&
  & YDMF_PHYS%OUT%FPLCL , YDMF_PHYS%OUT%FPLCN , YDMF_PHYS%OUT%FPLCG , YDMF_PHYS%OUT%FPLCH, YDMF_PHYS%OUT%FPLSL , YDMF_PHYS%OUT%FPLSN ,YDMF_PHYS%OUT%FPLSG ,&
  & YDMF_PHYS%OUT%FPLSH, YDMF_PHYS%OUT%FRSO  , YDMF_PHYS%OUT%FRSOC ,&
  & YDMF_PHYS%OUT%FRTH  , YDMF_PHYS%OUT%FRTHC , YDMF_PHYS%OUT%STRCU , YDMF_PHYS%OUT%STRCV , YDMF_PHYS%OUT%STRDU , YDMF_PHYS%OUT%STRDV , YDMF_PHYS%OUT%STRTU ,&
  & YDMF_PHYS%OUT%STRTV , YDMF_PHYS%OUT%STRMU , YDMF_PHYS%OUT%STRMV , YDMF_PHYS%OUT%FRMH  , ZFRMQ  ,&
  & YDMF_PHYS%OUT%DIFCQLC,YDMF_PHYS%OUT%DIFCQIC,YDMF_PHYS%OUT%FIMCC,&
  & YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC,YDMF_PHYS%OUT%FCNEGQIC,YDMF_PHYS%OUT%FCNEGQRC,YDMF_PHYS%OUT%FCNEGQSC,&
  & YDMF_PHYS%OUT%FCHOZ , ZCPS   , ZLHS   ,&
  & ZRS    , ZLH    , ZLSCPE , YDCPG_MISC%NEB   , YDCPG_MISC%QICE  , YDCPG_MISC%QLI   , ZQSAT  ,&
  & ZQW    , YDCPG_MISC%RH    , ZTW    , YDMF_PHYS%OUT%ALB , ZCD    , ZCDN   , ZCH    ,&
  & ZC1    , ZC2    , YDMF_PHYS%OUT%CT    , ZEMIS  , YDMF_PHYS%OUT%FCHSP , YDMF_PHYS%OUT%FCLL  , YDMF_PHYS%OUT%FCLN  ,&
  & YDMF_PHYS%OUT%FCS   , ZFEVI  , YDMF_PHYS%OUT%FEVL  , YDMF_PHYS%OUT%FEVN  , YDMF_PHYS%OUT%FEVV  , YDMF_PHYS%OUT%FLASH , YDMF_PHYS%OUT%FTR   , YDMF_PHYS%OUT%FLWSP ,&
  & YDMF_PHYS%OUT%FONTE , YDMF_PHYS%OUT%FGEL  , YDMF_PHYS%OUT%FGELS ,&
  & YDMF_PHYS%OUT%FRSGNI, YDMF_PHYS%OUT%FRSDNI, YDMF_PHYS%OUT%FRSODS, YDMF_PHYS%OUT%FRSOPS, YDMF_PHYS%OUT%FRSOPT, YDMF_PHYS%OUT%FRSOLU, YDMF_PHYS%OUT%FRTHDS,&
  & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN,YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN,YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPEVPSG,YDMF_PHYS%OUT%FPEVPCL,&
  & YDMF_PHYS%OUT%FPEVPCN,YDMF_PHYS%OUT%FPEVPCG,ZFTKE  , ZFTKEI , ZFEFB1 , ZFEFB2   , ZFEFB3,&
  & YDMF_PHYS%OUT%GZ0   , YDMF_PHYS%OUT%GZ0H  , ZNEIJ  , ZVEG   , YDCPG_MISC%QS    , ZQSATS , YDMF_PHYS%OUT%RUISL ,&
  & YDMF_PHYS%OUT%RUISP , YDMF_PHYS%OUT%RUISS , YDMF_PHYS%OUT%UCLS  , YDMF_PHYS%OUT%VCLS  , YDMF_PHYS%OUT%NUCLS , YDMF_PHYS%OUT%NVCLS , YDMF_PHYS%OUT%TCLS  , YDMF_PHYS%OUT%MRT,&
  & YDMF_PHYS%OUT%QCLS  , YDMF_PHYS%OUT%RHCLS , YDCPG_MISC%CLCT  , YDMF_PHYS%OUT%CLCH  , YDMF_PHYS%OUT%CLCM  , YDMF_PHYS%OUT%CLCL  , YDMF_PHYS%OUT%CLCC  ,&
  & YDMF_PHYS%OUT%CAPE  , YDMF_PHYS%OUT%CTOP  , ICLPH  , YDMF_PHYS%OUT%CLPH  , YDMF_PHYS%OUT%VEIN  , YDMF_PHYS%OUT%UGST  , YDMF_PHYS%OUT%VGST  ,&
  & YDMF_PHYS%OUT%DIAGH , ZEDR, YDMF_PHYS%OUT%VISICLD, YDMF_PHYS%OUT%VISIHYD, YDMF_PHYS%OUT%MXCLWC)

!*       2.    Complete physics.
!              -----------------

!        2.2  Complete physics.
!             -----------------

IF (LLDIAB.AND.(.NOT.LMPA)) THEN

  ! PAS DE TEMPS DE LA PHYSIQUE (/YOMPHY2/)
  ! Dans le cas des iterations de l'initialisation par modes normaux,
  ! le pas de temps pour la physique ne peux pas etre nul pour APLPAR et
  ! CPATY (par contre c'est bien PDTPHY qui est passe en argument aux autres
  ! sous-prog. de la physique). Ceci est du a l'impossibilite de prendre en
  ! compte des flux qui deviennent infinis pour TSPHY=0 (flux de masse du au
  ! reajustement des sursaturations par exemple...). Mais les tendances phys.
  ! sont bien nulles dans le cas de la configuration 'E' (Modes Normaux).
  ! PHYSICS TIME STEP (/YOMPHY2/)
  ! In case of normal mode initialisation iterations, the physics time
  ! step cannot be zero for APLPAR and CPATY (nevertheless it is PDTPHY
  ! which is passed to other physics subroutines). This is due to the
  ! impossibility to take into account fluxes which are infinite for TSPHY=0
  ! (e.g.: mass flux due to oversaturation...).

  TSPHY = MAX(PDTPHY,1.0_JPRB)

  ! CALL PARAMETERISATIONS

  DO JROF=KST,KEND
    YDMF_PHYS_SURF%GSD_VF%PLSM(JROF)=REAL(NINT(YDMF_PHYS_SURF%GSD_VF%PLSM(JROF)),JPRB)
  ENDDO
  
    IF (LTWOTL) THEN

      IF (LAJUCV) THEN
        DO JLEV=1,NFLEVG
          DO JROF=KST,KEND
            ZTAUX(JROF,JLEV)=YDVARS%T%T0(JROF,JLEV)
          ENDDO
        ENDDO
        CALL ACAJUCV(YDMODEL%YRML_PHY_MF%YRPHY0,KST,KEND,NPROMA,NTPLUI,NFLEVG,NTAJUC,&
         & YDCPG_PHY0%PREHYD,YDCPG_PHY0%XYB%ALPH,YDCPG_PHY0%XYB%DELP,&
         & YDCPG_PHY0%XYB%LNPR,YDVARS%T%T0)
        DO JLEV=1,NFLEVG
          DO JROF=KST,KEND
            ZDTAJU(JROF,JLEV)=YDVARS%T%T0(JROF,JLEV)-ZTAUX(JROF,JLEV)
          ENDDO
        ENDDO
      ENDIF

      ITDIA=1_JPIM
      CALL APLPAR(YDGEOMETRY,YDCPG_MISC,YDMF_PHYS,YDMF_PHYS_SURF,YDVARS,YDSURF, YDXFU, YDCFU, YDMODEL, KST    , KEND   , NPROMA ,&
       & ITDIA  , NFLEVG  , KSTGLO,&
       & NVCLIS , YSD_VVD%NUMFLDS ,&
       & NSTEP  ,&
       & NTSSG  , YSP_SBD%NLEVS   ,&
       & KBL    , KGPCOMP, YDCFU%NFRRC, PDTPHY,YDCSGEOM%RINDX,YDCSGEOM%RINDY, LLXFUMSE,&
       & YLAPLPAR%YCPG_DYN%PHI  , YLAPLPAR%YCPG_PHY%PREHYD  , YLAPLPAR%YCPG_DYN%PHIF , YLAPLPAR%YCPG_PHY%PREHYDF ,YLAPLPAR%YCPG_PHY%XYB%ALPH,&
       & YLAPLPAR%YCPG_PHY%XYB%DELP,&
       & YLAPLPAR%YCPG_PHY%XYB%LNPR,YLAPLPAR%YCPG_PHY%XYB%RDELP,&
       & YDGSGEOM%RCORI, YLAPLPAR%P1EXT,&
       & YLAPLPAR%U,YLAPLPAR%V,YLAPLPAR%T,&
       & YLAPLPAR%Q,YLAPLPAR%I,YLAPLPAR%L,&
       & YLAPLPAR%S,YLAPLPAR%R,YLAPLPAR%G,YLAPLPAR%TKE,&
       & YLAPLPAR%EFB1,YLAPLPAR%EFB2,YLAPLPAR%EFB3,&
       & YLAPLPAR%CVV,YLAPLPAR%O3,YLAPLPAR%P1CHEM,YLAPLPAR%P1NOGW,YLAPLPAR%P2NOGW, PGFL, &
       & YLAPLPAR%VOR,&
       & YLAPLPAR%YCPG_DYN%RCP%CP, ZCVGQ  ,YLAPLPAR%YCPG_DYN%RCP%R, PKOZO  , ZFPLCH , ZFPLSH ,&
       & YDCPG_DYN0%CTY%EVEL,&
       & PGPAR , &
       & YLAPLPAR%YGSP_SG%F , YLAPLPAR%YGSP_SG%A, YLAPLPAR%YGSP_SG%R,&
       & YLAPLPAR%YGSP_SB%T , YLAPLPAR%YGSP_RR%T , YLAPLPAR%YGSP_RR%FC ,&
       & YLAPLPAR%YGSP_SB%Q , YLAPLPAR%YGSP_SB%TL, YLAPLPAR%YGSP_RR%W  ,&
       & YLAPLPAR%YGSP_RR%IC,&
       & YDCPG_DYN0%CTY%VVEL(:,1:),&
       & ZMU0   , ZMU0LU , ZMU0M  ,ZMU0N,YDGSGEOM%GELAM,YDGSGEOM%GEMU,YDGSGEOM%GM,&
       & ZAC    , ZAC_HC , ZMCOR , ZMMU0  , &
       & ZMRAB3C,ZMRAB3N,&
       & ZMRAB4C,ZMRAB4N,&
       & ZMRAB6C,ZMRAB6N,&
       & ZMRAT1C,ZMRAT1N,&
       & ZMRAT2C,ZMRAT2N,&
       & ZMRAT3C,ZMRAT3N,&
       & ZMRAT4C,ZMRAT4N,&
       & ZMRAT5C,ZMRAT5N,&
       & YDOROG%OROG,&
       & YLAPLPAR%YCPG_PHY%W,YLAPLPAR%DIV,YDVARS%U%DL,YDVARS%V%DL,YLAPLPAR%YCPG_PHY%WL,YLAPLPAR%YCPG_PHY%WM,&
       & ZDIFEXT, &
       & ZFRMQ  ,&
       & ZCPS   , ZLHS   ,&
       & ZRS    , ZLH    , ZLSCPE , &
       & ZQRCONV, ZQSCONV,&
       & ZQSAT  ,&
       & ZQW    , ZTW    , ZCD    , ZCDN   ,&
       & ZCH    ,&
       & ZC1    , ZC2    , ZEMIS  , &
       & ZFEVI  , &
       & ZFTKE  , ZFTKEI,  ZFEFB1 , ZFEFB2 , ZFEFB3  ,&
       & ZNEIJ  , ZVEG   , &
       & ZQSATS , &
       & ICLPH  , &
       & ZDPRECIPS,ZDPRECIPS2,&
       & ZP1EZDIAG,ZTENDPTKE, ZKUROV_H, ZKTROV_H,&
       & ZTENDEXT_DEP, PTRAJ_PHYS,&
       & ZEDR, YDDDH, &
       & PFTCNS,&
       & ZGP2DSPP)

      IF (LAJUCV) THEN
        DO JLEV=1,NFLEVG
          DO JROF=KST,KEND
            YDVARS%T%T0(JROF,JLEV)=ZTAUX(JROF,JLEV)
          ENDDO
        ENDDO
        DO JLEV=1,NFLEVG
          DO JROF=KST,KEND
            PB1(JROF,ISLB1T9+JLEV-NFLSA)=PB1(JROF,ISLB1T9+JLEV-NFLSA)+ZDTAJU(JROF,JLEV)
          ENDDO
        ENDDO
      ENDIF

    ELSE

      ! IF (LAJUCV) THEN
      !   missing code under LAJUCV for leap-frog schemes.
      ! ENDIF

      ITDIA=1_JPIM
      CALL APLPAR(YDGEOMETRY,YDCPG_MISC,YDMF_PHYS,YDMF_PHYS_SURF,YDVARS,YDSURF,  YDXFU, YDCFU,YDMODEL,  KST,    KEND   , NPROMA ,&
       & ITDIA  , NFLEVG  , KSTGLO ,&
       & NVCLIS , YSD_VVD%NUMFLDS ,&
       & NSTEP  ,&
       & NTSSG  , YSP_SBD%NLEVS   ,&
       & KBL    , KGPCOMP, YDCFU%NFRRC, PDTPHY,YDCSGEOM%RINDX,YDCSGEOM%RINDY, LLXFUMSE,&
       & YLAPLPAR%YCPG_DYN%PHI  , YLAPLPAR%YCPG_PHY%PREHYD  , YLAPLPAR%YCPG_DYN%PHIF , YLAPLPAR%YCPG_PHY%PREHYDF , YLAPLPAR%YCPG_PHY%XYB%ALPH , &
       & YLAPLPAR%YCPG_PHY%XYB%DELP,&
       & YLAPLPAR%YCPG_PHY%XYB%LNPR,YLAPLPAR%YCPG_PHY%XYB%RDELP,&
       & YDGSGEOM%RCORI, YLAPLPAR%P1EXT,&
       & YLAPLPAR%U,YLAPLPAR%V,YLAPLPAR%T,&
       & YLAPLPAR%Q,YLAPLPAR%I,YLAPLPAR%L,&
       & YLAPLPAR%S,YLAPLPAR%R,YLAPLPAR%G,YLAPLPAR%TKE,&
       & YLAPLPAR%EFB1,YLAPLPAR%EFB2,YLAPLPAR%EFB3,&
       & YLAPLPAR%CVV,YLAPLPAR%O3,YLAPLPAR%P1CHEM,YLAPLPAR%P1NOGW,YLAPLPAR%P2NOGW, PGFL, &
       & YLAPLPAR%VOR,&
       & YLAPLPAR%YCPG_DYN%RCP%CP, ZCVGQ  ,YLAPLPAR%YCPG_DYN%RCP%R, PKOZO  , ZFPLCH, ZFPLSH  ,&
       & YDCPG_DYN0%CTY%EVEL,&
       & PGPAR , &
       & YLAPLPAR%YGSP_SG%F , YLAPLPAR%YGSP_SG%A, YLAPLPAR%YGSP_SG%R,&
       & YLAPLPAR%YGSP_SB%T , YLAPLPAR%YGSP_RR%T , YLAPLPAR%YGSP_RR%FC ,&
       & YLAPLPAR%YGSP_SB%Q , YLAPLPAR%YGSP_SB%TL, YLAPLPAR%YGSP_RR%W  ,&
       & YLAPLPAR%YGSP_RR%IC,&
       & YDCPG_DYN0%CTY%VVEL(:,1:),&
       & ZMU0   , ZMU0LU , ZMU0M  , ZMU0N,YDGSGEOM%GELAM,YDGSGEOM%GEMU,YDGSGEOM%GM,&
       & ZAC    , ZAC_HC , ZMCOR , ZMMU0  , &
       & ZMRAB3C,ZMRAB3N,&
       & ZMRAB4C,ZMRAB4N,&
       & ZMRAB6C,ZMRAB6N,&
       & ZMRAT1C,ZMRAT1N,&
       & ZMRAT2C,ZMRAT2N,&
       & ZMRAT3C,ZMRAT3N,&
       & ZMRAT4C,ZMRAT4N,&
       & ZMRAT5C,ZMRAT5N,&
       & YDOROG%OROG,&
       & YLAPLPAR%YCPG_PHY%W,YLAPLPAR%DIV,YDVARS%U%DL,YDVARS%V%DL,YLAPLPAR%YCPG_PHY%WL,YLAPLPAR%YCPG_PHY%WM,&
       & ZDIFEXT, &
       & ZFRMQ  ,&
       & ZCPS   , ZLHS   ,&
       & ZRS    , ZLH    , ZLSCPE , &
       & ZQRCONV, ZQSCONV,&
       & ZQSAT  ,&
       & ZQW    , ZTW    , ZCD    , ZCDN   ,&
       & ZCH    ,&
       & ZC1    , ZC2    , ZEMIS  , &
       & ZFEVI  , &
       & ZFTKE  , ZFTKEI,  ZFEFB1 , ZFEFB2 , ZFEFB3  ,&
       & ZNEIJ  , ZVEG   , &
       & ZQSATS , &
       & ICLPH  , &
       & ZDPRECIPS,ZDPRECIPS2,&
       & ZP1EZDIAG,ZTENDPTKE, ZKUROV_H, ZKTROV_H ,&
       & ZTENDEXT_DEP, PTRAJ_PHYS,&
       & ZEDR,YDDDH, &
       & PFTCNS,&
       & ZGP2DSPP)

      ! IF (LAJUCV) THEN
      !   missing code under LAJUCV for leap-frog schemes.
      ! ENDIF

    ENDIF

  !    convert to flexible interface structure
  IF (LINTFLEX) THEN
    CALL APLPAR2INTFLEX(YGFL,YDPHY,NPROMA, KST, KEND, NFLEVG,&
      & YDMF_PHYS%OUT%DIFCQ , YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS ,ZDIFEXT,&
      & YDMF_PHYS%OUT%DIFTQ , YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTS ,&
      & YDMF_PHYS%OUT%FCCQL , YDMF_PHYS%OUT%FCCQN , YDMF_PHYS%OUT%FCSQL , YDMF_PHYS%OUT%FCSQN ,&
      & YDMF_PHYS%OUT%FPLSL , YDMF_PHYS%OUT%FPLSN , YDMF_PHYS%OUT%FPLCL,  YDMF_PHYS%OUT%FPLCN,&
      & YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,&
      & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPCL ,YDMF_PHYS%OUT%FPFPCN ,&
      & YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG,&
      & YDMF_PHYS%OUT%FCQNG , YDMF_PHYS%OUT%FRMH  , ZFRMQ  , YDMF_PHYS%OUT%FRSO  , YDMF_PHYS%OUT%FRTH  ,&
      & YDMF_PHYS%OUT%STRCU , YDMF_PHYS%OUT%STRCV , YDMF_PHYS%OUT%STRDU , YDMF_PHYS%OUT%STRDV ,&
      & YDMF_PHYS%OUT%STRTU , YDMF_PHYS%OUT%STRTV , YDMF_PHYS%OUT%STRMU , YDMF_PHYS%OUT%STRMV ,&
      & YDMF_PHYS%OUT%DIFCQLC,YDMF_PHYS%OUT%DIFCQIC,YDMF_PHYS%OUT%FIMCC,&
      & YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC,YDMF_PHYS%OUT%FCNEGQIC,YDMF_PHYS%OUT%FCNEGQRC,YDMF_PHYS%OUT%FCNEGQSC,&
      & ZFTKE,&
      & ZTENDPTKE, ZTENDEXT, ZTENDEXT_DEP,&
      & YLPROCSET )
  ENDIF

ENDIF ! LLDIAB.AND..NOT.LMPA

!        2.3  Computes MOCON in the CLP.
!             --------------------------
YDMF_PHYS%OUT%MOCON(KST:KEND)=0.0_JPRB
IF (LLDIAB.AND.(.NOT.LMPA)) THEN

  IF (LTWOTL) THEN
    DO JLEV=1,NFLEVG
      DO JROF = KST, KEND
        YDMF_PHYS%OUT%MOCON(JROF) = YDMF_PHYS%OUT%MOCON(JROF)+(ZLCVQ(JROF,JLEV)-YDVARS%Q%T0(JROF,JLEV)*&
         & YDVARS%DIV%T0(JROF,JLEV))*YDCPG_PHY0%XYB%DELP(JROF,JLEV)&
         & *MAX(0,SIGN(1,JLEV-ICLPH(JROF)))  
      ENDDO
    ENDDO
    DO JROF = KST, KEND
      YDMF_PHYS%OUT%MOCON(JROF) = YDMF_PHYS%OUT%MOCON(JROF)/(YDCPG_PHY0%PREHYD(JROF,NFLEVG)-YDCPG_PHY0%PREHYD(JROF,ICLPH(JROF)-1))
    ENDDO
  ELSE
    DO JLEV=1,NFLEVG
      DO JROF = KST, KEND
        YDMF_PHYS%OUT%MOCON(JROF) = YDMF_PHYS%OUT%MOCON(JROF)+(ZLCVQ(JROF,JLEV)-YDVARS%Q%T9(JROF,JLEV)*&
         & YDVARS%DIV%T0(JROF,JLEV))*YDCPG_PHY9%XYB%DELP(JROF,JLEV)&
         & *MAX(0,SIGN(1,JLEV-ICLPH(JROF)))  
      ENDDO
    ENDDO
    DO JROF = KST, KEND
      YDMF_PHYS%OUT%MOCON(JROF) = YDMF_PHYS%OUT%MOCON(JROF)/(YDCPG_PHY9%PREHYD(JROF,NFLEVG)-YDCPG_PHY9%PREHYD(JROF,ICLPH(JROF)-1))
    ENDDO
  ENDIF
ENDIF

! Store surface water flux P and E for water conservation
IF (LCORWAT) THEN
  YDMF_PHYS_SURF%GSD_VH%PPSL(KST:KEND) = YDMF_PHYS%OUT%FPLSL(KST:KEND,NFLEVG)
  YDMF_PHYS_SURF%GSD_VH%PPCL(KST:KEND) = YDMF_PHYS%OUT%FPLCL(KST:KEND,NFLEVG)
  YDMF_PHYS_SURF%GSD_VH%PPSN(KST:KEND) = YDMF_PHYS%OUT%FPLSN(KST:KEND,NFLEVG)
  YDMF_PHYS_SURF%GSD_VH%PPCN(KST:KEND) = YDMF_PHYS%OUT%FPLCN(KST:KEND,NFLEVG)
  YDMF_PHYS_SURF%GSD_VH%PEVA(KST:KEND) = YDMF_PHYS%OUT%FEVN(KST:KEND,1)-YDMF_PHYS%OUT%FEVL(KST:KEND,1)
ENDIF

!        2.4  Stores radiation coefficients.
!             ------------------------------

! * writes grid-point transmission coefficients for simplified physics.

IF (LRCOEF.AND.(NSTEP == 1).AND.LLDIAB.AND.(.NOT.LMPA)) THEN
  IFIELDSS=NG3SR*NFLEVG
  CALL WRRADCOEF(YDGEOMETRY,YDRCOEF,KST,KEND,KSTGLO,IFIELDSS,ZRADTC,ZAC_HC)
ENDIF

!       2.5   Ozone
!             -----

IF (LLDIAB.AND.LOZONE.AND.(.NOT.LMPA)) THEN
  ! * Caution: this part has not been yet validated relative
  !   to the GFL implementation, and LOZONE (the setup of
  !   which has not yet been updated) can be true only if
  !   the GFL ozone is activated as a prognostic and advected
  !   variable.
  IPO3=(YO3%MP_SL1-1)*(NFLEVG+2*NFLSUL)
  IF (LSLAG) THEN
    IF (LTWOTL) THEN
      CALL CPOZO (NPROMA,KST,KEND,NFLEVG,PDTPHY,YDMF_PHYS%OUT%FCHOZ,&
       & PB1(1,ISLB1GFL9+IPO3+1-NFLSA),YDCPG_PHY0%XYB%RDELP)  
    ELSE
      CALL CPOZO (NPROMA,KST,KEND,NFLEVG,PDTPHY,YDMF_PHYS%OUT%FCHOZ,&
       & PB1(1,ISLB1GFL9+IPO3+1-NFLSA),YDCPG_PHY9%XYB%RDELP)  
    ENDIF
  ELSE
    CALL CPOZO (NPROMA,KST,KEND,NFLEVG,PDTPHY,YDMF_PHYS%OUT%FCHOZ,&
     & YDVARS%O3%T1,YDCPG_PHY9%XYB%RDELP)
  ENDIF
ENDIF

!        2.5.1 Chemical species   
!              ----------------
IF (LCHEM_ARPCLIM) THEN
   ! Processes described in my_phys ARPEGE-Climat 6.3 : to be added later here
   ! Modify also calls in CPTEND_NEW, etc.. as done ARPEGE-Climat 6.3.   
ENDIF   

!        2.6   surface specific humidity necessary to compute the vertical
!              advection of q in the case "delta m=1" (unlagged physics only).
!              ---------------------------------------------------------------

IF (LLDIAB.AND.(NDPSFI == 1)) THEN
  CALL CPQSOL(YDGEOMETRY%YRDIMV,YDPHY,NPROMA,KST,KEND,YDCPG_PHY0%PREHYD,YDMF_PHYS_SURF%GSP_RR%PT_T0,YDCPG_MISC%QS,ZQSATS,YDCPG_MISC%QSOL)
ENDIF

!        2.7  Computation of tendencies T,u,v and Q.
!             --------------------------------------

! Set GFL tendencies to 0

ZTENDGFL(:,:,:) = 0.0_JPRB

IF (LLDIAB.AND.(.NOT.LSIMPH).AND.(.NOT.LMPA)) THEN

  TSPHY = MAX(PDTPHY,1.0_JPRB)

  ! * CPTEND+CPUTQY = Old( CPATY + CPDUP + CPDTHP )
  ! Calcul des tendances de T , U et de Q et modifications
  ! eventuelles de W et de OMEGA/P

  IF (LTWOTL) THEN

    IF (LINTFLEX.AND.(.NOT.LDCONFX)) THEN
      CALL CPTEND_FLEX( YDLDDH,YDMDDH,YGFL,YDPHY,NPROMA, KST, KEND, NFLEVG,YDGSGEOM%GNORDL,YDGSGEOM%GNORDM,&
       & YDCPG_PHY0%XYB%DELP ,&
       & YDCPG_PHY0%XYB%RDELP, YDCPG_DYN0%RCP%CP,&
       & YDVARS%U%T0,YDVARS%V%T0,YDVARS%T%T0,YDMF_PHYS_SURF%GSP_RR%PT_T0,&
       & PGFL,&
       & YLPROCSET,&
       & YDMF_PHYS%OUT%TENDU , YDMF_PHYS%OUT%TENDV , ZTENDH , ZTENDGFL,&
       & YDMF_PHYS%OUT%FHSCL ,YDMF_PHYS%OUT%FHSCN,YDMF_PHYS%OUT%FHSSL,YDMF_PHYS%OUT%FHSSN,&
       & YDMF_PHYS%OUT%FHPCL ,YDMF_PHYS%OUT%FHPCN,YDMF_PHYS%OUT%FHPSL,YDMF_PHYS%OUT%FHPSN,&
       & ZFHP   ,ZFP   ,  YDMF_PHYS%OUT%FEPFP, YDMF_PHYS%OUT%FCMPCQ, YDMF_PHYS%OUT%FCMPSN, YDMF_PHYS%OUT%FCMPSL,YDDDH )
    ELSE
      CALL CPTEND_NEW( YDMODEL, NPROMA, KST, KEND, NFLEVG,YDGSGEOM%GNORDL,YDGSGEOM%GNORDM,&
       & YDMF_PHYS%OUT%DIFCQ , YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS ,ZDIFEXT,&
       & YDMF_PHYS%OUT%DIFTQ , YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTS ,&
       & YDMF_PHYS%OUT%FCCQL , YDMF_PHYS%OUT%FCCQN , YDMF_PHYS%OUT%FCSQL , YDMF_PHYS%OUT%FCSQN ,&
       & YDMF_PHYS%OUT%FPLSL , YDMF_PHYS%OUT%FPLSN , YDMF_PHYS%OUT%FPLSG , YDMF_PHYS%OUT%FPLCL,  YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLCG,&
       & YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPEVPSG,YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,YDMF_PHYS%OUT%FPEVPCG,&
       & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPFPCL ,YDMF_PHYS%OUT%FPFPCN ,&
       & YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FCQGNG ,&
       & YDMF_PHYS%OUT%FCQNG , YDMF_PHYS%OUT%FRMH  , ZFRMQ  , YDMF_PHYS%OUT%FRSO  , YDMF_PHYS%OUT%FRTH  ,&
       & YDMF_PHYS%OUT%STRCU , YDMF_PHYS%OUT%STRCV , YDMF_PHYS%OUT%STRDU , YDMF_PHYS%OUT%STRDV ,&
       & YDMF_PHYS%OUT%STRTU , YDMF_PHYS%OUT%STRTV , YDMF_PHYS%OUT%STRMU , YDMF_PHYS%OUT%STRMV ,&
       & YDMF_PHYS%OUT%DIFCQLC,YDMF_PHYS%OUT%DIFCQIC,YDMF_PHYS%OUT%FIMCC,&
       & YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC,YDMF_PHYS%OUT%FCNEGQIC,YDMF_PHYS%OUT%FCNEGQRC,YDMF_PHYS%OUT%FCNEGQSC,&
       & ZFTKE  , ZFTKEI,  ZFEFB1 , ZFEFB2 , ZFEFB3 ,YDCPG_PHY0%XYB%DELP ,&
       & YDCPG_PHY0%XYB%RDELP, YDCPG_DYN0%PHIF  , YDCPG_DYN0%RCP%CP,&
       & YDVARS%U%T0,YDVARS%V%T0,YDVARS%T%T0,&
       & YDVARS%Q%T0,YDVARS%I%T0,YDVARS%L%T0,&
       & YDVARS%LCONV%T0,YDVARS%ICONV%T0,YDVARS%RCONV%T0,YDVARS%SCONV%T0,&
       & YDVARS%R%T0,YDVARS%S%T0,YDVARS%G%T0,&
       & ZCPS   , YDMF_PHYS_SURF%GSP_RR%PT_T0  ,&
       & YDMF_PHYS%OUT%FHSCL ,YDMF_PHYS%OUT%FHSCN,YDMF_PHYS%OUT%FHSSL,YDMF_PHYS%OUT%FHSSN,YDMF_PHYS%OUT%FHSSG,&
       & YDMF_PHYS%OUT%FHPCL ,YDMF_PHYS%OUT%FHPCN,YDMF_PHYS%OUT%FHPCG,YDMF_PHYS%OUT%FHPSL,YDMF_PHYS%OUT%FHPSN,YDMF_PHYS%OUT%FHPSG,&
       & ZFHP   ,ZFP   ,  YDMF_PHYS%OUT%FEPFP, YDMF_PHYS%OUT%FCMPCQ, YDMF_PHYS%OUT%FCMPSN, YDMF_PHYS%OUT%FCMPSL,&
       & YDMF_PHYS%OUT%TENDU , YDMF_PHYS%OUT%TENDV , ZTENDU, ZTENDV, ZTENDH ,&
       & ZPTENDQ1,ZPTENDI1,ZPTENDL1,&
       & ZPTENDLCONV1,ZPTENDICONV1,&
       & ZPTENDRCONV1,ZPTENDSCONV1,&
       & ZPTENDR1,ZPTENDS1,ZPTENDG1,ZPTENDTKE1,&
       & ZPTENDEFB11,ZPTENDEFB21,ZPTENDEFB31,&
       & ZTENDEXT,YDDDH)
    ENDIF
  ELSE
    IF (LINTFLEX.AND.(.NOT.LDCONFX)) THEN
      CALL CPTEND_FLEX( YDLDDH,YDMDDH,YGFL,YDPHY,NPROMA, KST, KEND, NFLEVG,YDGSGEOM%GNORDL,YDGSGEOM%GNORDM,&
       & YDCPG_PHY9%XYB%DELP ,&
       & YDCPG_PHY9%XYB%RDELP, YDCPG_DYN9%RCP%CP,&
       & YDVARS%U%T9,YDVARS%V%T9,YDVARS%T%T9,YDMF_PHYS_SURF%GSP_RR%PT_T9,&
       & PGFL,&
       & YLPROCSET,&
       & YDMF_PHYS%OUT%TENDU , YDMF_PHYS%OUT%TENDV , ZTENDH , ZTENDGFL,&
       & YDMF_PHYS%OUT%FHSCL ,YDMF_PHYS%OUT%FHSCN,YDMF_PHYS%OUT%FHSSL,YDMF_PHYS%OUT%FHSSN,&
       & YDMF_PHYS%OUT%FHPCL ,YDMF_PHYS%OUT%FHPCN,YDMF_PHYS%OUT%FHPSL,YDMF_PHYS%OUT%FHPSN,&
       & ZFHP   ,ZFP   ,  YDMF_PHYS%OUT%FEPFP, YDMF_PHYS%OUT%FCMPCQ, YDMF_PHYS%OUT%FCMPSN, YDMF_PHYS%OUT%FCMPSL,YDDDH )
    ELSE
      CALL CPTEND_NEW( YDMODEL, NPROMA, KST, KEND, NFLEVG,YDGSGEOM%GNORDL,YDGSGEOM%GNORDM,&
       & YDMF_PHYS%OUT%DIFCQ , YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS ,ZDIFEXT,&
       & YDMF_PHYS%OUT%DIFTQ , YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTS ,&
       & YDMF_PHYS%OUT%FCCQL , YDMF_PHYS%OUT%FCCQN , YDMF_PHYS%OUT%FCSQL , YDMF_PHYS%OUT%FCSQN ,&
       & YDMF_PHYS%OUT%FPLSL , YDMF_PHYS%OUT%FPLSN , YDMF_PHYS%OUT%FPLSG , YDMF_PHYS%OUT%FPLCL,  YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLCG,&
      & YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPEVPSG,YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,YDMF_PHYS%OUT%FPEVPCG,&
       & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPFPCL ,YDMF_PHYS%OUT%FPFPCN,&
       & YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FCQGNG,&
       & YDMF_PHYS%OUT%FCQNG , YDMF_PHYS%OUT%FRMH  , ZFRMQ  , YDMF_PHYS%OUT%FRSO  , YDMF_PHYS%OUT%FRTH  ,&
       & YDMF_PHYS%OUT%STRCU , YDMF_PHYS%OUT%STRCV , YDMF_PHYS%OUT%STRDU , YDMF_PHYS%OUT%STRDV ,&
       & YDMF_PHYS%OUT%STRTU , YDMF_PHYS%OUT%STRTV , YDMF_PHYS%OUT%STRMU , YDMF_PHYS%OUT%STRMV ,&
       & YDMF_PHYS%OUT%DIFCQLC,YDMF_PHYS%OUT%DIFCQIC,YDMF_PHYS%OUT%FIMCC,&
       & YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC,YDMF_PHYS%OUT%FCNEGQIC,YDMF_PHYS%OUT%FCNEGQRC,YDMF_PHYS%OUT%FCNEGQSC,&
       & ZFTKE  , ZFTKEI,  ZFEFB1 , ZFEFB2 , ZFEFB3 , YDCPG_PHY9%XYB%DELP ,&
       & YDCPG_PHY9%XYB%RDELP, YDCPG_DYN9%PHIF  , YDCPG_DYN9%RCP%CP,&
       & YDVARS%U%T9,YDVARS%V%T9,YDVARS%T%T9,&
       & YDVARS%Q%T9,YDVARS%I%T9,YDVARS%L%T9,&
       & YDVARS%LCONV%T0,YDVARS%ICONV%T0,YDVARS%RCONV%T0,YDVARS%SCONV%T0,&
       & YDVARS%R%T9,YDVARS%S%T9,YDVARS%G%T9,&
       & ZCPS   , YDMF_PHYS_SURF%GSP_RR%PT_T9,&
       & YDMF_PHYS%OUT%FHSCL ,YDMF_PHYS%OUT%FHSCN,YDMF_PHYS%OUT%FHSSL,YDMF_PHYS%OUT%FHSSN,YDMF_PHYS%OUT%FHSSG,&
       & YDMF_PHYS%OUT%FHPCL ,YDMF_PHYS%OUT%FHPCN,YDMF_PHYS%OUT%FHPCG,YDMF_PHYS%OUT%FHPSL,YDMF_PHYS%OUT%FHPSN,YDMF_PHYS%OUT%FHPSG,&
       & ZFHP   ,ZFP   ,  YDMF_PHYS%OUT%FEPFP, YDMF_PHYS%OUT%FCMPCQ, YDMF_PHYS%OUT%FCMPSN, YDMF_PHYS%OUT%FCMPSL,&
       & YDMF_PHYS%OUT%TENDU , YDMF_PHYS%OUT%TENDV , ZTENDU, ZTENDV, ZTENDH ,&
       & ZPTENDQ1,ZPTENDI1,ZPTENDL1,&
       & ZPTENDLCONV1,ZPTENDICONV1,&
       & ZPTENDRCONV1,ZPTENDSCONV1,&
       & ZPTENDR1,ZPTENDS1,ZPTENDG1,ZPTENDTKE1,&
       & ZPTENDEFB11,ZPTENDEFB21,ZPTENDEFB31,&
       & ZTENDEXT,YDDDH)
    ENDIF
    
    IF ( L3MT.OR.LSTRAPRO.OR.(NDPSFI==1)) THEN
!     PFEPFP was ZFEPFP in CPTEND_NEW, before, ZFEPFP still in CPFHPFS
      DO JLEV= 0, NFLEVG 
        DO JROF = 1, NPROMA
          YDMF_PHYS%OUT%FEPFP(JROF,JLEV) = 0.0_JPRB
          YDMF_PHYS%OUT%FCMPCQ(JROF,JLEV) = 0.0_JPRB
          YDMF_PHYS%OUT%FCMPSN(JROF,JLEV) = 0.0_JPRB
          YDMF_PHYS%OUT%FCMPSL(JROF,JLEV) = 0.0_JPRB
        ENDDO
      ENDDO
    ENDIF

  ENDIF ! LTWOTL



!        2.7.1  Diagnostics on physical tendencies
!               ----------------------------------

  IF (.NOT.LDCONFX) THEN
    IF ((GCHETN%LFREQD).OR.(GCHETN%LCOORD).OR.(GCHETN%LPROFV)) THEN
      IF (LTWOTL) THEN
        CALL CPCHET(  YDRIP,YDPHY, NPROMA, KST, KEND, NFLEVG, NSTEP,&
        & YDMF_PHYS%OUT%FHSCL , YDMF_PHYS%OUT%FHSCN , YDMF_PHYS%OUT%FHSSL , YDMF_PHYS%OUT%FHSSN ,&
        & YDMF_PHYS%OUT%FHPCL , YDMF_PHYS%OUT%FHPCN , YDMF_PHYS%OUT%FHPSL , YDMF_PHYS%OUT%FHPSN ,&
        & YDMF_PHYS%OUT%DIFCQ , YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS ,&
        & YDMF_PHYS%OUT%DIFTQ , YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTS ,&
        & YDMF_PHYS%OUT%FCCQL , YDMF_PHYS%OUT%FCCQN , YDMF_PHYS%OUT%FCSQL , YDMF_PHYS%OUT%FCSQN ,&
        & YDMF_PHYS%OUT%FPLCL , YDMF_PHYS%OUT%FPLCN , YDMF_PHYS%OUT%FPLSL , YDMF_PHYS%OUT%FPLSN ,&
        & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN,&
        & YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,&
        & YDMF_PHYS%OUT%FRMH  , ZFRMQ  , YDMF_PHYS%OUT%FRSO  , YDMF_PHYS%OUT%FRTH  ,&
        & YDMF_PHYS%OUT%STRCU , YDMF_PHYS%OUT%STRCV , YDMF_PHYS%OUT%STRDU , YDMF_PHYS%OUT%STRDV ,&
        & YDMF_PHYS%OUT%STRTU , YDMF_PHYS%OUT%STRTV , YDMF_PHYS%OUT%STRMU , YDMF_PHYS%OUT%STRMV ,&
        & YDCPG_PHY0%XYB%RDELP, YDCPG_DYN0%RCP%CP,&
        & YDVARS%T%T0,YDVARS%Q%T0,YDVARS%I%T0,YDVARS%L%T0,&
        & ZCPS   , YDMF_PHYS_SURF%GSP_RR%PT_T0 , YDCPG_MISC%QS  ,&
        & YDMF_PHYS%OUT%TENDU , YDMF_PHYS%OUT%TENDV , ZTENDH ,&
        & ZPTENDQ1,ZPTENDI1,ZPTENDL1,&
        & ZPTENDR1,ZPTENDS1,&
        & YDGSGEOM%GEMU,YDGSGEOM%GELAM, YDCPG_DYN0%PHI   , YDCPG_DYN0%PHIF  )
      ELSE
        CALL CPCHET(  YDRIP,YDPHY, NPROMA, KST, KEND, NFLEVG, NSTEP,&
        & YDMF_PHYS%OUT%FHSCL , YDMF_PHYS%OUT%FHSCN , YDMF_PHYS%OUT%FHSSL , YDMF_PHYS%OUT%FHSSN ,&
        & YDMF_PHYS%OUT%FHPCL , YDMF_PHYS%OUT%FHPCN , YDMF_PHYS%OUT%FHPSL , YDMF_PHYS%OUT%FHPSN ,&
        & YDMF_PHYS%OUT%DIFCQ , YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS ,&
        & YDMF_PHYS%OUT%DIFTQ , YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTS ,&
        & YDMF_PHYS%OUT%FCCQL , YDMF_PHYS%OUT%FCCQN , YDMF_PHYS%OUT%FCSQL , YDMF_PHYS%OUT%FCSQN ,&
        & YDMF_PHYS%OUT%FPLCL , YDMF_PHYS%OUT%FPLCN , YDMF_PHYS%OUT%FPLSL , YDMF_PHYS%OUT%FPLSN ,&
        & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN,&
        & YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,&
        & YDMF_PHYS%OUT%FRMH  , ZFRMQ  , YDMF_PHYS%OUT%FRSO  , YDMF_PHYS%OUT%FRTH  ,&
        & YDMF_PHYS%OUT%STRCU , YDMF_PHYS%OUT%STRCV , YDMF_PHYS%OUT%STRDU , YDMF_PHYS%OUT%STRDV ,&
        & YDMF_PHYS%OUT%STRTU , YDMF_PHYS%OUT%STRTV , YDMF_PHYS%OUT%STRMU , YDMF_PHYS%OUT%STRMV ,&
        & YDCPG_PHY9%XYB%RDELP, YDCPG_DYN9%RCP%CP,&
        & YDVARS%T%T9,YDVARS%Q%T9,YDVARS%I%T9,YDVARS%L%T9,&
        & ZCPS   , YDMF_PHYS_SURF%GSP_RR%PT_T9  , YDCPG_MISC%QS    ,&
        & YDMF_PHYS%OUT%TENDU , YDMF_PHYS%OUT%TENDV , ZTENDH ,&
        & ZPTENDQ1,ZPTENDI1,ZPTENDL1,&
        & ZPTENDR1,ZPTENDS1,&
        & YDGSGEOM%GEMU,YDGSGEOM%GELAM, YDCPG_DYN9%PHI   , YDCPG_DYN9%PHIF  )
      ENDIF
    ENDIF
    
    IF (GCHETN%LPROFV)&
     & CALL PROFILECHET(YDGEOMETRY,YDSURF,YDDPHY,YDRIP,YDMODEL%YRML_PHY_MF,&
     & KEND,&
     & YDGSGEOM%GELAM,YDGSGEOM%GEMU,&
     & YDGSGEOM%GM,ZMU0,YDOROG%OROG,YDCPG_DYN0%OROGL,YDCPG_DYN0%OROGM,YDGSGEOM%RCORI,YDCSGEOM%RATATH,YDCSGEOM%RATATX,&
     & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PARG, YDMF_PHYS_SURF%GSD_VV%PSAB, YDMF_PHYS_SURF%GSD_VF%PALBF,&
     & YDMF_PHYS_SURF%GSD_VF%PALBSF, YDMF_PHYS_SURF%GSD_VF%PEMISF, YDMF_PHYS_SURF%GSD_VV%PD2, YDMF_PHYS_SURF%GSD_VV%PIVEG,&
     & YDMF_PHYS_SURF%GSD_VV%PLAI, YDMF_PHYS%OUT%CT, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H, YDMF_PHYS_SURF%GSD_VF%PZ0RLF,&
     & YDMF_PHYS_SURF%GSD_VF%PGETRL, YDMF_PHYS_SURF%GSD_VF%PVRLAN, YDMF_PHYS_SURF%GSD_VF%PVRLDI, YDMF_PHYS_SURF%GSD_VV%PRSMIN,&
     & YDMF_PHYS_SURF%GSD_VF%PVEG, YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VA%PSEA, YDMF_PHYS_SURF%GSD_VA%PLAN,&
     & YDMF_PHYS_SURF%GSD_VA%PSOO, YDMF_PHYS_SURF%GSD_VA%PDES, YDMF_PHYS_SURF%GSD_VC%PGROUP,&
     & YDVARS%SP%T0,YDVARS%SP%DL,YDVARS%SP%DM,&
     & YDVARS%T%T0,YDVARS%T%DL,YDVARS%T%DM,&
     & YDVARS%Q%T0,YDVARS%Q%DL,YDVARS%Q%DM,&
     & YDVARS%U%T0,YDVARS%V%T0,YDVARS%VOR%T0,YDVARS%DIV%T0, ZCVGQ, ZLCVQ,&
     & YDMF_PHYS_SURF%GSP_RR%PT_T9, YDMF_PHYS_SURF%GSP_SB%PT_T9, YDMF_PHYS_SURF%GSP_RR%PFC_T9, YDMF_PHYS_SURF%GSP_RR%PW_T9,&
     & YDMF_PHYS_SURF%GSP_RR%PIC_T9, YDMF_PHYS_SURF%GSP_SB%PQ_T9, YDMF_PHYS_SURF%GSP_SB%PTL_T9, YDMF_PHYS_SURF%GSP_SG%PF_T9,&
     & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL,&
     & YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%FCQNNG,&
     & YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH,&
     & YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, YDMF_PHYS%OUT%STRDU, YDMF_PHYS%OUT%STRDV, YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV, YDMF_PHYS%OUT%STRMU, YDMF_PHYS%OUT%STRMV,&
     & YDMF_PHYS%OUT%FRMH, YDMF_PHYS%OUT%FCHOZ, YDCPG_MISC%NEB, YDCPG_MISC%QICE, YDCPG_MISC%QLI, YDCPG_MISC%RH,&
     & YDMF_PHYS%OUT%FCS, YDMF_PHYS%OUT%FCLL, YDMF_PHYS%OUT%FCLN, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FTR, YDMF_PHYS%OUT%FLWSP, YDMF_PHYS%OUT%FONTE,&
     & YDMF_PHYS%OUT%FGEL, YDMF_PHYS%OUT%FGELS, YDMF_PHYS%OUT%FRSODS, YDMF_PHYS%OUT%FRSOPS, YDMF_PHYS%OUT%FRSOPT, YDMF_PHYS%OUT%FRTHDS,&
     & YDCPG_MISC%QS, YDMF_PHYS%OUT%RUISL, YDMF_PHYS%OUT%RUISP, YDMF_PHYS%OUT%RUISS,&
     & YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%QCLS, YDMF_PHYS%OUT%RHCLS,&
     & YDCPG_MISC%CLCT, YDMF_PHYS%OUT%CLCH, YDMF_PHYS%OUT%CLCM, YDMF_PHYS%OUT%CLCL, YDMF_PHYS%OUT%CLCC,&
     & ZFPLCH, ZFPLSH,&
     & YDMF_PHYS_SURF%GSD_VH%PTCCH,  YDMF_PHYS_SURF%GSD_VH%PSCCH, YDMF_PHYS_SURF%GSD_VH%PBCCH, YDMF_PHYS_SURF%GSD_VH%PPBLH )

  ENDIF

ENDIF

!        2.8  Modification of vertical velocities
!             by some physics output when required.
!             -------------------------------------

IF (LLDIAB.AND.(.NOT.LSIMPH)) THEN

  ! * MODIFICATION DE LA VITESSE VERTICALE ET DE LA TENDANCE DE
  ! PRESSION DE SURFACE SI NDPSFI=1 ( MASSE VARIABLE ).
  ! Ajout de la physique dans l'equation de continuite/Add physics
  ! in continuity equation.

  IF (NDPSFI == 1) THEN
    IF (LSLAG .AND. LTWOTL) THEN
      CALL CPMVVPS(YDVAB,NPROMA,KST,KEND,NFLEVG,PDTPHY,&
       & ZFP,YDCPG_PHY0%PREHYD(:,NFLEVG),YDMF_PHYS%OUT%FEVL,YDMF_PHYS%OUT%FEVN,&
       & YDCPG_DYN0%CTY%EVEL,YDCPG_DYN0%CTY%PSDVBC,PB1(1,MSLB1SP9))
    ELSEIF (LSLAG .AND. (.NOT.LTWOTL)) THEN
      CALL CPMVVPS(YDVAB,NPROMA,KST,KEND,NFLEVG,PDTPHY,&
       & ZFP,YDCPG_PHY9%PREHYD(:,NFLEVG),YDMF_PHYS%OUT%FEVL,YDMF_PHYS%OUT%FEVN,&
       & YDCPG_DYN0%CTY%EVEL,YDCPG_DYN0%CTY%PSDVBC,PB1(1,MSLB1SP9))
    ELSE
      CALL CPMVVPS(YDVAB,NPROMA,KST,KEND,NFLEVG,PDTPHY,&
       & ZFP,YDCPG_PHY9%PREHYD(:,NFLEVG),YDMF_PHYS%OUT%FEVL,YDMF_PHYS%OUT%FEVN,&
       & YDCPG_DYN0%CTY%EVEL,YDCPG_DYN0%CTY%PSDVBC,YDVARS%SP%T1)
    ENDIF
  ENDIF

ENDIF

!        2.9  Computation of evolution of T, u, v and Q.
!             ------------------------------------------

! * Calculation of IPGFL, since the old pointers
!   MSLB1[X]9 (=MSLB1GFL9+IP[X]) do not exist any longer in PTRSLB1.

! usefull pointer for new version of cputqy

DO JGFL=1,NUMFLDS
  IF ((YCOMP(JGFL)%MP1 > 0) .AND. (YCOMP(JGFL)%MP_SL1 > 0)) THEN
     IPGFL(YCOMP(JGFL)%MP1) = (YCOMP(JGFL)%MP_SL1-1)*(NFLEVG+2*NFLSUL)
  ENDIF   
ENDDO  

!  ALARO does not respect the coding rules, tendency of pseudo-TKE is computed in APLPAR and not
!  in CPTEND_NEW. To use the new version of cputqy it is then necessary to write it in GFL tendencies array.
! This memory transfer is not necessary, please respect coding rules to avoid it.

! Not necessary for intflex: already done in aplpar2intflex
IF (.NOT.(LINTFLEX.AND.(.NOT.LDCONFX))) THEN
  IF (LPTKE) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        ZPTENDTKE1(JROF,JLEV) = ZTENDPTKE(JROF,JLEV)
      ENDDO
    ENDDO    
  ENDIF
  ! Extra-GFL
  IF(LMDUST.AND.(NGFL_EXT/=0)) THEN
    DO JGFL=1, NGFL_EXT
      DO JLEV=1,NFLEVG
        DO JROF=KST,KEND
          ZTENDGFL(JROF,JLEV,YEXT(JGFL)%MP1) = ZTENDEXT(JROF,JLEV,JGFL)+&! turbulent tendency
                                             & ZTENDEXT_DEP(JROF,JLEV,JGFL) ! moist tendency
! moist tendency
! moist tendency
! moist tendency
! moist tendency
! moist tendency
! moist tendency
! moist tendency
! moist tendency
        ENDDO
      ENDDO 
    ENDDO   
  ENDIF
ENDIF

! ky: non-zero option not yet coded for the time being.
ZTENDD=0.0_JPRB

IF (LLDIAB.AND.(.NOT.LSIMPH).AND.(.NOT.LMPA)) THEN
  ! Calcul de T , Q et du Vent a l'instant 1

  CALL CPUTQY(YDGEOMETRY%YRDIMV,YDGMV,YGFL,YDPTRSLB1,YDPHY,NPROMA,KST,KEND,NFLEVG,PDTPHY,IPGFL,&
   & ISLB1T9,ISLB1U9,ISLB1V9,ISLB1VD9,ISLB1GFL9,&
   & ZTENDH, ZTENDT, YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDU, ZTENDV, ZTENDD, ZTENDGFL,&
   & YDCPG_DYN0%RCP%CP,YDCPG_PHY0%XYB%DELP,YDVARS%T%T0,YDVARS%U%T0,YDVARS%V%T0,&
   & YDCPG_DYN9%RCP%CP,YDCPG_PHY9%XYB%DELP,YDVARS%T%T9,YDVARS%U%T9,YDVARS%V%T9,&
   & PB1, PGMVT1, PGFLT1,&
   & YDMF_PHYS%OUT%FDIS)

ENDIF

!        2.9a Evolution of precipitation fluxes
!             ------------------------------------------

IF(LNEBN.OR.LNEBR.OR.LRRGUST) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      YDVARS%CPF%T1(JROF,JLEV)=ZFPLCH(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF
IF(LRRGUST) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      YDVARS%SPF%T1(JROF,JLEV)=ZFPLSH(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

!       2.9b Prognostic convection etc.
!            --------------------------

! TRANSFER NOT ADVECTED VARIABLES INTO PGFLT1
IF (LCVPRO.OR.LGPCMT) THEN
   IF (.NOT.YUAL%LADV) THEN
      YDVARS%UAL%T1(KST:KEND,1:NFLEVG)=YDVARS%UAL%T0(KST:KEND,1:NFLEVG)
      YDVARS%UOM%T1(KST:KEND,1:NFLEVG)=YDVARS%UOM%T0(KST:KEND,1:NFLEVG)
   ENDIF
   IF (LCDDPRO) THEN
    IF (.NOT.YDAL%LADV) THEN
      YDVARS%DAL%T1(KST:KEND,1:NFLEVG)=YDVARS%DAL%T0(KST:KEND,1:NFLEVG)
      YDVARS%DOM%T1(KST:KEND,1:NFLEVG)=YDVARS%DOM%T0(KST:KEND,1:NFLEVG)
    ENDIF
   ENDIF
ENDIF
IF (YUNEBH%LACTIVE.AND..NOT.YUNEBH%LADV)THEN
   YDVARS%UNEBH%T1(KST:KEND,1:NFLEVG)=YDVARS%UNEBH%T0(KST:KEND,1:NFLEVG)
ENDIF
IF (YUEN%LACTIVE.AND..NOT.YUEN%LADV) THEN
   YDVARS%UEN%T1(KST:KEND,1:NFLEVG)=YDVARS%UEN%T0(KST:KEND,1:NFLEVG)
ENDIF
IF (YTTE%LACTIVE.AND..NOT.YTTE%LADV) THEN
   YDVARS%TTE%T1(KST:KEND,1:NFLEVG)=YDVARS%TTE%T0(KST:KEND,1:NFLEVG)
ENDIF
IF (YMXL%LACTIVE.AND..NOT.YMXL%LADV) THEN
   YDVARS%MXL%T1(KST:KEND,1:NFLEVG)=YDVARS%MXL%T0(KST:KEND,1:NFLEVG)
ENDIF
IF (YSHTUR%LACTIVE.AND..NOT.YSHTUR%LADV) THEN
   YDVARS%SHTUR%T1(KST:KEND,1:NFLEVG)=YDVARS%SHTUR%T0(KST:KEND,1:NFLEVG)
ENDIF
IF (YFQTUR%LACTIVE.AND..NOT.YFQTUR%LADV) THEN
   YDVARS%FQTUR%T1(KST:KEND,1:NFLEVG)=YDVARS%FQTUR%T0(KST:KEND,1:NFLEVG)
ENDIF
IF (YFSTUR%LACTIVE.AND..NOT.YFSTUR%LADV) THEN
   YDVARS%FSTUR%T1(KST:KEND,1:NFLEVG)=YDVARS%FSTUR%T0(KST:KEND,1:NFLEVG)
ENDIF
IF (YRKTH%LACTIVE)THEN
   YDVARS%RKTH%T1(KST:KEND,1:NFLEVG)=YDVARS%RKTH%T0(KST:KEND,1:NFLEVG)
ENDIF
IF (YRKTQV%LACTIVE)THEN
   YDVARS%RKTQV%T1(KST:KEND,1:NFLEVG)=YDVARS%RKTQV%T0(KST:KEND,1:NFLEVG)
ENDIF
IF (YRKTQC%LACTIVE)THEN
   YDVARS%RKTQC%T1(KST:KEND,1:NFLEVG)=YDVARS%RKTQC%T0(KST:KEND,1:NFLEVG)
ENDIF

!        2.10  Surface variables.
!              ------------------

IF (LLDIAB.AND.LMPHYS.AND.(.NOT.LMPA).AND.(.NOT.LSFORCS)) THEN
  
  IF (.NOT.LMSE) THEN
    DO JLEV=0,NFLEVG
      DO JROF=KST,KEND
        ZFPLSN(JROF,JLEV)=YDMF_PHYS%OUT%FPLSN(JROF,JLEV)+YDMF_PHYS%OUT%FPLSG(JROF,JLEV)
      ENDDO
    ENDDO
    CALL CPTENDS( YDMODEL%YRML_PHY_MF, NPROMA, KST, KEND, NFLEVG, YSP_SBD%NLEVS, PDTPHY,&
     & YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLCN, ZFPLSN,&
     & YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH,&
     & YDMF_PHYS_SURF%GSP_SG%PA_T1,YDMF_PHYS%OUT%CT, ZC1, ZC2,&
     & YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FCLL, YDMF_PHYS%OUT%FCLN, YDMF_PHYS%OUT%FCS,&
     & ZFEVI,YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN,&
     & YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FGEL, YDMF_PHYS%OUT%FGELS, YDMF_PHYS%OUT%FLWSP, YDMF_PHYS%OUT%FONTE, YDMF_PHYS%OUT%FTR,&
     & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSP_SG%PR_T1,&
     & YDMF_PHYS%OUT%RUISL, YDMF_PHYS%OUT%RUISP, YDMF_PHYS%OUT%RUISS, YDMF_PHYS_SURF%GSP_SG%PF_T1, ZVEG,&
     & ZTDTS, ZTDTP, ZTDWS, ZTDWSI, ZTDWP, ZTDWPI, ZTDWL,&
     & ZTDSNS, ZTDALBNS, ZTDRHONS)  

    CALL CPWTS(YDSURF, YDMODEL%YRML_AOC%YRMCC,YDPHY,YDMODEL%YRML_PHY_MF%YRPHY1, NPROMA, KST,KEND, YSP_SBD%NLEVS, PDTPHY,&
     & ZTDTS, ZTDTP, ZTDWS, ZTDWSI, ZTDWP, ZTDWPI, ZTDWL,&
     & ZTDSNS, ZTDALBNS, ZTDRHONS,&
     & YDMF_PHYS_SURF%GSD_VP%PTPC,YDMF_PHYS_SURF%GSD_VP%PWPC,YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS_SURF%GSD_VV%PIVEG,&
     & YDMF_PHYS_SURF%GSP_RR%PT_T1,YDMF_PHYS_SURF%GSP_SB%PT_T1,YDMF_PHYS_SURF%GSP_RR%PW_T1,&
     & YDMF_PHYS_SURF%GSP_RR%PIC_T1,&
     & YDMF_PHYS_SURF%GSP_SB%PQ_T1,YDMF_PHYS_SURF%GSP_SB%PTL_T1,YDMF_PHYS_SURF%GSP_RR%PFC_T1,&
     & YDMF_PHYS_SURF%GSP_SG%PF_T1,YDMF_PHYS_SURF%GSP_SG%PA_T1,YDMF_PHYS_SURF%GSP_SG%PR_T1)  
  ELSE
    IF (LLXFUMSE) THEN
      DO JROF=KST,KEND
        YDMF_PHYS_SURF%GSP_RR%PT_T0(JROF)=PGPAR(JROF,MVTS)
      ENDDO
    ELSE
      DO JROF=KST,KEND
        YDMF_PHYS_SURF%GSP_RR%PT_T1(JROF)=PGPAR(JROF,MVTS)
      ENDDO
    ENDIF
  ENDIF
  IF(LNUDG)THEN

    ! * Calculation of IPQ since the old pointers
    !   MSLB1[X]9 (=MSLB1GFL9+IP[X]) do not exist any longer in PTRSLB1.
    IPQ=(YQ%MP_SL1-1)*(NFLEVG+2*NFLSUL)

    CALL CPNUDG ( NPROMA, KST,KEND, NFNUDG, NFLEVG, IBLK,&
     & XPNUDG,&
     & YDMF_PHYS_SURF%GSD_VF%PNUDM,&
     & YDMF_PHYS_SURF%GSP_RR%PT_T1,YDMF_PHYS_SURF%GSP_RR%PW_T1,&
     & YDMF_PHYS_SURF%GSP_SB%PQ_T1,YDMF_PHYS_SURF%GSP_SG%PF_T1,&
     & PB1(1,ISLB1T9+1-NFLSA),PB1(1,ISLB1GFL9+IPQ+1-NFLSA),&
     & PB1(1,ISLB1U9+1-NFLSA),PB1(1,ISLB1V9+1-NFLSA),&
     & PB1(1,MSLB1SP9),&
     & YDVARS%T%T0,YDVARS%Q%T0,YDVARS%U%T0,&
     & YDVARS%V%T0,YDCPG_PHY0%PREHYD(:,NFLEVG),YDGSGEOM%GM,YDMF_PHYS_SURF%GSD_VF%PLSM)
  ENDIF
ENDIF

!        2.11 Evolution of CVV (GY)
!             ---------------------

IF(LCVPGY) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      YDVARS%CVV%T1(JROF,JLEV)=YDVARS%CVV%T0(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

!    -------------------------------------------------------------------

!*       3.    Simplified physics.
!              -------------------

!        3.1  Preliminary calculations necessary for simplified physics.
!             ----------------------------------------------------------

IF (LSIMPH) THEN

  ! * read grid-point transmission coefficients for simplified physics.
  IF (LRAYSP.AND.LRCOEF.AND.(NSTEP > 1.OR.LTLADDIA)) THEN
    IFIELDSS=NG3SR*NFLEVG
    CALL RDRADCOEF(YDGEOMETRY,YDRCOEF,KST,KEND,KSTGLO,IFIELDSS,ZRADTC,ZAC_HC)
  ENDIF

ENDIF

!        3.2  Simplified physics.
!             -------------------

IF (LSIMPH) THEN

  TSPHY = MAX(PDTPHY,1.0_JPRB)

  IF (LTWOTL) THEN

    CALL APLPASSH (YDPHY,YDMODEL%YRML_PHY_MF%YRPHY1,KST,KEND,NPROMA,NFLEVG,&
     & YDCPG_PHY0%PREHYD,YDVARS%Q%T0,&
     & YDMF_PHYS_SURF%GSP_SG%PF_T1,YDMF_PHYS_SURF%GSP_RR%PT_T1,YDMF_PHYS_SURF%GSP_RR%PW_T1,&
     & YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS_SURF%GSD_VF%PVEG,&
     & ZQS1)  

    IF (.NOT.LLDIAB) THEN
      CALL APLPASSH (YDPHY,YDMODEL%YRML_PHY_MF%YRPHY1,KST,KEND,NPROMA,NFLEVG,&
       & YDCPG_PHY0%PREHYD,YDVARS%Q%T0,&
       & YDMF_PHYS_SURF%GSP_SG%PF_T0,YDMF_PHYS_SURF%GSP_RR%PT_T0,YDMF_PHYS_SURF%GSP_RR%PW_T0,&
       & YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS_SURF%GSD_VF%PVEG,&
       & YDCPG_MISC%QS)  
    ENDIF

    CALL APLPARS(YDGEOMETRY,YDRCOEF,YDMODEL%YRML_PHY_MF,KST,KEND,NPROMA,1,NFLEVG,NTSSG,NSTEP,&
     & YDCPG_DYN0%PHI,YDCPG_PHY0%PREHYD,YDCPG_DYN0%PHIF,YDCPG_PHY0%PREHYDF,YDCPG_PHY0%XYB%DELP,YDCPG_PHY0%XYB%RDELP,&
     & YDVARS%U%T0,YDVARS%V%T0,YDVARS%T%T0,YDVARS%Q%T0,YDCPG_DYN0%RCP%CP,ZCVGQ,&
     & YDMF_PHYS_SURF%GSP_SG%PF_T0,YDMF_PHYS_SURF%GSP_RR%PT_T0,YDCPG_MISC%QS,&
     & YDMF_PHYS_SURF%GSP_RR%PT_T1,ZQS1,&
     & YDMF_PHYS_SURF%GSD_VF%PGETRL,YDMF_PHYS_SURF%GSD_VF%PLSM,&
     & YDMF_PHYS_SURF%GSD_VF%PZ0F,YDMF_PHYS_SURF%GSD_VF%PVRLAN,YDMF_PHYS_SURF%GSD_VF%PVRLDI,&
     & YDMF_PHYS_SURF%GSD_VF%PALBF,ZMU0,YDGSGEOM%GM,&
     & ZAC_HC,ZMCOR,&
     & ZMRAB3C,ZMRAB3N,&
     & ZMRAB4C,ZMRAB4N,&
     & ZMRAB6C,ZMRAB6N,&
     & ZMRAT1C,ZMRAT1N,&
     & ZMRAT2C,ZMRAT2N,&
     & ZMRAT3C,ZMRAT3N,&
     & ZMRAT4C,ZMRAT4N,&
     & ZMRAT5C,ZMRAT5N,&
     & YDMF_PHYS%OUT%DIFCQ,YDMF_PHYS%OUT%DIFCS,YDMF_PHYS%OUT%DIFTQ,YDMF_PHYS%OUT%DIFTS,&
     & YDMF_PHYS%OUT%FCCQL,YDMF_PHYS%OUT%FCCQN,YDMF_PHYS%OUT%FCSQL,YDMF_PHYS%OUT%FCSQN,&
     & YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FPLSN,YDMF_PHYS%OUT%FRSO,YDMF_PHYS%OUT%FRTH,&
     & YDMF_PHYS%OUT%STRCU,YDMF_PHYS%OUT%STRCV,YDMF_PHYS%OUT%STRDU,YDMF_PHYS%OUT%STRDV,YDMF_PHYS%OUT%STRTU,YDMF_PHYS%OUT%STRTV,&
     & YDMF_PHYS%OUT%STRMU,YDMF_PHYS%OUT%STRMV,YDMF_PHYS%OUT%FRMH)

  ELSE

    CALL APLPASSH (YDPHY,YDMODEL%YRML_PHY_MF%YRPHY1,KST,KEND,NPROMA,NFLEVG,&
     & YDCPG_PHY0%PREHYD,YDVARS%Q%T0,&
     & YDMF_PHYS_SURF%GSP_SG%PF_T1,YDMF_PHYS_SURF%GSP_RR%PT_T1,YDMF_PHYS_SURF%GSP_RR%PW_T1,&
     & YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS_SURF%GSD_VF%PVEG,&
     & ZQS1)  

    IF (.NOT.LLDIAB) THEN
      CALL APLPASSH (YDPHY,YDMODEL%YRML_PHY_MF%YRPHY1,KST,KEND,NPROMA,NFLEVG,&
       & YDCPG_PHY9%PREHYD,YDVARS%Q%T9,&
       & YDMF_PHYS_SURF%GSP_SG%PF_T0,YDMF_PHYS_SURF%GSP_RR%PT_T9,YDMF_PHYS_SURF%GSP_RR%PW_T9,&
       & YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS_SURF%GSD_VF%PVEG,&
       & YDCPG_MISC%QS)  
    ENDIF

    CALL APLPARS(YDGEOMETRY,YDRCOEF,YDMODEL%YRML_PHY_MF,KST,KEND,NPROMA,1,NFLEVG,NTSSG,NSTEP,&
     & YDCPG_DYN9%PHI,YDCPG_PHY9%PREHYD,YDCPG_DYN9%PHIF,YDCPG_PHY9%PREHYDF,YDCPG_PHY9%XYB%DELP,YDCPG_PHY9%XYB%RDELP,&
     & YDVARS%U%T9,YDVARS%V%T9,YDVARS%T%T9,YDVARS%Q%T9,YDCPG_DYN9%RCP%CP,ZCVGQ,&
     & YDMF_PHYS_SURF%GSP_SG%PF_T9,YDMF_PHYS_SURF%GSP_RR%PT_T9,YDCPG_MISC%QS,&
     & YDMF_PHYS_SURF%GSP_RR%PT_T1,ZQS1,&
     & YDMF_PHYS_SURF%GSD_VF%PGETRL,YDMF_PHYS_SURF%GSD_VF%PLSM,&
     & YDMF_PHYS_SURF%GSD_VF%PZ0F,YDMF_PHYS_SURF%GSD_VF%PVRLAN,YDMF_PHYS_SURF%GSD_VF%PVRLDI,&
     & YDMF_PHYS_SURF%GSD_VF%PALBF,ZMU0,YDGSGEOM%GM,&
     & ZAC_HC,ZMCOR,&
     & ZMRAB3C,ZMRAB3N,&
     & ZMRAB4C,ZMRAB4N,&
     & ZMRAB6C,ZMRAB6N,&
     & ZMRAT1C,ZMRAT1N,&
     & ZMRAT2C,ZMRAT2N,&
     & ZMRAT3C,ZMRAT3N,&
     & ZMRAT4C,ZMRAT4N,&
     & ZMRAT5C,ZMRAT5N,&
     & YDMF_PHYS%OUT%DIFCQ,YDMF_PHYS%OUT%DIFCS,YDMF_PHYS%OUT%DIFTQ,YDMF_PHYS%OUT%DIFTS,&
     & YDMF_PHYS%OUT%FCCQL,YDMF_PHYS%OUT%FCCQN,YDMF_PHYS%OUT%FCSQL,YDMF_PHYS%OUT%FCSQN,&
     & YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FPLSN,YDMF_PHYS%OUT%FRSO,YDMF_PHYS%OUT%FRTH,&
     & YDMF_PHYS%OUT%STRCU,YDMF_PHYS%OUT%STRCV,YDMF_PHYS%OUT%STRDU,YDMF_PHYS%OUT%STRDV,YDMF_PHYS%OUT%STRTU,YDMF_PHYS%OUT%STRTV,&
     & YDMF_PHYS%OUT%STRMU,YDMF_PHYS%OUT%STRMV,YDMF_PHYS%OUT%FRMH)

  ENDIF

ENDIF

!        3.3  Store the model trajectory at t-dt (leap-frog) or t (sl2tl).
!             ------------------------------------------------------------

IF (LTRAJPS) THEN
  IF (LTWOTL) THEN
    PTRAJ_PHYS%PQSSMF(KST:KEND)=YDCPG_MISC%QS(KST:KEND)
    PTRAJ_PHYS%PTSMF(KST:KEND) =YDMF_PHYS_SURF%GSP_RR%PT_T0(KST:KEND)
    PTRAJ_PHYS%PSNSMF(KST:KEND)=YDMF_PHYS_SURF%GSP_SG%PF_T0(KST:KEND,1)
  ELSE
    CALL WRPHTRAJM(YDGEOMETRY,YDSIMPHL,KST,KEND,PTRAJ_PHYS,&
     & YDVARS%U%T9,YDVARS%V%T9,YDVARS%T%T9,&
     & YDVARS%Q%T9,YDVARS%L%T9,YDVARS%I%T9,YDVARS%SP%T9)  

    PTRAJ_PHYS%PQSSMF(KST:KEND)=YDCPG_MISC%QS(KST:KEND)
    PTRAJ_PHYS%PTSMF(KST:KEND) =YDMF_PHYS_SURF%GSP_RR%PT_T9(KST:KEND)
    PTRAJ_PHYS%PSNSMF(KST:KEND)=YDMF_PHYS_SURF%GSP_SG%PF_T9(KST:KEND,1)
  ENDIF
  IF (LPRTTRAJ.AND.PTRAJ_PHYS%LASTCHUNK) WRITE(NULOUT,*)'GREPTRAJ STORE TRAJ_PHYS in MF_PHYS'
ENDIF

!        3.4  Computation of tendencies T,u,v and Q.
!             --------------------------------------

IF (LSIMPH) THEN

  TSPHY = MAX(PDTPHY,1.0_JPRB)

  IF (LTWOTL) THEN
    CALL CPTENDSM (YDPHY,NPROMA,KST,KEND,NFLEVG,&
     & YDMF_PHYS%OUT%DIFCQ,YDMF_PHYS%OUT%DIFCS,YDMF_PHYS%OUT%DIFTQ,YDMF_PHYS%OUT%DIFTS,&
     & YDMF_PHYS%OUT%FCCQL,YDMF_PHYS%OUT%FCCQN,YDMF_PHYS%OUT%FCSQL,YDMF_PHYS%OUT%FCSQN,&
     & YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FPLSN,&
     & YDMF_PHYS%OUT%FRSO,YDMF_PHYS%OUT%FRTH,&
     & YDMF_PHYS%OUT%STRCU,YDMF_PHYS%OUT%STRCV,YDMF_PHYS%OUT%STRDU,YDMF_PHYS%OUT%STRDV,&
     & YDMF_PHYS%OUT%STRTU,YDMF_PHYS%OUT%STRTV,&
     & YDMF_PHYS%OUT%STRMU,YDMF_PHYS%OUT%STRMV,YDMF_PHYS%OUT%FRMH,&
     & YDCPG_PHY0%XYB%RDELP,YDCPG_DYN0%PHIF,&
     & YDVARS%U%T0,YDVARS%V%T0,YDVARS%T%T0,YDVARS%Q%T0,&
     & YDCPG_MISC%QS,YDMF_PHYS_SURF%GSP_RR%PT_T0,YDOROG%OROG,&
     & YDMF_PHYS%OUT%TENDU,YDMF_PHYS%OUT%TENDV,ZTENDH,ZTENDQ,&
     & YDMF_PHYS%OUT%FHPCL,YDMF_PHYS%OUT%FHPCN,YDMF_PHYS%OUT%FHPSL,YDMF_PHYS%OUT%FHPSN,&
     & YDMF_PHYS%OUT%FHSCL,YDMF_PHYS%OUT%FHSCN,YDMF_PHYS%OUT%FHSSL,YDMF_PHYS%OUT%FHSSN)  
  ELSE
    CALL CPTENDSM (YDPHY,NPROMA,KST,KEND,NFLEVG,&
     & YDMF_PHYS%OUT%DIFCQ,YDMF_PHYS%OUT%DIFCS,YDMF_PHYS%OUT%DIFTQ,YDMF_PHYS%OUT%DIFTS,&
     & YDMF_PHYS%OUT%FCCQL,YDMF_PHYS%OUT%FCCQN,YDMF_PHYS%OUT%FCSQL,YDMF_PHYS%OUT%FCSQN,&
     & YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FPLSN,&
     & YDMF_PHYS%OUT%FRSO,YDMF_PHYS%OUT%FRTH,&
     & YDMF_PHYS%OUT%STRCU,YDMF_PHYS%OUT%STRCV,YDMF_PHYS%OUT%STRDU,YDMF_PHYS%OUT%STRDV,&
     & YDMF_PHYS%OUT%STRTU,YDMF_PHYS%OUT%STRTV,&
     & YDMF_PHYS%OUT%STRMU,YDMF_PHYS%OUT%STRMV,YDMF_PHYS%OUT%FRMH,&
     & YDCPG_PHY9%XYB%RDELP,YDCPG_DYN9%PHIF,&
     & YDVARS%U%T9,YDVARS%V%T9,YDVARS%T%T9,YDVARS%Q%T9,&
     & YDCPG_MISC%QS,YDMF_PHYS_SURF%GSP_RR%PT_T9,YDOROG%OROG,&
     & YDMF_PHYS%OUT%TENDU,YDMF_PHYS%OUT%TENDV,ZTENDH,ZTENDQ,&
     & YDMF_PHYS%OUT%FHPCL,YDMF_PHYS%OUT%FHPCN,YDMF_PHYS%OUT%FHPSL,YDMF_PHYS%OUT%FHPSN,&
     & YDMF_PHYS%OUT%FHSCL,YDMF_PHYS%OUT%FHSCN,YDMF_PHYS%OUT%FHSSL,YDMF_PHYS%OUT%FHSSN)  
  ENDIF

ENDIF

!        3.5  Computation of evolution of T, u, v and Q.
!             ------------------------------------------

IF (LSIMPH) THEN

  CALL CPUTQYS(YDGEOMETRY%YRDIMV,YDGMV,YGFL,YDPTRSLB1,&
   & YDSTOPH, YDPHY2, &
   & NPROMA,KST,KEND,NFLEVG,PDTPHY,IPGFL,&
   & ISLB1U9,ISLB1V9,ISLB1T9,ISLB1GFL9,&
   & ZTENDH, ZTENDQ, YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV,&
   & YDMF_PHYS%FOR%U, YDMF_PHYS%FOR%V, YDMF_PHYS%FOR%T, YDMF_PHYS%FOR%Q, &
   & YDCPG_DYN0%RCP%CP,YDCPG_PHY0%XYB%DELP,YDVARS%T%T0,YDVARS%U%T0,YDVARS%V%T0,&
   & YDCPG_DYN9%RCP%CP,YDCPG_PHY9%XYB%DELP,YDVARS%T%T9,YDVARS%U%T9,YDVARS%V%T9,&
   & PB1, PGMVT1, PGFLT1,&
   & YDMF_PHYS%OUT%FDIS)

ENDIF

!     ------------------------------------------------------------------

!*       4.    AROME  physics.
!              ---------------

IF (LMPA) THEN

  !      4.1  CALL APL_AROME
  !           --------------

  IPTR(:) = 0 ! means no fields defined in ZGFLTENDR at start ; > 0 means defined.

  ! * ZTENDR     
  IRR=1
  IPTR_CONT = IRR
  IF (YQ%LACTIVE) THEN
    IPTR(YQ%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
  ENDIF   
  IF (YL%LACTIVE) THEN
    IPTR(YL%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
  ENDIF   
  IF (YR%LACTIVE) THEN
    IPTR(YR%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
  ENDIF   
  IF (YI%LACTIVE) THEN
    IPTR(YI%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
  ENDIF   
  IF (YS%LACTIVE) THEN
    IPTR(YS%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
  ENDIF   
  IF (YG%LACTIVE) THEN
    IPTR(YG%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
  ENDIF   
  IF (YH%LACTIVE) THEN
    IPTR(YH%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
  ENDIF   

  ! * ZTENDTKE   
  IF (YTKE%LACTIVE) THEN
    IPTR(YTKE%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
    IPTRTKE=IPTR(YTKE%MP1)
  ELSE
    IPTRTKE=0
  ENDIF
  ! * ZTENDEFB1 ZTENDEFB2 ZTENDEFB3
  IF (YEFB1%LACTIVE) THEN
    IPTR(YEFB1%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
    IEFB1 = IPTR(YEFB1%MP1)
  ELSE
    IEFB1 = 0
  ENDIF   
  IF (YEFB2%LACTIVE) THEN
    IPTR(YEFB2%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
    IEFB2 = IPTR(YEFB2%MP1)
  ELSE
    IEFB2 = 0
  ENDIF
  IF (YEFB3%LACTIVE) THEN   
    IPTR(YEFB3%MP1) = IPTR_CONT ; IPTR_CONT = IPTR_CONT+1
    IEFB3 = IPTR(YEFB3%MP1)
  ELSE
    IEFB3 = 0
  ENDIF

  ! * ZTENDEXT
  IF (NGFL_EXT > 0) THEN
    DO JGFL=1,NGFL_EXT
      IF (YEXT(JGFL)%LACTIVE) THEN   
        IPTR(YEXT(JGFL)%MP1)=IPTR_CONT
        IPTR_CONT = IPTR_CONT+1   
      ENDIF
    ENDDO  
    IPTREXT=IPTR(YEXT(1)%MP1)  
  ELSE
    IPTREXT=0
  ENDIF

  ! * LIMA
  IF (NLIMA > 0) THEN
    DO JGFL=1,NLIMA
      IF (YLIMA(JGFL)%LACTIVE) THEN   
        IPTR(YLIMA(JGFL)%MP1)=IPTR_CONT
        IPTR_CONT = IPTR_CONT+1   
      ENDIF
    ENDDO  
    IPTRLIMA=IPTR(YLIMA(1)%MP1)  
  ELSE
    IPTRLIMA=0
  ENDIF
 

  ! If an incorrect address is used, then the initialization below will detect it :
  ZTENDGFLR(:,:,0)=HUGE(1._JPRB)
  
  IF (LMFSHAL .AND. CMF_UPDRAFT=='DUAL') THEN
    IMAXDRAFT=3
  ELSE
    IMAXDRAFT=0
  ENDIF

  IF(LTWOTL) THEN
    CALL APL_AROME(YDGEOMETRY,YDCPG_MISC,YDMF_PHYS,YDMF_PHYS_SURF,YDVARS,YDSURF,YDCFU,YDXFU,YDMODEL, KBL, KGPCOMP, KST  , KEND   , NPROMA ,&
     & 1      , NFLEVG , NSTEP  ,&
     & IMAXDRAFT, NTSSG, YDCFU%NFRRC, PDTPHY, LLXFUMSE,YDCSGEOM%RINDX,YDCSGEOM%RINDY,&
     & YDGSGEOM%GEMU,YDGSGEOM%GELAM,YDOROG%OROG,YDGSGEOM%GM,&
     & ZMU0,ZMU0LU,ZMU0M,ZMU0N,&
     & YDGSGEOM%GECLO,YDGSGEOM%GESLO,&
     & PGP2DSDT, ZGP2DSPP, &
     & YDCPG_DYN0%PHI,YDCPG_DYN0%PHIF,YDCPG_PHY0%PRE,YDCPG_PHY0%PREF,YDCPG_PHY0%XYB%RDELP,YDCPG_PHY0%XYB%DELP,&
     & YDVARS%T%T0,YDVARS%Q%T0,&
     & YDCPG_DYN0%RCP%CP,YDCPG_DYN0%RCP%R,YDCPG_PHY0%XYB%ALPH,YDCPG_PHY0%XYB%LNPR,&
     & YDVARS%L%T0,YDVARS%I%T0,YDVARS%R%T0,YDVARS%S%T0,&
     & YDVARS%G%T0,YDVARS%H%T0,ZP1LIMA0,YDVARS%TKE%T0,&
     & YDVARS%EFB1%T0,YDVARS%EFB2%T0,YDVARS%EFB3%T0,&
     & YDVARS%SRC%T0,ZP1EXT0,&
     & YDVARS%U%T0,YDVARS%V%T0, YDCPG_PHY0%W,ZEDR,&
     & PGPAR,&
     ! outputs
     & ZTENDT, ZTENDGFLR(1,1,IRR), ZTENDW, ZTENDGFLR(1,1,IPTRLIMA), ZTENDGFLR(1,1,IPTRTKE),&
     & ZTENDGFLR(1,1,IEFB1),ZTENDGFLR(1,1,IEFB2),ZTENDGFLR(1,1,IEFB3),&
     & ZTENDGFLR(1,1,IPTREXT),&
     & ZDPRECIPS, ZDPRECIPS2,ZP1EZDIAG,&
     & YLPROCSET,YDDDH )
  ELSE    
    CALL APL_AROME(YDGEOMETRY,YDCPG_MISC,YDMF_PHYS,YDMF_PHYS_SURF,YDVARS,YDSURF,YDCFU,YDXFU,YDMODEL, KBL, KGPCOMP, KST  , KEND   , NPROMA ,&
     & 1      , NFLEVG , NSTEP  ,&
     & IMAXDRAFT, NTSSG, YDCFU%NFRRC, PDTPHY, LLXFUMSE,YDCSGEOM%RINDX,YDCSGEOM%RINDY,&
     & YDGSGEOM%GEMU,YDGSGEOM%GELAM,YDOROG%OROG,YDGSGEOM%GM,&
     & ZMU0,ZMU0LU,ZMU0M,ZMU0N,&
     & YDGSGEOM%GECLO,YDGSGEOM%GESLO,&
     & PGP2DSDT, ZGP2DSPP, &
     & YDCPG_DYN9%PHI,YDCPG_DYN9%PHIF,YDCPG_PHY9%PRE,YDCPG_PHY9%PREF,YDCPG_PHY9%XYB%RDELP,YDCPG_PHY9%XYB%DELP,&
     & YDVARS%T%T9,YDVARS%Q%T9,&
     & YDCPG_DYN9%RCP%CP,YDCPG_DYN9%RCP%R,YDCPG_PHY9%XYB%ALPH,YDCPG_PHY9%XYB%LNPR,&
     & YDVARS%L%T9,YDVARS%I%T9,YDVARS%R%T9,YDVARS%S%T9,&
     & YDVARS%G%T9,YDVARS%H%T9,ZP1LIMA9,YDVARS%TKE%T9,&
     & YDVARS%EFB1%T9,YDVARS%EFB2%T9,YDVARS%EFB3%T9,&
     & YDVARS%SRC%T9,ZP1EXT9,&
     & YDVARS%U%T9,YDVARS%V%T9, YDCPG_PHY9%W,ZEDR ,&
     & PGPAR,&
     ! outputs
     & ZTENDT, ZTENDGFLR(1,1,IRR), ZTENDW,  ZTENDGFLR(1,1,IPTRLIMA), ZTENDGFLR(1,1,IPTRTKE),&
     & ZTENDGFLR(1,1,IEFB1),ZTENDGFLR(1,1,IEFB2),ZTENDGFLR(1,1,IEFB3),&
     & ZTENDGFLR(1,1,IPTREXT),&
     & ZDPRECIPS, ZDPRECIPS2,ZP1EZDIAG,&
     & YLPROCSET,YDDDH )
  ENDIF 
  !Save surface temperature
  IF (LMSE.OR.LSFORCS) THEN
    IF (LLXFUMSE) THEN
      DO JROF=KST,KEND
        YDMF_PHYS_SURF%GSP_RR%PT_T0(JROF)=PGPAR(JROF,MVTS)
      ENDDO
    ELSE
      DO JROF=KST,KEND
        YDMF_PHYS_SURF%GSP_RR%PT_T1(JROF)=PGPAR(JROF,MVTS)
      ENDDO
    ENDIF 
  ENDIF  
  !      4.2  COMPUTE THE PHYS. TENDENCY FOR "T" AND "w"
  !           ------------------------------------------


  IF (LVERTFE.AND.LVFE_GWMPA) THEN
    ! * case LVFE_GWMPA not yet coded.
    !   (in this case ZGWT1 must be computed at full levels and
    !   not at half levels)
    CALL ABOR1(' MF_PHYS: case LVFE_GWMPA not yet coded if LMPA=T!')
  ENDIF

  ZDT = PDTPHY

  ZTENDD=0.0_JPRB

  ! * compute ZTT1:
  IF (LSLAG.AND.LTWOTL) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        ZTT1(JROF,JLEV)=YDVARS%T%T0(JROF,JLEV)+ZDT*ZTENDT(JROF,JLEV)
      ENDDO
    ENDDO
  ELSE
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        ZTT1(JROF,JLEV)=YDVARS%T%T9(JROF,JLEV)+ZDT*ZTENDT(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF

  ! * compute ZGWT1 = tendency of gw:
  IF (LNHDYN) THEN
    ! Valid for LVFE_GWMPA=F only; ZGWT1 assumed to be half level values.
    DO JLEV=1,NFLEVG-1
      DO JROF=KST,KEND
        ZGWT1(JROF,JLEV)=0.5_JPRB*RG*(ZTENDW(JROF,JLEV)+ZTENDW(JROF,JLEV+1))
      ENDDO
    ENDDO
    DO JROF=KST,KEND
      ZGWT1(JROF,NFLEVG)=0.0_JPRB
      ZGWT1(JROF,0)=0.0_JPRB
    ENDDO
  ENDIF

  ! * convert gw tendency in d tendency:
  IF(LNHDYN) THEN

    IF (LGWADV) THEN
      ZTENDD(KST:KEND,1:NFLEVG)=ZGWT1(KST:KEND,1:NFLEVG)
    ELSE

      ! * Provide the appropriate version of (RT) at t+dt for GNHGW2SVDAROME:
      IF (L_RDRY_VD) THEN
        ! Use Rd because "dver" is currently defined with Rd.
        ZRTT1(KST:KEND,1:NFLEVG)=RD*ZTT1(KST:KEND,1:NFLEVG)
      ELSE
        ! Use "moist R" because "dver" is defined with "moist R".
        ! Unfortunately, R(t+dt) is not yet available there, use R(t) instead.
        ! "Moist R" tendency is neglected in the below call to GNHGW2SVDAROME.
        DO JLEV=1,NFLEVG
          DO JROF=KST,KEND
            ZRTT1(JROF,JLEV)=YDCPG_DYN0%RCP%R(JROF,JLEV)*ZTT1(JROF,JLEV)
          ENDDO
        ENDDO
      ENDIF

      ! * Do conversion:
      IF (LSLAG.AND.LTWOTL) THEN
        CALL GNHGW2SVDAROME(YDGEOMETRY,KST,KEND,YDCPG_PHY0%PREHYDF,YDCPG_PHY0%XYB%LNPR,ZRTT1,YDCPG_PHY0%PREF,ZGWT1,&
         & ZTENDD)  
      ELSE
        CALL GNHGW2SVDAROME(YDGEOMETRY,KST,KEND,YDCPG_PHY9%PREHYDF,YDCPG_PHY9%XYB%LNPR,ZRTT1,YDCPG_PHY9%PREF,ZGWT1,&
         & ZTENDD)  
      ENDIF

    ENDIF
  ELSE
    ZTENDD=0.0_JPRB
  ENDIF

  !      4.3  PUT THE TENDENCIES IN PB1/GFLT1/GMVT1.
  !           --------------------------------------


  IF ( LINTFLEX ) THEN
    IF (LTWOTL) THEN
      CALL CPTEND_FLEX( YDLDDH,YDMDDH,YGFL,YDPHY,NPROMA, KST, KEND, NFLEVG,YDGSGEOM%GNORDL,YDGSGEOM%GNORDM,&
       & YDCPG_PHY0%XYB%DELP ,&
       & YDCPG_PHY0%XYB%RDELP, YDCPG_DYN0%RCP%CP,&
       & YDVARS%U%T0,YDVARS%V%T0,YDVARS%T%T0,YDMF_PHYS_SURF%GSP_RR%PT_T0,&
       & PGFL,&
       & YLPROCSET,&
       & YDMF_PHYS%OUT%TENDU , YDMF_PHYS%OUT%TENDV , ZTENDH , ZTENDGFL,&
       & YDMF_PHYS%OUT%FHSCL ,YDMF_PHYS%OUT%FHSCN,YDMF_PHYS%OUT%FHSSL,YDMF_PHYS%OUT%FHSSN,&
       & YDMF_PHYS%OUT%FHPCL ,YDMF_PHYS%OUT%FHPCN,YDMF_PHYS%OUT%FHPSL,YDMF_PHYS%OUT%FHPSN,&
       & ZFHP   ,ZFP   ,  YDMF_PHYS%OUT%FEPFP, YDMF_PHYS%OUT%FCMPCQ, YDMF_PHYS%OUT%FCMPSN, YDMF_PHYS%OUT%FCMPSL,YDDDH )
    ELSE
      CALL CPTEND_FLEX( YDLDDH,YDMDDH,YGFL,YDPHY,NPROMA, KST, KEND, NFLEVG,YDGSGEOM%GNORDL,YDGSGEOM%GNORDM,&
       & YDCPG_PHY9%XYB%DELP ,&
       & YDCPG_PHY9%XYB%RDELP, YDCPG_DYN9%RCP%CP,&
       & YDVARS%U%T9,YDVARS%V%T9,YDVARS%T%T9,YDMF_PHYS_SURF%GSP_RR%PT_T9,&
       & PGFL,&
       & YLPROCSET,&
       & YDMF_PHYS%OUT%TENDU , YDMF_PHYS%OUT%TENDV , ZTENDH , ZTENDGFL,&
       & YDMF_PHYS%OUT%FHSCL ,YDMF_PHYS%OUT%FHSCN,YDMF_PHYS%OUT%FHSSL,YDMF_PHYS%OUT%FHSSN,&
       & YDMF_PHYS%OUT%FHPCL ,YDMF_PHYS%OUT%FHPCN,YDMF_PHYS%OUT%FHPSL,YDMF_PHYS%OUT%FHPSN,&
       & ZFHP   ,ZFP   ,  YDMF_PHYS%OUT%FEPFP, YDMF_PHYS%OUT%FCMPCQ, YDMF_PHYS%OUT%FCMPSN, YDMF_PHYS%OUT%FCMPSL,YDDDH )      
    ENDIF
    
    CALL CPUTQY(YDGEOMETRY%YRDIMV,YDGMV,YGFL,YDPTRSLB1,YDPHY,NPROMA,KST,KEND,NFLEVG,PDTPHY,IPGFL,&
     & ISLB1T9,ISLB1U9,ISLB1V9,ISLB1VD9,ISLB1GFL9,&
     & ZTENDH, ZTENDT, YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDU, ZTENDV, ZTENDD, ZTENDGFL,&
     & YDCPG_DYN0%RCP%CP,YDCPG_PHY0%XYB%DELP,YDVARS%T%T0,YDVARS%U%T0,YDVARS%V%T0,&
     & YDCPG_DYN9%RCP%CP,YDCPG_PHY9%XYB%DELP,YDVARS%T%T9,YDVARS%U%T9,YDVARS%V%T9,&
     & PB1, PGMVT1, PGFLT1,&
     & YDMF_PHYS%OUT%FDIS)    
     
  ELSE
  
    ! start ZTENDGFLR at 1 because it is dimensionned (:,:,0:n)
    CALL CPUTQY_AROME(YDGEOMETRY%YRDIMV,YDGMV,YGFL,YDPTRSLB1,NPROMA,KST,KEND,NFLEVG,ZDT,IPGFL,IPTR,&
       & ISLB1U9,ISLB1V9,ISLB1T9,ISLB1GFL9,ISLB1VD9 ,&
       & ZTENDT, ZTENDGFLR(:,:,1), YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDD ,&
       & PB1, PGMVT1, PGFLT1)
  ENDIF
  
ENDIF

!     ------------------------------------------------------------------

!*       5.    Final calculations.
!              -------------------

! * Restore the initial value of some pseudo-historical surface buffers
!   if relevant.
IF (LL_SAVE_PHSURF) THEN
  IF(YSD_VV%YHV%LSET) YDMF_PHYS_SURF%GSD_VV%PHV(1:NPROMA)=ZHV(1:NPROMA)
  IF(YSD_VF%YZ0F%LSET) YDMF_PHYS_SURF%GSD_VF%PZ0F(1:NPROMA)=ZGZ0F(1:NPROMA)
  IF(YSD_VV%YZ0H%LSET) YDMF_PHYS_SURF%GSD_VV%PZ0H(1:NPROMA)=ZGZ0HF(1:NPROMA)
  IF(YSD_VH%YPBLH%LSET) YDMF_PHYS_SURF%GSD_VH%PPBLH(1:NPROMA)=ZPBLH(1:NPROMA)
  IF(YSD_VH%YSPSH%LSET) YDMF_PHYS_SURF%GSD_VH%PSPSH(1:NPROMA)=ZFHPS(1:NPROMA)
  IF(YSD_VH%YQSH%LSET)  YDMF_PHYS_SURF%GSD_VH%PQSH(1:NPROMA)=ZQSH(1:NPROMA)
  IF(YSD_VK%YUDGRO%LSET) YDMF_PHYS_SURF%GSD_VK%PUDGRO(1:NPROMA)=ZUDGRO(1:NPROMA)
  IF(LCVPRO.OR.LGPCMT) THEN
    YDVARS%UAL%T0(:,:)=ZUDAL(:,:)
    YDVARS%UOM%T0(:,:)=ZUDOM(:,:)
    IF(LCDDPRO) THEN
      YDVARS%DAL%T0(:,:)=ZDDAL(:,:)
      YDVARS%DOM%T0(:,:)=ZDDOM(:,:)
    ENDIF
  ENDIF
  IF(YUNEBH%LACTIVE) YDVARS%UNEBH%T0(:,:)=ZUNEBH(:,:)
  IF(YUEN%LACTIVE)   YDVARS%UEN%T0(:,:)=ZENTCH(:,:)
ENDIF

! Store horizontal exchange coefficients (3D turbulence) to SL2 buffers
IF (L3DTURB) THEN
  DO JLEV=1,NFLEVG
    PB2(KST:KEND,MSLB2KAPPAM+JLEV-1)=ZKUROV_H(KST:KEND,JLEV)
    PB2(KST:KEND,MSLB2KAPPAH+JLEV-1)=ZKTROV_H(KST:KEND,JLEV)
  ENDDO
ENDIF

!--------------------------------------------------------------------
! BAYRAD
! Fill convective hydrometeors mixing ratio in GFL
!--------------------------------------------------------------------
IF((.NOT.LGPCMT).AND.(.NOT.LAROME)) THEN
   IF(YRCONV%LACTIVE) THEN
     YDVARS%RCONV%T1(:,:) = ZQRCONV(:,:)
   ENDIF
   IF(YSCONV%LACTIVE) THEN
     YDVARS%SCONV%T1(:,:) = ZQSCONV(:,:)
   ENDIF
ENDIF



!-------------------------------------------------
! Extract Single Column Model profiles from 3D run or 
! write LFA file for MUSC (1D model)
!-------------------------------------------------
IF(LGSCM.OR.LMUSCLFA) THEN
  IF (LAROME) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        YDCPG_MISC%NEB(JROF,JLEV)=YDVARS%A%T1(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
  CALL WRITEPHYSIO(YDGEOMETRY,YDSURF,YDDPHY,YDRIP,YDMODEL%YRML_PHY_MF,&
       & KEND,&
       & KST, KGL1, KGL2, KSTGLO,&
       & NSTEP  , NTSSG  , YSP_SBD%NLEVS   ,&
       & YDGSGEOM%GELAM,YDGSGEOM%GEMU,&
       & YDGSGEOM%GM, ZMU0,YDOROG%OROG,YDCPG_DYN0%OROGL,YDCPG_DYN0%OROGM,YDGSGEOM%RCORI,YDCSGEOM%RATATH,YDCSGEOM%RATATX,&
       & YDCPG_DYN0%PHI  , YDCPG_PHY0%PREHYD  , YDCPG_DYN0%PHIF , YDCPG_PHY0%PREHYDF , YDCPG_PHY0%XYB%ALPH, YDCPG_PHY0%XYB%DELP,&
       & YDCPG_PHY0%XYB%LNPR, YDCPG_PHY0%XYB%RDELP,&
       & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PARG, YDMF_PHYS_SURF%GSD_VV%PSAB,&
       & YDMF_PHYS_SURF%GSD_VV%PD2, YDMF_PHYS_SURF%GSD_VV%PIVEG, YDMF_PHYS_SURF%GSD_VV%PLAI, YDMF_PHYS%OUT%CT,&
       & YDMF_PHYS_SURF%GSD_VV%PALV, YDMF_PHYS%OUT%ALB, YDMF_PHYS_SURF%GSD_VF%PALBF, YDMF_PHYS_SURF%GSD_VF%PALBSF, YDMF_PHYS_SURF%GSP_SG%PT_T1,&
       & YDMF_PHYS_SURF%GSD_VF%PEMISF,&
       & YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H, YDMF_PHYS_SURF%GSD_VF%PZ0RLF,&
       & YDMF_PHYS_SURF%GSD_VF%PGETRL, YDMF_PHYS_SURF%GSD_VF%PVRLAN, YDMF_PHYS_SURF%GSD_VF%PVRLDI, YDMF_PHYS_SURF%GSD_VV%PRSMIN,&
       & YDMF_PHYS_SURF%GSD_VF%PVEG, YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VA%PSEA, YDMF_PHYS_SURF%GSD_VA%PLAN,&
       & YDMF_PHYS_SURF%GSD_VA%PSOO, YDMF_PHYS_SURF%GSD_VA%PDES, YDMF_PHYS_SURF%GSD_VC%PGROUP,&
       & YDVARS%SP%T0,YDVARS%SP%DL,YDVARS%SP%DM,&
       & YDVARS%T%T0,YDVARS%T%DL,YDVARS%T%DM,&
       & YDVARS%Q%T0,YDVARS%Q%DL,YDVARS%Q%DM,&
       & YDVARS%I%T0,YDVARS%L%T0, YDVARS%S%T0,YDVARS%R%T0,YDVARS%G%T0,YDVARS%TKE%T0,&
       & YDVARS%EFB1%T0,YDVARS%EFB2%T0,YDVARS%EFB3%T0,&
       & YDCPG_DYN0%RCP%CP, YDCPG_DYN0%RCP%R,&
       & YDVARS%U%T0,YDVARS%V%T0,YDVARS%VOR%T0,YDVARS%DIV%T0, ZCVGQ, ZLCVQ,&
       & YDMF_PHYS_SURF%GSP_RR%PT_T9, YDMF_PHYS_SURF%GSP_SB%PT_T9, YDMF_PHYS_SURF%GSP_RR%PFC_T9, YDMF_PHYS_SURF%GSP_RR%PW_T9,&
       & YDMF_PHYS_SURF%GSP_RR%PIC_T9, YDMF_PHYS_SURF%GSP_SB%PQ_T9, YDMF_PHYS_SURF%GSP_SB%PTL_T9, YDMF_PHYS_SURF%GSP_SG%PF_T9,&
       & YDMF_PHYS_SURF%GSP_SG%PA_T0, YDMF_PHYS_SURF%GSP_SG%PR_T0, YDCPG_DYN0%CTY%VVEL(:,1:),&
       & YDMF_PHYS%RAD%EMTD , YDMF_PHYS%RAD%EMTU,  YDMF_PHYS%RAD%TRSW ,YDGSGEOM%GECLO,YDGSGEOM%GESLO,&
       & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL,&
       & YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%FCQNNG,&
       & YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLSG, YDMF_PHYS%OUT%FPLSH,&
       & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN,&
       & ZFTKE, ZFEFB1, ZFEFB2, ZFEFB3, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRSOC, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%FRTHC, YDMF_PHYS%OUT%FRSOLU,&
       & YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, YDMF_PHYS%OUT%STRDU, YDMF_PHYS%OUT%STRDV, YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV, YDMF_PHYS%OUT%STRMU, YDMF_PHYS%OUT%STRMV,&
       & YDMF_PHYS%OUT%FRMH, ZFRMQ, YDMF_PHYS%OUT%FCHOZ, YDCPG_MISC%NEB, YDCPG_MISC%QICE, YDCPG_MISC%QLI, YDCPG_MISC%RH,&
       & YDMF_PHYS%OUT%FCS, YDMF_PHYS%OUT%FCLL, YDMF_PHYS%OUT%FCLN, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, ZFEVI, YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FTR, YDMF_PHYS%OUT%FLWSP, YDMF_PHYS%OUT%FONTE,&
       & YDMF_PHYS%OUT%FGEL, YDMF_PHYS%OUT%FGELS, YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FRSODS, YDMF_PHYS%OUT%FRSOPS, YDMF_PHYS%OUT%FRSOPT, YDMF_PHYS%OUT%FRTHDS,&
       & ZCD, ZCDN, ZCH, ZC1, ZC2, ZEMIS, YDMF_PHYS%OUT%GZ0   , YDMF_PHYS%OUT%GZ0H  , ZNEIJ  , ZVEG,&
       & ZCPS, ZLHS, ZRS, ZLH, ZLSCPE, ZQSAT, ZQW, ZTW,&
       & YDCPG_MISC%QS, ZQSATS, YDMF_PHYS%OUT%RUISL, YDMF_PHYS%OUT%RUISP, YDMF_PHYS%OUT%RUISS,&
       & YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%QCLS, YDMF_PHYS%OUT%RHCLS,&
       & YDCPG_MISC%CLCT, YDMF_PHYS%OUT%CLCH, YDMF_PHYS%OUT%CLCM, YDMF_PHYS%OUT%CLCL, YDMF_PHYS%OUT%CLCC,&
       & YDMF_PHYS%OUT%CAPE  , YDMF_PHYS%OUT%CTOP, ICLPH  , YDMF_PHYS%OUT%CLPH  , YDMF_PHYS%OUT%UGST  , YDMF_PHYS%OUT%VGST,&
       & ZFPLCH, ZFPLSH,&
       & YDMF_PHYS_SURF%GSD_VH%PTCCH,  YDMF_PHYS_SURF%GSD_VH%PSCCH, YDMF_PHYS_SURF%GSD_VH%PBCCH, YDMF_PHYS_SURF%GSD_VH%PPBLH )
ENDIF

IF (LEDR) THEN
  YDMF_PHYS_SURF%GSD_DI%PXEDR(:,:)=ZEDR(:,:)
ENDIF

IF (LDPRECIPS) THEN
  YDMF_PHYS_SURF%GSD_XP%PPRECIP(KST:KEND,NDTPRECCUR)=ZDPRECIPS(KST:KEND,NDTPRECCUR)
ENDIF

IF (LDPRECIPS2) THEN
  YDMF_PHYS_SURF%GSD_XP2%PPRECIP2(KST:KEND,NDTPRECCUR2)=ZDPRECIPS2(KST:KEND,NDTPRECCUR2)
ENDIF

! Restore Tt and grad(Tt) for NHQE model.
IF (LNHQE) THEN
  ! At instant t (with the derivatives):
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      YDVARS%T%T0(JROF,JLEV)=ZTT0_SAVE(JROF,JLEV)
      YDVARS%T%DL(JROF,JLEV)=ZTT0L_SAVE(JROF,JLEV)
      YDVARS%T%DM(JROF,JLEV)=ZTT0M_SAVE(JROF,JLEV)
    ENDDO
  ENDDO
  ! At instant t-dt for leap-frog advections (without the derivatives):
  IF (.NOT.LTWOTL) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        YDVARS%T%T9(JROF,JLEV)=ZTT9_SAVE(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!       6. destructor for procset
IF (LINTFLEX) CALL CLEANINTPROCSET(YLPROCSET)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('MF_PHYS',1,ZHOOK_HANDLE)
END SUBROUTINE MF_PHYS
