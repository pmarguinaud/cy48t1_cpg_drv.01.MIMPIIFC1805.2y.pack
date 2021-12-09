#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE MF_PHYS_PAR(YDGEOMETRY, YDCPG_DIM, YDCPG_MISC, YDCPG_PHY0, YDCPG_PHY9, YDMF_PHYS, YDMF_PHYS_TMP, YDAPLPAR_TMP, &
& YDCPG_DYN0, YDCPG_DYN9, YDMF_PHYS_SURF, YDVARS, YDGMV, YDSURF, YDCFU, YDXFU, &
& YDMODEL, LDCONFX, PDTPHY, &
& PGFL, PKOZO, PGP2DSDT, PB1, PB2, PGMVT1, PGFLT1, PGPAR, PTRAJ_PHYS, YDDDH,   &
& PFTCNS)

!**** *APLPAR * - APPEL DES PARAMETRISATIONS PHYSIQUES.

!     Sujet.
!     ------
!     - APPEL DES SOUS-PROGRAMMES DE PARAMETRISATION
!       INTERFACE AVEC LES PARAMETRISATIONS PHYSIQUES (IALPP).
!     - CALL THE SUBROUTINES OF THE E.C.M.W.F. PHYSICS PACKAGE.

!**   Interface.
!     ----------
!        *CALL* *APLPAR*

!-----------------------------------------------------------------------

! -   ARGUMENTS D'ENTREE.
! -   INPUT ARGUMENTS.
!     -------------------

! - NOM DES PARAMETRES DE DIMENSIONNEMENT DE LA PHYSIQUE.
! - DIMENSIONS.

! KBL  : NUMERO DE BLOC NPROMA
! KBL  : NPROMA-PACKETS NUMBER
! KGPCOMP : NOMBRE TOTAL DE POINTS DE GRILLE SUR LE DOMAINE
! KGPCOMP : TOTAL GRID POINTS NUMBER IN THE DOMAIN
! KIDIA, KFDIA : BORNES BOUCLES HORIZONTALES   (IST,IEND DANS CPG).
! KIDIA, KFDIA : START/END OF HORIZONTAL LOOP  (IST,IEND IN *CPG*).
! KLON : DIMENSION HORIZONTALE                 (NPROMA DANS CPG).
! KLON : HORIZONTAL DIMENSION                  (NPROMA IN *CPG*).
! KLON: IDEM BUT FOR ARRAYS USED SOLELY BY DMN PHYSICS
! KLEV : FIN BOUCLE VERTICE ET DIMENSION VERTICALE (NFLEVG DANS CPG).
! KLEV : END OF VERTICAL LOOP AND VERTICAL DIMENSION(NFLEVG IN *CPG*).
! KSTGLO     : global offset.
! KSGST      : NOMBRE DE TEMPERATURES ET DE FLUX DE SURFACE SOUS-MAILLE
!                     (NTSSG DANS CPG)
! KSGST      : NUMBER OF SUBGRID SURFACE TEMPERATURES AND FLUXES
!                     (NTSSG IN *CPG*)
! KNFRRC     : FREQUENCE DE CALCUL DU RAYONNEMENT CLEAR SKY
! KNFRRC     : FREQUENCY FOR CLEAR SKY RADIATION CALCULATION
! PDT        : TIME STEP (in s)

! LOGICAL
! LDXFUMSE   : T if CDCONF=X in order not to increment surfex timer in that case

! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
!   CATEGORIE).

! - 2D (0:KLEV) .


! - 2D (1:KLEV) .

! PGFL       : GFL FIELDS
! PCVGQ      : CONVERGENCE D'HUMIDITE (CONDITION DE "KUO").
! PCVGQ      : CONVERGENCE OF HUMIDITY ("KUO" CONDITION).
! PKOZO      : CHAMPS POUR LA PHOTOCHIMIE DE L'OZONE (KVCLIS CHAMPS).
! PKOZO      : FIELDS FOR PHOTOCHEMISTERY OF OZONE   (KVCLIS FIELDS).

! - 1D (PROGNOSTIQUE) .
! - 1D (PROGNOSTIC QUANTITIES) .


! - 1D (GEOGRAPHIQUE) .
! - 1D (GEOGRAPHICAL DISTRIBUTION) .

! PALV       : ALBEDO DE LA VEGETATION.
! PALV       : VEGETATION ALBEDO.
! PALBF      : ALBEDO FIXE (SOL SANS NEIGE).
! PALBF      : BACKGROUND SURFACE SHORTWAVE ALBEDO (SNOW-FREE).
! PALBSF     : ALBEDO FIXE DU SOL (SOL SANS NEIGE).
! PALBSF     : BACKGROUND SOIL SHORTWAVE ALBEDO (SNOW-FREE).
! PARG       : POURCENTAGE D'ARGILE DANS LE SOL.
! PARG       : SILT PERCENTAGE WITHIN THE SOIL.
! PD2        : EPAISSEUR DU RESERVOIR PROFOND.
! PD2        : DEEP-LAYER OF THE OF THE PROFOUND WATER-TANK.
! PEMISF     : EMISSIVITE FIXE (SOL SANS NEIGE).
! PEMISF     : BACKGROUND SURFACE LONGWAVE EMISSIVITY (SNOW-FREE).
! PGETRL     : ECART TYPE DE L'OROGRAPHIE (EN J/KG).
! PGETRL     : STANDARD DEVIATION OF OROGRAPHY  (UNIT: J/KG).
! PLSM       : INDICE TERRE/MER.
! PLSM       : LAND/SEA MASK.
! PIVEG      : TYPE DE VEGETATION.
! PIVEG      : TYPE OF VEGETATION.
! PLAI       : INDICE FOLIAIRE.
! PLAI       : FOLIAIRE INDICE.
! PRCORI     : FACTEUR DE CORIOLIS.
! PRCORI     : CORIOLIS FACTOR.
! PRSMIN     : RESISTANCE STOMATIQUE MINIMALE.
! PRSMIN     : STOMATAL MINIMAL RESISTANCE.
! PSAB       : POURCENTAGE DE SABLE DANS LE SOL.
! PSAB       : PERCENTAGE OF SAND WITHIN THE SOIL.
! PVRLAN     : ANISOTROPIE DU RELIEF SOUS MAILLE.
! PVRLDI     : ANGLE DE LA DIRECTION DU RELIEF AVEC L'AXE DES X.
!              DUE AU RELIEF.

! - 1D (DIAGNOSTIQUE) .
! - 1D (DIAGNOSTIC QUANTITIES).

! PHV        : RESISTANCE A L'EVAPOTRANSPIRATION.
! PHV        : RESISTANCE TO EVAPOTRANSPIRATION.
! PMU0       : COSINUS LOCAL INSTANTANE DE L'ANGLE ZENITHAL SOLAIRE.
! PMU0LU     : COSINUS LOCAL INSTANTANE DE L'ANGLE ZENITHAL LUNAIRE.
! PMU0       : LOCAL COSINE OF INSTANTANEOUS SOLAR ZENITH ANGLE.
! PMU0LU     : LOCAL COSINE OF INSTANTANEOUS LUNAR ZENITH ANGLE.
! PMU0M      : COSINUS LOCAL MOYEN DE L'ANGLE ZENITHAL.
! PMU0M      : LOCAL COSINE OF AVERAGED SOLAR ZENITH ANGLE.
! PMU0N      : COSINUS LOCAL AU PAS DE TEMPS SUIVANT DE L'ANGLE ZENITHAL SOLAIRE.
! PMU0N      : NEXT TIME STEP COSINUS LOCAL INSTANTANE DE L'ANGLE ZENITHAL SOLAIRE.
! PFPLCH     : CHAMP "HISTORIQUE" POUR LES PRECIPITATIONS CONVECTIVES.
! PFPLCH     : HISTORICAL FIELD FOR CONVECTIVE PRECIPITATIONS.
! PFPLSH     : CHAMP "HISTORIQUE" POUR LES PRECIPITATIONS STRATIFORMES.
! PFPLSH     : HISTORICAL FIELD FOR STRATIFORM PRECIPITATIONS.
! PGPAR       : BUFFER FOR 2D FIELDS - CONTAINS PRECIP, ALBEDO, EMISS, TS
!             : SURFACE FLUXES

! PMMU0      : MEAN SOLAR ANGLE FOR GIVEN GRID POINT FOR SIMPL.
!              RADIATION SCHEME.
!              SIMPLIFIED RADIATION SCHEME.

! PTCCH      : PSEUDO-HISTORICAL ARRAY FOR TOTAL CONVECTIVE CLOUDINESS (1D).
! PSCCH      : PSEUDO-HISTORICAL ARRAY FOR CONVECTIVE CLOUD SUMMIT (1D).
! PBCCH      : PSEUDO-HISTORICAL ARRAY FOR CONVECTIVE CLOUD BASE (1D).
! PPBLH      : PSEUDO-HISTORICAL ARRAY FOR PBL HEIGHT (1D).
! PQSH       : PSEUDO-HISTORICAL ARRAY FOR SURFACE MOISTURE (1D).
! PUDGRO     : PSEUDO-HISTORICAL ARRAY FOR UPDRAFT RISING TOP (1D) - ACCSU
!              FRACTION OF PATH COMPLETED ABOVE LAST ACTIVE LEVEL AT T9

! - CONSTANTES

! POMPAC    : Horizontally variable Curtis matrix for thermal rad. (clear sky)
! PAC       : Horizontally cst Curtis matrix for thermal rad. (clear sky)



!-----------------------------------------------------------------------

! -   INPUT/OUTPUT ARGUMENTS.
!     -----------------------

! - 2D (ACRANEB2 INTERMITTENCY STORAGE)

!             LINEAR T_e CORRECTION
!             LINEAR T_e CORRECTION

! - 1D (ACRANEB2 INTERMITTENCY STORAGE)


!-----------------------------------------------------------------------

! -   ARGUMENTS DE SORTIE.
! -   OUTPUT ARGUMENTS.
!     --------------------

! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
!   CATEGORIE).

! - 2D (0:KLEV) .

! PDIFSV     : FLUX TURBULENT DES SCALAIRES PASSIFS.
! PDIFSV     : TURBULENT FLUX OF PASSIVE SCALAR.
! PFTKE      : FLUX DE TKE.
! PFTKE      : TKE FLUX.
! PFRMQ      : FLUX MESOSPHERIQUE D'humidite.
! PFRMQ      : MESOSPHERIC humidity flux

! - 2D (1:KLEV) .

! PLH        : CHALEUR LATENTE A LA TEMPERATURE DE L'AIR.
! PLH        : LATENT HEAT AT AIR TEMPERATURE.
! PLSCPE     : RAPPORT EFECTIF DES L ET CP EN CONDENSATION/EVAPORATION.
! PLSCPE     : EFFECTIVE RATIO OF L AND CP FOR CONDENSATION/EVAPORATION.
! PQSAT      : HUMIDITE SPECIFIQUE DE SATURATION.
! PQSAT      : SPECIFIC HUMIDITY AT SATURATION.
! PQW        : HUMIDITE SPECIFIQUE DU THERMOMETRE MOUILLE.
! PQW        : SPECIFIC HUMIDITY OF THE WET THERMOMETER.
! PTW        : TEMPERATURE DU THERMOMETRE MOUILLE.
! PTW        : TEMPERATURE OF THE WET THERMOMETER.
! PTENDPTKE  : INCREMENT DE LA ENERGIE CINETIQUE TURBULENTE (schema pTKE)
! PTENDPTKE  : TKE INCREMENT (USED BY PSEUDOPROGNOSTIC TKE SCHEME)
! PKUROV_H   : Coefficient d'echange horizontal de u et v
! PKUROV_H   : Horizontal exchange coefficient for momentum
! PKTROV_H   : Coefficient d'echange horizontal de T et q
! PKTROV_H   : Horizontal exchange coefficient for heat

! - 2D (CLOUD AND RADIATION I/O FOR ECMWF PHYSICS)

! ZTENT      : TENDENCY OF TEMPERATURE.

! - 2D (0:1)

! - 1D (DIAGNOSTIQUE) .

! PCD        : COEFFICIENT D'ECHANGE EN SURFACE POUR U ET V.
! PCD        : EXCHANGE COEFFICIENT AT SURFACE LEVEL FOR U AND V.
! PCDN       : COEFFICIENT NEUTRE D'ECHANGE EN SURFACE.
! PCDN       : EXCHANGE COEFF. AT SURFACE LEVEL IN NEUTRAL CONDITIONS.
! PCH        : COEFFICIENT D'ECHANGE EN SURFACE POUR T ET Q.
! PCH        : EXCHANGE COEFFICIENT AT SURFACE LEVEL FOR T AND Q.
! PCPS       : CHALEUR MASSIQUE DE L'AIR EN SURFACE.
! PC1        : COEFF. HYDRIQUE REPRESENTANT L'INTENSITE AVEC LAQUELLE
!              LES FLUX DE SURFACE PARTICIPENT A L'EVOLUTION DE WS.
! PC1        : HYDROLOGICAL COEFF. SHOWING THE CONTRIBUTION OF SURFACE
!              FLUXES IN THE WS EVOLUTION.
! PC2        : COEFF. HYDRIQUE TRADUISANT LA RAPIDITE DES TRANSFERTS
!              D'EAU ENTRE LES DEUX RESERVOIRS.
! PC2        : HYDROLOGICAL COEFFICIENT SHOWING THE QUICKNESS OF WATER
!              TRANSFERS BETWEEN BOTH TANKS.
! PEMIS      : EMISSIVITE DE SURFACE COURANTE.
! PEMIS      : MODEL SURFACE LONGWAVE EMISSIVITY.
! PFEVI      : FLUX DE VAPEUR D'EAU SUR SOL GELE.
! PFEVI      : WATER VAPOUR FLUX OVER FROZEN SOIL.
! PLHS       : CHALEUR LATENTE EN SURFACE.
! PNEIJ      : PROPORTION DE SOL ENNEIGE.
! PNEIJ      : FRACTION OF SOIL COVERED BY SNOW.
! PVEG       : FRACTION DE VEGETATION APPARENTE.
! PVEG       : FRACTIONAL COVER BY APPARENT VEGETATION.
! PQSATS     : HUMIDITE SPECIFIQUE DE SATURATION EN SURFACE.
! PQSATS     : SATURATED SPECIFIC HUMIDITY AT SURFACE LEVEL.
! PRS        : CONSTANTE DES GAZ DE L'AIR EN SURFACE.
! PCLCC      : SORTIE DIAGNOSTIQUE DE LA NEBULOSITE CONVECTIVE.
! PCLCC      : CONVECTIVE CLOUD COVER (DIAGNOSTIC).
! PDPRECIPS  : PRECIPITATION TYPE
! PDPRECIPS2 : PRECIPITATION TYPE FOR 2NDE PERIOD
! PTENDEXT_DEP:WET TENDENCY OF PASSIVES SCALAIRE
! PTENDEXT_DEP:TENDANCE HUMIDE DES SCALAR PASSIFS
! - INPUT/OUTPUT 1D
! YDDDH      : DDH superstructure
!-----------------------------------------------------------------------

! -   ARGUMENTS IMPLICITES.
! -   IMPLICIT ARGUMENTS.
!     ---------------------

!-----------------------------------------------------------------------

!     Externes.
!     ---------

!     Methode.
!     --------
!     - TERMINE LES INITIALISATIONS.
!     - APPELLE LES SS-PRGMS TAMPONS SUIVANT LA LOGIQUE TROUVEE
!        DANS /YOMPHY/. EUX MEMES VONT DECLARER LES TABLEAUX DE TRAVAIL
!        ET APPELER LES PARAMETRISATIONS ELLES MEMES.
!     - FINISH UP THE INITIALIZATION.
!     - CALL THE BUFFER SUBROUTINES FOLLOWING /YOEPHY/ REQUIREMENTS
!        WHICH IN TURN CALL THE ACTUAL PHYSICS SUBROUTINES
!        (THIS LAST POINT NOT PARTIALLY DONE)

!     Auteur.
!     -------
!     90-09-28: A. Joly, *CNRM*.

!     Modifications.
!     --------------
!     2007-02-01 M.Janousek : Introduction of 3MT routines
!     2007-02-19 R.Brozkova : Cleaning obsolet features (LSRCON, LSRCONT, LNEBT,
!                             pre-ISBA, modularisation and racionalisation.
!     2007-04-17 S.Ivatek-S : Over dimensioning of PGPAR (KLON,NGPAR+1) is used
!                             boundary checking bf
!     2007-05-10 E. Bazile  : Introduction of the AROME shallow convection (LCVPPKF)
!     2007-03-21 A. Alias   : Modifications for SURFEX (IGFL_EXT)
!     2007-05-07 F. Bouyssel: Several modifications for SURFEX
!     2007-05-07 F. Bouyssel: New argument in ACCOEFK
!     2007-06-27 A. Alias   : Use NGFL_EXT instead of IGFL_EXT
!     2008-02-18 F. Bouyssel: New acdifv1 & acdifv2 & arp_ground_param
!     2008-02-21 E. Bazile  : Cleaning for the call of the AROME shallow convection (LCVPPKF)
!     4-Mar-2008 Y. Seity : Cleaning IR and WV similated sat radiances
!                            (replaced by Fullpos Calculations)
!     2008-03-14 Y. Bouteloup: Store diffusion coefficients from non-linear model
!     2007-10-11 A. Alias   : New Call to ACHMT/ACNEBR/ACPBLH (P. Marquet/JF. Gueremy)
!     2008-02-01 P. Marquet : modify ZALBD/ZALBP and PFRSODS if LRAYFM15 (idem V4)
!     2008-03-26 F. Bouyssel: Intrduction of LACDIFUS
!     2008-04-28 E. Bazile  : Introduction of ZPROTH_CVPP for the TKE scheme
!     2008-05-09, J.F. Gueremy : Flux MEMO sur mer (ACFLUSO/LFLUSO) +
!            and  P. Marquet   : ZCEROV as new argument in ACDIFUS
!     2008-06-01 F. Bouyssel: Interface of radozc (ECMWF ozone)
!     2008-09-02 F. Vana  : Better split of ACPTKE and ACDIFV1 code
!     2008-10-01 F. Bouyssel: Call of radozcmf instead of radozc (ECMWF ozone)
!     2008-10-05 E. Bazile : Computation of the PBL height from the TKE.
!     03-Oct-2008 J. Masek    parameters for NER statistical model via namelist
!     2009-Jan-21 F. Vana : new mixing lengths for pTKE + few fixes
!     2008-11    C. Payan: Neutral Wind (new arg in the call of ACHMT)
!     2008-11-15 F. Bouyssel: Correction of negative humidities
!     2009-05-01 Y. Bouteloup : Store radiative cloud water and ice in GFL (YIRAD and YLRAD)
!     2009-05-25 F. Bouyssel: Cleaning
!     K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!     2009-08-07 A. Alias   : LCALLSFX introduced to call only once SURFEX at NSTEP=0
!                             Computation of ZRTI (INVERSE DE R*T) added (after ACHMTLS)-not done
!                             Negatives humidity correction PFCQVNG in acdifus (J.F. Gueremy)
!                             add ZAESUL/ZAEVOL to CALL RADAER
!     2009-09-21  D. Banciu: complete the cascade within 3MT frame;
!            prepare the environment for Rash and Kristjansson condensation (RK) scheme
!            remove some arguments of ACPUM and ACUPD (LUDEN option was removed)
!     12-Oct-2009 F. Vana : optimization + update of mixing lengths for p/eTKE
!     2009-10-15 Y. Bouteloup : Store radiative cloud water and ice in GFL (YIRAD and YLRAD)
!     2009-10-23 O.Riviere Intro. of LGWDSPNL for GWD in simpl. phys.
!     2010-01-21 S. Riette: PU, PV and ZDEPTH_HEIGHT in 3D for ARO_GROUND_DIAG
!     2010-05-11 F. Bouyssel : Use of PINDX, PINDY
!     2010-06-20 Y. Seity : Use AROCLDIA to compute PBLH
!     2010-10    A. Alias Compute Sunshine duration
!     2010-10    A. Alias modify ZALBD/ZALBP and PEMIS if LRAYFM for CLIMAT (JF Gueremy)
!     2010-08-10 P.marguinaud : Cleaning
!     2011-01-10 F. Bouyssel: Intro. of LADJCLD and some cleaning.
!     2010-12-01 E. Bazile: TKE1 for AROCLDIA and contributions terms of the
!          TKE equations for DDH.
!     2010-12 S. Riette: aro_ground_diag interface modified to add snow cover
!     2010-12    B. Decharme  : modify the radiative coupling with surfex (SW per band in ACRADIN and RADHEAT)
!     2011-02    A. Alias     : Computation of ZRTI (INVERSE DE R*T) added (after ACHMTLS)
!     2011-02    A. Voldoire : add ZAERINDS to CALL RADAER and ACRADIN
!                              for sulfate indirect effect computation
!     L. Bengtsson-Sedlar & F. Vana 18-Feb-2011 : CA scheme for convection
!     I. Bastak-Duran, F. Vana & R. Brozkova  16-Mar-2011: TOUCANS, version 0
!     2011-02-01 M. Mokhtari: Several modifications for aplpar and introduction of the key LMDUST
!                             (treatment of the desert aerosols)
!     2011-02-24 Y. Bouteloup : EDKF + Surface forcing for MUSC
!     2011-03-26 F. Bouyssel: Intro. of PSPSG (snow cover with surfex)
!     2011-09-07 J.M. Piriou: PCMT convection scheme.
!     2011-11-17 J.F. Gueremy: ZQLI_CVP diagnostic convective water content
!     2011-06: M. Jerczynski - some cleaning to meet norms
!     2011-11-29 K-I. Ivarsson, L. Bengtsson: RK-scheme modifications   
!     26-Jan-2012: F. Vana + I. Bastak-Duran - TOUCANS update + bugfixes 
!     2012-04-24 F. Bouyssel: Bug correction on surface water fluxes with surfex
!     2012-06-09 M. Mile: Bug correction for undefined z0;z0h at 0th step CALL ARO_GROUND_DIAG_Z0
!     2012-09-11 : P.Marguinaud : Add control threshold for
!     2013-06-17 J.M. Piriou: evaporation for PCMT scheme.
!     T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     2013-11-08 Y. Bouteloup New version of ACDIFV1 and ACDIFV2 for "full implicit PMMC09 scheme"
!     F. Vana  28-Nov-2013 : Redesigned trajectory handling.
!     2013-11, J. Masek: Introduction of ACRANEB2 scheme, externalized
!                        computation of direct albedo for ACRANEB/ACRANEB2.
!                        Phasing to cy40t1.
!     K. Yessad (July 2014): Move some variables.
!     2014-09, C. Wastl: Adaptations for orographic shadowing
!     2014-10, R. Brozkova: phasing TOUCANS.
!     2016-03, E. Bazile: phasing MUSC for surf_ideal_flux
!     2016-03, L. Gerard: LNSDO AND LCVCSD
!     2016-04, J. Masek: Exponential-random cloud overlap with variable
!                        decorrelation depth.
!     2016-09, J. Masek: Proper diagnostics of sunshine duration in ACRANEB2.
!     2016-09, M. Mokhtari & A. Ambar: preliminary calculation for passive scalar
!     2016-10, P. Marguinaud : Port to single precision
!     K. Yessad (Dec 2016): Prune obsolete options.
!     2016-06, F.Taillefer: add of aro_ground_diag_2isba call
!     2017-09, Y.Bouteloup: Phased Francoise's modification on cy45
!     2017-09, J. Masek: Fix for protected convective cloudiness,
!                        shifted dimensioning of PGMU0.
!     R. El Khatib 05-Feb-2018 fix bounds violations
!     2018-09, F. Duruisseau: Add PQRCONV1 and PQSCONV1 out (BAYRAD)
!     2018-09, D. St-Martin: Add non-orographic GWD scheme (ACNORGWD)
!     2018-09, R. Roehrig: add ACTKE input/output (ZQLC/ZQIC and ZKQROV/ZKQLROV) (from JF Guérémy)
!     2018-09, M. Michou: Add call to chem_main to activate ARPEGE-Climat chemistry scheme
!     2018-09, R. Brozkova: Fixes in thermodynamic adjustment - deep convective
!                           condensates protection. Passing of diagnostic hail.
!     2018-09, J. Masek: Calculation of snow fractions over bare ground and
!                        vegetation moved to ACSOL (case LVGSN=T). Coding of
!                        ALARO-1 fixes for LZ0HSREL=T. Diagnostics of global
!                        normal irradiance and mean radiant temperature.
!     2018-11, J.M. Piriou: correct 2010 historical bug about adding cloud sedimentation to resolved surface precipitation.
!     R. Hogan     24-Jan-2019 Removed radiation scheme from cycle 15
!     R. El Khatib 30-Apr-2019 fix uninitialized variable
!     2018-10, I. Etchevers : add Visibilities
!     2019-01, I. Etchevers, Y. Seity : add Precipitation Type
!     2019-05, J.M. Piriou: LCVRESDYN + LADJCLD.
!     2019-09, M. Hrastinski: Dataflow for TKE and TTE terms in ALARO DDH (PFTCNS).
!     2019-09, L. Gerard: Modified call to ACNSDO.
!     2019-09, R. Brozkova: Introduction of new NDIFFNEB options.
!     2019-09, J. Masek: Introduction of ETKE_MIN, efficient ACRANEB2 clearsky
!                        computations.
!    2019-10, I. Etchevers : Visibilities in ACVISIH, AROCLDIA=>ACCLDIA
!    2019-10, Y.Bouteloup : New anti-GPS in accvimp.F90 
!    2019-10, Y.Bouteloup and M. Bouzghaiam : Radiation modifications. Remove of FMR15, remove acradin.F90 direct
!                   call to recmwf.F90 and add interface to ecrad (in recmwf !)
!    2020-07, J.M. Piriou and O. Jaron: lightning flash density.
!    2020-10, J. Masek : modified call to ACCLDIA
!    2020-10, M. Hrastinski: Reorganized computation of the moist gustiness
!             correction. Modified call of ACMRIP and ACMIXELEN subroutines.
!    2020-11, Y.Bouteloup : Interface to IFS deep convection scheme under LCVTDK key
!    2020-12, U.Andrae : Introduce SPP for HARMONIE-AROME
! End Modifications
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!     ******************************************************************
!     ****** IDIOSYNCRASIES *** IDIOSYNCRASIES *** IDIOSYNCRASIES ******
!     ******************************************************************
!     ***  HEALTH WARNING:                                           ***
!     ***  ===============                                           ***
!     ***  NOTE THAT WITHIN THE E.C.M.W.F. PHYSICS HALF-LEVELS       ***
!     ***  ARE INDEXED FROM 1 TO NFLEVG+1 WHILE THEY ARE BETWEEN     ***
!     ***  0 AND NFLEVG IN THE REST OF THE MODEL. THE CHANGE IS TAKEN***
!     ***  CARE OF IN THE CALL TO THE VARIOUS SUBROUTINES OF THE     ***
!     ***  PHYSICS PACKAGE                                           ***
!     ***                                                            ***
!     ***    THIS IS SUPPOSED TO BE A "TEMPORARY" FEATURE TO BE      ***
!     ***    STRAIGHTENED OUT IN THE "NEAR" FUTURE                   ***
!     ******************************************************************

!     ------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE MF_PHYS_TYPE_MOD   , ONLY : MF_PHYS_TYPE, MF_PHYS_TMP_TYPE, APLPAR_TMP_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_MISC_TYPE, CPG_DYN_TYPE
USE CPG_DIM_TYPE_MOD   , ONLY : CPG_DIM_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD,ONLY : MF_PHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMXFU             , ONLY : TXFU
USE TYPE_MODEL         , ONLY : MODEL
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK    ,DR_HOOK
USE YOMMP0             , ONLY : NPRINTLEV
USE YOMVERT            , ONLY : VP00
USE YOMCST             , ONLY : RG       ,RSIGMA   ,RV       ,RD       ,&
                              & RCPV     ,RETV     ,RCW      ,RCS      ,RLVTT    ,&
                              & RLSTT    ,RTT      ,RALPW    ,RBETW    ,RGAMW    ,&
                              & RALPS    ,RBETS    ,RGAMS    ,RALPD    ,RBETD    ,&
                              & RGAMD    ,RCPD     ,RATM     ,RKAPPA   ,RLMLT
USE YOMCT0             , ONLY : LCALLSFX ,LSFORCS, LELAM
USE YOMDYNA            , ONLY : L3DTURB
USE YOMRIP0            , ONLY : NINDAT
USE DDH_MIX            , ONLY : ADD_FIELD_3D ,ADD_FIELD_2D ,NTOTSVAR , & 
                              & NTOTSURF ,NTOTSVFS, NEW_ADD_FIELD_3D, NEW_ADD_FIELD_2D, &
                              & TYP_DDH
USE YOMLUN             , ONLY : NULOUT
USE YOMLSFORC          , ONLY : LMUSCLFA,NMUSCLFA
USE YOMTRAJ            , ONLY : TRAJ_PHYS_TYPE
USE YOMCFU             , ONLY : TCFU !!! for parameters of FLASH
USE SPP_MOD  , ONLY : YSPP, YSPP_CONFIG
USE MF_PHYS_STATE_TYPE_MOD &
                       , ONLY : MF_PHYS_STATE_TYPE


USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE MF_PHYS_TYPE_MOD   , ONLY : MF_PHYS_TYPE, MF_PHYS_TMP_TYPE, APLPAR_TMP_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_DYN_TYPE, CPG_PHY_TYPE, &
                              & CPG_MISC_TYPE
USE CPG_DIM_TYPE_MOD   , ONLY : CPG_DIM_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD,ONLY : MF_PHYS_SURF_TYPE
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
USE MF_PHYS_STATE_TYPE_MOD &
                       , ONLY : MF_PHYS_STATE_TYPE
!     -------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_DIM_TYPE),INTENT(IN)    :: YDCPG_DIM
TYPE(CPG_MISC_TYPE),INTENT(INOUT):: YDCPG_MISC
TYPE(CPG_PHY_TYPE),INTENT(INOUT), TARGET :: YDCPG_PHY0
TYPE(CPG_PHY_TYPE),INTENT(INOUT), TARGET :: YDCPG_PHY9
TYPE(MF_PHYS_TYPE),INTENT(INOUT) :: YDMF_PHYS
TYPE(MF_PHYS_TMP_TYPE), INTENT(INOUT) :: YDMF_PHYS_TMP
TYPE(APLPAR_TMP_TYPE), INTENT(INOUT) :: YDAPLPAR_TMP
TYPE(CPG_DYN_TYPE),INTENT(INOUT), TARGET :: YDCPG_DYN0
TYPE(CPG_DYN_TYPE),INTENT(INOUT), TARGET :: YDCPG_DYN9
TYPE(MF_PHYS_SURF_TYPE),INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
TYPE(TCFU)        ,INTENT(INOUT) :: YDCFU
TYPE(TXFU)        ,INTENT(INOUT) :: YDXFU
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
LOGICAL           ,INTENT(IN)    :: LDCONFX
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTPHY 
REAL(KIND=JPRB)   ,INTENT(INOUT), TARGET :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM) 
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
INTEGER(KIND=JPIM) :: IPQ,IPO3
INTEGER(KIND=JPIM) :: IPGFL(YDMODEL%YRML_GCONF%YGFL%NUMFLDS)

INTEGER(KIND=JPIM) :: INSTEP_DEB,INSTEP_FIN
INTEGER(KIND=JPIM) :: JROF
INTEGER(KIND=JPIM) :: ISLB1U9  ,ISLB1V9  ,ISLB1T9  ,ISLB1GFL9, ISLB1VD9

!     --- UPPER AIR PHYSICAL TENDENCIES.
REAL(KIND=JPRB) :: ZTENDH(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)     ! Enthalpy tendency.
REAL(KIND=JPRB) :: ZTENDQ(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)     ! Moisture tendency.
REAL(KIND=JPRB) :: ZTENDPTKE(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG)  ! Pseudo progn. TKE

! GFL tendencies for APL_AROME (assumes YDMODEL%YRML_GCONF%YGFL%NUMFLDS>=YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
! for now, use Jovi's trick :
REAL(KIND=JPRB) :: ZTENDGFL(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NUMFLDS)   ! GFL tendencies

!     --- UPPER AIR PHYSICAL TENDENCIES FOR AROME.
!       (the previous one are not used in AROME)
REAL(KIND=JPRB) :: ZTENDT (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)        ! temperature tendency
REAL(KIND=JPRB) :: ZTENDD (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)        ! d  tendency
REAL(KIND=JPRB) :: ZTENDEXT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)      ! GFL EXTRA tendency
REAL(KIND=JPRB) :: ZTENDEXT_DEP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)  ! GFL EXTRA tendency
REAL(KIND=JPRB) :: ZDIFEXT(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)     ! Extra-GFL fluxes.

REAL(KIND=JPRB) :: ZTENDU (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! U tendency without deep convection contribution
REAL(KIND=JPRB) :: ZTENDV (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! V tendency without deep convection contribution

!     --- RADIATION COEFFICIENTS FOR SIMPLIFIED PHYSICS IN GRID-POINT ---
REAL(KIND=JPRB) :: ZAC(YDGEOMETRY%YRDIM%NPROMM,(YDGEOMETRY%YRDIMV%NFLEVG+1)*(YDGEOMETRY%YRDIMV%NFLEVG+1))   ! Curtis matrix.
REAL(KIND=JPRB) :: ZAC_HC(YDGEOMETRY%YRDIMV%NFLEVG+1,YDGEOMETRY%YRDIMV%NFLEVG+1)           ! horizontally-constant field for ZAC.



! required for INTFLEX
TYPE(TYPE_INTPROCSET) :: YLPROCSET

! SPP
REAL(KIND=JPRB) :: ZGP2DSPP(YDGEOMETRY%YRDIM%NPROMA,YSPP%N2D)


REAL(KIND=JPRB), POINTER :: ZPTENDEFB11(:,:), ZPTENDEFB21(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDEFB31(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDG1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDICONV1(:,:), ZPTENDI1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDLCONV1(:,:)
REAL(KIND=JPRB), POINTER :: ZP1EZDIAG(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDQ1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDRCONV1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDR1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDSCONV1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDS1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDTKE1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDL1(:,:)

TYPE (MF_PHYS_STATE_TYPE) :: YLMF_PHYS_STATE
LOGICAL :: LLNEW
CHARACTER*16 :: CLNEW

!     ------------------------------------------------------------------
!        ATTENTION SI KVCLIG < 7 LES CHAMPS SUIVANTS NE SONT
!        PAS REELLEMENT ALLOUES EN MEMOIRE.
!*
!     ------------------------------------------------------------------
!     DECLARATION DES TABLEAUX LOCAUX-GLOBAUX DES PARAMETRISATIONS
INTEGER(KIND=JPIM) :: INLAB(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), INLAB_CVPP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZVETAH(0:YDCPG_DIM%KFLEVG),ZNLAB(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZNLABCVP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZVETAF(YDCPG_DIM%KFLEVG)

!     ------------------------------------------------------------------
!     ARE DIMENSIONNED 0:KLEV ONLY IN ORDER TO KEEP IN MIND
!     THAT THEY ARE COMPUTED AT "HALF LEVELS".
!     THEY ARE USED HOWEVER FROM 1 TO KLEV.

REAL(KIND=JPRB) :: ZXTROV(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZXUROV(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZXPTKEROV(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZRRCOR(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZMRIPP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZMRIFPP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZBNEBCVPP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZBNEBQ(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZMN2PP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZMN2_ES(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZMN2_EQ(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZMN2_DS(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZMN2_DQ(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
! ZMRIMC     : M(C) from P. Marquet's moist Ri computation - for TKE correction after TOMs
! ZMRICTERM  : Rv/R.F(C)-1/M(C).T/Tv from P. Marquet's moist Ri computation - for TKE correction after TOMs
REAL(KIND=JPRB) :: ZMRIMC(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZMRICTERM(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTSTAR(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZTSTAR2(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) !diag 
REAL(KIND=JPRB) :: ZTSTARQ(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZTSTAR2Q(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)!diag
REAL(KIND=JPRB) :: ZTAU_TKE(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)!DISSIPATION TIME SCALE TAU  -FOR TOM's CALCULATION
REAL(KIND=JPRB) :: ZF_EPS(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)   !  Conversion function lm-L
REAL(KIND=JPRB) :: ZFUN_TTE(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)   !  Function in computation of tte_tilde
REAL(KIND=JPRB) :: ZKTROV(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZKUROV(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZNBVNO(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZKQROV(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZKQLROV(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZKNROV(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)

REAL(KIND=JPRB) :: ZFHORM(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZFHORH(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! arrays for 3D turb
! ZFMTKE - F_m function for static K_m computation
! ZFHTKE - F_h function for static K_h computation
REAL(KIND=JPRB) :: ZFMTKE(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZFTTKE(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZRHS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
! ZAUTKE - alpha_u for dry AF scheme
! ZATTKE - alpha_theta for dry AF scheme
REAL(KIND=JPRB) :: ZAUTKE(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZATTKE(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
!ZTH_FUN, ZWW_FUN - T_h, A_h, F_ww - stability functions for TOMs par.
REAL(KIND=JPRB) :: ZTH_FUN(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZWW_FUN(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
! ZFMGST - stability function F_m for moist gustiness correction
REAL(KIND=JPRB) :: ZFMGST(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)


!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZATSLC(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZNEBS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZQLIS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZNEBS0(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZQLIS0(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZNEBC0(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)   !Nebulosite convective radiative
REAL(KIND=JPRB) :: ZNEBDIFF(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) !Nebulosite: calcul de la diffusion
REAL(KIND=JPRB) :: ZNEBCH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)   !Nebulosite convective condensation
REAL(KIND=JPRB) :: ZUNEBH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)   !Nebulosite convective histo
REAL(KIND=JPRB) :: ZDETFI(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)   !fraction of instantaneous detrained air
REAL(KIND=JPRB) :: ZFPCOR(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZFHP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZLMT(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZZLMT(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZLMU(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZLMU2(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZLMT2(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! temporary storage of lm,lh
REAL(KIND=JPRB) :: ZLML(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! TKE type mixing length
REAL(KIND=JPRB) :: ZLMLTILD(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! 'STATIC' TKE type mixing length
REAL(KIND=JPRB) :: ZOME(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)   ! updraught envt vert vel*dt
REAL(KIND=JPRB) :: ZFALLR(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! fall velocity of rain
REAL(KIND=JPRB) :: ZFALLS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! fall velocity of snow
REAL(KIND=JPRB) :: ZFALLG(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! fall velocity of graupel
REAL(KIND=JPRB) :: ZICEFR1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)! Resolved Condensate ice fraction
REAL(KIND=JPRB) :: ZRHCRI(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Smith scheme critical RH
REAL(KIND=JPRB) :: ZRHDFDA(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)! RK scheme change in RH over cloud
REAL(KIND=JPRB) :: ZLHS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)   ! Sublimation latent heat
REAL(KIND=JPRB) :: ZLHV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)   ! Evaporation latent heat
REAL(KIND=JPRB) :: ZLH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! Temporar storage for updated PLH
REAL(KIND=JPRB) :: ZLSCPE(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Temporar storage for updated PLSCPE
REAL(KIND=JPRB) :: ZQSAT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)  ! Temporar storage for updated PQSAT
REAL(KIND=JPRB) :: ZQSATS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! QSAT of resolved cond./evap. scheme
REAL(KIND=JPRB) :: ZQW(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! Temporar storage for updated PQW
REAL(KIND=JPRB) :: ZRH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! Temporar storage for updated PRH
REAL(KIND=JPRB) :: ZTW(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! Temporar storage for updated PTW)
REAL(KIND=JPRB) :: ZDQ(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! Saturation departure for a given thermodynamic state
REAL(KIND=JPRB) :: ZDQM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)   ! maximum saturation departure
REAL(KIND=JPRB) :: ZPOID(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)  ! DP/(RG*DT) FOR A GIVEN LEVEL AND A GIVEN TIME STEP.
REAL(KIND=JPRB) :: ZIPOI(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)  ! INVERSE OF ZPOID.

REAL(KIND=JPRB) :: ZQU(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Updraught Specific moisture
REAL(KIND=JPRB) :: ZTU(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Updraught Temperature
REAL(KIND=JPRB) :: ZUU(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Updraught zonal wind
REAL(KIND=JPRB) :: ZVU(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Updraught merid. wind

REAL(KIND=JPRB) :: ZTMIC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Temperature for microphysics
REAL(KIND=JPRB) :: ZQMIC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Specific moisture for microphysics

REAL(KIND=JPRB) :: ZT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)     ! updated temperature T for cascading parameterization
REAL(KIND=JPRB) :: ZTCORR(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! temperature corr. for convective cloud
REAL(KIND=JPRB) :: ZMELNET(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! net melting (-freezing) rate of ice
REAL(KIND=JPRB) :: ZMELGET(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)! net melting (-freezing) rate of graupel
REAL(KIND=JPRB) :: ZU(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)     ! updated zonal velocity
REAL(KIND=JPRB) :: ZV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)     ! updated meridional velocity

REAL(KIND=JPRB) :: ZQV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! corrected (for negative values) vapour
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZQI(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! corrected (for negative values) cloud ice
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZQL(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! corrected (for negative values) cloud liquid
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZQR(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! corrected (for negative values) rain
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZQS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! corrected (for negative values) snow
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZCP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! new cp for turbulent diffusion
REAL(KIND=JPRB) :: ZQT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENHA(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENQVA(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZFCQVNG(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! correction flux increment for neg vapour
REAL(KIND=JPRB) :: ZFCQING(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! correction flux increment for neg ice
REAL(KIND=JPRB) :: ZFCQLNG(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! correction flux increment for neg liquid water
REAL(KIND=JPRB) :: ZFPLSL(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! total liquid water flux: diff+sedi+rain
REAL(KIND=JPRB) :: ZFPLSN(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! total solid water flux: diff+sedi+snow
REAL(KIND=JPRB) :: ZFCQL(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)   ! condensation flux(liquid)
REAL(KIND=JPRB) :: ZFCQI(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)   ! condensation flux(ice)
REAL(KIND=JPRB) :: ZDIFCQD(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! downdraft flux of specific humidity
REAL(KIND=JPRB) :: ZDIFCQLD(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! downdraft flux of liquid water
REAL(KIND=JPRB) :: ZDIFCQID(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! downdraft flux of  solid water
REAL(KIND=JPRB) :: ZSEDIQL(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! sedimentation flux of cloud liquid water
REAL(KIND=JPRB) :: ZSEDIQI(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! sedimentation flux of cloud ice water
REAL(KIND=JPRB) :: ZDIFCSD(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! downdraft entalphy flux
REAL(KIND=JPRB) :: ZSTRCUD(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! change in horizontal mom.
REAL(KIND=JPRB) :: ZSTRCVD(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! change in horizontal mom.
REAL(KIND=JPRB) :: ZRCVOTT(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! degree of inhomogeneity in precips.
REAL(KIND=JPRB) :: ZSIGPC(YDCPG_DIM%KLON)         ! Convective precipit mesh fraction
REAL(KIND=JPRB) :: ZSIGP(YDCPG_DIM%KLON)         ! Precipitation mesh fraction
REAL(KIND=JPRB) :: ZAUXPRC(YDCPG_DIM%KLON)        ! Precipitation auxilary
REAL(KIND=JPRB) :: ZDIFCVPPQ(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de CVPP (KFB or EDKF) sur Qv
REAL(KIND=JPRB) :: ZDIFCVPPS(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de CVPP (KFB or EDKF) sur CpT
REAL(KIND=JPRB) :: ZDIFCVTH(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de CV sur Theta air sec
REAL(KIND=JPRB) :: ZDIFCVPPU(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de CVPP (EDKF) sur U
REAL(KIND=JPRB) :: ZDIFCVPPV(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de CVPP (EDKF) sur V

REAL(KIND=JPRB) :: ZEDMFQ(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Mass flux part of EDMF flux for Qv
REAL(KIND=JPRB) :: ZEDMFS(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Mass flux part of EDMF flux for CpT
REAL(KIND=JPRB) :: ZEDMFU(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Mass flux part of EDMF flux for U
REAL(KIND=JPRB) :: ZEDMFV(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Mass flux part of EDMF flux for V
REAL(KIND=JPRB) :: ZMF_UP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Mass flux for implicit formulation of EDMF equation (LEDMFI)
REAL(KIND=JPRB) :: ZMU(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de masse (updraft) pour XIOS output
REAL(KIND=JPRB) :: ZMD(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de masse (downdraft) pour XIOS output

REAL(KIND=JPRB) :: ZCONDCVPPL(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de condensation liquide du a CVVPP (KFB)
REAL(KIND=JPRB) :: ZCONDCVPPI(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de condensation glace du a CVVPP (KFB)
REAL(KIND=JPRB) :: ZPRODTH_CVPP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de production thermique de TKE du a CVPP(KFB)
REAL(KIND=JPRB) :: ZDTRAD(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! radiation contribution to T tendency
REAL(KIND=JPRB) :: ZDQVDIFF(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! turtb.diff contribution to Qv tendency
REAL(KIND=JPRB) :: ZRKQCTEND(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Qc input for RK condensation scheme
REAL(KIND=JPRB) :: ZRKQVTEND(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Qv input for RK condensation scheme
REAL(KIND=JPRB) :: ZRKTTEND(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! T input for RK condensation scheme
REAL(KIND=JPRB) :: ZDQV, ZDQI, ZDQL, ZDQR, ZDQS, ZDQC, ZGDT, ZGDTI,&
                 & ZQV0, ZQX0, ZQX1,&
                 & ZCONVC, ZTOTC,ZDTURDIFF
REAL(KIND=JPRB) :: ZTMPAF(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! temporary array for Add_Field_3d.
REAL(KIND=JPRB) :: ZTMPAS(YDCPG_DIM%KLON)         ! temporary array for Add_Field_2d..
REAL(KIND=JPRB) :: ZTMPPRODTH(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! temporary array

!!for BAYRAD
REAL(KIND=JPRB) :: ZDE2MR(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! temporary array for conversion of density to mixing ratio
REAL(KIND=JPRB) :: ZCOEFRAIN(2) ! RTTOVSCATT coefficients to convert flux to density for rain
REAL(KIND=JPRB) :: ZCOEFSNOW(2) ! RTTOVSCATT coefficients to convert flux to density for snow
REAL(KIND=JPRB) :: ZRHORAIN ! RTTOVSCATT fixed density for rain
REAL(KIND=JPRB) :: ZRHOSNOW ! RTTOVSCATT fixed density of snow
!-----------------------------------------------------------------

! - 2D (0:KLEV) .

! ZKTROV     : COEFFICIENT D'ECHANGE VERTICAL DE T ET Q EN KG/(M*M*S).
! ZKUROV     : COEFFICIENT D'ECHANGE VERTICAL DE U ET V EN KG/(M*M*S).
! ZKNROV     : COEFFICIENT D'ECHANGE VERTICAL NEUTRE EN KG/(M*M*S).
! ZNBVNO     : FREQUENCE DE BRUNT-VAISALA DIVISEE PAR G FOIS LA DENSITE.

! ZNEBS      : NEBULOSITE STRATIFORME (SCHEMA STATISTIQUE DE NUAGES).
!            : STRATIFORM CLOUDINESS (STATISTICAL CLOUD SCHEME)
! ZQLIS      : QUANTITE D'EAU LIQUIDE STRATIFORME (SCHEMA STATISTIQUE).
!            : STRATIFORM LIQUID WATER (STATISTICAL CLOUD SCHEME)


! - 2D (1:KLEV) .

! INLAB      : INDICE D'INSTABILITE CONVECTIVE.

INTEGER(KIND=JPIM) :: INND(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZXDROV(YDCPG_DIM%KLON),ZXHROV(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZUGST(YDCPG_DIM%KLON),ZVGST(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZCDROV(YDCPG_DIM%KLON),ZCHROV(YDCPG_DIM%KLON),ZDQSTS(YDCPG_DIM%KLON),ZGWDCS(YDCPG_DIM%KLON),&
 & ZHQ(YDCPG_DIM%KLON),ZHU(YDCPG_DIM%KLON),ZHTR(YDCPG_DIM%KLON),ZCDNH(YDCPG_DIM%KLON),ZMOD(YDCPG_DIM%KLON),&
 & ZRTI(YDCPG_DIM%KLON),ZDPHI(YDCPG_DIM%KLON),ZPRS(YDCPG_DIM%KLON),ZSTAB(YDCPG_DIM%KLON),ZTAUX(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZWFC(YDCPG_DIM%KLON),ZWPMX(YDCPG_DIM%KLON),ZWLMX(YDCPG_DIM%KLON),ZWSEQ(YDCPG_DIM%KLON),&
 & ZWSMX(YDCPG_DIM%KLON),ZWWILT(YDCPG_DIM%KLON),&
 & ZC3(YDCPG_DIM%KLON),ZCG(YDCPG_DIM%KLON),ZCN(YDCPG_DIM%KLON),&
 & ZNEIJG(YDCPG_DIM%KLON),ZNEIJV(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZPCLS(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZPREN(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZFRSODS(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZCD(YDCPG_DIM%KLON)

! - 1D (DIAGNOSTIQUE) .

! ZCDROV     : PCD RENORME EN DENSITE FOIS VITESSE.
! ZCDNH      : COEFFICIENT NEUTRE D'ECHANGE EN SURFACE POUR LA CHALEUR.
! ZCDNH      : EXCHANGE COEFF. AT SURFACE LEVEL IN NEUTRAL CONDITIONS FOR HEAT.
! ZCG        : COEFFICIENT THERMIQUE DU SOL NU.
! ZCG        : THERMICAL COEFFICIENT OF BARE GROUND.
! ZCN        : COEFFICIENT THERMIQUE DE LA NEIGE.
! ZCN        : THERMICAL COEFFICIENT OF SNOW.
! ZCHROV     : PCH RENORME EN DENSITE FOIS VITESSE.
! ZC3        : COEFFICIENT UTILE POUR LE CALCUL DU DRAINAGE
! ZDQSTS     : DERIVEE DE PQSATS PAR RAPPORT A LA TEMPERATURE.
! ZGWDCS     : VARIABLE DE SURFACE POUR LE DRAG OROGRAPHIQUE (RHO*N0/G).
! ZHQ        : POIDS DE L'HUMIDITE DE L'AIR DANS L'HUMIDITE DE SURFACE.
! ZHTR       : RESISTANCE A LA TRANSPIRATION DU COUVERT VEGETAL.
! ZHTR       : FOLIAGE TRANSPIRATION RESISTANCE.
! ZHU        : POIDS DE L'HUMIDITE SATURANTE DANS L'HUMIDITE DE SURFACE.
! ZNEIJG     : FRACTION DE NEIGE RECOUVRANT LE SOL.
! ZNEIJV     : FRACTION DE NEIGE RECOUVRANT LA VEGETATION.
! ZRTI       : INVERSE DE R*T.
! ZDPHI      : EPAISSEUR EN GEOPOTENTIEL DU NIVEAU DE SURFACE.
! ZPRS       : CONSTANTE DES GAZ POUR L'AIR AU SOL.
! ZSTAB      : INDICE DE STABILITE A LA SURFACE.
! INND       : INDICE DE PRECIPITATIONS CONVECTIVES.
! ZWFC       : TENEUR EN EAU A LA CAPACITE AUX CHAMPS.
! ZWFC       : FIELD CAPACITY WATER CONTENT.
! ZWPMX      : TENEUR EN EAU MAXIMALE DU RESERVOIR PROFOND.
! ZWPMX      : MAXIMUM WATER CONTENT OF THE DEEP WATER-TANK.
! ZWLMX      : TENEUR EN EAU MAXIMALE DU RESERVOIR D'INTERCEPTION.
! ZWLMX      : MAXIMUM WATER CONTENT OF THE INTERCEPTION WATER-TANK.
! ZWSEQ      : TENEUR EN EAU A L'EQUILIBRE (EQUILIBRE ENTRE GRAVITE ET
!              CAPILLARITE) EN SURFACE.
! ZWSEQ      : SURFACE WATER CONTENT FOR THE BALANCE BETWEEN GRAVITY
!              AND CAPILLARITY
! ZWSMX      : TENEUR EN EAU MAXIMALE DU RESERVOIR SUPERFICIEL.
! ZWSMX      : MAXIMUM WATER CONTENT FOR THE SUPERFICIAL WATER-TANK.
! ZWWILT     : TENEUR EN EAU AU POINT DE FLETRISSEMENT.
! ZWWILT     : WATER CONTENT AT THE WILTING POINT.
! ZSSO_STDEV : OROGRAPHY STANDARD DEVIATION
! ZTWSNOW    : SNOW COVER FROM SURFEX
! ZTOWNS     : FRACTION OF TOWN FROM SURFEX
REAL(KIND=JPRB) :: ZDAER(YDCPG_DIM%KFLEVG), ZHUC(YDCPG_DIM%KFLEVG), ZBLH(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZQO3(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZAER(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,6)
REAL(KIND=JPRB) :: ZAERD(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZAERINDS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZQCO2(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZROZ(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZDPHIV(YDCPG_DIM%KLON),ZDPHIT(YDCPG_DIM%KLON)

REAL(KIND=JPRB) :: ZMAN(0:YDCPG_DIM%KFLEVG), ZMAK(0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZSSO_STDEV(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZTWSNOW(YDCPG_DIM%KLON),ZTOWNS(YDCPG_DIM%KLON)
! - (PROFILS DIAGNOSTIQUES)
! ZDAER     : EPAISSEUR OPTIQUE DES AEROSOLS STANDARDS DANS LA COUCHE
! ZDAER     : OPTICAL DEPTH OF STANDARD AEROSOLS IN THE LAYER.
! ZHUC      : HUMIDITE CRITIQUE PAR NIVEAU POUR LE CALCUL DE PNEB
! ZHUC      : CRITICAL MOISTURE FOR EACH LEVEL FOR PNEB CALCULATION.
!   ZQO3    : RAPPORT DE MELANGE MASSIQUE D'OZONE
!        (0):     "         "    MOYEN AU-DESSUS DU MODELE
!   ZQO3    : OZONE MIXING RATIO (MASS).
!        (0):AVERAGED-ABOVE      "         " .
!   ZQCO2   : RAPPORT MASSIQUE LOCAL DU CO2.
!   ZQCO2   : CO2 MIXING RATIO (MASS).

! IJN        : DIMENSION TABLEAUX ETENDUS POUR CYCLE DIURNE RAYONNEMENT
!              IJN AU PLUSL A KLON

!* INPUT ARGUMENTS FOR ACRADIN ( RAYT ECMWF POUR CLIMAT )

!            1-D ARRAYS
!            ----------

REAL(KIND=JPRB) :: ZTRSOD(YDCPG_DIM%KLON)

!            2-D ARRAYS
!            ----------

!* OUTPUT ARGUMENTS FOR THE ECMWF PHYSICS

!            0.2  LOCAL ARRAYS FOR ECMWF PHYSICS PACKAGE
!                 --------------------------------------

REAL(KIND=JPRB) :: ZCEMTR(YDCPG_DIM%KLON,0:1) , ZCTRSO(YDCPG_DIM%KLON,0:1)
REAL(KIND=JPRB) :: ZALBD(YDCPG_DIM%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW), ZALBP(YDCPG_DIM%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW),&
                 & ZALB(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZSFSWDIR (YDCPG_DIM%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW), ZSFSWDIF (YDCPG_DIM%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB) :: ZTRSODIR (YDCPG_DIM%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW), ZTRSODIF (YDCPG_DIM%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB) :: ZFSDNN(YDCPG_DIM%KLON),ZFSDNV(YDCPG_DIM%KLON)

!            1-D ARRAYS
!            ----------

REAL(KIND=JPRB) :: ZSUDU(YDCPG_DIM%KLON) , ZDSRP(YDCPG_DIM%KLON) , ZSDUR(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZTHETAVS(YDCPG_DIM%KLON), ZTHETAS(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZAESEA(YDCPG_DIM%KLON), ZAELAN(YDCPG_DIM%KLON), ZAESOO(YDCPG_DIM%KLON), ZAEDES(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZAESUL(YDCPG_DIM%KLON), ZAEVOL(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZMERL(YDCPG_DIM%KLON)

!            2-D ARRAYS
!            ----------

REAL(KIND=JPRB) :: ZTENT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZGEOSLC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTHETAV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

!            LOCAL ARRAYS FOR TKE
! ZCOEFN : COEFFICIENT STATISTIQUE POUR LES FLUX D'EAUX CONDENSEES.

REAL(KIND=JPRB) :: ZCOEFN(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

!            LOCAL ARRAYS FOR ACVPPKF
REAL(KIND=JPRB) :: ZQLI_CVPP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZNEB_CVPP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
!            LOCAL ARRAYS FOR EDKF
REAL(KIND=JPRB) :: ZIMPL

!           2-D ARRAY FOR SIMPL.RADIATION SCHEME

REAL(KIND=JPRB) :: ZZNEB(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

INTEGER(KIND=JPIM) :: IJN, IMLTYPE, JCHA, JLEV, JLON, JSG, JGFL, JDIAG
INTEGER(KIND=JPIM) :: ILONMNH, IKRR, ISPLITR ! useful size of klon arrays for mesonh physics
LOGICAL :: LLCLS, LLHMT, LLMAF, LLREDPR

REAL(KIND=JPRB) :: ZAEN, ZAEO, ZALBV, ZCARDI, ZEPS, ZEPS0, ZEPSNEB, ZEPSO3
REAL(KIND=JPRB) :: ZEPSA, ZALBPMER

!            2-D ARRAYS

! ZNEBC      : NEBULOSITE  CONVECTIVE A L'ECHELLE DE LA MAILLE.
! ZQLIC      : EAU LIQUIDE CONVECTIVE A L'ECHELLE DE LA MAILLE.
! ZQCL       : CONDENSAT STRATIFORME LIQUIDE
! ZQCI       : CONDENSAT STRATIFORME SOLIDE
! ZFHEVPPC   : FLUX DE CHALEUR DU A L'EVAPORATION DES PREC. CONVECTIVES.
! ZFHMLTSC   : FLUX DE CHALEUR DU A LA FONTE/GEL DES PREC. CONVECTIVES.
! ZFPEVPPC   : EVAPORATION DES PREC. CONVECTIVES.
! ICIS       : INDICE DE NIVEAU D'INSTABILITE SECHE.

REAL(KIND=JPRB) :: ZNEBC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZQLIC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZQCL(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZQCI(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZFHEVPPC(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZFHMLTSC(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)&
 & ,ZFPEVPPC(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
INTEGER(KIND=JPIM) :: ICIS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)


!           SURFEX local VARIABLES
!           ----------------------------

! Implicit coupling coefficients
INTEGER(KIND=JPIM) :: IRR
REAL(KIND=JPRB) :: ZDTMSE,ZRHGMT,ZSTATI
REAL(KIND=JPRB) :: ZCFAQ(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZCFAS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),&
 & ZCFATH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZCFAU(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZCFBQ(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),&
 & ZCFBS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZCFBTH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),&
 & ZCFBU(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZCFBV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),&
 & ZCFBU_G(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZCFBV_G(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),&
 & ZCFBS_G(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZCFBQ_G(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

REAL(KIND=JPRB) :: ZDSE(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZFEV(YDCPG_DIM%KLON),ZFMDU(YDCPG_DIM%KLON),ZFMDV(YDCPG_DIM%KLON),ZFEVS(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZSRAIN(YDCPG_DIM%KLON), ZSSNOW(YDCPG_DIM%KLON), ZSGROUPEL(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZSFCO2(YDCPG_DIM%KLON), ZRHODREFM(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZDEPTH_HEIGHT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZZS(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZTSN(YDCPG_DIM%KLON),ZTN(YDCPG_DIM%KLON)
REAL(KIND=JPRB)   :: ZBUDTH (YDCPG_DIM%KLON)
REAL(KIND=JPRB)   :: ZBUDSO (YDCPG_DIM%KLON)
REAL(KIND=JPRB)   :: ZFCLL  (YDCPG_DIM%KLON)
! FOR Hv
REAL(KIND=JPRB)   :: ZHV2(YDCPG_DIM%KLON)
! FOR DUST
REAL(KIND=JPRB), DIMENSION (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ::  ZQDM
REAL(KIND=JPRB) :: ZCFASV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,1:YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) ! SCOND MEMBRE POUR LES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZCFBSV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,1:YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) ! SCOND MEMBRE POUR LES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZSMOOTRAC(1:YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) ! SCOND MEMBRE POUR LES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZINVDT, ZINVG, ZRSCP, ZINVATM

REAL(KIND=JPRB),  DIMENSION(YDCPG_DIM%KLON,1,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT):: ZZI_SVM
REAL(KIND=JPRB),  DIMENSION (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG):: ZZI_PEZDIAG
REAL(KIND=JPRB), DIMENSION (:,:,:), ALLOCATABLE  :: ZSVM, ZPSV
REAL(KIND=JPRB), DIMENSION (:,:),   ALLOCATABLE  :: ZSFSV  ! passifs scalaires surf flux
! TRAITEMENT DES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZTM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB), DIMENSION (YDCPG_DIM%KLON,1,YDCPG_DIM%KFLEVG) :: ZZZ,ZDZZ,ZZI_PABSM, ZZI_THM,&
                & ZZI_EXNREFM, ZZI_RHODREFM,ZEVAP,ZZDEP,ZZI_RHO
REAL(KIND=JPRB) :: ZZI_APHI(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB), DIMENSION (YDCPG_DIM%KLON,1,YDCPG_DIM%KFLEVG,6) :: ZZI_RM
!            3-D ARRAYS
REAL(KIND=JPRB), DIMENSION(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_PHY_RAD%YRERAD%NSW):: ZPIZA_DST !Single scattering
                                             ! albedo of dust (points,lev,wvl)
REAL(KIND=JPRB), DIMENSION(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_PHY_RAD%YRERAD%NSW):: ZCGA_DST  !Assymetry factor
                                             ! for dust (points,lev,wvl)
REAL(KIND=JPRB), DIMENSION(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_PHY_RAD%YRERAD%NSW):: ZTAUREL_DST !tau/tau_{550}
                                             !dust (points,lev,wvl)

!         ACFLUSO (ECUME) local variable
!-------------------------------------------
REAL(KIND=JPRB) :: ZCE(YDCPG_DIM%KLON), ZCEROV(YDCPG_DIM%KLON), ZCRTI(YDCPG_DIM%KLON)


!        New ACDIFV1 local variable
!--------------------------------------------
REAL(KIND=JPRB)   :: ZXURO(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)   :: ZXQRO(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)   :: ZXTRO(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)

!        New ACNORGWD local variables
!--------------------------------------------
REAL(KIND=JPRB) :: ZD_U(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZD_V(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: Z_PP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) 
REAL(KIND=JPRB) :: Z_UU(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) 
REAL(KIND=JPRB) :: Z_VV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: Z_TT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: Z_VO(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)  
REAL(KIND=JPRB) :: ZFLX_LOTT_GWU(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG), ZFLX_LOTT_GWV(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZPRECGWD(YDCPG_DIM%KLON)

!    TKE+ for ACCLDIA
REAL(KIND=JPRB)   :: ZTKE1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)   :: ZTPRDY(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

!   For ACVISIH
REAL(KIND=JPRB)   :: ZQGM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)


!        New ARP_GROUND_PARAM local variable
!------------------------------------------------

REAL(KIND=JPRB)   :: ZALPHA1(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)   :: ZCOEFA (YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)   :: ZLVT   (YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)   :: ZQICE  (YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)   :: ZDIFWQ (YDCPG_DIM%KLON)
REAL(KIND=JPRB)   :: ZDIFWS (YDCPG_DIM%KLON)
REAL(KIND=JPRB)   :: ZSC_FEVI (YDCPG_DIM%KLON),ZSC_FEVN(YDCPG_DIM%KLON),ZSC_FCLL(YDCPG_DIM%KLON),ZSC_FCLN(YDCPG_DIM%KLON)

!           TRAJECTORY (For diffusion !) local VARIABLES
!           ----------------------------
REAL(KIND=JPRB) :: ZCDROV_SAVE(YDCPG_DIM%KLON),ZCHROV_SAVE(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZKTROV_SAVE(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZKUROV_SAVE(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTRAJGWD(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) !Traj buffer saved for TL/AD (YDMODEL%YRML_PHY_MF%YRSIMPHL%LGWDSPNL)

REAL(KIND=JPRB)    :: ZRVMD,ZDELTA
LOGICAL            :: LLAERO, LLAROME, LLCALLRAD
REAL(KIND=JPRB)    :: ZAIPCMT(YDCPG_DIM%KLON) ! Activity Index of PCMT: 1. if PCMT is active, 0. else case.
REAL(KIND=JPRB)    :: ZALF_CAPE(YDCPG_DIM%KLON)
REAL(KIND=JPRB)    :: ZALF_CVGQ(YDCPG_DIM%KLON)
REAL(KIND=JPRB)    :: ZQIC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZQLC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZQRC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZQSC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZQVI
REAL(KIND=JPRB)    :: ZQLI_CVP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)    :: ZQC_DET_PCMT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)    :: ZFPLS(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZFPLC(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZFPL(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)    :: ZCSGC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)    :: ZTZER
CHARACTER(LEN=200) :: CLERR

! Tracers: prognostique aerosols, passive scalars...
INTEGER(KIND=JPIM) :: INBTRA
INTEGER(KIND=JPIM) :: INBTRA_DEP
REAL(KIND=JPRB), ALLOCATABLE :: ZSTRCTRA(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZTRA(:,:,:)

!        New chemistry local variables
!-------------------------------------
INTEGER(KIND=JPIM)               :: IFLDX, IFLDX2, ILEVX 
INTEGER(KIND=JPIM), ALLOCATABLE  :: INDCHEM(:), IGPLAT(:)
REAL(KIND=JPRB), ALLOCATABLE     :: ZSD_XA(:,:,:), ZSD_X2(:,:), ZTENGFL(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE     :: ZCFLX(:,:), ZCFLXO(:,:), ZCHEMDV(:,:), ZAEROP(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE     :: ZTENC(:,:,:), ZDELP(:,:), ZWND(:), ZDUMMY1(:,:), ZGELAT(:)
REAL(KIND=JPRB), ALLOCATABLE     :: ZNEEFLX(:), ZCHEM2AER(:,:,:)
CHARACTER (LEN = 20) , POINTER ::   CGMIXLEN
REAL(KIND=JPRB) :: ZDCAPE(YDCPG_DIM%KLON) ! Descending CAPE for gusts.

! ACRANEB/ACRANEB2 local variables
! --------------------------------
REAL(KIND=JPRB) :: ZLAMB           ! proportion of Lambertian scattering
REAL(KIND=JPRB) :: ZALBDIR  (YDCPG_DIM%KLON) ! direct (parallel) surface albedo
REAL(KIND=JPRB) :: ZCLCT_RAD(YDCPG_DIM%KLON) ! total cloud cover for radiation
REAL(KIND=JPRB) :: ZDECRD   (YDCPG_DIM%KLON) ! decorrelation depth for cloud overlaps
REAL(KIND=JPRB) :: ZDECRD_MF(YDCPG_DIM%KLON) ! decorrelation depth for cloud overlaps
                                   ! in microphysics

REAL(KIND=JPRB) :: ZQG(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZDQG

! IFS deep convection scheme local variables
! --------------------------------
INTEGER(KIND=JPIM) :: ITOPC(YDCPG_DIM%KLON),IBASC(YDCPG_DIM%KLON),ITYPE(YDCPG_DIM%KLON),ISPPN2D
INTEGER(KIND=JPIM) :: ICBOT(YDCPG_DIM%KLON),ICTOP(YDCPG_DIM%KLON),IBOTSC(YDCPG_DIM%KLON)
LOGICAL :: LLDSLPHY,LLPTQ,LLLAND(YDCPG_DIM%KLON),LLCUM(YDCPG_DIM%KLON),LLSC(YDCPG_DIM%KLON),LLSHCV(YDCPG_DIM%KLON),LLLINOX(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZLIGH_CTG(YDCPG_DIM%KLON),ZCTOPH(YDCPG_DIM%KLON),ZPRECMX(YDCPG_DIM%KLON),ZICE(YDCPG_DIM%KLON),ZCDEPTH(YDCPG_DIM%KLON),ZWMFU(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZVERVEL(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZGEOM1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZGEOMH(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENQ(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZTENU(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZTENV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZLCRIT_AER(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZSNDE(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,2)
REAL(KIND=JPRB) :: ZCUCONVCA(YDCPG_DIM%KLON),ZGAW(YDCPG_DIM%KLON),ZVDIFTS,ZDXTDK(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZLU(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZLUDE(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZMFU(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZLISUM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZMFD(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZWMEAN(YDCPG_DIM%KLON),ZACPR(YDCPG_DIM%KLON),ZDIFF(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZMFUDE_RATE(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZMFDDE_RATE(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZFHPCL(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZFHPCN(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZFCQLF(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZFCQLI(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZFRSO(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZFRTH(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZGP2DSPPA(YDCPG_DIM%KLON,1),ZLUDELI(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,4),ZLRAIN(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZRSUD(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,2)
REAL(KIND=JPRB), ALLOCATABLE :: ZCEN(:,:,:),ZSCAV(:)

! Precipitation type diagnostics
!--------------------------------
REAL(KIND=JPRB)   :: ZFPLSG(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)

REAL(KIND=JPRB) :: ZENTCH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)


REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "aplpar.intfb.h"
#include "aplpar2intflex.intfb.h"
#include "cpchet.intfb.h"
#include "cpmvvps.intfb.h"
#include "cpnudg.intfb.h"
#include "cpozo.intfb.h"
#include "cpphinp.intfb.h"
#include "cpqsol.intfb.h"
#include "cptend_new.intfb.h"
#include "cptend_flex.intfb.h"
#include "cptends.intfb.h"
#include "cputqy.intfb.h"
#include "cpwts.intfb.h"
#include "cp_ptrslb1.intfb.h"
#include "mf_phys_init.intfb.h"
#include "aplpar_init.intfb.h"
#include "profilechet.intfb.h"
#include "writephysio.intfb.h"
#include "wrphtrajm.intfb.h"
#include "wrradcoef.intfb.h"
#include "acajucv.intfb.h"
#include "mf_phys_nhqe_part1.intfb.h"
#include "mf_phys_nhqe_part2.intfb.h"
#include "mf_phys_save_phsurf_part1.intfb.h"
#include "mf_phys_save_phsurf_part2.intfb.h"
#include "mf_phys_corwat.intfb.h"
#include "mf_phys_transfer.intfb.h"
#include "mf_phys_fpl_part1.intfb.h"
#include "mf_phys_fpl_part2.intfb.h"
#include "mf_phys_mocon.intfb.h"
#include "mf_phys_precips.intfb.h"
#include "mf_phys_bayrad.intfb.h"
#include "mf_phys_cvv.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MF_PHYS_PAR',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDVAB=>YDGEOMETRY%YRVAB, YDCSGEOM=>   &
& YDGEOMETRY%YRCSGEOM(YDCPG_DIM%KBL), YDGSGEOM=>YDGEOMETRY%YRGSGEOM(YDCPG_DIM%KBL), YDOROG=>YDGEOMETRY%YROROG(YDCPG_DIM%KBL),  &
& YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, YDPTRSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1, YDPTRSLB2=>YDMODEL%      &
& YRML_DYN%YRPTRSLB2, YDTOPH=>YDMODEL%YRML_PHY_MF%YRTOPH, YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL,   &
& YDRIP=>YDMODEL%YRML_GCONF%YRRIP, YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH, YDRCOEF=>YDMODEL%YRML_PHY_RAD% &
& YRRCOEF, YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY, YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY, YDLDDH=>YDMODEL%&
& YRML_DIAG%YRLDDH, YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, YGFL=>YDMODEL%YRML_GCONF%YGFL, YDEPHY=>     &
& YDMODEL%YRML_PHY_EC%YREPHY, YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR, YDPRECIPS=>YDMODEL%YRML_PHY_MF% &
& YRPHY%YRDPRECIPS, YDSTOPH=>YDMODEL%YRML_PHY_STOCH%YRSTOPH)

ASSOCIATE(MVTS=>YDPARAR%MVTS, CMF_UPDRAFT=>YDPARAR%CMF_UPDRAFT, TSPHY=>YDPHY2%TSPHY, NPROMA=>YDDIM%   &
& NPROMA, NTSSG=>YDDPHY%NTSSG, NVCLIS=>YDDPHY%NVCLIS, LMDUST=>YDARPHY%LMDUST, &
& LMSE=>YDARPHY%LMSE, LMFSHAL=>YDARPHY%LMFSHAL, YI=>YGFL%YI, YH=>YGFL%YH, YEZDIAG=>YGFL%YEZDIAG, YL   &
& =>YGFL%YL, YEXT=>YGFL%YEXT, YG=>YGFL%YG, YQ=>YGFL%YQ, YR=>YGFL%YR, YSCONV=>YGFL%YSCONV, YS=>YGFL%   &
& YS, YEFB3=>YGFL%YEFB3, YEFB2=>YGFL%YEFB2, YEFB1=>YGFL%YEFB1, YO3=>YGFL%YO3, YNOGW=>YGFL%YNOGW,      &
& LCHEM_ARPCLIM=>YDMODEL%YRML_CHEM%YRCHEM%LCHEM_ARPCLIM, YCHEM=>YGFL%YCHEM, NGFL_EXT=>YGFL%NGFL_EXT,  &
& YTKE=>YGFL%YTKE, YCOMP=>YGFL%YCOMP, YLCONV=>YGFL%YLCONV, YRCONV=>YGFL%YRCONV, YICONV=>YGFL%YICONV,  &
& YSP_SBD=>YDSURF%YSP_SBD, YSD_VVD=>YDSURF%YSD_VVD, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA,       &
& NFLSUL=>YDDIMV%NFLSUL, LTRAJPS=>YDSIMPHL%LTRAJPS, &
& LNEBR=>YDPHY%LNEBR, LNEBN=>YDPHY%LNEBN, LSTRAPRO=>YDPHY%LSTRAPRO, LPTKE=>     &
& YDPHY%LPTKE, NDPSFI=>YDPHY%NDPSFI, LOZONE=>YDPHY%LOZONE, L3MT=>YDPHY%L3MT, LGPCMT=>YDPHY%LGPCMT,    &
& LAJUCV=>YDPHY%LAJUCV, LCVPGY=>YDPHY%LCVPGY, LRRGUST=>YDPHY%LRRGUST, LEDR=>YDPHY%LEDR, NTAJUC=>      &
& YDTOPH%NTAJUC, NTPLUI=>YDTOPH%NTPLUI, NUMFLDS=>YGFL%NUMFLDS, LAGPHY=>YDEPHY%LAGPHY, &
& LDPRECIPS=>YDPHY%LDPRECIPS, LDPRECIPS2=>YDPHY%LDPRECIPS2, NDTPRECCUR=>YDPRECIPS%           &
& NDTPRECCUR, NDTPRECCUR2=>YDPRECIPS%NDTPRECCUR2, LSDDH=>YDLDDH%LSDDH, HDSF=>YDMDDH%HDSF, MSLB1SP9=>  &
& YDPTRSLB1%MSLB1SP9, MSLB2KAPPAH=>YDPTRSLB2%MSLB2KAPPAH, MSLB2KAPPAM=>YDPTRSLB2%MSLB2KAPPAM, LRCOEF  &
& =>YDRCOEF%LRCOEF, LTLADDIA=>YDRCOEF%LTLADDIA, NG3SR=>YDRCOEF%NG3SR, NLIMA=>YGFL%NLIMA, YLIMA=>YGFL  &
& %YLIMA)

CLNEW = ''
CALL GETENV ('APLPAR_NEW', CLNEW)
LLNEW = (CLNEW /= '') .AND. (CLNEW /= '0')

CALL SC2PRG(YEFB1%MP1 ,ZTENDGFL ,ZPTENDEFB11)
CALL SC2PRG(YEFB2%MP1 ,ZTENDGFL ,ZPTENDEFB21)
CALL SC2PRG(YEFB3%MP1 ,ZTENDGFL ,ZPTENDEFB31)
CALL SC2PRG(YG%MP1    ,ZTENDGFL ,ZPTENDG1)
CALL SC2PRG(YICONV%MP1,ZTENDGFL ,ZPTENDICONV1)
CALL SC2PRG(YI%MP1    ,ZTENDGFL ,ZPTENDI1)
CALL SC2PRG(YLCONV%MP1,ZTENDGFL ,ZPTENDLCONV1)
CALL SC2PRG(YL%MP1    ,ZTENDGFL ,ZPTENDL1)
CALL SC2PRG(YQ%MP1    ,ZTENDGFL ,ZPTENDQ1)
CALL SC2PRG(YRCONV%MP1,ZTENDGFL ,ZPTENDRCONV1)
CALL SC2PRG(YR%MP1    ,ZTENDGFL ,ZPTENDR1)
CALL SC2PRG(YSCONV%MP1,ZTENDGFL ,ZPTENDSCONV1)
CALL SC2PRG(YS%MP1    ,ZTENDGFL ,ZPTENDS1)
CALL SC2PRG(YTKE%MP1  ,ZTENDGFL ,ZPTENDTKE1)

CALL SC2PRG(1,YEZDIAG(:)%MP,PGFL,ZP1EZDIAG)

CALL YLMF_PHYS_STATE%INIT (LTWOTL, YDCPG_DYN0, YDCPG_DYN9, YDCPG_PHY0, YDCPG_PHY9, YDVARS, YDMF_PHYS_SURF, PGFL=PGFL, YDMODEL=YDMODEL)

!     ------------------------------------------------------------------

!        0.    constructor for procset
IF (LINTFLEX) YLPROCSET=NEWINTPROCSET()

!        1.    Preliminary calculations necessary
!              for all types of physics.
!              ------------------------------------


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
   ZGP2DSPP(:,JROF) = YSPP%GP_ARP(JROF)%GP2D(:,1,YDCPG_DIM%KBL)
 ENDDO
ENDIF

! Complete physics is called.
LLDIAB=.NOT.LAGPHY


! In the NHQE model, MF_PHYS_PAR enters with Tt and grad(Tt), where Tt = T * exp(-(R/cp) log(pre/prehyd)).
! But calculations of MF_PHYS_PAR must use T and grad(T).
! So we do a conversion Tt -> T.
IF (LNHQE) THEN
  CALL MF_PHYS_NHQE_PART1 (YDGEOMETRY, YDCPG_DIM, YDMF_PHYS_TMP, YDVARS, YDMODEL, PGFL)
ENDIF

IF (LLDIAB) THEN
  CALL CPPHINP(YDGEOMETRY, YDMODEL, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDGSGEOM%GEMU, YDGSGEOM%GELAM, YDVARS%U%T0, YDVARS%V% &
  & T0, YDVARS%Q%T0, YDVARS%Q%DL, YDVARS%Q%DM, YDVARS%CVGQ%DL, YDVARS%CVGQ%DM, YDCPG_PHY0%XYB%RDELP, &
  & YDCPG_DYN0%CTY%EVEL, YDVARS%CVGQ%T0, YDMF_PHYS_TMP%RDG%MU0, YDMF_PHYS_TMP%RDG%MU0LU, YDMF_PHYS_TMP%RDG%MU0M, YDMF_PHYS_TMP%RDG%MU0N, YDMF_PHYS_TMP%RDG%CVGQ)
  YDMF_PHYS_TMP%RDG%LCVQ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YDMF_PHYS_TMP%RDG%CVGQ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
ENDIF

IF (LLDIAB) THEN

  DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    YDAPLPAR_TMP%FLU%QSATS(JROF)=0.0_JPRB
  ENDDO

  CALL MF_PHYS_FPL_PART1 (YDCPG_DIM, YDAPLPAR_TMP, YDVARS, YDMODEL)

ENDIF

! * In some cases, some pseudo-historic surface buffers (like z0) should
!   not be modified between the entrance and the output of MF_PHYS_PAR
!   (this is the case for example if LDCONFX=T).
!   For the time being, we must save:
!   - HV (group VV) : resistance to evapotranspiration
!   - Z0F (group VD): gravity * surface roughness length
!   - Z0H (group VV): gravity * roughness length for heat
!   - PBLH (group VH): PBL height
!   - SPSH (group VH):
!   - QSH (group VH):

LL_SAVE_PHSURF = .FALSE.

IF (LLDIAB) THEN
  LL_SAVE_PHSURF=LLDIAB.AND.LDCONFX
  IF (LL_SAVE_PHSURF) THEN
    CALL MF_PHYS_SAVE_PHSURF_PART1 (YDCPG_DIM, YDMF_PHYS_TMP, YDMF_PHYS_SURF, YDVARS, YDSURF, YDMODEL)
  ENDIF
ENDIF


CALL MF_PHYS_INIT ( YGFL, YDARPHY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, NTSSG, YSP_SBD%NLEVS,&
  & YDMF_PHYS%OUT%DIFCQ , YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS , YDMF_PHYS%OUT%DIFTQ , YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL,&
  & YDMF_PHYS%OUT%DIFTS , YDMF_PHYS%OUT%FCCQL , YDMF_PHYS%OUT%FCCQN , YDMF_PHYS%OUT%FCSQL , YDMF_PHYS%OUT%FCSQN , YDMF_PHYS%OUT%FCQNG , YDMF_PHYS%OUT%FCQNNG,&
  & YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG,YDMF_PHYS%OUT%FCQGNG,&
  & YDMF_PHYS%OUT%FPLCL , YDMF_PHYS%OUT%FPLCN , YDMF_PHYS%OUT%FPLCG , YDMF_PHYS%OUT%FPLCH, YDMF_PHYS%OUT%FPLSL , YDMF_PHYS%OUT%FPLSN ,YDMF_PHYS%OUT%FPLSG ,&
  & YDMF_PHYS%OUT%FPLSH, YDMF_PHYS%OUT%FRSO  , YDMF_PHYS%OUT%FRSOC ,&
  & YDMF_PHYS%OUT%FRTH  , YDMF_PHYS%OUT%FRTHC , YDMF_PHYS%OUT%STRCU , YDMF_PHYS%OUT%STRCV , YDMF_PHYS%OUT%STRDU , YDMF_PHYS%OUT%STRDV , YDMF_PHYS%OUT%STRTU ,&
  & YDMF_PHYS%OUT%STRTV , YDMF_PHYS%OUT%STRMU , YDMF_PHYS%OUT%STRMV , YDMF_PHYS%OUT%FRMH  , &
  & YDMF_PHYS%OUT%DIFCQLC,YDMF_PHYS%OUT%DIFCQIC,YDMF_PHYS%OUT%FIMCC,&
  & YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC,YDMF_PHYS%OUT%FCNEGQIC,YDMF_PHYS%OUT%FCNEGQRC,YDMF_PHYS%OUT%FCNEGQSC,&
  & YDMF_PHYS%OUT%FCHOZ , &
  & YDCPG_MISC%NEB   , YDCPG_MISC%QICE  , YDCPG_MISC%QLI   , &
  & YDCPG_MISC%RH    , YDMF_PHYS%OUT%ALB , &
  & YDMF_PHYS%OUT%CT    , YDMF_PHYS%OUT%FCHSP , YDMF_PHYS%OUT%FCLL  , YDMF_PHYS%OUT%FCLN  ,&
  & YDMF_PHYS%OUT%FCS   , YDMF_PHYS%OUT%FEVL  , YDMF_PHYS%OUT%FEVN  , YDMF_PHYS%OUT%FEVV  , YDMF_PHYS%OUT%FLASH , YDMF_PHYS%OUT%FTR   , YDMF_PHYS%OUT%FLWSP ,&
  & YDMF_PHYS%OUT%FONTE , YDMF_PHYS%OUT%FGEL  , YDMF_PHYS%OUT%FGELS ,&
  & YDMF_PHYS%OUT%FRSGNI, YDMF_PHYS%OUT%FRSDNI, YDMF_PHYS%OUT%FRSODS, YDMF_PHYS%OUT%FRSOPS, YDMF_PHYS%OUT%FRSOPT, YDMF_PHYS%OUT%FRSOLU, YDMF_PHYS%OUT%FRTHDS,&
  & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN,YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN,YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPEVPSG,YDMF_PHYS%OUT%FPEVPCL,&
  & YDMF_PHYS%OUT%FPEVPCN,YDMF_PHYS%OUT%FPEVPCG, &
  & YDMF_PHYS%OUT%GZ0   , YDMF_PHYS%OUT%GZ0H  , YDCPG_MISC%QS    , YDMF_PHYS%OUT%RUISL ,&
  & YDMF_PHYS%OUT%RUISP , YDMF_PHYS%OUT%RUISS , YDMF_PHYS%OUT%UCLS  , YDMF_PHYS%OUT%VCLS  , YDMF_PHYS%OUT%NUCLS , YDMF_PHYS%OUT%NVCLS , YDMF_PHYS%OUT%TCLS  , YDMF_PHYS%OUT%MRT,&
  & YDMF_PHYS%OUT%QCLS  , YDMF_PHYS%OUT%RHCLS , YDCPG_MISC%CLCT  , YDMF_PHYS%OUT%CLCH  , YDMF_PHYS%OUT%CLCM  , YDMF_PHYS%OUT%CLCL  , YDMF_PHYS%OUT%CLCC  ,&
  & YDMF_PHYS%OUT%CAPE  , YDMF_PHYS%OUT%CTOP  , YDMF_PHYS%OUT%CLPH  , YDMF_PHYS%OUT%VEIN  , YDMF_PHYS%OUT%UGST  , YDMF_PHYS%OUT%VGST  ,&
  & YDMF_PHYS%OUT%DIAGH , YDMF_PHYS%OUT%EDR, YDMF_PHYS%OUT%VISICLD, YDMF_PHYS%OUT%VISIHYD, YDMF_PHYS%OUT%MXCLWC, YDMF_PHYS%OUT%MOCON)

IF (LMDUST) THEN
  ZDIFEXT (:,:,:) = 0.0_JPRB
ENDIF

CALL APLPAR_INIT ( YGFL, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, NTSSG, YSP_SBD%NLEVS,&
  & YDMF_PHYS_SURF%GSD_VF%PVEG ,&
  & YDAPLPAR_TMP%MSC%FRMQ  ,&
  & YDAPLPAR_TMP%DSA%CPS   , YDAPLPAR_TMP%DSA%LHS   ,&
  & YDAPLPAR_TMP%DSA%RS    , YDAPLPAR_TMP%MSC%LH    , YDAPLPAR_TMP%MSC%LSCPE , YDAPLPAR_TMP%FLU%QSAT  ,&
  & YDAPLPAR_TMP%MSC%QW    , YDAPLPAR_TMP%MSC%TW    , YDAPLPAR_TMP%FLU%CD    , YDAPLPAR_TMP%FLU%CDN   , YDAPLPAR_TMP%FLU%CH    ,&
  & YDAPLPAR_TMP%DSA%C1    , YDAPLPAR_TMP%DSA%C2    , YDAPLPAR_TMP%FLU%EMIS  , &
  & YDAPLPAR_TMP%FLU%FEVI  , &
  & YDAPLPAR_TMP%PFL%FTKE  , YDAPLPAR_TMP%PFL%FTKEI , YDAPLPAR_TMP%PFL%FEFB1 , YDAPLPAR_TMP%PFL%FEFB2   , YDAPLPAR_TMP%PFL%FEFB3,&
  & YDAPLPAR_TMP%FLU%NEIJ  , YDAPLPAR_TMP%FLU%VEG   , YDAPLPAR_TMP%FLU%QSATS , &
  & YDAPLPAR_TMP%MOC%CLPH)

!*       2.    Complete physics.
!              -----------------

!        2.2  Complete physics.
!             -----------------

IF (LLDIAB) THEN

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

  DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    YDMF_PHYS_SURF%GSD_VF%PLSM(JROF)=REAL(NINT(YDMF_PHYS_SURF%GSD_VF%PLSM(JROF)),JPRB)
  ENDDO
  
    IF (LTWOTL) THEN
      IF (LAJUCV) THEN
        DO JLEV=1,YDCPG_DIM%KFLEVG
          DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            YDAPLPAR_TMP%ADJ%TAUX(JROF,JLEV)=YDVARS%T%T0(JROF,JLEV)
          ENDDO
        ENDDO
        CALL ACAJUCV(YDMODEL%YRML_PHY_MF%YRPHY0,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,NTPLUI,YDCPG_DIM%KFLEVG,NTAJUC,&
         & YDCPG_PHY0%PREHYD,YDCPG_PHY0%XYB%ALPH,YDCPG_PHY0%XYB%DELP,&
         & YDCPG_PHY0%XYB%LNPR,YDVARS%T%T0)
        DO JLEV=1,YDCPG_DIM%KFLEVG
          DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            YDAPLPAR_TMP%ADJ%DTAJU(JROF,JLEV)=YDVARS%T%T0(JROF,JLEV)-YDAPLPAR_TMP%ADJ%TAUX(JROF,JLEV)
          ENDDO
        ENDDO
      ENDIF
   ELSE
      ! IF (LAJUCV) THEN
      !   missing code under LAJUCV for leap-frog schemes.
      ! ENDIF
   ENDIF


   IF (.NOT. LLNEW) THEN

   CALL APLPAR(YDGEOMETRY, YDCPG_DIM, YLMF_PHYS_STATE, YDCPG_DYN0, YDCPG_MISC, YDMF_PHYS, YDMF_PHYS_TMP, &
   & YDAPLPAR_TMP, YDMF_PHYS_SURF, YDVARS, YDSURF, YDXFU, YDCFU, YDMODEL, &
   & NTSSG, PDTPHY, &
   & LLXFUMSE, PGFL, PKOZO, PGPAR, &
   & ZAC, ZAC_HC, ZDIFEXT, ZP1EZDIAG, ZTENDPTKE, ZTENDEXT_DEP,       &
   & PTRAJ_PHYS, YDDDH, PFTCNS, ZGP2DSPP)

   ELSE

   ENDIF

    IF (LTWOTL) THEN
      IF (LSLAG) CALL CP_PTRSLB1(YDMODEL%YRML_DYN%YRDYN,YDPTRSLB1,ISLB1U9,ISLB1V9,ISLB1T9,ISLB1VD9,ISLB1GFL9)
      IF (LAJUCV) THEN
        DO JLEV=1,YDCPG_DIM%KFLEVG
          DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            YDVARS%T%T0(JROF,JLEV)=YDAPLPAR_TMP%ADJ%TAUX(JROF,JLEV)
          ENDDO
        ENDDO
        DO JLEV=1,YDCPG_DIM%KFLEVG
          DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            PB1(JROF,ISLB1T9+JLEV-NFLSA)=PB1(JROF,ISLB1T9+JLEV-NFLSA)+YDAPLPAR_TMP%ADJ%DTAJU(JROF,JLEV)
          ENDDO
        ENDDO
      ENDIF
    ELSE
      ! IF (LAJUCV) THEN
      !   missing code under LAJUCV for leap-frog schemes.
      ! ENDIF
    ENDIF

  !    convert to flexible interface structure
  IF (LINTFLEX) THEN
    CALL APLPAR2INTFLEX(YGFL, YDPHY, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%   &
    & DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS, ZDIFEXT, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT  &
    & %DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN,   &
    & YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%  &
    & OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%    &
    & FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%      &
    & FPFPCL, YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQRNG, &
    & YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%FRMH, YDAPLPAR_TMP%MSC%FRMQ,   &
    & YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, YDMF_PHYS%    &
    & OUT%STRDU, YDMF_PHYS%OUT%STRDV, YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV, YDMF_PHYS%OUT%STRMU,  &
    & YDMF_PHYS%OUT%STRMV, YDMF_PHYS%OUT%DIFCQLC, YDMF_PHYS%OUT%DIFCQIC, YDMF_PHYS%OUT%FIMCC,         &
    & YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC,         &
    & YDMF_PHYS%OUT%FCNEGQLC, YDMF_PHYS%OUT%FCNEGQIC, YDMF_PHYS%OUT%FCNEGQRC, YDMF_PHYS%OUT%FCNEGQSC, &
    & YDAPLPAR_TMP%PFL%FTKE, ZTENDPTKE, ZTENDEXT, ZTENDEXT_DEP, YLPROCSET )
  ENDIF

ENDIF ! LLDIAB

!        2.3  Computes MOCON in the CLP.
!             --------------------------
IF (LLDIAB) THEN
  CALL MF_PHYS_MOCON (YDCPG_DIM, YDMF_PHYS, YDMF_PHYS_TMP, YDAPLPAR_TMP, YLMF_PHYS_STATE)
ENDIF

! Store surface water flux P and E for water conservation
IF (LCORWAT) THEN
  CALL MF_PHYS_CORWAT (YDCPG_DIM, YDMF_PHYS, YDMF_PHYS_SURF)
ENDIF

!        2.4  Stores radiation coefficients.
!             ------------------------------

! * writes grid-point transmission coefficients for simplified physics.

IF (LLDIAB) THEN

  IF (LRCOEF.AND.(YDCPG_DIM%KSTEP == 1)) THEN
    IFIELDSS=NG3SR*YDCPG_DIM%KFLEVG
    CALL WRRADCOEF(YDGEOMETRY,YDRCOEF,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KSTGLO,IFIELDSS, &
  & YDAPLPAR_TMP%RDT%COR  , YDAPLPAR_TMP%RDT%RAB3C, YDAPLPAR_TMP%RDT%RAB3N, YDAPLPAR_TMP%RDT%RAB4C, &
  & YDAPLPAR_TMP%RDT%RAB4N, YDAPLPAR_TMP%RDT%RAB6C, YDAPLPAR_TMP%RDT%RAB6N, YDAPLPAR_TMP%RDT%RAT1C, &
  & YDAPLPAR_TMP%RDT%RAT1N, YDAPLPAR_TMP%RDT%RAT2C, YDAPLPAR_TMP%RDT%RAT2N, YDAPLPAR_TMP%RDT%RAT3C, &
  & YDAPLPAR_TMP%RDT%RAT3N, YDAPLPAR_TMP%RDT%RAT4C, YDAPLPAR_TMP%RDT%RAT4N, YDAPLPAR_TMP%RDT%RAT5C, &
  & YDAPLPAR_TMP%RDT%RAT5N, ZAC_HC)
  ENDIF

ENDIF

!       2.5   Ozone
!             -----

IF (LLDIAB) THEN

  IF (LOZONE) THEN
    ! * Caution: this part has not been yet validated relative
    !   to the GFL implementation, and LOZONE (the setup of
    !   which has not yet been updated) can be true only if
    !   the GFL ozone is activated as a prognostic and advected
    !   variable.
    IPO3=(YO3%MP_SL1-1)*(YDCPG_DIM%KFLEVG+2*NFLSUL)
    IF (LSLAG) THEN
      IF (LSLAG) CALL CP_PTRSLB1(YDMODEL%YRML_DYN%YRDYN,YDPTRSLB1,ISLB1U9,ISLB1V9,ISLB1T9,ISLB1VD9,ISLB1GFL9)
      IF (LTWOTL) THEN
        CALL CPOZO (YDCPG_DIM%KLON,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,PDTPHY,YDMF_PHYS%OUT%FCHOZ,&
         & PB1(1,ISLB1GFL9+IPO3+1-NFLSA),YDCPG_PHY0%XYB%RDELP)  
      ELSE
        CALL CPOZO (YDCPG_DIM%KLON,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,PDTPHY,YDMF_PHYS%OUT%FCHOZ,&
         & PB1(1,ISLB1GFL9+IPO3+1-NFLSA),YDCPG_PHY9%XYB%RDELP)  
      ENDIF
    ELSE
      CALL CPOZO (YDCPG_DIM%KLON,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,PDTPHY,YDMF_PHYS%OUT%FCHOZ,&
       & YDVARS%O3%T1,YDCPG_PHY9%XYB%RDELP)
    ENDIF
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
  CALL CPQSOL(YDGEOMETRY%YRDIMV,YDPHY,YDCPG_DIM%KLON,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_PHY0%PREHYD,YDMF_PHYS_SURF%GSP_RR%PT_T0,YDCPG_MISC%QS,YDAPLPAR_TMP%FLU%QSATS,YDCPG_MISC%QSOL)
ENDIF

!        2.7  Computation of tendencies T,u,v and Q.
!             --------------------------------------

! Set GFL tendencies to 0

ZTENDGFL(:,:,:) = 0.0_JPRB

IF (LLDIAB) THEN

  TSPHY = MAX(PDTPHY,1.0_JPRB)

  ! * CPTEND+CPUTQY = Old( CPATY + CPDUP + CPDTHP )
  ! Calcul des tendances de T , U et de Q et modifications
  ! eventuelles de W et de OMEGA/P

  IF (LINTFLEX.AND.(.NOT.LDCONFX)) THEN
    CALL CPTEND_FLEX( YDLDDH, YDMDDH, YGFL, YDPHY, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, YDGSGEOM%GNORDL,      &
    & YDGSGEOM%GNORDM, YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,       &
    & YLMF_PHYS_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_STATE%U, YLMF_PHYS_STATE%V, YLMF_PHYS_STATE%T,     &
    & YLMF_PHYS_STATE%YGSP_RR%T, PGFL, YLPROCSET, YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDH, &
    & ZTENDGFL, YDMF_PHYS%OUT%FHSCL, YDMF_PHYS%OUT%FHSCN, YDMF_PHYS%OUT%FHSSL, YDMF_PHYS%OUT%FHSSN, &
    & YDMF_PHYS%OUT%FHPCL, YDMF_PHYS%OUT%FHPCN, YDMF_PHYS%OUT%FHPSL, YDMF_PHYS%OUT%FHPSN, &
    & PFHP=YDAPLPAR_TMP%MSC%FHP, PFP=YDAPLPAR_TMP%PFL%FP, PFEPFP=YDMF_PHYS%OUT%FEPFP,&
    & PFCMPCQ=YDMF_PHYS%OUT%FCMPCQ, PFCMPSN=YDMF_PHYS%OUT%FCMPSN, PFCMPSL=YDMF_PHYS%OUT%FCMPSL,      &
    & YDDDH=YDDDH)
  ELSE
    CALL CPTEND_NEW( YDMODEL, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, YDGSGEOM%GNORDL,YDGSGEOM%GNORDM,&
     & YDMF_PHYS%OUT%DIFCQ , YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS ,ZDIFEXT,&
     & YDMF_PHYS%OUT%DIFTQ , YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTS ,&
     & YDMF_PHYS%OUT%FCCQL , YDMF_PHYS%OUT%FCCQN , YDMF_PHYS%OUT%FCSQL , YDMF_PHYS%OUT%FCSQN ,&
     & YDMF_PHYS%OUT%FPLSL , YDMF_PHYS%OUT%FPLSN , YDMF_PHYS%OUT%FPLSG , YDMF_PHYS%OUT%FPLCL,  YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLCG,&
     & YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPEVPSG,YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,YDMF_PHYS%OUT%FPEVPCG,&
     & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPFPCL ,YDMF_PHYS%OUT%FPFPCN ,&
     & YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FCQGNG ,&
     & YDMF_PHYS%OUT%FCQNG , YDMF_PHYS%OUT%FRMH  , YDAPLPAR_TMP%MSC%FRMQ  , YDMF_PHYS%OUT%FRSO  , YDMF_PHYS%OUT%FRTH  ,&
     & YDMF_PHYS%OUT%STRCU , YDMF_PHYS%OUT%STRCV , YDMF_PHYS%OUT%STRDU , YDMF_PHYS%OUT%STRDV ,&
     & YDMF_PHYS%OUT%STRTU , YDMF_PHYS%OUT%STRTV , YDMF_PHYS%OUT%STRMU , YDMF_PHYS%OUT%STRMV ,&
     & YDMF_PHYS%OUT%DIFCQLC,YDMF_PHYS%OUT%DIFCQIC,YDMF_PHYS%OUT%FIMCC,&
     & YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC,YDMF_PHYS%OUT%FCNEGQIC,YDMF_PHYS%OUT%FCNEGQRC,YDMF_PHYS%OUT%FCNEGQSC,&
     & YDAPLPAR_TMP%PFL%FTKE, YDAPLPAR_TMP%PFL%FTKEI,  YDAPLPAR_TMP%PFL%FEFB1 , YDAPLPAR_TMP%PFL%FEFB2, YDAPLPAR_TMP%PFL%FEFB3 ,YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP ,&
     & YLMF_PHYS_STATE%YCPG_PHY%XYB%RDELP, YLMF_PHYS_STATE%YCPG_DYN%PHIF  , YLMF_PHYS_STATE%YCPG_DYN%RCP%CP,&
     & YLMF_PHYS_STATE%U,YLMF_PHYS_STATE%V,YLMF_PHYS_STATE%T,&
     & YLMF_PHYS_STATE%Q,YLMF_PHYS_STATE%I,YLMF_PHYS_STATE%L,&
     & YDVARS%LCONV%T0,YDVARS%ICONV%T0,YDVARS%RCONV%T0,YDVARS%SCONV%T0,&
     & YLMF_PHYS_STATE%R,YLMF_PHYS_STATE%S,YLMF_PHYS_STATE%G,&
     & YDAPLPAR_TMP%DSA%CPS, YLMF_PHYS_STATE%YGSP_RR%T  ,&
     & YDMF_PHYS%OUT%FHSCL ,YDMF_PHYS%OUT%FHSCN,YDMF_PHYS%OUT%FHSSL,YDMF_PHYS%OUT%FHSSN,YDMF_PHYS%OUT%FHSSG,&
     & YDMF_PHYS%OUT%FHPCL ,YDMF_PHYS%OUT%FHPCN,YDMF_PHYS%OUT%FHPCG,YDMF_PHYS%OUT%FHPSL,YDMF_PHYS%OUT%FHPSN,YDMF_PHYS%OUT%FHPSG,&
     & YDAPLPAR_TMP%MSC%FHP   ,YDAPLPAR_TMP%PFL%FP   ,  YDMF_PHYS%OUT%FEPFP, YDMF_PHYS%OUT%FCMPCQ, YDMF_PHYS%OUT%FCMPSN, YDMF_PHYS%OUT%FCMPSL,&
     & YDMF_PHYS%OUT%TENDU , YDMF_PHYS%OUT%TENDV , ZTENDU, ZTENDV, ZTENDH ,&
     & ZPTENDQ1,ZPTENDI1,ZPTENDL1,&
     & ZPTENDLCONV1,ZPTENDICONV1,&
     & ZPTENDRCONV1,ZPTENDSCONV1,&
     & ZPTENDR1,ZPTENDS1,ZPTENDG1,ZPTENDTKE1,&
     & ZPTENDEFB11,ZPTENDEFB21,ZPTENDEFB31,&
     & ZTENDEXT,YDDDH)
  ENDIF

  IF (LTWOTL) THEN

  ELSE
    
    IF ( L3MT.OR.LSTRAPRO.OR.(NDPSFI==1)) THEN
!     PFEPFP was ZFEPFP in CPTEND_NEW, before, ZFEPFP still in CPFHPFS
      DO JLEV= 0, YDCPG_DIM%KFLEVG 
        DO JROF = 1, YDCPG_DIM%KLON
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
      CALL CPCHET (YDMF_PHYS, YLMF_PHYS_STATE, YDCPG_MISC, YDRIP, YDPHY, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA,        &
      & YDCPG_DIM%KFLEVG, YDCPG_DIM%KSTEP, YDAPLPAR_TMP%MSC%FRMQ, YDAPLPAR_TMP%DSA%CPS, ZTENDH, ZPTENDQ1, &
      & ZPTENDI1, ZPTENDL1, ZPTENDR1, ZPTENDS1, YDGSGEOM%GEMU, YDGSGEOM%GELAM)
    ENDIF
    
    IF (GCHETN%LPROFV)&
     & CALL PROFILECHET(YDGEOMETRY, YDCPG_MISC, YDMF_PHYS, YDMF_PHYS_TMP, YDAPLPAR_TMP, YDCPG_DYN0, YDMF_PHYS_SURF, YDVARS,     &
       & YDSURF, YDDPHY, YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KFDIA, YDGSGEOM%GELAM, YDGSGEOM%GEMU, YDGSGEOM%GM, &
       & YDOROG%OROG, YDGSGEOM%RCORI, YDCSGEOM%RATATH, YDCSGEOM%RATATX)

  ENDIF

ENDIF

!        2.8  Modification of vertical velocities
!             by some physics output when required.
!             -------------------------------------

IF (LLDIAB) THEN

  ! * MODIFICATION DE LA VITESSE VERTICALE ET DE LA TENDANCE DE
  ! PRESSION DE SURFACE SI NDPSFI=1 ( MASSE VARIABLE ).
  ! Ajout de la physique dans l'equation de continuite/Add physics
  ! in continuity equation.

  IF (NDPSFI == 1) THEN
    IF (LSLAG .AND. LTWOTL) THEN
      CALL CPMVVPS(YDVAB,YDCPG_DIM%KLON,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,PDTPHY,&
       & YDAPLPAR_TMP%PFL%FP,YDCPG_PHY0%PREHYD(:,YDCPG_DIM%KFLEVG),YDMF_PHYS%OUT%FEVL,YDMF_PHYS%OUT%FEVN,&
       & YDCPG_DYN0%CTY%EVEL,YDCPG_DYN0%CTY%PSDVBC,PB1(1,MSLB1SP9))
    ELSEIF (LSLAG .AND. (.NOT.LTWOTL)) THEN
      CALL CPMVVPS(YDVAB,YDCPG_DIM%KLON,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,PDTPHY,&
       & YDAPLPAR_TMP%PFL%FP,YDCPG_PHY9%PREHYD(:,YDCPG_DIM%KFLEVG),YDMF_PHYS%OUT%FEVL,YDMF_PHYS%OUT%FEVN,&
       & YDCPG_DYN0%CTY%EVEL,YDCPG_DYN0%CTY%PSDVBC,PB1(1,MSLB1SP9))
    ELSE
      CALL CPMVVPS(YDVAB,YDCPG_DIM%KLON,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,PDTPHY,&
       & YDAPLPAR_TMP%PFL%FP,YDCPG_PHY9%PREHYD(:,YDCPG_DIM%KFLEVG),YDMF_PHYS%OUT%FEVL,YDMF_PHYS%OUT%FEVN,&
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
     IPGFL(YCOMP(JGFL)%MP1) = (YCOMP(JGFL)%MP_SL1-1)*(YDCPG_DIM%KFLEVG+2*NFLSUL)
  ENDIF   
ENDDO  

!  ALARO does not respect the coding rules, tendency of pseudo-TKE is computed in APLPAR and not
!  in CPTEND_NEW. To use the new version of cputqy it is then necessary to write it in GFL tendencies array.
! This memory transfer is not necessary, please respect coding rules to avoid it.

! Not necessary for intflex: already done in aplpar2intflex
IF (.NOT.(LINTFLEX.AND.(.NOT.LDCONFX))) THEN
  IF (LPTKE) THEN
    DO JLEV=1,YDCPG_DIM%KFLEVG
      DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZPTENDTKE1(JROF,JLEV) = ZTENDPTKE(JROF,JLEV)
      ENDDO
    ENDDO    
  ENDIF
  ! Extra-GFL
  IF(LMDUST.AND.(NGFL_EXT/=0)) THEN
    DO JGFL=1, NGFL_EXT
      DO JLEV=1,YDCPG_DIM%KFLEVG
        DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZTENDGFL(JROF,JLEV,YEXT(JGFL)%MP1) = ZTENDEXT(JROF,JLEV,JGFL)+&! turbulent tendency
                                             & ZTENDEXT_DEP(JROF,JLEV,JGFL) ! moist tendency
        ENDDO
      ENDDO 
    ENDDO   
  ENDIF
ENDIF

! ky: non-zero option not yet coded for the time being.
ZTENDD=0.0_JPRB

IF (LLDIAB) THEN
  ! Calcul de T , Q et du Vent a l'instant 1

  IF (LSLAG) CALL CP_PTRSLB1(YDMODEL%YRML_DYN%YRDYN,YDPTRSLB1,ISLB1U9,ISLB1V9,ISLB1T9,ISLB1VD9,ISLB1GFL9)
  CALL CPUTQY(YDGEOMETRY%YRDIMV, YDGMV, YGFL, YDPTRSLB1, YDPHY, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, &
  & YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, PDTPHY,   &
  & IPGFL, ISLB1T9, ISLB1U9, ISLB1V9, ISLB1VD9, ISLB1GFL9, ZTENDH, ZTENDT, YDMF_PHYS%OUT%TENDU,      &
  & YDMF_PHYS%OUT%TENDV, ZTENDU, ZTENDV, ZTENDD, ZTENDGFL, YLMF_PHYS_STATE%YCPG_DYN%RCP%CP,          &
  & YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_STATE%T, YLMF_PHYS_STATE%U, YLMF_PHYS_STATE%V, PB1, &
  & PGMVT1, PGFLT1, YDMF_PHYS%OUT%FDIS)

ENDIF

IF (LLDIAB) THEN
  CALL MF_PHYS_FPL_PART2 (YDCPG_DIM, YDAPLPAR_TMP, YDVARS, YDMODEL)
ENDIF

!       2.9b Prognostic convection etc.
!            --------------------------

! TRANSFER NOT ADVECTED VARIABLES INTO PGFLT1
CALL MF_PHYS_TRANSFER (YDCPG_DIM, YDVARS, YDMODEL)

!        2.10  Surface variables.
!              ------------------

IF (LLDIAB.AND.(.NOT.LSFORCS)) THEN
  
  IF (.NOT.LMSE) THEN
    DO JLEV=0,YDCPG_DIM%KFLEVG
      DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        YDAPLPAR_TMP%PFL%FPLSN(JROF,JLEV)=YDMF_PHYS%OUT%FPLSN(JROF,JLEV)+YDMF_PHYS%OUT%FPLSG(JROF,JLEV)
      ENDDO
    ENDDO
    CALL CPTENDS( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, YSP_SBD%NLEVS, PDTPHY, YDMF_PHYS% &
    & OUT%FPLCL, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLCN, YDAPLPAR_TMP%PFL%FPLSN, YDMF_PHYS&
    & %OUT%FRSO, YDMF_PHYS%OUT%FRTH, YDMF_PHYS_SURF%GSP_SG%PA_T1, YDMF_PHYS%OUT%CT, YDAPLPAR_TMP%DSA%C1, &
    & YDAPLPAR_TMP%DSA%C2, YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FCLL, YDMF_PHYS&
    & %OUT%FCLN, YDMF_PHYS%OUT%FCS, YDAPLPAR_TMP%FLU%FEVI, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT&
    & %FEVN, YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FGEL, YDMF_PHYS%OUT%FGELS, YDMF_PHYS%OUT%FLWSP,      &
    & YDMF_PHYS%OUT%FONTE, YDMF_PHYS%OUT%FTR, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSP_SG%    &
    & PR_T1, YDMF_PHYS%OUT%RUISL, YDMF_PHYS%OUT%RUISP, YDMF_PHYS%OUT%RUISS, YDMF_PHYS_SURF%GSP_SG%  &
    & PF_T1, YDAPLPAR_TMP%FLU%VEG, YDAPLPAR_TMP%TDS%TDTS, YDAPLPAR_TMP%TDS% &
    & TDTP, YDAPLPAR_TMP%TDS%TDWS, YDAPLPAR_TMP%TDS%TDWSI, YDAPLPAR_TMP%TDS%&
    & TDWP, YDAPLPAR_TMP%TDS%TDWPI, YDAPLPAR_TMP%TDS%TDWL, YDAPLPAR_TMP%TDS%&
    & TDSNS, YDAPLPAR_TMP%TDS%TDALBNS, YDAPLPAR_TMP%TDS%TDRHONS)  

    CALL CPWTS(YDSURF, YDMODEL%YRML_AOC%YRMCC, YDPHY, YDMODEL%YRML_PHY_MF%YRPHY1, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA,  &
    & YSP_SBD%NLEVS, PDTPHY, YDAPLPAR_TMP%TDS%TDTS, YDAPLPAR_TMP%TDS%TDTP, YDAPLPAR_TMP%TDS%TDWS, YDAPLPAR_TMP%TDS%TDWSI, YDAPLPAR_TMP%TDS%TDWP, &
    & YDAPLPAR_TMP%TDS%TDWPI, YDAPLPAR_TMP%TDS%TDWL, YDAPLPAR_TMP%TDS%TDSNS,           &
    & YDAPLPAR_TMP%TDS%TDALBNS, YDAPLPAR_TMP%TDS%TDRHONS, YDMF_PHYS_SURF%GSD_VP%PTPC, &
    & YDMF_PHYS_SURF%GSD_VP%PWPC, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PIVEG,            &
    & YDMF_PHYS_SURF%GSP_RR%PT_T1, YDMF_PHYS_SURF%GSP_SB%PT_T1, YDMF_PHYS_SURF%GSP_RR%PW_T1,          &
    & YDMF_PHYS_SURF%GSP_RR%PIC_T1, YDMF_PHYS_SURF%GSP_SB%PQ_T1, YDMF_PHYS_SURF%GSP_SB%PTL_T1,        &
    & YDMF_PHYS_SURF%GSP_RR%PFC_T1, YDMF_PHYS_SURF%GSP_SG%PF_T1, YDMF_PHYS_SURF%GSP_SG%PA_T1,         &
    & YDMF_PHYS_SURF%GSP_SG%PR_T1)  
  ELSE
    IF (LLXFUMSE) THEN
      DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        YDMF_PHYS_SURF%GSP_RR%PT_T0(JROF)=PGPAR(JROF,MVTS)
      ENDDO
    ELSE
      DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        YDMF_PHYS_SURF%GSP_RR%PT_T1(JROF)=PGPAR(JROF,MVTS)
      ENDDO
    ENDIF
  ENDIF
  IF(LNUDG)THEN

    ! * Calculation of IPQ since the old pointers
    !   MSLB1[X]9 (=MSLB1GFL9+IP[X]) do not exist any longer in PTRSLB1.
    IPQ=(YQ%MP_SL1-1)*(YDCPG_DIM%KFLEVG+2*NFLSUL)

    IF (LSLAG) CALL CP_PTRSLB1(YDMODEL%YRML_DYN%YRDYN,YDPTRSLB1,ISLB1U9,ISLB1V9,ISLB1T9,ISLB1VD9,ISLB1GFL9)
    CALL CPNUDG ( YDCPG_DIM%KLON, YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA, NFNUDG, YDCPG_DIM%KFLEVG, YDCPG_DIM%KBL,&
     & XPNUDG,&
     & YDMF_PHYS_SURF%GSD_VF%PNUDM,&
     & YDMF_PHYS_SURF%GSP_RR%PT_T1,YDMF_PHYS_SURF%GSP_RR%PW_T1,&
     & YDMF_PHYS_SURF%GSP_SB%PQ_T1,YDMF_PHYS_SURF%GSP_SG%PF_T1,&
     & PB1(1,ISLB1T9+1-NFLSA),PB1(1,ISLB1GFL9+IPQ+1-NFLSA),&
     & PB1(1,ISLB1U9+1-NFLSA),PB1(1,ISLB1V9+1-NFLSA),&
     & PB1(1,MSLB1SP9),&
     & YDVARS%T%T0,YDVARS%Q%T0,YDVARS%U%T0,&
     & YDVARS%V%T0,YDCPG_PHY0%PREHYD(:,YDCPG_DIM%KFLEVG),YDGSGEOM%GM,YDMF_PHYS_SURF%GSD_VF%PLSM)
  ENDIF
ENDIF


CALL MF_PHYS_CVV (YDCPG_DIM, YDVARS, YDMODEL)

!        3.3  Store the model trajectory at t-dt (leap-frog) or t (sl2tl).
!             ------------------------------------------------------------

IF (LTRAJPS) THEN
  PTRAJ_PHYS%PQSSMF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_MISC%QS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
  PTRAJ_PHYS%PTSMF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA) =YLMF_PHYS_STATE%YGSP_RR%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
  PTRAJ_PHYS%PSNSMF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YLMF_PHYS_STATE%YGSP_SG%F(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)

  IF (.NOT. LTWOTL) THEN
    CALL WRPHTRAJM(YDGEOMETRY,YDSIMPHL,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,PTRAJ_PHYS,&
     & YDVARS%U%T9,YDVARS%V%T9,YDVARS%T%T9,&
     & YDVARS%Q%T9,YDVARS%L%T9,YDVARS%I%T9,YDVARS%SP%T9)  
  ENDIF

  IF (LPRTTRAJ.AND.PTRAJ_PHYS%LASTCHUNK) WRITE(NULOUT,*)'GREPTRAJ STORE TRAJ_PHYS in MF_PHYS_PAR'
ENDIF

!     ------------------------------------------------------------------

!*       5.    Final calculations.
!              -------------------

! * Restore the initial value of some pseudo-historical surface buffers
!   if relevant.
IF (LLDIAB) THEN
  IF (LL_SAVE_PHSURF) THEN
    CALL MF_PHYS_SAVE_PHSURF_PART2 (YDCPG_DIM, YDMF_PHYS_TMP, YDMF_PHYS_SURF, YDVARS, YDSURF, YDMODEL)
  ENDIF
ENDIF

IF (LLDIAB) THEN

  ! Store horizontal exchange coefficients (3D turbulence) to SL2 buffers
  IF (L3DTURB) THEN
    DO JLEV=1,YDCPG_DIM%KFLEVG
      PB2(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,MSLB2KAPPAM+JLEV-1)=YDAPLPAR_TMP%KUR%KUROV_H(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)
      PB2(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,MSLB2KAPPAH+JLEV-1)=YDAPLPAR_TMP%KUR%KTROV_H(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)
    ENDDO
  ENDIF

ENDIF

CALL MF_PHYS_BAYRAD (YDAPLPAR_TMP, YDVARS, YDMODEL, LAROME)


!-------------------------------------------------
! Extract Single Column Model profiles from 3D run or 
! write LFA file for MUSC (1D model)
!-------------------------------------------------
IF(LGSCM.OR.LMUSCLFA) THEN
  IF (LAROME) THEN
    DO JLEV=1,YDCPG_DIM%KFLEVG
      DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        YDCPG_MISC%NEB(JROF,JLEV)=YDVARS%A%T1(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
  CALL WRITEPHYSIO(YDGEOMETRY, YDCPG_MISC, YDCPG_PHY0, YDCPG_PHY9, YDMF_PHYS, YDCPG_DYN0, YDCPG_DYN9, &
  & YDMF_PHYS_SURF, YDVARS, YDSURF, YDDPHY, YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KFDIA, YDCPG_DIM%KIDIA, YDCPG_DIM%KGL1, YDCPG_DIM%KGL2,        &
  & YDCPG_DIM%KSTGLO, YDCPG_DIM%KSTEP, NTSSG, YSP_SBD%      NLEVS, YDGSGEOM%GELAM, YDGSGEOM%GEMU, YDGSGEOM%GM, YDOROG%    &
  & OROG, YDGSGEOM%RCORI, YDCSGEOM%RATATH, YDCSGEOM%RATATX, YDGSGEOM%       GECLO, YDGSGEOM%GESLO, YDMF_PHYS_TMP, YDAPLPAR_TMP  )
ENDIF

IF (LEDR) THEN
  YDMF_PHYS_SURF%GSD_DI%PXEDR(:,:)=YDMF_PHYS%OUT%EDR(:,:)
ENDIF

CALL MF_PHYS_PRECIPS (YDCPG_DIM, YDMF_PHYS_TMP, YDMF_PHYS_SURF, YDMODEL)

! Restore Tt and grad(Tt) for NHQE model.
IF (LNHQE) THEN
  CALL MF_PHYS_NHQE_PART2 (YDGEOMETRY, YDCPG_DIM, YDMF_PHYS_TMP, YDVARS)
ENDIF

!     ------------------------------------------------------------------

!       6. destructor for procset
IF (LINTFLEX) CALL CLEANINTPROCSET(YLPROCSET)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('MF_PHYS_PAR',1,ZHOOK_HANDLE)
END SUBROUTINE MF_PHYS_PAR
