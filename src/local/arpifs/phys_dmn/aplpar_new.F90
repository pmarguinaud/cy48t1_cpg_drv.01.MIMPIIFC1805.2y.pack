!OPTIONS XOPT(NOEVAL)
SUBROUTINE APLPAR_NEW(YDGEOMETRY,YDCPG_DIM,YDMF_PHYS_STATE,YDCPG_DYN0,YDCPG_MISC,YDMF_PHYS,YDMF_PHYS_SURF,YDVARS,YDSURF, YDXFU,YDCFU,YDMODEL,&
 & KTDIA  , &
 & KVCLIS , KVCLIV , &
 & KSGST  , KCSS   ,&
 & KNFRRC , PDT,&
 !---------------------------------------------------------------------
 ! - INPUT .
 & PINDX  , PINDY  ,&
 & LDXFUMSE,&
 & PRCORI ,&
 & PGFL, &
 & PKOZO  ,&
 & PGPAR  , &
 & PGELAM ,&
 & PGEMU  , PGM    , POMPAC , PAC    , &
 & POROG  ,&
 ! - OUTPUT .
 & PDIFSV , &
 & KCLPH  , &
 & PEZDIAG,PTENDPTKE,&
 & PTENDEXT_DEP, PTRAJ_PHYS,PTDISS,YDDDH,&
 & PFTCNS, &
 & PGP2DSPP )

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
! KTDIA : DEBUT BOUCLE VERTICALE DANS LA PHYSIQUE.
! KTDIA : START OF THE VERTICAL LOOP IN THE PHYSICS (IF SOME LEVELS ARE
!                     SKIPPED AT THE TOP OF THE MODEL).
! KLEV : FIN BOUCLE VERTICE ET DIMENSION VERTICALE (NFLEVG DANS CPG).
! KLEV : END OF VERTICAL LOOP AND VERTICAL DIMENSION(NFLEVG IN *CPG*).
! KSTGLO     : global offset.
! KVCLIS     : NOMBRE DE CHAMPS POUR LA PHOTO-CHIMIE DE L'OZONE.
!                     (NVCLIS DANS CPG).
! KVCLIS     : NUMBER OF FIELDS FOR OZONE PHOTO-CHEMISTRY.
!                     (NVCLIS IN CPG).
! KVCLIV     : NOMBRE DE CHAMPS POUR LA VEGETATION (NVCLIV DANS CPG).
! KVCLIV     : NUMBER OF FIELDS FOR VEGETATION (NVCLIV IN CPG).
! KSGST      : NOMBRE DE TEMPERATURES ET DE FLUX DE SURFACE SOUS-MAILLE
!                     (NTSSG DANS CPG)
! KSGST      : NUMBER OF SUBGRID SURFACE TEMPERATURES AND FLUXES
!                     (NTSSG IN *CPG*)
! KCSS       : NBRE DE NIVEAUX DANS LE SOL PROFOND
! KCSS       : NBR OF VERTICAL LAYERS IN THE DEEP SOIL
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
USE MF_PHYS_TYPE_MOD   , ONLY : MF_PHYS_TYPE
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

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_DIM_TYPE),INTENT(IN)    :: YDCPG_DIM
TYPE (MF_PHYS_STATE_TYPE),INTENT(INOUT):: YDMF_PHYS_STATE
TYPE(CPG_DYN_TYPE), INTENT(INOUT) :: YDCPG_DYN0
TYPE(CPG_MISC_TYPE),INTENT(INOUT):: YDCPG_MISC
TYPE(MF_PHYS_TYPE),INTENT(INOUT) :: YDMF_PHYS
TYPE(MF_PHYS_SURF_TYPE),INTENT(INOUT):: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
TYPE(TXFU)        ,INTENT(INOUT) :: YDXFU
TYPE(TCFU)        ,INTENT(INOUT) :: YDCFU
TYPE(MODEL)       ,INTENT(INOUT),TARGET :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KNFRRC
INTEGER(KIND=JPIM),INTENT(IN)    :: KVCLIS
INTEGER(KIND=JPIM),INTENT(IN)    :: KSGST
INTEGER(KIND=JPIM),INTENT(IN)    :: KCSS
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KVCLIV
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT
LOGICAL           ,INTENT(IN)    :: LDXFUMSE
REAL(KIND=JPRB)   ,INTENT(IN)    :: PINDX(YDCPG_DIM%KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PINDY(YDCPG_DIM%KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCORI(YDCPG_DIM%KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKOZO(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,KVCLIS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGPAR(YDCPG_DIM%KLON,YDMODEL%YRML_PHY_MF%YRPARAR%NGPAR+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(YDCPG_DIM%KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(YDCPG_DIM%KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGM(YDCPG_DIM%KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(YDCPG_DIM%KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POMPAC(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,0:YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAC(YDCPG_DIM%KFLEVG+1,YDCPG_DIM%KFLEVG+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDIFSV(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,1:YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)

! Prognostic convection variables (MUST BE 1:KLEV for GFL)
! Unless mf_phys has copied them into a 0:KLEV array
! Up to now, we leave 1:KLEV for pentch: see later.
!---------------------------------------------------

INTEGER(KIND=JPIM),INTENT(OUT)   :: KCLPH(YDCPG_DIM%KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFTCNS(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,6)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEZDIAG(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDPTKE(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENDEXT_DEP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTDISS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
TYPE (TRAJ_PHYS_TYPE), INTENT(INOUT) :: PTRAJ_PHYS
TYPE(TYP_DDH), INTENT(INOUT)     :: YDDDH
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGP2DSPP(YDCPG_DIM%KLON,YSPP%N2D)

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


 ! arrays for 3D turb
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
   ! updraught envt vert vel*dt
 ! fall velocity of rain
 ! fall velocity of snow
 ! fall velocity of graupel
REAL(KIND=JPRB) :: ZICEFR1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)! Resolved Condensate ice fraction
REAL(KIND=JPRB) :: ZRHCRI(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Smith scheme critical RH
REAL(KIND=JPRB) :: ZRHDFDA(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)! RK scheme change in RH over cloud
REAL(KIND=JPRB) :: ZLHS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)   ! Sublimation latent heat
REAL(KIND=JPRB) :: ZLHV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)   ! Evaporation latent heat
    ! Temporar storage for updated PLH
 ! Temporar storage for updated PLSCPE
  ! Temporar storage for updated PQSAT
REAL(KIND=JPRB) :: ZQSATS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! QSAT of resolved cond./evap. scheme
    ! Temporar storage for updated PQW
REAL(KIND=JPRB) :: ZRH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! Temporar storage for updated PRH
REAL(KIND=JPRB) :: ZTW(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! Temporar storage for updated PTW)
    ! Saturation departure for a given thermodynamic state
   ! maximum saturation departure
REAL(KIND=JPRB) :: ZPOID(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)  ! DP/(RG*DT) FOR A GIVEN LEVEL AND A GIVEN TIME STEP.
REAL(KIND=JPRB) :: ZIPOI(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)  ! INVERSE OF ZPOID.

 ! Updraught Specific moisture
 ! Updraught Temperature
 ! Updraught zonal wind
 ! Updraught merid. wind

REAL(KIND=JPRB) :: ZTMIC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Temperature for microphysics
REAL(KIND=JPRB) :: ZQMIC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Specific moisture for microphysics

     ! updated temperature T for cascading parameterization
 ! temperature corr. for convective cloud
 ! net melting (-freezing) rate of ice
! net melting (-freezing) rate of graupel
     ! updated zonal velocity
     ! updated meridional velocity

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
 ! downdraft flux of specific humidity
 ! downdraft flux of liquid water
 ! downdraft flux of  solid water
REAL(KIND=JPRB) :: ZSEDIQL(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! sedimentation flux of cloud liquid water
REAL(KIND=JPRB) :: ZSEDIQI(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! sedimentation flux of cloud ice water
 ! downdraft entalphy flux
 ! change in horizontal mom.
 ! change in horizontal mom.
 ! degree of inhomogeneity in precips.
         ! Convective precipit mesh fraction
         ! Precipitation mesh fraction
        ! Precipitation auxilary
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
 ! Flux de masse (updraft) pour XIOS output
 ! Flux de masse (downdraft) pour XIOS output

REAL(KIND=JPRB) :: ZCONDCVPPL(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de condensation liquide du a CVVPP (KFB)
REAL(KIND=JPRB) :: ZCONDCVPPI(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de condensation glace du a CVVPP (KFB)
REAL(KIND=JPRB) :: ZPRODTH_CVPP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de production thermique de TKE du a CVPP(KFB)
REAL(KIND=JPRB) :: ZDTRAD(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! radiation contribution to T tendency
REAL(KIND=JPRB) :: ZDQVDIFF(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! turtb.diff contribution to Qv tendency
REAL(KIND=JPRB) :: ZRKQCTEND(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Qc input for RK condensation scheme
REAL(KIND=JPRB) :: ZRKQVTEND(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Qv input for RK condensation scheme
REAL(KIND=JPRB) :: ZRKTTEND(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! T input for RK condensation scheme
REAL(KIND=JPRB) :: ZDQV, ZDQI, ZDQL, ZDQR, ZDQS, ZDQC, ZGDT, ZGDTI,&
                 & ZQV0
REAL(KIND=JPRB) :: ZTMPAF(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! temporary array for Add_Field_3d.
REAL(KIND=JPRB) :: ZTMPAS(YDCPG_DIM%KLON)         ! temporary array for Add_Field_2d..
 ! temporary array

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
 & ZHQ(YDCPG_DIM%KLON),ZHU(YDCPG_DIM%KLON),ZHTR(YDCPG_DIM%KLON),ZCDNH(YDCPG_DIM%KLON),&
 & ZRTI(YDCPG_DIM%KLON),ZDPHI(YDCPG_DIM%KLON),ZPRS(YDCPG_DIM%KLON),ZSTAB(YDCPG_DIM%KLON),ZTAUX(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZWFC(YDCPG_DIM%KLON),ZWPMX(YDCPG_DIM%KLON),ZWLMX(YDCPG_DIM%KLON),ZWSEQ(YDCPG_DIM%KLON),&
 & ZWSMX(YDCPG_DIM%KLON),ZWWILT(YDCPG_DIM%KLON),&
 & ZC3(YDCPG_DIM%KLON),ZCG(YDCPG_DIM%KLON),ZCN(YDCPG_DIM%KLON),&
 & ZNEIJG(YDCPG_DIM%KLON),ZNEIJV(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZPCLS(YDCPG_DIM%KLON)

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
REAL(KIND=JPRB) :: ZTHETAVS(YDCPG_DIM%KLON)
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


!           2-D ARRAY FOR SIMPL.RADIATION SCHEME



INTEGER(KIND=JPIM) :: IJN, JCHA, JLEV, JLON, JSG
 ! useful size of klon arrays for mesonh physics
LOGICAL :: LLCLS, LLHMT, LLREDPR

REAL(KIND=JPRB) :: ZAEN, ZAEO, ZCARDI, ZEPS, ZEPS0, ZEPSNEB, ZEPSO3
REAL(KIND=JPRB) :: ZEPSA

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
 & ZCFBU(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZCFBV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

REAL(KIND=JPRB) :: ZDSE(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZFEV(YDCPG_DIM%KLON),ZFMDU(YDCPG_DIM%KLON),ZFMDV(YDCPG_DIM%KLON),ZFEVS(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZSRAIN(YDCPG_DIM%KLON), ZSSNOW(YDCPG_DIM%KLON), ZSGROUPEL(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZSFCO2(YDCPG_DIM%KLON), ZRHODREFM(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZDEPTH_HEIGHT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZZS(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZTSN(YDCPG_DIM%KLON)
REAL(KIND=JPRB)   :: ZBUDTH (YDCPG_DIM%KLON)
REAL(KIND=JPRB)   :: ZBUDSO (YDCPG_DIM%KLON)
REAL(KIND=JPRB)   :: ZFCLL  (YDCPG_DIM%KLON)
! FOR Hv
REAL(KIND=JPRB)   :: ZHV2(YDCPG_DIM%KLON)
! FOR DUST

REAL(KIND=JPRB) :: ZCFASV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,1:YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) ! SCOND MEMBRE POUR LES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZCFBSV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,1:YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) ! SCOND MEMBRE POUR LES SCALAIRES PASSIFS
 ! SCOND MEMBRE POUR LES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZINVG



REAL(KIND=JPRB), DIMENSION (:,:,:), ALLOCATABLE  :: ZSVM, ZPSV
REAL(KIND=JPRB), DIMENSION (:,:),   ALLOCATABLE  :: ZSFSV  ! passifs scalaires surf flux
! TRAITEMENT DES SCALAIRES PASSIFS

REAL(KIND=JPRB), DIMENSION (YDCPG_DIM%KLON,1,YDCPG_DIM%KFLEVG) :: ZZZ,ZDZZ


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

 
 


  
REAL(KIND=JPRB) :: ZFLX_LOTT_GWU(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG), ZFLX_LOTT_GWV(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)


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

REAL(KIND=JPRB)    :: ZRVMD
LOGICAL            :: LLAERO, LLCALLRAD
REAL(KIND=JPRB)    :: ZAIPCMT(YDCPG_DIM%KLON) ! Activity Index of PCMT: 1. if PCMT is active, 0. else case.


REAL(KIND=JPRB)    :: ZQIC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZQLC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZQRC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZQSC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)    :: ZQLI_CVP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)    :: ZQC_DET_PCMT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)    :: ZFPLS(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZFPLC(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)




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
           ! proportion of Lambertian scattering
 ! direct (parallel) surface albedo
REAL(KIND=JPRB) :: ZCLCT_RAD(YDCPG_DIM%KLON) ! total cloud cover for radiation
REAL(KIND=JPRB) :: ZDECRD   (YDCPG_DIM%KLON) ! decorrelation depth for cloud overlaps
REAL(KIND=JPRB) :: ZDECRD_MF(YDCPG_DIM%KLON) ! decorrelation depth for cloud overlaps
                                   ! in microphysics

REAL(KIND=JPRB) :: ZQG(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)


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
REAL(KIND=JPRB) :: ZGP2DSPP(YDCPG_DIM%KLON,1),ZLUDELI(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,4),ZLRAIN(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZRSUD(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,2)
REAL(KIND=JPRB), ALLOCATABLE :: ZCEN(:,:,:),ZSCAV(:)

! Precipitation type diagnostics
!--------------------------------
REAL(KIND=JPRB)   :: ZFPLSG(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)



REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "acclph.intfb.h"




#include "cucalln_mf.intfb.h"


#include "acdnshf.intfb.h"


#include "acdrag.intfb.h"
#include "acdrme.intfb.h"
#include "acdrov.intfb.h"
#include "acfluso.intfb.h"
#include "achmt.intfb.h"




#include "acnebn.intfb.h"





#include "acpluis.intfb.h"
#include "acpluiz.intfb.h"

#include "recmwf.intfb.h"
#include "acdayd.intfb.h"
#include "acrso.intfb.h"


#include "acsol.intfb.h"
#include "actqsat.intfb.h"




#include "acuptq.intfb.h"
#include "acveg.intfb.h"

#include "qngcor.intfb.h"
#include "radaer.intfb.h"
#include "radheat.intfb.h"
#include "radozcmf.intfb.h"
#include "suozon.intfb.h"
#include "actke.intfb.h"
#include "acdifv1.intfb.h"
#include "acdifv2.intfb.h"

#include "achmtls.intfb.h"
#include "aro_ground_param.h"
#include "aro_ground_diag.h"
#include "aro_ground_diag_2isba.h"
#include "aro_ground_diag_z0.h"





#include "acnebcond.intfb.h"

#include "acnpart.intfb.h"
#include "acvppkf.intfb.h"

#include "arp_ground_param.intfb.h"

#include "accldia.intfb.h"










#include "checkmv.intfb.h"
#include "acaa1.intfb.h"
!include "chem_main.intfb.h"

#include "acevadcape.intfb.h"
#include "ppwetpoint.intfb.h"
#include "dprecips.intfb.h"
#include "acvisih.intfb.h"
#include "culight.intfb.h"

!     ------------------------------------------------------------------

#include "fcttrm.func.h"
!!#include "fctdoi.func.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('APLPAR',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDGEM=>YDGEOMETRY%YRGEM, &
 &   YDSTA=>YDGEOMETRY%YRSTA, &
 &  YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, &
 & YDTOPH=>YDMODEL%YRML_PHY_MF%YRTOPH,YDRCOEF=>YDMODEL%YRML_PHY_RAD%YRRCOEF, &
 & YDERDI=>YDMODEL%YRML_PHY_RAD%YRERDI, YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL, &
 & YDMCC=>YDMODEL%YRML_AOC%YRMCC,YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY,  &
 & YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD,YDPHY3=>YDMODEL%YRML_PHY_MF%YRPHY3,  &
 & YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2,YDPHY1=>YDMODEL%YRML_PHY_MF%YRPHY1,YDPHY0=>YDMODEL%YRML_PHY_MF%YRPHY0,  &
 & YGFL=>YDMODEL%YRML_GCONF%YGFL, &
 & YDNORGWD=>YDMODEL%YRML_PHY_MF%YRNORGWD, &
  & YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR,YDMSE=>YDMODEL%YRML_PHY_MF%YRMSE,  &
  & YDPHYDS=>YDMODEL%YRML_PHY_MF%YRPHYDS,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY)

ASSOCIATE(XMINLM=>YDPHY0%XMINLM, RTCAPE=>YDPHY0%RTCAPE, GRSO=>YDPHY0%GRSO, GCVTSMO=>YDPHY0%GCVTSMO, &
 & XKLM=>YDPHY0%XKLM, GAEPS=>YDPHY0%GAEPS, AERCS1=>YDPHY0%AERCS1, &
 & AERCS3=>YDPHY0%AERCS3, AERCS5=>YDPHY0%AERCS5, HUTIL2=>YDPHY0%HUTIL2, &
 & HUTIL1=>YDPHY0%HUTIL1, XMAXLM=>YDPHY0%XMAXLM, &
 & TEQK=>YDPHY0%TEQK, HUCOE=>YDPHY0%HUCOE, &
 & LCVNHD=>YDPHY0%LCVNHD, TEQC=>YDPHY0%TEQC, UHDIFV=>YDPHY0%UHDIFV, &
 & HUTIL=>YDPHY0%HUTIL, NPCLO1=>YDPHY0%NPCLO1, NPCLO2=>YDPHY0%NPCLO2, &
 & RDECRD=>YDPHY0%RDECRD, RDECRD1=>YDPHY0%RDECRD1, RDECRD2=>YDPHY0%RDECRD2, &
 & RDECRD3=>YDPHY0%RDECRD3, RDECRD4=>YDPHY0%RDECRD4, &
 & ETKE_MIN=>YDPHY0%ETKE_MIN, &
 & HSOLIWR=>YDPHY1%HSOLIWR, ALCRIN=>YDPHY1%ALCRIN, ALBMED=>YDPHY1%ALBMED, &
 & WSMX=>YDPHY1%WSMX, LALBMERCLIM=>YDPHY1%LALBMERCLIM, &
 & HSOLIT0=>YDPHY1%HSOLIT0, HSOL=>YDPHY1%HSOL, WPMX=>YDPHY1%WPMX, &
 & EMCRIN=>YDPHY1%EMCRIN, EMMMER=>YDPHY1%EMMMER, TMERGL=>YDPHY1%TMERGL, &
 & RTINER=>YDPHY1%RTINER, EMMGLA=>YDPHY1%EMMGLA, &
 & TSPHY=>YDPHY2%TSPHY, LRAFTKE=>YDPHY2%LRAFTKE, LRAFTUR=>YDPHY2%LRAFTUR, &
 & HVCLS=>YDPHY2%HVCLS, HTCLS=>YDPHY2%HTCLS, &
 & FSM_HH=>YDPHY3%FSM_HH, FSM_GG=>YDPHY3%FSM_GG, FSM_FF=>YDPHY3%FSM_FF, &
 & FSM_EE=>YDPHY3%FSM_EE, FSM_II=>YDPHY3%FSM_II, FSM_CC=>YDPHY3%FSM_CC, &
 & FSM_DD=>YDPHY3%FSM_DD, RLAMB_WATER=>YDPHY3%RLAMB_WATER, RII0=>YDPHY3%RII0, &
 & QCO2=>YDPHY3%QCO2, RLAMB_SOLID=>YDPHY3%RLAMB_SOLID, &
 & NDLUNG=>YDDIM%NDLUNG, NDGUNG=>YDDIM%NDGUNG, &
 & NDLUXG=>YDDIM%NDLUXG, NDGUXG=>YDDIM%NDGUXG, &
 & LRDEPOS=>YDARPHY%LRDEPOS, LMPA=>YDARPHY%LMPA, CCOUPLING=>YDARPHY%CCOUPLING, &
 & LMDUST=>YDARPHY%LMDUST, LMSE=>YDARPHY%LMSE, &
 & LSURFEX_KFROM=>YDARPHY%LSURFEX_KFROM, &
 & YA=>YGFL%YA, NGFL_EZDIAG=>YGFL%NGFL_EZDIAG, YFQTUR=>YGFL%YFQTUR, &
 & YFSTUR=>YGFL%YFSTUR, YIRAD=>YGFL%YIRAD, YLRAD=>YGFL%YLRAD, &
 & XZSEPS=>YDMSE%XZSEPS, &
 & LTRAJPS=>YDSIMPHL%LTRAJPS, LVDIFSPNL=>YDSIMPHL%LVDIFSPNL, &
 & LGWDSPNL=>YDSIMPHL%LGWDSPNL, LRAYSP=>YDSIMPHL%LRAYSP, &
 & LSTRA=>YDPHY%LSTRA, LAEROSOO=>YDPHY%LAEROSOO, LCDDPRO=>YDPHY%LCDDPRO, &
 & LRKCDEV=>YDPHY%LRKCDEV, LHUCN=>YDPHY%LHUCN, LCOEFK_TOMS=>YDPHY%LCOEFK_TOMS, &
 & LVDIF=>YDPHY%LVDIF, LRRMES=>YDPHY%LRRMES, LCVRA=>YDPHY%LCVRA, LCVTDK=>YDPHY%LCVTDK, &
 & LPTKE=>YDPHY%LPTKE, LCOEFK_RIS=>YDPHY%LCOEFK_RIS, LAEROLAN=>YDPHY%LAEROLAN, &
 & LNEBECT=>YDPHY%LNEBECT, LCVPGY=>YDPHY%LCVPGY, LAERODES=>YDPHY%LAERODES, &
 & LNEWSTAT=>YDPHY%LNEWSTAT, LTHERMO=>YDPHY%LTHERMO, LO3FL=>YDPHY%LO3FL, &
 & LRRGUST=>YDPHY%LRRGUST, LPHSPSH=>YDPHY%LPHSPSH, LSNV=>YDPHY%LSNV, &
 & LECSHAL=>YDPHY%LECSHAL, LECT=>YDPHY%LECT, LSTRAPRO=>YDPHY%LSTRAPRO, &
 & LDIFCONS=>YDPHY%LDIFCONS, LNODIFQC=>YDPHY%LNODIFQC, &
 & LAEROVOL=>YDPHY%LAEROVOL, LRSTAER=>YDPHY%LRSTAER, LOZONE=>YDPHY%LOZONE, &
 & NCALLRAD=>YDPHY%NCALLRAD, LNCVPGY=>YDPHY%LNCVPGY, LRAYLU=>YDPHY%LRAYLU, &
 & LAEROSUL=>YDPHY%LAEROSUL, LO3ABC=>YDPHY%LO3ABC, LSTRAS=>YDPHY%LSTRAS, &
 & LCOEFK_PTTE=>YDPHY%LCOEFK_PTTE, LSFHYD=>YDPHY%LSFHYD, &
 & LAEROSEA=>YDPHY%LAEROSEA, NDIFFNEB=>YDPHY%NDIFFNEB, LEDKF=>YDPHY%LEDKF, &
 & LRCVOTT=>YDPHY%LRCVOTT, LNEBR=>YDPHY%LNEBR, LNEBN=>YDPHY%LNEBN, &
 & LMPHYS=>YDPHY%LMPHYS, LZ0HSREL=>YDPHY%LZ0HSREL, &
 & LCAMOD=>YDPHY%LCAMOD, LCOMOD=>YDPHY%LCOMOD, LCONDWT=>YDPHY%LCONDWT, &
 & LCVCSD=>YDPHY%LCVCSD, LNSDO=>YDPHY%LNSDO, LUDEVOL=>YDPHY%LUDEVOL, &
 & LRAY=>YDPHY%LRAY, LGWD=>YDPHY%LGWD, LCOEFKTKE=>YDPHY%LCOEFKTKE, &
 & L3MT=>YDPHY%L3MT, LRAYFM=>YDPHY%LRAYFM, LECDEEP=>YDPHY%LECDEEP, &
 & LCVGQD=>YDPHY%LCVGQD, LGWDC=>YDPHY%LGWDC, LNORGWD=>YDPHY%LNORGWD, LFLUSO=>YDPHY%LFLUSO, &
 & LNEBCO=>YDPHY%LNEBCO, LNEBCV=>YDPHY%LNEBCV, LSOLV=>YDPHY%LSOLV, &
 & LCOEFKSURF=>YDPHY%LCOEFKSURF, LRNUEXP=>YDPHY%LRNUEXP, &
 & NRAY=>YDPHY%NRAY, LDAYD=>YDPHY%LDAYD, &
 & LEDMFI=>YDPHY%LEDMFI, LGPCMT=>YDPHY%LGPCMT, LFPCOR=>YDPHY%LFPCOR, &
 & LRPROX=>YDPHY%LRPROX, LPROCLD=>YDPHY%LPROCLD, &
 & LACDIFUS=>YDPHY%LACDIFUS, LCAPE=>YDPHY%LCAPE, LCVRAV3=>YDPHY%LCVRAV3, &
 & LHMTO=>YDPHY%LHMTO,  LVGSN=>YDPHY%LVGSN, &
 & LCVPPKF=>YDPHY%LCVPPKF, LCVPRO=>YDPHY%LCVPRO, LADJCLD=>YDPHY%LADJCLD, &
 & LGRAPRO=>YDPHY%LGRAPRO, RDECLI=>YDRIP%RDECLI, &
 & MVTS=>YDPARAR%MVTS, MALBSCA=>YDPARAR%MALBSCA, MVQS=>YDPARAR%MVQS, &
 & MRAIN=>YDPARAR%MRAIN, CMF_CLOUD=>YDPARAR%CMF_CLOUD, &
 & MALBDIR=>YDPARAR%MALBDIR, XSW_BANDS=>YDPARAR%XSW_BANDS, &
 & MVEMIS=>YDPARAR%MVEMIS, MGZ0H=>YDPARAR%MGZ0H, MGZ0=>YDPARAR%MGZ0, &
 & MCD=>YDPARAR%MCD, &
 & NSWB_MNH=>YDPARAR%NSWB_MNH, MSWDIR=>YDPARAR%MSWDIR, LMIXUV=>YDPARAR%LMIXUV, &
 & MSWDIF=>YDPARAR%MSWDIF, MSNOW=>YDPARAR%MSNOW, &
 & CMF_UPDRAFT=>YDPARAR%CMF_UPDRAFT, &
 & NTRADI=>YDTOPH%NTRADI, NTPLUI=>YDTOPH%NTPLUI, NTNEBU=>YDTOPH%NTNEBU, &
 & NTDIFU=>YDTOPH%NTDIFU, NTOZON=>YDTOPH%NTOZON, NTDRME=>YDTOPH%NTDRME, &
 & NTCVIM=>YDTOPH%NTCVIM, NTCOET=>YDTOPH%NTCOET, NTDRAG=>YDTOPH%NTDRAG, &
 & RMESOQ=>YDTOPH%RMESOQ, RMESOT=>YDTOPH%RMESOT, RMESOU=>YDTOPH%RMESOU, &
 & NTQSAT=>YDTOPH%NTQSAT, NTCOEF=>YDTOPH%NTCOEF, &
 & NAER=>YDERAD%NAER, NOZOCL=>YDERAD%NOZOCL, &
 & NRADFR=>YDERAD%NRADFR, NSW=>YDERAD%NSW, RCARDI=>YDERDI%RCARDI, &
 & RSUNDUR=>YDERDI%RSUNDUR, LXVISI=>YDXFU%LXVISI,&
 & LFLEXDIA=>YDLDDH%LFLEXDIA, LXVISI2=>YDXFU%LXVISI2,&
 & LDDH_OMP=>YDLDDH%LDDH_OMP, &
 & LMCC03=>YDMCC%LMCC03, &
 & LRCOEF=>YDRCOEF%LRCOEF, &
 & NSTOP=>YDRIP%NSTOP, &
 & RCODEC=>YDRIP%RCODEC, &
 & RHGMT=>YDRIP%RHGMT, &
 & RSIDEC=>YDRIP%RSIDEC, &
 & RSOVR=>YDRIP%RSOVR, &
 & RSTATI=>YDRIP%RSTATI, &
 & TSTEP=>YDRIP%TSTEP, &
 & LDPRECIPS=>YDPHY%LDPRECIPS,LDPRECIPS2=>YDPHY%LDPRECIPS2,&
 & NDTPREC=>YDPHY%YRDPRECIPS%NDTPREC,NDTPREC2=>YDPHY%YRDPRECIPS%NDTPREC2,&
 & NDTPRECCUR=>YDPHY%YRDPRECIPS%NDTPRECCUR,NDTPRECCUR2=>YDPHY%YRDPRECIPS%NDTPRECCUR2,&
 & STPRE=>YDSTA%STPRE,STPREH=>YDSTA%STPREH, STTEM=>YDSTA%STTEM, &
 & NORGWD_NNOVERDIF=>YDNORGWD%NORGWD_NNOVERDIF, &
 & LGCHECKMV=>YDPHY%LGCHECKMV, NAERO=>YGFL%NAERO, NGFL_EXT=>YGFL%NGFL_EXT, &
 & LCHEM_ARPCLIM=>YDMODEL%YRML_CHEM%YRCHEM%LCHEM_ARPCLIM,NCHEM=>YDMODEL%YRML_GCONF%YGFL%NCHEM, &
 & NACTAERO=>YGFL%NACTAERO, LXMRT=>YDXFU%LXMRT,RDELXN=>YDGEM%RDELXN,LFLASH =>YDCFU%LFLASH)
 CGMIXLEN=>YDMODEL%YRML_PHY_MF%YRPHY%CGMIXLEN
!
!-------------------------------------------------
! Check magnitude of model variables.
!-------------------------------------------------
!
IF(LGCHECKMV) CALL CHECKMV(YDRIP,YDPHY0,YDPHY2,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDCPG_DIM%KSTEP,YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF &
& ,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,PGELAM,PGEMU,YDMF_PHYS%TMP%RDG%MU0,YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%Q,YDMF_PHYS_STATE%YGSP_RR%T)
!     ------------------------------------------------------------------

LLREDPR=.FALSE.
ZRVMD=RV-RD
! SURFEX  and passive scalar
IF (LDXFUMSE) THEN
  ZDTMSE=0.01_JPRB
  ZSTATI=REAL(RSTATI,JPRB)-ZDTMSE/2._JPRB
ELSE
  ZDTMSE=TSPHY
  ZSTATI=REAL(RSTATI,JPRB)
ENDIF
ZRHGMT=REAL(RHGMT,JPRB)
ZAIPCMT(:)=0._JPRB

ALLOCATE(ZSVM   (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,NGFL_EXT))
ALLOCATE(ZSFSV  (YDCPG_DIM%KLON,NGFL_EXT))
ALLOCATE(ZPSV   (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,NGFL_EXT))

!     ------------------------------------------------------------------
!     1.- INITIALISATIONS COMPLEMENTAIRES
!     -----------------------------------
IJN = YDCPG_DIM%KLON

!     CALCUL FIN DE IJN
!      IJN = 1
!      ZMU0 = MAX(0.,-SIGN(1.,-PMU0(KIDIA)))
!      DO 1  JLON = KIDIA+1, KFDIA
!      ZMUN = MAX(0.,-SIGN(1.,-PMU0(JLON)))
!      IJN = IJN + IABS(NINT(ZMUN)-NINT(ZMU0))
!      ZMU0 = ZMUN
!   1  CONTINUE

!*        1.0 DECORRELATION DEPTH FOR CLOUD OVERLAPS
IF ( RDECRD <= 0._JPRB .OR. LRNUEXP ) THEN
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZDECRD(JLON)=RDECRD1+RDECRD2* &
     & EXP(-((ASIN(PGEMU(JLON))-RDECRD3*RDECLI)/RDECRD4)**2)
  ENDDO
ENDIF
IF ( RDECRD <= 0._JPRB ) THEN 
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZDECRD_MF(JLON)=ZDECRD(JLON)
  ENDDO
ELSE
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZDECRD_MF(JLON)=RDECRD
  ENDDO
ENDIF

!*        1.1 INITIALISATION DE L'OZONE


  IF (KVCLIS == 1) THEN
    ZEPSO3=1.E-11_JPRB
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZQO3(JLON,0)=1.E-9_JPRB
    ENDDO
    DO JLEV=KTDIA,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZQO3(JLON,JLEV) = PKOZO(JLON,JLEV,1)
      ENDDO
    ENDDO

  ELSEIF ((NOZOCL == 1)) THEN
    IF (MOD(YDCPG_DIM%KSTEP,NRADFR) == 0) THEN
      CALL RADOZCMF(YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,PGEMU,ZROZ)
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZQO3(JLON,0)=1.E-9_JPRB
      ENDDO
      DO JLEV=KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZQO3(JLON,JLEV)=ZROZ(JLON,JLEV)*YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

  ELSE
    CALL SUOZON(YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,ZQO3,.FALSE.,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,LO3ABC,YDMF_PHYS_SURF%GSD_VC%PGROUP)
  ENDIF

!     GAZ CARBONIQUE.

  DO JLEV=KTDIA,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZQCO2(JLON,JLEV)=QCO2
    ENDDO
  ENDDO


  ! INITIALISATION DE LA COORDONNEE ETA.
  ! INITIALISATION DE LA COORDONNEE ETA.

  ZVETAH(KTDIA-1)=STPREH(KTDIA-1)/VP00
  DO JLEV=KTDIA,YDCPG_DIM%KFLEVG
    ZVETAH(JLEV)=STPREH(JLEV)/VP00
    ZVETAF(JLEV)=STPRE (JLEV)/VP00
  ENDDO

!     EPAISSEUR STD AEROSOLS
  ZAEO = AERCS1*ZVETAH(KTDIA-1) + AERCS3*ZVETAH(KTDIA-1)**3&
   & + AERCS5*ZVETAH(KTDIA-1)**5
  DO JLEV = KTDIA, YDCPG_DIM%KFLEVG
    ZAEN = AERCS1*ZVETAH(JLEV) + AERCS3*ZVETAH(JLEV)**3&
     & + AERCS5*ZVETAH(JLEV)**5
    ZDAER(JLEV) = ZAEN - ZAEO
    ZAEO = ZAEN
  ENDDO

!     HUMIDITE CRITIQUE
  ZEPS=1.E-12_JPRB
  IF(LHUCN) THEN
    DO JLEV = KTDIA, YDCPG_DIM%KFLEVG
      ZHUC(JLEV)=1.0_JPRB-MAX( HUCOE*ZVETAF(JLEV)*(1.0_JPRB-ZVETAF(JLEV))/((&
      & 1.0_JPRB+HUTIL1*(ZVETAF(JLEV)-0.5_JPRB))*(&
      & 1.0_JPRB+HUTIL2*(ZVETAF(JLEV)-0.5_JPRB))),ZEPS)
    ENDDO
  ELSE
    DO JLEV = KTDIA, YDCPG_DIM%KFLEVG
      ZHUC(JLEV)=1.0_JPRB-MAX( HUCOE*ZVETAF(JLEV)**NPCLO1*(&
       & 1.0_JPRB-ZVETAF(JLEV))**NPCLO2*(&
       & 1.0_JPRB+SQRT(HUTIL)*(ZVETAF(JLEV)-0.5_JPRB)),ZEPS)
    ENDDO
  ENDIF

  

    DO JLEV = KTDIA, YDCPG_DIM%KFLEVG-1
      ZMAN(JLEV) = FSM_CC * TANH(FSM_DD*ZVETAH(JLEV))
      ZMAK(JLEV) = FSM_EE * ZVETAH(JLEV)**FSM_FF +&
       & FSM_GG * (1-ZVETAH(JLEV))**FSM_HH + FSM_II
    ENDDO

!     MATHEMATICAL FILTER FOR THE EDGES

    

     !LNEWSTAT

!3MT
!     INCREMENTAL CORRECTION FLUX FOR NEGAVTIVE HUMIDITY VALUES

  ZFCQVNG(:,:)=0.0_JPRB
  ZFCQING(:,:)=0.0_JPRB
  ZFCQLNG(:,:)=0.0_JPRB

  ZPRODTH_CVPP(:,:)=0.0_JPRB

  ZGDT=RG*TSPHY
  ZGDTI=1.0_JPRB/ZGDT

  DO JLEV=KTDIA,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA

!     LAYER WEIGHTS

      ZPOID(JLON,JLEV)=YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZGDTI
      ZIPOI(JLON,JLEV)=YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP(JLON,JLEV)*ZGDT

!     CALCULATION OF LATENT HEATS

      ZLHV(JLON,JLEV)=FOLH(YDMF_PHYS_STATE%T(JLON,JLEV),0.0_JPRB)
      ZLHS(JLON,JLEV)=FOLH(YDMF_PHYS_STATE%T(JLON,JLEV),1.0_JPRB)

    ENDDO
  ENDDO

!    ------------------------------------------------------------------
!     PROGNOSTIC GEMS/MACC AEROSOLS - INITIAL COMPUTATIONS
!     IMPORTANT for IFS: Tracer order is : CO2 - other tracers - react Gases - Aerosol - extra GFL
!    ------------------------------------------------------------------

  ! Preliminary for prog. aerosol or extra gfl
  INBTRA=0
  
  IF (NAERO>0)                  INBTRA=NAERO           ! the two cases exclude each other
  ALLOCATE(ZSTRCTRA(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,INBTRA)) ! to cover both prog aero and extra gfl cases
  ALLOCATE(ZTRA    (YDCPG_DIM%KLON,  YDCPG_DIM%KFLEVG,INBTRA))
  IF (INBTRA > 0) THEN
    ZSTRCTRA(:,:,:) = 0._JPRB
    ZTRA(:,:,:)     = 0._JPRB
  ENDIF 
  IF(INBTRA == 0) THEN
    INBTRA_DEP=0
  ELSE
    INBTRA_DEP=1
  ENDIF

!    ------------------------------------------------------------------
!      - CHANGEMENTS DE VARIABLES ET INVERSION DES NIVEAUX
!          POUR LE TRAITEMENT DES SCALAIRES PASSIFS
!     --------------------------------------------------------------------

  ZSFSV=0.0_JPRB    ! surf. flux of scalars
    ! ENDIF  (LMDUST & NGFL_EXT)
  !LMPHYS
!**
!     ------------------------------------------------------------------
!     2.- MISES A ZERO DE SECURITE EN CAS DE NON-APPEL DES PARAMETRIS.
!     ------------------------------------------------------------------
ZEPS0=1.E-12_JPRB
ZEPSNEB=1.E-10_JPRB

! To profitize from the vectorization collapsing the (:,:) form is preferable.
! (Even better would be to completely avoid any useless initialization.)

! arrays dimensioned from 0:KLEV (half level quantities)
ZFPCOR  (:,:) = 0.0_JPRB
ZFHP    (:,:) = 0.0_JPRB
ZXTROV  (:,:) = 1.0_JPRB
ZXUROV  (:,:) = 1.0_JPRB
ZLMT    (:,:) = 0.0_JPRB
ZZLMT   (:,:) = 0.0_JPRB
ZLMU    (:,:) = 0.0_JPRB
ZLMU2   (:,:) = 0.0_JPRB
ZLMT2   (:,:) = 0.0_JPRB
ZKTROV  (:,:) = 0.0_JPRB
ZKQROV  (:,:) = 0.0_JPRB
ZKQLROV (:,:) = 0.0_JPRB
ZKUROV  (:,:) = 0.0_JPRB
ZFHEVPPC(:,:) = 0.0_JPRB
ZFHMLTSC(:,:) = 0.0_JPRB
ZFPEVPPC(:,:) = 0.0_JPRB
ZFCQL   (:,:) = 0.0_JPRB
ZFCQI   (:,:) = 0.0_JPRB
ZDIFCVPPQ (:,:) = 0.0_JPRB
ZDIFCVPPS (:,:) = 0.0_JPRB
ZDIFCVTH (:,:) = 0.0_JPRB
ZDIFCVPPU (:,:) = 0.0_JPRB
ZDIFCVPPV (:,:) = 0.0_JPRB
ZCONDCVPPL(:,:) = 0.0_JPRB
ZCONDCVPPI(:,:) = 0.0_JPRB
ZSEDIQL(:,:) = 0.0_JPRB
ZSEDIQI(:,:) = 0.0_JPRB

ZXURO   (:,:) = 0.0_JPRB
ZXQRO   (:,:) = 0.0_JPRB
ZXTRO   (:,:) = 0.0_JPRB

ZALPHA1 (:,:) = 0.0_JPRB
ZCOEFA  (:,:) = 0.0_JPRB
ZLVT    (:,:) = 0.0_JPRB
ZQICE   (:,:) = 0.0_JPRB

ZF_EPS (:,:) = 1.0_JPRB
ZFUN_TTE (:,:) = 1.0_JPRB
ZMRIPP (:,:) = 1.E-12_JPRB
ZMRIMC  (:,:) = 1.0_JPRB
ZMRICTERM (:,:) = 1.0_JPRB
ZRRCOR (:,:) = 1.0_JPRB
ZTAU_TKE (:,:) = 0.0_JPRB
ZTH_FUN (:,:) = 1.0_JPRB
ZMRIFPP (:,:) = 1.E-12_JPRB
ZMN2PP (:,:) = 1.E-12_JPRB
ZMN2_ES (:,:) = 1.0_JPRB
ZMN2_EQ (:,:) = 1.0_JPRB
ZMN2_DS (:,:) = 1.0_JPRB
ZMN2_DQ (:,:) = 1.0_JPRB
ZTSTAR (:,:) = 1.E-12_JPRB
ZTSTAR2 (:,:) = 1.E-12_JPRB
ZTSTARQ (:,:) = 1.E-12_JPRB
ZTSTAR2Q (:,:) = 1.E-12_JPRB
ZFMGST (:,:) = 1.0_JPRB
ZFMTKE (:,:) = 1.0_JPRB
ZFTTKE (:,:) = 1.0_JPRB
ZAUTKE (:,:) = 1.0_JPRB
ZATTKE (:,:) = 1.0_JPRB
ZTH_FUN(:,:) = 1.0_JPRB
ZWW_FUN(:,:) = 1.0_JPRB
ZBNEBCVPP(:,:) = 0.0_JPRB
ZBNEBQ(:,:)   = 0.0_JPRB
ZRHS(:,:)     = 0.0_JPRB
ZLML(:,:)     = 1.0_JPRB
ZLMLTILD(:,:) = 1.0_JPRB

ZDIFWQ  (:) = 0.0_JPRB
ZDIFWS  (:) = 0.0_JPRB
ZSC_FEVI(:) = 1.0_JPRB     
ZSC_FEVN(:) = 1.0_JPRB     
ZSC_FCLL(:) = 1.0_JPRB     
ZSC_FCLN(:) = 1.0_JPRB
ZCDNH(:)    = 1.0_JPRB

! arrays dimensioned from 1:KLEV (full level quantities)
ZTENT   (:,:) = 0.0_JPRB
ZNEBS   (:,:) = ZEPS0
ZNEBC   (:,:) = ZEPS0
ZNEBS0  (:,:) = ZEPS0
ZNEBC0  (:,:) = ZEPS0
ZNEBCH  (:,:) = 0.0_JPRB
ZUNEBH  (:,:) = 0.0_JPRB
ZDETFI (:,:) = 0.0_JPRB
ZNEBDIFF(:,:) = 0.0_JPRB
ZQLIS   (:,:) = 0.0_JPRB
ZQLIS0  (:,:) = 0.0_JPRB
ZCFATH  (:,:) = 0.0_JPRB
ZCFAU   (:,:) = 0.0_JPRB
ZCFBTH  (:,:) = 0.0_JPRB
ZCFBU   (:,:) = 0.0_JPRB
ZCFBV   (:,:) = 0.0_JPRB
ZQLIC   (:,:) = 0.0_JPRB
INLAB   (:,:) = 0
INLAB_CVPP(:,:) = 0
ICIS    (:,:) = 1
ZQLI_CVPP(:,:) = 0.0_JPRB
ZNEB_CVPP(:,:) = ZEPS0

ZEDMFQ  (:,:)  = 0.0_JPRB
ZEDMFS  (:,:)  = 0.0_JPRB
ZEDMFU  (:,:)  = 0.0_JPRB
ZEDMFV  (:,:)  = 0.0_JPRB
ZMF_UP  (:,: ) = 0.0_JPRB
ZQLI_CVP(:,:) = 0.0_JPRB
ZQC_DET_PCMT(:,:) = 0.0_JPRB
ZTENHA(:,:)   = 0.0_JPRB
ZTENQVA(:,:)  = 0.0_JPRB
ZRHDFDA(:,:)  = 0.0_JPRB
ZQIC   (:,:)  = 0.0_JPRB
ZQLC   (:,:)  = 0.0_JPRB
ZQRC   (:,:)  = 0.0_JPRB
ZQSC   (:,:)  = 0.0_JPRB
ZQG    (:,:)  = 0.0_JPRB

!  ---------------------------------------------------------------------
!  Correction of negative advected humidity and precipitation values
!  ---------------------------------------------------------------------


  
    
      DO JLEV = 1, YDCPG_DIM%KFLEVG
        DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZQI(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_STATE%I(JLON,JLEV))
          ZQL(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_STATE%L(JLON,JLEV))
          ZQR(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_STATE%R(JLON,JLEV))
          ZQS(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_STATE%S(JLON,JLEV))

    ! CORRECTION OF NEGATIVE ADVECTED VALUES:
    !                         VAPOUR PUT IN PFCQVNG
    !                         LIQUID,ICE PUT IN PFCQL/ING
    !                         FOR RAIN/SNOW PUT IN PFCQR/SNG

          ZDQI=ZQI(JLON,JLEV)-YDMF_PHYS_STATE%I(JLON,JLEV)
          ZDQL=ZQL(JLON,JLEV)-YDMF_PHYS_STATE%L(JLON,JLEV)
          ZDQR=ZQR(JLON,JLEV)-YDMF_PHYS_STATE%R(JLON,JLEV)
          ZDQS=ZQS(JLON,JLEV)-YDMF_PHYS_STATE%S(JLON,JLEV)
          ZDQC=ZDQI+ZDQL+ZDQR+ZDQS

          ZQV0=YDMF_PHYS_STATE%Q(JLON,JLEV)-ZIPOI(JLON,JLEV)*(0.0_JPRB- YDMF_PHYS%OUT%FCQNG(JLON,JLEV-1)&
          & -YDMF_PHYS%OUT%FCQNNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQLNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQRNG(JLON,JLEV-1)&
          & -YDMF_PHYS%OUT%FCQSNG(JLON,JLEV-1))
          ZQV(JLON,JLEV)=MAX(0.0_JPRB,ZQV0-ZDQC)
          ZDQV=MAX(0.0_JPRB,ZQV0-ZDQC)-YDMF_PHYS_STATE%Q(JLON,JLEV)

          YDMF_PHYS%OUT%FCQNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNG(JLON,JLEV-1)-ZDQV*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQNNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNNG(JLON,JLEV-1)-ZDQI*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQLNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQLNG(JLON,JLEV-1)-ZDQL*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQRNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQRNG(JLON,JLEV-1)-ZDQR*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQSNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQSNG(JLON,JLEV-1)-ZDQS*ZPOID(JLON,JLEV)
        ENDDO
      ENDDO
    

    
  
 ! LCONDWT

DO JCHA = 1, 6
  DO JLEV = KTDIA, YDCPG_DIM%KFLEVG
    DO JLON = 1, YDCPG_DIM%KLON
      ZAER(JLON,JLEV,JCHA)=0.0_JPRB
    ENDDO
  ENDDO
ENDDO

DO JLEV = KTDIA, YDCPG_DIM%KFLEVG
  DO JLON = 1, YDCPG_DIM%KLON
    ZAERINDS(JLON,JLEV)=0.0_JPRB
  ENDDO
ENDDO

DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
  ZCEMTR  (JLON,0) = 0.0_JPRB
  ZCEMTR  (JLON,1) = 0.0_JPRB
  ZCTRSO  (JLON,0) = 0.0_JPRB
  ZCTRSO  (JLON,1) = 0.0_JPRB
  ZTRSOD  (JLON)   = 0.0_JPRB
  ZSUDU   (JLON)   = 0.0_JPRB
  ZXDROV  (JLON)   = 1.0_JPRB
  ZXHROV  (JLON)   = 1.0_JPRB
  ZTAUX   (JLON)   = ZEPS0
  KCLPH   (JLON)   = YDCPG_DIM%KFLEVG
ENDDO
DO JSG = 1, NSW
  DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZALBD   (JLON,JSG) = 0.0_JPRB
    ZALBP   (JLON,JSG) = 0.0_JPRB
    ZSFSWDIF(JLON,JSG) = 0.0_JPRB
    ZSFSWDIR(JLON,JSG) = 0.0_JPRB
    ZTRSODIF(JLON,JSG) = 0.0_JPRB
    ZTRSODIR(JLON,JSG) = 0.0_JPRB
  ENDDO
ENDDO

!  -------------------------------------------------------
!  Security values for pseudo-historical arrays at KSTEP=0
!  -------------------------------------------------------




IF((LNEBN.OR.LRRGUST).AND.YDCPG_DIM%KSTEP == 0) THEN
  DO JLEV=KTDIA-1,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%TMP%APLPAR%PFL%FPLCH(JLON,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF
IF(LRRGUST.AND.YDCPG_DIM%KSTEP == 0) THEN
  DO JLEV=KTDIA-1,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%TMP%APLPAR%PFL%FPLSH(JLON,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF




IF(LRKCDEV.AND.YDCPG_DIM%KSTEP == 0) THEN
  DO JLEV=KTDIA,YDCPG_DIM%KFLEVG
     DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZDTRAD(JLON,JLEV)=0.0_JPRB
        ZRKTTEND(JLON,JLEV)=0.0_JPRB
        ZRKQVTEND(JLON,JLEV)=0.0_JPRB
        ZRKQCTEND(JLON,JLEV)=0.0_JPRB
        ZDQVDIFF(JLON,JLEV)=0.0_JPRB
        YDVARS%RKTH%T0(JLON,JLEV) = 0.0_JPRB
        YDVARS%RKTQV%T0(JLON,JLEV)= 0.0_JPRB
        YDVARS%RKTQC%T0(JLON,JLEV)= 0.0_JPRB
     ENDDO
  ENDDO
ENDIF


!     ------------------------------------------------------------------
!  The LMPHYS and LEPHYS keys should in the following be
! dispached to individual parametrizations.

!*
!     ------------------------------------------------------------------
!     4.- CALCULS THERMODYNAMIQUES
!     ----------------------------
  
    CALL ACTQSAT ( YDPHY,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,NTQSAT,YDCPG_DIM%KFLEVG,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, ZQV, YDMF_PHYS_STATE%T,&
     & ZGEOSLC, YDMF_PHYS%TMP%APLPAR%MSC%LH, YDMF_PHYS%TMP%APLPAR%MSC%LSCPE, YDMF_PHYS%TMP%APLPAR%FLU%QSAT, YDMF_PHYS%TMP%APLPAR%MSC%QW, YDCPG_MISC%RH, YDMF_PHYS%TMP%APLPAR%MSC%TW)
  

!*
!     ------------------------------------------------------------------
!     4.BIS. COEFFICIENTS THERMO-HYDRIQUES DU SOL
!     -------------------------------------------

    ! End of surface forcing for 1D model MUSC

  IF ( .NOT.LMSE ) THEN
    IF ( LSOLV ) THEN
      LLHMT=.FALSE.
      CALL ACSOL ( YDPHY,YDPHY1,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,&
       & YDMF_PHYS_SURF%GSD_VV%PARG,YDMF_PHYS_SURF%GSD_VV%PD2,YDMF_PHYS_SURF%GSD_VF%PZ0F,YDMF_PHYS_SURF%GSD_VV%PZ0H,YDMF_PHYS_SURF%GSD_VF%PZ0RLF,YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS_SURF%GSD_VV%PIVEG,YDMF_PHYS_SURF%GSD_VV%PLAI,YDMF_PHYS_STATE%YGSP_SG%A,&
       & YDMF_PHYS_STATE%YGSP_SG%R,YDMF_PHYS_SURF%GSD_VV%PSAB,YDMF_PHYS_STATE%YGSP_SG%F,YDMF_PHYS_STATE%YGSP_RR%T,YDMF_PHYS_SURF%GSD_VF%PVEG,YDMF_PHYS_STATE%YGSP_SB%Q,YDMF_PHYS_STATE%YGSP_SB%TL,YDMF_PHYS_STATE%YGSP_RR%W,YDMF_PHYS_STATE%YGSP_RR%IC,&
       & LLHMT,&
       & YDMF_PHYS%TMP%APLPAR%DSA%C1,YDMF_PHYS%TMP%APLPAR%DSA%C2,ZC3,ZCG,ZCN,YDMF_PHYS%OUT%CT,ZNEIJG,ZNEIJV,&
       & ZWFC,ZWPMX,ZWSEQ,ZWSMX,ZWWILT)
    ELSE

!            INITIALISATION DE L'INERTIE THERMIQUE DU SOL.

!DEC$ IVDEP
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        YDMF_PHYS%OUT%CT(JLON)=HSOL /(&
         & 1.0_JPRB+HSOLIWR*(YDMF_PHYS_STATE%YGSP_RR%W(JLON)+YDMF_PHYS_STATE%YGSP_SB%Q(JLON,1))/(WSMX+WPMX)&
         & *EXP(-0.5_JPRB*(HSOLIT0*(YDMF_PHYS_STATE%YGSP_RR%T(JLON)-RTT))**2))
      ENDDO
    ENDIF
  ENDIF

!*
!     ------------------------------------------------------------------
!     4.TER.- INITIALISATIONS LIEES AU SCHEMA DE SURFACE EXTERNALISE
!     ------------------------------------------------------------------

  IF (LMSE) THEN

!     INITIALISATION DU SCHEMA DE SURFACE EXTERNALISE ET DES
!     VARIABLES PSEUDO-HISTORIQUES ASSOCIEES

    IF ( (NSWB_MNH /= NSW) ) CALL ABOR1('APLPAR: NSWB_MNH not = NSW')

    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%GZ0     (JLON) = PGPAR(JLON,MGZ0)
      YDMF_PHYS%OUT%GZ0H    (JLON) = PGPAR(JLON,MGZ0H)
      YDMF_PHYS%TMP%APLPAR%FLU%EMIS    (JLON) = PGPAR(JLON,MVEMIS)
      YDCPG_MISC%QS      (JLON) = PGPAR(JLON,MVQS)
      YDMF_PHYS_STATE%YGSP_RR%T      (JLON) = PGPAR(JLON,MVTS)
      ZSRAIN   (JLON) = PGPAR(JLON,MRAIN)
      ZSSNOW   (JLON) = PGPAR(JLON,MSNOW)
      ZSGROUPEL(JLON) = 0._JPRB
      ZTSN     (JLON) = YDMF_PHYS_STATE%YGSP_RR%T(JLON)
    ENDDO

    
      ! FMR  radiation => NSW solar bands
      DO JSG=1,NSW
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZALBP(JLON,JSG) = PGPAR(JLON,MALBDIR-1+JSG)
          ZALBD(JLON,JSG) = PGPAR(JLON,MALBSCA-1+JSG)
        ENDDO
      ENDDO
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        YDMF_PHYS%OUT%ALB(JLON)=0.0_JPRB
      ENDDO
      DO JSG=1,NSW
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          YDMF_PHYS%OUT%ALB(JLON)=YDMF_PHYS%OUT%ALB(JLON)+0.5*(ZALBP(JLON,JSG)+ZALBD(JLON,JSG))
        ENDDO
      ENDDO
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        YDMF_PHYS%OUT%ALB(JLON)=YDMF_PHYS%OUT%ALB(JLON)/REAL(NSW,JPRB)
      ENDDO
    

    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      IF (YDMF_PHYS%TMP%APLPAR%FLU%EMIS(JLON)==0._JPRB) THEN
        YDMF_PHYS%TMP%APLPAR%FLU%EMIS(JLON) = 0.99_JPRB
        YDMF_PHYS_STATE%YGSP_RR%T  (JLON) = 288.0_JPRB
        YDMF_PHYS%OUT%ALB (JLON) = 0.1_JPRB
      ENDIF
    ENDDO

  ENDIF  ! LMSE
!
! Define z0;z0h if it's necessary
!
  IF (LMSE.AND.YDCPG_DIM%KSTEP == 0) THEN
  CALL ARO_GROUND_DIAG_Z0( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP,&
          & YDCPG_DIM%KFDIA-YDCPG_DIM%KIDIA+1, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA,&
          & NDGUNG, NDGUXG, NDLUNG, NDLUXG,&
          & PINDX(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), PINDY(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & LSURFEX_KFROM,&
          & YDMF_PHYS%OUT%GZ0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%GZ0H(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA))
  ENDIF
!*
!     ------------------------------------------------------------------
!     5.- STRUCTURE ET CHAMPS DANS LA COUCHE LIMITE DE SURFACE
!     ------------------------------------------------------------------

!        INITIALISATION DES HAUTEURS "METEO".

  
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZDPHIV(JLON)=RG*HVCLS
      ZDPHIT(JLON)=RG*HTCLS
    ENDDO
  

  
    LLCLS=.TRUE.
    LLHMT=.TRUE.
    IF (LMSE) THEN
      
        CALL ACHMTLS ( YDMODEL%YRML_PHY_MF,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,&
          & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,&
          & YDMF_PHYS_STATE%YGSP_RR%T,YDCPG_MISC%QS,ZDPHIT,YDMF_PHYS%OUT%GZ0,YDMF_PHYS%OUT%GZ0H,LLCLS,&
          & ZNBVNO,ZMRIPP,YDMF_PHYS%TMP%APLPAR%DSA%CPS,ZGWDCS,YDMF_PHYS%TMP%APLPAR%DSA%LHS,ZPCLS,YDMF_PHYS%TMP%APLPAR%FLU%CD,YDMF_PHYS%TMP%APLPAR%FLU%CDN)
!       Computation of ZRTI
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZDPHI(JLON)=YDMF_PHYS_STATE%YCPG_DYN%PHIF(JLON,YDCPG_DIM%KFLEVG)-YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,YDCPG_DIM%KFLEVG)
          ZPRS(JLON)=RD+ZRVMD*YDCPG_MISC%QS(JLON)
          ZRTI(JLON)=2.0_JPRB/(YDMF_PHYS_STATE%YCPG_DYN%RCP%R(JLON,YDCPG_DIM%KFLEVG)*YDMF_PHYS_STATE%T(JLON,YDCPG_DIM%KFLEVG)+RKAPPA*ZDPHI(JLON)&
           & +ZPRS(JLON)*YDMF_PHYS_STATE%YGSP_RR%T(JLON))
        ENDDO
      
    ELSE
      
        CALL ACHMT ( YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,KVCLIV>=8.AND.LSOLV,&
         & YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, ZQV, YDMF_PHYS_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_STATE%T, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V,&
         & YDMF_PHYS%TMP%APLPAR%PFL%FPLSH,YDMF_PHYS%TMP%APLPAR%PFL%FPLCH,&
         & ZDPHIT, ZDPHIV, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H, YDMF_PHYS_SURF%GSD_VF%PZ0RLF, YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VF%PLSM,&
         & ZNEIJG, ZNEIJV, YDMF_PHYS_STATE%YGSP_SG%F, YDMF_PHYS_STATE%YGSP_RR%T, YDMF_PHYS_SURF%GSD_VF%PVEG, ZWFC, YDMF_PHYS_STATE%YGSP_RR%W, YDMF_PHYS_STATE%YGSP_RR%IC,&
         & LLCLS, LLHMT,&
         & ZNBVNO,ZMRIPP,&
         & YDMF_PHYS%TMP%APLPAR%FLU%CD, YDMF_PHYS%TMP%APLPAR%FLU%CDN, ZCDROV, YDMF_PHYS%TMP%APLPAR%FLU%CH, ZCHROV, YDMF_PHYS%TMP%APLPAR%DSA%CPS, ZDQSTS, ZGWDCS,&
         & YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, ZHQ, ZHU, YDMF_PHYS%TMP%APLPAR%FLU%NEIJ, YDMF_PHYS%OUT%QCLS, YDCPG_MISC%QS, YDMF_PHYS%TMP%APLPAR%FLU%QSATS, YDMF_PHYS%OUT%RHCLS,&
         & YDMF_PHYS%TMP%APLPAR%DSA%RS, ZRTI, ZSTAB, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS%OUT%NUCLS, YDMF_PHYS%OUT%NVCLS, ZPCLS, YDMF_PHYS%TMP%APLPAR%FLU%VEG,&
         & ZXDROV, ZXHROV,YDMF_PHYS%OUT%UGST,YDMF_PHYS%OUT%VGST)
      
    ENDIF

    
    

    
      ZCP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG) = YDMF_PHYS_STATE%YCPG_DYN%RCP%CP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
    

    

    ! COMPUTATION OF mixing lengths from  Ri*,Ri** - FIRST GUES for moist AF

    !---------------------------------------------------
    ! COMPUTATION OF 'DRY' mixing lengths : lm_d lh_d
    ! COMPUTATION OF ZPBLH - PBL HEIGHT

    
      DO JLEV=KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZTHETAV(JLON,JLEV)=YDMF_PHYS_STATE%T(JLON,JLEV)*(1.0_JPRB+RETV*ZQV(JLON,JLEV))&
           & *(RATM/YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(JLON,JLEV))**RKAPPA  
        ENDDO
      ENDDO
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZTHETAVS(JLON)=YDMF_PHYS_STATE%YGSP_RR%T(JLON)*(1.0_JPRB+RETV*YDCPG_MISC%QS(JLON))&
         & *(RATM/YDMF_PHYS_STATE%YCPG_PHY%PREHYD(JLON,YDCPG_DIM%KFLEVG))**RKAPPA  
      ENDDO
      CALL ACCLPH ( YDPHY0,YDPHY2,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,KTDIA,YDCPG_DIM%KFLEVG,&
       & ZTHETAV,YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,ZTHETAVS,&
       & KCLPH,YDMF_PHYS%OUT%CLPH,YDMF_PHYS%OUT%VEIN,ZUGST,ZVGST)
      ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS%OUT%CLPH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      
    

   ! end of LVDIF or LHMTO or LGWD

   ! (LVDIF or LGWD) and( not(LNEBR or LECT))

  

    IF ((.NOT.LMSE)) THEN
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZCRTI(JLON) = 1.0_JPRB/(YDMF_PHYS_STATE%YGSP_RR%T(JLON)*YDMF_PHYS%TMP%APLPAR%DSA%RS(JLON))
      ENDDO
      CALL ACFLUSO ( YDMODEL%YRML_PHY_MF,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%Q,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,&
       & ZDPHIT,ZDPHIV,YDMF_PHYS%OUT%GZ0,YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS%TMP%APLPAR%FLU%QSATS,ZCRTI,YDMF_PHYS_STATE%YGSP_RR%T,&
       & LLHMT,YDMF_PHYS%TMP%APLPAR%FLU%CD,YDMF_PHYS%TMP%APLPAR%FLU%CDN,ZCDROV,ZCE,ZCEROV,YDMF_PHYS%TMP%APLPAR%FLU%CH,ZCHROV,&
       & YDMF_PHYS%OUT%QCLS,YDMF_PHYS%OUT%RHCLS,YDMF_PHYS%OUT%TCLS,YDMF_PHYS%OUT%UCLS,YDMF_PHYS%OUT%VCLS,YDMF_PHYS%OUT%UGST,YDMF_PHYS%OUT%VGST)
    ELSE
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZCE   (JLON) = YDMF_PHYS%TMP%APLPAR%FLU%CH   (JLON)
        ZCEROV(JLON) = ZCHROV(JLON)
      ENDDO
    ENDIF

   ! LVDIF or LHMTO or LGWD

  IF (LRAFTKE) THEN
    YDMF_PHYS%OUT%CAPE(:)=0._JPRB
    ZDCAPE(:)=0._JPRB
    CALL ACCLDIA(YDXFU,YDPHY,YDPHY2,YDTOPH,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,&
      & YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, YDMF_PHYS%OUT%CAPE, ZDCAPE, YDMF_PHYS_STATE%TKE, YDMF_PHYS_STATE%YCPG_DYN%PHIF, POROG, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST,ZBLH,KCLPH)
  ENDIF

!**
!     ------------------------------------------------------------------
!     6.- TURBULENCE: COEFFICIENTS D'ECHANGE
!     ------------------------------------------------------------------

  

  
    CALL ACVPPKF( YDMODEL%YRML_PHY_MF,YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTCVIM, YDCPG_DIM%KFLEVG,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_STATE%T, ZQV,&
     & ZQL, ZQI, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, YDCPG_DYN0%CTY%VVEL(:,1:), YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_STATE%TKE,&
     & ZDIFCVPPQ, ZDIFCVPPS, ZCONDCVPPL, ZCONDCVPPI,&
     & ZPRODTH_CVPP, INLAB_CVPP, ZQLI_CVPP, ZNEB_CVPP, INND)
  

     !-------------------------------------------------
     !  Call to EDKF
     !-------------------------------------------------
  

  
    
       YDCPG_MISC%QICE(:,:)= ZQI(:,:)
       YDCPG_MISC%QLI(:,:) = ZQL(:,:)
    

! Computation of the 2 znlab used in acbl89
    
    
       ZNLABCVP(:,:) = 1.0_JPRB
    
    IF(LNEBN.OR.LRRGUST) THEN
      DO JLEV=KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZNLABCVP(JLON,JLEV) = ZNLABCVP(JLON,JLEV)&
         & *MAX(0.0_JPRB,SIGN(1.0_JPRB,YDMF_PHYS%TMP%APLPAR%PFL%FPLCH(JLON,JLEV)-YDMF_PHYS%TMP%APLPAR%PFL%FPLCH(JLON,JLEV-1)-ZEPS))
          ZNLAB(JLON,JLEV) = REAL(INLAB_CVPP(JLON,JLEV),JPRB)
        ENDDO
      ENDDO
    ENDIF
    DO JLEV=KTDIA,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZNLAB(JLON,JLEV) = REAL(INLAB_CVPP(JLON,JLEV),JPRB)
      ENDDO
    ENDDO

    CALL ACTKE ( YDLDDH,YDMODEL%YRML_DIAG%YRMDDH,YDMODEL%YRML_PHY_MF,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,NTCOEF,NTCOET,YDCPG_DIM%KFLEVG,&
     & YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,&
     & YDMF_PHYS_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_STATE%T, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, ZQV, ZQIC, ZQLC, YDMF_PHYS%TMP%APLPAR%MSC%LSCPE,&
     & YDMF_PHYS%TMP%APLPAR%FLU%CD, YDMF_PHYS%TMP%APLPAR%FLU%CH, YDMF_PHYS%OUT%GZ0, YDMF_PHYS_STATE%YGSP_RR%T, YDCPG_MISC%QS,&
     & YDCPG_MISC%QICE, YDCPG_MISC%QLI, YDMF_PHYS_STATE%TKE, ZPRODTH_CVPP, ZNLAB, ZNLABCVP,&
     & ZKTROV, ZKQROV, ZKQLROV, ZKUROV, ZXTROV, ZXUROV,&
     & ZNBVNO, ZNEBS, ZQLIS, ZNEBS0, ZQLIS0, ZCOEFN , YDMF_PHYS%TMP%APLPAR%PFL%FTKE, YDMF_PHYS%TMP%APLPAR%PFL%FTKEI, ZTKE1,ZTPRDY,&
     & PTDISS,YDDDH)
    YDMF_PHYS%OUT%CLPH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=MIN(XMAXLM,MAX(XMINLM,YDMF_PHYS%OUT%CLPH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)))
  


     !-------------------------------------------------
     ! Store diffusion coefficients in trajectory in temporary variables
     ! before final writing.
     !-------------------------------------------------
  ZKTROV_SAVE(:,:)=ZKTROV(:,:)
  ZKUROV_SAVE(:,:)=ZKUROV(:,:)
  ZCDROV_SAVE(:)=ZCDROV(:)
  ZCHROV_SAVE(:)=ZCHROV(:)

!**
!     ------------------------------------------------------------------
!     7.- RAYONNEMENT
!     ----------------
!     --------------------------------------------------------------------
!      - COMPUTE DUST PROPERTIES FOR RADIATION IF LMDUST=T
!     --------------------------------------------------------------------
  

!      7.1 Albedo et emissivite en presence de neige
!          Albedo and emissivity with snow

  IF (.NOT.LMSE) THEN
!DEC$ IVDEP
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      
        
          YDMF_PHYS%OUT%ALB(JLON)=YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)-YDMF_PHYS%TMP%APLPAR%FLU%NEIJ(JLON)*(YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)&
           & -MAX(YDMF_PHYS_SURF%GSD_VF%PALBF(JLON),ALCRIN))
          YDMF_PHYS%TMP%APLPAR%FLU%EMIS(JLON)=YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)-YDMF_PHYS%TMP%APLPAR%FLU%NEIJ(JLON)*(YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)-EMCRIN)
        
      
    ENDDO

    
      ! diffuse and direct (parallel) albedo in NSW solar intervals
      
!DEC$ IVDEP
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZMERL(JLON)=(1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PLSM(JLON))*(1.0_JPRB&
            & -MAX(0.0_JPRB,SIGN(1.0_JPRB,TMERGL-YDMF_PHYS_STATE%YGSP_RR%T(JLON))))
          YDMF_PHYS%TMP%APLPAR%FLU%EMIS(JLON)=YDMF_PHYS%TMP%APLPAR%FLU%EMIS(JLON)*YDMF_PHYS_SURF%GSD_VF%PLSM(JLON)+ZMERL(JLON)*EMMMER&
            & +(1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PLSM(JLON))*(1.0_JPRB-ZMERL(JLON))*EMMGLA
        ENDDO
        DO JSG=1,NSW
          DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            ZALBD(JLON,JSG)=(1.0_JPRB-ZMERL(JLON))*YDMF_PHYS%OUT%ALB(JLON)+&
                                    & ZMERL(JLON) *ALBMED
            ZALBP(JLON,JSG)=(1.0_JPRB-ZMERL(JLON))*YDMF_PHYS%OUT%ALB(JLON)+&
             & ZMERL(JLON) *&
             & MAX(0.037_JPRB/(1.1_JPRB*YDMF_PHYS%TMP%RDG%MU0(JLON)**1.4_JPRB+0.15_JPRB),ZEPS0)
          ENDDO
        ENDDO
      
    

  ENDIF  ! .NOT.LMSE

  ! Appel de la routine d'aerosols

  LLAERO=LAEROSEA.AND.LAEROLAN.AND.LAEROSOO.AND.LAERODES

  

  IF    (   ((MOD(YDCPG_DIM%KSTEP,NRADFR) == 0)) ) THEN

    IF (LLAERO) THEN
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZAESEA(JLON) = YDMF_PHYS_SURF%GSD_VA%PSEA(JLON)
        ZAELAN(JLON) = YDMF_PHYS_SURF%GSD_VA%PLAN(JLON)
        ZAESOO(JLON) = YDMF_PHYS_SURF%GSD_VA%PSOO(JLON)
        ZAEDES(JLON) = YDMF_PHYS_SURF%GSD_VA%PDES(JLON)
      ENDDO
    ELSE
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZAESEA(JLON) = 0.0_JPRB
        ZAELAN(JLON) = 0.0_JPRB
        ZAESOO(JLON) = 0.0_JPRB
        ZAEDES(JLON) = 0.0_JPRB
      ENDDO
    ENDIF
    
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZAESUL(JLON) = 0.0_JPRB
      ENDDO
    
    
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZAEVOL(JLON) = 0.0_JPRB
      ENDDO
    

    IF ( ( (NAER /= 0)).AND.LLAERO )  THEN
      CALL RADAER ( YDMODEL%YRML_PHY_RAD%YREAERD,YDERAD,YDPHY, YDCPG_DIM%KIDIA , YDCPG_DIM%KFDIA , YDCPG_DIM%KLON  , YDCPG_DIM%KFLEVG,&
       & YDMF_PHYS_STATE%YCPG_PHY%PREHYD , YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%T    , YDMF_PHYS_STATE%YGSP_RR%T,&
       & ZAESEA, ZAELAN, ZAESOO, ZAEDES, ZAESUL, ZAEVOL,&
       & ZAER,ZAERINDS                                 )
    ENDIF

  ENDIF ! FOR AEROSOLS

! We uses the extinction coefficient explicitely solved by ARO_MNHDUST
  


!      7.2 Flux radiatifs par ciel clair (Code Geleyn)
!          Clear sky radiative fluxes    (Geleyn's scheme)

  ! separate clearsky call is kept only for old ACRANEB; for ACRANEB2
  ! duplicit calculation of gaseous transmissions is avoided
  

!      7.3 Nebulosite et Convection
!          Cloud cover and Convection
!      7.3.1 Shallow + Deep convection

   !  (LCVPGY)

  

    

    CALL ACNEBCOND ( YDRIP,YDMODEL%YRML_PHY_MF,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,NTPLUI,YDCPG_DIM%KFLEVG,LLREDPR,&
     & ZHUC,ZVETAF,YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDCPG_MISC%RH,ZBLH,ZQV,ZQI,ZQL,&
     & YDMF_PHYS%TMP%APLPAR%MSC%QW,YDMF_PHYS_STATE%T,ZNEBCH,PGM,YDMF_PHYS_STATE%YGSP_RR%T,&
     & ZQLIS,ZNEBS,ZRHCRI,ZRH,ZQSATS,ZICEFR1,ZQLIS0,ZNEBS0)
     
     
    
    

    IF(LRKCDEV) THEN
! Rash-Kristiansson cloud water scheme - second part.
      DO JLEV=KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ! analytical solution dRH/dCLOUD:                                                
          ZRHDFDA(JLON,JLEV)=2._JPRB*(1._JPRB - ZNEBS(JLON,JLEV))&
               & *(1.0_JPRB-ZRHCRI(JLON,JLEV))
        ENDDO
      ENDDO
    ENDIF


   ! LCONDWT .AND. .NOT.LNEBECT

  !-------------------------------------------------
  ! PCMT convection scheme.
  !-------------------------------------------------
  

!         Appel du calcul de nebulosite.

  IF(LNEBN) THEN
    CALL ACNEBN ( YDMODEL%YRML_PHY_MF,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,NTNEBU,YDCPG_DIM%KFLEVG,&
     & YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,ZQV,ZQL,ZQI,YDMF_PHYS%TMP%APLPAR%FLU%QSAT,YDMF_PHYS_STATE%T,YDMF_PHYS%TMP%APLPAR%PFL%FPLCH,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,&
     & ZUNEBH,&
     & ZNEBS0,ZQLIS0,ZQLI_CVP,ZNEB_CVPP,ZQLI_CVPP,&
     & ZAIPCMT,YDCPG_MISC%NEB,ZNEBC0,YDCPG_MISC%QICE,YDCPG_MISC%QLI,&
     & ZHUC,ZVETAF)
    DO JLEV=KTDIA,YDCPG_DIM%KFLEVG
!DEC$ IVDEP
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        YDCPG_MISC%NEB(JLON,JLEV)=YDCPG_MISC%NEB(JLON,JLEV)*GAEPS
      ENDDO
    ENDDO
  ENDIF

  

!         Diagnostique de nebulosite partielle.
  CALL ACNPART(YDMODEL%YRML_PHY_MF,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,NTNEBU,YDCPG_DIM%KFLEVG,&
   & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,ZDECRD,YDCPG_MISC%NEB,&
   & YDMF_PHYS%OUT%CLCH,YDMF_PHYS%OUT%CLCM,YDMF_PHYS%OUT%CLCL,YDCPG_MISC%CLCT,ZCLCT_RAD,&
   ! optional arguments (convective cloud cover)
   & PCLCC=YDMF_PHYS%OUT%CLCC,PNEBC=ZNEBC0,PTOPC=YDMF_PHYS%OUT%CTOP)

!     7.3.5 Computation of the equivalent coefficients for simplified
!           radiation scheme

  

!     7.3.6 Module chimique - Chemistry module

  IF (LCHEM_ARPCLIM) THEN ! at this stage call when ARPEGE-Climat chemistry only

     ! initialisation below needs to be refined later for more general use
     IFLDX  = 1_JPIM
     IFLDX2 = 1_JPIM
     ILEVX  = 1_JPIM
     ALLOCATE(ZSD_XA(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,IFLDX), ZSD_X2(YDCPG_DIM%KLON,IFLDX2))
     ALLOCATE(INDCHEM(NCHEM),IGPLAT(YDCPG_DIM%KLON))
     ALLOCATE(ZCFLX(YDCPG_DIM%KLON,NCHEM), ZCFLXO(YDCPG_DIM%KLON,NCHEM), ZCHEMDV(YDCPG_DIM%KLON,0)) ! no species with dry deposition in ARPCLIM
     ALLOCATE(ZAEROP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,NACTAERO),ZTENC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,NCHEM))
     ALLOCATE(ZDELP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZWND(YDCPG_DIM%KLON),ZDUMMY1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZGELAT(YDCPG_DIM%KLON))
     ALLOCATE(ZNEEFLX(YDCPG_DIM%KLON),ZCHEM2AER(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,4))
     INDCHEM(:)     = 1_JPIM ! we should have here the indexes of the chemical species in the YCHEM array, to be implemented later
     IGPLAT (:)     = 1_JPIM
     ZSD_XA (:,:,:) = 0._JPRB
     ZSD_X2 (:,:)   = 0._JPRB
     DO JLEV=KTDIA,YDCPG_DIM%KFLEVG
       DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
         ZDELP(JLON,JLEV) = YDMF_PHYS_STATE%YCPG_PHY%PREHYD(JLON,JLEV) - YDMF_PHYS_STATE%YCPG_PHY%PREHYD(JLON,JLEV-1)
       ENDDO
     ENDDO
     ZWND  (:)      = 0._JPRB  ! not used in  ARPEGE-Climat chemistry
     ZDUMMY1 (:,:)  = 1.0E-18_JPRB
     ZGELAT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA) = ASIN(PGEMU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA))
     ZCFLX (:,:)        = 0._JPRB
     ZCFLXO (:,:)       = 0._JPRB
     ZCHEMDV (:,:)      = 0._JPRB
     ZAEROP (:,:,:)     = 0._JPRB  ! no interaction aerosol/chemistry in ARPCLIM
     ZTENGFL(:,:,:)     = 0._JPRB  ! only used later for diagnostics
     ZTENC (:,:,:)      = 0._JPRB  ! not used in  ARPEGE-Climat chemistry
     ZNEEFLX (:)        = 0._JPRB  ! not used in  ARPEGE-Climat chemistry
     ZCHEM2AER (:,:,:)  = 0._JPRB  ! not used in  ARPEGE-Climat chemistry
     
!#if false
! !   YDVAB YDDIMV not used for ARPEGE-Climat chemistry
!      CALL CHEM_MAIN &
!    &( YDVAB, YDDIMV, YDMODEL, KIDIA , KFDIA , KLON , KLEV, KVCLIS, NCHEM, INDCHEM,&
!    &  TSPHY , IGPLAT,  IFLDX  , IFLDX2 , ILEVX,&
!    &  ZSD_XA , ZSD_X2,  ZDELP, PAPRS, PAPRSF, PAPHI, PQ, PT,&
!    &  ZDUMMY1, ZDUMMY1, ZDUMMY1, ZDUMMY1,  PNEB,& 
!    &  PFPLCL, PFPLCN, PFPLSL, PFPLSN, ZDUMMY1,& 
!    &  PALB, ZWND, PLSM,& 
!    &  PMU0, ZGELAT, PGELAM, PGEMU, PKOZO, ZCFLX, ZCFLXO, ZCHEMDV, PGFL,&
!    &  ZAEROP, ZTENGFL, PCHEM, ZTENC, ZNEEFLX, ZCHEM2AER )
!#endif     
  ENDIF

  
!      7.4 Rayonnement Geleyn
!          Geleyn's radiation

  

    LLCALLRAD=(MOD(YDCPG_DIM%KSTEP,NRADFR) == 0 )
!    IF (NCALLRAD==1) ! <== not yet
    IF (NCALLRAD==2) LLCALLRAD=(LLCALLRAD.AND.(YDCPG_DIM%KSTEP<=NSTOP-1))
!    IF (NCALLRAD==3) ! <== not yet

    ! ---- Intermittent call to radiation scheme
    IF (LLCALLRAD) THEN
      CALL RECMWF(YDGEOMETRY%YRDIMV,YDMODEL,     &
       &  YDCPG_DIM%KIDIA , YDCPG_DIM%KFDIA, YDCPG_DIM%KLON  , YDCPG_DIM%KFLEVG   ,         &
       &  ZALBD , ZALBP, YDMF_PHYS_STATE%YCPG_PHY%PREHYD , YDMF_PHYS_STATE%YCPG_PHY%PREHYDF ,         &
       &  YDCPG_MISC%NEB  , ZQO3 , ZAER  , YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP  , YDMF_PHYS%TMP%APLPAR%FLU%EMIS , &
       &  YDMF_PHYS%TMP%RDG%MU0M , ZQV  , YDMF_PHYS%TMP%APLPAR%FLU%QSAT , YDCPG_MISC%QICE  , YDCPG_MISC%QLI  , & 
       &  ZQS   , ZQR  , YDMF_PHYS_SURF%GSD_VF%PLSM  , YDMF_PHYS_STATE%T     , YDMF_PHYS_STATE%YGSP_RR%T   , &
       &  PGP2DSPP, PEZDIAG, &
       &  YDMF_PHYS%RAD%EMTD , YDMF_PHYS%RAD%EMTU, YDMF_PHYS%RAD%TRSW ,                   &
       &  YDMF_PHYS%OUT%FRTHC, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%FRSOC   , YDMF_PHYS%OUT%FRSO ,&
       &  ZSFSWDIR     , ZSFSWDIF , ZFSDNN , ZFSDNV,&
       &  ZCTRSO,ZCEMTR, ZTRSOD ,&
       &  ZTRSODIR, ZTRSODIF,&
       &  ZPIZA_DST,ZCGA_DST,ZTAUREL_DST,ZAERINDS,&
       &  PGELAM, PGEMU, PGPAR, YDMF_PHYS%TMP%RDG%MU0LU , YDMF_PHYS%OUT%ALB, YDMF_PHYS%RAD%RMOON)
    ELSE
      IF (LMSE) THEN
      DO JSG=1,NSW
        ZTRSODIR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSG)=PGPAR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,MSWDIR+JSG-1)
        ZTRSODIF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSG)=PGPAR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,MSWDIF+JSG-1)
      ENDDO
      ENDIF
    ENDIF

    YDMF_PHYS%OUT%FRSOLU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS%RAD%RMOON(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)

    ! ---- Flux update and radiative heating rates
    CALL RADHEAT ( YDERAD,YDERDI,YDMODEL%YRML_PHY_MF, &
    &   YDCPG_DIM%KIDIA  , YDCPG_DIM%KFDIA  , YDCPG_DIM%KLON    , YDCPG_DIM%KFLEVG,&
    &   YDMF_PHYS_STATE%YCPG_PHY%PREHYD  , YDMF_PHYS%TMP%APLPAR%FLU%EMIS  , YDMF_PHYS%RAD%EMTD   , YDMF_PHYS%TMP%RDG%MU0, ZQV  ,&
    &   ZTENT  , YDMF_PHYS%RAD%TRSW   , ZTRSOD , YDMF_PHYS_STATE%YGSP_RR%T , TSPHY,&
    &   ZTRSODIR, ZTRSODIF, ZALBD , ZALBP,&
    &   YDMF_PHYS%OUT%FRSO  , YDMF_PHYS%OUT%FRTH  , YDMF_PHYS%OUT%FRSODS , YDMF_PHYS%OUT%FRTHDS,&
    &   ZCEMTR , ZCTRSO , YDMF_PHYS%OUT%FRSOC  , YDMF_PHYS%OUT%FRTHC,&
    &   ZSUDU  , ZSDUR  , ZDSRP  , ZSFSWDIR , ZSFSWDIF, YDMF_PHYS%OUT%FRSOPS, ZFRSODS, YDMF_PHYS%OUT%FRSOPT )


    ! ---- Take into account day duration depending on altitude.
    IF(LDAYD) CALL ACDAYD(YDRIP,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,KTDIA,KSGST,PGEMU ,YDMF_PHYS%TMP%RDG%MU0,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDMF_PHYS%OUT%FRSO)

    ! ---- Correct solar absorption as a function of pmu0.
    IF(GRSO < 1._JPRB) CALL ACRSO(YDPHY0,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,KTDIA,KSGST,PGEMU  ,YDMF_PHYS%TMP%RDG%MU0,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDMF_PHYS%OUT%FRSO)

    IF(.NOT.LMSE) THEN
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZALB(JLON)=0.0_JPRB
        DO JSG=1,NSW
          ZALB(JLON)=ZALB(JLON)+0.5_JPRB*(ZALBD(JLON,JSG)+ZALBP(JLON,JSG))
        ENDDO
        ZALB(JLON)=ZALB(JLON)/FLOAT(NSW)
        YDMF_PHYS%OUT%FRSODS(JLON)=YDMF_PHYS%OUT%FRSO(JLON,YDCPG_DIM%KFLEVG,1)/(1.0_JPRB-ZALB(JLON))
      ENDDO
    ENDIF
    ! Compute Sunshine Duration (in seconds)
!DEC$ IVDEP
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      IF(YDMF_PHYS%OUT%FRSODS(JLON) >= RSUNDUR) THEN
        YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)=YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)+1.0_JPRB*TSTEP
      ENDIF
    ENDDO

  

  

    ! Direct normal irradiance with securities
    DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%FRSDNI(JLON)=YDMF_PHYS%OUT%FRSOPS(JLON)
      IF (YDMF_PHYS%TMP%RDG%MU0(JLON) > 3.0E-02_JPRB) THEN
        YDMF_PHYS%OUT%FRSDNI(JLON)=YDMF_PHYS%OUT%FRSOPS(JLON)/YDMF_PHYS%TMP%RDG%MU0(JLON)
      ENDIF
      YDMF_PHYS%OUT%FRSDNI(JLON)=MAX(0.0_JPRB,YDMF_PHYS%OUT%FRSDNI(JLON))
    ENDDO
  

  

    ! global normal irradiance
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%FRSGNI(JLON)=YDMF_PHYS%OUT%FRSDNI(JLON)+0.5_JPRB*(                         &
       & (1.0_JPRB+YDMF_PHYS%TMP%RDG%MU0(JLON))*(YDMF_PHYS%OUT%FRSODS(JLON)-YDMF_PHYS%OUT%FRSOPS(JLON)       )+ &
       & (1.0_JPRB-YDMF_PHYS%TMP%RDG%MU0(JLON))*(YDMF_PHYS%OUT%FRSODS(JLON)-YDMF_PHYS%OUT%FRSO  (JLON,YDCPG_DIM%KFLEVG,1)))
    ENDDO

    ! mean radiant temperature
    

  

!*
!     ------------------------------------------------------------------

!     7.BIS. BILAN HYDRIQUE DU SOL
!     ----------------------------
!     CALCUL DES RESISTANCES A L'EVAPOTRANSPIRATION HV ET
!     A LA TRANSPIRATION
!     ------------------------------------------------------------------
!     HTR DU COUVERT VEGETAL
!     ----------------------

  IF (LSOLV.AND.(.NOT.LMSE)) THEN
    CALL ACVEG ( YDPHY,YDPHY1,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,&
     & YDMF_PHYS%OUT%FRSO,ZQV,YDMF_PHYS%TMP%APLPAR%FLU%QSAT,YDMF_PHYS_STATE%T,&
     & YDMF_PHYS_SURF%GSD_VV%PD2,YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS_SURF%GSD_VV%PIVEG,YDMF_PHYS_SURF%GSD_VV%PLAI,YDMF_PHYS%TMP%APLPAR%FLU%NEIJ,YDMF_PHYS%TMP%APLPAR%FLU%VEG,YDMF_PHYS_SURF%GSD_VV%PRSMIN,&
     & ZCHROV,ZGWDCS,ZWFC,YDMF_PHYS_STATE%YGSP_RR%FC,ZWWILT,YDMF_PHYS_STATE%YGSP_SB%Q,YDMF_PHYS%TMP%APLPAR%FLU%QSATS,&
     & ZHQ,ZHTR,ZHU,YDMF_PHYS_SURF%GSD_VV%PHV,ZWLMX)
  ENDIF

!     ------------------------------------------------------------------
!     8.- DIFFUSION VERTICALE TURBULENTE
!     ----------------------------------
  

!   Sauvegarde temporaire de l'ancien acdifus pour les besoins du Climat
    

      

      CALL ACDIFV1 ( YGFL,YDMODEL%YRML_PHY_MF,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,KTDIA,YDCPG_DIM%KFLEVG,&
        & YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,ZCP,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,ZKTROV,ZKQROV,ZKUROV,&
        & ZQV,ZQL,ZQI,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,ZSVM,ZXTROV,ZXUROV,&
        & ZXDROV,ZXHROV,&
        & ZEDMFS,ZEDMFQ,ZEDMFU,ZEDMFV,ZMF_UP, &
        & ZXURO,ZXQRO,ZXTRO,&
        & ZCFAQ,ZCFAS,ZCFATH,ZCFAU,ZCFASV,&
        & ZCFBQ,ZCFBS,ZCFBTH,ZCFBU,ZCFBV,ZCFBSV,&
        & ZDSE,ZQT)   
        

      IF ( LMSE.AND.LCALLSFX ) THEN

        
          ZCARDI=RCARDI
        

        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZRHODREFM(JLON)=YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(JLON,YDCPG_DIM%KFLEVG)/(YDMF_PHYS_STATE%T(JLON,YDCPG_DIM%KFLEVG)*YDMF_PHYS_STATE%YCPG_DYN%RCP%R(JLON,YDCPG_DIM%KFLEVG))
          ZDEPTH_HEIGHT(JLON,:)=(YDMF_PHYS_STATE%YCPG_DYN%PHIF(JLON,:)-YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,YDCPG_DIM%KFLEVG))/RG
          ZZS(JLON)=YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,YDCPG_DIM%KFLEVG)/RG
        ENDDO

        IRR=2

        CALL ARO_GROUND_PARAM( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP,&
          & YDCPG_DIM%KFDIA-YDCPG_DIM%KIDIA+1, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KSTEP,&
          & IRR, NSW, NGFL_EXT,NDGUNG, NDGUXG, NDLUNG, NDLUXG,LSURFEX_KFROM,&
          & LMPA, CCOUPLING,LDXFUMSE, &
          & NINDAT, ZRHGMT,ZSTATI,RSOVR,RCODEC,RSIDEC, &
          & PINDX(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), PINDY(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), &
          & YDMF_PHYS_STATE%U(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),YDMF_PHYS_STATE%V(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG), &
          & YDMF_PHYS_STATE%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),YDMF_PHYS_STATE%Q(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),&
          & ZSVM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG,1:NGFL_EXT),&
          & ZCARDI,&
          & ZRHODREFM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),YDMF_PHYS_STATE%YCPG_PHY%PREHYD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),&
          & ZDTMSE,ZDEPTH_HEIGHT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),ZZS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), XZSEPS,&
          & YDMF_PHYS%TMP%RDG%MU0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),YDMF_PHYS%TMP%RDG%MU0N(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),PGELAM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & PGEMU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),XSW_BANDS,&
          & ZSRAIN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),ZSSNOW(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & ZSGROUPEL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & YDMF_PHYS%OUT%FRTHDS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),ZSFSWDIF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:NSW),&
          & ZSFSWDIR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:NSW),&
          & ZCFAQ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG), ZCFATH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),&
          & ZCFAU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),&
          & ZCFBQ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG), ZCFBTH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),&
          & ZCFBU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),&
          & ZCFBV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),&
          & YDMF_PHYS%OUT%FCS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1),ZFEV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & ZSFSV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:NGFL_EXT),ZSFCO2(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & ZFMDU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),ZFMDV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & ZALBP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:NSW),ZALBD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:NSW),&
          & YDMF_PHYS%TMP%APLPAR%FLU%EMIS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),ZTSN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & YDMF_PHYS%OUT%FRTH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,1))               !orographic shadowing

! Opposite water vapor flux
        ZFEVS(:)=-ZFEV(:)

        CALL ARO_GROUND_DIAG( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP,&
          & YDCPG_DIM%KFDIA-YDCPG_DIM%KIDIA+1, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, -1,&
          & NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM,&
          & ZZS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),ZFEVS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & YDMF_PHYS_STATE%U(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG), YDMF_PHYS_STATE%V(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG),&
          & ZDEPTH_HEIGHT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG),&
          & YDMF_PHYS%OUT%FRTH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,1),YDMF_PHYS%OUT%FRSO(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,1),&
          & PINDX(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), PINDY(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & YDCPG_MISC%QS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%GZ0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & YDMF_PHYS%OUT%GZ0H(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%TCLS (YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & YDMF_PHYS%OUT%QCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%RHCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & YDMF_PHYS%OUT%UCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%VCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & YDMF_PHYS%OUT%NUCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%NVCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), &
          & YDMF_PHYS%OUT%FCLL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1), YDMF_PHYS%OUT%FCLN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1),&
          & YDMF_PHYS%OUT%FEVL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1), YDMF_PHYS%OUT%FEVN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1),&
          & ZSSO_STDEV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZTWSNOW(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & ZBUDTH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZBUDSO(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & ZFCLL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),ZTOWNS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),&
          & ZCD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)                        )

        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          PGPAR(JLON,MGZ0)   = YDMF_PHYS%OUT%GZ0 (JLON)
          PGPAR(JLON,MGZ0H)  = YDMF_PHYS%OUT%GZ0H(JLON)
          PGPAR(JLON,MVEMIS) = YDMF_PHYS%TMP%APLPAR%FLU%EMIS(JLON)
          PGPAR(JLON,MVQS)   = YDCPG_MISC%QS  (JLON)
          PGPAR(JLON,MCD)    = ZCD  (JLON)
        ENDDO

        DO JSG=1,NSW
          PGPAR(:,MALBSCA-1+JSG) = ZALBD(:,JSG)
          PGPAR(:,MALBDIR-1+JSG) = ZALBP(:,JSG)
        ENDDO

        DO JSG  = 1, KSGST+1
          DO JLEV = 0, YDCPG_DIM%KFLEVG
            DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
              YDMF_PHYS%OUT%FRTH(JLON,JLEV,JSG)=YDMF_PHYS%OUT%FRTH(JLON,JLEV,JSG)+ZBUDTH(JLON)
            ENDDO
          ENDDO
        ENDDO

! calculation of variables for the old "ISBA" atmosphere scheme
        
          CALL ARO_GROUND_DIAG_2ISBA( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP, &
            & YDCPG_DIM%KFDIA-YDCPG_DIM%KIDIA+1, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, &
            & NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM, &
            & PINDX(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), PINDY(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), &
            & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PARG, YDMF_PHYS_SURF%GSD_VV%PSAB, YDMF_PHYS_SURF%GSD_VV%PD2, ZTSN, ZTWSNOW, &
            & YDMF_PHYS_SURF%GSP_SB%PT_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1), YDMF_PHYS_SURF%GSP_RR%PW_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS_SURF%GSP_SB%PQ_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1), &
            & YDMF_PHYS_SURF%GSP_RR%PIC_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS_SURF%GSP_SB%PTL_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1), YDMF_PHYS_SURF%GSP_RR%PFC_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), &
            & YDMF_PHYS_SURF%GSP_SG%PA_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1), YDMF_PHYS_SURF%GSP_SG%PR_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1), ZHV2 )
          DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            YDMF_PHYS_SURF%GSD_VV%PHV(JLON) = ZHV2(JLON)
          ENDDO
        

        IF (.NOT.LNODIFQC) THEN
          CALL ACAA1 ( YDMODEL%YRML_PHY_MF,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,KTDIA,YDCPG_DIM%KFLEVG,&
            & YDMF_PHYS_STATE%YCPG_PHY%PREHYD,ZCOEFN,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,ZQL,ZQI,YDMF_PHYS_STATE%T,&
            & ZALPHA1,ZCOEFA,ZLVT,ZQICE)
        ENDIF

      ELSEIF (.NOT.LMSE) THEN

        IF(NDIFFNEB == 5) THEN
          ZNEBDIFF(:,:)=ZNEBS(:,:)
        ELSEIF(NDIFFNEB == 2) THEN
          ZNEBDIFF(:,:)=YDCPG_MISC%NEB(:,:)
        ELSEIF(NDIFFNEB == 3) THEN
          ZNEBDIFF(:,:)=ZNEBS(:,:)+(1.0_JPRB-ZNEBS(:,:))*ZNEBCH(:,:)
        ENDIF

        
             
          CALL ARP_GROUND_PARAM ( YDMCC,YDMODEL%YRML_PHY_MF,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,KTDIA,YDCPG_DIM%KFLEVG,KCSS,&
           & YDMF_PHYS_STATE%YCPG_DYN%PHI,ZCP,YDMF_PHYS%OUT%FRSO,&
           & ZQV,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,ZXURO,ZXQRO,ZXTRO,&
           & YDMF_PHYS_STATE%YCPG_PHY%PREHYD,ZQL,ZQI,ZNEBDIFF,&
           & YDMF_PHYS%TMP%APLPAR%FLU%CD,YDMF_PHYS%TMP%APLPAR%FLU%CDN,ZCDROV,ZCHROV,ZCEROV,YDMF_PHYS%TMP%APLPAR%DSA%CPS,YDMF_PHYS%OUT%CT,ZDQSTS,YDMF_PHYS%TMP%APLPAR%FLU%EMIS,YDMF_PHYS_SURF%GSD_VH%PSPSH,&
           & ZHQ,ZHTR,ZHU,YDMF_PHYS_SURF%GSD_VV%PHV,YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS_SURF%GSD_VV%PIVEG,YDMF_PHYS%TMP%APLPAR%FLU%NEIJ,YDCPG_MISC%QS,YDMF_PHYS%TMP%APLPAR%FLU%QSATS,&
           & YDMF_PHYS_STATE%YGSP_SB%T,YDMF_PHYS_STATE%YGSP_RR%T,YDMF_PHYS%TMP%APLPAR%FLU%VEG,ZXDROV,ZXHROV,YDMF_PHYS_STATE%YGSP_RR%W,YDMF_PHYS_STATE%YGSP_RR%IC,&
           & ZDSE,&
           & ZCFAS,ZCFAU,&
           & ZCFBS,ZCFBU,ZCFBV,&
           & ZCFBQ,ZCOEFA,ZALPHA1,ZLVT,ZQICE,&
           & ZDIFWQ,ZDIFWS,ZFMDU,ZFMDV,&
           & ZSC_FEVI,ZSC_FEVN,ZSC_FCLL,ZSC_FCLN,&
           & YDMF_PHYS%OUT%FCHSP,YDMF_PHYS%OUT%FCLL,YDMF_PHYS%OUT%FCLN,YDMF_PHYS%OUT%FCS,YDMF_PHYS%TMP%APLPAR%FLU%FEVI,YDMF_PHYS%OUT%FEVL,YDMF_PHYS%OUT%FEVN,&
           & YDMF_PHYS%OUT%FEVV,YDMF_PHYS%OUT%FTR,YDMF_PHYS%TMP%APLPAR%DSA%LHS,YDMF_PHYS_SURF%GSD_VH%PQSH,&
           & YDMF_PHYS%OUT%FRTH,&
           & YDMF_PHYS_SURF%GSD_VF%PZ0F,YDMF_PHYS_SURF%GSD_VV%PZ0H)

          DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            YDMF_PHYS%OUT%DIFTQ(JLON,YDCPG_DIM%KFLEVG)=ZDIFWQ(JLON)
            YDMF_PHYS%OUT%DIFTS(JLON,YDCPG_DIM%KFLEVG)=ZDIFWS(JLON)
          ENDDO

             

           ! <== LSFORCS

      ELSEIF (.NOT. LCALLSFX) THEN
        YDMF_PHYS%OUT%FCS(:,:) = 0.0_JPRB
        ZFEV(:)   = 0.0_JPRB

      ENDIF  !LMSE.AND.LCALLSFX

      CALL ACDIFV2 ( YGFL,YDMODEL%YRML_PHY_MF,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,KTDIA,YDCPG_DIM%KFLEVG,&
        & YDMF_PHYS_STATE%YCPG_PHY%PREHYD,ZCFAQ,ZCFAS,ZCFAU,ZCFASV,ZCFBQ,ZCFBS,ZCFBU,ZCFBV,ZCFBSV,&
        & ZKTROV,ZKQROV,ZKQLROV,ZKUROV,ZDSE,ZQT,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,ZPOID,&
        & YDMF_PHYS_STATE%T,ZQL,ZQI,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,&
        & ZCOEFA,ZALPHA1,ZLVT,ZQICE,&
        & ZSFSV,YDMF_PHYS%OUT%FCS,ZFEV,ZFMDU,ZFMDV,ZTSN,ZXHROV,&
        & PDIFSV,YDMF_PHYS%OUT%DIFTQ,YDMF_PHYS%OUT%DIFTS,YDMF_PHYS%OUT%STRTU,YDMF_PHYS%OUT%STRTV,YDMF_PHYS%OUT%DIFTQL,YDMF_PHYS%OUT%DIFTQN,YDMF_PHYS%OUT%DIFCQ,YDMF_PHYS%OUT%DIFCS,YDMF_PHYS%OUT%STRCU,YDMF_PHYS%OUT%STRCV,YDVARS%SHTUR%T0)

      

      IF ( LMSE ) THEN
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          PGPAR(JLON,MVTS) = ZTSN (JLON)
          YDMF_PHYS_SURF%GSP_SG%PF_T1(JLON,1) = ZTWSNOW (JLON)
        ENDDO
      ENDIF

      ! First compute horizontal exchange coefficients for momentum:
      !  (there's mo TOMs contribution, thus has to be done at latest here)
      

       ! LCOEFKTKE

      ! Now the heat coefficient can be completed by TKE+ containing
      !  the TOMs contribution. 
      

    

!-----------------------------------------------------------------------------
!   THE DEEP CONVECTION WILL SEE THE SHALLOW PART FROM KFB AS IT IS WITH LOUIS
!   SCHEME AND THE MODIFIED RI
!----------------------------------------------------------------------------
    
      DO JLEV=KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          YDMF_PHYS%OUT%DIFTQ   (JLON,JLEV) = YDMF_PHYS%OUT%DIFTQ(JLON,JLEV) + ZDIFCVPPQ(JLON,JLEV)
          YDMF_PHYS%OUT%DIFTS   (JLON,JLEV) = YDMF_PHYS%OUT%DIFTS(JLON,JLEV) + ZDIFCVPPS(JLON,JLEV)
          YDMF_PHYS%OUT%STRTU   (JLON,JLEV) = YDMF_PHYS%OUT%STRTU(JLON,JLEV) + ZDIFCVPPU(JLON,JLEV)
          YDMF_PHYS%OUT%STRTV   (JLON,JLEV) = YDMF_PHYS%OUT%STRTV(JLON,JLEV) + ZDIFCVPPV(JLON,JLEV)
        ENDDO
      ENDDO
    

     ! 3MT || LSTRAPRO

    ! ------------------------------------------------------
    ! DIAGNOSTIC DE LA HAUTEUR DE COUCHE LIMITE, SELON TROEN ET MAHRT,
    ! POUR USAGE AU PAS DE TEMPS SUIVANT.
    ! ------------------------------------------------------

    

    
    ! ------------------------------------------------------------------
    ! UPDATE PASSIFS SCALAIRS DUE TO THE TURBULENT DIFFUSION PROCESSES
    ! ------------------------------------------------------------------
      ! LMDUST
   ! LVDIF

!     DIAGNOSTIC SUPPLEMENTAIRE FLUX DE RAYONNEMENT (PFRTHDS ET PFRSOPT)
!     ADDITIONAL DIAGNOSTICS OF RADIATIVE FLUXES (PFRTHDS AND PFRSOPT)

!DEC$ IVDEP
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    IF(.NOT.LMSE)YDMF_PHYS%OUT%FRTHDS(JLON)=YDMF_PHYS%OUT%FRTH(JLON,YDCPG_DIM%KFLEVG,1)/YDMF_PHYS%TMP%APLPAR%FLU%EMIS(JLON)+RSIGMA*YDMF_PHYS_STATE%YGSP_RR%T(JLON)**4
    IF(YDMF_PHYS%TMP%RDG%MU0(JLON) <= 0.0_JPRB) THEN


      YDMF_PHYS%OUT%FRSOPT(JLON)=0.0_JPRB
    ELSE
      YDMF_PHYS%OUT%FRSOPT(JLON)=RII0*YDMF_PHYS%TMP%RDG%MU0(JLON)
    ENDIF
  ENDDO

! ADDITIONAL DIAGNOSTIC OF THE DERIVATIVE OF THE NON SOLAR SURFACE
! HEAT FLUX WITH RESPECT TO SURFACE TEMPERATURE (PDERNSHF)

  IF(LMCC03)THEN
    CALL ACDNSHF(YDPHY,YDPHY1,YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, KTDIA, YDCPG_DIM%KFLEVG,&
     & YDMF_PHYS%TMP%APLPAR%FLU%EMIS, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS%TMP%APLPAR%FLU%NEIJ, ZQV, YDCPG_MISC%QS, YDMF_PHYS_STATE%YGSP_RR%T, ZCHROV,  ZDQSTS,&
     & YDMF_PHYS%OUT%DRNSHF)
  ENDIF

!*
!     ------------------------------------------------------------------
!     9.- TRAINEE DES ONDES DE GRAVITE INDUITES PAR LE RELIEF
!     ------------------------------------------------------------------
  
    CALL ACDRAG ( YDMODEL%YRML_PHY_MF,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,NTDRAG,YDCPG_DIM%KFLEVG,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,&
     & YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP, ZNBVNO, YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V,&
     & PRCORI,YDMF_PHYS_SURF%GSD_VF%PGETRL, ZGWDCS,&
     & YDMF_PHYS_SURF%GSD_VF%PVRLAN , YDMF_PHYS_SURF%GSD_VF%PVRLDI,&
     & YDMF_PHYS%OUT%STRDU, YDMF_PHYS%OUT%STRDV,ZTRAJGWD)
  
 ! SAVE FOR TL/NL COEFS FROM VERT. DIFF AND GWD

  

!     ------------------------------------------------------------------
!     10.- PRECIPITATIONS STRATIFORMES.
!     ---------------------------------

  

  IF ( LSTRAS ) THEN
    CALL ACPLUIS ( YDMODEL%YRML_PHY_MF,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,NTPLUI,YDCPG_DIM%KFLEVG,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,ZNEBS,ZQV,ZQLIS,YDMF_PHYS%TMP%APLPAR%MSC%QW,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%T,&
     & YDMF_PHYS%OUT%FCSQL,YDMF_PHYS%OUT%FCSQN,YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FPLSN,YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPFPSL,YDMF_PHYS%OUT%FPFPSN)
  ENDIF

  

   ! L3MT

!     11.- PRECIPITATIONS SOUS-MAILLES.
!     ---------------------------------

  


   ! L3MT

! -----------------------------------------------------

       ! <== IFS deep convection scheme

    DO JLEV=1,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZGEOM1(JLON,JLEV)     = YDMF_PHYS_STATE%YCPG_DYN%PHIF(JLON,JLEV)-YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,YDCPG_DIM%KFLEVG)
        ZVERVEL(JLON,JLEV)    = YDCPG_DYN0%CTY%VVEL(JLON,JLEV)*YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(JLON,JLEV)
        ZLISUM(JLON,JLEV)     = 0.0_JPRB
        ZLCRIT_AER(JLON,JLEV) = 5.E-4_JPRB
      ENDDO
    ENDDO
    DO JLEV=0,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZGEOMH(JLON,JLEV)=YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,JLEV)-YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,YDCPG_DIM%KFLEVG)
        ZFCQLF(JLON,JLEV)=0.0_JPRB
        ZFCQLI(JLON,JLEV)=0.0_JPRB
      ENDDO
    ENDDO
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
       LLLAND(JLON)    = (YDMF_PHYS_SURF%GSD_VF%PLSM(JLON)>0.5_JPRB)
       LLSHCV(JLON)    = .FALSE.
       ZCUCONVCA(JLON) = 0.0_JPRB
       ZACPR(JLON)     = 0.0_JPRB
       ZDXTDK(JLON)    = RDELXN/PGM(JLON)
    ENDDO
    LLPTQ = .FALSE.
    CALL ACUPTQ ( YDCPG_DIM%KLON,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,LLPTQ,YDMF_PHYS%OUT%FRSO,YDMF_PHYS%OUT%FRTH,&
     & YDMF_PHYS%OUT%DIFCQ,YDMF_PHYS%OUT%DIFCS,YDMF_PHYS%OUT%DIFTQ,YDMF_PHYS%OUT%DIFTS,YDMF_PHYS%OUT%FCCQL,YDMF_PHYS%OUT%FCCQN,YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,&
     & YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,YDMF_PHYS%OUT%FPFPCL,YDMF_PHYS%OUT%FPFPCN,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%Q,YDMF_PHYS_STATE%YGSP_RR%T,&
     & ZTENHA,ZTENQVA )
    
    DO JLEV=1,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA 
         ZTENT(JLON,JLEV) = ZTENHA(JLON,JLEV)/YDMF_PHYS_STATE%YCPG_DYN%RCP%CP(JLON,JLEV)
         ZTENQ(JLON,JLEV) = ZTENQVA(JLON,JLEV)
         ZTENU(JLON,JLEV) = 0.0_JPRB
         ZTENV(JLON,JLEV) = 0.0_JPRB
      ENDDO
    ENDDO
    LLDSLPHY=.TRUE.
    ZVDIFTS = 0._JPRB
    ISPPN2D = 0
         
    CALL CUCALLN_MF &
     & (YDMODEL%YRML_PHY_RAD%YRERAD,YDMODEL%YRML_PHY_SLIN,YDMODEL%YRML_PHY_EC,YDMODEL%YRML_GCONF%YGFL,&
     & YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,0,YDCPG_DIM%KFLEVG,ZDXTDK,ISPPN2D,LLLAND,LLDSLPHY,TSPHY,ZVDIFTS,&
     & YDMF_PHYS_STATE%T,ZQV,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,ZLISUM,ZVERVEL,YDMF_PHYS%OUT%DIFTQ,YDMF_PHYS%OUT%DIFTS,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,ZGEOM1,ZGEOMH,PGM, &
     & ZCUCONVCA,ZGP2DSPP,ZTENT,ZTENQ,ZTENU,ZTENV,ZACPR,&
     & ITOPC,IBASC,ITYPE,ICBOT,ICTOP,IBOTSC,LLCUM,LLSC,LLSHCV,ZLCRIT_AER,&
     & ZLU,ZLUDE,ZLUDELI,ZSNDE,ZMFU,ZMFD,YDMF_PHYS%OUT%DIFCQ,YDMF_PHYS%OUT%DIFCS,ZFHPCL,ZFHPCN,&
     & YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,ZLRAIN,ZRSUD,YDMF_PHYS%OUT%STRCU,YDMF_PHYS%OUT%STRCV,ZFCQLF,ZFCQLI,&
     & ZMFUDE_RATE,ZMFDDE_RATE,YDMF_PHYS%OUT%CAPE,ZWMEAN,ZDIFF,0,ZCEN,ZTENC,ZSCAV)
    DO JLEV=0,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA 
         ZFPCOR  (JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)+YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
         YDMF_PHYS%OUT%FCCQL  (JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)
         YDMF_PHYS%OUT%FCCQN  (JLON,JLEV)=YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
         YDMF_PHYS%OUT%FPFPCL (JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)
         YDMF_PHYS%OUT%FPFPCN (JLON,JLEV)=YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
         YDMF_PHYS%OUT%FPEVPCL(JLON,JLEV)=0._JPRB
         YDMF_PHYS%OUT%FPEVPCN(JLON,JLEV)=0._JPRB
       ENDDO
    ENDDO

  
  IF(LFLASH) THEN
    ! LIGHTNING PARAMETERIZATION. OUTPUT IS PFLASH, TOTAL LIGHTNING FLASH RATES.
    ZGAW(:)=0._JPRB
    CALL CULIGHT &
      & ( YDEPHY,YGFL,YDCPG_DIM%KIDIA ,  YDCPG_DIM%KFDIA,   YDCPG_DIM%KLON,    YDCPG_DIM%KFLEVG,  ZGAW, ZGAW, &
      &   YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,     YDMF_PHYS_STATE%YCPG_PHY%PREHYD  ,  YDMF_PHYS_STATE%YCPG_DYN%PHI,   YDMF_PHYS_STATE%YCPG_DYN%PHIF,  LLLAND,&
      &   YDMF_PHYS_STATE%T    ,  ZLU   ,  ZMFU,    YDMF_PHYS%OUT%CAPE, &
      &   YDMF_PHYS%OUT%FPLCL,  YDMF_PHYS%OUT%FPLCN,&
      &   LLCUM ,  ICBOT ,  ICTOP,&
      ! Outputs
      &   LLLINOX, YDMF_PHYS%OUT%FLASH,  ZLIGH_CTG,  ZCTOPH, &
      &   ZPRECMX, ZICE,   ZCDEPTH,  ZWMFU) 
    ! LIGHTNING FLASH RATES ARE CONVERTED IN fl/km2/s BEFORE ENTERING CFU TIME ACCUMULATION.
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%FLASH(JLON)=YDMF_PHYS%OUT%FLASH(JLON)/86400._JPRB
    ENDDO
  ENDIF

  
    ZFRSO(:,:) = 0.0_JPRB
    ZFRTH(:,:) = 0.0_JPRB
    LLPTQ = .TRUE.
    CALL ACUPTQ ( YDCPG_DIM%KLON,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,LLPTQ,ZFRSO,ZFRTH,&
     & YDMF_PHYS%OUT%DIFCQ,YDMF_PHYS%OUT%DIFCS,YDMF_PHYS%OUT%DIFTQ,YDMF_PHYS%OUT%DIFTS,YDMF_PHYS%OUT%FCCQL,YDMF_PHYS%OUT%FCCQN,YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,&
     & YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,YDMF_PHYS%OUT%FPFPCL,YDMF_PHYS%OUT%FPFPCN,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%Q,YDMF_PHYS_STATE%YGSP_RR%T,&
     & ZTENHA,ZTENQVA )
  

  
    
      ! Microphysics occurs in the resolved state.
      ZTMIC(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,KTDIA:YDCPG_DIM%KFLEVG)=YDMF_PHYS_STATE%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,KTDIA:YDCPG_DIM%KFLEVG)
      ZQMIC(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,KTDIA:YDCPG_DIM%KFLEVG)=ZQV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,KTDIA:YDCPG_DIM%KFLEVG)
    

    CALL ACPLUIZ ( YDMODEL%YRML_PHY_MF,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,NTPLUI,YDCPG_DIM%KFLEVG,&
     & ZTMIC, ZQMIC, ZQL, ZQI, ZQR, ZQS,&
     & YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_STATE%YCPG_DYN%RCP%R,&
     & ZNEBS, ZQLIS, ZNEB_CVPP, ZQLI_CVPP,ZQC_DET_PCMT,&
     & ZTENHA, ZTENQVA, .TRUE., YDMF_PHYS_STATE%YCPG_DYN%PHI,&
     & YDMF_PHYS_STATE%YGSP_RR%T, YDMF_PHYS%TMP%APLPAR%FLU%NEIJ, YDMF_PHYS_SURF%GSD_VF%PLSM, PGM, ZVETAF,&
     & YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN,&
     & YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, ZSEDIQL, ZSEDIQI )

    
  

  
  

  ZFLX_LOTT_GWU(:, :) = 0.0_JPRB
  ZFLX_LOTT_GWV(:, :) = 0.0_JPRB

  

!*
!     ------------------------------------------------------------------
!         SAUVEGARDE DES FLUX DE PRECIPITATION CONVECTIVE ET STRATIFORME.

  IF(LNEBN.OR.LRRGUST) THEN
    
      
        DO JLEV=KTDIA-1,YDCPG_DIM%KFLEVG
          DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            YDMF_PHYS%TMP%APLPAR%PFL%FPLCH(JLON,JLEV)=ZFPCOR(JLON,JLEV)
          ENDDO
        ENDDO
      
    
  ENDIF

  IF(LRRGUST) THEN
    DO JLEV=KTDIA-1,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        YDMF_PHYS%TMP%APLPAR%PFL%FPLSH(JLON,JLEV)=YDMF_PHYS%OUT%FPLSL(JLON,JLEV)+YDMF_PHYS%OUT%FPLSN(JLON,JLEV)+YDMF_PHYS%OUT%FPLSG(JLON,JLEV)
      ENDDO
    ENDDO
  ENDIF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! UPDATE TRANSPORT FLUXES DUE TO SEDIMENTATION OF CLOUDS.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  DO JLEV=KTDIA,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)+ZSEDIQL(JLON,JLEV)
      YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)+ZSEDIQI(JLON,JLEV)
      ZFPLSL (JLON,JLEV)=YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)+YDMF_PHYS%OUT%FPLSL (JLON,JLEV)
      ZFPLSN (JLON,JLEV)=YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)+YDMF_PHYS%OUT%FPLSN (JLON,JLEV)
    ENDDO
  ENDDO
  

  ! - - - - - - - - - - - - - - - - -
  ! CORRECT NEGATIVE WATER CONTENTS.
  ! - - - - - - - - - - - - - - - - -

  
    CALL QNGCOR (YDPHY2,YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTPLUI, YDCPG_DIM%KFLEVG,&
     & ZQV, ZQL, ZQI, YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,&
     & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL,&
     & YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN,&
     & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN,&
     & YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FCQNG )

    ! Comment the following lines, which generate negative QR/QS values in CPTEND_NEW.
    !DO JLON=KIDIA,KFDIA
    !  PFPLSL(JLON,KLEV)=PFPLSL(JLON,KLEV)+PDIFTQL(JLON,KLEV)
    !  PFPLSN(JLON,KLEV)=PFPLSN(JLON,KLEV)+PDIFTQI(JLON,KLEV)
    !ENDDO

  

!*
!     ------------------------------------------------------------------
!     12. - BILAN HYDRIQUE DU SOL
!     ---------------------------
  IF ( LMSE ) THEN

    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      PGPAR(JLON,MRAIN)=ZFPLSL(JLON,YDCPG_DIM%KFLEVG)+YDMF_PHYS%OUT%FPLCL(JLON,YDCPG_DIM%KFLEVG)
      PGPAR(JLON,MSNOW)=ZFPLSN(JLON,YDCPG_DIM%KFLEVG)+YDMF_PHYS%OUT%FPLCN(JLON,YDCPG_DIM%KFLEVG)
    ENDDO

  ELSE

    IF ( LSOLV ) THEN
      CALL ACDROV ( YDMODEL%YRML_PHY_MF,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,KCSS,&
       & YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, ZFPLSL, ZFPLSN, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH,&
       & YDMF_PHYS%TMP%APLPAR%DSA%C1, YDMF_PHYS%TMP%APLPAR%DSA%C2, ZC3, ZCN, YDMF_PHYS%OUT%CT, YDMF_PHYS_SURF%GSD_VV%PD2, YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FTR, YDMF_PHYS_SURF%GSD_VV%PLAI, YDMF_PHYS%TMP%APLPAR%FLU%NEIJ,&
       & YDMF_PHYS%TMP%APLPAR%FLU%VEG, ZWFC, ZWPMX, YDMF_PHYS_STATE%YGSP_RR%FC, ZWLMX, ZWSEQ, ZWSMX,&
       & YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FCLL, YDMF_PHYS%OUT%FCLN, YDMF_PHYS%OUT%FCS, YDMF_PHYS%TMP%APLPAR%FLU%FEVI, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_STATE%YGSP_SG%F,&
       & YDMF_PHYS_STATE%YGSP_SB%T, YDMF_PHYS_STATE%YGSP_RR%T, YDMF_PHYS_STATE%YGSP_SB%Q, YDMF_PHYS_STATE%YGSP_SB%TL, YDMF_PHYS_STATE%YGSP_RR%W, YDMF_PHYS_STATE%YGSP_RR%IC,&
       & YDMF_PHYS%OUT%FGEL, YDMF_PHYS%OUT%FGELS,&
       & YDMF_PHYS%OUT%FLWSP, YDMF_PHYS%OUT%FONTE, YDMF_PHYS%OUT%RUISP, YDMF_PHYS%OUT%RUISL, YDMF_PHYS%OUT%RUISS)
    ENDIF

  ENDIF

!*
!-  --------------------------------------------------------------------
!     13.- DRAG MESOSPHERIQUE POUR UN MODELE POSSEDANT DES NIVEAUX
!               AU-DESSUS DE 50 KM (I.E. DANS LA MESOSPHERE)
!     ------------------------------------------------------------------
  
    CALL ACDRME ( YDPHY2,YDTOPH,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,NTDRME,YDCPG_DIM%KFLEVG,&
     & YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDMF_PHYS_STATE%T,ZQV,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,&
     & YDMF_PHYS%OUT%FRMH,YDMF_PHYS%TMP%APLPAR%MSC%FRMQ,YDMF_PHYS%OUT%STRMU,YDMF_PHYS%OUT%STRMV,&
     & RMESOT,RMESOQ,RMESOU,STTEM)
  

!*
!     ------------------------------------------------------------------
!     14.- FLUX PHOTO-CHIMIQUE D'OZONE
!     ------------------------------------------------------------------
  
!*
!     ------------------------------------------------------------------
!     15.- CALCUL DES FLUX D'ENTHALPIE ET DE CHALEUR SENSIBLE LIES AUX
!          PRECIPITATIONS EN FONCTION DES FLUX DE PRECIPITATION
!          ET DE CONDENSATION.
!     ------------------------------------------------------------------

  ! STORE THE PSEUDO-HISTORIC SURFACE PRECIPITATION SENSIBLE HEAT FLUX
  ! -------------------------------------------------------------------

   !LPHSPSH

! Store radiative cloudiness in GFL structure for ISP, Historical files or PostProcessing
  IF (YIRAD%LGP) YDVARS%IRAD%T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) = YDCPG_MISC%QICE(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
  IF (YLRAD%LGP) YDVARS%LRAD%T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) = YDCPG_MISC%QLI (YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
  IF (YA%LGP)    YDVARS%A%T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) = YDCPG_MISC%NEB (YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)


!*
!     -----------------------------------------------------------------------
!     16.- DDH FLEXIBLES POUR LES CHAMPS DE SURFACE VARIABLES/FLUX/TENDANCES.
!     -----------------------------------------------------------------------

  IF(LFLEXDIA) THEN
    IF (LDDH_OMP) THEN
      ! Store in DDH the clouds used for radiative calculations
      ZTMPAF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)*YDCPG_MISC%QICE(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAF,'VIR',YDDDH)
      ZTMPAF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)*YDCPG_MISC%QLI(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAF,'VLR',YDDDH)
      ZTMPAF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)*YDCPG_MISC%NEB(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAF,'VNT',YDDDH)

      ! surface variables
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS_SURF%GSD_VF%PLSM,'VLSM',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS_STATE%YGSP_RR%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)/YDMF_PHYS%OUT%CT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'VSTS',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*RTINER*YDMF_PHYS_STATE%YGSP_SB%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,KCSS)/YDMF_PHYS%OUT%CT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'VSTP',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS_STATE%YGSP_RR%W(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'VSWS',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS_STATE%YGSP_SB%Q(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'VSWP',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS_STATE%YGSP_SG%F(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'VSWN',YDDDH)

      ! surface free-style variables
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%TCLS,'VTCLS',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%QCLS,'VQCLS',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%UCLS,'VUCLS',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%VCLS,'VVCLS',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%NUCLS,'VNUCLS',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%NVCLS,'VNVCLS',YDDDH,CDTYPE='S')
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_STATE%YCPG_DYN%PHI(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'VSPHI',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%GZ0,'VSGZ0',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%GZ0H,'VSGZH',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%ALB,'VSALB',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%CLPH,'VSCPH',YDDDH,CDTYPE='S')

      ! surface fluxes free stlye or not
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FRSO(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSRAYSODW',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FRTH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSRAYTHUP',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FCLL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSCHLATLI',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FCLN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSCHLATNE',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FCS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSCHLATSF',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FCHSP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,KCSS)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSCHGSNPL',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FONTE(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FFONTESLI',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FPLSL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSPRELIGE',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FPLSN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSPRENEGE',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FPLCL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSPRELICO',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FPLCN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSPRENECO',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FEVL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSEVAPLIQ',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FEVN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSEVAPNEG',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FLWSP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSLIQSNPL',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%RUISS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FRUISSURF',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%RUISP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FRUISSPRF',YDDDH)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%FRTHDS,'FSRAYTHDS',YDDDH)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=-YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*RLMLT*YDMF_PHYS%OUT%FONTE(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FFONTESNE',YDDDH)

      ZEPSA = 1.E-4_JPRB
      DO JLON=YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
        ZTMPAS(JLON)=YDMF_PHYS%OUT%FRSO(JLON,YDCPG_DIM%KFLEVG,1)/MAX(1.0_JPRB-YDMF_PHYS%OUT%ALB(JLON),ZEPSA)
      ENDDO
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSRAYSOLR',YDDDH)
      IF (NPRINTLEV >= 2) THEN
      WRITE(NULOUT,*) 'LFLEXDIA ARPEGE WITH NTOTSURF = ',NTOTSURF,&
         & ' AND NTOTSVAR = ',NTOTSVAR, ' AND NTOTSVFS = ',NTOTSVFS
      ENDIF

    ELSE  

      ! Store in DDH the clouds used for radiative calculations
      ZTMPAF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)*YDCPG_MISC%QICE(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
      CALL ADD_FIELD_3D(YDLDDH,ZTMPAF,'VIR','V','ARP',.TRUE.,.TRUE.)
      ZTMPAF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)*YDCPG_MISC%QLI(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
      CALL ADD_FIELD_3D(YDLDDH,ZTMPAF,'VLR','V','ARP',.TRUE.,.TRUE.)
      ZTMPAF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)*YDCPG_MISC%NEB(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
      CALL ADD_FIELD_3D(YDLDDH,ZTMPAF,'VNT','V','ARP',.TRUE.,.TRUE.)

      ! surface variables
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS_SURF%GSD_VF%PLSM,'VLSM','V','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS_STATE%YGSP_RR%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)/YDMF_PHYS%OUT%CT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'VSTS','V','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*RTINER*YDMF_PHYS_STATE%YGSP_SB%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,KCSS)/YDMF_PHYS%OUT%CT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'VSTP','V','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS_STATE%YGSP_RR%W(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'VSWS','V','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS_STATE%YGSP_SB%Q(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'VSWP','V','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS_STATE%YGSP_SG%F(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'VSWN','V','ARP',.TRUE.,.TRUE.)

      ! surface free-style variables
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%TCLS,'VTCLS','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%QCLS,'VQCLS','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%UCLS,'VUCLS','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%VCLS,'VVCLS','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%NUCLS,'VNUCLS','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%NVCLS,'VNVCLS','S','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_STATE%YCPG_DYN%PHI(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'VSPHI','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%GZ0,'VSGZ0','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%GZ0H,'VSGZH','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%ALB,'VSALB','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%CLPH,'VSCPH','S','ARP',.TRUE.,.TRUE.)

      ! surface fluxes free stlye or not
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FRSO(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSRAYSODW','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FRTH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSRAYTHUP','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FCLL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSCHLATLI','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FCLN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSCHLATNE','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FCS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSCHLATSF','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FCHSP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,KCSS)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSCHGSNPL','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FONTE(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FFONTESLI','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FPLSL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSPRELIGE','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FPLSN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSPRENEGE','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FPLCL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSPRELICO','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FPLCN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSPRENECO','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FEVL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSEVAPLIQ','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FEVN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSEVAPNEG','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%FLWSP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSLIQSNPL','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%RUISS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FRUISSURF','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*YDMF_PHYS%OUT%RUISP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FRUISSPRF','F','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%FRTHDS,'FSRAYTHDS','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=-YDMF_PHYS_SURF%GSD_VF%PLSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)*RLMLT*YDMF_PHYS%OUT%FONTE(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FFONTESNE','F','ARP',.TRUE.,.TRUE.)
 
      ZEPSA = 1.E-4_JPRB
      DO JLON=YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
        ZTMPAS(JLON)=YDMF_PHYS%OUT%FRSO(JLON,YDCPG_DIM%KFLEVG,1)/MAX(1.0_JPRB-YDMF_PHYS%OUT%ALB(JLON),ZEPSA)
      ENDDO
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSRAYSOLR','F','ARP',.TRUE.,.TRUE.)
      IF (NPRINTLEV >= 2) THEN
      WRITE(NULOUT,*) 'LFLEXDIA ARPEGE WITH NTOTSURF = ',NTOTSURF,&
         & ' AND NTOTSVAR = ',NTOTSVAR, ' AND NTOTSVFS = ',NTOTSVFS
      ENDIF
    ENDIF
  ENDIF

  IF (LRAFTKE) THEN
     ! DCAPE due to precipitation evaporation.
     CALL ACEVADCAPE(YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,&
      & YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FPLSN,YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_DYN%PHI,ZDCAPE)
     ! Gusts.
     ZQGM=ZEPSNEB
     CALL ACCLDIA(YDXFU,YDPHY,YDPHY2,YDTOPH, YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,&
      & YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, YDMF_PHYS%OUT%CAPE, ZDCAPE, ZTKE1, YDMF_PHYS_STATE%YCPG_DYN%PHIF, POROG, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST, ZBLH, KCLPH)
     YDMF_PHYS%OUT%CLPH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=MIN(XMAXLM,MAX(XMINLM,ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)))
  ENDIF

  
     CALL ACVISIH(YDMODEL%YRML_PHY_MF%YRPHY%YRDVISI,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,KTDIA,YDCPG_DIM%KFLEVG,YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,&
      & YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDCPG_MISC%QLI,YDCPG_MISC%QICE,ZQR,ZQS,ZQGM,YDMF_PHYS%OUT%VISICLD, YDMF_PHYS%OUT%VISIHYD,YDMF_PHYS%OUT%MXCLWC)
  

 !LMPHYS
!----------------------------------------------------
!  CALCUL DE DEPOT HUMIDE POUR LES AEROSOLS DESERTIQUES
!----------------------------------------------------
 ! ENDIF  (LMDUST & NGFL_EXT & LRDEPOS)



!     ------------------------------------------------------------------

! BAYRAD
! Compute convective hydrometeors mixing ratio from diagnistic fluxes
!--------------------------------------------------------------------


   ! Convert from flux [kg/m2/s] to density [kg/m3] using old RTTOV-SCATT
   !  a    b         ! RR = a * LWC^b, [RR]=mm/h, [LWC]=g/m^3
   ! 20.89 1.15      ! rain
   ! 29.51 1.10      ! snow
   ZRHORAIN     = 1.0_JPRB
   ZRHOSNOW     = 0.1_JPRB

   ZCOEFRAIN(1) = ((1.0_JPRB / 20.89_JPRB) * 3600.0_JPRB) / ZRHORAIN
   ZCOEFRAIN(2) = (1.0_JPRB / 1.15_JPRB)
   ZCOEFSNOW(1) = ((1.0_JPRB / 29.51_JPRB) * 3600.0_JPRB) / ZRHOSNOW
   ZCOEFSNOW(2) = (1.0_JPRB / 1.10_JPRB)

   YDMF_PHYS%TMP%APLPAR%BAY%QRCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) =  0.0_JPRB
   YDMF_PHYS%TMP%APLPAR%BAY%QSCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) =  0.0_JPRB


   DO JLEV=1,YDCPG_DIM%KFLEVG
     DO JLON= YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
       YDMF_PHYS%TMP%APLPAR%BAY%QRCONV(JLON,JLEV)  = 0.001_JPRB * ( (ABS(YDMF_PHYS%OUT%FPLCL(JLON,JLEV)) * ZCOEFRAIN(1)) ** ZCOEFRAIN(2) )
       YDMF_PHYS%TMP%APLPAR%BAY%QSCONV(JLON,JLEV)  = 0.001_JPRB * ( (ABS(YDMF_PHYS%OUT%FPLCN(JLON,JLEV)) * ZCOEFSNOW(1)) ** ZCOEFSNOW(2) )    
     ENDDO
   ENDDO


   ! Convert density [kg/m3] to mixing ratio [kg/kg]
   ! R_dry (dry air constant)

   ZDE2MR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)  = RD * YDMF_PHYS_STATE%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) / YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
   YDMF_PHYS%TMP%APLPAR%BAY%QRCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) = YDMF_PHYS%TMP%APLPAR%BAY%QRCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) * ZDE2MR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
   YDMF_PHYS%TMP%APLPAR%BAY%QSCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) = YDMF_PHYS%TMP%APLPAR%BAY%QSCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) * ZDE2MR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)




! Precipitation Type

! Compute wet-bulb temperature at 2 meters (suppose homogeneity of qv/ql/qi )
!ZPCLS(KIDIA:KFDIA)=PAPRS(KIDIA:KFDIA,KLEV)-2._JPRB/ZZZF(KIDIA:KFDIA,1,KLEV)*&
!                 &(PAPRS(KIDIA:KFDIA,KLEV)-PAPRSF(KIDIA:KFDIA,KLEV))
CALL PPWETPOINT(YDPHY,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,ZPCLS,YDMF_PHYS%OUT%TCLS,&
  & YDMF_PHYS_STATE%Q(:,YDCPG_DIM%KFLEVG),YDMF_PHYS_STATE%L(:,YDCPG_DIM%KFLEVG),YDMF_PHYS_STATE%I(:,YDCPG_DIM%KFLEVG),YDMF_PHYS%OUT%TPWCLS)



   ! Defined precipitation type
   !
   ! Compute wet-bulb temperature
  DO JLEV=1,YDCPG_DIM%KFLEVG
    CALL PPWETPOINT(YDPHY,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(:,JLEV),YDMF_PHYS_STATE%T(:,JLEV),&
     & YDMF_PHYS_STATE%Q(:,JLEV),YDMF_PHYS_STATE%L(:,JLEV),YDMF_PHYS_STATE%I(:,JLEV),ZTW(:,JLEV))
  ENDDO

  DO JLON=1,YDCPG_DIM%KLON
      ZFPLS(JLON,YDCPG_DIM%KFLEVG)=YDMF_PHYS%OUT%FPLCN(JLON,YDCPG_DIM%KFLEVG)+YDMF_PHYS%OUT%FPLSN(JLON,YDCPG_DIM%KFLEVG)
      ZFPLC(JLON,YDCPG_DIM%KFLEVG)=YDMF_PHYS%OUT%FPLCL(JLON,YDCPG_DIM%KFLEVG)+YDMF_PHYS%OUT%FPLSL(JLON,YDCPG_DIM%KFLEVG)
      ZFPLSG(JLON,YDCPG_DIM%KFLEVG)=0._JPRB
  ENDDO

  !initialisation de ZZZ
  DO JLEV = 1,YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZZZ(JLON,1,JLEV)=YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO

  !initialisation de ZDZZ
  DO JLEV = 2, YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZDZZ(JLON,1,JLEV)=ZZZ(JLON,1,JLEV-1)-ZZZ(JLON,1,JLEV)
    ENDDO
  ENDDO
  DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZDZZ(JLON,1,1)=YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,0)*ZINVG-ZZZ(JLON,1,1)
  ENDDO

  

   NDTPRECCUR=INT(MOD(ZSTATI/TSTEP,REAL(NDTPREC)))+1_JPIM

   CALL DPRECIPS (YDPHY%YRDPRECIPS,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,POROG,YDMF_PHYS%OUT%TPWCLS,YDMF_PHYS%OUT%DIAGH,YDMF_PHYS_STATE%YCPG_DYN%PHIF,&
      & ZDZZ,ZTW,YDMF_PHYS_STATE%L,ZFPLC(:,YDCPG_DIM%KFLEVG),ZFPLS(:,YDCPG_DIM%KFLEVG),ZFPLSG(:,YDCPG_DIM%KFLEVG),YDMF_PHYS%TMP%PRC%DPRECIPS(:,NDTPRECCUR))

  ! WRITE(20,*)'sous aplpar ZDPRECIPS ', ZDPRECIPS(KIDIA:KFDIA,NDTPRECCUR),NDTPRECCUR

  

  

 !Idem for an other time step and an other period
   NDTPRECCUR2=INT(MOD(ZSTATI/TSTEP,REAL(NDTPREC2)))+1_JPIM

   CALL DPRECIPS(YDPHY%YRDPRECIPS,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,POROG,YDMF_PHYS%OUT%TPWCLS,YDMF_PHYS%OUT%DIAGH,YDMF_PHYS_STATE%YCPG_DYN%PHIF,&
      & ZDZZ,ZTW,YDMF_PHYS_STATE%L,ZFPLC(:,YDCPG_DIM%KFLEVG),ZFPLS(:,YDCPG_DIM%KFLEVG),ZFPLSG(:,YDCPG_DIM%KFLEVG),YDMF_PHYS%TMP%PRC%DPRECIPS2(:,NDTPRECCUR2))

  ! WRITE(20,*)'sous aplpar ZDPRECIPS2',ZDPRECIPS2(KIDIA:KFDIA,NDTPRECCUR2),NDTPRECCUR2

  
  

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('APLPAR',1,ZHOOK_HANDLE)
END SUBROUTINE APLPAR_NEW
