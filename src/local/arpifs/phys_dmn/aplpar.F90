!OPTIONS XOPT(NOEVAL)
SUBROUTINE APLPAR(YDGEOMETRY,YDMF_PHYS_STATE,YDCPG_DYN0,YDCPG_MISC,YDMF_PHYS,YDMF_PHYS_SURF,YDVARS,YDSURF, YDXFU,YDCFU,YDMODEL,KIDIA , KFDIA  , KLON      ,&
 & KTDIA  , KLEV   , KSTGLO ,&
 & KVCLIS , KVCLIV , KSTEP  ,&
 & KSGST  , KCSS   ,&
 & KBL    , KGPCOMP, KNFRRC , PDT,&
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
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KBL
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPCOMP
INTEGER(KIND=JPIM),INTENT(IN)    :: KNFRRC
INTEGER(KIND=JPIM),INTENT(IN)    :: KVCLIS
INTEGER(KIND=JPIM),INTENT(IN)    :: KSGST
INTEGER(KIND=JPIM),INTENT(IN)    :: KCSS
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTGLO
INTEGER(KIND=JPIM),INTENT(IN)    :: KVCLIV
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT
LOGICAL           ,INTENT(IN)    :: LDXFUMSE
REAL(KIND=JPRB)   ,INTENT(IN)    :: PINDX(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PINDY(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCORI(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKOZO(KLON,KLEV,KVCLIS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGPAR(KLON,YDMODEL%YRML_PHY_MF%YRPARAR%NGPAR+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POMPAC(KLON,0:KLEV,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAC(KLEV+1,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDIFSV(KLON,0:KLEV,1:YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)

! Prognostic convection variables (MUST BE 1:KLEV for GFL)
! Unless mf_phys has copied them into a 0:KLEV array
! Up to now, we leave 1:KLEV for pentch: see later.
!---------------------------------------------------

INTEGER(KIND=JPIM),INTENT(OUT)   :: KCLPH(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFTCNS(KLON,0:KLEV,6)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEZDIAG(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDPTKE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENDEXT_DEP(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTDISS(KLON,KLEV)
TYPE (TRAJ_PHYS_TYPE), INTENT(INOUT) :: PTRAJ_PHYS
TYPE(TYP_DDH), INTENT(INOUT)     :: YDDDH
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGP2DSPP(KLON,YSPP%N2D)

!     ------------------------------------------------------------------
!        ATTENTION SI KVCLIG < 7 LES CHAMPS SUIVANTS NE SONT
!        PAS REELLEMENT ALLOUES EN MEMOIRE.
!*
!     ------------------------------------------------------------------
!     DECLARATION DES TABLEAUX LOCAUX-GLOBAUX DES PARAMETRISATIONS
INTEGER(KIND=JPIM) :: INLAB(KLON,KLEV), INLAB_CVPP(KLON,KLEV)
REAL(KIND=JPRB) :: ZVETAH(0:KLEV),ZNLAB(KLON,KLEV),ZNLABCVP(KLON,KLEV)
REAL(KIND=JPRB) :: ZVETAF(KLEV)

!     ------------------------------------------------------------------
!     ARE DIMENSIONNED 0:KLEV ONLY IN ORDER TO KEEP IN MIND
!     THAT THEY ARE COMPUTED AT "HALF LEVELS".
!     THEY ARE USED HOWEVER FROM 1 TO KLEV.

REAL(KIND=JPRB) :: ZXTROV(KLON,0:KLEV),ZXUROV(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZXPTKEROV(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZRRCOR(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZMRIPP(KLON,0:KLEV),ZMRIFPP(KLON,0:KLEV),ZBNEBCVPP(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZBNEBQ(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZMN2PP(KLON,0:KLEV),ZMN2_ES(KLON,0:KLEV),ZMN2_EQ(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZMN2_DS(KLON,0:KLEV),ZMN2_DQ(KLON,0:KLEV)
! ZMRIMC     : M(C) from P. Marquet's moist Ri computation - for TKE correction after TOMs
! ZMRICTERM  : Rv/R.F(C)-1/M(C).T/Tv from P. Marquet's moist Ri computation - for TKE correction after TOMs
REAL(KIND=JPRB) :: ZMRIMC(KLON,0:KLEV),ZMRICTERM(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZTSTAR(KLON,KLEV),ZTSTAR2(KLON,0:KLEV) !diag 
REAL(KIND=JPRB) :: ZTSTARQ(KLON,KLEV),ZTSTAR2Q(KLON,0:KLEV)!diag
REAL(KIND=JPRB) :: ZTAU_TKE(KLON,0:KLEV)!DISSIPATION TIME SCALE TAU  -FOR TOM's CALCULATION
REAL(KIND=JPRB) :: ZF_EPS(KLON,0:KLEV)   !  Conversion function lm-L
REAL(KIND=JPRB) :: ZFUN_TTE(KLON,0:KLEV)   !  Function in computation of tte_tilde
REAL(KIND=JPRB) :: ZKTROV(KLON,0:KLEV),ZKUROV(KLON,0:KLEV),ZNBVNO(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZKQROV(KLON,0:KLEV),ZKQLROV(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZKNROV(KLON,0:KLEV)

REAL(KIND=JPRB) :: ZFHORM(KLON,0:KLEV),ZFHORH(KLON,0:KLEV) ! arrays for 3D turb
! ZFMTKE - F_m function for static K_m computation
! ZFHTKE - F_h function for static K_h computation
REAL(KIND=JPRB) :: ZFMTKE(KLON,0:KLEV),ZFTTKE(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZRHS(KLON,KLEV)
! ZAUTKE - alpha_u for dry AF scheme
! ZATTKE - alpha_theta for dry AF scheme
REAL(KIND=JPRB) :: ZAUTKE(KLON,0:KLEV),ZATTKE(KLON,0:KLEV)
!ZTH_FUN, ZWW_FUN - T_h, A_h, F_ww - stability functions for TOMs par.
REAL(KIND=JPRB) :: ZTH_FUN(KLON,0:KLEV),ZWW_FUN(KLON,0:KLEV)
! ZFMGST - stability function F_m for moist gustiness correction
REAL(KIND=JPRB) :: ZFMGST(KLON,0:KLEV)


!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZATSLC(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZNEBS(KLON,KLEV),ZQLIS(KLON,KLEV)
REAL(KIND=JPRB) :: ZNEBS0(KLON,KLEV),ZQLIS0(KLON,KLEV)
REAL(KIND=JPRB) :: ZNEBC0(KLON,KLEV)   !Nebulosite convective radiative
REAL(KIND=JPRB) :: ZNEBDIFF(KLON,KLEV) !Nebulosite: calcul de la diffusion
REAL(KIND=JPRB) :: ZNEBCH(KLON,KLEV)   !Nebulosite convective condensation
REAL(KIND=JPRB) :: ZUNEBH(KLON,KLEV)   !Nebulosite convective histo
REAL(KIND=JPRB) :: ZDETFI(KLON,KLEV)   !fraction of instantaneous detrained air
REAL(KIND=JPRB) :: ZFPCOR(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZFHP(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZLMT(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZZLMT(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZLMU(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZLMU2(KLON,0:KLEV),ZLMT2(KLON,0:KLEV) ! temporary storage of lm,lh
REAL(KIND=JPRB) :: ZLML(KLON,0:KLEV) ! TKE type mixing length
REAL(KIND=JPRB) :: ZLMLTILD(KLON,0:KLEV) ! 'STATIC' TKE type mixing length
REAL(KIND=JPRB) :: ZOME(KLON,KLEV)   ! updraught envt vert vel*dt
REAL(KIND=JPRB) :: ZFALLR(KLON,KLEV) ! fall velocity of rain
REAL(KIND=JPRB) :: ZFALLS(KLON,KLEV) ! fall velocity of snow
REAL(KIND=JPRB) :: ZFALLG(KLON,KLEV) ! fall velocity of graupel
REAL(KIND=JPRB) :: ZICEFR1(KLON,KLEV)! Resolved Condensate ice fraction
REAL(KIND=JPRB) :: ZRHCRI(KLON,KLEV) ! Smith scheme critical RH
REAL(KIND=JPRB) :: ZRHDFDA(KLON,KLEV)! RK scheme change in RH over cloud
REAL(KIND=JPRB) :: ZLHS(KLON,KLEV)   ! Sublimation latent heat
REAL(KIND=JPRB) :: ZLHV(KLON,KLEV)   ! Evaporation latent heat
REAL(KIND=JPRB) :: ZLH(KLON,KLEV)    ! Temporar storage for updated PLH
REAL(KIND=JPRB) :: ZLSCPE(KLON,KLEV) ! Temporar storage for updated PLSCPE
REAL(KIND=JPRB) :: ZQSAT(KLON,KLEV)  ! Temporar storage for updated PQSAT
REAL(KIND=JPRB) :: ZQSATS(KLON,KLEV) ! QSAT of resolved cond./evap. scheme
REAL(KIND=JPRB) :: ZQW(KLON,KLEV)    ! Temporar storage for updated PQW
REAL(KIND=JPRB) :: ZRH(KLON,KLEV)    ! Temporar storage for updated PRH
REAL(KIND=JPRB) :: ZTW(KLON,KLEV)    ! Temporar storage for updated PTW)
REAL(KIND=JPRB) :: ZDQ(KLON,KLEV)    ! Saturation departure for a given thermodynamic state
REAL(KIND=JPRB) :: ZDQM(KLON,KLEV)   ! maximum saturation departure
REAL(KIND=JPRB) :: ZPOID(KLON,KLEV)  ! DP/(RG*DT) FOR A GIVEN LEVEL AND A GIVEN TIME STEP.
REAL(KIND=JPRB) :: ZIPOI(KLON,KLEV)  ! INVERSE OF ZPOID.

REAL(KIND=JPRB) :: ZQU(KLON,KLEV) ! Updraught Specific moisture
REAL(KIND=JPRB) :: ZTU(KLON,KLEV) ! Updraught Temperature
REAL(KIND=JPRB) :: ZUU(KLON,KLEV) ! Updraught zonal wind
REAL(KIND=JPRB) :: ZVU(KLON,KLEV) ! Updraught merid. wind

REAL(KIND=JPRB) :: ZTMIC(KLON,KLEV) ! Temperature for microphysics
REAL(KIND=JPRB) :: ZQMIC(KLON,KLEV) ! Specific moisture for microphysics

REAL(KIND=JPRB) :: ZT(KLON,KLEV)     ! updated temperature T for cascading parameterization
REAL(KIND=JPRB) :: ZTCORR(KLON,KLEV) ! temperature corr. for convective cloud
REAL(KIND=JPRB) :: ZMELNET(KLON,KLEV) ! net melting (-freezing) rate of ice
REAL(KIND=JPRB) :: ZMELGET(KLON,KLEV)! net melting (-freezing) rate of graupel
REAL(KIND=JPRB) :: ZU(KLON,KLEV)     ! updated zonal velocity
REAL(KIND=JPRB) :: ZV(KLON,KLEV)     ! updated meridional velocity

REAL(KIND=JPRB) :: ZQV(KLON,KLEV)    ! corrected (for negative values) vapour
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZQI(KLON,KLEV)    ! corrected (for negative values) cloud ice
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZQL(KLON,KLEV)    ! corrected (for negative values) cloud liquid
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZQR(KLON,KLEV)    ! corrected (for negative values) rain
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZQS(KLON,KLEV)    ! corrected (for negative values) snow
                                     ! updated value for cascading parameterization
REAL(KIND=JPRB) :: ZCP(KLON,KLEV)    ! new cp for turbulent diffusion
REAL(KIND=JPRB) :: ZQT(KLON,KLEV)
REAL(KIND=JPRB) :: ZTENHA(KLON,KLEV)
REAL(KIND=JPRB) :: ZTENQVA(KLON,KLEV)
REAL(KIND=JPRB) :: ZFCQVNG(KLON,0:KLEV) ! correction flux increment for neg vapour
REAL(KIND=JPRB) :: ZFCQING(KLON,0:KLEV) ! correction flux increment for neg ice
REAL(KIND=JPRB) :: ZFCQLNG(KLON,0:KLEV) ! correction flux increment for neg liquid water
REAL(KIND=JPRB) :: ZFPLSL(KLON,0:KLEV) ! total liquid water flux: diff+sedi+rain
REAL(KIND=JPRB) :: ZFPLSN(KLON,0:KLEV) ! total solid water flux: diff+sedi+snow
REAL(KIND=JPRB) :: ZFCQL(KLON,0:KLEV)   ! condensation flux(liquid)
REAL(KIND=JPRB) :: ZFCQI(KLON,0:KLEV)   ! condensation flux(ice)
REAL(KIND=JPRB) :: ZDIFCQD(KLON,0:KLEV) ! downdraft flux of specific humidity
REAL(KIND=JPRB) :: ZDIFCQLD(KLON,0:KLEV) ! downdraft flux of liquid water
REAL(KIND=JPRB) :: ZDIFCQID(KLON,0:KLEV) ! downdraft flux of  solid water
REAL(KIND=JPRB) :: ZSEDIQL(KLON,0:KLEV) ! sedimentation flux of cloud liquid water
REAL(KIND=JPRB) :: ZSEDIQI(KLON,0:KLEV) ! sedimentation flux of cloud ice water
REAL(KIND=JPRB) :: ZDIFCSD(KLON,0:KLEV) ! downdraft entalphy flux
REAL(KIND=JPRB) :: ZSTRCUD(KLON,0:KLEV) ! change in horizontal mom.
REAL(KIND=JPRB) :: ZSTRCVD(KLON,0:KLEV) ! change in horizontal mom.
REAL(KIND=JPRB) :: ZRCVOTT(KLON,0:KLEV) ! degree of inhomogeneity in precips.
REAL(KIND=JPRB) :: ZSIGPC(KLON)         ! Convective precipit mesh fraction
REAL(KIND=JPRB) :: ZSIGP(KLON)         ! Precipitation mesh fraction
REAL(KIND=JPRB) :: ZAUXPRC(KLON)        ! Precipitation auxilary
REAL(KIND=JPRB) :: ZDIFCVPPQ(KLON,0:KLEV) ! Flux de CVPP (KFB or EDKF) sur Qv
REAL(KIND=JPRB) :: ZDIFCVPPS(KLON,0:KLEV) ! Flux de CVPP (KFB or EDKF) sur CpT
REAL(KIND=JPRB) :: ZDIFCVTH(KLON,0:KLEV) ! Flux de CV sur Theta air sec
REAL(KIND=JPRB) :: ZDIFCVPPU(KLON,0:KLEV) ! Flux de CVPP (EDKF) sur U
REAL(KIND=JPRB) :: ZDIFCVPPV(KLON,0:KLEV) ! Flux de CVPP (EDKF) sur V

REAL(KIND=JPRB) :: ZEDMFQ(KLON,0:KLEV) ! Mass flux part of EDMF flux for Qv
REAL(KIND=JPRB) :: ZEDMFS(KLON,0:KLEV) ! Mass flux part of EDMF flux for CpT
REAL(KIND=JPRB) :: ZEDMFU(KLON,0:KLEV) ! Mass flux part of EDMF flux for U
REAL(KIND=JPRB) :: ZEDMFV(KLON,0:KLEV) ! Mass flux part of EDMF flux for V
REAL(KIND=JPRB) :: ZMF_UP(KLON,0:KLEV) ! Mass flux for implicit formulation of EDMF equation (LEDMFI)
REAL(KIND=JPRB) :: ZMU(KLON,0:KLEV) ! Flux de masse (updraft) pour XIOS output
REAL(KIND=JPRB) :: ZMD(KLON,0:KLEV) ! Flux de masse (downdraft) pour XIOS output

REAL(KIND=JPRB) :: ZCONDCVPPL(KLON,0:KLEV) ! Flux de condensation liquide du a CVVPP (KFB)
REAL(KIND=JPRB) :: ZCONDCVPPI(KLON,0:KLEV) ! Flux de condensation glace du a CVVPP (KFB)
REAL(KIND=JPRB) :: ZPRODTH_CVPP(KLON,0:KLEV) ! Flux de production thermique de TKE du a CVPP(KFB)
REAL(KIND=JPRB) :: ZDTRAD(KLON,KLEV) ! radiation contribution to T tendency
REAL(KIND=JPRB) :: ZDQVDIFF(KLON,KLEV) ! turtb.diff contribution to Qv tendency
REAL(KIND=JPRB) :: ZRKQCTEND(KLON,KLEV) ! Qc input for RK condensation scheme
REAL(KIND=JPRB) :: ZRKQVTEND(KLON,KLEV) ! Qv input for RK condensation scheme
REAL(KIND=JPRB) :: ZRKTTEND(KLON,KLEV) ! T input for RK condensation scheme
REAL(KIND=JPRB) :: ZDQV, ZDQI, ZDQL, ZDQR, ZDQS, ZDQC, ZGDT, ZGDTI,&
                 & ZQV0, ZQX0, ZQX1,&
                 & ZCONVC, ZTOTC,ZDTURDIFF
REAL(KIND=JPRB) :: ZTMPAF(KLON,KLEV)    ! temporary array for Add_Field_3d.
REAL(KIND=JPRB) :: ZTMPAS(KLON)         ! temporary array for Add_Field_2d..
REAL(KIND=JPRB) :: ZTMPPRODTH(KLON,0:KLEV) ! temporary array

!!for BAYRAD
REAL(KIND=JPRB) :: ZDE2MR(KLON,KLEV) ! temporary array for conversion of density to mixing ratio
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

INTEGER(KIND=JPIM) :: INND(KLON)
REAL(KIND=JPRB) :: ZXDROV(KLON),ZXHROV(KLON)
REAL(KIND=JPRB) :: ZUGST(KLON),ZVGST(KLON)
REAL(KIND=JPRB) :: ZCDROV(KLON),ZCHROV(KLON),ZDQSTS(KLON),ZGWDCS(KLON),&
 & ZHQ(KLON),ZHU(KLON),ZHTR(KLON),ZCDNH(KLON),ZMOD(KLON),&
 & ZRTI(KLON),ZDPHI(KLON),ZPRS(KLON),ZSTAB(KLON),ZTAUX(KLON)
REAL(KIND=JPRB) :: ZWFC(KLON),ZWPMX(KLON),ZWLMX(KLON),ZWSEQ(KLON),&
 & ZWSMX(KLON),ZWWILT(KLON),&
 & ZC3(KLON),ZCG(KLON),ZCN(KLON),&
 & ZNEIJG(KLON),ZNEIJV(KLON)
REAL(KIND=JPRB) :: ZPCLS(KLON)
REAL(KIND=JPRB) :: ZPREN(KLON)
REAL(KIND=JPRB) :: ZFRSODS(KLON)
REAL(KIND=JPRB) :: ZCD(KLON)

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
REAL(KIND=JPRB) :: ZDAER(KLEV), ZHUC(KLEV), ZBLH(KLON)
REAL(KIND=JPRB) :: ZQO3(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZAER(KLON,KLEV,6)
REAL(KIND=JPRB) :: ZAERD(KLON,KLEV)
REAL(KIND=JPRB) :: ZAERINDS(KLON,KLEV)
REAL(KIND=JPRB) :: ZQCO2(KLON,KLEV)
REAL(KIND=JPRB) :: ZROZ(KLON,KLEV)
REAL(KIND=JPRB) :: ZDPHIV(KLON),ZDPHIT(KLON)

REAL(KIND=JPRB) :: ZMAN(0:KLEV), ZMAK(0:KLEV)
REAL(KIND=JPRB) :: ZSSO_STDEV(KLON)
REAL(KIND=JPRB) :: ZTWSNOW(KLON),ZTOWNS(KLON)
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

REAL(KIND=JPRB) :: ZTRSOD(KLON)

!            2-D ARRAYS
!            ----------

!* OUTPUT ARGUMENTS FOR THE ECMWF PHYSICS

!            0.2  LOCAL ARRAYS FOR ECMWF PHYSICS PACKAGE
!                 --------------------------------------

REAL(KIND=JPRB) :: ZCEMTR(KLON,0:1) , ZCTRSO(KLON,0:1)
REAL(KIND=JPRB) :: ZALBD(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW), ZALBP(KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW),&
                 & ZALB(KLON)
REAL(KIND=JPRB) :: ZSFSWDIR (KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW), ZSFSWDIF (KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB) :: ZTRSODIR (KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW), ZTRSODIF (KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB) :: ZFSDNN(KLON),ZFSDNV(KLON)

!            1-D ARRAYS
!            ----------

REAL(KIND=JPRB) :: ZSUDU(KLON) , ZDSRP(KLON) , ZSDUR(KLON)
REAL(KIND=JPRB) :: ZTHETAVS(KLON), ZTHETAS(KLON)
REAL(KIND=JPRB) :: ZAESEA(KLON), ZAELAN(KLON), ZAESOO(KLON), ZAEDES(KLON)
REAL(KIND=JPRB) :: ZAESUL(KLON), ZAEVOL(KLON)
REAL(KIND=JPRB) :: ZMERL(KLON)

!            2-D ARRAYS
!            ----------

REAL(KIND=JPRB) :: ZTENT(KLON,KLEV), ZGEOSLC(KLON,KLEV)
REAL(KIND=JPRB) :: ZTHETAV(KLON,KLEV)

!            LOCAL ARRAYS FOR TKE
! ZCOEFN : COEFFICIENT STATISTIQUE POUR LES FLUX D'EAUX CONDENSEES.

REAL(KIND=JPRB) :: ZCOEFN(KLON,KLEV)

!            LOCAL ARRAYS FOR ACVPPKF
REAL(KIND=JPRB) :: ZQLI_CVPP(KLON,KLEV)
REAL(KIND=JPRB) :: ZNEB_CVPP(KLON,KLEV)
!            LOCAL ARRAYS FOR EDKF
REAL(KIND=JPRB) :: ZIMPL

!           2-D ARRAY FOR SIMPL.RADIATION SCHEME

REAL(KIND=JPRB) :: ZZNEB(KLON,KLEV)

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

REAL(KIND=JPRB) :: ZNEBC(KLON,KLEV),ZQLIC(KLON,KLEV)
REAL(KIND=JPRB) :: ZQCL(KLON,KLEV),ZQCI(KLON,KLEV)
REAL(KIND=JPRB) :: ZFHEVPPC(KLON,0:KLEV),ZFHMLTSC(KLON,0:KLEV)&
 & ,ZFPEVPPC(KLON,0:KLEV)
INTEGER(KIND=JPIM) :: ICIS(KLON,KLEV)


!           SURFEX local VARIABLES
!           ----------------------------

! Implicit coupling coefficients
INTEGER(KIND=JPIM) :: IRR
REAL(KIND=JPRB) :: ZDTMSE,ZRHGMT,ZSTATI
REAL(KIND=JPRB) :: ZCFAQ(KLON,KLEV),ZCFAS(KLON,KLEV),&
 & ZCFATH(KLON,KLEV),ZCFAU(KLON,KLEV),ZCFBQ(KLON,KLEV),&
 & ZCFBS(KLON,KLEV),ZCFBTH(KLON,KLEV),&
 & ZCFBU(KLON,KLEV),ZCFBV(KLON,KLEV),&
 & ZCFBU_G(KLON,KLEV),ZCFBV_G(KLON,KLEV),&
 & ZCFBS_G(KLON,KLEV),ZCFBQ_G(KLON,KLEV)

REAL(KIND=JPRB) :: ZDSE(KLON,KLEV)
REAL(KIND=JPRB) :: ZFEV(KLON),ZFMDU(KLON),ZFMDV(KLON),ZFEVS(KLON)
REAL(KIND=JPRB) :: ZSRAIN(KLON), ZSSNOW(KLON), ZSGROUPEL(KLON)
REAL(KIND=JPRB) :: ZSFCO2(KLON), ZRHODREFM(KLON)
REAL(KIND=JPRB) :: ZDEPTH_HEIGHT(KLON,KLEV), ZZS(KLON)
REAL(KIND=JPRB) :: ZTSN(KLON),ZTN(KLON)
REAL(KIND=JPRB)   :: ZBUDTH (KLON)
REAL(KIND=JPRB)   :: ZBUDSO (KLON)
REAL(KIND=JPRB)   :: ZFCLL  (KLON)
! FOR Hv
REAL(KIND=JPRB)   :: ZHV2(KLON)
! FOR DUST
REAL(KIND=JPRB), DIMENSION (KLON,KLEV) ::  ZQDM
REAL(KIND=JPRB) :: ZCFASV(KLON,KLEV,1:YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) ! SCOND MEMBRE POUR LES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZCFBSV(KLON,KLEV,1:YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) ! SCOND MEMBRE POUR LES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZSMOOTRAC(1:YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) ! SCOND MEMBRE POUR LES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZINVDT, ZINVG, ZRSCP, ZINVATM

REAL(KIND=JPRB),  DIMENSION(KLON,1,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT):: ZZI_SVM
REAL(KIND=JPRB),  DIMENSION (KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG):: ZZI_PEZDIAG
REAL(KIND=JPRB), DIMENSION (:,:,:), ALLOCATABLE  :: ZSVM, ZPSV
REAL(KIND=JPRB), DIMENSION (:,:),   ALLOCATABLE  :: ZSFSV  ! passifs scalaires surf flux
! TRAITEMENT DES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZTM(KLON,KLEV)
REAL(KIND=JPRB), DIMENSION (KLON,1,KLEV) :: ZZZ,ZDZZ,ZZI_PABSM, ZZI_THM,&
                & ZZI_EXNREFM, ZZI_RHODREFM,ZEVAP,ZZDEP,ZZI_RHO
REAL(KIND=JPRB) :: ZZI_APHI(KLON,0:KLEV)
REAL(KIND=JPRB), DIMENSION (KLON,1,KLEV,6) :: ZZI_RM
!            3-D ARRAYS
REAL(KIND=JPRB), DIMENSION(KLON,KLEV,YDMODEL%YRML_PHY_RAD%YRERAD%NSW):: ZPIZA_DST !Single scattering
                                             ! albedo of dust (points,lev,wvl)
REAL(KIND=JPRB), DIMENSION(KLON,KLEV,YDMODEL%YRML_PHY_RAD%YRERAD%NSW):: ZCGA_DST  !Assymetry factor
                                             ! for dust (points,lev,wvl)
REAL(KIND=JPRB), DIMENSION(KLON,KLEV,YDMODEL%YRML_PHY_RAD%YRERAD%NSW):: ZTAUREL_DST !tau/tau_{550}
                                             !dust (points,lev,wvl)

!         ACFLUSO (ECUME) local variable
!-------------------------------------------
REAL(KIND=JPRB) :: ZCE(KLON), ZCEROV(KLON), ZCRTI(KLON)


!        New ACDIFV1 local variable
!--------------------------------------------
REAL(KIND=JPRB)   :: ZXURO(KLON,0:KLEV)
REAL(KIND=JPRB)   :: ZXQRO(KLON,0:KLEV)
REAL(KIND=JPRB)   :: ZXTRO(KLON,0:KLEV)

!        New ACNORGWD local variables
!--------------------------------------------
REAL(KIND=JPRB) :: ZD_U(KLON,KLEV), ZD_V(KLON,KLEV)
REAL(KIND=JPRB) :: Z_PP(KLON,KLEV) 
REAL(KIND=JPRB) :: Z_UU(KLON,KLEV) 
REAL(KIND=JPRB) :: Z_VV(KLON,KLEV)
REAL(KIND=JPRB) :: Z_TT(KLON,KLEV)
REAL(KIND=JPRB) :: Z_VO(KLON,KLEV)  
REAL(KIND=JPRB) :: ZFLX_LOTT_GWU(KLON,0:KLEV), ZFLX_LOTT_GWV(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZPRECGWD(KLON)

!    TKE+ for ACCLDIA
REAL(KIND=JPRB)   :: ZTKE1(KLON,KLEV)
REAL(KIND=JPRB)   :: ZTPRDY(KLON,KLEV)

!   For ACVISIH
REAL(KIND=JPRB)   :: ZQGM(KLON,KLEV)


!        New ARP_GROUND_PARAM local variable
!------------------------------------------------

REAL(KIND=JPRB)   :: ZALPHA1(KLON,0:KLEV)
REAL(KIND=JPRB)   :: ZCOEFA (KLON,0:KLEV)
REAL(KIND=JPRB)   :: ZLVT   (KLON,0:KLEV)
REAL(KIND=JPRB)   :: ZQICE  (KLON,0:KLEV)
REAL(KIND=JPRB)   :: ZDIFWQ (KLON)
REAL(KIND=JPRB)   :: ZDIFWS (KLON)
REAL(KIND=JPRB)   :: ZSC_FEVI (KLON),ZSC_FEVN(KLON),ZSC_FCLL(KLON),ZSC_FCLN(KLON)

!           TRAJECTORY (For diffusion !) local VARIABLES
!           ----------------------------
REAL(KIND=JPRB) :: ZCDROV_SAVE(KLON),ZCHROV_SAVE(KLON)
REAL(KIND=JPRB) :: ZKTROV_SAVE(KLON,0:KLEV),ZKUROV_SAVE(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZTRAJGWD(KLON,0:KLEV) !Traj buffer saved for TL/AD (YDMODEL%YRML_PHY_MF%YRSIMPHL%LGWDSPNL)

REAL(KIND=JPRB)    :: ZRVMD,ZDELTA
LOGICAL            :: LLAERO, LLAROME, LLCALLRAD
REAL(KIND=JPRB)    :: ZAIPCMT(KLON) ! Activity Index of PCMT: 1. if PCMT is active, 0. else case.
REAL(KIND=JPRB)    :: ZALF_CAPE(KLON)
REAL(KIND=JPRB)    :: ZALF_CVGQ(KLON)
REAL(KIND=JPRB)    :: ZQIC(KLON,KLEV),ZQLC(KLON,KLEV),ZQRC(KLON,KLEV),ZQSC(KLON,KLEV),ZQVI
REAL(KIND=JPRB)    :: ZQLI_CVP(KLON,KLEV)
REAL(KIND=JPRB)    :: ZQC_DET_PCMT(KLON,KLEV)
REAL(KIND=JPRB)    :: ZFPLS(KLON,0:KLEV),ZFPLC(KLON,0:KLEV),ZFPL(KLON,0:KLEV)
REAL(KIND=JPRB)    :: ZCSGC(KLON,KLEV)
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
REAL(KIND=JPRB) :: ZDCAPE(KLON) ! Descending CAPE for gusts.

! ACRANEB/ACRANEB2 local variables
! --------------------------------
REAL(KIND=JPRB) :: ZLAMB           ! proportion of Lambertian scattering
REAL(KIND=JPRB) :: ZALBDIR  (KLON) ! direct (parallel) surface albedo
REAL(KIND=JPRB) :: ZCLCT_RAD(KLON) ! total cloud cover for radiation
REAL(KIND=JPRB) :: ZDECRD   (KLON) ! decorrelation depth for cloud overlaps
REAL(KIND=JPRB) :: ZDECRD_MF(KLON) ! decorrelation depth for cloud overlaps
                                   ! in microphysics

REAL(KIND=JPRB) :: ZQG(KLON,KLEV)
REAL(KIND=JPRB) :: ZDQG

! IFS deep convection scheme local variables
! --------------------------------
INTEGER(KIND=JPIM) :: ITOPC(KLON),IBASC(KLON),ITYPE(KLON),ISPPN2D
INTEGER(KIND=JPIM) :: ICBOT(KLON),ICTOP(KLON),IBOTSC(KLON)
LOGICAL :: LLDSLPHY,LLPTQ,LLLAND(KLON),LLCUM(KLON),LLSC(KLON),LLSHCV(KLON),LLLINOX(KLON)
REAL(KIND=JPRB) :: ZLIGH_CTG(KLON),ZCTOPH(KLON),ZPRECMX(KLON),ZICE(KLON),ZCDEPTH(KLON),ZWMFU(KLON)
REAL(KIND=JPRB) :: ZVERVEL(KLON,KLEV), ZGEOM1(KLON,KLEV), ZGEOMH(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZTENQ(KLON,KLEV),ZTENU(KLON,KLEV),ZTENV(KLON,KLEV)
REAL(KIND=JPRB) :: ZLCRIT_AER(KLON,KLEV),ZSNDE(KLON,KLEV,2)
REAL(KIND=JPRB) :: ZCUCONVCA(KLON),ZGAW(KLON),ZVDIFTS,ZDXTDK(KLON)
REAL(KIND=JPRB) :: ZLU(KLON,KLEV),ZLUDE(KLON,KLEV),ZMFU(KLON,KLEV)
REAL(KIND=JPRB) :: ZLISUM(KLON,KLEV),ZMFD(KLON,KLEV)
REAL(KIND=JPRB) :: ZWMEAN(KLON),ZACPR(KLON),ZDIFF(KLON,KLEV)
REAL(KIND=JPRB) :: ZMFUDE_RATE(KLON,KLEV),ZMFDDE_RATE(KLON,KLEV)
REAL(KIND=JPRB) :: ZFHPCL(KLON,0:KLEV),ZFHPCN(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZFCQLF(KLON,0:KLEV),ZFCQLI(KLON,0:KLEV),ZFRSO(KLON,0:KLEV),ZFRTH(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZGP2DSPP(KLON,1),ZLUDELI(KLON,KLEV,4),ZLRAIN(KLON,KLEV),ZRSUD(KLON,KLEV,2)
REAL(KIND=JPRB), ALLOCATABLE :: ZCEN(:,:,:),ZSCAV(:)

! Precipitation type diagnostics
!--------------------------------
REAL(KIND=JPRB)   :: ZFPLSG(KLON,0:KLEV)

REAL(KIND=JPRB) :: ZENTCH(KLON,KLEV)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "acclph.intfb.h"
#include "accoefk.intfb.h"
#include "accvimp.intfb.h"
#include "accvimp_v3.intfb.h"
#include "accvimpgy.intfb.h"
#include "cucalln_mf.intfb.h"
#include "acdifoz.intfb.h"
#include "acdifus.intfb.h"
#include "acdnshf.intfb.h"
#include "acdrac.intfb.h"
#include "acnorgwd.intfb.h"
#include "acdrag.intfb.h"
#include "acdrme.intfb.h"
#include "acdrov.intfb.h"
#include "acfluso.intfb.h"
#include "achmt.intfb.h"
#include "acmixlenz.intfb.h"
#include "acmixlentm.intfb.h"
#include "acmixelen.intfb.h"
#include "acnebc.intfb.h"
#include "acnebn.intfb.h"
#include "acnebr.intfb.h"
#include "acozone.intfb.h"
#include "acpblh.intfb.h"
#include "acpblhtm.intfb.h"
#include "acpluie.intfb.h"
#include "acpluis.intfb.h"
#include "acpluiz.intfb.h"
#include "acradcoef.intfb.h"
#include "recmwf.intfb.h"
#include "acdayd.intfb.h"
#include "acrso.intfb.h"
#include "acraneb.intfb.h"
#include "acraneb2.intfb.h"
#include "acsol.intfb.h"
#include "actqsat.intfb.h"
#include "actqsats.intfb.h"
#include "acupu.intfb.h"
#include "acupm.intfb.h"
#include "acupd.intfb.h"
#include "acuptq.intfb.h"
#include "acveg.intfb.h"
#include "aplmphys.intfb.h"
#include "qngcor.intfb.h"
#include "radaer.intfb.h"
#include "radheat.intfb.h"
#include "radozcmf.intfb.h"
#include "suozon.intfb.h"
#include "actke.intfb.h"
#include "acdifv1.intfb.h"
#include "acdifv2.intfb.h"
#include "acdifv3.intfb.h"
#include "achmtls.intfb.h"
#include "aro_ground_param.h"
#include "aro_ground_diag.h"
#include "aro_ground_diag_2isba.h"
#include "aro_ground_diag_z0.h"
#include "accvud.intfb.h"
#include "accsu.intfb.h"
#include "acpcmt.intfb.h"
#include "acmodo.intfb.h"
#include "acnsdo.intfb.h"
#include "acnebcond.intfb.h"
#include "accdev.intfb.h"
#include "acnpart.intfb.h"
#include "acvppkf.intfb.h"
#include "acptke.intfb.h"
#include "arp_ground_param.intfb.h"
#include "wrphtrajtm_nl.intfb.h"
#include "accldia.intfb.h"
#include "surf_ideal_flux.intfb.h"
#include "actkecoefk.intfb.h"
#include "actkehmt.intfb.h"
#include "acmrip.intfb.h"
#include "acmris.intfb.h"
#include "acmriss.intfb.h"
#include "actkecoefkh.intfb.h"
#include "actkezotls.intfb.h"
#include "aro_mnhdust.h"
#include "wrscmr.intfb.h"
#include "checkmv.intfb.h"
#include "acaa1.intfb.h"
!include "chem_main.intfb.h"
#include "mean_rad_temp.intfb.h"
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
IF(LGCHECKMV) CALL CHECKMV(YDRIP,YDPHY0,YDPHY2,KIDIA,KFDIA,KLON,KLEV,KSTEP,YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF &
& ,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,PGELAM,PGEMU,YDMF_PHYS%TMP%RDG%MU0,YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%Q,YDMF_PHYS_STATE%YGSP_RR%T)
!     ------------------------------------------------------------------

LLREDPR=LCVCSD
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

ALLOCATE(ZSVM   (KLON,KLEV,NGFL_EXT))
ALLOCATE(ZSFSV  (KLON,NGFL_EXT))
ALLOCATE(ZPSV   (KLON,KLEV,NGFL_EXT))

!     ------------------------------------------------------------------
!     1.- INITIALISATIONS COMPLEMENTAIRES
!     -----------------------------------
IJN = KLON

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
  DO JLON=KIDIA,KFDIA
    ZDECRD(JLON)=RDECRD1+RDECRD2* &
     & EXP(-((ASIN(PGEMU(JLON))-RDECRD3*RDECLI)/RDECRD4)**2)
  ENDDO
ENDIF
IF ( RDECRD <= 0._JPRB ) THEN 
  DO JLON=KIDIA,KFDIA
    ZDECRD_MF(JLON)=ZDECRD(JLON)
  ENDDO
ELSE
  DO JLON=KIDIA,KFDIA
    ZDECRD_MF(JLON)=RDECRD
  ENDDO
ENDIF

!*        1.1 INITIALISATION DE L'OZONE

IF (LMPHYS) THEN
  IF (LOZONE) THEN

    ! L'ozone est initialise par PO3 (ozone GFL).
    ! Ozone is computed from PO3 (GFL ozone).

    ZEPSO3=1.E-11_JPRB
    DO JLON=KIDIA,KFDIA
      ZQO3(JLON,0)=1.E-9_JPRB
    ENDDO
    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA
        ZQO3(JLON,JLEV) = MAX(ZEPSO3,YDMF_PHYS_STATE%O3(JLON,JLEV))
      ENDDO
    ENDDO

  ELSEIF(KVCLIS == 1) THEN
    ZEPSO3=1.E-11_JPRB
    DO JLON=KIDIA,KFDIA
      ZQO3(JLON,0)=1.E-9_JPRB
    ENDDO
    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA
        ZQO3(JLON,JLEV) = PKOZO(JLON,JLEV,1)
      ENDDO
    ENDDO

  ELSEIF ((LO3FL).AND.(NOZOCL == 1).AND.(LRAYFM)) THEN
    IF (MOD(KSTEP,NRADFR) == 0) THEN
      CALL RADOZCMF(KIDIA,KFDIA,KLON,KLEV,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,PGEMU,ZROZ)
      DO JLON=KIDIA,KFDIA
        ZQO3(JLON,0)=1.E-9_JPRB
      ENDDO
      DO JLEV=KTDIA,KLEV
        DO JLON=KIDIA,KFDIA
          ZQO3(JLON,JLEV)=ZROZ(JLON,JLEV)*YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

  ELSE
    CALL SUOZON(KIDIA,KFDIA,KLON,KLEV,ZQO3,.FALSE.,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,LO3ABC,YDMF_PHYS_SURF%GSD_VC%PGROUP)
  ENDIF

!     GAZ CARBONIQUE.

  DO JLEV=KTDIA,KLEV
    DO JLON=KIDIA,KFDIA
      ZQCO2(JLON,JLEV)=QCO2
    ENDDO
  ENDDO


  ! INITIALISATION DE LA COORDONNEE ETA.
  ! INITIALISATION DE LA COORDONNEE ETA.

  ZVETAH(KTDIA-1)=STPREH(KTDIA-1)/VP00
  DO JLEV=KTDIA,KLEV
    ZVETAH(JLEV)=STPREH(JLEV)/VP00
    ZVETAF(JLEV)=STPRE (JLEV)/VP00
  ENDDO

!     EPAISSEUR STD AEROSOLS
  ZAEO = AERCS1*ZVETAH(KTDIA-1) + AERCS3*ZVETAH(KTDIA-1)**3&
   & + AERCS5*ZVETAH(KTDIA-1)**5
  DO JLEV = KTDIA, KLEV
    ZAEN = AERCS1*ZVETAH(JLEV) + AERCS3*ZVETAH(JLEV)**3&
     & + AERCS5*ZVETAH(JLEV)**5
    ZDAER(JLEV) = ZAEN - ZAEO
    ZAEO = ZAEN
  ENDDO

!     HUMIDITE CRITIQUE
  ZEPS=1.E-12_JPRB
  IF(LHUCN) THEN
    DO JLEV = KTDIA, KLEV
      ZHUC(JLEV)=1.0_JPRB-MAX( HUCOE*ZVETAF(JLEV)*(1.0_JPRB-ZVETAF(JLEV))/((&
      & 1.0_JPRB+HUTIL1*(ZVETAF(JLEV)-0.5_JPRB))*(&
      & 1.0_JPRB+HUTIL2*(ZVETAF(JLEV)-0.5_JPRB))),ZEPS)
    ENDDO
  ELSE
    DO JLEV = KTDIA, KLEV
      ZHUC(JLEV)=1.0_JPRB-MAX( HUCOE*ZVETAF(JLEV)**NPCLO1*(&
       & 1.0_JPRB-ZVETAF(JLEV))**NPCLO2*(&
       & 1.0_JPRB+SQRT(HUTIL)*(ZVETAF(JLEV)-0.5_JPRB)),ZEPS)
    ENDDO
  ENDIF

  IF ( LNEWSTAT ) THEN

    DO JLEV = KTDIA, KLEV-1
      ZMAN(JLEV) = FSM_CC * TANH(FSM_DD*ZVETAH(JLEV))
      ZMAK(JLEV) = FSM_EE * ZVETAH(JLEV)**FSM_FF +&
       & FSM_GG * (1-ZVETAH(JLEV))**FSM_HH + FSM_II
    ENDDO

!     MATHEMATICAL FILTER FOR THE EDGES

    IF ( LRPROX ) THEN
      DO JLEV = KTDIA, KTDIA+3
        ZMAK(JLEV) = ZMAK(JLEV) / 2**(KTDIA+4-JLEV)
      ENDDO
      DO JLEV = KLEV-4, KLEV-1
        ZMAK(JLEV) = ZMAK(JLEV) / 2**(5-KLEV+JLEV)
      ENDDO
    ENDIF

  ELSE
    DO JLEV = KTDIA, KLEV-1
      ZMAN(JLEV) = 0.3_JPRB*ZVETAH(JLEV)
      ZMAK(JLEV) = 0.1_JPRB
    ENDDO
  ENDIF   !LNEWSTAT

!3MT
!     INCREMENTAL CORRECTION FLUX FOR NEGAVTIVE HUMIDITY VALUES

  ZFCQVNG(:,:)=0.0_JPRB
  ZFCQING(:,:)=0.0_JPRB
  ZFCQLNG(:,:)=0.0_JPRB

  ZPRODTH_CVPP(:,:)=0.0_JPRB

  ZGDT=RG*TSPHY
  ZGDTI=1.0_JPRB/ZGDT

  DO JLEV=KTDIA,KLEV
    DO JLON=KIDIA,KFDIA

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
  IF (LMDUST.AND.(NGFL_EXT/=0)) INBTRA=NGFL_EXT
  IF (NAERO>0)                  INBTRA=NAERO           ! the two cases exclude each other
  ALLOCATE(ZSTRCTRA(KLON,0:KLEV,INBTRA)) ! to cover both prog aero and extra gfl cases
  ALLOCATE(ZTRA    (KLON,  KLEV,INBTRA))
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
  IF (LMDUST.AND.(NGFL_EXT/=0)) THEN
!  SIZE OF ARRAY FOR MSE
  ILONMNH=KFDIA-KIDIA+1
  ZINVDT=1/PDT
  ZINVG=1/RG
  ZZI_APHI=YDMF_PHYS_STATE%YCPG_DYN%PHI
  ZTM=YDMF_PHYS_STATE%T
!SETUP
  ZZI_SVM=0.0_JPRB
  ZZI_PEZDIAG=0.0_JPRB
  ZZI_PABSM=101325._JPRB
  ZZI_RHODREFM=1.0_JPRB
  ZZI_RHO=1.0_JPRB
  ZZZ=0.0_JPRB
  ZAERD=0.0_JPRB
  PEZDIAG=0.0_JPRB

!Initialisation des scalaires passifs pour aro_ground_param
  DO JGFL=1,NGFL_EXT
    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA
       ZSVM(JLON,JLEV,JGFL)=MAX(YDMF_PHYS_STATE%P1EXT(JLON,JLEV,JGFL),0.0_JPRB)
      ENDDO
    ENDDO
  ENDDO

  !initialisation de ZZZ
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      ZZZ(JLON,1,JLEV)=ZZI_APHI(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO

  !initialisation de ZDZZ
  DO JLEV = 2, KLEV
    DO JLON = KIDIA,KFDIA
      ZDZZ(JLON,1,JLEV)=ZZZ(JLON,1,JLEV-1)-ZZZ(JLON,1,JLEV)
    ENDDO
  ENDDO
  DO JLON = KIDIA,KFDIA
      ZDZZ(JLON,1,1)=ZZI_APHI(JLON,0)*ZINVG-ZZZ(JLON,1,1)
  ENDDO


!Initialisation de ZZI_RHODREFM
  DO JLEV = 1 , KLEV
    DO JLON=KIDIA,KFDIA
      ZZI_RHODREFM(JLON,1,JLEV)=YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(JLON,JLEV)/&
       & (YDMF_PHYS_STATE%T(JLON,JLEV)*YDMF_PHYS_STATE%YCPG_DYN%RCP%R(JLON,JLEV))
     ENDDO
  ENDDO

!Initialisation de ZZI_PABSM
  DO JLON = KIDIA,KFDIA
    DO JLEV = 1 , KLEV
      ZZI_PABSM(JLON,1,JLEV)=YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(JLON,JLEV)
    ENDDO
  ENDDO

!Initialisation de ZZI_EXNREFM
  ZRSCP=RD/RCPD
  ZINVATM=1/RATM
  DO JLEV = 1 , KLEV
    DO JLON = KIDIA,KFDIA
      ZZI_EXNREFM(JLON,1,JLEV)=(YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(JLON,JLEV)*ZINVATM)**(ZRSCP)
    ENDDO
  ENDDO

!Initialisation de ZZI_THM
  DO JLEV = 1 , KLEV
    DO JLON = KIDIA,KFDIA
      ZZI_THM(JLON,1,JLEV)=ZTM(JLON,JLEV)/ZZI_EXNREFM(JLON,1,JLEV)
    ENDDO
  ENDDO

!Initialisation des scalaires passifs pour aro_mnhdust (inversion des niveaux)
  DO JGFL=1,NGFL_EXT
    DO JLON=1,KLON
      DO JLEV=1,KLEV
        ZZI_SVM(JLON,1,JLEV,JGFL)=ZSVM(JLON,JLEV,JGFL)
      ENDDO
    ENDDO
  ENDDO

  ENDIF  ! ENDIF  (LMDUST & NGFL_EXT)
ENDIF  !LMPHYS
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

IF (LCONDWT) THEN
  IF (L3MT .OR. LSTRAPRO .OR. LPROCLD) THEN
    IF (LGRAPRO) THEN
!cdir unroll=8
      DO JLEV = 1, KLEV
        DO JLON = KIDIA,KFDIA
          ZQI(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_STATE%I(JLON,JLEV))
          ZQL(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_STATE%L(JLON,JLEV))
          ZQR(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_STATE%R(JLON,JLEV))
          ZQS(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_STATE%S(JLON,JLEV))
          ZQG(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_STATE%G(JLON,JLEV))
  ! CORRECTION OF NEGATIVE ADVECTED VALUES:
  !                         VAPOUR PUT IN PFCQVNG
  !                         LIQUID,ICE PUT IN PFCQL/ING
  !                         FOR RAIN/SNOW/GRAUPEL PUT IN PFCQR/S/GNG

          ZDQI=ZQI(JLON,JLEV)-YDMF_PHYS_STATE%I(JLON,JLEV)
          ZDQL=ZQL(JLON,JLEV)-YDMF_PHYS_STATE%L(JLON,JLEV)
          ZDQR=ZQR(JLON,JLEV)-YDMF_PHYS_STATE%R(JLON,JLEV)
          ZDQS=ZQS(JLON,JLEV)-YDMF_PHYS_STATE%S(JLON,JLEV)
          ZDQG=ZQG(JLON,JLEV)-YDMF_PHYS_STATE%G(JLON,JLEV)
          ZDQC=ZDQI+ZDQL+ZDQR+ZDQS+ZDQG

          ZQV0=YDMF_PHYS_STATE%Q(JLON,JLEV)-ZIPOI(JLON,JLEV)*(0.0_JPRB- YDMF_PHYS%OUT%FCQNG(JLON,JLEV-1) &
          & -YDMF_PHYS%OUT%FCQNNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQLNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQRNG(JLON,JLEV-1) &
          & -YDMF_PHYS%OUT%FCQSNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQGNG(JLON,JLEV-1))
          ZQV(JLON,JLEV)=MAX(0.0_JPRB,ZQV0-ZDQC)
          ZDQV=MAX(0.0_JPRB,ZQV0-ZDQC)-YDMF_PHYS_STATE%Q(JLON,JLEV)

          YDMF_PHYS%OUT%FCQNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNG(JLON,JLEV-1)-ZDQV*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQNNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNNG(JLON,JLEV-1)-ZDQI*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQLNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQLNG(JLON,JLEV-1)-ZDQL*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQRNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQRNG(JLON,JLEV-1)-ZDQR*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQSNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQSNG(JLON,JLEV-1)-ZDQS*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQGNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQGNG(JLON,JLEV-1)-ZDQG*ZPOID(JLON,JLEV)
        ENDDO
      ENDDO
    ELSE
      DO JLEV = 1, KLEV
        DO JLON = KIDIA,KFDIA
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
    ENDIF

    IF (LGPCMT) THEN
      DO JLEV = 1, KLEV
        DO JLON = KIDIA,KFDIA
          ZTZER=MAX(0.0_JPRB,-SIGN(1.0_JPRB,-(YDVARS%ICONV%T0(JLON,JLEV)-ZEPS0)+0.0_JPRB))
          ZQIC(JLON,JLEV)=ZTZER*YDVARS%ICONV%T0(JLON,JLEV)
          ZDQI=ZQIC(JLON,JLEV)-YDVARS%ICONV%T0(JLON,JLEV)

          ZTZER=MAX(0.0_JPRB,-SIGN(1.0_JPRB,-(YDVARS%LCONV%T0(JLON,JLEV)-ZEPS0)+0.0_JPRB))
          ZQLC(JLON,JLEV)=ZTZER*YDVARS%LCONV%T0(JLON,JLEV)
          ZDQL=ZQLC(JLON,JLEV)-YDVARS%LCONV%T0(JLON,JLEV)

          ZTZER=MAX(0.0_JPRB,-SIGN(1.0_JPRB,-(YDVARS%RCONV%T0(JLON,JLEV)-ZEPS0)+0.0_JPRB))
          ZQRC(JLON,JLEV)=ZTZER*YDVARS%RCONV%T0(JLON,JLEV)
          ZDQR=ZQRC(JLON,JLEV)-YDVARS%RCONV%T0(JLON,JLEV)

          ZTZER=MAX(0.0_JPRB,-SIGN(1.0_JPRB,-(YDVARS%SCONV%T0(JLON,JLEV)-ZEPS0)+0.0_JPRB))
          ZQSC(JLON,JLEV)=ZTZER*YDVARS%SCONV%T0(JLON,JLEV)
          ZDQS=ZQSC(JLON,JLEV)-YDVARS%SCONV%T0(JLON,JLEV)

          ZDQC=ZDQI+ZDQL+ZDQR+ZDQS

          ZQV0=ZQV(JLON,JLEV)-ZIPOI(JLON,JLEV)*(0.0_JPRB-ZFCQVNG(JLON,JLEV-1)&
          & -YDMF_PHYS%OUT%FCNEGQIC(JLON,JLEV-1)-YDMF_PHYS%OUT%FCNEGQLC(JLON,JLEV-1)&
          & -YDMF_PHYS%OUT%FCNEGQRC(JLON,JLEV-1)-YDMF_PHYS%OUT%FCNEGQSC(JLON,JLEV-1))
          ZQVI=MAX(0.0_JPRB,ZQV0-ZDQC)
          ZDQV=ZQVI-ZQV(JLON,JLEV)
          ZQV(JLON,JLEV)=ZQVI

          ZFCQVNG(JLON,JLEV)=ZFCQVNG(JLON,JLEV-1)-ZDQV*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNG(JLON,JLEV)+ZFCQVNG(JLON,JLEV)
          YDMF_PHYS%OUT%FCNEGQIC(JLON,JLEV)=YDMF_PHYS%OUT%FCNEGQIC(JLON,JLEV-1)-ZDQI*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCNEGQLC(JLON,JLEV)=YDMF_PHYS%OUT%FCNEGQLC(JLON,JLEV-1)-ZDQL*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCNEGQRC(JLON,JLEV)=YDMF_PHYS%OUT%FCNEGQRC(JLON,JLEV-1)-ZDQR*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCNEGQSC(JLON,JLEV)=YDMF_PHYS%OUT%FCNEGQSC(JLON,JLEV-1)-ZDQS*ZPOID(JLON,JLEV)
        ENDDO
      ENDDO
      ZFCQVNG(:,:)=0.0_JPRB
    ENDIF
  ELSE
    DO JLEV = 1, KLEV
      DO JLON = KIDIA,KFDIA
        ZQI(JLON,JLEV)=YDMF_PHYS_STATE%I(JLON,JLEV)
        ZQL(JLON,JLEV)=YDMF_PHYS_STATE%L(JLON,JLEV)
        ZQR(JLON,JLEV)=YDMF_PHYS_STATE%R(JLON,JLEV)
        ZQS(JLON,JLEV)=YDMF_PHYS_STATE%S(JLON,JLEV)
        IF (LGRAPRO) THEN
          ZQG(JLON,JLEV)=YDMF_PHYS_STATE%G(JLON,JLEV)
        ENDIF
        ZQV(JLON,JLEV)=YDMF_PHYS_STATE%Q(JLON,JLEV)
      ENDDO
    ENDDO
  ENDIF
ELSE
  DO JLEV = 1, KLEV
    DO JLON = KIDIA,KFDIA
      ZQI(JLON,JLEV)=0.0_JPRB
      ZQL(JLON,JLEV)=0.0_JPRB
      ZQV(JLON,JLEV)=YDMF_PHYS_STATE%Q(JLON,JLEV)
    ENDDO
  ENDDO
ENDIF ! LCONDWT

DO JCHA = 1, 6
  DO JLEV = KTDIA, KLEV
    DO JLON = 1, KLON
      ZAER(JLON,JLEV,JCHA)=0.0_JPRB
    ENDDO
  ENDDO
ENDDO

DO JLEV = KTDIA, KLEV
  DO JLON = 1, KLON
    ZAERINDS(JLON,JLEV)=0.0_JPRB
  ENDDO
ENDDO

DO JLON = KIDIA,KFDIA
  ZCEMTR  (JLON,0) = 0.0_JPRB
  ZCEMTR  (JLON,1) = 0.0_JPRB
  ZCTRSO  (JLON,0) = 0.0_JPRB
  ZCTRSO  (JLON,1) = 0.0_JPRB
  ZTRSOD  (JLON)   = 0.0_JPRB
  ZSUDU   (JLON)   = 0.0_JPRB
  ZXDROV  (JLON)   = 1.0_JPRB
  ZXHROV  (JLON)   = 1.0_JPRB
  ZTAUX   (JLON)   = ZEPS0
  KCLPH   (JLON)   = KLEV
ENDDO
DO JSG = 1, NSW
  DO JLON = KIDIA,KFDIA
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

IF((LNEBR.OR.(TRIM(CGMIXLEN) == 'TM')&
       & .OR.(TRIM(CGMIXLEN) == 'TMC')).AND.KSTEP == 0) THEN
!DEC$ IVDEP
  DO JLON=KIDIA,KFDIA
    YDMF_PHYS_SURF%GSD_VH%PPBLH(JLON)=(XMINLM+XMAXLM)*0.5_JPRB
  ENDDO
ENDIF
IF((LNEBCO.OR.LGWDC).AND.KSTEP == 0) THEN
  DO JLON=KIDIA,KFDIA
    YDMF_PHYS_SURF%GSD_VH%PTCCH(JLON)=0.0_JPRB
    YDMF_PHYS_SURF%GSD_VH%PSCCH(JLON)=0.0_JPRB
    YDMF_PHYS_SURF%GSD_VH%PBCCH(JLON)=0.0_JPRB
  ENDDO
ENDIF

IF((LNEBN.OR.LNEBR.OR.LRRGUST).AND.KSTEP == 0) THEN
  DO JLEV=KTDIA-1,KLEV
    DO JLON=KIDIA,KFDIA
      YDMF_PHYS%TMP%APLPAR%PFL%FPLCH(JLON,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF
IF(LRRGUST.AND.KSTEP == 0) THEN
  DO JLEV=KTDIA-1,KLEV
    DO JLON=KIDIA,KFDIA
      YDMF_PHYS%TMP%APLPAR%PFL%FPLSH(JLON,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF
IF(LCVPGY.AND.KSTEP == 0) THEN
  DO JLEV=KTDIA,KLEV
    DO JLON=KIDIA,KFDIA
      YDMF_PHYS_STATE%CVV(JLON,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF
IF(LPHSPSH.AND.KSTEP == 0) THEN
  YDMF_PHYS_SURF%GSD_VH%PSPSH(KIDIA:KFDIA)=0._JPRB
ENDIF
IF(LCOEFKTKE.AND.KSTEP == 0) THEN
  YDMF_PHYS_SURF%GSD_VH%PQSH(KIDIA:KFDIA)=YDMF_PHYS_STATE%Q(KIDIA:KFDIA,KLEV)
ENDIF
IF(LCVCSD.AND.LUDEVOL.AND.KSTEP==0) THEN
  YDMF_PHYS_SURF%GSD_VK%PUDGRO(KIDIA:KFDIA)=0._JPRB
ENDIF
IF(LRKCDEV.AND.KSTEP == 0) THEN
  DO JLEV=KTDIA,KLEV
     DO JLON=KIDIA,KFDIA
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
IF (L3MT) THEN
  IF (LCVPRO) THEN
    YDVARS%UNEBH%T0(KIDIA:KFDIA,:)=MAX(0._JPRB,YDVARS%UNEBH%T0(KIDIA:KFDIA,:))
    ZUNEBH(:,:)=MIN(1.0_JPRB-ZEPS0,YDVARS%UNEBH%T0(:,:)+MAX(0._JPRB,YDVARS%UAL%T0(:,:)))
  ELSE
    ZUNEBH(:,:)=YDVARS%UNEBH%T0(:,:)
  ENDIF
ENDIF

!     ------------------------------------------------------------------
!  The LMPHYS and LEPHYS keys should in the following be
! dispached to individual parametrizations.
IF(LMPHYS) THEN
!*
!     ------------------------------------------------------------------
!     4.- CALCULS THERMODYNAMIQUES
!     ----------------------------
  IF ( LTHERMO ) THEN
    CALL ACTQSAT ( YDPHY,KIDIA,KFDIA,KLON,NTQSAT,KLEV,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, ZQV, YDMF_PHYS_STATE%T,&
     & ZGEOSLC, YDMF_PHYS%TMP%APLPAR%MSC%LH, YDMF_PHYS%TMP%APLPAR%MSC%LSCPE, YDMF_PHYS%TMP%APLPAR%FLU%QSAT, YDMF_PHYS%TMP%APLPAR%MSC%QW, YDCPG_MISC%RH, YDMF_PHYS%TMP%APLPAR%MSC%TW)
  ENDIF

!*
!     ------------------------------------------------------------------
!     4.BIS. COEFFICIENTS THERMO-HYDRIQUES DU SOL
!     -------------------------------------------

  IF (LSFORCS) THEN  ! Surface forcing for 1D model MUSC
     DO JLON=KIDIA,KFDIA
        ZTSN(JLON)=YDMF_PHYS_STATE%YGSP_RR%T(JLON)
        ZTN(JLON) =YDMF_PHYS_STATE%T(JLON,KLEV)
     ENDDO
     DO JLON=KIDIA,KFDIA
        ZRHODREFM(JLON) = YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(JLON,KLEV)/(YDMF_PHYS_STATE%T(JLON,KLEV)*YDMF_PHYS_STATE%YCPG_DYN%RCP%R(JLON,KLEV))
        ZTHETAS(JLON)   = ZTSN(JLON)*(RATM/YDMF_PHYS_STATE%YCPG_PHY%PREHYD(JLON,KLEV))**RKAPPA
     ENDDO
     LLAROME = .FALSE.
     CALL SURF_IDEAL_FLUX(YDRIP,YDPHY0,YDPHYDS,LLAROME,KIDIA,KFDIA,KLON,YDMF_PHYS_STATE%YCPG_DYN%PHIF(:,KLEV),ZRHODREFM,YDMF_PHYS_SURF%GSD_SFO%PGROUP,&
        & ZTN,ZTSN,YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS_STATE%Q(:,KLEV),YDMF_PHYS_STATE%U(:,KLEV),YDMF_PHYS_STATE%V(:,KLEV),ZTHETAS,YDMF_PHYS%OUT%FCS(KIDIA:KFDIA,1),ZFEV,ZFMDU,ZFMDV)
     DO JLON=KIDIA,KFDIA
        YDMF_PHYS_STATE%YGSP_RR%T(JLON)=ZTSN(JLON)
     ENDDO
     DO JLON=KIDIA,KFDIA
        ZDELTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-YDMF_PHYS_STATE%YGSP_RR%T(JLON)))
        ! To be equivalent to surfex forcing
        YDMF_PHYS%OUT%FEVL(JLON,1)=ZFEV(JLON)*(1.0_JPRB-ZDELTA)
        YDMF_PHYS%OUT%FEVN(JLON,1)=ZFEV(JLON)*ZDELTA
        YDMF_PHYS%OUT%FCLL(JLON,1)=YDMF_PHYS%OUT%FEVL(JLON,1)*YDMF_PHYS%TMP%APLPAR%MSC%LH(JLON,KLEV)
        YDMF_PHYS%OUT%FCLN(JLON,1)=YDMF_PHYS%OUT%FEVN(JLON,1)*YDMF_PHYS%TMP%APLPAR%MSC%LH(JLON,KLEV)
        YDMF_PHYS%TMP%APLPAR%DSA%LHS(JLON)=FOLH(YDMF_PHYS_STATE%YGSP_RR%T(JLON),0.0_JPRB)
     ENDDO
  ENDIF  ! End of surface forcing for 1D model MUSC

  IF ( .NOT.LMSE ) THEN
    IF ( LSOLV ) THEN
      LLHMT=.FALSE.
      CALL ACSOL ( YDPHY,YDPHY1,KIDIA,KFDIA,KLON,&
       & YDMF_PHYS_SURF%GSD_VV%PARG,YDMF_PHYS_SURF%GSD_VV%PD2,YDMF_PHYS_SURF%GSD_VF%PZ0F,YDMF_PHYS_SURF%GSD_VV%PZ0H,YDMF_PHYS_SURF%GSD_VF%PZ0RLF,YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS_SURF%GSD_VV%PIVEG,YDMF_PHYS_SURF%GSD_VV%PLAI,YDMF_PHYS_STATE%YGSP_SG%A,&
       & YDMF_PHYS_STATE%YGSP_SG%R,YDMF_PHYS_SURF%GSD_VV%PSAB,YDMF_PHYS_STATE%YGSP_SG%F,YDMF_PHYS_STATE%YGSP_RR%T,YDMF_PHYS_SURF%GSD_VF%PVEG,YDMF_PHYS_STATE%YGSP_SB%Q,YDMF_PHYS_STATE%YGSP_SB%TL,YDMF_PHYS_STATE%YGSP_RR%W,YDMF_PHYS_STATE%YGSP_RR%IC,&
       & LLHMT,&
       & YDMF_PHYS%TMP%APLPAR%DSA%C1,YDMF_PHYS%TMP%APLPAR%DSA%C2,ZC3,ZCG,ZCN,YDMF_PHYS%OUT%CT,ZNEIJG,ZNEIJV,&
       & ZWFC,ZWPMX,ZWSEQ,ZWSMX,ZWWILT)
    ELSE

!            INITIALISATION DE L'INERTIE THERMIQUE DU SOL.

!DEC$ IVDEP
      DO JLON = KIDIA,KFDIA
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

    IF ( (NSWB_MNH /= NSW) .AND. LRAYFM ) CALL ABOR1('APLPAR: NSWB_MNH not = NSW')

    DO JLON=KIDIA,KFDIA
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

    IF ( LRAY ) THEN
      ! ACRANEB/ACRANEB2 radiation => one solar band
      DO JLON=KIDIA,KFDIA
        YDMF_PHYS%OUT%ALB   (JLON) = PGPAR(JLON,MALBSCA)
        ZALBDIR(JLON) = PGPAR(JLON,MALBDIR)
      ENDDO
    ELSE
      ! FMR  radiation => NSW solar bands
      DO JSG=1,NSW
        DO JLON=KIDIA,KFDIA
          ZALBP(JLON,JSG) = PGPAR(JLON,MALBDIR-1+JSG)
          ZALBD(JLON,JSG) = PGPAR(JLON,MALBSCA-1+JSG)
        ENDDO
      ENDDO
      DO JLON=KIDIA,KFDIA
        YDMF_PHYS%OUT%ALB(JLON)=0.0_JPRB
      ENDDO
      DO JSG=1,NSW
        DO JLON=KIDIA,KFDIA
          YDMF_PHYS%OUT%ALB(JLON)=YDMF_PHYS%OUT%ALB(JLON)+0.5*(ZALBP(JLON,JSG)+ZALBD(JLON,JSG))
        ENDDO
      ENDDO
      DO JLON=KIDIA,KFDIA
        YDMF_PHYS%OUT%ALB(JLON)=YDMF_PHYS%OUT%ALB(JLON)/REAL(NSW,JPRB)
      ENDDO
    ENDIF

    DO JLON=KIDIA,KFDIA
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
  IF (LMSE.AND.(.NOT.LCOEFKTKE).AND.(.NOT.LCOEFK_TOMS).AND.KSTEP == 0) THEN
  CALL ARO_GROUND_DIAG_Z0( KBL, KGPCOMP,&
          & KFDIA-KIDIA+1, KIDIA, KFDIA,&
          & NDGUNG, NDGUXG, NDLUNG, NDLUXG,&
          & PINDX(KIDIA:KFDIA), PINDY(KIDIA:KFDIA),&
          & LSURFEX_KFROM,&
          & YDMF_PHYS%OUT%GZ0(KIDIA:KFDIA), YDMF_PHYS%OUT%GZ0H(KIDIA:KFDIA))
  ENDIF
!*
!     ------------------------------------------------------------------
!     5.- STRUCTURE ET CHAMPS DANS LA COUCHE LIMITE DE SURFACE
!     ------------------------------------------------------------------

!        INITIALISATION DES HAUTEURS "METEO".

  IF ( LHMTO ) THEN
    DO JLON=KIDIA,KFDIA
      ZDPHIV(JLON)=RG*HVCLS
      ZDPHIT(JLON)=RG*HTCLS
    ENDDO
  ENDIF

  IF ( LVDIF.OR.LHMTO.OR.LGWD ) THEN
    LLCLS=LGWD.OR.LVDIF
    LLHMT=LHMTO
    IF (LMSE) THEN
      IF(LCOEFKTKE.AND.LCOEFKSURF) THEN

       IF (KSTEP == 0) THEN 
         CALL ARO_GROUND_DIAG_Z0( KBL, KGPCOMP, &
          & KFDIA-KIDIA+1, KIDIA, KFDIA, &
          & NDGUNG, NDGUXG, NDLUNG, NDLUXG, &
          & PINDX(KIDIA:KFDIA), PINDY(KIDIA:KFDIA), &
          & LSURFEX_KFROM, &
          & YDMF_PHYS%OUT%GZ0(KIDIA:KFDIA), YDMF_PHYS%OUT%GZ0H(KIDIA:KFDIA))
        ENDIF


        CALL ACTKEZOTLS ( YDPHY0,KIDIA,KFDIA,KLON,KLEV,&
         & YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_STATE%T,&
         & YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, YDMF_PHYS_STATE%YGSP_RR%T, YDCPG_MISC%QS,&
         & YDMF_PHYS%TMP%APLPAR%FLU%CDN, ZCDNH, YDMF_PHYS%TMP%APLPAR%DSA%CPS, ZRTI, YDMF_PHYS%TMP%APLPAR%DSA%RS)
         
      ELSE
        CALL ACHMTLS ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,KLEV,&
          & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,&
          & YDMF_PHYS_STATE%YGSP_RR%T,YDCPG_MISC%QS,ZDPHIT,YDMF_PHYS%OUT%GZ0,YDMF_PHYS%OUT%GZ0H,LLCLS,&
          & ZNBVNO,ZMRIPP,YDMF_PHYS%TMP%APLPAR%DSA%CPS,ZGWDCS,YDMF_PHYS%TMP%APLPAR%DSA%LHS,ZPCLS,YDMF_PHYS%TMP%APLPAR%FLU%CD,YDMF_PHYS%TMP%APLPAR%FLU%CDN)
!       Computation of ZRTI
        DO JLON=KIDIA,KFDIA
          ZDPHI(JLON)=YDMF_PHYS_STATE%YCPG_DYN%PHIF(JLON,KLEV)-YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,KLEV)
          ZPRS(JLON)=RD+ZRVMD*YDCPG_MISC%QS(JLON)
          ZRTI(JLON)=2.0_JPRB/(YDMF_PHYS_STATE%YCPG_DYN%RCP%R(JLON,KLEV)*YDMF_PHYS_STATE%T(JLON,KLEV)+RKAPPA*ZDPHI(JLON)&
           & +ZPRS(JLON)*YDMF_PHYS_STATE%YGSP_RR%T(JLON))
        ENDDO
      ENDIF
    ELSE
      IF (LCOEFKSURF) THEN
        CALL ACTKEHMT ( KIDIA,KFDIA,KLON,KLEV,KVCLIV>=8.AND.LSOLV,&
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
      ELSE
        CALL ACHMT ( KIDIA,KFDIA,KLON,KLEV,KVCLIV>=8.AND.LSOLV,&
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
    ENDIF

    IF (LPTKE) THEN
      YDMF_PHYS_STATE%TKE(KIDIA:KFDIA,1:KLEV) = MAX(YDMF_PHYS_STATE%TKE(KIDIA:KFDIA,1:KLEV),ETKE_MIN)
    ENDIF
    IF (LCOEFK_PTTE) THEN
      YDVARS%TTE%T0(KIDIA:KFDIA,1:KLEV) = MAX(YDVARS%TTE%T0(KIDIA:KFDIA,1:KLEV),ETKE_MIN)
    ENDIF

    IF(LCOEFKTKE) THEN
      ZCP(KIDIA:KFDIA,1:KLEV) = RCPD*(1.0_JPRB+(RCPV/RCPD-1.0_JPRB)*(&
        & ZQV(KIDIA:KFDIA,1:KLEV)+ZQI(KIDIA:KFDIA,1:KLEV)+&
        & ZQL(KIDIA:KFDIA,1:KLEV)))
    ELSE
      ZCP(KIDIA:KFDIA,1:KLEV) = YDMF_PHYS_STATE%YCPG_DYN%RCP%CP(KIDIA:KFDIA,1:KLEV)
    ENDIF

    IF(LCOEFK_RIS .AND. LCOEFKTKE) THEN
      !  computation of Ri*,Ri** for mixing lenth computation
      CALL ACMRISS ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTCOEF,KLEV,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, ZQV, ZQL, ZQI, YDMF_PHYS%TMP%APLPAR%FLU%QSAT,&
       & YDMF_PHYS_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_STATE%T, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, YDMF_PHYS%TMP%APLPAR%MSC%LSCPE,YDMF_PHYS%OUT%GZ0,&
       & ZMN2PP,ZMRIPP)
    ENDIF

    ! COMPUTATION OF mixing lengths from  Ri*,Ri** - FIRST GUES for moist AF

    !---------------------------------------------------
    ! COMPUTATION OF 'DRY' mixing lengths : lm_d lh_d
    ! COMPUTATION OF ZPBLH - PBL HEIGHT

    IF (CGMIXLEN == 'Z'  .OR. &
     &  CGMIXLEN == 'EL0'.OR. &
     &  CGMIXLEN == 'EL1'.OR. &
     &  CGMIXLEN == 'EL2'.OR. &
     &  CGMIXLEN == 'AY' .OR. &
     &  CGMIXLEN == 'AYC'.AND.(.NOT.LECT)) THEN
      DO JLEV=KTDIA,KLEV
        DO JLON=KIDIA,KFDIA
          ZTHETAV(JLON,JLEV)=YDMF_PHYS_STATE%T(JLON,JLEV)*(1.0_JPRB+RETV*ZQV(JLON,JLEV))&
           & *(RATM/YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(JLON,JLEV))**RKAPPA  
        ENDDO
      ENDDO
      DO JLON=KIDIA,KFDIA
        ZTHETAVS(JLON)=YDMF_PHYS_STATE%YGSP_RR%T(JLON)*(1.0_JPRB+RETV*YDCPG_MISC%QS(JLON))&
         & *(RATM/YDMF_PHYS_STATE%YCPG_PHY%PREHYD(JLON,KLEV))**RKAPPA  
      ENDDO
      CALL ACCLPH ( YDPHY0,YDPHY2,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
       & ZTHETAV,YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,ZTHETAVS,&
       & KCLPH,YDMF_PHYS%OUT%CLPH,YDMF_PHYS%OUT%VEIN,ZUGST,ZVGST)
      ZBLH(KIDIA:KFDIA)=YDMF_PHYS%OUT%CLPH(KIDIA:KFDIA)
      IF (.NOT.LRAFTUR) THEN
        YDMF_PHYS%OUT%UGST(KIDIA:KFDIA)=ZUGST(KIDIA:KFDIA)
        YDMF_PHYS%OUT%VGST(KIDIA:KFDIA)=ZVGST(KIDIA:KFDIA)
      ENDIF
    ELSE
      ZBLH(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VH%PPBLH(KIDIA:KFDIA)
    ENDIF

  ENDIF ! end of LVDIF or LHMTO or LGWD

  IF ( (LVDIF.OR.LGWD).AND.(.NOT.(LNEBR.OR.LECT)) ) THEN

    IF(TRIM(CGMIXLEN) == 'Z') THEN
      !-------------------------------------------------
      ! "z dependent" mixing length.
      !-------------------------------------------------
      ZBLH(KIDIA:KFDIA)=1.0_JPRB/UHDIFV
      CALL ACMIXLENZ ( YDPHY,YDPHY0,KIDIA,KFDIA,KLON,1,KLEV,.FALSE.,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,ZBLH,YDMF_PHYS%OUT%GZ0,YDMF_PHYS%OUT%GZ0H,ZLMU,ZLMT)

    ELSEIF((TRIM(CGMIXLEN) == 'TMC').OR.(TRIM(CGMIXLEN) == 'AYC')) THEN
      !     Cubique du climat
      ZBLH(KIDIA:KFDIA)=MIN(XMAXLM,MAX(XMINLM,ZBLH(KIDIA:KFDIA)))
      CALL ACMIXLENTM ( YDPHY0,KIDIA,KFDIA,KLON,KLEV, &
       & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS%OUT%GZ0,YDMF_PHYS%OUT%GZ0H,ZBLH,ZLMU,ZLMT)

    ELSEIF(TRIM(CGMIXLEN) == 'TM') THEN
      !     Ancienne formulation pour Lm
      ZBLH(KIDIA:KFDIA)=MIN(XMAXLM,MAX(XMINLM,ZBLH(KIDIA:KFDIA)))*XKLM
      CALL ACMIXLENZ ( YDPHY,YDPHY0,KIDIA,KFDIA,KLON,1,KLEV,.FALSE.,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,ZBLH,YDMF_PHYS%OUT%GZ0,YDMF_PHYS%OUT%GZ0H,ZLMU,ZLMT)

      !     Cubique du climat pour Lh
      ZBLH(KIDIA:KFDIA)=ZBLH(KIDIA:KFDIA)/XKLM
      CALL ACMIXLENTM ( YDPHY0,KIDIA,KFDIA,KLON,KLEV,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS%OUT%GZ0,YDMF_PHYS%OUT%GZ0H,ZBLH,ZZLMT,ZLMT)

    ELSEIF(TRIM(CGMIXLEN) == 'AY') THEN
      !     new Ayotte-Tudor ZBLH & mixing length
      CALL ACMIXLENZ ( YDPHY,YDPHY0,KIDIA,KFDIA,KLON,1,KLEV,.FALSE.,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,ZBLH,YDMF_PHYS%OUT%GZ0,YDMF_PHYS%OUT%GZ0H,ZLMU,ZLMT)

    ELSEIF((CGMIXLEN(1:2) == 'EL').AND.LPTKE) THEN
      !     e-type mixing length converted to Prandtl type
      CALL ACMIXLENZ(YDPHY,YDPHY0,KIDIA,KFDIA,KLON,1,KLEV,.TRUE.,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,ZBLH,YDMF_PHYS%OUT%GZ0,YDMF_PHYS%OUT%GZ0H,ZLMU,ZLMT)
      ZLMU2(:,:)=ZLMU(:,:)
      ZLMT2(:,:)=ZLMT(:,:)
      IF     (CGMIXLEN == 'EL0') THEN
        IMLTYPE=0
        ! to have identical mixing length like in pTKE
        ! e-type mixing length converted to Prandtl type
        CALL ACMIXLENZ(YDPHY,YDPHY0,KIDIA,KFDIA,KLON,1,KLEV,.FALSE.,&
         & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,ZBLH,YDMF_PHYS%OUT%GZ0,YDMF_PHYS%OUT%GZ0H,ZLMU,ZLMT)
        ZLMU2(:,:)=ZLMU(:,:)
        ZLMT2(:,:)=ZLMT(:,:)
      ELSEIF (CGMIXLEN == 'EL1') THEN
        IMLTYPE=1
      ELSEIF (CGMIXLEN == 'EL2') THEN
        IMLTYPE=2
      ELSE
        CLERR='APLPAR: UNEXPECTED VALUE FOR CGMIXLEN: '//TRIM(CGMIXLEN)
        CALL ABOR1(CLERR)
      ENDIF

      IF( LCOEFK_RIS) THEN
        LLMAF=.TRUE.
        CALL ACMIXELEN(YGFL,YDPHY,YDPHY0,&
         & KIDIA,KFDIA,KLON,KTDIA,KLEV,KSTEP,IMLTYPE,&
         & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%T,ZQV,ZQL,ZQI,YDMF_PHYS_STATE%R,YDMF_PHYS_STATE%S,YDMF_PHYS_STATE%TKE,&
         & YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH,YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR,ZMN2PP,ZFMGST,YDMF_PHYS%TMP%APLPAR%PFL%FPLSH,YDMF_PHYS%TMP%APLPAR%PFL%FPLCH,&
         & YDMF_PHYS%OUT%GZ0,YDMF_PHYS_STATE%YGSP_RR%T,YDMF_PHYS%TMP%APLPAR%FLU%CDN,ZBLH,ZLMU,ZLMT,YDVARS%MXL%T0,ZLML,ZLMLTILD,ZRRCOR,LLMAF)
      ENDIF

    ELSE
      CLERR='APLPAR: UNEXPECTED VALUE FOR CGMIXLEN: '//TRIM(CGMIXLEN)
      CALL ABOR1(CLERR)
    ENDIF

    IF(LCOEFKTKE) THEN

      ! ------------------------------------------------------------- 
      ! COMPUTATION OF Ri', NCVPP AND COEFFICIENT FOR MOIST GUSTINESS
      ! ------------------------------------------------------------- 
      IF(LCOEFK_RIS) THEN
        CALL ACMRIS ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTCOEF,KLEV,&
         & YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS%TMP%APLPAR%MSC%LSCPE,&
         & ZQV, ZQL, ZQI, YDMF_PHYS%TMP%APLPAR%FLU%QSAT,&
         & YDMF_PHYS_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_STATE%T, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, ZLMU, ZLMT,YDMF_PHYS%OUT%GZ0,&
         & ZMRIPP)
      ENDIF

      CALL ACMRIP(YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTCOEF,KLEV,KSTEP,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,&
       & ZQV, ZQL, ZQI,ZCP,&
       & YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH,YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR,YDMF_PHYS%TMP%APLPAR%FLU%QSAT,YDMF_PHYS%TMP%APLPAR%MSC%QW,YDMF_PHYS%TMP%APLPAR%MSC%TW,&
       & YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,&
       & YDVARS%FQTUR%T0,YDVARS%FSTUR%T0,YDVARS%SHTUR%T0,&
       & YDMF_PHYS_STATE%TKE,YDVARS%TTE%T0,&
       & YDMF_PHYS_STATE%YGSP_RR%T,YDMF_PHYS_SURF%GSD_VH%PQSH,YDMF_PHYS%TMP%APLPAR%DSA%RS,YDMF_PHYS%TMP%APLPAR%DSA%CPS,ZRTI,YDMF_PHYS%OUT%GZ0,&
       & LLCLS,&
       & ZMRIPP,ZMRIFPP,&
       & ZBNEBCVPP,ZBNEBQ,ZNBVNO,&
       & ZFMTKE,ZFTTKE,ZF_EPS,ZFUN_TTE,&
       & ZAUTKE,ZATTKE,&
       & ZFHORM,ZFHORH,&
       & ZTH_FUN,ZWW_FUN,&
       & ZMRIMC,ZMRICTERM,ZMN2PP,&
       & ZMN2_ES,ZMN2_EQ,ZMN2_DS,ZMN2_DQ,ZFMGST)

    ENDIF ! LCOEFTKE

    ! FINISHING MIXING LENGTH COMPUTATION
    IF((CGMIXLEN(1:2) == 'EL').AND.LPTKE) THEN
      ZLMU(:,:)=ZLMU2(:,:)
      ZLMT(:,:)=ZLMT2(:,:)
      IF     (CGMIXLEN == 'EL0') THEN
        IMLTYPE=0
      ELSEIF (CGMIXLEN == 'EL1') THEN
        IMLTYPE=1
      ELSEIF (CGMIXLEN == 'EL2') THEN
        IMLTYPE=2
      ENDIF
      LLMAF=.FALSE.
      CALL ACMIXELEN(YGFL,YDPHY,YDPHY0,&
       & KIDIA,KFDIA,KLON,KTDIA,KLEV,KSTEP,IMLTYPE,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%T,ZQV,ZQL,ZQI,YDMF_PHYS_STATE%R,YDMF_PHYS_STATE%S,YDMF_PHYS_STATE%TKE,&
       & YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH,YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR,ZMN2PP,ZFMGST,YDMF_PHYS%TMP%APLPAR%PFL%FPLSH,YDMF_PHYS%TMP%APLPAR%PFL%FPLCH,&
       & YDMF_PHYS%OUT%GZ0,YDMF_PHYS_STATE%YGSP_RR%T,YDMF_PHYS%TMP%APLPAR%FLU%CDN,ZBLH,ZLMU,ZLMT,YDVARS%MXL%T0,ZLML,ZLMLTILD,ZRRCOR,LLMAF)
    ENDIF

  ENDIF ! (LVDIF or LGWD) and( not(LNEBR or LECT))

  IF ( LVDIF.OR.LHMTO.OR.LGWD ) THEN

    IF (LFLUSO.AND.(.NOT.LMSE)) THEN
      DO JLON=KIDIA,KFDIA
        ZCRTI(JLON) = 1.0_JPRB/(YDMF_PHYS_STATE%YGSP_RR%T(JLON)*YDMF_PHYS%TMP%APLPAR%DSA%RS(JLON))
      ENDDO
      CALL ACFLUSO ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,KLEV,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%Q,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,&
       & ZDPHIT,ZDPHIV,YDMF_PHYS%OUT%GZ0,YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS%TMP%APLPAR%FLU%QSATS,ZCRTI,YDMF_PHYS_STATE%YGSP_RR%T,&
       & LLHMT,YDMF_PHYS%TMP%APLPAR%FLU%CD,YDMF_PHYS%TMP%APLPAR%FLU%CDN,ZCDROV,ZCE,ZCEROV,YDMF_PHYS%TMP%APLPAR%FLU%CH,ZCHROV,&
       & YDMF_PHYS%OUT%QCLS,YDMF_PHYS%OUT%RHCLS,YDMF_PHYS%OUT%TCLS,YDMF_PHYS%OUT%UCLS,YDMF_PHYS%OUT%VCLS,YDMF_PHYS%OUT%UGST,YDMF_PHYS%OUT%VGST)
    ELSE
      DO JLON=KIDIA,KFDIA
        ZCE   (JLON) = YDMF_PHYS%TMP%APLPAR%FLU%CH   (JLON)
        ZCEROV(JLON) = ZCHROV(JLON)
      ENDDO
    ENDIF

  ENDIF ! LVDIF or LHMTO or LGWD

  IF (LRAFTKE) THEN
    YDMF_PHYS%OUT%CAPE(:)=0._JPRB
    ZDCAPE(:)=0._JPRB
    CALL ACCLDIA(YDXFU,YDPHY,YDPHY2,YDTOPH,KIDIA,KFDIA, KLON, KLEV,&
      & YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, YDMF_PHYS%OUT%CAPE, ZDCAPE, YDMF_PHYS_STATE%TKE, YDMF_PHYS_STATE%YCPG_DYN%PHIF, POROG, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST,ZBLH,KCLPH)
  ENDIF

!**
!     ------------------------------------------------------------------
!     6.- TURBULENCE: COEFFICIENTS D'ECHANGE
!     ------------------------------------------------------------------

  IF ( (LVDIF.OR.LGWD).AND.(.NOT.(LNEBR.OR.LECT)) ) THEN
     !-------------------------------------------------
     ! Compute diffusion coefficients.
     !-------------------------------------------------

    IF(LCOEFKTKE) THEN
      CALL ACTKECOEFK ( YDPHY0,YDPHY2,KIDIA,KFDIA,KLON,NTCOEF,KLEV,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%PREHYD,&
       & ZFMTKE,ZFTTKE,ZAUTKE,ZATTKE,&
       & YDMF_PHYS_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_STATE%T, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, ZLMU, ZLMT, YDMF_PHYS%OUT%GZ0,&
       & ZKTROV, ZKUROV, ZKNROV, ZXTROV, ZXUROV, ZXPTKEROV)
    ELSE
      CALL ACCOEFK ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTCOEF,KLEV,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS%TMP%APLPAR%MSC%LSCPE, ZQV, ZQL, ZQI, YDMF_PHYS%TMP%APLPAR%FLU%QSAT,&
       & YDMF_PHYS_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_STATE%T, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, YDMF_PHYS%TMP%APLPAR%PFL%FPLSH, YDMF_PHYS%TMP%APLPAR%PFL%FPLCH, ZLMU, ZLMT, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, ZBLH,&
       & ZKTROV, ZKUROV, ZKNROV, ZNBVNO, ZXTROV, ZXUROV, ZXPTKEROV)
    ENDIF
  ENDIF

  IF (LCVPPKF) THEN
    CALL ACVPPKF( YDMODEL%YRML_PHY_MF,KIDIA, KFDIA, KLON, NTCVIM, KLEV,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_STATE%T, ZQV,&
     & ZQL, ZQI, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, YDCPG_DYN0%CTY%VVEL(:,1:), YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_STATE%TKE,&
     & ZDIFCVPPQ, ZDIFCVPPS, ZCONDCVPPL, ZCONDCVPPI,&
     & ZPRODTH_CVPP, INLAB_CVPP, ZQLI_CVPP, ZNEB_CVPP, INND)
  ENDIF

     !-------------------------------------------------
     !  Call to EDKF
     !-------------------------------------------------
  IF (LEDKF) THEN
    IF (LEDMFI) THEN
       ZIMPL=0._JPRB
    ELSE  
       ZIMPL=1._JPRB 
    ENDIF
      
    IF (KSTEP == 0) YDMF_PHYS_SURF%GSD_SFL%PGROUP(:,:) = 0.0_JPRB
    CALL ARP_SHALLOW_MF( KIDIA,KFDIA,KLON,KTDIA,KLEV,ZIMPL,TSPHY,YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, &
                      & CMF_UPDRAFT,CMF_CLOUD,LMIXUV, &
                      & YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,YDMF_PHYS_STATE%T,ZQV,ZQL,ZQI,ZQR,ZQS,YDMF_PHYS_STATE%TKE,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,          &
                      & ZEDMFQ,ZEDMFS,ZEDMFU,ZEDMFV,                &
                      & YDMF_PHYS_SURF%GSD_SFL%PGROUP(:,1),YDMF_PHYS_SURF%GSD_SFL%PGROUP(:,2),ZPRODTH_CVPP,ZQLI_CVPP,ZNEB_CVPP,INLAB_CVPP,ZMF_UP)

  ENDIF

  IF ( LVDIF.AND.LECT ) THEN
    IF ( LCONDWT ) THEN
       YDCPG_MISC%QICE(:,:)= ZQI(:,:)
       YDCPG_MISC%QLI(:,:) = ZQL(:,:)
    ELSE
       YDCPG_MISC%QICE(:,:)= 0.0_JPRB
       YDCPG_MISC%QLI(:,:) = 0.0_JPRB
    ENDIF

! Computation of the 2 znlab used in acbl89
    IF (.NOT. LECSHAL) INLAB_CVPP(:,:) = 0
    IF (LECDEEP) THEN
       ZNLABCVP(:,:) = 1.0_JPRB
    ELSE
       ZNLABCVP(:,:) = 0.0_JPRB
    ENDIF
    IF(LNEBN.OR.LNEBR.OR.LRRGUST) THEN
      DO JLEV=KTDIA,KLEV
        DO JLON=KIDIA,KFDIA
          ZNLABCVP(JLON,JLEV) = ZNLABCVP(JLON,JLEV)&
         & *MAX(0.0_JPRB,SIGN(1.0_JPRB,YDMF_PHYS%TMP%APLPAR%PFL%FPLCH(JLON,JLEV)-YDMF_PHYS%TMP%APLPAR%PFL%FPLCH(JLON,JLEV-1)-ZEPS))
          ZNLAB(JLON,JLEV) = REAL(INLAB_CVPP(JLON,JLEV),JPRB)
        ENDDO
      ENDDO
    ENDIF
    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA
        ZNLAB(JLON,JLEV) = REAL(INLAB_CVPP(JLON,JLEV),JPRB)
      ENDDO
    ENDDO

    CALL ACTKE ( YDLDDH,YDMODEL%YRML_DIAG%YRMDDH,YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTCOEF,NTCOET,KLEV,&
     & YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,&
     & YDMF_PHYS_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_STATE%T, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, ZQV, ZQIC, ZQLC, YDMF_PHYS%TMP%APLPAR%MSC%LSCPE,&
     & YDMF_PHYS%TMP%APLPAR%FLU%CD, YDMF_PHYS%TMP%APLPAR%FLU%CH, YDMF_PHYS%OUT%GZ0, YDMF_PHYS_STATE%YGSP_RR%T, YDCPG_MISC%QS,&
     & YDCPG_MISC%QICE, YDCPG_MISC%QLI, YDMF_PHYS_STATE%TKE, ZPRODTH_CVPP, ZNLAB, ZNLABCVP,&
     & ZKTROV, ZKQROV, ZKQLROV, ZKUROV, ZXTROV, ZXUROV,&
     & ZNBVNO, ZNEBS, ZQLIS, ZNEBS0, ZQLIS0, ZCOEFN , YDMF_PHYS%TMP%APLPAR%PFL%FTKE, YDMF_PHYS%TMP%APLPAR%PFL%FTKEI, ZTKE1,ZTPRDY,&
     & PTDISS,YDDDH)
    YDMF_PHYS%OUT%CLPH(KIDIA:KFDIA)=MIN(XMAXLM,MAX(XMINLM,YDMF_PHYS%OUT%CLPH(KIDIA:KFDIA)))
  ENDIF


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
  IF (LMDUST.AND.(NGFL_EXT/=0)) THEN
! input dust scalar concentration in ppp from

  CALL ARO_MNHDUST (1,ILONMNH,KLEV,NGFL_EXT, PDT,&
                 & ZZI_SVM(KIDIA:KFDIA,:,:,1:NGFL_EXT),&
                 & ZZZ(KIDIA:KFDIA,:,:),&
                 & ZDZZ(KIDIA:KFDIA,:,:),&
                 & ZZI_PABSM(KIDIA:KFDIA,:,:),&
                 & ZZI_THM(KIDIA:KFDIA,:,:),&
                 & ZZI_RHODREFM(KIDIA:KFDIA,:,:),&
                 & NSWB_MNH,&
                 & KSTEP+1,&
                 & ZZI_SVM(KIDIA:KFDIA,:,:,1:NGFL_EXT),&
                 & ZPIZA_DST(KIDIA:KFDIA,:,:),&
                 & ZCGA_DST(KIDIA:KFDIA,:,:),&
                 & ZTAUREL_DST(KIDIA:KFDIA,:,:),&
                 & ZAERD(KIDIA:KFDIA,:),&
                 & NGFL_EZDIAG,&
                 & ZZI_PEZDIAG(KIDIA:KFDIA,:,:)           )

  PEZDIAG(KIDIA:KFDIA,:,:)=ZZI_PEZDIAG(KIDIA:KFDIA,:,:)

! return to aladin environment (inversion des niveaux)
   DO JGFL=1,NGFL_EXT
     DO JLON=1,KLON
       DO JLEV=1,KLEV
         ZSVM(JLON,JLEV,JGFL)=ZZI_SVM(JLON,1,JLEV,JGFL)
       ENDDO
     ENDDO
   ENDDO
  ENDIF

!      7.1 Albedo et emissivite en presence de neige
!          Albedo and emissivity with snow

  IF (.NOT.LMSE) THEN
!DEC$ IVDEP
    DO JLON=KIDIA,KFDIA
      IF (LSNV) THEN
        IF ((YDMF_PHYS_SURF%GSD_VF%PVEG(JLON) < 0.01_JPRB).OR.(YDMF_PHYS_SURF%GSD_VF%PALBF(JLON) >= 0.60_JPRB)) THEN
          ZALBV=0.0_JPRB
          YDMF_PHYS_SURF%GSD_VF%PALBSF(JLON)=YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)
        ELSE
          ZALBV=(YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)-(1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PVEG(JLON))*YDMF_PHYS_SURF%GSD_VF%PALBSF(JLON))/YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)
        ENDIF
        YDMF_PHYS%OUT%ALB(JLON)= (1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PVEG(JLON))*(1.0_JPRB-ZNEIJG(JLON)) *&
         & YDMF_PHYS_SURF%GSD_VF%PALBSF(JLON)&
         & + (1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PVEG(JLON))*ZNEIJG(JLON) *&
         & MAX(YDMF_PHYS_SURF%GSD_VF%PALBSF(JLON),YDMF_PHYS_STATE%YGSP_SG%A(JLON,1))&
         & + YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)*ZNEIJV(JLON) *&
         & MAX(ZALBV,YDMF_PHYS_STATE%YGSP_SG%A(JLON,1))&
         & + YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)*(1.0_JPRB-ZNEIJV(JLON)) * ZALBV
        YDMF_PHYS%TMP%APLPAR%FLU%EMIS(JLON)= (1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PVEG(JLON))*(1.0_JPRB-ZNEIJG(JLON)) *&
         & YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)&
         & + (1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PVEG(JLON))*ZNEIJG(JLON) * EMCRIN&
         & + YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)*ZNEIJV(JLON) * EMCRIN&
         & + YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)*(1.0_JPRB-ZNEIJV(JLON)) * YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)
      ELSE
        IF (LVGSN) THEN
          IF (LZ0HSREL.AND.LCOEFKSURF) THEN
            ! new treatment, PNEIJ is gridbox snow fraction
            YDMF_PHYS%OUT%ALB(JLON)=(1.0_JPRB-YDMF_PHYS%TMP%APLPAR%FLU%VEG(JLON)-YDMF_PHYS%TMP%APLPAR%FLU%NEIJ(JLON))*YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)+ &
             & YDMF_PHYS%TMP%APLPAR%FLU%VEG(JLON)*YDMF_PHYS_SURF%GSD_VV%PALV(JLON)+YDMF_PHYS%TMP%APLPAR%FLU%NEIJ(JLON)*YDMF_PHYS_STATE%YGSP_SG%A(JLON,1)
          ELSE
            ! old treatment, PNEIJ is snow fraction for bare ground
            YDMF_PHYS%OUT%ALB(JLON)=YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)-YDMF_PHYS%TMP%APLPAR%FLU%NEIJ(JLON)*(YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)- &
             & YDMF_PHYS_STATE%YGSP_SG%A(JLON,1))+(YDMF_PHYS%TMP%APLPAR%FLU%NEIJ(JLON)-ZNEIJV(JLON))*     &
             & YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)*(YDMF_PHYS_SURF%GSD_VV%PALV(JLON)-YDMF_PHYS_STATE%YGSP_SG%A(JLON,1))
          ENDIF

          YDMF_PHYS%OUT%ALB(JLON)=MIN(ABS(YDMF_PHYS_SURF%GSD_VV%PIVEG(JLON)-2._JPRB),1.0_JPRB) * YDMF_PHYS%OUT%ALB(JLON) +(&
           & 1.0_JPRB-MIN(ABS(YDMF_PHYS_SURF%GSD_VV%PIVEG(JLON)-2._JPRB),1.0_JPRB))&
           & * MAX(ALCRIN,YDMF_PHYS%OUT%ALB(JLON))
          YDMF_PHYS_SURF%GSP_SG%PT_T1(JLON,1)=YDMF_PHYS%OUT%ALB(JLON)

          YDMF_PHYS%TMP%APLPAR%FLU%EMIS(JLON)=YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)-YDMF_PHYS%TMP%APLPAR%FLU%NEIJ(JLON)*(YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)-EMCRIN)

        ELSE
          YDMF_PHYS%OUT%ALB(JLON)=YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)-YDMF_PHYS%TMP%APLPAR%FLU%NEIJ(JLON)*(YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)&
           & -MAX(YDMF_PHYS_SURF%GSD_VF%PALBF(JLON),ALCRIN))
          YDMF_PHYS%TMP%APLPAR%FLU%EMIS(JLON)=YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)-YDMF_PHYS%TMP%APLPAR%FLU%NEIJ(JLON)*(YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)-EMCRIN)
        ENDIF
      ENDIF
    ENDDO

    IF (LRAYFM) THEN
      ! diffuse and direct (parallel) albedo in NSW solar intervals
      IF (LALBMERCLIM) THEN
        DO JSG=1,NSW
          DO JLON=KIDIA,KFDIA
            ZALBD(JLON,JSG)=YDMF_PHYS%OUT%ALB(JLON)
            ZALBPMER=(1.0_JPRB+&
             & 0.5_JPRB*YDMF_PHYS%TMP%RDG%MU0M(JLON)*(1.0_JPRB/YDMF_PHYS%OUT%ALB(JLON)-1.0_JPRB))/   (&
             & 1.0_JPRB+YDMF_PHYS%TMP%RDG%MU0M(JLON)*(1.0_JPRB/YDMF_PHYS%OUT%ALB(JLON)-1.0_JPRB))**2
            ZALBP(JLON,JSG)=YDMF_PHYS%OUT%ALB(JLON)*          YDMF_PHYS_SURF%GSD_VF%PLSM(JLON)+&
                          & ZALBPMER  *(1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PLSM(JLON))
          ENDDO
        ENDDO
      ELSE
!DEC$ IVDEP
        DO JLON=KIDIA,KFDIA
          ZMERL(JLON)=(1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PLSM(JLON))*(1.0_JPRB&
            & -MAX(0.0_JPRB,SIGN(1.0_JPRB,TMERGL-YDMF_PHYS_STATE%YGSP_RR%T(JLON))))
          YDMF_PHYS%TMP%APLPAR%FLU%EMIS(JLON)=YDMF_PHYS%TMP%APLPAR%FLU%EMIS(JLON)*YDMF_PHYS_SURF%GSD_VF%PLSM(JLON)+ZMERL(JLON)*EMMMER&
            & +(1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PLSM(JLON))*(1.0_JPRB-ZMERL(JLON))*EMMGLA
        ENDDO
        DO JSG=1,NSW
          DO JLON=KIDIA,KFDIA
            ZALBD(JLON,JSG)=(1.0_JPRB-ZMERL(JLON))*YDMF_PHYS%OUT%ALB(JLON)+&
                                    & ZMERL(JLON) *ALBMED
            ZALBP(JLON,JSG)=(1.0_JPRB-ZMERL(JLON))*YDMF_PHYS%OUT%ALB(JLON)+&
             & ZMERL(JLON) *&
             & MAX(0.037_JPRB/(1.1_JPRB*YDMF_PHYS%TMP%RDG%MU0(JLON)**1.4_JPRB+0.15_JPRB),ZEPS0)
          ENDDO
        ENDDO
      ENDIF
    ELSEIF (LRAY) THEN
      ! direct (parallel) albedo for ACRANEB/ACRANEB2, Geleyn's formula
      ! with given proportion of Lambertian scattering
      DO JLON=KIDIA,KFDIA
        IF ( YDMF_PHYS_SURF%GSD_VF%PLSM(JLON) < 0.5_JPRB .AND. YDMF_PHYS_STATE%YGSP_RR%T(JLON) >= TMERGL ) THEN
          ZLAMB=RLAMB_WATER  ! water surface (open sea)
        ELSE
          ZLAMB=RLAMB_SOLID  ! solid surface (frozen sea or land)
        ENDIF
        ZALBDIR(JLON)=ZLAMB*YDMF_PHYS%OUT%ALB(JLON)+(1._JPRB-ZLAMB)*(1._JPRB+&
         & 0.5_JPRB*YDMF_PHYS%TMP%RDG%MU0(JLON)*(1.0_JPRB/YDMF_PHYS%OUT%ALB(JLON)-1.0_JPRB))/   (&
         & 1.0_JPRB+YDMF_PHYS%TMP%RDG%MU0(JLON)*(1.0_JPRB/YDMF_PHYS%OUT%ALB(JLON)-1.0_JPRB))**2
      ENDDO
    ENDIF

  ENDIF  ! .NOT.LMSE

  ! Appel de la routine d'aerosols

  LLAERO=LAEROSEA.AND.LAEROLAN.AND.LAEROSOO.AND.LAERODES

  IF ( .not. LRAYFM) THEN
    NAER=0
    NRADFR=NSTOP+1
  ENDIF

  IF    (   (LRAYFM.AND.(MOD(KSTEP,NRADFR) == 0)) &
  & .OR.  ( (LRAY.OR.LRAYSP).AND.(.NOT.LRSTAER)) ) THEN

    IF (LLAERO) THEN
      DO JLON = KIDIA,KFDIA
        ZAESEA(JLON) = YDMF_PHYS_SURF%GSD_VA%PSEA(JLON)
        ZAELAN(JLON) = YDMF_PHYS_SURF%GSD_VA%PLAN(JLON)
        ZAESOO(JLON) = YDMF_PHYS_SURF%GSD_VA%PSOO(JLON)
        ZAEDES(JLON) = YDMF_PHYS_SURF%GSD_VA%PDES(JLON)
      ENDDO
    ELSE
      DO JLON = KIDIA,KFDIA
        ZAESEA(JLON) = 0.0_JPRB
        ZAELAN(JLON) = 0.0_JPRB
        ZAESOO(JLON) = 0.0_JPRB
        ZAEDES(JLON) = 0.0_JPRB
      ENDDO
    ENDIF
    IF (LAEROSUL) THEN
      DO JLON = KIDIA,KFDIA
        ZAESUL(JLON) = YDMF_PHYS_SURF%GSD_VA%PSUL(JLON)
      ENDDO
    ELSE
      DO JLON = KIDIA,KFDIA
        ZAESUL(JLON) = 0.0_JPRB
      ENDDO
    ENDIF
    IF (LAEROVOL) THEN
       DO JLON = KIDIA,KFDIA
        ZAEVOL(JLON) = YDMF_PHYS_SURF%GSD_VA%PVOL(JLON)
      ENDDO
    ELSE
      DO JLON = KIDIA,KFDIA
        ZAEVOL(JLON) = 0.0_JPRB
      ENDDO
    ENDIF

    IF ( ( (LRAYFM.AND.NAER /= 0) .OR.LRAY.OR.LRAYSP).AND.LLAERO )  THEN
      CALL RADAER ( YDMODEL%YRML_PHY_RAD%YREAERD,YDERAD,YDPHY, KIDIA , KFDIA , KLON  , KLEV,&
       & YDMF_PHYS_STATE%YCPG_PHY%PREHYD , YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%T    , YDMF_PHYS_STATE%YGSP_RR%T,&
       & ZAESEA, ZAELAN, ZAESOO, ZAEDES, ZAESUL, ZAEVOL,&
       & ZAER,ZAERINDS                                 )
    ENDIF

  ELSEIF ( (LRAY.OR.LRAYSP).AND.LRSTAER ) THEN

    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA
        ZAER(JLON,JLEV,1)=ZDAER(JLEV)
        ZAER(JLON,JLEV,2:6)=0._JPRB
      ENDDO
    ENDDO

  ENDIF ! FOR AEROSOLS

! We uses the extinction coefficient explicitely solved by ARO_MNHDUST
  IF (LMDUST.AND.(NGFL_EXT/=0)) THEN
         ZAER(KIDIA:KFDIA,:,3) = ZAERD(KIDIA:KFDIA,:)
  ENDIF


!      7.2 Flux radiatifs par ciel clair (Code Geleyn)
!          Clear sky radiative fluxes    (Geleyn's scheme)

  ! separate clearsky call is kept only for old ACRANEB; for ACRANEB2
  ! duplicit calculation of gaseous transmissions is avoided
  IF (LRAY.AND.NRAY == 1.AND.KNFRRC /= 0) THEN
    IF (MOD(KSTEP,KNFRRC) == 0) THEN
      CALL ACRANEB(YDRIP,YDMODEL%YRML_PHY_MF, &
       & KIDIA,KFDIA,KLON,NTRADI,KLEV,IJN,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDCPG_MISC%NEB,&
       & ZQV,ZQCO2,YDCPG_MISC%QICE,YDCPG_MISC%QLI,ZQO3,YDMF_PHYS_STATE%T,&
       & YDMF_PHYS%OUT%ALB,ZALBDIR,YDMF_PHYS%TMP%APLPAR%FLU%EMIS,YDMF_PHYS%TMP%RDG%MU0,PGEMU,PGELAM,YDMF_PHYS%TMP%RDG%MU0LU,YDMF_PHYS_STATE%YGSP_RR%T,&
       & YDMF_PHYS%OUT%FRSO,YDMF_PHYS%OUT%FRTH,&
       & ZFRSODS,YDMF_PHYS%OUT%FRSOPS,YDMF_PHYS%OUT%FRSOLU,YDMF_PHYS%OUT%FRTHDS,ZAER,&
       & ZMAK,ZMAN)
      DO JLON=KIDIA,KFDIA
        YDMF_PHYS%OUT%FRSOC(JLON,0)=YDMF_PHYS%OUT%FRSO(JLON,NTRADI-1,1)
        YDMF_PHYS%OUT%FRSOC(JLON,1)=YDMF_PHYS%OUT%FRSO(JLON,KLEV,1)
        YDMF_PHYS%OUT%FRTHC(JLON,0)=YDMF_PHYS%OUT%FRTH(JLON,NTRADI-1,1)
        YDMF_PHYS%OUT%FRTHC(JLON,1)=YDMF_PHYS%OUT%FRTH(JLON,KLEV,1)
      ENDDO
    ENDIF
  ENDIF

!      7.3 Nebulosite et Convection
!          Cloud cover and Convection
!      7.3.1 Shallow + Deep convection

  IF (LCVPGY) THEN
    ! Le schema de convection de J. F. Gueremy
    IF (LCONDWT) THEN
      DO JLEV=KTDIA,KLEV
        DO JLON=KIDIA,KFDIA
          ZQCL(JLON,JLEV)=ZQL(JLON,JLEV)
          ZQCI(JLON,JLEV)=ZQI(JLON,JLEV)
        ENDDO
      ENDDO
    ELSE
      DO JLEV=KTDIA,KLEV
        DO JLON=KIDIA,KFDIA
          ZQCL(JLON,JLEV)=0.0_JPRB
          ZQCI(JLON,JLEV)=0.0_JPRB
        ENDDO
      ENDDO
    ENDIF
    CALL ACCVIMPGY ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTCVIM,KLEV,&
     & YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,&
     & YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDMF_PHYS%TMP%APLPAR%MSC%LH,YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR,ZQV,ZQCI,ZQCL,ZQLIS,YDMF_PHYS%TMP%APLPAR%FLU%QSAT,YDMF_PHYS%TMP%APLPAR%MSC%QW,&
     & YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,YDMF_PHYS_STATE%T,YDMF_PHYS%TMP%APLPAR%MSC%TW, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V,&
     & YDMF_PHYS%TMP%APLPAR%DSA%CPS,PGM,YDMF_PHYS_STATE%YGSP_RR%T,&
     & YDMF_PHYS%OUT%DIFCQ,YDMF_PHYS%OUT%DIFCQL,YDMF_PHYS%OUT%DIFCQN,YDMF_PHYS%OUT%DIFCS,YDMF_PHYS%OUT%FCCQL,YDMF_PHYS%OUT%FCCQN,&
     & YDMF_PHYS%OUT%FPFPCL,YDMF_PHYS%OUT%FPFPCN,YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,&
     & ZFHMLTSC,ZFHEVPPC,ZFPEVPPC,&
     & YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,ZNEBC,ZQLIC,YDMF_PHYS%OUT%STRCU,YDMF_PHYS%OUT%STRCV,&
     & ICIS,INLAB,&
     & INND,&
     & YDMF_PHYS_STATE%CVV)

    DO JLEV = KTDIA, KLEV
      DO JLON = KIDIA, KFDIA
        YDMF_PHYS%OUT%DIFCS(JLON,JLEV) = YDMF_PHYS%OUT%DIFCS(JLON,JLEV) - ZFHEVPPC(JLON,JLEV)&
                                            & - ZFHMLTSC(JLON,JLEV)
        YDMF_PHYS%OUT%DIFCQ(JLON,JLEV) = YDMF_PHYS%OUT%DIFCQ(JLON,JLEV) + YDMF_PHYS%OUT%FPEVPCL(JLON,JLEV)&
                                            & + YDMF_PHYS%OUT%FPEVPCN(JLON,JLEV)
      ENDDO
    ENDDO
    YDMF_PHYS%OUT%FPEVPCL=0._JPRB
    YDMF_PHYS%OUT%FPEVPCN=0._JPRB
    IF (LGRAPRO) THEN
      YDMF_PHYS%OUT%FPEVPCG=0._JPRB
    ENDIF
    ! Prise en compte des nuages convectifs diagnostiques sortant d'ACMTUD
    IF (LNCVPGY) THEN
      DO JLEV=KTDIA,KLEV
        DO JLON=KIDIA,KFDIA
          ZNEBC0(JLON,JLEV)=ZNEBC(JLON,JLEV)
          ZQLI_CVP(JLON,JLEV)=ZQLIC(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

    ! Annulation possible des flux convectifs pour les eaux condensees.
    IF (.FALSE.) THEN
      DO JLEV=0,KLEV
        DO JLON=KIDIA,KFDIA
          YDMF_PHYS%OUT%DIFCQL(JLON,JLEV)=0.0_JPRB
          YDMF_PHYS%OUT%DIFCQN(JLON,JLEV)=0.0_JPRB
        ENDDO
      ENDDO
    ENDIF
  ENDIF !  (LCVPGY)

  IF ( LCONDWT.AND.(.NOT.LNEBECT)) THEN

    IF(LCVPRO.AND.LNEBCV) THEN
! convective cloudiness in case we need protection of convective cloud water.
      DO JLEV=KTDIA,KLEV
        DO JLON=KIDIA,KFDIA
          ZNEBCH(JLON,JLEV)=ZUNEBH(JLON,JLEV)
        ENDDO
      ENDDO
      IF (LCVCSD) THEN
        DO JLEV=KTDIA,KLEV
          DO JLON=KIDIA,KFDIA
            ZNEBC0(JLON,JLEV)=MAX(ZEPSNEB,&
                            & MIN(1._JPRB-ZEPSNEB,ZUNEBH(JLON,JLEV)))
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    CALL ACNEBCOND ( YDRIP,YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTPLUI,KLEV,LLREDPR,&
     & ZHUC,ZVETAF,YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDCPG_MISC%RH,ZBLH,ZQV,ZQI,ZQL,&
     & YDMF_PHYS%TMP%APLPAR%MSC%QW,YDMF_PHYS_STATE%T,ZNEBCH,PGM,YDMF_PHYS_STATE%YGSP_RR%T,&
     & ZQLIS,ZNEBS,ZRHCRI,ZRH,ZQSATS,ZICEFR1,ZQLIS0,ZNEBS0)
     
    IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'ZRHCRI',ZRHCRI,KLON,KLEV) 
    IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'ZQLIS0',ZQLIS0,KLON,KLEV)
    IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'ZQLIS',ZQLIS,KLON,KLEV)

    IF(LRKCDEV) THEN
! Rash-Kristiansson cloud water scheme - second part.
      DO JLEV=KTDIA,KLEV
        DO JLON=KIDIA,KFDIA
        ! analytical solution dRH/dCLOUD:                                                
          ZRHDFDA(JLON,JLEV)=2._JPRB*(1._JPRB - ZNEBS(JLON,JLEV))&
               & *(1.0_JPRB-ZRHCRI(JLON,JLEV))
        ENDDO
      ENDDO
    ENDIF


  ENDIF ! LCONDWT .AND. .NOT.LNEBECT

  !-------------------------------------------------
  ! PCMT convection scheme.
  !-------------------------------------------------
  IF(LGPCMT) THEN
    ZSMOOTRAC(1:INBTRA) = GCVTSMO ! as usual 

    IF(LEDMFI) THEN
      CALL ACPCMT(YDGEM,YDGEOMETRY%YRDIM,YDGEOMETRY%YREGEO,  YDLDDH,YDMODEL%YRML_DIAG%YRMDDH,YDRIP,YDMODEL%YRML_PHY_MF, &
      & KIDIA,KFDIA,KLON,NTCVIM,KLEV,INBTRA,&
      ! - INPUT  2D .
      & YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH,YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDCPG_DYN0%CTY%VVEL(:,1:),YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS%TMP%APLPAR%MSC%LH,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,&
      & YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR,&
      & YDMF_PHYS_STATE%Q,ZQI,ZQL,YDMF_PHYS_STATE%R,YDMF_PHYS_STATE%S,ZQLIS,YDMF_PHYS%TMP%APLPAR%FLU%QSAT,YDMF_PHYS%TMP%APLPAR%MSC%QW,&
      & YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,YDMF_PHYS%TMP%APLPAR%MSC%CVGQ,YDMF_PHYS_STATE%T,YDMF_PHYS%TMP%APLPAR%MSC%TW,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,YDMF_PHYS_STATE%TKE,ZTRA(:,:,INBTRA_DEP:INBTRA),ZSMOOTRAC(1:INBTRA),&
      & ZQLC, ZQIC, ZQRC, ZQSC, &
      ! - INPUT 1D .
      & PGM, YDMF_PHYS_STATE%YGSP_RR%T, YDCPG_MISC%QS, YDMF_PHYS%TMP%APLPAR%FLU%NEIJ, YDMF_PHYS_SURF%GSD_VF%PLSM, ZVETAF, &
      ! - OUTPUT 2D .
      & ZEDMFQ,YDMF_PHYS%OUT%DIFCQL,YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQLC, YDMF_PHYS%OUT%DIFCQIC,&
      & ZEDMFS,ZDIFCVTH,ZTMPPRODTH,&
      & YDMF_PHYS%OUT%FCCQL,YDMF_PHYS%OUT%FCCQN,ZEDMFU,ZEDMFV,ZSTRCTRA(:,:,INBTRA_DEP:INBTRA),&
      & YDMF_PHYS%OUT%FIMCC,YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,&
      & YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,ZMF_UP,ZMU,ZMD,&
      & YDMF_PHYS%OUT%FPFPCL,YDMF_PHYS%OUT%FPFPCN,&
      & YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC,&
      & YDMF_PHYS%OUT%FCNEGQLC,YDMF_PHYS%OUT%FCNEGQIC,YDMF_PHYS%OUT%FCNEGQRC,YDMF_PHYS%OUT%FCNEGQSC, &
      & INLAB,&
      & ZNEBC0,ZQLI_CVP,ZTU,ZQU,ZQC_DET_PCMT,ZCSGC,ZENTCH, &
      ! - OUTPUT 1D .
      & INND,YDMF_PHYS%OUT%CAPE,ZAIPCMT,ZALF_CAPE,ZALF_CVGQ,&
      ! - INPUT/OUTPUT 2D .
      & YDVARS%UAL%T0,YDVARS%UOM%T0,YDVARS%DAL%T0,YDVARS%DOM%T0,YDDDH)
    ELSE
      CALL ACPCMT(YDGEM,YDGEOMETRY%YRDIM,YDGEOMETRY%YREGEO,  YDLDDH,YDMODEL%YRML_DIAG%YRMDDH,YDRIP,YDMODEL%YRML_PHY_MF, &
      & KIDIA,KFDIA,KLON,NTCVIM,KLEV,INBTRA,&
      ! - INPUT  2D .
      & YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH,YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDCPG_DYN0%CTY%VVEL(:,1:),YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS%TMP%APLPAR%MSC%LH,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,&
      & YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR,&
      & YDMF_PHYS_STATE%Q,ZQI,ZQL,YDMF_PHYS_STATE%R,YDMF_PHYS_STATE%S,ZQLIS,YDMF_PHYS%TMP%APLPAR%FLU%QSAT,YDMF_PHYS%TMP%APLPAR%MSC%QW,&
      & YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,YDMF_PHYS%TMP%APLPAR%MSC%CVGQ,YDMF_PHYS_STATE%T,YDMF_PHYS%TMP%APLPAR%MSC%TW,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,YDMF_PHYS_STATE%TKE,ZTRA(:,:,INBTRA_DEP:INBTRA),ZSMOOTRAC(1:INBTRA),&
      & ZQLC, ZQIC, ZQRC, ZQSC, &
      ! - INPUT 1D .
      & PGM, YDMF_PHYS_STATE%YGSP_RR%T, YDCPG_MISC%QS, YDMF_PHYS%TMP%APLPAR%FLU%NEIJ, YDMF_PHYS_SURF%GSD_VF%PLSM, ZVETAF, &
      ! - OUTPUT 2D .
      & YDMF_PHYS%OUT%DIFCQ,YDMF_PHYS%OUT%DIFCQL,YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQLC, YDMF_PHYS%OUT%DIFCQIC,&
      & YDMF_PHYS%OUT%DIFCS,ZDIFCVTH,ZTMPPRODTH,&
      & YDMF_PHYS%OUT%FCCQL,YDMF_PHYS%OUT%FCCQN,YDMF_PHYS%OUT%STRCU,YDMF_PHYS%OUT%STRCV,ZSTRCTRA(:,:,INBTRA_DEP:INBTRA),&
      & YDMF_PHYS%OUT%FIMCC,YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,&
      & YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,ZMF_UP,ZMU,ZMD,&
      & YDMF_PHYS%OUT%FPFPCL,YDMF_PHYS%OUT%FPFPCN,&
      & YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC,&
      & YDMF_PHYS%OUT%FCNEGQLC,YDMF_PHYS%OUT%FCNEGQIC,YDMF_PHYS%OUT%FCNEGQRC,YDMF_PHYS%OUT%FCNEGQSC, &
      & INLAB,&
      & ZNEBC0,ZQLI_CVP,ZTU,ZQU,ZQC_DET_PCMT,ZCSGC,ZENTCH, &
      ! - OUTPUT 1D .
      & INND,YDMF_PHYS%OUT%CAPE,ZAIPCMT,ZALF_CAPE,ZALF_CVGQ,&
      ! - INPUT/OUTPUT 2D .
      & YDVARS%UAL%T0,YDVARS%UOM%T0,YDVARS%DAL%T0,YDVARS%DOM%T0,YDDDH)
    ENDIF

  ENDIF

!         Appel du calcul de nebulosite.

  IF(LNEBN) THEN
    CALL ACNEBN ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTNEBU,KLEV,&
     & YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,ZQV,ZQL,ZQI,YDMF_PHYS%TMP%APLPAR%FLU%QSAT,YDMF_PHYS_STATE%T,YDMF_PHYS%TMP%APLPAR%PFL%FPLCH,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,&
     & ZUNEBH,&
     & ZNEBS0,ZQLIS0,ZQLI_CVP,ZNEB_CVPP,ZQLI_CVPP,&
     & ZAIPCMT,YDCPG_MISC%NEB,ZNEBC0,YDCPG_MISC%QICE,YDCPG_MISC%QLI,&
     & ZHUC,ZVETAF)
    DO JLEV=KTDIA,KLEV
!DEC$ IVDEP
      DO JLON=KIDIA,KFDIA
        YDCPG_MISC%NEB(JLON,JLEV)=YDCPG_MISC%NEB(JLON,JLEV)*GAEPS
      ENDDO
    ENDDO
  ENDIF

  IF ( LNEBR ) THEN
    CALL ACNEBR ( YDERAD,YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTCOEF,NTNEBU,KLEV,&
     & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS%TMP%APLPAR%MSC%LSCPE,ZQV,&
     & YDMF_PHYS%TMP%APLPAR%FLU%QSAT,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,&
     & YDMF_PHYS%OUT%GZ0,YDMF_PHYS%TMP%APLPAR%PFL%FPLCH,YDMF_PHYS_SURF%GSD_VH%PPBLH,YDCPG_MISC%QS,YDMF_PHYS_STATE%YGSP_RR%T,&
     & ZKTROV,ZKUROV,ZNBVNO,YDCPG_MISC%NEB,ZNEBS,YDCPG_MISC%QICE,YDCPG_MISC%QLI,&
     & ZQLIS)
  ENDIF

!         Diagnostique de nebulosite partielle.
  CALL ACNPART(YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTNEBU,KLEV,&
   & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,ZDECRD,YDCPG_MISC%NEB,&
   & YDMF_PHYS%OUT%CLCH,YDMF_PHYS%OUT%CLCM,YDMF_PHYS%OUT%CLCL,YDCPG_MISC%CLCT,ZCLCT_RAD,&
   ! optional arguments (convective cloud cover)
   & PCLCC=YDMF_PHYS%OUT%CLCC,PNEBC=ZNEBC0,PTOPC=YDMF_PHYS%OUT%CTOP)

!     7.3.5 Computation of the equivalent coefficients for simplified
!           radiation scheme

  IF ( LRAYSP .AND.(KSTEP == 1).AND.LRCOEF ) THEN
    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA
        ZZNEB(JLON,JLEV)=0.0_JPRB
      ENDDO
    ENDDO

    CALL ACRADCOEF ( YDRCOEF,YDPHY3,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,ZZNEB,ZQV,ZQCO2,YDCPG_MISC%QICE,YDCPG_MISC%QLI,&
     & ZQO3,YDMF_PHYS_STATE%T,&
     & YDMF_PHYS%TMP%RDG%MMU0,YDCPG_MISC%DHSF,YDMF_PHYS%TMP%APLPAR%FLU%EMIS,&
     & ZAER,&
     & POMPAC,PAC,YDMF_PHYS%TMP%APLPAR%RDT%COR,&
     & YDMF_PHYS%TMP%APLPAR%RDT%RAB3C,YDMF_PHYS%TMP%APLPAR%RDT%RAB3N,YDMF_PHYS%TMP%APLPAR%RDT%RAB4C,YDMF_PHYS%TMP%APLPAR%RDT%RAB4N,YDMF_PHYS%TMP%APLPAR%RDT%RAB6C,YDMF_PHYS%TMP%APLPAR%RDT%RAB6N,&
     & YDMF_PHYS%TMP%APLPAR%RDT%RAT1C,YDMF_PHYS%TMP%APLPAR%RDT%RAT1N,YDMF_PHYS%TMP%APLPAR%RDT%RAT2C,YDMF_PHYS%TMP%APLPAR%RDT%RAT2N,YDMF_PHYS%TMP%APLPAR%RDT%RAT3C,YDMF_PHYS%TMP%APLPAR%RDT%RAT3N,&
     & YDMF_PHYS%TMP%APLPAR%RDT%RAT4C,YDMF_PHYS%TMP%APLPAR%RDT%RAT4N,YDMF_PHYS%TMP%APLPAR%RDT%RAT5C,YDMF_PHYS%TMP%APLPAR%RDT%RAT5N)
  ENDIF

!     7.3.6 Module chimique - Chemistry module

  IF (LCHEM_ARPCLIM) THEN ! at this stage call when ARPEGE-Climat chemistry only

     ! initialisation below needs to be refined later for more general use
     IFLDX  = 1_JPIM
     IFLDX2 = 1_JPIM
     ILEVX  = 1_JPIM
     ALLOCATE(ZSD_XA(KLON,KLEV,IFLDX), ZSD_X2(KLON,IFLDX2))
     ALLOCATE(INDCHEM(NCHEM),IGPLAT(KLON))
     ALLOCATE(ZCFLX(KLON,NCHEM), ZCFLXO(KLON,NCHEM), ZCHEMDV(KLON,0)) ! no species with dry deposition in ARPCLIM
     ALLOCATE(ZAEROP(KLON,KLEV,NACTAERO),ZTENC(KLON,KLEV,NCHEM))
     ALLOCATE(ZDELP(KLON,KLEV), ZWND(KLON),ZDUMMY1(KLON,KLEV), ZGELAT(KLON))
     ALLOCATE(ZNEEFLX(KLON),ZCHEM2AER(KLON,KLEV,4))
     INDCHEM(:)     = 1_JPIM ! we should have here the indexes of the chemical species in the YCHEM array, to be implemented later
     IGPLAT (:)     = 1_JPIM
     ZSD_XA (:,:,:) = 0._JPRB
     ZSD_X2 (:,:)   = 0._JPRB
     DO JLEV=KTDIA,KLEV
       DO JLON=KIDIA,KFDIA
         ZDELP(JLON,JLEV) = YDMF_PHYS_STATE%YCPG_PHY%PREHYD(JLON,JLEV) - YDMF_PHYS_STATE%YCPG_PHY%PREHYD(JLON,JLEV-1)
       ENDDO
     ENDDO
     ZWND  (:)      = 0._JPRB  ! not used in  ARPEGE-Climat chemistry
     ZDUMMY1 (:,:)  = 1.0E-18_JPRB
     ZGELAT(KIDIA:KFDIA) = ASIN(PGEMU(KIDIA:KFDIA))
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

  IF ( LRAY ) THEN
    SELECT CASE (NRAY)
      CASE(1)
        CALL ACRANEB(YDRIP,YDMODEL%YRML_PHY_MF,&
         & KIDIA,KFDIA,KLON,NTRADI,KLEV,IJN,&
         & YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDCPG_MISC%NEB,&
         & ZQV,ZQCO2,YDCPG_MISC%QICE,YDCPG_MISC%QLI,ZQO3,YDMF_PHYS_STATE%T,&
         & YDMF_PHYS%OUT%ALB,ZALBDIR,YDMF_PHYS%TMP%APLPAR%FLU%EMIS,YDMF_PHYS%TMP%RDG%MU0,PGEMU,PGELAM,YDMF_PHYS%TMP%RDG%MU0LU,YDMF_PHYS_STATE%YGSP_RR%T,&
         & YDMF_PHYS%OUT%FRSO,YDMF_PHYS%OUT%FRTH,&
         & ZFRSODS,YDMF_PHYS%OUT%FRSOPS,YDMF_PHYS%OUT%FRSOLU,YDMF_PHYS%OUT%FRTHDS,ZAER,&
         & ZMAK,ZMAN)

        ! update sunshine duration [s]
!DEC$ IVDEP
        DO JLON=KIDIA,KFDIA
          IF ( YDMF_PHYS%OUT%FRSOPS(JLON) > RSUNDUR*YDMF_PHYS%TMP%RDG%MU0(JLON) ) THEN
            YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)=YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)+TSTEP
          ENDIF
        ENDDO
      CASE(2)
        CALL ACRANEB2(YDERDI,YDRIP,YDMODEL%YRML_PHY_MF,&
         & KIDIA,KFDIA,KLON,NTRADI,KLEV,IJN,KSTEP,KNFRRC,&
         & YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDCPG_MISC%NEB,&
         & ZQV,ZQCO2,YDCPG_MISC%QICE,YDCPG_MISC%QLI,ZQO3,YDMF_PHYS_STATE%T,&
         & YDMF_PHYS%OUT%ALB,ZALBDIR,YDMF_PHYS%TMP%APLPAR%FLU%EMIS,PGELAM,PGEMU,YDMF_PHYS%TMP%RDG%MU0,YDMF_PHYS%TMP%RDG%MU0LU,YDMF_PHYS_STATE%YGSP_RR%T,ZDECRD,ZCLCT_RAD,&
         & YDMF_PHYS%OPT%GDEOSI,YDMF_PHYS%OPT%GUEOSI,YDMF_PHYS%OPT%GMU0,YDMF_PHYS%OPT%GMU0_MIN,YDMF_PHYS%OPT%GMU0_MAX,&
         & YDMF_PHYS%OPT%GDEOTI,YDMF_PHYS%OPT%GDEOTI2,YDMF_PHYS%OPT%GUEOTI,YDMF_PHYS%OPT%GUEOTI2,YDMF_PHYS%OPT%GEOLT,YDMF_PHYS%OPT%GEOXT,&
         & YDMF_PHYS%OPT%GRPROX,YDMF_PHYS%OPT%GMIXP,YDMF_PHYS%OPT%GFLUXC,YDMF_PHYS%OPT%GRSURF,YDMF_PHYS_SURF%GSD_VD%PSUND,&
         & YDMF_PHYS%OUT%FRSO,YDMF_PHYS%OUT%FRTH,&
         & YDMF_PHYS%OUT%FRSOC,YDMF_PHYS%OUT%FRTHC,ZFRSODS,YDMF_PHYS%OUT%FRSOPS,YDMF_PHYS%OUT%FRSOLU,YDMF_PHYS%OUT%FRTHDS,ZAER)
    ENDSELECT

    ! sum downward diffuse and direct solar radiation at surface
    DO JLON=KIDIA,KFDIA
      YDMF_PHYS%OUT%FRSODS(JLON)=ZFRSODS(JLON)+YDMF_PHYS%OUT%FRSOPS(JLON)
    ENDDO

!      7.5 Rayonnement Morcrette
!          Morcrette's radiation

  ELSEIF ( LRAYFM ) THEN

    LLCALLRAD=(MOD(KSTEP,NRADFR) == 0 )
!    IF (NCALLRAD==1) ! <== not yet
    IF (NCALLRAD==2) LLCALLRAD=(LLCALLRAD.AND.(KSTEP<=NSTOP-1))
!    IF (NCALLRAD==3) ! <== not yet

    ! ---- Intermittent call to radiation scheme
    IF (LLCALLRAD) THEN
      CALL RECMWF(YDGEOMETRY%YRDIMV,YDMODEL,     &
       &  KIDIA , KFDIA, KLON  , KLEV   ,         &
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
        ZTRSODIR(KIDIA:KFDIA,JSG)=PGPAR(KIDIA:KFDIA,MSWDIR+JSG-1)
        ZTRSODIF(KIDIA:KFDIA,JSG)=PGPAR(KIDIA:KFDIA,MSWDIF+JSG-1)
      ENDDO
      ENDIF
    ENDIF

    IF (LRAYLU)  YDMF_PHYS%OUT%FRSOLU(KIDIA:KFDIA)=YDMF_PHYS%RAD%RMOON(KIDIA:KFDIA)

    ! ---- Flux update and radiative heating rates
    CALL RADHEAT ( YDERAD,YDERDI,YDMODEL%YRML_PHY_MF, &
    &   KIDIA  , KFDIA  , KLON    , KLEV,&
    &   YDMF_PHYS_STATE%YCPG_PHY%PREHYD  , YDMF_PHYS%TMP%APLPAR%FLU%EMIS  , YDMF_PHYS%RAD%EMTD   , YDMF_PHYS%TMP%RDG%MU0, ZQV  ,&
    &   ZTENT  , YDMF_PHYS%RAD%TRSW   , ZTRSOD , YDMF_PHYS_STATE%YGSP_RR%T , TSPHY,&
    &   ZTRSODIR, ZTRSODIF, ZALBD , ZALBP,&
    &   YDMF_PHYS%OUT%FRSO  , YDMF_PHYS%OUT%FRTH  , YDMF_PHYS%OUT%FRSODS , YDMF_PHYS%OUT%FRTHDS,&
    &   ZCEMTR , ZCTRSO , YDMF_PHYS%OUT%FRSOC  , YDMF_PHYS%OUT%FRTHC,&
    &   ZSUDU  , ZSDUR  , ZDSRP  , ZSFSWDIR , ZSFSWDIF, YDMF_PHYS%OUT%FRSOPS, ZFRSODS, YDMF_PHYS%OUT%FRSOPT )


    ! ---- Take into account day duration depending on altitude.
    IF(LDAYD) CALL ACDAYD(YDRIP,KIDIA,KFDIA,KLON,KLEV,KTDIA,KSGST,PGEMU ,YDMF_PHYS%TMP%RDG%MU0,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDMF_PHYS%OUT%FRSO)

    ! ---- Correct solar absorption as a function of pmu0.
    IF(GRSO < 1._JPRB) CALL ACRSO(YDPHY0,KIDIA,KFDIA,KLON,KLEV,KTDIA,KSGST,PGEMU  ,YDMF_PHYS%TMP%RDG%MU0,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDMF_PHYS%OUT%FRSO)

    IF(.NOT.LMSE) THEN
      DO JLON=KIDIA,KFDIA
        ZALB(JLON)=0.0_JPRB
        DO JSG=1,NSW
          ZALB(JLON)=ZALB(JLON)+0.5_JPRB*(ZALBD(JLON,JSG)+ZALBP(JLON,JSG))
        ENDDO
        ZALB(JLON)=ZALB(JLON)/FLOAT(NSW)
        YDMF_PHYS%OUT%FRSODS(JLON)=YDMF_PHYS%OUT%FRSO(JLON,KLEV,1)/(1.0_JPRB-ZALB(JLON))
      ENDDO
    ENDIF
    ! Compute Sunshine Duration (in seconds)
!DEC$ IVDEP
    DO JLON=KIDIA,KFDIA
      IF(YDMF_PHYS%OUT%FRSODS(JLON) >= RSUNDUR) THEN
        YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)=YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)+1.0_JPRB*TSTEP
      ENDIF
    ENDDO

  ENDIF

  IF (LRAY.OR.LRAYFM ) THEN

    ! Direct normal irradiance with securities
    DO JLON = KIDIA, KFDIA
      YDMF_PHYS%OUT%FRSDNI(JLON)=YDMF_PHYS%OUT%FRSOPS(JLON)
      IF (YDMF_PHYS%TMP%RDG%MU0(JLON) > 3.0E-02_JPRB) THEN
        YDMF_PHYS%OUT%FRSDNI(JLON)=YDMF_PHYS%OUT%FRSOPS(JLON)/YDMF_PHYS%TMP%RDG%MU0(JLON)
      ENDIF
      YDMF_PHYS%OUT%FRSDNI(JLON)=MAX(0.0_JPRB,YDMF_PHYS%OUT%FRSDNI(JLON))
    ENDDO
  ENDIF

  IF (LRAY.OR.LRAYFM) THEN

    ! global normal irradiance
    DO JLON=KIDIA,KFDIA
      YDMF_PHYS%OUT%FRSGNI(JLON)=YDMF_PHYS%OUT%FRSDNI(JLON)+0.5_JPRB*(                         &
       & (1.0_JPRB+YDMF_PHYS%TMP%RDG%MU0(JLON))*(YDMF_PHYS%OUT%FRSODS(JLON)-YDMF_PHYS%OUT%FRSOPS(JLON)       )+ &
       & (1.0_JPRB-YDMF_PHYS%TMP%RDG%MU0(JLON))*(YDMF_PHYS%OUT%FRSODS(JLON)-YDMF_PHYS%OUT%FRSO  (JLON,KLEV,1)))
    ENDDO

    ! mean radiant temperature
    IF (LXMRT) THEN
      CALL MEAN_RAD_TEMP(KIDIA,KFDIA,KLON,                             &
       & YDMF_PHYS%TMP%RDG%MU0,YDMF_PHYS%OUT%FRSO(:,KLEV,1),YDMF_PHYS%OUT%FRSODS,YDMF_PHYS%OUT%FRSOPS,YDMF_PHYS%OUT%FRTH(:,KLEV,1),YDMF_PHYS%OUT%FRTHDS, &
       & YDMF_PHYS%OUT%MRT) 
    ENDIF

  ENDIF

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
    CALL ACVEG ( YDPHY,YDPHY1,KIDIA,KFDIA,KLON,KLEV,&
     & YDMF_PHYS%OUT%FRSO,ZQV,YDMF_PHYS%TMP%APLPAR%FLU%QSAT,YDMF_PHYS_STATE%T,&
     & YDMF_PHYS_SURF%GSD_VV%PD2,YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS_SURF%GSD_VV%PIVEG,YDMF_PHYS_SURF%GSD_VV%PLAI,YDMF_PHYS%TMP%APLPAR%FLU%NEIJ,YDMF_PHYS%TMP%APLPAR%FLU%VEG,YDMF_PHYS_SURF%GSD_VV%PRSMIN,&
     & ZCHROV,ZGWDCS,ZWFC,YDMF_PHYS_STATE%YGSP_RR%FC,ZWWILT,YDMF_PHYS_STATE%YGSP_SB%Q,YDMF_PHYS%TMP%APLPAR%FLU%QSATS,&
     & ZHQ,ZHTR,ZHU,YDMF_PHYS_SURF%GSD_VV%PHV,ZWLMX)
  ENDIF

!     ------------------------------------------------------------------
!     8.- DIFFUSION VERTICALE TURBULENTE
!     ----------------------------------
  IF ( LVDIF ) THEN

!   Sauvegarde temporaire de l'ancien acdifus pour les besoins du Climat
    IF ( LACDIFUS ) THEN

      IF(NDIFFNEB == 1) THEN
        ZNEBDIFF(:,:)=ZNEBS(:,:)
      ELSEIF(NDIFFNEB == 2) THEN
        ZNEBDIFF(:,:)=YDCPG_MISC%NEB(:,:)
      ELSEIF(NDIFFNEB == 3) THEN
        ZNEBDIFF(:,:)=ZNEBS(:,:)+(1.0_JPRB-ZNEBS(:,:))*ZNEBCH(:,:)
      ENDIF

      CALL ACDIFUS ( YDMCC,YGFL,YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTDIFU,KLEV,KCSS,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS%OUT%FRSO, ZKTROV, ZKUROV, ZKNROV,&
       & ZQV, YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS_STATE%T, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, ZXTROV, ZXUROV, ZXPTKEROV, YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH,&
       & YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_DYN%RCP%R, ZQL, ZQI, ZNEBDIFF, YDMF_PHYS_STATE%TKE, ZLMU, YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR,&
       & YDMF_PHYS%TMP%APLPAR%FLU%CD, YDMF_PHYS%TMP%APLPAR%FLU%CDN, ZCDROV, ZCHROV, ZCEROV, YDMF_PHYS%TMP%APLPAR%DSA%CPS, YDMF_PHYS%OUT%CT,ZDQSTS,&
       & YDMF_PHYS%TMP%APLPAR%FLU%EMIS, YDMF_PHYS_SURF%GSD_VH%PSPSH, ZHQ, ZHTR, ZHU, YDMF_PHYS_SURF%GSD_VV%PHV,&
       & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PIVEG, YDMF_PHYS%TMP%APLPAR%FLU%NEIJ, YDCPG_MISC%QS, YDMF_PHYS%TMP%APLPAR%FLU%QSATS, YDMF_PHYS_STATE%YGSP_SB%T, YDMF_PHYS_STATE%YGSP_RR%T, YDMF_PHYS%TMP%APLPAR%FLU%VEG,&
       & ZXDROV, ZXHROV, YDMF_PHYS_STATE%YGSP_RR%W, YDMF_PHYS_STATE%YGSP_RR%IC,&
       & YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV,&
       & YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTQN, PTENDPTKE,&
       & YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FCLL, YDMF_PHYS%OUT%FCLN, YDMF_PHYS%OUT%FCS, YDMF_PHYS%TMP%APLPAR%FLU%FEVI, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, YDMF_PHYS%OUT%FEVV,&
       & YDMF_PHYS%OUT%FTR, YDMF_PHYS%TMP%APLPAR%DSA%LHS, YDMF_PHYS%OUT%FRTH, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H)

    ELSE

      IF ( LPTKE ) THEN
        IF ( LMSE.AND.LCALLSFX ) THEN       
          IF (KSTEP == 0) THEN
            YDMF_PHYS%TMP%APLPAR%FLU%CD(:)=YDMF_PHYS%TMP%APLPAR%FLU%CDN(:)   ! very first approximation
          ELSE
            YDMF_PHYS%TMP%APLPAR%FLU%CD(:)=MAX(PGPAR(:,MCD),ZEPS0)
          ENDIF
        ENDIF
        CALL ACPTKE(YGFL,YDLDDH,YDMODEL%YRML_PHY_MF,&
          & KIDIA,KFDIA,KLON,KTDIA,KLEV,KSTEP,&
          & YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,ZKNROV,&
          & YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH,&
          & YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%TKE,ZLMU,YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR,ZF_EPS,ZFUN_TTE,&
          & ZMRIPP,ZMRIFPP,ZRHS,&
          & ZMN2_ES,ZMN2_EQ,ZRRCOR,YDVARS%FQTUR%T0,YDVARS%FSTUR%T0,YDVARS%SHTUR%T0,&
          & YDMF_PHYS%TMP%APLPAR%FLU%CD,YDMF_PHYS%OUT%GZ0,&
          & ZKUROV,ZKTROV,ZXPTKEROV,&
          & ZLMLTILD,YDVARS%MXL%T0,ZLML,&
          & ZTH_FUN,ZWW_FUN,&
          & PTENDPTKE,YDVARS%TTE%T0,PFTCNS)
      ENDIF

      CALL ACDIFV1 ( YGFL,YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
        & YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,ZCP,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,ZKTROV,ZKQROV,ZKUROV,&
        & ZQV,ZQL,ZQI,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,ZSVM,ZXTROV,ZXUROV,&
        & ZXDROV,ZXHROV,&
        & ZEDMFS,ZEDMFQ,ZEDMFU,ZEDMFV,ZMF_UP, &
        & ZXURO,ZXQRO,ZXTRO,&
        & ZCFAQ,ZCFAS,ZCFATH,ZCFAU,ZCFASV,&
        & ZCFBQ,ZCFBS,ZCFBTH,ZCFBU,ZCFBV,ZCFBSV,&
        & ZDSE,ZQT)   
        

      IF ( LMSE.AND.LCALLSFX ) THEN

        IF (LRAYFM) THEN
          ZCARDI=RCARDI
        ELSEIF (LRAY) THEN
          ZCARDI=QCO2
          DO JLON=KIDIA,KFDIA
            ZSFSWDIF(JLON,1)=ZFRSODS(JLON)
            ZSFSWDIR(JLON,1)=YDMF_PHYS%OUT%FRSOPS(JLON)
          ENDDO
        ENDIF

        DO JLON=KIDIA,KFDIA
          ZRHODREFM(JLON)=YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(JLON,KLEV)/(YDMF_PHYS_STATE%T(JLON,KLEV)*YDMF_PHYS_STATE%YCPG_DYN%RCP%R(JLON,KLEV))
          ZDEPTH_HEIGHT(JLON,:)=(YDMF_PHYS_STATE%YCPG_DYN%PHIF(JLON,:)-YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,KLEV))/RG
          ZZS(JLON)=YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,KLEV)/RG
        ENDDO

        IRR=2

        CALL ARO_GROUND_PARAM( KBL, KGPCOMP,&
          & KFDIA-KIDIA+1, KIDIA, KFDIA, KSTEP,&
          & IRR, NSW, NGFL_EXT,NDGUNG, NDGUXG, NDLUNG, NDLUXG,LSURFEX_KFROM,&
          & LMPA, CCOUPLING,LDXFUMSE, &
          & NINDAT, ZRHGMT,ZSTATI,RSOVR,RCODEC,RSIDEC, &
          & PINDX(KIDIA:KFDIA), PINDY(KIDIA:KFDIA), &
          & YDMF_PHYS_STATE%U(KIDIA:KFDIA,KLEV:KLEV),YDMF_PHYS_STATE%V(KIDIA:KFDIA,KLEV:KLEV), &
          & YDMF_PHYS_STATE%T(KIDIA:KFDIA,KLEV:KLEV),YDMF_PHYS_STATE%Q(KIDIA:KFDIA,KLEV:KLEV),&
          & ZSVM(KIDIA:KFDIA,KLEV:KLEV,1:NGFL_EXT),&
          & ZCARDI,&
          & ZRHODREFM(KIDIA:KFDIA),&
          & YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(KIDIA:KFDIA,KLEV:KLEV),YDMF_PHYS_STATE%YCPG_PHY%PREHYD(KIDIA:KFDIA,KLEV:KLEV),&
          & ZDTMSE,ZDEPTH_HEIGHT(KIDIA:KFDIA,KLEV),ZZS(KIDIA:KFDIA), XZSEPS,&
          & YDMF_PHYS%TMP%RDG%MU0(KIDIA:KFDIA),YDMF_PHYS%TMP%RDG%MU0N(KIDIA:KFDIA),PGELAM(KIDIA:KFDIA),&
          & PGEMU(KIDIA:KFDIA),XSW_BANDS,&
          & ZSRAIN(KIDIA:KFDIA),ZSSNOW(KIDIA:KFDIA),&
          & ZSGROUPEL(KIDIA:KFDIA),&
          & YDMF_PHYS%OUT%FRTHDS(KIDIA:KFDIA),ZSFSWDIF(KIDIA:KFDIA,1:NSW),&
          & ZSFSWDIR(KIDIA:KFDIA,1:NSW),&
          & ZCFAQ(KIDIA:KFDIA,KLEV:KLEV), ZCFATH(KIDIA:KFDIA,KLEV:KLEV),&
          & ZCFAU(KIDIA:KFDIA,KLEV:KLEV),&
          & ZCFBQ(KIDIA:KFDIA,KLEV:KLEV), ZCFBTH(KIDIA:KFDIA,KLEV:KLEV),&
          & ZCFBU(KIDIA:KFDIA,KLEV:KLEV),&
          & ZCFBV(KIDIA:KFDIA,KLEV:KLEV),&
          & YDMF_PHYS%OUT%FCS(KIDIA:KFDIA,1),ZFEV(KIDIA:KFDIA),&
          & ZSFSV(KIDIA:KFDIA,1:NGFL_EXT),ZSFCO2(KIDIA:KFDIA),&
          & ZFMDU(KIDIA:KFDIA),ZFMDV(KIDIA:KFDIA),&
          & ZALBP(KIDIA:KFDIA,1:NSW),ZALBD(KIDIA:KFDIA,1:NSW),&
          & YDMF_PHYS%TMP%APLPAR%FLU%EMIS(KIDIA:KFDIA),ZTSN(KIDIA:KFDIA),&
          & YDMF_PHYS%OUT%FRTH(KIDIA:KFDIA,KLEV,1))               !orographic shadowing

! Opposite water vapor flux
        ZFEVS(:)=-ZFEV(:)

        CALL ARO_GROUND_DIAG( KBL, KGPCOMP,&
          & KFDIA-KIDIA+1, KIDIA, KFDIA, KLEV, -1,&
          & NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM,&
          & ZZS(KIDIA:KFDIA),ZFEVS(KIDIA:KFDIA),&
          & YDMF_PHYS_STATE%U(KIDIA:KFDIA,1:KLEV), YDMF_PHYS_STATE%V(KIDIA:KFDIA,1:KLEV),&
          & ZDEPTH_HEIGHT(KIDIA:KFDIA,1:KLEV),&
          & YDMF_PHYS%OUT%FRTH(KIDIA:KFDIA,KLEV,1),YDMF_PHYS%OUT%FRSO(KIDIA:KFDIA,KLEV,1),&
          & PINDX(KIDIA:KFDIA), PINDY(KIDIA:KFDIA),&
          & YDCPG_MISC%QS(KIDIA:KFDIA), YDMF_PHYS%OUT%GZ0(KIDIA:KFDIA),&
          & YDMF_PHYS%OUT%GZ0H(KIDIA:KFDIA), YDMF_PHYS%OUT%TCLS (KIDIA:KFDIA),&
          & YDMF_PHYS%OUT%QCLS(KIDIA:KFDIA), YDMF_PHYS%OUT%RHCLS(KIDIA:KFDIA),&
          & YDMF_PHYS%OUT%UCLS(KIDIA:KFDIA), YDMF_PHYS%OUT%VCLS(KIDIA:KFDIA),&
          & YDMF_PHYS%OUT%NUCLS(KIDIA:KFDIA), YDMF_PHYS%OUT%NVCLS(KIDIA:KFDIA), &
          & YDMF_PHYS%OUT%FCLL(KIDIA:KFDIA,1), YDMF_PHYS%OUT%FCLN(KIDIA:KFDIA,1),&
          & YDMF_PHYS%OUT%FEVL(KIDIA:KFDIA,1), YDMF_PHYS%OUT%FEVN(KIDIA:KFDIA,1),&
          & ZSSO_STDEV(KIDIA:KFDIA), ZTWSNOW(KIDIA:KFDIA),&
          & ZBUDTH(KIDIA:KFDIA), ZBUDSO(KIDIA:KFDIA),&
          & ZFCLL(KIDIA:KFDIA),ZTOWNS(KIDIA:KFDIA),&
          & ZCD(KIDIA:KFDIA)                        )

        DO JLON=KIDIA,KFDIA
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
          DO JLEV = 0, KLEV
            DO JLON = KIDIA, KFDIA
              YDMF_PHYS%OUT%FRTH(JLON,JLEV,JSG)=YDMF_PHYS%OUT%FRTH(JLON,JLEV,JSG)+ZBUDTH(JLON)
            ENDDO
          ENDDO
        ENDDO

! calculation of variables for the old "ISBA" atmosphere scheme
        IF (.NOT. LELAM) THEN
          CALL ARO_GROUND_DIAG_2ISBA( KBL, KGPCOMP, &
            & KFDIA-KIDIA+1, KIDIA, KFDIA, &
            & NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM, &
            & PINDX(KIDIA:KFDIA), PINDY(KIDIA:KFDIA), &
            & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PARG, YDMF_PHYS_SURF%GSD_VV%PSAB, YDMF_PHYS_SURF%GSD_VV%PD2, ZTSN, ZTWSNOW, &
            & YDMF_PHYS_SURF%GSP_SB%PT_T1(KIDIA:KFDIA,1), YDMF_PHYS_SURF%GSP_RR%PW_T1(KIDIA:KFDIA), YDMF_PHYS_SURF%GSP_SB%PQ_T1(KIDIA:KFDIA,1), &
            & YDMF_PHYS_SURF%GSP_RR%PIC_T1(KIDIA:KFDIA), YDMF_PHYS_SURF%GSP_SB%PTL_T1(KIDIA:KFDIA,1), YDMF_PHYS_SURF%GSP_RR%PFC_T1(KIDIA:KFDIA), &
            & YDMF_PHYS_SURF%GSP_SG%PA_T1(KIDIA:KFDIA,1), YDMF_PHYS_SURF%GSP_SG%PR_T1(KIDIA:KFDIA,1), ZHV2 )
          DO JLON=KIDIA,KFDIA
            YDMF_PHYS_SURF%GSD_VV%PHV(JLON) = ZHV2(JLON)
          ENDDO
        ENDIF

        IF (LDIFCONS.AND..NOT.LNODIFQC) THEN
          CALL ACAA1 ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
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

        IF (.NOT. LSFORCS) THEN
          IF (LEDMFI) THEN
            ZCFBS_G(:,:) = ZCFBS(:,:)
            ZCFBQ_G(:,:) = ZCFBQ(:,:) 
            ZCFBU_G(:,:) = ZCFBU(:,:)
            ZCFBV_G(:,:) = ZCFBV(:,:)
          ENDIF   
          CALL ARP_GROUND_PARAM ( YDMCC,YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,KCSS,&
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

          DO JLON=KIDIA,KFDIA
            YDMF_PHYS%OUT%DIFTQ(JLON,KLEV)=ZDIFWQ(JLON)
            YDMF_PHYS%OUT%DIFTS(JLON,KLEV)=ZDIFWS(JLON)
          ENDDO

          IF (LEDMFI) THEN
            ZCFBQ(:,:)  = ZCFBQ_G(:,:)
            ZCFBS(:,:)  = ZCFBS_G(:,:) 
            ZCFBU(:,KLEV)  = ZCFBU_G(:,KLEV)
            ZCFBV(:,KLEV)  = ZCFBV_G(:,KLEV) 
          ENDIF   

        ENDIF   ! <== LSFORCS

      ELSEIF (.NOT. LCALLSFX) THEN
        YDMF_PHYS%OUT%FCS(:,:) = 0.0_JPRB
        ZFEV(:)   = 0.0_JPRB

      ENDIF  !LMSE.AND.LCALLSFX

      CALL ACDIFV2 ( YGFL,YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
        & YDMF_PHYS_STATE%YCPG_PHY%PREHYD,ZCFAQ,ZCFAS,ZCFAU,ZCFASV,ZCFBQ,ZCFBS,ZCFBU,ZCFBV,ZCFBSV,&
        & ZKTROV,ZKQROV,ZKQLROV,ZKUROV,ZDSE,ZQT,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,ZPOID,&
        & YDMF_PHYS_STATE%T,ZQL,ZQI,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,&
        & ZCOEFA,ZALPHA1,ZLVT,ZQICE,&
        & ZSFSV,YDMF_PHYS%OUT%FCS,ZFEV,ZFMDU,ZFMDV,ZTSN,ZXHROV,&
        & PDIFSV,YDMF_PHYS%OUT%DIFTQ,YDMF_PHYS%OUT%DIFTS,YDMF_PHYS%OUT%STRTU,YDMF_PHYS%OUT%STRTV,YDMF_PHYS%OUT%DIFTQL,YDMF_PHYS%OUT%DIFTQN,YDMF_PHYS%OUT%DIFCQ,YDMF_PHYS%OUT%DIFCS,YDMF_PHYS%OUT%STRCU,YDMF_PHYS%OUT%STRCV,YDVARS%SHTUR%T0)

      IF (LEDKF) THEN
        IF (LSFORCS) THEN
        DO JLON=KIDIA,KFDIA
           YDMF_PHYS_SURF%GSD_SFL%PGROUP(JLON,1) = YDMF_PHYS%OUT%FCS(JLON,1)
           YDMF_PHYS_SURF%GSD_SFL%PGROUP(JLON,2) = ZFEV(JLON)
        ENDDO
        ELSE
          DO JLON=KIDIA,KFDIA             
            YDMF_PHYS_SURF%GSD_SFL%PGROUP(JLON,1) = YDMF_PHYS%OUT%DIFTS(JLON,KLEV)
            YDMF_PHYS_SURF%GSD_SFL%PGROUP(JLON,2) = YDMF_PHYS%OUT%DIFTQ(JLON,KLEV)
          ENDDO
        ENDIF   
      ENDIF

      IF ( LMSE ) THEN
        DO JLON=KIDIA,KFDIA
          PGPAR(JLON,MVTS) = ZTSN (JLON)
          YDMF_PHYS_SURF%GSP_SG%PF_T1(JLON,1) = ZTWSNOW (JLON)
        ENDDO
      ENDIF

      ! First compute horizontal exchange coefficients for momentum:
      !  (there's mo TOMs contribution, thus has to be done at latest here)
      IF (L3DTURB) THEN
        CALL ACTKECOEFKH(YDRIP,YDMODEL%YRML_PHY_MF,YDGEOMETRY%YREGEO,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
          & YDMF_PHYS_STATE%TKE, PTENDPTKE, YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%T, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH,&
          & YDMF_PHYS_STATE%DIV, YDMF_PHYS_STATE%VOR, YDVARS%U%DL, YDVARS%V%DL, YDMF_PHYS_STATE%YCPG_PHY%W, YDMF_PHYS_STATE%YCPG_PHY%WL, YDMF_PHYS_STATE%YCPG_PHY%WM, YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR,&
          & ZLML, ZFHORM, ZFHORH, ZKUROV,&
          & ZRTI, YDMF_PHYS%TMP%APLPAR%FLU%CD, ZCDROV,&
          & LPTKE, YDMF_PHYS%TMP%APLPAR%KUR%KUROV_H, YDMF_PHYS%TMP%APLPAR%KUR%KTROV_H, ZRHS)
      ENDIF

      IF (LCOEFKTKE) THEN
        IF (NDIFFNEB == 1) THEN
          ZCOEFA(:,:) = ZBNEBQ(:,:)
        ELSEIF (NDIFFNEB == 4) THEN
          ZCOEFA(:,:) = ZBNEBCVPP(:,:)
        ENDIF
        CALL ACDIFV3 ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,KSTEP,&
          & YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,ZCP,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,ZKTROV,ZXTROV,&
          & YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,ZCHROV,ZXHROV,&
          & ZCOEFA,ZQV,ZQL,ZQI,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH,YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR,&
          & ZLML, ZTH_FUN,ZWW_FUN, ZF_EPS,&
          & YDMF_PHYS%OUT%DIFTS,YDMF_PHYS%OUT%DIFTQ,YDMF_PHYS_STATE%TKE,PTENDPTKE,&
          & ZMN2_ES,ZMN2_EQ,ZMN2_DS,ZMN2_DQ,&
          & YDMF_PHYS%OUT%DIFTQL,YDMF_PHYS%OUT%DIFTQN,&
          & ZDIFWQ,ZDIFWS,&
          & YDMF_PHYS%OUT%FCLL,YDMF_PHYS%OUT%FCLN,YDMF_PHYS%OUT%FCS,YDMF_PHYS%TMP%APLPAR%FLU%FEVI,YDMF_PHYS%OUT%FEVL,YDMF_PHYS%OUT%FEVN,&
          & YDMF_PHYS%OUT%GZ0,ZRTI,&
          & ZSC_FEVI,ZSC_FEVN,ZSC_FCLL,ZSC_FCLN,&
          & ZTSTAR,ZTSTAR2,ZTSTARQ,ZTSTAR2Q)

        !store fluxes and shear term
        !GFL fields are on full levels, fluxes on half levels
        IF (YFQTUR%LGP.AND.YFSTUR%LGP) THEN
          DO JLEV=KTDIA,KLEV
            DO JLON=KIDIA,KFDIA
              YDVARS%FQTUR%T0(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)+YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)&
               &                                 +YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)
              YDVARS%FSTUR%T0(JLON,JLEV)=YDMF_PHYS%OUT%DIFTS(JLON,JLEV)
            ENDDO
          ENDDO
        ENDIF
      ENDIF ! LCOEFKTKE

      ! Now the heat coefficient can be completed by TKE+ containing
      !  the TOMs contribution. 
      IF (L3DTURB) THEN
        CALL ACTKECOEFKH(YDRIP,YDMODEL%YRML_PHY_MF,YDGEOMETRY%YREGEO,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
          & YDMF_PHYS_STATE%TKE, PTENDPTKE, YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%T, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH,&
          & YDMF_PHYS_STATE%DIV, YDMF_PHYS_STATE%VOR, YDVARS%U%DL, YDVARS%V%DL, YDMF_PHYS_STATE%YCPG_PHY%W, YDMF_PHYS_STATE%YCPG_PHY%WL, YDMF_PHYS_STATE%YCPG_PHY%WM, YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR,&
          & ZLML, ZFHORM, ZFHORH, ZKUROV,&
          & ZRTI, YDMF_PHYS%TMP%APLPAR%FLU%CD, ZCDROV,&
          & LPTKE, YDMF_PHYS%TMP%APLPAR%KUR%KUROV_H, YDMF_PHYS%TMP%APLPAR%KUR%KTROV_H, ZRHS)
      ENDIF

    ENDIF

!-----------------------------------------------------------------------------
!   THE DEEP CONVECTION WILL SEE THE SHALLOW PART FROM KFB AS IT IS WITH LOUIS
!   SCHEME AND THE MODIFIED RI
!----------------------------------------------------------------------------
    IF (LCVPPKF.OR.(LEDKF .AND. .NOT. LEDMFI)) THEN
      DO JLEV=KTDIA,KLEV
        DO JLON=KIDIA,KFDIA
          YDMF_PHYS%OUT%DIFTQ   (JLON,JLEV) = YDMF_PHYS%OUT%DIFTQ(JLON,JLEV) + ZDIFCVPPQ(JLON,JLEV)
          YDMF_PHYS%OUT%DIFTS   (JLON,JLEV) = YDMF_PHYS%OUT%DIFTS(JLON,JLEV) + ZDIFCVPPS(JLON,JLEV)
          YDMF_PHYS%OUT%STRTU   (JLON,JLEV) = YDMF_PHYS%OUT%STRTU(JLON,JLEV) + ZDIFCVPPU(JLON,JLEV)
          YDMF_PHYS%OUT%STRTV   (JLON,JLEV) = YDMF_PHYS%OUT%STRTV(JLON,JLEV) + ZDIFCVPPV(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

    IF (L3MT) THEN

    ! ------------------------------------------------------------------
    ! UPDATE TEMPERATURE, LIQUID WATER AND ICE BY THE TURBULENT FLUXES
    ! INCLUDING CORRECTION OF NEGATIVE VALUES OF WATER SPECIES
    ! SAVE THE INCREMENTAL PART DUE TO THE TURBULENT DIFFUSION PROCESSES
    ! ------------------------------------------------------------------

!   setup of auxiliary variables for water vapour updating for RK condensation scheme

    IF(LRKCDEV) THEN
      ZDTURDIFF=1._JPRB
      IF (LCVGQD) ZDTURDIFF=0._JPRB
    ENDIF

!cdir unroll=8
    DO JLEV = KTDIA,KLEV
!DEC$ IVDEP
      DO JLON = KIDIA,KFDIA

        ZT(JLON,JLEV)=YDMF_PHYS_STATE%T(JLON,JLEV)-ZIPOI(JLON,JLEV)/YDMF_PHYS_STATE%YCPG_DYN%RCP%CP(JLON,JLEV)&
         & *(YDMF_PHYS%OUT%DIFTS(JLON,JLEV)-YDMF_PHYS%OUT%DIFTS(JLON,JLEV-1))
        ZQX0=ZQI(JLON,JLEV)
        ZQX1=ZQX0-ZIPOI(JLON,JLEV)*(YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)&
          & -YDMF_PHYS%OUT%DIFTQN(JLON,JLEV-1))
        ZQI(JLON,JLEV)=MAX(0.0_JPRB,ZQX1)
        ZDQI=MAX(0.0_JPRB,ZQX1)-ZQX1
        ZFCQING(JLON,JLEV)=ZFCQING(JLON,JLEV-1)-ZDQI*ZPOID(JLON,JLEV)
        YDMF_PHYS%OUT%FCQNNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNNG(JLON,JLEV)+ZFCQING(JLON,JLEV)
        ZQX0=ZQL(JLON,JLEV)
        ZQX1=ZQX0-ZIPOI(JLON,JLEV)*(YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)&
          & -YDMF_PHYS%OUT%DIFTQL(JLON,JLEV-1))
        ZQL(JLON,JLEV)=MAX(0.0_JPRB,ZQX1)
        ZDQL= MAX(0.0_JPRB,ZQX1)-ZQX1
        ZFCQLNG(JLON,JLEV)=ZFCQLNG(JLON,JLEV-1)-ZDQL*ZPOID(JLON,JLEV)
        YDMF_PHYS%OUT%FCQLNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQLNG(JLON,JLEV)+ZFCQLNG(JLON,JLEV)
        ZDQC=ZDQI+ZDQL

        ZQX0=ZQV(JLON,JLEV)
        ZQX1=ZQX0-ZIPOI(JLON,JLEV)*(YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)&
          & -YDMF_PHYS%OUT%DIFTQ(JLON,JLEV-1))
        ZQV0=ZQX1-ZIPOI(JLON,JLEV)*(0.0_JPRB&
          & -ZFCQVNG(JLON,JLEV-1)-ZFCQING(JLON,JLEV-1)&
          & -ZFCQLNG(JLON,JLEV-1))
        IF(LCVGQD) THEN
          ZQV(JLON,JLEV)=MAX(0.0_JPRB,ZQV0-ZDQC)
        ENDIF
        ZDQV=MAX(0.0_JPRB,ZQV0-ZDQC)-ZQX1
        ZDQVDIFF(JLON,JLEV)=ZDQV+ZQX1-ZQX0
        ZFCQVNG(JLON,JLEV)=ZFCQVNG(JLON,JLEV-1)-ZDQV*ZPOID(JLON,JLEV)
        YDMF_PHYS%OUT%FCQNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNG(JLON,JLEV)+ZFCQVNG(JLON,JLEV)
        IF(LRKCDEV.AND.KSTEP>0) THEN
          ZDTRAD(JLON,JLEV)=ZIPOI(JLON,JLEV)/YDMF_PHYS_STATE%YCPG_DYN%RCP%CP(JLON,JLEV)* (&
            & YDMF_PHYS%OUT%FRSO(JLON,JLEV-1,1)-YDMF_PHYS%OUT%FRSO(JLON,JLEV,1)&
            & +YDMF_PHYS%OUT%FRTH(JLON,JLEV-1,1)-YDMF_PHYS%OUT%FRTH(JLON,JLEV,1) )
          ZRKTTEND(JLON,JLEV)=(ZT(JLON,JLEV)-YDVARS%RKTH%T0(JLON,JLEV)&
             & +ZDTRAD(JLON,JLEV))/TSPHY
          ZRKQVTEND(JLON,JLEV)=(MAX(0._JPRB,ZQV(JLON,JLEV)&
             & +ZDQVDIFF(JLON,JLEV)*ZDTURDIFF)&
             & -MAX(0._JPRB,YDVARS%RKTQV%T0(JLON,JLEV)))/TSPHY
          ZRKQCTEND(JLON,JLEV)= (MAX(0._JPRB,ZQI(JLON,JLEV)+ZQL(JLON,JLEV))&
             & -MAX(0._JPRB,YDVARS%RKTQC%T0(JLON,JLEV)))/TSPHY
        ENDIF
      ENDDO
    ENDDO

    ELSEIF (LSTRAPRO) THEN

    ! ------------------------------------------------------------------
    ! UPDATE THE CORRECTION FLUXES FOR NEGATIVE VALUES OF WATER SPECIES
    ! SAVE THE INCREMENTAL PART DUE TO THE TURBULENT DIFFUSION PROCESSES
    ! ------------------------------------------------------------------
    DO JLEV = KTDIA, KLEV
      DO JLON = KIDIA,KFDIA

        ZT (JLON,JLEV)= YDMF_PHYS_STATE%T(JLON,JLEV)
        ZQX0=ZQI(JLON,JLEV)
        ZQX1=ZQX0-ZIPOI(JLON,JLEV)*(YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)&
          & -YDMF_PHYS%OUT%DIFTQN(JLON,JLEV-1))
        ZDQI= MAX(0.0_JPRB,ZQX1)-ZQX1
        ZFCQING(JLON,JLEV)=ZFCQING(JLON,JLEV-1)-ZDQI*ZPOID(JLON,JLEV)
        YDMF_PHYS%OUT%FCQNNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNNG(JLON,JLEV)+ZFCQING(JLON,JLEV)
        ZQX0=ZQL(JLON,JLEV)
        ZQX1=ZQX0-ZIPOI(JLON,JLEV)*(YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)&
          & -YDMF_PHYS%OUT%DIFTQL(JLON,JLEV-1))
        ZDQL= MAX(0.0_JPRB,ZQX1)-ZQX1
        ZFCQLNG(JLON,JLEV)=ZFCQLNG(JLON,JLEV-1)-ZDQL*ZPOID(JLON,JLEV)
        YDMF_PHYS%OUT%FCQLNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQLNG(JLON,JLEV)+ZFCQLNG(JLON,JLEV)
        ZDQC=ZDQI+ZDQL

        ZQX0=ZQV(JLON,JLEV)
        ZQX1=ZQX0-ZIPOI(JLON,JLEV)*(YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)&
          & -YDMF_PHYS%OUT%DIFTQ(JLON,JLEV-1))
        ZQV0=ZQX1-ZIPOI(JLON,JLEV)*(0.0_JPRB&
          & -ZFCQVNG(JLON,JLEV-1)-ZFCQING(JLON,JLEV-1)&
          & -ZFCQLNG(JLON,JLEV-1))
        ZDQV=MAX(0.0_JPRB,ZQV0-ZDQC)-ZQX1
        ZFCQVNG(JLON,JLEV)=ZFCQVNG(JLON,JLEV-1)-ZDQV*ZPOID(JLON,JLEV)
        YDMF_PHYS%OUT%FCQNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNG(JLON,JLEV)+ZFCQVNG(JLON,JLEV)

      ENDDO
    ENDDO

    ENDIF ! 3MT || LSTRAPRO

    ! ------------------------------------------------------
    ! DIAGNOSTIC DE LA HAUTEUR DE COUCHE LIMITE, SELON TROEN ET MAHRT,
    ! POUR USAGE AU PAS DE TEMPS SUIVANT.
    ! ------------------------------------------------------

    IF((TRIM(CGMIXLEN) == 'TM').OR.(TRIM(CGMIXLEN) == 'TMC')) THEN
      ZBLH(KIDIA:KFDIA)=MIN(XMAXLM,MAX(XMINLM,YDMF_PHYS_SURF%GSD_VH%PPBLH(KIDIA:KFDIA)))
      CALL ACPBLHTM ( YDPHY0,KIDIA,KFDIA,KLON,NTDIFU,KLEV,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,&
       & ZQV, YDMF_PHYS_STATE%T, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,&
       & YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV,&
       & ZPCLS, YDMF_PHYS%OUT%TCLS,YDMF_PHYS%OUT%QCLS,&
       & YDMF_PHYS%OUT%FCS,YDMF_PHYS%OUT%FCLL,YDMF_PHYS%OUT%FCLN,YDMF_PHYS%TMP%APLPAR%DSA%LHS, YDMF_PHYS%OUT%GZ0,&
       & KCLPH,ZBLH,ZUGST,ZVGST)
      YDMF_PHYS_SURF%GSD_VH%PPBLH(KIDIA:KFDIA)=ZBLH(KIDIA:KFDIA)
      YDMF_PHYS%OUT%CLPH(KIDIA:KFDIA)=ZBLH(KIDIA:KFDIA)
      IF (.NOT.LRAFTUR) THEN
        YDMF_PHYS%OUT%UGST(KIDIA:KFDIA)=ZUGST(KIDIA:KFDIA)
        YDMF_PHYS%OUT%VGST(KIDIA:KFDIA)=ZVGST(KIDIA:KFDIA)
      ENDIF
    ENDIF

    IF(LNEBCO.AND.LNEBR)THEN
      DO JLON = KIDIA,KFDIA
        ZPREN(JLON)=YDMF_PHYS_STATE%YCPG_PHY%PREHYD(JLON,KLEV)
      ENDDO
      CALL ACPBLH ( YDPHY0,KIDIA,KFDIA,KLON,NTDIFU,KLEV,&
       & YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,&
       & ZQV, YDMF_PHYS_STATE%T, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,&
       & YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV,&
       & ZPREN, YDMF_PHYS_STATE%YGSP_RR%T, YDCPG_MISC%QS,&
       & YDMF_PHYS%OUT%FCS,YDMF_PHYS%OUT%FCLL,YDMF_PHYS%OUT%FCLN,YDMF_PHYS%TMP%APLPAR%DSA%LHS,ZRTI,YDMF_PHYS_SURF%GSD_VH%PPBLH)
    ENDIF
    ! ------------------------------------------------------------------
    ! UPDATE PASSIFS SCALAIRS DUE TO THE TURBULENT DIFFUSION PROCESSES
    ! ------------------------------------------------------------------
    IF (LMDUST.AND.(NGFL_EXT/=0)) THEN
      DO JGFL=1,NGFL_EXT
        DO JLON=1,KLON
          DO JLEV=1,KLEV
            ZSVM(JLON,JLEV,JGFL)=ZSVM(JLON,JLEV,JGFL)-ZIPOI(JLON,JLEV)&
          & *(PDIFSV(JLON,JLEV,JGFL)-PDIFSV(JLON,JLEV-1,JGFL))
            ZSVM(JLON,JLEV,JGFL)=MAX(ZSVM(JLON,JLEV,JGFL),0.0_JPRB)
          ENDDO
        ENDDO
      ENDDO
    ENDIF  ! LMDUST
  ENDIF ! LVDIF

!     DIAGNOSTIC SUPPLEMENTAIRE FLUX DE RAYONNEMENT (PFRTHDS ET PFRSOPT)
!     ADDITIONAL DIAGNOSTICS OF RADIATIVE FLUXES (PFRTHDS AND PFRSOPT)

!DEC$ IVDEP
  DO JLON=KIDIA,KFDIA
    IF(.NOT.LMSE)YDMF_PHYS%OUT%FRTHDS(JLON)=YDMF_PHYS%OUT%FRTH(JLON,KLEV,1)/YDMF_PHYS%TMP%APLPAR%FLU%EMIS(JLON)+RSIGMA*YDMF_PHYS_STATE%YGSP_RR%T(JLON)**4
    IF(YDMF_PHYS%TMP%RDG%MU0(JLON) <= 0.0_JPRB) THEN


      YDMF_PHYS%OUT%FRSOPT(JLON)=0.0_JPRB
    ELSE
      YDMF_PHYS%OUT%FRSOPT(JLON)=RII0*YDMF_PHYS%TMP%RDG%MU0(JLON)
    ENDIF
  ENDDO

! ADDITIONAL DIAGNOSTIC OF THE DERIVATIVE OF THE NON SOLAR SURFACE
! HEAT FLUX WITH RESPECT TO SURFACE TEMPERATURE (PDERNSHF)

  IF(LMCC03)THEN
    CALL ACDNSHF(YDPHY,YDPHY1,KIDIA, KFDIA, KLON, KTDIA, KLEV,&
     & YDMF_PHYS%TMP%APLPAR%FLU%EMIS, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS%TMP%APLPAR%FLU%NEIJ, ZQV, YDCPG_MISC%QS, YDMF_PHYS_STATE%YGSP_RR%T, ZCHROV,  ZDQSTS,&
     & YDMF_PHYS%OUT%DRNSHF)
  ENDIF

!*
!     ------------------------------------------------------------------
!     9.- TRAINEE DES ONDES DE GRAVITE INDUITES PAR LE RELIEF
!     ------------------------------------------------------------------
  IF ( LGWD ) THEN
    CALL ACDRAG ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTDRAG,KLEV,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,&
     & YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP, ZNBVNO, YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V,&
     & PRCORI,YDMF_PHYS_SURF%GSD_VF%PGETRL, ZGWDCS,&
     & YDMF_PHYS_SURF%GSD_VF%PVRLAN , YDMF_PHYS_SURF%GSD_VF%PVRLDI,&
     & YDMF_PHYS%OUT%STRDU, YDMF_PHYS%OUT%STRDV,ZTRAJGWD)
  ENDIF
 ! SAVE FOR TL/NL COEFS FROM VERT. DIFF AND GWD

  IF (LTRAJPS.AND.(LVDIFSPNL.OR.LGWDSPNL)) THEN


     IF(.NOT.LVDIFSPNL) THEN
       ZKTROV_SAVE(:,:)=0.0_JPRB
       ZKUROV_SAVE(:,:)=0.0_JPRB
       ZCDROV_SAVE(:)=0.0_JPRB
       ZCHROV_SAVE(:)=0.0_JPRB
     ENDIF
     IF(.NOT. LGWDSPNL) THEN
       ZTRAJGWD(:,:)=0.0_JPRB
     ENDIF

     CALL WRPHTRAJTM_NL(YDGEOMETRY,YDSIMPHL,KIDIA,KFDIA,PTRAJ_PHYS,&
                & ZKTROV_SAVE,ZKUROV_SAVE,&
                & ZCDROV_SAVE,ZCHROV_SAVE,ZTRAJGWD,&
                & YDMF_PHYS_STATE%L,YDMF_PHYS_STATE%I,YDMF_PHYS_STATE%R,YDMF_PHYS_STATE%S,ZQLIS,ZNEBS)
  ENDIF

!     ------------------------------------------------------------------
!     10.- PRECIPITATIONS STRATIFORMES.
!     ---------------------------------

  IF ( LSTRA.AND.(.NOT.LSTRAPRO) ) THEN
    CALL ACPLUIE ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTPLUI,KLEV,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,&
     & ZQV, YDMF_PHYS%TMP%APLPAR%MSC%QW, YDMF_PHYS_STATE%T,&
     & YDMF_PHYS%OUT%FPFPSL,YDMF_PHYS%OUT%FPFPSN,YDMF_PHYS%OUT%FCSQL,YDMF_PHYS%OUT%FCSQN,YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FPLSN)
  ENDIF

  IF ( LSTRAS ) THEN
    CALL ACPLUIS ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTPLUI,KLEV,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,ZNEBS,ZQV,ZQLIS,YDMF_PHYS%TMP%APLPAR%MSC%QW,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%T,&
     & YDMF_PHYS%OUT%FCSQL,YDMF_PHYS%OUT%FCSQN,YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FPLSN,YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPFPSL,YDMF_PHYS%OUT%FPFPSN)
  ENDIF

  IF (L3MT.OR.LSTRAPRO) THEN

    IF (LCOEFKTKE.OR.LLREDPR) THEN
      CALL ACTQSATS ( YDPHY,YDPHY0,KIDIA,KFDIA,KLON,NTQSAT,KLEV,&
       & YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, ZQV, ZT,&
       & ZGEOSLC, ZLH, ZLSCPE, ZQSAT, ZQW, ZRH, ZTW)
      CALL ACNEBCOND ( YDRIP,YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTPLUI,KLEV,LLREDPR,&
       & ZHUC,ZVETAF,YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDCPG_MISC%RH,ZBLH,ZQV,ZQI,ZQL,&
       & ZQW,ZT,ZNEBCH,PGM,YDMF_PHYS_STATE%YGSP_RR%T,&
       & ZQLIS,ZNEBS,ZRHCRI,ZRH,ZQSATS,ZICEFR1,ZQLIS0,ZNEBS0)
      CALL ACCDEV ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTPLUI,KLEV,&
       & ZVETAF,YDMF_PHYS_STATE%YCPG_DYN%PHI,&
       & YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,ZQV,ZQI,ZQL,ZQS,ZQR,ZQG,ZQW,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,ZT,YDMF_PHYS%TMP%APLPAR%MSC%TW,ZNEBS,&
       & ZRHCRI,ZICEFR1,ZQSATS,ZNEBCH,ZRHDFDA,ZRKTTEND,ZRKQVTEND,ZRKQCTEND,&
       & ZQLI_CVPP,ZNEB_CVPP,YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS%TMP%APLPAR%FLU%NEIJ,PGM,YDMF_PHYS_STATE%YGSP_RR%T,ZDECRD_MF,&
       & YDMF_PHYS%OUT%FPFPSL,YDMF_PHYS%OUT%FPFPSN,YDMF_PHYS%OUT%FPFPSG,YDMF_PHYS%OUT%FCSQL,YDMF_PHYS%OUT%FCSQN,YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPEVPSG,YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FPLSN,YDMF_PHYS%OUT%FPLSG,&
       & ZSEDIQL,ZSEDIQI,YDMF_PHYS%OUT%DIAGH)
    ELSE
      CALL ACCDEV ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTPLUI,KLEV,&
       & ZVETAF,YDMF_PHYS_STATE%YCPG_DYN%PHI,&
       & YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,ZQV,ZQI,ZQL,ZQS,ZQR,ZQG,YDMF_PHYS%TMP%APLPAR%MSC%QW,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%T,YDMF_PHYS%TMP%APLPAR%MSC%TW,ZNEBS,&
       & ZRHCRI,ZICEFR1,ZQSATS,ZNEBCH,ZRHDFDA,ZRKTTEND,ZRKQVTEND,ZRKQCTEND,&
       & ZQLI_CVPP,ZNEB_CVPP,YDMF_PHYS_SURF%GSD_VF%PLSM,YDMF_PHYS%TMP%APLPAR%FLU%NEIJ,PGM,YDMF_PHYS_STATE%YGSP_RR%T,ZDECRD_MF,&
       & YDMF_PHYS%OUT%FPFPSL,YDMF_PHYS%OUT%FPFPSN,YDMF_PHYS%OUT%FPFPSG,YDMF_PHYS%OUT%FCSQL,YDMF_PHYS%OUT%FCSQN,YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPEVPSG,YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FPLSN,YDMF_PHYS%OUT%FPLSG,&
       & ZSEDIQL,ZSEDIQI,YDMF_PHYS%OUT%DIAGH)
    ENDIF

  ENDIF

  IF (L3MT) THEN


    !  ----------------------------------------------------------------
    !   UPDATE HUMIDITY VARIABLES BY RESOLVED CONDENSATION FLUXES
    !  ----------------------------------------------------------------
!cdir unroll=8
    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA
        ZDQL=ZIPOI(JLON,JLEV)* (&
          & YDMF_PHYS%OUT%FCSQL(JLON,JLEV)-YDMF_PHYS%OUT%FCSQL(JLON,JLEV-1))
        ZQL(JLON,JLEV)=ZQL(JLON,JLEV)+ZDQL
        ZDQI=ZIPOI(JLON,JLEV)* (&
          & YDMF_PHYS%OUT%FCSQN(JLON,JLEV)-YDMF_PHYS%OUT%FCSQN(JLON,JLEV-1))
        ZQI(JLON,JLEV)=ZQI(JLON,JLEV)+ZDQI
        ZQV(JLON,JLEV)=ZQV(JLON,JLEV)-ZDQI-ZDQL

        ZT(JLON,JLEV)=ZT(JLON,JLEV)+(ZLHV(JLON,JLEV)*ZDQL&
                   & +ZLHS(JLON,JLEV)*ZDQI)/YDMF_PHYS_STATE%YCPG_DYN%RCP%CP(JLON,JLEV)
      ENDDO
    ENDDO

!rb : Store values for next time-step: all cases
    IF(LRKCDEV) THEN
      DO JLEV=KTDIA,KLEV
        DO JLON=KIDIA,KFDIA
          YDVARS%RKTH%T0(JLON,JLEV)= ZT(JLON,JLEV)+ZDTRAD(JLON,JLEV)
          YDVARS%RKTQV%T0(JLON,JLEV)=ZQV(JLON,JLEV)+ZDQVDIFF(JLON,JLEV)*ZDTURDIFF
          YDVARS%RKTQC%T0(JLON,JLEV)= ZQI(JLON,JLEV)+ZQL(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

  ENDIF ! L3MT

!     11.- PRECIPITATIONS SOUS-MAILLES.
!     ---------------------------------

  IF (LVDIF.AND.(LCVRA.OR.L3MT.OR.LCVRAV3)) THEN
!           LA TENDANCE DYNAMIQUE DE Q EST MULTIPLIEE PAR UN
!           FACTEUR CORRECTIF FONCTION DE LA RESOLUTION LOCALE
!           DU MODELE PUIS AUGMENTEE DE LA CONTRIBUTION
!           DE L EVAPORATION DU SOL.
    IF(LCVCSD.OR..NOT.LCVGQD) THEN
      IF (LCOMOD) THEN
          DO JLON=KIDIA,KFDIA
            ZMOD(JLON)=1._JPRB/(1.0_JPRB+PGM(JLON)*TEQK)
         ENDDO
      ELSE
         ZMOD(KIDIA:KFDIA)=1._JPRB
      ENDIF
!cdir unroll=8
      DO JLEV = KTDIA, KLEV
        DO JLON = KIDIA, KFDIA
          YDMF_PHYS%TMP%APLPAR%MSC%CVGQ(JLON,JLEV) = YDMF_PHYS%TMP%APLPAR%MSC%CVGQ(JLON,JLEV)*ZMOD(JLON) &
           & -RG*YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP(JLON,JLEV)&
           & *(YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)-YDMF_PHYS%OUT%DIFTQ(JLON,JLEV-1)&
           & +ZFCQVNG(JLON,JLEV)-ZFCQVNG(JLON,JLEV-1))
        ENDDO
      ENDDO
    ENDIF
    IF (LCAPE) THEN
      IF (LCAMOD) THEN
         DO JLON=KIDIA,KFDIA
            ZTAUX(JLON)=RTCAPE*(PGM(JLON)*TEQC)
         ENDDO
      ELSE
         ZTAUX=RTCAPE
      ENDIF ! lcamod
    ENDIF
  ENDIF


  IF (L3MT) THEN
    ! - TEMPORAIRES
    ! ZSIGPC: CONVECTIVE FRACTION OF PRECIPITATION FLUX, USED A POSTERIORI 
    ! ZSIGP : PRECIPITATION MESH FRACTION
    ZSIGPC=0.0_JPRB
    ZSIGP=1.0_JPRB  ! used to limit/scale dd area
  IF (LCVPRO) THEN
  !  -------------------------
  !  UPDRAUGHT CONTRIBUTION
  !  -------------------------

    CALL ACTQSATS ( YDPHY,YDPHY0,KIDIA,KFDIA,KLON,NTQSAT,KLEV,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, ZQV, ZT,&
     & ZGEOSLC, ZLH, ZLSCPE, ZQSAT, ZQW, ZRH, ZTW)

    IF (.NOT.LCVCSD) THEN
      CALL ACCVUD ( YDMODEL%YRML_PHY_EC%YRECUCONVCA,YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
       & YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,&
       & YDMF_PHYS%TMP%APLPAR%MSC%CVGQ, YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR,&
       & ZQV, ZQI, ZQL, ZQSAT, ZQW,&
       & YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP, ZT, ZTW, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, YDMF_PHYS_STATE%VOR,&
       & YDVARS%DAL%T0,&
       & YDMF_PHYS_STATE%YGSP_RR%T, ZTAUX, PRCORI,PGM, ZDECRD_MF,&
       & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCS,&
       & YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, ZATSLC,&
       & INLAB, YDVARS%UNEBH%T0, ZDETFI, ZTU, ZQU, ZUU, ZVU,&
       & INND,YDMF_PHYS%OUT%CAPE,YDMF_PHYS%OUT%CUCONVCA,YDMF_PHYS%OUT%NLCONVCA,&
       & ZGEOSLC, YDVARS%UEN%T0, YDVARS%UAL%T0, YDVARS%UOM%T0)

    ! Initialize temperature correction
      ZTCORR(KIDIA:KFDIA,:)=ZTU(KIDIA:KFDIA,:)-ZT(KIDIA:KFDIA,:)
    ! Initialize Environment horizontal velocity
      ZU(KIDIA:KFDIA,:)=YDMF_PHYS_STATE%U(KIDIA:KFDIA,:)
      ZV(KIDIA:KFDIA,:)=YDMF_PHYS_STATE%V(KIDIA:KFDIA,:)
    ELSE
! PUT IN ZTCORR THE CONDENSATE USED IN TRIGGERING
      DO JLEV=KTDIA,KLEV
       DO JLON=KIDIA,KFDIA
          ZTCORR(JLON,JLEV)=MAX(0._JPRB,ZQI(JLON,JLEV)+ZQL(JLON,JLEV) -&
          &  MAX(0._JPRB,YDMF_PHYS_STATE%I(JLON,JLEV)+YDMF_PHYS_STATE%L(JLON,JLEV)) )
       ENDDO
      ENDDO
! PASS T9 CONDENSATES TO ACCSU
      CALL ACCSU ( YDMODEL%YRML_PHY_EC%YRECUCONVCA,YDMODEL%YRML_PHY_EC%YRECUMF,YDMODEL%YRML_PHY_MF, &
        & KIDIA,KFDIA,KLON,KTDIA,KLEV,&
        & YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,&
        & YDMF_PHYS%TMP%APLPAR%MSC%CVGQ, YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR,&
        & ZQV, ZTCORR, YDMF_PHYS_STATE%I, YDMF_PHYS_STATE%L, ZQSAT, ZQW,&
        & YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP, ZT, ZTW, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, YDMF_PHYS_STATE%VOR,&
        & YDCPG_DYN0%CTY%VVEL(:,1:), &
        & YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, &
        & YDMF_PHYS_STATE%YGSP_RR%T, ZTAUX, PRCORI,ZDECRD_MF,&
        & YDMF_PHYS_SURF%GSD_VK%PUDGRO, &
        & YDMF_PHYS%OUT%CUCONVCA,YDMF_PHYS%OUT%NLCONVCA,&
        & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCS,&
        & YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, ZATSLC,&
        & INLAB, YDVARS%UNEBH%T0, ZDETFI, ZTU, ZQU, ZUU, ZVU,&
        & INND,YDMF_PHYS%OUT%CAPE, &
        & ZGEOSLC, YDVARS%UEN%T0, YDVARS%UAL%T0, YDVARS%UOM%T0)

    ! Initialize temperature correction
      ZTCORR(:,:)=ZTU(:,:)
    ! Initialize Environment horizontal velocity
      ZU(KIDIA:KFDIA,:)=YDMF_PHYS_STATE%U(KIDIA:KFDIA,:)
      ZV(KIDIA:KFDIA,:)=YDMF_PHYS_STATE%V(KIDIA:KFDIA,:)
    ENDIF


    CALL ACUPU(YDPHY,YDPHY0,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
               & YDVARS%UAL%T0,YDVARS%UNEBH%T0,ZDETFI,&
               & ZPOID,ZIPOI,ZLHV,ZLHS,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,&
               & YDMF_PHYS%OUT%DIFCQ,YDMF_PHYS%OUT%DIFCQN,YDMF_PHYS%OUT%DIFCQL,YDMF_PHYS%OUT%DIFCS,&
               & YDMF_PHYS%OUT%FCCQL,YDMF_PHYS%OUT%FCCQN,YDMF_PHYS%OUT%FCSQL,YDMF_PHYS%OUT%FCSQN,&
! output  1d .
               & ZSIGP, ZSIGPC,&
!input/output  2d .
               & ZNEBS,ZQI,ZQL,ZQV,ZT,&
               & YDMF_PHYS%OUT%FCQNG,YDMF_PHYS%OUT%FCQNNG,YDMF_PHYS%OUT%FCQLNG)

!rb store "-B" state: no protection of Ncv
    IF(LRKCDEV.AND.(.NOT.LNEBCV)) THEN
      DO JLEV=KTDIA,KLEV
        DO JLON=KIDIA,KFDIA
          YDVARS%RKTH%T0(JLON,JLEV)=YDVARS%RKTH%T0(JLON,JLEV)-ZT(JLON,JLEV)
          YDVARS%RKTQV%T0(JLON,JLEV)=YDVARS%RKTQV%T0(JLON,JLEV)-ZQV(JLON,JLEV)
          YDVARS%RKTQC%T0(JLON,JLEV)=YDVARS%RKTQC%T0(JLON,JLEV)-ZQI(JLON,JLEV)&
        & -ZQL(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

  ENDIF !LCVPRO


  !  ---------------------------------------------------------
  !  MICROPHYSICS (AUTOCONVERSION, COLLECTION AND EVAPORATION)
  !        OF TOTAL CONDENSATE (RESOLVED + SUB-GRID)
  !  ---------------------------------------------------------

    ZAUXPRC(:)=0._JPRB
    ZRCVOTT(:,:)=0._JPRB

    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA

    !    COMPUTE TOTAL CONDENSATION FLUX
    !    -------------------------------
        ZFCQL(JLON,JLEV)=YDMF_PHYS%OUT%FCSQL(JLON,JLEV)+YDMF_PHYS%OUT%FCCQL(JLON,JLEV)
        ZFCQI(JLON,JLEV)=YDMF_PHYS%OUT%FCSQN(JLON,JLEV)+YDMF_PHYS%OUT%FCCQN(JLON,JLEV)
        IF (LRCVOTT) THEN
          ZCONVC=MAX(0.0_JPRB, YDMF_PHYS%OUT%FCCQL(JLON,JLEV)+YDMF_PHYS%OUT%FCCQN(JLON,JLEV))
          ZAUXPRC(JLON)=ZAUXPRC(JLON)+&
         & MAX(0.0_JPRB,YDMF_PHYS%OUT%FCSQL(JLON,JLEV)+YDMF_PHYS%OUT%FCSQN(JLON,JLEV)&
         & -YDMF_PHYS%OUT%FCSQL(JLON,JLEV-1)-YDMF_PHYS%OUT%FCSQN(JLON,JLEV-1))
          ZTOTC=ZCONVC+ZAUXPRC(JLON)
          ZRCVOTT(JLON,JLEV)=MAX(0.0_JPRB,SIGN(1.0_JPRB,ZTOTC-ZEPS0))&
         & *ZCONVC/MAX(ZTOTC,ZEPS0)
        ENDIF
      ENDDO
    ENDDO

    CALL ACTQSATS ( YDPHY,YDPHY0,KIDIA,KFDIA,KLON,NTQSAT,KLEV,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, ZQV, ZT,&
     & ZGEOSLC, ZLH, ZLSCPE, ZQSAT, ZQW, ZRH, ZTW)

    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA
        ZDQ(JLON,JLEV)=ZQW(JLON,JLEV)&
         & -(ZQV(JLON,JLEV)+ZQI(JLON,JLEV)+ZQL(JLON,JLEV))
        ZDQM(JLON,JLEV)=ZQW(JLON,JLEV)
      ENDDO
    ENDDO

    CALL APLMPHYS( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, ZQV, ZQI, ZQL, ZQS, ZQR, ZQG, YDMF_PHYS_STATE%YCPG_DYN%RCP%R, ZT,&
     & ZIPOI, ZDQ, ZDQM, ZLHS, ZLHV, ZNEBS, ZPOID, ZRCVOTT,&
     & ZTCORR,YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS%TMP%APLPAR%FLU%NEIJ, YDMF_PHYS_STATE%YGSP_RR%T, ZDECRD_MF,&
     & ZFCQL, ZFCQI,&
     & ZFALLR, ZFALLS,ZFALLG,&
     & YDMF_PHYS%OUT%FPFPSL,YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN,YDMF_PHYS%OUT%FPLSG,&
     & ZSEDIQL,ZSEDIQI,ZMELNET,ZMELGET,YDMF_PHYS%OUT%DIAGH)

    !  ------------------------------------------------------------
    !   UPDATE AFTER MICROPHYSICS - 3MT
    !  ------------------------------------------------------------

   CALL ACUPM( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
              & YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, ZLHS, ZLHV,&
              & YDCPG_DYN0%CTY%EVEL, ZIPOI, ZPOID, YDVARS%UAL%T0, YDVARS%UOM%T0,&
              & YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPEVPSG,YDMF_PHYS%OUT%FPFPSL,YDMF_PHYS%OUT%FPFPSN,YDMF_PHYS%OUT%FPFPSG,&
              & YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FPLSN,YDMF_PHYS%OUT%FPLSG,&
              & ZFCQL,ZFCQI,YDMF_PHYS%OUT%FCCQL,YDMF_PHYS%OUT%FCCQN,&
! output 2D
              & ZOME, ZFHP,&
! input/output  2d .
              & ZQV,ZQL,ZQI,ZQR,ZQS,ZQG,ZT,&
              & YDMF_PHYS%OUT%FCSQL,YDMF_PHYS%OUT%FCSQN,&
              & YDMF_PHYS%OUT%FCQNG,YDMF_PHYS%OUT%FCQRNG,YDMF_PHYS%OUT%FCQSNG,YDMF_PHYS%OUT%FCQGNG)

  IF (LCDDPRO) THEN
  !  -------------------------
  !  DOWNDRAUGHT CONTRIBUTION
  !  -------------------------

    ZDIFCQD(:,:)=0.0_JPRB
    ZDIFCQLD(:,:)=0.0_JPRB
    ZDIFCQID(:,:)=0.0_JPRB
    ZDIFCSD(:,:)=0.0_JPRB
    ZSTRCUD(:,:)=0.0_JPRB
    ZSTRCVD(:,:)=0.0_JPRB

    CALL ACTQSATS ( YDPHY,YDPHY0,KIDIA,KFDIA,KLON,NTQSAT,KLEV,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, ZQV, ZT,&
     & ZGEOSLC, ZLH, ZLSCPE, ZQSAT, ZQW, ZRH, ZTW)

    IF (LNSDO) THEN
      CALL ACNSDO(YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
       & YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,&
       & YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR,ZQV,ZQI,ZQL,ZQR,ZQS,ZQW,&
       & YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP, ZSIGP,&
       & ZT,ZTW,ZU,ZV,YDCPG_DYN0%CTY%VVEL(:,1:),&
       & ZATSLC,ZGEOSLC,&
       & YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FPLSN,YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN,&
       & ZDIFCQD,ZDIFCQLD,ZDIFCQID,ZDIFCSD,&
       & YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,&
       & ZSTRCUD,ZSTRCVD,&
       & YDVARS%DAL%T0,YDVARS%DOM%T0)
    ELSE
      CALL ACMODO(YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
       & YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,&
       & YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR,ZQV,ZQI,ZQL,ZQW,&
       & YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,ZSIGP,ZT,ZTW,ZU,ZV,YDCPG_DYN0%CTY%EVEL,ZOME,&
       & ZATSLC,ZGEOSLC,ZFHP,&
       & YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FPLSN,YDMF_PHYS%OUT%FPLSG,&
       & ZDIFCQD,ZDIFCQLD,ZDIFCQID,ZDIFCSD,&
       & YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,YDMF_PHYS%OUT%FPEVPCG,&
       & ZSTRCUD,ZSTRCVD,&
       & YDVARS%DAL%T0,YDVARS%DOM%T0)
    ENDIF
  !  ---------------------------------------------
  !  UPDATE VARIABLES  BY DOWNDRAUGHT CONTRIBUTION
  !  ---------------------------------------------


    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA

    !   UPDATE CONVECTIVE DIFFUSION AND EVAPORATION FLUXES
    !   --------------------------------------------------

        YDMF_PHYS%OUT%DIFCS(JLON,JLEV) =YDMF_PHYS%OUT%DIFCS(JLON,JLEV) +ZDIFCSD(JLON,JLEV)
        YDMF_PHYS%OUT%DIFCQ(JLON,JLEV) =YDMF_PHYS%OUT%DIFCQ(JLON,JLEV) +ZDIFCQD(JLON,JLEV)
        YDMF_PHYS%OUT%DIFCQL(JLON,JLEV)=YDMF_PHYS%OUT%DIFCQL(JLON,JLEV)+ZDIFCQLD(JLON,JLEV)
        YDMF_PHYS%OUT%DIFCQN(JLON,JLEV)=YDMF_PHYS%OUT%DIFCQN(JLON,JLEV)+ZDIFCQID(JLON,JLEV)
        YDMF_PHYS%OUT%STRCU(JLON,JLEV)=YDMF_PHYS%OUT%STRCU(JLON,JLEV)+ZSTRCUD(JLON,JLEV)
        YDMF_PHYS%OUT%STRCV(JLON,JLEV)=YDMF_PHYS%OUT%STRCV(JLON,JLEV)+ZSTRCVD(JLON,JLEV)
      ENDDO
    ENDDO

    CALL ACUPD(YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
! input  2d .
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,ZLHS,ZLHV,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,&
     & ZIPOI,ZPOID,&
     & ZFALLR,ZFALLS,ZFALLG,&
     & YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPEVPSG,YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,YDMF_PHYS%OUT%FPEVPCG,ZFHP,&
     & ZDIFCSD,ZDIFCQD,ZDIFCQLD,ZDIFCQID,&
! input/output  2d .
     & ZQV,ZQL,ZQI,ZQR,ZQS,ZQG,ZT,YDVARS%UEN%T0,&
     & YDMF_PHYS%OUT%FCQNG,YDMF_PHYS%OUT%FCQNNG,YDMF_PHYS%OUT%FCQLNG,YDMF_PHYS%OUT%FCQRNG,YDMF_PHYS%OUT%FCQSNG,YDMF_PHYS%OUT%FCQGNG,&
     & YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FPLSN,YDMF_PHYS%OUT%FPLSG)

!rb  store "+C": case with no protection of Ncv
    IF((LRKCDEV).AND.(.NOT.LNEBCV)) THEN
      DO JLEV=KTDIA,KLEV
        DO JLON=KIDIA,KFDIA
          YDVARS%RKTH%T0(JLON,JLEV)=YDVARS%RKTH%T0(JLON,JLEV)+ZT(JLON,JLEV)
          YDVARS%RKTQV%T0(JLON,JLEV)=YDVARS%RKTQV%T0(JLON,JLEV)+ZQV(JLON,JLEV)
          YDVARS%RKTQC%T0(JLON,JLEV)=YDVARS%RKTQC%T0(JLON,JLEV)&
          & +ZQI(JLON,JLEV)+ZQL(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

  !  PARTITION CONVECTIVE/STRATIFORM PRECIPITATION
  ! ---------------------------------------------
    DO JLEV=KTDIA-1,KLEV
      DO JLON=KIDIA, KFDIA
         YDMF_PHYS%OUT%FPLCL(JLON,JLEV)=ZSIGPC(JLON)*YDMF_PHYS%OUT%FPLSL(JLON,JLEV)
         YDMF_PHYS%OUT%FPLSL(JLON,JLEV)=YDMF_PHYS%OUT%FPLSL(JLON,JLEV)-YDMF_PHYS%OUT%FPLCL(JLON,JLEV)
         YDMF_PHYS%OUT%FPLCN(JLON,JLEV)=ZSIGPC(JLON)*YDMF_PHYS%OUT%FPLSN(JLON,JLEV)
         YDMF_PHYS%OUT%FPLSN(JLON,JLEV)=YDMF_PHYS%OUT%FPLSN(JLON,JLEV)-YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
      ENDDO
    ENDDO
    IF (LGRAPRO) THEN
      DO JLEV=KTDIA-1,KLEV
        DO JLON=KIDIA, KFDIA
           YDMF_PHYS%OUT%FPLCG(JLON,KLEV)=0.0_JPRB !ZSIGPC(JLON)*PFPLSG(JLON,JLEV)
           YDMF_PHYS%OUT%FPLSG(JLON,JLEV)=YDMF_PHYS%OUT%FPLSG(JLON,JLEV)-YDMF_PHYS%OUT%FPLCG(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF
  ENDIF ! LCDDPRO

  ENDIF ! L3MT

! -----------------------------------------------------

  IF ( LCVRA ) THEN

    IF (LSTRAPRO) THEN
      DO JLEV = KTDIA-1, KLEV
        DO JLON = KIDIA, KFDIA
          YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)+ZFCQVNG(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

    CALL ACCVIMP ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTCVIM,KLEV,&
     & YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS%TMP%APLPAR%MSC%CVGQ, YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,&
     & YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR, ZQV, YDMF_PHYS%TMP%APLPAR%FLU%QSAT, YDMF_PHYS%TMP%APLPAR%MSC%QW, YDMF_PHYS_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS_STATE%T, YDMF_PHYS%TMP%APLPAR%MSC%TW, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, YDMF_PHYS_STATE%VOR, YDCPG_DYN0%CTY%VVEL(:,1:),&
     & PGM, YDMF_PHYS_STATE%YGSP_RR%T, ZTAUX, PRCORI,&
     & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN,&
     & YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, ZFPCOR, INLAB, YDMF_PHYS%OUT%CAPE, INND,&
     & YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS, ZGEOSLC, PGEMU)

    IF (LSTRAPRO) THEN
      DO JLEV = KTDIA-1, KLEV
        DO JLON = KIDIA, KFDIA
          YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)-ZFCQVNG(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

  ELSEIF (LCVRAV3) THEN

    CALL ACCVIMP_V3 (YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTCVIM,KLEV,&
     & YDMF_PHYS_STATE%YCPG_PHY%XYB%ALPH, YDMF_PHYS_STATE%YCPG_DYN%PHI, YDMF_PHYS_STATE%YCPG_DYN%PHIF, YDMF_PHYS_STATE%YCPG_PHY%PREHYD, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS%TMP%APLPAR%MSC%CVGQ, YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,&
     & YDMF_PHYS_STATE%YCPG_PHY%XYB%LNPR, ZQV, YDMF_PHYS%TMP%APLPAR%MSC%QW, YDMF_PHYS_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS_STATE%T, YDMF_PHYS%TMP%APLPAR%MSC%TW, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V,&
     & YDMF_PHYS_STATE%YGSP_RR%T, YDMF_PHYS_SURF%GSD_VH%PPBLH,&
     & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN,&
     & YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, INLAB, INND,&
     & YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS)

  ELSEIF (LCVTDK) THEN     ! <== IFS deep convection scheme

    DO JLEV=1,KLEV
      DO JLON=KIDIA,KFDIA
        ZGEOM1(JLON,JLEV)     = YDMF_PHYS_STATE%YCPG_DYN%PHIF(JLON,JLEV)-YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,KLEV)
        ZVERVEL(JLON,JLEV)    = YDCPG_DYN0%CTY%VVEL(JLON,JLEV)*YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(JLON,JLEV)
        ZLISUM(JLON,JLEV)     = 0.0_JPRB
        ZLCRIT_AER(JLON,JLEV) = 5.E-4_JPRB
      ENDDO
    ENDDO
    DO JLEV=0,KLEV
      DO JLON=KIDIA,KFDIA
        ZGEOMH(JLON,JLEV)=YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,JLEV)-YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,KLEV)
        ZFCQLF(JLON,JLEV)=0.0_JPRB
        ZFCQLI(JLON,JLEV)=0.0_JPRB
      ENDDO
    ENDDO
    DO JLON=KIDIA,KFDIA
       LLLAND(JLON)    = (YDMF_PHYS_SURF%GSD_VF%PLSM(JLON)>0.5_JPRB)
       LLSHCV(JLON)    = .FALSE.
       ZCUCONVCA(JLON) = 0.0_JPRB
       ZACPR(JLON)     = 0.0_JPRB
       ZDXTDK(JLON)    = RDELXN/PGM(JLON)
    ENDDO
    LLPTQ = .FALSE.
    CALL ACUPTQ ( KLON,KIDIA,KFDIA,KLEV,LLPTQ,YDMF_PHYS%OUT%FRSO,YDMF_PHYS%OUT%FRTH,&
     & YDMF_PHYS%OUT%DIFCQ,YDMF_PHYS%OUT%DIFCS,YDMF_PHYS%OUT%DIFTQ,YDMF_PHYS%OUT%DIFTS,YDMF_PHYS%OUT%FCCQL,YDMF_PHYS%OUT%FCCQN,YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,&
     & YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,YDMF_PHYS%OUT%FPFPCL,YDMF_PHYS%OUT%FPFPCN,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%Q,YDMF_PHYS_STATE%YGSP_RR%T,&
     & ZTENHA,ZTENQVA )
    
    DO JLEV=1,KLEV
      DO JLON=KIDIA,KFDIA 
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
     & KIDIA,KFDIA,KLON,0,KLEV,ZDXTDK,ISPPN2D,LLLAND,LLDSLPHY,TSPHY,ZVDIFTS,&
     & YDMF_PHYS_STATE%T,ZQV,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,ZLISUM,ZVERVEL,YDMF_PHYS%OUT%DIFTQ,YDMF_PHYS%OUT%DIFTS,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,ZGEOM1,ZGEOMH,PGM, &
     & ZCUCONVCA,ZGP2DSPP,ZTENT,ZTENQ,ZTENU,ZTENV,ZACPR,&
     & ITOPC,IBASC,ITYPE,ICBOT,ICTOP,IBOTSC,LLCUM,LLSC,LLSHCV,ZLCRIT_AER,&
     & ZLU,ZLUDE,ZLUDELI,ZSNDE,ZMFU,ZMFD,YDMF_PHYS%OUT%DIFCQ,YDMF_PHYS%OUT%DIFCS,ZFHPCL,ZFHPCN,&
     & YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,ZLRAIN,ZRSUD,YDMF_PHYS%OUT%STRCU,YDMF_PHYS%OUT%STRCV,ZFCQLF,ZFCQLI,&
     & ZMFUDE_RATE,ZMFDDE_RATE,YDMF_PHYS%OUT%CAPE,ZWMEAN,ZDIFF,0,ZCEN,ZTENC,ZSCAV)
    DO JLEV=0,KLEV
      DO JLON=KIDIA,KFDIA 
         ZFPCOR  (JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)+YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
         YDMF_PHYS%OUT%FCCQL  (JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)
         YDMF_PHYS%OUT%FCCQN  (JLON,JLEV)=YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
         YDMF_PHYS%OUT%FPFPCL (JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)
         YDMF_PHYS%OUT%FPFPCN (JLON,JLEV)=YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
         YDMF_PHYS%OUT%FPEVPCL(JLON,JLEV)=0._JPRB
         YDMF_PHYS%OUT%FPEVPCN(JLON,JLEV)=0._JPRB
       ENDDO
    ENDDO

  ENDIF
  IF(LFLASH .AND. LCVTDK) THEN
    ! LIGHTNING PARAMETERIZATION. OUTPUT IS PFLASH, TOTAL LIGHTNING FLASH RATES.
    ZGAW(:)=0._JPRB
    CALL CULIGHT &
      & ( YDEPHY,YGFL,KIDIA ,  KFDIA,   KLON,    KLEV,  ZGAW, ZGAW, &
      &   YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,     YDMF_PHYS_STATE%YCPG_PHY%PREHYD  ,  YDMF_PHYS_STATE%YCPG_DYN%PHI,   YDMF_PHYS_STATE%YCPG_DYN%PHIF,  LLLAND,&
      &   YDMF_PHYS_STATE%T    ,  ZLU   ,  ZMFU,    YDMF_PHYS%OUT%CAPE, &
      &   YDMF_PHYS%OUT%FPLCL,  YDMF_PHYS%OUT%FPLCN,&
      &   LLCUM ,  ICBOT ,  ICTOP,&
      ! Outputs
      &   LLLINOX, YDMF_PHYS%OUT%FLASH,  ZLIGH_CTG,  ZCTOPH, &
      &   ZPRECMX, ZICE,   ZCDEPTH,  ZWMFU) 
    ! LIGHTNING FLASH RATES ARE CONVERTED IN fl/km2/s BEFORE ENTERING CFU TIME ACCUMULATION.
    DO JLON=KIDIA,KFDIA
      YDMF_PHYS%OUT%FLASH(JLON)=YDMF_PHYS%OUT%FLASH(JLON)/86400._JPRB
    ENDDO
  ENDIF

  IF ( LADJCLD ) THEN
    ZFRSO(:,:) = 0.0_JPRB
    ZFRTH(:,:) = 0.0_JPRB
    LLPTQ = .TRUE.
    CALL ACUPTQ ( KLON,KIDIA,KFDIA,KLEV,LLPTQ,ZFRSO,ZFRTH,&
     & YDMF_PHYS%OUT%DIFCQ,YDMF_PHYS%OUT%DIFCS,YDMF_PHYS%OUT%DIFTQ,YDMF_PHYS%OUT%DIFTS,YDMF_PHYS%OUT%FCCQL,YDMF_PHYS%OUT%FCCQN,YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,&
     & YDMF_PHYS%OUT%FPEVPCL,YDMF_PHYS%OUT%FPEVPCN,YDMF_PHYS%OUT%FPFPCL,YDMF_PHYS%OUT%FPFPCN,YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%Q,YDMF_PHYS_STATE%YGSP_RR%T,&
     & ZTENHA,ZTENQVA )
  ENDIF

  IF ( LPROCLD ) THEN
    IF(LGPCMT) THEN
      ! Microphysics occurs in a shell between updraft and its resolved environment.
      ZNEBS(KIDIA:KFDIA,KTDIA:KLEV)=1._JPRB-(1._JPRB-ZNEBS(KIDIA:KFDIA,KTDIA:KLEV))&
        & *(1._JPRB-ZCSGC(KIDIA:KFDIA,KTDIA:KLEV))
      ZTMIC(KIDIA:KFDIA,KTDIA:KLEV)=YDMF_PHYS_STATE%T(KIDIA:KFDIA,KTDIA:KLEV)
      ZQMIC(KIDIA:KFDIA,KTDIA:KLEV)=ZCSGC(KIDIA:KFDIA,KTDIA:KLEV)&
        & *YDMF_PHYS%TMP%APLPAR%FLU%QSAT(KIDIA:KFDIA,KTDIA:KLEV)&
        & +(1._JPRB-ZCSGC(KIDIA:KFDIA,KTDIA:KLEV))*ZQV(KIDIA:KFDIA,KTDIA:KLEV)
      ZQLIS(KIDIA:KFDIA,KTDIA:KLEV)=MAX(ZQLIS(KIDIA:KFDIA,KTDIA:KLEV),&
        & ZCSGC(KIDIA:KFDIA,KTDIA:KLEV)&
        & *(ZQL(KIDIA:KFDIA,KTDIA:KLEV)+ZQI(KIDIA:KFDIA,KTDIA:KLEV))&
        & +(1._JPRB-ZCSGC(KIDIA:KFDIA,KTDIA:KLEV))*ZQLIS(KIDIA:KFDIA,KTDIA:KLEV))
      ZQC_DET_PCMT(:,:)=0._JPRB
    ELSE
      ! Microphysics occurs in the resolved state.
      ZTMIC(KIDIA:KFDIA,KTDIA:KLEV)=YDMF_PHYS_STATE%T(KIDIA:KFDIA,KTDIA:KLEV)
      ZQMIC(KIDIA:KFDIA,KTDIA:KLEV)=ZQV(KIDIA:KFDIA,KTDIA:KLEV)
    ENDIF

    CALL ACPLUIZ ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,NTPLUI,KLEV,&
     & ZTMIC, ZQMIC, ZQL, ZQI, ZQR, ZQS,&
     & YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_STATE%YCPG_DYN%RCP%R,&
     & ZNEBS, ZQLIS, ZNEB_CVPP, ZQLI_CVPP,ZQC_DET_PCMT,&
     & ZTENHA, ZTENQVA, LADJCLD, YDMF_PHYS_STATE%YCPG_DYN%PHI,&
     & YDMF_PHYS_STATE%YGSP_RR%T, YDMF_PHYS%TMP%APLPAR%FLU%NEIJ, YDMF_PHYS_SURF%GSD_VF%PLSM, PGM, ZVETAF,&
     & YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN,&
     & YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, ZSEDIQL, ZSEDIQI )

    IF(LGPCMT.AND.LCVNHD) THEN
      ! Evaporation processes within convective environment
      ! will be used by the convection, to compute their feedback on convective updrafts. 
      ! This information is put here in PDDAL.
      YDVARS%DAL%T0(:,:)=0._JPRB
      DO JLEV=KTDIA,KLEV
        DO JLON=KIDIA,KFDIA
          YDVARS%DAL%T0(JLON,JLEV)=( &
            & +YDMF_PHYS%OUT%FPEVPSL(JLON,JLEV)-YDMF_PHYS%OUT%FPEVPSL(JLON,JLEV-1)   &
            & +YDMF_PHYS%OUT%FPEVPSN(JLON,JLEV)-YDMF_PHYS%OUT%FPEVPSN(JLON,JLEV-1)   &
            & +YDMF_PHYS%OUT%FPEVPCL(JLON,JLEV)-YDMF_PHYS%OUT%FPEVPCL(JLON,JLEV-1)   &
            & +YDMF_PHYS%OUT%FPEVPCN(JLON,JLEV)-YDMF_PHYS%OUT%FPEVPCN(JLON,JLEV-1)) &
            & /RG*YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP(JLON,JLEV)
        ENDDO
      ENDDO
      IF(LMUSCLFA) CALL WRSCMR(NMUSCLFA,'PFPEVPSL0',YDMF_PHYS%OUT%FPEVPSL,KLON,KLEV+1)
      IF(LMUSCLFA) CALL WRSCMR(NMUSCLFA,'PFPEVPSN0',YDMF_PHYS%OUT%FPEVPSN,KLON,KLEV+1)
      IF(LMUSCLFA) CALL WRSCMR(NMUSCLFA,'PFPEVPCL0',YDMF_PHYS%OUT%FPEVPCL,KLON,KLEV+1)
      IF(LMUSCLFA) CALL WRSCMR(NMUSCLFA,'PFPEVPCN0',YDMF_PHYS%OUT%FPEVPCN,KLON,KLEV+1)
    ENDIF
  ENDIF

  IF ( LNEBCO ) THEN
    CALL ACNEBC ( YDPHY0,KIDIA,KFDIA,KLON,NTNEBU,KLEV,&
     & YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,INLAB,&
     & INND,&
     & YDMF_PHYS_SURF%GSD_VH%PTCCH,YDMF_PHYS_SURF%GSD_VH%PSCCH,YDMF_PHYS_SURF%GSD_VH%PBCCH)
  ENDIF
  IF ( LGWDC ) THEN
    CALL ACDRAC ( YDPHY0,YDPHY2,KIDIA,KFDIA,KLON,NTDRAG,KLEV,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYD, ZNBVNO, YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V,&
     & YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS_SURF%GSD_VH%PSCCH, YDMF_PHYS_SURF%GSD_VH%PBCCH,&
     & YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV)
  ENDIF

  ZFLX_LOTT_GWU(:, :) = 0.0_JPRB
  ZFLX_LOTT_GWV(:, :) = 0.0_JPRB

  IF ( LNORGWD ) THEN

    ! Inversion du sens des niveaux
    DO JLON = KIDIA, KFDIA
      DO JLEV = 1, KLEV
        Z_PP(JLON, JLEV) = YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(JLON, KLEV - JLEV + 1)
        Z_UU(JLON, JLEV) = YDMF_PHYS_STATE%U(JLON, KLEV - JLEV + 1)
        Z_VV(JLON, JLEV) = YDMF_PHYS_STATE%V(JLON, KLEV - JLEV + 1)
        Z_TT(JLON, JLEV) = YDMF_PHYS_STATE%T(JLON, KLEV - JLEV + 1)
        Z_VO(JLON, JLEV) = YDMF_PHYS_STATE%VOR(JLON, KLEV - JLEV + 1)
        ZD_U(JLON, JLEV) = TSPHY * YDMF_PHYS_STATE%P1NOGW(JLON, KLEV - JLEV + 1)
        ZD_V(JLON, JLEV) = TSPHY * YDMF_PHYS_STATE%P2NOGW(JLON, KLEV - JLEV + 1)
      ENDDO
    ENDDO

    ZPRECGWD(KIDIA:KFDIA) = MAX(0.0_JPRB,                  &
      & YDMF_PHYS%OUT%FPLCL(KIDIA:KFDIA,KLEV) + YDMF_PHYS%OUT%FPLCN(KIDIA:KFDIA,KLEV))
    
    CALL ACNORGWD(YDNORGWD, KIDIA, KFDIA, KLON, KLEV, TSPHY, &
                    & Z_PP, PGEMU, Z_TT, Z_UU, Z_VV, Z_VO, ZPRECGWD, &
                    & ZD_U, ZD_V)
    
    DO JLON = KIDIA, KFDIA
      DO JLEV = 1, KLEV
        YDMF_PHYS_STATE%P1NOGW(JLON, JLEV) = ZD_U(JLON, KLEV - JLEV + 1) / TSPHY
        YDMF_PHYS_STATE%P2NOGW(JLON, JLEV) = ZD_V(JLON, KLEV - JLEV + 1) / TSPHY
      ENDDO
    ENDDO

    !-- CALCUL DU FLUX, PAR INTEGRATION DE LA TENDANCE DE HAUT EN BAS.
    !   LES FLUX SONT SUPPOSES NULS AU PREMIER NIVEAU (1 = TOP) DE CALCUL.
    DO JLEV = 1, KLEV
       DO JLON = KIDIA, KFDIA
          ZFLX_LOTT_GWU(JLON, JLEV) = ZFLX_LOTT_GWU(JLON, JLEV - 1) - YDMF_PHYS_STATE%P1NOGW(JLON, JLEV)*YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP(JLON, JLEV)/RG
          ZFLX_LOTT_GWV(JLON, JLEV) = ZFLX_LOTT_GWV(JLON, JLEV - 1) - YDMF_PHYS_STATE%P2NOGW(JLON, JLEV)*YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP(JLON, JLEV)/RG
          YDMF_PHYS%OUT%STRCU(JLON, JLEV) = YDMF_PHYS%OUT%STRCU(JLON, JLEV) + ZFLX_LOTT_GWU(JLON, JLEV)
          YDMF_PHYS%OUT%STRCV(JLON, JLEV) = YDMF_PHYS%OUT%STRCV(JLON, JLEV) + ZFLX_LOTT_GWV(JLON, JLEV)          
    !+ TODOLATER +!      PSTRNORGWDU(JLON, JLEV) = PSTRNORGWDU(JLON, JLEV) + ZFLX_LOTT_GWU(JLON, JLEV)
    !+ TODOLATER +!      PSTRNORGWDV(JLON, JLEV) = PSTRNORGWDV(JLON, JLEV) + ZFLX_LOTT_GWV(JLON, JLEV)          
       ENDDO
    ENDDO

    ! NO VERTICAL DIFFUSION IN THE MIDDLE ATMOSPHERE
    DO JLEV = 1, NORGWD_NNOVERDIF
       DO JLON = KIDIA, KFDIA
          YDMF_PHYS%OUT%STRTU(JLON, JLEV) = 0.0
          YDMF_PHYS%OUT%STRTV(JLON, JLEV) = 0.0          
       ENDDO
    ENDDO

    !+ TODOLATER +! IF (LNOWINDTEND) THEN
    !+ TODOLATER +!  DO JLEV = 0, KLEV
    !+ TODOLATER +!    DO JLON = KIDIA, KFDIA
    !+ TODOLATER +!      PSTRCU(JLON,JLEV) = 0._JPRB
    !+ TODOLATER +!      PSTRCV(JLON,JLEV) = 0._JPRB
    !+ TODOLATER +!    ENDDO
    !+ TODOLATER +!  ENDDO
    !+ TODOLATER +! ENDIF

  ENDIF

!*
!     ------------------------------------------------------------------
!         SAUVEGARDE DES FLUX DE PRECIPITATION CONVECTIVE ET STRATIFORME.

  IF(LNEBN.OR.LNEBR.OR.LRRGUST) THEN
    IF (LFPCOR) THEN
      IF (LCDDPRO) THEN
        DO JLEV=KTDIA-1,KLEV
          DO JLON=KIDIA,KFDIA
            YDMF_PHYS%TMP%APLPAR%PFL%FPLCH(JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)+YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
          ENDDO
        ENDDO
      ELSE
        DO JLEV=KTDIA-1,KLEV
          DO JLON=KIDIA,KFDIA
            YDMF_PHYS%TMP%APLPAR%PFL%FPLCH(JLON,JLEV)=ZFPCOR(JLON,JLEV)
          ENDDO
        ENDDO
      ENDIF
    ELSE
      DO JLEV=KTDIA-1,KLEV
        DO JLON=KIDIA,KFDIA
          YDMF_PHYS%TMP%APLPAR%PFL%FPLCH(JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)+YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  IF(LRRGUST) THEN
    DO JLEV=KTDIA-1,KLEV
      DO JLON=KIDIA,KFDIA
        YDMF_PHYS%TMP%APLPAR%PFL%FPLSH(JLON,JLEV)=YDMF_PHYS%OUT%FPLSL(JLON,JLEV)+YDMF_PHYS%OUT%FPLSN(JLON,JLEV)+YDMF_PHYS%OUT%FPLSG(JLON,JLEV)
      ENDDO
    ENDDO
  ENDIF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! UPDATE TRANSPORT FLUXES DUE TO SEDIMENTATION OF CLOUDS.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  DO JLEV=KTDIA,KLEV
    DO JLON=KIDIA,KFDIA
      YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)+ZSEDIQL(JLON,JLEV)
      YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)+ZSEDIQI(JLON,JLEV)
      ZFPLSL (JLON,JLEV)=YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)+YDMF_PHYS%OUT%FPLSL (JLON,JLEV)
      ZFPLSN (JLON,JLEV)=YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)+YDMF_PHYS%OUT%FPLSN (JLON,JLEV)
    ENDDO
  ENDDO
  IF (LGRAPRO) THEN
    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA
        ZFPLSN (JLON,JLEV)=ZFPLSN(JLON,JLEV)+YDMF_PHYS%OUT%FPLSG(JLON,JLEV)
      ENDDO
    ENDDO
  ENDIF

  ! - - - - - - - - - - - - - - - - -
  ! CORRECT NEGATIVE WATER CONTENTS.
  ! - - - - - - - - - - - - - - - - -

  IF ( LCONDWT.AND.LPROCLD.AND..NOT.LGPCMT ) THEN
    CALL QNGCOR (YDPHY2,KIDIA, KFDIA, KLON, NTPLUI, KLEV,&
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

  ENDIF

!*
!     ------------------------------------------------------------------
!     12. - BILAN HYDRIQUE DU SOL
!     ---------------------------
  IF ( LMSE ) THEN

    DO JLON=KIDIA,KFDIA
      PGPAR(JLON,MRAIN)=ZFPLSL(JLON,KLEV)+YDMF_PHYS%OUT%FPLCL(JLON,KLEV)
      PGPAR(JLON,MSNOW)=ZFPLSN(JLON,KLEV)+YDMF_PHYS%OUT%FPLCN(JLON,KLEV)
    ENDDO

  ELSE

    IF ( LSFHYD.AND.LSOLV ) THEN
      CALL ACDROV ( YDMODEL%YRML_PHY_MF,KIDIA,KFDIA,KLON,KLEV,KCSS,&
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
  IF ( LRRMES ) THEN
    CALL ACDRME ( YDPHY2,YDTOPH,KIDIA,KFDIA,KLON,NTDRME,KLEV,&
     & YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDMF_PHYS_STATE%T,ZQV,YDMF_PHYS_STATE%U,YDMF_PHYS_STATE%V,&
     & YDMF_PHYS%OUT%FRMH,YDMF_PHYS%TMP%APLPAR%MSC%FRMQ,YDMF_PHYS%OUT%STRMU,YDMF_PHYS%OUT%STRMV,&
     & RMESOT,RMESOQ,RMESOU,STTEM)
  ENDIF

!*
!     ------------------------------------------------------------------
!     14.- FLUX PHOTO-CHIMIQUE D'OZONE
!     ------------------------------------------------------------------
  IF ( LOZONE ) THEN
    CALL ACOZONE ( YDPHY2,YDTOPH,KIDIA,KFDIA,KLON,NTOZON,KLEV,KVCLIS,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,PKOZO,ZQO3(1,1),YDMF_PHYS_STATE%T,YDMF_PHYS%TMP%RDG%MU0,&
     & YDMF_PHYS%OUT%FCHOZ)

!           DIFFUSION TURBULENTE/DEPOT SEC DE L'OZONE

    CALL ACDIFOZ ( YDRIP,YDPHY2,YDMODEL%YRML_PHY_MF%YRVDOZ,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
     & YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP,YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FRSO,&
     & ZKTROV,ZQO3(1,1),YDMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,YDMF_PHYS_STATE%T,ZXTROV,&
     & YDMF_PHYS%TMP%APLPAR%FLU%NEIJ,PGEMU,YDMF_PHYS_SURF%GSD_VV%PIVEG,&
     & YDMF_PHYS%OUT%FCHOZ)
  ENDIF
!*
!     ------------------------------------------------------------------
!     15.- CALCUL DES FLUX D'ENTHALPIE ET DE CHALEUR SENSIBLE LIES AUX
!          PRECIPITATIONS EN FONCTION DES FLUX DE PRECIPITATION
!          ET DE CONDENSATION.
!     ------------------------------------------------------------------

  ! STORE THE PSEUDO-HISTORIC SURFACE PRECIPITATION SENSIBLE HEAT FLUX
  ! -------------------------------------------------------------------

  IF (LPHSPSH) THEN
    DO JLON=KIDIA, KFDIA
       YDMF_PHYS_SURF%GSD_VH%PSPSH(JLON)=(YDMF_PHYS_STATE%T(JLON,KLEV)-YDMF_PHYS_STATE%YGSP_RR%T(JLON)) * ((&
          & RCW-RCPD)*(YDMF_PHYS%OUT%FPLSL(JLON,KLEV)+YDMF_PHYS%OUT%FPLCL(JLON,KLEV))&
          & +(RCS-RCPD)*(YDMF_PHYS%OUT%FPLSN(JLON,KLEV)+YDMF_PHYS%OUT%FPLCN(JLON,KLEV)) )
    ENDDO
    IF (LGRAPRO) THEN
      YDMF_PHYS_SURF%GSD_VH%PSPSH(JLON)=YDMF_PHYS_SURF%GSD_VH%PSPSH(JLON)+(YDMF_PHYS_STATE%T(JLON,KLEV)-YDMF_PHYS_STATE%YGSP_RR%T(JLON)) * (&
        &(RCS-RCPD)*(YDMF_PHYS%OUT%FPLSG(JLON,KLEV)+YDMF_PHYS%OUT%FPLCG(JLON,KLEV)))
    ENDIF
  ENDIF !LPHSPSH

! Store radiative cloudiness in GFL structure for ISP, Historical files or PostProcessing
  IF (YIRAD%LGP) YDVARS%IRAD%T1(KIDIA:KFDIA,:) = YDCPG_MISC%QICE(KIDIA:KFDIA,:)
  IF (YLRAD%LGP) YDVARS%LRAD%T1(KIDIA:KFDIA,:) = YDCPG_MISC%QLI (KIDIA:KFDIA,:)
  IF (YA%LGP)    YDVARS%A%T1(KIDIA:KFDIA,:) = YDCPG_MISC%NEB (KIDIA:KFDIA,:)


!*
!     -----------------------------------------------------------------------
!     16.- DDH FLEXIBLES POUR LES CHAMPS DE SURFACE VARIABLES/FLUX/TENDANCES.
!     -----------------------------------------------------------------------

  IF(LFLEXDIA) THEN
    IF (LDDH_OMP) THEN
      ! Store in DDH the clouds used for radiative calculations
      ZTMPAF(KIDIA:KFDIA,1:KLEV)=YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP(KIDIA:KFDIA,1:KLEV)*YDCPG_MISC%QICE(KIDIA:KFDIA,1:KLEV)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAF,'VIR',YDDDH)
      ZTMPAF(KIDIA:KFDIA,1:KLEV)=YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP(KIDIA:KFDIA,1:KLEV)*YDCPG_MISC%QLI(KIDIA:KFDIA,1:KLEV)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAF,'VLR',YDDDH)
      ZTMPAF(KIDIA:KFDIA,1:KLEV)=YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP(KIDIA:KFDIA,1:KLEV)*YDCPG_MISC%NEB(KIDIA:KFDIA,1:KLEV)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAF,'VNT',YDDDH)

      ! surface variables
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS_SURF%GSD_VF%PLSM,'VLSM',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS_STATE%YGSP_RR%T(KIDIA:KFDIA)/YDMF_PHYS%OUT%CT(KIDIA:KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'VSTS',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*RTINER*YDMF_PHYS_STATE%YGSP_SB%T(KIDIA:KFDIA,KCSS)/YDMF_PHYS%OUT%CT(KIDIA:KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'VSTP',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS_STATE%YGSP_RR%W(KIDIA:KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'VSWS',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS_STATE%YGSP_SB%Q(KIDIA:KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'VSWP',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS_STATE%YGSP_SG%F(KIDIA:KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'VSWN',YDDDH)

      ! surface free-style variables
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%TCLS,'VTCLS',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%QCLS,'VQCLS',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%UCLS,'VUCLS',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%VCLS,'VVCLS',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%NUCLS,'VNUCLS',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%NVCLS,'VNVCLS',YDDDH,CDTYPE='S')
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_STATE%YCPG_DYN%PHI(KIDIA:KFDIA,KLEV)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'VSPHI',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%GZ0,'VSGZ0',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%GZ0H,'VSGZH',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%ALB,'VSALB',YDDDH,CDTYPE='S')
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%CLPH,'VSCPH',YDDDH,CDTYPE='S')

      ! surface fluxes free stlye or not
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FRSO(KIDIA:KFDIA,KLEV,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSRAYSODW',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FRTH(KIDIA:KFDIA,KLEV,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSRAYTHUP',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FCLL(KIDIA:KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSCHLATLI',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FCLN(KIDIA:KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSCHLATNE',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FCS(KIDIA:KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSCHLATSF',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FCHSP(KIDIA:KFDIA,KCSS)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSCHGSNPL',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FONTE(KIDIA:KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FFONTESLI',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FPLSL(KIDIA:KFDIA,KLEV)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSPRELIGE',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FPLSN(KIDIA:KFDIA,KLEV)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSPRENEGE',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FPLCL(KIDIA:KFDIA,KLEV)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSPRELICO',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FPLCN(KIDIA:KFDIA,KLEV)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSPRENECO',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FEVL(KIDIA:KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSEVAPLIQ',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FEVN(KIDIA:KFDIA,1)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSEVAPNEG',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FLWSP(KIDIA:KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSLIQSNPL',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%RUISS(KIDIA:KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FRUISSURF',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%RUISP(KIDIA:KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FRUISSPRF',YDDDH)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,YDMF_PHYS%OUT%FRTHDS,'FSRAYTHDS',YDDDH)
      ZTMPAS(KIDIA:KFDIA)=-YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*RLMLT*YDMF_PHYS%OUT%FONTE(KIDIA:KFDIA)
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FFONTESNE',YDDDH)

      ZEPSA = 1.E-4_JPRB
      DO JLON=KIDIA, KFDIA
        ZTMPAS(JLON)=YDMF_PHYS%OUT%FRSO(JLON,KLEV,1)/MAX(1.0_JPRB-YDMF_PHYS%OUT%ALB(JLON),ZEPSA)
      ENDDO
      CALL NEW_ADD_FIELD_2D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPAS,'FSRAYSOLR',YDDDH)
      IF (NPRINTLEV >= 2) THEN
      WRITE(NULOUT,*) 'LFLEXDIA ARPEGE WITH NTOTSURF = ',NTOTSURF,&
         & ' AND NTOTSVAR = ',NTOTSVAR, ' AND NTOTSVFS = ',NTOTSVFS
      ENDIF

    ELSE  

      ! Store in DDH the clouds used for radiative calculations
      ZTMPAF(KIDIA:KFDIA,1:KLEV)=YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP(KIDIA:KFDIA,1:KLEV)*YDCPG_MISC%QICE(KIDIA:KFDIA,1:KLEV)
      CALL ADD_FIELD_3D(YDLDDH,ZTMPAF,'VIR','V','ARP',.TRUE.,.TRUE.)
      ZTMPAF(KIDIA:KFDIA,1:KLEV)=YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP(KIDIA:KFDIA,1:KLEV)*YDCPG_MISC%QLI(KIDIA:KFDIA,1:KLEV)
      CALL ADD_FIELD_3D(YDLDDH,ZTMPAF,'VLR','V','ARP',.TRUE.,.TRUE.)
      ZTMPAF(KIDIA:KFDIA,1:KLEV)=YDMF_PHYS_STATE%YCPG_PHY%XYB%DELP(KIDIA:KFDIA,1:KLEV)*YDCPG_MISC%NEB(KIDIA:KFDIA,1:KLEV)
      CALL ADD_FIELD_3D(YDLDDH,ZTMPAF,'VNT','V','ARP',.TRUE.,.TRUE.)

      ! surface variables
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS_SURF%GSD_VF%PLSM,'VLSM','V','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS_STATE%YGSP_RR%T(KIDIA:KFDIA)/YDMF_PHYS%OUT%CT(KIDIA:KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'VSTS','V','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*RTINER*YDMF_PHYS_STATE%YGSP_SB%T(KIDIA:KFDIA,KCSS)/YDMF_PHYS%OUT%CT(KIDIA:KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'VSTP','V','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS_STATE%YGSP_RR%W(KIDIA:KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'VSWS','V','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS_STATE%YGSP_SB%Q(KIDIA:KFDIA,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'VSWP','V','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS_STATE%YGSP_SG%F(KIDIA:KFDIA,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'VSWN','V','ARP',.TRUE.,.TRUE.)

      ! surface free-style variables
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%TCLS,'VTCLS','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%QCLS,'VQCLS','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%UCLS,'VUCLS','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%VCLS,'VVCLS','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%NUCLS,'VNUCLS','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%NVCLS,'VNVCLS','S','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_STATE%YCPG_DYN%PHI(KIDIA:KFDIA,KLEV)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'VSPHI','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%GZ0,'VSGZ0','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%GZ0H,'VSGZH','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%ALB,'VSALB','S','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%CLPH,'VSCPH','S','ARP',.TRUE.,.TRUE.)

      ! surface fluxes free stlye or not
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FRSO(KIDIA:KFDIA,KLEV,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSRAYSODW','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FRTH(KIDIA:KFDIA,KLEV,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSRAYTHUP','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FCLL(KIDIA:KFDIA,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSCHLATLI','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FCLN(KIDIA:KFDIA,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSCHLATNE','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FCS(KIDIA:KFDIA,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSCHLATSF','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FCHSP(KIDIA:KFDIA,KCSS)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSCHGSNPL','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FONTE(KIDIA:KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FFONTESLI','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FPLSL(KIDIA:KFDIA,KLEV)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSPRELIGE','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FPLSN(KIDIA:KFDIA,KLEV)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSPRENEGE','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FPLCL(KIDIA:KFDIA,KLEV)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSPRELICO','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FPLCN(KIDIA:KFDIA,KLEV)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSPRENECO','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FEVL(KIDIA:KFDIA,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSEVAPLIQ','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FEVN(KIDIA:KFDIA,1)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSEVAPNEG','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%FLWSP(KIDIA:KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSLIQSNPL','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%RUISS(KIDIA:KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FRUISSURF','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*YDMF_PHYS%OUT%RUISP(KIDIA:KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FRUISSPRF','F','ARP',.TRUE.,.TRUE.)
      CALL ADD_FIELD_2D(YDLDDH,YDMF_PHYS%OUT%FRTHDS,'FSRAYTHDS','F','ARP',.TRUE.,.TRUE.)
      ZTMPAS(KIDIA:KFDIA)=-YDMF_PHYS_SURF%GSD_VF%PLSM(KIDIA:KFDIA)*RLMLT*YDMF_PHYS%OUT%FONTE(KIDIA:KFDIA)
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FFONTESNE','F','ARP',.TRUE.,.TRUE.)
 
      ZEPSA = 1.E-4_JPRB
      DO JLON=KIDIA, KFDIA
        ZTMPAS(JLON)=YDMF_PHYS%OUT%FRSO(JLON,KLEV,1)/MAX(1.0_JPRB-YDMF_PHYS%OUT%ALB(JLON),ZEPSA)
      ENDDO
      CALL ADD_FIELD_2D(YDLDDH,ZTMPAS,'FSRAYSOLR','F','ARP',.TRUE.,.TRUE.)
      IF (NPRINTLEV >= 2) THEN
      WRITE(NULOUT,*) 'LFLEXDIA ARPEGE WITH NTOTSURF = ',NTOTSURF,&
         & ' AND NTOTSVAR = ',NTOTSVAR, ' AND NTOTSVFS = ',NTOTSVFS
      ENDIF
    ENDIF
  ENDIF

  IF (LECT.AND.LRAFTKE) THEN
     ! DCAPE due to precipitation evaporation.
     CALL ACEVADCAPE(KIDIA,KFDIA, KLON, KLEV,&
      & YDMF_PHYS%OUT%FPLSL,YDMF_PHYS%OUT%FPLSN,YDMF_PHYS%OUT%FPLCL,YDMF_PHYS%OUT%FPLCN,YDMF_PHYS_STATE%YCPG_PHY%PREHYD,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%YCPG_DYN%RCP%CP,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_DYN%PHI,ZDCAPE)
     ! Gusts.
     ZQGM=ZEPSNEB
     CALL ACCLDIA(YDXFU,YDPHY,YDPHY2,YDTOPH, KIDIA,KFDIA, KLON, KLEV,&
      & YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS_STATE%U, YDMF_PHYS_STATE%V, YDMF_PHYS%OUT%CAPE, ZDCAPE, ZTKE1, YDMF_PHYS_STATE%YCPG_DYN%PHIF, POROG, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST, ZBLH, KCLPH)
     YDMF_PHYS%OUT%CLPH(KIDIA:KFDIA)=MIN(XMAXLM,MAX(XMINLM,ZBLH(KIDIA:KFDIA)))
  ENDIF

  IF (LXVISI.OR.LXVISI2) THEN
     CALL ACVISIH(YDMODEL%YRML_PHY_MF%YRPHY%YRDVISI,KIDIA,KFDIA,KLON,KTDIA,KLEV,YDMF_PHYS_STATE%YCPG_DYN%PHI,YDMF_PHYS_STATE%YCPG_DYN%PHIF,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF,&
      & YDMF_PHYS_STATE%T,YDMF_PHYS_STATE%YCPG_DYN%RCP%R,YDCPG_MISC%QLI,YDCPG_MISC%QICE,ZQR,ZQS,ZQGM,YDMF_PHYS%OUT%VISICLD, YDMF_PHYS%OUT%VISIHYD,YDMF_PHYS%OUT%MXCLWC)
  ENDIF

ENDIF !LMPHYS
!----------------------------------------------------
!  CALCUL DE DEPOT HUMIDE POUR LES AEROSOLS DESERTIQUES
!----------------------------------------------------
IF (LMDUST.AND.(NGFL_EXT/=0).AND.LRDEPOS) THEN

IKRR=6
ISPLITR=1
ZEVAP=0.0_JPRB
ZZDEP=0.0_JPRB

  DO JLEV = 1 , KLEV
    DO JLON = KIDIA,KFDIA
      ZZDEP(JLON,1,JLEV)=ZZI_APHI(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO


  DO JLEV=1,KLEV
    DO JLON= KIDIA, KFDIA
      ZQDM(JLON,JLEV)=1._JPRB-ZQV(JLON,JLEV)-ZQL(JLON,JLEV)-ZQR(JLON,JLEV) &
       & -ZQI(JLON,JLEV)-ZQS(JLON,JLEV)
    ENDDO
  ENDDO

  DO JLEV = 1, KLEV
    DO JLON = KIDIA,KFDIA
      ZZI_RM(JLON,1,JLEV,2)=ZQL(JLON,JLEV)&
       & /ZQDM(JLON,JLEV)
    ENDDO
  ENDDO

  DO JLEV = 1, KLEV
    DO JLON = KIDIA,KFDIA
      ZZI_RM(JLON,1,JLEV,3)=ZQR(JLON,JLEV)&
       & /ZQDM(JLON,JLEV)
    ENDDO
  ENDDO

  DO JLEV = 1, KLEV
    DO JLON = KIDIA,KFDIA
      ZEVAP(JLON,1,JLEV)=YDMF_PHYS%OUT%FPEVPSL(JLON,JLEV)&
       & +YDMF_PHYS%OUT%FPEVPCL(JLON,JLEV)
    ENDDO
  ENDDO

  DO JGFL=1,NGFL_EXT
    DO JLON=1,KLON
      DO JLEV=1,KLEV
        ZZI_SVM(JLON,1,JLEV,JGFL)=MAX(0._JPRB,ZSVM(JLON,JLEV,JGFL))
      ENDDO
    ENDDO
  ENDDO

    CALL ARO_WETDEP(ILONMNH,KLEV,NGFL_EXT,IKRR, PDT,&
                 & ZZI_SVM(KIDIA:KFDIA,:,:,1:NGFL_EXT),&
                 & ZZDEP(KIDIA:KFDIA,:,:),&
                 & ZZI_PABSM(KIDIA:KFDIA,:,:),&
                 & ZZI_THM(KIDIA:KFDIA,:,:),&
                 & ZZI_RHODREFM(KIDIA:KFDIA,:,:),&
                 & KSTEP+1,&
                 & ZZI_RM(KIDIA:KFDIA,:,:,:),&
                 & ZEVAP(KIDIA:KFDIA,:,:),&
                 & ISPLITR               )
! return to tendency
   DO JGFL=1,NGFL_EXT
     DO JLON=1,KLON
       DO JLEV=1,KLEV
         ZPSV(JLON,JLEV,JGFL)=ZZI_SVM(JLON,1,JLEV,JGFL)
       ENDDO
     ENDDO
   ENDDO

      DO JGFL=1,NGFL_EXT
        DO JLON=1,KLON
          DO JLEV=1,KLEV
            PTENDEXT_DEP(JLON,JLEV,JGFL)=(ZPSV(JLON,JLEV,JGFL)-ZSVM(JLON,JLEV,JGFL))*ZINVDT
          ENDDO
        ENDDO
      ENDDO

ENDIF ! ENDIF  (LMDUST & NGFL_EXT & LRDEPOS)

IF(LMUSCLFA) THEN
  DO JLEV=0,KLEV
    DO JLON=1,KLON
      ZFPLS(JLON,JLEV)=YDMF_PHYS%OUT%FPLSL(JLON,JLEV)+YDMF_PHYS%OUT%FPLSN(JLON,JLEV)
      ZFPLC(JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)+YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
      ZFPL (JLON,JLEV)=ZFPLC (JLON,JLEV)+ZFPLS (JLON,JLEV)
    ENDDO
  ENDDO
  CALL WRSCMR(NMUSCLFA,'ZFPLS',ZFPLS,KLON,KLEV+1)
  CALL WRSCMR(NMUSCLFA,'ZFPLC',ZFPLC,KLON,KLEV+1)
  CALL WRSCMR(NMUSCLFA,'ZFPL',ZFPL,KLON,KLEV+1)
ENDIF

!     ------------------------------------------------------------------

! BAYRAD
! Compute convective hydrometeors mixing ratio from diagnistic fluxes
!--------------------------------------------------------------------
IF ( (.NOT. LGPCMT) ) THEN

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

   YDMF_PHYS%TMP%APLPAR%BAY%QRCONV(KIDIA:KFDIA,:) =  0.0_JPRB
   YDMF_PHYS%TMP%APLPAR%BAY%QSCONV(KIDIA:KFDIA,:) =  0.0_JPRB


   DO JLEV=1,KLEV
     DO JLON= KIDIA, KFDIA
       YDMF_PHYS%TMP%APLPAR%BAY%QRCONV(JLON,JLEV)  = 0.001_JPRB * ( (ABS(YDMF_PHYS%OUT%FPLCL(JLON,JLEV)) * ZCOEFRAIN(1)) ** ZCOEFRAIN(2) )
       YDMF_PHYS%TMP%APLPAR%BAY%QSCONV(JLON,JLEV)  = 0.001_JPRB * ( (ABS(YDMF_PHYS%OUT%FPLCN(JLON,JLEV)) * ZCOEFSNOW(1)) ** ZCOEFSNOW(2) )    
     ENDDO
   ENDDO


   ! Convert density [kg/m3] to mixing ratio [kg/kg]
   ! R_dry (dry air constant)

   ZDE2MR(KIDIA:KFDIA,:)  = RD * YDMF_PHYS_STATE%T(KIDIA:KFDIA,:) / YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(KIDIA:KFDIA,:)
   YDMF_PHYS%TMP%APLPAR%BAY%QRCONV(KIDIA:KFDIA,:) = YDMF_PHYS%TMP%APLPAR%BAY%QRCONV(KIDIA:KFDIA,:) * ZDE2MR(KIDIA:KFDIA,:)
   YDMF_PHYS%TMP%APLPAR%BAY%QSCONV(KIDIA:KFDIA,:) = YDMF_PHYS%TMP%APLPAR%BAY%QSCONV(KIDIA:KFDIA,:) * ZDE2MR(KIDIA:KFDIA,:)

ENDIF


! Precipitation Type

! Compute wet-bulb temperature at 2 meters (suppose homogeneity of qv/ql/qi )
!ZPCLS(KIDIA:KFDIA)=PAPRS(KIDIA:KFDIA,KLEV)-2._JPRB/ZZZF(KIDIA:KFDIA,1,KLEV)*&
!                 &(PAPRS(KIDIA:KFDIA,KLEV)-PAPRSF(KIDIA:KFDIA,KLEV))
CALL PPWETPOINT(YDPHY,KIDIA,KFDIA,KLON,ZPCLS,YDMF_PHYS%OUT%TCLS,&
  & YDMF_PHYS_STATE%Q(:,KLEV),YDMF_PHYS_STATE%L(:,KLEV),YDMF_PHYS_STATE%I(:,KLEV),YDMF_PHYS%OUT%TPWCLS)


IF (LDPRECIPS .OR. LDPRECIPS2) THEN
   ! Defined precipitation type
   !
   ! Compute wet-bulb temperature
  DO JLEV=1,KLEV
    CALL PPWETPOINT(YDPHY,KIDIA,KFDIA,KLON,YDMF_PHYS_STATE%YCPG_PHY%PREHYDF(:,JLEV),YDMF_PHYS_STATE%T(:,JLEV),&
     & YDMF_PHYS_STATE%Q(:,JLEV),YDMF_PHYS_STATE%L(:,JLEV),YDMF_PHYS_STATE%I(:,JLEV),ZTW(:,JLEV))
  ENDDO

  DO JLON=1,KLON
      ZFPLS(JLON,KLEV)=YDMF_PHYS%OUT%FPLCN(JLON,KLEV)+YDMF_PHYS%OUT%FPLSN(JLON,KLEV)
      ZFPLC(JLON,KLEV)=YDMF_PHYS%OUT%FPLCL(JLON,KLEV)+YDMF_PHYS%OUT%FPLSL(JLON,KLEV)
      ZFPLSG(JLON,KLEV)=0._JPRB
  ENDDO

  !initialisation de ZZZ
  DO JLEV = 1,KLEV
    DO JLON = KIDIA,KFDIA
      ZZZ(JLON,1,JLEV)=YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO

  !initialisation de ZDZZ
  DO JLEV = 2, KLEV
    DO JLON = KIDIA,KFDIA
      ZDZZ(JLON,1,JLEV)=ZZZ(JLON,1,JLEV-1)-ZZZ(JLON,1,JLEV)
    ENDDO
  ENDDO
  DO JLON = KIDIA,KFDIA
      ZDZZ(JLON,1,1)=YDMF_PHYS_STATE%YCPG_DYN%PHI(JLON,0)*ZINVG-ZZZ(JLON,1,1)
  ENDDO

  IF (LDPRECIPS) THEN

   NDTPRECCUR=INT(MOD(ZSTATI/TSTEP,REAL(NDTPREC)))+1_JPIM

   CALL DPRECIPS (YDPHY%YRDPRECIPS,KIDIA,KFDIA,KLON,KLEV,POROG,YDMF_PHYS%OUT%TPWCLS,YDMF_PHYS%OUT%DIAGH,YDMF_PHYS_STATE%YCPG_DYN%PHIF,&
      & ZDZZ,ZTW,YDMF_PHYS_STATE%L,ZFPLC(:,KLEV),ZFPLS(:,KLEV),ZFPLSG(:,KLEV),YDMF_PHYS%TMP%PRC%DPRECIPS(:,NDTPRECCUR))

  ! WRITE(20,*)'sous aplpar ZDPRECIPS ', ZDPRECIPS(KIDIA:KFDIA,NDTPRECCUR),NDTPRECCUR

  ENDIF

  IF (LDPRECIPS2) THEN

 !Idem for an other time step and an other period
   NDTPRECCUR2=INT(MOD(ZSTATI/TSTEP,REAL(NDTPREC2)))+1_JPIM

   CALL DPRECIPS(YDPHY%YRDPRECIPS,KIDIA,KFDIA,KLON,KLEV,POROG,YDMF_PHYS%OUT%TPWCLS,YDMF_PHYS%OUT%DIAGH,YDMF_PHYS_STATE%YCPG_DYN%PHIF,&
      & ZDZZ,ZTW,YDMF_PHYS_STATE%L,ZFPLC(:,KLEV),ZFPLS(:,KLEV),ZFPLSG(:,KLEV),YDMF_PHYS%TMP%PRC%DPRECIPS2(:,NDTPRECCUR2))

  ! WRITE(20,*)'sous aplpar ZDPRECIPS2',ZDPRECIPS2(KIDIA:KFDIA,NDTPRECCUR2),NDTPRECCUR2

  ENDIF
ENDIF  

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('APLPAR',1,ZHOOK_HANDLE)
END SUBROUTINE APLPAR
