SUBROUTINE APL_ARPEGE(YDMF_PHYS_BASE_STATE, YDMF_PHYS_NEXT_STATE, YDGEOMETRY, YDCPG_DIM, YDCPG_MISC, &
& YDCPG_GPAR, YDCPG_PHY0, YDMF_PHYS, YDCPG_DYN0, YDMF_PHYS_SURF, YDCPG_SL2, YDVARS, YDSURF, YDCFU,   &
& YDXFU, YDMODEL, LDCONFX, PDTPHY, YDDDH)

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

! - 2D (1:KLEV) .

! PKOZO      : CHAMPS POUR LA PHOTOCHIMIE DE L'OZONE (KVCLIS CHAMPS).
! PKOZO      : FIELDS FOR PHOTOCHEMISTERY OF OZONE   (KVCLIS FIELDS).

! PGPAR       : BUFFER FOR 2D FIELDS - CONTAINS PRECIP, ALBEDO, EMISS, TS
!             : SURFACE FLUXES
! - INPUT/OUTPUT 1D
! YDDDH      : DDH superstructure

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
USE CPG_TYPE_MOD       , ONLY : CPG_MISC_TYPE, CPG_DYN_TYPE, &
                              & CPG_SL2_TYPE, CPG_GPAR_TYPE
USE CPG_DIM_TYPE_MOD   , ONLY : CPG_DIM_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD,ONLY : MF_PHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMXFU             , ONLY : TXFU
USE TYPE_MODEL         , ONLY : MODEL
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK    ,DR_HOOK
USE YOMVERT            , ONLY : VP00
USE YOMCST             , ONLY : RG       ,RSIGMA   ,RV       ,RD       ,&
                              & RCPV     ,RETV     ,RCW      ,RCS      ,RLVTT ,&
                              & RLSTT    ,RTT      ,RALPW    ,RBETW    ,RGAMW ,&
                              & RALPS    ,RBETS    ,RGAMS    ,RALPD    ,RBETD ,&
                              & RGAMD    ,RCPD     ,RATM     ,RKAPPA   
USE YOMCT0             , ONLY : LSFORCS, LELAM
USE YOMDYNA            , ONLY : L3DTURB
USE DDH_MIX            , ONLY : TYP_DDH
USE YOMCFU             , ONLY : TCFU !!! for parameters of FLASH
USE SPP_MOD  , ONLY : YSPP, YSPP_CONFIG
USE MF_PHYS_BASE_STATE_TYPE_MOD &
                       , ONLY : MF_PHYS_BASE_STATE_TYPE
USE MF_PHYS_NEXT_STATE_TYPE_MOD &
                       , ONLY : MF_PHYS_NEXT_STATE_TYPE


USE CPG_TYPE_MOD       , ONLY : CPG_PHY_TYPE
USE SC2PRG_MOD         , ONLY : SC2PRG

USE YOMCT0             , ONLY : LTWOTL, LAROME, LCORWAT
USE YOMNUD             , ONLY : NFNUDG   ,LNUDG 
USE YOMSNU             , ONLY : XPNUDG
USE YOMCHET            , ONLY : GCHETN
!     -------------------------------------------------------------------------

IMPLICIT NONE


TYPE (MF_PHYS_BASE_STATE_TYPE), INTENT(IN)    :: YDMF_PHYS_BASE_STATE
TYPE (MF_PHYS_NEXT_STATE_TYPE), INTENT(INOUT) :: YDMF_PHYS_NEXT_STATE
TYPE(GEOMETRY),                 INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_DIM_TYPE),             INTENT(IN)    :: YDCPG_DIM
TYPE(CPG_MISC_TYPE),            INTENT(INOUT) :: YDCPG_MISC
TYPE(CPG_GPAR_TYPE),            INTENT(INOUT) :: YDCPG_GPAR
TYPE(CPG_PHY_TYPE),             INTENT(IN)    :: YDCPG_PHY0
TYPE(MF_PHYS_TYPE),             INTENT(INOUT) :: YDMF_PHYS
TYPE(CPG_DYN_TYPE),             INTENT(IN)    :: YDCPG_DYN0
TYPE(MF_PHYS_SURF_TYPE),        INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(CPG_SL2_TYPE),             INTENT(INOUT) :: YDCPG_SL2
TYPE(FIELD_VARIABLES),          INTENT(INOUT) :: YDVARS
TYPE(TSURF),                    INTENT(IN)    :: YDSURF
TYPE(TCFU),                     INTENT(IN)    :: YDCFU
TYPE(TXFU),                     INTENT(IN)    :: YDXFU
TYPE(MODEL),                    INTENT(IN)    :: YDMODEL
LOGICAL,                        INTENT(IN)    :: LDCONFX
REAL(KIND=JPRB),                INTENT(IN)    :: PDTPHY 
TYPE(TYP_DDH),                  INTENT(INOUT) :: YDDDH


!     ------------------------------------------------------------------
LOGICAL :: LL_SAVE_PHSURF
INTEGER(KIND=JPIM) :: INSTEP_DEB,INSTEP_FIN
INTEGER(KIND=JPIM) :: JROF, JSPP

!     --- UPPER AIR PHYSICAL TENDENCIES.
REAL(KIND=JPRB) :: ZTENDH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)     ! Enthalpy tendency.
     ! Moisture tendency.
  ! Pseudo progn. TKE

! GFL tendencies for APL_AROME (assumes YDMODEL%YRML_GCONF%YGFL%NUMFLDS>=YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
! for now, use Jovi's trick :

!     --- UPPER AIR PHYSICAL TENDENCIES FOR AROME.
!       (the previous one are not used in AROME)
REAL(KIND=JPRB) :: ZTENDT (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)        ! temperature tendency
REAL(KIND=JPRB) :: ZTENDD (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)        ! d  tendency
REAL(KIND=JPRB) :: ZTENDEXT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)      ! GFL EXTRA tendency
  ! GFL EXTRA tendency
REAL(KIND=JPRB) :: ZDIFEXT(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)     ! Extra-GFL fluxes.

REAL(KIND=JPRB) :: ZTENDU (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! U tendency without deep convection contribution
REAL(KIND=JPRB) :: ZTENDV (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! V tendency without deep convection contribution

! SPP
REAL(KIND=JPRB) :: ZGP2DSPP(YDCPG_DIM%KLON,YSPP%N2D)

REAL(KIND=JPRB) :: ZTENDEFB1  (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENDEFB2  (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENDEFB3  (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENDG     (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENDICONV (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENDI     (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENDLCONV (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENDQ     (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENDRCONV (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENDR     (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENDSCONV (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENDS     (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENDTKE   (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENDL     (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

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
REAL(KIND=JPRB) :: ZICEFR1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)! Resolved Condensate ice fraction
REAL(KIND=JPRB) :: ZRHCRI(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! Smith scheme critical RH
REAL(KIND=JPRB) :: ZRHDFDA(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)! RK scheme change in RH over cloud
REAL(KIND=JPRB) :: ZLHS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)   ! Sublimation latent heat
REAL(KIND=JPRB) :: ZLHV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)   ! Evaporation latent heat
REAL(KIND=JPRB) :: ZQSATS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! QSAT of resolved cond./evap. scheme
REAL(KIND=JPRB) :: ZRH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! Temporar storage for updated PRH
REAL(KIND=JPRB) :: ZTW(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! Temporar storage for updated PTW)
REAL(KIND=JPRB) :: ZPOID(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)  ! DP/(RG*DT) FOR A GIVEN LEVEL AND A GIVEN TIME STEP.
REAL(KIND=JPRB) :: ZIPOI(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)  ! INVERSE OF ZPOID.
REAL(KIND=JPRB) :: ZQV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! corrected (for negative values) vapour
REAL(KIND=JPRB) :: ZQI(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! corrected (for negative values) cloud ice
REAL(KIND=JPRB) :: ZQL(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! corrected (for negative values) cloud liquid
REAL(KIND=JPRB) :: ZQR(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! corrected (for negative values) rain
REAL(KIND=JPRB) :: ZQS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! corrected (for negative values) snow
REAL(KIND=JPRB) :: ZTENHA(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENQVA(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZCP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! new cp for turbulent diffusion
REAL(KIND=JPRB) :: ZFCQVNG(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! correction flux increment for neg vapour
REAL(KIND=JPRB) :: ZFCQING(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! correction flux increment for neg ice
REAL(KIND=JPRB) :: ZFCQLNG(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! correction flux increment for neg liquid water
REAL(KIND=JPRB) :: ZFPLSL(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! total liquid water flux: diff+sedi+rain
REAL(KIND=JPRB) :: ZFPLSN(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! total solid water flux: diff+sedi+snow
REAL(KIND=JPRB) :: ZFCQL(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)   ! condensation flux(liquid)
REAL(KIND=JPRB) :: ZFCQI(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)   ! condensation flux(ice)
REAL(KIND=JPRB) :: ZSEDIQL(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! sedimentation flux of cloud liquid water
REAL(KIND=JPRB) :: ZSEDIQI(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! sedimentation flux of cloud ice water
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
REAL(KIND=JPRB) :: ZCONDCVPPL(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de condensation liquide du a CVVPP (KFB)
REAL(KIND=JPRB) :: ZCONDCVPPI(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de condensation glace du a CVVPP (KFB)
REAL(KIND=JPRB) :: ZPRODTH_CVPP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG) ! Flux de production thermique de TKE du a CVPP(KFB)
REAL(KIND=JPRB) :: ZDQV, ZDQI, ZDQL, ZDQR, ZDQS, ZDQC, ZGDT, ZGDTI, ZQV0

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
REAL(KIND=JPRB) :: ZALBD(YDCPG_DIM%KLON,YDCPG_DIM%KSW), ZALBP(YDCPG_DIM%KLON,YDCPG_DIM%KSW),&
                 & ZALB(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZSFSWDIR (YDCPG_DIM%KLON,YDCPG_DIM%KSW), ZSFSWDIF (YDCPG_DIM%KLON,YDCPG_DIM%KSW)
REAL(KIND=JPRB) :: ZTRSODIR (YDCPG_DIM%KLON,YDCPG_DIM%KSW), ZTRSODIF (YDCPG_DIM%KLON,YDCPG_DIM%KSW)
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

REAL(KIND=JPRB) :: ZGEOSLC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTHETAV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

!            LOCAL ARRAYS FOR TKE
! ZCOEFN : COEFFICIENT STATISTIQUE POUR LES FLUX D'EAUX CONDENSEES.

REAL(KIND=JPRB) :: ZCOEFN(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

!            LOCAL ARRAYS FOR ACVPPKF
REAL(KIND=JPRB) :: ZQLI_CVPP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZNEB_CVPP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
!            LOCAL ARRAYS FOR EDKF


!           2-D ARRAY FOR SIMPL.RADIATION SCHEME



INTEGER(KIND=JPIM) :: JCHA, JLEV, JLON, JSG
 ! useful size of klon arrays for mesonh physics
LOGICAL :: LLCLS, LLHMT, LLREDPR

REAL(KIND=JPRB) :: ZAEN, ZAEO, ZEPS, ZEPS0, ZEPSNEB, ZEPSO3


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
REAL(KIND=JPRB) :: ZCFBTH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZDTMSE,ZRHGMT,ZSTATI
REAL(KIND=JPRB) :: &
 & ZCFATH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZCFAU(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),&
 & ZCFBU(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZCFBV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

REAL(KIND=JPRB) :: ZSRAIN(YDCPG_DIM%KLON), ZSSNOW(YDCPG_DIM%KLON), ZSGROUPEL(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZTSN(YDCPG_DIM%KLON)
! FOR Hv
! FOR DUST

 ! SCOND MEMBRE POUR LES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZINVG

! TRAITEMENT DES SCALAIRES PASSIFS

REAL(KIND=JPRB) :: ZDZZ(YDCPG_DIM%KLON,1,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZZZ (YDCPG_DIM%KLON,1,YDCPG_DIM%KFLEVG)

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
LOGICAL            :: LLAERO
REAL(KIND=JPRB)    :: ZAIPCMT(YDCPG_DIM%KLON) ! Activity Index of PCMT: 1. if PCMT is active, 0. else case.


REAL(KIND=JPRB)    :: ZQIC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZQLC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZQRC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZQSC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)    :: ZQLI_CVP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)    :: ZQC_DET_PCMT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB)    :: ZFPLS(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG),ZFPLC(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)

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

! Precipitation type diagnostics
!--------------------------------
REAL(KIND=JPRB)   :: ZFPLSG(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)
INTEGER (KIND=JPIM)  :: IMOC_CLPH (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZBAY_QRCONV (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZBAY_QSCONV (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZDSA_C1 (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZDSA_C2 (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZDSA_CPS (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZDSA_LHS (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZDSA_RS (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZFLU_CDN (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZFLU_CD (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZFLU_CH (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZFLU_EMIS (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZFLU_FEVI (YDCPG_DIM%KLON, 1:YDCPG_DIM%KTSSG+1)
REAL (KIND=JPRB)     :: ZFLU_NEIJ (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZFLU_QSATS (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZFLU_QSAT (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZFLU_VEG (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZMSC_FHP (YDCPG_DIM%KLON, 0:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZMSC_FRMQ (YDCPG_DIM%KLON, 0:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZMSC_LH (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZMSC_LSCPE (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZMSC_QW (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZMSC_TW (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FEFB1 (YDCPG_DIM%KLON, 0:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FEFB2 (YDCPG_DIM%KLON, 0:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FEFB3 (YDCPG_DIM%KLON, 0:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FPLCH (YDCPG_DIM%KLON, 0:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FPLSH (YDCPG_DIM%KLON, 0:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FPLSN (YDCPG_DIM%KLON, 0:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FP (YDCPG_DIM%KLON, 0:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FTKEI (YDCPG_DIM%KLON, 0:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZPFL_FTKE (YDCPG_DIM%KLON, 0:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZTDS_TDALBNS (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDRHONS (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDSNS (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDTP (YDCPG_DIM%KLON, 1:YDCPG_DIM%YRSURF%YSP_SBD%NLEVS)
REAL (KIND=JPRB)     :: ZTDS_TDTS (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDWL (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDWPI (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDWP (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDWSI (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDWS (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZPRC_DPRECIPS2 (YDCPG_DIM%KLON, 1:YDCPG_DIM%KDTPREC2)
REAL (KIND=JPRB)     :: ZPRC_DPRECIPS (YDCPG_DIM%KLON, 1:YDCPG_DIM%KDTPREC)
REAL (KIND=JPRB)     :: ZRDG_CVGQ (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDG_LCVQ (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDG_MU0LU (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZRDG_MU0M (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZRDG_MU0N (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZRDG_MU0 (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZSAV_DDAL (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZSAV_DDOM (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZSAV_ENTCH (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZSAV_FHPS (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZSAV_GZ0F (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZSAV_GZ0HF (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZSAV_HV (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZSAV_PBLH (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZSAV_QSH  (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZSAV_UDAL (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZSAV_UDGRO (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZSAV_UDOM (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZSAV_UNEBH (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "accldia.intfb.h"
#include "acclph.intfb.h"
#include "acdnshf.intfb.h"
#include "acdrag.intfb.h"
#include "acdrme.intfb.h"
#include "acevadcape.intfb.h"
#include "achmt.intfb.h"
#include "achmtls.intfb.h"
#include "acpluis.intfb.h"
#include "acsol.intfb.h"
#include "actqsat.intfb.h"
#include "acvisih.intfb.h"
#include "aplpar_init.intfb.h"
#include "aro_ground_diag_z0.h"
#include "checkmv.intfb.h"
#include "cpchet.intfb.h"
#include "cpmvvps.intfb.h"
#include "cpnudg.intfb.h"
#include "cpphinp.intfb.h"
#include "cpqsol.intfb.h"
#include "cptend_new.intfb.h"
#include "cptends.intfb.h"
#include "cputqy_aplpar_expl.intfb.h"
#include "cpwts.intfb.h"
#include "dprecips.intfb.h"
#include "mf_phys_bayrad.intfb.h"
#include "mf_phys_corwat.intfb.h"
#include "mf_phys_cvv.intfb.h"
#include "mf_phys_fpl_part1.intfb.h"
#include "mf_phys_fpl_part2.intfb.h"
#include "mf_phys_mocon.intfb.h"
#include "mf_phys_precips.intfb.h"
#include "mf_phys_save_phsurf_part1.intfb.h"
#include "mf_phys_save_phsurf_part2.intfb.h"
#include "mf_phys_transfer.intfb.h"
#include "ppwetpoint.intfb.h"
#include "profilechet.intfb.h"
#include "qngcor.intfb.h"
#include "radozcmf.intfb.h"
#include "suozon.intfb.h"
#include "aplpar_flexdia.intfb.h"

!     ------------------------------------------------------------------

#include "fcttrm.func.h"
#include "apl_arpege_oceanic_fluxes.intfb.h"
#include "apl_wind_gust.intfb.h"
#include "apl_arpege_shallow_convection_and_turbulence.intfb.h"
#include "apl_arpege_albedo_computation.intfb.h"
#include "apl_arpege_aerosols_for_radiation.intfb.h"
#include "apl_arpege_cloudiness.intfb.h"
#include "apl_arpege_radiation.intfb.h"
#include "apl_arpege_soil_hydro.intfb.h"
#include "apl_arpege_surface.intfb.h"
#include "apl_arpege_deep_convection.intfb.h"
#include "apl_arpege_precipitation.intfb.h"
#include "apl_arpege_hydro_budget.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('APL_ARPEGE', 0, ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDVAB=>YDGEOMETRY%YRVAB, YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, YDTOPH=>YDMODEL%YRML_PHY_MF%YRTOPH, &
& YDRIP=>YDMODEL%YRML_GCONF%YRRIP, YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY, YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY,                       &
& YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, YGFL=>YDMODEL%YRML_GCONF%YGFL,                            &
& YDEPHY=> YDMODEL%YRML_PHY_EC%YREPHY, YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR, YDPRECIPS=>YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS,     &
& YDGEM=>YDGEOMETRY%YRGEM, YDSTA=>YDGEOMETRY%YRSTA, YDERDI=>YDMODEL%YRML_PHY_RAD%YRERDI, YDMCC=>YDMODEL%YRML_AOC%YRMCC,           &
& YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD, YDPHY3=>YDMODEL%YRML_PHY_MF%YRPHY3, YDPHY1=>YDMODEL%YRML_PHY_MF%YRPHY1,                    &
& YDPHY0=>YDMODEL%YRML_PHY_MF%YRPHY0, YDMSE=>YDMODEL%YRML_PHY_MF%YRMSE)


ASSOCIATE(TSPHY=>YDPHY2%TSPHY, NTSSG=>YDDPHY%NTSSG, LMSE=>YDARPHY%LMSE, NGFL_EXT=>YGFL%NGFL_EXT, YSP_SBD=>YDSURF%YSP_SBD,   &
& LNEBN=>YDPHY%LNEBN, NDPSFI=>YDPHY%NDPSFI, LRRGUST=>YDPHY%LRRGUST, LEDR=>YDPHY%LEDR, NTPLUI=>YDTOPH%NTPLUI,                &
& NDTPRECCUR=>YDPRECIPS%NDTPRECCUR, NDTPRECCUR2=>YDPRECIPS%NDTPRECCUR2, XMINLM=>YDPHY0%XMINLM, GRSO=>YDPHY0%GRSO,           &
& GAEPS=>YDPHY0%GAEPS, AERCS1=>YDPHY0%AERCS1, AERCS3=>YDPHY0%AERCS3, AERCS5=>YDPHY0%AERCS5, HUTIL2=>YDPHY0%HUTIL2,          &
& HUTIL1=>YDPHY0%HUTIL1, XMAXLM=>YDPHY0%XMAXLM, HUCOE=>YDPHY0%HUCOE, HUTIL=>YDPHY0%HUTIL, NPCLO1=>YDPHY0%NPCLO1,            &
& NPCLO2=>YDPHY0%NPCLO2, RDECRD=>YDPHY0%RDECRD, RDECRD1=>YDPHY0%RDECRD1, RDECRD2=>YDPHY0%RDECRD2, RDECRD3=>YDPHY0%RDECRD3,  &
& RDECRD4=>YDPHY0%RDECRD4, HSOLIWR=>YDPHY1%HSOLIWR, ALCRIN=>YDPHY1%ALCRIN, ALBMED=>YDPHY1%ALBMED, WSMX=>YDPHY1%WSMX,        &
& HSOLIT0=>YDPHY1%HSOLIT0, HSOL=>YDPHY1%HSOL, WPMX=>YDPHY1%WPMX, EMCRIN=>YDPHY1%EMCRIN, EMMMER=>YDPHY1%EMMMER,              &
& TMERGL=>YDPHY1%TMERGL, EMMGLA=>YDPHY1%EMMGLA, LRAFTKE=>YDPHY2%LRAFTKE, HVCLS=>YDPHY2%HVCLS, HTCLS=>YDPHY2%HTCLS,          &
& FSM_HH=>YDPHY3%FSM_HH, FSM_GG=>YDPHY3%FSM_GG, FSM_FF=>YDPHY3%FSM_FF, FSM_EE=>YDPHY3%FSM_EE, FSM_II=>YDPHY3%FSM_II,        &
& FSM_CC=>YDPHY3%FSM_CC, FSM_DD=>YDPHY3%FSM_DD, RII0=>YDPHY3%RII0, QCO2=>YDPHY3%QCO2, NDLUNG=>YDDIM%NDLUNG,                 &
& NDGUNG=>YDDIM%NDGUNG, NDLUXG=>YDDIM%NDLUXG, NDGUXG=>YDDIM%NDGUXG, LMPA=>YDARPHY%LMPA, CCOUPLING=>YDARPHY%CCOUPLING,       &
& LSURFEX_KFROM=>YDARPHY%LSURFEX_KFROM, YA=>YGFL%YA, YIRAD=>YGFL%YIRAD, YLRAD=>YGFL%YLRAD, XZSEPS=>YDMSE%XZSEPS,            &
& LAEROSOO=>YDPHY%LAEROSOO, LRKCDEV=>YDPHY%LRKCDEV, LHUCN=>YDPHY%LHUCN, LAEROLAN=>YDPHY%LAEROLAN, LAERODES=>YDPHY%LAERODES, &
& LNODIFQC=>YDPHY%LNODIFQC, NCALLRAD=>YDPHY%NCALLRAD, LO3ABC=>YDPHY%LO3ABC, LSTRAS=>YDPHY%LSTRAS, LAEROSEA=>YDPHY%LAEROSEA, &
& NDIFFNEB=>YDPHY%NDIFFNEB, LSOLV=>YDPHY%LSOLV, LRNUEXP=>YDPHY%LRNUEXP, LDAYD=>YDPHY%LDAYD, RDECLI=>YDRIP%RDECLI,           &
& XSW_BANDS=>YDPARAR%XSW_BANDS, NSWB_MNH=>YDPARAR%NSWB_MNH, NTNEBU=>YDTOPH%NTNEBU, NTDRME=>YDTOPH%NTDRME,                   &
& NTCVIM=>YDTOPH%NTCVIM, NTCOET=>YDTOPH%NTCOET, NTDRAG=>YDTOPH%NTDRAG, RMESOQ=>YDTOPH%RMESOQ, RMESOT=>YDTOPH%RMESOT,        &
& RMESOU=>YDTOPH%RMESOU, NTQSAT=>YDTOPH%NTQSAT, NTCOEF=>YDTOPH%NTCOEF, NAER=>YDERAD%NAER, NOZOCL=>YDERAD%NOZOCL,            &
& NRADFR=>YDERAD%NRADFR, NSW=>YDERAD%NSW, RCARDI=>YDERDI%RCARDI, RSUNDUR=>YDERDI%RSUNDUR, LMCC03=>YDMCC%LMCC03,             &
& NSTOP=>YDRIP%NSTOP, RCODEC=>YDRIP%RCODEC, RHGMT=>YDRIP%RHGMT, RSIDEC=>YDRIP%RSIDEC, RSOVR=>YDRIP%RSOVR,                   &
& RSTATI=>YDRIP%RSTATI, TSTEP=>YDRIP%TSTEP, STPRE=>YDSTA%STPRE, STPREH=>YDSTA%STPREH, STTEM=>YDSTA%STTEM,                   &
& LGCHECKMV=>YDPHY%LGCHECKMV, NAERO=>YGFL%NAERO, RDELXN=>YDGEM%RDELXN, LFLASH =>YDCFU%LFLASH)

ASSOCIATE(  PAPRSF=> YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, PAPRS => YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,                                  &
& PAPHIF=> YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, PAPHI => YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, PDELP => YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
& PR    => YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, PT    => YDMF_PHYS_BASE_STATE%T, PU    => YDMF_PHYS_BASE_STATE%U,                           &
& PV    => YDMF_PHYS_BASE_STATE%V, PCP   => YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, PTKE  => YDMF_PHYS_BASE_STATE%TKE,                        &
& PVERVEL => YDCPG_DYN0%CTY%VVEL  )


!     ------------------------------------------------------------------

!        0.    constructor for procset


!        1.    Preliminary calculations necessary
!              for all types of physics.
!              ------------------------------------


INSTEP_DEB=1
INSTEP_FIN=1

! SPP 
IF ( YSPP_CONFIG%LSPP ) THEN
 DO JSPP=1,YSPP%N2D
   ZGP2DSPP(:,JSPP) = YSPP%GP_ARP(JSPP)%GP2D(:,1,YDCPG_DIM%KBL)
 ENDDO
ENDIF

CALL CPPHINP(YDGEOMETRY, YDMODEL, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDVARS%GEOMETRY%GEMU%T0, YDVARS%GEOMETRY%GELAM%T0,   &
& YDVARS%U%T0, YDVARS%V%T0, YDVARS%Q%T0, YDVARS%Q%DL, YDVARS%Q%DM, YDVARS%CVGQ%DL, YDVARS%CVGQ%DM, YDCPG_PHY0%XYB%RDELP, &
& YDCPG_DYN0%CTY%EVEL, YDVARS%CVGQ%T0, ZRDG_MU0, ZRDG_MU0LU, ZRDG_MU0M, ZRDG_MU0N, ZRDG_CVGQ)
ZRDG_LCVQ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=ZRDG_CVGQ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)

DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
  ZFLU_QSATS(JROF)=0.0_JPRB
ENDDO

CALL MF_PHYS_FPL_PART1 (YDCPG_DIM, ZPFL_FPLCH, ZPFL_FPLSH, YDVARS%CPF%T0, YDVARS%SPF%T0, YDMODEL)


! * In some cases, some pseudo-historic surface buffers (like z0) should
!   not be modified between the entrance and the output of APLPAR
!   (this is the case for example if LDCONFX=T).
!   For the time being, we must save:
!   - HV (group VV) : resistance to evapotranspiration
!   - Z0F (group VD): gravity * surface roughness length
!   - Z0H (group VV): gravity * roughness length for heat
!   - PBLH (group VH): PBL height
!   - SPSH (group VH):
!   - QSH (group VH):

LL_SAVE_PHSURF=LDCONFX
IF (LL_SAVE_PHSURF) THEN
  CALL MF_PHYS_SAVE_PHSURF_PART1 (YDCPG_DIM, ZSAV_DDAL, ZSAV_DDOM, ZSAV_ENTCH, ZSAV_FHPS, ZSAV_GZ0F,                    &
  & ZSAV_GZ0HF, ZSAV_HV, ZSAV_PBLH, ZSAV_QSH, ZSAV_UDAL, ZSAV_UDGRO, ZSAV_UDOM, ZSAV_UNEBH, YDMF_PHYS_SURF%GSD_VF%PZ0F, &
  & YDMF_PHYS_SURF%GSD_VH%PPBLH, YDMF_PHYS_SURF%GSD_VH%PQSH, YDMF_PHYS_SURF%GSD_VH%PSPSH, YDMF_PHYS_SURF%GSD_VK%PUDGRO, &
  & YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VV%PZ0H, YDVARS%DAL%T0, YDVARS%DOM%T0, YDVARS%UAL%T0,                 &
  & YDVARS%UEN%T0, YDVARS%UNEBH%T0, YDVARS%UOM%T0, YDSURF, YDMODEL)
ENDIF




CALL APLPAR_INIT (YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, NTSSG, YSP_SBD%NLEVS, &
& YDMF_PHYS_SURF%GSD_VF%PVEG, ZMSC_FRMQ, ZDSA_CPS, ZDSA_LHS, ZDSA_RS, ZMSC_LH, ZMSC_LSCPE, ZFLU_QSAT,       &
& ZMSC_QW, ZMSC_TW, ZFLU_CD, ZFLU_CDN, ZFLU_CH, ZDSA_C1, ZDSA_C2, ZFLU_EMIS, ZFLU_FEVI, ZPFL_FTKE,          &
& ZPFL_FTKEI, ZPFL_FEFB1, ZPFL_FEFB2, ZPFL_FEFB3, ZFLU_NEIJ, ZFLU_VEG, ZFLU_QSATS, IMOC_CLPH)

!*       2.    Complete physics.
!              -----------------

!        2.2  Complete physics.
!             -----------------

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

! CALL PARAMETERISATIONS

DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
  YDMF_PHYS_SURF%GSD_VF%PLSM(JROF)=REAL(NINT(YDMF_PHYS_SURF%GSD_VF%PLSM(JROF)),JPRB)
ENDDO

IF (LTWOTL) THEN
  
ELSE
   ! IF (LAJUCV) THEN
   !   missing code under LAJUCV for leap-frog schemes.
   ! ENDIF
ENDIF


!
!-------------------------------------------------
! Check magnitude of model variables.
!-------------------------------------------------
!
IF(LGCHECKMV) CALL CHECKMV(YDRIP, YDPHY0, YDPHY2, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, &
              & YDCPG_DIM%KSTEP, PAPHI, PAPHIF, PAPRS, PAPRSF, YDVARS%GEOMETRY%GELAM%T0, YDVARS%GEOMETRY%GEMU%T0,     &
              & ZRDG_MU0, YDMF_PHYS_SURF%GSD_VF%PLSM, PT, YDMF_PHYS_BASE_STATE%Q, YDMF_PHYS_BASE_STATE%YGSP_RR%T      &
              &                       )
!     ------------------------------------------------------------------

LLREDPR=.FALSE.
ZRVMD=RV-RD
! SURFEX  and passive scalar
IF (LDCONFX) THEN
  ZDTMSE=0.01_JPRB
  ZSTATI=RSTATI-ZDTMSE/2._JPRB
ELSE
  ZDTMSE=TSPHY
  ZSTATI=RSTATI
ENDIF
ZRHGMT=REAL(RHGMT,JPRB)
ZAIPCMT(:)=0._JPRB


!*        1.0 DECORRELATION DEPTH FOR CLOUD OVERLAPS
IF ( RDECRD <= 0._JPRB .OR. LRNUEXP ) THEN
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZDECRD(JLON)=RDECRD1+RDECRD2* &
     & EXP(-((ASIN(YDVARS%GEOMETRY%GEMU%T0(JLON))-RDECRD3*RDECLI)/RDECRD4)**2)
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


  IF (YDMODEL%YRML_PHY_G%YRDPHY%NVCLIS == 1) THEN
    ZEPSO3=1.E-11_JPRB
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZQO3(JLON,0)=1.E-9_JPRB
    ENDDO
    DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZQO3(JLON,JLEV) = YDCPG_MISC%KOZO(JLON,JLEV,1)
      ENDDO
    ENDDO

  ELSEIF (NOZOCL == 1) THEN
    IF (MOD(YDCPG_DIM%KSTEP,NRADFR) == 0) THEN
      CALL RADOZCMF(YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, PAPRS, YDVARS%GEOMETRY%GEMU%T0, &
      & ZROZ)
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZQO3(JLON,0)=1.E-9_JPRB
      ENDDO
      DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZQO3(JLON,JLEV)=ZROZ(JLON,JLEV)*YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

  ELSE
    CALL SUOZON(YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, ZQO3, .FALSE., &
    & PAPRS, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, LO3ABC, YDMF_PHYS_SURF%GSD_VC%PGROUP)
  ENDIF

!     GAZ CARBONIQUE.

  DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZQCO2(JLON,JLEV)=QCO2
    ENDDO
  ENDDO


  ! INITIALISATION DE LA COORDONNEE ETA.
  ! INITIALISATION DE LA COORDONNEE ETA.

  ZVETAH(YDCPG_DIM%KTDIA-1)=STPREH(YDCPG_DIM%KTDIA-1)/VP00
  DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
    ZVETAH(JLEV)=STPREH(JLEV)/VP00
    ZVETAF(JLEV)=STPRE (JLEV)/VP00
  ENDDO

!     EPAISSEUR STD AEROSOLS
  ZAEO = AERCS1*ZVETAH(YDCPG_DIM%KTDIA-1) + AERCS3*ZVETAH(YDCPG_DIM%KTDIA-1)**3&
   & + AERCS5*ZVETAH(YDCPG_DIM%KTDIA-1)**5
  DO JLEV = YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG
    ZAEN = AERCS1*ZVETAH(JLEV) + AERCS3*ZVETAH(JLEV)**3&
     & + AERCS5*ZVETAH(JLEV)**5
    ZDAER(JLEV) = ZAEN - ZAEO
    ZAEO = ZAEN
  ENDDO

!     HUMIDITE CRITIQUE
  ZEPS=1.E-12_JPRB
  IF(LHUCN) THEN
    DO JLEV = YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG
      ZHUC(JLEV)=1.0_JPRB-MAX( HUCOE*ZVETAF(JLEV)*(1.0_JPRB-ZVETAF(JLEV))/((&
      & 1.0_JPRB+HUTIL1*(ZVETAF(JLEV)-0.5_JPRB))*(&
      & 1.0_JPRB+HUTIL2*(ZVETAF(JLEV)-0.5_JPRB))),ZEPS)
    ENDDO
  ELSE
    DO JLEV = YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG
      ZHUC(JLEV)=1.0_JPRB-MAX( HUCOE*ZVETAF(JLEV)**NPCLO1*(&
       & 1.0_JPRB-ZVETAF(JLEV))**NPCLO2*(&
       & 1.0_JPRB+SQRT(HUTIL)*(ZVETAF(JLEV)-0.5_JPRB)),ZEPS)
    ENDDO
  ENDIF

  

    DO JLEV = YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG-1
      ZMAN(JLEV) = FSM_CC * TANH(FSM_DD*ZVETAH(JLEV))
      ZMAK(JLEV) = FSM_EE * ZVETAH(JLEV)**FSM_FF +&
       & FSM_GG * (1-ZVETAH(JLEV))**FSM_HH + FSM_II
    ENDDO

!     INCREMENTAL CORRECTION FLUX FOR NEGAVTIVE HUMIDITY VALUES

  ZFCQVNG(:,:)=0.0_JPRB
  ZFCQING(:,:)=0.0_JPRB
  ZFCQLNG(:,:)=0.0_JPRB

  ZPRODTH_CVPP(:,:)=0.0_JPRB

  ZGDT=RG*TSPHY
  ZGDTI=1.0_JPRB/ZGDT

  DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA

!     LAYER WEIGHTS

      ZPOID(JLON,JLEV)=PDELP(JLON,JLEV)*ZGDTI
      ZIPOI(JLON,JLEV)=YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP(JLON,JLEV)*ZGDT

!     CALCULATION OF LATENT HEATS

      ZLHV(JLON,JLEV)=FOLH(PT(JLON,JLEV),0.0_JPRB)
      ZLHS(JLON,JLEV)=FOLH(PT(JLON,JLEV),1.0_JPRB)

    ENDDO
  ENDDO


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
    ZQI(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_BASE_STATE%I(JLON,JLEV))
    ZQL(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_BASE_STATE%L(JLON,JLEV))
    ZQR(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_BASE_STATE%R(JLON,JLEV))
    ZQS(JLON,JLEV)=MAX(0.0_JPRB,YDMF_PHYS_BASE_STATE%S(JLON,JLEV))

    ! CORRECTION OF NEGATIVE ADVECTED VALUES:
    !                         VAPOUR PUT IN PFCQVNG
    !                         LIQUID,ICE PUT IN PFCQL/ING
    !                         FOR RAIN/SNOW PUT IN PFCQR/SNG

    ZDQI=ZQI(JLON,JLEV)-YDMF_PHYS_BASE_STATE%I(JLON,JLEV)
    ZDQL=ZQL(JLON,JLEV)-YDMF_PHYS_BASE_STATE%L(JLON,JLEV)
    ZDQR=ZQR(JLON,JLEV)-YDMF_PHYS_BASE_STATE%R(JLON,JLEV)
    ZDQS=ZQS(JLON,JLEV)-YDMF_PHYS_BASE_STATE%S(JLON,JLEV)
    ZDQC=ZDQI+ZDQL+ZDQR+ZDQS

    ZQV0=YDMF_PHYS_BASE_STATE%Q(JLON,JLEV)-ZIPOI(JLON,JLEV)*(0.0_JPRB- YDMF_PHYS%OUT%FCQNG(JLON,JLEV-1)&
    & -YDMF_PHYS%OUT%FCQNNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQLNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQRNG(JLON,JLEV-1)&
    & -YDMF_PHYS%OUT%FCQSNG(JLON,JLEV-1))
    ZQV(JLON,JLEV)=MAX(0.0_JPRB,ZQV0-ZDQC)
    ZDQV=MAX(0.0_JPRB,ZQV0-ZDQC)-YDMF_PHYS_BASE_STATE%Q(JLON,JLEV)

    YDMF_PHYS%OUT%FCQNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNG(JLON,JLEV-1)-ZDQV*ZPOID(JLON,JLEV)
    YDMF_PHYS%OUT%FCQNNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNNG(JLON,JLEV-1)-ZDQI*ZPOID(JLON,JLEV)
    YDMF_PHYS%OUT%FCQLNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQLNG(JLON,JLEV-1)-ZDQL*ZPOID(JLON,JLEV)
    YDMF_PHYS%OUT%FCQRNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQRNG(JLON,JLEV-1)-ZDQR*ZPOID(JLON,JLEV)
    YDMF_PHYS%OUT%FCQSNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQSNG(JLON,JLEV-1)-ZDQS*ZPOID(JLON,JLEV)
  ENDDO
ENDDO
    
DO JCHA = 1, 6
  DO JLEV = YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG
    DO JLON = 1, YDCPG_DIM%KLON
      ZAER(JLON,JLEV,JCHA)=0.0_JPRB
    ENDDO
  ENDDO
ENDDO

DO JLEV = YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG
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
  IMOC_CLPH   (JLON)   = YDCPG_DIM%KFLEVG
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
  DO JLEV=YDCPG_DIM%KTDIA-1,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZPFL_FPLCH(JLON,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF
IF(LRRGUST.AND.YDCPG_DIM%KSTEP == 0) THEN
  DO JLEV=YDCPG_DIM%KTDIA-1,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZPFL_FPLSH(JLON,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF





!*
!     ------------------------------------------------------------------
!     4.- CALCULS THERMODYNAMIQUES
!     ----------------------------
  
CALL ACTQSAT ( YDPHY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTQSAT, YDCPG_DIM%KFLEVG, &
& PAPRSF, PCP, ZQV, PT, ZGEOSLC, ZMSC_LH, ZMSC_LSCPE, ZFLU_QSAT, ZMSC_QW, YDCPG_MISC%RH, ZMSC_TW)
  

!*
!     ------------------------------------------------------------------
!     4.BIS. COEFFICIENTS THERMO-HYDRIQUES DU SOL
!     -------------------------------------------

IF ( .NOT.LMSE ) THEN
  IF ( LSOLV ) THEN
    LLHMT=.FALSE.
    CALL ACSOL ( YDPHY, YDPHY1, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDMF_PHYS_SURF%GSD_VV%PARG,                      &
    & YDMF_PHYS_SURF%GSD_VV%PD2, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H, YDMF_PHYS_SURF%GSD_VF%PZ0RLF,             &
    & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PIVEG, YDMF_PHYS_SURF%GSD_VV%PLAI, YDMF_PHYS_BASE_STATE%YGSP_SG%A,         &
    & YDMF_PHYS_BASE_STATE%YGSP_SG%R, YDMF_PHYS_SURF%GSD_VV%PSAB, YDMF_PHYS_BASE_STATE%YGSP_SG%F, YDMF_PHYS_BASE_STATE%YGSP_RR%T,  &
    & YDMF_PHYS_SURF%GSD_VF%PVEG, YDMF_PHYS_BASE_STATE%YGSP_SB%Q, YDMF_PHYS_BASE_STATE%YGSP_SB%TL, YDMF_PHYS_BASE_STATE%YGSP_RR%W, &
    & YDMF_PHYS_BASE_STATE%YGSP_RR%IC, LLHMT, ZDSA_C1, ZDSA_C2, ZC3, ZCG, ZCN, YDMF_PHYS%OUT%CT,                                   &
    & ZNEIJG, ZNEIJV, ZWFC, ZWPMX, ZWSEQ, ZWSMX, ZWWILT)
  ELSE

!            INITIALISATION DE L'INERTIE THERMIQUE DU SOL.

!DEC$ IVDEP
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%CT(JLON)=HSOL /(&
      & 1.0_JPRB+HSOLIWR*(YDMF_PHYS_BASE_STATE%YGSP_RR%W(JLON)+YDMF_PHYS_BASE_STATE%YGSP_SB%Q(JLON,1))/(WSMX+WPMX)&
      & *EXP(-0.5_JPRB*(HSOLIT0*(YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)-RTT))**2))
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
    YDMF_PHYS%OUT%GZ0     (JLON) = YDCPG_GPAR%GZ0(JLON)
    YDMF_PHYS%OUT%GZ0H    (JLON) = YDCPG_GPAR%GZ0H(JLON)
    ZFLU_EMIS    (JLON) = YDCPG_GPAR%VEMIS(JLON)
    YDCPG_MISC%QS      (JLON) = YDCPG_GPAR%VQS(JLON)
    YDMF_PHYS_BASE_STATE%YGSP_RR%T      (JLON) = YDCPG_GPAR%VTS(JLON)
    ZSRAIN   (JLON) = YDCPG_GPAR%RAIN(JLON)
    ZSSNOW   (JLON) = YDCPG_GPAR%SNOW(JLON)
    ZSGROUPEL(JLON) = 0._JPRB
    ZTSN     (JLON) = YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)
  ENDDO

    
      ! FMR  radiation => NSW solar bands
  DO JSG=1,NSW
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZALBP(JLON,JSG) = YDCPG_GPAR%ALBDIR(JLON,JSG)
      ZALBD(JLON,JSG) = YDCPG_GPAR%ALBSCA(JLON,JSG)
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
    IF (ZFLU_EMIS(JLON)==0._JPRB) THEN
      ZFLU_EMIS(JLON) = 0.99_JPRB
      YDMF_PHYS_BASE_STATE%YGSP_RR%T  (JLON) = 288.0_JPRB
      YDMF_PHYS%OUT%ALB (JLON) = 0.1_JPRB
    ENDIF
  ENDDO

ENDIF  ! LMSE
!
! Define z0;z0h if it's necessary
!
IF (LMSE.AND.YDCPG_DIM%KSTEP == 0) THEN
  CALL ARO_GROUND_DIAG_Z0( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP, YDCPG_DIM%KFDIA-YDCPG_DIM%KIDIA+1, YDCPG_DIM%KIDIA,                  &
  & YDCPG_DIM%KFDIA, NDGUNG, NDGUXG, NDLUNG, NDLUXG, YDVARS%GEOMETRY%RINDX%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                   &
  & YDVARS%GEOMETRY%RINDY%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), LSURFEX_KFROM, YDMF_PHYS%OUT%GZ0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), &
  & YDMF_PHYS%OUT%GZ0H(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA))
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
  CALL ACHMTLS ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, &
  & PAPHI, PAPHIF, PAPRS, PAPRSF, PR, PT, PU, PV, YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDCPG_MISC%QS,          &
  & ZDPHIT, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, LLCLS, ZNBVNO, ZMRIPP, ZDSA_CPS, ZGWDCS, ZDSA_LHS,     &
  & ZPCLS, ZFLU_CD, ZFLU_CDN)
!       Computation of ZRTI
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZDPHI(JLON)=PAPHIF(JLON,YDCPG_DIM%KFLEVG)-PAPHI(JLON,YDCPG_DIM%KFLEVG)
    ZPRS(JLON)=RD+ZRVMD*YDCPG_MISC%QS(JLON)
    ZRTI(JLON)=2.0_JPRB/(PR(JLON,YDCPG_DIM%KFLEVG)*PT(JLON,YDCPG_DIM%KFLEVG)+RKAPPA*ZDPHI(JLON)&
    & +ZPRS(JLON)*YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON))
  ENDDO      
ELSE      
  CALL ACHMT ( YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YDSURF%YSD_VVD%NUMFLDS>=8.AND.LSOLV,         &
  & PAPHI, PAPHIF, PAPRS, PAPRSF, PCP, ZQV, PR, PT, PU, PV, ZPFL_FPLSH, ZPFL_FPLCH, ZDPHIT, ZDPHIV, YDMF_PHYS_SURF%GSD_VF%PZ0F, &
  & YDMF_PHYS_SURF%GSD_VV%PZ0H, YDMF_PHYS_SURF%GSD_VF%PZ0RLF, YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VF%PLSM,            &
  & ZNEIJG, ZNEIJV, YDMF_PHYS_BASE_STATE%YGSP_SG%F, YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS_SURF%GSD_VF%PVEG,                 &
  & ZWFC, YDMF_PHYS_BASE_STATE%YGSP_RR%W, YDMF_PHYS_BASE_STATE%YGSP_RR%IC, LLCLS, LLHMT, ZNBVNO,                                &
  & ZMRIPP, ZFLU_CD, ZFLU_CDN, ZCDROV, ZFLU_CH, ZCHROV, ZDSA_CPS, ZDQSTS, ZGWDCS, YDMF_PHYS%OUT%GZ0,                            &
  & YDMF_PHYS%OUT%GZ0H, ZHQ, ZHU, ZFLU_NEIJ, YDMF_PHYS%OUT%QCLS, YDCPG_MISC%QS, ZFLU_QSATS, YDMF_PHYS%OUT%RHCLS,                &
  & ZDSA_RS, ZRTI, ZSTAB, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS%OUT%NUCLS,                      &
  & YDMF_PHYS%OUT%NVCLS, ZPCLS, ZFLU_VEG, ZXDROV, ZXHROV, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST                                &
  &                              )      
ENDIF
    
ZCP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG) = PCP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)

    ! COMPUTATION OF 'DRY' mixing lengths : lm_d lh_d
    ! COMPUTATION OF ZPBLH - PBL HEIGHT

    
DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZTHETAV(JLON,JLEV)=PT(JLON,JLEV)*(1.0_JPRB+RETV*ZQV(JLON,JLEV))&
    & *(RATM/PAPRSF(JLON,JLEV))**RKAPPA  
  ENDDO
ENDDO
DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
  ZTHETAVS(JLON)=YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)*(1.0_JPRB+RETV*YDCPG_MISC%QS(JLON))&
  & *(RATM/PAPRS(JLON,YDCPG_DIM%KFLEVG))**RKAPPA  
ENDDO
CALL ACCLPH ( YDPHY0, YDPHY2, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG, &
& ZTHETAV, PAPHI, PAPHIF, PU, PV, ZTHETAVS, IMOC_CLPH, YDMF_PHYS%OUT%CLPH, YDMF_PHYS%OUT%VEIN, ZUGST,              &
& ZVGST)
ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS%OUT%CLPH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      
    

CALL APL_ARPEGE_OCEANIC_FLUXES (YDMF_PHYS_BASE_STATE, YDCPG_DIM, YDMF_PHYS, YDMF_PHYS_SURF, YDMODEL, LLHMT, ZCDROV, ZCE, &
                              & ZCEROV, ZCHROV, ZCRTI, ZDPHIT, ZDPHIV, ZDSA_RS, ZFLU_CD, ZFLU_CDN, ZFLU_CH, ZFLU_QSATS)


CALL APL_WIND_GUST (YDMF_PHYS_BASE_STATE, YDCPG_DIM, YDMF_PHYS, YDVARS, YDXFU, YDMODEL, IMOC_CLPH, ZBLH, ZDCAPE)


CALL APL_ARPEGE_SHALLOW_CONVECTION_AND_TURBULENCE (YDMF_PHYS_BASE_STATE, YDCPG_DIM, YDCPG_MISC, YDMF_PHYS, YDCPG_DYN0, &
& YDMODEL, YDDDH, INLAB_CVPP, INND, ZCDROV, ZCDROV_SAVE, ZCHROV, ZCHROV_SAVE, ZCOEFN, ZCONDCVPPI, ZCONDCVPPL, &
& ZDIFCVPPQ, ZDIFCVPPS, ZEPS, ZFLU_CD, ZFLU_CH, ZKQLROV, ZKQROV, ZKTROV, ZKTROV_SAVE, ZKUROV, ZKUROV_SAVE, ZMSC_LSCPE, &
& ZNBVNO, ZNEBS, ZNEBS0, ZNEB_CVPP, ZNLAB, ZNLABCVP, ZPFL_FPLCH, ZPFL_FTKE, ZPFL_FTKEI, ZPRODTH_CVPP, ZQI, ZQIC, ZQL, &
& ZQLC, ZQLIS, ZQLIS0, ZQLI_CVPP, ZQV, ZTKE1, ZTPRDY, ZXTROV, ZXUROV)

!**

!      7.1 Albedo et emissivite en presence de neige
!          Albedo and emissivity with snow

CALL APL_ARPEGE_ALBEDO_COMPUTATION (YDMF_PHYS_BASE_STATE, YDCPG_DIM, YDMF_PHYS, YDMF_PHYS_SURF, YDMODEL, ZALBD, ZALBP, ZEPS0, ZFLU_EMIS, ZFLU_NEIJ, ZMERL, ZRDG_MU0)

CALL APL_ARPEGE_AEROSOLS_FOR_RADIATION (YDMF_PHYS_BASE_STATE, YDCPG_DIM, YDMF_PHYS_SURF, YDMODEL, ZAEDES, ZAELAN, ZAER, ZAERINDS, ZAESEA, ZAESOO, ZAESUL, ZAEVOL)

CALL APL_ARPEGE_CLOUDINESS (YDMF_PHYS_BASE_STATE, YDCPG_DIM, YDCPG_MISC, YDMF_PHYS, YDVARS, &
& YDMODEL, LLREDPR, ZAIPCMT, ZBLH, ZCLCT_RAD, ZDECRD, ZFLU_QSAT, ZHUC, ZICEFR1, ZMSC_QW, &
& ZNEBC0, ZNEBCH, ZNEBS, ZNEBS0, ZNEB_CVPP, ZPFL_FPLCH, ZQI, ZQL, ZQLIS, ZQLIS0, ZQLI_CVP, &
& ZQLI_CVPP, ZQSATS, ZQV, ZRH, ZRHCRI, ZUNEBH, ZVETAF)
  
CALL APL_ARPEGE_RADIATION (YDMF_PHYS_BASE_STATE, YDGEOMETRY, YDCPG_DIM, YDCPG_MISC, YDCPG_GPAR, &
& YDMF_PHYS, YDMF_PHYS_SURF, YDVARS, YDMODEL, ZAER, ZAERINDS, ZALB, ZALBD, ZALBP, ZCEMTR, &
& ZCTRSO, ZDSRP, ZFLU_EMIS, ZFLU_QSAT, ZFRSODS, ZFSDNN, ZFSDNV, ZQO3, ZQR, ZQS, ZQV, &
& ZRDG_MU0, ZRDG_MU0LU, ZRDG_MU0M, ZSDUR, ZSFSWDIF, ZSFSWDIR, ZSUDU, ZTENT, ZTRSOD, &
& ZTRSODIF, ZTRSODIR)

!*
!     ------------------------------------------------------------------

!     7.BIS. BILAN HYDRIQUE DU SOL
!     ----------------------------
!     CALCUL DES RESISTANCES A L'EVAPOTRANSPIRATION HV ET
!     A LA TRANSPIRATION
!     ------------------------------------------------------------------
!     HTR DU COUVERT VEGETAL
!     ----------------------

CALL APL_ARPEGE_SOIL_HYDRO (YDMF_PHYS_BASE_STATE, YDCPG_DIM, YDMF_PHYS, YDMF_PHYS_SURF, YDMODEL, ZCHROV, ZFLU_NEIJ, ZFLU_QSAT, ZFLU_QSATS, ZFLU_VEG, ZGWDCS, ZHQ, ZHTR, ZHU, ZQV, ZWFC, ZWLMX, ZWWILT)

!     ------------------------------------------------------------------
!     8.- DIFFUSION VERTICALE TURBULENTE
!     ----------------------------------
  
CALL APL_ARPEGE_SURFACE (YDMF_PHYS_BASE_STATE, YDGEOMETRY, YDCPG_DIM, YDCPG_MISC, YDCPG_GPAR, YDMF_PHYS, &
& YDMF_PHYS_SURF, YDVARS, YDSURF, YDMODEL, LDCONFX, ZALBD, ZALBP, ZALPHA1, &
& ZCDROV, ZCEROV, ZCFATH, ZCFAU, ZCFBTH, &
& ZCFBU, ZCFBV, ZCHROV, ZCOEFA, ZCOEFN, ZCP, ZDIFEXT, ZDIFWQ, ZDIFWS, ZDQSTS, ZDSA_CPS, &
& ZDSA_LHS, ZDTMSE, ZEDMFQ, ZEDMFS, ZEDMFU, ZEDMFV, ZFLU_CD, ZFLU_CDN, &
& ZFLU_EMIS, ZFLU_FEVI, ZFLU_NEIJ, ZFLU_QSATS, ZFLU_VEG, ZHQ, ZHTR, ZHU, ZKQLROV, &
& ZKQROV, ZKTROV, ZKUROV, ZLVT, ZMF_UP, ZNEBCH, ZNEBDIFF, ZNEBS, ZPOID, ZQI, ZQICE, ZQL, ZQV, &
& ZRDG_MU0, ZRDG_MU0N, ZRHGMT, ZSC_FCLL, ZSC_FCLN, ZSC_FEVI, ZSC_FEVN, &
& ZSFSWDIF, ZSFSWDIR, ZSGROUPEL, ZSRAIN, ZSSNOW, ZSTATI, ZTSN, &
& ZXDROV, ZXHROV, ZXQRO, ZXTRO, ZXTROV, ZXURO, ZXUROV)
    

!-----------------------------------------------------------------------------
!   THE DEEP CONVECTION WILL SEE THE SHALLOW PART FROM KFB AS IT IS WITH LOUIS
!   SCHEME AND THE MODIFIED RI
!----------------------------------------------------------------------------
    
DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    YDMF_PHYS%OUT%DIFTQ   (JLON,JLEV) = YDMF_PHYS%OUT%DIFTQ(JLON,JLEV) + ZDIFCVPPQ(JLON,JLEV)
    YDMF_PHYS%OUT%DIFTS   (JLON,JLEV) = YDMF_PHYS%OUT%DIFTS(JLON,JLEV) + ZDIFCVPPS(JLON,JLEV)
    YDMF_PHYS%OUT%STRTU   (JLON,JLEV) = YDMF_PHYS%OUT%STRTU(JLON,JLEV) + ZDIFCVPPU(JLON,JLEV)
    YDMF_PHYS%OUT%STRTV   (JLON,JLEV) = YDMF_PHYS%OUT%STRTV(JLON,JLEV) + ZDIFCVPPV(JLON,JLEV)
  ENDDO
ENDDO
    


!     DIAGNOSTIC SUPPLEMENTAIRE FLUX DE RAYONNEMENT (PFRTHDS ET PFRSOPT)
!     ADDITIONAL DIAGNOSTICS OF RADIATIVE FLUXES (PFRTHDS AND PFRSOPT)

!DEC$ IVDEP
DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
  IF(.NOT.LMSE)YDMF_PHYS%OUT%FRTHDS(JLON)=YDMF_PHYS%OUT%FRTH(JLON,YDCPG_DIM%KFLEVG,1)/ZFLU_EMIS(JLON)+RSIGMA*YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)**4
  IF(ZRDG_MU0(JLON) <= 0.0_JPRB) THEN
    YDMF_PHYS%OUT%FRSOPT(JLON)=0.0_JPRB
  ELSE
    YDMF_PHYS%OUT%FRSOPT(JLON)=RII0*ZRDG_MU0(JLON)
  ENDIF
ENDDO

! ADDITIONAL DIAGNOSTIC OF THE DERIVATIVE OF THE NON SOLAR SURFACE
! HEAT FLUX WITH RESPECT TO SURFACE TEMPERATURE (PDERNSHF)

IF(LMCC03)THEN
  CALL ACDNSHF(YDPHY, YDPHY1, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG, &
  & ZFLU_EMIS, YDMF_PHYS_SURF%GSD_VF%PLSM, ZFLU_NEIJ, ZQV, YDCPG_MISC%QS, YDMF_PHYS_BASE_STATE%YGSP_RR%T,          &
  & ZCHROV, ZDQSTS, YDMF_PHYS%OUT%DRNSHF)
ENDIF

!*
!     ------------------------------------------------------------------
!     9.- TRAINEE DES ONDES DE GRAVITE INDUITES PAR LE RELIEF
!     ------------------------------------------------------------------

! GRAVITY WAVE DRAG  
CALL ACDRAG ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTDRAG, YDCPG_DIM%KFLEVG,           &
& PAPRS, PAPRSF, PDELP, ZNBVNO, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, PU, PV, YDVARS%GEOMETRY%RCORI%T0,               &
& YDMF_PHYS_SURF%GSD_VF%PGETRL, ZGWDCS, YDMF_PHYS_SURF%GSD_VF%PVRLAN, YDMF_PHYS_SURF%GSD_VF%PVRLDI, YDMF_PHYS%OUT%STRDU, &
& YDMF_PHYS%OUT%STRDV, ZTRAJGWD)
  
 ! SAVE FOR TL/NL COEFS FROM VERT. DIFF AND GWD

  

!     ------------------------------------------------------------------
!     10.- PRECIPITATIONS STRATIFORMES.
!     ---------------------------------

  
! ALARO PRECIPITATION SCHEME
IF ( LSTRAS ) THEN
  CALL ACPLUIS ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTPLUI, YDCPG_DIM%KFLEVG,     &
  & PAPRSF, PDELP, ZNEBS, ZQV, ZQLIS, ZMSC_QW, PR, PT, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLSL, &
  & YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN     &
  &        )
ENDIF
  
CALL APL_ARPEGE_DEEP_CONVECTION (YDMF_PHYS_BASE_STATE, YDGEOMETRY, YDCPG_DIM, YDMF_PHYS, YDCPG_DYN0, &
& YDMF_PHYS_SURF, YDVARS, YDCFU, YDMODEL, ZFPCOR, ZQV, ZTENHA, ZTENQVA, ZTENT)

CALL APL_ARPEGE_PRECIPITATION (YDMF_PHYS_BASE_STATE, YDCPG_DIM, YDMF_PHYS, YDMF_PHYS_SURF, YDVARS, YDMODEL, &
& ZFLU_NEIJ, ZNEBS, ZNEB_CVPP, ZQC_DET_PCMT, ZQI, ZQL, ZQLIS, ZQLI_CVPP, ZQR, ZQS, ZQV, ZSEDIQI, ZSEDIQL, &
& ZTENHA, ZTENQVA, ZVETAF)
    
ZFLX_LOTT_GWU(:, :) = 0.0_JPRB
ZFLX_LOTT_GWV(:, :) = 0.0_JPRB

!*
!     ------------------------------------------------------------------
!         SAUVEGARDE DES FLUX DE PRECIPITATION CONVECTIVE ET STRATIFORME.

IF(LNEBN.OR.LRRGUST) THEN
  DO JLEV=YDCPG_DIM%KTDIA-1,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZPFL_FPLCH(JLON,JLEV)=ZFPCOR(JLON,JLEV)
    ENDDO
  ENDDO    
ENDIF

IF(LRRGUST) THEN
  DO JLEV=YDCPG_DIM%KTDIA-1,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZPFL_FPLSH(JLON,JLEV)=YDMF_PHYS%OUT%FPLSL(JLON,JLEV)+YDMF_PHYS%OUT%FPLSN(JLON,JLEV)+YDMF_PHYS%OUT%FPLSG(JLON,JLEV)
    ENDDO
  ENDDO
ENDIF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! UPDATE TRANSPORT FLUXES DUE TO SEDIMENTATION OF CLOUDS.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -

DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
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

  
CALL QNGCOR (YDPHY2, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTPLUI, YDCPG_DIM%KFLEVG, ZQV,                 &
& ZQL, ZQI, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, &
& YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN,      &
& YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPCL,     &
& YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN,           &
& YDMF_PHYS%OUT%FCQNG )

!*
!     ------------------------------------------------------------------
!     12. - BILAN HYDRIQUE DU SOL
!     ---------------------------

CALL APL_ARPEGE_HYDRO_BUDGET (YDMF_PHYS_BASE_STATE, YDCPG_DIM, YDCPG_GPAR, YDMF_PHYS, YDMF_PHYS_SURF, YDSURF, YDMODEL, ZC3, ZCN, ZDSA_C1, ZDSA_C2, ZFLU_FEVI, ZFLU_NEIJ, ZFLU_VEG, ZFPLSL, ZFPLSN, ZWFC, ZWLMX, ZWPMX, ZWSEQ, ZWSMX)

!*
!-  --------------------------------------------------------------------
!     13.- DRAG MESOSPHERIQUE POUR UN MODELE POSSEDANT DES NIVEAUX
!               AU-DESSUS DE 50 KM (I.E. DANS LA MESOSPHERE)
!     ------------------------------------------------------------------

! 
! MESOSPHERIC DRAG  
CALL ACDRME ( YDMODEL%YRCST, YDPHY2, YDTOPH, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,               &
& NTDRME, YDCPG_DIM%KFLEVG, PCP, PDELP, PT, ZQV, PU, PV, YDMF_PHYS%OUT%FRMH, ZMSC_FRMQ, YDMF_PHYS%OUT%STRMU, &
& YDMF_PHYS%OUT%STRMV, RMESOT, RMESOQ, RMESOU, STTEM)
  
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

IF (YDMODEL%YRML_DIAG%YRLDDH%LFLEXDIA) THEN
  CALL APLPAR_FLEXDIA (YDCPG_DIM, YDCPG_MISC, YDMF_PHYS, YDMF_PHYS_SURF, YDSURF, YDMODEL, YDDDH, &
  & YDMF_PHYS_BASE_STATE)
ENDIF

IF (LRAFTKE) THEN
     ! DCAPE due to precipitation evaporation.
  CALL ACEVADCAPE(YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YDMF_PHYS%OUT%FPLSL, &
  & YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, PAPRS, PAPRSF, PT, PCP, PAPHIF,         &
  & PAPHI, ZDCAPE)
     ! Gusts.
  ZQGM=ZEPSNEB
  CALL ACCLDIA(YDXFU, YDPHY, YDPHY2, YDTOPH, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,        &
  & YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, PU, PV, YDMF_PHYS%OUT%CAPE, ZDCAPE, ZTKE1, PAPHIF, YDVARS%GEOMETRY%OROG%T0, &
  & YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST, ZBLH, IMOC_CLPH)
  YDMF_PHYS%OUT%CLPH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=MIN(XMAXLM,MAX(XMINLM,ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)))
ENDIF

  
CALL ACVISIH(YDMODEL%YRML_PHY_MF%YRPHY%YRDVISI, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA,         &
& YDCPG_DIM%KFLEVG, PAPHI, PAPHIF, PAPRSF, PT, PR, YDCPG_MISC%QLI, YDCPG_MISC%QICE, ZQR, ZQS, ZQGM, YDMF_PHYS%OUT%VISICLD, &
& YDMF_PHYS%OUT%VISIHYD, YDMF_PHYS%OUT%MXCLWC    )
  

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

ZBAY_QRCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) =  0.0_JPRB
ZBAY_QSCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) =  0.0_JPRB

DO JLEV=1,YDCPG_DIM%KFLEVG
  DO JLON= YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
    ZBAY_QRCONV(JLON,JLEV)  = 0.001_JPRB * ( (ABS(YDMF_PHYS%OUT%FPLCL(JLON,JLEV)) * ZCOEFRAIN(1)) ** ZCOEFRAIN(2) )
    ZBAY_QSCONV(JLON,JLEV)  = 0.001_JPRB * ( (ABS(YDMF_PHYS%OUT%FPLCN(JLON,JLEV)) * ZCOEFSNOW(1)) ** ZCOEFSNOW(2) )    
  ENDDO
ENDDO

   ! Convert density [kg/m3] to mixing ratio [kg/kg]
   ! R_dry (dry air constant)

ZDE2MR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)  = RD * PT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) / PAPRSF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
ZBAY_QRCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) = ZBAY_QRCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) * ZDE2MR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
ZBAY_QSCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) = ZBAY_QSCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) * ZDE2MR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)

! Precipitation Type

! Compute wet-bulb temperature at 2 meters (suppose homogeneity of qv/ql/qi )
!ZPCLS(KIDIA:KFDIA)=PAPRS(KIDIA:KFDIA,KLEV)-2._JPRB/ZZZF(KIDIA:KFDIA,1,KLEV)*&
!                 &(PAPRS(KIDIA:KFDIA,KLEV)-PAPRSF(KIDIA:KFDIA,KLEV))
CALL PPWETPOINT(YDPHY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, ZPCLS, YDMF_PHYS%OUT%TCLS, YDMF_PHYS_BASE_STATE%Q(:, YDCPG_DIM%KFLEVG), &
& YDMF_PHYS_BASE_STATE%L(:, YDCPG_DIM%KFLEVG), YDMF_PHYS_BASE_STATE%I(:, YDCPG_DIM%KFLEVG), YDMF_PHYS%OUT%TPWCLS                                 &
& )                            

   ! Defined precipitation type
   !
   ! Compute wet-bulb temperature
DO JLEV=1,YDCPG_DIM%KFLEVG
  CALL PPWETPOINT(YDPHY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, PAPRSF(:, JLEV), PT(:, JLEV), &
  & YDMF_PHYS_BASE_STATE%Q(:, JLEV), YDMF_PHYS_BASE_STATE%L(:, JLEV), YDMF_PHYS_BASE_STATE%I(:, JLEV),   &
  & ZTW(:, JLEV))
ENDDO

DO JLON=1,YDCPG_DIM%KLON
  ZFPLS(JLON,YDCPG_DIM%KFLEVG)=YDMF_PHYS%OUT%FPLCN(JLON,YDCPG_DIM%KFLEVG)+YDMF_PHYS%OUT%FPLSN(JLON,YDCPG_DIM%KFLEVG)
  ZFPLC(JLON,YDCPG_DIM%KFLEVG)=YDMF_PHYS%OUT%FPLCL(JLON,YDCPG_DIM%KFLEVG)+YDMF_PHYS%OUT%FPLSL(JLON,YDCPG_DIM%KFLEVG)
  ZFPLSG(JLON,YDCPG_DIM%KFLEVG)=0._JPRB
ENDDO

  !initialisation de ZZZ
DO JLEV = 1,YDCPG_DIM%KFLEVG
  DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZZZ(JLON,1,JLEV)=PAPHI(JLON,JLEV)*ZINVG
  ENDDO
ENDDO

  !initialisation de ZDZZ
DO JLEV = 2, YDCPG_DIM%KFLEVG
  DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZDZZ(JLON,1,JLEV)=ZZZ(JLON,1,JLEV-1)-ZZZ(JLON,1,JLEV)
  ENDDO
ENDDO
DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
  ZDZZ(JLON,1,1)=PAPHI(JLON,0)*ZINVG-ZZZ(JLON,1,1)
ENDDO

CALL DPRECIPS (YDPHY%YRDPRECIPS, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,               &
& YDVARS%GEOMETRY%OROG%T0, YDMF_PHYS%OUT%TPWCLS, YDMF_PHYS%OUT%DIAGH, PAPHIF, ZDZZ, ZTW, YDMF_PHYS_BASE_STATE%L,   &
& ZFPLC(:, YDCPG_DIM%KFLEVG), ZFPLS(:, YDCPG_DIM%KFLEVG), ZFPLSG(:, YDCPG_DIM%KFLEVG), ZPRC_DPRECIPS(:, NDTPRECCUR)&
& )
  
 !Idem for an other time step and an other period

CALL DPRECIPS(YDPHY%YRDPRECIPS, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YDVARS%GEOMETRY%OROG%T0, &
& YDMF_PHYS%OUT%TPWCLS, YDMF_PHYS%OUT%DIAGH, PAPHIF, ZDZZ, ZTW, YDMF_PHYS_BASE_STATE%L, ZFPLC(:, YDCPG_DIM%KFLEVG),          &
& ZFPLS(:, YDCPG_DIM%KFLEVG), ZFPLSG(:, YDCPG_DIM%KFLEVG), ZPRC_DPRECIPS2(:, NDTPRECCUR2))

!        2.3  Computes MOCON in the CLP.
!             --------------------------
CALL MF_PHYS_MOCON (YDCPG_DIM, YDMF_PHYS%OUT%MOCON, ZRDG_LCVQ, IMOC_CLPH, YDMF_PHYS_BASE_STATE)

! Store surface water flux P and E for water conservation
IF (LCORWAT) THEN
  CALL MF_PHYS_CORWAT (YDCPG_DIM, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, &
  & YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS_SURF%GSD_VH%PEVA, YDMF_PHYS_SURF%GSD_VH%PPCL,               &
  & YDMF_PHYS_SURF%GSD_VH%PPCN, YDMF_PHYS_SURF%GSD_VH%PPSL, YDMF_PHYS_SURF%GSD_VH%PPSN)
ENDIF

!        2.6   surface specific humidity necessary to compute the vertical
!              advection of q in the case "delta m=1" (unlagged physics only).
!              ---------------------------------------------------------------

IF (NDPSFI == 1) THEN
  CALL CPQSOL(YDGEOMETRY%YRDIMV, YDPHY, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_PHY0%PREHYD, &
  & YDMF_PHYS_SURF%GSP_RR%PT_T0, YDCPG_MISC%QS, ZFLU_QSATS, YDCPG_MISC%QSOL)
ENDIF

!        2.7  Computation of tendencies T,u,v and Q.
!             --------------------------------------

! Set GFL tendencies to 0

ZTENDEFB1   = 0._JPRB
ZTENDEFB2   = 0._JPRB
ZTENDEFB3   = 0._JPRB
ZTENDG      = 0._JPRB
ZTENDI      = 0._JPRB
ZTENDQ      = 0._JPRB
ZTENDR      = 0._JPRB
ZTENDS      = 0._JPRB
ZTENDL      = 0._JPRB
ZTENDICONV  = 0._JPRB
ZTENDLCONV  = 0._JPRB
ZTENDRCONV  = 0._JPRB
ZTENDSCONV  = 0._JPRB
ZTENDTKE    = 0._JPRB

! * CPTEND+CPUTQY = Old( CPATY + CPDUP + CPDTHP )
! Calcul des tendances de T , U et de Q et modifications
! eventuelles de W et de OMEGA/P


CALL CPTEND_NEW( YDMODEL, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, YDVARS%GEOMETRY%GNORDL%T0, &
& YDVARS%GEOMETRY%GNORDM%T0, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS,       &
& ZDIFEXT, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%FCCQL,    &
& YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN,               &
& YDMF_PHYS%OUT%FPLSG, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLCG, YDMF_PHYS%OUT%FPEVPSL,             &
& YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPEVPCG,     &
& YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN,          &
& YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FCQGNG,          &
& YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%FRMH, ZMSC_FRMQ, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%STRCU,       &
& YDMF_PHYS%OUT%STRCV, YDMF_PHYS%OUT%STRDU, YDMF_PHYS%OUT%STRDV, YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV,               &
& YDMF_PHYS%OUT%STRMU, YDMF_PHYS%OUT%STRMV, YDMF_PHYS%OUT%DIFCQLC, YDMF_PHYS%OUT%DIFCQIC, YDMF_PHYS%OUT%FIMCC,           &
& YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC,        &
& YDMF_PHYS%OUT%FCNEGQIC, YDMF_PHYS%OUT%FCNEGQRC, YDMF_PHYS%OUT%FCNEGQSC, ZPFL_FTKE, ZPFL_FTKEI,                         &
& ZPFL_FEFB1, ZPFL_FEFB2, ZPFL_FEFB3, PDELP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, PAPHIF, PCP,                       &
& PU, PV, PT, YDMF_PHYS_BASE_STATE%Q, YDMF_PHYS_BASE_STATE%I, YDMF_PHYS_BASE_STATE%L, YDVARS%LCONV%T0,                   &
& YDVARS%ICONV%T0, YDVARS%RCONV%T0, YDVARS%SCONV%T0, YDMF_PHYS_BASE_STATE%R, YDMF_PHYS_BASE_STATE%S,                     &
& YDMF_PHYS_BASE_STATE%G, ZDSA_CPS, YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS%OUT%FHSCL, YDMF_PHYS%OUT%FHSCN,            &
& YDMF_PHYS%OUT%FHSSL, YDMF_PHYS%OUT%FHSSN, YDMF_PHYS%OUT%FHSSG, YDMF_PHYS%OUT%FHPCL, YDMF_PHYS%OUT%FHPCN,               &
& YDMF_PHYS%OUT%FHPCG, YDMF_PHYS%OUT%FHPSL, YDMF_PHYS%OUT%FHPSN, YDMF_PHYS%OUT%FHPSG, ZMSC_FHP,                          &
& ZPFL_FP, YDMF_PHYS%OUT%FEPFP, YDMF_PHYS%OUT%FCMPCQ, YDMF_PHYS%OUT%FCMPSN, YDMF_PHYS%OUT%FCMPSL, YDMF_PHYS%OUT%TENDU,   &
& YDMF_PHYS%OUT%TENDV, ZTENDU, ZTENDV, ZTENDH, ZTENDQ, ZTENDI, ZTENDL, ZTENDLCONV, ZTENDICONV,                           &
& ZTENDRCONV, ZTENDSCONV, ZTENDR, ZTENDS, ZTENDG, ZTENDTKE, ZTENDEFB1, ZTENDEFB2, ZTENDEFB3,                             &
& ZTENDEXT, YDDDH)


IF (LTWOTL) THEN
ELSE    
  IF ( (NDPSFI==1)) THEN
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
ENDIF 

!        2.7.1  Diagnostics on physical tendencies
!               ----------------------------------

IF (.NOT.LDCONFX) THEN
  IF ((GCHETN%LFREQD).OR.(GCHETN%LCOORD).OR.(GCHETN%LPROFV)) THEN
    CALL CPCHET (YDMF_PHYS, YDMF_PHYS_BASE_STATE, YDCPG_MISC, YDRIP, YDPHY, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, &
    & YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, YDCPG_DIM%KSTEP, ZMSC_FRMQ, ZDSA_CPS, ZTENDH, ZTENDQ,               &
    & ZTENDI, ZTENDL, ZTENDR, ZTENDS, YDVARS%GEOMETRY%GEMU%T0, YDVARS%GEOMETRY%GELAM%T0)
  ENDIF
  
  IF (GCHETN%LPROFV)&
   & CALL PROFILECHET(YDGEOMETRY, YDCPG_MISC, YDMF_PHYS, ZRDG_CVGQ, ZRDG_LCVQ, ZRDG_MU0, YDCPG_DYN0,                  &
     & YDMF_PHYS_SURF, YDVARS, YDSURF, YDDPHY, YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KFDIA, YDVARS%GEOMETRY%GELAM%T0, &
     & YDVARS%GEOMETRY%GEMU%T0, YDVARS%GEOMETRY%GM%T0, YDVARS%GEOMETRY%OROG%T0, YDVARS%GEOMETRY%RCORI%T0,             &
     & YDVARS%GEOMETRY%RATATH%T0, YDVARS%GEOMETRY%RATATX%T0)
ENDIF

!        2.8  Modification of vertical velocities
!             by some physics output when required.
!             -------------------------------------


! * MODIFICATION DE LA VITESSE VERTICALE ET DE LA TENDANCE DE
! PRESSION DE SURFACE SI NDPSFI=1 ( MASSE VARIABLE ).
! Ajout de la physique dans l'equation de continuite/Add physics
! in continuity equation.

IF (NDPSFI == 1) THEN
  CALL CPMVVPS(YDVAB, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, PDTPHY,     &
  & ZPFL_FP, PAPRS(:, YDCPG_DIM%KFLEVG), YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, YDCPG_DYN0%CTY%EVEL, &
  & YDCPG_DYN0%CTY%PSDVBC, YDMF_PHYS_NEXT_STATE%SP)
ENDIF

!        2.9  Computation of evolution of T, u, v and Q.
!             ------------------------------------------

! ky: non-zero option not yet coded for the time being.
ZTENDD=0.0_JPRB

! Calcul de T , Q et du Vent a l'instant 1

CALL CPUTQY_APLPAR_EXPL(YDCPG_DIM, YDMF_PHYS_NEXT_STATE, YDMF_PHYS_BASE_STATE, YDVARS, YDPHY, PDTPHY, &
& ZTENDH, ZTENDT, YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDU, ZTENDV, ZTENDD, ZTENDEFB1,        &
& ZTENDEFB2, ZTENDEFB3, ZTENDG, ZTENDICONV, ZTENDI, ZTENDLCONV, ZTENDL, ZTENDQ, ZTENDRCONV, ZTENDR,   &
& ZTENDSCONV, ZTENDS, ZTENDTKE, YDMF_PHYS%OUT%FDIS)

CALL MF_PHYS_FPL_PART2 (YDCPG_DIM, ZPFL_FPLCH, ZPFL_FPLSH, YDVARS%CPF%T1, YDVARS%SPF%T1, YDMODEL)

!       2.9b Prognostic convection etc.
!            --------------------------

! TRANSFER NOT ADVECTED VARIABLES INTO PGFLT1
CALL MF_PHYS_TRANSFER (YDCPG_DIM, YDVARS, YDMODEL)

!        2.10  Surface variables.
!              ------------------

IF (.NOT.LMSE) THEN
  DO JLEV=0,YDCPG_DIM%KFLEVG
    DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZPFL_FPLSN(JROF,JLEV)=YDMF_PHYS%OUT%FPLSN(JROF,JLEV)+YDMF_PHYS%OUT%FPLSG(JROF,JLEV)
    ENDDO
  ENDDO
  CALL CPTENDS( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG,                  &
  & YSP_SBD%NLEVS, PDTPHY, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLCN, ZPFL_FPLSN,                     &
  & YDMF_PHYS  %OUT%FRSO, YDMF_PHYS%OUT%FRTH, YDMF_PHYS_SURF%GSP_SG%PA_T1, YDMF_PHYS%OUT%CT, ZDSA_C1,                     &
  & ZDSA_C2, YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FCLL, YDMF_PHYS  %OUT%FCLN, YDMF_PHYS%OUT%FCS,                            &
  & ZFLU_FEVI, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT  %FEVN, YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FGEL, YDMF_PHYS%OUT%FGELS,     &
  & YDMF_PHYS%OUT%FLWSP, YDMF_PHYS%OUT%FONTE, YDMF_PHYS%OUT%FTR, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSP_SG%PR_T1, &
  & YDMF_PHYS%OUT%RUISL, YDMF_PHYS%OUT%RUISP, YDMF_PHYS%OUT%RUISS, YDMF_PHYS_SURF%GSP_SG%PF_T1,                           &
  & ZFLU_VEG, ZTDS_TDTS, ZTDS_TDTP, ZTDS_TDWS, ZTDS_TDWSI, ZTDS_TDWP, ZTDS_TDWPI, ZTDS_TDWL,                              &
  & ZTDS_TDSNS, ZTDS_TDALBNS, ZTDS_TDRHONS)  

  CALL CPWTS(YDSURF, YDMODEL%YRML_AOC%YRMCC, YDPHY, YDMODEL%YRML_PHY_MF%YRPHY1, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA,           &
  & YDCPG_DIM%KFDIA, YSP_SBD%NLEVS, PDTPHY, ZTDS_TDTS, ZTDS_TDTP, ZTDS_TDWS, ZTDS_TDWSI, ZTDS_TDWP,                        &
  & ZTDS_TDWPI, ZTDS_TDWL, ZTDS_TDSNS, ZTDS_TDALBNS, ZTDS_TDRHONS, YDMF_PHYS_SURF%GSD_VP%PTPC, YDMF_PHYS_SURF%GSD_VP%PWPC, &
  & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PIVEG, YDMF_PHYS_SURF%GSP_RR%PT_T1, YDMF_PHYS_SURF%GSP_SB%PT_T1,     &
  & YDMF_PHYS_SURF%GSP_RR%PW_T1, YDMF_PHYS_SURF%GSP_RR%PIC_T1, YDMF_PHYS_SURF%GSP_SB%PQ_T1, YDMF_PHYS_SURF%GSP_SB%PTL_T1,  &
  & YDMF_PHYS_SURF%GSP_RR%PFC_T1, YDMF_PHYS_SURF%GSP_SG%PF_T1, YDMF_PHYS_SURF%GSP_SG%PA_T1, YDMF_PHYS_SURF%GSP_SG%PR_T1    &
  & )  
ELSE
  IF (LDCONFX) THEN
    DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS_SURF%GSP_RR%PT_T0(JROF)=YDCPG_GPAR%VTS(JROF)
    ENDDO
  ELSE
    DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS_SURF%GSP_RR%PT_T1(JROF)=YDCPG_GPAR%VTS(JROF)
    ENDDO
  ENDIF
ENDIF
IF(LNUDG)THEN
  CALL CPNUDG ( YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, NFNUDG, YDCPG_DIM%KFLEVG, YDCPG_DIM%KBL,                       &
  & XPNUDG, YDMF_PHYS_SURF%GSD_VF%PNUDM, YDMF_PHYS_SURF%GSP_RR%PT_T1, YDMF_PHYS_SURF%GSP_RR%PW_T1, YDMF_PHYS_SURF%GSP_SB%PQ_T1,  &
  & YDMF_PHYS_SURF%GSP_SG%PF_T1, YDMF_PHYS_NEXT_STATE%T (:, 1:YDCPG_DIM%KFLEVG), YDMF_PHYS_NEXT_STATE%Q (:, 1:YDCPG_DIM%KFLEVG), &
  & YDMF_PHYS_NEXT_STATE%U (:, 1:YDCPG_DIM%KFLEVG), YDMF_PHYS_NEXT_STATE%V (:, 1:YDCPG_DIM%KFLEVG), YDMF_PHYS_NEXT_STATE%SP,     &
  & YDVARS%T%T0, YDVARS%Q%T0, YDVARS%U%T0, YDVARS%V%T0, YDCPG_PHY0%PREHYD(:, YDCPG_DIM%KFLEVG), YDVARS%GEOMETRY%GM%T0,           &
  & YDMF_PHYS_SURF%GSD_VF%PLSM)
ENDIF

CALL MF_PHYS_CVV (YDCPG_DIM, YDVARS%CVV%T0, YDVARS%CVV%T1, YDMODEL)


!*       5.    Final calculations.
!              -------------------

! * Restore the initial value of some pseudo-historical surface buffers
!   if relevant.
IF (LL_SAVE_PHSURF) THEN
  CALL MF_PHYS_SAVE_PHSURF_PART2 (YDCPG_DIM, ZSAV_DDAL, ZSAV_DDOM, ZSAV_ENTCH, ZSAV_FHPS, ZSAV_GZ0F,                    &
  & ZSAV_GZ0HF, ZSAV_HV, ZSAV_PBLH, ZSAV_QSH, ZSAV_UDAL, ZSAV_UDGRO, ZSAV_UDOM, ZSAV_UNEBH, YDMF_PHYS_SURF%GSD_VF%PZ0F, &
  & YDMF_PHYS_SURF%GSD_VH%PPBLH, YDMF_PHYS_SURF%GSD_VH%PQSH, YDMF_PHYS_SURF%GSD_VH%PSPSH, YDMF_PHYS_SURF%GSD_VK%PUDGRO, &
  & YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VV%PZ0H, YDVARS%DAL%T0, YDVARS%DOM%T0, YDVARS%UAL%T0,                 &
  & YDVARS%UEN%T0, YDVARS%UNEBH%T0, YDVARS%UOM%T0, YDSURF, YDMODEL)
ENDIF

! Store horizontal exchange coefficients (3D turbulence) to SL2 buffers

CALL MF_PHYS_BAYRAD (YDCPG_DIM, ZBAY_QRCONV, ZBAY_QSCONV, YDVARS%RCONV%T1, YDVARS%SCONV%T1, YDMODEL, &
& LAROME)

IF (LEDR) THEN
  YDMF_PHYS_SURF%GSD_DI%PXEDR(:,:)=YDMF_PHYS%OUT%EDR(:,:)
ENDIF

CALL MF_PHYS_PRECIPS (YDCPG_DIM, ZPRC_DPRECIPS, ZPRC_DPRECIPS2, YDMF_PHYS_SURF%GSD_XP%PPRECIP, YDMF_PHYS_SURF%GSD_XP2%PPRECIP2, &
& YDMODEL)

END ASSOCIATE
END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('APL_ARPEGE', 1, ZHOOK_HANDLE)

END SUBROUTINE APL_ARPEGE
