SUBROUTINE APL_ARPEGE(YDMF_PHYS_BASE_STATE, YDMF_PHYS_NEXT_STATE, YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, &
& YDCPG_MISC, YDCPG_GPAR, YDCPG_PHY0, YDMF_PHYS, YDCPG_DYN0, YDMF_PHYS_SURF, YDCPG_SL2, YDVARS,       &
& YDCFU, YDMODEL, YDDDH)

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
                              & CPG_SL2_TYPE, CPG_GPAR_TYPE, &
                              & CPG_PHY_TYPE
USE CPG_OPTS_TYPE_MOD   , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD  &
                       , ONLY : MF_PHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE MF_PHYS_BASE_STATE_TYPE_MOD &
                       , ONLY : MF_PHYS_BASE_STATE_TYPE
USE MF_PHYS_NEXT_STATE_TYPE_MOD &
                       , ONLY : MF_PHYS_NEXT_STATE_TYPE
USE TYPE_MODEL         , ONLY : MODEL

USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK    ,DR_HOOK
USE YOMVERT            , ONLY : VP00
USE DDH_MIX            , ONLY : TYP_DDH

USE YOMCFU             , ONLY : TCFU !!! for parameters of FLASH

USE SPP_MOD            , ONLY : YSPP, YSPP_CONFIG

!     -------------------------------------------------------------------------

IMPLICIT NONE


TYPE (MF_PHYS_BASE_STATE_TYPE), INTENT(IN)    :: YDMF_PHYS_BASE_STATE
TYPE (MF_PHYS_NEXT_STATE_TYPE), INTENT(INOUT) :: YDMF_PHYS_NEXT_STATE
TYPE(GEOMETRY),                 INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_BNDS_TYPE),            INTENT(IN)    :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE),            INTENT(IN)    :: YDCPG_OPTS
TYPE(CPG_MISC_TYPE),            INTENT(INOUT) :: YDCPG_MISC
TYPE(CPG_GPAR_TYPE),            INTENT(INOUT) :: YDCPG_GPAR
TYPE(CPG_PHY_TYPE),             INTENT(IN)    :: YDCPG_PHY0
TYPE(MF_PHYS_TYPE),             INTENT(INOUT) :: YDMF_PHYS
TYPE(CPG_DYN_TYPE),             INTENT(IN)    :: YDCPG_DYN0
TYPE(MF_PHYS_SURF_TYPE),        INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(CPG_SL2_TYPE),             INTENT(INOUT) :: YDCPG_SL2
TYPE(FIELD_VARIABLES),          INTENT(INOUT) :: YDVARS
TYPE(TCFU),                     INTENT(IN)    :: YDCFU
TYPE(MODEL),                    INTENT(IN)    :: YDMODEL

 
TYPE(TYP_DDH),                  INTENT(INOUT) :: YDDDH


!     ------------------------------------------------------------------
LOGICAL :: LL_SAVE_PHSURF
INTEGER(KIND=JPIM) :: INSTEP_DEB,INSTEP_FIN
INTEGER(KIND=JPIM) :: JROF, JSPP

REAL(KIND=JPRB) :: ZDIFEXT(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)     ! Extra-GFL fluxes.

REAL(KIND=JPRB) :: ZGP2DSPP(YDCPG_OPTS%KLON,YSPP%N2D)

!     ------------------------------------------------------------------
!        ATTENTION SI KVCLIG < 7 LES CHAMPS SUIVANTS NE SONT
!        PAS REELLEMENT ALLOUES EN MEMOIRE.
!*
!     ------------------------------------------------------------------
!     DECLARATION DES TABLEAUX LOCAUX-GLOBAUX DES PARAMETRISATIONS
INTEGER(KIND=JPIM) :: INLAB_CVPP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZVETAH(0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZVETAF(YDCPG_OPTS%KFLEVG)

!     ------------------------------------------------------------------
!     ARE DIMENSIONNED 0:KLEV ONLY IN ORDER TO KEEP IN MIND
!     THAT THEY ARE COMPUTED AT "HALF LEVELS".
!     THEY ARE USED HOWEVER FROM 1 TO KLEV.

REAL(KIND=JPRB) :: ZXTROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZXUROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZMRIPP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)

! ZMRIMC     : M(C) from P. Marquet's moist Ri computation - for TKE correction after TOMs
! ZMRICTERM  : Rv/R.F(C)-1/M(C).T/Tv from P. Marquet's moist Ri computation - for TKE correction after TOMs

!DISSIPATION TIME SCALE TAU  -FOR TOM's CALCULATION
REAL(KIND=JPRB) :: ZKTROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZKUROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZNBVNO(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZKQROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZKQLROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZNEBS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZQLIS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZNEBS0(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZQLIS0(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZNEBC0(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)   !Nebulosite convective radiative
REAL(KIND=JPRB) :: ZNEBDIFF(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) !Nebulosite: calcul de la diffusion
REAL(KIND=JPRB) :: ZNEBCH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)   !Nebulosite convective condensation
REAL(KIND=JPRB) :: ZUNEBH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)   !Nebulosite convective histo
REAL(KIND=JPRB) :: ZFPCOR(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTW(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! Temporar storage for updated PTW)
REAL(KIND=JPRB) :: ZPOID(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)  ! DP/(YDMODEL%YRCST%RG*DT) FOR A GIVEN LEVEL AND A GIVEN TIME STEP.
REAL(KIND=JPRB) :: ZQV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! corrected (for negative values) vapour
REAL(KIND=JPRB) :: ZQI(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! corrected (for negative values) cloud ice
REAL(KIND=JPRB) :: ZQL(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! corrected (for negative values) cloud liquid
REAL(KIND=JPRB) :: ZQR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! corrected (for negative values) rain
REAL(KIND=JPRB) :: ZQS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! corrected (for negative values) snow
REAL(KIND=JPRB) :: ZTENHA(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENQVA(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZCP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)    ! new cp for turbulent diffusion
REAL(KIND=JPRB) :: ZFPLSL(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! total liquid water flux: diff+sedi+rain
REAL(KIND=JPRB) :: ZFPLSN(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! total solid water flux: diff+sedi+snow
REAL(KIND=JPRB) :: ZSEDIQL(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! sedimentation flux of cloud liquid water
REAL(KIND=JPRB) :: ZSEDIQI(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! sedimentation flux of cloud ice water
REAL(KIND=JPRB) :: ZDIFCVPPQ(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de CVPP (KFB or EDKF) sur Qv
REAL(KIND=JPRB) :: ZDIFCVPPS(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de CVPP (KFB or EDKF) sur CpT
REAL(KIND=JPRB) :: ZDIFCVPPU(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de CVPP (EDKF) sur U
REAL(KIND=JPRB) :: ZDIFCVPPV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de CVPP (EDKF) sur V
REAL(KIND=JPRB) :: ZEDMFQ(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Mass flux part of EDMF flux for Qv
REAL(KIND=JPRB) :: ZEDMFS(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Mass flux part of EDMF flux for CpT
REAL(KIND=JPRB) :: ZEDMFU(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Mass flux part of EDMF flux for U
REAL(KIND=JPRB) :: ZEDMFV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Mass flux part of EDMF flux for V
REAL(KIND=JPRB) :: ZMF_UP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Mass flux for implicit formulation of EDMF equation (LEDMFI)
REAL(KIND=JPRB) :: ZCONDCVPPL(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de condensation liquide du a CVVPP (KFB)
REAL(KIND=JPRB) :: ZCONDCVPPI(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de condensation glace du a CVVPP (KFB)
REAL(KIND=JPRB) :: ZPRODTH_CVPP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) ! Flux de production thermique de TKE du a CVPP(KFB)

!!for BAYRAD
REAL(KIND=JPRB) :: ZDE2MR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG) ! temporary array for conversion of density to mixing ratio
REAL(KIND=JPRB) :: ZCOEFRAIN(2) ! RTTOVSCATT coefficients to convert flux to density for rain
REAL(KIND=JPRB) :: ZCOEFSNOW(2) ! RTTOVSCATT coefficients to convert flux to density for snow
REAL(KIND=JPRB) :: ZRHORAIN ! RTTOVSCATT fixed density for rain
REAL(KIND=JPRB) :: ZRHOSNOW ! RTTOVSCATT fixed density of snow
!-----------------------------------------------------------------

REAL(KIND=JPRB) :: ZXDROV(YDCPG_OPTS%KLON),ZXHROV(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZUGST(YDCPG_OPTS%KLON),ZVGST(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZCDROV(YDCPG_OPTS%KLON),ZCHROV(YDCPG_OPTS%KLON),ZDQSTS(YDCPG_OPTS%KLON),ZGWDCS(YDCPG_OPTS%KLON),&
 & ZHQ(YDCPG_OPTS%KLON),ZHU(YDCPG_OPTS%KLON),ZHTR(YDCPG_OPTS%KLON),&
 & ZRTI(YDCPG_OPTS%KLON),ZDPHI(YDCPG_OPTS%KLON),ZPRS(YDCPG_OPTS%KLON),ZSTAB(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZWFC(YDCPG_OPTS%KLON),ZWPMX(YDCPG_OPTS%KLON),ZWLMX(YDCPG_OPTS%KLON),ZWSEQ(YDCPG_OPTS%KLON),&
 & ZWSMX(YDCPG_OPTS%KLON),ZWWILT(YDCPG_OPTS%KLON),&
 & ZC3(YDCPG_OPTS%KLON),ZCG(YDCPG_OPTS%KLON),ZCN(YDCPG_OPTS%KLON),&
 & ZNEIJG(YDCPG_OPTS%KLON),ZNEIJV(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZPCLS(YDCPG_OPTS%KLON)


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
REAL(KIND=JPRB) :: ZHUC(YDCPG_OPTS%KFLEVG), ZBLH(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZQO3(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZAER(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,6)
REAL(KIND=JPRB) :: ZAERINDS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZDPHIV(YDCPG_OPTS%KLON),ZDPHIT(YDCPG_OPTS%KLON)

! - (PROFILS DIAGNOSTIQUES)
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

REAL(KIND=JPRB) :: ZTRSOD(YDCPG_OPTS%KLON)

!            2-D ARRAYS
!            ----------

!* OUTPUT ARGUMENTS FOR THE ECMWF PHYSICS

!            0.2  LOCAL ARRAYS FOR ECMWF PHYSICS PACKAGE
!                 --------------------------------------

REAL(KIND=JPRB) :: ZCEMTR(YDCPG_OPTS%KLON,0:1) , ZCTRSO(YDCPG_OPTS%KLON,0:1)
REAL(KIND=JPRB) :: ZALBD(YDCPG_OPTS%KLON,YDCPG_OPTS%KSW), ZALBP(YDCPG_OPTS%KLON,YDCPG_OPTS%KSW)
REAL(KIND=JPRB) :: ZSFSWDIR (YDCPG_OPTS%KLON,YDCPG_OPTS%KSW), ZSFSWDIF (YDCPG_OPTS%KLON,YDCPG_OPTS%KSW)
REAL(KIND=JPRB) :: ZTRSODIR (YDCPG_OPTS%KLON,YDCPG_OPTS%KSW), ZTRSODIF (YDCPG_OPTS%KLON,YDCPG_OPTS%KSW)

!            1-D ARRAYS
!            ----------

REAL(KIND=JPRB) :: ZSUDU(YDCPG_OPTS%KLON) 
REAL(KIND=JPRB) :: ZTHETAVS(YDCPG_OPTS%KLON)

!            2-D ARRAYS
!            ----------

REAL(KIND=JPRB) :: ZGEOSLC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTHETAV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

!            LOCAL ARRAYS FOR TKE
! ZCOEFN : COEFFICIENT STATISTIQUE POUR LES FLUX D'EAUX CONDENSEES.

REAL(KIND=JPRB) :: ZCOEFN(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

!            LOCAL ARRAYS FOR ACVPPKF
REAL(KIND=JPRB) :: ZQLI_CVPP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZNEB_CVPP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

INTEGER(KIND=JPIM) :: JCHA, JLEV, JLON, JSG
 ! useful size of klon arrays for mesonh physics
LOGICAL :: LLCLS, LLHMT, LLREDPR

REAL(KIND=JPRB) :: ZAEN, ZEPS, ZEPS0, ZEPSNEB

! Implicit coupling coefficients
REAL(KIND=JPRB) :: ZCFBTH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZDTMSE,ZRHGMT,ZSTATI
REAL(KIND=JPRB) :: &
 & ZCFATH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZCFAU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),&
 & ZCFBU(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZCFBV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZSRAIN(YDCPG_OPTS%KLON), ZSSNOW(YDCPG_OPTS%KLON), ZSGROUPEL(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZTSN(YDCPG_OPTS%KLON)
! FOR Hv
! FOR DUST

 ! SCOND MEMBRE POUR LES SCALAIRES PASSIFS
REAL(KIND=JPRB) :: ZINVG

! TRAITEMENT DES SCALAIRES PASSIFS

REAL(KIND=JPRB) :: ZDZZ(YDCPG_OPTS%KLON,1,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZZZ (YDCPG_OPTS%KLON,1,YDCPG_OPTS%KFLEVG)

!         ACFLUSO (ECUME) local variable
!-------------------------------------------
REAL(KIND=JPRB) :: ZCE(YDCPG_OPTS%KLON), ZCEROV(YDCPG_OPTS%KLON), ZCRTI(YDCPG_OPTS%KLON)

!        New ACDIFV1 local variable
!--------------------------------------------
REAL(KIND=JPRB)   :: ZXURO(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)   :: ZXQRO(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)   :: ZXTRO(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)

!        New ACNORGWD local variables
!--------------------------------------------

REAL(KIND=JPRB) :: ZFLX_LOTT_GWU(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG), ZFLX_LOTT_GWV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)


!    TKE+ for ACCLDIA
REAL(KIND=JPRB)   :: ZTKE1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

!   For ACVISIH
REAL(KIND=JPRB)   :: ZQGM(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)


!        New ARP_GROUND_PARAM local variable
!------------------------------------------------

REAL(KIND=JPRB)   :: ZALPHA1(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)   :: ZCOEFA (YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)   :: ZLVT   (YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)   :: ZQICE  (YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)   :: ZDIFWQ (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)   :: ZDIFWS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)   :: ZSC_FEVI (YDCPG_OPTS%KLON),ZSC_FEVN(YDCPG_OPTS%KLON),ZSC_FCLL(YDCPG_OPTS%KLON),ZSC_FCLN(YDCPG_OPTS%KLON)

!           TRAJECTORY (For diffusion !) local VARIABLES
!           ----------------------------
REAL(KIND=JPRB) :: ZTRAJGWD(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG) !Traj buffer saved for TL/AD (YDMODEL%YRML_PHY_MF%YRSIMPHL%LGWDSPNL)

REAL(KIND=JPRB)    :: ZRVMD
LOGICAL            :: LLAERO
REAL(KIND=JPRB)    :: ZAIPCMT(YDCPG_OPTS%KLON) ! Activity Index of PCMT: 1. if PCMT is active, 0. else case.


REAL(KIND=JPRB)    :: ZQIC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG),ZQLC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)    :: ZQLI_CVP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)    :: ZQC_DET_PCMT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)    :: ZFPLS(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG),ZFPLC(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZDCAPE(YDCPG_OPTS%KLON) ! Descending CAPE for gusts.

! ACRANEB/ACRANEB2 local variables
! --------------------------------
REAL(KIND=JPRB) :: ZDECRD   (YDCPG_OPTS%KLON) ! decorrelation depth for cloud overlaps

! Precipitation type diagnostics
!--------------------------------
REAL(KIND=JPRB)   :: ZFPLSG(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
INTEGER (KIND=JPIM)  :: IMOC_CLPH (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZBAY_QRCONV (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZBAY_QSCONV (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZDSA_C1 (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZDSA_C2 (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZDSA_CPS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZDSA_LHS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZDSA_RS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZFLU_CDN (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZFLU_CD (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZFLU_CH (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZFLU_EMIS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZFLU_FEVI (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KTSSG+1)
REAL(KIND=JPRB)     :: ZFLU_NEIJ (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZFLU_QSATS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZFLU_QSAT (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZFLU_VEG (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZMSC_FRMQ (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZMSC_LH (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZMSC_LSCPE (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZMSC_QW (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZMSC_TW (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZPFL_FEFB1 (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZPFL_FEFB2 (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZPFL_FEFB3 (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZPFL_FPLCH (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZPFL_FPLSH (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZPFL_FTKEI (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZPFL_FTKE (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZPRC_DPRECIPS2 (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%NDTPREC2)
REAL(KIND=JPRB)     :: ZPRC_DPRECIPS (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%NDTPREC)
REAL(KIND=JPRB)     :: ZRDG_CVGQ (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZRDG_LCVQ (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZRDG_MU0LU (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZRDG_MU0M (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZRDG_MU0N (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZRDG_MU0 (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZSAV_DDAL (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZSAV_DDOM (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZSAV_ENTCH (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZSAV_FHPS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZSAV_GZ0F (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZSAV_GZ0HF (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZSAV_HV (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZSAV_PBLH (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZSAV_QSH  (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZSAV_UDAL (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZSAV_UDGRO (YDCPG_OPTS%KLON)
REAL(KIND=JPRB)     :: ZSAV_UDOM (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB)     :: ZSAV_UNEBH (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

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
#include "checkmv.intfb.h"
#include "cpphinp.intfb.h"
#include "cpqsol.intfb.h"
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
#include "qngcor.intfb.h"
#include "aplpar_flexdia.intfb.h"
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
#include "apl_arpege_surface_update.intfb.h"
#include "apl_arpege_atmosphere_update.intfb.h"
#include "apl_arpege_init.intfb.h"
#include "apl_arpege_init_surfex.intfb.h"

IF (LHOOK) CALL DR_HOOK('APL_ARPEGE', 0, ZHOOK_HANDLE)

ASSOCIATE (YDCST => YDMODEL%YRCST)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, YDTOPH=>YDMODEL%YRML_PHY_MF%YRTOPH,       &
& YDRIP=>YDMODEL%YRML_GCONF%YRRIP, YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY, YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY,    &
& YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, YGFL=>YDMODEL%YRML_GCONF%YGFL, YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR,     &
& YDPRECIPS=>YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS, YDSTA=>YDGEOMETRY%YRSTA, YDMCC=>YDMODEL%YRML_AOC%YRMCC,     &
& YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD, YDPHY3=>YDMODEL%YRML_PHY_MF%YRPHY3, YDPHY1=>YDMODEL%YRML_PHY_MF%YRPHY1, &
& YDPHY0=>YDMODEL%YRML_PHY_MF%YRPHY0)


ASSOCIATE(TSPHY=>YDPHY2%TSPHY, NTSSG=>YDDPHY%NTSSG, LMSE=>YDARPHY%LMSE, LNEBN=>YDPHY%LNEBN, NDPSFI=>YDPHY%NDPSFI,    &
& LRRGUST=>YDPHY%LRRGUST, LEDR=>YDPHY%LEDR, NTPLUI=>YDTOPH%NTPLUI, XMINLM=>YDPHY0%XMINLM, HUTIL2=>YDPHY0%HUTIL2,     &
& HUTIL1=>YDPHY0%HUTIL1, XMAXLM=>YDPHY0%XMAXLM, HUCOE=>YDPHY0%HUCOE, HUTIL=>YDPHY0%HUTIL, NPCLO1=>YDPHY0%NPCLO1,     &
& NPCLO2=>YDPHY0%NPCLO2, HSOLIWR=>YDPHY1%HSOLIWR, WSMX=>YDPHY1%WSMX, HSOLIT0=>YDPHY1%HSOLIT0, HSOL=>YDPHY1%HSOL,     &
& WPMX=>YDPHY1%WPMX, LRAFTKE=>YDPHY2%LRAFTKE, HVCLS=>YDPHY2%HVCLS, HTCLS=>YDPHY2%HTCLS, RII0=>YDPHY3%RII0,           &
& YA=>YGFL%YA, YIRAD=>YGFL%YIRAD, YLRAD=>YGFL%YLRAD, LHUCN=>YDPHY%LHUCN, LSTRAS=>YDPHY%LSTRAS, LSOLV=>YDPHY%LSOLV,   &
& NTDRME=>YDTOPH%NTDRME, NTDRAG=>YDTOPH%NTDRAG, RMESOQ=>YDTOPH%RMESOQ, RMESOT=>YDTOPH%RMESOT, RMESOU=>YDTOPH%RMESOU, &
& NTQSAT=>YDTOPH%NTQSAT, LMCC03=>YDMCC%LMCC03, RHGMT=>YDRIP%RHGMT, RSTATI=>YDRIP%RSTATI, STPRE=>YDSTA%STPRE,         &
& STPREH=>YDSTA%STPREH, STTEM=>YDSTA%STTEM, LGCHECKMV=>YDPHY%LGCHECKMV             )

ASSOCIATE(PAPRSF=> YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, PAPRS => YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,                                   &
& PAPHIF=> YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, PAPHI=> YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, PDELP => YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
& PR=>YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, PT=> YDMF_PHYS_BASE_STATE%T, PU=> YDMF_PHYS_BASE_STATE%U,                                       &
& PV=>YDMF_PHYS_BASE_STATE%V, PCP=>YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP)


INSTEP_DEB=1
INSTEP_FIN=1

! SPP 
IF ( YSPP_CONFIG%LSPP ) THEN
 DO JSPP=1,YSPP%N2D
   ZGP2DSPP(:,JSPP) = YSPP%GP_ARP(JSPP)%GP2D(:,1,YDCPG_BNDS%KBL)
 ENDDO
ENDIF

CALL CPPHINP(YDCPG_OPTS%LVERTFE, YDGEOMETRY, YDMODEL, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDVARS%GEOMETRY%GEMU%T0, YDVARS%GEOMETRY%GELAM%T0, &
& YDVARS%U%T0, YDVARS%V%T0, YDVARS%Q%T0, YDVARS%Q%DL, YDVARS%Q%DM, YDVARS%CVGQ%DL, YDVARS%CVGQ%DM, YDCPG_PHY0%XYB%RDELP,                     &
& YDCPG_DYN0%CTY%EVEL, YDVARS%CVGQ%T0, ZRDG_MU0, ZRDG_MU0LU, ZRDG_MU0M, ZRDG_MU0N, ZRDG_CVGQ)
ZRDG_LCVQ(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)=ZRDG_CVGQ(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)

DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
  ZFLU_QSATS(JROF)=0.0_JPRB
ENDDO

CALL MF_PHYS_FPL_PART1 (YDCPG_BNDS, YDCPG_OPTS, ZPFL_FPLCH, ZPFL_FPLSH, YDVARS%CPF%T0, YDVARS%SPF%T0, &
& YDMODEL)


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

LL_SAVE_PHSURF=YDCPG_OPTS%LCONFX
IF (LL_SAVE_PHSURF) THEN
  CALL MF_PHYS_SAVE_PHSURF_PART1 (YDCPG_BNDS, YDCPG_OPTS, ZSAV_DDAL, ZSAV_DDOM, ZSAV_ENTCH,                           &
  & ZSAV_FHPS, ZSAV_GZ0F, ZSAV_GZ0HF, ZSAV_HV, ZSAV_PBLH, ZSAV_QSH, ZSAV_UDAL, ZSAV_UDGRO, ZSAV_UDOM,                 &
  & ZSAV_UNEBH, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VH%PPBLH, YDMF_PHYS_SURF%GSD_VH%PQSH,                  &
  & YDMF_PHYS_SURF%GSD_VH%PSPSH, YDMF_PHYS_SURF%GSD_VK%PUDGRO, YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VV%PZ0H, &
  & YDVARS%DAL%T0, YDVARS%DOM%T0, YDVARS%UAL%T0, YDVARS%UEN%T0, YDVARS%UNEBH%T0, YDVARS%UOM%T0,                       &
  & YDMODEL)
ENDIF




CALL APLPAR_INIT (YDCPG_OPTS%LAROME, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, NTSSG,       &
& YDCPG_OPTS%YRSURF_DIMS%YSP_SBD%NLEVS, YDMF_PHYS_SURF%GSD_VF%PVEG, ZMSC_FRMQ, ZDSA_CPS, ZDSA_LHS, ZDSA_RS, ZMSC_LH,      &
& ZMSC_LSCPE, ZFLU_QSAT, ZMSC_QW, ZMSC_TW, ZFLU_CD, ZFLU_CDN, ZFLU_CH, ZDSA_C1, ZDSA_C2, ZFLU_EMIS, ZFLU_FEVI, ZPFL_FTKE, &
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

DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
  YDMF_PHYS_SURF%GSD_VF%PLSM(JROF)=REAL(NINT(YDMF_PHYS_SURF%GSD_VF%PLSM(JROF)),JPRB)
ENDDO

!
!-------------------------------------------------
! Check magnitude of model variables.
!-------------------------------------------------
!
IF(LGCHECKMV) CALL CHECKMV(YDCPG_OPTS%NINDAT, YDMODEL%YRCST, YDRIP, YDPHY0, YDPHY2, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA,       &
              & YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDCPG_OPTS%KSTEP, PAPHI, PAPHIF, PAPRS, PAPRSF, YDVARS%GEOMETRY%GELAM%T0, &
              & YDVARS%GEOMETRY%GEMU%T0, ZRDG_MU0, YDMF_PHYS_SURF%GSD_VF%PLSM, PT, YDMF_PHYS_BASE_STATE%Q,                    &
              & YDMF_PHYS_BASE_STATE%YGSP_RR%T                                                                                &
              &                                                                       )
!     ------------------------------------------------------------------

LLREDPR=.FALSE.
ZRVMD=YDMODEL%YRCST%RV-YDMODEL%YRCST%RD
! SURFEX  and passive scalar
IF (YDCPG_OPTS%LCONFX) THEN
  ZDTMSE=0.01_JPRB
  ZSTATI=RSTATI-ZDTMSE/2._JPRB
ELSE
  ZDTMSE=TSPHY
  ZSTATI=RSTATI
ENDIF
ZRHGMT=REAL(RHGMT,JPRB)


! INITIALISATION DE LA COORDONNEE ETA.

ZVETAH(YDCPG_OPTS%KTDIA-1)=STPREH(YDCPG_OPTS%KTDIA-1)/VP00
DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
  ZVETAH(JLEV)=STPREH(JLEV)/VP00
  ZVETAF(JLEV)=STPRE (JLEV)/VP00
ENDDO

!     HUMIDITE CRITIQUE
ZEPS=1.E-12_JPRB
IF(LHUCN) THEN
  DO JLEV = YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG
    ZHUC(JLEV)=1.0_JPRB-MAX( HUCOE*ZVETAF(JLEV)*(1.0_JPRB-ZVETAF(JLEV))/((&
    & 1.0_JPRB+HUTIL1*(ZVETAF(JLEV)-0.5_JPRB))*(&
    & 1.0_JPRB+HUTIL2*(ZVETAF(JLEV)-0.5_JPRB))),ZEPS)
  ENDDO
ELSE
  DO JLEV = YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG
    ZHUC(JLEV)=1.0_JPRB-MAX( HUCOE*ZVETAF(JLEV)**NPCLO1*(&
     & 1.0_JPRB-ZVETAF(JLEV))**NPCLO2*(&
     & 1.0_JPRB+SQRT(HUTIL)*(ZVETAF(JLEV)-0.5_JPRB)),ZEPS)
  ENDDO
ENDIF

CALL APL_ARPEGE_INIT (YDMODEL%YRCST, YDMF_PHYS_BASE_STATE, YDCPG_BNDS, YDCPG_OPTS, YDCPG_MISC, YDMF_PHYS, &
& YDMF_PHYS_SURF, YDVARS, YDMODEL, IMOC_CLPH, INLAB_CVPP, ZAER, ZAERINDS, ZAIPCMT, ZALBD, ZALBP,          &
& ZALPHA1, ZCEMTR, ZCFATH, ZCFAU, ZCFBTH, ZCFBU, ZCFBV, ZCOEFA, ZCONDCVPPI, ZCONDCVPPL, ZCTRSO,           &
& ZDECRD, ZDIFCVPPQ, ZDIFCVPPS, ZDIFCVPPU, ZDIFCVPPV, ZDIFWQ, ZDIFWS, ZEDMFQ, ZEDMFS, ZEDMFU, ZEDMFV,     &
& ZEPS0, ZEPSNEB, ZFPCOR, ZKQLROV, ZKQROV, ZKTROV, ZKUROV, ZLVT, ZMF_UP, ZMRIPP, ZNEBC0, ZNEBCH,          &
& ZNEBDIFF, ZNEBS, ZNEBS0, ZNEB_CVPP, ZPFL_FPLCH, ZPFL_FPLSH, ZPOID, ZPRODTH_CVPP, ZQC_DET_PCMT, ZQI,     &
& ZQIC, ZQICE, ZQL, ZQLC, ZQLIS, ZQLIS0, ZQLI_CVP, ZQLI_CVPP, ZQO3, ZQR, ZQS, ZQV, ZSC_FCLL,              &
& ZSC_FCLN, ZSC_FEVI, ZSC_FEVN, ZSEDIQI, ZSEDIQL, ZSFSWDIF, ZSFSWDIR, ZSUDU, ZTENHA, ZTENQVA, ZTENT,      &
& ZTRSOD, ZTRSODIF, ZTRSODIR, ZUNEBH, ZXDROV, ZXHROV, ZXQRO, ZXTRO, ZXTROV, ZXURO, ZXUROV)


!*
!     ------------------------------------------------------------------
!     4.- CALCULS THERMODYNAMIQUES
!     ----------------------------
  
CALL ACTQSAT (YDMODEL%YRCST, YDPHY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTQSAT, YDCPG_OPTS%KFLEVG, &
& PAPRSF, PCP, ZQV, PT, ZGEOSLC, ZMSC_LH, ZMSC_LSCPE, ZFLU_QSAT, ZMSC_QW, YDCPG_MISC%RH, ZMSC_TW)
  

!*
!     ------------------------------------------------------------------
!     4.BIS. COEFFICIENTS THERMO-HYDRIQUES DU SOL
!     -------------------------------------------

IF ( .NOT.LMSE ) THEN
  IF ( LSOLV ) THEN
    LLHMT=.FALSE.
    CALL ACSOL (YDCPG_OPTS%YRCLI, YDMODEL%YRCST, YDPHY, YDPHY1, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDMF_PHYS_SURF%GSD_VV%PARG, &
    & YDMF_PHYS_SURF%GSD_VV%PD2, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H, YDMF_PHYS_SURF%GSD_VF%PZ0RLF,                           &
    & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PIVEG, YDMF_PHYS_SURF%GSD_VV%PLAI, YDMF_PHYS_BASE_STATE%YGSP_SG%A,                       &
    & YDMF_PHYS_BASE_STATE%YGSP_SG%R, YDMF_PHYS_SURF%GSD_VV%PSAB, YDMF_PHYS_BASE_STATE%YGSP_SG%F, YDMF_PHYS_BASE_STATE%YGSP_RR%T,                &
    & YDMF_PHYS_SURF%GSD_VF%PVEG, YDMF_PHYS_BASE_STATE%YGSP_SB%Q, YDMF_PHYS_BASE_STATE%YGSP_SB%TL, YDMF_PHYS_BASE_STATE%YGSP_RR%W,               &
    & YDMF_PHYS_BASE_STATE%YGSP_RR%IC, LLHMT, ZDSA_C1, ZDSA_C2, ZC3, ZCG, ZCN, YDMF_PHYS%OUT%CT,                                                 &
    & ZNEIJG, ZNEIJV, ZWFC, ZWPMX, ZWSEQ, ZWSMX, ZWWILT)
  ELSE

!            INITIALISATION DE L'INERTIE THERMIQUE DU SOL.

!DEC$ IVDEP
    DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      YDMF_PHYS%OUT%CT(JLON)=HSOL /(&
      & 1.0_JPRB+HSOLIWR*(YDMF_PHYS_BASE_STATE%YGSP_RR%W(JLON)+YDMF_PHYS_BASE_STATE%YGSP_SB%Q(JLON,1))/(WSMX+WPMX)&
      & *EXP(-0.5_JPRB*(HSOLIT0*(YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)-YDMODEL%YRCST%RTT))**2))
    ENDDO
  ENDIF
ENDIF


CALL APL_ARPEGE_INIT_SURFEX (YDMF_PHYS_BASE_STATE, YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDCPG_MISC, &
& YDCPG_GPAR, YDMF_PHYS, YDVARS, YDMODEL, ZALBD, ZALBP, ZFLU_EMIS, ZSGROUPEL, ZSRAIN, ZSSNOW, ZTSN)

!*
!     ------------------------------------------------------------------
!     5.- STRUCTURE ET CHAMPS DANS LA COUCHE LIMITE DE SURFACE
!     ------------------------------------------------------------------

!        INITIALISATION DES HAUTEURS "METEO".

  
DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
  ZDPHIV(JLON)=YDMODEL%YRCST%RG*HVCLS
  ZDPHIT(JLON)=YDMODEL%YRCST%RG*HTCLS
ENDDO
  

  
LLCLS=.TRUE.
LLHMT=.TRUE.
IF (LMSE) THEN      
  CALL ACHMTLS (YDMODEL%YRCST, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, &
  & PAPHI, PAPHIF, PAPRS, PAPRSF, PR, PT, PU, PV, YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDCPG_MISC%QS,                            &
  & ZDPHIT, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, LLCLS, ZNBVNO, ZMRIPP, ZDSA_CPS, ZGWDCS, ZDSA_LHS,                       &
  & ZPCLS, ZFLU_CD, ZFLU_CDN)
!       Computation of ZRTI
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZDPHI(JLON)=PAPHIF(JLON,YDCPG_OPTS%KFLEVG)-PAPHI(JLON,YDCPG_OPTS%KFLEVG)
    ZPRS(JLON)=YDMODEL%YRCST%RD+ZRVMD*YDCPG_MISC%QS(JLON)
    ZRTI(JLON)=2.0_JPRB/(PR(JLON,YDCPG_OPTS%KFLEVG)*PT(JLON,YDCPG_OPTS%KFLEVG)+YDMODEL%YRCST%RKAPPA*ZDPHI(JLON)&
    & +ZPRS(JLON)*YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON))
  ENDDO      
ELSE      
  CALL ACHMT (YDCPG_OPTS%YRCLI, YDMODEL%YRML_PHY_MF%YRPHY, YDMODEL%YRML_PHY_MF%YRPHY0, YDMODEL%YRML_PHY_MF%YRPHY1, YDMODEL%YRML_PHY_MF%YRPHY2,             &
  & YDMODEL%YRCST, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDCPG_OPTS%YRSURF_DIMS%YSD_VVD%NUMFLDS>=8.AND.LSOLV,            &
  & PAPHI, PAPHIF, PAPRS, PAPRSF, PCP, ZQV, PR, PT, PU, PV, ZPFL_FPLSH, ZPFL_FPLCH, ZDPHIT, ZDPHIV, YDMF_PHYS_SURF%GSD_VF%PZ0F,                            &
  & YDMF_PHYS_SURF%GSD_VV%PZ0H, YDMF_PHYS_SURF%GSD_VF%PZ0RLF, YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VF%PLSM,                                       &
  & ZNEIJG, ZNEIJV, YDMF_PHYS_BASE_STATE%YGSP_SG%F, YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS_SURF%GSD_VF%PVEG,                                            &
  & ZWFC, YDMF_PHYS_BASE_STATE%YGSP_RR%W, YDMF_PHYS_BASE_STATE%YGSP_RR%IC, LLCLS, LLHMT, ZNBVNO,                                                           &
  & ZMRIPP, ZFLU_CD, ZFLU_CDN, ZCDROV, ZFLU_CH, ZCHROV, ZDSA_CPS, ZDQSTS, ZGWDCS, YDMF_PHYS%OUT%GZ0,                                                       &
  & YDMF_PHYS%OUT%GZ0H, ZHQ, ZHU, ZFLU_NEIJ, YDMF_PHYS%OUT%QCLS, YDCPG_MISC%QS, ZFLU_QSATS, YDMF_PHYS%OUT%RHCLS,                                           &
  & ZDSA_RS, ZRTI, ZSTAB, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS%OUT%NUCLS,                                                 &
  & YDMF_PHYS%OUT%NVCLS, ZPCLS, ZFLU_VEG, ZXDROV, ZXHROV, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST)
ENDIF
    
ZCP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG) = PCP(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,1:YDCPG_OPTS%KFLEVG)

    ! COMPUTATION OF 'DRY' mixing lengths : lm_d lh_d
    ! COMPUTATION OF ZPBLH - PBL HEIGHT

    
DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZTHETAV(JLON,JLEV)=PT(JLON,JLEV)*(1.0_JPRB+YDMODEL%YRCST%RETV*ZQV(JLON,JLEV))&
    & *(YDMODEL%YRCST%RATM/PAPRSF(JLON,JLEV))**YDMODEL%YRCST%RKAPPA  
  ENDDO
ENDDO
DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
  ZTHETAVS(JLON)=YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)*(1.0_JPRB+YDMODEL%YRCST%RETV*YDCPG_MISC%QS(JLON))&
  & *(YDMODEL%YRCST%RATM/PAPRS(JLON,YDCPG_OPTS%KFLEVG))**YDMODEL%YRCST%RKAPPA  
ENDDO

CALL ACCLPH (YDMODEL%YRCST, YDPHY0, YDPHY2, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA,&
& YDCPG_OPTS%KFLEVG, ZTHETAV, PAPHI, PAPHIF, PU, PV, ZTHETAVS, IMOC_CLPH, YDMF_PHYS%OUT%CLPH, YDMF_PHYS%OUT%VEIN, &
& ZUGST, ZVGST)

ZBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=YDMF_PHYS%OUT%CLPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)

CALL APL_ARPEGE_OCEANIC_FLUXES (YDMF_PHYS_BASE_STATE, YDCPG_BNDS, YDCPG_OPTS, YDMF_PHYS,              &
& YDMF_PHYS_SURF, YDMODEL, LLHMT, ZCDROV, ZCEROV, ZCHROV, ZDPHIT, ZDPHIV, ZDSA_RS, ZFLU_CD, ZFLU_CDN, &
& ZFLU_CH, ZFLU_QSATS)


CALL APL_WIND_GUST (YDMF_PHYS_BASE_STATE, YDCPG_BNDS, YDCPG_OPTS, YDMF_PHYS, YDVARS, YDMODEL, &
& IMOC_CLPH, ZBLH, ZDCAPE)


CALL APL_ARPEGE_SHALLOW_CONVECTION_AND_TURBULENCE (YDMF_PHYS_BASE_STATE, YDCPG_BNDS, YDCPG_OPTS,      &
& YDCPG_MISC, YDMF_PHYS, YDCPG_DYN0, YDMODEL, YDDDH, INLAB_CVPP, ZCDROV, ZCHROV, ZCOEFN, ZCONDCVPPI,  &
& ZCONDCVPPL, ZDIFCVPPQ, ZDIFCVPPS, ZEPS, ZFLU_CD, ZFLU_CH, ZKQLROV, ZKQROV, ZKTROV, ZKUROV,          &
& ZMSC_LSCPE, ZNBVNO, ZNEBS, ZNEBS0, ZNEB_CVPP, ZPFL_FPLCH, ZPFL_FTKE, ZPFL_FTKEI, ZPRODTH_CVPP, ZQI, &
& ZQIC, ZQL, ZQLC, ZQLIS, ZQLIS0, ZQLI_CVPP, ZQV, ZTKE1, ZXTROV, ZXUROV)

!**

!      7.1 Albedo et emissivite en presence de neige
!          Albedo and emissivity with snow

CALL APL_ARPEGE_ALBEDO_COMPUTATION (YDMF_PHYS_BASE_STATE, YDCPG_BNDS, YDCPG_OPTS, YDMF_PHYS, &
& YDMF_PHYS_SURF, YDMODEL, ZALBD, ZALBP, ZEPS0, ZFLU_EMIS, ZFLU_NEIJ, ZRDG_MU0)

CALL APL_ARPEGE_AEROSOLS_FOR_RADIATION (YDMF_PHYS_BASE_STATE, YDCPG_BNDS, YDCPG_OPTS, YDMF_PHYS_SURF, &
& YDMODEL, ZAER, ZAERINDS)

CALL APL_ARPEGE_CLOUDINESS (YDMF_PHYS_BASE_STATE, YDCPG_BNDS, YDCPG_OPTS, YDCPG_MISC, YDMF_PHYS,    &
& YDVARS, YDMODEL, LLREDPR, ZAIPCMT, ZBLH, ZDECRD, ZFLU_QSAT, ZHUC, ZMSC_QW, ZNEBC0, ZNEBCH, ZNEBS, &
& ZNEBS0, ZNEB_CVPP, ZPFL_FPLCH, ZQI, ZQL, ZQLIS, ZQLIS0, ZQLI_CVP, ZQLI_CVPP, ZQV, ZUNEBH, ZVETAF)
  
CALL APL_ARPEGE_RADIATION (YDMF_PHYS_BASE_STATE, YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDCPG_MISC, &
& YDCPG_GPAR, YDMF_PHYS, YDMF_PHYS_SURF, YDVARS, YDMODEL, ZAER, ZAERINDS, ZALBD, ZALBP, ZCEMTR,  &
& ZCTRSO, ZFLU_EMIS, ZFLU_QSAT, ZQO3, ZQR, ZQS, ZQV, ZRDG_MU0, ZRDG_MU0LU, ZRDG_MU0M, ZSFSWDIF,  &
& ZSFSWDIR, ZSUDU, ZTENT, ZTRSOD, ZTRSODIF, ZTRSODIR)

!*
!     ------------------------------------------------------------------

!     7.BIS. BILAN HYDRIQUE DU SOL
!     ----------------------------
!     CALCUL DES RESISTANCES A L'EVAPOTRANSPIRATION HV ET
!     A LA TRANSPIRATION
!     ------------------------------------------------------------------
!     HTR DU COUVERT VEGETAL
!     ----------------------

CALL APL_ARPEGE_SOIL_HYDRO (YDMF_PHYS_BASE_STATE, YDCPG_BNDS, YDCPG_OPTS, YDMF_PHYS, YDMF_PHYS_SURF, &
& YDMODEL, ZCHROV, ZFLU_NEIJ, ZFLU_QSAT, ZFLU_QSATS, ZFLU_VEG, ZGWDCS, ZHQ, ZHTR, ZHU, ZQV, ZWFC,    &
& ZWLMX, ZWWILT)

!     ------------------------------------------------------------------
!     8.- DIFFUSION VERTICALE TURBULENTE
!     ----------------------------------
  
CALL APL_ARPEGE_SURFACE (YDMF_PHYS_BASE_STATE, YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDCPG_MISC,       &
& YDCPG_GPAR, YDMF_PHYS, YDMF_PHYS_SURF, YDVARS, YDMODEL, YDCPG_OPTS%LCONFX, ZALBD, ZALBP, ZALPHA1, ZCDROV,    &
& ZCEROV, ZCFATH, ZCFAU, ZCFBTH, ZCFBU, ZCFBV, ZCHROV, ZCOEFA, ZCOEFN, ZCP, ZDIFEXT, ZDIFWQ, ZDIFWS, &
& ZDQSTS, ZDSA_CPS, ZDSA_LHS, ZDTMSE, ZEDMFQ, ZEDMFS, ZEDMFU, ZEDMFV, ZFLU_CD, ZFLU_CDN, ZFLU_EMIS,  &
& ZFLU_FEVI, ZFLU_NEIJ, ZFLU_QSATS, ZFLU_VEG, ZHQ, ZHTR, ZHU, ZKQLROV, ZKQROV, ZKTROV, ZKUROV, ZLVT, &
& ZMF_UP, ZNEBCH, ZNEBDIFF, ZNEBS, ZPOID, ZQI, ZQICE, ZQL, ZQV, ZRDG_MU0, ZRDG_MU0N, ZRHGMT,         &
& ZSC_FCLL, ZSC_FCLN, ZSC_FEVI, ZSC_FEVN, ZSFSWDIF, ZSFSWDIR, ZSGROUPEL, ZSRAIN, ZSSNOW, ZSTATI,     &
& ZTSN, ZXDROV, ZXHROV, ZXQRO, ZXTRO, ZXTROV, ZXURO, ZXUROV)
    

!-----------------------------------------------------------------------------
!   THE DEEP CONVECTION WILL SEE THE SHALLOW PART FROM KFB AS IT IS WITH LOUIS
!   SCHEME AND THE MODIFIED RI
!----------------------------------------------------------------------------
    
DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    YDMF_PHYS%OUT%DIFTQ   (JLON,JLEV) = YDMF_PHYS%OUT%DIFTQ(JLON,JLEV) + ZDIFCVPPQ(JLON,JLEV)
    YDMF_PHYS%OUT%DIFTS   (JLON,JLEV) = YDMF_PHYS%OUT%DIFTS(JLON,JLEV) + ZDIFCVPPS(JLON,JLEV)
    YDMF_PHYS%OUT%STRTU   (JLON,JLEV) = YDMF_PHYS%OUT%STRTU(JLON,JLEV) + ZDIFCVPPU(JLON,JLEV)
    YDMF_PHYS%OUT%STRTV   (JLON,JLEV) = YDMF_PHYS%OUT%STRTV(JLON,JLEV) + ZDIFCVPPV(JLON,JLEV)
  ENDDO
ENDDO
    


!     DIAGNOSTIC SUPPLEMENTAIRE FLUX DE RAYONNEMENT (PFRTHDS ET PFRSOPT)
!     ADDITIONAL DIAGNOSTICS OF RADIATIVE FLUXES (PFRTHDS AND PFRSOPT)

!DEC$ IVDEP
DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
  IF(.NOT.LMSE)YDMF_PHYS%OUT%FRTHDS(JLON)=YDMF_PHYS%OUT%FRTH(JLON,YDCPG_OPTS%KFLEVG,1)/ZFLU_EMIS(JLON)+YDMODEL%YRCST%RSIGMA*YDMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)**4
  IF(ZRDG_MU0(JLON) <= 0.0_JPRB) THEN
    YDMF_PHYS%OUT%FRSOPT(JLON)=0.0_JPRB
  ELSE
    YDMF_PHYS%OUT%FRSOPT(JLON)=RII0*ZRDG_MU0(JLON)
  ENDIF
ENDDO

! ADDITIONAL DIAGNOSTIC OF THE DERIVATIVE OF THE NON SOLAR SURFACE
! HEAT FLUX WITH RESPECT TO SURFACE TEMPERATURE (PDERNSHF)

IF(LMCC03)THEN
  CALL ACDNSHF(YDMODEL%YRCST, YDPHY, YDPHY1, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KTDIA,          &
  & YDCPG_OPTS%KFLEVG, ZFLU_EMIS, YDMF_PHYS_SURF%GSD_VF%PLSM, ZFLU_NEIJ, ZQV, YDCPG_MISC%QS, YDMF_PHYS_BASE_STATE%YGSP_RR%T, &
  & ZCHROV, ZDQSTS, YDMF_PHYS%OUT%DRNSHF)
ENDIF

!*
!     ------------------------------------------------------------------
!     9.- TRAINEE DES ONDES DE GRAVITE INDUITES PAR LE RELIEF
!     ------------------------------------------------------------------

! GRAVITY WAVE DRAG  
CALL ACDRAG (YDMODEL%YRCST, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTDRAG, YDCPG_OPTS%KFLEVG, &
& PAPRS, PAPRSF, PDELP, ZNBVNO, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, PU, PV, YDVARS%GEOMETRY%RCORI%T0,                       &
& YDMF_PHYS_SURF%GSD_VF%PGETRL, ZGWDCS, YDMF_PHYS_SURF%GSD_VF%PVRLAN, YDMF_PHYS_SURF%GSD_VF%PVRLDI, YDMF_PHYS%OUT%STRDU,         &
& YDMF_PHYS%OUT%STRDV, ZTRAJGWD)
  
 ! SAVE FOR TL/NL COEFS FROM VERT. DIFF AND GWD

  

!     ------------------------------------------------------------------
!     10.- PRECIPITATIONS STRATIFORMES.
!     ---------------------------------

  
! ALARO PRECIPITATION SCHEME
IF ( LSTRAS ) THEN
  CALL ACPLUIS (YDMODEL%YRCST, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTPLUI, YDCPG_OPTS%KFLEVG, &
  & PAPRSF, PDELP, ZNEBS, ZQV, ZQLIS, ZMSC_QW, PR, PT, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLSL,               &
  & YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN                   &
  &                                                  )
ENDIF
  
CALL APL_ARPEGE_DEEP_CONVECTION (YDMF_PHYS_BASE_STATE, YDGEOMETRY, YDCPG_BNDS, YDCPG_OPTS, YDMF_PHYS, &
& YDCPG_DYN0, YDMF_PHYS_SURF, YDVARS, YDCFU, YDMODEL, ZFPCOR, ZQV, ZTENHA, ZTENQVA, ZTENT)

CALL APL_ARPEGE_PRECIPITATION (YDMF_PHYS_BASE_STATE, YDCPG_BNDS, YDCPG_OPTS, YDMF_PHYS,        &
& YDMF_PHYS_SURF, YDVARS, YDMODEL, ZFLU_NEIJ, ZNEBS, ZNEB_CVPP, ZQC_DET_PCMT, ZQI, ZQL, ZQLIS, &
& ZQLI_CVPP, ZQR, ZQS, ZQV, ZSEDIQI, ZSEDIQL, ZTENHA, ZTENQVA, ZVETAF)
    
ZFLX_LOTT_GWU(:, :) = 0.0_JPRB
ZFLX_LOTT_GWV(:, :) = 0.0_JPRB

!*
!     ------------------------------------------------------------------
!         SAUVEGARDE DES FLUX DE PRECIPITATION CONVECTIVE ET STRATIFORME.

IF(LNEBN.OR.LRRGUST) THEN
  DO JLEV=YDCPG_OPTS%KTDIA-1,YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZPFL_FPLCH(JLON,JLEV)=ZFPCOR(JLON,JLEV)
    ENDDO
  ENDDO    
ENDIF

IF(LRRGUST) THEN
  DO JLEV=YDCPG_OPTS%KTDIA-1,YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZPFL_FPLSH(JLON,JLEV)=YDMF_PHYS%OUT%FPLSL(JLON,JLEV)+YDMF_PHYS%OUT%FPLSN(JLON,JLEV)+YDMF_PHYS%OUT%FPLSG(JLON,JLEV)
    ENDDO
  ENDDO
ENDIF

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! UPDATE TRANSPORT FLUXES DUE TO SEDIMENTATION OF CLOUDS.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - -

DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)+ZSEDIQL(JLON,JLEV)
    YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)+ZSEDIQI(JLON,JLEV)
    ZFPLSL (JLON,JLEV)=YDMF_PHYS%OUT%DIFTQL(JLON,JLEV)+YDMF_PHYS%OUT%FPLSL (JLON,JLEV)
    ZFPLSN (JLON,JLEV)=YDMF_PHYS%OUT%DIFTQN(JLON,JLEV)+YDMF_PHYS%OUT%FPLSN (JLON,JLEV)
  ENDDO
ENDDO
  

  ! - - - - - - - - - - - - - - - - -
  ! CORRECT NEGATIVE WATER CONTENTS.
  ! - - - - - - - - - - - - - - - - -

  
CALL QNGCOR (YDMODEL%YRCST, YDPHY2, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, NTPLUI, YDCPG_OPTS%KFLEVG,&
& ZQV, ZQL, ZQI, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN,               &
& YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%FPEVPSL,    &
& YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, &
& YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL,       &
& YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FCQNG )

!*
!     ------------------------------------------------------------------
!     12. - BILAN HYDRIQUE DU SOL
!     ---------------------------

CALL APL_ARPEGE_HYDRO_BUDGET (YDMF_PHYS_BASE_STATE, YDCPG_BNDS, YDCPG_OPTS, YDCPG_GPAR, YDMF_PHYS, &
& YDMF_PHYS_SURF, YDMODEL, ZC3, ZCN, ZDSA_C1, ZDSA_C2, ZFLU_FEVI, ZFLU_NEIJ, ZFLU_VEG, ZFPLSL,     &
& ZFPLSN, ZWFC, ZWLMX, ZWPMX, ZWSEQ, ZWSMX)

!*
!-  --------------------------------------------------------------------
!     13.- DRAG MESOSPHERIQUE POUR UN MODELE POSSEDANT DES NIVEAUX
!               AU-DESSUS DE 50 KM (I.E. DANS LA MESOSPHERE)
!     ------------------------------------------------------------------

! 
! MESOSPHERIC DRAG  
CALL ACDRME ( YDMODEL%YRCST, YDSTA, YDPHY2, YDTOPH, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,      &
& NTDRME, YDCPG_OPTS%KFLEVG, PCP, PDELP, PT, ZQV, PU, PV, YDMF_PHYS%OUT%FRMH, ZMSC_FRMQ, YDMF_PHYS%OUT%STRMU, &
& YDMF_PHYS%OUT%STRMV)
  
!     ------------------------------------------------------------------

  ! STORE THE PSEUDO-HISTORIC SURFACE PRECIPITATION SENSIBLE HEAT FLUX
  ! -------------------------------------------------------------------

   !LPHSPSH

! Store radiative cloudiness in GFL structure for ISP, Historical files or PostProcessing
IF (YIRAD%LGP) YDVARS%IRAD%T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) = YDCPG_MISC%QICE(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
IF (YLRAD%LGP) YDVARS%LRAD%T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) = YDCPG_MISC%QLI (YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
IF (YA%LGP)    YDVARS%A%T1(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) = YDCPG_MISC%NEB (YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)


!*
!     -----------------------------------------------------------------------
!     16.- DDH FLEXIBLES POUR LES CHAMPS DE SURFACE VARIABLES/FLUX/TENDANCES.
!     -----------------------------------------------------------------------

IF (YDMODEL%YRML_DIAG%YRLDDH%LFLEXDIA) THEN
  CALL APLPAR_FLEXDIA (YDCPG_BNDS, YDCPG_OPTS, YDCPG_MISC, YDMF_PHYS, YDMF_PHYS_SURF, YDMODEL, YDDDH, &
  & YDMF_PHYS_BASE_STATE)
ENDIF

IF (LRAFTKE) THEN
     ! DCAPE due to precipitation evaporation.
  CALL ACEVADCAPE(YDMODEL%YRML_PHY_MF%YRPHY2, YDMODEL%YRCST, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMF_PHYS%OUT%FPLSL, &
  & YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, PAPRS, PAPRSF, PT, PCP, PAPHIF,                            &
  & PAPHI, ZDCAPE)
     ! Gusts.
  ZQGM=ZEPSNEB
  CALL ACCLDIA(YDMODEL%YRCST, YDCPG_OPTS%LXCLP, YDCPG_OPTS%LXTGST, YDCPG_OPTS%LXXGST, YDPHY, YDPHY2, YDTOPH, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, &
  & YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG, YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, PU, PV, YDMF_PHYS%OUT%CAPE, ZDCAPE, ZTKE1, PAPHIF,               &
  & YDVARS%GEOMETRY%OROG%T0, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST, ZBLH, IMOC_CLPH)
  YDMF_PHYS%OUT%CLPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=MIN(XMAXLM,MAX(XMINLM,ZBLH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)))
ENDIF

  
CALL ACVISIH(YDMODEL%YRCST, YDMODEL%YRML_PHY_MF%YRPHY%YRDVISI, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON,   &
& YDCPG_OPTS%KTDIA, YDCPG_OPTS%KFLEVG, PAPHI, PAPHIF, PAPRSF, PT, PR, YDCPG_MISC%QLI, YDCPG_MISC%QICE,                &
& ZQR, ZQS, ZQGM, YDMF_PHYS%OUT%VISICLD, YDMF_PHYS%OUT%VISIHYD, YDMF_PHYS%OUT%MXCLWC    )
  

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

ZBAY_QRCONV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) =  0.0_JPRB
ZBAY_QSCONV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) =  0.0_JPRB

DO JLEV=1,YDCPG_OPTS%KFLEVG
  DO JLON= YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA
    ZBAY_QRCONV(JLON,JLEV)  = 0.001_JPRB * ( (ABS(YDMF_PHYS%OUT%FPLCL(JLON,JLEV)) * ZCOEFRAIN(1)) ** ZCOEFRAIN(2) )
    ZBAY_QSCONV(JLON,JLEV)  = 0.001_JPRB * ( (ABS(YDMF_PHYS%OUT%FPLCN(JLON,JLEV)) * ZCOEFSNOW(1)) ** ZCOEFSNOW(2) )    
  ENDDO
ENDDO

   ! Convert density [kg/m3] to mixing ratio [kg/kg]
   ! R_dry (dry air constant)

ZDE2MR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)  = YDMODEL%YRCST%RD * PT(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) / PAPRSF(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
ZBAY_QRCONV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) = ZBAY_QRCONV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) * ZDE2MR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)
ZBAY_QSCONV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) = ZBAY_QSCONV(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:) * ZDE2MR(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,:)

! Precipitation Type

! Compute wet-bulb temperature at 2 meters (suppose homogeneity of qv/ql/qi )
!ZPCLS(KIDIA:KFDIA)=PAPRS(KIDIA:KFDIA,KLEV)-2._JPRB/ZZZF(KIDIA:KFDIA,1,KLEV)*&
!                 &(PAPRS(KIDIA:KFDIA,KLEV)-PAPRSF(KIDIA:KFDIA,KLEV))
CALL PPWETPOINT(YDMODEL%YRCST, YDPHY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, ZPCLS, YDMF_PHYS%OUT%TCLS,                       &
& YDMF_PHYS_BASE_STATE%Q(:, YDCPG_OPTS%KFLEVG), YDMF_PHYS_BASE_STATE%L(:, YDCPG_OPTS%KFLEVG), YDMF_PHYS_BASE_STATE%I(:, YDCPG_OPTS%KFLEVG), &
& YDMF_PHYS%OUT%TPWCLS)

   ! Defined precipitation type
   !
   ! Compute wet-bulb temperature
DO JLEV=1,YDCPG_OPTS%KFLEVG
  CALL PPWETPOINT(YDMODEL%YRCST, YDPHY, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, PAPRSF(:, JLEV), PT(:, JLEV), &
  & YDMF_PHYS_BASE_STATE%Q(:, JLEV), YDMF_PHYS_BASE_STATE%L(:, JLEV), YDMF_PHYS_BASE_STATE%I(:, JLEV),      &
  & ZTW(:, JLEV))
ENDDO

DO JLON=1,YDCPG_OPTS%KLON
  ZFPLS(JLON,YDCPG_OPTS%KFLEVG)=YDMF_PHYS%OUT%FPLCN(JLON,YDCPG_OPTS%KFLEVG)+YDMF_PHYS%OUT%FPLSN(JLON,YDCPG_OPTS%KFLEVG)
  ZFPLC(JLON,YDCPG_OPTS%KFLEVG)=YDMF_PHYS%OUT%FPLCL(JLON,YDCPG_OPTS%KFLEVG)+YDMF_PHYS%OUT%FPLSL(JLON,YDCPG_OPTS%KFLEVG)
  ZFPLSG(JLON,YDCPG_OPTS%KFLEVG)=0._JPRB
ENDDO

  !initialisation de ZZZ
DO JLEV = 1,YDCPG_OPTS%KFLEVG
  DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZZZ(JLON,1,JLEV)=PAPHI(JLON,JLEV)*ZINVG
  ENDDO
ENDDO

  !initialisation de ZDZZ
DO JLEV = 2, YDCPG_OPTS%KFLEVG
  DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZDZZ(JLON,1,JLEV)=ZZZ(JLON,1,JLEV-1)-ZZZ(JLON,1,JLEV)
  ENDDO
ENDDO
DO JLON = YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
  ZDZZ(JLON,1,1)=PAPHI(JLON,0)*ZINVG-ZZZ(JLON,1,1)
ENDDO

CALL DPRECIPS (YDMODEL%YRCST, YDPHY%YRDPRECIPS, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,          &
& YDVARS%GEOMETRY%OROG%T0, YDMF_PHYS%OUT%TPWCLS, YDMF_PHYS%OUT%DIAGH, PAPHIF, ZDZZ, ZTW, YDMF_PHYS_BASE_STATE%L,                 &
& ZFPLC(:, YDCPG_OPTS%KFLEVG), ZFPLS(:, YDCPG_OPTS%KFLEVG), ZFPLSG(:, YDCPG_OPTS%KFLEVG), ZPRC_DPRECIPS(:, YDCPG_OPTS%NDTPRECCUR)&
&             )
  
 !Idem for an other time step and an other period

CALL DPRECIPS(YDMODEL%YRCST, YDPHY%YRDPRECIPS, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDCPG_OPTS%KFLEVG,             &
& YDVARS%GEOMETRY%OROG%T0, YDMF_PHYS%OUT%TPWCLS, YDMF_PHYS%OUT%DIAGH, PAPHIF, ZDZZ, ZTW, YDMF_PHYS_BASE_STATE%L,                   &
& ZFPLC(:, YDCPG_OPTS%KFLEVG), ZFPLS(:, YDCPG_OPTS%KFLEVG), ZFPLSG(:, YDCPG_OPTS%KFLEVG), ZPRC_DPRECIPS2(:, YDCPG_OPTS%NDTPRECCUR2)&
& )

!        2.3  Computes MOCON in the CLP.
!             --------------------------
CALL MF_PHYS_MOCON (YDCPG_BNDS, YDCPG_OPTS, YDMF_PHYS%OUT%MOCON, ZRDG_LCVQ, IMOC_CLPH, &
& YDMF_PHYS_BASE_STATE)

! Store surface water flux P and E for water conservation
IF (YDCPG_OPTS%LCORWAT) THEN
  CALL MF_PHYS_CORWAT (YDCPG_BNDS, YDCPG_OPTS, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, YDMF_PHYS%OUT%FPLCL,                &
  & YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS_SURF%GSD_VH%PEVA, YDMF_PHYS_SURF%GSD_VH%PPCL, &
  & YDMF_PHYS_SURF%GSD_VH%PPCN, YDMF_PHYS_SURF%GSD_VH%PPSL, YDMF_PHYS_SURF%GSD_VH%PPSN)
ENDIF

!        2.6   surface specific humidity necessary to compute the vertical
!              advection of q in the case "delta m=1" (unlagged physics only).
!              ---------------------------------------------------------------

IF (NDPSFI == 1) THEN
  CALL CPQSOL(YDMODEL%YRCST, YDGEOMETRY%YRDIMV, YDPHY, YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_PHY0%PREHYD, &
  & YDMF_PHYS_SURF%GSP_RR%PT_T0, YDCPG_MISC%QS, ZFLU_QSATS, YDCPG_MISC%QSOL)
ENDIF


CALL APL_ARPEGE_ATMOSPHERE_UPDATE (YDMF_PHYS_BASE_STATE, YDMF_PHYS_NEXT_STATE, YDGEOMETRY,            &
& YDCPG_BNDS, YDCPG_OPTS, YDCPG_MISC, YDMF_PHYS, YDCPG_DYN0, YDMF_PHYS_SURF, YDVARS, YDMODEL,         &
& YDCPG_OPTS%LCONFX, YDCPG_OPTS%ZDTPHY, YDDDH, ZDIFEXT, ZDSA_CPS, ZMSC_FRMQ, ZPFL_FEFB1, ZPFL_FEFB2, ZPFL_FEFB3,           &
& ZPFL_FTKE, ZPFL_FTKEI)

CALL MF_PHYS_FPL_PART2 (YDCPG_BNDS, YDCPG_OPTS, ZPFL_FPLCH, ZPFL_FPLSH, YDVARS%CPF%T1, YDVARS%SPF%T1, &
& YDMODEL)

!       2.9b Prognostic convection etc.
!            --------------------------

! TRANSFER NOT ADVECTED VARIABLES INTO PGFLT1
CALL MF_PHYS_TRANSFER (YDCPG_BNDS, YDCPG_OPTS, YDVARS, YDMODEL%YRML_PHY_MF%YRPHY, YDMODEL%YRML_GCONF%YGFL)

!        2.10  Surface variables.
!              ------------------

CALL APL_ARPEGE_SURFACE_UPDATE (YDCPG_BNDS, YDCPG_OPTS, YDCPG_GPAR, YDMF_PHYS, YDMF_PHYS_SURF, &
& YDMODEL, YDCPG_OPTS%LCONFX, YDCPG_OPTS%ZDTPHY, ZDSA_C1, ZDSA_C2, ZFLU_FEVI, ZFLU_VEG)

IF(YDMODEL%YRML_PHY_MF%YRPHY%LCVPGY) THEN
  CALL MF_PHYS_CVV (YDCPG_BNDS, YDCPG_OPTS, YDVARS%CVV%T0, YDVARS%CVV%T1)
ENDIF


!*       5.    Final calculations.
!              -------------------

! * Restore the initial value of some pseudo-historical surface buffers
!   if relevant.
IF (LL_SAVE_PHSURF) THEN
  CALL MF_PHYS_SAVE_PHSURF_PART2 (YDCPG_BNDS, YDCPG_OPTS, ZSAV_DDAL, ZSAV_DDOM, ZSAV_ENTCH,                           &
  & ZSAV_FHPS, ZSAV_GZ0F, ZSAV_GZ0HF, ZSAV_HV, ZSAV_PBLH, ZSAV_QSH, ZSAV_UDAL, ZSAV_UDGRO, ZSAV_UDOM,                 &
  & ZSAV_UNEBH, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VH%PPBLH, YDMF_PHYS_SURF%GSD_VH%PQSH,                  &
  & YDMF_PHYS_SURF%GSD_VH%PSPSH, YDMF_PHYS_SURF%GSD_VK%PUDGRO, YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VV%PZ0H, &
  & YDVARS%DAL%T0, YDVARS%DOM%T0, YDVARS%UAL%T0, YDVARS%UEN%T0, YDVARS%UNEBH%T0, YDVARS%UOM%T0,                       &
  & YDMODEL)
ENDIF

! Store horizontal exchange coefficients (3D turbulence) to SL2 buffers

CALL MF_PHYS_BAYRAD (YDCPG_BNDS, YDCPG_OPTS, ZBAY_QRCONV, ZBAY_QSCONV, YDVARS%RCONV%T1, YDVARS%SCONV%T1, &
& YDMODEL, YDCPG_OPTS%LAROME)

IF (LEDR) THEN
  YDMF_PHYS_SURF%GSD_DI%PXEDR(:,:)=YDMF_PHYS%OUT%EDR(:,:)
ENDIF

CALL MF_PHYS_PRECIPS (YDCPG_BNDS, YDCPG_OPTS, ZPRC_DPRECIPS, ZPRC_DPRECIPS2, YDMF_PHYS_SURF%GSD_XP%PPRECIP, &
& YDMF_PHYS_SURF%GSD_XP2%PPRECIP2, YDMODEL)

END ASSOCIATE
END ASSOCIATE
END ASSOCIATE

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('APL_ARPEGE', 1, ZHOOK_HANDLE)

END SUBROUTINE APL_ARPEGE
