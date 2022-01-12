#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE APLPAR(YDGEOMETRY, YDCPG_DIM, YDCPG_MISC, YDCPG_GPAR, YDCPG_PHY0, YDCPG_PHY9, YDMF_PHYS, &
& YDCPG_DYN0, YDCPG_DYN9, YDMF_PHYS_SURF, YDCPG_SL1, YDCPG_SL2, YDVARS, YDGMV, YDSURF, YDCFU, &
& YDXFU, YDMODEL, LDCONFX, PDTPHY, PGFL, PGMVT1, PGFLT1, PTRAJ_PHYS, YDDDH)

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

! PGFL       : GFL FIELDS
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
                              & CPG_SL1_TYPE, CPG_SL2_TYPE, CPG_GPAR_TYPE
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
USE YOMCT0             , ONLY : LCALLSFX ,LSFORCS, LELAM
USE YOMDYNA            , ONLY : L3DTURB
USE YOMRIP0            , ONLY : NINDAT
USE DDH_MIX            , ONLY : TYP_DDH
USE YOMLUN             , ONLY : NULOUT
USE YOMLSFORC          , ONLY : LMUSCLFA,NMUSCLFA
USE YOMTRAJ            , ONLY : TRAJ_PHYS_TYPE
USE YOMCFU             , ONLY : TCFU !!! for parameters of FLASH
USE SPP_MOD  , ONLY : YSPP, YSPP_CONFIG
USE MF_PHYS_BASE_STATE_TYPE_MOD &
                       , ONLY : MF_PHYS_BASE_STATE_TYPE
USE MF_PHYS_NEXT_STATE_TYPE_MOD &
                       , ONLY : MF_PHYS_NEXT_STATE_TYPE


USE CPG_TYPE_MOD       , ONLY : CPG_PHY_TYPE
USE YOMGMV             , ONLY : TGMV
USE SC2PRG_MOD         , ONLY : SC2PRG

USE YOMCT0             , ONLY : LTWOTL, LAROME, LCORWAT
USE YOMNUD             , ONLY : NFNUDG   ,LNUDG 
USE YOMSNU             , ONLY : XPNUDG
USE YOMSCM             , ONLY : LGSCM
USE YOMCHET            , ONLY : GCHETN
USE YOMTRAJ            , ONLY : LPRTTRAJ

USE INTFLEX_MOD        , ONLY : LINTFLEX, TYPE_INTPROCSET, NEWINTPROCSET, CLEANINTPROCSET
!     -------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_DIM_TYPE),INTENT(IN)    :: YDCPG_DIM
TYPE(CPG_MISC_TYPE),INTENT(INOUT):: YDCPG_MISC
TYPE(CPG_GPAR_TYPE),INTENT(INOUT):: YDCPG_GPAR
TYPE(CPG_PHY_TYPE),INTENT(INOUT), TARGET :: YDCPG_PHY0
TYPE(CPG_PHY_TYPE),INTENT(INOUT), TARGET :: YDCPG_PHY9
TYPE(MF_PHYS_TYPE),INTENT(INOUT) :: YDMF_PHYS
TYPE(CPG_DYN_TYPE),INTENT(INOUT), TARGET :: YDCPG_DYN0
TYPE(CPG_DYN_TYPE),INTENT(INOUT), TARGET :: YDCPG_DYN9
TYPE(MF_PHYS_SURF_TYPE),INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(CPG_SL1_TYPE),INTENT(INOUT), TARGET :: YDCPG_SL1
TYPE(CPG_SL2_TYPE),INTENT(INOUT) :: YDCPG_SL2
TYPE(FIELD_VARIABLES),INTENT(INOUT), TARGET :: YDVARS
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
TYPE(TCFU)        ,INTENT(INOUT) :: YDCFU
TYPE(TXFU)        ,INTENT(INOUT) :: YDXFU
TYPE(MODEL)       ,INTENT(INOUT), TARGET :: YDMODEL
LOGICAL           ,INTENT(IN)    :: LDCONFX
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTPHY 
REAL(KIND=JPRB)   ,INTENT(INOUT), TARGET :: PGFL(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDGMV%YT1%NDIM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLT1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1)
 

TYPE (TRAJ_PHYS_TYPE), INTENT(INOUT) :: PTRAJ_PHYS
TYPE(TYP_DDH)     ,INTENT(INOUT) :: YDDDH

!     ------------------------------------------------------------------
LOGICAL :: LL_SAVE_PHSURF
LOGICAL :: LLXFUMSE

INTEGER(KIND=JPIM) :: IFIELDSS

INTEGER(KIND=JPIM) :: INSTEP_DEB,INSTEP_FIN
INTEGER(KIND=JPIM) :: JROF
INTEGER(KIND=JPIM) :: ISLB1U9  ,ISLB1V9  ,ISLB1T9  ,ISLB1GFL9, ISLB1VD9

!     --- UPPER AIR PHYSICAL TENDENCIES.
REAL(KIND=JPRB) :: ZTENDH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)     ! Enthalpy tendency.
REAL(KIND=JPRB) :: ZTENDQ(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)     ! Moisture tendency.
REAL(KIND=JPRB) :: ZTENDPTKE(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)  ! Pseudo progn. TKE

! GFL tendencies for APL_AROME (assumes YDMODEL%YRML_GCONF%YGFL%NUMFLDS>=YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
! for now, use Jovi's trick :
REAL(KIND=JPRB) :: ZTENDGFL(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NUMFLDS)   ! GFL tendencies

!     --- UPPER AIR PHYSICAL TENDENCIES FOR AROME.
!       (the previous one are not used in AROME)
REAL(KIND=JPRB) :: ZTENDT (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)        ! temperature tendency
REAL(KIND=JPRB) :: ZTENDD (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)        ! d  tendency
REAL(KIND=JPRB) :: ZTENDEXT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)      ! GFL EXTRA tendency
REAL(KIND=JPRB) :: ZTENDEXT_DEP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)  ! GFL EXTRA tendency
REAL(KIND=JPRB) :: ZDIFEXT(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)     ! Extra-GFL fluxes.

REAL(KIND=JPRB) :: ZTENDU (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! U tendency without deep convection contribution
REAL(KIND=JPRB) :: ZTENDV (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! V tendency without deep convection contribution

!     --- RADIATION COEFFICIENTS FOR SIMPLIFIED PHYSICS IN GRID-POINT ---
REAL(KIND=JPRB) :: ZAC(YDCPG_DIM%KLON,(YDCPG_DIM%KFLEVG+1)*(YDCPG_DIM%KFLEVG+1))   ! Curtis matrix.
REAL(KIND=JPRB) :: ZAC_HC(YDCPG_DIM%KFLEVG+1,YDCPG_DIM%KFLEVG+1)           ! horizontally-constant field for ZAC.



! required for INTFLEX
TYPE(TYPE_INTPROCSET) :: YLPROCSET

! SPP
REAL(KIND=JPRB) :: ZGP2DSPP(YDCPG_DIM%KLON,YSPP%N2D)


REAL(KIND=JPRB), POINTER :: ZPTENDEFB11(:,:), ZPTENDEFB21(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDEFB31(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDG1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDICONV1(:,:), ZPTENDI1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDLCONV1(:,:)
REAL(KIND=JPRB), POINTER :: ZP1EZDIAG(:,:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDQ1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDRCONV1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDR1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDSCONV1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDS1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDTKE1(:,:)
REAL(KIND=JPRB), POINTER :: ZPTENDL1(:,:)

TYPE (MF_PHYS_BASE_STATE_TYPE) :: YLMF_PHYS_BASE_STATE
TYPE (MF_PHYS_NEXT_STATE_TYPE) :: YLMF_PHYS_NEXT_STATE

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
REAL(KIND=JPRB) :: ZALBD(YDCPG_DIM%KLON,YDCPG_DIM%KSW), ZALBP(YDCPG_DIM%KLON,YDCPG_DIM%KSW),&
                 & ZALB(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZSFSWDIR (YDCPG_DIM%KLON,YDCPG_DIM%KSW), ZSFSWDIF (YDCPG_DIM%KLON,YDCPG_DIM%KSW)
REAL(KIND=JPRB) :: ZTRSODIR (YDCPG_DIM%KLON,YDCPG_DIM%KSW), ZTRSODIF (YDCPG_DIM%KLON,YDCPG_DIM%KSW)
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
REAL(KIND=JPRB) :: ZALBPMER

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
REAL(KIND=JPRB), DIMENSION(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDCPG_DIM%KSW):: ZPIZA_DST !Single scattering
                                             ! albedo of dust (points,lev,wvl)
REAL(KIND=JPRB), DIMENSION(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDCPG_DIM%KSW):: ZCGA_DST  !Assymetry factor
                                             ! for dust (points,lev,wvl)
REAL(KIND=JPRB), DIMENSION(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDCPG_DIM%KSW):: ZTAUREL_DST !tau/tau_{550}
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

INTEGER (KIND=JPIM)  :: IMOC_CLPH (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZADJ_DTAJU (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZADJ_TAUX (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
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
REAL (KIND=JPRB)     :: ZFLU_QS1 (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZFLU_QSATS (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZFLU_QSAT (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZFLU_VEG (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZKUR_KTROV_H (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZKUR_KUROV_H (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
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
REAL (KIND=JPRB)     :: ZRDT_COR (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAB3C (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAB3N (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAB4C (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAB4N (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAB6C (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAB6N (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT1C (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT1N (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT2C (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT2N (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT3C (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT3N (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT4C (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT4N (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT5C (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZRDT_RAT5N (YDCPG_DIM%KLON, 1:YDCPG_DIM%KFLEVG)
REAL (KIND=JPRB)     :: ZTDS_TDALBNS (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDRHONS (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDSNS (YDCPG_DIM%KLON)
REAL (KIND=JPRB)     :: ZTDS_TDTP (YDCPG_DIM%KLON, 1:YDCPG_DIM%KSBDLEV)
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
REAL (KIND=JPRB)     :: ZRDG_MMU0 (YDCPG_DIM%KLON)
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
#include "acaa1.intfb.h"
#include "acajucv.intfb.h"
#include "accdev.intfb.h"
#include "accldia.intfb.h"
#include "acclph.intfb.h"
#include "accoefk.intfb.h"
#include "accsu.intfb.h"
#include "accvimpgy.intfb.h"
#include "accvimp.intfb.h"
#include "accvimp_v3.intfb.h"
#include "accvud.intfb.h"
#include "acdayd.intfb.h"
#include "acdifoz.intfb.h"
#include "acdifus.intfb.h"
#include "acdifv1.intfb.h"
#include "acdifv2.intfb.h"
#include "acdifv3.intfb.h"
#include "acdnshf.intfb.h"
#include "acdrac.intfb.h"
#include "acdrag.intfb.h"
#include "acdrme.intfb.h"
#include "acdrov.intfb.h"
#include "acevadcape.intfb.h"
#include "acfluso.intfb.h"
#include "achmt.intfb.h"
#include "achmtls.intfb.h"
#include "acmixelen.intfb.h"
#include "acmixlentm.intfb.h"
#include "acmixlenz.intfb.h"
#include "acmodo.intfb.h"
#include "acmrip.intfb.h"
#include "acmris.intfb.h"
#include "acmriss.intfb.h"
#include "acnebc.intfb.h"
#include "acnebcond.intfb.h"
#include "acnebn.intfb.h"
#include "acnebr.intfb.h"
#include "acnorgwd.intfb.h"
#include "acnpart.intfb.h"
#include "acnsdo.intfb.h"
#include "acozone.intfb.h"
#include "acpblh.intfb.h"
#include "acpblhtm.intfb.h"
#include "acpcmt.intfb.h"
#include "acpluie.intfb.h"
#include "acpluis.intfb.h"
#include "acpluiz.intfb.h"
#include "acptke.intfb.h"
#include "acradcoef.intfb.h"
#include "acraneb2.intfb.h"
#include "acraneb.intfb.h"
#include "acrso.intfb.h"
#include "acsol.intfb.h"
#include "actkecoefkh.intfb.h"
#include "actkecoefk.intfb.h"
#include "actkehmt.intfb.h"
#include "actke.intfb.h"
#include "actkezotls.intfb.h"
#include "actqsat.intfb.h"
#include "actqsats.intfb.h"
#include "acupd.intfb.h"
#include "acupm.intfb.h"
#include "acuptq.intfb.h"
#include "acupu.intfb.h"
#include "acveg.intfb.h"
#include "acvisih.intfb.h"
#include "acvppkf.intfb.h"
#include "aplmphys.intfb.h"
#include "aplpar2intflex.intfb.h"
#include "aplpar_init.intfb.h"
#include "aro_ground_diag_2isba.h"
#include "aro_ground_diag.h"
#include "aro_ground_diag_z0.h"
#include "aro_ground_param.h"
#include "aro_mnhdust.h"
#include "arp_ground_param.intfb.h"
#include "checkmv.intfb.h"
!include "chem_main.intfb.h"
#include "cpchet.intfb.h"
#include "cpmvvps.intfb.h"
#include "cpnudg.intfb.h"
#include "cpozo.intfb.h"
#include "cpphinp.intfb.h"
#include "cpqsol.intfb.h"
#include "cptend_flex.intfb.h"
#include "cptend_new.intfb.h"
#include "cptends.intfb.h"
#include "cputqy_aplpar_expl.intfb.h"
#include "cputqy_aplpar_loop.intfb.h"
#include "cpwts.intfb.h"
#include "cucalln_mf.intfb.h"
#include "culight.intfb.h"
#include "dprecips.intfb.h"
#include "mean_rad_temp.intfb.h"
#include "mf_phys_bayrad.intfb.h"
#include "mf_phys_corwat.intfb.h"
#include "mf_phys_cvv.intfb.h"
#include "mf_phys_fpl_part1.intfb.h"
#include "mf_phys_fpl_part2.intfb.h"
#include "mf_phys_init.intfb.h"
#include "mf_phys_mocon.intfb.h"
#include "mf_phys_precips.intfb.h"
#include "mf_phys_save_phsurf_part1.intfb.h"
#include "mf_phys_save_phsurf_part2.intfb.h"
#include "mf_phys_transfer.intfb.h"
#include "ppwetpoint.intfb.h"
#include "profilechet.intfb.h"
#include "qngcor.intfb.h"
#include "radaer.intfb.h"
#include "radheat.intfb.h"
#include "radozcmf.intfb.h"
#include "recmwf.intfb.h"
#include "suozon.intfb.h"
#include "surf_ideal_flux.intfb.h"
#include "writephysio.intfb.h"
#include "wrphtrajm.intfb.h"
#include "wrphtrajtm_nl.intfb.h"
#include "wrradcoef.intfb.h"
#include "wrscmr.intfb.h"
#include "aplpar_flexdia.intfb.h"

!     ------------------------------------------------------------------

#include "fcttrm.func.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('APLPAR', 0, ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDVAB=>YDGEOMETRY%YRVAB, &
& YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,           &
& YDPTRSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1, YDPTRSLB2=>YDMODEL%YRML_DYN%YRPTRSLB2, YDTOPH=>YDMODEL%YRML_PHY_MF%YRTOPH,                   &
& YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL, YDRIP=>YDMODEL%YRML_GCONF%YRRIP, YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH,                          &
& YDRCOEF=>YDMODEL%YRML_PHY_RAD%YRRCOEF, YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY, YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY,                     &
& YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, YGFL=>YDMODEL%YRML_GCONF%YGFL,                                &
& YDEPHY=> YDMODEL%YRML_PHY_EC%YREPHY, YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR, YDPRECIPS=>YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS,         &
& YDGEM=>YDGEOMETRY%YRGEM, YDSTA=>YDGEOMETRY%YRSTA, YDERDI=>YDMODEL%YRML_PHY_RAD%YRERDI, YDMCC=>YDMODEL%YRML_AOC%YRMCC,               &
& YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD, YDPHY3=>YDMODEL%YRML_PHY_MF%YRPHY3, YDPHY1=>YDMODEL%YRML_PHY_MF%YRPHY1,                        &
& YDPHY0=>YDMODEL%YRML_PHY_MF%YRPHY0, YDNORGWD=>YDMODEL%YRML_PHY_MF%YRNORGWD, YDMSE=>YDMODEL%YRML_PHY_MF%YRMSE,                       &
& YDPHYDS=>YDMODEL%YRML_PHY_MF%YRPHYDS    )


ASSOCIATE(CMF_UPDRAFT=>YDPARAR%CMF_UPDRAFT, TSPHY=>YDPHY2%TSPHY, NTSSG=>YDDPHY%NTSSG, LMDUST=>YDARPHY%LMDUST,                  &
& LMSE=>YDARPHY%LMSE, YI=>YGFL%YI, YEZDIAG=>YGFL%YEZDIAG, YL     =>YGFL%YL, YEXT=>YGFL%YEXT, YG=>YGFL%YG,                      &
& YQ=>YGFL%YQ, YR=>YGFL%YR, YSCONV=>YGFL%YSCONV, YS=>YGFL%YS, YEFB3=>YGFL%YEFB3, YEFB2=>YGFL%YEFB2, YEFB1=>YGFL%YEFB1,         &
& LCHEM_ARPCLIM=>YDMODEL%YRML_CHEM%YRCHEM%LCHEM_ARPCLIM, NGFL_EXT=>YGFL%NGFL_EXT, YTKE=>YGFL%YTKE, YLCONV=>YGFL%YLCONV,        &
& YRCONV=>YGFL%YRCONV, YICONV=>YGFL%YICONV, YSP_SBD=>YDSURF%YSP_SBD, LTRAJPS=>YDSIMPHL%LTRAJPS, LNEBR=>YDPHY%LNEBR,            &
& LNEBN=>YDPHY%LNEBN, LSTRAPRO=>YDPHY%LSTRAPRO, LPTKE=> YDPHY%LPTKE, NDPSFI=>YDPHY%NDPSFI, LOZONE=>YDPHY%LOZONE,               &
& L3MT=>YDPHY%L3MT, LGPCMT=>YDPHY%LGPCMT, LAJUCV=>YDPHY%LAJUCV, LCVPGY=>YDPHY%LCVPGY, LRRGUST=>YDPHY%LRRGUST,                  &
& LEDR=>YDPHY%LEDR, NTAJUC=> YDTOPH%NTAJUC, NTPLUI=>YDTOPH%NTPLUI, LDPRECIPS=>YDPHY%LDPRECIPS, LDPRECIPS2=>YDPHY%LDPRECIPS2,   &
& NDTPRECCUR=>YDPRECIPS%NDTPRECCUR, NDTPRECCUR2=>YDPRECIPS%NDTPRECCUR2, LRCOEF    =>YDRCOEF%LRCOEF, NG3SR=>YDRCOEF%NG3SR,      &
& XMINLM=>YDPHY0%XMINLM, RTCAPE=>YDPHY0%RTCAPE, GRSO=>YDPHY0%GRSO, GCVTSMO=>YDPHY0%GCVTSMO, XKLM=>YDPHY0%XKLM,                 &
& GAEPS=>YDPHY0%GAEPS, AERCS1=>YDPHY0%AERCS1, AERCS3=>YDPHY0%AERCS3, AERCS5=>YDPHY0%AERCS5, HUTIL2=>YDPHY0%HUTIL2,             &
& HUTIL1=>YDPHY0%HUTIL1, XMAXLM=>YDPHY0%XMAXLM, TEQK=>YDPHY0%TEQK, HUCOE=>YDPHY0%HUCOE, LCVNHD=>YDPHY0%LCVNHD,                 &
& TEQC=>YDPHY0%TEQC, UHDIFV=>YDPHY0%UHDIFV, HUTIL=>YDPHY0%HUTIL, NPCLO1=>YDPHY0%NPCLO1, NPCLO2=>YDPHY0%NPCLO2,                 &
& RDECRD=>YDPHY0%RDECRD, RDECRD1=>YDPHY0%RDECRD1, RDECRD2=>YDPHY0%RDECRD2, RDECRD3=>YDPHY0%RDECRD3, RDECRD4=>YDPHY0%RDECRD4,   &
& ETKE_MIN=>YDPHY0%ETKE_MIN, HSOLIWR=>YDPHY1%HSOLIWR, ALCRIN=>YDPHY1%ALCRIN, ALBMED=>YDPHY1%ALBMED, WSMX=>YDPHY1%WSMX,         &
& LALBMERCLIM=>YDPHY1%LALBMERCLIM, HSOLIT0=>YDPHY1%HSOLIT0, HSOL=>YDPHY1%HSOL, WPMX=>YDPHY1%WPMX, EMCRIN=>YDPHY1%EMCRIN,       &
& EMMMER=>YDPHY1%EMMMER, TMERGL=>YDPHY1%TMERGL, EMMGLA=>YDPHY1%EMMGLA, LRAFTKE=>YDPHY2%LRAFTKE, LRAFTUR=>YDPHY2%LRAFTUR,       &
& HVCLS=>YDPHY2%HVCLS, HTCLS=>YDPHY2%HTCLS, FSM_HH=>YDPHY3%FSM_HH, FSM_GG=>YDPHY3%FSM_GG, FSM_FF=>YDPHY3%FSM_FF,               &
& FSM_EE=>YDPHY3%FSM_EE, FSM_II=>YDPHY3%FSM_II, FSM_CC=>YDPHY3%FSM_CC, FSM_DD=>YDPHY3%FSM_DD, RLAMB_WATER=>YDPHY3%RLAMB_WATER, &
& RII0=>YDPHY3%RII0, QCO2=>YDPHY3%QCO2, RLAMB_SOLID=>YDPHY3%RLAMB_SOLID, NDLUNG=>YDDIM%NDLUNG, NDGUNG=>YDDIM%NDGUNG,           &
& NDLUXG=>YDDIM%NDLUXG, NDGUXG=>YDDIM%NDGUXG, LRDEPOS=>YDARPHY%LRDEPOS, LMPA=>YDARPHY%LMPA, CCOUPLING=>YDARPHY%CCOUPLING,      &
& LSURFEX_KFROM=>YDARPHY%LSURFEX_KFROM, YA=>YGFL%YA, NGFL_EZDIAG=>YGFL%NGFL_EZDIAG, YFQTUR=>YGFL%YFQTUR,                       &
& YFSTUR=>YGFL%YFSTUR, YIRAD=>YGFL%YIRAD, YLRAD=>YGFL%YLRAD, XZSEPS=>YDMSE%XZSEPS, LVDIFSPNL=>YDSIMPHL%LVDIFSPNL,              &
& LGWDSPNL=>YDSIMPHL%LGWDSPNL, LRAYSP=>YDSIMPHL%LRAYSP, LSTRA=>YDPHY%LSTRA, LAEROSOO=>YDPHY%LAEROSOO,                          &
& LCDDPRO=>YDPHY%LCDDPRO, LRKCDEV=>YDPHY%LRKCDEV, LHUCN=>YDPHY%LHUCN, LCOEFK_TOMS=>YDPHY%LCOEFK_TOMS,                          &
& LVDIF=>YDPHY%LVDIF, LRRMES=>YDPHY%LRRMES, LCVRA=>YDPHY%LCVRA, LCVTDK=>YDPHY%LCVTDK, LCOEFK_RIS=>YDPHY%LCOEFK_RIS,            &
& LAEROLAN=>YDPHY%LAEROLAN, LNEBECT=>YDPHY%LNEBECT, LAERODES=>YDPHY%LAERODES, LNEWSTAT=>YDPHY%LNEWSTAT,                        &
& LTHERMO=>YDPHY%LTHERMO, LO3FL=>YDPHY%LO3FL, LPHSPSH=>YDPHY%LPHSPSH, LSNV=>YDPHY%LSNV, LECSHAL=>YDPHY%LECSHAL,                &
& LECT=>YDPHY%LECT, LDIFCONS=>YDPHY%LDIFCONS, LNODIFQC=>YDPHY%LNODIFQC, LAEROVOL=>YDPHY%LAEROVOL, LRSTAER=>YDPHY%LRSTAER,      &
& NCALLRAD=>YDPHY%NCALLRAD, LNCVPGY=>YDPHY%LNCVPGY, LRAYLU=>YDPHY%LRAYLU, LAEROSUL=>YDPHY%LAEROSUL, LO3ABC=>YDPHY%LO3ABC,      &
& LSTRAS=>YDPHY%LSTRAS, LCOEFK_PTTE=>YDPHY%LCOEFK_PTTE, LSFHYD=>YDPHY%LSFHYD, LAEROSEA=>YDPHY%LAEROSEA,                        &
& NDIFFNEB=>YDPHY%NDIFFNEB, LEDKF=>YDPHY%LEDKF, LRCVOTT=>YDPHY%LRCVOTT, LMPHYS=>YDPHY%LMPHYS, LZ0HSREL=>YDPHY%LZ0HSREL,        &
& LCAMOD=>YDPHY%LCAMOD, LCOMOD=>YDPHY%LCOMOD, LCONDWT=>YDPHY%LCONDWT, LCVCSD=>YDPHY%LCVCSD, LNSDO=>YDPHY%LNSDO,                &
& LUDEVOL=>YDPHY%LUDEVOL, LRAY=>YDPHY%LRAY, LGWD=>YDPHY%LGWD, LCOEFKTKE=>YDPHY%LCOEFKTKE, LRAYFM=>YDPHY%LRAYFM,                &
& LECDEEP=>YDPHY%LECDEEP, LCVGQD=>YDPHY%LCVGQD, LGWDC=>YDPHY%LGWDC, LNORGWD=>YDPHY%LNORGWD, LFLUSO=>YDPHY%LFLUSO,              &
& LNEBCO=>YDPHY%LNEBCO, LNEBCV=>YDPHY%LNEBCV, LSOLV=>YDPHY%LSOLV, LCOEFKSURF=>YDPHY%LCOEFKSURF, LRNUEXP=>YDPHY%LRNUEXP,        &
& NRAY=>YDPHY%NRAY, LDAYD=>YDPHY%LDAYD, LEDMFI=>YDPHY%LEDMFI, LFPCOR=>YDPHY%LFPCOR, LRPROX=>YDPHY%LRPROX,                      &
& LPROCLD=>YDPHY%LPROCLD, LACDIFUS=>YDPHY%LACDIFUS, LCAPE=>YDPHY%LCAPE, LCVRAV3=>YDPHY%LCVRAV3, LHMTO=>YDPHY%LHMTO,            &
& LVGSN=>YDPHY%LVGSN, LCVPPKF=>YDPHY%LCVPPKF, LCVPRO=>YDPHY%LCVPRO, LADJCLD=>YDPHY%LADJCLD, LGRAPRO=>YDPHY%LGRAPRO,            &
& RDECLI=>YDRIP%RDECLI, CMF_CLOUD=>YDPARAR%CMF_CLOUD, XSW_BANDS=>YDPARAR%XSW_BANDS, NSWB_MNH=>YDPARAR%NSWB_MNH,                &
& LMIXUV=>YDPARAR%LMIXUV, NTRADI=>YDTOPH%NTRADI, NTNEBU=>YDTOPH%NTNEBU, NTDIFU=>YDTOPH%NTDIFU, NTOZON=>YDTOPH%NTOZON,          &
& NTDRME=>YDTOPH%NTDRME, NTCVIM=>YDTOPH%NTCVIM, NTCOET=>YDTOPH%NTCOET, NTDRAG=>YDTOPH%NTDRAG, RMESOQ=>YDTOPH%RMESOQ,           &
& RMESOT=>YDTOPH%RMESOT, RMESOU=>YDTOPH%RMESOU, NTQSAT=>YDTOPH%NTQSAT, NTCOEF=>YDTOPH%NTCOEF, NAER=>YDERAD%NAER,               &
& NOZOCL=>YDERAD%NOZOCL, NRADFR=>YDERAD%NRADFR, NSW=>YDERAD%NSW, RCARDI=>YDERDI%RCARDI, RSUNDUR=>YDERDI%RSUNDUR,               &
& LXVISI=>YDXFU%LXVISI, LXVISI2=>YDXFU%LXVISI2, LMCC03=>YDMCC%LMCC03, NSTOP=>YDRIP%NSTOP, RCODEC=>YDRIP%RCODEC,                &
& RHGMT=>YDRIP%RHGMT, RSIDEC=>YDRIP%RSIDEC, RSOVR=>YDRIP%RSOVR, RSTATI=>YDRIP%RSTATI, TSTEP=>YDRIP%TSTEP,                      &
& NDTPREC=>YDPHY%YRDPRECIPS%NDTPREC, NDTPREC2=>YDPHY%YRDPRECIPS%NDTPREC2, STPRE=>YDSTA%STPRE, STPREH=>YDSTA%STPREH,            &
& STTEM=>YDSTA%STTEM, NORGWD_NNOVERDIF=>YDNORGWD%NORGWD_NNOVERDIF, LGCHECKMV=>YDPHY%LGCHECKMV, NAERO=>YGFL%NAERO,              &
& NCHEM=>YDMODEL%YRML_GCONF%YGFL%NCHEM, NACTAERO=>YGFL%NACTAERO, LXMRT=>YDXFU%LXMRT, RDELXN=>YDGEM%RDELXN,                     &
& LFLASH =>YDCFU%LFLASH)

CALL SC2PRG(YEFB1%MP1, ZTENDGFL, ZPTENDEFB11)
CALL SC2PRG(YEFB2%MP1, ZTENDGFL, ZPTENDEFB21)
CALL SC2PRG(YEFB3%MP1, ZTENDGFL, ZPTENDEFB31)
CALL SC2PRG(YG%MP1, ZTENDGFL, ZPTENDG1)
CALL SC2PRG(YICONV%MP1, ZTENDGFL, ZPTENDICONV1)
CALL SC2PRG(YI%MP1, ZTENDGFL, ZPTENDI1)
CALL SC2PRG(YLCONV%MP1, ZTENDGFL, ZPTENDLCONV1)
CALL SC2PRG(YL%MP1, ZTENDGFL, ZPTENDL1)
CALL SC2PRG(YQ%MP1, ZTENDGFL, ZPTENDQ1)
CALL SC2PRG(YRCONV%MP1, ZTENDGFL, ZPTENDRCONV1)
CALL SC2PRG(YR%MP1, ZTENDGFL, ZPTENDR1)
CALL SC2PRG(YSCONV%MP1, ZTENDGFL, ZPTENDSCONV1)
CALL SC2PRG(YS%MP1, ZTENDGFL, ZPTENDS1)
CALL SC2PRG(YTKE%MP1, ZTENDGFL, ZPTENDTKE1)

CALL SC2PRG(1, YEZDIAG(:)%MP, YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG, PGFL, ZP1EZDIAG)

CALL YLMF_PHYS_BASE_STATE%INIT (LTWOTL, YDCPG_DYN0, YDCPG_DYN9, YDCPG_PHY0, YDCPG_PHY9, YDVARS, &
& YDMF_PHYS_SURF, PGFL=PGFL, YDMODEL=YDMODEL)

CALL YLMF_PHYS_NEXT_STATE%INIT (YDCPG_SL1, YDVARS, YDMODEL)

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

CALL CPPHINP(YDGEOMETRY, YDMODEL, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDVARS%GEOMETRY%GEMU%T0, YDVARS%GEOMETRY%GELAM%T0, YDVARS%U%T0, &
& YDVARS%V%   T0, YDVARS%Q%T0, YDVARS%Q%DL, YDVARS%Q%DM, YDVARS%CVGQ%DL, YDVARS%CVGQ%DM, YDCPG_PHY0%XYB%RDELP,  &
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


CALL MF_PHYS_INIT ( YGFL, YDARPHY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,                  &
& NTSSG, YSP_SBD%NLEVS, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS,           &
& YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%FCCQL,            &
& YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%FCQNNG,             &
& YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FCQGNG, YDMF_PHYS%OUT%FPLCL,          &
& YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLCG, YDMF_PHYS%OUT%FPLCH, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN,              &
& YDMF_PHYS%OUT%FPLSG, YDMF_PHYS%OUT%FPLSH, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRSOC, YDMF_PHYS%OUT%FRTH,                &
& YDMF_PHYS%OUT%FRTHC, YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, YDMF_PHYS%OUT%STRDU, YDMF_PHYS%OUT%STRDV,              &
& YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV, YDMF_PHYS%OUT%STRMU, YDMF_PHYS%OUT%STRMV, YDMF_PHYS%OUT%FRMH,               &
& YDMF_PHYS%OUT%DIFCQLC, YDMF_PHYS%OUT%DIFCQIC, YDMF_PHYS%OUT%FIMCC, YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC,        &
& YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC, YDMF_PHYS%OUT%FCNEGQIC, YDMF_PHYS%OUT%FCNEGQRC,   &
& YDMF_PHYS%OUT%FCNEGQSC, YDMF_PHYS%OUT%FCHOZ, YDCPG_MISC%NEB, YDCPG_MISC%QICE, YDCPG_MISC%QLI, YDCPG_MISC%RH,          &
& YDMF_PHYS%OUT%ALB, YDMF_PHYS%OUT%CT, YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FCLL, YDMF_PHYS%OUT%FCLN, YDMF_PHYS%OUT%FCS,  &
& YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FLASH, YDMF_PHYS%OUT%FTR,                   &
& YDMF_PHYS%OUT%FLWSP, YDMF_PHYS%OUT%FONTE, YDMF_PHYS%OUT%FGEL, YDMF_PHYS%OUT%FGELS, YDMF_PHYS%OUT%FRSGNI,              &
& YDMF_PHYS%OUT%FRSDNI, YDMF_PHYS%OUT%FRSODS, YDMF_PHYS%OUT%FRSOPS, YDMF_PHYS%OUT%FRSOPT, YDMF_PHYS%OUT%FRSOLU,         &
& YDMF_PHYS%OUT%FRTHDS, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPFPCL,         &
& YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG, YDMF_PHYS%OUT%FPEVPCL,     &
& YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPEVPCG, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, YDCPG_MISC%QS,                   &
& YDMF_PHYS%OUT%RUISL, YDMF_PHYS%OUT%RUISP, YDMF_PHYS%OUT%RUISS, YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS,                &
& YDMF_PHYS%OUT%NUCLS, YDMF_PHYS%OUT%NVCLS, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%MRT, YDMF_PHYS%OUT%QCLS,                  &
& YDMF_PHYS%OUT%RHCLS, YDCPG_MISC%CLCT, YDMF_PHYS%OUT%CLCH, YDMF_PHYS%OUT%CLCM, YDMF_PHYS%OUT%CLCL, YDMF_PHYS%OUT%CLCC, &
& YDMF_PHYS%OUT%CAPE, YDMF_PHYS%OUT%CTOP, YDMF_PHYS%OUT%CLPH, YDMF_PHYS%OUT%VEIN, YDMF_PHYS%OUT%UGST,                   &
& YDMF_PHYS%OUT%VGST, YDMF_PHYS%OUT%DIAGH, YDMF_PHYS%OUT%EDR, YDMF_PHYS%OUT%VISICLD, YDMF_PHYS%OUT%VISIHYD,             &
& YDMF_PHYS%OUT%MXCLWC, YDMF_PHYS%OUT%MOCON)

IF (LMDUST) THEN
  ZDIFEXT (:,:,:) = 0.0_JPRB
ENDIF

CALL APLPAR_INIT ( YGFL, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, NTSSG, YSP_SBD%NLEVS, &
& YDMF_PHYS_SURF%GSD_VF%PVEG, ZMSC_FRMQ, ZDSA_CPS, ZDSA_LHS, ZDSA_RS, ZMSC_LH, ZMSC_LSCPE, ZFLU_QSAT,              &
& ZMSC_QW, ZMSC_TW, ZFLU_CD, ZFLU_CDN, ZFLU_CH, ZDSA_C1, ZDSA_C2, ZFLU_EMIS, ZFLU_FEVI, ZPFL_FTKE,                 &
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

TSPHY = MAX(PDTPHY,1.0_JPRB)

! CALL PARAMETERISATIONS

DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
  YDMF_PHYS_SURF%GSD_VF%PLSM(JROF)=REAL(NINT(YDMF_PHYS_SURF%GSD_VF%PLSM(JROF)),JPRB)
ENDDO

IF (LTWOTL) THEN
  IF (LAJUCV) THEN
    DO JLEV=1,YDCPG_DIM%KFLEVG
      DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZADJ_TAUX(JROF,JLEV)=YDVARS%T%T0(JROF,JLEV)
      ENDDO
    ENDDO
    CALL ACAJUCV(YDMODEL%YRML_PHY_MF%YRPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,       &
    & NTPLUI, YDCPG_DIM%KFLEVG, NTAJUC, YDCPG_PHY0%PREHYD, YDCPG_PHY0%XYB%ALPH, YDCPG_PHY0%XYB%DELP, &
    & YDCPG_PHY0%XYB%LNPR, YDVARS%T%T0)
    DO JLEV=1,YDCPG_DIM%KFLEVG
      DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZADJ_DTAJU(JROF,JLEV)=YDVARS%T%T0(JROF,JLEV)-ZADJ_TAUX(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
ELSE
   ! IF (LAJUCV) THEN
   !   missing code under LAJUCV for leap-frog schemes.
   ! ENDIF
ENDIF


CGMIXLEN=>YDMODEL%YRML_PHY_MF%YRPHY%CGMIXLEN

!
!-------------------------------------------------
! Check magnitude of model variables.
!-------------------------------------------------
!
IF(LGCHECKMV) CALL CHECKMV(YDRIP, YDPHY0, YDPHY2, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, &
              & YDCPG_DIM%KSTEP, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,               &
              & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDVARS%GEOMETRY%GELAM%T0,          &
              & YDVARS%GEOMETRY%GEMU%T0, ZRDG_MU0, YDMF_PHYS_SURF%GSD_VF%PLSM, YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%Q,  &
              & YLMF_PHYS_BASE_STATE%YGSP_RR%T                      )
!     ------------------------------------------------------------------

LLREDPR=LCVCSD
ZRVMD=RV-RD
! SURFEX  and passive scalar
IF (LLXFUMSE) THEN
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

IF (LMPHYS) THEN
  IF (LOZONE) THEN

    ! L'ozone est initialise par PO3 (ozone GFL).
    ! Ozone is computed from PO3 (GFL ozone).

    ZEPSO3=1.E-11_JPRB
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZQO3(JLON,0)=1.E-9_JPRB
    ENDDO
    DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZQO3(JLON,JLEV) = MAX(ZEPSO3,YLMF_PHYS_BASE_STATE%O3(JLON,JLEV))
      ENDDO
    ENDDO

  ELSEIF(YDMODEL%YRML_PHY_G%YRDPHY%NVCLIS == 1) THEN
    ZEPSO3=1.E-11_JPRB
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZQO3(JLON,0)=1.E-9_JPRB
    ENDDO
    DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZQO3(JLON,JLEV) = YDCPG_MISC%KOZO(JLON,JLEV,1)
      ENDDO
    ENDDO

  ELSEIF ((LO3FL).AND.(NOZOCL == 1).AND.(LRAYFM)) THEN
    IF (MOD(YDCPG_DIM%KSTEP,NRADFR) == 0) THEN
      CALL RADOZCMF(YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, &
      & YDVARS%GEOMETRY%GEMU%T0, ZROZ)
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZQO3(JLON,0)=1.E-9_JPRB
      ENDDO
      DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZQO3(JLON,JLEV)=ZROZ(JLON,JLEV)*YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

  ELSE
    CALL SUOZON(YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, ZQO3, .FALSE., YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, LO3ABC, YDMF_PHYS_SURF%GSD_VC%PGROUP)
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

  IF ( LNEWSTAT ) THEN

    DO JLEV = YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG-1
      ZMAN(JLEV) = FSM_CC * TANH(FSM_DD*ZVETAH(JLEV))
      ZMAK(JLEV) = FSM_EE * ZVETAH(JLEV)**FSM_FF +&
       & FSM_GG * (1-ZVETAH(JLEV))**FSM_HH + FSM_II
    ENDDO

!     MATHEMATICAL FILTER FOR THE EDGES

    IF ( LRPROX ) THEN
      DO JLEV = YDCPG_DIM%KTDIA, YDCPG_DIM%KTDIA+3
        ZMAK(JLEV) = ZMAK(JLEV) / 2**(YDCPG_DIM%KTDIA+4-JLEV)
      ENDDO
      DO JLEV = YDCPG_DIM%KFLEVG-4, YDCPG_DIM%KFLEVG-1
        ZMAK(JLEV) = ZMAK(JLEV) / 2**(5-YDCPG_DIM%KFLEVG+JLEV)
      ENDDO
    ENDIF

  ELSE
    DO JLEV = YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG-1
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

  DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA

!     LAYER WEIGHTS

      ZPOID(JLON,JLEV)=YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZGDTI
      ZIPOI(JLON,JLEV)=YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP(JLON,JLEV)*ZGDT

!     CALCULATION OF LATENT HEATS

      ZLHV(JLON,JLEV)=FOLH(YLMF_PHYS_BASE_STATE%T(JLON,JLEV),0.0_JPRB)
      ZLHS(JLON,JLEV)=FOLH(YLMF_PHYS_BASE_STATE%T(JLON,JLEV),1.0_JPRB)

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
  IF (LMDUST.AND.(NGFL_EXT/=0)) THEN
!  SIZE OF ARRAY FOR MSE
  ILONMNH=YDCPG_DIM%KFDIA-YDCPG_DIM%KIDIA+1
  ZINVDT=1/PDTPHY
  ZINVG=1/RG
  ZZI_APHI=YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI
  ZTM=YLMF_PHYS_BASE_STATE%T
!SETUP
  ZZI_SVM=0.0_JPRB
  ZZI_PEZDIAG=0.0_JPRB
  ZZI_PABSM=101325._JPRB
  ZZI_RHODREFM=1.0_JPRB
  ZZI_RHO=1.0_JPRB
  ZZZ=0.0_JPRB
  ZAERD=0.0_JPRB
  ZP1EZDIAG=0.0_JPRB

!Initialisation des scalaires passifs pour aro_ground_param
  DO JGFL=1,NGFL_EXT
    DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
       ZSVM(JLON,JLEV,JGFL)=MAX(YLMF_PHYS_BASE_STATE%P1EXT(JLON,JLEV,JGFL),0.0_JPRB)
      ENDDO
    ENDDO
  ENDDO

  !initialisation de ZZZ
  DO JLEV = 1,YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZZZ(JLON,1,JLEV)=ZZI_APHI(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO

  !initialisation de ZDZZ
  DO JLEV = 2, YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZDZZ(JLON,1,JLEV)=ZZZ(JLON,1,JLEV-1)-ZZZ(JLON,1,JLEV)
    ENDDO
  ENDDO
  DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZDZZ(JLON,1,1)=ZZI_APHI(JLON,0)*ZINVG-ZZZ(JLON,1,1)
  ENDDO


!Initialisation de ZZI_RHODREFM
  DO JLEV = 1 , YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZZI_RHODREFM(JLON,1,JLEV)=YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON,JLEV)/&
       & (YLMF_PHYS_BASE_STATE%T(JLON,JLEV)*YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R(JLON,JLEV))
     ENDDO
  ENDDO

!Initialisation de ZZI_PABSM
  DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    DO JLEV = 1 , YDCPG_DIM%KFLEVG
      ZZI_PABSM(JLON,1,JLEV)=YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON,JLEV)
    ENDDO
  ENDDO

!Initialisation de ZZI_EXNREFM
  ZRSCP=RD/RCPD
  ZINVATM=1/RATM
  DO JLEV = 1 , YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZZI_EXNREFM(JLON,1,JLEV)=(YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON,JLEV)*ZINVATM)**(ZRSCP)
    ENDDO
  ENDDO

!Initialisation de ZZI_THM
  DO JLEV = 1 , YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZZI_THM(JLON,1,JLEV)=ZTM(JLON,JLEV)/ZZI_EXNREFM(JLON,1,JLEV)
    ENDDO
  ENDDO

!Initialisation des scalaires passifs pour aro_mnhdust (inversion des niveaux)
  DO JGFL=1,NGFL_EXT
    DO JLON=1,YDCPG_DIM%KLON
      DO JLEV=1,YDCPG_DIM%KFLEVG
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
      DO JLEV = 1, YDCPG_DIM%KFLEVG
        DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZQI(JLON,JLEV)=MAX(0.0_JPRB,YLMF_PHYS_BASE_STATE%I(JLON,JLEV))
          ZQL(JLON,JLEV)=MAX(0.0_JPRB,YLMF_PHYS_BASE_STATE%L(JLON,JLEV))
          ZQR(JLON,JLEV)=MAX(0.0_JPRB,YLMF_PHYS_BASE_STATE%R(JLON,JLEV))
          ZQS(JLON,JLEV)=MAX(0.0_JPRB,YLMF_PHYS_BASE_STATE%S(JLON,JLEV))
          ZQG(JLON,JLEV)=MAX(0.0_JPRB,YLMF_PHYS_BASE_STATE%G(JLON,JLEV))
  ! CORRECTION OF NEGATIVE ADVECTED VALUES:
  !                         VAPOUR PUT IN PFCQVNG
  !                         LIQUID,ICE PUT IN PFCQL/ING
  !                         FOR RAIN/SNOW/GRAUPEL PUT IN PFCQR/S/GNG

          ZDQI=ZQI(JLON,JLEV)-YLMF_PHYS_BASE_STATE%I(JLON,JLEV)
          ZDQL=ZQL(JLON,JLEV)-YLMF_PHYS_BASE_STATE%L(JLON,JLEV)
          ZDQR=ZQR(JLON,JLEV)-YLMF_PHYS_BASE_STATE%R(JLON,JLEV)
          ZDQS=ZQS(JLON,JLEV)-YLMF_PHYS_BASE_STATE%S(JLON,JLEV)
          ZDQG=ZQG(JLON,JLEV)-YLMF_PHYS_BASE_STATE%G(JLON,JLEV)
          ZDQC=ZDQI+ZDQL+ZDQR+ZDQS+ZDQG

          ZQV0=YLMF_PHYS_BASE_STATE%Q(JLON,JLEV)-ZIPOI(JLON,JLEV)*(0.0_JPRB- YDMF_PHYS%OUT%FCQNG(JLON,JLEV-1) &
          & -YDMF_PHYS%OUT%FCQNNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQLNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQRNG(JLON,JLEV-1) &
          & -YDMF_PHYS%OUT%FCQSNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQGNG(JLON,JLEV-1))
          ZQV(JLON,JLEV)=MAX(0.0_JPRB,ZQV0-ZDQC)
          ZDQV=MAX(0.0_JPRB,ZQV0-ZDQC)-YLMF_PHYS_BASE_STATE%Q(JLON,JLEV)

          YDMF_PHYS%OUT%FCQNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNG(JLON,JLEV-1)-ZDQV*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQNNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNNG(JLON,JLEV-1)-ZDQI*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQLNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQLNG(JLON,JLEV-1)-ZDQL*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQRNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQRNG(JLON,JLEV-1)-ZDQR*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQSNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQSNG(JLON,JLEV-1)-ZDQS*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQGNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQGNG(JLON,JLEV-1)-ZDQG*ZPOID(JLON,JLEV)
        ENDDO
      ENDDO
    ELSE
      DO JLEV = 1, YDCPG_DIM%KFLEVG
        DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZQI(JLON,JLEV)=MAX(0.0_JPRB,YLMF_PHYS_BASE_STATE%I(JLON,JLEV))
          ZQL(JLON,JLEV)=MAX(0.0_JPRB,YLMF_PHYS_BASE_STATE%L(JLON,JLEV))
          ZQR(JLON,JLEV)=MAX(0.0_JPRB,YLMF_PHYS_BASE_STATE%R(JLON,JLEV))
          ZQS(JLON,JLEV)=MAX(0.0_JPRB,YLMF_PHYS_BASE_STATE%S(JLON,JLEV))

    ! CORRECTION OF NEGATIVE ADVECTED VALUES:
    !                         VAPOUR PUT IN PFCQVNG
    !                         LIQUID,ICE PUT IN PFCQL/ING
    !                         FOR RAIN/SNOW PUT IN PFCQR/SNG

          ZDQI=ZQI(JLON,JLEV)-YLMF_PHYS_BASE_STATE%I(JLON,JLEV)
          ZDQL=ZQL(JLON,JLEV)-YLMF_PHYS_BASE_STATE%L(JLON,JLEV)
          ZDQR=ZQR(JLON,JLEV)-YLMF_PHYS_BASE_STATE%R(JLON,JLEV)
          ZDQS=ZQS(JLON,JLEV)-YLMF_PHYS_BASE_STATE%S(JLON,JLEV)
          ZDQC=ZDQI+ZDQL+ZDQR+ZDQS

          ZQV0=YLMF_PHYS_BASE_STATE%Q(JLON,JLEV)-ZIPOI(JLON,JLEV)*(0.0_JPRB- YDMF_PHYS%OUT%FCQNG(JLON,JLEV-1)&
          & -YDMF_PHYS%OUT%FCQNNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQLNG(JLON,JLEV-1)-YDMF_PHYS%OUT%FCQRNG(JLON,JLEV-1)&
          & -YDMF_PHYS%OUT%FCQSNG(JLON,JLEV-1))
          ZQV(JLON,JLEV)=MAX(0.0_JPRB,ZQV0-ZDQC)
          ZDQV=MAX(0.0_JPRB,ZQV0-ZDQC)-YLMF_PHYS_BASE_STATE%Q(JLON,JLEV)

          YDMF_PHYS%OUT%FCQNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNG(JLON,JLEV-1)-ZDQV*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQNNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQNNG(JLON,JLEV-1)-ZDQI*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQLNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQLNG(JLON,JLEV-1)-ZDQL*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQRNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQRNG(JLON,JLEV-1)-ZDQR*ZPOID(JLON,JLEV)
          YDMF_PHYS%OUT%FCQSNG(JLON,JLEV)=YDMF_PHYS%OUT%FCQSNG(JLON,JLEV-1)-ZDQS*ZPOID(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

    IF (LGPCMT) THEN
      DO JLEV = 1, YDCPG_DIM%KFLEVG
        DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
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
    DO JLEV = 1, YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZQI(JLON,JLEV)=YLMF_PHYS_BASE_STATE%I(JLON,JLEV)
        ZQL(JLON,JLEV)=YLMF_PHYS_BASE_STATE%L(JLON,JLEV)
        ZQR(JLON,JLEV)=YLMF_PHYS_BASE_STATE%R(JLON,JLEV)
        ZQS(JLON,JLEV)=YLMF_PHYS_BASE_STATE%S(JLON,JLEV)
        IF (LGRAPRO) THEN
          ZQG(JLON,JLEV)=YLMF_PHYS_BASE_STATE%G(JLON,JLEV)
        ENDIF
        ZQV(JLON,JLEV)=YLMF_PHYS_BASE_STATE%Q(JLON,JLEV)
      ENDDO
    ENDDO
  ENDIF
ELSE
  DO JLEV = 1, YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZQI(JLON,JLEV)=0.0_JPRB
      ZQL(JLON,JLEV)=0.0_JPRB
      ZQV(JLON,JLEV)=YLMF_PHYS_BASE_STATE%Q(JLON,JLEV)
    ENDDO
  ENDDO
ENDIF ! LCONDWT

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

IF((LNEBR.OR.(TRIM(CGMIXLEN) == 'TM')&
       & .OR.(TRIM(CGMIXLEN) == 'TMC')).AND.YDCPG_DIM%KSTEP == 0) THEN
!DEC$ IVDEP
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    YDMF_PHYS_SURF%GSD_VH%PPBLH(JLON)=(XMINLM+XMAXLM)*0.5_JPRB
  ENDDO
ENDIF
IF((LNEBCO.OR.LGWDC).AND.YDCPG_DIM%KSTEP == 0) THEN
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    YDMF_PHYS_SURF%GSD_VH%PTCCH(JLON)=0.0_JPRB
    YDMF_PHYS_SURF%GSD_VH%PSCCH(JLON)=0.0_JPRB
    YDMF_PHYS_SURF%GSD_VH%PBCCH(JLON)=0.0_JPRB
  ENDDO
ENDIF

IF((LNEBN.OR.LNEBR.OR.LRRGUST).AND.YDCPG_DIM%KSTEP == 0) THEN
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
IF(LCVPGY.AND.YDCPG_DIM%KSTEP == 0) THEN
  DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YLMF_PHYS_BASE_STATE%CVV(JLON,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF
IF(LPHSPSH.AND.YDCPG_DIM%KSTEP == 0) THEN
  YDMF_PHYS_SURF%GSD_VH%PSPSH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0._JPRB
ENDIF
IF(LCOEFKTKE.AND.YDCPG_DIM%KSTEP == 0) THEN
  YDMF_PHYS_SURF%GSD_VH%PQSH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YLMF_PHYS_BASE_STATE%Q(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
ENDIF
IF(LCVCSD.AND.LUDEVOL.AND.YDCPG_DIM%KSTEP==0) THEN
  YDMF_PHYS_SURF%GSD_VK%PUDGRO(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0._JPRB
ENDIF
IF(LRKCDEV.AND.YDCPG_DIM%KSTEP == 0) THEN
  DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
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
IF (L3MT) THEN
  IF (LCVPRO) THEN
    YDVARS%UNEBH%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)=MAX(0._JPRB,YDVARS%UNEBH%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:))
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
    CALL ACTQSAT ( YDPHY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTQSAT, YDCPG_DIM%KFLEVG,           &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, YLMF_PHYS_BASE_STATE%T, &
    & ZGEOSLC, ZMSC_LH, ZMSC_LSCPE, ZFLU_QSAT, ZMSC_QW, YDCPG_MISC%RH, ZMSC_TW)
  ENDIF

!*
!     ------------------------------------------------------------------
!     4.BIS. COEFFICIENTS THERMO-HYDRIQUES DU SOL
!     -------------------------------------------

  IF (LSFORCS) THEN  ! Surface forcing for 1D model MUSC
     DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZTSN(JLON)=YLMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)
        ZTN(JLON) =YLMF_PHYS_BASE_STATE%T(JLON,YDCPG_DIM%KFLEVG)
     ENDDO
     DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZRHODREFM(JLON) = YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON,YDCPG_DIM%KFLEVG)/(YLMF_PHYS_BASE_STATE%T(JLON,YDCPG_DIM%KFLEVG)*YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R(JLON,YDCPG_DIM%KFLEVG))
        ZTHETAS(JLON)   = ZTSN(JLON)*(RATM/YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD(JLON,YDCPG_DIM%KFLEVG))**RKAPPA
     ENDDO
     LLAROME = .FALSE.
     CALL SURF_IDEAL_FLUX(YDRIP, YDPHY0, YDPHYDS, LLAROME, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,                            &
     & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF(:, YDCPG_DIM%KFLEVG), ZRHODREFM, YDMF_PHYS_SURF%GSD_SFO%PGROUP,                               &
     & ZTN, ZTSN, YDMF_PHYS_SURF%GSD_VF%PLSM, YLMF_PHYS_BASE_STATE%Q(:, YDCPG_DIM%KFLEVG), YLMF_PHYS_BASE_STATE%U(:, YDCPG_DIM%KFLEVG), &
     & YLMF_PHYS_BASE_STATE%V(:, YDCPG_DIM%KFLEVG), ZTHETAS, YDMF_PHYS%OUT%FCS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1),                     &
     & ZFEV, ZFMDU, ZFMDV)
     DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        YLMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)=ZTSN(JLON)
     ENDDO
     DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZDELTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-YLMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)))
        ! To be equivalent to surfex forcing
        YDMF_PHYS%OUT%FEVL(JLON,1)=ZFEV(JLON)*(1.0_JPRB-ZDELTA)
        YDMF_PHYS%OUT%FEVN(JLON,1)=ZFEV(JLON)*ZDELTA
        YDMF_PHYS%OUT%FCLL(JLON,1)=YDMF_PHYS%OUT%FEVL(JLON,1)*ZMSC_LH(JLON,YDCPG_DIM%KFLEVG)
        YDMF_PHYS%OUT%FCLN(JLON,1)=YDMF_PHYS%OUT%FEVN(JLON,1)*ZMSC_LH(JLON,YDCPG_DIM%KFLEVG)
        ZDSA_LHS(JLON)=FOLH(YLMF_PHYS_BASE_STATE%YGSP_RR%T(JLON),0.0_JPRB)
     ENDDO
  ENDIF  ! End of surface forcing for 1D model MUSC

  IF ( .NOT.LMSE ) THEN
    IF ( LSOLV ) THEN
      LLHMT=.FALSE.
      CALL ACSOL ( YDPHY, YDPHY1, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDMF_PHYS_SURF%GSD_VV%PARG,                     &
      & YDMF_PHYS_SURF%GSD_VV%PD2, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H, YDMF_PHYS_SURF%GSD_VF%PZ0RLF,            &
      & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PIVEG, YDMF_PHYS_SURF%GSD_VV%PLAI, YLMF_PHYS_BASE_STATE%YGSP_SG%A,        &
      & YLMF_PHYS_BASE_STATE%YGSP_SG%R, YDMF_PHYS_SURF%GSD_VV%PSAB, YLMF_PHYS_BASE_STATE%YGSP_SG%F, YLMF_PHYS_BASE_STATE%YGSP_RR%T, &
      & YDMF_PHYS_SURF%GSD_VF%PVEG, YLMF_PHYS_BASE_STATE%YGSP_SB%Q, YLMF_PHYS_BASE_STATE%YGSP_SB%TL,                                &
      & YLMF_PHYS_BASE_STATE%YGSP_RR%W, YLMF_PHYS_BASE_STATE%YGSP_RR%IC, LLHMT, ZDSA_C1, ZDSA_C2,                                   &
      & ZC3, ZCG, ZCN, YDMF_PHYS%OUT%CT, ZNEIJG, ZNEIJV, ZWFC, ZWPMX, ZWSEQ, ZWSMX, ZWWILT)
    ELSE

!            INITIALISATION DE L'INERTIE THERMIQUE DU SOL.

!DEC$ IVDEP
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        YDMF_PHYS%OUT%CT(JLON)=HSOL /(&
         & 1.0_JPRB+HSOLIWR*(YLMF_PHYS_BASE_STATE%YGSP_RR%W(JLON)+YLMF_PHYS_BASE_STATE%YGSP_SB%Q(JLON,1))/(WSMX+WPMX)&
         & *EXP(-0.5_JPRB*(HSOLIT0*(YLMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)-RTT))**2))
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

    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%GZ0     (JLON) = YDCPG_GPAR%GZ0(JLON)
      YDMF_PHYS%OUT%GZ0H    (JLON) = YDCPG_GPAR%GZ0H(JLON)
      ZFLU_EMIS    (JLON) = YDCPG_GPAR%VEMIS(JLON)
      YDCPG_MISC%QS      (JLON) = YDCPG_GPAR%VQS(JLON)
      YLMF_PHYS_BASE_STATE%YGSP_RR%T      (JLON) = YDCPG_GPAR%VTS(JLON)
      ZSRAIN   (JLON) = YDCPG_GPAR%RAIN(JLON)
      ZSSNOW   (JLON) = YDCPG_GPAR%SNOW(JLON)
      ZSGROUPEL(JLON) = 0._JPRB
      ZTSN     (JLON) = YLMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)
    ENDDO

    IF ( LRAY ) THEN
      ! ACRANEB/ACRANEB2 radiation => one solar band
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        YDMF_PHYS%OUT%ALB   (JLON) = YDCPG_GPAR%ALBSCA(JLON,1)
        ZALBDIR(JLON) = YDCPG_GPAR%ALBDIR(JLON,1)
      ENDDO
    ELSE
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
    ENDIF

    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      IF (ZFLU_EMIS(JLON)==0._JPRB) THEN
        ZFLU_EMIS(JLON) = 0.99_JPRB
        YLMF_PHYS_BASE_STATE%YGSP_RR%T  (JLON) = 288.0_JPRB
        YDMF_PHYS%OUT%ALB (JLON) = 0.1_JPRB
      ENDIF
    ENDDO

  ENDIF  ! LMSE
!
! Define z0;z0h if it's necessary
!
  IF (LMSE.AND.(.NOT.LCOEFKTKE).AND.(.NOT.LCOEFK_TOMS).AND.YDCPG_DIM%KSTEP == 0) THEN
  CALL ARO_GROUND_DIAG_Z0( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP, YDCPG_DIM%KFDIA-YDCPG_DIM%KIDIA+1, YDCPG_DIM%KIDIA,        &
  & YDCPG_DIM%KFDIA, NDGUNG, NDGUXG, NDLUNG, NDLUXG, YDVARS%GEOMETRY%RINDX%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                   &
  & YDVARS%GEOMETRY%RINDY%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), LSURFEX_KFROM, YDMF_PHYS%OUT%GZ0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), &
  & YDMF_PHYS%OUT%GZ0H(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA))
  ENDIF
!*
!     ------------------------------------------------------------------
!     5.- STRUCTURE ET CHAMPS DANS LA COUCHE LIMITE DE SURFACE
!     ------------------------------------------------------------------

!        INITIALISATION DES HAUTEURS "METEO".

  IF ( LHMTO ) THEN
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZDPHIV(JLON)=RG*HVCLS
      ZDPHIT(JLON)=RG*HTCLS
    ENDDO
  ENDIF

  IF ( LVDIF.OR.LHMTO.OR.LGWD ) THEN
    LLCLS=LGWD.OR.LVDIF
    LLHMT=LHMTO
    IF (LMSE) THEN
      IF(LCOEFKTKE.AND.LCOEFKSURF) THEN

       IF (YDCPG_DIM%KSTEP == 0) THEN 
         CALL ARO_GROUND_DIAG_Z0( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP, YDCPG_DIM%KFDIA-YDCPG_DIM%KIDIA+1,                         &
         & YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, NDGUNG, NDGUXG, NDLUNG, NDLUXG, YDVARS%GEOMETRY%RINDX%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),  &
         & YDVARS%GEOMETRY%RINDY%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), LSURFEX_KFROM, YDMF_PHYS%OUT%GZ0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), &
         & YDMF_PHYS%OUT%GZ0H(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA))
        ENDIF


        CALL ACTKEZOTLS ( YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,                 &
        & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, &
        & YLMF_PHYS_BASE_STATE%T, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, YLMF_PHYS_BASE_STATE%YGSP_RR%T,              &
        & YDCPG_MISC%QS, ZFLU_CDN, ZCDNH, ZDSA_CPS, ZRTI, ZDSA_RS)
         
      ELSE
        CALL ACHMTLS ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,        &
        & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%T,          &
        & YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YLMF_PHYS_BASE_STATE%YGSP_RR%T, YDCPG_MISC%QS,               &
        & ZDPHIT, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, LLCLS, ZNBVNO, ZMRIPP, ZDSA_CPS, ZGWDCS,                      &
        & ZDSA_LHS, ZPCLS, ZFLU_CD, ZFLU_CDN)
!       Computation of ZRTI
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZDPHI(JLON)=YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF(JLON,YDCPG_DIM%KFLEVG)-YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,YDCPG_DIM%KFLEVG)
          ZPRS(JLON)=RD+ZRVMD*YDCPG_MISC%QS(JLON)
          ZRTI(JLON)=2.0_JPRB/(YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R(JLON,YDCPG_DIM%KFLEVG)*YLMF_PHYS_BASE_STATE%T(JLON,YDCPG_DIM%KFLEVG)+RKAPPA*ZDPHI(JLON)&
           & +ZPRS(JLON)*YLMF_PHYS_BASE_STATE%YGSP_RR%T(JLON))
        ENDDO
      ENDIF
    ELSE
      IF (LCOEFKSURF) THEN
        CALL ACTKEHMT ( YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YDSURF%YSD_VVD%NUMFLDS>=8.AND.LSOLV, &
        & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,           &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, &
        & YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, ZPFL_FPLSH,                                    &
        & ZPFL_FPLCH, ZDPHIT, ZDPHIV, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H, YDMF_PHYS_SURF%GSD_VF%PZ0RLF,      &
        & YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VF%PLSM, ZNEIJG, ZNEIJV, YLMF_PHYS_BASE_STATE%YGSP_SG%F,                 &
        & YLMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS_SURF%GSD_VF%PVEG, ZWFC, YLMF_PHYS_BASE_STATE%YGSP_RR%W,                      &
        & YLMF_PHYS_BASE_STATE%YGSP_RR%IC, LLCLS, LLHMT, ZNBVNO, ZMRIPP, ZFLU_CD, ZFLU_CDN, ZCDROV,                              &
        & ZFLU_CH, ZCHROV, ZDSA_CPS, ZDQSTS, ZGWDCS, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, ZHQ, ZHU,                            &
        & ZFLU_NEIJ, YDMF_PHYS%OUT%QCLS, YDCPG_MISC%QS, ZFLU_QSATS, YDMF_PHYS%OUT%RHCLS, ZDSA_RS,                                &
        & ZRTI, ZSTAB, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS%OUT%NUCLS,                          &
        & YDMF_PHYS%OUT%NVCLS, ZPCLS, ZFLU_VEG, ZXDROV, ZXHROV, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST                           &
        &                                                                                                                     )
      ELSE
        CALL ACHMT ( YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YDSURF%YSD_VVD%NUMFLDS>=8.AND.LSOLV,    &
        & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,           &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, &
        & YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, ZPFL_FPLSH,                                    &
        & ZPFL_FPLCH, ZDPHIT, ZDPHIV, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H, YDMF_PHYS_SURF%GSD_VF%PZ0RLF,      &
        & YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VF%PLSM, ZNEIJG, ZNEIJV, YLMF_PHYS_BASE_STATE%YGSP_SG%F,                 &
        & YLMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS_SURF%GSD_VF%PVEG, ZWFC, YLMF_PHYS_BASE_STATE%YGSP_RR%W,                      &
        & YLMF_PHYS_BASE_STATE%YGSP_RR%IC, LLCLS, LLHMT, ZNBVNO, ZMRIPP, ZFLU_CD, ZFLU_CDN, ZCDROV,                              &
        & ZFLU_CH, ZCHROV, ZDSA_CPS, ZDQSTS, ZGWDCS, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, ZHQ, ZHU,                            &
        & ZFLU_NEIJ, YDMF_PHYS%OUT%QCLS, YDCPG_MISC%QS, ZFLU_QSATS, YDMF_PHYS%OUT%RHCLS, ZDSA_RS,                                &
        & ZRTI, ZSTAB, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YDMF_PHYS%OUT%NUCLS,                          &
        & YDMF_PHYS%OUT%NVCLS, ZPCLS, ZFLU_VEG, ZXDROV, ZXHROV, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST                           &
        &                                                                                                                     )
      ENDIF
    ENDIF

    IF (LPTKE) THEN
      YLMF_PHYS_BASE_STATE%TKE(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG) = MAX(YLMF_PHYS_BASE_STATE%TKE(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG),ETKE_MIN)
    ENDIF
    IF (LCOEFK_PTTE) THEN
      YDVARS%TTE%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG) = MAX(YDVARS%TTE%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG),ETKE_MIN)
    ENDIF

    IF(LCOEFKTKE) THEN
      ZCP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG) = RCPD*(1.0_JPRB+(RCPV/RCPD-1.0_JPRB)*(&
        & ZQV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)+ZQI(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)+&
        & ZQL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)))
    ELSE
      ZCP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG) = YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
    ENDIF

    IF(LCOEFK_RIS .AND. LCOEFKTKE) THEN
      !  computation of Ri*,Ri** for mixing lenth computation
      CALL ACMRISS ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTCOEF, YDCPG_DIM%KFLEVG, &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,  &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, ZQL, ZQI,                   &
      & ZFLU_QSAT, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%U,               &
      & YLMF_PHYS_BASE_STATE%V, ZMSC_LSCPE, YDMF_PHYS%OUT%GZ0, ZMN2PP, ZMRIPP)
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
      DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZTHETAV(JLON,JLEV)=YLMF_PHYS_BASE_STATE%T(JLON,JLEV)*(1.0_JPRB+RETV*ZQV(JLON,JLEV))&
           & *(RATM/YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON,JLEV))**RKAPPA  
        ENDDO
      ENDDO
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZTHETAVS(JLON)=YLMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)*(1.0_JPRB+RETV*YDCPG_MISC%QS(JLON))&
         & *(RATM/YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD(JLON,YDCPG_DIM%KFLEVG))**RKAPPA  
      ENDDO
      CALL ACCLPH ( YDPHY0, YDPHY2, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA,               &
      & YDCPG_DIM%KFLEVG, ZTHETAV, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,            &
      & YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, ZTHETAVS, IMOC_CLPH, YDMF_PHYS%OUT%CLPH, YDMF_PHYS%OUT%VEIN, &
      & ZUGST, ZVGST)
      ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS%OUT%CLPH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      IF (.NOT.LRAFTUR) THEN
        YDMF_PHYS%OUT%UGST(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=ZUGST(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
        YDMF_PHYS%OUT%VGST(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=ZVGST(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      ENDIF
    ELSE
      ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS_SURF%GSD_VH%PPBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
    ENDIF

  ENDIF ! end of LVDIF or LHMTO or LGWD

  IF ( (LVDIF.OR.LGWD).AND.(.NOT.(LNEBR.OR.LECT)) ) THEN

    IF(TRIM(CGMIXLEN) == 'Z') THEN
      !-------------------------------------------------
      ! "z dependent" mixing length.
      !-------------------------------------------------
      ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=1.0_JPRB/UHDIFV
      CALL ACMIXLENZ ( YDPHY, YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, 1, YDCPG_DIM%KFLEVG,     &
      & .FALSE., YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, ZBLH, YDMF_PHYS%OUT%GZ0, &
      & YDMF_PHYS%OUT%GZ0H, ZLMU, ZLMT)

    ELSEIF((TRIM(CGMIXLEN) == 'TMC').OR.(TRIM(CGMIXLEN) == 'AYC')) THEN
      !     Cubique du climat
      ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=MIN(XMAXLM,MAX(XMINLM,ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)))
      CALL ACMIXLENTM ( YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, ZBLH, ZLMU, ZLMT)

    ELSEIF(TRIM(CGMIXLEN) == 'TM') THEN
      !     Ancienne formulation pour Lm
      ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=MIN(XMAXLM,MAX(XMINLM,ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)))*XKLM
      CALL ACMIXLENZ ( YDPHY, YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, 1, YDCPG_DIM%KFLEVG,     &
      & .FALSE., YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, ZBLH, YDMF_PHYS%OUT%GZ0, &
      & YDMF_PHYS%OUT%GZ0H, ZLMU, ZLMT)

      !     Cubique du climat pour Lh
      ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)/XKLM
      CALL ACMIXLENTM ( YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS%OUT%GZ0, YDMF_PHYS%OUT%GZ0H, ZBLH, ZZLMT, ZLMT                                   &
      &                                                                                                                                                     )

    ELSEIF(TRIM(CGMIXLEN) == 'AY') THEN
      !     new Ayotte-Tudor ZBLH & mixing length
      CALL ACMIXLENZ ( YDPHY, YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, 1, YDCPG_DIM%KFLEVG,     &
      & .FALSE., YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, ZBLH, YDMF_PHYS%OUT%GZ0, &
      & YDMF_PHYS%OUT%GZ0H, ZLMU, ZLMT)

    ELSEIF((CGMIXLEN(1:2) == 'EL').AND.LPTKE) THEN
      !     e-type mixing length converted to Prandtl type
      CALL ACMIXLENZ(YDPHY, YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, 1, YDCPG_DIM%KFLEVG,      &
      & .TRUE., YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, ZBLH, YDMF_PHYS%OUT%GZ0, &
      & YDMF_PHYS%OUT%GZ0H, ZLMU, ZLMT)
      ZLMU2(:,:)=ZLMU(:,:)
      ZLMT2(:,:)=ZLMT(:,:)
      IF     (CGMIXLEN == 'EL0') THEN
        IMLTYPE=0
        ! to have identical mixing length like in pTKE
        ! e-type mixing length converted to Prandtl type
        CALL ACMIXLENZ(YDPHY, YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, 1, YDCPG_DIM%KFLEVG,       &
        & .FALSE., YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, ZBLH, YDMF_PHYS%OUT%GZ0, &
        & YDMF_PHYS%OUT%GZ0H, ZLMU, ZLMT)
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
        CALL ACMIXELEN(YGFL, YDPHY, YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA,                            &
        & YDCPG_DIM%KFLEVG, YDCPG_DIM%KSTEP, IMLTYPE, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,              &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%T,                            &
        & ZQV, ZQL, ZQI, YLMF_PHYS_BASE_STATE%R, YLMF_PHYS_BASE_STATE%S, YLMF_PHYS_BASE_STATE%TKE, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,   &
        & YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, &
        & ZMN2PP, ZFMGST, ZPFL_FPLSH, ZPFL_FPLCH, YDMF_PHYS%OUT%GZ0, YLMF_PHYS_BASE_STATE%YGSP_RR%T,                                      &
        & ZFLU_CDN, ZBLH, ZLMU, ZLMT, YDVARS%MXL%T0, ZLML, ZLMLTILD, ZRRCOR, LLMAF)
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
        CALL ACMRIS ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTCOEF,                                     &
        & YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,           &
        & ZMSC_LSCPE, ZQV, ZQL, ZQI, ZFLU_QSAT, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%T,                             &
        & YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, ZLMU, ZLMT, YDMF_PHYS%OUT%GZ0, ZMRIPP)
      ENDIF

      CALL ACMRIP(YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTCOEF, YDCPG_DIM%KFLEVG,                    &
      & YDCPG_DIM%KSTEP, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, ZQV, ZQL, ZQI, ZCP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH,                            &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZFLU_QSAT, ZMSC_QW, ZMSC_TW, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,                     &
      & YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YDVARS%FQTUR%T0, YDVARS%FSTUR%T0,                     &
      & YDVARS%SHTUR%T0, YLMF_PHYS_BASE_STATE%TKE, YDVARS%TTE%T0, YLMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS_SURF%GSD_VH%PQSH,         &
      & ZDSA_RS, ZDSA_CPS, ZRTI, YDMF_PHYS%OUT%GZ0, LLCLS, ZMRIPP, ZMRIFPP, ZBNEBCVPP, ZBNEBQ,                                        &
      & ZNBVNO, ZFMTKE, ZFTTKE, ZF_EPS, ZFUN_TTE, ZAUTKE, ZATTKE, ZFHORM, ZFHORH, ZTH_FUN, ZWW_FUN,                                   &
      & ZMRIMC, ZMRICTERM, ZMN2PP, ZMN2_ES, ZMN2_EQ, ZMN2_DS, ZMN2_DQ, ZFMGST)

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
      CALL ACMIXELEN(YGFL, YDPHY, YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA,                            &
      & YDCPG_DIM%KFLEVG, YDCPG_DIM%KSTEP, IMLTYPE, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,              &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%T,                            &
      & ZQV, ZQL, ZQI, YLMF_PHYS_BASE_STATE%R, YLMF_PHYS_BASE_STATE%S, YLMF_PHYS_BASE_STATE%TKE, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,   &
      & YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, &
      & ZMN2PP, ZFMGST, ZPFL_FPLSH, ZPFL_FPLCH, YDMF_PHYS%OUT%GZ0, YLMF_PHYS_BASE_STATE%YGSP_RR%T,                                      &
      & ZFLU_CDN, ZBLH, ZLMU, ZLMT, YDVARS%MXL%T0, ZLML, ZLMLTILD, ZRRCOR, LLMAF)
    ENDIF

  ENDIF ! (LVDIF or LGWD) and( not(LNEBR or LECT))

  IF ( LVDIF.OR.LHMTO.OR.LGWD ) THEN

    IF (LFLUSO.AND.(.NOT.LMSE)) THEN
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZCRTI(JLON) = 1.0_JPRB/(YLMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)*ZDSA_RS(JLON))
      ENDDO
      CALL ACFLUSO ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,          &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,   &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%Q, YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%U, &
      & YLMF_PHYS_BASE_STATE%V, ZDPHIT, ZDPHIV, YDMF_PHYS%OUT%GZ0, YDMF_PHYS_SURF%GSD_VF%PLSM,                         &
      & ZFLU_QSATS, ZCRTI, YLMF_PHYS_BASE_STATE%YGSP_RR%T, LLHMT, ZFLU_CD, ZFLU_CDN, ZCDROV, ZCE,                      &
      & ZCEROV, ZFLU_CH, ZCHROV, YDMF_PHYS%OUT%QCLS, YDMF_PHYS%OUT%RHCLS, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%UCLS,      &
      & YDMF_PHYS%OUT%VCLS, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST)
    ELSE
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZCE   (JLON) = ZFLU_CH   (JLON)
        ZCEROV(JLON) = ZCHROV(JLON)
      ENDDO
    ENDIF

  ENDIF ! LVDIF or LHMTO or LGWD

  IF (LRAFTKE) THEN
    YDMF_PHYS%OUT%CAPE(:)=0._JPRB
    ZDCAPE(:)=0._JPRB
    CALL ACCLDIA(YDXFU, YDPHY, YDPHY2, YDTOPH, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, &
    & YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YDMF_PHYS%OUT%CAPE,  &
    & ZDCAPE, YLMF_PHYS_BASE_STATE%TKE, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDVARS%GEOMETRY%OROG%T0, YDMF_PHYS%OUT%UGST,       &
    & YDMF_PHYS%OUT%VGST, ZBLH, IMOC_CLPH)
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
      CALL ACTKECOEFK ( YDPHY0, YDPHY2, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTCOEF, YDCPG_DIM%KFLEVG,  &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, &
      & ZFMTKE, ZFTTKE, ZAUTKE, ZATTKE, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%T,                 &
      & YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, ZLMU, ZLMT, YDMF_PHYS%OUT%GZ0, ZKTROV,                       &
      & ZKUROV, ZKNROV, ZXTROV, ZXUROV, ZXPTKEROV)
    ELSE
      CALL ACCOEFK ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTCOEF, YDCPG_DIM%KFLEVG,        &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,         &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
      & ZMSC_LSCPE, ZQV, ZQL, ZQI, ZFLU_QSAT, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%T,                   &
      & YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, ZPFL_FPLSH, ZPFL_FPLCH, ZLMU, ZLMT, YDMF_PHYS%OUT%GZ0,               &
      & YDMF_PHYS%OUT%GZ0H, ZBLH, ZKTROV, ZKUROV, ZKNROV, ZNBVNO, ZXTROV, ZXUROV, ZXPTKEROV)
    ENDIF
  ENDIF

  IF (LCVPPKF) THEN
    CALL ACVPPKF( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTCVIM, YDCPG_DIM%KFLEVG,        &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,  &
    & YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%T, ZQV, ZQL, ZQI, YLMF_PHYS_BASE_STATE%U,                 &
    & YLMF_PHYS_BASE_STATE%V, YDCPG_DYN0%CTY%VVEL(:, 1:), YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%TKE, &
    & ZDIFCVPPQ, ZDIFCVPPS, ZCONDCVPPL, ZCONDCVPPI, ZPRODTH_CVPP, INLAB_CVPP, ZQLI_CVPP, ZNEB_CVPP,                       &
    & INND)
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
      
    IF (YDCPG_DIM%KSTEP == 0) YDMF_PHYS_SURF%GSD_SFL%PGROUP(:,:) = 0.0_JPRB
    CALL ARP_SHALLOW_MF( YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG,                   &
    & ZIMPL, TSPHY, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, &
    & YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, CMF_UPDRAFT, CMF_CLOUD, LMIXUV, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V,     &
    & YLMF_PHYS_BASE_STATE%T, ZQV, ZQL, ZQI, ZQR, ZQS, YLMF_PHYS_BASE_STATE%TKE, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF,         &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZEDMFQ, ZEDMFS, ZEDMFU, ZEDMFV, YDMF_PHYS_SURF%GSD_SFL%PGROUP(:, 1),              &
    & YDMF_PHYS_SURF%GSD_SFL%PGROUP(:, 2), ZPRODTH_CVPP, ZQLI_CVPP, ZNEB_CVPP, INLAB_CVPP, ZMF_UP)

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
      DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZNLABCVP(JLON,JLEV) = ZNLABCVP(JLON,JLEV)&
         & *MAX(0.0_JPRB,SIGN(1.0_JPRB,ZPFL_FPLCH(JLON,JLEV)-ZPFL_FPLCH(JLON,JLEV-1)-ZEPS))
          ZNLAB(JLON,JLEV) = REAL(INLAB_CVPP(JLON,JLEV),JPRB)
        ENDDO
      ENDDO
    ENDIF
    DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZNLAB(JLON,JLEV) = REAL(INLAB_CVPP(JLON,JLEV),JPRB)
      ENDDO
    ENDDO

    CALL ACTKE ( YDLDDH, YDMODEL%YRML_DIAG%YRMDDH, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA,                      &
    & YDCPG_DIM%KLON, NTCOEF, NTCOET, YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,     &
    & YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V,             &
    & ZQV, ZQIC, ZQLC, ZMSC_LSCPE, ZFLU_CD, ZFLU_CH, YDMF_PHYS%OUT%GZ0, YLMF_PHYS_BASE_STATE%YGSP_RR%T,                        &
    & YDCPG_MISC%QS, YDCPG_MISC%QICE, YDCPG_MISC%QLI, YLMF_PHYS_BASE_STATE%TKE, ZPRODTH_CVPP, ZNLAB,                           &
    & ZNLABCVP, ZKTROV, ZKQROV, ZKQLROV, ZKUROV, ZXTROV, ZXUROV, ZNBVNO, ZNEBS, ZQLIS, ZNEBS0,                                 &
    & ZQLIS0, ZCOEFN, ZPFL_FTKE, ZPFL_FTKEI, ZTKE1, ZTPRDY, YDMF_PHYS%OUT%EDR, YDDDH)
    YDMF_PHYS%OUT%CLPH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=MIN(XMAXLM,MAX(XMINLM,YDMF_PHYS%OUT%CLPH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)))
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

  CALL ARO_MNHDUST (1, ILONMNH, YDCPG_DIM%KFLEVG, NGFL_EXT, PDTPHY, ZZI_SVM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :, 1:NGFL_EXT),                &
  & ZZZ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :), ZDZZ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :), ZZI_PABSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :), &
  & ZZI_THM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :), ZZI_RHODREFM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :),                                       &
  & NSWB_MNH, YDCPG_DIM%KSTEP+1, ZZI_SVM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :, 1:NGFL_EXT), ZPIZA_DST(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :), &
  & ZCGA_DST(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :), ZTAUREL_DST(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :),                                       &
  & ZAERD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :), NGFL_EZDIAG, ZZI_PEZDIAG(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :)                                 &
  &                                                                                                                                                                                           )

  ZP1EZDIAG(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:,:)=ZZI_PEZDIAG(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:,:)

! return to aladin environment (inversion des niveaux)
   DO JGFL=1,NGFL_EXT
     DO JLON=1,YDCPG_DIM%KLON
       DO JLEV=1,YDCPG_DIM%KFLEVG
         ZSVM(JLON,JLEV,JGFL)=ZZI_SVM(JLON,1,JLEV,JGFL)
       ENDDO
     ENDDO
   ENDDO
  ENDIF

!      7.1 Albedo et emissivite en presence de neige
!          Albedo and emissivity with snow

  IF (.NOT.LMSE) THEN
!DEC$ IVDEP
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
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
         & MAX(YDMF_PHYS_SURF%GSD_VF%PALBSF(JLON),YLMF_PHYS_BASE_STATE%YGSP_SG%A(JLON,1))&
         & + YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)*ZNEIJV(JLON) *&
         & MAX(ZALBV,YLMF_PHYS_BASE_STATE%YGSP_SG%A(JLON,1))&
         & + YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)*(1.0_JPRB-ZNEIJV(JLON)) * ZALBV
        ZFLU_EMIS(JLON)= (1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PVEG(JLON))*(1.0_JPRB-ZNEIJG(JLON)) *&
         & YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)&
         & + (1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PVEG(JLON))*ZNEIJG(JLON) * EMCRIN&
         & + YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)*ZNEIJV(JLON) * EMCRIN&
         & + YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)*(1.0_JPRB-ZNEIJV(JLON)) * YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)
      ELSE
        IF (LVGSN) THEN
          IF (LZ0HSREL.AND.LCOEFKSURF) THEN
            ! new treatment, PNEIJ is gridbox snow fraction
            YDMF_PHYS%OUT%ALB(JLON)=(1.0_JPRB-ZFLU_VEG(JLON)-ZFLU_NEIJ(JLON))*YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)+ &
             & ZFLU_VEG(JLON)*YDMF_PHYS_SURF%GSD_VV%PALV(JLON)+ZFLU_NEIJ(JLON)*YLMF_PHYS_BASE_STATE%YGSP_SG%A(JLON,1)
          ELSE
            ! old treatment, PNEIJ is snow fraction for bare ground
            YDMF_PHYS%OUT%ALB(JLON)=YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)-ZFLU_NEIJ(JLON)*(YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)- &
             & YLMF_PHYS_BASE_STATE%YGSP_SG%A(JLON,1))+(ZFLU_NEIJ(JLON)-ZNEIJV(JLON))*     &
             & YDMF_PHYS_SURF%GSD_VF%PVEG(JLON)*(YDMF_PHYS_SURF%GSD_VV%PALV(JLON)-YLMF_PHYS_BASE_STATE%YGSP_SG%A(JLON,1))
          ENDIF

          YDMF_PHYS%OUT%ALB(JLON)=MIN(ABS(YDMF_PHYS_SURF%GSD_VV%PIVEG(JLON)-2._JPRB),1.0_JPRB) * YDMF_PHYS%OUT%ALB(JLON) +(&
           & 1.0_JPRB-MIN(ABS(YDMF_PHYS_SURF%GSD_VV%PIVEG(JLON)-2._JPRB),1.0_JPRB))&
           & * MAX(ALCRIN,YDMF_PHYS%OUT%ALB(JLON))
          YDMF_PHYS_SURF%GSP_SG%PT_T1(JLON,1)=YDMF_PHYS%OUT%ALB(JLON)

          ZFLU_EMIS(JLON)=YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)-ZFLU_NEIJ(JLON)*(YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)-EMCRIN)

        ELSE
          YDMF_PHYS%OUT%ALB(JLON)=YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)-ZFLU_NEIJ(JLON)*(YDMF_PHYS_SURF%GSD_VF%PALBF(JLON)&
           & -MAX(YDMF_PHYS_SURF%GSD_VF%PALBF(JLON),ALCRIN))
          ZFLU_EMIS(JLON)=YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)-ZFLU_NEIJ(JLON)*(YDMF_PHYS_SURF%GSD_VF%PEMISF(JLON)-EMCRIN)
        ENDIF
      ENDIF
    ENDDO

    IF (LRAYFM) THEN
      ! diffuse and direct (parallel) albedo in NSW solar intervals
      IF (LALBMERCLIM) THEN
        DO JSG=1,NSW
          DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            ZALBD(JLON,JSG)=YDMF_PHYS%OUT%ALB(JLON)
            ZALBPMER=(1.0_JPRB+&
             & 0.5_JPRB*ZRDG_MU0M(JLON)*(1.0_JPRB/YDMF_PHYS%OUT%ALB(JLON)-1.0_JPRB))/   (&
             & 1.0_JPRB+ZRDG_MU0M(JLON)*(1.0_JPRB/YDMF_PHYS%OUT%ALB(JLON)-1.0_JPRB))**2
            ZALBP(JLON,JSG)=YDMF_PHYS%OUT%ALB(JLON)*          YDMF_PHYS_SURF%GSD_VF%PLSM(JLON)+&
                          & ZALBPMER  *(1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PLSM(JLON))
          ENDDO
        ENDDO
      ELSE
!DEC$ IVDEP
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZMERL(JLON)=(1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PLSM(JLON))*(1.0_JPRB&
            & -MAX(0.0_JPRB,SIGN(1.0_JPRB,TMERGL-YLMF_PHYS_BASE_STATE%YGSP_RR%T(JLON))))
          ZFLU_EMIS(JLON)=ZFLU_EMIS(JLON)*YDMF_PHYS_SURF%GSD_VF%PLSM(JLON)+ZMERL(JLON)*EMMMER&
            & +(1.0_JPRB-YDMF_PHYS_SURF%GSD_VF%PLSM(JLON))*(1.0_JPRB-ZMERL(JLON))*EMMGLA
        ENDDO
        DO JSG=1,NSW
          DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            ZALBD(JLON,JSG)=(1.0_JPRB-ZMERL(JLON))*YDMF_PHYS%OUT%ALB(JLON)+&
                                    & ZMERL(JLON) *ALBMED
            ZALBP(JLON,JSG)=(1.0_JPRB-ZMERL(JLON))*YDMF_PHYS%OUT%ALB(JLON)+&
             & ZMERL(JLON) *&
             & MAX(0.037_JPRB/(1.1_JPRB*ZRDG_MU0(JLON)**1.4_JPRB+0.15_JPRB),ZEPS0)
          ENDDO
        ENDDO
      ENDIF
    ELSEIF (LRAY) THEN
      ! direct (parallel) albedo for ACRANEB/ACRANEB2, Geleyn's formula
      ! with given proportion of Lambertian scattering
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        IF ( YDMF_PHYS_SURF%GSD_VF%PLSM(JLON) < 0.5_JPRB .AND. YLMF_PHYS_BASE_STATE%YGSP_RR%T(JLON) >= TMERGL ) THEN
          ZLAMB=RLAMB_WATER  ! water surface (open sea)
        ELSE
          ZLAMB=RLAMB_SOLID  ! solid surface (frozen sea or land)
        ENDIF
        ZALBDIR(JLON)=ZLAMB*YDMF_PHYS%OUT%ALB(JLON)+(1._JPRB-ZLAMB)*(1._JPRB+&
         & 0.5_JPRB*ZRDG_MU0(JLON)*(1.0_JPRB/YDMF_PHYS%OUT%ALB(JLON)-1.0_JPRB))/   (&
         & 1.0_JPRB+ZRDG_MU0(JLON)*(1.0_JPRB/YDMF_PHYS%OUT%ALB(JLON)-1.0_JPRB))**2
      ENDDO
    ENDIF

  ENDIF  ! .NOT.LMSE

  ! Appel de la routine d'aerosols

  LLAERO=LAEROSEA.AND.LAEROLAN.AND.LAEROSOO.AND.LAERODES

  IF ( .not. LRAYFM) THEN
    NAER=0
    NRADFR=NSTOP+1
  ENDIF

  IF    (   (LRAYFM.AND.(MOD(YDCPG_DIM%KSTEP,NRADFR) == 0)) &
  & .OR.  ( (LRAY.OR.LRAYSP).AND.(.NOT.LRSTAER)) ) THEN

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
    IF (LAEROSUL) THEN
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZAESUL(JLON) = YDMF_PHYS_SURF%GSD_VA%PSUL(JLON)
      ENDDO
    ELSE
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZAESUL(JLON) = 0.0_JPRB
      ENDDO
    ENDIF
    IF (LAEROVOL) THEN
       DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZAEVOL(JLON) = YDMF_PHYS_SURF%GSD_VA%PVOL(JLON)
      ENDDO
    ELSE
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZAEVOL(JLON) = 0.0_JPRB
      ENDDO
    ENDIF

    IF ( ( (LRAYFM.AND.NAER /= 0) .OR.LRAY.OR.LRAYSP).AND.LLAERO )  THEN
      CALL RADAER ( YDMODEL%YRML_PHY_RAD%YREAERD, YDERAD, YDPHY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, &
      & YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF,             &
      & YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%YGSP_RR%T, ZAESEA, ZAELAN, ZAESOO, ZAEDES,                    &
      & ZAESUL, ZAEVOL, ZAER, ZAERINDS                                                     )
    ENDIF

  ELSEIF ( (LRAY.OR.LRAYSP).AND.LRSTAER ) THEN

    DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZAER(JLON,JLEV,1)=ZDAER(JLEV)
        ZAER(JLON,JLEV,2:6)=0._JPRB
      ENDDO
    ENDDO

  ENDIF ! FOR AEROSOLS

! We uses the extinction coefficient explicitely solved by ARO_MNHDUST
  IF (LMDUST.AND.(NGFL_EXT/=0)) THEN
         ZAER(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:,3) = ZAERD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
  ENDIF


!      7.2 Flux radiatifs par ciel clair (Code Geleyn)
!          Clear sky radiative fluxes    (Geleyn's scheme)

  ! separate clearsky call is kept only for old ACRANEB; for ACRANEB2
  ! duplicit calculation of gaseous transmissions is avoided
  IF (LRAY.AND.NRAY == 1.AND.YDCFU%NFRRC /= 0) THEN
    IF (MOD(YDCPG_DIM%KSTEP,YDCFU%NFRRC) == 0) THEN
      CALL ACRANEB(YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,                          &
      & NTRADI, YDCPG_DIM%KFLEVG, IJN, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,          &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDCPG_MISC%NEB, ZQV, ZQCO2, YDCPG_MISC%QICE, YDCPG_MISC%QLI,              &
      & ZQO3, YLMF_PHYS_BASE_STATE%T, YDMF_PHYS%OUT%ALB, ZALBDIR, ZFLU_EMIS, ZRDG_MU0, YDVARS%GEOMETRY%GEMU%T0,                     &
      & YDVARS%GEOMETRY%GELAM%T0, ZRDG_MU0LU, YLMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH,               &
      & ZFRSODS, YDMF_PHYS%OUT%FRSOPS, YDMF_PHYS%OUT%FRSOLU, YDMF_PHYS%OUT%FRTHDS, ZAER, ZMAK, ZMAN)
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        YDMF_PHYS%OUT%FRSOC(JLON,0)=YDMF_PHYS%OUT%FRSO(JLON,NTRADI-1,1)
        YDMF_PHYS%OUT%FRSOC(JLON,1)=YDMF_PHYS%OUT%FRSO(JLON,YDCPG_DIM%KFLEVG,1)
        YDMF_PHYS%OUT%FRTHC(JLON,0)=YDMF_PHYS%OUT%FRTH(JLON,NTRADI-1,1)
        YDMF_PHYS%OUT%FRTHC(JLON,1)=YDMF_PHYS%OUT%FRTH(JLON,YDCPG_DIM%KFLEVG,1)
      ENDDO
    ENDIF
  ENDIF

!      7.3 Nebulosite et Convection
!          Cloud cover and Convection
!      7.3.1 Shallow + Deep convection

  IF (LCVPGY) THEN
    ! Le schema de convection de J. F. Gueremy
    IF (LCONDWT) THEN
      DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZQCL(JLON,JLEV)=ZQL(JLON,JLEV)
          ZQCI(JLON,JLEV)=ZQI(JLON,JLEV)
        ENDDO
      ENDDO
    ELSE
      DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZQCL(JLON,JLEV)=0.0_JPRB
          ZQCI(JLON,JLEV)=0.0_JPRB
        ENDDO
      ENDDO
    ENDIF
    CALL ACCVIMPGY ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTCVIM, YDCPG_DIM%KFLEVG,                   &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,                 &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,              &
    & ZMSC_LH, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZQV, ZQCI, ZQCL, ZQLIS, ZFLU_QSAT, ZMSC_QW, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YLMF_PHYS_BASE_STATE%T, ZMSC_TW, YLMF_PHYS_BASE_STATE%U,                                 &
    & YLMF_PHYS_BASE_STATE%V, ZDSA_CPS, YDVARS%GEOMETRY%GM%T0, YLMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS%OUT%DIFCQ,                               &
    & YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN,                        &
    & YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN,                                         &
    & ZFHMLTSC, ZFHEVPPC, ZFPEVPPC, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, ZNEBC, ZQLIC, YDMF_PHYS%OUT%STRCU,                        &
    & YDMF_PHYS%OUT%STRCV, ICIS, INLAB, INND, YLMF_PHYS_BASE_STATE%CVV)

    DO JLEV = YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
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
      DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZNEBC0(JLON,JLEV)=ZNEBC(JLON,JLEV)
          ZQLI_CVP(JLON,JLEV)=ZQLIC(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

    ! Annulation possible des flux convectifs pour les eaux condensees.
    IF (.FALSE.) THEN
      DO JLEV=0,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          YDMF_PHYS%OUT%DIFCQL(JLON,JLEV)=0.0_JPRB
          YDMF_PHYS%OUT%DIFCQN(JLON,JLEV)=0.0_JPRB
        ENDDO
      ENDDO
    ENDIF
  ENDIF !  (LCVPGY)

  IF ( LCONDWT.AND.(.NOT.LNEBECT)) THEN

    IF(LCVPRO.AND.LNEBCV) THEN
! convective cloudiness in case we need protection of convective cloud water.
      DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZNEBCH(JLON,JLEV)=ZUNEBH(JLON,JLEV)
        ENDDO
      ENDDO
      IF (LCVCSD) THEN
        DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
          DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            ZNEBC0(JLON,JLEV)=MAX(ZEPSNEB,&
                            & MIN(1._JPRB-ZEPSNEB,ZUNEBH(JLON,JLEV)))
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    CALL ACNEBCOND ( YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,                            &
    & NTPLUI, YDCPG_DIM%KFLEVG, LLREDPR, ZHUC, ZVETAF, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,       &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDCPG_MISC%RH, ZBLH, ZQV, ZQI, ZQL, ZMSC_QW, YLMF_PHYS_BASE_STATE%T,            &
    & ZNEBCH, YDVARS%GEOMETRY%GM%T0, YLMF_PHYS_BASE_STATE%YGSP_RR%T, ZQLIS, ZNEBS, ZRHCRI, ZRH, ZQSATS,                                 &
    & ZICEFR1, ZQLIS0, ZNEBS0)
     
    IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA, 'ZRHCRI', ZRHCRI, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG) 
    IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA, 'ZQLIS0', ZQLIS0, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG)
    IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA, 'ZQLIS', ZQLIS, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG)

    IF(LRKCDEV) THEN
! Rash-Kristiansson cloud water scheme - second part.
      DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
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
      CALL ACPCMT(YDGEM, YDGEOMETRY%YRDIM, YDGEOMETRY%YREGEO, YDLDDH, YDMODEL%YRML_DIAG%YRMDDH,                                     &
      & YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTCVIM, YDCPG_DIM%KFLEVG,                     &
      & INBTRA, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,      &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDCPG_DYN0%CTY%VVEL(:, 1:),                    &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZMSC_LH, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, YLMF_PHYS_BASE_STATE%Q, ZQI, ZQL, YLMF_PHYS_BASE_STATE%R,                           &
      & YLMF_PHYS_BASE_STATE%S, ZQLIS, ZFLU_QSAT, ZMSC_QW, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,                                 &
      & ZRDG_CVGQ, YLMF_PHYS_BASE_STATE%T, ZMSC_TW, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YLMF_PHYS_BASE_STATE%TKE,       &
      & ZTRA(:, :, INBTRA_DEP:INBTRA), ZSMOOTRAC(1:INBTRA), ZQLC, ZQIC, ZQRC, ZQSC, YDVARS%GEOMETRY%GM%T0, YLMF_PHYS_BASE_STATE%YGSP_RR%T,    &
      & YDCPG_MISC%QS, ZFLU_NEIJ, YDMF_PHYS_SURF%GSD_VF%PLSM, ZVETAF, ZEDMFQ, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCQN,           &
      & YDMF_PHYS%OUT%DIFCQLC, YDMF_PHYS%OUT%DIFCQIC, ZEDMFS, ZDIFCVTH, ZTMPPRODTH, YDMF_PHYS%OUT%FCCQL,                            &
      & YDMF_PHYS%OUT%FCCQN, ZEDMFU, ZEDMFV, ZSTRCTRA(:, :, INBTRA_DEP:INBTRA), YDMF_PHYS%OUT%FIMCC,                                &
      & YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN,                                     &
      & ZMF_UP, ZMU, ZMD, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC,                   &
      & YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC, YDMF_PHYS%OUT%FCNEGQIC, YDMF_PHYS%OUT%FCNEGQRC,         &
      & YDMF_PHYS%OUT%FCNEGQSC, INLAB, ZNEBC0, ZQLI_CVP, ZTU, ZQU, ZQC_DET_PCMT, ZCSGC, ZENTCH, INND,                               &
      & YDMF_PHYS%OUT%CAPE, ZAIPCMT, ZALF_CAPE, ZALF_CVGQ, YDVARS%UAL%T0, YDVARS%UOM%T0, YDVARS%DAL%T0,                             &
      & YDVARS%DOM%T0, YDDDH)
    ELSE
      CALL ACPCMT(YDGEM, YDGEOMETRY%YRDIM, YDGEOMETRY%YREGEO, YDLDDH, YDMODEL%YRML_DIAG%YRMDDH,                                     &
      & YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTCVIM, YDCPG_DIM%KFLEVG,                     &
      & INBTRA, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,      &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDCPG_DYN0%CTY%VVEL(:, 1:),                    &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZMSC_LH, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, YLMF_PHYS_BASE_STATE%Q, ZQI, ZQL, YLMF_PHYS_BASE_STATE%R,                           &
      & YLMF_PHYS_BASE_STATE%S, ZQLIS, ZFLU_QSAT, ZMSC_QW, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,                                 &
      & ZRDG_CVGQ, YLMF_PHYS_BASE_STATE%T, ZMSC_TW, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YLMF_PHYS_BASE_STATE%TKE,       &
      & ZTRA(:, :, INBTRA_DEP:INBTRA), ZSMOOTRAC(1:INBTRA), ZQLC, ZQIC, ZQRC, ZQSC, YDVARS%GEOMETRY%GM%T0, YLMF_PHYS_BASE_STATE%YGSP_RR%T,    &
      & YDCPG_MISC%QS, ZFLU_NEIJ, YDMF_PHYS_SURF%GSD_VF%PLSM, ZVETAF, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQL,                    &
      & YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQLC, YDMF_PHYS%OUT%DIFCQIC, YDMF_PHYS%OUT%DIFCS,                                    &
      & ZDIFCVTH, ZTMPPRODTH, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV,                   &
      & ZSTRCTRA(:, :, INBTRA_DEP:INBTRA), YDMF_PHYS%OUT%FIMCC, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN,                       &
      & YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, ZMF_UP, ZMU, ZMD, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN,                     &
      & YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC,             &
      & YDMF_PHYS%OUT%FCNEGQIC, YDMF_PHYS%OUT%FCNEGQRC, YDMF_PHYS%OUT%FCNEGQSC, INLAB, ZNEBC0,                                      &
      & ZQLI_CVP, ZTU, ZQU, ZQC_DET_PCMT, ZCSGC, ZENTCH, INND, YDMF_PHYS%OUT%CAPE, ZAIPCMT,                                         &
      & ZALF_CAPE, ZALF_CVGQ, YDVARS%UAL%T0, YDVARS%UOM%T0, YDVARS%DAL%T0, YDVARS%DOM%T0, YDDDH)
    ENDIF

  ENDIF

!         Appel du calcul de nebulosite.

  IF(LNEBN) THEN
    CALL ACNEBN ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTNEBU, YDCPG_DIM%KFLEVG,                    &
    & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, ZQL, ZQI,                                        &
    & ZFLU_QSAT, YLMF_PHYS_BASE_STATE%T, ZPFL_FPLCH, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, ZUNEBH, ZNEBS0, ZQLIS0, ZQLI_CVP, ZNEB_CVPP, ZQLI_CVPP,                                  &
    & ZAIPCMT, YDCPG_MISC%NEB, ZNEBC0, YDCPG_MISC%QICE, YDCPG_MISC%QLI, ZHUC, ZVETAF)
    DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
!DEC$ IVDEP
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        YDCPG_MISC%NEB(JLON,JLEV)=YDCPG_MISC%NEB(JLON,JLEV)*GAEPS
      ENDDO
    ENDDO
  ENDIF

  IF ( LNEBR ) THEN
    CALL ACNEBR ( YDERAD, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,                 &
    & NTCOEF, NTNEBU, YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,   &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, ZMSC_LSCPE, ZQV,              &
    & ZFLU_QSAT, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%U,            &
    & YLMF_PHYS_BASE_STATE%V, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,   &
    & YDMF_PHYS%OUT%GZ0, ZPFL_FPLCH, YDMF_PHYS_SURF%GSD_VH%PPBLH, YDCPG_MISC%QS, YLMF_PHYS_BASE_STATE%YGSP_RR%T, &
    & ZKTROV, ZKUROV, ZNBVNO, YDCPG_MISC%NEB, ZNEBS, YDCPG_MISC%QICE, YDCPG_MISC%QLI, ZQLIS)
  ENDIF

!         Diagnostique de nebulosite partielle.
  CALL ACNPART(YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTNEBU, YDCPG_DIM%KFLEVG,   &
  & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, &
  & ZDECRD, YDCPG_MISC%NEB, YDMF_PHYS%OUT%CLCH, YDMF_PHYS%OUT%CLCM, YDMF_PHYS%OUT%CLCL, YDCPG_MISC%CLCT,          &
  & ZCLCT_RAD, PCLCC=YDMF_PHYS%OUT%CLCC, PNEBC=ZNEBC0, PTOPC=YDMF_PHYS%OUT%CTOP)

!     7.3.5 Computation of the equivalent coefficients for simplified
!           radiation scheme

  IF ( LRAYSP .AND.(YDCPG_DIM%KSTEP == 1).AND.LRCOEF ) THEN
    DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZZNEB(JLON,JLEV)=0.0_JPRB
      ENDDO
    ENDDO

    CALL ACRADCOEF ( YDRCOEF, YDPHY3, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA, &
    & YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF,     &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZZNEB, ZQV, ZQCO2, YDCPG_MISC%QICE, YDCPG_MISC%QLI,        &
    & ZQO3, YLMF_PHYS_BASE_STATE%T, ZRDG_MMU0, YDCPG_MISC%DHSF, ZFLU_EMIS, ZAER, ZAC, ZAC_HC,            &
    & ZRDT_COR, ZRDT_RAB3C, ZRDT_RAB3N, ZRDT_RAB4C, ZRDT_RAB4N, ZRDT_RAB6C, ZRDT_RAB6N, ZRDT_RAT1C,      &
    & ZRDT_RAT1N, ZRDT_RAT2C, ZRDT_RAT2N, ZRDT_RAT3C, ZRDT_RAT3N, ZRDT_RAT4C, ZRDT_RAT4N, ZRDT_RAT5C,    &
    & ZRDT_RAT5N)
  ENDIF

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
     DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
       DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
         ZDELP(JLON,JLEV) = YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD(JLON,JLEV) - YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD(JLON,JLEV-1)
       ENDDO
     ENDDO
     ZWND  (:)      = 0._JPRB  ! not used in  ARPEGE-Climat chemistry
     ZDUMMY1 (:,:)  = 1.0E-18_JPRB
     ZGELAT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA) = ASIN(YDVARS%GEOMETRY%GEMU%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA))
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
        CALL ACRANEB(YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,                          &
        & NTRADI, YDCPG_DIM%KFLEVG, IJN, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,          &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDCPG_MISC%NEB, ZQV, ZQCO2, YDCPG_MISC%QICE, YDCPG_MISC%QLI,              &
        & ZQO3, YLMF_PHYS_BASE_STATE%T, YDMF_PHYS%OUT%ALB, ZALBDIR, ZFLU_EMIS, ZRDG_MU0, YDVARS%GEOMETRY%GEMU%T0,                     &
        & YDVARS%GEOMETRY%GELAM%T0, ZRDG_MU0LU, YLMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH,               &
        & ZFRSODS, YDMF_PHYS%OUT%FRSOPS, YDMF_PHYS%OUT%FRSOLU, YDMF_PHYS%OUT%FRTHDS, ZAER, ZMAK,                            &
        & ZMAN)

        ! update sunshine duration [s]
!DEC$ IVDEP
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          IF ( YDMF_PHYS%OUT%FRSOPS(JLON) > RSUNDUR*ZRDG_MU0(JLON) ) THEN
            YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)=YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)+TSTEP
          ENDIF
        ENDDO
      CASE(2)
        CALL ACRANEB2(YDERDI, YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,                   &
        & NTRADI, YDCPG_DIM%KFLEVG, IJN, YDCPG_DIM%KSTEP, YDCFU%NFRRC, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,                  &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,   &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDCPG_MISC%NEB, ZQV, ZQCO2, YDCPG_MISC%QICE, YDCPG_MISC%QLI,                &
        & ZQO3, YLMF_PHYS_BASE_STATE%T, YDMF_PHYS%OUT%ALB, ZALBDIR, ZFLU_EMIS, YDVARS%GEOMETRY%GELAM%T0, YDVARS%GEOMETRY%GEMU%T0,                 &
        & ZRDG_MU0, ZRDG_MU0LU, YLMF_PHYS_BASE_STATE%YGSP_RR%T, ZDECRD, ZCLCT_RAD, YDMF_PHYS%OPT%GDEOSI,                      &
        & YDMF_PHYS%OPT%GUEOSI, YDMF_PHYS%OPT%GMU0, YDMF_PHYS%OPT%GMU0_MIN, YDMF_PHYS%OPT%GMU0_MAX, YDMF_PHYS%OPT%GDEOTI,     &
        & YDMF_PHYS%OPT%GDEOTI2, YDMF_PHYS%OPT%GUEOTI, YDMF_PHYS%OPT%GUEOTI2, YDMF_PHYS%OPT%GEOLT, YDMF_PHYS%OPT%GEOXT,       &
        & YDMF_PHYS%OPT%GRPROX, YDMF_PHYS%OPT%GMIXP, YDMF_PHYS%OPT%GFLUXC, YDMF_PHYS%OPT%GRSURF, YDMF_PHYS_SURF%GSD_VD%PSUND, &
        & YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%FRSOC, YDMF_PHYS%OUT%FRTHC, ZFRSODS,                          &
        & YDMF_PHYS%OUT%FRSOPS, YDMF_PHYS%OUT%FRSOLU, YDMF_PHYS%OUT%FRTHDS, ZAER)
    ENDSELECT

    ! sum downward diffuse and direct solar radiation at surface
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%FRSODS(JLON)=ZFRSODS(JLON)+YDMF_PHYS%OUT%FRSOPS(JLON)
    ENDDO

!      7.5 Rayonnement Morcrette
!          Morcrette's radiation

  ELSEIF ( LRAYFM ) THEN

    LLCALLRAD=(MOD(YDCPG_DIM%KSTEP,NRADFR) == 0 )
!    IF (NCALLRAD==1) ! <== not yet
    IF (NCALLRAD==2) LLCALLRAD=(LLCALLRAD.AND.(YDCPG_DIM%KSTEP<=NSTOP-1))
!    IF (NCALLRAD==3) ! <== not yet

    ! ---- Intermittent call to radiation scheme
    IF (LLCALLRAD) THEN
      CALL RECMWF(YDGEOMETRY%YRDIMV, YDMODEL, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,      &
      & YDCPG_DIM%KSW, ZALBD, ZALBP, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF,      &
      & YDCPG_MISC%NEB, ZQO3, ZAER, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZFLU_EMIS, ZRDG_MU0M,                      &
      & ZQV, ZFLU_QSAT, YDCPG_MISC%QICE, YDCPG_MISC%QLI, ZQS, ZQR, YDMF_PHYS_SURF%GSD_VF%PLSM, YLMF_PHYS_BASE_STATE%T, &
      & YLMF_PHYS_BASE_STATE%YGSP_RR%T, ZGP2DSPP, ZP1EZDIAG, YDMF_PHYS%RAD%EMTD, YDMF_PHYS%RAD%EMTU,                   &
      & YDMF_PHYS%RAD%TRSW, YDMF_PHYS%OUT%FRTHC, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%FRSOC, YDMF_PHYS%OUT%FRSO,          &
      & ZSFSWDIR, ZSFSWDIF, ZFSDNN, ZFSDNV, ZCTRSO, ZCEMTR, ZTRSOD, ZTRSODIR, ZTRSODIF, ZPIZA_DST,                     &
      & ZCGA_DST, ZTAUREL_DST, ZAERINDS, YDVARS%GEOMETRY%GELAM%T0, YDVARS%GEOMETRY%GEMU%T0, YDCPG_GPAR%SWDIR, YDCPG_GPAR%SWDIF,            &
      & ZRDG_MU0LU, YDMF_PHYS%OUT%ALB, YDMF_PHYS%RAD%RMOON)
    ELSE
      IF (LMSE) THEN
      DO JSG=1,NSW
        ZTRSODIR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSG)=YDCPG_GPAR%SWDIR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSG)
        ZTRSODIF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSG)=YDCPG_GPAR%SWDIF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSG)
      ENDDO
      ENDIF
    ENDIF

    IF (LRAYLU)  YDMF_PHYS%OUT%FRSOLU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS%RAD%RMOON(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)

    ! ---- Flux update and radiative heating rates
    CALL RADHEAT ( YDERAD, YDERDI, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, &
    & YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, ZFLU_EMIS, YDMF_PHYS%RAD%EMTD,              &
    & ZRDG_MU0, ZQV, ZTENT, YDMF_PHYS%RAD%TRSW, ZTRSOD, YLMF_PHYS_BASE_STATE%YGSP_RR%T, TSPHY,            &
    & ZTRSODIR, ZTRSODIF, ZALBD, ZALBP, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%FRSODS,     &
    & YDMF_PHYS%OUT%FRTHDS, ZCEMTR, ZCTRSO, YDMF_PHYS%OUT%FRSOC, YDMF_PHYS%OUT%FRTHC, ZSUDU, ZSDUR,       &
    & ZDSRP, ZSFSWDIR, ZSFSWDIF, YDMF_PHYS%OUT%FRSOPS, ZFRSODS, YDMF_PHYS%OUT%FRSOPT )


    ! ---- Take into account day duration depending on altitude.
    IF(LDAYD) CALL ACDAYD(YDRIP, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, &
              & YDCPG_DIM%KTDIA, NTSSG, YDVARS%GEOMETRY%GEMU%T0, ZRDG_MU0, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, &
              & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS%OUT%FRSO)

    ! ---- Correct solar absorption as a function of pmu0.
    IF(GRSO < 1._JPRB) CALL ACRSO(YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, &
                       & YDCPG_DIM%KTDIA, NTSSG, YDVARS%GEOMETRY%GEMU%T0, ZRDG_MU0, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, &
                       & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS%OUT%FRSO)

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

  ENDIF

  IF (LRAY.OR.LRAYFM ) THEN

    ! Direct normal irradiance with securities
    DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%FRSDNI(JLON)=YDMF_PHYS%OUT%FRSOPS(JLON)
      IF (ZRDG_MU0(JLON) > 3.0E-02_JPRB) THEN
        YDMF_PHYS%OUT%FRSDNI(JLON)=YDMF_PHYS%OUT%FRSOPS(JLON)/ZRDG_MU0(JLON)
      ENDIF
      YDMF_PHYS%OUT%FRSDNI(JLON)=MAX(0.0_JPRB,YDMF_PHYS%OUT%FRSDNI(JLON))
    ENDDO
  ENDIF

  IF (LRAY.OR.LRAYFM) THEN

    ! global normal irradiance
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%FRSGNI(JLON)=YDMF_PHYS%OUT%FRSDNI(JLON)+0.5_JPRB*(                         &
       & (1.0_JPRB+ZRDG_MU0(JLON))*(YDMF_PHYS%OUT%FRSODS(JLON)-YDMF_PHYS%OUT%FRSOPS(JLON)       )+ &
       & (1.0_JPRB-ZRDG_MU0(JLON))*(YDMF_PHYS%OUT%FRSODS(JLON)-YDMF_PHYS%OUT%FRSO  (JLON,YDCPG_DIM%KFLEVG,1)))
    ENDDO

    ! mean radiant temperature
    IF (LXMRT) THEN
      CALL MEAN_RAD_TEMP(YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, ZRDG_MU0, YDMF_PHYS%OUT%FRSO(:, YDCPG_DIM%KFLEVG, 1), &
      & YDMF_PHYS%OUT%FRSODS, YDMF_PHYS%OUT%FRSOPS, YDMF_PHYS%OUT%FRTH(:, YDCPG_DIM%KFLEVG, 1), YDMF_PHYS%OUT%FRTHDS,            &
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
    CALL ACVEG ( YDPHY, YDPHY1, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YDMF_PHYS%OUT%FRSO, &
    & ZQV, ZFLU_QSAT, YLMF_PHYS_BASE_STATE%T, YDMF_PHYS_SURF%GSD_VV%PD2, YDMF_PHYS_SURF%GSD_VF%PLSM,                    &
    & YDMF_PHYS_SURF%GSD_VV%PIVEG, YDMF_PHYS_SURF%GSD_VV%PLAI, ZFLU_NEIJ, ZFLU_VEG, YDMF_PHYS_SURF%GSD_VV%PRSMIN,       &
    & ZCHROV, ZGWDCS, ZWFC, YLMF_PHYS_BASE_STATE%YGSP_RR%FC, ZWWILT, YLMF_PHYS_BASE_STATE%YGSP_SB%Q,                    &
    & ZFLU_QSATS, ZHQ, ZHTR, ZHU, YDMF_PHYS_SURF%GSD_VV%PHV, ZWLMX)
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

      CALL ACDIFUS ( YDMCC, YGFL, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,                   &
      & NTDIFU, YDCPG_DIM%KFLEVG, YSP_SBD%NLEVS, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,    &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS%OUT%FRSO,                  &
      & ZKTROV, ZKUROV, ZKNROV, ZQV, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YLMF_PHYS_BASE_STATE%T,                      &
      & YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, ZXTROV, ZXUROV, ZXPTKEROV, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, ZQL, ZQI,                               &
      & ZNEBDIFF, YLMF_PHYS_BASE_STATE%TKE, ZLMU, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZFLU_CD,                         &
      & ZFLU_CDN, ZCDROV, ZCHROV, ZCEROV, ZDSA_CPS, YDMF_PHYS%OUT%CT, ZDQSTS, ZFLU_EMIS, YDMF_PHYS_SURF%GSD_VH%PSPSH,      &
      & ZHQ, ZHTR, ZHU, YDMF_PHYS_SURF%GSD_VV%PHV, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PIVEG,                &
      & ZFLU_NEIJ, YDCPG_MISC%QS, ZFLU_QSATS, YLMF_PHYS_BASE_STATE%YGSP_SB%T, YLMF_PHYS_BASE_STATE%YGSP_RR%T,              &
      & ZFLU_VEG, ZXDROV, ZXHROV, YLMF_PHYS_BASE_STATE%YGSP_RR%W, YLMF_PHYS_BASE_STATE%YGSP_RR%IC, YDMF_PHYS%OUT%DIFTQ,    &
      & YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV, YDMF_PHYS%OUT%DIFTQL,          &
      & YDMF_PHYS%OUT%DIFTQN, ZTENDPTKE, YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FCLL, YDMF_PHYS%OUT%FCLN,                      &
      & YDMF_PHYS%OUT%FCS, ZFLU_FEVI, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FTR,       &
      & ZDSA_LHS, YDMF_PHYS%OUT%FRTH, YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H)

    ELSE

      IF ( LPTKE ) THEN
        IF ( LMSE.AND.LCALLSFX ) THEN       
          IF (YDCPG_DIM%KSTEP == 0) THEN
            ZFLU_CD(:)=ZFLU_CDN(:)   ! very first approximation
          ELSE
            ZFLU_CD(:)=MAX(YDCPG_GPAR%CD(:),ZEPS0)
          ENDIF
        ENDIF
        CALL ACPTKE(YGFL, YDLDDH, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,                                &
        & YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG, YDCPG_DIM%KSTEP, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,    &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, ZKNROV, YLMF_PHYS_BASE_STATE%T,              &
        & YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%TKE,                         &
        & ZLMU, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZF_EPS, ZFUN_TTE, ZMRIPP, ZMRIFPP, ZRHS,                                        &
        & ZMN2_ES, ZMN2_EQ, ZRRCOR, YDVARS%FQTUR%T0, YDVARS%FSTUR%T0, YDVARS%SHTUR%T0, ZFLU_CD, YDMF_PHYS%OUT%GZ0,                      &
        & ZKUROV, ZKTROV, ZXPTKEROV, ZLMLTILD, YDVARS%MXL%T0, ZLML, ZTH_FUN, ZWW_FUN, ZTENDPTKE, YDVARS%TTE%T0,                         &
        & YDCPG_MISC%FTCNS)
      ENDIF

      CALL ACDIFV1 ( YGFL, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA,                   &
      & YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF,                                 &
      & ZCP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZKTROV, ZKQROV, ZKUROV, ZQV, ZQL, ZQI, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, &
      & YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, ZSVM, ZXTROV, ZXUROV,                                &
      & ZXDROV, ZXHROV, ZEDMFS, ZEDMFQ, ZEDMFU, ZEDMFV, ZMF_UP, ZXURO, ZXQRO, ZXTRO, ZCFAQ, ZCFAS,                                   &
      & ZCFATH, ZCFAU, ZCFASV, ZCFBQ, ZCFBS, ZCFBTH, ZCFBU, ZCFBV, ZCFBSV, ZDSE, ZQT)   
        

      IF ( LMSE.AND.LCALLSFX ) THEN

        IF (LRAYFM) THEN
          ZCARDI=RCARDI
        ELSEIF (LRAY) THEN
          ZCARDI=QCO2
          DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            ZSFSWDIF(JLON,1)=ZFRSODS(JLON)
            ZSFSWDIR(JLON,1)=YDMF_PHYS%OUT%FRSOPS(JLON)
          ENDDO
        ENDIF

        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZRHODREFM(JLON)=YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON,YDCPG_DIM%KFLEVG)/(YLMF_PHYS_BASE_STATE%T(JLON,YDCPG_DIM%KFLEVG)*YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R(JLON,YDCPG_DIM%KFLEVG))
          ZDEPTH_HEIGHT(JLON,:)=(YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF(JLON,:)-YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,YDCPG_DIM%KFLEVG))/RG
          ZZS(JLON)=YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,YDCPG_DIM%KFLEVG)/RG
        ENDDO

        IRR=2

        CALL ARO_GROUND_PARAM( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP, YDCPG_DIM%KFDIA-YDCPG_DIM%KIDIA+1, YDCPG_DIM%KIDIA,                                                     &
        & YDCPG_DIM%KFDIA, YDCPG_DIM%KSTEP, IRR, NSW, NGFL_EXT, NDGUNG, NDGUXG, NDLUNG, NDLUXG,                                                                          &
        & LSURFEX_KFROM, LMPA, CCOUPLING, LLXFUMSE, NINDAT, ZRHGMT, ZSTATI, RSOVR, RCODEC, RSIDEC, YDVARS%GEOMETRY%RINDX%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                      &
        & YDVARS%GEOMETRY%RINDY%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YLMF_PHYS_BASE_STATE%U(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),                   &
        & YLMF_PHYS_BASE_STATE%V(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),                                                                    &
        & YLMF_PHYS_BASE_STATE%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),                                                                    &
        & YLMF_PHYS_BASE_STATE%Q(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),                                                                    &
        & ZSVM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG, 1:NGFL_EXT),                                                                          &
        & ZCARDI, ZRHODREFM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG), &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),                                                      &
        & ZDTMSE, ZDEPTH_HEIGHT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG), ZZS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                                &
        & XZSEPS, ZRDG_MU0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZRDG_MU0N(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                                                 &
        & YDVARS%GEOMETRY%GELAM%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDVARS%GEOMETRY%GEMU%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                                               &
        & XSW_BANDS, ZSRAIN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZSSNOW(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                                                   &
        & ZSGROUPEL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%FRTHDS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                                             &
        & ZSFSWDIF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:NSW), ZSFSWDIR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:NSW),                                                            &
        & ZCFAQ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG), ZCFATH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),         &
        & ZCFAU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG), ZCFBQ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),          &
        & ZCFBTH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG), ZCFBU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG),         &
        & ZCFBV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG:YDCPG_DIM%KFLEVG), YDMF_PHYS%OUT%FCS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1),                              &
        & ZFEV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZSFSV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:NGFL_EXT),                                                                     &
        & ZSFCO2(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZFMDU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZFMDV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                       &
        & ZALBP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:NSW), ZALBD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:NSW),                                                                  &
        & ZFLU_EMIS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZTSN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%FRTH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, 1)    &
        &                         )               !orographic shadowing

! Opposite water vapor flux
        ZFEVS(:)=-ZFEV(:)

        CALL ARO_GROUND_DIAG( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP, YDCPG_DIM%KFDIA-YDCPG_DIM%KIDIA+1, YDCPG_DIM%KIDIA,                                           &
        & YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, -1, NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM, ZZS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                         &
        & ZFEVS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YLMF_PHYS_BASE_STATE%U(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG),                                &
        & YLMF_PHYS_BASE_STATE%V(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG), ZDEPTH_HEIGHT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG),    &
        & YDMF_PHYS%OUT%FRTH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, 1), YDMF_PHYS%OUT%FRSO(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, 1), &
        & YDVARS%GEOMETRY%RINDX%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDVARS%GEOMETRY%RINDY%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                                   &
        & YDCPG_MISC%QS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%GZ0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                                 &
        & YDMF_PHYS%OUT%GZ0H(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%TCLS (YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                          &
        & YDMF_PHYS%OUT%QCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%RHCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                          &
        & YDMF_PHYS%OUT%UCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%VCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                           &
        & YDMF_PHYS%OUT%NUCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%NVCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                         &
        & YDMF_PHYS%OUT%FCLL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1), YDMF_PHYS%OUT%FCLN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1),                                     &
        & YDMF_PHYS%OUT%FEVL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1), YDMF_PHYS%OUT%FEVN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1),                                     &
        & ZSSO_STDEV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZTWSNOW(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZBUDTH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                     &
        & ZBUDSO(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZFCLL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZTOWNS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                           &
        & ZCD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)                        )

        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          YDCPG_GPAR%GZ0(JLON)   = YDMF_PHYS%OUT%GZ0 (JLON)
          YDCPG_GPAR%GZ0H(JLON)  = YDMF_PHYS%OUT%GZ0H(JLON)
          YDCPG_GPAR%VEMIS(JLON) = ZFLU_EMIS(JLON)
          YDCPG_GPAR%VQS(JLON)   = YDCPG_MISC%QS  (JLON)
          YDCPG_GPAR%CD(JLON)    = ZCD  (JLON)
        ENDDO

        DO JSG=1,NSW
          YDCPG_GPAR%ALBSCA(:,JSG) = ZALBD(:,JSG)
          YDCPG_GPAR%ALBDIR(:,JSG) = ZALBP(:,JSG)
        ENDDO

        DO JSG  = 1, NTSSG+1
          DO JLEV = 0, YDCPG_DIM%KFLEVG
            DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
              YDMF_PHYS%OUT%FRTH(JLON,JLEV,JSG)=YDMF_PHYS%OUT%FRTH(JLON,JLEV,JSG)+ZBUDTH(JLON)
            ENDDO
          ENDDO
        ENDDO

! calculation of variables for the old "ISBA" atmosphere scheme
        IF (.NOT. LELAM) THEN
          CALL ARO_GROUND_DIAG_2ISBA( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP, YDCPG_DIM%KFDIA-YDCPG_DIM%KIDIA+1,                                         &
          & YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM, YDVARS%GEOMETRY%RINDX%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),      &
          & YDVARS%GEOMETRY%RINDY%T0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PARG,                               &
          & YDMF_PHYS_SURF%GSD_VV%PSAB, YDMF_PHYS_SURF%GSD_VV%PD2, ZTSN, ZTWSNOW, YDMF_PHYS_SURF%GSP_SB%PT_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1), &
          & YDMF_PHYS_SURF%GSP_RR%PW_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS_SURF%GSP_SB%PQ_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1),         &
          & YDMF_PHYS_SURF%GSP_RR%PIC_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS_SURF%GSP_SB%PTL_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1),       &
          & YDMF_PHYS_SURF%GSP_RR%PFC_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS_SURF%GSP_SG%PA_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1),        &
          & YDMF_PHYS_SURF%GSP_SG%PR_T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1), ZHV2 )
          DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            YDMF_PHYS_SURF%GSD_VV%PHV(JLON) = ZHV2(JLON)
          ENDDO
        ENDIF

        IF (LDIFCONS.AND..NOT.LNODIFQC) THEN
          CALL ACAA1 ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA,    &
          & YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, ZCOEFN, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, &
          & ZQL, ZQI, YLMF_PHYS_BASE_STATE%T, ZALPHA1, ZCOEFA, ZLVT, ZQICE)
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
          CALL ARP_GROUND_PARAM ( YDMCC, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,        &
          & YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG, YSP_SBD%NLEVS, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, ZCP,                  &
          & YDMF_PHYS%OUT%FRSO, ZQV, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YLMF_PHYS_BASE_STATE%T,                  &
          & YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, ZXURO, ZXQRO, ZXTRO, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, &
          & ZQL, ZQI, ZNEBDIFF, ZFLU_CD, ZFLU_CDN, ZCDROV, ZCHROV, ZCEROV, ZDSA_CPS, YDMF_PHYS%OUT%CT,                 &
          & ZDQSTS, ZFLU_EMIS, YDMF_PHYS_SURF%GSD_VH%PSPSH, ZHQ, ZHTR, ZHU, YDMF_PHYS_SURF%GSD_VV%PHV,                 &
          & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PIVEG, ZFLU_NEIJ, YDCPG_MISC%QS,                         &
          & ZFLU_QSATS, YLMF_PHYS_BASE_STATE%YGSP_SB%T, YLMF_PHYS_BASE_STATE%YGSP_RR%T, ZFLU_VEG,                      &
          & ZXDROV, ZXHROV, YLMF_PHYS_BASE_STATE%YGSP_RR%W, YLMF_PHYS_BASE_STATE%YGSP_RR%IC, ZDSE,                     &
          & ZCFAS, ZCFAU, ZCFBS, ZCFBU, ZCFBV, ZCFBQ, ZCOEFA, ZALPHA1, ZLVT, ZQICE, ZDIFWQ, ZDIFWS,                    &
          & ZFMDU, ZFMDV, ZSC_FEVI, ZSC_FEVN, ZSC_FCLL, ZSC_FCLN, YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FCLL,             &
          & YDMF_PHYS%OUT%FCLN, YDMF_PHYS%OUT%FCS, ZFLU_FEVI, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN,                  &
          & YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FTR, ZDSA_LHS, YDMF_PHYS_SURF%GSD_VH%PQSH, YDMF_PHYS%OUT%FRTH,           &
          & YDMF_PHYS_SURF%GSD_VF%PZ0F, YDMF_PHYS_SURF%GSD_VV%PZ0H)

          DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            YDMF_PHYS%OUT%DIFTQ(JLON,YDCPG_DIM%KFLEVG)=ZDIFWQ(JLON)
            YDMF_PHYS%OUT%DIFTS(JLON,YDCPG_DIM%KFLEVG)=ZDIFWS(JLON)
          ENDDO

          IF (LEDMFI) THEN
            ZCFBQ(:,:)  = ZCFBQ_G(:,:)
            ZCFBS(:,:)  = ZCFBS_G(:,:) 
            ZCFBU(:,YDCPG_DIM%KFLEVG)  = ZCFBU_G(:,YDCPG_DIM%KFLEVG)
            ZCFBV(:,YDCPG_DIM%KFLEVG)  = ZCFBV_G(:,YDCPG_DIM%KFLEVG) 
          ENDIF   

        ENDIF   ! <== LSFORCS

      ELSEIF (.NOT. LCALLSFX) THEN
        YDMF_PHYS%OUT%FCS(:,:) = 0.0_JPRB
        ZFEV(:)   = 0.0_JPRB

      ENDIF  !LMSE.AND.LCALLSFX

      CALL ACDIFV2 ( YGFL, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA, &
      & YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, ZCFAQ, ZCFAS, ZCFAU, ZCFASV, ZCFBQ,                &
      & ZCFBS, ZCFBU, ZCFBV, ZCFBSV, ZKTROV, ZKQROV, ZKQLROV, ZKUROV, ZDSE, ZQT, YLMF_PHYS_BASE_STATE%U,           &
      & YLMF_PHYS_BASE_STATE%V, ZPOID, YLMF_PHYS_BASE_STATE%T, ZQL, ZQI, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,  &
      & ZCOEFA, ZALPHA1, ZLVT, ZQICE, ZSFSV, YDMF_PHYS%OUT%FCS, ZFEV, ZFMDU, ZFMDV, ZTSN, ZXHROV,                  &
      & ZDIFEXT, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV,               &
      & YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%STRCU, &
      & YDMF_PHYS%OUT%STRCV, YDVARS%SHTUR%T0          )

      IF (LEDKF) THEN
        IF (LSFORCS) THEN
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
           YDMF_PHYS_SURF%GSD_SFL%PGROUP(JLON,1) = YDMF_PHYS%OUT%FCS(JLON,1)
           YDMF_PHYS_SURF%GSD_SFL%PGROUP(JLON,2) = ZFEV(JLON)
        ENDDO
        ELSE
          DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA             
            YDMF_PHYS_SURF%GSD_SFL%PGROUP(JLON,1) = YDMF_PHYS%OUT%DIFTS(JLON,YDCPG_DIM%KFLEVG)
            YDMF_PHYS_SURF%GSD_SFL%PGROUP(JLON,2) = YDMF_PHYS%OUT%DIFTQ(JLON,YDCPG_DIM%KFLEVG)
          ENDDO
        ENDIF   
      ENDIF

      IF ( LMSE ) THEN
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          YDCPG_GPAR%VTS(JLON) = ZTSN (JLON)
          YDMF_PHYS_SURF%GSP_SG%PF_T1(JLON,1) = ZTWSNOW (JLON)
        ENDDO
      ENDIF

      ! First compute horizontal exchange coefficients for momentum:
      !  (there's mo TOMs contribution, thus has to be done at latest here)
      IF (L3DTURB) THEN
        CALL ACTKECOEFKH(YDRIP, YDMODEL%YRML_PHY_MF, YDGEOMETRY%YREGEO, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA,                            &
        & YDCPG_DIM%KLON, YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%TKE, ZTENDPTKE, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, &
        & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V,                &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YLMF_PHYS_BASE_STATE%DIV, YLMF_PHYS_BASE_STATE%VOR,                                &
        & YDVARS%U%DL, YDVARS%V%DL, YLMF_PHYS_BASE_STATE%YCPG_PHY%W, YLMF_PHYS_BASE_STATE%YCPG_PHY%WL,                               &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%WM, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,               &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZLML, ZFHORM, ZFHORH, ZKUROV, ZRTI, ZFLU_CD,                                       &
        & ZCDROV, LPTKE, ZKUR_KUROV_H, ZKUR_KTROV_H, ZRHS)
      ENDIF

      IF (LCOEFKTKE) THEN
        IF (NDIFFNEB == 1) THEN
          ZCOEFA(:,:) = ZBNEBQ(:,:)
        ELSEIF (NDIFFNEB == 4) THEN
          ZCOEFA(:,:) = ZBNEBCVPP(:,:)
        ENDIF
        CALL ACDIFV3 ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA,           &
        & YDCPG_DIM%KFLEVG, YDCPG_DIM%KSTEP, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI,    &
        & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, ZCP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZKTROV,                       &
        & ZXTROV, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,                          &
        & ZCHROV, ZXHROV, ZCOEFA, ZQV, ZQL, ZQI, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YLMF_PHYS_BASE_STATE%T, &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZLML,                          &
        & ZTH_FUN, ZWW_FUN, ZF_EPS, YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%DIFTQ, YLMF_PHYS_BASE_STATE%TKE,                  &
        & ZTENDPTKE, ZMN2_ES, ZMN2_EQ, ZMN2_DS, ZMN2_DQ, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTQN,                     &
        & ZDIFWQ, ZDIFWS, YDMF_PHYS%OUT%FCLL, YDMF_PHYS%OUT%FCLN, YDMF_PHYS%OUT%FCS, ZFLU_FEVI, YDMF_PHYS%OUT%FEVL,      &
        & YDMF_PHYS%OUT%FEVN, YDMF_PHYS%OUT%GZ0, ZRTI, ZSC_FEVI, ZSC_FEVN, ZSC_FCLL, ZSC_FCLN,                           &
        & ZTSTAR, ZTSTAR2, ZTSTARQ, ZTSTAR2Q)

        !store fluxes and shear term
        !GFL fields are on full levels, fluxes on half levels
        IF (YFQTUR%LGP.AND.YFSTUR%LGP) THEN
          DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
            DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
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
        CALL ACTKECOEFKH(YDRIP, YDMODEL%YRML_PHY_MF, YDGEOMETRY%YREGEO, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA,                            &
        & YDCPG_DIM%KLON, YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%TKE, ZTENDPTKE, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, &
        & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V,                &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YLMF_PHYS_BASE_STATE%DIV, YLMF_PHYS_BASE_STATE%VOR,                                &
        & YDVARS%U%DL, YDVARS%V%DL, YLMF_PHYS_BASE_STATE%YCPG_PHY%W, YLMF_PHYS_BASE_STATE%YCPG_PHY%WL,                               &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%WM, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,               &
        & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZLML, ZFHORM, ZFHORH, ZKUROV, ZRTI, ZFLU_CD,                                       &
        & ZCDROV, LPTKE, ZKUR_KUROV_H, ZKUR_KTROV_H, ZRHS)
      ENDIF

    ENDIF

!-----------------------------------------------------------------------------
!   THE DEEP CONVECTION WILL SEE THE SHALLOW PART FROM KFB AS IT IS WITH LOUIS
!   SCHEME AND THE MODIFIED RI
!----------------------------------------------------------------------------
    IF (LCVPPKF.OR.(LEDKF .AND. .NOT. LEDMFI)) THEN
      DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
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
    DO JLEV = YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
!DEC$ IVDEP
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA

        ZT(JLON,JLEV)=YLMF_PHYS_BASE_STATE%T(JLON,JLEV)-ZIPOI(JLON,JLEV)/YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP(JLON,JLEV)&
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
        IF(LRKCDEV.AND.YDCPG_DIM%KSTEP>0) THEN
          ZDTRAD(JLON,JLEV)=ZIPOI(JLON,JLEV)/YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP(JLON,JLEV)* (&
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
    DO JLEV = YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA

        ZT (JLON,JLEV)= YLMF_PHYS_BASE_STATE%T(JLON,JLEV)
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
      ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=MIN(XMAXLM,MAX(XMINLM,YDMF_PHYS_SURF%GSD_VH%PPBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)))
      CALL ACPBLHTM ( YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTDIFU, YDCPG_DIM%KFLEVG,                  &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,     &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS%OUT%STRTU,                   &
      & YDMF_PHYS%OUT%STRTV, ZPCLS, YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%QCLS, YDMF_PHYS%OUT%FCS, YDMF_PHYS%OUT%FCLL,         &
      & YDMF_PHYS%OUT%FCLN, ZDSA_LHS, YDMF_PHYS%OUT%GZ0, IMOC_CLPH, ZBLH, ZUGST, ZVGST)
      YDMF_PHYS_SURF%GSD_VH%PPBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      YDMF_PHYS%OUT%CLPH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      IF (.NOT.LRAFTUR) THEN
        YDMF_PHYS%OUT%UGST(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=ZUGST(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
        YDMF_PHYS%OUT%VGST(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=ZVGST(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
      ENDIF
    ENDIF

    IF(LNEBCO.AND.LNEBR)THEN
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZPREN(JLON)=YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD(JLON,YDCPG_DIM%KFLEVG)
      ENDDO
      CALL ACPBLH ( YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTDIFU, YDCPG_DIM%KFLEVG,                     &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP,        &
      & ZQV, YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, &
      & YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV, ZPREN, YLMF_PHYS_BASE_STATE%YGSP_RR%T, YDCPG_MISC%QS,                     &
      & YDMF_PHYS%OUT%FCS, YDMF_PHYS%OUT%FCLL, YDMF_PHYS%OUT%FCLN, ZDSA_LHS, ZRTI, YDMF_PHYS_SURF%GSD_VH%PPBLH              &
      &                                                                                  )
    ENDIF
    ! ------------------------------------------------------------------
    ! UPDATE PASSIFS SCALAIRS DUE TO THE TURBULENT DIFFUSION PROCESSES
    ! ------------------------------------------------------------------
    IF (LMDUST.AND.(NGFL_EXT/=0)) THEN
      DO JGFL=1,NGFL_EXT
        DO JLON=1,YDCPG_DIM%KLON
          DO JLEV=1,YDCPG_DIM%KFLEVG
            ZSVM(JLON,JLEV,JGFL)=ZSVM(JLON,JLEV,JGFL)-ZIPOI(JLON,JLEV)&
          & *(ZDIFEXT(JLON,JLEV,JGFL)-ZDIFEXT(JLON,JLEV-1,JGFL))
            ZSVM(JLON,JLEV,JGFL)=MAX(ZSVM(JLON,JLEV,JGFL),0.0_JPRB)
          ENDDO
        ENDDO
      ENDDO
    ENDIF  ! LMDUST
  ENDIF ! LVDIF

!     DIAGNOSTIC SUPPLEMENTAIRE FLUX DE RAYONNEMENT (PFRTHDS ET PFRSOPT)
!     ADDITIONAL DIAGNOSTICS OF RADIATIVE FLUXES (PFRTHDS AND PFRSOPT)

!DEC$ IVDEP
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    IF(.NOT.LMSE)YDMF_PHYS%OUT%FRTHDS(JLON)=YDMF_PHYS%OUT%FRTH(JLON,YDCPG_DIM%KFLEVG,1)/ZFLU_EMIS(JLON)+RSIGMA*YLMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)**4
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
    & ZFLU_EMIS, YDMF_PHYS_SURF%GSD_VF%PLSM, ZFLU_NEIJ, ZQV, YDCPG_MISC%QS, YLMF_PHYS_BASE_STATE%YGSP_RR%T,          &
    & ZCHROV, ZDQSTS, YDMF_PHYS%OUT%DRNSHF)
  ENDIF

!*
!     ------------------------------------------------------------------
!     9.- TRAINEE DES ONDES DE GRAVITE INDUITES PAR LE RELIEF
!     ------------------------------------------------------------------
  IF ( LGWD ) THEN
    CALL ACDRAG ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTDRAG, YDCPG_DIM%KFLEVG,         &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
    & ZNBVNO, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V,                     &
    & YDVARS%GEOMETRY%RCORI%T0, YDMF_PHYS_SURF%GSD_VF%PGETRL, ZGWDCS, YDMF_PHYS_SURF%GSD_VF%PVRLAN, YDMF_PHYS_SURF%GSD_VF%PVRLDI,    &
    & YDMF_PHYS%OUT%STRDU, YDMF_PHYS%OUT%STRDV, ZTRAJGWD)
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

     CALL WRPHTRAJTM_NL(YDGEOMETRY, YDSIMPHL, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, PTRAJ_PHYS,                          &
     & ZKTROV_SAVE, ZKUROV_SAVE, ZCDROV_SAVE, ZCHROV_SAVE, ZTRAJGWD, YLMF_PHYS_BASE_STATE%L, YLMF_PHYS_BASE_STATE%I, &
     & YLMF_PHYS_BASE_STATE%R, YLMF_PHYS_BASE_STATE%S, ZQLIS, ZNEBS)
  ENDIF

!     ------------------------------------------------------------------
!     10.- PRECIPITATIONS STRATIFORMES.
!     ---------------------------------

  IF ( LSTRA.AND.(.NOT.LSTRAPRO) ) THEN
    CALL ACPLUIE ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTPLUI, YDCPG_DIM%KFLEVG,       &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZQV, ZMSC_QW, YLMF_PHYS_BASE_STATE%T, &
    & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPEVPSL,        &
    & YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN)
  ENDIF

  IF ( LSTRAS ) THEN
    CALL ACPLUIS ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTPLUI, YDCPG_DIM%KFLEVG, &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZNEBS, ZQV,                    &
    & ZQLIS, ZMSC_QW, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%T, YDMF_PHYS%OUT%FCSQL,             &
    & YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN,  &
    & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN       )
  ENDIF

  IF (L3MT.OR.LSTRAPRO) THEN

    IF (LCOEFKTKE.OR.LLREDPR) THEN
      CALL ACTQSATS ( YDPHY, YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTQSAT, YDCPG_DIM%KFLEVG, &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, ZT,                    &
      & ZGEOSLC, ZLH, ZLSCPE, ZQSAT, ZQW, ZRH, ZTW)
      CALL ACNEBCOND ( YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,                            &
      & NTPLUI, YDCPG_DIM%KFLEVG, LLREDPR, ZHUC, ZVETAF, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,       &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDCPG_MISC%RH, ZBLH, ZQV, ZQI, ZQL, ZQW, ZT, ZNEBCH,                            &
      & YDVARS%GEOMETRY%GM%T0, YLMF_PHYS_BASE_STATE%YGSP_RR%T, ZQLIS, ZNEBS, ZRHCRI, ZRH, ZQSATS, ZICEFR1,                                &
      & ZQLIS0, ZNEBS0)
      CALL ACCDEV ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTPLUI, YDCPG_DIM%KFLEVG,            &
      & ZVETAF, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZQV, ZQI, ZQL,                            &
      & ZQS, ZQR, ZQG, ZQW, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, ZT, ZMSC_TW, ZNEBS, ZRHCRI, ZICEFR1,                           &
      & ZQSATS, ZNEBCH, ZRHDFDA, ZRKTTEND, ZRKQVTEND, ZRKQCTEND, ZQLI_CVPP, ZNEB_CVPP, YDMF_PHYS_SURF%GSD_VF%PLSM,              &
      & ZFLU_NEIJ, YDVARS%GEOMETRY%GM%T0, YLMF_PHYS_BASE_STATE%YGSP_RR%T, ZDECRD_MF, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN,          &
      & YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN,           &
      & YDMF_PHYS%OUT%FPEVPSG, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLSG,                                   &
      & ZSEDIQL, ZSEDIQI, YDMF_PHYS%OUT%DIAGH)
    ELSE
      CALL ACCDEV ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTPLUI, YDCPG_DIM%KFLEVG,            &
      & ZVETAF, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, ZQV, ZQI, ZQL,                            &
      & ZQS, ZQR, ZQG, ZMSC_QW, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%T, ZMSC_TW,                           &
      & ZNEBS, ZRHCRI, ZICEFR1, ZQSATS, ZNEBCH, ZRHDFDA, ZRKTTEND, ZRKQVTEND, ZRKQCTEND, ZQLI_CVPP,                             &
      & ZNEB_CVPP, YDMF_PHYS_SURF%GSD_VF%PLSM, ZFLU_NEIJ, YDVARS%GEOMETRY%GM%T0, YLMF_PHYS_BASE_STATE%YGSP_RR%T,                          &
      & ZDECRD_MF, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FCSQL,                       &
      & YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG, YDMF_PHYS%OUT%FPLSL,          &
      & YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLSG, ZSEDIQL, ZSEDIQI, YDMF_PHYS%OUT%DIAGH)
    ENDIF

  ENDIF

  IF (L3MT) THEN


    !  ----------------------------------------------------------------
    !   UPDATE HUMIDITY VARIABLES BY RESOLVED CONDENSATION FLUXES
    !  ----------------------------------------------------------------
!cdir unroll=8
    DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZDQL=ZIPOI(JLON,JLEV)* (&
          & YDMF_PHYS%OUT%FCSQL(JLON,JLEV)-YDMF_PHYS%OUT%FCSQL(JLON,JLEV-1))
        ZQL(JLON,JLEV)=ZQL(JLON,JLEV)+ZDQL
        ZDQI=ZIPOI(JLON,JLEV)* (&
          & YDMF_PHYS%OUT%FCSQN(JLON,JLEV)-YDMF_PHYS%OUT%FCSQN(JLON,JLEV-1))
        ZQI(JLON,JLEV)=ZQI(JLON,JLEV)+ZDQI
        ZQV(JLON,JLEV)=ZQV(JLON,JLEV)-ZDQI-ZDQL

        ZT(JLON,JLEV)=ZT(JLON,JLEV)+(ZLHV(JLON,JLEV)*ZDQL&
                   & +ZLHS(JLON,JLEV)*ZDQI)/YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP(JLON,JLEV)
      ENDDO
    ENDDO

!rb : Store values for next time-step: all cases
    IF(LRKCDEV) THEN
      DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
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
          DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            ZMOD(JLON)=1._JPRB/(1.0_JPRB+YDVARS%GEOMETRY%GM%T0(JLON)*TEQK)
         ENDDO
      ELSE
         ZMOD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=1._JPRB
      ENDIF
!cdir unroll=8
      DO JLEV = YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG
        DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
          ZRDG_CVGQ(JLON,JLEV) = ZRDG_CVGQ(JLON,JLEV)*ZMOD(JLON) &
           & -RG*YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP(JLON,JLEV)&
           & *(YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)-YDMF_PHYS%OUT%DIFTQ(JLON,JLEV-1)&
           & +ZFCQVNG(JLON,JLEV)-ZFCQVNG(JLON,JLEV-1))
        ENDDO
      ENDDO
    ENDIF
    IF (LCAPE) THEN
      IF (LCAMOD) THEN
         DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            ZTAUX(JLON)=RTCAPE*(YDVARS%GEOMETRY%GM%T0(JLON)*TEQC)
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

    CALL ACTQSATS ( YDPHY, YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTQSAT, YDCPG_DIM%KFLEVG, &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, ZT, ZGEOSLC,           &
    & ZLH, ZLSCPE, ZQSAT, ZQW, ZRH, ZTW)

    IF (.NOT.LCVCSD) THEN
      CALL ACCVUD ( YDMODEL%YRML_PHY_EC%YRECUCONVCA, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA,                           &
      & YDCPG_DIM%KLON, YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF,              &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, ZRDG_CVGQ, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZQV, ZQI, ZQL, ZQSAT, ZQW, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,                   &
      & ZT, ZTW, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YLMF_PHYS_BASE_STATE%VOR, YDVARS%DAL%T0,                             &
      & YLMF_PHYS_BASE_STATE%YGSP_RR%T, ZTAUX, YDVARS%GEOMETRY%RCORI%T0, YDVARS%GEOMETRY%GM%T0, ZDECRD_MF, YDMF_PHYS%OUT%DIFCQ,                           &
      & YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN,                    &
      & YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, ZATSLC, INLAB, YDVARS%UNEBH%T0, ZDETFI, ZTU, ZQU,                                   &
      & ZUU, ZVU, INND, YDMF_PHYS%OUT%CAPE, YDMF_PHYS%OUT%CUCONVCA, YDMF_PHYS%OUT%NLCONVCA, ZGEOSLC,                                  &
      & YDVARS%UEN%T0, YDVARS%UAL%T0, YDVARS%UOM%T0)

    ! Initialize temperature correction
      ZTCORR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)=ZTU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)-ZT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
    ! Initialize Environment horizontal velocity
      ZU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)=YLMF_PHYS_BASE_STATE%U(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
      ZV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)=YLMF_PHYS_BASE_STATE%V(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
    ELSE
! PUT IN ZTCORR THE CONDENSATE USED IN TRIGGERING
      DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
       DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZTCORR(JLON,JLEV)=MAX(0._JPRB,ZQI(JLON,JLEV)+ZQL(JLON,JLEV) -&
          &  MAX(0._JPRB,YLMF_PHYS_BASE_STATE%I(JLON,JLEV)+YLMF_PHYS_BASE_STATE%L(JLON,JLEV)) )
       ENDDO
      ENDDO
! PASS T9 CONDENSATES TO ACCSU
      CALL ACCSU ( YDMODEL%YRML_PHY_EC%YRECUCONVCA, YDMODEL%YRML_PHY_EC%YRECUMF, YDMODEL%YRML_PHY_MF,                                &
      & YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, &
      & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,                 &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,            &
      & ZRDG_CVGQ, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR,                                   &
      & ZQV, ZTCORR, YLMF_PHYS_BASE_STATE%I, YLMF_PHYS_BASE_STATE%L, ZQSAT, ZQW, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,            &
      & ZT, ZTW, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YLMF_PHYS_BASE_STATE%VOR, YDCPG_DYN0%CTY%VVEL(:, 1:),               &
      & YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YLMF_PHYS_BASE_STATE%YGSP_RR%T, ZTAUX, YDVARS%GEOMETRY%RCORI%T0,                             &
      & ZDECRD_MF, YDMF_PHYS_SURF%GSD_VK%PUDGRO, YDMF_PHYS%OUT%CUCONVCA, YDMF_PHYS%OUT%NLCONVCA, YDMF_PHYS%OUT%DIFCQ,                &
      & YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN,                   &
      & YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, ZATSLC, INLAB, YDVARS%UNEBH%T0, ZDETFI, ZTU, ZQU,                                  &
      & ZUU, ZVU, INND, YDMF_PHYS%OUT%CAPE, ZGEOSLC, YDVARS%UEN%T0, YDVARS%UAL%T0, YDVARS%UOM%T0)

    ! Initialize temperature correction
      ZTCORR(:,:)=ZTU(:,:)
    ! Initialize Environment horizontal velocity
      ZU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)=YLMF_PHYS_BASE_STATE%U(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
      ZV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)=YLMF_PHYS_BASE_STATE%V(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
    ENDIF


    CALL ACUPU(YDPHY, YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG, &
    & YDVARS%UAL%T0, YDVARS%UNEBH%T0, ZDETFI, ZPOID, ZIPOI, ZLHV, ZLHS, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP,      &
    & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%FCCQL,   &
    & YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, ZSIGP, ZSIGPC, ZNEBS, ZQI, ZQL,               &
    & ZQV, ZT, YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQLNG)

!rb store "-B" state: no protection of Ncv
    IF(LRKCDEV.AND.(.NOT.LNEBCV)) THEN
      DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
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

    DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA

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

    CALL ACTQSATS ( YDPHY, YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTQSAT, YDCPG_DIM%KFLEVG, &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, ZT, ZGEOSLC,           &
    & ZLH, ZLSCPE, ZQSAT, ZQW, ZRH, ZTW)

    DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZDQ(JLON,JLEV)=ZQW(JLON,JLEV)&
         & -(ZQV(JLON,JLEV)+ZQI(JLON,JLEV)+ZQL(JLON,JLEV))
        ZDQM(JLON,JLEV)=ZQW(JLON,JLEV)
      ENDDO
    ENDDO

    CALL APLMPHYS( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA,              &
    & YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF,                    &
    & YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, ZQI, ZQL, ZQS, ZQR, ZQG, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,          &
    & ZT, ZIPOI, ZDQ, ZDQM, ZLHS, ZLHV, ZNEBS, ZPOID, ZRCVOTT, ZTCORR, YDMF_PHYS_SURF%GSD_VF%PLSM,                      &
    & ZFLU_NEIJ, YLMF_PHYS_BASE_STATE%YGSP_RR%T, ZDECRD_MF, ZFCQL, ZFCQI, ZFALLR, ZFALLS, ZFALLG, YDMF_PHYS%OUT%FPFPSL, &
    & YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG,  &
    & YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLSG, ZSEDIQL, ZSEDIQI, ZMELNET,                         &
    & ZMELGET, YDMF_PHYS%OUT%DIAGH)

    !  ------------------------------------------------------------
    !   UPDATE AFTER MICROPHYSICS - 3MT
    !  ------------------------------------------------------------

   CALL ACUPM( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA,                 &
   & YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZLHS, ZLHV, YDCPG_DYN0%CTY%EVEL, ZIPOI,                   &
   & ZPOID, YDVARS%UAL%T0, YDVARS%UOM%T0, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG,         &
   & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN,       &
   & YDMF_PHYS%OUT%FPLSG, ZFCQL, ZFCQI, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, ZOME, ZFHP, ZQV,                     &
   & ZQL, ZQI, ZQR, ZQS, ZQG, ZT, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%FCQRNG, &
   & YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FCQGNG)

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

    CALL ACTQSATS ( YDPHY, YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTQSAT, YDCPG_DIM%KFLEVG, &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, ZQV, ZT, ZGEOSLC,           &
    & ZLH, ZLSCPE, ZQSAT, ZQW, ZRH, ZTW)

    IF (LNSDO) THEN
      CALL ACNSDO(YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA,                    &
      & YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,                        &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZQV, ZQI, ZQL, ZQR, ZQS, ZQW, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,       &
      & ZSIGP, ZT, ZTW, ZU, ZV, YDCPG_DYN0%CTY%VVEL(:, 1:), ZATSLC, ZGEOSLC, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN,       &
      & YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN, ZDIFCQD, ZDIFCQLD, ZDIFCQID, ZDIFCSD, YDMF_PHYS%OUT%FPEVPCL,           &
      & YDMF_PHYS%OUT%FPEVPCN, ZSTRCUD, ZSTRCVD, YDVARS%DAL%T0, YDVARS%DOM%T0)
    ELSE
      CALL ACMODO(YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA,                  &
      & YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,                      &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, &
      & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZQV, ZQI,                          &
      & ZQL, ZQW, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,                            &
      & ZSIGP, ZT, ZTW, ZU, ZV, YDCPG_DYN0%CTY%EVEL, ZOME, ZATSLC, ZGEOSLC, ZFHP, YDMF_PHYS%OUT%FPLSL,                     &
      & YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLSG, ZDIFCQD, ZDIFCQLD, ZDIFCQID, ZDIFCSD, YDMF_PHYS%OUT%FPEVPCL,             &
      & YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPEVPCG, ZSTRCUD, ZSTRCVD, YDVARS%DAL%T0, YDVARS%DOM%T0                       &
      &                                                                                                                                      )
    ENDIF
  !  ---------------------------------------------
  !  UPDATE VARIABLES  BY DOWNDRAUGHT CONTRIBUTION
  !  ---------------------------------------------


    DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA

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

    CALL ACUPD(YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA,                   &
    & YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP,                     &
    & ZLHS, ZLHV, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, ZIPOI, ZPOID, ZFALLR, ZFALLS, ZFALLG, YDMF_PHYS%OUT%FPEVPSL,      &
    & YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPEVPCG, &
    & ZFHP, ZDIFCSD, ZDIFCQD, ZDIFCQLD, ZDIFCQID, ZQV, ZQL, ZQI, ZQR, ZQS, ZQG, ZT, YDVARS%UEN%T0, YDMF_PHYS%OUT%FCQNG,  &
    & YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FCQGNG,      &
    & YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLSG)

!rb  store "+C": case with no protection of Ncv
    IF((LRKCDEV).AND.(.NOT.LNEBCV)) THEN
      DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          YDVARS%RKTH%T0(JLON,JLEV)=YDVARS%RKTH%T0(JLON,JLEV)+ZT(JLON,JLEV)
          YDVARS%RKTQV%T0(JLON,JLEV)=YDVARS%RKTQV%T0(JLON,JLEV)+ZQV(JLON,JLEV)
          YDVARS%RKTQC%T0(JLON,JLEV)=YDVARS%RKTQC%T0(JLON,JLEV)&
          & +ZQI(JLON,JLEV)+ZQL(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

  !  PARTITION CONVECTIVE/STRATIFORM PRECIPITATION
  ! ---------------------------------------------
    DO JLEV=YDCPG_DIM%KTDIA-1,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
         YDMF_PHYS%OUT%FPLCL(JLON,JLEV)=ZSIGPC(JLON)*YDMF_PHYS%OUT%FPLSL(JLON,JLEV)
         YDMF_PHYS%OUT%FPLSL(JLON,JLEV)=YDMF_PHYS%OUT%FPLSL(JLON,JLEV)-YDMF_PHYS%OUT%FPLCL(JLON,JLEV)
         YDMF_PHYS%OUT%FPLCN(JLON,JLEV)=ZSIGPC(JLON)*YDMF_PHYS%OUT%FPLSN(JLON,JLEV)
         YDMF_PHYS%OUT%FPLSN(JLON,JLEV)=YDMF_PHYS%OUT%FPLSN(JLON,JLEV)-YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
      ENDDO
    ENDDO
    IF (LGRAPRO) THEN
      DO JLEV=YDCPG_DIM%KTDIA-1,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
           YDMF_PHYS%OUT%FPLCG(JLON,YDCPG_DIM%KFLEVG)=0.0_JPRB !ZSIGPC(JLON)*PFPLSG(JLON,JLEV)
           YDMF_PHYS%OUT%FPLSG(JLON,JLEV)=YDMF_PHYS%OUT%FPLSG(JLON,JLEV)-YDMF_PHYS%OUT%FPLCG(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF
  ENDIF ! LCDDPRO

  ENDIF ! L3MT

! -----------------------------------------------------

  IF ( LCVRA ) THEN

    IF (LSTRAPRO) THEN
      DO JLEV = YDCPG_DIM%KTDIA-1, YDCPG_DIM%KFLEVG
        DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
          YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)+ZFCQVNG(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

    CALL ACCVIMP ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTCVIM, YDCPG_DIM%KFLEVG,      &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,     &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, &
    & ZRDG_CVGQ, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZQV,                    &
    & ZFLU_QSAT, ZMSC_QW, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,                  &
    & YLMF_PHYS_BASE_STATE%T, ZMSC_TW, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YLMF_PHYS_BASE_STATE%VOR,         &
    & YDCPG_DYN0%CTY%VVEL(:, 1:), YDVARS%GEOMETRY%GM%T0, YLMF_PHYS_BASE_STATE%YGSP_RR%T, ZTAUX, YDVARS%GEOMETRY%RCORI%T0,                    &
    & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FPLCL,           &
    & YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN,     &
    & YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, ZFPCOR, INLAB, YDMF_PHYS%OUT%CAPE, INND, YDMF_PHYS%OUT%DIFTQ,            &
    & YDMF_PHYS%OUT%DIFTS, ZGEOSLC, YDVARS%GEOMETRY%GEMU%T0)

    IF (LSTRAPRO) THEN
      DO JLEV = YDCPG_DIM%KTDIA-1, YDCPG_DIM%KFLEVG
        DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
          YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)=YDMF_PHYS%OUT%DIFTQ(JLON,JLEV)-ZFCQVNG(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF

  ELSEIF (LCVRAV3) THEN

    CALL ACCVIMP_V3 (YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTCVIM, YDCPG_DIM%KFLEVG,       &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%ALPH, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,        &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP,    &
    & ZRDG_CVGQ, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%LNPR, ZQV,                       &
    & ZMSC_QW, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YLMF_PHYS_BASE_STATE%T,        &
    & ZMSC_TW, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YLMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS_SURF%GSD_VH%PPBLH, &
    & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FPLCL,              &
    & YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN,        &
    & YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV, INLAB, INND, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS                       &
    &                                                                                                                           )

  ELSEIF (LCVTDK) THEN     ! <== IFS deep convection scheme

    DO JLEV=1,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZGEOM1(JLON,JLEV)     = YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF(JLON,JLEV)-YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,YDCPG_DIM%KFLEVG)
        ZVERVEL(JLON,JLEV)    = YDCPG_DYN0%CTY%VVEL(JLON,JLEV)*YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON,JLEV)
        ZLISUM(JLON,JLEV)     = 0.0_JPRB
        ZLCRIT_AER(JLON,JLEV) = 5.E-4_JPRB
      ENDDO
    ENDDO
    DO JLEV=0,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZGEOMH(JLON,JLEV)=YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,JLEV)-YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,YDCPG_DIM%KFLEVG)
        ZFCQLF(JLON,JLEV)=0.0_JPRB
        ZFCQLI(JLON,JLEV)=0.0_JPRB
      ENDDO
    ENDDO
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
       LLLAND(JLON)    = (YDMF_PHYS_SURF%GSD_VF%PLSM(JLON)>0.5_JPRB)
       LLSHCV(JLON)    = .FALSE.
       ZCUCONVCA(JLON) = 0.0_JPRB
       ZACPR(JLON)     = 0.0_JPRB
       ZDXTDK(JLON)    = RDELXN/YDVARS%GEOMETRY%GM%T0(JLON)
    ENDDO
    LLPTQ = .FALSE.
    CALL ACUPTQ ( YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, LLPTQ, YDMF_PHYS%OUT%FRSO,  &
    & YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS,     &
    & YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPEVPCL,  &
    & YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, &
    & YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%Q, YLMF_PHYS_BASE_STATE%YGSP_RR%T, ZTENHA,                     &
    & ZTENQVA )
    
    DO JLEV=1,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA 
         ZTENT(JLON,JLEV) = ZTENHA(JLON,JLEV)/YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP(JLON,JLEV)
         ZTENQ(JLON,JLEV) = ZTENQVA(JLON,JLEV)
         ZTENU(JLON,JLEV) = 0.0_JPRB
         ZTENV(JLON,JLEV) = 0.0_JPRB
      ENDDO
    ENDDO
    LLDSLPHY=.TRUE.
    ZVDIFTS = 0._JPRB
    ISPPN2D = 0
         
    CALL CUCALLN_MF   (YDMODEL%YRML_PHY_RAD%YRERAD, YDMODEL%YRML_PHY_SLIN, YDMODEL%YRML_PHY_EC, YDMODEL%YRML_GCONF%YGFL, &
    & YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, 0, YDCPG_DIM%KFLEVG, ZDXTDK, ISPPN2D, LLLAND,                    &
    & LLDSLPHY, TSPHY, ZVDIFTS, YLMF_PHYS_BASE_STATE%T, ZQV, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V,             &
    & ZLISUM, ZVERVEL, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,                   &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, ZGEOM1, ZGEOMH, YDVARS%GEOMETRY%GM%T0,          &
    & ZCUCONVCA, ZGP2DSPPA, ZTENT, ZTENQ, ZTENU, ZTENV, ZACPR, ITOPC, IBASC, ITYPE, ICBOT, ICTOP,                        &
    & IBOTSC, LLCUM, LLSC, LLSHCV, ZLCRIT_AER, ZLU, ZLUDE, ZLUDELI, ZSNDE, ZMFU, ZMFD, YDMF_PHYS%OUT%DIFCQ,              &
    & YDMF_PHYS%OUT%DIFCS, ZFHPCL, ZFHPCN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, ZLRAIN, ZRSUD, YDMF_PHYS%OUT%STRCU, &
    & YDMF_PHYS%OUT%STRCV, ZFCQLF, ZFCQLI, ZMFUDE_RATE, ZMFDDE_RATE, YDMF_PHYS%OUT%CAPE, ZWMEAN,                         &
    & ZDIFF, 0, ZCEN, ZTENC, ZSCAV)
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

  ENDIF
  IF(LFLASH .AND. LCVTDK) THEN
    ! LIGHTNING PARAMETERIZATION. OUTPUT IS PFLASH, TOTAL LIGHTNING FLASH RATES.
    ZGAW(:)=0._JPRB
    CALL CULIGHT   ( YDEPHY, YGFL, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,                            &
    & ZGAW, ZGAW, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, &
    & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, LLLAND, YLMF_PHYS_BASE_STATE%T, ZLU, ZMFU, YDMF_PHYS%OUT%CAPE,                          &
    & YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, LLCUM, ICBOT, ICTOP, LLLINOX, YDMF_PHYS%OUT%FLASH,                                &
    & ZLIGH_CTG, ZCTOPH, ZPRECMX, ZICE, ZCDEPTH, ZWMFU) 
    ! LIGHTNING FLASH RATES ARE CONVERTED IN fl/km2/s BEFORE ENTERING CFU TIME ACCUMULATION.
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%FLASH(JLON)=YDMF_PHYS%OUT%FLASH(JLON)/86400._JPRB
    ENDDO
  ENDIF

  IF ( LADJCLD ) THEN
    ZFRSO(:,:) = 0.0_JPRB
    ZFRTH(:,:) = 0.0_JPRB
    LLPTQ = .TRUE.
    CALL ACUPTQ ( YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, LLPTQ, ZFRSO,                   &
    & ZFRTH, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCS, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%FCCQL, &
    & YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN,    &
    & YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YLMF_PHYS_BASE_STATE%T,    &
    & YLMF_PHYS_BASE_STATE%Q, YLMF_PHYS_BASE_STATE%YGSP_RR%T, ZTENHA, ZTENQVA )
  ENDIF

  IF ( LPROCLD ) THEN
    IF(LGPCMT) THEN
      ! Microphysics occurs in a shell between updraft and its resolved environment.
      ZNEBS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)=1._JPRB-(1._JPRB-ZNEBS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG))&
        & *(1._JPRB-ZCSGC(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG))
      ZTMIC(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)=YLMF_PHYS_BASE_STATE%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)
      ZQMIC(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)=ZCSGC(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)&
        & *ZFLU_QSAT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)&
        & +(1._JPRB-ZCSGC(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG))*ZQV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)
      ZQLIS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)=MAX(ZQLIS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG),&
        & ZCSGC(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)&
        & *(ZQL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)+ZQI(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG))&
        & +(1._JPRB-ZCSGC(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG))*ZQLIS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG))
      ZQC_DET_PCMT(:,:)=0._JPRB
    ELSE
      ! Microphysics occurs in the resolved state.
      ZTMIC(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)=YLMF_PHYS_BASE_STATE%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)
      ZQMIC(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)=ZQV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)
    ENDIF

    CALL ACPLUIZ ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTPLUI, YDCPG_DIM%KFLEVG,    &
    & ZTMIC, ZQMIC, ZQL, ZQI, ZQR, ZQS, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, &
    & YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, ZNEBS, ZQLIS,                         &
    & ZNEB_CVPP, ZQLI_CVPP, ZQC_DET_PCMT, ZTENHA, ZTENQVA, LADJCLD, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI,                 &
    & YLMF_PHYS_BASE_STATE%YGSP_RR%T, ZFLU_NEIJ, YDMF_PHYS_SURF%GSD_VF%PLSM, YDVARS%GEOMETRY%GM%T0, ZVETAF, YDMF_PHYS%OUT%FCSQL, &
    & YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN,     &
    & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, ZSEDIQL, ZSEDIQI )

    IF(LGPCMT.AND.LCVNHD) THEN
      ! Evaporation processes within convective environment
      ! will be used by the convection, to compute their feedback on convective updrafts. 
      ! This information is put here in PDDAL.
      YDVARS%DAL%T0(:,:)=0._JPRB
      DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          YDVARS%DAL%T0(JLON,JLEV)=( &
            & +YDMF_PHYS%OUT%FPEVPSL(JLON,JLEV)-YDMF_PHYS%OUT%FPEVPSL(JLON,JLEV-1)   &
            & +YDMF_PHYS%OUT%FPEVPSN(JLON,JLEV)-YDMF_PHYS%OUT%FPEVPSN(JLON,JLEV-1)   &
            & +YDMF_PHYS%OUT%FPEVPCL(JLON,JLEV)-YDMF_PHYS%OUT%FPEVPCL(JLON,JLEV-1)   &
            & +YDMF_PHYS%OUT%FPEVPCN(JLON,JLEV)-YDMF_PHYS%OUT%FPEVPCN(JLON,JLEV-1)) &
            & /RG*YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP(JLON,JLEV)
        ENDDO
      ENDDO
      IF(LMUSCLFA) CALL WRSCMR(NMUSCLFA, 'PFPEVPSL0', YDMF_PHYS%OUT%FPEVPSL, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG+1&
                   &           )
      IF(LMUSCLFA) CALL WRSCMR(NMUSCLFA, 'PFPEVPSN0', YDMF_PHYS%OUT%FPEVPSN, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG+1&
                   &           )
      IF(LMUSCLFA) CALL WRSCMR(NMUSCLFA, 'PFPEVPCL0', YDMF_PHYS%OUT%FPEVPCL, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG+1&
                   &           )
      IF(LMUSCLFA) CALL WRSCMR(NMUSCLFA, 'PFPEVPCN0', YDMF_PHYS%OUT%FPEVPCN, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG+1&
                   &           )
    ENDIF
  ENDIF

  IF ( LNEBCO ) THEN
    CALL ACNEBC ( YDPHY0, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTNEBU, YDCPG_DIM%KFLEVG,                  &
    & YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, INLAB, INND, YDMF_PHYS_SURF%GSD_VH%PTCCH, YDMF_PHYS_SURF%GSD_VH%PSCCH, &
    & YDMF_PHYS_SURF%GSD_VH%PBCCH)
  ENDIF
  IF ( LGWDC ) THEN
    CALL ACDRAC ( YDPHY0, YDPHY2, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTDRAG, YDCPG_DIM%KFLEVG,        &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, ZNBVNO, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YLMF_PHYS_BASE_STATE%U, &
    & YLMF_PHYS_BASE_STATE%V, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS_SURF%GSD_VH%PSCCH,                 &
    & YDMF_PHYS_SURF%GSD_VH%PBCCH, YDMF_PHYS%OUT%STRCU, YDMF_PHYS%OUT%STRCV)
  ENDIF

  ZFLX_LOTT_GWU(:, :) = 0.0_JPRB
  ZFLX_LOTT_GWV(:, :) = 0.0_JPRB

  IF ( LNORGWD ) THEN

    ! Inversion du sens des niveaux
    DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
      DO JLEV = 1, YDCPG_DIM%KFLEVG
        Z_PP(JLON, JLEV) = YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(JLON, YDCPG_DIM%KFLEVG - JLEV + 1)
        Z_UU(JLON, JLEV) = YLMF_PHYS_BASE_STATE%U(JLON, YDCPG_DIM%KFLEVG - JLEV + 1)
        Z_VV(JLON, JLEV) = YLMF_PHYS_BASE_STATE%V(JLON, YDCPG_DIM%KFLEVG - JLEV + 1)
        Z_TT(JLON, JLEV) = YLMF_PHYS_BASE_STATE%T(JLON, YDCPG_DIM%KFLEVG - JLEV + 1)
        Z_VO(JLON, JLEV) = YLMF_PHYS_BASE_STATE%VOR(JLON, YDCPG_DIM%KFLEVG - JLEV + 1)
        ZD_U(JLON, JLEV) = TSPHY * YLMF_PHYS_BASE_STATE%P1NOGW(JLON, YDCPG_DIM%KFLEVG - JLEV + 1)
        ZD_V(JLON, JLEV) = TSPHY * YLMF_PHYS_BASE_STATE%P2NOGW(JLON, YDCPG_DIM%KFLEVG - JLEV + 1)
      ENDDO
    ENDDO

    ZPRECGWD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA) = MAX(0.0_JPRB,                  &
      & YDMF_PHYS%OUT%FPLCL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG) + YDMF_PHYS%OUT%FPLCN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG))
    
    CALL ACNORGWD(YDNORGWD, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, &
    & TSPHY, Z_PP, YDVARS%GEOMETRY%GEMU%T0, Z_TT, Z_UU, Z_VV, Z_VO, ZPRECGWD, ZD_U, ZD_V)
    
    DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
      DO JLEV = 1, YDCPG_DIM%KFLEVG
        YLMF_PHYS_BASE_STATE%P1NOGW(JLON, JLEV) = ZD_U(JLON, YDCPG_DIM%KFLEVG - JLEV + 1) / TSPHY
        YLMF_PHYS_BASE_STATE%P2NOGW(JLON, JLEV) = ZD_V(JLON, YDCPG_DIM%KFLEVG - JLEV + 1) / TSPHY
      ENDDO
    ENDDO

    !-- CALCUL DU FLUX, PAR INTEGRATION DE LA TENDANCE DE HAUT EN BAS.
    !   LES FLUX SONT SUPPOSES NULS AU PREMIER NIVEAU (1 = TOP) DE CALCUL.
    DO JLEV = 1, YDCPG_DIM%KFLEVG
       DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
          ZFLX_LOTT_GWU(JLON, JLEV) = ZFLX_LOTT_GWU(JLON, JLEV - 1) - YLMF_PHYS_BASE_STATE%P1NOGW(JLON, JLEV)*YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON, JLEV)/RG
          ZFLX_LOTT_GWV(JLON, JLEV) = ZFLX_LOTT_GWV(JLON, JLEV - 1) - YLMF_PHYS_BASE_STATE%P2NOGW(JLON, JLEV)*YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP(JLON, JLEV)/RG
          YDMF_PHYS%OUT%STRCU(JLON, JLEV) = YDMF_PHYS%OUT%STRCU(JLON, JLEV) + ZFLX_LOTT_GWU(JLON, JLEV)
          YDMF_PHYS%OUT%STRCV(JLON, JLEV) = YDMF_PHYS%OUT%STRCV(JLON, JLEV) + ZFLX_LOTT_GWV(JLON, JLEV)          
    !+ TODOLATER +!      PSTRNORGWDU(JLON, JLEV) = PSTRNORGWDU(JLON, JLEV) + ZFLX_LOTT_GWU(JLON, JLEV)
    !+ TODOLATER +!      PSTRNORGWDV(JLON, JLEV) = PSTRNORGWDV(JLON, JLEV) + ZFLX_LOTT_GWV(JLON, JLEV)          
       ENDDO
    ENDDO

    ! NO VERTICAL DIFFUSION IN THE MIDDLE ATMOSPHERE
    DO JLEV = 1, NORGWD_NNOVERDIF
       DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
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
        DO JLEV=YDCPG_DIM%KTDIA-1,YDCPG_DIM%KFLEVG
          DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            ZPFL_FPLCH(JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)+YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
          ENDDO
        ENDDO
      ELSE
        DO JLEV=YDCPG_DIM%KTDIA-1,YDCPG_DIM%KFLEVG
          DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            ZPFL_FPLCH(JLON,JLEV)=ZFPCOR(JLON,JLEV)
          ENDDO
        ENDDO
      ENDIF
    ELSE
      DO JLEV=YDCPG_DIM%KTDIA-1,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZPFL_FPLCH(JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)+YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
        ENDDO
      ENDDO
    ENDIF
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
  IF (LGRAPRO) THEN
    DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZFPLSN (JLON,JLEV)=ZFPLSN(JLON,JLEV)+YDMF_PHYS%OUT%FPLSG(JLON,JLEV)
      ENDDO
    ENDDO
  ENDIF

  ! - - - - - - - - - - - - - - - - -
  ! CORRECT NEGATIVE WATER CONTENTS.
  ! - - - - - - - - - - - - - - - - -

  IF ( LCONDWT.AND.LPROCLD.AND..NOT.LGPCMT ) THEN
    CALL QNGCOR (YDPHY2, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTPLUI, YDCPG_DIM%KFLEVG,                   &
    & ZQV, ZQL, ZQI, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN,               &
    & YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%FPEVPSL,    &
    & YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, &
    & YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL,       &
    & YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FCQNG )

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

    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDCPG_GPAR%RAIN(JLON)=ZFPLSL(JLON,YDCPG_DIM%KFLEVG)+YDMF_PHYS%OUT%FPLCL(JLON,YDCPG_DIM%KFLEVG)
      YDCPG_GPAR%SNOW(JLON)=ZFPLSN(JLON,YDCPG_DIM%KFLEVG)+YDMF_PHYS%OUT%FPLCN(JLON,YDCPG_DIM%KFLEVG)
    ENDDO

  ELSE

    IF ( LSFHYD.AND.LSOLV ) THEN
      CALL ACDROV ( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,                  &
      & YSP_SBD%NLEVS, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, ZFPLSL, ZFPLSN, YDMF_PHYS%OUT%FRSO,                          &
      & YDMF_PHYS%OUT%FRTH, ZDSA_C1, ZDSA_C2, ZC3, ZCN, YDMF_PHYS%OUT%CT, YDMF_PHYS_SURF%GSD_VV%PD2,                          &
      & YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FTR, YDMF_PHYS_SURF%GSD_VV%PLAI, ZFLU_NEIJ, ZFLU_VEG, ZWFC,                         &
      & ZWPMX, YLMF_PHYS_BASE_STATE%YGSP_RR%FC, ZWLMX, ZWSEQ, ZWSMX, YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FCLL,                 &
      & YDMF_PHYS%OUT%FCLN, YDMF_PHYS%OUT%FCS, ZFLU_FEVI, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, YDMF_PHYS_SURF%GSD_VF%PLSM, &
      & YLMF_PHYS_BASE_STATE%YGSP_SG%F, YLMF_PHYS_BASE_STATE%YGSP_SB%T, YLMF_PHYS_BASE_STATE%YGSP_RR%T,                       &
      & YLMF_PHYS_BASE_STATE%YGSP_SB%Q, YLMF_PHYS_BASE_STATE%YGSP_SB%TL, YLMF_PHYS_BASE_STATE%YGSP_RR%W,                      &
      & YLMF_PHYS_BASE_STATE%YGSP_RR%IC, YDMF_PHYS%OUT%FGEL, YDMF_PHYS%OUT%FGELS, YDMF_PHYS%OUT%FLWSP,                        &
      & YDMF_PHYS%OUT%FONTE, YDMF_PHYS%OUT%RUISP, YDMF_PHYS%OUT%RUISL, YDMF_PHYS%OUT%RUISS)
    ENDIF

  ENDIF

!*
!-  --------------------------------------------------------------------
!     13.- DRAG MESOSPHERIQUE POUR UN MODELE POSSEDANT DES NIVEAUX
!               AU-DESSUS DE 50 KM (I.E. DANS LA MESOSPHERE)
!     ------------------------------------------------------------------
  IF ( LRRMES ) THEN
    CALL ACDRME ( YDPHY2, YDTOPH, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTDRME, YDCPG_DIM%KFLEVG,  &
    & YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_BASE_STATE%T,    &
    & ZQV, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YDMF_PHYS%OUT%FRMH, ZMSC_FRMQ, YDMF_PHYS%OUT%STRMU, &
    & YDMF_PHYS%OUT%STRMV, RMESOT, RMESOQ, RMESOU, STTEM)
  ENDIF

!*
!     ------------------------------------------------------------------
!     14.- FLUX PHOTO-CHIMIQUE D'OZONE
!     ------------------------------------------------------------------
  IF ( LOZONE ) THEN
    CALL ACOZONE ( YDPHY2, YDTOPH, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTOZON, YDCPG_DIM%KFLEVG,       &
    & YDMODEL%YRML_PHY_G%YRDPHY%NVCLIS, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, &
    & YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDCPG_MISC%KOZO, ZQO3(1, 1), YLMF_PHYS_BASE_STATE%T,                   &
    & ZRDG_MU0, YDMF_PHYS%OUT%FCHOZ)

!           DIFFUSION TURBULENTE/DEPOT SEC DE L'OZONE

    CALL ACDIFOZ ( YDRIP, YDPHY2, YDMODEL%YRML_PHY_MF%YRVDOZ, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,                  &
    & YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP,          &
    & YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FRSO, ZKTROV, ZQO3(1, 1), YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, &
    & YLMF_PHYS_BASE_STATE%T, ZXTROV, ZFLU_NEIJ, YDVARS%GEOMETRY%GEMU%T0, YDMF_PHYS_SURF%GSD_VV%PIVEG, YDMF_PHYS%OUT%FCHOZ                 &
    &                                                                                    )
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
    DO JLON=YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
       YDMF_PHYS_SURF%GSD_VH%PSPSH(JLON)=(YLMF_PHYS_BASE_STATE%T(JLON,YDCPG_DIM%KFLEVG)-YLMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)) * ((&
          & RCW-RCPD)*(YDMF_PHYS%OUT%FPLSL(JLON,YDCPG_DIM%KFLEVG)+YDMF_PHYS%OUT%FPLCL(JLON,YDCPG_DIM%KFLEVG))&
          & +(RCS-RCPD)*(YDMF_PHYS%OUT%FPLSN(JLON,YDCPG_DIM%KFLEVG)+YDMF_PHYS%OUT%FPLCN(JLON,YDCPG_DIM%KFLEVG)) )
    ENDDO
    IF (LGRAPRO) THEN
      YDMF_PHYS_SURF%GSD_VH%PSPSH(JLON)=YDMF_PHYS_SURF%GSD_VH%PSPSH(JLON)+(YLMF_PHYS_BASE_STATE%T(JLON,YDCPG_DIM%KFLEVG)-YLMF_PHYS_BASE_STATE%YGSP_RR%T(JLON)) * (&
        &(RCS-RCPD)*(YDMF_PHYS%OUT%FPLSG(JLON,YDCPG_DIM%KFLEVG)+YDMF_PHYS%OUT%FPLCG(JLON,YDCPG_DIM%KFLEVG)))
    ENDIF
  ENDIF !LPHSPSH

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
    & YLMF_PHYS_BASE_STATE)
  ENDIF

  IF (LECT.AND.LRAFTKE) THEN
     ! DCAPE due to precipitation evaporation.
     CALL ACEVADCAPE(YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YDMF_PHYS%OUT%FPLSL, &
     & YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD,   &
     & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP,   &
     & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, ZDCAPE)
     ! Gusts.
     ZQGM=ZEPSNEB
     CALL ACCLDIA(YDXFU, YDPHY, YDPHY2, YDTOPH, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, &
     & YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YDMF_PHYS%OUT%CAPE,  &
     & ZDCAPE, ZTKE1, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDVARS%GEOMETRY%OROG%T0, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST,      &
     & ZBLH, IMOC_CLPH)
     YDMF_PHYS%OUT%CLPH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=MIN(XMAXLM,MAX(XMINLM,ZBLH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)))
  ENDIF

  IF (LXVISI.OR.LXVISI2) THEN
     CALL ACVISIH(YDMODEL%YRML_PHY_MF%YRPHY%YRDVISI, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,           &
     & YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, &
     & YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,       &
     & YDCPG_MISC%QLI, YDCPG_MISC%QICE, ZQR, ZQS, ZQGM, YDMF_PHYS%OUT%VISICLD, YDMF_PHYS%OUT%VISIHYD,            &
     & YDMF_PHYS%OUT%MXCLWC    )
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

  DO JLEV = 1 , YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZZDEP(JLON,1,JLEV)=ZZI_APHI(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO


  DO JLEV=1,YDCPG_DIM%KFLEVG
    DO JLON= YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
      ZQDM(JLON,JLEV)=1._JPRB-ZQV(JLON,JLEV)-ZQL(JLON,JLEV)-ZQR(JLON,JLEV) &
       & -ZQI(JLON,JLEV)-ZQS(JLON,JLEV)
    ENDDO
  ENDDO

  DO JLEV = 1, YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZZI_RM(JLON,1,JLEV,2)=ZQL(JLON,JLEV)&
       & /ZQDM(JLON,JLEV)
    ENDDO
  ENDDO

  DO JLEV = 1, YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZZI_RM(JLON,1,JLEV,3)=ZQR(JLON,JLEV)&
       & /ZQDM(JLON,JLEV)
    ENDDO
  ENDDO

  DO JLEV = 1, YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZEVAP(JLON,1,JLEV)=YDMF_PHYS%OUT%FPEVPSL(JLON,JLEV)&
       & +YDMF_PHYS%OUT%FPEVPCL(JLON,JLEV)
    ENDDO
  ENDDO

  DO JGFL=1,NGFL_EXT
    DO JLON=1,YDCPG_DIM%KLON
      DO JLEV=1,YDCPG_DIM%KFLEVG
        ZZI_SVM(JLON,1,JLEV,JGFL)=MAX(0._JPRB,ZSVM(JLON,JLEV,JGFL))
      ENDDO
    ENDDO
  ENDDO

    CALL ARO_WETDEP(ILONMNH, YDCPG_DIM%KFLEVG, NGFL_EXT, IKRR, PDTPHY, ZZI_SVM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :, 1:NGFL_EXT), &
    & ZZDEP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :), ZZI_PABSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :),                              &
    & ZZI_THM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :), ZZI_RHODREFM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :),                         &
    & YDCPG_DIM%KSTEP+1, ZZI_RM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :, :), ZEVAP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :, :),           &
    & ISPLITR               )
! return to tendency
   DO JGFL=1,NGFL_EXT
     DO JLON=1,YDCPG_DIM%KLON
       DO JLEV=1,YDCPG_DIM%KFLEVG
         ZPSV(JLON,JLEV,JGFL)=ZZI_SVM(JLON,1,JLEV,JGFL)
       ENDDO
     ENDDO
   ENDDO

      DO JGFL=1,NGFL_EXT
        DO JLON=1,YDCPG_DIM%KLON
          DO JLEV=1,YDCPG_DIM%KFLEVG
            ZTENDEXT_DEP(JLON,JLEV,JGFL)=(ZPSV(JLON,JLEV,JGFL)-ZSVM(JLON,JLEV,JGFL))*ZINVDT
          ENDDO
        ENDDO
      ENDDO

ENDIF ! ENDIF  (LMDUST & NGFL_EXT & LRDEPOS)

IF(LMUSCLFA) THEN
  DO JLEV=0,YDCPG_DIM%KFLEVG
    DO JLON=1,YDCPG_DIM%KLON
      ZFPLS(JLON,JLEV)=YDMF_PHYS%OUT%FPLSL(JLON,JLEV)+YDMF_PHYS%OUT%FPLSN(JLON,JLEV)
      ZFPLC(JLON,JLEV)=YDMF_PHYS%OUT%FPLCL(JLON,JLEV)+YDMF_PHYS%OUT%FPLCN(JLON,JLEV)
      ZFPL (JLON,JLEV)=ZFPLC (JLON,JLEV)+ZFPLS (JLON,JLEV)
    ENDDO
  ENDDO
  CALL WRSCMR(NMUSCLFA, 'ZFPLS', ZFPLS, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG+1)
  CALL WRSCMR(NMUSCLFA, 'ZFPLC', ZFPLC, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG+1)
  CALL WRSCMR(NMUSCLFA, 'ZFPL', ZFPL, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG+1)
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

   ZDE2MR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)  = RD * YLMF_PHYS_BASE_STATE%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) / YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
   ZBAY_QRCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) = ZBAY_QRCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) * ZDE2MR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
   ZBAY_QSCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) = ZBAY_QSCONV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) * ZDE2MR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)

ENDIF


! Precipitation Type

! Compute wet-bulb temperature at 2 meters (suppose homogeneity of qv/ql/qi )
!ZPCLS(KIDIA:KFDIA)=PAPRS(KIDIA:KFDIA,KLEV)-2._JPRB/ZZZF(KIDIA:KFDIA,1,KLEV)*&
!                 &(PAPRS(KIDIA:KFDIA,KLEV)-PAPRSF(KIDIA:KFDIA,KLEV))
CALL PPWETPOINT(YDPHY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, ZPCLS, YDMF_PHYS%OUT%TCLS, YLMF_PHYS_BASE_STATE%Q(:, YDCPG_DIM%KFLEVG), &
& YLMF_PHYS_BASE_STATE%L(:, YDCPG_DIM%KFLEVG), YLMF_PHYS_BASE_STATE%I(:, YDCPG_DIM%KFLEVG), YDMF_PHYS%OUT%TPWCLS                                 &
&                                                                                                                                                                                     )


IF (LDPRECIPS .OR. LDPRECIPS2) THEN
   ! Defined precipitation type
   !
   ! Compute wet-bulb temperature
  DO JLEV=1,YDCPG_DIM%KFLEVG
    CALL PPWETPOINT(YDPHY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF(:, JLEV), &
    & YLMF_PHYS_BASE_STATE%T(:, JLEV), YLMF_PHYS_BASE_STATE%Q(:, JLEV), YLMF_PHYS_BASE_STATE%L(:, JLEV),                     &
    & YLMF_PHYS_BASE_STATE%I(:, JLEV), ZTW(:, JLEV))
  ENDDO

  DO JLON=1,YDCPG_DIM%KLON
      ZFPLS(JLON,YDCPG_DIM%KFLEVG)=YDMF_PHYS%OUT%FPLCN(JLON,YDCPG_DIM%KFLEVG)+YDMF_PHYS%OUT%FPLSN(JLON,YDCPG_DIM%KFLEVG)
      ZFPLC(JLON,YDCPG_DIM%KFLEVG)=YDMF_PHYS%OUT%FPLCL(JLON,YDCPG_DIM%KFLEVG)+YDMF_PHYS%OUT%FPLSL(JLON,YDCPG_DIM%KFLEVG)
      ZFPLSG(JLON,YDCPG_DIM%KFLEVG)=0._JPRB
  ENDDO

  !initialisation de ZZZ
  DO JLEV = 1,YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZZZ(JLON,1,JLEV)=YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO

  !initialisation de ZDZZ
  DO JLEV = 2, YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZDZZ(JLON,1,JLEV)=ZZZ(JLON,1,JLEV-1)-ZZZ(JLON,1,JLEV)
    ENDDO
  ENDDO
  DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZDZZ(JLON,1,1)=YLMF_PHYS_BASE_STATE%YCPG_DYN%PHI(JLON,0)*ZINVG-ZZZ(JLON,1,1)
  ENDDO

  IF (LDPRECIPS) THEN

   NDTPRECCUR=INT(MOD(ZSTATI/TSTEP,REAL(NDTPREC)))+1_JPIM

   CALL DPRECIPS (YDPHY%YRDPRECIPS, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,                      &
   & YDVARS%GEOMETRY%OROG%T0, YDMF_PHYS%OUT%TPWCLS, YDMF_PHYS%OUT%DIAGH, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,                             &
   & ZDZZ, ZTW, YLMF_PHYS_BASE_STATE%L, ZFPLC(:, YDCPG_DIM%KFLEVG), ZFPLS(:, YDCPG_DIM%KFLEVG), ZFPLSG(:, YDCPG_DIM%KFLEVG), &
   & ZPRC_DPRECIPS(:, NDTPRECCUR))

  ! WRITE(20,*)'sous aplpar ZDPRECIPS ', ZDPRECIPS(KIDIA:KFDIA,NDTPRECCUR),NDTPRECCUR

  ENDIF

  IF (LDPRECIPS2) THEN

 !Idem for an other time step and an other period
   NDTPRECCUR2=INT(MOD(ZSTATI/TSTEP,REAL(NDTPREC2)))+1_JPIM

   CALL DPRECIPS(YDPHY%YRDPRECIPS, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,                       &
   & YDVARS%GEOMETRY%OROG%T0, YDMF_PHYS%OUT%TPWCLS, YDMF_PHYS%OUT%DIAGH, YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF,                             &
   & ZDZZ, ZTW, YLMF_PHYS_BASE_STATE%L, ZFPLC(:, YDCPG_DIM%KFLEVG), ZFPLS(:, YDCPG_DIM%KFLEVG), ZFPLSG(:, YDCPG_DIM%KFLEVG), &
   & ZPRC_DPRECIPS2(:, NDTPRECCUR2))

  ! WRITE(20,*)'sous aplpar ZDPRECIPS2',ZDPRECIPS2(KIDIA:KFDIA,NDTPRECCUR2),NDTPRECCUR2

  ENDIF
ENDIF  

IF (LAJUCV) THEN
  DO JLEV=1,YDCPG_DIM%KFLEVG
    DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDVARS%T%T0(JROF,JLEV)=ZADJ_TAUX(JROF,JLEV)
    ENDDO
  ENDDO
  DO JLEV=1,YDCPG_DIM%KFLEVG
    DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YLMF_PHYS_NEXT_STATE%T (JROF, JLEV) = YLMF_PHYS_NEXT_STATE%T (JROF, JLEV) + ZADJ_DTAJU(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

!    convert to flexible interface structure
IF (LINTFLEX) THEN
  CALL APLPAR2INTFLEX(YGFL, YDPHY, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG,                            &
  & YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%     DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS,                                    &
  & ZDIFEXT, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT    %DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTS,                            &
  & YDMF_PHYS%OUT%FCCQL, YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLSL,                      &
  & YDMF_PHYS%OUT%FPLSN, YDMF_PHYS%    OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPEVPSL, YDMF_PHYS%OUT%FPEVPSN,              &
  & YDMF_PHYS%OUT%      FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%        FPFPCL, &
  & YDMF_PHYS%OUT%FPFPCN, YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG,                 &
  & YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%FRMH, ZMSC_FRMQ, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%STRCU,              &
  & YDMF_PHYS%OUT%STRCV, YDMF_PHYS%      OUT%STRDU, YDMF_PHYS%OUT%STRDV, YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV,                &
  & YDMF_PHYS%OUT%STRMU, YDMF_PHYS%OUT%STRMV, YDMF_PHYS%OUT%DIFCQLC, YDMF_PHYS%OUT%DIFCQIC, YDMF_PHYS%OUT%FIMCC,                  &
  & YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC,               &
  & YDMF_PHYS%OUT%FCNEGQIC, YDMF_PHYS%OUT%FCNEGQRC, YDMF_PHYS%OUT%FCNEGQSC, ZPFL_FTKE, ZTENDPTKE,                                 &
  & ZTENDEXT, ZTENDEXT_DEP, YLPROCSET )
ENDIF

!        2.3  Computes MOCON in the CLP.
!             --------------------------
CALL MF_PHYS_MOCON (YDCPG_DIM, YDMF_PHYS%OUT%MOCON, ZRDG_LCVQ, IMOC_CLPH, YLMF_PHYS_BASE_STATE)

! Store surface water flux P and E for water conservation
IF (LCORWAT) THEN
  CALL MF_PHYS_CORWAT (YDCPG_DIM, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, &
  & YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN, YDMF_PHYS_SURF%GSD_VH%PEVA, YDMF_PHYS_SURF%GSD_VH%PPCL,               &
  & YDMF_PHYS_SURF%GSD_VH%PPCN, YDMF_PHYS_SURF%GSD_VH%PPSL, YDMF_PHYS_SURF%GSD_VH%PPSN)
ENDIF

!        2.4  Stores radiation coefficients.
!             ------------------------------

! * writes grid-point transmission coefficients for simplified physics.


IF (LRCOEF.AND.(YDCPG_DIM%KSTEP == 1)) THEN
    IFIELDSS=NG3SR*YDCPG_DIM%KFLEVG
    CALL WRRADCOEF(YDGEOMETRY, YDRCOEF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KSTGLO, IFIELDSS, &
    & ZRDT_COR, ZRDT_RAB3C, ZRDT_RAB3N, ZRDT_RAB4C, ZRDT_RAB4N, ZRDT_RAB6C, ZRDT_RAB6N, ZRDT_RAT1C,   &
    & ZRDT_RAT1N, ZRDT_RAT2C, ZRDT_RAT2N, ZRDT_RAT3C, ZRDT_RAT3N, ZRDT_RAT4C, ZRDT_RAT4N, ZRDT_RAT5C, &
    & ZRDT_RAT5N, ZAC_HC)
ENDIF

!       2.5   Ozone
!             -----


IF (LOZONE) THEN
  ! * Caution: this part has not been yet validated relative
  !   to the GFL implementation, and LOZONE (the setup of
  !   which has not yet been updated) can be true only if
  !   the GFL ozone is activated as a prognostic and advected
  !   variable.
  CALL CPOZO (YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, PDTPHY, YDMF_PHYS%OUT%FCHOZ, &
  & YLMF_PHYS_NEXT_STATE%O3 (:, 1:YDCPG_DIM%KFLEVG), YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP)
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

IF (NDPSFI == 1) THEN
  CALL CPQSOL(YDGEOMETRY%YRDIMV, YDPHY, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_PHY0%PREHYD, &
  & YDMF_PHYS_SURF%GSP_RR%PT_T0, YDCPG_MISC%QS, ZFLU_QSATS, YDCPG_MISC%QSOL)
ENDIF

!        2.7  Computation of tendencies T,u,v and Q.
!             --------------------------------------

! Set GFL tendencies to 0

ZTENDGFL(:,:,:) = 0.0_JPRB

TSPHY = MAX(PDTPHY,1.0_JPRB)

! * CPTEND+CPUTQY = Old( CPATY + CPDUP + CPDTHP )
! Calcul des tendances de T , U et de Q et modifications
! eventuelles de W et de OMEGA/P

IF (LINTFLEX.AND.(.NOT.LDCONFX)) THEN
  CALL CPTEND_FLEX( YDLDDH, YDMDDH, YGFL, YDPHY, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG,   &
  & YDVARS%GEOMETRY%GNORDL%T0, YDVARS%GEOMETRY%GNORDM%T0, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP, &
  & YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%U, YLMF_PHYS_BASE_STATE%V, YLMF_PHYS_BASE_STATE%T,      &
  & YLMF_PHYS_BASE_STATE%YGSP_RR%T, PGFL, YLPROCSET, YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV,                         &
  & ZTENDH, ZTENDGFL, YDMF_PHYS%OUT%FHSCL, YDMF_PHYS%OUT%FHSCN, YDMF_PHYS%OUT%FHSSL, YDMF_PHYS%OUT%FHSSN,              &
  & YDMF_PHYS%OUT%FHPCL, YDMF_PHYS%OUT%FHPCN, YDMF_PHYS%OUT%FHPSL, YDMF_PHYS%OUT%FHPSN, PFHP=ZMSC_FHP,                 &
  & PFP=ZPFL_FP, PFEPFP=YDMF_PHYS%OUT%FEPFP, PFCMPCQ=YDMF_PHYS%OUT%FCMPCQ, PFCMPSN=YDMF_PHYS%OUT%FCMPSN,               &
  & PFCMPSL=YDMF_PHYS%OUT%FCMPSL, YDDDH=YDDDH)
ELSE
  CALL CPTEND_NEW( YDMODEL, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, YDVARS%GEOMETRY%GNORDL%T0,            &
  & YDVARS%GEOMETRY%GNORDM%T0, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS,                  &
  & ZDIFEXT, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%FCCQL,     &
  & YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN,                &
  & YDMF_PHYS%OUT%FPLSG, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLCG, YDMF_PHYS%OUT%FPEVPSL,              &
  & YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPEVPCG,      &
  & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN,           &
  & YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FCQGNG,           &
  & YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%FRMH, ZMSC_FRMQ, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%STRCU,        &
  & YDMF_PHYS%OUT%STRCV, YDMF_PHYS%OUT%STRDU, YDMF_PHYS%OUT%STRDV, YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV,                &
  & YDMF_PHYS%OUT%STRMU, YDMF_PHYS%OUT%STRMV, YDMF_PHYS%OUT%DIFCQLC, YDMF_PHYS%OUT%DIFCQIC, YDMF_PHYS%OUT%FIMCC,            &
  & YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC,         &
  & YDMF_PHYS%OUT%FCNEGQIC, YDMF_PHYS%OUT%FCNEGQRC, YDMF_PHYS%OUT%FCNEGQSC, ZPFL_FTKE, ZPFL_FTKEI,                          &
  & ZPFL_FEFB1, ZPFL_FEFB2, ZPFL_FEFB3, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,    &
  & YLMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YLMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_BASE_STATE%U,                       &
  & YLMF_PHYS_BASE_STATE%V, YLMF_PHYS_BASE_STATE%T, YLMF_PHYS_BASE_STATE%Q, YLMF_PHYS_BASE_STATE%I, YLMF_PHYS_BASE_STATE%L, &
  & YDVARS%LCONV%T0, YDVARS%ICONV%T0, YDVARS%RCONV%T0, YDVARS%SCONV%T0, YLMF_PHYS_BASE_STATE%R, YLMF_PHYS_BASE_STATE%S,     &
  & YLMF_PHYS_BASE_STATE%G, ZDSA_CPS, YLMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS%OUT%FHSCL, YDMF_PHYS%OUT%FHSCN,             &
  & YDMF_PHYS%OUT%FHSSL, YDMF_PHYS%OUT%FHSSN, YDMF_PHYS%OUT%FHSSG, YDMF_PHYS%OUT%FHPCL, YDMF_PHYS%OUT%FHPCN,                &
  & YDMF_PHYS%OUT%FHPCG, YDMF_PHYS%OUT%FHPSL, YDMF_PHYS%OUT%FHPSN, YDMF_PHYS%OUT%FHPSG, ZMSC_FHP,                           &
  & ZPFL_FP, YDMF_PHYS%OUT%FEPFP, YDMF_PHYS%OUT%FCMPCQ, YDMF_PHYS%OUT%FCMPSN, YDMF_PHYS%OUT%FCMPSL, YDMF_PHYS%OUT%TENDU,    &
  & YDMF_PHYS%OUT%TENDV, ZTENDU, ZTENDV, ZTENDH, ZPTENDQ1, ZPTENDI1, ZPTENDL1, ZPTENDLCONV1,                                &
  & ZPTENDICONV1, ZPTENDRCONV1, ZPTENDSCONV1, ZPTENDR1, ZPTENDS1, ZPTENDG1, ZPTENDTKE1, ZPTENDEFB11,                        &
  & ZPTENDEFB21, ZPTENDEFB31, ZTENDEXT, YDDDH)
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

ENDIF 



!        2.7.1  Diagnostics on physical tendencies
!               ----------------------------------

IF (.NOT.LDCONFX) THEN
  IF ((GCHETN%LFREQD).OR.(GCHETN%LCOORD).OR.(GCHETN%LPROFV)) THEN
    CALL CPCHET (YDMF_PHYS, YLMF_PHYS_BASE_STATE, YDCPG_MISC, YDRIP, YDPHY, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, &
    & YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, YDCPG_DIM%KSTEP, ZMSC_FRMQ, ZDSA_CPS, ZTENDH, ZPTENDQ1,             &
    & ZPTENDI1, ZPTENDL1, ZPTENDR1, ZPTENDS1, YDVARS%GEOMETRY%GEMU%T0, YDVARS%GEOMETRY%GELAM%T0)
  ENDIF
  
  IF (GCHETN%LPROFV)&
   & CALL PROFILECHET(YDGEOMETRY, YDCPG_MISC, YDMF_PHYS, ZRDG_CVGQ, ZRDG_LCVQ, ZRDG_MU0, YDCPG_DYN0,        &
     & YDMF_PHYS_SURF, YDVARS, YDSURF, YDDPHY, YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KFDIA, YDVARS%GEOMETRY%GELAM%T0, &
     & YDVARS%GEOMETRY%GEMU%T0, YDVARS%GEOMETRY%GM%T0, YDVARS%GEOMETRY%OROG%T0, YDVARS%GEOMETRY%RCORI%T0, YDVARS%GEOMETRY%RATATH%T0, YDVARS%GEOMETRY%RATATX%T0)

ENDIF

!        2.8  Modification of vertical velocities
!             by some physics output when required.
!             -------------------------------------


! * MODIFICATION DE LA VITESSE VERTICALE ET DE LA TENDANCE DE
! PRESSION DE SURFACE SI NDPSFI=1 ( MASSE VARIABLE ).
! Ajout de la physique dans l'equation de continuite/Add physics
! in continuity equation.

IF (NDPSFI == 1) THEN
  CALL CPMVVPS(YDVAB, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, PDTPHY,               &
  & ZPFL_FP, YLMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD(:, YDCPG_DIM%KFLEVG), YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT%FEVN, &
  & YDCPG_DYN0%CTY%EVEL, YDCPG_DYN0%CTY%PSDVBC, YLMF_PHYS_NEXT_STATE%SP)
ENDIF

!        2.9  Computation of evolution of T, u, v and Q.
!             ------------------------------------------

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

! Calcul de T , Q et du Vent a l'instant 1

CALL CPUTQY_APLPAR_EXPL(YDCPG_DIM, YLMF_PHYS_NEXT_STATE, YLMF_PHYS_BASE_STATE, YDVARS, YDPHY, PDTPHY, &
& ZTENDH, ZTENDT, YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDU, ZTENDV, ZTENDD, ZPTENDEFB11,      &
& ZPTENDEFB21, ZPTENDEFB31, ZPTENDG1, ZPTENDICONV1, ZPTENDI1, ZPTENDLCONV1, ZPTENDL1, ZPTENDQ1,       &
& ZPTENDRCONV1, ZPTENDR1, ZPTENDSCONV1, ZPTENDS1, ZPTENDTKE1, YDMF_PHYS%OUT%FDIS)

CALL CPUTQY_APLPAR_LOOP(YDMODEL%YRML_DYN%YRDYN, YDGEOMETRY%YRDIMV, YDGMV, YGFL, YDPTRSLB1, YDPHY, YDCPG_DIM%KLON, &
& YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, PDTPHY, ZTENDGFL, YDCPG_SL1%ZVIEW, PGMVT1,                  &
& PGFLT1)

CALL MF_PHYS_FPL_PART2 (YDCPG_DIM, ZPFL_FPLCH, ZPFL_FPLSH, YDVARS%CPF%T1, YDVARS%SPF%T1, YDMODEL)

!       2.9b Prognostic convection etc.
!            --------------------------

! TRANSFER NOT ADVECTED VARIABLES INTO PGFLT1
CALL MF_PHYS_TRANSFER (YDCPG_DIM, YDVARS%DAL%T0, YDVARS%DAL%T1, YDVARS%DOM%T0, YDVARS%DOM%T1, YDVARS%FQTUR%T0,     &
& YDVARS%FQTUR%T1, YDVARS%FSTUR%T0, YDVARS%FSTUR%T1, YDVARS%MXL%T0, YDVARS%MXL%T1, YDVARS%RKTH%T0, YDVARS%RKTH%T1, &
& YDVARS%RKTQC%T0, YDVARS%RKTQC%T1, YDVARS%RKTQV%T0, YDVARS%RKTQV%T1, YDVARS%SHTUR%T0, YDVARS%SHTUR%T1,            &
& YDVARS%TTE%T0, YDVARS%TTE%T1, YDVARS%UAL%T0, YDVARS%UAL%T1, YDVARS%UEN%T0, YDVARS%UEN%T1, YDVARS%UNEBH%T0,       &
& YDVARS%UNEBH%T1, YDVARS%UOM%T0, YDVARS%UOM%T1, YDMODEL)

!        2.10  Surface variables.
!              ------------------

IF ((.NOT.LSFORCS)) THEN
  
  IF (.NOT.LMSE) THEN
    DO JLEV=0,YDCPG_DIM%KFLEVG
      DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZPFL_FPLSN(JROF,JLEV)=YDMF_PHYS%OUT%FPLSN(JROF,JLEV)+YDMF_PHYS%OUT%FPLSG(JROF,JLEV)
      ENDDO
    ENDDO
    CALL CPTENDS( YDMODEL%YRML_PHY_MF, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG,                        &
    & YSP_SBD%NLEVS, PDTPHY, YDMF_PHYS%   OUT%FPLCL, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLCN,                                    &
    & ZPFL_FPLSN, YDMF_PHYS  %OUT%FRSO, YDMF_PHYS%OUT%FRTH, YDMF_PHYS_SURF%GSP_SG%PA_T1, YDMF_PHYS%OUT%CT,                        &
    & ZDSA_C1, ZDSA_C2, YDMF_PHYS%OUT%FCHSP, YDMF_PHYS%OUT%FCLL, YDMF_PHYS  %OUT%FCLN, YDMF_PHYS%OUT%FCS,                         &
    & ZFLU_FEVI, YDMF_PHYS%OUT%FEVL, YDMF_PHYS%OUT  %FEVN, YDMF_PHYS%OUT%FEVV, YDMF_PHYS%OUT%FGEL, YDMF_PHYS%OUT%FGELS,           &
    & YDMF_PHYS%OUT%FLWSP, YDMF_PHYS%OUT%FONTE, YDMF_PHYS%OUT%FTR, YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSP_SG%      PR_T1, &
    & YDMF_PHYS%OUT%RUISL, YDMF_PHYS%OUT%RUISP, YDMF_PHYS%OUT%RUISS, YDMF_PHYS_SURF%GSP_SG%    PF_T1,                             &
    & ZFLU_VEG, ZTDS_TDTS, ZTDS_TDTP, ZTDS_TDWS, ZTDS_TDWSI, ZTDS_TDWP, ZTDS_TDWPI, ZTDS_TDWL,                                    &
    & ZTDS_TDSNS, ZTDS_TDALBNS, ZTDS_TDRHONS)  

    CALL CPWTS(YDSURF, YDMODEL%YRML_AOC%YRMCC, YDPHY, YDMODEL%YRML_PHY_MF%YRPHY1, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA,           &
    & YDCPG_DIM%KFDIA, YSP_SBD%NLEVS, PDTPHY, ZTDS_TDTS, ZTDS_TDTP, ZTDS_TDWS, ZTDS_TDWSI, ZTDS_TDWP,                        &
    & ZTDS_TDWPI, ZTDS_TDWL, ZTDS_TDSNS, ZTDS_TDALBNS, ZTDS_TDRHONS, YDMF_PHYS_SURF%GSD_VP%PTPC, YDMF_PHYS_SURF%GSD_VP%PWPC, &
    & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS_SURF%GSD_VV%PIVEG, YDMF_PHYS_SURF%GSP_RR%PT_T1, YDMF_PHYS_SURF%GSP_SB%PT_T1,     &
    & YDMF_PHYS_SURF%GSP_RR%PW_T1, YDMF_PHYS_SURF%GSP_RR%PIC_T1, YDMF_PHYS_SURF%GSP_SB%PQ_T1, YDMF_PHYS_SURF%GSP_SB%PTL_T1,  &
    & YDMF_PHYS_SURF%GSP_RR%PFC_T1, YDMF_PHYS_SURF%GSP_SG%PF_T1, YDMF_PHYS_SURF%GSP_SG%PA_T1, YDMF_PHYS_SURF%GSP_SG%PR_T1    &
    &                         )  
  ELSE
    IF (LLXFUMSE) THEN
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
    CALL CPNUDG ( YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, NFNUDG, YDCPG_DIM%KFLEVG, YDCPG_DIM%KBL,    &
    & XPNUDG, YDMF_PHYS_SURF%GSD_VF%PNUDM, YDMF_PHYS_SURF%GSP_RR%PT_T1, YDMF_PHYS_SURF%GSP_RR%PW_T1,            &
    & YDMF_PHYS_SURF%GSP_SB%PQ_T1, YDMF_PHYS_SURF%GSP_SG%PF_T1, YLMF_PHYS_NEXT_STATE%T (:, 1:YDCPG_DIM%KFLEVG), &
    & YLMF_PHYS_NEXT_STATE%Q (:, 1:YDCPG_DIM%KFLEVG), YLMF_PHYS_NEXT_STATE%U (:, 1:YDCPG_DIM%KFLEVG),           &
    & YLMF_PHYS_NEXT_STATE%V (:, 1:YDCPG_DIM%KFLEVG), YLMF_PHYS_NEXT_STATE%SP, YDVARS%T%T0, YDVARS%Q%T0,        &
    & YDVARS%U%T0, YDVARS%V%T0, YDCPG_PHY0%PREHYD(:, YDCPG_DIM%KFLEVG), YDVARS%GEOMETRY%GM%T0, YDMF_PHYS_SURF%GSD_VF%PLSM &
    &          )
  ENDIF
ENDIF


CALL MF_PHYS_CVV (YDCPG_DIM, YDVARS%CVV%T0, YDVARS%CVV%T1, YDMODEL)

!        3.3  Store the model trajectory at t-dt (leap-frog) or t (sl2tl).
!             ------------------------------------------------------------

IF (LTRAJPS) THEN
  PTRAJ_PHYS%PQSSMF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_MISC%QS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
  PTRAJ_PHYS%PTSMF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA) =YLMF_PHYS_BASE_STATE%YGSP_RR%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
  PTRAJ_PHYS%PSNSMF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YLMF_PHYS_BASE_STATE%YGSP_SG%F(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)

  IF (.NOT. LTWOTL) THEN
    CALL WRPHTRAJM(YDGEOMETRY, YDSIMPHL, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, PTRAJ_PHYS, YDVARS%U%T9, YDVARS%V%T9, &
    & YDVARS%T%T9, YDVARS%Q%T9, YDVARS%L%T9, YDVARS%I%T9, YDVARS%SP%T9)  
  ENDIF

  IF (LPRTTRAJ.AND.PTRAJ_PHYS%LASTCHUNK) WRITE(NULOUT,*)'GREPTRAJ STORE TRAJ_PHYS in APLPAR'
ENDIF

!     ------------------------------------------------------------------

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
IF (L3DTURB) THEN
  DO JLEV=1,YDCPG_DIM%KFLEVG
    YDCPG_SL2%KAPPAM (YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, JLEV) = ZKUR_KUROV_H(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)
    YDCPG_SL2%KAPPAH (YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, JLEV) = ZKUR_KTROV_H(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)
  ENDDO
ENDIF

CALL MF_PHYS_BAYRAD (YDCPG_DIM, ZBAY_QRCONV, ZBAY_QSCONV, YDVARS%RCONV%T1, YDVARS%SCONV%T1, YDMODEL, &
& LAROME)


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
  CALL WRITEPHYSIO(YDGEOMETRY, YDCPG_MISC, YDCPG_PHY0, YDCPG_PHY9, YDMF_PHYS, YDCPG_DYN0, YDCPG_DYN9,        &
  & YDMF_PHYS_SURF, YDVARS, YDSURF, YDDPHY, YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KFDIA, YDCPG_DIM%KIDIA,    &
  & YDCPG_DIM%KGL1, YDCPG_DIM%KGL2, YDCPG_DIM%KSTGLO, YDCPG_DIM%KSTEP, NTSSG, YSP_SBD%NLEVS, YDVARS%GEOMETRY%GELAM%T0, &
  & YDVARS%GEOMETRY%GEMU%T0, YDVARS%GEOMETRY%GM%T0, YDVARS%GEOMETRY%OROG%T0, YDVARS%GEOMETRY%RCORI%T0, YDVARS%GEOMETRY%RATATH%T0, YDVARS%GEOMETRY%RATATX%T0,         &
  & YDVARS%GEOMETRY%GECLO%T0, YDVARS%GEOMETRY%GESLO%T0, ZRDG_CVGQ, ZRDG_LCVQ, ZRDG_MU0, ZDSA_C1, ZDSA_C2, ZDSA_CPS,              &
  & ZDSA_LHS, ZDSA_RS, ZFLU_CD, ZFLU_CDN, ZFLU_CH, ZFLU_EMIS, ZFLU_FEVI, ZFLU_NEIJ, ZFLU_QSAT,               &
  & ZFLU_QSATS, ZFLU_VEG, IMOC_CLPH, ZMSC_FRMQ, ZMSC_LH, ZMSC_LSCPE, ZMSC_QW, ZMSC_TW, ZPFL_FEFB1,           &
  & ZPFL_FEFB2, ZPFL_FEFB3, ZPFL_FPLCH, ZPFL_FPLSH, ZPFL_FTKE  )
ENDIF

IF (LEDR) THEN
  YDMF_PHYS_SURF%GSD_DI%PXEDR(:,:)=YDMF_PHYS%OUT%EDR(:,:)
ENDIF

CALL MF_PHYS_PRECIPS (YDCPG_DIM, ZPRC_DPRECIPS, ZPRC_DPRECIPS2, YDMF_PHYS_SURF%GSD_XP%PPRECIP, YDMF_PHYS_SURF%GSD_XP2%PPRECIP2, YDMODEL)

!     ------------------------------------------------------------------

!       6. destructor for procset
IF (LINTFLEX) CALL CLEANINTPROCSET(YLPROCSET)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('APLPAR', 1, ZHOOK_HANDLE)
END SUBROUTINE APLPAR
