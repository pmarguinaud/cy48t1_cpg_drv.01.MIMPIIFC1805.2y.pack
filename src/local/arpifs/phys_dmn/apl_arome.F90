#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE APL_AROME(YDGEOMETRY, YDCPG_DIM, YDCPG_MISC, YDCPG_GPAR, YDCPG_PHY0, YDCPG_PHY9, YDMF_PHYS,    &
& YDMF_PHYS_TMP, YDCPG_DYN0, YDCPG_DYN9, YDMF_PHYS_SURF, YDVARS, YDGMV, YDSURF, YDCFU, YDXFU,   &
& YDMODEL, LDCONFX, PDTPHY, PGFL, PKOZO, PGP2DSDT, PB1, PB2, PGMVT1, PGFLT1, PTRAJ_PHYS, &
& YDDDH, PFTCNS)

!**** *APL_AROME * - CALL OF PHYSICAL PARAMETERISATIONS FOR ALARO/AROME

!     Sujet.
!     ------
!     - APPEL DES SOUS-PROGRAMMES DE PARAMETRISATION

!**   Interface.
!     ----------
!        *CALL* *APL_AROME*

!-----------------------------------------------------------------------

! -   ARGUMENTS D'ENTREE.
! -   INPUT ARGUMENTS.
!     -------------------

! - NOM DES PARAMETRES DE DIMENSIONNEMENT DE LA PHYSIQUE.
! - DIMENSIONS.

! KMAXDRAFT : MAX NUMBER OF DRAFTS (FOR DIMENSIONNING)
! KSGST : NUMBER OF SUBGRID SURFACE TEMPERATURES AND FLUXES (NTSSG IN *CPG*)
! KNFRRC : FREQUENCY FOR CLEAR SKY RADIATION CALCULATION
! PDT : TIME STEP (in s) 
! LDXFUMSE : T if CDCONF=X in order not to increment surfex timer in that case
!-----------------------------------------------------------------------
! PGEMU      : SINE OF GEOGRAPHICAL LATITUDE
! PGELAM     :  LONGITUDE
! POROG      : g * OROGRAPHY
! PGM        : MAP FACTOR (used in ALARO convection only)
! PCLON      : cosine of geographical longitude.
! PSLON      : sine of geographical longitude.
! PGP2DSDT   : STOCHASTIC PHYSICS PATTERNS

! FIELDS WITH SUBSCRIPT M FOR TIME T-DT IN 3TL OR T IN 2TL

! PDELPM     : LAYER THICKNESS IN PRESSURE UNITS

! PTM        : TEMPERATURE.
! PCPM        : SPECIFIC HEAT AT CONSTANT PRESSURE FOR AIR
! PRM         : GAS CONSTANT FOR AIR

! PTKEM       : TURBULENT KINETIC ENERGY
! PSVM        : PASSIVE SCALARS
! PUM         : ZONAL WIND
! PVM         : MERIDIAN WIND
! PWM         : VERTICAL VELOCITY (m/s)

!-----------------------------------------------------------------------
! - INOUT

! PGPAR       : BUFFER FOR 2D FIELDS - CONTAINS PRECIP, ALBEDO, EMISS, TS 
!             : SURFACE FLUXES

! ACRANEB2 intermittency storage

!               LINEAR T_e CORRECTION
!               LINEAR T_e CORRECTION
!-----------------------------------------------------------------------
! - OUTPUT (SUBSCRIPT S FOR T+DT)

! PSIGS       : SIGMA FOR SUBGRIDCOND
! PTENDT      : TEMPERATURE TENDENCY
! PTENDR      : HYDROMETEORE TENDENCIES
! PTENDW      : VERTICAL VELOCITY TENDENCY
! PTENDTKE    : TKE TENDENCY
! PTENDEXT    : PASSIVE SCALARS TENDENCY
! PFRSO       : SHORTWAVE RADIATIVE FLUX
! - 2D (0:1)

! variables used in input for radiation in case no surface scheme is used 


! Part of GFL strcture dedicated to easy diagnostics (to be used as a print...)
! PEZDIAG    : MULPITPLE ARRAY TO BE FILLED BY THE USER BY 3D FIELDS
!              (NGFL_EZDIAG ONES)
! output for CFU XFU

!-----------------------------------------------------------------------

!     Externes.
!     ---------

!     Method
!     ------
!     - convert aladin variables into mesonh variables (level inversion 
!       and q to r, t to theta) 
!     - call mesoNH physics and ECMWF radiation scheme
!     - convert mesoNH tendencies to aladin tendencies

!     Auteur.
!     -------
!      S.Malardel et Y. Seity
!      10-03-03
!      big cleaning (18/06/04) S. Malardel and Y. Seity
!     externalisation of surface scheme call + small cleaning (20-07-04) Y.Seity
!     Modifications
!     -------------
!      G. Hello 04-02-06: Add the call of KFB-convection scheme 
!                         for future use in ALARO
!      T.Kovacic 04-05-05: Added ZCVTENDPR_ and ZCVTENDPRS_ 
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!      F.Bouyssel 04-05-05: New arguments in ACRADIN
!     Y. Seity 30-Sept-2005 Add MNH Chemistry scheme
!     R. Zaaboul 15-feb-2006 add surface scheme call
!     T.Kovacic  2006-03-23: calls to subroutines for budgets 
!                             and new arguments PFRTH and PFRSO
!     Y. Seity   2007-05-07: add CFU and XFU calculations 
!                           and call aro_ground_diag
!     S.Ivatek-S 2007-04-17: Over dimensioning of PGPAR by NGPAR+1 just 
!                            (KLON,NGPAR) is used boundary checking bf
!     T.Kovacic  2007-03-16: Fourth dim. in APFT
!     JJMorcrette, ECMWF, 20080325: dummy arguments for RADACT to allow for 
!                        using a new sulphate climatology in the ECMWF model
!     Y. Seity   2008-06-15: correct calculations of PFRTHDS, PFRSODS and PFCLL
!     Y. Seity   2008-09-29: phasing Chemistry corrections
!     O.Riviere  2008-10-01: introduction of new data flow for DDH in Arome
!     Y. Seity   2009-05-03: new version of EDKF and implementation of EDMF
!     Y. Seity   2009-10-03: add missed deallocations 
!     S. Riette  2009-03-25: Arguments modification for AROCLDIA to add HTKERAF
!     K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!     A. Alias   2009-09-01: Sulfate and Volcano aerosols added (call radaer)
!     S. Riette  2010-01-19: ZUM__, ZVM__ and ZDEPTH_HEIGHT_ are given
!                            ARO_GROUND_DIAG in 3D.                     
!     Y. Seity   2010-03-09: add PFEVN and PFEVL
!     Y. Bouteloup 2010-03-26 : Add PQLRAD et PQIRAD
!     Y. Seity : Test TKE > 0.
!     Y. Seity : Optimized version of EDKF + diag HCLS
!     Y. Seity : 2010-09 Save Ts at the end of apl_arome for ICMSH+0000
!     L. Bengtsson (2010): Introduce cloud diagnostics based on geop. 
!                               height (LWMOCLOUD), AND cloud-overlap assumptions 
!                               from C. Wittman 2009 (LACPANMX + WMXOV)
!     S. Riette: 2010-12 aro_ground_diag interface modified
!     Y. Seity: 2010-12 add hail diagnostic
!     R. El Khatib 30-Jun-2010 NEC directive noloopfusion to preserve critical regions
!     P.Marguinaud 2010-06-29 : KSURFEXCTL flag (disable SURFEX) 
!     2010-12    B. Decharme  : modify the radiative coupling with surfex (SW per band in ACRADIN and RADHEAT)
!     2011-02    A. Voldoire : add ZAERINDS to CALL RADAER and ACRADIN
!                              for sulfate indirect effect computation
!     2011-06: M. Jerczynski - some cleaning to meet norms
!     S. Riette: 2011-10 : Modifications for DUAL-MF scheme (according to Wim de Rooy's apl_arome version)
!                          Ice in EDKF
!     Y. Seity : 2012-03 : add LMNHLEV option to revert/or not arrays for MesoNH parameterisations
!     F. Bouttier: 2012-07 add SPPT stochastic physics
!     JJMorcrette, ECMWF, 20120815 additional dummy due to changes in RADACT
!     P. Marguinaud : 2012-09 : Add control threshold for orography
!     Y. Seity : 2013-01 Cleaning LMNHLEV and remove JPVEXT points
!     Y. Seity : 2013-02 Cleaning (add compute_neb)
!     L. Bengtsson: 2013-02: add LOLSMC and LOTOWNC options to compute (or not) cloud sedimentation
!                            using different cloud droplet number conc. depending on land/sea/town.
!     2013-11, D. Degrauwe: Introduction of radflex interface, export
!                           upper-air precipitation fluxes PFPR.
!     2013-11, J. Masek: Inclusion of ACRANEB2 radiation scheme.
!     S. Riette: 2013-11: subgrid precipitation
!     K. Yessad (July 2014): Move some variables.
!     2014-09, C. Wastl: Adaptations for orographic shadowing
!     2014-11, Y. Seity: add TKE budgets for DDH
!     2016-03, E. Bazile: Phasing MUSC for surf_ideal_flux
!     2016-04, J. Masek: LRNUEXP cloud overlap option (COMPUTE_NEB replaced
!                        by ACNPART), passing of sushine duration, fix of
!                        E. Gleeson for ACRANEB2 with SURFEX.
!     2016-09, J. Masek: Proper calculation of sunshine duration in ACRANEB2.
!     2016-10, P. Marguinaud : Port to single precision
!     S. Riette 2016-11: Changes in ICE3/ICE4
!     2018-09, E. Gleeson: Corrected misplaced arguments in ACRANEB2 call.
!     2019-09-24 J.M. Piriou arguments for convective gusts.
!     R. El Khatib 30-Oct-2018 substantial rewrite for optimization and coding standards respect.
!     2018-10, I. Etchevers : add Visibilities
!     2019-01, I. Etchevers, Y. Seity : add Precipitation Type
!     2019-09, J. Masek: Corrected dimensioning of dummy argument PGMU0.
!                        Modified call to ACRANEB2 (clearsky fluxes).
!     2019-10, I. Etchevers : Visibilities in ACVISIH, AROCLDIA=>ACCLDIA
!     2019-10, Y.Bouteloup and M. Bouzghaiam : Radiation modifications. Remove acradin.F90 direct
!              call to recmwf.F90 and add interface to ecrad (in recmwf !)
!     2020-10, J. Masek: Modified call to ACCLDIA.
!     2020-12, F. Meier add call to latent heat nudging if LNUDGLH is TRUE
!     2020-12, U. Andrae : Introduce SPP for HARMONIE-AROME
! End modifications
!-------------------------------------------------------------------------------


USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE MF_PHYS_TYPE_MOD   , ONLY : MF_PHYS_TYPE, MF_PHYS_TMP_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_MISC_TYPE, CPG_GPAR_TYPE
USE CPG_DIM_TYPE_MOD   , ONLY : CPG_DIM_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD,ONLY : MF_PHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMXFU             , ONLY : TXFU
USE YOMCFU             , ONLY : TCFU
USE TYPE_MODEL         , ONLY : MODEL
USE PARKIND1           , ONLY : JPIM     ,JPRB     ,JPRD
USE YOMHOOK            , ONLY : LHOOK    ,DR_HOOK
USE YOESW      , ONLY : RSUN2
USE YOMCST     , ONLY : RG       ,RCPD     ,RD       ,RATM     ,RTT      ,&
          &             RCW      ,RCPV     ,RLVTT    ,RCS      ,RLSTT    ,RGAMW    ,&
          &             RBETW    ,RALPW    ,RGAMS    ,RBETS    ,RALPS    ,RGAMD    ,&
          &             RBETD    ,RALPD    ,RETV     ,RV       ,RKAPPA   ,RHOUR
USE YOMLUN     , ONLY : NULOUT
USE YOMCT0     , ONLY : LTWOTL, LSFORCS
USE YOMVERT    , ONLY : VP00
USE YOMRIP0    , ONLY : NINDAT
USE YOMNUDGLH , ONLY :  LNUDGLH, NSTARTNUDGLH, NSTOPNUDGLH, NINTNUDGLH, NTAUNUDGLH, &
          &             RAMPLIFY,RMAXNUDGLH,RMINNUDGLH,LNUDGLHCOMPT,NTIMESPLITNUDGLH
USE YOMNSV     , ONLY : NSV_CO2
USE DDH_MIX    , ONLY : ADD_FIELD_3D, NEW_ADD_FIELD_3D, TYP_DDH ! for new diag data flow
USE YOMSPSDT   , ONLY : YSPPT_CONFIG, YSPPT
USE SPP_MOD    , ONLY : YSPP_CONFIG, YSPP
USE YOMLSFORC  , ONLY : LMUSCLFA, NMUSCLFA, REMIS_FORC, RALB_FORC
USE INTFLEX_MOD, ONLY : LINTFLEX, LRADFLEX,&
                      & TYPE_INTPROC, TYPE_INTPROCSET,&
                      & NEWINTFIELD, NEWINTPROC
USE YOMMP0     , ONLY : MYPROC     
USE MF_PHYS_BASE_STATE_TYPE_MOD &
               , ONLY : MF_PHYS_BASE_STATE_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_DYN_TYPE, CPG_PHY_TYPE
USE YOMGMV             , ONLY : TGMV
USE SC2PRG_MOD         , ONLY : SC2PRG

USE YOMCT0             , ONLY : LSLAG, LNHDYN, LAROME, LNHQE
USE YOMCVER            , ONLY : LVERTFE  ,LVFE_GWMPA 
USE YOMDYNA            , ONLY : LGWADV, L_RDRY_VD
USE YOMSCM             , ONLY : LGSCM
USE YOMTRAJ            , ONLY : TRAJ_PHYS_TYPE

USE INTFLEX_MOD        , ONLY : NEWINTPROCSET, CLEANINTPROCSET
!     -------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_DIM_TYPE),INTENT(IN)    :: YDCPG_DIM
TYPE(CPG_MISC_TYPE),INTENT(INOUT):: YDCPG_MISC
TYPE(CPG_GPAR_TYPE),INTENT(INOUT):: YDCPG_GPAR
TYPE(CPG_PHY_TYPE),INTENT(INOUT), TARGET :: YDCPG_PHY0
TYPE(CPG_PHY_TYPE),INTENT(INOUT), TARGET :: YDCPG_PHY9
TYPE(MF_PHYS_TYPE),INTENT(INOUT) :: YDMF_PHYS
TYPE(MF_PHYS_TMP_TYPE),INTENT(INOUT) :: YDMF_PHYS_TMP
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
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKOZO(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG*YDMODEL%YRML_PHY_G%YRDPHY%NVCLIS+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP2DSDT(YDCPG_DIM%KLON,YSPPT%YGPSDT(1)%NG2D,YSPPT%N2D)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1(YDCPG_DIM%KLON,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB2(YDCPG_DIM%KLON,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDGMV%YT1%NDIM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLT1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1)
 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFTCNS(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,6)
TYPE (TRAJ_PHYS_TYPE), INTENT(INOUT) :: PTRAJ_PHYS
TYPE(TYP_DDH)     ,INTENT(INOUT) :: YDDDH

!*
!     ------------------------------------------------------------------

! 3D arrays de reference dans mesoNH. En 1D, thetavref=thetavM, mais la question
! concernant la facon d initialiser cette variable dans le 3D reste ouverte (idem pour RHODREF)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                                READ ME, PLEASE !

!                         CODING CONVENTIONS FOR ARPEGE vs MNH PHYSICS

! The horizontal representation in MNH physics is such that 2 dimensions are needed, while a single dimension
! is used in ARPEGE and AROME. The remapping from 1 to 2 dimensions will be made inplicitly through the 
! subroutines interfaces (this is a fortran property). Therefore there is non need to add an explicit dimension
! sized to 1.

! Local 3D arrays with extra levels for Meso-NH turbulence scheme :
! - first dimension is KFDIA not KDLON in order to limit array copies
! - suffixed with two underscore to be easily identified
! These arrays are passed in argument as ZXXX__(:,1:KLEV) except for aro_turb_mnh where they are passed as ZXXX__.

! Local 3D arrays with regular number of levels for Meso-NH interfaces :
! - first dimension is KFDIA not KDLON in order to limit array copies due to non-contiguous data.
! - suffixed with one underscore to be easily identified.

! Local 4D arrays with regular number of levels for Meso-NH interfaces :
! - first dimension is KFDIA not KDLON in order to limit array copies
! - suffixed with one underscore to be easily identified

! Local 2D arrays for Meso-NH interfaces :
! - first dimension is KFDIA not KDLON as a convention as well (no arrays copies to fear)
! - suffixed with one underscore to be easily identified

! Other arrays, which can be dummy arguments, or local but used as argument to IFS/ARPEGE physics
! should remained dimensionned KDLON and should not be suffixed with undersores.

!                         DO NOT USE ARRAY SYNTAX FOR COMPUTATIONAL LOOPS !!

! - They make the code less performant because memory cache is poorly used
! - They can make the code even less readable if the indexes are removed

!                         AVOID ARRAYS COPIES, OR MAKE THEM FAST !

! - if you do need to initialize or copy an array, do it as follows with explicit array syntax in first dimension
! because the compiler will be able to use an optimized function to initialize/copy a segment of memory,
! and may be able to address simultaneously several cach lines :
! 1D array : 
!   Z(KIDIA:KFDIA)=value
! 2D arrays :
!  DO JLEV=1,KLEV
!     ZX(KIDIA:KFDIA,JLEV)=xval
!     ZY(KIDIA:KFDIA,JLEV)=yval
!  ENDDO

! - if you need the bakup of an array, use a swapp mechanism, as what is done here for instance for
! ZSVM_ and ZSVMIN_ : 
!  for the developer ZSVMIN_ is always the bakup and ZSVM_ is always the current array.
!  ZSVMSWAP_ or ZSVMSAVE_ are used for the swapp mechanism, they should never be used by the developer.

! - do not initialize if not necessary. To avoid useless initialization use the mechanizm 'INIT0' coded below :
! IF (INIT0 == 0) THEN
!   ZVALUE=HUGE(1._JPRB)
! ELSE
!   ZVALUE=default_value
! ENDIF
! IF (INIT0 >= 0) THEN
!   Z(:)=ZVALUE
! ENDIF
! INIT0= 0 : initialize to HUGE (testing/debugging)
! INIT0= 1 : initialize to realistic value (discouraged !)
! INIT0=-1 : no initialization (optimized code) - this is the default.

! - pointers can advantageously be used to avoid copies or to avoid the initialization to zero of a sum.
! Example :
!  DO JI=IFIRST,ILAST
!    IF (JI == IFIRST) THEN
!      ! Fill the sum at the first iteration
!      ZARG => ZSUM(:)
!    ELSE
!      ! increment
!      ZARG => ZINC(:)
!    ENDIF
!    CALL COMPUTE(ZARG)
!    IF (JI > IFIRST) THEN
!      ! Add increment
!      ZSUM(KIDIA:KFDIA)=ZSUM(KIDIA:KFDIA)+ZINC(KIDIA:KFDIA)
!    ENDIF
!  ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL(KIND=JPRB) :: ZRHODJM__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1),      ZRHODREFM__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1),   ZPABSM__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1)
REAL(KIND=JPRB) :: ZUM__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1),          ZVM__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1),         ZTHM__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1)
REAL(KIND=JPRB) :: ZUS__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1),          ZVS__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1),         ZWS__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1)
REAL(KIND=JPRB) :: ZTKES_OUT__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1),    ZMF_UP__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1),      ZTHVREFM__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1)  ! thetav de l etat
REAL(KIND=JPRB) :: ZTENDU_TURB__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1),  ZTENDV_TURB__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1), ZTENDTHL_TURB__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1)
REAL(KIND=JPRB) :: ZTENDRT_TURB__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1), ZTKEM__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1),       ZSRCS__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1) 
REAL(KIND=JPRB) :: ZSIGS__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1),        ZEDR__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1)
! THE DDH budgets
REAL(KIND=JPRB) :: ZDP__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1),          ZTP__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1),         ZTPMF__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1)
REAL(KIND=JPRB) :: ZTDIFF__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1),       ZTDISS__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1)
! length scales for momentum and heat for mnh level definitions in case LHARATU=TRUE
REAL(KIND=JPRB) :: ZLENGTHM__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1), ZLENGTHH__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1)

REAL(KIND=JPRB), POINTER :: ZTHS__(:,:)
! horizontal gradients and diagnostics
REAL(KIND=JPRB) :: ZTURB3D__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1,YDMODEL%YRML_PHY_MF%YRARPHY%NGRADIENTS)
! WARNING ! Don't use ZTHSWAP__ or ZTHSAVE__ below because they may be swapped !
! Use only the pointer ZTHS__, and possibly ZTHSIN_ if you need the backup of input data.
REAL(KIND=JPRB), TARGET  :: ZTHSWAP__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1),        ZTHSAVE__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1)
REAL(KIND=JPRB), TARGET  :: ZFLXZTHVMF_SUM__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1), ZWM__(YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG+1)


! Updraft characteristics for Meso-NH world (input of ARO_SHALLOW_MF)
REAL(KIND=JPRB) :: ZTHETAL_UP_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG), ZTHETAV_UP_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG), ZZFRAC_UP_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZRT_UP_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),     ZRC_UP_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),     ZRI_UP_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZZU_UP_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),     ZZV_UP_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),     ZZW_UP_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZZRV_UP_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),    ZTKES_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),      ZZZ_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZDZZ_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),       ZZZ_F_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),      ZDZZ_F_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZCIT_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),       ZMFM_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),       ZEXNREFM_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZSIGM_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),      ZNEBMNH_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),    ZEVAP_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
! additions for MF scheme (Pergaud et al)
REAL(KIND=JPRB) :: ZSIGMF_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),     ZRC_MF_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),     ZRI_MF_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZCF_MF_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),     ZAERD_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),      ZCVTENDT_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZCVTENDRV_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),  ZCVTENDRC_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),  ZCVTENDRI_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZMFS_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),       ZTHLS_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),      ZRTS_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZMFUS_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),      ZMFVS_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG),      ZDEPTH_HEIGHT_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)

REAL(KIND=JPRB), TARGET :: ZFLXZTHVMF_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB), POINTER :: ZARG_FLXZTHVMF_(:,:)


! WARNING ! Don't use ZRSWAP_ or ZRSAVE_ below because they may be swapped !
! Use only the pointer ZRS_, and possibly ZRSIN_ if you need the backup of input data.
REAL(KIND=JPRB), POINTER :: ZRSIN_(:,:,:), ZRS_(:,:,:)
REAL(KIND=JPRB), TARGET :: ZRSWAP_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
REAL(KIND=JPRB), TARGET :: ZRSAVE_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)

REAL(KIND=JPRB), POINTER  :: ZPTRWNU_(:,:), ZTHSIN_(:,:)
REAL(KIND=JPRB), TARGET :: ZWNU_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)

! WARNING ! Don't use ZSVSWAP_ or ZSVSAVE_ below because they may be swapped !
! Use only the pointer ZSVS_, and possibly ZSVSIN_ if you need the backup of input data.
REAL(KIND=JPRB), POINTER :: ZSVSIN_(:,:,:), ZSVS_(:,:,:)
REAL(KIND=JPRB), TARGET :: ZSVSWAP_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB), TARGET :: ZSVSAVE_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)

REAL(KIND=JPRB) :: ZSVXXX_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)

! WARNING ! Don't use ZSVMSWAP_ or ZSVMSAVE_ below because they may be swapped !
! Use only the pointer ZSVM_, and possibly ZSVMIN_ if you need the backup of input data.
REAL(KIND=JPRB), POINTER :: ZSVMIN_(:,:,:), ZSVM_(:,:,:)
REAL(KIND=JPRB), TARGET :: ZSVMSWAP_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) 
REAL(KIND=JPRB), TARGET :: ZSVMSAVE_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB) :: ZSVMB_(YDCPG_DIM%KFDIA,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)

! WARNING ! Don't use ZLIMASWAP_ or ZLIMASAVE_ below because they may be swapped !
! Use only the pointer ZLIMAS_, and possibly ZLIMASIN_ if you need the backup of input data.
REAL(KIND=JPRB), POINTER :: ZLIMAS_(:,:,:), ZLIMASIN_(:,:,:)
REAL(KIND=JPRB), TARGET :: ZLIMASWAP_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NLIMA)
REAL(KIND=JPRB), TARGET :: ZLIMASAVE_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NLIMA)

REAL(KIND=JPRB) :: ZLIMAM_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NLIMA)
!INTEGER(KIND=JPIM) :: KSV_TURB !CPtoclean?
!CPtoclean REAL(KIND=JPRB) :: ZTURBM(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT+YDMODEL%YRML_GCONF%YGFL%NLIMA)
!CPtoclean REAL(KIND=JPRB) :: ZTURBS(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT+YDMODEL%YRML_GCONF%YGFL%NLIMA)
!not (yet ?) used. REK
!REAL(KIND=JPRB) :: ZSFTURB(KLON,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT+YDMODEL%YRML_GCONF%YGFL%NLIMA) ! surf. flux of SV (=0)
!REAL(KIND=JPRB) :: ZTENDSV_TURB2(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT+YDMODEL%YRML_GCONF%YGFL%NLIMA)  ! SV (=0)
REAL(KIND=JPRB) :: ZSFSVLIMA_(YDCPG_DIM%KFDIA,YDMODEL%YRML_GCONF%YGFL%NLIMA) ! surf. flux of LIMA vars
REAL(KIND=JPRB) :: ZTENDSV_TURBLIMA_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NLIMA) ! LIMA

! For radiation scheme
REAL(KIND=JPRB) :: ZRM_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
REAL(KIND=JPRB) :: ZPFPR_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)

REAL(KIND=JPRB) :: ZPEZDIAG_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG)

REAL(KIND=JPRB) :: ZSFSV_(YDCPG_DIM%KFDIA,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT) ! surf. flux of scalars

! Single scattering albedo of dust (points,lev,wvl) :
REAL(KIND=JPRB) :: ZPIZA_DST_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NSWB_MNH)
! Assymetry factor for dust (points,lev,wvl) :
REAL(KIND=JPRB) :: ZCGA_DST_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NSWB_MNH)
! tau/tau_{550} dust (points,lev,wvl) :
REAL(KIND=JPRB) :: ZTAUREL_DST_(YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NSWB_MNH)


! surface flux of theta and surface flux of vapor ; surface flux of CO2
REAL(KIND=JPRB) :: ZSFTH_(YDCPG_DIM%KFDIA),  ZSFRV_(YDCPG_DIM%KFDIA),          ZSFCO2_(YDCPG_DIM%KFDIA)
REAL(KIND=JPRB) :: ZACPRG_(YDCPG_DIM%KFDIA), ZINPRG_NOTINCR_(YDCPG_DIM%KFDIA), ZINPRG_(YDCPG_DIM%KFDIA)
REAL(KIND=JPRB) :: ZACPRR_(YDCPG_DIM%KFDIA), ZINPRR_NOTINCR_(YDCPG_DIM%KFDIA), ZINPRR_(YDCPG_DIM%KFDIA)
REAL(KIND=JPRB) :: ZACPRS_(YDCPG_DIM%KFDIA), ZINPRS_NOTINCR_(YDCPG_DIM%KFDIA), ZINPRS_(YDCPG_DIM%KFDIA)
REAL(KIND=JPRB) :: ZCFBTH_(YDCPG_DIM%KFDIA), ZINPRH_NOTINCR_(YDCPG_DIM%KFDIA), ZINPRH_(YDCPG_DIM%KFDIA)
REAL(KIND=JPRB) :: ZZS_(YDCPG_DIM%KFDIA),    ZSSO_STDEV_(YDCPG_DIM%KFDIA),     ZALB_UV_(YDCPG_DIM%KIDIA)
REAL(KIND=JPRB) :: ZLAT_(YDCPG_DIM%KIDIA),   ZLON_(YDCPG_DIM%KIDIA),           ZZENITH_(YDCPG_DIM%KIDIA)
REAL(KIND=JPRB) :: ZGZ0_(YDCPG_DIM%KFDIA),   ZGZ0H_(YDCPG_DIM%KFDIA),          ZTOWNS_(YDCPG_DIM%KFDIA) 
REAL(KIND=JPRB) :: ZCFAQ_(YDCPG_DIM%KFDIA),  ZCFBQ_(YDCPG_DIM%KFDIA),          ZCFATH_(YDCPG_DIM%KFDIA)
REAL(KIND=JPRB) :: ZCFAU_(YDCPG_DIM%KFDIA),  ZCFBU_(YDCPG_DIM%KFDIA),          ZCFBV_(YDCPG_DIM%KFDIA)
REAL(KIND=JPRB) :: ZBUDTH_ (YDCPG_DIM%KFDIA),ZBUDSO_(YDCPG_DIM%KFDIA),         ZFCLL_(YDCPG_DIM%KFDIA)
REAL(KIND=JPRB) :: ZCD_(YDCPG_DIM%KFDIA),    ZSEA_(YDCPG_DIM%KFDIA),           ZTOWN_(YDCPG_DIM%KFDIA)
REAL(KIND=JPRB) :: ZZTOP_(YDCPG_DIM%KFDIA),  ZCVTENDPR_(YDCPG_DIM%KFDIA),      ZCVTENDPRS_(YDCPG_DIM%KFDIA)
! surface flux of x and y component of wind. are they really necessary ? REK
REAL(KIND=JPRB) :: ZSFU_(YDCPG_DIM%KFDIA),   ZSFV_(YDCPG_DIM%KFDIA)


! Arpege-style dimensionning :
! --------------------------

!Variables used in case LHARATU=TRUE
! length scales for momentum and heat and TKE
REAL(KIND=JPRB) :: ZLENGTH_M(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZLENGTH_H(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTKEEDMF(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTKEEDMFS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

REAL(KIND=JPRB) :: ZEMIS (YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZQICE(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZQLIQ(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

REAL(KIND=JPRB) :: ZAER(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,6)
REAL(KIND=JPRB) :: ZAERINDS(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

REAL(KIND=JPRB) :: ZQSAT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

REAL(KIND=JPRB) :: ZFRSOFS(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZLH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZLSCPE(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZGEOSLC(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZQDM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZQV(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZQCO2(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB), TARGET :: ZTKEM4SLDDH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB), POINTER :: ZTKEM(:,:)
REAL(KIND=JPRB) :: ZQW(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZTW(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

REAL(KIND=JPRB) :: ZTENT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTENDTT(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG) ! array to save heating profile for LHN
REAL(KIND=JPRB) :: ZMAXTEND,ZMINTEND
REAL(KIND=JPRB) :: ZDZZ(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZTPW(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

! POUR GROUND 
REAL(KIND=JPRB) :: ZZS_FSWDIR(YDCPG_DIM%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW), ZZS_FSWDIF(YDCPG_DIM%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB) :: ZTRSODIR (YDCPG_DIM%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW), ZTRSODIF (YDCPG_DIM%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)

REAL(KIND=JPRB) :: ZTRSOD(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZSUDU(YDCPG_DIM%KLON), ZSDUR(YDCPG_DIM%KLON), ZDSRP(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZCEMTR(YDCPG_DIM%KLON,2), ZCTRSO(YDCPG_DIM%KLON,2)

REAL(KIND=JPRB) :: ZALBD(YDCPG_DIM%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW), ZALBP(YDCPG_DIM%KLON,YDMODEL%YRML_PHY_RAD%YRERAD%NSW)
REAL(KIND=JPRB) :: ZALBD1(YDCPG_DIM%KLON), ZALBP1(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZAPHIM(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG), ZAPHIFM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

REAL(KIND=JPRB) :: ZTM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZQVM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZQIM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZQCM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG),ZQHM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZQHGM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZQRM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZQSM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZQGM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZUPGENL(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZUPGENN(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZCLFR(YDCPG_DIM%KLON)

REAL(KIND=JPRB) :: ZCPM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZRHM(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

! Variables concerning updraft rain/snow for EDMF
REAL(KIND=JPRB) :: ZTENDTUP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZTENDQVUP(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

! specific to new data flow for diagnostics
REAL(KIND=JPRB) :: ZTENDTBAK(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG), ZTENDRBAK(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
REAL(KIND=JPRB) :: ZTMPAF(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

! daand: radflex
REAL(KIND=JPRB)  :: ZFPR(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,YDMODEL%YRML_PHY_MF%YRPARAR%NRR)

! Target should not be necessary. REK
REAL(KIND=JPRB), TARGET :: ZCON1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB), TARGET :: ZCON2(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB), TARGET :: ZCON3(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

! Ajout pour MF Dual Scheme (KNMI et al)
! Updraft characteristics in Arpege/IFS world
REAL(KIND=JPRB) :: ZMF_UP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,YDCPG_DIM%KMAXDRAFT)
REAL(KIND=JPRB) :: ZTHETAL_UP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,YDCPG_DIM%KMAXDRAFT)
REAL(KIND=JPRB) :: ZQT_UP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,YDCPG_DIM%KMAXDRAFT)
REAL(KIND=JPRB) :: ZTHTV_UP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,YDCPG_DIM%KMAXDRAFT)
REAL(KIND=JPRB) :: ZQC_UP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,YDCPG_DIM%KMAXDRAFT)
REAL(KIND=JPRB) :: ZQI_UP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,YDCPG_DIM%KMAXDRAFT)
REAL(KIND=JPRB) :: ZU_UP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,YDCPG_DIM%KMAXDRAFT)
REAL(KIND=JPRB) :: ZV_UP(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,YDCPG_DIM%KMAXDRAFT)
REAL(KIND=JPRB) :: ZTSURF(YDCPG_DIM%KLON), ZTN(YDCPG_DIM%KLON), ZQS(YDCPG_DIM%KLON)

REAL(KIND=JPRB) :: ZFRSOLU(YDCPG_DIM%KLON), ZFRSODS(YDCPG_DIM%KLON)
REAL(KIND=JPRB) :: ZFSDNN(YDCPG_DIM%KLON), ZFSDNV(YDCPG_DIM%KLON)

REAL(KIND=JPRB) :: ZSURFPREP(YDCPG_DIM%KLON), ZSURFSNOW(YDCPG_DIM%KLON)

REAL(KIND=JPRB) :: ZQO3(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)

REAL(KIND=JPRB) :: ZZS_FTH_(YDCPG_DIM%KLON), ZZS_FRV_(YDCPG_DIM%KLON), ZZS_FU_(YDCPG_DIM%KLON), ZZS_FV_(YDCPG_DIM%KLON)

! Surface forcing arrays for MUSC
REAL(KIND=JPRB) :: ZRHODREFM(YDCPG_DIM%KLON), ZTHETAS(YDCPG_DIM%KLON)

! ACRANEB2 local variables
REAL(KIND=JPRB) :: ZNEB0    (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)  ! protected cloud fractions
REAL(KIND=JPRB) :: ZCLCT_RAD(YDCPG_DIM%KLON)       ! total cloud cover for radiation
REAL(KIND=JPRB) :: ZDECRD   (YDCPG_DIM%KLON)       ! decorrelation depth

! Stochastic physics pattern & dummy tendencies for calling sppten
! Bof. REK
REAL(KIND=JPRB) :: ZDUMMY(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
REAL(KIND=JPRB) :: ZDUMMY1(YDCPG_DIM%KLON)

REAL(KIND=JPRB) :: ZROZ(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

! Can we remove ? REK
REAL(KIND=JPRB) :: ZEPSM(0,0,0) ! Dissipation of TKE (eps) at time t-dt
REAL(KIND=JPRB) :: ZEPSS(0,0,0) ! Dissipation of TKE at time t+dt


!    Integers
INTEGER(KIND=JPIM) :: JLEV, JLON, JRR, JGFL, JGR, ISPLITR
INTEGER(KIND=JPIM) :: IJN  ! max. number of day/night slices within NRPOMA
INTEGER(KIND=JPIM) :: IKL  !ordering of vert levels 1:MNH -1:AROME
INTEGER(KIND=JPIM) :: IOFF_MFSHAL, IEZDIAG_CHEM
INTEGER(KIND=JPIM) :: IKA,IKB,IKU,IKT,IKTE,IKTB ! vertical points as in mpa
INTEGER(KIND=JPIM) :: JSG, JK, JR, JSW
INTEGER(KIND=JPIM) :: IDRAFT,JDRAFT,INDRAFT
INTEGER(KIND=JPIM) :: ISURFEX
INTEGER(KIND=JPIM) :: IDAY,IYEAR,IMONTH

INTEGER(KIND=JPIM) :: INIT0 ! Kind of safety/debugging initialization :
                            ! 0 = initialize to HUGE (debugging)
                            ! 1 = initialize to realistic value (discouraged)
                            ! -1 = no initialization (optimized code) - this is the default.

INTEGER(KIND=JPIM) :: ICLPH(YDCPG_DIM%KLON)             !PBL top level
INTEGER(KIND=JPIM) :: JLHSTEP,ISTEP

!       Real
REAL(KIND=JPRB) :: ZRHO
REAL(KIND=JPRB) :: ZAEO, ZAEN, ZSALBCOR
REAL(KIND=JPRB) :: ZDT, ZDT2, ZINVDT, ZINVG, ZRSCP, ZINVATM, Z_WMAX, Z_WMIN
 ! pas de temps pour la surface externalise
REAL(KIND=JPRB) :: ZDTMSE,ZRHGMT,ZSTATI
REAL(KIND=JPRB) :: ZDELTA
REAL(KIND=JPRB) :: ZEPSNEB

! default values for initialization :
REAL(KIND=JPRB) :: ZVALUE, ZVALUE_ONE, ZVALUE_T, ZVALUE_P, ZVALUE_L, ZVALUE_EPSILON

REAL(KIND=JPRB) :: ZVETAH(0:YDCPG_DIM%KFLEVG)

!       Boolean
LOGICAL :: LLMSE, LLMSE_PARAM, LLMSE_DIAG
LOGICAL :: LLAROME
LOGICAL :: LLRAD
LOGICAL :: LLSWAP_THS, LLSWAP_RS, LLSWAP_SVS, LLSWAP_SVM, LLSWAP_LIMAS ! logical to swap or not pointers in and out
LOGICAL :: LLHN(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)
LOGICAL :: LNUDGLHNREAD

!       Characters
CHARACTER(LEN=11) :: CLNAME
CHARACTER(LEN=2),DIMENSION(7):: CLVARNAME=(/"QV","QL","QR","QI","QS","QG","QH"/)

! daand: radflex
REAL(KIND=JPRB), POINTER :: ZFRSO(:,:), ZFRTH(:,:)
TYPE(TYPE_INTPROC), POINTER :: YLRADPROC
REAL(KIND=JPRB)   :: ZCAPE(YDCPG_DIM%KLON), ZDCAPE(YDCPG_DIM%KLON)

!
! Phaser team note from CY43T1:
! there was a USE MODD_CTURB for accessing XTKEMIN here, but that created a forbidden
! dependence of APL_AROME (in "ifsarp") to the Méso-NH/Arome interfaces (in "mpa").
! There should be no USE MODD_* in APL_*.
! We decided to change the variable here to a local one, with the classical initial value for TKE.
!
REAL(KIND=JPRB), PARAMETER :: PPTKEMIN = 1.E-6


! Perturbed radiation-cloud interaction coef
REAL(KIND=JPRB), DIMENSION (YDCPG_DIM%KLON) :: ZRADGR,ZRADSN
REAL(KIND=JPRB) :: ZMU,ZVAL
INTEGER(KIND=JPIM) :: JKO,JKE

!     ------------------------------------------------------------------
LOGICAL :: LLDIAB
LOGICAL :: LL_SAVE_PHSURF
LOGICAL :: LLXFUMSE

INTEGER(KIND=JPIM) :: IPTREXT,IEFB1,IEFB2,IEFB3
INTEGER(KIND=JPIM) :: IPTR(YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
INTEGER(KIND=JPIM) :: IPTRLIMA
INTEGER(KIND=JPIM) :: IRR ! pointer of 1st hydrometeors in ZTENDGFLR
INTEGER(KIND=JPIM) :: IPTRTKE ! pointer of TKE in ZTENDGFLR

INTEGER(KIND=JPIM) :: IPGFL(YDMODEL%YRML_GCONF%YGFL%NUMFLDS)

INTEGER(KIND=JPIM) :: INSTEP_DEB,INSTEP_FIN
INTEGER(KIND=JPIM) :: JROF
INTEGER(KIND=JPIM) :: ISLB1U9  ,ISLB1V9  ,ISLB1T9  ,ISLB1GFL9, ISLB1VD9

!     --- UPPER AIR PHYSICAL TENDENCIES.
REAL(KIND=JPRB) :: ZTENDH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)     ! Enthalpy tendency.
REAL(KIND=JPRB) :: ZTENDPTKE(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)  ! Pseudo progn. TKE

! GFL tendencies for APL_AROME (assumes YDMODEL%YRML_GCONF%YGFL%NUMFLDS>=YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
! for now, use Jovi's trick :
REAL(KIND=JPRB), TARGET :: ZTENDGFLR(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,0:YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
REAL(KIND=JPRB) :: ZTENDGFL(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NUMFLDS)   ! GFL tendencies

!     --- UPPER AIR PHYSICAL TENDENCIES FOR AROME.
!       (the previous one are not used in AROME)
REAL(KIND=JPRB) :: ZTENDT (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)        ! temperature tendency
REAL(KIND=JPRB) :: ZTENDW (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)        ! W  tendency
REAL(KIND=JPRB) :: ZTENDD (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)        ! d  tendency

REAL(KIND=JPRB) :: ZTENDU (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! U tendency without deep convection contribution
REAL(KIND=JPRB) :: ZTENDV (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! V tendency without deep convection contribution


!     ---FOR AROME PHYSICS  ---
REAL(KIND=JPRB) :: ZGWT1(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)     ! vertical velocity calculated by cputqy_arome before convertion in d
REAL(KIND=JPRB) :: ZTT1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)       ! Temperature at t1

! ZRTT1: appropriate version of R*T at t1 for gnhgw2svd
!  Version of R must be consistent with definition of vertical divergence.
REAL(KIND=JPRB) :: ZRTT1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

! required for INTFLEX
TYPE(TYPE_INTPROCSET) :: YLPROCSET

! SPP
REAL(KIND=JPRB) :: ZGP2DSPP(YDCPG_DIM%KLON,YSPP%N2D)

REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDR    (:,:,:) 
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDLIMA (:,:,:) 
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDTKE  (:,:) 
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDEFB1 (:,:) 
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDEFB2 (:,:)
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDEFB3 (:,:)
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDEXT  (:,:,:)

REAL(KIND=JPRB), POINTER :: ZP1EZDIAG(:,:,:)

TYPE (MF_PHYS_BASE_STATE_TYPE) :: YLMF_PHYS_STATE

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "cpphinp.intfb.h"
#include "cptend_flex.intfb.h"
#include "cputqy_arome.intfb.h"
#include "cputqy.intfb.h"
#include "cp_ptrslb1.intfb.h"
#include "gnhgw2svdarome.intfb.h"
#include "mf_phys_init.intfb.h"
#include "writephysio.intfb.h"
#include "mf_phys_nhqe_part1.intfb.h"
#include "mf_phys_nhqe_part2.intfb.h"
#include "mf_phys_save_phsurf_part1.intfb.h"
#include "mf_phys_save_phsurf_part2.intfb.h"
#include "mf_phys_transfer.intfb.h"
#include "mf_phys_precips.intfb.h"
#include "apl_arome_calc_iptr.intfb.h"
#include "apl_arome_calc_ipgfl.intfb.h"
#include "abor1.intfb.h"
#include "recmwf.intfb.h"
#include "acraneb2.intfb.h"
#include "actqsat.intfb.h"
#include "acnpart.intfb.h"
#include "bri2acconv.intfb.h"
#include "gpgeo.intfb.h"
#include "gprcp.intfb.h"
#include "radheat.intfb.h"
#include "suozon.intfb.h"
#include "radaer.intfb.h"
#include "radozc.intfb.h"
#include "accldia.intfb.h"
#include "vdfhghthl.intfb.h"
#include "sppten.intfb.h"
#include "surf_ideal_flux.intfb.h"
#include "ecr1d.intfb.h"
#include "apl_arome2intflex.intfb.h"
#include "aro_rain_ice.h"
#include "nudglhprecip.intfb.h"
#include "nudglh.intfb.h"
#include "nudglhclimprof.intfb.h"
#include "nudglhprep.intfb.h"
#include "aro_turb_mnh.h"
#include "aro_adjust.h"
#include "aro_mnhc.h"
#include "aro_mnhdust.h"
#include "aro_startbu.h"
#include "aro_convbu.h"
#include "aro_ground_param.h"
#include "aro_ground_diag.h"
#include "aro_shallow_mf.h"
#include "aro_rainaero.h"
#include "aro_lima.h"
#include "diagflash.intfb.h"
#include "dprecips.intfb.h"
#include "ppwetpoint.intfb.h"
#include "acvisih.intfb.h"
#include "aro_ground_diag_2isba.h"

!     ------------------------------------------------------------------

#include "fcttrm.func.h"


!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('APL_AROME', 0, ZHOOK_HANDLE)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDGEM=>YDGEOMETRY%YRGEM, YDCSGEOM=> YDGEOMETRY%YRCSGEOM(YDCPG_DIM%KBL),                         &
& YDGSGEOM=>YDGEOMETRY%YRGSGEOM(YDCPG_DIM%KBL), YDOROG=>YDGEOMETRY%YROROG(YDCPG_DIM%KBL), YDSTA=>YDGEOMETRY%YRSTA,                 &
& YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, YLDDH=>YDMODEL%YRML_DIAG%YRLDDH, YDTOPH=>YDMODEL%YRML_PHY_MF%YRTOPH,                           &
& YDERDI=>YDMODEL%YRML_PHY_RAD%YRERDI, YDRIP=>YDMODEL%YRML_GCONF%YRRIP, YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY,                      &
& YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD, YDPHY3=>YDMODEL%YRML_PHY_MF%YRPHY3, YDPHYDS=>YDMODEL%YRML_PHY_MF%YRPHYDS,                   &
& YDPHY0=>YDMODEL%YRML_PHY_MF%YRPHY0, YDVISI=>YDMODEL%YRML_PHY_MF%YRPHY%YRDVISI, YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2,               &
& YGFL=>YDMODEL%YRML_GCONF%YGFL, YDMSE=>YDMODEL%YRML_PHY_MF%YRMSE, YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR,                           &
& YDPRECIPS=>YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS, YDSTOPH=>YDMODEL%YRML_PHY_STOCH%YRSTOPH, YDPTRSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1, &
& YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH, YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY, YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH,                           &
& YDEPHY=> YDMODEL%YRML_PHY_EC%YREPHY)

ASSOCIATE(MINPRR=>YDPARAR%MINPRR, MINPRS=>YDPARAR%MINPRS, MVQS=>YDPARAR%MVQS, MINPRG=>YDPARAR%MINPRG,                                 &
& LOTOWNC=>YDPARAR%LOTOWNC, LFPREC3D=>YDPARAR%LFPREC3D, NRRI=>YDPARAR%NRRI, NRRL=>YDPARAR%NRRL, CSUBG_AUCV_RC=>YDPARAR%CSUBG_AUCV_RC, &
& LTOTPREC=>YDPARAR%LTOTPREC, NPRINTFR=>YDPARAR%NPRINTFR, CMF_CLOUD=>YDPARAR%CMF_CLOUD, MALBDIR=>YDPARAR%MALBDIR,                     &
& NSWB_MNH=>YDPARAR%NSWB_MNH, XSW_BANDS=>YDPARAR%XSW_BANDS, MACPRG=>YDPARAR%MACPRG, MSWDIR=>YDPARAR%MSWDIR,                           &
& LMIXUV=>YDPARAR%LMIXUV, MSWDIF=>YDPARAR%MSWDIF, LOLSMC=>YDPARAR%LOLSMC, NDIAGWMAX=>YDPARAR%NDIAGWMAX,                               &
& MACPRS=>YDPARAR%MACPRS, MACPRR=>YDPARAR%MACPRR, LSQUALL=>YDPARAR%LSQUALL, VSIGQSAT=>YDPARAR%VSIGQSAT,                               &
& MALBSCA=>YDPARAR%MALBSCA, RADSN=>YDPARAR%RADSN, LOSEDIC=>YDPARAR%LOSEDIC, LDIAGWMAX=>YDPARAR%LDIAGWMAX,                             &
& CSEDIM=>YDPARAR%CSEDIM, NPTP=>YDPARAR%NPTP, NSPLITR=>YDPARAR%NSPLITR, NSPLITG=>YDPARAR%NSPLITG, NSV=>YDPARAR%NSV,                   &
& CFRAC_ICE_SHALLOW_MF=>YDPARAR%CFRAC_ICE_SHALLOW_MF, CFRAC_ICE_ADJUST=>YDPARAR%CFRAC_ICE_ADJUST, MVTS=>YDPARAR%MVTS,                 &
& NREFROI2=>YDPARAR%NREFROI2, NREFROI1=>YDPARAR%NREFROI1, MVEMIS=>YDPARAR%MVEMIS, LOWARM=>YDPARAR%LOWARM,                             &
& LOCND2=>YDPARAR%LOCND2, LGRSN=>YDPARAR%LGRSN, LOSIGMAS=>YDPARAR%LOSIGMAS, NRR=>YDPARAR%NRR, LOSUBG_COND=>YDPARAR%LOSUBG_COND,       &
& RADGR=>YDPARAR%RADGR, CMF_UPDRAFT=>YDPARAR%CMF_UPDRAFT, LHARATU=>YDPARAR%LHARATU, XMINLM=>YDPHY0%XMINLM,                            &
& XMAXLM=>YDPHY0%XMAXLM, AERCS1=>YDPHY0%AERCS1, AERCS3=>YDPHY0%AERCS3, AERCS5=>YDPHY0%AERCS5, RDECRD1=>YDPHY0%RDECRD1,                &
& RDECRD2=>YDPHY0%RDECRD2, RDECRD3=>YDPHY0%RDECRD3, RDECRD4=>YDPHY0%RDECRD4, LMPA=>YDARPHY%LMPA, LUSECHEM=>YDARPHY%LUSECHEM,          &
& LKFBCONV=>YDARPHY%LKFBCONV, LMFSHAL=>YDARPHY%LMFSHAL, LMICRO=>YDARPHY%LMICRO, CCOUPLING=>YDARPHY%CCOUPLING,                         &
& LTURB=>YDARPHY%LTURB, LGRADHPHY=>YDARPHY%LGRADHPHY, NSURFEX_ITER=>YDARPHY%NSURFEX_ITER, LRDUST=>YDARPHY%LRDUST,                     &
& NGRADIENTS=>YDARPHY%NGRADIENTS, LRDEPOS=>YDARPHY%LRDEPOS, LSURFEX_CRITICAL=>YDARPHY%LSURFEX_CRITICAL,                               &
& LRCO2=>YDARPHY%LRCO2, LMSE=>YDARPHY%LMSE, LSURFEX_KFROM=>YDARPHY%LSURFEX_KFROM, NSURFEXCTL=>YDMSE%NSURFEXCTL,                       &
& XZSEPS=>YDMSE%XZSEPS, NDLUNG=>YDDIM%NDLUNG, NDGUNG=>YDDIM%NDGUNG, NPROMA=>YDDIM%NPROMA, NDLUXG=>YDDIM%NDLUXG,                       &
& NDGUXG=>YDDIM%NDGUXG, NGFL_EXT=>YGFL%NGFL_EXT, YLRAD=>YGFL%YLRAD, YIRAD=>YGFL%YIRAD, NGFL_EZDIAG=>YGFL%NGFL_EZDIAG,                 &
& NLIMA=>YGFL%NLIMA, CMICRO=>YDPARAR%CMICRO, YSD_VAD=>YDSURF%YSD_VAD, QCO2=>YDPHY3%QCO2, NRAY=>YDPHY%NRAY,                            &
& LRAYFM=>YDPHY%LRAYFM, LO3ABC=>YDPHY%LO3ABC, LRAY=>YDPHY%LRAY, LRSTAER=>YDPHY%LRSTAER, LRNUEXP=>YDPHY%LRNUEXP,                       &
& AMAGSTOPH_CASBS=> YDSTOPH%AMAGSTOPH_CASBS, LFORCENL=>YDSTOPH%LFORCENL, NFORCESTART=>YDSTOPH%NFORCESTART,                            &
& NFORCEEND=>YDSTOPH%NFORCEEND, NTRADI=>YDTOPH%NTRADI, NTQSAT=>YDTOPH%NTQSAT, NTNEBU=>YDTOPH%NTNEBU,                                  &
& NAER=>YDERAD%NAER, LHLRADUPD=>YDPHY%LHLRADUPD, TSPHY=>YDPHY2%TSPHY, NOZOCL=>YDERAD%NOZOCL, NRADFR=>YDERAD%NRADFR,                   &
& NSW=>YDERAD%NSW, RCARDI=>YDERDI%RCARDI, LFLEXDIA=>YLDDH%LFLEXDIA, LDDH_OMP=>YLDDH%LDDH_OMP, LRSLDDH=>YLDDH%LRSLDDH,                 &
& RDECLI=>YDRIP%RDECLI, RCODEC=>YDRIP%RCODEC, RHGMT=>YDRIP%RHGMT, RSIDEC=>YDRIP%RSIDEC, RSOVR=>YDRIP%RSOVR,                           &
& RSTATI=>YDRIP%RSTATI, TSTEP=>YDRIP%TSTEP, STPREH=>YDSTA%STPREH, LXXDIAGH=>YDXFU%LXXDIAGH, LFLASH =>YDCFU%LFLASH,                    &
& LDPRECIPS=>YDPHY%LDPRECIPS, LDPRECIPS2=>YDPHY%LDPRECIPS2, NDTPREC=>YDPRECIPS%NDTPREC, NDTPREC2=>YDPRECIPS%NDTPREC2,                 &
& NDTPRECCUR=>YDPRECIPS%NDTPRECCUR, NDTPRECCUR2=>YDPRECIPS%NDTPRECCUR2, NGPTOT=>YDGEM%NGPTOT, NGPBLKS=>YDDIM%NGPBLKS,                 &
& NTSSG=>YDDPHY%NTSSG, YEZDIAG=>YGFL%YEZDIAG, YEXT=>YGFL%YEXT, YNOGW=>YGFL%YNOGW, YCHEM=>YGFL%YCHEM,                                  &
& YSP_SBD=>YDSURF%YSP_SBD, LEDR=>YDPHY%LEDR, LAGPHY=>YDEPHY%LAGPHY, YLIMA=>YGFL%YLIMA)

CALL SC2PRG(1, YEZDIAG(:)%MP, YDMODEL%YRML_GCONF%YGFL%NGFL_EZDIAG, PGFL, ZP1EZDIAG)

CALL YLMF_PHYS_STATE%INIT (LTWOTL, YDCPG_DYN0, YDCPG_DYN9, YDCPG_PHY0, YDCPG_PHY9, YDVARS, &
& YDMF_PHYS_SURF, PGFL=PGFL, YDMODEL=YDMODEL)

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
LLDIAB=(.NOT.LAGPHY)

! In the NHQE model, APL_AROME enters with Tt and grad(Tt), where Tt = T * exp(-(R/cp) log(pre/prehyd)).
! But calculations of APL_AROME must use T and grad(T).
! So we do a conversion Tt -> T.
IF (LNHQE) THEN
  CALL MF_PHYS_NHQE_PART1 (YDGEOMETRY, YDCPG_DIM, YDMF_PHYS_TMP%NHQ%TT0, YDMF_PHYS_TMP%NHQ%TT0L, YDMF_PHYS_TMP%NHQ%TT0M, YDMF_PHYS_TMP%NHQ%TT9, YDVARS, YDMODEL, PGFL)
ENDIF

IF (LLDIAB) THEN
  CALL CPPHINP(YDGEOMETRY, YDMODEL, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDGSGEOM%GEMU, YDGSGEOM%GELAM,           &
  & YDVARS%U%T0, YDVARS%V%   T0, YDVARS%Q%T0, YDVARS%Q%DL, YDVARS%Q%DM, YDVARS%CVGQ%DL, YDVARS%CVGQ%DM,        &
  & YDCPG_PHY0%XYB%RDELP, YDCPG_DYN0%CTY%EVEL, YDVARS%CVGQ%T0, YDMF_PHYS_TMP%RDG%MU0, YDMF_PHYS_TMP%RDG%MU0LU, &
  & YDMF_PHYS_TMP%RDG%MU0M, YDMF_PHYS_TMP%RDG%MU0N, YDMF_PHYS_TMP%RDG%CVGQ)
  YDMF_PHYS_TMP%RDG%LCVQ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YDMF_PHYS_TMP%RDG%CVGQ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
ENDIF

! * In some cases, some pseudo-historic surface buffers (like z0) should
!   not be modified between the entrance and the output of APL_AROME
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
  LL_SAVE_PHSURF=LDCONFX
  IF (LL_SAVE_PHSURF) THEN
    CALL MF_PHYS_SAVE_PHSURF_PART1 (YDCPG_DIM, YDMF_PHYS_TMP%SAV%DDAL, YDMF_PHYS_TMP%SAV%DDOM, YDMF_PHYS_TMP%SAV%ENTCH, YDMF_PHYS_TMP%SAV%FHPS, YDMF_PHYS_TMP%SAV%GZ0F, &
                                  & YDMF_PHYS_TMP%SAV%GZ0HF, YDMF_PHYS_TMP%SAV%HV, YDMF_PHYS_TMP%SAV%PBLH, YDMF_PHYS_TMP%SAV%QSH, YDMF_PHYS_TMP%SAV%UDAL, &
                                  & YDMF_PHYS_TMP%SAV%UDGRO, YDMF_PHYS_TMP%SAV%UDOM, YDMF_PHYS_TMP%SAV%UNEBH, YDMF_PHYS_SURF, YDVARS, YDSURF, YDMODEL)
  ENDIF
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

CALL APL_AROME_CALC_IPGFL (YDGEOMETRY, YDCPG_DIM, YDMODEL, IPGFL)

CALL MF_PHYS_TRANSFER (YDCPG_DIM, YDVARS, YDMODEL)

CALL APL_AROME_CALC_IPTR (YDMODEL, IEFB1, IEFB2, IEFB3, IPTR, IPTREXT, IPTRLIMA, IPTRTKE, IRR)

! If an incorrect address is used, then the initialization below will detect it :
ZTENDGFLR(:,:,0)=HUGE(1._JPRB)

ZTENDR    => ZTENDGFLR (:, :, IRR:IRR+YDMODEL%YRML_PHY_MF%YRPARAR%NRR-1)
ZTENDLIMA => ZTENDGFLR (:, :, IPTRLIMA:IPTRLIMA+YDMODEL%YRML_GCONF%YGFL%NLIMA-1)
ZTENDTKE  => ZTENDGFLR (:, :, IPTRTKE)
ZTENDEFB1 => ZTENDGFLR (:, :, IEFB1)
ZTENDEFB2 => ZTENDGFLR (:, :, IEFB2)
ZTENDEFB3 => ZTENDGFLR (:, :, IEFB3)
ZTENDEXT  => ZTENDGFLR (:, :, IPTREXT:IPTREXT+YDMODEL%YRML_GCONF%YGFL%NGFL_EXT-1)


!     --------------------------------------------------------------------------

!    ------------------------------------------------------------------
!     1 - Initialisations
!    - --------------------------------------------------------------------

INIT0=-1

IF (INIT0 == 0) THEN
  ZVALUE=HUGE(1._JPRB)
  ZVALUE_ONE=HUGE(1._JPRB)
  ZVALUE_T=HUGE(1._JPRB)
  ZVALUE_P=HUGE(1._JPRB)
  ZVALUE_L=HUGE(1._JPRB)
  ZVALUE_EPSILON=HUGE(1._JPRB)
ELSE
  ZVALUE=0._JPRB
  ZVALUE_ONE=1._JPRB
  ZVALUE_T=293._JPRB
  ZVALUE_P=101325._JPRB
  ZVALUE_L=0.01_JPRB
  ZVALUE_EPSILON=1E-12_JPRB
ENDIF

LLSWAP_THS=.TRUE. ! it can be as well true as false, actually : this is just to start up the swapp process
LLSWAP_RS=.TRUE. ! it can be as well true as false, actually : this is just to start up the swapp process
LLSWAP_SVS=.TRUE. ! it can be as well true as false, actually : this is just to start up the swapp process
LLSWAP_SVM=.TRUE. ! it can be as well true as false, actually : this is just to start up the swapp process
LLSWAP_LIMAS=.TRUE. ! it can be as well true as false, actually : this is just to start up the swapp process

!        1.0 numerical safety

IF (JPRD == JPRB) THEN
  ZEPSNEB=1.E-12
ELSE
  ZEPSNEB=1.E-06
ENDIF

NSV=0

!         1.3 time step initialisation
!             the mesoNH physics (turb and microphysics) is written 
!             for leap frog scheme
!             !!! be carefull for 2TL or 3TL 

IF (LTWOTL) THEN
  ZDT=PDTPHY/2._JPRB
ELSE
  IF (YDCPG_DIM%KSTEP/=0) THEN
    ZDT=PDTPHY/2._JPRB
  ELSE
    ZDT=PDTPHY
  ENDIF
ENDIF

ZINVDT=1/PDTPHY

ZINVG=1._JPRB/RG 

! initialisation de ZDTMSE
IF (LLXFUMSE) THEN
  ZDTMSE=0.01_JPRB
  ZSTATI=REAL(RSTATI,JPRB)-ZDTMSE/2._JPRB
ELSE
  ZDTMSE=PDTPHY
  ZSTATI=REAL(RSTATI,JPRB)
ENDIF

IF(LTWOTL) THEN
  ZRHGMT=REAL(RHGMT,JPRB)-ZDTMSE/2._JPRB
ELSE
  ZRHGMT=REAL(RHGMT,JPRB)
ENDIF


LLMSE=LMSE.AND.(NSURFEXCTL >= 2)
LLMSE_PARAM=LLMSE
LLMSE_DIAG=LLMSE.AND.(NSURFEXCTL >= 3)


!  Vertical points
IKA=YDCPG_DIM%KFLEVG
IKB=YDCPG_DIM%KFLEVG
IKU=1
IKT=YDCPG_DIM%KFLEVG
IKTE=YDCPG_DIM%KFLEVG
IKTB=1
IKL=-1

!  SETUP

IF (INIT0 >= 0) THEN

  ZUM__(:,:)=ZVALUE
  ZUS__(:,:)=ZVALUE
  ZVM__(:,:)=ZVALUE
  ZVS__(:,:)=ZVALUE
  ZWM__(:,:)=ZVALUE
  ZTHSWAP__(:,:)=ZVALUE
  ZTHSAVE__(:,:)=ZVALUE
  ZSRCS__(:,:)=ZVALUE
  ZSIGS__(:,:)=ZVALUE
  ZTHM__(:,:)=ZVALUE_T
  ZRHODREFM__(:,:)=ZVALUE_ONE
  ZPABSM__(:,:)=ZVALUE_P
  ZTURB3D__(:,:,:)=ZVALUE
  ZLENGTHM__(:,:)=ZVALUE_L
  ZLENGTHH__(:,:)=ZVALUE_L

  ZZZ_(:,:)=ZVALUE
  ZMFM_(:,:)=ZVALUE
  ZSIGM_(:,:)=ZVALUE
  ZNEBMNH_(:,:)=ZVALUE
  ZEVAP_(:,:)=ZVALUE

  ZRSWAP_(:,:,:)=ZVALUE
  ZRSAVE_(:,:,:)=ZVALUE

  ZRM_(:,:,:)=ZVALUE

  ZLIMASWAP_(:,:,:)=ZVALUE
  ZLIMASAVE_(:,:,:)=ZVALUE

  ZLIMAM_(:,:,:)=ZVALUE

  ZFLXZTHVMF_(:,:)=ZVALUE
  ZSIGMF_(:,:)=ZVALUE
  ZRC_MF_(:,:)=ZVALUE
  ZRI_MF_(:,:)=ZVALUE
  ZCF_MF_(:,:)=ZVALUE

  ZSVSWAP_(:,:,:)=ZVALUE
  ZSVSAVE_(:,:,:)=ZVALUE

  ZSVMSWAP_(:,:,:)=ZVALUE
  ZSVMSAVE_(:,:,:)=ZVALUE
  ZSVMB_(:,:)=ZVALUE

  ZPIZA_DST_(:,:,:)  = ZVALUE
  ZCGA_DST_(:,:,:)   = ZVALUE
  ZTAUREL_DST_(:,:,:) = ZVALUE_EPSILON

  ZAERD_(:,:)=ZVALUE

  ZMFS_(:,:)=ZVALUE

  ZFPR(:,:,:)=ZVALUE

  IF(LKFBCONV) THEN
    ZCVTENDRV_(:,:)=ZVALUE
    ZCVTENDRC_(:,:)=ZVALUE
    ZCVTENDRI_(:,:)=ZVALUE
    ZCVTENDT_(:,:)=ZVALUE
  ENDIF

  ZACPRG_(:)=ZVALUE
  ZINPRR_NOTINCR_(:)=ZVALUE
  ZINPRS_NOTINCR_(:)=ZVALUE
  ZINPRG_NOTINCR_(:)=ZVALUE
  ZINPRH_NOTINCR_(:)=ZVALUE

  ZTHLS_(:,:)=ZVALUE
  ZMFUS_(:,:)=ZVALUE
  ZMFVS_(:,:)=ZVALUE

  ZTKEEDMF(:,:)=ZVALUE

  ZSFTH_(:)=ZVALUE
  ZSFRV_(:)=ZVALUE
  ZSFU_(:)=ZVALUE
  ZSFV_(:)=ZVALUE

  ZSFCO2_(:)=ZVALUE
  ZEMIS(:)=ZVALUE_ONE

  ZQICE(:,:)=ZVALUE
  ZQLIQ(:,:)=ZVALUE
  ZQO3(:,:)=ZVALUE

  ZAER(:,:,:)=ZVALUE
  ZAERINDS(:,:)=ZVALUE

  ZQSAT(:,:)=ZVALUE
  ZLH(:,:)=ZVALUE
  ZLSCPE(:,:)=ZVALUE
  ZGEOSLC(:,:)=ZVALUE

  ZFRSOFS(:)=ZVALUE

  ZQW(:,:)=ZVALUE
  YDCPG_MISC%RH(:,:)=ZVALUE
  ZTW(:,:)=ZVALUE

  ZTRSODIF(:,:)=ZVALUE
  ZTRSODIR(:,:)=ZVALUE
  ZZS_FSWDIR(:,:)=ZVALUE
  ZZS_FSWDIF(:,:)=ZVALUE

  ZSDUR(:)=ZVALUE
  ZDSRP(:)=ZVALUE
  ZALBD(:,:)=ZVALUE
  ZALBP(:,:)=ZVALUE
  ZALBD1(:)=ZVALUE
  ZALBP1(:)=ZVALUE

  ZQV(:,:)=ZVALUE
  ZSFSV_(:,:)=ZVALUE

  YDMF_PHYS%OUT%FRSOC(:,:)=ZVALUE
  YDMF_PHYS%OUT%FRTHC(:,:)=ZVALUE
  YDMF_PHYS%OUT%FRSOPS(:)=ZVALUE
  YDMF_PHYS%OUT%DIAGH(:)=ZVALUE

  ZTENDSV_TURBLIMA_(:,:,:)=ZVALUE

ENDIF

!  INITIALIZE (CUMULATED) TENDENCIES

DO JLEV=1,YDCPG_DIM%KFLEVG
  ZTENDT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=0.0_JPRB
  YDMF_PHYS%OUT%TENDU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=0.0_JPRB
  YDMF_PHYS%OUT%TENDV(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=0.0_JPRB
  ZTENDW(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=0.0_JPRB
  ZTENDTKE(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=0.0_JPRB
ENDDO
DO JRR=1,NRR
  DO JLEV=1,YDCPG_DIM%KFLEVG
    ZTENDR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV,JRR)=0.0_JPRB
  ENDDO
ENDDO
DO JGFL=1,NGFL_EXT
  DO JLEV=1,YDCPG_DIM%KFLEVG
    ZTENDEXT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV,JGFL)=0.0_JPRB
  ENDDO
ENDDO
DO JGFL=1,NLIMA
  DO JLEV=1,YDCPG_DIM%KFLEVG
    ZTENDLIMA(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV,JGFL)=0.0_JPRB
  ENDDO
ENDDO

!  INITIALIZE CUMULATED STUFF

! Small array, OK. REK
ZINPRH_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0._JPRB
ZINPRR_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0._JPRB
ZACPRR_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0._JPRB
ZINPRS_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0._JPRB
ZACPRS_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0._JPRB
ZINPRG_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0._JPRB


DO JLEV = 1,YDCPG_DIM%KFLEVG
  DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZZZ_F_(JLON,JLEV)=YLMF_PHYS_STATE%YCPG_DYN%PHIF(JLON,JLEV)*ZINVG
    ZTENDTT(JLON,JLEV)=0._JPRB
  ENDDO
ENDDO

! adhoc solution to avoid negative tke values
! when SL advective ddh is activated 
IF (LRSLDDH) THEN
  DO JLEV=1, YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
      ZTKEM4SLDDH(JLON,JLEV)=MAX(YLMF_PHYS_STATE%TKE(JLON,JLEV),PPTKEMIN)
    ENDDO
  ENDDO
  ZTKEM => ZTKEM4SLDDH(:,:)
ELSE
  ZTKEM => YLMF_PHYS_STATE%TKE(:,:)
  !test TKE > 0.
  IF (MINVAL(ZTKEM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)) <= 0._JPRB) THEN
    CALL ABOR1('TKE < 0 under APL_AROME check YTKE_NL%NREQIN')
  ENDIF
ENDIF


!test invalid combinations
IF (LHARATU .AND. CMF_UPDRAFT == 'EDKF') THEN
  CALL ABOR1('Combination LHARATU and EDKF not valid!')
ENDIF

!initialisation of first useful field for EZDIAG use in Chemistry/Dust
IOFF_MFSHAL=1
IF(LFPREC3D) IOFF_MFSHAL=2

!    ------------------------------------------------------------------
!     2 - CHANGEMENTS DE VARIABLES ET INVERSION DES NIVEAUX 
!     --------------------------------------------------------------------

IF (LMICRO.OR.LTURB.OR.LLMSE.OR.LKFBCONV) THEN

  ZRSCP=RD/RCPD
  ZINVATM=1/RATM
  !initialisation de ZZZ_
  DO JLEV = 1,YDCPG_DIM%KFLEVG
   !initialisation de qdm (utile localement pour calculer rho  
   !et convertir q en r 
    IF (NRR==7) THEN
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZQDM(JLON,JLEV)=1._JPRB-YLMF_PHYS_STATE%Q(JLON,JLEV)-YLMF_PHYS_STATE%L(JLON,JLEV)-YLMF_PHYS_STATE%R(JLON,JLEV)&
         & -YLMF_PHYS_STATE%I(JLON,JLEV)-YLMF_PHYS_STATE%S(JLON,JLEV)-YLMF_PHYS_STATE%G(JLON,JLEV)-YLMF_PHYS_STATE%H(JLON,JLEV)   
      ENDDO
    ELSE
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZQDM(JLON,JLEV)=1._JPRB-YLMF_PHYS_STATE%Q(JLON,JLEV)-YLMF_PHYS_STATE%L(JLON,JLEV)-YLMF_PHYS_STATE%R(JLON,JLEV)&
         & -YLMF_PHYS_STATE%I(JLON,JLEV)-YLMF_PHYS_STATE%S(JLON,JLEV)-YLMF_PHYS_STATE%G(JLON,JLEV)
      ENDDO
    ENDIF
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
   !initialisation de ZRHODREFM__ (=qd*zrho) 
      ZRHO=YLMF_PHYS_STATE%YCPG_PHY%PREF(JLON,JLEV)/(YLMF_PHYS_STATE%YCPG_DYN%RCP%R(JLON,JLEV)*YLMF_PHYS_STATE%T(JLON,JLEV))
      ZRHODREFM__(JLON,JLEV)=ZRHO*ZQDM(JLON,JLEV)
      ZRHODJM__(JLON,JLEV)=YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
  !initialisation de ZEXNREFM_
      ZEXNREFM_(JLON,JLEV)=(YLMF_PHYS_STATE%YCPG_PHY%PREF(JLON,JLEV)*ZINVATM)**(ZRSCP)
 ! vent horizontal et TKE
      ZPABSM__(JLON,JLEV)=YLMF_PHYS_STATE%YCPG_PHY%PREF(JLON,JLEV)
      ZUM__(JLON,JLEV)= YLMF_PHYS_STATE%U(JLON,JLEV)
      ZVM__(JLON,JLEV)= YLMF_PHYS_STATE%V(JLON,JLEV)
      ZWM__(JLON,JLEV)= YLMF_PHYS_STATE%YCPG_PHY%W(JLON,JLEV)
      ZTKEM__(JLON,JLEV)= ZTKEM(JLON,JLEV)
      ZZZ_(JLON,JLEV)=YLMF_PHYS_STATE%YCPG_DYN%PHI(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO
  !initialise sigma for subgrid condensation coming
  !from previous time step turbulence scheme
  IF (LOSIGMAS) THEN
    DO JLEV = 1, YDCPG_DIM%KFLEVG 
      ZSIGM_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)= YLMF_PHYS_STATE%SRC(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)
    ENDDO
  ENDIF
  !initialise convective mas flux for subgrid condensation coming 
  !from previous time step convection scheme
  IF ( LKFBCONV.AND.LOSUBG_COND.AND..NOT.LOSIGMAS) THEN
    DO JLEV = 1, YDCPG_DIM%KFLEVG 
      ZMFM_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=YLMF_PHYS_STATE%SRC(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV) 
    ENDDO
  ENDIF
!!! initialisation des variables d etat MNH �t

  !initialisation de ZRM_ pour les hydrometeores (ri=qi/qd)
  DO JLEV = 1, YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZTHM__(JLON,JLEV)=YLMF_PHYS_STATE%T(JLON,JLEV)/ZEXNREFM_(JLON,JLEV)
      ZRM_(JLON,JLEV,1)=YLMF_PHYS_STATE%Q(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ZRM_(JLON,JLEV,2)=YLMF_PHYS_STATE%L(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ZRM_(JLON,JLEV,3)=YLMF_PHYS_STATE%R(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ZRM_(JLON,JLEV,4)=YLMF_PHYS_STATE%I(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ZRM_(JLON,JLEV,5)=YLMF_PHYS_STATE%S(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ZRM_(JLON,JLEV,6)=YLMF_PHYS_STATE%G(JLON,JLEV)/ZQDM(JLON,JLEV)  
    ENDDO
  ENDDO
  
  IF (NRR==7) THEN
    DO JLEV = 1, YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZRM_(JLON,JLEV,7)=YLMF_PHYS_STATE%H(JLON,JLEV)/ZQDM(JLON,JLEV)  
      ENDDO
    ENDDO
  ENDIF
   
  IF (NRR==6) THEN
    !initialisation de ZTHVREFM__
    DO JLEV = 1, YDCPG_DIM%KFLEVG 
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZTHVREFM__(JLON,JLEV)=ZTHM__(JLON,JLEV)*&
         & (1._JPRB+ZRM_(JLON,JLEV,1)*(RV/RD))/&
         & (1._JPRB+ZRM_(JLON,JLEV,1)+ZRM_(JLON,JLEV,2) +&
         & ZRM_(JLON,JLEV,3)+ZRM_(JLON,JLEV,4)+&
         & ZRM_(JLON,JLEV,5)+ZRM_(JLON,JLEV,6))
      ENDDO
    ENDDO
  ELSEIF (NRR==7) THEN
    DO JLEV = 1, YDCPG_DIM%KFLEVG 
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZTHVREFM__(JLON,JLEV)=ZTHM__(JLON,JLEV)*&
         & (1._JPRB+ZRM_(JLON,JLEV,1)*(RV/RD))/&
         & (1._JPRB+ZRM_(JLON,JLEV,1)+ZRM_(JLON,JLEV,2) +&
         & ZRM_(JLON,JLEV,3)+ZRM_(JLON,JLEV,4)+&
         & ZRM_(JLON,JLEV,5)+ZRM_(JLON,JLEV,6)+&
         &  ZRM_(JLON,JLEV,7) )
      ENDDO
    ENDDO
  ENDIF

!!! initialisation des variables d etat MNH a t+dt
!!! division pas le pas de temps
!!!(la multiplication par rhodj est faite plus tard, si necessaire, 
!!! suivant les parametrisations)   

  ! initialise pointers :
  CALL SWAP_THS
  ! vent horizontal
  DO JLEV = 1, YDCPG_DIM%KFLEVG 
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZUS__(JLON,JLEV)= YLMF_PHYS_STATE%U(JLON,JLEV)*ZINVDT
      ZVS__(JLON,JLEV)= YLMF_PHYS_STATE%V(JLON,JLEV)*ZINVDT
      ZWS__(JLON,JLEV)= YLMF_PHYS_STATE%YCPG_PHY%W(JLON,JLEV)*ZINVDT
      ZTKES_(JLON,JLEV)= ZTKEM(JLON,JLEV)*ZINVDT
      ZTHS__(JLON,JLEV)=ZTHM__(JLON,JLEV)*ZINVDT
    ENDDO
  ENDDO

  !initialisation de ZRS_ pour les hydrometeores
  ! initialise pointers :
  CALL SWAP_RS
  DO JRR=1,NRR 
    DO JLEV = 1, YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZRS_(JLON,JLEV,JRR)=ZRM_(JLON,JLEV,JRR)*ZINVDT
      ENDDO
    ENDDO
  ENDDO

!!! Initialisations temporaires d'arguments non-utilises
  !initialisation de ZCIT_
  ZCIT_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:IKT)=0.0_JPRB
  
  !initialisation des tableaux de precipitations inst. and cumulated 
  !and surface fluxes for turbulence
  IF (LLMSE.OR.LSFORCS) THEN
    ZACPRR_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_GPAR%ACPRR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
    ZACPRS_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_GPAR%ACPRS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
    ZACPRG_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_GPAR%ACPRG(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
    ZINPRR_NOTINCR_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_GPAR%INPRR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
    ZINPRS_NOTINCR_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_GPAR%INPRS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
    ZINPRG_NOTINCR_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_GPAR%INPRG(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
  ENDIF

  !initialisation des scalaires passifs
  ! initialise pointers :
  CALL SWAP_SVM
  CALL SWAP_SVS
  DO JGFL=1,NGFL_EXT
    DO JLEV = 1, YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZSVM_(JLON,JLEV,JGFL)=YLMF_PHYS_STATE%P1EXT(JLON,JLEV,JGFL)
        ZSVS_(JLON,JLEV,JGFL)=YLMF_PHYS_STATE%P1EXT(JLON,JLEV,JGFL)*ZINVDT
      ENDDO
    ENDDO
  ENDDO

  !initialisation des concentrations LIMA
  ! initialise pointers :
  CALL SWAP_LIMAS
  DO JGFL=1,NLIMA
    DO JLEV = 1, YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZLIMAM_(JLON,JLEV,JGFL)=YLMF_PHYS_STATE%P1LIMA(JLON,JLEV,JGFL)
        ZLIMAS_(JLON,JLEV,JGFL)=YLMF_PHYS_STATE%P1LIMA(JLON,JLEV,JGFL)*ZINVDT
      ENDDO
    ENDDO
  ENDDO

ENDIF

! daand: radflex
ZFRSO => YDMF_PHYS%OUT%FRSO(:,:,1)
ZFRTH => YDMF_PHYS%OUT%FRTH(:,:,1)

!    ------------------------------------------------------------------
!     3 - PRINTS FOR DIAGNOSTICS IF NEEDED
!    ------------------------------------------------------------------
IF (LDIAGWMAX) THEN
  IF (MOD(YDCPG_DIM%KSTEP+1,NDIAGWMAX)==0) THEN
  ! calcul de wmax
    DO JLEV = 1 , YDCPG_DIM%KFLEVG
      Z_WMAX=0._JPRB
      Z_WMIN=0._JPRB
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        IF (YLMF_PHYS_STATE%YCPG_PHY%W(JLON,JLEV)>Z_WMAX) THEN
          Z_WMAX=YLMF_PHYS_STATE%YCPG_PHY%W(JLON,JLEV)
        ENDIF
        IF (YLMF_PHYS_STATE%YCPG_PHY%W(JLON,JLEV)<Z_WMIN) THEN
          Z_WMIN=YLMF_PHYS_STATE%YCPG_PHY%W(JLON,JLEV)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
ENDIF

IF (LFLEXDIA) THEN
  !save tendencies
  ZTENDTBAK(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=ZTENDT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  DO JR=1,NRR
    ZTENDRBAK(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,JR)=ZTENDR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,JR)
  ENDDO
ENDIF


!    ------------------------------------------------------------------
!     4 - ADJUSTMENT (CALLED IF THE MICROPHYSICS IS SWITCH ON) 
!    ------------------------------------------------------------------

IF (LMICRO) THEN

  ! Swap pointers because input values of THS and RS should be saved
  CALL SWAP_THS
  CALL SWAP_RS

  IF (LMFSHAL .AND. (CMF_CLOUD=='DIRE'.OR.CMF_CLOUD=='BIGA')) THEN
    IOFF_MFSHAL=IOFF_MFSHAL+3
    IF (YDCPG_DIM%KSTEP==0) THEN
      DO JLEV = 1, YDCPG_DIM%KFLEVG 
        ZRC_MF_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=0._JPRB
        ZRI_MF_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=0._JPRB
        ZCF_MF_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=0._JPRB
      ENDDO
    ELSE
      DO JLEV = 1, YDCPG_DIM%KFLEVG 
        ZRC_MF_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=ZP1EZDIAG(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV,1)
        ZRI_MF_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=ZP1EZDIAG(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV,3)
        ZCF_MF_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=ZP1EZDIAG(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV,2)
      ENDDO
    ENDIF
    ZP1EZDIAG(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,1:3)=0._JPRB
  ENDIF

  IF (MOD(YDCPG_DIM%KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'avant aro_adjust sous apl_arome'
    WRITE(NULOUT,*)'JLEV   ZZZ_F_   RHODJM   EXNREFM       PABSM       THM      SIGM         MFM    '
    DO JLEV=1,YDCPG_DIM%KFLEVG+1 
      WRITE(NULOUT,'(I2,X,7F10.3)') JLEV,ZZZ_F_(NPTP,JLEV),ZRHODJM__(NPTP,JLEV),&
       & ZEXNREFM_(NPTP,JLEV), ZPABSM__(NPTP,JLEV), ZTHM__(NPTP,JLEV), ZSIGM_(NPTP,JLEV), ZMFM_(NPTP,JLEV)
    ENDDO 
    WRITE (NULOUT,*)'JLEV  rhoQv  rhoQc   rhoQr   rhoQi   rhoQs   rhoQg'
    DO JLEV=1,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,'(I2,X,6E11.4)') JLEV,ZRM_(NPTP,JLEV,1),&
       & ZRM_(NPTP,JLEV,2), ZRM_(NPTP,JLEV,3),ZRM_(NPTP,JLEV,4),ZRM_(NPTP,JLEV,5), ZRM_(NPTP,JLEV,6)  
    ENDDO
    WRITE (NULOUT,*)'JLEV  ZRC_MF_  ZRI_MF_  ZCF_MF_ ZTHSIN_  ZSRCS__ ZNEBMNH_'
    DO JLEV=1,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,'(I2,X,6E11.4)') JLEV,ZRC_MF_(NPTP,JLEV),&
       & ZRI_MF_(NPTP,JLEV),ZCF_MF_(NPTP,JLEV), ZTHSIN_(NPTP,JLEV),ZSRCS__(NPTP,JLEV), ZNEBMNH_(NPTP,JLEV)
    ENDDO
  ENDIF

  ! for now a copy is needed (see below, inside). I don't like than :-( REK
  ZTHS__(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=ZTHSIN_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  ZRS_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,1:NRR)=ZRSIN_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,1:NRR)

  IF (CMICRO == 'LIMA') THEN

    CALL SWAP_LIMAS
    ! for now a copy is needed (see below, inside). I don't like than :-( REK
    ZLIMAS_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,1:NLIMA)=ZLIMASIN_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,1:NLIMA)

    CALL ARO_ADJUST_LIMA (YDCPG_DIM%KFLEVG, IKU, IKL, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, NRR, NLIMA,                 &
    & YDCPG_DIM%KSTEP+1, LOSUBG_COND, LOSIGMAS, LOCND2, ZDT, VSIGQSAT, ZZZ_F_, ZRHODJM__(:, 1:YDCPG_DIM%KFLEVG),     &
    & ZRHODREFM__(:, 1:YDCPG_DIM%KFLEVG), ZEXNREFM_, ZPABSM__(:, 1:YDCPG_DIM%KFLEVG), ZTHM__(:, 1:YDCPG_DIM%KFLEVG), &
    & ZRM_, ZLIMAM_, ZSIGM_, ZMFM_, ZRC_MF_, ZRI_MF_, ZCF_MF_, ZTHS__(:, 1:YDCPG_DIM%KFLEVG), ZRS_,                  &
    & ZLIMAS_, ZSRCS__(:, 1:YDCPG_DIM%KFLEVG), ZNEBMNH_, YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH   &
    & )
  ELSE

    CALL ARO_ADJUST (YDCPG_DIM%KFLEVG, IKU, IKL, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, NRR, NGFL_EZDIAG,                   &
    & CFRAC_ICE_ADJUST, LOSUBG_COND, LOSIGMAS, CMICRO, LOCND2, ZDT, VSIGQSAT, ZZZ_F_, ZRHODJM__(:, 1:YDCPG_DIM%KFLEVG), &
    & ZEXNREFM_, ZPABSM__(:, 1:YDCPG_DIM%KFLEVG), ZTHM__(:, 1:YDCPG_DIM%KFLEVG), ZRM_, ZSIGM_, ZMFM_,                   &
    & ZRC_MF_, ZRI_MF_, ZCF_MF_, ZTHS__(:, 1:YDCPG_DIM%KFLEVG), ZRS_, ZSRCS__(:, 1:YDCPG_DIM%KFLEVG),                   &
    & ZNEBMNH_, ZGP2DSPP, ZP1EZDIAG, YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH)

  ENDIF
  
  IF (MOD(YDCPG_DIM%KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'apres aro_adjust sous apl_arome'
    WRITE(NULOUT,*)'JLEV   ZZZ_F_   RHODJM   EXNREFM       PABSM       THM      SIGM         MFM    '   
    DO JLEV=1,YDCPG_DIM%KFLEVG+1 
      WRITE(NULOUT,'(I2,X,7F10.3)') JLEV,ZZZ_F_(NPTP,JLEV),ZRHODJM__(NPTP,JLEV),&
       & ZEXNREFM_(NPTP,JLEV), ZPABSM__(NPTP,JLEV), ZTHM__(NPTP,JLEV), ZSIGM_(NPTP,JLEV), ZMFM_(NPTP,JLEV)
    ENDDO 
    WRITE (NULOUT,*)'JLEV  rhoQv  rhoQc   rhoQr   rhoQi   rhoQs   rhoQg'
    DO JLEV=1,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,'(I2,X,6E11.4)') JLEV,ZRS_(NPTP,JLEV,1),&
       & ZRS_(NPTP,JLEV,2), ZRS_(NPTP,JLEV,3),ZRS_(NPTP,JLEV,4),ZRS_(NPTP,JLEV,5), ZRS_(NPTP,JLEV,6)  
    ENDDO
    WRITE (NULOUT,*)'JLEV  ZRC_MF_  ZRI_MF_  ZCF_MF_ ZTHS__ ZSRCS__ ZNEBMNH_'
    DO JLEV=1,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,'(I2,X,6E11.4)')JLEV,ZRC_MF_(NPTP,JLEV),&
       & ZRI_MF_(NPTP,JLEV),ZCF_MF_(NPTP,JLEV), ZTHS__(NPTP,JLEV),ZSRCS__(NPTP,JLEV), ZNEBMNH_(NPTP,JLEV)
    ENDDO
  ENDIF

  DO JLEV=1,YDCPG_DIM%KFLEVG
    YDVARS%A%T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=ZNEBMNH_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV) 
  ENDDO

  !adjusted zthm and zrm
  DO JLEV = 1, YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZTHM__(JLON,JLEV)=ZTHS__(JLON,JLEV)*PDTPHY
    ENDDO
  ENDDO

  DO JRR=1,NRR
    DO JLEV = 1, YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZRM_(JLON,JLEV,JRR)=ZRS_(JLON,JLEV,JRR)*PDTPHY
      ENDDO
    ENDDO
  ENDDO


  !initialisation de qdm utile pour 
  !convertir tendance de r en tendance de q 
  IF (NRR==6) THEN
    DO JLEV=1,YDCPG_DIM%KFLEVG
      DO JLON= YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
        ZQDM(JLON,JLEV)=1._JPRB/(1._JPRB+ZRM_(JLON,JLEV,1)+&
        &ZRM_(JLON,JLEV,2)+ZRM_(JLON,JLEV,3)+ZRM_(JLON,JLEV,4)+ZRM_(JLON,JLEV,5)+&
        &ZRM_(JLON,JLEV,6) )
      ENDDO
    ENDDO
  ELSEIF (NRR==7) THEN
    DO JLEV=1,YDCPG_DIM%KFLEVG
      DO JLON= YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
        ZQDM(JLON,JLEV)=1._JPRB/(1._JPRB+ZRM_(JLON,JLEV,1)+&
        &ZRM_(JLON,JLEV,2)+ZRM_(JLON,JLEV,3)+ZRM_(JLON,JLEV,4)+ZRM_(JLON,JLEV,5)+&
        &ZRM_(JLON,JLEV,6)+ZRM_(JLON,JLEV,7) )
      ENDDO
    ENDDO
  ENDIF 
  !reinitialisation des qi
  DO JLEV = 1, YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZQVM(JLON,JLEV)=ZRM_(JLON,JLEV,1)*ZQDM(JLON,JLEV)
      ZQCM(JLON,JLEV)=ZRM_(JLON,JLEV,2)*ZQDM(JLON,JLEV)
      ZQRM(JLON,JLEV)=ZRM_(JLON,JLEV,3)*ZQDM(JLON,JLEV)
      ZQIM(JLON,JLEV)=ZRM_(JLON,JLEV,4)*ZQDM(JLON,JLEV)
      ZQSM(JLON,JLEV)=ZRM_(JLON,JLEV,5)*ZQDM(JLON,JLEV)
      ZQGM(JLON,JLEV)=ZRM_(JLON,JLEV,6)*ZQDM(JLON,JLEV)
    ENDDO
  ENDDO

  IF (NRR==7) THEN
    DO JLEV = 1,YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZQHM(JLON,JLEV)=ZRM_(JLON,JLEV,7)*ZQDM(JLON,JLEV)
      ENDDO
    ENDDO
  ELSE
    ZQHM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=0._JPRB
  ENDIF

  ! Tendances des variables LIMA
  DO JGFL=1,NLIMA
    DO JLEV = 1, YDCPG_DIM%KFLEVG
      ! Réinitialisation des variables LIMA
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZLIMAM_(JLON,JLEV,JGFL)=ZLIMAS_(JLON,JLEV,JGFL)*PDTPHY
        ZTENDLIMA(JLON,JLEV,JGFL)=ZTENDLIMA(JLON,JLEV,JGFL)+(ZLIMAS_(JLON,JLEV,JGFL)-ZLIMASIN_(JLON,JLEV,JGFL))  
      ENDDO
    ENDDO
  ENDDO

  !modif de R et CP
  ZQHGM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)=ZQHM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)+ZQGM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
  CALL GPRCP(YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, PQ=ZQVM, PQI=ZQIM, PQL=ZQCM, &
  & PQR=ZQRM, PQS=ZQSM, PQG=ZQHGM, PCP=ZCPM, PR=ZRHM)  

  DO JLEV = 1,YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZTM(JLON,JLEV)=ZTHM__(JLON,JLEV)*ZEXNREFM_(JLON,JLEV)
  !reinitialisation de ZRHODREFM__ (=qd*zrho)
      ZRHO=YLMF_PHYS_STATE%YCPG_PHY%PREF(JLON,JLEV)/(ZRHM(JLON,JLEV)*ZTM(JLON,JLEV))
      ZRHODREFM__(JLON,JLEV)=ZRHO*ZQDM(JLON,JLEV)
    ENDDO
  ENDDO

  !geopotentiel calculation
 
  ZAPHIM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)=YLMF_PHYS_STATE%YCPG_DYN%PHI(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
  CALL GPGEO(YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, ZAPHIM, ZAPHIFM,          &
  & ZTM, ZRHM, YLMF_PHYS_STATE%YCPG_PHY%XYB%LNPR, YLMF_PHYS_STATE%YCPG_PHY%XYB%ALPH, YDGEOMETRY%YRVERT_GEOM&
  & )
   
  !calcul de l'altitude
  DO JLEV = 1, YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZZZ_(JLON,JLEV)=ZAPHIM(JLON,JLEV)*ZINVG
      !initialisation de ZZZ_F_
      ZZZ_F_(JLON,JLEV)=ZAPHIFM(JLON,JLEV)*ZINVG
      ! tendency of T
      ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV)+(ZTHS__(JLON,JLEV)-ZTHSIN_(JLON,JLEV))*ZEXNREFM_(JLON,JLEV)
      ZTENDTT(JLON,JLEV)=ZTHS__(JLON,JLEV)-ZTHSIN_(JLON,JLEV)
    ENDDO
  ENDDO
  !inversion niveaux tendances des ri et conversion en qi en multipliant par qd 
  DO JR=1,NRR
    DO JLEV=1,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZTENDR(JLON,JLEV,JR)=ZTENDR(JLON,JLEV,JR)+(ZRS_(JLON,JLEV,JR)-ZRSIN_(JLON,JLEV,JR))*ZQDM(JLON,JLEV)
      ENDDO
    ENDDO
  ENDDO
  !initialisation de ZDZZ_
  DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZDZZ_(JLON,1)=ZAPHIM(JLON,0)*ZINVG-ZZZ_(JLON,1)
  ENDDO
  DO JLEV = 2, YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZDZZ_(JLON,JLEV)=ZZZ_(JLON,JLEV+IKL)-ZZZ_(JLON,JLEV)
    ENDDO
  ENDDO

ELSE

  ZTM (YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YLMF_PHYS_STATE%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  ZRHM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YLMF_PHYS_STATE%YCPG_DYN%RCP%R(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  ZQVM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YLMF_PHYS_STATE%Q(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  ZQIM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YLMF_PHYS_STATE%I(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  ZQCM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YLMF_PHYS_STATE%L(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  ZQRM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YLMF_PHYS_STATE%R(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  ZQSM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YLMF_PHYS_STATE%S(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  ZQGM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YLMF_PHYS_STATE%G(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  IF (NRR==7) THEN
    ZQHM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YLMF_PHYS_STATE%H(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  ELSE
    ZQHM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=0._JPRB
  ENDIF
  ZCPM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)=YLMF_PHYS_STATE%YCPG_DYN%RCP%CP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KTDIA:YDCPG_DIM%KFLEVG)
  ZAPHIM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG)=YLMF_PHYS_STATE%YCPG_DYN%PHI(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,0:YDCPG_DIM%KFLEVG)
  ZAPHIFM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YLMF_PHYS_STATE%YCPG_DYN%PHIF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  ZZZ_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YLMF_PHYS_STATE%YCPG_DYN%PHI(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)*ZINVG
  !initialisation of PCLFS outside LMICRO to be zero in case LMICRO=F
  YDVARS%A%T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=0._JPRB

ENDIF ! ADJUSTMENT LMICRO

IF (LFLEXDIA) THEN
  DO JLEV = 1, YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA 
      ZTMPAF(JLON,JLEV)=(ZTENDT(JLON,JLEV)-ZTENDTBAK(JLON,JLEV))*YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG*ZCPM(JLON,JLEV)
    ENDDO
  ENDDO
  IF (LDDH_OMP) THEN
    CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZTMPAF, 'TCTADJU', YDDDH)
  ELSE
    CALL ADD_FIELD_3D(YLDDH, ZTMPAF, 'TCTADJU', 'T', 'ARP', .TRUE., .TRUE.)
  ENDIF
  DO JR=1,NRR
    CLNAME='T'//CLVARNAME(JR)//'ADJU'
    DO JLEV = 1, YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA 
        ZTMPAF(JLON,JLEV)=(ZTENDR(JLON,JLEV,JR)-ZTENDRBAK(JLON,JLEV,JR))*YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
      ENDDO
    ENDDO
    IF (LDDH_OMP) THEN
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZTMPAF, CLNAME, YDDDH)
    ELSE
      CALL ADD_FIELD_3D(YLDDH, ZTMPAF, CLNAME, 'T', 'ARP', .TRUE., .TRUE.)
    ENDIF
    DO JLEV = 1, YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA 
        ZTMPAF(JLON,JLEV)=YDVARS%A%T1(JLON,JLEV)*YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)
      ENDDO
    ENDDO
    IF (LDDH_OMP) THEN
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZTMPAF, 'VNT', YDDDH)
    ELSE
      CALL ADD_FIELD_3D(YLDDH, ZTMPAF, 'VNT', 'V', 'ARP', .TRUE., .TRUE.)
    ENDIF
  ENDDO
! specific to new data flow for diagnostics
  IF (LDDH_OMP) THEN
    DO JLEV = 1, YDCPG_DIM%KFLEVG
      ZCON1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV) = 1.0_JPRB
      ZCON2(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV) = ZQDM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)
    ENDDO
    DO JLEV = 1, YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA 
        ZCON3(JLON,JLEV) = YLMF_PHYS_STATE%YCPG_DYN%RCP%CP(JLON,JLEV)*ZEXNREFM_(JLON,JLEV)
      ENDDO
    ENDDO
    ! missing interface !!! REK
    CALL ARO_SUINTBUDGET_OMP(YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YDCPG_DIM%KSTEP, ZCON1, ZCON2, ZCON3, &
    & YDDDH)
  ELSE
    ! missing interface !!! REK
    CALL ARO_SUINTBUDGET(YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YDCPG_DIM%KSTEP, YDCPG_DIM%KFDIA, ZQDM, &
    & ZEXNREFM_, YLMF_PHYS_STATE%YCPG_DYN%RCP%CP)
  ENDIF
ENDIF


DO JLEV = 1, YDCPG_DIM%KFLEVG-1
  DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZDZZ_F_(JLON,JLEV)=ZZZ_F_(JLON,JLEV)-ZZZ_F_(JLON,JLEV-IKL)
  ENDDO
ENDDO
DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
  ZDZZ_F_(JLON,YDCPG_DIM%KFLEVG)=ZZZ_F_(JLON,YDCPG_DIM%KFLEVG)-YDOROG%OROG(JLON)*ZINVG
ENDDO


!     --------------------------------------------------------------------
!     5 - COMPUTE DUST PROPERTIES FOR RADIATION IF LRDUST=T
!     --------------------------------------------------------------------
IF (LRDUST) THEN
  ZP1EZDIAG(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,IOFF_MFSHAL:NGFL_EZDIAG)=0.0_JPRB
  ! input dust scalar concentration in ppp from
  IEZDIAG_CHEM=NGFL_EZDIAG-IOFF_MFSHAL+1
  ! Swapp because IN and OUT will be needed simultaneously
  CALL SWAP_SVM
  ! input dust scalar concentration in ppp from
  CALL ARO_MNHDUST (IKL, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, NGFL_EXT, PDTPHY, ZSVMIN_, ZZZ_, ZDZZ_, ZPABSM__(:, 1:YDCPG_DIM%KFLEVG), &
  & ZTHM__(:, 1:YDCPG_DIM%KFLEVG), ZRHODREFM__(:, 1:YDCPG_DIM%KFLEVG), NSWB_MNH, YDCPG_DIM%KSTEP+1,                                  &
  & ZSVM_, ZPIZA_DST_, ZCGA_DST_, ZTAUREL_DST_, ZAERD_, IEZDIAG_CHEM, ZPEZDIAG_(:, :, IOFF_MFSHAL:NGFL_EZDIAG)                       &
  &            )
  ZP1EZDIAG(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,IOFF_MFSHAL:NGFL_EZDIAG)=ZPEZDIAG_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,IOFF_MFSHAL:NGFL_EZDIAG)
! return to tendency
  DO JGFL=1, NGFL_EXT
    DO JLEV = 1,YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZTENDEXT(JLON,JLEV,JGFL)=ZTENDEXT(JLON,JLEV,JGFL)+(ZSVM_(JLON,JLEV,JGFL)-ZSVMIN_(JLON,JLEV,JGFL))*ZINVDT
      ENDDO
    ENDDO
  ENDDO
ENDIF ! LRDUST

IF (LSFORCS) THEN   ! <== Surface forcing for MUSC

  ZTSURF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA) = YLMF_PHYS_STATE%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
  ZTN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)    = YLMF_PHYS_STATE%T(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
  ZQS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)    = YLMF_PHYS_STATE%Q(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZRHODREFM(JLON) = YLMF_PHYS_STATE%YCPG_PHY%PREF(JLON,YDCPG_DIM%KFLEVG)/(YLMF_PHYS_STATE%T(JLON,YDCPG_DIM%KFLEVG)*YLMF_PHYS_STATE%YCPG_DYN%RCP%R(JLON,YDCPG_DIM%KFLEVG))
    ZTHETAS(JLON)   = ZTSURF(JLON)*(RATM/YLMF_PHYS_STATE%YCPG_PHY%PRE(JLON,YDCPG_DIM%KFLEVG))**RKAPPA
  ENDDO

  LLAROME=.TRUE.
  CALL SURF_IDEAL_FLUX(YDRIP, YDPHY0, YDPHYDS, LLAROME, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,                    &
  & YLMF_PHYS_STATE%YCPG_DYN%PHIF(:, YDCPG_DIM%KFLEVG), ZRHODREFM, YDMF_PHYS_SURF%GSD_SFO%PGROUP,                            &
  & ZTN, ZTSURF, YDMF_PHYS_SURF%GSD_VF%PLSM, YLMF_PHYS_STATE%Q(:, YDCPG_DIM%KFLEVG), YLMF_PHYS_STATE%U(:, YDCPG_DIM%KFLEVG), &
  & YLMF_PHYS_STATE%V(:, YDCPG_DIM%KFLEVG), ZTHETAS, ZSFTH_, ZSFRV_, ZSFU_, ZSFV_)

!* Compute PBL-diagnostics
   
   ZCAPE(:)=0._JPRB
   ZDCAPE(:)=0._JPRB
  CALL ACCLDIA(YDXFU, YDPHY, YDMODEL%YRML_PHY_MF%YRPHY2, YDTOPH, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,                               &
  & YDCPG_DIM%KFLEVG, YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YLMF_PHYS_STATE%U(:, 1:YDCPG_DIM%KFLEVG),                                          &
  & YLMF_PHYS_STATE%V(:, 1:YDCPG_DIM%KFLEVG), ZCAPE, ZDCAPE, ZTKEM(:, 1:YDCPG_DIM%KFLEVG), YLMF_PHYS_STATE%YCPG_DYN%PHIF(:, 1:YDCPG_DIM%KFLEVG), &
  & YDOROG%OROG, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST, YDMF_PHYS%OUT%CLPH, ICLPH)

  YDMF_PHYS%OUT%CLPH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=MIN(XMAXLM,MAX(XMINLM,YDMF_PHYS%OUT%CLPH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA))) 

ENDIF    ! <== End of surface forcing for MUSC

!     --------------------------------------------------------------------
!     6 - RADIATION LRAYFM (IFS) or LRAY (ACRANEB2) 
!     --------------------------------------------------------------------
IF (LRAYFM.OR.LRAY) THEN
  ! prepare some input for both radiation schemes at every time step

  ! test de coherence sur le nombre de bandes spectrales entre ce qui sort de
  ! la surface et ce qu'attend le rayonnement
  IF( NSWB_MNH /= NSW) THEN
    CALL ABOR1 (' NSWB_MNH must be equal to NSW !')
  ENDIF

  ! compute saturated specific humidity
  CALL ACTQSAT ( YDPHY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTQSAT, YDCPG_DIM%KFLEVG, YLMF_PHYS_STATE%YCPG_PHY%PREF, &
  & ZCPM, ZQVM, ZTM, ZGEOSLC, ZLH, ZLSCPE, ZQSAT, ZQW, YDCPG_MISC%RH, ZTW)  

  IF (YSPP_CONFIG%LSPP.AND.YSPP_CONFIG%LPERT_RADGR) THEN
   IF (YSPP_CONFIG%LLNN_MEAN1.OR.YSPP_CONFIG%LLNN_MEAN1_RADGR) THEN
    ZMU = -0.5_JPRB * (YSPP_CONFIG%CMPERT_RADGR * YSPP_CONFIG%SDEV)**2
   ELSE
    ZMU = 0._JPRB
   ENDIF
   DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZVAL = RADGR*EXP(ZMU+YSPP_CONFIG%CMPERT_RADGR*ZGP2DSPP(JLON,YSPP%MP_RADGR))
    ZRADGR(JLON) = MAX(YSPP_CONFIG%CLIP_RADGR(1),MIN(ZVAL,YSPP_CONFIG%CLIP_RADGR(2)))
   ENDDO
   IF ( YSPP_CONFIG%IEZDIAG_POS > 0 ) THEN
    JKO=2*YSPP%MP_RADGR-1
    JKE=2*YSPP%MP_RADGR
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
     ZP1EZDIAG(JLON,JKO,YSPP_CONFIG%IEZDIAG_POS) = ZGP2DSPP(JLON,YSPP%MP_RADGR)
     ZP1EZDIAG(JLON,JKE,YSPP_CONFIG%IEZDIAG_POS) = ZRADGR(JLON)
    ENDDO
   ENDIF
  ELSE
   DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZRADGR(JLON) = RADGR
   ENDDO
  ENDIF

  IF (YSPP_CONFIG%LSPP.AND.YSPP_CONFIG%LPERT_RADSN) THEN
   IF (YSPP_CONFIG%LLNN_MEAN1.OR.YSPP_CONFIG%LLNN_MEAN1_RADSN) THEN
    ZMU = -0.5_JPRB * (YSPP_CONFIG%CMPERT_RADSN * YSPP_CONFIG%SDEV)**2
   ELSE
    ZMU = 0._JPRB
   ENDIF
   DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZVAL = RADSN*EXP(ZMU+YSPP_CONFIG%CMPERT_RADSN*ZGP2DSPP(JLON,YSPP%MP_RADSN))
    ZRADSN(JLON) = MAX(YSPP_CONFIG%CLIP_RADSN(1),MIN(ZVAL,YSPP_CONFIG%CLIP_RADSN(2)))
   ENDDO
   IF ( YSPP_CONFIG%IEZDIAG_POS > 0 ) THEN
    JKO=2*YSPP%MP_RADSN-1
    JKE=2*YSPP%MP_RADSN
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
     ZP1EZDIAG(JLON,JKO,YSPP_CONFIG%IEZDIAG_POS) = ZGP2DSPP(JLON,YSPP%MP_RADSN)
     ZP1EZDIAG(JLON,JKE,YSPP_CONFIG%IEZDIAG_POS) = ZRADSN(JLON)
    ENDDO
   ENDIF
  ELSE
   DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZRADSN(JLON) = RADSN
   ENDDO
  ENDIF
  ! initialisation des humidite (dans le rayonnement, l'eau liquide nuageuse 
  ! et la glace sont donne par des hu par rapport au gaz.
  ! (qi/qa+qv pour ice par ex. C'est donc different de ri)
  DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA 
       ZQICE(JLON,JLEV)= MAX(0.0_JPRB,&
        & (ZQIM(JLON,JLEV) + ZQSM(JLON,JLEV)*ZRADSN(JLON) + ZQGM(JLON,JLEV)*ZRADGR(JLON))/&
        & (1.0_JPRB-ZQIM(JLON,JLEV)-ZQCM(JLON,JLEV)-ZQRM(JLON,JLEV)&
        & -ZQGM(JLON,JLEV)-ZQSM(JLON,JLEV)-ZQHM(JLON,JLEV)))
       ZQLIQ(JLON,JLEV)=MAX(0.0_JPRB, ZQCM(JLON,JLEV)/&
        & (1.0_JPRB-ZQIM(JLON,JLEV)-ZQCM(JLON,JLEV)-ZQRM(JLON,JLEV)&
        & -ZQGM(JLON,JLEV)-ZQSM(JLON,JLEV)-ZQHM(JLON,JLEV)))
       ZQV(JLON,JLEV)=MAX(0.0_JPRB, ZQVM(JLON,JLEV)/&
        & (1.0_JPRB-ZQIM(JLON,JLEV)-ZQCM(JLON,JLEV)-ZQRM(JLON,JLEV)&
        & -ZQGM(JLON,JLEV)-ZQSM(JLON,JLEV)-ZQHM(JLON,JLEV)))
    ENDDO
  ENDDO

  ! store cloud water content for RTTOV
  IF (YIRAD%LGP) YDVARS%IRAD%T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) = ZQICE(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
  IF (YLRAD%LGP) YDVARS%LRAD%T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) = ZQLIQ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)

  ! Hannu Savijarvi diffuse -> direct albedo correction from hlradia,
  ! Assuming that SURFEX does not make difference between  
  ! dir/dif albedo as surfex/SURFEX/albedo_from_nir_vis.F90 defines
  ! PSCA_ALB(:,:) = PDIR_ALB(:,:)
  
! Albedo dans les intervalles, direct (parallel) et diffus (diffuse).
  IF (NSW==6.OR.NSW==1) THEN
    IF (LLMSE) THEN
      DO JSW=1,NSW
        ZALBP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)=YDCPG_GPAR%ALBDIR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)
        ZALBD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)=YDCPG_GPAR%ALBSCA(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)
        IF (LHLRADUPD) THEN
          DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            ZSALBCOR=0.2_JPRB/(1._JPRB+YDMF_PHYS_TMP%RDG%MU0(JLON))-0.12_JPRB
            ZALBP(JLON,JSW)=ZALBD(JLON,JSW)+ZSALBCOR
          ENDDO
        ENDIF
      ENDDO
    ELSEIF (LSFORCS) THEN
      DO JSW=1,NSW
        ZALBP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)=RALB_FORC
        ZALBD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)=RALB_FORC
!  direct>diffuse correction might be applied to RALB_FORC,too:
!              ZALBP(JLON,JSW)=RALB_FORC+ZSALBCOR
      ENDDO
    ELSE
     !pour pouvoir tourner sans la surface
      DO JSW=1,NSW
        ZALBP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)=YDMF_PHYS_SURF%GSD_VF%PALBF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
        ZALBD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)=YDMF_PHYS_SURF%GSD_VF%PALBF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
!              ZALBP(JLON,JSW)=PALBIN(JLON)+ZSALBCOR
      ENDDO
    ENDIF

  ! Spectral average albedo done with RSUN2 weights, 
  ! to be applied for HLRADIA, ACRANEB2 which use a single solar spectral band
    IF (LHLRADUPD) THEN
      ZALBP1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0._JPRB
      ZALBD1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0._JPRB
      DO JSW=1,NSW
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZALBP1(JLON)=ZALBP1(JLON)+RSUN2(JSW)*ZALBP(JLON,JSW)
          ZALBD1(JLON)=ZALBD1(JLON)+RSUN2(JSW)*ZALBD(JLON,JSW)
        ENDDO
      ENDDO
    ELSE
       ZALBP1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_GPAR%ALBDIR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
       ZALBD1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_GPAR%ALBSCA(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
    ENDIF

  ELSE

    CALL ABOR1 ('ALBEDO FOR NSW/= 1 or 6 not defined in apl_arome')

  ENDIF

  ! all albedo operations

  IF (LLMSE) THEN
    ZEMIS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_GPAR%VEMIS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
    ZTSURF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_GPAR%VTS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
    ! protection for E Zone, Where surface scheme send back EMIS and T =0
    ! the protection in aro_ground_paramn is not sufficient !!! WHY ??
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      IF (ZEMIS(JLON)==0._JPRB) THEN
        ZEMIS(JLON)=1.0_JPRB
        ZTSURF(JLON)=288.0_JPRB
      ENDIF
    ENDDO
  ELSEIF (LSFORCS) THEN
    ZEMIS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=REMIS_FORC
  ELSE
    ZEMIS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0.5_JPRB ! value 0.5 is suspicious
    ZTSURF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=ZTM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,YDCPG_DIM%KFLEVG)
  ENDIF !LLMSE EMIS
  
  ! initialization of aerosols
  IF (LRAY.AND.NRAY == 2.AND.LRADFLEX.AND.LRSTAER) THEN

    ! old ("standard") aerosols for LRAY only
    ZAER(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KTDIA-1,1)=0._JPRB
    DO JLEV=YDCPG_DIM%KTDIA-1,YDCPG_DIM%KFLEVG
      ZVETAH(JLEV)=STPREH(JLEV)/VP00
    ENDDO
    ZAEO=AERCS1*ZVETAH(YDCPG_DIM%KTDIA-1)+AERCS3*ZVETAH(YDCPG_DIM%KTDIA-1)**3+AERCS5*ZVETAH(YDCPG_DIM%KTDIA-1)**5
    DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
      ZAEN=AERCS1*ZVETAH(JLEV)+AERCS3*ZVETAH(JLEV)**3+AERCS5*ZVETAH(JLEV)**5
      ZAER(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV,1)=ZAEN-ZAEO
      ZAEO=ZAEN
    ENDDO
    ZAER(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,2:6)=0._JPRB
  
  ELSE
    
    IF (NAER >= 1 ) THEN
      IF(YSD_VAD%NUMFLDS >= 4) THEN
        CALL RADAER ( YDMODEL%YRML_PHY_RAD%YREAERD, YDERAD, YDPHY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA,       &
        & YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YLMF_PHYS_STATE%YCPG_PHY%PRE, YLMF_PHYS_STATE%YCPG_PHY%PREF,   &
        & ZTM, ZTSURF, YDMF_PHYS_SURF%GSD_VA%PSEA, YDMF_PHYS_SURF%GSD_VA%PLAN, YDMF_PHYS_SURF%GSD_VA%PSOO, &
        & YDMF_PHYS_SURF%GSD_VA%PDES, YDMF_PHYS_SURF%GSD_VA%PSUL, YDMF_PHYS_SURF%GSD_VA%PVOL, ZAER,        &
        & ZAERINDS)
      ELSE
        WRITE(NULOUT,*) 'YSD_VAD%NUMFLDS SHOULD BE >= 4, IT IS: ',YSD_VAD%NUMFLDS
        CALL ABOR1('APL_AROME: PB AEROSOLS!')
        ! NB : this abort excludes the use of radact. REK.
      ENDIF
    ENDIF

    IF (LRDUST) THEN
      ! We use the extinction coefficient explicitly solved by ARO_MNHDUST
      ZAER(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,3) = ZAERD_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
    ENDIF

  ENDIF
  ! end of old or new aerosols

  ! initialization of ozone
  IF (NOZOCL == 1) THEN
    ! as in IFS
    CALL RADOZC(YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, 1, YDCPG_DIM%KLON, &
    & 0, YLMF_PHYS_STATE%YCPG_PHY%PRE, YDGSGEOM%GEMU, ZROZ)
    DO JK=1,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZQO3(JLON,JK)=ZROZ(JLON,JK)/YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP(JLON,JK)
      ENDDO
    ENDDO
  ELSEIF (NOZOCL == 2) THEN
    ! as in ARPEGE (from clim profiles)
    CALL SUOZON(YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, ZQO3, .FALSE., YLMF_PHYS_STATE%YCPG_PHY%PRE, &
    & YLMF_PHYS_STATE%YCPG_PHY%XYB%RDELP, LO3ABC, YDMF_PHYS_SURF%GSD_VC%PGROUP)
  ENDIF

ELSE

  DO JSW=1,NSW
    ZALBP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)=0._JPRB 
    ZALBD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)=0._JPRB
  ENDDO

ENDIF
 !of preparation of input for LRAYFM, LRAY at every time step
 
 IF (LRAYFM) THEN
   ! Intermittent call to radiation interface
   IF (MOD(YDCPG_DIM%KSTEP,NRADFR) == 0) THEN 
     CALL RECMWF (YDGEOMETRY%YRDIMV, YDMODEL, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,   &
     & YDCPG_DIM%KSW, ZALBD, ZALBP, YLMF_PHYS_STATE%YCPG_PHY%PRE, YLMF_PHYS_STATE%YCPG_PHY%PREF, YDVARS%A%T1, ZQO3, &
     & ZAER, YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP, ZEMIS, YDMF_PHYS_TMP%RDG%MU0M, ZQV, ZQSAT, ZQICE,                   &
     & ZQLIQ, ZQSM, ZQRM, YDMF_PHYS_SURF%GSD_VF%PLSM, ZTM, ZTSURF, ZGP2DSPP, ZP1EZDIAG, YDMF_PHYS%RAD%EMTD,         &
     & YDMF_PHYS%RAD%EMTU, YDMF_PHYS%RAD%TRSW, YDMF_PHYS%OUT%FRTHC, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%FRSOC,        &
     & YDMF_PHYS%OUT%FRSO, ZZS_FSWDIR, ZZS_FSWDIF, ZFSDNN, ZFSDNV, ZCTRSO, ZCEMTR, ZTRSOD, ZTRSODIR,                &
     & ZTRSODIF, ZPIZA_DST_, ZCGA_DST_, ZTAUREL_DST_, ZAERINDS, YDGSGEOM%GELAM, YDGSGEOM%GEMU, YDCPG_GPAR%SWDIR,    &
     & YDCPG_GPAR%SWDIF, YDMF_PHYS_TMP%RDG%MU0LU, ZALBD1, ZFRSOLU)
   ELSE
     IF (LLMSE) THEN 
       DO JSW=1,NSW 
         ZTRSODIR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)=YDCPG_GPAR%SWDIR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)
         ZTRSODIF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)=YDCPG_GPAR%SWDIF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)
       ENDDO 
     ENDIF
     ZCTRSO(:,:)=0._JPRB
   ENDIF
   ! daand: radflex
   IF (LRADFLEX) THEN
     YLRADPROC => NEWINTPROC(YLPROCSET,'Radiation')
     ZFRSO => NEWINTFIELD(YLRADPROC,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,'FRSO','H','F')
     ZFRTH => NEWINTFIELD(YLRADPROC,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,'FRTH','H','F')
   ENDIF

    DO JLEV=1,YDCPG_DIM%KFLEVG
      ZTENT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=0.0_JPRB
    ENDDO

   ZSUDU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0.0_JPRB

   CALL RADHEAT  (  YDERAD, YDERDI, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,               &
   & YDCPG_DIM%KFLEVG, YLMF_PHYS_STATE%YCPG_PHY%PRE, ZEMIS, YDMF_PHYS%RAD%EMTD, YDMF_PHYS_TMP%RDG%MU0,                   &
   & ZQVM, ZTENT, YDMF_PHYS%RAD%TRSW, ZTRSOD, ZTSURF, PDTPHY, ZTRSODIR, ZTRSODIF, ZALBD, ZALBP,                          &
   & ZFRSO, ZFRTH, YDMF_PHYS%OUT%FRSODS, YDMF_PHYS%OUT%FRTHDS, ZCEMTR, ZCTRSO, YDMF_PHYS%OUT%FRSOC, YDMF_PHYS%OUT%FRTHC, &
   & ZSUDU, ZSDUR, ZDSRP, ZZS_FSWDIR, ZZS_FSWDIF, YDMF_PHYS%OUT%FRSOPS, ZFRSOFS, YDMF_PHYS%OUT%FRSOPT                    &
   &  )

  ! daand: radflex
  IF (LRADFLEX) THEN
    ! store for further calculations and diagnostics
    ! warning : pointers. REK
    YDMF_PHYS%OUT%FRSO(:,:,1)=ZFRSO
    YDMF_PHYS%OUT%FRTH(:,:,1)=ZFRTH
  ELSE
  ! daand: if LRADFLEX, the contribution to temperature is done by
  ! cptend_flex/cputqy
    ! update temperature tendency by radiative contribution
    DO JLEV=1,YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
        ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV)+ZTENT(JLON,JLEV)
      ENDDO
    ENDDO
  ENDIF

  DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
    ! update sunshine duration [s]
    YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)=YDMF_PHYS_SURF%GSD_VD%PSUND(JLON)+ZSDUR(JLON)*TSTEP
    ! Estimate of the direct normal irradiance, with securities
    IF (YDMF_PHYS_TMP%RDG%MU0(JLON) > 3.0E-02_JPRB) THEN
      YDMF_PHYS%OUT%FRSDNI(JLON)=MAX(0.0_JPRB,YDMF_PHYS%OUT%FRSOPS(JLON)/YDMF_PHYS_TMP%RDG%MU0(JLON))
    ELSE
      YDMF_PHYS%OUT%FRSDNI(JLON)=MAX(0.0_JPRB,YDMF_PHYS%OUT%FRSOPS(JLON))
    ENDIF
  ENDDO 

  IF( MOD(YDCPG_DIM%KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'sous apl_arome apres rayonnement ZTENT=',ZTENT(NPTP,30:41)
    IF (LLMSE) THEN
      DO JSW=1, NSW
        WRITE(NULOUT,*)'ZSFSWDIR ZSFSWDIF ZFSDNN ZFSDNV PFRSO',&
         & ZZS_FSWDIR(NPTP,JSW),ZZS_FSWDIF(NPTP,JSW),ZFSDNN(NPTP), ZFSDNV(NPTP),YDMF_PHYS%OUT%FRSO(NPTP,YDCPG_DIM%KFLEVG,1)
        WRITE(NULOUT,*)'ZALBD ZALBP',ZALBD(NPTP,JSW),ZALBP(NPTP,JSW)
      ENDDO
    ENDIF
    WRITE(NULOUT,*)ZFSDNN(NPTP),ZFSDNV(NPTP)
    WRITE (NULOUT,*)'TSURF EMIS ZFRTH',ZTSURF(NPTP),ZEMIS(NPTP),YDMF_PHYS%OUT%FRTHDS(NPTP)
  ENDIF

  IF (LFLEXDIA) THEN
    IF (LDDH_OMP) THEN
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, YDMF_PHYS%OUT%FRSO(:, :, 1), 'FCTRAYSO', YDDDH&
      & )
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, YDMF_PHYS%OUT%FRTH(:, :, 1), 'FCTRAYTH', YDDDH&
      & )
    ELSE
      CALL ADD_FIELD_3D(YLDDH, YDMF_PHYS%OUT%FRSO(:, :, 1), 'FCTRAYSO', 'F', 'ARP', .TRUE., .TRUE.)
      CALL ADD_FIELD_3D(YLDDH, YDMF_PHYS%OUT%FRTH(:, :, 1), 'FCTRAYTH', 'F', 'ARP', .TRUE., .TRUE.)
    ENDIF
  ENDIF

ELSE

  YDMF_PHYS%OUT%FRSOC(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,0:1)=0.0_JPRB
  YDMF_PHYS%OUT%FRTHC(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,0:1)=0.0_JPRB

ENDIF  ! LRAYFM

!     ------------------------------------------------------------------
!     NEBULOSITE (CONVECTIVE+STRATIFORME) A TROIS NIVEAUX.
!     DIAGNOSTIC OF THREE LEVELS (CONVECTIVE+STRATIFORM) CLOUDINESS.

! protect cloudiness from being 0 or 1 (needed for ACRANEB2 and ACNPART)
DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZNEB0(JLON,JLEV)=MAX(ZEPSNEB,MIN(1._JPRB-ZEPSNEB,YDVARS%A%T1(JLON,JLEV)))
  ENDDO
ENDDO

! decorrelation depth for cloud overlaps

IF (LRNUEXP) THEN
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZDECRD(JLON)=RDECRD1+RDECRD2*EXP(-((ASIN(YDGSGEOM%GEMU(JLON))-RDECRD3*RDECLI)/RDECRD4)**2)
  ENDDO
ENDIF

! calculate high, medium, low and total cloud cover
CALL ACNPART(YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, NTNEBU, YDCPG_DIM%KFLEVG, &
& YLMF_PHYS_STATE%YCPG_DYN%PHI, YLMF_PHYS_STATE%YCPG_DYN%PHIF, YLMF_PHYS_STATE%YCPG_PHY%PREF, ZDECRD,         &
& ZNEB0, YDMF_PHYS%OUT%CLCH, YDMF_PHYS%OUT%CLCM, YDMF_PHYS%OUT%CLCL, YDCPG_MISC%CLCT, ZCLCT_RAD)

IF (LRAY.AND.NRAY == 2.AND.LRADFLEX) THEN

  ! -------------------------
  ! ACRANEB2 radiation scheme
  ! -------------------------

!+++ The next input preparations are redundant:

  ! initialization of cloud ice, cloud liquid and specific humidity
  ! (with respect to moist air, i.e. excluding hydrometeors)
  DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZQICE(JLON,JLEV)=MAX(0.0_JPRB, ZQIM(JLON,JLEV)/&
       & (1.0_JPRB-ZQIM(JLON,JLEV)-ZQCM(JLON,JLEV)-ZQRM(JLON,JLEV)&
       & -ZQGM(JLON,JLEV)-ZQSM(JLON,JLEV)-ZQHM(JLON,JLEV)))
      ZQLIQ(JLON,JLEV)=MAX(0.0_JPRB, ZQCM(JLON,JLEV)/&
       & (1.0_JPRB-ZQIM(JLON,JLEV)-ZQCM(JLON,JLEV)-ZQRM(JLON,JLEV)&
       & -ZQGM(JLON,JLEV)-ZQSM(JLON,JLEV)-ZQHM(JLON,JLEV)))
      ZQV(JLON,JLEV)=MAX(0.0_JPRB, ZQVM(JLON,JLEV)/&
       & (1.0_JPRB-ZQIM(JLON,JLEV)-ZQCM(JLON,JLEV)-ZQRM(JLON,JLEV)&
       & -ZQGM(JLON,JLEV)-ZQSM(JLON,JLEV)-ZQHM(JLON,JLEV)))
    ENDDO
  ENDDO

  ! store cloud water content for RTTOV
  IF (YIRAD%LGP) YDVARS%IRAD%T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) = ZQICE(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)
  IF (YLRAD%LGP) YDVARS%LRAD%T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:) = ZQLIQ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)

  ! initialization of ozone
  IF (NOZOCL == 1) THEN
    ! as in IFS
    CALL RADOZC(YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, 1, YDCPG_DIM%KLON, &
    & 0, YLMF_PHYS_STATE%YCPG_PHY%PRE, YDGSGEOM%GEMU, ZROZ)
    DO JK=1,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZQO3(JLON,JK)=ZROZ(JLON,JK)/YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP(JLON,JK)
      ENDDO
    ENDDO
  ELSEIF (NOZOCL == 2) THEN
    ! as in ARPEGE (from clim profiles)
    CALL SUOZON(YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, ZQO3, .FALSE., YLMF_PHYS_STATE%YCPG_PHY%PRE, &
    & YLMF_PHYS_STATE%YCPG_PHY%XYB%RDELP, LO3ABC, YDMF_PHYS_SURF%GSD_VC%PGROUP)
  ENDIF

  ! initialization of aerosols
  IF (LRAY.AND.NRAY == 2.AND.LRADFLEX.AND.LRSTAER) THEN

    ZAER(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KTDIA-1,1)=0._JPRB
    ! old ("standard") aerosols
    DO JLEV=YDCPG_DIM%KTDIA-1,YDCPG_DIM%KFLEVG
      ZVETAH(JLEV)=STPREH(JLEV)/VP00
    ENDDO
    ZAEO=AERCS1*ZVETAH(YDCPG_DIM%KTDIA-1)+AERCS3*ZVETAH(YDCPG_DIM%KTDIA-1)**3+AERCS5*ZVETAH(YDCPG_DIM%KTDIA-1)**5
    DO JLEV=YDCPG_DIM%KTDIA,YDCPG_DIM%KFLEVG
      ZAEN=AERCS1*ZVETAH(JLEV)+AERCS3*ZVETAH(JLEV)**3+AERCS5*ZVETAH(JLEV)**5
      ZAER(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV,1)=ZAEN-ZAEO
      ZAEO=ZAEN
    ENDDO
    ZAER(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,2:6)=0._JPRB
  
  ELSE

    IF (NAER >= 1) THEN
      IF (YSD_VAD%NUMFLDS >= 4) THEN
        ! initialisation of aerosols as in ARPEGE (from clim files)
        CALL RADAER (YDMODEL%YRML_PHY_RAD%YREAERD, YDERAD, YDPHY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,       &
        & YDCPG_DIM%KFLEVG, YLMF_PHYS_STATE%YCPG_PHY%PRE, YLMF_PHYS_STATE%YCPG_PHY%PREF, ZTM, ZTSURF,                     &
        & YDMF_PHYS_SURF%GSD_VA%PSEA, YDMF_PHYS_SURF%GSD_VA%PLAN, YDMF_PHYS_SURF%GSD_VA%PSOO, YDMF_PHYS_SURF%GSD_VA%PDES, &
        & YDMF_PHYS_SURF%GSD_VA%PSUL, YDMF_PHYS_SURF%GSD_VA%PVOL, ZAER, ZAERINDS)
      ELSE
        CALL ABOR1('APL_AROME: PB AEROSOLS!')
        ! NB : this abort excludes the use of radact. REK.
      ENDIF
    ENDIF

    IF (LRDUST) THEN
      ! We uses the extinction coefficient explicitely solved by ARO_MNHDUST
      ZAER(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,3) = ZAERD_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
    ENDIF

  ENDIF ! (LRAY.AND.NRAY == 2.AND.LRADFLEX.AND.LRSTAER) THEN

  ! get diffuse and direct surface albedo, emissivity and temperature
  IF (.NOT.LHLRADUPD) THEN
    ZALBD1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_GPAR%ALBSCA(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
    ZALBP1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_GPAR%ALBDIR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)
  ENDIF
  ZEMIS  (YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_GPAR%VEMIS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
  ZTSURF (YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDCPG_GPAR%VTS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ! protection of E-zone (not to have zero emissivity and T_surf there)
    IF (ZEMIS(JLON) == 0._JPRB) THEN
      ZEMIS (JLON)=  1._JPRB
      ZTSURF(JLON)=288._JPRB
    ENDIF
  ENDDO

!+++ End of redundant input preparations for ACRANEB

  ! initialization of CO2(+), differs from IFS radiation scheme!
  ZQCO2(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=QCO2

  ! daand: radflex
  YLRADPROC => NEWINTPROC(YLPROCSET,'Radiation')
  ZFRSO => NEWINTFIELD(YLRADPROC,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG, 'FRSO','H','F')
  ZFRTH => NEWINTFIELD(YLRADPROC,YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG, 'FRTH','H','F')

  ! call radiation scheme
  IJN=YDCPG_DIM%KLON
  CALL ACRANEB2(YDERDI, YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,                         &
  & NTRADI, YDCPG_DIM%KFLEVG, IJN, YDCPG_DIM%KSTEP, YDCFU%NFRRC, YLMF_PHYS_STATE%YCPG_PHY%PRE, YLMF_PHYS_STATE%YCPG_PHY%PREF, &
  & YLMF_PHYS_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_STATE%YCPG_DYN%RCP%R, YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP,                       &
  & ZNEB0, ZQV, ZQCO2, ZQICE, ZQLIQ, ZQO3, YLMF_PHYS_STATE%T, ZALBD1, ZALBP1, ZEMIS, YDGSGEOM%GELAM,                          &
  & YDGSGEOM%GEMU, YDMF_PHYS_TMP%RDG%MU0, YDMF_PHYS_TMP%RDG%MU0LU, ZTSURF, ZDECRD, ZCLCT_RAD, YDMF_PHYS%OPT%GDEOSI,           &
  & YDMF_PHYS%OPT%GUEOSI, YDMF_PHYS%OPT%GMU0, YDMF_PHYS%OPT%GMU0_MIN, YDMF_PHYS%OPT%GMU0_MAX, YDMF_PHYS%OPT%GDEOTI,           &
  & YDMF_PHYS%OPT%GDEOTI2, YDMF_PHYS%OPT%GUEOTI, YDMF_PHYS%OPT%GUEOTI2, YDMF_PHYS%OPT%GEOLT, YDMF_PHYS%OPT%GEOXT,             &
  & YDMF_PHYS%OPT%GRPROX, YDMF_PHYS%OPT%GMIXP, YDMF_PHYS%OPT%GFLUXC, YDMF_PHYS%OPT%GRSURF, YDMF_PHYS_SURF%GSD_VD%PSUND,       &
  & ZFRSO, ZFRTH, YDMF_PHYS%OUT%FRSOC, YDMF_PHYS%OUT%FRTHC, ZFRSODS, YDMF_PHYS%OUT%FRSOPS, ZFRSOLU, YDMF_PHYS%OUT%FRTHDS,     &
  & ZAER)

  ! daand: radflex
  ! store for further calculations and diagnostics
  ! warning : pointers. REK
  YDMF_PHYS%OUT%FRSO(:,:,1)=ZFRSO
  YDMF_PHYS%OUT%FRTH(:,:,1)=ZFRTH

  ! extract surface fluxes
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    YDMF_PHYS%OUT%FRSODS(JLON)=ZFRSODS(JLON)+YDMF_PHYS%OUT%FRSOPS(JLON)  ! downward surface sw flux
  ENDDO

  IF (LLMSE) THEN
    IF (LHLRADUPD) THEN
      DO JSW = 1,NSW
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZZS_FSWDIR(JLON,JSW) = YDMF_PHYS%OUT%FRSOPS(JLON)*RSUN2(JSW)
          ZZS_FSWDIF(JLON,JSW) = ZFRSODS(JLON)*RSUN2(JSW)
         ENDDO
      ENDDO
    ELSE
      ZZS_FSWDIR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)=YDMF_PHYS%OUT%FRSOPS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA) ! direct surface swdn flux
      ZZS_FSWDIF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1)=ZFRSODS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA) ! diffuse surface swdn flux
    ENDIF
  ENDIF

  ! Estimate of the direct normal irradiance, with securities
  YDMF_PHYS%OUT%FRSDNI(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=YDMF_PHYS%OUT%FRSOPS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
  DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
    IF (YDMF_PHYS_TMP%RDG%MU0(JLON) > 3.0E-02_JPRB) THEN
      YDMF_PHYS%OUT%FRSDNI(JLON)=YDMF_PHYS%OUT%FRSOPS(JLON)/YDMF_PHYS_TMP%RDG%MU0(JLON)
    ENDIF
  ENDDO
  DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
    YDMF_PHYS%OUT%FRSDNI(JLON)=MAX(0.0_JPRB,YDMF_PHYS%OUT%FRSDNI(JLON))
  ENDDO

  IF (LFLEXDIA) THEN
    IF (LDDH_OMP) THEN
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, YDMF_PHYS%OUT%FRSO(:, :, 1), 'FCTRAYSO', YDDDH&
      & )
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, YDMF_PHYS%OUT%FRTH(:, :, 1), 'FCTRAYSO', YDDDH&
      & )
    ELSE
      CALL ADD_FIELD_3D(YLDDH, YDMF_PHYS%OUT%FRSO(:, :, 1), 'FCTRAYSO', 'F', 'ARP', .TRUE., .TRUE.)
      CALL ADD_FIELD_3D(YLDDH, YDMF_PHYS%OUT%FRTH(:, :, 1), 'FCTRAYTH', 'F', 'ARP', .TRUE., .TRUE.)
    ENDIF
  ENDIF

ENDIF

IF (.NOT.(LRAY.AND.NRAY == 2.AND.LRADFLEX).AND..NOT.LRAYFM) THEN
  DO JSW = 1,NSW
    ZZS_FSWDIR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW) = 0._JPRB
    ZZS_FSWDIF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW) = 0._JPRB
  ENDDO
  YDMF_PHYS%OUT%FRSOPS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0._JPRB
ENDIF

IF (LFLEXDIA) THEN
  CALL ARO_STARTBU( YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, NRR, NGFL_EXT, ZRHODJM__(:, 1:YDCPG_DIM%KFLEVG),     &
  & ZUS__(:, 1:YDCPG_DIM%KFLEVG), ZVS__(:, 1:YDCPG_DIM%KFLEVG), ZWS__(:, 1:YDCPG_DIM%KFLEVG), ZTHS__(:, 1:YDCPG_DIM%KFLEVG), &
  & ZRS_, ZTKES_, YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH)
ENDIF


!    ------------------------------------------------------------------
!     7 - CONVECTION. 
!     --------------------------------------------------------------------

IF(LKFBCONV) THEN

  ! No swapp needed becaus IN and OUT are not needed simultaneously

  CALL BRI2ACCONV(YDMODEL%YRML_PHY_MF, YDGEOMETRY%YREGEO, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFDIA,                          &
  & YDCPG_DIM%KFLEVG, YDGSGEOM%GM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZPABSM__(:, 1:YDCPG_DIM%KFLEVG),                                  &
  & ZZZ_F_, ZTM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, :), ZRM_(:, :, 1), ZRM_(:, :, 2), ZRM_(:, :, 4), ZRHODREFM__(:, 1:YDCPG_DIM%KFLEVG), &
  & ZUM__(:, 1:YDCPG_DIM%KFLEVG), ZVM__(:, 1:YDCPG_DIM%KFLEVG), ZWM__(:, 1:YDCPG_DIM%KFLEVG), ZMFS_,                                  &
  & ZCVTENDT_, ZCVTENDRV_, ZCVTENDRC_, ZCVTENDRI_, ZCVTENDPR_, ZCVTENDPRS_   )

  IF(MOD(YDCPG_DIM%KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)"Pluie conv au sol", ZCVTENDPR_(NPTP), &
     & MAXVAL(ZCVTENDPR_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)) ,MINVAL(ZCVTENDPR_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA))
  ENDIF

  DO JLEV = 1,YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
      ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV) + ZCVTENDT_(JLON,JLEV)
      ZTENDR(JLON,JLEV,1) = ZTENDR(JLON,JLEV,1) + ZCVTENDRV_(JLON,JLEV)*ZQDM(JLON,JLEV)  
      ZTENDR(JLON,JLEV,2) = ZTENDR(JLON,JLEV,2) + ZCVTENDRC_(JLON,JLEV)*ZQDM(JLON,JLEV)  
      ZTENDR(JLON,JLEV,4) = ZTENDR(JLON,JLEV,4) + ZCVTENDRI_(JLON,JLEV)*ZQDM(JLON,JLEV)  
      ZRS_(JLON,JLEV,1)=ZRS_(JLON,JLEV,1)+ZCVTENDRV_(JLON,JLEV)
      ZRS_(JLON,JLEV,2)=ZRS_(JLON,JLEV,2)+ZCVTENDRC_(JLON,JLEV)
      ZRS_(JLON,JLEV,4)=ZRS_(JLON,JLEV,4)+ZCVTENDRI_(JLON,JLEV)
      ZTHS__(JLON,JLEV)=ZTHS__(JLON,JLEV)+ZCVTENDT_(JLON,JLEV)*(RATM/YLMF_PHYS_STATE%YCPG_PHY%PREF(JLON,JLEV))**(RD/RCPD)  
    ENDDO
  ENDDO
  DO JLON =YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
    ZINPRR_(JLON)=ZINPRR_(JLON)+ZCVTENDPR_(JLON)-ZCVTENDPRS_(JLON)
    ZACPRR_(JLON)=ZACPRR_(JLON)+(ZCVTENDPR_(JLON)-ZCVTENDPRS_(JLON))*PDTPHY
    ZINPRS_(JLON)=ZINPRS_(JLON)+ZCVTENDPRS_(JLON)
    ZACPRS_(JLON)=ZACPRS_(JLON)+ZCVTENDPRS_(JLON)*PDTPHY
  ENDDO
  ! avance temporelle et inversion niveau pour ZMFS_
  ! on utilise PSIGS pour le flux de masse pour la condensation sous maille 
  ! car PSIGS n est utilise que si LOSIGMAS=T
  IF (LOSUBG_COND.AND..NOT.LOSIGMAS) THEN
    YDVARS%SRC%T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=ZMFS_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  ENDIF
  IF(MOD(YDCPG_DIM%KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)"aps CONV, TENRV, TENRC, TENRI"
    DO JLEV=1,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,*)ZTENDR(NPTP,JLEV,1),ZTENDR(NPTP,JLEV,2),ZTENDR(NPTP,JLEV,4)
    ENDDO
  ENDIF
  CALL ARO_CONVBU(YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, NRR, ZRHODJM__(:, 1:YDCPG_DIM%KFLEVG), ZRS_, ZTHS__(:, 1:YDCPG_DIM%KFLEVG), &
  & YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH)

ENDIF

!    ------------------------------------------------------------------
!     8 - SURFACE. 
!     --------------------------------------------------------------------

IF (LLMSE) THEN
! A loop around SURFEX in order to test OpenMP

  SURFEX_LOOP : DO ISURFEX = 1, NSURFEX_ITER

! Initialisations 

  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZZS_(JLON)=YDOROG%OROG(JLON)*ZINVG 
  ENDDO
  DO JLEV = 1,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZDEPTH_HEIGHT_(JLON,JLEV)=ZZZ_F_(JLON,JLEV)-ZZS_(JLON)  
    ENDDO
  ENDDO
  IF (MINVAL(ZDEPTH_HEIGHT_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,IKB)) <= 0._JPRB) THEN
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      IF (ZDEPTH_HEIGHT_(JLON,IKB) <= 0._JPRB) THEN
        WRITE (NULOUT,*)'sous apl_arome pb height en', JLON,ZAPHIFM(JLON,YDCPG_DIM%KFLEVG),YDOROG%OROG(JLON)
      ENDIF
    ENDDO
  ENDIF
  ! Can't use a section of pointer. An explicit copy shows, by the way, that a copy is needed
  ! because data is not contiguous. REK
  ZSVMB_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:NGFL_EXT)=ZSVM_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,IKB,1:NGFL_EXT)

  IF (LSURFEX_CRITICAL) THEN

!$OMP CRITICAL (ARO_GROUND_PARAM_LOCK)

    IF (LLMSE_PARAM) THEN
      CALL ARO_GROUND_PARAM( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP, YDCPG_DIM%KFDIA, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA,                                                &
      & YDCPG_DIM%KSTEP, NRR, NSW, NGFL_EXT, NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM, LMPA,                                                                &
      & CCOUPLING, LLXFUMSE, NINDAT, ZRHGMT, ZSTATI, RSOVR, RCODEC, RSIDEC, YDCSGEOM%RINDX(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                     &
      & YDCSGEOM%RINDY(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZUM__(:, IKB), ZVM__(:, IKB), ZTM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG),                   &
      & ZRM_(:, IKB, 1), ZSVMB_, RCARDI, ZRHODREFM__(:, IKB), ZPABSM__(:, IKB), YLMF_PHYS_STATE%YCPG_PHY%PRE(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG), &
      & ZDTMSE, ZDEPTH_HEIGHT_(:, IKB), ZZS_, XZSEPS, YDMF_PHYS_TMP%RDG%MU0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                                    &
      & YDMF_PHYS_TMP%RDG%MU0N(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDGSGEOM%GELAM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                                &
      & YDGSGEOM%GEMU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), XSW_BANDS, ZINPRR_NOTINCR_, ZINPRS_NOTINCR_,                                                             &
      & ZINPRG_NOTINCR_, YDMF_PHYS%OUT%FRTHDS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZZS_FSWDIF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:NSW),                              &
      & ZZS_FSWDIR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:NSW), ZCFAQ_, ZCFATH_, ZCFAU_, ZCFBQ_, ZCFBTH_,                                                            &
      & ZCFBU_, ZCFBV_, ZSFTH_, ZSFRV_, ZSFSV_, ZSFCO2_, ZSFU_, ZSFV_, ZALBP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:NSW),                                            &
      & ZALBD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:NSW), ZEMIS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZTSURF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                          &
      & YDMF_PHYS%OUT%FRTH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, 1))

    ENDIF

    IF (LRCO2) THEN
      ZSFSV_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,NSV_CO2)= ZSFCO2_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
!print*,' FLUX CO2 =', MINVAL(ZSFSV_(KIDIA:KFDIA,NSV_CO2)),&
!                    & MAXVAL(ZSFSV_(KIDIA:KFDIA,NSV_CO2))
    ENDIF

!!!!! TEST DDH ATTENTION
!ZSFRV_(KIDIA:KFDIA) = 0._JPRB

    IF (LLMSE_DIAG) THEN

      CALL ARO_GROUND_DIAG( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP, YDCPG_DIM%KFDIA, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA,                     &
      & YDCPG_DIM%KFLEVG, IKL, NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM, ZZS_, ZSFRV_, ZUM__(:, IKTB:IKTE),                     &
      & ZVM__(:, IKTB:IKTE), ZDEPTH_HEIGHT_(:, IKTB:IKTE), YDMF_PHYS%OUT%FRTH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, 1), &
      & YDMF_PHYS%OUT%FRSO(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, 1), YDCSGEOM%RINDX(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),   &
      & YDCSGEOM%RINDY(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZQS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZGZ0_,                                &
      & ZGZ0H_, YDMF_PHYS%OUT%TCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%QCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),            &
      & YDMF_PHYS%OUT%RHCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%UCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                   &
      & YDMF_PHYS%OUT%VCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%NUCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                   &
      & YDMF_PHYS%OUT%NVCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%FCLL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1),                &
      & YDMF_PHYS%OUT%FCLN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1), YDMF_PHYS%OUT%FEVL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1),              &
      & YDMF_PHYS%OUT%FEVN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1), ZSSO_STDEV_, YDMF_PHYS_SURF%GSP_SG%PF_T1,                            &
      & ZBUDTH_, ZBUDSO_, ZFCLL_, ZTOWNS_, ZCD_                         )
      CALL ARO_GROUND_DIAG_2ISBA( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP, YDCPG_DIM%KFDIA, YDCPG_DIM%KIDIA,                                     &
      & YDCPG_DIM%KFDIA, NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM, YDCSGEOM%RINDX(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                  &
      & YDCSGEOM%RINDY(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS_SURF%GSD_VF%PLSM, ZDUMMY1,                                             &
      & ZDUMMY1, ZDUMMY1, ZTSURF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS_SURF%GSP_SG%PF_T1, ZDUMMY1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), &
      & ZDUMMY1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZDUMMY1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZDUMMY1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),     &
      & ZDUMMY1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZDUMMY1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZDUMMY1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),     &
      & YDMF_PHYS_SURF%GSP_SG%PR_T1, ZDUMMY1 )

    ENDIF
 
!$OMP END CRITICAL (ARO_GROUND_PARAM_LOCK)
  ELSE

    IF (LLMSE_PARAM) THEN

      CALL ARO_GROUND_PARAM( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP, YDCPG_DIM%KFDIA, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA,                                                &
      & YDCPG_DIM%KSTEP, NRR, NSW, NGFL_EXT, NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM, LMPA,                                                                &
      & CCOUPLING, LLXFUMSE, NINDAT, ZRHGMT, ZSTATI, RSOVR, RCODEC, RSIDEC, YDCSGEOM%RINDX(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                     &
      & YDCSGEOM%RINDY(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZUM__(:, IKB), ZVM__(:, IKB), ZTM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG),                   &
      & ZRM_(:, IKB, 1), ZSVMB_, RCARDI, ZRHODREFM__(:, IKB), ZPABSM__(:, IKB), YLMF_PHYS_STATE%YCPG_PHY%PRE(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG), &
      & ZDTMSE, ZDEPTH_HEIGHT_(:, IKB), ZZS_, XZSEPS, YDMF_PHYS_TMP%RDG%MU0(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                                    &
      & YDMF_PHYS_TMP%RDG%MU0N(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDGSGEOM%GELAM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                                                &
      & YDGSGEOM%GEMU(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), XSW_BANDS, ZINPRR_NOTINCR_, ZINPRS_NOTINCR_,                                                             &
      & ZINPRG_NOTINCR_, YDMF_PHYS%OUT%FRTHDS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZZS_FSWDIF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:NSW),                              &
      & ZZS_FSWDIR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:NSW), ZCFAQ_, ZCFATH_, ZCFAU_, ZCFBQ_, ZCFBTH_,                                                            &
      & ZCFBU_, ZCFBV_, ZSFTH_, ZSFRV_, ZSFSV_, ZSFCO2_, ZSFU_, ZSFV_, ZALBP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:NSW),                                            &
      & ZALBD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:NSW), ZEMIS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZTSURF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                          &
      & YDMF_PHYS%OUT%FRTH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, 1))

    ENDIF

    IF (LRCO2) THEN
      ZSFSV_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,NSV_CO2)= ZSFCO2_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
!print*,' FLUX CO2 =', MINVAL(ZSFSV_(KIDIA:KFDIA,NSV_CO2)),&
!                    & MAXVAL(ZSFSV_(KIDIA:KFDIA,NSV_CO2))
    ENDIF

!!!!! TEST DDH ATTENTION
!ZSFRV_(KIDIA:KFDIA) = 0._JPRB

    IF (LLMSE_DIAG) THEN

      CALL ARO_GROUND_DIAG( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP, YDCPG_DIM%KFDIA, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA,                     &
      & YDCPG_DIM%KFLEVG, IKL, NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM, ZZS_, ZSFRV_, ZUM__(:, IKTB:IKTE),                     &
      & ZVM__(:, IKTB:IKTE), ZDEPTH_HEIGHT_(:, IKTB:IKTE), YDMF_PHYS%OUT%FRTH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, 1), &
      & YDMF_PHYS%OUT%FRSO(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, 1), YDCSGEOM%RINDX(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),   &
      & YDCSGEOM%RINDY(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZQS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZGZ0_,                                &
      & ZGZ0H_, YDMF_PHYS%OUT%TCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%QCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),            &
      & YDMF_PHYS%OUT%RHCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%UCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                   &
      & YDMF_PHYS%OUT%VCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%NUCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                   &
      & YDMF_PHYS%OUT%NVCLS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS%OUT%FCLL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1),                &
      & YDMF_PHYS%OUT%FCLN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1), YDMF_PHYS%OUT%FEVL(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1),              &
      & YDMF_PHYS%OUT%FEVN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1), ZSSO_STDEV_, YDMF_PHYS_SURF%GSP_SG%PF_T1,                            &
      & ZBUDTH_, ZBUDSO_, ZFCLL_, ZTOWNS_, ZCD_                         )

      CALL ARO_GROUND_DIAG_2ISBA( YDCPG_DIM%KBL, YDCPG_DIM%KGPCOMP, YDCPG_DIM%KFDIA, YDCPG_DIM%KIDIA,                                     &
      & YDCPG_DIM%KFDIA, NDGUNG, NDGUXG, NDLUNG, NDLUXG, LSURFEX_KFROM, YDCSGEOM%RINDX(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                  &
      & YDCSGEOM%RINDY(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS_SURF%GSD_VF%PLSM, ZDUMMY1,                                             &
      & ZDUMMY1, ZDUMMY1, ZTSURF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDMF_PHYS_SURF%GSP_SG%PF_T1, ZDUMMY1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), &
      & ZDUMMY1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZDUMMY1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZDUMMY1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),     &
      & ZDUMMY1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZDUMMY1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZDUMMY1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),     &
      & YDMF_PHYS_SURF%GSP_SG%PR_T1, ZDUMMY1 )
 
    ENDIF

  ENDIF

  ENDDO SURFEX_LOOP

!* Compute PBL-diagnostics

   ZCAPE(:)=0._JPRB
   ZDCAPE(:)=0._JPRB   
   CALL ACCLDIA(YDXFU, YDPHY, YDMODEL%YRML_PHY_MF%YRPHY2, YDTOPH, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON,                               &
   & YDCPG_DIM%KFLEVG, YDMF_PHYS%OUT%UCLS, YDMF_PHYS%OUT%VCLS, YLMF_PHYS_STATE%U(:, 1:YDCPG_DIM%KFLEVG),                                          &
   & YLMF_PHYS_STATE%V(:, 1:YDCPG_DIM%KFLEVG), ZCAPE, ZDCAPE, ZTKEM(:, 1:YDCPG_DIM%KFLEVG), YLMF_PHYS_STATE%YCPG_DYN%PHIF(:, 1:YDCPG_DIM%KFLEVG), &
   & YDOROG%OROG, YDMF_PHYS%OUT%UGST, YDMF_PHYS%OUT%VGST, YDMF_PHYS%OUT%CLPH, ICLPH)

   YDMF_PHYS%OUT%CLPH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=MIN(XMAXLM,MAX(XMINLM,YDMF_PHYS%OUT%CLPH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA))) 

   CALL ACVISIH(YDVISI, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KTDIA, YDCPG_DIM%KFLEVG, &
   & YLMF_PHYS_STATE%YCPG_DYN%PHI, YLMF_PHYS_STATE%YCPG_DYN%PHIF, YLMF_PHYS_STATE%YCPG_PHY%PREF, ZTM,        &
   & ZRHM, ZQCM, ZQIM, ZQRM, ZQSM, ZQGM, YDMF_PHYS%OUT%VISICLD, YDMF_PHYS%OUT%VISIHYD, YDMF_PHYS%OUT%MXCLWC  &
   & )

ELSE

  ZSFSV_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,:)=0._JPRB

ENDIF    !  <== End block "IF (LMSE)"

!*            IDEALIZED TURBULENT SURFACE FLUXES FOR SQUALL LINE CASE
!                --------------------------------------------------------

IF (LSQUALL.AND.LTURB) THEN
  ! on n'a besoin que d'un flux sur V (U est nul). 
  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    IF (ABS(ZVM__(JLON,IKB)) <= 1.E-12) THEN
      ZSFV_(JLON)=0._JPRB
    ELSE
      ZSFV_(JLON)=-(ZVM__(JLON,IKB))**2 *&
       & (0.4_JPRB  /(LOG(ZZZ_F_(JLON,IKB)/0.2_JPRB) ) )**2&
       & *ZVM__(JLON,IKB)/ABS(ZVM__(JLON,IKB))  
    ENDIF
  ENDDO
ENDIF

!    ------------------------------------------------------------------
!    9.  Shallow Mass Flux Mixing
!    ------------------------------------------------------------------


IF (LMFSHAL) THEN
  IF (CMF_UPDRAFT=='DUAL') THEN
    ! Updraft computation from EDMF/ECMWF dual proposal
    ! Version May 2007
    !
    ! The following routine  are using arrays with the vertical Arpege/IFS fashion (as in the radiation scheme)

    IDRAFT = 2 ! beginning of the loop for MF tendency equation
               ! only 2 and 3 are used for tendency computation in ARO_SHALLOW_MF
    INDRAFT=3   ! 1 for test, 2 for dry, 3 for wet

    IF (YDCPG_DIM%KMAXDRAFT < INDRAFT) THEN
      CALL ABOR1('APL_AROME : KMAXDRAFT TOO SMALL !')
    ENDIF

    DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
      ZZS_FTH_(JLON)=-1._JPRB*ZSFTH_(JLON)*(YLMF_PHYS_STATE%YCPG_PHY%PRE(JLON,YDCPG_DIM%KFLEVG)*ZINVATM)**(ZRSCP)
      ZZS_FRV_(JLON)=-1._JPRB*ZSFRV_(JLON)
    ENDDO
    ZZS_FU_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=ZSFU_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
    ZZS_FV_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=ZSFV_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)

    !  IF LHARATU=TRUE then TKE at t-dt is needed as input for vdfexcuhl so fill ZTKEEDMF with t-1 value  from PTKEM

    IF (LHARATU) THEN
      DO JLEV=1,YDCPG_DIM%KFLEVG
        ZTKEEDMF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=ZTKEM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)
        ZLENGTH_M(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=0.01_JPRB
        ZLENGTH_H(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=0.01_JPRB
      ENDDO
      IF (MAXVAL(ZTKEM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)) > 3300._JPRB) THEN
        DO JLEV=1, YDCPG_DIM%KFLEVG
          DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
            IF (ZTKEM(JLON,JLEV) > 3300._JPRB) THEN
              WRITE (NULOUT,*) 'TKE > 3300 ! '
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    CALL VDFHGHTHL(YDMODEL%YRML_PHY_G%YRVDF, YDMODEL%YRML_PHY_SLIN%YREPHLI, YDMODEL%YRML_PHY_EC%YRECUMF,      &
    & YDMODEL%YRML_PHY_EC%YREPHY, YDPARAR, YDCPG_DIM%KSTEP, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, &
    & YDCPG_DIM%KFLEVG, INDRAFT, PDTPHY, YLMF_PHYS_STATE%U, YLMF_PHYS_STATE%V, ZTM, ZQVM, ZQCM, ZQIM,         &
    & YDVARS%A%T1, YLMF_PHYS_STATE%YCPG_PHY%PRE, YLMF_PHYS_STATE%YCPG_PHY%PREF, ZAPHIFM, ZAPHIM,              &
    & ZZS_FTH_, ZZS_FRV_, ZZS_FU_, ZZS_FV_, ZMF_UP, ZTHETAL_UP, ZQT_UP, ZTHTV_UP, ZQC_UP, ZQI_UP,             &
    & ZU_UP, ZV_UP, ZGP2DSPP, NGFL_EZDIAG, ZP1EZDIAG, ZTENDQVUP, ZTENDTUP, ZSURFPREP, ZSURFSNOW,              &
    & ZUPGENL, ZUPGENN, ZCLFR, ZLENGTH_M, ZLENGTH_H, ZTKEEDMF)

    !  tendtup, tendqvup  tendencies for non-conserved AROME
    !  variables due to updraft precipitation/snow (and its evaporation)
    DO JLEV = 2 ,YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV) + ZTENDTUP(JLON,JLEV)
        ZTENDR(JLON,JLEV,1)=ZTENDR(JLON,JLEV,1) + ZTENDQVUP(JLON,JLEV)
      ENDDO
    ENDDO
    
    IF(LTOTPREC)THEN
      !Add rain and snow tendencies from the sub-grid scheme to tendencies and sources,
      !at all vertical levels, instead of diagnosing only surface precip. 
      ZSURFPREP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0.0_JPRB
      ZSURFSNOW(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0.0_JPRB
      DO JLEV= 1, YDCPG_DIM%KFLEVG
        DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
          !Add rain and snow to sources:
          ZRS_(JLON,JLEV,3)=ZRS_(JLON,JLEV,3)+ZUPGENL(JLON,JLEV)
          ZRS_(JLON,JLEV,5)=ZRS_(JLON,JLEV,5)+ZUPGENN(JLON,JLEV)
          ZTHS__(JLON,JLEV)=ZTHS__(JLON,JLEV)+ZTENDTUP(JLON,JLEV)*(RATM/&
           & YLMF_PHYS_STATE%YCPG_PHY%PREF(JLON,JLEV))**(RD/RCPD)
          !Update rain/snow tendencies:
          ZTENDR(JLON,JLEV,3)=ZTENDR(JLON,JLEV,3)+ZUPGENL(JLON,JLEV)
          ZTENDR(JLON,JLEV,5)=ZTENDR(JLON,JLEV,5)+ZUPGENN(JLON,JLEV)
        ENDDO
      ENDDO 
    ENDIF

  ELSE
    IDRAFT=3 ! only a wet updraft
    INDRAFT=1
    ZSURFPREP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0._JPRB
    ZSURFSNOW(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0._JPRB
  ENDIF

  DO JDRAFT=IDRAFT,3

    ! No need to swapp because IN and OUT are never needed simultaneously

    !!! Call mass fluxes computations
    ! If CMF_UPDRAFT='DUAL', the updraft characteritics are already computed and will be passed as inputs of SHALLOW_MF
    ! if not, they will be computed in SHALLOW_MF itself (from Méso-NH type routines)

    ! JDRAFT=2 : dry updraft
    ! JDRAFT=3 : wet updraft

    IF (CMF_UPDRAFT=='DUAL') THEN
      ! Goes from one of the updraft from the IFS level world to the Méso-NH level world
      ! go from q to r)
      DO JLEV = 1,YDCPG_DIM%KFLEVG
        DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZMF_UP__(JLON,JLEV) = ZMF_UP(JLON,JLEV,JDRAFT)
          ZZU_UP_(JLON,JLEV) = ZU_UP(JLON,JLEV,JDRAFT)
          ZZV_UP_(JLON,JLEV) = ZV_UP(JLON,JLEV,JDRAFT)
          ZTHETAL_UP_(JLON,JLEV) = ZTHETAL_UP(JLON,JLEV,JDRAFT)
          ZTHETAV_UP_(JLON,JLEV) = ZTHTV_UP(JLON,JLEV,JDRAFT)
          ZRT_UP_(JLON,JLEV)  = ZQT_UP(JLON,JLEV,JDRAFT)/(1.-ZQT_UP(JLON,JLEV,JDRAFT))
          ZRC_UP_(JLON,JLEV)  = ZQC_UP(JLON,JLEV,JDRAFT)/(1.-ZQT_UP(JLON,JLEV,JDRAFT))
          ZRI_UP_(JLON,JLEV)  = ZQI_UP(JLON,JLEV,JDRAFT)/(1.-ZQT_UP(JLON,JLEV,JDRAFT))
        ENDDO
      ENDDO
      ZZW_UP_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:IKT)=0._JPRB
      ZZFRAC_UP_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:IKT)=0._JPRB
      IF (LHARATU) THEN
        DO JLEV = 1,YDCPG_DIM%KFLEVG
          DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
            ZLENGTHM__(JLON,JLEV) = MAX(0.01_JPRB,ZLENGTH_M(JLON,JLEV))
            ZLENGTHH__(JLON,JLEV) = MAX(0.01_JPRB,ZLENGTH_H(JLON,JLEV))
            ! TKE should be bigger than a minimum value:
            ZTKEEDMFS(JLON,JLEV) = MAX(ZTKEEDMF(JLON,JLEV),PPTKEMIN)*ZINVDT
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    IF(MOD(YDCPG_DIM%KSTEP+1,NPRINTFR)==0) THEN
      WRITE(NULOUT,*)"apres surface zsfth zsfrv",ZSFTH_(NPTP),ZSFRV_(NPTP)
    ENDIF

    DO JLEV = 1, YDCPG_DIM%KFLEVG 
      ZRC_MF_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=0._JPRB
      ZRI_MF_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=0._JPRB
      ZCF_MF_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=0._JPRB
    ENDDO

    IF (JDRAFT == IDRAFT) THEN
      ! Fill the sum at the first iteration
      ZARG_FLXZTHVMF_ => ZFLXZTHVMF_SUM__(:,1:YDCPG_DIM%KFLEVG)
    ELSE
      ! increment
      ZARG_FLXZTHVMF_ => ZFLXZTHVMF_(:,1:YDCPG_DIM%KFLEVG)
    ENDIF

    CALL ARO_SHALLOW_MF (KKL=IKL, KLON=YDCPG_DIM%KFDIA, KLEV=YDCPG_DIM%KFLEVG, KRR=NRR, KRRL=NRRL, KRRI=NRRI,                  &
    & KSV=NGFL_EXT, HMF_UPDRAFT=CMF_UPDRAFT, HMF_CLOUD=CMF_CLOUD, HFRAC_ICE=CFRAC_ICE_SHALLOW_MF, OMIXUV=LMIXUV,               &
    & ONOMIXLG=.FALSE., KSV_LGBEG=0, KSV_LGEND=0, KTCOUNT=YDCPG_DIM%KSTEP+1, PTSTEP=ZDT, PZZ=ZZZ_, PZZF=ZZZ_F_,                &
    & PDZZF=ZDZZ_F_, PRHODJ=ZRHODJM__(:, 1:YDCPG_DIM%KFLEVG), PRHODREF=ZRHODREFM__(:, 1:YDCPG_DIM%KFLEVG),                     &
    & PPABSM=ZPABSM__(:, 1:YDCPG_DIM%KFLEVG), PEXNM=ZEXNREFM_, PSFTH=ZSFTH_, PSFRV=ZSFRV_, PTHM=ZTHM__(:, 1:YDCPG_DIM%KFLEVG), &
    & PRM=ZRM_, PUM=ZUM__(:, 1:YDCPG_DIM%KFLEVG), PVM=ZVM__(:, 1:YDCPG_DIM%KFLEVG), PTKEM=ZTKEM__(:, 1:YDCPG_DIM%KFLEVG),      &
    & PSVM=ZSVM_, PDUDT_MF=ZMFUS_, PDVDT_MF=ZMFVS_, PDTHLDT_MF=ZTHLS_, PDRTDT_MF=ZRTS_, PDSVDT_MF=ZSVXXX_,                     &
    & PSIGMF=ZSIGMF_, PRC_MF=ZRC_MF_, PRI_MF=ZRI_MF_, PCF_MF=ZCF_MF_, PFLXZTHVMF=ZARG_FLXZTHVMF_, PTHL_UP=ZTHETAL_UP_,         &
    & PRT_UP= ZRT_UP_, PRV_UP=ZZRV_UP_, PRC_UP=ZRC_UP_, PRI_UP=ZRI_UP_, PU_UP=ZZU_UP_, PV_UP=ZZV_UP_,                          &
    & PTHV_UP=ZTHETAV_UP_, PW_UP=ZZW_UP_, PFRAC_UP=ZZFRAC_UP_, PEMF=ZMF_UP__(:, 1:YDCPG_DIM%KFLEVG))

    IF (JDRAFT > IDRAFT) THEN
      ! Add increment
      ZFLXZTHVMF_SUM__(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=ZFLXZTHVMF_SUM__(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)+ZFLXZTHVMF_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
    ENDIF

    ! traitement des sorties pour repasser dans le monde Aladin

    IF ((CMF_CLOUD=='DIRE'.OR.CMF_CLOUD=='BIGA').AND.JDRAFT==3) THEN
      ! sauvegarde pour le schema de nuage
      DO JLEV = 1,YDCPG_DIM%KFLEVG
        ZP1EZDIAG(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV,1)=ZRC_MF_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)
        ZP1EZDIAG(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV,3)=ZRI_MF_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)
        ZP1EZDIAG(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV,2)=ZCF_MF_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)
      ENDDO
    ENDIF
    DO JLEV = 1,YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZUS__(JLON,JLEV)=ZUS__(JLON,JLEV)+ZMFUS_(JLON,JLEV)
        ZVS__(JLON,JLEV)=ZVS__(JLON,JLEV)+ZMFVS_(JLON,JLEV)
        ZTHS__(JLON,JLEV)=ZTHS__(JLON,JLEV)+ZTHLS_(JLON,JLEV)
        ZRS_(JLON,JLEV,1)=ZRS_(JLON,JLEV,1)+ZRTS_(JLON,JLEV)
        !calcul de tendance et inversion des niveaux pour le vent horizontal
        YDMF_PHYS%OUT%TENDU(JLON,JLEV)=YDMF_PHYS%OUT%TENDU(JLON,JLEV)+ZMFUS_(JLON,JLEV)
        YDMF_PHYS%OUT%TENDV(JLON,JLEV)=YDMF_PHYS%OUT%TENDV(JLON,JLEV)+ZMFVS_(JLON,JLEV)
        !conversion de la tendance de theta en tendance de T et inversion niveau
        ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV)+ZTHLS_(JLON,JLEV)*ZEXNREFM_(JLON,JLEV)
        ZTENDTT(JLON,JLEV)=ZTENDTT(JLON,JLEV)+ZTHLS_(JLON,JLEV)
        !inversion niveaux tendances des ri et conversion en qi en multipliant par qd
        ZTENDR(JLON,JLEV,1) = ZTENDR(JLON,JLEV,1)+ZRTS_(JLON,JLEV)*ZQDM(JLON,JLEV)
      ENDDO
    ENDDO

  ENDDO ! JDRAFT

ENDIF ! LMFSHAL


!    ------------------------------------------------------------------
!     10 - TURBULENCE.
!     --------------------------------------------------------------------

IF (LTURB) THEN

  ! Swapp because IN and OUT might be needed simultaneously (though commented out)
  CALL SWAP_LIMAS

  ! Swapp because IN and OUT will be needed simultaneously
  CALL SWAP_SVS
  ! well let's keep the copy, though for now  OUT=IN anyway.
  IF (NGFL_EXT /=0 ) THEN
    ZSVSIN_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,1:NGFL_EXT)=ZSVS_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,1:NGFL_EXT)
  ENDIF

  !prints
  IF(MOD(YDCPG_DIM%KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'avant d entrer dans turb sous apl_arome U'
    WRITE(NULOUT,*)MAXVAL(ZUM__(:,IKB)), MINVAL(ZUM__(:,IKB))
    WRITE(NULOUT,*)'avant d entrer dans turb sous apl_arome V'
    WRITE(NULOUT,*)MAXVAL(ZVM__(:,IKB)), MINVAL(ZVM__(:,IKB))
    WRITE(NULOUT,*)'avant d entrer dans turb sous apl_arome W'
    WRITE(NULOUT,*)MAXVAL(ZWM__(:,IKB)), MINVAL(ZWM__(:,IKB))
    WRITE(NULOUT,*)'avant d entrer dans turb sous apl_arome TKE'
    WRITE(NULOUT,*)MAXVAL(ZTKEM__(:,IKB)), MINVAL(ZTKEM__(:,IKB))
    DO JLEV=1,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,*)JLEV,ZUM__(NPTP,JLEV),ZVM__(NPTP,JLEV),ZWM__(NPTP,JLEV),ZTKEM__(NPTP,JLEV)
    ENDDO
    WRITE(NULOUT,*)'u v w tke a S'
    DO JLEV=1,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,*)JLEV,ZUS__(NPTP,JLEV),ZVS__(NPTP,JLEV),ZWS__(NPTP,JLEV),ZTKES_(NPTP,JLEV)
    ENDDO
    WRITE(NULOUT,*)'ZTHS__ avant turb'
    DO JLEV=1,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,*)JLEV,ZTHS__(NPTP,JLEV)
    ENDDO
  ENDIF

!!$
!!$! Allocation des variables SV (NGFL_EXT + NLIMA)
!!$  KSV_TURB=NGFL_EXT+NLIMA
!!$!
!!$  IF (NGFL_EXT/=0) THEN
!!$     DO JGFL=1,NGFL_EXT
!!$        DO JLON=KIDIA,KFDIA
!!$           ZSFTURB(JLON,JGFL)=ZSFSV_(JLON,JGFL)
!!$           DO JLEV = 1, KLEV
!!$              ZTURBM(JLON,JLEV,JGFL)=ZSVM_(JLON,1,JLEV,JGFL)
!!$              ZTURBS(JLON,JLEV,JGFL)=ZSVSIN_(JLON,1,JLEV,JGFL)
!!$           ENDDO
!!$        ENDDO
!!$     ENDDO
!!$  ENDIF
!!$!
!!$  IF (NLIMA/=0) THEN
!!$     DO JGFL=1,NLIMA
!!$        DO JLON=KIDIA,KFDIA
!!$           ZSFTURB(JLON,NGFL_EXT+JGFL)=0.
!!$           DO JLEV = 1, KLEV
!!$              ZTURBM(JLON,JLEV,NGFL_EXT+JGFL)=ZLIMAM_(JLON,JLEV,JGFL)
!!$              ZTURBS(JLON,JLEV,NGFL_EXT+JGFL)=ZLIMASIN_(JLON,JLEV,JGFL)
!!$           ENDDO
!!$        ENDDO
!!$     ENDDO
!!$  ENDIF

  ! Input variable indeed. REK
  ZSFSVLIMA_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:NLIMA)=0._JPRB

  ! 10.2 calcul TURB
  ZZTOP_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=ZAPHIM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,0)*ZINVG

  IF (LGRADHPHY) THEN
  !   
    DO JLEV = 1,YDCPG_DIM%KFLEVG
      DO JGR=1,NGRADIENTS
        ZTURB3D__(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV,JGR)=YDMF_PHYS%GRA%G(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JGR,JLEV)
      ENDDO
    ENDDO
  
  ENDIF

! Appel avec les arguments modifiés pour variables LIMA :
! KSV_TURB, ZSFTURB, ZTURBM, ZTURBS, ZTENDSV_TURB
  CALL ARO_TURB_MNH(KKA=IKA, KKU=IKU, KKL=IKL, KLON=YDCPG_DIM%KFDIA, KLEV=YDCPG_DIM%KFLEVG, KRR=NRR,                 &
  & KRRL=NRRL, KRRI= NRRI, KSV=NLIMA, KTCOUNT=YDCPG_DIM%KSTEP+1, KGRADIENTS=NGRADIENTS, LDHARATU=LHARATU,            &
  & PTSTEP=ZDT, PZZ=ZZZ_, PZZF=ZZZ_F_, PZZTOP=ZZTOP_, PRHODJ=ZRHODJM__, PTHVREF=ZTHVREFM__, PRHODREF=ZRHODREFM__,    &
  & HINST_SFU='M', HMF_UPDRAFT=CMF_UPDRAFT, PSFTH=ZSFTH_, PSFRV=ZSFRV_, PSFSV=ZSFSVLIMA_, PSFU=ZSFU_,                &
  & PSFV=ZSFV_, PPABSM=ZPABSM__, PUM=ZUM__, PVM=ZVM__, PWM=ZWM__, PTKEM=ZTKEM__, PEPSM=ZEPSM, PSVM=ZLIMAM_,          &
  & PSRCM=ZSRCS__, PTHM=ZTHM__, PRM=ZRM_, PRUS=ZUS__, PRVS=ZVS__, PRWS=ZWS__, PRTHS=ZTHS__, PRRS=ZRS_,               &
  & PRSVSIN=ZLIMASIN_, PRSVS=ZLIMAS_, PRTKES=ZTKES_, PRTKES_OUT=ZTKES_OUT__, PREPSS=ZEPSS, PHGRAD=ZTURB3D__,         &
  & PSIGS=ZSIGS__, OSUBG_COND=LOSUBG_COND, PFLXZTHVMF=ZFLXZTHVMF_SUM__, PLENGTHM=ZLENGTHM__, PLENGTHH=ZLENGTHH__,    &
  & MFMOIST=ZMF_UP__, PDRUS_TURB=ZTENDU_TURB__, PDRVS_TURB=ZTENDV_TURB__, PDRTHLS_TURB=ZTENDTHL_TURB__,              &
  & PDRRTS_TURB=ZTENDRT_TURB__, PDRSVS_TURB=ZTENDSV_TURBLIMA_, PDP=ZDP__, PTP=ZTP__, PTPMF=ZTPMF__, PTDIFF=ZTDIFF__, &
  & PTDISS=ZTDISS__, PEDR=ZEDR__, YDDDH=YDDDH, YDLDDH=YDMODEL%YRML_DIAG%YRLDDH, YDMDDH=YDMODEL%YRML_DIAG%YRMDDH      &
  & )


! Séparation des variables SV (NGFL_EXT + NLIMA)
!!$  IF (NGFL_EXT/=0) THEN
!!$     DO JGFL=1,NGFL_EXT
!!$        DO JLON=KIDIA,KFDIA
!!$           ZSFSV_(JLON,JGFL)=ZSFTURB(JLON,JGFL)
!!$           DO JLEV = 1, KLEV
!!$              ZSVM_(JLON,1,JLEV,JGFL)=ZTURBM(JLON,JLEV,JGFL)
!!$              ZSVS_(JLON,1,JLEV,JGFL)=ZTURBS(JLON,JLEV,JGFL)
!!$           ENDDO
!!$        ENDDO
!!$     ENDDO
!!$  ENDIF
!!$!
!!$  IF (NLIMA/=0) THEN
!!$     DO JGFL=1,NLIMA
!!$        DO JLON=KIDIA,KFDIA
!!$           DO JLEV = 1, KLEV
!!$              ZLIMAM_(JLON,JLEV,JGFL)=ZTURBM(JLON,JLEV,NGFL_EXT+JGFL)
!!$              ZLIMAS_(JLON,JLEV,JGFL)=ZTURBS(JLON,JLEV,NGFL_EXT+JGFL)
!!$           ENDDO
!!$        ENDDO
!!$     ENDDO
!!$  ENDIF


  DO JLEV = 1 , YDCPG_DIM%KFLEVG
     YDMF_PHYS%OUT%EDR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=ZEDR__(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)
  ENDDO
   
  IF (LFLEXDIA) THEN
    DO JLEV = 1,YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZDP__(JLON,JLEV)=ZDP__(JLON,JLEV)*YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
        ZTP__(JLON,JLEV)=(ZTP__(JLON,JLEV)-ZTPMF__(JLON,JLEV))*YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
        ZTPMF__(JLON,JLEV)=ZTPMF__(JLON,JLEV)*YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
        ZTDIFF__(JLON,JLEV)=ZTDIFF__(JLON,JLEV)*YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
        ZTDISS__(JLON,JLEV)=ZTDISS__(JLON,JLEV)*YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
      ENDDO
    ENDDO
    IF (LDDH_OMP) THEN
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZDP__(:, 1:YDCPG_DIM%KFLEVG), 'TKEPRDY', YDDDH&
      & )
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZTP__(:, 1:YDCPG_DIM%KFLEVG), 'TKEPRTH', YDDDH&
      & )
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZTPMF__(:, 1:YDCPG_DIM%KFLEVG), 'TKEPRTHMF', &
      & YDDDH)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZTDIFF__(:, 1:YDCPG_DIM%KFLEVG), 'TKEDIFF', &
      & YDDDH)
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH, ZTDISS__(:, 1:YDCPG_DIM%KFLEVG), 'TKEDISS', &
      & YDDDH)
    ELSE
      CALL ADD_FIELD_3D(YLDDH, ZDP__(:, 1:YDCPG_DIM%KFLEVG), 'TKEPRDY', 'T', 'ARO', .TRUE., .TRUE.)
      CALL ADD_FIELD_3D(YLDDH, ZTP__(:, 1:YDCPG_DIM%KFLEVG), 'TKEPRTH', 'T', 'ARO', .TRUE., .TRUE.)
      CALL ADD_FIELD_3D(YLDDH, ZTPMF__(:, 1:YDCPG_DIM%KFLEVG), 'TKEPRTHMF', 'T', 'ARO', .TRUE., &
      & .TRUE.)
      CALL ADD_FIELD_3D(YLDDH, ZTDIFF__(:, 1:YDCPG_DIM%KFLEVG), 'TKEDIFF', 'T', 'ARO', .TRUE., &
      & .TRUE.)
      CALL ADD_FIELD_3D(YLDDH, ZTDISS__(:, 1:YDCPG_DIM%KFLEVG), 'TKEDISS', 'T', 'ARO', .TRUE., &
      & .TRUE.)
    ENDIF
  ENDIF 

  IF(MOD(YDCPG_DIM%KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'u v w a S apres turb'
    DO JLEV=1,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,*)JLEV,ZUS__(NPTP,JLEV),ZVS__(NPTP,JLEV),ZWS__(NPTP,JLEV),ZTKES_OUT__(NPTP,JLEV)
    ENDDO
    WRITE(NULOUT,*)'THS TKES SIGS apres turb'
    DO JLEV=1,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,*)JLEV,ZTHS__(NPTP,JLEV),ZTKES_OUT__(NPTP,JLEV),ZSIGS__(NPTP,JLEV)
    ENDDO
  ENDIF

  ! avance temporelle et inversion niveau pour ZSIGS__
  IF (LOSUBG_COND .AND. LOSIGMAS) THEN
    IF (CMF_CLOUD=='DIRE'.OR.CMF_CLOUD=='BIGA'.OR.CMF_CLOUD=='NONE') THEN
      DO JLEV = 1,YDCPG_DIM%KFLEVG
        YDVARS%SRC%T1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)=ZSIGS__(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV)
      ENDDO
    ELSEIF (CMF_CLOUD=='STAT') THEN
      DO JLEV = 1,YDCPG_DIM%KFLEVG
        DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          YDVARS%SRC%T1(JLON,JLEV)=SQRT(ZSIGS__(JLON,JLEV)**2+ZSIGMF_(JLON,JLEV)**2 )
        ENDDO
      ENDDO
    ENDIF
  ENDIF


  !10.3. traitement des sorties pour repasser dans le monde Aladin
  !calcul de tendance et inversion des niveaux pour le vent horizontal et la TKE

  DO JLEV = 1,YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%TENDU(JLON,JLEV)=YDMF_PHYS%OUT%TENDU(JLON,JLEV)+ZTENDU_TURB__(JLON,JLEV)
      YDMF_PHYS%OUT%TENDV(JLON,JLEV)=YDMF_PHYS%OUT%TENDV(JLON,JLEV)+ZTENDV_TURB__(JLON,JLEV)
      ! for the moment, turbulence do not compute w tendency:
      ZTENDW(JLON,JLEV)=0.0_JPRB
      ! PTENDW(JLON,JLEV)+(ZWS__(JLON,JLEV)-&
      ! & ZWS_AVE(JLON,1,JLEV))
      !conversion de la tendance de theta en tendance de T et inversion niveau
      ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV)+ZTENDTHL_TURB__(JLON,JLEV)*ZEXNREFM_(JLON,JLEV)
  !inversion niveaux tendances des rv et conversion en qv en multipliant par qd
      ZTENDR(JLON,JLEV,1)= ZTENDR(JLON,JLEV,1)+ZTENDRT_TURB__(JLON,JLEV)*ZQDM(JLON,JLEV)
    ENDDO
  ENDDO


  IF (LHARATU) THEN
    DO JLEV = 1,YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZTENDTKE(JLON,JLEV)=ZTENDTKE(JLON,JLEV)+(ZTKEEDMFS(JLON,JLEV)-ZTKES_(JLON,JLEV))
      ENDDO
    ENDDO
  ELSE
    DO JLEV = 1,YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
         ZTENDTKE(JLON,JLEV)=ZTENDTKE(JLON,JLEV)+(ZTKES_OUT__(JLON,JLEV)-ZTKES_(JLON,JLEV))
      ENDDO
    ENDDO
  ENDIF

  DO JGFL=1,NGFL_EXT
    DO JLEV = 1,YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZTENDEXT(JLON,JLEV,JGFL)=ZTENDEXT(JLON,JLEV,JGFL)+(ZSVS_(JLON,JLEV,JGFL)-ZSVSIN_(JLON,JLEV,JGFL))
      ENDDO
    ENDDO
  ENDDO

! Tendances LIMA
  DO JGFL=1,NLIMA
    DO JLEV = 1, YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZTENDLIMA(JLON,JLEV,JGFL)=ZTENDLIMA(JLON,JLEV,JGFL)+ZTENDSV_TURBLIMA_(JLON,JLEV,JGFL)
!        PTENDLIMA(JLON,JLEV,:)=PTENDLIMA(JLON,JLEV,:)+ (ZLIMAS_(JLON,JLEV,:)-ZLIMASIN_(JLON,JLEV,:))
      ENDDO
    ENDDO
  ENDDO

ENDIF
!    ------------------------------------------------------------------
!     11 - MICROPHYSIQUE. 
!     --------------------------------------------------------------------

IF (LMICRO) THEN

  ! Swap pointers because input values of THS and RS should be saved
  CALL SWAP_THS
  CALL SWAP_RS
  CALL SWAP_LIMAS

  ! for now a copy is needed (see below, inside). I don't like than :-( REK
  ZTHS__(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=ZTHSIN_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  ZRS_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,1:NRR)=ZRSIN_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,1:NRR)
  ! for now a copy is needed (see below, inside). I don't like than :-( REK
  ZLIMAS_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,1:NLIMA)=ZLIMASIN_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,1:NLIMA)
     
  !prints
  IF (MOD(YDCPG_DIM%KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'avant rain_ice sous apl_arome'
    WRITE(NULOUT,*)'JLEV   ZZZ_F_      ZZZ_      ZRHODREF',&
     & '    ZRHODJ      ZPABSM__        ZTHSIN_       ZTHM__      '   
    DO JLEV=1,YDCPG_DIM%KFLEVG+1 
      WRITE(NULOUT,'(I2,X,7F10.3)')JLEV,ZZZ_F_(NPTP,JLEV),ZZZ_(NPTP,JLEV), ZRHODREFM__(NPTP,JLEV),&
       & ZRHODJM__(NPTP,JLEV), ZPABSM__(NPTP,JLEV), ZTHSIN_(NPTP,JLEV), ZTHM__(NPTP,JLEV)  
    ENDDO 
    WRITE(NULOUT,*)'JLEV        PDELPM        ZPABSM__         ZEXNREF','          ZSIGS__'
    DO JLEV=2,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,'(I2,X,4f10.3)')JLEV, YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP(NPTP,JLEV),&
       & ZPABSM__(NPTP,JLEV),ZEXNREFM_(NPTP,JLEV),ZSIGS__(NPTP,JLEV)  
    ENDDO
    WRITE(NULOUT,*)'JLEV    PTM       PRM          PCPM'
    DO JLEV=1,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,'(I2,X,3f10.3)')JLEV,ZTM(NPTP,YDCPG_DIM%KFLEVG+1-JLEV), ZRHM(NPTP,YDCPG_DIM%KFLEVG+1-JLEV) ,ZCPM(NPTP,YDCPG_DIM%KFLEVG+1-JLEV)  
    ENDDO
    WRITE (NULOUT,*)'JLEV  rhoQv  rhoQc   rhoQr   rhoQi   rhoQs   rhoQg'
    DO JLEV=1,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,'(I2,X,6E11.4)')JLEV,ZRM_(NPTP,JLEV,1), ZRM_(NPTP,JLEV,2),&
       & ZRM_(NPTP,JLEV,3),ZRM_(NPTP,JLEV,4),ZRM_(NPTP,JLEV,5), ZRM_(NPTP,JLEV,6)  
    ENDDO
    WRITE (NULOUT,*)'JLEV  ZRSQv  ZRSQc   ZRSQr   ZRSQi   ZRSQs   ZRSQg'
    DO JLEV=1,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,'(I2,X,6E11.4)')JLEV,ZRS_(NPTP,JLEV,1), ZRS_(NPTP,JLEV,2),&
       & ZRSIN_(NPTP,JLEV,3),ZRSIN_(NPTP,JLEV,4),ZRSIN_(NPTP,JLEV,5), ZRSIN_(NPTP,JLEV,6)  
    ENDDO
    WRITE(NULOUT,*)'ZDT=',ZDT
    WRITE(NULOUT,*)'NRR and co',NRR,YDCPG_DIM%KSTEP+1,NSPLITR,LOSUBG_COND, LOSIGMAS, CSUBG_AUCV_RC,LOWARM  
  ENDIF
  

  ZSEA_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0.0_JPRB
  IF (LOLSMC) THEN
    DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
      IF (YDMF_PHYS_SURF%GSD_VF%PLSM(JLON) < 0.5) THEN
        ZSEA_(JLON) = 1.0_JPRB
      ENDIF
    ENDDO
  ENDIF
         
  IF (LOTOWNC) THEN
    ZTOWN_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA) = ZTOWNS_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
  ELSE
    ZTOWN_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0.0_JPRB
  ENDIF  

  IF (CMICRO == 'LIMA') THEN

    IF (LTURB) THEN
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        DO JLEV=1,YDCPG_DIM%KFLEVG
          ZWNU_(JLON,JLEV) = ZWM__(JLON,JLEV) + 0.66*SQRT(ZTKEM__(JLON,JLEV))
        ENDDO
      ENDDO
      ZPTRWNU_ => ZWNU_(1:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
    ELSE
      ZPTRWNU_ => ZWM__(1:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
    ENDIF
    CALL ARO_LIMA(YDCPG_DIM%KFLEVG, IKU, IKL, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, NRR, NLIMA, YDCPG_DIM%KSTEP+1,           &
    & NSPLITR, NSPLITG, ZDT, ZDZZ_, ZRHODJM__(:, 1:YDCPG_DIM%KFLEVG), ZRHODREFM__(:, 1:YDCPG_DIM%KFLEVG),                 &
    & ZEXNREFM_, ZPABSM__(:, 1:YDCPG_DIM%KFLEVG), ZPTRWNU_, ZTHM__(:, 1:YDCPG_DIM%KFLEVG), ZRM_,                          &
    & ZLIMAM_, ZTHS__(:, 1:YDCPG_DIM%KFLEVG), ZRS_, ZLIMAS_, ZEVAP_, ZINPRR_NOTINCR_,                                     &
    & ZINPRS_NOTINCR_, ZINPRG_NOTINCR_, ZINPRH_NOTINCR_, ZPFPR_, YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH&
    &  )
  ELSE
    CALL ARO_RAIN_ICE (YDCPG_DIM%KFLEVG, IKU, IKL, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, NRR, YDCPG_DIM%KSTEP+1,                  &
    & NSPLITR, NGFL_EZDIAG, LOSUBG_COND, CSUBG_AUCV_RC, LOSEDIC, CSEDIM, CMICRO, ZDT, ZDZZ_, ZRHODJM__(:, 1:YDCPG_DIM%KFLEVG), &
    & ZRHODREFM__(:, 1:YDCPG_DIM%KFLEVG), ZEXNREFM_, ZPABSM__(:, 1:YDCPG_DIM%KFLEVG), ZTHM__(:, 1:YDCPG_DIM%KFLEVG),           &
    & ZRM_, ZSIGS__(:, 1:YDCPG_DIM%KFLEVG), ZNEBMNH_, ZTHS__(:, 1:YDCPG_DIM%KFLEVG), ZRS_, ZEVAP_,                             &
    & ZCIT_, LOWARM, ZSEA_, ZTOWN_, LOCND2, LGRSN, ZINPRR_NOTINCR_, ZINPRS_NOTINCR_, ZINPRG_NOTINCR_,                          &
    & ZINPRH_NOTINCR_, ZPFPR_, ZGP2DSPP, ZP1EZDIAG, YDDDH, YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH                  &
    & )
  ENDIF

  DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZINPRR_(JLON)=ZINPRR_(JLON)+ZINPRR_NOTINCR_(JLON)
    ZINPRS_(JLON)=ZINPRS_(JLON)+ZINPRS_NOTINCR_(JLON)
    ZINPRG_(JLON)=ZINPRG_(JLON)+ZINPRG_NOTINCR_(JLON)
    ZINPRH_(JLON)=ZINPRH_(JLON)+ZINPRH_NOTINCR_(JLON)
  ENDDO

  !conversion de la tendance de theta en tendance de T et inversion niveau
  DO JLEV = 1,YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZTENDT(JLON,JLEV)= ZTENDT(JLON,JLEV)+(ZTHS__(JLON,JLEV)-ZTHSIN_(JLON,JLEV))*ZEXNREFM_(JLON,JLEV)  
      ZTENDTT(JLON,JLEV)= ZTENDTT(JLON,JLEV)+ZTHS__(JLON,JLEV)-ZTHSIN_(JLON,JLEV)
    ENDDO
  ENDDO
  
  !inversion niveaux tendances des ri et conversion en qi en multipliant par qd
  DO JR=1,NRR
    DO JLEV=1,YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZTENDR(JLON,JLEV,JR)=ZTENDR(JLON,JLEV,JR)+(ZRS_(JLON,JLEV,JR)-ZRSIN_(JLON,JLEV,JR))*ZQDM(JLON,JLEV)  
      ENDDO
    ENDDO
  ENDDO

  ! Tendances des variables LIMA
  DO JGFL=1,NLIMA
    DO JLEV = 1, YDCPG_DIM%KFLEVG
      DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZTENDLIMA(JLON,JLEV,JGFL)=ZTENDLIMA(JLON,JLEV,JGFL)+(ZLIMAS_(JLON,JLEV,JGFL)-ZLIMASIN_(JLON,JLEV,JGFL))  
      ENDDO
    ENDDO
  ENDDO

  IF (LINTFLEX) THEN
    !inversion of levels of upper-air precipitation
    DO JR=2,NRR ! no precip for qv
      ZFPR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,0,JR)=0._JPRB  ! zero precip at top of atmosphere
      DO JLEV=1,YDCPG_DIM%KFLEVG
        ZFPR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV,JR)=ZPFPR_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JLEV,JR)
      ENDDO
    ENDDO
  ENDIF

  !store cumulative 3D precipitations for mocage      
  IF (LFPREC3D) THEN
    DO JR=2,NRR ! no precip for qv
      DO JLEV=1,YDCPG_DIM%KFLEVG
        DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZP1EZDIAG(JLON,JLEV,4)=ZP1EZDIAG(JLON,JLEV,4)+ZPFPR_(JLON,JLEV,JR)*1000._JPRB*PDTPHY
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  !prints    
  IF(MOD(YDCPG_DIM%KSTEP+1,NPRINTFR)==0) THEN
    WRITE(NULOUT,*)'PTENDT en sortie de rain_ice'
    WRITE(NULOUT,*)'ZTHS__ en sortie de rain_ice'
    DO JLEV=1,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,*)ZTENDT(NPTP,JLEV),ZTHS__(NPTP,JLEV)
    ENDDO
    WRITE (NULOUT,*)'JLEV  ZTENDQv  ZTZNDQc   ZTENDQr   ZTENDQi' ,'ZTENDQs   ZTENDQg'  
    DO JLEV=1,YDCPG_DIM%KFLEVG
      WRITE(NULOUT,'(I2,X,6E11.4)')JLEV,ZTENDR(NPTP,JLEV,1),ZTENDR(NPTP,JLEV,2),&
       & ZTENDR(NPTP,JLEV,3),ZTENDR(NPTP,JLEV,4),ZTENDR(NPTP,JLEV,5),ZTENDR(NPTP,JLEV,6)  
    ENDDO
    WRITE (NULOUT,*) 'ZSRCS__ et ZNEBMNH_',MAXVAL(ZSRCS__),MAXVAL(ZNEBMNH_) 
  ENDIF
  
  IF (LRDEPOS) THEN
    ISPLITR=NSPLITR
    ! Swapp because IN and OUT will be needed simultaneously
    CALL SWAP_SVM
    CALL ARO_RAINAERO(YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, NGFL_EXT, NRR, PDTPHY, ZSVMIN_, ZZZ_, ZPABSM__(:, 1:YDCPG_DIM%KFLEVG), &
    & ZTHM__(:, 1:YDCPG_DIM%KFLEVG), ZRHODREFM__(:, 1:YDCPG_DIM%KFLEVG), YDCPG_DIM%KSTEP+1, ZRM_,                               &
    & ZEVAP_, ISPLITR, ZSVM_           )
    ! return to tendency
    DO JGFL=1,NGFL_EXT
      DO JLEV = 1,YDCPG_DIM%KFLEVG
        DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZTENDEXT(JLON,JLEV,JGFL)=ZTENDEXT(JLON,JLEV,JGFL)+(ZSVM_(JLON,JLEV,JGFL)-ZSVMIN_(JLON,JLEV,JGFL))*ZINVDT
        ENDDO
      ENDDO
    ENDDO
  ENDIF ! LRDEPOS

ENDIF ! LMICRO

! start LHN F.Meier 2020 ******

LNUDGLHNREAD=.TRUE.
IF(MYPROC==1.AND.YDCPG_DIM%KSTEP==1.AND.LNUDGLH)THEN
  CALL NUDGLHCLIMPROF(YDCPG_DIM%KFLEVG, LNUDGLHNREAD)
ENDIF
! save accumulated precipitation for LHN
IF (LNUDGLH.AND.YDCPG_DIM%KSTEP == NSTARTNUDGLH.AND.NSTARTNUDGLH > 0) THEN
  !IF(MYPROC==1) WRITE(NULOUT,*)'save precip for LHN - STEP:',KSTEP, &
  !  & 'NUDGINGINT:',NINTNUDGLH,'NSTARTNUDGLH:',NSTARTNUDGLH
  CALL NUDGLHPRECIP(YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, ZACPRR_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), &
  & ZACPRS_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZACPRG_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), YDCPG_DIM%KBL           &
  & )
ENDIF
ISTEP=YDCPG_DIM%KSTEP-NSTARTNUDGLH
! if LNUDGLH and KSTEP in nudging interval
IF (LNUDGLH.AND.YDCPG_DIM%KSTEP > NSTARTNUDGLH.AND.YDCPG_DIM%KSTEP <= NSTOPNUDGLH) THEN
  ! safe LH profile for step before LHN step
  LLHN=.FALSE.
  IF(MOD(ISTEP+1,NINTNUDGLH)==0) THEN
    CALL NUDGLHPREP(YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, ZTENDTT, YDCPG_DIM%KBL&
    & )
  ENDIF
  ! LHN step
  IF(MOD(ISTEP,NINTNUDGLH)==0) THEN
    !IF(MYPROC==1) WRITE(NULOUT,*)'LH nudging applied - STEP:',KSTEP, &
    !  & 'NUDGINGINT:',NINTNUDGLH
    ! get index for correctly reading observation from array
    ! first two indices are reserved for other LHN stuff
    JLHSTEP=NINT(1.0_JPRB*ISTEP/(NTIMESPLITNUDGLH*NINTNUDGLH))+2
    !IF(MYPROC==1) WRITE(NULOUT,*)'observation array:',JLHSTEP
    ! call nudging routine to modify LHN profile where necessary
    CALL NUDGLH(NGPTOT, NPROMA, NGPBLKS, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,                                                &
    & ZACPRR_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZACPRS_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZACPRG_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                         &
    & ZTENDTT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG), JLHSTEP, YDCPG_DIM%KBL, ZEXNREFM_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG), &
    & .TRUE., LLHN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG), ZPABSM__(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG),                     &
    & ZDT, ZTHM__(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG), ZRM_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG, :),                       &
    & ZQDM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG), ZTENDR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG, 1),                            &
    & NRR, LNUDGLHNREAD)
    !IF(MYPROC==1) WRITE(NULOUT,*)'calling LH successful - convert TH to T and
    !add temperature tendency'
    ! add LHN tendency to physics tendency, limit LHN tendency
    ZMAXTEND=0.0_JPRB
    ZMINTEND=0.0_JPRB
    DO JLEV = 1,YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        IF(LLHN(JLON,JLEV))THEN
          ZTENDTT(JLON,JLEV)=MAX(ZTENDTT(JLON,JLEV),RMINNUDGLH)
          ZTENDTT(JLON,JLEV)=MIN(ZTENDTT(JLON,JLEV),RMAXNUDGLH)
          ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV)+ZTENDTT(JLON,JLEV)*&
           & RAMPLIFY*ZEXNREFM_(JLON,JLEV)
          ZMINTEND=MIN(ZTENDTT(JLON,JLEV),ZMINTEND)
          ZMAXTEND=MAX(ZTENDTT(JLON,JLEV),ZMAXTEND)
          ! keep RH constant if LNUDGLHCOMPT=T
          IF(YLMF_PHYS_STATE%T(JLON,JLEV)>0.01_JPRB.AND.LNUDGLHCOMPT)THEN
            ZTENDR(JLON,JLEV,1)=ZTENDR(JLON,JLEV,1)+RLVTT/RV/((YLMF_PHYS_STATE%T(JLON,JLEV))**2._JPRB)* &
            & ZTENDTT(JLON,JLEV)*RAMPLIFY*ZEXNREFM_(JLON,JLEV)*ZQSAT(JLON,JLEV)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    !IF(ZMINTEND<-0.01)WRITE(*,*)'ZMINTEND',ZMINTEND
    !IF(ZMAXTEND>0.01)WRITE(*,*)'ZMAXTEND',ZMAXTEND
    ! write LH profiles to array to save it for next time step
    CALL NUDGLHPREP(YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, ZTENDTT, YDCPG_DIM%KBL&
    & )
    IF(MYPROC==1) WRITE(NULOUT,*)'calling LH successful finished'
    ! use LHN factor again on following time steps depending on NTAUNUDGLH
  ELSEIF(MOD(ISTEP,NINTNUDGLH)<NTAUNUDGLH.AND.MOD(ISTEP,NINTNUDGLH)>0 &
      & .AND.ISTEP>NINTNUDGLH) THEN
    IF(MYPROC==1)THEN
      WRITE(NULOUT,*)'LH nudging applied-STEP:',YDCPG_DIM%KSTEP,'NUDGINGINT:',NINTNUDGLH
      WRITE(NULOUT,*)'NTAUNUDGLH:',NTAUNUDGLH
    ENDIF
    ! get index for reading correctly most recent obs
    JLHSTEP=2+NINT(1.0_JPRB*(ISTEP-MOD(ISTEP,NINTNUDGLH))/(NTIMESPLITNUDGLH*NINTNUDGLH))
    !IF(MYPROC==1) WRITE(NULOUT,*)'observation array:',JLHSTEP
    ! call nudging routine to modify LHN profile where necessary
    ! LHN factor is not recalculated but might be damped by RDAMPNUDGLH
    CALL NUDGLH(NGPTOT, NPROMA, NGPBLKS, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,                                                &
    & ZACPRR_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZACPRS_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA), ZACPRG_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA),                         &
    & ZTENDTT(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG), JLHSTEP, YDCPG_DIM%KBL, ZEXNREFM_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG), &
    & .FALSE., LLHN(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG), ZPABSM__(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG),                    &
    & ZDT, ZTHM__(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG), ZRM_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG, :),                       &
    & ZQDM(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG), ZTENDR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA, 1:YDCPG_DIM%KFLEVG, 1),                            &
    & NRR, LNUDGLHNREAD)
    !IF(MYPROC==1) WRITE(NULOUT,*)'calling LH successful - convert TH to T and
    !add temperature tendency'
    ! add LHN tendency to physics tendency, limit LHN tendency
    DO JLEV = 1,YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        IF(LLHN(JLON,JLEV))THEN
          ZTENDTT(JLON,JLEV)=MAX(ZTENDTT(JLON,JLEV),RMINNUDGLH)
          ZTENDTT(JLON,JLEV)=MIN(ZTENDTT(JLON,JLEV),RMAXNUDGLH)
          ZTENDT(JLON,JLEV)= ZTENDT(JLON,JLEV)+ZTENDTT(JLON,JLEV)*&
           & RAMPLIFY*ZEXNREFM_(JLON,JLEV)
          ZMINTEND=MIN(ZTENDTT(JLON,JLEV),ZMINTEND)
          ZMAXTEND=MAX(ZTENDTT(JLON,JLEV),ZMAXTEND)
          ! keep RH constant if LNUDGLHCOMPT=T
          IF(YLMF_PHYS_STATE%T(JLON,JLEV)>0.01_JPRB.AND.LNUDGLHCOMPT)THEN
             ZTENDR(JLON,JLEV,1)=ZTENDR(JLON,JLEV,1)+RLVTT/RV/((YLMF_PHYS_STATE%T(JLON,JLEV))**2._JPRB)*&
              & ZTENDTT(JLON,JLEV)*RAMPLIFY*ZEXNREFM_(JLON,JLEV)*ZQSAT(JLON,JLEV)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    !IF(ZMAXTEND>0.01) WRITE(*,*)'ZMAXTEND',ZMAXTEND
    !IF(ZMINTEND<-0.01) WRITE(*,*)'ZMINTEND',ZMINTEND
    ! write LHN profiles to array for next timestep
    CALL NUDGLHPREP(YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, ZTENDTT, YDCPG_DIM%KBL&
    & )
    IF(MYPROC==1) WRITE(NULOUT,*)'calling LH successful finished'
  ENDIF
ENDIF
! **end latent heat nudging***********

    
!    ------------------------------------------------------------------
!     11 - SAVE FIELDS FOR EXT. SURFACE.
!     --------------------------------------------------------------------
!    Cette partie n'est plus necessaire apres branchement de la physique 
!    de surface sous apl_arome

!    ------------------------------------------------------------------
!     12 - CALL CHEMICAL SCHEME.
!     --------------------------------------------------------------------
IF (LUSECHEM) THEN

  ! ANNEE
  IYEAR = NINDAT / 10000
  ! MOIS
  IMONTH = (NINDAT - 10000*IYEAR ) / 100
  ! JOUR DU MOIS
  IDAY = NINDAT - 10000*IYEAR - 100*IMONTH

  DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZLAT_(JLON) = 180. * ASIN(YDGSGEOM%GEMU(JLON)) / (2.*ASIN(1.))
    ZLON_(JLON) = 180. * YDGSGEOM%GELAM(JLON) / (2.*ASIN(1.))
    ZZENITH_(JLON) = ACOS( YDMF_PHYS_TMP%RDG%MU0(JLON) )
    ZZS_(JLON)=YDOROG%OROG(JLON)/RG
    ZALB_UV_(JLON)=ZALBP(JLON,1)
  ENDDO

  ! Swapp because IN and OUT will be needed simultaneously
  CALL SWAP_SVS

  DO JGFL=1,NGFL_EXT
    DO JLEV=1,YDCPG_DIM%KFLEVG
      DO JLON= YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ! modify input
        ZSVSIN_(JLON,JLEV,JGFL)=MAX(0.0_JPRB, ZSVSIN_(JLON,JLEV,JGFL))
      ENDDO
    ENDDO
  ENDDO
  IEZDIAG_CHEM=NGFL_EZDIAG-IOFF_MFSHAL+1
  CALL ARO_MNHC(ZSVSIN_, ZRHODREFM__(:, 1:YDCPG_DIM%KFLEVG), PDTPHY, ZTHM__(:, 1:YDCPG_DIM%KFLEVG), ZPABSM__(:, 1:YDCPG_DIM%KFLEVG),     &
  & ZRM_, ZLAT_, ZLON_, ZALB_UV_, ZZS_, ZZENITH_, ZZZ_, IYEAR, IMONTH, IDAY, REAL(RHGMT, JPRB)+PDTPHY/2._JPRB,                           &
  & YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, NGFL_EXT, NRR, YDCPG_DIM%KSTEP+1, NULOUT, IEZDIAG_CHEM, ZPEZDIAG_(:, :, IOFF_MFSHAL:NGFL_EZDIAG), &
  & ZSVS_ )
 
  ZP1EZDIAG(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,IOFF_MFSHAL:NGFL_EZDIAG)=ZPEZDIAG_(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG,IOFF_MFSHAL:NGFL_EZDIAG)

  !inversion niveau de la tendance des scalaires passifs
  DO JGFL=1,NGFL_EXT
    DO JLEV = 1,YDCPG_DIM%KFLEVG
      DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        ZTENDEXT(JLON,JLEV,JGFL)=ZTENDEXT(JLON,JLEV,JGFL)+(ZSVS_(JLON,JLEV,JGFL)-ZSVSIN_(JLON,JLEV,JGFL))
      ENDDO
    ENDDO
  ENDDO

ENDIF ! LUSECHEM

!    ------------------------------------------------------------------
!     13 - STOCHASTIC PHYSICS : PERTURB TENDENCIES
!     -----------------------------------------------------------------

IF(YSPPT_CONFIG%LSPSDT) THEN

  ZDUMMY(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=0.0_JPRB               ! Dummy nonphys tendency for compatibility with ecmwf stochphy
  CALL SPPTEN (YDMODEL%YRML_PHY_EC%YRECLDP, YGFL, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG,          &
  & 1, PDTPHY, PTSL=YLMF_PHYS_STATE%T, PQSL=YLMF_PHYS_STATE%Q, PA=YDVARS%A%T1, PAP=YLMF_PHYS_STATE%YCPG_PHY%PREF,              &
  & PAPH=YLMF_PHYS_STATE%YCPG_PHY%PRE, PDYN_U=ZDUMMY, PDYN_V=ZDUMMY, PDYN_T=ZDUMMY, PDYN_Q=ZDUMMY, PUNP_U=YDMF_PHYS%OUT%TENDU, &
  & PUNP_V=YDMF_PHYS%OUT%TENDV, PUNP_T=ZTENDT, PUNP_Q=ZTENDR(:, :, 1), PMULNOISE=PGP2DSDT(1, 1, 1), PTENU=YDMF_PHYS%OUT%TENDU, &
  & PTENV=YDMF_PHYS%OUT%TENDV, PTENT=ZTENDT, PTENQ=ZTENDR(:, :, 1) )    ! Out: (u,v,t,qv) total perturbed tendencies
ENDIF

IF(LFORCENL.AND.(YDCPG_DIM%KSTEP*(TSPHY/RHOUR)>=NFORCESTART).AND.&
              & (YDCPG_DIM%KSTEP*(TSPHY/RHOUR)<=NFORCEEND)) THEN
  DO JLEV=1,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%TENDU(JLON,JLEV)=YDMF_PHYS%OUT%TENDU(JLON,JLEV)+AMAGSTOPH_CASBS*YDMF_PHYS%FOR%U(JLON,JLEV)
      YDMF_PHYS%OUT%TENDV(JLON,JLEV)=YDMF_PHYS%OUT%TENDV(JLON,JLEV)+AMAGSTOPH_CASBS*YDMF_PHYS%FOR%V(JLON,JLEV)
      ZTENDT(JLON,JLEV)=ZTENDT(JLON,JLEV)+AMAGSTOPH_CASBS*YDMF_PHYS%FOR%T(JLON,JLEV)
      ZTENDR(JLON,JLEV,1)=ZTENDR(JLON,JLEV,1)+AMAGSTOPH_CASBS*YDMF_PHYS%FOR%Q(JLON,JLEV)
    ENDDO
  ENDDO
ENDIF

!    ------------------------------------------------------------------
!     14 - FINAL CALCULATIONS.
!     --------------------------------------------------------------------

!forcage pour declencher la ligne de grain 
IF (LSQUALL) THEN
  IF (LTWOTL) THEN
    ZDT2=2*ZDT
  ELSE
    ZDT2=ZDT
  ENDIF
  IF((YDCPG_DIM%KSTEP+1)*ZDT2 < 600._JPRB) THEN
    WRITE(NULOUT, *)'refroidissement impose de',NREFROI1,' a ',NREFROI2
    DO JLEV=YDCPG_DIM%KFLEVG,YDCPG_DIM%KFLEVG-20,-1
      ZTENDT(NREFROI1:NREFROI2,JLEV)=-0.01_JPRB
    ENDDO
  ENDIF
ENDIF


!ecriture du buffer
IF(LLMSE.OR.LSFORCS) THEN
  DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    YDCPG_GPAR%INPRR(JLON)=ZINPRR_(JLON)+ZSURFPREP(JLON)/1000._JPRB
    YDCPG_GPAR%INPRS(JLON)=ZINPRS_(JLON)+ZSURFSNOW(JLON)/1000._JPRB
    YDCPG_GPAR%INPRG(JLON)=ZINPRG_(JLON)+ZINPRH_(JLON)
    YDCPG_GPAR%ACPRR(JLON)=YDCPG_GPAR%ACPRR(JLON)+(ZINPRR_(JLON)+ZSURFPREP(JLON)/1000._JPRB)*PDTPHY
    YDCPG_GPAR%ACPRS(JLON)=YDCPG_GPAR%ACPRS(JLON)+(ZINPRS_(JLON)+ZSURFSNOW(JLON)/1000._JPRB)*PDTPHY
    YDCPG_GPAR%ACPRG(JLON)=YDCPG_GPAR%ACPRG(JLON)+(ZINPRG_(JLON)+ZINPRH_(JLON))*PDTPHY
  ENDDO
  YDCPG_GPAR%VTS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=ZTSURF(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
  YDCPG_GPAR%VEMIS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=ZEMIS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
  YDCPG_GPAR%VQS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=ZQS(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)
  DO JSW=1,NSW
    YDCPG_GPAR%ALBDIR(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)=ZALBP(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)
    YDCPG_GPAR%ALBSCA(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)=ZALBD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,JSW)
  ENDDO
ENDIF

IF (LMUSCLFA) CALL ECR1D(NMUSCLFA, 'PCLCT_apl', YDCPG_MISC%CLCT, 1, YDCPG_DIM%KLON)
! initialisations for CFU for Rainfalls
DO JLEV = 0,YDCPG_DIM%KFLEVG
  DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ! conversion from m/s in mm/s
    YDMF_PHYS%OUT%FPLSL(JLON,JLEV)= ZINPRR_(JLON)*1000._JPRB+ZSURFPREP(JLON)
    YDMF_PHYS%OUT%FPLSN(JLON,JLEV)= ZINPRS_(JLON)*1000._JPRB+ZSURFSNOW(JLON)
    YDMF_PHYS%OUT%FPLSG(JLON,JLEV)= ZINPRG_(JLON)*1000._JPRB
    YDMF_PHYS%OUT%FPLSH(JLON,JLEV)= ZINPRH_(JLON)*1000._JPRB
    ! conversion in correct Unit for BADP (same as ALADIN)
    YDMF_PHYS%OUT%STRTU(JLON,JLEV)= ZSFU_(JLON)*ZRHODREFM__(JLON,IKB) 
    YDMF_PHYS%OUT%STRTV(JLON,JLEV)= ZSFV_(JLON)*ZRHODREFM__(JLON,IKB) 
  ENDDO
ENDDO
!Hail diagnostic
YDMF_PHYS%OUT%DIAGH(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA)=0._JPRB
IF (LXXDIAGH) THEN
  DO JLEV=1,YDCPG_DIM%KFLEVG
    DO JLON=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%DIAGH(JLON)=YDMF_PHYS%OUT%DIAGH(JLON)+ZQGM(JLON,JLEV)*YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP(JLON,JLEV)*ZINVG
    ENDDO
  ENDDO
ENDIF
! lightening density
IF (LFLASH) THEN
  IF (YDCPG_DIM%KSTEP==0) YDMF_PHYS%OUT%FLASH=0._JPRB

  CALL DIAGFLASH(YDCFU, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YDCPG_DIM%KSTEP, &
  & ZQCM, ZQIM, ZQRM, ZQSM, ZQGM, ZQHM, YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP, ZTM, YLMF_PHYS_STATE%YCPG_PHY%W,  &
  & YDMF_PHYS_SURF%GSD_VF%PLSM, YDMF_PHYS%OUT%FLASH)
ENDIF
!!! modif pour LMSE non activee
IF (LLMSE) THEN
  DO JLEV=1,NTSSG+1
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%FCS(JLON,JLEV)=-ZSFTH_(JLON)*ZRHODREFM__(JLON,IKB)*RCPD
      YDMF_PHYS%OUT%FCLL(JLON,JLEV) = YDMF_PHYS%OUT%FCLL(JLON,JLEV)*ZRHODREFM__(JLON,IKB)
      YDMF_PHYS%OUT%FCLN(JLON,JLEV) = YDMF_PHYS%OUT%FCLN(JLON,JLEV)*ZRHODREFM__(JLON,IKB)
      YDMF_PHYS%OUT%FEVL(JLON,JLEV) = YDMF_PHYS%OUT%FEVL(JLON,JLEV)*ZRHODREFM__(JLON,IKB)
      YDMF_PHYS%OUT%FEVN(JLON,JLEV) = YDMF_PHYS%OUT%FEVN(JLON,JLEV)*ZRHODREFM__(JLON,IKB)
    ENDDO
  ENDDO
ENDIF
IF (LSFORCS) THEN
  DO JLEV=1,NTSSG+1
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%FCS(JLON,JLEV)=-ZSFTH_(JLON)*ZRHODREFM__(JLON,IKB)*RCPD
      ZDELTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-ZTSURF(JLON)))
      YDMF_PHYS%OUT%FCLL(JLON,JLEV)=-ZSFRV_(JLON)*ZRHODREFM__(JLON,IKB)* FOLH (ZTSURF(JLON),0._JPRB)*(1.0_JPRB-ZDELTA)
      YDMF_PHYS%OUT%FCLN(JLON,JLEV)=-ZSFRV_(JLON)*ZRHODREFM__(JLON,IKB)* FOLH (ZTSURF(JLON),0._JPRB)*ZDELTA
    ENDDO
  ENDDO
ENDIF  

DO JSG  = 1, NTSSG+1
  DO JLEV = 0, YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA
      YDMF_PHYS%OUT%FRTH(JLON,JLEV,JSG)=YDMF_PHYS%OUT%FRTH(JLON,JLEV,JSG)+ZBUDTH_(JLON)
    ENDDO
  ENDDO
ENDDO

! daand: radflex
IF (LINTFLEX) THEN
  ! account for radiation separately
  LLRAD=.NOT.LRADFLEX
    
  CALL APL_AROME2INTFLEX(YGFL, YDPARAR, YDPHY, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, &
  & PDTPHY, YLMF_PHYS_STATE%YCPG_PHY%XYB%RDELP, YLMF_PHYS_STATE%U, YLMF_PHYS_STATE%V, YLMF_PHYS_STATE%T,           &
  & YDCPG_GPAR%VTS, YLMF_PHYS_STATE%YCPG_DYN%RCP%CP, ZFPR, LLRAD, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%FRSO,          &
  & YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDT, ZTENDR, ZTENDTKE, ZTENDEXT, YLPROCSET)
ENDIF


! Precipitation Type

! Compute wet-bulb temperature at 2 meters (suppose homogeneity of qv/ql/qi )
!ZPCLS(KIDIA:KFDIA)=PAPRSM(KIDIA:KFDIA,KLEV)-2._JPRB/ZZZF(KIDIA:KFDIA,1,KLEV)*&
!                 &(PAPRSM(KIDIA:KFDIA,KLEV)-PAPRSFM(KIDIA:KFDIA,KLEV))

CALL PPWETPOINT(YDPHY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YLMF_PHYS_STATE%YCPG_PHY%PRE(:, YDCPG_DIM%KFLEVG), &
& YDMF_PHYS%OUT%TCLS, YDMF_PHYS%OUT%QCLS, ZQCM(:, YDCPG_DIM%KFLEVG), ZQIM(:, YDCPG_DIM%KFLEVG), YDMF_PHYS%OUT%TPWCLS        &
& )

IF (LDPRECIPS.OR.LDPRECIPS2) THEN

  !initialisation de ZDZZ
  DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZDZZ(JLON,1)=ZAPHIM(JLON,0)*ZINVG-ZZZ_(JLON,1)
  ENDDO
  DO JLEV = 2, YDCPG_DIM%KFLEVG
    DO JLON = YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZDZZ(JLON,JLEV)=ZZZ_(JLON,JLEV+IKL)-ZZZ_(JLON,JLEV)
    ENDDO
  ENDDO


  ! Compute wet-bulb temperature
  DO JLEV=1,YDCPG_DIM%KFLEVG
      CALL PPWETPOINT(YDPHY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YLMF_PHYS_STATE%YCPG_PHY%PREF(:, JLEV), &
      & ZTM(:, JLEV), ZQVM(:, JLEV), ZQCM(:, JLEV), ZQIM(:, JLEV), ZTPW(:, JLEV))
  ENDDO

  IF (LDPRECIPS) THEN
   ! Defined precipitation type 
   !
   NDTPRECCUR=INT(MOD(ZSTATI/TSTEP,REAL(NDTPREC)))+1_JPIM
   !PDPRECIPS(:,NDTPRECCUR)=HUGE(1._JPRB)
   YDMF_PHYS_TMP%PRC%DPRECIPS(:,NDTPRECCUR)=0._JPRB

   !WRITE(NULOUT,*)'sous apl_arome NDTPRECCUR=',NDTPRECCUR,NDTPREC
   CALL DPRECIPS(YDPRECIPS, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YDOROG%OROG,                               &
   & YDMF_PHYS%OUT%TPWCLS, YDMF_PHYS%OUT%DIAGH, YLMF_PHYS_STATE%YCPG_DYN%PHIF, ZDZZ, ZTPW, ZQCM, YDMF_PHYS%OUT%FPLSL(:, YDCPG_DIM%KFLEVG), &
   & YDMF_PHYS%OUT%FPLSN(:, YDCPG_DIM%KFLEVG), YDMF_PHYS%OUT%FPLSG(:, YDCPG_DIM%KFLEVG), YDMF_PHYS_TMP%PRC%DPRECIPS(:, NDTPRECCUR)         &
   & )
  ENDIF

  IF (LDPRECIPS2) THEN

   !Idem for an other time step and an other period
   NDTPRECCUR2=INT(MOD(ZSTATI/TSTEP,REAL(NDTPREC2)))+1_JPIM
   YDMF_PHYS_TMP%PRC%DPRECIPS2(:,NDTPRECCUR2)=0._JPRB

   CALL DPRECIPS(YDPRECIPS, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, YDOROG%OROG,                               &
   & YDMF_PHYS%OUT%TPWCLS, YDMF_PHYS%OUT%DIAGH, YLMF_PHYS_STATE%YCPG_DYN%PHIF, ZDZZ, ZTPW, ZQCM, YDMF_PHYS%OUT%FPLSL(:, YDCPG_DIM%KFLEVG), &
   & YDMF_PHYS%OUT%FPLSN(:, YDCPG_DIM%KFLEVG), YDMF_PHYS%OUT%FPLSG(:, YDCPG_DIM%KFLEVG), YDMF_PHYS_TMP%PRC%DPRECIPS2(:, NDTPRECCUR2)       &
   & )

  ENDIF

ENDIF

!Save surface temperature
IF (LMSE.OR.LSFORCS) THEN
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
!      4.2  COMPUTE THE PHYS. TENDENCY FOR "T" AND "w"
!           ------------------------------------------

IF (LVERTFE.AND.LVFE_GWMPA) THEN
  ! * case LVFE_GWMPA not yet coded.
  !   (in this case ZGWT1 must be computed at full levels and
  !   not at half levels)
  CALL ABOR1(' APL_AROME: case LVFE_GWMPA not yet coded if LMPA=T!')
ENDIF

! * compute ZTT1:
IF (LSLAG.AND.LTWOTL) THEN
  DO JLEV=1,YDCPG_DIM%KFLEVG
    DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZTT1(JROF,JLEV)=YDVARS%T%T0(JROF,JLEV)+PDTPHY*ZTENDT(JROF,JLEV)
    ENDDO
  ENDDO
ELSE
  DO JLEV=1,YDCPG_DIM%KFLEVG
    DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZTT1(JROF,JLEV)=YDVARS%T%T9(JROF,JLEV)+PDTPHY*ZTENDT(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

! * compute ZGWT1 = tendency of gw:
IF (LNHDYN) THEN
  ! Valid for LVFE_GWMPA=F only; ZGWT1 assumed to be half level values.
  DO JLEV=1,YDCPG_DIM%KFLEVG-1
    DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZGWT1(JROF,JLEV)=0.5_JPRB*RG*(ZTENDW(JROF,JLEV)+ZTENDW(JROF,JLEV+1))
    ENDDO
  ENDDO
  DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZGWT1(JROF,YDCPG_DIM%KFLEVG)=0.0_JPRB
    ZGWT1(JROF,0)=0.0_JPRB
  ENDDO
ENDIF

! * convert gw tendency in d tendency:
IF(LNHDYN) THEN

  IF (LGWADV) THEN
    ZTENDD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=ZGWT1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  ELSE

    ! * Provide the appropriate version of (RT) at t+dt for GNHGW2SVDAROME:
    IF (L_RDRY_VD) THEN
      ! Use Rd because "dver" is currently defined with Rd.
      ZRTT1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=RD*ZTT1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
    ELSE
      ! Use "moist R" because "dver" is defined with "moist R".
      ! Unfortunately, R(t+dt) is not yet available there, use R(t) instead.
      ! "Moist R" tendency is neglected in the below call to GNHGW2SVDAROME.
      DO JLEV=1,YDCPG_DIM%KFLEVG
        DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZRTT1(JROF,JLEV)=YDCPG_DYN0%RCP%R(JROF,JLEV)*ZTT1(JROF,JLEV)
        ENDDO
      ENDDO
    ENDIF

    ! * Do conversion:
    IF (LSLAG.AND.LTWOTL) THEN
      CALL GNHGW2SVDAROME(YDGEOMETRY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_PHY0%PREHYDF, YDCPG_PHY0%XYB%LNPR, &
      & ZRTT1, YDCPG_PHY0%PREF, ZGWT1, ZTENDD)  
    ELSE
      CALL GNHGW2SVDAROME(YDGEOMETRY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_PHY9%PREHYDF, YDCPG_PHY9%XYB%LNPR, &
      & ZRTT1, YDCPG_PHY9%PREF, ZGWT1, ZTENDD)  
    ENDIF

  ENDIF
ELSE
  ZTENDD=0.0_JPRB
ENDIF

!      4.3  PUT THE TENDENCIES IN PB1/GFLT1/GMVT1.
!           --------------------------------------


IF (LSLAG) CALL CP_PTRSLB1(YDMODEL%YRML_DYN%YRDYN, YDPTRSLB1, ISLB1U9, ISLB1V9, ISLB1T9, ISLB1VD9, &
           & ISLB1GFL9)
IF ( LINTFLEX ) THEN

  ! Set GFL tendencies to 0
  ZTENDGFL(:,:,:) = 0.0_JPRB

  CALL CPTEND_FLEX( YDLDDH, YDMDDH, YGFL, YDPHY, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG,     &
  & YDGSGEOM%GNORDL, YDGSGEOM%GNORDM, YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,             &
  & YLMF_PHYS_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_STATE%U, YLMF_PHYS_STATE%V, YLMF_PHYS_STATE%T, YLMF_PHYS_STATE%YGSP_RR%T, &
  & PGFL, YLPROCSET, YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDH, ZTENDGFL, YDMF_PHYS%OUT%FHSCL,                    &
  & YDMF_PHYS%OUT%FHSCN, YDMF_PHYS%OUT%FHSSL, YDMF_PHYS%OUT%FHSSN, YDMF_PHYS%OUT%FHPCL, YDMF_PHYS%OUT%FHPCN,             &
  & YDMF_PHYS%OUT%FHPSL, YDMF_PHYS%OUT%FHPSN, PFEPFP =YDMF_PHYS%OUT%FEPFP, PFCMPCQ=YDMF_PHYS%OUT%FCMPCQ,                 &
  & PFCMPSN=YDMF_PHYS%OUT%FCMPSN, PFCMPSL=YDMF_PHYS%OUT%FCMPSL, YDDDH=YDDDH )
  
  CALL CPUTQY(YDGEOMETRY%YRDIMV, YDGMV, YGFL, YDPTRSLB1, YDPHY, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, &
  & YDCPG_DIM%KFLEVG, PDTPHY, IPGFL, ISLB1T9, ISLB1U9, ISLB1V9, ISLB1VD9, ISLB1GFL9, ZTENDH, ZTENDT,              &
  & YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDU, ZTENDV, ZTENDD, ZTENDGFL, YLMF_PHYS_STATE%YCPG_DYN%RCP%CP,  &
  & YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_STATE%T, YLMF_PHYS_STATE%U, YLMF_PHYS_STATE%V, PB1,              &
  & PGMVT1, PGFLT1, YDMF_PHYS%OUT%FDIS)    
   
ELSE

  ! start ZTENDGFLR at 1 because it is dimensionned (:,:,0:n)
  CALL CPUTQY_AROME(YDMODEL, YDGEOMETRY%YRDIMV, YDGMV, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA,          &
  & YDCPG_DIM%KFLEVG, PDTPHY, IPGFL, IPTR, ZTENDT, ZTENDGFLR(:, :, 1:), YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, &
  & ZTENDD, PB1, PGMVT1, PGFLT1)
ENDIF
  

!     ------------------------------------------------------------------ 
!*       5.    Final calculations.
!              -------------------

! * Restore the initial value of some pseudo-historical surface buffers
!   if relevant.
IF (LLDIAB) THEN
  IF (LL_SAVE_PHSURF) THEN
    CALL MF_PHYS_SAVE_PHSURF_PART2 (YDCPG_DIM, YDMF_PHYS_TMP%SAV%DDAL, YDMF_PHYS_TMP%SAV%DDOM, YDMF_PHYS_TMP%SAV%ENTCH, YDMF_PHYS_TMP%SAV%FHPS, &
                                  & YDMF_PHYS_TMP%SAV%GZ0F, YDMF_PHYS_TMP%SAV%GZ0HF, YDMF_PHYS_TMP%SAV%HV, YDMF_PHYS_TMP%SAV%PBLH, YDMF_PHYS_TMP%SAV%QSH, &
                                  & YDMF_PHYS_TMP%SAV%UDAL, YDMF_PHYS_TMP%SAV%UDGRO, YDMF_PHYS_TMP%SAV%UDOM, YDMF_PHYS_TMP%SAV%UNEBH, YDMF_PHYS_SURF, &
                                  & YDVARS, YDSURF, YDMODEL)
  ENDIF
ENDIF

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
  CALL WRITEPHYSIO(YDGEOMETRY, YDCPG_MISC, YDCPG_PHY0, YDCPG_PHY9, YDMF_PHYS, YDCPG_DYN0, YDCPG_DYN9,              &
  & YDMF_PHYS_SURF, YDVARS, YDSURF, YDDPHY, YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KFDIA, YDCPG_DIM%KIDIA,          &
  & YDCPG_DIM%KGL1, YDCPG_DIM%KGL2, YDCPG_DIM%KSTGLO, YDCPG_DIM%KSTEP, NTSSG, YSP_SBD%      NLEVS, YDGSGEOM%GELAM, &
  & YDGSGEOM%GEMU, YDGSGEOM%GM, YDOROG%      OROG, YDGSGEOM%RCORI, YDCSGEOM%RATATH, YDCSGEOM%RATATX,               &
  & YDGSGEOM%       GECLO, YDGSGEOM%GESLO, YDMF_PHYS_TMP%RDG%CVGQ, YDMF_PHYS_TMP%RDG%LCVQ, YDMF_PHYS_TMP%RDG%MU0)
ENDIF

IF (LEDR) THEN
  YDMF_PHYS_SURF%GSD_DI%PXEDR(:,:)=YDMF_PHYS%OUT%EDR(:,:)
ENDIF

CALL MF_PHYS_PRECIPS (YDCPG_DIM, YDMF_PHYS_TMP%PRC%DPRECIPS, YDMF_PHYS_TMP%PRC%DPRECIPS2, YDMF_PHYS_SURF, YDMODEL)

! Restore Tt and grad(Tt) for NHQE model.
IF (LNHQE) THEN
  CALL MF_PHYS_NHQE_PART2 (YDGEOMETRY, YDCPG_DIM, YDMF_PHYS_TMP%NHQ%TT0, YDMF_PHYS_TMP%NHQ%TT0L, YDMF_PHYS_TMP%NHQ%TT0M, YDMF_PHYS_TMP%NHQ%TT9, YDVARS)
ENDIF

!     ------------------------------------------------------------------

!       6. destructor for procset
IF (LINTFLEX) CALL CLEANINTPROCSET(YLPROCSET)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('APL_AROME', 1, ZHOOK_HANDLE)

CONTAINS

SUBROUTINE SWAP_THS
IF (LLSWAP_THS) THEN
  ZTHSIN_ => ZTHSAVE__(:,1:YDCPG_DIM%KFLEVG)
  ZTHS__ => ZTHSWAP__
ELSE
  ZTHSIN_ => ZTHSWAP__(:,1:YDCPG_DIM%KFLEVG)
  ZTHS__ => ZTHSAVE__
ENDIF
LLSWAP_THS=.NOT.LLSWAP_THS
END SUBROUTINE SWAP_THS

SUBROUTINE SWAP_RS
IF (LLSWAP_RS) THEN
  ZRSIN_ => ZRSAVE_
  ZRS_   => ZRSWAP_
ELSE
  ZRSIN_ => ZRSWAP_
  ZRS_   => ZRSAVE_
ENDIF
LLSWAP_RS=.NOT.LLSWAP_RS
END SUBROUTINE SWAP_RS

SUBROUTINE SWAP_SVS
IF (LLSWAP_SVS) THEN
  ZSVSIN_ => ZSVSAVE_
  ZSVS_   => ZSVSWAP_
ELSE
  ZSVSIN_ => ZSVSWAP_
  ZSVS_   => ZSVSAVE_
ENDIF
LLSWAP_SVS=.NOT.LLSWAP_SVS
END SUBROUTINE SWAP_SVS

SUBROUTINE SWAP_SVM
IF (LLSWAP_SVM) THEN
  ZSVMIN_ => ZSVMSAVE_
  ZSVM_   => ZSVMSWAP_
ELSE
  ZSVMIN_ => ZSVMSWAP_
  ZSVM_   => ZSVMSAVE_
ENDIF
LLSWAP_SVM=.NOT.LLSWAP_SVM
END SUBROUTINE SWAP_SVM

SUBROUTINE SWAP_LIMAS
IF (LLSWAP_LIMAS) THEN
  ZLIMASIN_ => ZLIMASAVE_
  ZLIMAS_   => ZLIMASWAP_
ELSE
  ZLIMASIN_ => ZLIMASWAP_
  ZLIMAS_   => ZLIMASAVE_
ENDIF
LLSWAP_LIMAS=.NOT.LLSWAP_LIMAS
END SUBROUTINE SWAP_LIMAS

END SUBROUTINE APL_AROME
