!OPTIONS XOPT(HSFUN)
#ifdef RS6K
@PROCESS HOT(NOVECTOR) NOSTRICT
#endif
SUBROUTINE CLOUDSC &
 !---input
 & (YDECLDP,YDECUMF,YDEPHLI,YDERAD, YDEPHY,KIDIA,    KFDIA,    KLON,    KLEV, KSPPN2D, &
 & PTSPHY, &
 & PT, PQ, &
 & TENDENCY_CML_T, TENDENCY_CML_Q, TENDENCY_CML_A, TENDENCY_CML_CLD, &
 & TENDENCY_TMP_T, TENDENCY_TMP_Q, TENDENCY_TMP_A, TENDENCY_TMP_CLD, &
 & TENDENCY_LOC_T, TENDENCY_LOC_Q, TENDENCY_LOC_A, TENDENCY_LOC_CLD, &
 & PVFA, PVFL, PVFI, PDYNA, PDYNL, PDYNI, &
 & PHRSW,    PHRLW, &
 & PVERVEL,  PAP,      PAPH, &
 & PLSM,     PGAW,     LDCUM,    KTYPE, &
 & PLU,      PLUDE,    PLUDELI,  PSNDE,    PMFU,     PMFD, PGP2DSPP, &
 & LDSLPHY, &
 !---prognostic fields
 & PA, &
 & PCLV,  &
 & PSUPSAT, &
!-- arrays for aerosol-cloud interactions
!!! & PQAER,    KAER, &
 & PLCRIT_AER,PICRIT_AER, &
 & PRE_ICE, &
 & PCCN,     PNICE, &
 !---diagnostic output
 & PCOVPTOT, PFSD, PRAINFRAC_TOPRFZ, &
 !---resulting fluxes
 & PFSQLF,   PFSQIF ,  PFCQNNG,  PFCQLNG, &
 & PFSQRF,   PFSQSF ,  PFCQRNG,  PFCQSNG, &
 & PFSQLTUR, PFSQITUR , &
 & PFPLSL,   PFPLSN,   PFHPSL,   PFHPSN, &
 & PEXTRA,   KFLDX)  

!===============================================================================
!**** *CLOUDSC* -  ROUTINE FOR PARAMETRIZATION OF CLOUD PROCESSES
!                  FOR PROGNOSTIC CLOUD SCHEME
!!
!     M.Tiedtke, C.Jakob, A.Tompkins, R.Forbes     (E.C.M.W.F.)
!!
!     PURPOSE
!     -------
!          THIS ROUTINE UPDATES THE CONV/STRAT CLOUD FIELDS.
!
!        1. Initial set up
!        2. Tidy up of input values
!        3. Subgrid sources/sinks (convection
!           3.1a Clipping of supersaturation
!           3.1b Addition of cloud from clipping at previous timestep
!           3.2 Detrainment of cloud water from convective updrafts
!           3.3 Vertical advection due to convective subsidence, and 
!               subsequent evaporation due to adiabatic warming 
!           3.4 Erosion at cloud edges by turbulent mixing of cloud air
!               with unsaturated environmental air
!           3.5 Evaporation/condensation of cloud water in connection
!               with heating/cooling such as by subsidence/ascent
!        4. Microphysical processes
!           4.1 Deposition onto ice when liquid water present (Bergeron-Findeison) 
!           4.2 Sedimentation of rain, snow and ice
!           4.3 Define subgrid precipiation fractions
!           4.4 Conversion of cloud ice to snow (aggregation)
!           4.5 Conversion of cloud water into rain (collision-coalescence)
!           4.6 Riming of snow
!           4.7 Melting of snow and ice
!           4.8 Freezing of rain
!           4.9 Freezing of cloud
!           4.10 Evaporation of rain
!           4.11 Sublimation of snow
!
!        5. Implicit solver for all processes
!
!        Note: Turbulent transports of s,q,u,v at cloud tops due to
!           buoyancy fluxes and lw radiative cooling are treated in 
!           the VDF scheme
!!
!     INTERFACE.
!     ----------
!     *CLOUDSC* IS CALLED FROM *CLOUD_LAYER* CALLED FROM *CALLPAR*
!     THE ROUTINE TAKES ITS INPUT FROM THE PROGNOSTIC VARIABLES:
!     T,Q,L,I,A AND DETRAINMENT OF CLOUD WATER FROM THE
!     CONVECTIVE CLOUDS (MASSFLUX CONVECTION SCHEME), BOUNDARY
!     LAYER TURBULENT FLUXES OF HEAT AND MOISTURE, RADIATIVE FLUXES,
!     OMEGA.
!     IT RETURNS ITS OUTPUT TO:
!      1.MODIFIED TENDENCIES OF MODEL VARIABLES T AND Q
!        AS WELL AS CLOUD LIQUID, ICE, RAIN, SNOW AND CLOUD FRACTION
!      2.GENERATES PRECIPITATION FLUXES FROM GRID-SCALE CLOUDS
!!
!     EXTERNALS.
!     ----------
!          NONE
!!
!     MODIFICATIONS.
!     -------------
!      M. TIEDTKE    E.C.M.W.F.     8/1988, 2/1990
!     CH. JAKOB      E.C.M.W.F.     2/1994 IMPLEMENTATION INTO IFS
!     A.TOMPKINS     E.C.M.W.F.     2002   NEW NUMERICS
!        01-05-22 : D.Salmond   Safety modifications
!        02-05-29 : D.Salmond   Optimisation
!        03-01-13 : J.Hague     MASS Vector Functions  J.Hague
!        03-10-01 : M.Hamrud    Cleaning
!        04-12-14 : A.Tompkins  New implicit solver and physics changes
!        04-12-03 : A.Tompkins & M.Ko"hler  moist PBL
!     G.Mozdzynski  09-Jan-2006  EXP security fix
!        19-01-09 : P.Bechtold  Changed increased RCLDIFF value for KTYPE=2
!        07-07-10 : A.Tompkins/R.Forbes  4-Phase flexible microphysics
!        01-03-11 : R.Forbes    Mixed phase changes and tidy up
!        01-10-11 : R.Forbes    Melt ice to rain, allow rain to freeze
!        01-10-11 : R.Forbes    Limit supersat to avoid excessive values
!        31-10-11 : M.Ahlgrimm  Add rain, snow and PEXTRA to DDH output
!        17-02-12 : F.Vana      Simplified/optimized LU factorization
!        18-05-12 : F.Vana      Cleaning + better support of sequential physics
!        N.Semane+P.Bechtold     04-10-2012 Add RVRFACTOR factor for small planet
!        01-02-13 : R.Forbes    New params of autoconv/acc,rain evap,snow riming
!        15-03-13 : F. Vana     New dataflow + more tendencies from the first call
!        K. Yessad (July 2014): Move some variables.
!        F. Vana  05-Mar-2015  Support for single precision
!        15-01-15 : R.Forbes    Added new options for snow evap & ice deposition
!        10-01-15 : R.Forbes    New physics for rain freezing
!        23-10-14 : P. Bechtold remove zeroing of convection arrays
!        15-12-2015 : R. Forbes Added inhomog option and variable rain fallspeed
!        27-01-2016 : M. Leutbecher & S.-J. Lock  Introduced SPP scheme (LSPP)
!        01-10-2016 : R. Forbes Tidy up routine
!        01-04-2017 : R. Forbes Modified numerics for rain/snow, 
!                               removed threshold for autoconv/accretion, 
!                               new turbulent erosion, ice evap, snow deposition
!        Oct-2017   : S.-J. Lock  Enabled options for new SPP microphysics perturbations
!        2017-11-11 M Ahlgrimm extent cloud heterogeneity for ice FSD
!        2019-01 : R.Forbes Passed in separate T/Q/A/CLD arrays for improved portability
!        2019-01 : R.Forbes Added additional diagnostics for cloud budget 
!        2019-01 : R.Forbes Added new parameters that can be set in the namelist
!        2019-01 : R.Forbes Tidy up and corrections to cloud budget diagnostics
!
!     REFERENCES.
!     ----------
!     Tietdke MWR 1993 - original description of the cloud parametrization
!     Jakob PhD 2000
!     Gregory et al. (2000) QJRMS
!     Tompkins el al. (2007) QJRMS - ice supersaturation parametrization
!     Forbes and Tompkins (2011) ECMWF Newsletter 129 - new prognostic liq/ice/rain/snow
!     Forbes et al. (2011) ECMWF Tech Memo 649 - new prognostic liq/ice/rain/snow
!     Forbes et al. (2014) ECMWF Newsletter 141 - freezing rain
!     Ahlgrimm and Forbes (2014) MWR - warm rain processes
!     Forbes and Ahlgrimm (2014) MWR - mixed-phase cloud processes
!     Ahlgrimm and Forbes (2015) MWR - subgrid heterogeneity of cloud condensate
!!
!===============================================================================

USE YOECLDP  , ONLY : TECLDP, NCLDQV, NCLDQL, NCLDQR, NCLDQI, NCLDQS, NCLV
USE YOEPHLI  , ONLY : TEPHLI
USE YOERAD   , ONLY : TERAD
USE YOEPHY   , ONLY : TEPHY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMMP0   , ONLY : LSCMEC
USE YOMCST   , ONLY : RG, RD, RCPD, RETV, RLVTT, RLSTT, RLMLT, RTT, RV, RA, RPI  , YRCST
USE YOETHF   , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 &                    R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE, RTICECU, &
 &                    RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2, YRTHF
USE YOECUMF  , ONLY : TECUMF
USE YOMPHYDER, ONLY : STATE_TYPE
USE SPP_MOD  , ONLY : YSPP_CONFIG, YSPP

IMPLICIT NONE

!-------------------------------------------------------------------------------
!                 Declare input/output arguments
!-------------------------------------------------------------------------------
 
TYPE(TECLDP)      ,INTENT(INOUT) :: YDECLDP
TYPE(TECUMF)      ,INTENT(INOUT) :: YDECUMF
TYPE(TEPHLI)      ,INTENT(INOUT) :: YDEPHLI
TYPE(TERAD)       ,INTENT(INOUT) :: YDERAD
TYPE(TEPHY)       ,INTENT(INOUT) :: YDEPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON             ! Number of grid points
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV             ! Number of levels
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPPN2D          ! Number of 2D patterns in SPP scheme
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY           ! Physics timestep
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)    ! T at start of callpar
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV)    ! Q at start of callpar
REAL(KIND=JPRB)   ,INTENT(IN)    :: TENDENCY_CML_T(KLON,KLEV)   ! T cumulative tendency
REAL(KIND=JPRB)   ,INTENT(IN)    :: TENDENCY_CML_Q(KLON,KLEV)   ! Q cumulative tendency
REAL(KIND=JPRB)   ,INTENT(IN)    :: TENDENCY_CML_A(KLON,KLEV)   ! A cumulative tendency
REAL(KIND=JPRB)   ,INTENT(IN)    :: TENDENCY_CML_CLD(KLON,KLEV,NCLV) ! CLD cumulative tendency
REAL(KIND=JPRB)   ,INTENT(OUT)   :: TENDENCY_LOC_T(KLON,KLEV)   ! T local output tendency
REAL(KIND=JPRB)   ,INTENT(OUT)   :: TENDENCY_LOC_Q(KLON,KLEV)   ! Q local output tendency
REAL(KIND=JPRB)   ,INTENT(OUT)   :: TENDENCY_LOC_A(KLON,KLEV)   ! A local output tendency
REAL(KIND=JPRB)   ,INTENT(OUT)   :: TENDENCY_LOC_CLD(KLON,KLEV,NCLV) ! CLD local output tendency
REAL(KIND=JPRB)   ,INTENT(IN)    :: TENDENCY_TMP_T(KLON,KLEV)   ! T tmp output tendency
REAL(KIND=JPRB)   ,INTENT(IN)    :: TENDENCY_TMP_Q(KLON,KLEV)   ! Q tmp output tendency
REAL(KIND=JPRB)   ,INTENT(IN)    :: TENDENCY_TMP_A(KLON,KLEV)   ! A tmp output tendency
REAL(KIND=JPRB)   ,INTENT(IN)    :: TENDENCY_TMP_CLD(KLON,KLEV,NCLV) ! CLD tmp output tendency
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFA(KLON,KLEV)  ! CC from VDF scheme
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFL(KLON,KLEV)  ! Liq from VDF scheme
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFI(KLON,KLEV)  ! Ice from VDF scheme
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNA(KLON,KLEV) ! CC from Dynamics
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNL(KLON,KLEV) ! Liq from Dynamics
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNI(KLON,KLEV) ! Liq from Dynamics
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRSW(KLON,KLEV) ! Short-wave heating rate
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRLW(KLON,KLEV) ! Long-wave heating rate
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) ! Vertical velocity
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)   ! Pressure on full levels
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1)! Pressure on half levels
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON)       ! Land fraction (0-1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAW(KLON)       ! Grid area=PGAW*4*RPI*RA**2
! From convection parametrization 
LOGICAL           ,INTENT(IN)    :: LDCUM(KLON)      ! Convection active
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYPE(KLON)      ! Convection type 0,1,2
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLU(KLON,KLEV)   ! Conv. condensate
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLUDE(KLON,KLEV) ! Conv. detrained water 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLUDELI(KLON,KLEV,2)! Conv. detrained droplets/ice
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNDE(KLON,KLEV,2)! Conv. detrained snow/rain
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV)  ! Conv. mass flux up
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFD(KLON,KLEV)  ! Conv. mass flux down
! Options and cloud prognostic variables
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP2DSPP(KLON,KSPPN2D) ! perturbation pattern
LOGICAL           ,INTENT(IN)    :: LDSLPHY          ! True if semi-lag physics
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV)    ! Original Cloud fraction (t)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLV(KLON,KLEV,NCLV) ! Cloud/precip prognostics
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSUPSAT(KLON,KLEV) ! Supersat clipped at previous time level in SLTEND
! Prognostic aerosol
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLCRIT_AER(KLON,KLEV) ! critical liquid mmr for rain autoconversion process
REAL(KIND=JPRB)   ,INTENT(IN)    :: PICRIT_AER(KLON,KLEV) ! critical liquid mmr for snow autoconversion process
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE_ICE(KLON,KLEV)    ! ice effective radius
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCCN(KLON,KLEV)       ! liquid cloud condensation nuclei
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNICE(KLON,KLEV)      ! ice number concentration (cf. CCN) 
! Precipitation related
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOVPTOT(KLON,KLEV)   ! Precip fraction
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAINFRAC_TOPRFZ(KLON)! Rain/snow fraction at top of refreezing layer 
! Flux diagnostics for DDH budget
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQLF(KLON,KLEV+1)  ! Flux of liquid
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQIF(KLON,KLEV+1)  ! Flux of ice
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQLNG(KLON,KLEV+1) ! -ve corr for liq
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQNNG(KLON,KLEV+1) ! -ve corr for ice
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQRF(KLON,KLEV+1)  ! Flux diagnostics
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQSF(KLON,KLEV+1)  !    for DDH, generic
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQRNG(KLON,KLEV+1) ! rain
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQSNG(KLON,KLEV+1) ! snow
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQLTUR(KLON,KLEV+1) ! liquid flux due to VDF
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSQITUR(KLON,KLEV+1) ! ice flux due to VDF
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLSL(KLON,KLEV+1) ! liq+rain sedim flux
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLSN(KLON,KLEV+1) ! ice+snow sedim flux
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPSL(KLON,KLEV+1) ! Enthalpy flux for liq
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPSN(KLON,KLEV+1) ! Enthalp flux for ice
! Extra fields for diagnostics
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTRA(KLON,KLEV,KFLDX) ! extra fields
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX ! Number of extra fields
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFSD(KLON,KLEV) ! cloud condensate fractional standard deviation

!-------------------------------------------------------------------------------
!                       Declare local variables
!-------------------------------------------------------------------------------

!  condensation and evaporation terms
REAL(KIND=JPRB) :: ZLCOND1(KLON), ZLCOND2(KLON)
REAL(KIND=JPRB) :: ZLEVAP,        ZLEROS
REAL(KIND=JPRB) :: ZLEVAPL(KLON), ZLEVAPI(KLON)
! autoconversion terms
REAL(KIND=JPRB) :: ZRAINAUT(KLON), ZSNOWAUT(KLON)
REAL(KIND=JPRB) :: ZLIQCLD(KLON),  ZICECLD(KLON)

REAL(KIND=JPRB) :: ZFOKOOP(KLON), ZFOEALFA(KLON,KLEV+1)
REAL(KIND=JPRB) :: ZICENUCLEI(KLON) ! number concentration of ice nuclei
REAL(KIND=JPRB) :: ZLICLD(KLON)
REAL(KIND=JPRB) :: ZACOND
REAL(KIND=JPRB) :: ZAEROS
REAL(KIND=JPRB) :: ZLFINALSUM(KLON)
REAL(KIND=JPRB) :: ZDQS(KLON)
REAL(KIND=JPRB) :: ZTOLD(KLON)
REAL(KIND=JPRB) :: ZQOLD(KLON)  
REAL(KIND=JPRB) :: ZDTGDP(KLON) 
REAL(KIND=JPRB) :: ZRDTGDP(KLON)  
REAL(KIND=JPRB) :: ZTRPAUS(KLON)
REAL(KIND=JPRB) :: ZCOVPCLR(KLON)   
REAL(KIND=JPRB) :: ZPRECLR
REAL(KIND=JPRB) :: ZCOVPTOT(KLON)    
REAL(KIND=JPRB) :: ZCOVPMAX(KLON)
REAL(KIND=JPRB) :: ZQPRETOT(KLON)
REAL(KIND=JPRB) :: ZDPEVAP
REAL(KIND=JPRB) :: ZDTFORC
REAL(KIND=JPRB) :: ZDTDIAB
REAL(KIND=JPRB) :: ZTP1(KLON,KLEV)   
REAL(KIND=JPRB) :: ZLDEFR(KLON)
REAL(KIND=JPRB) :: ZLDIFDT(KLON)
REAL(KIND=JPRB) :: ZDTGDPF(KLON)
REAL(KIND=JPRB) :: ZLCUST(KLON,NCLV)
REAL(KIND=JPRB) :: ZACUST(KLON)
REAL(KIND=JPRB) :: ZMF(KLON) 

REAL(KIND=JPRB) :: ZRHO(KLON)
REAL(KIND=JPRB) :: ZGDP(KLON)

! Accumulators of A,B,and C factors for cloud equations
REAL(KIND=JPRB) :: ZSOLAB(KLON) ! -ve implicit CC
REAL(KIND=JPRB) :: ZSOLAC(KLON) ! linear CC
REAL(KIND=JPRB) :: ZANEW
REAL(KIND=JPRB) :: ZANEWM1(KLON) 

REAL(KIND=JPRB) :: ZDA(KLON)
REAL(KIND=JPRB) :: ZLI(KLON,KLEV)
REAL(KIND=JPRB) :: ZA(KLON,KLEV)
REAL(KIND=JPRB) :: ZAORIG(KLON,KLEV) ! start of scheme value for CC

LOGICAL :: LLFLAG(KLON)
LOGICAL :: LLO1
LOGICAL :: LLRAMID, LLRCLDIFF, LLRCLCRIT, LLRLCRITSNOW
LOGICAL :: LLRAINEVAP, LLSNOWSUBLIM, LLQSATVERVEL

INTEGER(KIND=JPIM) :: ICALL, IK, JK, JL, JM, JN, JO, IS
INTEGER(KIND=JPIM) :: IPRAMID, IPRCLDIFF, IPRCLCRIT, IPRLCRITSNOW
INTEGER(KIND=JPIM) :: IPRAINEVAP, IPSNOWSUBLIM, IPQSATVERVEL

REAL(KIND=JPRB) :: ZDP(KLON), ZPAPHD(KLON)

REAL(KIND=JPRB) :: ZALFA
REAL(KIND=JPRB) :: ZALFAW
REAL(KIND=JPRB) :: ZBETA,ZBETA1
!REAL(KIND=JPRB) :: ZBOTT
REAL(KIND=JPRB) :: ZCFPR
REAL(KIND=JPRB) :: ZCOR
REAL(KIND=JPRB) :: ZCDMAX
REAL(KIND=JPRB) :: ZLCONDLIM
REAL(KIND=JPRB) :: ZDENOM
REAL(KIND=JPRB) :: ZDPMXDT
REAL(KIND=JPRB) :: ZDPR
REAL(KIND=JPRB) :: ZDTDP
REAL(KIND=JPRB) :: ZE
REAL(KIND=JPRB) :: ZEPSEC
REAL(KIND=JPRB) :: ZFAC, ZFACI, ZFACW
REAL(KIND=JPRB) :: ZGDCP
REAL(KIND=JPRB) :: ZINEW
REAL(KIND=JPRB) :: ZLCRIT
REAL(KIND=JPRB) :: ZMFDN
REAL(KIND=JPRB) :: ZPRECIP
REAL(KIND=JPRB) :: ZQE
REAL(KIND=JPRB) :: ZQSAT, ZQTMST, ZRDCP
REAL(KIND=JPRB) :: ZRHC, ZSIG, ZSIGK
REAL(KIND=JPRB) :: ZWTOT
REAL(KIND=JPRB) :: ZZCO, ZZDL, ZZRH, ZZZDT, ZQADJ
REAL(KIND=JPRB) :: ZQNEW, ZTNEW
REAL(KIND=JPRB) :: ZRG_R,ZGDPH_R,ZCONS1,ZCOND,ZCONS1A
REAL(KIND=JPRB) :: ZLFINAL
REAL(KIND=JPRB) :: ZMELT
REAL(KIND=JPRB) :: ZEVAP
REAL(KIND=JPRB) :: ZFRZ
REAL(KIND=JPRB) :: ZVPLIQ, ZVPICE
REAL(KIND=JPRB) :: ZADD, ZBDD, ZCVDS, ZICE0, ZDEPOS
REAL(KIND=JPRB) :: ZSUPSAT(KLON)
REAL(KIND=JPRB) :: ZFALL
REAL(KIND=JPRB) :: ZRE_ICE
REAL(KIND=JPRB) :: ZRLDCP
REAL(KIND=JPRB) :: ZQP1ENV
REAL(KIND=JPRB) :: ZDZ
REAL(KIND=JPRB) :: ZMU_RAMID,   ZXRAMID
REAL(KIND=JPRB) :: ZMU_RCLDIFF
REAL(KIND=JPRB) :: ZMU_RCLCRITL, ZMU_RCLCRITS
REAL(KIND=JPRB) :: ZMU_RLCRITSNOW
REAL(KIND=JPRB) :: ZMU_RAINEVAP
REAL(KIND=JPRB) :: ZMU_SNOWSUBLIM
REAL(KIND=JPRB) :: ZMU_QSATVERVEL, ZPVERVEL(KLON,KLEV)

INTEGER(KIND=JPIM) :: IPHASE(NCLV) ! marker for water phase of each species
                                   ! 0=vapour, 1=liquid, 2=ice

INTEGER(KIND=JPIM) :: IMELT(NCLV)  ! marks melting linkage for ice categories
                                   ! ice->liquid, snow->rain

LOGICAL :: LLFALL(NCLV)      ! marks falling species
                             ! LLFALL=0, cloud cover must > 0 for zqx > 0
                             ! LLFALL=1, no cloud needed, zqx can evaporate

REAL(KIND=JPRB) :: ZLIQFRAC(KLON,KLEV)  ! cloud liquid water fraction: ql/(ql+qi)
REAL(KIND=JPRB) :: ZICEFRAC(KLON,KLEV)  ! cloud ice water fraction: qi/(ql+qi)
REAL(KIND=JPRB) :: ZQX(KLON,KLEV,NCLV)  ! water variables
REAL(KIND=JPRB) :: ZQX0(KLON,KLEV,NCLV) ! water variables at start of scheme
REAL(KIND=JPRB) :: ZQXN(KLON,NCLV)      ! new values for zqx at time+1
REAL(KIND=JPRB) :: ZQXFG(KLON,NCLV)     ! first guess values including precip
REAL(KIND=JPRB) :: ZQXNM1(KLON,NCLV)    ! new values for zqx at time+1 at level above
REAL(KIND=JPRB) :: ZFLUXQ(KLON,NCLV)    ! fluxes convergence of species (needed?)

REAL(KIND=JPRB) :: ZPFPLSX(KLON,KLEV+1,NCLV) ! generalized precipitation flux
REAL(KIND=JPRB) :: ZLNEG(KLON,KLEV,NCLV)     ! for negative correction diagnostics
REAL(KIND=JPRB) :: ZMELTMAX(KLON)
REAL(KIND=JPRB) :: ZFRZMAX(KLON)
REAL(KIND=JPRB) :: ZICETOT(KLON)

REAL(KIND=JPRB) :: ZQXN2D(KLON,KLEV,NCLV)   ! water variables store

REAL(KIND=JPRB) :: ZQSMIX(KLON,KLEV) ! diagnostic mixed phase saturation 
REAL(KIND=JPRB) :: ZQSLIQ(KLON,KLEV) ! liquid water saturation
REAL(KIND=JPRB) :: ZQSICE(KLON,KLEV) ! ice water saturation

!REAL(KIND=JPRB) :: ZRHM(KLON,KLEV) ! diagnostic mixed phase RH
!REAL(KIND=JPRB) :: ZRHL(KLON,KLEV) ! RH wrt liq
!REAL(KIND=JPRB) :: ZRHI(KLON,KLEV) ! RH wrt ice

REAL(KIND=JPRB) :: ZFOEEWMT(KLON,KLEV)
REAL(KIND=JPRB) :: ZFOEEW(KLON,KLEV)
REAL(KIND=JPRB) :: ZFOEELIQT(KLON,KLEV)

REAL(KIND=JPRB) :: ZDQSLIQDT(KLON), ZDQSICEDT(KLON), ZDQSMIXDT(KLON)
REAL(KIND=JPRB) :: ZCORQSLIQ(KLON)
REAL(KIND=JPRB) :: ZCORQSICE(KLON) 
REAL(KIND=JPRB) :: ZCORQSMIX(KLON)
REAL(KIND=JPRB) :: ZEVAPLIMLIQ(KLON), ZEVAPLIMICE(KLON), ZEVAPLIMMIX(KLON)
!-------------------------------------------------------
! SOURCE/SINK array for implicit and explicit terms
!-------------------------------------------------------
! a POSITIVE value entered into the arrays is a...
!            Source of this variable
!            |
!            |   Sink of this variable
!            |   |
!            V   V
! ZSOLQA(JL,IQa,IQb)  = explicit terms
! ZSOLQB(JL,IQa,IQb)  = implicit terms
! Thus if ZSOLAB(JL,NCLDQL,IQV)=K where K>0 then this is 
! a source of NCLDQL and a sink of IQV
! put 'external' source terms such as PLUDE from 
! detrainment into explicit source/sink array diagnognal
! ZSOLQA(NCLDQL,NCLDQL)= -PLUDE
! i.e. a positive value is a sink! 
!-------------------------------------------------------

REAL(KIND=JPRB) :: ZSOLQA(KLON,NCLV,NCLV) ! explicit sources and sinks
REAL(KIND=JPRB) :: ZSOLQB(KLON,NCLV,NCLV) ! implicit sources and sinks
                        ! e.g. microphysical pathways between ice variables.
REAL(KIND=JPRB) :: ZQLHS(KLON,NCLV,NCLV)  ! n x n matrix storing the LHS of implicit solver
REAL(KIND=JPRB) :: ZVQX(NCLV)        ! fall speeds of three categories
REAL(KIND=JPRB) :: ZEXPLICIT

! REAL(KIND=JPRB) :: ZSINKSUM(KLON,NCLV)

! for sedimentation source/sink terms
REAL(KIND=JPRB) :: ZFALLSINK(KLON,NCLV)
REAL(KIND=JPRB) :: ZFALLSRCE(KLON,NCLV)

! for convection detrainment source and subsidence source/sink terms
REAL(KIND=JPRB) :: ZCONVSRCE(KLON,NCLV)
REAL(KIND=JPRB) :: ZCONVSINK(KLON,NCLV)
REAL(KIND=JPRB) :: ZADVW(KLON), ZADVWD(KLON)

! for supersaturation source term from previous timestep
REAL(KIND=JPRB) :: ZPSUPSATSRCE(KLON,NCLV)

! Numerical fit to wet bulb temperature
REAL(KIND=JPRB),PARAMETER :: ZTW1 = 1329.31_JPRB
REAL(KIND=JPRB),PARAMETER :: ZTW2 = 0.0074615_JPRB
REAL(KIND=JPRB),PARAMETER :: ZTW3 = 0.85E5_JPRB
REAL(KIND=JPRB),PARAMETER :: ZTW4 = 40.637_JPRB
REAL(KIND=JPRB),PARAMETER :: ZTW5 = 275.0_JPRB

REAL(KIND=JPRB) :: ZSUBSAT  ! Subsaturation for snow melting term         
REAL(KIND=JPRB) :: ZTDMTW0  ! Diff between dry-bulb temperature and 
                            ! temperature when wet-bulb = 0degC 

! Variables for deposition term
REAL(KIND=JPRB) :: ZTCG ! Temperature dependent function for ice PSD
REAL(KIND=JPRB) :: ZFACX1I, ZFACX1S! PSD correction factor
REAL(KIND=JPRB) :: ZAPLUSB,ZCORRFAC,ZCORRFAC2,ZPR02,ZTERM1,ZTERM2 ! for ice dep
REAL(KIND=JPRB) :: ZCLDTOPDIST(KLON) ! Distance from cloud top
REAL(KIND=JPRB) :: ZINFACTOR         ! No. of ice nuclei factor for deposition
REAL(KIND=JPRB) :: ZOVERLAP_LIQICE   ! Overlap fraction between SLW and ice

! Option control variables
INTEGER(KIND=JPIM) :: IWARMRAIN
INTEGER(KIND=JPIM) :: IRAINACC
INTEGER(KIND=JPIM) :: IEVAPSNOW
INTEGER(KIND=JPIM) :: IEVAPICE
INTEGER(KIND=JPIM) :: IDEPICE
INTEGER(KIND=JPIM) :: IVARFALL
INTEGER(KIND=JPIM) :: ITURBEROSION

! Autoconversion/accretion/riming/evaporation
REAL(KIND=JPRB) :: ZRAINACC(KLON)
REAL(KIND=JPRB) :: ZRAINCLD(KLON)
REAL(KIND=JPRB) :: ZRAINCLDM1(KLON)
REAL(KIND=JPRB) :: ZSNOWRIME(KLON)
REAL(KIND=JPRB) :: ZSNOWCLD(KLON)
REAL(KIND=JPRB) :: ZSNOWCLDM1(KLON)
REAL(KIND=JPRB) :: ZESATLIQ
REAL(KIND=JPRB) :: ZFALLCORR
REAL(KIND=JPRB) :: ZLAMBDA
REAL(KIND=JPRB) :: ZEVAP_DENOM
REAL(KIND=JPRB) :: ZCORR2
REAL(KIND=JPRB) :: ZKA
REAL(KIND=JPRB) :: ZCONST
REAL(KIND=JPRB) :: ZTEMP

! Cloud and precipitation inhomogeneity
REAL(KIND=JPRB) :: ZGRIDLEN, ZCLDRAINCORR
REAL(KIND=JPRB) :: ZPHIC, ZFRACSDC, ZEAUT
REAL(KIND=JPRB) :: ZPHIR, ZFRACSDR, ZEACC
REAL(KIND=JPRB) :: ZPHIP1, ZPHIP2, ZPHIP3, ZQTOT

! Rain freezing
LOGICAL :: LLRAINLIQ(KLON)  ! True if majority of raindrops are liquid (no ice core)

! SCM budget statistics 
REAL(KIND=JPRB), ALLOCATABLE :: ZSUMQ0(:,:),  ZSUMQ1(:,:) , ZERRORQ(:,:), &
                               &ZSUMH0(:,:),  ZSUMH1(:,:) , ZERRORH(:,:)
REAL(KIND=JPRB) :: ZRAIN

! Cloud budget
REAL(KIND=JPRB) :: ZBUDCC(KLON,11) ! cloud fraction budget array
REAL(KIND=JPRB) :: ZBUDL(KLON,21)  ! cloud liquid budget array
REAL(KIND=JPRB) :: ZBUDI(KLON,17)  ! cloud ice budget array

! Miscellaneous
REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZTMPL,ZTMPI,ZTMPA
REAL(KIND=JPRB) :: ZMM,ZRR

! Cloud heterogeneity variables
! 'FSD' is the fractional standard deviation= standard deviation/mean
! of cloud condensate
REAL(KIND=JPRB) :: ZDELZ              ! Layer thickness in km
REAL(KIND=JPRB) :: ZPHI               ! Factor to account for layer thickness and total water in ice FSD calculation
REAL(KIND=JPRB) :: ZRATFSD(KLON,KLEV) ! detrainment ratio: detrained mass from convection scheme
                                      ! divided by total cloud mass liq/ice in the grid box
REAL(KIND=JPRB) :: ZLIQF(KLON,KLEV)   ! Final cloud liquid at end of cloudsc
REAL(KIND=JPRB) :: ZICEF(KLON,KLEV)   ! Final cloud ice at end of cloudsc
REAL(KIND=JPRB) :: ZLFSD(KLON,KLEV)   ! Liquid condensate FSD
REAL(KIND=JPRB) :: ZIFSD(KLON,KLEV)   ! Ice condensate FSD
REAL(KIND=JPRB) :: ZFSD(KLON,KLEV)    ! Merged FSD for liquid and ice, passed to radiation scheme
REAL(KIND=JPRB) :: ZANEWP(KLON,KLEV)  ! Modified cloud fraction at end of time step, taking into account precip fraction, for ice FSD
REAL(KIND=JPRB) :: ZANEW2(KLON,KLEV)  ! Modified cloud fraction for liquid FSD (avoiding very small cloud fractions)
REAL(KIND=JPRB) :: ZIFSDBACK          ! Background ice FSD dependent on model grid scale

REAL(KIND=JPRB), PARAMETER :: ZR12=1.3_JPRB ! Factor accounting for parameterized along-track variability being an underestimate of the true area variability. 1.3 taken from Hill et al. 2015
!updated Oct 27th 2016
!REAL(KIND=JPRB), PARAMETER :: PPu1=0.517863_JPRB !Parameters used in background ice FSD calculation
!REAL(KIND=JPRB), PARAMETER :: ZU2=0.519894_JPRB 
!REAL(KIND=JPRB), PARAMETER :: PPu3=0.416821_JPRB
!REAL(KIND=JPRB), PARAMETER :: ZU4=-0.723216_JPRB

!updated Sep 28th 2017
REAL(KIND=JPRB), PARAMETER :: ZU1=0.446585_JPRB !Parameters used in background ice FSD calculation
REAL(KIND=JPRB), PARAMETER :: ZU2=0.308061_JPRB 
REAL(KIND=JPRB), PARAMETER :: ZU3=0.395736_JPRB
REAL(KIND=JPRB), PARAMETER :: ZU4=-0.744527_JPRB

REAL(KIND=JPRB) :: ZPAP(KLON)

LOGICAL            :: LLINDEX3(KLON,NCLV) ! index variable
INTEGER(KIND=JPIM) :: IORDER(KLON,NCLV), IORDV(NCLV) ! arrays for sorting explicit terms
REAL(KIND=JPRB) :: ZSINKSUM(KLON) 
REAL(KIND=JPRB) :: ZRATIO(KLON,NCLV), ZZRATIO, ZRAT, ZMAX
REAL(KIND=JPRB) :: ZEPSILON

#include "abor1.intfb.h"
#include "cuadjtq.intfb.h"

#include "fcttre.func.h"
#include "fccld.func.h"

!===============================================================================
IF (LHOOK) CALL DR_HOOK('CLOUDSC',0,ZHOOK_HANDLE)

ASSOCIATE(LAERICEAUTO=>YDECLDP%LAERICEAUTO, LAERICESED=>YDECLDP%LAERICESED, &
 & LAERLIQAUTOLSP=>YDECLDP%LAERLIQAUTOLSP, LAERLIQCOLL=>YDECLDP%LAERLIQCOLL, &
 & LCLDBUDC=>YDECLDP%LCLDBUDC, LCLDBUDL=>YDECLDP%LCLDBUDL, &
 & LCLDBUDI=>YDECLDP%LCLDBUDI, LCLDBUDT=>YDECLDP%LCLDBUDT, &
 & LCLDBUD_VERTINT=>YDECLDP%LCLDBUD_VERTINT, &
 & LCLDBUD_TIMEINT=>YDECLDP%LCLDBUD_TIMEINT, &
 & LCLDBUDGET=>YDECLDP%LCLDBUDGET, NCLDTOP=>YDECLDP%NCLDTOP, &
 & NSSOPT=>YDECLDP%NSSOPT, RAMID=>YDECLDP%RAMID, RAMIN=>YDECLDP%RAMIN, &
 & RCCN=>YDECLDP%RCCN, RCLCRIT_LAND=>YDECLDP%RCLCRIT_LAND, &
 & RCLCRIT_SEA=>YDECLDP%RCLCRIT_SEA, RCLDIFF=>YDECLDP%RCLDIFF, &
 & RCLDIFF_CONVI=>YDECLDP%RCLDIFF_CONVI, RCLDTOPCF=>YDECLDP%RCLDTOPCF, &
 & RCL_APB1=>YDECLDP%RCL_APB1, RCL_APB2=>YDECLDP%RCL_APB2, &
 & RCL_APB3=>YDECLDP%RCL_APB3, RCL_CDENOM1=>YDECLDP%RCL_CDENOM1, &
 & RCL_CDENOM2=>YDECLDP%RCL_CDENOM2, RCL_CDENOM3=>YDECLDP%RCL_CDENOM3, &
 & RCL_CONST1I=>YDECLDP%RCL_CONST1I, RCL_CONST1R=>YDECLDP%RCL_CONST1R, &
 & RCL_CONST1S=>YDECLDP%RCL_CONST1S, RCL_CONST2I=>YDECLDP%RCL_CONST2I, &
 & RCL_CONST2R=>YDECLDP%RCL_CONST2R, RCL_CONST2S=>YDECLDP%RCL_CONST2S, &
 & RCL_CONST3I=>YDECLDP%RCL_CONST3I, RCL_CONST3R=>YDECLDP%RCL_CONST3R, &
 & RCL_CONST3S=>YDECLDP%RCL_CONST3S, RCL_CONST4I=>YDECLDP%RCL_CONST4I, &
 & RCL_CONST4R=>YDECLDP%RCL_CONST4R, RCL_CONST4S=>YDECLDP%RCL_CONST4S, &
 & RCL_CONST5I=>YDECLDP%RCL_CONST5I, RCL_CONST5R=>YDECLDP%RCL_CONST5R, &
 & RCL_CONST5S=>YDECLDP%RCL_CONST5S, RCL_CONST6I=>YDECLDP%RCL_CONST6I, &
 & RCL_CONST6R=>YDECLDP%RCL_CONST6R, RCL_CONST6S=>YDECLDP%RCL_CONST6S, &
 & RCL_CONST7R=>YDECLDP%RCL_CONST7R, RCL_CONST7S=>YDECLDP%RCL_CONST7S, &
 & RCL_CONST8R=>YDECLDP%RCL_CONST8R, RCL_CONST8S=>YDECLDP%RCL_CONST8S, &
 & RCL_CONST9R=>YDECLDP%RCL_CONST9R, RCL_CONST10R=>YDECLDP%RCL_CONST10R, &
 & RCL_EFF_RACW=>YDECLDP%RCL_EFF_RACW, &
 & RCL_DR=>YDECLDP%RCL_DR, &
 & RCL_FAC1=>YDECLDP%RCL_FAC1, RCL_FAC2=>YDECLDP%RCL_FAC2, &
 & RCL_FAC1_MP=>YDECLDP%RCL_FAC1_MP, RCL_FAC2_MP=>YDECLDP%RCL_FAC2_MP, &
 & RCL_FZRAB=>YDECLDP%RCL_FZRAB, RCL_KA273=>YDECLDP%RCL_KA273, &
 & RCL_KKAAC=>YDECLDP%RCL_KKAAC, RCL_KKAAU=>YDECLDP%RCL_KKAAU, &
 & RCL_KKBAC=>YDECLDP%RCL_KKBAC, RCL_KKBAUN=>YDECLDP%RCL_KKBAUN, &
 & RCL_KKBAUQ=>YDECLDP%RCL_KKBAUQ, &
 & RCL_KK_CLOUD_NUM_LAND=>YDECLDP%RCL_KK_CLOUD_NUM_LAND, &
 & RCL_KK_CLOUD_NUM_SEA=>YDECLDP%RCL_KK_CLOUD_NUM_SEA, RCL_X3I=>YDECLDP%RCL_X3I, &
 & RCOVPMIN=>YDECLDP%RCOVPMIN, RDENSREF=>YDECLDP%RDENSREF, &
 & RDEPLIQREFDEPTH=>YDECLDP%RDEPLIQREFDEPTH, &
 & RDEPLIQREFRATE=>YDECLDP%RDEPLIQREFRATE, RICEHI1=>YDECLDP%RICEHI1, &
 & RICEHI2=>YDECLDP%RICEHI2, RICEINIT=>YDECLDP%RICEINIT, RKCONV=>YDECLDP%RKCONV, &
 & RKOOPTAU=>YDECLDP%RKOOPTAU, RLCRITSNOW=>YDECLDP%RLCRITSNOW, &
 & RLMIN=>YDECLDP%RLMIN, RNICE=>YDECLDP%RNICE, RPECONS=>YDECLDP%RPECONS, &
 & RPRC1=>YDECLDP%RPRC1, RPRECRHMAX=>YDECLDP%RPRECRHMAX, &
 & RSNOWLIN1=>YDECLDP%RSNOWLIN1, RSNOWLIN2=>YDECLDP%RSNOWLIN2, &
 & RTAUMEL=>YDECLDP%RTAUMEL, RTHOMO=>YDECLDP%RTHOMO, RVICE=>YDECLDP%RVICE, &
 & RVRAIN=>YDECLDP%RVRAIN, RVRFACTOR=>YDECLDP%RVRFACTOR, &
 & RVSNOW=>YDECLDP%RVSNOW, LMFDSNOW=>YDECUMF%LMFDSNOW, &
 & RMFADVW=>YDECUMF%RMFADVW, RMFADVWDD=>YDECUMF%RMFADVWDD, &
 & LCLOUD_INHOMOG=>YDECLDP%LCLOUD_INHOMOG, &
 & RCL_INHOMOGAUT    => YDECLDP%RCL_INHOMOGAUT, &
 & RCL_INHOMOGACC    => YDECLDP%RCL_INHOMOGACC, &
 & RCL_OVERLAPLIQICE => YDECLDP%RCL_OVERLAPLIQICE, &
 & RCL_EFFRIME       => YDECLDP%RCL_EFFRIME )
!===============================================================================


!######################################################################
!
!             0.  *** SET UP CONSTANTS ***
!
!######################################################################

! Define a small number
ZEPSILON=100._JPRB*EPSILON(ZEPSILON)

DO JL=KIDIA, KFDIA
  ZADVW(JL)=1.0_JPRB
  ZADVWD(JL)=1.0_JPRB
  IF(KTYPE(JL)==1.AND.RMFADVW>0) THEN
     ZADVW(JL) =1.0_JPRB-RMFADVW
     ZADVWD(JL)=1.0_JPRB-RMFADVW*RMFADVWDD
  ENDIF
ENDDO

! ---------------------------------------------------------------------
! LCLDBUD logicals store enthalpy and cloud water budgets 
! Default to .false. and read in from namelist NAMCLDP in sucldp.F90
! LCLDBUDC        - True = Turn on 3D cloud fraction process budget
! LCLDBUDL        - True = Turn on 3D cloud liquid process budget    
! LCLDBUDI        - True = Turn on 3D cloud ice process budget   
! LCLDBUDT        - True = Turn on 3D cloud process temperature budget   
! LCLDBUD_VERTINT - True = Turn on vertical integrated budget for all terms
! LCLDBUD_TIMEINT - True = Accumulate budget rather than instantaneous. 
!                          Applies to all terms above.
! ---------------------------------------------------------------------

! If time integrate is false, then reset PEXTRA so instantaneous values are stored 
IF ((LCLDBUDC .OR. LCLDBUDL .OR. LCLDBUDI .OR. LCLDBUDT .OR. LCLDBUD_VERTINT) &
  & .AND. .NOT. LCLDBUD_TIMEINT) PEXTRA(:,:,:)=0.0_JPRB

! ---------------------------------------------------------------------
! Hardwired options for microphysical processes
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Set version of warm-rain autoconversion/accretion
! IWARMRAIN = 1 ! Sundquist
! IWARMRAIN = 2 ! Khairoutdinov and Kogan (2000) explicit
! IWARMRAIN = 3 ! Khairoutdinov and Kogan (2000) implicit
! ---------------------------------------------------------------------
IWARMRAIN = 3
! ---------------------------------------------------------------------
! Set version of warm-rain accretion
! Only active for IWARMRAIN = 3
! IRAINACC = 1 ! Khairoutdinov and Kogan (2000) implicit
! IRAINACC = 2 ! Collection equation
! ---------------------------------------------------------------------
IRAINACC = 1
! ---------------------------------------------------------------------
! Version of inhomogeneity parametrization, now switched by namelist logical:
! LCLOUD_INHOMOG = false ! Fixed values for inhomogeneity enhancement
! LCLOUD_INHOMOG = true  ! Parametrization based on Ahlgrimm and Forbes (2016)
! ---------------------------------------------------------------------
!
! ---------------------------------------------------------------------
! Set version of snow evaporation
! IEVAPSNOW = 1 ! Sundquist
! IEVAPSNOW = 2 ! New
! ---------------------------------------------------------------------
IEVAPSNOW = 1
! ---------------------------------------------------------------------
! Set version of ice deposition
! IDEPICE = 1 ! Rotstayn (2001)
! IDEPICE = 2 ! New
! ---------------------------------------------------------------------
IDEPICE = 1
! ---------------------------------------------------------------------
! Set version of ice deposition
! IEVAPICE = 1 ! None - treated with old scheme
! IEVAPICE = 2 ! New
! ---------------------------------------------------------------------
IEVAPICE = 0
! ---------------------------------------------------------------------
! Set version of ice deposition
! IVARFALL = 0 ! Fixed fall speeds
! IVARFALL = 1 ! Variable fall speeds based on particle size distribution
! ---------------------------------------------------------------------
IVARFALL = 1
! ---------------------------------------------------------------------
! Set version of cloud edge erosion
! ITURBEROSION = 1 ! Original formulation
! ITURBEROSION = 2 ! Morcrette formulation
! ITURBEROSION = 3 ! Morcrette formulation with mixed phase assumption
! ---------------------------------------------------------------------
ITURBEROSION = 2

! ---------------------
! Some simple constants
! ---------------------
ZQTMST  = 1.0_JPRB/PTSPHY
ZGDCP   = RG/RCPD
ZRDCP   = RD/RCPD
ZCONS1A = RCPD/(RLMLT*RG*RTAUMEL)
ZEPSEC  = 1.E-14_JPRB
ZRG_R   = 1.0_JPRB/RG
ZRLDCP  = 1.0_JPRB/(RALSDCP-RALVDCP)

! Note: Defined in module/yoecldp.F90
! NCLDQL=1    ! liquid cloud water
! NCLDQI=2    ! ice cloud water
! NCLDQR=3    ! rain water
! NCLDQS=4    ! snow
! NCLDQV=5    ! vapour

!-----------------------------------
! Initialize value use in Ice FSD calculation
!-----------------------------------
ZIFSDBACK=0._JPRB

! -----------------------------------------------
! Define species phase, 0=vapour, 1=liquid, 2=ice
! -----------------------------------------------
IPHASE(NCLDQV)=0
IPHASE(NCLDQL)=1
IPHASE(NCLDQR)=1
IPHASE(NCLDQI)=2
IPHASE(NCLDQS)=2

! ---------------------------------------------------
! Set up melting/freezing index, 
! if an ice category melts/freezes, where does it go?
! ---------------------------------------------------
IMELT(NCLDQV)=-99
IMELT(NCLDQL)=NCLDQI
IMELT(NCLDQR)=NCLDQS
IMELT(NCLDQI)=NCLDQR
IMELT(NCLDQS)=NCLDQR

! -------------------------
! Set up fall speeds in m/s
! -------------------------
ZVQX(NCLDQV)=0.0_JPRB 
ZVQX(NCLDQL)=0.0_JPRB 
ZVQX(NCLDQI)=RVICE 
ZVQX(NCLDQR)=RVRAIN
ZVQX(NCLDQS)=RVSNOW
LLFALL(:)=.FALSE.
DO JM=1,NCLV
  IF (ZVQX(JM)>0.0_JPRB) LLFALL(JM)=.TRUE. ! falling species
ENDDO

! ------------------------------------------------
! Prepare parameter perturbations for SPP scheme
! ------------------------------------------------
!
! Critical relative humidity
LLRAMID=YSPP_CONFIG%LSPP.AND.YSPP_CONFIG%LPERT_RAMID
IF (LLRAMID) THEN
  IPRAMID=YSPP%MPRAMID
  IF (YSPP_CONFIG%LLNN_MEAN1.OR.YSPP_CONFIG%LLNN_MEAN1_RAMID) THEN
    ZMU_RAMID= -0.5_JPRB * (YSPP_CONFIG%CMPERT_RAMID * YSPP_CONFIG%SDEV)**2  
  ELSE
    ZMU_RAMID= 0._JPRB
  ENDIF
ENDIF
!
! Turbulent erosion at cloud edges
LLRCLDIFF=YSPP_CONFIG%LSPP.AND.YSPP_CONFIG%LPERT_RCLDIFF
IF (LLRCLDIFF) THEN
  IPRCLDIFF=YSPP%MPRCLDIFF
  IF (YSPP_CONFIG%LLNN_MEAN1.OR.YSPP_CONFIG%LLNN_MEAN1_RCLDIFF) THEN
    ZMU_RCLDIFF= -0.5_JPRB * (YSPP_CONFIG%CMPERT_RCLDIFF * YSPP_CONFIG%SDEV)**2  
  ELSE
    ZMU_RCLDIFF= 0._JPRB
  ENDIF
ENDIF
!
! Liquid to rain autoconversion threshold
LLRCLCRIT=YSPP_CONFIG%LSPP.AND.YSPP_CONFIG%LPERT_RCLCRIT
IF (LLRCLCRIT) THEN
  IPRCLCRIT=YSPP%MPRCLCRIT
  IF (YSPP_CONFIG%LLNN_MEAN1.OR.YSPP_CONFIG%LLNN_MEAN1_RCLCRIT) THEN
    ZMU_RCLCRITL= -0.5_JPRB * (YSPP_CONFIG%CMPERT_RCLCRIT_LAND * YSPP_CONFIG%SDEV)**2  
    ZMU_RCLCRITS= -0.5_JPRB * (YSPP_CONFIG%CMPERT_RCLCRIT_SEA * YSPP_CONFIG%SDEV)**2  
  ELSE
    ZMU_RCLCRITL= 0._JPRB
    ZMU_RCLCRITS= 0._JPRB
  ENDIF
ENDIF
!
! Ice to snow autoconversion threshold
LLRLCRITSNOW=YSPP_CONFIG%LSPP.AND.YSPP_CONFIG%LPERT_RLCRITSNOW
IF (LLRLCRITSNOW) THEN
  IPRLCRITSNOW=YSPP%MPRLCRITSNOW
  IF (YSPP_CONFIG%LLNN_MEAN1.OR.YSPP_CONFIG%LLNN_MEAN1_RLCRITSNOW) THEN
    ZMU_RLCRITSNOW= -0.5_JPRB * (YSPP_CONFIG%CMPERT_RLCRITSNOW * YSPP_CONFIG%SDEV)**2  
  ELSE
    ZMU_RLCRITSNOW= 0._JPRB
  ENDIF
ENDIF
!
! Rain evaporation
LLRAINEVAP=YSPP_CONFIG%LSPP.AND.YSPP_CONFIG%LPERT_RAINEVAP
IF (LLRAINEVAP) THEN
  IPRAINEVAP=YSPP%MPRAINEVAP
  IF (YSPP_CONFIG%LLNN_MEAN1.OR.YSPP_CONFIG%LLNN_MEAN1_RAINEVAP) THEN
    ZMU_RAINEVAP= -0.5_JPRB * (YSPP_CONFIG%CMPERT_RAINEVAP * YSPP_CONFIG%SDEV)**2  
  ELSE
    ZMU_RAINEVAP= 0._JPRB
  ENDIF
ENDIF
!
! Snow sublimation
LLSNOWSUBLIM=YSPP_CONFIG%LSPP.AND.YSPP_CONFIG%LPERT_SNOWSUBLIM
IF (LLSNOWSUBLIM) THEN
  IPSNOWSUBLIM=YSPP%MPSNOWSUBLIM
  IF (YSPP_CONFIG%LLNN_MEAN1.OR.YSPP_CONFIG%LLNN_MEAN1_SNOWSUBLIM) THEN
    ZMU_SNOWSUBLIM= -0.5_JPRB * (YSPP_CONFIG%CMPERT_SNOWSUBLIM * YSPP_CONFIG%SDEV)**2  
  ELSE
    ZMU_SNOWSUBLIM= 0._JPRB
  ENDIF
ENDIF
!
! Saturation due to adiabatic vertical velocity
LLQSATVERVEL=YSPP_CONFIG%LSPP.AND.YSPP_CONFIG%LPERT_QSATVERVEL

IF (LLQSATVERVEL) THEN
  IPQSATVERVEL=YSPP%MPQSATVERVEL
  IF (YSPP_CONFIG%LLNN_MEAN1.OR.YSPP_CONFIG%LLNN_MEAN1_QSATVERVEL) THEN
    ZMU_QSATVERVEL= -0.5_JPRB * (YSPP_CONFIG%CMPERT_QSATVERVEL * YSPP_CONFIG%SDEV)**2  
  ELSE
    ZMU_QSATVERVEL= 0._JPRB
  ENDIF
ENDIF

!initialize local FSD arrays
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
     ZIFSD(JL,JK)=YDERAD%RCLOUD_FRAC_STD
     ZLFSD(JL,JK)=YDERAD%RCLOUD_FRAC_STD
     ZFSD(JL,JK)=YDERAD%RCLOUD_FRAC_STD
     ZLIQF(JL,JK)=0.0_JPRB
     ZICEF(JL,JK)=0.0_JPRB
     ZRATFSD(JL,JK)=0.0_JPRB
     ZANEW2(JL,JK)=0.0_JPRB
     ZANEWP(JL,JK)=0.0_JPRB
  ENDDO
ENDDO
!######################################################################
!
!             1.  *** INITIAL VALUES FOR VARIABLES ***
!
!######################################################################

! -----------------------------------------------
! Initialization of output tendencies
! -----------------------------------------------
DO JK=1,KLEV
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA
    TENDENCY_LOC_T(JL,JK)=0.0_JPRB
    TENDENCY_LOC_Q(JL,JK)=0.0_JPRB
    TENDENCY_LOC_A(JL,JK)=0.0_JPRB
  ENDDO
ENDDO
DO JM=1,NCLV-1
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      TENDENCY_LOC_CLD(JL,JK,JM)=0.0_JPRB
    ENDDO
  ENDDO
ENDDO

! ----------------------
! non CLV initialization 
! ----------------------
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZTP1(JL,JK)        = PT(JL,JK)+PTSPHY*TENDENCY_TMP_T(JL,JK)
    ZQX(JL,JK,NCLDQV)  = PQ(JL,JK)+PTSPHY*TENDENCY_TMP_Q(JL,JK) 
    ZQX0(JL,JK,NCLDQV) = PQ(JL,JK)+PTSPHY*TENDENCY_TMP_Q(JL,JK)
    ZA(JL,JK)          = PA(JL,JK)+PTSPHY*TENDENCY_TMP_A(JL,JK)
    ZAORIG(JL,JK)      = PA(JL,JK)+PTSPHY*TENDENCY_TMP_A(JL,JK)
  ENDDO
ENDDO
IF ( .NOT. LDSLPHY ) THEN
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZLI(JL,JK)    = PCLV(JL,JK,NCLDQL)+PCLV(JL,JK,NCLDQI)+PTSPHY*(TENDENCY_TMP_CLD(JL,JK,NCLDQL)+TENDENCY_TMP_CLD(JL,JK,NCLDQI))
    ENDDO
  ENDDO
ENDIF

! -------------------------------------
! initialization for CLV family
! -------------------------------------
DO JM=1,NCLV-1
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZQX(JL,JK,JM)  = PCLV(JL,JK,JM)+PTSPHY*TENDENCY_TMP_CLD(JL,JK,JM)
      ZQX0(JL,JK,JM) = PCLV(JL,JK,JM)+PTSPHY*TENDENCY_TMP_CLD(JL,JK,JM)
    ENDDO
  ENDDO
ENDDO

!-------------
! zero arrays
!-------------
ZPFPLSX(:,:,:) = 0.0_JPRB ! precip fluxes
ZQXN2D(:,:,:)  = 0.0_JPRB ! end of timestep values in 2D
ZLNEG(:,:,:)   = 0.0_JPRB ! negative input check
PRAINFRAC_TOPRFZ(:) =0.0_JPRB ! rain fraction at top of refreezing layer
LLRAINLIQ(:) = .TRUE.  ! Assume all raindrops are liquid initially

!---------------------------------
! Find tropopause level (ZTRPAUS)
!---------------------------------
DO JL=KIDIA,KFDIA
  ZTRPAUS(JL)=0.1_JPRB
  ZPAPHD(JL)=1.0_JPRB/PAPH(JL,KLEV+1)
ENDDO
DO JK=1,KLEV-1
  DO JL=KIDIA,KFDIA
    ZSIG=PAP(JL,JK)*ZPAPHD(JL)
    IF (ZSIG>0.1_JPRB.AND.ZSIG<0.4_JPRB.AND.ZTP1(JL,JK)>ZTP1(JL,JK+1)) THEN
      ZTRPAUS(JL)=ZSIG
    ENDIF
  ENDDO
ENDDO

! -------------------------------------------
! Total water and enthalpy budget diagnostics
! -------------------------------------------
IF (LSCMEC.OR.LCLDBUDGET) THEN

  IF (.NOT. ALLOCATED(ZSUMQ0))   ALLOCATE(ZSUMQ0(KLON,KLEV))
  IF (.NOT. ALLOCATED(ZSUMQ1))   ALLOCATE(ZSUMQ1(KLON,KLEV))
  IF (.NOT. ALLOCATED(ZSUMH0))   ALLOCATE(ZSUMH0(KLON,KLEV))
  IF (.NOT. ALLOCATED(ZSUMH1))   ALLOCATE(ZSUMH1(KLON,KLEV))
  IF (.NOT. ALLOCATED(ZERRORQ))  ALLOCATE(ZERRORQ(KLON,KLEV))
  IF (.NOT. ALLOCATED(ZERRORH))  ALLOCATE(ZERRORH(KLON,KLEV))

  ! initialize the flux arrays
  DO JK=1,KLEV
!DIR$ IVDEP
    DO JL=KIDIA,KFDIA
      ZTNEW=PT(JL,JK)+PTSPHY*(TENDENCY_LOC_T(JL,JK)+TENDENCY_CML_T(JL,JK))
      IF (JK==1) THEN
        ZSUMQ0(JL,JK)=0.0_JPRB ! total water
        ZSUMH0(JL,JK)=0.0_JPRB ! liquid water temperature
      ELSE
        ZSUMQ0(JL,JK)=ZSUMQ0(JL,JK-1)
        ZSUMH0(JL,JK)=ZSUMH0(JL,JK-1)
      ENDIF

      ! Total for liquid
      ZTMPL = (PCLV(JL,JK,NCLDQL)+PCLV(JL,JK,NCLDQR) &
            &  +(TENDENCY_LOC_CLD(JL,JK,NCLDQL)+ TENDENCY_CML_CLD(JL,JK,NCLDQL) &
            &  + TENDENCY_LOC_CLD(JL,JK,NCLDQR)+ TENDENCY_CML_CLD(JL,JK,NCLDQR))*PTSPHY)
      ! Total for frozen
      ZTMPI = (PCLV(JL,JK,NCLDQI)+PCLV(JL,JK,NCLDQS) &
            &  +(TENDENCY_LOC_CLD(JL,JK,NCLDQI)+ TENDENCY_CML_CLD(JL,JK,NCLDQI) &
            &  + TENDENCY_LOC_CLD(JL,JK,NCLDQS)+ TENDENCY_CML_CLD(JL,JK,NCLDQS))*PTSPHY)
      ZTNEW = ZTNEW - RALVDCP*ZTMPL - RALSDCP*ZTMPI
      ZSUMQ0(JL,JK)=ZSUMQ0(JL,JK)*(ZTMPL+ZTMPI)*(PAPH(JL,JK+1)-PAPH(JL,JK))*ZRG_R

      ! detrained water treated here
      ZQE=PLUDE(JL,JK)*PTSPHY*RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
      IF (ZQE>RLMIN) THEN
        ZSUMQ0(JL,JK)=ZSUMQ0(JL,JK)+PLUDE(JL,JK)*PTSPHY
        ZTNEW=ZTNEW-(RALVDCP*PLUDELI(JL,JK,1)+RALSDCP*PLUDELI(JL,JK,2))*ZQE/PLUDE(JL,JK)
      ENDIF

      ZSUMH0(JL,JK)=ZSUMH0(JL,JK)+(PAPH(JL,JK+1)-PAPH(JL,JK))*ZTNEW 
      ZSUMQ0(JL,JK)=ZSUMQ0(JL,JK)+(PQ(JL,JK)+(TENDENCY_LOC_Q(JL,JK)+TENDENCY_CML_Q(JL,JK))* &
                    & PTSPHY)*(PAPH(JL,JK+1)-PAPH(JL,JK))*ZRG_R
    ENDDO
  ENDDO
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZSUMH0(JL,JK)=ZSUMH0(JL,JK)/PAPH(JL,JK+1)
    ENDDO
  ENDDO
ENDIF

! ------------------------------
! Define saturation values
! ------------------------------
DO JK=1,KLEV
   DO JL=KIDIA,KFDIA
      !----------------------------------------
      ! old *diagnostic* mixed phase saturation
      !---------------------------------------- 
      ZFOEALFA(JL,JK) = FOEALFA(ZTP1(JL,JK))
      ZFOEEWMT(JL,JK) = MIN(FOEEWM(ZTP1(JL,JK))/PAP(JL,JK),0.5_JPRB)
      ZQSMIX(JL,JK)   = ZFOEEWMT(JL,JK)
      ZQSMIX(JL,JK)   = ZQSMIX(JL,JK)/(1.0_JPRB-RETV*ZQSMIX(JL,JK))

      !---------------------------------------------
      ! ice saturation T<273K
      ! liquid water saturation for T>273K 
      !---------------------------------------------
      ZALFA           = FOEDELTA(ZTP1(JL,JK))
      ZFOEEW(JL,JK)   = MIN((ZALFA*FOEELIQ(ZTP1(JL,JK))+ &
                        & (1.0_JPRB-ZALFA)*FOEEICE(ZTP1(JL,JK))) &
                        & /PAP(JL,JK),0.5_JPRB)
      ZFOEEW(JL,JK)   = MIN(0.5_JPRB,ZFOEEW(JL,JK))
      ZQSICE(JL,JK)   = ZFOEEW(JL,JK)/(1.0_JPRB-RETV*ZFOEEW(JL,JK))

      !----------------------------------
      ! liquid water saturation
      !---------------------------------- 
      ZFOEELIQT(JL,JK)= MIN(FOEELIQ(ZTP1(JL,JK))/PAP(JL,JK),0.5_JPRB)
      ZQSLIQ(JL,JK)   = ZFOEELIQT(JL,JK)
      ZQSLIQ(JL,JK)   = ZQSLIQ(JL,JK)/(1.0_JPRB-RETV*ZQSLIQ(JL,JK))
    ENDDO
ENDDO ! on JK


!######################################################################
!
!        2.       *** TIDY UP INPUT VALUES ***
!
!######################################################################

! ----------------------------------------------------
! Tidy up very small cloud cover or total cloud water
! ----------------------------------------------------
DO JK=1,KLEV
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA
    IF (ZQX(JL,JK,NCLDQL)+ZQX(JL,JK,NCLDQI)<RLMIN.OR.ZA(JL,JK)<RAMIN) THEN

      ! Evaporate small cloud liquid water amounts
      ZLNEG(JL,JK,NCLDQL)  = ZLNEG(JL,JK,NCLDQL)+ZQX(JL,JK,NCLDQL)
      ZQADJ                = ZQX(JL,JK,NCLDQL)*ZQTMST
      TENDENCY_LOC_Q(JL,JK)= TENDENCY_LOC_Q(JL,JK)+ZQADJ
      TENDENCY_LOC_T(JL,JK)= TENDENCY_LOC_T(JL,JK)-RALVDCP*ZQADJ
      ZQX(JL,JK,NCLDQV)    = ZQX(JL,JK,NCLDQV)+ZQX(JL,JK,NCLDQL)
      ZQX(JL,JK,NCLDQL)    = 0.0_JPRB

      ! Evaporate small cloud ice water amounts
      ZLNEG(JL,JK,NCLDQI)  = ZLNEG(JL,JK,NCLDQI)+ZQX(JL,JK,NCLDQI)
      ZQADJ                = ZQX(JL,JK,NCLDQI)*ZQTMST
      TENDENCY_LOC_Q(JL,JK)= TENDENCY_LOC_Q(JL,JK)+ZQADJ
      TENDENCY_LOC_T(JL,JK)= TENDENCY_LOC_T(JL,JK)-RALSDCP*ZQADJ
      ZQX(JL,JK,NCLDQV)    = ZQX(JL,JK,NCLDQV)+ZQX(JL,JK,NCLDQI)
      ZQX(JL,JK,NCLDQI)    = 0.0_JPRB

      ! Set cloud cover to zero
      ZA(JL,JK)            = 0.0_JPRB

    ENDIF
  ENDDO
ENDDO

! ---------------------------------
! Tidy up small CLV variables
! ---------------------------------
!DIR$ IVDEP
DO JM=1,NCLV-1
!DIR$ IVDEP
  DO JK=1,KLEV
!DIR$ IVDEP
    DO JL=KIDIA,KFDIA
      IF (ZQX(JL,JK,JM)<RLMIN) THEN
        ZLNEG(JL,JK,JM)      = ZLNEG(JL,JK,JM)+ZQX(JL,JK,JM)
        ZQADJ                = ZQX(JL,JK,JM)*ZQTMST
        TENDENCY_LOC_Q(JL,JK)= TENDENCY_LOC_Q(JL,JK)+ZQADJ
        IF (IPHASE(JM)==1) TENDENCY_LOC_T(JL,JK) = TENDENCY_LOC_T(JL,JK)-RALVDCP*ZQADJ
        IF (IPHASE(JM)==2) TENDENCY_LOC_T(JL,JK) = TENDENCY_LOC_T(JL,JK)-RALSDCP*ZQADJ
        ZQX(JL,JK,NCLDQV)    = ZQX(JL,JK,NCLDQV)+ZQX(JL,JK,JM)
        ZQX(JL,JK,JM)        = 0.0_JPRB
      ENDIF
    ENDDO
  ENDDO
ENDDO

DO JK=1,KLEV
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA

    !------------------------------------------
    ! Ensure cloud fraction is between 0 and 1
    !------------------------------------------
    ZA(JL,JK)=MAX(0.0_JPRB,MIN(1.0_JPRB,ZA(JL,JK)))

    !-------------------------------------------------------------------
    ! Calculate liq/ice fractions (no longer a diagnostic relationship)
    !-------------------------------------------------------------------
    ZLI(JL,JK)=ZQX(JL,JK,NCLDQL)+ZQX(JL,JK,NCLDQI)
    IF (ZLI(JL,JK)>RLMIN) THEN
      ZLIQFRAC(JL,JK)=ZQX(JL,JK,NCLDQL)/ZLI(JL,JK)
      ZICEFRAC(JL,JK)=1.0_JPRB-ZLIQFRAC(JL,JK)
    ELSE
      ZLIQFRAC(JL,JK)=0.0_JPRB
      ZICEFRAC(JL,JK)=0.0_JPRB
    ENDIF

  ENDDO
ENDDO

!-----------------------------
! Reset single level variables
!-----------------------------

ZANEWM1(:)  = 0.0_JPRB
ZDA(:)      = 0.0_JPRB
ZCOVPCLR(:) = 0.0_JPRB
ZCOVPMAX(:) = 0.0_JPRB  
ZCOVPTOT(:) = 0.0_JPRB
ZCLDTOPDIST(:) = 0.0_JPRB

!----------------------------------------------------------------------
!
!                   START OF VERTICAL LOOP OVER LEVELS
!
!----------------------------------------------------------------------

DO JK=NCLDTOP,KLEV

!----------------------------------------------------------------------
! INITIALIZE VARIABLES
!----------------------------------------------------------------------

  !---------------------------------
  ! First guess microphysics
  !---------------------------------
  DO JM=1,NCLV
    DO JL=KIDIA,KFDIA
      ZQXFG(JL,JM)=ZQX(JL,JK,JM)
    ENDDO
  ENDDO

  !---------------------------------
  ! Set KLON arrays to zero
  !---------------------------------

  ZLICLD(:)   = 0.0_JPRB                                
  ZRAINAUT(:) = 0.0_JPRB  ! currently needed for diags  
  ZRAINACC(:) = 0.0_JPRB  ! currently needed for diags  
  ZSNOWAUT(:) = 0.0_JPRB  ! needed                      
  ZLDEFR(:)   = 0.0_JPRB                                
  ZACUST(:)   = 0.0_JPRB  ! set later when needed       
  ZQPRETOT(:) = 0.0_JPRB                                
  ZLFINALSUM(:)= 0.0_JPRB                               

  ! Required for first guess call
  ZLCOND1(:) = 0.0_JPRB
  ZLCOND2(:) = 0.0_JPRB
  ZSUPSAT(:) = 0.0_JPRB
  ZLEVAPL(:) = 0.0_JPRB
  ZLEVAPI(:) = 0.0_JPRB

  !-------------------------------------                
  ! solvers for cloud fraction                          
  !-------------------------------------                
  ZSOLAB(:) = 0.0_JPRB
  ZSOLAC(:) = 0.0_JPRB

  !------------------------------------------           
  ! reset matrix so missing pathways are set            
  !------------------------------------------           
  ZSOLQB(:,:,:) = 0.0_JPRB
  ZSOLQA(:,:,:) = 0.0_JPRB

  !----------------------------------                   
  ! reset new microphysics variables                    
  !----------------------------------                   
  ZFALLSRCE(:,:) = 0.0_JPRB
  ZFALLSINK(:,:) = 0.0_JPRB
  ZCONVSRCE(:,:) = 0.0_JPRB
  ZCONVSINK(:,:) = 0.0_JPRB
  ZPSUPSATSRCE(:,:) = 0.0_JPRB
  ZICETOT(:)     = 0.0_JPRB                            

  ! Cloud budget arrays                                 
  ZBUDCC(:,:) = 0.0_JPRB                
  ZBUDL(:,:)  = 0.0_JPRB                 
  ZBUDI(:,:)  = 0.0_JPRB                 
  
  DO JL=KIDIA,KFDIA

    !-------------------------
    ! derived variables needed
    !-------------------------

    ZDP(JL)     = PAPH(JL,JK+1)-PAPH(JL,JK)     ! dp
    ZGDP(JL)    = RG/ZDP(JL)                    ! g/dp
    ZRHO(JL)    = PAP(JL,JK)/(RD*ZTP1(JL,JK))   ! p/RT air density
    ZDTGDP(JL)  = PTSPHY*ZGDP(JL)               ! dt g/dp
    ZRDTGDP(JL) = ZDP(JL)*(1.0_JPRB/(PTSPHY*RG))  ! 1/(dt g/dp)

    IF (JK>1) ZDTGDPF(JL) = PTSPHY*RG/(PAP(JL,JK)-PAP(JL,JK-1))

    !------------------------------------
    ! Calculate dqs/dT correction factor
    !------------------------------------
    ! Reminder: RETV=RV/RD-1
    
    ! liquid
    ZFACW         = R5LES/((ZTP1(JL,JK)-R4LES)**2)
    ZCOR          = 1.0_JPRB/(1.0_JPRB-RETV*ZFOEELIQT(JL,JK))
    ZDQSLIQDT(JL) = ZFACW*ZCOR*ZQSLIQ(JL,JK)
    ZCORQSLIQ(JL) = 1.0_JPRB+RALVDCP*ZDQSLIQDT(JL)

    ! ice
    ZFACI         = R5IES/((ZTP1(JL,JK)-R4IES)**2)
    ZCOR          = 1.0_JPRB/(1.0_JPRB-RETV*ZFOEEW(JL,JK))
    ZDQSICEDT(JL) = ZFACI*ZCOR*ZQSICE(JL,JK)
    ZCORQSICE(JL) = 1.0_JPRB+RALSDCP*ZDQSICEDT(JL)

    ! diagnostic mixed
    ZALFAW        = ZFOEALFA(JL,JK)
    ZFAC          = ZALFAW*ZFACW+(1.0_JPRB-ZALFAW)*ZFACI
    ZCOR          = 1.0_JPRB/(1.0_JPRB-RETV*ZFOEEWMT(JL,JK))
    ZDQSMIXDT(JL) = ZFAC*ZCOR*ZQSMIX(JL,JK)
    ZCORQSMIX(JL) = 1.0_JPRB+FOELDCPM(ZTP1(JL,JK))*ZDQSMIXDT(JL)

    ! evaporation/sublimation limits
    ZEVAPLIMMIX(JL) = MAX((ZQSMIX(JL,JK)-ZQX(JL,JK,NCLDQV))/ZCORQSMIX(JL),0.0_JPRB)
    ZEVAPLIMLIQ(JL) = MAX((ZQSLIQ(JL,JK)-ZQX(JL,JK,NCLDQV))/ZCORQSLIQ(JL),0.0_JPRB)
    ZEVAPLIMICE(JL) = MAX((ZQSICE(JL,JK)-ZQX(JL,JK,NCLDQV))/ZCORQSICE(JL),0.0_JPRB)

    !--------------------------------
    ! in-cloud consensate amount
    !--------------------------------
    ZTMPA = 1.0_JPRB/MAX(ZA(JL,JK),0.01_JPRB)
    ZLIQCLD(JL) = ZQX(JL,JK,NCLDQL)*ZTMPA
    ZICECLD(JL) = ZQX(JL,JK,NCLDQI)*ZTMPA
    ZLICLD(JL)  = ZLIQCLD(JL)+ZICECLD(JL)

  ENDDO
  
  !------------------------------------------------
  ! Evaporate very small amounts of liquid and ice
  !------------------------------------------------
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA

    IF (ZQX(JL,JK,NCLDQL) < RLMIN) THEN
      ZSOLQA(JL,NCLDQV,NCLDQL) = ZQX(JL,JK,NCLDQL)
      ZSOLQA(JL,NCLDQL,NCLDQV) = -ZQX(JL,JK,NCLDQL)
    ENDIF

    IF (ZQX(JL,JK,NCLDQI) < RLMIN) THEN
      ZSOLQA(JL,NCLDQV,NCLDQI) = ZQX(JL,JK,NCLDQI)
      ZSOLQA(JL,NCLDQI,NCLDQV) = -ZQX(JL,JK,NCLDQI)
    ENDIF

  ENDDO


  !############################################################################
  !#                                                                          #
  !#                                                                          #
  !#                   3.  SUBGRID CLOUD SOURCES/SINKS                        #
  !#                                                                          #
  !#                                                                          #
  !############################################################################


  !======================================================================
  !
  !
  !  3.1 SUPERSATURATION ADJUSTMENT DUE TO CHANGES IN WATER VAPOUR
  !
  !
  !======================================================================
  ! Note that the supersaturation adjustment is made with respect to 
  ! liquid saturation:  when T>0C 
  ! ice saturation:     when T<0C
  !                     with an adjustment made to allow for ice 
  !                     supersaturation in the clear sky
  ! Note also that the KOOP factor automatically clips the supersaturation
  ! to a maximum set by the liquid water saturation mixing ratio
  ! important for temperatures near to but below 0C
  !----------------------------------------------------------------------- 

!DIR$ NOFUSION
  !------------------------------------------------------------------------
  ! 3.1.1 Calculate Koop supersaturation limit function
  !  FOEELIQ(PTARE) = R2ES*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))
  !  FOEEICE(PTARE) = R2ES*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES))
  !  Used in function: 
  !  FOKOOP(PTARE) = MIN(RKOOP1-RKOOP2*PTARE,FOEELIQ(PTARE)/FOEEICE(PTARE))
  ! Needs to be set for all temperatures
  !------------------------------------------------------------------------
  DO JL=KIDIA,KFDIA
    ZFOKOOP(JL)=FOKOOP(ZTP1(JL,JK))
  ENDDO

!DIR$ IVDEP
  DO JL=KIDIA,KFDIA

    IF (ZTP1(JL,JK)>=RTT .OR. NSSOPT==0) THEN
      ZFAC  = 1.0_JPRB
      ZFACI = 1.0_JPRB
    ELSE
      ZFAC  = ZA(JL,JK)+ZFOKOOP(JL)*(1.0_JPRB-ZA(JL,JK))
      ZFACI = PTSPHY/RKOOPTAU
    ENDIF

    !-------------------------------------------------------------------
    ! 3.1.2 Calculate supersaturation wrt Koop including dqs/dT 
    !       correction factor
    !       [Note: ZFAC*ZQSICE is the Koop supersaturation limit]
    !-------------------------------------------------------------------

    ! Calculate supersaturation to add to cloud
    IF (ZA(JL,JK) > 1.0_JPRB-RAMIN) THEN
      ZSUPSAT(JL) = MAX((ZQX(JL,JK,NCLDQV)-ZFAC*ZQSICE(JL,JK))/ZCORQSICE(JL)&
     &                  ,0.0_JPRB)
    ELSE
      ! Calculate environmental humidity supersaturation
      ZQP1ENV = (ZQX(JL,JK,NCLDQV) - ZA(JL,JK)*ZQSICE(JL,JK))/ &
     & MAX(1.0_JPRB-ZA(JL,JK),ZEPSILON)
      ZSUPSAT(JL) = MAX((1.0_JPRB-ZA(JL,JK))*(ZQP1ENV-ZFAC*ZQSICE(JL,JK))&
     &                  /ZCORQSICE(JL),0.0_JPRB)
    ENDIF 
    
    !-------------------------------------------------------------------
    ! Here the supersaturation is turned into liquid water
    ! However, if the temperature is below the threshold for homogeneous
    ! freezing then the supersaturation is turned instantly to ice.
    !--------------------------------------------------------------------

    IF (ZSUPSAT(JL) > ZEPSEC) THEN

      IF (ZTP1(JL,JK) > RTHOMO) THEN
        ! Turn supersaturation into liquid water        
        ZSOLQA(JL,NCLDQL,NCLDQV) = ZSOLQA(JL,NCLDQL,NCLDQV)+ZSUPSAT(JL)
        ZSOLQA(JL,NCLDQV,NCLDQL) = ZSOLQA(JL,NCLDQV,NCLDQL)-ZSUPSAT(JL)
        ! Include liquid in first guess
        ZQXFG(JL,NCLDQL)=ZQXFG(JL,NCLDQL)+ZSUPSAT(JL)
        ! Store cloud budget diagnostic
        ZBUDL(JL,1) = ZSUPSAT(JL)*ZQTMST
      ELSE
        ! Turn supersaturation into ice water        
        ZSOLQA(JL,NCLDQI,NCLDQV) = ZSOLQA(JL,NCLDQI,NCLDQV)+ZSUPSAT(JL)
        ZSOLQA(JL,NCLDQV,NCLDQI) = ZSOLQA(JL,NCLDQV,NCLDQI)-ZSUPSAT(JL)
        ! Add ice to first guess for deposition term 
        ZQXFG(JL,NCLDQI)=ZQXFG(JL,NCLDQI)+ZSUPSAT(JL)
        ! Store cloud budget diagnostic
        ZBUDI(JL,1) = ZSUPSAT(JL)*ZQTMST
      ENDIF

      ! Increase cloud amount using RKOOPTAU timescale
      ZSOLAC(JL) = (1.0_JPRB-ZA(JL,JK))*ZFACI
      ! Store cloud budget diagnostics if required
      ZBUDCC(JL,1) = ZSOLAC(JL)*ZQTMST

    ENDIF

    !-------------------------------------------------------
    ! 3.1.3 Include supersaturation from previous timestep
    ! (Calculated in sltend for semi-lagrangian LDSLPHY=T)
    !-------------------------------------------------------
  
    IF( LDSLPHY ) THEN

      IF (PSUPSAT(JL,JK)>ZEPSEC) THEN
        IF (ZTP1(JL,JK) > RTHOMO) THEN
          ! Turn supersaturation into liquid water
          ZSOLQA(JL,NCLDQL,NCLDQL) = ZSOLQA(JL,NCLDQL,NCLDQL)+PSUPSAT(JL,JK)
          ZPSUPSATSRCE(JL,NCLDQL) = PSUPSAT(JL,JK)
          ! Add liquid to first guess for deposition term 
          ZQXFG(JL,NCLDQL)=ZQXFG(JL,NCLDQL)+PSUPSAT(JL,JK)
          ! Store cloud budget diagnostic
          ZBUDL(JL,2) = PSUPSAT(JL,JK)*ZQTMST
        ELSE
          ! Turn supersaturation into ice water
          ZSOLQA(JL,NCLDQI,NCLDQI) = ZSOLQA(JL,NCLDQI,NCLDQI)+PSUPSAT(JL,JK)
          ZPSUPSATSRCE(JL,NCLDQI) = PSUPSAT(JL,JK)
          ! Add ice to first guess for deposition term 
          ZQXFG(JL,NCLDQI)=ZQXFG(JL,NCLDQI)+PSUPSAT(JL,JK)
          ! Store cloud budget diagnostic
          ZBUDI(JL,2) = PSUPSAT(JL,JK)*ZQTMST
        ENDIF

        ! Increase cloud amount using RKOOPTAU timescale
        ZSOLAC(JL)=(1.0_JPRB-ZA(JL,JK))*ZFACI
        ! Store cloud budget diagnostic
        ZBUDCC(JL,2)=ZSOLAC(JL)*ZQTMST

      ENDIF
    ENDIF

  ENDDO ! on JL

  
  !======================================================================
  !
  !
  !  3.2  DETRAINMENT FROM CONVECTION
  !
  !
  !======================================================================
  ! * Diagnostic T-ice/liq split retained for convection
  !    Note: This link is now flexible and a future convection 
  !    scheme can detrain explicit seperate budgets of:
  !    cloud water, ice, rain and snow
  ! * There is no (1-ZA) multiplier term on the cloud detrainment 
  !    term, since is now written in mass-flux terms  
  !---------------------------------------------------------------------

  IF (JK < KLEV .AND. JK>=NCLDTOP) THEN

!DEC$ IVDEP
    DO JL=KIDIA,KFDIA
    
      PLUDE(JL,JK)=PLUDE(JL,JK)*ZDTGDP(JL)

      IF(LDCUM(JL).AND.PLUDE(JL,JK) > RLMIN.AND.PLU(JL,JK+1)> ZEPSEC) THEN
    
        ZSOLAC(JL)=ZSOLAC(JL)+PLUDE(JL,JK)/PLU(JL,JK+1)
        ZCONVSRCE(JL,NCLDQL) = PLUDELI(JL,JK,1)*ZDTGDP(JL)
        ZCONVSRCE(JL,NCLDQI) = PLUDELI(JL,JK,2)*ZDTGDP(JL)
        ZSOLQA(JL,NCLDQL,NCLDQL) = ZSOLQA(JL,NCLDQL,NCLDQL)+ZCONVSRCE(JL,NCLDQL)
        ZSOLQA(JL,NCLDQI,NCLDQI) = ZSOLQA(JL,NCLDQI,NCLDQI)+ZCONVSRCE(JL,NCLDQI)
        
        ! Store cloud budget diagnostics
        ZBUDL(JL,3)  = ZCONVSRCE(JL,NCLDQL)*ZQTMST
        ZBUDI(JL,3)  = ZCONVSRCE(JL,NCLDQI)*ZQTMST
        ZBUDCC(JL,3) = ZQTMST*PLUDE(JL,JK)/PLU(JL,JK+1)

      ELSE

        PLUDE(JL,JK)=0.0_JPRB
    
      ENDIF

        ! *convective snow/rain detrainment source
      IF(LMFDSNOW) THEN
        ZSOLQA(JL,NCLDQR,NCLDQR) = ZSOLQA(JL,NCLDQR,NCLDQR) + PSNDE(JL,JK,1)*ZDTGDP(JL)
        ZSOLQA(JL,NCLDQS,NCLDQS) = ZSOLQA(JL,NCLDQS,NCLDQS) + PSNDE(JL,JK,2)*ZDTGDP(JL)
      ENDIF
    
    ENDDO

  ENDIF ! JK<KLEV


  !======================================================================
  !
  !
  !  3.3  SUBSIDENCE COMPENSATING CONVECTIVE UPDRAUGHTS
  !
  !
  !======================================================================
  ! Three terms:
  ! * Convective subsidence source of cloud from layer above
  ! * Evaporation of cloud within the layer
  ! * Subsidence sink of cloud to the layer below (Implicit solution)
  !---------------------------------------------------------------------

  !-----------------------------------------------
  ! Subsidence source from layer above
  !               and 
  ! Evaporation of cloud within the layer
  !-----------------------------------------------
  IF (JK > NCLDTOP) THEN

    DO JL=KIDIA,KFDIA
      ZMF(JL)=MAX(0.0_JPRB,(PMFU(JL,JK)*ZADVW(JL)+PMFD(JL,JK)*ZADVWD(JL))*ZDTGDP(JL) )
!      ZMF(JL)=MAX(0.0_JPRB,(PMFU(JL,JK)+PMFD(JL,JK))*ZDTGDP(JL) )
      ZACUST(JL)=ZMF(JL)*ZANEWM1(JL)
    ENDDO

    DO JL=KIDIA,KFDIA
      ZLCUST(JL,NCLDQL) = ZMF(JL)*ZQXNM1(JL,NCLDQL)
      ZLCUST(JL,NCLDQI) = ZMF(JL)*ZQXNM1(JL,NCLDQI)
      ! record total flux for enthalpy budget:
      ZCONVSRCE(JL,NCLDQL) = ZCONVSRCE(JL,NCLDQL)+ZLCUST(JL,NCLDQL)
      ZCONVSRCE(JL,NCLDQI) = ZCONVSRCE(JL,NCLDQI)+ZLCUST(JL,NCLDQI)
    ENDDO

    ! Now have to work out how much liquid evaporates at arrival point 
    ! since there is no prognostic memory for in-cloud humidity, i.e. 
    ! we always assume cloud is saturated. 

    DO JL=KIDIA,KFDIA
      ZDTDP=ZRDCP*0.5_JPRB*(ZTP1(JL,JK-1)+ZTP1(JL,JK))/PAPH(JL,JK)
      ZDTFORC = ZDTDP*(PAP(JL,JK)-PAP(JL,JK-1))
      ![#Note: Diagnostic mixed phase should be replaced below]
      ZDQS(JL)=ZANEWM1(JL)*ZDTFORC*ZDQSMIXDT(JL)
    ENDDO

    ! Cloud liquid (NCLDQL)
    DO JL=KIDIA,KFDIA
      ZLFINAL=MAX(0.0_JPRB,ZLCUST(JL,NCLDQL)-ZDQS(JL)) !lim to zero
      ! no supersaturation allowed incloud ---V
      ZEVAP=MIN((ZLCUST(JL,NCLDQL)-ZLFINAL),ZEVAPLIMMIX(JL)) 
      ZLFINAL=ZLCUST(JL,NCLDQL)-ZEVAP 
      ZLFINALSUM(JL)=ZLFINALSUM(JL)+ZLFINAL ! sum 

      ZSOLQA(JL,NCLDQL,NCLDQL) = ZSOLQA(JL,NCLDQL,NCLDQL)+ZLCUST(JL,NCLDQL)
      ZSOLQA(JL,NCLDQV,NCLDQL) = ZSOLQA(JL,NCLDQV,NCLDQL)+ZEVAP
      ZSOLQA(JL,NCLDQL,NCLDQV) = ZSOLQA(JL,NCLDQL,NCLDQV)-ZEVAP
      ! Store cloud budget diagnostic
      ZBUDL(JL,4) = ZLCUST(JL,NCLDQL)*ZQTMST
      ZBUDL(JL,5) = -ZEVAP*ZQTMST
    ENDDO

    ! Cloud ice (NCLDQI)
    DO JL=KIDIA,KFDIA
      ZLFINAL=MAX(0.0_JPRB,ZLCUST(JL,NCLDQI)-ZDQS(JL)) !lim to zero
      ! no supersaturation allowed incloud ---V
      ZEVAP=MIN((ZLCUST(JL,NCLDQI)-ZLFINAL),ZEVAPLIMMIX(JL)) 
      ZLFINAL=ZLCUST(JL,NCLDQI)-ZEVAP 
      ZLFINALSUM(JL)=ZLFINALSUM(JL)+ZLFINAL ! sum 

      ZSOLQA(JL,NCLDQI,NCLDQI) = ZSOLQA(JL,NCLDQI,NCLDQI)+ZLCUST(JL,NCLDQI)
      ZSOLQA(JL,NCLDQV,NCLDQI) = ZSOLQA(JL,NCLDQV,NCLDQI)+ZEVAP
      ZSOLQA(JL,NCLDQI,NCLDQV) = ZSOLQA(JL,NCLDQI,NCLDQV)-ZEVAP
      ! Store cloud budget diagnostic
      ZBUDI(JL,4) = ZLCUST(JL,NCLDQI)*ZQTMST
      ZBUDI(JL,5) = -ZEVAP*ZQTMST 
    ENDDO
    
    !  Reset the cloud contribution if no cloud water survives to this level:
    DO JL=KIDIA,KFDIA
      IF (ZLFINALSUM(JL)<ZEPSEC) ZACUST(JL)=0.0_JPRB
      ZSOLAC(JL)=ZSOLAC(JL)+ZACUST(JL)
      ! Store cloud budget diagnostic
      ZBUDCC(JL,4) = ZACUST(JL)*ZQTMST
      !ZBUDCC(JL,5) isn't included as only reduced if cloud->zero
    ENDDO

  ENDIF ! on  JK>NCLDTOP

  !---------------------------------------------------------------------
  ! Subsidence sink of cloud to the layer below 
  ! (Implicit - re. CFL limit on convective mass flux)
  !---------------------------------------------------------------------

  DO JL=KIDIA,KFDIA

    IF(JK<KLEV) THEN

      ZMFDN=MAX(0.0_JPRB,(PMFU(JL,JK+1)*ZADVW(JL)+PMFD(JL,JK+1)*ZADVWD(JL))*ZDTGDP(JL) )
!      ZMFDN=MAX(0.0_JPRB,(PMFU(JL,JK+1)+PMFD(JL,JK+1))*ZDTGDP(JL) )
      ZSOLAB(JL)=ZSOLAB(JL)+ZMFDN
      ZSOLQB(JL,NCLDQL,NCLDQL)=ZSOLQB(JL,NCLDQL,NCLDQL)+ZMFDN
      ZSOLQB(JL,NCLDQI,NCLDQI)=ZSOLQB(JL,NCLDQI,NCLDQI)+ZMFDN

      ! Record sink for cloud budget and enthalpy budget diagnostics
      ZCONVSINK(JL,NCLDQL) = ZMFDN
      ZCONVSINK(JL,NCLDQI) = ZMFDN

    ENDIF

  ENDDO


  !======================================================================
  !
  !
  ! 3.4  EROSION OF CLOUDS BY TURBULENT MIXING
  !
  !
  !======================================================================
  ! NOTE: This process decreases the cloud area 
  !       but leaves the specific cloud water content
  !       *within clouds* unchanged
  !----------------------------------------------------------------------

 IF (ITURBEROSION == 1) THEN

  ! ------------------------------
  ! Define turbulent erosion rate
  ! ------------------------------
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA
    !original version (possibly perturbed by SPP)
    IF (LLRCLDIFF) THEN !Apply SPP perturbations
      ZLDIFDT(JL)=RCLDIFF*PTSPHY*EXP(ZMU_RCLDIFF+YSPP_CONFIG%CMPERT_RCLDIFF*PGP2DSPP(JL, IPRCLDIFF))
    ELSE
      ZLDIFDT(JL)=RCLDIFF*PTSPHY ! (unperturbed)
    ENDIF
    !Increase by factor of 5 for convective points
    IF(KTYPE(JL) > 0 .AND. PLUDE(JL,JK) > ZEPSEC) THEN
      IF(.NOT.(KTYPE(JL) >= 2)) &
       & ZLDIFDT(JL)=RCLDIFF_CONVI*ZLDIFDT(JL)  
    ENDIF
  ENDDO

  ! At the moment, works on mixed RH profile and partitioned ice/liq fraction
  ! so that it is similar to previous scheme
  ! Should apply RHw for liquid cloud and RHi for ice cloud separately 
  DO JL=KIDIA,KFDIA
    IF(ZLI(JL,JK) > ZEPSEC) THEN
      ! Calculate environmental humidity
!      ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSMIX(JL,JK))/&
!    &      MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))  
!      ZE=ZLDIFDT(JL)*MAX(ZQSMIX(JL,JK)-ZQE,0.0_JPRB)
      ZE=ZLDIFDT(JL)*MAX(ZQSMIX(JL,JK)-ZQX(JL,JK,NCLDQV),0.0_JPRB)
      ZLEROS=ZA(JL,JK)*ZE
      ZLEROS=MIN(ZLEROS,ZEVAPLIMMIX(JL))
      ZLEROS=MIN(ZLEROS,ZLI(JL,JK))
      ZAEROS=ZLEROS/ZLICLD(JL)  !if linear term

      ! Erosion is -ve LINEAR in L,A
      ZSOLAC(JL)=ZSOLAC(JL)-ZAEROS !linear

      ZSOLQA(JL,NCLDQV,NCLDQL) = ZSOLQA(JL,NCLDQV,NCLDQL)+ZLIQFRAC(JL,JK)*ZLEROS
      ZSOLQA(JL,NCLDQL,NCLDQV) = ZSOLQA(JL,NCLDQL,NCLDQV)-ZLIQFRAC(JL,JK)*ZLEROS
      ZSOLQA(JL,NCLDQV,NCLDQI) = ZSOLQA(JL,NCLDQV,NCLDQI)+ZICEFRAC(JL,JK)*ZLEROS
      ZSOLQA(JL,NCLDQI,NCLDQV) = ZSOLQA(JL,NCLDQI,NCLDQV)-ZICEFRAC(JL,JK)*ZLEROS

      ! Store cloud budget diagnostics if required
      ZBUDL(JL,7)  = -ZLIQFRAC(JL,JK)*ZLEROS*ZQTMST
      ZBUDI(JL,7)  = -ZICEFRAC(JL,JK)*ZLEROS*ZQTMST
      ZBUDCC(JL,7) = -ZAEROS*ZQTMST

    ENDIF
  ENDDO

 ELSEIF (ITURBEROSION == 2) THEN

  ! ------------------------------
  ! Define turbulent erosion rate
  ! Based on Morcrette
  ! Todo:
  ! To reduce calculations could remove duplication of some code for liq/ice
  ! Could also remove if tests as there is a limiter
  ! ZLICLD(JL) is LIQ+ICE/ZA
  ! Should this use ZCORQSLIQ(JL)?
  ! Have to use right subsaturation calculation otherwise duplicate evap
  ! ------------------------------
  DO JL=KIDIA,KFDIA

    ZLDIFDT(JL) = RCLDIFF*PTSPHY

    ! Increase erosion rate if convective
    !  KTYPE=1  Penetrative (deep) convection
    !  KTYPE=2  Shallow convection
    !  KTYPE=3  Mid-level convection
    ! Increase erosion for shallow/mid-level convection only
    IF(KTYPE(JL) >= 2 .AND. PLUDE(JL,JK) > ZEPSEC) THEN
       ZLDIFDT(JL) = RCLDIFF_CONVI*ZLDIFDT(JL)  
    ENDIF

  ENDDO

  ! Turbulent erosion of liquid
  DO JL=KIDIA,KFDIA
    ! If cloud liquid water is present
    IF (ZQX(JL,JK,NCLDQL) > ZEPSEC) THEN

      ! Calculate saturation deficit (wrt water)
      ZE = MAX(ZQSLIQ(JL,JK)-ZQX(JL,JK,NCLDQV),0.0_JPRB)

      ! Following Morcrette (2012)
      ZLEROS = 0.333_JPRB * ZA(JL,JK)*(1._JPRB-ZA(JL,JK))*ZLDIFDT(JL)*ZE
      
      ! Limiter taking account of evaporative cooling reducing saturation
      ZLEROS = MIN(ZLEROS,ZEVAPLIMLIQ(JL))
      
      ! Limiter to not remove more liquid than is present
      ZLEROS = MIN(ZLEROS,ZQX(JL,JK,NCLDQL))
      
      ! Cloud fraction decrease linearly proportional to liq+ice water decrease
      ZAEROS = ZLEROS/ZLICLD(JL)  !if linear term

      ! Erosion is -ve LINEAR in L,A
      ZSOLAC(JL) = ZSOLAC(JL)-ZAEROS !linear

      ! Update source/sink terms
      ZSOLQA(JL,NCLDQV,NCLDQL) = ZSOLQA(JL,NCLDQV,NCLDQL) + ZLEROS
      ZSOLQA(JL,NCLDQL,NCLDQV) = ZSOLQA(JL,NCLDQL,NCLDQV) - ZLEROS

      ! Store cloud budget diagnostics
      ZBUDL(JL,7)  = -ZLEROS*ZQTMST
      ZBUDCC(JL,7) = -ZAEROS*ZQTMST

    ENDIF
  ENDDO
  
  ! Turbulent erosion of ice
  DO JL=KIDIA,KFDIA
    ! If cloud ice is present
    IF (ZQX(JL,JK,NCLDQI) > ZEPSEC) THEN

      ! Calculate saturation deficit (wrt water)
      ZE = MAX(ZQSICE(JL,JK)-ZQX(JL,JK,NCLDQV),0.0_JPRB)

      ! Following Morcrette (2012)
      ZLEROS = 0.333_JPRB * ZA(JL,JK)*(1._JPRB-ZA(JL,JK))*ZLDIFDT(JL)*ZE
      
      ! Limiter taking account of evaporative cooling reducing saturation
      ZLEROS = MIN(ZLEROS,ZEVAPLIMICE(JL))
      
      ! Limiter to not remove more liquid than is present
      ZLEROS = MIN(ZLEROS,ZQX(JL,JK,NCLDQI))
      
      ! Cloud fraction decrease linearly proportional to liq+ice water decrease
      ZAEROS = ZLEROS/ZLICLD(JL)  !if linear term

      ! Erosion is -ve LINEAR in L,A
      ZSOLAC(JL) = ZSOLAC(JL)-ZAEROS !linear

      ! Update source/sink terms
      ZSOLQA(JL,NCLDQV,NCLDQI) = ZSOLQA(JL,NCLDQV,NCLDQI) + ZLEROS
      ZSOLQA(JL,NCLDQI,NCLDQV) = ZSOLQA(JL,NCLDQI,NCLDQV) - ZLEROS

      ! Store cloud budget diagnostics
      ZBUDI(JL,7)  = -ZLEROS*ZQTMST
      ZBUDCC(JL,7) = ZBUDCC(JL,7)-ZAEROS*ZQTMST

    ENDIF
  ENDDO
  
 ELSEIF (ITURBEROSION == 3) THEN

  ! ------------------------------
  ! Define turbulent erosion rate
  ! Based on Morcrette
  ! Todo:
  ! To reduce calculations could remove duplication of some code for liq/ice
  ! Could also remove if tests as there is a limiter
  ! ZLICLD(JL) is LIQ+ICE/ZA
  ! Should this use ZCORQSLIQ(JL)? No because it is the saturation deficit
  ! but do want to limit to max evaporation possible. Hence ZEVAPLIMLIQ
  ! Have to use right subsaturation calculation otherwise duplicate evap
  ! This version works on mixed phase
  ! ------------------------------
  DO JL=KIDIA,KFDIA

  !original version (possibly perturbed by SPP)
    IF (LLRCLDIFF) THEN
      ZLDIFDT(JL)=RCLDIFF*PTSPHY*EXP(ZMU_RCLDIFF+YSPP_CONFIG%CMPERT_RCLDIFF*PGP2DSPP(JL, IPRCLDIFF))
    ELSE
      ZLDIFDT(JL)=RCLDIFF*PTSPHY !original version (unperturbed)
    ENDIF

    ! Increase erosion rate if convective
    !  KTYPE=1  Penetrative (deep) convection
    !  KTYPE=2  Shallow convection
    !  KTYPE=3  Mid-level convection
    ! Increase erosion for shallow/mid-level convection only
    IF(KTYPE(JL) >= 2 .AND. PLUDE(JL,JK) > ZEPSEC) THEN
       ZLDIFDT(JL) = RCLDIFF_CONVI*ZLDIFDT(JL)  
    ENDIF

  ENDDO

  ! Turbulent erosion of liquid
  DO JL=KIDIA,KFDIA

    ! If cloud liquid water is present
    IF(ZLI(JL,JK) > ZEPSEC) THEN

      ! Calculate environmental humidity
      ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSMIX(JL,JK))/ &
    &      MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))

      ! Calculate saturation deficit (wrt mixed phase)
      ZE = MAX(ZQSMIX(JL,JK)-ZQE,0.0_JPRB)
            
      ! Following Morcrette (2012)
      ZLEROS = 0.333_JPRB * ZA(JL,JK)*(1._JPRB-ZA(JL,JK))*ZLDIFDT(JL)*ZE
      
      ! Limiter taking account of evaporative cooling reducing saturation
      ZLEROS = MIN(ZLEROS,ZEVAPLIMMIX(JL))
      
      ! Limiter to not remove more condensate than is present
      ZLEROS = MIN(ZLEROS,ZLI(JL,JK))
         
      ! Cloud fraction decrease linearly proportional to liq+ice water decrease
      ZAEROS = ZLEROS/ZLICLD(JL)  !if linear term
       
      ! Erosion is -ve LINEAR in L,A
      ZSOLAC(JL) = ZSOLAC(JL)-ZAEROS !linear

      ! Update source/sink terms
      ZSOLQA(JL,NCLDQV,NCLDQL) = ZSOLQA(JL,NCLDQV,NCLDQL)+ZLIQFRAC(JL,JK)*ZLEROS
      ZSOLQA(JL,NCLDQL,NCLDQV) = ZSOLQA(JL,NCLDQL,NCLDQV)-ZLIQFRAC(JL,JK)*ZLEROS
      ZSOLQA(JL,NCLDQV,NCLDQI) = ZSOLQA(JL,NCLDQV,NCLDQI)+ZICEFRAC(JL,JK)*ZLEROS
      ZSOLQA(JL,NCLDQI,NCLDQV) = ZSOLQA(JL,NCLDQI,NCLDQV)-ZICEFRAC(JL,JK)*ZLEROS

      ! Store cloud budget diagnostics if required
      ZBUDL(JL,7)  = -ZLIQFRAC(JL,JK)*ZLEROS*ZQTMST
      ZBUDI(JL,7)  = -ZICEFRAC(JL,JK)*ZLEROS*ZQTMST
      ZBUDCC(JL,7) = -ZAEROS*ZQTMST

    ENDIF
  ENDDO

 ENDIF ! on ITURBEROSION


  !======================================================================
  !
  !
  ! 3.5  CONDENSATION/EVAPORATION DUE TO DQSAT/DT
  !
  !
  !======================================================================
  !  calculate dqs/dt
  !  Note: For the separate prognostic Qi and Ql, one would ideally use
  !  Qsat/DT wrt liquid/Koop here, since the physics is that new clouds
  !  forms by liquid droplets [liq] or when aqueous aerosols [Koop] form.
  !  These would then instantaneous freeze if T<-38C or lead to ice growth 
  !  by deposition in warmer mixed phase clouds.  However, since we do 
  !  not have a separate prognostic equation for in-cloud humidity or a 
  !  statistical scheme approach in place, the depositional growth of ice 
  !  in the mixed phase can not be modelled and we resort to supersaturation  
  !  wrt ice instanteously converting to ice over one timestep 
  !  (see Tompkins et al. QJRMS 2007 for details)
  !  Thus for the initial implementation the diagnostic mixed phase is 
  !  retained for the moment, and the level of approximation noted.  
  !----------------------------------------------------------------------

  DO JL=KIDIA,KFDIA
    IF (LLQSATVERVEL) THEN  ! Compute SPP perturbation:
      ZPVERVEL(JL,JK)=PVERVEL(JL,JK)*EXP(ZMU_QSATVERVEL+ &
              & YSPP_CONFIG%CMPERT_QSATVERVEL*PGP2DSPP(JL,IPQSATVERVEL))
    ELSE
      ZPVERVEL(JL,JK)=PVERVEL(JL,JK)  ! (unperturbed)
    ENDIF
  ENDDO
  
  DO JL=KIDIA,KFDIA
    ZDTDP   = ZRDCP*ZTP1(JL,JK)/PAP(JL,JK)
    ZDPMXDT = ZDP(JL)*ZQTMST
    ZMFDN   = 0.0_JPRB
    IF(JK < KLEV) ZMFDN=PMFU(JL,JK+1)*ZADVW(JL)+PMFD(JL,JK+1)*ZADVWD(JL)
    ZWTOT   = ZPVERVEL(JL,JK)+(0.5_JPRB*RG*(PMFU(JL,JK)*ZADVW(JL)+PMFD(JL,JK)*ZADVWD(JL)+ZMFDN))
!    ZWTOT   = ZPVERVEL(JL,JK)+0.5_JPRB*RG*(PMFU(JL,JK)+PMFD(JL,JK)+ZMFDN)
    ZWTOT   = MIN(ZDPMXDT,MAX(-ZDPMXDT,ZWTOT))
    ZZZDT   = PHRSW(JL,JK)+PHRLW(JL,JK)
    ZDTDIAB = MIN(ZDPMXDT*ZDTDP,MAX(-ZDPMXDT*ZDTDP,ZZZDT))&
                    & *PTSPHY+RALFDCP*ZLDEFR(JL)  
! Note: ZLDEFR should be set to the difference between the mixed phase functions
! in the convection and cloud scheme, but this is not calculated, so is zero and
! the functions must be the same
    ZDTFORC = ZDTDP*ZWTOT*PTSPHY+ZDTDIAB
    ZQOLD(JL)   = ZQSMIX(JL,JK)
    ZTOLD(JL)   = ZTP1(JL,JK)
    ZTP1(JL,JK) = ZTP1(JL,JK)+ZDTFORC
    ZTP1(JL,JK) = MAX(ZTP1(JL,JK),160.0_JPRB)
    LLFLAG(JL)  = .TRUE.
  ENDDO

  IK=JK
  ICALL=5
  ZPAP(1:KLON) = PAP(KIDIA:(KIDIA+KLON-1),JK)
  CALL CUADJTQ &
   & (YRTHF, YRCST, YDEPHLI, KIDIA, KFDIA, KLON, KLEV, IK,&
   & ZPAP, ZTP1,    ZQSMIX,   LLFLAG,  ICALL)  

  DO JL=KIDIA,KFDIA
    ZDQS(JL)      = ZQSMIX(JL,JK)-ZQOLD(JL)
    ZQSMIX(JL,JK) = ZQOLD(JL)
    ZTP1(JL,JK)   = ZTOLD(JL)
  ENDDO

  !----------------------------------------------------------------------
  ! 3.5a  ZDQS(JL) > 0:  EVAPORATION OF CLOUDS
  ! ----------------------------------------------------------------------
  ! Assumes a delta function of cloud condensate, no change to cloud cover

  DO JL=KIDIA,KFDIA

   IF (ZDQS(JL) > 0.0_JPRB) THEN
    ZLEVAP = ZA(JL,JK)*MIN(ZDQS(JL),ZLICLD(JL))
    ZLEVAP = MIN(ZLEVAP,ZEVAPLIMMIX(JL))
    ZLEVAP = MIN(ZLEVAP,MAX(ZQSMIX(JL,JK)-ZQX(JL,JK,NCLDQV),0.0_JPRB))

    ! For first guess call
    ZLEVAPL(JL) = ZLIQFRAC(JL,JK)*ZLEVAP
    ZLEVAPI(JL) = ZICEFRAC(JL,JK)*ZLEVAP

    ZSOLQA(JL,NCLDQV,NCLDQL) = ZSOLQA(JL,NCLDQV,NCLDQL)+ZLIQFRAC(JL,JK)*ZLEVAP
    ZSOLQA(JL,NCLDQL,NCLDQV) = ZSOLQA(JL,NCLDQL,NCLDQV)-ZLIQFRAC(JL,JK)*ZLEVAP

    ZSOLQA(JL,NCLDQV,NCLDQI) = ZSOLQA(JL,NCLDQV,NCLDQI)+ZICEFRAC(JL,JK)*ZLEVAP
    ZSOLQA(JL,NCLDQI,NCLDQV) = ZSOLQA(JL,NCLDQI,NCLDQV)-ZICEFRAC(JL,JK)*ZLEVAP

    ! Store cloud budget diagnostics if required
    ZBUDL(JL,8)=-ZLIQFRAC(JL,JK)*ZLEVAP*ZQTMST
    ZBUDI(JL,8)=-ZICEFRAC(JL,JK)*ZLEVAP*ZQTMST

   ENDIF

  ENDDO

  !----------------------------------------------------------------------
  ! 3.5b ZDQS(JL) < 0: FORMATION OF CLOUDS
  !----------------------------------------------------------------------
  ! (1) Increase of cloud water in existing clouds
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA
    IF(ZA(JL,JK) > ZEPSEC.AND.ZDQS(JL) <= -RLMIN) THEN

      ZLCOND1(JL)=MAX(-ZDQS(JL),0.0_JPRB) !new limiter

!old limiter (significantly improves upper tropospheric humidity rms)
      IF(ZA(JL,JK) > 0.99_JPRB) THEN
        ZCOR=1.0_JPRB/(1.0_JPRB-RETV*ZQSMIX(JL,JK))
        ZCDMAX=(ZQX(JL,JK,NCLDQV)-ZQSMIX(JL,JK))/&
         & (1.0_JPRB+ZCOR*ZQSMIX(JL,JK)*FOEDEM(ZTP1(JL,JK)))  
      ELSE
        ZCDMAX=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSMIX(JL,JK))/ZA(JL,JK)
      ENDIF
      ZLCOND1(JL)=MAX(MIN(ZLCOND1(JL),ZCDMAX),0.0_JPRB)
! end old limiter
      
      ZLCOND1(JL)=ZA(JL,JK)*ZLCOND1(JL)
      IF(ZLCOND1(JL) < RLMIN) ZLCOND1(JL)=0.0_JPRB
      
      !-------------------------------------------------------------------------
      ! All increase goes into liquid unless so cold cloud homogeneously freezes
      ! Include new liquid formation in first guess value, otherwise liquid 
      ! remains at cold temperatures until next timestep.
      !-------------------------------------------------------------------------
      IF (ZTP1(JL,JK)>RTHOMO) THEN
        ZSOLQA(JL,NCLDQL,NCLDQV)=ZSOLQA(JL,NCLDQL,NCLDQV)+ZLCOND1(JL)
        ZSOLQA(JL,NCLDQV,NCLDQL)=ZSOLQA(JL,NCLDQV,NCLDQL)-ZLCOND1(JL)
        ZQXFG(JL,NCLDQL)=ZQXFG(JL,NCLDQL)+ZLCOND1(JL)
        ! Store cloud liquid diagnostic if required
        ZBUDL(JL,9)=ZLCOND1(JL)*ZQTMST
      ELSE
        ZSOLQA(JL,NCLDQI,NCLDQV)=ZSOLQA(JL,NCLDQI,NCLDQV)+ZLCOND1(JL)
        ZSOLQA(JL,NCLDQV,NCLDQI)=ZSOLQA(JL,NCLDQV,NCLDQI)-ZLCOND1(JL)
        ZQXFG(JL,NCLDQI)=ZQXFG(JL,NCLDQI)+ZLCOND1(JL)
        ! Store cloud ice diagnostic if required
        ZBUDI(JL,9)=ZLCOND1(JL)*ZQTMST
      ENDIF
    ENDIF
  ENDDO

  ! (2) Generation of new clouds (da/dt>0)
  
  DO JL=KIDIA,KFDIA

    IF(ZDQS(JL) <= -RLMIN .AND. ZA(JL,JK)<1.0_JPRB-ZEPSEC) THEN

      !---------------------------
      ! Critical relative humidity
      !---------------------------
      IF (LLRAMID) THEN !Apply SPP perturbations
        ZXRAMID= RAMID*EXP(ZMU_RAMID+YSPP_CONFIG%CMPERT_RAMID*PGP2DSPP(JL, IPRAMID))
        IF (YSPP_CONFIG%LRAMIDLIMIT1) THEN
          ZXRAMID= MIN(1.0_JPRB, ZXRAMID)
        ENDIF
      ELSE
        ZXRAMID= RAMID  ! (unperturbed)
      ENDIF
      
      ZRHC=ZXRAMID
      ZSIGK=PAP(JL,JK)/PAPH(JL,KLEV+1)
      ! Increase RHcrit to 1.0 towards the surface (eta>0.8)
      IF(ZSIGK > 0.8_JPRB) THEN
        ZRHC=ZXRAMID+(1.0_JPRB-ZXRAMID)*((ZSIGK-0.8_JPRB)/0.2_JPRB)**2
      ENDIF

! Commented out for CY37R1 to reduce humidity in high trop and strat
!      ! Increase RHcrit to 1.0 towards the tropopause (trop-0.2) and above
!      ZBOTT=ZTRPAUS(JL)+0.2_JPRB
!      IF(ZSIGK < ZBOTT) THEN
!        ZRHC=ZXRAMID+(1.0_JPRB-ZXRAMID)*MIN(((ZBOTT-ZSIGK)/0.2_JPRB)**2,1.0_JPRB)
!      ENDIF

      !---------------------------
      ! Supersaturation options
      !---------------------------      
      IF (NSSOPT==0) THEN 
        ! No scheme
        ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSICE(JL,JK))/&
            & MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))  
        ZQE=MAX(0.0_JPRB,ZQE)
      ELSEIF (NSSOPT==1) THEN 
        ! Tompkins 
        ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSICE(JL,JK))/&
            & MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))  
        ZQE=MAX(0.0_JPRB,ZQE)
      ELSEIF (NSSOPT==2) THEN 
        ! Lohmann and Karcher
        ZQE=ZQX(JL,JK,NCLDQV)  
      ELSEIF (NSSOPT==3) THEN 
        ! Gierens
        ZQE=ZQX(JL,JK,NCLDQV)+ZLI(JL,JK)
      ENDIF

      IF (ZTP1(JL,JK)>=RTT .OR. NSSOPT==0) THEN 
        ! No ice supersaturation allowed
        ZFAC=1.0_JPRB        
      ELSE
        ! Ice supersaturation
        ZFAC=ZFOKOOP(JL)
      ENDIF

      IF(ZQE >= ZRHC*ZQSICE(JL,JK)*ZFAC.AND.ZQE<ZQSICE(JL,JK)*ZFAC) THEN
        ! note: not **2 on 1-a term if ZQE is used. 
        ! Added correction term ZFAC to numerator 15/03/2010
        ZACOND=-(1.0_JPRB-ZA(JL,JK))*ZFAC*ZDQS(JL)/&
         &MAX(2.0_JPRB*(ZFAC*ZQSICE(JL,JK)-ZQE),ZEPSEC)

        ZACOND=MIN(ZACOND,1.0_JPRB-ZA(JL,JK))  !PUT THE LIMITER BACK

        ! Linear term:
        ! Added correction term ZFAC 15/03/2010
        ZLCOND2(JL)=-ZFAC*ZDQS(JL)*0.5_JPRB*ZACOND !mine linear

        ! new limiter formulation
        ZZDL=2.0_JPRB*(ZFAC*ZQSICE(JL,JK)-ZQE)/MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))
        ! Added correction term ZFAC 15/03/2010
        IF (ZFAC*ZDQS(JL)<-ZZDL) THEN
          ! ZLCONDLIM=(ZA(JL,JK)-1.0_JPRB)*ZDQS(JL)-ZQSICE(JL,JK)+ZQX(JL,JK,NCLDQV)
          ZLCONDLIM=(ZA(JL,JK)-1.0_JPRB)*ZFAC*ZDQS(JL)- &
     &               ZFAC*ZQSICE(JL,JK)+ZQX(JL,JK,NCLDQV)
          ZLCOND2(JL)=MIN(ZLCOND2(JL),ZLCONDLIM)
        ENDIF
        ZLCOND2(JL)=MAX(ZLCOND2(JL),0.0_JPRB)

        IF(ZLCOND2(JL) < RLMIN .OR. (1.0_JPRB-ZA(JL,JK))<ZEPSEC ) THEN
          ZLCOND2(JL) = 0.0_JPRB
          ZACOND      = 0.0_JPRB
        ENDIF
        IF(ZLCOND2(JL) == 0.0_JPRB) ZACOND=0.0_JPRB

        ! Large-scale generation is LINEAR in A and LINEAR in L
        ZSOLAC(JL) = ZSOLAC(JL)+ZACOND !linear
        
        ! Store cloud fraction diagnostic if required
        ZBUDCC(JL,10)=ZACOND*ZQTMST

        !------------------------------------------------------------------------
        ! All increase goes into liquid unless so cold cloud homogeneously freezes
        ! Include new liquid formation in first guess value, otherwise liquid 
        ! remains at cold temperatures until next timestep.
        !------------------------------------------------------------------------
        IF (ZTP1(JL,JK)>RTHOMO) THEN
          ZSOLQA(JL,NCLDQL,NCLDQV)=ZSOLQA(JL,NCLDQL,NCLDQV)+ZLCOND2(JL)
          ZSOLQA(JL,NCLDQV,NCLDQL)=ZSOLQA(JL,NCLDQV,NCLDQL)-ZLCOND2(JL)
          ZQXFG(JL,NCLDQL)=ZQXFG(JL,NCLDQL)+ZLCOND2(JL)
          ! Store cloud liquid diagnostic if required
          ZBUDL(JL,10)=ZLCOND2(JL)*ZQTMST
        ELSE ! homogeneous freezing
          ZSOLQA(JL,NCLDQI,NCLDQV)=ZSOLQA(JL,NCLDQI,NCLDQV)+ZLCOND2(JL)
          ZSOLQA(JL,NCLDQV,NCLDQI)=ZSOLQA(JL,NCLDQV,NCLDQI)-ZLCOND2(JL)
          ZQXFG(JL,NCLDQI)=ZQXFG(JL,NCLDQI)+ZLCOND2(JL)
          ! Store cloud ice diagnostic if required
          ZBUDI(JL,10)=ZLCOND2(JL)*ZQTMST
        ENDIF

      ENDIF
    ENDIF
  ENDDO


  !############################################################################
  !#                                                                          #
  !#                                                                          #
  !#                     4.  MICROPHYSICAL PROCESSES                          #
  !#                                                                          #
  !#                                                                          #
  !############################################################################


  !======================================================================
  !
  !
  ! 4.1 GROWTH OF ICE BY VAPOUR DEPOSITION 
  !
  !
  !======================================================================
  ! does not use the ice nuclei number from cloudaer.F90
  ! but rather a simple Meyers et al. 1992 form based on the 
  ! supersaturation and assuming clouds are saturated with 
  ! respect to liquid water (well mixed), (or Koop adjustment)
  ! Growth considered as sink of liquid water if present
  !----------------------------------------------------------------------

  !--------------------------------------------------------
  !-
  !- Ice deposition following Rotstayn et al. (2001)
  !-  (monodisperse ice particle size distribution)
  !-
  !--------------------------------------------------------
  IF (IDEPICE == 1) THEN

!DIR$ IVDEP  
  DO JL=KIDIA,KFDIA

    !--------------------------------------------------------------
    ! Calculate distance from cloud top 
    ! defined by cloudy layer below a layer with cloud frac <0.01
    ! ZDZ = ZDP(JL)/(ZRHO(JL)*RG)
    !--------------------------------------------------------------
      
    IF (ZA(JL,JK-1) < RCLDTOPCF .AND. ZA(JL,JK) >= RCLDTOPCF) THEN
      ZCLDTOPDIST(JL) = 0.0_JPRB
    ELSE
      ZCLDTOPDIST(JL) = ZCLDTOPDIST(JL) + ZDP(JL)/(ZRHO(JL)*RG)
    ENDIF

    !--------------------------------------------------------------
    ! Set subgrid overlap fraction of supercooled liquid and ice
    ! Reduce in shallow convection because assume SLW in active 
    ! updraught is less overlapped with ice in less active part
    !--------------------------------------------------------------
    ZOVERLAP_LIQICE = RCL_OVERLAPLIQICE
    
    IF (KTYPE(JL) > 0 .AND. PLUDE(JL,JK) > ZEPSEC) THEN
      ZOVERLAP_LIQICE = 0.1_JPRB
    ENDIF

    !--------------------------------------------------------------
    ! only treat depositional growth if liquid present. due to fact 
    ! that can not model ice growth from vapour without additional 
    ! in-cloud water vapour variable
    !--------------------------------------------------------------
    IF (ZTP1(JL,JK)<RTT .AND. ZQXFG(JL,NCLDQL)>RLMIN) THEN  ! T<273K

      ZVPICE=FOEEICE(ZTP1(JL,JK))*RV/RD
      ZVPLIQ=ZVPICE*ZFOKOOP(JL) 
      ZICENUCLEI(JL)=1000.0_JPRB*EXP(12.96_JPRB*(ZVPLIQ-ZVPICE)/ZVPLIQ-0.639_JPRB)

      !------------------------------------------------
      !   2.4e-2 is conductivity of air
      !   8.8 = 700**1/3 = density of ice to the third
      !------------------------------------------------
      ZADD=RLSTT*(RLSTT/(RV*ZTP1(JL,JK))-1.0_JPRB)/(2.4E-2_JPRB*ZTP1(JL,JK))
      ZBDD=RV*ZTP1(JL,JK)*PAP(JL,JK)/(2.21_JPRB*ZVPICE)
      ZCVDS=7.8_JPRB*(ZICENUCLEI(JL)/ZRHO(JL))**0.666_JPRB*(ZVPLIQ-ZVPICE) / &
         & (8.87_JPRB*(ZADD+ZBDD)*ZVPICE)

      !-----------------------------------------------------
      ! RICEINIT=1.E-12_JPRB is initial mass of ice particle
      !-----------------------------------------------------
      ZICE0=MAX(ZICECLD(JL), ZICENUCLEI(JL)*RICEINIT/ZRHO(JL))

      !------------------
      ! new value of ice:
      !------------------
      ZINEW=(0.666_JPRB*ZCVDS*PTSPHY+ZICE0**0.666_JPRB)**1.5_JPRB

      !---------------------------
      ! grid-mean deposition rate:
      !--------------------------- 
      ZDEPOS=MAX(ZOVERLAP_LIQICE*ZA(JL,JK)*(ZINEW-ZICE0),0.0_JPRB)

      !--------------------------------------------------------------------
      ! Limit deposition to liquid water amount
      ! If liquid is all frozen, ice would use up reservoir of water 
      ! vapour in excess of ice saturation mixing ratio - However this 
      ! can not be represented without a in-cloud humidity variable. Using 
      ! the grid-mean humidity would imply a large artificial horizontal 
      ! flux from the clear sky to the cloudy area. We thus rely on the 
      ! supersaturation check to clean up any remaining supersaturation
      !--------------------------------------------------------------------
      ZDEPOS=MIN(ZDEPOS,ZQXFG(JL,NCLDQL)) ! limit to liquid water amount
      
      !--------------------------------------------------------------------
      ! At top of cloud, reduce deposition rate near cloud top to account for
      ! small scale turbulent processes, limited ice nucleation and ice fallout 
      !--------------------------------------------------------------------
      ! Include dependence on ice nuclei concentration
      ! to increase deposition rate with decreasing temperatures 
      ZINFACTOR = MIN(ZICENUCLEI(JL)/15000._JPRB, 1.0_JPRB)
      ZDEPOS = ZDEPOS*MIN(ZINFACTOR + (1.0_JPRB-ZINFACTOR)* &
                  & (RDEPLIQREFRATE+ZCLDTOPDIST(JL)/RDEPLIQREFDEPTH),1.0_JPRB)

      !--------------
      ! add to matrix 
      !--------------
      ZSOLQA(JL,NCLDQI,NCLDQL)=ZSOLQA(JL,NCLDQI,NCLDQL)+ZDEPOS
      ZSOLQA(JL,NCLDQL,NCLDQI)=ZSOLQA(JL,NCLDQL,NCLDQI)-ZDEPOS
      ZQXFG(JL,NCLDQI)=ZQXFG(JL,NCLDQI)+ZDEPOS
      ZQXFG(JL,NCLDQL)=ZQXFG(JL,NCLDQL)-ZDEPOS
      ! Store cloud budget diagnostics if required
      ZBUDL(JL,11) = -ZDEPOS*ZQTMST
      ZBUDI(JL,11) = ZDEPOS*ZQTMST

    ENDIF
  ENDDO

  !--------------------------------------------------------
  !-
  !- Ice deposition assuming ice PSD
  !-
  !--------------------------------------------------------
  ELSEIF (IDEPICE == 2) THEN

    DO JL=KIDIA,KFDIA

      !--------------------------------------------------------------
      ! Calculate distance from cloud top 
      ! defined by cloudy layer below a layer with cloud frac <0.01
      ! ZDZ = ZDP(JL)/(ZRHO(JL)*RG)
      !--------------------------------------------------------------

      IF (ZA(JL,JK-1) < RCLDTOPCF .AND. ZA(JL,JK) >= RCLDTOPCF) THEN
        ZCLDTOPDIST(JL) = 0.0_JPRB
      ELSE
        ZCLDTOPDIST(JL) = ZCLDTOPDIST(JL) + ZDP(JL)/(ZRHO(JL)*RG)
      ENDIF

      !--------------------------------------------------------------
      ! only treat depositional growth if liquid present. due to fact 
      ! that can not model ice growth from vapour without additional 
      ! in-cloud water vapour variable
      !--------------------------------------------------------------
      IF (ZTP1(JL,JK)<RTT .AND. ZQXFG(JL,NCLDQL)>RLMIN) THEN  ! T<273K
      
        ZVPICE = FOEEICE(ZTP1(JL,JK))*RV/RD
        ZVPLIQ = ZVPICE*ZFOKOOP(JL) 
        ZICENUCLEI(JL)=1000.0_JPRB*EXP(12.96_JPRB*(ZVPLIQ-ZVPICE)/ZVPLIQ-0.639_JPRB)

        !-----------------------------------------------------
        ! RICEINIT=1.E-12_JPRB is initial mass of ice particle
        !-----------------------------------------------------
        ZICE0=MAX(ZICECLD(JL), ZICENUCLEI(JL)*RICEINIT/ZRHO(JL))
        
        ! Particle size distribution
        ZTCG    = 1.0_JPRB
        ZFACX1I = 1.0_JPRB

        ZAPLUSB   = RCL_APB1*ZVPICE-RCL_APB2*ZVPICE*ZTP1(JL,JK)+ &
       &             PAP(JL,JK)*RCL_APB3*ZTP1(JL,JK)**3._JPRB
        ZCORRFAC  = (1.0_JPRB/ZRHO(JL))**0.5_JPRB
        ZCORRFAC2 = ((ZTP1(JL,JK)/273.0_JPRB)**1.5_JPRB) &
       &             *(393.0_JPRB/(ZTP1(JL,JK)+120.0_JPRB))

        ZPR02  = ZRHO(JL)*ZICE0*RCL_CONST1I/(ZTCG*ZFACX1I)

        ZTERM1 = (ZVPLIQ-ZVPICE)*ZTP1(JL,JK)**2.0_JPRB*ZVPICE*ZCORRFAC2*ZTCG* &
       &          RCL_CONST2I*ZFACX1I/(ZRHO(JL)*ZAPLUSB*ZVPICE)
        ZTERM2 = 0.65_JPRB*RCL_CONST6I*ZPR02**RCL_CONST4I+RCL_CONST3I &
       &          *ZCORRFAC**0.5_JPRB*ZRHO(JL)**0.5_JPRB &
       &          *ZPR02**RCL_CONST5I/ZCORRFAC2**0.5_JPRB

        ZDEPOS = MAX(ZA(JL,JK)*ZTERM1*ZTERM2*PTSPHY,0.0_JPRB)

        !--------------------------------------------------------------------
        ! Limit deposition to liquid water amount
        ! If liquid is all frozen, ice would use up reservoir of water 
        ! vapour in excess of ice saturation mixing ratio - However this 
        ! can not be represented without a in-cloud humidity variable. Using 
        ! the grid-mean humidity would imply a large artificial horizontal 
        ! flux from the clear sky to the cloudy area. We thus rely on the 
        ! supersaturation check to clean up any remaining supersaturation
        !--------------------------------------------------------------------
        ZDEPOS=MIN(ZDEPOS,ZQXFG(JL,NCLDQL)) ! limit to liquid water amount

        !--------------------------------------------------------------------
        ! At top of cloud, reduce deposition rate near cloud top to account for
        ! small scale turbulent processes, limited ice nucleation and ice fallout 
        !--------------------------------------------------------------------
        ! Change to include dependence on ice nuclei concentration
        ! to increase deposition rate with decreasing temperatures 
        ZINFACTOR = MIN(ZICENUCLEI(JL)/15000._JPRB, 1.0_JPRB)
        ZDEPOS = ZDEPOS*MIN(ZINFACTOR + (1.0_JPRB-ZINFACTOR)* &
                    & (RDEPLIQREFRATE+ZCLDTOPDIST(JL)/RDEPLIQREFDEPTH),1.0_JPRB)

        !--------------
        ! add to matrix 
        !--------------
        ZSOLQA(JL,NCLDQI,NCLDQL) = ZSOLQA(JL,NCLDQI,NCLDQL)+ZDEPOS
        ZSOLQA(JL,NCLDQL,NCLDQI) = ZSOLQA(JL,NCLDQL,NCLDQI)-ZDEPOS
        ZQXFG(JL,NCLDQI) = ZQXFG(JL,NCLDQI)+ZDEPOS
        ZQXFG(JL,NCLDQL) = ZQXFG(JL,NCLDQL)-ZDEPOS
        ! Store cloud budget diagnostics if required
        ZBUDL(JL,11) = -ZDEPOS*ZQTMST
        ZBUDI(JL,11) = ZDEPOS*ZQTMST

      ENDIF
    ENDDO

  ENDIF ! on IDEPICE

  !--------------------------------------------------------
  !-
  !- Ice sublimation assuming ice PSD
  !-
  !--------------------------------------------------------
  IF (IEVAPICE == 1) THEN

    DO JL=KIDIA,KFDIA

      !-----------------------------------------------------------------------
      ! Calculate relative humidity limit for snow evaporation 
      !-----------------------------------------------------------------------
      ZZRH=RPRECRHMAX+(1.0_JPRB-RPRECRHMAX)*ZCOVPMAX(JL)/MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))
      ZZRH=MIN(MAX(ZZRH,RPRECRHMAX),1.0_JPRB)
!      ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSICE(JL,JK))/ &
!      & MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))  
      ZQE=ZQX(JL,JK,NCLDQV)  

      !---------------------------------------------
      ! humidity in moistest ZCOVPCLR part of domain
      !---------------------------------------------
!      ZQE=MAX(0.0_JPRB,MIN(ZQE,ZQSICE(JL,JK)))
      LLO1=ZA(JL,JK) > ZEPSEC .AND. ZQX(JL,JK,NCLDQI)>ZEPSEC &
     &  .AND. ZQE<ZQSICE(JL,JK)
!     &  .AND. ZQE<ZZRH*ZQSICE(JL,JK)

      IF(LLO1) THEN
             
        ! Calculate local in-cloud ice (kg/kg)
        ZICE0 = ZQX(JL,JK,NCLDQI)/ZA(JL,JK)
        ! Calculate saturation wrt ice
        ZVPICE = FOEEICE(ZTP1(JL,JK))*RV/RD
        
        ! Particle size distribution
        ZTCG    = 1.0_JPRB
        ZFACX1I = 1.0_JPRB

        ZAPLUSB   = RCL_APB1*ZVPICE-RCL_APB2*ZVPICE*ZTP1(JL,JK)+ &
       &             PAP(JL,JK)*RCL_APB3*ZTP1(JL,JK)**3._JPRB
        ZCORRFAC  = (1.0_JPRB/ZRHO(JL))**0.5_JPRB
        ZCORRFAC2 = ((ZTP1(JL,JK)/273.0_JPRB)**1.5_JPRB) &
       &             *(393.0_JPRB/(ZTP1(JL,JK)+120.0_JPRB))

        ZPR02  = ZRHO(JL)*ZICE0*RCL_CONST1I/(ZTCG*ZFACX1I)

        ZTERM1 = (ZQSICE(JL,JK)-ZQE)*ZTP1(JL,JK)**2.0_JPRB*ZVPICE*ZCORRFAC2*ZTCG* &
       &          RCL_CONST2I*ZFACX1I/(ZRHO(JL)*ZAPLUSB*ZVPICE)
        ZTERM2 = 0.65_JPRB*RCL_CONST6I*ZPR02**RCL_CONST4I+RCL_CONST3I &
       &          *ZCORRFAC**0.5_JPRB*ZRHO(JL)**0.5_JPRB &
       &          *ZPR02**RCL_CONST5I/ZCORRFAC2**0.5_JPRB

        ZDPEVAP = MAX(ZA(JL,JK)*ZTERM1*ZTERM2*PTSPHY,0.0_JPRB)

        !--------------------------------------------------------------------
        ! Limit sublimation to ice water amount
        !--------------------------------------------------------------------
        ZEVAP = MIN(ZDPEVAP,ZEVAPLIMICE(JL))
        ZEVAP = MIN(ZEVAP,ZQX(JL,JK,NCLDQI))

        !--------------
        ! add to matrix 
        !--------------
        ZSOLQA(JL,NCLDQV,NCLDQI) = ZSOLQA(JL,NCLDQV,NCLDQI)+ZEVAP
        ZSOLQA(JL,NCLDQI,NCLDQV) = ZSOLQA(JL,NCLDQI,NCLDQV)-ZEVAP
        ZQXFG(JL,NCLDQI) = ZQXFG(JL,NCLDQI)-ZEVAP
        ! Store cloud budget diagnostic
        ZBUDI(JL,8) = -ZEVAP*ZQTMST

        ! Decrease cloud amount using RKOOPTAU timescale
        ZFACI = PTSPHY/(2.0_JPRB*RKOOPTAU)
        ZSOLAC(JL) = ZSOLAC(JL)-ZA(JL,JK)*ZFACI

      ENDIF
    ENDDO

  ENDIF ! on IEVAPICE

  
  !----------------------------------
  ! revise in-cloud consensate amount
  !----------------------------------
  DO JL=KIDIA,KFDIA
    ZTMPA = 1.0_JPRB/MAX(ZA(JL,JK),0.01_JPRB)
    ZLIQCLD(JL) = ZQXFG(JL,NCLDQL)*ZTMPA
    ZICECLD(JL) = ZQXFG(JL,NCLDQI)*ZTMPA
    ZLICLD(JL)  = ZLIQCLD(JL)+ZICECLD(JL)
  ENDDO


  !======================================================================
  !
  !
  ! 4.2 SEDIMENTATION/FALLING OF HYDROMETEORS
  !
  !
  !======================================================================
  ! All hydrometeors fall downwards (can be carried upwards by dynamics)
  !----------------------------------------------------------------------
  
  ! Loop over sedimenting hydrometeors
  DO JM = 1,NCLV
    IF (LLFALL(JM)) THEN
      DO JL=KIDIA,KFDIA
        !------------------------
        ! source from layer above 
        !------------------------
        IF (JK > NCLDTOP) THEN
          ZFALLSRCE(JL,JM) = ZPFPLSX(JL,JK,JM)*ZDTGDP(JL) 
          ZSOLQA(JL,JM,JM) = ZSOLQA(JL,JM,JM)+ZFALLSRCE(JL,JM)
          ZQXFG(JL,JM)     = ZQXFG(JL,JM)+ZFALLSRCE(JL,JM)
          ! use first guess precip----------V
          ZQPRETOT(JL)     = ZQPRETOT(JL)+ZQXFG(JL,JM) 
          IF (JM == NCLDQI) THEN
            ZBUDI(JL,12)=ZFALLSRCE(JL,JM)*ZQTMST
          ENDIF
        ENDIF

        !-------------------------------------------------
        ! Calculate fall speeds
        !-------------------------------------------------

        ! Default to constant fall speed
        ZFALL = ZVQX(JM)

        !-----------
        ! Cloud ice 
        !-----------
        ! if aerosol effect then override 
        !  note that for T>233K this is the same as above.
        IF (LAERICESED .AND. JM == NCLDQI) THEN
          ZRE_ICE=PRE_ICE(JL,JK) 
          ! The exponent value is from 
          ! Morrison et al. JAS 2005 Appendix
          ZFALL = 0.002_JPRB*ZRE_ICE**1.0_JPRB
        ENDIF

        ! Option to modify as fn(p,T) Heymsfield and Iaquinta JAS 2000
        !  ZFALL = ZFALL*((PAP(JL,JK)*RICEHI1)**(-0.178_JPRB)) &
        !             &*((ZTP1(JL,JK)*RICEHI2)**(-0.394_JPRB))
    
        !---------------------------------
        ! Rain (if variable fallspeed on)
        !---------------------------------
        IF (JM == NCLDQR .AND. IVARFALL == 1) THEN

          ! Fallspeed air density correction 
          ZFALLCORR = (RDENSREF/ZRHO(JL))**0.4

          ! Rain water content (kg/kg) simple average between this layer and layer above
          ZTEMP = 0.5*(ZQX(JL,JK,NCLDQR)+ZQX(JL,JK-1,NCLDQR))
          
          IF (ZTEMP > ZEPSEC) THEN
            ! Slope of particle size distribution
            ZLAMBDA = (RCL_FAC1/(ZRHO(JL)*ZTEMP))**RCL_FAC2 
  
            ! Calculate fallspeed
            ZFALL = ZFALLCORR*RCL_CONST7R*ZLAMBDA**(-1._JPRB*RCL_DR)

          ENDIF
        
        ENDIF
  
        !-------------------------------------------------
        ! Calculate sink to next layer (implicit)
        !-------------------------------------------------
        ZFALLSINK(JL,JM)=ZDTGDP(JL)*ZRHO(JL)*ZFALL
        ! Cloud budget diagnostic stored at end as implicit
      
      ENDDO ! jl  
    ENDIF ! LLFALL
  ENDDO ! End loop over hydrometeor type


  !======================================================================
  !
  !
  ! 4.3 DEFINE PRECIPITATION GRIDBOX FRACTION
  !
  !
  !======================================================================
  ! Although precipitation (rain/snow) are prognostic variables
  ! precipitation fraction is diagnostic, so needs to be calculated 
  ! using the prognostic cloud cover working down through each grid 
  ! column every timestep. Maximum-random overlap is assumed
  ! Since precipitation may be advected into a column with no cloud above
  ! an arbitrary minimum coverage if precip>0 is assumed (RCOVPMIN).
  ! Since there is no memory of the clear sky precip fraction, 
  ! the precipitation cover ZCOVPTOT, which has the memory in the column,
  ! is reduced proportionally with the precip evaporation rate.
  !---------------------------------------------------------------

  DO JL=KIDIA,KFDIA
    IF (ZQPRETOT(JL) > ZEPSEC) THEN
      ZCOVPTOT(JL)   = 1.0_JPRB - ((1.0_JPRB-ZCOVPTOT(JL))* &
       &              (1.0_JPRB - MAX(ZA(JL,JK),ZA(JL,JK-1)))/ &
       &              (1.0_JPRB - MIN(ZA(JL,JK-1),1.0_JPRB-1.E-06_JPRB)) )  
      ZCOVPTOT(JL)   = MAX(ZCOVPTOT(JL),RCOVPMIN)
      ZCOVPCLR(JL)   = MAX(0.0_JPRB,ZCOVPTOT(JL)-ZA(JL,JK)) ! clear sky proportion
      ZRAINCLD(JL)   = 0.5*(ZQX(JL,JK,NCLDQR)+ZQX(JL,JK-1,NCLDQR))/ZCOVPTOT(JL)
      ZRAINCLDM1(JL) = 0.5*(ZQX(JL,JK,NCLDQR)+ZQXNM1(JL,NCLDQR))/ZCOVPTOT(JL)
      ZSNOWCLD(JL)   = 0.5*(ZQX(JL,JK,NCLDQS)+ZQX(JL,JK-1,NCLDQS))/ZCOVPTOT(JL)
      ZSNOWCLDM1(JL) = 0.5*(ZQX(JL,JK,NCLDQS)+ZQXNM1(JL,NCLDQS))/ZCOVPTOT(JL)
      ZCOVPMAX(JL)   = MAX(ZCOVPTOT(JL),ZCOVPMAX(JL))
    ELSE
      ZRAINCLD(JL)   = 0.0_JPRB ! no precip
      ZRAINCLDM1(JL) = 0.0_JPRB ! no precip
      ZSNOWCLD(JL)   = 0.0_JPRB ! no precip
      ZSNOWCLDM1(JL) = 0.0_JPRB ! no precip
      ZCOVPTOT(JL)   = 0.0_JPRB ! no flux - reset cover
      ZCOVPCLR(JL)   = 0.0_JPRB ! reset clear sky proportion 
      ZCOVPMAX(JL)   = 0.0_JPRB ! reset max cover for ZZRH calc 
    ENDIF
  ENDDO


  !======================================================================
  !
  !
  ! 4.4  AUTOCONVERSION OF ICE TO SNOW
  !
  !
  !======================================================================
  ! Formulation follows Lin et al. (1983)
  ! Sink of ice, source of snow
  ! Implicit in ice water content 
  !----------------------------------------------------------------------

  DO JL=KIDIA,KFDIA
 
    IF(ZTP1(JL,JK) <= RTT) THEN
      IF (ZICECLD(JL)>ZEPSEC) THEN

        ZZCO=PTSPHY*RSNOWLIN1*EXP(RSNOWLIN2*(ZTP1(JL,JK)-RTT))

        IF (LAERICEAUTO) THEN
          ZLCRIT=PICRIT_AER(JL,JK)
          ! 0.3 = N**0.333 with N=0.027 
          ZZCO=ZZCO*(RNICE/PNICE(JL,JK))**0.333_JPRB
        ELSE
          IF (LLRLCRITSNOW) THEN !Apply SPP perturbations
            ZLCRIT= RLCRITSNOW*EXP(ZMU_RLCRITSNOW+YSPP_CONFIG%CMPERT_RLCRITSNOW*PGP2DSPP(JL, IPRLCRITSNOW))
          ELSE
            ZLCRIT=RLCRITSNOW    ! (unperturbed)
          ENDIF
        ENDIF

        ZSNOWAUT(JL)=ZZCO*(1.0_JPRB-EXP(-(ZICECLD(JL)/ZLCRIT)**2))
        ZSOLQB(JL,NCLDQS,NCLDQI)=ZSOLQB(JL,NCLDQS,NCLDQI)+ZSNOWAUT(JL)

      ENDIF
    ENDIF 


  !======================================================================
  !
  !
  ! 4.5  AUTOCONVERSION OF LIQUID TO RAIN
  !
  !
  !======================================================================

   IF (ZLIQCLD(JL) > ZEPSEC) THEN

      !----------------------------------------
      ! Calculate inhomogeneity 
      !----------------------------------------
      IF (LCLOUD_INHOMOG .OR. YDEPHY%LRAD_CLOUD_INHOMOG) THEN

        ! Total water (vapour+liquid) in g/kg
        ZQTOT = MIN((ZQX(JL,JK,NCLDQV)+ZQX(JL,JK,NCLDQL))*1000._JPRB,30._JPRB)

        ! Representative grid box length (km)
        ! PGAW = normalised gaussian quadrature weight / no. longitude pts
        ZGRIDLEN = 2*RA*SQRT(RPI*PGAW(JL))*0.001_JPRB

        ! Correlation between cloud and rain
        ZCLDRAINCORR = 1._JPRB-0.8_JPRB*ZA(JL,JK)

        !Maike's new parameters for modified Boutle parameterization
        ZPHIP1=.2_JPRB+.01_JPRB*ZQTOT+.0027_JPRB*ZQTOT**2-.00008_JPRB*ZQTOT**3
        ZPHIP2 = ZPHIP1-0.2_JPRB
        ZPHIP3 = 0.123_JPRB*EXP(-ZPHIP1*0.55_JPRB)
        !cut off very low cloud fraction at 0.1  to avoid very low FSD values
        ZANEW2(JL,JK)=MAX(ZA(JL,JK),0.1_JPRB) ! cut off very low cloud fraction to avoid

        ! Autoconversion enhancement factor from Boutle et al. 2013
        ZPHIC = (ZGRIDLEN*ZANEW2(JL,JK))**(1._JPRB/3._JPRB)* &
     &      ((ZPHIP3*ZGRIDLEN*ZANEW2(JL,JK))**1.5_JPRB + 3._JPRB*ZPHIP3)**(-0.17_JPRB)

        ! Fractional standard deviation for cloud condensate
        ZFRACSDC = (ZPHIP1 - ZPHIP2*ZANEW2(JL,JK))*ZPHIC

        ! Unique value when essentially no cloud edges
        IF (ZA(JL,JK) > 0.95_JPRB) ZFRACSDC = 0.17_JPRB*ZPHIC

        ! multiply by 1D-to-2D variability enhancement factor
        ZLFSD(JL,JK)=ZR12*ZFRACSDC
       ELSE
          ! Default FSD value
          ZLFSD(JL,JK)=YDERAD%RCLOUD_FRAC_STD
       ENDIF 


    !--------------------------------------------------------
    !-
    !- Warm-rain process follow Sundqvist (1989)
    !-
    !--------------------------------------------------------
    ! Implicit in liquid water content 
    
    IF (IWARMRAIN == 1) THEN

      ZZCO=RKCONV*PTSPHY

      IF (LAERLIQAUTOLSP) THEN
        ZLCRIT=PLCRIT_AER(JL,JK)
        ! 0.3 = N**0.333 with N=125 cm-3 
        ZZCO=ZZCO*(RCCN/PCCN(JL,JK))**0.333_JPRB
      ELSE
        ! Modify autoconversion threshold dependent on: 
        !  land (polluted, high CCN, smaller droplets, higher threshold)
        !  sea  (clean, low CCN, larger droplets, lower threshold)
        IF (LLRCLCRIT) THEN  !Apply SPP perturbations
          IF (PLSM(JL) > 0.5_JPRB) THEN
             ! perturbed land value of RCLCRIT
            ZLCRIT = RCLCRIT_LAND*EXP(ZMU_RCLCRITL+YSPP_CONFIG%CMPERT_RCLCRIT_LAND*PGP2DSPP(JL, IPRCLCRIT))
          ELSE
            ! perturbed ocean value of RCLCRIT
            ZLCRIT = RCLCRIT_SEA *EXP(ZMU_RCLCRITS+YSPP_CONFIG%CMPERT_RCLCRIT_SEA*PGP2DSPP(JL, IPRCLCRIT))
          ENDIF
        ELSE
          IF (PLSM(JL) > 0.5_JPRB) THEN
            ZLCRIT = RCLCRIT_LAND ! land  (unperturbed)
          ELSE
            ZLCRIT = RCLCRIT_SEA  ! ocean (unperturbed)
          ENDIF
        ENDIF
      ENDIF 

      !------------------------------------------------------------------
      ! Parameters for cloud collection by rain and snow.
      ! Note that with new prognostic variable it is now possible 
      ! to REPLACE this with an explicit collection parametrization
      !------------------------------------------------------------------   
      ZPRECIP=(ZPFPLSX(JL,JK,NCLDQS)+ZPFPLSX(JL,JK,NCLDQR))/MAX(ZEPSEC,ZCOVPTOT(JL))
      ZCFPR=1.0_JPRB + RPRC1*SQRT(MAX(ZPRECIP,0.0_JPRB))
!      ZCFPR=1.0_JPRB + RPRC1*SQRT(MAX(ZPRECIP,0.0_JPRB))* &
!       &ZCOVPTOT(JL)/(MAX(ZA(JL,JK),ZEPSEC))

      IF (LAERLIQCOLL) THEN 
        ! 5.0 = N**0.333 with N=125 cm-3 
        ZCFPR=ZCFPR*(RCCN/PCCN(JL,JK))**0.333_JPRB
      ENDIF

      ZZCO=ZZCO*ZCFPR
      ZLCRIT=ZLCRIT/MAX(ZCFPR,ZEPSEC)
  
      IF(ZLIQCLD(JL)/ZLCRIT < 20.0_JPRB )THEN ! Security for exp for some compilers
        ZRAINAUT(JL)=ZZCO*(1.0_JPRB-EXP(-(ZLIQCLD(JL)/ZLCRIT)**2))
      ELSE
        ZRAINAUT(JL)=ZZCO
      ENDIF

      ! rain freezes instantly
      IF(ZTP1(JL,JK) <= RTT) THEN
        ZSOLQB(JL,NCLDQS,NCLDQL)=ZSOLQB(JL,NCLDQS,NCLDQL)+ZRAINAUT(JL)
      ELSE
        ZSOLQB(JL,NCLDQR,NCLDQL)=ZSOLQB(JL,NCLDQR,NCLDQL)+ZRAINAUT(JL)
      ENDIF

    !--------------------------------------------------------
    !-
    !- Warm-rain process follow Khairoutdinov and Kogan (2000)
    !-
    !--------------------------------------------------------
    ELSEIF (IWARMRAIN == 2 .OR. IWARMRAIN == 3) THEN

      IF (LLRCLCRIT) THEN  !Apply SPP perturbations
        IF (PLSM(JL) > 0.5_JPRB) THEN
          ZCONST = RCL_KK_CLOUD_NUM_LAND
          ! perturbed land value of RCLCRIT
          ZLCRIT = RCLCRIT_LAND*EXP(ZMU_RCLCRITL+YSPP_CONFIG%CMPERT_RCLCRIT_LAND*PGP2DSPP(JL, IPRCLCRIT))
        ELSE
          ZCONST = RCL_KK_CLOUD_NUM_SEA
          ! perturbed ocean value of RCLCRIT
          ZLCRIT = RCLCRIT_SEA *EXP(ZMU_RCLCRITS+YSPP_CONFIG%CMPERT_RCLCRIT_SEA*PGP2DSPP(JL, IPRCLCRIT))
        ENDIF
      ELSE
        IF (PLSM(JL) > 0.5_JPRB) THEN ! land  (unperturbed)
          ZCONST = RCL_KK_CLOUD_NUM_LAND
          ZLCRIT = RCLCRIT_LAND
        ELSE                          ! ocean (unperturbed)
          ZCONST = RCL_KK_CLOUD_NUM_SEA
          ZLCRIT = RCLCRIT_SEA
        ENDIF
      ENDIF
 
      !----------------------------------------
      ! Calculate inhomogeneity 
      !----------------------------------------
      IF (LCLOUD_INHOMOG) THEN

         ZFRACSDC=ZLFSD(JL,JK)

         ZEAUT = (1._JPRB+ZFRACSDC**2._JPRB)**((RCL_KKBAUQ-1._JPRB) &
              &       *RCL_KKBAUQ/2._JPRB)

        ! Accretion enhancement factor from Boutle et al. 2013
        ZPHIR = (ZGRIDLEN*ZA(JL,JK))**(1._JPRB/3._JPRB)* &
              &       ((0.11_JPRB*ZGRIDLEN*ZA(JL,JK))**1.14_JPRB + 1._JPRB)**(-0.22_JPRB)
         
        ! Fractional standard deviation for rain condensate
         ZFRACSDR = (1.1_JPRB - 0.8_JPRB*ZA(JL,JK))*ZPHIR

        ! Unique value when essentially no cloud edges
         IF (ZA(JL,JK) > 0.95_JPRB) ZFRACSDR = 0.3_JPRB*ZPHIR

         ZEACC = &
     &    (1._JPRB+ZFRACSDC**2._JPRB)**((RCL_KKBAC-1._JPRB)*RCL_KKBAC/2._JPRB) &
     &   *(1._JPRB+ZFRACSDR**2._JPRB)**((RCL_KKBAC-1._JPRB)*RCL_KKBAC/2._JPRB) &
     &   * EXP(ZCLDRAINCORR*RCL_KKBAC*RCL_KKBAC* &
     &   SQRT(LOG(1._JPRB+ZFRACSDC**2._JPRB)*LOG(1._JPRB+ZFRACSDR**2._JPRB)))

       ELSE

        ! Simple constant multiplier to take account of inhomogeneity
        ZEAUT = RCL_INHOMOGAUT
        ZEACC = RCL_INHOMOGACC

       ENDIF
        
      !-------------------------------------------------------------------------
      ! Calculate autoconversion of cloud liquid droplets to rain
      ! Calculate accretion of cloud liquid droplets by rain
      !-------------------------------------------------------------------------

      !----------------------
      ! Explicit formulation
      !----------------------
      IF (IWARMRAIN == 2) THEN       
        ! Explicit formulation

        ! Autoconversion
        ZRAINAUT(JL)  = ZEAUT*ZA(JL,JK)*PTSPHY* &
     &                  RCL_KKAAU * ZLIQCLD(JL)**RCL_KKBAUQ * ZCONST**RCL_KKBAUN

        ZRAINAUT(JL) = MIN(ZRAINAUT(JL),ZQXFG(JL,NCLDQL))
        IF (ZRAINAUT(JL) < ZEPSEC) ZRAINAUT(JL) = 0.0_JPRB


        ! Accretion
        ZRAINACC(JL) = ZEACC*ZA(JL,JK)*PTSPHY* &
     &                 RCL_KKAAC * (ZLIQCLD(JL)*ZRAINCLD(JL))**RCL_KKBAC

        ZRAINACC(JL) = MIN(ZRAINACC(JL),ZQXFG(JL,NCLDQL))
        IF (ZRAINACC(JL) < ZEPSEC) ZRAINACC(JL) = 0.0_JPRB

        ! If temperature < 0, then autoconversion produces snow rather than rain
        ! Explicit
        IF(ZTP1(JL,JK) <= RTT) THEN
          ZSOLQA(JL,NCLDQS,NCLDQL)=ZSOLQA(JL,NCLDQS,NCLDQL)+ZRAINAUT(JL)
          ZSOLQA(JL,NCLDQS,NCLDQL)=ZSOLQA(JL,NCLDQS,NCLDQL)+ZRAINACC(JL)
          ZSOLQA(JL,NCLDQL,NCLDQS)=ZSOLQA(JL,NCLDQL,NCLDQS)-ZRAINAUT(JL)
          ZSOLQA(JL,NCLDQL,NCLDQS)=ZSOLQA(JL,NCLDQL,NCLDQS)-ZRAINACC(JL)
          ! Store cloud budget diagnostics
          ZBUDL(JL,12) = -ZRAINAUT(JL)*ZQTMST
          ZBUDL(JL,13) = -ZRAINACC(JL)*ZQTMST
        ELSE
          ZSOLQA(JL,NCLDQR,NCLDQL)=ZSOLQA(JL,NCLDQR,NCLDQL)+ZRAINAUT(JL)
          ZSOLQA(JL,NCLDQR,NCLDQL)=ZSOLQA(JL,NCLDQR,NCLDQL)+ZRAINACC(JL)
          ZSOLQA(JL,NCLDQL,NCLDQR)=ZSOLQA(JL,NCLDQL,NCLDQR)-ZRAINAUT(JL)
          ZSOLQA(JL,NCLDQL,NCLDQR)=ZSOLQA(JL,NCLDQL,NCLDQR)-ZRAINACC(JL)
          ! Store cloud budget diagnostics
          ZBUDL(JL,14) = -ZRAINAUT(JL)*ZQTMST
          ZBUDL(JL,15) = -ZRAINACC(JL)*ZQTMST
        ENDIF

      
      !----------------------
      ! Implicit formulation
      !----------------------
      ELSEIF (IWARMRAIN == 3) THEN
        
        ! (zqxfg taken out and zliqcld=zqxfg/za, so multiply by 1/za
        !  and za cancels out) 
        ZRAINAUT(JL) = ZEAUT*PTSPHY* &
     &                 RCL_KKAAU * ZLIQCLD(JL)**(RCL_KKBAUQ-1.0_JPRB) &
     &                 * ZCONST**RCL_KKBAUN

        ! If temperature < 0, then autoconversion produces snow rather than rain
        ! Implicit
        IF(ZTP1(JL,JK) <= RTT) THEN
          ZSOLQB(JL,NCLDQS,NCLDQL)=ZSOLQB(JL,NCLDQS,NCLDQL)+ZRAINAUT(JL)
        ELSE
          ZSOLQB(JL,NCLDQR,NCLDQL)=ZSOLQB(JL,NCLDQR,NCLDQL)+ZRAINAUT(JL)
        ENDIF

        IF (ZRAINCLD(JL) > ZEPSEC) THEN

          IF (IRAINACC == 1) THEN

            ! Khairoutdinov and Kogan
            ! (zqxfg taken out and zliqcld=zqxfg/za, so multiply by 1/za
            !  and za cancels out) 
            ZRAINACC(JL) = ZEACC*PTSPHY* &
       &                   RCL_KKAAC * ZLIQCLD(JL)**(RCL_KKBAC-1.0_JPRB) &
       &                   * ZRAINCLD(JL)**RCL_KKBAC
          ELSE

            ! Sweep out
            ! (zqxfg taken out and zliqcld=zqxfg/za, so multiply by 1/za
            !  and za cancels out) - implicit in lwc
          
            ! Fallspeed air density correction 
            ZFALLCORR = (RDENSREF/ZRHO(JL))**0.4

            ! Slope of particle size distribution
            ZLAMBDA = (RCL_FAC1/(ZRHO(JL)*ZRAINCLD(JL)))**RCL_FAC2 
                    
            ! Calculate accretion term
            ! Factor of liq water taken out because implicit
            ZRAINACC(JL) = ZEACC*PTSPHY*RCL_EFF_RACW &
                         & *RCL_CONST9R*ZFALLCORR/(ZLAMBDA**RCL_CONST10R)

            ! Limit rain accretion term - needed?
            ZRAINACC(JL)=MIN(ZRAINACC(JL),1.0_JPRB)
          
          ENDIF
          
        ENDIF

        ! If temperature < 0, then autoconversion produces snow rather than rain
        ! Implicit
        IF(ZTP1(JL,JK) <= RTT) THEN
          ZSOLQB(JL,NCLDQS,NCLDQL)=ZSOLQB(JL,NCLDQS,NCLDQL)+ZRAINACC(JL)
        ELSE
          ZSOLQB(JL,NCLDQR,NCLDQL)=ZSOLQB(JL,NCLDQR,NCLDQL)+ZRAINACC(JL)
        ENDIF

      ENDIF ! on IWARMRAIN = 2 or 3
    
    ENDIF ! on IWARMRAIN

   ENDIF ! on ZLIQCLD > ZEPSEC
  ENDDO


  !======================================================================
  !
  !
  ! 4.6 RIMING - COLLECTION OF CLOUD LIQUID DROPS BY SNOW AND ICE
  !
  !
  !======================================================================
  ! Only active if T<0degC and supercooled liquid water is present
  ! AND if not Sundquist autoconversion (as this includes riming)
  !----------------------------------------------------------------------
  
  IF (IWARMRAIN > 1) THEN

  DO JL=KIDIA,KFDIA
    
    ZSNOWRIME(JL) = 0.0_JPRB
     
    IF(ZTP1(JL,JK) <= RTT .AND. ZLIQCLD(JL)>ZEPSEC) THEN

      ! Fallspeed air density correction 
      ZFALLCORR = (RDENSREF/ZRHO(JL))**0.4

      !------------------------------------------------------------------
      ! Riming of snow by cloud water - implicit in lwc
      !------------------------------------------------------------------
      IF (ZSNOWCLD(JL)>ZEPSEC .AND. ZCOVPTOT(JL)>0.01_JPRB) THEN

        ! Calculate riming term
        ! Factor of liq water taken out because implicit
        ZSNOWRIME(JL) = RCL_EFFRIME*ZCOVPTOT(JL)*PTSPHY*RCL_CONST7S*ZFALLCORR &
     &                  *(ZRHO(JL)*ZSNOWCLD(JL)*RCL_CONST1S)**RCL_CONST8S

        ! Limit snow riming term
        ZSNOWRIME(JL)=MIN(ZSNOWRIME(JL),1.0_JPRB)

        ZSOLQB(JL,NCLDQS,NCLDQL) = ZSOLQB(JL,NCLDQS,NCLDQL) + ZSNOWRIME(JL)

      ENDIF

      !------------------------------------------------------------------
      ! Riming of ice by cloud water - implicit in lwc
      ! NOT YET ACTIVE
      !------------------------------------------------------------------
!      IF (ZICECLD(JL)>ZEPSEC .AND. ZA(JL,JK)>0.01_JPRB) THEN
!
!        ! Calculate riming term
!        ! Factor of liq water taken out because implicit
!        ZSNOWRIME(JL) = ZA(JL,JK)*PTSPHY*RCL_CONST7S*ZFALLCORR &
!     &                  *(ZRHO(JL)*ZICECLD(JL)*RCL_CONST1S)**RCL_CONST8S
!
!        ! Limit ice riming term
!        ZSNOWRIME(JL)=MIN(ZSNOWRIME(JL),1.0_JPRB)
!
!        ZSOLQB(JL,NCLDQI,NCLDQL) = ZSOLQB(JL,NCLDQI,NCLDQL) + ZSNOWRIME(JL)
!
!      ENDIF
    ENDIF
  ENDDO
  
  ENDIF ! on IWARMRAIN > 1


  !======================================================================
  !
  !
  ! 4.7  MELTING OF SNOW AND ICE
  !
  !
  !======================================================================
  ! With implicit solver this also has to treat snow or ice
  ! precipitating from the level above... i.e. local ice AND flux.
  ! in situ ice and snow: could arise from LS advection or warming
  ! falling ice and snow: arrives by precipitation process
  !----------------------------------------------------------------------
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA

    ZICETOT(JL)=ZQXFG(JL,NCLDQI)+ZQXFG(JL,NCLDQS)
    ZMELTMAX(JL) = 0.0_JPRB

    ! If there are frozen hydrometeors present and dry-bulb temperature > 0degC
    IF(ZICETOT(JL) > ZEPSEC .AND. ZTP1(JL,JK) > RTT) THEN

      ! Calculate subsaturation
      ZSUBSAT = MAX(ZQSICE(JL,JK)-ZQX(JL,JK,NCLDQV),0.0_JPRB)
      
      ! Calculate difference between dry-bulb (ZTP1) and the temperature 
      ! at which the wet-bulb=0degC (RTT-ZSUBSAT*....) using an approx.
      ! Melting only occurs if the wet-bulb temperature >0
      ! i.e. warming of ice particle due to melting > cooling 
      ! due to evaporation.
      ZTDMTW0 = ZTP1(JL,JK)-RTT-ZSUBSAT* &
                & (ZTW1+ZTW2*(PAP(JL,JK)-ZTW3)-ZTW4*(ZTP1(JL,JK)-ZTW5))
      ! Not implicit yet... 
      ! Ensure ZCONS1 is positive so that ZMELTMAX=0 if ZTDMTW0<0
      ZCONS1 = ABS(PTSPHY*(1.0_JPRB+0.5_JPRB*ZTDMTW0)/RTAUMEL)
      ZMELTMAX(JL) = MAX(ZTDMTW0*ZCONS1*ZRLDCP,0.0_JPRB)
    ENDIF
  ENDDO

  ! Loop over frozen hydrometeors (ice, snow)
  DO JM=1,NCLV
   IF (IPHASE(JM) == 2) THEN
    JN = IMELT(JM)
    DO JL=KIDIA,KFDIA
      IF(ZMELTMAX(JL)>ZEPSEC .AND. ZICETOT(JL)>ZEPSEC) THEN
        ! Apply melting in same proportion as frozen hydrometeor fractions 
        ZALFA = ZQXFG(JL,JM)/ZICETOT(JL)
        ZMELT = MIN(ZQXFG(JL,JM),ZALFA*ZMELTMAX(JL))
        ! needed in first guess
        ! This implies that zqpretot has to be recalculated below
        ! since is not conserved here if ice falls and liquid doesn't
        ZQXFG(JL,JM)     = ZQXFG(JL,JM)-ZMELT
        ZQXFG(JL,JN)     = ZQXFG(JL,JN)+ZMELT
        ZSOLQA(JL,JN,JM) = ZSOLQA(JL,JN,JM)+ZMELT
        ZSOLQA(JL,JM,JN) = ZSOLQA(JL,JM,JN)-ZMELT
        IF (JM==NCLDQI) ZBUDI(JL,15) = -ZMELT*ZQTMST
        IF (JM==NCLDQI) ZBUDL(JL,17) =  ZMELT*ZQTMST
        IF (JM==NCLDQS) ZBUDI(JL,16) = -ZMELT*ZQTMST
      ENDIF
    ENDDO
   ENDIF
  ENDDO

  
  !======================================================================
  !
  !
  ! 4.8  FREEZING OF RAIN
  !
  !
  !======================================================================
  ! Rain drop freezing rate based on Bigg(1953) and Wisner(1972)
  ! Rain -> Snow
  !----------------------------------------------------------------------
  
!DEC$ IVDEP
  DO JL=KIDIA,KFDIA 

    ! If rain present
    IF (ZQX(JL,JK,NCLDQR) > ZEPSEC) THEN

      IF (ZTP1(JL,JK) <= RTT .AND. ZTP1(JL,JK-1) > RTT) THEN
        ! Base of melting layer/top of refreezing layer so
        ! store rain/snow fraction for precip type diagnosis
        ! If mostly rain, then supercooled rain slow to freeze
        ! otherwise faster to freeze (snow or ice pellets)
        ZQPRETOT(JL) = MAX(ZQX(JL,JK,NCLDQS)+ZQX(JL,JK,NCLDQR),ZEPSEC)
        PRAINFRAC_TOPRFZ(JL) = ZQX(JL,JK,NCLDQR)/ZQPRETOT(JL)
        IF (PRAINFRAC_TOPRFZ(JL) > 0.8) THEN 
          LLRAINLIQ(JL) = .TRUE.
        ELSE
          LLRAINLIQ(JL) = .FALSE.
        ENDIF
      ENDIF
    
      ! If temperature less than zero
      IF (ZTP1(JL,JK) < RTT) THEN

        IF (LLRAINLIQ(JL)) THEN 

          ! Majority of raindrops completely melted
          ! Refreezing is by slow heterogeneous freezing
          
          ! Slope of rain particle size distribution
          ZLAMBDA = (RCL_FAC1/(ZRHO(JL)*ZQX(JL,JK,NCLDQR)))**RCL_FAC2

          ! Calculate freezing rate based on Bigg(1953) and Wisner(1972)
          ZTEMP = RCL_FZRAB * (ZTP1(JL,JK)-RTT)
          ZFRZ  = PTSPHY * (RCL_CONST5R/ZRHO(JL)) * (EXP(ZTEMP)-1._JPRB) &
                  & * ZLAMBDA**RCL_CONST6R
          ZFRZMAX(JL) = MAX(ZFRZ,0.0_JPRB)

        ELSE

          ! Majority of raindrops only partially melted 
          ! Refreeze with a shorter timescale (reverse of melting...for now)
          
          ZCONS1 = ABS(PTSPHY*(1.0_JPRB+0.5_JPRB*(RTT-ZTP1(JL,JK)))/RTAUMEL)
          ZFRZMAX(JL) = MAX((RTT-ZTP1(JL,JK))*ZCONS1*ZRLDCP,0.0_JPRB)

        ENDIF

        IF(ZFRZMAX(JL)>ZEPSEC) THEN
          ZFRZ = MIN(ZQX(JL,JK,NCLDQR),ZFRZMAX(JL))
          ZSOLQA(JL,NCLDQS,NCLDQR) = ZSOLQA(JL,NCLDQS,NCLDQR)+ZFRZ
          ZSOLQA(JL,NCLDQR,NCLDQS) = ZSOLQA(JL,NCLDQR,NCLDQS)-ZFRZ
          ZBUDL(JL,18) = ZFRZ*ZQTMST
        ENDIF
      ENDIF

    ENDIF

  ENDDO


  !======================================================================
  !
  !
  ! 4.9   FREEZING OF CLOUD LIQUID 
  !
  !
  !======================================================================
  ! All liquid cloud drops assumed to freeze instantaneously to ice crystals
  ! below the homogeneous freezing temperature (-38degC)
  ! Liquid -> Ice
  !----------------------------------------------------------------------
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA 
    ! not implicit yet... 
    ZFRZMAX(JL)=MAX((RTHOMO-ZTP1(JL,JK))*ZRLDCP,0.0_JPRB)
  ENDDO

  DO JL=KIDIA,KFDIA
    IF(ZFRZMAX(JL)>ZEPSEC .AND. ZQXFG(JL,NCLDQL)>ZEPSEC) THEN
      ZFRZ = MIN(ZQXFG(JL,NCLDQL),ZFRZMAX(JL))
      ZSOLQA(JL,NCLDQI,NCLDQL) = ZSOLQA(JL,NCLDQI,NCLDQL)+ZFRZ
      ZSOLQA(JL,NCLDQL,NCLDQI) = ZSOLQA(JL,NCLDQL,NCLDQI)-ZFRZ
      ZBUDL(JL,19) = -ZFRZ*ZQTMST
    ENDIF
  ENDDO


  !======================================================================
  !
  !
  ! 4.10  EVAPORATION OF RAIN
  !
  !
  !======================================================================
  ! Rain -> Vapour
  !----------------------------------------------------------------------

  DO JL=KIDIA,KFDIA

    !-----------------------------------------------------------------------
    ! Calculate relative humidity limit for rain evaporation 
    ! to avoid cloud formation and saturation of the grid box
    !-----------------------------------------------------------------------
    ! Limit RH for rain evaporation dependent on precipitation fraction 
    ZZRH=RPRECRHMAX+(1.0_JPRB-RPRECRHMAX)*ZCOVPMAX(JL)/MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))
    ZZRH=MIN(MAX(ZZRH,RPRECRHMAX),1.0_JPRB)

    ! Critical relative humidity
    IF (LLRAMID) THEN !Apply SPP perturbations
      ZXRAMID= RAMID*EXP(ZMU_RAMID+YSPP_CONFIG%CMPERT_RAMID*PGP2DSPP(JL, IPRAMID))
      IF (YSPP_CONFIG%LRAMIDLIMIT1) THEN
        ZXRAMID= MIN(1.0_JPRB, ZXRAMID)
      ENDIF
    ELSE
      ZXRAMID= RAMID ! (unperturbed)
    ENDIF
    ZRHC=ZXRAMID
    ZSIGK=PAP(JL,JK)/PAPH(JL,KLEV+1)
    ! Increase RHcrit to 1.0 towards the surface (eta>0.8)
    IF(ZSIGK > 0.8_JPRB) THEN
      ZRHC=ZXRAMID+(1.0_JPRB-ZXRAMID)*((ZSIGK-0.8_JPRB)/0.2_JPRB)**2
    ENDIF
    
    ! Limit evaporation to RHcrit threshold
    ZZRH = MIN(ZRHC,ZZRH)
    ! Limit further to not go above 90%
    ZZRH = MIN(0.9_JPRB,ZZRH)

  
    ZQE=MAX(0.0_JPRB,MIN(ZQX(JL,JK,NCLDQV),ZQSLIQ(JL,JK)))

    ! If there is precipitation in clear sky, and there is rain present 
    ! and the humidity is less than the threshold defined above then... 
    LLO1=ZCOVPCLR(JL)>ZEPSEC .AND. &
       & ZRAINCLDM1(JL)>ZEPSEC .AND. & 
       & ZQE<ZZRH*ZQSLIQ(JL,JK)

    IF(LLO1) THEN

      !-------------------------------------------
      ! Evaporation
      !-------------------------------------------
      ! Calculate local precipitation (kg/kg)
      ZPRECLR = ZRAINCLDM1(JL)

      ! Fallspeed air density correction 
      ZFALLCORR = (RDENSREF/ZRHO(JL))**0.4

      ! Saturation vapour pressure with respect to liquid phase
      ZESATLIQ = RV/RD*FOEELIQ(ZTP1(JL,JK))

      ! Slope of particle size distribution
      ZLAMBDA = (RCL_FAC1/(ZRHO(JL)*ZPRECLR))**RCL_FAC2 ! ZPRECLR=kg/kg

      ZEVAP_DENOM = RCL_CDENOM1*ZESATLIQ - RCL_CDENOM2*ZTP1(JL,JK)*ZESATLIQ &
              & + RCL_CDENOM3*ZTP1(JL,JK)**3._JPRB*PAP(JL,JK)

      ! Temperature dependent conductivity
      ZCORR2= (ZTP1(JL,JK)/273._JPRB)**1.5_JPRB*393._JPRB/(ZTP1(JL,JK)+120._JPRB)
      ZKA = RCL_KA273*ZCORR2

      ZSUBSAT = MAX(ZZRH*ZQSLIQ(JL,JK)-ZQE,0.0_JPRB)

      ZBETA = (0.15_JPRB/ZQSLIQ(JL,JK))*ZTP1(JL,JK)**2._JPRB*ZESATLIQ* &
     & RCL_CONST1R*(ZCORR2/ZEVAP_DENOM)*(0.78_JPRB/(ZLAMBDA**RCL_CONST4R)+ &
     & RCL_CONST2R*(ZRHO(JL)*ZFALLCORR)**0.5_JPRB/ &
     & (ZCORR2**0.5_JPRB*ZLAMBDA**RCL_CONST3R))
     
      ZDENOM  = 1.0_JPRB+ZBETA*PTSPHY*ZCORQSLIQ(JL)
      ZDPEVAP = ZCOVPCLR(JL)*ZBETA*PTSPHY*ZSUBSAT/ZDENOM

      !Apply SPP perturbations
      IF (LLRAINEVAP) THEN
        ZDPEVAP = ZDPEVAP*EXP(ZMU_RAINEVAP+YSPP_CONFIG%CMPERT_RAINEVAP*PGP2DSPP(JL, IPRAINEVAP))
      ENDIF

      !---------------------------------------------------------
      ! Add evaporation term to explicit sink.
      ! this has to be explicit since if treated in the implicit
      ! term evaporation can not reduce rain to zero and model
      ! produces small amounts of rainfall everywhere. 
      !---------------------------------------------------------
      
      ! Limit rain evaporation
      ZEVAP = MIN(ZDPEVAP,ZQXFG(JL,NCLDQR))

      ZSOLQA(JL,NCLDQV,NCLDQR) = ZSOLQA(JL,NCLDQV,NCLDQR)+ZEVAP
      ZSOLQA(JL,NCLDQR,NCLDQV) = ZSOLQA(JL,NCLDQR,NCLDQV)-ZEVAP

      ZBUDL(JL,20) = -ZEVAP*ZQTMST

      !-------------------------------------------------------------
      ! Reduce the total precip coverage proportional to evaporation
      ! to mimic the previous scheme which had a diagnostic
      ! 2-flux treatment, abandoned due to the new prognostic precip
      !-------------------------------------------------------------
      ! Comment out reduction of precipitation fraction with evaporation 
      !ZCOVPTOT(JL) = MAX(RCOVPMIN,ZCOVPTOT(JL)-MAX(0.0_JPRB, &
      ! &            (ZCOVPTOT(JL)-ZA(JL,JK))*ZEVAP/ZQXFG(JL,NCLDQR)))

      ! Update fg field 
      ZQXFG(JL,NCLDQR) = ZQXFG(JL,NCLDQR)-ZEVAP
    
    ENDIF
  ENDDO

  
  !======================================================================
  !
  !
  ! 4.11  SUBLIMATION OF SNOW
  !
  !
  !======================================================================
  ! Snow -> Vapour
  !----------------------------------------------------------------------

 IF (IEVAPSNOW == 1) THEN
  
  DO JL=KIDIA,KFDIA
    ZZRH=RPRECRHMAX+(1.0_JPRB-RPRECRHMAX)*ZCOVPMAX(JL)/MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))
    ZZRH=MIN(MAX(ZZRH,RPRECRHMAX),1.0_JPRB)
    ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSICE(JL,JK))/ &
    & MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))  

    !---------------------------------------------
    ! humidity in moistest ZCOVPCLR part of domain
    !---------------------------------------------
    ZQE=MAX(0.0_JPRB,MIN(ZQE,ZQSICE(JL,JK)))
    LLO1=ZCOVPCLR(JL)>ZEPSEC .AND. &
       & ZQXFG(JL,NCLDQS)>ZEPSEC .AND. &
       & ZQE<ZZRH*ZQSICE(JL,JK)

    IF(LLO1) THEN
      ! note: zpreclr is a rain flux a
      ZPRECLR=ZQXFG(JL,NCLDQS)*ZCOVPCLR(JL)/ &
       & SIGN(MAX(ABS(ZCOVPTOT(JL)*ZDTGDP(JL)),ZEPSILON),ZCOVPTOT(JL)*ZDTGDP(JL))

      !--------------------------------------
      ! actual microphysics formula in zbeta
      !--------------------------------------

      ZBETA1=SQRT(PAP(JL,JK)/ &
       & PAPH(JL,KLEV+1))/RVRFACTOR*ZPRECLR/ &
       & MAX(ZCOVPCLR(JL),ZEPSEC)

      ZBETA=RG*RPECONS*(ZBETA1)**0.5777_JPRB  

      ZDENOM=1.0_JPRB+ZBETA*PTSPHY*ZCORQSICE(JL)
      ZDPR = ZCOVPCLR(JL)*ZBETA*(ZQSICE(JL,JK)-ZQE)/ZDENOM*ZDP(JL)*ZRG_R
      ZDPEVAP=ZDPR*ZDTGDP(JL)

      !Apply SPP perturbations
      IF (LLSNOWSUBLIM) THEN
        ZDPEVAP = ZDPEVAP*EXP(ZMU_SNOWSUBLIM+YSPP_CONFIG%CMPERT_SNOWSUBLIM*PGP2DSPP(JL, IPSNOWSUBLIM))
      ENDIF

      !---------------------------------------------------------
      ! add evaporation term to explicit sink.
      ! this has to be explicit since if treated in the implicit
      ! term evaporation can not reduce snow to zero and model
      ! produces small amounts of snowfall everywhere. 
      !---------------------------------------------------------
      
      ! Evaporate snow
      ZEVAP = MIN(ZDPEVAP,ZQXFG(JL,NCLDQS))

      ZSOLQA(JL,NCLDQV,NCLDQS) = ZSOLQA(JL,NCLDQV,NCLDQS)+ZEVAP
      ZSOLQA(JL,NCLDQS,NCLDQV) = ZSOLQA(JL,NCLDQS,NCLDQV)-ZEVAP
      ZBUDI(JL,17) = -ZEVAP*ZQTMST
      
      !-------------------------------------------------------------
      ! Reduce the total precip coverage proportional to evaporation
      ! to mimic the previous scheme which had a diagnostic
      ! 2-flux treatment, abandoned due to the new prognostic precip
      !-------------------------------------------------------------
      ZCOVPTOT(JL) = MAX(RCOVPMIN,ZCOVPTOT(JL)-MAX(0.0_JPRB, &
     &              (ZCOVPTOT(JL)-ZA(JL,JK))*ZEVAP/ZQXFG(JL,NCLDQS)))
      
      !Update first guess field
      ZQXFG(JL,NCLDQS) = ZQXFG(JL,NCLDQS)-ZEVAP

    ENDIF
  ENDDO

  !---------------------------------------------------------
  ELSEIF (IEVAPSNOW == 2) THEN

 
   DO JL=KIDIA,KFDIA

    !-----------------------------------------------------------------------
    ! Calculate relative humidity limit for snow evaporation 
    !-----------------------------------------------------------------------
    ZZRH=RPRECRHMAX+(1.0_JPRB-RPRECRHMAX)*ZCOVPMAX(JL)/MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))
    ZZRH=MIN(MAX(ZZRH,RPRECRHMAX),1.0_JPRB)

    ZQE=(ZQX(JL,JK,NCLDQV)-ZA(JL,JK)*ZQSICE(JL,JK))/ &
    & MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))  
     
    !---------------------------------------------
    ! humidity in moistest ZCOVPCLR part of domain
    !---------------------------------------------
    ZQE=MAX(0.0_JPRB,MIN(ZQE,ZQSICE(JL,JK)))

    LLO1=ZCOVPCLR(JL)>ZEPSEC .AND. &
       & ZSNOWCLDM1(JL)>ZEPSEC .AND. & 
       & ZQE<ZZRH*ZQSICE(JL,JK)

    IF(LLO1) THEN
      
      ! Calculate local precipitation (kg/kg)
      ZPRECLR = ZSNOWCLDM1(JL)
     
      ! Saturation vapour pressure with respect to ice phase
      ZVPICE = RV/RD*FOEEICE(ZTP1(JL,JK))

      ! Particle size distribution
      ! ZTCG increases Ni with colder temperatures - essentially a 
      ! Fletcher or Meyers scheme? 
      ZTCG=1.0_JPRB !v1 EXP(RCL_X3I*(273.15_JPRB-ZTP1(JL,JK))/8.18_JPRB)
      ! ZFACX1I modification is based on Andrew Barrett's results
      ZFACX1S = 1.0_JPRB !v1 (ZICE0/1.E-5_JPRB)**0.627_JPRB

      ZAPLUSB   = RCL_APB1*ZVPICE-RCL_APB2*ZVPICE*ZTP1(JL,JK)+ &
     &             PAP(JL,JK)*RCL_APB3*ZTP1(JL,JK)**3
      ZCORRFAC  = (1.0/ZRHO(JL))**0.5
      ZCORRFAC2 = ((ZTP1(JL,JK)/273.0)**1.5)*(393.0/(ZTP1(JL,JK)+120.0))

      ZPR02 = ZRHO(JL)*ZPRECLR*RCL_CONST1S/(ZTCG*ZFACX1S)

      ZTERM1 = (ZQSICE(JL,JK)-ZQE)*ZTP1(JL,JK)**2*ZVPICE*ZCORRFAC2*ZTCG* &
     &          RCL_CONST2S*ZFACX1S/(ZRHO(JL)*ZAPLUSB*ZQSICE(JL,JK))
      ZTERM2 = 0.65*RCL_CONST6S*ZPR02**RCL_CONST4S+RCL_CONST3S*ZCORRFAC**0.5 &
     &          *ZRHO(JL)**0.5*ZPR02**RCL_CONST5S/ZCORRFAC2**0.5

      ZDPEVAP = MAX(ZCOVPCLR(JL)*ZTERM1*ZTERM2*PTSPHY,0.0_JPRB)
 
      !--------------------------------------------------------------------
      ! Limit evaporation to snow amount
      !--------------------------------------------------------------------
      ZEVAP = MIN(ZDPEVAP,ZEVAPLIMICE(JL))
      ZEVAP = MIN(ZEVAP,ZQXFG(JL,NCLDQS))
            
      ZSOLQA(JL,NCLDQV,NCLDQS) = ZSOLQA(JL,NCLDQV,NCLDQS)+ZEVAP
      ZSOLQA(JL,NCLDQS,NCLDQV) = ZSOLQA(JL,NCLDQS,NCLDQV)-ZEVAP
      ZBUDI(JL,17) = -ZEVAP*ZQTMST
      
      !-------------------------------------------------------------
      ! Reduce the total precip coverage proportional to evaporation
      ! to mimic the previous scheme which had a diagnostic
      ! 2-flux treatment, abandoned due to the new prognostic precip
      !-------------------------------------------------------------
      ZCOVPTOT(JL) = MAX(RCOVPMIN,ZCOVPTOT(JL)-MAX(0.0_JPRB, &
     &              (ZCOVPTOT(JL)-ZA(JL,JK))*ZEVAP/ZQXFG(JL,NCLDQS)))
      
      !Update first guess field
      ZQXFG(JL,NCLDQS) = ZQXFG(JL,NCLDQS)-ZEVAP

    ENDIF    
  ENDDO
     
ENDIF ! on IEVAPSNOW

  !--------------------------------------
  ! Evaporate small precipitation amounts
  !--------------------------------------
  DO JM=1,NCLV
   IF (LLFALL(JM)) THEN 
!DIR$ IVDEP
    DO JL=KIDIA,KFDIA
      IF (ZQXFG(JL,JM)<RLMIN) THEN
        ZSOLQA(JL,NCLDQV,JM) = ZSOLQA(JL,NCLDQV,JM)+ZQXFG(JL,JM)
        ZSOLQA(JL,JM,NCLDQV) = ZSOLQA(JL,JM,NCLDQV)-ZQXFG(JL,JM)
      ENDIF
    ENDDO
   ENDIF
  ENDDO

  
  !######################################################################
  !
  !            5.  *** SOLVERS FOR A AND L ***
  !
  ! Use an implicit solution rather than exact solution.
  ! Solver is forward in time, upstream difference for advection.
  !######################################################################

  !======================================================================
  !
  ! 5.1 Solver for cloud cover
  !
  !======================================================================
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA
    ZANEW=(ZA(JL,JK)+ZSOLAC(JL))/(1.0_JPRB+ZSOLAB(JL))
    ZANEW=MIN(ZANEW,1.0_JPRB)
    IF (ZANEW<RAMIN) ZANEW=0.0_JPRB
    ZDA(JL)=ZANEW-ZAORIG(JL,JK)
    !---------------------------------
    ! variables needed for next level
    !---------------------------------
    ZANEWM1(JL)=ZANEW
  ENDDO

  !======================================================================
  !
  ! 5.2 Solver for the microphysics
  !
  !======================================================================

  !--------------------------------------------------------------
  ! 5.2.1 Truncate explicit sinks to avoid negatives 
  ! Note: Species are treated in the order in which they run out
  ! since the clipping will alter the balance for the other vars
  !--------------------------------------------------------------

  !----------------------------
  ! compute sink terms
  !----------------------------
  DO JM=1,NCLV
     DO JL=KIDIA,KFDIA
         ZSINKSUM(JL)=0.0_JPRB
     ENDDO
     DO JN=1,NCLV
        DO JL=KIDIA,KFDIA
           ZSINKSUM(JL)=ZSINKSUM(JL)-ZSOLQA(JL,JM,JN) ! +ve total is bad
        ENDDO
     ENDDO

     !----------------------------------------------
     ! calculate overshoot and scaling factor
     ! if ZSINKSUM is -ve, no overshoot and ZRATIO=1
     !----------------------------------------------
     DO JL=KIDIA,KFDIA
        ZMAX=MAX(ZQX(JL,JK,JM),ZEPSEC)
        ZRAT=MAX(ZSINKSUM(JL),ZMAX)
        ZRATIO(JL,JM)=ZMAX/ZRAT
     ENDDO
  ENDDO

  !--------------------------------------------------------
  ! sort zratio to find out which species run out first
  !--------------------------------------------------------
  DO JL=KIDIA,KFDIA
     IORDV(1)=1

     DO JM=2,NCLV
        ! Make room to move ZRV(JM) to its final place
        ! in the sorted sequence 
        ! ZRATIO(JL,IORDV(1)) ... ZRATIO(JL,IORDV(JM-1))
        DO JN=JM-1,1,-1
           IF (ZRATIO(JL,IORDV(JN))<=ZRATIO(JL,JM)) EXIT
           IORDV(JN+1)=IORDV(JN)
        ENDDO
        
        IORDV(JN+1)=JM
     ENDDO

     IORDER(JL,1:NCLV)=IORDV(1:NCLV)
  ENDDO

  !----------------
  ! recalculate sum
  !----------------
  DO JM=1,NCLV
     DO JL=KIDIA,KFDIA
        ZSINKSUM(JL)=0.0_JPRB
     ENDDO

     DO JN=1,NCLV
        DO JL=KIDIA,KFDIA
           JO=IORDER(JL,JM)
           LLINDEX3(JL,JN)=ZSOLQA(JL,JO,JN)<0.0_JPRB
           ZSINKSUM(JL)=ZSINKSUM(JL)-ZSOLQA(JL,JO,JN) ! +ve total is bad
        ENDDO
     ENDDO
     !---------------------------
     ! recalculate scaling factor
     !---------------------------
     DO JL=KIDIA,KFDIA
        JO=IORDER(JL,JM)
        ZMM=MAX(ZQX(JL,JK,JO),ZEPSEC)
        ZRR=MAX(ZSINKSUM(JL),ZMM)
        ZRATIO(JL,1)=ZMM/ZRR
     ENDDO
     !------
     ! scale
     !------
     DO JL=KIDIA,KFDIA
        JO=IORDER(JL,JM)
        ZZRATIO=ZRATIO(JL,1)
        !DIR$ IVDEP
        !DIR$ PREFERVECTOR
        DO JN=1,NCLV
           IF (LLINDEX3(JL,JN)) THEN
              ZSOLQA(JL,JO,JN)=ZSOLQA(JL,JO,JN)*ZZRATIO
              ZSOLQA(JL,JN,JO)=ZSOLQA(JL,JN,JO)*ZZRATIO
           ENDIF
        ENDDO
     ENDDO
  ENDDO


  !--------------------------------------------------------------
  ! 5.2.2 Solver
  !------------------------

  !------------------------
  ! set the LHS of equation  
  !------------------------
  DO JM=1,NCLV
     DO JN=1,NCLV
        !----------------------------------------------
        ! diagonals: microphysical sink terms+transport
        !----------------------------------------------
        IF (JN==JM) THEN
           DO JL=KIDIA,KFDIA
              ZQLHS(JL,JN,JM)=1.0_JPRB + ZFALLSINK(JL,JM)
           ENDDO
           !$OMP SIMD PRIVATE(JO)
           DO JL=KIDIA,KFDIA
!DIR$ UNROLL
              DO JO=1,NCLV
                 ZQLHS(JL,JN,JM)=ZQLHS(JL,JN,JM) + ZSOLQB(JL,JO,JN)
              ENDDO
           ENDDO
           !------------------------------------------
           ! non-diagonals: microphysical source terms
           !------------------------------------------
        ELSE
           DO JL=KIDIA,KFDIA
              ZQLHS(JL,JN,JM)= -ZSOLQB(JL,JN,JM) ! here is the delta T - missing from doc.
           ENDDO
        ENDIF
     ENDDO
  ENDDO

  !------------------------
  ! set the RHS of equation  
  !------------------------
  DO JM=1,NCLV
    DO JL=KIDIA,KFDIA
      !---------------------------------
      ! sum the explicit source and sink
      !---------------------------------
      ZEXPLICIT=0.0_JPRB
      DO JN=1,NCLV
        ZEXPLICIT=ZEXPLICIT+ZSOLQA(JL,JM,JN) ! sum over middle index
      ENDDO
      ZQXN(JL,JM)=ZQX(JL,JK,JM)+ZEXPLICIT
    ENDDO
  ENDDO

  !-----------------------------------
  ! *** solve by LU decomposition: ***
  !-----------------------------------
  ! Note: This fast way of solving NCLVxNCLV system
  !       assumes a good behaviour (i.e. non-zero diagonal
  !       terms with comparable orders) of the matrix stored
  !       in ZQLHS. For the moment this is the case but
  !       be aware to preserve it when doing eventual 
  !       modifications.

  ! Non pivoting recursive factorization 
  DO JN = 1, NCLV-1  ! number of steps
    DO JM = JN+1,NCLV ! row index
      ZQLHS(KIDIA:KFDIA,JM,JN)=ZQLHS(KIDIA:KFDIA,JM,JN) &
       &                     / ZQLHS(KIDIA:KFDIA,JN,JN)
      DO IK=JN+1,NCLV ! column index
        DO JL=KIDIA,KFDIA
          ZQLHS(JL,JM,IK)=ZQLHS(JL,JM,IK)-ZQLHS(JL,JM,JN)*ZQLHS(JL,JN,IK)
        ENDDO
      ENDDO
    ENDDO
  ENDDO        

  ! Backsubstitution 
  !  step 1 
  DO JN=2,NCLV
    DO JM = 1,JN-1
      ZQXN(KIDIA:KFDIA,JN)=ZQXN(KIDIA:KFDIA,JN)-ZQLHS(KIDIA:KFDIA,JN,JM) &
       &  *ZQXN(KIDIA:KFDIA,JM)
    ENDDO
  ENDDO
  !  step 2
  ZQXN(KIDIA:KFDIA,NCLV)=ZQXN(KIDIA:KFDIA,NCLV)/ZQLHS(KIDIA:KFDIA,NCLV,NCLV)
  DO JN=NCLV-1,1,-1
    DO JM = JN+1,NCLV
      ZQXN(KIDIA:KFDIA,JN)=ZQXN(KIDIA:KFDIA,JN)-ZQLHS(KIDIA:KFDIA,JN,JM) &
       &  *ZQXN(KIDIA:KFDIA,JM)
    ENDDO
    ZQXN(KIDIA:KFDIA,JN)=ZQXN(KIDIA:KFDIA,JN)/ZQLHS(KIDIA:KFDIA,JN,JN)
  ENDDO

  ! Ensure no small values (including negatives) remain in cloud variables nor
  ! precipitation rates.
  ! Evaporate l,i,r,s to water vapour. Latent heating taken into account below
  DO JN=1,NCLV-1
    DO JL=KIDIA,KFDIA
      IF (ZQXN(JL,JN) < ZEPSEC) THEN
        ZQXN(JL,NCLDQV) = ZQXN(JL,NCLDQV)+ZQXN(JL,JN)
        ZQXN(JL,JN)     = 0.0_JPRB
      ENDIF
    ENDDO
  ENDDO

  !--------------------------------
  ! variables needed for next level
  !--------------------------------
  DO JM=1,NCLV
    DO JL=KIDIA,KFDIA
      ZQXNM1(JL,JM)    = ZQXN(JL,JM)
      ZQXN2D(JL,JK,JM) = ZQXN(JL,JM)
    ENDDO
  ENDDO

  !------------------------------------------------------------------------
  ! 5.3 Precipitation/sedimentation fluxes to next level
  !     diagnostic precipitation fluxes
  !     It is this scaled flux that must be used for source to next layer
  !------------------------------------------------------------------------

  DO JM=1,NCLV
    DO JL=KIDIA,KFDIA
      ZPFPLSX(JL,JK+1,JM) = ZFALLSINK(JL,JM)*ZQXN(JL,JM)*ZRDTGDP(JL)
    ENDDO
  ENDDO

  ! Ensure precipitation fraction is zero if no precipitation
  DO JL=KIDIA,KFDIA
    ZQPRETOT(JL) =ZPFPLSX(JL,JK+1,NCLDQS)+ZPFPLSX(JL,JK+1,NCLDQR)
  ENDDO
  DO JL=KIDIA,KFDIA
    IF (ZQPRETOT(JL)<ZEPSEC) THEN
      ZCOVPTOT(JL)=0.0_JPRB
    ENDIF
  ENDDO
  
  ! Calculate diagnosed process rates for implicit quantities
  DO JL=KIDIA,KFDIA
    IF (IWARMRAIN == 3) THEN
      IF(ZTP1(JL,JK) <= RTT) THEN
        ZBUDL(JL,12) = -ZRAINAUT(JL)*ZQXN(JL,NCLDQL)*ZQTMST
        ZBUDL(JL,13) = -ZRAINACC(JL)*ZQXN(JL,NCLDQL)*ZQTMST
      ELSE
        ZBUDL(JL,14) = -ZRAINAUT(JL)*ZQXN(JL,NCLDQL)*ZQTMST
        ZBUDL(JL,15) = -ZRAINACC(JL)*ZQXN(JL,NCLDQL)*ZQTMST
      ENDIF
    ENDIF
    ZBUDL(JL,6) = -ZCONVSINK(JL,NCLDQL)*ZQXN(JL,NCLDQL)*ZQTMST
    ZBUDI(JL,6) = -ZCONVSINK(JL,NCLDQI)*ZQXN(JL,NCLDQI)*ZQTMST
    ZBUDI(JL,13)= -ZFALLSINK(JL,NCLDQI)*ZQXN(JL,NCLDQI)*ZQTMST
    ZBUDL(JL,21)= -ZSNOWRIME(JL)*ZQXN(JL,NCLDQL)*ZQTMST
    ZBUDI(JL,14)= -ZSNOWAUT(JL)*ZQXN(JL,NCLDQI)*ZQTMST
  ENDDO


  !######################################################################
  !
  !              6.  *** UPDATE TENDANCIES ***
  !
  !######################################################################

  !----------------------------------------------
  ! 6.1 Temperature and cloud condensate budgets 
  !----------------------------------------------

  DO JM=1,NCLV-1
    DO JL=KIDIA,KFDIA

      ! calculate fluxes in and out of box for conservation of TL
      ZFLUXQ(JL,JM)=ZPSUPSATSRCE(JL,JM)+ZCONVSRCE(JL,JM)+ZFALLSRCE(JL,JM)-&
                    & (ZFALLSINK(JL,JM)+ZCONVSINK(JL,JM))*ZQXN(JL,JM)
    ENDDO

    IF (IPHASE(JM)==1) THEN
      DO JL=KIDIA,KFDIA
        TENDENCY_LOC_T(JL,JK)=TENDENCY_LOC_T(JL,JK)+ &
          & RALVDCP*(ZQXN(JL,JM)-ZQX(JL,JK,JM)-ZFLUXQ(JL,JM))*ZQTMST
      ENDDO
    ENDIF

    IF (IPHASE(JM)==2) THEN
      DO JL=KIDIA,KFDIA
        TENDENCY_LOC_T(JL,JK)=TENDENCY_LOC_T(JL,JK)+ &
          & RALSDCP*(ZQXN(JL,JM)-ZQX(JL,JK,JM)-ZFLUXQ(JL,JM))*ZQTMST
      ENDDO
    ENDIF

      !----------------------------------------------------------------------
      ! New prognostic tendencies - ice,liquid rain,snow 
      ! Note: CLV arrays use PCLV in calculation of tendency while humidity
      !       uses ZQX. This is due to clipping at start of cloudsc which
      !       include the tendency already in TENDENCY_LOC_T and TENDENCY_LOC_q. ZQX was reset
      !----------------------------------------------------------------------
    DO JL=KIDIA,KFDIA
      TENDENCY_LOC_CLD(JL,JK,JM)=TENDENCY_LOC_CLD(JL,JK,JM)+(ZQXN(JL,JM)-ZQX0(JL,JK,JM))*ZQTMST
    ENDDO

  ENDDO
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA
    !----------------------
    ! 6.2 Humidity budget
    !----------------------
    TENDENCY_LOC_Q(JL,JK)=TENDENCY_LOC_Q(JL,JK)+(ZQXN(JL,NCLDQV)-ZQX(JL,JK,NCLDQV))*ZQTMST

    !-------------------
    ! 6.3 cloud cover 
    !-----------------------
    TENDENCY_LOC_A(JL,JK)=TENDENCY_LOC_A(JL,JK)+ZDA(JL)*ZQTMST

    !------------------------------------------------------------------
    ! check for supersaturation when Semi-Langrangian advection is off
    !------------------------------------------------------------------
    IF( .NOT. LDSLPHY ) THEN
      ZQNEW=PQ(JL,JK)+PTSPHY*(TENDENCY_LOC_Q(JL,JK)+TENDENCY_CML_Q(JL,JK))
      ZTNEW=PT(JL,JK)+PTSPHY*(TENDENCY_LOC_T(JL,JK)+TENDENCY_CML_T(JL,JK))
      ZANEW=PA(JL,JK)+PTSPHY*(TENDENCY_LOC_A(JL,JK)+TENDENCY_CML_A(JL,JK))
      ZALFAW=FOEALFA(ZTNEW)
      IF (ZTNEW>=RTT .OR. NSSOPT==0) THEN
        ZFAC=1.0_JPRB
        ZFACI=ZQTMST
      ELSE
        ZFAC=ZANEW+FOKOOP(ZTNEW)*(1.0_JPRB-ZANEW)
        ZFACI=1.0_JPRB/RKOOPTAU
      ENDIF
      ZQSAT=FOEEWM(ZTNEW)/PAP(JL,JK)
      ZQSAT=MIN(0.5_JPRB,ZQSAT)
      ZQSAT=ZQSAT/(1.0_JPRB-RETV*ZQSAT)
      IF (ZQNEW>ZFAC*ZQSAT) THEN 
        ZCOND=(ZQNEW-ZFAC*ZQSAT)/ZCORQSMIX(JL)
        TENDENCY_LOC_Q(JL,JK)=TENDENCY_LOC_Q(JL,JK)-ZCOND*ZQTMST
        TENDENCY_LOC_T(JL,JK)=TENDENCY_LOC_T(JL,JK)+FOELDCPM(ZTP1(JL,JK))*ZCOND*ZQTMST
        TENDENCY_LOC_A(JL,JK)=TENDENCY_LOC_A(JL,JK)+(1.0_JPRB-ZANEW)*ZFACI
        TENDENCY_LOC_CLD(JL,JK,NCLDQL)=TENDENCY_LOC_CLD(JL,JK,NCLDQL)+ZCOND*ZALFAW*ZQTMST
        TENDENCY_LOC_CLD(JL,JK,NCLDQI)=TENDENCY_LOC_CLD(JL,JK,NCLDQI)+ZCOND*(1.0_JPRB-ZALFAW)*ZQTMST
      ENDIF
    ENDIF
  ENDDO


!--------------------------------------------------
! Copy precipitation fraction into output variable
!-------------------------------------------------
 DO JL=KIDIA,KFDIA
    PCOVPTOT(JL,JK) = ZCOVPTOT(JL)
 ENDDO

  ! cloud fraction/precip fraction at the end of the routine
  DO JL=KIDIA,KFDIA
     IF (PA(JL,JK)+ZDA(JL) > 0.001_JPRB) THEN !same cloud fraction threshold as in radiation
     ZANEWP(JL,JK) = MAX(MAX(PA(JL,JK)+ZDA(JL),PCOVPTOT(JL,JK)),0.1_JPRB)
     ENDIF
  ENDDO
!######################################################################
!
!              7.  *** CLOUD BUDGET DIAGNOSTICS ***
!
!######################################################################

  IS = 0  ! Set extra diagnostic counter. Note these diagnostics currently 
          ! won't work with LBUD23 also turned on
          
  !-----------------------------------------------------------------
  ! Vertical integral of all cloud process terms in one 3D field 
  ! Requires certain number of levels.
  ! At some point need to move this to individual 2D diagnostic fields
  !-----------------------------------------------------------------
  IF (LCLDBUD_VERTINT) THEN

   IF (KLEV < 60) CALL ABOR1('CLOUDSC ERROR: Not enough levels for cloud vertical integral budget.')
   
   DO JL=KIDIA,KFDIA

      ! Layer depth (m)
      ZDZ = ZDP(JL)/(ZRHO(JL)*RG)
            
      PEXTRA(JL,1,1)  = PEXTRA(JL,1,1) + ZBUDCC(JL,1)*ZDZ  ! + Supersat clipping so far this timestep                          
      PEXTRA(JL,2,1)  = PEXTRA(JL,2,1) + ZBUDCC(JL,2)*ZDZ  ! + Supersat clipping from t-1 sltend (PSUPSAT)                     
      PEXTRA(JL,3,1)  = PEXTRA(JL,3,1) + ZBUDCC(JL,3)*ZDZ  ! + Convective detrainment                                          
      PEXTRA(JL,4,1)  = PEXTRA(JL,4,1) + (ZBUDCC(JL,4)+ZBUDCC(JL,6))*ZDZ ! +- Convective subsidence source and sink
      PEXTRA(JL,5,1)  = PEXTRA(JL,5,1) + ZBUDCC(JL,7)*ZDZ  ! - Turbulent erosion                                               
      PEXTRA(JL,6,1)  = PEXTRA(JL,6,1) + ZBUDCC(JL,10)*ZDZ ! + Condensation of new cloud      
      PEXTRA(JL,7,1)  = PEXTRA(JL,7,1) + ZBUDCC(JL,11)*ZDZ ! Tidy up       
      PEXTRA(JL,8,1)  = PEXTRA(JL,8,1) + PVFA(JL,JK)*ZDZ   ! Vertical diffusion
      PEXTRA(JL,9,1)  = PEXTRA(JL,9,1) + PDYNA(JL,JK)*ZDZ  ! Advection from dynamics    
      
      PEXTRA(JL,11,1) = PEXTRA(JL,11,1) + ZBUDL(JL,1)*ZDZ ! + Supersat clipping so far this timestep                          
      PEXTRA(JL,12,1) = PEXTRA(JL,12,1) + ZBUDL(JL,2)*ZDZ ! + Supersat clipping from t-1 sltend (PSUPSAT)                     
      PEXTRA(JL,13,1) = PEXTRA(JL,13,1) + ZBUDL(JL,3)*ZDZ ! + Convective detrainment                                          
      PEXTRA(JL,14,1) = PEXTRA(JL,14,1) + (ZBUDL(JL,4)+ZBUDL(JL,6))*ZDZ ! +- Convective subsidence source and sink
      PEXTRA(JL,15,1) = PEXTRA(JL,15,1) + ZBUDL(JL,5)*ZDZ ! +- Evaporation due to convective subsidence 
      PEXTRA(JL,16,1) = PEXTRA(JL,16,1) + ZBUDL(JL,7)*ZDZ ! - Turbulent erosion                                               
      PEXTRA(JL,17,1) = PEXTRA(JL,17,1) + ZBUDL(JL,8)*ZDZ ! - Evaporation of existing cloud (dqs increasing = subsat)         
      PEXTRA(JL,18,1) = PEXTRA(JL,18,1) + ZBUDL(JL,9)*ZDZ  ! + Condensation of existing cloud (dqs decreasing = supersat)      
      PEXTRA(JL,19,1) = PEXTRA(JL,19,1) + ZBUDL(JL,10)*ZDZ ! + Condensation of new cloud (dqs decreasing = supersat)           
      PEXTRA(JL,20,1) = PEXTRA(JL,20,1) + ZBUDL(JL,11)*ZDZ! - Deposition of liquid to ice                                     
      PEXTRA(JL,21,1) = PEXTRA(JL,21,1) + ZBUDL(JL,12)*ZDZ! - Autoconversion to rain+freezing->snow (IMPLICIT)(ZRAINAUT) 
      PEXTRA(JL,22,1) = PEXTRA(JL,22,1) + ZBUDL(JL,13)*ZDZ! - Accretion of cloud to rain+freezing->snow (IMPLICIT)(ZRAINACC)      
      PEXTRA(JL,23,1) = PEXTRA(JL,23,1) + ZBUDL(JL,14)*ZDZ! - Autoconversion to rain (IMPLICIT)(ZRAINAUT)                     
      PEXTRA(JL,24,1) = PEXTRA(JL,24,1) + ZBUDL(JL,15)*ZDZ! - Accretion of cloud to rain (IMPLICIT)(ZRAINACC)      
      PEXTRA(JL,25,1) = PEXTRA(JL,25,1) + ZBUDL(JL,17)*ZDZ! + Melting of ice to liquid                                        
      PEXTRA(JL,26,1) = PEXTRA(JL,26,1) + (ZBUDL(JL,18)+ZBUDL(JL,19))*ZDZ ! - Freezing of rain-to-snow, liq-to-ice
      PEXTRA(JL,27,1) = PEXTRA(JL,27,1) + ZBUDL(JL,20)*ZDZ ! - Evaporation of rain
      PEXTRA(JL,28,1) = PEXTRA(JL,28,1) + ZBUDL(JL,21)*ZDZ ! - Riming of cloud liquid to snow
      PEXTRA(JL,29,1) = PEXTRA(JL,29,1) + PVFL(JL,JK)*ZDZ  ! +- Vertical diffusion
      PEXTRA(JL,30,1) = PEXTRA(JL,30,1) + PDYNL(JL,JK)*ZDZ  ! +- Advection from dynamics

      PEXTRA(JL,41,1) = PEXTRA(JL,41,1) + ZBUDI(JL,1)*ZDZ  ! + Supersat clipping so far this timestep                          
      PEXTRA(JL,42,1) = PEXTRA(JL,42,1) + ZBUDI(JL,2)*ZDZ  ! + Supersat clipping from t-1 sltend (PSUPSAT)                     
      PEXTRA(JL,43,1) = PEXTRA(JL,43,1) + ZBUDI(JL,3)*ZDZ  ! + Convective detrainment                                          
      PEXTRA(JL,44,1) = PEXTRA(JL,44,1) + (ZBUDI(JL,4)+ZBUDI(JL,6))*ZDZ ! +- Convective subsidence source and sink
      PEXTRA(JL,45,1) = PEXTRA(JL,45,1) + ZBUDI(JL,5)*ZDZ  ! +- Evaporation due to convective subsidence 
      PEXTRA(JL,46,1) = PEXTRA(JL,46,1) + ZBUDI(JL,7)*ZDZ  ! - Turbulent erosion                                               
      PEXTRA(JL,47,1) = PEXTRA(JL,47,1) + ZBUDI(JL,8)*ZDZ  ! - Evaporation of existing cloud (dqs increasing = subsat)         
      PEXTRA(JL,48,1) = PEXTRA(JL,48,1) + ZBUDI(JL,9)*ZDZ  ! + Condensation of existing cloud (dqs decreasing = supersat)      
      PEXTRA(JL,49,1) = PEXTRA(JL,49,1) + ZBUDI(JL,10)*ZDZ ! + Condensation of new cloud (dqs decreasing = supersat)           
      PEXTRA(JL,50,1) = PEXTRA(JL,50,1) + ZBUDI(JL,11)*ZDZ ! + Deposition of liquid to ice                                     
      PEXTRA(JL,51,1) = PEXTRA(JL,51,1) + (ZBUDI(JL,12)+ZBUDI(JL,13))*ZDZ  ! Ice sedimentation
      PEXTRA(JL,52,1) = PEXTRA(JL,52,1) + ZBUDI(JL,14)*ZDZ ! - Autoconversion to snow (IMPLICIT)(ZSNOWAUT)
      PEXTRA(JL,53,1) = PEXTRA(JL,53,1) + (ZBUDI(JL,15)+ZBUDI(JL,16))*ZDZ ! - Melting of ice/snow to rain
      PEXTRA(JL,54,1) = PEXTRA(JL,54,1) + ZBUDI(JL,17)*ZDZ ! - Evaporation of rain
      PEXTRA(JL,55,1) = PEXTRA(JL,55,1) + PVFI(JL,JK)*ZDZ  ! +- Vertical diffusion
      PEXTRA(JL,56,1) = PEXTRA(JL,56,1) + PDYNI(JL,JK)*ZDZ ! +- Advection from dynamics
    
    ENDDO
    IS = IS + 1
    IF (KFLDX < IS) CALL ABOR1('CLOUDSC ERROR: Not enough PEXTRA variables for cloud vertical integral budget.')
  ENDIF

  !-----------------------------------------------------------------
  ! Cloud fraction budget 
  !-----------------------------------------------------------------
  IF (LCLDBUDC) THEN
    DO JL=KIDIA,KFDIA
      PEXTRA(JL,JK,IS+1)  = PEXTRA(JL,JK,IS+1) + PA(JL,JK)     ! Initial total cloud fraction
      PEXTRA(JL,JK,IS+2)  = PEXTRA(JL,JK,IS+2) + ZBUDCC(JL,1)  ! Supersat clipping
      PEXTRA(JL,JK,IS+3)  = PEXTRA(JL,JK,IS+3) + ZBUDCC(JL,2)  ! Supersat clipping from t-1 sltend
      PEXTRA(JL,JK,IS+4)  = PEXTRA(JL,JK,IS+4) + ZBUDCC(JL,3)  ! Convective detrainment
      PEXTRA(JL,JK,IS+5)  = PEXTRA(JL,JK,IS+5) + ZBUDCC(JL,4)+ZBUDCC(JL,6)  ! Convective mass flux
      PEXTRA(JL,JK,IS+6)  = PEXTRA(JL,JK,IS+6) + ZBUDCC(JL,7)  ! Turbulent erosion
      PEXTRA(JL,JK,IS+7)  = PEXTRA(JL,JK,IS+7) + ZBUDCC(JL,10) ! Condensation of new cloud
      PEXTRA(JL,JK,IS+8)  = PEXTRA(JL,JK,IS+8) + ZBUDCC(JL,11) ! Small value tidy up
      PEXTRA(JL,JK,IS+9)  = PEXTRA(JL,JK,IS+9) + PVFA(JL,JK)   ! Vertical diffusion
      PEXTRA(JL,JK,IS+10) = PEXTRA(JL,JK,IS+10) + PDYNA(JL,JK)  ! Advection from dynamics
      PEXTRA(JL,JK,IS+11) = PEXTRA(JL,JK,IS+11) + PA(JL,JK)+ZDA(JL) ! Final cloud frac
    ENDDO
    IS = IS + 11
    IF (KFLDX < IS) CALL ABOR1('CLOUDSC ERROR: Not enough PEXTRA variables for cloud fraction budget.')
  ENDIF

  !-----------------------------------------------------------------
  ! Cloud liquid condensate budget 
  !-----------------------------------------------------------------
  IF (LCLDBUDL) THEN
    DO JL=KIDIA,KFDIA
      PEXTRA(JL,JK,IS+1)  = PEXTRA(JL,JK,IS+1) + ZQX0(JL,JK,NCLDQL) ! Initial condensate
      PEXTRA(JL,JK,IS+2)  = PEXTRA(JL,JK,IS+2) + ZBUDL(JL,1) ! + Supersat clipping so far this timestep
      PEXTRA(JL,JK,IS+3)  = PEXTRA(JL,JK,IS+3) + ZBUDL(JL,2) ! + Supersat clipping from t-1 sltend (PSUPSAT)
      PEXTRA(JL,JK,IS+4)  = PEXTRA(JL,JK,IS+4) + ZBUDL(JL,3) ! + Convective detrainment
      PEXTRA(JL,JK,IS+5)  = PEXTRA(JL,JK,IS+5) + ZBUDL(JL,4)+ZBUDL(JL,5)+ZBUDL(JL,6) ! +- Convective subsidence
       ! ZBUDL(JL,4) + Convective subsidence source from layer above
       ! ZBUDL(JL,5) - Convective subsidence source evaporation in layer
       ! ZBUDL(JL,6) - Convective subsidence sink to layer below (IMPLICIT) (ZCONVSINK)
      PEXTRA(JL,JK,IS+6)  = PEXTRA(JL,JK,IS+6) + ZBUDL(JL,7) ! - Turbulent erosion
      PEXTRA(JL,JK,IS+7)  = PEXTRA(JL,JK,IS+7) + ZBUDL(JL,8) ! - Evaporation of existing cloud (dqs increasing = subsat)
      PEXTRA(JL,JK,IS+8)  = PEXTRA(JL,JK,IS+8) + ZBUDL(JL,9)  ! + Condensation of existing cloud (dqs decreasing = supersat)
      PEXTRA(JL,JK,IS+9)  = PEXTRA(JL,JK,IS+9) + ZBUDL(JL,10) ! + Condensation of new cloud (dqs decreasing = supersat)
      PEXTRA(JL,JK,IS+10) = PEXTRA(JL,JK,IS+10) + ZBUDL(JL,11)! - Deposition of liquid to ice
      PEXTRA(JL,JK,IS+11) = PEXTRA(JL,JK,IS+11) + ZBUDL(JL,12)! - Autoconversion to rain+freezing->snow (IMPLICIT)(ZRAINAUT) 
      PEXTRA(JL,JK,IS+12) = PEXTRA(JL,JK,IS+12) + ZBUDL(JL,13)! - Accretion of cloud to rain+freezing->snow (IMPLICIT)(ZRAINACC)
      PEXTRA(JL,JK,IS+13) = PEXTRA(JL,JK,IS+13) + ZBUDL(JL,14)! - Autoconversion to rain (IMPLICIT)(ZRAINAUT)
      PEXTRA(JL,JK,IS+14) = PEXTRA(JL,JK,IS+14) + ZBUDL(JL,15)! - Accretion of cloud to rain (IMPLICIT)(ZRAINACC)
      PEXTRA(JL,JK,IS+15) = PEXTRA(JL,JK,IS+15) + ZBUDL(JL,17)! + Melting of ice to liquid
      PEXTRA(JL,JK,IS+16) = PEXTRA(JL,JK,IS+16) + ZBUDL(JL,18)+ZBUDL(JL,19) ! - Freezing of rain/liq
       ! ZBUDL(JL,18) - Freezing of rain to snow
       ! ZBUDL(JL,19) - Freezing of liquid to ice
      PEXTRA(JL,JK,IS+17) = PEXTRA(JL,JK,IS+17) + ZBUDL(JL,20) ! - Evaporation of rain
      PEXTRA(JL,JK,IS+18) = PEXTRA(JL,JK,IS+18) + ZBUDL(JL,21) ! - Riming of cloud liquid to snow
      PEXTRA(JL,JK,IS+19) = PEXTRA(JL,JK,IS+19) + PVFL(JL,JK)  ! +- Vertical diffusion
      PEXTRA(JL,JK,IS+20) = PEXTRA(JL,JK,IS+20) + PDYNL(JL,JK) ! +- Advection from dynamics
      PEXTRA(JL,JK,IS+21) = PEXTRA(JL,JK,IS+21) + ZQXN(JL,NCLDQL)  ! Final condensate
    ENDDO
    IS = IS + 21
    IF (KFLDX < IS) CALL ABOR1('CLOUDSC ERROR: Not enough PEXTRA variables for cloud liquid budget.')
  ENDIF

  !-----------------------------------------------------------------
  ! Cloud ice condensate budget 
  !-----------------------------------------------------------------
  IF (LCLDBUDI) THEN
    DO JL=KIDIA,KFDIA
      PEXTRA(JL,JK,IS+1)  = PEXTRA(JL,JK,IS+1) + ZQX0(JL,JK,NCLDQI) ! Initial condensate
      PEXTRA(JL,JK,IS+2)  = PEXTRA(JL,JK,IS+2) + ZBUDI(JL,1)  ! + Supersat clipping so far this timestep 
      PEXTRA(JL,JK,IS+3)  = PEXTRA(JL,JK,IS+3) + ZBUDI(JL,2)  ! + Supersat clipping from t-1 sltend (PSUPSAT)
      PEXTRA(JL,JK,IS+4)  = PEXTRA(JL,JK,IS+4) + ZBUDI(JL,3)  ! + Convective detrainment
      PEXTRA(JL,JK,IS+5)  = PEXTRA(JL,JK,IS+5) + ZBUDI(JL,4)+ZBUDI(JL,5)+ZBUDI(JL,6)! +- Convective subsidence
       ! ZBUDI(JL,4) + Convective subsidence source from layer above
       ! ZBUDI(JL,5) - Convective subsidence source evaporation in layer
       ! ZBUDI(JL,6) - Convective subsidence sink to layer below (IMPLICIT) (ZCONVSINK)
      PEXTRA(JL,JK,IS+6)  = PEXTRA(JL,JK,IS+6) + ZBUDI(JL,7)  ! - Turbulent erosion
      PEXTRA(JL,JK,IS+7)  = PEXTRA(JL,JK,IS+7) + ZBUDI(JL,8)  ! - Evaporation of existing cloud (dqs increasing = subsat)
      PEXTRA(JL,JK,IS+8)  = PEXTRA(JL,JK,IS+8) + ZBUDI(JL,9)  ! + Condensation of existing cloud (dqs decreasing = supersat)
      PEXTRA(JL,JK,IS+9)  = PEXTRA(JL,JK,IS+9) + ZBUDI(JL,10) ! + Condensation of new cloud (dqs decreasing = supersat)
      PEXTRA(JL,JK,IS+10) = PEXTRA(JL,JK,IS+10) + ZBUDI(JL,11)! + Deposition of liquid to ice
      PEXTRA(JL,JK,IS+11) = PEXTRA(JL,JK,IS+11) + ZBUDI(JL,12)+ZBUDI(JL,13)! + Ice sedimentation
       ! ZBUDI(JL,12) + Ice sedimentation source from above
       ! ZBUDI(JL,13) - Ice sedimentation sink to below (IMPLICIT)(ZFALLSINK)
      PEXTRA(JL,JK,IS+12) = PEXTRA(JL,JK,IS+12) + ZBUDI(JL,14) ! - Autoconversion to snow (IMPLICIT) (ZSNOWAUT)
      PEXTRA(JL,JK,IS+13) = PEXTRA(JL,JK,IS+13) + ZBUDI(JL,15)+ZBUDI(JL,16) ! - Melting of ice/snow to rain
       ! ZBUDI(JL,15)! - Melting of ice to rain
       ! ZBUDI(JL,16)! - Melting of snow to rain
      PEXTRA(JL,JK,IS+14) = PEXTRA(JL,JK,IS+14) + ZBUDI(JL,17)! - Evaporation of snow
      PEXTRA(JL,JK,IS+15) = PEXTRA(JL,JK,IS+15) + PVFI(JL,JK) ! Vertical diffusion
      PEXTRA(JL,JK,IS+16) = PEXTRA(JL,JK,IS+16) + PDYNI(JL,JK)! Advection from dynamics
      PEXTRA(JL,JK,IS+17) = PEXTRA(JL,JK,IS+17) + ZQXN(JL,NCLDQI)  ! Final condensate
    ENDDO
    IS = IS + 17
    IF (KFLDX < IS) CALL ABOR1('CLOUDSC ERROR: Not enough PEXTRA variables for cloud ice budget.')
  ENDIF

  !-----------------------------------------------------------------
  ! Cloud processes temperature budget 
  !-----------------------------------------------------------------
  ! RALVDCP latent heat of condensation (vapour to liquid)
  ! RALSDCP latent heat of sublimation (vapour to solid)
  ! RALFDCP latent heat of melting/freezing (liquid to solid)
  !-----------------------------------------------------------------
  IF (LCLDBUDT) THEN
    DO JL=KIDIA,KFDIA
      ! Note, PEXTRA(:,:,1) is set to convective T tendency in callpar, so start at 2 here

      ! Radiative heating rates (shortwave and longwave)
      PEXTRA(JL,JK,IS+1) = PEXTRA(JL,JK,IS+1) + PHRSW(JL,JK)
      PEXTRA(JL,JK,IS+2) = PEXTRA(JL,JK,IS+2) + PHRLW(JL,JK)

      ! Condensation (vapour to liquid) (heating +ve RALVDCP) 
      !    = "Supersat adjust cloudsc" + "Supersat adjust sltend t-1" 
      !    + "Existing cloud" + "New cloud"

      PEXTRA(JL,JK,IS+3)  = PEXTRA(JL,JK,IS+3) + (ZBUDL(JL,1) + ZBUDL(JL,2) &
     &                     + ZBUDL(JL,9) + ZBUDL(JL,10))*RALVDCP&
     &                     -ZBUDI(JL,11)*RALVDCP

      ! Condensation (vapour to ice) (heating +ve RALSDCP) 
      !    = "Supersat adjust cloudsc" + "Supersat adjust sltend t-1" 
      !    + "Existing cloud" + "New cloud"
      PEXTRA(JL,JK,IS+4)  = PEXTRA(JL,JK,IS+4) + (ZBUDI(JL,1) + ZBUDI(JL,2) &
     &                     + ZBUDI(JL,9) + ZBUDI(JL,10)+ZBUDI(JL,11))*RALSDCP

      ! Deposition of liquid to ice  (Bergeron Findeisen liquid-vapour-ice) (heating +ve RALFDCP)        
      PEXTRA(JL,JK,IS+5) = PEXTRA(JL,JK,IS+5) + ZBUDI(JL,11) * RALFDCP                                      
      	
      ! Evaporation (liquid to vapour) (cooling -ve RALVDCP)
      !    = "Turbulent erosion" + "Evaporation of existing cloud" 
      !    + "Convective subsidence evap" 
      PEXTRA(JL,JK,IS+6)  = PEXTRA(JL,JK,IS+6) + (ZBUDL(JL,7) + ZBUDL(JL,8) + ZBUDL(JL,5))*RALVDCP

      ! Evaporation (ice to vapour) (cooling -ve RALSDCP)
      !    = "Turbulent erosion" + "Evaporation of existing cloud" 
      !    + "Convective subsidence evap" 
      PEXTRA(JL,JK,IS+7)  = PEXTRA(JL,JK,IS+7) + (ZBUDI(JL,7) + ZBUDI(JL,8) + ZBUDI(JL,5))*RALSDCP

      ! Evaporation of rain to vapour (cooling -ve RALVDCP)
      PEXTRA(JL,JK,IS+8) = PEXTRA(JL,JK,IS+8) + ZBUDL(JL,20) * RALVDCP

      ! Evaporation of snow to vapour (cooling -ve RALSDCP)
      PEXTRA(JL,JK,IS+9) = PEXTRA(JL,JK,IS+9) + ZBUDI(JL,17) * RALSDCP 

      ! Melting of ice to liquid (cooling -ve RALFDCP)
      PEXTRA(JL,JK,IS+10) = PEXTRA(JL,JK,IS+10) + ZBUDI(JL,15) * RALFDCP                                       

      ! Melting of snow to rain (cooling -ve RALFDCP)
      PEXTRA(JL,JK,IS+11) = PEXTRA(JL,JK,IS+11) + ZBUDI(JL,16) * RALFDCP                                       
             
      ! Freezing rate of liquid to ice (heating +ve RALFDCP)
      !    "Autoconversion & Accretion to rain+freezing->snow" 
      !    + "Freezing of rain to snow" + "Freezing of liquid to ice"
      PEXTRA(JL,JK,IS+12) = PEXTRA(JL,JK,IS+12) + (ZBUDL(JL,12) + ZBUDL(JL,13) + ZBUDL(JL,18) + &
     &                       ZBUDL(JL,19)) * RALFDCP 

      ! Riming of liquid to snow
      PEXTRA(JL,JK,IS+13) = PEXTRA(JL,JK,IS+13) + ZBUDL(JL,21) * RALFDCP                                       
      
    ENDDO
    IS = IS + 13
    IF (KFLDX < IS) CALL ABOR1('CLOUDSC ERROR: Not enough PEXTRA variables for cloud temperature budget.')
  ENDIF

  IF (YDEPHY%LRAD_CLOUD_INHOMOG) THEN
     DO JL=KIDIA,KFDIA
        ! for ice FSD
        !liquid, ice and detrained condensate at the beginning of timestep
        ! calculate ratio of detrained condensate to all condensate
        IF (ZQX0(JL,JK,NCLDQL)+ZQX0(JL,JK,NCLDQI)+PLUDE(JL,JK) > RLMIN) THEN
           ZRATFSD(JL,JK)=PLUDE(JL,JK)/(PLUDE(JL,JK)+ZQX0(JL,JK,NCLDQL)+ZQX0(JL,JK,NCLDQI))
           !safety - ratio between 0 and 1
           ZRATFSD(JL,JK)=MAX(0.0_JPRB,MIN(1.0_JPRB,ZRATFSD(JL,JK)))
        ENDIF
        ! final liquid and ice condensate
        ZICEF(JL,JK)=ZQXN(JL,NCLDQI)
        ZLIQF(JL,JK)=ZQXN(JL,NCLDQL)
     ENDDO
  ENDIF

ENDDO ! on vertical level JK
!----------------------------------------------------------------------
!
!                  END OF VERTICAL LOOP OVER LEVELS
!
!----------------------------------------------------------------------



!######################################################################
!
!              8.  *** FLUX/DIAGNOSTICS COMPUTATIONS ***
!
!######################################################################

!-------------------------------------
! Enthalpy and total water diagnostics (LSCMEC and LCLDBUDGET normally false)
!-------------------------------------
IF (LSCMEC.OR.LCLDBUDGET) THEN

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZTNEW=PT(JL,JK)+PTSPHY*(TENDENCY_LOC_T(JL,JK)+TENDENCY_CML_T(JL,JK))
      IF (JK==1) THEN
        ZSUMQ1(JL,JK)=0.0_JPRB
        ZSUMH1(JL,JK)=0.0_JPRB
      ELSE
        ZSUMQ1(JL,JK)=ZSUMQ1(JL,JK-1)
        ZSUMH1(JL,JK)=ZSUMH1(JL,JK-1)
      ENDIF

      ! cld vars
      DO JM=1,NCLV-1
        IF (IPHASE(JM)==1) ZTNEW=ZTNEW-RALVDCP*(PCLV(JL,JK,JM)+ &
          & (TENDENCY_LOC_CLD(JL,JK,JM)+TENDENCY_CML_CLD(JL,JK,JM))*PTSPHY)
        IF (IPHASE(JM)==2) ZTNEW=ZTNEW-RALSDCP*(PCLV(JL,JK,JM)+ &
          & (TENDENCY_LOC_CLD(JL,JK,JM)+TENDENCY_CML_CLD(JL,JK,JM))*PTSPHY)
        ZSUMQ1(JL,JK)=ZSUMQ1(JL,JK)+ &
        & (PCLV(JL,JK,JM)+(TENDENCY_LOC_CLD(JL,JK,JM)+TENDENCY_CML_CLD(JL,JK,JM))*PTSPHY)* &
        & (PAPH(JL,JK+1)-PAPH(JL,JK))*ZRG_R
      ENDDO
      ZSUMH1(JL,JK)=ZSUMH1(JL,JK)+(PAPH(JL,JK+1)-PAPH(JL,JK))*ZTNEW 

      ! humidity
      ZSUMQ1(JL,JK)=ZSUMQ1(JL,JK)+ &
        &(PQ(JL,JK)+(TENDENCY_LOC_Q(JL,JK)+TENDENCY_CML_Q(JL,JK))*PTSPHY)*(PAPH(JL,JK+1)-PAPH(JL,JK))*ZRG_R

      ZRAIN=0.0_JPRB
      DO JM=1,NCLV
        ZRAIN=ZRAIN+PTSPHY*ZPFPLSX(JL,JK+1,JM)
      ENDDO
      ZERRORQ(JL,JK)=ZSUMQ1(JL,JK)+ZRAIN-ZSUMQ0(JL,JK)
    ENDDO
  ENDDO

  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
      ZDTGDP(JL)=PTSPHY*RG/(PAPH(JL,JK+1)-PAPH(JL,JK))
      ZRAIN=0.0_JPRB
      DO JM=1,NCLV
        IF (IPHASE(JM)==1) ZRAIN=ZRAIN+RALVDCP*ZDTGDP(JL)*ZPFPLSX(JL,JK+1,JM)* &
                           & (PAPH(JL,JK+1)-PAPH(JL,JK))
        IF (IPHASE(JM)==2) ZRAIN=ZRAIN+RALSDCP*ZDTGDP(JL)*ZPFPLSX(JL,JK+1,JM)* &
                           & (PAPH(JL,JK+1)-PAPH(JL,JK))
      ENDDO
      ZSUMH1(JL,JK)=(ZSUMH1(JL,JK)-ZRAIN)/PAPH(JL,JK+1)
      ZERRORH(JL,JK)=ZSUMH1(JL,JK)-ZSUMH0(JL,JK)
    ENDDO
  ENDDO

  DO JL=KIDIA,KFDIA
    IF (ABS(ZERRORQ(JL,KLEV))>1.E-13_JPRB.OR.ABS(ZERRORH(JL,KLEV))>1.E-13_JPRB) THEN
      ZQADJ=0.0_JPRB ! dummy statement
                     ! place totalview break here to catch non-conservation
    ENDIF
  ENDDO

  IF (ALLOCATED(ZSUMQ0))  DEALLOCATE(ZSUMQ0)
  IF (ALLOCATED(ZSUMQ1))  DEALLOCATE(ZSUMQ1)
  IF (ALLOCATED(ZSUMH0))  DEALLOCATE(ZSUMH0)
  IF (ALLOCATED(ZSUMH1))  DEALLOCATE(ZSUMH1)
  IF (ALLOCATED(ZERRORQ)) DEALLOCATE(ZERRORQ)
  IF (ALLOCATED(ZERRORH)) DEALLOCATE(ZERRORH)

ENDIF

!--------------------------------------------------------------------
! Copy general precip arrays back into PFP arrays for GRIB archiving
! Add rain and liquid fluxes, ice and snow fluxes
!--------------------------------------------------------------------
DO JK=1,KLEV+1
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA
    PFPLSL(JL,JK) = ZPFPLSX(JL,JK,NCLDQR)+ZPFPLSX(JL,JK,NCLDQL)
    PFPLSN(JL,JK) = ZPFPLSX(JL,JK,NCLDQS)+ZPFPLSX(JL,JK,NCLDQI)
  ENDDO
ENDDO

!--------
! Fluxes:
!--------
!DIR$ IVDEP
DO JL=KIDIA,KFDIA
  PFSQLF(JL,1)  = 0.0_JPRB
  PFSQIF(JL,1)  = 0.0_JPRB
  PFSQRF(JL,1)  = 0.0_JPRB
  PFSQSF(JL,1)  = 0.0_JPRB
  PFCQLNG(JL,1) = 0.0_JPRB
  PFCQNNG(JL,1) = 0.0_JPRB
  PFCQRNG(JL,1) = 0.0_JPRB !rain
  PFCQSNG(JL,1) = 0.0_JPRB !snow
! fluxes due to turbulence
  PFSQLTUR(JL,1) = 0.0_JPRB
  PFSQITUR(JL,1) = 0.0_JPRB
ENDDO

DO JK=1,KLEV
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA

    ZGDPH_R = -ZRG_R*(PAPH(JL,JK+1)-PAPH(JL,JK))*ZQTMST
    PFSQLF(JL,JK+1)  = PFSQLF(JL,JK)
    PFSQIF(JL,JK+1)  = PFSQIF(JL,JK)
    PFSQRF(JL,JK+1)  = PFSQLF(JL,JK)
    PFSQSF(JL,JK+1)  = PFSQIF(JL,JK)
    PFCQLNG(JL,JK+1) = PFCQLNG(JL,JK)
    PFCQNNG(JL,JK+1) = PFCQNNG(JL,JK)
    PFCQRNG(JL,JK+1) = PFCQLNG(JL,JK)
    PFCQSNG(JL,JK+1) = PFCQNNG(JL,JK)
    PFSQLTUR(JL,JK+1) = PFSQLTUR(JL,JK)
    PFSQITUR(JL,JK+1) = PFSQITUR(JL,JK)

    ZALFAW=ZFOEALFA(JL,JK)

    ! Liquid , LS scheme minus detrainment
    PFSQLF(JL,JK+1)=PFSQLF(JL,JK+1)+ &
   ! &(ZQXN2D(JL,JK,NCLDQL)-ZQX0(JL,JK,NCLDQL)+PVFL(JL,JK)*PTSPHY-PLUDELI(JL,JK,1))*ZGDPH_R
     &(ZQXN2D(JL,JK,NCLDQL)-ZQX0(JL,JK,NCLDQL)+PVFL(JL,JK)*PTSPHY)*ZGDPH_R+PLUDELI(JL,JK,1)
    ! liquid, negative numbers 
    PFCQLNG(JL,JK+1)=PFCQLNG(JL,JK+1)+ZLNEG(JL,JK,NCLDQL)*ZGDPH_R

    ! liquid, vertical diffusion
    PFSQLTUR(JL,JK+1)=PFSQLTUR(JL,JK+1)+PVFL(JL,JK)*PTSPHY*ZGDPH_R

    ! Rain, LS scheme 
    PFSQRF(JL,JK+1)=PFSQRF(JL,JK+1)+(ZQXN2D(JL,JK,NCLDQR)-ZQX0(JL,JK,NCLDQR))*ZGDPH_R 
    ! rain, negative numbers
    PFCQRNG(JL,JK+1)=PFCQRNG(JL,JK+1)+ZLNEG(JL,JK,NCLDQR)*ZGDPH_R

    ! Ice , LS scheme minus detrainment
    PFSQIF(JL,JK+1)=PFSQIF(JL,JK+1)+ &
   ! & (ZQXN2D(JL,JK,NCLDQI)-ZQX0(JL,JK,NCLDQI)+PVFI(JL,JK)*PTSPHY-PLUDELI(JL,JK,2))*ZGDPH_R
     & (ZQXN2D(JL,JK,NCLDQI)-ZQX0(JL,JK,NCLDQI)+PVFI(JL,JK)*PTSPHY)*ZGDPH_R+PLUDELI(JL,JK,2)
     ! ice, negative numbers
    PFCQNNG(JL,JK+1)=PFCQNNG(JL,JK+1)+ZLNEG(JL,JK,NCLDQI)*ZGDPH_R

    ! ice, vertical diffusion
    PFSQITUR(JL,JK+1)=PFSQITUR(JL,JK+1)+PVFI(JL,JK)*PTSPHY*ZGDPH_R

    ! snow, LS scheme
    PFSQSF(JL,JK+1)=PFSQSF(JL,JK+1)+(ZQXN2D(JL,JK,NCLDQS)-ZQX0(JL,JK,NCLDQS))*ZGDPH_R 
    ! snow, negative numbers
    PFCQSNG(JL,JK+1)=PFCQSNG(JL,JK+1)+ZLNEG(JL,JK,NCLDQS)*ZGDPH_R
  ENDDO
ENDDO

!-----------------------------------
! enthalpy flux due to precipitation
!-----------------------------------
DO JK=1,KLEV+1
!DIR$ IVDEP
  DO JL=KIDIA,KFDIA
    PFHPSL(JL,JK) = -RLVTT*PFPLSL(JL,JK)
    PFHPSN(JL,JK) = -RLSTT*PFPLSN(JL,JK)
  ENDDO
ENDDO

! Ice FSD calculation
!-----------------------------------
IF (YDEPHY%LRAD_CLOUD_INHOMOG) THEN
  DO JK=1,KLEV
    DO JL=KIDIA,KFDIA
       ! Representative grid box length (km)
       ! PGAW = normalised gaussian quadrature weight / no. longitude pts
       ZGRIDLEN = 2*RA*SQRT(RPI*PGAW(JL))*0.001_JPRB
       ZIFSDBACK=ZU1*ZGRIDLEN**(1._JPRB/3._JPRB)*((ZU2*ZGRIDLEN)**ZU3+1._JPRB)**ZU4
       
     !param updated in Sept 2017.
     ZQTOT = max(min((ZQX(JL,JK,NCLDQV)+ZQX(JL,JK,NCLDQI))*1000._JPRB,30._JPRB),0.00001)
     ZDP(JL)     = PAPH(JL,JK+1)-PAPH(JL,JK)     ! dp
     ZRHO(JL)    = PAP(JL,JK)/(RD*ZTP1(JL,JK))   ! p/RT air density
     ZDELZ =1._JPRB/(RG*ZRHO(JL))*ZDP(JL)*0.001_JPRB  !layer thickness in km

     ZDELZ=max(min(ZDELZ,.5),.001)
     ZPHI=(ZDELZ/.24)**.11*(ZQTOT/10.)**.03

     IF (ZANEWP(JL,JK) > 0.95_JPRB) THEN !treat as overcast
        ZIFSD(JL,JK)=ZIFSDBACK*ZPHI
        !add detrainment enhancement, multiply by 1D-to-2D enhancement factor
        ZIFSD(JL,JK)=ZR12*(ZIFSD(JL,JK)+ZRATFSD(JL,JK)*1.5_JPRB)
     ELSEIF (ZANEWP(JL,JK) > 0.001_JPRB .AND. ZANEWP(JL,JK) <= 0.95_JPRB) THEN 
        ZIFSD(JL,JK)=ZIFSDBACK*ZPHI*1.5_JPRB
        !add detrainment enhancement, multiply by 1D-to-2D enhancement factor
        ZIFSD(JL,JK)=ZR12*(ZIFSD(JL,JK)+ZRATFSD(JL,JK)*1.5_JPRB)
     ENDIF 

     !assign liquid or ice fsd based on liquid fraction. 
     !Global default value set to 1. all other cases
     ZFSD(JL,JK)=YDERAD%RCLOUD_FRAC_STD
     IF (ZICEF(JL,JK) > ZEPSEC ) THEN
        ZFSD(JL,JK)=ZIFSD(JL,JK)
     ENDIF
     IF (ZLIQF(JL,JK) > 0._JPRB .AND. ZICEF(JL,JK) <= ZEPSEC) THEN
        ZFSD(JL,JK)=ZLFSD(JL,JK) 
     ENDIF

     ! consistent with limits of possible FSD [0.1,3.575] in 
     ! lookup tables for Gamma/Log-normal functions
     ZFSD(JL,JK)=MIN(MAX(0.1_JPRB,ZFSD(JL,JK)),3.575_JPRB)

    ENDDO
  ENDDO
ENDIF ! LRAD_CLOUD_INHOMOG

! assign FSD to PFSD variable
IF (YDEPHY%LRAD_CLOUD_INHOMOG) THEN
   DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
         PFSD(JL,JK)=ZFSD(JL,JK)
      ENDDO
   ENDDO
ENDIF

!===============================================================================
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CLOUDSC',1,ZHOOK_HANDLE)

END SUBROUTINE CLOUDSC
