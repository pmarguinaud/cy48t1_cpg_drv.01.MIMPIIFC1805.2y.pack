SUBROUTINE VDFEXCUHL(YDVDF,YDEPHY,KIDIA  , KFDIA  , KLON   , KLEV   , PTMST  , &
                  &PUM1   , PVM1   , PTM1   , PQM1   , PLM1   , PIM1   , &
                  &PSLGM1 , PQTM1  ,  &
                  &PKMFL  , PKHFL  , PKQFL  ,  &
                  &PAPHM1 , PAPM1  , PGEOM1 , PGEOH, &
                  &PZINV  ,  KPBLTYPE , KDRAFT , &
                  &PQSVAR , PDQSDTEMP , &
                  &ZBUOY_COR, ZEN_CASC , ZWU ,&
                  &PGP2DSPP, KGFL_EZDIAG, PEZDIAG, &
!
! OUTPUT
                  &PTKE   , PMFLX , &
                  &ZLENGTH_M, ZLENGTH_H )

!
!     ------------------------------------------------------------------

!**   *VDFEXCUHL* - DETERMINES THE EXCHANGE COEFFICIENTS BETWEEN THE
!                 UPPER MODEL LEVELS WITH STABILITY AS A FUNCTION OF
!                 OBUKHOV-L


!     Original  A.C.M. BELJAARS       E.C.M.W.F.    26/03/90.
!     Modified  A.C.M. BELJAARS   26/03/99   Tiling of the land surface.
!     Modified  Geert Lenderink 2005   TKE scheme
!     Modified  Stephan de Roode 1 September 2005  Extension to moist convection
!     Modified  Wim de Rooy Implementation LHARATU in Harmonie June 2015
!     Modified  Ulf Andrae, Introduce SPP in HARMONIE-AROME

!     PURPOSE
!     -------

!     DETERMINE EXCHANGE COEFFICIENTS BETWEEN THE UPPER MODEL LEVELS

!     INTERFACE
!     ---------

!     *VDFEXCU* IS CALLED BY *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLON*         NUMBER OF GRID POINTS PER PACKET
!     *KLEV*         NUMBER OF LEVELS
!     *KVARTOP*      LEVEL OF MIXED LAYER TOP OR CUMULUS CLOUD TOP
!     *KHPBL*        PBL HEIGHT INDEX
!     *KPBLTYPE*     0=SBL,1=CBL,2=SC,3=SHALLOW CU,4=DEEP CU
!     *KDRAFT*       NUMBER OF EXPLICITLY MODELED DRAFTS - CURRENTLY SET TO 3 IN VDFMAIN
!     *KFLDX2*       NUMBER OF VARIABLES IN EXTRA DIAGNOSTIC FIELDS
!     *KLEVX*        NUMBER OF LEVELS IN EXTRA MULTI-LEVEL DIAGNOSTIC FIELDS
!     *KFLDX*        NUMBER OF VARIABLES IN EXTRA DIAGNOSTIC FIELDS


!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)
!     *PUM1*         X-VELOCITY COMPONENT          AT T-1
!     *PVM1*         Y-VELOCITY COMPONENT          AT T-1
!     *PTM1*         TEMPERATURE                   AT T-1
!     *PQM1*         SPECIFIC HUMIDITY             AT T-1
!     *PLM1*         SPECIFIC CLOUD LIQUID WATER   AT T-1                KG/KG
!     *PIM1*         SPECIFIC CLOUD ICE            AT T-1                KG/KG
!     *PSLGM1*       LIQUID WATER STATIC ENERGY    AT T-1
!     *PQTM1*        TOTAL SPECIFIC HUMIDITY       AT T-1 (INCLUDING ICE AND LIQUID)
!     *PAM1*         CLOUD FRACTION                AT T-1  
!     *PKMFL*        KINEMATIC MOMENTUM FLUX                [#]
!     *PKHFL*        KINEMATIC HEAT FLUX                    [#]
!     *PKQFL*        KINEMATIC MOISTURE FLUX                [#]
!     *PMCU*         CLOUD-DEPTH AVERAGE MOIST UPDRAFT MASS FLUX 
!     *PAPHM1*       HALF LEVEL PRESSURE           AT T-1
!     *PAPM1*        FULL LEVEL PRESSURE           AT T-1
!     *PRICUI*       RECIPROCAL CUMULUS INVERSION RICHARDSON NUMBER 
!     *PGEOM1*       GEOPOTENTIAL ON FULL LEVELS   AT T-1
!     *PGEOH *       GEOPOTENTIAL ON HALF LEVELS   AT T-1
!     *PZINV*        PBL HEIGHT (MOIST PARCEL, NOT FOR STABLE PBL) M 
!     *PQSVAR*       SATURATION SPECIFIC HUMIDITY (FROM TEMPERATURE AFTER DYNAMICS AND RADIATION)
!     *PDQSDTEMP*    QSAT DERIVATIVE WRT TEMPERATURE
!     *ZBUOY_COR*    STABILITY CORRECTION PARAMETER COMPUTED IN VDFHGHTN
!     *ZEN_CASC*     MASS FLUX TO TKE ENERGY CASCADE
!     *ZWU*          VERTICAL VELOCITY OF SECOND UPDRAFT



!     OUTPUT PARAMETERS (REAL):

!     *PCFM*         PROP. TO EXCH. COEFF. FOR MOMENTUM (C-STAR IN DOC.)
!     *PCFH*         PROP. TO EXCH. COEFF. FOR HEAT     (C-STAR IN DOC.)
!                    (ONLY PCFM(*,1:KLEV-1) AND
!                          PCFH(*,1:KLEV-1) ARE COMPUTED)
!     *PKH*          TURB. DIFF. COEFF. FOR HEAT ABOVE SURF. LAY.  (M2/S)
!     *PTKE*         TURBULENT KINETIC ENERGY AT T-1, UPDATED IN THIS ROUTINE  M^2/S^2
!     *PMFLX*        CLOUD MASS FLUX
!     *PEXTR2*       EXTRA DIAGNOSTIC FIELDS


!     ------
!     METHOD
!     ------

!     SEE  LENDERINK AND HOLTSLAG, 2004, QJRMS 

!     ------------------------------------------------------------------

USE PARKIND1,  ONLY : JPIM, JPRB

USE YOMCST   , ONLY : &  !setup/sucst.F90
   & RG      , &  ! = 9.80665_JPR
   & RD      , &  ! = 1000._JPRB*R/RMD =287.0597 
   & RCPD    , &  ! = 3.5_JPRB*RD = 1004.709
   & RETV    , &  ! = RV/RD-1.0_JPRB = 0.6077667
   & RLVTT        

USE YOEPHY   , ONLY : TEPHY

USE YOETHF   , ONLY :&!source1c/suphec1c.F90, setup/sucst.F90
   & RVTMP2       ! = RCPV/RCPD-1.0_JPRB
                  !   RCPV=4._JPRB *RV = 1846.100

USE PARPHY   , ONLY : RKAP,REPDU2
USE YOEVDF   , ONLY : TVDF

USE SPP_MOD  , ONLY : YSPP,YSPP_CONFIG


IMPLICIT NONE


!*     0.1    GLOBAL VARIABLES 

TYPE(TVDF)        ,INTENT(IN) :: YDVDF
TYPE(TEPHY)       ,INTENT(IN) :: YDEPHY
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
REAL(KIND=JPRB)   ,INTENT(IN) :: PTMST
REAL(KIND=JPRB)   ,INTENT(IN) :: PUM1      (KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PVM1      (KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PTM1      (KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PQM1      (KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PLM1      (KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PIM1      (KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PSLGM1    (KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PQTM1     (KLON,KLEV)
!REAL(KIND=JPRB)   ,INTENT(IN) :: PAM1      (KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PKMFL     (KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: PKHFL     (KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: PKQFL     (KLON)
!REAL(KIND=JPRB)   ,INTENT(IN) :: PMCU      (KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: PAPHM1    (KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN) :: PAPM1     (KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PGEOM1    (KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PGEOH    (KLON,0:KLEV)
!INTEGER(KIND=JPIM),INTENT(IN) :: KVARTOP   (KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: PZINV     (KLON)
!INTEGER(KIND=JPIM),INTENT(IN) :: KHPBL     (KLON)
INTEGER(KIND=JPIM),INTENT(IN) :: KPBLTYPE  (KLON)
INTEGER(KIND=JPIM),INTENT(IN) :: KDRAFT
!REAL(KIND=JPRB)   ,INTENT(IN) :: PRICUI    (KLON)
REAL(KIND=JPRB)   ,INTENT(IN) :: PQSVAR    (KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PDQSDTEMP (KLON,KLEV)
!INTEGER(KIND=JPIM),INTENT(IN) :: KFLDX2
!INTEGER(KIND=JPIM),INTENT(IN) :: KLEVX
!INTEGER(KIND=JPIM),INTENT(IN) :: KFLDX
!cstep/GL --- variables associated with dualM-TKE scheme
REAL(KIND=JPRB)   ,INTENT(IN) ::  ZBUOY_COR (KLON,0:KLEV)  
REAL(KIND=JPRB)   ,INTENT(IN) ::  ZEN_CASC  (KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) ::  ZWU       (KLON,0:KLEV)
!cstep/GL ---------------------------------------------

REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP2DSPP(KLON,YSPP%N2D)
INTEGER(KIND=JPIM),INTENT(IN)    :: KGFL_EZDIAG
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEZDIAG(KLON,KLEV,KGFL_EZDIAG) 

! DIAGNOSTIC OUTPUT
!REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCFM   (KLON,KLEV)
!REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCFH   (KLON,KLEV)
!REAL(KIND=JPRB)   ,INTENT(INOUT) :: PKH    (KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTKE   (KLON,KLEV)
!NB PMFLX is IN only
REAL(KIND=JPRB)   ,INTENT(IN) :: PMFLX(KLON,0:KLEV,KDRAFT)
REAL(KIND=JPRB)   ,INTENT(OUT) :: ZLENGTH_M(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT) :: ZLENGTH_H(KLON,KLEV)
!REAL(KIND=JPRB)   ,INTENT(INOUT) ::  PEXTR2(KLON,KFLDX2)



!*            LOCAL STORAGE
!             ----- -------

REAL(KIND=JPRB) ::&
           & ZDU2   (KLON) ,&! DELTA U^2 + DELTA V^2
           & ZMGEOM (KLON) ,&! G DELTA HEIGHT
           & ZUST   (KLON) ,&! FRICTION VELOCITY
           & ZDSL          ,&! DELTA LIQUID STATIC ENERGY
           & ZDQT          ,&! DELTA TOTAL SPECIFIC HUMIDITY 
           & ZKHVFL (KLON) ,&! SURFACE WTHV FLUX
           & ZWS    (KLON)     ! CHARACTERISTIC TURBULENT VELOCITY SCALE W_TURB

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: ILEVM1        ,&! HELP VARIABLE = KLEV -1       
          &  JK            ,&! COUNTER FOR VERTICAL LEVELS
          &  JL                ! COUNTER FOR HORIZONTAL POSITION

!     LOCAL REAL SCALARS
REAL(KIND=JPRB) ::    Z2GEOMF       ,&! HELP VARIABLE TO INTERPOLATE GEOPOTENTIAL HEIGHT
           & ZCFNC1        ,&! DIFFUSION COEFFICIENT
           & ZRG           ,&! = 1.0_JPRB/RG
           & ZCONS2        ,&! ZCONS2  = 0.5_JPRB * RKAP / RG
           & ZCONS10       ,&! = ZTPFAC1 * PTMST * RG**2 / ( 0.5_JPRB * RD )
           & ZCONS13       ,&! = 1/3, POWER IN TURB CHAR VEL SCALE
           & ZCONS23       ,&! = 0.5_JPRB * RKAP * RLAM / RG 
           & ZEPS          ,&! = 1.E-10_JPRB     !MINIMUM HEAT FLUX VALUE
           & ZTPFAC1           !FACTOR FOR TIME STEP WEIGHTING IN *VDF....*
                               !phys_ec/suvdf.F90:RVDIFTS=1.5_JPRB

#include "surf_inq.h"



REAL(KIND=JPRB) ::      ZTKEVN (KLON,KLEV)   !UPDATED TKE BEFORE TKE DIFFUSION
REAL(KIND=JPRB) ::      ZTKEOLD (KLON,KLEV)   !input tke

!Cgeert 20070802
REAL(KIND=JPRB) ::      ZCORRI (KLON,KLEV)   !CORRECTION  FACTOR BEFORE RI used in integral length scale


 
REAL(KIND=JPRB) ::  ZC_0    ,&! = 3.75 ,              C_0, EQ. (4)
         & ZC_D    ,&! = 1.0_JPRB/(ZC_0**2),    DISSIPATION LENGTH SCALE , BELOW EQ. (5)
         & ZC_N    ,&! = 1.0_JPRB/(ZC_0**0.5),  C_N = C_D^0.25 = 1/C_0^0.5, EQ. (A9), 
                          !                       SCALING NEUTRAL LENGTH SCALE AT SURFACE
         & ZWSTF   ,&! = 0.2                 FACTOR IN EQ. (4)
!        & ZL_INF  ,&! = 75.                 L_INFINITY, EQ. (B.3)
         & ZTKEMIN ,&! 1e-3                  MINIMUM THRESHOLD FOR TKE
!
         & ZHEIGHT    ,&!      HEIGHT AT HALF LEVELS
         & ZNBRUNT    ,&!      D/DZ (G/THV) THETAV = N_BV**2 , Brunt-Vaisala freq. squared
         & ZTKESQ     ,&!      SQUARE ROOT TKE
!
         & ZMIXCH  ,&!      SEE EQS (B1) and (B2)
         & ZMIXCM  ,&!      SEE EQS (B1) and (B2)
         & ZKTEST  ,&!      TEST IF TKE (T+1) > 0
         & ZDISL   ,&!      MODIFIED DISSIPATION TO SOLVE TKE
!
!        next block, see Eqs. (2) and (3) LH04
!
!         & ZLENGTH_M (KLON,KLEV)   , & ! LENGTH SCALE FOR MOMENTUM, l_m 
!         & ZLENGTH_H (KLON,KLEV)   , & ! LENGTH SCALE FOR HEAT, l_h
         & ZBUOY     (KLON,KLEV)   ,&! VERTICAL GRADIENT: D/DZ (G/THV) THETAV
! Cgeert testing array
         & ZBUOY2     (KLON,KLEV)   ,&! VERTICAL GRADIENT: D/DZ (G/THV) THETAV
         & ZSHEAR    (KLON,KLEV)   ,&! VERTICAL GRADIENT: D/DZ^2 (U**2 + V**2)
         & ZKH       (KLON,KLEV)   ,&! EDDY DIFFUSIVITY FOR HEAT, K_h
         & ZKM       (KLON,KLEV)   ,&! EDDY DIFFUSIVITY FOR MOMENTUM, K_m
         & ZCDUM     (KLON,KLEV)   ,&! MODIFIED EDDY DIFFUSIVITY FOR TKE, PROP to 2 * K_m
!
         & ZTHVM1    (KLON,KLEV)   ,&! VIRTUAL POTENTIAL TEMPERATURE AT T-1
         & ZTVIRDIF                ,&!
         & ZTVIRDIFCHECK(KLON,KLEV),&! bewaar voor check
         & ZKAPPA                  ,&! !R_d/c_pd  , NEEDED FOR EXNER FUNCTION
         & ZTVM1     (KLON,KLEV)   ,&! VIRTUAL TEMPERATURE
         & ZHH       (KLON,KLEV)   ,&! HALF LEVEL HEIGHT
         & ZRI2      (KLON,KLEV)   ,&! 
!CGeert
         & ZCLOUDTKE  (KLON,KLEV)  ,&
         & ZSTRESFAC (KLON) 

REAL(KIND=JPRB) ::  ZMIXUPH   (KLON,KLEV) , ZMIXDWH   (KLON,KLEV) , ZMIXUPM   (KLON,KLEV) , & 
             & ZMIXDWM   (KLON,KLEV) , ZMIXQUADM (KLON,KLEV) , ZMIXQUADH (KLON,KLEV), &
             &  PCFM(KLON,KLEV) , PCFH(KLON,KLEV) , PKH(KLON,KLEV)

INTEGER ::     ITOP    , ITOPP1
REAL(KIND=JPRB) ::      ZEDIF     (KLON,KLEV) , ZTCOE (KLON)
REAL(KIND=JPRB) ::      ZEBSM     (KLON,KLEV)
REAL(KIND=JPRB) ::      ZTPFAC2 , ZQDP     , ZDISC    , ZFAC    , ZTPFAC3  , &
             & ZZB 

REAL(KIND=JPRB) ::      ZPRANDT , PI2      , ZALPHA_C , ZALPHA_N, RHO      , &
             & DZ      , ZRI_FUL  , XAR      , X       , ZB       , DZM , X2, &
             & DZH     , DZM_INT  , DZH_INT  , ZAVE 
REAL(KIND=JPRB) ::      ZLENGTHRH,ZLENGTHRM, ZZ2GEO   , ZMIX    , &
!            & ZC_H     , & ! FACTOR Ch in Eq. (B1)
             & ZC_M     , & ! FACTOR Cm in Eq. (B2)
             & ZMIX2    , &
             & ZMIXH   , ZMIXM    , RPOW     , ZMIXH_STABLE, &
             & ZMIXM_STABLE , ZC_M_MAX
	
	     !CGeert
REAL(KIND=JPRB) ::     ZDQDVAR
	

REAL(KIND=JPRB) ::      ZHFAC, ZIFAC, ZTMP
REAL(KIND=JPRB) ::      ZTHM1                 !POTENTIAL TEMPERATURE 
REAL(KIND=JPRB) ::      ZVIRT                 !WATER PART OF VIRTUAL (POTENTIAL) TEMPERATURE
REAL(KIND=JPRB) ::      ZPRES_0               !REFERENCE PRESSURE FOR EXNER FUNCTION
REAL(KIND=JPRB) ::      ZCFRAC
REAL(KIND=JPRB) ::      ZALFA_TM, ZALFA_QM  , ZALFA_TD, ZALFA_QD , ZALFA_T , ZALFA_Q

LOGICAL ::         LLRICU
LOGICAL ::         EN_CASC

LOGICAL ::         LWINDADJ
REAL(KIND=JPRB) :: ZREPUST,  ZCFHNEW
REAL(KIND=JPRB) :: ZWECUTOP(KLON)

! SPP
REAL(KIND=JPRB) :: ZFAC_TWO_COEF(KLON),ZC_H(KLON),ZL_INF(KLON),ZMU,ZVAL
INTEGER(KIND=JPIM) :: JKO,JKE

ASSOCIATE(YSURF=>YDEPHY%YSURF,RLAM=>YDVDF%RLAM,RVDIFTS=>YDVDF%RVDIFTS, &
  & RFAC_TWO_COEF=>YDVDF%RFAC_TWO_COEF,RZC_H=>YDVDF%RZC_H,RZL_INF=>YDVDF%RZL_INF)

CALL SURF_INQ(YSURF,PREPUST=ZREPUST)



!     ------------------------------------------------------------------
!*       1.     INITIALIZE CONSTANTS
!               ---------- ----------

LWINDADJ = .TRUE.
! LWINDADJ = .FALSE.

ZRG       = 1.0_JPRB/RG

ZTPFAC1 = RVDIFTS
ZTPFAC2 = 1.0_JPRB / ZTPFAC1  ! = 2/3
ZTPFAC3 = 1.0_JPRB - ZTPFAC2  ! = 1/3

ZCONS2  = 0.5_JPRB * RKAP / RG
ZCONS10 = ZTPFAC1 * PTMST * RG**2 / ( 0.5_JPRB * RD )
ZCONS13 = 1.0_JPRB / 3._JPRB
ZCONS23 = ZCONS2 * RLAM

ZEPS    = 1.E-10_JPRB   
ZKAPPA  = RD / RCPD   
ZPRES_0 = 100000._JPRB

ILEVM1  = KLEV - 1

! constants for tke scheme, described by Lenderink and Holtslag QJRMS 2004

ZC_0    = 3.75_JPRB
ZC_D    = 1.0_JPRB/(ZC_0**2)
ZC_N    = 1.0_JPRB/(ZC_0**0.5_JPRB)

IF (LWINDADJ) THEN
  ZC_N    = 1.2_JPRB*ZC_N
ENDIF

ZWSTF   = 0.2_JPRB
ZTKEMIN = 1E-3_JPRB

ZPRANDT = 0.6_JPRB                             ! see LH04 page 3411, from Cuxart et al. 2000, 
PI2     = 2.0_JPRB * ATAN(1._JPRB)             ! 0.5 * pi
ZALPHA_N = ZC_N *RKAP                          ! alpha_n = c_n * kappa, Eq. (11)
ZALPHA_C = ZC_N *RKAP * 3._JPRB / ZPRANDT      ! alpha_c = 3 * c_n * kappa / Pr , p. 3411 
ZB       = 4._JPRB                             !factor b, Eq. (10)

!CGL different values
!ZC_H     = 0.20_JPRB                    !c_h, Eq. (B.1)
ZC_H     = 0.15_JPRB                    !c_h, Eq. (B.1)
ZC_M_MAX = 2._JPRB                     !limiter for c_m

IF (LWINDADJ) THEN
  ZC_M_MAX = 3._JPRB 
ENDIF

!ZL_INF   = 75._JPRB                    != l_infinity, Eq. (B.3) 
!ZL_INF   = 100._JPRB                    != l_infinity, Eq. (B.3) 
                                 
!FAC_TWO_COEF = 2.0_JPRB            !Eddy diffusivity for TKE diffusion is 2K_m

! TKE is defined at half levels.
! note I picked definition of levels from ECHAM4 code. This means that
! TKE(klev) is the surface value (usual half flux level klev+1) as also
! is done in ECHAM4. So be carefull !
! anyway, in ECMWF PCFM etc is also defined at the TKE half (1) leveling.
! So I take this numbering for the length scales, buoyancy, shear etc.
!
!     (1)        (2)      (3)
!     TKE         M     Flux levels
!
!     k-1                  k         k-1/2
!                 k                   k
!      k                  k+1        k+1/2
!                k+1
!

LLRICU = .TRUE.   ! switch for top-entrainment efficiency closure using Ri^cu at cumulus PBL top
LLRICU = .FALSE.
EN_CASC = .TRUE.

! Init SPP
IF (YSPP_CONFIG%LSPP.AND.YSPP_CONFIG%LPERT_RFAC_TWOC) THEN
  IF (YSPP_CONFIG%LLNN_MEAN1.OR.YSPP_CONFIG%LLNN_MEAN1_RFAC_TWOC) THEN
    ZMU = -0.5_JPRB * (YSPP_CONFIG%CMPERT_RFAC_TWOC * YSPP_CONFIG%SDEV)**2
  ELSE
    ZMU = 0._JPRB
  ENDIF
  DO JL=KIDIA,KFDIA
    ZVAL = RFAC_TWO_COEF*EXP(ZMU+YSPP_CONFIG%CMPERT_RFAC_TWOC*PGP2DSPP(JL,YSPP%MP_RFAC_TWOC))
    ZFAC_TWO_COEF(JL) = MAX(YSPP_CONFIG%CLIP_RFAC_TWOC(1),MIN(ZVAL,YSPP_CONFIG%CLIP_RFAC_TWOC(2)))
  ENDDO
  IF ( YSPP_CONFIG%IEZDIAG_POS > 0 ) THEN
   JKO=2*YSPP%MP_RFAC_TWOC-1
   JKE=2*YSPP%MP_RFAC_TWOC
   DO JL=KIDIA,KFDIA
    PEZDIAG(JL,JKO,YSPP_CONFIG%IEZDIAG_POS) = PGP2DSPP(JL,YSPP%MP_RFAC_TWOC)
    PEZDIAG(JL,JKE,YSPP_CONFIG%IEZDIAG_POS) = ZFAC_TWO_COEF(JL)
   ENDDO
  ENDIF
ELSE
  DO JL=KIDIA,KFDIA
    ZFAC_TWO_COEF(JL) = RFAC_TWO_COEF
  ENDDO
ENDIF

IF (YSPP_CONFIG%LSPP.AND.YSPP_CONFIG%LPERT_RZC_H) THEN
  IF (YSPP_CONFIG%LLNN_MEAN1.OR.YSPP_CONFIG%LLNN_MEAN1_RZC_H) THEN
    ZMU = -0.5_JPRB * (YSPP_CONFIG%CMPERT_RZC_H * YSPP_CONFIG%SDEV)**2
  ELSE
    ZMU = 0._JPRB
  ENDIF
  DO JL=KIDIA,KFDIA
    ZVAL = RZC_H*EXP(ZMU+YSPP_CONFIG%CMPERT_RZC_H*PGP2DSPP(JL,YSPP%MP_RZC_H))
    ZC_H(JL) = MAX(YSPP_CONFIG%CLIP_RZC_H(1),MIN(ZVAL,YSPP_CONFIG%CLIP_RZC_H(2)))
  ENDDO
  IF ( YSPP_CONFIG%IEZDIAG_POS > 0 ) THEN
   JKO=2*YSPP%MP_RZC_H-1
   JKE=2*YSPP%MP_RZC_H
   DO JL=KIDIA,KFDIA
    PEZDIAG(JL,JKO,YSPP_CONFIG%IEZDIAG_POS) = PGP2DSPP(JL,YSPP%MP_RZC_H)
    PEZDIAG(JL,JKE,YSPP_CONFIG%IEZDIAG_POS) = ZC_H(JL)
   ENDDO
  ENDIF
ELSE
  DO JL=KIDIA,KFDIA
    ZC_H(JL) = RZC_H
  ENDDO
ENDIF

IF (YSPP_CONFIG%LSPP.AND.YSPP_CONFIG%LPERT_RZL_INF) THEN
  IF (YSPP_CONFIG%LLNN_MEAN1.OR.YSPP_CONFIG%LLNN_MEAN1_RZL_INF) THEN
    ZMU = -0.5_JPRB * (YSPP_CONFIG%CMPERT_RZL_INF * YSPP_CONFIG%SDEV)**2
  ELSE
    ZMU = 0._JPRB
  ENDIF
  DO JL=KIDIA,KFDIA
    ZVAL = RZL_INF*EXP(ZMU+YSPP_CONFIG%CMPERT_RZL_INF*PGP2DSPP(JL,YSPP%MP_RZL_INF))
    ZL_INF(JL) = MAX(YSPP_CONFIG%CLIP_RZL_INF(1),MIN(ZVAL,YSPP_CONFIG%CLIP_RZL_INF(2)))
  ENDDO
  IF ( YSPP_CONFIG%IEZDIAG_POS > 0 ) THEN
   JKO=2*YSPP%MP_RZL_INF-1
   JKE=2*YSPP%MP_RZL_INF
   DO JL=KIDIA,KFDIA
    PEZDIAG(JL,JKO,YSPP_CONFIG%IEZDIAG_POS) = PGP2DSPP(JL,YSPP%MP_RZL_INF)
    PEZDIAG(JL,JKE,YSPP_CONFIG%IEZDIAG_POS) = ZL_INF(JL)
   ENDDO
  ENDIF
ELSE
  DO JL=KIDIA,KFDIA
    ZL_INF(JL) = RZL_INF
  ENDDO
ENDIF

!        1.1    PREPARE FOR MIXED LAYER DEPTH
!               ------- --- ----- ----- -----

DO JL = KIDIA,KFDIA
  ZUST   (JL)       = MAX ( SQRT (PKMFL(JL)) , ZREPUST )
  ZKHVFL (JL)       = PKHFL (JL) + RETV * PTM1 (JL,KLEV) * PKQFL(JL)
ENDDO



!        1.2    PREPARE SCALING COEFFICIENTS FOR MIXED LAYER
!               ------- ------- ------------ --- ----- -----

DO JL=KIDIA,KFDIA
  IF (ZKHVFL (JL)  <  0.0_JPRB) THEN   !UNSTABLE STRATIFICATION
    ZWS    (JL) = ( - ZKHVFL (JL) * RG / PTM1 (JL,KLEV) * PZINV (JL)) ** ZCONS13
  ELSE
    ZWS    (JL) = 0.0_JPRB
  ENDIF

  DO JK=1,KLEV
    ZKM       (JL,JK) = 0.0_JPRB
    ZKH       (JL,JK) = 0.0_JPRB
    ZBUOY     (JL,JK) = 0.0_JPRB
    ZSHEAR    (JL,JK) = 0.0_JPRB
    ZLENGTH_M (JL,JK) = 0.0_JPRB
    ZLENGTH_H (JL,JK) = 0.0_JPRB
    ZTKEOLD   (JL,JK) = PTKE(JL,JK)
    ZMIXQUADM   (JL,JK) = 0.0_JPRB
    ZMIXQUADH   (JL,JK) = 0.0_JPRB
  ENDDO

! set new surface boundary condition TKE

  PTKE   (JL,KLEV) = ZC_0 *ZUST (JL)** 2 + ZWSTF * ZWS(JL)**2
  ZTKEVN (JL,KLEV) = PTKE (JL,KLEV)


! note: global values PCFM, PCFH,PKH already computed on lowest model level

  DO JK = 1 , KLEV-1
    PCFM (JL,JK) = 0.0_JPRB
    PCFH (JL,JK) = 0.0_JPRB
    PKH  (JL,JK) = 0.0_JPRB
  ENDDO

ENDDO

  
!  compute shear and stability to be used in other parts of the code
       
DO JK = KLEV , 1 , -1
  DO JL = KIDIA , KFDIA
    ZVIRT          = (1.0_JPRB + RETV * PQM1 (JL,JK) - PLM1 (JL,JK) - PIM1(JL,JK))
    ZTHM1          = PTM1  (JL,JK) * (ZPRES_0 / PAPM1 (JL,JK) ) ** ZKAPPA
    ZTVM1  (JL,JK) = PTM1  (JL,JK) * ZVIRT
    ZTHVM1 (JL,JK) = ZTHM1 * ZVIRT
!wc
! watch out geopotential corrected for geopotential surface (not msl)
    Z2GEOMF        = (PGEOM1 (JL,JK)-PGEOH(JL,KLEV)) + (PGEOM1 (JL,MIN(JK+1,KLEV))-PGEOH(JL,KLEV))
    ZHH    (JL,JK) = Z2GEOMF / (RG * 2.0_JPRB) 
  ENDDO
ENDDO  

DO JL = KIDIA , KFDIA
  ZHH (JL,KLEV) = ZEPS
ENDDO 

    !CGeert : recompute cloud fraction for testing purposes

DO JK = 1, KLEV
  DO JL = KIDIA , KFDIA
    ZDQDVAR = PQM1(JL,JK) + PLM1 (JL,JK) + PIM1 (JL,JK) - PQSVAR (JL,JK)
    ZDQDVAR =  ZDQDVAR / (MAX(0.04_JPRB*PQSVAR (JL,JK),1E-8_JPRB))    
    ZCLOUDTKE(JL,JK) = 0.5_JPRB + 0.36_JPRB * ATAN( 1.55_JPRB * ZDQDVAR ) 
    ZCLOUDTKE(JL,JK) = MIN(MAX(ZCLOUDTKE(JL,JK),0.0_JPRB),1.0_JPRB)
!    write (188,'(10f10.3)')  PGEOM1(JL,JK)*ZRG, ZCLOUDTKE(JL,JK), ZDQDVAR, 1000*PQSVAR (JL,JK),1000*( PQM1(JL,JK) + PLM1 (JL,JK) + PIM1 (JL,JK) - PQSVAR (JL,JK))
    IF ( PLM1 (JL,JK) + PIM1 (JL,JK) .LE. 0.0_JPRB )  ZCLOUDTKE(JL,JK) = 0._JPRB
!    write (911,'(8f10.2)') PGEOM1 (JL,JK)/10., ZCLOUDTKE(JL,JK), 1000*(PQM1(JL,JK) + PLM1 (JL,JK) + PIM1 (JL,JK)), 1000*PQSVAR (JL,JK), ZDQDVAR
  ENDDO
ENDDO    

!    write (188,*)  
!    write (188,*)   

! CGL 20151217
! this determines the stress dependent modification of length scales   
DO JL = KIDIA , KFDIA
  ZSTRESFAC(JL) = (1.0_JPRB -EXP(-(ZUST (JL)/1.0_JPRB)**4.0_JPRB ) ) 
ENDDO 
  
DO JK = KLEV - 1 , 1 , -1
  DO JL = KIDIA , KFDIA

    ZDU2 (JL)  = MAX (REPDU2 , ( PUM1 (JL,JK) - PUM1 (JL,JK+1) )**2 +&
                       &       ( PVM1 (JL,JK) - PVM1 (JL,JK+1) )**2 )
       
    ZMGEOM(JL) = PGEOM1 (JL,JK) - PGEOM1 (JL,JK+1)
  
    ZDSL    = (PSLGM1(JL,JK) - PSLGM1(JL,JK+1)) / (ZMGEOM (JL) / RG)
    ZDQT    = (PQTM1(JL,JK) - PQTM1(JL,JK+1)) / (ZMGEOM (JL) / RG)

!
!  Compute thermodynamic coefficients, neglect difference full/half levels
!  See B. Stevens, Entrainment in stratocumulus mixed layers, QJRMS
!
    ZALFA_TM = ( 1.0_JPRB + (1.0_JPRB+RETV) * PQSVAR(JL,JK) - PQTM1(JL,JK)&
         &    + (1.0_JPRB+RETV) * PTM1(JL,JK) * PDQSDTEMP(JL,JK)      )&
         & / ( 1.0_JPRB + RLVTT/RCPD * PDQSDTEMP(JL,JK) )
    ZALFA_QM = RLVTT * ZALFA_TM / (RCPD * PTM1(JL,JK)) - 1.0_JPRB 

    ZALFA_TD = 1.0_JPRB + RETV *  PQTM1(JL,JK)
    ZALFA_QD = RETV 
    
    !ZCFRAC  = (PAM1(JL,JK) + PAM1(JL,JK+1)) / 2.0_JPRB 
    !CGeert : recompute cloud fraction 
    ZCFRAC  = (ZCLOUDTKE(JL,JK) +ZCLOUDTKE(JL,JK+1)) / 2.0_JPRB
    
    ZALFA_T = (ZALFA_TD * (1.0_JPRB - ZCFRAC) + ZCFRAC * ZALFA_TM) / RCPD
    ZALFA_Q = (ZALFA_QD * (1.0_JPRB - ZCFRAC) + ZCFRAC * ZALFA_QM) * PTM1(JL,JK)

    ZTVIRDIF       = ZALFA_T * ZDSL + ZALFA_Q * ZDQT
    ZTVIRDIFCHECK(JL,JK)=ZTVIRDIF
    ZBUOY  (JL,JK) = 2.0_JPRB* ZTVIRDIF * RG / (ZTHVM1 (JL,JK) + ZTHVM1 (JL,JK+1))
    ZSHEAR (JL,JK) =  ZDU2 (JL) / (ZMGEOM(JL) / RG)**2
    ZRI2   (JL,JK) = ZBUOY (JL,JK) / ZSHEAR(JL,JK)

    ! CGL re-compute stability parameters
    ! ZBUOY2 is only used in stability parameters for length scale formulation !
    ! correction is computer is vdfhghtn
    ZBUOY2(JL,JK) = ZBUOY(JL,JK)   - ZBUOY_COR (JL,JK)  * 1.
!CGL check wether jk is from 0:klev

     ! CGL added some shear to ZRI2 formulation !
      !cstep    write (6,*) 'cstep',JL,JK,ZSHEAR(JL,JK), ZWU(JL,JK)
    ZRI2   (JL,JK) = ZBUOY2 (JL,JK) / ( ZSHEAR(JL,JK)**2 + (ZWU(JL,JK)/RLAM)**4 )**0.5_JPRB

    IF (LWINDADJ) THEN
!     do a simple contribution of TKE to windshear to compute RI number
      Z2GEOMF = (PGEOM1 (JL,JK)-PGEOH(JL,KLEV)) + (PGEOM1 (JL,JK+1) - PGEOH(JL,KLEV))
      ZHEIGHT = Z2GEOMF / (RG * 2.0_JPRB)
      ZMIX2   = ZL_INF(JL) / (1.0_JPRB + ZL_INF(JL) / ( ZC_N * RKAP * ZHEIGHT))
      ZRI2   (JL,JK) = ZBUOY2 (JL,JK) /&
     &    ( ZSHEAR(JL,JK)**2 + (ZWU(JL,JK)/RLAM)**4 + ZSTRESFAC(JL)*PTKE(JL,JK)**2/ZMIX2**4)**0.5_JPRB

    ENDIF

  ENDDO
ENDDO

DO JL = KIDIA , KFDIA
  ZRI2 (JL,KLEV) = ZRI2 (JL,KLEV-1)
ENDDO


! new length scale formulation !

! bottom-up length scale
DO JL = KIDIA , KFDIA
  ZMIXUPH (JL, KLEV) = 0.0_JPRB
  ZMIXUPM (JL, KLEV) = 0.0_JPRB
ENDDO

XAR = ZB * PI2 * ZC_N * RKAP  / (ZALPHA_C * ZPRANDT - ZALPHA_N)

DO JK = KLEV-1, 2, -1             ! VERTICAL LOOP BOTTOM-UP
  DO JL=KIDIA,KFDIA

    RHO = 0.5_JPRB*(PAPHM1(JL,JK+1)+PAPHM1(JL,JK+2))&
          &      /(RD*ZTVM1(JL,JK+1))
    DZ = (PAPHM1(JL,JK+2) - PAPHM1(JL,JK+1))/(RG*RHO)  !check half levels
                                   !preshalf array index differs from ptke in this routine 
       
    ZRI_FUL = 0.5_JPRB * (ZRI2(JL,JK) + ZRI2(JL,MIN(JK+1,KLEV-1))) 
    X = XAR * ZRI_FUL

    DZM_INT = ( ZALPHA_C * ZPRANDT - ZALPHA_N ) / PI2 * DZ
    DZH_INT = ( ZALPHA_C - ZALPHA_N ) / PI2 * DZ

    IF (X.GT.0._JPRB) THEN
      DZM = ZALPHA_N*DZ - DZM_INT * X
      DZH = ZALPHA_N*DZ - DZH_INT * X 
    ELSE
      DZM = ZALPHA_N*DZ - DZM_INT * ATAN(X)
      DZH = ZALPHA_N*DZ - DZH_INT * ATAN(X)
    ENDIF

    ZMIXUPH(JL, JK) =  ZMIXUPH (JL, JK + 1) + DZH
    ZMIXUPH(JL, JK) = MAX(ZMIXUPH(JL, JK) , ZEPS )

    ZMIXUPM(JL, JK) =  ZMIXUPM (JL, JK + 1) + DZM
    ZMIXUPM(JL, JK) = MAX(ZMIXUPM(JL, JK) , ZEPS)

  ENDDO
ENDDO

!top-down length scale
DO JL=KIDIA,KFDIA                   
  ZMIXDWH(JL, 1) = 0.0_JPRB 
  ZMIXDWM(JL, 1) = 0.0_JPRB
ENDDO
!

DO JK = 2, KLEV-1            ! VERTICAL LOOP TOP, DOWN
  DO JL = KIDIA,KFDIA     

    RHO = 0.5_JPRB *(PAPHM1(JL,JK+1)+PAPHM1(JL,JK))&
             &    /(RD*(ZTVM1(JL,JK)))
    DZ = (PAPHM1(JL,JK+1) - PAPHM1(JL,JK))/(RG*RHO)
                   
    ZRI_FUL = 0.5_JPRB *(ZRI2(JL,JK) + ZRI2(JL,JK-1)) 

    X = XAR * ZRI_FUL
    X2 = X
    IF (LWINDADJ) THEN
      X2 = XAR * (ZRI_FUL - 0.3_JPRB*ZSTRESFAC(JL) )
    ENDIF

    DZM_INT = ( ZALPHA_C * ZPRANDT - ZALPHA_N ) / PI2 * DZ
    DZH_INT = ( ZALPHA_C - ZALPHA_N ) / PI2 * DZ

    IF (X.GT.0._JPRB) THEN
      DZH = ZALPHA_N*DZ - DZH_INT * X 
    ELSE
      DZH = ZALPHA_N*DZ - DZH_INT * ATAN(X) 
    ENDIF

    IF (X2.GT.0._JPRB) THEN
      DZM = ZALPHA_N*DZ - DZM_INT * X2
    ELSE
      DZM = ZALPHA_N*DZ - DZM_INT * ATAN(X2) 
    ENDIF

    ZMIXDWH(JL,JK) =  ZMIXDWH (JL,JK - 1)  +  DZH             
    ZMIXDWH(JL,JK) = MAX (ZMIXDWH(JL,JK),ZEPS) 

    ZMIXDWM(JL,JK) =  ZMIXDWM (JL,JK - 1)  +  DZM
    ZMIXDWM(JL,JK) = MAX (ZMIXDWM(JL,JK),ZEPS)    

!    COMPOSED LENGTH SCALE OF TOP-DOWN LS AND BOTTOM-UP LS

    ZAVE = 1.0_JPRB
!cstep : prevent division by zero

    ZLENGTHRH = 1.0_JPRB / ZMIXDWH(JL,JK)**ZAVE&
          &      + 1.0_JPRB / ZMIXUPH(JL,JK)**ZAVE
    ZMIXQUADH(JL,JK) = 1.0_JPRB / ZLENGTHRH**(1.0_JPRB/ZAVE)
    ZLENGTHRM = 1.0_JPRB / ZMIXDWM(JL,JK)**ZAVE&
          &   + 1.0_JPRB / ZMIXUPM(JL,JK)**ZAVE
    ZMIXQUADM(JL,JK) = 1.0_JPRB / ZLENGTHRM**(1.0_JPRB/ZAVE)

  ENDDO
ENDDO                                  ! END VERTICAL LOOP


!   END COMPUTATION FREE TURBULENCE LENGTH SCALE 
!   LENGTHSCALES ARE IN ZMIXQUADH (HEAT) AND ZMIXQUADM (MOMENTUM)
     

DO JK = 1, KLEV-1            
  DO JL = KIDIA,KFDIA    
!wc 
! watch out geopotential corrected for geopotential surface (not msl)

    Z2GEOMF = (PGEOM1 (JL,JK)-PGEOH(JL,KLEV)) + (PGEOM1 (JL,JK+1) - PGEOH(JL,KLEV))
    ZZ2GEO  = ZCONS2 * Z2GEOMF
    ZMIX    = ZZ2GEO / (1.0_JPRB + ZCONS23 * Z2GEOMF)

    ZC_M  = ZC_H(JL) * MAX(1.0_JPRB,MIN(1.0_JPRB +ZRI2(JL,JK),ZC_M_MAX))
    ZHEIGHT = Z2GEOMF / (RG * 2.0_JPRB)

    IF (LWINDADJ) THEN

     ! some adjustment to increase momentum mixing in strong wind conditions
     ! prandtl nr 1 -> 2
     ! enhanced mxixing (max. fact. ZIFAC) as decreasing function of height, height = poormans monin obhukov length ZTMP
     ! do only for strong wind conditions, with hight u* ( 1- exp(-(u*/1.0)^4)
     
     ZC_M  = ZC_H(JL) * MAX(1.0_JPRB,MIN(1.0_JPRB +2.0_JPRB*ZRI2(JL,JK),ZC_M_MAX))    
     ZDU2 (JL)  = MAX (REPDU2 , ( PUM1 (JL,KLEV-1) - PUM1 (JL,KLEV) )**2 +&
                       &        ( PVM1 (JL,KLEV-1) - PVM1 (JL,KLEV) )**2 )
     ZTMP = ZUST(JL)/MAX(ABS(ZBUOY(JL,KLEV-1)),0.00001_JPRB)*ZDU2(JL)**0.5/&
                       &       ((PGEOM1 (JL,KLEV-1) - PGEOM1 (JL,KLEV))/RG)   
     ZTMP = MIN(MAX(ZTMP,1._JPRB),1000._JPRB)   
     ZHFAC = EXP(-ZHEIGHT/ZTMP)
     ZC_M  = ZC_M * ( 1.0_JPRB + 2.0_JPRB*ZSTRESFAC(JL)*ZHFAC )
    ENDIF

    ZHEIGHT = Z2GEOMF / (RG * 2.0_JPRB) 

! limited to half the neutral value
    ZMIX2   = ZL_INF(JL) / (1.0_JPRB + ZL_INF(JL) / (0.5_JPRB * ZC_N *RKAP * ZHEIGHT))

    RPOW = 2.0_JPRB      ! NO PHYS.
                 ! ADDED FOR CONTINUITY BETWEEN UNSTABLE AND STABLE
    ZMIXH =  (ZMIX2**RPOW + ZMIXQUADH(JL,JK)**RPOW)**(1.0_JPRB/RPOW)
    ZMIXM =  (ZMIX2**RPOW + ZMIXQUADM(JL,JK)**RPOW)**(1.0_JPRB/RPOW)
	
    ZNBRUNT = ZBUOY(JL,JK)
   !cgeert
    ZNBRUNT = ZBUOY2(JL,JK)

    ZTKESQ  = SQRT(MAX(PTKE(JL,JK),ZTKEMIN))

    IF (ZNBRUNT.GT.0._JPRB) THEN
      ZMIXCH = 1.0_JPRB +  ZMIXH * SQRT ( ZNBRUNT ) / ( ZC_H(JL) * ZTKESQ )     
      ZMIXCM = 1.0_JPRB +  ZMIXM * SQRT ( ZNBRUNT ) / ( ZC_M * ZTKESQ )     
    ELSE
      ZMIXCH = 1.0_JPRB
      ZMIXCM = 1.0_JPRB 
    ENDIF        
        
    ZLENGTH_M (JL,JK) = ZMIXM / ZMIXCM
    ZLENGTH_H (JL,JK) = ZMIXH / ZMIXCH

!  IF ( JL .EQ. 1  .AND.  JK .GT. 35 ) THEN
!    write (186,*) ZHEIGHT, ZLENGTH_M (JL,JK), ZLENGTH_H (JL,JK), ZMIX2, ZC_H(JL) * ZTKESQ/SQRT( MAX(ZNBRUNT,1E-5) ), &
!   &    ZC_M * ZTKESQ/SQRT( MAX(ZNBRUNT,1E-5)), ZTKESQ, ZNBRUNT, ZMIXM, ZMIXCM, ZMIXDWM(JL,JK), ZRI2(JL,JK), ZMIXQUADM(JL,JK), &
!   &    (PUM1(JL,JK+1)**2 + PVM1(JL,JK+1)**2)**0.5, ZC_M, ZUST (JL), ZIFAC*ZHFAC, "V10 cleanup version test"
!  ENDIF


  ENDDO 
ENDDO

! compute K

DO JK=KLEV-1,1,-1
  DO JL=KIDIA,KFDIA
    
! ZCFNC1 is needed to get is fit into the original code. 
! preferably this should be done in the diffusion solver, and not here !
! but somehow this solver expects K/ZCFNC1 
 
! recompute these because vertical loop has been broken (needed for length scale)
    ZDU2   (JL) = MAX (REPDU2, (PUM1(JL,JK) - PUM1(JL,JK+1))**2 +&
                      &        (PVM1(JL,JK) - PVM1(JL,JK+1))**2 )
    ZMGEOM (JL) = PGEOM1 (JL,JK) - PGEOM1(JL,JK+1)
              
    ZCFNC1 = ZCONS10 * PAPHM1 (JL,JK+1)&
       & / ( ZMGEOM(JL)&
       & * ( ZTVM1 (JL,JK) + ZTVM1 (JL,JK+1)))
       
    ZTKESQ     = SQRT (MAX ( PTKE(JL,JK) , ZTKEMIN))
    ZKH(JL,JK) = ZTKESQ * ZLENGTH_H (JL,JK)
    ZKM(JL,JK) = ZTKESQ * ZLENGTH_M (JL,JK)
    
!CGL  ! add some mixing to help subcloud-cloud interaction.
!     100. is a lengthscale; typically 0.1 * cloud depth 
    ZKH(JL,JK) =  ZKH(JL,JK) +  50._JPRB*PMFLX(JL,JK,3)  
    ZKM(JL,JK) =  ZKM(JL,JK) +  50._JPRB*PMFLX(JL,JK,3)  
            
!    PCFM(JL,JK) = ZKM(JL,JK) * ZCFNC1   
!    PCFH(JL,JK) = ZKH(JL,JK) * ZCFNC1 
    PCFM(JL,JK) = ZKM(JL,JK)    
    PCFH(JL,JK) = ZKH(JL,JK) 
  
   
      
! a factor two to compensate for the pE transport term of TKE
! I am not entirely sure whether ZCDUM is similar to PCFM 
! check this in the diffusion solver

! picked from ECHAM4 code
!    ZCDUM (JL,JK) = ZFAC_TWO_COEF(JL) * PCFM (JL,JK)
    ZCDUM (JL,JK) = ZFAC_TWO_COEF(JL) * PCFM (JL,JK) * ZCFNC1

   PKH (JL,JK) = PCFH (JL,JK)
  ENDDO
ENDDO

!    write (182,*)  
!    write (182,*) 

     
!    
! do tke integration in time
!

! compute bouyancy flux, shear production, 
! and do first implicit timestepping (diffusion equation without diffusion solved implicitly
! as described in Brinkop and Roeckner as in ECHAM4). 
! stable, but not entirely optimal from a numerical point
! of view. explicit may be better, but more prone to instability
! a fractional timestepping might also improve things
! 
!
! cstep: in the following, BrinkopRoeckner (BrRo) E^t-1 is replaced by PTKE = E^t
!      : as a consequence, 2Del_t factor in BrRo is replaced by 1 Del_t!
!

DO JK=KLEV-1,1,-1  
  DO JL=KIDIA,KFDIA
    ZTKESQ  = SQRT (MAX ( PTKE(JL,JK) , ZTKEMIN ))        
    IF (EN_CASC) THEN
       ZZB     = ZSHEAR(JL,JK)* ZLENGTH_M (JL,JK)  - ZBUOY(JL,JK) * ZLENGTH_H (JL,JK)&
            & + ZEN_CASC(JL,JK)/ZTKESQ 
    ELSE
       ZZB     = ZSHEAR(JL,JK)* ZLENGTH_M (JL,JK)  - ZBUOY(JL,JK) * ZLENGTH_H (JL,JK)
    ENDIF

    ZDISL   =  ZLENGTH_M (JL,JK) / (PTMST * ZC_D)

    ZKTEST= 1.0_JPRB + ( ZZB * PTMST + SQRT ( PTKE (JL,JK)) * 2.0_JPRB) / ZDISL
    IF (ZKTEST.LE. 1.0_JPRB) THEN
      ZTKEVN(JL,JK) = ZTKEMIN
    ELSE
      ZTKEVN(JL,JK) = MAX (ZTKEMIN , (ZDISL * (SQRT (ZKTEST) - 1.0_JPRB )) **2 )
    ENDIF
!    write (186,'(12f10.3)') PGEOM1(JL,JK)*ZRG, ZBUOY(JL,JK)*PCFH(JL,JK)*1e4, 1e4*ZEN_CASC(JL,JK), PCFH(JL,JK), ZCLOUDTKE(JL,JK), 1000*ZBUOY(JL,JK), 1000*ZBUOY2(JL,JK)
	
  ENDDO
ENDDO

!    write (186,'(10f8.3)') 
!    write (186,'(10f8.3)')  

! finally do the diffusion of tke according to ECHAM4 code

ITOP   = 1
ITOPP1 = ITOP + 1
      
DO JK = ITOP, KLEV
  DO JL = KIDIA, KFDIA
    ZEDIF (JL,JK) = ZTPFAC2 * ZTKEVN (JL,JK)
  ENDDO
ENDDO

DO JL = KIDIA , KFDIA
  ZCDUM(JL,KLEV) =  ZCDUM(JL,KLEV-1)
  ZQDP            = 1.0_JPRB / (PAPM1 (JL,ITOPP1) - PAPM1 (JL,ITOP))
  ZTCOE (JL)      = (ZCDUM (JL,ITOP) + ZCDUM (JL,ITOPP1)) * 0.5_JPRB
  ZDISC           = 1.0_JPRB / ( 1.0_JPRB + ZTCOE (JL) * ZQDP)
  ZEBSM (JL,ITOP) = ZDISC * ZTCOE (JL) * ZQDP
  ZEDIF (JL,ITOP) = ZDISC * ZEDIF (JL,ITOP)
ENDDO
 
DO JK = ITOPP1 , KLEV-2
  DO JL = KIDIA , KFDIA
    ZQDP         = 1.0_JPRB / (PAPM1 (JL,JK+1) - PAPM1 (JL,JK))
    ZFAC         = ZTCOE (JL) * ZQDP
    ZTCOE(JL)    = (ZCDUM (JL,JK+1) + ZCDUM (JL,JK)) * 0.5_JPRB
    ZDISC        = 1.0_JPRB / (1.0_JPRB + ZFAC * ( 1.0_JPRB - ZEBSM (JL,JK-1)) + ZTCOE(JL) * ZQDP)
    ZEBSM(JL,JK) = ZDISC * ZTCOE (JL) * ZQDP
    ZEDIF(JL,JK) = ZDISC * (ZEDIF(JL,JK) + ZFAC * ZEDIF (JL,JK-1))
  ENDDO
ENDDO
 
DO JL=KIDIA,KFDIA
  ZQDP              = 1.0_JPRB / ( PAPM1(JL,KLEV) - PAPM1 (JL,KLEV-1) )
  ZFAC              = ZTCOE (JL) * ZQDP
  ZTCOE(JL)         = (ZCDUM (JL,KLEV) + ZCDUM (JL,KLEV-1) ) * 0.5_JPRB
  ZDISC             = 1.0_JPRB / (1.0_JPRB + ZFAC * ( 1.0_JPRB - ZEBSM (JL,KLEV-2)) +&
                    &     ZTCOE(JL) *ZQDP)
  ZEDIF (JL,KLEV-1) = ZDISC * (ZTCOE (JL)&
     &   * ZQDP * ZEDIF (JL,KLEV) + ZEDIF (JL,KLEV-1) + ZFAC * ZEDIF (JL,KLEV-2))
ENDDO
 
DO JK=KLEV-2,ITOP,-1
  DO JL=KIDIA,KFDIA
    ZEDIF (JL,JK) = ZEDIF (JL,JK) + ZEBSM (JL,JK) * ZEDIF (JL,JK+1)
  ENDDO
ENDDO

!     --------------------------------------------------------------------
!     TIME INTEGRATION OF TURBULENT KINETIC ENERGY AND CHECK
!     --------------------------------------------------------------------

DO JK = ITOP, KLEV-1
  DO JL = KIDIA,KFDIA
    PTKE (JL,JK) = ZEDIF (JL,JK) + ZTPFAC3 * ZTKEVN (JL,JK)
    IF (PTKE (JL,JK) <= 0.0_JPRB ) THEN
      PRINT *, 'TKE IS NEGATIVE = ', PTKE(JL,JK)
      PTKE(JL,JK) = ZTKEMIN
    ENDIF
  ENDDO
ENDDO


DO JK = KLEV,1,-1
  DO JL = KIDIA,KFDIA
    ZTKESQ  = SQRT(MAX(PTKE(JL,JK),ZTKEMIN))
  ENDDO	
ENDDO

END ASSOCIATE
RETURN
END SUBROUTINE VDFEXCUHL
