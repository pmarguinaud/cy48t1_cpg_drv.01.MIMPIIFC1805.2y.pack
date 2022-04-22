SUBROUTINE AOD_GLOMAP_TL&
 & (KFLEVG, KDLEN  , KLEN  , KMXCOUNT, &
 & YDGP_TL,&
 & YDGP5,&
 & PAEROD_PRIME)

!*** *AOD_GLOMAP_TL* Tangent Linear Operator routine for aerosol optical depth observations with the GLOMAP model.
!
!**   INTERFACE.
!     ----------
!
!     WHERE KFLEVG    = Number of vertical levels
!           KDLEN     = First dimension: Max number of observations, this set (INPUT)
!           KLEN      = Number of aerosol optical depth obs (INPUT)
!           KMXCOUNT  = Number of channels (INPUT)
!           YDGP_TL   = Perturbs in variables at observation points, model levels (INPUT)
!           YDGP5     = Nonlinear trajectory values at observation points, model levels (INPUT)
!
!           PAEROD_PRIME      = Aerosol optical depth perturbation at different wavelengths (OUTPUT)
!
!
!**   AUTHOR.
!     -------
!        Melanie Ades        *ECMWF*
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2018-08-01
!
!
!     PURPOSE.
!     --------
!     - derive tangent linear aerosol optical thickness from profiles of mode number concentration and component mass mixing ratio 
!
!-----------------------------------------------------------------------

USE PARKIND1       , ONLY : JPIM, JPRB
USE YOMHOOK        , ONLY : LHOOK, DR_HOOK
USE YOMCST         , ONLY : RG, YRCST
USE YOETHF   , ONLY : YRTHF
USE GOM_PLUS       , ONLY : TYPE_GOM_PLUS, IH
USE YOM_GRIB_CODES, ONLY : NGRBAERLG
USE UKCA_MODE_SETUP, ONLY: NMODES, NCP, COMPONENT, MODE
USE UKCA_RADAER_MODETYPE, ONLY: IP_UKCA_MODE_NUCLEATION, IP_UKCA_MODE_AITKEN, IP_UKCA_MODE_ACCUM, &
    & IP_UKCA_MODE_COARSE
USE UKCA_CONSTANTS, ONLY: ZBOLTZ
USE UKCA_RADAER_COMPTYPE, ONLY: IP_UKCA_SULPHATE, IP_UKCA_BLACKCARBON, &
 & IP_UKCA_ORGANICCARBON, IP_UKCA_SEASALT, IP_UKCA_DUST, &
 & IP_UKCA_SECONDORGANIC, IP_UKCA_NITRATE, IP_UKCA_AMMONIUM

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEVG, KDLEN, KLEN, KMXCOUNT
TYPE (TYPE_GOM_PLUS),INTENT(IN)  :: YDGP_TL
TYPE (TYPE_GOM_PLUS),INTENT(IN)  :: YDGP5
REAL(KIND=JPRB)   ,INTENT(OUT) :: PAEROD_PRIME(KDLEN,KMXCOUNT) ! AEROSOL OPT DEPTH PERTURBATION 

!*       0.5   LOCAL VARIABLES
!              ---------------
INTEGER(KIND=JPIM) :: IFLAG
INTEGER(KIND=JPIM) :: IWAVLINDEX
INTEGER(KIND=JPIM) :: JLEN, JLEV, JK, JL, JF, JJ
INTEGER(KIND=JPIM) :: JMODE, JCP !Glomap - mode identifier and component identifier
INTEGER(KIND=JPIM), DIMENSION(NMODES) :: IMTYPE !Glomap: type of each mode from ukca_radaer_modetype.F90.
INTEGER(KIND=JPIM), DIMENSION(NMODES,NCP) :: ICPTYPE !Glomap: type of each component in each mode from ukca_radaer_comptype.F90.

!Variables that describe the mode
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES) :: ZNPM5 !Glomap: number of particles per mole of air in each mode 
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES,NCP) :: ZCPMMR5 !Glomap: mass mixing ratio of each component in each mode 
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES) :: ZMODAL_MMR5 !Glomap: mass mixing ratio of each mode 
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES) :: ZMODAL_MMRTL !Glomap: mass mixing ratio of each mode 
REAL(KIND=JPRB) :: ZMODAL_NBR5 !Glomap: particles per m3 
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG)   :: ZTMMR5 ! total aer mixing ratios nonlinear
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES):: ZFRACMMR5 ! fraction of mixing ratios wrt total mass
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG)   :: ZTMMRTL ! total aer mixing ratios pert

!Variables for calcuting the AOD
REAL(KIND=JPRB), DIMENSION(KLEN,KMXCOUNT) :: ZAEROD5
REAL(KIND=JPRB), DIMENSION(KLEN,NMODES)   :: ZAERTAU5
REAL(KIND=JPRB), DIMENSION(KLEN,NMODES)   :: ZAERTAU_TL
REAL(KIND=JPRB)    :: ZAEREXT5, ZAEREXT_TL, ZDP5, ZDP_TL

!Meteorologoical variables
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG)   :: ZQSAT5 
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG)   :: ZT5 !Temperature at observation points, model levels (INPUT through YDGP5)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG)   :: ZPF5 !Full level pressure values at obs points, model levels (INPUT through YDGP5)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG)   :: ZQF5 !Specific humidity at observation points, model levels (INPUT through YDGP5)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG) :: ZRH5 !Glomap: relative humidity on a scale of 0-1.
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG) :: ZAIRD !Glomap - number density of air per cm3

!Variables needed for coefficient calculation
REAL(KIND=JPRB) :: ZLOC_ABS5, ZLOC_SCA5 !Glomap: Absorption and scattering coefficients for a single mode and one model level
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES) :: ZWET_DIAM5 !Glomap: Avg wet diameter of size mode (m)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES) :: ZDRY_DIAM5 !Glomap: Median particle dry diameter for each mode (m)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES,NCP) :: ZCOMP_VOL5 !Glomap: Partial volumes of each component in each mode (m3)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES) :: ZMODAL_WTV5 !Glomap: Partial volume of water in each mode (m3)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES) :: ZMODAL_RHO5 !Glomap: Particle density [incl. H2O & insoluble cpts] (kgm^-3)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES) :: ZMODAL_VOL5 !Glomap: Avg volume of size mode (m3)

REAL(KIND=JPRB), PARAMETER :: ZTHRESHOLD_MMR = 0.0_JPRB, & ! kg/kg
   &                 ZTHRESHOLD_NBR = 1.0E+04_JPRB    ! m-3

LOGICAL, DIMENSION(NMODES) :: LLMSOL !Glomap: soluble or not.

REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL :: LLPHYLIN 

#include "satur.intfb.h"
#include "aod_glomap_convert_mdnd.intfb.h"
#include "aod_glomap_calc_coeff.intfb.h"
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AOD_GLOMAP_TL',0,ZHOOK_HANDLE)

!Start by copying the number density (particles per molecule of air) of each mode and 
!the mass mixing ratio of each component into local arrays
ZNPM5(:,:,:)=0._JPRB !(needed since no secondary organics?)
ZCPMMR5(:,:,:,:)=0._JPRB !(needed since no secondary organics?)
ZMODAL_MMRTL(:,:,:)=0._JPRB
ICPTYPE(:,:)=0_JPIM 

DO JLEN=1,KLEN
!    !*****************How are the GFL fields setup?**********************************
!    !Check that this structure still matches ifs/module/ukca_mode_setup.F90
    DO JF=1,YDGP5%NGEMS
        !MODES
        IF(YDGP5%GEMS_IGRIB(JF) == 212020) THEN !Nucleation soluble
            ZNPM5(JLEN,:,1) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            IMTYPE(1) = IP_UKCA_MODE_NUCLEATION
            LLMSOL(1) = .TRUE.
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212021) THEN !Aitken soluble
            ZNPM5(JLEN,:,2) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            IMTYPE(2) = IP_UKCA_MODE_AITKEN
            LLMSOL(2) = .TRUE.
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212022) THEN !Accumulation Soluble
            ZNPM5(JLEN,:,3) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            IMTYPE(3) = IP_UKCA_MODE_ACCUM
            LLMSOL(3) = .TRUE.
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212023) THEN !Coarse Soluble
            ZNPM5(JLEN,:,4) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            IMTYPE(4) = IP_UKCA_MODE_COARSE
            LLMSOL(4) = .TRUE.
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212024) THEN !Aitken Insoluble
            ZNPM5(JLEN,:,5) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            IMTYPE(5) = IP_UKCA_MODE_AITKEN
            LLMSOL(5) = .FALSE.
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212025) THEN !Accumulation Insoluble
            ZNPM5(JLEN,:,6) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            IMTYPE(6) = IP_UKCA_MODE_ACCUM
            LLMSOL(6) = .FALSE.
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212026) THEN !Coarse Insoluble
            ZNPM5(JLEN,:,7) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            IMTYPE(7) = IP_UKCA_MODE_COARSE
            LLMSOL(7) = .FALSE.
        ENDIF
        
     
        !COMPONENTS
        IF(YDGP5%GEMS_IGRIB(JF) == 212001) THEN  !Nucleation soluble SO4
            ZCPMMR5(JLEN,:,1,1) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(1,1) = IP_UKCA_SULPHATE
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212019) THEN  !Nucleation soluble OM
            ZCPMMR5(JLEN,:,1,3) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(1,3) = IP_UKCA_ORGANICCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212002) THEN  !Aitken soluble SO4
            ZCPMMR5(JLEN,:,2,1) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(2,1) = IP_UKCA_SULPHATE
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212005) THEN  !Aitken soluble BC
            ZCPMMR5(JLEN,:,2,2) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(2,2) = IP_UKCA_BLACKCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212009) THEN  !Aitken soluble OM
            ZCPMMR5(JLEN,:,2,3) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(2,3) = IP_UKCA_ORGANICCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212003) THEN  !Accumulation Soluble SO4
            ZCPMMR5(JLEN,:,3,1) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(3,1) = IP_UKCA_SULPHATE
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212006) THEN  !Accumulation Soluble BC
            ZCPMMR5(JLEN,:,3,2) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(3,2) = IP_UKCA_BLACKCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212010) THEN  !Accumulation Soluble OM
            ZCPMMR5(JLEN,:,3,3) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(3,3) = IP_UKCA_ORGANICCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212013) THEN  !Accumulation Soluble SS
            ZCPMMR5(JLEN,:,3,4) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(3,4) = IP_UKCA_SEASALT
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212015) THEN  !Accumulation Soluble DU
            ZCPMMR5(JLEN,:,3,5) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(3,5) = IP_UKCA_SEASALT
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212004) THEN  !Coarse Soluble SO4
            ZCPMMR5(JLEN,:,4,1) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(4,1) = IP_UKCA_SULPHATE
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212007) THEN  !Coarse Soluble BC
            ZCPMMR5(JLEN,:,4,2) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(4,2) = IP_UKCA_BLACKCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212011) THEN  !Coarse Soluble OM
            ZCPMMR5(JLEN,:,4,3) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(4,3) = IP_UKCA_ORGANICCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212014) THEN  !Coarse Soluble SS
            ZCPMMR5(JLEN,:,4,4) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(4,4) = IP_UKCA_SEASALT
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212016) THEN  !Coarse Soluble DU
            ZCPMMR5(JLEN,:,4,5) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(4,5) = IP_UKCA_DUST
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212008) THEN  !Aitken Insoluble BC
            ZCPMMR5(JLEN,:,5,2) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(5,2) = IP_UKCA_BLACKCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212012) THEN  !Aitken Insoluble OM
            ZCPMMR5(JLEN,:,5,3) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(5,3) = IP_UKCA_ORGANICCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212017) THEN  !Accumulation Insoluble DU
            ZCPMMR5(JLEN,:,6,5) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(6,5) = IP_UKCA_DUST
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212018) THEN  !Coarse Insoluble DU
            ZCPMMR5(JLEN,:,7,5) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(7,5) = IP_UKCA_DUST
        ENDIF
        !!WARNING: NITRATE AND AMMONIUM not included
        
        IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERLG) THEN ! Total aerosol mixing ratio and pert
           ZTMMRTL(JLEN,:)=YDGP_TL%GEMS(JLEN,1:KFLEVG,IH,JF)
        ENDIF
    ENDDO
    
    !Gather the component mass mixing ratios together to give the modal mass mixing ratios
    ZMODAL_MMR5(JLEN,:,:) = 0.0_JPRB
    DO JMODE=1,NMODES
      IF(MODE(JMODE)) THEN
         DO JCP=1,NCP
            IF(COMPONENT(JMODE,JCP)) THEN
              ZMODAL_MMR5(JLEN,:,JMODE) = ZMODAL_MMR5(JLEN,:,JMODE) + ZCPMMR5(JLEN,:,JMODE,JCP)
            ENDIF
         ENDDO
      ENDIF
    ENDDO

ENDDO !End observations loop

! Compute total mixing ratio by summing single contributions from all but the nucleation mode
ZTMMR5(:,:) = 0.0_JPRB
DO JMODE=2,NMODES
  IF(MODE(JMODE)) THEN
      ZTMMR5(:,:) = ZTMMR5(:,:) + ZMODAL_MMR5(:,:,JMODE)
  ENDIF
ENDDO

! The perturbations on the various species are computed from the perturbation
! on the total mass mixing ratio using the mass fraction assumed constant 
ZFRACMMR5(:,:,:) = 0.0_JPRB
DO JLEN=1,KLEN
  DO JMODE=2,NMODES
    IF(MODE(JMODE)) THEN
      DO JLEV=1,KFLEVG
        IF(ZTMMR5(JLEN,JLEV) /= 0.0_JPRB) THEN
          ZFRACMMR5(JLEN,JLEV,JMODE) = ZMODAL_MMR5(JLEN,JLEV,JMODE)/ZTMMR5(JLEN,JLEV)
        ENDIF
        ZMODAL_MMRTL(JLEN,JLEV,JMODE)=ZTMMRTL(JLEN,JLEV) * ZFRACMMR5(JLEN,JLEV,JMODE)
      ENDDO
    ENDIF
  ENDDO
ENDDO

!Meteorological variables needed for calculation of extinction coefficient
ZT5(1:KLEN,:) = YDGP5%TF(1:KLEN,1:KFLEVG,IH) !Temperature at obs locations
ZPF5(1:KLEN,:)= YDGP5%PRESF(1:KLEN,:,IH) !Pressure at obs locations
ZQF5(1:KLEN,:) = YDGP5%QF(1:KLEN,1:KFLEVG,IH) !Specific humidity at obs locations
CALL SATUR (YRTHF, YRCST, 1 , KLEN , KLEN  , 1 , KFLEVG, LLPHYLIN, &
  & ZPF5, ZT5 , ZQSAT5 , IFLAG)
ZRH5(1:KLEN,:) = ZQF5(1:KLEN,:)/ZQSAT5(1:KLEN,1:KFLEVG) !Compute relative humidity at observation points, on model levels
ZAIRD(1:KLEN,:) = ZPF5(1:KLEN,:)/(ZT5(1:KLEN,:)*ZBOLTZ*1.0E6_JPRB) !number denisty of air per cm3 

!Aerosol mode/component information directly needed for calculation of extinction coefficient
! initialise to zero
ZWET_DIAM5(:,:,:)=0._JPRB
ZDRY_DIAM5(:,:,:) = 0._JPRB
ZCOMP_VOL5(:,:,:,:) = 0._JPRB
ZMODAL_WTV5(:,:,:) = 0._JPRB
ZMODAL_RHO5(:,:,:) = 0._JPRB
ZMODAL_VOL5(:,:,:) = 0._JPRB

!Convert the mass mixing ratio and number concentrations into the quantaties needed
!calculate the extinction (absorption and scattering) coefficient.
CALL AOD_GLOMAP_CONVERT_MDND(KLEN, KFLEVG, ZNPM5, ZCPMMR5, LLMSOL, ZT5, ZPF5, ZQF5, ZRH5, ZAIRD, &
   & ZWET_DIAM5, ZDRY_DIAM5, ZCOMP_VOL5, ZMODAL_WTV5, ZMODAL_RHO5, ZMODAL_VOL5)
!If an individual component is less than zero then the ZCOMP_VOL is set to zero and 
!the component won't count towards the extinction coefficient. If all components in
!a mode are zero then the ZMODAL_MMR will be less than the ZTHRESHOLD_MMR and the
!extinction coefficient will not be calculated and remain zero below.

! initialise perturbations to zero 
PAEROD_PRIME(:,:) = 0.0_JPRB 
! initialise trajectory 
ZAEROD5(:,:) = 0.0_JPRB 

!Wavelengths available in GLOMAP are currently
!3.80000E-07, 4.40000E-07, 5.50000E-07, 6.70000E-07, 8.70000E-07, 1.02000E-06
!Taken from pcalc_hadgem_v3.ukca (search 'aerosol optical depth')
IF(KMXCOUNT ==1) IWAVLINDEX=3   ! 550nm
DO JJ=1,KMXCOUNT !Wavelengths
!
! *** to be activated, for now assume AOD is only defined at 550nm , i.e. IWAVLINDEX=3 always
!     IWAVLINDEX=IWAVL(JJ)

ZAERTAU5(:,:)=0.0_JPRB
ZAERTAU_TL(:,:)=0.0_JPRB

  DO JL=1,KLEN !Observation
    DO JMODE=2,NMODES !Modes - exclude nucleation mode
 
       DO JK=1,KFLEVG !Model level
!    
         IF(MODE(JMODE)) THEN
         
            ZMODAL_NBR5 = ZNPM5(JL,JK,JMODE)*ZAIRD(JL,JK) * 1E6 ! m-3

            !Calculate the scattering and absorption coefficients
            ZLOC_ABS5 = 0.0_JPRB
            ZLOC_SCA5 = 0.0_JPRB
            IF (ZMODAL_MMR5(JL,JK,JMODE) > ZTHRESHOLD_MMR .AND. ZMODAL_NBR5 > ZTHRESHOLD_NBR) THEN
               CALL AOD_GLOMAP_CALC_COEFF(IWAVLINDEX, JMODE, IMTYPE(JMODE), ICPTYPE(JMODE,1:NCP), LLMSOL(JMODE), &
                 &  ZWET_DIAM5(JL,JK,JMODE), ZDRY_DIAM5(JL,JK,JMODE), ZCOMP_VOL5(JL,JK,JMODE,1:NCP), ZMODAL_WTV5(JL,JK,JMODE), &
                 &  ZMODAL_VOL5(JL,JK,JMODE), ZMODAL_RHO5(JL,JK,JMODE), ZLOC_ABS5, ZLOC_SCA5)
            ENDIF
 
            !Calculate mass mixing ratio of the mode times the extinction coefficient
            ZAEREXT5=ZMODAL_MMR5(JL,JK,JMODE)*(ZLOC_ABS5+ZLOC_SCA5)
            ZAEREXT_TL=ZMODAL_MMRTL(JL,JK,JMODE)*(ZLOC_ABS5+ZLOC_SCA5)
             
            !Calculate the thickness
            ZDP5=YDGP5%PRESH(JL,JK,IH)-YDGP5%PRESH(JL,JK-1,IH)
            ZDP_TL=YDGP_TL%PRESH(JL,JK,IH)-YDGP_TL%PRESH(JL,JK-1,IH)
             
            !Calculate the AOD for the mode across all the model levels
            ZAERTAU5(JL,JMODE)= ZAERTAU5(JL,JMODE) + ZAEREXT5*(ZDP5/RG)
            ZAERTAU_TL(JL,JMODE) = ZAERTAU_TL(JL,JMODE) + ZAEREXT_TL*(ZDP5/RG) + ZAEREXT5*(ZDP_TL/RG)
             
         ENDIF !end test of mode
 
       ENDDO ! loop over model vertical levels
      
      !Sum over MODES to get total optical depth as function of horizontal position and channel
      ZAEROD5(JL,JJ) = ZAEROD5(JL,JJ) + ZAERTAU5(JL,JMODE)
      PAEROD_PRIME(JL,JJ) = PAEROD_PRIME(JL,JJ) + ZAERTAU_TL(JL,JMODE)
      
    ENDDO ! loop over aerosol modes
  ENDDO ! loop over horizontal observational points
ENDDO ! loop over MODIS channels

! Security check for negative optical depths 
DO JJ = 1,KMXCOUNT
  DO JL=1, KLEN
    IF(ZAEROD5(JL,JJ) < 0.0_JPRB) THEN
      PAEROD_PRIME(JL,JJ)= 0.0_JPRB
    ENDIF
  ENDDO
ENDDO

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AOD_GLOMAP_TL',1,ZHOOK_HANDLE)
END SUBROUTINE AOD_GLOMAP_TL
