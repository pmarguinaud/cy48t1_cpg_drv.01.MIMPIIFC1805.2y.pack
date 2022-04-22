SUBROUTINE AOD_GLOMAP_OP&
 & (KFLEVG, KDLEN  , KLEN  , KMXCOUNT, &
 & YDGP5,&
 & PAEROD )

!*** *AOD_GLOMAP_OP* Operator routine for aerosol optical depth observations with the GLOMAP model.
!
!**   INTERFACE.
!     ----------
!          *AOD_GLOMAP_OP* IS CALLED FROM *HOP*|*OBSOP_COMPOSITION* (oops design).
!
!     WHERE KFLEVG    = Number of vertical levels
!           KDLEN     = First dimension: Max number of observations, this set (INPUT)
!           KLEN      = Number of aerosol optical depth obs (INPUT)
!           KMXCOUNT  = Number of channels (INPUT)
!           YDGP5     = Full nonlinear trajectory fields at observation locations (INPUT)
!           PAEROD    = Aerosol optical depth at different wavelengths (OUTPUT)
!
!**   AUTHOR.
!     -------
!        Melanie Ades        *ECMWF*
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2018-08-01 (adapted from aod_op.F90 by A. Benedetti and ukca_radaer_ifs.F90/ukca_radaer_compute_aod.F90)
!
!
!     PURPOSE.
!     --------
!     - derive aerosol optical thickness from profiles of mode number density and component mass mixing ratio
!-----------------------------------------------------------------------

USE PARKIND1       , ONLY : JPIM, JPRB
USE YOMHOOK        , ONLY : LHOOK, DR_HOOK
USE YOMCST         , ONLY : RG, YRCST
USE YOETHF         , ONLY : YRTHF
USE YOMLUN         , ONLY : NULOUT
USE GOM_PLUS       , ONLY : TYPE_GOM_PLUS, IH
USE YOEPHLI        , ONLY : YREPHLI
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
TYPE (TYPE_GOM_PLUS),INTENT(IN)  :: YDGP5
REAL(KIND=JPRB)   ,INTENT(OUT) :: PAEROD(KDLEN,KMXCOUNT) ! AEROSOL OPT DEPTH (FUNCTION OF HORIZONTAL POSITION AND CHANNEL #)

!*       0.5   LOCAL VARIABLES
!              ---------------
INTEGER(KIND=JPIM) :: IFLAG
INTEGER(KIND=JPIM) :: IWAVLINDEX
INTEGER(KIND=JPIM) :: JLEN, JK, JL, JF, JJ
INTEGER(KIND=JPIM), DIMENSION(NMODES) :: IMTYPE !Glomap: type of each mode from ukca_radaer_modetype.F90.
INTEGER(KIND=JPIM), DIMENSION(NMODES,NCP) :: ICPTYPE !Glomap: type of each component in each mode from ukca_radaer_comptype.F90.

REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG)   :: ZQSAT 
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG)   :: ZT !Temperature at observation points, model levels (INPUT through YDGP5)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG)   :: ZPF !Full level pressure values at obs points, model levels (INPUT through YDGP5)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG)   :: ZQF !Specific humidity at observation points, model levels (INPUT through YDGP5)
REAL(KIND=JPRB), DIMENSION(KLEN,NMODES)   :: ZAERTAU
REAL(KIND=JPRB)    :: ZAEREXT, ZDP

REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES) :: ZNPM !Glomap: number of particles per mole of air in each mode 
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES,NCP) :: ZCPMMR !Glomap: mass mixing ratio of each component in each mode 
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES) :: ZMODAL_MMR !Glomap: mass mixing ratio of each mode 
REAL(KIND=JPRB) :: ZMODAL_NBR !Glomap: particles per m3 
REAL(KIND=JPRB) :: ZLOC_ABS, ZLOC_SCA !Glomap: Absorption and scattering coefficients for a single mode and one model level
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG) :: ZAIRD !Glomap - number density of air per cm3
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG) :: ZRH !Glomap: relative humidity on a scale of 0-1.
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES) :: ZWET_DIAM !Glomap: Avg wet diameter of size mode (m)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES) :: ZDRY_DIAM !Glomap: Median particle dry diameter for each mode (m)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES,NCP) :: ZCOMP_VOL !Glomap: Partial volumes of each component in each mode (m3)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES) :: ZMODAL_WTV !Glomap: Partial volume of water in each mode (m3)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES) :: ZMODAL_RHO !Glomap: Particle density [incl. H2O & insoluble cpts] (kgm^-3)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG,NMODES) :: ZMODAL_VOL !Glomap: Avg volume of size mode (m3)

REAL(KIND=JPRB), PARAMETER :: ZTHRESHOLD_MMR = 0.0_JPRB, & ! kg/kg
   &                 ZTHRESHOLD_NBR = 1.0E+04_JPRB    ! m-3

LOGICAL, DIMENSION(NMODES) :: LLMSOL !Glomap: soluble or not.

INTEGER(KIND=JPIM) :: JMODE, JCP !Glomap - mode identifier and component identifier

REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL :: LLPHYLIN 

#include "satur.intfb.h"
#include "aod_glomap_convert_mdnd.intfb.h"
#include "aod_glomap_calc_coeff.intfb.h"
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AOD_GLOMAP_OP',0,ZHOOK_HANDLE)

LLPHYLIN = YREPHLI%LPHYLIN

WRITE(NULOUT,*) 'MA: Made it into Glomap!'

!Start by copying the number density (particles per molecule of air) of each mode and 
!the mass mixing ratio of each component into local arrays
ZNPM(:,:,:)=0._JPRB !(needed since no secondary organics?)
ZCPMMR(:,:,:,:)=0._JPRB !(needed since no secondary organics?)
ICPTYPE(:,:)=0_JPIM 
DO JLEN=1,KLEN
    !*****************How are the GFL fields setup?**********************************
    !Check that this structure still matches ifs/module/ukca_mode_setup.F90
    DO JF=1,YDGP5%NGEMS
        !MODES
        IF(YDGP5%GEMS_IGRIB(JF) == 212020) THEN !Nucleation soluble
            ZNPM(JLEN,:,1) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            IMTYPE(1) = IP_UKCA_MODE_NUCLEATION
            LLMSOL(1) = .TRUE.
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212021) THEN !Aitken soluble
            ZNPM(JLEN,:,2) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            IMTYPE(2) = IP_UKCA_MODE_AITKEN
            LLMSOL(2) = .TRUE.
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212022) THEN !Accumulation Soluble
            ZNPM(JLEN,:,3) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            IMTYPE(3) = IP_UKCA_MODE_ACCUM
            LLMSOL(3) = .TRUE.
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212023) THEN !Coarse Soluble
            ZNPM(JLEN,:,4) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            IMTYPE(4) = IP_UKCA_MODE_COARSE
            LLMSOL(4) = .TRUE.
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212024) THEN !Aitken Insoluble
            ZNPM(JLEN,:,5) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            IMTYPE(5) = IP_UKCA_MODE_AITKEN
            LLMSOL(5) = .FALSE.
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212025) THEN !Accumulation Insoluble
            ZNPM(JLEN,:,6) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            IMTYPE(6) = IP_UKCA_MODE_ACCUM
            LLMSOL(6) = .FALSE.
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212026) THEN !Coarse Insoluble
            ZNPM(JLEN,:,7) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            IMTYPE(7) = IP_UKCA_MODE_COARSE
            LLMSOL(7) = .FALSE.
        ENDIF
        
     
        !COMPONENTS
        IF(YDGP5%GEMS_IGRIB(JF) == 212001) THEN  !Nucleation soluble SO4
            ZCPMMR(JLEN,:,1,1) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(1,1) = IP_UKCA_SULPHATE
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212019) THEN  !Nucleation soluble OM
            ZCPMMR(JLEN,:,1,3) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(1,3) = IP_UKCA_ORGANICCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212002) THEN  !Aitken soluble SO4
            ZCPMMR(JLEN,:,2,1) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(2,1) = IP_UKCA_SULPHATE
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212005) THEN  !Aitken soluble BC
            ZCPMMR(JLEN,:,2,2) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(2,2) = IP_UKCA_BLACKCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212009) THEN  !Aitken soluble OM
            ZCPMMR(JLEN,:,2,3) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(2,3) = IP_UKCA_ORGANICCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212003) THEN  !Accumulation Soluble SO4
            ZCPMMR(JLEN,:,3,1) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(3,1) = IP_UKCA_SULPHATE
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212006) THEN  !Accumulation Soluble BC
            ZCPMMR(JLEN,:,3,2) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(3,2) = IP_UKCA_BLACKCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212010) THEN  !Accumulation Soluble OM
            ZCPMMR(JLEN,:,3,3) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(3,3) = IP_UKCA_ORGANICCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212013) THEN  !Accumulation Soluble SS
            ZCPMMR(JLEN,:,3,4) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(3,4) = IP_UKCA_SEASALT
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212015) THEN  !Accumulation Soluble DU
            ZCPMMR(JLEN,:,3,5) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(3,5) = IP_UKCA_SEASALT
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212004) THEN  !Coarse Soluble SO4
            ZCPMMR(JLEN,:,4,1) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(4,1) = IP_UKCA_SULPHATE
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212007) THEN  !Coarse Soluble BC
            ZCPMMR(JLEN,:,4,2) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(4,2) = IP_UKCA_BLACKCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212011) THEN  !Coarse Soluble OM
            ZCPMMR(JLEN,:,4,3) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(4,3) = IP_UKCA_ORGANICCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212014) THEN  !Coarse Soluble SS
            ZCPMMR(JLEN,:,4,4) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(4,4) = IP_UKCA_SEASALT
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212016) THEN  !Coarse Soluble DU
            ZCPMMR(JLEN,:,4,5) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(4,5) = IP_UKCA_DUST
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212008) THEN  !Aitken Insoluble BC
            ZCPMMR(JLEN,:,5,2) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(5,2) = IP_UKCA_BLACKCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212012) THEN  !Aitken Insoluble OM
            ZCPMMR(JLEN,:,5,3) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(5,3) = IP_UKCA_ORGANICCARBON
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212017) THEN  !Accumulation Insoluble DU
            ZCPMMR(JLEN,:,6,5) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(6,5) = IP_UKCA_DUST
        ENDIF
        IF(YDGP5%GEMS_IGRIB(JF) == 212018) THEN  !Coarse Insoluble DU
            ZCPMMR(JLEN,:,7,5) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
            ICPTYPE(7,5) = IP_UKCA_DUST
        ENDIF
        !!WARNING: NITRATE AND AMMONIUM not included
    ENDDO
    
    !Gather the component mass mixing ratios together to give the modal mass mixing ratios
    ZMODAL_MMR(JLEN,:,:) = 0.0_JPRB
    DO JMODE=1,NMODES
      IF(MODE(JMODE)) THEN
         DO JCP=1,NCP
            IF(COMPONENT(JMODE,JCP)) THEN
              ZMODAL_MMR(JLEN,:,JMODE) = ZMODAL_MMR(JLEN,:,JMODE) + ZCPMMR(JLEN,:,JMODE,JCP)
            ENDIF
         ENDDO
      ENDIF
    ENDDO

ENDDO !End observations loop

!Meteorological variables needed for calculation of extinction coefficient
ZT(1:KLEN,:) = YDGP5%TF(1:KLEN,1:KFLEVG,IH) !Temperature at obs locations
ZPF(1:KLEN,:)= YDGP5%PRESF(1:KLEN,:,IH) !Pressure at obs locations
ZQF(1:KLEN,:) = YDGP5%QF(1:KLEN,1:KFLEVG,IH) !Specific humidity at obs locations
CALL SATUR (YRTHF, YRCST, 1 , KLEN , KLEN  , 1 , KFLEVG, LLPHYLIN, &
  & ZPF, ZT , ZQSAT , IFLAG)
ZRH(1:KLEN,:) = ZQF(1:KLEN,:)/ZQSAT(1:KLEN,1:KFLEVG) !Compute relative humidity at observation points, on model levels
ZAIRD(1:KLEN,:) = ZPF(1:KLEN,:)/(ZT(1:KLEN,:)*ZBOLTZ*1.0E6_JPRB) !number denisty of air per cm3 

!Aerosol mode/component information directly needed for calculation of extinction coefficient
! initialise to zero
ZWET_DIAM(:,:,:)=0._JPRB
ZDRY_DIAM(:,:,:) = 0._JPRB
ZCOMP_VOL(:,:,:,:) = 0._JPRB
ZMODAL_WTV(:,:,:) = 0._JPRB
ZMODAL_RHO(:,:,:) = 0._JPRB
ZMODAL_VOL(:,:,:) = 0._JPRB

!Convert the mass mixing ratio and number concentrations into the quantaties needed
!calculate the extinction (absorption and scattering) coefficient.
CALL AOD_GLOMAP_CONVERT_MDND(KLEN, KFLEVG, ZNPM, ZCPMMR, LLMSOL, ZT, ZPF, ZQF, ZRH, ZAIRD, &
   & ZWET_DIAM, ZDRY_DIAM, ZCOMP_VOL, ZMODAL_WTV, ZMODAL_RHO, ZMODAL_VOL)
!If an individual component is less than zero then the ZCOMP_VOL is set to zero and 
!the component won't count towards the extinction coefficient. If all components in
!a mode are zero then the ZMODAL_MMR will be less than the ZTHRESHOLD_MMR and the
!extinction coefficient will not be calculated and remain zero below.

! initialise optical depth 
PAEROD(:,:) = 0.0_JPRB 

!Wavelengths available in GLOMAP are currently
!3.80000E-07, 4.40000E-07, 5.50000E-07, 6.70000E-07, 8.70000E-07, 1.02000E-06
!Taken from pcalc_hadgem_v3.ukca (search 'aerosol optical depth')
IF(KMXCOUNT ==1) IWAVLINDEX=3   ! 550nm
DO JJ=1,KMXCOUNT !Wavelengths
!
! *** to be activated, for now assume AOD is only defined at 550nm , i.e. IWAVLINDEX=3 always
!     IWAVLINDEX=IWAVL(JJ)
ZAERTAU(:,:)=0.0_JPRB

  DO JL=1,KLEN !Observation
    DO JMODE=2,NMODES !Modes - exclude nucleation mode
      DO JK=1,KFLEVG !Model level
   
        IF(MODE(JMODE)) THEN

            ZMODAL_NBR = ZNPM(JL,JK,JMODE)*ZAIRD(JL,JK) * 1E6 ! m-3

            !Calculate the scattering and absorption coefficients
            ZLOC_ABS = 0.0_JPRB
            ZLOC_SCA = 0.0_JPRB
            IF (ZMODAL_MMR(JL,JK,JMODE) > ZTHRESHOLD_MMR .AND. ZMODAL_NBR > ZTHRESHOLD_NBR) THEN
               CALL AOD_GLOMAP_CALC_COEFF(IWAVLINDEX, JMODE, IMTYPE(JMODE), ICPTYPE(JMODE,1:NCP), LLMSOL(JMODE), &
                 &  ZWET_DIAM(JL,JK,JMODE), ZDRY_DIAM(JL,JK,JMODE), ZCOMP_VOL(JL,JK,JMODE,1:NCP), ZMODAL_WTV(JL,JK,JMODE), &
                 &  ZMODAL_VOL(JL,JK,JMODE), ZMODAL_RHO(JL,JK,JMODE), ZLOC_ABS, ZLOC_SCA)
            ENDIF

            !Calculate mass mixing ratio of the mode times the extinction coefficient
            ZAEREXT=ZMODAL_MMR(JL,JK,JMODE)*(ZLOC_ABS+ZLOC_SCA)
            
            !Calculate the thickness
            ZDP=YDGP5%PRESH(JL,JK,IH)-YDGP5%PRESH(JL,JK-1,IH)
            
            !Calculate the AOD for the mode across all the model levels
            ZAERTAU(JL,JMODE)=ZAERTAU(JL,JMODE)+ ZAEREXT*(ZDP/RG)
        ENDIF

      ENDDO ! loop over model vertical levels
      
      !Sum over MODES to get total optical depth as function of horizontal position and channel
      PAEROD(JL,JJ)=PAEROD(JL,JJ)+ZAERTAU(JL,JMODE)
    ENDDO ! loop over aerosol modes
  ENDDO ! loop over horizontal observational points
ENDDO ! loop over MODIS channels

! Security check for negative optical depths
DO  JJ = 1,KMXCOUNT
  DO JL=1, KLEN
    IF(PAEROD(JL,JJ) < 0.0_JPRB) THEN
       PAEROD(JL,JJ)= 0.0_JPRB
    ENDIF
  ENDDO
ENDDO


!-----------------------------------------------------------------------
!END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AOD_GLOMAP_OP',1,ZHOOK_HANDLE)
END SUBROUTINE AOD_GLOMAP_OP

