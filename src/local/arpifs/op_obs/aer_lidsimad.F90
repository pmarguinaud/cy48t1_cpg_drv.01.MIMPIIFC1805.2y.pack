SUBROUTINE AER_LIDSIMAD &
  &( KDLEN, KLEN , KLEV , KDFS, LDPHYLIN,  &
  &  PLIDPRES , & 
  &  PCO215 , PNO215, PGEOP15, PTQO15  ,PAELIDPROF5,  &
  &  PCO21  , PNO21 , PGEOP1 , PTQO1   ,PAELIDPROF,  &
  &  YDGP5, &
  &  PLISIS, PLISIT, PULISIS,PULISIT,PEXTT, PEXTS &
  & )

!**** *AER_LIDSIMAD* - ADJOINT VERSION: COMPUTES GRADIENTS FOR THE LIDAR SIGNAL OUT OF THE ATMOSPHERE
!                    ASSUMED CLOUDLESS BUT INCLUDING THE GEMS-AEROSOLS.

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE TOTAL BACK-SCATTERING COEFFICIENT 
!          ACCOUNTING RAYLEIGH SCATTERING, AEROSOL SCATTERING AND GASEOUS 
!          ABSORPTION FOR THE THREE USUAL WAVELENGTHS (355, 532 AND 1064 NM)
!          FROM THE STANDARD "LIDAR EQUATION".

!**   INTERFACE.
!     ----------

!          *AER_LIDSIMAD* IS CALLED FROM *HOPAD*

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
! PLIDPRES (KDLEN,KLEV) : LIDAR LEVEL PRESSURE (Pa)
! PCO215   (KDLEN,KLEV)   : CO2 CONCENTRATION (kg/kg) from traj
! PNO215   (KDLEN,KLEV)   : NO2 CONCENTRATION (kg/kg) from traj
! PTQO15   (KDLEN,KLEV)   : TEMPERATURE (K), SPECIFIC HUMIDITY (kg/kg) and OZONE from traj
! PGEOP15  (KDLEN,KLEV)   : GEOPOTENTIAL from traj
! PCO21    (KDLEN,KLEV)   : Gradient for CO2 CONCENTRATION (kg/kg)
! PNO21    (KDLEN,KLEV)   : Gradient for NO2 CONCENTRATION (kg/kg)
! PTQO1    (KDLEN,KLEV,KDFS)   : Gradient for TEMPERATURE (K), SPECIFIC HUMIDITY (kg/kg) and OZONE
! PGEOP1   (KDLEN,KLEV)   : Gradient for GEOPOTENTIAL

!     ==== OUTPUTS ===
! PLISIS (KDLEN,NWLID,0:KLEV) : Gradient for BACK-SCATTERING SIGNAL FOR SURFACE LIDAR
! PLISIT (KDLEN,NWLID,0:KLEV) : Gradient for BACK-SCATTERING SIGNAL FOR SATELLITE LIDAR
! PULISIS (KDLEN,NWLID,0:KLEV): Gradient for UNATTENUATED BACK-SCATTERING SIGNAL FOR SURFACE LIDAR
! PULISIT (KDLEN,NWLID,0:KLEV): Gradient for UNATTENUATED BACK-SCATTERING SIGNAL FOR SATELLITE LIDAR
! PEXTT (KDLEN,NWLID,0:KLEV) :  Gradient for EXTINCTION FOR SATELLITE LIDAR
! PEXTS (KDLEN,NWLID,0:KLEV) :  Gradient for EXTINCTION FOR SURFACE LIDAR

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 20090702 

!     Angela Benedetti: 6-Sept-2010 Modifications to call it from hopad.F90
!     A. Geer           29-Dec-2015 OOPS cleaning (GOM_PLUS added; YGFL and GEMSPROF removed)
!        S. Remy     16-Sep-2016    Add nitrate and ammonium
!        S. Remy     13-Nov-2017    Add SOA
!     A Benedetti    05-Mar-2018    Add unattenuated backscatter
!     ------------------------------------------------------------------


USE PARKIND1      , ONLY : JPIM, JPRB
USE YOMHOOK       , ONLY : LHOOK, DR_HOOK
USE YOEAERLID     , ONLY : YREAERLID
USE YOEAEROP      , ONLY : ALF_SU, OMG_SU, ALF_OM, OMG_OM, ALF_DD, OMG_DD, ALF_SS,OMG_SS, &
 &                         ALF_BC, OMG_BC,  RALI_BC,RALI_DD,RALI_OM,RALI_SU,RALI_SS, ALF_NI, OMG_NI, &
 &                         RALI_NI, ALF_AM, OMG_AM, RALI_AM, ALF_SOA, OMG_SOA,RALI_SOA
USE YOEAERATM     , ONLY : YREAERATM
USE YOEAERSNK     , ONLY : YREAERSNK
USE YOEAERVOL     , ONLY : YREAERVOL
USE YOMCST        , ONLY : RG, RMD, RMCO2, RMNO2, RMO3, RNAVO , RPI, YRCST
USE YOETHF        , ONLY : YRTHF
USE GOM_PLUS      , ONLY : TYPE_GOM_PLUS
USE YOM_GRIB_CODES, ONLY : NGRBAERMR01,NGRBAERMR02, NGRBAERMR03, NGRBAERMR04, NGRBAERMR05, &
 &                         NGRBAERMR06, NGRBAERMR07, NGRBAERMR08, NGRBAERMR09, NGRBAERMR10, NGRBAERMR11, NGRBAERMR12, &
 &                         NGRBAERMR13, NGRBAERMR14, NGRBAERMR15, NGRBAERMR16, NGRBAERMR17, NGRBAERMR18, &
 &                         NGRBAERMR19, NGRBAERMR20, NGRBAERMR21, NGRBAERLG

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KDLEN 
INTEGER(KIND=JPIM),INTENT(IN) :: KLEN   ! Horizontal points
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV   ! Lidar vertical levels
INTEGER(KIND=JPIM),INTENT(IN) :: KDFS   !Size of T q O3 array from interpolation
LOGICAL           ,INTENT(IN) :: LDPHYLIN

TYPE(TYPE_GOM_PLUS), INTENT(IN) :: YDGP5 ! Model variables at observation locations

REAL(KIND=JPRB)   ,INTENT(IN) :: PLIDPRES(KDLEN,KLEV)  ! Fixed lidar pressure levels
REAL(KIND=JPRB)   ,INTENT(IN) :: PCO215(KDLEN,KLEV) ! CO2 profile from traj
REAL(KIND=JPRB)   ,INTENT(INOUT):: PCO21(KDLEN,KLEV)  ! Gradients on CO2
REAL(KIND=JPRB)   ,INTENT(IN) :: PNO215(KDLEN,KLEV) ! NO2 from traj
REAL(KIND=JPRB)   ,INTENT(INOUT):: PNO21(KDLEN,KLEV)  ! Gradients on NO2
REAL(KIND=JPRB)   ,INTENT(IN) :: PGEOP15(KDLEN,KLEV) ! Geopotential
REAL(KIND=JPRB)   ,INTENT(INOUT):: PGEOP1(KDLEN,KLEV)  ! Gradients on geopotential
REAL(KIND=JPRB)   ,INTENT(IN) :: PTQO15(KDLEN ,KLEV,KDFS) ! Interpolated profiles of T, q, O3
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTQO1(KDLEN ,KLEV,KDFS) ! Gradients on T, q, O3
REAL(KIND=JPRB)   ,INTENT(IN) :: PAELIDPROF5(KDLEN,KLEV,YDGP5%NGEMS) !GEMS profiles including aerosols from traj
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAELIDPROF(KDLEN,KLEV,YDGP5%NGEMS) !Gradients on GEMS profiles including aerosols

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLISIS(KDLEN,YREAERLID%NWLID,0:KLEV), PLISIT(KDLEN,YREAERLID%NWLID,0:KLEV), &
     &                            PULISIS(KDLEN,YREAERLID%NWLID,0:KLEV), PULISIT(KDLEN,YREAERLID%NWLID,0:KLEV), &
     &                            PEXTS(KDLEN,YREAERLID%NWLID,0:KLEV), PEXTT(KDLEN,YREAERLID%NWLID,0:KLEV) ! Gradients for lidar variables

!     ------------------------------------------------------------------

!*       0.1   LOCAL ARRAYS
!              ------------

INTEGER(KIND=JPIM) :: IAER, IBIN, IK, ITYP, IWAVL,IEFRH
INTEGER(KIND=JPIM) :: JAER, JL, JK, JTAB, JWL, JLEV, JP, JF
INTEGER(KIND=JPIM) :: IFLAG, IRH(KDLEN,KLEV)
INTEGER(KIND=JPIM), PARAMETER  :: JTABMX=12
INTEGER(KIND=JPIM) :: INACTAERO

LOGICAL :: LLPRINT

REAL(KIND=JPRB) :: ZAERF15(KDLEN,KLEV,YDGP5%NGEMS), ZAERF1(KDLEN,KLEV,YDGP5%NGEMS)

REAL(KIND=JPRB), DIMENSION(KLEN,KLEV) :: ZAERFT15 ! sum of single contributions of aer mixing ratio
REAL(KIND=JPRB), DIMENSION(KLEN,KLEV) :: ZAERTMR15 ! total aer mixing ratios from traj
REAL(KIND=JPRB), DIMENSION(KLEN,KLEV) :: ZAERTMR1 ! total aer mixing ratios grad 
REAL(KIND=JPRB), DIMENSION(KLEN,KLEV,YDGP5%NGEMS) :: ZFRACAER15 ! fraction of mixing ratios wrt total mass
REAL(KIND=JPRB), DIMENSION(KLEN,KLEV)   :: ZQSAT15, ZTF15, ZPF15 ! Temperature, pressure, q sat from traj
REAL(KIND=JPRB), DIMENSION(KLEN,KLEV)   :: ZTF1 , ZPF1 ! Gradients for temperature, pressure, q sat
REAL(KIND=JPRB), DIMENSION(KDLEN,0:KLEV):: ZGEOMH15, ZGEOMH1
REAL(KIND=JPRB), DIMENSION(KLEN,KLEV)   :: ZO3F15, ZQF15 ! O3 and q from traj
REAL(KIND=JPRB), DIMENSION(KLEN,KLEV)   :: ZO3F1, ZQF1 ! Gradients of O3 and q from traj
REAL(KIND=JPRB) :: ZDP15(KDLEN, KLEV), ZDZ15(KDLEN,KLEV)
REAL(KIND=JPRB) :: ZDP1(KDLEN, KLEV), ZDZ1(KDLEN,KLEV)
REAL(KIND=JPRB) :: ZRHCL15

REAL(KIND=JPRB) :: ZLISIS5(KDLEN,YREAERLID%NWLID,0:KLEV), ZLISIT5(KDLEN,YREAERLID%NWLID,0:KLEV) 

REAL(KIND=JPRB) :: ZOTRAY15(KDLEN,YREAERLID%NWLID,1:KLEV), ZOTCO215(KDLEN,YREAERLID%NWLID,1:KLEV), &
                  & ZOTO215(KDLEN,YREAERLID%NWLID,1:KLEV)
REAL(KIND=JPRB) :: ZOTRAY1 (KDLEN,YREAERLID%NWLID,1:KLEV), ZOTCO21 (KDLEN,YREAERLID%NWLID,1:KLEV), &
                  & ZOTO21 (KDLEN,YREAERLID%NWLID,1:KLEV)
REAL(KIND=JPRB) :: ZOTNO215(KDLEN,YREAERLID%NWLID,1:KLEV), ZOTO315(KDLEN,YREAERLID%NWLID,1:KLEV) , &
                  & ZPATHT15(KDLEN,YREAERLID%NWLID,0:KLEV)
REAL(KIND=JPRB) :: ZOTNO21 (KDLEN,YREAERLID%NWLID,1:KLEV), ZOTO31 (KDLEN,YREAERLID%NWLID,1:KLEV) , &
                  & ZPATHT1 (KDLEN,YREAERLID%NWLID,0:KLEV)
REAL(KIND=JPRB) :: ZPATHS15(KDLEN,YREAERLID%NWLID,0:KLEV)
REAL(KIND=JPRB) :: ZPATHS1 (KDLEN,YREAERLID%NWLID,0:KLEV)

REAL(KIND=JPRB) :: ZAERMSS15(KDLEN,KLEV,YDGP5%NGEMS)
REAL(KIND=JPRB) :: ZAERMSS1 (KDLEN,KLEV,YDGP5%NGEMS), ZALF(YDGP5%NGEMS), ZOMG(YDGP5%NGEMS), ZLIR(YDGP5%NGEMS)
REAL(KIND=JPRB) :: ZBSC15(KDLEN,YREAERLID%NWLID,KLEV)      , ZAEROD15(KDLEN,YREAERLID%NWLID,KLEV)
REAL(KIND=JPRB) :: ZBSC1 (KDLEN,YREAERLID%NWLID,KLEV)      , ZAEROD1 (KDLEN,YREAERLID%NWLID,KLEV)
REAL(KIND=JPRB) :: ZBSCMOL15(KDLEN,YREAERLID%NWLID,KLEV)   , ZBSCMOL1 (KDLEN,YREAERLID%NWLID,KLEV)
REAL(KIND=JPRB) :: ZBOLTZ, ZCONST

REAL(KIND=JPRB) ::                                          ZDPG15,          ZFACT15
REAL(KIND=JPRB) :: ZA1, ZA2, ZA3, ZA4, ZAA, ZAN, ZAN2, ZDPG1 , ZEPSAER, ZFAC, ZFACT1 , ZMD, ZLIDRAT

REAL(KIND=JPRB) :: ZNO2MOL15, ZNO2TAU15, ZNO2V15, ZO3MOL15, ZO3TAU15, ZO3V15
REAL(KIND=JPRB) :: ZNO2MOL1 , ZNO2TAU1 , ZNO2V1 , ZO3MOL1 , ZO3TAU1 , ZO3V1 ,  ZREF, ZSIGMA

REAL(KIND=JPRB) ::                             ZTAUR15, ZTE15, ZTE215
REAL(KIND=JPRB) :: ZTAUR1 , ZTE1 , ZTE21

REAL(KIND=JPRB) :: ZTRSOT15, ZTRTOT15, ZUSOT15, ZUTOT15
REAL(KIND=JPRB) :: ZTRSOT1 , ZTRTOT1 , ZUSOT1 , ZUTOT1 , ZWLUM, ZWLCM, ZWL2CM, ZWL4CM, ZWN, ZWN2

REAL(KIND=JPRB) :: ZGCO215,                       ZPHIC15, ZPSIC15, ZUCO215, ZUPCO215, ZVCO215
REAL(KIND=JPRB) :: ZGCO21 , ZMCO2,ZPHIC1 , ZPSIC1 , ZUCO21 , ZUPCO21 , ZVCO21

REAL(KIND=JPRB) ::       ZGO215,                        ZPHIO15, ZPSIO15, ZUO215, ZUPO215
REAL(KIND=JPRB) :: ZO2 , ZGO21 , ZMO2 , ZPHIO1 , ZPSIO1 , ZUO21 , ZUPO21

REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "satur.intfb.h"

!***********************************************************************************
IF (LHOOK) CALL DR_HOOK('AER_LIDSIMAD',0,ZHOOK_HANDLE)
ASSOCIATE( &
 & RSS_RH80_MASSFAC=>YREAERATM%RSS_RH80_MASSFAC, &
 & NVOLOPTP=>YREAERVOL%NVOLOPTP, &
 & JWLID=>YREAERLID%JWLID, NWLID=>YREAERLID%NWLID, RLICLS=>YREAERLID%RLICLS, &
 & RLICO2=>YREAERLID%RLICO2, RLIDELT=>YREAERLID%RLIDELT, &
 & RLINO2=>YREAERLID%RLINO2, RLINS=>YREAERLID%RLINS, RLIO2=>YREAERLID%RLIO2, &
 & RLIO3=>YREAERLID%RLIO3, RLIPREF=>YREAERLID%RLIPREF, RLIT0=>YREAERLID%RLIT0, &
 & RLITREF=>YREAERLID%RLITREF, RWLID=>YREAERLID%RWLID, &
 & RRHTAB=>YREAERSNK%RRHTAB, &
 & YAERO_DESC=>YREAERATM%YAERO_DESC)
!***********************************************************************************


!*         0.     SPECTRAL INFORMATION OF RELEVANCE
!                 ---------------------------------

! Rayleigh |  355nm / 28169cm-1 |  532nm / 18797cm-1 | 1064nm / 9398.5cm-1 |
! Index    |          2         |          8         |         16          ! 
! CO2      |         no         |         no         |        yes          |
! O2       |         no         |         no         |        yes          |
! O3       |        yes         |        yes         |        yes          |
! NO2      |        yes         |        yes         |        yes          |



!*         1.     INITIALIZATION AND PREPARATORY WORK
!                 -----------------------------------

!************************************************************************************
! (0)  COMPUTATIONS of LINEARISATION STATE VARIABLES
!************************************************************************************
ZREF=RLITREF/RLIPREF
ZO2=RMO3*2._JPRB/3._JPRB
ZEPSAER=1.E-15_JPRB
LLPRINT=.FALSE.
!-- the effective relative hunidity is the low value (20%) assumed for
!hydrophobic component of OM
IEFRH=3

! Extract aerosol mixing ratios from GEMS profiles structure (trajectory)
ZAERF15(:,:,:)=0.0_JPRB
DO JL=1,KLEN
  IAER=0
  DO JF=1,YDGP5%NGEMS
    JP=0
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR01 .OR.  &
     &   YDGP5%GEMS_IGRIB(JF) == NGRBAERMR02 .OR.  &
     &   YDGP5%GEMS_IGRIB(JF) == NGRBAERMR03) JP=1 !SS 3 bins
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR04 .OR.  &
     &   YDGP5%GEMS_IGRIB(JF) == NGRBAERMR05 .OR.  &
     &   YDGP5%GEMS_IGRIB(JF) == NGRBAERMR06 ) JP=2 !DD 3 bins
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR07 .OR. &
     &   YDGP5%GEMS_IGRIB(JF) == NGRBAERMR08) JP=3 ! ORGANIC MATTER - 2 bins
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR09 .OR. &
     &   YDGP5%GEMS_IGRIB(JF) == NGRBAERMR10) JP=4 ! Black Carbon - 2 bins
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR11 .OR. & 
     &   YDGP5%GEMS_IGRIB(JF) == NGRBAERMR12)  JP=5 ! SO4 & SO2
    IF((YDGP5%GEMS_IGRIB(JF) == NGRBAERMR16 .OR. &
     & YDGP5%GEMS_IGRIB(JF) == NGRBAERMR17)) JP=6 ! NO3
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR18) JP=7 ! NH4
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR19 .OR. &
     & YDGP5%GEMS_IGRIB(JF) == NGRBAERMR20 .OR. &
     & YDGP5%GEMS_IGRIB(JF) == NGRBAERMR21) JP=8 ! SOA
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR13) JP=9 ! Volcanic ashes 
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR14 .OR. &
     & YDGP5%GEMS_IGRIB(JF) == NGRBAERMR15)  JP=10 ! Volcanic SO4 and SO2

    IF(JP > 0) THEN
      IAER=IAER+1
      ZAERF15(JL,:,IAER)=PAELIDPROF5(JL,1:KLEV,JF)
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERLG) THEN ! Total aerosol mixing ratio
      ZAERTMR15(JL,:)= PAELIDPROF5(JL,1:KLEV,JF)
    ENDIF
    ! Zero out adjoint variables
    PAELIDPROF(JL,1:KLEV,JF) = 0.0_JPRB
  ENDDO
ENDDO
INACTAERO=IAER

! Compute total mixing ratio by summing single contributions
    ZAERFT15(:,:) = 0.0_JPRB
    DO JAER=1,INACTAERO 
      IF ((YAERO_DESC(JAER)%NTYP /= 5 .AND. YAERO_DESC(JAER)%NTYP /= 9) .OR. YAERO_DESC(JAER)%NBIN /= 2) THEN
        ZAERFT15(:,:)=ZAERFT15(:,:)+ ZAERF15(:,:,JAER)
      ENDIF
    ENDDO

! BEN the perturbations on the various species are computed from the perturbation
! on the total mass mixing ratio using the mass fraction assumed constant

ZFRACAER15(:,:,:) = 0.0_JPRB
DO JL=1,KLEN
  DO JAER=1,INACTAERO 
    IF ((YAERO_DESC(JAER)%NTYP /= 5 .AND. YAERO_DESC(JAER)%NTYP /= 9) .OR. YAERO_DESC(JAER)%NBIN /= 2) THEN
      DO JLEV=1,KLEV
        IF(ZAERFT15(JL,JLEV) /= 0.0_JPRB) THEN
          ZFRACAER15(JL,JLEV,JAER) = ZAERF15(JL,JLEV,JAER)/ZAERFT15(JL,JLEV)
        ENDIF
      ENDDO
    ENDIF
  ENDDO
ENDDO

! Save T, Q, and O3 profiles in separate arrays
ZTF15 (1:KLEN,:) = PTQO15(1:KLEN,:,1)
ZQF15 (1:KLEN,:) = PTQO15(1:KLEN,:,2)
ZO3F15(1:KLEN,:) = PTQO15(1:KLEN,:,3)

ZPF15(1:KLEN,:) = PLIDPRES(1:KLEN,:)

! COMPUTE THICKNESS

ZGEOMH15(:,0) = PGEOP15(:,1)  ! Top layer is assumed to have zero thickness...
ZGEOMH15(:,1:KLEV) = PGEOP15(:,1:KLEV)

DO JL=1,KLEN
 DO JK=1,KLEV
    ZDZ15(JL,JK) = (ZGEOMH15(JL,JK-1)-ZGEOMH15(JL,JK)) /RG
  ENDDO
    ZDZ15(JL,1) = 60.0_JPRB ! set top layer equal to 60 meters
ENDDO

!     INTEGRATE HYDROSTATIC EQUATION TO COMPUTE ZDP from ZDZ

DO JK=1,KLEV
  DO JL=1,KLEN
     ZDP15(JL,JK) = RG * ZDZ15(JL,JK) / ( 287.0_JPRB*ZTF15(JL,JK)* &
     & (1.0_JPRB+0.608_JPRB*ZQF15(JL,JK)) ) * ZPF15(JL,JK)
  ENDDO
ENDDO

IFLAG=2
CALL SATUR (YRTHF, YRCST, 1 , KLEN , KLEN  , 1 , KLEV, LDPHYLIN, &
  & ZPF15, ZTF15 , ZQSAT15 , IFLAG)

!-- define RH index from "clear-sky" relative humidity
IRH(:,:)=1
ZRHCL15= 0.0_JPRB
DO JK=1,KLEV
  DO JL=1,KLEN
    ZRHCL15=(ZQF15(JL,JK)/ZQSAT15(JL,JK))*100._JPRB
    DO JTAB=1,JTABMX
      IF (ZRHCL15 > RRHTAB(JTAB)) THEN
        IRH(JL,JK)=JTAB
      ENDIF
    ENDDO
  ENDDO
ENDDO

ZBSC15(1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZBSCMOL15(1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZAEROD15(1:KLEN,1:NWLID,1:KLEV)=0._JPRB

ZOTRAY15(1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZOTNO215(1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZOTO315 (1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZOTO215 (1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZOTCO215(1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZPATHT15(1:KLEN,1:NWLID,0:KLEV)=0._JPRB

ZPATHS15(1:KLEN,1:NWLID,0:KLEV)=0._JPRB



!*         2.      RAYLEIGH AND GASEOUS OPTICAL THICKNESSES
!                  ----------------------------------------

DO JWL=1,NWLID
!-- RWLID in m, ZWLUM in um, ZWLCM in cm!!!!
  ZWLUM=RWLID(JWL)*1.0E+06_JPRB
  ZWLCM=RWLID(JWL)*1.0E+02_JPRB
  ZWL2CM=ZWLCM*ZWLCM
  ZWL4CM=ZWL2CM*ZWL2CM

  ZWN=1._JPRB/ZWLUM
  ZWN2=ZWN*ZWN

!* Air Refractive Index: Edlen, 1966, Metrologia, 2, 71-80 w pw=0
  ZA1=130._JPRB-ZWN2
  ZA2=38.9_JPRB-ZWN2
  ZA3=2406030._JPRB/ZA1
  ZA4=15997._JPRB/ZA2
  ZAN=1._JPRB+(8342.13_JPRB+ZA3+ZA4)*1.E-08_JPRB
  ZAN2=ZAN*ZAN
  ZAA=(24._JPRB*RPI**3) * ((ZAN2-1._JPRB)/(ZAN2+2._JPRB))**2 * (6._JPRB+3._JPRB*RLIDELT)/(6._JPRB-7._JPRB*RLIDELT)

!-- According to Bohdaine et al. (1999), T and p dependence of the molecular density Ns and of the function of 
!   the refractive index (ns**2-1)/(ns**2+2) cancels out, and Rayleigh scattering coefficient ZSIGMA is independent 
!   of T and p; wavelength MUST be in cm to get ZSIGMA in cm2. 
!   ZTAUR for the entire atmosphere should be between 0.59 at 355 nm and < 0.0007 at 1064 nm.

  ZSIGMA=ZAA/(ZWL4CM*RLINS*RLINS)

  ZOTRAY15(1:KLEN,JWL,1:KLEV)=0._JPRB
  ZOTNO215(1:KLEN,JWL,1:KLEV)=0._JPRB
  ZOTO315 (1:KLEN,JWL,1:KLEV)=0._JPRB

  DO JK=1,KLEV
    DO JL=1,KLEN

!-- ZFACT is the number of gas molecules in a layer of pressure ZDP at temperature T
      ZFACT15=RNAVO*(ZDP15(JL,JK)*RLITREF)/(ZTF15(JL,JK)*RLIPREF)

!-- Rayleigh: Bodhaine et al., 1999: JAOT, 16, 1854-1861.   tau= sigma*DelP*Avo / (mmw_air * g)
      ZVCO215=PCO215(JL,JK)*RMD/RMCO2
      ZMD = RMD
      ZTAUR15=ZSIGMA*ZDP15(JL,JK)*RNAVO/(ZMD*RG)/10._JPRB
      ZOTRAY15(JL,JWL,JK)= ZTAUR15

!-- O3 
      ZO3V15=ZO3F15(JL,JK)*RMD/RMO3                              ! volume mixing ratio in layer
      ZO3MOL15=ZO3V15*ZFACT15                                     ! number of O3 molecules in layer (pressure unit: kg m-1 s-2)
      ZO3TAU15=1.E-04_JPRB*RLIO3(3,JWL)*ZO3MOL15/(RMO3*RG)    ! cross-section from cm2 molec-1 to m-2 molec-1
      ZOTO315(JL,JWL,JK)=ZO3TAU15

!-- NO2
      ZNO2V15=PNO215(JL,JK)*RMD/RMNO2                           ! volume mixing ratio in layer
      ZNO2MOL15=ZNO2V15*ZFACT15                                   ! number of NO2 molecules in layer (pressure unit: kg m-1 s-2)
      ZNO2TAU15=1.E-04_JPRB*RLINO2(3,JWL)*ZNO2MOL15/(RMNO2*RG)  ! cross-section from cm2 molec-1 to m-2 molec-1
      ZOTNO215(JL,JWL,JK)=ZNO2TAU15

    ENDDO
  ENDDO
ENDDO

!-- to start with (and to be replaced by more proper quantities once debugged!)
ZMCO2= 357.E-06_JPRB*RMCO2/RMD
ZMO2 = 0.20947_JPRB *ZO2  /RMD

ZOTCO215(1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZOTO215 (1:KLEN,1:NWLID,1:KLEV)=0._JPRB
DO JWL=1,NWLID
 DO JK=1,KLEV
  DO JL=1,KLEN
    ZTE15=ZTF15(JL,JK)-RLIT0
    ZTE215=ZTE15*ZTE15
    ZDPG15=ZDP15(JL,JK)/RG

!-- CO2 -------------------------------
    ZPHIC15=EXP(RLICO2(3)*ZTE15+RLICO2(4)*ZTE215)
    ZPSIC15=EXP(RLICO2(5)*ZTE15+RLICO2(6)*ZTE215)
    ZUCO215 =ZPHIC15*ZMCO2*ZDPG15
    ZUPCO215=ZPSIC15*ZMCO2*ZDPG15*ZPF15(JL,JK)/RLIPREF

!-- Goody model for CO2
    ZGCO215=RLICO2(1)*ZUCO215/SQRT(1._JPRB+RLICO2(1)*ZUCO215*ZUCO215/(RLICO2(2)*ZUPCO215))

    ZOTCO215(JL,JWL,JK)=ZGCO215

!-- O2 --------------------------------
    ZPHIO15=EXP(RLIO2(3) *ZTE15+RLIO2(4) *ZTE215)
    ZPSIO15=EXP(RLIO2(5) *ZTE15+RLIO2(6) *ZTE215)
    ZUO215 =ZPHIO15*ZMO2*ZDPG15
    ZUPO215=ZPSIO15*ZMO2*ZDPG15*ZPF15(JL,JK)/RLIPREF
!-- Goody model for O2
    ZGO215=RLIO2(1)*ZUO215/SQRT(1._JPRB+RLIO2(1)*ZUO215*ZUO215/(RLIO2(2)*ZUPO215))

    ZOTO215(JL,JWL,JK)=ZGO215

  ENDDO
 ENDDO
ENDDO

!*         3.      AEROSOL-RELATED QUANTITIES
!                  --------------------------
ZAERMSS15(:,:,:) = 0.0_JPRB
DO JAER=1,INACTAERO
  DO JK=1,KLEV
    DO JL=1,KLEN
      IF (ZAERF15(JL,JK,JAER) > ZEPSAER) THEN
        ZAERMSS15(JL,JK,JAER)= ZAERF15(JL,JK,JAER)*ZDP15(JL,JK)/RG
      ENDIF
    ENDDO
  ENDDO
ENDDO


!-- path from TOA
ZPATHS15(1:KLEN,1:NWLID,0:KLEV)=0._JPRB
ZPATHT15(1:KLEN,1:NWLID,0:KLEV)=0._JPRB


DO JWL=1,NWLID
  IWAVL=JWLID(JWL)

  DO JAER=1,INACTAERO

    ITYP=YAERO_DESC(JAER)%NTYP
    IBIN=YAERO_DESC(JAER)%NBIN

!-- ITYP is the aerosol type 1:SS,   2:DD,   3:OM    4:BC,   5:SU,   6:NI,   7:AM,   8=SOA,  9=VA, 10=VS
!   IBIN is the bin index: 1-3:SS, 1-3:DD,   2:OM,   2:BC,   2:SU,   1-2:NI, 1:AM,   1-3:SOA,1:VA, 1-2:VS
!   N.B.: extinction coefficients are in m2 g-1

    DO JK=1,KLEV
      DO JL=1,KLEN
        ZFAC = 1.0_JPRB
        IF (ITYP == 1) THEN
          ZALF(JAER)=ALF_SS(IRH(JL,JK),IWAVL,IBIN)
          ZOMG(JAER)=OMG_SS(IRH(JL,JK),IWAVL,IBIN)
          ZLIR(JAER)=RALI_SS(IRH(JL,JK),IWAVL,IBIN)
          ZFAC = RSS_RH80_MASSFAC
        ELSEIF (ITYP == 2) THEN
          ZALF(JAER)=ALF_DD(IBIN,IWAVL)
          ZOMG(JAER)=OMG_DD(IBIN,IWAVL)
          ZLIR(JAER)=RALI_DD(IBIN,IWAVL)
        ELSEIF (ITYP == 3) THEN
          ZALF(JAER)=ALF_OM(IRH(JL,JK),IWAVL)
          ZOMG(JAER)=OMG_OM(IRH(JL,JK),IWAVL)
          ZLIR(JAER)=RALI_OM(IRH(JL,JK),IWAVL)
        ELSEIF (ITYP == 4) THEN
          ZALF(JAER)=ALF_BC(IWAVL)
          ZOMG(JAER)=OMG_BC(IWAVL)
          ZLIR(JAER)=RALI_BC(IWAVL)
        ELSEIF (ITYP == 5 .OR. ITYP == 10) THEN
          ZALF(JAER)=ALF_SU(IRH(JL,JK),IWAVL)
          ZOMG(JAER)=OMG_SU(IRH(JL,JK),IWAVL)
          ZLIR(JAER)=RALI_SU(IRH(JL,JK),IWAVL)
!-- SO2 does not contribute to optical depth, only SO4 does.
          IF (IBIN == 2) THEN
            ZALF(JAER)=0._JPRB
            ZOMG(JAER)=0._JPRB
            ZLIR(JAER)=0._JPRB
          ENDIF
        ELSEIF (ITYP == 6) THEN
          ZALF(JAER)=ALF_NI(IRH(JL,JK),IWAVL,IBIN)
          ZOMG(JAER)=OMG_NI(IRH(JL,JK),IWAVL,IBIN)
          ZLIR(JAER)=RALI_NI(IRH(JL,JK),IWAVL,IBIN)
        ELSEIF (ITYP == 7) THEN
          ZALF(JAER)=ALF_AM(IRH(JL,JK),IWAVL)
          ZOMG(JAER)=OMG_AM(IRH(JL,JK),IWAVL)
          ZLIR(JAER)=RALI_AM(IRH(JL,JK),IWAVL)
        ELSEIF (ITYP == 8) THEN
          ZALF(JAER)=ALF_SOA(IRH(JL,JK),IWAVL,IBIN)
          ZOMG(JAER)=OMG_SOA(IRH(JL,JK),IWAVL,IBIN)
          ZLIR(JAER)=RALI_SOA(IRH(JL,JK),IWAVL,IBIN)
        ELSEIF (ITYP ==9) THEN
          IF (NVOLOPTP == 1) THEN
            ZALF(JAER)= ALF_SU(IEFRH, IWAVL)
            ZOMG(JAER)= OMG_SU(IEFRH, IWAVL)
            ZLIR(JAER)=RALI_SU(IEFRH, IWAVL)
          ELSEIF (NVOLOPTP == 2) THEN
            ZALF(JAER)= ALF_BC(IWAVL)
            ZOMG(JAER)= OMG_BC(IWAVL)
            ZLIR(JAER)=RALI_BC(IWAVL)
          ELSEIF (NVOLOPTP == 3) THEN
            ZALF(JAER)= ALF_DD(3,IWAVL)
            ZOMG(JAER)= OMG_DD(3,IWAVL)
            ZLIR(JAER)=RALI_DD(3,IWAVL)
          ENDIF
        ENDIF
        IF (ZLIR(JAER) /= 0._JPRB) THEN
          ZLIDRAT=1._JPRB/ZLIR(JAER)
        ELSE
          ZLIDRAT=ZLIR(JAER)
        ENDIF

!-- extinction within a layer (sum on aerosol contribution)
        ZAEROD15(JL,JWL,JK) = ZAEROD15(JL,JWL,JK) + ZAERMSS15(JL,JK,JAER) * ZFAC * 1.E+03_JPRB * ZALF(JAER)

!-- backscatter within a layer
        ZBSC15(JL,JWL,JK) = ZBSC15(JL,JWL,JK) + ZAERMSS15(JL,JK,JAER) * ZFAC * 1.E+03_JPRB * ZALF(JAER) * ZLIDRAT

      ENDDO
    ENDDO
  ENDDO

ZBOLTZ = 1.38E-23_JPRB ! (J/K, Boltzman constant)
! Calculation of molecular backscatter
 ZWLUM=RWLID(JWL)*1.0E+06_JPRB ! wavelenght in micrometers
  DO JK=1,KLEV
   DO JL=1,KLEN
! Molecular backscatter from Collis and Russel (1976)
! Trajectory
      ZBSCMOL15(JL,JWL,JK) = 5.45E-32_JPRB * ZPF15(JL,JK)/(ZBOLTZ * ZTF15(JL,JK)) * (ZWLUM/0.55_JPRB)**(-4.09_JPRB)
    ENDDO
   ENDDO

  ZPATHS15(1:KLEN,JWL,KLEV)=0._JPRB
  ZPATHT15(1:KLEN,JWL,0   )=0._JPRB

  DO JK=1,KLEV
    IK=KLEV-JK
    DO JL=1,KLEN
!-- lidar signal seen from TOA
      ZPATHT15(JL,JWL,JK)=ZPATHT15(JL,JWL,JK-1)+2._JPRB*(ZAEROD15(JL,JWL,JK) &
      & +ZOTRAY15(JL,JWL,JK)+ZOTO315(JL,JWL,JK)+ZOTNO215(JL,JWL,JK)+ZOTO215(JL,JWL,JK)+ZOTCO215(JL,JWL,JK))
      ZUTOT15=200._JPRB
      IF(ZPATHT15(JL,JWL,JK) < 200._JPRB) THEN
        ZUTOT15=ZPATHT15(JL,JWL,JK)
      ENDIF
      ZTRTOT15=EXP(-ZUTOT15)
      ZLISIT5(JL,JWL,JK)=RLICLS* (ZBSC15(JL,JWL,JK)/ZDZ15(JL,JK) + ZBSCMOL15(JL,JWL,JK) ) * ZTRTOT15
!-- lidar signal seen from the surface
      ZPATHS15(JL,JWL,IK)=ZPATHS15(JL,JWL,IK+1)+2._JPRB*(ZAEROD15(JL,JWL,IK+1) &
      & +ZOTRAY15(JL,JWL,IK+1)+ZOTO315(JL,JWL,IK+1)+ZOTNO215(JL,JWL,IK+1)+ZOTO215(JL,JWL,IK+1)+ZOTCO215(JL,JWL,IK+1))

      ZUSOT15=200._JPRB
      IF(ZPATHS15(JL,JWL,IK) < 200._JPRB) THEN
        ZUSOT15=ZPATHS15(JL,JWL,IK)
      ENDIF
      ZTRSOT15=EXP(-ZUSOT15)
      ZLISIS5(JL,JWL,IK)=RLICLS* (ZBSC15(JL,JWL,IK+1)/ZDZ15(JL,IK+1) + ZBSCMOL15(JL,JWL,IK+1) ) * ZTRTOT15
    ENDDO
  ENDDO
ENDDO

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! (1) STARTING PERTURBATION COMPUTATIONS
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!========================================================================
! (1a) setting all gradients to zero
!========================================================================

ZAERF1(:,:,:) = 0._JPRB
ZDP1(:,:) = 0._JPRB
ZDZ1(:,:) = 0._JPRB
ZO3F1(:,:)  = 0._JPRB
ZTF1(:,:) = 0._JPRB
ZQF1(:,:) = 0._JPRB
ZPF1(:,:) = 0._JPRB
ZGEOMH1(:,:) = 0._JPRB

ZOTRAY1(:,:,:)= 0._JPRB
ZOTCO21(:,:,:)= 0._JPRB
ZOTO21(:,:,:) = 0._JPRB
ZOTNO21(:,:,:)= 0._JPRB
ZOTO31(:,:,:) = 0._JPRB

ZPATHT1(:,:,:)= 0._JPRB
ZPATHS1(:,:,:) = 0._JPRB

ZAERMSS1(:,:,:) = 0._JPRB

ZBSC1(:,:,:)      =  0._JPRB
ZBSCMOL1(:,:,:)      =  0._JPRB

ZAEROD1(:,:,:) = 0._JPRB
ZDPG1= 0._JPRB
ZFACT1 = 0._JPRB

ZNO2MOL1= 0._JPRB
ZNO2TAU1= 0._JPRB
ZNO2V1= 0._JPRB
ZO3MOL1= 0._JPRB
ZO3TAU1= 0._JPRB
ZO3V1 = 0._JPRB

ZTAUR1= 0._JPRB
ZTE1= 0._JPRB
ZTE21 = 0._JPRB

ZTRSOT1= 0._JPRB
ZTRTOT1= 0._JPRB
ZUSOT1= 0._JPRB
ZUTOT1 = 0._JPRB

ZGCO21= 0._JPRB
ZPHIC1= 0._JPRB
ZPSIC1= 0._JPRB
ZUCO21= 0._JPRB
ZUPCO21= 0._JPRB
ZVCO21 = 0._JPRB

ZGO21= 0._JPRB
ZPHIO1= 0._JPRB
ZPSIO1= 0._JPRB
ZUO21= 0._JPRB
ZUPO21 = 0._JPRB


!======================================================================================================
! (1b) main part of perturbation computations
!======================================================================================================
DO JWL=1,NWLID
  ZPATHS1(1:KLEN,JWL,KLEV)=0._JPRB
  ZPATHT1(1:KLEN,JWL,0   )=0._JPRB
  DO JK=KLEV,1,-1
    IK=KLEV-JK
    DO JL=1,KLEN
      !========================================================================
      ! local linearisation state variables
      !========================================================================
      ZUSOT15=200._JPRB
      IF(ZPATHS15(JL,JWL,IK) < 200._JPRB) THEN
        ZUSOT15=ZPATHS15(JL,JWL,IK)
      ENDIF
      ZTRSOT15=EXP(-ZUSOT15)
      !========================================================================
      ! perturbation computations
      !========================================================================
      ZTRSOT1=ZTRSOT1+PEXTS(JL,JWL,IK)  ! extinction contribution 
      
      ZBSC1(JL,JWL,IK+1)=ZBSC1(JL,JWL,IK+1)+RLICLS/ZDZ15(JL,IK+1)                             *PULISIS (JL,JWL,IK) ! unattenuated backscatter
      ZDZ1(JL,IK+1) = ZDZ1(JL,IK+1) - RLICLS/ZDZ15(JL,IK+1)**2 * ZBSC15(JL,JWL,IK+1)          *PULISIS (JL,JWL,IK)

      ZTRSOT1=ZTRSOT1+RLICLS * (ZBSC15(JL,JWL,IK+1)/ZDZ15(JL,IK+1) + ZBSCMOL15(JL,JWL,IK+1) ) *PLISIS (JL,JWL,IK)  ! attenuated backscatter
      ZBSC1(JL,JWL,IK+1)=ZBSC1(JL,JWL,IK+1)+RLICLS/ZDZ15(JL,IK+1)*ZTRSOT15            *PLISIS (JL,JWL,IK)
      ZDZ1(JL,IK+1) = ZDZ1(JL,IK+1) - RLICLS/ZDZ15(JL,IK+1)**2 * ZBSC15(JL,JWL,IK+1)*ZTRSOT15 *PLISIS (JL,JWL,IK)
! No perturbations in molecolar backscatter
!      ZBSCMOL1(JL,JWL,IK+1) =  ZBSCMOL1(JL,JWL,IK+1) + RLICLS* ZTRSOT15 * PLISIS (JL,JWL,IK)

      ZUSOT1 = ZUSOT1 -ZTRSOT15 * ZTRSOT1
      ZTRSOT1= 0.0_JPRB 

      IF(ZPATHS15(JL,JWL,IK) < 200._JPRB) THEN
        ZPATHS1(JL,JWL,IK)= ZPATHS1(JL,JWL,IK) + ZUSOT1
        ZUSOT1 = 0.0_JPRB
      ENDIF
      ZUSOT1 = 0.0_JPRB

      ZPATHS1 (JL,JWL,IK+1)  = ZPATHS1(JL,JWL,IK+1)  +        ZPATHS1(JL,JWL,IK)
      ZAEROD1 (JL,JWL,IK+1)  = ZAEROD1(JL,JWL,IK+1)  +2._JPRB*ZPATHS1(JL,JWL,IK)
      ZOTRAY1 (JL,JWL,IK+1)  = ZOTRAY1(JL,JWL,IK+1)  +2._JPRB*ZPATHS1(JL,JWL,IK)
      ZOTO31  (JL,JWL,IK+1)  = ZOTO31 (JL,JWL,IK+1)  +2._JPRB*ZPATHS1(JL,JWL,IK)
      ZOTNO21 (JL,JWL,IK+1)  = ZOTNO21(JL,JWL,IK+1)  +2._JPRB*ZPATHS1(JL,JWL,IK)
      ZOTO21  (JL,JWL,IK+1)  = ZOTO21 (JL,JWL,IK+1)  +2._JPRB*ZPATHS1(JL,JWL,IK)
      ZOTCO21 (JL,JWL,IK+1)  = ZOTCO21(JL,JWL,IK+1)  +2._JPRB*ZPATHS1(JL,JWL,IK)
      ZPATHS1 (JL,JWL,IK)   = 0.0_JPRB
      
!--  end of adjoint calculations for lidar signal seen from the surface
    !-------------------------------------------------------------------------------------
      !========================================================================
      ! local linearisation state variables
      !========================================================================
      ZUTOT15=200._JPRB
      IF(ZPATHT15(JL,JWL,JK) < 200._JPRB) THEN
        ZUTOT15=ZPATHT15(JL,JWL,JK)
      ENDIF
      ZTRTOT15=EXP(-ZUTOT15)
      !========================================================================
      ! perturbation computations
      !========================================================================
      ZTRTOT1           =ZTRTOT1 + PEXTT(JL,JWL,JK) !extinction from TOA 
  
      ZBSC1(JL,JWL,JK)=ZBSC1(JL,JWL,JK)+RLICLS/ZDZ15(JL,JK)                           *PULISIT (JL,JWL,JK) !unattenuated backscatter
      ZDZ1(JL,JK) = ZDZ1(JL,JK) - RLICLS/ZDZ15(JL,JK)**2 * ZBSC15(JL,JWL,JK)          *PULISIT (JL,JWL,JK)
      ZTRTOT1           =ZTRTOT1           +RLICLS * (ZBSC15(JL,JWL,JK)/ZDZ15(JL,JK) + ZBSCMOL15(JL,JWL,JK) ) *PLISIT (JL,JWL,JK)
      ZBSC1(JL,JWL,JK)=ZBSC1(JL,JWL,JK)+RLICLS/ZDZ15(JL,JK)*ZTRTOT15            *PLISIT (JL,JWL,JK)
      ZDZ1(JL,JK) = ZDZ1(JL,JK) - RLICLS/ZDZ15(JL,JK)**2 * ZBSC15(JL,JWL,JK)*ZTRTOT15 *PLISIT (JL,JWL,JK)
! No perturbations in molecular backscatter
!      ZBSCMOL1(JL,JWL,JK) =  ZBSCMOL1(JL,JWL,JK) + RLICLS* ZTRTOT15 * PLISIT (JL,JWL,JK)

      ZUTOT1 = ZUTOT1 -ZTRTOT15 * ZTRTOT1
      ZTRTOT1= 0._JPRB

      IF(ZPATHT15(JL,JWL,JK) < 200._JPRB) THEN
        ZPATHT1(JL,JWL,JK)= ZPATHT1(JL,JWL,JK) + ZUTOT1
        ZUTOT1= 0._JPRB
      ENDIF
      ZUTOT1 =0._JPRB

      ZPATHT1(JL,JWL,JK-1) =ZPATHT1 (JL,JWL,JK-1)+        ZPATHT1(JL,JWL,JK)
      ZAEROD1 (JL,JWL,JK)    = ZAEROD1  (JL,JWL,JK)  +2._JPRB*ZPATHT1(JL,JWL,JK)
      ZOTRAY1 (JL,JWL,JK)  = ZOTRAY1(JL,JWL,JK)  +2._JPRB*ZPATHT1(JL,JWL,JK)
      ZOTO31  (JL,JWL,JK)  = ZOTO31 (JL,JWL,JK)  +2._JPRB*ZPATHT1(JL,JWL,JK)
      ZOTNO21 (JL,JWL,JK)  = ZOTNO21(JL,JWL,JK)  +2._JPRB*ZPATHT1(JL,JWL,JK)
      ZOTO21  (JL,JWL,JK)  = ZOTO21 (JL,JWL,JK)  +2._JPRB*ZPATHT1(JL,JWL,JK)
      ZOTCO21 (JL,JWL,JK)  = ZOTCO21(JL,JWL,JK)  +2._JPRB*ZPATHT1(JL,JWL,JK)
      ZPATHT1(JL,JWL,JK)   = 0._JPRB
    !-------------------------------------------------------------------------------------
!-- end of adjoint calculations for lidar signal seen from TOA 

    ENDDO
  ENDDO

! Molecular backscattering 
  ZBOLTZ = 1.38E-23_JPRB ! (J/K, Boltzman constant)
  ZWLUM=RWLID(JWL)*1.0E+06_JPRB ! wavelenght in micrometers
  ZCONST= 5.45E-32_JPRB /ZBOLTZ * (ZWLUM/0.55_JPRB)**(-4.09_JPRB) 
    DO JK=KLEV,1,-1
      DO JL=1,KLEN

        ZTF1(JL,JK) = ZTF1(JL,JK) + ZCONST * 1._JPRB/ZTF15(JL,JK) *  ZBSCMOL1(JL,JWL,JK)
        ZPF1(JL,JK) = ZPF1(JL,JK) - ZCONST * ZPF15(JL,JK)/ZTF15(JL,JK)**2 * ZBSCMOL1(JL,JWL,JK)  

      ENDDO
    ENDDO

  IWAVL=JWLID(JWL)

  DO JAER=1,INACTAERO

    ITYP=YAERO_DESC(JAER)%NTYP
    IBIN=YAERO_DESC(JAER)%NBIN

!-- ITYP is the aerosol type 1:SS,   2:DD,   3:OM    4:BC,   5:SU,   6:FA,   7:BS,   8=VS,
!   IBIN is the bin index: 1-3:SS, 1-3:DD,   2:OM,   2:BC,   2:SU,   1:FA,   1:SB    1=VS,
!   N.B.: extinction coefficients are in m2 g-1

    DO JK=KLEV,1,-1
      DO JL=1,KLEN
        ZFAC = 1.0_JPRB
        IF (ITYP == 1) THEN
          ZALF(JAER)=ALF_SS(IRH(JL,JK),IWAVL,IBIN)
          ZOMG(JAER)=OMG_SS(IRH(JL,JK),IWAVL,IBIN)
          ZLIR(JAER)=RALI_SS(IRH(JL,JK),IWAVL,IBIN)
          ZFAC = RSS_RH80_MASSFAC
        ELSEIF (ITYP == 2) THEN
          ZALF(JAER)=ALF_DD(IBIN,IWAVL)
          ZOMG(JAER)=OMG_DD(IBIN,IWAVL)
          ZLIR(JAER)=RALI_DD(IBIN,IWAVL)
        ELSEIF (ITYP == 3) THEN
          ZALF(JAER)=ALF_OM(IRH(JL,JK),IWAVL)
          ZOMG(JAER)=OMG_OM(IRH(JL,JK),IWAVL)
          ZLIR(JAER)=RALI_OM(IRH(JL,JK),IWAVL)
        ELSEIF (ITYP == 4) THEN
          ZALF(JAER)=ALF_BC(IWAVL)
          ZOMG(JAER)=OMG_BC(IWAVL)
          ZLIR(JAER)=RALI_BC(IWAVL)
        ELSEIF (ITYP == 5 .OR. ITYP == 10) THEN
          ZALF(JAER)=ALF_SU(IRH(JL,JK),IWAVL)
          ZOMG(JAER)=OMG_SU(IRH(JL,JK),IWAVL)
          ZLIR(JAER)=RALI_SU(IRH(JL,JK),IWAVL)
!-- SO2 does not contribute to optical depth, only SO4 does.
          IF (IBIN == 2) THEN
            ZALF(JAER)=0._JPRB
            ZOMG(JAER)=0._JPRB
            ZLIR(JAER)=0._JPRB
          ENDIF
        ELSEIF (ITYP == 6) THEN
          ZALF(JAER)=ALF_NI(IRH(JL,JK),IWAVL,IBIN)
          ZOMG(JAER)=OMG_NI(IRH(JL,JK),IWAVL,IBIN)
          ZLIR(JAER)=RALI_NI(IRH(JL,JK),IWAVL,IBIN)
        ELSEIF (ITYP == 7) THEN
          ZALF(JAER)=ALF_AM(IRH(JL,JK),IWAVL)
          ZOMG(JAER)=OMG_AM(IRH(JL,JK),IWAVL)
          ZLIR(JAER)=RALI_AM(IRH(JL,JK),IWAVL)
        ELSEIF (ITYP == 8) THEN
          ZALF(JAER)=ALF_SOA(IRH(JL,JK),IWAVL,IBIN)
          ZOMG(JAER)=OMG_SOA(IRH(JL,JK),IWAVL,IBIN)
          ZLIR(JAER)=RALI_SOA(IRH(JL,JK),IWAVL,IBIN)
        ELSEIF (ITYP ==9) THEN
          IF (NVOLOPTP == 1) THEN
            ZALF(JAER)= ALF_SU(IEFRH, IWAVL)
            ZOMG(JAER)= OMG_SU(IEFRH, IWAVL)
            ZLIR(JAER)=RALI_SU(IEFRH, IWAVL)
          ELSEIF (NVOLOPTP == 2) THEN
            ZALF(JAER)= ALF_BC(IWAVL)
            ZOMG(JAER)= OMG_BC(IWAVL)
            ZLIR(JAER)=RALI_BC(IWAVL)
          ELSEIF (NVOLOPTP == 3) THEN
            ZALF(JAER)= ALF_DD(3,IWAVL)
            ZOMG(JAER)= OMG_DD(3,IWAVL)
            ZLIR(JAER)=RALI_DD(3,IWAVL)
          ENDIF
        ELSE
          WRITE(*,*) 'ITYP=',ITYP,JAER,'aer_lidsimad.F90:XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        ENDIF
        IF (ZLIR(JAER) /= 0._JPRB) THEN
          ZLIDRAT=1._JPRB/ZLIR(JAER)
        ELSE
          ZLIDRAT=ZLIR(JAER)
        ENDIF

!-- extinction within a layer (sum on aerosol contribution)
       !----------
        ZAERMSS1 (JL,JK,JAER) = ZAERMSS1 (JL,JK,JAER) + 1.E+03_JPRB*ZALF(JAER)*ZFAC*ZAEROD1 (JL,JWL,JK)
       !----------

!-- back-scatter within a layer
       !----------
        ZAERMSS1(JL,JK,JAER) = ZAERMSS1 (JL,JK,JAER) + 1.E+03_JPRB*ZALF(JAER)*ZFAC*ZLIDRAT * ZBSC1 (JL,JWL,JK)
       !----------

      ENDDO
    ENDDO
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!*         3.      AEROSOL-RELATED QUANTITIES
!                  --------------------------

DO JAER=1,INACTAERO
  DO JK=KLEV,1,-1
    DO JL=1,KLEN
      IF (ZAERF15(JL,JK,JAER) > ZEPSAER) THEN 
        ZAERF1 (JL,JK,JAER) = ZAERF1(JL,JK,JAER) + ZDP15(JL,JK)/RG        * ZAERMSS1(JL,JK,JAER)
        ZDP1(JL,JK) =         ZDP1(JL,JK)        + ZAERF15(JL,JK,JAER)/RG *  ZAERMSS1(JL,JK,JAER)
        ZAERMSS1(JL,JK,JAER)= 0._JPRB
      ENDIF
    ENDDO
  ENDDO
ENDDO

ZMCO2= 357.E-06_JPRB*RMCO2/RMD
ZMO2 = 0.20947_JPRB *ZO2  /RMD

DO JWL=1,NWLID
DO JK=KLEV,1,-1
  DO JL=1,KLEN
   !========================================================================
   ! local linearisation state variables
   !========================================================================
    ZTE15=ZTF15(JL,JK)-RLIT0
    ZTE215=ZTE15*ZTE15
    ZDPG15=ZDP15(JL,JK)/RG

!-- CO2 -------------------------------
    ZPHIC15=EXP(RLICO2(3)*ZTE15+RLICO2(4)*ZTE215)
    ZPSIC15=EXP(RLICO2(5)*ZTE15+RLICO2(6)*ZTE215)

    ZUCO215 =ZPHIC15*ZMCO2*ZDPG15
    ZUPCO215=ZPSIC15*ZMCO2*ZDPG15*ZPF15(JL,JK)/RLIPREF

!-- Goody model for CO2
    ZGCO215=RLICO2(1)*ZUCO215/SQRT(1._JPRB+RLICO2(1)*ZUCO215*ZUCO215/(RLICO2(2)*ZUPCO215))
!-- O2 --------------------------------
    ZPHIO15=EXP(RLIO2(3) *ZTE15+RLIO2(4) *ZTE215)
    ZPSIO15=EXP(RLIO2(5) *ZTE15+RLIO2(6) *ZTE215)

    ZUO215 =ZPHIO15*ZMO2*ZDPG15
    ZUPO215=ZPSIO15*ZMO2*ZDPG15*ZPF15(JL,JK)/RLIPREF
!-- Goody model for O2
    ZGO215=RLIO2(1)*ZUO215/SQRT(1._JPRB+RLIO2(1)*ZUO215*ZUO215/(RLIO2(2)*ZUPO215))

   !========================================================================
   ! perturbation computations
   !========================================================================

    ZGO21               = ZGO21              +ZOTO21 (JL,JWL,JK)
    ZOTO21 (JL,JWL,JK)  =0.0_JPRB

    ZUO21 = ZUO21 +        ZGO215/ZUO215                                                       *ZGO21
    ZUO21 = ZUO21-         ZGO215/ZUO215 /((RLIO2(2)*ZUPO215)/(RLIO2(1)*ZUO215*ZUO215)+1._JPRB)*ZGO21
    ZUPO21=ZUPO21+0.5_JPRB*ZGO215/ZUPO215/((RLIO2(2)*ZUPO215)/(RLIO2(1)*ZUO215*ZUO215)+1._JPRB)*ZGO21
    ZGO21 =0.0_JPRB
                 
    ZPSIO1      =ZPSIO1      + ZUPO215/ZPSIO15      * ZUPO21
    ZDPG1       = ZDPG1      + ZUPO215/ZDPG15       * ZUPO21
    ZUPO21      =0.0_JPRB

    ZPHIO1 = ZPHIO1 + ZUO215/ZPHIO15  * ZUO21
    ZDPG1  = ZDPG1  + ZUO215/ZDPG15   * ZUO21
    ZUO21  =0.0_JPRB

    ZTE1 = ZTE1 + ZPSIO15*RLIO2(5) * ZPSIO1
    ZTE21=ZTE21 + ZPSIO15*RLIO2(6) * ZPSIO1
    ZPSIO1 =0.0_JPRB

    ZTE1 = ZTE1 + ZPHIO15*RLIO2(3) * ZPHIO1
    ZTE21=ZTE21 + ZPHIO15*RLIO2(4) * ZPHIO1
    ZPHIO1 =0.0_JPRB

    ZGCO21                = ZGCO21                + ZOTCO21 (JL,JWL,JK)
    ZOTCO21 (JL,JWL,JK)   =0.0_JPRB

    ZUCO21 = ZUCO21 +         ZGCO215/ZUCO215                                                             *ZGCO21 
    ZUCO21 = ZUCO21 -         ZGCO215/ZUCO215  /((RLICO2(2)*ZUPCO215)/(RLICO2(1)*ZUCO215*ZUCO215)+1._JPRB)*ZGCO21 
    ZUPCO21=ZUPCO21 +0.5_JPRB*ZGCO215/ZUPCO215 /((RLICO2(2)*ZUPCO215)/(RLICO2(1)*ZUCO215*ZUCO215)+1._JPRB)*ZGCO21
    ZGCO21 =0.0_JPRB

    ZPSIC1     = ZPSIC1    + ZUPCO215/ZPSIC15      * ZUPCO21
    ZDPG1      = ZDPG1     + ZUPCO215/ZDPG15       * ZUPCO21
    ZPF1(JL,JK)=ZPF1(JL,JK)+ ZUPCO215/ZPF15(JL,JK) * ZUPCO21
    ZUPCO21    =0.0_JPRB

    ZPHIC1 = ZPHIC1 + ZUCO215/ZPHIC15 * ZUCO21
    ZDPG1  = ZDPG1  + ZUCO215/ZDPG15  * ZUCO21
    ZUCO21 =0.0_JPRB

    ZTE21= ZTE21 +ZPSIC15*RLICO2(6)* ZPSIC1
    ZTE1 = ZTE1  +ZPSIC15*RLICO2(5)* ZPSIC1
    ZPSIC1=0.0_JPRB
    
    ZTE21= ZTE21 +ZPHIC15*RLICO2(4)* ZPHIC1
    ZTE1 = ZTE1  +ZPHIC15*RLICO2(3)* ZPHIC1
    ZPHIC1=0.0_JPRB

    ZDP1 (JL,JK)=ZDP1 (JL,JK)+ZDPG1/RG
    ZDPG1       =0.0_JPRB

    ZTE1=ZTE1+2._JPRB*ZTE15*ZTE21
    ZTE21=0.0_JPRB

    ZTF1(JL,JK)=ZTF1(JL,JK)+ZTE1
    ZTE1=0.0_JPRB

  ENDDO
 ENDDO
ENDDO

ZOTCO21 (1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZOTO21  (1:KLEN,1:NWLID,1:KLEV)=0._JPRB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!*         2.      RAYLEIGH AND GASEOUS OPTICAL THICKNESSES
!                  ----------------------------------------

DO JWL=1,NWLID
!-- RWLID in m, ZWLUM in um, ZWLCM in cm!!!!
  ZWLUM=RWLID(JWL)*1.0E+06_JPRB
  ZWLCM=RWLID(JWL)*1.0E+02_JPRB
  ZWL2CM=ZWLCM*ZWLCM
  ZWL4CM=ZWL2CM*ZWL2CM

  ZWN=1._JPRB/ZWLUM
  ZWN2=ZWN*ZWN

!* Air Refractive Index: Edlen, 1966, Metrologia, 2, 71-80 w pw=0
  ZA1=130._JPRB-ZWN2
  ZA2=38.9_JPRB-ZWN2
  ZA3=2406030._JPRB/ZA1
  ZA4=15997._JPRB/ZA2
  ZAN=1._JPRB+(8342.13_JPRB+ZA3+ZA4)*1.E-08_JPRB
  ZAN2=ZAN*ZAN
  ZAA=(24._JPRB*RPI**3) * ((ZAN2-1._JPRB)/(ZAN2+2._JPRB))**2 * (6._JPRB+3._JPRB*RLIDELT)/(6._JPRB-7._JPRB*RLIDELT)
  ZSIGMA=ZAA/(ZWL4CM*RLINS*RLINS)

  DO JK=KLEV,1,-1
    DO JL=1,KLEN
     !========================================================================
     ! local linearisation state variables
     !========================================================================
      ZNO2V15=PNO215(JL,JK)*RMD/RMNO2                           ! volume mixing ratio in layer
      ZNO2MOL15=ZNO2V15*ZFACT15                                   ! number of NO2 molecules in layer (pressure unit: kg m-1 s-2)
      ZNO2TAU15=1.E-04_JPRB*RLINO2(3,JWL)*ZNO2MOL15/(RMNO2*RG)  ! cross-section from cm2 molec-1 to m-2 molec-1
      ZFACT15=RNAVO*(ZDP15(JL,JK)*RLITREF)/(ZTF15(JL,JK)*RLIPREF)
!-- NO2
      ZO3V15=ZO3F15(JL,JK)*RMD/RMO3                              ! volume mixing ratio in layer
      ZO3MOL15=ZO3V15*ZFACT15                                     ! number of O3 molecules in layer (pressure unit: kg m-1 s-2)
      ZO3TAU15=1.E-04_JPRB*RLIO3(3,JWL)*ZO3MOL15/(RMO3*RG)    ! cross-section from cm2 molec-1 to m-2 molec-1
!-- O3 
      ZVCO215=PCO215(JL,JK)*RMD/RMCO2
      ZMD = RMD
      ZTAUR15=ZSIGMA*ZDP15(JL,JK)*RNAVO/(ZMD*RG)/10._JPRB


     !========================================================================
     ! perturbation computations
     !========================================================================

      ZNO2TAU1             =ZNO2TAU1             +ZOTNO21 (JL,JWL,JK)
      ZOTNO21 (JL,JWL,JK)  =0.0_JPRB
  
      ZNO2MOL1 =ZNO2MOL1 + 1.E-04_JPRB*RLINO2(3,JWL)/(RMNO2*RG) * ZNO2TAU1
      ZNO2TAU1 =0.0_JPRB

      ZFACT1   =ZFACT1 + ZNO2V15*ZNO2MOL1
      ZNO2V1   =ZNO2V1 + ZFACT15*ZNO2MOL1
      ZNO2MOL1 =0.0_JPRB

      PNO21 (JL,JK)=PNO21 (JL,JK) + RMD/RMNO2 * ZNO2V1
      ZNO2V1 =0.0_JPRB

      ZO3TAU1             =ZO3TAU1             +ZOTO31(JL,JWL,JK)
      ZOTO31(JL,JWL,JK)   =0.0_JPRB

      ZO3MOL1 = ZO3MOL1 + 1.E-04_JPRB*RLIO3(3,JWL)/(RMO3*RG) * ZO3TAU1
      ZO3TAU1 = 0.0_JPRB

      ZO3V1   = ZO3V1   + ZFACT15 * ZO3MOL1
      ZFACT1  = ZFACT1  + ZO3V15  * ZO3MOL1
      ZO3MOL1 = 0.0_JPRB

      ZO3F1 (JL,JK) = ZO3F1 (JL,JK) + RMD/RMO3 * ZO3V1
      ZO3V1        = 0.0_JPRB

      ZTAUR1               = ZTAUR1                + ZOTRAY1 (JL,JWL,JK)
      ZOTRAY1 (JL,JWL,JK)  = 0.0_JPRB
      
      ZDP1 (JL,JK) = ZDP1 (JL,JK) + ZSIGMA*RNAVO/(ZMD*RG)/10._JPRB * ZTAUR1
      ZTAUR1       = 0.0_JPRB

      PCO21 (JL,JK) = PCO21 (JL,JK) + RMD/RMCO2 * ZVCO21
      ZVCO21        = 0.0_JPRB

!-- ZFACT is the number of gas molecules in a layer of pressure ZDP at temperature T

      ZTF1(JL,JK) = ZTF1(JL,JK) - ZFACT15/ZTF15(JL,JK) * ZFACT1
      ZDP1(JL,JK)=ZDP1(JL,JK) + ZFACT15/ZDP15(JL,JK)* ZFACT1
      ZFACT1     = 0.0_JPRB

    ENDDO
  ENDDO
ENDDO

! BEN - neglect contribution to gradients from computation of saturation (for now)
!   CALL SATURAD (1 , KLEN , KLEN  , 1 , KLEV,&
!  &              ZPF15, ZTF15 , ZQSAT15,
!  &              ZPF1,  ZTF1 ,  ZQSAT1)

!     INTEGRATE BACKWARD HYDROSTATIC EQUATION TO COMPUTE gradients for ZDZ, T, and q  from ZDP
DO JK=KLEV,1,-1
  DO JL=1,KLEN

     ZDZ1(JL,JK) = ZDZ1(JL,JK) + ZDP1(JL,JK)/ ( 287.0 *ZTF15(JL,JK)* &
     &                         (1.0 +0.608 *ZQF15(JL,JK)) ) * ZPF15(JL,JK)* RG

     ZTF1(JL,JK)  = ZTF1(JL,JK) - ZDP1(JL,JK) * ZDZ15(JL,JK)/ ( 287.0 *ZTF15(JL,JK)**2 * &
     &                         (1.0 +0.608 *ZQF15(JL,JK)) ) * ZPF15(JL,JK)* RG

     ZQF1(JL,JK)  = ZQF1(JL,JK) - ZDP1(JL,JK) * ZDZ15(JL,JK) * 0.608 / ( 287.0 *ZTF15(JL,JK)* &
     &                            (1+ 0.608 *ZQF15(JL,JK)) **2   )   * ZPF15(JL,JK)* RG


     ZDP1(JL,JK) =  0.0_JPRB
  ENDDO
ENDDO

DO JK=KLEV,1,-1
  DO JL=1,KLEN
    ZGEOMH1 (JL,JK-1) = ZGEOMH1 (JL,JK-1) + ZDZ1 (JL,JK) /RG  
    ZGEOMH1 (JL,JK)   = ZGEOMH1 (JL,JK)   - ZDZ1 (JL,JK) /RG 
    ZDZ1 (JL,JK) = 0.0_JPRB
  ENDDO
ENDDO

! Setting to zero adjoint variables

    PGEOP1(:,:) = 0.0_JPRB
    PTQO1(:,:,:)= 0.0_JPRB

    PGEOP1(:,1:KLEV) = PGEOP1(:,1:KLEV) + ZGEOMH1(:,1:KLEV)
    PGEOP1(:,1)      = PGEOP1(:,1) + ZGEOMH1(:,0)
    ZGEOMH1(:,0:KLEV) = 0.0_JPRB

    ZPF1 (1:KLEN,:) = 0.0_JPRB ! no variations of pressure

    PTQO1(1:KLEN,:,3) = PTQO1(1:KLEN,:,3) + ZO3F1(1:KLEN,:)
    ZO3F1(1:KLEN,:) = 0.0_JPRB
    PTQO1(1:KLEN,:,2) = PTQO1(1:KLEN,:,2) + ZQF1 (1:KLEN,:)
    ZQF1 (1:KLEN,:) = 0.0_JPRB
    PTQO1(1:KLEN,:,1) = PTQO1(1:KLEN,:,1) + ZTF1 (1:KLEN,:)
    ZTF1 (1:KLEN,:) = 0.0_JPRB
    
ZAERTMR1(:,:)=0.0_JPRB
DO JL=1,KLEN
  DO JAER=1,INACTAERO
    ZAERTMR1(JL,:) = ZAERTMR1(JL,:) + ZAERF1(JL,:,JAER)* ZFRACAER15(JL,:,JAER)
    ZAERF1(JL,:,JAER)=0.0_JPRB
  ENDDO
ENDDO
! Put back gradient of cost function with respect to total aerosol mixing ratio
DO JL=1,KLEN
  DO JF=1,YDGP5%NGEMS
    IF(YDGP5%GEMS_IGRIB(JF) == 210048) THEN
      PAELIDPROF(JL,1:KLEV,JF) = PAELIDPROF(JL,1:KLEV,JF) + ZAERTMR1(JL,:)
      ZAERTMR1(JL,:)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_LIDSIMAD',1,ZHOOK_HANDLE)
END SUBROUTINE AER_LIDSIMAD
