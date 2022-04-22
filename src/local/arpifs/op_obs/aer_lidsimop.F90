SUBROUTINE AER_LIDSIMOP &
  &( KDLEN, KLEN , KLEV , KDFS, &
  &  PLIDPRES , & 
  &  PCO2 , PNO2, PGEOP, PTQO, PAELIDPROF, YDGP5,  &
  &  PLISIS, PLISIT, PULISIS,PULISIT,PEXTT, PEXTS &
  & )

!**** *AER_LIDSIMOP* - SIMULATES THE LIDAR SIGNAL OUT OF THE ATMOSPHERE
!                    ASSUMED CLOUDLESS BUT INCLUDING THE GEMS-AEROSOLS.

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE TOTAL BACK-SCATTERING COEFFICIENT 
!          ACCOUNTING RAYLEIGH SCATTERING, AEROSOL SCATTERING AND GASEOUS 
!          ABSORPTION FOR THE THREE USUAL WAVELENGTHS (355, 532 AND 1064 NM)
!          FROM THE STANDARD "LIDAR EQUATION".

!**   INTERFACE.
!     ----------

!          *AER_LIDSIMOP* IS CALLED FROM *HOP*

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
! PLIDPRES (KDLEN,KLEV) : LIDAR LEVEL PRESSURE (Pa)
! PCO2   (KDLEN,KLEV)   : CO2 CONCENTRATION (kg/kg)
! PNO2   (KDLEN,KLEV)   : NO2 CONCENTRATION (kg/kg)
! PTQO   (KDLEN,KLEV)   : TEMPERATURE (K), SPECIFIC HUMIDITY (kg/kg) and OZONE
! PGEOP  (KDLEN,KLEV)   : GEOPOTENTIAL 

!     ==== OUTPUTS ===
! PLISIS (KDLEN,NWLID,0:KLEV) : BACK-SCATTERING SIGNAL FOR SURFACE LIDAR
! PLISIT (KDLEN,NWLID,0:KLEV) : BACK-SCATTERING SIGNAL FOR SATELLITE LIDAR
! PULISIS (KDLEN,NWLID,0:KLEV) : UNATTENUATED BACK-SCATTERING SIGNAL FOR SURFACE LIDAR
! PULISIT (KDLEN,NWLID,0:KLEV) : UNATTENUATED BACK-SCATTERING SIGNAL FOR SATELLITE LIDAR
! PEXTT (KDLEN,NWLID,0:KLEV) : EXTINCTION FOR SATELLITE LIDAR
! PEXTS (KDLEN,NWLID,0:KLEV) : EXTINCTION FOR SURFACE LIDAR

!     METHOD.
!     -------
!         from Ackerman, J., 1998, J.Atmos.Ocean.Technol., 15, 1043-1050
!         and Huneeus, N. and O. Boucher, 2007, J.Geophys.Res. 

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
! 
!     Angela Benedetti: 24-Aug-2010 Modifications to call it from hop.F90
!     A. Geer           29-Dec-2015 OOPS cleaning (GOM_PLUS added; GEMSPROF removed)
!        S. Remy     16-Sep-2016    Add nitrate and ammonium
!     A Benedetti :  05-Mar-2018    Add unattenuated backscatter
!     ------------------------------------------------------------------

USE PARKIND1      , ONLY : JPIM, JPRB
USE YOMHOOK       , ONLY : LHOOK, DR_HOOK
USE YOEAERLID     , ONLY : YREAERLID
USE YOEAEROP      , ONLY : ALF_SU, OMG_SU, ALF_OM, OMG_OM, ALF_DD, OMG_DD, ALF_SS, OMG_SS, &
 &                         ALF_BC, OMG_BC,  RALI_BC,RALI_DD,RALI_OM,RALI_SU,RALI_SS, ALF_NI, OMG_NI, &
 &                         RALI_NI, ALF_AM, OMG_AM, RALI_AM, ALF_SOA, OMG_SOA,RALI_SOA
USE YOEAERSNK     , ONLY : YREAERSNK
USE YOEAERATM     , ONLY : YREAERATM
USE YOEAERVOL     , ONLY : YREAERVOL
USE YOMCST        , ONLY : RG, RMD, RMCO2, RMNO2, RMO3, RNAVO, RPI, YRCST
USE YOETHF        , ONLY : YRTHF
USE GOM_PLUS      , ONLY : TYPE_GOM_PLUS
USE YOM_GRIB_CODES, ONLY : NGRBAERMR01,NGRBAERMR02, NGRBAERMR03, NGRBAERMR04, NGRBAERMR05, &
 & NGRBAERMR06, NGRBAERMR07, NGRBAERMR08, NGRBAERMR09, NGRBAERMR10, NGRBAERMR11, NGRBAERMR12, &
 & NGRBAERMR13, NGRBAERMR14, NGRBAERMR15, NGRBAERMR16, NGRBAERMR17, NGRBAERMR18, &
 & NGRBAERMR19, NGRBAERMR20, NGRBAERMR21, NGRBAERLG
USE YOEPHLI       , ONLY : YREPHLI

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KDLEN  ! Total number of observations
INTEGER(KIND=JPIM),INTENT(IN) :: KLEN   ! Current number of observations
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV   ! Number of lidar vertical levels <= YREAERLID%NAELID
INTEGER(KIND=JPIM),INTENT(IN) :: KDFS   ! Dimension of T, q, O3 array 

TYPE(TYPE_GOM_PLUS), INTENT(IN) :: YDGP5 ! Model variables at observation locations
REAL(KIND=JPRB)   ,INTENT(IN) :: PLIDPRES(KDLEN,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN) :: PCO2(KDLEN,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN) :: PNO2(KDLEN,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN) :: PGEOP(KDLEN,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN) :: PTQO(KDLEN,KLEV,KDFS) ! Interpolated profiles of temperature, specific
                                                        ! humidity and ozone
REAL(KIND=JPRB)   ,INTENT(IN) :: PAELIDPROF(KDLEN,KLEV,YDGP5%NGEMS) ! Interpolated composition profiles

REAL(KIND=JPRB)   ,INTENT(OUT) :: PLISIS(KDLEN,YREAERLID%NWLID,0:KLEV), PLISIT(KDLEN,YREAERLID%NWLID,0:KLEV), & 
     &                            PULISIS(KDLEN,YREAERLID%NWLID,0:KLEV), PULISIT(KDLEN,YREAERLID%NWLID,0:KLEV), &
     &                            PEXTS(KDLEN,YREAERLID%NWLID,0:KLEV), PEXTT(KDLEN,YREAERLID%NWLID,0:KLEV)

!     ------------------------------------------------------------------

!*       0.1   LOCAL ARRAYS
!              ------------

INTEGER(KIND=JPIM) :: IAER, IBIN, IK, ITYP, IWAVL, IFLAG, IEFRH
INTEGER(KIND=JPIM), PARAMETER  :: JTABMX=12
INTEGER(KIND=JPIM) :: JAER, JL, JK, JP, JTAB, JWL, JF, JLEN
INTEGER(KIND=JPIM) :: IRH(KLEN,KLEV)
INTEGER(KIND=JPIM) :: INACTAERO

LOGICAL :: LLPRINT

REAL(KIND=JPRB) :: ZAERF(KDLEN,KLEV,YDGP5%NGEMS)
REAL(KIND=JPRB) :: ZTF(KDLEN, KLEV), ZQF(KDLEN,KLEV),  ZO3F(KDLEN,KLEV), ZPF(KDLEN,KLEV)
REAL(KIND=JPRB), DIMENSION(KLEN,KLEV)   :: ZQSAT
REAL(KIND=JPRB) :: ZRHCL
REAL(KIND=JPRB) :: ZDP(KDLEN, KLEV), ZDZ(KDLEN,KLEV)
REAL(KIND=JPRB) :: ZGEOMH(KDLEN,0:KLEV)

REAL(KIND=JPRB) :: ZOTRAY(KDLEN,YREAERLID%NWLID,KLEV), ZOTCO2(KDLEN,YREAERLID%NWLID,KLEV), ZOTO2(KDLEN,YREAERLID%NWLID,KLEV)
REAL(KIND=JPRB) :: ZOTNO2(KDLEN,YREAERLID%NWLID,KLEV), ZOTO3(KDLEN,YREAERLID%NWLID,KLEV) , ZPATHT(KDLEN,YREAERLID%NWLID,0:KLEV)
REAL(KIND=JPRB) :: ZPATHS(KDLEN,YREAERLID%NWLID,0:KLEV)

REAL(KIND=JPRB) :: ZAERMSS(KDLEN,KLEV,YDGP5%NGEMS), ZALF(YDGP5%NGEMS), ZOMG(YDGP5%NGEMS), ZLIR(YDGP5%NGEMS)
REAL(KIND=JPRB) :: ZAERTMR(KDLEN,KLEV)  ! Total aerosol mixing ratio (control variable) - not needed
REAL(KIND=JPRB) :: ZBSC(KDLEN,YREAERLID%NWLID,KLEV), ZAEROD(KDLEN,YREAERLID%NWLID,KLEV) 
REAL(KIND=JPRB) :: ZBSCMOL(KDLEN,YREAERLID%NWLID,KLEV), ZBOLTZ

REAL(KIND=JPRB) :: ZA1, ZA2, ZA3, ZA4, ZAA, ZAN, ZAN2, ZDPG, ZEPSAER, ZFAC, ZFACT, ZMD, ZLIDRAT, &
  & ZNO2MOL, ZNO2TAU, ZNO2V, ZO3MOL, ZO3TAU, ZO3V, ZREF, ZSIGMA, &
  & ZTAUR, ZTE, ZTE2, &
  & ZTRSOT, ZTRTOT, ZUSOT, ZUTOT, ZWLUM, ZWLCM, ZWL2CM, ZWL4CM, ZWN, ZWN2

REAL(KIND=JPRB) :: ZGCO2, ZMCO2, ZPHIC, ZPSIC, ZUCO2, ZUPCO2, ZVCO2
REAL(KIND=JPRB) :: ZO2 , ZGO2 , ZMO2 , ZPHIO, ZPSIO, ZUO2 , ZUPO2

REAL(KIND=JPRB) :: ZCAER(KDLEN), ZCRAY(KDLEN), ZCOZO(KDLEN), ZCOXY(KDLEN), ZCCO2(KDLEN), ZCNO2(KDLEN)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL :: LLPHYLIN 

#include "satur.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_LIDSIMOP',0,ZHOOK_HANDLE)
ASSOCIATE( &
 & RSS_RH80_MASSFAC=>YREAERATM%RSS_RH80_MASSFAC, &
 & JWLID=>YREAERLID%JWLID, NAELID=>YREAERLID%NAELID, NWLID=>YREAERLID%NWLID, &
 & RLICLS=>YREAERLID%RLICLS, RLICO2=>YREAERLID%RLICO2, &
 & RLIDELT=>YREAERLID%RLIDELT, RLINO2=>YREAERLID%RLINO2, RLINS=>YREAERLID%RLINS, &
 & RLIO2=>YREAERLID%RLIO2, RLIO3=>YREAERLID%RLIO3, RLIPREF=>YREAERLID%RLIPREF, &
 & RLIT0=>YREAERLID%RLIT0, RLITREF=>YREAERLID%RLITREF, RWLID=>YREAERLID%RWLID, &
 & RRHTAB=>YREAERSNK%RRHTAB, &
 & NVOLOPTP=>YREAERVOL%NVOLOPTP, &
 & YAERO_DESC=>YREAERATM%YAERO_DESC)


LLPHYLIN = YREPHLI%LPHYLIN

!*         0.     SPECTRAL INFORMATION OF RELEVANCE
!                 ---------------------------------

! Rayleigh |  355nm / 28169cm-1 |  532nm / 18797cm-1 | 1064nm / 9398.5cm-1 |
! Index    |          2         |          8         |         16          ! 
! CO2      |         no         |         no         |        yes          |
! O2       |         no         |         no         |        yes          |
! O3       |        yes         |        yes         |        yes          |
! NO2      |        yes         |        yes         |      yes          |

!*    --------------------------------------------------------------------------

!*         1.     INITIALIZATION AND PREPARATORY WORK
!                 -----------------------------------

ZREF=RLITREF/RLIPREF   
ZO2=RMO3*2._JPRB/3._JPRB
ZEPSAER=1.E-15_JPRB
LLPRINT=.FALSE.
!-- the effective relative hunidity is the low value (20%) assumed for
!hydrophobic component of OM
IEFRH=3


! Extract aerosol mixing ratios from GEMS profiles structure
DO JLEN=1,KLEN
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
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR16 .OR. &
     & YDGP5%GEMS_IGRIB(JF) == NGRBAERMR17) JP=6 ! NO3
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR18) JP=7 ! NH4
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR19 .OR. &
     & YDGP5%GEMS_IGRIB(JF) == NGRBAERMR20 .OR. &
     & YDGP5%GEMS_IGRIB(JF) == NGRBAERMR21) JP=8 ! SOA
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR13) JP=9 ! Volcanic ashes 
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR14 .OR. &
     & YDGP5%GEMS_IGRIB(JF) == NGRBAERMR15)  JP=10 ! Volcanic SO4 and SO2


    IF(JP > 0) THEN
      IAER=IAER+1
      ZAERF(JLEN,:,IAER)=PAELIDPROF(JLEN,1:KLEV,JF)
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERLG) THEN ! Total aerosol mixing ratio
      ZAERTMR(JLEN,:)=PAELIDPROF(JLEN,1:KLEV,JF)
    ENDIF
  ENDDO
ENDDO
INACTAERO=IAER

! Save T, Q, and O3 profiles in separate arrays
ZTF (1:KLEN,:) = PTQO(1:KLEN,:,1)
ZQF (1:KLEN,:) = PTQO(1:KLEN,:,2)
ZO3F(1:KLEN,:) = PTQO(1:KLEN,:,3)

ZPF(1:KLEN,:)= PLIDPRES(1:KLEN,:)

! COMPUTE THICKNESS

ZGEOMH(:,0) = PGEOP(:,1) 
ZGEOMH(:,1:KLEV) = PGEOP(:,1:KLEV)

DO JL=1,KLEN
 DO JK=1,KLEV
    ZDZ(JL,JK) = (ZGEOMH(JL,JK-1)-ZGEOMH(JL,JK)) /RG
 ENDDO
    ZDZ(JL,1) = 60.0_JPRB ! set top layer equal to 60 meters
ENDDO


!     INTEGRATE HYDROSTATIC EQUATION TO COMPUTE ZDP from ZDZ

DO JK=1,KLEV
  DO JL=1,KLEN
     ZDP(JL,JK) = RG* ZDZ(JL,JK) / ( 287.0_JPRB*ZTF(JL,JK)* &
     & (1.0_JPRB+0.608_JPRB*ZQF(JL,JK)) ) * ZPF(JL,JK) 
  ENDDO
ENDDO

IFLAG=2
CALL SATUR (YRTHF, YRCST, 1 , KLEN , KLEN  , 1 , KLEV, LLPHYLIN, &
  & ZPF, ZTF , ZQSAT , IFLAG)

!-- define RH index from "clear-sky" relative humidity
IRH(:,:)=1
ZRHCL=0.0_JPRB
DO JK=1,KLEV
  DO JL=1,KLEN
    ZRHCL=(ZQF(JL,JK)/ZQSAT(JL,JK))*100._JPRB
    DO JTAB=1,JTABMX
      IF (ZRHCL > RRHTAB(JTAB)) THEN
        IRH(JL,JK)=JTAB
      ENDIF
    ENDDO
  ENDDO
ENDDO

ZBSC  (1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZAEROD(1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZBSCMOL(1:KLEN,1:NWLID,1:KLEV)=0._JPRB

ZOTRAY(1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZOTNO2(1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZOTO3 (1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZOTO2 (1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZOTCO2(1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZPATHT(1:KLEN,1:NWLID,0:KLEV)=0._JPRB
ZPATHS(1:KLEN,1:NWLID,0:KLEV)=0._JPRB

!*    --------------------------------------------------------------------------

!*         2.      RAYLEIGH AND GASEOUS OPTICAL THICKNESSES LAYER-BY-LAYER
!                  -------------------------------------------------------

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
  ZTAUR =ZSIGMA*101325._JPRB*RNAVO/(RMD*RG)/10._JPRB

  ZOTRAY(1:KLEN,JWL,1:KLEV)=0._JPRB
  ZOTNO2(1:KLEN,JWL,1:KLEV)=0._JPRB
  ZOTO3 (1:KLEN,JWL,1:KLEV)=0._JPRB

  DO JK=1,KLEV
    DO JL=1,KLEN

!-- ZFACT is the number of gas molecules in a layer of pressure ZDP at temperature T
      ZFACT=RNAVO*(ZDP(JL,JK)*RLITREF)/(ZTF(JL,JK)*RLIPREF)

!-- Rayleigh: Bodhaine et al., 1999: JAOT, 16, 1854-1861.   tau= sigma*DelP*Avo / (mmw_air * g)
      ZVCO2=PCO2(JL,JK)*RMD/RMCO2
      ZMD = 28.9595_JPRB + 15.0556_JPRB * ZVCO2
      ZMD = RMD
      ZTAUR=ZSIGMA*ZDP(JL,JK)*RNAVO/(ZMD*RG)/10._JPRB
      ZOTRAY(JL,JWL,JK) = ZTAUR

!-- O3 
      ZO3V=ZO3F(JL,JK)*RMD/RMO3                             ! volume mixing ratio in layer
      ZO3MOL=ZO3V*ZFACT                                     ! number of O3 molecules in layer (pressure unit: kg m-1 s-2)
      ZO3TAU=1.E-04_JPRB*RLIO3(3,JWL)*ZO3MOL/(RMO3*RG)    ! cross-section from cm2 molec-1 to m-2 molec-1
      ZOTO3(JL,JWL,JK) = ZO3TAU

!-- NO2
      ZNO2V=PNO2(JL,JK)*RMD/RMNO2                           ! volume mixing ratio in layer
      ZNO2MOL=ZNO2V*ZFACT                                   ! number of NO2 molecules in layer (pressure unit: kg m-1 s-2)
      ZNO2TAU=1.E-04_JPRB*RLINO2(3,JWL)*ZNO2MOL/(RMNO2*RG)  ! cross-section from cm2 molec-1 to m-2 molec-1
      ZOTNO2(JL,JWL,JK) = ZNO2TAU
    ENDDO
  ENDDO
ENDDO


!-- to start with (and to be replaced by more proper quantities once debugged!)
ZMCO2= 357.E-06_JPRB*RMCO2/RMD
ZMO2 = 0.20947_JPRB *ZO2  /RMD
   
ZOTCO2(1:KLEN,1:NWLID,1:KLEV)=0._JPRB
ZOTO2 (1:KLEN,1:NWLID,1:KLEV)=0._JPRB
DO JWL=1,NWLID
DO JK=1,KLEV
  DO JL=1,KLEN
    ZTE=ZTF(JL,JK)-RLIT0
    ZTE2=ZTE*ZTE
    ZDPG=ZDP(JL,JK)/RG

!-- CO2 -------------------------------
    ZPHIC=EXP(RLICO2(3)*ZTE+RLICO2(4)*ZTE2)
    ZPSIC=EXP(RLICO2(5)*ZTE+RLICO2(6)*ZTE2)    
    ZUCO2 =ZPHIC*ZMCO2*ZDPG
    ZUPCO2=ZPSIC*ZMCO2*ZDPG*ZPF(JL,JK)/RLIPREF

!-- Goody model for CO2
    ZGCO2=RLICO2(1)*ZUCO2/SQRT(1._JPRB+RLICO2(1)*ZUCO2*ZUCO2/(RLICO2(2)*ZUPCO2))

    ZOTCO2(JL,JWL,JK) = ZGCO2

!-- O2 --------------------------------
    ZPHIO=EXP(RLIO2(3) *ZTE+RLIO2(4) *ZTE2)
    ZPSIO=EXP(RLIO2(5) *ZTE+RLIO2(6) *ZTE2)
    ZUO2 =ZPHIO*ZMO2*ZDPG
    ZUPO2=ZPSIO*ZMO2*ZDPG*ZPF(JL,JK)/RLIPREF
!-- Goody model for O2
    ZGO2=RLIO2(1)*ZUO2/SQRT(1._JPRB+RLIO2(1)*ZUO2*ZUO2/(RLIO2(2)*ZUPO2))

    ZOTO2(JL,JWL,JK) = ZGO2
  ENDDO
ENDDO
ENDDO

!*    --------------------------------------------------------------------------

!*         3.      AEROSOL-RELATED QUANTITIES
!                  --------------------------

ZAERMSS(:,:,:) = 0.0_JPRB
DO JAER=1,INACTAERO
  DO JK=1,KLEV
    DO JL=1,KLEN
      IF (ZAERF(JL,JK,JAER) > ZEPSAER) THEN
        ZAERMSS(JL,JK,JAER)= ZAERF(JL,JK,JAER)*ZDP(JL,JK)/RG
      ENDIF
    ENDDO
  ENDDO
ENDDO


!-- path from TOA
ZPATHS(1:KLEN,1:NWLID,0:KLEV)=0._JPRB
ZPATHT(1:KLEN,1:NWLID,0:KLEV)=0._JPRB

PLISIT(:,:,:)=0._JPRB
PLISIS(:,:,:)=0._JPRB
PULISIT(:,:,:)=0._JPRB
PULISIS(:,:,:)=0._JPRB
PEXTT(:,:,:)=0._JPRB
PEXTS(:,:,:)=0._JPRB

DO JWL=1,NWLID
  IWAVL=JWLID(JWL)

  DO JAER=1,INACTAERO
    ITYP=YAERO_DESC(JAER)%NTYP
    IBIN=YAERO_DESC(JAER)%NBIN

!-- ITYP is the aerosol type 1:SS,   2:DD,   3:OM    4:BC,   5:SU,   6:FA,   7:BS,   8=VS,
!   IBIN is the bin index: 1-3:SS, 1-3:DD,   2:OM,   2:BC,   2:SU,   1:FA,   1:SB    1=VS,
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

!-- total aerosol optical depth within a layer (sum on aerosol contribution)
        ZAEROD(JL,JWL,JK) = ZAEROD(JL,JWL,JK) + ZAERMSS(JL,JK,JAER) * ZFAC * 1.E+03_JPRB * ZALF(JAER)

!-- back-scatter coefficient computed from:
!    1/ the extinction coefficient multiplied by the lidar ratio  (= backscatter / extinction)
!    2/ this extinction coefficient being a layer integral will have to be divided by the geometrical thickness
        ZBSC(JL,JWL,JK) = ZBSC(JL,JWL,JK) + ZAERMSS(JL,JK,JAER) * ZFAC * 1.E+03_JPRB * ZALF(JAER) * ZLIDRAT
      ENDDO
    ENDDO
  ENDDO

ZBOLTZ = 1.38E-23_JPRB ! (J/K, Boltzman constant)
!DO JWL=1,NWLID
! Calculation of molecular backscatter
  ZWLUM=RWLID(JWL)*1.0E+06_JPRB ! wavelenght in micrometers
    DO JK=1,KLEV
      DO JL=1,KLEN
! Molecular backscatter from Collis and Russel (1976)
        ZBSCMOL(JL,JWL,JK) = 5.45E-32_JPRB * ZPF(JL,JK)/(ZBOLTZ * ZTF(JL,JK)) * (ZWLUM/0.55_JPRB)**(-4.09_JPRB)
      ENDDO
    ENDDO
!*    --------------------------------------------------------------------------

!*         4.      COMPUTATION OF LIDAR SIGNAL
!                  ---------------------------


  ZPATHS(1:KLEN,JWL,KLEV)=0._JPRB
  ZPATHT(1:KLEN,JWL,0   )=0._JPRB


!*         4.1     TOP-DOWN (SATELLITE VIEW)
!                  -------------------------

  ZCAER(:)=0._JPRB
  ZCRAY(:)=0._JPRB
  ZCOZO(:)=0._JPRB
  ZCNO2(:)=0._JPRB
  ZCOXY(:)=0._JPRB
  ZCCO2(:)=0._JPRB
  DO JK=1,KLEV
    IK=KLEV-JK
    DO JL=1,KLEN
      ZCAER(JL)=ZCAER(JL)+ZAEROD(JL,JWL,JK)
      ZCRAY(JL)=ZCRAY(JL)+ZOTRAY(JL,JWL,JK)
      ZCOZO(JL)=ZCOZO(JL)+ZOTO3 (JL,JWL,JK)
      ZCNO2(JL)=ZCNO2(JL)+ZOTNO2(JL,JWL,JK)
      ZCOXY(JL)=ZCOXY(JL)+ZOTO2 (JL,JWL,JK)
      ZCCO2(JL)=ZCCO2(JL)+ZOTCO2(JL,JWL,JK)

!-- total optical depth (aerosols + gases) from TOA to level JK
      ZPATHT(JL,JWL,JK) = ZPATHT(JL,JWL,JK-1) + 2._JPRB * ( ZAEROD(JL,JWL,JK) &
      & +ZOTRAY(JL,JWL,JK)+ZOTO3(JL,JWL,JK)+ZOTNO2(JL,JWL,JK)+ZOTO2(JL,JWL,JK)+ZOTCO2(JL,JWL,JK) )
      ZUTOT=MIN(200._JPRB,ZPATHT(JL,JWL,JK))
!-- attenuation seen from TOA to level JK
      ZTRTOT=EXP(-ZUTOT)

!-- lidar signal seen from TOA
      PLISIT(JL,JWL,JK)=RLICLS* (ZBSC(JL,JWL,JK)/ZDZ(JL,JK) + ZBSCMOL(JL,JWL,JK) ) * ZTRTOT
!-- unattenuated lidar signal seen from TOA  (without molecular backscatter)
      PULISIT(JL,JWL,JK)=RLICLS* ZBSC(JL,JWL,JK)/ZDZ(JL,JK) 
!--extinction from TOA
      PEXTT(JL,JWL,JK)= ZTRTOT

!*         4.2     BOTTOM-UP (GROUND STATION VIEW)
!                  -------------------------------

!-- lidar signal seen from the surface
      ZPATHS(JL,JWL,IK) = ZPATHS(JL,JWL,IK+1) + 2._JPRB * (ZAEROD(JL,JWL,IK+1) &
      & +ZOTRAY(JL,JWL,IK+1)+ZOTO3(JL,JWL,IK+1)+ZOTNO2(JL,JWL,IK+1)+ZOTO2(JL,JWL,IK+1)+ZOTCO2(JL,JWL,IK+1))
      ZUSOT=MIN(200._JPRB,ZPATHS(JL,JWL,IK))

!-- attenuation seen from surface to level JK
      ZTRSOT=EXP(-ZUSOT)
      PLISIS(JL,JWL,IK)=RLICLS*(ZBSC(JL,JWL,IK+1)/ZDZ(JL,IK+1) + ZBSCMOL(JL,JWL,IK+1) )*ZTRSOT 
!-- unattenuated lidar signal from ground-based lidar  (without molecular backscatter)
      PULISIS(JL,JWL,IK)=RLICLS*ZBSC(JL,JWL,IK+1)/ZDZ(JL,IK+1) 
!--extinction from ground-based lidar
      PEXTS(JL,JWL,IK)= ZTRSOT


    ENDDO
  ENDDO
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_LIDSIMOP',1,ZHOOK_HANDLE)
END SUBROUTINE AER_LIDSIMOP
