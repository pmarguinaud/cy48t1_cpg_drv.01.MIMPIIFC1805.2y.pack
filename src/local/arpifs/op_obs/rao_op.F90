SUBROUTINE RAO_OP&
 & (KLEVG, KDLEN  , KLEN  , KMXCOUNT, &
 & YDGP5,&
 & PAEROD )

!*** *RAO_OP* Operator routine for aerosol optical depth observations.

!**   INTERFACE.
!     ----------

!     WHERE KDLEN      = First dimension: Max number of observations, this set (INPUT)
!           KLEN       = Number of aerosol optical depth obs (INPUT)
!           IMXCOUNT   = Number of channels (INPUT)
!           ZQF       =  Specific humidity at observation points, model levels (INPUT)
!           ZTF       =  Temperature at observation points, model levels (INPUT)
!           ZPRESH,ZPRESF =   Half/Full level pressure values at obs points, model levels (INPUT)

!           ZXPP       = Aerosol optical depth at different wavelengths (OUTPUT)



!**   AUTHOR.
!     -------
!        ANGELA BENEDETTI        *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2005-06-17 (adapted from aer_bdgmass.F90 by J.-J. Morcrette)
!        A. Benedetti: 2009-02-25. Readaptation of code to allow for two sulphate variables
!                                   with SO2 not contributing to the total aerosol optical depth
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     A. Geer     29-Dec-2015    OOPS cleaning: remove GEMSPROF and GFL references
!        S. Remy     16-Sep-2016    Add nitrate and ammonium


!     PURPOSE.
!     --------
!     - derive aerosol optical thickness from profiles of mass mixing ratio 

! Optical properties in su_aerop.F90 avaliable at following wavelenghts
! (340)355(380)400 440 469 500 532 550 645 670 800 858 865(1020)1064 1240 1640 2130
! corresponding IWAVL 
!   1   2   3   4   5   6   7   8   9  10   11  12  13  14  15  16   17    18   19
!                      ---         ---     ---         ---           ---- ---- ---- MODIS
!              ***                             ***                                  SEVIRI
!  ---     ---     xxx     ---             xxx         xxx xxxx                     AERONET
!      000                     000                               000                EARLINET/CALIPSO
! 

!     'MODIS Channels' (channels from other instruments can be added) 
!      -         470 nm, 550 nm, 660 nm, 870 nm, 1240 nm, 1630 nm, 2130 nm
! Corresponding IWAVL
!      IWAVL = [ 6, 9, 11, 14, 17, 18, 19] 
! **** DEPENDS ON J.-J.'s choice of numbering in su_aerop.F90 !!!!!! CHECK PERIODICALLY !!!


!-----------------------------------------------------------------------

USE PARKIND1      , ONLY : JPIM, JPRB
USE YOMHOOK       , ONLY : LHOOK, DR_HOOK
USE YOMDB
USE YOMCST        , ONLY : RG, YRCST
USE YOETHF        , ONLY : YRTHF
USE YOEAEROP      , ONLY : ALF_SU, ALF_OM, ALF_DD, ALF_SS, ALF_BC, ALF_NI, ALF_AM, ALF_SOA
USE YOEAERATM     , ONLY : YREAERATM
USE YOEAERSNK     , ONLY : YREAERSNK
USE GOM_PLUS      , ONLY : TYPE_GOM_PLUS, IH
USE YOM_GRIB_CODES, ONLY : NGRBAERMR01,NGRBAERMR04, &
 & NGRBAERMR07, NGRBAERMR08, NGRBAERMR09, NGRBAERMR10, NGRBAERMR11, &
 & NGRBAERMR16, NGRBAERMR18, NGRBAERMR19, NGRBAERMR20, NGRBAERSM
USE YOEPHLI       , ONLY : YREPHLI

IMPLICIT NONE


!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVG, KDLEN, KLEN, KMXCOUNT

TYPE (TYPE_GOM_PLUS),INTENT(IN)  :: YDGP5

REAL(KIND=JPRB)   ,INTENT(OUT) :: PAEROD(KDLEN,KMXCOUNT) ! AEROSOL OPT DEPTH (FUNCTION OF HORIZONTAL POSITION AND CHANNEL #)



!*       0.5   LOCAL VARIABLES
!              ---------------
! INTEGER(KIND=JPIM), PARAMETER :: IWAVLMX=7 ! Max number of MODIS channels
!INTEGER(KIND=JPIM), DIMENSION(IWAVLMX) :: IWAVL=(/6, 9, 11, 14, 17, 18, 19/) 
                                         ! Index for MODIS channels - based on latest su_aerop.F90
INTEGER(KIND=JPIM) :: IAER, ITYP, IFLAG, IWAVLINDEX
INTEGER(KIND=JPIM) :: JAER, JK, JL, JJ, JF, JP, JLEN, JTAB
INTEGER(KIND=JPIM), PARAMETER  :: JTABMX=12
INTEGER(KIND=JPIM), ALLOCATABLE :: ITYPE(:) !index to determine optically active aerosol species
INTEGER(KIND=JPIM), DIMENSION(KLEN,KLEVG) :: IRH

INTEGER(KIND=JPIM) :: INACTAERO

REAL(KIND=JPRB), DIMENSION(KLEN,KLEVG)   :: ZQSAT, ZT, ZPF
REAL(KIND=JPRB), ALLOCATABLE   :: ZAERTAU(:,:)
REAL(KIND=JPRB), ALLOCATABLE   :: ZAERF(:,:,:) ! aer mixing ratios from traj (interpolated to obs locations)
REAL(KIND=JPRB), DIMENSION(KLEN,KLEVG) :: ZAERFMR ! fine aer mixing ratios from traj
REAL(KIND=JPRB), ALLOCATABLE :: ZALF(:)
REAL(KIND=JPRB)    :: ZRHCL(KLEN,KLEVG), ZAEREXT, ZDP, ZFAC

REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL :: LLPHYLIN 

#include "satur.intfb.h"
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RAO_OP',0,ZHOOK_HANDLE)
ASSOCIATE(&
 & RSS_RH80_MASSFAC=>YREAERATM%RSS_RH80_MASSFAC, &
 & RRHTAB=>YREAERSNK%RRHTAB, &
 & NMAXTAER=>YREAERATM%NMAXTAER, NTYPAER=>YREAERATM%NTYPAER)

LLPHYLIN = YREPHLI%LPHYLIN

! Remember to put in a check for maximum IBIN allowed per aerosol type


ALLOCATE(ZAERF(KLEN,KLEVG,YDGP5%NGEMS))
ALLOCATE(ITYPE(YDGP5%NGEMS))
ALLOCATE(ZAERTAU(KLEN,YDGP5%NGEMS))
ALLOCATE(ZALF(YDGP5%NGEMS))

ZAERF(:,:,:)=0.0
ITYPE(:)=0
! Extract aerosol mixing ratios from GEMS profiles structure
DO JLEN=1,KLEN
  IAER=0
  DO JF=1,YDGP5%NGEMS
    JP=0
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR01) JP=1 !SS only first bin
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR04) JP=2 !DD only first bin
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR07 .OR. &
     &   YDGP5%GEMS_IGRIB(JF) == NGRBAERMR08) JP=3 ! ORGANIC MATTER - 2 bins
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR09 .OR. &
     &   YDGP5%GEMS_IGRIB(JF) == NGRBAERMR10) JP=4 ! Black Carbon - 2 bins
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR11)  JP=5 ! SO4 
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR16) JP=6 ! NO3
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR18 .OR. &
     &   YDGP5%GEMS_IGRIB(JF) == NGRBAERMR19 .OR. &
     &   YDGP5%GEMS_IGRIB(JF) == NGRBAERMR20) JP=8 ! SOA - 2 bins
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR18) JP=7 ! NH4

    IF(JP > 0) THEN
      IAER=IAER+1
      ZAERF(JLEN,:,IAER)=YDGP5%GEMS(JLEN,1:KLEVG,IH,JF)
      ITYPE(IAER)=JP
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERSM) THEN ! Fine mode aerosol mixing ratio
      ZAERFMR(JLEN,:)=YDGP5%GEMS(JLEN,1:KLEVG,IH,JF)
    ENDIF
  ENDDO
ENDDO
INACTAERO=IAER

! initialise optical depth 
PAEROD(:,:) = 0.0_JPRB 

! Compute relative humidity at observation points, on model levels

ZT(1:KLEN,:) = YDGP5%TF(1:KLEN,1:KLEVG,IH)
ZPF(1:KLEN,:)= YDGP5%PRESF(1:KLEN,:,IH)
IFLAG=2
CALL SATUR (YRTHF, YRCST, 1 , KLEN , KLEN  , 1 , KLEVG, LLPHYLIN, &
  & ZPF, ZT , ZQSAT , IFLAG)

!-- define RH index from "clear-sky" relative humidity
IRH(:,:)=1
DO JK=1,KLEVG
  DO JL=1,KLEN
    ZRHCL(JL,JK)=(YDGP5%QF(JL,JK,IH)/ZQSAT(JL,JK))*100._JPRB
    DO JTAB=1,JTABMX
      IF (ZRHCL(JL,JK) > RRHTAB(JTAB)) THEN
        IRH(JL,JK)=JTAB
      ENDIF
    ENDDO
  ENDDO
ENDDO
  
IF(KMXCOUNT ==1) IWAVLINDEX=9   ! 550nm

DO JJ=1,KMXCOUNT

!   IWAVLINDEX is the wavelength index: for now the choice of channels is tailored to MODIS optical
!         depth retrievals (Remer et al. 2005  and Benedetti et al, 2009, ECMWF Tech. Memo 571)
! *** to be activated, for now assume AOD is only defined at 550nm , i.e. IWAVLINDEX=9 always
!     IWAVLINDEX=IWAVL(JJ)

  DO JL=1,KLEN
   DO JAER=1,INACTAERO 

!-- ITYP is the aerosol type 1:SeaSalt,   2:DesertDust,   3:OrganicMatter, 4: BlackCarbon,
!                            5:Sulphates, 6:FlyingAshes
!   IBIN is the bin index: 1-3:SS, 1-3:DD, 1:OM, 1:BC, 2:SU, 1:FA 
!   Only one "bin" is used for sulphates (SO4) as SO2 does not contribute to the aerosol optical depth. 
! This indeces are defined in su_aerw.F90 

      ITYP=ITYPE(JAER)

      ZAERTAU(JL,JAER)=0.0_JPRB
        
      DO JK=1,KLEVG

! Always double check dimensioning of arrays ALF_** in yoeaerop.F90 !!!!!!!!!
! It has changed a number of times!!!!!!!!
         ZFAC=1.0_JPRB
         IF (ITYP == 1) THEN
          ZALF(JAER)=ALF_SS(IRH(JL,JK),IWAVLINDEX,1) !only one bin of SS contributes
          ZFAC=RSS_RH80_MASSFAC
         ELSEIF (ITYP == 2) THEN                     
          ZALF(JAER)=ALF_DD(1,IWAVLINDEX)            !only one bin of DD contributes
         ELSEIF (ITYP == 3) THEN
          ZALF(JAER)=ALF_OM(IRH(JL,JK),IWAVLINDEX)
         ELSEIF (ITYP == 4) THEN
          ZALF(JAER)=ALF_BC(IWAVLINDEX)
          ELSEIF (ITYP == 5 .OR. ITYP == 10) THEN
          ZALF(JAER)=ALF_SU(IRH(JL,JK),IWAVLINDEX)
         ELSEIF (ITYP == 6) THEN ! Nitrate
            ZALF(JAER)=ALF_NI(IRH(JL,JK),IWAVLINDEX,1)
         ELSEIF (ITYP == 7) THEN ! Ammonium
            ZALF(JAER)=ALF_AM(IRH(JL,JK),IWAVLINDEX)
         ELSEIF (ITYP == 8) THEN ! SOA
            ZALF(JAER)=ALF_SOA(IRH(JL,JK),IWAVLINDEX,1)
         ELSEIF (ITYP == 9) THEN ! volcanic ashes
!-- use dust for 0.9-20 um bin
            ZALF(JAER)=ALF_DD(3,IWAVLINDEX)
         ENDIF
               
! convert units from m^2/g to m^2/kg
            ZALF(JAER)=ZALF(JAER)*1000._JPRB ! conversion from m^2/g to m^2/kg

            ZDP=YDGP5%PRESH(JL,JK,IH)-YDGP5%PRESH(JL,JK-1,IH)
            ZAEREXT=ZAERF(JL,JK,JAER)*ZFAC*ZALF(JAER)

            ZAERTAU(JL,JAER)=ZAERTAU(JL,JAER)+ ZAEREXT*(ZDP/RG)

      ENDDO ! loop over model vertical levels

! Sum over JAER to get total optical depth as function of horizontal position and channel
      PAEROD(JL,JJ)=PAEROD(JL,JJ)+ZAERTAU(JL,JAER)
    ENDDO ! loop over aerosol type

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

DEALLOCATE(ZAERF)
DEALLOCATE(ITYPE)
DEALLOCATE(ZAERTAU)
DEALLOCATE(ZALF)

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RAO_OP',1,ZHOOK_HANDLE)
END SUBROUTINE RAO_OP

