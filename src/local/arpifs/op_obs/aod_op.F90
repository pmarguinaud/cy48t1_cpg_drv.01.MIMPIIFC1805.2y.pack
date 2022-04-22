SUBROUTINE AOD_OP&
 & (KFLEVG, KDLEN  , KLEN  , KMXCOUNT, &
 & YDGP5,&
 & PAEROD )

!*** *AOD_OP* Operator routine for aerosol optical depth observations.

!**   INTERFACE.
!     ----------
!          *AOD_OP* IS CALLED FROM *HOP*|*OBSOP_COMPOSITION* (oops design).

!     WHERE KFLEVG    = Number of vertical levels
!           KDLEN     = First dimension: Max number of observations, this set (INPUT)
!           KLEN      = Number of aerosol optical depth obs (INPUT)
!           KMXCOUNT  = Number of channels (INPUT)
!           PQF       = Specific humidity at observation points, model levels (INPUT)
!           PTF       = Temperature at observation points, model levels (INPUT)
!           PPRESH,PPRESF = Half/Full level pressure values at obs points, model levels (INPUT)
!           PGEMSPROF = GEMS variables (INPUT)
!           PAEROD    = Aerosol optical depth at different wavelengths (OUTPUT)



!**   AUTHOR.
!     -------
!        ANGELA BENEDETTI        *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2005-06-17 (adapted from aer_bdgmass.F90 by J.-J. Morcrette)

!        A. Benedetti: 2009-02-25. Readaptation of code to allow for two sulphate variables
!                                   with SO2 not contributing to the total aerosol optical depth
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!        C. Payan: Nov 2015 CY42 oops!
!        A. Geer     29-Dec-2015    OOPS cleaning: remove GEMSPROF and GFL references
!        S. Remy     16-Sep-2016    Add nitrate and ammonium

!     PURPOSE.
!     --------
!     - derive aerosol optical thickness from profiles of mass mixing ratio 

! Optical properties in su_aerop.F90 avaliable at following wavelenghts
! (340)355(380)400 440 469 500 532 550 645 670 800 858 865(1020)1064 1240 1640 2130 10000
! corresponding IWAVL 
!   1   2   3   4   5   6   7   8   9  10   11  12  13  14  15  16   17    18   19   20
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

USE PARKIND1       , ONLY : JPIM, JPRB
USE YOMHOOK        , ONLY : LHOOK, DR_HOOK
USE YOMDB
USE YOMCST         , ONLY : RG, YRCST
USE YOETHF         , ONLY : YRTHF
USE YOMLUN         , ONLY : NULERR
USE YOEAEROP       , ONLY : ALF_SU, ALF_OM, ALF_DD, ALF_SS, ALF_BC, ALF_NI, ALF_AM, ALF_SOA
USE YOEAERATM      , ONLY : YREAERATM
USE YOEAERSNK      , ONLY : YREAERSNK
USE GOM_PLUS       , ONLY : TYPE_GOM_PLUS, IH
USE YOM_GRIB_CODES , ONLY : NGRBAERMR01,NGRBAERMR02, NGRBAERMR03, NGRBAERMR04, NGRBAERMR05, &
 & NGRBAERMR06, NGRBAERMR07, NGRBAERMR08, NGRBAERMR09, NGRBAERMR10, NGRBAERMR11, NGRBAERMR12, &
 & NGRBAERMR13, NGRBAERMR14, NGRBAERMR15, NGRBAERMR16, NGRBAERMR17, NGRBAERMR18, &
 & NGRBAERMR19, NGRBAERMR20, NGRBAERMR21, NGRBAERLG
USE YOEPHLI        , ONLY : YREPHLI

IMPLICIT NONE


!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEVG, KDLEN, KLEN, KMXCOUNT

TYPE (TYPE_GOM_PLUS),INTENT(IN)  :: YDGP5

REAL(KIND=JPRB)   ,INTENT(OUT) :: PAEROD(KDLEN,KMXCOUNT) ! AEROSOL OPT DEPTH (FUNCTION OF HORIZONTAL POSITION AND CHANNEL #)



!*       0.5   LOCAL VARIABLES
!              ---------------
!INTEGER(KIND=JPIM), PARAMETER :: IWAVLMX=7 ! Max number of MODIS channels
! IWAVL: Index for MODIS channels - based on latest su_aerop.F90
!INTEGER(KIND=JPIM), DIMENSION(IWAVLMX) :: IWAVL=(/6, 9, 11, 14, 17, 18, 19/) 
INTEGER(KIND=JPIM) :: IAER, ITYP, IBIN, IFLAG, IWAVLINDEX
INTEGER(KIND=JPIM) :: JAER, JK, JL, JJ, JF, JP, JLEN, JTAB
INTEGER(KIND=JPIM), PARAMETER  :: JTABMX=12
INTEGER(KIND=JPIM), ALLOCATABLE :: ITYPE(:) !index to determine optically active aerosol species
INTEGER(KIND=JPIM), DIMENSION(KLEN,KFLEVG) :: IRH

REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG)   :: ZQSAT, ZT, ZPF
REAL(KIND=JPRB), ALLOCATABLE   :: ZAERTAU(:,:)
REAL(KIND=JPRB), ALLOCATABLE   :: ZAERF(:,:,:) ! aer mixing ratios from traj (interpolated to obs locations)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG) :: ZAERTMR ! total aer mixing ratios from traj
REAL(KIND=JPRB), ALLOCATABLE :: ZALF(:)
REAL(KIND=JPRB)    :: ZRHCL(KLEN,KFLEVG), ZAEREXT, ZDP, ZFAC
INTEGER(KIND=JPIM) :: IIRH, IEFRH

INTEGER(KIND=JPIM) :: INACTAERO

REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL :: LLPHYLIN 

#include "abor1.intfb.h"
#include "satur.intfb.h"
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AOD_OP',0,ZHOOK_HANDLE)
ASSOCIATE(&
 & RSS_RH80_MASSFAC=>YREAERATM%RSS_RH80_MASSFAC, &
 & RRHTAB=>YREAERSNK%RRHTAB, &
 & RSSGROWTH_RHTAB=>YREAERSNK%RSSGROWTH_RHTAB, &
 & RSSDENS_RHTAB=>YREAERSNK%RSSDENS_RHTAB,&
 & NMAXTAER=>YREAERATM%NMAXTAER, NTYPAER=>YREAERATM%NTYPAER, &
 & YAERO_DESC=>YREAERATM%YAERO_DESC)

LLPHYLIN = YREPHLI%LPHYLIN

!-----------------------------------------------------------------------
! Remember to put in a check for maximum IBIN allowed per aerosol type

ALLOCATE(ZAERF(KLEN,KFLEVG,YDGP5%NGEMS))
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
     & YDGP5%GEMS_IGRIB(JF) == NGRBAERMR12)  JP=5 ! SO4 and SO2
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
      ZAERF(JLEN,:,IAER)=YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
      ITYPE(IAER)=JP
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERLG) THEN ! Total aerosol mixing ratio
      ZAERTMR(JLEN,:)=YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
    ENDIF
  ENDDO
ENDDO
INACTAERO=IAER

! initialise optical depth 
PAEROD(:,:) = 0.0_JPRB 

! Compute relative humidity at observation points, on model levels

ZT(1:KLEN,:) = YDGP5%TF(1:KLEN,1:KFLEVG,IH)
ZPF(1:KLEN,:)= YDGP5%PRESF(1:KLEN,:,IH)
IFLAG=2
CALL SATUR (YRTHF, YRCST, 1 , KLEN , KLEN  , 1 , KFLEVG, LLPHYLIN, &
  & ZPF, ZT , ZQSAT , IFLAG)

!-- define RH index from "clear-sky" relative humidity
IRH(:,:)=1
DO JK=1,KFLEVG
  DO JL=1,KLEN
    ZRHCL(JL,JK)=(YDGP5%QF(JL,JK,IH)/ZQSAT(JL,JK))*100._JPRB
    DO JTAB=1,JTABMX
      IF (ZRHCL(JL,JK) > RRHTAB(JTAB)) THEN
        IRH(JL,JK)=JTAB
      ENDIF
    ENDDO
  ENDDO
ENDDO

!-- the effective relative humidity is the low value (20%) assumed for hydrophobic component of OM
IEFRH=3
  
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
! This indeces are defined in su_aerw.F90 

         ITYP=YAERO_DESC(JAER)%NTYP
         IBIN=YAERO_DESC(JAER)%NBIN
         IF(ITYP /= 0 .AND. IBIN == 0) THEN
          WRITE(NULERR,*)'AOD_OP:AT LEAST ONE BIN NEEDS TO BE ALLOCATED FOR TYPE ', ITYP
          WRITE(NULERR,*) 'PROGRAM WILL STOP'
          CALL ABOR1('') 
         ENDIF
         IF(ITYP /= ITYPE(JAER)) THEN
          WRITE(NULERR,*)'AOD_OP: AEROSOL TYPE IS NOT MATCHING', ITYP, ITYPE(JAER)
          WRITE(NULERR,*) 'PROGRAM WILL STOP'
          CALL ABOR1('')
         ENDIF

      ZAERTAU(JL,JAER)=0.0_JPRB
        
      DO JK=1,KFLEVG

! Always double check dimensioning of arrays ALF_** in yoeaerop.F90 !!!!!!!!!
! It has changed a number of times!!!!!!!!
      
         IIRH=IRH(JL,JK)
         ZFAC=1.0_JPRB

         IF (ITYP == 1) THEN
          ZALF(JAER)=ALF_SS(IIRH,IWAVLINDEX,IBIN)
          ZFAC=RSS_RH80_MASSFAC
         ELSEIF (ITYP == 2) THEN
          ZALF(JAER)=ALF_DD(IBIN,IWAVLINDEX)
         ELSEIF (ITYP == 3) THEN
      !-- for bin 2 (hydrophobic), use the 20% value of the OM optical properties
          IF (IBIN == 2) IIRH=IEFRH
          ZALF(JAER)=ALF_OM(IIRH,IWAVLINDEX)
         ELSEIF (ITYP == 4) THEN
          ZALF(JAER)=ALF_BC(IWAVLINDEX)
         ELSEIF (ITYP == 4) THEN
          ZALF(JAER)=ALF_BC(IWAVLINDEX)
         ELSEIF (ITYP == 5 .OR. ITYP == 10) THEN
          ZALF(JAER)=ALF_SU(IIRH,IWAVLINDEX)
          !-- SO2 does not contribute to optical depth, only SO4 does.
          IF (IBIN == 2) THEN
            ZALF(JAER)=0._JPRB
          ENDIF
         ELSEIF (ITYP == 6) THEN ! Nitrate
            ZALF(JAER)=ALF_NI(IIRH,IWAVLINDEX,IBIN)
         ELSEIF (ITYP == 7) THEN ! Ammonium
            ZALF(JAER)=ALF_AM(IIRH,IWAVLINDEX)
         ELSEIF (ITYP == 8) THEN ! SOA
            ZALF(JAER)=ALF_SOA(IIRH,IWAVLINDEX,IBIN)
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
IF (LHOOK) CALL DR_HOOK('AOD_OP',1,ZHOOK_HANDLE)
END SUBROUTINE AOD_OP

