SUBROUTINE AOD_TL&
 & (KFLEVG, KDLEN  , KLEN  , KMXCOUNT, &
 & YDGP_TL,&
 & YDGP5,&
 & PAEROD_PRIME)

!*** *AOD_TL* Tangent Linear Operator routine for aerosol optical depth observations.

!**   INTERFACE.
!     ----------

!     WHERE KFLEVG    = Number of vertical levels
!           KDLEN     = First dimension: Max number of observations, this set (INPUT)
!           KLEN      = Number of aerosol optical depth obs (INPUT)
!           IMXCOUNT  = Number of channels (INPUT)
!           ZQF       =  Perturb in specific humidity at observation points, model levels (INPUT)
!           ZTF       =  Perturb in temperature at observation points, model levels (INPUT)
!           ZPRESH,ZPRESF = Perturb in Half/Full level pressure values at obs points, model levels (INPUT)

!           ZQF5      =  Specific humidity at observation points, model levels (INPUT)
!           ZTF5      =  Temperature at observation points, model levels (INPUT)
!           ZPRESH5,ZPRESF5 = Half/Full level pressure values at obs points, model levels (INPUT)


!           ZXPP      = Aerosol optical depth at different wavelengths (OUTPUT)



!**   AUTHOR.
!     -------
!        ANGELA BENEDETTI        *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2005-06-17
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     A. Geer     29-Dec-2015    OOPS cleaning: remove GEMSPROF and GFL references
!        S. Remy     16-Sep-2016    Add nitrate and ammonium


!     PURPOSE.
!     --------
!     - derive tangent linear aerosol optical thickness from profiles of mass mixing ratio 

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
USE YOMLUN        , ONLY : NULERR
USE YOEAEROP      , ONLY : ALF_SU, ALF_OM, ALF_DD, ALF_SS, ALF_BC, ALF_NI, ALF_AM, ALF_SOA
USE YOEAERATM     , ONLY : YREAERATM
USE YOEAERSNK     , ONLY : YREAERSNK
USE GOM_PLUS      , ONLY : TYPE_GOM_PLUS, IH
USE YOM_GRIB_CODES, ONLY : NGRBAERMR01,NGRBAERMR02, NGRBAERMR03, NGRBAERMR04, NGRBAERMR05, &
 & NGRBAERMR06, NGRBAERMR07, NGRBAERMR08, NGRBAERMR09, NGRBAERMR10, NGRBAERMR11, NGRBAERMR12, &
 & NGRBAERMR13, NGRBAERMR14, NGRBAERMR15, NGRBAERMR16, NGRBAERMR17, NGRBAERMR18, &
 & NGRBAERMR19, NGRBAERMR20, NGRBAERMR21, NGRBAERLG

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

!INTEGER(KIND=JPIM), PARAMETER :: IWAVLMX=7 ! Max number of MODIS channels
!INTEGER(KIND=JPIM), DIMENSION(IWAVLMX) :: IWAVL=(/6, 9, 11, 14, 17, 18, 19/) ! Index for MODIS channels
INTEGER(KIND=JPIM) :: IAER, ITYP, IBIN, ISUM, IWAVLINDEX
INTEGER(KIND=JPIM) :: JAER, JK, JL, JJ, JF, JP, JLEN, JLEV, JTAB
INTEGER(KIND=JPIM) :: IFLAG 
INTEGER(KIND=JPIM), PARAMETER  :: JTABMX=12
INTEGER(KIND=JPIM), ALLOCATABLE :: ITYPE(:) !index to determine optically active aerosol species
INTEGER(KIND=JPIM), DIMENSION(KLEN,KFLEVG) :: IRH
INTEGER(KIND=JPIM):: IIRH, IEFRH

REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG)   :: ZQSAT5, ZT5, ZPF5
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG)   :: ZT, ZPF
REAL(KIND=JPRB), ALLOCATABLE   :: ZAERTAU5(:,:)
REAL(KIND=JPRB), ALLOCATABLE   :: ZAERTAU(:,:)
REAL(KIND=JPRB), ALLOCATABLE   :: ZAERF5(:,:,:) ! aer mixing ratios from traj (interpolated to obs locations)
REAL(KIND=JPRB), ALLOCATABLE   :: ZAERF(:,:,:) ! pert aer mixing ratios from traj (interpolated to obs locations)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG) :: ZAERFT5 ! sum of single contributions of aer mixing ratio
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG) :: ZAERTMR5 ! total aer mixing ratios from traj
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG) :: ZAERTMR ! total aer mixing ratios pert

REAL(KIND=JPRB), ALLOCATABLE :: ZFRACAER5(:,:,:) ! fraction of mixing ratios wrt total mass
REAL(KIND=JPRB) :: ZAEROD5(KDLEN,KMXCOUNT) ! AEROSOL OPT DEPTH 

REAL(KIND=JPRB), ALLOCATABLE :: ZALF5(:)

REAL(KIND=JPRB)    :: ZRHCL5, ZAEREXT5, ZDP5,ZFAC
REAL(KIND=JPRB)    :: ZAEREXT, ZDP

INTEGER(KIND=JPIM) :: INACTAERO

REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL :: LLPHYLIN

#include "abor1.intfb.h"
#include "satur.intfb.h"
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AOD_TL',0,ZHOOK_HANDLE)
ASSOCIATE(&
 & RSS_RH80_MASSFAC=>YREAERATM%RSS_RH80_MASSFAC, &
 & RRHTAB=>YREAERSNK%RRHTAB, &
 & RSSGROWTH_RHTAB=>YREAERSNK%RSSGROWTH_RHTAB, &
 & RSSDENS_RHTAB=>YREAERSNK%RSSDENS_RHTAB, &
 & NMAXTAER=>YREAERATM%NMAXTAER, NTYPAER=>YREAERATM%NTYPAER, &
 & YAERO_DESC=>YREAERATM%YAERO_DESC)

! Check that total number of aerosol bins is not larger than NACTAERO
ISUM=0
DO IAER=1,NMAXTAER
  ISUM=ISUM+NTYPAER(IAER)
ENDDO

! Remember to put in a check for maximum IBIN allowed per aerosol type

ALLOCATE(ZAERF5(KLEN,KFLEVG,YDGP5%NGEMS))
ALLOCATE(ZAERF(KLEN,KFLEVG,YDGP5%NGEMS))
ALLOCATE(ITYPE(YDGP5%NGEMS))
ALLOCATE(ZAERTAU(KLEN,YDGP5%NGEMS))
ALLOCATE(ZAERTAU5(KLEN,YDGP5%NGEMS))
ALLOCATE(ZFRACAER5(KLEN,KFLEVG,YDGP5%NGEMS))
ALLOCATE(ZALF5(YDGP5%NGEMS))

ZAERF5(:,:,:)=0.0_JPRB
ZAERF(:,:,:)=0.0_JPRB

! Initialize ITYPE
ITYPE(:)=0

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
     & YDGP5%GEMS_IGRIB(JF) == NGRBAERMR21) JP=8 ! NO3
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR13) JP=9 ! Volcanic ashes 
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR14 .OR. &
     & YDGP5%GEMS_IGRIB(JF) == NGRBAERMR15)  JP=10 ! Volcanic SO4 and SO2


    IF(JP > 0) THEN
      IAER=IAER+1
      ZAERF5(JLEN,:,IAER)=YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
      ITYPE(IAER)=JP

    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERLG) THEN ! Total aerosol mixing ratio and pert
      ZAERTMR5(JLEN,:)=YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
      ZAERTMR(JLEN,:)=YDGP_TL%GEMS(JLEN,1:KFLEVG,IH,JF)
    ENDIF
  ENDDO ! loop over GEMS species
ENDDO

INACTAERO=IAER
IF (ISUM > INACTAERO) THEN
  CALL ABOR1('AOD_TL:TOTAL NUMBER OF REQUESTED BINS PER TYPE EXCEEDS NAERO')
ENDIF

! Compute total mixing ratio by summing single contributions 
    ZAERFT5(:,:) = 0.0_JPRB
    DO JAER=1,INACTAERO
      IF ((YAERO_DESC(JAER)%NTYP /= 5 .AND. YAERO_DESC(JAER)%NTYP /= 9) .OR. YAERO_DESC(JAER)%NBIN /= 2) THEN
        ZAERFT5(:,:)=ZAERFT5(:,:)+ ZAERF5(:,:,JAER)
      ENDIF
    ENDDO

! BEN the perturbations on the various species are computed from the perturbation
! on the total mass mixing ratio using the mass fraction assumed constant 

ZFRACAER5(:,:,:) = 0.0_JPRB
DO JLEN=1,KLEN
  DO JAER=1,INACTAERO
    IF ((YAERO_DESC(JAER)%NTYP /= 5 .AND. YAERO_DESC(JAER)%NTYP /= 9) .OR. YAERO_DESC(JAER)%NBIN /= 2) THEN
      DO JLEV=1,KFLEVG
        IF(ZAERFT5(JLEN,JLEV) /= 0.0_JPRB) THEN
          ZFRACAER5(JLEN,JLEV,JAER) = ZAERF5(JLEN,JLEV,JAER)/ZAERFT5(JLEN,JLEV)
        ENDIF
        ZAERF(JLEN,JLEV,JAER)=ZAERTMR(JLEN,JLEV) * ZFRACAER5(JLEN,JLEV,JAER)
      ENDDO
    ENDIF
  ENDDO
ENDDO

!! BEN the perturbations on the various species are computed from the perturbation
!! on the total mass mixing ratio using the mass fraction assumed constant (horrible!!!!!)

!DO JLEN=1,KLEN
!  DO JAER=1,INACTAERO
!    DO JLEV=1,KFLEVG
!      IF(ZAERTMR5(JLEN,JLEV) /= 0.0_JPRB) THEN
!        ZFRACAER5(JLEN,JLEV,JAER) = ZAERF5(JLEN,JLEV,JAER)/ZAERTMR5(JLEN,JLEV) 
!      ELSE
!        ZFRACAER5(JLEN,JLEV,JAER) = 0.0_JPRB
!      ENDIF
!        ZAERF(JLEN,JLEV,JAER)=ZAERTMR(JLEN,JLEV) * ZFRACAER5(JLEN,JLEV,JAER)
!    ENDDO
!  ENDDO
!ENDDO

! initialise perturbations to zero 
PAEROD_PRIME(:,:) = 0.0_JPRB 
! initialise trajectory 
ZAEROD5(:,:) = 0.0_JPRB 

! Compute relative humidity at observation points, on model levels

ZT5(1:KLEN,:) = YDGP5%TF(1:KLEN,1:KFLEVG,IH) ! trajectory
ZPF5(1:KLEN,:)= YDGP5%PRESF(1:KLEN,:,IH)     ! trajectory
ZT(1:KLEN,:) = YDGP_TL%TF(1:KLEN,1:KFLEVG,IH)   ! perturb
ZPF(1:KLEN,:)= YDGP_TL%PRESF(1:KLEN,:,IH)       ! perturb
!CALL SATURTL (1 , KLEN , KLEN  , 1 , KFLEVG,&
!  & ZPF5, ZT5 , ZQSAT5,&  
!  & ZPF, ZT , ZQSAT )

IFLAG=2
CALL SATUR (YRTHF, YRCST, 1 , KLEN , KLEN  , 1 , KFLEVG, LLPHYLIN, &
  & ZPF5, ZT5 , ZQSAT5 , IFLAG)

!-- define RH index from "clear-sky" relative humidity
IRH(:,:)=1
DO JK=1,KFLEVG
  DO JL=1,KLEN
    ZRHCL5=(YDGP5%QF(JL,JK,IH)/ZQSAT5(JL,JK))*100._JPRB
    DO JTAB=1,JTABMX
      IF (ZRHCL5 > RRHTAB(JTAB)) THEN
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
!         depth retrievals (Remer et al. 2005 and Benedetti et al., 2009, ECMWF Tech. Memo 571).
!
!      IWAVLINDEX=IWAVL(JJ)

  DO JL=1,KLEN
    DO JAER=1,INACTAERO

!-- ITYP is the aerosol type 1:SeaSalt,   2:DesertDust,   3:OrganicMatter, 4: BlackCarbon,
!                            5:SUlphates, 6:FlyingAshes
!   IBIN is the bin index: 1-3:SS, 1-3:DD, 1:OM, 1:BC, 1:SU, 1:FA

     ITYP=YAERO_DESC(JAER)%NTYP
     IBIN=YAERO_DESC(JAER)%NBIN
     IF(ITYP /= 0 .AND. IBIN == 0) THEN
      WRITE(NULERR,*) 'AT LEAST ONE BIN NEEDS TO BE ALLOCATED FOR TYPE ', ITYP
      WRITE(NULERR,*) 'PROGRAM WILL STOP'
      CALL ABOR1('') 
     ENDIF
     IF(ITYP /= ITYPE(JAER)) THEN
      WRITE(NULERR,*)'AOD_TL: AEROSOL TYPE IS NOT MATCHING', ITYP
      WRITE(NULERR,*) 'PROGRAM WILL STOP'
      CALL ABOR1('')
     ENDIF

     ZAERTAU(JL,JAER)=0.0_JPRB
     ZAERTAU5(JL,JAER)=0.0_JPRB

        
      DO JK=1,KFLEVG

         IIRH=IRH(JL,JK)
         ZFAC=1.0_JPRB

         IF (ITYP == 1) THEN
          ZALF5(JAER)=ALF_SS(IIRH,IWAVLINDEX,IBIN)
          ZFAC=RSS_RH80_MASSFAC
         ELSEIF (ITYP == 2) THEN
          ZALF5(JAER)=ALF_DD(IBIN,IWAVLINDEX)
         ELSEIF (ITYP == 3) THEN
      !-- for bin 2 (hydrophobic), use the 20% value of the OM optical properties
          IF (IBIN == 2) IIRH=IEFRH
          ZALF5(JAER)=ALF_OM(IIRH,IWAVLINDEX)
         ELSEIF (ITYP == 4) THEN
          ZALF5(JAER)=ALF_BC(IWAVLINDEX)
          ELSEIF (ITYP == 5 .OR. ITYP == 10) THEN
          ZALF5(JAER)=ALF_SU(IIRH,IWAVLINDEX)
          !-- SO2 does not contribute to optical depth, only SO4 does.
          IF (IBIN == 2) THEN
            ZALF5(JAER)=0._JPRB
          ENDIF
         ELSEIF (ITYP == 6) THEN ! Nitrate
            ZALF5(JAER)=ALF_NI(IIRH,IWAVLINDEX,IBIN)
         ELSEIF (ITYP == 7) THEN ! Ammonium
            ZALF5(JAER)=ALF_AM(IIRH,IWAVLINDEX)
         ELSEIF (ITYP == 8) THEN ! SOA
            ZALF5(JAER)=ALF_SOA(IIRH,IWAVLINDEX,IBIN)
         ELSEIF (ITYP == 9) THEN ! volcanic ashes
!-- use dust for 0.9-20 um bin
            ZALF5(JAER)=ALF_DD(3,IWAVLINDEX)
         ENDIF


            ZALF5(JAER)=ZALF5(JAER)*1000._JPRB ! conversion from m^2/g to m^2/kg

            ZDP5=YDGP5%PRESH(JL,JK,IH)-YDGP5%PRESH(JL,JK-1,IH)
            ZDP=YDGP_TL%PRESH(JL,JK,IH)-YDGP_TL%PRESH(JL,JK-1,IH)

        ZAEREXT5=ZAERF5(JL,JK,JAER)*ZFAC*ZALF5(JAER)
        ZAEREXT=ZAERF(JL,JK,JAER)*ZFAC*ZALF5(JAER)

        ZAERTAU5(JL,JAER)  =ZAERTAU5(JL,JAER)+ZAEREXT5 * (ZDP5/RG)
        ZAERTAU (JL,JAER)  =ZAERTAU (JL,JAER)+ZAEREXT * (ZDP5/RG) + ZAEREXT5 *(ZDP/RG)

     ENDDO ! loop over model vertical levels

      ZAEROD5(JL,JJ)=ZAEROD5(JL,JJ)+ZAERTAU5(JL,JAER)
      PAEROD_PRIME(JL,JJ)=PAEROD_PRIME(JL,JJ)+ZAERTAU(JL,JAER)
      
    ENDDO ! loop over aerosol type
  ENDDO ! loop over horizontal observational points
ENDDO ! loop over MODIS channels

! Security check for negative optical depths 
   DO  JJ = 1,KMXCOUNT
    DO JL=1, KLEN
     IF(ZAEROD5(JL,JJ) < 0.0_JPRB) THEN
      PAEROD_PRIME(JL,JJ)= 0.0_JPRB
     ENDIF
    ENDDO
   ENDDO

DEALLOCATE(ZAERF5)
DEALLOCATE(ZAERF)
DEALLOCATE(ITYPE)
DEALLOCATE(ZAERTAU)
DEALLOCATE(ZAERTAU5)
DEALLOCATE(ZFRACAER5)
DEALLOCATE(ZALF5)
!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AOD_TL',1,ZHOOK_HANDLE)
END SUBROUTINE AOD_TL
