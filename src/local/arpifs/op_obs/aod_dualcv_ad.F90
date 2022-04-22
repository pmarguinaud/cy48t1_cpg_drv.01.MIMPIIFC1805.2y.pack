SUBROUTINE AOD_DUALCV_AD&
 & (KFLEVG, KDLEN  , KLEN  , KMXCOUNT, &
 & YDGP_AD,&
 & YDGP5,&
 & PAEROD_HAT)

!*** *AOD_DUALCV_AD* Adjoint Operator routine for aerosol optical depth observations.

!**   INTERFACE.
!     ----------

!     WHERE KFLEVG    = Number of vertical levels
!           KDLEN      = First dimension: Max number of observations, this set (INPUT)
!           KLEN       = Number of aerosol optical depth obs (INPUT)
!           IMXCOUNT   = Number of channels (INPUT)
!           ZQF       =  Perturb in specific humidity at observation points, model levels (INPUT)
!           ZTF       =  Perturb in temperature at observation points, model levels (INPUT)
!           ZPRESH,ZPRESF = Perturb in Half/Full level pressure values at obs points, model levels (INPUT)
!           ZQF5       =  Specific humidity at observation points, model levels (INPUT)
!           ZTF5       =  Temperature at observation points, model levels (INPUT)
!           ZPRESH5,ZPRESF5 =   Half/Full level pressure values at obs points, model levels (INPUT)

!           ZXPP       = Aerosol optical depth at different wavelengths (OUTPUT)



!**   AUTHOR.
!     -------
!        ANGELA BENEDETTI        *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2005-06-17
! O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     A. Geer     29-Dec-2015    OOPS cleaning: remove GEMSPROF and GFL references
!        S. Remy     16-Sep-2016    Add nitrate and ammonium


!     PURPOSE.
!     --------
!     - derive adjoint aerosol optical thickness from profiles of mass mixing ratio 

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
USE YOEAERSNK     , ONLY : YREAERSNK
USE YOEAERATM     , ONLY : YREAERATM
USE GOM_PLUS      , ONLY : TYPE_GOM_PLUS, IH
USE YOM_GRIB_CODES, ONLY : NGRBAERMR01,NGRBAERMR02, NGRBAERMR03, NGRBAERMR04, NGRBAERMR05, &
 & NGRBAERMR06, NGRBAERMR07, NGRBAERMR08, NGRBAERMR09, NGRBAERMR10, NGRBAERMR11, &
 & NGRBAERMR16, NGRBAERMR17, NGRBAERMR18, NGRBAERMR19, NGRBAERMR20, &
 & NGRBAERLG, NGRBAERSM
USE YOEPHLI       , ONLY : YREPHLI

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KFLEVG, KDLEN, KLEN, KMXCOUNT

TYPE (TYPE_GOM_PLUS),INTENT(INOUT) :: YDGP_AD
TYPE (TYPE_GOM_PLUS),INTENT(IN)    :: YDGP5

REAL(KIND=JPRB)     ,INTENT(INOUT) :: PAEROD_HAT(KDLEN,KMXCOUNT) ! AEROSOL OPT DEPTH GRADIENT 

!*       0.5   LOCAL VARIABLES
!              ---------------
!INTEGER(KIND=JPIM), PARAMETER :: IWAVLMX=7 ! Max number of MODIS channels
!INTEGER(KIND=JPIM), DIMENSION(IWAVLMX) :: IWAVL=(/6, 9, 11, 14, 17, 18, 19/) ! Index for MODIS channels
INTEGER(KIND=JPIM) :: IACTAEROP ! Number of optically active aerosol species (YGFL%NACTAERO-1, for now)
INTEGER(KIND=JPIM) :: IAER, ITYP, IBIN, ISUM, IWAVLINDEX
INTEGER(KIND=JPIM) :: JAER, JK, JL, JJ, JF, JLEN, JLEV, JTAB
INTEGER(KIND=JPIM) :: IFLAG 
INTEGER(KIND=JPIM), PARAMETER  :: JTABMX=12
INTEGER(KIND=JPIM), ALLOCATABLE :: ITYPE(:) !index to determine optically active aerosol species
INTEGER(KIND=JPIM), DIMENSION(KLEN,KFLEVG) :: IRH

REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG)   :: ZQSAT5, ZT5, ZPF5
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG)   :: ZT, ZPF

REAL(KIND=JPRB), ALLOCATABLE   :: ZAERTAU5(:,:)
REAL(KIND=JPRB), ALLOCATABLE   :: ZAERTAU(:,:)

REAL(KIND=JPRB), ALLOCATABLE   :: ZAER5(:,:,:) ! single-species aer mixing ratios from traj (interpolated to obs locations)
REAL(KIND=JPRB), ALLOCATABLE   :: ZAER(:,:,:) ! single-species pert aer mixing ratios from traj (interpolated to obs locations)
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG) :: ZAERFT5 ! sum of single contributions of fine aer mixing ratio
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG) :: ZAERCT5 ! sum of single contributions of coarse aer mixing ratio
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG) :: ZAERSMMR5 ! fine aer mixing ratios from traj
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG) :: ZAERSMMR ! fine aer mixing ratios pert
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG) :: ZAERLGMR5 ! coarse aer mixing ratios from traj
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG) :: ZAERLGMR !  coarse aer mixing ratios pert


REAL(KIND=JPRB), ALLOCATABLE :: ZFRAERSM5(:,:,:) ! fraction of mixing ratios wrt fine mass
REAL(KIND=JPRB), ALLOCATABLE :: ZFRAERLG5(:,:,:) ! fraction of mixing ratios wrt fine mass

INTEGER(KIND=JPIM)              :: JPACTAERSM=7, JPACTAERLG=4
INTEGER(KIND=JPIM), ALLOCATABLE :: IACTAERSM(:), IACTAERLG(:)

REAL(KIND=JPRB), DIMENSION(KDLEN,KMXCOUNT) :: ZAEROD5 ! aer optical depth
REAL(KIND=JPRB), ALLOCATABLE :: ZALF5(:)
REAL(KIND=JPRB)   :: ZRHCL5, ZAEREXT5
REAL(KIND=JPRB), DIMENSION(KLEN,KFLEVG):: ZDP5 ! pressure difference from trajectory
REAL(KIND=JPRB)   :: ZAEREXT, ZDP, ZFAC
INTEGER(KIND=JPIM) :: INACTAERO, INAER
REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL :: LLPHYLIN 

#include "satur.intfb.h"
#include "abor1.intfb.h"
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AOD_DUALCV_AD',0,ZHOOK_HANDLE)
ASSOCIATE(&
 & RSS_RH80_MASSFAC=>YREAERATM%RSS_RH80_MASSFAC, &
 & RRHTAB=>YREAERSNK%RRHTAB, &
 & NMAXTAER=>YREAERATM%NMAXTAER, NTYPAER=>YREAERATM%NTYPAER, &
 & YAERO_DESC=>YREAERATM%YAERO_DESC)

LLPHYLIN = YREPHLI%LPHYLIN

! Check that total number of aerosol bins is not larger than NACTAERO
ISUM=0
DO IAER=1,NMAXTAER
  ISUM=ISUM+NTYPAER(IAER)
ENDDO


ALLOCATE(ZAER5(KLEN,KFLEVG,YDGP5%NGEMS))
ALLOCATE(ZAER(KLEN,KFLEVG,YDGP5%NGEMS))
ALLOCATE(ITYPE(YDGP5%NGEMS))
ALLOCATE(ZAERTAU(KLEN,YDGP5%NGEMS))
ALLOCATE(ZAERTAU5(KLEN,YDGP5%NGEMS))
ALLOCATE(ZFRAERSM5(KLEN,KFLEVG,YDGP5%NGEMS))
ALLOCATE(ZFRAERLG5(KLEN,KFLEVG,YDGP5%NGEMS))
ALLOCATE(ZALF5(YDGP5%NGEMS))

ZAER5(:,:,:)=0.0_JPRB
ZAER(:,:,:)=0.0_JPRB
ITYPE(:)=0

! Extract aerosol mixing ratios from GEMS profiles structure
DO JLEN=1,KLEN
  INAER=0
  DO JF=1,YDGP5%NGEMS
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR01) THEN ! Sea-salt first bin
     IAER=1
     ITYPE(IAER)=1   
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR02) THEN ! Sea-salt second bin
     IAER=2
     ITYPE(IAER)=1
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR03) THEN ! Sea-salt third bin
     IAER=3
     ITYPE(IAER)=1
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR04) THEN ! Desert dust first bin
     IAER=4
     ITYPE(IAER)=2
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR05) THEN ! Desert dust second bin
     IAER=5
     ITYPE(IAER)=2
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR06) THEN ! Desert dust third bin
     IAER=6
     ITYPE(IAER)=2
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR07) THEN ! Organic matter hydrophobic 
     IAER=7
     ITYPE(IAER)=3
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR08) THEN ! Organic matter hydrophilic
     IAER=8
     ITYPE(IAER)=3
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR09) THEN ! Black carbon hydrophobic
     IAER=9
     ITYPE(IAER)=4
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR10) THEN ! Black carbon hydrophilic
     IAER=10
     ITYPE(IAER)=4
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR11) THEN ! Sulphate  (SO4)
     IAER=11
     ITYPE(IAER)=5
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR16) THEN ! Nitrate bin 1
     IAER=12
     ITYPE(IAER)=6
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR17) THEN ! Nitrate bin 2
     IAER=13
     ITYPE(IAER)=6
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR18) THEN ! Ammonium
     IAER=14
     ITYPE(IAER)=7
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR19) THEN ! SOA Biogenic
     IAER=15
     ITYPE(IAER)=8
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR20) THEN ! SOA Anthro 1
     IAER=16
     ITYPE(IAER)=8
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR20) THEN ! SOA Anthro 2
     IAER=17
     ITYPE(IAER)=8
     ZAER5(JLEN,:,IAER) = YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
     INAER=INAER+1
    ENDIF

    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERLG) THEN ! Coarse aerosol mixing ratio and pert
      ZAERLGMR5(JLEN,:)=YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERSM) THEN ! Fine aerosol mixing ratio and pert
      ZAERSMMR5(JLEN,:)=YDGP5%GEMS(JLEN,1:KFLEVG,IH,JF)
    ENDIF
  ENDDO ! loop over GEMS species
ENDDO

INACTAERO=INAER
IF (ISUM > INACTAERO) THEN
  CALL ABOR1('AOD_DUALCV_TL:TOTAL NUMBER OF REQUESTED BINS PER TYPE EXCEEDS NAERO')
ENDIF

IACTAEROP=INACTAERO-1 ! Number of optically active aerosol species

! Check that sum of active fine mode and coarse mode does not exceed total number of 
! optically active species

 IF(( JPACTAERSM + JPACTAERLG) > IACTAEROP)  THEN
  CALL ABOR1('AOD_DUALCV_TL:SUM OF FINE AND COARSE EXCEEDS IACTAEROP')
 ENDIF

! Compute fine mode mixing ratio by summing single contributions
    ZAERFT5(:,:) = 0.0_JPRB
    ZAERFT5(:,:) = ZAER5(:,:,1)+ ZAER5(:,:,4)  + ZAER5(:,:,7) +ZAER5(:,:,8)+ &
    &              ZAER5(:,:,9)+ ZAER5(:,:,10) + ZAER5(:,:,11)
! Nitrates
IF (NTYPAER(6)+NTYPAER(7) > 0) THEN
ZAERFT5(:,:) = ZAERFT5(:,:) + ZAER5(:,:,12)+ZAER5(:,:,14)
JPACTAERSM=JPACTAERSM+2
ENDIF
! SOA
IF (NTYPAER(8) > 0) THEN
ZAERFT5(:,:) = ZAERFT5(:,:) + ZAER5(:,:,15)+ZAER5(:,:,16)+ZAER5(:,:,17)
JPACTAERSM=JPACTAERSM+3
ENDIF

! Compute coarse mode mixing ratio by summing single contributions
    ZAERCT5(:,:) = 0.0_JPRB
    ZAERCT5(:,:) = ZAER5(:,:,2)+ ZAER5(:,:,3)  + ZAER5(:,:,5) +ZAER5(:,:,6)

! Nitrates
IF (NTYPAER(6)+NTYPAER(7) > 0) THEN
ZAERCT5(:,:) = ZAERCT5(:,:) + ZAER5(:,:,13)
JPACTAERLG=JPACTAERLG+1
ENDIF

! BEN the perturbations on the various species are computed from the perturbation
! on the total mass mixing ratio using the mass fraction assumed constant

ZFRAERSM5(:,:,:) = 0.0_JPRB
ZFRAERLG5(:,:,:) = 0.0_JPRB

ALLOCATE(IACTAERSM(JPACTAERSM))
ALLOCATE(IACTAERLG(JPACTAERLG))

IACTAERSM(1)=1
IACTAERSM(2)=4
IACTAERSM(3)=7
IACTAERSM(4)=8
IACTAERSM(5)=9
IACTAERSM(6)=10
IACTAERSM(7)=11
! Nitrates
IF (NTYPAER(6)+NTYPAER(7) > 0) THEN
  IACTAERSM(8)=12
  IACTAERSM(9)=14
  IF (NTYPAER(8) > 0) THEN
    ! SOA
    IACTAERSM(10)=12
    IACTAERSM(11)=14
   ENDIF
ELSEIF (NTYPAER(8) > 0) THEN
  ! SOA
  IACTAERSM(8)=12
  IACTAERSM(9)=14
ENDIF

IACTAERLG(1)=2
IACTAERLG(2)=3
IACTAERLG(3)=5
IACTAERLG(4)=6
! Nitrates
IF (NTYPAER(6)+NTYPAER(7) > 0) THEN
  IACTAERLG(5)=13
ENDIF

DO JLEN=1,KLEN
  DO JLEV=1,KFLEVG

! Fine mode
    DO JAER=1,JPACTAERSM
      IF(ZAERFT5(JLEN,JLEV) /= 0.0_JPRB) THEN
        ZFRAERSM5(JLEN,JLEV,IACTAERSM(JAER)) = ZAER5(JLEN,JLEV,IACTAERSM(JAER))/ZAERFT5(JLEN,JLEV)
      ENDIF
    ENDDO
! Coarse mode
    DO JAER=1,JPACTAERLG
      IF(ZAERCT5(JLEN,JLEV) /= 0.0_JPRB) THEN
        ZFRAERLG5(JLEN,JLEV,IACTAERLG(JAER)) = ZAER5(JLEN,JLEV,IACTAERLG(JAER))/ZAERCT5(JLEN,JLEV)
      ENDIF
    ENDDO

  ENDDO
ENDDO


!! Compute fraction of given species (and bin) wrt total mass mixing ratio

!DO JLEN=1,KLEN
!  DO JAER=1,IACTAEROP
!    DO JLEV=1,KFLEVG
!      IF(ZAERTMR5(JLEN,JLEV) /= 0.0_JPRB) THEN
!        ZFRACAER5(JLEN,JLEV,JAER) = ZAER5(JLEN,JLEV,JAER)/ZAERTMR5(JLEN,JLEV)
!      ELSE
!        ZFRACAER5(JLEN,JLEV,JAER) = 0.0_JPRB
!      ENDIF
!    ENDDO
!  ENDDO
!ENDDO

! Compute relative humidity at observation points, on model levels
ZT5(1:KLEN,:) = YDGP5%TF(1:KLEN,1:KFLEVG,IH) ! trajectory
ZPF5(1:KLEN,:)= YDGP5%PRESF(1:KLEN,:,IH)     ! trajectory
ZT(1:KLEN,:) = 0.0_JPRB       ! perturb
ZPF(1:KLEN,:)= 0.0_JPRB       ! perturb
!CALL SATURAD (1 , KLEN , KLEN  , 1 , KFLEVG,&
!  & ZPF5, ZT5 , ZQSAT5, &
!  & ZPF,  ZT ,  ZQSAT)

IFLAG=2
CALL SATUR (YRTHF, YRCST, 1 , KLEN , KLEN  , 1 , KFLEVG, LLPHYLIN, &
  & ZPF5, ZT5 , ZQSAT5 , IFLAG)


!PPRESF(1:KLEN,:)=PPRESF(1:KLEN,:)+ZPF(1:KLEN,:)     ! perturb in pressure - do not include
!PTF(1:KLEN,:)= PTF(1:KLEN,:) + ZT(1:KLEN,:)         ! perturb in temperature - do not include 
!ZT(1:KLEN,:) = 0.0_JPRB       
!ZPF(1:KLEN,:)= 0.0_JPRB   

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
! Save the pressure difference from trajectory 
    ZDP5(JL,JK)=YDGP5%PRESH(JL,JK,IH)-YDGP5%PRESH(JL,JK-1,IH)
  ENDDO
ENDDO

    IF(KMXCOUNT ==1) IWAVLINDEX=9   ! 550nm

 ZAEROD5(:,:)=0.0_JPRB
! Trajectory computations
 DO JJ=1,KMXCOUNT

!   IWAVLINDEX is the wavelength index: for now the choice of channels is tailored to MODIS optical
!         depth retrievals (Remer et al. 2005 and Benedetti et al, 2009, ECMWF Tech. Memo 571).
!
!      IWAVLINDEX=IWAVL(JJ)

  DO JL=1,KLEN
   DO JAER=1,IACTAEROP

      ZAERTAU5(JL,JAER)= 0.0_JPRB

      ITYP=YAERO_DESC(JAER)%NTYP
      IBIN=YAERO_DESC(JAER)%NBIN
      IF(ITYP /= 0 .AND. IBIN == 0) THEN
       WRITE(NULERR,*) 'AT LEAST ONE BIN NEEDS TO BE ALLOCATED FOR TYPE ', ITYP,'!'
       WRITE(NULERR,*) 'PROGRAM WILL STOP'
       CALL ABOR1('') 
      ENDIF
     IF(ITYP /= ITYPE(JAER)) THEN
      WRITE(NULERR,*)'AOD_DUALCV_AD: AEROSOL TYPE IS NOT MATCHING', ITYP
      WRITE(NULERR,*) 'PROGRAM WILL STOP'
      CALL ABOR1('')
     ENDIF


      DO JK=1,KFLEVG
      
         ZFAC=1.0_JPRB
         IF (ITYP == 1) THEN
          ZALF5(JAER)=ALF_SS(IRH(JL,JK),IWAVLINDEX,IBIN)
          ZFAC=RSS_RH80_MASSFAC
         ELSEIF (ITYP == 2) THEN
          ZALF5(JAER)=ALF_DD(IBIN,IWAVLINDEX)
         ELSEIF (ITYP == 3) THEN
          ZALF5(JAER)=ALF_OM(IRH(JL,JK),IWAVLINDEX)
         ELSEIF (ITYP == 4) THEN
          ZALF5(JAER)=ALF_BC(IWAVLINDEX)
         ELSEIF (ITYP == 5 .OR. ITYP == 10) THEN
          ZALF5(JAER)=ALF_SU(IRH(JL,JK),IWAVLINDEX)
          !-- SO2 does not contribute to optical depth, only SO4 does.
          IF (IBIN == 2) THEN
            ZALF5(JAER)=0._JPRB
          ENDIF
         ELSEIF (ITYP == 6) THEN ! Nitrate
            ZALF5(JAER)=ALF_NI(IRH(JL,JK),IWAVLINDEX,IBIN)
         ELSEIF (ITYP == 7) THEN ! Ammonium
            ZALF5(JAER)=ALF_AM(IRH(JL,JK),IWAVLINDEX)
         ELSEIF (ITYP == 8) THEN ! SOA
            ZALF5(JAER)=ALF_SOA(IRH(JL,JK),IWAVLINDEX,IBIN)
         ELSEIF (ITYP == 9) THEN ! volcanic ashes
!-- use dust for 0.9-20 um bin
            ZALF5(JAER)=ALF_DD(3,IWAVLINDEX)
         ENDIF

! convert units from m^2/g to m^2/kg
         ZALF5(JAER)=ZALF5(JAER)*1000._JPRB

         ZAEREXT5=ZAER5(JL,JK,JAER)*ZFAC*ZALF5(JAER)
         ZAERTAU5(JL,JAER)  =ZAERTAU5(JL,JAER)+ZAEREXT5 * (ZDP5(JL,JK)/RG)
 
     ENDDO  ! loop over vertical levels
      ZAEROD5(JL,JJ)=ZAEROD5(JL,JJ)+ZAERTAU5(JL,JAER)
    ENDDO ! loop over aerosol type
   ENDDO ! loop over horizontal observational points
 ENDDO ! loop over MODIS channels


! Security check for negative optical depths (adjoint)
   DO  JJ = 1,KMXCOUNT
    DO JL=1, KLEN
     IF(ZAEROD5(JL,JJ) < 0.0_JPRB) THEN
      PAEROD_HAT(JL,JJ)= 0.0_JPRB
     ENDIF
    ENDDO
   ENDDO

!---------------------------------------------------------
! Adjoint calculations
!---------------------------------------------------------
DO JJ=1,KMXCOUNT

!   IWAVLINDEX is the wavelength index: for now the choice of channels is tailored to MODIS optical
!         depth retrievals (Remer et al. 2005 and Benedetti et al, 2009, ECMWF Tech. Memo 571)
!
!      IWAVLINDEX=IWAVL(JJ)

  DO JL=1,KLEN
    DO JAER=IACTAEROP,1,-1

! Initialize local variables
      ZDP     = 0.0_JPRB
      ZAEREXT = 0.0_JPRB
      ZAERTAU(JL,JAER)=0.0_JPRB

      ZAERTAU(JL,JAER)= ZAERTAU(JL,JAER)+ PAEROD_HAT(JL,JJ)

!-- ITYP is the aerosol type 1:SeaSalt,   2:DesertDust,   3:OrganicMatter, 4: BlackCarbon,
!                            7:Strat background
!   IBIN is the bin index: 1-3:SS, 1-3:DD, 1-2:OM, 1-2:BC, 1:STRAT BACKGROUND

      ITYP=YAERO_DESC(JAER)%NTYP
      IBIN=YAERO_DESC(JAER)%NBIN
      IF(ITYP /= 0 .AND. IBIN == 0) THEN
        WRITE(NULERR,*) 'AT LEAST ONE BIN NEEDS TO BE ALLOCATED FOR TYPE ', ITYP,'!'
        WRITE(NULERR,*) 'PROGRAM WILL STOP'
        CALL ABOR1('') 
      ENDIF
      IF(ITYP /= ITYPE(JAER)) THEN
       WRITE(NULERR,*)'AOD_DUALCV_AD: AEROSOL TYPE IS NOT MATCHING', ITYP
       WRITE(NULERR,*) 'PROGRAM WILL STOP'
       CALL ABOR1('')
      ENDIF

      DO JK=KFLEVG,1,-1
      
         ZFAC=1.0_JPRB
         IF (ITYP == 1) THEN
          ZALF5(JAER)=ALF_SS(IRH(JL,JK),IWAVLINDEX,IBIN)
          ZFAC=RSS_RH80_MASSFAC
         ELSEIF (ITYP == 2) THEN
          ZALF5(JAER)=ALF_DD(IBIN,IWAVLINDEX)
         ELSEIF (ITYP == 3) THEN
          ZALF5(JAER)=ALF_OM(IRH(JL,JK),IWAVLINDEX)
         ELSEIF (ITYP == 4) THEN
          ZALF5(JAER)=ALF_BC(IWAVLINDEX)
         ELSEIF (ITYP == 5 .OR. ITYP == 10) THEN
          ZALF5(JAER)=ALF_SU(IRH(JL,JK),IWAVLINDEX)
          !-- SO2 does not contribute to optical depth, only SO4 does.
          IF (IBIN == 2) THEN
            ZALF5(JAER)=0._JPRB
          ENDIF
         ELSEIF (ITYP == 6) THEN ! Nitrate
            ZALF5(JAER)=ALF_NI(IRH(JL,JK),IWAVLINDEX,IBIN)
         ELSEIF (ITYP == 7) THEN ! Ammonium
            ZALF5(JAER)=ALF_AM(IRH(JL,JK),IWAVLINDEX)
         ELSEIF (ITYP == 8) THEN ! SOA
            ZALF5(JAER)=ALF_SOA(IRH(JL,JK),IWAVLINDEX,IBIN)
         ELSEIF (ITYP == 9) THEN ! volcanic ashes
!-- use dust for 0.9-20 um bin
            ZALF5(JAER)=ALF_DD(3,IWAVLINDEX)
         ENDIF

! convert units from m^2/g to m^2/kg
            ZALF5(JAER)=ZALF5(JAER)*1000._JPRB

            ZAEREXT5=ZAER5(JL,JK,JAER)*ZFAC*ZALF5(JAER)
 
            ZAEREXT= ZAEREXT + (ZDP5(JL,JK)/RG)* ZAERTAU(JL,JAER)
            ZDP    = ZDP + ZAEREXT5/RG * ZAERTAU (JL,JAER)  
! For now set ZDP to zero here and avoid increments in surface pressure
!         ZDP= 0.0_JPRB
   
            ZAER(JL,JK,JAER) = ZAER(JL,JK,JAER) + ZALF5(JAER)*ZFAC*ZAEREXT
            ZAEREXT= 0.0_JPRB

            YDGP_AD%PRESH(JL,JK,IH)   = YDGP_AD%PRESH(JL,JK,IH)  + ZDP
            YDGP_AD%PRESH(JL,JK-1,IH) = YDGP_AD%PRESH(JL,JK-1,IH)- ZDP
            ZDP =0.0_JPRB
  
      ENDDO  ! loop over vertical levels
      ZAERTAU(JL,JAER)=0.0_JPRB 
    ENDDO
    PAEROD_HAT(JL,JJ)=0.0_JPRB
  ENDDO
ENDDO

!! BEN the gradient of the total mass mixing ratio are computed from the gradient of
!  the various species using the mass fraction assumed constant 

ZAERLGMR(:,:)=0.0_JPRB
DO JLEN=1,KLEN
  DO JAER=JPACTAERLG,1,-1
    ZAERLGMR(JLEN,:) = ZAERLGMR(JLEN,:) + ZAER(JLEN,:,IACTAERLG(JAER))* ZFRAERLG5(JLEN,:,IACTAERLG(JAER))
    ZAER(JLEN,:,IACTAERLG(JAER))=0.0_JPRB
  ENDDO
ENDDO
ZAERSMMR(:,:)=0.0_JPRB
DO JLEN=1,KLEN
  DO JAER=JPACTAERSM,1,-1
    ZAERSMMR(JLEN,:) = ZAERSMMR(JLEN,:) + ZAER(JLEN,:,IACTAERSM(JAER))* ZFRAERSM5(JLEN,:,IACTAERSM(JAER))
    ZAER(JLEN,:,IACTAERSM(JAER))=0.0_JPRB
  ENDDO
ENDDO

! Put back gradient of cost function with respect to total aerosol mixing ratio
DO JLEN=1,KLEN
  DO JF=1,YDGP5%NGEMS
    IF(YDGP5%GEMS_IGRIB(JF) == 210047) THEN 
      YDGP_AD%GEMS(JLEN,1:KFLEVG,IH,JF) = YDGP_AD%GEMS(JLEN,1:KFLEVG,IH,JF) + ZAERSMMR(JLEN,:)
      ZAERSMMR(JLEN,:)=0.0_JPRB
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == 210048) THEN
      YDGP_AD%GEMS(JLEN,1:KFLEVG,IH,JF) = YDGP_AD%GEMS(JLEN,1:KFLEVG,IH,JF) + ZAERLGMR(JLEN,:)
      ZAERLGMR(JLEN,:)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

DEALLOCATE(ZAER5)
DEALLOCATE(ZAER)
DEALLOCATE(ITYPE)
DEALLOCATE(ZAERTAU)
DEALLOCATE(ZAERTAU5)
DEALLOCATE(ZFRAERSM5)
DEALLOCATE(ZFRAERLG5)
DEALLOCATE(ZALF5)
!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AOD_DUALCV_AD',1,ZHOOK_HANDLE)
END SUBROUTINE AOD_DUALCV_AD
