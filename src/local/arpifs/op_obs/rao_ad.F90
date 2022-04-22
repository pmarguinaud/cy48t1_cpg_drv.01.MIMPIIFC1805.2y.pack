SUBROUTINE RAO_AD&
 & (KLEVG, KDLEN  , KLEN  , KMXCOUNT, &
 & YDGP_AD,&
 & YDGP5,&
 & PAEROD_HAT)

!*** *RAO_AD* Adjoint Operator routine for aerosol optical depth observations.

!**   INTERFACE.
!     ----------

!     WHERE KDLEN      = First dimension: Max number of observations, this set (INPUT)
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
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!        S. Remy     16-Sep-2016    Add nitrate and ammonium
!        S. Remy     13-Nov-2017    Add SOA


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
USE YOEAERATM     , ONLY : YREAERATM
USE YOEAEROP      , ONLY : ALF_SU, ALF_OM, ALF_DD, ALF_SS, ALF_BC, ALF_NI, ALF_AM, ALF_SOA
USE YOEAERSNK     , ONLY : YREAERSNK
USE GOM_PLUS      , ONLY : TYPE_GOM_PLUS, IH
USE YOM_GRIB_CODES, ONLY : NGRBAERMR01,NGRBAERMR04, &
 & NGRBAERMR07, NGRBAERMR08, NGRBAERMR09, NGRBAERMR10, NGRBAERMR11, &
 & NGRBAERMR16, NGRBAERMR18, NGRBAERMR19, NGRBAERMR20, NGRBAERMR21, NGRBAERSM
USE YOEPHLI       , ONLY : YREPHLI

IMPLICIT NONE

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KLEVG, KDLEN, KLEN, KMXCOUNT

TYPE (TYPE_GOM_PLUS),INTENT(INOUT) :: YDGP_AD
TYPE (TYPE_GOM_PLUS),INTENT(IN)    :: YDGP5

REAL(KIND=JPRB)     ,INTENT(INOUT) :: PAEROD_HAT(KDLEN,KMXCOUNT) ! AEROSOL OPT DEPTH GRADIENT 

!*       0.5   LOCAL VARIABLES
!              ---------------
!INTEGER(KIND=JPIM), PARAMETER :: IWAVLMX=7 ! Max number of MODIS channels
!INTEGER(KIND=JPIM), DIMENSION(IWAVLMX) :: IWAVL=(/6, 9, 11, 14, 17, 18, 19/) ! Index for MODIS channels
INTEGER(KIND=JPIM) :: IAER, ITYP, IWAVLINDEX
INTEGER(KIND=JPIM) :: JAER, JK, JL, JJ, JF, JP, JLEN, JLEV, JTAB
INTEGER(KIND=JPIM) :: IFLAG 
INTEGER(KIND=JPIM), PARAMETER  :: JTABMX=12
INTEGER(KIND=JPIM), ALLOCATABLE :: ITYPE(:) !index to determine optically active aerosol species
INTEGER(KIND=JPIM), DIMENSION(KLEN,KLEVG) :: IRH

REAL(KIND=JPRB), DIMENSION(KLEN,KLEVG)   :: ZQSAT5, ZT5, ZPF5
REAL(KIND=JPRB), DIMENSION(KLEN,KLEVG)   :: ZT, ZPF

REAL(KIND=JPRB), ALLOCATABLE   :: ZAERTAU5(:,:)
REAL(KIND=JPRB), ALLOCATABLE   :: ZAERTAU(:,:)
REAL(KIND=JPRB), ALLOCATABLE   :: ZAERF5(:,:,:) ! aer mixing ratios from traj (interpolated to obs locations)
REAL(KIND=JPRB), ALLOCATABLE   :: ZAERF(:,:,:) ! pert aer mixing ratios from traj (interpolated to obs locations)

REAL(KIND=JPRB), DIMENSION(KLEN,KLEVG) :: ZAERFT5 !
REAL(KIND=JPRB), DIMENSION(KLEN,KLEVG) :: ZAERFMR5 ! total aer mixing ratios from traj
REAL(KIND=JPRB), DIMENSION(KLEN,KLEVG) :: ZAERFMR ! total aer mixing ratios pert

REAL(KIND=JPRB), ALLOCATABLE :: ZFRACAER5(:,:,:) ! fraction of mixing ratios wrt total mass

REAL(KIND=JPRB), DIMENSION(KDLEN,KMXCOUNT) :: ZAEROD5 ! aer optical depth
REAL(KIND=JPRB), ALLOCATABLE :: ZALF5(:)
REAL(KIND=JPRB)   :: ZRHCL5, ZAEREXT5
REAL(KIND=JPRB), DIMENSION(KLEN,KLEVG):: ZDP5 ! pressure difference from trajectory
REAL(KIND=JPRB)   :: ZAEREXT, ZDP, ZFAC

INTEGER(KIND=JPIM) :: INACTAERO

REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL :: LLPHYLIN

#include "satur.intfb.h"
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RAO_AD',0,ZHOOK_HANDLE)
ASSOCIATE(&
 & RSS_RH80_MASSFAC=>YREAERATM%RSS_RH80_MASSFAC, &
 & RRHTAB=>YREAERSNK%RRHTAB, &
 & NMAXTAER=>YREAERATM%NMAXTAER, NTYPAER=>YREAERATM%NTYPAER)

LLPHYLIN = YREPHLI%LPHYLIN

ALLOCATE(ZAERF5(KLEN,KLEVG,YDGP5%NGEMS))
ALLOCATE(ZAERF(KLEN,KLEVG,YDGP5%NGEMS))
ALLOCATE(ITYPE(YDGP5%NGEMS))
ALLOCATE(ZAERTAU(KLEN,YDGP5%NGEMS))
ALLOCATE(ZAERTAU5(KLEN,YDGP5%NGEMS))
ALLOCATE(ZFRACAER5(KLEN,KLEVG,YDGP5%NGEMS))
ALLOCATE(ZALF5(YDGP5%NGEMS))

ZAERF5(:,:,:)=0.0_JPRB
ZAERF(:,:,:)=0.0_JPRB
ITYPE(:)=0

! Extract aerosol mixing ratios from GEMS profiles structure
DO JLEN=1,KLEN
  IAER=0
  DO JF=1,YDGP5%NGEMS
    JP=0
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR01) JP=1 ! SS 1 bin
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR04) JP=2 ! DD 1 bin 
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR07 .OR. & 
     &   YDGP5%GEMS_IGRIB(JF) == NGRBAERMR08) JP=3 ! ORGANIC MATTER - 2 bins
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR09 .OR. &
     &   YDGP5%GEMS_IGRIB(JF) == NGRBAERMR10) JP=4 ! Black Carbon - 2 bins
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR11) JP=5 ! SO4
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR16) JP=6 ! NO3
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR18) JP=7 ! NH4
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERMR19 .OR. & 
    &   YDGP5%GEMS_IGRIB(JF) == NGRBAERMR20 .OR. &
     &   YDGP5%GEMS_IGRIB(JF) == NGRBAERMR21) JP=8 ! SOA - 3 bins

    IF(JP > 0) THEN
      IAER=IAER+1
      ZAERF5(JLEN,:,IAER)=YDGP5%GEMS(JLEN,1:KLEVG,IH,JF)
      ITYPE(IAER)=JP
    ENDIF
    IF(YDGP5%GEMS_IGRIB(JF) == NGRBAERSM) THEN ! Fine mode aerosol mixing ratio
      ZAERFMR5(JLEN,:)=YDGP5%GEMS(JLEN,1:KLEVG,IH,JF)
    ENDIF
  ENDDO
ENDDO
INACTAERO=IAER

! Compute total mixing ratio by summing single contributions
    ZAERFT5(:,:) = 0.0_JPRB
    DO JAER=1,INACTAERO
      ZAERFT5(:,:)=ZAERFT5(:,:)+ ZAERF5(:,:,JAER)
    ENDDO

! BEN the perturbations on the various species are computed from the perturbation
! on the total mass mixing ratio using the mass fraction assumed constant

ZFRACAER5(:,:,:) = 0.0_JPRB
DO JLEN=1,KLEN
  DO JAER=1,INACTAERO
    DO JLEV=1,KLEVG
      IF(ZAERFT5(JLEN,JLEV) /= 0.0_JPRB) THEN
        ZFRACAER5(JLEN,JLEV,JAER) = ZAERF5(JLEN,JLEV,JAER)/ZAERFT5(JLEN,JLEV)
      ENDIF
    ENDDO
  ENDDO
ENDDO

!! Compute fraction of given species (and bin) wrt total mass mixing ratio

!DO JLEN=1,KLEN
!  DO JAER=1,INACTAERO
!    DO JLEV=1,KLEVG
!      IF(ZAERFMR5(JLEN,JLEV) /= 0.0_JPRB) THEN
!        ZFRACAER5(JLEN,JLEV,JAER) = ZAERF5(JLEN,JLEV,JAER)/ZAERFMR5(JLEN,JLEV)
!      ELSE
!        ZFRACAER5(JLEN,JLEV,JAER) = 0.0_JPRB
!      ENDIF
!    ENDDO
!  ENDDO
!ENDDO

! Compute relative humidity at observation points, on model levels
ZT5(1:KLEN,:) = YDGP5%TF(1:KLEN,1:KLEVG,IH) ! trajectory
ZPF5(1:KLEN,:)= YDGP5%PRESF(1:KLEN,:,IH)     ! trajectory
ZT(1:KLEN,:) = 0.0_JPRB       ! perturb
ZPF(1:KLEN,:)= 0.0_JPRB       ! perturb
!CALL SATURAD (1 , KLEN , KLEN  , 1 , KLEVG,&
!  & ZPF5, ZT5 , ZQSAT5, &
!  & ZPF,  ZT ,  ZQSAT)

IFLAG=2
CALL SATUR (YRTHF, YRCST, 1 , KLEN , KLEN  , 1 , KLEVG, LLPHYLIN, &
  & ZPF5, ZT5 , ZQSAT5 , IFLAG)


!PPRESF(1:KLEN,:)=PPRESF(1:KLEN,:)+ZPF(1:KLEN,:)     ! perturb in pressure - do not include
!PTF(1:KLEN,:)= PTF(1:KLEN,:) + ZT(1:KLEN,:)         ! perturb in temperature - do not include 
!ZT(1:KLEN,:) = 0.0_JPRB       
!ZPF(1:KLEN,:)= 0.0_JPRB   

!-- define RH index from "clear-sky" relative humidity
IRH(:,:)=1
DO JK=1,KLEVG
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
   DO JAER=1,INACTAERO

      ZAERTAU5(JL,JAER)= 0.0_JPRB

      ITYP=ITYPE(JAER)

      DO JK=1,KLEVG
         ZFAC=1.0_JPRB
      
         IF (ITYP == 1) THEN
          ZALF5(JAER)=ALF_SS(IRH(JL,JK),IWAVLINDEX,1)
          ZFAC=RSS_RH80_MASSFAC
         ELSEIF (ITYP == 2) THEN
          ZALF5(JAER)=ALF_DD(1,IWAVLINDEX)
         ELSEIF (ITYP == 3) THEN
          ZALF5(JAER)=ALF_OM(IRH(JL,JK),IWAVLINDEX)
         ELSEIF (ITYP == 4) THEN
          ZALF5(JAER)=ALF_BC(IWAVLINDEX)
         ELSEIF (ITYP == 5 .OR. ITYP == 10) THEN
          ZALF5(JAER)=ALF_SU(IRH(JL,JK),IWAVLINDEX)
         ELSEIF (ITYP == 6) THEN ! Nitrate
            ZALF5(JAER)=ALF_NI(IRH(JL,JK),IWAVLINDEX,1)
         ELSEIF (ITYP == 7) THEN ! Ammonium
            ZALF5(JAER)=ALF_AM(IRH(JL,JK),IWAVLINDEX)
         ELSEIF (ITYP == 8) THEN ! SOA
            ZALF5(JAER)=ALF_SOA(IRH(JL,JK),IWAVLINDEX,1)
         ELSEIF (ITYP == 9) THEN ! volcanic ashes
!-- use dust for 0.9-20 um bin
            ZALF5(JAER)=ALF_DD(3,IWAVLINDEX)
         ENDIF
! convert units from m^2/g to m^2/kg
         ZALF5(JAER)=ZALF5(JAER)*1000._JPRB

         ZAEREXT5=ZAERF5(JL,JK,JAER)*ZFAC*ZALF5(JAER)
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
    DO JAER=INACTAERO,1,-1

! Initialize local variables
      ZDP     = 0.0_JPRB
      ZAEREXT = 0.0_JPRB
      ZAERTAU(JL,JAER)=0.0_JPRB

      ZAERTAU(JL,JAER)= ZAERTAU(JL,JAER)+ PAEROD_HAT(JL,JJ)

!-- ITYP is the aerosol type 1:SeaSalt,   2:DesertDust,   3:OrganicMatter, 4: BlackCarbon, 5: Sulphate, 6: Nitrate
!                            7:Ammonium, 8: SOA

      ITYP=ITYPE(JAER)

      DO JK=KLEVG,1,-1
         ZFAC=1.0_JPRB
      
         IF (ITYP == 1) THEN
          ZALF5(JAER)=ALF_SS(IRH(JL,JK),IWAVLINDEX,1)
          ZFAC=RSS_RH80_MASSFAC
         ELSEIF (ITYP == 2) THEN
          ZALF5(JAER)=ALF_DD(1,IWAVLINDEX)
         ELSEIF (ITYP == 3) THEN
          ZALF5(JAER)=ALF_OM(IRH(JL,JK),IWAVLINDEX)
         ELSEIF (ITYP == 4) THEN
          ZALF5(JAER)=ALF_BC(IWAVLINDEX)
         ELSEIF (ITYP == 5 .OR. ITYP == 10) THEN
          ZALF5(JAER)=ALF_SU(IRH(JL,JK),IWAVLINDEX)
         ELSEIF (ITYP == 6) THEN ! Nitrate
            ZALF5(JAER)=ALF_NI(IRH(JL,JK),IWAVLINDEX,1)
         ELSEIF (ITYP == 7) THEN ! Ammonium
            ZALF5(JAER)=ALF_AM(IRH(JL,JK),IWAVLINDEX)
         ELSEIF (ITYP == 8) THEN ! Ammonium
            ZALF5(JAER)=ALF_SOA(IRH(JL,JK),IWAVLINDEX,1)
         ELSEIF (ITYP == 9) THEN ! volcanic ashes
!-- use dust for 0.9-20 um bin
            ZALF5(JAER)=ALF_DD(3,IWAVLINDEX)
         ENDIF
! convert units from m^2/g to m^2/kg
            ZALF5(JAER)=ZALF5(JAER)*1000._JPRB

            ZAEREXT5=ZAERF5(JL,JK,JAER)*ZFAC*ZALF5(JAER)
 
            ZAEREXT= ZAEREXT + (ZDP5(JL,JK)/RG)* ZAERTAU(JL,JAER)
            ZDP    = ZDP + ZAEREXT5/RG * ZAERTAU (JL,JAER)  
! For now set ZDP to zero here and avoid increments in surface pressure
!         ZDP= 0.0_JPRB
   
            ZAERF(JL,JK,JAER) = ZAERF(JL,JK,JAER) + ZALF5(JAER)*ZFAC*ZAEREXT
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

ZAERFMR(:,:)=0.0_JPRB
DO JLEN=1,KLEN
  DO JAER=INACTAERO,1,-1
    ZAERFMR(JLEN,:) = ZAERFMR(JLEN,:) + ZAERF(JLEN,:,JAER)* ZFRACAER5(JLEN,:,JAER)
    ZAERF(JLEN,:,JAER)=0.0_JPRB
  ENDDO
ENDDO
! Put back gradient of cost function with respect to total aerosol mixing ratio
DO JLEN=1,KLEN
  DO JF=1,YDGP5%NGEMS
    IF(YDGP5%GEMS_IGRIB(JF) == 210047) THEN 
      YDGP_AD%GEMS(JLEN,1:KLEVG,IH,JF) = YDGP_AD%GEMS(JLEN,1:KLEVG,IH,JF) + ZAERFMR(JLEN,:)
      ZAERFMR(JLEN,:)=0.0_JPRB
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
IF (LHOOK) CALL DR_HOOK('RAO_AD',1,ZHOOK_HANDLE)
END SUBROUTINE RAO_AD
