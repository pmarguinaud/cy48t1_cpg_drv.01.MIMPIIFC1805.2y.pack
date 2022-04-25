SUBROUTINE HSEAICE(YDGEOMETRY)

!**** *HSEAICE*

!     PURPOSE.
!     --------
!       Creation of an FA file containing a reference sea surface temperature.

!**   INTERFACE.
!     ----------
!       CALL HSEAICE

!     METHOD.
!     -------
!       Reading of an external SST field and interpolation over the required AAA domain
!         (the input field may be a NESDIS or OSTIA daily analysis).
!       Creation and filling of a FA output "clim" SST file.

!     EXTERNALS.
!     ----------
!       FAITOU FANDAR FAUTIF FACILE FAIENC FAIRME FACIES LFIMST ELECI INCLITC

!     AUTHORS.
!     --------
!       A. Napoly 08-2020

!     MODIFICATIONS.
!     --------------
!
!-------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMLUN   , ONLY : NULOUT, NINISH
USE YOMOPH0  , ONLY : CNMCA
USE YOMCT0   , ONLY : CNMEXP
USE YOMRIP0  , ONLY : NINDAT, NSSSSS
USE YOMCLTC  , ONLY : NLONTC, NLATTC,&
                     & LPERTSST, NSEEDSST, ZSIGMA_A_SST, ZSIGMA_A_SST_TH,&
                     & ZINFL_PERT_SST
USE YOMCLI   , ONLY : YRCLI
USE NETCDF
USE LECECR_NETCDF_MOD
!-------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN)   :: YDGEOMETRY
INTEGER(KIND=JPIM) :: IDATEF(11)

REAL(KIND=JPRB),ALLOCATABLE :: ZS(:),ZLSM(:)

CHARACTER :: CLNOMF*13
CHARACTER (LEN = 200) ::  CLDATE
INTEGER(KIND=JPIM) :: ILEN, INBARI, INBARP, INIMES, IUNIT,&
 & IREP, JJ, IA1, IM1, IJ1
INTEGER(KIND=JPIM) :: INCIDHN ! nc id: logical unit of NETCDF file.
INTEGER(KIND=JPIM) ::ISTATUS
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------

#include "inclitc_hice.intfb.h"

!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('HSEAICE',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NDGLG=>YDDIM%NDGLG, NDGUNG=>YDDIM%NDGUNG, NDGUXG=>YDDIM%NDGUXG, &
 & NDLON=>YDDIM%NDLON, NDLUNG=>YDDIM%NDLUNG, NDLUXG=>YDDIM%NDLUXG)
!-------------------------------------------------------------------------------
!     1. READ NAMELIST.
!        -------------

NLONTC=4320!360 ANTMPTEST
NLATTC=2041!180 ANTMPTEST
!-------------------------------------------------------------------------------
!     2. SET INITIAL VALUES.
!        ------------------

YRCLI%LGLOBE=.TRUE.

IREP=0
INBARI=0
INIMES=1
INBARP=1
IUNIT=15
CLNOMF='ICMSH'//CNMEXP(1:4)//'HICE'

ILEN=NDGLG*NDLON

ALLOCATE(ZLSM(ILEN))
ALLOCATE(ZS(ILEN))

!     3. OUTPUT GRID LAND/SEA MASK READING.
!        ---------------------------------
CALL FAITOU (IREP,NINISH,.TRUE.,'Const.Clim.sfx','UNKNOWN',.TRUE.,.TRUE.,&
 & INIMES,INBARP,INBARI,'CADRE LECTURE   ')
CALL FACILE (IREP,NINISH,'SFX',0,'.FRAC_SEA',ZLSM,.FALSE.)
CALL FAIRME (IREP,NINISH,'UNKNOWN')
DO JJ=1,ILEN
  IF (ZLSM(JJ)>0) ZLSM(JJ)=1.0   !0.=land ans 1.0=sea (opposite of the usual ITM)
ENDDO


!-------------------------------------------------------------------------------
!     4. INPUT DATA TREATMENT.
!        --------------------
! Open NETCDF file.
ISTATUS=NF90_OPEN('OCEANFILE.nc',NF90_NOWRITE,INCIDHN)
IF (ISTATUS /= NF90_NOERR) THEN
    !CALL HANDLE_ERR(ISTATUS)
ELSE
    ! Date.
    ISTATUS=NF90_GET_ATT(INCIDHN,NF90_GLOBAL,"field_date", CLDATE)
    IF (ISTATUS /= NF90_NOERR) CALL HANDLE_ERR(ISTATUS)
    READ(CLDATE,FMT='(I4,1X,I2,1X,I2)') IA1,IM1,IJ1

    !     4.2  First field interpolation

    CALL INCLITC_HICE(YDGEOMETRY,NLONTC,NLATTC,ILEN,ZLSM,ZS)

    IF (ALLOCATED(ZLSM)) DEALLOCATE(ZLSM)

    !-------------------------------------------------------------------------------
    !     5. OUTPUT FILE WRITING.
    !        -------------------
    
    !     5.1 File creation
    CALL FAITOU (IREP,IUNIT,.TRUE.,CLNOMF,'NEW',.TRUE.,.TRUE.,INIMES,INBARP,INBARI,CNMCA)  
    CALL LFIMST (IREP,IUNIT,.FALSE.)
    CALL FAUTIF (IREP,IUNIT,'SEA ICE DEPTH')
    
    !     5.2 Date setting
    WRITE(NULOUT,'(/,''  date of the model: '',I8,I8)') NINDAT,NSSSSS
    WRITE(NULOUT,'(''  date of the data:  '',I4,I2.2,I2.2,/)') IA1,IM1,IJ1
    IDATEF(1)= IA1
    IDATEF(2)= IM1
    IDATEF(3)= IJ1
    IDATEF(4)=0
    IDATEF(5)=0
    IDATEF(6)=1
    IDATEF(7:11)=0
    CALL FANDAR (IREP,IUNIT,IDATEF)
    
    !     5.3 Field writing
    CALL FAIENC (IREP,IUNIT,'SURF',0,'SEAICE.THICK',ZS,.FALSE.)
    CALL FAIRME (IREP,IUNIT,'UNKNOWN')
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('HSEAICE',1,ZHOOK_HANDLE)
END SUBROUTINE HSEAICE
