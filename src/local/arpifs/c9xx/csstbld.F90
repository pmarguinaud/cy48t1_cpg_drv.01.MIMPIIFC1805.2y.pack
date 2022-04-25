SUBROUTINE CSSTBLD(YDGEOMETRY)

!**** *CSSTBLD*

!     PURPOSE.
!     --------
!       Creation of an FA file containing a reference sea surface temperature.

!**   INTERFACE.
!     ----------
!       CALL CSSTBLD

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
!       F. TAILLEFER    APRIL 2010  (from INCLITC)

!     MODIFICATIONS.
!     --------------
!       Y. Michel    06/12 : add option for explicit perturbation
!                    06/13 : add LELAM case
!       K. Yessad (July 2014): Move some variables.
!       2019-04-11 Jean-Marcel Piriou and Adrien Napoly : read NETCDF SST file instead unformatted file (LLNETCDF).
!-------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMLUN   , ONLY : NULOUT, NULNAM, NINISH
USE YOMOPH0  , ONLY : CNMCA
USE YOMCT0   , ONLY : CNMEXP, LELAM, NCONF
USE YOMRIP0  , ONLY : NINDAT, NSSSSS
USE YOMCLTC  , ONLY : NLONTC, NLATTC, NLONTCLK, NLATTCLK, LAKEFIELD,&
                     & LPERTSST, NSEEDSST, ZSIGMA_A_SST, ZSIGMA_A_SST_TH,&
                     & ZLENGTH_SCALE_SST, ZINFL_PERT_SST, ZGAMMA_HP,     &
                     & LLDEBUGPERTSST
USE YOMCLI   , ONLY : YRCLI
USE NETCDF
USE LECECR_NETCDF_MOD
USE IOGRIDA_MOD, ONLY : RUNDEF
!-------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN)   :: YDGEOMETRY
INTEGER(KIND=JPIM) :: IDATEF(11),ILMO(12)

REAL(KIND=JPRB),ALLOCATABLE :: ZS(:),ZLSM(:),ZLSM2(:),ZS2(:)

CHARACTER :: CLNOMF*12, CLNOMF2*16
CHARACTER (LEN = 200) ::  CLDATE
INTEGER(KIND=JPIM) :: ILEN, INBARI, INBARP, INIMES, INIQ, INJQ, IUNIT,&
 & IREP, JI, JJ, IA1, IM1, IJ1, IA2, IM2, IJ2, IA3, IM3, IJ3, ID1, ID2, ICOUNT
INTEGER(KIND=JPIM) :: INCIDHN ! nc id: logical unit of NETCDF file.
INTEGER(KIND=JPIM) ::ISTATUS
REAL(KIND=JPRB) :: ZSEUIL,ZVALAND,ZVALMSK,ZLIMASK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
LOGICAL :: LLNETCDF
!-------------------------------------------------------------------------------
#include "faieno.h"
#include "posnam.intfb.h"
#include "eleci.intfb.h"
#include "inclitc.intfb.h"
#include "updcal.intfb.h"
#include "abor1.intfb.h"
#include "namcltc.nam.h"

!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CSSTBLD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NDGLG=>YDDIM%NDGLG, NDGUNG=>YDDIM%NDGUNG, NDGUXG=>YDDIM%NDGUXG, &
 & NDLON=>YDDIM%NDLON, NDLUNG=>YDDIM%NDLUNG, NDLUXG=>YDDIM%NDLUXG)
!-------------------------------------------------------------------------------
!     1. READ NAMELIST.
!        -------------
NLONTC=360
NLATTC=180
LAKEFIELD=.FALSE.
LLNETCDF=.FALSE.
IF (NCONF==933) LLNETCDF=.TRUE.
NLONTCLK=1800
NLATTCLK=900
! By default, do not perturb SST
LPERTSST=.FALSE.
ZINFL_PERT_SST=1.0_JPRB
NSEEDSST=0
! Constants adapted from from Donlon et.al., 2012
ZSIGMA_A_SST=0.60_JPRB
ZSIGMA_A_SST_TH=2.00_JPRB
! WARNING: This lengthscale is specific to the correlation model of gen_corr_pert
IF (LELAM) THEN
  ZLENGTH_SCALE_SST=200.0E3_JPRB
ELSE
  ZLENGTH_SCALE_SST=250.0E3_JPRB
ENDIF

CALL POSNAM (NULNAM,'NAMCLTC')
READ (NULNAM,NAMCLTC)

WRITE (UNIT=NULOUT,&
 & FMT='(/,'' input main SST grid size : NLONTC = '',I4,''   NLATTC = '',I4)')&
 & NLONTC,NLATTC
WRITE (NULOUT,*) '    '
WRITE (NULOUT,*) ' second field for lakes updating : LAKEFIELD=',LAKEFIELD
IF (LAKEFIELD) WRITE (UNIT=NULOUT,&
 & FMT='(/,'' input second SST grid (for lakes) : NLONTCLK = '',I4,&
 & ''    NLATTCLK = '',I4)') NLONTCLK,NLATTCLK
IF (LPERTSST) THEN
  WRITE (NULOUT,*) 'LPERTSST -> EXPLICIT PERTURBATION OF SST VALUES WITH SEED=',NSEEDSST
  WRITE (NULOUT,'(A,F12.4,A,F12.4,A,F12.4)') 'LPERTSST -> LENGTHSCALE=',ZLENGTH_SCALE_SST,&
 & ' SIGMA_A=',ZSIGMA_A_SST, ' INFLATION FACTOR=',ZINFL_PERT_SST
  IF (ZLENGTH_SCALE_SST<0.0) CALL ABOR1('CSSTBLD: LENGTHSCALE IS NEGATIVE')
  IF (ZSIGMA_A_SST<0.0)  CALL ABOR1('CSSTBLD: SIGMA_A IS NEGATIVE')
  IF (ZSIGMA_A_SST_TH<ZSIGMA_A_SST)&
    & CALL ABOR1('CSSTBLD: THRESHOLD ON SIGMA_A IS LOWER THAN SIGMA_A')
ENDIF
!-------------------------------------------------------------------------------
!     2. SET INITIAL VALUES.
!        ------------------

YRCLI%LGLOBE=.TRUE.

IREP=0
INBARI=0
INIMES=1
INBARP=1
IUNIT=15
CLNOMF='ICMSH'//CNMEXP(1:4)//'SST'
CLNOMF2='ICMSH'//CNMEXP(1:4)//'INIT'

IF (LELAM) THEN
  INJQ=NDGUXG-NDGUNG+1
  INIQ=NDLUXG-NDLUNG+1
  ILEN=INIQ*INJQ
ELSE
  ILEN=NDGLG*NDLON
ENDIF

ALLOCATE(ZLSM(ILEN))
ALLOCATE(ZLSM2(ILEN))
ALLOCATE(ZS(ILEN))
IF (LAKEFIELD) ALLOCATE(ZS2(ILEN))

ZLIMASK=0.5_JPRB
ZVALMSK=271.15_JPRB
ZVALAND=288.15_JPRB
ZSEUIL=350._JPRB

!-------------------------------------------------------------------------------
!     3. OUTPUT GRID LAND/SEA MASK READING.
!        ---------------------------------

! This field is read in the initial conditions FA file (climatology)

CALL FAITOU (IREP,NINISH,.TRUE.,CLNOMF2,'UNKNOWN',.TRUE.,.TRUE.,&
 & INIMES,INBARP,INBARI,'CADRE LECTURE   ')  
CALL LFIMST (IREP,NINISH,.FALSE.)
IF (LELAM) THEN
  CALL ELECI(YDGEOMETRY%YRDIM,IREP,NINISH,'SURF',0,'IND.TERREMER',ZLSM,NULOUT)
ELSE
  CALL FACILE (IREP,NINISH,'SURF',0,'IND.TERREMER',ZLSM,.FALSE.)
ENDIF
CALL FAIRME (IREP,NINISH,'UNKNOWN')

ZLSM(:)=1.0_JPRB-ZLSM(:)


!Add sea information from PGD field to the LSM
CALL FAITOU (IREP,NINISH,.TRUE.,'Const.Clim.sfx','UNKNOWN',.TRUE.,.TRUE.,&
             & INIMES,INBARP,INBARI,'CADRE LECTURE   ')
CALL FACILE (IREP,NINISH,'SFX',0,'.FRAC_SEA',ZLSM2,.FALSE.)
CALL FAIRME (IREP,NINISH,'UNKNOWN')
DO JJ=1,ILEN
  IF (ZLSM2(JJ)>0) ZLSM(JJ)=1.0   !0.=land ans 1.0=sea (opposite of the usual ITM)
ENDDO



!-------------------------------------------------------------------------------
!     4. INPUT DATA TREATMENT.
!        --------------------
IF (LLNETCDF) THEN
  ! Open NETCDF file.
  ISTATUS=NF90_OPEN('OSTIA.nc',NF90_NOWRITE,INCIDHN)
  IF (ISTATUS /= NF90_NOERR) CALL HANDLE_ERR(ISTATUS)
  ! Date.
  ISTATUS=NF90_GET_ATT(INCIDHN,NF90_GLOBAL,"start_date", CLDATE)
  IF (ISTATUS /= NF90_NOERR) CALL HANDLE_ERR(ISTATUS)
  READ(CLDATE,FMT='(I4,1X,I2,1X,I2)') IA1,IM1,IJ1
ELSE
!     4.1  Date control
  OPEN(UNIT=3,FILE='sst.reference',FORM='UNFORMATTED')
  READ(3) IA1,IM1,IJ1
  CLOSE(3)

  IF (LAKEFIELD) THEN
    OPEN(UNIT=4,FILE='sst.lake',FORM='UNFORMATTED')
    READ(4) IA2,IM2,IJ2
    CLOSE(4)
    IF (IA2/=IA1.OR.IM2/=IM1.OR.IJ2/=IJ1) THEN
      ID1=IA1*10000+IM1*100+IJ1
      ID2=IA2*10000+IM2*100+IJ2
      WRITE(NULOUT,*) '   '
      WRITE(NULOUT,*) ' dates file 1 and 2 are different :  ',ID1,ID2
  !    10 days old authorized for the lake field
      CALL UPDCAL(IJ1,IM1,IA1,-10,IJ3,IM3,IA3,ILMO,99)
      ID1=IA3*10000+IM3*100+IJ3
      IF (ID1>ID2) THEN
        WRITE(NULOUT,*) '    file 2 is too old  --> no use of the second field'
        LAKEFIELD=.FALSE.
      ENDIF
    ENDIF
  ENDIF
ENDIF
!     4.2  First field interpolation

IREP=1
CALL INCLITC(YDGEOMETRY,IREP,NLONTC,NLATTC,ILEN,ZLSM,ZS)

!     4.3  Second field interpolation and final SST updating

IF (LAKEFIELD) THEN
  ICOUNT=0
  IREP=2
  CALL INCLITC(YDGEOMETRY,IREP,NLONTCLK,NLATTCLK,ILEN,ZLSM,ZS2)
  DO JJ = 1,ILEN
    IF (ZS(JJ) > ZSEUIL .AND. ZS2(JJ) < ZSEUIL ) THEN
      ICOUNT=ICOUNT+1
      ZS(JJ)=ZS2(JJ)
    ENDIF
  ENDDO
  WRITE(NULOUT,*) '   '
  WRITE(NULOUT,*) '    final SST updating (borders/lakes) on ',ICOUNT,'points'
ENDIF

IF (ALLOCATED(ZS2)) DEALLOCATE(ZS2)

!     4.4  Final field last checking

IA2=0
DO JI = 1,ILEN
  IF (ZLSM(JI) < ZLIMASK) THEN
    ZS(JI)=RUNDEF
  ELSE
    IF (ZS(JI) < 0. .OR. ZS(JI) > ZSEUIL) THEN
      IA2=IA2+1
      ZS(JI)=RUNDEF
    ENDIF
   ENDIF
ENDDO

IF (ALLOCATED(ZLSM)) DEALLOCATE(ZLSM)
IF (ALLOCATED(ZLSM2)) DEALLOCATE(ZLSM2)

WRITE(NULOUT,*) '  '
WRITE(NULOUT,*) '    il y a ',IA2,'points indefinis sur mer.Ils sont mis a UNDEF '

!-------------------------------------------------------------------------------
!     5. OUTPUT FILE WRITING.
!        -------------------

!     5.1 File creation

CALL FAITOU (IREP,IUNIT,.TRUE.,CLNOMF,'NEW',.TRUE.,.TRUE.,INIMES,&
 & INBARP,INBARI,CNMCA)  
CALL LFIMST (IREP,IUNIT,.FALSE.)
CALL FAUTIF (IREP,IUNIT,'REFERENCE SST')

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

IF (LELAM) THEN
  IA1=0
  IM1=0
  ALLOCATE(ZS2(NDLON*NDGLG))
  ZS2(:)=0.
  DO JJ=1,NDGUXG
    ZSEUIL=0._JPRB
    DO JI=1,NDLUXG
      IA1=IA1+1
      IM1=IM1+1
      ZS2(IA1)=ZS(IM1)
      ZSEUIL=ZSEUIL+ZS(IM1)
    ENDDO
    ZSEUIL=ZSEUIL/FLOAT(NDLUXG)
    DO JI=NDLUXG+1,NDLON
      IA1=IA1+1
      ZS2(IA1)=ZSEUIL
    ENDDO
  ENDDO
  DO JJ=NDGUXG+1,NDGLG
    DO JI=1,NDLON
      IA1=IA1+1
      ZS2(IA1)=ZSEUIL
    ENDDO
  ENDDO
  WHERE(ZS2>1000. .and. ZS2/=RUNDEF)ZS2=RUNDEF ! valeurs qui sont devenues plus grandes que RUNDEF
  CALL FAIENO (IREP,IUNIT,'SURF',0,'SST.CLIM.',ZS2, &
              & .FALSE., LDUNDF=.TRUE., PUNDF=RUNDEF)
ELSE
  CALL FAIENO (IREP,IUNIT,'SURF',0,'SST.CLIM.',ZS, &
              & .FALSE., LDUNDF=.TRUE., PUNDF=RUNDEF)
ENDIF

CALL FAIRME (IREP,IUNIT,'UNKNOWN')

WRITE(NULOUT,*) ' END of configuration ',NCONF
WRITE(NULOUT,*) '    '

!-------------------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CSSTBLD',1,ZHOOK_HANDLE)
END SUBROUTINE CSSTBLD
