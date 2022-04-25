SUBROUTINE INCLITC(YDGEOMETRY,KF,KLON,KLAT,KLEN,PLSM,PS)

!**** *INCLITC* ***

!     PURPOSE.
!     --------
!       Interpolation of a field from a global lat/lon grid to an Arp/Ald/Aro grid.

!**   INTERFACE.
!     ----------
!       CALL INCLITC(KF,KLON,KLAT,KLEN,PLSM,PS)
!              KF   (IN) : 1 for SST, 2 for lakes
!              KLON (IN) : number of longitudes of the input grid
!              KLAT (IN) : number of latitudes of the input grid
!              KLEN (IN) : output grid vector length
!              PLSM (IN) : output grid land/sea mask
!              PS  (OUT) : input SST interpolated over the output grid

!     METHOD.
!     -------
!       The field is read on unit 3 (prescribed format).
!       It's interpolated on the gaussian grid for ARPEGE or on the specified ALADIN grid.

!     EXTERNALS.
!     ----------
!       INTER0 GEO923 INTER6 EINTER0 EINTER6 EGEO923

!     AUTHORS.
!     --------
!       F. TAILLEFER    JULY 99

!     MODIFICATIONS.
!     --------------
!       F. Taillefer 11/01 : moving input grid dimensionning from parameter to namelist
!       F. Taillefer 02/04 : changing INTER2 to INTER6 (use of land/sea mask)
!       M.Hamrud   01-Oct-2003 CY28 Cleaning
!       F. Taillefer 01/08 : add LAM case
!       F. Taillefer 04/10 : add lakes case (read from a second finer field) / CSSTBLD
!       T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!       Y. Michel    06/12 : add option for explicit perturbation
!                    06/13 : add LELAM case
!       2019-04-11 Jean-Marcel Piriou and Adrien Napoly : read NETCDF SST file instead of unformatted file (LLNETCDF).
!-------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMLUN   , ONLY : NULOUT
USE YOMCT0   , ONLY : LELAM, NCONF
USE YOMCLI   , ONLY : YRCLI
USE YOMCST   , ONLY : RTT
USE NETCDF
USE LECECR_NETCDF_MOD

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEN,KLON,KLAT,KF
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLEN)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PS(KLEN)

! JPBX/JPBY : number of extra longitudes/latitudes on each side of the data
!             grid required by interpolation (4 or 12 points -> 2/2)

INTEGER(KIND=JPIM) :: JPBX,JPBY
PARAMETER(JPBX=2,JPBY=2)

REAL(KIND=JPRB) :: ZBO(YDGEOMETRY%YRDIM%NDGLG+1),ZLONG(YDGEOMETRY%YRDIM%NDGLG),ZMU(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON)

REAL(KIND=JPRB),ALLOCATABLE :: ZFLD(:,:),ZITM(:,:),ZSLA(:),ZSLO(:,:),ZCLO(:,:),ZS2(:)
REAL(KIND=JPRB),ALLOCATABLE :: ZFLDTMP(:,:)

INTEGER(KIND=JPIM) :: INTLO, INTLA
INTEGER(KIND=JPIM),ALLOCATABLE :: ILIND(:)

INTEGER(KIND=JPIM) :: ICO,INIQ,INJQ,INS,IPNIV,ITOGOS,JI,JJ,INPNPP,IL1,IL2,ILIMGP,ISTATUS
INTEGER(KIND=JPIM) :: INCIDHN ! nc id: logical unit of NETCDF file.
INTEGER(KIND=JPIM) :: J2,J1,ILON,ILAT

REAL(KIND=JPRB) :: ZLIMASK, ZSEUIL
REAL(KIND=JPRB) :: ZHOOK_HANDLE

CHARACTER :: CLNOMF*20

REAL(KIND=JPRB), DIMENSION(:,:,:), POINTER :: ZSST
LOGICAL :: LLNETCDF

#include "geo923.intfb.h"
#include "inter0.intfb.h"
#include "inter6.intfb.h"
#include "abor1.intfb.h"
#include "egeo923.intfb.h"
#include "einter0.intfb.h"
#include "einter6.intfb.h"
!-------------------------------------------------------------------------------
!     1. SET INITIAL VALUES.
!        ------------------
LLNETCDF=.FALSE.
IF (NCONF==933) LLNETCDF=.TRUE.

IF (LHOOK) CALL DR_HOOK('INCLITC',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDEGEO=>YDGEOMETRY%YREGEO)
ASSOCIATE(NDGLG=>YDDIM%NDGLG, NDGUNG=>YDDIM%NDGUNG, NDGUXG=>YDDIM%NDGUXG, &
 & NDLON=>YDDIM%NDLON, NDLUNG=>YDDIM%NDLUNG, NDLUXG=>YDDIM%NDLUXG, &
 & EDELX=>YDEGEO%EDELX, EDELY=>YDEGEO%EDELY, &
 & NLOENG=>YDGEM%NLOENG, RSTRET=>YDGEM%RSTRET)

IF (KF==1) THEN
  CLNOMF='sst.reference'
ELSEIF(KF==2) THEN
  CLNOMF='sst.lake'
ENDIF

WRITE(NULOUT,*) '   '
WRITE(NULOUT,*) ' call to INCLITC - field',KF,'- file : ',CLNOMF

YRCLI%NGLOBX=KLON
YRCLI%NGLOBY=KLAT

YRCLI%NDATX=KLON
YRCLI%NDATY=KLAT

IF (LELAM) THEN
  CALL EINTER0(YDGEOMETRY,KLON,KLAT,INPNPP)
  INS=-1
  IL1=KLAT+JPBY
  IL2=1+JPBY
  INJQ=NDGUXG-NDGUNG+1
  INIQ=NDLUXG-NDLUNG+1
ELSE
  CALL INTER0 (NDGLG,NDLON,KLAT,KLON,RSTRET,INPNPP)
  INS=1
  IL1=1+JPBY
  IL2=KLAT+JPBY
  INJQ=NDGLG*INPNPP
  INIQ=NDLON*INPNPP
ENDIF

INTLO=KLON+2*JPBX
INTLA=KLAT+2*JPBY

!         Allocate temporary space

ALLOCATE(ZFLD(INTLO,INTLA)) 
ALLOCATE(ZFLDTMP(KLON,KLAT))
ALLOCATE(ZITM(INTLO,INTLA))
ALLOCATE(ILIND(KLEN))

ZSEUIL=10000._JPRB
ZLIMASK=0.5_JPRB
ILIMGP=1000000

!-------------------------------------------------------------------------------
!     2. INPUT DATA UPDATING.
!        -------------------

!     2.1  Input SST file reading

IF (LLNETCDF) THEN
  ! Open NETCDF file.
  ISTATUS=NF90_OPEN('OSTIA.nc',NF90_NOWRITE,INCIDHN)
  IF (ISTATUS /= NF90_NOERR) CALL HANDLE_ERR(ISTATUS)
  CALL LECNETCDF(INCIDHN,'analysed_sst',ZSST)

  !OSTIA SST is in hundredth of degrees celsius. Convert to K.
  ZSST(:,:,:)=ZSST(:,:,:)*0.01_JPRB+RTT 

  ILON=SIZE(ZSST,1)
  ILAT=SIZE(ZSST,2)

  ! FORCE JLON=1 OVER GREENWICH MERIDIAN, RATHER THAN DATE CHANGE MERIDIAN.
  DO J2=1,ILAT
    DO J1=1,ILON
      ZFLDTMP(J1,J2)=ZSST(MODULO(J1+ILON/2-1,ILON)+1,ILAT-J2+1,1)
    ENDDO
  ENDDO


  ! COPY INTO A LARGER ARRAY.
  DO JJ=IL1,IL2,INS
    DO JI=1+JPBX,ILON+JPBX
      ZFLD(JI,JJ)=ZFLDTMP(JI-JPBX,JJ-IL1+1)
    ENDDO
  ENDDO
  ! 2.2  Input land/sea mask
  ZITM(:,:)=1.0_JPRB
  WHERE ( ZFLD < 0.0_JPRB ) ZITM=0.0_JPRB !over land, ZFLD <0.

ELSE
  OPEN(UNIT=3,FILE=CLNOMF,FORM='UNFORMATTED')
  READ(3) ICO,IPNIV,ITOGOS
  READ(3) ((ZFLD(JI,JJ),JI=1+JPBX,KLON+JPBX),JJ=IL1,IL2,INS)
  CLOSE(3)

  ! 2.2  Input land/sea mask
  ZITM(:,:)=1.0_JPRB
  WHERE ( ZFLD > ZSEUIL ) ZITM=0.0_JPRB

ENDIF

IF ( ANY (ZITM < ZLIMASK) ) THEN
  WRITE(NULOUT,*) '    masked input SST '
ELSE
  WRITE(NULOUT,*) '    no mask in the input SST '
ENDIF
WRITE(NULOUT,*) '  '

!-------------------------------------------------------------------------------
!     3. DEFINITION OF THE SUBGRID (ARPEGE CASE MAINLY).
!        -------------------------

!     ZBO = Sine of the boundary latitudes of the NDGLG boxes.
!     ZSLA = Sine of the latitudes of the subgrid.
!     ZSLO = Sine of the longitudes of the subgrid.
!     ZCLO = Cosine of the longitudes of the subgrid.

IF (.NOT.LELAM) THEN
  ALLOCATE(ZSLA(INJQ))
  ALLOCATE(ZSLO(INIQ,NDGLG))
  ALLOCATE(ZCLO(INIQ,NDGLG))
  ITOGOS=0
  CALL GEO923(YDGEOMETRY,INPNPP,ITOGOS,ZBO,ZLONG,ZMU,ZSLA,ZSLO,ZCLO)
ELSE
  ITOGOS=INIQ*INJQ
  IF (ITOGOS /= KLEN) THEN
    WRITE(NULOUT,*) 'PB IN ALD ARRAYS : ITOGOS=',ITOGOS,'  KLEN=',KLEN,'  field ',KF
    CALL ABOR1 ('PB INCLITC')
  ENDIF
  CALL EGEO923(YDGEOMETRY%YRGEM)
ENDIF

!-------------------------------------------------------------------------------
!     4. INPUT FIELD INTERPOLATION.
!        -------------------------

IPNIV=1
ICO=INPNPP*INPNPP
IF (LELAM) THEN
  CALL EINTER6 (YDGEOMETRY%YREGEO,ZFLD,INTLA,INTLO,IPNIV,PS,INJQ,INIQ,KLEN,EDELX,EDELY,INPNPP,&
 & ICO,NULOUT,ZITM,PLSM,ZLIMASK)
ELSE
  CALL INTER6 (NDGLG,NDLON,INTLA,INTLO,IPNIV,INPNPP,ITOGOS,ICO,INJQ,INIQ,&
 & NLOENG(1:NDGLG),KLEN,PS,ZFLD,ZSLA,ZSLO,ZCLO,ZITM,PLSM,ZLIMASK)
  DEALLOCATE(ZSLA)
  DEALLOCATE(ZSLO)
  DEALLOCATE(ZCLO)
ENDIF



!-------------------------------------------------------------------------------
!     5. OUTPUT FILE CHECKING.
!        --------------------

!     8.1  First field checking

INS=0

ZSEUIL=350._JPRB

DO JI = 1,ITOGOS
  IF (PLSM(JI) > ZLIMASK .AND. PS(JI) > ZSEUIL) THEN
    INS=INS+1
    ILIND(INS)=JI
  ENDIF
ENDDO
WRITE(NULOUT,*) ' resultat interpolation :',INS,'points indefinis sur mer '

!     8.2  Second loop with a bigger INPNPP

IF (((.NOT.LELAM .AND. KLON<2*NDLON).OR.(LELAM .AND. KLEN<ILIMGP)) .AND. INS /= 0) THEN
  INPNPP=INPNPP*4+1
  ICO=INPNPP*INPNPP
  IPNIV=1
  ALLOCATE(ZS2(KLEN))
  IF (LELAM) THEN
    CALL EINTER6 (YDGEOMETRY%YREGEO,ZFLD,INTLA,INTLO,IPNIV,ZS2,INJQ,INIQ,KLEN,EDELX,EDELY,INPNPP,&
   & ICO,NULOUT,ZITM,PLSM,ZLIMASK)
  ELSE
    INJQ=NDGLG*INPNPP
    INIQ=NDLON*INPNPP
    ALLOCATE(ZSLA(INJQ))
    ALLOCATE(ZSLO(INIQ,NDGLG))
    ALLOCATE(ZCLO(INIQ,NDGLG))
    ITOGOS=0
    CALL GEO923(YDGEOMETRY,INPNPP,ITOGOS,ZBO,ZLONG,ZMU,ZSLA,ZSLO,ZCLO)
    CALL INTER6 (NDGLG,NDLON,INTLA,INTLO,IPNIV,INPNPP,ITOGOS,ICO,INJQ,INIQ,&
   & NLOENG(1:NDGLG),KLEN,ZS2,ZFLD,ZSLA,ZSLO,ZCLO,ZITM,PLSM,ZLIMASK)
  ENDIF
  DO JI = 1,INS
    PS(ILIND(JI))=ZS2(ILIND(JI))
  ENDDO
  ICO=0
  DO JI = 1,ITOGOS
    IF (PLSM(JI) > ZLIMASK .AND. PS(JI) > ZSEUIL) ICO=ICO+1
  ENDDO
  WRITE(NULOUT,*) '         ',ICO,'points indefinis sur mer apres 2e passage '
ENDIF

!     8.3  Deallocate temporary space

IF (ALLOCATED(ZFLD))    DEALLOCATE(ZFLD)
IF (ALLOCATED(ZFLDTMP)) DEALLOCATE(ZFLDTMP)
IF (ALLOCATED(ZSLA))    DEALLOCATE(ZSLA)
IF (ALLOCATED(ZSLO))    DEALLOCATE(ZSLO)
IF (ALLOCATED(ZCLO))    DEALLOCATE(ZCLO)
IF (ALLOCATED(ZS2))     DEALLOCATE(ZS2)

!-------------------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INCLITC',1,ZHOOK_HANDLE)
END SUBROUTINE INCLITC
