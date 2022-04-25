SUBROUTINE EINCLI7(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF)

!**** *EINCLI7*

!     PURPOSE
!     -------
!     This routine modifies climatological fields describing water surfaces
!     on a limited rectangular domain of the globe, in 1 climatological file.
!     Fields computed in EINCLI2 and EINCLI6 and not constant for water may 
!     be changed, if the corresponding high resolution data are available.
!     This routine may be used to introduce a sea/lake contrast.

!**   INTERFACE
!     ---------

!     CALL EINCLI7

!     METHOD
!     ------
!     This routine modifies up to 7 fields wich characterize seas and lakes :
!      - land use type
!      - useful soil depth     
!      - albedo and emissivity (according to climatological temperature)
!      - maximum soil depth
!      - climatological temperatures (superficial, mean)
!     The target grid is a regular grid (geographical or in plane projection).
!     The source grid can be any regular rectangular latitude by longitude
!     grid, the first point being at the NW edge, longitudes going eastwards,
!     and latitude going southwards (the NE edge is before the SW edge).
!     The earth is supposed to be flat, i.e. the periodicity of the longitudes
!     as well as the symmetry of the latitudes is ignored; this program can be
!     used for a global source grid (although it has not been designed for
!     this task), but the interpolation at the boundaries will not be optimal.
!     Inside the input domain, the variables are averaged in a rectangular
!     latitude x longitude box around each point of the target grid; the size
!     of the box is approximately twice the distance between two points in the
!     target grid.
!     The missing data have a negative value, and it is assumed that the
!     location of missing data is the same for each field. Outside the domain
!     and over parts where data are missing, the original value is kept. At
!     the boundaries of the domain, a linear combination of the old and new
!     values is performed. Values are not modified if the corresponding file
!     is not available.
!     Input data are units 10 to 15 :
!     -mask describing missing data              (unit 10,file msk_HR)
!       -9999. -> missing  1. -> available
!       values over land are expected to be missing
!     -dominant land use type                    (unit 11,file itp_HR)
!       only relevant values are expected : NTPMER, NTPLAC or -9999.
!     -g*height (useful to modify temperatures)  (unit 12,file rel_HR)
!     -superficial climatological temperature    (unit 13,file tsl_HR)
!     -mean        climatological temperature    (unit 14,file tpl_HR)
!     -maximum soil depth                        (unit 15,file dps_HR)
!     The output fields are added to or modified in the ARPEGE file.

!     EXTERNALS
!     ---------
!     ABOR1
!     CCHIEN
!     ELECI
!     EBICLI
!     ECHK923
!     FA-LFI package (FAITOU,LFILAF,FACILE,FAIRME)

!     AUTHORS
!     -------
!     S. Kertesz 16/12/1999 from EINCLI5

!     MODIFICATIONS
!     -------------
!     R. El Khatib : 01-12-06 Cleaning sm variables
!     M.Hamrud    03-10-01  CY28 Cleaning
!     D. Giard    04-09-15  cleaning
!     D. Giard    05-04-07 : no call to EBICLI
!     O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE YOEPHY       , ONLY : TEPHY
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCLI   , ONLY : YRCLI
USE YOMCST   , ONLY : RPI      ,RG
USE YOMLUN   , ONLY : NULOUT
USE YOMSTA   , ONLY : RDTDZ1
USE YOM_YGFL , ONLY : TYPE_GFLD

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT)   :: YDGEOMETRY
TYPE(TYPE_GFLD)    ,INTENT(INOUT):: YDGFL
TYPE(TEPHY)    ,INTENT(INOUT)   :: YDEPHY
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
INTEGER(KIND=JPIM) :: JPTYVE
REAL(KIND=JPRB) :: PPONDER

PARAMETER (PPONDER=.25_JPRB,JPTYVE=5)

! Input grid arrays:
REAL(KIND=JPRB) :: ZDPS0(YRCLI%NDATX,YRCLI%NDATY),ZITP0(YRCLI%NDATX,YRCLI%NDATY),ZMSK0(YRCLI%NDATX,YRCLI%NDATY)&
 & ,ZREL0(YRCLI%NDATX,YRCLI%NDATY),ZTPL0(YRCLI%NDATX,YRCLI%NDATY),ZTSL0(YRCLI%NDATX,YRCLI%NDATY)   
REAL(KIND=JPRD) ,ALLOCATABLE :: ZAUX(:)
! Output grid arrays (C+I):
REAL(KIND=JPRB) :: ZALB1((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1))&
 & ,ZDEP1((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1))&
 & ,ZDPS1((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1))&
 & ,ZEMI1((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1))&
 & ,ZLSM1((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1))&
 & ,ZITP1((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1))&
 & ,ZLND1((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1))&
 & ,ZREL1((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1))&
 & ,ZTPL1((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1))&
 & ,ZTSL1((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1))  
! Writing arrays (C+I+E)
REAL(KIND=JPRB) :: ZEXT(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON,0:7)
REAL(KIND=JPRB) :: ZITP(JPTYVE)
REAL(KIND=JPRB) :: ZCV, ZD, ZDLA, ZDLO, ZDPS, ZEPS, ZGST, ZITPM, ZITPN,&
 & ZLAT, ZLIS, ZLON, ZMIS, ZMSK, ZPTL, ZTPL, ZTSL, ZW1, ZW2, ZZDL
REAL(KIND=JPRB) :: ZHOOK_HANDLE

CHARACTER :: CLFORM*12,CLNOMC*16,CLNOMF*10
CHARACTER :: CLPREF(0:7)*8,CLSUFF(0:7)*12

INTEGER(KIND=JPIM) :: INIVL(0:7)
INTEGER(KIND=JPIM) :: IADL(YDGEOMETRY%YRDIM%NDGLG)
INTEGER(KIND=JPIM) :: I, I1LA, I1LO, I2LA, I2LO, IARI, IARP, IC, IDLA,&
 & IDLO, IDPS, IFLD, IINF, IITP, ILAT, ILON, IM,&
 & IMES, INDEX1, INDEX2, INIV, INUM, IOS, IPTL, IREP,&
 & IT, ITFING, ITPL, ITSL, IXFING, IYFING, IZ,&
 & J, JC, JJ, JLA, JLO, JT, JX, JY  

LOGICAL :: LLBIP(0:7),LLWRI(0:7),LLPAC(0:7)
LOGICAL :: LLCOSP, LLDPS, LLIMST, LLITP, LLMSK, LLOPEN,&
 & LLREL, LLTPL, LLTSL  

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "cchien.intfb.h"
#include "ebicli.intfb.h"
#include "echk923.intfb.h"
#include "eleci.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('EINCLI7',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDPHY1=>YDML_PHY_MF%YRPHY1)
ASSOCIATE(TMERGL=>YDPHY1%TMERGL, &
 & NTSTAGP=>YDGEM%NTSTAGP, &
 & YRGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, &
 & NDGLG=>YDDIM%NDGLG, NDLUXG=>YDDIM%NDLUXG, NDGUNG=>YDDIM%NDGUNG, &
 & NDGUXG=>YDDIM%NDGUXG, NDLUNG=>YDDIM%NDLUNG, NDLON=>YDDIM%NDLON)
!     ------------------------------------------------------------------

!     1. SET INITIAL VALUES.
!        -------------------

!     1.1 Miscellaneous

! Final grid geometry:

IXFING=NDLUXG-NDLUNG+1
IYFING=NDGUXG-NDGUNG+1
ITFING=IXFING*IYFING

ZCV= 180._JPRB/RPI
ZEPS= 1.E-10_JPRB
ZGST= RDTDZ1/RG
ZLIS= 1.0_JPRB/REAL(YRCLI%NPINT,JPRB)
ZMIS= YRCLI%SMANQ - 1.0_JPRB

IADL(1:NDGLG)=NTSTAGP(1:NDGLG)-1

!     1.2 Type of data files

IF (YRCLI%LIEEE) THEN
  CLFORM='UNFORMATTED'
ELSE
  CLFORM='FORMATTED'
ENDIF

!     1.3 Missing data

!  Default is missing data
LLMSK=.FALSE.
LLDPS=.FALSE.
LLITP=.FALSE.
LLREL=.FALSE.
LLTPL=.FALSE.
LLTSL=.FALSE.
IDPS= 0
IITP= 0
ITPL= 0
ITSL= 0
ZMSK0(:,:)= ZMIS
ZDPS0(:,:)= ZMIS
ZITP0(:,:)= ZMIS
ZREL0(:,:)= ZMIS
ZTPL0(:,:)= ZMIS
ZTSL0(:,:)= ZMIS

!     1.4 Final biperiodization and writing

DO J=0,7
  LLBIP(J)=.TRUE.
  LLWRI(J)=.TRUE.
  LLPAC(J)=.TRUE.
  INIVL(J)=0
ENDDO
LLBIP(0)=.FALSE.
LLBIP(1)=.FALSE.
LLWRI(0)=.FALSE.
INIVL(6)=1
INIVL(7)=1

DO JC=0,7
  DO J=1,NDGLG*NDLON
    ZEXT(J,JC)=0.0_JPRB
  ENDDO
ENDDO

DO J=0,6
  CLPREF(J)='SURF'
ENDDO
CLPREF(7)='PROF'

CLSUFF( 0)='IND.TERREMER'
CLSUFF( 1)='IND.VEG.DOMI'
CLSUFF( 2)='EPAIS.SOL   '
CLSUFF( 3)='EPAI.SOL.MAX'
CLSUFF( 4)='ALBEDO      '
CLSUFF( 5)='EMISSIVITE  '
CLSUFF( 6)='TEMPERATURE '
CLSUFF( 7)='TEMPERATURE '

!     ------------------------------------------------------------------

!     2. CHECK AND READ DATASETS.
!        ------------------------

!     2.1 Check which files are available

!  Common mask
LLOPEN=.FALSE.
IOS= 0
INQUIRE(FILE='msk_HR',IOSTAT=IOS,EXIST=LLMSK,OPENED=LLOPEN)
LLMSK= LLMSK .AND. (IOS == 0) .AND. .NOT.LLOPEN
IF (.NOT.LLMSK) THEN
  WRITE(NULOUT,'('' ERROR IN EINCLI7 :'',&
   & '' THE COMMON MASK (DEFINING MISSING DATA) IS NOT GIVEN !'')')  
  CALL ABOR1('EINCLI7 : THE COMMON MASK (DEFINING MISSING DATA) IS NOT GIVEN !')
ENDIF

!  Land use type
INQUIRE(FILE='itp_HR',IOSTAT=IOS,EXIST=LLITP,OPENED=LLOPEN)
LLITP= LLITP .AND. (IOS == 0) .AND. .NOT.LLOPEN
! Orography
INQUIRE(FILE='rel_HR',IOSTAT=IOS,EXIST=LLREL,OPENED=LLOPEN)
LLREL= LLREL .AND. (IOS == 0) .AND. .NOT.LLOPEN
! Superficial temperature
INQUIRE(FILE='tsl_HR',IOSTAT=IOS,EXIST=LLTSL,OPENED=LLOPEN)
LLTSL= LLTSL .AND. (IOS == 0) .AND. .NOT.LLOPEN
! Mean temperature
INQUIRE(FILE='tpl_HR',IOSTAT=IOS,EXIST=LLTPL,OPENED=LLOPEN)
LLTPL= LLTPL .AND. (IOS == 0) .AND. .NOT.LLOPEN
!  Maximum soil depth
INQUIRE(FILE='dps_HR',IOSTAT=IOS,EXIST=LLDPS,OPENED=LLOPEN)
LLDPS= LLDPS .AND. (IOS == 0) .AND. .NOT.LLOPEN

!     2.2 Check which fields will be modified

IF (.NOT.(LLDPS.OR.LLITP.OR.LLTSL.OR.LLTPL))&
 & CALL ABOR1('ERROR IN EINCLI7 : NO FIELDS TO MODIFY !')  

IF ((LLTSL.OR.LLTPL) .AND. .NOT.LLREL)&
 & CALL ABOR1(' ERROR IN EINCLI7 : NO HEIGHT CORR. FOR TEMPERATURE !')  

LLREL= LLREL.AND.(LLTSL.OR.LLTPL)

!     2.3 Read and check the mask (missing)

IF (YRCLI%LIEEE) THEN
  ALLOCATE (ZAUX(YRCLI%NDATX*YRCLI%NDATY))
ENDIF

OPEN(UNIT=10,FILE='msk_HR',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  READ(10) ZAUX
  DO JJ=YRCLI%NDATY,1,-1
    DO J=1,YRCLI%NDATX
      IZ=J+(YRCLI%NDATY-JJ)*YRCLI%NDATX
      ZMSK0(J,JJ)=ZAUX(IZ)
    ENDDO
  ENDDO
ELSE
  READ(10,*) ((ZMSK0(J,JJ),J=1,YRCLI%NDATX),JJ=YRCLI%NDATY,1,-1)
ENDIF
CLOSE(10)
DO JJ=1,YRCLI%NDATY
  DO J=1,YRCLI%NDATX
    IF ((ZMSK0(J,JJ) > YRCLI%SMANQ) .AND. (ABS(ZMSK0(J,JJ)-1.0_JPRB) > ZEPS))&
     & CALL ABOR1(' ERROR IN EINCLI7 : UNEXPECTED MASK VALUE !')  
  ENDDO
ENDDO

!     2.4 Read new data

!  Land use type
IF (LLITP) THEN
  OPEN(UNIT=11,FILE='itp_HR',FORM=CLFORM)
  IF (YRCLI%LIEEE) THEN
    READ(11) ZAUX
    DO JJ=YRCLI%NDATY,1,-1
      DO J=1,YRCLI%NDATX
        IZ=J+(YRCLI%NDATY-JJ)*YRCLI%NDATX
        ZITP0(J,JJ)=ZAUX(IZ)
      ENDDO
    ENDDO
  ELSE
    READ(11,*) ((ZITP0(J,JJ),J=1,YRCLI%NDATX),JJ=YRCLI%NDATY,1,-1)
  ENDIF
  CLOSE(11)
  IITP= 1
ENDIF

!  g*Orography
IF (LLREL) THEN
  OPEN(UNIT=12,FILE='rel_HR',FORM=CLFORM)
  IF (YRCLI%LIEEE) THEN
    READ(12) ZAUX
    DO JJ=YRCLI%NDATY,1,-1
      DO J=1,YRCLI%NDATX
        IZ=J+(YRCLI%NDATY-JJ)*YRCLI%NDATX
        ZREL0(J,JJ)=ZAUX(IZ)
      ENDDO
    ENDDO
  ELSE
    READ(12,*) ((ZREL0(J,JJ),J=1,YRCLI%NDATX),JJ=YRCLI%NDATY,1,-1)
  ENDIF
  CLOSE(12)
ENDIF

!  Superficial temperature
IF (LLTSL) THEN
  OPEN(UNIT=13,FILE='tsl_HR',FORM=CLFORM)
  IF (YRCLI%LIEEE) THEN
    READ(13) ZAUX
    DO JJ=YRCLI%NDATY,1,-1
      DO J=1,YRCLI%NDATX
        IZ=J+(YRCLI%NDATY-JJ)*YRCLI%NDATX
        ZTSL0(J,JJ)=ZAUX(IZ)
      ENDDO
    ENDDO
  ELSE
    READ(13,*) ((ZTSL0(J,JJ),J=1,YRCLI%NDATX),JJ=YRCLI%NDATY,1,-1)
  ENDIF
  CLOSE(13)
  ITSL= 1
  DO JJ=1,YRCLI%NDATY
    DO J=1,YRCLI%NDATX
      ZTSL0(J,JJ)= ZTSL0(J,JJ) - ZGST*ZREL0(J,JJ)
    ENDDO
  ENDDO
ENDIF

!  Mean temperature
IF (LLTPL) THEN
  OPEN(UNIT=14,FILE='tpl_HR',FORM=CLFORM)
  IF (YRCLI%LIEEE) THEN
    READ(14) ZAUX
    DO JJ=YRCLI%NDATY,1,-1
      DO J=1,YRCLI%NDATX
        IZ=J+(YRCLI%NDATY-JJ)*YRCLI%NDATX
        ZTPL0(J,JJ)=ZAUX(IZ)
      ENDDO
    ENDDO
  ELSE
    READ(14,*) ((ZTPL0(J,JJ),J=1,YRCLI%NDATX),JJ=YRCLI%NDATY,1,-1)
  ENDIF
  CLOSE(14)
  ITPL= 1
  DO JJ=1,YRCLI%NDATY
    DO J=1,YRCLI%NDATX
      ZTPL0(J,JJ)= ZTPL0(J,JJ) - ZGST*ZREL0(J,JJ)
    ENDDO
  ENDDO
ENDIF

!  Maximum soil depth
IF (LLDPS) THEN
  OPEN(UNIT=15,FILE='dps_HR',FORM=CLFORM)
  IF (YRCLI%LIEEE) THEN
    READ(15) ZAUX
    DO JJ=YRCLI%NDATY,1,-1
      DO J=1,YRCLI%NDATX
        IZ=J+(YRCLI%NDATY-JJ)*YRCLI%NDATX
        ZDPS0(J,JJ)=ZAUX(IZ)
      ENDDO
    ENDDO
  ELSE
    READ(15,*) ((ZDPS0(J,JJ),J=1,YRCLI%NDATX),JJ=YRCLI%NDATY,1,-1)
  ENDIF
  CLOSE(15)
  IDPS= 1
ENDIF

!  Deallocate temporary space
IF (YRCLI%LIEEE) THEN
  DEALLOCATE (ZAUX)
ENDIF

IC= IDPS + IITP + ITPL + ITSL

!     ------------------------------------------------------------------

!     3. READ ARPEGE FIELDS.
!        -------------------

INUM=3
IREP=0
IARI=0
IARP=31
IMES=1
CLNOMF='Const.Clim'
CLNOMC='Const.Clim.Surfa'
LLIMST=.TRUE.

IINF=-1
LLCOSP=.FALSE.

CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'OLD',.TRUE.,LLIMST,&
 & IMES,IARP,IARI,CLNOMC)  
CALL CCHIEN(YDGEOMETRY,CLNOMC,INUM,IINF)

INIV=0
CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'SURF',INIV,'IND.TERREMER',ZLSM1,NULOUT)
CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'SURF',INIV,'GEOPOTENTIEL',ZREL1,NULOUT)
CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'SURF',INIV,'PROP.TERRE  ',ZLND1,NULOUT)
CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'SURF',INIV,'IND.VEG.DOMI',ZITP1,NULOUT)
CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'SURF',INIV,'EPAIS.SOL   ',ZDEP1,NULOUT)
CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'SURF',INIV,'EPAI.SOL.MAX',ZDPS1,NULOUT)
CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'SURF',INIV,'ALBEDO      ',ZALB1,NULOUT)
CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'SURF',INIV,'EMISSIVITE  ',ZEMI1,NULOUT)

INIV=1
CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'SURF',INIV,'TEMPERATURE ',ZTSL1,NULOUT)
CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'PROF',INIV,'TEMPERATURE ',ZTPL1,NULOUT)

ZEXT(1:ITFING,0)=ZLSM1(1:ITFING)

!     ------------------------------------------------------------------

!     4. LOOP ON THE ALADIN GRID.
!       ------------------------

! ALADIN GEOMETRY: zone C+I
!   NDGUNG to NDGUXG, NDLUNG to NDLUXG
!   I = index of the element in the final grid

I=0
IM=0
IT=0
DO JJ=1,IYFING
  DO J=1,IXFING
    I=I+1

!  Only water points are considered here :
    IF (ZLSM1(I) > YRCLI%SMASK) GO TO 410
    IT=IT+1

!  Check the position of the point relative to the dataset and define the
!  interpolation box if required. Only points inside the domain or close to
!  its boundary are considered here.
!   ILAT/ILON : latitude/longitude of the nearest point in the input grid
!   IDLA/IDLO : latitude/longitude radius of the box
!   I1LA/I2LA/I1LO/I2LO : north/south/west/east edges of the box
!  The size of the box inside which input data are considered and averaged 
!  is by default twice the model mesh-size if each direction. This can be 
!  reduced (to keep sharp gradients) by increasing the namelist variable 
!  NPINT (1 by default), i.e. decreasing ZLIS.

!  Check the latitude
    ZLAT= YRGSGEOM_NB%GELAT(J+IADL(JJ))*ZCV
    ILAT= NINT(0.5_JPRB+(ZLAT-YRCLI%ELATSW)/YRCLI%EDLAT)
!  Latitude mesh size evaluation (including projected grids) (in degrees):
    IF (JJ < IYFING.AND.JJ > 1) THEN
      ZZDL= (YRGSGEOM_NB%GELAT(J+IADL(JJ+1))-YRGSGEOM_NB%GELAT(J+IADL(JJ-1)))*ZCV/2.0_JPRB
    ELSEIF(JJ < IYFING) THEN
      ZZDL= (YRGSGEOM_NB%GELAT(J+IADL(JJ+1))-YRGSGEOM_NB%GELAT(J+IADL(JJ  )))*ZCV
    ELSE
      ZZDL= (YRGSGEOM_NB%GELAT(J+IADL(JJ  ))-YRGSGEOM_NB%GELAT(J+IADL(JJ-1)))*ZCV
    ENDIF
    ZDLA= ZLIS*ZZDL/YRCLI%EDLAT
    IDLA= INT(ZDLA)
    IF (ILAT < (1-IDLA) .OR. ILAT > (YRCLI%NDATY+IDLA)) GO TO 410
    I1LA= MAX(1    ,ILAT-IDLA)
    I2LA= MIN(YRCLI%NDATY,ILAT+IDLA)
!  Check the longitude
    ZLON= YRGSGEOM_NB%GELAM(J+IADL(JJ))*ZCV
    ZZDL=MOD(ZLON-YRCLI%ELONSW,360._JPRB)
    IF (ZZDL < 0) ZZDL=ZZDL+360._JPRB
    ILON= NINT(0.5_JPRB+(ZZDL)/YRCLI%EDLON)

    IF (ILON <= 0.0_JPRB.OR.ILON > YRCLI%NDATX) GOTO 410
! Longitude mesh size evaluation (including projected grids) (in degrees):
    IF (J < IXFING.AND.J > 1) THEN
      ZZDL= (YRGSGEOM_NB%GELAM(J+1+IADL(JJ))-YRGSGEOM_NB%GELAM(J-1+IADL(JJ)))*ZCV/2.0_JPRB
      IF (ZZDL < 0.0_JPRB) ZZDL=ZZDL+180._JPRB
    ELSE
      IF(J < IXFING) THEN
        ZZDL= (YRGSGEOM_NB%GELAM(J+1+IADL(JJ))-YRGSGEOM_NB%GELAM(  J+IADL(JJ)))*ZCV
      ELSE
        ZZDL= (YRGSGEOM_NB%GELAM(  J+IADL(JJ))-YRGSGEOM_NB%GELAM(J-1+IADL(JJ)))*ZCV
      ENDIF
! ZZDL MUST be >0:
      ZZDL=MOD(ZZDL,360._JPRB)
      IF (ZZDL < 0.0_JPRB) ZZDL=ZZDL+360._JPRB
    ENDIF
    ZDLO= ZLIS*ZZDL/YRCLI%EDLON
    IDLO= MIN(INT(ZDLO),YRCLI%NDATX/4)
    I1LO=MAX(ILON-IDLO,1)
    I2LO=MIN(ILON+IDLO,YRCLI%NDATX)

!  Initialize new fields
    ZMSK= 0.0_JPRB
    ZDPS= 0.0_JPRB
    ZTSL= 0.0_JPRB
    ZTPL= 0.0_JPRB
    DO JT=1,JPTYVE
      ZITP(JT)= 0.0_JPRB
    ENDDO

!  Compute new values
    IPTL=0
    DO JLA=I1LA,I2LA
      DO JLO=I1LO,I2LO
!  Mask
        ZMSK= ZMSK + MAX(0.0_JPRB,ZMSK0(JLO,JLA))
!  Land use type
        DO JT=1,JPTYVE
          IF (NINT(ZITP0(JLO,JLA)) == JT) ZITP(JT)= ZITP(JT) + 1.0_JPRB
        ENDDO
!  Soil depth        
        ZDPS= ZDPS + MAX(0.0_JPRB,ZDPS0(JLO,JLA))        
!  Superficial temperature 
        ZTSL= ZTSL + MAX(0.0_JPRB,ZTSL0(JLO,JLA))
!  Mean temperature               
        ZTPL= ZTPL + MAX(0.0_JPRB,ZTPL0(JLO,JLA))
! Number of points in the box:
        IPTL=IPTL+1
      ENDDO
    ENDDO
!  End of the two loops in the rectangular latitude x longitude box
!  ZMSK : Number of non-missing data in the box
!  ZPTL : Maximum number of water points in the box
!         (including the part of the box outside the domain)
    IF (ZMSK < ZEPS) GO TO 410
    IM=IM+1
!c This cannot work with the borders' cut
!c      ZPTL=REAL((2*IDLA+1)*(2*IDLO+1),KIND(ZPTL))*ZLND1(I)
    ZPTL=IPTL*(1.0_JPRB-ZLND1(I))
    ZPTL=MAX(ZMSK,ZPTL)

!  Compute new fields
!  Choice of the dominant type of vegetation. The dominant vegetation must
!  exceed 25%, otherwise the old type is kept.
    ZITPN= ZITP1(I)
    ZITPM= PPONDER*ZPTL
    DO JT=1,JPTYVE
      IF ((JT == YRCLI%NTPMER.OR.JT == YRCLI%NTPLAC) .AND. (ZITP(JT) >= ZITPM)) THEN      
        ZITPN= REAL(JT,JPRB)
        ZITPM= ZITP(JT)
      ENDIF
    ENDDO
    ZITP1(I)= ZITPN*IITP + ZITP1(I)*(1-IITP)
!  For the other parameters, a linear combination is done with the weights
!  ZW1 for the input data and ZW2 for the initial data, provided data are not
!  missing (Ixxx=1). Else (Ixxx=0), the previous value is kept.
    ZW1=1.0_JPRB/ZPTL
    ZW2=1.0_JPRB-ZMSK*ZW1
!  Superficial and mean temperature (with height correction)
    IF (LLTPL.AND..NOT.LLTSL) ZTSL=ZTPL
    IF (LLTSL.AND..NOT.LLTPL) ZTPL=ZTSL
    ZTSL= ZTSL + ZMSK*ZGST*ZREL1(I)
    ZTPL= ZTPL + ZMSK*ZGST*ZREL1(I)
    ZTSL1(I)= (ZW1*ZTSL + ZW2*ZTSL1(I))*ITSL + ZTSL1(I)*(1-ITSL)
    ZTPL1(I)= (ZW1*ZTPL + ZW2*ZTPL1(I))*ITPL + ZTPL1(I)*(1-ITPL)
!  Soil  depth
    ZDPS1(I)= (ZW1*ZDPS + ZW2*ZDPS1(I))*IDPS + ZDPS1(I)*(1-IDPS)
    ZDPS1(I)= MIN(YRCLI%SDEPX,MAX(YRCLI%SDEPN,ZDPS1(I)))
    ZDEP1(I)= ZDPS1(I)
!  Correction of albedo and emissivity
    IF (ZTSL1(I) <= TMERGL) THEN            
      ZALB1(I)= YRCLI%SALBB                        
      ZEMI1(I)= YRCLI%SEMIB
    ELSE              
      ZALB1(I)= YRCLI%SALBM
      ZEMI1(I)= YRCLI%SEMIM
    ENDIF   
        
410 CONTINUE
  ENDDO
ENDDO

IF (IM == 0) CALL ABOR1(' ERROR IN EINCLI7 : NO MODIFICATION !')
WRITE(NULOUT,'(I6,''/'',I6,'' POINTS MODIFIED BY EINCLI7'')')IM,IT
WRITE(NULOUT,'(I2,'' FIELDS MODIFIED BY EINCLI7'')') IC

!     ------------------------------------------------------------------

!     5. BIPER + WRITING THE CLIM FILE.
!        ------------------------------

! Special case for land-use type: no biperiodization
! Find out the Globally dominant type:
DO JT=1,JPTYVE
  ZITP(JT)=0.0_JPRB
ENDDO
DO JJ=1,ITFING
  IT=NINT(ZITP1(JJ))
  ZITP(IT)=ZITP(IT)+1
ENDDO
IT=0
ZD=0.0_JPRB
DO JT=1,JPTYVE
  IF (ZD < ZITP(JT))THEN
    IT=JT
    ZD=ZITP(JT)
  ENDIF
ENDDO
! Now use this type IT to fill the Extension zone
INDEX1=0
INDEX2=0
DO JY=1,IYFING
  DO JX=1,IXFING
    ZEXT(JX+INDEX2,1)=ZITP1(JX+INDEX1)
  ENDDO
  DO JX=IXFING+1,NDLON
    ZEXT(JX+INDEX2,1)=IT
  ENDDO
  INDEX1=INDEX1+IXFING
  INDEX2=INDEX2+NDLON
ENDDO
DO JY=IYFING+1,NDGLG
  DO JX=1,NDLON
    ZEXT(JX+INDEX2,1)=IT
  ENDDO
  INDEX2=INDEX2+NDLON
ENDDO

! Biperiodization and storing for other fields:

DO J=1,ITFING 
  ZEXT(J, 2)=ZDEP1(J)
  ZEXT(J, 3)=ZDPS1(J)  
  ZEXT(J, 4)=ZALB1(J)
  ZEXT(J, 5)=ZEMI1(J)
  ZEXT(J, 6)=ZTSL1(J)
  ZEXT(J, 7)=ZTPL1(J)   
ENDDO
IFLD=7+1
CALL EBICLI(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF,IFLD,INIVL,CLPREF,CLSUFF,INUM,ZEXT,LLBIP,LLWRI,LLPAC)

CALL ECHK923(YDGEOMETRY%YRDIM,INUM)

CALL LFILAF(IREP,INUM,.FALSE.)
CALL FAIRME(IREP,INUM,'KEEP')

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EINCLI7',1,ZHOOK_HANDLE)
END SUBROUTINE EINCLI7
