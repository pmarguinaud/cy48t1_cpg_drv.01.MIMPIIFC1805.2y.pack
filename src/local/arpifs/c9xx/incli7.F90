SUBROUTINE INCLI7(YDGEOMETRY,YDPHY1)

!**** *INCLI7*

!     PURPOSE
!     -------
!     This routine modifies climatological fields describing water surfaces
!     on a limited rectangular domain of the globe, in 1 climatological file.
!     Fields computed in INCLI2 and INCLI6 and not constant for water may 
!     be changed, if the corresponding high resolution data are available.
!     This routine may be used to introduce a sea/lake contrast.

!**   INTERFACE
!     ---------

!     CALL INCLI7

!     METHOD
!     ------
!     This routine modifies up to 7 fields wich characterize seas and lakes :
!      - land use type
!      - useful soil depth     
!      - albedo and emissivity (according to climatological temperature)
!      - maximum soil depth
!      - climatological temperatures (superficial, mean)
!     The target grid is a Gaussian grid with rotation of the pole, "Schmidt"
!     stretching as well as reduction of the number of grid points towards the
!     pole of the final representation.
!     The source grid can be any regular rectangular latitude by longitude 
!     grid, the first point being at the NW edge, longitudes going eastwards,
!     and latitude going northwards (the NE edge is before the SW edge).
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

!     FA-LFI package (FAITOU,FACILE,FAIENC,LFILAF,FAIRME)
!     CHIEN
!     ABOR1
!     CHK923

!     AUTHORS
!     -------
!      D. Giard    00-01-31  from EINCLI7

!     MODIFICATIONS
!     -------------
!      M.Hamrud    03-10-01  CY28 Cleaning
!      D. Giard    04-09-15  cleaning
!      K. Yessad (Jan 2010): externalisation of group EGGX in XRD/IFSAUX
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM and TCSGLEG
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCT0   , ONLY : NQUAD
USE YOMLUN   , ONLY : NULOUT
USE YOMCST   , ONLY : RPI, RG
USE YOMVERT  , ONLY : VP00
USE YOMSTA   , ONLY : RDTDZ1
USE YOMPHY1  , ONLY : TPHY1
USE YOMCLI   , ONLY : YRCLI  

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN)   :: YDGEOMETRY
TYPE(TPHY1)    ,INTENT(INOUT):: YDPHY1
INTEGER(KIND=JPIM) :: JPTYVE
REAL(KIND=JPRB) :: PPONDER

PARAMETER (PPONDER=.25_JPRB,JPTYVE=5)

REAL(KIND=JPRB) :: ZDPS0(YRCLI%NDATX,YRCLI%NDATY),ZITP0(YRCLI%NDATX,YRCLI%NDATY),ZMSK0(YRCLI%NDATX,YRCLI%NDATY)&
 & ,ZREL0(YRCLI%NDATX,YRCLI%NDATY),ZTPL0(YRCLI%NDATX,YRCLI%NDATY),ZTSL0(YRCLI%NDATX,YRCLI%NDATY)  
REAL(KIND=JPRB) :: ZALB1(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON),ZDEP1(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON),&
 & ZDPS1(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON)&
 & ,ZEMI1(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON),ZLSM1(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON),&
 & ZITP1(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON)&
 & ,ZLND1(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON),ZREL1(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON),&
 & ZTPL1(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON)&
 & ,ZTSL1(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON)  
REAL(KIND=JPRB) :: ZITP(JPTYVE)
REAL(KIND=JPRB) ,ALLOCATABLE :: ZAUX(:)
REAL(KIND=JPRB) :: ZCV, ZDL, ZDLA, ZDLO, ZDPS, ZEPS, ZITPM, ZITPN, ZLAT, ZLON,&
 & ZMSK, ZPTL, ZTPL, ZTSL, ZW1, ZW2, ZZDL, ZLIS ,ZGST ,ZMIS  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

CHARACTER :: CLNOMC*16,CLNOMF*10,CLFORM*12

INTEGER(KIND=JPIM) :: I, I1LA, I1LO, I2LA, I2LO, IA, IARI, IARP, IC, ID, IDLA,&
 & IDLO, IDPS, IINF, IITP, ILAT, ILON, IM, IMES, INIV, INUM,&
 & IOS, IREP, IT, ITPL, ITSL, IZ, J, JJ, JLA, JLO, JT  

LOGICAL :: LLCOSP, LLDPS, LLIMST, LLITP, LLMSK, LLOPEN,&
 & LLPOLE, LLREL, LLTPL, LLTSL  

!     ------------------------------------------------------------------

#include "chien.h"

#include "abor1.intfb.h"
#include "chk923.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INCLI7',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDCSGLEG=>YDGEOMETRY%YRCSGLEG,  &
 & YDSTA=>YDGEOMETRY%YRSTA, YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE,  &
 & YDLAP=>YDGEOMETRY%YRLAP, YDVSPLIP=>YDGEOMETRY%YRVSPLIP,  &
 & YDVSLETA=>YDGEOMETRY%YRVSLETA, &
  & YDHSLMER=>YDGEOMETRY%YRHSLMER, YDCSGEOM=>YDGEOMETRY%YRCSGEOM, YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,  &
  & YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(TMERGL=>YDPHY1%TMERGL, &
 & NDGENG=>YDDIM%NDGENG, NDGLG=>YDDIM%NDGLG, NDGSAG=>YDDIM%NDGSAG, &
 & NDLON=>YDDIM%NDLON, NSMAX=>YDDIM%NSMAX, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NHTYP=>YDGEM%NHTYP, NLOENG=>YDGEM%NLOENG, NMENG=>YDGEM%NMENG, &
 & NSTAGP=>YDGEM%NSTAGP, NSTTYP=>YDGEM%NSTTYP, RLOCEN=>YDGEM%RLOCEN, &
 & RMUCEN=>YDGEM%RMUCEN, RSTRET=>YDGEM%RSTRET)
!     ------------------------------------------------------------------

!     1. SET INITIAL VALUES.
!        -------------------

!     1.1 Miscellaneous

ZCV= 180._JPRB/RPI
ZDL= 180._JPRB/REAL(NDGLG,JPRB)
ZEPS= 1.E-10_JPRB
ZGST= RDTDZ1/RG
ZLIS= 1.0_JPRB/REAL(YRCLI%NPINT,JPRB)
ZMIS= YRCLI%SMANQ - 1.0_JPRB

ID= 1

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
  WRITE(NULOUT,'('' ERROR IN INCLI7 :'',&
   & '' THE COMMON MASK (DEFINING MISSING DATA) IS NOT GIVEN !'')')  
  CALL ABOR1('INCLI5')
ENDIF

!  Land use type
INQUIRE(FILE='itp_HR',IOSTAT=IOS,EXIST=LLITP,OPENED=LLOPEN)
LLITP= LLITP .AND. (IOS == 0) .AND. .NOT.LLOPEN
! Orography
INQUIRE(FILE='rel_HR',IOSTAT=IOS,EXIST=LLTPL,OPENED=LLOPEN)
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
 & CALL ABOR1('ERROR IN INCLI7 : NO FIELDS TO MODIFY !')  

IF ((LLTSL.OR.LLTPL) .AND. .NOT.LLREL)&
 & CALL ABOR1(' ERROR IN INCLI7 : NO HEIGHT CORR. FOR TEMPERATURE !')  

LLREL= LLREL.AND.(LLTSL.OR.LLTPL)

!     2.3 Read and check the mask (missing)

IF (YRCLI%LIEEE) THEN
  ALLOCATE (ZAUX(YRCLI%NDATX*YRCLI%NDATY))
ENDIF

OPEN(UNIT=10,FILE='msk_HR',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  READ(10) ZAUX
  DO JJ=1,YRCLI%NDATY
    DO J=1,YRCLI%NDATX
      IZ=J+(JJ-1)*YRCLI%NDATX
      ZMSK0(J,JJ)=ZAUX(IZ)
    ENDDO
  ENDDO
ELSE
  READ(10,*) ((ZMSK0(J,JJ),J=1,YRCLI%NDATX),JJ=1,YRCLI%NDATY)
ENDIF
CLOSE(10)
DO JJ=1,YRCLI%NDATY
  DO J=1,YRCLI%NDATX
    IF ((ZMSK0(J,JJ) > YRCLI%SMANQ) .AND. (ABS(ZMSK0(J,JJ)-1.0_JPRB) > ZEPS))&
     & CALL ABOR1(' ERROR IN INCLI7 : UNEXPECTED MASK VALUE !')           
  ENDDO
ENDDO

!     2.4 Read new data

!  Land use type
IF (LLITP) THEN
  OPEN(UNIT=11,FILE='itp_HR',FORM=CLFORM)
  IF (YRCLI%LIEEE) THEN
    READ(11) ZAUX
    DO JJ=1,YRCLI%NDATY
      DO J=1,YRCLI%NDATX
        IZ=J+(JJ-1)*YRCLI%NDATX
        ZITP0(J,JJ)=ZAUX(IZ)
      ENDDO
    ENDDO
  ELSE
    READ(11,*) ((ZITP0(J,JJ),J=1,YRCLI%NDATX),JJ=1,YRCLI%NDATY)
  ENDIF
  CLOSE(11)
  IITP= 1
ENDIF

!  g*Orography
IF (LLREL) THEN
  OPEN(UNIT=12,FILE='rel_HR',FORM=CLFORM)
  IF (YRCLI%LIEEE) THEN
    READ(12) ZAUX
    DO JJ=1,YRCLI%NDATY
      DO J=1,YRCLI%NDATX
        IZ=J+(JJ-1)*YRCLI%NDATX
        ZREL0(J,JJ)=ZAUX(IZ)
      ENDDO
    ENDDO
  ELSE
    READ(12,*) ((ZREL0(J,JJ),J=1,YRCLI%NDATX),JJ=1,YRCLI%NDATY)
  ENDIF
  CLOSE(12)
ENDIF

!  Superficial temperature
IF (LLTSL) THEN
  OPEN(UNIT=13,FILE='tsl_HR',FORM=CLFORM)
  IF (YRCLI%LIEEE) THEN
    READ(13) ZAUX
    DO JJ=1,YRCLI%NDATY
      DO J=1,YRCLI%NDATX
        IZ=J+(JJ-1)*YRCLI%NDATX
        ZTSL0(J,JJ)=ZAUX(IZ)
      ENDDO
    ENDDO
  ELSE
    READ(13,*) ((ZTSL0(J,JJ),J=1,YRCLI%NDATX),JJ=1,YRCLI%NDATY)
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
    DO JJ=1,YRCLI%NDATY
      DO J=1,YRCLI%NDATX
        IZ=J+(JJ-1)*YRCLI%NDATX
        ZTPL0(J,JJ)=ZAUX(IZ)
      ENDDO
    ENDDO
  ELSE
    READ(14,*) ((ZTPL0(J,JJ),J=1,YRCLI%NDATX),JJ=1,YRCLI%NDATY)
  ENDIF
  CLOSE(14)
  ITSL= 1
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
    READ(15)ZAUX
    DO JJ=1,YRCLI%NDATY
      DO J=1,YRCLI%NDATX
        IZ=J+(JJ-1)*YRCLI%NDATX
        ZDPS0(J,JJ)=ZAUX(IZ)
      ENDDO
    ENDDO
  ELSE
    READ(15,*) ((ZDPS0(J,JJ),J=1,YRCLI%NDATX),JJ=1,YRCLI%NDATY)
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
LLPOLE=.TRUE.
LLCOSP=.FALSE.

CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'OLD',.TRUE.,LLIMST,&
 & IMES,IARP,IARI,CLNOMC)  
CALL CHIEN(CLNOMC,NSTTYP,RMUCEN,RLOCEN,RSTRET,NSMAX,&
 & NDGLG,NDLON,NLOENG,NMENG,NHTYP,NFLEVG,VP00,YDVAB%VALH,YDVAB%VBH,&
 & NQUAD,IINF,NDGSAG,NDGENG,ZEPS,LLPOLE,NULOUT)  
IF (LLPOLE) CALL ABOR1(' CLIM. FILES MUST NOT HAVE POLES !')

INIV=0
CALL FACILE(IREP,INUM,'SURF',INIV,'IND.TERREMER',ZLSM1,LLCOSP)
CALL FACILE(IREP,INUM,'SURF',INIV,'GEOPOTENTIEL',ZREL1,LLCOSP)
CALL FACILE(IREP,INUM,'SURF',INIV,'PROP.TERRE  ',ZLND1,LLCOSP)
CALL FACILE(IREP,INUM,'SURF',INIV,'IND.VEG.DOMI',ZITP1,LLCOSP)
CALL FACILE(IREP,INUM,'SURF',INIV,'EPAIS.SOL   ',ZDEP1,LLCOSP)
CALL FACILE(IREP,INUM,'SURF',INIV,'EPAI.SOL.MAX',ZDPS1,LLCOSP)
CALL FACILE(IREP,INUM,'SURF',INIV,'ALBEDO      ',ZALB1,LLCOSP)
CALL FACILE(IREP,INUM,'SURF',INIV,'EMISSIVITE  ',ZEMI1,LLCOSP)
INIV=1
CALL FACILE(IREP,INUM,'SURF',INIV,'TEMPERATURE ',ZTSL1,LLCOSP)
CALL FACILE(IREP,INUM,'PROF',INIV,'TEMPERATURE ',ZTPL1,LLCOSP)

!     ------------------------------------------------------------------

!     4. LOOP ON THE ARPEGE GRID.
!        ------------------------

!     4.1 Standard grid

I=0
IM=0
IT=0
DO JJ=1,NDGLG
  DO J=1,NLOENG(JJ)
!  Addresses
    I=I+1
    IA=J+NSTAGP(JJ)-ID

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
    ZZDL= ZDL/YDGSGEOM_NB%GM(IA)
    ZLAT= YDGSGEOM_NB%GELAT(IA)*ZCV
    ILAT= NINT((YRCLI%ELATNE-ZLAT)/YRCLI%EDLAT)+1
    ZDLA= ZLIS*ZZDL/YRCLI%EDLAT
    IDLA= INT(ZDLA)
    IF (ILAT < (1-IDLA) .OR. ILAT > (YRCLI%NDATY+IDLA)) GO TO 410
    I1LA= MAX(1    ,ILAT-IDLA)
    I2LA= MIN(YRCLI%NDATY,ILAT+IDLA)

!  Check the longitude
    ZLON= YDGSGEOM_NB%GELAM(IA)*ZCV
    ILON= NINT((ZLON-YRCLI%ELONSW)/YRCLI%EDLON)+1
    ZDLO= ZLIS*ZZDL/YRCLI%EDLON/COS(YDGSGEOM_NB%GELAT(IA))
    IDLO= MIN(INT(ZDLO),YRCLI%NDATX/4)
    I1LO= MOD(ILON-IDLO-1+2*YRCLI%NGLOBX,YRCLI%NGLOBX)+1
    I2LO= MOD(ILON+IDLO-1+2*YRCLI%NGLOBX,YRCLI%NGLOBX)+1
    IF (I1LO > YRCLI%NDATX .AND. I2LO >= I1LO) GO TO 410
    IF (I1LO > YRCLI%NDATX) I1LO=1
    I2LO=MIN(YRCLI%NDATX,I2LO)

!  Initialize new fields
    ZMSK= 0.0_JPRB
    ZDPS= 0.0_JPRB
    ZTSL= 0.0_JPRB
    ZTPL= 0.0_JPRB
    DO JT=1,JPTYVE
      ZITP(JT)= 0.0_JPRB
    ENDDO

!  Compute new values
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
      ENDDO
    ENDDO
!  End of the two loops in the rectangular latitude x longitude box
!  ZMSK : Number of non-missing data in the box
!  ZPTL : Maximum number of land points in the box 
!         (including the part of the box outside the domain)
    IF (ZMSK < ZEPS) GO TO 410
    IM=IM+1
    ZPTL=REAL((2*IDLA+1)*(2*IDLO+1),JPRB)*(1.0_JPRB-ZLND1(I))
    ZPTL=MAX(ZMSK,ZPTL)

!  Compute new fields
!  Choice of the dominant type of vegetation. The dominant vegetation must
!  exceed 25%, otherwise the old type is kept.
    ZITPN= ZITP1(I)
    ZITPM= PPONDER*ZPTL
    DO JT=1,JPTYVE
      IF ((JT == YRCLI%NTPMER.OR.JT == YRCLI%NTPLAC).AND.(ZITP(JT) >= ZITPM)) THEN
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
!  Soil depth
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

IF (IM == 0) CALL ABOR1(' ERROR IN INCLI7 : NO MODIFICATION !')
WRITE(NULOUT,'(I6,''/'',I6,'' POINTS MODIFIED BY INCLI7'')') IM,IT
WRITE(NULOUT,'(I2,'' FIELDS MODIFIED BY INCLI7'')') IC

!     ------------------------------------------------------------------

!     5. WRITING THE ARPEGE FILE.
!        ------------------------

INIV=0
CALL FAIENC(IREP,INUM,'SURF',INIV,'IND.VEG.DOMI',ZITP1,LLCOSP)
CALL FAIENC(IREP,INUM,'SURF',INIV,'EPAIS.SOL   ',ZDEP1,LLCOSP)
CALL FAIENC(IREP,INUM,'SURF',INIV,'EPAI.SOL.MAX',ZDPS1,LLCOSP)
CALL FAIENC(IREP,INUM,'SURF',INIV,'ALBEDO      ',ZALB1,LLCOSP)
CALL FAIENC(IREP,INUM,'SURF',INIV,'EMISSIVITE  ',ZEMI1,LLCOSP)
INIV=1
CALL FAIENC(IREP,INUM,'SURF',INIV,'TEMPERATURE ',ZTSL1,LLCOSP)
CALL FAIENC(IREP,INUM,'PROF',INIV,'TEMPERATURE ',ZTPL1,LLCOSP)

CALL CHK923(YDGEOMETRY,INUM)

CALL LFILAF(IREP,INUM,.FALSE.)
CALL FAIRME(IREP,INUM,'KEEP')

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INCLI7',1,ZHOOK_HANDLE)
END SUBROUTINE INCLI7
