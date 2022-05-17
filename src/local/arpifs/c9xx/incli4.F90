SUBROUTINE INCLI4(YDGEOMETRY)

!**** *INCLI4*

!     PURPOSE
!     -------

!     This routine calculates monthly climatological fields.

!**   INTERFACE
!     ---------

!     CALL INCLI4

!     METHOD
!     ------

!     This routine calculates 8 monthly fields :
!     - fraction of vegetation cover (1)
!     - g*cinetic roughness length   (2)
!     - albedo                       (3)
!     - leaf area index              (4)
!     - minimum surface resistance   (5)
!     - g*thermal roughness length   (6)
!     - g*cinetic roughness length : contribution from vegetation (7)
!     - albedo : contribution from vegetation                     (8)
!     The target grid is a Gaussian grid with rotation of the pole and "Schmidt"
!     stretching as well as -if asked- reduction of the number of grid points
!     towards the pole of the final representation.
!     Input data are global, and coming from units 11 to 15 :
!     - fraction of vegetation cover                (file veg_GL, unit 11)
!     - roughness length of vegetation              (file z0v_Gl, unit 12)
!     - albedo of vegetation                        (file alv_Gl, unit 13)
!     - leaf area index                             (file lai_GL, unit 14)
!     - minimum surface resistance                  (file rsm_GL, unit 15)
!     or from fixed climatological fields :
!     - maximum vegetation fraction        (initial value of 1)
!     - bare ground albedo (complementary) (------- ----- -- 2)
!     - g*roughness length of topography   (------- ----- -- 3)
!     - land use type
!     Fields are calculated with a mask, taking into account missing data.
!     A 4-point interpolation operator is used for all fields. The values on
!     the Gaussian grid are the averages of the values in the boxes (a box is
!     a part of the subgrid corresponding to a point of the Gaussian grid).
!     Ln(z0) and LAI/Rsmin are interpolated rather than z0 and Rsmin.
!     The fraction of vegetation cover is weighted by the fraction of land.
!     Vegetation characteristics are taken into account only as long as the
!     fraction of vegetation cover keeps significant.
!     The ratio of thermal to cinetic roughness length is fixed to 1/10 .
!     Coherence between fields computed in INCLI1,INCLI2,INCLI3 and INCLI4 is
!     checked. The output fields are added to or modified in the ARPEGE file.

!     EXTERNALS
!     ---------

!     GEO923
!     INTER8
!     INTER10
!     FA-LFI package (FAITOU,FACILE,FAIENC,LFILAF,FAIRME)
!     CHIEN
!     ABOR1
!     CHK923

!     AUTHORS
!     -------
!      D. Giard 97-03-13  from INCLIC and INCLICM

!     MODIFICATIONS
!     -------------
!      D. Giard 01-11-15 reduced memory cost; better consistency between masks
!      D. Giard 01-11-23 add lake type
!      R. El Khatib 03-08-18 Roughness lengths not packed
!      M.Hamrud 01-10-03 CY28 Cleaning
!      D. Giard 04-09-15 cleaning
!      F. Taillefer 09-06-02 add LZ0THER
!      K. Yessad (Jan 2010): externalisation of group EGGX in XRD/IFSAUX
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCT0   , ONLY : NQUAD
USE YOMLUN   , ONLY : NULOUT
USE YOMCST   , ONLY : RG
USE YOMVERT  , ONLY : VP00
USE YOMCLI   , ONLY : YRCLI  

IMPLICIT NONE

! JPNVAR    : number of input datasets
! JPNVAR + 3: number of fields to be derived (monthly surface characteristics)
! JPBX/JPBY : number of extra longitudes/latitudes on each side of the data
!          grid required by interpolation (4 or 12 points -> 2/2)

TYPE(GEOMETRY), INTENT(IN)   :: YDGEOMETRY
INTEGER(KIND=JPIM) :: JPBX
INTEGER(KIND=JPIM) :: JPBY
INTEGER(KIND=JPIM) :: JPNVAR

PARAMETER(JPNVAR=5,JPBX=2,JPBY=2)

REAL(KIND=JPRB) :: ZDTA(YRCLI%NDATX,YRCLI%NDATY,2),ZDES(YRCLI%NDATX,YRCLI%NDATY)
REAL(KIND=JPRB) :: ZFLD(YRCLI%NDATX+2*JPBX,YRCLI%NDATY+2*JPBY,0:2)
REAL(KIND=JPRB) :: ZRES(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON,-2:JPNVAR+3),&
 & ZCMP(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON,0:JPNVAR)&
 & ,ZLSM(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON),ZMSK(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON)&
 & ,ZBO(YDGEOMETRY%YRDIM%NDGLG+1),ZLONG(YDGEOMETRY%YRDIM%NDGLG),ZMU(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON),&
 & ZSLA(YDGEOMETRY%YRDIM%NDGLG*YRCLI%NPINT)&
 & ,ZSLO(YDGEOMETRY%YRDIM%NDLON*YRCLI%NPINT,YDGEOMETRY%YRDIM%NDGLG),ZCLO(YDGEOMETRY%YRDIM%NDLON*YRCLI%NPINT,YDGEOMETRY%YRDIM%NDGLG)  
REAL(KIND=JPRB) ,ALLOCATABLE :: ZAUX(:)
REAL(KIND=JPRB) ,ALLOCATABLE :: ZCMS(:,:)
REAL(KIND=JPRB) :: Z0, ZA, ZEPS, ZFZ, ZL, ZMIS, ZR, ZU, ZV, ZVI
REAL(KIND=JPRB) :: ZHOOK_HANDLE

CHARACTER :: CLPREF(-2:JPNVAR+3)*8,CLSUFF(-2:JPNVAR+3)*12
CHARACTER :: CLNOMC*16,CLNOMF*10,CLFORM*12

INTEGER(KIND=JPIM) :: IARI, IARP, IC, ICF, ICN, ICO, IDATX, IDATY, IINF,&
 & IIQ, IJQ, ILENE, ILENT, IMES, INIV,&
 & INUM, IREP, IV, IX, IY, IZ, J, JF, JJ  
INTEGER(KIND=JPIM) :: INGRIG,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL

LOGICAL :: LLCOSP, LLIMST, LLPOLE
LOGICAL :: LLPACK(-2:JPNVAR+3)

!     ------------------------------------------------------------------

#include "chien.h"

#include "abor1.intfb.h"
#include "chk923.intfb.h"
#include "geo923.intfb.h"
#include "inter10.intfb.h"
#include "inter8.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INCLI4',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDVAB=>YDGEOMETRY%YRVAB&
& )
ASSOCIATE(NDGENG=>YDDIM%NDGENG, NDGLG=>YDDIM%NDGLG, NDGSAG=>YDDIM%NDGSAG,   NDLON=>YDDIM%NDLON, NSMAX=>YDDIM%NSMAX,   &
& NFLEVG=>YDDIMV%NFLEVG,   NHTYP=>YDGEM%NHTYP, NLOENG=>YDGEM%NLOENG, NMENG=>YDGEM%NMENG,   NSTTYP=>YDGEM%NSTTYP,      &
& RLOCEN=>YDGEM%RLOCEN, RMUCEN=>YDGEM%RMUCEN,   RSTRET=>YDGEM%RSTRET)
!     ------------------------------------------------------------------

!     1. SET INITIAL VALUES.
!        -------------------

!     1.1 Dimensions for data and interpolation.

ICO=YRCLI%NPINT*YRCLI%NPINT
IJQ=NDGLG*YRCLI%NPINT
IIQ=NDLON*YRCLI%NPINT
ILENT=NDGLG*NDLON
IDATX=YRCLI%NDATX+2*JPBX
IDATY=YRCLI%NDATY+2*JPBY

!     1.2 Type of data files.

IF (YRCLI%LIEEE) THEN
  CLFORM='UNFORMATTED'
ELSE
  CLFORM='FORMATTED'
ENDIF

!     1.3 Model grid.

ILENE=0
CALL GEO923(YDGEOMETRY,YRCLI%NPINT,ILENE,ZBO,ZLONG,ZMU,ZSLA,ZSLO,ZCLO)

!     1.4 Miscellaneous.

ZFZ= YRCLI%SFCZ0*RG
ZMIS= YRCLI%SMANQ - 1.0_JPRB
LLPACK(:)=.TRUE.

!     ------------------------------------------------------------------

!     2. READING INITIAL ARPEGE CLIM FILE.
!        ---------------------------------

!    2.1 Setup

INUM=3
IREP=0
IARI=0
IARP=27
IMES=1
CLNOMF='Const.Clim'
CLNOMC='Const.Clim.Surfa'
LLIMST=.TRUE.

IINF=-1
ZEPS=1.E-10_JPRB
LLPOLE=.TRUE.

INIV=0
LLCOSP=.FALSE.

!    2.2 Open file

CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'OLD',.TRUE.,LLIMST,&
 & IMES,IARP,IARI,CLNOMC)  
CALL CHIEN(CLNOMC,NSTTYP,RMUCEN,RLOCEN,RSTRET,NSMAX,&
 & NDGLG,NDLON,NLOENG,NMENG,NHTYP,NFLEVG,VP00,YDVAB%VALH,YDVAB%VBH,&
 & NQUAD,IINF,NDGSAG,NDGENG,ZEPS,LLPOLE,NULOUT)  
IF (LLPOLE) CALL ABOR1(' CLIM. FILES MUST NOT HAVE POLES !')

!    2.3 Read useful fields, including initial values for some ones

CLSUFF(-2)='PROP.URBANIS'
CLSUFF(-1)='PROP.TERRE  '
CLSUFF( 0)='IND.VEG.DOMI'
CLSUFF( 1)='PROP.VEGETAT'
CLSUFF( 2)='Z0.FOIS.G   '
CLSUFF( 3)='ALBEDO      '
LLPACK( 2)=.FALSE.
DO J=-2,3
  CLPREF(J)='SURF    '
  CALL FACILE(IREP,INUM,CLPREF(J),INIV,CLSUFF(J),ZRES(1,J),LLCOSP)
ENDDO
CALL FACILE(IREP,INUM,'SURF',INIV,'IND.TERREMER',ZLSM,LLCOSP)

!     ------------------------------------------------------------------

!     3. READ AND INTERPOLATE DATA : VEGETATION FRACTION.
!        ------------------------------------------------

ICF=1
ICN=1

!    3.1 Read input data

IF (YRCLI%LIEEE) THEN
  ALLOCATE ( ZAUX(YRCLI%NDATX*YRCLI%NDATY) )
ENDIF

OPEN(UNIT=11,FILE='veg_GL',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  READ(11) ZAUX
  DO JJ=1,YRCLI%NDATY
    DO J=1,YRCLI%NDATX
      IZ=J+(JJ-1)*YRCLI%NDATX
      ZDTA(J,JJ,1)=ZAUX(IZ)
    ENDDO
  ENDDO
ELSE
  READ(11,*) ((ZDTA(J,JJ,1),J=1,YRCLI%NDATX),JJ=1,YRCLI%NDATY)
ENDIF
CLOSE(11)

IF (YRCLI%LIEEE) THEN
  DEALLOCATE ( ZAUX )
ENDIF

!    3.2 Define masks

! Input grid (missing data + desert)
DO JJ=JPBY+1,JPBY+YRCLI%NDATY
  DO J=JPBX+1,JPBX+YRCLI%NDATX
    IX=J -JPBX
    IY=JJ-JPBY
    IF (ZDTA(IX,IY,1) <= YRCLI%SMANQ) THEN
      ZFLD(J,JJ,0)= 1.0_JPRB
    ELSE
      ZFLD(J,JJ,0)= 0.0_JPRB
    ENDIF
    ZFLD(J,JJ,1)= ZDTA(IX,IY,1)
    IF (ZDTA(IX,IY,1) < YRCLI%SVEG) THEN
      ZDES(IX,IY)= 1.0_JPRB
    ELSE
      ZDES(IX,IY)= 0.0_JPRB
    ENDIF
  ENDDO
ENDDO

! Output grid (sea mask)
DO J=1,ILENE
  ZMSK(J)= 1.0_JPRB-ZLSM(J)
ENDDO

!     3.3 Interpolations

! Interpolation of mask (missing or suspect data)
CALL INTER8(NDGLG,NDLON,IDATY,IDATX,1,YRCLI%NPINT,ILENE,ICO,IJQ,IIQ,&
 & NLOENG(1:),ILENT,ZCMP(1,0),ZFLD(1,1,0),&
 & ZSLA,ZSLO,ZCLO)  

! Interpolation of fields, with mask
ALLOCATE ( ZCMS(ILENT,0:ICN) )
CALL INTER10(NDGLG,NDLON,IDATY,IDATX,ICN,YRCLI%NPINT,ILENE,ICO,IJQ,IIQ,&
 & NLOENG(1:),ILENT,ZCMP(1,ICF),ZCMS(1,1),ZFLD(1,1,1),&
 & ZSLA,ZSLO,ZCLO,&
 & ZFLD(1,1,0),ZCMP(1,0),ZMSK,ZCMS(1,0),YRCLI%SMASK)  
DO JF=ICF,ICF+ICN-1
  IC=JF-ICF+1
  DO J=1,ILENE
    IF (ZCMP(J,JF) <= YRCLI%SMANQ) ZCMP(J,JF)= ZCMS(J,IC)
  ENDDO
ENDDO
DEALLOCATE ( ZCMS )

!     ------------------------------------------------------------------

!     4. READ AND INTERPOLATE DATA : ROUGHNESS LENGTH.
!        ---------------------------------------------

ICF=2
ICN=1

!    4.1 Read input data

IF (YRCLI%LIEEE) THEN
  ALLOCATE ( ZAUX(YRCLI%NDATX*YRCLI%NDATY) )
ENDIF

OPEN(UNIT=12,FILE='z0v_GL',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  READ(12) ZAUX
  DO JJ=1,YRCLI%NDATY
    DO J=1,YRCLI%NDATX
      IZ=J+(JJ-1)*YRCLI%NDATX
      ZDTA(J,JJ,1)=ZAUX(IZ)
    ENDDO
  ENDDO
ELSE
  READ(12,*) ((ZDTA(J,JJ,1),J=1,YRCLI%NDATX),JJ=1,YRCLI%NDATY)
ENDIF
CLOSE(12)

IF (YRCLI%LIEEE) THEN
  DEALLOCATE ( ZAUX )
ENDIF

!    4.2 Define masks

! Input grid (missing data + desert)
DO JJ=JPBY+1,JPBY+YRCLI%NDATY
  DO J=JPBX+1,JPBX+YRCLI%NDATX
    IX=J -JPBX
    IY=JJ-JPBY
    IF (ZDTA(IX,IY,1) < ZEPS .OR. ZDES(IX,IY) >= YRCLI%SMASK) THEN
      ZFLD(J,JJ,0)= 1.0_JPRB
      ZFLD(J,JJ,1)= ZMIS
    ELSE
      ZFLD(J,JJ,0)= 0.0_JPRB
      ZFLD(J,JJ,1)= LOG(ZDTA(IX,IY,1))
    ENDIF
  ENDDO
ENDDO

! Output grid (sea + desert mask)
DO J=1,ILENE
  ZMSK(J)= 1.0_JPRB-ZLSM(J)*MAX(0.0_JPRB,SIGN(1.0_JPRB,ZRES(J,1)-YRCLI%SVEG))
ENDDO

!     4.3 Interpolations

! Interpolation of mask (missing or suspect data)
CALL INTER8(NDGLG,NDLON,IDATY,IDATX,1,YRCLI%NPINT,ILENE,ICO,IJQ,IIQ,&
 & NLOENG(1:),ILENT,ZCMP(1,0),ZFLD(1,1,0),&
 & ZSLA,ZSLO,ZCLO)  

! Interpolation of fields, with mask
ALLOCATE ( ZCMS(ILENT,0:ICN) )
CALL INTER10(NDGLG,NDLON,IDATY,IDATX,ICN,YRCLI%NPINT,ILENE,ICO,IJQ,IIQ,&
 & NLOENG(1:),ILENT,ZCMP(1,ICF),ZCMS(1,1),ZFLD(1,1,1),&
 & ZSLA,ZSLO,ZCLO,&
 & ZFLD(1,1,0),ZCMP(1,0),ZMSK,ZCMS(1,0),YRCLI%SMASK)  
DO JF=ICF,ICF+ICN-1
  IC=JF-ICF+1
  DO J=1,ILENE
    IF (ZCMP(J,JF) <= YRCLI%SMANQ) ZCMP(J,JF)= ZCMS(J,IC)
  ENDDO
ENDDO
DEALLOCATE ( ZCMS )

!     ------------------------------------------------------------------

!     5. READ AND INTERPOLATE DATA : ALBEDO.
!        -----------------------------------

ICF=3
ICN=1

!    5.1 Read input data

IF (YRCLI%LIEEE) THEN
  ALLOCATE ( ZAUX(YRCLI%NDATX*YRCLI%NDATY) )
ENDIF

OPEN(UNIT=13,FILE='alv_GL',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  READ(13) ZAUX
  DO JJ=1,YRCLI%NDATY
    DO J=1,YRCLI%NDATX
      IZ=J+(JJ-1)*YRCLI%NDATX
      ZDTA(J,JJ,1)=ZAUX(IZ)
    ENDDO
  ENDDO
ELSE
  READ(13,*) ((ZDTA(J,JJ,1),J=1,YRCLI%NDATX),JJ=1,YRCLI%NDATY)
ENDIF
CLOSE(13)

IF (YRCLI%LIEEE) THEN
  DEALLOCATE ( ZAUX )
ENDIF

!    5.2 Define masks

! Input grid (missing data + desert)
DO JJ=JPBY+1,JPBY+YRCLI%NDATY
  DO J=JPBX+1,JPBX+YRCLI%NDATX
    IX=J -JPBX
    IY=JJ-JPBY
    IF (ZDTA(IX,IY,1) <= YRCLI%SMANQ .OR. ZDES(IX,IY) >= YRCLI%SMASK) THEN
      ZFLD(J,JJ,0)= 1.0_JPRB
    ELSE
      ZFLD(J,JJ,0)= 0.0_JPRB
    ENDIF
    ZFLD(J,JJ,1)= ZDTA(IX,IY,1)
  ENDDO
ENDDO

! Output grid (sea + desert mask, no change)

!     5.3 Interpolations

! Interpolation of mask (missing or suspect data)
CALL INTER8(NDGLG,NDLON,IDATY,IDATX,1,YRCLI%NPINT,ILENE,ICO,IJQ,IIQ,&
 & NLOENG(1:),ILENT,ZCMP(1,0),ZFLD(1,1,0),&
 & ZSLA,ZSLO,ZCLO)  

! Interpolation of fields, with mask
ALLOCATE ( ZCMS(ILENT,0:ICN) )
CALL INTER10(NDGLG,NDLON,IDATY,IDATX,ICN,YRCLI%NPINT,ILENE,ICO,IJQ,IIQ,&
 & NLOENG(1:),ILENT,ZCMP(1,ICF),ZCMS(1,1),ZFLD(1,1,1),&
 & ZSLA,ZSLO,ZCLO,&
 & ZFLD(1,1,0),ZCMP(1,0),ZMSK,ZCMS(1,0),YRCLI%SMASK)  
DO JF=ICF,ICF+ICN-1
  IC=JF-ICF+1
  DO J=1,ILENE
    IF (ZCMP(J,JF) <= YRCLI%SMANQ) ZCMP(J,JF)= ZCMS(J,IC)
  ENDDO
ENDDO
DEALLOCATE ( ZCMS )

!     ------------------------------------------------------------------

!     6. READ AND INTERPOLATE DATA : LEAF AREA INDEX, MIN. SURF. RESISTANCE.
!        -------------------------------------------------------------------

ICF=4
ICN=2

!    6.1 Read input data

IF (YRCLI%LIEEE) THEN
  ALLOCATE ( ZAUX(YRCLI%NDATX*YRCLI%NDATY) )
ENDIF

OPEN(UNIT=14,FILE='lai_GL',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  READ(14) ZAUX
  DO JJ=1,YRCLI%NDATY
    DO J=1,YRCLI%NDATX
      IZ=J+(JJ-1)*YRCLI%NDATX
      ZDTA(J,JJ,1)=ZAUX(IZ)
    ENDDO
  ENDDO
ELSE
  READ(14,*) ((ZDTA(J,JJ,1),J=1,YRCLI%NDATX),JJ=1,YRCLI%NDATY)
ENDIF
CLOSE(14)

OPEN(UNIT=15,FILE='rsm_GL',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  READ(15) ZAUX
  DO JJ=1,YRCLI%NDATY
    DO J=1,YRCLI%NDATX
      IZ=J+(JJ-1)*YRCLI%NDATX
      ZDTA(J,JJ,2)=ZAUX(IZ)
    ENDDO
  ENDDO
ELSE
  READ(15,*) ((ZDTA(J,JJ,2),J=1,YRCLI%NDATX),JJ=1,YRCLI%NDATY)
ENDIF
CLOSE(15)

IF (YRCLI%LIEEE) THEN
  DEALLOCATE ( ZAUX )
ENDIF

!    6.2 Define masks

! Input grid (missing data + desert)
DO JJ=JPBY+1,JPBY+YRCLI%NDATY
  DO J=JPBX+1,JPBX+YRCLI%NDATX
    IX=J -JPBX
    IY=JJ-JPBY
    IF (ZDTA(IX,IY,1) <= YRCLI%SMANQ .OR. ZDTA(IX,IY,2) <= ZEPS&
       & .OR. ZDES(IX,IY) >= YRCLI%SMASK) THEN  
      ZFLD(J,JJ,0)= 1.0_JPRB
      ZFLD(J,JJ,1)= ZMIS
      ZFLD(J,JJ,2)= ZMIS
    ELSE
      ZFLD(J,JJ,0)= 0.0_JPRB
      ZFLD(J,JJ,1)= ZDTA(IX,IY,1)
      ZFLD(J,JJ,2)= ZDTA(IX,IY,1)/ZDTA(IX,IY,2)
    ENDIF
  ENDDO
ENDDO

! Output grid (sea + desert masks, no change in ZMSK)

!    6.3 Interpolations

! Interpolation of mask (missing or suspect data)
CALL INTER8(NDGLG,NDLON,IDATY,IDATX,1,YRCLI%NPINT,ILENE,ICO,IJQ,IIQ,&
 & NLOENG(1:),ILENT,ZCMP(1,0),ZFLD(1,1,0),&
 & ZSLA,ZSLO,ZCLO)  

! Interpolation of fields, with mask
ALLOCATE ( ZCMS(ILENT,0:ICN) )
CALL INTER10(NDGLG,NDLON,IDATY,IDATX,ICN,YRCLI%NPINT,ILENE,ICO,IJQ,IIQ,&
 & NLOENG(1:),ILENT,ZCMP(1,ICF),ZCMS(1,1),ZFLD(1,1,1),&
 & ZSLA,ZSLO,ZCLO,&
 & ZFLD(1,1,0),ZCMP(1,0),ZMSK,ZCMS(1,0),YRCLI%SMASK)  
DO JF=ICF,ICF+ICN-1
  IC=JF-ICF+1
  DO J=1,ILENE
    IF (ZCMP(J,JF) <= YRCLI%SMANQ) ZCMP(J,JF)= ZCMS(J,IC)
  ENDDO
ENDDO
DEALLOCATE ( ZCMS )

!     ------------------------------------------------------------------

!     7. MODIFY (1-3) OR CREATE (4-8) OF ARPEGE FIELDS.
!        ----------------------------------------------

DO J=1,ILENE

!    7.1 Initial and limit values

!  Vegetation fraction (corrected with land fraction)
  ZVI= ZCMP(J,1)
  ZVI= MAX(0.0_JPRB,MIN(1.0_JPRB,ZVI))
  ZV= MIN(ZRES(J,1),ZVI*ZRES(J,-1))

!  Roughness length of vegetation
  Z0= ZCMP(J,2)
  Z0= EXP(Z0)

!  Albedo of vegetation
  ZA= ZCMP(J,3)
  ZA= MIN(YRCLI%SALBX,MAX(YRCLI%SALBN,ZA))

!  Leaf area index
  ZL= ZCMP(J,4)
  ZL= MAX(0.0_JPRB,ZL)

!  Minimum surface resistance
  ZR= ZCMP(J,5)
  ZR= ZL/MAX(ZR,ZEPS)

!    7.2 Corrections according to land-use type
 
!  Set values over sea (no vegetation, albedo and roughness unchanged)
  IV= NINT(ZRES(J,0))
  IF (IV == YRCLI%NTPMER .OR. IV == YRCLI%NTPLAC) THEN
    ZRES(J,1)= 0.0_JPRB
    ZRES(J,4)= 0.0_JPRB
    ZRES(J,5)= YRCLI%SRSMX
    ZRES(J,6)= ZRES(J,2)
    ZRES(J,7)= 0.0_JPRB
    ZRES(J,8)= YRCLI%SALBN
  ELSE

!  Correct if vegetation cover is not significant
    IF (IV == YRCLI%NTPGLA .OR. ZVI < YRCLI%SVEG) THEN
      ZV = 0.0_JPRB
      Z0 = YRCLI%SZZ0D
      ZA = YRCLI%SALBN
      ZL = 0.0_JPRB
      ZR = YRCLI%SRSMD
    ENDIF

!  Roughness length :
!  contributions from vegetation, urban areas and desert ; scaling
    ZU = ZRES(J,-2)
    Z0 = ZFZ*(ZV*Z0 + ZU*YRCLI%SZZ0U + (1.0_JPRB-ZU-ZV)*YRCLI%SZZ0D)

!  Final values
    ZRES(J,1)= ZV
    ZRES(J,2)= SQRT(ZRES(J,2)**2 + Z0**2)
    ZRES(J,3)= MIN(YRCLI%SALBX,MAX(YRCLI%SALBN,(1.0_JPRB-ZV)*ZRES(J,3)+ZV*ZA))
    ZRES(J,4)= ZL
    ZRES(J,5)= MIN(YRCLI%SRSMX,MAX(YRCLI%SRSMN,ZR))
    IF (YRCLI%LZ0THER) THEN
      ZRES(J,6)= ZRES(J,2)*YRCLI%STHER
    ELSE
      ZRES(J,6)= Z0*YRCLI%STHER
    ENDIF
    ZRES(J,7)= Z0
    ZRES(J,8)= ZA
  ENDIF

ENDDO

!     ------------------------------------------------------------------

!     4. WRITING FIELDS IN ARPEGE CLIM FILE, THEN CHECKING.
!     -----------------------------------------------------

CLPREF(4)='SURF    '
CLPREF(5)='SURF    '
CLPREF(6)='SURF    '
CLPREF(7)='SURF    '
CLPREF(8)='SURF    '
CLSUFF(4)='IND.FOLIAIRE'
CLSUFF(5)='RESI.STO.MIN'
CLSUFF(6)='GZ0.THERM   '
CLSUFF(7)='Z0VEG.FOIS.G'
CLSUFF(8)='ALBEDO.VEG  '

LLPACK(6)=.FALSE.
LLPACK(7)=.FALSE.

DO J=1,JPNVAR+3
  IF (.NOT.LLPACK(J)) THEN
!   Do not pack roughness lengths
    CALL FAVEUR(IREP,INUM,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
    INGRIG=0
    CALL FAGOTE(IREP,INUM,INGRIG,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  ENDIF
  CALL FAIENC(IREP,INUM,CLPREF(J),INIV,CLSUFF(J),ZRES(1,J),LLCOSP)
  IF (.NOT.LLPACK(J)) THEN
!   Reset packing after the treatment of roughness lengths
    CALL FAGOTE(IREP,INUM,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  ENDIF
ENDDO

CALL CHK923(YDGEOMETRY,INUM)

CALL LFILAF(IREP,INUM,.TRUE.)
CALL FAIRME(IREP,INUM,'KEEP')

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INCLI4',1,ZHOOK_HANDLE)
END SUBROUTINE INCLI4
