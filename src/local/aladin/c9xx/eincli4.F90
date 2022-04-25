SUBROUTINE EINCLI4(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF)

!**** *EINCLI4*

!     PURPOSE
!     -------

!     This routine calculates monthly climatological fields.

!**   INTERFACE
!     ---------

!     CALL EINCLI4

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
!     The target grid is a regular grid (geographical or in plane projection).

!     Input data are global, and coming from units 11 to 15 :
!     ALL input fields are organized from 0 to 360(E) deg in longitude
!                                    from North to South in latitude
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
!     Coherence between fields computed in EINCLI1,EINCLI2,EINCLI3 and EINCLI4
!     is checked. The output fields are added to or modified in the ALADIN file.

!     EXTERNALS
!     ---------

!     CCHIEN
!     EGEO923
!     EINTER8
!     EINTER10
!     EBICLI (previously EBIEN)
!     ELECI
!     FA-LFI package (FAITOU,LFILAF,FAIRME)
!     ECHK923

!     AUTHORS
!     -------
!     L. Gerard 3/06/1997 from INCLI4 and EINCLIC

!     MODIFICATIONS
!     -------------
!     R. El Khatib : 01-12-06 Cleaning sm variables
!     D. Giard 01-11-23 add lake type
!     D. Giard 01-11-30 corrections for defaults
!     M.Hamrud 03-10-01 CY28 Cleaning
!     D. Giard 02-11-29 case of no water; mask checking depending on resolution
!     D. Giard 04-09-15 update
!     D. Giard 05-04-07 : no packing for roughness lengths, new EBICLI
!     F. Taillefer 09-06-02 add LZ0THER
!     O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!------------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE YOEPHY       , ONLY : TEPHY
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCLI   , ONLY : YRCLI 
USE YOMCST   , ONLY : RG
USE YOMLUN   , ONLY : NULOUT
USE YOM_YGFL , ONLY : TYPE_GFLD

!------------------------------------------------------------------------

IMPLICIT NONE

! JPNVAR    : number of input datasets
! JPNVAR + 3: number of fields to be derived (monthly surface characteristics)
! JPBX/JPBY : number of extra longitudes/latitudes on each side of the data
!          grid required by interpolation (4 or 12 points -> 2/2)

TYPE(GEOMETRY), INTENT(INOUT)   :: YDGEOMETRY
TYPE(TYPE_GFLD)    ,INTENT(INOUT):: YDGFL
TYPE(TEPHY)    ,INTENT(INOUT)   :: YDEPHY
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
INTEGER(KIND=JPIM) :: JPBX
INTEGER(KIND=JPIM) :: JPBY
INTEGER(KIND=JPIM) :: JPNVAR

PARAMETER(JPNVAR=5,JPBX=2,JPBY=2)

! Input grid geometry:
! ZFLD needs 2 empty rows / columns at each border for interpolation routines:
REAL(KIND=JPRB) :: ZDTA(YRCLI%NDATX,YRCLI%NDATY,2),ZDES(YRCLI%NDATX,YRCLI%NDATY)
REAL(KIND=JPRB) :: ZFLD(YRCLI%NDATX+2*JPBX,YRCLI%NDATY+2*JPBY,0:2)
REAL(KIND=JPRD) ,ALLOCATABLE :: ZAUX(:)
! Final grid geometry (C+I):
REAL(KIND=JPRB) :: ZRES((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1),&
 & -3:JPNVAR+3)&
 & ,ZCMP((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1),0:JPNVAR)&
 & ,ZMSK((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1))  
REAL(KIND=JPRB) ,ALLOCATABLE :: ZCMS(:,:)
! Total grid (C+I+E) for EBICLI
REAL(KIND=JPRB) :: ZEXT(YDGEOMETRY%YRDIM%NDLON*YDGEOMETRY%YRDIM%NDGLG,0:JPNVAR+3)
REAL(KIND=JPRB) :: ZDEF(JPNVAR)
REAL(KIND=JPRB) :: Z0, ZA, ZEPS, ZEPX, ZFZ, ZL, ZMIS, ZR, ZU, ZV, ZVI
REAL(KIND=JPRB) :: ZHOOK_HANDLE

CHARACTER :: CLPREF(-3:JPNVAR+3)*8,CLSUFF(-3:JPNVAR+3)*12
CHARACTER :: CLFORM*12,CLNOMC*16,CLNOMF*10

INTEGER(KIND=JPIM) :: INIVL(0:JPNVAR+3)
INTEGER(KIND=JPIM) :: IARI, IARP, IC, ICF, ICN, ICO, IDATX, IDATY, IINF,&
 & IMES, IMSKI, IMSKR, INIV, INUM, IREP, ITFING, IV,&
 & IX, IXFING, IY, IYFING, IZ, J, JF, JJ  

LOGICAL :: LLBIP(0:JPNVAR+3),LLWRI(0:JPNVAR+3),LLPAC(0:JPNVAR+3)
LOGICAL :: LLIMST

!------------------------------------------------------------------------

#include "cchien.intfb.h"
#include "ebicli.intfb.h"
#include "echk923.intfb.h"
#include "egeo923.intfb.h"
#include "einter8.intfb.h"
#include "einter10.intfb.h"
#include "eleci.intfb.h"

!------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EINCLI4',0,ZHOOK_HANDLE)

ASSOCIATE( &
 & NDGUNG=>YDGEOMETRY%YRDIM%NDGUNG, NDGUXG=>YDGEOMETRY%YRDIM%NDGUXG, NDLUNG=>YDGEOMETRY%YRDIM%NDLUNG, &
 & NDLUXG=>YDGEOMETRY%YRDIM%NDLUXG, NDGLG=>YDGEOMETRY%YRDIM%NDGLG, NDLON=>YDGEOMETRY%YRDIM%NDLON, &
 & EDELX=>YDGEOMETRY%YREGEO%EDELX, EDELY=>YDGEOMETRY%YREGEO%EDELY)
!------------------------------------------------------------------------

!     1. SET INITIAL VALUES.
!        -------------------

!     1.1 Dimensions for data and interpolation.

ICO=YRCLI%NPINT*YRCLI%NPINT
! initial grid:
IDATX=YRCLI%NDATX+2*JPBX
IDATY=YRCLI%NDATY+2*JPBY
! final grid:
IXFING=NDLUXG-NDLUNG+1
IYFING=NDGUXG-NDGUNG+1
ITFING=IXFING*IYFING

!     1.2 Type of data files.

IF (YRCLI%LIEEE) THEN
  CLFORM='UNFORMATTED'
ELSE
  CLFORM='FORMATTED'
ENDIF

!     1.3 Geometry : initialize YEMCLI and YOMDIL

CALL EGEO923(YDGEOMETRY%YRGEM)

!     1.4 Final biperiodization and writing

DO J=0,JPNVAR+3
  LLBIP(J)=.TRUE.
  LLWRI(J)=.TRUE.
  LLPAC(J)=.TRUE.
  INIVL(J)=0
ENDDO
LLBIP(0)=.FALSE.
LLWRI(0)=.FALSE.

DO JF=0,JPNVAR+3
  DO J=1,NDGLG*NDLON
    ZEXT(J,JF)=0.0_JPRB
  ENDDO
ENDDO

!     1.5 Default values for interpolated fields (desert)

! Vegetation fraction
ZDEF(1)=0.0_JPRB
! Roughness length of vegetation (Log(z) interpolated)
ZDEF(2)=LOG(YRCLI%SZZ0D)
! Albedo of vegetation
ZDEF(3)=YRCLI%SALBD
! Leaf Area Index
ZDEF(4)=0.0_JPRB
! Minimum surface resistance (LAI/Rsm interpolated)
ZDEF(5)=0.0_JPRB

!     1.6 Miscellaneous.

ZFZ= YRCLI%SFCZ0*RG
ZMIS= YRCLI%SMANQ - 1.0_JPRB
ZEPS=0.1_JPRB * MAX(YRCLI%EDLAT,YRCLI%EDLON)
ZEPX=1.E-10_JPRB

!     ------------------------------------------------------------------

!     2. READ INITIAL ALADIN CLIM FILE.
!        ------------------------------

! Setup
INUM=3
IREP=0
IARI=0
IARP=27
IMES=1
CLNOMF='Const.Clim'
CLNOMC='Const.Clim.Surfa'
LLIMST=.TRUE.

IINF=-1
INIV=0

! Open file
CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'OLD',.TRUE.,LLIMST,&
 & IMES,IARP,IARI,CLNOMC)  
CALL CCHIEN(YDGEOMETRY,CLNOMC,INUM,IINF)

! Read useful fields, including initial values for some ones
CLSUFF(-3)='PROP.URBANIS'
CLSUFF(-2)='PROP.TERRE  '
CLSUFF(-1)='IND.VEG.DOMI'
CLSUFF( 0)='IND.TERREMER'
CLSUFF( 1)='PROP.VEGETAT'
CLSUFF( 2)='Z0.FOIS.G   '
CLSUFF( 3)='ALBEDO      '

! Reload fields on C+I
DO J=-3,3
  CLPREF(J)='SURF    '
  CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,CLPREF(J),INIV,CLSUFF(J),ZRES(1,J),NULOUT)
ENDDO

! Put land-sea mask on C+I domain for EBICLI
ZEXT(1:ITFING,0)=ZRES(1:ITFING,0)

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
  DO JJ=YRCLI%NDATY,1,-1
    DO J=1,YRCLI%NDATX
      IZ=J+(YRCLI%NDATY-JJ)*YRCLI%NDATX
      ZDTA(J,JJ,1)= ZAUX(IZ)
    ENDDO
  ENDDO
ELSE
  READ(11,*) ((ZDTA(J,JJ,1),J=1,YRCLI%NDATX),JJ=YRCLI%NDATY,1,-1)
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
DO J=1,ITFING
  ZMSK(J)= 1.0_JPRB-ZRES(J,0)
ENDDO

!    3.3 Interpolations

! Interpolation of mask (missing or suspect data)
CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,0),IDATY,IDATX,1,ZCMP(1,0),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

! Intersection between mask and domain
IMSKR=0
IMSKI=0
DO J=1,ITFING
  IF (ZCMP(J,0) > ZEPS) IMSKR= IMSKR+1
  IF (ZCMP(J,0) > 1.0_JPRB-ZEPS) IMSKI= IMSKI+1
ENDDO

! No mask over the ALADIN domain : interpolation of data without mask
IF (IMSKR == 0) THEN
  WRITE(NULOUT,'('' NO EFFECTIVE MASK - VEG'')')
  CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,1),IDATY,IDATX,ICN,ZCMP(1,ICF),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

! No available data over the ALADIN domain
ELSEIF (IMSKI == ITFING) THEN
  WRITE(NULOUT,'('' NO AVAILABLE DATA FOR THIS AREA ! - VEG'')')
  DO JF=ICF,ICF+ICN-1
    DO J=1,ITFING
      ZCMP(J,JF)= ZDEF(JF)
    ENDDO
  ENDDO

! Standard case : interpolation of data with a mask
ELSE
  ALLOCATE ( ZCMS(ITFING,0:ICN) )
  CALL EINTER10(YDGEOMETRY%YREGEO,ZFLD(1,1,1),IDATY,IDATX,ICN,ZCMP(1,ICF),ZCMS(1,1),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT,&
   & ZFLD(1,1,0),ZCMP(1,0),ZMSK,ZCMS(1,0),YRCLI%SMASK)  
  DO JF=ICF,ICF+ICN-1
    IC=JF-ICF+1
    DO J=1,ITFING
      IF (ZCMP(J,JF) <= YRCLI%SMANQ) ZCMP(J,JF)= ZCMS(J,IC)
    ENDDO
  ENDDO
  DEALLOCATE ( ZCMS )
ENDIF

!     ------------------------------------------------------------------

!     4. READ AND INTERPOLATE DATA : LOG OF ROUGHNESS LENGTH.
!        ----------------------------------------------------

ICF=2
ICN=1

!    4.1 Read input data

IF (YRCLI%LIEEE) THEN
  ALLOCATE ( ZAUX(YRCLI%NDATX*YRCLI%NDATY) )
ENDIF

OPEN(UNIT=12,FILE='z0v_GL',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  READ(12) ZAUX
  DO JJ=YRCLI%NDATY,1,-1
    DO J=1,YRCLI%NDATX
      IZ=J+(YRCLI%NDATY-JJ)*YRCLI%NDATX
      ZDTA(J,JJ,1)= ZAUX(IZ)
    ENDDO
  ENDDO
ELSE
  READ(12,*) ((ZDTA(J,JJ,1),J=1,YRCLI%NDATX),JJ=YRCLI%NDATY,1,-1)
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
    IF (ZDTA(IX,IY,1) < ZEPX .OR. ZDES(IX,IY) >= YRCLI%SMASK) THEN
      ZFLD(J,JJ,0)= 1.0_JPRB
      ZFLD(J,JJ,1)= ZMIS
    ELSE
      ZFLD(J,JJ,0)= 0.0_JPRB
      ZFLD(J,JJ,1)= LOG(ZDTA(IX,IY,1))
    ENDIF
  ENDDO
ENDDO

! Output grid (sea + desert masks)
DO J=1,ITFING
  ZMSK(J)= 1.0_JPRB-ZRES(J,0)*MAX(0.0_JPRB,SIGN(1.0_JPRB,ZRES(J,1)-YRCLI%SVEG))
ENDDO

!    4.3 Interpolations

! Interpolation of mask (missing or suspect data)
CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,0),IDATY,IDATX,1,ZCMP(1,0),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

! Intersection between mask and domain
IMSKR=0
IMSKI=0
DO J=1,ITFING
  IF (ZCMP(J,0) > ZEPS) IMSKR= IMSKR+1
  IF (ZCMP(J,0) > 1.0_JPRB-ZEPS) IMSKI= IMSKI+1
ENDDO

! No mask over the ALADIN domain : interpolation of data without mask
IF (IMSKR == 0) THEN
  WRITE(NULOUT,'('' NO EFFECTIVE MASK - Z0V'')')
  CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,1),IDATY,IDATX,ICN,ZCMP(1,ICF),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

! No available data over the ALADIN domain
ELSEIF (IMSKI == ITFING) THEN
  WRITE(NULOUT,'('' NO AVAILABLE DATA FOR THIS AREA ! - Z0V'')')
  DO JF=ICF,ICF+ICN-1
    DO J=1,ITFING
      ZCMP(J,JF)= ZDEF(JF)
    ENDDO
  ENDDO

! Standard case : interpolation of data with a mask
ELSE
  ALLOCATE ( ZCMS(ITFING,0:ICN) )
  CALL EINTER10(YDGEOMETRY%YREGEO,ZFLD(1,1,1),IDATY,IDATX,ICN,ZCMP(1,ICF),ZCMS(1,1),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT,&
   & ZFLD(1,1,0),ZCMP(1,0),ZMSK,ZCMS(1,0),YRCLI%SMASK)  
  DO JF=ICF,ICF+ICN-1
    IC=JF-ICF+1
    DO J=1,ITFING
      IF (ZCMP(J,JF) <= YRCLI%SMANQ) ZCMP(J,JF)= ZCMS(J,IC)
    ENDDO
  ENDDO
  DEALLOCATE ( ZCMS )
ENDIF
                               
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
  DO JJ=YRCLI%NDATY,1,-1
    DO J=1,YRCLI%NDATX
      IZ=J+(YRCLI%NDATY-JJ)*YRCLI%NDATX
      ZDTA(J,JJ,1)= ZAUX(IZ)
    ENDDO
  ENDDO
ELSE
  READ(13,*) ((ZDTA(J,JJ,1),J=1,YRCLI%NDATX),JJ=YRCLI%NDATY,1,-1)
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

! Output grid (sea + desert masks, no change in ZMSK)

!    5.3 Interpolations

! Interpolation of mask (missing or suspect data)
CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,0),IDATY,IDATX,1,ZCMP(1,0),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

! Intersection between mask and domain
IMSKR=0
IMSKI=0
DO J=1,ITFING
  IF (ZCMP(J,0) > ZEPS) IMSKR= IMSKR+1
  IF (ZCMP(J,0) > 1.0_JPRB-ZEPS) IMSKI= IMSKI+1
ENDDO

! No mask over the ALADIN domain : interpolation of data without mask
IF (IMSKR == 0) THEN
  WRITE(NULOUT,'('' NO EFFECTIVE MASK - ALV'')')
  CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,1),IDATY,IDATX,ICN,ZCMP(1,ICF),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

! No available data over the ALADIN domain
ELSEIF (IMSKI == ITFING) THEN
  WRITE(NULOUT,'('' NO AVAILABLE DATA FOR THIS AREA ! - ALV'')')
  DO JF=ICF,ICF+ICN-1
    DO J=1,ITFING
      ZCMP(J,JF)= ZDEF(JF)
    ENDDO
  ENDDO

! Standard case : interpolation of data with a mask
ELSE
  ALLOCATE ( ZCMS(ITFING,0:ICN) )
  CALL EINTER10(YDGEOMETRY%YREGEO,ZFLD(1,1,1),IDATY,IDATX,ICN,ZCMP(1,ICF),ZCMS(1,1),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT,&
   & ZFLD(1,1,0),ZCMP(1,0),ZMSK,ZCMS(1,0),YRCLI%SMASK)  
  DO JF=ICF,ICF+ICN-1
    IC=JF-ICF+1
    DO J=1,ITFING
      IF (ZCMP(J,JF) <= YRCLI%SMANQ) ZCMP(J,JF)= ZCMS(J,IC)
    ENDDO
  ENDDO
  DEALLOCATE ( ZCMS )
ENDIF
                               
!     ------------------------------------------------------------------

!     6. READ AND INTERPOLATE DATA : LEAF AREA INDEX, MIN. SURF. RESISTANCE.
!        -------------------------------------------------------------------

ICF=4
ICN=2

!    6.1 Read input data

IF (YRCLI%LIEEE) THEN
  ALLOCATE ( ZAUX(YRCLI%NDATX*YRCLI%NDATY) )
ENDIF

! Leaf area index
OPEN(UNIT=14,FILE='lai_GL',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  READ(14) ZAUX
  DO JJ=YRCLI%NDATY,1,-1
    DO J=1,YRCLI%NDATX
      IZ=J+(YRCLI%NDATY-JJ)*YRCLI%NDATX
      ZDTA(J,JJ,1)= ZAUX(IZ)
    ENDDO
  ENDDO
ELSE
  READ(14,*) ((ZDTA(J,JJ,1),J=1,YRCLI%NDATX),JJ=YRCLI%NDATY,1,-1)
ENDIF
CLOSE(14)

! Minimum surface resistance
OPEN(UNIT=15,FILE='rsm_GL',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  READ(15) ZAUX
  DO JJ=YRCLI%NDATY,1,-1
    DO J=1,YRCLI%NDATX
      IZ=J+(YRCLI%NDATY-JJ)*YRCLI%NDATX
      ZDTA(J,JJ,2)= ZAUX(IZ)
    ENDDO
  ENDDO
ELSE
  READ(15,*) ((ZDTA(J,JJ,2),J=1,YRCLI%NDATX),JJ=YRCLI%NDATY,1,-1)
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
    IF (ZDTA(IX,IY,1) <= YRCLI%SMANQ .OR. ZDTA(IX,IY,2) <= ZEPX&
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
CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,0),IDATY,IDATX,1,ZCMP(1,0),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

! Intersection between mask and domain
IMSKR=0
IMSKI=0
DO J=1,ITFING
  IF (ZCMP(J,0) > ZEPS) IMSKR= IMSKR+1
  IF (ZCMP(J,0) > 1.0_JPRB-ZEPS) IMSKI= IMSKI+1
ENDDO

! No mask over the ALADIN domain : interpolation of data without mask
IF (IMSKR == 0) THEN
  WRITE(NULOUT,'('' NO EFFECTIVE MASK - LAI,RSM'')')
  CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,1),IDATY,IDATX,ICN,ZCMP(1,ICF),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

! No available data over the ALADIN domain
ELSEIF (IMSKI == ITFING) THEN
  WRITE(NULOUT,'('' NO AVAILABLE DATA FOR THIS AREA ! - LAI,RSM'')')
  DO JF=ICF,ICF+ICN-1
    DO J=1,ITFING
      ZCMP(J,JF)= ZDEF(JF)
    ENDDO
  ENDDO

! Standard case : interpolation of data with a mask
ELSE
  ALLOCATE ( ZCMS(ITFING,0:ICN) )
  CALL EINTER10(YDGEOMETRY%YREGEO,ZFLD(1,1,1),IDATY,IDATX,ICN,ZCMP(1,ICF),ZCMS(1,1),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT,&
   & ZFLD(1,1,0),ZCMP(1,0),ZMSK,ZCMS(1,0),YRCLI%SMASK)  
  DO JF=ICF,ICF+ICN-1
    IC=JF-ICF+1
    DO J=1,ITFING
      IF (ZCMP(J,JF) <= YRCLI%SMANQ) ZCMP(J,JF)= ZCMS(J,IC)
    ENDDO
  ENDDO
  DEALLOCATE ( ZCMS )
ENDIF

!     ------------------------------------------------------------------

!     7. MODIFY (1-3) OR CREATE (4-8) OF ALADIN FIELDS.
!        ----------------------------------------------

DO J=1,ITFING

!    7.1 Initial and limit values

!  Vegetation fraction (corrected with land fraction)
  ZVI= ZCMP(J,1)
  ZVI= MAX(0.0_JPRB,MIN(1.0_JPRB,ZVI))
  ZV = MIN(ZRES(J,1),ZVI*ZRES(J,-2))

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
  ZR= ZL/MAX(ZR,ZEPX)

!    7.2 Corrections according to land-use type
 
!  Set values over sea (no vegetation, albedo and roughness unchanged)
  IV= NINT(ZRES(J,-1))
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
    ZU = ZRES(J,-3)
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

!     8. BIPERIODIZE AND WRITE FIELDS IN CLIM FILE
!     ------------------------------------------------

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

LLPAC(2)=.FALSE.
LLPAC(6)=.FALSE.
LLPAC(7)=.FALSE.

DO JF=1,JPNVAR+3
  DO J=1,ITFING
    ZEXT(J,JF)=ZRES(J,JF)
  ENDDO
ENDDO

ICN=JPNVAR+4
CALL EBICLI(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF,ICN,INIVL,CLPREF(0),CLSUFF(0),INUM,ZEXT,LLBIP,LLWRI,LLPAC)

CALL ECHK923(YDGEOMETRY%YRDIM,INUM)

CALL LFILAF(IREP,INUM,.TRUE.)
CALL FAIRME(IREP,INUM,'KEEP')

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EINCLI4',1,ZHOOK_HANDLE)
END SUBROUTINE EINCLI4
