SUBROUTINE EINCLI2(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF)

!**** *EINCLI2*

!     PURPOSE
!     -------

!     This routine calculates fixed climatological fields.

!**   INTERFACE
!     ---------

!     CALL EINCLI2

!     METHOD
!     ------

!     This routine calculates 8 fixed fields :
!     - dominant land use type
!     - bare ground albedo
!     - emissivity
!     - maximum depth of the soil column
!     - percentage of clay
!     - percentage of sand
!     - useful depth of the soil column
!     - maximum vegetation fraction

!     The target grid is a regular grid (geographical or in plane projection).
!     Input data are global, and coming from units 11 to 18 :
!     ALL input fields are organized from 0 to 360(E) deg in longitude
!                                    from North to South in latitude
!     - dominant land use type                          (file itp_GL, unit 11)
!     - bare ground albedo                              (file alb_GL, unit 12)
!     - emissivity                                      (file emi_GL, unit 13)
!     - hydrological depth of soil, corresponding percentages of clay and sand
!                        (files dps_GL, arg_GL, sab_GL units 14, 15, 16 resp.)
!     - max. vegetation cover along the annual cycle    (file vgx_GL, unit 17)
!     - root depth                                      (file dpr_GL, unit 18)

!     Fields are calculated with a mask, taking into account missing data. The
!     same mask is imposed for the following fields:
!     - bare ground albedo and emissivity
!     - hydrological depth of soil, percentages of clay and sand
!     Data for sea, lake and ice-cap are assumed to be missing (undefined) for 
!     soil depth and texture as well as for vegetation characteristics.
!     The case of no missing data is considered for bare ground albedo and 
!     emissivity (no mask).
!     Information from the fraction of vegetation cover is used for root depth.
!     The fraction of vegetation cover is weighted by the fraction of land.
!     A 4-point interpolation operator is used for all fields. The values on
!     the Gaussian grid are the averages of the values in the boxes (a box is
!     a part of the subgrid corresponding to a point of the Gaussian grid).
!     The output fields are biperiodicised and added to the ARPEGE file.

!     When old-fashioned fields are required (case LSOLV=.FALSE.), only the
!     following fields are computed or modified :
!     - albedo
!     - emissivity
!     - vegetation cover
!     - g*roughness length
!       (adding the contributions from urbanization and vegetation)

!     EXTERNALS
!     ---------

!     ABOR1
!     CCHIEN
!     EGEO923
!     EINTER8
!     EINTER10
!     EBICLI (previously EBIEN)
!     ELECI
!     FA-LFI package (FAITOU,LFILAF,FAIENC,FAIRME)

!     AUTHORS
!     -------
!     L. Gerard 30/05/97 from INCLI2, EINCLIA and EINCLIC

!     MODIFICATIONS
!     -------------
!     R. El Khatib : 01-12-06 Cleaning sm variables
!     D. Giard 01-11-13 more allocatable arrays to save memory
!     M.Hamrud 03-10-01 CY28 Cleaning
!     D. Giard 02-11-29 case of no water; mask checking depending on resolution
!     D. Giard 04-09-15 update
!     D. Giard 05-04-07 : no packing for roughness lengths, new EBICLI
!     O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE YOEPHY       , ONLY : TEPHY
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCLI   , ONLY : YRCLI
USE YOMCST   , ONLY : RG
USE YOMLUN   , ONLY : NULOUT
USE YOM_YGFL , ONLY : TYPE_GFLD

!     ------------------------------------------------------------------

IMPLICIT NONE

! JPTYVE : number of land-use types, at least 3 :
!   1 -> sea , 2 -> ice , 3 -> desert or low vegetation (lakes are 1 or 5)
! JPNFIX : number of fixed fields to be derived
! JPNF   : maximum number of fields to be interpolated together
! JPBX/JPBY : number of extra longitudes/latitudes on each side of the data
!          grid required by interpolation (4 or 12 points -> 2/2)

TYPE(GEOMETRY), INTENT(INOUT)   :: YDGEOMETRY
TYPE(TYPE_GFLD)    ,INTENT(INOUT):: YDGFL
TYPE(TEPHY)    ,INTENT(INOUT)   :: YDEPHY
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
INTEGER(KIND=JPIM) :: JPBX
INTEGER(KIND=JPIM) :: JPBY
INTEGER(KIND=JPIM) :: JPNF
INTEGER(KIND=JPIM) :: JPNFIX
INTEGER(KIND=JPIM) :: JPTYVE

PARAMETER(JPTYVE=5,JPNFIX=8,JPNF=MAX(JPTYVE,4),JPBX=2,JPBY=2)

REAL(KIND=JPRB),ALLOCATABLE :: ZFLD (:,:,:)
REAL(KIND=JPRD),ALLOCATABLE :: ZAUX1(:),ZAUX2(:),ZAUX3(:),ZITP(:)
REAL(KIND=JPRB),ALLOCATABLE :: ZIVD (:),ZIVS(:) 
REAL(KIND=JPRB) :: ZRES((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1),&
 & -1:JPNFIX+1)
REAL(KIND=JPRB) :: ZCMP((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1),&
 & JPNF)&
 & ,ZCMS((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1),JPNF)&
 & ,ZLS((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1))  
REAL(KIND=JPRB) :: ZEXT(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON,-1:JPNFIX+4)
REAL(KIND=JPRB) :: ZD, ZEPS, ZFZ0, ZMAX1, ZMAX2, ZSUM, ZV, ZZ0R, ZZARG, ZZSAB
REAL(KIND=JPRB) :: ZHOOK_HANDLE

CHARACTER :: CLPREF(-1:JPNFIX+4)*8,CLSUFF(-1:JPNFIX+4)*12
CHARACTER :: CLFORM*12,CLNOMC*16,CLNOMF*10

LOGICAL :: LLBIP(0:JPNFIX+4),LLWRI(0:JPNFIX+4),LLPAC(0:JPNFIX+4) 

INTEGER(KIND=JPIM) :: INIVL(0:JPNFIX+4)
INTEGER(KIND=JPIM) :: IARI, IARP, IC, ICO, IDATX, IDATY, IFLD, IINF,&
 & IK, IMES, IMSK, IMSKI, IMSKR, INDEX1, INDEX2, INIV, INUM, IREP,&
 & ITFING, ILSM, ITP, IV, IVD, IVN, IVS, IX, IXFING, IY, IYFING, IZ,&
 & J, JC, JJ, JK, JT, JX, JY  

LOGICAL :: LLCOSP, LLIMST

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "cchien.intfb.h"
#include "ebicli.intfb.h"
#include "egeo923.intfb.h"
#include "einter8.intfb.h"
#include "einter10.intfb.h"
#include "eleci.intfb.h"

!     ------------------------------------------------------------------
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EINCLI2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDPHY=>YDML_PHY_MF%YRPHY)
ASSOCIATE(EDELY=>YDGEOMETRY%YREGEO%EDELY, EDELX=>YDGEOMETRY%YREGEO%EDELX, &
 & NDGLG=>YDDIM%NDGLG, NDLUXG=>YDDIM%NDLUXG, NDGUNG=>YDDIM%NDGUNG, &
 & NDGUXG=>YDDIM%NDGUXG, NDLUNG=>YDDIM%NDLUNG, NDLON=>YDDIM%NDLON, &
 & LSOLV=>YDPHY%LSOLV)
!     1. SET INITIAL VALUES.
!        -------------------

!     1.1 Dimensions for data and interpolation.

IDATX=YRCLI%NDATX+2*JPBX
IDATY=YRCLI%NDATY+2*JPBY

IXFING=NDLUXG-NDLUNG+1
IYFING=NDGUXG-NDGUNG+1
ITFING=IXFING*IYFING
ICO=YRCLI%NPINT*YRCLI%NPINT

ZEPS=0.1_JPRB * MAX(YRCLI%EDLAT,YRCLI%EDLON)

!     1.2 Type of data files.

IF (YRCLI%LIEEE) THEN
  CLFORM='UNFORMATTED'
ELSE
  CLFORM='FORMATTED'
ENDIF

!     1.3 Geometry : initialize YEMCLI and YOMDIL

CALL EGEO923(YDGEOMETRY%YRGEM)

!     1.4 Names of ARPEGE fields

DO J=-1,JPNFIX+4
  CLPREF(J)='SURF'
ENDDO
CLSUFF(-1)='PROP.TERRE  '
CLSUFF( 0)='IND.TERREMER'
CLSUFF( 1)='IND.VEG.DOMI'
CLSUFF( 2)='ALBEDO.SOLNU'
CLSUFF( 3)='EMISSIVITE  '
CLSUFF( 4)='EPAI.SOL.MAX'
CLSUFF( 5)='PROP.ARGILE '
CLSUFF( 6)='PROP.SABLE  '
CLSUFF( 7)='PROP.VEG.MAX'
CLSUFF( 8)='EPAIS.SOL   '
CLSUFF( 9)='PROP.URBANIS'
CLSUFF(10)='ALBEDO      '
CLSUFF(11)='PROP.VEGETAT'
CLSUFF(12)='ALBEDO.COMPL'
IF (.NOT. LSOLV) THEN
  CLSUFF( 5)='Z0.FOIS.G   '
  CLSUFF( 6)='PROP.URBANIS'
ENDIF

!     1.5 Arrays for final biperiodization and writing

DO J=0,JPNFIX+4
  LLBIP(J)=.TRUE.
  LLWRI(J)=.TRUE.
  LLPAC(J)=.TRUE.
  INIVL(J)=0
ENDDO
LLBIP(0)=.FALSE.
LLBIP(1)=.FALSE.
LLWRI(0)=.FALSE.
IF (.NOT. LSOLV) LLPAC(5)=.FALSE.

DO JC=0,JPNFIX+4
  DO J=1,NDGLG*NDLON
    ZEXT(J,JC)=0.0_JPRB
  ENDDO
ENDDO

!     1.6 Open ARPEGE File and read land-sea mask and land fraction

INUM=3
IREP=0
IARI=0
IARP=10
IMES=1
CLNOMF='Const.Clim'
CLNOMC='Const.Clim.Surfa'
LLIMST=.TRUE.
IINF=-1

INIV=0
LLCOSP=.FALSE.

CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'OLD',.TRUE.,LLIMST,&
 & IMES,IARP,IARI,CLNOMC)  

! Control CADRE compatibility:
CALL CCHIEN(YDGEOMETRY,CLNOMC,INUM,IINF)

! Reload the 2 needed fields: on C+I domain only
DO J=-1,0
  CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,CLPREF(J),INIV,CLSUFF(J),ZRES(1,J),NULOUT)
ENDDO
! ZLS fraction of sea for interpolation routines
! ZEXT(.,0) land sea mask for EBICLI
DO J=1,ITFING
  ZLS(J)=1.0_JPRB-ZRES(J,0)
  ZEXT(J,0)=ZRES(J,0)
ENDDO

! Case of OLD_FASHIONED stuff
IF (.NOT. LSOLV) THEN
  CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,CLPREF(5),INIV,CLSUFF(5),ZRES(1,5),NULOUT)
  CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,CLPREF(6),INIV,CLSUFF(6),ZRES(1,6),NULOUT)
  GOTO 300
ELSE
  CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,CLPREF(9),INIV,CLSUFF(9),ZRES(1,9),NULOUT)
ENDIF

!     ------------------------------------------------------------------

!     2. DOMINANT LAND USE TYPE.
!        -----------------------

!     2.1 Allocate temporary space

ALLOCATE ( ZFLD (IDATX,IDATY,JPTYVE) )
IF (YRCLI%LIEEE) THEN
  ALLOCATE ( ZITP (YRCLI%NDATX*YRCLI%NDATY) )
ELSE
  ALLOCATE ( ZITP (YRCLI%NDATX) )
ENDIF

!     2.2 Read data

DO JJ=JPBY+1,JPBY+YRCLI%NDATY
  DO J=JPBX+1,JPBX+YRCLI%NDATX
    DO JT=1,JPTYVE
      ZFLD(J,JJ,JT)=0.0_JPRB
    ENDDO
  ENDDO
ENDDO

OPEN(UNIT=11,FILE='itp_GL',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  READ(11)ZITP
  DO JK=1,YRCLI%NDATY
    JJ=YRCLI%NDATY+JPBY-JK+1
    DO J=1,YRCLI%NDATX
      IX=J+JPBX
      IZ=J+(JK-1)*YRCLI%NDATX
      ITP=NINT(ZITP(IZ))
      IF ((ITP <= 0) .OR. (ITP > JPTYVE))&
       & CALL ABOR1(' EINCLI2 : INVALID LAND-USE TYPE ! -1-')  
      ZFLD(IX,JJ,ITP)=1.0_JPRB
    ENDDO
  ENDDO
ELSE
  DO JJ=YRCLI%NDATY+JPBY,1+JPBY,-1
    READ(11,*) (ZITP(J),J=1,YRCLI%NDATX)
    DO J=1,YRCLI%NDATX
      IX=J+JPBX
      ITP=NINT(ZITP(J))
      IF ((ITP <= 0) .OR. (ITP > JPTYVE))&
       & CALL ABOR1(' EINCLI2 : INVALID LAND-USE TYPE ! -1-')  
      ZFLD(IX,JJ,ITP)=1.0_JPRB
    ENDDO
  ENDDO
ENDIF
CLOSE(11)

!     2.3 Interpolation -> proportion of each type

IFLD=JPTYVE
CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,1),IDATY,IDATX,IFLD,ZCMP(1,1),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!     2.4 Deallocate/allocate temporary space

DEALLOCATE ( ZITP )
DEALLOCATE ( ZFLD )

ALLOCATE ( ZIVD (ITFING) )
ALLOCATE ( ZIVS (ITFING) )

!     2.5 Determination of dominant and secondary types

DO J=1,ITFING
  ZMAX1=0.0_JPRB
  ZMAX2=0.0_JPRB
  ZIVD(J)=0.0_JPRB
  ZIVS(J)=0.0_JPRB
  DO JT=1,JPTYVE
    IF (ZCMP(J,JT) > ZMAX1) THEN
      ZMAX1=ZCMP(J,JT)
      ZIVD(J)=REAL(JT,JPRB)
    ENDIF
  ENDDO
  DO JT=1,JPTYVE
    IF ((ZCMP(J,JT) < ZMAX1).AND.(ZCMP(J,JT) > ZMAX2)) THEN
      ZMAX2=ZCMP(J,JT)
      ZIVS(J)=REAL(JT,JPRB)
    ENDIF
  ENDDO
ENDDO

!     2.6 Adjustment of the dominant type of vegetation and the land/sea index

DO J=1,ITFING
  IVD=NINT(ZIVD(J))
  IVS=NINT(ZIVS(J))
  ILSM=NINT(ZRES(J,0))
! Water -> Dominant type must be 1 (sea) or 5  (lake)
  IF (ILSM == 0 .AND. IVD /= YRCLI%NTPLAC) IVD=YRCLI%NTPMER
! Land / Dominant type is 1 or 5 -> take secondary type
  IF ((IVD == YRCLI%NTPMER .OR. IVD == YRCLI%NTPLAC) .AND. (ILSM == 1)) THEN
    IVD=IVS
  ENDIF
! Land / -> research of a neighbour if
! No dominant type (0) or dominant type is 1 and no secondary type (0)
  IF (IVD == 0) THEN
    IC=0
    IX=1
    IY=J+1
    IK=MAX(1,MIN(ITFING,IY))
    IVN=NINT(ZIVD(IK))
262 CONTINUE 
    IF ((IVN /= YRCLI%NTPMER).AND.(IVN /= YRCLI%NTPLAC).AND.(IVN /= 0)) GO TO 263
    IC=IC+1
    IF (IC > ITFING) GO TO 264
    IX=SIGN(IABS(IX)+1,-IX)
    IY=IY+IX
    IK=MAX(1,MIN(ITFING,IY))
    IVN=NINT(ZIVD(IK))
    GOTO 262
263 CONTINUE
    IVD=IVN
! Desert taken as default value if the search is unsuccessful
264 CONTINUE
    IF (IVD == 0) IVD=YRCLI%NTPDES
    WRITE(NULOUT,FMT=&
     & '('' DOMINANT TYPE OF VEGET. FORCED TO'',I2,'', J='',I6)')&
     & IVD,J  
  ENDIF
  ZIVD(J)=REAL(IVD,JPRB)
ENDDO

!     2.8 ALADIN field

DO J=1,ITFING
  ZRES(J,1)=ZIVD(J)
ENDDO

!     2.9 Compute globally dominant land use type; recycle array zivd for this

DO JT=1,JPTYVE
  ZIVD(JT)=0.0_JPRB
ENDDO
DO JT=1,JPTYVE
  DO J=1,ITFING
    ZIVD(JT)=ZIVD(JT)+ZCMP(J,JT)
  ENDDO
ENDDO
IVD=0
ZD=0.0_JPRB
DO JT=1,JPTYVE
  IF (ZD < ZIVD(JT)) THEN
    IVD=JT
    ZD=ZIVD(JT)
  ENDIF
ENDDO
WRITE(NULOUT,*) 'GLOBALLY DOMINANT LAND USE IS ',IVD

! Now this type IVD has to be used to fill the Extension zone

INDEX1=0
INDEX2=0
DO JY=1,IYFING
  DO JX=1,IXFING
    ZEXT(JX+INDEX2,1)=ZRES(JX+INDEX1,1)
  ENDDO
  DO JX=IXFING+1,NDLON
    ZEXT(JX+INDEX2,1)=IVD
  ENDDO
  INDEX1=INDEX1+IXFING
  INDEX2=INDEX2+NDLON
ENDDO
DO JY=IYFING+1,NDGLG
  DO JX=1,NDLON
    ZEXT(JX+INDEX2,1)=IVD
  ENDDO
  INDEX2=INDEX2+NDLON
ENDDO

!     2.10 Deallocate temporary space

DEALLOCATE ( ZIVS )
DEALLOCATE ( ZIVD )

!     ------------------------------------------------------------------

!     3. BARE GROUND ALBEDO AND EMISSIVITY.
!        ----------------------------------

300 CONTINUE

!     3.1 Allocate temporary space

ALLOCATE ( ZFLD (IDATX,IDATY,3) )
IF (YRCLI%LIEEE) THEN
  ALLOCATE (ZAUX1(YRCLI%NDATX*YRCLI%NDATY))
  ALLOCATE (ZAUX2(YRCLI%NDATX*YRCLI%NDATY))
ENDIF

!     3.2 Read data

OPEN(UNIT=12,FILE='alb_GL',FORM=CLFORM)
OPEN(UNIT=13,FILE='emi_GL',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  READ(12)ZAUX1
  READ(13)ZAUX2
  DO JJ=YRCLI%NDATY+JPBY,1+JPBY,-1
    DO J=1,YRCLI%NDATX
      IX=J+JPBX
      IZ=J+(YRCLI%NDATY+JPBY-JJ)*YRCLI%NDATX
      ZFLD(IX,JJ,2)=ZAUX1(IZ)
      ZFLD(IX,JJ,3)=ZAUX2(IZ)
    ENDDO
  ENDDO
ELSE
  READ(12,*)((ZFLD(J,JJ,2),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=YRCLI%NDATY+JPBY,1+JPBY,-1)
  READ(13,*)((ZFLD(J,JJ,3),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=YRCLI%NDATY+JPBY,1+JPBY,-1)
ENDIF
CLOSE(12)
CLOSE(13)

!     3.3 Define the mask (missing data)

IMSK=0
DO JJ=JPBY+1,JPBY+YRCLI%NDATY
  DO J=JPBX+1,JPBX+YRCLI%NDATX
    IF (MIN(ZFLD(J,JJ,2),ZFLD(J,JJ,3)) <= YRCLI%SMANQ) THEN
      ZFLD(J,JJ,1)=1.0_JPRB
      IMSK=IMSK+1
    ELSE
      ZFLD(J,JJ,1)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

!     3.4 Interpolate mask and fields

!  No mask in the global dataset : interpolation of data without mask
IF (IMSK == 0) THEN
  WRITE(NULOUT,'('' NO MASK FOR BARE GROUND ALBEDO AND EMISSIVITY !'')')
  IFLD=2
  CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,2),IDATY,IDATX,IFLD,ZCMP(1,2),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  Only missing data
ELSEIF (IMSK == YRCLI%NDATX*YRCLI%NDATY) THEN
  CALL ABOR1(' EINCLI2 : ALBEDO OR EMISSIVITY MISSING !')

!  Standard case : interpolation of mask
ELSE
  IFLD=1
  CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,1),IDATY,IDATX,IFLD,ZCMP(1,1),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  Intersection between mask and domain
  IMSKR=0
  IMSKI=0
  DO J=1,ITFING
    IF (ZCMP(J,1) > ZEPS) IMSKR=IMSKR+1
    IF (ZCMP(J,1) > 1.0_JPRB-ZEPS) IMSKI=IMSKI+1
  ENDDO

!  No mask over the ALADIN domain : interpolation of data without mask
  IF (IMSKR == 0) THEN
    WRITE(NULOUT,&
     & '('' NO EFFECTIVE MASK - BARE GROUND ALBEDO AND EMISSIVITY'')')  
    IFLD=2
    CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,2),IDATY,IDATX,IFLD,ZCMP(1,2),&
     & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  No available data over the ALADIN domain : default values (desert) used
  ELSEIF (IMSKI == ITFING) THEN
    WRITE(NULOUT,'('' ALBEDO AND EMISSIVITY NOT AVAILABLE !'')')
    DO J=1,ITFING
      ZCMP(J,2)= YRCLI%SALBD
      ZCMS(J,2)= YRCLI%SALBD
      ZCMP(J,3)= YRCLI%SEMID
      ZCMS(J,3)= YRCLI%SEMID
    ENDDO

!  Standard case : interpolation of data with a mask
  ELSE
    IFLD=2
    CALL EINTER10(YDGEOMETRY%YREGEO,ZFLD(1,1,2),IDATY,IDATX,IFLD,ZCMP(1,2),ZCMS(1,2),&
     & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT,&
     & ZFLD(1,1,1),ZCMP(1,1),ZLS,ZCMS(1,1),YRCLI%SMASK)  
  ENDIF
ENDIF

!     3.5 Deallocate temporary space

DEALLOCATE ( ZFLD )
IF (YRCLI%LIEEE) THEN
  DEALLOCATE ( ZAUX1 )
  DEALLOCATE ( ZAUX2 )
ENDIF

!     3.7 Move to ARPEGE fields

DO J=1,ITFING
  IF (MIN(ZCMP(J,2),ZCMP(J,3)) <= YRCLI%SMANQ) THEN
    ZRES(J,2)=ZCMS(J,2)
    ZRES(J,3)=ZCMS(J,3)
  ELSE
    ZRES(J,2)=ZCMP(J,2)
    ZRES(J,3)=ZCMP(J,3)
  ENDIF
ENDDO

!     ------------------------------------------------------------------

!     4. DEPTH OF SOIL, PERCENTAGES OF CLAY AND SAND.
!        --------------------------------------------

IF (.NOT. LSOLV) GO TO 500

!     4.1 Allocate temporary space

ALLOCATE ( ZFLD (IDATX,IDATY,4) )
IF (YRCLI%LIEEE) THEN
  ALLOCATE (ZAUX1(YRCLI%NDATX*YRCLI%NDATY))
  ALLOCATE (ZAUX2(YRCLI%NDATX*YRCLI%NDATY))
  ALLOCATE (ZAUX3(YRCLI%NDATX*YRCLI%NDATY))
ENDIF

!     4.2 Read data

OPEN(UNIT=14,FILE='dps_GL',FORM=CLFORM)
OPEN(UNIT=15,FILE='arg_GL',FORM=CLFORM)
OPEN(UNIT=16,FILE='sab_GL',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  READ(14)ZAUX1
  READ(15)ZAUX2
  READ(16)ZAUX3
  DO JJ=YRCLI%NDATY+JPBY,1+JPBY,-1
    DO J=1,YRCLI%NDATX
      IX=J+JPBX
      IZ=J+(YRCLI%NDATY+JPBY-JJ)*YRCLI%NDATX
      ZFLD(IX,JJ,2)=ZAUX1(IZ)
      ZFLD(IX,JJ,3)=ZAUX2(IZ)
      ZFLD(IX,JJ,4)=ZAUX3(IZ)
    ENDDO
  ENDDO
ELSE
  READ(14,*)((ZFLD(J,JJ,2),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=YRCLI%NDATY+JPBY,1+JPBY,-1)
  READ(15,*)((ZFLD(J,JJ,3),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=YRCLI%NDATY+JPBY,1+JPBY,-1)
  READ(16,*)((ZFLD(J,JJ,4),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=YRCLI%NDATY+JPBY,1+JPBY,-1)
ENDIF
CLOSE(14)
CLOSE(15)
CLOSE(16)

!     4.3 Define the mask (missing data)

DO JJ=JPBY+1,JPBY+YRCLI%NDATY
  DO J=JPBX+1,JPBX+YRCLI%NDATX
    IF (MIN(ZFLD(J,JJ,2),ZFLD(J,JJ,3),ZFLD(J,JJ,4)) <= YRCLI%SMANQ)THEN
      ZFLD(J,JJ,1)=1.0_JPRB
    ELSE
      ZFLD(J,JJ,1)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

!     4.4 Interpolate mask and fields

!  Interpolation of mask
IFLD=1
CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,1),IDATY,IDATX,IFLD,ZCMP(1,1),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  Intersection between mask and domain
IMSKR=0
IMSKI=0
DO J=1,ITFING
  IF (ZCMP(J,1) > ZEPS) IMSKR=IMSKR+1
  IF (ZCMP(J,1) > 1.0_JPRB-ZEPS) IMSKI=IMSKI+1
ENDDO

!  No mask over the ALADIN domain : interpolation of data without mask
IF (IMSKR == 0) THEN
  WRITE(NULOUT,'('' NO EFFECTIVE MASK - SOIL DEPTH AND TEXTURE'')')
  IFLD=3
  CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,2),IDATY,IDATX,IFLD,ZCMP(1,2),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  No available data over the ALADIN domain : default values (desert) used

ELSEIF (IMSKI == ITFING) THEN
  WRITE(NULOUT,'('' SOIL DEPTH AND TEXTURE NOT AVAILABLE !'')')
  DO J=1,ITFING
    ZCMP(J,2)= YRCLI%SDEPD
    ZCMS(J,2)= YRCLI%SDEPD
    ZCMP(J,3)= YRCLI%SARGD
    ZCMS(J,3)= YRCLI%SARGD
    ZCMP(J,4)= YRCLI%SSABD
    ZCMS(J,4)= YRCLI%SSABD
  ENDDO

!  Standard case : interpolation of data with a mask
ELSE
  IFLD=3
  CALL EINTER10(YDGEOMETRY%YREGEO,ZFLD(1,1,2),IDATY,IDATX,IFLD,ZCMP(1,2),ZCMS(1,2),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT,&
   & ZFLD(1,1,1),ZCMP(1,1),ZLS,ZCMS(1,1),YRCLI%SMASK)  
ENDIF

!     4.5 Deallocate temporary space

DEALLOCATE ( ZFLD )
IF (YRCLI%LIEEE) THEN
  DEALLOCATE (ZAUX1)
  DEALLOCATE (ZAUX2)
  DEALLOCATE (ZAUX3)
ENDIF

!     4.7 Move to ARPEGE fields

DO J=1,ITFING
  IF (MIN(ZCMP(J,2),ZCMP(J,3),ZCMP(J,4)) <= YRCLI%SMANQ) THEN
    ZRES(J,4)=ZCMS(J,2)
    ZRES(J,5)=ZCMS(J,3)
    ZRES(J,6)=ZCMS(J,4)
  ELSE
    ZRES(J,4)=ZCMP(J,2)
    ZRES(J,5)=ZCMP(J,3)
    ZRES(J,6)=ZCMP(J,4)
  ENDIF
ENDDO

!     ------------------------------------------------------------------

!     5. MAXIMUM VEGETATION COVER AND ROOT DEPTH / ROUGHNESS LENGTH.
!        -----------------------------------------------------------

500 CONTINUE

!     5.1 Allocate temporary space

ALLOCATE ( ZFLD (IDATX,IDATY,4) )
IF (YRCLI%LIEEE) THEN
  ALLOCATE ( ZAUX1(YRCLI%NDATX*YRCLI%NDATY) )
  ALLOCATE ( ZAUX2(YRCLI%NDATX*YRCLI%NDATY) )
ENDIF

!     5.2 Read data

OPEN(UNIT=17,FILE='vgx_GL',FORM=CLFORM)
IF (LSOLV) THEN
  OPEN(UNIT=18,FILE='dpr_GL',FORM=CLFORM)
ELSE
  OPEN(UNIT=18,FILE='z0v_GL',FORM=CLFORM)
ENDIF
IF (YRCLI%LIEEE) THEN
  READ(17)ZAUX1
  READ(18)ZAUX2
  DO JJ=YRCLI%NDATY+JPBY,1+JPBY,-1
    DO J=1,YRCLI%NDATX
      IX=J+JPBX
      IZ=J+(YRCLI%NDATY+JPBY-JJ)*YRCLI%NDATX
      ZFLD(IX,JJ,2)=ZAUX1(IZ)
      ZFLD(IX,JJ,4)=ZAUX2(IZ)
    ENDDO
  ENDDO
ELSE
  READ(17,*)((ZFLD(J,JJ,2),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=YRCLI%NDATY+JPBY,1+JPBY,-1)
  READ(18,*)((ZFLD(J,JJ,4),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=YRCLI%NDATY+JPBY,1+JPBY,-1)
ENDIF
CLOSE(17)
CLOSE(18)

!     5.3 Define the masks (missing data)

DO JJ=JPBY+1,JPBY+YRCLI%NDATY
  DO J=JPBX+1,JPBX+YRCLI%NDATX
    IF (ZFLD(J,JJ,2) <= YRCLI%SMANQ) THEN
      ZFLD(J,JJ,1)=1.0_JPRB
    ELSE
      ZFLD(J,JJ,1)=0.0_JPRB
    ENDIF
    IF ((ZFLD(J,JJ,4) <= YRCLI%SMANQ) .OR. (ZFLD(J,JJ,2) < YRCLI%SVEG)) THEN
      ZFLD(J,JJ,3)=1.0_JPRB
    ELSE
      ZFLD(J,JJ,3)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

!     5.4 Interpolate masks and fields

!  Interpolation of mask for vegetation cover
IFLD=1
CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,1),IDATY,IDATX,IFLD,ZCMP(1,1),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  Intersection between mask and domain
IMSKR=0
IMSKI=0
DO J=1,ITFING
  IF (ZCMP(J,1) > ZEPS) IMSKR=IMSKR+1
  IF (ZCMP(J,1) > 1.0_JPRB- ZEPS) IMSKI=IMSKI+1
ENDDO

!  No mask over the ALADIN domain : interpolation of data without mask
IF (IMSKR == 0) THEN
  WRITE(NULOUT,'('' NO EFFECTIVE MASK - VEGETATION COVER'')')
  IFLD=1
  CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,2),IDATY,IDATX,IFLD,ZCMP(1,2),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  No available data over the ALADIN domain : default values (desert) used
ELSEIF (IMSKI == ITFING) THEN
  WRITE(NULOUT,'('' VEGETATION COVER NOT AVAILABLE !'')')
  DO J=1,ITFING
    ZCMP(J,2)= 0.0_JPRB
    ZCMS(J,2)= 0.0_JPRB
  ENDDO

!  Standard case : interpolation of data with a mask
ELSE
  IFLD=1
  CALL EINTER10(YDGEOMETRY%YREGEO,ZFLD(1,1,2),IDATY,IDATX,IFLD,ZCMP(1,2),ZCMS(1,2),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT,&
   & ZFLD(1,1,1),ZCMP(1,1),ZLS,ZCMS(1,1),YRCLI%SMASK)  
ENDIF

!  Interpolation of mask for vegetation characteristics
IFLD=1
CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,3),IDATY,IDATX,IFLD,ZCMP(1,3),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  Intersection between mask and domain
IMSKR=0
IMSKI=0
DO J=1,ITFING
  IF (ZCMP(J,3) > 0.0_JPRB) IMSKR=IMSKR+1
  IF (ZCMP(J,3) > 1.0_JPRB- ZEPS) IMSKI=IMSKI+1
ENDDO

!  No mask over the ALADIN domain : interpolation of data without mask
IF (IMSKR == 0) THEN
  WRITE(NULOUT,'('' NO EFFECTIVE MASK - VEGETATION CHARACTERISTICS'')')
  IFLD=1
  CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,4),IDATY,IDATX,IFLD,ZCMP(1,4),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  No available data over the ALADIN domain : default values (desert) used
ELSEIF (IMSKI == ITFING) THEN
  WRITE(NULOUT,'('' VEG. CHARACTERISTICS NOT AVAILABLE !'')')
  IF (LSOLV) THEN
    ZD=YRCLI%SDEPD
  ELSE
    ZD=YRCLI%SZZ0D
  ENDIF
  DO J=1,ITFING
    ZCMP(J,4)= ZD
    ZCMS(J,4)= ZD
  ENDDO

!  Standard case : interpolation of data with a mask
ELSE
  IFLD=1
  CALL EINTER10(YDGEOMETRY%YREGEO,ZFLD(1,1,4),IDATY,IDATX,IFLD,ZCMP(1,4),ZCMS(1,4),&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT,&
   & ZFLD(1,1,3),ZCMP(1,3),ZLS,ZCMS(1,3),YRCLI%SMASK)  
ENDIF

!     5.5 Deallocate temporary space

DEALLOCATE ( ZFLD )
IF (YRCLI%LIEEE) THEN
  DEALLOCATE (ZAUX1)
  DEALLOCATE (ZAUX2)
ENDIF

!     5.7 Move to ARPEGE field (vegetation cover) and modify useful soil depth

DO J=1,ITFING
  IF (ZCMP(J,2) <= YRCLI%SMANQ) THEN
    ZV=ZCMS(J,2)
  ELSE
    ZV=ZCMP(J,2)
  ENDIF
  IF (ZCMP(J,4) <= YRCLI%SMANQ) THEN
    ZD=ZCMS(J,4)
  ELSE
    ZD=ZCMP(J,4)
  ENDIF
  ZV=MAX(0.0_JPRB,MIN(1.0_JPRB,ZV))
  IF (LSOLV) THEN
    IF (ZV < YRCLI%SVEG) ZV= 0.0_JPRB
    ZRES(J,8)=ZV*MIN(ZD,ZRES(J,4))+(1.0_JPRB-ZV)*ZRES(J,4)
  ELSE
    ZRES(J,4)=MAX(0.0_JPRB,ZD)
  ENDIF
  ZRES(J,7)=ZV
ENDDO

!     ------------------------------------------------------------------

!     6. FINAL CALCULATIONS AND WRITING.
!     ----------------------------------

!     6.1 Physical computations on the final grid

!  New fields

IF (LSOLV) THEN
  DO J=1,ITFING
    IV=NINT(ZRES(J,1))

!  Sea or lake
    IF (IV == YRCLI%NTPMER .OR. IV == YRCLI%NTPLAC) THEN
      ZRES(J,2)=YRCLI%SALBM
      ZRES(J,3)=YRCLI%SEMIM
      ZRES(J,4)=YRCLI%SDEPX
      ZRES(J,5)=YRCLI%SARGN
      ZRES(J,6)=YRCLI%SSABN
      ZRES(J,7)=0.0_JPRB
      ZRES(J,8)=YRCLI%SDEPX
      ZRES(J,9)=0.0_JPRB

!  Ice
    ELSEIF (IV == YRCLI%NTPGLA) THEN
      ZRES(J,2)=YRCLI%SALBG
      ZRES(J,3)=YRCLI%SEMIG
      ZRES(J,4)=YRCLI%SDEPN
      ZRES(J,5)=YRCLI%SARGN
      ZRES(J,6)=YRCLI%SSABN
      ZRES(J,7)=0.0_JPRB
      ZRES(J,8)=YRCLI%SDEPN
      ZRES(J,9)=0.0_JPRB

!  Land
    ELSEIF ( (IV > 0) .AND. (IV <= JPTYVE) ) THEN
      ZRES(J,2)=MAX(YRCLI%SALBN,MIN(YRCLI%SALBX,ZRES(J,2)))
      ZRES(J,3)=MAX(YRCLI%SEMIN,MIN(YRCLI%SEMIX,ZRES(J,3)))
      ZRES(J,4)=MAX(YRCLI%SDEPN,MIN(YRCLI%SDEPX,ZRES(J,4)))
      ZZARG=MAX(YRCLI%SARGN,MIN(YRCLI%SARGX,ZRES(J,5)))
      ZZSAB=MAX(YRCLI%SSABN,MIN(YRCLI%SSABX,ZRES(J,6)))
      ZSUM=MAX(100._JPRB,ZZARG+ZZSAB)/100._JPRB
      ZRES(J,5)=ZZARG/ZSUM
      ZRES(J,6)=ZZSAB/ZSUM
      ZRES(J,7)=ZRES(J,7)*ZRES(J,-1)
      ZRES(J,8)=MAX(YRCLI%SDEPN,MIN(ZRES(J,4),ZRES(J,8)))
      ZRES(J,9)=MIN(1.0_JPRB-ZRES(J,7),ZRES(J,9))

    ELSE
      CALL ABOR1(' EINCLI2 : INVALID LAND USE TYPE ! -2-')
    ENDIF
  ENDDO

!  Old fields

ELSE
  ZFZ0=(YRCLI%SFCZ0*RG)**2
  DO J=1,ITFING

!  Sea or lake
    IF ( NINT(ZRES(J,0)) == 0 ) THEN
      IF (ZRES(J,2) <= 0.5_JPRB*(YRCLI%SALBM+YRCLI%SALBB)) THEN
        ZRES(J,2)=YRCLI%SALBM
        ZRES(J,3)=YRCLI%SEMIM
        ZRES(J,5)=YRCLI%SZZ0M*RG
      ELSE
        ZRES(J,2)=YRCLI%SALBB
        ZRES(J,3)=YRCLI%SEMIB
        ZRES(J,5)=YRCLI%SZZ0B*RG
      ENDIF
      ZRES(J,7)=0.0_JPRB

!  Land
    ELSEIF ( NINT(ZRES(J,0)) == 1 ) THEN
      ZRES(J,2)=MAX(YRCLI%SALBN,MIN(YRCLI%SALBX,ZRES(J,2)))
      ZRES(J,3)=MAX(YRCLI%SEMIN,MIN(YRCLI%SEMIX,ZRES(J,3)))
      ZZ0R=ZRES(J,7)*ZRES(J,4)**2 + ZRES(J,6)*YRCLI%SZZ0U**2
      ZRES(J,5)=SQRT(ZRES(J,5)**2 + ZFZ0*ZZ0R)
    ELSE
      CALL ABOR1(' EINCLI2 : INVALID LAND SEA MASK !')
    ENDIF
  ENDDO
ENDIF

!     6.2 Final biperiodicisation and writing on Arpege File

! Albedo (monthly) and bare ground albedo,
! vegetation fraction (monthly) and maximum vegetation fraction,
! are written as equal at this stage.

! New fields
IF (LSOLV) THEN
  IFLD=JPNFIX+4

  DO JC=2,JPNFIX+1
    DO  J=1,ITFING
      ZEXT(J,JC)=ZRES(J,JC)
    ENDDO
  ENDDO
  DO  J=1,ITFING
    ZEXT(J,JPNFIX+2)=ZRES(J,2)
    ZEXT(J,JPNFIX+3)=ZRES(J,7)
    ZEXT(J,JPNFIX+4)=ZRES(J,2)
  ENDDO

  CALL EBICLI(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF,IFLD+1,INIVL(0),CLPREF(0),CLSUFF(0),INUM,ZEXT(1,0),&
   & LLBIP(0),LLWRI(0),LLPAC(0))

! Old fields
ELSE
  IFLD=4

  DO  J=1,ITFING
    ZEXT(J,1)=ZRES(J,2)
    ZEXT(J,2)=ZRES(J,3)
    ZEXT(J,3)=ZRES(J,5)
    ZEXT(J,4)=ZRES(J,7)
  ENDDO
  CLSUFF(1)=CLSUFF(10)
  CLSUFF(2)=CLSUFF( 3)
  CLSUFF(3)=CLSUFF( 5)
  CLSUFF(4)=CLSUFF(11)

  CALL EBICLI(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF,IFLD+1,INIVL(0),CLPREF(0),CLSUFF(0),INUM,ZEXT(1,0),&
   & LLBIP(0),LLWRI(0),LLPAC(0))

ENDIF

CALL LFILAF(IREP,INUM,.TRUE.)
CALL FAIRME(IREP,INUM,'KEEP')

!     ------------------------------------------------------------------


END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EINCLI2',1,ZHOOK_HANDLE)

END SUBROUTINE EINCLI2
