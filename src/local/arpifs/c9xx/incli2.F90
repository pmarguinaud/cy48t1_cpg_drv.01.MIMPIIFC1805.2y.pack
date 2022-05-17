SUBROUTINE INCLI2(YDGEOMETRY,YDPHY)

!**** *INCLI2*

!     PURPOSE
!     -------

!     This routine calculates fixed climatological fields.

!**   INTERFACE
!     ---------

!     CALL INCLI2

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
!     The target grid is a Gaussian grid with rotation of the pole and "Schmidt"
!     stretching as well as -if asked- reduction of the number of grid points
!     towards the pole of the final representation.
!     Input data are global, and coming from units 11 to 18 :
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
!     Data over sea are assumed to be missing (undefined) for soil depth and
!     texture as well as for vegetation characteristics. The case of no missing
!     data (no mask) is considered for bare ground albedo and emissivity.
!     Information from the fraction of vegetation cover is used for root depth.
!     A 4-point interpolation operator is used for all fields. The values on
!     the Gaussian grid are the averages of the values in the boxes (a box is
!     a part of the subgrid corresponding to a point of the Gaussian grid).
!     The fraction of vegetation cover is weighted by the fraction of land.
!     The output fields are added to the ARPEGE file.

!     When old-fashioned fields are required (case LSOLV=.FALSE.), only the
!     following fields are computed or modified :
!     - albedo
!     - emissivity
!     - vegetation cover
!     - g*roughness length
!       (adding the contributions from urbanization and vegetation)

!     EXTERNALS
!     ---------

!     GEO923
!     INTER8
!     INTER10
!     FA-LFI package (FAITOU,FACILE,FAIENC,LFILAF,FAIRME)
!     CHIEN
!     ABOR1

!     AUTHORS
!     -------
!      D. Giard 97-03-13  from INCLIA and INCLIC

!     MODIFICATIONS
!     -------------
!      D. Giard 01-11-14 use of global variables for DM
!                        more allocatable arrays to save memory
!      R. El Khatib 03-08-18 Roughness lengths not packed
!      M.Hamrud 01-10-03 CY28 Cleaning
!      D. Giard 04-09-15 cosmetic changes after automatic ones 
!      K. Yessad (Jan 2010): externalisation of group EGGX in XRD/IFSAUX
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMMP0   , ONLY : NPRINTLEV
USE YOMCT0   , ONLY : NQUAD
USE YOMLUN   , ONLY : NULOUT
USE YOMPHY   , ONLY : TPHY
USE YOMCST   , ONLY : RG
USE YOMVERT  , ONLY : VP00
USE YOMCLI   , ONLY : YRCLI  

!     ------------------------------------------------------------------

IMPLICIT NONE

! JPTYVE : number of land-use types, at least 3 :
!   1 -> sea , 2 -> ice , 3 -> desert or low vegetation (lakes are 1 or 5)
! JPNFIX : number of fixed fields to be derived
! JPNF   : maximum number of fields to be interpolated together
! JPBX/JPBY : number of extra longitudes/latitudes on each side of the data
!          grid required by interpolation (4 or 12 points -> 2/2)

TYPE(GEOMETRY), INTENT(IN)   :: YDGEOMETRY
TYPE(TPHY)     ,INTENT(INOUT):: YDPHY
INTEGER(KIND=JPIM) :: JPBX
INTEGER(KIND=JPIM) :: JPBY
INTEGER(KIND=JPIM) :: JPNF
INTEGER(KIND=JPIM) :: JPNFIX
INTEGER(KIND=JPIM) :: JPTYVE

PARAMETER(JPTYVE=5,JPNFIX=8,JPNF=MAX(JPTYVE,4),JPBX=2,JPBY=2)

REAL(KIND=JPRB),ALLOCATABLE :: ZFLD (:,:,:),ZITP(:)
REAL(KIND=JPRB),ALLOCATABLE :: ZAUX1(:),ZAUX2(:),ZAUX3(:)
REAL(KIND=JPRB),ALLOCATABLE :: ZIVD (:),ZIVS(:) 
REAL(KIND=JPRB) :: ZRES(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON,-1:JPNFIX+1)
REAL(KIND=JPRB) :: ZLSR(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON)
REAL(KIND=JPRB) :: ZCMP(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON,JPNF),&
 & ZCMS(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON,JPNF)&
 & ,ZBO(YDGEOMETRY%YRDIM%NDGLG+1),ZLONG(YDGEOMETRY%YRDIM%NDGLG),ZMU(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON),&
 & ZSLA(YDGEOMETRY%YRDIM%NDGLG*YRCLI%NPINT)&
 & ,ZSLO(YDGEOMETRY%YRDIM%NDLON*YRCLI%NPINT,YDGEOMETRY%YRDIM%NDGLG),ZCLO(YDGEOMETRY%YRDIM%NDLON*YRCLI%NPINT,YDGEOMETRY%YRDIM%NDGLG)  
REAL(KIND=JPRB) :: ZD, ZEPS, ZFZ0, ZMAX1, ZMAX2, ZSUM, ZV, ZZ0R, ZZARG, ZZSAB
REAL(KIND=JPRB) :: ZHOOK_HANDLE

CHARACTER :: CLPREF(-1:JPNFIX+4)*8,CLSUFF(-1:JPNFIX+4)*12
CHARACTER :: CLNOMC*16,CLNOMF*10,CLFORM*12

INTEGER(KIND=JPIM) :: IARI, IARP, ICO, IDATX, IDATY, IFLD, IINF,&
 & IK, ILENE, ILENT, IMES, IMSK, INIQ, INIV, INJQ, INUM, IREP,&
 & ILSM, ITP, IV, IVD, IVN, IVS, IX, IY, IZ, J, JJ, JT  
INTEGER(KIND=JPIM) :: INGRIG,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL

LOGICAL :: LLCOSP, LLIMST, LLPOLE

!     ------------------------------------------------------------------

#include "chien.h"

#include "abor1.intfb.h"
#include "geo923.intfb.h"
#include "inter10.intfb.h"
#include "inter8.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INCLI2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDVAB=>YDGEOMETRY%YRVAB&
& )
ASSOCIATE(NDGENG=>YDDIM%NDGENG, NDGLG=>YDDIM%NDGLG, NDGSAG=>YDDIM%NDGSAG,   NDLON=>YDDIM%NDLON, NSMAX=>YDDIM%NSMAX,   &
& NFLEVG=>YDDIMV%NFLEVG,   NHTYP=>YDGEM%NHTYP, NLOENG=>YDGEM%NLOENG, NMENG=>YDGEM%NMENG,   NSTTYP=>YDGEM%NSTTYP,      &
& RLOCEN=>YDGEM%RLOCEN, RMUCEN=>YDGEM%RMUCEN,   RSTRET=>YDGEM%RSTRET, LSOLV=>YDPHY%LSOLV)
!     ------------------------------------------------------------------

!     1. SET INITIAL VALUES.
!        -------------------

!     1.1 Dimensions for data and interpolation.

ICO=YRCLI%NPINT*YRCLI%NPINT
INJQ=NDGLG*YRCLI%NPINT
INIQ=NDLON*YRCLI%NPINT
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

!     1.5 Open ARPEGE File and read land-sea mask and land fraction

INUM=3
IREP=0
IARI=0
IARP=10
IMES=1
CLNOMF='Const.Clim'
CLNOMC='Const.Clim.Surfa'
LLIMST=.TRUE.

IINF=-1
ZEPS=1.E-10_JPRB
LLPOLE=.TRUE.

INIV=0
LLCOSP=.FALSE.

CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'OLD',.TRUE.,LLIMST,&
 & IMES,IARP,IARI,CLNOMC)  
CALL CHIEN(CLNOMC,NSTTYP,RMUCEN,RLOCEN,RSTRET,NSMAX,&
 & NDGLG,NDLON,NLOENG,NMENG,NHTYP,NFLEVG,VP00,YDVAB%VALH,YDVAB%VBH,&
 & NQUAD,IINF,NDGSAG,NDGENG,ZEPS,LLPOLE,NULOUT)  
IF (LLPOLE) CALL ABOR1(' CLIM. FILES MUST NOT HAVE POLES !')

DO J=-1,0
  CALL FACILE(IREP,INUM,CLPREF(J),INIV,CLSUFF(J),ZRES(1,J),LLCOSP)
ENDDO
DO J=1,ILENE
  ZLSR(J)=1.0_JPRB-ZRES(J,0)
ENDDO
IF (.NOT. LSOLV) THEN
  CALL FACILE(IREP,INUM,CLPREF(5),INIV,CLSUFF(5),ZRES(1,5),LLCOSP)
  CALL FACILE(IREP,INUM,CLPREF(6),INIV,CLSUFF(6),ZRES(1,6),LLCOSP)
ELSE
  CALL FACILE(IREP,INUM,CLPREF(9),INIV,CLSUFF(9),ZRES(1,9),LLCOSP)
ENDIF

!     ------------------------------------------------------------------

!     2. DOMINANT LAND USE TYPE.
!        -----------------------

IF (.NOT. LSOLV) GO TO 300

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
  DO JJ=1,YRCLI%NDATY
    IY=JJ+JPBY
    DO J=1,YRCLI%NDATX
      IX=J+JPBX
      IZ=J+(JJ-1)*YRCLI%NDATX
      ITP=NINT(ZITP(IZ))
      IF ((ITP <= 0) .OR. (ITP > JPTYVE))&
       & CALL ABOR1(' INCLI2 : INVALID LAND-USE TYPE ! -1-')  
      ZFLD(IX,IY,ITP)=1.0_JPRB
    ENDDO
  ENDDO
ELSE
  DO JJ=1,YRCLI%NDATY
    IY=JJ+JPBY
    READ(11,*) (ZITP(J),J=1,YRCLI%NDATX)
    DO J=1,YRCLI%NDATX
      IX=J+JPBX
      ITP=NINT(ZITP(J))
      IF ((ITP <= 0) .OR. (ITP > JPTYVE))&
       & CALL ABOR1(' INCLI2 : INVALID LAND-USE TYPE ! -1-')  
      ZFLD(IX,IY,ITP)=1.0_JPRB
    ENDDO
  ENDDO
ENDIF
CLOSE(11)

!     2.3 Interpolation -> proportion of each type

IFLD=JPTYVE
CALL INTER8(NDGLG,NDLON,IDATY,IDATX,IFLD,YRCLI%NPINT,ILENE,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZCMP(1,1),ZFLD(1,1,1),&
 & ZSLA,ZSLO,ZCLO)  

!     2.4 Deallocate/allocate temporary space

DEALLOCATE ( ZITP )
DEALLOCATE ( ZFLD )

ALLOCATE ( ZIVD (ILENT) )
ALLOCATE ( ZIVS (ILENT) )

!     2.5 Determination of dominant and secondary types

DO J=1,ILENE
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

DO J=1,ILENE
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
    IX=1
    IY=J+1
    IK=MAX(1,MIN(ILENE,IY))
    IVN=NINT(ZIVD(IK))
262 CONTINUE
    IF ((IVN /= YRCLI%NTPMER).AND.(IVN /= YRCLI%NTPLAC).AND.(IVN /= 0)) GO TO 263
    IX=SIGN(IABS(IX)+1,-IX)
    IY=IY+IX
    IK=MAX(1,MIN(ILENE,IY))
    IVN=NINT(ZIVD(IK))
    GOTO 262
263 CONTINUE
    IVD=IVN
    IF (NPRINTLEV >= 1) WRITE(NULOUT,FMT='(&
     & '' DOMINANT TYPE OF VEGET. FORCED TO'',I2,'', J='',I8)') IVD,J  
  ENDIF
  ZIVD(J)=REAL(IVD,JPRB)
ENDDO

!     2.7 ARPEGE field

DO J=1,ILENE
  ZRES(J,1)=ZIVD(J)
ENDDO

!     2.8 Deallocate temporary space

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
  DO JJ=1+JPBY,YRCLI%NDATY+JPBY
    DO J=1,YRCLI%NDATX
      IX=J+JPBX
      IZ=J+(JJ-1-JPBY)*YRCLI%NDATX
      ZFLD(IX,JJ,2)=ZAUX1(IZ)
      ZFLD(IX,JJ,3)=ZAUX2(IZ)
    ENDDO
  ENDDO
ELSE
  READ(12,*)((ZFLD(J,JJ,2),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=1+JPBY,YRCLI%NDATY+JPBY)
  READ(13,*)((ZFLD(J,JJ,3),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=1+JPBY,YRCLI%NDATY+JPBY)
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

IF (IMSK == 0) THEN
  WRITE(NULOUT,'('' NO MASK FOR BARE GROUND ALBEDO AND EMISSIVITY !'')')
  IFLD=2
  CALL INTER8(NDGLG,NDLON,IDATY,IDATX,IFLD,YRCLI%NPINT,ILENE,&
   & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZCMP(1,2),ZFLD(1,1,2),&
   & ZSLA,ZSLO,ZCLO)  
ELSEIF (IMSK == YRCLI%NDATX*YRCLI%NDATY) THEN
  CALL ABOR1(' INCLI2 : ALBEDO OR EMISSIVITY MISSING !')
ELSE
  IFLD=1
  CALL INTER8(NDGLG,NDLON,IDATY,IDATX,IFLD,YRCLI%NPINT,ILENE,&
   & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZCMP(1,1),ZFLD(1,1,1),&
   & ZSLA,ZSLO,ZCLO)  
  IFLD=2
  CALL INTER10(NDGLG,NDLON,IDATY,IDATX,IFLD,YRCLI%NPINT,ILENE,&
   & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZCMP(1,2),ZCMS(1,2),&
   & ZFLD(1,1,2),ZSLA,ZSLO,ZCLO,&
   & ZFLD(1,1,1),ZCMP(1,1),ZLSR(1),ZCMS(1,1),YRCLI%SMASK)  
ENDIF

!     3.5 Deallocate temporary space

DEALLOCATE ( ZFLD )
IF (YRCLI%LIEEE) THEN
  DEALLOCATE (ZAUX1)
  DEALLOCATE (ZAUX2)
ENDIF

!     3.6 Move to ARPEGE fields

DO J=1,ILENE
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
  DO JJ=1+JPBY,YRCLI%NDATY+JPBY
    DO J=1,YRCLI%NDATX
      IX=J+JPBX
      IZ=J+(JJ-1-JPBY)*YRCLI%NDATX
      ZFLD(IX,JJ,2)=ZAUX1(IZ)
      ZFLD(IX,JJ,3)=ZAUX2(IZ)
      ZFLD(IX,JJ,4)=ZAUX3(IZ)
    ENDDO
  ENDDO
ELSE
  READ(14,*)((ZFLD(J,JJ,2),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=1+JPBY,YRCLI%NDATY+JPBY)
  READ(15,*)((ZFLD(J,JJ,3),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=1+JPBY,YRCLI%NDATY+JPBY)
  READ(16,*)((ZFLD(J,JJ,4),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=1+JPBY,YRCLI%NDATY+JPBY)
ENDIF
CLOSE(14)
CLOSE(15)
CLOSE(16)

!     4.3 Define the mask (missing data)

DO JJ=JPBY+1,JPBY+YRCLI%NDATY
  DO J=JPBX+1,JPBX+YRCLI%NDATX
    IF (MIN(ZFLD(J,JJ,2),ZFLD(J,JJ,3),ZFLD(J,JJ,4)) <= YRCLI%SMANQ) THEN
      ZFLD(J,JJ,1)=1.0_JPRB
    ELSE
      ZFLD(J,JJ,1)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

!     4.4 Interpolate mask and fields
IFLD=1
CALL INTER8(NDGLG,NDLON,IDATY,IDATX,IFLD,YRCLI%NPINT,ILENE,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZCMP(1,1),ZFLD(1,1,1),&
 & ZSLA,ZSLO,ZCLO)  

IFLD=3
CALL INTER10(NDGLG,NDLON,IDATY,IDATX,IFLD,YRCLI%NPINT,ILENE,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZCMP(1,2),ZCMS(1,2),&
 & ZFLD(1,1,2),ZSLA,ZSLO,ZCLO,&
 & ZFLD(1,1,1),ZCMP(1,1),ZLSR(1),ZCMS(1,1),YRCLI%SMASK)  

!     4.5 Deallocate temporary space

DEALLOCATE ( ZFLD )
IF (YRCLI%LIEEE) THEN
  DEALLOCATE (ZAUX1)
  DEALLOCATE (ZAUX2)
  DEALLOCATE (ZAUX3)
ENDIF

!     4.6 Move to ARPEGE fields

DO J=1,ILENE
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
  DO JJ=1+JPBY,YRCLI%NDATY+JPBY
    DO J=1,YRCLI%NDATX
      IX=J+JPBX
      IZ=J+(JJ-1-JPBY)*YRCLI%NDATX
      ZFLD(IX,JJ,2)=ZAUX1(IZ)
      ZFLD(IX,JJ,4)=ZAUX2(IZ)
    ENDDO
  ENDDO
ELSE
  READ(17,*)((ZFLD(J,JJ,2),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=1+JPBY,YRCLI%NDATY+JPBY)
  READ(18,*)((ZFLD(J,JJ,4),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=1+JPBY,YRCLI%NDATY+JPBY)
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
IFLD=1
CALL INTER8(NDGLG,NDLON,IDATY,IDATX,IFLD,YRCLI%NPINT,ILENE,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZCMP(1,1),ZFLD(1,1,1),&
 & ZSLA,ZSLO,ZCLO)  

IFLD=1
CALL INTER10(NDGLG,NDLON,IDATY,IDATX,IFLD,YRCLI%NPINT,ILENE,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZCMP(1,2),ZCMS(1,2),&
 & ZFLD(1,1,2),ZSLA,ZSLO,ZCLO,&
 & ZFLD(1,1,1),ZCMP(1,1),ZLSR(1),ZCMS(1,1),YRCLI%SMASK)  

IFLD=1
CALL INTER8(NDGLG,NDLON,IDATY,IDATX,IFLD,YRCLI%NPINT,ILENE,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZCMP(1,3),ZFLD(1,1,3),&
 & ZSLA,ZSLO,ZCLO)  

IFLD=1
CALL INTER10(NDGLG,NDLON,IDATY,IDATX,IFLD,YRCLI%NPINT,ILENE,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZCMP(1,4),ZCMS(1,4),&
 & ZFLD(1,1,4),ZSLA,ZSLO,ZCLO,&
 & ZFLD(1,1,3),ZCMP(1,3),ZLSR(1),ZCMS(1,3),YRCLI%SMASK)  

!     5.5 Deallocate temporary space

DEALLOCATE ( ZFLD )
IF (YRCLI%LIEEE) THEN
  DEALLOCATE (ZAUX1)
  DEALLOCATE (ZAUX2)
ENDIF

!     5.6 Move to ARPEGE field (vegetation cover) and modify useful soil depth

DO J=1,ILENE
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

!     6.1 Physical computations on the reduced Gaussian grid

!  New fields

IF (LSOLV) THEN
  DO J=1,ILENE
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
      CALL ABOR1(' INCLI2 : INVALID LAND USE TYPE ! -2-')
    ENDIF
  ENDDO

!  Old fields

ELSE
  ZFZ0=(YRCLI%SFCZ0*RG)**2
  DO J=1,ILENE
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
      CALL ABOR1(' INCLI2 : INVALID LAND USE TYPE ! -2-')
    ENDIF
  ENDDO
ENDIF

!     6.2 Final writing on Arpege File

! Albedo (monthly) and bare ground albedo,
! vegetation fraction (monthly) and maximum vegetation fraction,
! are written as equal at this stage.

!  New fields
IF (LSOLV) THEN
  DO J=1,JPNFIX+1
    CALL FAIENC(IREP,INUM,CLPREF(J),INIV,CLSUFF(J),ZRES(1,J),LLCOSP)
  ENDDO
  CALL FAIENC(IREP,INUM,CLPREF(10),INIV,CLSUFF(10),ZRES(1,2),LLCOSP)
  CALL FAIENC(IREP,INUM,CLPREF(11),INIV,CLSUFF(11),ZRES(1,7),LLCOSP)
  CALL FAIENC(IREP,INUM,CLPREF(12),INIV,CLSUFF(12),ZRES(1,2),LLCOSP)
!  Old fields
ELSE
  CALL FAIENC(IREP,INUM,CLPREF(10),INIV,CLSUFF(10),ZRES(1,2),LLCOSP)
  CALL FAIENC(IREP,INUM,CLPREF( 3),INIV,CLSUFF( 3),ZRES(1,3),LLCOSP)
! Do not pack roughness lengths
  CALL FAVEUR(IREP,INUM,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  INGRIG=0
  CALL FAGOTE(IREP,INUM,INGRIG,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  CALL FAIENC(IREP,INUM,CLPREF( 5),INIV,CLSUFF( 5),ZRES(1,5),LLCOSP)
! Reset packing after the treatment of roughness lengths
  CALL FAGOTE(IREP,INUM,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  CALL FAIENC(IREP,INUM,CLPREF(11),INIV,CLSUFF(11),ZRES(1,7),LLCOSP)
ENDIF

CALL LFILAF(IREP,INUM,.TRUE.)
CALL FAIRME(IREP,INUM,'KEEP')

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INCLI2',1,ZHOOK_HANDLE)
END SUBROUTINE INCLI2
