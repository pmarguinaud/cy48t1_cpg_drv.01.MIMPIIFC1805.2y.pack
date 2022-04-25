SUBROUTINE EINCLI1(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF)

!**** *EINCLI1*

!     PURPOSE.
!     --------

!     This routine calculates first 7 fixed fields :
!      .land(1)/water(0) mask,
!      .surface geopotential (grid point and spectral)
!        -including enveloppe effect, if asked-,
!      .G times the standard deviation of the orography,
!      .the anisotropy coefficient
!      .the direction of principal axis of topography (in radian),
!      .G times the roughness length of bare surface
!        -with the possibility of scaling it by an arbitrary parameter-,
!      .the fraction of land.
!      .the fraction of urbanization.
!       The target grid is a regular grid(geographical or in plan-projection)
!     with rotation of the pole(if asked)
!     Input data are coming from the *MANU* type files, based on GLOBE data
!     ALL input fields are organized from 0 to 360 deg in longitude
!                                    from N to S in latitude

!**   INTERFACE.
!     ----------

!     CALL EINCLI1

!     METHOD.
!     -------

!     The *MANU* dataset is read on units 31-39.
!     A subgrid of the final grid is introduced, and a value on this subgrid
!     is the value of the nearest point on the initial grid.
!     The values on the final grid are the averages of the values in the
!     boxes (a box is the part of the subgrid corresponding to a point of the
!     final grid).
!     The 7 "fixed" output fields are written on unit n=41 (one output file).

!     EXTERNALS.
!     ----------
!      ABOR1
!      EGEO923
!      EINTER1
!      ELISLAP
!      EBICLI (previously ETIBIHI and EBIEN)
!      GTOPTX2
!      GTOPTXY
!      GTOPTY2
!      LOCMAXI (LOCMAXP)
!      RELNEW
!      RELSPE
!      SPREORD
!      FA-LFI package (FACADE,FAITOU,FANDAR,FAVEUR,FAGOTE,FAIENC,LFILAF,FAIRME)

!     AUTHORS.
!     --------
!      L. GERARD 97-05-26 from INCLI1 and E. Cordoneanu's EINCLIA

!     MODIFICATIONS.
!     --------------
!      M. Szczech 01-03-14 : Phasing ->AL15 (no more lmessp)
!      R. El Khatib : 01-12-06 Cleaning sm variables
!      R. El Khatib : 01-08-07 Pruning options
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R.ZAABOUL : 08-Dec-2003 new arguments for EBICLI
!      D. Paradis & R. El Khatib : 04-07-22 GRIBEX
!      J. Woyciechowska 03-12 all spectral transforms moved to ERELSPE
!                             initializations moved to EINCLI1
!      D. Giard 04-01-07 : diagnostics moved to EINCLI1, simpler code
!      D. Giard 05-04-04 : update and cleaning, semi-envelope suppressed
!      D. Giard 05-04-07 : no packing for roughness lengths, new EBICLI
!      JD. Gril 05-05-12 : add LMRT=true => NROTEQ/ZSINLA(1) = -2
!      JD. Gril : 01-02-2006 Test for LMRT=.T. and NCADFORM=0 => abort
!      A. Trojakova : 30/05/2006 new argument for ELISLAP 
!                     due to Bf for map factor addressing on multiproc
!      S. Ivatek-Sahdan : 06-07-06 map factor for ELISLAP just for (C+I)
!      K.Essaouini  Jan. 2010 : Add envelop if orography read from a PGD file (LIPGD)
!      A. Bogatchev 12-04-2013 YRVAB => YOMVERT
!      P.Marguinaud 10-10-2014 Use FAIENO
!      B. Bochenek (Apr 2015): Phasing: move some variables.
!      O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!      F. Taillefer : March 2017 bf (replace call faieno by faienc)
!----------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE YOEPHY       , ONLY : TEPHY
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPIS
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : RPI      ,RA       ,RG
USE YOMCLA   , ONLY : NLISSZ   ,NLISSR   ,FACZ0    ,FENVN    ,&
 & LNORO    ,LNLSM    ,LNEWORO  ,LNEWORO2 ,LKEYF    ,LIPGD   ,&
 & LSPSMORO
USE YOMCLI   , ONLY : YRCLI
USE YEMCLI   , ONLY : NJLCLI1  ,NJLCLI2
USE YOMVERT  , ONLY : VP00 
USE YOMLUN   , ONLY : NULOUT
USE PTRSPOR  , ONLY : XRNM     ,XREV
USE YOMOPH0  , ONLY : NCADFORM
USE YOM_YGFL , ONLY : TYPE_GFLD

!----------------------------------------------------------------------

IMPLICIT NONE

! Input geometry:
TYPE(GEOMETRY), INTENT(INOUT)   :: YDGEOMETRY
TYPE(TYPE_GFLD)    ,INTENT(INOUT):: YDGFL
TYPE(TEPHY)    ,INTENT(INOUT)   :: YDEPHY
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
REAL(KIND=JPRB),ALLOCATABLE :: ZVA (:,:)
INTEGER(KIND=JPIS) :: IARRAY16(YRCLI%NDATX)
INTEGER(KIND=JPIM) :: IARRAY32(YRCLI%NDATX)
CHARACTER (LEN = 1) ::  CL_CARRAY(YRCLI%NDATX)
! Output geometry (C+I):
REAL(KIND=JPRB) :: ZS((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1),14)
REAL(KIND=JPRB) :: ZGM((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1))
! Output array (C+I+E)
REAL(KIND=JPRB) :: ZEXT(YDGEOMETRY%YRDIM%NDLON*YDGEOMETRY%YRDIM%NDGLG,8),ZWR(YDGEOMETRY%YRDIM%NDLON*YDGEOMETRY%YRDIM%NDGLG,1),&
 & ZZWI(YDGEOMETRY%YRDIM%NDLON*YDGEOMETRY%YRDIM%NDGLG,1)
! Cadre and Writing Setup
REAL(KIND=JPRB) :: ZSINLA(YDGEOMETRY%YRDIM%NDGNH+1)
INTEGER(KIND=JPIM) :: INIP(YDGEOMETRY%YRDIM%NDGLG+2),INOZPA(YDGEOMETRY%YRDIM%NDGNH+1),IDATEF(11),&
 & INIVL(8),INIV(1)
INTEGER(KIND=JPIM) :: IADL(YDGEOMETRY%YRDIM%NDGLG)
CHARACTER :: CLNOMC*16,CLNOMF*10,CLPREF(0:9)*8,CLSUFF(0:9)*12
LOGICAL :: LLBIP(8),LLWRI(8),LLPAC(8),LLBI(1),LLWR(1),LLPA(1)

INTEGER(KIND=JPIM) :: I8, IBCSP, IBPDG, IC, ICO, IADR, IDATX, IDATY,&
 & IDMO, IFLD, IGRIB, IGRIG, IJJ1,&
 & IJJ2, ILAP, IMES, INBARI, INBARP,&
 & INUM, IREP, ITFING, ITN, ITRON, ITRONC,&
 & ITYPTR, IXFING, IYFING, J, JA, JC, JI, JJ,&
 & JLOOP, JS, JX, JY, ITR,&
 & ISTARTCIE, ISTARTCI, IENDCIE, IENDCI

LOGICAL :: LLCOSP, LLERFA, LLFIT, LLIMST, LLKEEP, LLNOMM, LLOPT

REAL(KIND=JPRB) :: Z1, Z2, ZCM, ZCODIL, ZCV, ZEPS, ZPOR, ZPOZ, ZSURF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!----------------------------------------------------------------------

#include "abor1.intfb.h"
#include "ebicli.intfb.h"
#include "eganiso.intfb.h"
#include "egeo923.intfb.h"
#include "einter1.intfb.h"
#include "elislap.intfb.h"
#include "gtoptx2.intfb.h"
#include "gtoptxy.intfb.h"
#include "gtopty2.intfb.h"
#include "inipz.intfb.h"
#include "inirp.intfb.h"
#include "locmaxi.intfb.h"
#include "relnew.intfb.h"
#include "relspe.intfb.h"
#include "faieno.h"

!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EINCLI1',0,ZHOOK_HANDLE)

ASSOCIATE( YDVAB=>YDGEOMETRY%YRVAB, &
 & NBZONG=>YDGEOMETRY%YREDIM%NBZONG, NBZONL=>YDGEOMETRY%YREDIM%NBZONL, &
 & NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & NDGLG=>YDGEOMETRY%YRDIM%NDGLG, NDGNH=>YDGEOMETRY%YRDIM%NDGNH, NDGUNG=>YDGEOMETRY%YRDIM%NDGUNG, &
 & NDGUXG=>YDGEOMETRY%YRDIM%NDGUXG, NDLON=>YDGEOMETRY%YRDIM%NDLON, NDLUNG=>YDGEOMETRY%YRDIM%NDLUNG, &
 & NDLUXG=>YDGEOMETRY%YRDIM%NDLUXG, NMSMAX=>YDGEOMETRY%YRDIM%NMSMAX, NSMAX=>YDGEOMETRY%YRDIM%NSMAX, &
 & NSEFRE=>YDGEOMETRY%YRDIM%NSEFRE, NSPEC2G=>YDGEOMETRY%YRDIM%NSPEC2G, &
 & NGPTOT=>YDGEOMETRY%YRGEM%NGPTOT, NTSTAGP=>YDGEOMETRY%YRGEM%NTSTAGP, &
 & ERPK=>YDGEOMETRY%YREGEO%ERPK, ELONC=>YDGEOMETRY%YREGEO%ELONC, ELATC=>YDGEOMETRY%YREGEO%ELATC, &
 & ELON1=>YDGEOMETRY%YREGEO%ELON1, ELAT1=>YDGEOMETRY%YREGEO%ELAT1, ELON2=>YDGEOMETRY%YREGEO%ELON2, &
 & ELAT2=>YDGEOMETRY%YREGEO%ELAT2, ELON0=>YDGEOMETRY%YREGEO%ELON0, ELAT0=>YDGEOMETRY%YREGEO%ELAT0, &
 & EDELX=>YDGEOMETRY%YREGEO%EDELX, EDELY=>YDGEOMETRY%YREGEO%EDELY, ELX=>YDGEOMETRY%YREGEO%ELX, ELY=>YDGEOMETRY%YREGEO%ELY, &
 & EXWN=>YDGEOMETRY%YREGEO%EXWN, EYWN=>YDGEOMETRY%YREGEO%EYWN, LMAP=>YDGEOMETRY%YREGEO%LMAP, LMRT=>YDGEOMETRY%YREGEO%LMRT )
!----------------------------------------------------------------------

!     1. SET INITIAL VALUES
!        ------------------

!     1.1 Constants, initial setup of arrays.

ZEPS=1.E-06_JPRB

IADL(1:NDGLG)=NTSTAGP(1:NDGLG)-1

! Characteristic lengths for smoothing (z0, h)
ZPOR=4000._JPRB
ZPOZ=4000._JPRB
 
! Orography set to 0 at first
DO J=1,NGPTOT
  XREV(J)= 0.0_JPRB
ENDDO
DO JS=1,NSEFRE
  XRNM(JS)= 0.0_JPRB
ENDDO

! Same default
DO JC=1,8
  DO J=1,NDGLG*NDLON
    ZEXT(J,JC)= 0.0_JPRB
  ENDDO
ENDDO
DO J=1,NDGLG*NDLON
  ZWR(J,1)= 0.0_JPRB
ENDDO

!     Smoothing is allways done on just one CPU
!     reordering of in this case GM(C+I+E) to ZGM(C+I)
DO JY=1,NDGUXG-NDGUNG+1
  ISTARTCIE=(JY+NDGUNG-2)*NDLON+NDLUNG
  IENDCIE=(JY+NDGUNG-2)*NDLON+NDLUXG
  ISTARTCI=(JY-1)*(NDLUXG-NDLUNG+1)+1
  IENDCI=JY*(NDLUXG-NDLUNG+1)
  ZGM(ISTARTCI:IENDCI)=YDGEOMETRY%YRGSGEOM_NB%GM(ISTARTCIE:IENDCIE)
ENDDO

!     1.2 Dimensions for data interpolation

! Original grid: NDATX x NDATY
IDATX=YRCLI%NDATX
IDATY=YRCLI%NDATY
ITN=YRCLI%NDATX*YRCLI%NDATY
! Final grid:
IXFING=NDLUXG-NDLUNG+1
IYFING=NDGUXG-NDGUNG+1
ITFING=IXFING*IYFING
! number of fields to interpolate:
IFLD=1
! NPINT: size of the interpolation boxes
ICO=YRCLI%NPINT*YRCLI%NPINT

!     1.3 Geometry : initialize YEMCLI and YOMDIL

CALL EGEO923(YDGEOMETRY%YRGEM)

!     1.4 Fields description for biperiodization and writing

INUM=41
CLNOMF='Const.Clim'

IJJ1=NJLCLI1
IJJ2=NJLCLI2

CLPREF( 0)='SPECSURF'
CLPREF( 1)='SURF    '
CLPREF( 2)='SURF    '
CLPREF( 3)='SURF    '
CLPREF( 4)='SURF    '
CLPREF( 5)='SURF    '
CLPREF( 6)='SURF    '
CLPREF( 7)='SURF    '
CLPREF( 8)='SURF    '
CLPREF( 9)='SURF    '
CLSUFF( 0)='GEOPOTEN    '
CLSUFF( 1)='IND.TERREMER'
CLSUFF( 2)='ET.GEOPOTENT'
CLSUFF( 3)='VAR.GEOP.ANI'
CLSUFF( 4)='VAR.GEOP.DIR'
CLSUFF( 5)='Z0REL.FOIS.G'
CLSUFF( 6)='PROP.TERRE  '
CLSUFF( 7)='PROP.URBANIS'
CLSUFF( 8)='Z0.FOIS.G   '
CLSUFF( 9)='GEOPOTENTIEL'

DO J=1,8
  LLBIP(J)=.TRUE.
  LLWRI(J)=.TRUE.
  LLPAC(J)=.TRUE.
  INIVL(J)=0
ENDDO
LLPAC(5)=.FALSE.
LLPAC(8)=.FALSE.

LLBI(1)=.TRUE.
LLWR(1)=.TRUE.
LLPA(1)=.FALSE.
INIV(1)=0

!     ------------------------------------------------------------------
!*
!     2. READING DATA.
!        -------------

!     2.1 Read *MANU* tape.

OPEN(31, FILE='Water_Percentage',&
 & FORM='UNFORMATTED', ACCESS='DIRECT', RECL=YRCLI%NDATX)  

OPEN(32, FILE='Oro_Mean',&
 & FORM='UNFORMATTED', ACCESS='DIRECT', RECL=2*YRCLI%NDATX)  

OPEN(33, FILE='Sigma',&
 & FORM='UNFORMATTED', ACCESS='DIRECT', RECL=2*YRCLI%NDATX)  

OPEN(34, FILE='Nb_Peaks',&
 & FORM='UNFORMATTED', ACCESS='DIRECT', RECL=YRCLI%NDATX)  

OPEN(35, FILE='Urbanisation',&
 & FORM='UNFORMATTED', ACCESS='DIRECT', RECL=YRCLI%NDATX)  

OPEN(36, FILE='Dh_over_Dx_Dh_over_Dy',&
 & FORM='UNFORMATTED', ACCESS='DIRECT', RECL=4*YRCLI%NDATX)  

OPEN(37, FILE='Dh_over_Dx_square',&
 & FORM='UNFORMATTED', ACCESS='DIRECT', RECL=4*YRCLI%NDATX)  

OPEN(38, FILE='Dh_over_Dy_square',&
 & FORM='UNFORMATTED', ACCESS='DIRECT', RECL=4*YRCLI%NDATX)  

OPEN(39, FILE='Hmax-HxH-Hmin_ov4',&
 & FORM='UNFORMATTED', ACCESS='DIRECT', RECL=4*YRCLI%NDATX)  

!     2.2 Allocate temporary space.

ALLOCATE ( ZVA (IDATX,IDATY) )

!     ------------------------------------------------------------------

!     3. DATA READING and INTERPOLATION IN MODEL GRID POINTS
!        ---------------------------------------------------

! REMARK: When reading *MANU* fields, we do a mirror transformation along Y

!    ORIGINAL ATTRIBUTIONS of LOCAL ARRAY ZS:
!    ----------------------------------------
!    ZS1: %w water covered fraction
!    ZS2: mo: mean orography
!    ZS3: mo2: mean squared orography
!    ZS4: sigma**2: small scale (0.5 mn <resol< 2.5 mn) orography variance
!    ZS5: Z0a = N/s: resolved peak density
!    ZS6: z01r small scale (0.5 mn <resol< 2.5 mn) peaks roughness
!    ZS7: v1: tiny scale (resol< 0.5 mn) orography variance
!    ZS8: urb: urbanisation fraction
! Resolved orography covariance tensor:
!    ZS9:  V11: TGCT-11a
!    ZS10: V22: TGCT-22a
!    ZS11: V12: TGCT-12a
! Unresolved orography covariance tensor:
!    ZS12: nu11: TGCT-11b
!    ZS13: nu22: TGCT-22b
!    ZS14: nu12: TGCT-12b

!   1 fraction of water covered area: ZS1=%w

DO JJ=YRCLI%NDATY,1,-1
  JA=YRCLI%NDATY-JJ+1
  READ(31,REC=JA)CL_CARRAY
  DO JI=1,YRCLI%NDATX
    I8=ICHAR(CL_CARRAY(JI))
    ZVA(JI,JJ)=REAL(I8,JPRB)/100._JPRB
  ENDDO
ENDDO
CALL EINTER1(YDGEOMETRY%YREGEO,ZVA,IDATY,IDATX,ITN,IFLD,ZS(1,1),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!   2 mean orography: ZS2=mo

DO JJ=YRCLI%NDATY,1,-1
  JA=YRCLI%NDATY-JJ+1
  READ(32,REC=JA)IARRAY16
  DO JI=1,YRCLI%NDATX
    ZVA(JI,JJ)=REAL(IARRAY16(JI),JPRB)
  ENDDO
ENDDO

CALL EINTER1(YDGEOMETRY%YREGEO,ZVA,IDATY,IDATX,ITN,IFLD,ZS(1,2),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

! Now compute derived fields from mean orography:

!   3 square mean orography: ZS3=mo2

DO JJ=1,YRCLI%NDATY
  DO JI=1,YRCLI%NDATX
    ZVA(JI,JJ)=ZVA(JI,JJ)*ZVA(JI,JJ)
  ENDDO
ENDDO

CALL EINTER1(YDGEOMETRY%YREGEO,ZVA,IDATY,IDATX,ITN,IFLD,ZS(1,3),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!   5 orographical roughness length (part one): ZS5=Z0a
!      N/s=original grid resolved peaks density in finer grid mesh
!  Search for local maxima: ZVA becomes 1 at local max, 0 elsewhere.

CALL LOCMAXI(IDATY,IDATX,ZVA)
!     CALL LOCMAXP(IDATY,IDATX,ZVA)

ZCV=RPI/REAL(YRCLI%NGLOBY,JPRB)
ZCM=2*RPI*RA**2/REAL(YRCLI%NGLOBX,JPRB)
DO JJ=IJJ1,IJJ2
! ZSURF is the surface of a mesh in the original (long-lat) grid:
  ZSURF=ZCM*(COS(ZCV*(JJ-1))-COS(ZCV*JJ))
  DO JI=1,YRCLI%NDATX
    ZVA(JI,JJ-IJJ1+1)=ZVA(JI,JJ-IJJ1+1)/ZSURF
  ENDDO
ENDDO

CALL EINTER1(YDGEOMETRY%YREGEO,ZVA,IDATY,IDATX,ITN,IFLD,ZS(1,5),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!   4 orography standard deviation (sigma, part two): ZS4=sigma**2
!     variance of orography for 2.5>resol>0.5

DO JJ=YRCLI%NDATY,1,-1
  JA=YRCLI%NDATY-JJ+1
  READ(33,REC=JA)IARRAY16
  DO JI=1,YRCLI%NDATX
    ZVA(JI,JJ)=REAL(IARRAY16(JI),JPRB)
    ZVA(JI,JJ)=ZVA(JI,JJ)*ZVA(JI,JJ)
  ENDDO
ENDDO

CALL EINTER1(YDGEOMETRY%YREGEO,ZVA,IDATY,IDATX,ITN,IFLD,ZS(1,4),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!   6 n=nb_peaks: ZS6=z01r= sqrt(n/s) * sigma**2

ZCV=RPI/REAL(YRCLI%NGLOBY,JPRB)
ZCM=2*RPI*RA**2/REAL(YRCLI%NGLOBX,JPRB)
DO JJ=YRCLI%NDATY,1,-1
  JA=YRCLI%NDATY-JJ+1
  JS=JJ+IJJ1-1
  ZSURF=ZCM*(COS(ZCV*(JS-1))-COS(ZCV*JS))
  READ(34,REC=JA)CL_CARRAY
  DO JI=1,YRCLI%NDATX
    I8=ICHAR(CL_CARRAY(JI))
    ZVA(JI,JJ)=ZVA(JI,JJ)*SQRT(REAL(I8,JPRB)/ZSURF)
  ENDDO
ENDDO

CALL EINTER1(YDGEOMETRY%YREGEO,ZVA,IDATY,IDATX,ITN,IFLD,ZS(1,6),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  7 small scale orography variance ZS7=v1=(hmax-h)(h-hmin)/4
!    variance du relief non resolu: relief de resolution <0.5 degres
!    UP TO NOW plays a role only when computing the total orographical
!    variance.
DO JJ=YRCLI%NDATY,1,-1
  JA=YRCLI%NDATY-JJ+1
  READ(39,REC=JA)IARRAY32
  DO JI=1,YRCLI%NDATX
    ZVA(JI,JJ)=REAL(IARRAY32(JI),JPRB)
  ENDDO
ENDDO

CALL EINTER1(YDGEOMETRY%YREGEO,ZVA,IDATY,IDATX,ITN,IFLD,ZS(1,7),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!   8 urbanisation fraction: ZS8=urb

DO JJ=YRCLI%NDATY,1,-1
  JA=YRCLI%NDATY-JJ+1
  READ(35,REC=JA)CL_CARRAY
  DO JI=1,YRCLI%NDATX
    I8=ICHAR(CL_CARRAY(JI))
    ZVA(JI,JJ)=REAL(I8,JPRB)/100._JPRB
  ENDDO
ENDDO

CALL EINTER1(YDGEOMETRY%YREGEO,ZVA,IDATY,IDATX,ITN,IFLD,ZS(1,8),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!   9 TGCT-11 (part one): (dh/dx)^2: ZS9=V11

DO JJ=YRCLI%NDATY,1,-1
  JA=YRCLI%NDATY-JJ+1
  READ(32,REC=JA)IARRAY16
  DO JI=1,YRCLI%NDATX
    ZVA(JI,JJ)=REAL(IARRAY16(JI),JPRB)
  ENDDO
ENDDO

CALL GTOPTX2(IDATY,IDATX,IJJ1,ZVA)
CALL EINTER1(YDGEOMETRY%YREGEO,ZVA,IDATY,IDATX,ITN,IFLD,ZS(1,9),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!   10 TGCT-22 (part one): (dh/dy)^2: ZS10=V22

DO JJ=YRCLI%NDATY,1,-1
  JA=YRCLI%NDATY-JJ+1
  READ(32,REC=JA)IARRAY16
  DO JI=1,YRCLI%NDATX
    ZVA(JI,JJ)=REAL(IARRAY16(JI),JPRB)
  ENDDO
ENDDO

CALL GTOPTY2(IDATY,IDATX,ZVA)
CALL EINTER1(YDGEOMETRY%YREGEO,ZVA,IDATY,IDATX,ITN,IFLD,ZS(1,10),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  11 TGCT-12 (part one): (dh/dx)(dh/dy): ZS11=V12

DO JJ=YRCLI%NDATY,1,-1
  JA=YRCLI%NDATY-JJ+1
  READ(32,REC=JA)IARRAY16
  DO JI=1,YRCLI%NDATX
    ZVA(JI,JJ)=REAL(IARRAY16(JI),JPRB)
  ENDDO
ENDDO

CALL GTOPTXY(IDATY,IDATX,IJJ1,IJJ2,ZVA)
CALL EINTER1(YDGEOMETRY%YREGEO,ZVA,IDATY,IDATX,ITN,IFLD,ZS(1,11),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  12 TGCT-11 (part two): ZS12=nu11

DO JJ=YRCLI%NDATY,1,-1
  JA=YRCLI%NDATY-JJ+1
  READ(37,REC=JA)IARRAY32
  DO JI=1,YRCLI%NDATX
    ZVA(JI,JJ)=REAL(IARRAY32(JI),JPRB)/1.E06_JPRB
  ENDDO
ENDDO

CALL EINTER1(YDGEOMETRY%YREGEO,ZVA,IDATY,IDATX,ITN,IFLD,ZS(1,12),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  13 TGCT-22 (part two): ZS13=nu22

DO JJ=YRCLI%NDATY,1,-1
  JA=YRCLI%NDATY-JJ+1
  READ(38,REC=JA)IARRAY32
  DO JI=1,YRCLI%NDATX
    ZVA(JI,JJ)=REAL(IARRAY32(JI),JPRB)/1.E06_JPRB
  ENDDO
ENDDO

CALL EINTER1(YDGEOMETRY%YREGEO,ZVA,IDATY,IDATX,ITN,IFLD,ZS(1,13),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  14 TGCT-12 (part two): ZS14=nu12

DO JJ=YRCLI%NDATY,1,-1
  JA=YRCLI%NDATY-JJ+1
  READ(36,REC=JA)IARRAY32
  DO JI=1,YRCLI%NDATX
    ZVA(JI,JJ)=REAL(IARRAY32(JI),JPRB)/1.E06_JPRB
  ENDDO
ENDDO

CALL EINTER1(YDGEOMETRY%YREGEO,ZVA,IDATY,IDATX,ITN,IFLD,ZS(1,14),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!     ------------------------------------------------------------------

!     4. CLOSE FILES AND DEALLOCATE SPACE.
!        ---------------------------------

CLOSE(31)
CLOSE(32)
CLOSE(33)
CLOSE(34)
CLOSE(35)
CLOSE(36)
CLOSE(37)
CLOSE(38)
CLOSE(39)

DEALLOCATE ( ZVA  )

!     ------------------------------------------------------------------

!     5. FINAL CALCULATIONS AND WRITING.
!        -------------------------------

! FINAL ATTRIBUTIONS of LOCAL ARRAY ZS:
!    ZS1: land-sea mask at the end
!    ZS2: mo: mean orography with possible envelope (*g at the end)
!    ZS3: standard deviation of orography (*g at the end)
!    ZS4: gamma returned by EGANISO (anisotropy coefficient)
!    ZS5: theta returned by EGANISO (principal direction)
!    ZS6: Z0R orographical roughness length (*g at the end)
!    ZS7: land fraction (1-%w)
!    ZS8: urb: urbanisation fraction
! INTERMEDIATE ATTRIBUTIONS
! Resolved orography covariance tensor:
!    ZS9:  V11: TGCT-11a
!    ZS10: V22: TGCT-22a
!    ZS11: V12: TGT-12a
! Unresolved orography covariance tensor:
!    ZS12: nu11: TGTC11b
!    ZS13: nu22: TGTC22b
!    ZS14: nu12: TGTC12b

DO J=1,ITFING

!     5.1 Calculate V=mo2-mo^2 resolved orography variance

  ZS(J,3)=ZS(J,3) - ZS(J,2)**2

!     5.2 Calculate TGC Tensor Main dir and anisotropy

! Z1=V11+V22, Z2=nu11+nu22
  Z1=ZS(J, 9)+ZS(J, 10)
  Z2=ZS(J,12)+ZS(J,13)
  IF (Z1 < ZEPS) THEN
    ZS(J, 9)=0.0_JPRB
    ZS(J,10)=0.0_JPRB
    ZS(J,11)=0.0_JPRB
  ELSE
! +=V/Z1
    ZS(J,9 )=ZS(J, 9)*ZS(J,3)/Z1
    ZS(J,10)=ZS(J,10)*ZS(J,3)/Z1
    ZS(J,11)=ZS(J,11)*ZS(J,3)/Z1
  ENDIF
  IF (Z2 >= ZEPS) THEN
! +=sigma/Z2
    ZS(J,9 )=ZS(J, 9)+ZS(J,4)*ZS(J,12)/Z2
    ZS(J,10)=ZS(J,10)+ZS(J,4)*ZS(J,13)/Z2
    ZS(J,11)=ZS(J,11)+ZS(J,4)*ZS(J,14)/Z2
  ENDIF

!     5.3 Calculate the rest

! ZS6= Z0R = z01r + sqrt(N/s) * V
  ZS(J,6)=ZS(J,6)+SQRT(MAX(0.0_JPRB,ZS(J,5)))*ZS(J,3)

! ZS3=sqrt(V + sigma + v1)=orography std dev
  ZS(J,3)=SQRT(MAX(0.0_JPRB,ZS(J,3)+ZS(J,4)+ZS(J,7)))
  
! save the orography without envelope
  ZZWI(J,1)=ZS(J,2)

! possible addition to yield envelope orography:
  ZS(J,2)=ZS(J,2)+FENVN*ZS(J,3)*(1.0_JPRB-ZS(J,1))

! ZS7 = land fraction
  ZS(J,7)=MAX(0.0_JPRB,MIN(1.0_JPRB,1.0_JPRB-ZS(J,1)))

! ZS8 = urbanisation fraction
  ZS(J,8)=MAX(0.0_JPRB,MIN(ZS(J,7),ZS(J,8)))

! Apply Land-Sea mask on Z0R:
  IF (ZS(J,1) <= YRCLI%SMASK) THEN
    ZS(J,6)=MAX(FACZ0*ZS(J,6),YRCLI%SZZ0N)
  ELSE
    ZS(J,6)=YRCLI%SZZ0M
  ENDIF

ENDDO

!     5.4 Calculate anisotropy parameters, smooth roughness length

CALL EGANISO(YDGEOMETRY,ITFING,ZS(1,9),ZS(1,10),ZS(1,11),ZS(1,4),ZS(1,5))

!  Smoothing of roughness lengths.
!  The ZS6 field, still in meters, will be overwritten by its smoothed
!  equivalent. Sea values values will be preserved during the process.

DO J=1,ITFING
  ZS(J,6)=LOG(ZS(J,6))
ENDDO
DO JLOOP=1,NLISSZ
  CALL ELISLAP(YDGEOMETRY%YREGEO,ZS(1,6),IYFING,IXFING,ITFING,ZS(1,1),ZGM,YRCLI%SMASK,ZPOZ)
ENDDO
DO J=1,ITFING
  ZS(J,6)=EXP(ZS(J,6))
ENDDO

!     5.5 Fit of the orography

! Smooth orography if asked
DO JLOOP=1,NLISSR
  CALL ELISLAP(YDGEOMETRY%YREGEO,ZS(1,2),IYFING,IXFING,ITFING,ZS(1,1),ZGM,YRCLI%SMASK,ZPOR)
ENDDO

IF (LNLSM.AND..NOT.LNORO) THEN
  WRITE(NULOUT,*) ' IF YOU WANT TO IMPORT A NEW LAND-SEA MASK,'
  WRITE(NULOUT,*) ' ASK FOR A NEW OROGRAPHY !'
  CALL ABOR1('EINCLI1: NEW LAND-SEA MASK MUST BE USED WITH A NEW OROGRAPHY !')
ENDIF

LLFIT=.TRUE.

! Read another orography and a new land-sea mask -if asked- on a ALADIN file.
IF (LNORO) THEN
  CALL RELNEW(YDGEOMETRY,LNLSM,LLFIT,ITFING,ZS(1,2),ZS(1,1))
ENDIF
! Add envelop if imported PGD orography
IF (LIPGD) THEN
  DO J=1,ITFING
    ZS(J,2)=ZS(J,2)+FENVN*ZS(J,3)*(1.0_JPRB-ZS(J,1))
  ENDDO
ENDIF

! Check the number of points over land, set orography to 0 if sea everywhere
ITR=0
DO J=1,ITFING
  IF (ZS(J,1) <= YRCLI%SMASK) ITR= ITR+1
ENDDO 
IF (ITR == 0) THEN
  LLFIT=.FALSE.
! Gridpoint and spectral orography already set to 0
ENDIF

! Spectral fit to the orography, unless  :
!   - either a new orography is read as spectral coefficients,
!   - or the domains covers sea only.
! The spectrally fitted equivalent of the ZS(.,2) field, still in meters,
! is stored in *XREV*. The spectral coefficients are in *XRNM*

IF (LLFIT) THEN

! Biperiodization
  DO J=1,ITFING
    ZWR(J,1)=ZS(J,2)
  ENDDO

  LLBI(1)=.TRUE.
  LLWR(1)=.FALSE.
  CALL EBICLI(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF,1,INIV,CLPREF(9),CLSUFF(9),INUM,ZWR,LLBI,LLWR,LLPA)

  DO JY=1,NDGLG
    DO JX=1,NDLON
      IADR=JX+IADL(JY)
      XREV(IADR)=ZWR((JY-1)*NDLON+JX,1)
    ENDDO
  ENDDO
! DO J=1,NGPTOT
!   XREV(J)=ZWR(J,1)
! ENDDO

!* Spectral fit.
  IF (LKEYF) THEN

! Weights for the cost functions
    CALL INIPZ(YDGEOMETRY,ZS(1,1),YRCLI%SMASK,ITFING)
    CALL INIRP(YDGEOMETRY%YRLAP,YDGEOMETRY%YRELAP,YDGEOMETRY%YRDIM)

! Raw spectral fit or optimization, according to LLOPT
    WRITE(NULOUT,*) ' SPECTRAL FIT'
    LLOPT= (LNEWORO .OR. LNEWORO2) .AND. .NOT.LSPSMORO
    CALL RELSPE(YDGEOMETRY,LLOPT)

  ENDIF

ENDIF

!     5.6 Translation into final units.

DO J=1,ITFING
! Land-Sea mask:
  ZS(J,1)=MAX(0.0_JPRB,SIGN(1.0_JPRB,YRCLI%SMASK-ZS(J,1)))
! Std dev of orography
  ZS(J,3)=RG*ZS(J,3)
! Roughness length of bare surface:
  ZS(J,6)=RG*ZS(J,6)
ENDDO

! Spectral orography:
IF (LKEYF) THEN
  DO JS=1,NSEFRE
    XRNM(JS)=RG*XRNM(JS)
  ENDDO
ELSE
! Spectral orography written as missing data
  XRNM(:)=YRCLI%SMANQ - 1.0_JPRB
ENDIF

! Grid-point orography, scaling and reordering:
DO J=1,NGPTOT
  XREV(J)=RG*XREV(J)
ENDDO
DO JY=1,NDGLG
  DO JX=1,NDLON
    IADR=JX+IADL(JY)
    ZWR((JY-1)*NDLON+JX,1)=XREV(IADR)
  ENDDO
ENDDO
!DO J=1,NGPTOT
!  ZWR(J,1)=XREV(J)
!ENDDO

!     5.7 Final writing on ALADIN file as result.

! Set-up for ALADIN file

! CADRE
CLNOMC='Const.Clim.Surfa'
ITYPTR=-NMSMAX
ITRONC=NSMAX
LLKEEP=.TRUE.

ZCODIL = 0.0_JPRB
IF (.NOT.LMAP) ZCODIL = -1.0_JPRB

DO J=1,NDGNH+1
  INOZPA(J)=0
  ZSINLA(J)=0.0_JPRB
ENDDO

DO J=1,NDGLG+2
  INIP(J) = 0
ENDDO

IF(NCADFORM == 0) THEN
  IF (LMRT) THEN
    CALL ABOR1('SUFRAME : Mercator Rot/Til and OLD CADRE not allowed !')
  ENDIF
  ZSINLA(1)  = 0.0_JPRB              ! ex-NROTEQ
  ZSINLA(2)  = 0.0_JPRB              ! ex-ELONR
  ZSINLA(3)  = 0.0_JPRB              ! ex-ELATR
  ZSINLA(4)  = ELON1
  ZSINLA(5)  = ELAT1
  ZSINLA(6)  = ELON2
  ZSINLA(7)  = ELAT2
  ZSINLA(8)  = ELON0
  ZSINLA(9)  = ELAT0
  ZSINLA(10) = ERPK
  ZSINLA(11) = 0.0_JPRB              ! ex-NSOTRP
  ZSINLA(12) = 0.0_JPRB              ! ex-NGIV0
  ZSINLA(13) = ELX
  ZSINLA(14) = ELY
  ZSINLA(15) = EDELX
  ZSINLA(16) = EDELY
  ZSINLA(17) = EXWN
  ZSINLA(18) = EYWN
ELSE
  ZSINLA(1)  = -1.0_JPRB
  IF (LMRT) ZSINLA(1)  = ZSINLA(1) - 1.0_JPRB ! Mercator Rot/Tilted = -2
  ZSINLA(2)  = ERPK
  ZSINLA(3)  = ELON0
  ZSINLA(4)  = ELAT0
  ZSINLA(5)  = ELONC
  ZSINLA(6)  = ELATC
  ZSINLA(7)  = EDELX
  ZSINLA(8)  = EDELY
  ZSINLA(9)  = ELX
  ZSINLA(10) = ELY
  ZSINLA(11) = EXWN
  ZSINLA(12) = EYWN
  ZSINLA(13) = ELON1
  ZSINLA(14) = ELAT1
  ZSINLA(15) = ELON2
  ZSINLA(16) = ELAT2
  ZSINLA(17) = 0.0_JPRB              ! free
  ZSINLA(18) = 0.0_JPRB              ! free
ENDIF

INIP(1) = 10
INIP(2) = 1
IF (.NOT.LKEYF) INIP(2) = -1
INIP(3) = NDLUNG
INIP(4) = NDLUXG
INIP(5) = NDGUNG
INIP(6) = NDGUXG
INIP(7) = NBZONL
INIP(8) = NBZONG

CALL FACADE (CLNOMC,ITYPTR,0.0_JPRB,0.0_JPRB,0.0_JPRB,ZCODIL,&
 & ITRONC,NDGLG,NDLON,INIP(1),INOZPA(1),ZSINLA(1),&
 & NFLEVG,VP00,YDVAB%VALH,YDVAB%VBH,LLKEEP)  

! Open file
IREP=0
IMES=1
INBARI=0
INBARP=10
LLNOMM=.TRUE.
LLERFA=.TRUE.
LLIMST=.TRUE.

CALL FAITOU(IREP,INUM,LLNOMM,CLNOMF,'NEW',LLERFA,LLIMST,IMES,&
 & INBARP,INBARI,CLNOMC)  

! Set date
DO J=1,11
  IDATEF(J)=0
ENDDO
IDATEF(1)=1
IDATEF(2)=1
IDATEF(3)=15
IDATEF(6)=1

CALL FANDAR(IREP,INUM,IDATEF)

! Get GRIB options
IGRIB=0
IGRIG=0
IBPDG=0
IBCSP=0
ITRON=0
ILAP=0
IDMO=0

CALL FAVEUR(IREP,INUM,IGRIB,IBPDG,IBCSP,ITRON,ILAP,IDMO)
! new spectral ordering only if packing type is -1 or 3
IF (IGRIB == -1 .OR. IGRIB == 3) THEN
  IGRIG=-1
ELSE
  IGRIG=0
ENDIF

! Write spectral and gridpoint orography (2 fields)

! Change GRIB level of coding to "none"
CALL FAGOTE(IREP,INUM,IGRIG,IBPDG,IBCSP,ITRON,ILAP,IDMO)
! Grid-point orography
LLBI(1)=.FALSE.
LLWR(1)=.TRUE.
CALL EBICLI(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF,1,INIV,CLPREF(9),CLSUFF(9),INUM,ZWR,LLBI,LLWR,LLPA)
! Spectral orography
LLCOSP=.TRUE.
! FT: replace FAIENO by FAIENC here : FAIENO reorders the spectral arrays,
! and as it is already done in 923 it leads to a desordering 
! but Francoise said that it could also be done with a clean use of FAIENO
! (will be necessary to be able to produce grib2 files)
CALL FAIENC(IREP,INUM,CLPREF(0),INIV(1),CLSUFF(0),XRNM,LLCOSP)
! Reset previous GRIB options
CALL FAGOTE(IREP,INUM,IGRIB,IBPDG,IBCSP,ITRON,ILAP,IDMO)

! Write other gridpoint fields (8)

! Transfer to C+I+E arrays
DO JC=1,8
  IC=JC+1
  IF (JC == 1) IC=1
  IF (JC == 8) IC=6
  DO J=1,ITFING
    ZEXT(J,JC)=ZS(J,IC)
  ENDDO
ENDDO
! Biperiodization and writing
CALL EBICLI(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF,8,INIVL,CLPREF(1),CLSUFF(1),INUM,ZEXT,LLBIP,LLWRI,LLPAC)

! Close file
CALL LFILAF(IREP,INUM,.TRUE.)
CALL FAIRME(IREP,INUM,'KEEP')

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EINCLI1',1,ZHOOK_HANDLE)
END SUBROUTINE EINCLI1

