SUBROUTINE INCLI1(YDGEOMETRY)

!**** *INCLI1*

!     PURPOSE.
!     --------

!     This routine calculates first 7 fixed fields :
!      .land(1)/water(0) mask,
!      .surface geopotential (grid point and spectral)
!        -including enveloppe effect, if asked-,
!      .G times the standard deviation of the orography,
!      .the anisotropy coefficient
!      .the direction of principal axis of topography (in radian),
!      .G times the roughness length of topography
!        -with the possibility of scaling it by an arbitrary parameter-,
!      .the fraction of land.
!      .the fraction of urbanization.
!     The target grid is a Gaussian grid with rotation of the pole, and
!     "Schmidt"
!     stretching as well as -if asked- reduction of the number of grid points
!     towards the poles of the final representation.
!     Input data are coming from the *MANU* type files, based on GLOBE data

!**   INTERFACE.
!     ----------

!     CALL INCLI1

!     METHOD.
!     -------

!     The *MANU* dataset is read on units 31-39.
!     A subgrid of the gaussian grid is introduced, and a value on this subgrid
!     is the value of the nearest point on the initial grid.
!     The values on the gaussian grid are the averages of the values in the 
!     boxes (a box is the part of the subgrid corresponding to a point of the
!     gaussian grid).
!     The 7 "fixed" output fields are written on unit n=41 (one output file).

!     EXTERNALS.
!     ----------

!      GEO923
!      INTER1
!      LOCMAXI (or LOCMAXP)
!      GTOPTX2
!      GTOPTXY
!      GTOPTY2
!      GANISO
!      RELNEW
!      RELSPE
!      INIRP
!      LISLAP
!      SPREORD
!      FA-LFI package (FACADE,FAITOU,FANDAR,FAVEUR,FAGOTE,FAIENC,LFILAF,FAIRME)

!     AUTHORS.
!     --------
!      E. Cordoneanu 96-12-15 from INCLIA

!     MODIFICATIONS.
!     --------------
!      D. Giard 01-11-19 : message passing
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib : 03-08-18 Roughness lengths not packed
!      D. Paradis & R. El Khatib : 04-07-22 GRIBEX
!      D. Giard 04-01-07 : ARPEGE-ALADIN consistency for spectral fit
!      D. Giard 05-04-04 : update and cleaning
!                          semi-envelope suppressed
!                          raw or no spectral fit added
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB, JPIS
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMDIL   , ONLY : SLAPO, SLOPO, GLOPO, FACDI
USE YOMCST   , ONLY : RPI, RA, RG
USE YOMLUN   , ONLY : NULOUT
USE YOMVERT  , ONLY : VP00
USE YOMCLI   , ONLY : YRCLI
USE YOMCLA   , ONLY : NLISSZ, FACZ0, FENVN, FENVS, LNORO, LNLSM, LNEWORO,&
 & LNEWORO2, LKEYF
USE PTRSPOR  , ONLY : XRNM, XREV

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT)   :: YDGEOMETRY
REAL(KIND=JPRB) :: ZS(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON,14),ZSINLA(YDGEOMETRY%YRDIM%NDGNH)&
 & ,ZBO(YDGEOMETRY%YRDIM%NDGLG+1),ZLONG(YDGEOMETRY%YRDIM%NDGLG),ZMU(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON),&
 & ZSLA(YDGEOMETRY%YRDIM%NDGLG*YRCLI%NPINT)&
 & ,ZSLO(YDGEOMETRY%YRDIM%NDLON*YRCLI%NPINT,YDGEOMETRY%YRDIM%NDGLG),ZCLO(YDGEOMETRY%YRDIM%NDLON*YRCLI%NPINT,YDGEOMETRY%YRDIM%NDGLG)  
REAL(KIND=JPRB),ALLOCATABLE :: ZVA(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZRNM(:)
REAL(KIND=JPRB) :: Z, Z1, Z2, ZCM, ZCV, ZDENV, ZEPS, ZMENV, ZPOZ, ZSURF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: INLOPA(YDGEOMETRY%YRDIM%NDGNH),INOZPA(YDGEOMETRY%YRDIM%NDGNH),IDATEF(11),IARRAY32(YRCLI%NDATX)
INTEGER(KIND=JPIS) :: IARRAY16(YRCLI%NDATX)
INTEGER(KIND=JPIM) :: I8, IADR, IARI, IARP, IBCSP, IBPDG, ICO, IFLD,&
 & IGRIBC, IGRIBO, ILAP, ILEND, ILENE, ILENT,&
 & IMES, IMOD, INDC, INDK, INIQ, INIV, INJQ, INSIGN,&
 & INUM, IREP, ITRON, J, JJ  

CHARACTER :: CLNOMC*16,CLNOMF*10,CLPREF(0:9)*8,CLSUFF(0:9)*12
CHARACTER (LEN = 1) ::  CLARRAY(YRCLI%NDATX)

LOGICAL :: LLCOSP, LLFIT, LLIMST, LLKEEP, LLOPT
LOGICAL :: LLPACK(0:9)

#include "abor1.intfb.h"
#include "ganiso.intfb.h"
#include "geo923.intfb.h"
#include "gtoptx2.intfb.h"
#include "gtoptxy.intfb.h"
#include "gtopty2.intfb.h"
#include "inipz.intfb.h"
#include "inirp.intfb.h"
#include "inter1.intfb.h"
#include "lislap.intfb.h"
#include "locmaxi.intfb.h"
#include "relnew.intfb.h"
#include "relspe.intfb.h"
#include "spreord.intfb.h"

!     ------------------------------------------------------------------

!     1. SET INITIAL VALUES AND DEFINE THE SUBGRID.
!        ------------------------------------------

IF (LHOOK) CALL DR_HOOK('INCLI1',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDLAP=>YDGEOMETRY%YRLAP, &
& YDCSGLEG=>YDGEOMETRY%YRCSGLEG, YDVAB=>YDGEOMETRY%YRVAB, YDELAP=>YDGEOMETRY%YRELAP)
ASSOCIATE(NDGLG=>YDDIM%NDGLG, NDGNH=>YDDIM%NDGNH, NDLON=>YDDIM%NDLON,   NSEFRE=>YDDIM%NSEFRE, NSMAX=>YDDIM%NSMAX,      &
& NSPEC2G=>YDDIM%NSPEC2G,   NFLEVG=>YDDIMV%NFLEVG,   NLOENG=>YDGEM%NLOENG, NMENG=>YDGEM%NMENG, NSTAGP=>YDGEM%NSTAGP,   &
& NSTTYP=>YDGEM%NSTTYP)

ZEPS=1.E-6_JPRB

! Characteristic length for smoothing z0
ZPOZ=4000._JPRB

!     1.1 Dimensions for data and interpolation.

ICO=YRCLI%NPINT*YRCLI%NPINT
INJQ=NDGLG*YRCLI%NPINT
INIQ=NDLON*YRCLI%NPINT
ILENT=NDGLG*NDLON
ILEND=YRCLI%NDATX*YRCLI%NDATY
IFLD=1
ILENE=0

!     1.2 Model grid.

CALL GEO923(YDGEOMETRY,YRCLI%NPINT,ILENE,ZBO,ZLONG,ZMU,ZSLA,ZSLO,ZCLO)

!     ------------------------------------------------------------------

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

ALLOCATE ( ZVA (YRCLI%NDATX,YRCLI%NDATY) )

!     ------------------------------------------------------------------

!     3. INTERPOLATION IN MODEL GRID POINTS.
!        -----------------------------------

!   1 fraction of water covered area
!   2 mean orography
!   3 orographycal and urbanisation roughness length (part 3)
!   4 orography standard deviation                   (part 2)
!   5 orographical and urbanisation roughness length (part 1)
!   6 orographical and urbanisation roughness length (part 2)
!   7 square mean orography
!  Topographic Gradient Correlation Tensor(TGCT)
!   8 TGCT-11                                        (part 1)
!   9 TGCT-22                                        (part 1)
!  10 TGCT-12                                        (part 1)
!  11 TGCT-11                                        (part 2)
!  12 TGCT-22                                        (part 2)
!  13 TGCT-12                                        (part 2)
!  14 orography variance                             (part 3)

!   1 fraction of water covered area

DO JJ=1,YRCLI%NDATY
  READ(31,REC=JJ) CLARRAY
  DO J=1,YRCLI%NDATX
    I8=ICHAR(CLARRAY(J))
    ZVA(J,JJ)=REAL(I8,JPRB)/100._JPRB
  ENDDO
ENDDO
CALL INTER1(NDGLG,NDLON,YRCLI%NDATY,YRCLI%NDATX,IFLD,YRCLI%NPINT,ILENE,ILEND,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,1),ZVA,ZSLA,ZSLO,ZCLO)  

!   2 mean orography

DO JJ=1,YRCLI%NDATY
  READ(32,REC=JJ) IARRAY16
  DO J=1,YRCLI%NDATX
    ZVA(J,JJ)=REAL(IARRAY16(J),JPRB)
  ENDDO
ENDDO
CALL INTER1(NDGLG,NDLON,YRCLI%NDATY,YRCLI%NDATX,IFLD,YRCLI%NPINT,ILENE,ILEND,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,2),ZVA,ZSLA,ZSLO,ZCLO)  

!   7 square mean orography

DO JJ=1,YRCLI%NDATY
  DO J=1,YRCLI%NDATX
    ZVA(J,JJ)=ZVA(J,JJ)*ZVA(J,JJ)
  ENDDO
ENDDO
CALL INTER1(NDGLG,NDLON,YRCLI%NDATY,YRCLI%NDATX,IFLD,YRCLI%NPINT,ILENE,ILEND,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,7),ZVA,ZSLA,ZSLO,ZCLO)  

!   5 orographical and urbanisation roughness length (part one)

!  Search for local maxima.
CALL LOCMAXI(YRCLI%NDATY,YRCLI%NDATX,ZVA)
!     CALL LOCMAXP(NDATY,NDATX,ZVA)

ZCV=RPI/REAL(YRCLI%NDATY,JPRB)
ZCM=2*RPI*RA**2/REAL(YRCLI%NDATX,JPRB)
DO JJ=1,YRCLI%NDATY
  ZSURF=ZCM*(COS(ZCV*(JJ-1))-COS(ZCV*JJ))
  DO J=1,YRCLI%NDATX
    ZVA(J,JJ)=ZVA(J,JJ)/ZSURF
  ENDDO
ENDDO
CALL INTER1(NDGLG,NDLON,YRCLI%NDATY,YRCLI%NDATX,IFLD,YRCLI%NPINT,ILENE,ILEND,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,5),ZVA,ZSLA,ZSLO,ZCLO)  

!   8 TGCT-11 (part one)

DO JJ=1,YRCLI%NDATY
  READ(32,REC=JJ) IARRAY16
  DO J=1,YRCLI%NDATX
    ZVA(J,JJ)=REAL(IARRAY16(J),JPRB)
  ENDDO
ENDDO
CALL GTOPTX2(YRCLI%NDATY,YRCLI%NDATX,1,ZVA)
CALL INTER1(NDGLG,NDLON,YRCLI%NDATY,YRCLI%NDATX,IFLD,YRCLI%NPINT,ILENE,ILEND,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,8),ZVA,ZSLA,ZSLO,ZCLO)  

!   9 TGCT-22 (part one)

DO JJ=1,YRCLI%NDATY
  READ(32,REC=JJ) IARRAY16
  DO J=1,YRCLI%NDATX
    ZVA(J,JJ)=REAL(IARRAY16(J),JPRB)
  ENDDO
ENDDO
CALL GTOPTY2(YRCLI%NDATY,YRCLI%NDATX,ZVA)
CALL INTER1(NDGLG,NDLON,YRCLI%NDATY,YRCLI%NDATX,IFLD,YRCLI%NPINT,ILENE,ILEND,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,9),ZVA,ZSLA,ZSLO,ZCLO)  

!  10 TGCT-12 (part one)

DO JJ=1,YRCLI%NDATY
  READ(32,REC=JJ) IARRAY16
  DO J=1,YRCLI%NDATX
    ZVA(J,JJ)=REAL(IARRAY16(J),JPRB)
  ENDDO
ENDDO
CALL GTOPTXY(YRCLI%NDATY,YRCLI%NDATX,1,YRCLI%NDATY,ZVA)
CALL INTER1(NDGLG,NDLON,YRCLI%NDATY,YRCLI%NDATX,IFLD,YRCLI%NPINT,ILENE,ILEND,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,10),ZVA,ZSLA,ZSLO,ZCLO)  

!   4 orography standard deviation (part two)

DO JJ=1,YRCLI%NDATY
  READ(33,REC=JJ) IARRAY16
  DO J=1,YRCLI%NDATX
    ZVA(J,JJ)=REAL(IARRAY16(J),JPRB)
    ZVA(J,JJ)=ZVA(J,JJ)*ZVA(J,JJ)
  ENDDO
ENDDO
CALL INTER1(NDGLG,NDLON,YRCLI%NDATY,YRCLI%NDATX,IFLD,YRCLI%NPINT,ILENE,ILEND,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,4),ZVA,ZSLA,ZSLO,ZCLO)  

!   6 orographical and urbanisation roughness length (part two)

ZCV=RPI/REAL(YRCLI%NDATY,JPRB)
ZCM=2*RPI*RA**2/REAL(YRCLI%NDATX,JPRB)
DO JJ=1,YRCLI%NDATY
  ZSURF=ZCM*(COS(ZCV*(JJ-1))-COS(ZCV*JJ))
  READ(34,REC=JJ) CLARRAY
  DO J=1,YRCLI%NDATX
    I8=ICHAR(CLARRAY(J))
    ZVA(J,JJ)=ZVA(J,JJ)*SQRT(REAL(I8,JPRB)/ZSURF)
  ENDDO
ENDDO
CALL INTER1(NDGLG,NDLON,YRCLI%NDATY,YRCLI%NDATX,IFLD,YRCLI%NPINT,ILENE,ILEND,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,6),ZVA,ZSLA,ZSLO,ZCLO)  

!   3 orographical and urbanisation roughness length (part three)

DO JJ=1,YRCLI%NDATY
  READ(35,REC=JJ) CLARRAY
  DO J=1,YRCLI%NDATX
    I8=ICHAR(CLARRAY(J))
    ZVA(J,JJ)=REAL(I8,JPRB)/100._JPRB
  ENDDO
ENDDO
CALL INTER1(NDGLG,NDLON,YRCLI%NDATY,YRCLI%NDATX,IFLD,YRCLI%NPINT,ILENE,ILEND,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,3),ZVA,ZSLA,ZSLO,ZCLO)  

!  13 TGCT-12 (part two)

DO JJ=1,YRCLI%NDATY
  READ(36,REC=JJ) IARRAY32
  DO J=1,YRCLI%NDATX
    ZVA(J,JJ)=REAL(IARRAY32(J),JPRB)/1.E06_JPRB
  ENDDO
ENDDO
CALL INTER1(NDGLG,NDLON,YRCLI%NDATY,YRCLI%NDATX,IFLD,YRCLI%NPINT,ILENE,ILEND,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,13),ZVA,ZSLA,ZSLO,ZCLO)  

!  11 TGCT-11 (part two)

DO JJ=1,YRCLI%NDATY
  READ(37,REC=JJ) IARRAY32
  DO J=1,YRCLI%NDATX
    ZVA(J,JJ)=REAL(IARRAY32(J),JPRB)/1.E06_JPRB
  ENDDO
ENDDO
CALL INTER1(NDGLG,NDLON,YRCLI%NDATY,YRCLI%NDATX,IFLD,YRCLI%NPINT,ILENE,ILEND,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,11),ZVA,ZSLA,ZSLO,ZCLO)  

!  12 TGCT-22 (part two)

DO JJ=1,YRCLI%NDATY
  READ(38,REC=JJ) IARRAY32
  DO J=1,YRCLI%NDATX
    ZVA(J,JJ)=REAL(IARRAY32(J),JPRB)/1.E06_JPRB
  ENDDO
ENDDO
CALL INTER1(NDGLG,NDLON,YRCLI%NDATY,YRCLI%NDATX,IFLD,YRCLI%NPINT,ILENE,ILEND,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,12),ZVA,ZSLA,ZSLO,ZCLO)  

!  14 orography variance (part three) 

DO JJ=1,YRCLI%NDATY
  READ(39,REC=JJ) IARRAY32
  DO J=1,YRCLI%NDATX
    ZVA(J,JJ)=MAX(REAL(IARRAY32(J),JPRB),0.0_JPRB)
  ENDDO
ENDDO
CALL INTER1(NDGLG,NDLON,YRCLI%NDATY,YRCLI%NDATX,IFLD,YRCLI%NPINT,ILENE,ILEND,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZS(1,14),ZVA,ZSLA,ZSLO,ZCLO)  

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

!     5.1 Physical computations on the reduced gaussian grid.

!   1 Land/Sea index (1=land, but only at the end of the routine)
!   2 Orography (times g, idem)
!   3 Standard deviation of orography (times g, idem)
!   4 Anisotropy coefficient of topography
!   5 Direction of the pricipal axis of the topography
!   6 Roughness length of topography (times g, idem)    
!   7 Fraction of land
!   8 Fraction of urbanization (at the end of the routine)
!   9 Mean orography (without envelope)

DO J=1,ILENE
  Z=ZS(J,3)
  ZS(J,3)=ZS(J,7) - ZS(J,2)**2
  ZS(J,7)=Z
  Z1=ZS(J, 8)+ZS(J, 9)
  Z2=ZS(J,11)+ZS(J,12)
  IF (Z1 < ZEPS) THEN
    ZS(J, 8)=0.0_JPRB
    ZS(J, 9)=0.0_JPRB
    ZS(J,10)=0.0_JPRB
  ELSE
    ZS(J,8 )=ZS(J, 8)*ZS(J,3)/Z1
    ZS(J,9 )=ZS(J, 9)*ZS(J,3)/Z1
    ZS(J,10)=ZS(J,10)*ZS(J,3)/Z1
  ENDIF
  IF (Z2 >= ZEPS) THEN
    ZS(J,8 )=ZS(J, 8)+ZS(J,4)*ZS(J,11)/Z2
    ZS(J,9 )=ZS(J, 9)+ZS(J,4)*ZS(J,12)/Z2
    ZS(J,10)=ZS(J,10)+ZS(J,4)*ZS(J,13)/Z2
  ENDIF
ENDDO

!  Calculation of the degree of anisotropy and principal axis of topography

CALL GANISO(YDGEOMETRY,ILENE,ZS(1,8),ZS(1,9),ZS(1,10),ZS(1,12),ZS(1,13))

ZDENV=0.5_JPRB*(FENVN-FENVS)
ZMENV=0.5_JPRB*(FENVN+FENVS)
DO J=1,ILENE
  ZS(J,6)=ZS(J,6) + SQRT( MAX(0.0_JPRB,ZS(J,5)) )*ZS(J,3)
  ZS(J,3)=SQRT( MAX(0.0_JPRB,ZS(J,3)+ZS(J,4)+ZS(J,14)) )
  ZS(J,9)=ZS(J,2)
  ZS(J,2)=ZS(J,2) + (ZDENV*ZMU(J)+ZMENV)*ZS(J,3)*(1.0_JPRB-ZS(J,1))
  ZS(J,8)=ZS(J,7)
  ZS(J,7)=ZS(J,1)
  IF (ZS(J,1) <= YRCLI%SMASK) THEN
    ZS(J,6)=MAX(FACZ0*ZS(J,6),YRCLI%SZZ0N)
  ELSE
    ZS(J,6)=YRCLI%SZZ0M
  ENDIF
  ZS(J,4)=ZS(J,12)
  ZS(J,5)=ZS(J,13)
ENDDO

!  Smoothing of roughness length.
!  The ZS(.,6) field, still in meters, will be overwritten by its smoothed
!  equivalent. Sea values will be preserved during the process.
!  ZBO is overwritten here.
DO JJ=1,NDGLG
  INDC=MIN(JJ,NDGLG+1-JJ)
  INSIGN=SIGN(1,NDGNH-JJ)
  ZBO(JJ)=REAL(INSIGN,JPRB)*ASIN(YDCSGLEG%RMU(INDC))
ENDDO

DO J=1,ILENE
  ZS(J,6)=LOG(ZS(J,6))
ENDDO
DO J=1,NLISSZ
  CALL LISLAP(ZS(1,6),NLOENG(1:),NDGLG,NDLON,ILENE,ZLONG,ZBO,FACDI,&
   & ZS(1,1),YRCLI%SMASK,ZPOZ)  
ENDDO
DO J=1,ILENE
  ZS(J,6)=EXP(ZS(J,6))
ENDDO

! Default values for spectral orography

XRNM(:)=0.0_JPRB
XREV(:)=0.0_JPRB

IF (LNORO) THEN

!  In this case another orography and a new land-sea mask -if asked- is used.
!  It is read on an ARPEGE file and, after checking the *CADRE*, it undergoes
!  a spectral truncation. Final fields are stored in XRNM and XREV.
  LLFIT=.TRUE.
  CALL RELNEW(YDGEOMETRY,LNLSM,LLFIT,ILENE,ZS(1,2),ZS(1,1))

ELSE 

! Initialization of gridpoint orography, stored in XREV
  INDK=0
  DO JJ=1,NDGLG
    DO J=1,NLOENG(JJ)
      INDC=INDK+J
      IADR=J+NSTAGP(JJ)-1
      XREV(IADR)=ZS(INDC,2)
    ENDDO
    INDK=INDK+NLOENG(JJ)
  ENDDO

  IF (LKEYF) THEN
! Computation and optimization of spectral orography
! The spectrally fitted equivalent of the ZS(.,2) field, still in meters,
! is stored in *XREV*. The spectral coefficients are in *XRNM*.
! The weights for the cost-functions are computed by INIPZ and INIRP.
! Spectral fit is performed in RELSPE, with or without optimization.

    CALL INIPZ(YDGEOMETRY,ZS(1,1),YRCLI%SMASK,ILENE)
    CALL INIRP(YDLAP,YDELAP,YDGEOMETRY%YRDIM)

    WRITE(NULOUT,FMT='('' SPECTRAL FIT'')')
    LLOPT= LNEWORO .OR. LNEWORO2
    CALL RELSPE(YDGEOMETRY,LLOPT)

    INDK=0
    DO JJ=1,NDGLG
      DO J=1,NLOENG(JJ)
        INDC=INDK+J
        IADR=J+NSTAGP(JJ)-1
        ZS(INDC,2)=XREV(IADR)
      ENDDO
      INDK=INDK+NLOENG(JJ)
    ENDDO
  ENDIF

ENDIF

!     5.2 Translation into final units.

DO J=1,ILENE
  ZS(J,1)=MAX(0.0_JPRB,SIGN(1.0_JPRB,YRCLI%SMASK-ZS(J,1)))
  ZS(J,2)=RG*ZS(J,2)
  ZS(J,3)=RG*ZS(J,3)
  ZS(J,6)=RG*ZS(J,6)
  ZS(J,7)=MAX(0.0_JPRB,MIN(1.0_JPRB,1.0_JPRB-ZS(J,7)))
  ZS(J,8)=MAX(0.0_JPRB,MIN(ZS(J,7),ZS(J,8)))
  ZS(J,9)=ZS(J,6)
ENDDO
IF (LKEYF) THEN
  DO J=1,NSEFRE
    XRNM(J)=RG*XRNM(J)
  ENDDO
ELSE
! Spectral orography written as missing data
  XRNM(:)=YRCLI%SMANQ - 1.0_JPRB
ENDIF

!     5.3 Final writing on ARPEGE files as results.

!  Field identifiers

CLPREF( 1)='SURF    '
CLPREF( 2)='SURF    '
CLPREF( 3)='SURF    '
CLPREF( 4)='SURF    '
CLPREF( 5)='SURF    '
CLPREF( 6)='SURF    '
CLPREF( 7)='SURF    '
CLPREF( 8)='SURF    '
CLPREF( 9)='SURF    '
CLPREF( 0)='SPECSURF'
CLSUFF( 1)='IND.TERREMER'
CLSUFF( 2)='GEOPOTENTIEL'
CLSUFF( 3)='ET.GEOPOTENT'
CLSUFF( 4)='VAR.GEOP.ANI'
CLSUFF( 5)='VAR.GEOP.DIR'
CLSUFF( 6)='Z0REL.FOIS.G'
CLSUFF( 7)='PROP.TERRE  '
CLSUFF( 8)='PROP.URBANIS'
CLSUFF( 9)='Z0.FOIS.G   '
CLSUFF( 0)='GEOPOTEN    '

!  Set-up for ARPEGE files

!  CADRE
CLNOMC='Const.Clim.Surfa'
LLKEEP=.TRUE.

DO J=1,NDGNH
  INLOPA(J)=NLOENG(J)
  INOZPA(J)=NMENG(J)
  ZSINLA(J)=YDCSGLEG%RMU(J)
ENDDO

!  Date
DO J=1,11
  IDATEF(J)=0
ENDDO
IDATEF(1)=1
IDATEF(2)=1
IDATEF(3)=15
IDATEF(6)=1

!  Output file
INUM=41
CLNOMF='Const.Clim'
INIV=0
IMES=1
IARP=10
LLIMST=.TRUE.
IREP=0
IARI=0

!  GRIB coding
IGRIBC=0
IGRIBO=0
IBPDG=0
IBCSP=0
ITRON=0
ILAP=0
IMOD=0
LLCOSP=.FALSE.

! Orography and roughness length are not packed at all
LLPACK(:)=.TRUE.
LLPACK( 0)=.FALSE.
LLPACK( 2)=.FALSE.
LLPACK( 6)=.FALSE.
LLPACK( 9)=.FALSE.

!  Define the CADRE

CALL FACADE(CLNOMC,NSTTYP,SLAPO,GLOPO,SLOPO,FACDI,&
 & NSMAX,NDGLG,NDLON,INLOPA,INOZPA,ZSINLA,&
 & NFLEVG,VP00,YDVAB%VALH,YDVAB%VBH,LLKEEP)  

!  Open the file and write the date

CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'NEW',.TRUE.,LLIMST,&
 & IMES,IARP,IARI,CLNOMC)  
CALL FANDAR(IREP,INUM,IDATEF)

!  Get GRIB options

CALL FAVEUR(IREP,INUM,IGRIBC,IBPDG,IBCSP,ITRON,ILAP,IMOD)
! new spectral ordering only if packing type is -1 or 3
IF (IGRIBC == -1 .OR. IGRIBC == 3) THEN
  IGRIBO=-1
ELSE
  IGRIBO=0
ENDIF

DO J=1,9
  IF (.NOT.LLPACK(J)) THEN
!   Change GRIB level of coding to "none"
    CALL FAGOTE(IREP,INUM,IGRIBO,IBPDG,IBCSP,ITRON,ILAP,IMOD)
  ENDIF
  CALL FAIENC(IREP,INUM,CLPREF(J),INIV,CLSUFF(J),ZS(1,J),LLCOSP)
! Orography : writing spectral components too
  IF (CLPREF(J) == 'SURF'.AND.CLSUFF(J) == 'GEOPOTENTIEL') THEN
    LLCOSP=.TRUE.
    SELECT CASE (IGRIBO)
    CASE (0)
      CALL FAIENC(IREP,INUM,CLPREF(0),INIV,CLSUFF(0),XRNM,LLCOSP)
    CASE (-1)
      ALLOCATE ( ZRNM (NSPEC2G) )
      ZRNM(:)=0.0_JPRB
      CALL SPREORD(YDGEOMETRY%YRDIM,YDGEOMETRY%YREDIM,YDGEOMETRY%YRELAP,1,XRNM,ZRNM,.TRUE.)
      CALL FAIENC(IREP,INUM,CLPREF(0),INIV,CLSUFF(0),ZRNM,LLCOSP)
      DEALLOCATE ( ZRNM )
    CASE DEFAULT
      CALL ABOR1('INCLI1: WRONG VALUE IGRIBO')
    END SELECT
    LLCOSP=.FALSE.
  ENDIF
! Reset previous GRIB options
  IF (.NOT.LLPACK(J)) THEN
    CALL FAGOTE(IREP,INUM,IGRIBC,IBPDG,IBCSP,ITRON,ILAP,IMOD)
  ENDIF
ENDDO

CALL LFILAF(IREP,INUM,.TRUE.)
CALL FAIRME(IREP,INUM,'KEEP')

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INCLI1',1,ZHOOK_HANDLE)
END SUBROUTINE INCLI1
