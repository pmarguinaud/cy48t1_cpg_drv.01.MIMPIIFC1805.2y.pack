SUBROUTINE LADINETL(YDGEOMETRY,YDRIP,YDDYN,KSTART,KPROF,YDSL,KPROMA,KFLDN,KFLDX,&
 & KIBL,&
 & PURL0,PVRL0,PURL,PVRL,&
 & PUL0,PVL0,PCL0,PUL9,PVL9,PCL9,&
 & PURL05,PVRL05,PURL5,PVRL5,&
 & PUL05,PVL05,PCL05,PUL95,PVL95,PCL95,&
 & PUSI,PVSI,PUSI5,PVSI5,PUT1,PVT1,PSPT1,PUT15,PVT15 )  

!**** *LADINETL* - semi-LAgrangian scheme:   (tangent-linear version)
!                Interface subroutine for interpolations and lagrangian
!                trends. (Programme INterface d'Ensemble).
!                2D Shallow water.

!     Purpose.
!     --------
!           Grid point calculations in dynamics.

!           Computes the medium point M and origin point O.
!           Does interpolations at O.
!           Computes lagrangian trends.

!**   Interface.
!     ----------
!        *CALL* *LADINETL(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KSTART  - first element of arrays where
!                    computations are performed.
!          KPROF   - depth of work.
!          YDSL    - SL_STRUCT definition.
!          KPROMA  - horizontal dimension.
!          KFLDN   - number of the first field for semi-lag variables.
!          KFLDX   - number of the last field for semi-lag variables.
!          KIBL    - index into YRGSGEOM/YRCSGEOM instances in YDGEOMETRY
!          PURL0   - U-wind, to be interpolated.
!          PVRL0   - V-wind, to be interpolated.
!          PURL    - U-wind for research of medium point.
!          PVRL    - V-wind for research of medium point.
!          PUL0    - Quantity to interpolate linearly at O in the
!                    U-wind equation if NWLAG=3.
!          PVL0    - Quantity to interpolate linearly at O in the
!                    V-wind equation if NWLAG=3.
!          PCL0    - Quantity to interpolate linearly at O in the
!                    continuity equation if NVLAG=3.
!          PUL9    - Quantity to interpolate at O in the U-wind equation.
!          PVL9    - Quantity to interpolate at O in the V-wind equation.
!          PCL9    - Quantity to interpolate at O in the continuity equation.
!          PUSI    - Semi-implicit term in the U-wind equation  ) needed for 
!          PVSI    - Semi-implicit term in the V-wind equation  ) L2TLFF
! ---------------------------- trajectory variables ----------------------
!          PURL05  - U-component of the wind, to be interpolated.
!          PVRL05  - V-component of the wind, to be interpolated.
!          PURL5   - U-component of the wind for research of medium point.
!          PVRL5   - V-component of the wind for research of medium point.
!          PUL05   - Quantity to interpolate linearly at O in the
!                    U-wind equation if NWLAG=3.
!          PVL05   - Quantity to interpolate linearly at O in the
!                    V-wind equation if NWLAG=3.
!          PCL05   - Quantity to interpolate linearly at O in the
!                    continuity equation if NVLAG=3.
!          PUL95   - Quantity to interpolate at O in the U-wind equation.
!          PVL95   - Quantity to interpolate at O in the V-wind equation.
!          PCL95   - Quantity to interpolate at O in the continuity equation.
!          PUSI5   - Semi-implicit term in the U-equation  ) needed for
!          PVSI5   - Semi-implicit term in the V-equation  ) L2TLFF
!          PUT15   - t+dt U-wind                           )
!          PVT15   - t+dt V-wind                           )
! ------------------------------------------------------------------------

!        INPUT/OUTPUT:
!          PUT1    - t+dt U-wind.
!          PVT1    - t+dt V-wind.
!          PSPT1   - t+dt equivalent height.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE documentation vol 2 ch 1 and vol 3 ch 6

!     Author.
!     -------
!      C. Temperton (ECMWF)
!      Original : 98/06/15

!     Modifications.
!     --------------
!   Modified 01-08-30 by K. YESSAD: pruning and some other cleanings.
!   R. El Khatib : 01-08-07 Pruning options
!   M.Hamrud      01-Oct-2003 CY28 Cleaning
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!   K. Yessad (Jan 2011): call LARCHE2TL instead of LARCHETL.
!   G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!   G. Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM, TCSGEOM and TCSGLEG
!   G. Mozdzynski (May 2012): further cleaning
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCT0       , ONLY : LRPLANE
USE YOMDYNA      , ONLY : YRDYNA
USE YOMDYN       , ONLY : TDYN
USE YOMCST       , ONLY : RPI, ROMEGA, RA
USE YOMRIP       , ONLY : TRIP
USE EINT_MOD     , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(IN)    :: PURL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PURL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUL9(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVL9(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCL9(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PURL05(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRL05(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PURL5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRL5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUL05(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVL05(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCL05(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUL95(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVL95(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCL95(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSI5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSI5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUT1(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVT1(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPT1(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT15(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVT15(KPROMA) 

!     ------------------------------------------------------------------

! - trigonometric arrays.
REAL(KIND=JPRB) :: ZSINCO(KPROMA)
REAL(KIND=JPRB) :: ZCOSCO(KPROMA)
REAL(KIND=JPRB) :: ZSINLA(KPROMA)
REAL(KIND=JPRB) :: ZCOPHI(KPROMA)
REAL(KIND=JPRB) :: ZP    (KPROMA)
REAL(KIND=JPRB) :: ZQ    (KPROMA)
REAL(KIND=JPRB) :: ZR    (KPROMA)
REAL(KIND=JPRB) :: ZS    (KPROMA)
REAL(KIND=JPRB) :: ZSINCO5(KPROMA)
REAL(KIND=JPRB) :: ZCOSCO5(KPROMA)
REAL(KIND=JPRB) :: ZSINLA5(KPROMA)
REAL(KIND=JPRB) :: ZCOPHI5(KPROMA)
REAL(KIND=JPRB) :: ZP5    (KPROMA)
REAL(KIND=JPRB) :: ZQ5    (KPROMA)
REAL(KIND=JPRB) :: ZR5    (KPROMA)
REAL(KIND=JPRB) :: ZS5    (KPROMA)
REAL(KIND=JPRB) :: ZLON   (KPROMA)
REAL(KIND=JPRB) :: ZLAT   (KPROMA)
REAL(KIND=JPRB) :: ZLON5  (KPROMA)
REAL(KIND=JPRB) :: ZLAT5  (KPROMA)
REAL(KIND=JPRB) :: ZLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)
! - interpolated quantities arrays.
REAL(KIND=JPRB) :: ZU9(KPROMA)
REAL(KIND=JPRB) :: ZV9(KPROMA)
REAL(KIND=JPRB) :: ZC9(KPROMA)
REAL(KIND=JPRB) :: ZZUZ9(KPROMA)         ! To allow YRDYN%L2TLFF option
REAL(KIND=JPRB) :: ZZVZ9(KPROMA)         ! To allow YRDYN%L2TLFF option
REAL(KIND=JPRB) :: ZU95(KPROMA)
REAL(KIND=JPRB) :: ZV95(KPROMA)
REAL(KIND=JPRB) :: ZC95(KPROMA)
REAL(KIND=JPRB) :: ZZUZ95(KPROMA)         ! To allow YRDYN%L2TLFF option
REAL(KIND=JPRB) :: ZZVZ95(KPROMA)         ! To allow YRDYN%L2TLFF option
REAL(KIND=JPRB) :: ZUT15(KPROMA)
REAL(KIND=JPRB) :: ZVT15(KPROMA)

INTEGER(KIND=JPIM) :: ISTALAT(YDSL%NDGSAH:YDSL%NDGENH)

INTEGER(KIND=JPIM) :: IGLGLO, JGL, JROF

LOGICAL :: LLSETTLST

REAL(KIND=JPRB) :: ZCOPH, ZDEPI, ZDSTRET, ZPHI, ZPIS2, ZPP,&
 & ZPP5, ZQQ, ZQQ5, ZSINX, ZZUT1, ZZVT1,&
 & ZPP9, ZQQ9, ZPP1, ZQQ1, ZPP95, ZQQ95,&
 & ZZUT15, ZZVT15, ZPHI5, ZCOPH5, ZSINX5  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "laconetl.intfb.h"
#include "lainor2tl.intfb.h"
#include "larche2tl.intfb.h"
#include "larmes2tl.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LADINETL',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG,YDVSPLIP=>YDGEOMETRY%YRVSPLIP,YDVSLETA=>YDGEOMETRY%YRVSLETA, &
 & YDHSLMER=>YDGEOMETRY%YRHSLMER, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL),YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(L2TLFF=>YDDYN%L2TLFF, &
 & LADVF=>YDDYN%LADVF, LQMHP=>YDDYN%LQMHP, LQMHW=>YDDYN%LQMHW, &
 & NITMP=>YDDYN%NITMP, NVLAG=>YDDYN%NVLAG, NWLAG=>YDDYN%NWLAG, &
 & NSTTYP=>YDGEOMETRY%YRGEM%NSTTYP, &
 & R4JP=>YDGEOMETRY%YRGEM%R4JP, RC2M1=>YDGEOMETRY%YRGEM%RC2M1, RC2P1=>YDGEOMETRY%YRGEM%RC2P1, &
 & RLOCEN=>YDGEOMETRY%YRGEM%RLOCEN, RMUCEN=>YDGEOMETRY%YRGEM%RMUCEN, RSTRET=>YDGEOMETRY%YRGEM%RSTRET, &
 & RDTS22=>YDRIP%RDTS22, &
 & RDTS62=>YDRIP%RDTS62, RDTSA=>YDRIP%RDTSA, RTDT=>YDRIP%RTDT, &
 & YDGEM=>YDGEOMETRY%YRGEM)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS:
!              ----------------------------

ZDEPI  =2.0_JPRB*RPI
ZPIS2  =0.5_JPRB*RPI

ZDSTRET=2.0_JPRB*RSTRET
DO JGL=MAX(YDSL%NDGSAG,YDSL%NDGSAL+YDSL%NFRSTLOFF-YDSL%NSLWIDE)-YDSL%NFRSTLOFF,&
   & MIN(YDSL%NDGENG,YDSL%NDGENL+YDSL%NFRSTLOFF+YDSL%NSLWIDE)-YDSL%NFRSTLOFF  
  IGLGLO=JGL+YDSL%NFRSTLOFF
  ZLSDEPI(JGL)=REAL(YDSL%NLOENG(IGLGLO),JPRB)/ZDEPI
ENDDO

DO JGL=MAX(YDSL%NDGSAH,LBOUND(YDSL%NSLOFF,1)),MIN(YDSL%NDGENH,UBOUND(YDSL%NSLOFF,1))
  ISTALAT(JGL)=YDSL%NSLOFF(JGL)
ENDDO

!     ------------------------------------------------------------------

!*       2.    COMPUTATION OF THE MEDIUM POINT AND
!              INTERPOLATIONS AT THIS POINT.
!              ------------------------------------

LLSETTLST=YRDYNA%LSETTLST.AND.(.NOT.YRDYNA%LELTRA)
CALL LARMES2TL(YDGEOMETRY%YRHSLMER,YDDYN,YDSL,KPROMA,KSTART,KPROF,KFLDN,KFLDX,&
 & ISTALAT,&
 & NVLAG,NWLAG,NITMP,&
 & LRPLANE,LLSETTLST,YRDYNA%LELTRA,&
 & PURL0,PVRL0,PURL,PVRL,&
 & PURL05,PVRL05,PURL5,PVRL5,&
 & NSTTYP,ZDSTRET,RC2M1,RC2P1,R4JP,RPI,ZDEPI,ZPIS2,&
 & RDTSA,RDTS62,RDTS22,&
 & RLOCEN,RMUCEN,ZLSDEPI,YDCSGLEG%RLATI(YDSL%NDGSAH:),&
 & YDGSGEOM,YDCSGEOM,&
 & ZLON,ZLAT,&
 & ZCOSCO,ZSINCO,ZSINLA,ZCOPHI,&
 & ZR,ZS,&
 & ZLON5,ZLAT5,&
 & ZCOSCO5,ZSINCO5,ZSINLA5,ZCOPHI5,&
 & ZR5,ZS5)

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF THE ORIGIN POINT AND
!              INTERPOLATIONS AT THIS POINT.
!              ------------------------------------

CALL LAINOR2TL(YDGEOMETRY%YRHSLMER,YDDYN,YDSL,KPROMA,KSTART,KPROF,KFLDN,KFLDX,&
 & ISTALAT,&
 & NVLAG,NWLAG,&
 & LRPLANE,LQMHW,LQMHP,&
 & PUL9,PVL9,PCL9,PUL0,PVL0,PCL0,&
 & PUL95,PVL95,PCL95,PUL05,PVL05,PCL05,&
 & YDGSGEOM,YDCSGEOM,&
 & NSTTYP,ZDSTRET,RC2M1,RC2P1,R4JP,RPI,ZDEPI,ZPIS2,&
 & RLOCEN,RMUCEN,ZLSDEPI,YDCSGLEG%RLATI(YDSL%NDGSAH:),&
 & ZCOSCO,ZSINCO,ZSINLA,ZCOPHI,&
 & ZCOSCO5,ZSINCO5,ZSINLA5,ZCOPHI5,&
 & ZLON,ZLAT,&
 & ZP,ZQ,&
 & ZU9,ZV9,ZC9,&
 & ZLON5,ZLAT5,&
 & ZP5,ZQ5,&
 & ZU95,ZV95,ZC95,&
 & ZZUZ9,ZZVZ9,ZZUZ95,ZZVZ95)

!     ------------------------------------------------------------------

!*       4.    COMPUTATION OF LAGRANGIAN TRENDS.
!              ---------------------------------

!*       4.1   MOMENTUM EQUATION.

DO JROF=KSTART,KPROF
  ZPP5=ZP5(JROF)*YDGSGEOM%GNORDM(JROF)-ZQ5(JROF)*YDGSGEOM%GNORDL(JROF)
  ZQQ5=ZQ5(JROF)*YDGSGEOM%GNORDM(JROF)+ZP5(JROF)*YDGSGEOM%GNORDL(JROF)
  ZPP=ZP(JROF)*YDGSGEOM%GNORDM(JROF)-ZQ(JROF)*YDGSGEOM%GNORDL(JROF)
  ZQQ=ZQ(JROF)*YDGSGEOM%GNORDM(JROF)+ZP(JROF)*YDGSGEOM%GNORDL(JROF)
  ZUT15(JROF) = PUT15(JROF) + ZPP5*ZU95(JROF) + ZQQ5*ZV95(JROF)
  ZVT15(JROF) = PVT15(JROF) + ZPP5*ZV95(JROF) - ZQQ5*ZU95(JROF)
  PUT1(JROF) = PUT1(JROF) + ZPP5*ZU9(JROF) + ZQQ5*ZV9(JROF)&
   & + ZPP*ZU95(JROF) + ZQQ*ZV95(JROF)  
  PVT1(JROF) = PVT1(JROF) - ZQQ5*ZU9(JROF) + ZPP5*ZV9(JROF)&
   & - ZQQ*ZU95(JROF) + ZPP*ZV95(JROF)  
ENDDO

IF (LADVF) THEN     !! LRPLANE version not done !!

  CALL LACONETL(KPROMA,KSTART,KPROF,1,NSTTYP,&
   & ROMEGA,ZDSTRET,RC2M1,RC2P1,RMUCEN,&
   & RA,RPI,ZLON5,ZLAT5,ZU95,ZV95,ZLON,ZLAT,ZU9,ZV9)  

  DO JROF=KSTART,KPROF
    ZPP5=ZP5(JROF)*YDGSGEOM%GNORDM(JROF)-ZQ5(JROF)*YDGSGEOM%GNORDL(JROF)
    ZQQ5=ZQ5(JROF)*YDGSGEOM%GNORDM(JROF)+ZP5(JROF)*YDGSGEOM%GNORDL(JROF)
    ZPP=ZP(JROF)*YDGSGEOM%GNORDM(JROF)-ZQ(JROF)*YDGSGEOM%GNORDL(JROF)
    ZQQ=ZQ(JROF)*YDGSGEOM%GNORDM(JROF)+ZP(JROF)*YDGSGEOM%GNORDL(JROF)
    ZUT15(JROF) = ZUT15(JROF) + ZPP5*ZU95(JROF) + ZQQ5*ZV95(JROF)
    ZVT15(JROF) = ZVT15(JROF) - ZQQ5*ZU95(JROF) + ZPP5*ZV95(JROF)
    PUT1(JROF) = PUT1(JROF) + ZPP5*ZU9(JROF) + ZQQ5*ZV9(JROF)&
     & + ZPP*ZU95(JROF) + ZQQ*ZV95(JROF)  
    PVT1(JROF) = PVT1(JROF) - ZQQ5*ZU9(JROF) + ZPP5*ZV9(JROF)&
     & - ZQQ*ZU95(JROF) + ZPP*ZV95(JROF)  
  ENDDO

ENDIF

IF (L2TLFF) THEN  !! LRPLANE version not done !!

! * compute new origin point:
  DO JROF=KSTART,KPROF
!   * old rotation and vector deplacement:
    ZPP5=ZP5(JROF)*YDGSGEOM%GNORDM(JROF)-ZQ5(JROF)*YDGSGEOM%GNORDL(JROF)
    ZQQ5=ZQ5(JROF)*YDGSGEOM%GNORDM(JROF)+ZP5(JROF)*YDGSGEOM%GNORDL(JROF)
    ZPP95=ZP5(JROF)
    ZQQ95=ZQ5(JROF)
    ZPP=ZP(JROF)*YDGSGEOM%GNORDM(JROF)-ZQ(JROF)*YDGSGEOM%GNORDL(JROF)
    ZQQ=ZQ(JROF)*YDGSGEOM%GNORDM(JROF)+ZP(JROF)*YDGSGEOM%GNORDL(JROF)
    ZPP9=ZP(JROF)
    ZQQ9=ZQ(JROF)
    ZPP1=YDGSGEOM%GNORDM(JROF)
    ZQQ1=-YDGSGEOM%GNORDL(JROF)
!   * new estimate of trajectory wind:
    ZZUT15=0.5_JPRB*(&
     & ZPP1*(ZUT15(JROF)-PUSI5(JROF))+ZQQ1*(ZVT15(JROF)-PVSI5(JROF))&
     & +ZPP95*ZZUZ95(JROF)+ZQQ95*ZZVZ95(JROF))  
    ZZVT15=0.5_JPRB*(&
     & -ZQQ1*(ZUT15(JROF)-PUSI5(JROF))+ZPP1*(ZVT15(JROF)-PVSI5(JROF))&
     & -ZQQ95*ZZUZ95(JROF)+ZPP95*ZZVZ95(JROF))  
    ZZUT1=0.5_JPRB*(&
     & ZPP1*(PUT1(JROF)-PUSI(JROF))+ZQQ1*(PVT1(JROF)-PVSI(JROF))&
     & +ZPP9*ZZUZ95(JROF)+ZQQ9*ZZVZ95(JROF)&
     & +ZPP95*ZZUZ9(JROF)+ZQQ95*ZZVZ9(JROF))  
    ZZVT1=0.5_JPRB*(&
     & -ZQQ1*(PUT1(JROF)-PUSI(JROF))+ZPP1*(PVT1(JROF)-PVSI(JROF))&
     & -ZQQ9*ZZUZ95(JROF)+ZPP9*ZZVZ95(JROF)&
     & -ZQQ95*ZZUZ9(JROF)+ZPP95*ZZVZ9(JROF))  
!   * subtract old estimate of advected Coriolis term:
    PUT1(JROF) = PUT1(JROF) - ZPP5*ZU9(JROF) - ZQQ5*ZV9(JROF)&
     & - ZPP*ZU95(JROF) - ZQQ*ZV95(JROF)  
    PVT1(JROF) = PVT1(JROF) + ZQQ5*ZU9(JROF) - ZPP5*ZV9(JROF)&
     & + ZQQ*ZU95(JROF) - ZPP*ZV95(JROF)  
!   * compute new origin point:
    ZPHI5=0.5_JPRB*RTDT*SQRT(ZZUT15*ZZUT15+ZZVT15*ZZVT15)/RA
    ZCOPH5=1.0_JPRB-0.5_JPRB*ZPHI5*ZPHI5
    ZSINX5=ZCOPH5*RTDT*(1.0_JPRB-ZPHI5*ZPHI5/6.0_JPRB)/RA
    ZSINLA5(JROF)=YDGSGEOM%GEMU(JROF)*(2.0_JPRB*ZCOPH5*ZCOPH5-1.0_JPRB)&
     & -ZZVT15*YDGSGEOM%GSQM2(JROF)*ZSINX5  
    ZCOSCO5(JROF)=YDGSGEOM%GSQM2(JROF)*(2.0_JPRB*ZCOPH5*ZCOPH5-1.0_JPRB)&
     & +ZZVT15*YDGSGEOM%GEMU(JROF)*ZSINX5  
    ZSINCO5(JROF)=-ZZUT15*ZSINX5
    ZCOPHI5(JROF)=2.0_JPRB*ZCOPH5*ZCOPH5-1.0_JPRB
    ZPHI=0.5_JPRB*RTDT*(ZZUT1*ZZUT15+ZZVT1*ZZVT15)&
     & /(RA*SQRT(ZZUT15*ZZUT15+ZZVT15*ZZVT15))  
    ZCOPH=-ZPHI*ZPHI5
    ZSINX=RTDT*(ZCOPH*(1.0_JPRB-ZPHI5*ZPHI5/6.0_JPRB)&
     & -ZCOPH5*ZPHI*ZPHI5/3.0_JPRB)/RA  
    ZSINLA(JROF)=YDGSGEOM%GEMU(JROF)*4.0_JPRB*ZCOPH*ZCOPH5&
     & -YDGSGEOM%GSQM2(JROF)*(ZZVT1*ZSINX5+ZZVT15*ZSINX)  
    ZCOSCO(JROF)=YDGSGEOM%GSQM2(JROF)*4.0_JPRB*ZCOPH*ZCOPH5&
     & +YDGSGEOM%GEMU(JROF)*(ZZVT1*ZSINX5+ZZVT15*ZSINX)  
    ZSINCO(JROF)=-ZZUT1*ZSINX5-ZZUT15*ZSINX
    ZCOPHI(JROF)=4.0_JPRB*ZCOPH*ZCOPH5
  ENDDO
! * compute rotation matrix and Coriolis term at the new origin point.
  CALL LARCHE2TL(KPROMA,KSTART,KPROF,NSTTYP,ZDSTRET,RC2M1,RC2P1,&
   & RPI,ZDEPI,RLOCEN,RMUCEN,YDGSGEOM,YDCSGEOM,&
   & ZCOSCO,ZSINCO,ZSINLA,ZCOPHI,ZCOSCO5,ZSINCO5,&
   & ZSINLA5,ZCOPHI5,1,ZLON,ZLAT,ZP,ZQ,ZLON5,ZLAT5,ZP5,ZQ5)  
  CALL LACONETL(KPROMA,KSTART,KPROF,1,NSTTYP,&
   & ROMEGA,ZDSTRET,RC2M1,RC2P1,RMUCEN,&
   & RA,RPI,ZLON5,ZLAT5,ZU95,ZV95,ZLON,ZLAT,ZU9,ZV9)  
! * add new estimate of advected Coriolis term:
  DO JROF=KSTART,KPROF
    ZPP5=ZP5(JROF)*YDGSGEOM%GNORDM(JROF)-ZQ5(JROF)*YDGSGEOM%GNORDL(JROF)
    ZQQ5=ZQ5(JROF)*YDGSGEOM%GNORDM(JROF)+ZP5(JROF)*YDGSGEOM%GNORDL(JROF)
    ZPP=ZP(JROF)*YDGSGEOM%GNORDM(JROF)-ZQ(JROF)*YDGSGEOM%GNORDL(JROF)
    ZQQ=ZQ(JROF)*YDGSGEOM%GNORDM(JROF)+ZP(JROF)*YDGSGEOM%GNORDL(JROF)
    PUT1(JROF) = PUT1(JROF) + ZPP5*ZU9(JROF) + ZQQ5*ZV9(JROF)&
     & + ZPP*ZU95(JROF) + ZQQ*ZV95(JROF)  
    PVT1(JROF) = PVT1(JROF) - ZQQ5*ZU9(JROF) + ZPP5*ZV9(JROF)&
     & - ZQQ*ZU95(JROF) + ZPP*ZV95(JROF)  
  ENDDO

ENDIF

!*       4.2   CONTINUITY EQUATION.

DO JROF=KSTART,KPROF
  PSPT1(JROF)=PSPT1(JROF)+ZC9(JROF)
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LADINETL',1,ZHOOK_HANDLE)
END SUBROUTINE LADINETL
