SUBROUTINE LARCIN2TL(YDHSLMER,YDDYN,YDSL,KPROMA,KSTART,KPROF,KFLDN,KFLDX,&
 & KSTABUF,&
 & KWISA,KXLAG,KROT,LDPLANE,&
 & LDQMHW,LDQMHP,&
 & KSTTYP,PDSTRET,PC2M1,PC2P1,P4JP,PI,PDEPI,PIS2,&
 & PLOCEN,PMUCEN,PLSDEPI,PLATI,&
 & YDGSGEOM,YDCSGEOM,&
 & PCOSCO,PSINCO,PSINLA,PCOPHI,&
 & PCOSCO5,PSINCO5,PSINLA5,PCOPHI5,&
 & PURL0,PVRL0,&
 & PUSL,PVSL,PCSL,PUSL2,PVSL2,PCSL2,&
 & PURL05,PVRL05,&
 & PUSL5,PVSL5,PCSL5,PUSL25,PVSL25,PCSL25,&
 & PLON,PLAT,&
 & PO,PQ,&
 & PUF0,PVF0,&
 & PUF,PVF,PCF,&
 & PLON5,PLAT5,&
 & PO5,PQ5,&
 & PUF05,PVF05,&
 & PUF5,PVF5,PCF5,&
 & PUFZ,PVFZ,PUFZ5,PVFZ5)  

!**** *LARCIN2TL -semi-LAgrangian scheme:   (tangent-linear version)
!                 Research of the Coordinates (of the medium or origin
!                 point) and INterpolations (2-d model).

!     Purpose.
!     --------
!       Computes the longitude and latitude of the interpolation
!       point from its cartesian coordinates.
!       Then computes the vector displacement matrix
!                        I po pq I
!                        I       I
!                        I-pq po I
!       from the interpolation point to the grid point.
!       At last determines the interpolation grid:
!       - computation of the latitude and the longitude of the
!         point situated at the upper left corner of the 16 points square,
!         and of the interpolation point.
!       - interpolations.

!**   Interface.
!     ----------
!        *CALL* *LARCINTL(......)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          YDSL    - SL_STRUCT definition
!          KPROMA  - horizontal dimension.
!          KSTART  - first element of arrays where
!                    computations are performed.
!          KPROF   - depth of work.
!          KFLDN   - number of the first field.
!          KFLDX   - number of the last field.
!          KSTABUF - for a latitude IGL, KSTABUF(IGL) is the
!                    address of the element corresponding to
!                    (ILON=1,IGL) in the KPROMA arrays.
!          KWISA   - switch :
!                    KWISA=1: interpolations for wind to find trajectory.
!                    KWISA=3: origin point interp.
!          KXLAG   - array containing switches concerning discretisations
!                    of equations (NVLAG,NTLAG,NWLAG,NSVDLAG,NSPDLAG).
!          KROT    - KROT=1: computation of the elements po and pq
!                    of the wind displacement matrix.
!                    KROT=0: no computation.
!          LDPLANE - switch: .T. = plane geometry; .F. = spherical geometry.
!          LDQMH[X]- Use quasi-monotone interpolation in the horizontal
!                    for field X.

!          KSTTYP  - 1: Not tilted pole;  2: Tilted pole.
!          PDSTRET - 2*c (where c is the stretching factor).
!          PC2M1   - c*c-1.
!          PC2P1   - c*c+1.
!          P4JP    - Approximative inverse of the differences of latitudes.
!          PI      - number PI
!          PDEPI   - 2 * PI
!          PIS2    - PI / 2
!          PLOCEN  - geographical longitude of the stretching pole.
!          PMUCEN  - sine of the geographical latitude of the stretching pole.
!          PLSDEPI - (Number of points by latitude) / (2 * PI) .
!          PLATI   - latitude on the computational sphere.
!          YDGSGEOM- grid point geometry.
!          YDCSGEOM- computational sphere geometry.
!          PCOSCO  - cos(Longitude-Longitude(grid-point))*cos(Latitude)
!                    of the interpolation point (geographical longitude
!                    and latitude).
!                    ! if LDPLANE: x - coordinate (fractional system).
!          PSINCO  - sin(Longitude-Longitude(grid-point))*cos(Latitude)
!                    of the interpolation point (geographical longitude
!                    and latitude).
!                    ! if LDPLANE: y - coordinate (fractional system).
!          PSINLA  - sine of the interpolation point geographical latitude.
!          PCOPHI  - cosine of the geographical angle between the
!                    interpolation point and the grid-point.
! ------------------------ trajectory variables ---------------------------
!          PCOSCO5 - cos(Longitude-Longitude(grid-point))*cos(Latitude)
!                    of the interpolation point (geographical longitude
!                    and latitude).
!                    ! if LDPLANE: x - coordinate (fractional system).
!          PSINCO5 - sin(Longitude-Longitude(grid-point))*cos(Latitude)
!                    of the interpolation point (geographical longitude
!                    and latitude).
!                    ! if LDPLANE: y - coordinate (fractional system).
!          PSINLA5 - sine of the interpolation point geographical latitude.
!          PCOPHI5 - cosine of the geographical angle between the
!                    interpolation point and the grid-point.
! -------------------------------------------------------------------------
!          PURL0   - U-component of the wind.
!          PVRL0   - V-component of the wind.
!          PUSL    - Quantity to interpolate at O in the U-wind equation.
!          PVSL    - Quantity to interpolate at O in the V-wind equation.
!          PCSL    - Quantity to interpolate at O in the continuity equation.
!          PUSL2   - Second quantity to interpolate at O
!                    (if NWLAG=3) in the U-component of the wind equation.
!          PVSL2   - Second quantity to interpolate at O
!                    (if NWLAG=3) in the V-component of the wind equation.
!          PCSL2   - Second quantity to interpolate at O
!                    (if NVLAG=3) in the continuity equation.
! ------------------------ trajectory variables ---------------------------
!          PURL05  - U-component of the wind.
!          PVRL05  - V-component of the wind.
!          PUSL5   - Quantity to interpolate at O in the U-wind equation.
!          PVSL5   - Quantity to interpolate at O in the V-wind equation.
!          PCSL5   - Quantity to interpolate at O in the continuity equation.
!          PUSL25  - Second quantity to interpolate at O
!                    (if NWLAG=3) in the U-component of the wind equation.
!          PVSL25  - Second quantity to interpolate at O
!                    (if NWLAG=3) in the V-component of the wind equation.
!          PCSL25  - Second quantity to interpolate at O
!                    (if NVLAG=3) in the continuity equation.
! -------------------------------------------------------------------------

!        OUTPUT:
!          PLON    - computational sphere longitude of interpolation point.
!          PLAT    - computational sphere latitude of interpolation point.
!          PO      - first element of the wind displacement matrix.
!          PQ      - second element of the wind displacement matrix.
!          PUF0    - Interpolated U-wind at the medium point.
!          PVF0    - Interpolated V-wind at the medium point.
!          PUF     - Interpolated quantity at O in the U-wind equation.
!          PVF     - Interpolated quantity at O in the V-wind equation.
!          PCF     - Interpolated quantity at O in the continuity equation.
! ------------------------ trajectory variables ---------------------------
!          PLON5   - computational sphere longitude of interpolation point.
!          PLAT5   - computational sphere latitude of interpolation point.
!          PO5     - first element of the wind displacement matrix.
!          PQ5     - second element of the wind displacement matrix.
!          PUF05   - Interpolated U-wind at the medium point.
!          PVF05   - Interpolated V-wind at the medium point.
!          PUF5    - Interpolated quantity at O in the U-wind equation.
!          PVF5    - Interpolated quantity at O in the V-wind equation.
!          PCF5    - Interpolated quantity at O in the continuity equation.
! -------------------------------------------------------------------------
!          PUFZ    - Interpolated quantity at O in the U-wind equation
!                    in the SL2TL scheme for case L2TLFF=.T.
!          PVFZ    - Interpolated quantity at O in the V-wind equation
!                    in the SL2TL scheme for case L2TLFF=.T.
! ------------------------ trajectory variables ---------------------------
!          PUFZ5   - Interpolated quantity at O in the U-wind equation
!                    in the SL2TL scheme for case L2TLFF=.T.
!          PVFZ5   - Interpolated quantity at O in the V-wind equation
!                    in the SL2TL scheme for case L2TLFF=.T.
! -------------------------------------------------------------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!      C. Temperton (ECMWF)
!      Original : 98/06/15

!     Modifications.
!     --------------
!      Modified 01-08-30 by K. YESSAD: pruning and some other cleanings.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Modified 30-Jun-2008 by F. Vana: new interpolation data-flow
!      K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!      K. Yessad (Aug 2009): use RIPI, RSLD
!      K. Yessad (Jan 2011): call LARCHE2TL instead of LARCHETL.
!      G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
!      G. Mozdzynski (May 2012): further cleaning
!      F. Vana 13-Feb-2014  Update of arguments for LASCAWTL
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      F. Vana    21-Nov-2017: Option LHOISLT
!     ------------------------------------------------------------------

USE YOMDYN    , ONLY : TDYN
USE YOMHSLMER , ONLY : THSLMER
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK
USE EINT_MOD  , ONLY : SL_STRUCT
USE YOMDYNA   , ONLY : YRDYNA
USE YOMCSGEOM , ONLY : TCSGEOM
USE YOMGSGEOM , ONLY : TGSGEOM

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(THSLMER)     ,INTENT(IN)    :: YDHSLMER
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTABUF(YDSL%NDGSAH:YDSL%NDGENH) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KWISA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KXLAG(6) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KROT 
LOGICAL           ,INTENT(IN)    :: LDPLANE 
LOGICAL           ,INTENT(IN)    :: LDQMHW 
LOGICAL           ,INTENT(IN)    :: LDQMHP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTTYP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDSTRET
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC2M1 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC2P1
REAL(KIND=JPRB)   ,INTENT(IN)    :: P4JP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDEPI 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLOCEN
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMUCEN
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLATI(YDSL%NDGSAH:YDSL%NDGENH) 
TYPE(TGSGEOM)     ,INTENT(IN)    :: YDGSGEOM
TYPE(TCSGEOM)     ,INTENT(IN)    :: YDCSGEOM
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOSCO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSINCO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSINLA(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOPHI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOSCO5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSINCO5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSINLA5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOPHI5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PURL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSL(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSL(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSL(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSL2(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSL2(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSL2(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PURL05(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRL05(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSL5(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSL5(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSL5(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSL25(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSL25(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSL25(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLON(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLAT(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQ(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUF0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVF0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUF(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVF(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCF(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLON5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLAT5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PO5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQ5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUF05(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVF05(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUF5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVF5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCF5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUFZ(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVFZ(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUFZ5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVFZ5(KPROMA) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZDLAT(KPROMA)
REAL(KIND=JPRB) :: ZCLA (KPROMA,1,3)
REAL(KIND=JPRB) :: ZDLO (KPROMA,1,0:3)
REAL(KIND=JPRB) :: ZCLO (KPROMA,1,3,2)
REAL(KIND=JPRB) :: ZDLAT5(KPROMA)
REAL(KIND=JPRB) :: ZCLA5 (KPROMA,1,3)
REAL(KIND=JPRB) :: ZDLO5 (KPROMA,1,0:3)
REAL(KIND=JPRB) :: ZCLO5 (KPROMA,1,3,2)
REAL(KIND=JPRB) :: ZUF     (KPROMA)
REAL(KIND=JPRB) :: ZVF     (KPROMA)
REAL(KIND=JPRB) :: ZTF     (KPROMA)
REAL(KIND=JPRB) :: ZUF5    (KPROMA)
REAL(KIND=JPRB) :: ZVF5    (KPROMA)
REAL(KIND=JPRB) :: ZTF5    (KPROMA)

INTEGER(KIND=JPIM) :: IL0 (KPROMA,1,0:3)
INTEGER(KIND=JPIM) :: ILH0(KPROMA,1,0:3)

INTEGER(KIND=JPIM) :: IHOR, IVLAG, IWIS, IWLAG, JROF

! unused arguments in call to LASCAW
! a) input - array dimensions omitted
INTEGER(KIND=JPIM) :: IUN_VAUTF(1)
REAL(KIND=JPRB) :: ZUN_SLD(1,3)
REAL(KIND=JPRB) :: ZUN_SLDW(1)
REAL(KIND=JPRB) :: ZUN_LEV(1)
REAL(KIND=JPRB) :: ZUN_VETAF(1)
REAL(KIND=JPRB) :: ZUN_VCUICO(1)
REAL(KIND=JPRB) :: ZUN_VSLD(1)
REAL(KIND=JPRB) :: ZUN_VSLDW(1)
REAL(KIND=JPRB) :: ZUN_KAPPA(1),ZUN_KAPPAT(1)
REAL(KIND=JPRB) :: ZUN_GAMMA_WENO(1)
! b) output - array dimensions kept for safety
INTEGER(KIND=JPIM) :: IUN_LEV(KPROMA,1)
INTEGER(KIND=JPIM) :: IUN_NOWENO(KPROMA,1)
REAL(KIND=JPRB) :: ZUN_CLASLD(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_CLOSLD(KPROMA,1,3,2)
REAL(KIND=JPRB) :: ZUN_CLASLT(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_CLOSLT(KPROMA,1,3,2)
REAL(KIND=JPRB) :: ZUN_DVER(KPROMA,1)
REAL(KIND=JPRB) :: ZUN_VINTW(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_VINTWSLD(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_VINTWSLT(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_CLASLD5(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_CLOSLD5(KPROMA,1,3,2)
REAL(KIND=JPRB) :: ZUN_CLASLT5(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_CLOSLT5(KPROMA,1,3,2)
REAL(KIND=JPRB) :: ZUN_DVER5(KPROMA,1)
REAL(KIND=JPRB) :: ZUN_VINTW5(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_VINTWSLD5(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_VINTWSLT5(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_CW(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_CW5(KPROMA,1,3)

LOGICAL         :: LLSLHD,LLSLHDQUAD,LLSLHD_OLD

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "laidditl.intfb.h"
#include "laidlitl.intfb.h"
#include "larche2tl.intfb.h"
#include "lascawtl.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LARCIN2TL',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

! Switching off the SLHD related weights computation (not coded for 2D model)
LLSLHD=.FALSE.
LLSLHDQUAD=.FALSE.
LLSLHD_OLD=.FALSE.

IHOR=0

! * I[X]LAG:
IVLAG=KXLAG(1)
IWLAG=KXLAG(4)

! * Input variable IWIS for LASCAW.
IF (KWISA == 1) THEN
  ! * trajectory research.
  IF (YRDYNA%LHOISLT) THEN
    IWIS=202
  ELSE
    IWIS=201
  ENDIF
ELSEIF (KWISA == 3) THEN
  ! * origin point interpolations.
  IWIS=203
ELSE
  IWIS=-999
  CALL ABOR1('LARCIN2TL: WRONG VALUE FOR IWIS')
ENDIF
!     ------------------------------------------------------------------

!*       2.    COMPUTATION OF LAT LON OF THE INTERPOLATION POINT.
!              IF KROT=1 COMPUTATION OF THE WIND DISPLACEMENT MATRIX
!              FROM THE INTERPOLATION POINT TO THE FINAL POINT.
!              ( T FOR LATITUDE THETA, L FOR LONGITUDE LAMBDA).
!              PO = ( 1 / (1+cos(PHI)) )
!                  *( cos(TG)*cos(T) + (1+sin(TG)*sin(T))*cos(L-LG) )
!              PQ = (-1 / (1+cos(PHI)) )
!                  *( sin(TG)+sin(T) ) * sin(L-LG)

!     ------------------------------------------------------------------

IF (LDPLANE) THEN                  !! TL not done yet !!

  ! Not coded

ELSE

  CALL LARCHE2TL(KPROMA,KSTART,KPROF,&
   & KSTTYP,PDSTRET,PC2M1,PC2P1,PI,PDEPI,&
   & PLOCEN,PMUCEN,YDGSGEOM,YDCSGEOM,&
   & PCOSCO,PSINCO,PSINLA,PCOPHI,&
   & PCOSCO5,PSINCO5,PSINLA5,PCOPHI5,&
   & KROT,PLON,PLAT,PO,PQ,&
   & PLON5,PLAT5,PO5,PQ5)  

  CALL LASCAWTL(YDSL,KPROMA,KSTART,KPROF,1,&
   & KFLDN,KSTABUF,IWIS,IHOR,1,&
   & LLSLHD,LLSLHDQUAD,LLSLHD_OLD,YDDYN%LSLHDHEAT, &
   & P4JP,PIS2,PLSDEPI,PLATI,&
   & YDHSLMER%RIPI(YDSL%NDGSAH:YDSL%NDGENH,1:3),ZUN_SLD,ZUN_SLDW,&
   & PLON,PLAT,ZUN_LEV,PLON5,PLAT5,ZUN_LEV,&
   & ZUN_VETAF,IUN_VAUTF,&
   & ZUN_VCUICO,ZUN_VSLD,ZUN_VSLDW,ZUN_GAMMA_WENO,1,1.0_JPRB,ZUN_KAPPA,ZUN_KAPPA,ZUN_KAPPAT,ZUN_KAPPAT,&
   & ZDLAT,ZCLA,ZUN_CLASLD,ZUN_CLASLT,ZDLO,ZCLO,ZUN_CLOSLD,ZUN_CLOSLT,IL0(1,1,0),ILH0(1,1,0),IUN_LEV,&
   & IUN_NOWENO,ZUN_CW,ZUN_DVER,ZUN_VINTW,ZUN_VINTWSLD,ZUN_VINTWSLT,&
   & ZDLAT5,ZCLA5,ZUN_CLASLD5,ZUN_CLASLT5,ZDLO5,ZCLO5,ZUN_CLOSLD5,ZUN_CLOSLT5,&
   & ZUN_CW5,ZUN_DVER5,ZUN_VINTW5,ZUN_VINTWSLD5,ZUN_VINTWSLT5)

ENDIF

!     ------------------------------------------------------------------

!*       3.    MEDIUM POINT INTERPOLATIONS OF WIND
!              FOR COMPUTATION OF MEDIUM POINT.
!              -----------------------------------

IF (KWISA == 1) THEN

  IF (YRDYNA%LHOISLT) THEN
    ! Not coded
    call ABOR1('LARCIN2TL: LHOISLT=true is not coded for shallow water. Sorry.')
  ELSE
  ! * Bilinear interpolations.
  CALL LAIDLITL(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
   & KFLDX,ZDLAT,ZDLO(1,1,1),ZDLAT5,ZDLO5(1,1,1),&
   & IL0(1,1,1),PURL0,PURL05,PUF0,PUF05)  
  CALL LAIDLITL(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
   & KFLDX,ZDLAT,ZDLO(1,1,1),ZDLAT5,ZDLO5(1,1,1),&
   & IL0(1,1,1),PVRL0,PVRL05,PVF0,PVF05)  
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       7.    ORIGIN POINT INTERPOLATIONS.
!              ----------------------------

IF (KWISA == 3) THEN

  IF (LDQMHW.OR.LDQMHP) THEN
    CALL ABOR1('LARCIN2TL: QUASI-MONOTONE OPTIONS NOT ALLOWED')
  ENDIF

  CALL LAIDDITL(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
   & KFLDX,ZCLA,ZDLO,ZCLO,ZCLA5,ZDLO5,ZCLO5,IL0(1,1,0),&
   & PUSL,PUSL5,PUF,PUF5)  
  CALL LAIDDITL(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
   & KFLDX,ZCLA,ZDLO,ZCLO,ZCLA5,ZDLO5,ZCLO5,IL0(1,1,0),&
   & PVSL,PVSL5,PVF,PVF5)  
  CALL LAIDDITL(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
   & KFLDX,ZCLA,ZDLO,ZCLO,ZCLA5,ZDLO5,ZCLO5,IL0(1,1,0),&
   & PCSL,PCSL5,PCF,PCF5)  

  IF (IWLAG == 3) THEN
!   Save PUF,PVF in case of refined treatment of
!   Coriolis term in 2TL scheme:
    DO JROF=KSTART,KPROF
      PUFZ(JROF)=PUF(JROF)
      PVFZ(JROF)=PVF(JROF)
      PUFZ5(JROF)=PUF5(JROF)
      PVFZ5(JROF)=PVF5(JROF)
    ENDDO
    CALL LAIDLITL(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
     & KFLDX,ZDLAT,ZDLO(1,1,1),ZDLAT5,ZDLO5(1,1,1),&
     & IL0(1,1,1),PUSL2,PUSL25,ZUF,ZUF5)  
    CALL LAIDLITL(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
     & KFLDX,ZDLAT,ZDLO(1,1,1),ZDLAT5,ZDLO5(1,1,1),&
     & IL0(1,1,1),PVSL2,PVSL25,ZVF,ZVF5)  
    DO JROF=KSTART,KPROF
      PUF(JROF)=PUF(JROF)+ZUF(JROF)
      PVF(JROF)=PVF(JROF)+ZVF(JROF)
      PUF5(JROF)=PUF5(JROF)+ZUF5(JROF)
      PVF5(JROF)=PVF5(JROF)+ZVF5(JROF)
    ENDDO
  ENDIF

  IF (IVLAG == 3) THEN
    CALL LAIDLITL(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
     & KFLDX,ZDLAT,ZDLO(1,1,1),ZDLAT5,ZDLO5(1,1,1),&
     & IL0(1,1,1),PCSL2,PCSL25,ZTF,ZTF5)  
    DO JROF=KSTART,KPROF
      PCF(JROF)=PCF(JROF)+ZTF(JROF)
      PCF5(JROF)=PCF5(JROF)+ZTF5(JROF)
    ENDDO
  ENDIF

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LARCIN2TL',1,ZHOOK_HANDLE)
END SUBROUTINE LARCIN2TL

