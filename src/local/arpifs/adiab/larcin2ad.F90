SUBROUTINE LARCIN2AD(YDGEOMETRY,YDML_DYN,YDSL,KPROMA,KSTART,KPROF,KFLDN,KFLDX,&
 & KSTABUF,&
 & KWISA,KXLAG,KROT,LDPLANE,&
 & KDIM,KFLDSLB1,PSLBUF1,&
 & KSTTYP,PDSTRET,PC2M1,PC2P1,P4JP,PI,PDEPI,PIS2,&
 & PLOCEN,PMUCEN,PLSDEPI,PLATI,&
 & KIBL,&
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
 & PUFZ,PVFZ)  

!**** *LARCIN2AD -semi-LAgrangian scheme:   (adjoint version)
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
!        *CALL* *LARCIN2AD(......)

!        Explicit arguments :
!        --------------------

!        INPUT or INPUT/OUTPUT matching with INPUT in the TL code:
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
!          KDIM    - number of independent updates
!          KFLDSLB1- second dimension of PSLBUF1

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
!          KIBL    - index into YRCSGEOM/YRGSGEOM instances in YDGEOMETRY
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

!        INPUT, INPUT/OUTPUT and OUTPUT, matching with OUTPUT in the TL code:
!          PSLBUF1 - set of semi-Lagrangian buffers for interpolation.
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
!      Original : 98/10/20

!     Modifications.
!     --------------
!   Modified 28-Aug-2007 by F. Vana: update of arguments for LASCAW
!   30-Jun-2008 J. Masek     Updated arguments of LASCAW.
!   K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!   K. Yessad (Aug 2009): use RIPI, RSLD
!   K. Yessad (Nov 2009): routine renaming (LAI..TLAD -> LAI..AD).
!   K. Yessad (Jan 2011): LARCHE2.. instead of LARCHE...
!   F. Vana 21-Feb-2011: update of arguments for LASCAW
!   G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!   G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
!   G. Mozdzynski (May 2012): further cleaning
!   F. Vana 13-Feb-2014: update of arguments for LASCAW
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   S. Malardel (Nov 2013): COMAD weights for SL interpolations
!   K. Yessad (March 2017): simplify level numbering in interpolator.
!   F. Vana    21-Nov-2017: Option LHOISLT
! End Modifications
!-------------------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE YOMVAR             , ONLY : LVECADIN, LSLADREP

USE EINT_MOD           , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDSLB1 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTABUF(YDSL%NDGSAH:YDSL%NDGENH) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KWISA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KXLAG(6) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KROT 
LOGICAL           ,INTENT(IN)    :: LDPLANE 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIM 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSLBUF1(YDSL%NASLB1*KFLDSLB1) 
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
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCOSCO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSINCO(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSINLA(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCOPHI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOSCO5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSINCO5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSINLA5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOPHI5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PURL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVRL0(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSL(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSL(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCSL(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSL2(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSL2(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCSL2(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PURL05(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRL05(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSL5(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSL5(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSL5(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUSL25(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSL25(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCSL25(YDSL%NASLB1,KFLDN:KFLDX) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLON(KPROMA,1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLAT(KPROMA,1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PO(KPROMA,1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQ(KPROMA,1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUF0(KPROMA,1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVF0(KPROMA,1) 
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUFZ(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFZ(KPROMA) 

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
REAL(KIND=JPRB) :: ZINC    (KDIM+1,KPROMA)
REAL(KIND=JPRB) :: ZSIGN4  (KPROMA,4 ,1,2)
REAL(KIND=JPRB) :: ZSIGN12 (KPROMA,12,1,2)

INTEGER(KIND=JPIM) :: IL0  (KPROMA,1,0:3)
INTEGER(KIND=JPIM) :: ILH0 (KPROMA,1,0:3)
INTEGER(KIND=JPIM) :: INC  (KDIM+1,KPROMA)
INTEGER(KIND=JPIM) :: IMAP4 (KPROMA,4 ,1)
INTEGER(KIND=JPIM) :: IMAP12(KPROMA,12,1)
INTEGER(KIND=JPIM) :: ISIGN(KFLDSLB1)

INTEGER(KIND=JPIM) :: ISPLTHOI, IDIMK
INTEGER(KIND=JPIM) :: IHOR, INCA, ISEP, IVLAG, IWIS, IWLAG,&
 & JINC, JROF, JFLD, IDEP(1)  

REAL(KIND=JPRB) :: ZLEV(1)
REAL(KIND=JPRB) :: ZHOOK_HANDLE

! unused arguments in call to LASCAW
! a) input - array dimensions omitted
INTEGER(KIND=JPIM) :: IUN_VAUTF(1)
REAL(KIND=JPRB) :: ZUN_SLD(1,3)
REAL(KIND=JPRB) :: ZUN_SLDW(1)
REAL(KIND=JPRB) :: ZUN_3DTW(1)
REAL(KIND=JPRB) :: ZUN_LEV(1)
REAL(KIND=JPRB) :: ZUN_LEV5(1)
REAL(KIND=JPRB) :: ZUN_VETAF(1)
REAL(KIND=JPRB) :: ZUN_VCUICO(1)
REAL(KIND=JPRB) :: ZUN_VSLD(1)
REAL(KIND=JPRB) :: ZUN_VSLDW(1)
REAL(KIND=JPRB) :: ZUN_KAPPA(KPROMA,1) ! I/O in adjoint
REAL(KIND=JPRB) :: ZUN_KAPPAT(KPROMA,1) ! I/O in adjoint
REAL(KIND=JPRB) :: ZUN_KAPPAM(1),ZUN_KAPPAH(1)
REAL(KIND=JPRB) :: ZUN_KAPPA5(1),ZUN_KAPPAT5(1)
REAL(KIND=JPRB) :: ZUN_MADU(1),ZUN_MADV(1),ZUN_MADW(1)
REAL(KIND=JPRB) :: ZUN_GAMMA_WENO(1)
REAL(KIND=JPRB) :: ZUN_DVER5(1)
! b) output - array dimensions kept for safety
INTEGER(KIND=JPIM) :: IUN_LEV(KPROMA,1)
INTEGER(KIND=JPIM) :: IUN_NOWENO(KPROMA,1)
REAL(KIND=JPRB) :: ZUN_DLAMAD(KPROMA)
REAL(KIND=JPRB) :: ZUN_CLASLD(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_CLAMAD(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_DLOMAD (KPROMA,1,0:3)
REAL(KIND=JPRB) :: ZUN_CLOSLD(KPROMA,1,3,2)
REAL(KIND=JPRB) :: ZUN_CLOMAD(KPROMA,1,3,2)
REAL(KIND=JPRB) :: ZUN_CLASLT(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_CLOSLT(KPROMA,1,3,2)
REAL(KIND=JPRB) :: ZUN_DVER(KPROMA,1)
REAL(KIND=JPRB) :: ZUN_DVERMAD(KPROMA,1)
REAL(KIND=JPRB) :: ZUN_VINTW(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_VINTWSLD(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_VINTWMAD(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_VINTWSLT(KPROMA,1,3)
REAL(KIND=JPRB) :: ZUN_VINTWS(KPROMA,1,1:4)
REAL(KIND=JPRB) :: ZUN_VDERW(KPROMA,1,0,0)  ! case KHVI=0
REAL(KIND=JPRB) :: ZUN_HVW(KPROMA,1,0)      ! case KHVI=0
REAL(KIND=JPRB) :: ZUN_CW(KPROMA,1,3)

LOGICAL         :: LLSLHD,LLSLHDQUAD,LL3DTURB
LOGICAL         :: LLCOMAD, LLCOMADH, LLCOMADV
!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "laiddi_init.intfb.h"
#include "laiddiad.intfb.h"
#include "laidli_init.intfb.h"
#include "laidliad.intfb.h"
#include "larche2.intfb.h"
#include "larche2ad.intfb.h"
#include "lascaw.intfb.h"
#include "lascawad.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LARCIN2AD',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB, &
 & YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP, &
 & YDCSGLEG=>YDGEOMETRY%YRCSGLEG,YDVSPLIP=>YDGEOMETRY%YRVSPLIP,YDVSLETA=>YDGEOMETRY%YRVSLETA, &
 & YDHSLMER=>YDGEOMETRY%YRHSLMER, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL), &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDSPGEOM=>YDGEOMETRY%YSPGEOM, YDPTRSLB1=>YDML_DYN%YRPTRSLB1)

ASSOCIATE(MSLB1SP0=>YDPTRSLB1%MSLB1SP0, MSLB1SP9=>YDPTRSLB1%MSLB1SP9, &
 & MSLB1U0=>YDPTRSLB1%MSLB1U0, MSLB1U9=>YDPTRSLB1%MSLB1U9, &
 & MSLB1UR0=>YDPTRSLB1%MSLB1UR0, MSLB1V0=>YDPTRSLB1%MSLB1V0, &
 & MSLB1V9=>YDPTRSLB1%MSLB1V9, MSLB1VR0=>YDPTRSLB1%MSLB1VR0, &
 & RPARSL1=>YDPTRSLB1%RPARSL1)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

! Switching off the SLHD/COMAD related weights computation
LLSLHD=.FALSE.
LLSLHDQUAD=.FALSE.
LLCOMAD=.FALSE.
LLCOMADH=.FALSE.
LLCOMADV=.FALSE.

ISPLTHOI=0
LL3DTURB=.FALSE.
IDIMK=1

IHOR=0
ISEP=1

! * I[X]LAG:
IVLAG=KXLAG(1)
IWLAG=KXLAG(4)

DO JFLD=1,KFLDSLB1
  IF(     RPARSL1(JFLD) ==  1.0_JPRB )THEN
    ISIGN(JFLD)=1
  ELSEIF( RPARSL1(JFLD) == -1.0_JPRB )THEN
    ISIGN(JFLD)=2
  ELSE
    CALL ABOR1("LARCIN2AD: RPARSL1 invalid sign")
  ENDIF
ENDDO

! * Input variable IWIS for LASCAW.
IF (KWISA == 1) THEN
  ! * trajectory research.
  IF (YDML_DYN%YRDYNA%LHOISLT) THEN
    IWIS=202
  ELSE
    IWIS=201
  ENDIF
ELSEIF (KWISA == 3) THEN
  ! * origin point interpolations.
  IWIS=203
ELSE
  IWIS=-999
  CALL ABOR1('LARCIN2AD: WRONG VALUE FOR IWIS')
ENDIF

! * ZERO LOCAL ARRAYS:
ZDLAT(1:KPROMA)=0.0_JPRB
ZCLA(1:KPROMA,1,1:3)=0.0_JPRB
ZDLO(1:KPROMA,1,0:3)=0.0_JPRB
ZCLO(1:KPROMA,1,1:3,1:2)=0.0_JPRB

!     ------------------------------------------------------------------

!*       2.    COMPUTATION OF LAT LON OF THE INTERPOLATION POINT.
!              IF KROT=1 COMPUTATION OF THE WIND DISPLACEMENT MATRIX
!              FROM THE INTERPOLATION POINT TO THE FINAL POINT.
!              ( T FOR LATITUDE THETA, L FOR LONGITUDE LAMBDA).
!              PO = ( 1 / (1+cos(PHI)) )
!                  *( cos(TG)*cos(T) + (1+sin(TG)*sin(T))*cos(L-LG) )
!              PQ = (-1 / (1+cos(PHI)) )
!                  *( sin(TG)+sin(T) ) * sin(L-LG)
!        (This is a repeat of the TRAJECTORY calculation)

!     ------------------------------------------------------------------

IF (LDPLANE) THEN    !! AD not yet coded !!

  ! Not coded

ELSE

  CALL LARCHE2(KPROMA,KSTART,KPROF,&
   & KSTTYP,PDSTRET,PC2M1,PC2P1,PI,PDEPI,&
   & PLOCEN,PMUCEN,YDGSGEOM,YDCSGEOM,&
   & PCOSCO5,PSINCO5,PSINLA5,PCOPHI5,&
   & KROT,PLON5,PLAT5,PO5,PQ5)  

  CALL LASCAW(YDGEOMETRY%YRVSPLIP,YDSL,KPROMA,IDIMK,KSTART,KPROF,1,&
   & KFLDN,KSTABUF,IWIS,IHOR,1,0,&
   & LLSLHD,LLSLHDQUAD,YDML_DYN%YRDYNA%LSLHD_OLD,YDML_DYN%YRDYN%LSLHDHEAT,LL3DTURB,&
   & LLCOMAD,LLCOMADH,LLCOMADV,ISPLTHOI,&
   & P4JP,PIS2,PLSDEPI,PLATI,&
   & YDHSLMER%RIPI(YDSL%NDGSAH:YDSL%NDGENH,1),&
   & YDHSLMER%RIPI(YDSL%NDGSAH:YDSL%NDGENH,2),&
   & YDHSLMER%RIPI(YDSL%NDGSAH:YDSL%NDGENH,3),&
   & ZUN_SLD(1,1),ZUN_SLD(1,2),ZUN_SLD(1,3),ZUN_SLDW,ZUN_3DTW,&
   & PLON5,PLAT5,ZUN_LEV,&
   & ZUN_VETAF,IUN_VAUTF,&
   & ZUN_VCUICO,ZUN_VSLD,ZUN_VSLDW,ZUN_GAMMA_WENO,1,1.0_JPRB,&
   & ZUN_KAPPA,ZUN_KAPPAT,ZUN_KAPPAM,ZUN_KAPPAH,&
   & ZUN_MADU,ZUN_MADV,ZUN_MADW,&
   & ZDLAT5,ZUN_DLAMAD,ZCLA5(1,1,1),ZUN_CLASLD,ZUN_CLASLT,ZUN_CLAMAD,&
   & ZDLO5(1,1,0),ZUN_DLOMAD,ZCLO5,ZUN_CLOSLD,ZUN_CLOSLT,ZUN_CLOMAD,&
   & IL0(1,1,0),ILH0(1,1,0),IUN_LEV,IUN_NOWENO,ZUN_CW,&
   & ZUN_DVER,ZUN_DVERMAD,ZUN_VINTW,ZUN_VINTWSLD,ZUN_VINTWSLT,ZUN_VINTWMAD,ZUN_VINTWS,&
   & ZUN_VDERW,ZUN_HVW)

ENDIF

IF( LSLADREP )THEN
  IF( IWIS == 201 )THEN
    CALL LAIDLI_INIT(YDML_DYN%YRSLREP,YDSL%NASLB1,KPROMA,KSTART,KPROF,1,IL0(1,1,1),IMAP4,ZSIGN4)
  ELSEIF( IWIS == 202 .OR. IWIS == 203 )THEN
    CALL LAIDLI_INIT(YDML_DYN%YRSLREP,YDSL%NASLB1,KPROMA,KSTART,KPROF,1,IL0(1,1,1),IMAP4,ZSIGN4)
    CALL LAIDDI_INIT(YDML_DYN%YRSLREP,YDSL%NASLB1,KPROMA,KSTART,KPROF,1,IL0(1,1,0),IMAP12,ZSIGN12)
  ELSE
    CALL ABOR1('LARCIN2AD: INVALID IWIS')
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       3.    MEDIUM POINT INTERPOLATIONS OF WIND
!              FOR COMPUTATION OF MEDIUM POINT.
!              -----------------------------------

IF (KWISA == 1) THEN

  IF (YDML_DYN%YRDYNA%LHOISLT) THEN
    ! Not coded
    CALL ABOR1('LARCIN2AD: LHOISLT=true is not coded for shallow water. Sorry.')
  ELSE
  !* Bilinear interpolation
  INCA=1
  CALL LAIDLIAD(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
   & KFLDX,MSLB1UR0,ISEP,IMAP4,ZSIGN4(1,1,1,ISIGN(MSLB1UR0)),&
   & LVECADIN,&
   & ZDLAT,ZDLO(1,1,1),ZDLAT5,ZDLO5(1,1,1),&
   & IL0(1,1,1),INC(INCA,1),ZINC(INCA,1),KDIM,PURL0,PURL05,PUF0,PUF05)  
  INCA=INCA+4
  CALL LAIDLIAD(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
   & KFLDX,MSLB1VR0,ISEP,IMAP4,ZSIGN4(1,1,1,ISIGN(MSLB1VR0)),&
   & LVECADIN,&
   & ZDLAT,ZDLO(1,1,1),ZDLAT5,ZDLO5(1,1,1),&
   & IL0(1,1,1),INC(INCA,1),ZINC(INCA,1),KDIM,PVRL0,PVRL05,PVF0,PVF05)  
  INCA=INCA+4
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       7.    ORIGIN POINT INTERPOLATIONS. 2D MODEL.
!              --------------------------------------

IF (KWISA == 3) THEN

  INCA=1

  IF (IWLAG == 3) THEN
    DO JROF=KSTART,KPROF
      ZUF(JROF)=PUF(JROF)
      ZVF(JROF)=PVF(JROF)
    ENDDO
    CALL LAIDLIAD(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
     & KFLDX,MSLB1U0,ISEP,IMAP4,ZSIGN4(1,1,1,ISIGN(MSLB1U0)),&
     & LVECADIN,&
     & ZDLAT,ZDLO(1,1,1),ZDLAT5,ZDLO5(1,1,1),&
     & IL0(1,1,1),INC(INCA,1),ZINC(INCA,1),KDIM,PUSL2,PUSL25,ZUF,ZUF5)  
    INCA=INCA+4
    CALL LAIDLIAD(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
     & KFLDX,MSLB1V0,ISEP,IMAP4,ZSIGN4(1,1,1,ISIGN(MSLB1V0)),&
     & LVECADIN,&
     & ZDLAT,ZDLO(1,1,1),ZDLAT5,ZDLO5(1,1,1),&
     & IL0(1,1,1),INC(INCA,1),ZINC(INCA,1),KDIM,PVSL2,PVSL25,ZVF,ZVF5)  
    INCA=INCA+4
    ! Save PUF,PVF in case of refined treatment of
    ! Coriolis term in 2TL scheme:
    DO JROF=KSTART,KPROF
      PUF(JROF)=PUF(JROF)+PUFZ(JROF)
      PVF(JROF)=PVF(JROF)+PVFZ(JROF)
    ENDDO
  ENDIF

  IF (IVLAG == 3) THEN
    DO JROF=KSTART,KPROF
      ZTF(JROF)=PCF(JROF)
    ENDDO
    CALL LAIDLIAD(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
     & KFLDX,MSLB1SP0,ISEP,IMAP4,ZSIGN4(1,1,1,ISIGN(MSLB1SP0)),&
     & LVECADIN,&
     & ZDLAT,ZDLO(1,1,1),ZDLAT5,ZDLO5(1,1,1),&
     & IL0(1,1,1),INC(INCA,1),ZINC(INCA,1),KDIM,PCSL2,PCSL25,ZTF,ZTF5)  
    INCA=INCA+4
  ENDIF

  CALL LAIDDIAD(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
   & KFLDX,MSLB1U9,ISEP,IMAP12,ZSIGN12(1,1,1,ISIGN(MSLB1U9)),&
   & IDEP,LVECADIN,LDPLANE,&
   & ZCLA,ZDLO,ZCLO,&
   & ZCLA5,ZDLO5,ZCLO5,IL0(1,1,0),&
   & INC(INCA,1),ZINC(INCA,1),KDIM,PUSL,PUSL5,PUF)  
  INCA=INCA+12
  CALL LAIDDIAD(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
   & KFLDX,MSLB1V9,ISEP,IMAP12,ZSIGN12(1,1,1,ISIGN(MSLB1V9)),&
   & IDEP,LVECADIN,LDPLANE,&
   & ZCLA,ZDLO,ZCLO,&
   & ZCLA5,ZDLO5,ZCLO5,IL0(1,1,0),&
   & INC(INCA,1),ZINC(INCA,1),KDIM,PVSL,PVSL5,PVF)  
  INCA=INCA+12
  CALL LAIDDIAD(YDSL%NASLB1,KPROMA,KSTART,KPROF,1,KFLDN,&
   & KFLDX,MSLB1SP9,ISEP,IMAP12,ZSIGN12(1,1,1,ISIGN(MSLB1SP9)),&
   & IDEP,LVECADIN,LDPLANE,&
   & ZCLA,ZDLO,ZCLO,&
   & ZCLA5,ZDLO5,ZCLO5,IL0(1,1,0),&
   & INC(INCA,1),ZINC(INCA,1),KDIM,PCSL,PCSL5,PCF)  
  INCA=INCA+12
  ! - for completeness!   (not really necessary)
  IF (IWLAG == 3) THEN
    DO JROF=KSTART,KPROF
      PUF5(JROF)=PUF5(JROF)+ZUF5(JROF)
      PVF5(JROF)=PVF5(JROF)+ZVF5(JROF)
    ENDDO
  ENDIF
  IF (IVLAG == 3) THEN
    DO JROF=KSTART,KPROF
      PCF5(JROF)=PCF5(JROF)+ZTF5(JROF)
    ENDDO
  ENDIF

ENDIF

! Finally, do all the updates together:

!!! NOTE: On no account should the following loop ordering be changed!
!!! (loop on JINC must be innermost for vectorization;
!!! outer loop must be on JROF for reproducibility)

INCA=INCA-1
DO JROF=KSTART,KPROF
#if defined(VPP)
!OCL NOVREC
#elif defined(NECSX)
!CDIR VWORK=STACK
!CDIR VWORKSZ=4M
!CDIR UNROLL=4
!CDIR NODEP
!CDIR GTHREORDER
!CDIR VOVERTAKE(PSLBUF1),VOB
#endif
  DO JINC=1,INCA
    PSLBUF1(INC(JINC,JROF))=PSLBUF1(INC(JINC,JROF))&
     & +ZINC(JINC,JROF)  
  ENDDO
ENDDO

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

IF (LDPLANE) THEN                  !! AD not done yet !!

  ! Not coded

ELSE

  CALL LASCAWAD(YDSL,KPROMA,KSTART,KPROF,1,&
   & KFLDN,IWIS,IHOR,1,&
   & LLSLHD,LLSLHDQUAD,YDML_DYN%YRDYNA%LSLHD_OLD,YDML_DYN%YRDYN%LSLHDHEAT,&
   & P4JP,PIS2,PLSDEPI,PLATI,&
   & YDHSLMER%RIPI(YDSL%NDGSAH:YDSL%NDGENH,1:3),ZUN_SLD,ZUN_SLDW,ZUN_GAMMA_WENO,&
   & PLON,PLAT,ZLEV,ZUN_KAPPA,ZUN_KAPPAT,&
   & PLON5,PLAT5,ZUN_LEV5,ZUN_KAPPA5,ZUN_KAPPAT5,ZUN_DVER5,ZUN_VETAF,&
   & IUN_VAUTF,ZUN_VCUICO,ZUN_VSLD,ZUN_VSLDW,1,1.0_JPRB,&
   & ZDLAT,ZCLA,ZUN_CLASLD,ZUN_CLASLT,ZDLO,ZCLO,ZUN_CLOSLD,ZUN_CLOSLT,&
   & ZUN_CW,ZUN_DVER,ZUN_VINTW,ZUN_VINTWSLD,ZUN_VINTWSLT)

  CALL LARCHE2AD(KPROMA,KSTART,KPROF,&
   & KSTTYP,PDSTRET,PC2M1,PC2P1,PI,PDEPI,&
   & PLOCEN,PMUCEN,YDGSGEOM,YDCSGEOM,&
   & PCOSCO,PSINCO,PSINLA,PCOPHI,&
   & PCOSCO5,PSINCO5,PSINLA5,PCOPHI5,&
   & KROT,PLON,PLAT,PO,PQ)

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LARCIN2AD',1,ZHOOK_HANDLE)
END SUBROUTINE LARCIN2AD
