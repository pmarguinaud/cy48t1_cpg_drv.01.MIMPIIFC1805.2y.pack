SUBROUTINE LASCAWTL(YDSL,KPROMB,KST,KPROF,KFLEV,&
 & KFLDN,KSTABUF,KWIS,KHOR,KWENO,&
 & LDSLHD,LDSLHDQUAD,LDSLHD_OLD,LDSLHDHEAT,&
 & P4JP,PIS2,PLSDEPI,PLATI,&
 & PIPI,PSLD,PSLDW,PLON,PLAT,PLEV,&
 & PLON5,PLAT5,PLEV5,&
 & PVETA,KVAUT,&
 & PVCUICO,PVSLD,PVSLDW,PGAMMA_WENO,KRLEVX,PVRLEVX,PKAPPA,PKAPPA5,PKAPPAT,PKAPPAT5,&
 & PDLAT,PCLA,PCLASLD,PCLASLT,PDLO,PCLO,PCLOSLD,PCLOSLT,KL0,KLH0,KLEV,&
 & KNOWENO,PCW,PDVER,PVINTW,PVINTWSLD,PVINTWSLT,&
 & PDLAT5,PCLA5,PCLASLD5,PCLASLT5,PDLO5,PCLO5,PCLOSLD5,PCLOSLT5,&
 & PCW5,PDVER5,PVINTW5,PVINTWSLD5,PVINTWSLT5)  

!     ------------------------------------------------------------------

!**** *LASCAWTL - Externalisable interpolator:   (tangent-linear version)
!                 Storage of Coordinates And Weights.
!                 Spherical geometry version (Gaussian grids)

!     Purpose.
!     --------
!       Determines the interpolation grid:
!       - computation of the coordinates of the
!         point situated at the upper left corner of the 16 points
!         square, and of the interpolation point.
!       - computation of weights.
!       Storage of coordinates and weights.

!**   Interface.
!     ----------
!        *CALL* *LASCAWTL(.........)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          YDSL    - SL_STRUCT definition
!          KPROMB  - horizontal dimension for interpolation point quantities.
!          KST     - first element of arrays where computations are performed.
!          KPROF   - depth of work.
!          KFLEV   - vertical dimension.
!          KFLDN   - number of the first field.
!          KSTABUF - for a latitude IGL, KSTABUF(IGL) is the
!                    address of the element corresponding to
!                    (ILON=0,IGL) in the KPROMB arrays.
!          KWIS    - kind of interpolation.
!          KHOR    - 0: Horizontal interpolation for each level
!                       of a 3D variable.
!                    1: Interpolation for all origin points corresponding
!                       to a final point of a 2D variable.
!          KWENO   - 1/3 = off/on WENO: extra dimension used for WENO
!          LDSLHD  - key activating SLHD weights precomputation.
!          LDSLHDQUAD - key activating quadratic weights precomputation.
!          LDSLHD_OLD - use old SLHD interpolator
!          LDSLHDHEAT - If true, the triggering function for heat variables differs from the one for momentum variables
!          P4JP    - Approximative inverse of the differences of latitudes.
!          PIS2    - PI / 2
!          PLSDEPI - (Number of points by latitude) / (2 * PI) .
!          PLATI   - latitude on the computational sphere.
!          PIPI    - coefficients for the bicubic interpolations.
!          PSLD    - auxiliary quantity for SLHD interpolation in latitude.
!          PSLDW   - weights for SLHD Laplacian smoother in latitude.
!          PLON    - Interpolation point longitude on the work sphere.
!          PLAT    - Interpolation point latitude on the work sphere.
!          PLEV    - vertical coordinate of the interpolation point.
! ---------------------- Trajectory variables ------------------------
!          PLON5   - Interpolation point longitude on the work sphere.
!          PLAT5   - Interpolation point latitude on the work sphere.
!          PLEV5   - vertical coordinate of the interpolation point.
! --------------------------------------------------------------------
!          PVETA   - Values of ETA.
!          KVAUT   - Help table for vertical box search: gives the number
!                    of the layer immediately above eta.
!          PVCUICO - Denominators of the vertical cubic interpolation
!                    coefficients
!          PVSLD   - auxiliary quantities for vertical SLHD interpolation
!          PVSLDW  - weights for SLHD vertical Laplacian smoother
!          PGAMMA_WENO - weights for vertical WENO interpolation
!          KRLEVX  - Dimension of KVAUT
!          PVRLEVX - REAL(KRLEVX).
!          PKAPPA  - kappa function ("coefficient of SLHD") based on the
!                    rescaled horizontal deformation of the flow evaluated
!                    at instant "t" for the final point F
!          PKAPPAT - kappa function for heat variables
! ---------------------- trajectory variable -------------------------
!          PKAPPA5 - kappa function ("coefficient of SLHD") based on the
!                    rescaled horizontal deformation of the flow evaluated
!                    at instant "t" for the final point F
!          PKAPPAT5 - kappa function for heat variables
! --------------------------------------------------------------------

!        OUTPUT:
!          PDLAT     - distance for horizontal linear interpolations
!                      in latitude.
!          PCLA      - weights for horizontal cubic interpolations
!                      in latitude.
!          PCLASLD   - weights for horizontal cubic interpolations
!                      in latitude, SLHD case
!          PCLASLT   - weights for horizontal cubic interpolations
!                      in latitude, SLHD case on T
!          PDLO      - distances for horizontal linear interpolations
!                      in longitude (latitude rows 0, 1, 2, 3)
!          PCLO      - weights for horizontal cubic interpolations in
!                      longitude (latitude rows 1, 2)
!          PCLOSLD   - weights for horizontal cubic interpolations in
!                      longitude, SLHD case (latitude rows 1, 2)
!          PCLOSLT   - weights for horizontal cubic interpolations in
!                      longitude, SLHD case on T (latitude rows 1, 2)
!          KL0       - index of the four western points
!                      of the 16 points interpolation grid.
!          KLH0      - second value of index of the four western points
!                      of the 16 points interpolation grid if needed.
!          KLEV      - lower level of the vertical interpolation
!                      grid needed for vertical interpolations.
!          KNOWENO   - special boundary treatment for WENO
!          PCW       - C_k weights for the vertical WENO interpolation
!          PDVER     - distance for vertical linear interpolation
!          PVINTW    - vertical cubic interpolation weights
!          PVINTWSLD - vertical cubic interpolation weights, SLHD case
!          PVINTWSLT - vertical cubic interpolation weights, SLHD case on T

!  ------------------- Trajectory variables ------------------------
!          PDLAT5    - distance for horizontal linear interpolations
!                      in latitude.
!          PCLA5     - weights for horizontal cubic interpolations
!                      in latitude.
!          PCLASLD5  - weights for horizontal cubic interpolations
!                      in latitude, SLHD case
!          PCLASLT5  - weights for horizontal cubic interpolations
!                      in latitude, SLHD case on T
!          PDLO5     - distances for horizontal linear interpolations
!                      in longitude (latitude rows 0, 1, 2, 3)
!          PCLO5     - weights for horizontal cubic interpolations in
!                      longitude (latitude rows 1, 2)
!          PCLOSLD5  - weights for horizontal cubic interpolations in
!                      longitude, SLHD case (latitude rows 1, 2)
!          PCLOSLT5  - weights for horizontal cubic interpolations in
!                      longitude, SLHD case (latitude rows 1, 2) on T
!          PCW5      - C_k weights for the vertical WENO interpolation
!          PDVER5    - distance for vertical linear interpolation
!          PVINTW5   - vertical cubic interpolation weights
!          PVINTWSLD5- vertical cubic interpolation weights, SLHD case
!          PVINTWSLT5- vertical cubic interpolation weights, SLHD case on T
! --------------------------------------------------------------------

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
!      C. Temperton, ECMWF
!      Original : 98/06/12.

!     Modifications.
!     --------------
!      Modified 26-May-2008 : F. Vana - weights driven interpolation fully
!                              projected to TL (still excluding P/C and half
!                              level interpolation) + vectorization support
!      K. Yessad (Dec 2008): remove useless dummy arguments
!      K. Yessad (Feb 2009): split loops, rewrite in a shorter way.
!      K. Yessad (Aug 2009): use RIPI, RSLD
!      G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      G.Mozdzynski (Aug 2011): support higher order interpolation
!      F. Vana 13-feb-2014 SLHD weights for heat variables
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      K. Yessad (March 2017): simplify level numbering.
!      F. Vana    21-Nov-2017: Option LHOISLT
!      F. Vana    25-Feb-2019: quintic vertical interpolation
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK

! arp/ifs dependencies to be solved later.
USE YOMMP0   , ONLY : LOPT_SCALAR, NPROC
USE YOMCT0   , ONLY : LREGETA
USE YOMDYNA  , ONLY : SLHDKMIN, HOISLTH, HOISLTV, LHOISLT, LSLTVWENO

USE EINT_MOD , ONLY : SL_STRUCT,JPDUP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(SL_STRUCT),   INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMB 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KRLEVX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTABUF(YDSL%NDGSAH:YDSL%NDGENH) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KWIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KHOR 
INTEGER(KIND=JPIM),INTENT(IN)    :: KWENO
LOGICAL           ,INTENT(IN)    :: LDSLHD
LOGICAL           ,INTENT(IN)    :: LDSLHDQUAD
LOGICAL           ,INTENT(IN)    :: LDSLHD_OLD
LOGICAL           ,INTENT(IN)    :: LDSLHDHEAT
REAL(KIND=JPRB)   ,INTENT(IN)    :: P4JP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLATI(YDSL%NDGSAH:YDSL%NDGENH) 
REAL(KIND=JPRD)   ,INTENT(IN)    :: PIPI(YDSL%NDGSAH:YDSL%NDGENH,3) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLD(YDSL%NDGSAH:YDSL%NDGENH,3) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLDW(3,3,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLON(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLEV(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLON5(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT5(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLEV5(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVETA(0:KFLEV+1) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVAUT(0:KRLEVX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCUICO(4,0:KFLEV-1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSLD(3,0:KFLEV-1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSLDW(3,3,0:KFLEV-1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAMMA_WENO(KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRLEVX 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPPA(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPPA5(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPPAT(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPPAT5(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLAT(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLA(KPROMB,KFLEV,3) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLASLD(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLASLT(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLO(KPROMB,KFLEV,0:3) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLO(KPROMB,KFLEV,3,2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLOSLD(KPROMB,KFLEV,3,2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLOSLT(KPROMB,KFLEV,3,2)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KL0(KPROMB,KFLEV,0:3) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLH0(KPROMB,KFLEV,0:3) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLEV(KPROMB,KFLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KNOWENO(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCW(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDVER(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTW(KPROMB,KFLEV,3*KWENO) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTWSLD(KPROMB,KFLEV,3) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTWSLT(KPROMB,KFLEV,3) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLAT5(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLA5(KPROMB,KFLEV,3) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLASLD5(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLASLT5(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLO5(KPROMB,KFLEV,0:3) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLO5(KPROMB,KFLEV,3,2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLOSLD5(KPROMB,KFLEV,3,2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLOSLT5(KPROMB,KFLEV,3,2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCW5(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDVER5(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTW5(KPROMB,KFLEV,3*KWENO) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTWSLD5(KPROMB,KFLEV,3) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTWSLT5(KPROMB,KFLEV,3) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IADDR(YDSL%NDGSAH:YDSL%NDGENH)

INTEGER(KIND=JPIM) :: IFLVM2, ILA, ILEV, ILEVV, IBCLIM, &
 & ILO, ILO1, ILO2, ILO3, IZLAT, JLAT, JLEV, JROF, IJ_, J_, JJ
INTEGER(KIND=JPIM) :: IILA(KPROMB,KFLEV)
INTEGER(KIND=JPIM) :: IWLEV(KPROMB,KFLEV)

REAL(KIND=JPRB) :: ZDA, ZDB, ZDC, ZDD, ZFAC, ZLO, ZLO1, ZLO2, ZLO3, ZEPS

LOGICAL :: LLT_SLHD(3),LLSLHD,LLSLHDQUAD,LLSLHD_OLD

! Duplicata of some dummy or local arrays for optimisation on NEC platform.
INTEGER(KIND=JPIM) :: IADDR_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB) :: ZVCUICO_(JPDUP,4,0:KFLEV-1)
REAL(KIND=JPRB) :: ZVSLD_(JPDUP,3,0:KFLEV-1)
REAL(KIND=JPRB) :: ZVSLDW_(JPDUP,3,3,0:KFLEV-1)
REAL(KIND=JPRB) :: ZLSDEPI_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB) :: ZRLATI_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB) :: ZRIPI0_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB) :: ZRIPI1_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB) :: ZRIPI2_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB) :: ZSLD1_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB) :: ZSLD2_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB) :: ZSLD3_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB) :: ZSLDW_(JPDUP,3,3,YDSL%NDGSAH:YDSL%NDGENH)

REAL(KIND=JPRB) :: ZZWH(KPROMB,3,KFLEV)
REAL(KIND=JPRB) :: ZZWH5(KPROMB,3,KFLEV)

REAL(KIND=JPRB) :: ZSLHDKMINH,ZSLHDKMINV
REAL(KIND=JPRB), PARAMETER :: ZSLHDKMINV_WENO=0._JPRB   ! WENO only runs with Lagrangian cubic!!!

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "lascaw_cla_tl.intfb.h"
#include "lascaw_clo_tl.intfb.h"
#include "lascaw_vintw_tl.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LASCAWTL',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       0.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

ZEPS=100.0_JPRB*TINY(1.0_JPRB)

! cases relevant for SLHD scheme (switches LDSLHD, LDSLHDQUAD are
! deactivated during computation of medium points in LAPINEA)
LLSLHD=LDSLHD.AND.(KWIS==103.OR.KWIS==106.OR.KWIS==203)
LLSLHDQUAD=(LDSLHDQUAD.AND.(KWIS==103.OR.KWIS==106.OR.KWIS==203)).OR. &
   & ((KWIS==102.OR.KWIS==202).AND.LHOISLT.AND.((HOISLTV/=0_JPRB).OR.(HOISLTH/=0_JPRB)))

! switch for old SLHD scheme
LLSLHD_OLD=LLSLHD.AND.LDSLHD_OLD

LLT_SLHD(1)=LLSLHD
LLT_SLHD(2)=LLSLHDQUAD
LLT_SLHD(3)=LLSLHD_OLD

ZSLHDKMINH=SLHDKMIN
ZSLHDKMINV=SLHDKMIN

! Modify previous defaults for high order trajectory research
IF (KWIS==102.OR.KWIS==202) THEN
  ZSLHDKMINH=HOISLTH
  ZSLHDKMINV=HOISLTV
ENDIF

DO JLAT=YDSL%NDGSAH,YDSL%NDGENH
  IADDR(JLAT)=KSTABUF(JLAT)+YDSL%NASLB1*(0-KFLDN)
ENDDO

IFLVM2=KFLEV-2

DO J_ = 1, JPDUP
  IADDR_(J_,YDSL%NDGSAH:YDSL%NDGENH)=IADDR(YDSL%NDGSAH:YDSL%NDGENH)
ENDDO
DO J_ = 1, JPDUP
  ZLSDEPI_(J_,YDSL%NDGSAH:YDSL%NDGENH)=PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)
  ZRLATI_(J_,YDSL%NDGSAH:YDSL%NDGENH)=PLATI(YDSL%NDGSAH:YDSL%NDGENH)
  ZRIPI0_(J_,YDSL%NDGSAH:YDSL%NDGENH)=PIPI(YDSL%NDGSAH:YDSL%NDGENH,1)
  ZRIPI1_(J_,YDSL%NDGSAH:YDSL%NDGENH)=PIPI(YDSL%NDGSAH:YDSL%NDGENH,2)
  ZRIPI2_(J_,YDSL%NDGSAH:YDSL%NDGENH)=PIPI(YDSL%NDGSAH:YDSL%NDGENH,3)
  IF (LOPT_SCALAR) THEN
    ZVCUICO_(J_,1:4,0:KFLEV-1)=PVCUICO(1:4,0:KFLEV-1)
  ELSE
!CDIR NOUNROLL
    DO JLEV=1,4*KFLEV
      ZVCUICO_(J_,JLEV,0)=PVCUICO(JLEV,0)
    ENDDO
  ENDIF
ENDDO
IF (LLSLHDQUAD) THEN
  DO J_ = 1, JPDUP
    ZSLD1_(J_,YDSL%NDGSAH:YDSL%NDGENH)=PSLD(YDSL%NDGSAH:YDSL%NDGENH,1)
    ZSLD2_(J_,YDSL%NDGSAH:YDSL%NDGENH)=PSLD(YDSL%NDGSAH:YDSL%NDGENH,2)
    ZSLD3_(J_,YDSL%NDGSAH:YDSL%NDGENH)=PSLD(YDSL%NDGSAH:YDSL%NDGENH,3)
    IF (LOPT_SCALAR) THEN
      ZVSLD_(J_,1:3,0:KFLEV-1)=PVSLD(1:3,0:KFLEV-1)
    ELSE
!CDIR NOUNROLL
      DO JLEV=1,3*KFLEV
        ZVSLD_(J_,JLEV,0)=PVSLD(JLEV,0)
      ENDDO
    ENDIF
  ENDDO
ENDIF
IF ( LLSLHD ) THEN
  DO J_ = 1, JPDUP
    ZSLDW_(J_,1:3,1:3,YDSL%NDGSAH:YDSL%NDGENH)=PSLDW(1:3,1:3,YDSL%NDGSAH:YDSL%NDGENH)
    IF (LOPT_SCALAR) THEN
      ZVSLDW_(J_,1:3,1:3,0:KFLEV-1)=PVSLDW(1:3,1:3,0:KFLEV-1)
    ELSE
!CDIR NOUNROLL
      DO JLEV=1,9*KFLEV
        ZVSLDW_(J_,JLEV,1,0)=PVSLDW(JLEV,1,0)
      ENDDO
    ENDIF
  ENDDO
ENDIF

! Update the LSLHDQUAD for horizontal interpolation
IF (KWIS == 102 .OR. KWIS == 202) LLT_SLHD(2)= (HOISLTH /= 0._JPRB)

!     ------------------------------------------------------------------

!*       1.    3D MODEL.
!              ---------

!        1.01  Coordinates and weights for trilinear interpolations.

IF (KWIS == 101) THEN

  ZFAC=PVRLEVX/(PVETA(KFLEV+1)-PVETA(0))

  DO JLEV=1,KFLEV

    ! * Calculation of linear weights, KL0.
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1

      IZLAT =INT(P4JP*(PIS2-PLAT5(JROF,JLEV))+0.75_JPRB+ZEPS)-YDSL%NFRSTLOFF
      ILA   =IZLAT+NINT(SIGN(0.5_JPRB,ZRLATI_(IJ_,IZLAT)-PLAT5(JROF,JLEV))-1.5_JPRB)
      PDLAT5(JROF,JLEV)=(PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+1))&
       & /(ZRLATI_(IJ_,ILA+2)-ZRLATI_(IJ_,ILA+1)+ZEPS)  
      PDLAT(JROF,JLEV)=PLAT(JROF,JLEV)/(ZRLATI_(IJ_,ILA+2)-ZRLATI_(IJ_,ILA+1)+ZEPS)

      ZLO1  =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+1)
      ILO1  =INT(ZLO1)
      PDLO5(JROF,JLEV,1)=ZLO1-REAL(ILO1,JPRB)
      PDLO(JROF,JLEV,1)=PLON(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+1)
      ZLO2  =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+2)
      ILO2  =INT(ZLO2)
      PDLO5(JROF,JLEV,2)=ZLO2-REAL(ILO2,JPRB)
      PDLO(JROF,JLEV,2)=PLON(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+2)

      ILEV  =KVAUT(INT(PLEV5(JROF,JLEV)*ZFAC))-1
      IF(ILEV < IFLVM2.AND.&
       & (PLEV5(JROF,JLEV)-PVETA(ILEV+2)) > 0.0_JPRB) ILEV=ILEV+1  

      KLEV(JROF,JLEV)=ILEV

      PDVER5(JROF,JLEV)=(PLEV5(JROF,JLEV)-PVETA(ILEV+1))/&
       & (PVETA(ILEV+2)-PVETA(ILEV+1))  
      PDVER(JROF,JLEV)=PLEV(JROF,JLEV)/(PVETA(ILEV+2)-PVETA(ILEV+1))

      KL0(JROF,JLEV,1)=IADDR_(IJ_,ILA+1)+YDSL%NSLEXT(ILO1,ILA+1)
      KL0(JROF,JLEV,2)=IADDR_(IJ_,ILA+2)+YDSL%NSLEXT(ILO2,ILA+2)
    ENDDO

    ! * Mask calculation for on-demand communications:
    IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
!CDIR NODEP
      DO JROF=KST,KPROF
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)  )=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)+1)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)+2)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)  )=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)+1)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)+2)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)+3)=1
      ENDDO
    ENDIF

    DO JROF=KST,KPROF
      KL0(JROF,JLEV,1)=KL0(JROF,JLEV,1)+YDSL%NASLB1*KLEV(JROF,JLEV)
      KL0(JROF,JLEV,2)=KL0(JROF,JLEV,2)+YDSL%NASLB1*KLEV(JROF,JLEV)
    ENDDO

  ENDDO

ENDIF

!        1.03  Coordinates and weights for ( horizontal 12 points
!              + vertical cubic + 32 points interpolations ) or
!              ( horizontal 12 points + 32 points interpolations ).

IF (KWIS == 102 .OR. KWIS == 103 .OR. KWIS == 104 .OR. KWIS == 105 .OR. KWIS == 106) THEN

  ZFAC=PVRLEVX/(PVETA(KFLEV+1)-PVETA(0))

  DO JLEV=1,KFLEV
    IF (KHOR == 0) ILEV=JLEV
    IF (KHOR == 1) ILEV=KFLDN

    ! * Calculation of linear weights, KL0, KLH0.
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1

      IZLAT =INT(P4JP*(PIS2-PLAT5(JROF,JLEV))+0.75_JPRB+ZEPS)-YDSL%NFRSTLOFF
      ILA   =IZLAT+NINT(SIGN(0.5_JPRB,ZRLATI_(IJ_,IZLAT)-PLAT5(JROF,JLEV))-1.5_JPRB)
      PDLAT5(JROF,JLEV)=(PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+1))&
       & /(ZRLATI_(IJ_,ILA+2)-ZRLATI_(IJ_,ILA+1)+ZEPS)  
      PDLAT(JROF,JLEV)=PLAT(JROF,JLEV)/(ZRLATI_(IJ_,ILA+2)-ZRLATI_(IJ_,ILA+1)+ZEPS)

      ZDA   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA)
      ZDB   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+1)
      ZDC   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+2)
      ZDD   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+3)

      ! for input to LASCAW_CLA_TL.
      IILA(JROF,JLEV)=ILA
      ZZWH5(JROF,1,JLEV)=(ZDA*ZDC)*ZDD*ZRIPI0_(IJ_,ILA+1)
      ZZWH5(JROF,2,JLEV)=(ZDA*ZDB)*ZDD*ZRIPI1_(IJ_,ILA+1)
      ZZWH5(JROF,3,JLEV)=(ZDA*ZDB)*ZDC*ZRIPI2_(IJ_,ILA+1)
      ZZWH(JROF,1,JLEV)=(ZDA*(ZDC+ZDD)+ZDC*ZDD)*PLAT(JROF,JLEV)*ZRIPI0_(IJ_,ILA+1)
      ZZWH(JROF,2,JLEV)=(ZDA*(ZDB+ZDD)+ZDB*ZDD)*PLAT(JROF,JLEV)*ZRIPI1_(IJ_,ILA+1)
      ZZWH(JROF,3,JLEV)=(ZDA*(ZDB+ZDC)+ZDB*ZDC)*PLAT(JROF,JLEV)*ZRIPI2_(IJ_,ILA+1)

      ZLO   =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA  )
      ILO   =INT(ZLO )
      PDLO5(JROF,JLEV,0)=ZLO -REAL(ILO ,JPRB)
      PDLO(JROF,JLEV,0)=PLON(JROF,JLEV)*ZLSDEPI_(IJ_,ILA  )
      ZLO1  =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+1)
      ILO1  =INT(ZLO1)
      PDLO5(JROF,JLEV,1)=ZLO1-REAL(ILO1,JPRB)
      PDLO(JROF,JLEV,1)=PLON(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+1)
      ZLO2  =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+2)
      ILO2  =INT(ZLO2)
      PDLO5(JROF,JLEV,2)=ZLO2-REAL(ILO2,JPRB)
      PDLO(JROF,JLEV,2)=PLON(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+2)
      ZLO3  =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+3)
      ILO3  =INT(ZLO3)
      PDLO5(JROF,JLEV,3)=ZLO3-REAL(ILO3,JPRB)
      PDLO(JROF,JLEV,3)=PLON(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+3)

      ILEVV=KVAUT(INT(PLEV5(JROF,JLEV)*ZFAC))-1
      IF(ILEVV < IFLVM2.AND.&
       & (PLEV5(JROF,JLEV)-PVETA(ILEVV+2)) > 0.0_JPRB) ILEVV=ILEVV+1  

      KLEV(JROF,JLEV)=ILEVV
      PDVER5(JROF,JLEV)=(PLEV5(JROF,JLEV)-PVETA(ILEVV+1))/&
       & (PVETA(ILEVV+2)-PVETA(ILEVV+1))  
      PDVER(JROF,JLEV)=PLEV(JROF,JLEV)/(PVETA(ILEVV+2)-PVETA(ILEVV+1))

      KL0(JROF,JLEV,0)=IADDR_(IJ_,ILA  )+YDSL%NSLEXT(ILO ,ILA  )
      KL0(JROF,JLEV,1)=IADDR_(IJ_,ILA+1)+YDSL%NSLEXT(ILO1,ILA+1)
      KL0(JROF,JLEV,2)=IADDR_(IJ_,ILA+2)+YDSL%NSLEXT(ILO2,ILA+2)
      KL0(JROF,JLEV,3)=IADDR_(IJ_,ILA+3)+YDSL%NSLEXT(ILO3,ILA+3)

      KLH0(JROF,JLEV,0)=IADDR_(IJ_,ILA  )+YDSL%NSLEXT(ILO ,ILA  )+YDSL%NASLB1*ILEV
      KLH0(JROF,JLEV,1)=IADDR_(IJ_,ILA+1)+YDSL%NSLEXT(ILO1,ILA+1)+YDSL%NASLB1*ILEV
      KLH0(JROF,JLEV,2)=IADDR_(IJ_,ILA+2)+YDSL%NSLEXT(ILO2,ILA+2)+YDSL%NASLB1*ILEV
      KLH0(JROF,JLEV,3)=IADDR_(IJ_,ILA+3)+YDSL%NSLEXT(ILO3,ILA+3)+YDSL%NASLB1*ILEV

    ENDDO

    ! * Mask calculation for on-demand communications:
    IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
!CDIR NODEP
      DO JROF=KST,KPROF
        YDSL%MASK_SL2(KL0(JROF,JLEV,0)  )=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,0)+1)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,0)+2)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,0)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)  )=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)+1)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)+2)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)  )=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)+1)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)+2)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,3)  )=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,3)+1)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,3)+2)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,3)+3)=1
      ENDDO
    ENDIF

    DO JROF=KST,KPROF
      KL0(JROF,JLEV,0)=KL0(JROF,JLEV,0)+YDSL%NASLB1*KLEV(JROF,JLEV)
      KL0(JROF,JLEV,1)=KL0(JROF,JLEV,1)+YDSL%NASLB1*KLEV(JROF,JLEV)
      KL0(JROF,JLEV,2)=KL0(JROF,JLEV,2)+YDSL%NASLB1*KLEV(JROF,JLEV)
      KL0(JROF,JLEV,3)=KL0(JROF,JLEV,3)+YDSL%NASLB1*KLEV(JROF,JLEV)
    ENDDO

  ENDDO

  ! * Calculation of PCLA and PCLASLD
  !   + calculation of PCLA5 and PCLASLD5:
  CALL LASCAW_CLA_TL(YDSL,KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINH,&
   & IILA,ZZWH,PDLAT,PKAPPA,PKAPPAT,ZSLD1_,ZSLD2_,ZSLD3_,ZSLDW_,&
   & ZZWH5,PDLAT5,PKAPPA5,PKAPPAT5,&
   & PCLA,PCLASLD,PCLASLT,&
   & PCLA5,PCLASLD5,PCLASLT5)

  ! * Calculation of PCLO and PCLOSLD
  !   + calculation of PCLO5 and PCLOSLD5:
  CALL LASCAW_CLO_TL(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINH,&
   & PDLO(:,:,1),PKAPPA,PKAPPAT,&
   & PDLO5(:,:,1),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,1),PCLOSLD(:,:,:,1),PCLOSLT(:,:,:,1),&
   & PCLO5(:,:,:,1),PCLOSLD5(:,:,:,1),PCLOSLT5(:,:,:,1))
  CALL LASCAW_CLO_TL(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINH,&
   & PDLO(:,:,2),PKAPPA,PKAPPAT,&
   & PDLO5(:,:,2),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,2),PCLOSLD(:,:,:,2),PCLOSLT(:,:,:,2),&
   & PCLO5(:,:,:,2),PCLOSLD5(:,:,:,2),PCLOSLT5(:,:,:,2))

   ! * Calculation of PVINTW and PVINTWSLD
  !   + calculation of PVINTW5 and PVINTWSLD5:
  IF (KWIS == 102) LLT_SLHD(2) = (HOISLTV /= 0._JPRB)

  IF (((KWIS == 102).AND.LSLTVWENO) .OR. (KWIS == 106)) THEN

    ! Set value for boundary offset
    IF (LREGETA) THEN
      IBCLIM=0
    ELSE
      IBCLIM=1
    ENDIF

    ! WENO computation
    KNOWENO(KST:KPROF,1:KFLEV) = 0
    DO JJ=3,1,-1

      SELECT CASE (JJ)
        CASE (1)
          IWLEV(KST:KPROF,1:KFLEV)=KLEV(KST:KPROF,1:KFLEV)
        CASE (2)
          IWLEV(KST:KPROF,1:KFLEV)=MIN(IFLVM2-IBCLIM,KLEV(KST:KPROF,1:KFLEV)+1) 
          KNOWENO(KST:KPROF,1:KFLEV)=KNOWENO(KST:KPROF,1:KFLEV) &
           & + IWLEV(KST:KPROF,1:KFLEV)-KLEV(KST:KPROF,1:KFLEV) - 1
        CASE (3)
          ! can't be  KSLEV-1 as there is no half level on -1
          IWLEV(KST:KPROF,1:KFLEV)=MAX(IBCLIM,KLEV(KST:KPROF,1:KFLEV)-1)
          KNOWENO(KST:KPROF,1:KFLEV)=KNOWENO(KST:KPROF,1:KFLEV) &
           & + IWLEV(KST:KPROF,1:KFLEV)-KLEV(KST:KPROF,1:KFLEV) + 1
        CASE DEFAULT
          CALL ABOR1(' LASCAWTL: WENO PROBLEM')
      END SELECT    
        
      CALL LASCAW_VINTW_TL(KPROMB,KFLEV,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINV_WENO,IWLEV,&
       & PLEV,PDVER,PKAPPA,PKAPPAT,PVETA,ZVCUICO_,ZVSLD_,ZVSLDW_,&
       & PLEV5,PDVER5,PKAPPA5,PKAPPAT5,&
       & PVINTW(:,:,3*(JJ-1)+1),PVINTWSLD,PVINTWSLT,&
       & PVINTW5(:,:,3*(JJ-1)+1),PVINTWSLD5,PVINTWSLT5)

    ENDDO
    ! make sure it only keeps -1,0,+1 values
    if ((maxval(KNOWENO(KST:KPROF,1:KFLEV)) > 1+IBCLIM) .or. &
     &  (minval(KNOWENO(KST:KPROF,1:KFLEV)) <-1-IBCLIM))     &
     &  call abor1(' LASCAWTL: Something strange is happenig about level shifts.')

    ! C_k functions 
    IF (LREGETA) THEN
      DO JLEV=1,KFLEV
        DO JROF=KST,KPROF
          ! regular mesh (LREGETA=.t. case) but convenient for the BCs
          PCW5(JROF,JLEV,1)=0.10_JPRB*(2.0_JPRB+PDVER5(JROF,JLEV))*(3.0_JPRB-PDVER5(JROF,JLEV))   ! central
          PCW5(JROF,JLEV,2)=0.05_JPRB*(2.0_JPRB+PDVER5(JROF,JLEV))*(1.0_JPRB+PDVER5(JROF,JLEV))   ! lower
          PCW5(JROF,JLEV,3)=0.05_JPRB*(2.0_JPRB-PDVER5(JROF,JLEV))*(3.0_JPRB-PDVER5(JROF,JLEV))   ! upper

          PCW (JROF,JLEV,1)=0.10_JPRB*( 1.0_JPRB-2.0_JPRB*PDVER5(JROF,JLEV))*PDVER(JROF,JLEV)     ! central
          PCW (JROF,JLEV,2)=0.05_JPRB*( 3.0_JPRB+2.0_JPRB*PDVER5(JROF,JLEV))*PDVER(JROF,JLEV)     ! lower
          PCW (JROF,JLEV,3)=0.05_JPRB*(-5.0_JPRB+2.0_JPRB*PDVER5(JROF,JLEV))*PDVER(JROF,JLEV)     ! upper
        ENDDO
      ENDDO
    ELSE
      DO JLEV=1,KFLEV
        DO JROF=KST,KPROF
          ! general form
          ILEVV=KLEV(JROF,JLEV)
          IF ((ILEVV > 1) .AND. (ILEVV < KFLEV-3)) THEN
            PCW5(JROF,JLEV,1)=PGAMMA_WENO(ILEVV,1) &
             & *(PLEV5(JROF,JLEV)-PVETA(ILEVV-1))*(PLEV5(JROF,JLEV)-PVETA(ILEVV+4))  ! central
            PCW5(JROF,JLEV,2)=PGAMMA_WENO(ILEVV,2) &
             & *(PLEV5(JROF,JLEV)-PVETA(ILEVV-1))*(PLEV5(JROF,JLEV)-PVETA(ILEVV  ))  ! lower
            PCW5(JROF,JLEV,3)=PGAMMA_WENO(ILEVV,3) &
             & *(PLEV5(JROF,JLEV)-PVETA(ILEVV+3))*(PLEV5(JROF,JLEV)-PVETA(ILEVV+4))  ! upper

            PCW(JROF,JLEV,1)=PGAMMA_WENO(ILEVV,1) &
             & *( 2._JPRB*PLEV5(JROF,JLEV)-PVETA(ILEVV+4)-PVETA(ILEVV-1))*PLEV(JROF,JLEV)  ! central
            PCW(JROF,JLEV,2)=PGAMMA_WENO(ILEVV,2) &
             & *( 2._JPRB*PLEV5(JROF,JLEV)-PVETA(ILEVV  )-PVETA(ILEVV-1))*PLEV(JROF,JLEV)  ! lower
            PCW(JROF,JLEV,3)=PGAMMA_WENO(ILEVV,3) &
             & *( 2._JPRB*PLEV5(JROF,JLEV)-PVETA(ILEVV+4)-PVETA(ILEVV+3))*PLEV(JROF,JLEV)  ! upper
          ENDIF
        ENDDO
      ENDDO
    ENDIF

  ELSE
    CALL LASCAW_VINTW_TL(KPROMB,KFLEV,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINV,KLEV,&
     & PLEV,PDVER,PKAPPA,PKAPPAT,PVETA,ZVCUICO_,ZVSLD_,ZVSLDW_,&
     & PLEV5,PDVER5,PKAPPA5,PKAPPAT5,&
     & PVINTW,PVINTWSLD,PVINTWSLT,&
     & PVINTW5,PVINTWSLD5,PVINTWSLT5)
  ENDIF

  IF (KWIS == 104) THEN !!! TL not supported !!!
    CALL ABOR1('LASCAWTL: The option KWIS=104 is not supported by TL code')
  ENDIF

  IF (KWIS == 105) THEN !!! TL not supported !!!
    CALL ABOR1('LASCAWTL: The option KWIS=105 is not supported by TL code')
  ENDIF

ENDIF

!     ----------------------------------------------------------------

!*       2.    2D MODEL AND CASES IN THE 3D MODEL WHERE ONLY
!              2D INTERPOLATIONS ARE NEEDED.
!              ---------------------------------------------

!        2.01  Coordinates and weights for bilinear interpolations.

IF (KWIS == 201) THEN

  DO JLEV=1,KFLEV
    IF (KHOR == 0) ILEV=JLEV
    IF (KHOR == 1) ILEV=KFLDN

    ! * Calculation of linear weights, KL0.
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1

      IZLAT =INT(P4JP*(PIS2-PLAT5(JROF,JLEV))+0.75_JPRB)-YDSL%NFRSTLOFF
      ILA   =IZLAT+NINT(SIGN(0.5_JPRB,ZRLATI_(IJ_,IZLAT)-PLAT5(JROF,JLEV))-1.5_JPRB)
      PDLAT5(JROF,JLEV)=(PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+1))&
       & /(ZRLATI_(IJ_,ILA+2)-ZRLATI_(IJ_,ILA+1))  
      PDLAT(JROF,JLEV)=PLAT(JROF,JLEV)  /(ZRLATI_(IJ_,ILA+2)-ZRLATI_(IJ_,ILA+1))

      ZLO1  =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+1)
      ILO1  =INT(ZLO1)
      PDLO5(JROF,JLEV,1)=ZLO1-REAL(ILO1,JPRB)
      PDLO(JROF,JLEV,1)=PLON(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+1)
      ZLO2  =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+2)
      ILO2  =INT(ZLO2)
      PDLO5(JROF,JLEV,2)=ZLO2-REAL(ILO2,JPRB)
      PDLO(JROF,JLEV,2)=PLON(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+2)

      KL0(JROF,JLEV,1)=IADDR_(IJ_,ILA+1)+YDSL%NSLEXT(ILO1,ILA+1)+YDSL%NASLB1*ILEV
      KL0(JROF,JLEV,2)=IADDR_(IJ_,ILA+2)+YDSL%NSLEXT(ILO2,ILA+2)+YDSL%NASLB1*ILEV

    ENDDO

    ! * Mask calculation for on-demand communications:
    IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
!CDIR NODEP
      DO JROF=KST,KPROF
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)  )=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)+1)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)+2)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)  )=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)+1)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)+2)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)+3)=1
      ENDDO
    ENDIF

  ENDDO

ENDIF

!        2.03  Coordinates and weights for 12 points interpolations.

IF (KWIS == 202 .OR. KWIS == 203) THEN

  DO JLEV=1,KFLEV
    IF (KHOR == 0) ILEV=JLEV
    IF (KHOR == 1) ILEV=KFLDN

    ! * Calculation of linear weights, KL0.
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1

      IZLAT =INT(P4JP*(PIS2-PLAT5(JROF,JLEV))+0.75_JPRB)-YDSL%NFRSTLOFF
      ILA   =IZLAT+NINT(SIGN(0.5_JPRB,ZRLATI_(IJ_,IZLAT)-PLAT5(JROF,JLEV))-1.5_JPRB)
      PDLAT5(JROF,JLEV)=(PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+1))&
       & /(ZRLATI_(IJ_,ILA+2)-ZRLATI_(IJ_,ILA+1))  
      PDLAT(JROF,JLEV)=PLAT(JROF,JLEV) /(ZRLATI_(IJ_,ILA+2)-ZRLATI_(IJ_,ILA+1))

      ZDA   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA)
      ZDB   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+1)
      ZDC   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+2)
      ZDD   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+3)

      ! for input to LASCAW_CLA_TL.
      IILA(JROF,JLEV)=ILA
      ZZWH5(JROF,1,JLEV)=(ZDA*ZDC)*ZDD*ZRIPI0_(IJ_,ILA+1)
      ZZWH5(JROF,2,JLEV)=(ZDA*ZDB)*ZDD*ZRIPI1_(IJ_,ILA+1)
      ZZWH5(JROF,3,JLEV)=(ZDA*ZDB)*ZDC*ZRIPI2_(IJ_,ILA+1)
      ZZWH(JROF,1,JLEV)=(ZDA*(ZDC+ZDD)+ZDC*ZDD)*PLAT(JROF,JLEV)*ZRIPI0_(IJ_,ILA+1)
      ZZWH(JROF,2,JLEV)=(ZDA*(ZDB+ZDD)+ZDB*ZDD)*PLAT(JROF,JLEV)*ZRIPI1_(IJ_,ILA+1)
      ZZWH(JROF,3,JLEV)=(ZDA*(ZDB+ZDC)+ZDB*ZDC)*PLAT(JROF,JLEV)*ZRIPI2_(IJ_,ILA+1)

      ZLO   =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA  )
      ILO   =INT(ZLO )
      PDLO5(JROF,JLEV,0)=ZLO -REAL(ILO ,JPRB)
      PDLO(JROF,JLEV,0)=PLON(JROF,JLEV)*ZLSDEPI_(IJ_,ILA  )
      ZLO1  =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+1)
      ILO1  =INT(ZLO1)
      PDLO5(JROF,JLEV,1)=ZLO1-REAL(ILO1,JPRB)
      PDLO(JROF,JLEV,1)=PLON(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+1)
      ZLO2  =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+2)
      ILO2  =INT(ZLO2)
      PDLO5(JROF,JLEV,2)=ZLO2-REAL(ILO2,JPRB)
      PDLO(JROF,JLEV,2)=PLON(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+2)
      ZLO3  =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+3)
      ILO3  =INT(ZLO3)
      PDLO5(JROF,JLEV,3)=ZLO3-REAL(ILO3,JPRB)
      PDLO(JROF,JLEV,3)=PLON(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+3)

      KL0(JROF,JLEV,0)=IADDR_(IJ_,ILA  )+YDSL%NSLEXT(ILO ,ILA  )+YDSL%NASLB1*ILEV
      KL0(JROF,JLEV,1)=IADDR_(IJ_,ILA+1)+YDSL%NSLEXT(ILO1,ILA+1)+YDSL%NASLB1*ILEV
      KL0(JROF,JLEV,2)=IADDR_(IJ_,ILA+2)+YDSL%NSLEXT(ILO2,ILA+2)+YDSL%NASLB1*ILEV
      KL0(JROF,JLEV,3)=IADDR_(IJ_,ILA+3)+YDSL%NSLEXT(ILO3,ILA+3)+YDSL%NASLB1*ILEV

    ENDDO

    ! * Mask calculation for on-demand communications:
    IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
!CDIR NODEP
      DO JROF=KST,KPROF
        YDSL%MASK_SL2(KL0(JROF,JLEV,0)  )=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,0)+1)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,0)+2)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,0)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)  )=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)+1)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)+2)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)  )=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)+1)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)+2)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,3)  )=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,3)+1)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,3)+2)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,3)+3)=1
      ENDDO
    ENDIF

  ENDDO

  ! * Calculation of PCLA and PCLASLD
  !   + calculation of PCLA5 and PCLASLD5:
  CALL LASCAW_CLA_TL(YDSL,KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINH,&
   & IILA,ZZWH,PDLAT,PKAPPA,PKAPPAT,ZSLD1_,ZSLD2_,ZSLD3_,ZSLDW_,&
   & ZZWH5,PDLAT5,PKAPPA5,PKAPPAT5,&
   & PCLA,PCLASLD,PCLASLT,&
   & PCLA5,PCLASLD5,PCLASLT5)

  ! * Calculation of PCLO and PCLOSLD
  !   + calculation of PCLO5 and PCLOSLD5:
  CALL LASCAW_CLO_TL(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINH,&
   & PDLO(:,:,1),PKAPPA,PKAPPAT,&
   & PDLO5(:,:,1),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,1),PCLOSLD(:,:,:,1),PCLOSLT(:,:,:,1),&
   & PCLO5(:,:,:,1),PCLOSLD5(:,:,:,1),PCLOSLT5(:,:,:,1))
  CALL LASCAW_CLO_TL(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINH,&
   & PDLO(:,:,2),PKAPPA,PKAPPAT,&
   & PDLO5(:,:,2),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,2),PCLOSLD(:,:,:,2),PCLOSLT(:,:,:,2),&
   & PCLO5(:,:,:,2),PCLOSLD5(:,:,:,2),PCLOSLT5(:,:,:,2))

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LASCAWTL',1,ZHOOK_HANDLE)
END SUBROUTINE LASCAWTL
