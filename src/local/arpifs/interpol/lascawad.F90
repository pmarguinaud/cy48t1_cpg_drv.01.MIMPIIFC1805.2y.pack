!ocl list_copy(32,PLSDEPI)
!ocl list_copy(32,PLATI)
!ocl list_copy(32,PIPI)
!ocl list_copy(32,PVETA)
!ocl list_copy(32,PVCUICO)
SUBROUTINE LASCAWAD(YDSL,KPROMB,KST,KPROF,KFLEV,&
 & KFLDN,KWIS,KHOR,KWENO,&
 & LDSLHD,LDSLHDQUAD,LDSLHD_OLD,LDSLHDHEAT,&
 & P4JP,PIS2,PLSDEPI,PLATI,&
 & PIPI,PSLD,PSLDW,PGAMMA_WENO,&
 & PLON,PLAT,PLEV,PKAPPA,PKAPPAT,&
 & PLON5,PLAT5,PLEV5,PKAPPA5,PKAPPAT5,PDVER5,&
 & PVETA,KVAUT,&
 & PVCUICO,PVSLD,PVSLDW,KRLEVX,PVRLEVX,&
 & PDLAT,PCLA,PCLASLD,PCLASLT,PDLO,PCLO,PCLOSLD,PCLOSLT,PCW,PDVER,PVINTW,PVINTWSLD,PVINTWSLT)  

!**** *LASCAWAD - Externalisable interpolator:   (adjoint version)
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
!        *CALL* *LASCAWAD(.........)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          YDSL    - SL_STRUCT definition
!          KPROMB  - horizontal dimension for interpolation point quantities.
!          KST     - first element of arrays where computations are performed.
!          KPROF   - depth of work.
!          KFLEV   - vertical dimension.
!          KFLDN   - number of the first field.
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
!          PGAMMA_WENO - weights for vertical WENO interpolation
! ------------------ input variables (output in TL) ----------------
!          PDLAT     - distance for horizontal linear interpolations
!                      in latitude.
!          PCLA      - weights for horizontal cubic interpolations
!                      in latitude.
!          PCLASLD   - weights for horizontal cubic interpolations
!                      in latitude, SLHD case
!          PCLASLT   - weights for horizontal cubic interpolations
!                      in latitude, SLHD case on T
!          PCW       - C_k weights for the vertical WENO interpolation
!          PDLO      - distances for horizontal linear interpolations
!                      in longitude (latitude rows 0, 1, 2, 3)
!          PCLO      - weights for horizontal cubic interpolations in
!                      longitude (latitude rows 1, 2)
!          PCLOSLD   - weights for horizontal cubic interpolations in
!                      longitude, SLHD case (latitude rows 1, 2)
!          PCLOSLT   - weights for horizontal cubic interpolations in
!                      longitude, SLHD case (latitude rows 1, 2) on T
!          PDVER     - distance for vertical linear interpolation
!          PVINTW    - vertical cubic interpolation weights
!          PVINTWSLD - vertical cubic interpolation weights, SLHD case
!          PVINTWSLT - vertical cubic interpolation weights, SLHD case on T
! ---------------------- trajectory variables ------------------------
!          PLON5   - Interpolation point longitude on the work sphere.
!          PLAT5   - Interpolation point latitude on the work sphere.
!          PLEV5   - vertical coordinate of the interpolation point.
!          PKAPPA(T)5 - kappa function ("coefficient of SLHD") based on the
!                    rescaled horizontal deformation of the flow evaluated
!                    at instant "t" for the final point F for momentum and heat
!          PDVER5    - distance for vertical linear interpolation
! --------------------------------------------------------------------
!          PVETA   - Values of ETA.
!          KVAUT   - Help table for vertical box search: gives the number
!                    of the layer immediately above eta.
!          PVCUICO - Denominators of the vertical cubic interpolation
!                    coefficients
!          PVSLD   - auxiliary quantities for vertical SLHD interpolation
!          PVSLDW  - weights for SLHD vertical Laplacian smoother
!          KRLEVX  - Dimension of KVAUT
!          PVRLEVX - REAL(KRLEVX).

! ------------ output variables (input in TL) ----------------------
!          PLON    - Interpolation point longitude on the work sphere.
!          PLAT    - Interpolation point latitude on the work sphere.
!          PLEV    - vertical coordinate of the interpolation point.
!          PKAPPA  - kappa function ("coefficient of SLHD") based on the
!                    rescaled horizontal deformation of the flow evaluated
!                    at instant "t" for the final point F
!          PKAPPAT - kappa function for heat (scalar) variables

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
!      Original : 98/10/19.

!     Modifications.
!     --------------
!      Modified 04-Aug-2008 : F. Vana - weights driven interpolation fully
!                            projected to TL (still excluding P/C and half
!                            level interpolation) + optimization for vect.
!      K. Yessad (Dec 2008): remove useless dummy arguments
!      F. Vana  13-Jan-2009: fixed bug in PKAPPA computation
!      K. Yessad + F. Vana (Feb 2009): split loops, rewrite in a shorter way.
!      K. Yessad (Aug 2009): use RIPI, RSLD
!      G. Mozdzynski (May 2012): further cleaning
!      K. Yessad (July 2014): Move some variables.
!      F. Vana 13-feb-2014 SLHD weights for heat variables
!      K. Yessad (March 2017): simplify level numbering.
!      F. Vana    21-Nov-2017: Option LHOISLT
!      F. Vana  26-Feb-2019: quintic vertical interpolation
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK

! arp/ifs dependencies to be solved later.
USE YOMMP0   , ONLY : LOPT_SCALAR
USE YOMCT0   , ONLY : LREGETA
USE YOMDYNA  , ONLY : YRDYNA

USE EINT_MOD , ONLY : SL_STRUCT, JPDUP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(SL_STRUCT),   INTENT(IN)    :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMB 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KRLEVX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAMMA_WENO(KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLON(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLAT(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLEV(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PKAPPA(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PKAPPAT(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLON5(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT5(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLEV5(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPPA5(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPPAT5(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDVER5(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVETA(0:KFLEV+1) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KVAUT(0:KRLEVX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCUICO(4,0:KFLEV-1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSLD(3,0:KFLEV-1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSLDW(3,3,0:KFLEV-1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRLEVX 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDLAT(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLA(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLASLD(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLASLT(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDLO(KPROMB,KFLEV,0:3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLO(KPROMB,KFLEV,3,2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLOSLD(KPROMB,KFLEV,3,2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLOSLT(KPROMB,KFLEV,3,2)
REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PCW(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDVER(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVINTW(KPROMB,KFLEV,3*KWENO)
REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PVINTWSLD(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PVINTWSLT(KPROMB,KFLEV,3)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IFLVM2, ILA, ILEV, ILEVV, IBCLIM,&
 & ILO1, ILO2, IZLAT, JLEV, JROF, IJ_, J_, JJ
INTEGER(KIND=JPIM) :: IILA(KPROMB,KFLEV), IILEV(KPROMB,KFLEV)

REAL(KIND=JPRB) ::  ZDA, ZDB, ZDC, ZDD, ZFAC, ZLO1, ZLO2, ZEPS
REAL(KIND=JPRB) ::  ZDLO5(KPROMB,KFLEV,2), ZDLAT5(KPROMB,KFLEV), ZDVER5(KPROMB,KFLEV)
INTEGER(KIND=JPIM) :: IWLEV(KPROMB,KFLEV)

LOGICAL :: LLT_SLHD(3),LLSLHD,LLSLHDQUAD,LLSLHD_OLD

! Duplicata of some dummy or local arrays for optimisation on NEC platform.
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
#include "lascaw_cla_ad.intfb.h"
#include "lascaw_clo_ad.intfb.h"
#include "lascaw_vintw_ad.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LASCAWAD',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       0.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

ZEPS=100.0_JPRB*TINY(1.0_JPRB)

! cases relevant for SLHD scheme (switches LDSLHD, LDSLHDQUAD are
! deactivated during computation of medium points in LAPINEA)
LLSLHD=LDSLHD.AND.(KWIS==103.OR.KWIS==106.OR.KWIS==203)
LLSLHDQUAD=(LDSLHDQUAD.AND.(KWIS==103.OR.KWIS==106.OR.KWIS==203)).OR. &
   & ((KWIS==102.OR.KWIS==202).AND.YRDYNA%LHOISLT.AND.((YRDYNA%HOISLTV/=0_JPRB).OR.(YRDYNA%HOISLTH/=0_JPRB)))

! switch for old SLHD scheme
LLSLHD_OLD=LLSLHD.AND.LDSLHD_OLD

LLT_SLHD(1)=LLSLHD
LLT_SLHD(2)=LLSLHDQUAD
LLT_SLHD(3)=LLSLHD_OLD

ZSLHDKMINH=YRDYNA%SLHDKMIN
ZSLHDKMINV=YRDYNA%SLHDKMIN

! Modify previous defaults for high order trajectory research
IF (KWIS==102.OR.KWIS==202) THEN
  ZSLHDKMINH=YRDYNA%HOISLTH
  ZSLHDKMINV=YRDYNA%HOISLTV  ! useless for shallow-water
ENDIF

IFLVM2=KFLEV-2

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

!     ------------------------------------------------------------------

!*       1.    3D MODEL.
!              ---------

!        1.01  Coordinates and weights for trilinear interpolations.

IF (KWIS == 101) THEN

  ZFAC=PVRLEVX/(PVETA(KFLEV+1)-PVETA(0))

  DO JLEV=1,KFLEV

    ! * AD calculation of linear weights.
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1

      IZLAT =INT(P4JP*(PIS2-PLAT5(JROF,JLEV))+0.75_JPRB+ZEPS)-YDSL%NFRSTLOFF
      ILA   =IZLAT+NINT(SIGN(0.5_JPRB,ZRLATI_(IJ_,IZLAT)-PLAT5(JROF,JLEV))-1.5_JPRB)

      PLAT(JROF,JLEV)=PLAT(JROF,JLEV)+PDLAT(JROF,JLEV)&
       & /(ZRLATI_(IJ_,ILA+2)-ZRLATI_(IJ_,ILA+1)+ZEPS)  

      PLON(JROF,JLEV)=PLON(JROF,JLEV)+&
       & PDLO(JROF,JLEV,1)*ZLSDEPI_(IJ_,ILA+1)+PDLO(JROF,JLEV,2)*ZLSDEPI_(IJ_,ILA+2)  

      ILEV  =KVAUT(INT(PLEV5(JROF,JLEV)*ZFAC))-1
      IF(ILEV < IFLVM2.AND.&
       & (PLEV5(JROF,JLEV)-PVETA(ILEV+2)) > 0.0_JPRB) ILEV=ILEV+1  

      PLEV(JROF,JLEV)=PLEV(JROF,JLEV)+PDVER(JROF,JLEV)/&
       & (PVETA(ILEV+2)-PVETA(ILEV+1))  

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

    ! * AD of calculation of linear weights: trajectory.
!CDIR NOLSTVAL
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1

      IZLAT =INT(P4JP*(PIS2-PLAT5(JROF,JLEV))+0.75_JPRB+ZEPS)-YDSL%NFRSTLOFF
      ILA   =IZLAT+NINT(SIGN(0.5_JPRB,ZRLATI_(IJ_,IZLAT)-PLAT5(JROF,JLEV))-1.5_JPRB)
      ZDLAT5(JROF,JLEV)=(PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+1))&
       & /(ZRLATI_(IJ_,ILA+2)-ZRLATI_(IJ_,ILA+1)+ZEPS)

      ZDA   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA)
      ZDB   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+1)
      ZDC   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+2)
      ZDD   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+3)

      ! for input to LASCAW_CLA_AD.
      IILA(JROF,JLEV)=ILA
      ZZWH5(JROF,1,JLEV)=(ZDA*ZDC)*ZDD*ZRIPI0_(IJ_,ILA+1)
      ZZWH5(JROF,2,JLEV)=(ZDA*ZDB)*ZDD*ZRIPI1_(IJ_,ILA+1)
      ZZWH5(JROF,3,JLEV)=(ZDA*ZDB)*ZDC*ZRIPI2_(IJ_,ILA+1)

      ZLO1  =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+1)
      ILO1  =INT(ZLO1)
      ZDLO5(JROF,JLEV,1)=ZLO1-REAL(ILO1,JPRB)
      ZLO2  =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+2)
      ILO2  =INT(ZLO2)
      ZDLO5(JROF,JLEV,2)=ZLO2-REAL(ILO2,JPRB)

      ILEVV=KVAUT(INT(PLEV5(JROF,JLEV)*ZFAC))-1
      IF(ILEVV < IFLVM2.AND.&
       & (PLEV5(JROF,JLEV)-PVETA(ILEVV+2)) > 0.0_JPRB) ILEVV=ILEVV+1  
      IILEV(JROF,JLEV)=ILEVV

      ZDVER5(JROF,JLEV)=(PLEV5(JROF,JLEV)-PVETA(ILEVV+1))/&
       & (PVETA(ILEVV+2)-PVETA(ILEVV+1)) 

    ENDDO
  ENDDO

  ! Modify previous defaults for high order trajectory research
  IF (KWIS==102) LLT_SLHD(2)= (YRDYNA%HOISLTV /= 0._JPRB)

  ! * AD of calculation of PVINTW and PVINTWSLD:
  IF (((KWIS == 102) .AND. YRDYNA%LSLTVWENO) .OR. (KWIS == 106)) THEN

    ! Set value for boundary offset
    IF (LREGETA) THEN
      IBCLIM=0
    ELSE
      IBCLIM=1
    ENDIF

    ! WENO computation

    ! C_k functions 
    IF (LREGETA) THEN
      DO JLEV=1,KFLEV
        DO JROF=KST,KPROF
          ! regular mesh 
          PDVER(JROF,JLEV)=PDVER(JROF,JLEV) &
           & + 0.10_JPRB*( 1.0_JPRB-2.0_JPRB*PDVER5(JROF,JLEV))*PCW(JROF,JLEV,1) &
           & + 0.05_JPRB*( 3.0_JPRB+2.0_JPRB*PDVER5(JROF,JLEV))*PCW(JROF,JLEV,2) &
           & + 0.05_JPRB*(-5.0_JPRB+2.0_JPRB*PDVER5(JROF,JLEV))*PCW(JROF,JLEV,3)
        ENDDO
      ENDDO
    ELSE
      DO JLEV=1,KFLEV
        DO JROF=KST,KPROF
          ! general form
          ILEVV=IILEV(JROF,JLEV)
          IF ((ILEVV > 1) .AND. (ILEVV < KFLEV-3)) PLEV(JROF,JLEV)=PLEV(JROF,JLEV) &
             & + PGAMMA_WENO(ILEVV,1)*(2._JPRB*PLEV5(JROF,JLEV)-PVETA(ILEVV+4)-PVETA(ILEVV-1))*PCW(JROF,JLEV,1) &
             & + PGAMMA_WENO(ILEVV,2)*(2._JPRB*PLEV5(JROF,JLEV)-PVETA(ILEVV  )-PVETA(ILEVV-1))*PCW(JROF,JLEV,2) &
             & + PGAMMA_WENO(ILEVV,3)*(2._JPRB*PLEV5(JROF,JLEV)-PVETA(ILEVV+4)-PVETA(ILEVV+3))*PCW(JROF,JLEV,3)
        ENDDO
      ENDDO
    ENDIF
    PCW(:,:,:)=0._JPRB

    DO JJ=1,3

      SELECT CASE (JJ)
        CASE (1)
          IWLEV(KST:KPROF,1:KFLEV)=IILEV(KST:KPROF,1:KFLEV)
        CASE (2)
          IWLEV(KST:KPROF,1:KFLEV)=MIN(IFLVM2-IBCLIM,IILEV(KST:KPROF,1:KFLEV)+1) 
        CASE (3)
          ! can't be  KSLEV-1 as there is no half level on -1
          IWLEV(KST:KPROF,1:KFLEV)=MAX(IBCLIM,IILEV(KST:KPROF,1:KFLEV)-1)
        CASE DEFAULT
          CALL ABOR1(' LASCAWAD: WENO PROBLEM')
      END SELECT    
        
      CALL LASCAW_VINTW_AD(KPROMB,KFLEV,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINV_WENO,IWLEV,&
       & PLEV,PDVER,PKAPPA,PKAPPAT,PVETA,ZVCUICO_,ZVSLD_,ZVSLDW_,&
       & PLEV5,ZDVER5,PKAPPA5,PKAPPAT5,&
       & PVINTW(:,:,3*(JJ-1)+1),PVINTWSLD,PVINTWSLT)
    ENDDO

  ELSE
    ! Standard high order interpolation
    CALL LASCAW_VINTW_AD(KPROMB,KFLEV,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINV,IILEV,&
     & PLEV,PDVER,PKAPPA,PKAPPAT,PVETA,ZVCUICO_,ZVSLD_,ZVSLDW_,&
     & PLEV5,ZDVER5,PKAPPA5,PKAPPAT5,&
     & PVINTW,PVINTWSLD,PVINTWSLT)
  ENDIF

  ! Modify previous defaults for high order trajectory research
  IF (KWIS==102) LLT_SLHD(2)= (YRDYNA%HOISLTH /= 0._JPRB)

    ! * AD of calculation of PCLO and PCLOSLD:
  CALL LASCAW_CLO_AD(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINH,&
   & PDLO(:,:,2),PKAPPA,PKAPPAT,&
   & ZDLO5(:,:,2),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,2),PCLOSLD(:,:,:,2),PCLOSLT(:,:,:,2))
  CALL LASCAW_CLO_AD(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINH,&
   & PDLO(:,:,1),PKAPPA,PKAPPAT,&
   & ZDLO5(:,:,1),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,1),PCLOSLD(:,:,:,1),PCLOSLT(:,:,:,1))

  ! * AD of calculation of PCLA and PCLASLD:
  CALL LASCAW_CLA_AD(YDSL,KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINH,&
   & IILA,ZZWH,PDLAT,PKAPPA,PKAPPAT,ZSLD1_,ZSLD2_,ZSLD3_,ZSLDW_,&
   & ZZWH5,ZDLAT5,PKAPPA5,PKAPPAT5,&
   & PCLA,PCLASLD,PCLASLT)
    

  DO JLEV=1,KFLEV
    ! * AD of calculation of linear weights: increments.
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1
      ILA=IILA(JROF,JLEV)

      ZDA   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA)
      ZDB   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+1)
      ZDC   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+2)
      ZDD   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+3)

      PLEV(JROF,JLEV)=PLEV(JROF,JLEV)+PDVER(JROF,JLEV)&
       & /(PVETA(IILEV(JROF,JLEV)+2)-PVETA(IILEV(JROF,JLEV)+1))

      PLON(JROF,JLEV)=PLON(JROF,JLEV)&
       & +PDLO(JROF,JLEV,0)*ZLSDEPI_(IJ_,ILA  )&
       & +PDLO(JROF,JLEV,1)*ZLSDEPI_(IJ_,ILA+1)&
       & +PDLO(JROF,JLEV,2)*ZLSDEPI_(IJ_,ILA+2)&
       & +PDLO(JROF,JLEV,3)*ZLSDEPI_(IJ_,ILA+3)

      PLAT(JROF,JLEV)=PLAT(JROF,JLEV)&
       & +(ZDA*(ZDC+ZDD)+ZDC*ZDD)*ZZWH(JROF,1,JLEV)*ZRIPI0_(IJ_,ILA+1)&
       & +(ZDA*(ZDB+ZDD)+ZDB*ZDD)*ZZWH(JROF,2,JLEV)*ZRIPI1_(IJ_,ILA+1)&
       & +(ZDA*(ZDB+ZDC)+ZDB*ZDC)*ZZWH(JROF,3,JLEV)*ZRIPI2_(IJ_,ILA+1)

      PLAT(JROF,JLEV)=PLAT(JROF,JLEV)+&
       & PDLAT(JROF,JLEV)/(ZRLATI_(IJ_,ILA+2)-ZRLATI_(IJ_,ILA+1)+ZEPS)

    ENDDO
  ENDDO

  IF (KWIS == 104) THEN     !!! AD not supported !!!
    CALL ABOR1('LASCAWAD: The option KWIS=104 is not supported by AD code')
  ENDIF

  IF (KWIS == 105) THEN !!! AD not supported !!!
    CALL ABOR1('LASCAWAD: The option KWIS=105 is not supported by AD code')
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

    ! * AD calculation of linear weights.
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1

      IZLAT =INT(P4JP*(PIS2-PLAT5(JROF,JLEV))+0.75_JPRB)-YDSL%NFRSTLOFF
      ILA   =IZLAT+NINT(SIGN(0.5_JPRB,ZRLATI_(IJ_,IZLAT)-PLAT5(JROF,JLEV))-1.5_JPRB)

      PLAT(JROF,JLEV)=PLAT(JROF,JLEV) +&
       & PDLAT(JROF,JLEV) / (ZRLATI_(IJ_,ILA+2)-ZRLATI_(IJ_,ILA+1))  

      PLON(JROF,JLEV)=PLON(JROF,JLEV)&
       & +PDLO(JROF,JLEV,1)*ZLSDEPI_(IJ_,ILA+1)&
       & +PDLO(JROF,JLEV,2)*ZLSDEPI_(IJ_,ILA+2)  

    ENDDO
  ENDDO

ENDIF

!        2.03  Coordinates and weights for 12 points interpolations.

IF (KWIS == 202 .OR. KWIS == 203) THEN

  DO JLEV=1,KFLEV
    IF (KHOR == 0) ILEV=JLEV
    IF (KHOR == 1) ILEV=KFLDN

    ! * AD of calculation of linear weights: trajectory.
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1

      IZLAT =INT(P4JP*(PIS2-PLAT5(JROF,JLEV))+0.75_JPRB)-YDSL%NFRSTLOFF
      ILA   =IZLAT+NINT(SIGN(0.5_JPRB,ZRLATI_(IJ_,IZLAT)-PLAT5(JROF,JLEV))-1.5_JPRB)
      ZDLAT5(JROF,JLEV)=(PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+1))&
       & /(ZRLATI_(IJ_,ILA+2)-ZRLATI_(IJ_,ILA+1))

      ZDA   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA)
      ZDB   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+1)
      ZDC   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+2)
      ZDD   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+3)

      ! for input to LASCAW_CLA_AD.
      IILA(JROF,JLEV)=ILA
      ZZWH5(JROF,1,JLEV)=(ZDA*ZDC)*ZDD*ZRIPI0_(IJ_,ILA+1)
      ZZWH5(JROF,2,JLEV)=(ZDA*ZDB)*ZDD*ZRIPI1_(IJ_,ILA+1)
      ZZWH5(JROF,3,JLEV)=(ZDA*ZDB)*ZDC*ZRIPI2_(IJ_,ILA+1)

      ZLO1  =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+1)
      ILO1  =INT(ZLO1)
      ZDLO5(JROF,JLEV,1)=ZLO1-REAL(ILO1,JPRB)
      ZLO2  =PLON5(JROF,JLEV)*ZLSDEPI_(IJ_,ILA+2)
      ILO2  =INT(ZLO2)
      ZDLO5(JROF,JLEV,2)=ZLO2-REAL(ILO2,JPRB)

    ENDDO
  ENDDO

  ! Update the LSLHDQUAD for horizontal interpolation
  IF (KWIS == 202) LLT_SLHD(2)= (YRDYNA%HOISLTH /= 0._JPRB)

  ! * AD of calculation of PCLO and PCLOSLD:
  CALL LASCAW_CLO_AD(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINH,&
   & PDLO(:,:,2),PKAPPA,PKAPPAT,&
   & ZDLO5(:,:,2),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,2),PCLOSLD(:,:,:,2),PCLOSLT(:,:,:,2))
  CALL LASCAW_CLO_AD(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINH,&
   & PDLO(:,:,1),PKAPPA,PKAPPAT,&
   & ZDLO5(:,:,1),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,1),PCLOSLD(:,:,:,1),PCLOSLT(:,:,:,1))

  ! * AD of calculation of PCLA and PCLASLD:
  CALL LASCAW_CLA_AD(YDSL,KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINH,&
   & IILA,ZZWH,PDLAT,PKAPPA,PKAPPAT,ZSLD1_,ZSLD2_,ZSLD3_,ZSLDW_,&
   & ZZWH5,ZDLAT5,PKAPPA5,PKAPPAT5,&
   & PCLA,PCLASLD,PCLASLT)

  DO JLEV=1,KFLEV
    ! * AD of calculation of linear weights: increments.
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1
      ILA=IILA(JROF,JLEV)

      ZDA   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA)
      ZDB   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+1)
      ZDC   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+2)
      ZDD   =PLAT5(JROF,JLEV)-ZRLATI_(IJ_,ILA+3)

      PLON(JROF,JLEV)=PLON(JROF,JLEV)&
       & +PDLO(JROF,JLEV,0)*ZLSDEPI_(IJ_,ILA  )&
       & +PDLO(JROF,JLEV,1)*ZLSDEPI_(IJ_,ILA+1)&
       & +PDLO(JROF,JLEV,2)*ZLSDEPI_(IJ_,ILA+2)&
       & +PDLO(JROF,JLEV,3)*ZLSDEPI_(IJ_,ILA+3)

      PLAT(JROF,JLEV)=PLAT(JROF,JLEV)&
       & +(ZDA*(ZDC+ZDD)+ZDC*ZDD)*ZZWH(JROF,1,JLEV)*ZRIPI0_(IJ_,ILA+1)&
       & +(ZDA*(ZDB+ZDD)+ZDB*ZDD)*ZZWH(JROF,2,JLEV)*ZRIPI1_(IJ_,ILA+1)&
       & +(ZDA*(ZDB+ZDC)+ZDB*ZDC)*ZZWH(JROF,3,JLEV)*ZRIPI2_(IJ_,ILA+1)

      PLAT(JROF,JLEV)=PLAT(JROF,JLEV)+&
       & PDLAT(JROF,JLEV)/(ZRLATI_(IJ_,ILA+2)-ZRLATI_(IJ_,ILA+1))

    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LASCAWAD',1,ZHOOK_HANDLE)
END SUBROUTINE LASCAWAD
