SUBROUTINE ELASCAWTL(&
 ! --- INPUT -------------------------------------------------
 & YDDYN,YDSL,KPROMB,KST,KPROF,KFLEV,&
 & KFLDN,KSTABUF,KWIS,KHOR,&
 & LDSLHD,LDSLHDQUAD,LDSLHD_OLD,&
 & PLON,PLAT,PLEV,&
 & PLON5,PLAT5,PLEV5,&
 & PVETA,KVAUT,&
 & PVCUICO,PVSLD,PVSLDW,KRLEVX,PVRLEVX,PKAPPA,PKAPPA5,PKAPPAT,PKAPPAT5,&
 ! --- OUTPUT ------------------------------------------------
 & PDLAT,PCLA,PCLASLD,PCLASLT,PDLO,PCLO,PCLOSLD,PCLOSLT,&
 & KL0,KLH0,KLEV,PDVER,PVINTW,PVINTWSLD,PVINTWSLT,&
 & PDLAT5,PCLA5,PCLASLD5,PCLASLT5,PDLO5,PCLO5,PCLOSLD5,PCLOSLT5,&
 & PDVER5,PVINTW5,PVINTWSLD5,PVINTWSLT5)  

!**** *ELASCAWTL  -  Externalisable interpolator:   (TL code)
!                 Storage of Coordinates And Weights.
!                 Plane geometry version

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
!        *CALL* *ELASCAWTL( ... )

!        Explicit arguments :
!        --------------------

!        INPUT:
!          YDDYN   - structure containing dynamics
!          YDSL    - SL_STRUCT definition
!          KPROMB  - horizontal dimension for interpolation point quantities.
!          KST     - first element of arrays where computations are performed.
!          KPROF   - depth of work.
!          KFLEV   - vertical dimension.
!          KFLDN   - number of the first field.
!          KSTABUF - for a latitude IGL, KSTABUF(IGL) is the
!                    address of the element corresponding to
!                    (ILON=1,IGL) in the KPROMB arrays.
!          KWIS    - kind of interpolation.
!          KHOR    - 0: Horizontal interpolation for each level
!                       of a 3D variable.
!                    1: Interpolation for all origin points corresponding
!                       to a final point of a 2D variable.
!          LDSLHD  - key activating SLHD weights precomputation.
!          LDSLHDQUAD - key activating quadratic weights precomputation.
!          LDSLHD_OLD - use old SLHD interpolator
!          PLON    - x-coordinate of the interpolation point
!                    (in the fractional system <NDLUNG;NDLON>)
!          PLAT    - y-coordinate of the interpolation point
!                    (in the fractional system <NDGSAG;NDGENG>)
!          PLEV    - vertical coordinate of the interpolation point.
!  ------------------- Trajectory variables ------------------------
!          PLON5   - x-coordinate of the interpolation point
!                    (in the fractional system <NDLUNG;NDLON>)
!          PLAT5   - y-coordinate of the interpolation point
!                    (in the fractional system <NDGSAG;NDGENG>)
!          PLEV5   - vertical coordinate of the interpolation point.
!  -----------------------------------------------------------------
!          PVETA   - Values of ETA.
!          KVAUT   - Help table for vertical box search: gives the number
!                    of the layer or interlayer immediately above eta.
!          PVCUICO - Denominators of the vertical cubic interpolation coefficients
!          PVSLD   - auxiliary quantities for vertical SLHD interpolation
!          PVSLDW  - weights for SLHD vertical Laplacian smoother
!          KRLEVX  - Dimension of KVAUT
!          PVRLEVX - REAL(KRLEVX).
!          PKAPPA  - kappa function ("coefficient of SLHD") based on the
!                    rescaled horizontal deformation of the flow evaluated
!                    at instant "t" for the final point F
! ---------------------- trajectory variable -------------------------
!          PKAPPA5 - kappa function ("coefficient of SLHD") based on the
!                    rescaled horizontal deformation of the flow evaluated
!                    at instant "t" for the final point F
! --------------------------------------------------------------------

!        OUTPUT:
!          PDLAT     - distance for horizontal linear interpolations in latitude
!          PCLA      - weights for horizontal cubic interpolations in latitude
!          PCLASLD   - weights for horizontal cubic interpolations in latitude, SLHD case
!          PDLO      - distances for horizontal linear interpolations
!                      in longitude (latitude rows 0, 1, 2, 3)
!          PCLO      - weights for horizontal cubic interpolations in
!                      longitude (latitude rows 1, 2)
!          PCLOSLD   - weights for horizontal cubic interpolations in
!                      longitude, SLHD case (latitude rows 1, 2)
!          KL0       - index of the four western points
!                      of the 16 points interpolation grid.
!          KLH0      - second value of index of the four western points
!                      of the 16 points interpolation grid if needed.
!          KLEV      - lower level of the vertical interpolation
!                      grid needed for vertical interpolations.
!          PDVER     - distance for vertical linear interpolation
!          PVINTW    - vertical cubic interpolation weights
!          PVINTWSLD - vertical cubic interpolation weights, SLHD case

!  ------------------- Trajectory variables ------------------------
!          PDLAT5    - distance for horizontal linear interpolations in latitude.
!          PCLA5     - weights for horizontal cubic interpolations in latitude.
!          PCLASLD5  - weights for horizontal cubic interpolations in latitude, SLHD case
!          PDLO5     - distances for horizontal linear interpolations
!                      in longitude (latitude rows 0, 1, 2, 3)
!          PCLO5     - weights for horizontal cubic interpolations in
!                      longitude (latitude rows 1, 2)
!          PCLOSLD5  - weights for horizontal cubic interpolations in
!                      longitude, SLHD case (latitude rows 1, 2)
!          PDVER5    - distance for vertical linear interpolation
!          PVINTW5   - vertical cubic interpolation weights
!          PVINTWSLD5- vertical cubic interpolation weights, SLHD case
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
!        F. VANA, being inspired by LASCAWTL
!        * ALADIN/LACE *
!        Original : JUNE 2006.

!     Modifications.
!     --------------
!        Modified 02-Jun-2008 : A. Bogatchev - USE YOMLASCAW et c.
!        Modified 13-May-2008 : F. Vana - weights driven interpolation
!        K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!        K. Yessad (Feb 2009): split loops, rewrite in a shorter way.
!        G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!        B. Bochenek (Apr 2015): Phasing: move some variables and update
!        K. Yessad (March 2017): simplify level numbering.
!        R. El Khatib 22-May-2019 fix the "IFL OUT OF BOUNDS" issue
!     ------------------------------------------------------------------

USE YOMDYN   , ONLY : TDYN
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

! arp/ifs dependencies to be solved later.
USE YOMMP0   , ONLY : LOPT_SCALAR, NPROC
USE YOMDYNA  , ONLY : SLHDKMIN

USE EINT_MOD , ONLY : SL_STRUCT,JPDUP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
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
LOGICAL           ,INTENT(IN)    :: LDSLHD
LOGICAL           ,INTENT(IN)    :: LDSLHDQUAD
LOGICAL           ,INTENT(IN)    :: LDSLHD_OLD
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
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDVER(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTW(KPROMB,KFLEV,3) 
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
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDVER5(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTW5(KPROMB,KFLEV,3) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTWSLD5(KPROMB,KFLEV,3) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTWSLT5(KPROMB,KFLEV,3) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IADDR(YDSL%NDGSAH:YDSL%NDGENH)
INTEGER(KIND=JPIM) :: IDGENH
INTEGER(KIND=JPIM) :: IDLUN1, IDLUX1, IFLVM2, ILA, ILA1, ILA2, ILA3,&
 & ILAG, ILEV, ILEVV, ILO, JLAT, JLEV, JROF, IJ_, J_

REAL(KIND=JPRB) :: ZFAC

LOGICAL :: LLT_SLHD(3),LLSLHD,LLSLHDQUAD,LLSLHD_OLD,LLSLHDHEAT

INTEGER(KIND=JPIM) :: IADDR_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)    :: ZVETA_(JPDUP,0:KFLEV+1) 
REAL(KIND=JPRB)    :: ZVCUICO_(JPDUP,4,0:KFLEV-1)
REAL(KIND=JPRB)    :: ZVSLD_(JPDUP,3,0:KFLEV-1)
REAL(KIND=JPRB)    :: ZVSLDW_(JPDUP,3,3,0:KFLEV-1)

REAL(KIND=JPRB)    :: ZSLHDKMINH,ZSLHDKMINV


REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "lascaw_clo_tl.intfb.h"
#include "lascaw_vintw_tl.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ELASCAWTL',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       0.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

! cases relevant for SLHD scheme (switches LDSLHD, LDSLHDQUAD are
! deactivated during computation of medium points in LAPINEA)
LLSLHD=LDSLHD.AND.(KWIS==103.OR.KWIS==203)
LLSLHDQUAD=LDSLHDQUAD.AND.(KWIS==103.OR.KWIS==203)

! switch for old SLHD scheme
LLSLHD_OLD=LLSLHD.AND.LDSLHD_OLD

LLT_SLHD(1)=LLSLHD
LLT_SLHD(2)=LLSLHDQUAD
LLT_SLHD(3)=LLSLHD_OLD

LLSLHDHEAT = YDDYN%LSLHDHEAT

ZSLHDKMINH=SLHDKMIN
ZSLHDKMINV=SLHDKMIN

DO JLAT=YDSL%NDGSAH,YDSL%NDGENH
  IADDR(JLAT)=KSTABUF(JLAT)+YDSL%NASLB1*(0-KFLDN)
ENDDO

DO J_ = 1, JPDUP
  IADDR_(J_,YDSL%NDGSAH:YDSL%NDGENH)=IADDR(YDSL%NDGSAH:YDSL%NDGENH)
ENDDO
DO J_ = 1, JPDUP
  ZVETA_ (J_,0:KFLEV+1)  =PVETA (0:KFLEV+1)
  ZVCUICO_(J_,2,0:KFLEV-1)=PVCUICO(2,0:KFLEV-1)
  ZVCUICO_(J_,3,0:KFLEV-1)=PVCUICO(3,0:KFLEV-1)
  ZVCUICO_(J_,4,0:KFLEV-1)=PVCUICO(4,0:KFLEV-1)
ENDDO
IF (LLSLHDQUAD) THEN
  DO J_ = 1, JPDUP
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

IFLVM2=KFLEV-2

!     ----------------------------------------------------------------

!*       1.    3D MODEL.
!              ---------

IDLUN1=YDSL%NDLUNG-1
IDLUX1=YDSL%NDLUXG-1

! Constraint to enable halo fully inside E zone,
! while limiting halo on C+I whenever possible for tasks in E zone
IDGENH=MAX(YDSL%NDGSAH,MIN(YDSL%NDGUXG-YDSL%NFRSTLOFF,YDSL%NDGENH))

!        1.01  Coordinates and weights for trilinear interpolations.

IF (KWIS == 101) THEN

  ZFAC=PVRLEVX/(PVETA(KFLEV+1)-PVETA(0))

  DO JLEV=1,KFLEV

    ! * Calculation of linear weights, KL0.
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1
      ILAG=INT(PLAT5(JROF,JLEV))-1
      ILA=ILAG - YDSL%NFRSTLOFF
      PDLAT5(JROF,JLEV)=PLAT5(JROF,JLEV)-REAL(ILAG+1,JPRB)
      PDLAT(JROF,JLEV)=PLAT(JROF,JLEV)
      ILO=INT(PLON5(JROF,JLEV))-1
      PDLO5(JROF,JLEV,0)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO5(JROF,JLEV,1)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO5(JROF,JLEV,2)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO5(JROF,JLEV,3)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO(JROF,JLEV,0)=PLON(JROF,JLEV)
      PDLO(JROF,JLEV,1)=PLON(JROF,JLEV)
      PDLO(JROF,JLEV,2)=PLON(JROF,JLEV)
      PDLO(JROF,JLEV,3)=PLON(JROF,JLEV)

      ILEV  =KVAUT(INT(PLEV5(JROF,JLEV)*ZFAC))-1
      IF(ILEV < IFLVM2.AND.&
       & (PLEV5(JROF,JLEV)-ZVETA_(IJ_,ILEV+2)) > 0.0_JPRB) ILEV=ILEV+1  

      KLEV(JROF,JLEV)=ILEV

      PDVER5(JROF,JLEV)=(PLEV5(JROF,JLEV)-ZVETA_(IJ_,ILEV+1))/&
       & (ZVETA_(IJ_,ILEV+2)-ZVETA_(IJ_,ILEV+1))  
      PDVER(JROF,JLEV)=PLEV(JROF,JLEV)/(ZVETA_(IJ_,ILEV+2)-ZVETA_(IJ_,ILEV+1))  

      ILA1=ILA+1
      ILA1=MIN(MAX(ILA1,YDSL%NDGSAH),IDGENH)
      ILA2=ILA+2
      ILA2=MIN(MAX(ILA2,YDSL%NDGSAH),IDGENH)
      ILO =MIN(MAX(ILO,IDLUN1),IDLUX1)
      KL0(JROF,JLEV,1)=IADDR_(IJ_,ILA1)+YDSL%NSLEXT(ILO,ILA1)
      KL0(JROF,JLEV,2)=IADDR_(IJ_,ILA2)+YDSL%NSLEXT(ILO,ILA2)
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

IF (KWIS == 103 .OR. KWIS == 104 .OR. KWIS == 105) THEN

  ZFAC=PVRLEVX/(PVETA(KFLEV+1)-PVETA(0))

  DO JLEV=1,KFLEV
    IF (KHOR == 0) ILEV=JLEV
    IF (KHOR == 1) ILEV=KFLDN

    ! * Calculation of linear weights, KL0, KLH0.
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1
      ILAG=INT(PLAT5(JROF,JLEV))-1
      ILA=ILAG - YDSL%NFRSTLOFF
      PDLAT5(JROF,JLEV)=PLAT5(JROF,JLEV)-REAL(ILAG+1,JPRB)
      PDLAT(JROF,JLEV)=PLAT(JROF,JLEV)

      ILO=INT(PLON5(JROF,JLEV))-1
      PDLO5(JROF,JLEV,0)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO5(JROF,JLEV,1)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO5(JROF,JLEV,2)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO5(JROF,JLEV,3)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO(JROF,JLEV,0)=PLON(JROF,JLEV)
      PDLO(JROF,JLEV,1)=PLON(JROF,JLEV)
      PDLO(JROF,JLEV,2)=PLON(JROF,JLEV)
      PDLO(JROF,JLEV,3)=PLON(JROF,JLEV)

      ILEVV=KVAUT(INT(PLEV5(JROF,JLEV)*ZFAC))-1
      IF(ILEVV < IFLVM2.AND.&
       & (PLEV5(JROF,JLEV)-ZVETA_(IJ_,ILEVV+2)) > 0.0_JPRB) ILEVV=ILEVV+1  

      KLEV(JROF,JLEV)=ILEVV
      PDVER5(JROF,JLEV)=(PLEV5(JROF,JLEV)-ZVETA_(IJ_,ILEVV+1))/&
       & (ZVETA_(IJ_,ILEVV+2)-ZVETA_(IJ_,ILEVV+1))  
      PDVER(JROF,JLEV)=PLEV(JROF,JLEV)/(ZVETA_(IJ_,ILEVV+2)-ZVETA_(IJ_,ILEVV+1))  

      ILA1=ILA+1
      ILA1=MIN(MAX(ILA1,YDSL%NDGSAH),IDGENH)
      ILA2=ILA+2
      ILA2=MIN(MAX(ILA2,YDSL%NDGSAH),IDGENH)
      ILA3=ILA+3
      ILA3=MIN(MAX(ILA3,YDSL%NDGSAH),IDGENH)
      ILA =MIN(MAX(ILA ,YDSL%NDGSAH),IDGENH)
      ILO =MIN(MAX(ILO,IDLUN1),IDLUX1)
      KL0(JROF,JLEV,0)=IADDR_(IJ_,ILA )+YDSL%NSLEXT(ILO,ILA)
      KL0(JROF,JLEV,1)=IADDR_(IJ_,ILA1)+YDSL%NSLEXT(ILO,ILA1)
      KL0(JROF,JLEV,2)=IADDR_(IJ_,ILA2)+YDSL%NSLEXT(ILO,ILA2)
      KL0(JROF,JLEV,3)=IADDR_(IJ_,ILA3)+YDSL%NSLEXT(ILO,ILA3)

      KLH0(JROF,JLEV,0)=IADDR_(IJ_,ILA )+YDSL%NSLEXT(ILO,ILA)+YDSL%NASLB1*ILEV
      KLH0(JROF,JLEV,1)=IADDR_(IJ_,ILA1)+YDSL%NSLEXT(ILO,ILA1)+YDSL%NASLB1*ILEV
      KLH0(JROF,JLEV,2)=IADDR_(IJ_,ILA2)+YDSL%NSLEXT(ILO,ILA2)+YDSL%NASLB1*ILEV
      KLH0(JROF,JLEV,3)=IADDR_(IJ_,ILA3)+YDSL%NSLEXT(ILO,ILA3)+YDSL%NASLB1*ILEV
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
  CALL LASCAW_CLO_TL(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LLSLHDHEAT,&
   & ZSLHDKMINH,PDLAT,PKAPPA,PKAPPAT,PDLAT5,PKAPPA5,PKAPPAT5,&
   & PCLA,PCLASLD,PCLASLT,PCLA5,PCLASLD5,PCLASLT5)

  ! * Calculation of PCLO and PCLOSLD
  !   + calculation of PCLO5 and PCLOSLD5:
  CALL LASCAW_CLO_TL(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LLSLHDHEAT,&
   & ZSLHDKMINH,PDLO(:,:,1),PKAPPA,PKAPPAT,&
   & PDLO5(:,:,1),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,1),PCLOSLD(:,:,:,1),PCLOSLT(:,:,:,1),&
   & PCLO5(:,:,:,1),PCLOSLD5(:,:,:,1),PCLOSLT5(:,:,:,1))
  CALL LASCAW_CLO_TL(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LLSLHDHEAT,&
   & ZSLHDKMINH,PDLO(:,:,2),PKAPPA,PKAPPAT,&
   & PDLO5(:,:,2),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,2),PCLOSLD(:,:,:,2),PCLOSLT(:,:,:,2),&
   & PCLO5(:,:,:,2),PCLOSLD5(:,:,:,2),PCLOSLT5(:,:,:,2))
  ! * Calculation of PVINTW and PVINTWSLD
  !   + calculation of PVINTW5 and PVINTWSLD5:
  CALL LASCAW_VINTW_TL(KPROMB,KFLEV,KST,KPROF,LLT_SLHD,LLSLHDHEAT,ZSLHDKMINV, &
   & KLEV,PLEV,PDVER,PKAPPA,PKAPPAT,PVETA,ZVCUICO_,ZVSLD_,ZVSLDW_,&
   & PLEV5,PDVER5,PKAPPA5,PKAPPAT5,PVINTW,PVINTWSLD,PVINTWSLT,&
   & PVINTW5,PVINTWSLD5,PVINTWSLT5)

  IF (KWIS == 104) THEN !!! TL not supported !!!
    CALL ABOR1('ELASCAWTL: The option KWIS=104 is not supported by TL code')
  ENDIF

  IF (KWIS == 105) THEN !!! TL not supported !!!
    CALL ABOR1('ELASCAWTL: The option KWIS=105 is not supported by TL code')
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

      ILAG=INT(PLAT5(JROF,JLEV))-1
      ILA=ILAG - YDSL%NFRSTLOFF
      PDLAT5(JROF,JLEV)=PLAT5(JROF,JLEV)-REAL(ILAG+1,JPRB)
      PDLAT(JROF,JLEV)=PLAT(JROF,JLEV)
      ILO=INT(PLON5(JROF,JLEV))-1
      PDLO5(JROF,JLEV,1)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO5(JROF,JLEV,2)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO(JROF,JLEV,1)=PLON(JROF,JLEV)
      PDLO(JROF,JLEV,2)=PLON(JROF,JLEV)

      ILA1=ILA+1
      ILA1=MIN(MAX(ILA1,YDSL%NDGSAH),IDGENH)
      ILA2=ILA+2
      ILA2=MIN(MAX(ILA2,YDSL%NDGSAH),IDGENH)
      ILO =MIN(MAX(ILO,IDLUN1),IDLUX1)
      KL0(JROF,JLEV,1)=IADDR_(IJ_,ILA1)+YDSL%NSLEXT(ILO,ILA1)+YDSL%NASLB1*ILEV
      KL0(JROF,JLEV,2)=IADDR_(IJ_,ILA2)+YDSL%NSLEXT(ILO,ILA2)+YDSL%NASLB1*ILEV

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

IF (KWIS == 203) THEN

  DO JLEV=1,KFLEV
    IF (KHOR == 0) ILEV=JLEV
    IF (KHOR == 1) ILEV=KFLDN

    ! * Calculation of linear weights, KL0.
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1

      ILAG=INT(PLAT5(JROF,JLEV))-1
      ILA=ILAG - YDSL%NFRSTLOFF
      PDLAT5(JROF,JLEV)=PLAT5(JROF,JLEV)-REAL(ILAG+1,JPRB)
      PDLAT(JROF,JLEV)=PLAT(JROF,JLEV)

      ILO =INT(PLON5(JROF,JLEV))-1
      PDLO5(JROF,JLEV,0)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO5(JROF,JLEV,1)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO5(JROF,JLEV,2)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO5(JROF,JLEV,3)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO(JROF,JLEV,0)=PLON(JROF,JLEV)
      PDLO(JROF,JLEV,1)=PLON(JROF,JLEV)
      PDLO(JROF,JLEV,2)=PLON(JROF,JLEV)
      PDLO(JROF,JLEV,3)=PLON(JROF,JLEV)

      ILA1=ILA+1
      ILA1=MIN(MAX(ILA1,YDSL%NDGSAH),IDGENH)
      ILA2=ILA+2
      ILA2=MIN(MAX(ILA2,YDSL%NDGSAH),IDGENH)
      ILA3=ILA+3
      ILA3=MIN(MAX(ILA3,YDSL%NDGSAH),IDGENH)
      ILA =MIN(MAX(ILA ,YDSL%NDGSAH),IDGENH)
      ILO =MIN(MAX(ILO,IDLUN1),IDLUX1)
      KL0(JROF,JLEV,0)=IADDR_(IJ_,ILA )+YDSL%NSLEXT(ILO,ILA)
      KL0(JROF,JLEV,1)=IADDR_(IJ_,ILA1)+YDSL%NSLEXT(ILO,ILA1)
      KL0(JROF,JLEV,2)=IADDR_(IJ_,ILA2)+YDSL%NSLEXT(ILO,ILA2)
      KL0(JROF,JLEV,3)=IADDR_(IJ_,ILA3)+YDSL%NSLEXT(ILO,ILA3)

      KL0(JROF,JLEV,0)=KL0(JROF,JLEV,0)+YDSL%NASLB1*ILEV
      KL0(JROF,JLEV,1)=KL0(JROF,JLEV,1)+YDSL%NASLB1*ILEV
      KL0(JROF,JLEV,2)=KL0(JROF,JLEV,2)+YDSL%NASLB1*ILEV
      KL0(JROF,JLEV,3)=KL0(JROF,JLEV,3)+YDSL%NASLB1*ILEV
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
  CALL LASCAW_CLO_TL(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LLSLHDHEAT,&
   & ZSLHDKMINH,PDLAT,PKAPPA,PKAPPAT,PDLAT5,PKAPPA5,PKAPPAT5,&
   & PCLA,PCLASLD,PCLASLT,PCLA5,PCLASLD5,PCLASLT5)

  ! * Calculation of PCLO and PCLOSLD
  !   + calculation of PCLO5 and PCLOSLD5:
  CALL LASCAW_CLO_TL(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LLSLHDHEAT,&
   & ZSLHDKMINH,PDLO(:,:,1),PKAPPA,PKAPPAT,&
   & PDLO5(:,:,1),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,1),PCLOSLD(:,:,:,1),PCLOSLT(:,:,:,1),&
   & PCLO5(:,:,:,1),PCLOSLD5(:,:,:,1),PCLOSLT5(:,:,:,1))
  CALL LASCAW_CLO_TL(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,LLSLHDHEAT,&
   & ZSLHDKMINH,PDLO(:,:,2),PKAPPA,PKAPPAT,&
   & PDLO5(:,:,2),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,2),PCLOSLD(:,:,:,2),PCLOSLT(:,:,:,2),&
   & PCLO5(:,:,:,2),PCLOSLD5(:,:,:,2),PCLOSLT5(:,:,:,2))


ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ELASCAWTL',1,ZHOOK_HANDLE)
END SUBROUTINE ELASCAWTL
