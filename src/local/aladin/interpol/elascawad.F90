!ocl list_copy(32,PVETA)
!ocl list_copy(32,PVCUICO)
SUBROUTINE ELASCAWAD(&
 ! --- INPUT -------------------------------------------------
 & YDDYN,YDSL,KPROMB,KST,KPROF,KFLEV,&
 & KFLDN,KWIS,KHOR,&
 & LDSLHD,LDSLHDQUAD,LDSLHD_OLD,&
 ! --- OUTPUT (input in TL) ----------------------------------
 & PLON,PLAT,PLEV,PKAPPA,PKAPPAT,&
 ! --- INPUT -------------------------------------------------
 & PLON5,PLAT5,PLEV5,PKAPPA5,PKAPPAT5,&
 & PVETA,KVAUT,&
 & PVCUICO,PVSLD,PVSLDW,KRLEVX,PVRLEVX,&
 ! --- INPUT/OUTPUT (output in TL) ---------------------------
 & PDLAT,PCLA,PCLASLD,PCLASLT,PDLO,PCLO,PCLOSLD,PCLOSLT,PDVER,PVINTW,PVINTWSLD,PVINTWSLT)

!**** *ELASCAWAD  -  Externalisable interpolator:   (adjoint code)
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
!        *CALL* *ELASCAWAD( ... )

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
!          KWIS    - kind of interpolation.
!          KHOR    - 0: Horizontal interpolation for each level
!                       of a 3D variable.
!                    1: Interpolation for all origin points corresponding
!                       to a final point of a 2D variable.
!          LDSLHD  - key activating SLHD weights precomputation.
!          LDSLHDQUAD - key activating quadratic weights precomputation.
!          LDSLHD_OLD - use old SLHD interpolator
!  ------------------- Trajectory variables ------------------------
!          PLON5   - x-coordinate of the interpolation point
!                    (in the fractional system <NDLUNG;NDLON>)
!          PLAT5   - y-coordinate of the interpolation point
!                    (in the fractional system <NDGSAG;NDGENG>)
!          PLEV5   - vertical coordinate of the interpolation point.
!          PKAPPA5 - kappa function ("coefficient of SLHD") based on the
!                    rescaled horizontal deformation of the flow evaluated
!                    at instant "t" for the final point F
!  -----------------------------------------------------------------
!          PVETA   - Values of ETA.
!          KVAUT   - Help table for vertical box search: gives the number
!                    of the layer or interlayer immediately above eta.
!          PVCUICO - Denominators of the vertical cubic interpolation coef.
!          PVSLD   - auxiliary quantities for vertical SLHD interpolation
!          PVSLDW  - weights for SLHD vertical Laplacian smoother
!          KRLEVX  - Dimension of KVAUT.
!          PVRLEVX - REAL(KRLEVX).

! ----------------- INOUT (OUTPUT of TL) --------------------------
!          PDLAT     - distance for horizontal linear interpolations in latitude.
!          PCLA      - weights for horizontal cubic interpolations in latitude.
!          PCLASLD   - weights for horizontal cubic interpolations in latitude, SLHD case
!          PDLO      - distances for horizontal linear interpolations
!                      in longitude (latitude rows 0, 1, 2, 3)
!          PCLO      - weights for horizontal cubic interpolations in
!                      longitude (latitude rows 1, 2)
!          PCLOSLD   - weights for horizontal cubic interpolations in
!                      longitude, SLHD case (latitude rows 1, 2)
!          PDVER     - distance for vertical linear interpolation
!          PVINTW    - vertical cubic interpolation weights
!          PVINTWSLD - vertical cubic interpolation weights, SLHD case
! ----------------------------------------------------------------
!        OUTPUT: (input of TL)
!          PLON    - x-coordinate of the interpolation point
!                    (in the fractional system <NDLUNG;NDLON>)
!          PLAT    - y-coordinate of the interpolation point
!                    (in the fractional system <NDGSAG;NDGENG>)
!          PLEV    - vertical coordinate of the interpolation point.
!          PKAPPA  - kappa function ("coefficient of SLHD") based on the
!                    rescaled horizontal deformation of the flow evaluated
!                    at instant "t" for the final point F

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
!        F. VANA, being inspired by LASCAWTL/AD
!        * ALADIN/LACE *
!        Original : SEPTEMBER 2006.

!     Modifications.
!     --------------
!        02-Jun-2008 A. Bogatchev - USE YOMLASCAW et c.
!        F. Vana  23-Jul-2008 : weights driven interpolation (including SLHD)
!        K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!        K. Yessad + F. Vana (Feb 2009): split loops, rewrite in a shorter way.
!        B. Bochenek (Apr 2015): Phasing: move some variables and update
!        K. Yessad (March 2017): simplify level numbering.
!     ------------------------------------------------------------------

USE YOMDYN    ,ONLY : TDYN
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

! arp/ifs dependencies to be solved later.
USE YOMMP0   , ONLY : LOPT_SCALAR
USE YOMDYNA  , ONLY : SLHDKMIN
USE EINT_MOD , ONLY : SL_STRUCT, JPDUP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(SL_STRUCT),   INTENT(IN)    :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMB 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KRLEVX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KWIS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KHOR 
LOGICAL           ,INTENT(IN)    :: LDSLHD
LOGICAL           ,INTENT(IN)    :: LDSLHDQUAD
LOGICAL           ,INTENT(IN)    :: LDSLHD_OLD
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
REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PCLO(KPROMB,KFLEV,3,2)
REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PCLOSLD(KPROMB,KFLEV,3,2)
REAL(KIND=JPRB)   ,INTENT(INOUT)    :: PCLOSLT(KPROMB,KFLEV,3,2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDVER(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVINTW(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVINTWSLD(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVINTWSLT(KPROMB,KFLEV,3)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IDLUN1, IDLUX1, IFLVM2, ILA, ILO,&
 & ILAG, ILEV, ILEVV, JLEV, JROF, IJ_, J_
INTEGER(KIND=JPIM) :: IILEV(KPROMB,KFLEV)

REAL(KIND=JPRB) :: ZFAC
REAL(KIND=JPRB) :: ZDLAT5(KPROMB,KFLEV),ZDLO5(KPROMB,KFLEV),ZDVER5(KPROMB,KFLEV)

LOGICAL :: LLT_SLHD(3),LLSLHD,LLSLHDQUAD,LLSLHD_OLD,LLSLHDHEAT

REAL(KIND=JPRB)    :: ZVETA_(JPDUP,0:KFLEV+1) 
REAL(KIND=JPRB)    :: ZVCUICO_(JPDUP,4,0:KFLEV-1)
REAL(KIND=JPRB)    :: ZVSLD_(JPDUP,3,0:KFLEV-1)
REAL(KIND=JPRB)    :: ZVSLDW_(JPDUP,3,3,0:KFLEV-1)

REAL(KIND=JPRB)    :: ZSLHDKMINH,ZSLHDKMINV

REAL(KIND=JPRB)    :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "lascaw_clo_ad.intfb.h"
#include "lascaw_vintw_ad.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ELASCAWAD',0,ZHOOK_HANDLE)

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
ZSLHDKMINH = SLHDKMIN
ZSLHDKMINV = SLHDKMIN

DO J_ = 1, JPDUP 
  ZVETA_ (J_,  0:KFLEV+1)=PVETA (0:KFLEV+1)
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
IF (LLSLHD) THEN
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

!        1.01  Coordinates and weights for trilinear interpolations.

IF (KWIS == 101) THEN

  ZFAC=PVRLEVX/(PVETA(KFLEV+1)-PVETA(0))

  DO JLEV=1,KFLEV

    ! * AD calculation of linear weights.
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1
      ILAG=INT(PLAT5(JROF,JLEV))-1
      ILA=ILAG - YDSL%NFRSTLOFF
      
      PLAT(JROF,JLEV)=PLAT(JROF,JLEV)+PDLAT(JROF,JLEV)

      PLON(JROF,JLEV)=PLON(JROF,JLEV)+PDLO(JROF,JLEV,0)+&
       & PDLO(JROF,JLEV,1)+PDLO(JROF,JLEV,2)+PDLO(JROF,JLEV,3)

      ILEV  =KVAUT(INT(PLEV5(JROF,JLEV)*ZFAC))-1
      IF(ILEV < IFLVM2.AND.&
       & (PLEV5(JROF,JLEV)-ZVETA_(IJ_,ILEV+2)) > 0.0_JPRB) ILEV=ILEV+1  

      PLEV(JROF,JLEV) = PLEV(JROF,JLEV) +&
       & PDVER(JROF,JLEV)/(ZVETA_(IJ_,ILEV+2)-ZVETA_(IJ_,ILEV+1))  

    ENDDO
  ENDDO

ENDIF

!        1.03  Coordinates and weights for ( horizontal 12 points
!              + vertical cubic + 32 points interpolations ) or
!              ( horizontal 12 points + 32 points interpolations ).

IF (KWIS == 103 .OR. KWIS == 104 .OR. KWIS == 105) THEN

  ZFAC=PVRLEVX/(PVETA(KFLEV+1)-PVETA(0))

  ! Following loop has been splitted for efficiency reasons (mainly
  ! to ease vectorization). The SLHD case is more complex and requires
  ! extra trajectories computation.

  DO JLEV=1,KFLEV
    IF (KHOR == 0) ILEV=JLEV
    IF (KHOR == 1) ILEV=KFLDN

    ! * AD of calculation of linear weights: trajectory.
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1
      ILAG=INT(PLAT5(JROF,JLEV))-1
      ILA=ILAG - YDSL%NFRSTLOFF
      ZDLAT5(JROF,JLEV)=PLAT5(JROF,JLEV)-REAL(ILAG+1,JPRB)

      ILO=INT(PLON5(JROF,JLEV))-1
      ZDLO5(JROF,JLEV)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)

      ILEVV=KVAUT(INT(PLEV5(JROF,JLEV)*ZFAC))-1
      IF(ILEVV < IFLVM2.AND.&
       & (PLEV5(JROF,JLEV)-ZVETA_(IJ_,ILEVV+2)) > 0.0_JPRB) ILEVV=ILEVV+1
      IILEV(JROF,JLEV)=ILEVV

      ZDVER5(JROF,JLEV)=(PLEV5(JROF,JLEV)-ZVETA_(IJ_,ILEVV+1))/&
       & (ZVETA_(IJ_,ILEVV+2)-ZVETA_(IJ_,ILEVV+1))

    ENDDO
  ENDDO

  ! * AD of calculation of PVINTW and PVINTWSLD:
  CALL LASCAW_VINTW_AD(KPROMB,KFLEV,KST,KPROF,LLT_SLHD,LLSLHDHEAT,ZSLHDKMINV,&
   & IILEV,PLEV,PDVER,PKAPPA,PKAPPAT,PVETA,ZVCUICO_,ZVSLD_,ZVSLDW_,&
   & PLEV5,ZDVER5,PKAPPA5,PKAPPAT5,&
   & PVINTW,PVINTWSLD,PVINTWSLT)

    ! * AD of calculation of PCLO and PCLOSLD:
  CALL LASCAW_CLO_AD(KFLEV,KPROMB,KST,KPROF,LLT_SLHD, &
   & LLSLHDHEAT,ZSLHDKMINH,PDLO(:,:,2),PKAPPA,PKAPPAT,&
   & ZDLO5(:,:),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,2),PCLOSLD(:,:,:,2),PCLOSLT(:,:,:,2))
  CALL LASCAW_CLO_AD(KFLEV,KPROMB,KST,KPROF,LLT_SLHD, &
   & LLSLHDHEAT,ZSLHDKMINH,PDLO(:,:,1),PKAPPA,PKAPPAT,&
   & ZDLO5(:,:),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,1),PCLOSLD(:,:,:,1),PCLOSLT(:,:,:,1))

  ! * AD of calculation of PCLA and PCLASLD:
  CALL LASCAW_CLO_AD(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,  &
   & LLSLHDHEAT,ZSLHDKMINH,PDLAT,PKAPPA,PKAPPAT,ZDLAT5,&
   & PKAPPA5,PKAPPAT5,PCLA,PCLASLD,PCLASLT)

  ! * AD of calculation of linear weights: increment.
  DO JLEV=1,KFLEV
!CDIR NODEP
    DO JROF=KST,KPROF

      IJ_ = MOD(JROF+1-KST,JPDUP)+1
      PLEV(JROF,JLEV)=PLEV(JROF,JLEV)+PDVER(JROF,JLEV)/&
       & (ZVETA_(IJ_,IILEV(JROF,JLEV)+2)-ZVETA_(IJ_,IILEV(JROF,JLEV)+1))  

      PLON(JROF,JLEV)=PLON(JROF,JLEV)+PDLO(JROF,JLEV,0)+&
       & PDLO(JROF,JLEV,1)+PDLO(JROF,JLEV,2)+PDLO(JROF,JLEV,3)

      PLAT(JROF,JLEV)=PLAT(JROF,JLEV)+PDLAT(JROF,JLEV)

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

      ILAG=INT(PLAT5(JROF,JLEV))-1
      ILA=ILAG - YDSL%NFRSTLOFF
      PLAT(JROF,JLEV)=PLAT(JROF,JLEV)+PDLAT(JROF,JLEV)

      PLON(JROF,JLEV)=PLON(JROF,JLEV)+PDLO(JROF,JLEV,1)+PDLO(JROF,JLEV,2)

    ENDDO
  ENDDO

ENDIF

!        2.03  Coordinates and weights for 12 points interpolations.

IF (KWIS == 203) THEN

  ! Following loop has been splitted for efficiency reasons (mainly
  ! to ease vectorization). The SLHD case is more complex and requires
  ! extra trajectories computation.

  DO JLEV=1,KFLEV
    IF (KHOR == 0) ILEV=JLEV
    IF (KHOR == 1) ILEV=KFLDN

    ! * AD of calculation of linear weights: trajectory.
!CDIR NODEP
    DO JROF=KST,KPROF

      ILAG=INT(PLAT5(JROF,JLEV))-1
      ILA=ILAG - YDSL%NFRSTLOFF
      ZDLAT5(JROF,JLEV)=PLAT5(JROF,JLEV)-REAL(ILAG+1,JPRB)

      ILO=INT(PLON5(JROF,JLEV))-1
      ZDLO5(JROF,JLEV)=PLON5(JROF,JLEV)-REAL(ILO+1,JPRB)

    ENDDO
  ENDDO

  ! * AD of calculation of PCLO and PCLOSLD:
  CALL LASCAW_CLO_AD(KFLEV,KPROMB,KST,KPROF,LLT_SLHD, &
   & LLSLHDHEAT,ZSLHDKMINH,PDLO(:,:,2),PKAPPA,PKAPPAT,&
   & ZDLO5(:,:),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,2),PCLOSLD(:,:,:,2),PCLOSLT(:,:,:,2))
  CALL LASCAW_CLO_AD(KFLEV,KPROMB,KST,KPROF,LLT_SLHD, &
   & LLSLHDHEAT,ZSLHDKMINH,PDLO(:,:,1),PKAPPA,PKAPPAT,&
   & ZDLO5(:,:),PKAPPA5,PKAPPAT5,&
   & PCLO(:,:,:,1),PCLOSLD(:,:,:,1),PCLOSLT(:,:,:,1))

  ! * AD of calculation of PCLA and PCLASLD:
  CALL LASCAW_CLO_AD(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,  &
   & LLSLHDHEAT,ZSLHDKMINH,PDLAT,PKAPPA,PKAPPAT,ZDLAT5,&
   & PKAPPA5,PKAPPAT5,PCLA,PCLASLD,PCLASLT)

  ! * AD of calculation of linear weights: increment.
  DO JLEV=1,KFLEV
!CDIR NODEP
    DO JROF=KST,KPROF

      PLON(JROF,JLEV)=PLON(JROF,JLEV)+PDLO(JROF,JLEV,0)+&
       & PDLO(JROF,JLEV,1)+PDLO(JROF,JLEV,2)+PDLO(JROF,JLEV,3)

      PLAT(JROF,JLEV)=PLAT(JROF,JLEV)+PDLAT(JROF,JLEV)

    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ELASCAWAD',1,ZHOOK_HANDLE)
END SUBROUTINE ELASCAWAD
