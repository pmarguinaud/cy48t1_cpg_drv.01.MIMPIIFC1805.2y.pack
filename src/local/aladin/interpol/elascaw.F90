!option! -O extendreorder
SUBROUTINE ELASCAW(&
 ! --- INPUT -------------------------------------------------
 & YDVSPLIP,LDSLHDHEAT,YDSL,KPROMB,KDIMK,KST,KPROF,KFLEV,&
 & KFLDN,KSTABUF,KWIS,KHOR,KWENO,KHVI,&
 & LDSLHD,LDSLHDQUAD,LDSLHD_OLD,LD3DTURB,&
 & LDCOMAD,LDCOMADH,LDCOMADV,KSPLTHOI,&
 & PLON,PLAT,PLEV,&
 & PVETA,KVAUT,&
 & PVCUICO,PVSLD,PVSLDW,PGAMMA_WENO,KRLEVX,PVRLEVX,&
 & PKAPPA,PKAPPAT,PKAPPAM,PKAPPAH,&
 & PSTDDISU,PSTDDISV,PSTDDISW,&
 ! --- OUTPUT ------------------------------------------------
 & PDLAT,PDLAMAD,PCLA,PCLASLD,PCLASLT,PCLAMAD,&
 & PDLO ,PDLOMAD,PCLO,PCLOSLD,PCLOSLT,PCLOMAD,&
 & KL0,KLH0,KLEV,KNOWENO,PCW,&
 & PDVER,PDVERMAD,PVINTW,PVINTWSLD,PVINTWSLT,PVINTWMAD,PVINTWS,&
 & PVDERW,PHVW,KDEP)

!**** *ELASCAW  -  Externalisable interpolator:
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

!       Note that this routine should not know if levels are half levels or full levels;
!       this information must remain in the caller.

!**   Interface.
!     ----------
!        *CALL* *ELASCAW( ... )

!        Explicit arguments :
!        --------------------

!        INPUT:
!          YDSL    - SL_STRUCT definition
!          YDDYN   - structure containing dynamics.
!          KPROMB  - horizontal dimension for interpolation point quantities.
!          KDIMK   - last dimension for some non-linear weights.
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
!          KHVI    - 1/0: filling weights arrays PVDERW and PHVW is necessary/not necessary.
!          LDSLHD  - key activating SLHD weights precomputation
!          LDSLHDQUAD - key activating quadratic weights precomputation
!          LDSLHD_OLD - use old SLHD interpolator
!          LD3DTURB- key activating 3D turbulence weights precomputation
!          LDCOMAD -  key activating COMAD weight computation
!          LDCOMADH-  key activating hor. COMAD
!          LDCOMADV-  key activating ver. COMAD
!          KSPLTHOI- controls additional weights precomputation
!          PLON    - x-coordinate of the interpolation point
!                    (in the fractional system <NDLUNG;NDLON>)
!          PLAT    - y-coordinate of the interpolation point
!                    (in the fractional system <NDGSAG;NDGENG>)
!          PLEV    - vertical coordinate of the interpolation point.
!          PVETA   - Values of ETA.
!          KVAUT   - Help table for vertical box search: gives the number
!                    of the layer immediately above eta.
!          PVCUICO - Denominators of the vertical cubic interpolation coefficients
!          PVSLD   - auxiliary quantities for vertical SLHD interpolation
!          PVSLDW  - weights for SLHD vertical Laplacian smoother
!          KRLEVX  - Dimension of KVAUT
!          PVRLEVX - REAL(KRLEVX).
!          PKAPPA  - kappa function ("coefficient of SLHD") based on the
!                    rescaled horizontal deformation of the flow evaluated
!                    at instant "t" for the final point F
!          PKAPPAT - PKAPPA for heat variable
!          PKAPPAM - horizontal exchange coefficient for momentum in 3D turb. 
!          PKAPPAH - horizontal exchange coefficient for heat in 3D turb.
!          PSTDDISU- COMAD correction coefficient based on estimated flow deformation
!                    along zonal direction of the trajectory but computed
!                    with wind derivatives at instant "t" at the final point F
!          PSTDDISV- COMAD correction coefficient based on estimated flow deformation
!                    along meridional direction of the trajectory but computed
!                    with wind derivatives at instant "t" at the final point F
!          PSTDDISW- COMAD correction coefficient based on estimated flow deformation
!                    along vertical direction of the trajectory but computed
!                    with wind derivatives at instant "t" at the final point F
!        OUTPUT:
!          PDLAT     - distance for horizontal linear interpolations in latitude
!          PDLAMAD   - PDLAT, COMAD case
!          PCLA      - weights for horizontal cubic interpolations in latitude
!          PCLASLD   - weights for horizontal cubic interpolations in latitude, SLHD case
!          PCLAMAD   - PCLA, COMAD case
!          PCLASLT   - weights for horizontal cubic interpolations in latitude, SLHD case on T
!          PDLO      - distances for horizontal linear interpolations
!                      in longitude (latitude rows 0, 1, 2, 3)
!          PDLOMAD   - PDLO, COMAD case
!          PCLO      - weights for horizontal cubic interpolations in
!                      longitude (latitude rows 1, 2)
!          PCLOSLD   - weights for horizontal cubic interpolations in
!                      longitude, SLHD case (latitude rows 1, 2)
!          PCLOMAD   - PCLO, COMAD case
!          PCLOSLT   - weights for horizontal cubic interpolations in
!                      longitude, SLHD case (latitude rows 1, 2) on T
!          KL0       - index of the four western points
!                      of the 16 points interpolation grid.
!          KLH0      - second value of index of the four western points
!                      of the 16 points interpolation grid if needed.
!          KLEV      - lower level of the vertical interpolation
!                      grid needed for vertical interpolations.
!          PDVER     - distance for vertical linear interpolation
!          PDVERMAD  - PDVER, COMAD case
!          PVINTW    - vertical cubic interpolation weights
!          PVINTWSLD - vertical cubic interpolation weights, SLHD case
!          PVINTWMAD - PVINTW, COMAD case
!          PVINTWSLT - vertical cubic interpolation weights, SLHD case on T
!          PVINTWS   - Vertical spline interpolation weights.
!          PVDERW    - weights to compute vertical derivatives necessary for
!                      Hermite cubic vertical interpolation.
!          PHVW      - Hermite vertical cubic interpolation weights.
!          KDEP      - indication of the interpolation stencil latitudial
!                      dependences (used for LVECADIN=.T. option in adjoint)

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
!        M. JANOUSEK, after the subroutine LASCAW written by K. YESSAD
!        METEO-FRANCE, CNRM/GMAP/EXT
!        Original : MARCH 1994.

!     Modifications.
!     --------------
!        Modified 28-Aug-2007 F. Vana : removing distances for 4p-splines
!        07-Nov-2007 J. Masek   New weights for SLHD interpolators.
!        Modified 02-Jun-2008 A. Bogatchev - use YOMLASCAW et c.
!        F. Vana  23-Jul-2008 support for vectorization
!        K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!        K. Yessad (Feb 2009): split loops, rewrite in a shorter way.
!        R. El Khatib 07-08-2009 Optimisation directive for NEC
!        F. Vana  23-Feb-2011: Horizontal turbulence + diff of phys. tend.
!        G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!        Y. Seity  20-Jan-2012: activation of SLHD in Preditor in case LPC_CHEAP
!        S. Malardel and D. Ricard (Nov 2013): COMAD weights for SL interpolations
!        B. Bochenek (Apr 2015): Phasing: move some variables and update
!        K. Yessad (March 2017): simplify level numbering.
!        R. El Khatib 07-May-2019 Remove bounds violations
!        R. El Khatib 22-May-2019 fix the "IFL OUT OF BOUNDS" issue
!        H Petithomme (Dec 2020): optimisation, COMAD bugfix
! End Modifications
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM,JPIA,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

! arp/ifs dependencies to be solved later.
USE YOMCT0    , ONLY : LREGETA
USE YOMDYNA  , ONLY : YRDYNA
USE YOMMP0   , ONLY : LOPT_SCALAR,NPROC

USE EINT_MOD , ONLY : SL_STRUCT,JPDUP
USE YOMVSPLIP , ONLY : TVSPLIP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVSPLIP)     ,INTENT(IN)    :: YDVSPLIP
LOGICAL           ,INTENT(IN)    :: LDSLHDHEAT
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMB
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIMK
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTABUF(YDSL%NDGSAH:YDSL%NDGENH)
INTEGER(KIND=JPIM),INTENT(IN)    :: KWIS
INTEGER(KIND=JPIM),INTENT(IN)    :: KHOR
INTEGER(KIND=JPIM),INTENT(IN)    :: KWENO
INTEGER(KIND=JPIM),INTENT(IN)    :: KHVI
LOGICAL           ,INTENT(IN)    :: LDSLHD
LOGICAL           ,INTENT(IN)    :: LDSLHDQUAD
LOGICAL           ,INTENT(IN)    :: LDSLHD_OLD
LOGICAL           ,INTENT(IN)    :: LD3DTURB
LOGICAL           ,INTENT(IN)    :: LDCOMAD
LOGICAL           ,INTENT(IN)    :: LDCOMADH
LOGICAL           ,INTENT(IN)    :: LDCOMADV
INTEGER(KIND=JPIM),INTENT(IN)    :: KSPLTHOI
INTEGER(KIND=JPIM),INTENT(IN)    :: KRLEVX
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLON(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLEV(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVETA(0:KFLEV+1)
INTEGER(KIND=JPIM),INTENT(IN)    :: KVAUT(0:KRLEVX)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVCUICO(4,0:KFLEV-1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSLD(3,0:KFLEV-1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVSLDW(3,3,0:KFLEV-1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAMMA_WENO(KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRLEVX
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPPA(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPPAT(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPPAM(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPPAH(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTDDISU(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTDDISV(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSTDDISW(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLAT(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLAMAD(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLA(KPROMB,KFLEV,3,KDIMK)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLASLD(KPROMB,KFLEV,3,KDIMK)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLASLT(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLAMAD(KPROMB,KFLEV,3,KDIMK)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLO(KPROMB,KFLEV,0:3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDLOMAD(KPROMB,KFLEV,0:3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLO(KPROMB,KFLEV,3,2,KDIMK)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLOSLD(KPROMB,KFLEV,3,2,KDIMK)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLOSLT(KPROMB,KFLEV,3,2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLOMAD(KPROMB,KFLEV,3,2,KDIMK)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KL0(KPROMB,KFLEV,0:3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLH0(KPROMB,KFLEV,0:3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLEV(KPROMB,KFLEV)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KNOWENO(KPROMB,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCW(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDVER(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDVERMAD(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTW(KPROMB,KFLEV,3*KWENO)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTWSLD(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTWMAD(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTWSLT(KPROMB,KFLEV,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVINTWS(KPROMB,KFLEV,1:4)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVDERW(KPROMB,KFLEV,2*KHVI,2*KHVI)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHVW(KPROMB,KFLEV,4*KHVI)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KDEP(KPROMB,KFLEV)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIA) :: ILEV64
INTEGER(KIND=JPIM) :: IDGENH,IOFF,JOFF
INTEGER(KIND=JPIM) :: IDLUN1, IDLUX1, IFLVM2, ILA, ILA1, ILA2, ILA3,&
 & ILAG, ILEV, ILO, JLAT, JLEV, JROF, JJ, IJ_
INTEGER(KIND=JPIM) :: IILEV(KPROMB,KFLEV),IBCLIM

REAL(KIND=JPRB) :: PD, ZD2, ZDVER, ZFAC, ZNUM,ZW1,ZW2,Z1,Z2,Z3,Z4
REAL(KIND=JPRB) :: ZKHTURB(KPROMB,KFLEV,KDIMK)
LOGICAL         :: LLT_SLHD(4),LLT_PHYS(4),LLSLHD,LLSLHDQUAD,LLSLHD_OLD,LL3DTURB
LOGICAL         :: LLCOMAD,LLCOMADH,LLCOMADV
REAL(KIND=JPRB) :: ZCLA(KPROMB,KFLEV,3) 
REAL(KIND=JPRB) :: ZCLO(KPROMB,KFLEV,3,2)
REAL(KIND=JPRB) :: ZSLHDKMINH,ZSLHDKMINV
REAL(KIND=JPRB), PARAMETER :: ZSLHDKMINV_WENO=0._JPRB ! WENO only runs with Lagrangian cubic!!!

! Duplicata of some dummy or local arrays for optimisation on NEC platform.
INTEGER(KIND=JPIM) :: IADDR_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB) :: ZVETA_(JPDUP,0:KFLEV+1)
REAL(KIND=JPRB) :: ZVCUICO_(JPDUP,4,0:KFLEV-1)
REAL(KIND=JPRB) :: ZVSLD_(JPDUP,3,0:KFLEV-1)
REAL(KIND=JPRB) :: ZVSLDW_(JPDUP,3,3,0:KFLEV-1)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
! functions

REAL(KIND=JPRB) :: FHLO1, FHLO2, FHLO3, FHLO4

! auxiliary functions for Hermite cubic interpolation
FHLO1(PD)= (1.0_JPRB-PD)*(1.0_JPRB-PD)*(1.0_JPRB+2.0_JPRB*PD)
FHLO2(PD)= PD*PD*(3._JPRB-2.0_JPRB*PD)
FHLO3(PD)= PD*(1.0_JPRB-PD)*(1.0_JPRB-PD)
FHLO4(PD)=-PD*PD*(1.0_JPRB-PD)

!     ------------------------------------------------------------------

#include "lascaw_clo.intfb.h"
#include "lascaw_cloturb.intfb.h"
#include "lascaw_vintw.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ELASCAW',0,ZHOOK_HANDLE)
ASSOCIATE(RFVV=>YDVSPLIP%RFVV)
!     ------------------------------------------------------------------

!*       0.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

! cases relevant for SLHD scheme (switches LDSLHD, LDSLHDQUAD are
! deactivated during computation of medium points in LAPINEA)
LLSLHD=LDSLHD.AND.(KWIS==103.OR.KWIS==104.OR.KWIS==105.OR.KWIS==106.OR.KWIS==203)
LLSLHDQUAD=LDSLHDQUAD.AND.(KWIS==103.OR.KWIS==104.OR.KWIS==105.OR.KWIS==106.OR.KWIS==203)
LL3DTURB=LD3DTURB.AND.(KWIS==103.OR.KWIS==104.OR.KWIS==105.OR.KWIS==106)

! switch for old SLHD scheme
LLSLHD_OLD=LLSLHD.AND.LDSLHD_OLD

LLT_SLHD(1)=LLSLHD
LLT_SLHD(2)=LLSLHDQUAD
LLT_SLHD(3)=LLSLHD_OLD
LLT_SLHD(4)=.FALSE.

ZSLHDKMINH=YRDYNA%SLHDKMIN
ZSLHDKMINV=YRDYNA%SLHDKMIN

! cases relevant for COMAD scheme (switches LDCOMADH and LDCOMADV  are
! deactivated during computation of interpolation points in LAPINEA)
LLCOMAD =LDCOMAD.AND.(KWIS==103.OR.KWIS==104.OR.KWIS==105.OR.KWIS==106.OR.KWIS==203)
LLCOMADH=LLCOMAD.AND.LDCOMADH
LLCOMADV=LLCOMAD.AND.LDCOMADV

! switches for interpolation of physics 
! It holds the same value for every iteration step (in ICI scheme). 
LLT_PHYS(1)=YRDYNA%LSLHD
LLT_PHYS(2)=YRDYNA%LSLHDQUAD
LLT_PHYS(3)=LDSLHD_OLD
LLT_PHYS(4)=.FALSE.

IF (KWIS==106) THEN
  ! Just make sure there is no change from default 3rd order Lagrangian cubic
  LLT_PHYS(2)=.FALSE.
ENDIF 

DO JJ = 1, JPDUP
  IADDR_(JJ,YDSL%NDGSAH:YDSL%NDGENH)=KSTABUF(YDSL%NDGSAH:YDSL%NDGENH)-YDSL%NASLB1*KFLDN
ENDDO
DO JJ = 1, JPDUP
  ZVETA_(JJ,:)=PVETA (:)
!CDIR NOUNROLL
  ZVCUICO_(JJ,:,:)=PVCUICO(:,:)
ENDDO
IF (LLSLHDQUAD) THEN
  DO JJ = 1, JPDUP
    IF (LOPT_SCALAR) THEN
      ZVSLD_(JJ,1:3,0:KFLEV-1)=PVSLD(1:3,0:KFLEV-1)
    ELSE
!CDIR NOUNROLL
      DO JLEV=1,3*KFLEV
        ZVSLD_(JJ,JLEV,0)=PVSLD(JLEV,0)
      ENDDO
    ENDIF
  ENDDO
ENDIF
IF (LLSLHD) THEN
  DO JJ = 1, JPDUP
    IF (LOPT_SCALAR) THEN
      ZVSLDW_(JJ,1:3,1:3,0:KFLEV-1)=PVSLDW(1:3,1:3,0:KFLEV-1)
    ELSE
!CDIR NOUNROLL
      DO JLEV=1,9*KFLEV
        ZVSLDW_(JJ,JLEV,1,0)=PVSLDW(JLEV,1,0)
      ENDDO
    ENDIF
  ENDDO
ENDIF
IF (LL3DTURB) THEN
  ! copy from kappa
  ZKHTURB(KST:KPROF,1:KFLEV,2)=PKAPPAM(KST:KPROF,1:KFLEV)
  ZKHTURB(KST:KPROF,1:KFLEV,3)=PKAPPAH(KST:KPROF,1:KFLEV)
ENDIF

IF (KSPLTHOI == 1) THEN
  ! In this case the ZKHTURB is set to static mode with maximum diffusion.
  ZKHTURB(KST:KPROF,1:KFLEV,KDIMK)=1._JPRB
ENDIF

IFLVM2=KFLEV-2

!     ------------------------------------------------------------------

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
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      IJ_ = MOD(JROF+1-KST,JPDUP)+1
      ILAG=INT(PLAT(JROF,JLEV))-1
      ILA=ILAG - YDSL%NFRSTLOFF
      PDLAT(JROF,JLEV)=PLAT(JROF,JLEV)-REAL(ILAG+1,JPRB)

      ILO=INT(PLON(JROF,JLEV))-1
      PDLO(JROF,JLEV,1)=PLON(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO(JROF,JLEV,2)=PDLO(JROF,JLEV,1)
      PDLO(JROF,JLEV,3)=PDLO(JROF,JLEV,1)
      PDLO(JROF,JLEV,0)=PDLO(JROF,JLEV,1)

      ! note: case 101 only sets linear weights, needing only 2 central lats
      ILA1=MIN(MAX(ILA+1,YDSL%NDGSAH),IDGENH)
      ILA2=MIN(MAX(ILA+2,YDSL%NDGSAH),IDGENH)
      ILO =MIN(MAX(ILO,IDLUN1),IDLUX1)
      KL0(JROF,JLEV,1)=IADDR_(IJ_,ILA1)+YDSL%NSLEXT(ILO,ILA1)
      KL0(JROF,JLEV,2)=IADDR_(IJ_,ILA2)+YDSL%NSLEXT(ILO,ILA2)
      KDEP(JROF,JLEV)= (ILA1 -ILA2 +1)*7
    ENDDO

    ! optim: does not vectorize
    DO JROF=KST,KPROF
      ILEV  =KVAUT(INT(PLEV(JROF,JLEV)*ZFAC))-1
      IF(ILEV < IFLVM2.AND.PLEV(JROF,JLEV) > PVETA(ILEV+2)) ILEV=ILEV+1

      KLEV(JROF,JLEV)=ILEV

      Z1 = PVETA(ILEV+1)
      PDVER(JROF,JLEV)=(PLEV(JROF,JLEV)-Z1)/(PVETA(ILEV+2)-Z1)  
    ENDDO

    ! * Mask calculation for on-demand communications:
    IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
!CDIR NODEP
!DIR$ PREFERVECTOR
!DIR IVDEP
      DO JROF=KST,KPROF
#ifdef __INTEL_COMPILER
        CALL MM_PREFETCH(YDSL%MASK_SL2(KL0(JROF,JLEV,1)),3)
        CALL MM_PREFETCH(YDSL%MASK_SL2(KL0(JROF,JLEV,2)),3)
#endif
        YDSL%MASK_SL2(KL0(JROF,JLEV,1):KL0(JROF,JLEV,1)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2):KL0(JROF,JLEV,2)+3)=1
      ENDDO
    ENDIF

    DO JROF=KST,KPROF
      KL0(JROF,JLEV,1)=KL0(JROF,JLEV,1)+YDSL%NASLB1*KLEV(JROF,JLEV)
      KL0(JROF,JLEV,2)=KL0(JROF,JLEV,2)+YDSL%NASLB1*KLEV(JROF,JLEV)
    ENDDO

    ! note: mirror setting for case 101 since extreme lats are not used by linear weights
    ! this is probably not useless (cf slcomm/slextpol)
    DO JROF=KST,KPROF
      KL0(JROF,JLEV,0)=KL0(JROF,JLEV,1)
      KL0(JROF,JLEV,3)=KL0(JROF,JLEV,2)
    ENDDO
  ENDDO
ENDIF

!        1.03  Coordinates and weights for ( horizontal 12 points
!              + vertical cubic + 32 points interpolations ) or
!              ( horizontal 12 points + 32 points interpolations ).
!              Optionally, Hermite cubic or cubic B-spline or WENO vertical interpolations weights are computed.

IF (KWIS == 103 .OR. KWIS == 104 .OR. KWIS == 105 .OR. KWIS == 106) THEN
  ZFAC=PVRLEVX/(PVETA(KFLEV+1)-PVETA(0))

  IF (KHOR == 1) IOFF=YDSL%NASLB1*KFLDN

  DO JLEV=1,KFLEV
    ! * Calculation of linear weights, KL0, KLH0.
!CDIR NODEP
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      IJ_ = MOD(JROF+1-KST,JPDUP)+1
      ILAG=INT(PLAT(JROF,JLEV))-1
      ILA=ILAG - YDSL%NFRSTLOFF

      ! meridional interpolation: linear weights (regular grid)
      ! general case 
      PDLAT(JROF,JLEV)=PLAT(JROF,JLEV)-REAL(ILAG+1,JPRB)

      ! zonal interpolation: linear weights for 4 lat. lines
      ! as the grid is regular in the zonal direction,
      ! the cubic weight computation does not need 
      ! other input than linear weights (LASCAW_CLO)
      ! general case
      ILO=INT(PLON(JROF,JLEV))-1
      PDLO(JROF,JLEV,0)=PLON(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO(JROF,JLEV,1)=PDLO(JROF,JLEV,0)
      PDLO(JROF,JLEV,2)=PDLO(JROF,JLEV,0)
      PDLO(JROF,JLEV,3)=PDLO(JROF,JLEV,0)

      ILA1=MIN(MAX(ILA+1,YDSL%NDGSAH),IDGENH)
      ILA2=MIN(MAX(ILA+2,YDSL%NDGSAH),IDGENH)
      ILA3=MIN(MAX(ILA+3,YDSL%NDGSAH),IDGENH)
      ILA =MIN(MAX(ILA ,YDSL%NDGSAH),IDGENH)
      ILO =MIN(MAX(ILO,IDLUN1),IDLUX1)
      KL0(JROF,JLEV,0)=IADDR_(IJ_,ILA )+YDSL%NSLEXT(ILO,ILA)
      KL0(JROF,JLEV,1)=IADDR_(IJ_,ILA1)+YDSL%NSLEXT(ILO,ILA1)
      KL0(JROF,JLEV,2)=IADDR_(IJ_,ILA2)+YDSL%NSLEXT(ILO,ILA2)
      KL0(JROF,JLEV,3)=IADDR_(IJ_,ILA3)+YDSL%NSLEXT(ILO,ILA3)
      KDEP(JROF,JLEV)=(ILA-ILA1+1) + (ILA1-ILA2+1)*2 + (ILA2-ILA3+1)*4
    ENDDO

    ! optim: does not vectorize
    DO JROF=KST,KPROF
      ! vertical interpolation: linear weights
      ! the cubic weight computation are done in 
      ! LASCAW_VINTW (including terms for grid irregularity)
      ILEV=KVAUT(INT(PLEV(JROF,JLEV)*ZFAC))-1
      IF(ILEV < IFLVM2.AND.PLEV(JROF,JLEV) > PVETA(ILEV+2)) ILEV=ILEV+1

      KLEV(JROF,JLEV)=ILEV

      Z1 = PVETA(ILEV+1)
      PDVER(JROF,JLEV)=(PLEV(JROF,JLEV)-Z1)/(PVETA(ILEV+2)-Z1)  
    ENDDO

    ! COMAD meridional and zonal interpolation 
    IF (LLCOMADH) THEN
      DO JROF=KST,KPROF
        PDLAMAD(JROF,JLEV)=PDLAT(JROF,JLEV)*PSTDDISV(JROF,JLEV)+0.5_JPRB*(1._JPRB-PSTDDISV(JROF,JLEV))
        PDLOMAD(JROF,JLEV,0)=PDLO(JROF,JLEV,0)*PSTDDISU(JROF,JLEV)+0.5_JPRB*(1._JPRB-PSTDDISU(JROF,JLEV))
        PDLOMAD(JROF,JLEV,1)=PDLO(JROF,JLEV,1)*PSTDDISU(JROF,JLEV)+0.5_JPRB*(1._JPRB-PSTDDISU(JROF,JLEV))
        PDLOMAD(JROF,JLEV,2)=PDLO(JROF,JLEV,2)*PSTDDISU(JROF,JLEV)+0.5_JPRB*(1._JPRB-PSTDDISU(JROF,JLEV))
        PDLOMAD(JROF,JLEV,3)=PDLO(JROF,JLEV,3)*PSTDDISU(JROF,JLEV)+0.5_JPRB*(1._JPRB-PSTDDISU(JROF,JLEV))
      ENDDO
    ELSE
      DO JROF=KST,KPROF
        PDLAMAD(JROF,JLEV)=PDLAT(JROF,JLEV)
        PDLOMAD(JROF,JLEV,0)=PDLO(JROF,JLEV,0)
        PDLOMAD(JROF,JLEV,1)=PDLO(JROF,JLEV,1)
        PDLOMAD(JROF,JLEV,2)=PDLO(JROF,JLEV,2)
        PDLOMAD(JROF,JLEV,3)=PDLO(JROF,JLEV,3)
      ENDDO
    ENDIF

    ! COMAD vertical interpolation 
    IF (LLCOMADV) THEN
      DO JROF=KST,KPROF
        PDVERMAD(JROF,JLEV)=PDVER(JROF,JLEV)*PSTDDISW(JROF,JLEV)+0.5_JPRB*(1._JPRB-PSTDDISW(JROF,JLEV))
      ENDDO
    ELSE
      DO JROF=KST,KPROF
        PDVERMAD(JROF,JLEV)=PDVER(JROF,JLEV)
      ENDDO
    ENDIF

    ! * Mask calculation for on-demand communications:
    IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
!CDIR NODEP
!DIR$ PREFERVECTOR
      DO JROF=KST,KPROF
#ifdef __INTEL_COMPILER
        CALL MM_PREFETCH(YDSL%MASK_SL2(KL0(JROF,JLEV,0)),3)
        CALL MM_PREFETCH(YDSL%MASK_SL2(KL0(JROF,JLEV,1)),3)
        CALL MM_PREFETCH(YDSL%MASK_SL2(KL0(JROF,JLEV,2)),3)
        CALL MM_PREFETCH(YDSL%MASK_SL2(KL0(JROF,JLEV,3)),3)
#endif
        YDSL%MASK_SL2(KL0(JROF,JLEV,0):KL0(JROF,JLEV,0)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1):KL0(JROF,JLEV,1)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2):KL0(JROF,JLEV,2)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,3):KL0(JROF,JLEV,3)+3)=1
      ENDDO
    ENDIF

    IF (KHOR == 0) IOFF=YDSL%NASLB1*JLEV

    DO JROF=KST,KPROF
      ! note: mind the different offsets JOFF/IOFF and the order
      KLH0(JROF,JLEV,0)=KL0(JROF,JLEV,0)+IOFF
      KLH0(JROF,JLEV,1)=KL0(JROF,JLEV,1)+IOFF
      KLH0(JROF,JLEV,2)=KL0(JROF,JLEV,2)+IOFF
      KLH0(JROF,JLEV,3)=KL0(JROF,JLEV,3)+IOFF

      JOFF = YDSL%NASLB1*KLEV(JROF,JLEV)
      KL0(JROF,JLEV,0)=KL0(JROF,JLEV,0)+JOFF
      KL0(JROF,JLEV,1)=KL0(JROF,JLEV,1)+JOFF
      KL0(JROF,JLEV,2)=KL0(JROF,JLEV,2)+JOFF
      KL0(JROF,JLEV,3)=KL0(JROF,JLEV,3)+JOFF
    ENDDO
  ENDDO

  IF (LLSLHD.AND.LDSLHDHEAT) THEN
    ! Computes the weights for heat fields affected by SLHD
    !  all the rest is recomputed once again bellow.

    ! * Calculation of PCLA and PCLASLD:
    CALL LASCAW_CLO(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,PDLAT,PDLAMAD,PKAPPAT,&
     & PCLA,PCLAMAD(:,:,:,1),PCLASLT)

    ! * Calculation of PCLO and PCLOSLD:
    CALL LASCAW_CLO(KFLEV,&
     & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,PDLO(:,:,1),PDLOMAD(:,:,1),PKAPPAT,&
     & PCLO(:,:,:,1,1),PCLOMAD(:,:,:,1,1),PCLOSLT(:,:,:,1))
    CALL LASCAW_CLO(KFLEV,&
     & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,PDLO(:,:,2),PDLOMAD(:,:,2),PKAPPAT,&
     & PCLO(:,:,:,2,1),PCLOMAD(:,:,:,2,1),PCLOSLT(:,:,:,2))
  ENDIF

  ! Loop over all horiz. weights for 3Dturb (computing in addition
  !                         two sets with and without SLHD)
  DO JJ=1,KDIMK
    IF ((KSPLTHOI == 1).AND.(JJ == KDIMK)) THEN
      ! Bit specific case computing diffusive weights for physical tendencies.
      ! In this case SLHD weights are of no use.

      ! * Calculation of PCLA and PCLASLD:
      CALL LASCAW_CLO(KFLEV,&
       & KPROMB,KST,KPROF,LLT_PHYS,ZSLHDKMINH,PDLAT,PDLAMAD,ZKHTURB(:,:,JJ),&
       & ZCLA,PCLAMAD(:,:,:,JJ),PCLA(:,:,:,JJ))

      ! * Calculation of PCLO and PCLOSLD:
      CALL LASCAW_CLO(KFLEV,&
       & KPROMB,KST,KPROF,LLT_PHYS,ZSLHDKMINH,PDLO(:,:,1),PDLOMAD(:,:,1),ZKHTURB(:,:,JJ),&
       & ZCLO,PCLOMAD(:,:,:,1,JJ),PCLO(:,:,:,1,JJ))
      CALL LASCAW_CLO(KFLEV,&
       & KPROMB,KST,KPROF,LLT_PHYS,ZSLHDKMINH,PDLO(:,:,2),PDLOMAD(:,:,2),ZKHTURB(:,:,JJ),&
       & ZCLO,PCLOMAD(:,:,:,2,JJ),PCLO(:,:,:,2,JJ))
    ELSE
      ! * Calculation of PCLA, PCLAMAD and PCLASLD:
      CALL LASCAW_CLO(KFLEV,&
       & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,PDLAT,PDLAMAD,PKAPPA,&
       & PCLA(:,:,:,JJ),PCLAMAD(:,:,:,JJ),PCLASLD(:,:,:,JJ))

      ! * Calculation of PCLO and PCLOSLD for central lat 1 and 2
      ! (linear int. only for lat 0 and 3)
      CALL LASCAW_CLO(KFLEV,&
       & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,PDLO(:,:,1),PDLOMAD(:,:,1),PKAPPA,&
       & PCLO(:,:,:,1,JJ),PCLOMAD(:,:,:,1,JJ),PCLOSLD(:,:,:,1,JJ))
      CALL LASCAW_CLO(KFLEV,&
       & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,PDLO(:,:,2),PDLOMAD(:,:,2),PKAPPA,&
       & PCLO(:,:,:,2,JJ),PCLOMAD(:,:,:,2,JJ),PCLOSLD(:,:,:,2,JJ))

      IF (JJ > 1.AND.LL3DTURB) THEN
        CALL LASCAW_CLOTURB(KFLEV,KPROMB,KST,KPROF,ZKHTURB(:,:,JJ),PCLA(:,:,:,JJ),&
         & PCLASLD(:,:,:,JJ))
        CALL LASCAW_CLOTURB(KFLEV,KPROMB,KST,KPROF,ZKHTURB(:,:,JJ),PCLO(:,:,:,1,JJ),&
         & PCLOSLD(:,:,:,1,JJ))
        CALL LASCAW_CLOTURB(KFLEV,KPROMB,KST,KPROF,ZKHTURB(:,:,JJ),PCLO(:,:,:,2,JJ),&
         & PCLOSLD(:,:,:,2,JJ))
      ENDIF
    ENDIF
  ENDDO

  ! note: case 102 not coded in LAM but same test as in lascaw
  IF (KWIS == 102.AND.YRDYNA%LSLTVWENO.OR.KWIS == 106) THEN
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
          IILEV(KST:KPROF,1:KFLEV)=KLEV(KST:KPROF,1:KFLEV)
        CASE (2)
          IILEV(KST:KPROF,1:KFLEV)=MIN(IFLVM2-IBCLIM,KLEV(KST:KPROF,1:KFLEV)+1)
          KNOWENO(KST:KPROF,1:KFLEV)=KNOWENO(KST:KPROF,1:KFLEV) &
           & + IILEV(KST:KPROF,1:KFLEV)-KLEV(KST:KPROF,1:KFLEV) - 1
        CASE (3)
          ! can't be  KSLEV-1 as there is no half level on -1
          IILEV(KST:KPROF,1:KFLEV)=MAX(IBCLIM,KLEV(KST:KPROF,1:KFLEV)-1)
          KNOWENO(KST:KPROF,1:KFLEV)=KNOWENO(KST:KPROF,1:KFLEV) &
           & + IILEV(KST:KPROF,1:KFLEV)-KLEV(KST:KPROF,1:KFLEV) + 1
        CASE DEFAULT
          CALL ABOR1(' LASCAW: WENO PROBLEM')
      END SELECT

      CALL LASCAW_VINTW(&
       & KPROMB,KFLEV,KST,KPROF,LLCOMADV,LLT_SLHD,LDSLHDHEAT,ZSLHDKMINV_WENO,IILEV,&
       & PLEV,PDVER,PDVERMAD,PSTDDISW,PKAPPA,PKAPPAT,PVETA,PVCUICO,PVSLD,PVSLDW,&
       & PVINTW(:,:,3*(JJ-1)+1),PVINTWMAD,PVINTWSLD,PVINTWSLT)
    ENDDO

    ! make sure it only keeps -1,0,+1 values
    IF ((MAXVAL(KNOWENO(KST:KPROF,1:KFLEV)) > 1+IBCLIM) .OR. &
     &  (MINVAL(KNOWENO(KST:KPROF,1:KFLEV)) <-1-IBCLIM))     &
     &  CALL ABOR1(' LASCAW: Something strange is happening about level shifts.')

    ! C_k functions 
    IF (LREGETA) THEN
      DO JLEV=1,KFLEV
        DO JROF=KST,KPROF
          ! regular mesh (LREGETA=.t. case)
          ! Note: This code doesn't seem to work for irregular vertical spacing.
          !       Hence the smart LWENOBC code is only allowed with LREGETA=t.
          PCW(JROF,JLEV,1)=0.10_JPRB*(2.0_JPRB+PDVER(JROF,JLEV))*(3.0_JPRB-PDVER(JROF,JLEV))   ! central
          PCW(JROF,JLEV,2)=0.05_JPRB*(2.0_JPRB+PDVER(JROF,JLEV))*(1.0_JPRB+PDVER(JROF,JLEV))   ! lower
          PCW(JROF,JLEV,3)=0.05_JPRB*(2.0_JPRB-PDVER(JROF,JLEV))*(3.0_JPRB-PDVER(JROF,JLEV))   ! upper
        ENDDO
      ENDDO
    ELSE
      DO JLEV=1,KFLEV
        ! optim: does not vectorize, intermediate variables lower gathers
        DO JROF=KST,KPROF
          ! general form
          ILEV=KLEV(JROF,JLEV)
          IF (ILEV <= 1.OR.ILEV >= KFLEV-3) CYCLE

          Z1 = PLEV(JROF,JLEV)-PVETA(ILEV-1)
          Z4 = PLEV(JROF,JLEV)-PVETA(ILEV+4)
          PCW(JROF,JLEV,1)=PGAMMA_WENO(ILEV,1)*Z1*Z4 ! central
          PCW(JROF,JLEV,2)=PGAMMA_WENO(ILEV,2)*Z1*(PLEV(JROF,JLEV)-PVETA(ILEV)) ! lower
          PCW(JROF,JLEV,3)=PGAMMA_WENO(ILEV,3)*(PLEV(JROF,JLEV)-PVETA(ILEV+3))*Z4 ! upper
        ENDDO
      ENDDO
    ENDIF

  ELSE
  ! * Calculation of PVINTW and PVINTWSLD:
    CALL LASCAW_VINTW(&
     & KPROMB,KFLEV,KST,KPROF,LLCOMADV,LLT_SLHD(1:3),LDSLHDHEAT,ZSLHDKMINV,KLEV,&
     & PLEV,PDVER,PDVERMAD,PSTDDISW,PKAPPA,PKAPPAT,PVETA,ZVCUICO_,ZVSLD_,ZVSLDW_,&
     & PVINTW,PVINTWMAD,PVINTWSLD,PVINTWSLT)
  ENDIF

  IF (KWIS == 104) THEN
    DO JLEV=1,KFLEV
      ! * Calculation of PHVW:
      DO JROF=KST,KPROF
        ZDVER=PDVER(JROF,JLEV)
        PHVW(JROF,JLEV,1)=FHLO1(ZDVER)
        PHVW(JROF,JLEV,2)=FHLO2(ZDVER)
        PHVW(JROF,JLEV,3)=FHLO3(ZDVER)
        PHVW(JROF,JLEV,4)=FHLO4(ZDVER)
      ENDDO
      ! * Calculation of PVDERW:
      ! optim: better vectorization with intermediate variables (lower gathers)
      DO JROF=KST,KPROF
        IJ_ = MOD(JROF+1-KST,JPDUP)+1
        ILEV=KLEV(JROF,JLEV)
        Z1 = ZVETA_(IJ_,ILEV+1)
        Z2 = ZVETA_(IJ_,ILEV+2)
        ZW1=(Z2-Z1)/(Z2-ZVETA_(IJ_,ILEV))
        ZW2=(Z2-Z1)/(ZVETA_(IJ_,ILEV+3)-Z1)
        IF(ILEV >= 1.AND.ILEV <= KFLEV-3) THEN
          PVDERW(JROF,JLEV,1,1)=ZW1
          PVDERW(JROF,JLEV,2,1)=ZW1
          PVDERW(JROF,JLEV,1,2)=ZW2
          PVDERW(JROF,JLEV,2,2)=ZW2
        ELSEIF (ILEV == 0) THEN
          PVDERW(JROF,JLEV,1,1)=0.0_JPRB
          PVDERW(JROF,JLEV,2,1)=2._JPRB*ZW1
          PVDERW(JROF,JLEV,1,2)=ZW2
          PVDERW(JROF,JLEV,2,2)=ZW2
        ELSEIF (ILEV == KFLEV-2) THEN
          PVDERW(JROF,JLEV,1,1)=ZW1
          PVDERW(JROF,JLEV,2,1)=ZW1
          PVDERW(JROF,JLEV,1,2)=2._JPRB*ZW2
          PVDERW(JROF,JLEV,2,2)=0.0_JPRB
        ENDIF
      ENDDO
    ENDDO
  ENDIF

  IF (KWIS == 105) THEN
    ! * Calculation of PVINTWS (weights for cubic spline interpolation).
    ! optim: does not vectorize, use of 64-bit indexing
    DO JLEV=1,KFLEV
      DO JROF=KST,KPROF
        IJ_ = MOD(JROF+1-KST,JPDUP)+1
        ILEV64=KLEV(JROF,JLEV)
        ZD2=PLEV(JROF,JLEV)-ZVETA_(IJ_,ILEV64+1)
        PVINTWS(JROF,JLEV,1)=RFVV(4,ILEV64,1)+ZD2*( RFVV(4,ILEV64,2) +&
         & ZD2*(RFVV(4,ILEV64,3) + ZD2*RFVV(4,ILEV64,4) ) )  
        PVINTWS(JROF,JLEV,2)=RFVV(3,ILEV64+1,1)+ZD2*( RFVV(3,ILEV64+1,2) +&
         & ZD2*( RFVV(3,ILEV64+1,3) + ZD2*RFVV(3,ILEV64+1,4) ) )  
        PVINTWS(JROF,JLEV,3)=RFVV(2,ILEV64+2,1)+ZD2*( RFVV(2,ILEV64+2,2) +&
         & ZD2*( RFVV(2,ILEV64+2,3) + ZD2*RFVV(2,ILEV64+2,4) ) )  
        PVINTWS(JROF,JLEV,4)=RFVV(1,ILEV64+3,1)+ZD2*( RFVV(1,ILEV64+3,2) +&
         & ZD2*( RFVV(1,ILEV64+3,3) + ZD2*RFVV(1,ILEV64+3,4) ) )  
      ENDDO
    ENDDO
  ENDIF

ENDIF

!     ----------------------------------------------------------------

!*       2.    2D MODEL AND CASES IN THE 3D MODEL WHERE ONLY
!              2D INTERPOLATIONS ARE NEEDED.
!              ---------------------------------------------

!        2.01  Coordinates and weights for bilinear interpolations.

IF (KWIS == 201) THEN
  IF (KHOR == 1) IOFF=YDSL%NASLB1*KFLDN

  DO JLEV=1,KFLEV
    IF (KHOR == 0) IOFF=YDSL%NASLB1*JLEV

    ! * Calculation of linear weights, KL0.
!CDIR NODEP
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      IJ_ = MOD(JROF+1-KST,JPDUP)+1

      ILAG=INT(PLAT(JROF,JLEV))-1
      ILA=ILAG - YDSL%NFRSTLOFF
      PDLAT(JROF,JLEV)=PLAT(JROF,JLEV)-REAL(ILAG+1,JPRB)
      ILO=INT(PLON(JROF,JLEV))-1
      PDLO(JROF,JLEV,1)=PLON(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO(JROF,JLEV,2)=PDLO(JROF,JLEV,1)

      ILA1=MIN(MAX(ILA+1,YDSL%NDGSAH),IDGENH)
      ILA2=MIN(MAX(ILA+2,YDSL%NDGSAH),IDGENH)
      ILO =MIN(MAX(ILO,IDLUN1),IDLUX1)
      KL0(JROF,JLEV,1)=IADDR_(IJ_,ILA1)+YDSL%NSLEXT(ILO,ILA1)+IOFF
      KL0(JROF,JLEV,2)=IADDR_(IJ_,ILA2)+YDSL%NSLEXT(ILO,ILA2)+IOFF
      KDEP(JROF,JLEV)= (ILA1 -ILA2 +1)*7
    ENDDO

    ! * Mask calculation for on-demand communications:
    IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
!CDIR NODEP
!DIR$ PREFERVECTOR
      DO JROF=KST,KPROF
#ifdef __INTEL_COMPILER
        CALL MM_PREFETCH(YDSL%MASK_SL2(KL0(JROF,JLEV,1)),3)
        CALL MM_PREFETCH(YDSL%MASK_SL2(KL0(JROF,JLEV,2)),3)
#endif
        YDSL%MASK_SL2(KL0(JROF,JLEV,1):KL0(JROF,JLEV,1)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2):KL0(JROF,JLEV,2)+3)=1
      ENDDO
    ENDIF
  ENDDO
ENDIF

!        2.03  Coordinates and weights for 12 points interpolations.

IF (KWIS == 203) THEN
  IF (KHOR == 1) IOFF=YDSL%NASLB1*KFLDN

  DO JLEV=1,KFLEV
    IF (KHOR == 0) IOFF=YDSL%NASLB1*JLEV

    ! * Calculation of linear weights, KL0.
!CDIR NODEP
!DIR$ PREFERVECTOR
    DO JROF=KST,KPROF
      IJ_ = MOD(JROF+1-KST,JPDUP)+1

      ILAG=INT(PLAT(JROF,JLEV))-1
      ILA=ILAG - YDSL%NFRSTLOFF

      ! meridional interpolation: linear weights (regular grid)
      ! general case 
      PDLAT(JROF,JLEV)=PLAT(JROF,JLEV)-REAL(ILAG+1,JPRB)

      ! zonal interpolation: linear weights for 4 lat. lines
      ! as the grid is regular in the zonal direction,
      ! the cubic weight computation does not need 
      ! other input than linear weights (LASCAW_CLO)
      ! general case
      ILO =INT(PLON(JROF,JLEV))-1
      PDLO(JROF,JLEV,0)=PLON(JROF,JLEV)-REAL(ILO+1,JPRB)
      PDLO(JROF,JLEV,1)=PDLO(JROF,JLEV,0)
      PDLO(JROF,JLEV,2)=PDLO(JROF,JLEV,0)
      PDLO(JROF,JLEV,3)=PDLO(JROF,JLEV,0)

      ILA1=MIN(MAX(ILA+1,YDSL%NDGSAH),IDGENH)
      ILA2=MIN(MAX(ILA+2,YDSL%NDGSAH),IDGENH)
      ILA3=MIN(MAX(ILA+3,YDSL%NDGSAH),IDGENH)
      ILA =MIN(MAX(ILA ,YDSL%NDGSAH),IDGENH)
      ILO =MIN(MAX(ILO,IDLUN1),IDLUX1)
      KL0(JROF,JLEV,0)=IADDR_(IJ_,ILA )+YDSL%NSLEXT(ILO,ILA)+IOFF
      KL0(JROF,JLEV,1)=IADDR_(IJ_,ILA1)+YDSL%NSLEXT(ILO,ILA1)+IOFF
      KL0(JROF,JLEV,2)=IADDR_(IJ_,ILA2)+YDSL%NSLEXT(ILO,ILA2)+IOFF
      KL0(JROF,JLEV,3)=IADDR_(IJ_,ILA3)+YDSL%NSLEXT(ILO,ILA3)+IOFF
      KDEP(JROF,JLEV)=(ILA  -ILA1 +1) + (ILA1 -ILA2 +1)*2 + (ILA2 -ILA3 +1)*4
    ENDDO

    ! * Mask calculation for on-demand communications:
    IF(NPROC > 1.AND.YDSL%LSLONDEM_ACTIVE)THEN
!CDIR NODEP
!DIR$ PREFERVECTOR
      DO JROF=KST,KPROF
#ifdef __INTEL_COMPILER
        CALL MM_PREFETCH(YDSL%MASK_SL2(KL0(JROF,JLEV,0)),3)
        CALL MM_PREFETCH(YDSL%MASK_SL2(KL0(JROF,JLEV,1)),3)
        CALL MM_PREFETCH(YDSL%MASK_SL2(KL0(JROF,JLEV,2)),3)
        CALL MM_PREFETCH(YDSL%MASK_SL2(KL0(JROF,JLEV,3)),3)
#endif
        YDSL%MASK_SL2(KL0(JROF,JLEV,0):KL0(JROF,JLEV,0)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,1):KL0(JROF,JLEV,1)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,2):KL0(JROF,JLEV,2)+3)=1
        YDSL%MASK_SL2(KL0(JROF,JLEV,3):KL0(JROF,JLEV,3)+3)=1
      ENDDO
    ENDIF

    ! COMAD meridional and zonal interpolation 
    IF (LLCOMADH) THEN
      DO JROF=KST,KPROF
        PDLAMAD(JROF,JLEV)=0.5_JPRB+(PDLAT(JROF,JLEV)-0.5_JPRB)*PSTDDISV(JROF,JLEV)
        PDLOMAD(JROF,JLEV,0)=0.5_JPRB+(PDLO(JROF,JLEV,0)-0.5_JPRB)*PSTDDISU(JROF,JLEV)
        PDLOMAD(JROF,JLEV,1)=0.5_JPRB+(PDLO(JROF,JLEV,1)-0.5_JPRB)*PSTDDISU(JROF,JLEV)
        PDLOMAD(JROF,JLEV,2)=0.5_JPRB+(PDLO(JROF,JLEV,2)-0.5_JPRB)*PSTDDISU(JROF,JLEV)
        PDLOMAD(JROF,JLEV,3)=0.5_JPRB+(PDLO(JROF,JLEV,3)-0.5_JPRB)*PSTDDISU(JROF,JLEV)
      ENDDO
    ELSE
      DO JROF=KST,KPROF
        PDLAMAD(JROF,JLEV)=PDLAT(JROF,JLEV)
        PDLOMAD(JROF,JLEV,0)=PDLO(JROF,JLEV,0)
        PDLOMAD(JROF,JLEV,1)=PDLO(JROF,JLEV,1)
        PDLOMAD(JROF,JLEV,2)=PDLO(JROF,JLEV,2)
        PDLOMAD(JROF,JLEV,3)=PDLO(JROF,JLEV,3)
      ENDDO
    ENDIF
  ENDDO

  ! * Calculation of PCLA and PCLASLD:
  CALL LASCAW_CLO(KFLEV,KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,PDLAT,PDLAMAD,PKAPPA,&
   & PCLA(:,:,:,1),PCLAMAD(:,:,:,1),PCLASLD(:,:,:,1))

  ! * Calculation of PCLO and PCLOSLD:
  CALL LASCAW_CLO(KFLEV,&
   & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,PDLO(:,:,1),PDLOMAD(:,:,1),PKAPPA,&
   & PCLO(:,:,:,1,1),PCLOMAD(:,:,:,1,1),PCLOSLD(:,:,:,1,1))
  CALL LASCAW_CLO(KFLEV,&
   & KPROMB,KST,KPROF,LLT_SLHD,ZSLHDKMINH,PDLO(:,:,2),PDLOMAD(:,:,2),PKAPPA,&
   & PCLO(:,:,:,2,1),PCLOMAD(:,:,:,2,1),PCLOSLD(:,:,:,2,1))
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ELASCAW',1,ZHOOK_HANDLE)
END SUBROUTINE ELASCAW
