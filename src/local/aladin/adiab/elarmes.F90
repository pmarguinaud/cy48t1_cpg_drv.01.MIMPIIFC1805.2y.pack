#ifdef RS6K
@PROCESS NOCHECK
#endif
!option! -O nomove
SUBROUTINE ELARMES(YDGEOMETRY,YDML_DYN,YDRIP,KST,KPROF,YDSL,KSTABUF,PGMDTX,PGMDTY,PB1,PB2,&
 & PLSDEPI,KIBL,&
 & KVSEPC,KVSEPL,&
 & PSCO,PLEV,PCCO,PUF,PVF,&
 & KL0,KLH0,KLEV,PLSCAW,PRSCAW)  

!     ------------------------------------------------------------------
!**** *ELARMES - semi-LAgrangian scheme:
!                Research of the origin point on the Sphere.

!     Purpose.
!     --------

!      The computation of the location of the interpolation point of
!     the lagrangian trajectory is performed by an iterative
!     method described by Robert and adapted to the sphere by M. Rochas.
!     Trajectories are great circles for spherical geometry.
!      In ELARMES the trajectories are straight lines on the plane.
!     Finally we find the departure (origin) point "O".

!**   Interface.
!     ----------
!        *CALL* *ELARMES(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST      - first element of arrays where computations are performed.
!          KPROF    - depth of work.
!          YDSL     - SL_STRUCT definition.
!          KSTABUF  - for a latitude IGL, KSTABUF(IGL) is the
!                     address of the element corresponding to
!                     (ILON=1,IGL) in the NPROMA arrays.
!          PGMDTX   - m * DELTA t / DELTA x.
!          PGMDTY   - m * DELTA t / DELTA y.
!          PB1      - SLBUF1-buffer for interpolations.
!          PB2      - SLBUF2-buf to communicate info from non lag. to lag. dyn.
!          PLSDEPI  - (Number of points by latitude) / (2 * PI) .
!          KIBL     - index into YRCSGEOM/YRGSGEOM instances in YDGEOMETRY

!        INPUT/OUTPUT:
!          KVSEPC   - vertical separation (used in S/L adjoint, cubic interp.)
!          KVSEPL   - vertical separation (used in S/L adjoint, linear interp.)

!        OUTPUT:
!          PSCO     - information about geographic position of interpol. point.
!          PLEV     - vertical coordinate of the interpolation point.
!          PCCO     - information about comput. space position of interpol. point.
!          PUF      - U-comp of wind necessary to
!                     find the position of the origin point,
!                     in a local repere linked to computational sphere.
!          PVF      - V-comp of wind necessary to
!                     find the position of the origin point,
!                     in a local repere linked to computational sphere.
!          KL0      - index of the four western points
!                     of the 16 points interpolation grid.
!          KLH0     - second value of index of the four western points
!                     of the 16 points interpolation grid if needed.
!          KLEV     - lower level of the vertical interpolation
!                     grid needed for vertical interpolations.
!          PLSCAW   - linear weights (distances) for interpolations.
!          PRSCAW   - non-linear weights for interpolations.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------
!        Calls  LARCINA.
!        Is called by LAPINEA (3D model)

!     Reference.
!     ----------

!     Author.
!     -------
!        R. Bubnova after LARMES (by Karim Yessad) and developments
!        by Martin Janousek      CNRM/GMAP/EXT
!        Original : JANUARY 1995

!     Modifications.
!     --------------
!        K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!        K. Yessad (Aug 2009): always use root (QX,QY) for (p,q) variables names
!        K. Yessad (Nov 2009): keep LLO.OR.LELTRA=T only in SL2TL.
!        K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!        F. Vana      : 23-Feb-2011  LARCINA arguments update
!        G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!        G. Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
!        K. Yessad (Nov 2011): introduce LRALTVDISP; do same printings as in LARMES.
!        S. Malardel and D. Ricard (Nov 2013): COMAD weights for SL interpolations
!        B. Bochenek (Apr 2015): Phasing: move some variables.
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!        K. Yessad (Feb 2018): remove deep-layer formulations.
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
!        R. El Khatib 22-May-2019 fix the "IFL OUT OF BOUNDS" issue
! End Modifications
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK
USE SC2PRG_MOD, ONLY : SC2PRG

USE YOMCST    , ONLY : RPI
USE YOMCT0    , ONLY : NCONF, LTWOTL
USE YOMDYNA   , ONLY : YRDYNA
USE YOMLUN    , ONLY : NULERR
USE YOMRIP    , ONLY : TRIP
USE EINT_MOD  , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN)       :: YDGEOMETRY
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
TYPE(TRIP)     ,INTENT(IN)       :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
TYPE(SL_STRUCT),   INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTABUF(YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMDTX(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMDTY(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB1(YDSL%NASLB1,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
INTEGER(KIND=JPIM),INTENT(INOUT) :: KVSEPC
INTEGER(KIND=JPIM),INTENT(INOUT) :: KVSEPL
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTSCO%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLH0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IHVI
INTEGER(KIND=JPIM) :: ILEVEXP
INTEGER(KIND=JPIM) :: IROFEXP
INTEGER(KIND=JPIM) :: IROT
INTEGER(KIND=JPIM) :: ISTESB(YDGEOMETRY%YRDIM%NPROMA)
INTEGER(KIND=JPIM) :: ISTEST(YDGEOMETRY%YRDIM%NPROMA)
INTEGER(KIND=JPIM) :: ISTEST0(YDGEOMETRY%YRDIM%NPROMA)
INTEGER(KIND=JPIM) :: ITIP
INTEGER(KIND=JPIM) :: JITER
INTEGER(KIND=JPIM) :: JLEV
INTEGER(KIND=JPIM) :: JROF
INTEGER(KIND=JPIM) :: IL0A (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
INTEGER(KIND=JPIM) :: IDEP (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: INOWENO (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

LOGICAL :: LLINTV
LOGICAL :: LLO
LOGICAL :: LLSLHD, LLSLHDQUAD

REAL(KIND=JPRB) :: ZDGUN(YDGEOMETRY%YRDIM%NPROMA), ZDGUX(YDGEOMETRY%YRDIM%NPROMA), ZDLUN, ZDLUX
REAL(KIND=JPRB) :: ZEW, ZEWX
REAL(KIND=JPRB) :: ZINDX, ZINDY, ZINEZ, ZINEZV(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZLEVB, ZLEVO, ZLEVT
REAL(KIND=JPRB) :: ZNOR, ZNORX
REAL(KIND=JPRB) :: ZPU, ZPV
REAL(KIND=JPRB) :: ZTXO, ZTYO
REAL(KIND=JPRB) :: ZVMAX1, ZVMAX2
REAL(KIND=JPRB) :: ZEPS
REAL(KIND=JPRB) :: ZVETAON, ZVETAOX

REAL(KIND=JPRB) :: ZZF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZWF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZWFSM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZWFASM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTSCO%NDIM)
REAL(KIND=JPRB) :: ZPLEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB), POINTER :: ZSLB1UR0(:), ZSLB1VR0(:), ZSLB1WR0(:), ZSLB1WRA(:)
REAL(KIND=JPRB), POINTER :: ZSLB1ZR0(:), ZSLB2STDDISU(:), ZSLB2STDDISV(:), ZSLB2STDDISW(:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! unused arguments in call to LARCINA:
! a) input arrays - will not be used since LSLHD was set to .FALSE.
!    in LAPINEA before call to LARMES/ELARMES => dimensions can be
!    contracted
REAL(KIND=JPRB) :: ZUN_KAPPA(1),ZUN_KAPPAT(1),ZUN_KAPPAM(1),ZUN_KAPPAH(1)

REAL(KIND=JPRB) :: ZVDISP_1,ZZ,ZVDISP_2,ZVDISP
REAL(KIND=JPRB) :: ZWF_2(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZWO_2(YDGEOMETRY%YRDIM%NPROMA, YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZHVETAON,ZHVETAOX

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "larcina.intfb.h"
#include "laismoa.intfb.h"

!     ------------------------------------------------------------------


IF (LHOOK) CALL DR_HOOK('ELARMES',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDMP=>YDGEOMETRY%YRMP, YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL), YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL) , &
 & YDPTRSLB2=>YDML_DYN%YRPTRSLB2, &
 & YDPTRSLB1=>YDML_DYN%YRPTRSLB1,YDDYN=>YDML_DYN%YRDYN, &
 & YDVETA=>YDGEOMETRY%YRVETA)

ASSOCIATE( &
 & NFLEVG=>YDDIMV%NFLEVG, NFLEN=>YDDIMV%NFLEN, NFLSA=>YDDIMV%NFLSA, &
 & NPROMA=>YDDIM%NPROMA, &
 & NITMP=>YDDYN%NITMP, VMAX1=>YDDYN%VMAX1, VMAX2=>YDDYN%VMAX2, &
 & VETAON=>YDDYN%VETAON, VETAOX=>YDDYN%VETAOX, &
 & LSVTSM=>YDDYN%LSVTSM, LFINDVSEP=>YDDYN%LFINDVSEP, &
 & RPRES_SVTSM=>YDDYN%RPRES_SVTSM, &
 & MSLB1UR0=>YDPTRSLB1%MSLB1UR0, MSLB1VR0=>YDPTRSLB1%MSLB1VR0, &
 & MSLB1ZR0=>YDPTRSLB1%MSLB1ZR0, &
 & MSLB1WR0=>YDPTRSLB1%MSLB1WR0, MSLB1WRA=>YDPTRSLB1%MSLB1WRA, &
 & MSLB2URL=>YDPTRSLB2%MSLB2URL, MSLB2VRL=>YDPTRSLB2%MSLB2VRL, &
 & MSLB2WRL=>YDPTRSLB2%MSLB2WRL, MSLB2STDDISU=>YDPTRSLB2%MSLB2STDDISU, &
 & MSLB2STDDISV=>YDPTRSLB2%MSLB2STDDISV, &
 & MSLB2STDDISW=>YDPTRSLB2%MSLB2STDDISW, RTDT=>YDRIP%RTDT, &
 & STPREH=>YDGEOMETRY%YRSTA%STPREH)

CALL SC2PRG(MSLB1UR0  ,PB1      ,ZSLB1UR0)
CALL SC2PRG(MSLB1VR0  ,PB1      ,ZSLB1VR0)
CALL SC2PRG(MSLB1WR0  ,PB1      ,ZSLB1WR0)
CALL SC2PRG(MSLB1WRA  ,PB1      ,ZSLB1WRA)
CALL SC2PRG(MSLB1ZR0  ,PB1      ,ZSLB1ZR0)
CALL SC2PRG(MSLB2STDDISU,PB2    ,ZSLB2STDDISU)
CALL SC2PRG(MSLB2STDDISV,PB2    ,ZSLB2STDDISV)
CALL SC2PRG(MSLB2STDDISW,PB2    ,ZSLB2STDDISW)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS AND TESTS.
!              --------------------------------------

!*       1.0   Set bounds of latitudes in such a way that 

! Set bounds of Y rows in such a way that :
! - gridpoints in C+I will remain in C+I wherever the origin point of the trajectory is
! - gridpoint in upper part of E zone will stay along their Y row.
! This will force the upper part of the E-zone to remain within its SL halo
! independently of any C+I consideration.
DO JROF=KST,KPROF
  IF (YDGEOMETRY%YRGSGEOM_NB%NGPLAT(JROF) > YDGEOMETRY%YRDIM%NDGUXG) THEN
    ZDGUN(JROF)=REAL((YDGEOMETRY%YRGSGEOM_NB%NGPLAT(JROF)-YDSL%NFRSTLOFF),JPRB)
    ZDGUX(JROF)=ZDGUN(JROF)
  ELSE
    ZDGUN(JROF)=REAL(YDSL%NDGUNG,JPRB)
    ZDGUX(JROF)=REAL(YDSL%NDGUXG,JPRB)
  ENDIF
ENDDO

!*       1.1   Test that wind is not too strong.

ZNOR=0.0_JPRB
ZEW=0.0_JPRB
ZDLUN=REAL(YDSL%NDLUNG,JPRB)
ZDLUX=REAL(YDSL%NDLUXG,JPRB)

ZVMAX1=VMAX1*VMAX1
ZVMAX2=VMAX2*VMAX2
ZEPS=1.E-6_JPRB
DO JROF=KST,KPROF
  ZINEZV(JROF)=MAX(0.0_JPRB,SIGN(1.0_JPRB,(YDCSGEOM%RINDX(JROF)-ZDLUN)*&
   & (ZDLUX-YDCSGEOM%RINDX(JROF))))*&
   & MAX(0.0_JPRB,SIGN(1.0_JPRB,(YDCSGEOM%RINDY(JROF)-ZDGUN(JROF))*(ZDGUX(JROF)-YDCSGEOM%RINDY(JROF))))  
ENDDO

DO JLEV=1,NFLEVG
  DO JROF=KST,KPROF
    ZNOR=MAX(PB2(JROF,MSLB2VRL+JLEV-1)*PB2(JROF,MSLB2VRL+JLEV-1),ZNOR)
    ZEW =MAX(PB2(JROF,MSLB2URL+JLEV-1)*PB2(JROF,MSLB2URL+JLEV-1),ZEW )
    PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RQX)=1.0_JPRB
    PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RQY)=1.0_JPRB
  ENDDO
ENDDO

IF (ZNOR > ZVMAX1) THEN
  WRITE(NULERR,*) ' MAX V WIND=',SQRT(ZNOR)
ENDIF
IF (ZEW > ZVMAX1) THEN
  WRITE(NULERR,*) ' MAX U WIND=',SQRT(ZEW)
ENDIF
IF (ZNOR > ZVMAX2) THEN
  ZNOR=0.0_JPRB
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      ZNORX=PB2(JROF,MSLB2VRL+JLEV-1)*PB2(JROF,MSLB2VRL+JLEV-1)
      ZNOR=MAX(ZNORX,ZNOR)
      IF(ZNOR == ZNORX) THEN
        ILEVEXP=JLEV
        IROFEXP=JROF
      ENDIF
    ENDDO
  ENDDO
  WRITE(NULERR,*) ' V WIND =',SQRT(ZNOR),' IS TOO STRONG, EXPLOSION.'
  WRITE(NULERR,*) ' LEVEL= ',ILEVEXP,' POINT= ',IROFEXP
  WRITE(NULERR,*) ' LON  = ',ACOS(YDCSGEOM%RCOLON(IROFEXP))*180._JPRB/RPI,' degrees'
  WRITE(NULERR,*) ' LAT  = ',ASIN(YDGSGEOM%GEMU(IROFEXP))*180._JPRB/RPI,' degrees'
  CALL ABOR1(' !V WIND TOO STRONG, EXPLOSION!!!')
ENDIF
IF (ZEW > ZVMAX2) THEN
  ZEW=0.0_JPRB
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      ZEWX=PB2(JROF,MSLB2URL+JLEV-1)*PB2(JROF,MSLB2URL+JLEV-1)
      ZEW=MAX(ZEWX,ZEW)
      IF(ZEW == ZEWX) THEN
        ILEVEXP=JLEV
        IROFEXP=JROF
      ENDIF
    ENDDO
  ENDDO
  WRITE(NULERR,*) ' U WIND =',SQRT(ZEW),' IS TOO STRONG, EXPLOSION.'
  WRITE(NULERR,*) ' LEVEL= ',ILEVEXP,' POINT= ',IROFEXP
  WRITE(NULERR,*) ' LON  = ',ACOS(YDCSGEOM%RCOLON(IROFEXP))*180._JPRB/RPI,' degrees'
  WRITE(NULERR,*) ' LAT  = ',ASIN(YDGSGEOM%GEMU(IROFEXP))*180._JPRB/RPI,' degrees'
  CALL ABOR1(' !U WIND TOO STRONG, EXPLOSION!!!')
ENDIF

!*       1.2   Miscellaneous preliminary initialisations.

! in practical LLO.OR.LELTRA should now be always T for SL2TL.
IF (LTWOTL) THEN
  LLO=.NOT.YRDYNA%LELTRA
ELSE
  LLO=.FALSE.
ENDIF

! deactivate computation of SLHD weights
LLSLHD     =.FALSE.
LLSLHDQUAD =.FALSE.

!*       1.3   Computation of weights for smooth interpolation at arrival point.

IF(LSVTSM) THEN
  !*     SMOOTH INTERPOLATION AT ARRIVAL POINT
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      ZSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_SINLA)=YDGSGEOM%GEMU(JROF)-ZEPS
      ZSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO)=YDGSGEOM%GSQM2(JROF)
      ZSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO)=ZEPS
      ZSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI)=1.0_JPRB-ZEPS
      ZPLEV(JROF,JLEV)=YDVETA%VETAF(JLEV)
    ENDDO
  ENDDO
  IHVI=0
  IROT=0
  ITIP=1
  LLINTV=.FALSE.
  CALL LARCINA(YDGEOMETRY,YDML_DYN,KST,KPROF,YDSL,IHVI,KSTABUF,LFINDVSEP,LLSLHD,LLSLHDQUAD,LLINTV,&
   & ITIP,IROT,.TRUE.,PLSDEPI,&
   & KIBL,ZSCO,ZPLEV,&
   & ZUN_KAPPA,ZUN_KAPPAT,ZUN_KAPPAM,ZUN_KAPPAH,&
   & ZSLB2STDDISU,ZSLB2STDDISV,ZSLB2STDDISW,&
   & ZSLB1UR0,ZSLB1VR0,ZSLB1ZR0,ZSLB1WRA,&
   & KVSEPC,KVSEPL,PCCO,&
   & PUF,PVF,ZZF,ZWF,ZWFASM,&
   & IL0A,KLH0,KLEV,PLSCAW,PRSCAW,IDEP,INOWENO)
  CALL LAISMOA(YDSL%NASLB1,NPROMA,KST,KPROF,NFLEVG,NFLSA,&
   & NFLEN,PLSCAW(1,1,YDML_DYN%YYTLSCAW%M_WDLO),IL0A,PLSCAW(1,1,YDML_DYN%YYTLSCAW%M_WDVER),ZSLB1WRA,ZWFASM)  
ENDIF

!     ------------------------------------------------------------------

!*       2.    ITERATIONS.
!              -----------

ISTEST0(:)=0
DO JITER=1,NITMP

  !*       2.1   DETERMINATION OF THE MEDIUM POINT "M" OR THE ORIGIN POINT "O".

  ! Computation of the coordinates of the medium or origin point.
  ! If (LLO=T or LELTRA=T) the origin point "O" is computed
  ! instead of the medium point "M".

  ZLEVT=YDVETA%VETAF(0)
  ZLEVB=YDVETA%VETAF(NFLEVG+1)
  ZVETAON=(1.0_JPRB-VETAON)*YDVETA%VETAH(0)+VETAON*YDVETA%VETAF(1)
  ZVETAOX=(1.0_JPRB-VETAOX)*YDVETA%VETAH(NFLEVG)+VETAOX*YDVETA%VETAF(NFLEVG)
  IF (YRDYNA%LRALTVDISP) THEN
    ZHVETAON=0.5_JPRB*ZVETAON
    ZHVETAOX=1.0_JPRB-0.5_JPRB*(1.0_JPRB-ZVETAOX)
  ENDIF
  ISTEST(:)=0
  ISTESB(:)=0

  IF (JITER == 1) THEN

    DO JLEV=1,NFLEVG
      DO JROF=KST,KPROF

        ZINEZ=ZINEZV(JROF)

        ! * computations on horizontal plans.

        ZINDX=YDCSGEOM%RINDX(JROF)
        ZINDY=YDCSGEOM%RINDY(JROF)

        !   - Compute the relative coordinates of departure point of trajectory

        PUF(JROF,JLEV)=PB2(JROF,MSLB2URL+JLEV-1)
        PVF(JROF,JLEV)=PB2(JROF,MSLB2VRL+JLEV-1)
        ZTXO = ZINDX-2.0_JPRB*PB2(JROF,MSLB2URL+JLEV-1)*PGMDTX(JROF)
        ZTYO = ZINDY-2.0_JPRB*PB2(JROF,MSLB2VRL+JLEV-1)*PGMDTY(JROF)
        PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLON) = ZTXO*ZINEZ +ZINDX*(1.0_JPRB-ZINEZ)
        PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLAT) = ZTYO*ZINEZ +ZINDY*(1.0_JPRB-ZINEZ)

        !   - Return back the departure point if it is out of the C+I zone

        ZTXO = MIN(MAX(ZTXO,ZDLUN),ZDLUX)
        ZTYO = MIN(MAX(ZTYO,ZDGUN(JROF)),ZDGUX(JROF))

        !   - Fill the array elements by coordinates of the origin point
        !     Set the origin point to the arrival point if it is left or right
        !     out of C+I zone

        IF(LLO.OR.YRDYNA%LELTRA) THEN
          PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO) = ZTXO*ZINEZ+ZINDX*(1.0_JPRB-ZINEZ)
          PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO) = ZTYO*ZINEZ+ZINDY*(1.0_JPRB-ZINEZ)
        ELSE
          PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO) = 0.5_JPRB*(ZTXO+ZINDX)*ZINEZ&
           & + ZINDX*(1.0_JPRB-ZINEZ)
          PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO) = 0.5_JPRB*(ZTYO+ZINDY)*ZINEZ&
           & + ZINDY*(1.0_JPRB-ZINEZ)
        ENDIF

        ! * computations on a vertical.

        IF(STPREH(JLEV) > RPRES_SVTSM .OR. (NCONF /= 1 .AND. NCONF /= 302) .OR..NOT.LSVTSM) THEN
          ZWF(JROF,JLEV)=PB2(JROF,MSLB2WRL+JLEV-1)
        ELSE
          ZWF(JROF,JLEV)=ZWFASM(JROF,JLEV)
        ENDIF
        IF (YRDYNA%LRALTVDISP) THEN
          ZVDISP_1=-RTDT*ZINEZ*ZWF(JROF,JLEV)
          ZZ=EXP(-RTDT*ZINEZ*ZWF(JROF,JLEV)*(1.0_JPRB/(YDVETA%VETAF(JLEV)-ZHVETAON)&
           & +1.0_JPRB/(ZHVETAOX-YDVETA%VETAF(JLEV))))&
           & *(YDVETA%VETAF(JLEV)-ZHVETAON)/(ZHVETAOX-YDVETA%VETAF(JLEV))
          ZVDISP_2=(ZHVETAON+ZHVETAOX*ZZ)/(1.0_JPRB+ZZ)-YDVETA%VETAF(JLEV)
          ZVDISP=SIGN(1.0_JPRB,ZVDISP_1)*MIN(ABS(ZVDISP_1),ABS(ZVDISP_2))
          ZLEVO=YDVETA%VETAF(JLEV)+ZVDISP
        ELSE
          ZLEVO=YDVETA%VETAF(JLEV)-RTDT*ZINEZ*ZWF(JROF,JLEV)
        ENDIF
        ISTEST(JROF)=ISTEST(JROF)-MIN(0,MAX(-1,NINT(ZLEVO-ZLEVT-0.5_JPRB)))
        ISTESB(JROF)=ISTESB(JROF)-MIN(0,MAX(-1,NINT(ZLEVB-ZLEVO-0.5_JPRB)))
        ZLEVO=MIN(ZVETAOX,MAX(ZVETAON,ZLEVO))
        IF(LLO.OR.YRDYNA%LELTRA) THEN
          PLEV(JROF,JLEV)=ZLEVO
        ELSE
          PLEV(JROF,JLEV)=0.5_JPRB*(ZLEVO+YDVETA%VETAF(JLEV))
        ENDIF

      ENDDO
    ENDDO

  ELSE

    DO JLEV=1,NFLEVG
      DO JROF=KST,KPROF

        ZINEZ=ZINEZV(JROF)

        ! * computations on horizontal plans.

        ZINDX=YDCSGEOM%RINDX(JROF)
        ZINDY=YDCSGEOM%RINDY(JROF)

        !   - Compute the relative coordinates of departure point of trajectory
        !     ZPU,ZPV are the coordinates of VM in the local repere related to F

        IF(YRDYNA%LELTRA) THEN
          ZPU=PUF(JROF,JLEV)
          ZPV=PVF(JROF,JLEV)
        ELSEIF(LLO) THEN
          ZPU=0.5_JPRB*(PUF(JROF,JLEV)+PB2(JROF,MSLB2URL+JLEV-1))
          ZPV=0.5_JPRB*(PVF(JROF,JLEV)+PB2(JROF,MSLB2VRL+JLEV-1))
        ELSE
          ZPU=PUF(JROF,JLEV)
          ZPV=PVF(JROF,JLEV)
        ENDIF

        ZTXO = ZINDX-2.0_JPRB*ZPU*PGMDTX(JROF)
        ZTYO = ZINDY-2.0_JPRB*ZPV*PGMDTY(JROF)

        PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLON) = ZTXO*ZINEZ +ZINDX*(1.0_JPRB-ZINEZ)
        PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLAT) = ZTYO*ZINEZ +ZINDY*(1.0_JPRB-ZINEZ)

        !   - Return back the departure point if it is out of the C+I zone

        ZTXO = MIN(MAX(ZTXO,ZDLUN),ZDLUX)
        ZTYO = MIN(MAX(ZTYO,ZDGUN(JROF)),ZDGUX(JROF))

        !   - Fill the array elements by coordinates of the origin point
        !     Set the origin point to the arrival point if it is left or right
        !     out of C+I zone

        IF(LLO.OR.YRDYNA%LELTRA) THEN
          PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO) = ZTXO*ZINEZ+ZINDX*(1.0_JPRB-ZINEZ)
          PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO) = ZTYO*ZINEZ+ZINDY*(1.0_JPRB-ZINEZ)
        ELSE
          PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO) = 0.5_JPRB*(ZTXO+ZINDX)*ZINEZ&
           & + ZINDX*(1.0_JPRB-ZINEZ)
          PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO) = 0.5_JPRB*(ZTYO+ZINDY)*ZINEZ&
           & + ZINDY*(1.0_JPRB-ZINEZ)
        ENDIF
        IF(JITER == NITMP) THEN
          ! * save (zpu,zpv) in (puf,pvf) with a rotation equal to identity
          !   in plane geometry.
          PUF(JROF,JLEV)=ZPU
          PVF(JROF,JLEV)=ZPV
        ENDIF

        ! * computations on a vertical.

        IF(LLO) THEN
          IF(STPREH(JLEV) > RPRES_SVTSM .OR. (NCONF /= 1 .AND. NCONF /= 302) .OR..NOT.LSVTSM) THEN
            IF (YRDYNA%LRALTVDISP) THEN
              ZWF_2(JROF,JLEV)=ZWF(JROF,JLEV)
              ZWO_2(JROF,JLEV)=PB2(JROF,MSLB2WRL+JLEV-1)
            ENDIF
            ZWF(JROF,JLEV)=0.5_JPRB*(ZWF(JROF,JLEV)+PB2(JROF,MSLB2WRL+JLEV-1)*ZINEZ)
          ELSE
            IF (YRDYNA%LRALTVDISP) THEN
              ZWF_2(JROF,JLEV)=ZWFASM(JROF,JLEV)
              ZWO_2(JROF,JLEV)=ZWFSM(JROF,JLEV)
            ENDIF
            ZWF(JROF,JLEV)=0.5_JPRB*(ZWFSM(JROF,JLEV)+ZWFASM(JROF,JLEV)*ZINEZ)
          ENDIF
        ENDIF

        IF (YRDYNA%LRALTVDISP) THEN
          ZVDISP_1=-RTDT*ZINEZ*ZWF(JROF,JLEV)
          IF (LLO) THEN
            ZZ=EXP(-0.5_JPRB*RTDT*ZINEZ*ZWF_2(JROF,JLEV)*(1.0_JPRB/(YDVETA%VETAF(JLEV)-ZHVETAON)&
             & +1.0_JPRB/(ZHVETAOX-YDVETA%VETAF(JLEV)))&
             & -0.5_JPRB*RTDT*ZINEZ*ZWO_2(JROF,JLEV)*(1.0_JPRB/(PLEV(JROF,JLEV)-ZHVETAON)&
             & +1.0_JPRB/(ZHVETAOX-PLEV(JROF,JLEV))))&
             & *(YDVETA%VETAF(JLEV)-ZHVETAON)/(ZHVETAOX-YDVETA%VETAF(JLEV))
          ELSE
            ZZ=EXP(-RTDT*ZINEZ*ZWF(JROF,JLEV)*(1.0_JPRB/(PLEV(JROF,JLEV)-ZHVETAON)&
             & +1.0_JPRB/(ZHVETAOX-PLEV(JROF,JLEV))))&
             & *(YDVETA%VETAF(JLEV)-ZHVETAON)/(ZHVETAOX-YDVETA%VETAF(JLEV))
          ENDIF
          ZVDISP_2=(ZHVETAON+ZHVETAOX*ZZ)/(1.0_JPRB+ZZ)-YDVETA%VETAF(JLEV)
          ZVDISP=SIGN(1.0_JPRB,ZVDISP_1)*MIN(ABS(ZVDISP_1),ABS(ZVDISP_2))
          ZLEVO=YDVETA%VETAF(JLEV)+ZVDISP
        ELSE
          ZLEVO=YDVETA%VETAF(JLEV)-RTDT*ZINEZ*ZWF(JROF,JLEV)
        ENDIF
        ISTEST(JROF)=ISTEST(JROF)-MIN(0,MAX(-1,NINT(ZLEVO-ZLEVT-0.5_JPRB)))
        ISTESB(JROF)=ISTESB(JROF)-MIN(0,MAX(-1,NINT(ZLEVB-ZLEVO-0.5_JPRB)))
        ZLEVO=MIN(ZVETAOX,MAX(ZVETAON,ZLEVO))
        IF(LLO.OR.YRDYNA%LELTRA) THEN
          PLEV(JROF,JLEV)=ZLEVO
        ELSE
          PLEV(JROF,JLEV)=0.5_JPRB*(ZLEVO+YDVETA%VETAF(JLEV))
        ENDIF

        !later IF( LSLDIA.AND.(JITER==NITMP) ) THEN
        !later   ! trajectory vertical velocity
        !later   PWF(JROF,JLEV)=ZWF(JROF,JLEV)
        !later ENDIF

      ENDDO
    ENDDO

  ENDIF

  IF( ANY( ISTEST /= ISTEST0 ) .AND. JITER == NITMP) THEN
    WRITE(NULERR,'(A,I6,A)') ' SMILAG TRAJECTORY OUT OF ATM ',SUM(ISTEST(KST:KPROF)),' TIMES.'
    ! print statistics
    IF (YRDYNA%LRPRSLTRJ) THEN
      DO JROF=KST,KPROF
        IF( ISTEST(JROF) /= 0 ) THEN
          WRITE(NULERR,*) ' POINT= ',JROF,' MAX ETADOT VERTICAL VEL.= ',MAXVAL(ZWF(JROF,:))
          WRITE(NULERR,*) ' LON  = ',ACOS(YDCSGEOM%RCOLON(JROF))*180._JPRB/RPI,' degrees'
          WRITE(NULERR,*) ' LAT  = ',ASIN(YDGSGEOM%GEMU(JROF))*180._JPRB/RPI,' degrees'
        ENDIF
      ENDDO
    ENDIF
  ENDIF
  IF( ANY( ISTESB /= ISTEST0 ) .AND. JITER == NITMP) THEN
    WRITE(NULERR,'(A,I6,A)') ' SMILAG TRAJECTORY UNDERGROUND ',SUM(ISTESB(KST:KPROF)),' TIMES.'
    ! print statistics
    IF (YRDYNA%LRPRSLTRJ) THEN
      DO JROF=KST,KPROF
        IF(ISTESB(JROF) /= 0 ) THEN
          WRITE(NULERR,*) ' POINT= ',JROF,' MAX ETADOT VERTICAL VEL.= ',MAXVAL(ZWF(JROF,:))
          WRITE(NULERR,*) ' LON  = ',ACOS(YDCSGEOM%RCOLON(JROF))*180._JPRB/RPI,' degrees'
          WRITE(NULERR,*) ' LAT  = ',ASIN(YDGSGEOM%GEMU(JROF))*180._JPRB/RPI,' degrees'
        ENDIF
      ENDDO
    ENDIF
  ENDIF

  !*       2.2   DETERMINATION OF THE wind AT "M" OR "O".

  IF(JITER /= NITMP) THEN

    ! If (LLO=T or LELTRA=T)
    ! and JITER < NITMP wind is interpolated at
    ! the origin point "O" instead of at the medium point "M".
    ! In the other cases the wind is interpolated at "M".

    IHVI=0
    IROT=0
    ITIP=1
    LLINTV=.TRUE.

    CALL LARCINA(YDGEOMETRY,YDML_DYN,KST,KPROF,YDSL,IHVI,KSTABUF,LFINDVSEP,LLSLHD,LLSLHDQUAD,LLINTV,&
     & ITIP,IROT,.TRUE.,PLSDEPI,&
     & KIBL,PSCO,PLEV,&
     & ZUN_KAPPA,ZUN_KAPPAT,ZUN_KAPPAM,ZUN_KAPPAH,&
     & ZSLB2STDDISU,ZSLB2STDDISV,ZSLB2STDDISW,&
     & ZSLB1UR0,ZSLB1VR0,ZSLB1ZR0,ZSLB1WR0,&
     & KVSEPC,KVSEPL,PCCO,&
     & PUF,PVF,ZZF,ZWF,ZWFSM,&
     & KL0,KLH0,KLEV,PLSCAW,PRSCAW,IDEP,INOWENO)

  ENDIF

ENDDO

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF THE ORIGIN POINT "O" COORDINATES.
!              ------------------------------------------------

IF (.NOT.(LLO.OR.YRDYNA%LELTRA)) THEN
  ! this case may occur only for SL3TL scheme.
  IF (LTWOTL) CALL ABOR1(' ELARMES 3 ')

  ZVETAON=(1.0_JPRB-VETAON)*YDVETA%VETAH(0)+VETAON*YDVETA%VETAF(1)
  ZVETAOX=(1.0_JPRB-VETAOX)*YDVETA%VETAH(NFLEVG)+VETAOX*YDVETA%VETAF(NFLEVG)

  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF

      ! computations on horizontal plans.
      PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO)=2.0_JPRB*PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO)-YDCSGEOM%RINDY(JROF)
      PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO)=2.0_JPRB*PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO)-YDCSGEOM%RINDX(JROF)

      ! computations on a vertical.

      PLEV(JROF,JLEV)=2.0_JPRB*PLEV(JROF,JLEV)-YDVETA%VETAF(JLEV)
      PLEV(JROF,JLEV)=MAX(ZVETAON,MIN(ZVETAOX,PLEV(JROF,JLEV)))

    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ELARMES',1,ZHOOK_HANDLE)
END SUBROUTINE ELARMES
