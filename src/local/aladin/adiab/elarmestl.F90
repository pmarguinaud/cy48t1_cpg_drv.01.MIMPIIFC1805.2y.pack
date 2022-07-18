SUBROUTINE ELARMESTL(YDGEOMETRY,YDML_DYN,YDRIP,KST,KPROF,YDSL,KSTABUF,PGMDTX,PGMDTY,PB1,PB15,PB2,&
 & PLSDEPI,KIBL,&
 & PSCO,PLEV,PCCO,&
 & KL0,KLH0,KLEV,PLSCAW,PRSCAW,&
 & PSCO5,PLEV5,PCCO5,PLSCAW5,PRSCAW5,POUT5)

!**** *ELARMESTL - semi-LAgrangian scheme:
!                Research of the origin point on the Sphere.

!     Purpose.
!     --------

!      The computation of the location of the origin point "O" of
!     the lagrangian trajectory is performed by an iterative
!     method described by Robert and adapted to the sphere by M. Rochas.
!     Trajectories are great circles for spherical geometry.
!      In ELARMESTL the trajectories are straight lines on the plane.

!**   Interface.
!     ----------
!        *CALL* *ELARMESTL(...)

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
!          PB15     - SLBUF1-buffer for interpolations (trajectory).
!          PB2      - SLBUF2-buf to communicate info from non lag. to lag. dyn.
!          PLSDEPI  - (Number of points by latitude) / (2 * PI) .
!          KIBL     - index into YRCSGEOM/YRGSGEOM instances in YDGEOMETRY
!          YDGSGEOM - grid point geometry.
!          YDCSGEOM - computational sphere geometry.

!        OUTPUT:
!          PSCO     - information about geographic position of interpol. point.
!          PLEV     - vertical coordinate of the interpolation point.
!          PCCO     - information about comput. space position of interpol. point.
!          KL0      - index of the four western points
!                     of the 16 points interpolation grid.
!          KLH0     - second value of index of the four western points
!                     of the 16 points interpolation grid if needed.
!          KLEV     - lower level of the vertical interpolation
!                     grid needed for vertical interpolations.
!          PLSCAW   - linear weights (distances) for interpolations.
!          PRSCAW   - non-linear weights for interpolations.

!        TRAJECTORY:
!          PSCO5    - cf. PSCO                                       (in)
!          PLEV5    - cf. PLEV                                       (in)
!          PCCO5    - cf. PCCO                                       (out)
!          PLSCAW5  - linear weights (distances) for interpolations  (out)
!          PRSCAW5  - non-linear weights for interpolations.         (out)
!          POUT5    - equal to zero if O inside the computational    (in)
!                     domain, otherwise contains:
!                      +1 increment if O is out of the C+I zone in x direction,
!                      +2 increment if O is out of the C+I zone in y direction,
!                      +4 increment if O is bellow or above the model atmosphere.
!                     (As derivatives of POUT are zero, its TL/AD is not coded.)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------
!        Calls  LARCINATL.
!        Is called by LAPINEATL (3D model)
!        Attention: some input dummies apparently relevant to the spherical
!        geometry are still present due to the call to LARCINATL, even if they
!        are not really used.

!     Reference.
!     ----------

!     Author.
!     -------
!        F. Vana

!     Modifications.
!     --------------
!        Original : MAY 2006
!        Modified 08-Jun-27 by F. Vana: update of arguments for LARCINATL
!        K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!        K. Yessad (Aug 2009): always use root (QX,QY) for (p,q) variables names
!        K. Yessad (Nov 2009): keep LLO.OR.LELTRA=T only.
!        K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!        G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!        G. Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
!        B. Bochenek (Apr 2015): Phasing: move some variables.
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
! End Modifications
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK
USE SC2PRG_MOD, ONLY : SC2PRG

USE YOMDYNA   , ONLY : YRDYNA
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB15(YDSL%NASLB1,YDML_DYN%YRPTRSLB15%NFLDSLB15)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTSCO%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3,YDML_DYN%YRDYN%NITMP)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLH0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTSCO%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLEV5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POUT5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IROT, ITIP, JITER, JLEV, JROF

LOGICAL :: LLO
LOGICAL :: LLSLHD,LLSLHDQUAD,LLSLHD_OLD

REAL(KIND=JPRB) :: ZDGUN, ZDGUX, ZDLUN, ZDLUX
REAL(KIND=JPRB) :: ZINDX, ZINDY, ZINEZ, ZINEZV(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZLEVO
REAL(KIND=JPRB) :: ZPU, ZPV
REAL(KIND=JPRB) :: ZTXO, ZTYO

REAL(KIND=JPRB) :: ZEPS
REAL(KIND=JPRB) :: ZVETAON, ZVETAOX

REAL(KIND=JPRB) :: ZUF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZVF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZWF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZZF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

! unused arguments to LARCINATL (they can be passed to LAPINEATL
!  but there's no reason for it, as the linear interpolation only
!  is concerned at this stage)
REAL(KIND=JPRB) :: ZUN_KAPPA(1),ZUN_KAPPAT(1),ZUN_KAPPAT5(1)
REAL(KIND=JPRB) :: ZUN_KAPPA5(1)

REAL(KIND=JPRB), POINTER :: ZSLB1UR0(:), ZSLB1VR0(:), ZSLB1WR0(:), ZSLB1ZR0(:)
REAL(KIND=JPRB), POINTER :: ZSLB1UR05(:), ZSLB1VR05(:), ZSLB1WR05(:), ZSLB1ZR05(:)

!SM: temporary fix I hope... (for call to LARCINATL)
INTEGER(KIND=JPIM) :: INOWENO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "larcinatl.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ELARMESTL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL),YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDPTRSLB15=>YDML_DYN%YRPTRSLB15, &
 & YDPTRSLB2=>YDML_DYN%YRPTRSLB2, &
 & YDPTRSLB1=>YDML_DYN%YRPTRSLB1,YDDYN=>YDML_DYN%YRDYN, &
 & YDVETA=>YDGEOMETRY%YRVETA)

ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, NPROMA=>YDDIM%NPROMA, NITMP=>YDDYN%NITMP, &
        & VETAON=>YDDYN%VETAON, VETAOX=>YDDYN%VETAOX, &
        & MSLB1UR0=>YDPTRSLB1%MSLB1UR0, MSLB1VR0=>YDPTRSLB1%MSLB1VR0, &
        & MSLB1WR0=>YDPTRSLB1%MSLB1WR0, MSLB1UR05=>YDPTRSLB15%MSLB1UR05, &
        & MSLB1VR05=>YDPTRSLB15%MSLB1VR05, MSLB1WR05=>YDPTRSLB15%MSLB1WR05, &
        & MSLB1ZR0=>YDPTRSLB1%MSLB1ZR0, MSLB1ZR05=>YDPTRSLB15%MSLB1ZR05, &
        & MSLB2URL=>YDPTRSLB2%MSLB2URL, MSLB2VRL=>YDPTRSLB2%MSLB2VRL, &
        & MSLB2WRL=>YDPTRSLB2%MSLB2WRL, RTDT=>YDRIP%RTDT)

CALL SC2PRG(MSLB1UR05 ,PB15     ,ZSLB1UR05)
CALL SC2PRG(MSLB1VR05 ,PB15     ,ZSLB1VR05)
CALL SC2PRG(MSLB1WR05 ,PB15     ,ZSLB1WR05)
CALL SC2PRG(MSLB1ZR05 ,PB15     ,ZSLB1ZR05)
CALL SC2PRG(MSLB1UR0  ,PB1      ,ZSLB1UR0)
CALL SC2PRG(MSLB1VR0  ,PB1      ,ZSLB1VR0)
CALL SC2PRG(MSLB1WR0  ,PB1      ,ZSLB1WR0)
CALL SC2PRG(MSLB1ZR0  ,PB1      ,ZSLB1ZR0)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS AND TESTS.
!              --------------------------------------

!*       1.2   Miscellaneous preliminary initialisations.

! in practical LLO.OR.LELTRA is now always T
LLO=.NOT.YRDYNA%LELTRA

IF(.NOT.(LLO.OR.YRDYNA%LELTRA)) THEN
  CALL ABOR1(' ELARMESTL: LLO.OR.LELTRA=F no longer supported')
ENDIF

! deactivate computation of SLHD weights
LLSLHD     =.FALSE.
LLSLHDQUAD =.FALSE.
LLSLHD_OLD =.FALSE.

ZDLUN=REAL(YDSL%NDLUNG,JPRB)
ZDLUX=REAL(YDSL%NDLUXG,JPRB)
ZDGUN=REAL(YDSL%NDGUNG,JPRB)
ZDGUX=REAL(YDSL%NDGUXG,JPRB)
ZEPS=1.E-6_JPRB
DO JROF=KST,KPROF
  ZINEZV(JROF)=MAX(0.0_JPRB,SIGN(1.0_JPRB,(YDCSGEOM%RINDX(JROF)-ZDLUN)*&
   & (ZDLUX-YDCSGEOM%RINDX(JROF))))*&
   & MAX(0.0_JPRB,SIGN(1.0_JPRB,(YDCSGEOM%RINDY(JROF)-ZDGUN)*(ZDGUX-YDCSGEOM%RINDY(JROF))))  
ENDDO

DO JLEV=1,NFLEVG
  DO JROF=KST,KPROF
    PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RQX)=0.0_JPRB
    PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RQY)=0.0_JPRB
  ENDDO
ENDDO

!SM
INOWENO(:,:)=0

!     ------------------------------------------------------------------

!*       2.    ITERATIONS.
!              -----------

DO JITER=1,NITMP

  !*       2.1   DETERMINATION OF THE ORIGIN POINT "O".

  ZVETAON=(1.0_JPRB-VETAON)*YDVETA%VETAH(0)+VETAON*YDVETA%VETAF(1)
  ZVETAOX=(1.0_JPRB-VETAOX)*YDVETA%VETAH(NFLEVG)+VETAOX*YDVETA%VETAF(NFLEVG)

  IF (JITER == 1) THEN

    DO JLEV=1,NFLEVG
      DO JROF=KST,KPROF

        ! * computations on horizontal plans.

        ZINDX=YDCSGEOM%RINDX(JROF)
        ZINDY=YDCSGEOM%RINDY(JROF)
        ZINEZ=ZINEZV(JROF)

        !   - Compute the relative coordinates of departure point of trajectory

        ZUF(JROF,JLEV)=PB2(JROF,MSLB2URL+JLEV-1)
        ZVF(JROF,JLEV)=PB2(JROF,MSLB2VRL+JLEV-1)
        ZTXO = -2.0_JPRB*PB2(JROF,MSLB2URL+JLEV-1)*PGMDTX(JROF)
        ZTYO = -2.0_JPRB*PB2(JROF,MSLB2VRL+JLEV-1)*PGMDTY(JROF)
        PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLON) = ZTXO*ZINEZ
        PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLAT) = ZTYO*ZINEZ

        !   - Return back the departure point if it is out of the C+I zone
        !     Note: TL of POUT is always 0. Hence not coded

        IF (MOD(POUT5(JROF,JLEV,1),2.0_JPRB) == 1.0_JPRB) THEN
          ZTXO=0.0_JPRB
        ENDIF
        IF (MOD(REAL(INT(POUT5(JROF,JLEV,1)/2.0_JPRB,JPIM),JPRB),2.0_JPRB)&
         & == 1.0_JPRB) THEN
          ZTYO=0.0_JPRB
        ENDIF

        !   - Fill the array elements by coordinates of the origin point
        !     Set the origin point to the arrival point if it is left or right
        !     out of C+I zone

        IF(LLO.OR.YRDYNA%LELTRA) THEN
          PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO) = ZTXO*ZINEZ
          PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO) = ZTYO*ZINEZ
        ENDIF

        ! * computations on a vertical.

        ZWF(JROF,JLEV)=PB2(JROF,MSLB2WRL+JLEV-1)
        ZLEVO=-RTDT*ZWF(JROF,JLEV)*ZINEZ
        IF (POUT5(JROF,JLEV,1) >= 4._JPRB ) ZLEVO=0.0_JPRB
        IF(LLO.OR.YRDYNA%LELTRA) THEN
          PLEV(JROF,JLEV)=ZLEVO
        ENDIF

      ENDDO
    ENDDO

  ELSE

    DO JLEV=1,NFLEVG
      DO JROF=KST,KPROF

        ! * ZPU,ZPV are the coordinates of VM in the local repere
        !   related to F.

        ZINDX=YDCSGEOM%RINDX(JROF)
        ZINDY=YDCSGEOM%RINDY(JROF)
        ZINEZ=ZINEZV(JROF)

        !   - Compute the relative coordinates of departure point of trajectory

        IF(YRDYNA%LELTRA) THEN
          ZPU=ZUF(JROF,JLEV)
          ZPV=ZVF(JROF,JLEV)
        ELSEIF(LLO) THEN
          ZPU=0.5_JPRB*(ZUF(JROF,JLEV)+PB2(JROF,MSLB2URL+JLEV-1))
          ZPV=0.5_JPRB*(ZVF(JROF,JLEV)+PB2(JROF,MSLB2VRL+JLEV-1))
          ZWF(JROF,JLEV)=&
           & 0.5_JPRB*(ZWF(JROF,JLEV)+PB2(JROF,MSLB2WRL+JLEV-1)*ZINEZ)
        ENDIF

        ZTXO = -2.0_JPRB*ZPU*PGMDTX(JROF)
        ZTYO = -2.0_JPRB*ZPV*PGMDTY(JROF)

        PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLON) = ZTXO*ZINEZ
        PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLAT) = ZTYO*ZINEZ

        !   - Return back the departure point if it is out of the C+I zone

        IF (MOD(POUT5(JROF,JLEV,JITER),2.0_JPRB) == 1.0_JPRB) THEN
          ZTXO=0.0_JPRB
        ENDIF
        IF (MOD(REAL(INT(POUT5(JROF,JLEV,JITER)/2.0_JPRB,JPIM),JPRB),&
         & 2.0_JPRB) == 1.0_JPRB) THEN
          ZTYO=0.0_JPRB
        ENDIF

        !   - Fill the array elements by coordinates of the origin point
        !     Set the origin point to the arrival point if it is left or right
        !     out of C+I zone

        IF(LLO.OR.YRDYNA%LELTRA) THEN
          PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO) = ZTXO*ZINEZ
          PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO) = ZTYO*ZINEZ
        ENDIF

        ! * save (zpu,zpv) in (puf,pvf) with a rotation equal to identity
        !   in plane geometry.
        ! Note: Since not needed afterwards this is useless for TL code

        ! * computations on a vertical.

        ZLEVO=-RTDT*ZWF(JROF,JLEV)*ZINEZ
        IF (POUT5(JROF,JLEV,JITER) >= 4._JPRB ) ZLEVO=0.0_JPRB
        IF(LLO.OR.YRDYNA%LELTRA) THEN
          PLEV(JROF,JLEV)=ZLEVO
        ENDIF

      ENDDO
    ENDDO

  ENDIF

  !*       2.2   DETERMINATION OF THE "(a/rs)*wind" AT "O".

  IF(JITER /= NITMP) THEN

    IROT=0
    ITIP=1

    CALL LARCINATL(YDGEOMETRY,YDML_DYN,KST,KPROF,YDSL,KSTABUF,LLSLHD,LLSLHDQUAD,LLSLHD_OLD,ITIP,IROT,.TRUE.,&
     & PLSDEPI,KIBL,&
     & PSCO,PLEV,PSCO5(1,1,1,JITER),PLEV5(1,1,JITER),&
     & ZUN_KAPPA,ZUN_KAPPA5,ZUN_KAPPAT,ZUN_KAPPAT5,&
     & ZSLB1UR0,ZSLB1VR0,ZSLB1ZR0,ZSLB1WR0,&
     & ZSLB1UR05,ZSLB1VR05,ZSLB1ZR05,ZSLB1WR05,&
     & PCCO,ZUF,ZVF,ZZF,ZWF,PCCO5(1,1,1,JITER),&
     & KL0(1,1,0,JITER),KLH0,KLEV,PLSCAW,PRSCAW,&
     & PLSCAW5(1,1,1,JITER),PRSCAW5,INOWENO)

  ENDIF

ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ELARMESTL',1,ZHOOK_HANDLE)
END SUBROUTINE ELARMESTL
