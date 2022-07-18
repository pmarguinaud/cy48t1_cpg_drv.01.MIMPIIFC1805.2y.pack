SUBROUTINE ELARMESAD(YDGEOMETRY,YDRIP,YDML_DYN,KSTGLO,KST,KPROF,YDSL,KSTABUF,PGMDTX,PGMDTY,PB1,PB15,PB2,&
 & PINC,KINC,KDIM1,KDIM2,&
 & PLSDEPI,KIBL,KDEP,KNOWENO,&
 & PSCO,PLEV,PCCO,KL0,KLOCK,&
 & PSCO5,PLEV5,PCCO5,PLSCAW5,PRSCAW5,POUT5)

!**** *ELARMESAD - semi-LAgrangian scheme (adjoint code):
!                Research of the origin point on the Sphere.

!     Purpose.
!     --------

!      The computation of the location of the origin point "O" of
!     the lagrangian trajectory is performed by an iterative
!     method described by Robert and adapted to the sphere by M. Rochas.
!     Trajectories are great circles for spherical geometry.
!      In ELARMES(AD) the trajectories are straight lines on the plane.

!**   Interface.
!     ----------
!        *CALL* *ELARMESAD(...)

!        Explicit arguments :
!        --------------------

!          KSTGLO   - global offset into nproma buffer.
!          KST      - first element of arrays where computations are performed.
!          KPROF    - depth of work.
!          YDSL     - SL_STRUCT definition
!          KSTABUF  - for a latitude IGL, KSTABUF(IGL) is the
!                     address of the element corresponding to
!                     (ILON=1,IGL) in the NPROMA arrays.
!          PGMDTX   - m * DELTA t / DELTA x.
!          PGMDTY   - m * DELTA t / DELTA y.
!          PB1      - SLBUF1-buffer for interpolations.
!          PB15     - SLBUF1-buffer for interpolations (trajectory).
!          PB2      - SLBUF2-buf to communicate info from non lag. to lag. dyn.
!          KDIM1,KDIM2 - dimensions of KINC, PINC
!          KINC     - addresses of interpolation increments in global buffer
!          PINC     - interpolation increments to global buffer (vector only)
!          PLSDEPI  - (Number of points by latitude) / (2 * PI) .
!          KIBL     - index into YRGSGEOM/YRCSGEOM instances in YDGEOMETRY
!          KDEP     - dependency indicator for the southernost or northermost latitudes (LAM only)
!          KNOWENO  - vertical boundary indicator for the WENO interpolation.
!          PSCO     - information about geographic position of interpol. point.
!          PLEV     - vertical coordinate of the interpolation point.
!          PCCO     - information about comput. space position of interpol. point.
!          KL0      - index of the four western points
!                     of the 16 points interpolation grid.
!          KLOCK    - ??? (used for OpenMP purposes).

!        TRAJECTORY:
!          PSCO5    - cf. PSCO                                       (in)
!          PLEV5    - cf. PLEV                                       (in)
!          PCCO5    - cf. PCCO                                       (in)
!          PLSCAW5  - cf. PLSCAW                                     (in)
!          PRSCAW5  - cf. PRSCAW

!          POUT5    -  equal to zero if O inside the computational   (in)
!                      domain, otherwise contains:
!                      +1 increment if O is out of the C+I zone in x direction,
!                      +2 increment if O is out of the C+I zone in y direction,
!                      +4 increment if O is bellow or above the model atmosphere.
!                     Note: TL/AD od POUT is always zero. Hence not coded.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------
!        Calls  LARCINAAD.
!        Is called by LAPINEAAD (3D model)
!        Attention: some input dummies apparently relevant to the spherical
!        geometry are still present due to the call to LARCINAAD, even if they
!        are not really used.

!     Reference.
!     ----------

!     Author.
!     -------
!        F. Vana

!     Modifications.
!     --------------
!        Original : SEPTEMBER 2006
!        Modifications:
!        F. Vana 08-Jan-2007 KDEP array for vectorized adjoint
!        F. Vana 14-Aug-2008 arguments update for LARCINAAD
!        K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!        F. Vana :     11-Dec-2008 OpenMP for vector platforms
!        K. Yessad (Aug 2009): always use root (QX,QY) for (p,q) variables names
!        K. Yessad (Nov 2009): keep LLO.OR.LELTRA=T only.
!        K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!        G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!        G. Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
!        B. Bochenek (Apr 2015): Phasing: move some variables.
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!        R. El Khatib 22-May-2019 fix the "IFL OUT OF BOUNDS" issue
!        R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
!        P. Sekula 17-02-2020 Monkey business - change of size array KNOWENO
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK
USE SC2PRG_MOD, ONLY : SC2PRG
USE YOMRIP    , ONLY : TRIP

USE YOMDYNA   , ONLY : YRDYNA
USE EINT_MOD  , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN)       :: YDGEOMETRY
TYPE(TRIP)    , INTENT(IN)       :: YDRIP
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTGLO
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
TYPE(SL_STRUCT),   INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTABUF(YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMDTX(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMDTY(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1(YDSL%NASLB1,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB15(YDSL%NASLB1,YDML_DYN%YRPTRSLB15%NFLDSLB15)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIM1
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIM2
INTEGER(KIND=JPIM),INTENT(INOUT) :: KINC(KDIM1,KDIM2,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PINC(KDIM1,KDIM2,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
INTEGER(KIND=JPIM),INTENT(IN)    :: KDEP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)
INTEGER(KIND=JPIM),INTENT(IN)    :: KNOWENO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTSCO%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM)
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3,YDML_DYN%YRDYN%NITMP)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KLOCK(8)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTSCO%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLEV5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM,&
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POUT5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IROT, ITIP, JITER, JLEV, JROF, IINC

LOGICAL :: LLO
LOGICAL :: LLSLHD,LLSLHDQUAD,LLSLHD_OLD
LOGICAL :: LLOCK_UP

REAL(KIND=JPRB) :: ZDGUN(YDGEOMETRY%YRDIM%NPROMA), ZDGUX(YDGEOMETRY%YRDIM%NPROMA), ZDLUN, ZDLUX
REAL(KIND=JPRB) :: ZINDX, ZINDY, ZINEZ, ZINEZV(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZLEVO
REAL(KIND=JPRB) :: ZPU, ZPV
REAL(KIND=JPRB) :: ZTXO, ZTYO

REAL(KIND=JPRB) :: ZEPS
REAL(KIND=JPRB) :: ZVETAON, ZVETAOX

REAL(KIND=JPRB) :: ZUF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZVF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZZF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZWF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZLSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM)
REAL(KIND=JPRB) :: ZRSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM)

! unused arguments to LARCINAAD (they can be passed to LAPINEAAD
!  but there's no reason for it, as only the linear interpolation
!  is concerned at this stage)
! a) input arrays
REAL(KIND=JPRB) :: ZUN_KAPPA5(1), ZUN_KAPPAT5(1)
! b) output arrays - dimensions kept for safety
REAL(KIND=JPRB) :: ZUN_KAPPA(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZUN_KAPPAT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB), POINTER :: ZSLB1UR0(:), ZSLB1VR0(:), ZSLB1WR0(:), ZSLB1ZR0(:)
REAL(KIND=JPRB), POINTER :: ZSLB1UR05(:), ZSLB1VR05(:), ZSLB1WR05(:), ZSLB1ZR05(:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "larcinaad.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ELARMESAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 &  YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL), &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDPTRSLB15=>YDML_DYN%YRPTRSLB15, &
 & YDPTRSLB2=>YDML_DYN%YRPTRSLB2, &
 & YDPTRSLB1=>YDML_DYN%YRPTRSLB1,YDDYN=>YDML_DYN%YRDYN, &
 & YDVETA=>YDGEOMETRY%YRVETA)

ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, NPROMA=>YDDIM%NPROMA,NITMP=>YDDYN%NITMP, &
        & VETAON=>YDDYN%VETAON,VETAOX=>YDDYN%VETAOX, &
        & NVSEPC=>YDDYN%NVSEPC,NVSEPL=>YDDYN%NVSEPL, &
        & MSLB1UR0=>YDPTRSLB1%MSLB1UR0, MSLB1VR0=>YDPTRSLB1%MSLB1VR0, &
        & MSLB1WR0=>YDPTRSLB1%MSLB1WR0, MSLB1UR05=>YDPTRSLB15%MSLB1UR05, &
        & MSLB1VR05=>YDPTRSLB15%MSLB1VR05, MSLB1WR05=>YDPTRSLB15%MSLB1WR05, &
        & MSLB1ZR05=>YDPTRSLB15%MSLB1ZR05, MSLB1ZR0=>YDPTRSLB1%MSLB1ZR0, &
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

! in practical LLO.OR.LELTRA is now always T
LLO=.NOT.YRDYNA%LELTRA

IF(.NOT.(LLO.OR.YRDYNA%LELTRA)) THEN
  CALL ABOR1(' ELARMESAD: LLO.OR.LELTRA=F no longer supported')
ENDIF

ZDLUN=REAL(YDSL%NDLUNG,JPRB)
ZDLUX=REAL(YDSL%NDLUXG,JPRB)
ZEPS=1.E-6_JPRB
DO JROF=KST,KPROF
  ZINEZV(JROF)=MAX(0.0_JPRB,SIGN(1.0_JPRB,(YDCSGEOM%RINDX(JROF)-ZDLUN)*&
   & (ZDLUX-YDCSGEOM%RINDX(JROF))))*&
   & MAX(0.0_JPRB,SIGN(1.0_JPRB,(YDCSGEOM%RINDY(JROF)-ZDGUN(JROF))*(ZDGUX(JROF)-YDCSGEOM%RINDY(JROF))))  
ENDDO

DO JLEV=1,NFLEVG
  DO JROF=KST,KPROF
    ZUF(JROF,JLEV)=0.0_JPRB
    ZVF(JROF,JLEV)=0.0_JPRB
    ZZF(JROF,JLEV)=0.0_JPRB
    ZWF(JROF,JLEV)=0.0_JPRB
    PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RQX)=0.0_JPRB
    PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RQY)=0.0_JPRB
  ENDDO
ENDDO

ZVETAON=(1.0_JPRB-VETAON)*YDVETA%VETAH(0)+VETAON*YDVETA%VETAF(1)
ZVETAOX=(1.0_JPRB-VETAOX)*YDVETA%VETAH(NFLEVG)+VETAOX*YDVETA%VETAF(NFLEVG)

! deactivate computation of SLHD weights
LLSLHD     =.FALSE.
LLSLHDQUAD =.FALSE.
LLSLHD_OLD =.FALSE.

!     ------------------------------------------------------------------

!*       2.    ITERATIONS.       (Adjoint)
!              -----------

DO JITER=NITMP,1,-1

  !*       2.2   DETERMINATION OF THE "(a/rs)*wind" AT "O".

  IF(JITER /= NITMP) THEN
    LLOCK_UP=(JITER == 1)

    IROT=0
    ITIP=1

    ! Compute number of updates for LARCINAAD:
    ! 3D-linear interpolations, 8 each :
    IINC=24

    CALL LARCINAAD(YDGEOMETRY,YDML_DYN,IINC,KSTGLO,KST,KPROF,YDSL,&
     & KSTABUF,LLSLHD,LLSLHDQUAD,LLSLHD_OLD,ITIP,IROT,.TRUE.,LLOCK_UP,PLSDEPI,&
     & KIBL,KDEP(1,1,JITER),KNOWENO(1,1,JITER),&
     & PSCO,PLEV,PSCO5(1,1,1,JITER),PLEV5(1,1,JITER),&
     & ZUN_KAPPA,ZUN_KAPPA5,ZUN_KAPPAT,ZUN_KAPPAT5,&
     & ZSLB1UR0,ZSLB1VR0,ZSLB1ZR0,ZSLB1WR0,&
     & ZSLB1UR05,ZSLB1VR05,ZSLB1ZR05,ZSLB1WR05,&
     & NVSEPC,NVSEPL,KLOCK,&
     & PINC(1,1,JITER),KINC(1,1,JITER),KDIM1,KDIM2,&
     & PCCO,ZUF,ZVF,ZZF,ZWF,PCCO5(1,1,1,JITER),&
     & KL0(1,1,0,JITER),ZLSCAW,ZRSCAW,PLSCAW5(1,1,1,JITER),PRSCAW5(1,1,1,JITER))

  ENDIF

  !*       2.1   DETERMINATION OF THE ORIGIN POINT "O".

  IF (JITER == 1) THEN

    DO JLEV=1,NFLEVG
!DEC$ IVDEP
      DO JROF=KST,KPROF

        ZINDX=YDCSGEOM%RINDX(JROF)
        ZINDY=YDCSGEOM%RINDY(JROF)
        ZINEZ=ZINEZV(JROF)

        ! * computations on a vertical.

        IF(LLO.OR.YRDYNA%LELTRA) THEN
          ZLEVO=PLEV(JROF,JLEV)
        ENDIF
        PLEV(JROF,JLEV)=0.0_JPRB
        IF (POUT5(JROF,JLEV,1) >= 4._JPRB ) ZLEVO=0.0_JPRB
        PB2(JROF,MSLB2WRL+JLEV-1)=PB2(JROF,MSLB2WRL+JLEV-1)-RTDT*ZLEVO*ZINEZ
        ZLEVO=0.0_JPRB

        ! * computations on horizontal plans.

        !   - Fill the array elements by coordinates of the origin point
        !     Set the origin point to the arrival point if it is left or right
        !     out of C+I zone

        IF(LLO.OR.YRDYNA%LELTRA) THEN
          ZTXO=PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO)*ZINEZ
          ZTYO=PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO)*ZINEZ
        ENDIF
        PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO)=0.0_JPRB
        PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO)=0.0_JPRB

        !   - Return back the departure point if it is out of the C+I zone

        IF (MOD(POUT5(JROF,JLEV,1),2.0_JPRB) == 1.0_JPRB) THEN
          ZTXO=0.0_JPRB
        ENDIF
        IF (MOD(REAL(INT(POUT5(JROF,JLEV,1)/2.0_JPRB,JPIM),JPRB),2.0_JPRB)&
         & == 1.0_JPRB) THEN
          ZTYO=0.0_JPRB
        ENDIF

        !   - Compute the relative coordinates of departure point of trajectory

        ZTXO=ZTXO+PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLON)*ZINEZ
        ZTYO=ZTYO+PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLAT)*ZINEZ
        PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLON)=0.0_JPRB
        PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLAT)=0.0_JPRB
        PB2(JROF,MSLB2URL+JLEV-1)=PB2(JROF,MSLB2URL+JLEV-1)&
         & -2.0_JPRB*ZTXO*PGMDTX(JROF)
        PB2(JROF,MSLB2VRL+JLEV-1)=PB2(JROF,MSLB2VRL+JLEV-1)&
         & -2.0_JPRB*ZTYO*PGMDTY(JROF)
        ZTXO=0.0_JPRB
        ZTYO=0.0_JPRB

      ENDDO
    ENDDO

  ELSE

    DO JLEV=1,NFLEVG
!DEC$ IVDEP
      DO JROF=KST,KPROF

        ! * ZPU,ZPV are the coordinates of VM in the local repere
        !   related to F.

        ZINDX=YDCSGEOM%RINDX(JROF)
        ZINDY=YDCSGEOM%RINDY(JROF)
        ZINEZ=ZINEZV(JROF)

        ! * computations on a vertical.

        IF(LLO.OR.YRDYNA%LELTRA) THEN
          ZLEVO=PLEV(JROF,JLEV)
        ENDIF
        PLEV(JROF,JLEV)=0.0_JPRB
        IF (POUT5(JROF,JLEV,JITER) >= 4._JPRB ) ZLEVO=0.0_JPRB
        ZWF(JROF,JLEV)=-RTDT*ZLEVO*ZINEZ
        ZLEVO=0.0_JPRB
        IF(LLO) THEN
          PB2(JROF,MSLB2WRL+JLEV-1)=PB2(JROF,MSLB2WRL+JLEV-1)&
           & +0.5_JPRB*ZWF(JROF,JLEV)*ZINEZ
          ZWF(JROF,JLEV)=ZWF(JROF,JLEV)*0.5_JPRB
        ENDIF

        ! * horizontal computation

        !   - Fill the array elements by coordinates of the origin point
        !     Set the origin point to the arrival point if it is left or right
        !     out of C+I zone

        IF(LLO.OR.YRDYNA%LELTRA) THEN
          ZTXO=PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO)*ZINEZ
          ZTYO=PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO)*ZINEZ
        ENDIF
        PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO)=0.0_JPRB
        PSCO(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO)=0.0_JPRB

        !   - Return back the departure point if it is out of the C+I zone

        IF (MOD(POUT5(JROF,JLEV,JITER),2.0_JPRB) == 1.0_JPRB) THEN
          ZTXO=0.0_JPRB
        ENDIF
        IF (MOD(REAL(INT(POUT5(JROF,JLEV,JITER)/2.0_JPRB,JPIM),JPRB),&
         & 2.0_JPRB) == 1.0_JPRB) THEN
          ZTYO=0.0_JPRB
        ENDIF

        !   - Compute the relative coordinates of departure point of trajectory

        ZTXO=ZTXO+PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLON)*ZINEZ
        ZTYO=ZTYO+PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLAT)*ZINEZ
        PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLON)= 0.0_JPRB 
        PCCO(JROF,JLEV,YDML_DYN%YYTCCO%M_RLAT)= 0.0_JPRB

        ZPU=-2.0_JPRB*ZTXO*PGMDTX(JROF)
        ZPV=-2.0_JPRB*ZTYO*PGMDTY(JROF)
        ZTXO = 0.0_JPRB
        ZTYO = 0.0_JPRB

        IF(YRDYNA%LELTRA) THEN
          ZUF(JROF,JLEV)=ZPU
          ZVF(JROF,JLEV)=ZPV
        ELSEIF(LLO) THEN
          ZUF(JROF,JLEV)=0.5_JPRB*ZPU
          PB2(JROF,MSLB2URL+JLEV-1)=PB2(JROF,MSLB2URL+JLEV-1)+0.5_JPRB*ZPU
          ZVF(JROF,JLEV)=0.5_JPRB*ZPV
          PB2(JROF,MSLB2VRL+JLEV-1)=PB2(JROF,MSLB2VRL+JLEV-1)+0.5_JPRB*ZPV
        ENDIF
        ZPU= 0.0_JPRB
        ZPV= 0.0_JPRB

      ENDDO
    ENDDO

  ENDIF

ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ELARMESAD',1,ZHOOK_HANDLE)
END SUBROUTINE ELARMESAD
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
