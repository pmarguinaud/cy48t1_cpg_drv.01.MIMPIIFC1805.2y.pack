SUBROUTINE ELARMES5(YDGEOMETRY,YDML_DYN,YDRIP,KST,KPROF,YDSL,KSTABUF,PGMDTX,PGMDTY,PB15,PB2,&
 & PLSDEPI,KIBL,KDEP,KNOWENO,&
 & PSCO5,PLEV5,PCCO5,PUF5,PVF5,KL0,KLH0,KLEV,PLSCAW5,POUT5)

!**** *ELARMES5 - semi-LAgrangian scheme:
!                Research of the origin point on the Sphere.
!                Trajectory computation for TL/AD models.

!     Purpose.
!     --------

!      The computation of the location of the interpolation point of
!     the lagrangian trajectory is performed by an iterative
!     method described by Robert and adapted to the sphere by M. Rochas.
!     Trajectories are great circles for spherical geometry.
!      In ELARMES5 the trajectories are straight lines on the plane.

!**   Interface.
!     ----------
!        *CALL* *ELARMES5(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST       - first element of arrays where computations are performed.
!          KPROF     - depth of work.
!          YDSL      - SL_STRUCT definition.
!          KSTABUF   - for a latitude IGL, KSTABUF(IGL) is the
!                      address of the element corresponding to
!                      (ILON=1,IGL) in the NPROMA arrays.
!          PGMDTX    - m * DELTA t / DELTA x.
!          PGMDTY    - m * DELTA t / DELTA y.
!          PB15      - SLBUF1-buffer for interpolations (trajectory).
!          PB2       - SLBUF2-buf to communicate info from non lag. to lag. dyn.
!          PLSDEPI   - (Number of points by latitude) / (2 * PI) .
!          KIBL      - index into YRGSGEOM/YRCSGEOM instances in YDGEOMETRY

!        OUTPUT:
!          KDEP      - indication of the interpolation stencil latitudial
!                      dependences (used for LVECADIN=.T. option in adjoint)
!          KNOWENO   - indication of the specific treatment of vertical
!                      boundary used by WENO scheme (only usefull to pass upwards in adjoint code).
!          PSCO5     - information about geographic position of interpol. point.
!          PLEV5     - vertical coordinate of the interpolation point.
!          PCCO5     - information about comput. space position of interpol. point.
!          PUF5,PVF5 - U-comp and V-comp of "(a/rs)*wind" necessary to
!                      find the position of the origin point,
!                      in a local repere linked to computational sphere.
!          KL0       - index of the four western points
!                      of the 16 points interpolation grid.
!          KLH0      - second value of index of the four western points
!                      of the 16 points interpolation grid if needed.
!          KLEV      - lower level of the vertical interpolation
!                      grid needed for vertical interpolations.
!          PLSCAW5   - linear weights (distances) for interpolations.
!          POUT5     - equal zero if O inside the computational domain,
!                      otherwise contains:
!                      +1 increment if O is out of the C+I zone in x direction,
!                      +2 increment if O is out of the C+I zone in y direction,
!                      +4 increment if O is above or bellow atmosphere.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------
!       Called by LAPINEA5 

!     Reference.
!     ----------

!     Author.
!     -------
!        F. Vana after ELARMES and TL/AD code of IFS/ARPEGE

!     Modifications.
!     --------------
!        Original : MAY 2006
!        Modifications:
!        F. Vana  08-Jan-2007  new array KDEP to be used in adjoint
!        F. Vana  28-Aug-2007  one less dummy array to LARCINA
!        07-Nov-2007 J. Masek    Updated arguments of LARCINA.
!        F. Vana  14-Aug-2008  optimization of KDEP (not usefull
!                              for trajectory research)
!        K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!        K. Yessad (Aug 2009): always use root (QX,QY) for (p,q) variables names
!        K. Yessad (Nov 2009): keep LLO.OR.LELTRA=T only.
!        K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!        F. Vana  23-Feb-2011  update of arguments to LARCINA
!        G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!        G. Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
!        S. Malardel and D. Ricard (Nov 2013): COMAD weights for SL interpolations
!        B. Bochenek (Apr 2015): Phasing: move some variables.
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!        R. El Khatib 22-May-2019 fix the "IFL OUT OF BOUNDS" issue
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK
USE SC2PRG_MOD, ONLY : SC2PRG

USE YOMDYNA   , ONLY : YRDYNA
USE YOMCSGEOM , ONLY : TCSGEOM
USE YOMGSGEOM , ONLY : TGSGEOM
USE YOMRIP    , ONLY : TRIP
USE EINT_MOD  , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN)       :: YDGEOMETRY
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
TYPE(TRIP)     ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
TYPE(SL_STRUCT),   INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTABUF(YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMDTX(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMDTY(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB15(YDSL%NASLB1,YDML_DYN%YRPTRSLB15%NFLDSLB15)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
INTEGER(KIND=JPIM),INTENT(OUT)    :: KDEP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)
INTEGER(KIND=JPIM),INTENT(OUT)    :: KNOWENO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTSCO%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLEV5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP-1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP-1)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3,YDML_DYN%YRDYN%NITMP)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLH0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POUT5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ITIP
INTEGER(KIND=JPIM) :: IROT
INTEGER(KIND=JPIM) :: IHVI
INTEGER(KIND=JPIM) :: JITER
INTEGER(KIND=JPIM) :: JLEV
INTEGER(KIND=JPIM) :: JROF

LOGICAL :: LLFINDVSEP
LOGICAL :: LLINTV
LOGICAL :: LLO
LOGICAL :: LLSLHD, LLSLHDQUAD

REAL(KIND=JPRB) :: ZDGUN(YDGEOMETRY%YRDIM%NPROMA), ZDGUX(YDGEOMETRY%YRDIM%NPROMA), ZDLUN, ZDLUX
REAL(KIND=JPRB) :: ZINDX, ZINDY, ZINEZ, ZINEZV(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZLEVO5
REAL(KIND=JPRB) :: ZPU5, ZPV5
REAL(KIND=JPRB) :: ZTXO5, ZTYO5
REAL(KIND=JPRB) :: ZXPM, ZYPM, ZWPM

REAL(KIND=JPRB) :: ZEPS
REAL(KIND=JPRB) :: ZVETAON, ZVETAOX

REAL(KIND=JPRB) :: ZWF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
               &   ZWFSM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZZ05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

! unused arguments in call to LARCINA:
! a) input arrays - will not be used since LSLHD was set to .FALSE.
!    in LAPINEA before call to LARMES/ELARMES => dimensions can be
!    contracted
REAL(KIND=JPRB) :: ZUN_KAPPA(1),ZUN_KAPPAM(1),ZUN_KAPPAH(1),ZUN_KAPPAT(1)
REAL(KIND=JPRB) :: ZUN_STDDISU(1),ZUN_STDDISV(1),ZUN_STDDISW(1)
! b) input/output
INTEGER(KIND=JPIM) :: IUN_VSEPC
INTEGER(KIND=JPIM) :: IUN_VSEPL
! c) output arrays - dimensions kept for safety
REAL(KIND=JPRB) :: ZUN_RSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM)

REAL(KIND=JPRB), POINTER :: ZSLB1UR05(:), ZSLB1VR05(:), ZSLB1WR05(:), ZSLB1ZR05(:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "larcina.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ELARMES5',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDMP=>YDGEOMETRY%YRMP, YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL),YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL) , &
 & YDPTRSLB15=>YDML_DYN%YRPTRSLB15, &
 & YDPTRSLB2=>YDML_DYN%YRPTRSLB2,YDDYN=>YDML_DYN%YRDYN, &
 & YDVETA=>YDGEOMETRY%YRVETA)

ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, NPROMA=>YDDIM%NPROMA, NITMP=>YDDYN%NITMP, &
        & VETAON=>YDDYN%VETAON, VETAOX=>YDDYN%VETAOX, &
        & MSLB1UR05=>YDPTRSLB15%MSLB1UR05, MSLB1VR05=>YDPTRSLB15%MSLB1VR05, &
        & MSLB1ZR05=>YDPTRSLB15%MSLB1ZR05, &
        & MSLB1WR05=>YDPTRSLB15%MSLB1WR05, MSLB2URL5=>YDPTRSLB2%MSLB2URL5, &
        & MSLB2VRL5=>YDPTRSLB2%MSLB2VRL5, MSLB2WRL5=>YDPTRSLB2%MSLB2WRL5, &
        & RTDT=>YDRIP%RTDT)

CALL SC2PRG(MSLB1UR05 ,PB15     ,ZSLB1UR05)
CALL SC2PRG(MSLB1VR05 ,PB15     ,ZSLB1VR05)
CALL SC2PRG(MSLB1WR05 ,PB15     ,ZSLB1WR05)
CALL SC2PRG(MSLB1ZR05 ,PB15     ,ZSLB1ZR05)
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

!*       1.2   Miscellaneous preliminary initialisations.

POUT5(:,:,:)=0.0_JPRB

ZDLUN=REAL(YDSL%NDLUNG,JPRB)
ZDLUX=REAL(YDSL%NDLUXG,JPRB)
ZEPS=1.E-6_JPRB
DO JROF=KST,KPROF
  ZINEZV(JROF)=MAX(0.0_JPRB,SIGN(1.0_JPRB,(YDCSGEOM%RINDX(JROF)-ZDLUN)*&
   & (ZDLUX-YDCSGEOM%RINDX(JROF))))*&
   & MAX(0.0_JPRB,SIGN(1.0_JPRB,(YDCSGEOM%RINDY(JROF)-ZDGUN(JROF))*(ZDGUX(JROF)-YDCSGEOM%RINDY(JROF))))  
ENDDO

! in practical LLO.OR.LELTRA is now always T
LLO=.NOT.YRDYNA%LELTRA

IF(.NOT.(LLO.OR.YRDYNA%LELTRA)) THEN
  CALL ABOR1(' ELARMES5: LLO.OR.LELTRA=F no longer supported')
ENDIF

! deactivate computation of SLHD weights
LLSLHD     =.FALSE.
LLSLHDQUAD =.FALSE.

!     ------------------------------------------------------------------

!*       2.    ITERATIONS.
!              -----------

DO JITER=1,NITMP

  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      PCCO5(JROF,JLEV,YDML_DYN%YYTCCO%M_RQX,JITER)=1.0_JPRB
      PCCO5(JROF,JLEV,YDML_DYN%YYTCCO%M_RQY,JITER)=0.0_JPRB
    ENDDO
  ENDDO

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

        PUF5(JROF,JLEV,1)=PB2(JROF,MSLB2URL5+JLEV-1)
        PVF5(JROF,JLEV,1)=PB2(JROF,MSLB2VRL5+JLEV-1)
        ZTXO5 = ZINDX-2.0_JPRB*PB2(JROF,MSLB2URL5+JLEV-1)*PGMDTX(JROF)
        ZTYO5 = ZINDY-2.0_JPRB*PB2(JROF,MSLB2VRL5+JLEV-1)*PGMDTY(JROF)
        PCCO5(JROF,JLEV,YDML_DYN%YYTCCO%M_RLON,1) = ZTXO5*ZINEZ +ZINDX*(1.0_JPRB-ZINEZ)
        PCCO5(JROF,JLEV,YDML_DYN%YYTCCO%M_RLAT,1) = ZTYO5*ZINEZ +ZINDY*(1.0_JPRB-ZINEZ)

        !   - Return back the departure point if it is out of the C+I zone

        ZXPM = ZTXO5
        ZYPM = ZTYO5
        ZTXO5= MIN(MAX(ZTXO5,ZDLUN),ZDLUX)
        ZTYO5= MIN(MAX(ZTYO5,ZDGUN(JROF)),ZDGUX(JROF))
        IF(ZXPM /= ZTXO5) THEN
          POUT5(JROF,JLEV,1)=1.0_JPRB + POUT5(JROF,JLEV,1)
        ENDIF
        IF(ZYPM /= ZTYO5) THEN
          POUT5(JROF,JLEV,1)=2.0_JPRB + POUT5(JROF,JLEV,1)
        ENDIF

        !   - Fill the array elements by coordinates of the origin point
        !     Set the origin point to the arrival point if it is left or right
        !     out of C+I zone

        IF(LLO.OR.YRDYNA%LELTRA) THEN
          PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO,1) = ZTXO5*ZINEZ+ZINDX*(1.0_JPRB-ZINEZ)
          PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO,1) = ZTYO5*ZINEZ+ZINDY*(1.0_JPRB-ZINEZ)
        ENDIF

        ! * computations on a vertical.

        ZWF5(JROF,JLEV)=PB2(JROF,MSLB2WRL5+JLEV-1)
        ZLEVO5=YDVETA%VETAF(JLEV)-RTDT*ZWF5(JROF,JLEV)*ZINEZ
        ZWPM=ZLEVO5
        ZLEVO5=MIN(ZVETAOX,MAX(ZVETAON,ZLEVO5))
        IF(ZWPM /= ZLEVO5) POUT5(JROF,JLEV,1)=POUT5(JROF,JLEV,1)+4._JPRB
        IF(LLO.OR.YRDYNA%LELTRA) THEN
          PLEV5(JROF,JLEV,1)=ZLEVO5
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
          ZPU5=PUF5(JROF,JLEV,JITER-1)
          ZPV5=PVF5(JROF,JLEV,JITER-1)
        ELSEIF(LLO) THEN
          ZPU5=0.5_JPRB*(PUF5(JROF,JLEV,JITER-1)+PB2(JROF,MSLB2URL5+JLEV-1))
          ZPV5=0.5_JPRB*(PVF5(JROF,JLEV,JITER-1)+PB2(JROF,MSLB2VRL5+JLEV-1))
          ZWF5(JROF,JLEV)=&
           & 0.5_JPRB*(ZWF5(JROF,JLEV)+PB2(JROF,MSLB2WRL5+JLEV-1)*ZINEZ)
        ENDIF

        ZTXO5 = ZINDX-2.0_JPRB*ZPU5*PGMDTX(JROF)
        ZTYO5 = ZINDY-2.0_JPRB*ZPV5*PGMDTY(JROF)

        PCCO5(JROF,JLEV,YDML_DYN%YYTCCO%M_RLON,JITER) = ZTXO5*ZINEZ +ZINDX*(1.0_JPRB-ZINEZ)
        PCCO5(JROF,JLEV,YDML_DYN%YYTCCO%M_RLAT,JITER) = ZTYO5*ZINEZ +ZINDY*(1.0_JPRB-ZINEZ)

        !   - Return back the departure point if it is out of the C+I zone

        ZXPM = ZTXO5
        ZYPM = ZTYO5
        ZTXO5 = MIN(MAX(ZTXO5,ZDLUN),ZDLUX)
        ZTYO5 = MIN(MAX(ZTYO5,ZDGUN(JROF)),ZDGUX(JROF))
        IF(ZXPM /= ZTXO5) THEN
          POUT5(JROF,JLEV,JITER)=1.0_JPRB + POUT5(JROF,JLEV,JITER)
        ENDIF
        IF(ZYPM /= ZTYO5) THEN
          POUT5(JROF,JLEV,JITER)=2.0_JPRB + POUT5(JROF,JLEV,JITER)
        ENDIF

        !   - Fill the array elements by coordinates of the origin point
        !     Set the origin point to the arrival point if it is left or right
        !     out of C+I zone

        IF(LLO.OR.YRDYNA%LELTRA) THEN
          PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO,JITER) = ZTXO5*ZINEZ+ZINDX*(1.0_JPRB-ZINEZ)
          PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO,JITER) = ZTYO5*ZINEZ+ZINDY*(1.0_JPRB-ZINEZ)
        ENDIF

        ! * computations on a vertical.

        ZLEVO5=YDVETA%VETAF(JLEV)-RTDT*ZWF5(JROF,JLEV)*ZINEZ
        ZWPM=ZLEVO5
        ZLEVO5=MIN(ZVETAOX,MAX(ZVETAON,ZLEVO5))
        IF(ZWPM /= ZLEVO5) POUT5(JROF,JLEV,JITER)=POUT5(JROF,JLEV,JITER)+4._JPRB
        IF(LLO.OR.YRDYNA%LELTRA) THEN
          PLEV5(JROF,JLEV,JITER)=ZLEVO5
        ENDIF

      ENDDO
    ENDDO

  ENDIF

  !*       2.2   DETERMINATION OF THE "(a/rs)*wind" AT "O".

  IF(JITER /= NITMP) THEN

    IHVI=0
    IROT=0
    ITIP=1
    LLFINDVSEP=.FALSE.
    LLINTV=.TRUE.

    CALL LARCINA(&
     & YDGEOMETRY,YDML_DYN,KST,KPROF,YDSL,IHVI,KSTABUF,LLFINDVSEP,LLSLHD,LLSLHDQUAD,LLINTV,&
     & ITIP,IROT,.TRUE.,PLSDEPI,&
     & KIBL,&
     & PSCO5(1,1,1,JITER),PLEV5(1,1,JITER),&
     & ZUN_KAPPA,ZUN_KAPPAT,ZUN_KAPPAM,ZUN_KAPPAH,&
     & ZUN_STDDISU,ZUN_STDDISV,ZUN_STDDISW,&
     & ZSLB1UR05,ZSLB1VR05,ZSLB1ZR05,ZSLB1WR05,&
     & IUN_VSEPC,IUN_VSEPL,PCCO5(1,1,1,JITER),&
     & PUF5(1,1,JITER),PVF5(1,1,JITER),ZZ05,ZWF5,ZWFSM,&
     & KL0(1,1,0,JITER),KLH0(1,1,0),KLEV,PLSCAW5(1,1,1,JITER),ZUN_RSCAW5, & 
     & KDEP(1,1,JITER),KNOWENO(1,1,JITER))

  ENDIF

ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ELARMES5',1,ZHOOK_HANDLE)
END SUBROUTINE ELARMES5
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
