SUBROUTINE LAPINEA5(YDGEOMETRY,YDRIP,YDML_DYN,KST,KPROF,YDSL,LD_AD,&
 & KIBL,KSTAGE2,KSTAGE3,&
 & PB15,PB2,PCCO5,&
 & PUF5,PVF5,PZF5,PUF05,PVF05,PZF05,PWF05,PUS5,PVS5,PZS5,PWS5,&
 & KL0,KLH0,PLSCAW5,PRSCAW5,PSCO5,PLEV5,POUT5,KDEP,KNOWENO)

!**** *LAPINEA5* - semi-LAgrangian scheme:(Trajectory)
!                Interface subroutine for interpolations and lagrangian
!                trends. (Programme INterface d'Ensemble).

!     Purpose.
!     --------
!           Computes the semi-Lagrangian trajectory.
!           Computes the weights ready for interpolations at origin point.
!           Trajectory version of LAPINEA.

!**   Interface.
!     ----------
!        *CALL* *LAPINEA5(.....)

!        Explicit arguments :
!        --------------------
!        INPUT:
!          KST      - first element of arrays where computations are performed.
!          KPROF    - depth of work.
!          YDSL     - SL_STRUCT definition.
!          LD_AD    - T: LAPINEA5 is called in the adjoint model.
!                     F: LAPINEA5 is called in the linear tangent model.
!          KIBL     - index into YDGSGEOM/YDCSGEOM instances in YDGEOMETRY
!          KSTAGE2/KSTAGE3 - dimension for 2 and 3 stage variables (used by RK4)
!          PB15     - SLBUF1-buf traj for interpolations
!          PB2      - SLBUF2-buf to communicate info from non lag. to lag. dyn.

!        OUTPUT:

!          PCCO5    - information about comput. space position of interpol. point.
!          PUF5     - U-comp of "(a/rs)*wind" necessary to
!                     find the position of the origin point,
!                     in a local repere linked to computational sphere.
!          PVF5     - V-comp of "(a/rs)*wind" necessary to
!                     find the position of the origin point,
!                     in a local repere linked to computational sphere.
!          PZF5     - Z-comp of "(a/rs)*wind" (if required) necessary to
!                     find the position of the origin point,
!                     in a local repere linked to computational sphere.
!          NOTE: To distinguish between default 2nd order scheme and 4th order RK4 exclusive output
!          variables the stage dimension is swapped with number of iterations for the latter. 
!          PUF05,PVF05,PZF05,PWF05 - interpolated quantities for RK4 scheme
!          PUS5,PVS5,PZS5,PWS5 - stage winds for RK4 scheme
!          KL0      - index of the four western points
!                     of the 16 points interpolation grid.
!          KLH0     - second value of index of the four western points
!                     of the 16 points interpolation grid if needed.
!          PLSCAW5  - linear weights (distances) for interpolations (trajectory).
!          PRSCAW5  - non-linear weights for interpolations.
!          PSCO5    - information about geographic position of interpol. point.
!          PLEV5    - vertical coordinate of the interpolation point.
!          POUT5    - equal to 0 if O lays inside the computational area,
!                     otherwise contains sum of increments for each
!                     direction for which it is out of domain.
!          KDEP     - indication of the interpolation stencil latitudial
!                     dependences (used for LVECADIN=.T. option in adjoint)
!          KNOWENO  - vertical boundary indicator for the WENO interpolation.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!       See includes below.
!       Called by CALL_SL_AD and CALL_SL_TL.

!     Reference.
!     ----------
!        ARPEGE documentation about semi-Lagrangian scheme.

!     Author.
!     -------
!      C. Temperton (ECMWF)
!      Original : 99/10/15

!     Modifications.
!     --------------
!      Modified 01-06-26 by K. YESSAD: update dummy arguments of LAIDEP.
!      Modified 01-08-30 by K. YESSAD: pruning and some other cleanings.
!      R. El Khatib : 01-08-07 Pruning options
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Modified 05-01-11 by C. Temperton: streamlining.
!      D.Salmond     08-Feb-2005 SLADJTR replaced by LARMES5
!      F. Vana       Nov-2004  one more dummy to LAIDEP
!      K. Yessad     07-Sep-2005 Now calls directly LARCINA.
!      F. Vana       01-Jun-2006 LRPLANE=.T.
!      F. Vana       08-Jan-2007  fix + vectorization in LAM adjoint 
!      F. Vana       Nov-2004  one less dummy to LARCINA
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      30-Jun-2008 J. Masek New arg PCLO5, updated arguments of LARCINA.
!      F. Vana       14-Aug-2008 Weights driven interpolation + KDEP optim.
!      K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!      F. Vana       13-Jan-2009: TL/AD of SLHD
!      K. Yessad (Aug 2009): always use root (QX,QY) for (p,q) variables names
!      K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!      F. Vana       21-Feb-2011: Extended dimensions of weights
!      G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      G. Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
!      G. Mozdzynski (May 2012): further cleaning
!      F. Vana  13-Feb-2014  KappaT for heat quantities.
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      S. Malardel (Nov 2013): COMAD weights for SL interpolations
!      F. Vana    21-Nov-2017: Options LSLDP_CURV and LHOISLT
!      F. Vana July 2018: RK4 scheme for trajectory research.
!      F. Vana October 2018: Extended LSLDP_CURV.
!      F. Vana 26-Feb-2019:  Vertical quintic interpolation
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
!      F. Vana 18-Mar-2019: Implicit RK4 for trajectory research.
! End Modifications
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE SC2PRG_MOD         , ONLY : SC2PRG
USE YOMCT0             , ONLY : LRPLANE
USE YOMDYNA            , ONLY : YRDYNA
USE YOMCST             , ONLY : RPI
USE YOMRIP             , ONLY : TRIP
USE EINT_MOD           , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
LOGICAL           ,INTENT(IN)    :: LD_AD
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTAGE2
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTAGE3
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB15(YDSL%NASLB1,YDML_DYN%YRPTRSLB15%NFLDSLB15)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM, &
 &                                        YDML_DYN%YRDYN%NITMP,KSTAGE3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP-1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP-1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP-1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUF05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KSTAGE2,&
  &                                       YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVF05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KSTAGE2,&
  &                                       YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZF05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KSTAGE2,&
  &                                       YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWF05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KSTAGE2,&
  &                                       YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUS5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       KSTAGE3,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVS5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       KSTAGE3,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZS5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       KSTAGE3,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWS5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       KSTAGE3,YDML_DYN%YRDYN%NITMP)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3,YDML_DYN%YRDYN%NITMP, &
 & KSTAGE3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLH0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM, &
 & YDML_DYN%YRDYN%NITMP,KSTAGE3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM, &
 & YDML_DYN%YRDYN%NITMP,KSTAGE3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTSCO%NDIM, &
 & YDML_DYN%YRDYN%NITMP,KSTAGE2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLEV5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP, &
 & KSTAGE2)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POUT5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KDEP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP, &
 & KSTAGE3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KNOWENO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP, &
 & KSTAGE3)


!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZGMDTX(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZGMDTY(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)
INTEGER(KIND=JPIM) :: ILEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: ISTALAT(YDSL%NDGSAH:YDSL%NDGENH)

INTEGER(KIND=JPIM) :: IGLGLO, IHVI, IVSEPC, IVSEPL, JGL
INTEGER(KIND=JPIM) :: IROT
INTEGER(KIND=JPIM) :: ITIP
LOGICAL :: LLVSEP
LOGICAL :: LLINTV

REAL(KIND=JPRB) :: ZUF0(1), ZURL0(1), ZVF0(1), ZZRL0(1)
REAL(KIND=JPRB) :: ZVRL0(1), ZZF0(1), ZWF0(1), ZWFSM(1), ZWRL0(1)

REAL(KIND=JPRB) :: ZUN_STDDISU(1),ZUN_STDDISV(1),ZUN_STDDISW(1)

REAL(KIND=JPRB) :: ZDEPI

REAL(KIND=JPRB), POINTER :: ZSLB2KAPPA5(:), ZSLB2KAPPAT5(:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "larcina.intfb.h"
#include "elarmes5.intfb.h"
#include "larmes5.intfb.h"
#include "larmes_rk5.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAPINEA5',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDMP=>YDGEOMETRY%YRMP, YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL),YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), &
 & YDEGEO=>YDGEOMETRY%YREGEO, &
 & YDDYN=>YDML_DYN%YRDYN,YDPTRSLB2=>YDML_DYN%YRPTRSLB2, &
 & YDPTRSLB15=>YDML_DYN%YRPTRSLB15)

ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & LADVF=>YDDYN%LADVF, NITMP=>YDDYN%NITMP, VETAON=>YDDYN%VETAON, &
 & LSLDP_RK=>YDDYN%LSLDP_RK, VETAOX=>YDDYN%VETAOX, &
 & EDELX=>YDEGEO%EDELX, EDELY=>YDEGEO%EDELY, &
 & NFLDSLB15=>YDPTRSLB15%NFLDSLB15, &
 & MSLB2KAPPA5=>YDPTRSLB2%MSLB2KAPPA5, MSLB2KAPPAT5=>YDPTRSLB2%MSLB2KAPPAT5, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2, &
 & RTDT=>YDRIP%RTDT)

CALL SC2PRG(MSLB2KAPPA5,PB2      ,ZSLB2KAPPA5)
CALL SC2PRG(MSLB2KAPPAT5,PB2      ,ZSLB2KAPPAT5)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS:

ZDEPI=2.0_JPRB*RPI
DO JGL=MAX(YDSL%NDGSAG,YDSL%NDGSAL+YDSL%NFRSTLOFF-YDSL%NSLWIDE)-YDSL%NFRSTLOFF,&
   & MIN(YDSL%NDGENG,YDSL%NDGENL+YDSL%NFRSTLOFF+YDSL%NSLWIDE)-YDSL%NFRSTLOFF  
  IGLGLO=JGL+YDSL%NFRSTLOFF
  ZLSDEPI(JGL)=REAL(YDSL%NLOENG(IGLGLO),JPRB)/ZDEPI
ENDDO

IF(LRPLANE) THEN
  ZGMDTX (KST:KPROF)=0.5_JPRB*YDGSGEOM%GM(KST:KPROF)*RTDT/EDELX
  ZGMDTY (KST:KPROF)=0.5_JPRB*YDGSGEOM%GM(KST:KPROF)*RTDT/EDELY
ENDIF

DO JGL=MAX(YDSL%NDGSAH,LBOUND(YDSL%NSLOFF,1)),MIN(YDSL%NDGENH,UBOUND(YDSL%NSLOFF,1))
  ISTALAT(JGL)=YDSL%NSLOFF(JGL)
ENDDO

!     ------------------------------------------------------------------

!*       2.    COMPUTATION OF THE SL TRAJECTORY.
!              ---------------------------------

! LAM only variables. Better to protect those for the LELAM=F case.
POUT5(KST:KPROF,1:NFLEVG,1:NITMP)=0.0_JPRB  ! When RK4 implemented to LAM, this should be 4D array
KDEP(KST:KPROF,1:NFLEVG,1:NITMP,1:KSTAGE3) =0

IF (LRPLANE) THEN
  CALL ELARMES5(YDGEOMETRY,YDML_DYN,YDRIP,KST,KPROF,YDSL,ISTALAT,ZGMDTX,ZGMDTY,PB15,PB2,&
   & ZLSDEPI,KIBL,KDEP,KNOWENO,&  ! LAM people add also KDEP,KNOWENO here. They are required for LHOISLT=true
   & PSCO5(:,:,:,:,KSTAGE2),PLEV5(:,:,:,KSTAGE2),PCCO5(:,:,:,:,KSTAGE3),PUF5,PVF5,&
   & KL0(:,:,:,:,KSTAGE3),KLH0,ILEV,PLSCAW5(:,:,:,:,KSTAGE3),POUT5)  
ELSE
  IF (LSLDP_RK) THEN
    CALL LARMES_RK5(YDGEOMETRY,YDRIP,YDML_DYN,KST,KPROF,YDSL,ISTALAT,PB15,PB2,&
     & ZLSDEPI,KIBL,KDEP,KNOWENO,&
     & PSCO5,PLEV5,PCCO5,PUF05,PVF05,PZF05,PWF05,PUS5,PVS5,PZS5,PWS5,&
     & KL0,KLH0,ILEV,PLSCAW5,PRSCAW5)  
  ELSE
    CALL LARMES5(YDGEOMETRY,YDRIP,YDML_DYN,KST,KPROF,YDSL,ISTALAT,PB15,PB2,&
     & ZLSDEPI,KIBL,KDEP(:,:,:,KSTAGE3),KNOWENO(:,:,:,KSTAGE3),&
     & PSCO5(:,:,:,:,KSTAGE2),PLEV5(:,:,:,KSTAGE2),PCCO5(:,:,:,:,KSTAGE3),PUF5,PVF5,PZF5,&
     & KL0(:,:,:,:,KSTAGE3),KLH0,ILEV,PLSCAW5(:,:,:,:,KSTAGE3),PRSCAW5(:,:,:,:,KSTAGE3))
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF WEIGHTS AND INDICES ARRAYS FOR INTERPOLATIONS.
!              -------------------------------------------------------------

!*       3.1   Full-level Quantities.
!              ---------------------

IHVI=0
LLVSEP=.FALSE.

IF(LD_AD)THEN

  IF (LRPLANE.AND.LADVF) THEN
    !*Use PSCO5 (SINLA and COPHI to transfer origin point coordinates to LAPINEAAD
    ! before they are recalculated by "2 Omega vectorial a k" recalculation.
    PSCO5(KST:KPROF,1:NFLEVG,YDML_DYN%YYTSCO%M_SINLA,NITMP,KSTAGE2)= &
     &   PCCO5(KST:KPROF,1:NFLEVG,YDML_DYN%YYTCCO%M_RLON,NITMP,KSTAGE3)
    PSCO5(KST:KPROF,1:NFLEVG,YDML_DYN%YYTSCO%M_COPHI,NITMP,KSTAGE2)= &
     &   PCCO5(KST:KPROF,1:NFLEVG,YDML_DYN%YYTCCO%M_RLAT,NITMP,KSTAGE3)
  ENDIF

  IROT=1
  ITIP=3

  LLINTV=.FALSE.

  !  Warning: ZURL0,ZVRL0,ZZRL0,ZWRL0,ZUF0,ZVF0,ZZF0,ZWF0,ZWFSM are not
  !           used in this call to LARCINA.

  CALL LARCINA(YDGEOMETRY,YDML_DYN,KST,KPROF,YDSL,IHVI,ISTALAT,LLVSEP,YRDYNA%LSLHD,YRDYNA%LSLHDQUAD,LLINTV,&
   & ITIP,IROT,LRPLANE,ZLSDEPI,&
   & KIBL,&
   & PSCO5(1,1,1,NITMP,KSTAGE2),PLEV5(1,1,NITMP,KSTAGE2),&
   & ZSLB2KAPPA5,ZSLB2KAPPAT5,ZSLB2KAPPA5,ZSLB2KAPPA5,&
   & ZUN_STDDISU,ZUN_STDDISV,ZUN_STDDISW,&
   & ZURL0,ZVRL0,ZZRL0,ZWRL0,&
   & IVSEPC,IVSEPL,PCCO5(1,1,1,NITMP,KSTAGE3),&
   & ZUF0,ZVF0,ZZF0,ZWF0,ZWFSM,&
   & KL0(1,1,0,NITMP,KSTAGE3),KLH0,ILEV,PLSCAW5(1,1,1,NITMP,KSTAGE3),PRSCAW5(1,1,1,NITMP,KSTAGE3),&
   & KDEP(1,1,NITMP,KSTAGE3),KNOWENO(1,1,NITMP,KSTAGE3))

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LAPINEA5',1,ZHOOK_HANDLE)
END SUBROUTINE LAPINEA5
