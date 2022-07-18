SUBROUTINE LAPINEATL(YDGEOMETRY,YDRIP,YDML_DYN,KST,KPROF,YDSL,&
 & KIBL,KSTAGE2,KSTAGE3,&
 & PB1,PB2,PB15,&
 & PCCO,PCCO5,&
 & KL0,KLH0,PLSCAW,PRSCAW,PLSCAW5,PRSCAW5,&
 & PSCO,PLEV,PSCO5,PLEV5,&
 & PUF5,PVF5,PZF5,PUS5,PVS5,PZS5,POUT5,KNOWENO)

!**** *LAPINEATL* - semi-LAgrangian scheme:(Trajectory)
!                                         (tangent-linear version)
!                Interface subroutine for interpolations and lagrangian
!                trends. (Programme INterface d'Ensemble).

!     Purpose.
!     --------
!           Grid point calculations in dynamics.

!           Computes the semi-Lagrangian trajectory.
!           Computes the weights ready for interpolations at origin point.

!**   Interface.
!     ----------
!        *CALL* *LAPINEATL(.....)

!        Explicit arguments :
!        --------------------
!        INPUT:
!          KST      - first element of arrays where computations are performed.
!          KPROF    - depth of work.
!          YDSL     - SL_STRUCT definition.
!          KIBL     - index into YRCSGEOM/YRGSGEOM instances in YDGEOMETRY
!          KSTAGE2/KSTAGE3 - dimension for 2 and 3 stage variables (used by RK4)
!          PB1      - SLBUF1-buffer for interpolations.
!          PB2      - SLBUF2-buf to communicate info from non lag. to lag. dyn.
!          PB15     - SLBUF1-buf traj for interpolations.

!        OUTPUT AND TRAJECTORY-OUTPUT:
!          PCCO     - information about comput. space position of interpol. point.
!          PCCO5    - cf. PCCO (trajectory)
!          KL0      - index of the four western points
!                     of the 16 points interpolation grid.
!          KLH0     - second value of index of the four western points
!                     of the 16 points interpolation grid if needed.
!          PLSCAW   - linear weights (distances) for interpolations.
!          PRSCAW   - non-linear weights for interpolations.
!          PLSCAW5  - linear weights (distances) for interpolations (trajectory).
!          PRSCAW5  - non-linear weights for interpolations (trajectory).
!          PSCO     - information about geographic position of interpol. point.
!          PLEV     - vertical coordinate of the interpolation point.
!          PSCO5    - cf. PSCO (trajectory)
!          PLEV5    - cf. PLEV (trajectory)
!          KNOWENO  - indication of the specific treatment of vertical boundary used
!                       by WENO scheme 

!        TRAJECTORY-INPUT:
!          PUF5     - interpolated U-wind.                             (in)
!          PVF5     - interpolated V-wind.                             (in)
!          PZF5     - interpolated Z-wind.                             (in)
!          PUS5     - stage U-wind.                                    (in)
!          PVS5     - stage V-wind.                                    (in)
!          PZS5     - stage Z-wind.                                    (in)
!          POUT5    - indicator of O point of the SL trajectory.       (in)
!                     If O lays inside the domain, the value is zero.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!       See includes below.
!       Called by CALL_SL_TL.

!     Reference.
!     ----------
!        ARPEGE documentation about semi-Lagrangian scheme.

!     Author.
!     -------
!      C. Temperton (ECMWF)
!      Original : 99/07/13

!     Modifications.
!     --------------
!      Modified 01-08-30 by K. YESSAD: pruning and some other cleanings.
!      R. El Khatib : 01-08-07 Pruning options
!      Modified 03-04-09 by C. Temperton: use SLADJTR to repeat    
!                                         trajectory calculations.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      D.Salmond     22-Dec-2003 Karim's Larcin cleaning
!      Modified 05-01-11 by C. Temperton: streamlining.
!      Modified 05-01-14 by C. Temperton: transfer computations from
!                              start of LAIDEP to end of LARMES.
!      D.Salmond     08-Feb-2005 SLADJTR replaced by LARMES5
!      K. Yessad     07-Sep-2005 Now calls directly LARCINATL.
!      F. Vana       01-Jun-2006 LRPLANE=.T.
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      F. Vana       27-Jun-2008 Weights driven interpolation
!      K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!      F. Vana 13-Jan-2009: TL/AD of SLHD
!      K. Yessad (Aug 2009): always use root (QX,QY) for (p,q) variables names
!      K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!      G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      G. Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
!      G. Mozdzynski (May 2012): further cleaning
!      F. Vana  13-Feb-2014  KappaT for heat quantities.
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      F. Vana    21-Nov-2017: Options LSLDP_CURV and LHOISLT
!      F. Vana    July 2018: RK4 scheme for trajectory research.
!      F. Vana October 2018: Extended LSLDP_CURV.
!      F. Vana    26-Feb-2019: Vertical quintic interpolation.
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
!      F. Vana    18-Mar-2019: Implicit RK4 for trajectory research.
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
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTAGE2
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTAGE3
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB1(YDSL%NASLB1,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB15(YDSL%NASLB1,YDML_DYN%YRPTRSLB15%NFLDSLB15)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM, &
 & YDML_DYN%YRDYN%NITMP,KSTAGE3)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3,YDML_DYN%YRDYN%NITMP,KSTAGE3)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KLH0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM, &
 & YDML_DYN%YRDYN%NITMP,KSTAGE3)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM, &
 & YDML_DYN%YRDYN%NITMP,KSTAGE3)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTSCO%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTSCO%NDIM, &
 & YDML_DYN%YRDYN%NITMP,KSTAGE2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLEV5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP,KSTAGE2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP-1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP-1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP-1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUS5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KSTAGE3,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVS5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KSTAGE3,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZS5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KSTAGE3,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POUT5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KNOWENO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP,KSTAGE3)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZGMDTX(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZGMDTY(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)

INTEGER(KIND=JPIM) :: ISTALAT(YDSL%NDGSAH:YDSL%NDGENH)
INTEGER(KIND=JPIM) :: ILEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IGLGLO, JGL
INTEGER(KIND=JPIM) :: ITIP
INTEGER(KIND=JPIM) :: IROT

REAL(KIND=JPRB) :: ZUF0(1), ZVF0(1), ZZF0(1), ZWF0(1), ZWRL0(1),&
 & ZWRL05(1), ZURL(1),ZURL5(1),ZVRL(1),ZVRL5(1), ZZRL(1), ZZRL5(1)

REAL(KIND=JPRB) :: ZDEPI

REAL(KIND=JPRB), POINTER :: ZSLB2KAPPA(:), ZSLB2KAPPA5(:), ZSLB2KAPPAT(:), ZSLB2KAPPAT5(:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "larcinatl.intfb.h"
#include "elarmestl.intfb.h"
#include "larmestl.intfb.h"
#include "larmes_rk_tl.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAPINEATL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL), YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDEGEO=>YDGEOMETRY%YREGEO, &
 & YDDYN=>YDML_DYN%YRDYN, &
 & YDPTRSLB1=>YDML_DYN%YRPTRSLB1,YDPTRSLB2=>YDML_DYN%YRPTRSLB2, &
 & YDPTRSLB15=>YDML_DYN%YRPTRSLB15)

ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & LADVF=>YDDYN%LADVF, NITMP=>YDDYN%NITMP, VETAON=>YDDYN%VETAON, &
 & LSLDP_RK=>YDDYN%LSLDP_RK, VETAOX=>YDDYN%VETAOX, &
 & EDELX=>YDEGEO%EDELX, EDELY=>YDEGEO%EDELY, &
 & NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & NFLDSLB15=>YDPTRSLB15%NFLDSLB15, &
 & MSLB2KAPPA=>YDPTRSLB2%MSLB2KAPPA, MSLB2KAPPA5=>YDPTRSLB2%MSLB2KAPPA5, &
 & MSLB2KAPPAT=>YDPTRSLB2%MSLB2KAPPAT, MSLB2KAPPAT5=>YDPTRSLB2%MSLB2KAPPAT5, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2, &
 & RTDT=>YDRIP%RTDT)

CALL SC2PRG(MSLB2KAPPA,PB2      ,ZSLB2KAPPA)
CALL SC2PRG(MSLB2KAPPA5,PB2     ,ZSLB2KAPPA5)
CALL SC2PRG(MSLB2KAPPAT,PB2     ,ZSLB2KAPPAT)
CALL SC2PRG(MSLB2KAPPAT5,PB2    ,ZSLB2KAPPAT5)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS:
!              ----------------------------

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

IF (LRPLANE) THEN
  CALL ELARMESTL(YDGEOMETRY,YDML_DYN,YDRIP,KST,KPROF,YDSL,ISTALAT,ZGMDTX,ZGMDTY,PB1,PB15,PB2,&
   & ZLSDEPI,KIBL,&
   & PSCO,PLEV,PCCO,&
   & KL0(:,:,:,:,KSTAGE3),KLH0,ILEV,PLSCAW,PRSCAW,&
   & PSCO5(:,:,:,:,KSTAGE2),PLEV5(:,:,:,KSTAGE2),PCCO5(:,:,:,:,KSTAGE3),&
   & PLSCAW5(:,:,:,:,KSTAGE3),PRSCAW5(:,:,:,:,KSTAGE3),POUT5)
ELSE
  IF (LSLDP_RK) THEN
    CALL LARMES_RK_TL(YDGEOMETRY,YDRIP,YDML_DYN,KST,KPROF,YDSL,ISTALAT,PB1,PB15,PB2,&
     & ZLSDEPI,KIBL,&
     & PSCO,PLEV,PCCO,&
     & KL0,KLH0,ILEV,PLSCAW,PRSCAW,&
     & PSCO5,PLEV5,PCCO5,PUS5,PVS5,PZS5,PLSCAW5,PRSCAW5)
  ELSE
    CALL LARMESTL(YDGEOMETRY,YDRIP,YDML_DYN,KST,KPROF,YDSL,ISTALAT,PB1,PB15,PB2,&
     & ZLSDEPI,KIBL,&
     & PSCO,PLEV,PCCO,&
     & KL0(:,:,:,:,KSTAGE3),KLH0,ILEV,PLSCAW,PRSCAW,&
     & PSCO5(:,:,:,:,KSTAGE2),PLEV5(:,:,:,KSTAGE2),PCCO5(:,:,:,:,KSTAGE3),&
     & PUF5,PVF5,PZF5,PLSCAW5(:,:,:,:,KSTAGE3),PRSCAW5(:,:,:,:,KSTAGE3))
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF WEIGHTS AND INDICES ARRAYS FOR INTERPOLATIONS.
!              -------------------------------------------------------------

ITIP=3
IROT=1

! Remark: ZURL,ZVRL,ZZRL,ZWRL0,ZURL5,ZVRL5,ZZRL5,ZWRL05,ZUF0,ZVF0,ZZF0,ZWF0
! are not used in this call of LARCINATL.

CALL LARCINATL(YDGEOMETRY,YDML_DYN,KST,KPROF,YDSL,ISTALAT,YRDYNA%LSLHD,YRDYNA%LSLHDQUAD,YRDYNA%LSLHD_OLD,ITIP,IROT,LRPLANE,&
 & ZLSDEPI,KIBL,&
 & PSCO,PLEV,&
 & PSCO5(1,1,1,NITMP,KSTAGE2),PLEV5(1,1,NITMP,KSTAGE2),&
 & ZSLB2KAPPA,ZSLB2KAPPA5,&
 & ZSLB2KAPPAT,ZSLB2KAPPAT5,&
 & ZURL,ZVRL,ZZRL,ZWRL0,ZURL5,ZVRL5,ZZRL5,ZWRL05,&
 & PCCO,&
 & ZUF0,ZVF0,ZZF0,ZWF0,&
 & PCCO5(1,1,1,NITMP,KSTAGE3),&
 & KL0(1,1,0,NITMP,KSTAGE3),KLH0,ILEV,PLSCAW,PRSCAW,&
 & PLSCAW5(1,1,1,NITMP,KSTAGE3),PRSCAW5(1,1,1,NITMP,KSTAGE3),KNOWENO(1,1,NITMP,KSTAGE3))

IF (LRPLANE) THEN
! * Use PSCO (SINLA and COPHI) to transfer origin point coordinates to LAPINEB
!   for "2 Omega vectorial a k" recalculation.
  IF(LADVF) THEN
    PSCO5(KST:KPROF,1:NFLEVG,YDML_DYN%YYTSCO%M_SINLA,NITMP,KSTAGE2)= &
     &   PCCO5(KST:KPROF,1:NFLEVG,YDML_DYN%YYTCCO%M_RLON,NITMP,KSTAGE3)
    PSCO5(KST:KPROF,1:NFLEVG,YDML_DYN%YYTSCO%M_COPHI,NITMP,KSTAGE2)= &
     &   PCCO5(KST:KPROF,1:NFLEVG,YDML_DYN%YYTCCO%M_RLAT,NITMP,KSTAGE3)
    PSCO(KST:KPROF,1:NFLEVG,YDML_DYN%YYTSCO%M_SINLA)=PCCO(KST:KPROF,1:NFLEVG,YDML_DYN%YYTCCO%M_RLON)
    PSCO(KST:KPROF,1:NFLEVG,YDML_DYN%YYTSCO%M_COPHI)=PCCO(KST:KPROF,1:NFLEVG,YDML_DYN%YYTCCO%M_RLAT)
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LAPINEATL',1,ZHOOK_HANDLE)
END SUBROUTINE LAPINEATL
