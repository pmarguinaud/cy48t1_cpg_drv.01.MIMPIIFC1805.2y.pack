SUBROUTINE LAPINEAAD(YDGEOMETRY,YDRIP,YDML_DYN,KST,KPROF,YDSL,KSTGLO,&
 & KIBL,KSTAGE2,KSTAGE3,KDEP,KNOWENO,&
 & PINC,KINC,KDIM1,KDIM2,KDIM3,&
 & PB1,PB2,PB15,&
 & PCCO,PCCO5,&
 & KL0,PLSCAW,PRSCAW,PLSCAW5,PRSCAW5,&
 & PSCO,PLEV,PSCO5,PLEV5,&
 & PUF5,PVF5,PZF5,PUF05,PVF05,PZF05,PWF05,PUS5,PVS5,PZS5,PWS5,POUT5,KLOCK)

!**** *LAPINEAAD* - semi-LAgrangian scheme:(Interpolation - Part A)
!                                               (adjoint version)
!                Interface subroutine for interpolations and lagrangian
!                trends. (Programme INterface d'Ensemble).

!     Purpose.
!     --------
!           Grid point calculations in dynamics.

!           Computes the semi-Lagrangian trajectory.
!           Computes the weights ready for interpolations at origin point.

!**   Interface.
!     ----------
!        *CALL* *LAPINEAAD(.....)

!        Explicit arguments :
!        --------------------
!          KST      - first element of arrays where computations are performed.
!          KPROF    - depth of work.
!          YDSL     - SL_STRUCT definition.
!          KSTGLO   - global offset into nproma buffer.
!          KIBL     - index into YRGSGEOM/YRCSGEOM instances in YDGEOMETRY
!          KSTAGE2/KSTAGE3 - dimension for 2 and 3 stage variables (used by RK4)
!          KDEP    - dependency indicator for the southernost or northermost latitudes (LAM only)
!          KNOWENO - vertical boundary indicator for the WENO/quintic interpolation.
!          KDIM1,KDIM2,KDIM3 - dimensions of KINC, PINC
!          KINC     - addresses of interpolation increments in global buffer
!          PINC     - interpolation increments to global buffer (vector only)
!          PB1      - SLBUF1-buffer for interpolations.
!          PB2      - SLBUF2-buf to communicate info from non lag. to lag. dyn.
!          PB15     - SLBUF1-buf traj for interpolations.

!          PCCO     - information about comput. space position of interpol. point.
!          PCCO5    - cf. PCCO (trajectory)
!          KL0      - index of the four western points
!                     of the 16 points interpolation grid.
!          PLSCAW   - linear weights (distances) for interpolations.
!          PRSCAW   - non-linear weights for interpolations.
!          PLSCAW5  - linear weights (distances) for interpolations (trajectory).
!          PRSCAW5  - non-linear weights for interpolations (trajectory).
!          PSCO     - information about geographic position of interpol. point.
!          PLEV     - vertical coordinate of the interpolation point.
!          PSCO5    - cf. PSCO (trajectory)
!          PLEV5    - cf. PLEV (trajectory)
!          PUF5     - interpolated U-wind (trajectory)
!          PVF5     - interpolated V-wind (trajectory)
!          PZF5     - interpolated Z-wind (trajectory)
!          PUF05    - U-comp of "(a/rs)*wind" interpolated quantity for stages 1 & 2
!          PVF05    - V-comp of "(a/rs)*wind" interpolated quantity for stages 1 & 2
!          PZF05    - Z-comp of "(a/rs)*wind" interpolated quantity for stages 1 & 2
!          PWF05    - Z-comp of "(a/rs)*wind" interpolated quantity for stages 1 & 2
!          PUS5,PVS5,PZS5,PWS5 - stage wind in Cartesian space
!          POUT5    - indicator of O point of the SL trajectory (for LAM only).
!                     If O lays inside the domain, the value is zero.
!          KLOCK    - Semaphor counter (used for OpenMP synchronization).

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!       See includes below.
!       Called by CALL_SL_AD.

!     Reference.
!     ----------
!        ARPEGE documentation about semi-Lagrangian scheme.

!     Author.
!     -------
!      C. Temperton (ECMWF)
!      Original : 99/10/15

! Modifications (list from LAPINEAD)
! -------------
!   Modified Jul 2001 by K. YESSAD: update dummy arguments of LAINOR.
!   Modified 01-08-30 by K. YESSAD: pruning and some other cleanings.
!   R. El Khatib : 01-08-07 Pruning options
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!   16-Jun-2004 J. Masek   NH cleaning (LPC_NOTR)
!   Modified 05-01-11 by C. Temperton: streamlining.
!   Modified 05-02-02 by C. Temperton: transfer computations from
!                            end of LAIDEPAD to start of LARMESAD.
!   D.Salmond     08-Feb-2005 New LAPINEAAD
!   K. Yessad     07-Sep-2005 Now calls directly LARCINAAD
!   Modified 21-06-2006 by F. Vana: LRPLANE=.T.
!   F. Vana       08-Jan-2007 KDEP for vectorized (LAM) adjoint
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   F. Vana       14-Aug-2008 Weights driven interpolation + optimization.
!   K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!   F. Vana :     11-Dec-2008 OpenMP for vector platforms
!   F. Vana     13-Jan-2009: TL/AD of SLHD
!   K. Yessad (Aug 2009): always use root (QX,QY) for (p,q) variables names
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!   G. Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
!   G. Mozdzynski (May 2012): further cleaning
!   F. Vana  13-Feb-2014 KappaT for heat quentities.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   F. Vana    21-Nov-2017: Options LSLDP_CURV and LHOISLT
!   F. Vana July 2018: RK4 scheme for trajectory research.
!   F. Vana October 2018: Extended LSLDP_CURV.
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
!   F. Vana 18-Mar-2019: Implicit RK4 following Lobatto IIIA 
! End Modifications
!------------------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD     , ONLY : MODEL_DYNAMICS_TYPE
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK
USE SC2PRG_MOD             , ONLY : SC2PRG
USE YOMCT0                 , ONLY : LRPLANE
USE YOMDYNA                , ONLY : YRDYNA
USE YOMCST                 , ONLY : RPI
USE YOMRIP                 , ONLY : TRIP
USE EINT_MOD               , ONLY : SL_STRUCT

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTGLO
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTAGE2
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTAGE3
INTEGER(KIND=JPIM),INTENT(IN)    :: KDEP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP,KSTAGE3)
INTEGER(KIND=JPIM),INTENT(IN)    :: KNOWENO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP,KSTAGE3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KINC(KDIM1,KDIM2,KDIM3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PINC(KDIM1,KDIM2,KDIM3)
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIM1
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIM2
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIM3
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1(YDSL%NASLB1,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB15(YDSL%NASLB1,YDML_DYN%YRPTRSLB15%NFLDSLB15)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM, &
 & YDML_DYN%YRDYN%NITMP,KSTAGE3)
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3,YDML_DYN%YRDYN%NITMP,KSTAGE3)
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUF05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KSTAGE2,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVF05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KSTAGE2,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZF05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KSTAGE2,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWF05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KSTAGE2,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUS5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KSTAGE3,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVS5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KSTAGE3,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZS5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KSTAGE3,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWS5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KSTAGE3,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POUT5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KLOCK(:)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)

REAL(KIND=JPRB) :: ZGMDTX(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZGMDTY(YDGEOMETRY%YRDIM%NPROMA)

INTEGER(KIND=JPIM) :: ISTALAT(YDSL%NDGSAH:YDSL%NDGENH)
INTEGER(KIND=JPIM) :: IGLGLO, JGL
INTEGER(KIND=JPIM) :: ITIP, IDIM, IROT, IL0(1,1,0:3)

LOGICAL :: LLOCK_UP

REAL(KIND=JPRB) :: ZDEPI

INTEGER(KIND=JPIM) :: IDEP(1),IINC(1)
REAL(KIND=JPRB) :: ZUF0(1),ZVF0(1),ZZF0(1),ZWF0(1),ZZRL0(1),ZZRL05(1),ZWRL0(1),ZWRL05(1),ZINC(1)
REAL(KIND=JPRB) :: ZLSCAW5(1,1,YDML_DYN%YYTLSCAW%NDIM)
REAL(KIND=JPRB) :: ZRSCAW5(1,1,YDML_DYN%YYTLSCAW%NDIM)

REAL(KIND=JPRB), POINTER :: ZSLB1UR9(:), ZSLB1VR9(:), ZSLB2KAPPA(:), ZSLB2KAPPA5(:)
REAL(KIND=JPRB), POINTER :: ZSLB2KAPPAT(:), ZSLB2KAPPAT5(:)
REAL(KIND=JPRB), POINTER :: ZSLB1UR95(:), ZSLB1VR95(:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "larcinaad.intfb.h"
#include "larmesad.intfb.h"
#include "larmes_rk_ad.intfb.h"
#include "elarmesad.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAPINEAAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 &  YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL), YDEGEO=>YDGEOMETRY%YREGEO, YDDYN=>YDML_DYN%YRDYN, &
 & YDPTRSLB1=>YDML_DYN%YRPTRSLB1, YDPTRSLB2=>YDML_DYN%YRPTRSLB2, &
 & YDPTRSLB15=>YDML_DYN%YRPTRSLB15)

ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & LADVF=>YDDYN%LADVF, NITMP=>YDDYN%NITMP, NVSEPC=>YDDYN%NVSEPC, &
 & NVSEPL=>YDDYN%NVSEPL, VETAON=>YDDYN%VETAON, VETAOX=>YDDYN%VETAOX, &
 & LSLDP_RK=>YDDYN%LSLDP_RK, EDELX=>YDEGEO%EDELX, EDELY=>YDEGEO%EDELY, &
 & MSLB1UR9=>YDPTRSLB1%MSLB1UR9, MSLB1VR9=>YDPTRSLB1%MSLB1VR9, &
 & NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & MSLB1UR95=>YDPTRSLB15%MSLB1UR95, MSLB1VR95=>YDPTRSLB15%MSLB1VR95, &
 & NFLDSLB15=>YDPTRSLB15%NFLDSLB15, &
 & MSLB2KAPPA=>YDPTRSLB2%MSLB2KAPPA, MSLB2KAPPA5=>YDPTRSLB2%MSLB2KAPPA5, &
 & MSLB2KAPPAT=>YDPTRSLB2%MSLB2KAPPAT, MSLB2KAPPAT5=>YDPTRSLB2%MSLB2KAPPAT5, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2, &
 & RTDT=>YDRIP%RTDT)

CALL SC2PRG(MSLB1UR95 ,PB15     ,ZSLB1UR95)
CALL SC2PRG(MSLB1VR95 ,PB15     ,ZSLB1VR95)
CALL SC2PRG(MSLB1UR9  ,PB1      ,ZSLB1UR9)
CALL SC2PRG(MSLB1VR9  ,PB1      ,ZSLB1VR9)
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

!*       3.    COMPUTATION OF WEIGHTS AND INDICES ARRAYS FOR INTERPOLATIONS.
!              -------------------------------------------------------------

IF (LRPLANE) THEN
! * Use PSCO (SINLA and COPHI) to transfer origin point coordinates from LAPINEA5
!   before they are overwritten by "2 Omega vectorial a k". The perturbations
!   are stored there in LARCINB in order to be consistent with the NL code.
  IF(LADVF) THEN
    PCCO5(KST:KPROF,1:NFLEVG,YDML_DYN%YYTCCO%M_RLON,NITMP,KSTAGE3)= &
     & PSCO5(KST:KPROF,1:NFLEVG,YDML_DYN%YYTSCO%M_SINLA,NITMP,KSTAGE2)
    PCCO5(KST:KPROF,1:NFLEVG,YDML_DYN%YYTCCO%M_RLAT,NITMP,KSTAGE3)= &
     & PSCO5(KST:KPROF,1:NFLEVG,YDML_DYN%YYTSCO%M_COPHI,NITMP,KSTAGE2)
    PCCO(KST:KPROF,1:NFLEVG,YDML_DYN%YYTCCO%M_RLON)=PSCO(KST:KPROF,1:NFLEVG,YDML_DYN%YYTSCO%M_SINLA)
    PCCO(KST:KPROF,1:NFLEVG,YDML_DYN%YYTCCO%M_RLAT)=PSCO(KST:KPROF,1:NFLEVG,YDML_DYN%YYTSCO%M_COPHI)
    PSCO(KST:KPROF,1:NFLEVG,YDML_DYN%YYTSCO%M_SINLA)=0.0_JPRB
    PSCO(KST:KPROF,1:NFLEVG,YDML_DYN%YYTSCO%M_COPHI)=0.0_JPRB
  ENDIF
ENDIF

IROT=1
IDIM=1
LLOCK_UP=.FALSE.
ITIP=3

! * Remark: IL0,ZLSCAW5,ZRSCAW5,ZZRL0,ZZRL05,ZWRL0,ZWRL05,ZUF0,ZVF0,ZZF0,ZWF0,IDEP,ZINC,IINC
!   are not used in this call of LARCINAAD.

CALL LARCINAAD(YDGEOMETRY,YDML_DYN,IDIM,KSTGLO,KST,KPROF,YDSL,&
 & ISTALAT,YRDYNA%LSLHD,YRDYNA%LSLHDQUAD,YRDYNA%LSLHD_OLD,ITIP,IROT,LRPLANE,LLOCK_UP,ZLSDEPI,&
 & KIBL,IDEP,KNOWENO,&
 & PSCO,PLEV,&
 & PSCO5(1,1,1,NITMP,KSTAGE2),PLEV5(1,1,NITMP,KSTAGE2),&
 & ZSLB2KAPPA,ZSLB2KAPPA5,&
 & ZSLB2KAPPAT,ZSLB2KAPPAT5,&
 & ZSLB1UR9,ZSLB1VR9,ZZRL0,ZWRL0,&
 & ZSLB1UR95,ZSLB1VR95,ZZRL05,ZWRL05,&
 & NVSEPC,NVSEPL,KLOCK,&
 & ZINC,IINC,1,1,&
 & PCCO,ZUF0,ZVF0,ZZF0,ZWF0,&
 & PCCO5(1,1,1,NITMP,KSTAGE3),&
 & IL0,PLSCAW,PRSCAW,ZLSCAW5,ZRSCAW5)

!     ------------------------------------------------------------------

!*       2.    COMPUTATION OF THE SL TRAJECTORY.
!              ---------------------------------

IF (LRPLANE) THEN
  CALL ELARMESAD(YDGEOMETRY,YDRIP,YDML_DYN,KSTGLO,KST,KPROF,YDSL,ISTALAT,ZGMDTX,ZGMDTY,PB1,PB15,PB2,&
   & PINC,KINC,KDIM1,KDIM2,&
   & ZLSDEPI,KIBL,KDEP,KNOWENO,&
   & PSCO,PLEV,PCCO,&
   & KL0(:,:,:,:,KSTAGE3),KLOCK,PSCO5(:,:,:,:,KSTAGE2),PLEV5(:,:,:,KSTAGE2),&
   & PCCO5(:,:,:,:,KSTAGE3),PLSCAW5(:,:,:,:,KSTAGE3),PRSCAW5(:,:,:,:,KSTAGE3),POUT5)
ELSE
  IF (LSLDP_RK) THEN
    CALL LARMES_RK_AD(YDGEOMETRY,YDRIP,YDML_DYN,KSTGLO,KST,KPROF,YDSL,ISTALAT,PB1,PB15,PB2,&
     & PINC,KINC,KDIM1,KDIM2,&
     & ZLSDEPI,KIBL,KDEP,KNOWENO,&
     & PSCO,PLEV,PCCO,&
     & KL0,KLOCK,&
     & PSCO5,PLEV5,PCCO5,PUF05,PVF05,PZF05,PZF05,PUS5,PVS5,PZS5,PWS5,PLSCAW5,PRSCAW5)  
  ELSE
    CALL LARMESAD(YDGEOMETRY,YDRIP,YDML_DYN,KSTGLO,KST,KPROF,YDSL,ISTALAT,PB1,PB15,PB2,&
     & PINC,KINC,KDIM1,KDIM2,&
     & ZLSDEPI,KIBL,KDEP,KNOWENO,&
     & PSCO,PLEV,PCCO,&
     & KL0(:,:,:,:,KSTAGE3),KLOCK,PSCO5(:,:,:,:,KSTAGE2),PLEV5(:,:,:,KSTAGE2),&
     & PCCO5(:,:,:,:,KSTAGE3),PUF5,PVF5,PZF5,PLSCAW5(:,:,:,:,KSTAGE3),PRSCAW5(:,:,:,:,KSTAGE3))
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LAPINEAAD',1,ZHOOK_HANDLE)
END SUBROUTINE LAPINEAAD
