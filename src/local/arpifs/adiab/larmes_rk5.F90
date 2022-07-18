!option! -O nomove
SUBROUTINE LARMES_RK5(YDGEOMETRY,YDRIP,YDML_DYN,KST,KPROF,YDSL,KSTABUF,PB15,PB2,&
 & PLSDEPI,KIBL,KDEP,KNOWENO,&
 & PSCO5,PLEV5,PCCO5,PUF05,PVF05,PZF05,PWF05,PUS5,PVS5,PZS5,PWS5,&
 & KL0,KLH0,KLEV,PLSCAW5,PRSCAW5)

!---- should be compiled for Cray with -hcontiguous----

!**** *LARMES_RK5 - semi-LAgrangian scheme:
!                Research of the origin point on the Sphere by 4th order Runge-Kutta method.

!     Purpose.
!     --------

!      The computation of the location of the interpolation point of
!     the lagrangian trajectory is performed by an 4th order Runge-Kutta 
!     method. 
!     Note the final trajectory is composite of four different trajectories on great circle.
!     Finally we find the departure (origin) point "O".

!**   Interface.
!     ----------
!        *CALL* *LARMES_RK5(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST      - first element of arrays where computations are performed.
!          KPROF    - depth of work.
!          YDSL     - SL_STRUCT definition
!          KSTABUF  - for a latitude IGL, KSTABUF(IGL) is the
!                     address of the element corresponding to
!                     (ILON=1,IGL) in the NPROMA arrays.
!          PB15     - SLBUF1-buffer for interpolations (trajectory).
!          PB2      - SLBUF2-buf to communicate info from non lag. to lag. dyn.
!          PLSDEPI  - (Number of points by latitude) / (2 * PI) .
!          KIBL     - index into YRCSGEOM/YRGSGEOM instances in YDGEOMETRY





!        OUTPUT:
!          KDEP     - indicate dependency (used in adjoint of LAM for vectorization)
!          KNOWENO  - lower boundary indicator used in vertical WENO interpolation.
!          PSCO5    - information about geographic position of interpol. point. (traj)
!          PLEV5    - vertical coordinate of the interpolation point. (traj)
!          PCCO5    - information about comput. space position of interpol. point. (traj)
!          PUF05    - U-comp of "(a/rs)*wind" interpolated quantity for stages 1 & 2
!          PVF05    - V-comp of "(a/rs)*wind" interpolated quantity for stages 1 & 2
!          PZF05    - Z-comp of "(a/rs)*wind" interpolated quantity for stages 1 & 2
!          PWF05    - Z-comp of "(a/rs)*wind" interpolated quantity for stages 1 & 2
!          PUS5,PVS5,PZS5,PWS5 - stage wind in Cartesian space 
!          KL0      - index of the four western points
!                     of the 16 points interpolation grid.
!          KLH0     - second value of index of the four western points
!                     of the 16 points interpolation grid if needed.
!          KLEV     - lower level of the vertical interpolation
!                     grid needed for vertical interpolations.
!          PLSCAW5  - linear weights (distances) for interpolations.
!          PRSCAW5  - non-linear weights for interpolations.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------
!        Calls  LARCINA
!        Is called by LAPINEA (3D model)

!     Reference.
!     ----------

!     Author.
!     -------
!      F. VANA  after LARMES5 and LARMES_RK
!      Original : MAY 2018

!     Modifications.
!     --------------
!   F. Vana October 2018: Extended LSLDP_CURV.
!   F. Vana  26-Feb-2019: Vertical quintic interpolation.
!   F. Vana  05-Mar-2019: Make the scheme implicit following Lobatto IIIA type

! End Modifications
!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMCT0    , ONLY : LTWOTL

USE YOMRIP    , ONLY : TRIP
USE EINT_MOD  , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
TYPE(SL_STRUCT)   ,INTENT(INOUT) :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTABUF(YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB15(YDSL%NASLB1,YDML_DYN%YRPTRSLB15%NFLDSLB15)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       YDML_DYN%YYTSCO%NDIM,YDML_DYN%YRDYN%NITMP,2:3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLEV5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       YDML_DYN%YRDYN%NITMP,2:3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       YDML_DYN%YYTCCO%NDIM,YDML_DYN%YRDYN%NITMP,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUF05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       2,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVF05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       2,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZF05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       2,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWF05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       2,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUS5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       3,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVS5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       3,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZS5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       3,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWS5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       3,YDML_DYN%YRDYN%NITMP)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3, &
  &                                       YDML_DYN%YRDYN%NITMP,3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLH0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       YDML_DYN%YYTLSCAW%NDIM,YDML_DYN%YRDYN%NITMP,3)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       YDML_DYN%YYTRSCAW%NDIM,YDML_DYN%YRDYN%NITMP,3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KDEP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       YDML_DYN%YRDYN%NITMP,3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KNOWENO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG, &
  &                                       YDML_DYN%YRDYN%NITMP,3)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IHVI
INTEGER(KIND=JPIM) :: IROT
INTEGER(KIND=JPIM) :: ITIP
INTEGER(KIND=JPIM) :: JITER
INTEGER(KIND=JPIM) :: JLEV
INTEGER(KIND=JPIM) :: JROF

LOGICAL :: LLINTV
LOGICAL :: LLO
LOGICAL :: LLSLHD, LLSLHDQUAD

REAL(KIND=JPRB) :: ZDTS22, ZDTS62, ZDTSA, ZINT, ZLEVO, ZNOR2, ZSPHSV

REAL(KIND=JPRB) :: ZVETAON, ZVETAOX

REAL(KIND=JPRB) :: ZPU, ZPV
REAL(KIND=JPRB) :: ZPU0, ZPV0, ZPZ0, ZPW0

REAL(KIND=JPRB) :: ZU0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZV0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZZ0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! unused arguments in call to LARCINA:
! a) input arrays - will not be used since LSLHD was set to .FALSE.
!    in LAPINEA before call to LARMES_RK5/ELARMES_RK5 => dimensions can be
!    contracted
REAL(KIND=JPRB) :: ZUN_KAPPA(1),ZUN_KAPPAT(1),ZUN_KAPPAM(1),ZUN_KAPPAH(1)
REAL(KIND=JPRB) :: ZUN_STDDISU(1),ZUN_STDDISV(1),ZUN_STDDISW(1)
! b) input/output
INTEGER(KIND=JPIM) :: IUN_VSEPC
INTEGER(KIND=JPIM) :: IUN_VSEPL
! c) output arrays - dimensions kept for safety
REAL(KIND=JPRB) :: ZUN_WFSM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "larcina.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LARMES_RK5',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDMP=>YDGEOMETRY%YRMP, YDSTA=>YDGEOMETRY%YRSTA, YDVAB=>YDGEOMETRY%YRVAB, &
 & YDVETA=>YDGEOMETRY%YRVETA, YDLAP=>YDGEOMETRY%YRLAP, YDCSGLEG=>YDGEOMETRY%YRCSGLEG,  &
 & YDVSPLIP=>YDGEOMETRY%YRVSPLIP, YDVSLETA=>YDGEOMETRY%YRVSLETA, &
 & YDHSLMER=>YDGEOMETRY%YRHSLMER, YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL), &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL),  &
 & YDSPGEOM=>YDGEOMETRY%YSPGEOM,YDDYN=>YDML_DYN%YRDYN,YDPTRSLB15=>YDML_DYN%YRPTRSLB15,  &
 & YDPTRSLB2=>YDML_DYN%YRPTRSLB2)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & LFINDVSEP=>YDDYN%LFINDVSEP, LSVTSM=>YDDYN%LSVTSM, NITMP=>YDDYN%NITMP, &
 & RPRES_SVTSM=>YDDYN%RPRES_SVTSM, VETAON=>YDDYN%VETAON, VETAOX=>YDDYN%VETAOX, &
 & LSLDP_CURV=>YDDYN%LSLDP_CURV, &
 & MSLB1UR05=>YDPTRSLB15%MSLB1UR05, MSLB1VR05=>YDPTRSLB15%MSLB1VR05, &
 & MSLB1ZR05=>YDPTRSLB15%MSLB1ZR05, MSLB1WR05=>YDPTRSLB15%MSLB1WR05, &
 & MSLB1UR005=>YDPTRSLB15%MSLB1UR005, MSLB1VR005=>YDPTRSLB15%MSLB1VR005, &
 & MSLB1ZR005=>YDPTRSLB15%MSLB1ZR005, MSLB1WR005=>YDPTRSLB15%MSLB1WR005, &
 & NFLDSLB15=>YDPTRSLB15%NFLDSLB15, &
 & MSLB2URL5=>YDPTRSLB2%MSLB2URL5, &
 & MSLB2VRL5=>YDPTRSLB2%MSLB2VRL5, MSLB2WRL5=>YDPTRSLB2%MSLB2WRL5, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2, &
 & RDTS22=>YDRIP%RDTS22, RDTS62=>YDRIP%RDTS62, RDTSA=>YDRIP%RDTSA, &
 & RTDT=>YDRIP%RTDT)

!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS AND TESTS.
!              --------------------------------------

!*       1.2   Miscellaneous preliminary initialisations.

! LLO.OR.LELTRA should now be always T for SL2TL.
IF (LTWOTL) THEN
  LLO=.NOT.YDML_DYN%YRDYNA%LELTRA
ELSE
  LLO=.FALSE.
ENDIF

IF(LLO.OR.YDML_DYN%YRDYNA%LELTRA) THEN
  ZDTS22=RDTS22*4._JPRB
  ZDTSA=RDTSA*2._JPRB
  ZDTS62=RDTS62*4._JPRB
ELSE
  ! this case may occur only for SL3TL scheme.
  IF (LTWOTL) CALL ABOR1(' LARMES_RK5 1.2 ')
  ZDTS22=RDTS22
  ZDTSA=RDTSA
  ZDTS62=RDTS62
ENDIF

! deactivate computation of SLHD weights
LLSLHD     =.FALSE.
LLSLHDQUAD =.FALSE.

!*       1.3   Computation of weights for smooth interpolation at arrival point.

IF(LSVTSM .OR. YDML_DYN%YRDYNA%LRALTVDISP .OR. YDML_DYN%YRDYNA%LELTRA) THEN
  CALL ABOR1('NOT CODED')
ENDIF

!     ------------------------------------------------------------------

!*       2.    ITERATIONS.
!              -----------


ZVETAON=(1.0_JPRB-VETAON)*YDVETA%VETAH(0)+VETAON*YDVETA%VETAF(1)
ZVETAOX=(1.0_JPRB-VETAOX)*YDVETA%VETAH(NFLEVG)+VETAOX*YDVETA%VETAF(NFLEVG)

! Preparatory stage setting values for 1-th iteration
JITER=1

DO JLEV=1,NFLEVG
  DO JROF=KST,KPROF

    ! * computations on horizontal plans.

    ! convert arrival point wind to Cartesian space for later averaging
    ! Arrival point wind
    ZPU0= +PB2(JROF,MSLB2URL5+JLEV-1)*YDGSGEOM%GNORDM(JROF)&
     &    -PB2(JROF,MSLB2VRL5+JLEV-1)*YDGSGEOM%GNORDL(JROF)  
    ZPV0= +PB2(JROF,MSLB2URL5+JLEV-1)*YDGSGEOM%GNORDL(JROF)&
     &    +PB2(JROF,MSLB2VRL5+JLEV-1)*YDGSGEOM%GNORDM(JROF)
    ! Arrival point wind converted to Cartesian space
    ZU0(JROF,JLEV)=-ZPU0*YDGSGEOM%GESLO(JROF)             &
     &             -ZPV0*YDGSGEOM%GECLO(JROF)*YDGSGEOM%GEMU(JROF)
    ZV0(JROF,JLEV)= ZPU0*YDGSGEOM%GECLO(JROF)             &
     &             -ZPV0*YDGSGEOM%GESLO(JROF)*YDGSGEOM%GEMU(JROF)
    ZZ0(JROF,JLEV)= ZPV0*YDGSGEOM%GSQM2(JROF)

    ZPU=PB2(JROF,MSLB2URL5+JLEV-1)
    ZPV=PB2(JROF,MSLB2VRL5+JLEV-1)

    ZNOR2=( ZPU*ZPU + ZPV*ZPV )
    PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI,JITER,3)=1.0_JPRB-ZDTS22*ZNOR2
    ZSPHSV=ZDTSA*(1.0_JPRB-ZDTS62*ZNOR2)
    ZINT=(ZPU*YDGSGEOM%GNORDL(JROF)+ZPV*YDGSGEOM%GNORDM(JROF))*ZSPHSV  
    PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_SINLA,JITER,3)= &
     & YDGSGEOM%GEMU(JROF)*PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI,JITER,3)-ZINT*YDGSGEOM%GSQM2(JROF)
    PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO,JITER,3)= &
     & YDGSGEOM%GSQM2(JROF)*PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI,JITER,3)+ZINT*YDGSGEOM%GEMU(JROF)
    PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO,JITER,3)=-( ZPU*YDGSGEOM%GNORDM(JROF)&
     & -ZPV*YDGSGEOM%GNORDL(JROF) )*ZSPHSV  

    ! * computations on a vertical.
    ZPW0=PB2(JROF,MSLB2WRL5+JLEV-1)
    ZLEVO=YDVETA%VETAF(JLEV)-RTDT*ZPW0
    ZLEVO=MIN(ZVETAOX,MAX(ZVETAON,ZLEVO))
    PLEV5(JROF,JLEV,JITER,3)=ZLEVO

  ENDDO
ENDDO

! Interpolation

IHVI=0
ITIP=1
LLINTV=.TRUE.
IROT = 0 

! Extra interpolation for t+dt extrapolated wind.
! Used for stages 1 (half) and 3 (full)
CALL LARCINA(YDGEOMETRY,YDML_DYN,KST,KPROF,YDSL,IHVI,KSTABUF,LFINDVSEP,LLSLHD,LLSLHDQUAD,LLINTV,&
 & ITIP,IROT,.FALSE.,PLSDEPI,&
 & KIBL,PSCO5(:,:,:,JITER,3),PLEV5(:,:,JITER,3),&
 & ZUN_KAPPA,ZUN_KAPPAT,ZUN_KAPPAM,ZUN_KAPPAH,&
 & ZUN_STDDISU,ZUN_STDDISV,ZUN_STDDISW,&
 & PB15(1,MSLB1UR005),PB15(1,MSLB1VR005),PB15(1,MSLB1ZR005),PB15(1,MSLB1WR005),&
 & IUN_VSEPC,IUN_VSEPL,PCCO5(:,:,:,JITER,3),&
 & PUS5(:,:,3,JITER),PVS5(:,:,3,JITER),PZS5(:,:,3,JITER),PWS5(:,:,3,JITER),ZUN_WFSM,&
 & KL0(:,:,:,JITER,3),KLH0,KLEV,PLSCAW5(:,:,:,JITER,3),PRSCAW5(:,:,:,JITER,3),&
 & KDEP(:,:,JITER,3),KNOWENO(:,:,JITER,3))  


!  Iterations starts


DO JITER=2,NITMP

  !*       2.1   DETERMINATION OF THE MEDIUM POINT "M" OR THE ORIGIN POINT "O".

  ! Computation of the norm of the real "(a/rs)*wind" vector.
  ! Computation of the angle (PHI=DT . NORM OF V ON RADIUS)**2
  ! then computation of the coordinates of the medium point "M".
  ! If (LLO=T or LELTRA=T) the origin point "O" is computed
  ! instead of the medium point "M".


  ! STAGE 1 (arrival point, t+dt) ============================================================

  ! Second interpolation to complete extrapolated wind by SETTLS
  CALL LARCINA(YDGEOMETRY,YDML_DYN,KST,KPROF,YDSL,IHVI,KSTABUF,LFINDVSEP,LLSLHD,LLSLHDQUAD,LLINTV,&
   & ITIP,IROT,.FALSE.,PLSDEPI,&
   & KIBL,PSCO5(:,:,:,JITER-1,3),PLEV5(:,:,JITER-1,3),&
   & ZUN_KAPPA,ZUN_KAPPAT,ZUN_KAPPAM,ZUN_KAPPAH,&
   & ZUN_STDDISU,ZUN_STDDISV,ZUN_STDDISW,&
   & PB15(1,MSLB1UR05),PB15(1,MSLB1VR05),PB15(1,MSLB1ZR05),PB15(1,MSLB1WR05),&
   & IUN_VSEPC,IUN_VSEPL,PCCO5(:,:,:,JITER,1),&
   & PUF05(:,:,1,JITER),PVF05(:,:,1,JITER),PZF05(:,:,1,JITER),PWF05(:,:,1,JITER),ZUN_WFSM,&
   & KL0(:,:,:,JITER,1),KLH0,KLEV,PLSCAW5(:,:,:,JITER,1),PRSCAW5(:,:,:,JITER,1),&
   & KDEP(:,:,JITER,1),KNOWENO(:,:,JITER,1))  

  IF (JITER == 2) THEN
    ! Stage 2 (middle point, t+0.5dt)
    ! Averaging wind by SETTLS to arrival point in Cartesian space
    PUS5(KST:KPROF,1:NFLEVG,2,JITER-1)=0.5_JPRB*(PUF05(KST:KPROF,1:NFLEVG,1,JITER) &
      &                                          +ZU0(KST:KPROF,1:NFLEVG))
    PVS5(KST:KPROF,1:NFLEVG,2,JITER-1)=0.5_JPRB*(PVF05(KST:KPROF,1:NFLEVG,1,JITER) &
      &                                          +ZV0(KST:KPROF,1:NFLEVG))
    PZS5(KST:KPROF,1:NFLEVG,2,JITER-1)=0.5_JPRB*(PZF05(KST:KPROF,1:NFLEVG,1,JITER) &
      &                                          +ZZ0(KST:KPROF,1:NFLEVG))
    PWS5(KST:KPROF,1:NFLEVG,2,JITER-1)=0.5_JPRB*(PWF05(KST:KPROF,1:NFLEVG,1,JITER) &
      &                                          +PB2(KST:KPROF,MSLB2WRL5:MSLB2WRL5+NFLEVG-1))
  ENDIF

  ! Extrapolating wind by SETTLS to arrival point in Cartesian space
  PUS5(KST:KPROF,1:NFLEVG,1,JITER)=PUF05(KST:KPROF,1:NFLEVG,1,JITER)-PUS5(KST:KPROF,1:NFLEVG,3,JITER-1) &
    &                              +ZU0(KST:KPROF,1:NFLEVG)
  PVS5(KST:KPROF,1:NFLEVG,1,JITER)=PVF05(KST:KPROF,1:NFLEVG,1,JITER)-PVS5(KST:KPROF,1:NFLEVG,3,JITER-1) &
    &                              +ZV0(KST:KPROF,1:NFLEVG)
  PZS5(KST:KPROF,1:NFLEVG,1,JITER)=PZF05(KST:KPROF,1:NFLEVG,1,JITER)-PZS5(KST:KPROF,1:NFLEVG,3,JITER-1) &
    &                              +ZZ0(KST:KPROF,1:NFLEVG)
  PWS5(KST:KPROF,1:NFLEVG,1,JITER)=PWF05(KST:KPROF,1:NFLEVG,1,JITER)-PWS5(KST:KPROF,1:NFLEVG,3,JITER-1) &
    &                              +PB2(KST:KPROF,MSLB2WRL5:MSLB2WRL5+NFLEVG-1)


  ! STAGE 2 (midpoint, t+0.5dt) ==============================================================
  ! Here we however compute the departure point position using then obvious SETTLS averaging into M

  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF

      ! * computations on horizontal plans.
      ZPU0=(5._JPRB*PUS5(JROF,JLEV,1,JITER)+8._JPRB*PUS5(JROF,JLEV,2,JITER-1) &
       &           -PUS5(JROF,JLEV,3,JITER-1))/12._JPRB
      ZPV0=(5._JPRB*PVS5(JROF,JLEV,1,JITER)+8._JPRB*PVS5(JROF,JLEV,2,JITER-1) &
       &           -PVS5(JROF,JLEV,3,JITER-1))/12._JPRB
      ZPZ0=(5._JPRB*PZS5(JROF,JLEV,1,JITER)+8._JPRB*PZS5(JROF,JLEV,2,JITER-1) &
       &           -PZS5(JROF,JLEV,3,JITER-1))/12._JPRB
      ! Convert wind back to lat,lon at Arrival point
      ZPU = -ZPU0*YDGSGEOM%GESLO(JROF)                     &
       &    +ZPV0*YDGSGEOM%GECLO(JROF)
      ZPV = -ZPU0*YDGSGEOM%GECLO(JROF)*YDGSGEOM%GEMU(JROF) &
       &    -ZPV0*YDGSGEOM%GESLO(JROF)*YDGSGEOM%GEMU(JROF) &
       &    +ZPZ0*YDGSGEOM%GSQM2(JROF)

      ZNOR2             =(ZPU*ZPU+ZPV*ZPV)
      PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI,JITER,2) =1.0_JPRB-ZDTS22*ZNOR2
      ZSPHSV            =ZDTSA*(1.0_JPRB-ZDTS62*ZNOR2)
      ZINT              =ZPV*ZSPHSV
      PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_SINLA,JITER,2) = &
       &  YDGSGEOM%GEMU(JROF)*PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI,JITER,2)-ZINT*YDGSGEOM%GSQM2(JROF)
      PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO,JITER,2) = &
       &  YDGSGEOM%GSQM2(JROF)*PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI,JITER,2)+ZINT*YDGSGEOM%GEMU(JROF)
      PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO,JITER,2) =-ZPU*ZSPHSV

      ! * computations on a vertical.
      ZPW0 =(5._JPRB*PWS5(JROF,JLEV,1,JITER)+8._JPRB*PWS5(JROF,JLEV,2,JITER-1) &
       &                      -PWS5(JROF,JLEV,3,JITER-1))/12._JPRB
      ZLEVO=YDVETA%VETAF(JLEV)-RTDT*ZPW0
      ZLEVO=MIN(ZVETAOX,MAX(ZVETAON,ZLEVO))
      PLEV5(JROF,JLEV,JITER,2)=ZLEVO

    ENDDO
  ENDDO

  ! Interpolation for stage 2
  CALL LARCINA(YDGEOMETRY,YDML_DYN,KST,KPROF,YDSL,IHVI,KSTABUF,LFINDVSEP,LLSLHD,LLSLHDQUAD,LLINTV,&
   & ITIP,IROT,.FALSE.,PLSDEPI,&
   & KIBL,PSCO5(:,:,:,JITER,2),PLEV5(:,:,JITER,2),&
   & ZUN_KAPPA,ZUN_KAPPAT,ZUN_KAPPAM,ZUN_KAPPAH,&
   & ZUN_STDDISU,ZUN_STDDISV,ZUN_STDDISW,&
   & PB15(1,MSLB1UR05),PB15(1,MSLB1VR05),PB15(1,MSLB1ZR05),PB15(1,MSLB1WR05),&
   & IUN_VSEPC,IUN_VSEPL,PCCO5(:,:,:,JITER,2),&
   & PUF05(:,:,2,JITER),PVF05(:,:,2,JITER),PZF05(:,:,2,JITER),PWF05(:,:,2,JITER),ZUN_WFSM,&
   & KL0(:,:,:,JITER,2),KLH0,KLEV,PLSCAW5(:,:,:,JITER,2),PRSCAW5(:,:,:,JITER,2),&
   & KDEP(:,:,JITER,2),KNOWENO(:,:,JITER,2))  
  ! New wind stage 2
  ! Stage 2 (middle point, t+0.5dt)
  ! Averaging wind by SETTLS to arrival point in Cartesian space
  PUS5(KST:KPROF,1:NFLEVG,2,JITER)=0.5_JPRB*(PUF05(KST:KPROF,1:NFLEVG,2,JITER) &
   &                                         +ZU0(KST:KPROF,1:NFLEVG))
  PVS5(KST:KPROF,1:NFLEVG,2,JITER)=0.5_JPRB*(PVF05(KST:KPROF,1:NFLEVG,2,JITER) &
   &                                         +ZV0(KST:KPROF,1:NFLEVG))
  PZS5(KST:KPROF,1:NFLEVG,2,JITER)=0.5_JPRB*(PZF05(KST:KPROF,1:NFLEVG,2,JITER) &
   &                                         +ZZ0(KST:KPROF,1:NFLEVG))
  PWS5(KST:KPROF,1:NFLEVG,2,JITER)=0.5_JPRB*(PWF05(KST:KPROF,1:NFLEVG,2,JITER) &
   &                                         +PB2(KST:KPROF,MSLB2WRL5:MSLB2WRL5+NFLEVG-1))


  ! STAGE 3 (departure point, t) =============================================================
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF

      ! * computations on horizontal plans.
      ZPU0=(PUS5(JROF,JLEV,1,JITER)+4._JPRB*PUS5(JROF,JLEV,2,JITER) &
       &   +PUS5(JROF,JLEV,3,JITER-1))/6._JPRB
      ZPV0=(PVS5(JROF,JLEV,1,JITER)+4._JPRB*PVS5(JROF,JLEV,2,JITER) &
       &   +PVS5(JROF,JLEV,3,JITER-1))/6._JPRB
      ZPZ0=(PZS5(JROF,JLEV,1,JITER)+4._JPRB*PZS5(JROF,JLEV,2,JITER) &
       &   +PZS5(JROF,JLEV,3,JITER-1))/6._JPRB
      ! Convert wind back to lat,lon at Arrival point
      ZPU = -ZPU0*YDGSGEOM%GESLO(JROF)                     &
       &    +ZPV0*YDGSGEOM%GECLO(JROF)
      ZPV = -ZPU0*YDGSGEOM%GECLO(JROF)*YDGSGEOM%GEMU(JROF) &
       &    -ZPV0*YDGSGEOM%GESLO(JROF)*YDGSGEOM%GEMU(JROF) &
       &    +ZPZ0*YDGSGEOM%GSQM2(JROF)

      ZNOR2             =(ZPU*ZPU+ZPV*ZPV)
      PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI,JITER,3) =1.0_JPRB-ZDTS22*ZNOR2
      ZSPHSV            =ZDTSA*(1.0_JPRB-ZDTS62*ZNOR2)
      ZINT              =ZPV*ZSPHSV
      PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_SINLA,JITER,3) = &
       &  YDGSGEOM%GEMU(JROF)*PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI,JITER,3)-ZINT*YDGSGEOM%GSQM2(JROF)
      PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO,JITER,3) = &
       &  YDGSGEOM%GSQM2(JROF)*PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI,JITER,3)+ZINT*YDGSGEOM%GEMU(JROF)
      PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO,JITER,3) =-ZPU*ZSPHSV

      ! * computations on a vertical.
      ZPW0=(PWS5(JROF,JLEV,1,JITER)+4._JPRB*PWS5(JROF,JLEV,2,JITER) &
       &              +PWS5(JROF,JLEV,3,JITER-1))/6._JPRB
      ZLEVO=YDVETA%VETAF(JLEV)-RTDT*ZPW0
      ZLEVO=MIN(ZVETAOX,MAX(ZVETAON,ZLEVO))
      PLEV5(JROF,JLEV,JITER,3)=ZLEVO

    ENDDO
  ENDDO

  ! Interpolation for stage 3
  IF (JITER /= NITMP) THEN
    CALL LARCINA(YDGEOMETRY,YDML_DYN,KST,KPROF,YDSL,IHVI,KSTABUF,LFINDVSEP,LLSLHD,LLSLHDQUAD,LLINTV,&
     & ITIP,IROT,.FALSE.,PLSDEPI,&
     & KIBL,PSCO5(:,:,:,JITER,3),PLEV5(:,:,JITER,3),&
     & ZUN_KAPPA,ZUN_KAPPAT,ZUN_KAPPAM,ZUN_KAPPAH,&
     & ZUN_STDDISU,ZUN_STDDISV,ZUN_STDDISW,&
     & PB15(1,MSLB1UR005),PB15(1,MSLB1VR005),PB15(1,MSLB1ZR005),PB15(1,MSLB1WR005),&
     & IUN_VSEPC,IUN_VSEPL,PCCO5(:,:,:,JITER,3),&
     & PUS5(:,:,3,JITER),PVS5(:,:,3,JITER),PZS5(:,:,3,JITER),PWS5(:,:,3,JITER),ZUN_WFSM,&
     & KL0(:,:,:,JITER,3),KLH0,KLEV,PLSCAW5(:,:,:,JITER,3),PRSCAW5(:,:,:,JITER,3),&
     & KDEP(:,:,JITER,3),KNOWENO(:,:,JITER,3))  
    ! New wind stage 3 overwrites Z[x]S
  ENDIF


ENDDO   ! end of JITER

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF THE ORIGIN POINT "O" COORDINATES.
!              ------------------------------------------------

IF (.NOT.(LLO.OR.YDML_DYN%YRDYNA%LELTRA)) THEN

  ! this case may occur only for SL3TL scheme.
  IF (LTWOTL) CALL ABOR1(' LARMES_RK5 3 ')
  CALL ABOR1(' LARMES_RK5: 3TL scheme not coded for RK4 trajectory research. ')

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LARMES_RK5',1,ZHOOK_HANDLE)
END SUBROUTINE LARMES_RK5
