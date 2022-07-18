SUBROUTINE LARMESTL(YDGEOMETRY,YDRIP,YDML_DYN,KST,KPROF,YDSL,KSTABUF,PB1,PB15,PB2,&
 & PLSDEPI,KIBL,&
 & PSCO,PLEV,PCCO,&
 & KL0,KLH0,KLEV,PLSCAW,PRSCAW,&
 & PSCO5,PLEV5,PCCO5,PUF5,PVF5,PZF5,PLSCAW5,PRSCAW5)  

!**** *LARMESTL - semi-LAgrangian scheme:  (tangent-linear version)
!                Research of the origin point on the Sphere.

!     Purpose.
!     --------

!      The computation of the location of the origin point "O" of
!     the lagrangian trajectory is performed by an iterative
!     method described by Robert and adapted to the sphere by M. Rochas.
!     Trajectories are great circles.

!**   Interface.
!     ----------
!        *CALL* *LARMESTL(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST      - first element of arrays where computations are performed.
!          KPROF    - depth of work.
!          YDSL     - SL_STRUCT definition
!          KSTABUF  - for a latitude IGL, KSTABUF(IGL) is the
!                     address of the element corresponding to
!                     (ILON=1,IGL) in the NPROMA arrays.
!          PB1      - SLBUF1-buffer for interpolations.
!          PB15     - SLBUF1-buffer for interpolations (trajectory).
!          PB2      - SLBUF2-buf to communicate info from non lag. to lag. dyn.
!          PLSDEPI  - (Number of points by latitude) / (2 * PI) .
!          KIBL     - index into YRCSGEOM/YRGSGEOM instances in YDGEOMETRY

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
!          PUF5     - interpolated U-wind.                           (in)
!          PVF5     - interpolated V-wind.                           (in)
!          PZF5     - interpolated Z-wind.                           (in)
!          PLSCAW5  - linear weights (distances) for interpolations  (out)
!          PRSCAW5  - non-linear weights for interpolations.         (out)


!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------
!        Calls  LARCINATL.
!        Is called by LAPINEATL (3D model) 

!     Reference.
!     ----------

!     Author.
!     -------
!      C. Temperton (ECMWF)
!      Original : 99/07/12

!     Modifications.
!     --------------
!      Modified 01-08-30 by K. YESSAD: pruning and some other cleanings.
!      Modified 03-04-09 by C. Temperton: some trajectory recomputations
!                                         already done by calling SLADJTR
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Modified 05-01-11 by C. Temperton: streamlining.
!      Modified 05-01-14 by C. Temperton: transfer computations from
!                              start of LAIDEP to end of LARMES.
!      Modified 06-09-04 by F. Vana: - bugfix of vertical bounding when
!                                      LSETTLS=LELTRA=.F.
!                                    - removing useless LDPLANE blocks
!      Modified 08-Jun-27 by F. Vana: update of arguments for LARCINATL
!      K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!      K. Yessad (Aug 2009): always use root (QX,QY) for (p,q) variables names
!      K. Yessad (Nov 2009): keep LLO.OR.LELTRA=T only.
!      K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!      G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      G. Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
!      F. Vana  13-Feb-2014  Updated arguments list for LARCINATL
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      F. Vana    21-Nov-2017: Options LSLDP_CURV and LHOISLT
!      F. Vana October 2018: Extended LSLDP_CURV.
!      F. Vana    26-Feb-2019: Vertical quintic interpolation
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
! End Modifications
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE SC2PRG_MOD         , ONLY : SC2PRG
USE YOMCT0             , ONLY : LTWOTL

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
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTABUF(YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB1(YDSL%NASLB1,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB15(YDSL%NASLB1,YDML_DYN%YRPTRSLB15%NFLDSLB15)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSDEPI(YDSL%NDGSAH:YDSL%NDGENH)
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTSCO%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCCO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3,YDML_DYN%YRDYN%NITMP)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KLH0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
INTEGER(KIND=JPIM),INTENT(INOUT) :: KLEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRSCAW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTSCO%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLEV5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP-1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP-1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP-1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM, &
 & YDML_DYN%YRDYN%NITMP)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IROT, ITIP, JITER, JLEV, JROF

LOGICAL :: LLO
LOGICAL :: LLSLHD,LLSLHDQUAD,LLSLHD_OLD

REAL(KIND=JPRB) :: ZDTS22, ZDTS62, ZDTSA, ZINT,&
 & ZLEVO, ZNOR2, ZPU, ZPV, ZSPHSV,&
 & ZVETAON, ZVETAOX, ZNOR25, ZSPHSV5, ZLEVO5, ZPU5, ZPV5

REAL(KIND=JPRB) :: ZUF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZVF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZZF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZWF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZU0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZV0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZZ0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPU0, ZPV0

! unused arguments to LARCINATL 
!  a/ input variables
REAL(KIND=JPRB) :: ZUN_KAPPA(1), ZUN_KAPPAT(1)
REAL(KIND=JPRB) :: ZUN_KAPPA5(1), ZUN_KAPPAT5(1)
!  b/ output variable
INTEGER(KIND=JPIM) :: INOWENO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB), POINTER :: ZSLB1UR0(:), ZSLB1VR0(:), ZSLB1WR0(:), ZSLB1ZR0(:)
REAL(KIND=JPRB), POINTER :: ZSLB1UR05(:), ZSLB1VR05(:), ZSLB1WR05(:), ZSLB1ZR05(:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "larcinatl.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LARMESTL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL),YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDDYN=>YDML_DYN%YRDYN, &
 & YDPTRSLB1=>YDML_DYN%YRPTRSLB1, YDPTRSLB2=>YDML_DYN%YRPTRSLB2,YDPTRSLB15=>YDML_DYN%YRPTRSLB15, &
 & YDTCCO=>YDML_DYN%YYTCCO, YDTSCO=>YDML_DYN%YYTSCO)

ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & NITMP=>YDDYN%NITMP, VETAON=>YDDYN%VETAON, VETAOX=>YDDYN%VETAOX, &
 & LSLDP_CURV=>YDDYN%LSLDP_CURV, &
 & MSLB1UR0=>YDPTRSLB1%MSLB1UR0, MSLB1VR0=>YDPTRSLB1%MSLB1VR0, &
 & MSLB1ZR0=>YDPTRSLB1%MSLB1ZR0, MSLB1ZR05=>YDPTRSLB15%MSLB1ZR05, &
 & MSLB1WR0=>YDPTRSLB1%MSLB1WR0, NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & MSLB1UR05=>YDPTRSLB15%MSLB1UR05, MSLB1VR05=>YDPTRSLB15%MSLB1VR05, &
 & MSLB1WR05=>YDPTRSLB15%MSLB1WR05, NFLDSLB15=>YDPTRSLB15%NFLDSLB15, &
 & MSLB2URL=>YDPTRSLB2%MSLB2URL, MSLB2URL5=>YDPTRSLB2%MSLB2URL5, &
 & MSLB2VRL=>YDPTRSLB2%MSLB2VRL, MSLB2VRL5=>YDPTRSLB2%MSLB2VRL5, &
 & MSLB2WRL=>YDPTRSLB2%MSLB2WRL, NFLDSLB2=>YDPTRSLB2%NFLDSLB2, &
 & RDTS22=>YDRIP%RDTS22, RDTS62=>YDRIP%RDTS62, RDTSA=>YDRIP%RDTSA, &
 & RTDT=>YDRIP%RTDT)

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
IF (LTWOTL) THEN
  LLO=.NOT.YDML_DYN%YRDYNA%LELTRA
ELSE
  LLO=.FALSE.
ENDIF

IF(LLO.OR.YDML_DYN%YRDYNA%LELTRA) THEN
  ZDTS22=RDTS22*4
  ZDTSA=RDTSA*2
  ZDTS62=RDTS62*4
ELSE
  CALL ABOR1(' LARMESTL: LLO.OR.LELTRA=F no longer supported')
ENDIF

! deactivate computation of SLHD weights
LLSLHD     =.FALSE.
LLSLHDQUAD =.FALSE.
LLSLHD_OLD =.FALSE.

!     ------------------------------------------------------------------

!*       2.    ITERATIONS.
!              -----------

ZVETAON=(1.0_JPRB-VETAON)*YDVETA%VETAH(0)+VETAON*YDVETA%VETAF(1)
ZVETAOX=(1.0_JPRB-VETAOX)*YDVETA%VETAH(NFLEVG)+VETAOX*YDVETA%VETAF(NFLEVG)

DO JITER=1,NITMP

  !*       2.1   DETERMINATION OF THE ORIGIN POINT "O".

  ! Computation of the norm of the real "(a/rs)*wind" vector.
  ! Computation of the angle (PHI=DT . NORM OF V ON RADIUS)**2
  ! then computation of the coordinates of the origin point "O".

  IF (JITER == 1) THEN

    DO JLEV=1,NFLEVG
      DO JROF=KST,KPROF

        ! * computations on horizontal plans.

        ! When relevant convert arrival point wind to cartesian space
        IF (LSLDP_CURV.AND.(.NOT.YDML_DYN%YRDYNA%LELTRA)) THEN  
          ! Arrival point wind
          ZPU0= +PB2(JROF,MSLB2URL+JLEV-1)*YDGSGEOM%GNORDM(JROF)&
           &    -PB2(JROF,MSLB2VRL+JLEV-1)*YDGSGEOM%GNORDL(JROF)  
          ZPV0= +PB2(JROF,MSLB2URL+JLEV-1)*YDGSGEOM%GNORDL(JROF)&
           &    +PB2(JROF,MSLB2VRL+JLEV-1)*YDGSGEOM%GNORDM(JROF)
          ! Arrival point wind converted to Cartesian space
          ZU0(JROF,JLEV)=-ZPU0*YDGSGEOM%GESLO(JROF)             &
           &             -ZPV0*YDGSGEOM%GECLO(JROF)*YDGSGEOM%GEMU(JROF)
          ZV0(JROF,JLEV)= ZPU0*YDGSGEOM%GECLO(JROF)             &
           &             -ZPV0*YDGSGEOM%GESLO(JROF)*YDGSGEOM%GEMU(JROF)
          ZZ0(JROF,JLEV)= ZPV0*YDGSGEOM%GSQM2(JROF)
        ENDIF

        ZNOR25=( PB2(JROF,MSLB2URL5+JLEV-1)*PB2(JROF,MSLB2URL5+JLEV-1)&
         & +PB2(JROF,MSLB2VRL5+JLEV-1)*PB2(JROF,MSLB2VRL5+JLEV-1) )  
        ZSPHSV5=ZDTSA*(1.0_JPRB-ZDTS62*ZNOR25)
        ZNOR2=2.0_JPRB*(PB2(JROF,MSLB2URL+JLEV-1)*PB2(JROF,MSLB2URL5+JLEV-1)+&
         & PB2(JROF,MSLB2VRL+JLEV-1)*PB2(JROF,MSLB2VRL5+JLEV-1))  
        PSCO(JROF,JLEV,YDTSCO%M_COPHI)=-ZDTS22*ZNOR2
        ZSPHSV=-ZDTSA*ZDTS62*ZNOR2
        ZINT=( PB2(JROF,MSLB2URL+JLEV-1)*YDGSGEOM%GNORDL(JROF)&
         & +PB2(JROF,MSLB2VRL+JLEV-1)*YDGSGEOM%GNORDM(JROF))*ZSPHSV5 +&
         & ( PB2(JROF,MSLB2URL5+JLEV-1)*YDGSGEOM%GNORDL(JROF)&
         & +PB2(JROF,MSLB2VRL5+JLEV-1)*YDGSGEOM%GNORDM(JROF))*ZSPHSV  
        PSCO(JROF,JLEV,YDTSCO%M_SINLA)= YDGSGEOM%GEMU(JROF)*PSCO(JROF,JLEV,YDTSCO%M_COPHI) &
         &                              -ZINT*YDGSGEOM%GSQM2(JROF)
        PSCO(JROF,JLEV,YDTSCO%M_COSCO)= YDGSGEOM%GSQM2(JROF)*PSCO(JROF,JLEV,YDTSCO%M_COPHI) &
         &                              +ZINT* YDGSGEOM%GEMU(JROF)
        PSCO(JROF,JLEV,YDTSCO%M_SINCO)=-( PB2(JROF,MSLB2URL+JLEV-1)*YDGSGEOM%GNORDM(JROF)&
         & -PB2(JROF,MSLB2VRL+JLEV-1)*YDGSGEOM%GNORDL(JROF) )*ZSPHSV5 -&
         & ( PB2(JROF,MSLB2URL5+JLEV-1)*YDGSGEOM%GNORDM(JROF)&
         & -PB2(JROF,MSLB2VRL5+JLEV-1)*YDGSGEOM%GNORDL(JROF) )*ZSPHSV  

        ! * computations on a vertical.

        ZWF(JROF,JLEV)=PB2(JROF,MSLB2WRL+JLEV-1)
        ZLEVO=-RTDT*PB2(JROF,MSLB2WRL+JLEV-1)
        IF(LLO.OR.YDML_DYN%YRDYNA%LELTRA) THEN
          ZLEVO5=PLEV5(JROF,JLEV,JITER)
          IF (ZLEVO5 <= ZVETAON) THEN
            ZLEVO=0.0_JPRB
          ELSEIF (ZLEVO5 >= ZVETAOX) THEN
            ZLEVO=0.0_JPRB
          ENDIF
          PLEV(JROF,JLEV)=ZLEVO
        ENDIF

      ENDDO
    ENDDO

  ELSE

    DO JLEV=1,NFLEVG
      DO JROF=KST,KPROF

        ! * ZPU,ZPV are the coordinates of VM in the local repere
        !   related to F.

        IF(LSLDP_CURV.AND.YDML_DYN%YRDYNA%LELTRA) THEN
          ! Convert Departure point wind back to lat,lon at Arrival point
          ZPU5 = -PUF5(JROF,JLEV,JITER-1)*YDGSGEOM%GESLO(JROF)                      &
           &     +PVF5(JROF,JLEV,JITER-1)*YDGSGEOM%GECLO(JROF)
          ZPV5 = -PUF5(JROF,JLEV,JITER-1)*YDGSGEOM%GECLO(JROF)*YDGSGEOM%GEMU(JROF)  &
           &     -PVF5(JROF,JLEV,JITER-1)*YDGSGEOM%GESLO(JROF)*YDGSGEOM%GEMU(JROF)  &
           &     +PZF5(JROF,JLEV,JITER-1)*YDGSGEOM%GSQM2(JROF)

          ZPU =  -ZUF(JROF,JLEV)*YDGSGEOM%GESLO(JROF)                      &
           &     +ZVF(JROF,JLEV)*YDGSGEOM%GECLO(JROF)
          ZPV =  -ZUF(JROF,JLEV)*YDGSGEOM%GECLO(JROF)*YDGSGEOM%GEMU(JROF)  &
           &     -ZVF(JROF,JLEV)*YDGSGEOM%GESLO(JROF)*YDGSGEOM%GEMU(JROF)  &
           &     +ZZF(JROF,JLEV)*YDGSGEOM%GSQM2(JROF)
        ELSEIF(YDML_DYN%YRDYNA%LELTRA) THEN
          ZPU5= PCCO5(JROF,JLEV,YDTCCO%M_RQX,JITER-1)*PUF5(JROF,JLEV,JITER-1)&
           & +PCCO5(JROF,JLEV,YDTCCO%M_RQY,JITER-1)*PVF5(JROF,JLEV,JITER-1)  
          ZPU= PCCO(JROF,JLEV,YDTCCO%M_RQX)*PUF5(JROF,JLEV,JITER-1)&
           & +PCCO5(JROF,JLEV,YDTCCO%M_RQX,JITER-1)*ZUF(JROF,JLEV)&
           & +PCCO(JROF,JLEV,YDTCCO%M_RQY)*PVF5(JROF,JLEV,JITER-1)&
           & +PCCO5(JROF,JLEV,YDTCCO%M_RQY,JITER-1)*ZVF(JROF,JLEV)  
          ZPV5=-PCCO5(JROF,JLEV,YDTCCO%M_RQY,JITER-1)*PUF5(JROF,JLEV,JITER-1)&
           & +PCCO5(JROF,JLEV,YDTCCO%M_RQX,JITER-1)*PVF5(JROF,JLEV,JITER-1)  
          ZPV=-PCCO(JROF,JLEV,YDTCCO%M_RQY)*PUF5(JROF,JLEV,JITER-1)&
           & -PCCO5(JROF,JLEV,YDTCCO%M_RQY,JITER-1)*ZUF(JROF,JLEV)&
           & +PCCO(JROF,JLEV,YDTCCO%M_RQX)*PVF5(JROF,JLEV,JITER-1)&
           & +PCCO5(JROF,JLEV,YDTCCO%M_RQX,JITER-1)*ZVF(JROF,JLEV)  
        ELSEIF(LSLDP_CURV.AND.LLO) THEN
          ! Averaging wind to Midpoint in Cartesian space
          ZUF(JROF,JLEV)=0.5_JPRB*(ZUF(JROF,JLEV)+ZU0(JROF,JLEV))
          ZVF(JROF,JLEV)=0.5_JPRB*(ZVF(JROF,JLEV)+ZV0(JROF,JLEV))
          ZZF(JROF,JLEV)=0.5_JPRB*(ZZF(JROF,JLEV)+ZZ0(JROF,JLEV))

          ! Convert Midpoint wind back to lat,lon at Arrival point
          ZPU5 = -PUF5(JROF,JLEV,JITER-1)*YDGSGEOM%GESLO(JROF)                      &
           &     +PVF5(JROF,JLEV,JITER-1)*YDGSGEOM%GECLO(JROF)
          ZPV5 = -PUF5(JROF,JLEV,JITER-1)*YDGSGEOM%GECLO(JROF)*YDGSGEOM%GEMU(JROF)  &
           &     -PVF5(JROF,JLEV,JITER-1)*YDGSGEOM%GESLO(JROF)*YDGSGEOM%GEMU(JROF)  &
           &     +PZF5(JROF,JLEV,JITER-1)*YDGSGEOM%GSQM2(JROF)

          ZPU =  -ZUF(JROF,JLEV)*YDGSGEOM%GESLO(JROF)                      &
           &     +ZVF(JROF,JLEV)*YDGSGEOM%GECLO(JROF)
          ZPV =  -ZUF(JROF,JLEV)*YDGSGEOM%GECLO(JROF)*YDGSGEOM%GEMU(JROF)  &
           &     -ZVF(JROF,JLEV)*YDGSGEOM%GESLO(JROF)*YDGSGEOM%GEMU(JROF)  &
           &     +ZZF(JROF,JLEV)*YDGSGEOM%GSQM2(JROF)
          ZWF(JROF,JLEV)=0.5_JPRB*(ZWF(JROF,JLEV)+PB2(JROF,MSLB2WRL+JLEV-1))
        ELSEIF(LLO) THEN
          ZPU5= 0.5_JPRB*(PCCO5(JROF,JLEV,YDTCCO%M_RQX,JITER-1)*PUF5(JROF,JLEV,JITER-1)&
           & +PCCO5(JROF,JLEV,YDTCCO%M_RQY,JITER-1)*PVF5(JROF,JLEV,JITER-1)&
           & +PB2(JROF,MSLB2URL5+JLEV-1)*YDGSGEOM%GNORDM(JROF)&
           & -PB2(JROF,MSLB2VRL5+JLEV-1)*YDGSGEOM%GNORDL(JROF))  
          ZPU= 0.5_JPRB*(PCCO(JROF,JLEV,YDTCCO%M_RQX)*PUF5(JROF,JLEV,JITER-1)&
           & +PCCO5(JROF,JLEV,YDTCCO%M_RQX,JITER-1)*ZUF(JROF,JLEV)&
           & +PCCO(JROF,JLEV,YDTCCO%M_RQY)*PVF5(JROF,JLEV,JITER-1)&
           & +PCCO5(JROF,JLEV,YDTCCO%M_RQY,JITER-1)*ZVF(JROF,JLEV)&
           & +PB2(JROF,MSLB2URL+JLEV-1)*YDGSGEOM%GNORDM(JROF)&
           & -PB2(JROF,MSLB2VRL+JLEV-1)*YDGSGEOM%GNORDL(JROF))  
          ZPV5=0.5_JPRB*(-PCCO5(JROF,JLEV,YDTCCO%M_RQY,JITER-1)*PUF5(JROF,JLEV,JITER-1)&
           & +PCCO5(JROF,JLEV,YDTCCO%M_RQX,JITER-1)*PVF5(JROF,JLEV,JITER-1)&
           & +PB2(JROF,MSLB2URL5+JLEV-1)*YDGSGEOM%GNORDL(JROF)&
           & +PB2(JROF,MSLB2VRL5+JLEV-1)*YDGSGEOM%GNORDM(JROF))  
          ZPV=0.5_JPRB*(-PCCO(JROF,JLEV,YDTCCO%M_RQY)*PUF5(JROF,JLEV,JITER-1)&
           & -PCCO5(JROF,JLEV,YDTCCO%M_RQY,JITER-1)*ZUF(JROF,JLEV)&
           & +PCCO(JROF,JLEV,YDTCCO%M_RQX)*PVF5(JROF,JLEV,JITER-1)&
           & +PCCO5(JROF,JLEV,YDTCCO%M_RQX,JITER-1)*ZVF(JROF,JLEV)&
           & +PB2(JROF,MSLB2URL+JLEV-1)*YDGSGEOM%GNORDL(JROF)&
           & +PB2(JROF,MSLB2VRL+JLEV-1)*YDGSGEOM%GNORDM(JROF))  
          ZWF(JROF,JLEV)=0.5_JPRB*(ZWF(JROF,JLEV)+PB2(JROF,MSLB2WRL+JLEV-1))
        ENDIF
        ZNOR25            =(ZPU5*ZPU5+ZPV5*ZPV5)
        ZNOR2             =2.0_JPRB*(ZPU*ZPU5+ZPV*ZPV5)
        PSCO(JROF,JLEV,YDTSCO%M_COPHI) =-ZDTS22*ZNOR2
        ZSPHSV5           =ZDTSA*(1.0_JPRB-ZDTS62*ZNOR25)
        ZSPHSV            =-ZDTSA*ZDTS62*ZNOR2
        ZINT              =ZPV*ZSPHSV5+ZPV5*ZSPHSV
        PSCO(JROF,JLEV,YDTSCO%M_SINLA) = YDGSGEOM%GEMU(JROF)*PSCO(JROF,JLEV,YDTSCO%M_COPHI)-ZINT*YDGSGEOM%GSQM2(JROF)
        PSCO(JROF,JLEV,YDTSCO%M_COSCO) = YDGSGEOM%GSQM2(JROF)*PSCO(JROF,JLEV,YDTSCO%M_COPHI)+ZINT*YDGSGEOM%GEMU(JROF)
        PSCO(JROF,JLEV,YDTSCO%M_SINCO) =-ZPU*ZSPHSV5-ZPU5*ZSPHSV

        ! * computations on a vertical.

        ZLEVO=-RTDT*ZWF(JROF,JLEV)
        IF(LLO.OR.YDML_DYN%YRDYNA%LELTRA) THEN
          ZLEVO5=PLEV5(JROF,JLEV,JITER)
          IF (ZLEVO5 <= ZVETAON) THEN
            ZLEVO=0.0_JPRB
          ELSEIF (ZLEVO5 >= ZVETAOX) THEN
            ZLEVO=0.0_JPRB
          ENDIF
          PLEV(JROF,JLEV)=ZLEVO
        ENDIF

      ENDDO
    ENDDO

  ENDIF

  !*       2.2   DETERMINATION OF THE "(a/rs)*wind" AT "O".

  IF(JITER /= NITMP) THEN

    IF (LSLDP_CURV) THEN
      IROT=0
    ELSE
      IROT=1
    ENDIF
    ITIP=1

    CALL LARCINATL(YDGEOMETRY,YDML_DYN,KST,KPROF,YDSL,KSTABUF,LLSLHD,LLSLHDQUAD,LLSLHD_OLD,ITIP,IROT,.FALSE.,&
     & PLSDEPI,KIBL,&
     & PSCO,PLEV,PSCO5(1,1,1,JITER),PLEV5(1,1,JITER),&
     & ZUN_KAPPA,ZUN_KAPPA5,ZUN_KAPPAT,ZUN_KAPPAT5,&
     & ZSLB1UR0,ZSLB1VR0,ZSLB1ZR0,ZSLB1WR0,&
     & ZSLB1UR05,ZSLB1VR05,ZSLB1ZR05,ZSLB1WR05,&
     & PCCO,ZUF,ZVF,ZZF,ZWF,PCCO5(1,1,1,JITER),&
     & KL0(1,1,0,JITER),KLH0,KLEV,PLSCAW,PRSCAW,&
     & PLSCAW5(1,1,1,JITER),PRSCAW5(1,1,1,JITER),INOWENO)

  ENDIF

ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LARMESTL',1,ZHOOK_HANDLE)
END SUBROUTINE LARMESTL
