SUBROUTINE LARMES5(YDGEOMETRY,YDRIP,YDML_DYN,KST,KPROF,YDSL,KSTABUF,PB15,PB2,&
 & PLSDEPI,KIBL,KDEP,KNOWENO,&
 & PSCO5,PLEV5,PCCO5,PUF5,PVF5,PZF5,KL0,KLH0,KLEV,PLSCAW5,PRSCAW5)

!**** *LARMES5 - semi-LAgrangian scheme:
!                Research of the origin point on the Sphere.
!                Trajectory computation for TL/AD models.

!     Purpose.
!     --------

!      The computation of the location of the interpolation point of
!     the lagrangian trajectory is performed by an iterative
!     method described by Robert and adapted to the sphere by M. Rochas.
!     Trajectories are great circles.
!     (Repeat of computations in LARMES, but useful results are saved
!     at each iteration - used by adjoint of semi-Lagrangian scheme)

!**   Interface.
!     ----------
!        *CALL* *LARMES5(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST       - first element of arrays where computations are performed.
!          KPROF     - depth of work.
!          YDSL      - SL_STRUCT definition.
!          KSTABUF   - for a latitude IGL, KSTABUF(IGL) is the
!                      address of the element corresponding to
!                      (ILON=1,IGL) in the NPROMA arrays.
!          PB15      - SLBUF1-buffer for interpolations (trajectory).
!          PB2       - SLBUF2-buf to communicate info from non lag. to lag. dyn.
!          PLSDEPI   - (Number of points by latitude) / (2 * PI) .
!          KIBL      - index into YRGSGEOM/YRCSGEOM instances in YDGEOMETRY

!        OUTPUT:
!          KDEP      - indication of the interpolation stencil latitudial
!                      dependences (used for LVECADIN=.T. option in adjoint)
!          KNOWENO   - indication of the specific treatment of vertical boundary used
!                       by WENO scheme (only usefull to pass upwards in adjoint code).
!          PSCO5     - information about geographic position of interpol. point.
!          PLEV5     - vertical coordinate of the interpolation point.
!          PCCO5     - information about comput. space position of interpol. point.
!          PUF5,PVF5,PZF5 - U,V&Z-comp of "(a/rs)*wind" necessary to
!                      find the position of the origin point,
!                      in a local repere linked to computational sphere.
!          KL0       - index of the four western points
!                      of the 16 points interpolation grid.
!          KLH0      - second value of index of the four western points
!                      of the 16 points interpolation grid if needed.
!          KLEV      - lower level of the vertical interpolation
!                      grid needed for vertical interpolations.
!          PLSCAW5   - linear weights (distances) for interpolations.
!          PRSCAW5   - high order weights (distances) for interpolations.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------
!        Calls  LARCINA.
!        Is called by LAPINEA5 (3D model) 

!     Reference.
!     ----------

!     Author.
!     -------
!      C. Temperton (ECMWF)
!      Original : 99/10/09

!     Modifications.
!     --------------
!   Modified 01-08-30 by K. YESSAD: pruning and some other cleanings.
!   01-Oct-2003 M. Hamrud     CY28 Cleaning
!   14-Jan-2005 C. Temperton  transfer computations from start of LAIDEP
!                             to end of LARMES.
!   08-Feb-2005 D. Salmond    SLADJTR replaced by LARMES5
!   F.Vana  09-Jan-2007 update of arguments for LARCINA
!   F.Vana  28-Aug-2007 update of arguments for LARCINA
!   30-Jun-2008 J. Masek      Updated arguments of LARCINA.
!   K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!   K. Yessad (Aug 2009): always use root (QX,QY) for (p,q) variables names
!   K. Yessad (Nov 2009): keep LLO.OR.LELTRA=T only.
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   F. Vana 21-Feb-2011: update of weights dimensions (hor. turbulence)
!   G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!   G. Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TGSGEOM and TCSGEOM
!   F. Vana 13-Feb-2014       Updated arguments of LARCINA.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   S. Malardel (Nov 2013): COMAD weights for SL interpolations
!   K. Yessad (July 2014): Move some variables.
!   F. Vana    21-Nov-2017: Options LSLDP_CURV and LHOISLT
!   F. Vana    July 2018 : Beauty content & optimization
!   F. Vana October 2018: Extended LSLDP_CURV.
!   F. Vana    26-Feb-2019: Vertical quintic interpolation
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
! End Modifications
!-------------------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB 
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE SC2PRG_MOD         , ONLY : SC2PRG
USE YOMDYNA            , ONLY : YRDYNA
USE YOMRIP             , ONLY : TRIP
USE EINT_MOD           , ONLY : SL_STRUCT

!-------------------------------------------------------------------------------

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
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTSCO%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLEV5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCCO5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTCCO%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP-1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP-1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP-1)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KL0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3,YDML_DYN%YRDYN%NITMP)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLH0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,0:3)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTLSCAW%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRSCAW5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YYTRSCAW%NDIM, &
 & YDML_DYN%YRDYN%NITMP)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KDEP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KNOWENO(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_DYN%YRDYN%NITMP)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ITIP

INTEGER(KIND=JPIM) :: IROT, JITER, JLEV, JROF, IHVI

REAL(KIND=JPRB)    :: ZPU5,ZLEVO5,ZINT5,ZPW5,ZPV5,ZSPHSV5,&
 & ZDTS62,ZDTSA,ZDTS22,ZNOR25

LOGICAL :: LLFINDVSEP
LOGICAL :: LLINTV
LOGICAL :: LLO
LOGICAL :: LLSLHD, LLSLHDQUAD

REAL(KIND=JPRB) :: ZVETAON, ZVETAOX

REAL(KIND=JPRB) :: ZWF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZU05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZV05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZZ05(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPU05, ZPV05

! unused arguments in call to LARCINA:
! a) input arrays - will not be used since LSLHD was set to .FALSE.
!    in LAPINEA before call to LARMES/ELARMES => dimensions can be
!    contracted
REAL(KIND=JPRB) :: ZUN_KAPPA(1),ZUN_KAPPAT(1),ZUN_KAPPAM(1),ZUN_KAPPAH(1)
REAL(KIND=JPRB) :: ZUN_STDDISU(1),ZUN_STDDISV(1),ZUN_STDDISW(1)
! b) input/output
INTEGER(KIND=JPIM) :: IUN_VSEPC
INTEGER(KIND=JPIM) :: IUN_VSEPL
! c) output arrays - dimensions kept for safety
REAL(KIND=JPRB) :: ZUN_WFSM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB), POINTER :: ZSLB1UR05(:), ZSLB1VR05(:), ZSLB1WR05(:), ZSLB1ZR05(:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "larcina.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LARMES5',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDMP=>YDGEOMETRY%YRMP, YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, &
 & YDVFE=>YDGEOMETRY%YRVFE, YDCSGEOM=>YDGEOMETRY%YRCSGEOM(KIBL),YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL) , &
 & YDDYN=>YDML_DYN%YRDYN, &
 & YDPTRSLB2=>YDML_DYN%YRPTRSLB2,YDPTRSLB15=>YDML_DYN%YRPTRSLB15)

ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & NITMP=>YDDYN%NITMP, VETAON=>YDDYN%VETAON, VETAOX=>YDDYN%VETAOX, &
 & LSLDP_CURV=>YDDYN%LSLDP_CURV, &
 & MSLB1UR05=>YDPTRSLB15%MSLB1UR05, MSLB1VR05=>YDPTRSLB15%MSLB1VR05, &
 & MSLB1ZR05=>YDPTRSLB15%MSLB1ZR05, &
 & MSLB1WR05=>YDPTRSLB15%MSLB1WR05, NFLDSLB15=>YDPTRSLB15%NFLDSLB15, &
 & MSLB2URL5=>YDPTRSLB2%MSLB2URL5, MSLB2VRL5=>YDPTRSLB2%MSLB2VRL5, &
 & MSLB2WRL5=>YDPTRSLB2%MSLB2WRL5, NFLDSLB2=>YDPTRSLB2%NFLDSLB2, &
 & RDTS22=>YDRIP%RDTS22, RDTS62=>YDRIP%RDTS62, RDTSA=>YDRIP%RDTSA, &
 & RTDT=>YDRIP%RTDT)

CALL SC2PRG(MSLB1UR05 ,PB15     ,ZSLB1UR05)
CALL SC2PRG(MSLB1VR05 ,PB15     ,ZSLB1VR05)
CALL SC2PRG(MSLB1WR05 ,PB15     ,ZSLB1WR05)
CALL SC2PRG(MSLB1ZR05 ,PB15     ,ZSLB1ZR05)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS AND TESTS.
!              --------------------------------------

!*       1.2   Miscellaneous preliminary initialisations.

! in practical LLO.OR.LELTRA is now always T
LLO=.NOT.YRDYNA%LELTRA

IF(LLO.OR.YRDYNA%LELTRA) THEN
  ZDTS22=RDTS22*4._JPRB
  ZDTSA=RDTSA*2._JPRB
  ZDTS62=RDTS62*4._JPRB
ELSE
  CALL ABOR1(' LARMES5: LLO.OR.LELTRA=F no longer supported')
ENDIF

! deactivate computation of SLHD weights
LLSLHD     =.FALSE.
LLSLHDQUAD =.FALSE.

!     ------------------------------------------------------------------

!*       2.    ITERATIONS.
!              -----------

ZVETAON=(1.0_JPRB-VETAON)*YDVETA%VETAH(0)+VETAON*YDVETA%VETAF(1)
ZVETAOX=(1.0_JPRB-VETAOX)*YDVETA%VETAH(NFLEVG)+VETAOX*YDVETA%VETAF(NFLEVG)

DO JITER=1,NITMP

  !*       2.1   DETERMINATION OF THE ORIGIN POINT "O".

  ! Computation of the norm of the real wind vector.
  ! Computation of the angle (PHI=DT . NORM OF V ON RADIUS)**2
  ! then computation of the coordinates of the origin point "O".

  IF (JITER == 1) THEN

    DO JLEV=1,NFLEVG
      DO JROF=KST,KPROF

        ! * computations on horizontal plans.

        ! When relevant, convert arrival point wind to cartesian space
        IF (LSLDP_CURV.AND.(.NOT.YRDYNA%LELTRA)) THEN  
          ! Arrival point wind
          ZPU05= +PB2(JROF,MSLB2URL5+JLEV-1)*YDGSGEOM%GNORDM(JROF)&
           &     -PB2(JROF,MSLB2VRL5+JLEV-1)*YDGSGEOM%GNORDL(JROF)  
          ZPV05= +PB2(JROF,MSLB2URL5+JLEV-1)*YDGSGEOM%GNORDL(JROF)&
           &     +PB2(JROF,MSLB2VRL5+JLEV-1)*YDGSGEOM%GNORDM(JROF)
          ! Arrival point wind converted to Cartesian space
          ZU05(JROF,JLEV)=-ZPU05*YDGSGEOM%GESLO(JROF)             &
           &              -ZPV05*YDGSGEOM%GECLO(JROF)*YDGSGEOM%GEMU(JROF)
          ZV05(JROF,JLEV)= ZPU05*YDGSGEOM%GECLO(JROF)             &
           &              -ZPV05*YDGSGEOM%GESLO(JROF)*YDGSGEOM%GEMU(JROF)
          ZZ05(JROF,JLEV)= ZPV05*YDGSGEOM%GSQM2(JROF)
        ENDIF

        ZNOR25=( PB2(JROF,MSLB2URL5+JLEV-1)*PB2(JROF,MSLB2URL5+JLEV-1)&
         & +PB2(JROF,MSLB2VRL5+JLEV-1)*PB2(JROF,MSLB2VRL5+JLEV-1) )  
        PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI,1)=1.0_JPRB-ZDTS22*ZNOR25
        ZSPHSV5=ZDTSA*(1.0_JPRB-ZDTS62*ZNOR25)
        ZINT5=( PB2(JROF,MSLB2URL5+JLEV-1)*YDGSGEOM%GNORDL(JROF)&
         & +PB2(JROF,MSLB2VRL5+JLEV-1)*YDGSGEOM%GNORDM(JROF))*ZSPHSV5  
        PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_SINLA,1)=&
         & YDGSGEOM%GEMU(JROF)*PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI,1)-ZINT5*YDGSGEOM%GSQM2(JROF)  
        PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO,1)=&
         & YDGSGEOM%GSQM2(JROF)*PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI,1)+ZINT5*YDGSGEOM%GEMU(JROF)  
        PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO,1)=-( PB2(JROF,MSLB2URL5+JLEV-1)*YDGSGEOM%GNORDM(JROF)&
         & -PB2(JROF,MSLB2VRL5+JLEV-1)*YDGSGEOM%GNORDL(JROF) )*ZSPHSV5  

        ! * computations on a vertical.

        ZLEVO5=YDVETA%VETAF(JLEV)-RTDT*PB2(JROF,MSLB2WRL5+JLEV-1)
        ZLEVO5=MIN(ZVETAOX,MAX(ZVETAON,ZLEVO5))
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

        IF(LSLDP_CURV.AND.YRDYNA%LELTRA) THEN
          ! Convert Departure wind back to lat,lon at Arrival point
          !  (trajectory values must remain in Cartesian space)
          ZPU5 = -PUF5(JROF,JLEV,JITER-1)*YDGSGEOM%GESLO(JROF)                      &
           &     +PVF5(JROF,JLEV,JITER-1)*YDGSGEOM%GECLO(JROF)
          ZPV5 = -PUF5(JROF,JLEV,JITER-1)*YDGSGEOM%GECLO(JROF)*YDGSGEOM%GEMU(JROF)  &
           &     -PVF5(JROF,JLEV,JITER-1)*YDGSGEOM%GESLO(JROF)*YDGSGEOM%GEMU(JROF)  &
           &     +PZF5(JROF,JLEV,JITER-1)*YDGSGEOM%GSQM2(JROF)
        ELSEIF(YRDYNA%LELTRA) THEN
          ZPU5=PCCO5(JROF,JLEV,YDML_DYN%YYTCCO%M_RQX,JITER-1)*PUF5(JROF,JLEV,JITER-1)&
           & +PCCO5(JROF,JLEV,YDML_DYN%YYTCCO%M_RQY,JITER-1)*PVF5(JROF,JLEV,JITER-1)  
          ZPV5=-PCCO5(JROF,JLEV,YDML_DYN%YYTCCO%M_RQY,JITER-1)*PUF5(JROF,JLEV,JITER-1)&
           & +PCCO5(JROF,JLEV,YDML_DYN%YYTCCO%M_RQX,JITER-1)*PVF5(JROF,JLEV,JITER-1)  
          ZPW5=ZWF5(JROF,JLEV)
        ELSEIF(LSLDP_CURV.AND.LLO) THEN
          ! Averaging wind to Midpoint in Cartesian space
          PUF5(JROF,JLEV,JITER-1)=0.5_JPRB*(PUF5(JROF,JLEV,JITER-1)+ZU05(JROF,JLEV))
          PVF5(JROF,JLEV,JITER-1)=0.5_JPRB*(PVF5(JROF,JLEV,JITER-1)+ZV05(JROF,JLEV))
          PZF5(JROF,JLEV,JITER-1)=0.5_JPRB*(PZF5(JROF,JLEV,JITER-1)+ZZ05(JROF,JLEV))

          ! Convert Midpoint wind back to lat,lon at Arrival point
          !  (trajectory values must remain in Cartesian space)
          ZPU5 = -PUF5(JROF,JLEV,JITER-1)*YDGSGEOM%GESLO(JROF)                      &
           &     +PVF5(JROF,JLEV,JITER-1)*YDGSGEOM%GECLO(JROF)
          ZPV5 = -PUF5(JROF,JLEV,JITER-1)*YDGSGEOM%GECLO(JROF)*YDGSGEOM%GEMU(JROF)  &
           &     -PVF5(JROF,JLEV,JITER-1)*YDGSGEOM%GESLO(JROF)*YDGSGEOM%GEMU(JROF)  &
           &     +PZF5(JROF,JLEV,JITER-1)*YDGSGEOM%GSQM2(JROF)
          ZPW5=0.5_JPRB*(ZWF5(JROF,JLEV)+PB2(JROF,MSLB2WRL5+JLEV-1))
        ELSEIF(LLO) THEN
          ZPU5=0.5_JPRB*(PCCO5(JROF,JLEV,YDML_DYN%YYTCCO%M_RQX,JITER-1)*PUF5(JROF,JLEV,JITER-1)&
           & +PCCO5(JROF,JLEV,YDML_DYN%YYTCCO%M_RQY,JITER-1)*PVF5(JROF,JLEV,JITER-1)&
           & +PB2(JROF,MSLB2URL5+JLEV-1)*YDGSGEOM%GNORDM(JROF)&
           & -PB2(JROF,MSLB2VRL5+JLEV-1)*YDGSGEOM%GNORDL(JROF))  
          ZPV5=0.5_JPRB*(-PCCO5(JROF,JLEV,YDML_DYN%YYTCCO%M_RQY,JITER-1)*PUF5(JROF,JLEV,JITER-1)&
           & +PCCO5(JROF,JLEV,YDML_DYN%YYTCCO%M_RQX,JITER-1)*PVF5(JROF,JLEV,JITER-1)&
           & +PB2(JROF,MSLB2URL5+JLEV-1)*YDGSGEOM%GNORDL(JROF)&
           & +PB2(JROF,MSLB2VRL5+JLEV-1)*YDGSGEOM%GNORDM(JROF))  
          ZPW5=0.5_JPRB*(ZWF5(JROF,JLEV)+PB2(JROF,MSLB2WRL5+JLEV-1))
        ENDIF
        ZNOR25            =(ZPU5*ZPU5+ZPV5*ZPV5)
        PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI,JITER) =1.0_JPRB-ZDTS22*ZNOR25
        ZSPHSV5            =ZDTSA*(1.0_JPRB-ZDTS62*ZNOR25)
        ZINT5             =ZPV5*ZSPHSV5
        PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_SINLA,JITER) =&
         & YDGSGEOM%GEMU(JROF)*PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI,JITER)-ZINT5*YDGSGEOM%GSQM2(JROF)  
        PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COSCO,JITER) =&
         & YDGSGEOM%GSQM2(JROF)*PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_COPHI,JITER)+ZINT5*YDGSGEOM%GEMU(JROF)  
        PSCO5(JROF,JLEV,YDML_DYN%YYTSCO%M_SINCO,JITER) =-ZPU5*ZSPHSV5

        ! * computations on a vertical.

        ZLEVO5=YDVETA%VETAF(JLEV)-RTDT*ZPW5
        ZLEVO5=MIN(ZVETAOX,MAX(ZVETAON,ZLEVO5))
        IF(LLO.OR.YRDYNA%LELTRA) THEN
          PLEV5(JROF,JLEV,JITER)=ZLEVO5
        ENDIF

      ENDDO
    ENDDO

  ENDIF

  !*       2.2   DETERMINATION OF THE "(a/rs)*wind" AT "O".

  IF(JITER /= NITMP) THEN

    IHVI=0
    IF (LSLDP_CURV) THEN
      IROT=0
    ELSE
      IROT=1
    ENDIF
    ITIP=1
    LLFINDVSEP=.FALSE.
    LLINTV=.TRUE.

    CALL LARCINA(&
     & YDGEOMETRY,YDML_DYN,KST,KPROF,YDSL,IHVI,KSTABUF,LLFINDVSEP,LLSLHD,LLSLHDQUAD,LLINTV,&
     & ITIP,IROT,.FALSE.,PLSDEPI,&
     & KIBL,&
     & PSCO5(:,:,:,JITER),PLEV5(:,:,JITER),&
     & ZUN_KAPPA,ZUN_KAPPAT,ZUN_KAPPAM,ZUN_KAPPAH,&
     & ZUN_STDDISU,ZUN_STDDISV,ZUN_STDDISW,&
     & ZSLB1UR05,ZSLB1VR05,ZSLB1ZR05,ZSLB1WR05,&
     & IUN_VSEPC,IUN_VSEPL,PCCO5(:,:,:,JITER),&
     & PUF5(:,:,JITER),PVF5(:,:,JITER),PZF5(:,:,JITER),ZWF5,ZUN_WFSM,&
     & KL0(:,:,:,JITER),KLH0,KLEV,PLSCAW5(:,:,:,JITER),PRSCAW5(:,:,:,JITER),&
     & KDEP(:,:,JITER),KNOWENO(:,:,JITER))

  ENDIF

ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LARMES5',1,ZHOOK_HANDLE)
END SUBROUTINE LARMES5
