SUBROUTINE STEPOTL(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,CDCONF,PII0,YDGOM5,YDGOM,YDACV)

!**** *STEPOTL*  - Controls tangent linear integration at lowest level

!     Purpose.
!     --------
!     Controls the tangent linear integration.

!**   Interface.
!     ----------
!        *CALL* *STEPOTL(CDCONF)

!        Explicit arguments :
!        --------------------
!     CHARACTER*9 CDCONF                (in)
!   1   : configuration of IOPACK        IO handling
!   2+3 : configuration of (E)TRANSINVH  inverse transforms
!   4   : configuration of SCAN2MTL      grid point computations   (under SCAN2M*)
!   6   : configuration of OBS           comparison to observations(under SCAN2M*)
!   7   : configuration of ECOUPL1       coupling for LAM models
!   8   : configuration of (E)TRANSDIRH  direct transforms
!   9   : configuration of (E)SPC..      spectral space computations

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.    see includes below.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!       Mats Hamrud and Philippe Courtier  *ECMWF*
!       Original    : 87-10-15

! Modifications
! -------------
!   J.Vivoda    : 03-2002  PC schemes for NH dynamics (LPC_XXXX keys)
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!   10-Jun-2004 J. Masek   NH cleaning (LPC_NOTR)
!   01-Jul-2008 :B. Chapnik LL_CPL consistent between stepotl and ad
!   F. Vana     : 13-Jan-2009 Remove useless LSPC_FROM_DI key
!   K. Yessad (Sep 2010): organigramme simplification.
!   K. Yessad (Jan 2011): new architecture for LBC modules and set-up.
!   K. Yessad (Dec 2011): various contributions.
!   F. Vana  28-Nov-2013: Redesigned trajectory handling
!   K. Yessad (July 2014): Move some variables.
!   F. Vana  05-Mar-2015  Support for single precision
!      R. El khatib 22-Feb-2016 NSTEP passed to opdis
!   O. Marsden  Aug 2016  Removed explicit YDSPEC argument, using that in YDFIELDS instead
!   K. Yessad (June 2017): Introduce NHQE model.
!   S. Massart 19-Feb-2019 : Solar constant optimisation
! End Modifications
!------------------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE FIELDS_MOD         , ONLY : FIELDS
USE MTRAJ_MOD          , ONLY : MTRAJ
USE PARKIND1           , ONLY : JPRD, JPIM, JPRB
USE YOMCT0             , ONLY : LR2D, NCONF, LOPDIS, LELAM
USE YOMDYNA            , ONLY : LRUBC
USE YOMCT3             , ONLY : NSTEP
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE YOMLUN             , ONLY : NULOUT
USE YOMTIM             , ONLY : RSTART, RVSTART, RTIMEF
USE YEMLBC_INIT         , ONLY : NECRIPL, NECOTL, LUNBC
USE YOMTRAJ            , ONLY : TRAJEC
USE SUPERGOM_CLASS     , ONLY : CLASS_SUPERGOM
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)
USE TYPE_ACV           , ONLY : ACV_CONTAINER

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      , INTENT(INOUT) :: YDGEOMETRY !! INOUT needed for call to IOPACK
TYPE(FIELDS)        , INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ)         , INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL)         , INTENT(INOUT) :: YDMODEL
CHARACTER(LEN=9)    , INTENT(IN)    :: CDCONF 
REAL(KIND=JPRB)     , INTENT(IN), OPTIONAL    :: PII0
TYPE(CLASS_SUPERGOM), INTENT(INOUT), OPTIONAL :: YDGOM5
TYPE(CLASS_SUPERGOM), INTENT(INOUT), OPTIONAL :: YDGOM
TYPE(ACV_CONTAINER) , INTENT(INOUT), OPTIONAL :: YDACV
!     ------------------------------------------------------------------
CHARACTER(LEN=1)   :: CLCONF_TRANSDIR
INTEGER(KIND=JPIM) :: ISPEC2, ISPEC2V
INTEGER(KIND=JPIM) :: IJUM

LOGICAL :: LLESPCL      ! do spectral nudging in LAM models
LOGICAL :: LLCPL        ! do coupling in LAM models
LOGICAL :: LL_TST_GPGFL ! do timestepping on grid-point GFL
LOGICAL :: LL_DFISTEP
LOGICAL :: LLIAU        ! do IAU in LAM models

REAL(KIND=JPRB) :: ZJUM, ZII0
REAL(KIND=JPRD) :: ZCT, ZVT, ZWT
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "user_clock.h"

#include "abor1.intfb.h"
#include "ecoupl1.intfb.h"
#include "ecoupl2.intfb.h"
#include "espcm.intfb.h"
#include "etransdirh.intfb.h"
#include "etransinvh.intfb.h"
#include "iopack.intfb.h"
#include "opdis.intfb.h"
#include "scan2mtl.intfb.h"
#include "spcm.intfb.h"
#include "spc2m.intfb.h"
#include "transdirh.intfb.h"
#include "transinvh.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('STEPOTL',0,ZHOOK_HANDLE)
ASSOCIATE(YDGFL5=>YDMTRAJ%YRGFL5, &
 & YDGMV5=>YDMTRAJ%YRGMV5, YDSPEC=>YDFIELDS%YRSPEC, &
 & YDGFL=>YDFIELDS%YRGFL,YDGMV=>YDFIELDS%YRGMV, YDSURF=>YDFIELDS%YRSURF, YDDIM=>YDGEOMETRY%YRDIM, &
 & YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, &
 & YDMP=>YDGEOMETRY%YRMP, YDLAP=>YDGEOMETRY%YRLAP, YDRIP=>YDMODEL%YRML_GCONF%YRRIP)

ASSOCIATE(NSTART=>YDRIP%NSTART, NSTOP=>YDRIP%NSTOP)
!     ------------------------------------------------------------------

!*       0.    VARIOUS INITIALIZATIONS
!              -----------------------

! Timestepping on grid-point GFL is done when both:
! - grid-point calculations are required
! - direct spectral transforms are required
LL_TST_GPGFL=(CDCONF(7:8) /= '00')

CALL GSTATS(1889,0)
CALL USER_CLOCK(PELAPSED_TIME=ZWT,PVECTOR_CP=ZVT,PTOTAL_CP=ZCT)
ZCT=ZCT-RSTART
ZVT=ZVT-RVSTART
ZWT=ZWT-RTIMEF
RSTART=RSTART+ZCT
RVSTART=RVSTART+ZVT
RTIMEF=RTIMEF+ZWT
ZJUM=10._JPRB**(INT(LOG10(REAL(MAX(NSTOP-NSTART,1),JPRB)))+1)
IJUM=NINT(MAX(ZJUM/100._JPRB,1.0_JPRB))
IF(NSTEP-NSTART <= 10.OR.MOD(NSTEP,IJUM) == 0)THEN
  WRITE(NULOUT,'('' NSTEP ='',I6,'' STEPOTL '',A9)') NSTEP,CDCONF
ENDIF
CALL GSTATS(1889,1)

IF (LOPDIS) THEN
  CALL OPDIS(CDCONF,'STEPOTL',ZCT,ZVT,ZWT,RSTART,RTIMEF,NSTEP)
ENDIF

!        0.1   Define configurations for coupling LBC in LAM

IF (LELAM) THEN
  LLCPL = (YDMODEL%YRML_LBC%LECOBI.AND.NECRIPL == 1.AND..NOT.LRUBC) .AND.((&
   & CDCONF(3:3) == 'A'.OR.CDCONF(3:3) == 'B'.OR.&
   & CDCONF(6:6) == '1')&
   & .AND.CDCONF(5:5) == '0').AND.(&
   & NCONF/100 == 0.OR.NCONF/100 == 2.OR.&
   & NCONF/100 == 7.OR.NCONF == 801.OR. NECOTL < 0)
ELSE
  LLCPL=.FALSE.
ENDIF

IF (CDCONF(7:7)=='D') THEN
  LL_DFISTEP=.TRUE.
ELSE
  LL_DFISTEP=.FALSE.
ENDIF

!     ------------------------------------------------------------------

!*       1.    READ TRAJECTORY (NOT NORMALLY DONE HERE ANY LONGER).
!              ----------------------------------------------------

IF(CDCONF(1:1) /= '0')THEN
  CALL IOPACK(YDGEOMETRY,YDFIELDS,YDGFL5,YDMODEL,CDCONF(1:1),PTRAJEC=TRAJEC(NSTEP))
ENDIF
!     ------------------------------------------------------------------

!*       2.    INVERSE TRANSFORMS.
!              -------------------

IF(CDCONF(2:2) /= '0')THEN

  IF (LELAM) THEN
    CALL ETRANSINVH(YDGEOMETRY,YDFIELDS%YRGFL,YDFIELDS%YRGMV,YDMODEL%YRML_DIAG,YDMODEL%YRML_GCONF,CDCONF, YDSPEC,YDMTRAJ)
  ELSE
    CALL TRANSINVH(YDGEOMETRY,YDFIELDS%YRGFL,YDFIELDS%YRGMV,YDMODEL%YRML_DIAG,YDMODEL%YRML_GCONF, &
 &                 YDMODEL%YRML_DYN%YRDYN%LRFRIC,YDMODEL%YRML_DYN%YRDYN%RKRF,CDCONF,YDSPEC,YDMTRAJ=YDMTRAJ)
  ENDIF

ENDIF
!     ------------------------------------------------------------------

!*       3.    GRIDPOINT COMPUTATIONS.
!             ------------------------

IF(CDCONF(3:8) /= '000000') THEN
  ZII0 = 0.0_JPRB
  IF (PRESENT(PII0)) ZII0 = PII0
  CALL SCAN2MTL(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDGOM5,YDGOM,CDCONF,&
 &              LL_TST_GPGFL,LL_DFISTEP,TRAJEC(NSTEP),ZII0,YDACV=YDACV)
ENDIF

!     ------------------------------------------------------------------

!*       4.    JO COMPUTATION
!              --------------

!Moved to CNT4TL
!     ------------------------------------------------------------------

!*       8.    DIRECT TRANSFORMS.
!              ------------------

IF(CDCONF(7:8) /= '00')THEN
  CLCONF_TRANSDIR=CDCONF(8:8)
  IF (LELAM) THEN
    IF (LLCPL) THEN
      CALL ECOUPL1(YDGEOMETRY,YDMODEL,YDFIELDS,LL_DFISTEP)
      IF(LUNBC) CALL ECOUPL2(YDGEOMETRY,YDMODEL,YDFIELDS,LL_DFISTEP)
    ENDIF
    CALL ETRANSDIRH(YDGEOMETRY,YDFIELDS%YRGFL,YDFIELDS%YRGMV,CLCONF_TRANSDIR,YDMODEL%YRML_GCONF%YRDIMF%NFTHER,YDSPEC)
  ELSE
    CALL TRANSDIRH(YDGEOMETRY,YDFIELDS%YRGFL,YDFIELDS%YRGMV,CLCONF_TRANSDIR,YDMODEL%YRML_GCONF%YRDIMF%NFTHER,YDSPEC)
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       9.    SPECTRAL COMPUTATIONS.
!              ----------------------

IF(CDCONF(9:9) /= '0')THEN

  IF (LELAM) THEN
    IF (LR2D) THEN
      CALL ABOR1('STEPOTL : ESPC2M DOES NOT EXIST')
    ELSE
      ! no spectral nudging in TL code for the time being
      LLESPCL=.FALSE.
      ! no IAU in TL code 
      LLIAU=.FALSE.
      ISPEC2=0
      ISPEC2V=0
      CALL ESPCM(YDGEOMETRY,YDMODEL,YDFIELDS,CDCONF(9:9),LLESPCL,LLIAU,ISPEC2,ISPEC2V)
    ENDIF
  ELSE
    IF (LR2D) THEN
      CALL SPC2M(YDGEOMETRY,YDRIP,YDMODEL%YRML_DYN%YRDYN,CDCONF(9:9),YDSPEC)
    ELSE
      CALL SPCM(YDGEOMETRY,YDMODEL,CDCONF(9:9),YDSPEC,YDGMV)
    ENDIF
  ENDIF

ENDIF
!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('STEPOTL',1,ZHOOK_HANDLE)
END SUBROUTINE STEPOTL
