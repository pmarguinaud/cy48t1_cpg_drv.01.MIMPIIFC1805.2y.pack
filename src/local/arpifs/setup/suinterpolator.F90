SUBROUTINE SUINTERPOLATOR(YDGEOMETRY)

!**** *SUVERT*  - Routine to initialize vertical interpolator

!     Purpose.
!     --------
!           Initialize vertical interpolator

!**   Interface.
!     ----------

!     *CALL* SUINTERPOLATOR
!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!        see the modules used above.

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Tomas Wilhelmsson from SUVERT.
!      Original : 2013-08-22

!     Modifications.
!     --------------
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!-------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCT0   , ONLY : LRPLANE
USE YOMMP0   , ONLY : NPRINTLEV, LOUTPUT

!-------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
LOGICAL :: LLVERBOSE,LLEQMER
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "suvsleta.intfb.h"
#include "suvsplip.intfb.h"
#include "suhslmer.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUINTERPOLATOR',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV)
ASSOCIATE(NDGENH=>YDDIM%NDGENH, NDGSAH=>YDDIM%NDGSAH, NDGSUR=>YDDIM%NDGSUR,   NFLEVG=>YDDIMV%NFLEVG)
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------

!*            1. SET UP VERTICAL INTERPOLATOR.
!             --------------------------------

IF (NFLEVG > 1) THEN
  ! * externalisable part of interpolator: computes RVSPTRI, RVSPC, RFVV.
  CALL SUVSPLIP(YDGEOMETRY)

  ! * externalisable part of interpolator: attributes of YRVSLETA.
  LLVERBOSE=LOUTPUT.AND.NPRINTLEV >= 1
  CALL SUVSLETA(YDGEOMETRY,LLVERBOSE)
ENDIF

! * externalisable part of interpolator: intermediate quantities for meridian cubic interpolations.
IF (LRPLANE .OR. (NDGSUR>=2 .AND. .NOT.LRPLANE)) THEN
  LLEQMER=LRPLANE
  CALL SUHSLMER(YDGEOMETRY,NDGSAH,NDGENH,LLEQMER)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUINTERPOLATOR',1,ZHOOK_HANDLE)
END SUBROUTINE SUINTERPOLATOR
