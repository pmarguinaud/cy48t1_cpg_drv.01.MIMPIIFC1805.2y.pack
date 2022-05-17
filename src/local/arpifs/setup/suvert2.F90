SUBROUTINE SUVERT2(YDGEOMETRY)

!**** *SUVERT2*  - Routine to initialize vertical coordinate
!                  Simplified version of SUVERT for 2D models and more generally
!                  for all configurations using NIOLEVG=NFLEVG=1.

!     Purpose.
!     --------
!           Initialize the hybrid-cordinate system of the model in 2D models.

!**   Interface.
!     ----------

!     *CALL* SUVERT2

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
!      K. Yessad (July 2012) after SUVERT

!     Modifications.
!     --------------
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
! End Modifications
!-------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK
USE YOMVV1   , ONLY : DVALH, DVBH

!-------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVERT2',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,   NIOLEVG=>YDGEOMETRY%YRDIMV%NIOLEVG)!     ------------------------------------------------------------------

IF (NIOLEVG /= 1 .OR. NFLEVG /= 1) THEN
  CALL ABOR1('SUVERT2 may be called only if NIOLEVG=1 and NFLEVG=1.')
ENDIF

! * fill DVALH and DVBH.
DVALH(0)=0.0_JPRB
DVALH(1)=0.0_JPRB
DVBH(0)=0.0_JPRB
DVBH(1)=1.0_JPRB

! * fill YRVETA.
YDVETA%VETAF(0)=0._JPRB
YDVETA%VETAF(1)=0.5_JPRB
YDVETA%VETAF(2)=1.0_JPRB
YDVETA%VETAH(0)=0._JPRB
YDVETA%VETAH(1)=1.0_JPRB

! * fill YRVAB (attributes VC, VDELA, VDELB, VAF, VBF are not used in this case).
YDVAB%VAH(0)=0._JPRB
YDVAB%VAH(1)=0._JPRB
YDVAB%VALH(0)=0._JPRB
YDVAB%VALH(1)=0._JPRB
YDVAB%VBH(0)=0._JPRB
YDVAB%VBH(1)=1._JPRB
YDVAB%VRATH(0)=0.0_JPRB
YDVAB%VRATH(1)=1.0_JPRB
YDVAB%VRATF(1)=0.5_JPRB

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVERT2',1,ZHOOK_HANDLE)
END SUBROUTINE SUVERT2
