SUBROUTINE SUVERTFE(YDGEOMETRY)

!**** *SUVERTFE*  - Setup VERTical Finite Element scheme

!     Purpose.
!     --------
!           Call of initialisation routines for the different 
!           versions of the vertical finite element scheme

!**   Interface.
!     ----------

!     *CALL* SUVERTFE

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!        see below

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mariano Hortal ECMWF
!      Original : 2000-05

!     Modifications.
!     --------------
!      K.Yessad and J.Vivoda: 28-08-2007 Set-up VFE for NH model.
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      J. Vivoda (Oct 2013): new options for VFE-NH
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!      P.Smolikova (Sep 2020): VFE pruning.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPRB ,JPIM
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCT0   , ONLY : LNHDYN
USE YOMCVER  , ONLY : LVERTFE, LVFE_ECMWF, LVFE_GW, LVFE_GW_HALF, NVFE_TYPE
USE YOMDYNA  , ONLY : LGWADV
USE YOMLUN   , ONLY : NULERR

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
REAL(KIND=JPRB)    :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "sunh_vertfe1d.intfb.h"
#include "sunh_vertfe1dd.intfb.h"
#include "sunh_vertfe3d.intfb.h"
#include "sunh_vertfe3dbc.intfb.h"
#include "sunh_vertfe3dd.intfb.h"
#include "suvertfe1.intfb.h"
#include "suvertfe3.intfb.h"
#include "suvertfe3d.intfb.h"
#include "sunh_vertfespline.intfb.h"
#include "sunh_vertfespline_half.intfb.h"
#include "sunh_vertfespline_inv.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUVERTFE',0,ZHOOK_HANDLE)
ASSOCIATE(YDVFE=>YDGEOMETRY%YRVFE)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)

IF(LVERTFE) THEN
  IF (LVFE_ECMWF) THEN

    IF(NVFE_TYPE == 1) THEN

      ! setup linear finite element scheme
      CALL SUVERTFE1(YDGEOMETRY)
      IF (LNHDYN) THEN
        CALL SUNH_VERTFE1D(YDGEOMETRY)
        ! ----- SUNH_VERTFE1DBC not yet coded ----------------------
        ! CALL SUNH_VERTFE1DBC
        ! Provisionally RDERB=RDERI with additional zeroed columns.
        YDVFE%RDERB(1:NFLEVG,1)=0.0_JPRB
        YDVFE%RDERB(1:NFLEVG,2:NFLEVG+1)=YDVFE%RDERI(1:NFLEVG,1:NFLEVG)
        YDVFE%RDERB(1:NFLEVG,NFLEVG+2)=0.0_JPRB
        ! ----------------------------------------------------------
        CALL SUNH_VERTFE1DD(YDGEOMETRY)
      ENDIF

    ELSEIF (NVFE_TYPE == 3) THEN

      ! setup cubic finite element scheme
      CALL SUVERTFE3(YDGEOMETRY)
      IF (LNHDYN) THEN
        CALL SUNH_VERTFE3D(YDGEOMETRY)
        CALL SUNH_VERTFE3DBC(YDGEOMETRY)
        CALL SUNH_VERTFE3DD(YDGEOMETRY)
      ELSE
        CALL SUVERTFE3D(YDGEOMETRY)
      ENDIF
    ELSE
      WRITE(NULERR,*) ' SUVERTFE: INVALID VALUE OF NVFE_TYPE : ',NVFE_TYPE
      WRITE(NULERR,*) ' FOR LVFE_ECMWF.'
    ENDIF

  ELSE

    ! setup B-spline finite element scheme
    CALL SUNH_VERTFESPLINE(YDGEOMETRY)
    IF (LNHDYN) THEN
      CALL SUNH_VERTFESPLINE_HALF(YDGEOMETRY)
      IF (LGWADV.AND.(LVFE_GW.OR.LVFE_GW_HALF)) THEN
        CALL SUNH_VERTFESPLINE_INV(YDGEOMETRY)
      ENDIF
    ENDIF

  ENDIF

ELSE
  CALL ABOR1(' IN SUVERTFE LVERTFE SHOULD BE TRUE')
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVERTFE',1,ZHOOK_HANDLE)
END SUBROUTINE SUVERTFE
