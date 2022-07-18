SUBROUTINE TRANSINVHAD(YDGEOMETRY,YDGFL,YDGMV,YDML_GCONF,CDCONF,YDSP)

!**** *TRANSINVHAD * - Inverse transforms - adjoint

!     Purpose.  Perform inverse transform (spectral to gridpoint)
!     --------

!     Explicit arguments :
!     --------------------
!        CDCONF     - configuration of work
!                     Current values in use for CDCONF(2:3): AA, GB

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-10-25
!        Y.Tremolet: 01-02-28 cleanup
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGFL   , ONLY : TGFL
USE YOMGMV   , ONLY : TGMV
USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMSP    , ONLY : SPVOR_FLT, SPDIV_FLT  
USE YOMDYNA   ,ONLY : YRDYNA
!USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)
USE SPECTRAL_FIELDS_MOD

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL)          ,INTENT(INOUT) :: YDGFL
TYPE(TGMV)          ,INTENT(INOUT) :: YDGMV
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
CHARACTER(LEN=9)    ,INTENT(IN)    :: CDCONF
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSP
!     ------------------------------------------------------------------
LOGICAL :: LLDERR,LLMODE,LLFSCOMP,LLVOR
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "transinv_mdlad.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TRANSINVHAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDDIMF=>YDML_GCONF%YRDIMF)
ASSOCIATE(LADER=>YDDIMF%LADER, LVOR=>YDDIMF%LVOR, &
 & GFL=>YDGFL%GFL, &
 & GMV=>YDGMV%GMV, GMVS=>YDGMV%GMVS)
!     ------------------------------------------------------------------

LLDERR = LLE(CDCONF(2:2),'F') .AND. LADER 
LLMODE = (CDCONF(2:3)=='AA').OR.(CDCONF(2:3)=='GB')
LLFSCOMP = .FALSE.
LLVOR=LVOR.AND.LLDERR

IF(LLMODE) THEN
  IF (YRDYNA%LGRADSP) THEN
    CALL TRANSINV_MDLAD(YDGEOMETRY,YDGMV,YDML_GCONF,YDSP%VOR,YDSP%DIV,YDSP%SP,YDSP%HV,YDSP%GFL,&
     & GMV,GMVS,GFL,LLDERR,LLVOR,LLFSCOMP,PSPVOR_FLT=SPVOR_FLT,PSPDIV_FLT=SPDIV_FLT)
  ELSE
    CALL TRANSINV_MDLAD(YDGEOMETRY,YDGMV,YDML_GCONF,YDSP%VOR,YDSP%DIV,YDSP%SP,YDSP%HV,YDSP%GFL,&
     & GMV,GMVS,GFL,LLDERR,LLVOR,LLFSCOMP)
  ENDIF
ELSE
  CALL ABOR1(' TRANSINVHAD - CONFIGURATION NOT SUPPORTED ')
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TRANSINVHAD',1,ZHOOK_HANDLE)
END SUBROUTINE TRANSINVHAD
