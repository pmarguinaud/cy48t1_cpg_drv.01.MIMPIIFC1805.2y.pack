SUBROUTINE GPMPFC5(YDGMV5,YDML_GCONF,KPROMA,KFLEV,KST,KEN,KFLAG,PGM,PGMV5,PGMV5S,PGFL5)

!**** *GPMPFC5* - Apply map factor to convert 
!                 reduced variables -> geographical variables if kflag=0
!                 or
!                 geographical variables -> reduced variables if kflag=1
!                 Version of GPMPFC for trajectory.

!     Purpose.
!     --------
!           Multiply or divide by map factor.

!**   Interface.
!     ----------
!        *CALL* *GPMPFC5(...)

!        Explicit arguments :
!        --------------------

!      INPUT:
!      ------
!       KPROMA    - horizontal dimensioning
!       KFLEV     - number of layers
!       KST       - start of work
!       KEN       - depth of work
!       KFLAG     - 0 -> multiply, 1-> divide
!       PGM       - map factor

!      INPUT/OUTPUT:
!      -------------
!       PGMV5     - GMV variables
!       PGMV5S    - GMVS variables
!       PGFL5     - GFL variables

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      K. Yessad (June 2011) after GPMPFC and GPMPFC_GMVS

! Modifications
! -------------
! End Modifications
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMGMV                 , ONLY : TGMV
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK
USE YOMCT0                 , ONLY : LNHDYN
USE YOMDYNA                , ONLY : YRDYNA

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV5
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLAG 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGM(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV5(KPROMA,KFLEV,YDGMV5%YT5%NDIM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV5S(KPROMA,YDGMV5%YT5%NDIMS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL5(KPROMA,KFLEV,YDML_GCONF%YGFL%NDIM5) 
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZGM(KPROMA)

INTEGER(KIND=JPIM) :: JLEV,JGFL,JROF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPMPFC5',0,ZHOOK_HANDLE)
ASSOCIATE(NDIM5=>YDML_GCONF%YGFL%NDIM5, NUMFLDS=>YDML_GCONF%YGFL%NUMFLDS, YCOMP=>YDML_GCONF%YGFL%YCOMP, &
 & LUVDER=>YDML_GCONF%YRDIMF%LUVDER, LVOR=>YDML_GCONF%YRDIMF%LVOR, &
 & YT5=>YDGMV5%YT5)
!     ------------------------------------------------------------------

!*       1. APPLY MAP FACTOR.
!           -----------------

IF(KFLAG == 0) THEN
  ZGM(KST:KEN) = PGM(KST:KEN)
ELSEIF(KFLAG == 1) THEN
  ZGM(KST:KEN) = 1.0_JPRB/PGM(KST:KEN)
ELSE
  CALL ABOR1('GPMPFC5: ILLEGAL KFLAG')
ENDIF

! * GMV5:

DO JLEV=1,KFLEV
  DO JROF=KST,KEN
    PGMV5(JROF,JLEV,YT5%MU)=PGMV5(JROF,JLEV,YT5%MU)*ZGM(JROF)
    PGMV5(JROF,JLEV,YT5%MV)=PGMV5(JROF,JLEV,YT5%MV)*ZGM(JROF)
    PGMV5(JROF,JLEV,YT5%MDIV)=PGMV5(JROF,JLEV,YT5%MDIV)*ZGM(JROF)*ZGM(JROF)
    PGMV5(JROF,JLEV,YT5%MTL)=PGMV5(JROF,JLEV,YT5%MTL)*ZGM(JROF)
    PGMV5(JROF,JLEV,YT5%MTM)=PGMV5(JROF,JLEV,YT5%MTM)*ZGM(JROF)
    IF(LUVDER) THEN
      PGMV5(JROF,JLEV,YT5%MUL)=PGMV5(JROF,JLEV,YT5%MUL)*ZGM(JROF)*ZGM(JROF)
      PGMV5(JROF,JLEV,YT5%MVL)=PGMV5(JROF,JLEV,YT5%MVL)*ZGM(JROF)*ZGM(JROF)
    ENDIF
    IF(LVOR) THEN
      PGMV5(JROF,JLEV,YT5%MVOR)=PGMV5(JROF,JLEV,YT5%MVOR)*ZGM(JROF)*ZGM(JROF)
    ENDIF
  ENDDO
ENDDO

IF(LNHDYN) THEN
  DO JLEV=1,KFLEV
    DO JROF=KST,KEN
      PGMV5(JROF,JLEV,YT5%MSPDL)=PGMV5(JROF,JLEV,YT5%MSPDL)*ZGM(JROF)
      PGMV5(JROF,JLEV,YT5%MSPDM)=PGMV5(JROF,JLEV,YT5%MSPDM)*ZGM(JROF)
      PGMV5(JROF,JLEV,YT5%MSVDL)=PGMV5(JROF,JLEV,YT5%MSVDL)*ZGM(JROF)
      PGMV5(JROF,JLEV,YT5%MSVDM)=PGMV5(JROF,JLEV,YT5%MSVDM)*ZGM(JROF)
      IF (YRDYNA%LNHXDER) THEN
        PGMV5(JROF,JLEV,YT5%MNHXL)=PGMV5(JROF,JLEV,YT5%MNHXL)*ZGM(JROF)
        PGMV5(JROF,JLEV,YT5%MNHXM)=PGMV5(JROF,JLEV,YT5%MNHXM)*ZGM(JROF)
      ENDIF
    ENDDO
  ENDDO
ENDIF

! * GMV5S:

DO JROF=KST,KEN
  PGMV5S(JROF,YT5%MSPL)=PGMV5S(JROF,YT5%MSPL)*ZGM(JROF)
  PGMV5S(JROF,YT5%MSPM)=PGMV5S(JROF,YT5%MSPM)*ZGM(JROF)
ENDDO

! * GFL5:

DO JGFL=1,NUMFLDS
  IF(YCOMP(JGFL)%LCDERS) THEN
    DO JLEV=1,KFLEV
      DO JROF=KST,KEN
        PGFL5(JROF,JLEV,YCOMP(JGFL)%MP5L)=PGFL5(JROF,JLEV,YCOMP(JGFL)%MP5L)*ZGM(JROF)  
        PGFL5(JROF,JLEV,YCOMP(JGFL)%MP5M)=PGFL5(JROF,JLEV,YCOMP(JGFL)%MP5M)*ZGM(JROF)  
      ENDDO
    ENDDO
  ENDIF
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPMPFC5',1,ZHOOK_HANDLE)
END SUBROUTINE GPMPFC5
