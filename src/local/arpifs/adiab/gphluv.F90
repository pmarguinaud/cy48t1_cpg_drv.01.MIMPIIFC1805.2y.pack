SUBROUTINE GPHLUV(YDDIMV,KPROMA,KSTART,KPROF,PU,PV,PUVH)

!**** *GPHLUV* - wind components calculation in half-levels

USE YOMDIMV   , ONLY : TDIMV
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK
USE INTDYN_MOD, ONLY : YYTHW

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIMV)        ,INTENT(IN)     :: YDDIMV
INTEGER(KIND=JPIM) ,INTENT(IN)     :: KPROMA 
INTEGER(KIND=JPIM) ,INTENT(IN)     :: KSTART 
INTEGER(KIND=JPIM) ,INTENT(IN)     :: KPROF 
REAL(KIND=JPRB)    ,INTENT(IN)     :: PU(KPROMA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)     :: PV(KPROMA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(INOUT)  :: PUVH(KPROMA,0:YDDIMV%NFLEVG,YYTHW%NDIM) 

#include "gphluv_expl.intfb.h"

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPHLUV',0,ZHOOK_HANDLE)

CALL GPHLUV_EXPL (YDDIMV,KPROMA,KSTART,KPROF,PU,PV,PUVH(:,:,YYTHW%M_WWI),&
                & PUVH(:,:,YYTHW%M_UH),PUVH(:,:,YYTHW%M_VH))

IF (LHOOK) CALL DR_HOOK('GPHLUV',1,ZHOOK_HANDLE)

END SUBROUTINE GPHLUV
