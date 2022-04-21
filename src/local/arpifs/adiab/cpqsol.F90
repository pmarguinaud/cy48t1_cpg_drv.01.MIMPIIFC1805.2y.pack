SUBROUTINE CPQSOL(YDCST, YDDIMV,YDPHY,KPROMA,KSTART,KPROF,PRES,PTS0,PQS,PQSATS,PQSOL)

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCST   , ONLY : TCST  
USE YOMPHY   , ONLY : TPHY

!**** *CPQSOL* - computation of the surface specific humidity.

!     Purpose.
!     --------
!           Computes the surface specific humidity necessary to compute
!           the vertical advection of humidity in the case "delta m = 1".
!           Assumes that the relative humidity is constant between the
!           layer number NFLEVG and the ground.

!**   Interface.
!     ----------
!        *CALL* *CPQSOL(KPROMA,KSTART,KPROF
!                      ,PRES,PTS0,PQS,PQSATS,PQSOL)

!        Explicit arguments :
!        --------------------
!          KPROMA  - horizontal dimension.                   (input)
!          KSTART  - first element of work.                  (input)
!          KPROF   - depth of work.                          (input)
!          PRES    - pressure on an interlayer.              (input)
!          PTS0    - surface temperature.                    (input)
!          PQS     - initial surface specific humidity,
!                    output of t-dt physics.                 (input)
!          PQSATS  - initial surface saturated specific
!                    humidity, output of t-dt physics.       (input)
!          PQSOL   - final surface specific humidity.        (output)

!        Implicit arguments : none.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!     Reference.
!     ----------
!        ARPEGE documentation

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : K. YESSAD (after CPDYN), 11 DEC 1992.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TPHY)        ,INTENT(IN)    :: YDPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRES(KPROMA,0:YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSATS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQSOL(KPROMA) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF

REAL(KIND=JPRB) :: ZDELT0, ZES0, ZHU, ZQSATS0
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "fcttrm.ycst.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CPQSOL',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & LNEIGE=>YDPHY%LNEIGE)
!     ------------------------------------------------------------------

DO JROF=KSTART,KPROF

  ZHU=PQS(JROF)*(1.0_JPRB+YDCST%RETV*PQSATS(JROF))/(PQSATS(JROF)*&
   & (1.0_JPRB+YDCST%RETV*PQS(JROF)))  

  IF (LNEIGE) THEN
    ZDELT0 = MAX(0.0_JPRB,SIGN(1.0_JPRB,YDCST%RTT-PTS0(JROF)))
  ELSE
    ZDELT0 = 0.0_JPRB
  ENDIF

  ZES0 = FOEW(PTS0(JROF),ZDELT0)
  ZQSATS0=FOQS(ZES0/PRES(JROF,NFLEVG))
  PQSOL(JROF)=ZHU*ZQSATS0/(1.0_JPRB+YDCST%RETV*ZQSATS0*(1.0_JPRB-ZHU))

ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPQSOL',1,ZHOOK_HANDLE)
END SUBROUTINE CPQSOL
