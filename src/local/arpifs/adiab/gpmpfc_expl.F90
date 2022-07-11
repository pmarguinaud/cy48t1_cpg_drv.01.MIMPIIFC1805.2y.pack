SUBROUTINE GPMPFC_EXPL (LDNHDYN, LDNHXDER, LDPC_FULL, LDTWOTL, YDVARS, YDML_GCONF, YDDYN, KPROMA, &
& KFLEV, KST, KEN, KFLAG, PGM)

!**** *GPMPFC_EXPL* - Apply map factor to convert 
!                reduced variables -> geographical variables if kflag=0
!                or
!                geographical variables -> reduced variables if kflag=1

!     Purpose.
!     --------
!           Multiply or divide by map factor.

!**   Interface.
!     ----------
!        *CALL* *GPMPFC_EXPL(...)

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
!       PGMV      - GMV variables
!       PGMVS     - GMVS variables
!       PGFL      - GFL variables

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
!      Mats Hamrud  *ECMWF*
!      Original : 1994-01-18

! Modifications
! -------------
!   Modified 2002-07-02 by C. Fischer  : rename NHS variables T0/T9
!   Modified 2002-11-13 by K. YESSAD   : some cleanings + improve vectorization
!   Modified 2003-08    by M. HAMRUD   : GFL
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!   08-Jun-2004 J. Masek   NH cleaning (LVSLWBC)
!   01-Jul-2004 K. Yessad  Make clearer the tests for PC scheme.
!   09-Feb-2006 M. Deque : Dasux compilance
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (Jun 2011): new dataflow, GMVS too.
!   K. Yessad (Nov 2012): simplify testings.
!   K. Yessad (July 2014): Move some variables.
!   H. Petithomme (Dec 2020): optimisation
! End Modifications
!     ------------------------------------------------------------------
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK


USE YOMDYN                 , ONLY : TDYN
USE FIELD_VARIABLES_MOD    , ONLY : FIELD_VARIABLES

!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL                       ,INTENT(IN)             :: LDNHDYN
LOGICAL                       ,INTENT(IN)             :: LDNHXDER
LOGICAL                       ,INTENT(IN)             :: LDPC_FULL
LOGICAL                       ,INTENT(IN)             :: LDTWOTL
TYPE(FIELD_VARIABLES)         ,INTENT(INOUT)          :: YDVARS
TYPE(TDYN)                    ,INTENT(IN)             :: YDDYN
TYPE(MODEL_GENERAL_CONF_TYPE) ,INTENT(IN)             :: YDML_GCONF
INTEGER(KIND=JPIM)            ,INTENT(IN)             :: KPROMA
INTEGER(KIND=JPIM)            ,INTENT(IN)             :: KFLEV
INTEGER(KIND=JPIM)            ,INTENT(IN)             :: KST 
INTEGER(KIND=JPIM)            ,INTENT(IN)             :: KEN 
INTEGER(KIND=JPIM)            ,INTENT(IN)             :: KFLAG 
REAL(KIND=JPRB)               ,INTENT(IN)    ,TARGET  :: PGM(KPROMA)


#include "abor1.intfb.h"
#include "gpmpfc_expl_part1.intfb.h"
#include "gpmpfc_expl_part2.intfb.h"

REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZGM(:)
REAL(KIND=JPRB),TARGET :: ZGM0(KPROMA)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('GPMPFC_EXPL', 0, ZHOOK_HANDLE)

!*       1. APPLY MAP FACTOR.
!           -----------------

IF(KFLAG == 0) THEN
  ZGM => PGM(:)
ELSEIF(KFLAG == 1) THEN
  ZGM0(KST:KEN) = 1.0_JPRB/PGM(KST:KEN)
  ZGM => ZGM0(:)
ELSE
  CALL ABOR1('GPMPFC_EXPL: ILLEGAL KFLAG')
ENDIF

CALL GPMPFC_EXPL_PART1 (LDNHDYN, LDNHXDER, LDPC_FULL, LDTWOTL, YDML_GCONF, YDDYN, KPROMA, KFLEV, KST,              &
& KEN, ZGM, P0U=YDVARS%U%T0, P0V=YDVARS%V%T0, P0DIV=YDVARS%DIV%T0, P0TL=YDVARS%T%DL, P0TM=YDVARS%T%DM,             &
& P9U=YDVARS%U%T9, P9V=YDVARS%V%T9, P0UL=YDVARS%U%DL, P0VL=YDVARS%V%DL, P0VOR=YDVARS%VOR%T0, P0SPDL=YDVARS%SPD%DL, &
& P0SPDM=YDVARS%SPD%DM, P0SVDL=YDVARS%SVD%DL, P0SVDM=YDVARS%SVD%DM, P0NHXL=YDVARS%NHX%DL, P0NHXM=YDVARS%NHX%DM,    &
& P9DIV=YDVARS%DIV%T9, P9TL=YDVARS%T%DL9, P9TM=YDVARS%T%DM9, P9SPDL=YDVARS%SPD%DL9, P9SPDM=YDVARS%SPD%DM9,         &
& P9SVDL=YDVARS%SVD%DL9, P9SVDM=YDVARS%SVD%DM9, P0SPL=YDVARS%SP%DL, P0SPM=YDVARS%SP%DM, P9SPL=YDVARS%SP%DL9,       &
& P9SPM=YDVARS%SP%DM9)

CALL GPMPFC_EXPL_PART2 (YDVARS, KPROMA, KFLEV, KST, KEN, PGM)

IF (LHOOK) CALL DR_HOOK('GPMPFC_EXPL', 1, ZHOOK_HANDLE)

END SUBROUTINE GPMPFC_EXPL

