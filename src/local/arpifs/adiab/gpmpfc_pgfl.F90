SUBROUTINE GPMPFC_PGFL(YDGMV,YDML_GCONF,YDDYN,KPROMA,KFLEV,KST,KEN,KFLAG,PGM,PGMV,PGMVS,PGFL)

!**** *GPMPFC_PGFL* - Apply map factor to convert 
!                reduced variables -> geographical variables if kflag=0
!                or
!                geographical variables -> reduced variables if kflag=1

!     Purpose.
!     --------
!           Multiply or divide by map factor.

!**   Interface.
!     ----------
!        *CALL* *GPMPFC_PGFL(...)

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
USE YOMGMV                 , ONLY : TGMV
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK
USE YOMCT0                 , ONLY : LTWOTL, LNHDYN
USE YOMDYNA                , ONLY : LNHXDER, LPC_FULL
USE YOMDYN                 , ONLY : TDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TGMV) , INTENT(INOUT) :: YDGMV
TYPE(TDYN)  ,INTENT(IN)    :: YDDYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLAG 
REAL(KIND=JPRB),TARGET,INTENT(IN) :: PGM(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(KPROMA,KFLEV,YDML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB),TARGET,INTENT(INOUT) :: PGMV(KPROMA,KFLEV,YDGMV%NDIMGMV)
REAL(KIND=JPRB),TARGET,INTENT(INOUT) :: PGMVS(KPROMA,YDGMV%NDIMGMVS)

#include "abor1.intfb.h"
#include "gpmpfc_expl_part1.intfb.h"
#include "gpmpfc_pgfl_part2.intfb.h"

REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z0U, Z0V, Z0DIV, Z0TL, Z0TM, Z9U, Z9V, Z0UL, Z0VL, Z0VOR, Z0SPDL, Z0SPDM
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z0SVDL, Z0SVDM, Z0NHXL, Z0NHXM, Z9DIV, Z9TL, Z9TM, Z9SPDL, Z9SPDM, Z9SVDL, Z9SVDM
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:) :: Z0SPL, Z0SPM, Z9SPL, Z9SPM

REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZGM(:)
REAL(KIND=JPRB),TARGET :: ZGM0(KPROMA)


REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('GPMPFC_PGFL',0,ZHOOK_HANDLE)

ASSOCIATE(NDIM=>YDML_GCONF%YGFL%NDIM, NUMFLDS=>YDML_GCONF%YGFL%NUMFLDS, YCOMP=>YDML_GCONF%YGFL%YCOMP, &
 & LUVDER=>YDML_GCONF%YRDIMF%LUVDER, LVOR=>YDML_GCONF%YRDIMF%LVOR, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0, &
 & YT9=>YDGMV%YT9)

Z0U    => NULL ()
Z0V    => NULL ()
Z0DIV  => NULL ()
Z0TL   => NULL ()
Z0TM   => NULL ()
Z9U    => NULL ()
Z9V    => NULL ()
Z0UL   => NULL ()
Z0VL   => NULL ()
Z0VOR  => NULL ()
Z0SPDL => NULL ()
Z0SPDM => NULL ()
Z0SVDL => NULL ()
Z0SVDM => NULL ()
Z0NHXL => NULL ()
Z0NHXM => NULL ()
Z9DIV  => NULL ()
Z9TL   => NULL ()
Z9TM   => NULL ()
Z9SPDL => NULL ()
Z9SPDM => NULL ()
Z9SVDL => NULL ()
Z9SVDM => NULL ()
Z0SPL  => NULL ()
Z0SPM  => NULL ()
Z9SPL  => NULL ()
Z9SPM  => NULL ()

IF (YT0%MU    > 0) Z0U    => PGMV(:,:,YT0%MU)
IF (YT0%MV    > 0) Z0V    => PGMV(:,:,YT0%MV)
IF (YT0%MDIV  > 0) Z0DIV  => PGMV(:,:,YT0%MDIV)
IF (YT0%MTL   > 0) Z0TL   => PGMV(:,:,YT0%MTL)
IF (YT0%MTM   > 0) Z0TM   => PGMV(:,:,YT0%MTM)
IF (YT9%MU    > 0) Z9U    => PGMV(:,:,YT9%MU)
IF (YT9%MV    > 0) Z9V    => PGMV(:,:,YT9%MV)
IF (YT0%MUL   > 0) Z0UL   => PGMV(:,:,YT0%MUL)
IF (YT0%MVL   > 0) Z0VL   => PGMV(:,:,YT0%MVL)
IF (YT0%MVOR  > 0) Z0VOR  => PGMV(:,:,YT0%MVOR)
IF (YT0%MSPDL > 0) Z0SPDL => PGMV(:,:,YT0%MSPDL)
IF (YT0%MSPDM > 0) Z0SPDM => PGMV(:,:,YT0%MSPDM)
IF (YT0%MSVDL > 0) Z0SVDL => PGMV(:,:,YT0%MSVDL)
IF (YT0%MSVDM > 0) Z0SVDM => PGMV(:,:,YT0%MSVDM)
IF (YT0%MNHXL > 0) Z0NHXL => PGMV(:,:,YT0%MNHXL)
IF (YT0%MNHXM > 0) Z0NHXM => PGMV(:,:,YT0%MNHXM)
IF (YT9%MDIV  > 0) Z9DIV  => PGMV(:,:,YT9%MDIV)
IF (YT9%MTL   > 0) Z9TL   => PGMV(:,:,YT9%MTL)
IF (YT9%MTM   > 0) Z9TM   => PGMV(:,:,YT9%MTM)
IF (YT9%MSPDL > 0) Z9SPDL => PGMV(:,:,YT9%MSPDL)
IF (YT9%MSPDM > 0) Z9SPDM => PGMV(:,:,YT9%MSPDM)
IF (YT9%MSVDL > 0) Z9SVDL => PGMV(:,:,YT9%MSVDL)
IF (YT9%MSVDM > 0) Z9SVDM => PGMV(:,:,YT9%MSVDM)
IF (YT0%MSPL  > 0) Z0SPL  => PGMVS(:,YT0%MSPL)
IF (YT0%MSPM  > 0) Z0SPM  => PGMVS(:,YT0%MSPM)
IF (YT9%MSPL  > 0) Z9SPL  => PGMVS(:,YT9%MSPL)
IF (YT9%MSPM  > 0) Z9SPM  => PGMVS(:,YT9%MSPM)

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

CALL GPMPFC_EXPL_PART1 (LNHDYN, LNHXDER, LPC_FULL, LTWOTL, YDML_GCONF, YDDYN, KPROMA, KFLEV, KST, KEN, ZGM, &
               & Z0U, Z0V, Z0DIV, Z0TL, Z0TM, Z9U, Z9V, Z0UL, Z0VL, Z0VOR, Z0SPDL, Z0SPDM, &
               & Z0SVDL, Z0SVDM, Z0NHXL, Z0NHXM, Z9DIV, Z9TL, Z9TM, Z9SPDL, Z9SPDM, Z9SVDL,&
               & Z9SVDM, Z0SPL, Z0SPM, Z9SPL, Z9SPM)

CALL GPMPFC_PGFL_PART2 (YDML_GCONF, KPROMA, KFLEV, KST, KEN, PGM, PGFL)

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('GPMPFC_PGFL',1,ZHOOK_HANDLE)

END SUBROUTINE GPMPFC_PGFL

