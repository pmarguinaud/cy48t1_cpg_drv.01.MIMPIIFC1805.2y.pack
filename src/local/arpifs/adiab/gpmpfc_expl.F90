SUBROUTINE GPMPFC_EXPL(YDGMV,YDML_GCONF,YDDYN,KPROMA,KFLEV,KST,KEN,KFLAG,PGM,PGFL,&
                     & P0U, P0V, P0DIV, P0TL, P0TM, P9U, P9V, P0UL, P0VL, P0VOR, P0SPDL, P0SPDM, &
                     & P0SVDL, P0SVDM, P0NHXL, P0NHXM, P9DIV, P9TL, P9TM, P9SPDL, P9SPDM, P9SVDL,&
                     & P9SVDM, P0SPL, P0SPM, P9SPL, P9SPM)

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
REAL(KIND=JPRB),INTENT(INOUT) :: PGFL(KPROMA,KFLEV,YDML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB),OPTIONAL,TARGET,INTENT(INOUT),DIMENSION(KPROMA,KFLEV) :: P0U, P0V, P0DIV, P0TL, P0TM, P9U, P9V, P0UL, P0VL, P0VOR, P0SPDL, P0SPDM
REAL(KIND=JPRB),OPTIONAL,TARGET,INTENT(INOUT),DIMENSION(KPROMA,KFLEV) :: P0SVDL, P0SVDM, P0NHXL, P0NHXM, P9DIV, P9TL, P9TM, P9SPDL, P9SPDM, P9SVDL, P9SVDM
REAL(KIND=JPRB),OPTIONAL,TARGET,INTENT(INOUT),DIMENSION(KPROMA) :: P0SPL, P0SPM, P9SPL, P9SPM
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JLEV,JGFL,JROF
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZGM(:)
REAL(KIND=JPRB),TARGET :: ZGM0(KPROMA)
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPMPFC_EXPL',0,ZHOOK_HANDLE)
ASSOCIATE(NDIM=>YDML_GCONF%YGFL%NDIM, NUMFLDS=>YDML_GCONF%YGFL%NUMFLDS, YCOMP=>YDML_GCONF%YGFL%YCOMP, &
 & LUVDER=>YDML_GCONF%YRDIMF%LUVDER, LVOR=>YDML_GCONF%YRDIMF%LVOR, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0, &
 & YT9=>YDGMV%YT9)
!     ------------------------------------------------------------------

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

! * GMV (unconditional t0 and t9 merged):

DO JLEV=1,KFLEV
  DO JROF=KST,KEN
    ! mapping of gmv t0
    P0U(JROF,JLEV)=P0U(JROF,JLEV)*ZGM(JROF)
    P0V(JROF,JLEV)=P0V(JROF,JLEV)*ZGM(JROF)
    P0DIV(JROF,JLEV)=P0DIV(JROF,JLEV)*ZGM(JROF)*ZGM(JROF)
    P0TL(JROF,JLEV)=P0TL(JROF,JLEV)*ZGM(JROF)
    P0TM(JROF,JLEV)=P0TM(JROF,JLEV)*ZGM(JROF)

    ! mapping of gmv t9
    P9U(JROF,JLEV)=P9U(JROF,JLEV)*ZGM(JROF)
    P9V(JROF,JLEV)=P9V(JROF,JLEV)*ZGM(JROF)
  ENDDO

  IF(LUVDER) THEN
    DO JROF=KST,KEN
      P0UL(JROF,JLEV)=P0UL(JROF,JLEV)*ZGM(JROF)*ZGM(JROF)
      P0VL(JROF,JLEV)=P0VL(JROF,JLEV)*ZGM(JROF)*ZGM(JROF)
    ENDDO
  ENDIF

  IF(LVOR) THEN
    DO JROF=KST,KEN
      P0VOR(JROF,JLEV)=P0VOR(JROF,JLEV)*ZGM(JROF)*ZGM(JROF)
    ENDDO
  ENDIF

  IF(LNHDYN) THEN
    DO JROF=KST,KEN
      P0SPDL(JROF,JLEV)=P0SPDL(JROF,JLEV)*ZGM(JROF)
      P0SPDM(JROF,JLEV)=P0SPDM(JROF,JLEV)*ZGM(JROF)
      P0SVDL(JROF,JLEV)=P0SVDL(JROF,JLEV)*ZGM(JROF)
      P0SVDM(JROF,JLEV)=P0SVDM(JROF,JLEV)*ZGM(JROF)
    ENDDO

    IF (LNHXDER) THEN
      DO JROF=KST,KEN
        P0NHXL(JROF,JLEV)=P0NHXL(JROF,JLEV)*ZGM(JROF)
        P0NHXM(JROF,JLEV)=P0NHXM(JROF,JLEV)*ZGM(JROF)
      ENDDO
    ENDIF
  ENDIF
ENDDO

IF(.NOT.LTWOTL) THEN
  DO JLEV=1,KFLEV
    DO JROF=KST,KEN
      P9DIV(JROF,JLEV)=P9DIV(JROF,JLEV)*ZGM(JROF)*ZGM(JROF)
      P9TL(JROF,JLEV)=P9TL(JROF,JLEV)*ZGM(JROF)
      P9TM(JROF,JLEV)=P9TM(JROF,JLEV)*ZGM(JROF)
    ENDDO

    IF(LNHDYN) THEN
      DO JROF=KST,KEN
        P9SPDL(JROF,JLEV)=P9SPDL(JROF,JLEV)*ZGM(JROF)
        P9SPDM(JROF,JLEV)=P9SPDM(JROF,JLEV)*ZGM(JROF)
        P9SVDL(JROF,JLEV)=P9SVDL(JROF,JLEV)*ZGM(JROF)
        P9SVDM(JROF,JLEV)=P9SVDM(JROF,JLEV)*ZGM(JROF)
      ENDDO
    ENDIF
  ENDDO
ENDIF

! * GMVS:

DO JROF=KST,KEN
  P0SPL(JROF)=P0SPL(JROF)*ZGM(JROF)
  P0SPM(JROF)=P0SPM(JROF)*ZGM(JROF)
ENDDO

IF (.NOT.LTWOTL.OR.NCURRENT_ITER > 0.AND.LPC_FULL) THEN
  DO JROF=KST,KEN
    P9SPL(JROF)=P9SPL(JROF)*ZGM(JROF)
    P9SPM(JROF)=P9SPM(JROF)*ZGM(JROF)
  ENDDO
ENDIF

! * GFL:

DO JGFL=1,NUMFLDS
  IF(YCOMP(JGFL)%LCDERS) THEN
    DO JLEV=1,KFLEV
      DO JROF=KST,KEN
        PGFL(JROF,JLEV,YCOMP(JGFL)%MPL)=PGFL(JROF,JLEV,YCOMP(JGFL)%MPL)*ZGM(JROF)  
        PGFL(JROF,JLEV,YCOMP(JGFL)%MPM)=PGFL(JROF,JLEV,YCOMP(JGFL)%MPM)*ZGM(JROF)  
      ENDDO
    ENDDO
  ENDIF
ENDDO
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('GPMPFC_EXPL',1,ZHOOK_HANDLE)
END SUBROUTINE GPMPFC_EXPL

