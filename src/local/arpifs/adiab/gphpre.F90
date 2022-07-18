SUBROUTINE GPHPRE(KPROMA,KFLEV,K1,K2,VAB,PRESH,PXYB,PRESF,LHSET,LDELP,LALPHA,LRTGR,LRPP)

!**** *GPHPRE* - Computes half and full level pressure
!                Modern version of former GPPRE.
!                Modern version of former GPPREH+GPXYB+GPPREF

!     Purpose.
!     --------
!           Computes pressures at half and full model levels.

!**   Interface.
!     ----------
!        *CALL* *GPHPRE(...)

!        Explicit arguments :
!        --------------------

!          KPROMA    : horizontal dimensioning                                (in)
!          KFLEV     : vertical dimensioning                                  (in)
!          KSTART    : start of work                                          (in)
!          KPROF     : depth of work                                          (in)
!          YDVAB     : contains information about hybrid vertical coordinate  (in)
!          PRESH     : half level pressure                                    (inout)
!          PXYB      : contains pressure depth, "delta", "alpha"              (opt out)
!          PRESF     : full level pressure                                    (opt out)
!          LDELP,LALPHA,... : activation keys for partial computations        (opt in)

!        Implicit arguments :  NONE.
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
!      K. YESSAD (Sep 2011) after GPPRE, GPPREH, GPXYB and GPPREF.

!     Modifications.
!     --------------
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (Mar 2017): Introduce NDLNPR=2 for NHQE model.
!   H Petithomme (Dec 2020): add options, use of pointers, group VFE tests
!     ------------------------------------------------------------------

USE PARKIND1,ONLY: JPIM,JPRB
USE YOMHOOK,ONLY: LHOOK,DR_HOOK
USE YOMVERT,ONLY: TVAB,TOPPRES,LVERTFE
USE YOMDYNA, ONLY : YRDYNA
USE INTDYN_MOD,ONLY: YYTXYB
USE YOMCST,ONLY: YRCST

IMPLICIT NONE

TYPE(TVAB),                     INTENT(IN)    :: VAB
INTEGER(KIND=JPIM),             INTENT(IN)    :: KPROMA,KFLEV,K1,K2
LOGICAL,OPTIONAL,               INTENT(IN)    :: LHSET,LDELP,LALPHA,LRTGR,LRPP
REAL(KIND=JPRB),                INTENT(INOUT) :: PRESH(KPROMA,0:KFLEV)
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(OUT)   :: PXYB(KPROMA,KFLEV,YYTXYB%NDIM)
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(OUT)   :: PRESF(KPROMA,KFLEV)

#include "gphpre_expl.intfb.h"

REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('GPHPRE',0,ZHOOK_HANDLE)

IF (PRESENT (PXYB)) THEN
  CALL GPHPRE_EXPL (YRDYNA%LAPRXPK, LVERTFE, YRDYNA%NDLNPR, YRDYNA%RHYDR0, TOPPRES, YRCST,KPROMA,KFLEV,K1,K2,VAB,PRESH,PRESF,&
                  & LHSET,LDELP,LALPHA,LRTGR,LRPP,&
                  & PXYB (:,:,YYTXYB%M_DELP ), PXYB (:,:,YYTXYB%M_LNPR ), PXYB (:,:,YYTXYB%M_RDELP), &
                  & PXYB (:,:,YYTXYB%M_ALPH ), PXYB (:,:,YYTXYB%M_RTGR ), PXYB (:,:,YYTXYB%M_RPRE ), &
                  & PXYB (:,:,YYTXYB%M_RPP  ))
ELSE
  CALL GPHPRE_EXPL (YRDYNA%LAPRXPK, LVERTFE, YRDYNA%NDLNPR, YRDYNA%RHYDR0, TOPPRES, YRCST,KPROMA,KFLEV,K1,K2,VAB,PRESH,&
                  & PRESF,LHSET,LDELP,LALPHA,LRTGR,LRPP)
ENDIF


IF (LHOOK) CALL DR_HOOK('GPHPRE',1,ZHOOK_HANDLE)

END SUBROUTINE

