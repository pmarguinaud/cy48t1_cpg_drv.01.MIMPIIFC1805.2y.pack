SUBROUTINE GPHLUV_EXPL (YDDIMV,KPROMA,KSTART,KPROF,PU,PV,PWWI,PUH,PVH)

!**** *GPHLUV_EXPL* - wind components calculation in half-levels

!     Purpose.
!     --------
!           Compute wind components in half-levels

!**   Interface.
!     ----------
!        *CALL* *GPHLUV_EXPL()

!        Explicit arguments :
!        --------------------
!        * INPUT:
!          KPROMA  - length of work
!          KSTART  - start of work
!          KPROF   - end of work
!          PU      - U-wind at full levels
!          PV      - V-wind at full levels

!        * IN/OUT:
!          PUVH    - horizontal wind and weights at half levels
!                    (IN) for weights, (OUT) for half level wind

!        Implicit arguments : none.
!        --------------------

!     Method.
!     -------
!        Interpolation from full-levels using logaritmic pressure profile
!        Then modify values on the top and bottom using boundary condition
!        (at the bottom free-slip condition)
!        Store also the weights of vertical interpolation

!     Externals.
!     ----------
!     Reference.
!     ----------
!        ARPEGE/ALADIN documentation

!     Author.
!     -------
!      Radmila Bubnova & Martin Janousek,  CNRM/GMAP/EXT
!      Original : February 1996

! Modifications
! -------------
!   Modified 02-07-02 by C. Fischer : intents for dummy arrays
!   Modified 08-2002 C. Smith and K. YESSAD:
!    - make optional the vertical averaging.
!    - some cleanings and optimisations.
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!   09-Jun-2004 J. Masek   NH cleaning (LVSLWBC)
!   21-Jan-2005 K. Yessad  Remove useless dummy arguments.
!   07-Mar-2007 K. Yessad  Remove LVERAVE_HLUV.
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
! End Modifications
!------------------------------------------------------------------

USE YOMDIMV   , ONLY : TDIMV
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIMV),       INTENT(IN)    :: YDDIMV
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KPROMA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KPROMA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWWI(KPROMA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUH(KPROMA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVH(KPROMA,0:YDDIMV%NFLEVG)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF, JLEV
REAL(KIND=JPRB) :: ZWW

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPHLUV_EXPL',0,ZHOOK_HANDLE)

ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG)

!     ------------------------------------------------------------------

!*    1. General case:

DO JLEV=1, NFLEVG - 1
  DO JROF=KSTART,KPROF
    ZWW=PWWI(JROF,JLEV)
    PUH(JROF,JLEV)=ZWW*PU(JROF,JLEV)+(1.0_JPRB-ZWW)*PU(JROF,JLEV+1)  
    PVH(JROF,JLEV)=ZWW*PV(JROF,JLEV)+(1.0_JPRB-ZWW)*PV(JROF,JLEV+1)  
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*    2. Top.

PUH(KSTART:KPROF,0) = PU(KSTART:KPROF,1)
PVH(KSTART:KPROF,0) = PV(KSTART:KPROF,1)

!     ------------------------------------------------------------------

!*    3. Surface.

PUH(KSTART:KPROF,NFLEVG) = PU(KSTART:KPROF,NFLEVG)
PVH(KSTART:KPROF,NFLEVG) = PV(KSTART:KPROF,NFLEVG)

!     ------------------------------------------------------------------

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('GPHLUV_EXPL',1,ZHOOK_HANDLE)

END SUBROUTINE GPHLUV_EXPL
