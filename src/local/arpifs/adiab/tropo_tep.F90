SUBROUTINE TROPO_TEP(PTOPPRES, PRPSTRA, PRPTROP, LDVERTFE, YDDYNA, YDCST,YDGEOMETRY,KPROMA,KSTART,KPROF,KFLEV,&
  & PSPT0,PTT0,KTROP)

!**** *TROPO_TEP* input for SLHD diffusion.
!                Computation of the model level above which
!                temperature is inversed. In a way searching
!                for tropopause.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *TROPO_TEP(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KPROMA      - horizontal dimension.
!          KSTART      - first element of work.
!          KPROF       - depth of work.
!          KFLEV       - number of layers.
!          PSPT0       - "ln(prehyds)" at t.
!          PTT0        - temperature at time t.

!        OUTPUT:
!          KTROP       - the lowermost stratospheric level.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           none
!           Called by LATTE_KAPPA.

!     Reference.
!     ----------
!             IFS documentation about semi-Lagrangian scheme.

!  Author.
!  -------
!    Original F. VANA : October 2013.

!  Modifications.
!  --------------

! End Modifications
!     ------------------------------------------------------------------
  
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM     ,JPRB
USE YOMHOOK      , ONLY : LHOOK    ,DR_HOOK

USE YOMVERT,ONLY: TVAB

USE YOMDYNA, ONLY : TDYNA
USE YOMCST,ONLY: TCST

!     ------------------------------------------------------------------

IMPLICIT NONE

REAL (KIND=JPRB), INTENT (IN) :: PTOPPRES
REAL (KIND=JPRB), INTENT (IN) :: PRPTROP
REAL (KIND=JPRB), INTENT (IN) :: PRPSTRA
LOGICAL, INTENT (IN) :: LDVERTFE
TYPE (TDYNA), INTENT (IN) :: YDDYNA
TYPE(TCST),        INTENT(IN)    :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT0(KPROMA,KFLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KTROP(KPROMA) 


!     -------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPRB) ::   ZPDIFFL, ZPDIFFU

REAL(KIND=JPRB) :: ZPREF(KPROMA,KFLEV),ZPREH (KPROMA,0:KFLEV)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     -------------------------------------------------------
#include "gphpre_expl.intfb.h"
!     -------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TROPO_TEP',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB)
!    --------------------------------------------------------

!*       0.   Basic settings

! Maximum diff in Pa between the ICAO troposfere and the diagnosed one 
ZPDIFFL=18000._JPRB  ! lower limit
ZPDIFFU=10000._JPRB  ! upper limit


!*       1.   Tropopause level detection

! Alternatives to this  computation could be found
! in vdfouter (by P. Bechtold or E. Holm). 

!*         Compute pressure
ZPREH (KSTART:KPROF,KFLEV) = EXP(PSPT0(KSTART:KPROF))
CALL GPHPRE_EXPL (YDDYNA%LAPRXPK, LDVERTFE, YDDYNA%NDLNPR, YDDYNA%RHYDR0, PTOPPRES, YDCST,KPROMA,KFLEV,KSTART,KPROF,YDVAB,ZPREH,PRESF=ZPREF)

!*         Find the lowest stratospheric level
KTROP(KSTART:KPROF)=1
DO JLEV=2,KFLEV
  DO JROF=KSTART,KPROF
    IF (ZPREF(JROF,KTROP(JROF)) < PRPSTRA-ZPDIFFU) KTROP(JROF)=JLEV
    IF ( (PTT0(JROF,JLEV) <= PTT0(JROF,JLEV-1)).AND. &
     &  (ZPREF(JROF,JLEV) < PRPTROP+ZPDIFFL))         &
     &  KTROP(JROF)=MAX(KTROP(JROF),JLEV)
  ENDDO
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TROPO_TEP',1,ZHOOK_HANDLE)
END SUBROUTINE TROPO_TEP
