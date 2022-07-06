SUBROUTINE GNHPREH(&
 ! --- INPUT -------------------------------------------
 & LDVERTFE, LDVFE_ECMWF, KPDVAR, YDGEOMETRY, KPROMA,KFLEV,KST,KEND,PSPD,PREH,PDELP,PLNPR,PNHPPI,PNHPRE,&
 ! --- OUTPUT ------------------------------------------
 & PNHPREH,PDELNHPRE&
 & )

! GNHPREH - Computation of the total pressure "pre" at half levels,
!           from the "pressure departure" prognostic variable.

! Purpose
! -------
!   Computes:
!    * the total pressure "pre" at half levels.
!    * the total pressure depth "Delta pre" at full levels.
!   The assumptions currently done are:
!    1/ LVERTFE=F: 
!       the pressure departure prognostic variable at half
!        levels is the average of those at full levels.
!    2/ LVERTFE=T: 
!       the total pressure "pre" at half levels is computed from the full
!       level total pressure depth "Delta pre" (first computed).
!    3/ The pressure departure "pre - prehyd" is zero at the top of the model.

! Interface
! ---------
!   INPUT:
!    KPROMA    - horizontal dimension.
!    KFLEV     - number of levels.
!    KST       - start of work.
!    KEND      - end of work.
!    PSPD      - pressure departure prognostic variable at full levels.
!    PREH      - hydrostatic pressure "prehyd" at half levels.
!    PDELP     - "delta prehyd" at full levels.
!    PLNPR     - delta = delta(log(prehyd)) at full levels.
!    PNHPPI    - ratio pre/prehyd at full levels.
!    PNHPRE    - NH pressure at full levels.

!   OUTPUT:
!    PNHPREH   - total pressure "pre" at half levels.
!    PDELNHPRE - total pressure depth "Delta pre" at full levels.

! Externals
! ---------

! Method
! ------

! Reference
! ---------

! Author
! ------
!   K. YESSAD (MF/CNRM/GMAP), Dec 2004.

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (Dec 2016): Prune obsolete options.
!   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   J. Vivoda and P. Smolikova (Sep 2020): VFE pruning.
! -----------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB


USE YOMHOOK      , ONLY : LHOOK, DR_HOOK

! -----------------------------------------------------------------------------

IMPLICIT NONE

LOGICAL, INTENT (IN) :: LDVERTFE
LOGICAL, INTENT (IN) :: LDVFE_ECMWF
INTEGER (KIND=JPIM), INTENT (IN) :: KPDVAR
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPD(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELP(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLNPR(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHPPI(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHPRE(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNHPREH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDELNHPRE(KPROMA,KFLEV) 

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZIN(KPROMA,0:KFLEV+1) 
REAL(KIND=JPRB) :: ZDX(KPROMA,KFLEV) 
REAL(KIND=JPRB) :: ZSPDH
REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

#include "verdisint.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHPREH',0,ZHOOK_HANDLE)

! -----------------------------------------------------------------------------

!*     1. FINITE DIFFERENCES VERTICAL DISCRETISATION
!      ---------------------------------------------

! Half level pressure "pre" is computed first; full level pressure depths
!  are then computed.

! * Computation of total pressure "pre" at half levels (PNHPREH):
!   Required for the time being as input for the AROME physics, but
!    in VFE discretisation, the use of this quantity is not recommanded
!    (in this case we use the FD discretisation, which does not take
!    account of the VFE calculation of PDELNHPRE, and we should notice
!    that the sum of PDELNHPRE*(delta eta) on a whole column does not
!    give exactly PNHPREH(surf)-PNHPREH(top)).

IF ( KPDVAR == 2 ) THEN

  ! * Top (the pressure departure "pre - prehyd" is assumed to be zero):
  PNHPREH(KST:KEND,0)=PREH(KST:KEND,0)

  ! * Half levels 1 to KFLEV-1:
  DO JLEV=1,KFLEV-1
    DO JROF=KST,KEND
      ZSPDH=0.5_JPRB*(PSPD(JROF,JLEV)+PSPD(JROF,JLEV+1))
      PNHPREH(JROF,JLEV)=EXP(ZSPDH)*PREH(JROF,JLEV)
    ENDDO
  ENDDO

  ! * Surface:
  !   The pressure departure variable is assumed to be constant
  !   under the full level l=nflevg.
  DO JROF=KST,KEND
    ZSPDH=PSPD(JROF,KFLEV)
    PNHPREH(JROF,KFLEV)=EXP(ZSPDH)*PREH(JROF,KFLEV)
  ENDDO

ENDIF

! * Computation of the total pressure depth "Delta pre" at full levels.
IF (.NOT.LDVERTFE) THEN
  DO JLEV=1,KFLEV
    DO JROF=KST,KEND
      PDELNHPRE(JROF,JLEV)=PNHPREH(JROF,JLEV)-PNHPREH(JROF,JLEV-1)
    ENDDO
  ENDDO
ENDIF

! -----------------------------------------------------------------------------

!*     2. VFE VERTICAL DISCRETISATION
!      ------------------------------

! Full level pressure depths "Delta pre" are computed first;
!  half level pressure "pre" is then computed, starting from the
!  upper boundary condition "pre - prehyd=0".

IF (LDVERTFE) THEN

  IF (LDVFE_ECMWF) THEN

    ! * Computation of the total pressure depth "Delta pre" at full levels.
    !   Calculation is done using the identity:
    !    d pre / d prehyd = pre/prehyd + d ((pre-prehyd)/prehyd)/ d (log(prehyd))
    !   then:
    !    Delta pre = Delta prehyd * [d pre / d prehyd]

    DO JLEV=1,KFLEV
      DO JROF=KST,KEND
        ZIN(JROF,JLEV)=PNHPPI(JROF,JLEV)-1.0_JPRB
      ENDDO
    ENDDO
    ! * Top: "[pre/prehyd]_top-1=0":
    ZIN(KST:KEND,0)=0.0_JPRB
    ! * Bottom: "[pre/prehyd]_surf=[pre/prehyd]_(l=L)":
    ZIN(KST:KEND,KFLEV+1)=ZIN(KST:KEND,KFLEV)
    CALL VERDISINT(YDGEOMETRY%YRVFE,'FDER','00',KPROMA,KST,KEND,KFLEV,ZIN,ZDX)
    DO JLEV=1,KFLEV
      DO JROF=KST,KEND
        PDELNHPRE(JROF,JLEV)=PDELP(JROF,JLEV)*(PNHPPI(JROF,JLEV)&
         & +ZDX(JROF,JLEV)/(PLNPR(JROF,JLEV)*YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)))
      ENDDO
    ENDDO

  ELSE

    IF( KPDVAR == 2 )THEN
      ! * Computation of the total pressure depth "Delta pre" at full levels.
      !   Calculation is done using the identity:
      !    Delta pre = Delta prehyd * (pre/prehyd) + pre * Delta (log (pre/prehyd))
      DO JLEV=1,KFLEV
        DO JROF=KST,KEND
          ZIN(JROF,JLEV)=PSPD(JROF,JLEV)
        ENDDO
      ENDDO
      ZIN(KST:KEND,0)=0.0_JPRB
      ZIN(KST:KEND,KFLEV+1)=0.0_JPRB
      CALL VERDISINT(YDGEOMETRY%YRVFE,'FDER','01',KPROMA,KST,KEND,KFLEV,ZIN,ZDX)
      DO JLEV=1,KFLEV
        DO JROF=KST,KEND
          PDELNHPRE(JROF,JLEV)=PDELP(JROF,JLEV)*PNHPPI(JROF,JLEV) &
           & + PNHPRE(JROF,JLEV)*ZDX(JROF,JLEV)/YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)
        ENDDO
      ENDDO
    ENDIF

  ENDIF

ENDIF

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHPREH',1,ZHOOK_HANDLE)

END SUBROUTINE GNHPREH

