SUBROUTINE GNHQE_PREH(   LDVERTFE, KPDVAR, YDGEOMETRY, KST, KEND, PREH, PQCHAF, PNHPREH, PREF, &
& PEQCHAF, PKAPH, PEQCHAH, PDEL_QCHAF)

! GNHQE_PREH - Computation of quantities linked to the total pressure "pre" at half levels,
!              from the "pressure departure" prognostic variable Qcha=log(pre/prehyd).
!              NHQE model.
!              Calculations are different from those of the NHEE model,
!              this is why the code is in a different routine.

! Purpose
! -------
!   Computation of quantities linked to the total pressure "pre" at half levels.
!   In the NHQE model we need at least:
!   * half level values of exp((R/Cp) Qcha) (=(pre/prehyd)**(R/Cp))
!   * full level values of [Delta exp((R/Cp) Qcha)].
!   * full level values of [Delta Qcha].
!   VFD and VFE discretisations must be taken into consideration.
!   The assumptions currently done are:
!    1/ VFD discretisation: 
!       the pressure departure prognostic variable at half
!        levels is the average of those at full levels.
!    2/ VFE discretisation: 
!       To be studied later; a "draft" code has been proposed.
!    3/ Qcha=log(pre/prehyd), and the pressure departure "pre - prehyd" are zero at the top of the model.

! Interface
! ---------
!   INPUT:
!    YDGEOMETRY     : structure containing geometry
!    KST            : start of work.
!    KEND           : end of work.
!    PREH           : hydrostatic pressure "prehyd" at half levels.
!    PQCHAF         : full level Qcha.

!   OUTPUT:
!    PNHPREH        : half level total pressure.

!   OPTIONAL INPUT:
!    PREF           : hydrostatic pressure "prehyd" at full levels.
!    PEQCHAF        : full level exp((R/Cp) Qcha).
!    PKAPH          : half level (R/Cp).

!   OPTIONAL OUTPUT:
!    PEQCHAH        : half level exp((R/Cp) Qcha) for horizontal pressure gradient term.
!    PDEL_QCHAF     : full level [Delta Qcha].

! Externals
! ---------

! Method
! ------

! Reference
! ---------

! Author
! ------
!   K. YESSAD (MF/CNRM/GMAP).
!   Original: March 2017

! Modifications
! -------------
!      K. Yessad (March 2018): rewrite Laplacian term in NHQE model.
! -----------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB


USE YOMHOOK      , ONLY : LHOOK, DR_HOOK

! -----------------------------------------------------------------------------

IMPLICIT NONE

LOGICAL            ,INTENT(IN)             :: LDVERTFE
INTEGER(KIND=JPIM) ,INTENT(IN)             :: KPDVAR
TYPE(GEOMETRY)     ,INTENT(IN)             :: YDGEOMETRY
INTEGER(KIND=JPIM) ,INTENT(IN)             :: KST 
INTEGER(KIND=JPIM) ,INTENT(IN)             :: KEND 
REAL(KIND=JPRB)    ,INTENT(IN)             :: PREH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)             :: PQCHAF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(OUT)            :: PNHPREH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)  ,OPTIONAL  :: PREF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)  ,OPTIONAL  :: PEQCHAF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)  ,OPTIONAL  :: PKAPH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(OUT) ,OPTIONAL  :: PEQCHAH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(OUT) ,OPTIONAL  :: PDEL_QCHAF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZQCHAH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

#include "abor1.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHQE_PREH', 0, ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)

! -----------------------------------------------------------------------------

IF ( KPDVAR /= 2 ) THEN
  ! This routine is coded only for NPDVAR=2 (prognostic variable = Qcha = log(pre/prehyd))
  CALL ABOR1(' GNHQE_PREH 0.1: NPDVAR should be equal to 2.')
ENDIF

! -----------------------------------------------------------------------------

!*     1. HALF LEVEL Qcha and total pressure.
!      --------------------------------------

! * Top (the pressure departure "pre - prehyd" is assumed to be zero):
ZQCHAH(KST:KEND,0)=0.0_JPRB

! * Half levels 1 to NFLEVG-1:
DO JLEV=1,NFLEVG-1
  DO JROF=KST,KEND
    ZQCHAH(JROF,JLEV)=0.5_JPRB*(PQCHAF(JROF,JLEV)+PQCHAF(JROF,JLEV+1))
  ENDDO
ENDDO

! * Surface: Qcha is assumed to be constant under the full level l=nflevg.
ZQCHAH(KST:KEND,NFLEVG)=PQCHAF(KST:KEND,NFLEVG)

! * Half level total pressure:
DO JLEV=0,NFLEVG
  DO JROF=KST,KEND
    PNHPREH(JROF,JLEV)=PREH(JROF,JLEV)*EXP(ZQCHAH(JROF,JLEV))
  ENDDO
ENDDO

! -----------------------------------------------------------------------------

!*     2. FULL LEVEL [Delta Qcha]
!      --------------------------

IF (PRESENT(PDEL_QCHAF)) THEN

  IF (LDVERTFE) THEN

    ! * VFE treatment of [Delta Qcha]:
    CALL ABOR1(' GNHQE_PREH: option not yet coded!')

  ELSE

    ! * VFD treatment of [Delta Qcha]:
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        PDEL_QCHAF(JROF,JLEV)=ZQCHAH(JROF,JLEV)-ZQCHAH(JROF,JLEV-1)
      ENDDO
    ENDDO

  ENDIF

ENDIF

! -----------------------------------------------------------------------------

!*     3. HALF LEVEL exp((R/Cp) Qcha).
!      -------------------------------

! For horizontal pressure gradient term:

IF (PRESENT(PEQCHAH).AND.PRESENT(PKAPH).AND.PRESENT(PEQCHAF)) THEN

  ! * Top (the pressure departure "pre - prehyd" is assumed to be zero):
  DO JROF=KST,KEND
    PEQCHAH(JROF,0)=1.0_JPRB
  ENDDO

  ! * Half levels 1 to NFLEVG-1:
  DO JLEV=1,NFLEVG-1
    DO JROF=KST,KEND
      PEQCHAH(JROF,JLEV)=EXP(PKAPH(JROF,JLEV)*ZQCHAH(JROF,JLEV))
    ENDDO
  ENDDO

  ! * Surface: Qcha is assumed to be constant under the full level l=nflevg.
  DO JROF=KST,KEND
    PEQCHAH(JROF,NFLEVG)=PEQCHAF(JROF,NFLEVG)
  ENDDO

ENDIF

! -----------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHQE_PREH', 1, ZHOOK_HANDLE)
END SUBROUTINE GNHQE_PREH

