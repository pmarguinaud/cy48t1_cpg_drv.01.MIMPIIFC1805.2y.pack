SUBROUTINE GNHGRPRE(&
 ! --- INPUT -------------------------------------------
 & YDGEOMETRY, KPROMA,KFLEV,KSTART,KEND,PRTGR,PREL,PREM,PNHPREF,PSPDL,PSPDM,&
 ! --- OUTPUT ------------------------------------------
 & PNHPREFL,PNHPREFM,PLNNHPREFL,PLNNHPREFM,PQCHAL,PQCHAM)

! GNHGRPRE - Computation of the horizontal gradient of the
!            total pressure "pre" at full levels.

! Purpose
! -------
!   Computes on model full levels.
!    * "grad pre" at full levels.
!    * the ratio [grad pre]/pre at full levels.
!    * the difference [grad pre]/pre - [grad prehyd]/prehyd
!      (which is equal to "grad(log(pre/prehyd))") at full levels.

!   Abbreviation "pdv" means "pressure departure variable"
!   (this is Qcha = log(pre/prehyd) for NPDVAR=2)

! Interface
! ---------
!   INPUT:
!    KPROMA     - horizontal dimension.
!    KFLEV      - number of levels.
!    KSTART     - start of work.
!    KEND       - end of work.
!    PRTGR      - the ratio [grad prehyd/prehyd]/[grad prehyds/prehyds]
!                 at full levels.
!    PREL       - zonal component of "grad (prehyds)".
!    PREM       - meridian component of "grad (prehyds)".
!    PNHPREF    - total pressure "pre" at full levels.
!    PSPDL      - zonal component of "grad pdv" at full levels.
!    PSPDM      - meridian component of "grad pdv" at full levels.

!   OUTPUT:
!    PNHPREFL   - zonal component of "grad pre" at full levels.
!    PNHPREFM   - merid component of "grad pre" at full levels.
!    PLNNHPREFL - zonal component of "[grad pre]/pre" at full levels.
!    PLNNHPREFM - merid component of "[grad pre]/pre" at full levels.
!    PQCHAL     - zonal comp of "grad(log(pre/prehyd))" at full levels.
!    PQCHAM     - merid comp of "grad(log(pre/prehyd))" at full levels.

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
!   K. Yessad (Dec 2016): Prune obsolete options.
!------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMDYNA      , ONLY : NPDVAR

! -----------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY),    INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRTGR(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREM(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHPREF(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDL(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPDM(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNHPREFL(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNHPREFM(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLNNHPREFL(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLNNHPREFM(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQCHAL(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQCHAM(KPROMA,KFLEV) 

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHGRPRE',0,ZHOOK_HANDLE)

! -----------------------------------------------------------------------------

IF ( NPDVAR == 2 ) THEN

  DO JLEV=1,KFLEV
    DO JROF=KSTART,KEND
      PLNNHPREFL(JROF,JLEV)=(PSPDL(JROF,JLEV)+PRTGR(JROF,JLEV)*PREL(JROF))
      PLNNHPREFM(JROF,JLEV)=(PSPDM(JROF,JLEV)+PRTGR(JROF,JLEV)*PREM(JROF))
      PNHPREFL(JROF,JLEV)=PNHPREF(JROF,JLEV)*PLNNHPREFL(JROF,JLEV)
      PNHPREFM(JROF,JLEV)=PNHPREF(JROF,JLEV)*PLNNHPREFM(JROF,JLEV)
      PQCHAL(JROF,JLEV)=PSPDL(JROF,JLEV)
      PQCHAM(JROF,JLEV)=PSPDM(JROF,JLEV)
    ENDDO
  ENDDO

ENDIF

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHGRPRE',1,ZHOOK_HANDLE)

END SUBROUTINE GNHGRPRE

