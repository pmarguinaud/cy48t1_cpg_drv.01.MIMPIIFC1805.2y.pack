!OCL  NOEVAL
SUBROUTINE GNHPRE(   KPDVAR, YDGEOMETRY, KPROMA, KFLEV, KSTART, KEND, PSPD, PREF, PKAP, PNHPREF, &
& PNHPPI, PRNHPPI, PDEP, PEQCHAF, PEIQCHAF)

! GNHPRE - Computation of the total pressure "pre" at full levels,
!          from the 'PD' prognostic variable.

! Purpose
! -------
!   Computes at model full levels.
!    * the total pressure "pre".
!    * the pressure departure "pre - prehyd".
!    * the ratios pre/prehyd and prehyd/pre.
!    * exp((R/cp) log(pre/prehyd)) and its inverse.

! Interface
! ---------
!   INPUT:
!    YDGEOMETRY   : structure containing all geometry.
!    KPROMA       : horizontal dimension.
!    KFLEV        : number of levels.
!    KSTART       : start of work.
!    KEND         : end of work.
!    PSPD         : NH pressure departure prognostic variable.
!    PREF         : hydrostatic pressure "prehyd" at full levels.

!   OPTIONAL INPUT:
!    PKAP         : kappa = R/cp (required for NHQE).

!   OPTIONAL OUTPUT:
!    PNHPREF      : total pressure "pre" at full levels.
!    PNHPPI       : ratio pre/prehyd at full levels (required for NHEE).
!    PRNHPPI      : ratio prehyd/pre at full levels (required for NHEE).
!    PDEP         : pressure departure "pre - prehyd" at full levels (required for NHEE).
!    PEQCHAF      : exp((R/cp) log(pre/prehyd)) at full levels (required for NHQE).
!    PEIQCHAF     : exp(-(R/cp) log(pre/prehyd)) at full levels (required for NHQE).

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
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Introduce NHQE model.
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   H Petithomme (Dec 2020): use of pointer
!------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK


! -----------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM) ,INTENT(IN)                     :: KPDVAR
TYPE(GEOMETRY)     ,INTENT(IN)                     :: YDGEOMETRY
INTEGER(KIND=JPIM) ,INTENT(IN)                     :: KPROMA
INTEGER(KIND=JPIM) ,INTENT(IN)                     :: KFLEV
INTEGER(KIND=JPIM) ,INTENT(IN)                     :: KSTART 
INTEGER(KIND=JPIM) ,INTENT(IN)                     :: KEND 
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PSPD(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PREF(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)  ,OPTIONAL          :: PKAP(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(OUT) ,OPTIONAL          :: PNHPREF(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(OUT) ,OPTIONAL ,TARGET  :: PNHPPI(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(OUT) ,OPTIONAL          :: PRNHPPI(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(OUT) ,OPTIONAL          :: PDEP(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(OUT) ,OPTIONAL          :: PEQCHAF(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(OUT) ,OPTIONAL          :: PEIQCHAF(KPROMA,KFLEV) 

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPRB),TARGET :: ZNHPPI0(KPROMA,KFLEV) 
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZNHPPI(:,:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

IF ( KPDVAR /= 2 ) RETURN

IF (LHOOK) CALL DR_HOOK('GNHPRE', 0, ZHOOK_HANDLE)

! -----------------------------------------------------------------------------

IF(PRESENT(PNHPPI)) THEN
  ZNHPPI => PNHPPI(:,:)
ELSE
  ZNHPPI => ZNHPPI0(:,:)
ENDIF

! valid for npdvar=2 only
DO JLEV=1,KFLEV
  DO JROF=KSTART,KEND
    ZNHPPI(JROF,JLEV)=EXP(PSPD(JROF,JLEV))
  ENDDO
ENDDO

! valid for any npdvar (but only npdvar=2 is treated here)
IF (PRESENT(PNHPREF)) THEN
  DO JLEV=1,KFLEV
    DO JROF=KSTART,KEND
      PNHPREF(JROF,JLEV)=ZNHPPI(JROF,JLEV)*PREF(JROF,JLEV)  
    ENDDO
  ENDDO
ENDIF

IF (PRESENT(PRNHPPI)) THEN
  DO JLEV=1,KFLEV
    DO JROF=KSTART,KEND
      PRNHPPI(JROF,JLEV)=1.0_JPRB/ZNHPPI(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

IF (PRESENT(PDEP)) THEN
  DO JLEV=1,KFLEV
    DO JROF=KSTART,KEND
      PDEP(JROF,JLEV)=( ZNHPPI(JROF,JLEV) - 1.0_JPRB )*PREF(JROF,JLEV)  
    ENDDO
  ENDDO
ENDIF

! valid for npdvar=2 only
IF (PRESENT(PKAP)) THEN
  IF (PRESENT(PEQCHAF)) THEN
    ! * useful for NHQE.
    DO JLEV=1,KFLEV
      DO JROF=KSTART,KEND
        PEQCHAF(JROF,JLEV)=EXP( PKAP(JROF,JLEV)*PSPD(JROF,JLEV) )
      ENDDO
    ENDDO
  ENDIF

  IF (PRESENT(PEIQCHAF)) THEN
    ! * useful for NHQE.
    IF (PRESENT(PEQCHAF)) THEN
      DO JLEV=1,KFLEV
        DO JROF=KSTART,KEND
          PEIQCHAF(JROF,JLEV)=1._JPRB/PEQCHAF(JROF,JLEV)
        ENDDO
      ENDDO
    ELSE
      DO JLEV=1,KFLEV
        DO JROF=KSTART,KEND
          PEIQCHAF(JROF,JLEV)=EXP(-PKAP(JROF,JLEV)*PSPD(JROF,JLEV))
        ENDDO
      ENDDO
    ENDIF
  ENDIF
ENDIF

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHPRE', 1, ZHOOK_HANDLE)

END SUBROUTINE GNHPRE

