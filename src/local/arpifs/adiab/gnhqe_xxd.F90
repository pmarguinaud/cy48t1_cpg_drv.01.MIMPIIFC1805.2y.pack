SUBROUTINE GNHQE_XXD(   LDRUBC, LDVERTFE, YDGEOMETRY, KFLEV, KPROMA, KST, KEND, PU, PV, PLNPR, &
& PRDELP, PALPH, PRTGR, PSPL, PSPM, PKAP, PNHXD)

! GNHQE_XXD - Diagnose "divergence part" of the NHX-term which appears in the NHQE model.
!             NHX writes NHX=NHX_s+NHX_d
!             In the NHQE model, d4 = dver + NHX, but D3 = D + d4 - NHX_d = D + dver + NHX_s

! Purpose
! -------
!   Diagnose "divergence part" of the NHX-term which appears in the NHQE model.
!   NHX writes NHX = NHX_s + NHX_d
!   NHX_s has an expression close to the NHEE-model total NHX-term, and is linked to the 
!    interaction between wind shear and orography gradient.
!   In the NHQE model, d4 = dver + NHX, but D3=D+d4-NHX_d=D+dver+NHX_s:
!    we need a separate evaluation of NHX_d and NHX_s for the RHS of dver equation.
!   It involves the calculation of a vertical integral belonging to the following family:
!    (1/Pi) int_(eta'=0 to eta'=eta) (dB/deta') Z deta'


! Interface
! ---------
!   * INPUT:
!   YDGEOMETRY   : structure containing all geometry
!   KFLEV        : number of levels.
!   KPROMA       : length of work
!   KST          : start of work
!   KEND         : end of work
!   PU           : U component of the wind, at full levels
!   PV           : V component of the wind, at full levels
!   PXYB         : contains pressure depth, "delta", "alpha".
!   PSPL         : zonal component of "grad prehyds"
!   PSPM         : meridian component of "grad prehyds"
!   PKAP         : R/Cp at full levels.

!   * OUTPUT:
!   PNHXD        : NHX_d-term at full levels.

! Externals
! ---------

! Method
! ------

! Reference
! ---------

! Author
! ------
!   K. Yessad and F. Voitus (July 2017)

! Modifications
! -------------
!   H Petithomme (Dec 2020): use ZOUT instead of Z_VERINT, merge VFD loops
! End Modifications
!---------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK



! -----------------------------------------------------------------------------

IMPLICIT NONE

LOGICAL            ,INTENT(IN)   :: LDRUBC
LOGICAL            ,INTENT(IN)   :: LDVERTFE
TYPE(GEOMETRY)     ,INTENT(IN)   :: YDGEOMETRY
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KFLEV
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KPROMA 
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KST
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KEND 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PU(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PV(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PLNPR(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PRDELP(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PALPH(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PRTGR(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PSPL(KPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PSPM(KPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PKAP(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(OUT)  :: PNHXD(KPROMA,KFLEV) 

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPRB) :: ZVP(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZOUT(KPROMA,KFLEV+1)
REAL(KIND=JPRB) :: ZIN(KPROMA,0:KFLEV+1)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

#include "abor1.intfb.h"
#include "verdisint.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHQE_XXD', 0, ZHOOK_HANDLE)

! -----------------------------------------------------------------------------

ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA)

IF (LDRUBC) THEN
  ! Input [etadot (d prehyd / d eta)]_top is not always available there.
  CALL ABOR1(' GNHQE_XXD: LRUBC not coded')
ENDIF

! * Computes "vec(V) * grad prehyds" at full levels.
DO JLEV=1,KFLEV
  DO JROF=KST,KEND
    ZVP(JROF,JLEV)=PU(JROF,JLEV)*PSPL(JROF)+PV(JROF,JLEV)*PSPM(JROF)
  ENDDO
ENDDO

! * Term containing vertical integral
!   (1/Pi) int_(eta'=0 to eta'=eta) (dB/deta') [vec(V) * grad prehyds] deta'
IF( LDVERTFE )THEN
  ZIN(KST:KEND,0)=0.0_JPRB
  ZIN(KST:KEND,KFLEV+1)=0.0_JPRB
  DO JLEV=1,KFLEV
    DO JROF=KST,KEND
      ZIN(JROF,JLEV)=YDVAB%VDELB(JLEV)*ZVP(JROF,JLEV)*YDVETA%VFE_RDETAH(JLEV)
    ENDDO
  ENDDO
  CALL VERDISINT(YDGEOMETRY%YRVFE, 'ITOP', '11', KPROMA, KST, KEND, KFLEV, ZIN, ZOUT)

  ! note that LNPR*RDELP is exactly [1/Pi] at full levels in this case.
  DO JLEV=1,KFLEV
    DO JROF=KST,KEND
      ZOUT(JROF,JLEV)=(PLNPR(JROF,JLEV)*PRDELP(JROF,JLEV))*ZOUT(JROF,JLEV)
    ENDDO
  ENDDO
ELSE
  ZOUT(KST:KEND,1)=0.0_JPRB
  DO JLEV=1,KFLEV
    DO JROF=KST,KEND
      ZOUT(JROF,JLEV+1)=ZOUT(JROF,JLEV)+YDVAB%VDELB(JLEV)*ZVP(JROF,JLEV)
      ZOUT(JROF,JLEV)=(PLNPR(JROF,JLEV)*PRDELP(JROF,JLEV)*ZOUT(JROF,JLEV) &
       & + PALPH(JROF,JLEV)*PRDELP(JROF,JLEV)*YDVAB%VDELB(JLEV)*ZVP(JROF,JLEV))
    ENDDO
  ENDDO
ENDIF

! * NHX_d
DO JLEV=1,KFLEV
  DO JROF=KST,KEND
    PNHXD(JROF,JLEV)=(1.0_JPRB - PKAP(JROF,JLEV))*(PRTGR(JROF,JLEV)*ZVP(JROF,JLEV)-ZOUT(JROF,JLEV))
  ENDDO
ENDDO
END ASSOCIATE

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHQE_XXD', 1, ZHOOK_HANDLE)
END SUBROUTINE GNHQE_XXD

