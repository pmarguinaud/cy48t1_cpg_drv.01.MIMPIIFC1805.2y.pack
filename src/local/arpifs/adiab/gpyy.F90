SUBROUTINE GPYY(&
 ! --- INPUT -----------------------------------------------------------------
 & LDVERTFE, LDVFE_GW, YDGEOMETRY,KFLEV,KPROMA,KST,KEND,POROGL,POROGM,PUH,PVH,&
 ! --- OUTPUT -------------------------------------------------------
 & PNHY)

! GPYY  - compute Y-term in NHEE model
!          This routine is called for NVDVAR=5

! Purpose
! -------

!   Retained VFD discretisation is:
!    [Y]_[lbar] = -[S_[lbar] V_[lbar] grad[Phi_s]

! Remarks:

! Interface
! ---------
!   * INPUT:
!   YDGEOMETRY   : structure containing all geometry.
!   KFLEV        : number of levels.
!   KPROMA       : length of work.
!   KST          : start of work.
!   KEND         : end of work.
!   POROGL       : zonal component of "grad[Phi_s]".
!   POROGM       : meridian component of "grad[Phi_s]".
!   PUH          : U-wind at half levels.
!   PVH          : V-wind at half levels.

!   * OUTPUT
!   PNHY         : -S_[lbar] V_[lbar].grad[Phi_s]

! Externals
! ---------

! Method
! ------

! Reference
! ---------

! Author
! ------
!   F. Voitus (August 2018)

! Modifications
! -------------
! End Modifications
!---------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK


! -----------------------------------------------------------------------------

IMPLICIT NONE

LOGICAL            ,INTENT(IN)   :: LDVERTFE
LOGICAL            ,INTENT(IN)   :: LDVFE_GW
TYPE(GEOMETRY)     ,INTENT(IN)   :: YDGEOMETRY
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KFLEV
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KPROMA 
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KST
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KEND 
REAL(KIND=JPRB)    ,INTENT(IN)   :: POROGL(KPROMA) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: POROGM(KPROMA) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PUH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PVH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)    ,INTENT(OUT)  :: PNHY(KPROMA,0:KFLEV) 

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

#include "abor1.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPYY',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB)

! -----------------------------------------------------------------------------


IF( LDVERTFE .AND. LDVFE_GW )THEN

  ! LVFE_GW not yet coded
  ! Calculation of ZDGWINCR13 requires a VFE vertical derivative, which uses full-level values of (U,V).
  CALL ABOR1(' GPYY: case LVFE_GW=T is not yet coded ')

ELSE

  ! VFD discretisation.
  ! Remark ky: VFD discretisation is not consistent
  !  with LVFE_X_TERM=T discretisation of term NHXS (better to keep LVFE_X_TERM=F).

  DO JLEV=0,KFLEV
    DO JROF=KST,KEND
      PNHY(JROF,JLEV) = -YDVAB%VRATH(JLEV)*(PUH(JROF,JLEV)*POROGL(JROF)+PVH(JROF,JLEV)*POROGM(JROF)) 
    ENDDO
  ENDDO

ENDIF


! -----------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPYY',1,ZHOOK_HANDLE)

END SUBROUTINE GPYY

