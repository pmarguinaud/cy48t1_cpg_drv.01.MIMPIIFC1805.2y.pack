SUBROUTINE GNHEE_SVDINCR13(&
 ! --- INPUT -----------------------------------------------------------------
 & LDVERTFE, LDVFE_GW, YDGEOMETRY,KFLEV,KPROMA,KST,KEND,POROGL,POROGM,PUH,PVH,&
 ! --- OPTIONAL INPUT -------------------------------------------------------
 & PLNPR,PRT,PNHPPI,&
 ! --- OPTIONAL OUTPUT -------------------------------------------------------
 & PDGWINCR13,PSVDINCR13)

! GNHEE_SVDINCR13 - compute d13 - dver in NHEE model
!                   This routine is called for NVDVAR=13 or 14

! Purpose
! -------
!   Diagnose (d13 - dver) in NHEE model for options NVDVAR=13 or 14, at full levels.
!    d13 - dver = (pre/(prehyd RT)) [d(S V)/d log(prehyd)] grad[Phi_s]
!   That looks like the term NHXS but replacing [dV/d log(prehyd)] grad(Phi)
!    by [d(S V)/d log(prehyd)] grad[Phi_s].
!   Only the VFD discretisation of this term is coded (assumes LVFE_GW=F,
!    and not very consistent with LVFE_X_TERM=T).

!   Retained VFD discretisation is:
!    [d13 - dver]_[l] = [pre/prehyd]_[l] * [ 1/(R_[l] T_[l] delta_[l]) ] *
!     [S_[lbar] V_[lbar] - S_[lbar-1] V_[lbar-1]] grad[Phi_s]

! Remarks:
! - to avoid inconsistencies in parts doing conversions between d13 and [gW],
!   the increment (d13 - dver) is defined with the same value of R as 'dver',
!   according to the value of L_RDRY_VD (input array PRT).

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

!   * OPTIONAL INPUT (used to compute PSVDINCR13):
!   PLNPR        : "delta" at full levels.
!   PRT          : (R*Temperature) at full levels.
!   PNHPPI       : [pre/prehyd] at full levels.

!   * OPTIONAL OUTPUT:
!   PDGWINCR13   : (d(gw) - d(gW)) at full levels.
!   PSVDINCR13   : (d13 - dver) at full levels.

! Externals
! ---------

! Method
! ------

! Reference
! ---------

! Author
! ------
!   K. Yessad (August 2018)

! Modifications
! -------------
! End Modifications
!---------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK


! -----------------------------------------------------------------------------

IMPLICIT NONE

LOGICAL            ,INTENT(IN)             :: LDVERTFE
LOGICAL            ,INTENT(IN)             :: LDVFE_GW
TYPE(GEOMETRY)     ,INTENT(IN)             :: YDGEOMETRY
INTEGER(KIND=JPIM) ,INTENT(IN)             :: KFLEV
INTEGER(KIND=JPIM) ,INTENT(IN)             :: KPROMA 
INTEGER(KIND=JPIM) ,INTENT(IN)             :: KST
INTEGER(KIND=JPIM) ,INTENT(IN)             :: KEND 
REAL(KIND=JPRB)    ,INTENT(IN)             :: POROGL(KPROMA) 
REAL(KIND=JPRB)    ,INTENT(IN)             :: POROGM(KPROMA) 
REAL(KIND=JPRB)    ,INTENT(IN)             :: PUH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)             :: PVH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)  ,OPTIONAL  :: PLNPR(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)  ,OPTIONAL  :: PRT(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)  ,OPTIONAL  :: PNHPPI(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(OUT) ,OPTIONAL  :: PDGWINCR13(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(OUT) ,OPTIONAL  :: PSVDINCR13(KPROMA,KFLEV) 

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPRB) :: ZSVDINCR13(KPROMA,KFLEV) 
REAL(KIND=JPRB) :: ZDGWINCR13(KPROMA,KFLEV) 

LOGICAL :: LL_PSVDINCR13

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

#include "abor1.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHEE_SVDINCR13',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB)

! -----------------------------------------------------------------------------

LL_PSVDINCR13=PRESENT(PSVDINCR13).AND.PRESENT(PLNPR).AND.PRESENT(PRT).AND.PRESENT(PNHPPI)

IF( LDVERTFE .AND. LDVFE_GW )THEN

  ! LVFE_GW not yet coded
  ! Calculation of ZDGWINCR13 requires a VFE vertical derivative, which uses full-level values of (U,V).
  CALL ABOR1(' GNHEE_SVDINCR13: case LVFE_GW=T is not yet coded ')

ELSE

  ! VFD discretisation.
  ! Remark ky: VFD discretisation is not consistent
  !  with LVFE_X_TERM=T discretisation of term NHXS (better to keep LVFE_X_TERM=F).

  ! compute [Delta(S V)] grad[Phi_s] at full levels:
  DO JLEV=1,KFLEV
    DO JROF=KST,KEND
      ZDGWINCR13(JROF,JLEV) = &
       & YDVAB%VRATH(JLEV)*(PUH(JROF,JLEV)*POROGL(JROF)+PVH(JROF,JLEV)*POROGM(JROF)) &
       & - YDVAB%VRATH(JLEV-1)*(PUH(JROF,JLEV-1)*POROGL(JROF)+PVH(JROF,JLEV-1)*POROGM(JROF))
    ENDDO
  ENDDO

  ! multiply by (pre/prehyd)*(1/(RT [delta log(prehyd)])):
  IF (LL_PSVDINCR13) THEN
    DO JLEV=1,KFLEV
      DO JROF=KST,KEND
        ZSVDINCR13(JROF,JLEV)=ZDGWINCR13(JROF,JLEV)*PNHPPI(JROF,JLEV)/(PRT(JROF,JLEV)*PLNPR(JROF,JLEV))
      ENDDO
    ENDDO
  ENDIF

ENDIF

IF (PRESENT(PDGWINCR13)) THEN
  PDGWINCR13(KST:KEND,1:KFLEV)=ZDGWINCR13(KST:KEND,1:KFLEV)
ENDIF

IF (LL_PSVDINCR13) THEN
  PSVDINCR13(KST:KEND,1:KFLEV)=ZSVDINCR13(KST:KEND,1:KFLEV)
ENDIF


! -----------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHEE_SVDINCR13',1,ZHOOK_HANDLE)

END SUBROUTINE GNHEE_SVDINCR13

