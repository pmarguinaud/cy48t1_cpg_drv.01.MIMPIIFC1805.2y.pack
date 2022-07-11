!OCL  NOEVAL
SUBROUTINE GNHEE_LAPL(&
 ! --- INPUT ---------------------------------------------------------------------
 & LDVFE_LAPL_BC,YDGEOMETRY,KSTART,KPROF,PIN,PINS,&
 & PLNPR,PALPH,PREF,PREH,&
 ! --- OUTPUT --------------------------------------------------------------------
 & POUT)  

!**** *GNHEE_LAPL* - Computation of [LAPL Z] at full levels, for NHEE model.

!     Purpose.
!     --------

!     Computes:

!                d [ dZ/d(log(prehyd)) + Z ]   d [ d(prehyd Z)/d(prehyd) ]
!     [LAPL Z] = --------------------------- = ---------------------------
!                       d(log(prehyd))               d(log(prehyd))

!     For example Z may be [exp(Qcha) - 1]=(pre-prehyd)/prehyd
!     Top value of Z is assumed to be 0.

!     LAPL is the non-linear counterpart of operator LLsstar computed in SISEVE

!     VFD general expression of [LAPL Z](l) is:
!      [LAPL Z](l) = CA(l) [Z(l-1)-Z(l)] + CC(l) [Z(l+1)-Z(l)]
!     where:
!      CA(l)=prehyd(l-1)/[delta(l) (prehyd(l)-prehyd(l-1))
!      CC(l)=prehyd(l+1)/[delta(l) (prehyd(l+1)-prehyd(l))
!     Expressions are slightly different for l=1 and l=NFLEVG
!     General usage is Z=[exp(Qcha) - 1] at full levels.

!     At the surface we add ZS/delta(l=L) to [LAPL Z](l=L); general usage is ZS=[D w_surf/Dt]/g

!**   Interface.
!     ----------
!        *CALL* *GNHEE_LAPL(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!           YDGEOMETRY        : structure containing geometry
!           KSTART            : start of work.
!           KPROF             : working length.
!           PLNPR             : "delta" at full levels.
!           PALPH             : "alpha" at full levels.
!           PREF              : "prehyd" at full levels.
!           PREH              : "prehyd" at half levels.
!           PIN               : input quantity containing Z at full levels.
!           PINS              : input quantity containing surface condition ZS (matches with (1/G) [D w_surf/Dt]).

!         * OUTPUT:
!           POUT              : output quantity containing [LAPL Z] at full levels.

!        Implicit arguments :   None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    None.
!     ----------

!     Reference.
!     ----------
!        See documentations Arpege about the NHEE model.

!     Author.
!     -------
!        K. YESSAD and F. VOITUS.
!        Original : July 2018

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL,           INTENT(IN)    :: LDVFE_LAPL_BC
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLNPR(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PINS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POUT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,JLEV
REAL(KIND=JPRB) :: ZA(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZC(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

#include "abor1.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHEE_LAPL',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA,NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)

IF (LDVFE_LAPL_BC) CALL ABOR1(' GNHEE_LAPL: VFE Laplacian not coded!')

!     ------------------------------------------------------------------

!*       1.    COMPUTATION FOR l=1
!              -------------------

DO JROF=KSTART,KPROF
  ZA(JROF,1)=PREH(JROF,0)/(PLNPR(JROF,1)*(PREF(JROF,1)-PREH(JROF,0)))
  ZC(JROF,1)=PREF(JROF,2)/(PLNPR(JROF,1)*(PREF(JROF,2)-PREF(JROF,1)))
  POUT(JROF,1)=-ZA(JROF,1)*PIN(JROF,1)+ZC(JROF,1)*(PIN(JROF,2)-PIN(JROF,1))
ENDDO

!     ------------------------------------------------------------------

!*       2.    COMPUTATION FOR l=2 to NFLEVG-1
!              ------------------------------

DO JLEV=2,NFLEVG-1
  DO JROF=KSTART,KPROF
    ZA(JROF,JLEV)=PREF(JROF,JLEV-1)/(PLNPR(JROF,JLEV)*(PREF(JROF,JLEV)-PREF(JROF,JLEV-1)))
    ZC(JROF,JLEV)=PREF(JROF,JLEV+1)/(PLNPR(JROF,JLEV)*(PREF(JROF,JLEV+1)-PREF(JROF,JLEV)))
    POUT(JROF,JLEV)=ZA(JROF,JLEV)*(PIN(JROF,JLEV-1)-PIN(JROF,JLEV))+&
     & ZC(JROF,JLEV)*(PIN(JROF,JLEV+1)-PIN(JROF,JLEV))
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*       3.    COMPUTATION FOR l=NFLEVG
!              -----------------------

! PIN contribution:
DO JROF=KSTART,KPROF
  ZA(JROF,NFLEVG)=PREF(JROF,NFLEVG-1)/(PLNPR(JROF,NFLEVG)*(PREF(JROF,NFLEVG)-PREF(JROF,NFLEVG-1)))
  ZC(JROF,NFLEVG)=1._JPRB/PLNPR(JROF,NFLEVG)
  POUT(JROF,NFLEVG)=ZA(JROF,NFLEVG)*(PIN(JROF,NFLEVG-1)-PIN(JROF,NFLEVG))-ZC(JROF,NFLEVG)*PIN(JROF,NFLEVG)
ENDDO

! PINS contribution:
DO JROF=KSTART,KPROF
  POUT(JROF,NFLEVG)=POUT(JROF,NFLEVG)+PINS(JROF)/PLNPR(JROF,NFLEVG)
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHEE_LAPL',1,ZHOOK_HANDLE)
END SUBROUTINE GNHEE_LAPL
