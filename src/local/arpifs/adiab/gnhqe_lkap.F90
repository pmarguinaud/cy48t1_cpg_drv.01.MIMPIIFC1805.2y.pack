!OCL  NOEVAL
SUBROUTINE GNHQE_LKAP(&
 ! --- INPUT ---------------------------------------------------------------------
 & LDVFE_LAPL_BC,YDGEOMETRY,KSTART,KPROF,PIN,PINS,&
 & PLNPR,PALPH,PREF,PREH,PKAPF,&
 ! --- OUTPUT --------------------------------------------------------------------
 & POUT)  

!**** *GNHQE_LKAP* - Computation of [LKAP Z] at full levels, for NHQE model.

!     Purpose.
!     --------

!     Computes:

!                d [ dZ/d(log(prehyd)) + Kap Z ]   d [ d(prehyd Z)/d(prehyd) + (Kap -1) Z ]
!     [LKAP Z] = ------------------------------- = ----------------------------------------
!                       d(log(prehyd))                        d(log(prehyd))

!     where Kap=R/cp and EQ=exp(Kap*Qcha).
!     For example Z may be [exp(Kap*Qcha) - 1]

!     LKAP is the non-linear counterpart of operator LLsstar_kap computed in SILKAP

!     VFD general expression of [LKAP Z](l) is:
!      [LKAP Z](l) = CA(l) [1+(Kap(l)-1)(alpha(l)-delta(l))] [Z(l-1)-Z(l)] + CC(l) [1+(Kap(l)-1) alpha(l)] [Z(l+1)-Z(l)]
!     where:
!      CA(l)=prehyd(l-1)/[delta(l) (prehyd(l)-prehyd(l-1))
!      CC(l)=prehyd(l+1)/[delta(l) (prehyd(l+1)-prehyd(l))
!     Expressions are slightly different for l=1 and l=NFLEVG
!     General usage is Z=[exp(Kap*Qcha) - 1] at full levels.

!     At the surface we add ZS/delta(l=L) to [LKAP Z](l=L); general usage is ZS=[D w_surf/Dt]/g

!**   Interface.
!     ----------
!        *CALL* *GNHQE_LKAP(...)

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
!           PKAPF             : "Kap=R/Cp" at full levels.
!           PIN               : input quantity containing Z at full levels.
!           PINS              : input quantity containing surface condition ZS (matches with (1/G) [D w_surf/Dt]).

!         * OUTPUT:
!           POUT              : output quantity containing [LKAP Z] at full levels.

!        Implicit arguments :   None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    None.
!     ----------

!     Reference.
!     ----------
!        See documentations Arpege about the NH model.

!     Author.
!     -------
!        K. YESSAD and F. VOITUS.
!        Original : March 2018

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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PINS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POUT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,JLEV
REAL(KIND=JPRB) :: ZA(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZAA(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZC(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZCC(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

#include "abor1.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHQE_LKAP',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA,NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)

IF (LDVFE_LAPL_BC) CALL ABOR1(' GNHQE_LKAP: VFE Laplacian not coded!')

!     ------------------------------------------------------------------

!*       1.    COMPUTATION FOR l=1
!              -------------------

DO JROF=KSTART,KPROF
  ZA(JROF,1)=PREH(JROF,0)/(PLNPR(JROF,1)*(PREF(JROF,1)-PREH(JROF,0)))
  ZC(JROF,1)=PREF(JROF,2)/(PLNPR(JROF,1)*(PREF(JROF,2)/PLNPR(JROF,1)-PREF(JROF,1)*(PALPH(JROF,2)-PLNPR(JROF,2))))
  ZAA(JROF,1)=ZA(JROF,1)*(1.0_JPRB+(PKAPF(JROF,1)-1.0_JPRB)*(PALPH(JROF,1)-PLNPR(JROF,1)))
  ZCC(JROF,1)=ZC(JROF,1)*(1.0_JPRB+(PKAPF(JROF,1)-1.0_JPRB)*PALPH(JROF,1))
  POUT(JROF,1)=-ZAA(JROF,1)*PIN(JROF,1)+ZCC(JROF,1)*(PIN(JROF,2)-PIN(JROF,1))
ENDDO

!     ------------------------------------------------------------------

!*       2.    COMPUTATION FOR l=2 to NFLEVG-1
!              ------------------------------

DO JLEV=2,NFLEVG-1
  DO JROF=KSTART,KPROF
    ZA(JROF,JLEV)=PREF(JROF,JLEV-1)/(PLNPR(JROF,JLEV)*(PREF(JROF,JLEV)-PREF(JROF,JLEV-1)))
    ZC(JROF,JLEV)=PREF(JROF,JLEV+1)/(PLNPR(JROF,JLEV)*(PREF(JROF,JLEV+1)-PREF(JROF,JLEV)))
    ZAA(JROF,JLEV)=ZA(JROF,JLEV)*(1.0_JPRB+(PKAPF(JROF,JLEV)-1.0_JPRB)*(PALPH(JROF,JLEV)-PLNPR(JROF,JLEV)))
    ZCC(JROF,JLEV)=ZC(JROF,JLEV)*(1.0_JPRB+(PKAPF(JROF,JLEV)-1.0_JPRB)*PALPH(JROF,JLEV))
    POUT(JROF,JLEV)=ZAA(JROF,JLEV)*(PIN(JROF,JLEV-1)-PIN(JROF,JLEV))+&
     & ZCC(JROF,JLEV)*(PIN(JROF,JLEV+1)-PIN(JROF,JLEV))
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*       3.    COMPUTATION FOR l=NFLEVG
!              -----------------------

! PIN contribution:
DO JROF=KSTART,KPROF
  ZA(JROF,NFLEVG)=PREF(JROF,NFLEVG-1)/(PLNPR(JROF,NFLEVG)*(PREF(JROF,NFLEVG)-PREF(JROF,NFLEVG-1)))
  ZC(JROF,NFLEVG)=PKAPF(JROF,NFLEVG)*PREH(JROF,NFLEVG) &
   & /(PLNPR(JROF,NFLEVG)*(PKAPF(JROF,NFLEVG)*PREH(JROF,NFLEVG)-(PKAPF(JROF,NFLEVG)-1.0_JPRB)*PREF(JROF,NFLEVG)))
  ZAA(JROF,NFLEVG)=ZA(JROF,NFLEVG)*(1.0_JPRB+(PKAPF(JROF,NFLEVG)-1.0_JPRB)*(PALPH(JROF,NFLEVG)-PLNPR(JROF,NFLEVG)))
  ZCC(JROF,NFLEVG)=ZC(JROF,NFLEVG)*(1.0_JPRB+(PKAPF(JROF,NFLEVG)-1.0_JPRB)*PALPH(JROF,NFLEVG))
  POUT(JROF,NFLEVG)=ZAA(JROF,NFLEVG)*(PIN(JROF,NFLEVG-1)-PIN(JROF,NFLEVG))-ZCC(JROF,NFLEVG)*PIN(JROF,NFLEVG)
ENDDO

! PINS contribution:
DO JROF=KSTART,KPROF
  POUT(JROF,NFLEVG)=POUT(JROF,NFLEVG)+PKAPF(JROF,NFLEVG)*PINS(JROF)/PLNPR(JROF,NFLEVG)
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHQE_LKAP',1,ZHOOK_HANDLE)
END SUBROUTINE GNHQE_LKAP
