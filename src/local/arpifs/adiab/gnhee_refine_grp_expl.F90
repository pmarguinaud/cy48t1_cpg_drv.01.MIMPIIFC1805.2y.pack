!OCL  NOEVAL
SUBROUTINE GNHEE_REFINE_GRP_EXPL(&
 ! --- INPUT -----------------------------------------------------------------
 & LDVERTFE, YDGEOMETRY,KST,KEND,&
 & PRT,PRTL,PRTM,PREL,PREM,PRTGR,PALPH,&
 & PRNHPPI,PQCHAL,PQCHAM,PHIHL,PHIHM,   &
 & PDEP,PREF,PREH,PWWH2F,&
 ! --- OUTPUT ----------------------------------------------------------------
 & PSGRTL,PSGRTM, &
 ! --- INPUT -----------------------------------------------------------------
 & PALPHPLL_DER, PALPHPLM_DER)

!**** *GNHEE_REFINE_GRP_EXPL* - Computation of the pressure gradient force term used in the
!                   RHS of the horizontal wind equation in NHEE model.
!                   This term is computed at full levels.

!     Purpose.
!     --------

!        The pressure gradient force term writes, for NHEE model:
!          - (d pre/d prehyd) grad[gz] - RT grad[pre]/pre

!**   Interface.
!     ----------
!        *CALL* *GNHEE_GRP(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!           KST         : start of work
!           KEND        : working length
!           PRT         : "R(air)*temperature" at full levels
!           PRTL        : zonal component of "grad (RT)" at full levels
!           PRTM        : meridian component of "grad (RT)" at full levels
!           PREL        : zonal component of "grad (prehyds)"
!           PREM        : meridian component of "grad (prehyds)"
!           PXYB        : contains pressure depth, "delta", "alpha".
!           PXYBDER     : contains grad(delta), grad(alpha), grad(alpha + log prehyd)
!           PRNHPPI     : (prehyd/pre) at full levels
!           PQCHAL      : zonal comp of "grad (log(pre/prehyd))" at full levels
!           PQCHAM      : merid comp of "grad (log(pre/prehyd))" at full levels
!           PHIHL       : zonal component of "grad[gz]" at half levels
!           PHIHM       : merid component of "grad[gz]" at half levels
!           PDEP        : (pre - prehyd) at full levels.
!           PREF        : "prehyd" at full levels.
!           PREH        : "prehyd" at half levels.
!           PWH2F       : interpolation coefficients at full levels.

!         * OUTPUT:
!           PSGRTL      : zonal comp. of pressure gradient force
!           PSGRTM      : merid comp. of pressure gradient force

!        Implicit arguments :   None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    None.
!     ----------

!     Reference.
!     ----------
!        See documentation Arpege

!     Author.
!     -------
!      Fabrice Voitus, December (2019)

!     Modifications.
!     --------------

!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK


!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL, INTENT (IN) :: LDVERTFE
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRTL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRTM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREL(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREM(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRTGR(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRNHPPI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQCHAL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQCHAM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHIHL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHIHM(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDEP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWWH2F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSGRTL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSGRTM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPHPLL_DER(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPHPLM_DER(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

REAL(KIND=JPRB)    :: ZDELPSDPREHYD(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    :: ZDEPH
INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPRB)    :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHEE_REFINE_GRP_EXPL',0,ZHOOK_HANDLE)

ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,NPROMA=>YDGEOMETRY%YRDIM%NPROMA)

!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF THE PRESSURE GRADIENT TERM.
!              ------------------------------------------
IF (LDVERTFE) THEN
  CALL ABOR1(' GNHEE_REFINE_GRP_EXPL 3.1.12: option not yet coded!')
ELSE
  ! * SOME PRE-CALCULATIONS:

  ! * Calculation of "d (pre-prehyd) / d prehyd" at half levels:
  DO JLEV=1,NFLEVG-1
    DO JROF=KST,KEND
      ZDELPSDPREHYD(JROF,JLEV)=(PDEP(JROF,JLEV+1)-PDEP(JROF,JLEV)) &
                         & /(PREF(JROF,JLEV+1)-PREF(JROF,JLEV))
    ENDDO
  ENDDO
  ! * Top of the model.
  !   [d (pre-prehyd)/d prehyd]_[l=0] is computed by
  !   [ (pre-prehyd)_[l=1]-(pre-prehyd)_[top] ]/[prehyd_[l=1]-prehyd_[top]]
  !   and (pre-prehyd)_[top] is assumed to be equal to zero.
  DO JROF=KST,KEND
    ZDELPSDPREHYD(JROF,0)=PDEP(JROF,1)/(PREF(JROF,1)-PREH(JROF,0))
  ENDDO
  ! * Surface.
  !   [d (pre-prehyd)/d prehyd]_[l=Lbar] is computed by
  !   [ (pre-prehyd)_[Lbar]-(pre-prehyd)_[L] ]/[prehyd_[Lbar] - prehyd_[L]]
  !   and [(pre-prehyd)]_[Lbar]=prehyd_[Lbar]*[(pre-prehyd)/ prehyd]_[L] is assumed.
  DO JROF=KST,KEND
    ZDELPSDPREHYD(JROF,NFLEVG)=PDEP(JROF,NFLEVG)*((PREH(JROF,NFLEVG)/PREF(JROF,NFLEVG))-1.0_JPRB)&
                            & /(PREH(JROF,NFLEVG)-PREF(JROF,NFLEVG))
  ENDDO


! * WE NOW COMPUTE THE PRESSURE GRADIENT TERM:
  ! * The pressure gradient term calculation requires some calculations on half levels.

  !   In the NH model, this term writes:
  !   (delta pre/delta prehyd) grap(Phi) + RT grad(pre) / pre
  !   For its discretisation we prefer to rewrite it as follows:
  !   [ grad(Phi) + (prehyd/pre) RT grad(log(prehyd)) ]
  !   + [ ((delta pre / delta prehyd) - 1) * grad(Phi) ]
  !   + [ RT grad(log(pre/prehyd)) ]
  !   + [ RT (1 - prehyd/pre) grad(log(prehyd)) ]
  !   The three last terms are non-zero only in the NH model.
  !   The first term [ grad(Phi) + (prehyd/pre) RT grad(log(prehyd)) ]
  !   looks like the hydrostatic pressure gradient term, but with an
  !   additional factor "prehyd/pre" and an additional term
  !   containing " (prehyd/pre) * grad(log(pre/prehyd)) ".
  !   This way of writing avoids the quasi cancellation of terms
  !   which have the same magnitude.

  ! * First get "grad[gz]" on half levels.
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      PSGRTL(JROF,JLEV)=PHIHL(JROF,JLEV)
      PSGRTM(JROF,JLEV)=PHIHM(JROF,JLEV)
    ENDDO
  ENDDO

  ! * Calculation of [ grad(Phi) + (prehyd/pre) RT grad(log(prehyd)) ]
  !   at full levels (cf. what is done in GPGRGEO to compute "grad(Phi)"
  !   at full levels, but replacing (alphl,alphm) by (alphpll,alphplm)):

  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      PSGRTL(JROF,JLEV)=PSGRTL(JROF,JLEV)&
       & +PRNHPPI(JROF,JLEV)*&
       & (PALPH(JROF,JLEV)*PRTL(JROF,JLEV)&
       & +PALPHPLL_DER(JROF,JLEV)*PRT(JROF,JLEV))
      PSGRTM(JROF,JLEV)=PSGRTM(JROF,JLEV)&
       & +PRNHPPI(JROF,JLEV)*&
       & (PALPH(JROF,JLEV)*PRTM(JROF,JLEV)&
       & +PALPHPLM_DER(JROF,JLEV)*PRT(JROF,JLEV))
    ENDDO
  ENDDO

  DO JLEV=1,NFLEVG 
    DO JROF=KST,KEND
      PSGRTL(JROF,JLEV)=PSGRTL(JROF,JLEV)&
       & -PRNHPPI(JROF,JLEV)*&
       & PALPH(JROF,JLEV)*PRT(JROF,JLEV)*PQCHAL(JROF,JLEV)
      PSGRTM(JROF,JLEV)=PSGRTM(JROF,JLEV)&
       & -PRNHPPI(JROF,JLEV)*&
       & PALPH(JROF,JLEV)*PRT(JROF,JLEV)*PQCHAM(JROF,JLEV)
    ENDDO
  ENDDO

  ! * Calculation of the additional anhydrostatic contributions at full levels:

  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      PSGRTL(JROF,JLEV)=PSGRTL(JROF,JLEV)&
       & + PWWH2F(JROF,JLEV,1)*ZDELPSDPREHYD(JROF,JLEV-1)*PHIHL(JROF,JLEV-1) &
       & + PWWH2F(JROF,JLEV,2)*ZDELPSDPREHYD(JROF,JLEV)*PHIHL(JROF,JLEV) &
       & + PRT(JROF,JLEV)*PQCHAL(JROF,JLEV) &
       & + PRT(JROF,JLEV)*(1.0_JPRB-PRNHPPI(JROF,JLEV))*PRTGR(JROF,JLEV) &
       & *PREL(JROF)
      PSGRTM(JROF,JLEV)=PSGRTM(JROF,JLEV)&
       & + PWWH2F(JROF,JLEV,1)*ZDELPSDPREHYD(JROF,JLEV-1)*PHIHM(JROF,JLEV-1) &
       & + PWWH2F(JROF,JLEV,2)*ZDELPSDPREHYD(JROF,JLEV)*PHIHM(JROF,JLEV) &
       & + PRT(JROF,JLEV)*PQCHAM(JROF,JLEV) &
       & + PRT(JROF,JLEV)*(1.0_JPRB-PRNHPPI(JROF,JLEV))*PRTGR(JROF,JLEV) &
       & *PREM(JROF)
    ENDDO
  ENDDO

ENDIF
!     ------------------------------------------------------------------

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('GNHEE_REFINE_GRP_EXPL',1,ZHOOK_HANDLE)

END SUBROUTINE GNHEE_REFINE_GRP_EXPL
