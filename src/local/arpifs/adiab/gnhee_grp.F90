!OCL  NOEVAL
SUBROUTINE GNHEE_GRP(&
 ! --- INPUT -----------------------------------------------------------------
 & LDVERTFE, YDGEOMETRY,KST,KEND,&
 & PRT,PRTL,PRTM,PREL,PREM,PRDELP,PRTGR,PALPH,&
 & PDELNHPRE,PLNNHPREFL,PLNNHPREFM,&
 & PRNHPPI,PQCHAL,PQCHAM,&
 & PHIHL,PHIHM,PHIFL,PHIFM,&
 ! --- OUTPUT ----------------------------------------------------------------
 & PSGRTL,PSGRTM, &
 & PALPHPLL_DER, PALPHPLM_DER)

!**** *GNHEE_GRP* - Computation of the pressure gradient force term used in the
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
!           PDELNHPRE   : [Delta pre] at full levels
!           PLNNHPREFL  : zonal component of "[grad pre]/pre" at full levels
!           PLNNHPREFM  : merid component of "[grad pre]/pre" at full levels
!           PRNHPPI     : (prehyd/pre) at full levels
!           PQCHAL      : zonal comp of "grad (log(pre/prehyd))" at full levels
!           PQCHAM      : merid comp of "grad (log(pre/prehyd))" at full levels
!           PHIHL       : zonal component of "grad[gz]" at half levels
!           PHIHM       : merid component of "grad[gz]" at half levels
!           PHIFL       : zonal component of "grad[gz]" at full levels
!           PHIFM       : merid component of "grad[gz]" at full levels

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
!      Karim YESSAD.
!      Original : DEC 2004 (in GPGRP); MAR 2017 (in GNHEE_GRP)

!     Modifications.
!     --------------
!      K. Yessad (Feb 2018): remove deep-layer formulations.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK


!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL            ,INTENT(IN)   :: LDVERTFE
TYPE(GEOMETRY)     ,INTENT(IN)   :: YDGEOMETRY
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KST 
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KEND 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PRT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PRTL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PRTM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PREL(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PREM(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PRDELP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PRTGR(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PALPH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PDELNHPRE(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PLNNHPREFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PLNNHPREFM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PRNHPPI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PQCHAL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PQCHAM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PHIHL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PHIHM(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PHIFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PHIFM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(OUT)  :: PSGRTL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(OUT)  :: PSGRTM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PALPHPLL_DER(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PALPHPLM_DER(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZDPRESDPREHYD(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHEE_GRP',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,NPROMA=>YDGEOMETRY%YRDIM%NPROMA)

!     ------------------------------------------------------------------


!*       1.    COMPUTATION OF THE PRESSURE GRADIENT TERM.
!              ------------------------------------------

! * SOME PRE-CALCULATIONS:

! * Calculation of "d pre / d prehyd" at full levels:
DO JLEV=1,NFLEVG
  DO JROF=KST,KEND
    ZDPRESDPREHYD(JROF,JLEV)=PDELNHPRE(JROF,JLEV)*PRDELP(JROF,JLEV)
  ENDDO
ENDDO

! * WE NOW COMPUTE THE PRESSURE GRADIENT TERM:

IF (LDVERTFE) THEN

  ! * The pressure gradient term is directly computed at full levels in
  !   this case and no quantity at half layers is computed.

  ! * First transfer " (d pre/d prehyd) grad[gz] + RT [grad pre]/pre " in (PSGRTL,PSGRTM).
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      PSGRTL(JROF,JLEV)=ZDPRESDPREHYD(JROF,JLEV)*PHIFL(JROF,JLEV) &
       & +PRT(JROF,JLEV)*PLNNHPREFL(JROF,JLEV)
      PSGRTM(JROF,JLEV)=ZDPRESDPREHYD(JROF,JLEV)*PHIFM(JROF,JLEV) &
       & +PRT(JROF,JLEV)*PLNNHPREFM(JROF,JLEV)
    ENDDO
  ENDDO

ELSE

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
       & + (ZDPRESDPREHYD(JROF,JLEV)-1.0_JPRB)*PHIFL(JROF,JLEV) &
       & + PRT(JROF,JLEV)*PQCHAL(JROF,JLEV) &
       & + PRT(JROF,JLEV)*(1.0_JPRB-PRNHPPI(JROF,JLEV))*PRTGR(JROF,JLEV) &
       & *PREL(JROF)
      PSGRTM(JROF,JLEV)=PSGRTM(JROF,JLEV)&
       & + (ZDPRESDPREHYD(JROF,JLEV)-1.0_JPRB)*PHIFM(JROF,JLEV) &
       & + PRT(JROF,JLEV)*PQCHAM(JROF,JLEV) &
       & + PRT(JROF,JLEV)*(1.0_JPRB-PRNHPPI(JROF,JLEV))*PRTGR(JROF,JLEV) &
       & *PREM(JROF)
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHEE_GRP',1,ZHOOK_HANDLE)
END SUBROUTINE GNHEE_GRP
