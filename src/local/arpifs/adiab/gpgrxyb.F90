SUBROUTINE GPGRXYB(KPROMA,KD,KF,KFLEV,LDCOEF,YDVAB,PREL,PREM,PXYB,PXYBDER, &
                 & PDELP, PLNPR, PRDELP, PALPH, PRTGR, PRPRE, PRPP)

!**** *GPGRXYB* - Complement to routine "GPXYB".
!                 Computation of the horizontal gradient of quantities
!                 "alpha" and "delta" at model levels.

!     Purpose.
!     --------

!     "alpha" and "delta" are computed at model levels in routine "GPXYB",
!     but not their horizontal gradient. So this routine provides the
!     horizontal gradients at full levels. Quantity
!     "(grad(alpha)) + (grad(prehyd)/prehyd)"
!     is also provided separately (for case LVERTFE=.F.
!     its "NDLNPR=0" expression is simpler than the expressions
!     of "grad(alpha)" and "grad(prehyd)/prehyd").
!     Discretisation depends on variables "NDLNPR" and "LVERTFE".

!     For LVERTFE=.F., NDLNPR=0, discretisations are:

!      (grad(delta))[l] =
!      - (A[lbar]*B[lbar-1]-A[lbar-1]*B[lbar])/(prehyd[lbar]*prehyd[lbar-1])
!      * (grad prehyds)

!      (grad(alpha))[l] + (grad(prehyd)/prehyd)[l] =
!      B[lbar]/prehyd[lbar] * (grad prehyds)

!      Quantity "(grad(alpha))[l]" is computed by substracting
!      "(grad(prehyd)/prehyd)[l]" from
!      "(grad(alpha))[l] + (grad(prehyd)/prehyd)[l]"

!     For LVERTFE=.F., NDLNPR=1 or 2, discretisations are:

!      (grad(delta))[l] =
!      - delta[l] * (A[lbar]*B[lbar-1]-A[lbar-1]*B[lbar])
!      * (1/sqrt(prehyd[lbar]*prehyd[lbar-1])) * (1/(delta prehyd[l]))
!      * (grad prehyds)

!      (grad(alpha))[l] =
!      - alpha[l] * (A[lbar]*B[lbar-1]-A[lbar-1]*B[lbar])
!      * (1/sqrt(prehyd[lbar]*prehyd[lbar-1])) * (1/(delta prehyd[l]))
!      * (grad prehyds)

!      (grad(prehyd)/prehyd)[l] = prtgr[l] * (grad prehyds)
!      where "prtgr[l]" is computed in routine "gpxyb" as:
!      prtgr[l] = { (delta B)[l]
!      + delta[l] * (A[lbar]*B[lbar-1]-A[lbar-1]*B[lbar])/(delta prehyd[l]) }
!      * { 1/(delta prehyd[l]) }

!      In this case "(grad(alpha))[l]" is computed prior to
!      "(grad(alpha))[l] + (grad(prehyd)/prehyd)[l]"

!     For LVERTFE=.T., NDLNPR=0, discretisations are:

!      (grad(delta))[l] =
!      delta[l] * ((Delta B)[l]/(Delta prehyd)[l] - B[l]/prehyd[l])
!      * (grad prehyds)

!      grad(alpha) is useless in this case.

!     Notations:
!      - "grad" is the horizontal gradient operator.
!        (grad X = vnabla X = M vnabla' X)
!      - "prehyd" is the hydrostatic pressure.
!      - "prehyds" is the surface hydrostatic pressure.

!**   Interface.
!     ----------
!        *CALL* *GPGRXYB(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!           KPROMA       : horizontal dimension
!           KD           : start of work
!           KF           : working length
!           KFLEV        : number of levels
!           LDCOEF       : if T, stores ZCOEFD, ZCOEFA, ZCOEFAPL in PXYBDER.
!           YDVAB        : contains information about hybrid vertical coordinate
!           PREL         : zonal component of "grad prehyds"
!           PREM         : meridian component of "grad prehyds"
!           PXYB         : contains pressure depth, "delta", "alpha".

!         * OUTPUT:
!           PXYBDER      : contains grad(delta), grad(alpha), grad(alpha + log prehyd)

!        Implicit arguments :   None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    None.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        K. YESSAD
!        Original : 00-08-11

!     Modifications.
!     --------------
!        K. Yessad (Dec 2008): remove dummy CDLOCK
!        K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!        K. Yessad (Dec 2011): use YDVAB.
!        K. Yessad (June 2017): introduce NDLNPR=2 (for NHQE model).
!        H. Petithomme (Dec 2020): use of pointers for optimisation
!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK
USE YOMDYNA   , ONLY : NDLNPR
USE YOMCVER   , ONLY : LVERTFE 
USE YOMVERT   , ONLY : TVAB
USE INTDYN_MOD, ONLY : YYTXYB, YYTXYBDER

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KD 
INTEGER(KIND=JPIM),INTENT(IN)    :: KF 
LOGICAL           ,INTENT(IN)    :: LDCOEF
TYPE(TVAB)        ,INTENT(IN)    :: YDVAB
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREM(KPROMA) 
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(IN) :: PXYB(KPROMA,KFLEV,YYTXYB%NDIM) 
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(OUT):: PXYBDER(KPROMA,KFLEV,YYTXYBDER%NDIM) 
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(IN) :: PDELP(KPROMA,KFLEV), PLNPR(KPROMA,KFLEV), PRDELP(KPROMA,KFLEV), PALPH(KPROMA,KFLEV)
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(IN) :: PRTGR(KPROMA,KFLEV), PRPRE(KPROMA,KFLEV), PRPP(KPROMA,KFLEV)

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPRB),TARGET :: ZCOEFA0(KPROMA),ZCOEFAPL0(KPROMA),ZCOEFD0(KPROMA)
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZCOEFA(:),ZCOEFAPL(:),ZCOEFD(:)
REAL(KIND=JPRB), POINTER :: ZDELP (:,:), ZLNPR (:,:), ZRDELP (:,:), ZALPH (:,:), ZRTGR (:,:), ZRPRE (:,:), ZRPP (:,:)  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPGRXYB',0,ZHOOK_HANDLE)

IF (PRESENT (PXYB)) THEN
  ZDELP  => PXYB (:,:,YYTXYB%M_DELP )
  ZLNPR  => PXYB (:,:,YYTXYB%M_LNPR )
  ZRDELP => PXYB (:,:,YYTXYB%M_RDELP)
  ZALPH  => PXYB (:,:,YYTXYB%M_ALPH )
  ZRTGR  => PXYB (:,:,YYTXYB%M_RTGR )
  ZRPRE  => PXYB (:,:,YYTXYB%M_RPRE )
  ZRPP   => PXYB (:,:,YYTXYB%M_RPP  )
ELSEIF (PRESENT (PDELP) .AND.  PRESENT (PLNPR) .AND.  PRESENT (PRDELP) .AND.  PRESENT (PALPH) .AND. &
     &  PRESENT (PRTGR) .AND.  PRESENT (PRPRE) .AND.  PRESENT (PRPP)) THEN
  ZDELP  => PDELP 
  ZLNPR  => PLNPR 
  ZRDELP => PRDELP
  ZALPH  => PALPH 
  ZRTGR  => PRTGR 
  ZRPRE  => PRPRE 
  ZRPP   => PRPP  
ELSE
  CALL ABOR1 ('GPGRXYB: PXYB OR PDELP/PLNPR/PRDELP/PALPH/PRTGR/PRPRE/PRPP REQUIRED')
ENDIF

!     ------------------------------------------------------------------

!*    1/ Calculation of "grad delta" at full levels.

IF (.NOT.LDCOEF) ZCOEFD => ZCOEFD0(:)

! optim: compilers may report dependence between pxyb and ydvab, ignored with ivdep/nodep

IF(LVERTFE) THEN
  DO JLEV=1,KFLEV
    IF (LDCOEF) ZCOEFD => PXYBDER(:,JLEV,YYTXYBDER%M_COEFD)
    !DIR$ IVDEP
    !CDIR NODEP
    DO JROF=KD,KF
      ZCOEFD(JROF)=(YDVAB%VDELB(JLEV)*ZRDELP(JROF,JLEV)&
       & -ZRTGR(JROF,JLEV))*ZLNPR(JROF,JLEV)  
      PXYBDER(JROF,JLEV,YYTXYBDER%M_LNPRL)=ZCOEFD(JROF)*PREL(JROF)
      PXYBDER(JROF,JLEV,YYTXYBDER%M_LNPRM)=ZCOEFD(JROF)*PREM(JROF)
    ENDDO
  ENDDO
ELSE
  DO JLEV=1,KFLEV
    IF (LDCOEF) ZCOEFD => PXYBDER(:,JLEV,YYTXYBDER%M_COEFD)
    !DIR$ IVDEP
    !CDIR NODEP
    DO JROF=KD,KF
      ZCOEFD(JROF)=-YDVAB%VC(JLEV)*ZRPP(JROF,JLEV)
      PXYBDER(JROF,JLEV,YYTXYBDER%M_LNPRL)=ZCOEFD(JROF)*PREL(JROF)
      PXYBDER(JROF,JLEV,YYTXYBDER%M_LNPRM)=ZCOEFD(JROF)*PREM(JROF)
    ENDDO
  ENDDO

!*    2/ Calculation of "grad (alpha + log prehyd)" at full levels
!        discretised as "grad alpha + (grad prehyd) / prehyd ",
!        and calculation of "grad alpha" at full levels.

  IF (.NOT.LDCOEF) THEN
   ZCOEFA => ZCOEFA0(:)
   ZCOEFAPL => ZCOEFAPL0(:)
  ENDIF

  IF(NDLNPR == 0) THEN
    DO JLEV=1,KFLEV
      IF (LDCOEF) THEN
        ZCOEFA => PXYBDER(:,JLEV,YYTXYBDER%M_COEFA)
        ZCOEFAPL => PXYBDER(:,JLEV,YYTXYBDER%M_COEFAPL)
      END IF

      !DIR$ IVDEP
      !CDIR NODEP
      DO JROF=KD,KF
        ZCOEFAPL(JROF)=YDVAB%VBH(JLEV)*ZRPRE(JROF,JLEV)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHPLL)=ZCOEFAPL(JROF)*PREL(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHPLM)=ZCOEFAPL(JROF)*PREM(JROF)

        ZCOEFA(JROF)=ZCOEFAPL(JROF)-ZRTGR(JROF,JLEV)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHL)=ZCOEFA(JROF)*PREL(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHM)=ZCOEFA(JROF)*PREM(JROF)
      ENDDO
    ENDDO
  ELSEIF(NDLNPR == 1 .OR. NDLNPR == 2) THEN
    DO JLEV=1,KFLEV
      IF (LDCOEF) THEN
        ZCOEFA => PXYBDER(:,JLEV,YYTXYBDER%M_COEFA)
        ZCOEFAPL => PXYBDER(:,JLEV,YYTXYBDER%M_COEFAPL)
      END IF

      !DIR$ IVDEP
      !CDIR NODEP
      DO JROF=KD,KF
        ZCOEFA(JROF)=-YDVAB%VC(JLEV)*ZRPP(JROF,JLEV)&
          & *ZALPH(JROF,JLEV)/ZLNPR(JROF,JLEV)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHL)=ZCOEFA(JROF)*PREL(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHM)=ZCOEFA(JROF)*PREM(JROF)

        ZCOEFAPL(JROF)=ZCOEFA(JROF)+ZRTGR(JROF,JLEV)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHPLL)=ZCOEFAPL(JROF)*PREL(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHPLM)=ZCOEFAPL(JROF)*PREM(JROF)
      ENDDO
    ENDDO
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('GPGRXYB',1,ZHOOK_HANDLE)
END SUBROUTINE GPGRXYB
