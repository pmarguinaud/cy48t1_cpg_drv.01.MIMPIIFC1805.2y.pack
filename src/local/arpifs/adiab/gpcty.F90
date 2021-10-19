!OCL  NOEVAL
SUBROUTINE GPCTY(YDVFE,KPROMA,KD,KF,KFLEV,LDRUBC,YDVAB,YDVETA,&
 & PU,PV,PD,PEVT,PXYB,PSPL,PSPM,PRPREF,PCTY,PDPHYCTY,&
 & PDELP, PLNPR, PRDELP, PALPH, PRTGR, PRPRE, PRPP)

!**** *GPCTY* - Computes vertical velocities.

!     Purpose.
!     --------

!      Computes the vertical velocities
!       "etadot (d prehyd / d eta)" at half levels
!       "(omega / prehyd)" at full levels

!      ----- The following discretisations are valid for lvertfe=.false. -----

!      Omitting the "delta m=1" flux precipitation terms, which will be later
!      added in routine "cpmvvps", the discretised expression of
!      "etadot (d prehyd / d eta)" on the half level number "lbar" is:

!       etadot (d prehyd / d eta) [lbar] =
!        B[lbar] * { sum[k=1 to L]
!        (Delta B)[k] * vec(V)[k] * (grad prehyds) }
!        + B[lbar] * { sum[k=1 to L] (Delta prehyd)[k] * (grad vec(V)[k]) }
!        - sum[k=1 to l] (Delta B)[k] * vec(V)[k] * (grad prehyds)
!        - sum[k=1 to l] (Delta prehyd)[k] * (grad vec(V)[k])
!        + { 1 - B[lbar] } * (etadot (d prehyd / d eta))[top]

!      where:
!       - "vec(V)" is the horizontal wind.
!       - "prehyds" is the surface hydrostatic pressure.
!       - "grad" is the "horizontal gradient" operator:
!         grad X = vnabla X = M vnabla' X

!      Omitting the "delta m=1" flux precipitation terms, which will be later
!      added in routine "cpmvvps", the discretised expression of
!      "(omega / prehyd)" on the full level number "l" is:

!       (omega / prehyd)[l] =
!        vec(V)[l] (grad prehyd / prehyd)[l]
!        - delta[l]/(Delta prehyd)[l] * { sum[k=1 to l-1]
!        (Delta B)[k] * vec(V)[k] * (grad prehyds) }
!        - delta[l]/(Delta prehyd)[l] * { sum[k=1 to l-1]
!        (Delta prehyd)[k] * (grad vec(V)[k]) }
!        - alpha[l]/(Delta prehyd)[l] * (Delta B)[l] * vec(V)[l] * (grad prehyds)
!        - alpha[l] * (grad vec(V)[l])
!        + delta[l]/(Delta prehyd)[l] * (etadot (d prehyd / d eta))[top]

!      This routine stores additional quantities:

!      * vertical integral of divergence without the "lrubc" contribution:
!        for the half level number "lbar" its discretised expression is:

!        psdiv[lbar] = sum[k=1 to l]
!        { (Delta B)[k] * vec(V)[k] * (grad prehyds)
!        + (Delta prehyd)[k] * (grad vec(V)[k]) }

!      * vertical integral of divergence with the "lrubc" contribution:
!        for the half level number "lbar" its discretised expression is:

!        psdvbc[lbar] = sum[k=1 to l]
!        { (Delta B)[k] * vec(V)[k] * (grad prehyds)
!        + (Delta prehyd)[k] * (grad vec(V)[k]) }
!        - (etadot (d prehyd / d eta))[top]

!      * divergence term
!        "pdivdp=grad( vec(V) * (Delta prehyd) )"
!        at full levels: for the full level number "l" its discretised expression is: 

!        pdivdp[l] = 
!        (Delta B)[l] * vec(V)[l] * (grad prehyds)
!        + (Delta prehyd)[l] * (grad vec(V)[l])

!      -----------------------------------------------------------------------

!**   Interface.
!     ----------
!        *CALL* *GPCTY(...)

!        Explicit arguments :
!        --------------------
!        * INPUT:
!          KPROMA       : horizontal dimensioning
!          KD           : start of work
!          KF           : depth of work
!          KFLEV        : number of layers
!          LDRUBC       : upper boundary condition switch
!          YDVAB        : contains information about hybrid vertical coordinate
!          YDVETA       : contains information about hybrid vertical coordinate: "eta"
!          PU           : U component of the wind, at full levels
!          PV           : V component of the wind, at full levels
!          PD           : horizontal divergence of the wind, at full levels
!          PEVT         : "etadot (d prehyd / d eta)" at the top
!          PXYB         : contains pressure depth, "delta", "alpha".
!          PSPL         : zonal component of "grad prehyds"
!          PSPM         : meridian component of "grad prehyds"
!          PRPREF       : inverse of full level pressures.

!        * OUTPUT:
!          PCTY         : contains vertical velocities, vertical integral of divergence.

!        * OPTIONAL INPUT:
!          PDPHYCTY     : mass source/sink from physics (previous dt in LAG physics)

!        Implicit arguments :   None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 88-02-04

!     Modifications.
!     --------------
!   K. Yessad (Jan 2008): complete (LVERCOR,LVERTFE)=(T,T).
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   K. Yessad (Dec 2011): use YDVAB.
!   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   H Petithomme (Dec 2020): use of pointer, add directives for dependencies
!     ------------------------------------------------------------------

USE YOMVERT   , ONLY : TVFE, TVAB, TVETA
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK
USE YOMCVER   , ONLY : LVERTFE
USE INTDYN_MOD, ONLY : YYTCTY, YYTXYB

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVFE)        ,INTENT(IN)    :: YDVFE
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KD
INTEGER(KIND=JPIM),INTENT(IN)    :: KF
LOGICAL           ,INTENT(IN)    :: LDRUBC
TYPE(TVAB)        ,INTENT(IN)    :: YDVAB
TYPE(TVETA)       ,INTENT(IN)    :: YDVETA
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PD(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVT(KPROMA)
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(IN) :: PXYB(KPROMA,KFLEV,YYTXYB%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPL(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPM(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPREF(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,TARGET,INTENT(OUT)   :: PCTY(KPROMA,0:KFLEV,YYTCTY%NDIM)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN) :: PDPHYCTY(KPROMA,KFLEV)
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(IN) :: PDELP(KPROMA,KFLEV), PLNPR(KPROMA,KFLEV), PRDELP(KPROMA,KFLEV), PALPH(KPROMA,KFLEV)
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(IN) :: PRTGR(KPROMA,KFLEV), PRPRE(KPROMA,KFLEV), PRPP(KPROMA,KFLEV)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSDIV(KPROMA,0:KFLEV+1)
REAL(KIND=JPRB) :: ZPSDIVFE(KPROMA,KFLEV+1)
REAL(KIND=JPRB) :: ZVP(KPROMA,KFLEV)
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZEVEL(:,:),ZVVEL(:,:),ZPSDIV(:,:),ZPSDVBC(:,:),ZDIVDP(:,:)
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZDELP (:,:), ZLNPR (:,:), ZRDELP (:,:), ZALPH (:,:), ZRTGR (:,:), ZRPRE (:,:), ZRPP (:,:)  

INTEGER(KIND=JPIM) :: JROF,JLEV
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPCTY',0,ZHOOK_HANDLE)

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
  CALL ABOR1 ('GPCTY: PXYB OR PDELP/PLNPR/PRDELP/PALPH/PRTGR/PRPRE/PRPP REQUIRED')
ENDIF



!     ------------------------------------------------------------------

!     check for non-supported configurations

IF (LDRUBC.AND.LVERTFE) CALL ABOR1('GPCTY: BAD OPTIONS')

!     ------------------------------------------------------------------

!*    1. Computes "vec(V) * grad prehyds" at full levels.
!     ---------------------------------------------------

DO JLEV=1,KFLEV
  ZVP(KD:KF,JLEV)=PU(KD:KF,JLEV)*PSPL(KD:KF)+PV(KD:KF,JLEV)*PSPM(KD:KF)
ENDDO

!     ------------------------------------------------------------------

!*    2. Sum divergence.
!     ------------------

ZEVEL(1:KPROMA,0:KFLEV) => PCTY(:,:,YYTCTY%M_EVEL)
ZVVEL(1:KPROMA,0:KFLEV) => PCTY(:,:,YYTCTY%M_VVEL)
ZPSDIV(1:KPROMA,0:KFLEV) => PCTY(:,:,YYTCTY%M_PSDIV)
ZPSDVBC(1:KPROMA,0:KFLEV) => PCTY(:,:,YYTCTY%M_PSDVBC)
ZDIVDP(1:KPROMA,0:KFLEV) => PCTY(:,:,YYTCTY%M_DIVDP)

IF (PRESENT(PDPHYCTY)) THEN
  DO JLEV=1,KFLEV
    ZDIVDP(KD:KF,JLEV)=PD(KD:KF,JLEV)*ZDELP(KD:KF,JLEV)&
     & +YDVAB%VDELB(JLEV)*ZVP(KD:KF,JLEV)&
     & -PDPHYCTY(KD:KF,JLEV)*ZDELP(KD:KF,JLEV)
  ENDDO
ELSE
  DO JLEV=1,KFLEV
    ZDIVDP(KD:KF,JLEV)=PD(KD:KF,JLEV)*ZDELP(KD:KF,JLEV)&
     & +YDVAB%VDELB(JLEV)*ZVP(KD:KF,JLEV)  
  ENDDO
ENDIF

ZPSDIV(KD:KF,0)=0.0_JPRB

IF(LVERTFE) THEN
  DO JLEV=1,KFLEV
    ZSDIV(KD:KF,JLEV)=ZDIVDP(KD:KF,JLEV)*YDVETA%VFE_RDETAH(JLEV)
  ENDDO

  ZSDIV(KD:KF,0)=0.0_JPRB
  ZSDIV(KD:KF,KFLEV+1)=0.0_JPRB
  CALL VERDISINT(YDVFE,'ITOP','11',KPROMA,KD,KF,KFLEV,ZSDIV,ZPSDIVFE)

  DO JLEV=1,KFLEV
    DO JROF=KD,KF
      ZPSDIV(JROF,JLEV)=ZPSDIVFE(JROF,JLEV)
    ENDDO
  ENDDO
ELSE
  DO JLEV=1,KFLEV
    !dir$ ivdep
    DO JROF=KD,KF
      ZPSDIV(JROF,JLEV)=ZPSDIV(JROF,JLEV-1)+ZDIVDP(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

IF (LDRUBC) THEN
  IF(LVERTFE) THEN
    DO JROF=KD,KF
      ZPSDVBC(JROF,KFLEV)=ZPSDIVFE(JROF,KFLEV+1)-PEVT(JROF)
    ENDDO
  ELSE
    DO JLEV=0,KFLEV
      !dir$ ivdep
      DO JROF=KD,KF
        ZPSDVBC(JROF,JLEV)=ZPSDIV(JROF,JLEV)-PEVT(JROF)
      ENDDO
    ENDDO
  ENDIF
ELSE
  IF(LVERTFE) THEN
    DO JROF=KD,KF
      ZPSDVBC(JROF,KFLEV)=ZPSDIVFE(JROF,KFLEV+1)
    ENDDO
  ELSE
    DO JLEV=0,KFLEV
      !dir$ ivdep
      DO JROF=KD,KF
        ZPSDVBC(JROF,JLEV)=ZPSDIV(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*    3. Computes "etadot (d prehyd / d eta)".
!     ----------------------------------------

! "etadot (d prehyd / d eta)" is computed at full levels if LVERTFE=T,
!  at half levels otherwise.

IF(LVERTFE) THEN
  DO JLEV=1,KFLEV
    !dir$ ivdep
    DO JROF=KD,KF
      ZEVEL(JROF,JLEV)=YDVAB%VBF(JLEV)*ZPSDIVFE(JROF,KFLEV+1)-ZPSDIV(JROF,JLEV)
    ENDDO
  ENDDO
ELSE
  DO JLEV=1,KFLEV-1
    !dir$ ivdep
    DO JROF=KD,KF
      ZEVEL(JROF,JLEV)=YDVAB%VBH(JLEV)*ZPSDIV(JROF,KFLEV)-ZPSDIV(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

ZEVEL(KD:KF,0)=0.0_JPRB
IF(.NOT.LVERTFE) ZEVEL(KD:KF,KFLEV)=0.0_JPRB

IF (LDRUBC) THEN
  IF (LVERTFE) THEN
    CALL ABOR1(' GPCTY: VFE not coded for LRUBC.')
  ELSE
    DO JLEV=0,KFLEV
      ZEVEL(KD:KF,JLEV)=ZEVEL(KD:KF,JLEV)+(1.0_JPRB-YDVAB%VBH(JLEV))*PEVT(KD:KF)
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*    4. Computes "(omega/prehyd)" at full levels.
!     --------------------------------------------

IF(LVERTFE) THEN
  DO JLEV=1,KFLEV
    !dir$ ivdep
    DO JROF=KD,KF
      ZVVEL(JROF,JLEV)=(YDVAB%VBF(JLEV)*ZVP(JROF,JLEV)-ZPSDIV(JROF,JLEV))*PRPREF(JROF,JLEV)
    ENDDO
  ENDDO
ELSE
  DO JLEV=1,KFLEV
    !dir$ ivdep
    DO JROF=KD,KF
      ZVVEL(JROF,JLEV)=ZRTGR(JROF,JLEV)*ZVP(JROF,JLEV)&
      & -ZRDELP(JROF,JLEV)*ZALPH(JROF,JLEV)*YDVAB%VDELB(JLEV)*ZVP(JROF,JLEV)&
      & -ZALPH(JROF,JLEV)*PD(JROF,JLEV)& 
      & -ZRDELP(JROF,JLEV)*ZPSDVBC(JROF,JLEV-1)*ZLNPR(JROF,JLEV)  
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPCTY',1,ZHOOK_HANDLE)
END SUBROUTINE GPCTY

