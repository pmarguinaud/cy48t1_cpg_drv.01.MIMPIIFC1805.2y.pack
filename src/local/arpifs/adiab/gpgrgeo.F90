!OCL  NOEVAL
SUBROUTINE GPGRGEO(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY,KPROMA,KD,KF,KFLEV,&
 & PRT,PRTL,PRTM,&
 & PLNPR,PALPH,PXYBDER,&
 & POROGL,POROGM,&
 ! --- OUTPUT ----------------------------------------------------------------
 & PHIFL,PHIFM,PHIHL,PHIHM,&
 ! --- OPTIONAL INPUT --------------------------------------------------------
 & LDNHEE,PRNHPPI,PQCHAL,PQCHAM,&
 ! --- OPTIONAL OUTPUT -------------------------------------------------------
 & PNH1L,PNH1M,&
 ! --- OPTIONAL INPUT --------------------------------------------------------
 & PLNPRL_DER, PLNPRM_DER, PALPHL_DER, PALPHM_DER)  

!**** *GPGRGEO* - Computes half and full level gradient of geopotential height "gz".

!     Purpose.
!     --------

!      Expression of this term is:

!       grad (gz) = grad Phi_s
!       + grad {int[prehyd'=prehyds to prehyd] (-RT/pre) d prehyd'}

!      where: 
!       - "Phi_s = g z[surf]" is the surface orography.
!       - "prehyd" is the hydrostatic pressure.
!       - "pre" is the total pressure including non-hydrostatic effects.
!       - "prehyds" is the surface hydrostatic pressure.
!       - "R" is the air constant (including water contribution).
!       - "T" is the temperature.

!      For the NHQE and the HYD model, pre is replaced by prehyd, some terms disappear.
!      For the NHQE model, T is a modified temperature

!      Discretisation of the gradient of geopotential height yields:

!      * "grad (gz)" at half level "lbar":

!        (grad (gz))[lbar] = (grad Phi_s)
!        + sum[k=L to l+1] (prehyd/pre)[k] R[k] T[k] (grad (delta))[k]
!        + sum[k=L to l+1] (prehyd/pre)[k] (grad (RT))[k] delta[k]
!        + sum[k=L to l+1] (prehyd/pre)[k] R[k] T[k] delta[k]
!          { (grad(prehyd)/prehyd)[k] - (grad(pre)/pre)[k] }

!      * "grad (gz)" at full level "l":

!        (grad (gz))[l] = (grad (gz))[lbar]
!        + (prehyd/pre)[l] R[l] T[l] (grad (alpha))[l]
!        + (prehyd/pre)[l] (grad (RT))[l] alpha[l]
!        + (prehyd/pre)[l] R[l] T[l] alpha[l]
!          { (grad(prehyd)/prehyd)[l] - (grad(pre)/pre)[l] }

!**   Interface.
!     ----------
!        *CALL* *GPGRGEO(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!           YDGEOMETRY   : structure containing all geometry
!           KPROMA       : horizontal dimension
!           KD           : start of work
!           KF           : working length
!           KFLEV        : number of levels
!           PRT          : (RT) at full levels
!           PRTL         : zonal component of "grad (RT)" at full levels
!           PRTM         : meridian component of "grad (RT)" at full levels
!           PLNPR        : "delta" at full levels
!           PALPH        : "alpha" at full levels
!           PXYBDER      : contains grad(delta), grad(alpha), grad(alpha + log prehyd)
!           POROGL       : zonal component of "grad(surf orography)"
!           POROGM       : meridian component of "grad(surf orography)"

!         * OUTPUT:
!           PHIFL        : zonal component of "grad (gz)" at full levels
!           PHIFM        : merid component of "grad (gz)" at full levels
!           PHIHL        : zonal component of "grad (gz)" at half levels
!           PHIHM        : merid component of "grad (gz)" at half levels

!         * INPUT OPTIONAL:
!           LDNHEE       : .T.: fully elastic non hydrostatic (NHEE) model.
!                          .F.: hydrostatic or NHQE model.
!           PRNHPPI      : "prehyd/pre" at full levels.
!           PQCHAL,PQCHAM: zonal and meridian components at full levels of
!                          "grad(log(pre/prehyd))=(grad pre)/pre - (grad(prehyd))/prehyd"

!         * OUTPUT OPTIONAL:
!           PNH1L        : zonal comp of RT grad(log(prehyd/pre)) at full levels
!           PNH1M        : merid comp of RT grad(log(prehyd/pre)) at full levels

!        Implicit arguments :   None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    None.
!     ----------

!     Reference.
!     ----------
!        See  documentation Arpege ALGORITHME CHAPITRE 6 paragraphe 6

!     Author.
!     -------
!      K. YESSAD
!      Original : 2000-08-11

!     Modifications.
!     --------------
!      K. Yessad (Dec 2008): remove dummy CDLOCK
!      K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!      K. Yessad (June 2017): Introduce NHQE model.
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!      H. Petithomme (Nov 2020): use of pointers for avoiding array copies
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCVER      , ONLY : LVERTFE
USE INTDYN_MOD   , ONLY : YYTXYBDER

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY),    INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KD 
INTEGER(KIND=JPIM),INTENT(IN)    :: KF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRT(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRTL(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRTM(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLNPR(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPH(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN),TARGET :: PXYBDER(KPROMA,KFLEV,YYTXYBDER%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHIFL(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHIFM(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHIHL(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHIHM(KPROMA,0:KFLEV) 
LOGICAL,OPTIONAL,INTENT(IN)          :: LDNHEE
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)  :: PRNHPPI(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)  :: PQCHAL(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)  :: PQCHAM(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PNH1L(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PNH1M(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN),TARGET :: PLNPRL_DER (KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN),TARGET :: PLNPRM_DER (KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN),TARGET :: PALPHL_DER (KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN),TARGET :: PALPHM_DER (KPROMA,KFLEV)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF
LOGICAL         :: LLNHEE
REAL(KIND=JPRB) :: ZNH1L(KPROMA,KFLEV), ZNH1M(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZGPHL(KPROMA,KFLEV+1),ZGPHM(KPROMA,KFLEV+1)
REAL(KIND=JPRB) :: ZINL(KPROMA,0:KFLEV+1),ZINM(KPROMA,0:KFLEV+1)
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZLNPRL_DER (:,:), ZLNPRM_DER (:,:), ZALPHL_DER (:,:), ZALPHM_DER (:,:)
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPGRGEO',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

IF (PRESENT(LDNHEE)) THEN
  LLNHEE=LDNHEE
ELSE
  LLNHEE=.FALSE.
ENDIF

ZLNPRL_DER => NULL ()
ZLNPRM_DER => NULL ()
ZALPHL_DER => NULL ()
ZALPHM_DER => NULL ()

IF (PRESENT (PXYBDER)) THEN
  ZLNPRL_DER => PXYBDER(:,:,YYTXYBDER%M_LNPRL)
  ZLNPRM_DER => PXYBDER(:,:,YYTXYBDER%M_LNPRM)
  ZALPHL_DER => PXYBDER(:,:,YYTXYBDER%M_ALPHL)
  ZALPHM_DER => PXYBDER(:,:,YYTXYBDER%M_ALPHM)
ENDIF

IF (PRESENT (PLNPRL_DER)) ZLNPRL_DER => PLNPRL_DER
IF (PRESENT (PLNPRM_DER)) ZLNPRM_DER => PLNPRM_DER
IF (PRESENT (PALPHL_DER)) ZALPHL_DER => PALPHL_DER
IF (PRESENT (PALPHM_DER)) ZALPHM_DER => PALPHM_DER

!     ------------------------------------------------------------------

!*    1. Computation of term:
!        R[l] T[l] {(grad(prehyd)/prehyd)[l] - (grad(pre)/pre)[l]}
!        according to "llnhee".
!*    2.1 Calculation of "grad (gz)" at half levels.
!         ("delta" and "grad delta" terms contributions.

DO JROF=KD,KF
  PHIHL(JROF,KFLEV)=POROGL(JROF)
  PHIHM(JROF,KFLEV)=POROGM(JROF)
ENDDO

IF (LLNHEE) THEN
  IF (.NOT.(PRESENT(PRNHPPI).AND.PRESENT(PQCHAL).AND.PRESENT(PQCHAM)))&
    CALL ABOR1(' GPGRGEO: missing input PRNHPPI, PQCHAL, PQCHAM !!!')

  DO JLEV=1,KFLEV
    ZNH1L(KD:KF,JLEV)=-PRT(KD:KF,JLEV)*PQCHAL(KD:KF,JLEV)
    ZNH1M(KD:KF,JLEV)=-PRT(KD:KF,JLEV)*PQCHAM(KD:KF,JLEV)
  ENDDO

  DO JLEV=KFLEV,1,-1
    DO JROF=KD,KF
      PHIHL(JROF,JLEV-1)= PHIHL(JROF,JLEV)&
       & +PLNPR(JROF,JLEV)*PRNHPPI(JROF,JLEV)*(PRTL(JROF,JLEV)+ZNH1L(JROF,JLEV))&
       & +ZLNPRL_DER(JROF,JLEV)*PRNHPPI(JROF,JLEV)*PRT(JROF,JLEV)
      PHIHM(JROF,JLEV-1)= PHIHM(JROF,JLEV)&
       & +PLNPR(JROF,JLEV)*PRNHPPI(JROF,JLEV)*(PRTM(JROF,JLEV)+ZNH1M(JROF,JLEV))&
       & +ZLNPRM_DER(JROF,JLEV)*PRNHPPI(JROF,JLEV)*PRT(JROF,JLEV)
    ENDDO
  ENDDO

  ! note: merge tests since both or none present (only present in cpg5_cp)
  IF (PRESENT(PNH1L).AND.PRESENT(PNH1M)) THEN
    PNH1L(KD:KF,1:KFLEV)=ZNH1l(KD:KF,1:KFLEV)
    PNH1M(KD:KF,1:KFLEV)=ZNH1M(KD:KF,1:KFLEV)
  END IF
ELSE
  ! if not present, these arrays are unused
  IF (PRESENT(PNH1L).AND.PRESENT(PNH1M)) THEN
    PNH1L(KD:KF,1:KFLEV)=0.0_JPRB 
    PNH1M(KD:KF,1:KFLEV)=0.0_JPRB
  ENDIF

  DO JLEV=KFLEV,1,-1
    DO JROF=KD,KF
      PHIHL(JROF,JLEV-1)=PHIHL(JROF,JLEV)+PLNPR(JROF,JLEV)*PRTL(JROF,JLEV)&
       & +ZLNPRL_DER(JROF,JLEV)*PRT(JROF,JLEV)
      PHIHM(JROF,JLEV-1)=PHIHM(JROF,JLEV)+PLNPR(JROF,JLEV)*PRTM(JROF,JLEV)&
       & +ZLNPRM_DER(JROF,JLEV)*PRT(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

!*    2.2 Calculation of "grad (gz)" at half levels.
!         "alpha" and "grad alpha" terms contributions.

IF(LVERTFE) THEN
  IF (LLNHEE) THEN
    DO JLEV=1,KFLEV
      DO JROF=KD,KF
        ZINL(JROF,JLEV)=(-PLNPR(JROF,JLEV)*PRNHPPI(JROF,JLEV)*&
         & (PRTL(JROF,JLEV)+ZNH1L(JROF,JLEV))&
         & -ZLNPRL_DER(JROF,JLEV)*PRNHPPI(JROF,JLEV)*PRT(JROF,JLEV))&
         & *YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)  
        ZINM(JROF,JLEV)=(-PLNPR(JROF,JLEV)*PRNHPPI(JROF,JLEV)*&
         & (PRTM(JROF,JLEV)+ZNH1M(JROF,jLEV))&
         & -ZLNPRM_DER(JROF,JLEV)*PRNHPPI(JROF,JLEV)*PRT(JROF,JLEV))&
         & *YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)  
      ENDDO
    ENDDO
  ELSE
    DO JLEV=1,KFLEV
      DO JROF=KD,KF
        ZINL(JROF,JLEV)=(-PLNPR(JROF,JLEV)*PRTL(JROF,JLEV)&
         & -ZLNPRL_DER(JROF,JLEV)*PRT(JROF,JLEV))&
         & *YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)
        ZINM(JROF,JLEV)=(-PLNPR(JROF,JLEV)*PRTM(JROF,JLEV)&
         & -ZLNPRM_DER(JROF,JLEV)*PRT(JROF,JLEV))&
         & *YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)
      ENDDO
    ENDDO
  ENDIF

  ZINL(KD:KF,0)=0.0_JPRB
  ZINL(KD:KF,KFLEV+1)=0.0_JPRB
  ZINM(KD:KF,0)=0.0_JPRB
  ZINM(KD:KF,KFLEV+1)=0.0_JPRB
  CALL VERDISINT(YDGEOMETRY%YRVFE,'IBOT','11',KPROMA,KD,KF,KFLEV,ZINL,ZGPHL)
  CALL VERDISINT(YDGEOMETRY%YRVFE,'IBOT','11',KPROMA,KD,KF,KFLEV,ZINM,ZGPHM)

  DO JLEV=1,KFLEV
    DO JROF=KD,KF
      PHIFL(JROF,JLEV)=ZGPHL(JROF,JLEV)+POROGL(JROF)
      PHIFM(JROF,JLEV)=ZGPHM(JROF,JLEV)+POROGM(JROF)
    ENDDO  
  ENDDO  
ELSE
  IF (LLNHEE) THEN
    DO JLEV=1,KFLEV
      DO JROF=KD,KF
        PHIFL(JROF,JLEV)=PHIHL(JROF,JLEV)&
         & +PALPH(JROF,JLEV)*PRNHPPI(JROF,JLEV)*(PRTL(JROF,JLEV)+ZNH1L(JROF,JLEV))&
         & +ZALPHL_DER(JROF,JLEV)*PRNHPPI(JROF,JLEV)*PRT(JROF,JLEV)
        PHIFM(JROF,JLEV)=PHIHM(JROF,JLEV)&
         & +PALPH(JROF,JLEV)*PRNHPPI(JROF,JLEV)*(PRTM(JROF,JLEV)+ZNH1M(JROF,JLEV))&
         & +ZALPHM_DER(JROF,JLEV)*PRNHPPI(JROF,JLEV)*PRT(JROF,JLEV)
      ENDDO
    ENDDO
  ELSE
    DO JLEV=1,KFLEV
      DO JROF=KD,KF
        PHIFL(JROF,JLEV)=PHIHL(JROF,JLEV)+PALPH(JROF,JLEV)*PRTL(JROF,JLEV)&
         & +ZALPHL_DER(JROF,JLEV)*PRT(JROF,JLEV)
        PHIFM(JROF,JLEV)=PHIHM(JROF,JLEV)+PALPH(JROF,JLEV)*PRTM(JROF,JLEV)&
         & +ZALPHM_DER(JROF,JLEV)*PRT(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPGRGEO',1,ZHOOK_HANDLE)
END SUBROUTINE GPGRGEO
