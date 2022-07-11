!OCL  NOEVAL
SUBROUTINE GNHGRGW(&
 ! --- INPUT -----------------------------------------------------------------
 & LDVERTFE, YDGEOMETRY,KPROMA,KFLEV,KD,KF,LDGWFDER,&
 & PRT,PRTL,PRTM,PDVER,PDVERL,PDVERM,&
 & PLNPR,PALPH,&
 & PUS,PVS,PUS_L,PVS_L,PUS_M,PVS_M,&
 & POROGL,POROGM,POROGLM,POROGLL,POROGMM,&
 ! --- OUTPUT ----------------------------------------------------------------
 & PGWFL,PGWFM,PGWHL,PGWHM,&
 ! --- OPTIONAL INPUT --------------------------------------------------------
 & LDNHEE,PRNHPPI,PQCHAL,PQCHAM,&
 & PLNPRL_DER, PLNPRM_DER, PALPHL_DER, PALPHM_DER)

!**** *GNHGRGW* - Computes half and full level gradient of "gw"
!                 ("g" times the vertical velocity "w").

!     Purpose.
!     --------
!     Diagnoses "grad(gw)" from the vertical divergence "dver" and its gradient.
!     For the finite element vertical discretization (lvertfe=T), "grad(gw)" is
!      computed only at full levels.
!     For the finite difference vertical discretization (lvertfe=F), "grad(gw)"
!      is computed only at both half and full levels.
!     Calculation is done by applying the gradient operator to the
!      vertical integration of the formula:
!      dver = - (g pre)/(Rd T prehyd) (d w / d log(prehyd))
!     with the following bottom condition:
!      w_surf = V_surf grad[Phi_s].

!     For the NHQE, pre is replaced by prehyd, some terms disappear.
!     For the NHQE model, T is a modified temperature
!     This routine is useless and should not be called in the hydrostatic model.

!**   Interface.
!     ----------
!        *CALL* *GNHGRGW(...)

!     Explicit arguments :
!     --------------------
!      * INPUT:
!        YDGEOMETRY   : structure containing all geometry.
!        KPROMA       : horizontal dimension.
!        KFLEV        : number of levels.
!        KD           : start of work.
!        KF           : working length.
!        LDGWFDER     : compute (PGWFL,PGWFM) in all cases if T
!        PRT          : (RT) at full levels, with the version of R used to define vertical divergence "dver".
!                       R may be Rdry or Rmoist according to definition of vertical divergence "dver".
!        PRTL         : zonal component of "grad (RT)" at full levels.
!        PRTM         : meridian component of "grad (RT)" at full levels.
!        PDVER        : vertical divergence "dver" at full levels.
!        PDVERL       : zonal comp. of grad(dver) at full levels.
!        PDVERM       : merid comp. of grad(dver) at full levels.
!        PLNPR        : "delta" at full levels.
!        PALPH        : "alpha" at full levels.
!        PXYBDER      : contains grad(delta), grad(alpha), grad(alpha + log prehyd)
!        PUS          : U-surface wind.
!        PVS          : V-surface wind.
!        PUS_L        : zonal comp. of grad(U-surface wind).
!        PVS_L        : zonal comp. of grad(V-surface wind).
!        PUS_M        : merid comp. of grad(U-surface wind).
!        PVS_M        : merid comp. of grad(V-surface wind).
!        POROGL       : zonal comp. of grad(surf orography).
!        POROGM       : merid comp. of grad(surf orography).
!        POROGLM      : -I
!        POROGLL      :  I- second order derivatives of "surf orography".
!        POROGMM      : -I

!      * OUTPUT:
!        PGWFL        : zonal comp. of "grad(gw)" at full levels.
!        PGWFM        : merid comp. of "grad(gw)" at full levels.
!        PGWHL        : zonal comp. of "grad(gw)" at half levels (lvertfe=F only).
!        PGWHM        : merid comp. of "grad(gw)" at half levels (lvertfe=F only).

!      * INPUT OPTIONAL:
!        LDNHEE       : .T.: fully elastic non hydrostatic (NHEE) model.
!                       .F.: hydrostatic or NHQE model.
!        PRNHPPI      : "prehyd/pre" at full levels.
!        PQCHAL,PQCHAM: zonal and meridian components at full levels of
!                       "grad(log(pre/prehyd))=(grad pre)/pre - (grad(prehyd))/prehyd"

!     Implicit arguments :   None.
!     --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    None.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!       K. YESSAD (after routines GNHPDVD, GNHGRP and GPGRGEO).
!       Original : 07-Dec-2004

!     Modifications.
!     --------------
!       K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!       F. Vana 21-Feb-2011: FL computation of grad(gw) when L3DTURB
!       K. Yessad (June 2017): Introduce NHQE model.
!       J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!       J. Vivoda and P. Smolikova (Sep 2020): VFE pruning.
!       H Petithomme (Dec 2020): use of pointers
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK


!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL            ,INTENT(IN)                     :: LDVERTFE
TYPE(GEOMETRY)     ,INTENT(IN)                     :: YDGEOMETRY
INTEGER(KIND=JPIM) ,INTENT(IN)                     :: KPROMA
INTEGER(KIND=JPIM) ,INTENT(IN)                     :: KFLEV
INTEGER(KIND=JPIM) ,INTENT(IN)                     :: KD
INTEGER(KIND=JPIM) ,INTENT(IN)                     :: KF
LOGICAL            ,INTENT(IN)                     :: LDGWFDER
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PRT(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PRTL(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PRTM(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PDVER(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PDVERL(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PDVERM(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PLNPR(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PALPH(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PUS(KPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PVS(KPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PUS_L(KPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PVS_L(KPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PUS_M(KPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PVS_M(KPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: POROGL(KPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: POROGM(KPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: POROGLM(KPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: POROGLL(KPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: POROGMM(KPROMA)
REAL(KIND=JPRB)    ,INTENT(OUT)                    :: PGWFL(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(OUT)                    :: PGWFM(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(OUT)           ,TARGET  :: PGWHL(KPROMA,0:KFLEV)
REAL(KIND=JPRB)    ,INTENT(OUT)           ,TARGET  :: PGWHM(KPROMA,0:KFLEV)
LOGICAL            ,INTENT(IN)  ,OPTIONAL          :: LDNHEE
REAL(KIND=JPRB)    ,INTENT(IN)  ,OPTIONAL          :: PRNHPPI(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)  ,OPTIONAL          :: PQCHAL(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)  ,OPTIONAL          :: PQCHAM(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PLNPRL_DER(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PLNPRM_DER(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PALPHL_DER(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)                     :: PALPHM_DER(KPROMA,KFLEV)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZGPHL(KPROMA,0:KFLEV),ZGPHM(KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZRTPI,ZRTV
REAL(KIND=JPRB) :: ZINL(KPROMA,0:KFLEV+1),ZINM(KPROMA,0:KFLEV+1)
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZGWSL(:),ZGWSM(:)

LOGICAL :: LLNHEE
INTEGER(KIND=JPIM) :: JL, JROF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHGRGW',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

IF (PRESENT(LDNHEE)) THEN
  LLNHEE=LDNHEE
ELSE
  LLNHEE=.FALSE.
ENDIF

! Remarks:
!  if LLNHEE=T: PRNHPPI, PQCHAL and PQCHAM are used and should be present
!  if LLNHEE=F: PRNHPPI, PQCHAL and PQCHAM should not be used
IF ( LLNHEE.AND. .NOT.(PRESENT(PRNHPPI).AND.PRESENT(PQCHAL).AND.PRESENT(PQCHAM)) ) THEN
  CALL ABOR1(' GNHGRGW: missing input PRNHPPI, PQCHAL, PQCHAM !!!')
ENDIF

!     ------------------------------------------------------------------

!*    1.1 Calculation of "grad (gw)" at the surface.
! optim: use pointers on last level
ZGWSL => PGWHL(:,KFLEV)
ZGWSM => PGWHM(:,KFLEV)

DO JROF=KD,KF
  ZGWSL(JROF)=PUS_L(JROF)*POROGL(JROF)+PUS(JROF)*POROGLL(JROF) &
   & +PVS_L(JROF)*POROGM(JROF)+PVS(JROF)*POROGLM(JROF)
  ZGWSM(JROF)=PUS_M(JROF)*POROGL(JROF)+PUS(JROF)*POROGLM(JROF) &
   & +PVS_M(JROF)*POROGM(JROF)+PVS(JROF)*POROGMM(JROF)
ENDDO

!*    1.2 Calculation of "grad (gw)" for the Z-term.

IF (LDVERTFE) THEN

  ! * Calculation of "grad (gw)" at full levels.

  IF (LLNHEE) THEN
    DO JL=1,KFLEV
      DO JROF=KD,KF
        ZRTPI = PRT(JROF,JL)*PRNHPPI(JROF,JL)
        ZRTV = ZRTPI*PDVER(JROF,JL)
        ZINL(JROF,JL)=YDGEOMETRY%YRVETA%VFE_RDETAH(JL)*(&
         & (-PRTL(JROF,JL)*PRNHPPI(JROF,JL)*PDVER(JROF,JL)&
         & + ZRTPI*PQCHAL(JROF,JL)*PDVER(JROF,JL)&
         & - ZRTPI*PDVERL(JROF,JL) )*PLNPR(JROF,JL)&
         & - ZRTV*PLNPRL_DER(JROF,JL) )
        ZINM(JROF,JL)=YDGEOMETRY%YRVETA%VFE_RDETAH(JL)*(&
         & (-PRTM(JROF,JL)*PRNHPPI(JROF,JL)*PDVER(JROF,JL)&
         & + ZRTPI*PQCHAM(JROF,JL)*PDVER(JROF,JL)&
         & - ZRTPI*PDVERM(JROF,JL) )*PLNPR(JROF,JL)&
         & - ZRTV*PLNPRM_DER(JROF,JL) )
      ENDDO
    ENDDO
  ELSE
    DO JL=1,KFLEV
      DO JROF=KD,KF
        ZRTV = PRT(JROF,JL)*PDVER(JROF,JL)
        ZINL(JROF,JL)=YDGEOMETRY%YRVETA%VFE_RDETAH(JL)*(&
         & (-PRTL(JROF,JL)*PDVER(JROF,JL)-PRT(JROF,JL)*PDVERL(JROF,JL))*PLNPR(JROF,JL)&
         & - ZRTV*PLNPRL_DER(JROF,JL) )
        ZINM(JROF,JL)=YDGEOMETRY%YRVETA%VFE_RDETAH(JL)*(&
         & (-PRTM(JROF,JL)*PDVER(JROF,JL)-PRT(JROF,JL)*PDVERM(JROF,JL))*PLNPR(JROF,JL)&
         & - ZRTV*PLNPRM_DER(JROF,JL) )
      ENDDO
    ENDDO
  ENDIF

  DO JROF=KD,KF
    ZINL(JROF,0)=0.0_JPRB
    ZINL(JROF,KFLEV+1)=0.0_JPRB
    ZINM(JROF,0)=0.0_JPRB
    ZINM(JROF,KFLEV+1)=0.0_JPRB
  ENDDO
  CALL VERDISINT(YDGEOMETRY%YRVFE,'IBOT','11',KPROMA,KD,KF,KFLEV,ZINL,ZGPHL)
  CALL VERDISINT(YDGEOMETRY%YRVFE,'IBOT','11',KPROMA,KD,KF,KFLEV,ZINM,ZGPHM)

  DO JL=1,KFLEV
    DO JROF=KD,KF
      PGWFL(JROF,JL)=ZGPHL(JROF,JL-1)+ZGWSL(JROF)
      PGWFM(JROF,JL)=ZGPHM(JROF,JL-1)+ZGWSM(JROF)
    ENDDO
  ENDDO

ELSE

  ! * Calculation of "grad (gw)" at half levels.
  !   "delta" and "grad delta" terms contributions.

  IF (LLNHEE) THEN
    DO JL=KFLEV,1,-1
      !DIR$ IVDEP
      !CDIR NODEP
      DO JROF=KD,KF
        ZRTPI = PRT(JROF,JL)*PRNHPPI(JROF,JL)
        ZRTV = ZRTPI*PDVER(JROF,JL)
        PGWHL(JROF,JL-1)=PGWHL(JROF,JL) &
         & +(PRTL(JROF,JL)*PRNHPPI(JROF,JL)*PDVER(JROF,JL) &
         & - ZRTPI*PQCHAL(JROF,JL)*PDVER(JROF,JL) &
         & + ZRTPI*PDVERL(JROF,JL) )*PLNPR(JROF,JL) &
         & + ZRTV*PLNPRL_DER(JROF,JL)
        PGWHM(JROF,JL-1)=PGWHM(JROF,JL) &
         & +(PRTM(JROF,JL)*PRNHPPI(JROF,JL)*PDVER(JROF,JL) &
         & - ZRTPI*PQCHAM(JROF,JL)*PDVER(JROF,JL) &
         & + ZRTPI*PDVERM(JROF,JL) )*PLNPR(JROF,JL) &
         & + ZRTV*PLNPRM_DER(JROF,JL)
      ENDDO
    ENDDO
  ELSE
    DO JL=KFLEV,1,-1
      !DIR$ IVDEP
      !CDIR NODEP
      DO JROF=KD,KF
        ZRTV = PRT(JROF,JL)*PDVER(JROF,JL)
        PGWHL(JROF,JL-1)=PGWHL(JROF,JL) &
         & +(PRTL(JROF,JL)*PDVER(JROF,JL)+PRT(JROF,JL)*PDVERL(JROF,JL))*PLNPR(JROF,JL) &
         & + ZRTV*PLNPRL_DER(JROF,JL)
        PGWHM(JROF,JL-1)=PGWHM(JROF,JL) &
         & +(PRTM(JROF,JL)*PDVER(JROF,JL)+PRT(JROF,JL)*PDVERM(JROF,JL))*PLNPR(JROF,JL) &
         & + ZRTV*PLNPRM_DER(JROF,JL)
      ENDDO
    ENDDO
  ENDIF

  ! * Calculation of "grad (gw)" at full levels.
  !   "alpha" and "grad alpha" terms contributions.

  IF (LDGWFDER) THEN
    IF (LLNHEE) THEN
      DO JL=1,KFLEV
        DO JROF=KD,KF
          ZRTPI = PRT(JROF,JL)*PRNHPPI(JROF,JL)
          ZRTV = ZRTPI*PDVER(JROF,JL)
          PGWFL(JROF,JL)=PGWHL(JROF,JL) &
           & +(PRTL(JROF,JL)*PRNHPPI(JROF,JL)*PDVER(JROF,JL) &
           & - ZRTPI*PQCHAL(JROF,JL)*PDVER(JROF,JL) &
           & + ZRTPI*PDVERL(JROF,JL) )*PALPH(JROF,JL) &
           & + ZRTV*PALPHL_DER(JROF,JL)
          PGWFM(JROF,JL)=PGWHM(JROF,JL) &
           & +(PRTM(JROF,JL)*PRNHPPI(JROF,JL)*PDVER(JROF,JL) &
           & - ZRTPI*PQCHAM(JROF,JL)*PDVER(JROF,JL) &
           & + ZRTPI*PDVERM(JROF,JL) )*PALPH(JROF,JL) &
           & + ZRTV*PALPHM_DER(JROF,JL)
        ENDDO
      ENDDO
    ELSE
      DO JL=1,KFLEV
        DO JROF=KD,KF
          ZRTV = PRT(JROF,JL)*PDVER(JROF,JL)
          PGWFL(JROF,JL)=PGWHL(JROF,JL) &
           & +(PRTL(JROF,JL)*PDVER(JROF,JL)+PRT(JROF,JL)*PDVERL(JROF,JL))*PALPH(JROF,JL) &
           & + ZRTV*PALPHL_DER(JROF,JL)
          PGWFM(JROF,JL)=PGWHM(JROF,JL) &
           & +(PRTM(JROF,JL)*PDVER(JROF,JL)+PRT(JROF,JL)*PDVERM(JROF,JL))*PALPH(JROF,JL) &
           & + ZRTV*PALPHM_DER(JROF,JL)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHGRGW',1,ZHOOK_HANDLE)
END SUBROUTINE GNHGRGW
