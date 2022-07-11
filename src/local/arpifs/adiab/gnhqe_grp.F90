!OCL  NOEVAL
SUBROUTINE GNHQE_GRP(&
 ! --- INPUT -----------------------------------------------------------------
 & LDVERTFE, YDGEOMETRY,KST,KEND,&
 & PRT,PRTL,PRTM,PREL,PREM,PRDELP,PRTGR,PALPH,PREHYDH,&
 & PQCHAF,PQCHAL,PQCHAM,PEQCHAF,PEQCHAH,PKAPF,PKAPH,&
 & PHIHL,PHIHM,PHIFL,PHIFM,&
 ! --- OUTPUT ----------------------------------------------------------------
 & PSGRTL,PSGRTM, &
 ! --- INPUT -----------------------------------------------------------------
 & PALPHPLL_DER, PALPHPLM_DER)

!**** *GNHQE_GRP* - Computation of the pressure gradient force term used in the
!                   RHS of the horizontal wind equation in the NHQE model.
!                   This term is computed at full levels.

!     Purpose.
!     --------
!     Computation of the pressure gradient force term, in the NHQE model.
!     Expression is: GRP_HYD + GRP_ADD
!     where:
!      * GRP_HYD looks like the hydrostatic part of the pressure gradient force term
!        but some NH effects are hidden in the modified temperature Tt.
!        GRP_HYD = grad(Phi) + R Tt grad(prehyd)/prehyd
!        This term is treated as closely as in the hydrostatic model.
!      * GRP_ADD is an additional contribution containing some NH effects.
!        GRP_ADD = R Tt [exp((R/Cp) Qcha)-1] grad(prehyd)/prehyd
!                  + R Tt exp((R/Cp) Qcha) grad(Qcha)
!                  + (Cp/R) d (prehyd [exp((R/Cp) Qcha)-1])/d prehyd grad(Phi)
!                  - (Cp/R) (1 - R/Cp) [exp((R/Cp) Qcha)-1] grad(Phi)
!                  - (Cp/R) Qcha exp((R/Cp) Qcha) [prehyd d (R/Cp)/d prehyd] grad(Phi)
!        Term containing [prehyd d (R/Cp)/d prehyd] is neglected and omitted for the time being.
!     This routine uses a modified temperature Tt = T exp(-(R/Cp) Qcha).

!**   Interface.
!     ----------
!        *CALL* *GNHQE_GRP(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!           YDGEOMETRY  : structure containing geometry
!           KST         : start of work
!           KEND        : working length
!           PRT         : "R Tt" at full levels
!           PRTL        : zonal component of "grad (R Tt)" at full levels
!           PRTM        : meridian component of "grad (R Tt)" at full levels
!           PREL        : zonal component of "grad (prehyds)"
!           PREM        : meridian component of "grad (prehyds)"
!           PXYB        : contains pressure depth, "delta", "alpha".
!           PXYBDER     : contains grad(delta), grad(alpha), grad(alpha + log prehyd)
!           PREHYDH     : hydrostatic pressure at half levels.
!           PQCHAF      : Qcha=(log(pre/prehyd)) at full levels
!           PQCHAL      : zonal comp of "grad (log(pre/prehyd))" at full levels
!           PQCHAM      : merid comp of "grad (log(pre/prehyd))" at full levels
!           PEQCHAF     : exp((R/Cp) Qcha) at full levels
!           PEQCHAH     : exp((R/Cp) Qcha) at half levels
!           PKAPF       : (R/Cp) at full levels
!           PKAPH       : (R/Cp) at half levels
!           PHIHL       : zonal component of "grad[Phi]" at half levels (LVERTFE=F only)
!           PHIHM       : merid component of "grad[Phi]" at half levels (LVERTFE=F only)
!           PHIFL       : zonal component of "grad[Phi]" at full levels
!           PHIFM       : merid component of "grad[Phi]" at full levels

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

!     Author.
!     -------
!      K. Yessad (METEO-FRANCE/CNRM/GMAP), after GPGRP.
!      Original : March 2017

!     Modifications.
!     --------------
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
REAL(KIND=JPRB)    ,INTENT(IN)   :: PREHYDH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PQCHAF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PQCHAL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PQCHAM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PEQCHAF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PEQCHAH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PKAPF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PKAPH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PHIHL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PHIHM(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PHIFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PHIFM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(OUT)  :: PSGRTL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(OUT)  :: PSGRTM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PALPHPLL_DER(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PALPHPLM_DER(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSGRTL_HYD(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSGRTM_HYD(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSGRTL_ADD(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSGRTM_ADD(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZUSKAPF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZVDERI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!!! #include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHQE_GRP',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,NPROMA=>YDGEOMETRY%YRDIM%NPROMA)

!     ------------------------------------------------------------------

!*       1.    PSEUDO-HYDROSTATIC PART GRP_HYD.
!              --------------------------------

IF (.NOT.LDVERTFE) THEN

  ! * First get "grad[Phi]" at half levels.
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZSGRTL_HYD(JROF,JLEV)=PHIHL(JROF,JLEV)
      ZSGRTM_HYD(JROF,JLEV)=PHIHM(JROF,JLEV)
    ENDDO
  ENDDO

  ! * Calculation of [ grad(Phi) + R Tt grad(log(prehyd)) ] at full levels.
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZSGRTL_HYD(JROF,JLEV)=ZSGRTL_HYD(JROF,JLEV)&
       & +PALPH(JROF,JLEV)*PRTL(JROF,JLEV)&
       & +PALPHPLL_DER(JROF,JLEV)*PRT(JROF,JLEV)
      ZSGRTM_HYD(JROF,JLEV)=ZSGRTM_HYD(JROF,JLEV)&
       & +PALPH(JROF,JLEV)*PRTM(JROF,JLEV)&
       & +PALPHPLM_DER(JROF,JLEV)*PRT(JROF,JLEV)
    ENDDO
  ENDDO

ELSEIF (LDVERTFE) THEN

  ! * Use (PHIFL,PHIFM) because (PHIHL,PHIHM) is not available in this case.
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZSGRTL_HYD(JROF,JLEV)=PHIFL(JROF,JLEV)+PRT(JROF,JLEV)*PRTGR(JROF,JLEV)*PREL(JROF)
      ZSGRTM_HYD(JROF,JLEV)=PHIFM(JROF,JLEV)+PRT(JROF,JLEV)*PRTGR(JROF,JLEV)*PREM(JROF)
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

!*       2.    ADDITIONAL CONTRIBUTION GRP_ADD.
!              --------------------------------

! * Pre-compute (Cp/R) at full levels:
DO JLEV=1,NFLEVG
  DO JROF=KST,KEND
    ZUSKAPF(JROF,JLEV)=1.0_JPRB/PKAPF(JROF,JLEV)
  ENDDO
ENDDO

! * Term "R Tt [exp((R/Cp) Qcha)-1] grad(prehyd)/prehyd" at full levels:
DO JLEV=1,NFLEVG
  DO JROF=KST,KEND
    ZSGRTL_ADD(JROF,JLEV)=PRT(JROF,JLEV)*(PEQCHAF(JROF,JLEV)-1.0_JPRB)*PRTGR(JROF,JLEV)*PREL(JROF)
    ZSGRTM_ADD(JROF,JLEV)=PRT(JROF,JLEV)*(PEQCHAF(JROF,JLEV)-1.0_JPRB)*PRTGR(JROF,JLEV)*PREM(JROF)
  ENDDO
ENDDO

! * Term "R Tt exp((R/Cp) Qcha) grad(Qcha)" at full levels:
DO JLEV=1,NFLEVG
  DO JROF=KST,KEND
    ZSGRTL_ADD(JROF,JLEV)=ZSGRTL_ADD(JROF,JLEV)+PRT(JROF,JLEV)*PEQCHAF(JROF,JLEV)*PQCHAL(JROF,JLEV)
    ZSGRTM_ADD(JROF,JLEV)=ZSGRTM_ADD(JROF,JLEV)+PRT(JROF,JLEV)*PEQCHAF(JROF,JLEV)*PQCHAM(JROF,JLEV)
  ENDDO
ENDDO

! * Term "(Cp/R) d (prehyd [exp((R/Cp) Qcha)-1])/d prehyd grad(Phi)" at full levels:
!   This piece of code is valid for VFD and VFE discretisation of grad(Phi), but
!   !!!!!! VFE discretisation of vertical derivative "ZVDERI = d (prehyd [exp((R/Cp) Qcha)-1])/d prehyd" is not yet coded !!!!!!
DO JLEV=1,NFLEVG
  DO JROF=KST,KEND
    ZVDERI(JROF,JLEV)=PRDELP(JROF,JLEV)*&
     & (PREHYDH(JROF,JLEV)*(PEQCHAH(JROF,JLEV)-1.0_JPRB) - PREHYDH(JROF,JLEV-1)*(PEQCHAH(JROF,JLEV-1)-1.0_JPRB))
    ZSGRTL_ADD(JROF,JLEV)=ZSGRTL_ADD(JROF,JLEV)+ZUSKAPF(JROF,JLEV)*ZVDERI(JROF,JLEV)*PHIFL(JROF,JLEV)
    ZSGRTM_ADD(JROF,JLEV)=ZSGRTM_ADD(JROF,JLEV)+ZUSKAPF(JROF,JLEV)*ZVDERI(JROF,JLEV)*PHIFM(JROF,JLEV)
  ENDDO
ENDDO

! * Term "- (Cp/R) (1 - R/Cp) [exp((R/Cp) Qcha)-1] grad(Phi)" at full levels:
!   This piece of code is valid for VFD and VFE discretisation of grad(Phi).
DO JLEV=1,NFLEVG
  DO JROF=KST,KEND
    ZSGRTL_ADD(JROF,JLEV)=ZSGRTL_ADD(JROF,JLEV)+(1.0_JPRB-ZUSKAPF(JROF,JLEV))*(PEQCHAF(JROF,JLEV)-1.0_JPRB)*PHIFL(JROF,JLEV)
    ZSGRTM_ADD(JROF,JLEV)=ZSGRTM_ADD(JROF,JLEV)+(1.0_JPRB-ZUSKAPF(JROF,JLEV))*(PEQCHAF(JROF,JLEV)-1.0_JPRB)*PHIFM(JROF,JLEV)
  ENDDO
ENDDO

! * Term "- (Cp/R) Qcha exp((R/Cp) Qcha) [prehyd d (R/Cp)/d prehyd] grad(Phi)" at full levels:
!   This piece of code is valid for VFD and VFE discretisation of grad(Phi), but
!   !!!!!! VFE discretisation of vertical derivative "ZVDERI = [prehyd d (R/Cp)/d prehyd]" is not yet coded !!!!!!
!   !!!!!! Currently neglected and omitted !!!!!!
!!! DO JLEV=1,NFLEVG
!!!   DO JROF=KST,KEND
!!!     ZVDERI(JROF,JLEV)=(PKAPH(JROF,JLEV)-PKAPH(JROF,JLEV-1))/PLNPR(JROF,JLEV)
!!!     ZSGRTL_ADD(JROF,JLEV)=ZSGRTL_ADD(JROF,JLEV)&
!!!      & -ZUSKAPF(JROF,JLEV)*PQCHAF(JROF,JLEV)*PEQCHAF(JROF,JLEV)*ZVDERI(JROF,JLEV)*PHIFL(JROF,JLEV)
!!!     ZSGRTM_ADD(JROF,JLEV)=ZSGRTM_ADD(JROF,JLEV)&
!!!      & -ZUSKAPF(JROF,JLEV)*PQCHAF(JROF,JLEV)*PEQCHAF(JROF,JLEV)*ZVDERI(JROF,JLEV)*PHIFM(JROF,JLEV)
!!!   ENDDO
!!! ENDDO

!     ------------------------------------------------------------------

!*       3.    FINAL CALCULATION.
!              ------------------

DO JLEV=1,NFLEVG
  DO JROF=KST,KEND
    PSGRTL(JROF,JLEV)=ZSGRTL_HYD(JROF,JLEV)+ZSGRTL_ADD(JROF,JLEV)
    PSGRTM(JROF,JLEV)=ZSGRTM_HYD(JROF,JLEV)+ZSGRTM_ADD(JROF,JLEV)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHQE_GRP',1,ZHOOK_HANDLE)
END SUBROUTINE GNHQE_GRP
