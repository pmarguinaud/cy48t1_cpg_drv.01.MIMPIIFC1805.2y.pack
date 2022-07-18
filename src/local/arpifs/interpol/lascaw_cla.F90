SUBROUTINE LASCAW_CLA(YDSL,KFLEV,KPROM,KST,KPROF,LDT,PSLHDKMIN,KILA,PDLAT,PKAPPA,&
 & PSLD1,PSLD2,PSLD3,PSLDW,PCLA,PCLASLD)

!     ------------------------------------------------------------------

!**** *LASCAW_CLA  -  Weights for semi-LAgrangian interpolator:
!                     Computes PCLA and PCLASLD for one layer
!                     (high-order zonal weights)
!      Spherical geometry only.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LASCAW_CLA( ... )

!        Explicit arguments :
!        --------------------

!        INPUT:
!          YDSL     - SL_STRUCT definition
!          KFLEV    - Vertical dimension
!          KPROM    - horizontal dimension.
!          KST      - first element of arrays where computations are performed.
!          KPROF    - depth of work.
!          LDT      - key for SLHD and horizontal turbulence.
!          KILA     - cf. ILA in LASCAW.
!          PSLHDKMIN - either HOISLH or SLHDKMIN
!          PDLAT    - distance for horizontal linear interpolations in latitude
!          PKAPPA   - kappa function ("coefficient of SLHD").
!          PSLD1   - auxiliary quantity for SLHD interpolation in latitude
!          PSLD2   - auxiliary quantity for SLHD interpolation in latitude
!          PSLD3   - auxiliary quantity for SLHD interpolation in latitude
!          PSLDW   - weights for SLHD Laplacian smoother in latitude

!        OUTPUT:
!          PCLA    - weights for horizontal cubic interpolations in latitude.
!          PCLASLD - cf. PCLA, SLHD case.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------

!        No external.
!        Called by LASCAW.

!     Reference.
!     ----------

!     Author.
!     -------
!        K. YESSAD, after former LASCAW code (JAN 2009).
!        METEO-FRANCE, CNRM/GMAP.

!     Modifications.
!     --------------
!      F. Vana  22-Feb-2011: horizontal turbulence and phys tendencies diff
!      G. Mozdzynski (May 2012): further cleaning
!      S. Malardel (Nov 2013): COMAD weights for SL interpolations
!      F. Vana 13-feb-2014 SLHD weights for heat variables
!      K. Yessad (March 2017): simplify level numbering.
!      F. Vana    21-Nov-2017: Option LHOISLT
!      H Petithomme (Dec 2020): use PCLA already set, option 3DTURB moved out
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIA,JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK 

! arp/ifs dependencies to be solved later.
USE YOMDYNA  , ONLY : YRDYNA
USE EINT_MOD , ONLY : SL_STRUCT,JPDUP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(SL_STRUCT),    INTENT(IN)  :: YDSL
INTEGER(KIND=JPIM), INTENT(IN)  :: KFLEV
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROM
INTEGER(KIND=JPIM), INTENT(IN)  :: KST
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROF
LOGICAL           , INTENT(IN)  :: LDT(4)
INTEGER(KIND=JPIM), INTENT(IN)  :: KILA(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLHDKMIN
REAL(KIND=JPRB)   , INTENT(IN)  :: PDLAT(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKAPPA(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLD1(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLD2(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLD3(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLDW(JPDUP,3,3,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PCLA(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLASLD(KPROM,KFLEV,3)

! note: JPDUP=1, except for NEC machines
#ifndef NECSX
INTEGER(KIND=JPIA),PARAMETER :: I64=1
INTEGER(KIND=JPIM),PARAMETER :: IJ=1
#endif
INTEGER(KIND=JPIA) :: ILA64
INTEGER(KIND=JPIM) :: ILA,JROF,JLEV
LOGICAL :: LLSLHD,LLSLHDQUAD,LLSLHD_OLD
REAL(KIND=JPRB) :: ZWD1,ZWD2,ZWD3,ZWDS1,ZWDS2,ZWDS3,ZWL1,ZWL2,ZWL3
REAL(KIND=JPRB) :: ZSLHDKMIN,ZKAP
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('LASCAW_CLA',0,ZHOOK_HANDLE)

LLSLHD=LDT(1)
LLSLHDQUAD=LDT(2)
LLSLHD_OLD=LDT(3)

! * Calculation of PCLA and PCLASLD:
IF (LLSLHDQUAD.OR.LLSLHD) THEN
  DO JLEV=1,KFLEV
    DO JROF=KST,KPROF
#ifdef NECSX
      IJ = MOD(JROF+1-KST,JPDUP)+1
      I64 = IJ
#endif

      ! warning: LLSLHD and LLSLHDQUAD can both be true, so no test simplification
      ! or reordering. On the contrary, LSLHDQUAD and LSLHD_OLD can not both be true.
      ! These simpler if/else tests help optimization (fewer predicates)
      IF (LLSLHDQUAD) THEN
        ILA=KILA(JROF,JLEV)+1
        ZWL2=PDLAT(JROF,JLEV)
        ZWL1=1.0_JPRB-ZWL2
        ZWL3=PSLD3(IJ,ILA)*ZWL1*ZWL2
        ZWL1=ZWL1+PSLD1(IJ,ILA)*ZWL3
        ZWL2=ZWL2+PSLD2(IJ,ILA)*ZWL3
      ENDIF

      IF (LLSLHD) THEN
        IF (LLSLHD_OLD) THEN
          ZWL2=PDLAT(JROF,JLEV)
          ZWL1=1.0_JPRB-ZWL2
          ZWL3=0.0_JPRB
        ENDIF

        ! optim: 64-bit indexing in PSLDW
        ILA64=KILA(JROF,JLEV)+1
#ifdef __INTEL_COMPILER
        CALL MM_PREFETCH(PSLDW(I64,1,1,ILA64),2)
#endif
        ! note: mind the order of settings for ZWDi/PCLA
        ZWD1=PCLA(JROF,JLEV,1)+YRDYNA%SLHDKMAX*(ZWL1-PCLA(JROF,JLEV,1))
        ZWD2=PCLA(JROF,JLEV,2)+YRDYNA%SLHDKMAX*(ZWL2-PCLA(JROF,JLEV,2))
        ZWD3=PCLA(JROF,JLEV,3)+YRDYNA%SLHDKMAX*(ZWL3-PCLA(JROF,JLEV,3))

        ZWDS1=PSLDW(I64,1,1,ILA64)*ZWD1+PSLDW(I64,1,2,ILA64)*ZWD2+PSLDW(I64,1,3,ILA64)*ZWD3
        ZWDS2=PSLDW(I64,2,1,ILA64)*ZWD1+PSLDW(I64,2,2,ILA64)*ZWD2+PSLDW(I64,2,3,ILA64)*ZWD3
        ZWDS3=PSLDW(I64,3,1,ILA64)*ZWD1+PSLDW(I64,3,2,ILA64)*ZWD2+PSLDW(I64,3,3,ILA64)*ZWD3

        IF (PKAPPA(JROF,JLEV) < 0._JPRB) THEN
          ZKAP=-PKAPPA(JROF,JLEV)
          ZSLHDKMIN=YRDYNA%SLHDKREF
        ELSE
          ZKAP=PKAPPA(JROF,JLEV)
          ZSLHDKMIN=PSLHDKMIN
        ENDIF

        PCLA(JROF,JLEV,1)=PCLA(JROF,JLEV,1)+ZSLHDKMIN*(ZWL1-PCLA(JROF,JLEV,1))
        PCLA(JROF,JLEV,2)=PCLA(JROF,JLEV,2)+ZSLHDKMIN*(ZWL2-PCLA(JROF,JLEV,2))
        PCLA(JROF,JLEV,3)=PCLA(JROF,JLEV,3)+ZSLHDKMIN*(ZWL3-PCLA(JROF,JLEV,3))

        PCLASLD(JROF,JLEV,1)=PCLA(JROF,JLEV,1)+ZKAP*(ZWDS1-PCLA(JROF,JLEV,1))
        PCLASLD(JROF,JLEV,2)=PCLA(JROF,JLEV,2)+ZKAP*(ZWDS2-PCLA(JROF,JLEV,2))
        PCLASLD(JROF,JLEV,3)=PCLA(JROF,JLEV,3)+ZKAP*(ZWDS3-PCLA(JROF,JLEV,3))
      ELSE
        ! note: mind the order of settings for PCLA/PCLASLD
        PCLA(JROF,JLEV,1)=PCLA(JROF,JLEV,1)+PSLHDKMIN*(ZWL1-PCLA(JROF,JLEV,1))
        PCLA(JROF,JLEV,2)=PCLA(JROF,JLEV,2)+PSLHDKMIN*(ZWL2-PCLA(JROF,JLEV,2))
        PCLA(JROF,JLEV,3)=PCLA(JROF,JLEV,3)+PSLHDKMIN*(ZWL3-PCLA(JROF,JLEV,3))
        PCLASLD(JROF,JLEV,1)=PCLA(JROF,JLEV,1)
        PCLASLD(JROF,JLEV,2)=PCLA(JROF,JLEV,2)
        PCLASLD(JROF,JLEV,3)=PCLA(JROF,JLEV,3)
      ENDIF
    ENDDO
  ENDDO
ELSE
  ! optim: separate case since no reference to KILA needed
  DO JLEV=1,KFLEV
    DO JROF=KST,KPROF
      PCLASLD(JROF,JLEV,1)=PCLA(JROF,JLEV,1)
      PCLASLD(JROF,JLEV,2)=PCLA(JROF,JLEV,2)
      PCLASLD(JROF,JLEV,3)=PCLA(JROF,JLEV,3)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_CLA',1,ZHOOK_HANDLE)
END SUBROUTINE LASCAW_CLA
