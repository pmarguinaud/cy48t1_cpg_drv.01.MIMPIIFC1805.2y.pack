SUBROUTINE LASCAW_VINTW(KPROM,KFLEV,KST,KPROF,&
 & LDCOMADV,LDT_SLHD,LDSLHDHEAT,PSLHDKMIN,&
 & KLEV,PLEV,PDVER,PDVERMAD,PSTDDISW,PKAPPA,PKAPPAT,PVETA,PVCUICO,PVSLD,PVSLDW,&
 & PVINTW,PVINTWMAD,PVINTWSLD,PVINTWSLT)

!     ------------------------------------------------------------------

!**** *LASCAW_VINTW  -  Weights for semi-LAgrangian interpolator:
!                       Computes PVINTW, PVINTWMAD and PVINTWSLD/T for one layer
!                       (high-order vertical weights)

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LASCAW_VINTW( ... )

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KPROM    - horizontal dimension.
!          KFLEV    - vertical dimension.
!          KST      - first element of arrays where computations are performed.
!          KPROF    - depth of work.
!          LDCOMADV - key for COMAD.
!          LDT_SLHD - key for SLHD.
!          PSLHDKMIN - either HOISLV or SLHDKMIN
!          KLEV     - lower level of the vertical interpolation
!                     grid needed for vertical interpolations.
!          PLEV     - vertical coordinate of the interpolation point.
!          PDVER    - distance for vertical linear interpolation.
!          PDVERMAD - PDVER for COMAD
!          PSTDDISW - STDDISW correction coef. for vertical COMAD
!          PKAPPA   - kappa function ("coefficient of SLHD").
!          PKAPPAT  - kappa function ("coefficient of SLHD") on T.
!          PVETA    - Values of ETA.
!          PVCUICO - Denominators of the vertical cubic interpolation coef.
!          PVSLD   - auxiliary quantities for vertical SLHD interpolation.
!          PVSLDW  - weights for SLHD vertical Laplacian smoother.

!        OUTPUT:
!          PVINTW    - vertical cubic interpolation weights.
!          PVINTWMAD - vertical cubic interpolation weights, COMAD case.
!          PVINTWSLD - vertical cubic interpolation weights, SLHD case.
!          PVINTWSLT - vertical cubic interpolation weights, SLHD case on T.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------

!        No external.
!        Called by some (E)LASCAW.. routines.

!     Reference.
!     ----------

!     Author.
!     -------
!        K. YESSAD, after former LASCAW code (JAN 2009).
!        METEO-FRANCE, CNRM/GMAP.

!     Modifications.
!     --------------
!        S. Malardel (Nov 2013): COMAD weights for SL interpolations
!        F. Vana 13-feb-2014 SLHD weights for heat variables
!        K. Yessad (March 2017): simplify level numbering.
!        F. Vana    21-Nov-2017: Option LHOISLT
!        H Petithomme (March 2021): optimization (test hoisting, access lowering)
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM,JPIA,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK 

! arp/ifs dependencies to be solved later.
USE YOMDYNA  , ONLY : SLHDKMAX,SLHDKREF
USE EINT_MOD , ONLY : JPDUP

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KPROM
INTEGER(KIND=JPIM), INTENT(IN)  :: KFLEV
INTEGER(KIND=JPIM), INTENT(IN)  :: KST
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROF
LOGICAL           , INTENT(IN)  :: LDCOMADV
LOGICAL           , INTENT(IN)  :: LDT_SLHD(3)
LOGICAL           , INTENT(IN)  :: LDSLHDHEAT
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLHDKMIN
INTEGER(KIND=JPIM), INTENT(IN)  :: KLEV(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PLEV(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PDVER(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PDVERMAD(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKAPPA(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSTDDISW(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKAPPAT(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PVETA(0:KFLEV+1)
REAL(KIND=JPRB)   , INTENT(IN)  :: PVCUICO(JPDUP,4,0:KFLEV-1)
REAL(KIND=JPRB)   , INTENT(IN)  :: PVSLD(JPDUP,3,0:KFLEV-1)
REAL(KIND=JPRB)   , INTENT(IN)  :: PVSLDW(JPDUP,3,3,0:KFLEV-1)
REAL(KIND=JPRB)   , INTENT(OUT) :: PVINTW(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PVINTWMAD(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PVINTWSLD(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PVINTWSLT(KPROM,KFLEV,3)

! optim: specific option for NEC machines, constant dim#1 (=1) otherwise
#ifndef NECSX
INTEGER(KIND=JPIA),PARAMETER :: I64=1
INTEGER(KIND=JPIM),PARAMETER :: IJ_=1
#endif
INTEGER(KIND=JPIA) :: ILEV64
INTEGER(KIND=JPIM) :: JROF,ILEV,JLEV
LOGICAL :: LLSLHD,LLSLHDQUAD
REAL(KIND=JPRB) :: ZWH1,ZWH2,ZWH3,ZWD1,ZWD2,ZWD3,ZWDS1,ZWDS2,ZWDS3,ZWL1,ZWL2,ZWL3
REAL(KIND=JPRB) :: ZD0,ZD1,ZD2,ZD3,ZDV,Z1,Z2,ZSLHDKMIN,ZKAP,ZKAPT
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('LASCAW_VINTW',0,ZHOOK_HANDLE)

LLSLHD=LDT_SLHD(1)
LLSLHDQUAD=LDT_SLHD(2)

! * Calculation of PVINTW and PVINTWSLD/T:
DO JLEV=1,KFLEV
  ! optim: manual hoisting of SLHD and COMAD options helps compiler optimizer
  IF (LLSLHD) THEN
    ! branch SLHD, old or quad scheme + varying SLD/SLT coefs
    ! optim: does not vectorize
    !CDIR NODEP
    DO JROF=KST,KPROF
      ILEV=KLEV(JROF,JLEV)
      IF (ILEV >= 1.AND.ILEV <= KFLEV-3) THEN
#ifdef NECSX
        IJ_ = MOD(JROF+1-KST,JPDUP)+1
        I64 = IJ_
#endif
        Z1=PVETA(ILEV+1)
        Z2=PVETA(ILEV+2)
        ZD0=PLEV(JROF,JLEV)-PVETA(ILEV)
        ZD1=PLEV(JROF,JLEV)-Z1
        ZD2=PLEV(JROF,JLEV)-Z2
        ZD3=PLEV(JROF,JLEV)-PVETA(ILEV+3)
        ILEV64=ILEV
        ZDV=ZD0*ZD1
        ZWH1=ZD0*ZD2*ZD3*PVCUICO(I64,2,ILEV64)
        ZWH2=ZDV*ZD3*PVCUICO(I64,3,ILEV64)
        ZWH3=ZDV*ZD2*PVCUICO(I64,4,ILEV64)

        ! warning: LSLHD and LSLHDQUAD can both be true, so no test simplification
        ! in this branch LLSLHD, alternative means old (linear) SLHD scheme
        IF (LLSLHDQUAD) THEN
          ZWL2=PDVER(JROF,JLEV)
          ZWL1=1.0_JPRB-ZWL2
          ZWL3=PVSLD(IJ_,3,ILEV)*ZWL1*ZWL2
          ZWL1=ZWL1+PVSLD(IJ_,1,ILEV)*ZWL3
          ZWL2=ZWL2+PVSLD(IJ_,2,ILEV)*ZWL3
        ELSE
          ZWL2=PDVER(JROF,JLEV)
          ZWL1=1.0_JPRB-ZWL2
          ZWL3=0.0_JPRB
        ENDIF

        ZWD1=ZWH1+SLHDKMAX*(ZWL1-ZWH1)
        ZWD2=ZWH2+SLHDKMAX*(ZWL2-ZWH2)
        ZWD3=ZWH3+SLHDKMAX*(ZWL3-ZWH3)
        ZWDS1=PVSLDW(I64,1,1,ILEV64)*ZWD1+PVSLDW(I64,1,2,ILEV64)*ZWD2+PVSLDW(I64,1,3,ILEV64)*ZWD3
        ZWDS2=PVSLDW(I64,2,1,ILEV64)*ZWD1+PVSLDW(I64,2,2,ILEV64)*ZWD2+PVSLDW(I64,2,3,ILEV64)*ZWD3
        ZWDS3=PVSLDW(I64,3,1,ILEV64)*ZWD1+PVSLDW(I64,3,2,ILEV64)*ZWD2+PVSLDW(I64,3,3,ILEV64)*ZWD3

        ! optim: no use of sign/abs functions anymore
        IF (PKAPPA(JROF,JLEV) < 0._JPRB) THEN
          ZKAP=-PKAPPA(JROF,JLEV)
          ZSLHDKMIN=SLHDKREF
        ELSE
          ZKAP=PKAPPA(JROF,JLEV)
          ZSLHDKMIN=PSLHDKMIN
        ENDIF

        ZWH1=ZWH1+ZSLHDKMIN*(ZWL1-ZWH1)
        ZWH2=ZWH2+ZSLHDKMIN*(ZWL2-ZWH2)
        ZWH3=ZWH3+ZSLHDKMIN*(ZWL3-ZWH3)
        PVINTW(JROF,JLEV,1)=ZWH1
        PVINTW(JROF,JLEV,2)=ZWH2
        PVINTW(JROF,JLEV,3)=ZWH3
        PVINTWSLD(JROF,JLEV,1)=ZWH1+ZKAP*(ZWDS1-ZWH1)
        PVINTWSLD(JROF,JLEV,2)=ZWH2+ZKAP*(ZWDS2-ZWH2)
        PVINTWSLD(JROF,JLEV,3)=ZWH3+ZKAP*(ZWDS3-ZWH3)

        IF (LDSLHDHEAT) THEN
          ZKAPT = ABS(PKAPPAT(JROF,JLEV))
          PVINTWSLT(JROF,JLEV,1)=ZWH1+ZKAPT*(ZWDS1-ZWH1)
          PVINTWSLT(JROF,JLEV,2)=ZWH2+ZKAPT*(ZWDS2-ZWH2)
          PVINTWSLT(JROF,JLEV,3)=ZWH3+ZKAPT*(ZWDS3-ZWH3)
        ENDIF
      ELSE
        PVINTW(JROF,JLEV,1)=1.0_JPRB-PDVER(JROF,JLEV)
        PVINTW(JROF,JLEV,2)=PDVER(JROF,JLEV)
        PVINTW(JROF,JLEV,3)=0.0_JPRB
        PVINTWSLD(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)
        PVINTWSLD(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)
        PVINTWSLD(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)

        IF (LDSLHDHEAT) THEN
          PVINTWSLT(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)
          PVINTWSLT(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)
          PVINTWSLT(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)
        ENDIF
      ENDIF
    ENDDO
  ELSE
    ! branch no SLHD or LLSLHDQUAD + constant SLD/SLT coefs
    ! optim: does not vectorize
    !CDIR NODEP
    DO JROF=KST,KPROF
      ILEV=KLEV(JROF,JLEV)
      IF (ILEV >= 1.AND.ILEV <= KFLEV-3) THEN
#ifdef NECSX
        IJ_ = MOD(JROF+1-KST,JPDUP)+1
        I64 = IJ_
#endif
        Z1=PVETA(ILEV+1)
        Z2=PVETA(ILEV+2)
        ZD0=PLEV(JROF,JLEV)-PVETA(ILEV)
        ZD1=PLEV(JROF,JLEV)-Z1
        ZD2=PLEV(JROF,JLEV)-Z2
        ZD3=PLEV(JROF,JLEV)-PVETA(ILEV+3)
        ILEV64=ILEV
        ZDV=ZD0*ZD1
        ZWH1=ZD0*ZD2*ZD3*PVCUICO(I64,2,ILEV64)
        ZWH2=ZDV*ZD3*PVCUICO(I64,3,ILEV64)
        ZWH3=ZDV*ZD2*PVCUICO(I64,4,ILEV64)

        IF (LLSLHDQUAD) THEN
          ZWL2=PDVER(JROF,JLEV)
          ZWL1=1.0_JPRB-ZWL2
          ZWL3=PVSLD(IJ_,3,ILEV)*ZWL1*ZWL2
          ZWL1=ZWL1+PVSLD(IJ_,1,ILEV)*ZWL3
          ZWL2=ZWL2+PVSLD(IJ_,2,ILEV)*ZWL3

          ZWH1=ZWH1+PSLHDKMIN*(ZWL1-ZWH1)
          ZWH2=ZWH2+PSLHDKMIN*(ZWL2-ZWH2)
          ZWH3=ZWH3+PSLHDKMIN*(ZWL3-ZWH3)
        ENDIF

        PVINTW(JROF,JLEV,1)=ZWH1
        PVINTW(JROF,JLEV,2)=ZWH2
        PVINTW(JROF,JLEV,3)=ZWH3
      ELSE
        PVINTW(JROF,JLEV,1)=1.0_JPRB-PDVER(JROF,JLEV)
        PVINTW(JROF,JLEV,2)=PDVER(JROF,JLEV)
        PVINTW(JROF,JLEV,3)=0.0_JPRB
      ENDIF
    ENDDO

    IF (LDSLHDHEAT) THEN
      !CDIR NODEP
      DO JROF=KST,KPROF
        PVINTWSLD(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)
        PVINTWSLD(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)
        PVINTWSLD(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)
        PVINTWSLT(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)
        PVINTWSLT(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)
        PVINTWSLT(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)
      ENDDO
    ELSE
      !CDIR NODEP
      DO JROF=KST,KPROF
        PVINTWSLD(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)
        PVINTWSLD(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)
        PVINTWSLD(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)
      ENDDO
    ENDIF
  ENDIF

  IF (LDCOMADV) THEN
    ! optim: does not vectorize
    !CDIR NODEP
    DO JROF=KST,KPROF
      ILEV=KLEV(JROF,JLEV)
      IF (ILEV >= 1.AND.ILEV <= KFLEV-3) THEN
        Z1 = PVETA(ILEV+1)
        Z2 = PVETA(ILEV+2)
        ZDV = 0.5_JPRB*(Z2-Z1)*(1._JPRB-PSTDDISW(JROF,JLEV))
        ZD1 = (PLEV(JROF,JLEV)-Z1)*PSTDDISW(JROF,JLEV)+ZDV
        ZD2 = (PLEV(JROF,JLEV)-Z2)*PSTDDISW(JROF,JLEV)-ZDV
        ZD0 = Z1-PVETA(ILEV)+ZD1
        ZD3 = Z2-PVETA(ILEV+3)+ZD2

        ILEV64 = ILEV
        ZDV=ZD0*ZD1
        PVINTWMAD(JROF,JLEV,1) = ZD0*ZD2*ZD3*PVCUICO(I64,2,ILEV64)
        PVINTWMAD(JROF,JLEV,2) = ZDV*ZD3*PVCUICO(I64,3,ILEV64)
        PVINTWMAD(JROF,JLEV,3) = ZDV*ZD2*PVCUICO(I64,4,ILEV64)
      ELSE
        PVINTWMAD(JROF,JLEV,1) = 1.0_JPRB-PDVERMAD(JROF,JLEV)
        PVINTWMAD(JROF,JLEV,2) = PDVERMAD(JROF,JLEV)
        PVINTWMAD(JROF,JLEV,3) = 0.0_JPRB
      ENDIF
    ENDDO
  ELSE
    !CDIR NODEP
    DO JROF=KST,KPROF
      PVINTWMAD(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)
      PVINTWMAD(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)
      PVINTWMAD(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)
    ENDDO
  ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('LASCAW_VINTW',1,ZHOOK_HANDLE)
END SUBROUTINE LASCAW_VINTW
