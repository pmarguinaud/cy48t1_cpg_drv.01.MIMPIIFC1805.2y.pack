SUBROUTINE LASCAW_CLO(KFLEV,KPROM,KST,KPROF,LDT,PSLHDKMIN,PDLO,PDLOMAD,PKAPPA,&
 & PCLO,PCLOMAD,PCLOSLD)

!     ------------------------------------------------------------------

!**** *LASCAW_CLO  -  Weights for semi-LAgrangian interpolator:
!                     Computes PCLO, PCLOMAD and PCLOSLD for one layer
!                     (high-order meridian weights)

!      Can be also called for zonal high-order weights if plane geometry
!       (no need to code a specific ELASCAW_CLA).

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LASCAW_CLO( ... )

!        Explicit arguments :
!        --------------------

!        * INPUT:
!        KFLEV    - Vertical dimension
!        KPROM    - horizontal dimension.
!        KST      - first element of arrays where computations are performed.
!        KPROF    - depth of work.
!        LDT      - keys for SLHD and horizontal turbulence.
!        PSLHDKMIN - either HOISLH or SLHDKMIN
!        PDLO     - distances for horizontal linear interpolations in longitude.
!        PDLOMAD  -  PDLO for COMAD
!        PKAPPA   - kappa function ("coefficient of SLHD").

!        * OUTPUT:
!        PCLO     - weights for horizontal cubic interpolations in longitude.
!        PCLOMAD  - cf. PCLO, COMAD case.
!        PCLOSLD  - cf. PCLO, SLHD case.

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
!     F. Vana 22-Feb-2011: Horiz. turbulence and diff on phys. tendencies
!     S. Malardel (Nov 2013): COMAD weights for SL interpolations
!     F. Vana 13-feb-2014 SLHD weights for heat variables
!     K. Yessad (March 2017): simplify level numbering.
!     F. Vana    21-Nov-2017: Option LHOISLT
!     H Petithomme (Dec 2020): optimization, option 3DTURB moved out
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK

! arp/ifs dependencies to be solved later.
USE YOMDYNA  , ONLY : YRDYNA

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KFLEV
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROM
INTEGER(KIND=JPIM), INTENT(IN)  :: KST
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROF
LOGICAL           , INTENT(IN)  :: LDT(4)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLHDKMIN
REAL(KIND=JPRB)   , INTENT(IN)  :: PDLO(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PDLOMAD(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKAPPA(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLO(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLOMAD(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLOSLD(KPROM,KFLEV,3)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,JLEV
REAL(KIND=JPRB) :: ZWH1,ZWH2,ZWH3,ZWD1,ZWD2,ZWD3,ZWDS1,ZWDS2,ZWDS3,ZWL1,ZWL2,ZWL3
REAL(KIND=JPRB) :: Z1M2EPSH,ZSIGN, ZSLHDKMIN, ZKAP
LOGICAL :: LLSLHD,LLSLHDQUAD,LLSLHD_OLD
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

REAL(KIND=JPRB),PARAMETER :: PP6_R=1.0_JPRB/6.0_JPRB
REAL(KIND=JPRB) :: PD,PD12
REAL(KIND=JPRB) :: FLAG1, FLAG2, FLAG3,F12
REAL(KIND=JPRB) :: FQUAD1, FQUAD2, FQUAD3

! weights for cubic Lagrange interpolation (regular nodes)
! optim: intermediate computation PD12=F12(PD) gains 1 multiply (fully reproducible form)
F12(PD) = 0.5_JPRB*(PD+1.0_JPRB)
FLAG1(PD,PD12)=PD12*(PD-1.0_JPRB)*(PD-2.0_JPRB)
FLAG2(PD,PD12)=-PD12*PD*(PD-2.0_JPRB)
FLAG3(PD)= PP6_R   *(PD+1.0_JPRB)*PD*(PD-1.0_JPRB)

! weights for quadratic SLHD interpolation (regular nodes)
FQUAD1(PD)=(1.0_JPRB-PD)*(1.0_JPRB+0.25_JPRB*PD)
FQUAD2(PD)=PD*(1.25_JPRB-0.25_JPRB*PD)
FQUAD3(PD)=0.25_JPRB*PD*(PD-1.0_JPRB)

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_CLO',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

LLSLHD=LDT(1)
LLSLHDQUAD=LDT(2)
LLSLHD_OLD=LDT(3)

! * Auxiliary quantity for Laplacian smoother:
Z1M2EPSH=1.0_JPRB-2.0_JPRB*YRDYNA%SLHDEPSH

! * Calculation of PCLO and PCLOSLD:
DO JLEV=1,KFLEV
  !CDIR NODEP
  DO JROF=KST,KPROF
    PD12=F12(PDLO(JROF,JLEV))
    ZWH1=FLAG1(PDLO(JROF,JLEV),PD12)
    ZWH2=FLAG2(PDLO(JROF,JLEV),PD12)
    ZWH3=FLAG3(PDLO(JROF,JLEV))

    ! warning: LSLHD and LSLHDQUAD can both be true, so no test simplification
    ! or reordering. On the contrary, LSLHDQUAD and LSLHD_OLD can not both be true.
    IF (LLSLHDQUAD) THEN
      ZWL1=FQUAD1(PDLO(JROF,JLEV))
      ZWL2=FQUAD2(PDLO(JROF,JLEV))
      ZWL3=FQUAD3(PDLO(JROF,JLEV))
    ENDIF

    IF (LLSLHD) THEN
      IF (LLSLHD_OLD) THEN
        ZWL2=PDLO(JROF,JLEV)
        ZWL1=1.0_JPRB-ZWL2
        ZWL3=0.0_JPRB
      ENDIF

      ZWD1=ZWH1+YRDYNA%SLHDKMAX*(ZWL1-ZWH1)
      ZWD2=ZWH2+YRDYNA%SLHDKMAX*(ZWL2-ZWH2)
      ZWD3=ZWH3+YRDYNA%SLHDKMAX*(ZWL3-ZWH3)
      ZWDS1=Z1M2EPSH*ZWD1+YRDYNA%SLHDEPSH*ZWD2
      ZWDS2=YRDYNA%SLHDEPSH*ZWD1+Z1M2EPSH*ZWD2
      ZWDS3=YRDYNA%SLHDEPSH*ZWD2+ZWD3

      IF (PKAPPA(JROF,JLEV) < 0._JPRB) THEN
        ZSLHDKMIN=YRDYNA%SLHDKREF
        ZKAP=-PKAPPA(JROF,JLEV)
      ELSE
        ZSLHDKMIN=PSLHDKMIN
        ZKAP=PKAPPA(JROF,JLEV)
      ENDIF

      ZWH1=ZWH1+ZSLHDKMIN*(ZWL1-ZWH1)
      ZWH2=ZWH2+ZSLHDKMIN*(ZWL2-ZWH2)
      ZWH3=ZWH3+ZSLHDKMIN*(ZWL3-ZWH3)

      PCLO(JROF,JLEV,1)=ZWH1
      PCLO(JROF,JLEV,2)=ZWH2
      PCLO(JROF,JLEV,3)=ZWH3
      PCLOSLD(JROF,JLEV,1)=ZWH1+ZKAP*(ZWDS1-ZWH1)
      PCLOSLD(JROF,JLEV,2)=ZWH2+ZKAP*(ZWDS2-ZWH2)
      PCLOSLD(JROF,JLEV,3)=ZWH3+ZKAP*(ZWDS3-ZWH3)
    ELSEIF (LLSLHDQUAD) THEN
      ZWH1=ZWH1+PSLHDKMIN*(ZWL1-ZWH1)
      ZWH2=ZWH2+PSLHDKMIN*(ZWL2-ZWH2)
      ZWH3=ZWH3+PSLHDKMIN*(ZWL3-ZWH3)
      PCLO(JROF,JLEV,1)=ZWH1
      PCLO(JROF,JLEV,2)=ZWH2
      PCLO(JROF,JLEV,3)=ZWH3
      PCLOSLD(JROF,JLEV,1)=ZWH1
      PCLOSLD(JROF,JLEV,2)=ZWH2
      PCLOSLD(JROF,JLEV,3)=ZWH3
    ELSE
      PCLO(JROF,JLEV,1)=ZWH1
      PCLO(JROF,JLEV,2)=ZWH2
      PCLO(JROF,JLEV,3)=ZWH3
      PCLOSLD(JROF,JLEV,1)=ZWH1
      PCLOSLD(JROF,JLEV,2)=ZWH2
      PCLOSLD(JROF,JLEV,3)=ZWH3
    ENDIF
  ENDDO

  !CDIR NODEP
  DO JROF=KST,KPROF
    PD12=F12(PDLOMAD(JROF,JLEV))
    PCLOMAD(JROF,JLEV,1)=FLAG1(PDLOMAD(JROF,JLEV),PD12)
    PCLOMAD(JROF,JLEV,2)=FLAG2(PDLOMAD(JROF,JLEV),PD12)
    PCLOMAD(JROF,JLEV,3)=FLAG3(PDLOMAD(JROF,JLEV))
  ENDDO
ENDDO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_CLO',1,ZHOOK_HANDLE)
END SUBROUTINE LASCAW_CLO
