SUBROUTINE LASCAW_CLO_AD(KFLEV,KPROM,KST,KPROF,LDT_SLHD,LDSLHDHEAT,PSLHDKMIN,&
 & PDLO,PKAPPA,PKAPPAT,&
 & PDLO5,PKAPPA5,PKAPPAT5,PCLO,PCLOSLD,PCLOSLT)

!     ------------------------------------------------------------------

!**** *LASCAW_CLO_AD  -  Weights for semi-LAgrangian interpolator:
!                        Computes PCLO and PCLOSLD/T for one layer
!                        (high-order meridian weights)
!                        AD + trajectory code

!      Can be also called for zonal high-order weights if plane geometry
!       (no need to code a specific ELASCAW_CLA_AD).

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LASCAW_CLO_AD( ... )

!        Explicit arguments :
!        --------------------

!        KFLEV    - Vertical dimension
!        KPROM    - horizontal dimension.
!        KST      - first element of arrays where computations are performed.
!        KPROF    - depth of work.
!        LDT_SLHD - keys for SLHD.
!        LDSLHDHEAT   - If true, the triggering function for heat variables differs from the one for momentum variables
!        PSLHDKMIN - either HOISLH or SLHDKMIN
!        PDLO     - distances for horizontal linear interpolations in longitude.
!        PKAPPA   - kappa function ("coefficient of SLHD").
!        PKAPPAT  - kappa function ("coefficient of SLHD") on T.
!        PDLO5    - cf. PDLO (trajectory).
!        PKAPPA5  - cf. PKAPPA (trajectory).
!        PKAPPAT5 - cf. PKAPPAT (trajectory).

!        PCLO     - weights for horizontal cubic interpolations in longitude.
!        PCLOSLD  - cf. PCLO, SLHD case.
!        PCLOSLT  - cf. PCLO, SLHD case on T.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------

!        No external.
!        Called by some (E)LASCAWAD routines.

!     Reference.
!     ----------

!     Author.
!     -------
!        K. YESSAD and F. VANA, after former LASCAWAD code (JAN 2009).
!        METEO-FRANCE, CNRM/GMAP.

!     Modifications.
!     --------------
!      F. Vana 13-feb-2014 SLHD weights for heat variables
!      F. Vana 21-Nov-2017: Option LHOISLT
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK 

! arp/ifs dependencies to be solved later.
USE YOMDYNA  , ONLY : SLHDKMAX ,SLHDEPSH,SLHDKREF

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM), INTENT(IN)    :: KPROM
INTEGER(KIND=JPIM), INTENT(IN)    :: KST
INTEGER(KIND=JPIM), INTENT(IN)    :: KPROF
LOGICAL           , INTENT(IN)    :: LDT_SLHD(3)
LOGICAL           , INTENT(IN)    :: LDSLHDHEAT
REAL(KIND=JPRB)   , INTENT(IN)    :: PSLHDKMIN
REAL(KIND=JPRB)   , INTENT(INOUT) :: PDLO(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PKAPPA(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PKAPPAT(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)    :: PDLO5(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)    :: PKAPPA5(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)    :: PKAPPAT5(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)    :: PCLO(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(IN)    :: PCLOSLD(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(IN)    :: PCLOSLT(KPROM,KFLEV,3)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,JLEV
REAL(KIND=JPRB) :: ZWA1,ZWA2,ZWA3,ZWD1,ZWD2,ZWD3,ZWH1,ZWH2,ZWH3
REAL(KIND=JPRB) :: ZWDS1,ZWDS2,ZWDS3,ZWL1,ZWL2,ZWL3
REAL(KIND=JPRB) :: ZWA15,ZWA25,ZWA35,ZWD15,ZWD25,ZWD35,ZWH15,ZWH25,ZWH35
REAL(KIND=JPRB) :: ZWDS15,ZWDS25,ZWDS35,ZWL15,ZWL25,ZWL35
REAL(KIND=JPRB) :: Z1M2EPSH, ZSIGN, ZSLHDKMIN
LOGICAL :: LLSLHD,LLSLHDQUAD,LLSLHD_OLD
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

REAL(KIND=JPRB),PARAMETER :: PP6_R=1.0_JPRB/6.0_JPRB
REAL(KIND=JPRB) :: PD
REAL(KIND=JPRB) :: FLAG1, FLAG2, FLAG3, FDLAG1, FDLAG2, FDLAG3
REAL(KIND=JPRB) :: FQUAD1, FQUAD2, FQUAD3, FDQUAD1, FDQUAD2, FDQUAD3

! weights for cubic Lagrange interpolation (regular nodes)
FLAG1(PD)= 0.5_JPRB*(PD+1.0_JPRB)   *(PD-1.0_JPRB)*(PD-2.0_JPRB)
FLAG2(PD)=-0.5_JPRB*(PD+1.0_JPRB)*PD              *(PD-2.0_JPRB)
FLAG3(PD)= PP6_R   *(PD+1.0_JPRB)*PD*(PD-1.0_JPRB)
FDLAG1(PD)= 0.5_JPRB*(3.0_JPRB*PD*PD-4.0_JPRB*PD-1.0_JPRB)
FDLAG2(PD)=-0.5_JPRB*(3.0_JPRB*PD*PD-2.0_JPRB*PD-2.0_JPRB)
FDLAG3(PD)= PP6_R   *(3.0_JPRB*PD*PD-1.0_JPRB)

! weights for quadratic SLHD interpolation (regular nodes)
FQUAD1(PD)=(1.0_JPRB-PD)*(1.0_JPRB+0.25_JPRB*PD)
FQUAD2(PD)=PD*(1.25_JPRB-0.25_JPRB*PD)
FQUAD3(PD)=0.25_JPRB*PD*(PD-1.0_JPRB)
FDQUAD1(PD)=-0.5_JPRB*PD - 0.75_JPRB
FDQUAD2(PD)=-0.5_JPRB*PD + 1.25_JPRB
FDQUAD3(PD)= 0.5_JPRB*PD - 0.25_JPRB

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_CLO_AD',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

LLSLHD=LDT_SLHD(1)
LLSLHDQUAD=LDT_SLHD(2)
LLSLHD_OLD=LDT_SLHD(3)

! * Auxiliary quantity for Laplacian smoother:
Z1M2EPSH=1.0_JPRB-2.0_JPRB*SLHDEPSH

! * Calculation of PCLO and PCLOSLD/T:

IF (LLSLHD) THEN ! LLSLHD

  DO JLEV=1,KFLEV
!CDIR NODEP
  DO JROF=KST,KPROF

    ! a/ Trajectory code:

    ! trajectory
    ZWH15=FLAG1(PDLO5(JROF,JLEV))
    ZWH25=FLAG2(PDLO5(JROF,JLEV))
    ZWH35=FLAG3(PDLO5(JROF,JLEV))
    IF (LLSLHDQUAD) THEN
      ZWL15=FQUAD1(PDLO5(JROF,JLEV))
      ZWL25=FQUAD2(PDLO5(JROF,JLEV))
      ZWL35=FQUAD3(PDLO5(JROF,JLEV))
    ELSEIF (LLSLHD_OLD) THEN
      ZWL25=PDLO5(JROF,JLEV)
      ZWL15=1.0_JPRB-ZWL25
      ZWL35=0.0_JPRB
    ENDIF

    ZSIGN=SIGN(0.5_JPRB,PKAPPA5(JROF,JLEV))
    ZSLHDKMIN=(0.5_JPRB+ZSIGN)*PSLHDKMIN - (ZSIGN-0.5_JPRB)*SLHDKREF
    ZWA15=ZWH15+ZSLHDKMIN*(ZWL15-ZWH15)
    ZWA25=ZWH25+ZSLHDKMIN*(ZWL25-ZWH25)
    ZWA35=ZWH35+ZSLHDKMIN*(ZWL35-ZWH35)
    ZWD15=ZWH15+SLHDKMAX*(ZWL15-ZWH15)
    ZWD25=ZWH25+SLHDKMAX*(ZWL25-ZWH25)
    ZWD35=ZWH35+SLHDKMAX*(ZWL35-ZWH35)
    ZWDS15=Z1M2EPSH*ZWD15+SLHDEPSH*ZWD25
    ZWDS25=SLHDEPSH*ZWD15+Z1M2EPSH*ZWD25
    ZWDS35=SLHDEPSH*ZWD25+ZWD35

    ! b/ AD code:

    PKAPPA(JROF,JLEV)=PKAPPA(JROF,JLEV)&
     & +(ZWDS15-ZWA15)*PCLOSLD(JROF,JLEV,1)&
     & +(ZWDS25-ZWA25)*PCLOSLD(JROF,JLEV,2)&
     & +(ZWDS35-ZWA35)*PCLOSLD(JROF,JLEV,3)
    ZWDS1=ABS(PKAPPA5(JROF,JLEV))*PCLOSLD(JROF,JLEV,1)
    ZWDS2=ABS(PKAPPA5(JROF,JLEV))*PCLOSLD(JROF,JLEV,2)
    ZWDS3=ABS(PKAPPA5(JROF,JLEV))*PCLOSLD(JROF,JLEV,3)
    ZWA1=(1._JPRB-ABS(PKAPPA5(JROF,JLEV)))*PCLOSLD(JROF,JLEV,1)+PCLO(JROF,JLEV,1)
    ZWA2=(1._JPRB-ABS(PKAPPA5(JROF,JLEV)))*PCLOSLD(JROF,JLEV,2)+PCLO(JROF,JLEV,2)
    ZWA3=(1._JPRB-ABS(PKAPPA5(JROF,JLEV)))*PCLOSLD(JROF,JLEV,3)+PCLO(JROF,JLEV,3)
    IF (LDSLHDHEAT) THEN
      PKAPPAT(JROF,JLEV)=PKAPPAT(JROF,JLEV)&
       & +(ZWDS15-ZWA15)*PCLOSLT(JROF,JLEV,1)&
       & +(ZWDS25-ZWA25)*PCLOSLT(JROF,JLEV,2)&
       & +(ZWDS35-ZWA35)*PCLOSLT(JROF,JLEV,3)
      ZWDS1=ZWDS1+ABS(PKAPPAT5(JROF,JLEV))*PCLOSLT(JROF,JLEV,1)
      ZWDS2=ZWDS2+ABS(PKAPPAT5(JROF,JLEV))*PCLOSLT(JROF,JLEV,2)
      ZWDS3=ZWDS3+ABS(PKAPPAT5(JROF,JLEV))*PCLOSLT(JROF,JLEV,3)
      ZWA1=ZWA1 + (1._JPRB-ABS(PKAPPAT5(JROF,JLEV)))*PCLOSLT(JROF,JLEV,1)
      ZWA2=ZWA2 + (1._JPRB-ABS(PKAPPAT5(JROF,JLEV)))*PCLOSLT(JROF,JLEV,2)
      ZWA3=ZWA3 + (1._JPRB-ABS(PKAPPAT5(JROF,JLEV)))*PCLOSLT(JROF,JLEV,3)
    ENDIF
    ZWD1=Z1M2EPSH*ZWDS1+SLHDEPSH*ZWDS2
    ZWD2=SLHDEPSH*(ZWDS1+ZWDS3)+Z1M2EPSH*ZWDS2
    ZWD3=ZWDS3
    ZWH1=(1._JPRB-ZSLHDKMIN)*ZWA1+(1._JPRB-SLHDKMAX)*ZWD1
    ZWH2=(1._JPRB-ZSLHDKMIN)*ZWA2+(1._JPRB-SLHDKMAX)*ZWD2
    ZWH3=(1._JPRB-ZSLHDKMIN)*ZWA3+(1._JPRB-SLHDKMAX)*ZWD3
    ZWL1=ZSLHDKMIN*ZWA1+SLHDKMAX*ZWD1
    ZWL2=ZSLHDKMIN*ZWA2+SLHDKMAX*ZWD2
    ZWL3=ZSLHDKMIN*ZWA3+SLHDKMAX*ZWD3

    IF (LLSLHDQUAD) THEN
      PDLO(JROF,JLEV)=PDLO(JROF,JLEV)+FDQUAD1(PDLO5(JROF,JLEV))*ZWL1&
       & +FDQUAD2(PDLO5(JROF,JLEV))*ZWL2+FDQUAD3(PDLO5(JROF,JLEV))*ZWL3
    ELSEIF (LLSLHD_OLD) THEN
      PDLO(JROF,JLEV)=PDLO(JROF,JLEV)-ZWL1+ZWL2
    ENDIF
    PDLO(JROF,JLEV)=PDLO(JROF,JLEV)+FDLAG1(PDLO5(JROF,JLEV))*ZWH1&
     & +FDLAG2(PDLO5(JROF,JLEV))*ZWH2+FDLAG3(PDLO5(JROF,JLEV))*ZWH3

  ENDDO
  ENDDO

ELSE ! .NOT.LLSLHD

  DO JLEV=1,KFLEV
!CDIR NODEP
  DO JROF=KST,KPROF

    ! b/ AD code:

    IF (LLSLHDQUAD) THEN
      ZWA1=PCLO(JROF,JLEV,1)
      ZWA2=PCLO(JROF,JLEV,2)
      ZWA3=PCLO(JROF,JLEV,3)
      ZWH1=(1._JPRB-PSLHDKMIN)*ZWA1
      ZWH2=(1._JPRB-PSLHDKMIN)*ZWA2
      ZWH3=(1._JPRB-PSLHDKMIN)*ZWA3
      ZWL1=PSLHDKMIN*ZWA1
      ZWL2=PSLHDKMIN*ZWA2
      ZWL3=PSLHDKMIN*ZWA3
    ELSE
      ZWH1=PCLO(JROF,JLEV,1)
      ZWH2=PCLO(JROF,JLEV,2)
      ZWH3=PCLO(JROF,JLEV,3)
    ENDIF

    IF (LLSLHDQUAD) THEN
      PDLO(JROF,JLEV)=PDLO(JROF,JLEV)+FDQUAD1(PDLO5(JROF,JLEV))*ZWL1&
       & +FDQUAD2(PDLO5(JROF,JLEV))*ZWL2+FDQUAD3(PDLO5(JROF,JLEV))*ZWL3
    ENDIF
    PDLO(JROF,JLEV)=PDLO(JROF,JLEV)+FDLAG1(PDLO5(JROF,JLEV))*ZWH1&
     & +FDLAG2(PDLO5(JROF,JLEV))*ZWH2+FDLAG3(PDLO5(JROF,JLEV))*ZWH3

  ENDDO
  ENDDO

ENDIF ! Test on LLSLHD

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_CLO_AD',1,ZHOOK_HANDLE)
END SUBROUTINE LASCAW_CLO_AD
