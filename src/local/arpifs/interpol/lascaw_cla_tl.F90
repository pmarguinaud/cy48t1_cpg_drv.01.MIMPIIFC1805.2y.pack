SUBROUTINE LASCAW_CLA_TL(YDSL,KFLEV,KPROM,KST,KPROF,LDT_SLHD,LDSLHDHEAT,&
 & PSLHDKMIN,&
 & KILA,PWH,PDLAT,PKAPPA,PKAPPAT,PSLD1_,PSLD2_,PSLD3_,PSLDW_,&
 & PWH5,PDLAT5,PKAPPA5,PKAPPAT5,&
 & PCLA,PCLASLD,PCLASLT,PCLA5,PCLASLD5,PCLASLT5)

!     ------------------------------------------------------------------

!**** *LASCAW_CLA_TL  -  Weights for semi-LAgrangian interpolator:
!                        Computes PCLA and PCLASLD for one layer
!                        (high-order zonal weights)
!                        TL + trajectory code
!      Spherical geometry only.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LASCAW_CLA_TL( ... )

!        Explicit arguments :
!        --------------------

!        INPUT:
!          YDSL     - SL_STRUCT definition
!          KFLEV    - Vertical dimension
!          KPROM    - horizontal dimension.
!          KST      - first element of arrays where computations are performed.
!          KPROF    - depth of work.
!          LDT_SLHD - keys for SLHD.
!          LDSLHDHEAT   - If true, the triggering function for heat variables differs from the one for momentum variables
!          PSLHDKMIN - either HOISLH or SLHDKMIN
!          KILA     - cf. ILA in LASCAW.
!          PWH      - cf. ZZWH in LASCAW.
!          PDLAT    - distance for horizontal linear interpolations in latitude
!          PKAPPA   - kappa function ("coefficient of SLHD").
!          PSLD1_   - auxiliary quantity for SLHD interpolation in latitude
!          PSLD2_   - auxiliary quantity for SLHD interpolation in latitude
!          PSLD3_   - auxiliary quantity for SLHD interpolation in latitude
!          PSLDW_   - weights for SLHD Laplacian smoother in latitude
!          PWH5     - cf. PWH (trajectory).
!          PDLAT5   - cf. PDLAT (trajectory).
!          PKAPPA5  - cf. PKAPPA (trajectory).

!        OUTPUT:
!          PCLA    - weights for horizontal cubic interpolations in latitude.
!          PCLASLD - cf. PCLA, SLHD case.
!          PCLASLT - cf. PCLA, SLHD case on T
!          PCLA5   - cf. PCLA (trajectory).
!          PCLASLD5- cf. PCLASLD (trajectory).
!          PCLASLT5- cf. PCLASLT (trajectory).

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------

!        No external.
!        Called by LASCAWTL.

!     Reference.
!     ----------

!     Author.
!     -------
!        K. YESSAD, after former LASCAWTL code (JAN 2009).
!        METEO-FRANCE, CNRM/GMAP.

!     Modifications.
!     --------------
!      G. Mozdzynski (May 2012): further cleaning
!      F. Vana 13-feb-2014 SLHD weights for heat variables
!      F. Vana 21-Nov-2017: Option LHOISLT
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK 

! arp/ifs dependencies to be solved later.
USE YOMDYNA  , ONLY : SLHDKMAX,SLHDKREF

USE EINT_MOD , ONLY : SL_STRUCT,JPDUP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(SL_STRUCT),    INTENT(IN)  :: YDSL
INTEGER(KIND=JPIM), INTENT(IN)  :: KFLEV
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROM
INTEGER(KIND=JPIM), INTENT(IN)  :: KST
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROF
LOGICAL           , INTENT(IN)  :: LDT_SLHD(3)
LOGICAL           , INTENT(IN)  :: LDSLHDHEAT
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLHDKMIN
INTEGER(KIND=JPIM), INTENT(IN)  :: KILA(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PWH(KPROM,3,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PDLAT(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKAPPA(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKAPPAT(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLD1_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLD2_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLD3_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   , INTENT(IN)  :: PSLDW_(JPDUP,3,3,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   , INTENT(IN)  :: PWH5(KPROM,3,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PDLAT5(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKAPPA5(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)  :: PKAPPAT5(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLA(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLASLD(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLASLT(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLA5(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLASLD5(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(OUT) :: PCLASLT5(KPROM,KFLEV,3)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,IJ_,ILA,JLEV
REAL(KIND=JPRB) :: ZWA1,ZWA2,ZWA3,ZWD1,ZWD2,ZWD3,ZWH1,ZWH2,ZWH3
REAL(KIND=JPRB) :: ZWDS1,ZWDS2,ZWDS3,ZWL1,ZWL2,ZWL3
REAL(KIND=JPRB) :: ZWA15,ZWA25,ZWA35,ZWD15,ZWD25,ZWD35,ZWH15,ZWH25,ZWH35
REAL(KIND=JPRB) :: ZWDS15,ZWDS25,ZWDS35,ZWL15,ZWL25,ZWL35
REAL(KIND=JPRB) :: ZSIGN, ZSLHDKMIN
LOGICAL :: LLSLHD,LLSLHDQUAD,LLSLHD_OLD
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_CLA_TL',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

LLSLHD=LDT_SLHD(1)
LLSLHDQUAD=LDT_SLHD(2)
LLSLHD_OLD=LDT_SLHD(3)

! * Calculation of PCLA and PCLASLD/T, PCLA5 and PCLASL[D/T]5:
DO JLEV=1,KFLEV
!CDIR NODEP
DO JROF=KST,KPROF

  ! a/ Preliminary calculations:

  IJ_ = MOD(JROF+1-KST,JPDUP)+1
  ILA=KILA(JROF,JLEV)

  ! b/ Trajectory code (in particuliar computes PCLA5 and PCLASLD/T5):

#if defined(NECSX)
  ! scalar variable initialization (to meet vectorization condition)
  ZWA15=0._JPRB
  ZWA25=0._JPRB
  ZWA35=0._JPRB
  ZWD15=0._JPRB
  ZWD25=0._JPRB
  ZWD35=0._JPRB
  ZWDS15=0._JPRB
  ZWDS25=0._JPRB
  ZWDS35=0._JPRB
#endif

  ZWH15=PWH5(JROF,1,JLEV)
  ZWH25=PWH5(JROF,2,JLEV)
  ZWH35=PWH5(JROF,3,JLEV)
  IF (LLSLHDQUAD) THEN
    ZWL25=PDLAT5(JROF,JLEV)
    ZWL15=1.0_JPRB-ZWL25
    ZWL35=PSLD3_(IJ_,ILA+1)*ZWL15*ZWL25
    ZWL15=ZWL15+PSLD1_(IJ_,ILA+1)*ZWL35
    ZWL25=ZWL25+PSLD2_(IJ_,ILA+1)*ZWL35
  ELSEIF (LLSLHD_OLD) THEN
    ZWL25=PDLAT5(JROF,JLEV)
    ZWL15=1.0_JPRB-ZWL25
    ZWL35=0.0_JPRB
  ENDIF

  IF (LLSLHD) THEN
    ZSIGN=SIGN(0.5_JPRB,PKAPPA5(JROF,JLEV))
    ZSLHDKMIN=(0.5_JPRB+ZSIGN)*PSLHDKMIN - (ZSIGN-0.5_JPRB)*SLHDKREF
    ZWA15=ZWH15+ZSLHDKMIN*(ZWL15-ZWH15)
    ZWA25=ZWH25+ZSLHDKMIN*(ZWL25-ZWH25)
    ZWA35=ZWH35+ZSLHDKMIN*(ZWL35-ZWH35)
    ZWD15=ZWH15+SLHDKMAX*(ZWL15-ZWH15)
    ZWD25=ZWH25+SLHDKMAX*(ZWL25-ZWH25)
    ZWD35=ZWH35+SLHDKMAX*(ZWL35-ZWH35)
    ZWDS15=PSLDW_(IJ_,1,1,ILA+1)*ZWD15+PSLDW_(IJ_,1,2,ILA+1)*ZWD25+&
     & PSLDW_(IJ_,1,3,ILA+1)*ZWD35
    ZWDS25=PSLDW_(IJ_,2,1,ILA+1)*ZWD15+PSLDW_(IJ_,2,2,ILA+1)*ZWD25+&
     & PSLDW_(IJ_,2,3,ILA+1)*ZWD35
    ZWDS35=PSLDW_(IJ_,3,1,ILA+1)*ZWD15+PSLDW_(IJ_,3,2,ILA+1)*ZWD25+&
     & PSLDW_(IJ_,3,3,ILA+1)*ZWD35
    PCLA5(JROF,JLEV,1)=ZWA15
    PCLA5(JROF,JLEV,2)=ZWA25
    PCLA5(JROF,JLEV,3)=ZWA35
    PCLASLD5(JROF,JLEV,1)=ZWA15+ABS(PKAPPA5(JROF,JLEV))*(ZWDS15-ZWA15)
    PCLASLD5(JROF,JLEV,2)=ZWA25+ABS(PKAPPA5(JROF,JLEV))*(ZWDS25-ZWA25)
    PCLASLD5(JROF,JLEV,3)=ZWA35+ABS(PKAPPA5(JROF,JLEV))*(ZWDS35-ZWA35)
    IF (LDSLHDHEAT) THEN
      PCLASLT5(JROF,JLEV,1)=ZWA15+ABS(PKAPPAT5(JROF,JLEV))*(ZWDS15-ZWA15)
      PCLASLT5(JROF,JLEV,2)=ZWA25+ABS(PKAPPAT5(JROF,JLEV))*(ZWDS25-ZWA25)
      PCLASLT5(JROF,JLEV,3)=ZWA35+ABS(PKAPPAT5(JROF,JLEV))*(ZWDS35-ZWA35)
    ENDIF
  ELSEIF (LLSLHDQUAD) THEN
    ZWA15=ZWH15+PSLHDKMIN*(ZWL15-ZWH15)
    ZWA25=ZWH25+PSLHDKMIN*(ZWL25-ZWH25)
    ZWA35=ZWH35+PSLHDKMIN*(ZWL35-ZWH35)
    PCLA5(JROF,JLEV,1)=ZWA15
    PCLA5(JROF,JLEV,2)=ZWA25
    PCLA5(JROF,JLEV,3)=ZWA35
  ELSE
    PCLA5(JROF,JLEV,1)=ZWH15
    PCLA5(JROF,JLEV,2)=ZWH25
    PCLA5(JROF,JLEV,3)=ZWH35
  ENDIF

  ! c/ TL (in particuliar computes PCLA and PCLASLT/D):

#if defined(NECSX)
  ! scalar variable initialization (to meet vectorization condition)
  ZWA1=0._JPRB
  ZWA2=0._JPRB
  ZWA3=0._JPRB
  ZWD1=0._JPRB
  ZWD2=0._JPRB
  ZWD3=0._JPRB
  ZWDS1=0._JPRB
  ZWDS2=0._JPRB
  ZWDS3=0._JPRB
#endif

  ZWH1=PWH(JROF,1,JLEV)
  ZWH2=PWH(JROF,2,JLEV)
  ZWH3=PWH(JROF,3,JLEV)
  IF (LLSLHDQUAD) THEN
    ZWL2=PDLAT(JROF,JLEV)
    ZWL1=-ZWL2
    ZWL3=PSLD3_(IJ_,ILA+1)*((ZWL1-ZWL2)*PDLAT5(JROF,JLEV) + ZWL2)
    ZWL1=ZWL1+PSLD1_(IJ_,ILA+1)*ZWL3
    ZWL2=ZWL2+PSLD2_(IJ_,ILA+1)*ZWL3
  ELSEIF (LLSLHD_OLD) THEN
    ZWL2=PDLAT(JROF,JLEV)
    ZWL1=-ZWL2
    ZWL3=0.0_JPRB
  ENDIF

  IF (LLSLHD) THEN
    ZWA1=ZWH1+ZSLHDKMIN*(ZWL1-ZWH1)
    ZWA2=ZWH2+ZSLHDKMIN*(ZWL2-ZWH2)
    ZWA3=ZWH3+ZSLHDKMIN*(ZWL3-ZWH3)
    ZWD1=ZWH1+SLHDKMAX*(ZWL1-ZWH1)
    ZWD2=ZWH2+SLHDKMAX*(ZWL2-ZWH2)
    ZWD3=ZWH3+SLHDKMAX*(ZWL3-ZWH3)
    ZWDS1=PSLDW_(IJ_,1,1,ILA+1)*ZWD1+PSLDW_(IJ_,1,2,ILA+1)*ZWD2+&
     & PSLDW_(IJ_,1,3,ILA+1)*ZWD3
    ZWDS2=PSLDW_(IJ_,2,1,ILA+1)*ZWD1+PSLDW_(IJ_,2,2,ILA+1)*ZWD2+&
     & PSLDW_(IJ_,2,3,ILA+1)*ZWD3
    ZWDS3=PSLDW_(IJ_,3,1,ILA+1)*ZWD1+PSLDW_(IJ_,3,2,ILA+1)*ZWD2+&
     & PSLDW_(IJ_,3,3,ILA+1)*ZWD3
    PCLA(JROF,JLEV,1)=ZWA1
    PCLA(JROF,JLEV,2)=ZWA2
    PCLA(JROF,JLEV,3)=ZWA3
    PCLASLD(JROF,JLEV,1)=ZWA1+ABS(PKAPPA5(JROF,JLEV))*(ZWDS1-ZWA1)&
     & +PKAPPA(JROF,JLEV)*(ZWDS15-ZWA15)
    PCLASLD(JROF,JLEV,2)=ZWA2+ABS(PKAPPA5(JROF,JLEV))*(ZWDS2-ZWA2)&
     & +PKAPPA(JROF,JLEV)*(ZWDS25-ZWA25)
    PCLASLD(JROF,JLEV,3)=ZWA3+ABS(PKAPPA5(JROF,JLEV))*(ZWDS3-ZWA3)&
     & +PKAPPA(JROF,JLEV)*(ZWDS35-ZWA35)
    IF (LDSLHDHEAT) THEN
      PCLASLT(JROF,JLEV,1)=ZWA1+ABS(PKAPPAT5(JROF,JLEV))*(ZWDS1-ZWA1)&
       & +PKAPPAT(JROF,JLEV)*(ZWDS15-ZWA15)
      PCLASLT(JROF,JLEV,2)=ZWA2+ABS(PKAPPAT5(JROF,JLEV))*(ZWDS2-ZWA2)&
       & +PKAPPAT(JROF,JLEV)*(ZWDS25-ZWA25)
      PCLASLT(JROF,JLEV,3)=ZWA3+ABS(PKAPPAT5(JROF,JLEV))*(ZWDS3-ZWA3)&
       & +PKAPPAT(JROF,JLEV)*(ZWDS35-ZWA35)
    ENDIF
  ELSEIF (LLSLHDQUAD) THEN
    ZWA1=ZWH1+PSLHDKMIN*(ZWL1-ZWH1)
    ZWA2=ZWH2+PSLHDKMIN*(ZWL2-ZWH2)
    ZWA3=ZWH3+PSLHDKMIN*(ZWL3-ZWH3)
    PCLA(JROF,JLEV,1)=ZWA1
    PCLA(JROF,JLEV,2)=ZWA2
    PCLA(JROF,JLEV,3)=ZWA3
  ELSE
    PCLA(JROF,JLEV,1)=ZWH1
    PCLA(JROF,JLEV,2)=ZWH2
    PCLA(JROF,JLEV,3)=ZWH3
  ENDIF

ENDDO
ENDDO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_CLA_TL',1,ZHOOK_HANDLE)
END SUBROUTINE LASCAW_CLA_TL
