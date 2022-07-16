SUBROUTINE LASCAW_CLA_AD(YDSL,KFLEV,KPROM,KST,KPROF,LDT_SLHD,LDSLHDHEAT,&
 & PSLHDKMIN,&
 & KILA,PWH,PDLAT,PKAPPA,PKAPPAT,PSLD1_,PSLD2_,PSLD3_,PSLDW_,&
 & PWH5,PDLAT5,PKAPPA5,PKAPPAT5,&
 & PCLA,PCLASLD,PCLASLT)

!     ------------------------------------------------------------------

!**** *LASCAW_CLA_AD  -  Weights for semi-LAgrangian interpolator:
!                        Computes PCLA and PCLASLD/T for one layer
!                        (high-order zonal weights)
!                        AD + trajectory code
!      Spherical geometry only.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LASCAW_CLA_AD( ... )

!        Explicit arguments :
!        --------------------

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
!          PKAPPAT  - kappa function ("coefficient of SLHD") for T.
!          PSLD1_   - auxiliary quantity for SLHD interpolation in latitude
!          PSLD2_   - auxiliary quantity for SLHD interpolation in latitude
!          PSLD3_   - auxiliary quantity for SLHD interpolation in latitude
!          PSLDW_   - weights for SLHD Laplacian smoother in latitude
!          PWH5     - cf. PWH (trajectory).
!          PDLAT5   - cf. PDLAT (trajectory).
!          PKAPPA5  - cf. PKAPPA (trajectory).
!          PKAPPAT5 - cf. PKAPPAT (trajectory).

!          PCLA    - weights for horizontal cubic interpolations in latitude.
!          PCLASLD - cf. PCLA, SLHD case.
!          PCLASLT - cf. PCLA, SLHD case on T.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about semi-Lagrangian scheme.

!     Externals.
!     ----------

!        No external.
!        Called by LASCAWAD.

!     Reference.
!     ----------

!     Author.
!     -------
!        K. YESSAD and F. VANA, after former LASCAWAD code (JAN 2009).
!        METEO-FRANCE, CNRM/GMAP.

!     Modifications.
!     --------------
!        G. Mozdzynski (May 2012): further cleaning
!      F. Vana 13-feb-2014 SLHD weights for heat variables
!      F. Vana 21-Nov-2017: Option LHOISLT
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK 

! arp/ifs dependencies to be solved later.
USE YOMDYNA  , ONLY : SLHDKMAX,SLHDKREF

USE EINT_MOD , ONLY : SL_STRUCT,JPDUP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(SL_STRUCT),    INTENT(IN)    :: YDSL
INTEGER(KIND=JPIM), INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM), INTENT(IN)    :: KPROM
INTEGER(KIND=JPIM), INTENT(IN)    :: KST
INTEGER(KIND=JPIM), INTENT(IN)    :: KPROF
LOGICAL           , INTENT(IN)    :: LDT_SLHD(3)
LOGICAL           , INTENT(IN)    :: LDSLHDHEAT
REAL(KIND=JPRB)   , INTENT(IN)    :: PSLHDKMIN
INTEGER(KIND=JPIM), INTENT(IN)    :: KILA(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(OUT)   :: PWH(KPROM,3,KFLEV)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PDLAT(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PKAPPA(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PKAPPAT(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)    :: PSLD1_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   , INTENT(IN)    :: PSLD2_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   , INTENT(IN)    :: PSLD3_(JPDUP,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   , INTENT(IN)    :: PSLDW_(JPDUP,3,3,YDSL%NDGSAH:YDSL%NDGENH)
REAL(KIND=JPRB)   , INTENT(IN)    :: PWH5(KPROM,3,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)    :: PDLAT5(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)    :: PKAPPA5(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)    :: PKAPPAT5(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)    :: PCLA(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(IN)    :: PCLASLD(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(IN)    :: PCLASLT(KPROM,KFLEV,3)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,IJ_,ILA,JLEV
REAL(KIND=JPRB) :: ZWA1,ZWA2,ZWA3,ZWD1,ZWD2,ZWD3
REAL(KIND=JPRB) :: ZWDS1,ZWDS2,ZWDS3,ZWL1,ZWL2,ZWL3
REAL(KIND=JPRB) :: ZWA15,ZWA25,ZWA35,ZWD15,ZWD25,ZWD35
REAL(KIND=JPRB) :: ZWDS15,ZWDS25,ZWDS35,ZWL15,ZWL25,ZWL35
REAL(KIND=JPRB) :: ZSIGN, ZSLHDKMIN
LOGICAL :: LLSLHD,LLSLHDQUAD,LLSLHD_OLD
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_CLA_AD',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

LLSLHD=LDT_SLHD(1)
LLSLHDQUAD=LDT_SLHD(2)
LLSLHD_OLD=LDT_SLHD(3)

! * AD of calculation of PCLA and PCLASLD/T:

IF (LLSLHD) THEN ! LLSLHD

  DO JLEV=1,KFLEV
!CDIR NODEP
  DO JROF=KST,KPROF

    ! a/ Preliminary calculations:

    IJ_ = MOD(JROF+1-KST,JPDUP)+1
    ILA=KILA(JROF,JLEV)

    ! b/ Trajectory code:

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

    ZSIGN=SIGN(0.5_JPRB,PKAPPA5(JROF,JLEV))
    ZSLHDKMIN=(0.5_JPRB+ZSIGN)*PSLHDKMIN - (ZSIGN-0.5_JPRB)*SLHDKREF
    ZWA15=PWH5(JROF,1,JLEV)+ZSLHDKMIN*(ZWL15-PWH5(JROF,1,JLEV))
    ZWA25=PWH5(JROF,2,JLEV)+ZSLHDKMIN*(ZWL25-PWH5(JROF,2,JLEV))
    ZWA35=PWH5(JROF,3,JLEV)+ZSLHDKMIN*(ZWL35-PWH5(JROF,3,JLEV))
    ZWD15=PWH5(JROF,1,JLEV)+SLHDKMAX*(ZWL15-PWH5(JROF,1,JLEV))
    ZWD25=PWH5(JROF,2,JLEV)+SLHDKMAX*(ZWL25-PWH5(JROF,2,JLEV))
    ZWD35=PWH5(JROF,3,JLEV)+SLHDKMAX*(ZWL35-PWH5(JROF,3,JLEV))
    ZWDS15=PSLDW_(IJ_,1,1,ILA+1)*ZWD15+PSLDW_(IJ_,1,2,ILA+1)*ZWD25+&
     & PSLDW_(IJ_,1,3,ILA+1)*ZWD35
    ZWDS25=PSLDW_(IJ_,2,1,ILA+1)*ZWD15+PSLDW_(IJ_,2,2,ILA+1)*ZWD25+&
     & PSLDW_(IJ_,2,3,ILA+1)*ZWD35
    ZWDS35=PSLDW_(IJ_,3,1,ILA+1)*ZWD15+PSLDW_(IJ_,3,2,ILA+1)*ZWD25+&
     & PSLDW_(IJ_,3,3,ILA+1)*ZWD35

    ! c/ AD code:

    PKAPPA(JROF,JLEV)=PKAPPA(JROF,JLEV)&
     & + (ZWDS15-ZWA15)*PCLASLD(JROF,JLEV,1)&
     & + (ZWDS25-ZWA25)*PCLASLD(JROF,JLEV,2)&
     & + (ZWDS35-ZWA35)*PCLASLD(JROF,JLEV,3)
    ZWA1=(1._JPRB-ABS(PKAPPA5(JROF,JLEV)))*PCLASLD(JROF,JLEV,1)+PCLA(JROF,JLEV,1) 
    ZWA2=(1._JPRB-ABS(PKAPPA5(JROF,JLEV)))*PCLASLD(JROF,JLEV,2)+PCLA(JROF,JLEV,2)
    ZWA3=(1._JPRB-ABS(PKAPPA5(JROF,JLEV)))*PCLASLD(JROF,JLEV,3)+PCLA(JROF,JLEV,3)
    ZWDS1=ABS(PKAPPA5(JROF,JLEV))*PCLASLD(JROF,JLEV,1)
    ZWDS2=ABS(PKAPPA5(JROF,JLEV))*PCLASLD(JROF,JLEV,2)
    ZWDS3=ABS(PKAPPA5(JROF,JLEV))*PCLASLD(JROF,JLEV,3)
    IF (LDSLHDHEAT) THEN
      PKAPPAT(JROF,JLEV)=PKAPPAT(JROF,JLEV)&
       & + (ZWDS15-ZWA15)*PCLASLT(JROF,JLEV,1)&
       & + (ZWDS25-ZWA25)*PCLASLT(JROF,JLEV,2)&
       & + (ZWDS35-ZWA35)*PCLASLT(JROF,JLEV,3)
      ZWA1=ZWA1 + (1._JPRB-ABS(PKAPPAT5(JROF,JLEV)))*PCLASLT(JROF,JLEV,1)
      ZWA2=ZWA2 + (1._JPRB-ABS(PKAPPAT5(JROF,JLEV)))*PCLASLT(JROF,JLEV,2)
      ZWA3=ZWA3 + (1._JPRB-ABS(PKAPPAT5(JROF,JLEV)))*PCLASLT(JROF,JLEV,3)
      ZWDS1=ZWDS1+ABS(PKAPPAT5(JROF,JLEV))*PCLASLT(JROF,JLEV,1)
      ZWDS2=ZWDS2+ABS(PKAPPAT5(JROF,JLEV))*PCLASLT(JROF,JLEV,2)
      ZWDS3=ZWDS3+ABS(PKAPPAT5(JROF,JLEV))*PCLASLT(JROF,JLEV,3)
    ENDIF
    ZWD1=PSLDW_(IJ_,1,1,ILA+1)*ZWDS1+PSLDW_(IJ_,2,1,ILA+1)*ZWDS2&
     & +PSLDW_(IJ_,3,1,ILA+1)*ZWDS3
    ZWD2=PSLDW_(IJ_,1,2,ILA+1)*ZWDS1+PSLDW_(IJ_,2,2,ILA+1)*ZWDS2&
     & +PSLDW_(IJ_,3,2,ILA+1)*ZWDS3
    ZWD3=PSLDW_(IJ_,1,3,ILA+1)*ZWDS1+PSLDW_(IJ_,2,3,ILA+1)*ZWDS2&
     & +PSLDW_(IJ_,3,3,ILA+1)*ZWDS3
    PWH(JROF,1,JLEV)=(1._JPRB-ZSLHDKMIN)*ZWA1+(1._JPRB-SLHDKMAX)*ZWD1
    PWH(JROF,2,JLEV)=(1._JPRB-ZSLHDKMIN)*ZWA2+(1._JPRB-SLHDKMAX)*ZWD2
    PWH(JROF,3,JLEV)=(1._JPRB-ZSLHDKMIN)*ZWA3+(1._JPRB-SLHDKMAX)*ZWD3
    ZWL1=ZSLHDKMIN*ZWA1+SLHDKMAX*ZWD1
    ZWL2=ZSLHDKMIN*ZWA2+SLHDKMAX*ZWD2
    ZWL3=ZSLHDKMIN*ZWA3+SLHDKMAX*ZWD3

    IF (LLSLHDQUAD) THEN
      ZWL3=ZWL3+PSLD1_(IJ_,ILA+1)*ZWL1+PSLD2_(IJ_,ILA+1)*ZWL2
      ZWL1=ZWL1+PSLD3_(IJ_,ILA+1)*PDLAT5(JROF,JLEV)*ZWL3
      ZWL2=ZWL2+PSLD3_(IJ_,ILA+1)*(1._JPRB-PDLAT5(JROF,JLEV))*ZWL3
      PDLAT(JROF,JLEV)=PDLAT(JROF,JLEV)-ZWL1+ZWL2
    ELSEIF (LLSLHD_OLD) THEN
      PDLAT(JROF,JLEV)=PDLAT(JROF,JLEV)-ZWL1+ZWL2
    ENDIF

  ENDDO
  ENDDO

ELSE ! .NOT.LLSLHD

  DO JLEV=1,KFLEV
!CDIR NODEP
  DO JROF=KST,KPROF

    ! a/ Preliminary calculations:
    
    IJ_ = MOD(JROF+1-KST,JPDUP)+1
    ILA=KILA(JROF,JLEV)

    ! c/ AD code:

    IF (LLSLHDQUAD) THEN
      ZWA1=PCLA(JROF,JLEV,1)
      ZWA2=PCLA(JROF,JLEV,2)
      ZWA3=PCLA(JROF,JLEV,3)
      PWH(JROF,1,JLEV)=(1._JPRB-PSLHDKMIN)*ZWA1
      PWH(JROF,2,JLEV)=(1._JPRB-PSLHDKMIN)*ZWA2
      PWH(JROF,3,JLEV)=(1._JPRB-PSLHDKMIN)*ZWA3
      ZWL1=PSLHDKMIN*ZWA1
      ZWL2=PSLHDKMIN*ZWA2
      ZWL3=PSLHDKMIN*ZWA3
      ZWL3=ZWL3+PSLD1_(IJ_,ILA+1)*ZWL1+PSLD2_(IJ_,ILA+1)*ZWL2
      ZWL1=ZWL1+PSLD3_(IJ_,ILA+1)*PDLAT5(JROF,JLEV)*ZWL3
      ZWL2=ZWL2+PSLD3_(IJ_,ILA+1)*(1._JPRB-PDLAT5(JROF,JLEV))*ZWL3
      PDLAT(JROF,JLEV)=PDLAT(JROF,JLEV)-ZWL1+ZWL2
    ELSE
      PWH(JROF,1,JLEV)=PCLA(JROF,JLEV,1)
      PWH(JROF,2,JLEV)=PCLA(JROF,JLEV,2)
      PWH(JROF,3,JLEV)=PCLA(JROF,JLEV,3)
    ENDIF

  ENDDO
  ENDDO

ENDIF ! Test on LLSLHD

!     -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_CLA_AD',1,ZHOOK_HANDLE)
END SUBROUTINE LASCAW_CLA_AD
