SUBROUTINE LASCAW_VINTW_AD(KPROM,KFLEV,KST,KPROF,LDT_SLHD,LDSLHDHEAT,PSLHDKMIN,&
 & KLEV,PLEV,PDVER,PKAPPA,PKAPPAT,PVETA,PVCUICO_,PVSLD_,PVSLDW_,&
 & PLEV5,PDVER5,PKAPPA5,PKAPPAT5,PVINTW,PVINTWSLD,PVINTWSLT)

!     ------------------------------------------------------------------

!**** *LASCAW_VINTW_AD-  Weights for semi-LAgrangian interpolator:
!                        Computes PVINTW and PVINTWSLD/T for one layer
!                        (high-order vertical weights):
!                        AD + trajectory code

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LASCAW_VINTW_AD( ... )

!        Explicit arguments :
!        --------------------

!          KPROM    - horizontal dimension.
!          KFLEV    - vertical dimension.
!          KST      - first element of arrays where computations are performed.
!          KPROF    - depth of work.
!          LDT_SLHD - keys for SLHD.
!          LDSLHDHEAT   - If true, the triggering function for heat variables differs from the one for momentum variables
!          PSLHDKMIN - either HOISLV or SLHDKMIN
!          KLEV     - lower level of the vertical interpolation
!                     grid needed for vertical interpolations.
!          PLEV     - vertical coordinate of the interpolation point.
!          PDVER    - distance for vertical linear interpolation.
!          PKAPPA   - kappa function ("coefficient of SLHD").
!          PKAPPAT  - kappa function ("coefficient of SLHD") on T.
!          PVETA    - Values of ETA.
!          PVCUICO_ - Denominators of the vertical cubic interpolation coef.
!          PVSLD_   - auxiliary quantities for vertical SLHD interpolation.
!          PVSLDW_  - weights for SLHD vertical Laplacian smoother.
!          PLEV5    - cf. PLEV (trajectory).
!          PDVER5   - cf. PDVER (trajectory).
!          PKAPPA5  - cf. PKAPPA (trajectory).
!          PKAPPAT5 - cf. PKAPPAT (trajectory).

!          PVINTW    - vertical cubic interpolation weights.
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

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK 

! arp/ifs dependencies to be solved later.
USE YOMDYNA  , ONLY : YRDYNA

USE EINT_MOD , ONLY : JPDUP

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KPROM
INTEGER(KIND=JPIM), INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM), INTENT(IN)    :: KST
INTEGER(KIND=JPIM), INTENT(IN)    :: KPROF
LOGICAL           , INTENT(IN)    :: LDT_SLHD(3)
LOGICAL           , INTENT(IN)    :: LDSLHDHEAT
REAL(KIND=JPRB)   , INTENT(IN)    :: PSLHDKMIN
INTEGER(KIND=JPIM), INTENT(IN)    :: KLEV(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PLEV(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PDVER(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PKAPPA(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PKAPPAT(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)    :: PVETA(0:KFLEV+1)
REAL(KIND=JPRB)   , INTENT(IN)    :: PVCUICO_(JPDUP,4,0:KFLEV-1)
REAL(KIND=JPRB)   , INTENT(IN)    :: PVSLD_(JPDUP,3,0:KFLEV-1)
REAL(KIND=JPRB)   , INTENT(IN)    :: PVSLDW_(JPDUP,3,3,0:KFLEV-1)
REAL(KIND=JPRB)   , INTENT(IN)    :: PLEV5(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)    :: PDVER5(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)    :: PKAPPA5(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(IN)    :: PKAPPAT5(KPROM,KFLEV)
REAL(KIND=JPRB)   , INTENT(INOUT) :: PVINTW(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(INOUT)    :: PVINTWSLD(KPROM,KFLEV,3)
REAL(KIND=JPRB)   , INTENT(INOUT)    :: PVINTWSLT(KPROM,KFLEV,3)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,IJ_,ILEVV,JLEV
REAL(KIND=JPRB) :: ZWA1,ZWA2,ZWA3,ZWD1,ZWD2,ZWD3,ZWH1,ZWH2,ZWH3
REAL(KIND=JPRB) :: ZWDS1,ZWDS2,ZWDS3,ZWL1,ZWL2,ZWL3
REAL(KIND=JPRB) :: ZWA15,ZWA25,ZWA35,ZWD15,ZWD25,ZWD35,ZWH15,ZWH25,ZWH35
REAL(KIND=JPRB) :: ZWDS15,ZWDS25,ZWDS35,ZWL15,ZWL25,ZWL35
REAL(KIND=JPRB) :: ZD05,ZD15,ZD25,ZD35, ZSIGN, ZSLHDKMIN
LOGICAL :: LLSLHD,LLSLHDQUAD,LLSLHD_OLD
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_VINTW_AD',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

LLSLHD=LDT_SLHD(1)
LLSLHDQUAD=LDT_SLHD(2)
LLSLHD_OLD=LDT_SLHD(3)

! * AD of calculation of PVINTW and PVINTWSLD/T:

IF (LLSLHD) THEN ! LLSLHD

  DO JLEV=1,KFLEV
!CDIR NODEP
  DO JROF=KST,KPROF

    ! a/ Preliminary calculations:
    IJ_ = MOD(JROF+1-KST,JPDUP)+1
    ILEVV=KLEV(JROF,JLEV)

    ! b/ Trajectory code:
    
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
    ZWH15=0._JPRB
    ZWH25=0._JPRB
    ZWH35=0._JPRB
    ZD05=0._JPRB
    ZD15=0._JPRB
    ZD25=0._JPRB
    ZD35=0._JPRB
#endif

    IF (ILEVV >= 1.AND.ILEVV <= KFLEV-3) THEN
      ZD05=PLEV5(JROF,JLEV)-PVETA(ILEVV  )
      ZD15=PLEV5(JROF,JLEV)-PVETA(ILEVV+1)
      ZD25=PLEV5(JROF,JLEV)-PVETA(ILEVV+2)
      ZD35=PLEV5(JROF,JLEV)-PVETA(ILEVV+3)
      ZWH15=ZD05     *ZD25*ZD35*PVCUICO_(IJ_,2,ILEVV)
      ZWH25=ZD05*ZD15     *ZD35*PVCUICO_(IJ_,3,ILEVV)
      ZWH35=ZD05*ZD15*ZD25     *PVCUICO_(IJ_,4,ILEVV)
      IF (LLSLHDQUAD) THEN
        ZWL25=PDVER5(JROF,JLEV)
        ZWL15=1.0_JPRB-ZWL25
        ZWL35=PVSLD_(IJ_,3,ILEVV)*ZWL15*ZWL25
        ZWL15=ZWL15+PVSLD_(IJ_,1,ILEVV)*ZWL35
        ZWL25=ZWL25+PVSLD_(IJ_,2,ILEVV)*ZWL35
      ELSEIF (LLSLHD_OLD) THEN
        ZWL25=PDVER5(JROF,JLEV)
        ZWL15=1.0_JPRB-ZWL25
        ZWL35=0.0_JPRB
      ENDIF

      ZSIGN=SIGN(0.5_JPRB,PKAPPA5(JROF,JLEV))
      ZSLHDKMIN=(0.5_JPRB+ZSIGN)*PSLHDKMIN - (ZSIGN-0.5_JPRB)*YRDYNA%SLHDKREF
      ZWA15=ZWH15+ZSLHDKMIN*(ZWL15-ZWH15)
      ZWA25=ZWH25+ZSLHDKMIN*(ZWL25-ZWH25)
      ZWA35=ZWH35+ZSLHDKMIN*(ZWL35-ZWH35)
      ZWD15=ZWH15+YRDYNA%SLHDKMAX*(ZWL15-ZWH15)
      ZWD25=ZWH25+YRDYNA%SLHDKMAX*(ZWL25-ZWH25)
      ZWD35=ZWH35+YRDYNA%SLHDKMAX*(ZWL35-ZWH35)
      ZWDS15=PVSLDW_(IJ_,1,1,ILEVV)*ZWD15+PVSLDW_(IJ_,1,2,ILEVV)*ZWD25+&
       & PVSLDW_(IJ_,1,3,ILEVV)*ZWD35
      ZWDS25=PVSLDW_(IJ_,2,1,ILEVV)*ZWD15+PVSLDW_(IJ_,2,2,ILEVV)*ZWD25+&
       & PVSLDW_(IJ_,2,3,ILEVV)*ZWD35
      ZWDS35=PVSLDW_(IJ_,3,1,ILEVV)*ZWD15+PVSLDW_(IJ_,3,2,ILEVV)*ZWD25+&
       & PVSLDW_(IJ_,3,3,ILEVV)*ZWD35
    ENDIF

    ! c/ AD code:

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
    ZWH1=0._JPRB
    ZWH2=0._JPRB
    ZWH3=0._JPRB
#endif

    IF (ILEVV >= 1.AND.ILEVV <= KFLEV-3) THEN
      PKAPPA(JROF,JLEV)=PKAPPA(JROF,JLEV)&
       & +(ZWDS15-ZWA15)*PVINTWSLD(JROF,JLEV,1)&
       & +(ZWDS25-ZWA25)*PVINTWSLD(JROF,JLEV,2)&
       & +(ZWDS35-ZWA35)*PVINTWSLD(JROF,JLEV,3)
      ZWDS1=ABS(PKAPPA5(JROF,JLEV))*PVINTWSLD(JROF,JLEV,1)
      ZWDS2=ABS(PKAPPA5(JROF,JLEV))*PVINTWSLD(JROF,JLEV,2)
      ZWDS3=ABS(PKAPPA5(JROF,JLEV))*PVINTWSLD(JROF,JLEV,3)
      ZWA1=(1._JPRB-ABS(PKAPPA5(JROF,JLEV)))*PVINTWSLD(JROF,JLEV,1)+PVINTW(JROF,JLEV,1)
      ZWA2=(1._JPRB-ABS(PKAPPA5(JROF,JLEV)))*PVINTWSLD(JROF,JLEV,2)+PVINTW(JROF,JLEV,2)
      ZWA3=(1._JPRB-ABS(PKAPPA5(JROF,JLEV)))*PVINTWSLD(JROF,JLEV,3)+PVINTW(JROF,JLEV,3)
      PVINTWSLD(JROF,JLEV,1:3)=0._JPRB
      PVINTW(JROF,JLEV,1:3)=0._JPRB
      IF (LDSLHDHEAT) THEN
        PKAPPAT(JROF,JLEV)=PKAPPAT(JROF,JLEV)&
         & +(ZWDS15-ZWA15)*PVINTWSLT(JROF,JLEV,1)&
         & +(ZWDS25-ZWA25)*PVINTWSLT(JROF,JLEV,2)&
         & +(ZWDS35-ZWA35)*PVINTWSLT(JROF,JLEV,3)
        ZWDS1=ZWDS1+ABS(PKAPPAT5(JROF,JLEV))*PVINTWSLT(JROF,JLEV,1)
        ZWDS2=ZWDS2+ABS(PKAPPAT5(JROF,JLEV))*PVINTWSLT(JROF,JLEV,2)
        ZWDS3=ZWDS3+ABS(PKAPPAT5(JROF,JLEV))*PVINTWSLT(JROF,JLEV,3)
        ZWA1=ZWA1 + (1._JPRB-ABS(PKAPPAT5(JROF,JLEV)))*PVINTWSLT(JROF,JLEV,1)
        ZWA2=ZWA2 + (1._JPRB-ABS(PKAPPAT5(JROF,JLEV)))*PVINTWSLT(JROF,JLEV,2)
        ZWA3=ZWA3 + (1._JPRB-ABS(PKAPPAT5(JROF,JLEV)))*PVINTWSLT(JROF,JLEV,3)
        PVINTWSLT(JROF,JLEV,1:3)=0._JPRB
      ENDIF
      ZWD1=PVSLDW_(IJ_,1,1,ILEVV)*ZWDS1+PVSLDW_(IJ_,2,1,ILEVV)*ZWDS2&
       & +PVSLDW_(IJ_,3,1,ILEVV)*ZWDS3
      ZWD2=PVSLDW_(IJ_,1,2,ILEVV)*ZWDS1+PVSLDW_(IJ_,2,2,ILEVV)*ZWDS2&
       & +PVSLDW_(IJ_,3,2,ILEVV)*ZWDS3
      ZWD3=PVSLDW_(IJ_,1,3,ILEVV)*ZWDS1+PVSLDW_(IJ_,2,3,ILEVV)*ZWDS2&
       & +PVSLDW_(IJ_,3,3,ILEVV)*ZWDS3
      ZWH1=(1._JPRB-ZSLHDKMIN)*ZWA1+(1._JPRB-YRDYNA%SLHDKMAX)*ZWD1
      ZWH2=(1._JPRB-ZSLHDKMIN)*ZWA2+(1._JPRB-YRDYNA%SLHDKMAX)*ZWD2
      ZWH3=(1._JPRB-ZSLHDKMIN)*ZWA3+(1._JPRB-YRDYNA%SLHDKMAX)*ZWD3
      ZWL1=ZSLHDKMIN*ZWA1+YRDYNA%SLHDKMAX*ZWD1
      ZWL2=ZSLHDKMIN*ZWA2+YRDYNA%SLHDKMAX*ZWD2
      ZWL3=ZSLHDKMIN*ZWA3+YRDYNA%SLHDKMAX*ZWD3

      IF (LLSLHDQUAD) THEN
        ZWL3=ZWL3+PVSLD_(IJ_,1,ILEVV)*ZWL1+PVSLD_(IJ_,2,ILEVV)*ZWL2
        ZWL1=ZWL1+PVSLD_(IJ_,3,ILEVV)*PDVER5(JROF,JLEV)*ZWL3
        ZWL2=ZWL2+PVSLD_(IJ_,3,ILEVV)*(1._JPRB-PDVER5(JROF,JLEV))*ZWL3
        PDVER(JROF,JLEV)=PDVER(JROF,JLEV)-ZWL1+ZWL2
      ELSEIF (LLSLHD_OLD) THEN
        PDVER(JROF,JLEV)=PDVER(JROF,JLEV)-ZWL1+ZWL2
      ENDIF
      PLEV(JROF,JLEV)=PLEV(JROF,JLEV)&
       & +(ZD05*(ZD25+ZD35)+ZD25*ZD35)*ZWH1*PVCUICO_(IJ_,2,ILEVV)&
       & +(ZD05*(ZD15+ZD35)+ZD15*ZD35)*ZWH2*PVCUICO_(IJ_,3,ILEVV)&
       & +(ZD05*(ZD15+ZD25)+ZD15*ZD25)*ZWH3*PVCUICO_(IJ_,4,ILEVV)
    ELSE
      PVINTW(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)+PVINTWSLD(JROF,JLEV,1)
      PVINTW(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)+PVINTWSLD(JROF,JLEV,2)
      PVINTW(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)+PVINTWSLD(JROF,JLEV,3)
      PVINTWSLD(JROF,JLEV,1:3)=0._JPRB
      IF (LDSLHDHEAT) THEN
        PVINTW(JROF,JLEV,1)=PVINTW(JROF,JLEV,1)+PVINTWSLT(JROF,JLEV,1)
        PVINTW(JROF,JLEV,2)=PVINTW(JROF,JLEV,2)+PVINTWSLT(JROF,JLEV,2)
        PVINTW(JROF,JLEV,3)=PVINTW(JROF,JLEV,3)+PVINTWSLT(JROF,JLEV,3)
        PVINTWSLT(JROF,JLEV,1:3)=0._JPRB
      ENDIF
      PDVER(JROF,JLEV)=PDVER(JROF,JLEV)-PVINTW(JROF,JLEV,1)+PVINTW(JROF,JLEV,2)
      PVINTW(JROF,JLEV,1:3)=0._JPRB
    ENDIF

  ENDDO
  ENDDO

ELSE ! .NOT.LLSLHD

  DO JLEV=1,KFLEV
!CDIR NODEP
  DO JROF=KST,KPROF

    ! a/ Preliminary calculations:
    IJ_ = MOD(JROF+1-KST,JPDUP)+1
    ILEVV=KLEV(JROF,JLEV)

    ! b/ Trajectory code:

#if defined(NECSX)
    ! scalar variable initialization (to meet vectorization condition)
    ZD05=0._JPRB
    ZD15=0._JPRB
    ZD25=0._JPRB
    ZD35=0._JPRB
#endif

    IF (ILEVV >= 1.AND.ILEVV <= KFLEV-3) THEN
      ZD05=PLEV5(JROF,JLEV)-PVETA(ILEVV  )
      ZD15=PLEV5(JROF,JLEV)-PVETA(ILEVV+1)
      ZD25=PLEV5(JROF,JLEV)-PVETA(ILEVV+2)
      ZD35=PLEV5(JROF,JLEV)-PVETA(ILEVV+3)
    ENDIF

    ! c/ AD code:

#if defined(NECSX)
    ! scalar variable initialization (to meet vectorization condition)
    ZWA1=0._JPRB
    ZWA2=0._JPRB
    ZWA3=0._JPRB
    ZWL1=0._JPRB
    ZWL2=0._JPRB
    ZWL3=0._JPRB
    ZWH1=0._JPRB
    ZWH2=0._JPRB
    ZWH3=0._JPRB
#endif

    IF (ILEVV >= 1.AND.ILEVV <= KFLEV-3) THEN
      IF (LLSLHDQUAD) THEN
        ZWA1=PVINTW(JROF,JLEV,1)
        ZWA2=PVINTW(JROF,JLEV,2)
        ZWA3=PVINTW(JROF,JLEV,3)
        ZWH1=(1._JPRB-PSLHDKMIN)*ZWA1
        ZWH2=(1._JPRB-PSLHDKMIN)*ZWA2
        ZWH3=(1._JPRB-PSLHDKMIN)*ZWA3
        ZWL1=PSLHDKMIN*ZWA1
        ZWL2=PSLHDKMIN*ZWA2
        ZWL3=PSLHDKMIN*ZWA3
        ZWL3=ZWL3+PVSLD_(IJ_,1,ILEVV)*ZWL1+PVSLD_(IJ_,2,ILEVV)*ZWL2
        ZWL1=ZWL1+PVSLD_(IJ_,3,ILEVV)*PDVER5(JROF,JLEV)*ZWL3
        ZWL2=ZWL2+PVSLD_(IJ_,3,ILEVV)*(1._JPRB-PDVER5(JROF,JLEV))*ZWL3
        PDVER(JROF,JLEV)=PDVER(JROF,JLEV)-ZWL1+ZWL2
      ELSE
        ZWH1=PVINTW(JROF,JLEV,1)
        ZWH2=PVINTW(JROF,JLEV,2)
        ZWH3=PVINTW(JROF,JLEV,3)
      ENDIF
      PVINTW(JROF,JLEV,1:3)=0._JPRB

      PLEV(JROF,JLEV)=PLEV(JROF,JLEV)&
       & +(ZD05*(ZD25+ZD35)+ZD25*ZD35)*ZWH1*PVCUICO_(IJ_,2,ILEVV)&
       & +(ZD05*(ZD15+ZD35)+ZD15*ZD35)*ZWH2*PVCUICO_(IJ_,3,ILEVV)&
       & +(ZD05*(ZD15+ZD25)+ZD15*ZD25)*ZWH3*PVCUICO_(IJ_,4,ILEVV)
    ELSE
      PDVER(JROF,JLEV)=PDVER(JROF,JLEV)-PVINTW(JROF,JLEV,1)+PVINTW(JROF,JLEV,2)
      PVINTW(JROF,JLEV,1:3)=0._JPRB
    ENDIF

  ENDDO
  ENDDO
ENDIF ! Test on LLSLHD

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LASCAW_VINTW_AD',1,ZHOOK_HANDLE)
END SUBROUTINE LASCAW_VINTW_AD
