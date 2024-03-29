!OPTIONS XOPT(NOEVAL)
SUBROUTINE ACNPART(YDML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
 & PAPHI,PAPHIF,PAPRSF,PDECRD,PNEB,&
 & PCLCH,PCLCM,PCLCL,PCLCT,PCLCT_RAD,&
 ! optional arguments (convective cloud cover)
 & PCLCC,PNEBC,PTOPC)

! Purpose:
! --------
!   ACNPART - computes high/medium/low, convective and total cloud cover.
!   Several overlap options are implemented.

! Interface:
! ----------
! INPUT:
!   KIDIA     - initial index for horizontal loops
!   KFDIA     - final index for horizontal loops
!   KLON      - horizontal dimension of arrays
!   KTDIA     - initial index for vertical loops (usually 1)
!   KLEV      - vertical dimension of full level arrays
!   PAPHI     - half level geopotential
!   PAPHIF    - full level geopotential
!   PAPRSF    - full level pressure
!   PDECRD    - decorrelation depth for cloud overlaps [Pa]
!   PNEB      - total cloud cover on levels (protected from 0 and 1)

! OUTPUT:
!   PCLCH     - high cloud cover
!   PCLCM     - medium cloud cover
!   PCLCL     - low cloud cover
!   PCLCT     - total cloud cover
!   PCLCT_RAD - total cloud cover for radiation

! INPUT, OPTIONAL:
!   PNEBC     - convective cloud cover on levels (protected from 0 and 1,
!               missing in AROME)

! OUTPUT, OPTIONAL:
!   PCLCC     - convective cloud cover (missing in AROME)
!   PTOPC     - TOP of convective cloud  [Pa] (missing in AROME)


! Externals:
! ----------

! Method:
! -------

! Reference:
! ----------

! Author:
! -------
!   2007-02, R. Brozkova

! Modifications:
! --------------
!   2009-03, C. Wittmann
!   Introduction of LACPANMX and WMXOV.
!
!   2009-07, K. Yessad
!   Remove CDLOCK + some cleaning.
!
!   2009-10, L. Bengtsson
!   Introduction of LWMOCLOUD.
!
!   2016-04, J. Masek
!   Introduction of LRNUEXP (exponential-random overlap), fix for LWMOCLOUD,
!   modularization, reordering of arguments, PCLCC optional (missing in AROME).
!
!   2016-09, J. Masek
!   Introduction of radiative cloud cover PCLCT_RAD, needed for consistent
!   calculation of sunshine duration.
!
!   2018-09, J. Masek
!   Fix of convective cloud cover when WMXOV or RDECRDRED differ from 1.
!
!   2018-07, O. Jaron
!   Introduction of CONV_BASE_TOP to diagnose base and top of convective
!   clouds. (Base desactivated)
!
! End Modifications
!-------------------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK


USE YOMCST   , ONLY : RG

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHI(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHIF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDECRD(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNEB(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLCH(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLCM(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLCL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLCT(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLCT_RAD(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT),OPTIONAL :: PCLCC(KLON)
REAL(KIND=JPRB)   ,INTENT(IN), OPTIONAL :: PNEBC(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT),OPTIONAL :: PTOPC(KLON)

!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZDECRDRED,ZWMXOV
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('ACNPART',0,ZHOOK_HANDLE)
ASSOCIATE(LACPANMX =>YDML_PHY_MF%YRPHY%LACPANMX,   &
 &        LRNUMX   =>YDML_PHY_MF%YRPHY%LRNUMX,     &
 &        LRNUEXP  =>YDML_PHY_MF%YRPHY%LRNUEXP,    &
 &        RDECRDRED=>YDML_PHY_MF%YRPHY0%RDECRDRED, &
 &        WMXOV    =>YDML_PHY_MF%YRPHY0%WMXOV,     &
 &        LWMOCLOUD=>YDML_PHY_MF%YRPHY2%LWMOCLOUD, &
 &        NTSML    =>YDML_PHY_MF%YRPHY2%NTSML,     &
 &        NTSHM    =>YDML_PHY_MF%YRPHY2%NTSHM,     &
 &        LPTOPC   =>YDML_PHY_MF%YRPHY%LPTOPC)
! settings for diagnostic cloud cover
ZDECRDRED=RDECRDRED
ZWMXOV   =WMXOV

IF ( LWMOCLOUD ) THEN

  ! high/medium/low cloud cover according to WMO heights
  CALL CLOUD_COVER_WMO(PNEB,PAPHI,PAPRSF,PCLCH,PCLCM,PCLCL)

ELSE

  ! high/medium/low cloud cover according to fixed model levels
  CALL CLOUD_COVER(KTDIA  ,NTSHM,PNEB,PAPRSF,PCLCH)  ! high
  CALL CLOUD_COVER(NTSHM+1,NTSML,PNEB,PAPRSF,PCLCM)  ! medium
  CALL CLOUD_COVER(NTSML+1,KLEV ,PNEB,PAPRSF,PCLCL)  ! low

ENDIF

! total cloud cover
CALL CLOUD_COVER(KTDIA,KLEV,PNEB,PAPRSF,PCLCT)

! convective cloud cover
IF ( PRESENT(PCLCC) ) THEN
  CALL CLOUD_COVER(KTDIA,KLEV,PNEBC,PAPRSF,PCLCC)
ENDIF

! total cloud cover for radiation
IF ( LRNUMX.AND.(LACPANMX.OR.LRNUEXP) ) THEN
  ZDECRDRED=1._JPRB  ! do not reduce decorrelation depth
  ZWMXOV   =1._JPRB  ! ignore LACPANMX
  CALL CLOUD_COVER(KTDIA,KLEV,PNEB ,PAPRSF,PCLCT_RAD)
ELSE
  PCLCT_RAD(:)=PCLCT(:)
ENDIF

! convective top
IF (PRESENT(PTOPC).AND. LPTOPC) THEN
  CALL CONV_BASE_TOP(PNEBC,PAPRSF,PAPHIF,PTOPC)
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACNPART',1,ZHOOK_HANDLE)

! ------------------------------------------------------------------------------

CONTAINS

! Private subroutines/functions
! -----------------------------

! cloud cover between model levels KL1 and KL2
SUBROUTINE CLOUD_COVER(KL1,KL2,PNEB,PAPRSF,PCLC)

! Interface:
! ----------
! INPUT:
!   KL1    - initial level
!   KL2    - final level
!   PNEB   - cloud cover on levels
!   PAPRSF - full level pressure

! OUTPUT:
!   PCLC   - cloud cover between model levels KL1 and KL2

INTEGER(KIND=JPIM),INTENT(IN) :: KL1
INTEGER(KIND=JPIM),INTENT(IN) :: KL2

REAL(KIND=JPRB),INTENT(IN)  :: PNEB  (KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)  :: PAPRSF(KLON,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PCLC  (KLON)

INTEGER(KIND=JPIM) :: JLON,JLEV

REAL(KIND=JPRB) :: ZCLOV
REAL(KIND=JPRB) :: ZHOOK_HANDLE

REAL(KIND=JPRB) :: ZOVLP(KLON)

IF (LHOOK) CALL DR_HOOK('ACNPART:CLOUD_COVER',0,ZHOOK_HANDLE)
ASSOCIATE(LACPANMX=>YDML_PHY_MF%YRPHY%LACPANMX, &
 &        LRNUMX  =>YDML_PHY_MF%YRPHY%LRNUMX,   &
 &        LRNUEXP =>YDML_PHY_MF%YRPHY%LRNUEXP)

! initialize cloud cover or its complement
IF (LRNUMX.AND.LACPANMX) THEN
  DO JLON=KIDIA,KFDIA
    PCLC(JLON)=1._JPRB-PNEB(JLON,KL1)
  ENDDO
ELSE
  DO JLON=KIDIA,KFDIA
    PCLC(JLON)=PNEB(JLON,KL1)
  ENDDO
ENDIF
 
! compute cloud cover using actual overlap hypothesis
DO JLEV=KL1+1,KL2

  IF (LRNUMX.AND.LACPANMX) THEN

    ! nearly maximum-random overlap
    DO JLON=KIDIA,KFDIA
      ZOVLP(JLON)=(MIN(1.0_JPRB-PNEB(JLON,JLEV),1.0_JPRB-ZWMXOV*&
       & PNEB(JLON,JLEV-1)))/(1.0_JPRB-ZWMXOV*PNEB(JLON,JLEV-1))
    ENDDO

  ELSEIF (LRNUMX.AND.LRNUEXP) THEN

    ! exponential-random overlap
    DO JLON=KIDIA,KFDIA
      ZCLOV      =EXP((PAPRSF(JLON,JLEV-1)-PAPRSF(JLON,JLEV))/&
       & (ZDECRDRED*PDECRD(JLON)))
      ZOVLP(JLON)=(1._JPRB-ZCLOV)*PNEB(JLON,JLEV)*PNEB(JLON,JLEV-1)+&
       & ZCLOV*MIN(PNEB(JLON,JLEV),PNEB(JLON,JLEV-1))
    ENDDO

  ELSEIF (LRNUMX) THEN

    ! maximum-random overlap
    DO JLON=KIDIA,KFDIA
      ZOVLP(JLON)=MIN(PNEB(JLON,JLEV),PNEB(JLON,JLEV-1))
    ENDDO

  ELSE

    ! random overlap
    DO JLON=KIDIA,KFDIA
      ZOVLP(JLON)=PNEB(JLON,JLEV)*PNEB(JLON,JLEV-1)
    ENDDO

  ENDIF

  ! update cloud cover
  IF (LRNUMX.AND.LACPANMX) THEN
    DO JLON=KIDIA,KFDIA
      PCLC(JLON)=PCLC(JLON)*ZOVLP(JLON)
    ENDDO
  ELSE
    DO JLON=KIDIA,KFDIA
      PCLC(JLON)=PCLC(JLON)+(PNEB(JLON,JLEV)-ZOVLP(JLON))*&
       & (1._JPRB-PCLC(JLON))/(1._JPRB-PNEB(JLON,JLEV-1))
    ENDDO
  ENDIF

ENDDO

! convert complement of cloud cover to cloud cover
IF (LRNUMX.AND.LACPANMX) THEN
  DO JLON=KIDIA,KFDIA
    PCLC(JLON)=1._JPRB-PCLC(JLON)
  ENDDO
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACNPART:CLOUD_COVER',1,ZHOOK_HANDLE)

END SUBROUTINE CLOUD_COVER

! -----

! high/medium/low cloud cover according to WMO heights
SUBROUTINE CLOUD_COVER_WMO(PNEB,PAPHI,PAPRSF,PCLCH,PCLCM,PCLCL)

! Interface:
! ----------
! INPUT:
!   PNEB   - cloud cover on levels
!   PAPHI  - half level geopotential
!   PAPRSF - full level pressure

! OUTPUT:
!   PCLCH  - high cloud cover
!   PCLCM  - medium cloud cover
!   PCLCL  - low cloud cover

REAL(KIND=JPRB),INTENT(IN)  :: PNEB  (KLON,  KLEV)
REAL(KIND=JPRB),INTENT(IN)  :: PAPHI (KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(IN)  :: PAPRSF(KLON,  KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PCLCH (KLON)
REAL(KIND=JPRB),INTENT(OUT) :: PCLCM (KLON)
REAL(KIND=JPRB),INTENT(OUT) :: PCLCL (KLON)

INTEGER(KIND=JPIM) :: JLON,JLEV

REAL(KIND=JPRB) :: ZCLOV,ZH,ZL
REAL(KIND=JPRB) :: ZHOOK_HANDLE

REAL(KIND=JPRB) :: ZOVLP(KLON)

IF (LHOOK) CALL DR_HOOK('ACNPART:CLOUD_COVER_WMO',0,ZHOOK_HANDLE)
ASSOCIATE(HWMOLOW =>YDML_PHY_MF%YRPHY2%HWMOLOW,  &
 &        HWMOHIGH=>YDML_PHY_MF%YRPHY2%HWMOHIGH, &
 &        LACPANMX=>YDML_PHY_MF%YRPHY%LACPANMX,  &
 &        LRNUMX  =>YDML_PHY_MF%YRPHY%LRNUMX,    &
 &        LRNUEXP =>YDML_PHY_MF%YRPHY%LRNUEXP)

! initialize high cloud cover or its complement
IF (LRNUMX.AND.LACPANMX) THEN
  DO JLON=KIDIA,KFDIA
    PCLCH(JLON)=1._JPRB-PNEB(JLON,KTDIA)
  ENDDO
ELSE
  DO JLON=KIDIA,KFDIA
    PCLCH(JLON)=PNEB(JLON,KTDIA)
  ENDDO
ENDIF

! compute cloud cover using actual overlap hypothesis
DO JLEV=KTDIA+1,KLEV
  IF (LRNUMX.AND.LACPANMX) THEN

    ! nearly maximum-random overlap
    DO JLON=KIDIA,KFDIA
      ZOVLP(JLON)=(MIN(1.0_JPRB-PNEB(JLON,JLEV),1.0_JPRB-ZWMXOV*&
       & PNEB(JLON,JLEV-1)))/(1.0_JPRB-ZWMXOV*PNEB(JLON,JLEV-1))
    ENDDO

  ELSEIF (LRNUMX.AND.LRNUEXP) THEN

    ! exponential-random overlap
    DO JLON=KIDIA,KFDIA
      ZCLOV      =EXP((PAPRSF(JLON,JLEV-1)-PAPRSF(JLON,JLEV))/&
       & (ZDECRDRED*PDECRD(JLON)))
      ZOVLP(JLON)=(1._JPRB-ZCLOV)*PNEB(JLON,JLEV)*PNEB(JLON,JLEV-1)+&
       & ZCLOV*MIN(PNEB(JLON,JLEV),PNEB(JLON,JLEV-1))
    ENDDO

  ELSEIF (LRNUMX) THEN

    ! maximum-random overlap
    DO JLON=KIDIA,KFDIA
      ZOVLP(JLON)=MIN(PNEB(JLON,JLEV),PNEB(JLON,JLEV-1))
    ENDDO

  ELSE

    ! random overlap
    DO JLON=KIDIA,KFDIA
      ZOVLP(JLON)=PNEB(JLON,JLEV)*PNEB(JLON,JLEV-1)
    ENDDO

  ENDIF

  ! update cloud cover
  IF (LRNUMX.AND.LACPANMX) THEN
    DO JLON=KIDIA,KFDIA
      ZH=(PAPHI(JLON,JLEV-1)-PAPHI(JLON,KLEV))/RG
      ZL=(PAPHI(JLON,JLEV  )-PAPHI(JLON,KLEV))/RG
      IF (ZL >= HWMOHIGH) THEN
        PCLCH(JLON)=PCLCH(JLON)*ZOVLP(JLON)
      ELSEIF (ZH >= HWMOHIGH) THEN
        PCLCM(JLON)=1._JPRB-PNEB(JLON,JLEV)  ! complement of high cloud cover
      ELSEIF (ZL >= HWMOLOW) THEN
        PCLCM(JLON)=PCLCM(JLON)*ZOVLP(JLON)
      ELSEIF (ZH >= HWMOLOW) THEN
        PCLCL(JLON)=1._JPRB-PNEB(JLON,JLEV)  ! complement of low cloud cover
      ELSE
        PCLCL(JLON)=PCLCL(JLON)*ZOVLP(JLON)
      ENDIF
    ENDDO
  ELSE
    DO JLON=KIDIA,KFDIA
      ZH=(PAPHI(JLON,JLEV-1)-PAPHI(JLON,KLEV))/RG
      ZL=(PAPHI(JLON,JLEV  )-PAPHI(JLON,KLEV))/RG
      IF (ZL >= HWMOHIGH) THEN
        PCLCH(JLON)=PCLCH(JLON)+(PNEB(JLON,JLEV)-ZOVLP(JLON))*&
         & (1._JPRB-PCLCH(JLON))/(1._JPRB-PNEB(JLON,JLEV-1))
      ELSEIF (ZH >= HWMOHIGH) THEN
        PCLCM(JLON)=PNEB(JLON,JLEV)  ! initialize medium cloud cover
      ELSEIF (ZL >= HWMOLOW) THEN
        PCLCM(JLON)=PCLCM(JLON)+(PNEB(JLON,JLEV)-ZOVLP(JLON))*&
         & (1._JPRB-PCLCM(JLON))/(1._JPRB-PNEB(JLON,JLEV-1))
      ELSEIF (ZH >= HWMOLOW) THEN
        PCLCL(JLON)=PNEB(JLON,JLEV)  ! initialize low cloud cover
      ELSE
        PCLCL(JLON)=PCLCL(JLON)+(PNEB(JLON,JLEV)-ZOVLP(JLON))*&
         & (1._JPRB-PCLCL(JLON))/(1._JPRB-PNEB(JLON,JLEV-1))
      ENDIF
    ENDDO
  ENDIF

ENDDO

! convert complement of cloud cover to cloud cover
IF (LRNUMX.AND.LACPANMX) THEN
  DO JLON=KIDIA,KFDIA
    PCLCH(JLON)=1._JPRB-PCLCH(JLON)
    PCLCM(JLON)=1._JPRB-PCLCM(JLON)
    PCLCL(JLON)=1._JPRB-PCLCL(JLON)
  ENDDO
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACNPART:CLOUD_COVER_WMO',1,ZHOOK_HANDLE)

END SUBROUTINE CLOUD_COVER_WMO

SUBROUTINE CONV_BASE_TOP(PNEBC,PAPRSF,PAPHIF,PTOPC)
 ! TOP  of convective cloud
 ! If BASE is needed, add PBASEC in the arguments of this subroutine.

! Interface:
! ----------
! INPUT:
!   PNEBC     - convective cloud cover on levels (protected from 0 and 1,
!               missing in AROME)
!   PAPRSF    - full level pressure
! OUTPUT :
!!   PBASEC    - BASE of convective cloud [Pa] 
!   PTOPC     - TOP of convective cloud  [Pa]

REAL(KIND=JPRB)   ,INTENT(IN)    :: PNEBC(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHIF(KLON,KLEV)
!REAL(KIND=JPRB)  ,INTENT(OUT)   :: PBASEC(KLON) ! Base desactivated
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTOPC(KLON)

REAL(KIND=JPRB) :: ZNCCRIT, ZALPHA(KLON),ZBETA(KLON),ZABINV(KLON)
INTEGER(KIND=JPIM) :: JLON,JLEV
INTEGER(KIND=JPIM) :: JCTOP(KLON), JCBASE(KLON)

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ACNPART:CONV_BASE_TOP',0,ZHOOK_HANDLE)

ASSOCIATE(GTOPDEPTH=>YDML_PHY_MF%YRPHY0%GTOPDEPTH)

ZNCCRIT=0.1_JPRB
JCTOP(KIDIA:KFDIA) = KLEV
JCBASE(KIDIA:KFDIA)= KLEV
! Recherche de la base et du sommet
DO JLEV=KLEV-1,1,-1
    DO JLON=KIDIA,KFDIA
      !On renverse la coordonnees verticale (KLEV-JLEV), 
      JCTOP(JLON)=KLEV-MAX((KLEV-JLEV)*INT(SIGN(1._JPRB,PNEBC(JLON,JLEV) &
            & -ZNCCRIT)),(KLEV-JCTOP(JLON)))
      JCBASE(JLON)=KLEV-MAX((KLEV-JLEV)*SIGN(1_JPIM,JCTOP(JLON)-KLEV), &
                     & KLEV-JCBASE(JLON))
    ENDDO
ENDDO
! Critere epaisseur 
!If NEBUL < 10%, JCBASE=1 et JCTOP=KLEV
DO JLON=KIDIA,KFDIA
  JCTOP(JLON)= KLEV-MAX((KLEV-JCTOP(JLON)) &
   & *INT(SIGN(1._JPRB,((PAPHIF(JLON,JCTOP(JLON))- & 
   & PAPHIF(JLON,JCBASE(JLON)))/RG)-GTOPDEPTH)),0_JPIM)
ENDDO

DO JLON=KIDIA,KFDIA
  !Mise à KLEV si pas de base convective
  !JCBASE(JLON)=MAX(SIGN(1_JPIM,JCBASE(JLON)-2_JPIM),0_JPIM)*JCBASE(JLON)&
  !          & +MAX(SIGN(1_JPIM,1_JPIM-JCBASE(JLON)),0_JPIM)*KLEV

  ! Interpolation lineaire sur Z de la base
  !ZALPHA(JLON)=ABS(PNEBC(JLON,JCBASE(JLON)-1_JPIM)-ZNCCRIT) !danger
  !ZBETA(JLON)= ABS(ZNCCRIT-PNEBC(JLON,JCBASE(JLON)))
  !ZABINV(JLON)=1._JPRB/(ZALPHA(JLON)+ZBETA(JLON))
  !PBASEC(JLON)=PAPRSF(JLON,JCBASE(JLON))*ZALPHA(JLON)*ZABINV(JLON) & 
  !     & + PAPRSF(JLON,JCBASE(JLON)-1_JPIM)*ZBETA(JLON)*ZABINV(JLON)

  ! Interpolation lineaire sur P du sommet 
  ZALPHA(JLON)=ABS(PNEBC(JLON,JCTOP(JLON))-ZNCCRIT)
  ZBETA(JLON)= ABS(ZNCCRIT-PNEBC(JLON,JCTOP(JLON)-1_JPIM))
  ZABINV(JLON)=1._JPRB/(ZALPHA(JLON)+ZBETA(JLON))
  PTOPC(JLON) = PAPRSF(JLON,JCTOP(JLON))*ZBETA(JLON)*ZABINV(JLON) &
           & +PAPRSF(JLON,JCTOP(JLON)-1_JPIM)*ZALPHA(JLON)*ZABINV(JLON)

  PTOPC(JLON) =  PAPRSF(JLON,KLEV) * FLOAT(MAX(0_JPIM, SIGN(1_JPIM, JCTOP(JLON)-KLEV)))&
    &         +  PTOPC(JLON)       * FLOAT(MAX(0_JPIM, SIGN(1_JPIM, KLEV-JCTOP(JLON)-1_JPIM)))
ENDDO



END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACNPART:CONV_BASE_TOP',1,ZHOOK_HANDLE)
END SUBROUTINE CONV_BASE_TOP

END SUBROUTINE ACNPART
