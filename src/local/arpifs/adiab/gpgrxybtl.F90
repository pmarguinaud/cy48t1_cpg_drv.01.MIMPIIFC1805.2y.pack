SUBROUTINE GPGRXYBTL(&
 ! --- INPUT --------------------------------------
 & KPROMA,KD,KF,KFLEV,YDVAB,PREL,PREM,PXYB,&
 ! --- OUTPUT -------------------------------------
 & PXYBDER,&
 ! --- INPUT (TRAJECTORY) -------------------------
 & PRE5L,PRE5M,PXYB5,PXYBDER5)  

!**** *GPGRXYBTL* - TL of routine "GPGRXYB".
!                 Computation of the horizontal gradient of quantities
!                 "alpha" and "delta" at model levels.

!     Purpose.
!     --------

!      See GPGRXYB

!**   Interface.
!     ----------
!        *CALL* *GPGRXYBTL(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!           KPROMA       : horizontal dimension
!           KD           : start of work
!           KF           : working length
!           KFLEV        : number of levels
!           YDVAB        : contains information about hybrid vertical coordinate
!           PREL         : zonal component of "grad prehyds"
!           PREM         : meridian component of "grad prehyds"
!           PXYB         : contains pressure depth, "delta", "alpha".

!         * OUTPUT:
!           PXYBDER      : contains grad(delta), grad(alpha), grad(alpha + log prehyd)

!         * INPUT-TRAJECTORY:
!           PRE5L,PRE5M,PXYB5,PXYBDER5

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
!        K. YESSAD

!     Modifications.
!     --------------
!        Original : March 2006
!        K. Yessad (Dec 2008): remove dummy CDLOCK
!        K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!        K. Yessad (Dec 2011): use YDVAB.
!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK
USE YOMDYNA   , ONLY : NDLNPR
USE YOMCVER   , ONLY : LVERTFE 
USE YOMVERT   , ONLY : TVAB
USE INTDYN_MOD, ONLY : YYTXYB, YYTXYBT, YYTXYBDER, YYTXYBDERT

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KD 
INTEGER(KIND=JPIM),INTENT(IN)    :: KF 
TYPE(TVAB)        ,INTENT(IN)    :: YDVAB
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREM(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXYB(KPROMA,KFLEV,YYTXYB%NDIM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXYBDER(KPROMA,KFLEV,YYTXYBDER%NDIM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE5L(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE5M(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXYB5(KPROMA,KFLEV,YYTXYBT%NDIM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXYBDER5(KPROMA,KFLEV,YYTXYBDERT%NDIM) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZCOEFA(KPROMA), ZCOEFAPL(KPROMA), ZCOEFD(KPROMA)

INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPGRXYBTL',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*    1/ Calculation of "grad delta" at full levels.

IF(LVERTFE) THEN
  DO JLEV=1,KFLEV
    DO JROF=KD,KF
      ZCOEFD(JROF)=&
       & (YDVAB%VDELB(JLEV)*PXYB(JROF,JLEV,YYTXYB%M_RDELP)-PXYB(JROF,JLEV,YYTXYB%M_RTGR))*PXYB5(JROF,JLEV,YYTXYBT%M_LNPR) &
       & +(YDVAB%VDELB(JLEV)*PXYB5(JROF,JLEV,YYTXYBT%M_RDELP)-PXYB5(JROF,JLEV,YYTXYBT%M_RTGR))*PXYB(JROF,JLEV,YYTXYB%M_LNPR)
      PXYBDER(JROF,JLEV,YYTXYBDER%M_LNPRL)=ZCOEFD(JROF)*PRE5L(JROF)+PXYBDER5(JROF,JLEV,YYTXYBDERT%M_COEFD)*PREL(JROF)
      PXYBDER(JROF,JLEV,YYTXYBDER%M_LNPRM)=ZCOEFD(JROF)*PRE5M(JROF)+PXYBDER5(JROF,JLEV,YYTXYBDERT%M_COEFD)*PREM(JROF)
    ENDDO
  ENDDO
ELSE
  DO JLEV=1,KFLEV
    DO JROF=KD,KF
      ZCOEFD(JROF)=-YDVAB%VC(JLEV)*PXYB(JROF,JLEV,YYTXYB%M_RPP)
      PXYBDER(JROF,JLEV,YYTXYBDER%M_LNPRL)=ZCOEFD(JROF)*PRE5L(JROF)+PXYBDER5(JROF,JLEV,YYTXYBDERT%M_COEFD)*PREL(JROF)
      PXYBDER(JROF,JLEV,YYTXYBDER%M_LNPRM)=ZCOEFD(JROF)*PRE5M(JROF)+PXYBDER5(JROF,JLEV,YYTXYBDERT%M_COEFD)*PREM(JROF)
    ENDDO
  ENDDO
ENDIF

!*    2/ Calculation of "grad (alpha + log prehyd)" at full levels
!        discretised as "grad alpha + (grad prehyd) / prehyd ",
!        and calculation of "grad alpha" at full levels.

IF(.NOT.LVERTFE) THEN
  IF(NDLNPR == 0) THEN
    DO JLEV=1,KFLEV
      DO JROF=KD,KF
        ZCOEFAPL(JROF)=YDVAB%VBH(JLEV)*PXYB(JROF,JLEV,YYTXYB%M_RPRE)
        ZCOEFA(JROF)=ZCOEFAPL(JROF)-PXYB(JROF,JLEV,YYTXYB%M_RTGR)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHPLL)=ZCOEFAPL(JROF)*PRE5L(JROF) &
         & + PXYBDER5(JROF,JLEV,YYTXYBDERT%M_COEFAPL)*PREL(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHPLM)=ZCOEFAPL(JROF)*PRE5M(JROF) &
         & + PXYBDER5(JROF,JLEV,YYTXYBDERT%M_COEFAPL)*PREM(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHL)=ZCOEFA(JROF)*PRE5L(JROF) &
         & + PXYBDER5(JROF,JLEV,YYTXYBDERT%M_COEFA)*PREL(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHM)=ZCOEFA(JROF)*PRE5M(JROF) &
         & + PXYBDER5(JROF,JLEV,YYTXYBDERT%M_COEFA)*PREM(JROF)
      ENDDO
    ENDDO
  ELSEIF(NDLNPR == 1) THEN
    DO JLEV=1,KFLEV
      DO JROF=KD,KF
        ZCOEFA(JROF)=-YDVAB%VC(JLEV)*( &
         & PXYB(JROF,JLEV,YYTXYB%M_RPP)*PXYB5(JROF,JLEV,YYTXYBT%M_ALPH)/PXYB5(JROF,JLEV,YYTXYBT%M_LNPR) &
         & + PXYB5(JROF,JLEV,YYTXYBT%M_RPP)*PXYB(JROF,JLEV,YYTXYB%M_ALPH)/PXYB5(JROF,JLEV,YYTXYBT%M_LNPR) &
         & - PXYB5(JROF,JLEV,YYTXYBT%M_RPP)*PXYB5(JROF,JLEV,YYTXYBT%M_ALPH)*PXYB(JROF,JLEV,YYTXYB%M_LNPR) &
         & /(PXYB5(JROF,JLEV,YYTXYBT%M_LNPR)*PXYB5(JROF,JLEV,YYTXYBT%M_LNPR)) )
        ZCOEFAPL(JROF)=ZCOEFA(JROF)+PXYB(JROF,JLEV,YYTXYB%M_RTGR)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHPLL)=ZCOEFAPL(JROF)*PRE5L(JROF) &
         & + PXYBDER5(JROF,JLEV,YYTXYBDERT%M_COEFAPL)*PREL(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHPLM)=ZCOEFAPL(JROF)*PRE5M(JROF) &
         & + PXYBDER5(JROF,JLEV,YYTXYBDERT%M_COEFAPL)*PREM(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHL)=ZCOEFA(JROF)*PRE5L(JROF) &
         & + PXYBDER5(JROF,JLEV,YYTXYBDERT%M_COEFA)*PREL(JROF)
        PXYBDER(JROF,JLEV,YYTXYBDER%M_ALPHM)=ZCOEFA(JROF)*PRE5M(JROF) &
         & + PXYBDER5(JROF,JLEV,YYTXYBDERT%M_COEFA)*PREM(JROF)
      ENDDO
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPGRXYBTL',1,ZHOOK_HANDLE)
END SUBROUTINE GPGRXYBTL
