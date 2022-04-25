SUBROUTINE EGEO923(YDGEM)

!**** *EGEO923*

!     PURPOSE.
!     --------
!      Initialize YEMCLI and YOMDIL in configuration 923 (ALADIN).

!     INTERFACE.
!     ----------
!      CALL EGEO923

!     AUTHORS.
!     --------
!      D. Giard 97-11-04 from INCLI0 (Luc Gerard) and GEO923

!     MODIFICATIONS.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     ------------------------------------------------------------------

USE YOMGEM   , ONLY : TGEM
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMLUN   , ONLY : NULOUT
USE YOMCLI   , ONLY : YRCLI  
USE YEMCLI   , ONLY : NILCLI1  ,NILCLI2  ,NJLCLI1  ,NJLCLI2
USE YOMDIL   , ONLY : NPOTYP   ,SLAPO    ,GLAPO    ,SLOPO    ,GLOPO    ,FACDI

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TGEM) , INTENT(IN) :: YDGEM
INTEGER(KIND=JPIM) :: II1, II2, IJ1, IJ2
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
#include "abor1.intfb.h"
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('EGEO923',0,ZHOOK_HANDLE)
ASSOCIATE(NSTTYP=>YDGEM%NSTTYP, RMUCEN=>YDGEM%RMUCEN, RLOCEN=>YDGEM%RLOCEN, RSTRET=>YDGEM%RSTRET)
!     ------------------------------------------------------------------

!     1. YEMCLI

! Set default values.
NILCLI1= 1
NILCLI2= YRCLI%NDATX
NJLCLI1= 1
NJLCLI2= YRCLI%NDATY

!  Compute and check corners coordinates in global grid.
IF (.NOT.YRCLI%LGLOBE) THEN
  IJ1= NINT((90._JPRB+YRCLI%ELATSW)/YRCLI%EDLAT+0.5_JPRB)
  IJ2= NINT((90._JPRB+YRCLI%ELATNE)/YRCLI%EDLAT-0.5_JPRB)
  II1= NINT((180._JPRB*(1.0_JPRB-SIGN(1.0_JPRB,YRCLI%ELONSW))+YRCLI%ELONSW)/YRCLI%EDLON+0.5_JPRB)
  II2= NINT((180._JPRB*(1.0_JPRB-SIGN(1.0_JPRB,YRCLI%ELONNE))+YRCLI%ELONNE)/YRCLI%EDLON-0.5_JPRB)

  IF ((IJ2-IJ1+1) /= YRCLI%NDATY) THEN
    WRITE(UNIT=NULOUT,FMT=&
     & '(''NO AGREEMENT BETWEEN ELATSW, ELATNE AND NDATY'')')  
    CALL ABOR1('EGEO923: NO AGREEMENT BETWEEN ELATSW, ELATNE AND NDATY')
  ENDIF
  IF ((II2-II1+1 /= YRCLI%NDATX).AND.(II2-II1+1+YRCLI%NGLOBX /= YRCLI%NDATX)) THEN
    WRITE(UNIT=NULOUT,FMT=&
     & '(''NO AGREEMENT BETWEEN ELONSW, ELONNE AND NDATX'')')  
    CALL ABOR1('EGEO923: NO AGREEMENT BETWEEN ELONSW, ELONNE AND NDATX')
  ENDIF

  NILCLI1=II1
  NILCLI2=II2
  NJLCLI1=IJ1
  NJLCLI2=IJ2
ENDIF

! Print.

WRITE(UNIT=NULOUT,FMT='(&
 & '' NILCLI1='',I5,'' NILCLI2='',I5,&
 & '' NJLCLI1='',I5,'' NJLCLI2='',I5)')&
 & NILCLI1,NILCLI2,NJLCLI1,NJLCLI2  

!     2. YOMDIL (used in part 1 of E923, by ERELSPE,...)

SLAPO = RMUCEN
GLAPO = SQRT(1.0_JPRB-RMUCEN**2)
SLOPO = SIN(RLOCEN)
GLOPO = COS(RLOCEN)
FACDI = RSTRET
NPOTYP= NSTTYP

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EGEO923',1,ZHOOK_HANDLE)
END SUBROUTINE EGEO923
