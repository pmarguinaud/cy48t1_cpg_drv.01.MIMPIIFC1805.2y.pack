SUBROUTINE COORD_DTS(KX,KY,PLON,PLAT)

!**** *COORD_DTS*

!     PURPOSE.
!     --------

!     Setting the latitudes and longitudes of the dataset grid.

!**   INTERFACE.
!     ----------

!      CALL COORD_DTS(KX,KY,PLON,PLAT)

!        KX   = nb of longitudes in the initial grid
!        KY   = nb of latitudes in the initial grid
!        PLON = longitudes of the initial grid
!        PLAT = latitudes of the initial grid

!     AUTHOR.
!     -------

!         F. TAILLEFER  02/12/2005

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCT0   , ONLY : LELAM
USE YOMCST   , ONLY : RPI
USE YOMCLI   , ONLY : YRCLI

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) ::KX,KY

REAL(KIND=JPRB),INTENT(OUT) :: PLON(KX,KY),PLAT(KX,KY)

INTEGER(KIND=JPIM) ::JX,JY

REAL(KIND=JPRB) :: ZCONV
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('COORD_DTS',0,ZHOOK_HANDLE)

IF (KX /= YRCLI%NDATX .OR. KY /= YRCLI%NDATY) CALL ABOR1('COORD_DTS : PB in dataset dimension')

IF (.NOT.YRCLI%LGLOBE) CALL ABOR1(' COORD_DTS : only for global dataset !')

ZCONV=RPI/180._JPRB
DO JY = 1,YRCLI%NDATY
  DO JX = 1,YRCLI%NDATX
    PLON(JX,JY) = (YRCLI%ELONSW + (0.5_JPRB*YRCLI%EDLON) + (JX-1)*YRCLI%EDLON)*ZCONV
    IF (LELAM) THEN
      PLAT(JX,YRCLI%NDATY+1-JY) = (YRCLI%ELATNE - (0.5_JPRB*YRCLI%EDLAT) - (JY-1)*YRCLI%EDLAT)*ZCONV
    ELSE
      PLAT(JX,JY) = (YRCLI%ELATNE - (0.5_JPRB*YRCLI%EDLAT) - (JY-1)*YRCLI%EDLAT)*ZCONV
    ENDIF
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('COORD_DTS',1,ZHOOK_HANDLE)
END SUBROUTINE COORD_DTS
