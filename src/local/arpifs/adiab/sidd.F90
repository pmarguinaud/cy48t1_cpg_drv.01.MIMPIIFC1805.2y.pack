SUBROUTINE SIDD(LDVERTFE, YDDYNA, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,PDH,PDV,PRNH,PT,PSP,KNLON)

!**** * SIDD* - Provide 3D divergence increment in semi-implicit
!               scheme for the case of nonhydrostatic dynamics

!     Purpose.
!     --------
!           Operators gamma, h, L.

!**   Interface.
!     ----------
!        *CALL* *SIDD(...)

!        Explicit arguments :
!        --------------------
!        KLEV   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE VERTICAL
!        KLON   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE SAME LEVEL

!           TYPICAL VALUES ARE  NDLSUR,1  FOR GRID POINT ARRAY
!                               1,NFLSUR  FOR SPECTRAL ARRAY

!        PDH   : HORIZONTAL DIVERGENCE                       - output
!        PDV   : VERTICAL DIVERGENCE VARIABLE                - output
!        PRNH  : NONHYDROSTATIC PRESSURE DEPARTURE VARIABLE  - input
!        PT    : TEMPERATURE                                 - input
!        PSP   : SURFACE PRESSURE                            - input
!        KNLON : NUMBER OF VERTICAL COLUMNS TREATED

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Arpege/Aladin Documentation

!     Author.
!     -------
!      Radmila Bubnova  *GMAP/COMPAS - stage MRE*
!      Original : 93-03-24

!     Modifications.
!     --------------
!      K. Yessad (Sep 2008): update comments + cleanings.
!      K. Yessad (June 2017): Vertical-dependent SITRA.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : TCST
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : TDYNA
USE YOMCVER      , ONLY : LVFE_LAPL_BC


!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL, INTENT (IN) :: LDVERTFE
TYPE(TDYNA), INTENT (IN) :: YDDYNA
TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNLON 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDH(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDV(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRNH(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP(KNLON) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZW3D(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) :: ZW2D(KNLON)

INTEGER(KIND=JPIM) :: ILEV, ILON, INLON, JAD, JLEV, JLON

REAL(KIND=JPRB) :: ZDIM
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "siseve.intfb.h"
#include "sigam.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SIDD',0,ZHOOK_HANDLE)

ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,   SITR=>YDDYN%SITR)
!     ------------------------------------------------------------------

!*       1.    SIGAM ON PT AND PSP: AS IN HYDROSTATIC.
!              ---------------------------------------

ILEV = KLEV
ILON = KLON
INLON= KNLON

CALL SIGAM(LDVERTFE, YDCST, YDGEOMETRY,YDDYN,ILEV,ILON,PDH,PT,PSP,INLON,NFLEVG)

!     ------------------------------------------------------------------

!*       2.    SIGAM ON PRNH, GAMMA APPLIED ONLY, FINISH PDH.
!              ----------------------------------------------

DO JLON = 1, INLON
  ZW2D(JLON)  = 0.0_JPRB
ENDDO

CALL SIGAM(LDVERTFE, YDCST, YDGEOMETRY,YDDYN,ILEV,ILON,ZW3D,PRNH,ZW2D,INLON,NFLEVG)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JAD)
DO JLEV = 1, NFLEVG
  DO JAD=1+(JLEV-1)*ILEV,1+(JLEV-1)*ILEV+(INLON-1)*ILON,ILON
    PDH(JAD) = PDH(JAD) - SITR*ZW3D(JAD) + YDCST%RD*SITR*PRNH(JAD)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

!     ------------------------------------------------------------------

!*       3.    SISEVE ON PRNH, FINISH PDV.
!              --------------------------

CALL SISEVE(LDVERTFE, LVFE_LAPL_BC, YDDYNA,YDGEOMETRY,YDDYN,ILEV,ILON,PRNH,ZW3D,INLON)

ZDIM = YDCST%RG*YDCST%RG/(YDCST%RD*SITR)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JAD)
DO JLEV = 1, NFLEVG
  DO JAD=1+(JLEV-1)*ILEV,1+(JLEV-1)*ILEV+(INLON-1)*ILON,ILON
    PDV(JAD) = ZDIM*ZW3D(JAD)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

!     ------------------------------------------------------------------
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('SIDD',1,ZHOOK_HANDLE)
END SUBROUTINE SIDD
