SUBROUTINE SIPTP(YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,PDH,PDV,PRNH,PT,PSP,KNLON)

!**** *SIPTP* - Counterpart of SITNU in nonhydrostatic

!     Purpose.
!     --------
!           Provide corrections of temperature, surface
!           pressure and nonhydrostatic perturbation from
!           3D divergence in semiimplicit

!**   Interface.
!     ----------
!        *CALL* *SIPTP(...)

!        Explicit arguments :
!        --------------------
!        KLEV   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE VERTICAL
!        KLON   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE SAME LEVEL

!           TYPICAL VALUES ARE  NDLSUR,1  FOR GRID POINT ARRAY
!                               1,NFLSUR  FOR SPECTRAL ARRAY

!        PDH   : HORIZONTAL DIVERGENCE                       - input
!        PDV   : VERTICAL DIVERGENCE VARIABLE                - input
!        PRNH  : NONHYDROSTATIC PRESSURE DEPARTURE VARIABLE  - output
!        PT    : TEMPERATURE                                 - output
!        PSP   : LOG(SURFACE PRESSURE)                       - output
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
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad (Sep 2008): update comments + cleanings.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : TCST
USE YOMDYN       , ONLY : TDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNLON 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDH(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDV(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRNH(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PT(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSP(KNLON) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ILEV, ILON, INLON, JAD, JLEV

REAL(KIND=JPRB) :: Z3DIVC
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "sitnu.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SIPTP',0,ZHOOK_HANDLE)

ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,   SITR=>YDDYN%SITR)
!     ------------------------------------------------------------------

!*       1.    SITNU ON PDH TO GET PSP AND PART OF PRNH
!              ----------------------------------------

ILEV = KLEV
ILON = KLON
INLON= KNLON

CALL SITNU(YDCST, YDGEOMETRY,YDDYN,ILEV,ILON,PDH,PRNH,PSP,INLON)

!     ------------------------------------------------------------------

!*       2.    FINISH PT AND PRNH CALCULATION.
!              -------------------------------

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JAD,Z3DIVC)
DO JLEV = 1, NFLEVG
!DEC$ IVDEP
  DO JAD=1+(JLEV-1)*ILEV,1+(JLEV-1)*ILEV+(INLON-1)*ILON,ILON
    Z3DIVC = (PDH(JAD) + PDV(JAD))/YDCST%RCVD
    PT(JAD) = YDCST%RD*SITR*Z3DIVC
    PRNH(JAD) = YDCST%RCPD*Z3DIVC - YDCST%RCPD*PRNH(JAD)/(YDCST%RD*SITR)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

!     ------------------------------------------------------------------
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('SIPTP',1,ZHOOK_HANDLE)
END SUBROUTINE SIPTP
