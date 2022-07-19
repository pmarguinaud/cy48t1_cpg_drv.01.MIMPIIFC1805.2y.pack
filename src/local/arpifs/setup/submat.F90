SUBROUTINE SUBMAT(YDCST,YDGEOMETRY,YDDYN)

!**** *SUBMAT*  Initialize vertical structure matrix

!     Purpose.
!     --------
!           Initialize vertical structure matrix used for semi-implicit
!         scheme and the normal mode initialization.

!**   Interface.
!     ----------
!        *CALL* *SUBMAT

!        Explicit arguments :  None
!        --------------------
!        Implicit arguments :
!        --------------------
!           Matrix SIB in YOMDYN
!     Method.
!     -------
!        See documentation

!     Externals.    SIGAM - gamma operator of the semi-implicit
!     ----------    SITNU - tau and nu operators of the semi-implicit

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-10-15
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMDYN   , ONLY : TDYN
USE YOMCST   , ONLY : TCST
USE YOMCVER      , ONLY : LVERTFE

IMPLICIT NONE

!     ZD   : DIVERGENCE ARRAY
!     ZT   : TEMPERATURE ARRAY
!     ZSP  : SURFACE PRESSURE ARRAY

TYPE(TCST),     INTENT(IN)    :: YDCST
TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT) :: YDDYN
REAL(KIND=JPRB) :: ZD(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZT(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSP(YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: JLEV, JLON
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "sigam.intfb.h"
#include "sitnu.intfb.h"

!      ----------------------------------------------------------------

!*       1.    INITIALIZE.
!              -----------

IF (LHOOK) CALL DR_HOOK('SUBMAT',0,ZHOOK_HANDLE)

ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,   SIB=>YDDYN%SIB)
DO JLEV=1,NFLEVG
  DO JLON=1,NFLEVG
    ZD(JLON,JLEV)=1.0_JPRB-MIN(ABS(REAL(JLON-JLEV,JPRB)),1.0_JPRB)
  ENDDO
ENDDO

!      ---------------------------------------------------------------

!*       2.    COMPUTES THE MATRIX.
!              --------------------

CALL SITNU(LVERTFE, YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZD,ZT,ZSP,NFLEVG)
CALL SIGAM(LVERTFE, YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZD,ZT,ZSP,NFLEVG,NFLEVG)
!      ----------------------------------------------------------------

!*       3.    STORAGE OF MATRIX.
!              ------------------

DO JLEV=1,NFLEVG
  DO JLON=1,NFLEVG
    SIB(JLON,JLEV)=ZD(JLON,JLEV)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('SUBMAT',1,ZHOOK_HANDLE)
END SUBROUTINE SUBMAT
