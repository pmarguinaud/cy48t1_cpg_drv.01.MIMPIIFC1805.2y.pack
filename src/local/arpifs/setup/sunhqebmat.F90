SUBROUTINE SUNHQEBMAT(YDCST, YDGEOMETRY,YDDYN)

!**** *SUNHQEBMAT*  Initialize vertical structure matrix

!     Purpose.
!     --------
!           Initialize vertical structure matrix used for NHQE model semi-implicit scheme.

!**   Interface.
!     ----------
!        *CALL* *SUNHQEBMAT

!     Explicit arguments :
!     --------------------
!      * INPUT:
!        YDGEOMETRY   : structure containing geometry

!      * INOUT:
!        YDDYN        : structure containing dynamics

!        Implicit arguments :
!        --------------------
!           Matrix YDDYN%SIB
!           Matrices YDDYN%SI_ILAPKSSI(.,.,1:2)

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation (IDSI)

!     Author.
!     -------
!        K. YESSAD (after SUSI and SUNHSI).
!        Original : march 2017

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYN       , ONLY : TDYN
USE YOMCST       , ONLY : TCST

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCVER      , ONLY : LVERTFE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)     :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(INOUT)  :: YDDYN

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZID(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZT(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSP(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSIB_HYD(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSIB_ADD(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZH2,ZBT2

INTEGER(KIND=JPIM) :: IOPT1, IOPT2, JLEV, JLON
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "sigam.intfb.h"
#include "sitnu.intfb.h"
#include "silkapi.intfb.h"

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUNHQEBMAT',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & RBTS2=>YDDYN%RBTS2,SITIME=>YDDYN%SITIME, &
 & SIB=>YDDYN%SIB,SI_ILAPKSSI=>YDDYN%SI_ILAPKSSI,SITR=>YDDYN%SITR)

!      ----------------------------------------------------------------

!*       1.    INITIALIZE.
!              -----------

IOPT1=1
IOPT2=2

! * ZID contains the identity matrix "I":
DO JLEV=1,NFLEVG
  DO JLON=1,NFLEVG
    ZID(JLON,JLEV)=1.0_JPRB-MIN(ABS(REAL(JLON-JLEV,JPRB)),1.0_JPRB)
  ENDDO
ENDDO

! * Compute the constant H**2:
ZH2=((YDCST%RD*SITR)/YDCST%RG)**2

! * Compute ZBT2 = (1+epsilon_uncentering)**2 * BETADT**2 * TSTEP**2:
ZBT2=(RBTS2*SITIME)**2

!      ---------------------------------------------------------------

!*       2.    COMPUTES THE MATRIX.
!              --------------------

! * 2.1  Computes hydrostatic part of SIB (cf. SUSI)

CALL SITNU(LVERTFE, YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZID,ZT,ZSP,NFLEVG)
CALL SIGAM(LVERTFE, YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZSIB_HYD,ZT,ZSP,NFLEVG,NFLEVG)

! * 2.2  Computes "LLsstar_kap**-1" and "LLsstar_kap**-1 Sstar_kap"

CALL SILKAPI(YDCST, YDGEOMETRY,YDDYN,IOPT1,1,NFLEVG,NFLEVG,ZID,SI_ILAPKSSI(1,1,1))
CALL SILKAPI(YDCST, YDGEOMETRY,YDDYN,IOPT2,1,NFLEVG,NFLEVG,ZID,SI_ILAPKSSI(1,1,2))

! * 2.3  Computes "anhydrostatic quasi-elastic" part of SIB

DO JLEV=1,NFLEVG
  DO JLON=1,NFLEVG
    ZSIB_ADD(JLON,JLEV)=-(ZH2/ZBT2)*SI_ILAPKSSI(JLON,JLEV,2)
  ENDDO
ENDDO

!      ----------------------------------------------------------------

!*       3.    STORAGE OF MATRIX.
!              ------------------

DO JLEV=1,NFLEVG
  DO JLON=1,NFLEVG
    SIB(JLON,JLEV)=ZSIB_HYD(JLON,JLEV)+ZSIB_ADD(JLON,JLEV)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUNHQEBMAT',1,ZHOOK_HANDLE)
END SUBROUTINE SUNHQEBMAT
