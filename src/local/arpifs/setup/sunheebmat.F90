SUBROUTINE SUNHEEBMAT(YDCST, YDGEOMETRY,YDDYN)

!**** *SUNHEEBMAT*  Initialize vertical structure matrix

!     Purpose.
!     --------
!           Initialize vertical structure matrix used for NHEE model semi-implicit scheme.
!           Unknown in Helmholtz equation is assumed to be horizontal divergence. 

!**   Interface.
!     ----------
!        *CALL* *SUNHEEBMAT_NEW

!     Explicit arguments :
!     --------------------
!      * INPUT:
!        YDGEOMETRY   : structure containing geometry

!      * INOUT:
!        YDDYN        : structure containing dynamics

!        Implicit arguments :
!        --------------------
!           Matrix YDDYN%SIB
!           Matrices YDDYN%SIFAC and YDDYN%SIFACI

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
!        K. YESSAD & F. VOITUS (after SUNHBMAT and SUNHQEBMAT).
!        Original : january 2018

!     Modifications.
!     -------------- 
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYN       , ONLY : TDYN
USE YOMCST       , ONLY : TCST
USE YOMDYNA      , ONLY : YRDYNA

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)     :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(INOUT)  :: YDDYN

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZID(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZT(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSP(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZZERO(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSIB_HYD(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSIB_ADD(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZZ1(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZZ2(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZZ3(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZC2,ZH2,ZBT2,ZCOEF

INTEGER(KIND=JPIM) :: JLEV, JLON
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "minv_caller.h"
#include "sigam.intfb.h"
#include "sitnu.intfb.h"
#include "siseve.intfb.h"

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUNHEEBMAT',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & RBTS2=>YDDYN%RBTS2,SITIME=>YDDYN%SITIME, &
 & SIB=>YDDYN%SIB,SIFAC=>YDDYN%SIFAC,SIFACI=>YDDYN%SIFACI,SITR=>YDDYN%SITR)

!      ----------------------------------------------------------------

!*       1.    INITIALIZE.
!              -----------

! * ZID contains the identity matrix "I":
DO JLEV=1,NFLEVG
  DO JLON=1,NFLEVG
    ZID(JLON,JLEV)=1.0_JPRB-MIN(ABS(REAL(JLON-JLEV,JPRB)),1.0_JPRB)
  ENDDO
ENDDO

! * Compute the constants C**2, H**2:
ZC2=YDCST%RD*SITR*YDCST%RCPD/YDCST%RCVD
ZH2=((YDCST%RD*SITR)/YDCST%RG)**2

! * Compute ZBT2 = (1+epsilon_uncentering)**2 * BETADT**2 * TSTEP**2:
ZBT2=(RBTS2*SITIME)**2

!      ---------------------------------------------------------------

!*       2.    COMPUTES THE MATRIX.
!              --------------------

! * 2.1  Computes hydrostatic part of SIB (cf. SUSI)

CALL SITNU(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZID,ZT,ZSP,NFLEVG)
CALL SIGAM(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZSIB_HYD,ZT,ZSP,NFLEVG,NFLEVG)

! * 2.2  Computes SIFAC and SIFACI

! * Apply the "Laplacian" operator LLstarstar to ZID.
CALL SISEVE(YRDYNA,YDGEOMETRY,YDDYN,1,NFLEVG,ZID,SIFAC,NFLEVG)

! * Multiply by
! "-(1+epsilon_uncentering)**2 beta**2 (Delta t)**2 C**2 (1/H**2)"
! then add the identity matrix:
ZCOEF=ZBT2*ZC2/ZH2
DO JLEV=1,NFLEVG
  DO JLON=1,NFLEVG
    SIFAC(JLON,JLEV)=ZID(JLON,JLEV)-ZCOEF*SIFAC(JLON,JLEV)
  ENDDO
ENDDO

! * Inverts matrix SIFAC:
CALL MINV_CALLER(.TRUE.,NFLEVG,SIFAC,SIFACI)

! * 2.3  Computes "anhydrostatic fully-elastic" part of SIB

! * compute ZT=tau*ZID
CALL SITNU(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZID,ZT,ZSP,NFLEVG)

! * compute ZZ1 = C**2 * ZID - Cpdry * [tau*ZID]
DO JLEV=1,NFLEVG
  DO JLON=1,NFLEVG
    ZZ1(JLON,JLEV)=ZC2*ZID(JLON,JLEV)-YDCST%RCPD*ZT(JLON,JLEV)
  ENDDO
ENDDO

! * compute ZZ2 = SIFACI * ZZ1
CALL MXMAOP(SIFACI,1,NFLEVG,ZZ1,1,NFLEVG,ZZ2,1,NFLEVG,NFLEVG,NFLEVG,NFLEVG)

! * compute ZZ3 = gamma*ZZ2
ZZERO(1:NFLEVG)=0.0_JPRB
CALL SIGAM(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZZ3,ZZ2,ZZERO,NFLEVG,NFLEVG)

! Compute ZSIB_ADD = ZZ2 - [SITR/C**2] * ZZ3
DO JLEV=1,NFLEVG
  DO JLON=1,NFLEVG
    ZSIB_ADD(JLON,JLEV)=ZZ2(JLON,JLEV)-(SITR/ZC2)*ZZ3(JLON,JLEV)
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
IF (LHOOK) CALL DR_HOOK('SUNHEEBMAT',1,ZHOOK_HANDLE)
END SUBROUTINE SUNHEEBMAT
