SUBROUTINE SI_CCCOR(YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,KNLON,PIN,POU,POU2)

!**** *SI_CCCOR* - Semi-implicit scheme in the NH model.
!                  Does a multiplication by C**2 * COR when the constraint
!                  C1 is not matched.
!                  Optionally, does a multiplication by (g**2/(N**2 C**2)) Lv* A2
!                  for the case when the constraint C2 is not matched.

!     Purpose.
!     --------
!       This routine is called in the SI scheme (NH model)

!       Action 1/
!        It computes C**2 * COR where:
!         C = sqrt( R_dry Trsi Cp_dry/Cv_dry )
!         C**2 is the acoustic phase speed
!        and:
!         COR = (Cv_dry/(R_d**2 Trsi)) Gamma Tau
!               - (Cv_dry/(R_dry*Cp_dry)) Gamma
!               - (Cv_dry/(R_dry*Trsi)) Tau
!               + (Cv_dry/Cv_dry) Nu
!        This is equivalent to write:
!         C**2 * COR = (Cp_dry/R_dry) Gamma Tau - Trsi Gamma
!          - Cp_dry Tau + (R_d*Trsi) Nu
!        It does a multiplication by C**2 * COR.

!       Action 2/ (when argument POU2 is present)
!        It computes (g**2/(N**2 C**2)) Lv* A2 where:
!         A2 = S*G* - (Cp_dry/Cv_dry) S* - (Cp_dry/Cv_dry) G*
!         S* is (Cp_dry/ (R_dry Trsi)) Tau
!         G* is (1/R_dry) Gamma
!         g is the acceleration of gravity
!         N**2 = g**2/(Cp_dry Trsi) is the Brunt-Vaisala frequency

!       Remarks:
!        COR=0 in the continuous equations and in the discretisations
!         matching the constraint C1 (ex: NDLNPR=1 in finite differences).
!        (g**2/(N**2 C**2)) Lv* A2 = Id in the continuous case.

!**   Interface.
!     ----------
!        *CALL* *SI_CCCOR(...)

!        Explicit arguments :
!        --------------------

!         KLEV     :  distance in memory between values of PIN      (input)
!                     at the same vertical.
!         KLON:    :  distance in memory between values of PIN      (input)
!                     at the same level.
!         For (KLEV,KLON), typical values are:
!                     NDLSUR,1 for grid-point array.
!                     1,NFLSUR for spectral array.
!         KNLON    :  number of vertical columns treated            (input)
!         PIN      :  known input vector.                           (input)
!         POU      :  unknown output vector linked to COR.          (output)
!         POU2     :  unknown output vector linked to Laplacian.    (opt output)

!        Implicit arguments :
!        --------------------
!        none.

!     Method.
!     -------
!        See documentation about semi-implicite scheme

!     Externals.
!     ----------
!        None.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        ARPEGE documentation about semi-implicit scheme

!     Author.
!     -------
!        K. YESSAD (CNRM/GMAP)

!     Modifications.
!     --------------
!      P Smolikova and J Vivoda (Oct 2013): add POU2
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : TCST
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : YRDYNA
USE YOMCVER      , ONLY : LVERTFE

!      ----------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KNLON
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIN(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POU(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)   :: POU2(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 

!      ----------------------------------------------------------------

INTEGER(KIND=JPIM) :: IDT,JLEV,JLON

REAL(KIND=JPRB) :: ZGAM(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) :: ZTAU(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) :: ZGAMTAU(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) :: ZTAUGAM(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) :: ZAUX(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) :: ZNU(KNLON)
REAL(KIND=JPRB) :: ZSP(KNLON)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

#include "sigam.intfb.h"
#include "sitnu.intfb.h"
#include "siseve.intfb.h"

!      ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SI_CCCOR',0,ZHOOK_HANDLE)

ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,   SITR=>YDDYN%SITR)
!      ----------------------------------------------------------------

!*       1.    COMPUTE C2*COR AND DO MULTIPLICATION
!              ------------------------------------

! * Computes ZTAU = Cp_dry * Tau * PIN and ZNU = Nu * PIN
CALL SITNU(LVERTFE, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,PIN,ZTAU,ZNU,KNLON)
ZTAU(:)=ZTAU(:)*YDCST%RCPD

! * Computes ZGAM = Gamma * PIN
ZSP(1:KNLON)=0.0_JPRB
CALL SIGAM(LVERTFE, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,ZGAM,PIN,ZSP,KNLON,NFLEVG)

! * Computes ZGAMTAU = Cp_dry * Gamma * (Tau * PIN) = Gamma * ZTAU
ZSP(1:KNLON)=0.0_JPRB
CALL SIGAM(LVERTFE, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,ZGAMTAU,ZTAU,ZSP,KNLON,NFLEVG)

! * Computes POU:
DO JLEV=1,NFLEVG
  DO JLON=1,KNLON
    IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
    POU(IDT)=ZGAMTAU(IDT)/YDCST%RD-SITR*ZGAM(IDT)-ZTAU(IDT)+YDCST%RD*SITR*ZNU(JLON)
  ENDDO
ENDDO

!      ----------------------------------------------------------------

!*       2.    COMPUTE POU2
!              ------------

IF (PRESENT(POU2)) THEN

  ! * Computes ZTAUGAM = Tau * ZGAM
  ZSP(1:KNLON)=0.0_JPRB
  CALL SITNU(LVERTFE, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,ZGAM,ZTAUGAM,ZSP,KNLON)

  ! * Computes (g**2/(N**2 C**2)) A2 * PIN:
  DO JLEV=1,NFLEVG
    DO JLON=1,KNLON
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      ZAUX(IDT)=((YDCST%RCVD/YDCST%RD)*ZTAUGAM(IDT)-ZTAU(IDT)-SITR*ZGAM(IDT))*YDCST%RCPD/(YDCST%RD*YDCST%RD*SITR)
    ENDDO
  ENDDO

  ! * Apply the tridiagonal operator Lv* to get POU2
  CALL SISEVE(YRDYNA,YDGEOMETRY,YDDYN,KLEV,KLON,ZAUX,POU2,KNLON)

ENDIF

!       ----------------------------------------------------------------

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('SI_CCCOR',1,ZHOOK_HANDLE)
END SUBROUTINE SI_CCCOR

