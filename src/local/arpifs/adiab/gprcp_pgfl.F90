SUBROUTINE GPRCP_PGFL(KPROMA,KST,KPROF,KFLEV,PCP,PR,PKAP,PGFL,KGFLTYP)

!**** *GPRCP_PGFL* - Computes Cp, R and R/Cp from Q

!     Purpose.
!     --------
!        Computes Cp, R and R/Cp from Q

!**   Interface.
!     ----------
!        *CALL* *GPRCP_PGFL(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KPROMA               - dimensioning.
!          KSTART               - start of work.
!          KPROF                - depth of work.
!          KFLEV                - number of layers.
!          PQ(KPROMA,KFLEV)     - specific humidity.
!          PQI(KPROMA,KFLEV)    - ice.
!          PQL(KPROMA,KFLEV)    - liquid water.
!          PQR(KPROMA,KFLEV)    - rain.
!          PQS(KPROMA,KFLEV)    - snow.
!          PQG(KPROMA,KFLEV)    - graupel.

!        OUTPUT:
!          PCP(KPROMA,KFLEV)    - CP
!          PR(KPROMA,KFLEV)     - R
!          PKAP(KPROMA,KFLEV)   - KAPPA

!        Implicit arguments :  Physical constants from YOMCST
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 88-02-04

!     Modifications.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Y.Seity  04-02-13 (Rain, Snow and Graupel)
!      M.Hamrud  15-Jan-2006  Revised GPRCP
!      K. Yessad (Jan 2011): more compact rewriting.
!      R. El Khatib 28-Aug-2014 Optimizations :
!       - compute R or CP only if required
!       - loop collapsing whenever possible, through pure array syntax
!      A. Geer      01-Oct-2015    For independence of observation operator in OOPS, 
!                                  allow calls without YGFL initialised. Removal
!                                  of all YGFL references will have to wait.
!      H Petithomme (Dec 2020): general rewrite for optimization
!     ------------------------------------------------------------------

USE PARKIND1,ONLY: JPIM,JPRB
USE YOMHOOK,ONLY: LHOOK,DR_HOOK
USE YOMCST,ONLY: RD,RV,RCPD,RCPV,RCW,RCS
USE YOM_YGFL,ONLY: TYPE_GFLD,TYPE_GFL_COMP,YGFL
USE CODETOOLS,ONLY: SETDEFAULTL,SETDEFAULTI

IMPLICIT NONE

#include "abor1.intfb.h"

INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA,KFLEV,KST,KPROF
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KGFLTYP
REAL(KIND=JPRB),OPTIONAL,INTENT(IN) :: PGFL(KPROMA,KFLEV,YGFL%NDIM)
REAL(KIND=JPRB),OPTIONAL,TARGET,INTENT(OUT) :: PCP(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,TARGET,INTENT(OUT) :: PR(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT) :: PKAP(KPROMA,KFLEV)

REAL(KIND=JPRB),TARGET :: ZCP0(KPROMA,KFLEV),ZR0(KPROMA,KFLEV)
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZCP(:,:),ZR(:,:)
REAL(KIND=JPRB) :: ZHOOK

IF (LHOOK) CALL DR_HOOK("GPRCP_PGFL",0,ZHOOK)

IF (PRESENT(PR)) THEN
  ZR => PR(:,:)
ELSEIF (PRESENT(PKAP)) THEN
  ZR => ZR0(:,:)
ENDIF

IF (PRESENT(PCP)) THEN
  ZCP => PCP(:,:)
ELSEIF (PRESENT(PKAP)) THEN
  ZCP => ZCP0(:,:)
ELSE
  ZCP => NULL()
ENDIF

!DIR$ FORCEINLINE
CALL GPRCPGFL(YGFL,KPROMA,KST,KPROF,KFLEV,PGFL,ZR,ZCP,KGFLTYP)

IF (PRESENT(PKAP)) PKAP(KST:KPROF,1:KFLEV)=ZR(KST:KPROF,1:KFLEV)/ZCP(KST:KPROF,1:KFLEV)

IF (LHOOK) CALL DR_HOOK("GPRCP_PGFL",1,ZHOOK)

CONTAINS

SUBROUTINE GPRCPGFL(YGFL,KPROMA,KST,KPROF,KFLEV,PGFL,ZR,ZCP,KGFLTYP)
  TYPE(TYPE_GFLD),INTENT(IN) :: YGFL
  INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA,KFLEV,KST,KPROF
  INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KGFLTYP
  REAL(KIND=JPRB),INTENT(IN) :: PGFL(KPROMA,KFLEV,YGFL%NDIM)
  REAL(KIND=JPRB),CONTIGUOUS,POINTER,INTENT(INOUT) :: ZR(:,:)
  REAL(KIND=JPRB),CONTIGUOUS,POINTER,INTENT(INOUT) :: ZCP(:,:)

  INTEGER(KIND=JPIM) :: IGFLTYP,JL,JP,JGFL,I,NACT,IPT(YGFL%NUMFLDS)
  REAL(KIND=JPRB) :: YR(YGFL%NUMFLDS),YCP(YGFL%NUMFLDS)

  ASSOCIATE(YCOMP=>YGFL%YCOMP)

  IGFLTYP = SETDEFAULTI(0,KGFLTYP)

  NACT = 0
  IF (IGFLTYP == 0) THEN
    DO JGFL=1,YGFL%NUMFLDS
      IF (.NOT.YCOMP(JGFL)%LTHERMACT) CYCLE

      NACT = NACT+1
      IPT(NACT) = YCOMP(JGFL)%MP
      YR(NACT) = YCOMP(JGFL)%R
      YCP(NACT) = YCOMP(JGFL)%RCP
    ENDDO
  ELSEIF (IGFLTYP == 1) THEN
    DO JGFL=1,YGFL%NUMFLDS
      IF (.NOT.YCOMP(JGFL)%LTHERMACT) CYCLE

      NACT = NACT+1
      IPT(NACT) = YCOMP(JGFL)%MP1
      YR(NACT) = YCOMP(JGFL)%R
      YCP(NACT) = YCOMP(JGFL)%RCP
    ENDDO
  ELSEIF (IGFLTYP == 5) THEN
    DO JGFL=1,YGFL%NUMFLDS
      IF (.NOT.YCOMP(JGFL)%LTHERMACT) CYCLE

      NACT = NACT+1
      IPT(NACT) = YCOMP(JGFL)%MP5
      YR(NACT) = YCOMP(JGFL)%R
      YCP(NACT) = YCOMP(JGFL)%RCP
    ENDDO
  ELSEIF (IGFLTYP == 9) THEN
    DO JGFL=1,YGFL%NUMFLDS
      IF (.NOT.YCOMP(JGFL)%LTHERMACT) CYCLE

      NACT = NACT+1
      IPT(NACT) = YCOMP(JGFL)%MP9_PH
      YR(NACT) = YCOMP(JGFL)%R
      YCP(NACT) = YCOMP(JGFL)%RCP
    ENDDO
  ELSEIF (IGFLTYP == 101) THEN
    DO JGFL=1,YGFL%NUMFLDS
      IF (.NOT.YCOMP(JGFL)%LTHERMACT) CYCLE

      NACT = NACT+1
      IPT(NACT) = YCOMP(JGFL)%MP_SL1
      YR(NACT) = YCOMP(JGFL)%R
      YCP(NACT) = YCOMP(JGFL)%RCP
    ENDDO
  ELSE
    CALL ABOR1("UNKNOWN GFL TYPE")
  ENDIF

  IF (ASSOCIATED(ZR)) THEN
    ! note: mind the brackets on PGFL, in order to add RD at last (precise computes)
    ! optim: explicit treatment of cases 0, 1, 2 and 3 (fewer loops)
    IF (NACT == 0) THEN
      ZR(KST:KPROF,1:KFLEV) = RD
    ELSE IF (NACT == 1) THEN
      ZR(KST:KPROF,1:KFLEV) = RD+(YR(1)-RD)*PGFL(KST:KPROF,1:KFLEV,IPT(1))
    ELSE IF (NACT == 2) THEN
      ZR(KST:KPROF,1:KFLEV) = RD+((YR(1)-RD)*PGFL(KST:KPROF,1:KFLEV,IPT(1))+&
       & (YR(2)-RD)*PGFL(KST:KPROF,1:KFLEV,IPT(2)))
    ELSE IF (NACT == 3) THEN
      ZR(KST:KPROF,1:KFLEV) = RD+((YR(1)-RD)*PGFL(KST:KPROF,1:KFLEV,IPT(1))+&
       & (YR(2)-RD)*PGFL(KST:KPROF,1:KFLEV,IPT(2))+&
       & (YR(3)-RD)*PGFL(KST:KPROF,1:KFLEV,IPT(3)))
    ELSE
      ZR(KST:KPROF,1:KFLEV) = (YR(1)-RD)*PGFL(KST:KPROF,1:KFLEV,IPT(1))+&
       & (YR(2)-RD)*PGFL(KST:KPROF,1:KFLEV,IPT(2))+&
       & (YR(3)-RD)*PGFL(KST:KPROF,1:KFLEV,IPT(3))

      DO I=4,NACT-1
        ZR(KST:KPROF,1:KFLEV) = ZR(KST:KPROF,1:KFLEV)+&
          (YR(I)-RD)*PGFL(KST:KPROF,1:KFLEV,IPT(I))
      ENDDO

      I=NACT
      ZR(KST:KPROF,1:KFLEV) = RD+(ZR(KST:KPROF,1:KFLEV)+&
          (YR(I)-RD)*PGFL(KST:KPROF,1:KFLEV,IPT(I)))
    ENDIF
  ENDIF

  IF (ASSOCIATED(ZCP)) THEN
    ! note: mind the brackets on PGFL, in order to add RCPD at last (precise computes)
    ! optim: explicit treatment of cases 0, 1, 2 and 3 (fewer loops)
    IF (NACT == 0) THEN
      ZCP(KST:KPROF,1:KFLEV) = RCPD
    ELSE IF (NACT == 1) THEN
      ZCP(KST:KPROF,1:KFLEV) = RCPD+(YCP(1)-RCPD)*PGFL(KST:KPROF,1:KFLEV,IPT(1))
    ELSE IF (NACT == 2) THEN
      ZCP(KST:KPROF,1:KFLEV) = RCPD+((YCP(1)-RCPD)*PGFL(KST:KPROF,1:KFLEV,IPT(1))+&
       & (YCP(2)-RCPD)*PGFL(KST:KPROF,1:KFLEV,IPT(2)))
    ELSE IF (NACT == 3) THEN
      ZCP(KST:KPROF,1:KFLEV) = RCPD+((YCP(1)-RCPD)*PGFL(KST:KPROF,1:KFLEV,IPT(1))+&
       & (YCP(2)-RCPD)*PGFL(KST:KPROF,1:KFLEV,IPT(2))+&
       & (YCP(3)-RCPD)*PGFL(KST:KPROF,1:KFLEV,IPT(3)))
    ELSE
      ZCP(KST:KPROF,1:KFLEV) = (YCP(1)-RCPD)*PGFL(KST:KPROF,1:KFLEV,IPT(1))+&
       & (YCP(2)-RCPD)*PGFL(KST:KPROF,1:KFLEV,IPT(2))+&
       & (YCP(3)-RCPD)*PGFL(KST:KPROF,1:KFLEV,IPT(3))

      DO I=4,NACT-1
        ZCP(KST:KPROF,1:KFLEV) = ZCP(KST:KPROF,1:KFLEV)+(YCP(I)-RCPD)*&
         & PGFL(KST:KPROF,1:KFLEV,IPT(I))
      ENDDO

      I=NACT
      ZCP(KST:KPROF,1:KFLEV) = RCPD+(ZCP(KST:KPROF,1:KFLEV)+(YCP(I)-RCPD)*&
       & PGFL(KST:KPROF,1:KFLEV,IPT(I)))
    ENDIF
  ENDIF
  END ASSOCIATE
END SUBROUTINE

END SUBROUTINE

