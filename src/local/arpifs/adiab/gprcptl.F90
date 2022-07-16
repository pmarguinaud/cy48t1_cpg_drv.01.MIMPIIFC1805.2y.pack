SUBROUTINE GPRCPTL(KPROMA,KSTART,KPROF,KFLEV,PQ,PQI,PQL,PQR,PQS,PQG,&
 & PCP,PR,PKAP,PCP5,PR5,PGFL,KGFLTYP,LDTHERMACT)  

!**** *GPRCPTL* - Computes Cp, R and R/Cp from Q, tangent linear

!     Purpose.
!     --------
!           COMPUTES CP AND R  AND R/CP FROM Q

!**   Interface.
!     ----------
!        *CALL* *GPRCPTL(...)*

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

!        INPUT (trajectory):
!          PCP5(KPROMA,KFLEV)   - CP (trajectory)
!          PR5(KPROMA,KFLEV)    - R (trajectory)

!        Implicit arguments :  Physical constants from YOMCST
!        --------------------

!     Method.
!     -------
!        See documentation
!        NOTE! This is NOT an strict linearization of the corresponding
!        non-linear routine, some q-related terms have been supressed

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
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Y.Seity : 04-02-18 : add rain, graupel and snow 
!      L.Isaksen     15-Jun-2005 Supress R dep. on q variables
!      M.Hamrud      15-Jan-2006  Revised GPRCP
!      Y.Bouteloup   06-04-2006 Re-introduction of R dep. on q variables under key
!      M.Jidane      08-04-2006 CY31 Phasing
!      K. Yessad (Jan 2011): more compact rewriting.
!      A. Geer      11-Apr-2016    For independence of observation operator in OOPS, 
!                                  allow calls without YGFL initialised. Removal
!                                  of all YGFL references will have to wait.

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCST   , ONLY : RD, RV, RCPD, RCPV, RCW, RCS
USE YOM_YGFL , ONLY : YGFL
USE YOMDYNA  , ONLY : LDRY_ECMWF

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PQ(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PQI(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PQL(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PQR(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PQS(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PQG(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PGFL(KPROMA,KFLEV,YGFL%NDIM) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(OUT)   :: PCP(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(OUT)   :: PR(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(OUT)   :: PKAP(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PCP5(KPROMA,KFLEV) 
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PR5(KPROMA,KFLEV) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)    :: KGFLTYP
LOGICAL,OPTIONAL,INTENT(IN) :: LDTHERMACT   ! To allow calls independent of YGFL 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZCP(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZR(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZKAP(KPROMA,KFLEV)

INTEGER(KIND=JPIM) :: JLEV, JLON, JGFL
INTEGER(KIND=JPIM) :: IACT(YGFL%NUMFLDS),IPT(YGFL%NUMFLDS),INUMACT,IGFL,IGFLTYP
LOGICAL :: LLGFL,LLQ,LLQL,LLQI,LLQR,LLQS,LLQG,LLKAP
LOGICAL :: LLQ_THERMACT,LLQL_THERMACT,LLQI_THERMACT,LLQR_THERMACT,LLQS_THERMACT,LLQG_THERMACT
REAL(KIND=JPRB) :: ZPCP_R5,ZPCP_R2
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPRCPTL',0,ZHOOK_HANDLE)
ASSOCIATE(NDIM=>YGFL%NDIM, NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & YG=>YGFL%YG, YI=>YGFL%YI, YL=>YGFL%YL, YQ=>YGFL%YQ, YR=>YGFL%YR, YS=>YGFL%YS)

!     ------------------------------------------------------------------

!*       1.    COMPUTES R AND CP AND KAPPA.
!              ----------------------------

LLGFL = PRESENT(PGFL)
IF(.NOT. LLGFL) THEN
  LLQ  = PRESENT(PQ)
  LLQL = PRESENT(PQL)
  LLQI = PRESENT(PQI)
  LLQR = PRESENT(PQR)
  LLQS = PRESENT(PQS)
  LLQG = PRESENT(PQG)
  IF(PRESENT(LDTHERMACT)) THEN
    LLQ_THERMACT=LDTHERMACT
    LLQL_THERMACT=LDTHERMACT
    LLQI_THERMACT=LDTHERMACT
    LLQR_THERMACT=LDTHERMACT
    LLQS_THERMACT=LDTHERMACT
    LLQG_THERMACT=LDTHERMACT
  ELSE
    LLQ_THERMACT=YQ%LTHERMACT
    LLQL_THERMACT=YL%LTHERMACT
    LLQI_THERMACT=YI%LTHERMACT
    LLQR_THERMACT=YR%LTHERMACT
    LLQS_THERMACT=YS%LTHERMACT
    LLQG_THERMACT=YG%LTHERMACT
  ENDIF
ENDIF

LLKAP=PRESENT(PKAP).AND.PRESENT(PR5).AND.PRESENT(PCP5)

! * compute IPT:
IF(LLGFL) THEN
  IGFLTYP = 0
  IF(PRESENT(KGFLTYP)) IGFLTYP=KGFLTYP
  INUMACT = 0
  DO JGFL=1,NUMFLDS
    IF(YCOMP(JGFL)%LTHERMACT) THEN
      INUMACT = INUMACT+1
      IACT(INUMACT) = JGFL
      IF(IGFLTYP == 0) THEN
        IPT(INUMACT) = YCOMP(JGFL)%MP
      ELSEIF(IGFLTYP == 1) THEN
        IPT(INUMACT) = YCOMP(JGFL)%MP1
      ELSEIF(IGFLTYP == 5) THEN
        IPT(INUMACT) = YCOMP(JGFL)%MP5
      ELSEIF(IGFLTYP == 9) THEN
        IPT(INUMACT) = YCOMP(JGFL)%MP9
      ELSEIF(IGFLTYP == 101) THEN
        IPT(INUMACT) = YCOMP(JGFL)%MP_SL1
      ELSE
        CALL ABOR1('GPRCPTL:UNKNOWN GFL TYPE')
      ENDIF
    ENDIF
  ENDDO
ENDIF

! * compute ZR,ZCP,ZKAP:
IF(LLGFL) THEN
  DO JLEV=1,KFLEV
    IF (LDRY_ECMWF) THEN
      ZR(KSTART:KPROF,JLEV)=0.0_JPRB
    ELSE
      ! take account of active GFL; not yet coded
      CALL ABOR1(' GPRCPTL: case not yet coded')
    ENDIF
    IF(INUMACT == 0) THEN
      DO JLON=KSTART,KPROF
        ZCP(JLON,JLEV) = 0.0_JPRB
      ENDDO
    ENDIF
    DO JGFL=1,INUMACT
      IGFL = IACT(JGFL)
      IF(JGFL == 1) THEN
        DO JLON=KSTART,KPROF
          ZCP(JLON,JLEV) = &
           & (YCOMP(IGFL)%RCP-RCPD)*PGFL(JLON,JLEV,IPT(JGFL))
        ENDDO
      ELSE
        DO JLON=KSTART,KPROF
          ZCP(JLON,JLEV) = ZCP(JLON,JLEV)+&
           & (YCOMP(IGFL)%RCP-RCPD)*PGFL(JLON,JLEV,IPT(JGFL))
        ENDDO
      ENDIF
    ENDDO
    IF (LLKAP) THEN
      DO JLON=KSTART,KPROF
        ZPCP_R5=1.0_JPRB/PCP5(JLON,JLEV)
        ZPCP_R2=ZPCP_R5*ZPCP_R5
        ZKAP(JLON,JLEV)= -PR5(JLON,JLEV)*ZPCP_R2*ZCP(JLON,JLEV)  
      ENDDO
    ENDIF
  ENDDO
ELSE
  DO JLEV=1,KFLEV
    DO JLON=KSTART,KPROF
      ZR(JLON,JLEV)  = 0.0_JPRB
      ZCP(JLON,JLEV) = 0.0_JPRB
    ENDDO
    IF(LLQ .AND. LLQ_THERMACT) THEN
      DO JLON=KSTART,KPROF
        IF (.NOT. LDRY_ECMWF) ZR(JLON,JLEV)=ZR(JLON,JLEV)+(RV-RD)*PQ(JLON,JLEV)
        ZCP(JLON,JLEV) = ZCP(JLON,JLEV)+(RCPV-RCPD)*PQ(JLON,JLEV)
      ENDDO
    ENDIF
    IF(LLQL .AND. LLQL_THERMACT) THEN
      DO JLON=KSTART,KPROF
        IF (.NOT. LDRY_ECMWF) ZR(JLON,JLEV)=ZR(JLON,JLEV)-RD*PQL(JLON,JLEV)
        ZCP(JLON,JLEV) = ZCP(JLON,JLEV)+(RCW-RCPD)*PQL(JLON,JLEV)
      ENDDO
    ENDIF
    IF(LLQI .AND. LLQI_THERMACT) THEN
      DO JLON=KSTART,KPROF
        IF (.NOT. LDRY_ECMWF) ZR(JLON,JLEV)=ZR(JLON,JLEV)-RD*PQI(JLON,JLEV)
        ZCP(JLON,JLEV) = ZCP(JLON,JLEV)+(RCS-RCPD)*PQI(JLON,JLEV)
      ENDDO
    ENDIF
    IF(LLQR .AND. LLQR_THERMACT) THEN
      DO JLON=KSTART,KPROF
        IF (.NOT. LDRY_ECMWF) ZR(JLON,JLEV)=ZR(JLON,JLEV)-RD*PQR(JLON,JLEV)
        ZCP(JLON,JLEV) = ZCP(JLON,JLEV)+(RCW-RCPD)*PQR(JLON,JLEV)
      ENDDO
    ENDIF
    IF(LLQS .AND. LLQS_THERMACT) THEN
      DO JLON=KSTART,KPROF
        IF (.NOT. LDRY_ECMWF) ZR(JLON,JLEV)=ZR(JLON,JLEV)-RD*PQS(JLON,JLEV)
        ZCP(JLON,JLEV) = ZCP(JLON,JLEV)+(RCS-RCPD)*PQS(JLON,JLEV)
      ENDDO
    ENDIF
    IF(LLQG .AND. LLQG_THERMACT) THEN
      DO JLON=KSTART,KPROF
        IF (.NOT. LDRY_ECMWF) ZR(JLON,JLEV)=ZR(JLON,JLEV)-RD*PQG(JLON,JLEV)
        ZCP(JLON,JLEV) = ZCP(JLON,JLEV)+(RCS-RCPD)*PQG(JLON,JLEV)
      ENDDO
    ENDIF
    IF (LLKAP) THEN
      DO JLON=KSTART,KPROF
        ZPCP_R5=1.0_JPRB/PCP5(JLON,JLEV)
        ZPCP_R2=ZPCP_R5*ZPCP_R5
        ZKAP(JLON,JLEV)= -PR5(JLON,JLEV)*ZPCP_R2*ZCP(JLON,JLEV)  
        IF (.NOT. LDRY_ECMWF) ZKAP(JLON,JLEV)=ZKAP(JLON,JLEV)+ZR(JLON,JLEV)*ZPCP_R5
      ENDDO
    ENDIF
  ENDDO
ENDIF

! * fill PR,PCP,PKAP:
IF (PRESENT(PR)) PR(KSTART:KPROF,1:KFLEV)=ZR(KSTART:KPROF,1:KFLEV)
IF (PRESENT(PCP)) PCP(KSTART:KPROF,1:KFLEV)=ZCP(KSTART:KPROF,1:KFLEV)
IF (LLKAP) PKAP(KSTART:KPROF,1:KFLEV)=ZKAP(KSTART:KPROF,1:KFLEV)

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPRCPTL',1,ZHOOK_HANDLE)
END SUBROUTINE GPRCPTL
