SUBROUTINE SIGAM(LDVERTFE, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,PD,PT,PSP,KNLON,KFLEVG)

!**** *SIGAM* - Solve hydrostatic operator in semi-implicit

!     Purpose.
!     --------
!           Operator gamma to compute p.

!**   Interface.
!     ----------
!        *CALL* *SIGAM(...)

!        Explicit arguments :
!        --------------------
!        KLEV   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE VERTICAL
!        KLON   : DISTANCE IN MEMORY BETWEEN VALUES OF THE DIVERGENCE
!                OR TEMPERATURE AT THE SAME LEVEL

!           TYPICAL VALUES ARE  NDLSUR,1  FOR GRID POINT ARRAY
!                               1,NFLSUR  FOR SPECTRAL ARRAY

!        PD    : DIVERGENCE       (output)
!        PT    : TEMPERATURE      (input)
!        PSP   : SURFACE PRESSURE (input)
!        KNLON : NUMBER OF VERTICAL COLUMNS TREATED
!        KFLEVG: NUMBER OF ELEMENTS IN A VERTICAL COLUMN

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      Modified : 09-Oct-2007 by K. YESSAD: possibility to have a specific
!                 value of LVERTFE in the SI linear model.
!      F. Vana + NEC 28-Apr-2009: OpenMP
!      P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!      G. Mozdzynski Oct 2012: OpenMP optimization
!      K. Yessad (Dec 2016): Prune obsolete options.
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : TCST
USE YOMDYN       , ONLY : TDYN


!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL, INTENT (IN) :: LDVERTFE
TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEVG
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PD(1+(KFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(1+(KFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP(KNLON) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSPHI(KNLON,0:KFLEVG+1)
REAL(KIND=JPRB) :: ZOUT(KNLON,0:KFLEVG)
REAL(KIND=JPRB) :: ZSPHIX(0:KFLEVG,KNLON)
INTEGER(KIND=JPIM) :: IDT, JLEV, JLON
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SIGAM',0,ZHOOK_HANDLE)
ASSOCIATE(YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE)
ASSOCIATE(SIALPH=>YDDYN%SIALPH, SILNPR=>YDDYN%SILNPR, SIRPRG=>YDDYN%SIRPRG)
!     ------------------------------------------------------------------

!*       1.    SUM GEOPOTENTIAL, COMPUTES P AND PUT IT IN PD.
!              ----------------------------------------------

IF(LDVERTFE) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JLON,IDT)
  DO JLON=1,KNLON
    DO JLEV=1,KFLEVG
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      ZSPHI(JLON,JLEV)=-YDCST%RD*PT(IDT)*SILNPR(JLEV)*YDVETA%VFE_RDETAH(JLEV)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

  ZSPHI(1:KNLON,0)=0.0_JPRB
  ZSPHI(1:KNLON,KFLEVG+1)=0.0_JPRB
  CALL VERDISINT(YDVFE,'IBOT','11',KNLON,1,KNLON,KFLEVG,ZSPHI,ZOUT,KCHUNK=YDGEOMETRY%YRDIM%NPROMA)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JLON,IDT)
  DO JLON=1,KNLON
    DO JLEV=1,KFLEVG
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      PD(IDT)=ZOUT(JLON,JLEV-1)+PSP(JLON)*SIRPRG
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLON,JLEV,IDT)
  DO JLON=1,KNLON
    ZSPHIX(KFLEVG,JLON)=0.0_JPRB
    DO JLEV=KFLEVG,1,-1
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      ZSPHIX(JLEV-1,JLON)=ZSPHIX(JLEV,JLON)+YDCST%RD*PT(IDT)*SILNPR(JLEV)
      PD(IDT)=ZSPHIX(JLEV,JLON)+SIALPH(JLEV)*YDCST%RD*PT(IDT)+PSP(JLON)*SIRPRG
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDIF

!      -----------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SIGAM',1,ZHOOK_HANDLE)
END SUBROUTINE SIGAM
