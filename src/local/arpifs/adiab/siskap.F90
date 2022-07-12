SUBROUTINE SISKAP(YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,KNLON,PIN,POUT)

!**** *SISKAP* - Semi-implicit scheme: compute Sstar_kap X

!     Purpose.
!     --------
!           Semi-implicit scheme: compute Sstar_kap X
!           This routine is called in the linear system of NHQE model

!           Sstar_kap Z_l = Z_l - [ ((Cpd/Rd)-1)/Tstar ] (tau Z)_l

!**   Interface.
!     ----------
!        *CALL* *SISKAP(...)

!        Explicit arguments :
!        --------------------
!        * INPUT:
!          YDGEOMETRY    : structure containing geometry
!          YDDYN         : structure containing dynamics
!          KLEV,KLON     : couple of dimensions; ex: (NPROMA,1) for grid-point array; (1,NFLEVG) for spectral array
!          KNLON         : number of vertical columns treated
!          PIN           : input array

!        * OUTPUT:
!          POUT          : output array

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
!        Documentation (IDSI)

!     Author.
!     -------
!      K. Yessad (METEO-FRANCE/CNRM/GMAP)
!      Original : march 2017

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYN       , ONLY : TDYN
USE YOMCST       , ONLY : TCST
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KNLON
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIN(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POUT(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZTAUZ(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB) :: ZSP(KNLON)
REAL(KIND=JPRB) :: ZFAC
INTEGER(KIND=JPIM) :: IDT, JLEV, JLON
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "sitnu.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SISKAP',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

ZFAC=((YDCST%RCPD/YDCST%RD)-1.0_JPRB)/YDDYN%SITR

! * Compute "tau Z" (output ZSP is not used):
CALL SITNU(YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,PIN,ZTAUZ,ZSP,KNLON)

! * Compute "Sstar_kap Z":
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLON,JLEV,IDT)
  DO JLON=1,KNLON
    DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      POUT(IDT)=PIN(IDT)-ZFAC*ZTAUZ(IDT)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

!      -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SISKAP',1,ZHOOK_HANDLE)
END SUBROUTINE SISKAP
