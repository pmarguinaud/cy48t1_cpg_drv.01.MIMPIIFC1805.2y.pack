SUBROUTINE SILKAP(LDVERTFE, LDVFE_LAPL_BC, YDCST, YDGEOMETRY,YDDYN,LDNHQE_C2,KLEV,KLON,KNLON,PIN,POUT,PMULFAC)

!**** *SILKAP* - Semi-implicit scheme: compute LLsstar_kap Z

!     Purpose.
!     --------
!           Semi-implicit scheme: compute LLsstar_kap Z
!           This routine is called in the linear system of NHQE model

!**   Interface.
!     ----------
!        *CALL* *SILKAP(...)

!        Explicit arguments :
!        --------------------
!        * INPUT:
!          YDGEOMETRY    : structure containing geometry
!          YDDYN         : structure containing dynamics
!          LDNHQE_C2     : T vs. F: constraint C2 is (vs. is not) matched
!          KLEV,KLON     : couple of dimensions; ex: (NPROMA,1) for grid-point array; (1,NFLEVG) for spectral array
!          KNLON         : number of vertical columns treated
!          PIN           : input array

!        * OUTPUT:
!          POUT          : output array

!        * OPTIONAL INPUT:
!          PMULFAC       : multiplication factor

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
!      J. Vivoda and P. Smolikova (Sep 2020): VFE pruning.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYN       , ONLY : TDYN
USE YOMCST       , ONLY : TCST

USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL, INTENT (IN) :: LDVFE_LAPL_BC
LOGICAL, INTENT (IN) :: LDVERTFE
TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)          ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)              ,INTENT(IN)    :: YDDYN
LOGICAL                 ,INTENT(IN)    :: LDNHQE_C2
INTEGER(KIND=JPIM)      ,INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM)      ,INTENT(IN)    :: KLON
INTEGER(KIND=JPIM)      ,INTENT(IN)    :: KNLON
REAL(KIND=JPRB)         ,INTENT(IN)    :: PIN(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB)         ,INTENT(OUT)   :: POUT(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PMULFAC

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZZ1(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) :: ZZ2(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) :: ZFACVAR(YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: IDT, JLEV, JLON
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "sinhqe_seve.intfb.h"
#include "sinhqe_vderi.intfb.h"
#include "siskap.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SILKAP',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

! * Preliminary computations:
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV)
DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
  ZFACVAR(JLEV)=YDDYN%SITR/YDDYN%SITRAM(JLEV)
ENDDO
!$OMP END PARALLEL DO

! * VFE treatment of LLsstar_kap, constraint C2 is not ensured
IF (LDVFE_LAPL_BC) CALL ABOR1(' SILKAP: case not yet coded!!!')

IF (LDNHQE_C2) THEN

  ! * VFD treatment of LLsstar_kap, constraint C2 is ensured
  !   LLsstar_kap Z = (Tstar/Tastar) Sstar_kap LLstar Z

  ! * Multiply by LLstar: ZZ1 = LLstar Z
  CALL SINHQE_SEVE(LDVFE_LAPL_BC, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,PIN,ZZ1,KNLON,LD_LSTAR=.TRUE.)

  ! * Multiply by Sstar_kap: ZZ2 = Sstar_kap ZZ1
  CALL SISKAP(LDVERTFE, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,KNLON,ZZ1,ZZ2)

  ! * Multiply by (Tstar/Tastar)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLON,JLEV,IDT)
  DO JLON=1,KNLON
    DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      POUT(IDT)=ZFACVAR(JLEV)*ZZ2(IDT)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

ELSEIF (.NOT.LDNHQE_C2) THEN
  
  ! * VFD treatment of LLsstar_kap, constraint C2 is not ensured
  !   LLsstar_kap Z = (Tstar/Tastar) [ LLstar Z + (Rd/Cpd - 1) partial_star Z ]

  ! * Multiply by LLstar: ZZ1 = LLstar Z
  CALL SINHQE_SEVE(LDVFE_LAPL_BC, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,PIN,ZZ1,KNLON,LD_LSTAR=.TRUE.)

  ! * Compute ZZ2 = partial_star Z
  CALL SINHQE_VDERI(YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,PIN,ZZ2,KNLON,LD_LSTAR=.TRUE.)

  ! * Final calculation
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLON,JLEV,IDT)
  DO JLON=1,KNLON
    DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      POUT(IDT)=ZFACVAR(JLEV)*(ZZ1(IDT)+(YDCST%RKAPPA-1.0_JPRB)*ZZ2(IDT))
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

ENDIF

! * Additional multiplication by a constant factor:
IF (PRESENT(PMULFAC)) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLON,JLEV,IDT)
  DO JLON=1,KNLON
    DO JLEV=1,YDGEOMETRY%YRDIMV%NFLEVG
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      POUT(IDT)=PMULFAC*POUT(IDT)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDIF

!      -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SILKAP',1,ZHOOK_HANDLE)
END SUBROUTINE SILKAP
