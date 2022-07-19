SUBROUTINE SILKAPI(YDCST, YDGEOMETRY,YDDYN,KOPT,KLEV,KLON,KNLON,PIN,POUT)

!**** *SILKAPI* - Semi-implicit scheme: compute LLsstar_kap**-1 Sstar_kap Z or LLsstar_kap**-1 Z

!     Purpose.
!     --------
!           Semi-implicit scheme:
!           compute LLsstar_kap**-1 Sstar_kap Z or LLsstar_kap**-1 Z
!           according to value of KOPT.
!           This routine is called in the linear system of NHQE model

!**   Interface.
!     ----------
!        *CALL* *SILKAPI(...)

!        Explicit arguments :
!        --------------------
!        * INPUT:
!          YDGEOMETRY    : structure containing geometry
!          YDDYN         : structure containing dynamics
!          KOPT          : 0: compute LLstar**-1 Z
!                          1: compute LLsstar_kap**-1 Z
!                          2: compute [LLsstar_kap**-1 Sstar_kap] Z
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
!      J. Vivoda and P. Smolikova (Sep 2020): VFE pruning.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : YRDYNA
USE YOMCVER      , ONLY : LVFE_LAPL_BC, LVERTFE
USE YOMCST       , ONLY : TCST
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KOPT
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KNLON
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIN(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POUT(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZZ0(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) :: ZZ1(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) :: ZZ2(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) :: ZZ3(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) :: ZZ4(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) :: ZZ5(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON)
REAL(KIND=JPRB) :: ZSP(KNLON)
REAL(KIND=JPRB) :: ZFACVAR(YDGEOMETRY%YRDIMV%NFLEVG),ZFACCST,ZCOEF
REAL(KIND=JPRB) :: ZID(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZL1(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZL2(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_LAPKSSI_TRID(1,YDGEOMETRY%YRDIMV%NFLEVG,-1:1)
REAL(KIND=JPRB) :: Z_ILAPKSSI(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IDT, JLEV, JLON
LOGICAL :: LLCSTSITRA, LL_TRIDIA
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "siskapi.intfb.h"
#include "siskap.intfb.h"
#include "siqq.intfb.h"
#include "sigam.intfb.h"
#include "sitnu.intfb.h"
#include "sinhqe_seve.intfb.h"
#include "sinhqe_vderi.intfb.h"
#include "minv_caller.h"
#include "tridia.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SILKAPI',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)
!     ------------------------------------------------------------------

! * VFE treatment of LLsstar_kap**-1, constraint C2 is not ensured
IF (LVFE_LAPL_BC) CALL ABOR1(' SILKAPI: case LVFE_LAPL_BC not yet coded!!!')

! * Preliminary computations:
!!! could be rewritten LLCSTSITRA=YDDYN%NOPT_SITRA == 0
LLCSTSITRA=MAXVAL(YDDYN%SITRAM(1:NFLEVG)) == MINVAL(YDDYN%SITRAM(1:NFLEVG))

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV)
DO JLEV=1,NFLEVG
  ZFACVAR(JLEV)=YDDYN%SITRAM(JLEV)/YDDYN%SITR
ENDDO
!$OMP END PARALLEL DO

ZCOEF=YDCST%RCPD/(YDCST%RD*YDCST%RD*YDDYN%SITR)

IF (YRDYNA%LNHQE_C2) THEN

  ! * VFD treatment of LLsstar_kap**-1, constraint C2 is ensured
  !   LLsstar_kap**-1 Z = - ( Cpd/(Rd*Rd*Tstar) (gamma tau + Rd Tstar nu) + QQstar ) Sstar_kap**-1 (Tastar/Tstar) Z

  IF (KOPT == 0) THEN

    ! * Multiply by QQstar: ZZ3 = QQstar Z
    CALL SIQQ(YDGEOMETRY,YDDYN,KLEV,KLON,KNLON,PIN,ZZ3)

    ! * Compute "ZZ4=tau Z" and "ZSP=nu Z"
    CALL SITNU(LVERTFE, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,PIN,ZZ4,ZSP,KNLON)

    ! * Compute "ZZ5 = gamma ZZ4 + Rd Tstar ZSP": ZZ5 = [gamma tau + Rd Tstar nu] Z
    CALL SIGAM(LVERTFE, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,ZZ5,ZZ4,ZSP,KNLON,NFLEVG)

    ! * Compute [- ( Cpd/(Rd*Rd*Tstar) ZZ5 + ZZ3 )]
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLON,JLEV,IDT)
    DO JLON=1,KNLON
      DO JLEV=1,NFLEVG
        IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
        POUT(IDT)=-ZCOEF*ZZ5(IDT)-ZZ3(IDT)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO

  ELSEIF (KOPT == 1) THEN

    ! * Multiply by (Tastar/Tstar)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLON,JLEV,IDT)
    DO JLON=1,KNLON
      DO JLEV=1,NFLEVG
        IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
        ZZ1(IDT)=ZFACVAR(JLEV)*PIN(IDT)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! * Multiply by Sstar_kap**-1: ZZ2 = Sstar_kap**-1 (Tastar/Tstar) Z
    CALL SISKAPI(YDGEOMETRY,YDDYN,KLEV,KLON,KNLON,ZZ1,ZZ2)

    ! * Multiply by QQstar: ZZ3 = QQstar Sstar_kap**-1 (Tastar/Tstar) Z
    CALL SIQQ(YDGEOMETRY,YDDYN,KLEV,KLON,KNLON,ZZ2,ZZ3)

    ! * Compute "ZZ4=tau ZZ2" and "ZSP=nu ZZ2"
    CALL SITNU(LVERTFE, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,ZZ2,ZZ4,ZSP,KNLON)

    ! * Compute "ZZ5 = gamma ZZ4 + Rd Tstar ZSP": ZZ5 = [gamma tau + Rd Tstar nu] Sstar_kap**-1 (Tastar/Tstar) Z
    CALL SIGAM(LVERTFE, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,ZZ5,ZZ4,ZSP,KNLON,NFLEVG)

    ! * Compute [- ( Cpd/(Rd*Rd*Tstar) ZZ5 + ZZ3 )]
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLON,JLEV,IDT)
    DO JLON=1,KNLON
      DO JLEV=1,NFLEVG
        IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
        POUT(IDT)=-ZCOEF*ZZ5(IDT)-ZZ3(IDT)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO

  ELSEIF (KOPT == 2 .AND. .NOT.LLCSTSITRA) THEN

    ! * Multiply by Sstar_kap: ZZ0 = Sstar_kap Z
    CALL SISKAP(YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,KNLON,PIN,ZZ0)

    ! * Multiply by (Tastar/Tstar)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLON,JLEV,IDT)
    DO JLON=1,KNLON
      DO JLEV=1,NFLEVG
        IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
        ZZ1(IDT)=ZFACVAR(JLEV)*ZZ0(IDT)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! * Multiply by Sstar_kap**-1: ZZ2 = Sstar_kap**-1 (Tastar/Tstar) ZZ0
    CALL SISKAPI(YDGEOMETRY,YDDYN,KLEV,KLON,KNLON,ZZ1,ZZ2)

    ! * Multiply by QQstar: ZZ3 = QQstar Sstar_kap**-1 (Tastar/Tstar) ZZ0
    CALL SIQQ(YDGEOMETRY,YDDYN,KLEV,KLON,KNLON,ZZ2,ZZ3)

    ! * Compute "ZZ4=tau ZZ2" and "ZSP=nu ZZ2"
    CALL SITNU(LVERTFE, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,ZZ2,ZZ4,ZSP,KNLON)

    ! * Compute "ZZ5 = gamma ZZ4 + Rd Tstar ZSP": ZZ5 = [gamma tau + Rd Tstar nu] Sstar_kap**-1 (Tastar/Tstar) ZZ0
    CALL SIGAM(LVERTFE, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,ZZ5,ZZ4,ZSP,KNLON,NFLEVG)

    ! * Compute [- ( Cpd/(Rd*Rd*Tstar) ZZ5 + ZZ3 )]
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLON,JLEV,IDT)
    DO JLON=1,KNLON
      DO JLEV=1,NFLEVG
        IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
        POUT(IDT)=-ZCOEF*ZZ5(IDT)-ZZ3(IDT)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO

  ELSEIF (KOPT == 2 .AND. LLCSTSITRA) THEN
    ! * Simplification: Sstar_kap**-1 (Tastar/Tstar) Sstar_kap = (Tastar/Tstar)

    ZFACCST=ZFACVAR(1)

    ! * Multiply by (Tastar/Tstar)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLON,JLEV,IDT)
    DO JLON=1,KNLON
      DO JLEV=1,NFLEVG
        IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
        ZZ2(IDT)=ZFACCST*PIN(IDT)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! * Multiply by QQstar: ZZ3 = QQstar (Tastar/Tstar) Z
    CALL SIQQ(YDGEOMETRY,YDDYN,KLEV,KLON,KNLON,ZZ2,ZZ3)

    ! * Compute "ZZ4=tau ZZ2" and "ZSP=nu ZZ2"
    CALL SITNU(LVERTFE, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,ZZ2,ZZ4,ZSP,KNLON)

    ! * Compute "ZZ5 = gamma ZZ4 + Rd Tstar ZSP": ZZ5 = [gamma tau + Rd Tstar nu] (Tastar/Tstar) Z
    CALL SIGAM(LVERTFE, YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,ZZ5,ZZ4,ZSP,KNLON,NFLEVG)

    ! * Compute [- ( Cpd/(Rd*Rd*Tstar) ZZ5 + ZZ3 )]
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLON,JLEV,IDT)
    DO JLON=1,KNLON
      DO JLEV=1,NFLEVG
        IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
        POUT(IDT)=-ZCOEF*ZZ5(IDT)-ZZ3(IDT)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO

  ENDIF

ELSE

  ! * VFD treatment of LLsstar_kap**-1, constraint C2 is not ensured
  !   LLsstar_kap is a tridiagonal matrix.

  IF (KOPT == 0) CALL ABOR1(' SILKAPI: case LNHQE_C2=F with KOPT=0 not coded!!!')

  ! * ZID contains the identity matrix "I":
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLON,JLEV)
  DO JLEV=1,NFLEVG
    DO JLON=1,NFLEVG
      ZID(JLON,JLEV)=1.0_JPRB-MIN(ABS(REAL(JLON-JLEV,JPRB)),1.0_JPRB)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

  ! * Multiply by LLstar: ZL1 = LLstar ZID
  CALL SINHQE_SEVE(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZID,ZL1,NFLEVG,LD_LSTAR=.TRUE.)

  ! * Compute ZL2 = partial_star ZID
  CALL SINHQE_VDERI(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZID,ZL2,NFLEVG,LD_LSTAR=.TRUE.)

  ! * Put LLsstar_kap ZID = (Tstar/Tastar) [LLstar + (RKAPPA-1) partial_star] ZID in ZL1
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLON,JLEV)
  DO JLEV=1,NFLEVG
    DO JLON=1,NFLEVG
      ZL1(JLON,JLEV)=(ZL1(JLON,JLEV)+(YDCST%RKAPPA-1.0_JPRB)*ZL2(JLON,JLEV))/ZFACVAR(JLON)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

  LL_TRIDIA=.FALSE.

  IF(LL_TRIDIA) THEN

    !! ky: still remains bugs.

    ! upper diagonal
    DO JLEV=2,NFLEVG
      ! ???? is it (1,JLEV-1,1) or (1,JLEV-1,-1)?
      Z_LAPKSSI_TRID(1,JLEV-1,1)=ZL1(JLEV-1,JLEV)
    ENDDO
    ! central diagonal
    DO JLEV=1,NFLEVG
      Z_LAPKSSI_TRID(1,JLEV,0)=ZL1(JLEV,JLEV)
    ENDDO
    ! lower diagonal
    DO JLEV=1,NFLEVG-1
      ! ???? is it (1,JLEV+1,-1) or (1,JLEV+1,1)?
      Z_LAPKSSI_TRID(1,JLEV+1,-1)=ZL1(JLEV+1,JLEV)
    ENDDO

    ! * Invert LLsstar_kap:
    !   (!!! there is maybe a better algorithm to invert this tridiagonal matrix)
    CALL TRIDIA(NFLEVG,NFLEVG,1,NFLEVG,1,Z_LAPKSSI_TRID,ZID,Z_ILAPKSSI)

  ELSE

    !! ky: not optimal way to invert a tridiagonal matrix, but it works!

    ! * Invert LLsstar_kap:
    CALL MINV_CALLER(.TRUE.,NFLEVG,ZL1,Z_ILAPKSSI)

  ENDIF

  IF (KOPT == 1) THEN
    ! * POUT = LLsstar_kap**-1 PIN
    CALL MXMAOP(Z_ILAPKSSI,1,NFLEVG,PIN,1,NFLEVG,POUT,1,NFLEVG,NFLEVG,NFLEVG,KNLON)
  ELSEIF (KOPT == 2) THEN
    ! * POUT = LLsstar_kap**-1 Sstar_kap PIN
    CALL SISKAP(YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,KNLON,PIN,ZZ0)
    CALL MXMAOP(Z_ILAPKSSI,1,NFLEVG,ZZ0,1,NFLEVG,POUT,1,NFLEVG,NFLEVG,NFLEVG,KNLON)
  ENDIF
ENDIF

!      -----------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SILKAPI',1,ZHOOK_HANDLE)
END SUBROUTINE SILKAPI
