SUBROUTINE LANHQESI(YDCST,YDGEOMETRY,YDVARS,YDCPG_BNDS,YDCPG_OPTS,YDGMV,YDGFL,YDDYN,LDGWADV,PRCORI,PGMV,PGMVS,PGFL,&
 & PSDIV0,PSDIV9,PTOD0,PTOD9,PGAGT0L,PGAGT0M,PGAGT9L,PGAGT9M,PLPD0,PLPD9)

!**** *LANHQESI*   Semi-Lagrangian scheme.
!                  Computation of linear terms used in the semi-implicit scheme. 
!                  NHQE model.

!     Purpose.
!     --------
!        Computation of linear terms used in the semi-implicit scheme for NHQE model:
!        Nabla(Gamma*T+Mu*Pi), (Tau*D), (Nu*D), and some additional anhydrostatic contributions.

!**   Interface.
!     ----------
!        *CALL* *LANHQESI(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          YDGEOMETRY   : structure containing geometry.
!          YDGMV        : structure containing GMV.
!          YDGFL        : structure containing GFL.
!          YDDYN        : structure containing dynamics.
!          KST          : first element of work.
!          KPROF        : depth of work.
!          LDGWADV      : T/F: SVD var is 'gw' or vertical divergence
!          PRCORI       : Coriolis parameter 2.OMEGA .
!          PGMV         : GMV variables at t-dt and t.
!          PGMVS        : GMVS variables at t-dt and t.
!          PGFL         : unified_treatment grid-point (GFL) fields.

!        OUTPUT:
!          PSDIV0       : linear term at time t for continuity equation.
!          PSDIV9       : linear term at time t-dt for continuity equation.
!          PTOD0        : linear term at time t for temperature equation.
!          PTOD9        : linear term at time t-dt for temperature equation.
!          PGAGT0L      : linear term at time t for U-wind equation.
!          PGAGT0M      : linear term at time t for V-wind equation.
!          PGAGT9L      : linear term at time t-dt for U-wind equation.
!          PGAGT9M      : linear term at time t-dt for V-wind equation.
!          PLPD0        : linear term at time t for SVD-var equation.
!          PLPD9        : linear term at time t-dt for SVD-var equation.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation.

!     Externals.
!     ----------
!        Called by LACDYN.

!     Reference.
!     ----------
!        Arpege documentation about semi-Lagrangian scheme (IDSL).
!        Documentation (IDSI).

!     Author.
!     -------
!      K. YESSAD (METEO FRANCE/CNRM/GMAP).
!      Original : March 2017.

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOM_YGFL     , ONLY : TYPE_GFLD
USE YOMCST       , ONLY : TCST, RD, RG
USE YOMCT0       , ONLY : LSPRT, LTWOTL
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : LNHQE_SIHYD, LNHQE_C2
USE FIELD_VARIABLES_MOD,ONLY : FIELD_VARIABLES
USE CPG_OPTS_TYPE_MOD, ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCST),        INTENT(IN)    :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(CPG_BNDS_TYPE)   ,INTENT(IN)   :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE)   ,INTENT(IN)   :: YDCPG_OPTS
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YDGFL
 
 
LOGICAL           ,INTENT(IN)    :: LDGWADV
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCORI(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGFL%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSDIV0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSDIV9(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTOD0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTOD9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAGT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAGT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAGT9L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAGT9M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLPD0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLPD9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZR9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZF

INTEGER(KIND=JPIM) :: IPROFS, JLEV, JROF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gprcp_pgfl.intfb.h"
#include "sigam.intfb.h"
#include "sitnu.intfb.h"
#include "silkap.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LANHQESI',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA,NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & LIMPF=>YDDYN%LIMPF,SITR=>YDDYN%SITR, &
 & YT0=>YDGMV%YT0,YT9=>YDGMV%YT9)

!     ------------------------------------------------------------------

IPROFS=YDCPG_BNDS%KFDIA-YDCPG_BNDS%KIDIA+1
ZF=RG*RG/(RD*SITR)

! * Variables at time 0 (SL3TL and SL2TL).

IF (LDGWADV) THEN

  DO JLEV=1,NFLEVG
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      PGAGT0L(JROF,JLEV)=0.0_JPRB
      PGAGT0M(JROF,JLEV)=0.0_JPRB
      PTOD0  (JROF,JLEV)=0.0_JPRB
      PLPD0  (JROF,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
  PSDIV0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0.0_JPRB

ELSEIF (.NOT.LDGWADV) THEN

  ! - Computation of Nu*D (SI term for continuity equation) and Tau*D (SI term for temperature equation).
  CALL SITNU(YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGMV(YDCPG_BNDS%KIDIA,1,YT0%MDIV),PTOD0(YDCPG_BNDS%KIDIA,1),PSDIV0(YDCPG_BNDS%KIDIA),IPROFS)
  ! - Computation of Nabla(Gamma*T+Mu*Pi) (SI term for momentum equation).
  CALL SIGAM(YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGAGT0L(YDCPG_BNDS%KIDIA,1),PGMV(YDCPG_BNDS%KIDIA,1,YT0%MTL),PGMVS(YDCPG_BNDS%KIDIA,YT0%MSPL),IPROFS,NFLEVG)
  CALL SIGAM(YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGAGT0M(YDCPG_BNDS%KIDIA,1),PGMV(YDCPG_BNDS%KIDIA,1,YT0%MTM),PGMVS(YDCPG_BNDS%KIDIA,YT0%MSPM),IPROFS,NFLEVG)
  ! - Add NH contributions to (PGAGT0L,PGAGT0M).
  IF (.NOT.LNHQE_SIHYD) THEN
    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        PGAGT0L(JROF,JLEV)=PGAGT0L(JROF,JLEV)+RD*SITR*PGMV(JROF,JLEV,YT0%MSPDL)
        PGAGT0M(JROF,JLEV)=PGAGT0M(JROF,JLEV)+RD*SITR*PGMV(JROF,JLEV,YT0%MSPDM)
      ENDDO
    ENDDO
    ! - Computation of linear terms for vertical divergence variable.
    CALL SILKAP(YDCST,YDGEOMETRY,YDDYN,LNHQE_C2,NPROMA,1,IPROFS,PGMV(YDCPG_BNDS%KIDIA,1,YT0%MSPD),PLPD0(YDCPG_BNDS%KIDIA,1),PMULFAC=ZF)
  ELSE
    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        PLPD0  (JROF,JLEV)=0.0_JPRB
      ENDDO
    ENDDO
  ENDIF

ENDIF

!     ------------------------------------------------------------------

! * Variables at time 9 (SL3TL only).

IF (.NOT.LTWOTL) THEN

  IF (LDGWADV) THEN

    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        PGAGT9L(JROF,JLEV)=0.0_JPRB
        PGAGT9M(JROF,JLEV)=0.0_JPRB
        PTOD9  (JROF,JLEV)=0.0_JPRB
        PLPD9  (JROF,JLEV)=0.0_JPRB
      ENDDO
    ENDDO
    PSDIV9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0.0_JPRB

  ELSEIF (.NOT.LDGWADV) THEN

    ! - Computation of Nu*D (SI term for continuity equation) and Tau*D (SI term for temperature equation).
    CALL SITNU(YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGMV(YDCPG_BNDS%KIDIA,1,YT9%MDIV),PTOD9(YDCPG_BNDS%KIDIA,1),PSDIV9(YDCPG_BNDS%KIDIA),IPROFS)
    ! - Computation of Nabla(Gamma*T+Mu*Pi) (SI term for momentum equation).
    CALL SIGAM(YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGAGT9L(YDCPG_BNDS%KIDIA,1),PGMV(YDCPG_BNDS%KIDIA,1,YT9%MTL),PGMVS(YDCPG_BNDS%KIDIA,YT9%MSPL),IPROFS,NFLEVG)
    CALL SIGAM(YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGAGT9M(YDCPG_BNDS%KIDIA,1),PGMV(YDCPG_BNDS%KIDIA,1,YT9%MTM),PGMVS(YDCPG_BNDS%KIDIA,YT9%MSPM),IPROFS,NFLEVG)
    ! - Add NH contributions to (PGAGT0L,PGAGT0M).
    IF (.NOT.LNHQE_SIHYD) THEN
      DO JLEV=1,NFLEVG
        DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          PGAGT9L(JROF,JLEV)=PGAGT9L(JROF,JLEV)+RD*SITR*PGMV(JROF,JLEV,YT9%MSPDL)
          PGAGT9M(JROF,JLEV)=PGAGT9M(JROF,JLEV)+RD*SITR*PGMV(JROF,JLEV,YT9%MSPDM)
        ENDDO
      ENDDO
      ! - Computation of linear terms for vertical divergence variable.
      CALL SILKAP(YDCST,YDGEOMETRY,YDDYN,LNHQE_C2,NPROMA,1,IPROFS,PGMV(YDCPG_BNDS%KIDIA,1,YT9%MSPD),PLPD9(YDCPG_BNDS%KIDIA,1),PMULFAC=ZF)
    ELSE
      DO JLEV=1,NFLEVG
        DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
          PLPD9  (JROF,JLEV)=0.0_JPRB
        ENDDO
      ENDDO
    ENDIF

  ENDIF

ENDIF

!     ------------------------------------------------------------------

! * For "spectral RT" option, adjust semi-implicit term in
!   T-equation to compensate for later multiplication by R/Rd

IF (LSPRT) THEN
  IF (LTWOTL) THEN
    ! remark for lpc_full:
    !  predictor: treatment of "t" data.
    !  corrector: treatment of provisional "t+dt" data.
    ! So in this case this is always the Y[X]%MP data which are used.
    CALL GPRCP_PGFL(NPROMA,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,NFLEVG,PGFL=PGFL,PR=ZR9)
    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        PTOD0(JROF,JLEV)=RD*PTOD0(JROF,JLEV)/ZR9(JROF,JLEV)
      ENDDO
    ENDDO
  ELSE
    CALL GPRCP_PGFL(NPROMA,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,NFLEVG,PGFL=PGFL,KGFLTYP=9,PR=ZR9)
    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        PTOD0(JROF,JLEV)=RD*PTOD0(JROF,JLEV)/ZR9(JROF,JLEV)
        PTOD9(JROF,JLEV)=RD*PTOD9(JROF,JLEV)/ZR9(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
ENDIF

! * Add semi-implicit Coriolis terms to momentum eqn SI terms if required (LIMPF=.T.).

IF (LIMPF .AND. (.NOT.LDGWADV)) THEN
  DO JLEV=1,NFLEVG
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      PGAGT0L(JROF,JLEV)=PGAGT0L(JROF,JLEV)-PRCORI(JROF)*PGMV(JROF,JLEV,YT0%MV)
      PGAGT0M(JROF,JLEV)=PGAGT0M(JROF,JLEV)+PRCORI(JROF)*PGMV(JROF,JLEV,YT0%MU)
    ENDDO
  ENDDO
  IF (.NOT.LTWOTL) THEN
    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        PGAGT9L(JROF,JLEV)=PGAGT9L(JROF,JLEV)-PRCORI(JROF)*PGMV(JROF,JLEV,YT9%MV)
        PGAGT9M(JROF,JLEV)=PGAGT9M(JROF,JLEV)+PRCORI(JROF)*PGMV(JROF,JLEV,YT9%MU)
      ENDDO
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LANHQESI',1,ZHOOK_HANDLE)
END SUBROUTINE LANHQESI

