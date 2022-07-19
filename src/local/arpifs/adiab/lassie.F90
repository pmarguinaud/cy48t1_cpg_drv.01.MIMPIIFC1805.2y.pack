SUBROUTINE LASSIE(YDCST,YDGEOMETRY,YDVARS,YDCPG_BNDS,YDCPG_OPTS,YDDYN,PRCORI,&
 & PSDIV0,PSDIV9,PTOD0,PTOD9,PGAGT0L,PGAGT0M,PGAGT9L,PGAGT9M)  

!**** *LASSIE*   Semi-Lagrangian scheme.
!                Computation of linear terms used in the semi-implicit scheme. 

!     Purpose.
!     --------
!        Computation of linear terms used in the semi-implicit scheme:
!        Nabla(Gamma*T+Mu*Pi), (Tau*D) and (Nu*D).

!**   Interface.
!     ----------
!        *CALL* *LASSIE(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST      - first element of work.
!          KPROF    - depth of work.
!          PRCORI   - Coriolis parameter 2.OMEGA .
!          PGMV     - GMV variables at t-dt and t.
!          PGMVS    - GMVS variables at t-dt and t.
!          PGFL     - unified_treatment grid-point (GFL) fields.

!        OUTPUT:
!          PSDIV0   - semi-implicit term at time t for continuity equation
!                     (Nu*D).
!          PSDIV9   - semi-implicit term at time t-dt for continuity equation
!                     (Nu*D).
!          PTOD0    - semi-implicit term at time t for temperature equation
!                     (Tau*D).
!          PTOD9    - semi-implicit term at time t-dt for temperature equation
!                     (Tau*D).
!          PGAGT0L  - semi-implicit term at time t for U-wind equation
!                     (zonal component of Nabla(Gamma*T+Mu*Pi)).
!          PGAGT0M  - semi-implicit term at time t for V-wind equation
!                     (meridian component of Nabla(Gamma*T+Mu*Pi)).
!          PGAGT9L  - semi-implicit term at time t-dt for U-wind equation
!                     (zonal component of Nabla(Gamma*T+Mu*Pi)).
!          PGAGT9M  - semi-implicit term at time t-dt for V-wind equation
!                     (meridian component of Nabla(Gamma*T+Mu*Pi)).

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           Called by LACDYN.

!     Reference.
!     ----------
!             Arpege documentation about semi-lagrangian scheme.

!     Author.
!     -------
!      K. YESSAD (METEO FRANCE/CNRM/GMAP) after part 3.1 of LACDYN.
!      Loops have been recoded according to F90 norms.
!      Original : AUGUST 1995.

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Modified 04-02-06: Y. Seity : new arguments to GPRCP (rain, snow and graupel)
!      M.Hamrud  15-Jan-2006  Revised GPRCP
!      K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : TCST
USE YOMCT0       , ONLY : LTWOTL, LSPRT
USE YOMDYN       , ONLY : TDYN
USE FIELD_VARIABLES_MOD,ONLY : FIELD_VARIABLES
USE CPG_OPTS_TYPE_MOD, ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCST),        INTENT(IN)    :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(CPG_BNDS_TYPE)   ,INTENT(IN)   :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE)   ,INTENT(IN)   :: YDCPG_OPTS
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCORI(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSDIV0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSDIV9(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTOD0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTOD9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAGT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAGT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAGT9L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAGT9M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZR9  (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IPROFS, JLEV, JROF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gprcp_expl.intfb.h"
#include "sigam.intfb.h"
#include "sitnu.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LASSIE',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV)
ASSOCIATE(NPROMA=>YDDIM%NPROMA,   NFLEVG=>YDDIMV%NFLEVG,   LIMPF=>YDDYN%LIMPF  )
!     ------------------------------------------------------------------

IPROFS=YDCPG_BNDS%KFDIA-YDCPG_BNDS%KIDIA+1

! * Variables at time 0 (SL3TL and SL2TL).

!   - Computation of Nu*D (SI term for continuity equation)
!     and Tau*D (SI term for temperature equation).
CALL SITNU(YDCPG_OPTS%LVERTFE, YDCST,YDGEOMETRY,YDDYN,NPROMA,1,YDVARS%DIV%T0,PTOD0(YDCPG_BNDS%KIDIA,1),PSDIV0(YDCPG_BNDS%KIDIA),IPROFS)
!   - Computation of Nabla(Gamma*T+Mu*Pi) (SI term for momentum equation).
CALL SIGAM(YDCPG_OPTS%LVERTFE, YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGAGT0L(YDCPG_BNDS%KIDIA,1),YDVARS%T%DL,YDVARS%SP%DL,IPROFS,NFLEVG)
CALL SIGAM(YDCPG_OPTS%LVERTFE, YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGAGT0M(YDCPG_BNDS%KIDIA,1),YDVARS%T%DM,YDVARS%SP%DM,IPROFS,NFLEVG)

! * Variables at time 9 (SL3TL only).

IF (.NOT.LTWOTL) THEN
  ! - Computation of Nu*D (SI term for continuity equation)
  !   and Tau*D (SI term for temperature equation).
  CALL SITNU(YDCPG_OPTS%LVERTFE, YDCST,YDGEOMETRY,YDDYN,NPROMA,1,YDVARS%DIV%T9,PTOD9(YDCPG_BNDS%KIDIA,1),PSDIV9(YDCPG_BNDS%KIDIA),IPROFS)
  ! - Computation of Nabla(Gamma*T+Mu*Pi) (SI term for momentum equation).
  CALL SIGAM(YDCPG_OPTS%LVERTFE, YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGAGT9L(YDCPG_BNDS%KIDIA,1),YDVARS%T%DL9,YDVARS%SP%DL9,IPROFS,NFLEVG)
  CALL SIGAM(YDCPG_OPTS%LVERTFE, YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGAGT9M(YDCPG_BNDS%KIDIA,1),YDVARS%T%DM9,YDVARS%SP%DM9,IPROFS,NFLEVG)
ENDIF

! * For "spectral RT" option, adjust semi-implicit term (Tau*D) in
!   T-equation to compensate for later multiplication by R/Rd

IF (LSPRT) THEN
  IF (LTWOTL) THEN
    ! remark for lpc_full:
    !  predictor: treatment of "t" data.
    !  corrector: treatment of provisional "t+dt" data.
    ! So in this case this is always the Y[X]%MP data which are used.
    CALL GPRCP_EXPL(YDCST,NPROMA,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,NFLEVG,YDVARS=YDVARS,PR=ZR9)  
    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        PTOD0(JROF,JLEV)=YDCST%RD*PTOD0(JROF,JLEV)/ZR9(JROF,JLEV)
      ENDDO
    ENDDO
  ELSE
    CALL GPRCP_EXPL(YDCST,NPROMA,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,NFLEVG,YDVARS=YDVARS,PR=ZR9,KGFLTYP=9)  
    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        PTOD0(JROF,JLEV)=YDCST%RD*PTOD0(JROF,JLEV)/ZR9(JROF,JLEV)
        PTOD9(JROF,JLEV)=YDCST%RD*PTOD9(JROF,JLEV)/ZR9(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
ENDIF

! * Add semi-implicit Coriolis terms to Nabla(Gamma*T+Mu*Pi)
!   if required (LIMPF=.T.).

IF (LIMPF) THEN
  DO JLEV=1,NFLEVG
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      PGAGT0L(JROF,JLEV)=PGAGT0L(JROF,JLEV)-PRCORI(JROF)*YDVARS%V%T0(JROF,JLEV)
      PGAGT0M(JROF,JLEV)=PGAGT0M(JROF,JLEV)+PRCORI(JROF)*YDVARS%U%T0(JROF,JLEV)
    ENDDO
  ENDDO
  IF (.NOT.LTWOTL) THEN
    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        PGAGT9L(JROF,JLEV)=PGAGT9L(JROF,JLEV)&
         & -PRCORI(JROF)*YDVARS%V%T9(JROF,JLEV)
        PGAGT9M(JROF,JLEV)=PGAGT9M(JROF,JLEV)&
         & +PRCORI(JROF)*YDVARS%U%T9(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LASSIE',1,ZHOOK_HANDLE)
END SUBROUTINE LASSIE

