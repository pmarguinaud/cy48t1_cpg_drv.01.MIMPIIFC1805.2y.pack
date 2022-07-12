SUBROUTINE LANHSI(YDCST,YDGEOMETRY,YDVARS,YDCPG_BNDS,YDCPG_OPTS,YDGMV,YDGFL,YDDYN,LDGWADV,PRCORI,PREDIV,PGMV,PGMVS,PGFL,&
 & PSDIV0,PSDIV9,PTOD0,PTOD9,PGAGT0L,PGAGT0M,PGAGT9L,PGAGT9M,PSPDS0,PSPDS9,PLPD0,PLPD9)

!**** *LANHSI*   Semi-Lagrangian scheme.
!                Computation of linear terms used in the semi-implicit scheme. 
!                NHEE model.

!     Purpose.
!     --------
!        Computation of linear terms used in the semi-implicit scheme for NHEE model:
!        for U,V (momentum) equation, T (temperature) equation,
!        ps or ln(ps) (continuity) equation, NH variables equations.
!        NH model with (log(pre/prehyd),d3 or d4) NH variables. 

!**   Interface.
!     ----------
!        *CALL* *LANHSI(..)

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
!          PREDIV       : (c**2)/(M**2).
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
!          PSPDS0       : linear term at time t for SPD-var equation.
!          PSPDS9       : linear term at time t-dt for SPD-var equation.
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
!        Arpege documentation about semi-implicit scheme (IDSI).

!     Author.
!     -------
!      K. YESSAD (METEO FRANCE/CNRM/GMAP) after routine LASSIE and
!      code introduced by M. JANOUSEK and R. BUBNOVA.
!      Original : APRIL 1996.

!     Modifications.
!     --------------
!   K. Yessad 07-03-2007: Remove useless (gw)_surf interpol. in NH+LGWADV.
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!   K. Yessad (June 2017): Vertical-dependent SITRA.
!  End Modifications
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOM_YGFL     , ONLY : TYPE_GFLD
USE YOMCST       , ONLY : TCST
USE YOMCT0       , ONLY : LSPRT, LTWOTL
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
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YDGFL
 
 
LOGICAL           ,INTENT(IN)    :: LDGWADV
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCORI(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREDIV(YDGEOMETRY%YRDIM%NPROMA) 
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
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDS0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDS9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLPD0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLPD9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZR9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSVDT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDUM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IPROFS, JLEV, JROF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gprcp_pgfl.intfb.h"
#include "siseve.intfb.h"
#include "sidd.intfb.h"
#include "siptp.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LANHSI',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA,NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & LIMPF=>YDDYN%LIMPF,SITR=>YDDYN%SITR, &
 & YT0=>YDGMV%YT0,YT9=>YDGMV%YT9)

!     ------------------------------------------------------------------

IPROFS=YDCPG_BNDS%KFDIA-YDCPG_BNDS%KIDIA+1

!     ------------------------------------------------------------------

! * Variables at time 0 (SL3TL and SL2TL).

IF (LDGWADV) THEN

  DO JLEV=1,NFLEVG
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      PGAGT0L(JROF,JLEV)=0.0_JPRB
      PGAGT0M(JROF,JLEV)=0.0_JPRB
      PTOD0  (JROF,JLEV)=0.0_JPRB
      PSPDS0 (JROF,JLEV)=0.0_JPRB
      PLPD0  (JROF,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
  PSDIV0(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0.0_JPRB

ELSEIF (.NOT.LDGWADV) THEN

  ! - Computation of SI terms for continuity, T and Pcha equations.
  DO JLEV=1,NFLEVG
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZSVDT(JROF,JLEV)=PGMV(JROF,JLEV,YT0%MSVD)/PREDIV(JROF)
    ENDDO
  ENDDO
  CALL SIPTP(YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGMV(YDCPG_BNDS%KIDIA,1,YT0%MDIV),ZSVDT(YDCPG_BNDS%KIDIA,1),&
   & PSPDS0(YDCPG_BNDS%KIDIA,1),PTOD0(YDCPG_BNDS%KIDIA,1),PSDIV0(YDCPG_BNDS%KIDIA),IPROFS)
  ! - Computation of SI terms for momentum equation.
  CALL SIDD(YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGAGT0L(YDCPG_BNDS%KIDIA,1),ZDUM(YDCPG_BNDS%KIDIA,1),&
   & PGMV(YDCPG_BNDS%KIDIA,1,YT0%MSPDL),PGMV(YDCPG_BNDS%KIDIA,1,YT0%MTL),PGMVS(YDCPG_BNDS%KIDIA,YT0%MSPL),IPROFS)
  CALL SIDD(YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGAGT0M(YDCPG_BNDS%KIDIA,1),ZDUM(YDCPG_BNDS%KIDIA,1),&
   & PGMV(YDCPG_BNDS%KIDIA,1,YT0%MSPDM),PGMV(YDCPG_BNDS%KIDIA,1,YT0%MTM),PGMVS(YDCPG_BNDS%KIDIA,YT0%MSPM),IPROFS)
  ! - Computation of SI terms for rescaled "pseudo" vertical divergence eqn.
  CALL SISEVE(YDGEOMETRY,YDDYN,NPROMA,1,PGMV(YDCPG_BNDS%KIDIA,1,YT0%MSPD),PLPD0(YDCPG_BNDS%KIDIA,1),IPROFS)
  DO JLEV=1,NFLEVG
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      PLPD0(JROF,JLEV)=(YDCST%RG*YDCST%RG/(YDCST%RD*SITR))*PLPD0(JROF,JLEV)
    ENDDO
  ENDDO

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
        PSPDS9 (JROF,JLEV)=0.0_JPRB
        PLPD9  (JROF,JLEV)=0.0_JPRB
      ENDDO
    ENDDO
    PSDIV9(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=0.0_JPRB

  ELSEIF (.NOT.LDGWADV) THEN

    ! - Computation of SI terms for continuity, T and Pcha equations.
    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        ZSVDT(JROF,JLEV)=PGMV(JROF,JLEV,YT9%MSVD)/PREDIV(JROF)
      ENDDO
    ENDDO
    CALL SIPTP(YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGMV(YDCPG_BNDS%KIDIA,1,YT9%MDIV),ZSVDT(YDCPG_BNDS%KIDIA,1),&
     & PSPDS9(YDCPG_BNDS%KIDIA,1),PTOD9(YDCPG_BNDS%KIDIA,1),PSDIV9(YDCPG_BNDS%KIDIA),IPROFS)
    ! - Computation of SI terms for momentum equation.
    CALL SIDD(YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGAGT9L(YDCPG_BNDS%KIDIA,1),ZDUM(YDCPG_BNDS%KIDIA,1),&
     & PGMV(YDCPG_BNDS%KIDIA,1,YT9%MSPDL),PGMV(YDCPG_BNDS%KIDIA,1,YT9%MTL),PGMVS(YDCPG_BNDS%KIDIA,YT9%MSPL),IPROFS)
    CALL SIDD(YDCST,YDGEOMETRY,YDDYN,NPROMA,1,PGAGT9M(YDCPG_BNDS%KIDIA,1),ZDUM(YDCPG_BNDS%KIDIA,1),&
     & PGMV(YDCPG_BNDS%KIDIA,1,YT9%MSPDM),PGMV(YDCPG_BNDS%KIDIA,1,YT9%MTM),PGMVS(YDCPG_BNDS%KIDIA,YT9%MSPM),IPROFS)
    ! - Computation of SI terms for rescaled "pseudo" vertical divergence eqn.
    CALL SISEVE(YDGEOMETRY,YDDYN,NPROMA,1,PGMV(YDCPG_BNDS%KIDIA,1,YT9%MSPD),PLPD9(YDCPG_BNDS%KIDIA,1),IPROFS)
    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        PLPD9(JROF,JLEV)=(YDCST%RG*YDCST%RG/(YDCST%RD*SITR))*PLPD9(JROF,JLEV)
      ENDDO
    ENDDO

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
        PTOD0(JROF,JLEV)=YDCST%RD*PTOD0(JROF,JLEV)/ZR9(JROF,JLEV)
      ENDDO
    ENDDO
  ELSE
    CALL GPRCP_PGFL(NPROMA,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,NFLEVG,PGFL=PGFL,KGFLTYP=9,PR=ZR9)
    DO JLEV=1,NFLEVG
      DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
        PTOD0(JROF,JLEV)=YDCST%RD*PTOD0(JROF,JLEV)/ZR9(JROF,JLEV)
        PTOD9(JROF,JLEV)=YDCST%RD*PTOD9(JROF,JLEV)/ZR9(JROF,JLEV)
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
IF (LHOOK) CALL DR_HOOK('LANHSI',1,ZHOOK_HANDLE)
END SUBROUTINE LANHSI

