SUBROUTINE LATTE_KAPPA(YDDYNA, YDCST, YDGEOMETRY,YDVARS,YDCPG_SL1,YDCPG_SL2, YDCPG_BNDS, YDCPG_OPTS,YDDYN,PDTS2,PSLHDA,PSLHDD0,&
 & PUVH0,PRES0,PRDELP)

!**** *LATTE_KAPPA* input for SLHD diffusion.
!                   Computation of the kappa function ("coefficient
!                   of diffusion") based on rescaled horizontal
!                   deformation of the flow to be evaluated at F, and put
!                   it in PB2.
!                   Used by SLHD scheme.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LATTE_KAPPA(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST     - first element of work.
!          KPROF   - depth of work.
!          PDTS2   - 0.5*"pdt", where "pdt" =
!                    time step for the first time-integration step of
!                    a leap-frog scheme or all time-integration steps of
!                    a two-time level scheme; 2*time step for the following
!                    time-integration steps of a leap-frog scheme.
!          PSLHDA  - Scaling factor of the deformation in f(d) function
!                    (including the model resolution correction)
!          PSLHDD0 - Treshold for deformation tensor enhancement
!          PGMV    - GMV variables at t-dt and t.
!          PGMVS   - GMVS variables at t-dt and t.
!          PUVH0 - horizontal wind at time t at half levels.
!          PRES0   - prehyds at t.
!          PRDELP  - 1/(pressure depth of layers) at t.

!        INPUT/OUTPUT:
!          PB2     - "SLBUF2" buffer.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Called by LACDYN.

!     Reference.
!     ----------
!        Arpege documentation about semi-Lagrangian scheme.
!        Arpege documentation about horizontal diffusion.

!     Author.
!     -------
!        Original (F. VANA and K. YESSAD): AUGUST 2003.

!     Modifications.
!     --------------
!        Modified 03-08-19 by F. Vana : horizontal flow deformation for SLHD
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        F. Vana  04-10-25 : New content of output + DFI fix -> renamed
!        F. Vana  06-06-08 : SLHDA,SLHDD0 changed to arays P...
!        K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!         + put calculation of Kappa in GP_KAPPA.
!        K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!        F. Vana  13-feb-2014  kappaT for heat variables, special treatment above tropopause
!        F. Vana  Jul 2018 : SLHD acting like a sponge
!        F. Vana  21-May-2019 : Computation limited to levels where required.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE INTDYN_MOD   , ONLY : YYTHW0
USE YOMDYNA      , ONLY : TDYNA
USE YOMDYN       , ONLY : TDYN
USE YOMCST       , ONLY : TCST
USE FIELD_VARIABLES_MOD,ONLY : FIELD_VARIABLES
USE CPG_OPTS_TYPE_MOD, ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE CPG_TYPE_MOD,ONLY : CPG_SL1_TYPE, CPG_SL2_TYPE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TDYNA), INTENT (IN) :: YDDYNA
TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(CPG_SL1_TYPE)   ,INTENT(INOUT) :: YDCPG_SL1
TYPE(CPG_SL2_TYPE)   ,INTENT(INOUT) :: YDCPG_SL2
TYPE(CPG_BNDS_TYPE)   ,INTENT(IN)     :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE)   ,INTENT(IN)     :: YDCPG_OPTS
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDA(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDD0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUVH0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTHW0%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRES0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
!     -------------------------------------------------------
INTEGER(KIND=JPIM) :: JLEV, JROF
INTEGER(KIND=JPIM) :: ITROP(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZKAPPA(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZKAPPAT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     -------------------------------------------------------

#include "tropo_tep.intfb.h"
#include "gp_kappa.intfb.h"
#include "gp_kappat.intfb.h"

!     -------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LATTE_KAPPA',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDVAB=>YDGEOMETRY%YRVAB)
ASSOCIATE(NPROMA=>YDDIM%NPROMA,   NFLEVG=>YDDIMV%NFLEVG, NLEV_SPONGE=>YDDYN%NLEV_SPONGE,   LSLHDHEAT=>YDDYN%LSLHDHEAT,      &
& SLHD_MASK_U=>YDDYN%SLHD_MASK_U, SLHD_MASK_T=>YDDYN%SLHD_MASK_T)
!    --------------------------------------------------------


!*       1.   Kappa functions computation

! Momentum
CALL GP_KAPPA(YDDYNA, YDVAB,YDDYN,NPROMA,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,NFLEVG,PDTS2,&
 & YDVARS%U%DL,YDVARS%V%DL,YDVARS%VOR%T0,YDVARS%DIV%T0,&
 & PUVH0(1,0,YYTHW0%M_UH),PUVH0(1,0,YYTHW0%M_VH),PRES0,YDVARS%SP%DL,YDVARS%SP%DM,&
 & PRDELP,PSLHDA,PSLHDD0,ZKAPPA)

! Heat
IF (LSLHDHEAT) THEN
  CALL GP_KAPPAT(YDCPG_OPTS%LVERTFE, YDCPG_OPTS%TOPPRES, YDDYNA, YDCST, YDGEOMETRY,YDDYN,NPROMA,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,NFLEVG,PDTS2,&
   & YDVARS%T%T0,YDVARS%T%DL,YDVARS%T%DM,&
   & YDVARS%SP%T0,YDVARS%SP%DL,YDVARS%SP%DM,&
   & PSLHDA,PSLHDD0,ZKAPPAT)
!ELSE
! ZKAPPAT(KST:KPROF,1:NFLEVG)=ZKAPPA(KST:KPROF,1:NFLEVG)
ENDIF

! Specific treatment above tropopause
IF (YDDYNA%SLHDKMIN /= YDDYNA%SLHDKREF) THEN

  !  Detection of tropopause level
  CALL TROPO_TEP(YDCPG_OPTS%RPSTRA, YDCPG_OPTS%RPTROP, YDCPG_OPTS%LVERTFE, YDDYNA, YDCST,YDGEOMETRY,NPROMA,YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA,NFLEVG,&
   & YDVARS%SP%T0,YDVARS%T%T0,ITROP)

  !  Tag areas above troposphere
  DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    DO JLEV=1,MIN(ITROP(JROF),NLEV_SPONGE)
      ZKAPPA(JROF,JLEV)=-ZKAPPA(JROF,JLEV)
    ENDDO
  ENDDO
  IF (LSLHDHEAT) THEN
    DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      DO JLEV=1,MIN(ITROP(JROF),NLEV_SPONGE)
        ZKAPPAT(JROF,JLEV)=-ZKAPPAT(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF

ENDIF

!*       2.   Fill PB2 with Kappa + masking out areas of no use

DO JLEV=1,NLEV_SPONGE
  YDCPG_SL2%KAPPA(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=SLHD_MASK_U(JLEV)*ZKAPPA(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
  IF (LSLHDHEAT) YDCPG_SL2%KAPPAT(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)= &
   &   SLHD_MASK_T(JLEV)*ZKAPPAT(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)
ENDDO
DO JLEV=NLEV_SPONGE+1,NFLEVG
  YDCPG_SL2%KAPPA(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)=0._JPRB
  IF (LSLHDHEAT) YDCPG_SL2%KAPPAT(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA,JLEV)= 0._JPRB
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LATTE_KAPPA',1,ZHOOK_HANDLE)
END SUBROUTINE LATTE_KAPPA
