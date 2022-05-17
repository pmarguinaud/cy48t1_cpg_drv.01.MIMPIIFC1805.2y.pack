SUBROUTINE LATTE_KAPPAAD(YDGEOMETRY,YDGMV,YDGMV5,YDDYN,YDPTRSLB2,KST,KPROF,PDTS2,PSLHDA,PSLHDD0,&
 & PGMV,PGMVS,PGMV5,PGMV5S,PB2)

!**** *LATTE_KAPPAAD* input for SLHD diffusion. (adjoint model version)
!                   Computation of the kappa function ("coefficient
!                   of diffusion") based on rescaled horizontal
!                   deformation of the flow to be evaluated at F, and put
!                   it in PB2.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LATTE_KAPPAAD(..)

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
!          PGMV5   - GMV variables (trajectory) at t-dt and t.
!          PGMV5S  - GMVS variables (trajectory) at t-dt and t.

!        INPUT/OUTPUT:
!          PGMV    - GMV variables at t-dt and t.
!          PGMVS   - GMVS variables at t-dt and t.
!          PB2     - "SLBUF2" buffer

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Called by LACDYNAD.

!     Reference.
!     ----------

!     Author.
!     -------
!        Original F. VANA : January 2009

!     Modifications.
!     --------------
!     F. Vana  13-Feb-2014  kappaT for heat variables
!     F. Vana   Jul 2018    SLHD acting like a sponge
!     F. Vana  21-May-2019 : Computation limited to levels where required.

!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM     ,JPRB
USE YOMHOOK      , ONLY : LHOOK    ,DR_HOOK

USE YOMDYN       , ONLY : TDYN
USE PTRSLB2      , ONLY : TPTRSLB2

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV5
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TPTRSLB2)    ,INTENT(IN)    :: YDPTRSLB2
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDA(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDD0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMV5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV5%YT5%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMV5S(YDGEOMETRY%YRDIM%NPROMA,YDGMV5%YT5%NDIMS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDPTRSLB2%NFLDSLB2)
!     -------------------------------------------------------
INTEGER(KIND=JPIM) :: JLEV

REAL(KIND=JPRB) :: ZKAPPA(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB) :: ZKAPPAT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     -------------------------------------------------------

#include "gp_kappaad.intfb.h"
#include "gp_kappatad.intfb.h"

!     -------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LATTE_KAPPAAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV)
ASSOCIATE(NPROMA=>YDDIM%NPROMA,   NFLEVG=>YDDIMV%NFLEVG, NLEV_SPONGE=>YDDYN%NLEV_SPONGE,   LSLHDHEAT=>YDDYN%LSLHDHEAT,   &
& SLHD_MASK_U=>YDDYN%SLHD_MASK_U, SLHD_MASK_T=>YDDYN%SLHD_MASK_T,   YT0=>YDGMV%YT0,   YT5=>YDGMV5%YT5,                   &
& MSLB2KAPPA=>YDPTRSLB2%MSLB2KAPPA, MSLB2KAPPAT=>YDPTRSLB2%MSLB2KAPPAT)
!    --------------------------------------------------------

!*       2.   Get kappa from PB2 + masking out areas of no use

DO JLEV=1,NLEV_SPONGE
  IF (LSLHDHEAT) ZKAPPAT(KST:KPROF,JLEV)= &
   & SLHD_MASK_T(JLEV) * PB2(KST:KPROF,MSLB2KAPPAT+JLEV-1)
  ZKAPPA(KST:KPROF,JLEV) = SLHD_MASK_U(JLEV) * PB2(KST:KPROF,MSLB2KAPPA+JLEV-1)
ENDDO

!*       1.   Incrementing GMVs according kappa function

!   Heat
IF (LSLHDHEAT) THEN
  CALL GP_KAPPATAD(YDGEOMETRY,YDDYN,NPROMA,KST,KPROF,NFLEVG,PDTS2,&
   & PGMV5(1,1,YT5%MT),PGMV5(1,1,YT5%MTL),PGMV5(1,1,YT5%MTM),&
   & PGMV5S(1,YT5%MSP),PGMV5S(1,YT5%MSPL),PGMV5S(1,YT5%MSPM),&
   & PGMV(1,1,YT0%MT),PGMV(1,1,YT0%MTL),PGMV(1,1,YT0%MTM),&
   & PGMVS(1,YT0%MSP),PGMVS(1,YT0%MSPL),PGMVS(1,YT0%MSPM),&
   & PSLHDA,PSLHDD0,ZKAPPAT)
!ELSE
!  ZKAPPA (KST:KPROF,1:NFLEVG)=ZKAPPA(KST:KPROF,1:NFLEVG) &
!   & + ZKAPPAT(KST:KPROF,1:NFLEVG)
!  ZKAPPAT(KST:KPROF,1:NFLEVG)=0.0_JPRB
ENDIF

!   Momentum
CALL GP_KAPPAAD(YDDYN,NPROMA,KST,KPROF,NFLEVG,PDTS2,&
 & PGMV5(1,1,YT5%MUL),PGMV5(1,1,YT5%MVL),&
 & PGMV5(1,1,YT5%MVOR),PGMV5(1,1,YT5%MDIV),&
 & PGMV(1,1,YT0%MUL),PGMV(1,1,YT0%MVL),&
 & PGMV(1,1,YT0%MVOR),PGMV(1,1,YT0%MDIV),&
 & PSLHDA,PSLHDD0,ZKAPPA)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LATTE_KAPPAAD',1,ZHOOK_HANDLE)
END SUBROUTINE LATTE_KAPPAAD
