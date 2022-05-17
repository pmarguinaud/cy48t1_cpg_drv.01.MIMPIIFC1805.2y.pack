SUBROUTINE LATTE_KAPPA(YDGEOMETRY,YDGMV,YDDYN,YDPTRSLB2,KST,KPROF,PDTS2,PSLHDA,PSLHDD0,PGMV,PGMVS,&
 & PUVH0,PRES0,PRDELP,PB2)

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
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE PTRSLB2      , ONLY : TPTRSLB2
USE INTDYN_MOD   , ONLY : YYTHW0
USE YOMDYNA      , ONLY : SLHDKMIN,SLHDKREF
USE YOMDYN       , ONLY : TDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TPTRSLB2)    ,INTENT(IN)    :: YDPTRSLB2
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDA(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDD0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUVH0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTHW0%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRES0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDPTRSLB2%NFLDSLB2)
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
ASSOCIATE(NPROMA=>YDDIM%NPROMA,   NFLEVG=>YDDIMV%NFLEVG, NLEV_SPONGE=>YDDYN%NLEV_SPONGE,   LSLHDHEAT=>YDDYN%LSLHDHEAT,    &
& SLHD_MASK_U=>YDDYN%SLHD_MASK_U, SLHD_MASK_T=>YDDYN%SLHD_MASK_T,   YT0=>YDGMV%YT0,   MSLB2KAPPA=>YDPTRSLB2%MSLB2KAPPA,   &
& MSLB2KAPPAT=>YDPTRSLB2%MSLB2KAPPAT)
!    --------------------------------------------------------


!*       1.   Kappa functions computation

! Momentum
CALL GP_KAPPA(YDVAB,YDDYN,NPROMA,KST,KPROF,NFLEVG,PDTS2,&
 & PGMV(1,1,YT0%MUL),PGMV(1,1,YT0%MVL),PGMV(1,1,YT0%MVOR),PGMV(1,1,YT0%MDIV),&
 & PUVH0(1,0,YYTHW0%M_UH),PUVH0(1,0,YYTHW0%M_VH),PRES0,PGMVS(1,YT0%MSPL),PGMVS(1,YT0%MSPM),&
 & PRDELP,PSLHDA,PSLHDD0,ZKAPPA)

! Heat
IF (LSLHDHEAT) THEN
  CALL GP_KAPPAT(YDGEOMETRY,YDDYN,NPROMA,KST,KPROF,NFLEVG,PDTS2,&
   & PGMV(1,1,YT0%MT),PGMV(1,1,YT0%MTL),PGMV(1,1,YT0%MTM),&
   & PGMVS(1,YT0%MSP),PGMVS(1,YT0%MSPL),PGMVS(1,YT0%MSPM),&
   & PSLHDA,PSLHDD0,ZKAPPAT)
!ELSE
! ZKAPPAT(KST:KPROF,1:NFLEVG)=ZKAPPA(KST:KPROF,1:NFLEVG)
ENDIF

! Specific treatment above tropopause
IF (SLHDKMIN /= SLHDKREF) THEN

  !  Detection of tropopause level
  CALL TROPO_TEP(YDGEOMETRY,NPROMA,KST,KPROF,NFLEVG,&
   & PGMVS(1,YT0%MSP),PGMV(1,1,YT0%MT),ITROP)

  !  Tag areas above troposphere
  DO JROF=KST,KPROF
    DO JLEV=1,MIN(ITROP(JROF),NLEV_SPONGE)
      ZKAPPA(JROF,JLEV)=-ZKAPPA(JROF,JLEV)
    ENDDO
  ENDDO
  IF (LSLHDHEAT) THEN
    DO JROF=KST,KPROF
      DO JLEV=1,MIN(ITROP(JROF),NLEV_SPONGE)
        ZKAPPAT(JROF,JLEV)=-ZKAPPAT(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF

ENDIF

!*       2.   Fill PB2 with Kappa + masking out areas of no use

DO JLEV=1,NLEV_SPONGE
  PB2(KST:KPROF,MSLB2KAPPA+JLEV-1)=SLHD_MASK_U(JLEV)*ZKAPPA(KST:KPROF,JLEV)
  IF (LSLHDHEAT) PB2(KST:KPROF,MSLB2KAPPAT+JLEV-1)= &
   &   SLHD_MASK_T(JLEV)*ZKAPPAT(KST:KPROF,JLEV)
ENDDO
DO JLEV=NLEV_SPONGE+1,NFLEVG
  PB2(KST:KPROF,MSLB2KAPPA+JLEV-1)=0._JPRB
  IF (LSLHDHEAT) PB2(KST:KPROF,MSLB2KAPPAT+JLEV-1)= 0._JPRB
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LATTE_KAPPA',1,ZHOOK_HANDLE)
END SUBROUTINE LATTE_KAPPA
