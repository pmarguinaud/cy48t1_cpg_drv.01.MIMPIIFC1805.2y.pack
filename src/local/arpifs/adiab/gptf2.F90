SUBROUTINE GPTF2(YDGEOMETRY,YDGMV,&
 ! --- INPUT ---------------------------------------------------------
 & YDML_GCONF,YDDYN,KST,KEN,LDFSTEP,&
 ! --- INPUT-OUTPUT --------------------------------------------------
 & PGMV,PGMVS,PGFL)

!**** *GPTF2* - Timefilter part 2

!     Purpose.
!     --------
!           Performs part 2 of the time-filtering.
!           - leap frog + ldfstep=true : pxt9. = pxt0.
!           - leap frog + ldfstep=false: pxt9. = pxt9. + eps2*pxt0.
!           - sl2tl     + ldfstep=true : put9=put0 and pvt9=pvt0 only.
!           - sl2tl     + ldfstep=false: nothing.

!**   Interface.
!     ----------
!        *CALL* *GPTF2(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!         KST          : start of work
!         KEN          : depth of work
!         LDLSTEP      : check on the last time step.

!        INPUT/OUTPUT:
!         PGMV         : "t" and "t-dt" upper air GMV variables.
!         PGMVS        : "t" and "t-dt" surface GMV variables.
!         PGFL         : "t" and "t-dt" GFL variables.

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud  *ECMWF*
!      Original : 92-02-01

! Modifications.
! --------------
!   Modified 02-09-30 by P. Smolikova (variable d4 in NH)
!   Modified 13-11-02 K. YESSAD : cleanings + improve vectorization.
!   Modified 2003-07-17 C. Fischer - psvdauxt0* come from pgmv
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!   08-Jun-2004 J. Masek   NH cleaning (LVSLWBC)
!   Modified Nov 2007 N. Wedi: bug correction
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2008): remove dummy CDLOCK
! End Modifications
!------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCT0       , ONLY : LTWOTL, LNHDYN
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : NVDVAR

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEN 
LOGICAL           ,INTENT(IN)    :: LDFSTEP 
REAL(KIND=JPRB),TARGET,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV) 
REAL(KIND=JPRB),TARGET,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM) 

#include "gptf2_expl_2tl.intfb.h"
#include "gptf2_expl_3tl.intfb.h"

!     ------------------------------------------------------------------
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:)   :: Z0SP  , Z0SPL , Z0SPM , Z9SP  , Z9SPL , Z9SPM 
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z0DIV , Z0NHX , Z0SPD , Z0SPDL, Z0SPDM, Z0SVD , Z0SVDL
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z0SVDM, Z0T   , Z0TL  , Z0TM  , Z0U   , Z0V   , Z9DIV 
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z9NHX , Z9SPD , Z9SPDL, Z9SPDM, Z9SVD , Z9SVDL, Z9SVDM
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z9T   , Z9TL  , Z9TM  , Z9U   , Z9V   

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPTF2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YGFL=>YDML_GCONF%YGFL,YDDIMF=>YDML_GCONF%YRDIMF)
ASSOCIATE(NDIM=>YGFL%NDIM, NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFTHER=>YDDIMF%NFTHER, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & REPS2=>YDDYN%REPS2, REPSM2=>YDDYN%REPSM2, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0, &
 & YT9=>YDGMV%YT9)
!     ------------------------------------------------------------------

Z0SP    => NULL ()
Z0SPL   => NULL ()
Z0SPM   => NULL ()
Z9SP    => NULL ()
Z9SPL   => NULL ()
Z9SPM   => NULL ()

Z0DIV   => NULL ()
Z0NHX   => NULL ()
Z0SPD   => NULL ()
Z0SPDL  => NULL ()
Z0SPDM  => NULL ()
Z0SVD   => NULL ()
Z0SVDL  => NULL ()
Z0SVDM  => NULL ()
Z0T     => NULL ()
Z0TL    => NULL ()
Z0TM    => NULL ()
Z0U     => NULL ()
Z0V     => NULL ()
Z9DIV   => NULL ()
Z9NHX   => NULL ()
Z9SPD   => NULL ()
Z9SPDL  => NULL ()
Z9SPDM  => NULL ()
Z9SVD   => NULL ()
Z9SVDL  => NULL ()
Z9SVDM  => NULL ()
Z9T     => NULL ()
Z9TL    => NULL ()
Z9TM    => NULL ()
Z9U     => NULL ()
Z9V     => NULL ()
           
IF (YT0%MSP  > 0) Z0SP  => PGMVS(:,YT0%MSP)                       
IF (YT0%MSPL > 0) Z0SPL => PGMVS(:,YT0%MSPL)                       
IF (YT0%MSPM > 0) Z0SPM => PGMVS(:,YT0%MSPM)                       
IF (YT9%MSP  > 0) Z9SP  => PGMVS(:,YT9%MSP)                       
IF (YT9%MSPL > 0) Z9SPL => PGMVS(:,YT9%MSPL)                       
IF (YT9%MSPM > 0) Z9SPM => PGMVS(:,YT9%MSPM)                       

IF (YT0%MDIV  > 0)  Z0DIV  => PGMV(:,:,YT0%MDIV)                       
IF (YT0%MNHX  > 0)  Z0NHX  => PGMV(:,:,YT0%MNHX)                       
IF (YT0%MSPD  > 0)  Z0SPD  => PGMV(:,:,YT0%MSPD)                       
IF (YT0%MSPDL > 0)  Z0SPDL => PGMV(:,:,YT0%MSPDL)                       
IF (YT0%MSPDM > 0)  Z0SPDM => PGMV(:,:,YT0%MSPDM)                       
IF (YT0%MSVD  > 0)  Z0SVD  => PGMV(:,:,YT0%MSVD)                       
IF (YT0%MSVDL > 0)  Z0SVDL => PGMV(:,:,YT0%MSVDL)                       
IF (YT0%MSVDM > 0)  Z0SVDM => PGMV(:,:,YT0%MSVDM)                       
IF (YT0%MT    > 0)  Z0T    => PGMV(:,:,YT0%MT)                       
IF (YT0%MTL   > 0)  Z0TL   => PGMV(:,:,YT0%MTL)                       
IF (YT0%MTM   > 0)  Z0TM   => PGMV(:,:,YT0%MTM)                       
IF (YT0%MU    > 0)  Z0U    => PGMV(:,:,YT0%MU)                       
IF (YT0%MV    > 0)  Z0V    => PGMV(:,:,YT0%MV)                       
IF (YT9%MDIV  > 0)  Z9DIV  => PGMV(:,:,YT9%MDIV)                       
IF (YT9%MNHX  > 0)  Z9NHX  => PGMV(:,:,YT9%MNHX)                       
IF (YT9%MSPD  > 0)  Z9SPD  => PGMV(:,:,YT9%MSPD)                       
IF (YT9%MSPDL > 0)  Z9SPDL => PGMV(:,:,YT9%MSPDL)                       
IF (YT9%MSPDM > 0)  Z9SPDM => PGMV(:,:,YT9%MSPDM)                       
IF (YT9%MSVD  > 0)  Z9SVD  => PGMV(:,:,YT9%MSVD)                       
IF (YT9%MSVDL > 0)  Z9SVDL => PGMV(:,:,YT9%MSVDL)                       
IF (YT9%MSVDM > 0)  Z9SVDM => PGMV(:,:,YT9%MSVDM)                       
IF (YT9%MT    > 0)  Z9T    => PGMV(:,:,YT9%MT)                       
IF (YT9%MTL   > 0)  Z9TL   => PGMV(:,:,YT9%MTL)                       
IF (YT9%MTM   > 0)  Z9TM   => PGMV(:,:,YT9%MTM)                       
IF (YT9%MU    > 0)  Z9U    => PGMV(:,:,YT9%MU)                       
IF (YT9%MV    > 0)  Z9V    => PGMV(:,:,YT9%MV)                       

IF (LTWOTL) THEN
  CALL GPTF2_EXPL_2TL (YDGEOMETRY, KST, KEN, LDFSTEP, Z0U, Z0V, Z9U, Z9V)
ELSE
  CALL GPTF2_EXPL_3TL (YDGEOMETRY, YDML_GCONF, YDDYN, KST, KEN, LDFSTEP, LNHDYN, NVDVAR, &
 & PGFL, &
 & Z0SP  , Z0SPL , Z0SPM , Z9SP  , Z9SPL , Z9SPM, &
 & Z0DIV , Z0NHX , Z0SPD , Z0SPDL, Z0SPDM, Z0SVD , Z0SVDL, &
 & Z0SVDM, Z0T   , Z0TL  , Z0TM  , Z0U   , Z0V   , Z9DIV,  &
 & Z9NHX , Z9SPD , Z9SPDL, Z9SPDM, Z9SVD , Z9SVDL, Z9SVDM, &
 & Z9T   , Z9TL  , Z9TM  , Z9U   , Z9V)
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPTF2',1,ZHOOK_HANDLE)
END SUBROUTINE GPTF2
