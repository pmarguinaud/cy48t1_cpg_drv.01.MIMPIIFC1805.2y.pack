SUBROUTINE SIMREL(YDGEOMETRY,KINDIC,KN,PX,PF,PG,KITER)

!**** *SIMREL* ****

!     PURPOSE.
!     --------

!       External of the minimisation routine *M1QN3R* : calculates the function
!     and its gradient.

!**   INTERFACE.
!     ----------

!     *CALL* *SIMREL(...)*

!         *KINDIC* - Required by the minimization package
!         *KN*     - Size of the argument (input)
!         *PX*     - Argument (input)
!         *PF*     - Value of the function
!         *PG*     - Value of the gradient
!         *KITER*  - Required by the minimization package

!     METHOD.
!     -------

!       Transforms into grid point , calculates over the grid the quadratic
!     function and its derivatives, then returns to spectral coefficients.

!     EXTERNALS.
!     ----------

!           SPREORD
!           REESPE
!           SPEREE
!           ESPERAD
!           ESPEREE

!     AUTHOR.
!     -------
!      Y. BOUTELOUP and M. DEQUE   1 FEB 91.

!     MODIFICATIONS.
!     --------------
!      Modified : 01-03-08 S. Alexandru, management of the extension zone
!      Modified : 01-11-23 D. Giard, clean call to ESPERAD
!      C. Fischer:03-04-07 implement control_vector type
!      Modified : 02-12-05 D. Giard, cleaner call to ESPERAD
!      Modified : 02-12-16 D. Giard, LNEWORO3-> LNEWENV
!      Modified : 03-01-13 D. Giard, spectral smoothing and cleaning
!                          (from K. Stadlbacher)
!      Modified : 05-03-25 F. Taillefer, back to an old minimization code
!      Modified : 05-04-04 D. Giard, update, semi-envelope and old ARPEGE 
!                          cost functions, adding a polynomial one
!      Modified : Feb 2011 G.Mozdzynski OOPS cleaning, use of derived type TCSGLEG
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE PTRSPOR  , ONLY : XREV     ,XPDS     ,XPDL     ,XEXT
USE YOMCLA   , ONLY : HDIM     ,QPOWER   ,QCONST   ,FACE     ,&
 & LNEWORO  , LNEWORO2 ,NLISSP
USE YOMCT0   , ONLY : LELAM

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(INOUT) :: KINDIC
INTEGER(KIND=JPIM),INTENT(IN)    :: KN
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PX(KN) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PF 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PG(KN) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KITER

REAL(KIND=JPRB) :: ZS(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZX1(YDGEOMETRY%YRDIM%NSPEC2),ZG1(YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) :: ZPX(YDGEOMETRY%YRDIM%NSEFRE),ZPG(YDGEOMETRY%YRDIM%NSEFRE)
REAL(KIND=JPRB),ALLOCATABLE :: ZST(:,:,:),ZGT(:,:)
REAL(KIND=JPRB) :: ZCST, ZDF, ZDS, ZW
REAL(KIND=JPRB) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: IADL(YDGEOMETRY%YRDIM%NDGLG)
INTEGER(KIND=JPIM) :: IADR, JI, JK, JM, JN

#include "abor1.intfb.h"
#include "esperad.intfb.h"
#include "esperee.intfb.h"
#include "reespe.intfb.h"
#include "speree.intfb.h"
#include "spreord.intfb.h"

IF (LHOOK) CALL DR_HOOK('SIMREL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDGEM=>YDGEOMETRY%YRGEM, YDLAP=>YDGEOMETRY%YRLAP, YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
& YDELAP=>YDGEOMETRY%YRELAP)
ASSOCIATE(NDGLG=>YDDIM%NDGLG, NMSMAX=>YDDIM%NMSMAX, NSEFRE=>YDDIM%NSEFRE,   NSMAX=>YDDIM%NSMAX, NSPEC2=>YDDIM%NSPEC2,   &
& NGPTOT=>YDGEM%NGPTOT, NLOENG=>YDGEM%NLOENG, NSTAGP=>YDGEM%NSTAGP,   NTSTAGP=>YDGEM%NTSTAGP,   NASN0=>YDLAP%NASN0      &
& )

!     ------------------------------------------------------------------

!*    1. SPECTRAL - GRID POINT TRANSFORM OF THE ARGUMENT
!        -----------------------------------------------

IF (KN /= NSEFRE) CALL ABOR1('SIMREL : DIMENSION MISMATCH !')

IF (LELAM) THEN
  IADL(1:NDGLG)=NTSTAGP(1:NDGLG)-1
ELSE
  IADL(1:NDGLG)=NSTAGP(1:NDGLG)-1
ENDIF  

ZPX(:)=PX(:)
CALL SPREORD(YDGEOMETRY%YRDIM,YDGEOMETRY%YREDIM,YDGEOMETRY%YRELAP,1,ZPX,ZX1,.TRUE.)
IF (LELAM) THEN
  CALL ESPEREE(YDGEOMETRY,1,1,ZX1,ZS)
ELSE
  CALL SPEREE(YDGEOMETRY,1,1,ZX1,ZS)
ENDIF

!     ------------------------------------------------------------------

!*    2. CALCULATION OF THE FUNCTION AND ITS GRADIENT.
!        ---------------------------------------------

PF=0.0_JPRB

IF (LNEWORO) THEN
  DO JK=1,NDGLG
    ZW=YDCSGLEG%RW(JK)/REAL(NLOENG(JK),JPRB)
    DO JI=1,NLOENG(JK)
      IADR=JI+IADL(JK)
      ZDF=(ABS(ZS(IADR)-XREV(IADR))/HDIM)**XPDS(IADR)
      ZDS=XPDS(IADR)*(ZS(IADR)-XREV(IADR))*&
       & (ABS(ZS(IADR)-XREV(IADR))/HDIM)**(XPDS(IADR)-2.0_JPRB)/HDIM**2  
      PF=PF+XEXT(JI,JK)*ZW*ZDF
      ZS(IADR)=XEXT(JI,JK)*ZDS
    ENDDO
  ENDDO

ELSEIF (LNEWORO2) THEN
  ZCST=QCONST**QPOWER
  DO JK=1,NDGLG
    ZW=YDCSGLEG%RW(JK)/REAL(NLOENG(JK),JPRB)
    DO JI=1,NLOENG(JK)
      IADR=JI+IADL(JK)
      ZDF=FACE*XPDS(IADR)*(ZS(IADR)-XREV(IADR))**2&
       & +ZCST*ABS(ZS(IADR)-XREV(IADR))**QPOWER  
      ZDS=2.0_JPRB*FACE*XPDS(IADR)*(ZS(IADR)-XREV(IADR))&
       & +ZCST*QPOWER*(ZS(IADR)-XREV(IADR))&
       & *ABS(ZS(IADR)-XREV(IADR))**(QPOWER-2.0_JPRB)  
      PF=PF+XEXT(JI,JK)*ZW*ZDF
      ZS(IADR)=XEXT(JI,JK)*ZDS
    ENDDO
  ENDDO
ENDIF

IF (LELAM) THEN
  DO JK=1,NDGLG
    ZW=YDCSGLEG%RW(JK)/REAL(NLOENG(JK),JPRB)
    DO JI=1,NLOENG(JK)
      IADR=JI+IADL(JK)
      ZS(IADR)=ZS(IADR)*ZW
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*    3. GRID POINT - SPECTRAL TRANSFORM OF THE GRADIENT
!        -----------------------------------------------

IF (LELAM) THEN
  ALLOCATE ( ZST (NGPTOT,1,1) )
  ALLOCATE ( ZGT (1,NSPEC2) )
  ZST(:,:,:)=0.0_JPRB
  ZGT(:,:)  =0.0_JPRB
  ZST(:,1,1)=ZS(:)
  CALL ESPERAD(YDGEOMETRY,1,1,ZGT,ZST)
  ZG1(:)=ZGT(1,:)
  DEALLOCATE ( ZST )
  DEALLOCATE ( ZGT )
ELSE
  CALL REESPE(YDGEOMETRY,1,1,ZG1,ZS)
ENDIF
CALL SPREORD(YDGEOMETRY%YRDIM,YDGEOMETRY%YREDIM,YDGEOMETRY%YRELAP,1,ZPG,ZG1,.FALSE.)

!     ------------------------------------------------------------------

!*    4. ADDITION OF A SPECTRAL TERM TO THE FUNCTION AND THE GRADIENT
!        ------------------------------------------------------------

IF (NLISSP == 2) THEN
 IF (LELAM) THEN
  DO JM=0,NMSMAX
    DO JN=0,YDELAP%NCPL4M(JM)-1
      IADR=YDELAP%NESM0(JM)+JN
      PF=PF+2*XPDL(IADR)*ZPX(IADR)**2
      ZPG(IADR)=ZPG(IADR)+4*XPDL(IADR)*ZPX(IADR)
    ENDDO
  ENDDO
 ELSE
  DO JN=0,NSMAX
    DO JM=-JN,JN
      IADR=NASN0(JN)+JM
      IF (JM == 0) THEN
        PF=PF+XPDL(IADR)*ZPX(IADR)**2
        ZPG(IADR)=ZPG(IADR)+2*XPDL(IADR)*ZPX(IADR)
      ELSE
        PF=PF+2*XPDL(IADR)*ZPX(IADR)**2
        ZPG(IADR)=ZPG(IADR)+4*XPDL(IADR)*ZPX(IADR)
      ENDIF
    ENDDO
  ENDDO
 ENDIF 
ENDIF

!     ------------------------------------------------------------------

PG(:)=ZPG(:)
PX(:)=ZPX(:)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SIMREL',1,ZHOOK_HANDLE)
END SUBROUTINE SIMREL
