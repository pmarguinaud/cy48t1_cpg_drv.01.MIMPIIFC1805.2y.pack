SUBROUTINE CPG_END_TL(YDGEOMETRY,YDGMV, &
 !---------------------------------------------------------------------
 ! - INPUT .
 & YDML_GCONF,YDDYN,KST,KEND,PGM,&
 !---------------------------------------------------------------------
 ! - INPUT/OUTPUT .
 & PGFL,PGMV,PGMVS)

!**** *CPG_END_TL* - Grid point calculations: end of non lagged part, TL code.

!     Purpose.
!     --------
!           Grid point calculations: end of non lagged part, TL code.
!           Note that this routine looks like a subset of CPG_END:
!           - part 3 is quasi-identical to part 3 of CPG_END.
!           - part 5 is a subset of the part 5 of CPG_END.

!**   Interface.
!     ----------
!        *CALL* *CPG_END_TL(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        KST       : first element of work.
!        KEND      : last element of work.
!        PGM       : mapping factor.

!     INPUT/OUTPUT:
!     -------------
!        PGFL      : GFL variables at time t-dt and t.
!        PGMV      : upper air GMV variables at time t and t-dt.
!        PGMVS     : surface GMV variables at time t and t-dt.

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

!     Author.
!     -------
!        K. Yessad, after CPGTL (old part 4 of CPGTL)

! Modifications
! -------------
!   Original     : 12 Jul 2004.
!   17-Apr-2007 S.Ivatek-S: Over dimensioning of PGPNH to NFGPNH+1, boundary
!                           checking problem if NFGPNH=0 bf
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (Aug 2009): remove LPC_OLD in TL and AD codes.
!   K. Yessad (Nov 2011): new GPMPFC.
! End Modifications
!------------------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMDYN       , ONLY : TDYN
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV   , ONLY : TGMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TDYN)        ,INTENT(INOUT) :: YDDYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS)
!     ------------------------------------------------------------------
LOGICAL :: LLSTR

INTEGER(KIND=JPIM) :: IFLAG

REAL(KIND=JPRB) :: ZEPS
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gpmpfc_pgfl.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPG_END_TL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YGFL=>YDML_GCONF%YGFL)
ASSOCIATE(NDIM=>YGFL%NDIM, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & RSTRET=>YDGEM%RSTRET, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0)
!     ------------------------------------------------------------------

!*    1.    PRELIMINARY INITIALISATIONS.
!           ----------------------------

ZEPS=100.0_JPRB*TINY(1.0_JPRB)
LLSTR=(ABS(RSTRET-1.0_JPRB)>ZEPS)

!     ------------------------------------------------------------------

!*    3.    DIVIDE BY MAPPING FACTOR.
!           -------------------------

IF(LLSTR) THEN
  IFLAG=1
  CALL GPMPFC_PGFL(YDGMV,YDML_GCONF,YDDYN,NPROMA,NFLEVG,KST,KEND,IFLAG,PGM,PGMV,PGMVS,PGFL)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG_END_TL',1,ZHOOK_HANDLE)
END SUBROUTINE CPG_END_TL
