SUBROUTINE ETRANSDIR_MDL_FROM_T0(YDGEOMETRY,YDGFL,YDGMV,CDCONF,KNFTHER,YDSP)

!**** *ETRANSDIR_MDL_FROM_T0 * - Direct transforms for model

!     Purpose.  Perform direct transform (gridpoint to spectral)
!     --------

!**   Interface.  CALL ETRANSDIR_MDL_FROM_T0(...)
!     ----------

!     Explicit arguments :
!     --------------------
!        CDCONF - configuration of work

!        Called by TRANSDIRH

!     Externals.
!     ----------
!     DIR_TRANS - inverse transform (TRANS library)

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-10-25
!     Modified : 02-09-30 P. Smolikova (EBIPAUX called for biperiodic. if d4 in NH)
!     R. El Khatib 03-05-07 Conf 'U' where geographical wind is transformed the reduced
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

!        Y. Seity and G. Radnoti 03-09-29 phasing for AL27 (new data flow)
!        G. Hello 03-11 bug correction for auxillary variable (LVDAUX)
!        O. Spaniel Oct-2004 phasing for AL29
!        R. Brozkova Jul-2005 VDAUX structure => standard GMV structures:
!                             cleaning the specific VDAUX structure
!        K. Yessad (Aug 2009): cleaning
!        B. Bochenek 15-04-2013 - Phasing cy40, coherence with modified modules
!        B. Bochenek (Apr 2015): Phasing: update
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGFL   , ONLY : TGFL
USE YOMGMV   , ONLY : TGMV
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE YOMSP    , ONLY : SPVOR_FLT, SPDIV_FLT, SPUB_FLT, SPVB_FLT
USE YOMMP0    ,ONLY : MYSETV
USE YOMDYNA   ,ONLY : YRDYNA
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TGFL) , INTENT(INOUT) :: YDGFL
TYPE(TGMV) , INTENT(INOUT) :: YDGMV
CHARACTER(LEN=1),INTENT(IN) :: CDCONF
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KNFTHER
TYPE(SPECTRAL_FIELD), INTENT(INOUT) :: YDSP
!     ------------------------------------------------------------------
REAL(KIND=JPRB),       ALLOCATABLE :: ZSPU(:,:), ZSPV(:,:)
REAL(KIND=JPRB)    :: ZSP(1,YDGEOMETRY%YRDIM%NSPEC2)
INTEGER(KIND=JPIM) :: IVSETSC(1)
INTEGER(KIND=JPIM) :: IST,IOFF,IDIM3,IDIMT1
LOGICAL :: LLOK, LLG, LLU
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: ISTUV,IENUV,ISTGMV,IENGMV
!     ------------------------------------------------------------------

#include "edir_trans.h"

#include "abor1.intfb.h"
#include "euvgeovd.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ETRANSDIR_MDL_FROM_T0',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
  & YGFL=>YDGFL%YGFL)
ASSOCIATE(NUMFLDS1=>YGFL%NUMFLDS1, NUMSPFLDS1=>YGFL%NUMSPFLDS1, &
 & NSPEC2=>YDDIM%NSPEC2, NPROMA=>YDDIM%NPROMA, NRESOL=>YDDIM%NRESOL, &
 & NBSETLEV=>YDMP%NBSETLEV, NBSETSP=>YDMP%NBSETSP, &
 & GMVT1=>YDGMV%GMVT1, GMVT1S=>YDGMV%GMVT1S, YT1=>YDGMV%YT1, &
 & GMV=>YDGMV%GMV, GMVS=>YDGMV%GMVS, YT0=>YDGMV%YT0, &
 & NFLEVL=>YDDIMV%NFLEVL, NFLSUR=>YDDIMV%NFLSUR)
!     ------------------------------------------------------------------

LLOK = CDCONF == 'A'.OR.CDCONF == 'G'.OR.CDCONF == 'U'
IF(.NOT. LLOK)THEN
  CALL ABOR1('TRANSDIR_MDL - UNKNOWN CONFIGURATION ')
ENDIF
LLG = CDCONF == 'G'
LLU = CDCONF == 'U'

IVSETSC(1) = NBSETSP
IF(NBSETSP == MYSETV) THEN
  IST = 1
ELSE
  IST = 0
ENDIF

IF(NUMSPFLDS1 > 0) THEN
  IOFF  = NUMFLDS1-NUMSPFLDS1
  IDIM3 = NUMSPFLDS1
ENDIF

IDIMT1=YT1%NDIM
ISTGMV = YT0%NDIMUV+1
IENGMV = YT0%NDIMUV+KNFTHER
ISTUV = YT0%MU
IENUV = YT0%MV

IF(LLG.OR.LLU) THEN
! Configuration where grid-point VOR,DIV transformed to spectral (for Jb calc.)
  IF(NUMSPFLDS1>0) THEN
    CALL EDIR_TRANS(PSPSC2=ZSP(1:IST,:),PSPSC3A=YDSP%SP3D(:,:,1:2+KNFTHER),&
     & PSPSC3B=YDSP%GFL,&
     & KRESOL=NRESOL,KPROMA=NPROMA,KVSETSC2=IVSETSC(1:1),&
     & KVSETSC3A=NBSETLEV,KVSETSC3B=NBSETLEV,&
     & PGP2=GMVS(:,1:YT1%NDIMS,:),&
     & PGP3A=GMV(:,:,1:IDIMT1,:),PGP3B=YDGFL%GFL(:,:,IOFF+1:IOFF+IDIM3,:))
  ELSE
    CALL EDIR_TRANS(PSPSC2=ZSP(1:IST,:),PSPSC3A=YDSP%SP3D(:,:,1:2+KNFTHER),&
     & KRESOL=NRESOL,KPROMA=NPROMA,KVSETSC2=IVSETSC(1:1),&
     & KVSETSC3A=NBSETLEV,&
     & PGP2=GMVS(:,1:YT1%NDIMS,:),&
     & PGP3A=GMV(:,:,1:IDIMT1,:))
  ENDIF
ELSEIF(YT0%NDIM>2) THEN
  IF( YRDYNA%LGRADSP ) THEN
    CALL EDIR_TRANS(PSPVOR=SPVOR_FLT,PSPDIV=SPDIV_FLT,KVSETUV=NBSETLEV,&
     & KRESOL=NRESOL,KPROMA=NPROMA,PGPUV=GMV(:,:,YT0%MSGRTL:YT0%MSGRTM,:),&
     & PMEANU=SPUB_FLT,PMEANV=SPVB_FLT)
  ENDIF
! Normal configuration (u,v grid-point to vor,div spectral)
  IF(NUMSPFLDS1>0) THEN
    CALL EDIR_TRANS(PSPVOR=YDSP%VOR,PSPDIV=YDSP%DIV,PSPSC2=ZSP(1:IST,:),PSPSC3A=YDSP%HV,&
     & PSPSC3B=YDSP%GFL,&
     & KRESOL=NRESOL,KPROMA=NPROMA,KVSETUV=NBSETLEV,KVSETSC2=IVSETSC(1:1),&
     & KVSETSC3A=NBSETLEV,KVSETSC3B=NBSETLEV,&
     & PGPUV=GMV(:,:,ISTUV:IENUV,:),PGP2=GMVS(:,1:YT1%NDIMS,:),&
     & PGP3A=GMV(:,:,ISTGMV:IENGMV,:),PGP3B=YDGFL%GFL(:,:,IOFF+1:IOFF+IDIM3,:),&
     & PMEANU=YDSP%UB,PMEANV=YDSP%VB)
  ELSE
    CALL EDIR_TRANS(PSPVOR=YDSP%VOR,PSPDIV=YDSP%DIV,PSPSC2=ZSP(1:IST,:),PSPSC3A=YDSP%HV,&
     & KRESOL=NRESOL,KPROMA=NPROMA,KVSETUV=NBSETLEV,KVSETSC2=IVSETSC(1:1),&
     & KVSETSC3A=NBSETLEV,&
     & PGPUV=GMV(:,:,ISTUV:IENUV,:),PGP2=GMVS(:,1:YT1%NDIMS,:),&
     & PGP3A=GMV(:,:,ISTGMV:IENGMV,:),PMEANU=YDSP%UB,PMEANV=YDSP%VB)
  ENDIF
ELSE
! Normal configuration (u,v grid-point to vor,div spectral)
  IF(NUMSPFLDS1>0) THEN
    CALL EDIR_TRANS(PSPVOR=YDSP%VOR,PSPDIV=YDSP%DIV,PSPSC2=ZSP(1:IST,:),&
     & PSPSC3B=YDSP%GFL,&
     & KRESOL=NRESOL,KPROMA=NPROMA,KVSETUV=NBSETLEV,KVSETSC2=IVSETSC(1:1),&
     & KVSETSC3B=NBSETLEV,&
     & PGPUV=GMV(:,:,ISTUV:IENUV,:),PGP2=GMVS(:,1:YT1%NDIMS,:),&
     & PGP3A=GMV(:,:,ISTGMV:IENGMV,:),PGP3B=YDGFL%GFL(:,:,IOFF+1:IOFF+IDIM3,:),&
     & PMEANU=YDSP%UB,PMEANV=YDSP%VB)
  ELSE
    CALL EDIR_TRANS(PSPVOR=YDSP%VOR,PSPDIV=YDSP%DIV,PSPSC2=ZSP(1:IST,:),&
     & KRESOL=NRESOL,KPROMA=NPROMA,KVSETUV=NBSETLEV,KVSETSC2=IVSETSC(1:1),&
     & PGPUV=GMV(:,:,ISTUV:IENUV,:),PGP2=GMVS(:,1:YT1%NDIMS,:),&
     & PMEANU=YDSP%UB,PMEANV=YDSP%VB)
  ENDIF
ENDIF

IST=1

IF(NFLEVL > 0 .AND. LLU) THEN
  ALLOCATE(ZSPU(NFLSUR,NSPEC2))
  ALLOCATE(ZSPV(NFLSUR,NSPEC2))
  ! For this specific case we use YT1%MU. An other solution
  ! would have been to use SPVOR and SPDIV
  ZSPU(1:NFLEVL,:) = YDSP%SP3D(1:NFLEVL,:,YT1%MU)
  ZSPV(1:NFLEVL,:) = YDSP%SP3D(1:NFLEVL,:,YT1%MV)
  CALL EUVGEOVD(YDGEOMETRY,ZSPU,ZSPV,YDSP%VOR,YDSP%DIV,YDSP%UB,YDSP%VB)
  DEALLOCATE(ZSPU)
  DEALLOCATE(ZSPV)
ENDIF

CALL GSTATS(1859,0)
IF(NBSETSP == MYSETV) THEN
  YDSP%SP(:) = ZSP(IST,:)
ENDIF
CALL GSTATS(1859,1)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ETRANSDIR_MDL_FROM_T0',1,ZHOOK_HANDLE)
END SUBROUTINE ETRANSDIR_MDL_FROM_T0
