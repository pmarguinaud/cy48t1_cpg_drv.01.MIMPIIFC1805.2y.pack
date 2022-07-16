SUBROUTINE TRANSDIR_MDL_FROM_T0(YDGEOMETRY,YDGFL,YDGMV,CDCONF,KNFTHER,YDSP)

!**** *TRANSDIR_MDL * - Direct transforms for model

!     Purpose.  Perform direct transform (gridpoint to spectral)
!     --------

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
!        Modified : 03-08-01 M.Hamrud - GFL introduction
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        C.Temperton 04-01-29: Fix for SWE (is there a cleaner way?)
!        R. El Khatib 13-Apr-2018 : fix and replace hard-coded dimensions by soft ones for GMV

!     ------------------------------------------------------------------

USE PARKIND1           , ONLY : JPIM, JPRB
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE YOMGFL             , ONLY : TGFL
USE YOMGMV             , ONLY : TGMV
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE YOMMP0             , ONLY : MYSETV
USE YOMSP              , ONLY : SPVOR_FLT, SPDIV_FLT
USE YOMDYNA            , ONLY : LGRADSP
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)
USE YOMLUN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL)          ,INTENT(INOUT) :: YDGFL
TYPE(TGMV)          ,INTENT(INOUT) :: YDGMV
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KNFTHER 
CHARACTER(LEN=1)    ,INTENT(IN)    :: CDCONF
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSP
!     ------------------------------------------------------------------
REAL(KIND=JPRB)    :: ZSP(1,YDGEOMETRY%YRDIM%NSPEC2)
INTEGER(KIND=JPIM) :: IVSETSC(1)
INTEGER(KIND=JPIM) :: IST,IOFF,IDIM3
LOGICAL :: LLOK, LLG
REAL(KIND=JPRB) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: ISTUV,IENUV,ISTGMV,IENGMV
!     ------------------------------------------------------------------

#include "dir_trans.h"

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TRANSDIR_MDL_FROM_T0',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
 & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YGFL=>YDGFL%YGFL)
ASSOCIATE(NUMFLDS1=>YGFL%NUMFLDS1, NUMSPFLDS1=>YGFL%NUMSPFLDS1, &
 & NPROMA=>YDDIM%NPROMA, NRESOL=>YDDIM%NRESOL, NSPEC2=>YDDIM%NSPEC2, &
 & GFLT1=>YDGFL%GFLT1, GFL=>YDGFL%GFL, &
 & GMV=>YDGMV%GMV, GMVT1=>YDGMV%GMVT1, GMVT1S=>YDGMV%GMVT1S, GMVS=>YDGMV%GMVS, YT0=>YDGMV%YT0, &
 & YT1=>YDGMV%YT1, &
 & NBSETLEV=>YDMP%NBSETLEV, NBSETSP=>YDMP%NBSETSP)
!     ------------------------------------------------------------------

LLOK = CDCONF == 'A'.OR.CDCONF == 'G'
IF(.NOT. LLOK)THEN
  CALL ABOR1('TRANSDIR_MDL - UNKNOWN CONFIGURATION ')
ENDIF
LLG = CDCONF == 'G'

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

!! new defs taken from transinv_mdl
ISTGMV = YT0%NDIMUV+1
IENGMV = YT0%NDIMUV+KNFTHER
ISTUV = YT0%MU
IENUV = YT0%MV

IF(LLG) THEN
! Configuration where grid-point VOR,DIV transformed to spectral (for Jb calc.)
  IF(NUMSPFLDS1>0) THEN
    CALL DIR_TRANS(PSPSC2=ZSP(1:IST,:),PSPSC3A=YDSP%SP3D(:,:,1:2+KNFTHER),&
     & PSPSC3B=YDSP%GFL,&
     & KRESOL=NRESOL,KPROMA=NPROMA,KVSETSC2=IVSETSC(1:1),&
     & KVSETSC3A=NBSETLEV,KVSETSC3B=NBSETLEV,&
     & PGP2=GMVS(:,1:YT1%NDIMS,:),&
     & PGP3A=GMV(:,:,1:YT1%NDIM,:),PGP3B=GFL(:,:,IOFF+1:IOFF+IDIM3,:))  
  ELSE
    CALL DIR_TRANS(PSPSC2=ZSP(1:IST,:),PSPSC3A=YDSP%SP3D(:,:,1:2+KNFTHER),&
     & KRESOL=NRESOL,KPROMA=NPROMA,KVSETSC2=IVSETSC(1:1),&
     & KVSETSC3A=NBSETLEV,&
     & PGP2=GMVS(:,1:YT1%NDIMS,:),&
     & PGP3A=GMV(:,:,1:YT1%NDIM,:))  
  ENDIF
ELSEIF (YT0%NDIM>2) THEN

  IF( LGRADSP ) THEN
    CALL DIR_TRANS(PSPVOR=SPVOR_FLT,PSPDIV=SPDIV_FLT,KVSETUV=NBSETLEV,&
     & KRESOL=NRESOL,KPROMA=NPROMA,PGPUV=GMV(:,:,YT0%MSGRTL:YT0%MSGRTM,:))
  ENDIF
  
! Normal configuration (u,v grid-point to vor,div spectral)
  IF(NUMSPFLDS1>0) THEN
    CALL DIR_TRANS(PSPVOR=YDSP%VOR,PSPDIV=YDSP%DIV,&
     & PSPSC2=ZSP(1:IST,:),PSPSC3A=YDSP%HV,&
     & PSPSC3B=YDSP%GFL,&
     & KRESOL=NRESOL,KPROMA=NPROMA,KVSETUV=NBSETLEV,KVSETSC2=IVSETSC(1:1),&
     & KVSETSC3A=NBSETLEV,KVSETSC3B=NBSETLEV,&
     & PGPUV=GMV(:,:,ISTUV:IENUV,:),PGP2=GMVS(:,1:YT1%NDIMS,:),&
     & PGP3A=GMV(:,:,ISTGMV:IENGMV,:),PGP3B=GFL(:,:,IOFF+1:IOFF+IDIM3,:))  
  ELSE
    WRITE(NULOUT,*) 'think I want to be calling this dir_trans ...'
    CALL DIR_TRANS(PSPVOR=YDSP%VOR,PSPDIV=YDSP%DIV,&
     & PSPSC2=ZSP(1:IST,:),PSPSC3A=YDSP%HV,&
     & KRESOL=NRESOL,KPROMA=NPROMA,KVSETUV=NBSETLEV,KVSETSC2=IVSETSC(1:1),&
     & KVSETSC3A=NBSETLEV,&
     & PGPUV=GMV(:,:,ISTUV:IENUV,:),PGP2=GMVS(:,1:YT1%NDIMS,:),&
     & PGP3A=GMV(:,:,ISTGMV:IENGMV,:))  
  ENDIF
ELSE
! Normal configuration (u,v grid-point to vor,div spectral)
  IF(NUMSPFLDS1>0) THEN
    CALL DIR_TRANS(PSPVOR=YDSP%VOR,PSPDIV=YDSP%DIV,PSPSC2=ZSP(1:IST,:),&
     & PSPSC3B=YDSP%GFL,&
     & KRESOL=NRESOL,KPROMA=NPROMA,KVSETUV=NBSETLEV,KVSETSC2=IVSETSC(1:1),&
     & KVSETSC3B=NBSETLEV,&
     & PGPUV=GMV(:,:,ISTUV:IENUV,:),PGP2=GMVS(:,1:YT1%NDIMS,:),&
     & PGP3A=GMV(:,:,ISTGMV:IENGMV,:),PGP3B=GFL(:,:,IOFF+1:IOFF+IDIM3,:))  
  ELSE
    CALL DIR_TRANS(PSPVOR=YDSP%VOR,PSPDIV=YDSP%DIV,PSPSC2=ZSP(1:IST,:),&
     & KRESOL=NRESOL,KPROMA=NPROMA,KVSETUV=NBSETLEV,KVSETSC2=IVSETSC(1:1),&
     & PGPUV=GMV(:,:,ISTUV:IENUV,:),PGP2=GMVS(:,1:YT1%NDIMS,:))
  ENDIF
ENDIF

CALL GSTATS(1036,0)
IF(NBSETSP == MYSETV) THEN
  YDSP%SP(:) = ZSP(IST,:)
ENDIF
CALL GSTATS(1036,1)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TRANSDIR_MDL_FROM_T0',1,ZHOOK_HANDLE)
END SUBROUTINE TRANSDIR_MDL_FROM_T0
