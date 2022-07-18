!OCL NOEVAL
SUBROUTINE TRANSDIR_NHCONV(YDGEOMETRY,LDMODEL_TO_FILE,PSPD,PSVD,PSNHX,YDSP)

!* TRANSDIR_NHCONV - Conversion of NH variables (model to file and vice-versa)
!                    Direct spectral transforms.

! Purpose
! -------

! Interface
! ---------
!  * INPUT: 
!    LDMODEL_TO_FILE : switch for model to file / file to model conversion
!    PSPD            : NH pressure departure variable
!    PSVD            : NH vertical divergence variable
!    PSNHX           : NH "X" part of the total divergence.     
!
!  * OUTPUT: 
!    YDSP            : spectral variable in which transformed information is written

! Externals
! ---------

! Method
! ------
!   See documentation

! Reference
! ---------

! Author
! ------
!   27-Jan-2005 K. Yessad (after old part 5 of GNHPDVDCONV)

! Modifications
! -------------
!   18-Jul-2005 R. Brozkova: old VDAUX => "X" term         
!   13-Oct-2005 K. Yessad  : B-level parallelisation.
!      Aug 2016 O. Marsden : Removed use of SPA3, replaced by explicit spectral argument
! End Modifications
!------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMDYNA  , ONLY : YRDYNA
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
LOGICAL             ,INTENT(IN)    :: LDMODEL_TO_FILE 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PSPD(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)     ,INTENT(IN)    :: PSVD(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)     ,INTENT(IN)    :: PSNHX(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSP

!------------------------------------------------------------------------------

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------

#include "reespe.intfb.h"

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TRANSDIR_NHCONV',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, NFLSUR=>YDDIMV%NFLSUR, &
 & NGPTOT=>YDGEM%NGPTOT)

!------------------------------------------------------------------------------

!* 1. TRANSFORM NH FIELDS BACK TO SPECTRAL SPACE
!     ------------------------------------------

! * pressure departure variable
CALL REESPE(YDGEOMETRY,NFLSUR,NFLEVG,YDSP%SPD,PSPD)

! * pseudo vertical divergence variable
CALL REESPE(YDGEOMETRY,NFLSUR,NFLEVG,YDSP%SVD,PSVD)

! * transform only if transformation
!   from file to model (to initialize NHX variable)
IF ( YRDYNA%LNHX.AND.(.NOT.LDMODEL_TO_FILE) ) THEN
  ! "X" part of the total divergence
  CALL REESPE(YDGEOMETRY,NFLSUR,NFLEVG,YDSP%NHX,PSNHX)
ENDIF

!------------------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TRANSDIR_NHCONV',1,ZHOOK_HANDLE)

END SUBROUTINE TRANSDIR_NHCONV

