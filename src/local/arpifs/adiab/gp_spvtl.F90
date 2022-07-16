SUBROUTINE GP_SPVTL(&
 ! - INPUT --------------------------------------------------------------
 & YDGEOMETRY, YDDYN,YDSIMPHL,KST,KEND,&
 & PSPT0,PSPT0L,PSPT0M,PSPT9,&
 ! - OUTPUT -------------------------------------------------------------
 & PRE0,PRE0L,PRE0M,PRE9,&
 ! - INPUT-TRAJECTORY ---------------------------------------------------
 & PSPT5L,PSPT5M,PRE5,PRE95)

!------------------------------------------------------------------------------
!**** *GP_SPVTL* - 

! Purpose
! -------

! Interface
! ---------
!   abbreviation "prehyds" means "hydrostatic surface pressure".

!   INPUT:
!   ------
!   KST                : start of work
!   KEND               : depth of work
!   PSPT0,PSPT0L,PSPT0M: log(prehyds) and derivatives at t.
!   PSPT9              : log(prehyds) at t-dt.

!   OUTPUT:
!   -------
!   PRE0,PRE0L,PRE0M   : prehyds and derivatives at t.
!   PRE9               : prehyds at t-dt.

!   INPUT-TRAJECTORY:
!   -----------------
!   PSPT5L,PSPT5M      : cf. PSPT0L,PSPT0M.
!   PRE5               : cf. PRE0.
!   PRE95              : cf. PRE9.

! Externals
! ---------

! Method
! ------

! Reference
! ---------
       
! Author
! ------
!   27-Jun-2002 C. Fischer (METEO-FRANCE)

! Modifications
! -------------
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!      Dec 2003 K. Yessad  multiplication by GM has moved in GPMPFC_GMVS.
!   09-Jun-2004 J. Masek   NH cleaning (LFULLIMP)
!   01-Jul-2004 K. Yessad  Make clearer the tests for PC scheme.
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (Nov 2012): simplify testings.
!   K. Yessad (July 2014): Move some variables.
! End Modifications
!------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCT0       , ONLY : LTWOTL, LNHDYN
USE YOMDYNA      , ONLY : LPC_FULL
USE YOMDYN       , ONLY : TDYN
USE YOMSIMPHL    , ONLY : TSIMPHL

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TSIMPHL)     ,INTENT(IN)    :: YDSIMPHL
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0L(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0M(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT9(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0L(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0M(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT5L(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT5M(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE5(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE95(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF
LOGICAL :: LLPRE9
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GP_SPVTL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
 & LSIMPH=>YDSIMPHL%LSIMPH)
!     ------------------------------------------------------------------

IF (.NOT.(LSIMPH.OR.LNHDYN)) THEN
  LLPRE9=.FALSE.
ELSE
  IF (LTWOTL) THEN
    LLPRE9=(NCURRENT_ITER > 0 .AND. LPC_FULL)
  ELSE
    LLPRE9=.TRUE.
  ENDIF
ENDIF

DO JROF=KST,KEND
  PRE0 (JROF,NFLEVG)=PRE5 (JROF,NFLEVG)*PSPT0(JROF)
  PRE0L(JROF       )=PSPT5L(JROF)*PRE0(JROF,NFLEVG)&
   & + PRE5(JROF,NFLEVG)*PSPT0L(JROF)
  PRE0M(JROF       )=PSPT5M(JROF)*PRE0(JROF,NFLEVG)&
   & + PRE5(JROF,NFLEVG)*PSPT0M(JROF)
ENDDO
IF (LLPRE9) THEN
  DO JROF=KST,KEND
    PRE9 (JROF,NFLEVG)=PRE95(JROF,NFLEVG)*PSPT9(JROF)
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GP_SPVTL',1,ZHOOK_HANDLE)
END SUBROUTINE GP_SPVTL
