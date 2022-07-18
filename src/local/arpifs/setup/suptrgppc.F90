SUBROUTINE SUPTRGPPC(YDDIMV,YGFL,LDSLPHY,KGPPC,&
 & KGPPCF_U ,KGPPCF_V,KGPPCF_T,KGPPCF_PD,KGPPCF_VD,KGPPCF_SP,KGPPCF_CP,&
 & KGPPCF_NHX,KGPPCF_UP,KGPPCF_VP,KGPPCF_TP,KGPPCF_GFLP,KGPPCF_BBC,&
 & KFGPPC)

! SUPTRGPPC - Pointers of PTRGPPC.
!             See PTRGPPC for comments of dummy arguments.
!             KGPPCF_.. are used for LPC_CHEAP.
!             KFGPPC is the number of 2D fields contained in GPPCBUF.
!             (matches with NFGPPC of PTRGPPC).

!             These pointers must later to be put in GMV+GMVS (new attributes
!             to create), and also GFL for KGPPCF_GFLP.
!             They are linked to the different PC schemes and they are no
!             longer purely NH pointers, so they have been renamed in Jan 2008.

! Purpose
! -------

! Interface
! ---------

! Externals
! ---------

! Method
! ------
!   Pointers for GPPCBUF buffer

! Reference
! ---------

! Author
! ------
!      Nov-2001 J. VIVODA (with name ald/adiab/ESTR_POINTERS).

! Modifications
! -------------
!   08-Jan-2008 N. Wedi and K. Yessad: different dev for NH model and PC scheme
!   23-Sep-2008 K. Yessad: remove lpc_xidt
!   K. Yessad (Nov 2009): prune lpc_old.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (July 2014): Move some variables.
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Introduce NHQE model.
! End Modifications
!------------------------------------------------------------------------------

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMCT0   , ONLY : LNHDYN, NUNDEFLD, LNHEE
USE YOMDYNA  , ONLY : YRDYNA
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TYPE_GFLD)   ,INTENT(INOUT) :: YGFL
LOGICAL           ,INTENT(IN)    :: LDSLPHY
INTEGER(KIND=JPIM),INTENT(OUT)   :: KGPPC
INTEGER(KIND=JPIM),INTENT(OUT)   :: KGPPCF_U
INTEGER(KIND=JPIM),INTENT(OUT)   :: KGPPCF_V
INTEGER(KIND=JPIM),INTENT(OUT)   :: KGPPCF_T
INTEGER(KIND=JPIM),INTENT(OUT)   :: KGPPCF_PD
INTEGER(KIND=JPIM),INTENT(OUT)   :: KGPPCF_VD
INTEGER(KIND=JPIM),INTENT(OUT)   :: KGPPCF_SP
INTEGER(KIND=JPIM),INTENT(OUT)   :: KGPPCF_CP
INTEGER(KIND=JPIM),INTENT(OUT)   :: KGPPCF_NHX
INTEGER(KIND=JPIM),INTENT(OUT)   :: KGPPCF_UP
INTEGER(KIND=JPIM),INTENT(OUT)   :: KGPPCF_VP
INTEGER(KIND=JPIM),INTENT(OUT)   :: KGPPCF_TP
INTEGER(KIND=JPIM),INTENT(OUT)   :: KGPPCF_GFLP
INTEGER(KIND=JPIM),INTENT(OUT)   :: KGPPCF_BBC
INTEGER(KIND=JPIM),INTENT(OUT)   :: KFGPPC

!------------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IDEFAULT, ITMP

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUPTRGPPC',0,ZHOOK_HANDLE)
ASSOCIATE(NUMFLDS=>YGFL%NUMFLDS, &
 & NFLEVG=>YDDIMV%NFLEVG)
!------------------------------------------------------------------------------

IDEFAULT = NUNDEFLD

! 0. Set default values
KGPPC           = IDEFAULT
KGPPCF_U        = IDEFAULT
KGPPCF_V        = IDEFAULT
KGPPCF_T        = IDEFAULT
KGPPCF_PD       = IDEFAULT
KGPPCF_VD       = IDEFAULT
KGPPCF_SP       = IDEFAULT
KGPPCF_CP       = IDEFAULT
KGPPCF_NHX      = IDEFAULT
KGPPCF_UP       = IDEFAULT
KGPPCF_VP       = IDEFAULT
KGPPCF_TP       = IDEFAULT
KGPPCF_GFLP     = IDEFAULT
KGPPCF_BBC      = IDEFAULT

! 1. Set_pointers

KGPPC= 1
ITMP = 1
KFGPPC=0

IF( YRDYNA%LPC_CHEAP ) THEN

  KGPPCF_U = ITMP;    ITMP = ITMP + NFLEVG
  KGPPCF_V = ITMP;    ITMP = ITMP + NFLEVG
  KGPPCF_T = ITMP;    ITMP = ITMP + NFLEVG
  KGPPCF_SP= ITMP;    ITMP = ITMP + NFLEVG
  KGPPCF_CP= ITMP;    ITMP = ITMP + NFLEVG
  KFGPPC = KFGPPC+5*NFLEVG

  IF( LDSLPHY ) THEN
    ! GFL provisionally put here (to be put in GFL(PC) later).
    KGPPCF_UP  = ITMP;    ITMP = ITMP + NFLEVG
    KGPPCF_VP  = ITMP;    ITMP = ITMP + NFLEVG
    KGPPCF_TP  = ITMP;    ITMP = ITMP + NFLEVG
    KGPPCF_GFLP= ITMP;    ITMP = ITMP + NUMFLDS*NFLEVG
    KFGPPC = KFGPPC+3*NFLEVG+NUMFLDS*NFLEVG
  ENDIF

  IF( LNHEE )THEN
    KGPPCF_PD = ITMP; ITMP = ITMP + NFLEVG
    KFGPPC = KFGPPC+NFLEVG
  ENDIF
  IF( LNHDYN )THEN
    KGPPCF_VD = ITMP; ITMP = ITMP + NFLEVG
    KFGPPC = KFGPPC+NFLEVG
    IF (YRDYNA%LNHX) THEN
      KGPPCF_NHX  = ITMP; ITMP = ITMP + NFLEVG
      KFGPPC = KFGPPC+NFLEVG
    ENDIF
    IF( YRDYNA%LRDBBC ) THEN
      KGPPCF_BBC = ITMP; ITMP = ITMP + 3
      KFGPPC = KFGPPC+3      
    ENDIF
  ENDIF

ENDIF

!------------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUPTRGPPC',1,ZHOOK_HANDLE)
END SUBROUTINE SUPTRGPPC

