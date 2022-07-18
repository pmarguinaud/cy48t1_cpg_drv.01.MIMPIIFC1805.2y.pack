SUBROUTINE CPG_ZERO_AD(&
 ! ---------------------------------------------------------------------------
 ! - INPUT:
 & YDGEOMETRY, YDDYN,YDPTRGPPC,YDSIMPHL,KST,KEND,&
 ! ---------------------------------------------------------------------------
 ! - OUTPUT:
 & PRT0L,PRT0M,PRE0F,PRE0,PRE0L,PRE0M,PRCP0,PHIF0,PHI0,PKENE0,&
 & PXYB0,PCTY0,PUVH0,&
 & PRE9F,PRE9,PRCP9,PHIF9,PHI9,PXYB9,&
 & PGPPC,PATND)

!**** *CPG_ZERO_AD* - Adjoint grid point calculations.
!                     Set to zero some sensitivities.

!     Purpose.
!     --------
!           Adjoint grid point calculations.
!           Set to zero some sensitivities.

!           Notice that this routine has no TL counterpart.
!           These settings to zero are not done in CPG_END_AD in order
!            to keep the maximum of consistency between CPG_END_AD and
!            CPG_END_TL.

!**   Interface.
!     ----------
!        *CALL* *CPG_ZERO_AD*

!        Explicit arguments :
!        --------------------

!*    INPUT:
!     ------
!        KST       : first element of work.
!        KEND      : last element of work.

!*    OUTPUT:
!     -------
!        PRT0L to PATND: see GPG_END_AD, CPG_DYN_AD or CPG_GP_AD for comments.

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
!        K. Yessad (Jul 2004) after part 7 of CPGAD.

!     Modifications.
!     --------------
!   Original : 12-Jul-2004
!   K. Yessad: 10-Apr-2006 Move some zeroing from CPG_END_AD.
!   S. Ivatek-S: 17-Apr-2007 Over dimensioning of PGPNH to NFGPNH+1,
!                           boundary checking problem if NFGPNH=0 bf
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2008): remove dummy CDLOCK
!   K. Yessad (Aug 2009): remove LPC_OLD in TL and AD codes.
!   K. Yessad (Nov 2009): cleanings, DT/Dt now pre-computed in CPG_GP_AD.
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   K. Yessad (July 2014): Move some variables.
! End Modifications
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCT0   , ONLY : LNHDYN, LTWOTL
USE YOMDYNA  , ONLY : YRDYNA
USE PTRGPPC  , ONLY : TPTRGPPC
USE YOMDYN   , ONLY : TDYN
USE YOMSIMPHL, ONLY : TSIMPHL
USE INTDYN_MOD,ONLY : YYTTND, YYTHW0, YYTCTY0, YYTRCP0, YYTRCP9, YYTXYB0, YYTXYB9

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(INOUT) :: YDDYN
TYPE(TPTRGPPC)    ,INTENT(INOUT) :: YDPTRGPPC
TYPE(TSIMPHL)     ,INTENT(INOUT) :: YDSIMPHL
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0L(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0M(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRCP0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP0%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHIF0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHI0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKENE0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXYB0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB0%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCTY0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTCTY0%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUVH0(YDGEOMETRY%YRDIM%NPROMNH,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTHW0%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE9F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRCP9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP9%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHIF9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHI9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXYB9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB9%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGPPC(YDGEOMETRY%YRDIM%NPROMA,YDPTRGPPC%NFGPPC+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PATND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTTND%NDIM)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPG_ZERO_AD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, NPROMNH=>YDDIM%NPROMNH, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LSIMPH=>YDSIMPHL%LSIMPH, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
 & NFGPPC=>YDPTRGPPC%NFGPPC)
!     ------------------------------------------------------------------

!*       1.    Set local sensitivities to zero for the start for var at t
!              ----------------------------------------------------------

PRT0L(KST:KEND,1:NFLEVG)=0.0_JPRB
PRT0M(KST:KEND,1:NFLEVG)=0.0_JPRB
PRE0F(KST:KEND,1:NFLEVG)=0.0_JPRB
PHIF0(KST:KEND,1:NFLEVG)=0.0_JPRB
PKENE0(KST:KEND,1:NFLEVG)=0.0_JPRB
PRE0(KST:KEND,0:NFLEVG)=0.0_JPRB
PHI0(KST:KEND,0:NFLEVG)=0.0_JPRB
PXYB0(KST:KEND,1:NFLEVG,1:YYTXYB0%NDIM)=0.0_JPRB
PCTY0(KST:KEND,0:NFLEVG,1:YYTCTY0%NDIM)=0.0_JPRB
PRCP0(KST:KEND,1:NFLEVG,1:YYTRCP0%NDIM)=0.0_JPRB

IF (LNHDYN) THEN
  PUVH0(KST:KEND,0:NFLEVG,1:YYTHW0%NDIM)=0.0_JPRB
ENDIF

PGPPC(KST:KEND,1:NFGPPC+1)=0.0_JPRB
PRE0L(KST:KEND)=0.0_JPRB
PRE0M(KST:KEND)=0.0_JPRB

PATND(KST:KEND,1:NFLEVG,1:YYTTND%NDIM)=0.0_JPRB

!     ------------------------------------------------------------------

!*       2.    Set local sensitivities to zero for the start for var at t-dt
!              -------------------------------------------------------------

IF ((.NOT.LTWOTL) .AND. (NCURRENT_ITER == 0) .AND. (LSIMPH.OR.LNHDYN) ) THEN
  PRE9F(KST:KEND,1:NFLEVG)=0.0_JPRB
  PXYB9(KST:KEND,1:NFLEVG,1:YYTXYB9%NDIM)=0.0_JPRB
ENDIF

IF ((.NOT.LTWOTL).AND.LSIMPH) THEN
  PRCP9(KST:KEND,1:NFLEVG,1:YYTRCP9%NDIM)=0.0_JPRB
  PHI9(KST:KEND,0:NFLEVG)=0.0_JPRB
  PHIF9(KST:KEND,1:NFLEVG)=0.0_JPRB
ENDIF

IF ((LTWOTL.AND.YRDYNA%LPC_FULL).OR.(.NOT.LTWOTL)) THEN
  PRE9(KST:KEND,0:NFLEVG)=0.0_JPRB
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG_ZERO_AD',1,ZHOOK_HANDLE)
END SUBROUTINE CPG_ZERO_AD
