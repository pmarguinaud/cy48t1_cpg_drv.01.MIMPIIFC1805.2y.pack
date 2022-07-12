SUBROUTINE CPCLS_ASSIM(YDGEOMETRY,KST,KND, &
 & PUCLS1, PVCLS1, PNUCLS1, PNVCLS1, PTCLS1, PHUCLS1, &
 & PXUCLS, PXVCLS, PXNUCLS, PXNVCLS, PXTCLS, PXRHCLS )

!**** *CPCLS_ASSIM* - INTERFACE FOR CLS FIELDS

!     Purpose.
!     --------
!           DIAGNOSTICS OF PHYSICAL FLUXES IN CLS ARRAYS

!**   Interface.
!     ----------
!        *CALL* *CPCLS_ASSIM*

!        Explicit arguments :
!        --------------------

!       NPROMA                 - HORIZONTAL DIMENSION                 (INPUT)
!       KST to KND             - NB OF POINTS                         (INPUT)
!       FLUXES COMING FROM THE PHYSICAL PARAMETERIZATIONS             (INPUT)

!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!      F. Taillefer
!      Original : 06/2016  from cpxfu

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KND 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUCLS1(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVCLS1(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNUCLS1(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNVCLS1(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTCLS1(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHUCLS1(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXUCLS(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXVCLS(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXNUCLS(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXNVCLS(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXTCLS(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXRHCLS(YDGEOMETRY%YRDIM%NPROMA) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPCLS_ASSIM',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

DO JROF = KST,KND
  PUCLS1(JROF)=PXUCLS(JROF)
  PVCLS1(JROF)=PXVCLS(JROF)
  PNUCLS1(JROF)=PXNUCLS(JROF)
  PNVCLS1(JROF)=PXNVCLS(JROF)
  PTCLS1(JROF)=PXTCLS(JROF)
  PHUCLS1(JROF)=MAX(0.0_JPRB,MIN(1.0_JPRB,PXRHCLS(JROF)))
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPCLS_ASSIM',1,ZHOOK_HANDLE)

END SUBROUTINE CPCLS_ASSIM
