SUBROUTINE CP_PTRSLB1(YDDYN,YDPTRSLB1,KSLB1U9,KSLB1V9,KSLB1T9,KSLB1VD9,KSLB1GFL9)

!**** *CP_PTRSLB1*

!     Purpose.
!     --------
!      Compute local versions of pointers for buffer "SLBUF1"

!**   Interface.
!     ----------
!        *CALL* *CP_PTRSLB1(...)*

!        Explicit arguments :
!        --------------------
!         KSLB1[X]9: output pointer for variable [X]   (output)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!      K. Yessad after code in CPG (Feb 2011)

!     Modifications
!     -------------
!     End Modifications
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMDYN   , ONLY : TDYN
USE PTRSLB1  , ONLY : TPTRSLB1

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TPTRSLB1)    ,INTENT(IN)    :: YDPTRSLB1
INTEGER(KIND=JPIM),INTENT(OUT)   :: KSLB1U9
INTEGER(KIND=JPIM),INTENT(OUT)   :: KSLB1V9
INTEGER(KIND=JPIM),INTENT(OUT)   :: KSLB1T9
INTEGER(KIND=JPIM),INTENT(OUT)   :: KSLB1VD9
INTEGER(KIND=JPIM),INTENT(OUT)   :: KSLB1GFL9

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CP_PTRSLB1',0,ZHOOK_HANDLE)
ASSOCIATE(LSPLTHOIGFL=>YDDYN%LSPLTHOIGFL, NSPLTHOI=>YDDYN%NSPLTHOI, &
 & NSVDLAG=>YDDYN%NSVDLAG, NTLAG=>YDDYN%NTLAG, NWLAG=>YDDYN%NWLAG, &
 & MSLB1GFL9=>YDPTRSLB1%MSLB1GFL9, MSLB1GFLF9=>YDPTRSLB1%MSLB1GFLF9, &
 & MSLB1T0=>YDPTRSLB1%MSLB1T0, MSLB1T9=>YDPTRSLB1%MSLB1T9, &
 & MSLB1TF9=>YDPTRSLB1%MSLB1TF9, MSLB1U0=>YDPTRSLB1%MSLB1U0, &
 & MSLB1U9=>YDPTRSLB1%MSLB1U9, MSLB1UF9=>YDPTRSLB1%MSLB1UF9, &
 & MSLB1V0=>YDPTRSLB1%MSLB1V0, MSLB1V9=>YDPTRSLB1%MSLB1V9, &
 & MSLB1VD0=>YDPTRSLB1%MSLB1VD0, MSLB1VD9=>YDPTRSLB1%MSLB1VD9, &
 & MSLB1VDF9=>YDPTRSLB1%MSLB1VDF9, MSLB1VF9=>YDPTRSLB1%MSLB1VF9)
!     ------------------------------------------------------------------

! * U,V components
IF (NWLAG == 4) THEN
  KSLB1U9  = MSLB1U0
  KSLB1V9  = MSLB1V0
ELSE
  IF (NSPLTHOI /= 0) THEN
    KSLB1U9  = MSLB1UF9
    KSLB1V9  = MSLB1VF9
  ELSE
    KSLB1U9  = MSLB1U9
    KSLB1V9  = MSLB1V9
  ENDIF
ENDIF
! * temperature
IF (NTLAG == 4) THEN
  KSLB1T9  = MSLB1T0
ELSE
  IF (NSPLTHOI /= 0) THEN
    KSLB1T9  = MSLB1TF9
  ELSE
    KSLB1T9  = MSLB1T9
  ENDIF
ENDIF
! * vertical divergence
IF (NSVDLAG == 4) THEN
  KSLB1VD9 = MSLB1VD0
ELSE
  IF (NSPLTHOI /= 0) THEN
    KSLB1VD9 = MSLB1VDF9
  ELSE
    KSLB1VD9 = MSLB1VD9
  ENDIF
ENDIF
! * GFL arrays
IF (LSPLTHOIGFL.OR.(NSPLTHOI /= 0)) THEN
  KSLB1GFL9= MSLB1GFLF9
ELSE
  KSLB1GFL9= MSLB1GFL9
ENDIF  

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CP_PTRSLB1',1,ZHOOK_HANDLE)
END SUBROUTINE CP_PTRSLB1
