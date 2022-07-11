SUBROUTINE GNHD3(    LD_RDRY_VD, YDCST, KFLEV, KPROMA, KSTART, KEND, PDIV, PNHX, PDVER, PR, PMAX_D3, &
& LDBOUND_D3, P3DIVG)

! GNHD3 - Diagnose D3-term in the NH model.

! Purpose
! -------
!   Diagnose D3-term in the NH model from D, NHX, dver, R.
!    D3 = D + NHX + dver (Rd/R) = D + dqua + dver (Rd/R - 1)
!   Remember that:
!    dver = - (g pre)/(Rd T) (d w/d prehyd)
!    dver (Rd/R) = - (g pre)/(RT) (d w/d prehyd)
!    dqua = dver + NHX

! Interface
! ---------
!   * INPUT:
!   KFLEV   - number of levels.
!   KPROMA  - length of work.
!   KSTART  - start of work.
!   KEND    - end of work.
!   PDIV    - Horizontal divergence D = grad V at full levels.
!   PNHX    - NHX-term at full levels.
!   PDVER   - vertical divergence at full levels.
!   PR      - R at full levels.
!   PMAX_D3 - maximum value allowed for absolute value of D3 when LDBOUND_D3=.T.
!   LDBOUND_D3 - bound absolute value of D3


!   * OUTPUT:
!   P3DIVG  - D3-term at full levels.

! Externals
! ---------

! Method
! ------

! Reference
! ---------

! Author
! ------
!   06 Dec 2004 K. Yessad (after GNHPDVD).

! Modifications
! -------------
!   M. Hortal (May 2015): Limit for increased stability
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   K. Yessad (Apr 2018): introduce key L_RDRY_VD (ensure consistent definition of "dver" everywhere).
! End Modifications
!---------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMCST   , ONLY : TCST

USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK

! -----------------------------------------------------------------------------

IMPLICIT NONE

LOGICAL            ,INTENT(IN)   :: LD_RDRY_VD
TYPE(TCST)         ,INTENT(IN)   :: YDCST
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KFLEV
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KPROMA 
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KSTART 
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KEND 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PDIV(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PNHX(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PDVER(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PR(KPROMA,KFLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PMAX_D3
LOGICAL            ,INTENT(IN)   :: LDBOUND_D3
REAL(KIND=JPRB)    ,INTENT(OUT)  :: P3DIVG(KPROMA,KFLEV) 

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHD3', 0, ZHOOK_HANDLE)

! -----------------------------------------------------------------------------

IF (LD_RDRY_VD) THEN
  ! Ratio RD/PR is due to current definition of PDVER (with RD)
  DO JLEV=1,KFLEV
    DO JROF=KSTART,KEND
      P3DIVG(JROF,JLEV)=PDIV(JROF,JLEV)+PNHX(JROF,JLEV)+PDVER(JROF,JLEV)*(YDCST%RD/PR(JROF,JLEV))
    ENDDO
  ENDDO
ELSE
  DO JLEV=1,KFLEV
    DO JROF=KSTART,KEND
      P3DIVG(JROF,JLEV)=PDIV(JROF,JLEV)+PNHX(JROF,JLEV)+PDVER(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

IF (LDBOUND_D3) THEN
  DO JLEV=1,KFLEV
    DO JROF=KSTART,KEND
      P3DIVG(JROF,JLEV)=SIGN(MIN(ABS(P3DIVG(JROF,JLEV)),PMAX_D3),P3DIVG(JROF,JLEV))
    ENDDO
  ENDDO
ENDIF

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHD3', 1, ZHOOK_HANDLE)

END SUBROUTINE GNHD3

