SUBROUTINE GTOPTX2(KLAT,KLON,KJJ1,PREL)

!     PURPOSE
!     -------

!     This routine calculates the component 11 of the topographic gradient
!     correlation tensor(TGCT)-part one

!     INTERFACE
!     ---------

!     CALL GTOPTX2(...)
!        KLAT-number of latitudes
!        KLON-number of longitudes
!        PREL-height of orography (component 11 of symmetric TGCT)

!     METHOD
!     ------

!     The derivatives are calculated using the values from two neighbouring
!     points : LON+N and LON-N for X-derivatives (N is such that the step 
!     is isotropic) 

!     EXTERNALS : None
!     ---------

!     AUTHORS
!     -------
!      E. Cordoneanu 96-12-15 from GTOPT

!     MODIFICATIONS
!     -------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : RPI      ,RA
USE YOMCLI   , ONLY : YRCLI

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLAT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KJJ1 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PREL(KLON,KLAT) 
REAL(KIND=JPRB) ::ZPX2(KLON)

INTEGER(KIND=JPIM) :: IIM, IIP, IX, IXM, IXP, J, JI, JJ

REAL(KIND=JPRB) :: ZC, ZIX, ZSIN, ZX, ZY
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('GTOPTX2',0,ZHOOK_HANDLE)
ZC=RPI/REAL(YRCLI%NGLOBY,JPRB)
ZY=2.0_JPRB*RA*ZC
IF(YRCLI%LGLOBE) THEN
  DO JJ=1,KLAT
    ZSIN=SIN((JJ-0.5_JPRB)*ZC)
    ZIX=1.0_JPRB/ZSIN
    IX=MIN(KLON/4,INT(ZIX))
    ZX=ZY*ZSIN*IX
    IIP=IX-1
    IIM=KLON-IX-1
    DO J=1,KLON
      IXP=1+MOD(J+IIP,KLON)
      IXM=1+MOD(J+IIM,KLON)
      ZPX2(J)=(PREL(IXP,JJ)-PREL(IXM,JJ))/ZX
      ZPX2(J)=ZPX2(J)**2
    ENDDO
    DO J=1,KLON
      PREL(J,JJ)=ZPX2(J)
    ENDDO
  ENDDO
ELSE
  DO JJ=1,KLAT
    ZSIN=SIN((JJ+(KJJ1-1)-0.5_JPRB)*ZC)
    ZIX=1.0_JPRB/ZSIN
    IX=MIN(YRCLI%NGLOBX/4,INT(ZIX))
    ZX=2*RPI*RA*ZSIN/YRCLI%NGLOBX
    DO JI=1,KLON
      IXP=MIN(JI+IX,KLON)
      IXM=MAX(JI-IX,1)
      ZPX2(JI)=(PREL(IXP,JJ)-PREL(IXM,JJ))/(ZX*(IXP-IXM))
      ZPX2(JI)=ZPX2(JI)*ZPX2(JI)
    ENDDO
    DO J=1,KLON
      PREL(J,JJ)=ZPX2(J)
    ENDDO
  ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('GTOPTX2',1,ZHOOK_HANDLE)
END SUBROUTINE GTOPTX2
