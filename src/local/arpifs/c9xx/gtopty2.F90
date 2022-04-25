SUBROUTINE GTOPTY2(KLAT,KLON,PREL)

!     PURPOSE
!     -------

!     This routine calculates the component 22 of the topographic gradient
!     correlation tensor(TGCT)-part one

!     INTERFACE
!     ---------

!     CALL GTOPTY2(...)
!        KLAT-number of latitudes
!        KLON-number of longitudes
!        PREL-height of orography (component 22 of symmetric TGCT)

!     METHOD
!     ------

!     The derivatives are calculated using the values from two neighbouring
!     points : LAT+1 and LAT-1 for Y-derivatives

!     EXTERNALS : None
!     ---------

!     AUTHORS
!     -------
!      E. Cordoneanu 96-12-15 from GTOPT

!     MODIFICATIONS
!     ------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : RPI      ,RA
USE YOMCLI   , ONLY : YRCLI

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLAT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PREL(KLON,KLAT) 
REAL(KIND=JPRB) ::ZPY2(KLON,2)

INTEGER(KIND=JPIM) :: IL, J, JI, JJ

REAL(KIND=JPRB) :: ZC, ZY
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('GTOPTY2',0,ZHOOK_HANDLE)
ZC=RPI/REAL(YRCLI%NGLOBY,JPRB)
ZY=2.0_JPRB*RA*ZC

!  First latitude circle
IF (YRCLI%LGLOBE) THEN
  DO J=1,KLON
    IL=1+MOD(J-1+KLON/2,KLON)
    ZPY2(J,1)=(PREL(IL,1)-PREL(J,2))/ZY
    ZPY2(J,1)=ZPY2(J,1)**2
  ENDDO
ELSE
  DO JI=1,KLON
    ZPY2(JI,1)=(PREL(JI,1)-PREL(JI,2))/ZY/2.0_JPRB
    ZPY2(JI,1)=ZPY2(JI,1)**2
  ENDDO
ENDIF

!  Internal latitude circles

DO JJ=2,KLAT-1
  DO J=1,KLON
    ZPY2(J,2)=(PREL(J,JJ-1)-PREL(J,JJ+1))/ZY
    ZPY2(J,2)=ZPY2(J,2)**2
  ENDDO
  DO J=1,KLON
    PREL(J,JJ-1)=ZPY2(J,1)
    ZPY2(J,1)=ZPY2(J,2)
  ENDDO
ENDDO

!  Last latitude circle
IF (YRCLI%LGLOBE) THEN
  DO J=1,KLON
    IL=1+MOD(J-1+KLON/2,KLON)
    ZPY2(J,2)=(PREL(J,KLAT-1)-PREL(IL,KLAT))/ZY
    ZPY2(J,2)=ZPY2(J,2)**2
  ENDDO
ELSE
  DO J=1,KLON
    ZPY2(J,2)=(PREL(J,KLAT-1)-PREL(J,KLAT))/ZY/2.0_JPRB
    ZPY2(J,2)=ZPY2(J,2)**2
  ENDDO
ENDIF

DO J=1,KLON
  PREL(J,KLAT-1)=ZPY2(J,1)
  PREL(J,KLAT  )=ZPY2(J,2)
ENDDO

IF (LHOOK) CALL DR_HOOK('GTOPTY2',1,ZHOOK_HANDLE)
END SUBROUTINE GTOPTY2
