SUBROUTINE GTOPTXY(KLAT,KLON,KJJ1,KJJ2,PREL)

!     PURPOSE
!     -------

!     This routine calculates the component 12 of the topographic gradient
!     correlation tensor(TGCT)-part one

!     INTERFACE
!     ---------

!     CALL GTOPTXY(...)
!        KLAT-number of latitudes
!        KLON-number of longitudes
!        PREL-height of orography (component 12 of symmetric TGCT)

!     METHOD
!     ------

!     The derivatives are calculated using the values from two neighbouring
!     points : LAT+1 and LAT-1 for Y-derivatives, LON+N and LON-N for
!     X-derivatives (N is such that the step is isotropic)

!     EXTERNALS : None
!     ---------

!     AUTHORS
!     -------
!      E. Cordoneanu 96-12-15 from GTOPT

!     MODIFICATIONS:
!     --------------
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
INTEGER(KIND=JPIM),INTENT(IN)    :: KJJ2 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PREL(KLON,KLAT) 
REAL(KIND=JPRB) ::ZPXY(KLON,2)

INTEGER(KIND=JPIM) :: IIM, IIMPOL, IIP, IIPPOL, IX, IXL, IXM, IXP, IXPOL, J, JI, JJ

REAL(KIND=JPRB) :: ZC, ZIX, ZIXPO, ZSIN, ZSPOL, ZX, ZXPOL, ZY
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('GTOPTXY',0,ZHOOK_HANDLE)
ZC=RPI/REAL(YRCLI%NGLOBY,JPRB)
ZY=2.0_JPRB*RA*ZC
IF(YRCLI%LGLOBE) THEN
  ZSPOL=SIN(0.5_JPRB*ZC)
  ZIXPO=1.0_JPRB/ZSPOL
  IXPOL=MIN(KLON/4,INT(ZIXPO))
  ZXPOL=ZY*ZSPOL*IXPOL
  IIPPOL=IXPOL-1
  IIMPOL=KLON-IXPOL-1

!  First latitude circle

  ZX=ZXPOL
  DO J=1,KLON
    IXP=1+MOD(J+IIPPOL  ,KLON)
    IXM=1+MOD(J+IIMPOL  ,KLON)
    IXL=1+MOD(J-1+KLON/2,KLON)
    ZPXY(J,1)=(PREL(IXP,1)-PREL(IXM,1))/ZX*(PREL(IXL,1)-PREL(  J,2))/ZY
  ENDDO

!  Internal latitude circles

  DO JJ=2,KLAT-1
    ZSIN=SIN((JJ-0.5_JPRB)*ZC)
    ZIX=1.0_JPRB/ZSIN
    IX=MIN(KLON/4,INT(ZIX))
    ZX=ZY*ZSIN*IX
    IIP=IX-1
    IIM=KLON-IX-1
    DO J=1,KLON
      IXP=1+MOD(J+IIP,KLON)
      IXM=1+MOD(J+IIM,KLON)
      ZPXY(J,2)=(PREL(IXP,JJ  )-PREL(IXM,JJ  ))/ZX &
       & *(PREL(  J,JJ-1)-PREL(  J,JJ+1))/ZY  
    ENDDO
    DO J=1,KLON
      PREL(J,JJ-1)=ZPXY(J,1)
      ZPXY(J,1)=ZPXY(J,2)
    ENDDO
  ENDDO

!  Last latitude circle

  ZX=ZXPOL
  DO J=1,KLON
    IXP=1+MOD(J+IIPPOL,KLON)
    IXM=1+MOD(J+IIMPOL,KLON)
    IXL=1+MOD(J-1+KLON/2,KLON)
    ZPXY(J,2)=(PREL(IXP,KLAT  )-PREL(IXM,KLAT))/ZX &
     & *(PREL(  J,KLAT-1)-PREL(IXL,KLAT))/ZY  
  ENDDO
ELSE
  ZSIN=SIN((KJJ1-0.5_JPRB)*ZC)
  IX=INT(1.0_JPRB/ZSIN)
  IX=MIN(YRCLI%NGLOBX/4,IX)
  ZX=2*RPI*RA*ZSIN/YRCLI%NGLOBX
  DO JI=1,KLON
    IXP=MIN((JI+IX),KLON)
    IXM=MAX((JI-IX),1)
    ZPXY(JI,1)=(PREL(IXP,1)-PREL(IXM,1))/(ZX*(IXP-IXM))*&
     & (PREL(JI,1)-PREL(JI,2))/ZY  
  ENDDO
!     Internal latitude circles

  DO JJ=2,KLAT-1
    ZSIN=SIN((JJ+(KJJ1-1)-0.5_JPRB)*ZC)
    IX=INT(1.0_JPRB/ZSIN)
    IX=MIN(YRCLI%NGLOBX/4,IX)
    ZX=2*RPI*RA*ZSIN/YRCLI%NGLOBX
    DO JI=1,KLON
      IXP=MIN(JI+IX,KLON)
      IXM=MAX(JI-IX,1)
      ZPXY(JI,2)=(PREL(IXP,JJ)-PREL(IXM,JJ))/(ZX*(IXP-IXM))*&
       & (PREL(JI,JJ-1)-PREL(JI,JJ+1))/ZY  
    ENDDO
    DO JI=1,KLON
      PREL(JI,JJ-1)=ZPXY(JI,1)
      ZPXY(JI,1)=ZPXY(JI,2)
    ENDDO
  ENDDO
!     Last latitude circle

  ZSIN=SIN((KJJ2-0.5_JPRB)*ZC)
  IX=INT(1.0_JPRB/ZSIN)
  IX=MIN(YRCLI%NGLOBX/4,IX)
  ZX=2*RPI*RA*ZSIN/YRCLI%NGLOBX
  DO JI=1,KLON
    IXP=MIN((JI+IX),KLON)
    IXM=MAX((JI-IX),1)
    ZPXY(JI,2)= (PREL(IXP,KLAT)-PREL(IXM,KLAT))/(ZX*(IXP-IXM))*&
     & (PREL(JI,KLAT-1)-PREL(JI,KLAT))/ZY  
  ENDDO
ENDIF

DO J=1,KLON
  PREL(J,KLAT-1)=ZPXY(J,1)
  PREL(J,KLAT)=ZPXY(J,2)
ENDDO

IF (LHOOK) CALL DR_HOOK('GTOPTXY',1,ZHOOK_HANDLE)
END SUBROUTINE GTOPTXY
