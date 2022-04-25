SUBROUTINE LOCMAXI(KLAT,KLON,PREL)

!     PURPOSE.
!     --------

!     This routine calculates the occurrence of a local maximum in a field
!     Simplified version.

!**   INTERFACE.
!     ----------

!     CALL LOCMAXI(...)
!        KLAT  = number of latitudes
!        KLON  = number of longitudes
!        PREL  = input/output field
!        output is 0 if not a maximum, 1 if a maximum

!     METHOD.
!     -------

!     The value of the field is compared with the 8 surrounding values
!     Maxima in the last and first latitude rows are not detected.

!     EXTERNALS : None.
!     ----------

!     AUTHORS.
!     --------
!     E. Cordoneanu 96-12-15 and D. Giard 97-03-13 from LOCMAX

!     MODIFICATIONS.
!     --------------
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCLI   , ONLY : YRCLI

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLAT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PREL(KLON,KLAT) 

REAL(KIND=JPRB) :: ZMAX(KLON,3)
REAL(KIND=JPRB) :: ZUMZEPS, ZEPS
REAL(KIND=JPRB) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: IC, IIC, IIM, IIP, IJI1, IJI2, IJJ1, IJJ2, IM, IP, J, JJ


IF (LHOOK) CALL DR_HOOK('LOCMAXI',0,ZHOOK_HANDLE)

! ZEPS controls a very slight recursive smoothing which is applied to get
! a number of peaks different from zero when the orography is duplicated.
! With ZEPS=1.E-6 and a 64 bit machine, the program is able to count the
! maxima even when the values are repeated four times in both directions.
! But with five times or more, ZEPS must be increased since 1.E-18 + 1.=1.
! The degradation of the data is negligible as long as ZEPS is less than
! 1.E-3 (a few meters of error on the mountain elevation)

ZEPS=1.E-3_JPRB
ZUMZEPS=1.0_JPRB-4._JPRB*ZEPS

DO JJ=1,3
  DO J=1,KLON
    ZMAX(J,JJ)=0.0_JPRB
  ENDDO
ENDDO
IF(YRCLI%LGLOBE) THEN
  IJJ1=1
  IJJ2=KLAT
  IJI1=1
  IJI2=KLON
ELSE
  IJJ1=2
  IJJ2=KLAT-1
  IJI1=2
  IJI2=KLON-1
ENDIF

DO JJ=IJJ1,IJJ2

  IP=MIN(JJ+1,KLAT)
  IC=    JJ
  IM=MAX(JJ-1,   1)
!DIR$ IVDEP
!OCL NOVREC
  DO J=IJI1,IJI2
    IIP=1+MOD(J       ,KLON)
    IIC=J
    IIM=1+MOD(J-2+KLON,KLON)
    PREL(IIC,IC)= PREL(IIC,IC)*ZUMZEPS &
     & + (PREL(IIC,IP)+PREL(IIC,IM)+ PREL(IIP,IC)+PREL(IIM,IC))*ZEPS  
  ENDDO

  DO J=IJI1,IJI2
    IIP=1+MOD(J       ,KLON)
    IIC=J
    IIM=1+MOD(J-2+KLON,KLON)
    IF (MAX(PREL(IIM,IM),PREL(IIC,IM),PREL(IIP,IM),&
       & PREL(IIM,IC),             PREL(IIP,IC),&
       & PREL(IIM,IP),PREL(IIC,IP),PREL(IIP,IP))&
       & <              PREL(IIC,IC)) THEN  
      ZMAX(J,3)=1.0_JPRB
    ELSE
      ZMAX(J,3)=0.0_JPRB
    ENDIF
  ENDDO

  IF (JJ > 2) THEN
    DO J=1,KLON
      PREL(J,JJ-2)=ZMAX(J,1)
    ENDDO
  ENDIF

  DO J=1,KLON
    ZMAX(J,1)=ZMAX(J,2)
    ZMAX(J,2)=ZMAX(J,3)
  ENDDO

ENDDO

DO J=1,KLON
  PREL(J,KLAT-1)=ZMAX(J,1)
  PREL(J,KLAT  )=ZMAX(J,2)
ENDDO

IF (LHOOK) CALL DR_HOOK('LOCMAXI',1,ZHOOK_HANDLE)
END SUBROUTINE LOCMAXI
