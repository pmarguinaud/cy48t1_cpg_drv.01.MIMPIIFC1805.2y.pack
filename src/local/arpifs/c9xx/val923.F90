SUBROUTINE VAL923(LDNEW)

!**** *GEO923*

!     PURPOSE.
!     --------
!      Compute the constants (YOMCLI) which are used by configuration 923.

!     INTERFACE.
!     ----------
!      CALL VAL923(...)
!        LDNEW = .FALSE. if old fields required
!        Results in YOMCLI.

!     AUTHORS.
!     --------
!      D. Giard 97-05-06

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
! ---------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCLI   , ONLY : YRCLI  
USE YOMLUN   , ONLY : NULOUT

IMPLICIT NONE

LOGICAL           ,INTENT(IN)    :: LDNEW 
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('VAL923',0,ZHOOK_HANDLE)

!  Threshold defining the mask
YRCLI%SMASK= 0.5_JPRB
!  Value for missing data + 1
YRCLI%SMANQ=-9998._JPRB
!  Land-use types for sea, ice-cap, desert, lakes
YRCLI%NTPMER= 1
YRCLI%NTPGLA= 2
YRCLI%NTPDES= 3
YRCLI%NTPLAC= 5
!  Roughness length : minimum, sea, sea-ice, urban areas, desert
YRCLI%SZZ0N= 0.001_JPRB
YRCLI%SZZ0M= 0.001_JPRB
YRCLI%SZZ0B= 0.001_JPRB
YRCLI%SZZ0U= 2.500_JPRB
YRCLI%SZZ0D= 0.001_JPRB
!  Ration of thermal to kinetic roughness length
YRCLI%STHER= 0.10_JPRB
!  Albedo : minimum, maximum, sea, ice-cap, sea-ice, desert
IF (LDNEW) THEN
  YRCLI%SALBN= 0.05_JPRB
  YRCLI%SALBX= 0.80_JPRB
ELSE
  YRCLI%SALBN= 0.07_JPRB
  YRCLI%SALBX= 0.70_JPRB
ENDIF
YRCLI%SALBM= 0.07_JPRB
YRCLI%SALBG= 0.75_JPRB
YRCLI%SALBB= 0.65_JPRB
YRCLI%SALBD= 0.10_JPRB
!  Emissivity : minimum, maximum, sea, ice-cap, sea-ice, desert
YRCLI%SEMIN= 0.90_JPRB
YRCLI%SEMIX= 1.00_JPRB
YRCLI%SEMIM= 0.96_JPRB
YRCLI%SEMIG= 0.98_JPRB
YRCLI%SEMIB= 0.97_JPRB
YRCLI%SEMID= 0.943_JPRB
!  Soil depth : minimum, maximum, desert
YRCLI%SDEPN= 0.10_JPRB
YRCLI%SDEPX= 8.00_JPRB
YRCLI%SDEPD= 0.10_JPRB
!  Percentage of clay : minimum, maximum, desert
YRCLI%SARGN=  3._JPRB
YRCLI%SARGX= 58._JPRB
YRCLI%SARGD=  3._JPRB
!  Percentage of sand : minimum, maximum, desert
YRCLI%SSABN=  6._JPRB
YRCLI%SSABX= 92._JPRB
YRCLI%SSABD= 92._JPRB
!  Minimum surface resistance : minimum, maximum, desert
YRCLI%SRSMX=5000._JPRB
YRCLI%SRSMN=   1.0_JPRB
YRCLI%SRSMD=5000._JPRB

WRITE(UNIT=NULOUT,FMT=111) YRCLI%SMASK,YRCLI%SMANQ,YRCLI%STHER,&
 & YRCLI%NTPMER,YRCLI%NTPGLA,YRCLI%NTPDES,YRCLI%NTPLAC  
WRITE(UNIT=NULOUT,FMT=112) YRCLI%SZZ0N,YRCLI%SZZ0M,YRCLI%SZZ0B,YRCLI%SZZ0U,YRCLI%SZZ0D
WRITE(UNIT=NULOUT,FMT=113) YRCLI%SALBN,YRCLI%SALBX,YRCLI%SALBM,YRCLI%SALBG,YRCLI%SALBB,YRCLI%SALBD,&
 & YRCLI%SEMIN,YRCLI%SEMIX,YRCLI%SEMIM,YRCLI%SEMIG,YRCLI%SEMIB,YRCLI%SEMID  
WRITE(UNIT=NULOUT,FMT=114) YRCLI%SDEPN,YRCLI%SDEPX,YRCLI%SDEPD,YRCLI%SARGN,YRCLI%SARGX,YRCLI%SARGD,&
 & YRCLI%SSABN,YRCLI%SSABX,YRCLI%SSABD,YRCLI%SRSMN,YRCLI%SRSMX,YRCLI%SRSMD  
111 FORMAT(' COMMON YOMCLI',/,&
 & ' SMASK=',F4.2,' SMANQ=',F6.0,' STHER=',F4.2,/&
 & ' NTPMER=',I2,' NTPGLA=',I2,' NTPDES=',I2,' NTPLAC=',I2)  
112 FORMAT(' LONGUEUR DE RUGOSITE :',/,&
 & ' minimum    mer    banquise  villes   desert ',&
 & /,5F9.3)  
113 FORMAT(' ALBEDO ET EMISSIVITE :',/,&
 & ' minimum  maximum    mer    glacier  banquise  desert ',&
 & 2(/,6F9.3))  
114 FORMAT(' PROFONDEUR, % ARGILE, % SABLE, RESIS. MIN. :',/,&
 & ' minimum  maximum   desert ',&
 & 4(/,3F9.3))  

IF (LHOOK) CALL DR_HOOK('VAL923',1,ZHOOK_HANDLE)
END SUBROUTINE VAL923
