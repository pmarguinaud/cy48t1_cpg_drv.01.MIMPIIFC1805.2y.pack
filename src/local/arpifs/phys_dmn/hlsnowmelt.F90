      SUBROUTINE HLSNOWMELT(YDPHY2,PSNOWAB, PRAINAB,              & ! INPUT AND OUTPUT
     & PT,PTW,PPF,PDZ,PDPF,PQ,PCOV2D,PVSNOW,               & ! INPUT
     & PNEWRAIN)                                             ! OUTPUT
!     -----------
!     Purpose:  (*)

!     Parameterization of snow melting used for the RK-condensations scheme.

!     -----------
!     Method:
     
!     Based on the method by Sundqvist. With a modification by using wet bulb temp. 
!     instead of temperature. 
!     Output is the new rain and snow amount. Note that they are also input.
!     PNREWRAIN is pure output and is the amount of snow that has melt.

!     -----------
!     Author: 
!     K-I Ivarsson, 09-2011                                                         

!     -------------------------------------------------------------------
!     INTERFACE  :  SUBROUTINE 'HLSNOWMELT' IS CALLED 
!     ------------  FROM SUBROUTINE 'HLRKCOND', BUT IT IS POSSIBLE TO USE 
!                   IN OTHER CONDENSATION SCHEMES AS WELL 
!     INPUT-OUTPUT  ARGUMENTS  (ARGUMENTS D'ENTREE/SORTIE)
!     ----------------------------------------------------- 
!     PSNOWAB   : SNOW INTENSITY (KG/MS^2)
!     PRAINAB   : RAIN INTENSITY (KG/MS^2)
       
!     INPUT  ARGUMENTS  (ARGUMENTS D'ENTREE)
!     ----------------------------------------------------- 
!     PT        : TEMPERATURE (K)
!     PTW       : WET BULB TEMPERATURE (K)
!     PPF       : PRESSURE  (PA)
!     PDZ       : LAYER THICKNESS (M)
!     PDPF      : LAYER THICKNESS (PA)
!     PQ        : SPECIFIC HUMIDITY (KG/KG)
!     PCOV2D    : MAX CLOUD COVER OVER GRIDSQUARE (FRACTION)

!     OUTPUT  ARGUMENTS  (ARGUMENTS D'SORTIE)
!     ----------------------------------------------------- 
!     PNEWRAIN  : SNOW MELTING TO RAIN AS INTENSITY (KG/MS^2)

!     WORK  VARIABLES  :
!     ----------------------------------------------------- 
!     ZTIW: WET BULB. TEMPERATURE (K)
!     ZFALLTIME : THE TIME FOR PRECIPITATION TO FALL TROUGH A MODEL LAYER (S)
!     ZNEWSNOW     : SNOW INTENSITY AFTER MELTING (KG/MS^2)

!     1. DECLARATIONS. 
!     ==================================================================

!     1.1 MODULES USED
USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY  : LHOOK,   DR_HOOK
USE YOMCST   , ONLY : RG        ,&
 & RLVTT    ,RLSTT    ,RTT      ,&
 & RCPD
USE YOMPHY2  , ONLY : TPHY2

!---------------------------------------------------------------------
!1.2 ARGUMENT DECLARATIONS

!  INPUT-OUPUT ARGUMENTS :
IMPLICIT NONE
TYPE(TPHY2)       ,INTENT(IN)        :: YDPHY2
REAL(KIND=JPRB) ,  INTENT(INOUT)     ::  PSNOWAB,PRAINAB
!  INPUT ARGUMENTS :
REAL(KIND=JPRB) ,  INTENT(IN)        ::  PT,PTW,PPF,PDZ,PDPF,PQ,&
& PCOV2D,PVSNOW
!  OUTPUT ARGUMENTS :
REAL(KIND=JPRB) ,  INTENT(OUT)       ::  PNEWRAIN
!  WORKING VARIBLES : 
REAL(KIND=JPRB) ::  ZFALLTIME,ZNEWSNOW,ZKMELT
!  EXTERNAL FUNCTION:
REAL(KIND=JPRB)    ::  ZHOOK_HANDLE
IF(LHOOK) CALL DR_HOOK('HLSNOWMELT',0,ZHOOK_HANDLE)
ASSOCIATE(TSPHY=>YDPHY2%TSPHY)
 ZKMELT=2.E-4_JPRB 
 PNEWRAIN=0._JPRB

      IF(PSNOWAB > 0._JPRB .AND.PTW > RTT.AND.PT > RTT)THEN  

         ZFALLTIME = MIN(PDZ/PVSNOW,TSPHY)

         ZNEWSNOW = SQRT(MAX(0._JPRB,PSNOWAB)) -&
     &        0.5_JPRB*ZKMELT*(PTW-RTT)*SQRT(PVSNOW)*&
     &        SQRT(MAX(0.1_JPRB,PCOV2D))*&
     &        ZFALLTIME

         ZNEWSNOW = (MAX(0._JPRB,ZNEWSNOW))**2
!     Avoid melting below freezing:
         ZNEWSNOW = MAX (ZNEWSNOW,&
     &        PSNOWAB -&
     &        PDPF / RG*(PT-RTT)*RCPD/(RLSTT-RLVTT)/ TSPHY)
 

!     Final result
         PNEWRAIN = PSNOWAB - ZNEWSNOW
         PRAINAB  = PRAINAB + PNEWRAIN
         PSNOWAB = ZNEWSNOW 

      ENDIF
 END ASSOCIATE
 IF(LHOOK) CALL DR_HOOK('HLSNOWMELT',1,ZHOOK_HANDLE)
      END SUBROUTINE HLSNOWMELT
