!OPTIONS XOPT(NOEVAL)
SUBROUTINE HLEVAPPREC(YDPHY0,YDPHY2,K,PSNOWAB,PRAINAB,                & !INPUT
  & PPF,PT,PQS,                                            & 
  & PDPF,PQ,PCOV2D,                                        & 
  & PSNOWSHAPE,PRAINSHAPE,                                 & 
!   OUTPUT
  & PSNOWRES,PRAINRES,PDZ,                                 & 
  & PEVAPS,PEVAPR,PCEVAP)                                       

!----------------------------------------------------------------------- 
!     Purpose:
!     ---------
!     Parameterization of evaporation of precipitation, used for the RK-
!     condensation scheme.

!     ---------
!     Externes.                                                         
!     ---------

!     Method.

!     Computations are partly based on the Sundqvist scheme, with some
!     extentions from the STRACO scheme. The pararmaterization include
!     depositional growth of snow in case of supersaturation with
!     respect to snow.  

!     --------
!     Author.

!     K-I. Ivarsson 09-2011

!     -------------------------------------------------------------------
!     INTERFACE  :  SUBROUTINE 'HLEVAPPREC' IS CALLED 
!     ------------  FROM SUBROUTINE 'HLRKCOND', OR 'PCOND' BUT IT IS POSSIBLE 
!                   TO USE IN OTHER CONDENSATION SCHEMES AS WELL
!                   
!     INPUT  ARGUMENTS  (ARGUMENTS D'ENTREE)
!     ----------------------------------------------------- 
!     K,        : MODEL LEVEL
!     PSNOWAB   : SNOW INTENSITY (KG/MS^2)
!     PRAINAB   : RAIN INTENSITY (KG/MS^2)
!     PPF       : PRESSURE  (PA)
!     PT        : TEMPERATURE (K)
!     PQS       : SATURATION SPEC HUM. OF WATER VAPOR (MIXED WATER -ICE)
!     PDPF      : LAYER THICKNESS (PA)
!     PQ        : SPECIFIC HUMIDITY (KG/KG)
!     PCOV2D    : MAX CLOUD COVER OVER GRIDSQUARE (FRACTION)
!     PSNOWSHAPE: DEGREE OF SPHERICALNESS OF SNOW 1= PERFECT SPHERES (~GRAUPEL)
!                 > 1 MORE LIKE NEEDELS, DENTRITES ETC.
!     PRAINSHAPE: DEGREE OF (NO) DRIZZLENESS : 0 MORE OF DRIZZLE TYPE. 1:
!                 NOT OF DRIZZLE TYPE
!     OUTPUT  ARGUMENTS  (ARGUMENTS D'SORTIE)
!     ----------------------------------------------------- 
!     PEVAPS,  : RATE OF EVAP OF FALLIN SNOW  (1/S) 
!     PEVAPR,  : RATE OF EVAP OF FALLING RAIN  (1/S) 
!     PCEVAP,  : RATE OF EVAP OF FALLING PRECIP  IN THE CLODY AREA (1/S)
!                SHOULD BE ZERO OR NEGATIVE. ( NEGATIVE FOR SNOW IF QSI > 1 )
!     PDZ      : THICKNESS OF MODEL LEVEL IN METRES  
!     PSNOWRES : SNOW INTENSITY (KG/MS^2) AFTER EVAPOATION, DEPOSITION
!     PRAINRES : RAIN INTENSITY (KG/MS^2) AFTER EVAPOATION


!     WORK  VARIABLES  :
!     ----------------------------------------------------- 
!     ZQSI      : SATURATION SPEC HUM. OF WATER VAPOR OVER ICE
!     ZQSW      : SATURATION SPEC HUM. OF WATER VAPOR OVER WATER
!     ZFALLTIME: THE TIME FOR PRECIPITATION TO FALL TROUGH A MODEL LAYER (S)
!     ZHKEVAP_MIX : HKEVAP MODYFIED DUE TO ASSUMED LAGRER SURFACE
!                  OF SNOWFLAKES THAN OF RAIN FOR THE SAME PRECIP AMOUNT. 


!     ==================================================================
!     1. DECLARATIONS. 
!     ==================================================================
!     1.1 MODULES USED   
!-----------------------------------------------------------------------

USE YOMHOOK   ,ONLY  : LHOOK,   DR_HOOK
USE YOMCST   , ONLY : RG        ,RDT      ,&
 & RV       ,RGAMD    ,RCPV     ,RETV     ,RCW      ,&
 & RCS      ,RLVTT    ,RLSTT    ,RTT      ,RALPW    ,&
 & RBETW    ,RGAMW    ,RALPS    ,RBETS    ,RGAMS    ,&
 & RALPD    ,RBETD    ,RMD      ,RMV
USE YOMPHY0  , ONLY : TPHY0
USE YOMPHY2  , ONLY : TPHY2
USE PARKIND1  ,ONLY : JPRB,JPIM

!
!---------------------------------------------------------------------
!  1.2 ARGUMENT DECLARATIONS
IMPLICIT NONE 
!  INPUT ARGUMENTS :
TYPE(TPHY0)         ,INTENT(IN) :: YDPHY0
TYPE(TPHY2)         ,INTENT(IN) :: YDPHY2
INTEGER(KIND=JPIM),  INTENT(IN) ::   K
REAL(KIND=JPRB) ,    INTENT(IN) ::   PSNOWAB,PRAINAB,PPF,&
& PT,PQS,PDPF,PQ,PCOV2D,PSNOWSHAPE,PRAINSHAPE  
!  OUTPUT ARGUMENTS :
REAL(KIND=JPRB) ,  INTENT(OUT)  ::   PEVAPS,PEVAPR,PCEVAP,&
& PSNOWRES,PRAINRES,PDZ

!  WORK VARIABLES :
REAL(KIND=JPRB)                 ::   ZQSI,ZQSW,&
& ZFALLTIME,ZHKEVAP_MIX,ZHVSNOW,ZHVRAIN,ZVDRIZ,&
& ZEPSILO,ZEW,ZEI,ZKEVAP,ZRHOE,ZRAINAB,ZSNOWAB
!  EXTERNAL FUNCTIONS:
REAL(KIND=JPRB)                 :: ZHOOK_HANDLE

#include "fctdoi.func.h"
#include "fcttrm.func.h"

IF (LHOOK) CALL DR_HOOK('HLEVAPPREC',0,ZHOOK_HANDLE)
ASSOCIATE(RDTFAC=>YDPHY0%RDTFAC, &
 & TSPHY=>YDPHY2%TSPHY)
 PCEVAP  = 0._JPRB
 PEVAPR = 0._JPRB
 PEVAPS  = 0._JPRB
 ZHVRAIN = 5._JPRB
 ZHVSNOW = 1._JPRB
 ZVDRIZ = 2._JPRB
 ZKEVAP = 1.2E-2_JPRB
 ZEPSILO=0.622_JPRB  !RMV/RMD
 PRAINRES = 0._JPRB
 PSNOWRES = 0._JPRB

      IF(PSNOWAB+PRAINAB > 0._JPRB)THEN
         ZRHOE=PPF/PT/287._JPRB
!        conversion from specific precipitation to kg/m S^2 (grid mean)
         ZSNOWAB = PSNOWAB*ZRHOE*ZHVSNOW
         ZRAINAB = PRAINAB*ZRHOE*ZHVRAIN
         PDZ = PDPF / RG / ZRHOE

         IF(ZRAINAB > 1.0E-20_JPRB)THEN
            ZEW= FOEW (PT,0._JPRB)
            ZQSW = MIN(ZEPSILO*ZEW/(PPF - ZEW),1._JPRB)
            ZHVRAIN = ZHVRAIN*PRAINSHAPE + ZVDRIZ*(1._JPRB-PRAINSHAPE)
            ZFALLTIME = MIN(PDZ/ZHVRAIN,TSPHY)
            PEVAPR = SQRT(ZRAINAB) -&
     &        0.5_JPRB*ZKEVAP*(ZQSW-PQ)*&
     &        SQRT(MAX(0.1_JPRB,PCOV2D)/ZHVRAIN)*&
     &        ZFALLTIME

            PEVAPR = MAX(0._JPRB,PEVAPR)**2
!           Residual rain in kg/m^2
            PRAINRES = MAX(0._JPRB,ZRAINAB-PEVAPR)

            PEVAPR =  (ZRAINAB - PEVAPR)*RG / PDPF

         ENDIF

         IF(ZSNOWAB > 1.0E-20_JPRB)THEN
            ZEI= FOEW (PT,1._JPRB)
            ZQSI = MIN(ZEPSILO*ZEI/(PPF - ZEI),1._JPRB)
            ZFALLTIME = MIN(PDZ/ZHVSNOW,TSPHY)
            ZHKEVAP_MIX = ZKEVAP*PSNOWSHAPE         
            PEVAPS = SQRT(ZSNOWAB) -&
     &        0.5_JPRB*ZHKEVAP_MIX*(ZQSI-PQ)*&
     &        SQRT(MAX(0.1_JPRB,PCOV2D)/ZHVSNOW)*&
     &        ZFALLTIME

            PEVAPS = MAX(0._JPRB,PEVAPS)**2

!           Residual snow in kg/m s^2

            PSNOWRES = MAX(0._JPRB,ZSNOWAB-PEVAPS)
 
            PEVAPS =  (ZSNOWAB - PEVAPS)*RG / PDPF 

!     Note that negative evaporation (depositional snow growth)
!     is allowed here, since it is an important part of the BF-process

!     For correct implementation in the CAM3 code, estimate snow
!     deposition in the cloudy area

            IF(PT < RTT)THEN                                             
               PCEVAP= SQRT(ZSNOWAB) - 0.5_JPRB*ZHKEVAP_MIX*(ZQSI-PQS)*&
     &            SQRT(MAX(0.1_JPRB,PCOV2D)/ZHVSNOW)*ZFALLTIME
               PCEVAP = MAX(0._JPRB,PCEVAP)**2
               PCEVAP =MIN(0._JPRB,(ZSNOWAB - PCEVAP)*RG/PDPF)
            ENDIF
         ENDIF
      ENDIF
  END ASSOCIATE
  IF(LHOOK) CALL DR_HOOK('HLEVAPPREC',1,ZHOOK_HANDLE)
END SUBROUTINE HLEVAPPREC
