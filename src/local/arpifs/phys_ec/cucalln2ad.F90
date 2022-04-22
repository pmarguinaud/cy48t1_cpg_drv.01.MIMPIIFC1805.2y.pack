SUBROUTINE CUCALLN2AD &
 & (YDTHF, YDCST, YDML_PHY_SLIN,     YDML_PHY_EC,&
 & KIDIA,    KFDIA,    KLON,   KLEV,&
 & LDLAND, LDRAIN1D, &
 & PTSPHY,&
 & PTM15,    PQM15,    PUM15,    PVM15,&
 & PVERVEL5, PQHFL5,   PAHFS5,   PAPHM15,&
 & PAP5,     PAPH5,    PGEO5,    PGEOH5, PGAW,&
 & PTENT5,   PGTENT5,  PTENQ5,   PGTENQ5,&
 & PTENU5,   PGTENU5,  PTENV5,   PGTENV5,&
 & PTENT25,  PTENQ25,  PTENU25,  PTENV25,&
 & KTYPE,&
 & KCBOT,    KCTOP,    KBOTSC,   LDCUM,   LDSC,&
 & PLU5,     PLUDE5,   PMFU5,    PMFD5,&
 & PFPLCL5,  PFPLCN5,&
 & KTRAC,    PCM15,    PTENC5,   PSCAV5,&
 & PTM1,     PQM1,     PUM1,     PVM1,&
 & PVERVEL,  PQHFL,    PAHFS,    PAPHM1,&
 & PAP,      PAPH,     PGEO,     PGEOH,&
 & PTENT ,   PGTENT,   PTENQ,    PGTENQ,&
 & PTENU,    PGTENU,   PTENV,    PGTENV,&
 & PLU,      PLUDE,    PMFU,     PMFD,&
 & PDIFCQ,   PDIFCS,   PFHPCL,   PFHPCN,&
 & PFPLCL,   PFPLCN,   PSTRCU,   PSTRCV,  PFCQLF ,PFCQIF,&
 & PCAPE,&
 & PCM1,     PTENC )  

!          *CUCALLN2AD* - MASTER ROUTINE - PROVIDES INTERFACE FOR:
!                        *CUMASTRN2AD* (CUMULUS PARAMETERIZATION)
!                        *CUCCDIAAD* (CUMULUS CLOUDS FOR RADIATION)
!                        *CUSTRATAD* (PBL_STRATOCUMULUS)
!                        (Adjoint)

!           ***  SIMPLIFIED CONVECTION SCHEME ***

!           P. LOPEZ     E.C.M.W.F.     15/12/2003

!           Adapted from 

!           M.TIEDTKE      E.C.M.W.F.     12/1989

!**   PURPOSE.
!     --------

!          *CUCALLN2AD* - INTERFACE FOR *CUMASTRN2AD*,*CUCCDIAAD* AND *CUSTRATAD*:
!                         PROVIDES INPUT FOR CUMASTRN2AD, CUCCDIANAD AND CUSTRATAD.
!                         RECEIVES UPDATED TENDENCIES, PRECIPITATION
!                         AND CONVECTIVE CLOUD PARAMETERS FOR RADIATION.

!**   INTERFACE.
!     ----------

!          *CUCALLN2AD* IS CALLED FROM *CALLPARAD*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS
!    *KTRAC*        NUMBER OF TRACERS

!     INPUT PARAMETERS (LOGICAL)

!    *LDLAND*       LAND SEA MASK (.TRUE. FOR LAND)

!     INPUT PARAMETERS (REAL)

!    *PTSPHY*       TIME STEP FOR THE PHYSICS                     S
!    *PTM15*        TEMPERATURE (T-1)                           (Trajectory)
!    *PQM15*        SPECIFIC HUMIDITY (T-1)                     (Trajectory)
!    *PUM15*        X-VELOCITY COMPONENT (T-1)                  (Trajectory)
!    *PVM15*        Y-VELOCITY COMPONENT (T-1)                  (Trajectory)
!    *PCM15*        CHEMICAL TRACERS                            (Trajectory)
!    *PVERVEL5*     VERTICAL VELOCITY                           (Trajectory)
!    *PQHFL5*       MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)      (Trajectory)
!    *PAHFS5*       SENSIBLE HEAT FLUX                          (Trajectory)
!    *PAPHM15*      PRESSURE ON HALF LEVELS                     (Trajectory)
!    *PAP5*         PROVISIONAL PRESSURE ON FULL LEVELS         (Trajectory)
!    *PAPH5*        PROVISIONAL PRESSURE ON HALF LEVELS         (Trajectory)
!    *PGEO5*        GEOPOTENTIAL                                (Trajectory)
!    *PGEOH5*       GEOPOTENTIAL ON HALF LEVELS                 (Trajectory)
!    *PSCAV5*       SCAVENGING COEFFICIENT                      (Trajectory)

!    *PTM1*         TEMPERATURE (T-1)                             K
!    *PQM1*         SPECIFIC HUMIDITY (T-1)                       KG/KG
!    *PUM1*         X-VELOCITY COMPONENT (T-1)                    M/S
!    *PVM1*         Y-VELOCITY COMPONENT (T-1)                    M/S
!    *PCM1*         CHEMICAL TRACERS                            
!    *PVERVEL*      VERTICAL VELOCITY                             PA/S
!    *PQHFL*        MOISTURE FLUX (EXCEPT FROM SNOW EVAP.)     KG/(SM2)
!    *PAHFS*        SENSIBLE HEAT FLUX                            W/M2
!    *PAPHM1*       PRESSURE ON HALF LEVELS                       PA
!    *PAP*          PROVISIONAL PRESSURE ON FULL LEVELS           PA 
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS           PA
!    *PGEO*         GEOPOTENTIAL                                  M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2

!    UPDATED PARAMETERS (REAL):

!    *PTENT5*       TEMPERATURE TENDENCY                       (Trajectory)
!    *PTENQ5*       MOISTURE TENDENCY                          (Trajectory)
!    *PTENU5*       TENDENCY OF U-COMP. OF WIND                (Trajectory)
!    *PTENV5*       TENDENCY OF V-COMP. OF WIND                (Trajectory)
!    *PTENC5*       TRACER TENDENCIES                          (Trajectory)

!    *PTENT*        TEMPERATURE TENDENCY                         K/S 
!    *PTENQ*        MOISTURE TENDENCY                           KG/(KG S)
!    *PTENU*        TENDENCY OF U-COMP. OF WIND                  M/S2
!    *PTENV*        TENDENCY OF V-COMP. OF WIND                  M/S2 
!    *PTENC*        TRACER TENDENCIES                          

!    OUTPUT PARAMETERS (INTEGER):

!    *KTYPE*        TYPE OF CONVECTION
!                       1 = PENETRATIVE CONVECTION
!                       2 = SHALLOW CONVECTION
!                       3 = MIDLEVEL CONVECTION 
!    *KCBOT*        CLOUD BASE LEVEL
!    *KCTOP*        CLOUD TOP LEVEL
!    *KBOTSC*       CLOUD BASE LEVEL FOR SC-CLOUDS

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDSC*         FLAG: .TRUE. FOR SC-POINTS

!    OUTPUT PARAMETERS (REAL):

!    *PLU5*         LIQUID WATER CONTENT IN UPDRAFTS          (Trajectory)
!    *PLUDE5*       DETRAINED LIQUID WATER                    (Trajectory) 
!    *PMFU5*        MASSFLUX UPDRAFTS                         (Trajectory) 
!    *PMFD5*        MASSFLUX DOWNDRAFTS                       (Trajectory) 
!    *PFPLCL5*      CONVECTIVE RAIN FLUX                      (Trajectory) 
!    *PFPLCN5*      CONVECTIVE SNOW FLUX                      (Trajectory) 
!    *PSTRTU5*      CONVECTIVE FLUX OF U-MOMEMTUM             (Trajectory)
!    *PSTRTV5*      CONVECTIVE FLUX OF V-MOMEMTUM             (Trajectory)
!    *PFCQLF5*      CONVECTIVE CONDENSATION FLUX LIQUID       (Trajectory) 
!    *PFCQIF5*      CONVECTIVE CONDENSATION FLUX ICE          (Trajectory) 

!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS            KG/KG
!    *PLUDE*        DETRAINED LIQUID WATER                      KG/(M3*S)
!    *PMFU*         MASSFLUX UPDRAFTS                           KG/(M2*S)
!    *PMFD*         MASSFLUX DOWNDRAFTS                         KG/(M2*S)
!    *PDIFCQ*       CONVECTIVE FLUX OF SPECIFIC HUMIDITY        KG/(M2*S)
!    *PDIFCS*       CONVECTIVE FLUX OF HEAT                      J/(M2*S)
!    *PFHPCL*       ENTHALPY FLUX DUE TO CONVECTIVE RAIN         J/(M2*S)
!    *PFHPCN*       ENTHALPY FLUX DUE TO CONVECTIVE SNOW         J/(M2*S)
!    *PFPLCL*       CONVECTIVE RAIN FLUX                        KG/(M2*S)
!    *PFPLCN*       CONVECTIVE SNOW FLUX                        KG/(M2*S)
!    *PSTRTU*       CONVECTIVE FLUX OF U-MOMEMTUM         (M/S)*KG/(M2*S)
!    *PSTRTV*       CONVECTIVE FLUX OF V-MOMEMTUM         (M/S)*KG/(M2*S)
!    *PFCQLF*       CONVECTIVE CONDENSATION FLUX LIQUID         KG/(M2*S)
!    *PFCQIF*       CONVECTIVE CONDENSATION FLUX ICE            KG/(M2*S)
!    *PCAPE*        CONVECTVE AVAILABLE POTENTIAL ENERGY         J/KG

!     EXTERNALS.
!     ----------

!          CUMASTRN2
!          CUCCDIA
!          CUSTRAT
!          SATUR
!          SATUR_1D
!          CUMASTRN2AD
!          CUCCDIAAD
!          CUSTRATAD
!          SATURAD

!          MODIFICATIONS
!          -------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        P. Lopez      14-02-2006  Added transport of tracers by mass flux
!        P. Lopez      11-01-2007  Removed switch LPHYLIN and added
!                                  routine SATUR_1D for 1D-Var rain
!        A. Geer       08-15-2008  Allow adjoint sensitivity to rain/snow flux
!                                  and LDRAIN1D name change to reflect usage
!        P. Lopez      17-09-2015  CAPE passed as argument.
!        P. Lopez      10-05-2016  Changed tendency computations to match full scheme

!----------------------------------------------------------------------

USE MODEL_PHYSICS_ECMWF_MOD      , ONLY : MODEL_PHYSICS_ECMWF_TYPE
USE MODEL_PHYSICS_SIMPLINEAR_MOD , ONLY : MODEL_PHYSICS_SIMPLINEAR_TYPE
USE PARKIND1                     , ONLY : JPIM     ,JPRB
USE YOMHOOK                      , ONLY : LHOOK,   DR_HOOK

USE YOMCST                       , ONLY : TCST  
USE YOETHF                       , ONLY : TTHF  
!USE YOMCT3                      , ONLY : NSTEP

IMPLICIT NONE

TYPE(TTHF)        ,INTENT(IN)    :: YDTHF
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(MODEL_PHYSICS_ECMWF_TYPE),INTENT(INOUT):: YDML_PHY_EC
TYPE(MODEL_PHYSICS_SIMPLINEAR_TYPE),INTENT(INOUT):: YDML_PHY_SLIN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
LOGICAL           ,INTENT(IN)    :: LDLAND(KLON)
LOGICAL           ,INTENT(IN)    :: LDRAIN1D 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM15(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM15(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUM15(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVM15(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCM15(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSCAV5(KTRAC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQHFL5(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAHFS5(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM15(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH5(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEO5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH5(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAW(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENT5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENQ5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENU5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENV5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGTENT5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGTENQ5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGTENU5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGTENV5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC5(KLON,KLEV,KTRAC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENT25(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENQ25(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENU25(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENV25(KLON,KLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KTYPE(KLON)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCBOT(KLON)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCTOP(KLON)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KBOTSC(KLON)
LOGICAL           ,INTENT(OUT)   :: LDCUM(KLON)
LOGICAL           ,INTENT(OUT)   :: LDSC(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLU5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLUDE5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFU5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFD5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLCL5(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLCN5(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCM1(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVERVEL(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQHFL(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFS(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAPHM1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAPH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGEO(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGEOH(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGTENT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGTENQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGTENU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGTENV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLUDE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMFU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PMFD(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFCQ(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFCS(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPCL(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPCN(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLCL(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLCN(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRCU(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRCV(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQLF(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCQIF(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCAPE(KLON) 

REAL(KIND=JPRB) ::     ZTP15(KLON,KLEV),       ZQP15(KLON,KLEV),&
 & ZUP15(KLON,KLEV),       ZVP15(KLON,KLEV),&
 & ZTU5(KLON,KLEV),        ZQU5(KLON,KLEV),&
 & ZQSAT5(KLON,KLEV),      ZRAIN5(KLON)   

REAL(KIND=JPRB) ::     ZTP1(KLON,KLEV),        ZQP1(KLON,KLEV),&
 & ZUP1(KLON,KLEV),        ZVP1(KLON,KLEV),&
 & ZTU(KLON,KLEV),         ZQU(KLON,KLEV),&
 & ZQSAT(KLON,KLEV),       ZRAIN(KLON)         
!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZENTHD5(KLON,KLEV),ZENTHS5(KLON,KLEV)  
REAL(KIND=JPRB) :: ZENTHD(KLON,KLEV),ZENTHS(KLON,KLEV)  
REAL(KIND=JPRB) :: ZCONDFLL(KLON),ZCONDFLN(KLON)

REAL(KIND=JPRB) :: ZTENQ55A(KLON,KLEV),         ZTENT55A(KLON,KLEV),&
 & ZTENU55A(KLON,KLEV),         ZTENV55A(KLON,KLEV),&
 & ZQSAT55A(KLON,KLEV)          

REAL(KIND=JPRB) :: ZTENQ55B(KLON,KLEV),         ZTENT55B(KLON,KLEV)

! Tracers
REAL(KIND=JPRB) :: ZCP15(KLON,KLEV,KTRAC), ZCP1(KLON,KLEV,KTRAC),&
                  &ZTENC55A(KLON,KLEV,KTRAC)

INTEGER(KIND=JPIM) :: IFLAG, JK, JL, JN

REAL(KIND=JPRB) :: ZALFAW5(KLON,KLEV), ZALFAW15(KLON,KLEV)
REAL(KIND=JPRB) :: ZALFAW
REAL(KIND=JPRB) :: ZCP5, ZGDPH5, Z1GDPH5
REAL(KIND=JPRB) :: ZCP , ZGDPH
LOGICAL :: LLTEST
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "cumastrn2ad.intfb.h"
#include "satur.intfb.h"
#include "satur_1d.intfb.h"
#include "saturad.intfb.h"

!DIR$ VFUNCTION EXPHF

!-----------------------------------------------------------------------

!              --------  TRAJECTORY COMPUTATIONS  --------

!-----------------------------------------------------------------------

!*    0.1          STORE T,Q,X,Y TENDENCIES FOR FLUX COMPUTATIONS
!*                 ----------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUCALLN2AD',0,ZHOOK_HANDLE)
ASSOCIATE(NJKT22=>YDML_PHY_SLIN%YRECUMF2%NJKT22, &
 & RLPTRC=>YDML_PHY_SLIN%YREPHLI%RLPTRC, &
 & LMFTRAC=>YDML_PHY_EC%YREPHY%LMFTRAC, &
 & LENCLD2=>YDML_PHY_SLIN%YRPHNC%LENCLD2)
! Setup of tendencies
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PTENQ5(JL,JK)=PGTENQ5(JL,JK)
    PTENT5(JL,JK)=PGTENT5(JL,JK)
    PTENV5(JL,JK)=PGTENV5(JL,JK)
    PTENU5(JL,JK)=PGTENU5(JL,JK)
    ZENTHD5(JL,JK)=0.0_JPRB
    ZENTHS5(JL,JK)=0.0_JPRB
  ENDDO
ENDDO
DO JK=1,KLEV+1
  DO JL=KIDIA,KFDIA
    PFPLCL5(JL,JK)=0.0_JPRB
    PFPLCN5(JL,JK)=0.0_JPRB
  ENDDO
ENDDO

!-----------------------------------------------------------------------

!*    1.           CALCULATE T,Q AND QS AT MAIN LEVELS
!*                 -----------------------------------

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZTP15(JL,JK)=PTM15(JL,JK)+PTENT5(JL,JK)*PTSPHY
    ZQP15(JL,JK)=PQM15(JL,JK)+PTENQ5(JL,JK)*PTSPHY
    ZUP15(JL,JK)=PUM15(JL,JK)+PTENU5(JL,JK)*PTSPHY
    ZVP15(JL,JK)=PVM15(JL,JK)+PTENV5(JL,JK)*PTSPHY
    ZQSAT5(JL,JK)=ZQP15(JL,JK)
  ENDDO
ENDDO
IF(KTRAC>0 .AND. LMFTRAC) THEN
  DO JN=1,KTRAC
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        ZCP15(JL,JK,JN)=PCM15(JL,JK,JN)+PTENC5(JL,JK,JN)*PTSPHY
      ENDDO
    ENDDO
  ENDDO
ENDIF

IFLAG=1
IF (LDRAIN1D) THEN
  CALL SATUR_1D (KIDIA , KFDIA , KLON  , NJKT22 , KLEV,&
   & PAP5  , ZTP15 , ZQSAT5 ) 
ELSE   
  CALL SATUR (YDTHF, YDCST, KIDIA , KFDIA , KLON  , NJKT22 , KLEV, YDML_PHY_SLIN%YREPHLI%LPHYLIN, &
   & PAP5  , ZTP15 , ZQSAT5, IFLAG  ) 
ENDIF   

DO JL=KIDIA,KFDIA
  ZRAIN5(JL)=0.0_JPRB
  LDCUM(JL)=.FALSE.
ENDDO

!-----------------------------------------------------------------------

!*    2.     CALL 'CUMASTRN2'(MASTER-ROUTINE FOR CUMULUS PARAMETERIZATION) 
!*           ------------------------------------------------------------ 

!        **   Storage of trajectory tendencies  **

ZTENT55A(:,:)=PTENT5(:,:)
ZTENQ55A(:,:)=PTENQ5(:,:)
ZTENU55A(:,:)=PTENU5(:,:)
ZTENV55A(:,:)=PTENV5(:,:)
IF(KTRAC>0 .AND. LMFTRAC) THEN
  ZTENC55A(:,:,:)=PTENC5(:,:,:)
ENDIF

!   Save specific humidity at saturation 
!   (possibly recomputed in CUFLXN)

ZQSAT55A(:,:)=ZQSAT5(:,:)

!      CALL CUMASTRN2
!     *  (  KIDIA,    KFDIA,    KLON,    KLEV 
!     *  ,  LDLAND,   LDRAIN1D, PTSPHY 
!     *  ,  ZTP15,    ZQP15,    ZUP15,    ZVP15 
!     *  ,  PVERVEL5, ZQSAT5,   PQHFL5,   PAHFS5 
!     *  ,  PAP5,     PAPH5,    PGEO5,    PGEOH5 
!     *  ,  PTENT5,   PTENQ5,   PTENU5,   PTENV5 
!     *  ,  LDCUM,    KTYPE,    KCBOT,    KCTOP
!     *  ,  KBOTSC,   LDSC
!     *  ,  ZTU5,     ZQU5,     PLU5,     PLUDE5
!     *  ,  ZENTHD5,  PFPLCL5,  PFPLCN5,  ZRAIN5
!     *  ,  PMFU5,    PMFD5
!     *  ,  KTRAC,    ZCP15,    PTENC5 )

!----------------------------------------------------------------------

! SECTIONS 3 AND 4 ARE ONLY REQUIRED IF THE DIAGNOSTIC CLOUD SCHEME
! IS USED  : LEPCLD=.FALSE.

!----------------------------------------------------------------------

LLTEST = .NOT.LENCLD2.AND..NOT.LDRAIN1D

IF (LLTEST) THEN
!-----------------------------------------------------------------------

!*    3.0       CALL 'CUCCDIA' TO UPDATE CLOUD PARAMETERS FOR RADIATION
!               -------------------------------------------------------

!      CALL CUCCDIA
!     *  (  KIDIA,    KFDIA,    KLON,   KLEV
!     *  ,  NSTEP,    KCBOT,    KCTOP
!     *  ,  LDCUM,    ZQU5,     PLU5,     PMFU5,   ZRAIN5
!     *  ,  PARPRC5,  KTOPC,    KBASEC  )

!---------------------------------------------------------------------

!*    4.0          CALL 'CUSTRAT' FOR PARAMETERIZATION OF PBL-CLOUDS
!                  -------------------------------------------------

!        **   Storage of trajectory tendencies  **

  ZTENT55B(:,:)=PTENT5(:,:)
  ZTENQ55B(:,:)=PTENQ5(:,:)

!         CALL CUSTRAT
!     *     (  KIDIA,    KFDIA,    KLON,   KLEV
!     *     ,  LDCUM,    PTSPHY
!     *     ,  PAP5,     PAPH5,    PGEO5
!     *     ,  ZTP15,    ZQP15,    ZQSAT5,  ZENTHS5
!     *     ,  PTENT5,   PTENQ5)

ENDIF

!---------------------------------------------------------------------

! Converting outputs to process tendencies
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PTENQ5(JL,JK)=PTENQ5(JL,JK)-PGTENQ5(JL,JK)
    PTENT5(JL,JK)=PTENT5(JL,JK)-PGTENT5(JL,JK)
    PTENV5(JL,JK)=PTENV5(JL,JK)-PGTENV5(JL,JK)
    PTENU5(JL,JK)=PTENU5(JL,JK)-PGTENU5(JL,JK)
  ENDDO
ENDDO

!-----------------------------------------------------------------------

!*    5.           FLUX COMPUTATIONS
!                  -----------------

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZALFAW5(JL,JK)=0.545_JPRB*(TANH(0.17_JPRB*(ZTP15(JL,JK)-RLPTRC))+1.0_JPRB)
    ZALFAW15(JL,JK)=ZALFAW5(JL,JK)
    IF (ZALFAW15(JL,JK) > 1.0_JPRB) THEN
      ZALFAW5(JL,JK) = 1.0_JPRB
    ENDIF
  ENDDO
ENDDO

!------------------------------------------------------------------------

!                --------  ADJOINT COMPUTATIONS  --------

!------------------------------------------------------------------------

!*                           INITIALISATIONS
!                            ---------------

ZTP1(:,:)  = 0.0_JPRB
ZQP1(:,:)  = 0.0_JPRB
ZUP1(:,:)  = 0.0_JPRB
ZVP1(:,:)  = 0.0_JPRB
ZTU(:,:)   = 0.0_JPRB
ZQU(:,:)   = 0.0_JPRB
ZQSAT(:,:) = 0.0_JPRB
ZRAIN(:)   = 0.0_JPRB
ZENTHD(:,:)= 0.0_JPRB
ZENTHS(:,:)= 0.0_JPRB
ZCONDFLL(:)= 0.0_JPRB
ZCONDFLN(:)= 0.0_JPRB
IF(KTRAC>0 .AND. LMFTRAC) THEN
  ZCP1(:,:,:)  = 0.0_JPRB
ENDIF

PDIFCQ(:,:)=0.0_JPRB
PDIFCS(:,:)=0.0_JPRB
PSTRCU(:,:)=0.0_JPRB
PSTRCV(:,:)=0.0_JPRB
PFCQLF(:,:)=0.0_JPRB
PFCQIF(:,:)=0.0_JPRB
PFHPCL(:,:)=0.0_JPRB
PFHPCN(:,:)=0.0_JPRB

DO JK=KLEV,1,-1
!DEC$ IVDEP
  DO JL=KIDIA,KFDIA

    ZALFAW=0.0_JPRB

!... enthalpy flux due to precipitations
    
    ZCONDFLN(JL)=ZCONDFLN(JL)+PFCQIF(JL,JK+1)
    PFCQIF(JL,JK+1)=0.0_JPRB

    ZCONDFLL(JL)=ZCONDFLL(JL)+PFCQLF(JL,JK+1)
    PFCQLF(JL,JK+1)=0.0_JPRB

    PFPLCL(JL,JK+1)=PFPLCL(JL,JK+1)-PDIFCQ(JL,JK+1)
    PFPLCN(JL,JK+1)=PFPLCN(JL,JK+1)-PDIFCQ(JL,JK+1)
    ZCONDFLL(JL)=ZCONDFLL(JL)-PDIFCQ(JL,JK+1)
    ZCONDFLN(JL)=ZCONDFLN(JL)-PDIFCQ(JL,JK+1)

    PFHPCL(JL,JK+1)=PFHPCL(JL,JK+1)-PDIFCS(JL,JK+1)
    PFHPCN(JL,JK+1)=PFHPCN(JL,JK+1)-PDIFCS(JL,JK+1)
    ZCONDFLL(JL)=ZCONDFLL(JL)+YDCST%RLVTT*PDIFCS(JL,JK+1)
    ZCONDFLN(JL)=ZCONDFLN(JL)+YDCST%RLSTT*PDIFCS(JL,JK+1)

    PLUDE(JL,JK)=PLUDE(JL,JK)+(1.0_JPRB-ZALFAW5(JL,JK))*ZCONDFLN(JL)
    ZALFAW=ZALFAW-PLUDE5(JL,JK)*ZCONDFLN(JL)

    PLUDE(JL,JK)=PLUDE(JL,JK)+ZALFAW5(JL,JK)*ZCONDFLL(JL)
    ZALFAW=ZALFAW+PLUDE5(JL,JK)*ZCONDFLL(JL)

    IF (ZALFAW15(JL,JK) > 1.0_JPRB) THEN
      ZALFAW=0.0_JPRB
    ENDIF

    ZTP1(JL,JK)=ZTP1(JL,JK)+ ZALFAW*0.17_JPRB*0.545_JPRB &
     & *(1.0_JPRB-TANH(0.17_JPRB*(ZTP15(JL,JK)-RLPTRC))**2)  
    ZALFAW=0.0_JPRB

    PFPLCN(JL,JK+1)=PFPLCN(JL,JK+1)-PFHPCN(JL,JK+1)*YDCST%RLSTT
    PFHPCN(JL,JK+1)=0.0_JPRB

    PFPLCL(JL,JK+1)=PFPLCL(JL,JK+1)-PFHPCL(JL,JK+1)*YDCST%RLVTT
    PFHPCL(JL,JK+1)=0.0_JPRB

  ENDDO
ENDDO

DO JK=KLEV,1,-1
  DO JL=KIDIA,KFDIA

    ZGDPH5  = -YDCST%RG/( PAPHM15(JL,JK+1)-PAPHM15(JL,JK) )
    Z1GDPH5 = 1.0_JPRB/ZGDPH5
    ZCP5 = YDCST%RCPD*(1+YDTHF%RVTMP2*PQM15(JL,JK))

    ZGDPH=0.0_JPRB
    ZCP=0.0_JPRB

!...increments in U,V
    PSTRCV(JL,JK)=PSTRCV(JL,JK)+PSTRCV(JL,JK+1)
    PTENV(JL,JK)=PTENV(JL,JK)+Z1GDPH5*PSTRCV(JL,JK+1)
    ZGDPH=ZGDPH-PTENV25(JL,JK)*Z1GDPH5**2 *PSTRCV(JL,JK+1)  
    PSTRCV(JL,JK+1)=0.0_JPRB

    PSTRCU(JL,JK)=PSTRCU(JL,JK)+PSTRCU(JL,JK+1)
    PTENU(JL,JK)=PTENU(JL,JK)+Z1GDPH5*PSTRCU(JL,JK+1)
    ZGDPH=ZGDPH-PTENU25(JL,JK)*Z1GDPH5**2 *PSTRCU(JL,JK+1)  
    PSTRCU(JL,JK+1)=0.0_JPRB

!...increment in Q
    PDIFCQ(JL,JK)=PDIFCQ(JL,JK)+PDIFCQ(JL,JK+1)
    PTENQ(JL,JK)=PTENQ(JL,JK)+Z1GDPH5*PDIFCQ(JL,JK+1)
    ZGDPH=ZGDPH-PTENQ25(JL,JK)*Z1GDPH5**2 *PDIFCQ(JL,JK+1)  
    PDIFCQ(JL,JK+1)=0.0_JPRB
    
!...increment in dry static energy is converted to flux of d.s.e.
    PDIFCS(JL,JK)=PDIFCS(JL,JK)+PDIFCS(JL,JK+1)
    PTENT(JL,JK)=PTENT(JL,JK)+ZCP5*Z1GDPH5*PDIFCS(JL,JK+1)
    ZCP=ZCP+PTENT25(JL,JK)*Z1GDPH5*PDIFCS(JL,JK+1)
    ZGDPH=ZGDPH-PTENT25(JL,JK)*(Z1GDPH5**2)*ZCP5*PDIFCS(JL,JK+1)  
    PDIFCS(JL,JK+1)=0.0_JPRB

    PQM1(JL,JK)=PQM1(JL,JK)+YDCST%RCPD*YDTHF%RVTMP2*ZCP
    ZCP=0.0_JPRB

    PAPHM1(JL,JK+1)=PAPHM1(JL,JK+1)+ZGDPH*(ZGDPH5**2)/YDCST%RG
    PAPHM1(JL,JK)=PAPHM1(JL,JK)-ZGDPH*(ZGDPH5**2)/YDCST%RG
    ZGDPH=0.0_JPRB

  ENDDO
ENDDO

DO JL=KIDIA,KFDIA
  PDIFCQ (JL,1)=0.0_JPRB
  PDIFCS (JL,1)=0.0_JPRB
  PSTRCU (JL,1)=0.0_JPRB
  PSTRCV (JL,1)=0.0_JPRB
  PFCQLF (JL,1)=0.0_JPRB
  PFCQIF (JL,1)=0.0_JPRB
  PFHPCL (JL,1)=0.0_JPRB
  PFHPCN (JL,1)=0.0_JPRB
  ZCONDFLL (JL)=0.0_JPRB
  ZCONDFLN (JL)=0.0_JPRB
ENDDO

!-----------------------------------------------------------------------

! Converting outputs to process tendencies
DO JK=KLEV,1,-1
  DO JL=KIDIA,KFDIA
    PGTENQ(JL,JK) =PGTENQ(JL,JK)-PTENQ(JL,JK)
    PGTENT(JL,JK) =PGTENT(JL,JK)-PTENT(JL,JK)
    PGTENV(JL,JK) =PGTENV(JL,JK)-PTENV(JL,JK)
    PGTENU(JL,JK) =PGTENU(JL,JK)-PTENU(JL,JK)
  ENDDO
ENDDO

!-----------------------------------------------------------------------

!*    4.0          CALL 'CUSTRATAD' FOR PARAMETERIZATION OF PBL-CLOUDS
!                  ---------------------------------------------------

IF (LLTEST) THEN

!         CALL CUSTRATAD
!     *     (  KIDIA,    KFDIA,    KLON,   KLEV
!     *     ,  LDCUM,    PTSPHY
!     *     ,  PAP5,     PAPH5,    PGEO5
!     *     ,  ZTP15,    ZQP15,    ZQSAT5,  ZENTHS5
!     *     ,  ZTENT55B, ZTENQ55B
!     *     ,  PAP,      PAPH,     PGEO
!     *     ,  ZTP1,     ZQP1,     ZQSAT,   ZENTHS
!     *     ,  PTENT,    PTENQ)

ENDIF

!-----------------------------------------------------------------------

!*    2.     CALL 'CUMASTRN2AD'(MASTER-ROUTINE FOR SIMPLIFIED
!*                              CUMULUS PARAMETERIZATION) 
!*           ------------------------------------------------ 

CALL CUMASTRN2AD &
 & (YDTHF, YDCST, YDML_PHY_SLIN,     YDML_PHY_EC,&
 & KIDIA,    KFDIA,    KLON,    KLEV,&
 & LDLAND,   LDRAIN1D, PTSPHY,&
 & ZTP15,    ZQP15,    ZUP15,    ZVP15,&
 & PVERVEL5, ZQSAT55A, PQHFL5,   PAHFS5,&
 & PAP5,     PAPH5,    PGEO5,    PGEOH5, PGAW,&
 & ZTENT55A, ZTENQ55A, ZTENU55A, ZTENV55A,&
 & LDCUM,    KTYPE,    KCBOT,    KCTOP,&
 & KBOTSC,   LDSC,&
 & ZTU5,     ZQU5,     PLU5,     PLUDE5,&
 & ZENTHD5,  PFPLCL5,  PFPLCN5,  ZRAIN5,&
 & PMFU5,    PMFD5,&
 & KTRAC,    ZCP15,    ZTENC55A, PSCAV5,&
 & ZTP1,     ZQP1,     ZUP1,     ZVP1,&
 & PVERVEL,  ZQSAT,    PQHFL,    PAHFS,&
 & PAP,      PAPH,     PGEO,     PGEOH,&
 & PTENT,    PTENQ,    PTENU,    PTENV,&
 & ZTU,      ZQU,      PLU,      PLUDE,&
 & ZENTHD,   PFPLCL,   PFPLCN,   ZRAIN,&
 & PMFU,     PMFD,     PCAPE,&
 & ZCP1,     PTENC )  

DO JL=KIDIA,KFDIA
  ZRAIN(JL)=0.0_JPRB
ENDDO

!-----------------------------------------------------------------------

!*    1.           CALCULATE T,Q AND QS AT MAIN LEVELS
!*                 -----------------------------------

CALL SATURAD (KIDIA , KFDIA , KLON  , NJKT22 , KLEV,&
 & PAP5  , ZTP15 , ZQSAT5, &
 & PAP   , ZTP1  , ZQSAT   )  

IF(KTRAC>0 .AND. LMFTRAC) THEN
  DO JN=1,KTRAC
    DO JK=1,KLEV
      DO JL=KIDIA,KFDIA
        PTENC (JL,JK,JN) = PTENC(JL,JK,JN) + ZCP1(JL,JK,JN)*PTSPHY
        PCM1 (JL,JK,JN) = PCM1(JL,JK,JN) + ZCP1(JL,JK,JN)
        ZCP1 (JL,JK,JN) = 0.0_JPRB
      ENDDO
    ENDDO
  ENDDO
ENDIF
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZQP1 (JL,JK) = ZQP1 (JL,JK) + ZQSAT (JL,JK)
    ZQSAT (JL,JK) = 0.0_JPRB
    PTENV (JL,JK) = PTENV (JL,JK) + ZVP1 (JL,JK)*PTSPHY
    PVM1 (JL,JK) = PVM1 (JL,JK) + ZVP1 (JL,JK)
    ZVP1 (JL,JK) = 0.0_JPRB
    PTENU (JL,JK) = PTENU (JL,JK) + ZUP1 (JL,JK)*PTSPHY
    PUM1 (JL,JK) = PUM1 (JL,JK) + ZUP1 (JL,JK)
    ZUP1 (JL,JK) = 0.0_JPRB
    PTENQ (JL,JK) = PTENQ (JL,JK) + ZQP1 (JL,JK)*PTSPHY
    PQM1 (JL,JK) = PQM1 (JL,JK) + ZQP1 (JL,JK)
    ZQP1 (JL,JK) = 0.0_JPRB
    PTENT (JL,JK) = PTENT (JL,JK) + ZTP1 (JL,JK)*PTSPHY
    PTM1 (JL,JK) = PTM1 (JL,JK) + ZTP1 (JL,JK)
    ZTP1 (JL,JK) = 0.0_JPRB
  ENDDO
ENDDO

!-----------------------------------------------------------------------

!*    0.1          STORE T,Q,X,Y TENDENCIES FOR FLUX COMPUTATIONS
!*                 ----------------------------------------------

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PGTENQ(JL,JK)=PGTENQ(JL,JK)+PTENQ(JL,JK)
    PGTENT(JL,JK)=PGTENT(JL,JK)+PTENT(JL,JK)
    PGTENU(JL,JK)=PGTENU(JL,JK)+PTENU(JL,JK)
    PGTENV(JL,JK)=PGTENV(JL,JK)+PTENV(JL,JK)
    PTENQ(JL,JK)=0.0_JPRB
    PTENT(JL,JK)=0.0_JPRB
    PTENU(JL,JK)=0.0_JPRB
    PTENV(JL,JK)=0.0_JPRB
  ENDDO
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUCALLN2AD',1,ZHOOK_HANDLE)
END SUBROUTINE CUCALLN2AD
