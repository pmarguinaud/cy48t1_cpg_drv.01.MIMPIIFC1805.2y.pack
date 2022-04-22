SUBROUTINE CUPDRAAD &
 & ( YDCUMFS,YDECUMF2, YDECLDP,&
 & KIDIA,    KFDIA,    KLON,     KLEV,&
 & K_JKSTART,  LD_LLGO_ON,&
 & PQENH5,   PGEO5,    PGEOH5,   PAPH5,    PQSEN5,&
 & PTU5,     PQU5,     PLU5,&
 & PSENH5,   PWU2H5,   PSUH5,&
 & PBUOH5,   PLGLAC5,  PQPRCV5,  PCAPE5,&
 & KLAB,     LDCUM,    LDSC,     KCBOT,&
 & KBOTSC,   KCTOP,    KSTUP,&
 & PQENH,    PGEO,     PGEOH,    PAPH,     PQSEN,&
 & PTU,      PQU,      PLU,&
 & PSENH,    PWU2H,    PSUH,&
 & PBUOH,    PLGLAC,   PQPRCV,   PCAPE )  

!          THIS ROUTINE CALCULATES CHARACTERISTICS OF 
!          THE CUMULUS UPDRAFTS, CLOUD BASE HEIGHT AND 
!          CLOUD TOP HEIGHT (ADJOINT VERSION).

!          Ph LOPEZ (ECMWF) (07/2002) 
!          inspired from A. Pier Siebesma (KNMI)
!          and C. Jakob (ECMWF) (01/2001)

!          PURPOSE.
!          --------
!          TO PRODUCE UPDRAFT CALCULATION, CLOUD BASE AND CLOUD
!          TOP VALUES FOR SIMPLIFIED CONVECTIVE PARAMETRIZATION

!          INTERFACE
!          ---------
!          THIS ROUTINE IS CALLED FROM *CUBASEN2AD*.
!          INPUT ARE ENVIRONM. VALUES OF T,Q,P,PHI AT HALF LEVELS.
!          IT RETURNS CLOUD FIELDS VALUES AND FLAGS AS FOLLOWS;
!                 KLAB=0 FOR STABLE LAYERS
!                 KLAB=1 FOR SUBCLOUD LEVELS
!                 KLAB=2 FOR CLOUD LEVELS LEVEL

!          METHOD.
!          --------
!          LIFT SURFACE AIR DRY-ADIABATICALLY TO CLOUD TOP.
!          ENTRAINING PLUME WITH ENTRAINMENT PROPORTIONAL TO 1/Z
!          FOR SHALLOW CONVECTION (STARTING FROM LOWEST MODEL LEVEL)
!          AND FIXED OTHERWISE.

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!    *KIDIA*        START POINT
!    *KFDIA*        END POINT
!    *KLON*         NUMBER OF GRID POINTS PER PACKET
!    *KLEV*         NUMBER OF LEVELS

!    INPUT PARAMETERS (REAL):

!    *PQENH5*       ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      (Trajectory)
!    *PGEO5*        GEOPOTENTIAL ON FULL LEVELS                   (Trajectory)
!    *PGEOH5*       GEOPOTENTIAL ON HALF LEVELS                   (Trajectory)
!    *PAPH5*        PROVISIONAL PRESSURE ON HALF LEVELS           (Trajectory)
!    *PQSEN5*       PROVISIONAL ENVIRONMENT SATU. HUMIDITY (T+1)  (Trajectory)

!    *PQENH*        ENV. SPEC. HUMIDITY (T+1) ON HALF LEVELS      KG/KG
!    *PGEO*         GEOPOTENTIAL ON FULL LEVELS                   M2/S2
!    *PGEOH*        GEOPOTENTIAL ON HALF LEVELS                   M2/S2
!    *PAPH*         PROVISIONAL PRESSURE ON HALF LEVELS           PA
!    *PQSEN*        PROVISIONAL ENVIRONMENT SATU. HUMIDITY (T+1)  KG/KG

!    UPDATED PARAMETERS (REAL):

!    *PTU5*         TEMPERATURE IN UPDRAFTS                       (Trajectory)
!    *PQU5*         SPEC. HUMIDITY IN UPDRAFTS                    (Trajectory)
!    *PLU5*         LIQUID WATER CONTENT IN UPDRAFTS              (Trajectory)
!    *PBUOH5*       BUOYANCY ON HALF LEVELS                       (Trajectory)

!    *PTU*          TEMPERATURE IN UPDRAFTS                       K
!    *PQU*          SPEC. HUMIDITY IN UPDRAFTS                    KG/KG
!    *PLU*          LIQUID WATER CONTENT IN UPDRAFTS              KG/KG
!    *PBUOH*        BUOYANCY ON HALF LEVELS                       M/S2

!    UPDATED PARAMETERS (INTEGER):

!    *KLAB*         FLAG KLAB=1 FOR SUBCLOUD LEVELS
!                        KLAB=2 FOR CLOUD LEVELS

!    OUTPUT PARAMETERS (LOGICAL):

!    *LDCUM*        FLAG: .TRUE. FOR CONVECTIVE POINTS 
!    *LDSC*         FLAG: .TRUE. IF BL-CLOUDS EXIST

!    OUTPUT PARAMETERS (INTEGER):

!    *KSTUP*       LEVEL AT WHICH THE UPDRAFT STARTS 
!    *KCBOT*       CLOUD BASE LEVEL     
!    *KCTOP*       CLOUD TOP LEVEL = HEIGHEST HALF LEVEL 
!                  WITH A NON-ZERO CLOUD UPDRAFT.
!    *KBOTSC*      CLOUD BASE LEVEL OF BL-CLOUDS

!    OUTPUT PARAMETERS (REAL):

!    *PLGLAC5*     FROZEN CLOUD WATER CONTENT                   (Trajectory)
!    *PQPRCV5*     CONVECTIVE PRECIPITATION CONTENT             (Trajectory)
!    *PCAPE5*      CONVECTIVE AVAILABLE POTENTIAL ENERGY        (Trajectory)

!    *PLGLAC*      FROZEN CLOUD WATER CONTENT                   KG/KG 
!    *PQPRCV*      CONVECTIVE PRECIPITATION CONTENT             KG/KG 
!    *PCAPE*       CONVECTIVE AVAILABLE POTENTIAL ENERGY        J/KG 

!    EXTERNALS
!    ---------
!    *CUADJTQS*   FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT
!    *CUADJTQSAD* FOR ADJUSTING T AND Q DUE TO CONDENSATION IN ASCENT

!    MODIFICATIONS
!    -------------
!      P.Lopez    01-Oct-2005  New version of linearized convection
!      P.Lopez    21-Mar-2007  Modified regularization of ZWU for SVs
!      P.Lopez    04-Sep-2007  Revised version to improve match to Tiedtke scheme
!      P.Lopez    09-Nov-2012  Revision to improve match to non-linear version
!      P.Lopez    15-Oct-2015  Added CAPE computation (excl. liquid water loading)
!      P.Lopez    30-Jan-2018  Changed in buoyancy computation to reflect full NL change
!      P.Lopez    28-Feb-2019  Changed entrainment and buoyancy computation to match ref NL

!----------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK
USE YOMCST    , ONLY : RETV, RD, RG, RCPD, RLVTT, RLSTT, RTT, YRCST
USE YOECUMF2  , ONLY : TECUMF2
USE YOECLDP   , ONLY : TECLDP
USE YOETHF    , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 &                     R5ALVCP,  R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE, RTICECU, &
 &                     RTWAT_RTICE_R, RTWAT_RTICECU_R, YRTHF
USE YOMCUMFS  , ONLY : TCUMFS
USE YOMDYNCORE, ONLY : RPLRG

IMPLICIT NONE

TYPE(TCUMFS)      ,INTENT(INOUT) :: YDCUMFS
TYPE(TECLDP)      ,INTENT(INOUT) :: YDECLDP
TYPE(TECUMF2)     ,INTENT(INOUT) :: YDECUMF2
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JKSTART 
LOGICAL           ,INTENT(INOUT) :: LD_LLGO_ON(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQENH5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEO5(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH5(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH5(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSEN5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTU5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQU5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLU5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSENH5(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWU2H5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSUH5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBUOH5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLGLAC5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQPRCV5(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCAPE5(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KLAB(KLON,KLEV) 
LOGICAL           ,INTENT(INOUT) :: LDCUM(KLON) 
LOGICAL           ,INTENT(INOUT) :: LDSC(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCBOT(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KBOTSC(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCTOP(KLON) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KSTUP(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQENH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGEO(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGEOH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAPH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQSEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSENH(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWU2H(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSUH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PBUOH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLGLAC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQPRCV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCAPE(KLON) 

!             LOCAL STORAGE
!             ----- -------

INTEGER(KIND=JPIM) ::  I_KKBASE(KLON)

LOGICAL :: LLGO_ON2(KLON,KLEV)

INTEGER(KIND=JPIM) :: JOK(KLEV)

INTEGER(KIND=JPIM) :: ICALL, IK, JK, JL, JKT2

REAL(KIND=JPRB) :: ZQOLD5(KLON,KLEV),  ZPH5(KLON)
REAL(KIND=JPRB) :: ZMIX5(KLON,KLEV),   ZDMIX5(KLON,KLEV)
REAL(KIND=JPRB) :: ZDZ5(KLON,KLEV)

REAL(KIND=JPRB) :: ZQOLD(KLON),        ZPH(KLON)
REAL(KIND=JPRB) :: ZMIX(KLON)       
REAL(KIND=JPRB) :: ZDZ(KLON),          ZCBASE(KLON)

REAL(KIND=JPRB) :: ZBUOF5(KLON,KLEV),  ZLUOLD5(KLON,KLEV)  
REAL(KIND=JPRB) :: ZFACI5(KLON),       ZFACW5(KLON),        ZFAC5(KLON)
REAL(KIND=JPRB) :: ZCOND5(KLON,KLEV),  ZLNEW5(KLON,KLEV)
REAL(KIND=JPRB) :: ZLMAX5(KLON,KLEV),  ZDQ5(KLON)
REAL(KIND=JPRB) :: ZTVENH5(KLON,KLEV), ZTVUH5(KLON,KLEV),   ZQSU15(KLON),&
                 & ZQSU25(KLON),       ZQSU35(KLON),        ZCOR15(KLON),&
                 & ZALFAW5(KLON),      ZDQSDT5(KLON),       ZDTDP5(KLON),&
                 & ZDP5(KLON),         ZPDIFFTOP5(KLON),    ZPDIFFBOT5(KLON),&
                 & ZSF5(KLON,KLEV),    ZQF5(KLON,KLEV),&
                 & ZCOR25(KLON),       ZEPS5(KLON,KLEV)

REAL(KIND=JPRB) :: ZTU155(KLON,KLEV), ZQU155(KLON,KLEV)
REAL(KIND=JPRB) :: ZQU15(KLON,KLEV),  ZSUH15(KLON,KLEV)
REAL(KIND=JPRB) :: ZLU15(KLON,KLEV),  ZLU25(KLON,KLEV)  

REAL(KIND=JPRB) :: ZBUOF, ZLUOLD, ZCOND, ZLNEW, ZLMAX
REAL(KIND=JPRB) :: ZEPS      ! FRACTIONAL ENTRAINMENT RATE   [M^-1]
REAL(KIND=JPRB) :: ZTVENH    ! ENVIRONMENT VIRTUAL TEMPERATURE AT HALF LEVELS (K)  
REAL(KIND=JPRB) :: ZTVUH     ! UPDRAFT VIRTUAL TEMPERATURE AT HALF LEVELS     (K)
REAL(KIND=JPRB) :: ZQSU, ZCOR, ZDQ, ZFACI, ZFACW, ZFAC, &
                 & ZESDP,ZDQSDT,ZDTDP,ZDP,ZPDIFFTOP,ZPDIFFBOT,&
                 & ZSF,ZQF,ZAW,ZBW,ZALFAW,ZALFAD,ZFACT1,ZFACT3
REAL(KIND=JPRB) :: ZESDP5,ZDMIX25,Z1DPH5,ZCBASE5
REAL(KIND=JPRB) :: ZRG, ZRCPD 
REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL :: LLTEST1(KLON,KLEV)

#include "cuadjtqs.intfb.h"
#include "cuadjtqsad.intfb.h"

!DIR$ VFUNCTION EXPHF
#include "fcttre.func.h"
#include "fcttread.func.h"

!----------------------------------------------------------------------
!     0.           INITIALIZE CONSTANTS AND FIELDS
!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CUPDRAAD',0,ZHOOK_HANDLE)
ASSOCIATE(LREGCV=>YDCUMFS%LREGCV, &
 & RLMIN=>YDECLDP%RLMIN, &
 & ENTRORG2=>YDECUMF2%ENTRORG2, ENTSTPC12=>YDECUMF2%ENTSTPC12, &
 & ENTSTPC22=>YDECUMF2%ENTSTPC22, NJKT22=>YDECUMF2%NJKT22)
ZAW=1.0_JPRB
ZBW=1.0_JPRB
ZRG=1.0_JPRB/RG
ZRCPD=1.0_JPRB/RCPD

JOK(:)=0
JKT2=NJKT22

ZTU155(:,:)=PTU5(:,:)
ZQU155(:,:)=PQU5(:,:)
ZQU15(:,:)=PQU5(:,:)
ZLU15(:,:)=0.0_JPRB
ZLU25(:,:)=0.0_JPRB
ZSUH15(:,:)=PSUH5(:,:)
I_KKBASE(:)=-1
LLGO_ON2(:,:)=.FALSE.

DO JL=KIDIA,KFDIA
  PCAPE5(JL) = 0.0_JPRB
ENDDO

!----------------------------------------------------------------------

!     1.0          DO ASCENT IN SUBCLOUD AND LAYER,
!                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
!                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!                  CHECK FOR BUOYANCY AND SET FLAGS
!                  -------------------------------------
!       ------------------------------------------------------------
!        1.1  DO THE VERTICAL ASCENT UNTIL VELOCITY BECOMES NEGATIVE
!       ------------------------------------------------------------

KLAB(:,1:K_JKSTART)=0

LLGO_ON2(:,K_JKSTART)=LD_LLGO_ON(:)

DO JK=K_JKSTART,JKT2,-1
  DO JL=KIDIA,KFDIA
    IF (LD_LLGO_ON(JL)) THEN
      JOK(JK)=JOK(JK)+1
      ZDZ5(JL,JK)= (PGEOH5(JL,JK) - PGEOH5(JL,JK+1))*ZRG
      IF (K_JKSTART == KLEV-1) THEN
        ZEPS5(JL,JK) = ENTSTPC12 / ((PGEO5(JL,JK)-PGEOH5(JL,KLEV+1))*ZRG*RPLRG) + ENTSTPC22
        ZMIX5(JL,JK)= 0.5_JPRB*RPLRG*ZDZ5(JL,JK)*ZEPS5(JL,JK)
        LLTEST1(JL,JK) = (ZMIX5(JL,JK) > 1.0_JPRB)
        IF (LLTEST1(JL,JK)) ZMIX5(JL,JK) = 1.0_JPRB
        ZDMIX5(JL,JK)= 1.0_JPRB/(1.0_JPRB+ZMIX5(JL,JK))
        ZQF5(JL,JK) = (PQENH5(JL,JK+1) + PQENH5(JL,JK))*0.5_JPRB
        ZSF5(JL,JK) = (PSENH5(JL,JK+1) + PSENH5(JL,JK))*0.5_JPRB
        PQU5(JL,JK) = (PQU5(JL,JK+1)*(1.0_JPRB-ZMIX5(JL,JK))&
                  & + 2.0_JPRB*ZQF5(JL,JK)*ZMIX5(JL,JK))*ZDMIX5(JL,JK)  
        PSUH5(JL,JK)= (PSUH5(JL,JK+1)*(1.0_JPRB-ZMIX5(JL,JK))&
                  & + 2.0_JPRB*ZSF5(JL,JK)*ZMIX5(JL,JK))*ZDMIX5(JL,JK)  
      ELSE
        ZEPS5(JL,JK) = 0.4_JPRB*ENTRORG2*MIN(1.0_JPRB,(PQSEN5(JL,JK)/PQSEN5(JL,KLEV))**3)
        ZMIX5(JL,JK)= ZDZ5(JL,JK)*ZEPS5(JL,JK)
        LLTEST1(JL,JK) = (ZMIX5(JL,JK) > 1.0_JPRB)
        IF (LLTEST1(JL,JK)) ZMIX5(JL,JK) = 1.0_JPRB
        ZQF5(JL,JK) = (PQENH5(JL,JK+1) + PQENH5(JL,JK))*0.5_JPRB
        ZSF5(JL,JK) = (PSENH5(JL,JK+1) + PSENH5(JL,JK))*0.5_JPRB
        PQU5(JL,JK)  = PQU5(JL,JK+1) *(1.0_JPRB-ZMIX5(JL,JK))&
                   & + ZQF5(JL,JK) * ZMIX5(JL,JK)
        PSUH5(JL,JK) = PSUH5(JL,JK+1) *(1.0_JPRB-ZMIX5(JL,JK))&
                   & + ZSF5(JL,JK) * ZMIX5(JL,JK) 
      ENDIF

      ZQOLD5(JL,JK)  = PQU5(JL,JK)
      PTU5(JL,JK) = (PSUH5(JL,JK)-PGEOH5(JL,JK))*ZRCPD
      ZPH5 (JL)   = PAPH5(JL,JK)
    ENDIF
  ENDDO
  
  DO JL=KIDIA,KFDIA
    ZQU15(JL,JK)=PQU5(JL,JK)
    ZSUH15(JL,JK)=PSUH5(JL,JK)
    ZTU155(JL,JK)=PTU5(JL,JK)
    ZQU155(JL,JK)=PQU5(JL,JK)
  ENDDO

  IK=JK
  ICALL=1
  
  CALL CUADJTQS &
   & ( YRTHF, YRCST, KIDIA,    KFDIA,    KLON,    KLEV,&
   & IK,&
   & ZPH5,     PTU5,     PQU5,     LD_LLGO_ON,  ICALL)  

!DIR$ IVDEP
!OCL NOVREC
  DO JL=KIDIA,KFDIA
    IF(LD_LLGO_ON(JL)) THEN

      ZLU15(JL,JK+1)=PLU5(JL,JK+1)

! condensation 

      ZCOND5(JL,JK)=MAX(ZQOLD5(JL,JK)-PQU5(JL,JK),0.0_JPRB)

! microphysics

      ZLUOLD5(JL,JK)=PLU5(JL,JK+1)
      ZLMAX5(JL,JK)=MAX(ZLUOLD5(JL,JK)+ZCOND5(JL,JK),0.0_JPRB)

! Simple microphysics at this stage (to mimic Tiedtke scheme)

      IF (K_JKSTART == KLEV-1) THEN  ! No precip for shallow convection
        ZLNEW5(JL,JK)=MIN(ZLMAX5(JL,JK),5.E-3_JPRB)
      ELSE 
        ZLNEW5(JL,JK)=0.5_JPRB*ZLMAX5(JL,JK)
      ENDIF
      PQPRCV5(JL,JK)=MAX(0.0_JPRB,ZLMAX5(JL,JK)-ZLNEW5(JL,JK))
      PLU5(JL,JK)=ZLNEW5(JL,JK)
      ZLU25(JL,JK)=PLU5(JL,JK)

! frozen water

      PLGLAC5(JL,JK)=ZCOND5(JL,JK)*((1.0_JPRB-FOEALFCU(PTU5(JL,JK))) &
       & -(1.0_JPRB-FOEALFCU(PTU5(JL,JK+1))))
         
! update dry static energy after condensation + freezing

      PSUH5(JL,JK) = RCPD*(PTU5(JL,JK)+RALFDCP*PLGLAC5(JL,JK))+PGEOH5(JL,JK)
      
! Buoyancy on half and full levels
         
      ZTVUH5(JL,JK) = (1.0_JPRB+RETV*PQU5(JL,JK)-PLU5(JL,JK))*PTU5(JL,JK) &
                  & + RALFDCP*PLGLAC5(JL,JK)
      ZTVENH5(JL,JK) = (1.0_JPRB+RETV*PQENH5(JL,JK)) &
                   & * (PSENH5(JL,JK)-PGEOH5(JL,JK))*ZRCPD  
      PBUOH5(JL,JK) = RG*(ZTVUH5(JL,JK)-ZTVENH5(JL,JK))/ZTVENH5(JL,JK)
      ZBUOF5(JL,JK) = (PBUOH5(JL,JK) + PBUOH5(JL,JK+1))*0.5_JPRB

      IF (ZBUOF5(JL,JK) > 0.0_JPRB) THEN
        PCAPE5(JL) = PCAPE5(JL) + ZBUOF5(JL,JK)*ZDZ5(JL,JK)
      ENDIF

! solve kinetic energy equation

      ZDMIX25=1.0_JPRB/(1.0_JPRB+2.0_JPRB*ZBW*ZMIX5(JL,JK))
      PWU2H5(JL,JK) = (PWU2H5(JL,JK+1)*(1.0_JPRB-2.0_JPRB*ZBW*ZMIX5(JL,JK)) &
       & + 2.0_JPRB*ZAW*ZBUOF5(JL,JK)*ZDZ5(JL,JK))*ZDMIX25  

! first layer with liquid water - find exact cloud base

      IF (PLU5(JL,JK) >0.0_JPRB.AND.KLAB(JL,JK+1)==1) THEN

        I_KKBASE(JL)=JK

        IK=JK+1

        ZALFAW5(JL)=FOEALFA(PTU5(JL,IK))
        Z1DPH5=1.0_JPRB/PAPH5(JL,IK)
        ZQSU15(JL)=FOEEWM(PTU5(JL,IK))*Z1DPH5
        ZQSU25(JL)=MIN(0.5_JPRB,ZQSU15(JL))
        ZCOR15(JL)=1.0_JPRB/(1.0_JPRB-RETV*ZQSU25(JL))
        ZQSU35(JL)=ZQSU25(JL)*ZCOR15(JL)
        ZDQ5(JL)=PQU5(JL,IK)-ZQSU35(JL)

        ZFACW5(JL)=R5LES/((PTU5(JL,IK)-R4LES)**2)
        ZFACI5(JL)=R5IES/((PTU5(JL,IK)-R4IES)**2)
        ZFAC5(JL)=ZALFAW5(JL)*ZFACW5(JL) &
         & +(1.-ZALFAW5(JL))*ZFACI5(JL)  

        ZESDP5=FOEEWM(PTU5(JL,IK))*Z1DPH5
        ZCOR25(JL)=1.0_JPRB/(1.0_JPRB-RETV*ZESDP5)
        ZDQSDT5(JL)=ZFAC5(JL)*ZCOR25(JL)*ZQSU35(JL)
        ZDTDP5(JL)=RD*PTU5(JL,IK)*ZRCPD*Z1DPH5
        ZFACT3=1.0_JPRB/(ZDQSDT5(JL)*ZDTDP5(JL))
        ZDP5(JL)=ZDQ5(JL)*ZFACT3
        ZCBASE5=PAPH5(JL,IK)+ZDP5(JL)
        
! chose nearest half level as cloud base

        ZPDIFFTOP5(JL)=ZCBASE5-PAPH5(JL,JK)
        ZPDIFFBOT5(JL)=PAPH5(JL,JK+1)-ZCBASE5
        
        IF(ZPDIFFTOP5(JL) > ZPDIFFBOT5(JL) .AND. &
           & PWU2H5(JL,JK+1) > 0.0_JPRB &
           & .AND. JK < KLEV-1) THEN  
          KLAB(JL,JK)=2
          KLAB(JL,JK+1)=2
          LDSC(JL)   =.TRUE.
          KBOTSC(JL) =JK+1
          KCBOT(JL)  =JK+1
          PLU5(JL,JK+1)= RLMIN
        ELSEIF((ZPDIFFTOP5(JL) <= ZPDIFFBOT5(JL) .OR. &
           & (ZPDIFFTOP5(JL) > ZPDIFFBOT5(JL) .AND. &
           & JK == KLEV-1)) .AND. &
           & PWU2H5(JL,JK) > 0.0_JPRB) THEN  
          KLAB(JL,JK)=2
          LDSC(JL)   =.TRUE.
          KBOTSC(JL) =JK
          KCBOT(JL)  =JK
        ENDIF
      ENDIF

! decide on presence of convection, cloud base and cloud top based on
! kinetic energy

      IF (PWU2H5(JL,JK) < 0.0_JPRB) THEN
        LD_LLGO_ON(JL) = .FALSE.
        IF (PLU5(JL,JK+1) > 0.0_JPRB .AND. JK < K_JKSTART) THEN
          KCTOP(JL) = JK
          KSTUP(JL) = K_JKSTART+1
          LDCUM(JL) = .TRUE.
        ELSE
          LDCUM(JL) = .FALSE.
          KCBOT(JL) = -1
          KCTOP(JL) = -1
        ENDIF
      ELSE
        IF (PLU5(JL,JK) > 0.0_JPRB) THEN
          KLAB(JL,JK) = 2
        ELSE
          KLAB(JL,JK) = 1
        ENDIF
      ENDIF

      LLGO_ON2(JL,JK-1)=LD_LLGO_ON(JL)

    ENDIF
  ENDDO

  IF (JOK(JK) == 0) EXIT 

ENDDO

!     -----------------------------------------------------------
!     -----------------------------------------------------------

!     *****               ADJOINT CALCULATIONS             ******

!     -----------------------------------------------------------
!     -----------------------------------------------------------

!*         0.     INITIALIZATIONS
!                 ---------------

!----------------------------------------------------------------------

!     1.0          DO ASCENT IN SUBCLOUD AND LAYER,
!                  CHECK FOR EXISTENCE OF CONDENSATION LEVEL,
!                  ADJUST T,Q AND L ACCORDINGLY IN *CUADJTQ*,
!                  CHECK FOR BUOYANCY AND SET FLAGS
!                  -------------------------------------
!       ------------------------------------------------------------
!        1.1  DO THE VERTICAL ASCENT UNTIL VELOCITY BECOMES NEGATIVE
!       ------------------------------------------------------------

DO JK=JKT2,K_JKSTART
  IF (JOK(JK) == 0) THEN
    CYCLE
  ELSE
    
    ZCBASE(:)=0.0_JPRB
    ZMIX(:)=0.0_JPRB
    ZDZ(:)=0.0_JPRB
    ZQOLD(:)=0.0_JPRB
    ZPH(:)=0.0_JPRB

!DIR$ IVDEP
!OCL NOVREC

    DO JL=KIDIA,KFDIA

      IF(LLGO_ON2(JL,JK)) THEN
    
        ZLNEW=0.0_JPRB
        ZLMAX=0.0_JPRB
        ZLUOLD=0.0_JPRB
        ZCOND=0.0_JPRB

! first layer with liquid water - find exact cloud base

        IF (JK == I_KKBASE(JL)) THEN

          IK=JK+1
          Z1DPH5=1.0_JPRB/PAPH5(JL,IK)

          ZDP=0.0_JPRB
          ZDQ=0.0_JPRB
          ZDQSDT=0.0_JPRB
          ZDTDP=0.0_JPRB
          ZQSU=0.0_JPRB
          ZCOR=0.0_JPRB
          ZFAC=0.0_JPRB
          ZESDP=0.0_JPRB
          ZFACW=0.0_JPRB
          ZFACI=0.0_JPRB
          ZALFAW=0.0_JPRB
          ZPDIFFBOT=0.0_JPRB
          ZPDIFFTOP=0.0_JPRB

! chose nearest half level as cloud base

          IF (ZPDIFFTOP5(JL) > ZPDIFFBOT5(JL) .AND. &
             & PWU2H5(JL,JK+1) > 0.0_JPRB .AND. JK < KLEV-1) THEN  
            PLU(JL,JK+1)=0.0_JPRB
          ENDIF

          PAPH(JL,JK+1)=PAPH(JL,JK+1)+ZPDIFFBOT
          ZCBASE(JL)=ZCBASE(JL)-ZPDIFFBOT
          ZPDIFFBOT=0.0_JPRB

          ZCBASE(JL)=ZCBASE(JL)+ZPDIFFTOP
          PAPH(JL,JK)=PAPH(JL,JK)-ZPDIFFTOP
          ZPDIFFTOP=0.0_JPRB

          PAPH(JL,IK)=PAPH(JL,IK)+ZCBASE(JL)
          ZDP=ZDP+ZCBASE(JL) 
          ZCBASE(JL)=0.0_JPRB

          ZFACT3=1.0_JPRB/(ZDQSDT5(JL)*ZDTDP5(JL))
          ZDQ   =ZDQ   +ZDP*ZFACT3
          ZDQSDT=ZDQSDT-ZDP*ZDQ5(JL)*ZDTDP5(JL)*(ZFACT3**2)
          ZDTDP =ZDTDP -ZDP*ZDQ5(JL)*ZDQSDT5(JL)*(ZFACT3**2)
          ZDP=0.0_JPRB

          PTU(JL,IK)=PTU(JL,IK)+ZDTDP*ZDTDP5(JL)/PTU5(JL,IK)
          PAPH(JL,IK)=PAPH(JL,IK)-ZDTDP*ZDTDP5(JL)*Z1DPH5
          ZDTDP=0.0_JPRB

          ZQSU=ZQSU+ZDQSDT*ZFAC5(JL)*ZCOR25(JL)
          ZCOR=ZCOR+ZDQSDT*ZFAC5(JL)*ZQSU35(JL)
          ZFAC=ZFAC+ZDQSDT*ZCOR25(JL)*ZQSU35(JL)
          ZDQSDT=0.0_JPRB

          ZESDP=ZESDP+ZCOR*RETV*(ZCOR25(JL)**2) 
          ZCOR=0.0_JPRB

          ZALFAD=FOEALFAAD(ZESDP,PTU5(JL,IK))
          PTU(JL,IK)=PTU(JL,IK)+Z1DPH5 &
           & *FOEEWMAD(ZESDP,PTU5(JL,IK),ZALFAD)  
          PAPH(JL,IK)=PAPH(JL,IK)-ZESDP*FOEEWM(PTU5(JL,IK)) &
           & *(Z1DPH5**2)  
          ZESDP=0.0_JPRB

          ZFACW  = ZFACW  + ZFAC* ZALFAW5(JL)
          ZALFAW = ZALFAW + ZFAC*(ZFACW5(JL)-ZFACI5(JL))
          ZFACI  = ZFACI  + ZFAC* (1.-ZALFAW5(JL))
          ZFAC=0.0_JPRB
         
          PTU(JL,IK)=PTU(JL,IK)-ZFACI*2.0_JPRB*ZFACI5(JL) &
           & /(PTU5(JL,IK)-R4IES)   
          ZFACI=0.0_JPRB
        
          PTU(JL,IK)=PTU(JL,IK)-ZFACW*2.0_JPRB*ZFACW5(JL) &
           & /(PTU5(JL,IK)-R4LES)   
          ZFACW=0.0_JPRB
        
          PQU(JL,IK)=PQU(JL,IK)+ZDQ
          ZQSU=ZQSU-ZDQ
          ZDQ=0.0_JPRB

          ZCOR=ZCOR+ZQSU*ZQSU25(JL)
          ZQSU=ZQSU*ZCOR15(JL)
        
          ZQSU=ZQSU+ZCOR*RETV*(ZCOR15(JL)**2)
          ZCOR=0.0_JPRB
        
          IF (ZQSU15(JL) > 0.5_JPRB) ZQSU=0.0_JPRB
        
          ZALFAD=FOEALFAAD(ZQSU,PTU5(JL,IK))
          PTU(JL,IK)=PTU(JL,IK)+Z1DPH5*FOEEWMAD(ZQSU,PTU5(JL,IK),ZALFAD) 
          PAPH(JL,IK)=PAPH(JL,IK)-ZQSU*FOEEWM(PTU5(JL,IK)) &
           & *(Z1DPH5**2)  
          ZQSU=0.0_JPRB

          PTU(JL,IK)=PTU(JL,IK)+FOEALFAAD(ZALFAW,PTU5(JL,IK))
          ZALFAW=0.0_JPRB

        ENDIF
        
! solve kinetic energy equation

        ZDMIX25=1.0_JPRB/(1.0_JPRB+2.0_JPRB*ZBW*ZMIX5(JL,JK))

        ZBUOF=0.0_JPRB

        PWU2H(JL,JK+1)=PWU2H(JL,JK+1)+PWU2H(JL,JK) &
         & *(1.0_JPRB-2.0_JPRB*ZBW*ZMIX5(JL,JK))*ZDMIX25  
        ZMIX(JL)=ZMIX(JL)-PWU2H(JL,JK)*2.0_JPRB*ZBW &
         & *(PWU2H5(JL,JK+1)+PWU2H5(JL,JK))*ZDMIX25  
        ZDZ(JL)=ZDZ(JL)+PWU2H(JL,JK)*2.0_JPRB*ZAW &
         & *ZBUOF5(JL,JK)*ZDMIX25  
        ZBUOF=ZBUOF+PWU2H(JL,JK)*2.0_JPRB*ZAW &
         & *ZDZ5(JL,JK)*ZDMIX25  
        PWU2H(JL,JK)=0.0_JPRB
      
        IF (ZBUOF5(JL,JK) > 0.0_JPRB) THEN
          ZBUOF   = ZBUOF   + PCAPE(JL) * ZDZ5(JL,JK)
          ZDZ(JL) = ZDZ(JL) + PCAPE(JL) * ZBUOF5(JL,JK)
        ENDIF

! Reduce buoyancy perturbation if required
        IF (LREGCV) THEN
          ZBUOF=ZBUOF*0.33_JPRB
        ENDIF

! Buoyancy on half and full levels

        PBUOH(JL,JK  )=PBUOH(JL,JK  )+0.5_JPRB*ZBUOF
        PBUOH(JL,JK+1)=PBUOH(JL,JK+1)+0.5_JPRB*ZBUOF
        ZBUOF=0.0_JPRB
  
        ZTVUH=0.0_JPRB
        ZTVENH=0.0_JPRB

        ZTVUH  = ZTVUH  + PBUOH(JL,JK) * RG / ZTVENH5(JL,JK)
        ZTVENH = ZTVENH - PBUOH(JL,JK) * RG * ZTVUH5(JL,JK) / (ZTVENH5(JL,JK)**2)  
        PBUOH(JL,JK)=0.0_JPRB

        PSENH(JL,JK)=PSENH(JL,JK)+ZTVENH &
         & *(1.0_JPRB+RETV*PQENH5(JL,JK))*ZRCPD  
        PGEOH(JL,JK)=PGEOH(JL,JK)-ZTVENH &
         & *(1.0_JPRB+RETV*PQENH5(JL,JK))*ZRCPD  
        PQENH(JL,JK)=PQENH(JL,JK)+ZTVENH &
         & *RETV*(PSENH5(JL,JK)-PGEOH5(JL,JK))*ZRCPD  
        ZTVENH=0.0_JPRB

        PQU(JL,JK)=PQU(JL,JK)+ZTVUH*RETV*PTU5(JL,JK)
        PLU(JL,JK)=PLU(JL,JK)-ZTVUH*PTU5(JL,JK)
        PTU(JL,JK)=PTU(JL,JK)+ZTVUH*(1.0_JPRB+RETV*PQU5(JL,JK)-ZLU25(JL,JK))
        PLGLAC(JL,JK)=PLGLAC(JL,JK)+RALFDCP*ZTVUH
        ZTVUH=0.0_JPRB

! update dry static energy after condensation + freezing

        PTU(JL,JK)=PTU(JL,JK)+PSUH(JL,JK)*RCPD
        PLGLAC(JL,JK)=PLGLAC(JL,JK)+PSUH(JL,JK)*RCPD*RALFDCP
        PGEOH(JL,JK)=PGEOH(JL,JK)+PSUH(JL,JK)
        PSUH(JL,JK)=0.0_JPRB
  
! frozen water

        PTU(JL,JK)=PTU(JL,JK)-ZCOND5(JL,JK) &
         & *FOEALFCUAD(PLGLAC(JL,JK),PTU5(JL,JK))   
        PTU(JL,JK+1)=PTU(JL,JK+1)+ZCOND5(JL,JK) &
         & *FOEALFCUAD(PLGLAC(JL,JK),PTU5(JL,JK+1))  
        ZCOND=ZCOND+PLGLAC(JL,JK) &
         & *( (1.0_JPRB-FOEALFCU(PTU5(JL,JK))) &
         & -  (1.0_JPRB-FOEALFCU(PTU5(JL,JK+1))) )  
        PLGLAC(JL,JK)=0.0_JPRB  
    
! microphysics

        ZLNEW=ZLNEW+PLU(JL,JK)
        PLU(JL,JK)=0.0_JPRB

        IF (ZLMAX5(JL,JK)-ZLNEW5(JL,JK) > 0.0_JPRB) THEN      
          ZLMAX=ZLMAX+PQPRCV(JL,JK)
          ZLNEW=ZLNEW-PQPRCV(JL,JK)
        ENDIF
        PQPRCV(JL,JK)=0.0_JPRB
      
! Simple microphysics at this stage (to mimic Tiedtke scheme)

        IF (K_JKSTART == KLEV-1) THEN  ! No precip for shallow convection
          IF (ZLMAX5(JL,JK) < 5.E-3_JPRB) THEN
            ZLMAX = ZLMAX + ZLNEW
          ENDIF
        ELSE 
          ZLMAX = ZLMAX + 0.5_JPRB*ZLNEW
        ENDIF
        ZLNEW=0.0_JPRB

        IF (ZLUOLD5(JL,JK)+ZCOND5(JL,JK) > 0.0_JPRB) THEN
          ZLUOLD=ZLUOLD+ZLMAX
          ZCOND=ZCOND+ZLMAX
        ENDIF
        ZLMAX=0.0_JPRB
      
        PLU(JL,JK+1)=PLU(JL,JK+1)+ZLUOLD
        ZLUOLD=0.0_JPRB

! condensation

        IF (ZQOLD5(JL,JK)-PQU5(JL,JK) > 0.0_JPRB) THEN
          ZQOLD(JL)=ZQOLD(JL)+ZCOND
          PQU(JL,JK)=PQU(JL,JK)-ZCOND
        ENDIF
        ZCOND=0.0_JPRB

! ENDIF LLGO_ON2(JL,JK)
      ENDIF

      ZPH5(JL)= PAPH5(JL,JK)
      LD_LLGO_ON(JL)=LLGO_ON2(JL,JK)

! LOOP OVER JL
    ENDDO

    IK=JK
    ICALL=1

    CALL CUADJTQSAD &
     & ( KIDIA,    KFDIA,    KLON,    KLEV,&
     & IK,&
     & ZPH5,     ZTU155,   ZQU155,&
     & ZPH,      PTU,      PQU,      LD_LLGO_ON,   ICALL)  

    DO JL=KIDIA,KFDIA
      IF (LLGO_ON2(JL,JK)) THEN
        ZSF=0.0_JPRB
        ZQF=0.0_JPRB
        ZEPS=0.0_JPRB

        PAPH(JL,JK)=PAPH(JL,JK)+ZPH(JL) 
        ZPH(JL)=0.0_JPRB
 
        PSUH(JL,JK)=PSUH(JL,JK)+PTU(JL,JK)*ZRCPD
        PGEOH(JL,JK)=PGEOH(JL,JK)-PTU(JL,JK)*ZRCPD
        PTU(JL,JK)=0.0_JPRB
    
        PQU(JL,JK)=PQU(JL,JK)+ZQOLD(JL)
        ZQOLD(JL)=0.0_JPRB

        IF (K_JKSTART == KLEV-1) THEN
          PSUH(JL,JK+1)=PSUH(JL,JK+1)+PSUH(JL,JK) &
                     & *(1.0_JPRB-ZMIX5(JL,JK))*ZDMIX5(JL,JK)   
          ZSF=ZSF+PSUH(JL,JK)*2.0_JPRB*ZMIX5(JL,JK)*ZDMIX5(JL,JK) 
          ZMIX(JL)=ZMIX(JL)+PSUH(JL,JK)*(2.0_JPRB*ZSF5(JL,JK) &
                & -PSUH5(JL,JK+1)-ZSUH15(JL,JK))*ZDMIX5(JL,JK)   
          PSUH(JL,JK)=0.0_JPRB
    
          PQU(JL,JK+1)=PQU(JL,JK+1)+PQU(JL,JK) &
                    & *(1.0_JPRB-ZMIX5(JL,JK))*ZDMIX5(JL,JK)   
          ZMIX(JL)=ZMIX(JL)+PQU(JL,JK)*(2.0_JPRB*ZQF5(JL,JK) &
                & -PQU5(JL,JK+1)-ZQU15(JL,JK))*ZDMIX5(JL,JK)   
          ZQF=ZQF+PQU(JL,JK)*2.0_JPRB*ZMIX5(JL,JK)*ZDMIX5(JL,JK) 
          PQU(JL,JK)=0.0_JPRB
    
          PSENH(JL,JK+1)=PSENH(JL,JK+1)+ZSF*0.5_JPRB
          PSENH(JL,JK  )=PSENH(JL,JK  )+ZSF*0.5_JPRB
          ZSF=0.0_JPRB
  
          PQENH(JL,JK+1)=PQENH(JL,JK+1)+ZQF*0.5_JPRB
          PQENH(JL,JK  )=PQENH(JL,JK  )+ZQF*0.5_JPRB
          ZQF=0.0_JPRB
  
          IF (LLTEST1(JL,JK)) ZMIX(JL) = 0.0_JPRB

          ZEPS   =ZEPS   +ZMIX(JL)*0.5_JPRB*RPLRG*ZDZ5(JL,JK)
          ZDZ(JL)=ZDZ(JL)+ZMIX(JL)*0.5_JPRB*RPLRG*ZEPS5(JL,JK)
          ZMIX(JL)=0.0_JPRB
 
          ZFACT1 = ENTSTPC12 / ((PGEO5(JL,JK) - PGEOH5(JL,KLEV+1))**2 * ZRG * RPLRG)
          PGEO(JL,JK)      = PGEO(JL,JK)      - ZEPS * ZFACT1
          PGEOH(JL,KLEV+1) = PGEOH(JL,KLEV+1) + ZEPS * ZFACT1
          ZEPS=0.0_JPRB
        ELSE
          PSUH(JL,JK+1)=PSUH(JL,JK+1)+PSUH(JL,JK)*(1.0_JPRB-ZMIX5(JL,JK))
          ZSF=ZSF+PSUH(JL,JK)*ZMIX5(JL,JK)
          ZMIX(JL)=ZMIX(JL)+PSUH(JL,JK)*(ZSF5(JL,JK)-PSUH5(JL,JK+1))   
          PSUH(JL,JK)=0.0_JPRB

          PQU(JL,JK+1)=PQU(JL,JK+1)+PQU(JL,JK)*(1.0_JPRB-ZMIX5(JL,JK))  
          ZQF=ZQF+PQU(JL,JK)*ZMIX5(JL,JK)
          ZMIX(JL)=ZMIX(JL)+PQU(JL,JK)*(ZQF5(JL,JK)-PQU5(JL,JK+1))
          PQU(JL,JK)=0.0_JPRB
    
          PSENH(JL,JK+1)=PSENH(JL,JK+1)+ZSF*0.5_JPRB
          PSENH(JL,JK  )=PSENH(JL,JK  )+ZSF*0.5_JPRB
          ZSF=0.0_JPRB
  
          PQENH(JL,JK+1)=PQENH(JL,JK+1)+ZQF*0.5_JPRB
          PQENH(JL,JK  )=PQENH(JL,JK  )+ZQF*0.5_JPRB
          ZQF=0.0_JPRB

          IF (LLTEST1(JL,JK)) ZMIX(JL) = 0.0_JPRB

          ZEPS=ZEPS+ZMIX(JL)*ZDZ5(JL,JK)
          ZDZ(JL)=ZDZ(JL)+ZMIX(JL)*ZEPS5(JL,JK)
          ZMIX(JL)=0.0_JPRB

          IF (PQSEN5(JL,JK)/PQSEN5(JL,KLEV) < 1.0_JPRB) THEN 
            PQSEN(JL,JK)  =PQSEN(JL,JK)  +ZEPS*3.0_JPRB*ZEPS5(JL,JK)/PQSEN5(JL,JK)
            PQSEN(JL,KLEV)=PQSEN(JL,KLEV)-ZEPS*3.0_JPRB*ZEPS5(JL,JK)/PQSEN5(JL,KLEV)
          ENDIF     
          ZEPS=0.0_JPRB
        ENDIF

        PGEOH(JL,JK  )=PGEOH(JL,JK  )+ZDZ(JL)*ZRG
        PGEOH(JL,JK+1)=PGEOH(JL,JK+1)-ZDZ(JL)*ZRG
        ZDZ(JL)=0.0_JPRB
              
! ENDIF LLGO_ON2(JL,JK)
      ENDIF
! LOOP OVER JL
    ENDDO

! ENDIF JOK(JK)=0
  ENDIF

! LOOP OVER JK
ENDDO

DO JL=KIDIA,KFDIA
  PCAPE(JL) = 0.0_JPRB
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CUPDRAAD',1,ZHOOK_HANDLE)
END SUBROUTINE CUPDRAAD
