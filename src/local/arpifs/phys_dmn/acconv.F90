!OPTIONS XOPT(NOEVAL)
SUBROUTINE ACCONV (YDCST,YDML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
 !-----------------------------------------------------------------------
 ! - INPUT  2D .
 & PAPHIF,PAPRS,PAPRSF,PCVGQ,&
 & PDELP,PQ,PRDELP,PT,PU,PV,&
 ! - INPUT  1D .
 & PTS,&
 ! - OUTPUT 2D .
 & PDIFCQ,PDIFCS,PFPLCL,PFPLCN,PSTRCU,PSTRCV)

!**** *ACCONV * - MASS-FLUX DEEP CONVECTION SCHEME.

!     Subject.
!     -------
!     - MASS-FLUX DEEP CONVECTION SCHEME, COMPUTING CONVECTIVE
!       PRECIPITATIONS AND "SUB-GRID SCALE" ASSOCIATED FLUXES .

!     - COMPUTATION OF ENTHALPY FLUXES LINKED TO CONVECT. PRECIPITATIONS
!              (WATER AND SNOW) .

!**   Interface.
!     ----------
!        *CALL* *ACCONV*

! -   ARGUMENTS FOR INPUT.
!     --------------------

! - DIMENSIONAL PARAMETERS OF PHYSICS.

! KIDIA      : START OF HORIZONTAL LOOP (IST in CPG).
! KFDIA      : END OF HORIZONTAL LOOP (IEND in CPG).
! KLON       : HORIZONTAL DIMENSIONAL.
! KTDIA      : START OF THE VERTICAL LOOP IN THE PHYSIC (1 IN GENERAL).
! KLEV       : END OF VERTICAL LOOP AND VERTICAL DIMENSION "FULL LEVEL".
!              (NFLEVG in CPG).

! - THE NAME OF THE PHYSICAL VARIABLES

! - 2D (0:KLEV) .

! PAPRS      : PRESSURE ON HALF LEVELS.

! - 2D (1:KLEV) .

! PAPHIF     : GEOPOTENTIAL ON FULL LEVELS.
! PAPRSF     : PRESSURE ON FULL LEVELS.
! PCVGQ      : CONVERGENCE D'HUMIDITE (CONDITION DE "KUO").
! PDELP      : LAYER THICKNESS IN PRESSURE UNITS.
! PQ         : SPECIFIC HUMIDITY OF WATER VAPOUR.
! PRDELP     : INVERSE OF LAYER THICKNESS IN PRESSURE UNITS.
! PT         : TEMPERATURE.
! PU         : X COMPONENT OF WIND.
! PV         : Y COMPONENT OF WIND.

! - 1D

! PTS        : SURFACE TEMPERATURE.

!-----------------------------------------------------------------------

! -   ARGUMENTS FOR OUTPUT.
!     ---------------------

! - 2D (0:KLEV) .

! PDIFCQ     : CONVECTIVE FLUX OF SPECIFIC HUMIDITY (NOT RAIN/SNOW).
! PDIFCS     : CONVECTIVE FLUX OF ENTHALPY (NOT RAIN/SNOW).
! PFPLCL     : CONVECTIVE PRECIPITATION AS RAIN.
! PFPLCN     : CONVECTIVE PRECIPITATION AS SNOW.
! PSTRCU     : CONVECTIVE FLUX OF MOMENTUM "U".
! PSTRCV     : CONVECTIVE FLUX OF MOMENTUM "V".

!-----------------------------------------------------------------------

! -   IMPLICIT ARGUMENTS.
!     ---------------------

!-----------------------------------------------------------------------

!     Externals.
!     ----------

!     Method.
!     -------

!     Author
!     -------
!      97-01, M. Janiskova.

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!-----------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE YOMCST   , ONLY : TCST

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVGQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFCQ(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFCS(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLCL(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFPLCN(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRCU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRCV(KLON,0:KLEV) 

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: I_KNLAB(KLON,KLEV),I_KNND(KLON)

REAL(KIND=JPRB) :: ZSNP(KLON,0:KLEV),ZLHN(KLON,0:KLEV),ZLHL(KLON,0:KLEV)&
 & ,ZFORM(KLON,0:KLEV),ZZLHP(KLON,0:KLEV)  
REAL(KIND=JPRB) :: ZCP(KLON,KLEV),ZLHE(KLON,KLEV),ZDSE(KLON,KLEV)&
 & ,ZQW(KLON,KLEV),ZTW(KLON,KLEV),ZMSE(KLON,KLEV)&
 & ,ZMSEU(KLON,KLEV),ZF(KLON,KLEV),ZMFC(KLON,KLEV)&
 & ,ZCORU(KLON,KLEV),ZCORV(KLON,KLEV),ZCORQ(KLON,KLEV)&
 & ,ZCORS(KLON,KLEV),ZFQSC(KLON,KLEV)&
 & ,ZFUSC(KLON,KLEV),ZFVSC(KLON,KLEV),ZFSSC(KLON,KLEV)&
 & ,ZQDN(KLON,KLEV),ZPOID(KLON,KLEV),ZIPOI(KLON,KLEV)&
 & ,ZMSEC(KLON,KLEV)  
REAL(KIND=JPRB) :: ZLHSB(KLON),ZTN(KLON),ZQN(KLON),ZS1(KLON),ZS2(KLON)&
 & ,ZPOII(KLON),ZALF(KLON),ZLN(KLON),ZUN(KLON),ZVN(KLON)&
 & ,ZMSEN(KLON),ZDFMDQ(KLON),ZDFMDU(KLON),ZDFMDV(KLON)&
 & ,ZDFMDS(KLON),ZDSEN(KLON),ZENTR(KLON),ZDETR(KLON)&
 & ,ZPOIL(KLON)  

INTEGER(KIND=JPIM) :: ISUM, ITOP, JIT, JLEV, JLON

LOGICAL :: LL_LCLOSN

REAL(KIND=JPRB) :: ZCORIN, ZCPSMD, ZCPVMD, ZCPVMS, ZCPVMW, ZCPWMD,&
 & ZCVMCD, ZDELQ, ZDELT, ZDELTA, ZDETC, ZDQW, &
 & ZEPS1, ZEPS2, ZEPS3, ZEPS4, ZEPS5, ZEPS6, &
 & ZEPS7, ZESP, ZEW, ZFCORQ, ZFMDQ, ZFMDS, ZFMDU, &
 & ZFMDV, ZFTOTQ, ZFTOTS, ZGDT, ZGDTI, ZINSTA, &
 & ZIV, ZIV1, ZMSEP, ZMSEPR, ZRVMD, ZTI, ZTQ, &
 & ZZQW, ZZVAL  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "fcttrm.ycst.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACCONV',0,ZHOOK_HANDLE)
ASSOCIATE(SCO=>YDML_PHY_MF%YRPHY0%SCO, USDMLT=>YDML_PHY_MF%YRPHY0%USDMLT, ECMNP=>YDML_PHY_MF%YRPHY0%ECMNP, &
 & TSPHY=>YDML_PHY_MF%YRPHY2%TSPHY, &
 & NBITER=>YDML_PHY_MF%YRPHY%NBITER, LNEIGE=>YDML_PHY_MF%YRPHY%LNEIGE)
!-----------------------------------------------------------------------

!*
!     ------------------------------------------------------------------
!     I - AUXILIARY CONSTANTS.

LL_LCLOSN=.FALSE.

ZRVMD=YDCST%RV-YDCST%RD
ZCPVMD=YDCST%RCPV-YDCST%RCPD
ZCPVMW=YDCST%RCPV-YDCST%RCW
ZCPVMS=YDCST%RCPV-YDCST%RCS
ZCPWMD=YDCST%RCW-YDCST%RCPD
ZCPSMD=YDCST%RCS-YDCST%RCPD
ZCVMCD=YDCST%RCPV-YDCST%RCPD

!*
!     ------------------------------------------------------------------
!     II - COMPUTATION OF DERIVED PARAMETERS, SECURITY CONSTANTS
!     AND INVERSE OF A CRITICAL ARGUMENT (TO AVOID "UNDERFLOW" IN
!     LIQUID (OR ICE) WATER CALCULATIONS FOR SUSTENTATION IN CLOUDY
!     ASCENTS).

ZGDT=YDCST%RG*TSPHY
ZGDTI=1.0_JPRB/ZGDT

ZEPS1=1.E-01_JPRB
ZEPS2=1.E+01_JPRB
ZEPS3=1.E-01_JPRB
ZEPS4=1.E-04_JPRB

ZEPS5=1.E-06_JPRB
ZEPS6=1.E-07_JPRB
ZEPS7=1.E-05_JPRB

ZDETC=0.0_JPRB
ZCORIN=0.0_JPRB

!*
!     ------------------------------------------------------------------
!     III - COMPUTATION OF VALUE OF ATMOSPHERIC Cp,
!     SETTING ALL FLUXES TO ZERO (FOR THE CASE WITHOUT CONVECTION)
!     AND PRELIMINARY CALCULATIONS OF DRY STATIC ENERGY, OF EFFECTIVE
!     LATENT HEAT (DEPENDING UPON THE OPTION WHETHER SURFACE PRESSURE IS
!     VARYING OR NOT WITH THE SURFACE WATER FLUXES), OF LEVEL PRESSURE
!     THICKNESSES DIVIDED BY G*DT AND OF THEIR INVERSES AT ALL LEVELS AS
!     WELL AS OF IMPOSED PROPORTIONS OF SNOW PRECIPITATION IN FUNCTION
!     OF TEMPERATURE AND PRESSURE.

! - TEMPORARY(S) 2D (1:KLEV) .

! ZCP        : VALUE OF ATMOSPHERIC Cp
! ZPOID     : DP/(RG*DT) FOR A GIVEN LEVEL AND A GIVEN TIME STEP.
! ZIPOI     : INVERSE OF ZPOID.
! ZDSE      : LOCAL ARRAY FOR DRY STATIC ENERGY.
! ZMSE      : LOCAL ARRAY FOR MOIST STATIC ENERGY.
! ZLHE      : LOCAL ARRAY FOR EFFECTIVE LATENT HEATS.

! - TEMPORARY(S) 2D (0:KLEV) .

! ZSNP      : PROPORTION OF SNOW PRECIPITATION.
! ZLHL      : MINUS THE LATENT HEAT OF CONSERVATION OF LIQUID ENTHALPY.
! ZLHN      : MINUS THE LATENT HEAT OF CONSERVATION OF SOLID ENTHALPY.

!     COMPUTATION OF VALUE OF ATMOSPHERIC Cp

DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA

    ZCP(JLON,JLEV)=YDCST%RCPD*(1.0_JPRB-PQ(JLON,JLEV))+YDCST%RCPV*PQ(JLON,JLEV)
  ENDDO
ENDDO

!     SNOW OPTION DEPENDENT CALCULATIONS.

DO JLON=KIDIA,KFDIA

  IF (LNEIGE) THEN
    ZSNP(JLON,KTDIA-1)=1.0_JPRB-MIN(1.0_JPRB,USDMLT*MAX(0.0_JPRB,PT(JLON,KTDIA)&
     & -YDCST%RTT)**2/MAX(1.E-03_JPRB,PAPRS(JLON,KTDIA-1)))  
  ELSE
    ZSNP(JLON,KTDIA-1)=0.0_JPRB
  ENDIF

ENDDO

!     SNOW OPTION DEPENDENT CALCULATIONS.

DO JLEV=KTDIA,KLEV-1
  DO JLON=KIDIA,KFDIA

    IF (LNEIGE) THEN
      ZSNP(JLON,JLEV)=1.0_JPRB-MIN(1.0_JPRB,USDMLT*MAX(0.0_JPRB,0.5_JPRB*(PT(JLON,JLEV)&
       & +PT(JLON,JLEV+1))-YDCST%RTT)**2/PAPRS(JLON,JLEV))  
    ELSE
      ZSNP(JLON,JLEV)=0.0_JPRB
    ENDIF

  ENDDO
ENDDO

!     SNOW OPTION DEPENDENT CALCULATIONS.

DO JLON=KIDIA,KFDIA

  IF (LNEIGE) THEN
    ZSNP(JLON,KLEV)=1.0_JPRB-MIN(1.0_JPRB,USDMLT*MAX(0.0_JPRB,PTS(JLON)-YDCST%RTT)**2&
     & /PAPRS(JLON,KLEV))  
  ELSE
    ZSNP(JLON,KLEV)=0.0_JPRB
  ENDIF

ENDDO

!     "FLUX LATENT HEAT" ARRAYS INITIALIZATION.

! --- COMPUTATION AT THE TOP
!     ----------------------
DO JLON=KIDIA,KFDIA
  ZTQ=ZCVMCD*PT(JLON,KTDIA)
  ZLHL(JLON,0)=-(FOLH(PT(JLON,KTDIA),0.0_JPRB)-ZTQ)
  ZLHN(JLON,0)=-(FOLH(PT(JLON,KTDIA),1.0_JPRB)-ZTQ)
ENDDO

! --- COMPUTATION FOR LEVELS
!     ----------------------
DO JLEV = KTDIA , KLEV-1
  DO JLON = KIDIA , KFDIA
    ZTI=0.5_JPRB*(PT(JLON,JLEV)+PT(JLON,JLEV+1))
    ZTQ=ZCVMCD*ZTI
    ZLHL(JLON,JLEV)=-(FOLH(ZTI,0.0_JPRB)-ZTQ)
    ZLHN(JLON,JLEV)=-(FOLH(ZTI,1.0_JPRB)-ZTQ)
  ENDDO
ENDDO

! --- SURFACE COMPUTATION
!     --------------------
DO JLON=KIDIA,KFDIA
  ZTQ=ZCVMCD*PTS(JLON)
  ZLHL(JLON,KLEV)=-(FOLH(PTS(JLON),0.0_JPRB)-ZTQ)
  ZLHN(JLON,KLEV)=-(FOLH(PTS(JLON),1.0_JPRB)-ZTQ)
ENDDO

DO JLON=KIDIA,KFDIA

!     CALL TO THE DEFINED FUNCTION TO COMPUTE THE REFERENCE OF EFFECTIVE
!     LATENT HEAT AT THE SURFACE.

! - TEMPORARY(S) 1D .

! ZLHSB     : BUDGET LATENT HEAT AS COMPUTED AT THE SURFACE.

  ZLHSB(JLON)=-ZLHL(JLON,KLEV)-ZSNP(JLON,KLEV)*(ZLHN(JLON,KLEV)&
   & -ZLHL(JLON,KLEV))  

!     SETTING TO ZERO FLUXES AT THE TOP.

  PFPLCL(JLON,KTDIA-1)=0.0_JPRB
  PFPLCN(JLON,KTDIA-1)=0.0_JPRB
  PDIFCQ(JLON,KTDIA-1)=0.0_JPRB
  PDIFCS(JLON,KTDIA-1)=0.0_JPRB
  PSTRCU(JLON,KTDIA-1)=0.0_JPRB
  PSTRCV(JLON,KTDIA-1)=0.0_JPRB
ENDDO

!     COMPUTATIONS AT ALL LEVELS. BEWARE, THE DEFINED FUNCTION FOLH WILL
!     BE CALLED WITH VALUES THAT ARE NOT NECESSARILY 0 OR 1 FOR ITS
!     SECOND ARGUMENT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    ZPOID(JLON,JLEV)=PDELP(JLON,JLEV)*ZGDTI
    ZIPOI(JLON,JLEV)=PRDELP(JLON,JLEV)*ZGDT
    ZDSE(JLON,JLEV)=ZCP(JLON,JLEV)*PT(JLON,JLEV)+PAPHIF(JLON,JLEV)
    ZDELTA=0.5_JPRB*(ZSNP(JLON,JLEV-1)+ZSNP(JLON,JLEV))
    ZLHE(JLON,JLEV)= FOLH (PT(JLON,JLEV),ZDELTA)+(ZCPWMD+ZDELTA&
     & *(ZCPSMD-ZCPWMD))*(PT(JLON,JLEV)&
     & -PTS(JLON))  

    ZMSE(JLON,JLEV)=ZCP(JLON,JLEV)*PT(JLON,JLEV)&
     & +ZLHE(JLON,JLEV)*PQ(JLON,JLEV)+PAPHIF(JLON,JLEV)  
    ZMSEC(JLON,JLEV)=ZCP(JLON,JLEV)*PT(JLON,JLEV)&
     & +ZLHSB(JLON)*PQ(JLON,JLEV)+PAPHIF(JLON,JLEV)  

!     SETTING FLUXES TO ZERO.

    PFPLCL(JLON,JLEV)=0.0_JPRB
    PFPLCN(JLON,JLEV)=0.0_JPRB
    PDIFCQ(JLON,JLEV)=0.0_JPRB
    PDIFCS(JLON,JLEV)=0.0_JPRB
    PSTRCU(JLON,JLEV)=0.0_JPRB
    PSTRCV(JLON,JLEV)=0.0_JPRB
  ENDDO
ENDDO

!*
!     ------------------------------------------------------------------
!     IV - COMPUTATION OF PROFILE OF MASS FLUX AND NORMALIZATION
!     FACTOR FOR THIS PROFILE CONNECTED WITH THE "KUO-TYPE" CLOSURE

!     ACTIVE CHARACTERISTICS OF THE CLOUD AT THE BOTTOM

! - TEMPORARY(S) 2D (0:KLEV) .

! ZFORM     : PROFILE OF MASS FLUX AND LATER MASS FLUX TIME DT.

! - TEMPORARY(S) 1D .

! ZS1       : INTEGRAL OF HUMIDITY CONVERGENCE
! ZS2       : INTEGRAL OF BUOYANCY

IF (.NOT.LL_LCLOSN) THEN
  DO JLON=KIDIA,KFDIA
    ZMSEU(JLON,KLEV)=ZMSE(JLON,KLEV)
    ZFORM(JLON,KLEV)=0.0_JPRB

    ZS1(JLON)=0.0_JPRB
    ZS2(JLON)=0.0_JPRB
  ENDDO

!       COMPUTATION FOR LEVELS

  DO JLEV=KLEV-1,KTDIA,-1
    DO JLON=KIDIA,KFDIA
      ZINSTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,ZMSEU(JLON,JLEV+1)-ZMSE(JLON,JLEV)&
       & -ZCORIN))  
      ZMSEU(JLON,JLEV)=ZINSTA*ZMSEU(JLON,JLEV+1)&
       & +(1.0_JPRB-ZINSTA)*ZMSE(JLON,JLEV)  
      ZFORM(JLON,JLEV)=ZINSTA*SQRT(MAX(0.0_JPRB,ZMSEU(JLON,JLEV)&
       & +ZMSEU(JLON,JLEV+1)&
       & -ZMSE(JLON,JLEV)-ZMSE(JLON,JLEV+1))&
       & *PAPRS(JLON,JLEV)&
       & /(PT(JLON,JLEV)+PT(JLON,JLEV+1)))  
      ZPOII(JLON)=MIN(1.0_JPRB,MAX(0.0_JPRB,(1.0_JPRB+ZDETC*PDELP(JLON,JLEV+1))&
       & *ZFORM(JLON,JLEV)-ZFORM(JLON,JLEV+1))&
       & /MAX(ZEPS7,ZFORM(JLON,JLEV)))  
      ZMSEP=ZPOII(JLON)*ZMSE(JLON,JLEV+1)&
       & +(1.0_JPRB-ZPOII(JLON))*ZMSEU(JLON,JLEV+1)  
      ZINSTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,ZMSEP-ZMSE(JLON,JLEV)-ZCORIN))
      ZMSEU(JLON,JLEV)=ZINSTA*ZMSEP+(1.0_JPRB-ZINSTA)*ZMSE(JLON,JLEV)
      ZFORM(JLON,JLEV)=ZINSTA*SQRT(MAX(0.0_JPRB,ZMSEU(JLON,JLEV)&
       & +ZMSEU(JLON,JLEV+1)&
       & -ZMSE(JLON,JLEV)-ZMSE(JLON,JLEV+1))&
       & *PAPRS(JLON,JLEV)&
       & /(PT(JLON,JLEV)+PT(JLON,JLEV+1)))  
      ZS1(JLON)=ZS1(JLON)+PCVGQ(JLON,JLEV+1)*PDELP(JLON,JLEV+1)*ZINSTA
      ZS2(JLON)=ZS2(JLON)+(PQ(JLON,JLEV+1)-PQ(JLON,JLEV))*ZFORM(JLON,JLEV)
    ENDDO
  ENDDO

  DO JLON=KIDIA,KFDIA
    ZALF(JLON)=(MAX(0.0_JPRB,ZS1(JLON))/MAX(ZEPS6,ZS2(JLON)))&
     & *MAX(0.0_JPRB,SIGN(1.0_JPRB,ZS2(JLON)-ZEPS6))  
  ENDDO
ENDIF

!*
!     ------------------------------------------------------------------
!     V - BOTTOM INITIALISATION. THE CLOUD PROFILE WILL EVENTUALY BE
!     CHARACTERISED BY ITS TOTAL HUMIDITY AND ITS DETRAINING DRY STATIC
!     ENERGY AS WELL AS BY A STABILITY INDEX (ONE AT BOTH ENDS OF ANY
!     CONVECTIVE SLAB (FROM MODEL LAYER TO THE ADJACENT ONE) AND ZERO
!     ELSEWHERE).

DO JLON=KIDIA,KFDIA

!     INITIALIZATION OF THE VARIABLES IN THE NEWTON LOOP, 
!     THE ITERATIVE TW AND QW HAVING TO EVOLVE

! - TEMPORARY 2D (1:KLEV)

! ZQW        : SPECIFIC WET BULB HUMIDITY.
! ZTW        : WET BULB TEMPERATURE.

!  505 CONTINUE

  ZQW(JLON,KLEV)=PQ(JLON,KLEV)
  ZTW(JLON,KLEV)=PT(JLON,KLEV)
ENDDO

!     NEWTON LOOP FOR THE "WET BULB" POINT.

DO JIT=1,NBITER
  DO JLON=KIDIA,KFDIA

!     SNOW OPTION DEPENDENT COMPUTATIONS.

    IF (LNEIGE) THEN
      ZDELTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,YDCST%RTT-ZTW(JLON,KLEV)))
    ELSE
      ZDELTA=0.0_JPRB
    ENDIF

!     SATURATION CALCULATIONS USING THE DEFINED FUNCTIONS

    ZEW= FOEW (ZTW(JLON,KLEV),ZDELTA)
    ZESP=ZEW/PAPRSF(JLON,KLEV)
    ZZQW= FOQS (ZESP)
    ZDQW= FODQS (ZZQW,ZESP, FODLEW (ZTW(JLON,KLEV),ZDELTA))

!     INCREMENTATIONS.

    ZDELQ=(ZZQW-ZQW(JLON,KLEV))*ZCP(JLON,KLEV)&
     & /(ZCP(JLON,KLEV)+ZLHSB(JLON)*ZDQW)  
    ZDELT=-ZDELQ*ZLHSB(JLON)/ZCP(JLON,KLEV)
    ZQW(JLON,KLEV)=ZQW(JLON,KLEV)+ZDELQ
    ZTW(JLON,KLEV)=ZTW(JLON,KLEV)+ZDELT
  ENDDO

ENDDO

DO JLON=KIDIA,KFDIA

!     THERMODYNAMIC CHARACTERISTICS OF THE CLOUD.

! - TEMPORARY(S) 1D .

! ZTN       : TEMPERATURE OF THE CLOUDY ASCENT.
! ZQN       : WATER VAPOUR SPECIFIC HUMIDITY OF THE CLOUDY ASCENT.
! ZLN       : CONDENSED SPECIFIC HUMIDITY OF THE CLOUDY ASCENT.
! ZUN       : U-WIND COMPONENT OF THE CLOUDY ASCENT.
! ZVN       : V-WIND COMPONENT OF THE CLOUDY ASCENT.
! ZMSEN     : MOIST STATIC ENERGY OF THE CLOUDY ASCENT.
! ZDSEN     : DRY STATIC ENERGY OF THE CLOUDY ASCENT.

  ZTN(JLON)=ZTW(JLON,KLEV)
  ZQN(JLON)=ZQW(JLON,KLEV)
  ZUN(JLON)=PU(JLON,KLEV)
  ZVN(JLON)=PV(JLON,KLEV)
  ZMSEN(JLON)=ZCP(JLON,KLEV)*ZTN(JLON)&
   & +ZLHSB(JLON)*ZQN(JLON)+PAPHIF(JLON,KLEV)  
  ZDSEN(JLON)=ZCP(JLON,KLEV)*ZTN(JLON)+PAPHIF(JLON,KLEV)
  ZLN(JLON)=0.0_JPRB

!     ACTIVE OTHER CHARACTERISTICS OF THE CLOUD.

! - TEMPORARY(S) 2D (1:KLEV)

! ZMCF      : MASS FLUX IN CLOUD
! ZQDN      : TOTAL CLOUD HUMIDITY

  ZMFC(JLON,KLEV)=0.0_JPRB
  ZQDN(JLON,KLEV)=ZQN(JLON)
  I_KNLAB(JLON,KLEV)=0
  I_KNND(JLON)=0

ENDDO

!*
!     ------------------------------------------------------------------
!     VI - FIRST VERTICAL LOOP (UPWARDS) INCLUDING THE SATURATED
!     ADIABATIC TYPE CALCULATION OF THE CLOUD PROFILE.

!     INITIALIZATION OF THE PARAMETER THAT SHALL EVENTUALLY HELP
!     AVOIDING UNNECESSARY COMPUTATIONS IN THE UPPER PART OF THE
!     ATMOSPHERE AND START OF THE VERTICAL LOOP.

ITOP=KLEV+1

!     OPTIONAL COMPUTATION ACCORDING CHOSEN CLOSURE ASSUMPTION
!     COMPUTATION OF VALUES OF DETRAINMENT AND ENTRAINMENT, TEST
!     OF STABILITY, COMPUTATION OF THE SECOND TERM OF THE BASIC EQUATION,
!     COMPUTATION OF VALUES u,v,q,h IN CLOUD

! - TEMPORARY(S) 2D (1:KLEV)

!   ZCORU   : THE SECOND TERM OF THE BASIC EQUATION FOR U-WIND COMP.
!   ZCORV   : THE SECOND TERM OF THE BASIC EQUATION FOR V-WIND COMP.
!   ZCORQ   : THE SECOND TERM OF THE BASIC EQUATION FOR SPEC.HUMIDITY
!   ZCORS   : THE SECOND TERM OF THE BASIC EQUATION FOR DRY STATIC
!             ENERGY

DO JLEV=KLEV-1,KTDIA,-1

  DO JLON=KIDIA,KFDIA

!     OPTIONAL COMPUTATION ACCORDING CHOSEN CLOSURE ASSUMPTION
!     IF LCLOSN=TRUE NEW CLOSURE ASSUMPTION OF SIMPLIFIED PHYSICS
!     WILL BE USED

    IF (LL_LCLOSN) THEN
      ZF(JLON,JLEV+1)=PCVGQ(JLON,JLEV+1)/MAX(ZEPS5,PQ(JLON,JLEV+1))
    ELSE
      ZF(JLON,JLEV+1)=ZALF(JLON)*(ZFORM(JLON,JLEV)&
       & -ZFORM(JLON,JLEV+1))/PDELP(JLON,JLEV+1)  
    ENDIF

    ZF(JLON,JLEV+1)=ZF(JLON,JLEV+1)/(1.0_JPRB+MAX(0.0_JPRB,ZF(JLON,JLEV+1))*TSPHY)
    ZMFC(JLON,JLEV)=MAX(0.0_JPRB,ZMFC(JLON,JLEV+1)&
     & +ZF(JLON,JLEV+1)*PDELP(JLON,JLEV+1))  

!     COMPUTATION OF DETRAINMENT AND ENTAINMENT

    ZENTR(JLON)=(MAX(0.0_JPRB,ZF(JLON,JLEV+1)+ZDETC*ZMFC(JLON,JLEV)))&
     & *PDELP(JLON,JLEV+1)  
    ZENTR(JLON)=MIN(ZENTR(JLON),ZMFC(JLON,JLEV))
    ZDETR(JLON)=ZENTR(JLON)+ZMFC(JLON,JLEV+1)-ZMFC(JLON,JLEV)

    ZPOIL(JLON)=ZENTR(JLON)/MAX(ZEPS4,ZMFC(JLON,JLEV))

!     STABILTY TEST

    ZMSEPR=ZPOIL(JLON)*ZMSEC(JLON,JLEV+1)+(1.0_JPRB-ZPOIL(JLON))*ZMSEN(JLON)

    I_KNLAB(JLON,JLEV)=NINT(MAX(0.0_JPRB,-SIGN(1.0_JPRB,&
     & -(ZMSEPR-ZMSEC(JLON,JLEV))))&
     & *MAX(MAX(0.0_JPRB,SIGN(1.0_JPRB,&
     & ZDSE(JLON,JLEV)-ZDSE(JLON,JLEV+1))),&
     & MAX(0.0_JPRB,-SIGN(1.0_JPRB,0.0_JPRB-ZMFC(JLON,JLEV+1)))))        

    I_KNLAB(JLON,JLEV+1)=MAX(I_KNLAB(JLON,JLEV+1),I_KNLAB(JLON,JLEV))
    I_KNND(JLON)=MAX(I_KNND(JLON),I_KNLAB(JLON,JLEV))

    ZMFC(JLON,JLEV)=I_KNLAB(JLON,JLEV)*ZMFC(JLON,JLEV)
    ZDETR(JLON)=I_KNLAB(JLON,JLEV)*ZDETR(JLON)&
     & +(1-I_KNLAB(JLON,JLEV))*ZMFC(JLON,JLEV+1)  
    ZENTR(JLON)=I_KNLAB(JLON,JLEV)*ZENTR(JLON)

!     COMPUTATION OF THE SECOND TERM OF THE BASIC EQUATION

    ZCORU(JLON,JLEV+1)=ZENTR(JLON)*PU(JLON,JLEV+1)-ZDETR(JLON)*ZUN(JLON)
    ZCORV(JLON,JLEV+1)=ZENTR(JLON)*PV(JLON,JLEV+1)-ZDETR(JLON)*ZVN(JLON)
    ZCORQ(JLON,JLEV+1)=ZENTR(JLON)*PQ(JLON,JLEV+1)&
     & -ZDETR(JLON)*(ZQN(JLON)+ZLN(JLON))  
    ZCORS(JLON,JLEV+1)=ZENTR(JLON)*ZDSE(JLON,JLEV+1)&
     & -ZDETR(JLON)&
     & *(ZDSEN(JLON)-ZLHSB(JLON)*ZLN(JLON))  

!     COMPUTATION OF VALUES u,v,h,q,T IN THE CLOUD

    ZUN(JLON)=I_KNLAB(JLON,JLEV)*(ZPOIL(JLON)*PU(JLON,JLEV+1)&
     & +(1.0_JPRB-ZPOIL(JLON))*ZUN(JLON))&
     & +(1-I_KNLAB(JLON,JLEV))*PU(JLON,JLEV)  
    ZVN(JLON)=I_KNLAB(JLON,JLEV)*(ZPOIL(JLON)*PV(JLON,JLEV+1)&
     & +(1.0_JPRB-ZPOIL(JLON))*ZVN(JLON))&
     & +(1-I_KNLAB(JLON,JLEV))*PV(JLON,JLEV)  
    ZQW(JLON,JLEV)=I_KNLAB(JLON,JLEV)*(ZPOIL(JLON)&
     & *PQ(JLON,JLEV+1)&
     & +(1.0_JPRB-ZPOIL(JLON))*ZQN(JLON))&
     & +(1-I_KNLAB(JLON,JLEV))*PQ(JLON,JLEV)  
    ZDSEN(JLON)=I_KNLAB(JLON,JLEV)*(ZPOIL(JLON)&
     & *ZDSE(JLON,JLEV+1)&
     & +(1.0_JPRB-ZPOIL(JLON))*ZDSEN(JLON))&
     & +(1-I_KNLAB(JLON,JLEV))*ZDSE(JLON,JLEV)  
    ZMSEN(JLON)=I_KNLAB(JLON,JLEV)*ZMSEPR+(1-I_KNLAB(JLON,JLEV))*ZMSEC(JLON,JLEV)

!     INITIALIZATION OF THE VARIABLES IN THE NEWTON LOOP, 
!     THE ITERATIVE TW AND QW HAVING TO EVOLVE.

    ZTW(JLON,JLEV)=(ZDSEN(JLON)-PAPHIF(JLON,JLEV))/ZCP(JLON,JLEV)

  ENDDO

!     NEWTON LOOP FOR THE "WET BULB" POINT.

  DO JIT=1,NBITER
    DO JLON=KIDIA,KFDIA

!     SNOW OPTION DEPENDENT COMPUTATIONS.

      IF (LNEIGE) THEN
        ZDELTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,YDCST%RTT-ZTW(JLON,JLEV)))
      ELSE
        ZDELTA=0.0_JPRB
      ENDIF

!     SATURATION CALCULATIONS USING THE DEFINED FUNCTIONS AND TEMPORARY

! - TEMPORARY 1D .

      ZEW= FOEW (ZTW(JLON,JLEV),ZDELTA)
      ZESP=ZEW/PAPRSF(JLON,JLEV)
      ZZQW= FOQS (ZESP)
      ZDQW= FODQS (ZZQW,ZESP, FODLEW (ZTW(JLON,JLEV),ZDELTA))

!     INCREMENTATIONS.

      ZDELQ=(ZZQW-ZQW(JLON,JLEV))*ZCP(JLON,JLEV)&
       & /(ZCP(JLON,JLEV)+ZLHSB(JLON)*ZDQW)  
      ZDELT=-ZDELQ*ZLHSB(JLON)/ZCP(JLON,JLEV)
      ZQW(JLON,JLEV)=ZQW(JLON,JLEV)+ZDELQ
      ZTW(JLON,JLEV)=ZTW(JLON,JLEV)+ZDELT
    ENDDO

  ENDDO

  DO JLON=KIDIA,KFDIA

    ZLN(JLON)=MAX(0.0_JPRB,ZLN(JLON)&
     & *(1.0_JPRB-ZPOIL(JLON)/(1.0_JPRB+ZPOIL(JLON)))&
     & *(1.0_JPRB-(PAPHIF(JLON,JLEV)-PAPHIF(JLON,JLEV+1))/ECMNP)&
     & -(ZQW(JLON,JLEV)-ZQN(JLON)&
     & -(ZPOIL(JLON)/(1.0_JPRB+ZPOIL(JLON)))&
     & *(PQ(JLON,JLEV+1)-ZQN(JLON))))*I_KNLAB(JLON,JLEV)  
    ZQN(JLON)=ZQW(JLON,JLEV)
    ZQDN(JLON,JLEV)=ZQN(JLON)+ZLN(JLON)
    ZTN(JLON)=ZTW(JLON,JLEV)

  ENDDO

!     VERIFICATION OF THE PRESENCE OF AT LEAST ONE CLOUD ALONG THE
!     "JLON" VECTOR AND MODIFICATION OF ITOP IF NECESSARY.

  ISUM=0
  DO JLON=KIDIA,KFDIA
    ISUM=ISUM+I_KNLAB(JLON,JLEV)+I_KNLAB(JLON,JLEV+1)
  ENDDO

  IF (ISUM /= 0) THEN
    ITOP=JLEV
  ENDIF

ENDDO

!     CALCULATION AT THE TOP OF THE ATMOSPHERE

DO JLON=KIDIA,KFDIA

  ZCORU(JLON,KTDIA)=-ZMFC(JLON,KTDIA)*ZUN(JLON)
  ZCORV(JLON,KTDIA)=-ZMFC(JLON,KTDIA)*ZVN(JLON)
  ZCORQ(JLON,KTDIA)=-ZMFC(JLON,KTDIA)*(ZQN(JLON)+ZLN(JLON))
  ZCORS(JLON,KTDIA)=-ZMFC(JLON,KTDIA)*(ZDSEN(JLON)-ZLHSB(JLON)*ZLN(JLON))

ENDDO

! ----------------------------------------------------------------------
!     QUITTING IN CASE OF NO CONVECTION EVERYWHERE. THE "IF THEN /
!     ENDIF" ENCOMPASSED SEQUENCE WILL NOT BE INDENTED GIVEN ITS LENGTH
!     AND ITS SPECIAL CHARACTER (PSEUDO "RETURN" STATEMENT).

IF (ITOP < KLEV+1) THEN
! ----------------------------------------------------------------------

!*
!     ------------------------------------------------------------------
!     VII - SECOND VERTICAL LOOP (DOWNWARDS) TO SOLVE THE IMPLICIT
!     ALGORITHM FOR COMPENSATING SUBSIDENCE ADVECTION. ONE WILL AT THE
!     SAME TIME STORE PROVISIONAL FINAL VALUES FOR HUMIDITY, MOIST
!     STATIC ENERGY AND FOR THE TWO WIND COMPONENTS.
!     STRAIGHTFORWARD ELIMINATION-TYPE SOLUTION OF THE LINEAR SYSTEM
!     (BI-DIAGONAL MATRIX).

! - TEMPORARY(S) 2D (1:KLEV) .

! ZFQSC     : FINAL Q VALUE (AT FIRST ONLY COMPENSATING SUBSIDENCE).
! ZFSSC     : FINAL S VALUE (AT FIRST ONLY COMPENSATING SUBSIDENCE).
! ZFUSC     : FINAL U VALUE (AT FIRST ONLY COMPENSATING SUBSIDENCE).
! ZFVSC     : FINAL V VALUE (AT FIRST ONLY COMPENSATING SUBSIDENCE).

! - TEMPORARY(S) 1D .

! ZDFMDS   :  HALF THE MASS FLUX TIME THE JUMP IN S.
! ZDFMDU   :  HALF THE MASS FLUX TIME THE JUMP IN U.
! ZDFMDV   :  HALF THE MASS FLUX TIME THE JUMP IN V.
! ZDFMDQ   :  HALF THE MASS FLUX TIME THE JUMP IN Q.

!     SPECIAL CASE OF THE FIRST LEVEL.

  DO JLON=KIDIA,KFDIA

    ZIV1=PDELP(JLON,ITOP)/TSPHY
    ZIV=1.0_JPRB/(ZIV1+ZMFC(JLON,ITOP))
    ZDFMDU(JLON)=0.5_JPRB*ZMFC(JLON,ITOP)*(PU(JLON,ITOP+1)-PU(JLON,ITOP))
    ZDFMDV(JLON)=0.5_JPRB*ZMFC(JLON,ITOP)*(PV(JLON,ITOP+1)-PV(JLON,ITOP))
    ZDFMDQ(JLON)=0.5_JPRB*ZMFC(JLON,ITOP)*(PQ(JLON,ITOP+1)-PQ(JLON,ITOP))
    ZDFMDS(JLON)=0.5_JPRB*ZMFC(JLON,ITOP)*(ZDSE(JLON,ITOP+1)-ZDSE(JLON,ITOP))

    ZFUSC(JLON,ITOP)=(PU(JLON,ITOP)*ZIV1-ZCORU(JLON,ITOP)-ZDFMDU(JLON))*ZIV
    ZFVSC(JLON,ITOP)=(PV(JLON,ITOP)*ZIV1-ZCORV(JLON,ITOP)-ZDFMDV(JLON))*ZIV
    ZFQSC(JLON,ITOP)=(PQ(JLON,ITOP)*ZIV1-ZCORQ(JLON,ITOP)-ZDFMDQ(JLON))*ZIV
    ZFSSC(JLON,ITOP)=(ZDSE(JLON,ITOP)*ZIV1-ZCORS(JLON,ITOP)-ZDFMDS(JLON))*ZIV

  ENDDO

!     VERTICAL LOOP WITH USING THE RESULTS FROM PREVIOUS UPPER LEVEL.

  DO JLEV=ITOP+1,KLEV-1
    DO JLON=KIDIA,KFDIA

      ZIV1=PDELP(JLON,JLEV)/TSPHY
      ZIV=1.0_JPRB/(ZIV1+ZMFC(JLON,JLEV))
      ZFMDU=0.5_JPRB*ZMFC(JLON,JLEV)*(PU(JLON,JLEV+1)-PU(JLON,JLEV))
      ZFMDV=0.5_JPRB*ZMFC(JLON,JLEV)*(PV(JLON,JLEV+1)-PV(JLON,JLEV))
      ZFMDQ=0.5_JPRB*ZMFC(JLON,JLEV)*(PQ(JLON,JLEV+1)-PQ(JLON,JLEV))
      ZFMDS=0.5_JPRB*ZMFC(JLON,JLEV)*(ZDSE(JLON,JLEV+1)-ZDSE(JLON,JLEV))

      ZFUSC(JLON,JLEV)=(PU(JLON,JLEV)*ZIV1-ZCORU(JLON,JLEV)&
       & +ZMFC(JLON,JLEV-1)*ZFUSC(JLON,JLEV-1)&
       & -(ZFMDU-ZDFMDU(JLON)))*ZIV  
      ZFVSC(JLON,JLEV)=(PV(JLON,JLEV)*ZIV1-ZCORV(JLON,JLEV)&
       & +ZMFC(JLON,JLEV-1)*ZFVSC(JLON,JLEV-1)&
       & -(ZFMDV-ZDFMDV(JLON)))*ZIV  
      ZFQSC(JLON,JLEV)=(PQ(JLON,JLEV)*ZIV1-ZCORQ(JLON,JLEV)&
       & +ZMFC(JLON,JLEV-1)*ZFQSC(JLON,JLEV-1)&
       & -(ZFMDQ-ZDFMDQ(JLON)))*ZIV  
      ZFSSC(JLON,JLEV)=(ZDSE(JLON,JLEV)*ZIV1-ZCORS(JLON,JLEV)&
       & +ZMFC(JLON,JLEV-1)*ZFSSC(JLON,JLEV-1)&
       & -(ZFMDS-ZDFMDS(JLON)))*ZIV  

      ZDFMDU(JLON)=ZFMDU
      ZDFMDV(JLON)=ZFMDV
      ZDFMDQ(JLON)=ZFMDQ
      ZDFMDS(JLON)=ZFMDS

    ENDDO
  ENDDO

!     SPECIAL CASE OF THE LAST LEVEL.

  IF (ITOP /= KLEV) THEN
    DO JLON=KIDIA,KFDIA
      ZIV=1.0_JPRB/(PDELP(JLON,KLEV)/TSPHY)
      ZFUSC(JLON,KLEV)=PU(JLON,KLEV)+(ZMFC(JLON,KLEV-1)&
       & *ZFUSC(JLON,KLEV-1)-ZCORU(JLON,KLEV)&
       & +ZDFMDU(JLON))*ZIV  
      ZFVSC(JLON,KLEV)=PV(JLON,KLEV)+(ZMFC(JLON,KLEV-1)&
       & *ZFVSC(JLON,KLEV-1)-ZCORV(JLON,KLEV)&
       & +ZDFMDV(JLON))*ZIV  
      ZFQSC(JLON,KLEV)=PQ(JLON,KLEV)+(ZMFC(JLON,KLEV-1)&
       & *ZFQSC(JLON,KLEV-1)-ZCORQ(JLON,KLEV)&
       & +ZDFMDQ(JLON))*ZIV  
      ZFSSC(JLON,KLEV)=ZDSE(JLON,KLEV)+(ZMFC(JLON,KLEV-1)&
       & *ZFSSC(JLON,KLEV-1)-ZCORS(JLON,KLEV)&
       & +ZDFMDS(JLON))*ZIV  

    ENDDO
  ENDIF

!*
!     ------------------------------------------------------------------
!     VIII - FLUX COMPUTATION

!     DIFFUSIVE HUMIDITY FLUX AND INTERPOLATED LATENT HEAT OF FLUXES.

! - TEMPORARY(S) 2D (0:KLEV) .

! ZZLHP     : FICTITIOUS PRECIPITATION FLUX LATENT HEAT.

  DO JLON=KIDIA,KFDIA
    ZZLHP(JLON,ITOP-1)=-ZLHL(JLON,ITOP-1)-ZSNP(JLON,ITOP-1)&
     & *(ZLHN(JLON,ITOP-1)-ZLHL(JLON,ITOP-1))  
  ENDDO
  DO JLEV=ITOP,KLEV-1
    DO JLON=KIDIA,KFDIA
      PDIFCQ(JLON,JLEV)=ZMFC(JLON,JLEV)&
       & *0.5_JPRB*(ZFQSC(JLON,JLEV)&
       & +ZFQSC(JLON,JLEV+1)-ZQDN(JLON,JLEV)&
       & -ZQDN(JLON,JLEV+1))/YDCST%RG  
      ZZLHP(JLON,JLEV)=-ZLHL(JLON,JLEV)-ZSNP(JLON,JLEV)&
       & *(ZLHN(JLON,JLEV)-ZLHL(JLON,JLEV))  
    ENDDO
  ENDDO

  DO JLON=KIDIA,KFDIA
    ZZLHP(JLON,KLEV)=ZLHSB(JLON)
  ENDDO

!     FINAL COMPUTATION OF FLUXES WITH SECURITY AGAINST NEGATIVE
!     PRECIPITATION.

  DO JLEV=ITOP,KLEV
    DO JLON=KIDIA,KFDIA

!     DIFFUSIVE MOMENTUM FLUXES.

      PSTRCU(JLON,JLEV)=PSTRCU(JLON,JLEV-1)&
       & -ZPOID(JLON,JLEV)*(ZFUSC(JLON,JLEV)-PU(JLON,JLEV))  
      PSTRCV(JLON,JLEV)=PSTRCV(JLON,JLEV-1)&
       & -ZPOID(JLON,JLEV)*(ZFVSC(JLON,JLEV)-PV(JLON,JLEV))  

!     HUMIDITY FLUXES (TOTAL AND PRECIPITATING).

      ZFTOTQ=PDIFCQ(JLON,JLEV-1)+PFPLCL(JLON,JLEV-1)&
       & +PFPLCN(JLON,JLEV-1)&
       & -(ZPOID(JLON,JLEV)*(ZFQSC(JLON,JLEV)-PQ(JLON,JLEV)) )  
      ZFCORQ=MAX(0.0_JPRB,PDIFCQ(JLON,JLEV)-ZFTOTQ)
      PFPLCL(JLON,JLEV)=(ZFCORQ-(PDIFCQ(JLON,JLEV)-ZFTOTQ))&
       & *(1.0_JPRB-ZSNP(JLON,JLEV))  
      PFPLCN(JLON,JLEV)=(ZFCORQ-(PDIFCQ(JLON,JLEV)-ZFTOTQ))*ZSNP(JLON,JLEV)

!     ENTHALPY FLUXES (TOTAL AND DIFFUSIVE).

      ZFTOTS=PDIFCS(JLON,JLEV-1)-ZZLHP(JLON,JLEV-1)&
       & *(PFPLCL(JLON,JLEV-1)+PFPLCN(JLON,JLEV-1))&
       & -ZPOID(JLON,JLEV)*(ZFSSC(JLON,JLEV)-ZDSE(JLON,JLEV))  
      PDIFCS(JLON,JLEV)=ZFTOTS-ZLHSB(JLON)*ZFCORQ&
       & +ZZLHP(JLON,JLEV)*(PFPLCL(JLON,JLEV)+PFPLCN(JLON,JLEV))  

    ENDDO
  ENDDO

!     SETTING ALL RESULTS TO ZERO IN CASE OF NEGATIVE PRECIPITATIONS AT
!     THE SURFACE .

  DO JLON=KIDIA,KFDIA
    ZZVAL=PFPLCL(JLON,KLEV)+PFPLCN(JLON,KLEV)-SCO
    I_KNND(JLON)=I_KNND(JLON)*MAX(0.0_JPRB,-SIGN(1.0_JPRB,0.0_JPRB-ZZVAL))
  ENDDO
  DO JLEV=ITOP,KLEV
    DO JLON=KIDIA,KFDIA
      PFPLCL(JLON,JLEV)=I_KNND(JLON)*PFPLCL(JLON,JLEV)
      PFPLCN(JLON,JLEV)=I_KNND(JLON)*PFPLCN(JLON,JLEV)
      PDIFCQ(JLON,JLEV)=I_KNND(JLON)*PDIFCQ(JLON,JLEV)
      PDIFCS(JLON,JLEV)=I_KNND(JLON)*PDIFCS(JLON,JLEV)
      PSTRCU(JLON,JLEV)=I_KNND(JLON)*PSTRCU(JLON,JLEV)
      PSTRCV(JLON,JLEV)=I_KNND(JLON)*PSTRCV(JLON,JLEV)
    ENDDO
  ENDDO

!*
!     RETURN FROM THE PSEUDO-RETURN SITUATED AFTER LOOP 670.

! -------------
ENDIF
! -------------

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACCONV',1,ZHOOK_HANDLE)
END SUBROUTINE ACCONV
