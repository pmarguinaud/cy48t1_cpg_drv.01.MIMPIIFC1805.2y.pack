!OPTIONS XOPT(NOEVAL)
SUBROUTINE ACDIFSP (YDCST,YDML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
 !-----------------------------------------------------------------------
 ! - INPUT  2D .
 & PAPHI,PAPHIF,PDELP,PRDELP,PT,PQ,&
 & PU,PV,PKTROV,PKUROV,PXTROV,PXUROV,&
 ! - INPUT  1D .
 & PCDROV,PCHROV,PXDROV,PXHROV,&
 & PLSM,PNEIJ,PQS,PTS,PNTS,PTS1,PQS1,&
 ! - OUTPUT 2D .
 & PDIFTQ,PDIFTS,PSTRTU,PSTRTV,&
 ! - OUTPUT 1D .
 & PFCLL,PFCLN,PFCS,PFEVL,PFEVN,PLHS)

!**** *ACDIFSP * - VERTICAL TURBULENT DIFFUSION - SIMPLIFIED

!     Subject.
!     --------
!     - COMPUTATION OF VERTICAL TURBULENT DIFFUSION INCLUDING A NEGATIVE
!       HUMIDITY CORRECTION : VERTICAL TURBULENT FLUXES, CALL TO THE
!       CHARNOCK FORMULA, THERMAL RADIATION FLUXES' CORRECTION AND ALL
!       CORRESPONDING FLUXES AT THE SURFACE OR IN THE SOIL .
!       SIMPLIFICATION OF ORIGINAL ACDIFUS IS BASED ON ASSUMPTION
!       OF CONSTANT FIELDS FOR Ts and Qs.

!**   Interface.
!     ----------
!        *CALL* *ACDIFSP*

!-----------------------------------------------------------------------

! -   ARGUMENTS FOR INPUT
!     -------------------

! - DIMENSIONAL PARAMETERS OF PHYSICS.

! KIDIA      : START OF HORIZONTAL LOOP (IST in CPG)
! KFDIA      : END OF HORIZONTAL LOOP.
! KLON       : HORIZONTAL DIMENSION.
! KTDIA      : START OF THE VERTICAL LOOP IN THE PHYSICS (1 IN GENERAL)
! KLEV       : END OF VERTICAL LOOP AND VERTCAL DIMENSION "FULL LEVEL"
!              (NFLEVG in CPG)

! - THE NAME OF THE PHYSICAL VARIABLES.

! - 2D (0:KLEV) .

! PAPHI      : GEOPOTENTIAL ON HALF LEVELS.
! PKTROV     : VERTICAL EXCHANGE COEFFICIENT OF T AND Q IN KG/(M*M*S).
! PKUROV     : VERTICAL EXCHANGE COEFFICIENT OF U AND V IN KG/(M*M*S).
! PXTROV     : "ANTI-FIBRILLATION" MULTIPLICATOR OF PKTROV.
! PXUROV     : "ANTI-FIBRILLATION" MULTIPLICATOR OF PKUROV.

! - 2D (1:KLEV) .

! PAPHIF     : GEOPOTENTIAL ON FULL LEVELS.
! PDELP      : LAYER THICKNESS IN PRESSURE UNITS.
! PQ         : SPECIFIC HUMIDITY OF WATER VAPOR.
! PRDELP     : INVERSE VALUE OF LAYER THICKNESS (OF PDELP).
! PT         : TEMPERATURE.
! PU         : X - COMPONENT OF WIND.
! PV         : Y - COMPONENT OF WIND.

! - 1D (PROGNOSTIC) .

! PTS        : SURFACE TEMPERATURE.
! PNTS       : SURFACE TEMPERATURE AT TIME t+dt.
! PQS        : SURFACE SPECIFIC HUMIDITY.
! PQS1       : SURFACE SPECIFIC HUMIDITY AT TIME t+dt.

! - 1D (GEOGRAPHIQUE) .

! PLSM       : LAND/SEA MASK.

! - 1D (DIAGNOSTIQUE) .

! PCDROV     : PCD RENORME EN DENSITE FOIS VITESSE.
! PCHROV     : PCH RENORME EN DENSITE FOIS VITESSE.
! PNEIJ      : FRACTION OF SOIL COVERED BY SNOW.
! PQS        : SURFACE SPECIFIC HUMIDITY.
! PXDROV     : "ANTI-FIBRILLATION" MULTIPLICATOR OF PCDROV.
! PXHROV     : "ANTI-FIBRILLATION" MULTIPLICATOR OF PCHROV.

!-----------------------------------------------------------------------

! -   ARGUMENTS FOR OUTPUT.
!     ---------------------

! - 2D (0:KLEV) .

! PDIFTQ     : TURBULENT FLUX (INC. Q NEGATIVE) OF SPECIFIC HUMIDITY.
! PDIFTS     : TURBULENT FLUX OF ENTHALPY (OR DRY STATIC ENERGY).
! PSTRTU     : TURBULENT FLUX OF MOMENTUM "U".
! PSTRTV     : TURBULENT FLUX OF MOMENTUM "V".

! - 1D (DIAGNOSTIC) .

! PFCLL      : LATENT HEAT FLUX OVER LIQUID WATER (OR WET SOIL).
! PFCLN      : LATENT HEAT FLUX OVER SNOW (OR ICE).
! PFCS       : SENSIBLE HEAT FLUX AT SURFACE LEVEL.
! PFEVL      : WATER VAPOR FLUX OVER LIQUID WATER (OR WET SOIL).
! PFEVN      : WATER VAPOR FLUX OVER SNOW (OR ICE).
! PLHS       : LATENT HEAT AT THE SURFACE.

!-----------------------------------------------------------------------

! -   IMPLICIT ARGUMENTS
!     -------------------

!-----------------------------------------------------------------------

!     Externals.
!     ----------

!     Methode.
!     --------

!     Authors.
!     --------
!      89-12, J.F. Geleyn.

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
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN) :: YDML_PHY_MF
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHI(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKTROV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKUROV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXTROV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXUROV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCDROV(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCHROV(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXDROV(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXHROV(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNEIJ(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNTS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS1(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQS1(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTQ(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIFTS(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRTV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCLL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCLN(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCS(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFEVL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFEVN(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLHS(KLON) 

!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZKSRV(KLON,0:KLEV),ZXTRO(KLON,0:KLEV),ZXURO(KLON,0:KLEV)&
 & ,ZFCQNG(KLON,0:KLEV)&
 & ,ZDSE(KLON,KLEV),ZIPOI(KLON,KLEV),ZPOID(KLON,KLEV)&
 & ,ZN1(KLON,KLEV),ZN2(KLON,KLEV),ZN3(KLON,KLEV),ZN4(KLON,KLEV)&
 & ,ZQCOR(KLON,KLEV),ZSUB1(KLON,KLEV),ZSUB2(KLON,KLEV)&
 & ,ZCP(KLON,KLEV)&
 & ,ZCHRS(KLON),ZDSES1(KLON),ZDLLS(KLON),ZDLLW(KLON)&
 & ,ZDLS(KLON),ZLHS(KLON),ZLHW(KLON),ZPN(KLON)&
 & ,ZXDRO(KLON),ZXSRO(KLON),ZCPS1(KLON)  

INTEGER(KIND=JPIM) :: JLEV, JLON

REAL(KIND=JPRB) :: ZCPVMD, ZCPVMS, ZCPVMW, ZDLDQ, ZELIM1, ZELIM2,&
 & ZEPS1, ZFENT, ZFQ, ZGDT, ZGDTI, ZMUL1, ZMUL2, &
 & ZPUL, ZQCORP  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "fcttrm.ycst.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACDIFSP',0,ZHOOK_HANDLE)
ASSOCIATE(TMERGL=>YDML_PHY_MF%YRPHY1%TMERGL, &
 & TSPHY=>YDML_PHY_MF%YRPHY2%TSPHY, &
 & LNEIGE=>YDML_PHY_MF%YRPHY%LNEIGE)
!-----------------------------------------------------------------------

!*
!     ------------------------------------------------------------------
!     I - AUXILIARY CONSTANTS.

ZCPVMD=YDCST%RCPV-YDCST%RCPD
ZCPVMW=YDCST%RCPV-YDCST%RCW
ZCPVMS=YDCST%RCPV-YDCST%RCS

!*
!     ------------------------------------------------------------------
!     II - COMPUTATION OF DERIVED PARAMETERS AND SECURITY CONSTANT (FOR
!     CSDT).

ZGDT=YDCST%RG*TSPHY
ZGDTI=1.0_JPRB/ZGDT
ZPUL=2.0_JPRB*YDCST%RPI/YDCST%RDAY

ZEPS1=1.E-08_JPRB

!*
!     ------------------------------------------------------------------
!     III - COMPUTATION OF SPECIFIC HEAT AT CONSTANT PRESSURE FOR AIR
!     Cp AND AT THE SURFACE Cps

! - TEMPORARY 2D.

! ZCP      : VALUE OF ATMOSPHERIC Cp

! - TEMPORARY 1D.

! ZCPS1    : VALUE OF Cp AT THE SURFACE AT t+dt

DO JLON=KIDIA,KFDIA
  ZCPS1(JLON)=YDCST%RCPD+ZCPVMD*PQS1(JLON)
ENDDO

DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    ZCP(JLON,JLEV)=YDCST%RCPD*(1.0_JPRB-PQ(JLON,JLEV))+YDCST%RCPV*PQ(JLON,JLEV)
  ENDDO
ENDDO

!*
!     ------------------------------------------------------------------
!     IV - PRELIMINARY COMPUTATIONS OF DRY STATIC ENERGY, SPECIFIC
!     HUMIDITY CORRECTED OF ITS NEGATIVE VALUES THROUGH UPWARDS PUMPING
!     (EXCEPT IN THE CASE OF THE LOWEST LEVEL ASSUMED TO POSSES A
!     FICTITIOUS SOURCE AT THE SURFACE) AND OF LAYERS' CHARACTERISTICS
!     IN THE UNITS OF THE NORMALIZED EXCHANGE COEFFICIENTS (KG/M**2/S)
!     FOR "POID" AND ITS INVERSE FOR "IPOI".

! - TEMPORARY(S) 2D (0:KLEV) .

! ZXURO     : PKUROV MODIFIED TO ATTENUATE FIBRILLATIONS.
! ZFCQNG    : PSEUDO-FLUX OF WATER TO CORRECT FOR Q<0.

! - TEMPORARY(S) 2D (1:KLEV) .

! ZPOID     : DP/(RG*DT) FOR A GIVEN LAYER AND A GIVEN TIME STEP.
! ZIPOI     : INVERSE OF ZPOID.
! ZDSE      : LOCAL ARRAY FOR DRY STATIC ENERGIES.
! ZQCOR     : LOCAL ARRAY FOR Q AFTER CORRECTION OF NEGATIVE VALUES.

! - TEMPORARY 1D.

! ZXDRO     : PCDROV MODIFIED TO ATTENUATE FIBRILLATIONS.
! ZDSES1    : LOCAL ARRAY FOR DRY STATIC ENERGIES AT THE SURFACE 
!             AT t+dt.

DO JLON=KIDIA,KFDIA
  ZFCQNG(JLON,KTDIA-1)=0.0_JPRB
  ZDSES1(JLON)=ZCPS1(JLON)*PTS1(JLON)+PAPHI(JLON,KLEV)
ENDDO

DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    ZPOID(JLON,JLEV)=PDELP(JLON,JLEV)*ZGDTI
    ZIPOI(JLON,JLEV)=PRDELP(JLON,JLEV)*ZGDT
    ZDSE(JLON,JLEV)=ZCP(JLON,JLEV)*PT(JLON,JLEV)+PAPHIF(JLON,JLEV)
    ZQCORP=PQ(JLON,JLEV)+ZFCQNG(JLON,JLEV-1)*ZIPOI(JLON,JLEV)
    ZFCQNG(JLON,JLEV)=-MAX(0.0_JPRB,-ZQCORP)*ZPOID(JLON,JLEV)
    ZQCOR(JLON,JLEV)=ZQCORP+MAX(0.0_JPRB,-ZQCORP)
  ENDDO
ENDDO

!     POTENTIAL "ANTI-FIBRILLATION" AMPLIFICATION OF THE EXCHANGE
!     COEFFICIENTS.

DO JLON=KIDIA,KFDIA
  ZXDRO(JLON)=PCDROV(JLON)*PXDROV(JLON)
ENDDO

DO JLEV=KTDIA,KLEV-1
  DO JLON=KIDIA,KFDIA
    ZXURO(JLON,JLEV)=PKUROV(JLON,JLEV)*PXUROV(JLON,JLEV)
  ENDDO
ENDDO

!*

!     ------------------------------------------------------------------
!     V - COMPUTATIONS FOR MOMENTUM / ELIMINATION PART. IN ALL THE
!     FOLLOWING ONE WILL USE ZELIM1 AND ZELIM2 (1 FOR THE COUPLE U,V
!     AND 2 FOR THE COUPLE S,Q) TO NAME THE SUBDIAGONAL ELEMENT BY
!     WHICH ONE SHOULD MULTIPLY THE PRECEEDING LINE TO ELIMINATE IT BY
!     SUBSTRACTION AND ZMUL1 AND ZMUL2 (1 FOR THE COUPLE U,V AND
!     2 FOR THE COUPLE S,Q) TO NAME THE MULTIPLYING FACTOR TO PUT THE
!     MAIN DIAGONAL TO "1". ZSUB1 AND ZSUB2 WILL CONTAIN THE MODIFIED
!     SUPERDIAGONAL COEFFICIENTS THAT WILL BE USED IN THE SUBSTITUTION
!     (1 FOR THE COUPLE U,V AND 2 FOR Q,S) AND ZN1, ZN2, ZN3 AND ZN4
!     WILL CONTAIN THE RIGHT HAND SIDES DUE TO BECOME THE FINAL VALUES
!     (1 FOR U, 2 FOR V, 3 FOR S AND 4 FOR Q).

! - TEMPORARY(S) 2D (1:KLEV) .

! ZSUB1     : ARRAY OF COEFFICIENTS FOR THE U,V SUBSTITUTION.
! ZSUB2     : ARRAY OF COEFFICIENTS FOR THE S,Q SUBSTITUTION.
! ZN1       : ARRAY OF RIGH HAND SIDES FOR U.
! ZN2       : ARRAY OF RIGHT HAND SIDES FOR V.
! ZN3       : ARRAY OF RIGH HAND SIDES FOR S (S=CP*T+PHI).
! ZN4       : ARRAY OF RIGHT HAND SIDES FOR Q.

!     ELIMINATION AT THE TOP.

DO JLON=KIDIA,KFDIA
  ZMUL1=1.0_JPRB/(1.0_JPRB+ZXURO(JLON,KTDIA)*ZIPOI(JLON,KTDIA))
  ZSUB1(JLON,KTDIA)=ZMUL1*ZXURO(JLON,KTDIA)*ZIPOI(JLON,KTDIA)
  ZN1(JLON,KTDIA)=ZMUL1*PU(JLON,KTDIA)
  ZN2(JLON,KTDIA)=ZMUL1*PV(JLON,KTDIA)
ENDDO

!     ELIMINATION FOR A STANDART LEVEL.

DO JLEV=KTDIA+1,KLEV-1
  DO JLON=KIDIA,KFDIA
    ZELIM1=ZXURO(JLON,JLEV-1)*ZIPOI(JLON,JLEV)
    ZMUL1=1.0_JPRB/(1.0_JPRB+ZELIM1*(1.0_JPRB-ZSUB1(JLON,JLEV-1))+ZXURO(JLON,JLEV)&
     & *ZIPOI(JLON,JLEV))  
    ZSUB1(JLON,JLEV)=ZMUL1*ZXURO(JLON,JLEV)*ZIPOI(JLON,JLEV)
    ZN1(JLON,JLEV)=ZMUL1*(PU(JLON,JLEV)+ZELIM1*ZN1(JLON,JLEV-1))
    ZN2(JLON,JLEV)=ZMUL1*(PV(JLON,JLEV)+ZELIM1*ZN2(JLON,JLEV-1))
  ENDDO
ENDDO

!*
!     ------------------------------------------------------------------
!     VI - COMPUTATION FOR MOMENTUM: SURFACE CONDITION, "CHARNOCK-TYPE"
!     CALCULATION THEN SUBSTITUTION WITH FLUXES' COMPUTATION AT ALL
!     LEVELS.

!     SURFACE CALCULATIONS.

DO JLON=KIDIA,KFDIA
  ZELIM1=ZXURO(JLON,KLEV-1)*ZIPOI(JLON,KLEV)
  ZMUL1=1.0_JPRB/(1.0_JPRB+ZELIM1*(1.0_JPRB-ZSUB1(JLON,KLEV-1))+ZXDRO(JLON)&
   & *ZIPOI(JLON,KLEV))  
  ZN1(JLON,KLEV)=ZMUL1*(PU(JLON,KLEV)+ZELIM1*ZN1(JLON,KLEV-1))
  ZN2(JLON,KLEV)=ZMUL1*(PV(JLON,KLEV)+ZELIM1*ZN2(JLON,KLEV-1))
  PSTRTU(JLON,KLEV)=PCDROV(JLON)*ZN1(JLON,KLEV)
  PSTRTV(JLON,KLEV)=PCDROV(JLON)*ZN2(JLON,KLEV)
ENDDO

!     BACK-SUBSTITUTION FOR A STANDART LAYER AND AT THE TOP.

DO JLEV=KLEV-1,KTDIA,-1
  DO JLON=KIDIA,KFDIA
    ZN1(JLON,JLEV)=ZN1(JLON,JLEV)+ZSUB1(JLON,JLEV)*ZN1(JLON,JLEV+1)
    ZN2(JLON,JLEV)=ZN2(JLON,JLEV)+ZSUB1(JLON,JLEV)*ZN2(JLON,JLEV+1)
    PSTRTU(JLON,JLEV)=PKUROV(JLON,JLEV)*(ZN1(JLON,JLEV)-ZN1(JLON,JLEV+1))
    PSTRTV(JLON,JLEV)=PKUROV(JLON,JLEV)*(ZN2(JLON,JLEV)-ZN2(JLON,JLEV+1))
  ENDDO
ENDDO

!     UPPER BOUNDARY CONDITION.

DO JLON=KIDIA,KFDIA
  PSTRTU(JLON,KTDIA-1)=0.0_JPRB
  PSTRTV(JLON,KTDIA-1)=0.0_JPRB
ENDDO

!*
!     ------------------------------------------------------------------
!     VII - PRELIMINARY COMPUTATIONS AT THE SURFACE FOR THE COUPLING OF
!     THE IMPLICIT SOLUTIONS IN S AND Q THROUGH THE TS VALUE (ONE WILL
!     HAVE TWO VALUES, A PRIORI IDENTICAL, OF THE NORMALIZED CH
!     COEFFICIENT TO ALLOW SENSITIVITY STUDIES FOR SENSIBLE/LATENT
!     HEATS).

DO JLON=KIDIA,KFDIA

!     COMPUTATIONS CORRESPONDING TO EVAPORATION/SUBLIMATION.

! - TEMPORARY(S) 1D .

! ZPN       : PROPORTION OF THE SURFACE EVAPORATING IN ICE PHASE.

  IF (LNEIGE) THEN
    ZPN(JLON)=PLSM(JLON)*PNEIJ(JLON)+(1.0_JPRB-PLSM(JLON))&
     & *MAX(0.0_JPRB,SIGN(1.0_JPRB,TMERGL-PTS(JLON)))  
  ELSE
    ZPN(JLON)=0.0_JPRB
  ENDIF

!     COMPUTATION OF THE WEIGHTED LATENT HEAT AND OF ITS TWO COMPONENTS
!     THROUGH CALLS TO THE DEFINED FUNCTION.

! - TEMPORARY(S) 1D .

! ZLHW      : CONTRIBUTION OF THE LIQUID PHASE TO THE LATENT HEAT.
! ZLHS      : CONTRIBUTION OF THE ICE PHASE TO THE LATENT HEAT.

  ZLHW(JLON)=(1.0_JPRB-ZPN(JLON))* FOLH (PTS(JLON),0.0_JPRB)
  ZLHS(JLON)=ZPN(JLON)* FOLH (PTS(JLON),1.0_JPRB)
  PLHS(JLON)=ZLHW(JLON)+ZLHS(JLON)

!     DIFFERENTIATION BETWEEN THE SURFACE EXCHANGE COEFFICIENTS FOR
!     LATENT/SENSIBLE HEAT.

!     POTENTIAL "ANTI-FIBRILLATION" AMPLIFICATION OF THE EXCHANGE
!     COEFFICIENTS.

! - TEMPORARY(S) 1D .

! ZCHRS     : EXCHANGE COEFFICIENT PCHROV FOR SENSIBLE HEAT.
! ZXSRO     : ZCHRS MODIFIED TO ATTENUATE FIBRILLATIONS.

  ZCHRS(JLON)=PCHROV(JLON)
  ZXSRO(JLON)=ZCHRS(JLON)*PXHROV(JLON)

!     CONTRIBUTION TO THE COEFFICIENTS OF THE TRIDIAGONAL MATRIX.

! - TEMPORARY(S) 1D .

! ZDLS      : TEMPERATURE DERIVATIVE OF THE TOTAL FLUX ; SENSIBLE PART.
! ZDLLW     : AS ZDLS ; MINUS LATENT PART, LIQUID PHASE.
! ZDLLS     : AS ZDLS ; MINUS LATENT PART, ICE PHASE.

  ZDLS(JLON)=(PQS(JLON)-PQ(JLON,KLEV))*ZCHRS(JLON)*ZCPVMD
  ZDLLW(JLON)=(PQS(JLON)-PQ(JLON,KLEV))*ZCHRS(JLON)*ZCPVMW*(1.0_JPRB-ZPN(JLON))
  ZDLLS(JLON)=(PQS(JLON)-PQ(JLON,KLEV))*ZCHRS(JLON)*ZCPVMS*ZPN(JLON)
  ZDLDQ=ZDLLW(JLON)+ZDLLS(JLON)-ZDLS(JLON)*PTS(JLON)
ENDDO

!*
!     ------------------------------------------------------------------
!     VIII - COMPUTATION FOR TEMPERATURE AND SPECIFIC HUMIDITY :
!     ELIMINATION PART. LIKE AT THE SURFACE ONE DIFFERENTIATES THE
!     (A PRIORI EQUAL) NORMALIZED COEFFICIENTS FOR S AND Q.

! - TEMPORARY(S) 2D (0:KLEV) .

! ZKSRV     : VALUE (DUPLICATED FROM PKTROV) FOR THE S DIFFUSION.
! ZXTRO     : ZKSRV MODIFIED TO ATTENUATE FIBRILLATIONS.

!     THE VERTICAL EXCHANGE COEFFICIENT FOR S AND Q.

!     POTENTIAL "ANTI-FIBRILLATION" AMPLIFICATION OF THE EXCHANGE
!     COEFFICIENT.

DO JLEV=KTDIA,KLEV-1
  DO JLON=KIDIA,KFDIA
    ZKSRV(JLON,JLEV)=PKTROV(JLON,JLEV)
    ZXTRO(JLON,JLEV)=ZKSRV(JLON,JLEV)*PXTROV(JLON,JLEV)
  ENDDO
ENDDO

!     TOP ELIMINATION FOR THE DRY STATIC ENERGY, S AND FOR THE SPECIFIC
!     HUMIDITY CORRECTED FROM NEGATIVE VALUES, Q.

DO JLON=KIDIA,KFDIA
  ZMUL2=1.0_JPRB/(1.0_JPRB+ZXTRO(JLON,KTDIA)*ZIPOI(JLON,KTDIA))
  ZSUB2(JLON,KTDIA)=ZMUL2*ZXTRO(JLON,KTDIA)*ZIPOI(JLON,KTDIA)
  ZN3(JLON,KTDIA)=ZMUL2*ZDSE(JLON,KTDIA)
  ZN4(JLON,KTDIA)=ZMUL2*ZQCOR(JLON,KTDIA)
ENDDO

!     ELIMINATION FOR A STANDART LAYER FOR S AND Q.

DO JLEV=KTDIA+1,KLEV-1
  DO JLON=KIDIA,KFDIA
    ZELIM2=ZXTRO(JLON,JLEV-1)*ZIPOI(JLON,JLEV)
    ZMUL2=1.0_JPRB/(1.0_JPRB+ZELIM2*(1.0_JPRB-ZSUB2(JLON,JLEV-1))+ZXTRO(JLON,JLEV)&
     & *ZIPOI(JLON,JLEV))  
    ZSUB2(JLON,JLEV)=ZMUL2*ZXTRO(JLON,JLEV)*ZIPOI(JLON,JLEV)
    ZN3(JLON,JLEV)=ZMUL2*(ZDSE(JLON,JLEV)+ZELIM2*ZN3(JLON,JLEV-1))
    ZN4(JLON,JLEV)=ZMUL2*(ZQCOR(JLON,JLEV)+ZELIM2*ZN4(JLON,JLEV-1))
  ENDDO
ENDDO

!     SURFACE ELIMINATION CALCULATION FOR S AND Q.

DO JLON=KIDIA,KFDIA
  ZELIM2=ZXTRO(JLON,KLEV-1)*ZIPOI(JLON,KLEV)
  ZMUL2=1.0_JPRB/(1.0_JPRB+ZELIM2*(1.0_JPRB-ZSUB2(JLON,KLEV-1))&
   & +ZXSRO(JLON)*ZIPOI(JLON,KLEV))  
  ZN3(JLON,KLEV)=ZMUL2*(ZDSE(JLON,KLEV)&
   & +ZXSRO(JLON)*ZIPOI(JLON,KLEV)*ZDSES1(JLON)&
   & +ZELIM2*ZN3(JLON,KLEV-1))  
  ZN4(JLON,KLEV)=ZMUL2*(ZQCOR(JLON,KLEV)&
   & +ZXSRO(JLON)*ZIPOI(JLON,KLEV)*PQS1(JLON)&
   & +ZELIM2*ZN4(JLON,KLEV-1))  
ENDDO

!*
!     ------------------------------------------------------------------
!     IX - COMPUTATIONS FOR TEMPERATURE AND SPECIFIC HUMIDITY :
!     BACK-SUBSTITUTION PART AND CALCULATION OF FLUXES AT ALL LEVELS AS
!     WELL AS AT THE SURFACE AND IN THE SOIL.

!     SURFACE AND SOIL FLUXES.

DO JLON=KIDIA,KFDIA
  ZFQ=ZCHRS(JLON)*(ZN4(JLON,KLEV)-PQS1(JLON))
  PFEVL(JLON)=ZFQ*(1.0_JPRB-ZPN(JLON))
  PFEVN(JLON)=ZFQ*ZPN(JLON)
  PFCLL(JLON)=ZLHW(JLON)*ZFQ-ZDLLW(JLON)*(PNTS(JLON)-PTS(JLON))
  PFCLN(JLON)=ZLHS(JLON)*ZFQ-ZDLLS(JLON)*(PNTS(JLON)-PTS(JLON))
  ZFENT=ZCHRS(JLON)*(ZN3(JLON,KLEV)-ZDSES1(JLON))
  PFCS(JLON)=ZFENT-ZCPVMD*PTS(JLON)*ZFQ+ZDLS(JLON)*(PNTS(JLON)-PTS(JLON))
  PDIFTQ(JLON,KLEV)=ZFQ+ZFCQNG(JLON,KLEV)
  PDIFTS(JLON,KLEV)=ZFENT

ENDDO

!     BACK-SUBSTITUTION FOR A STANDART LEVEL.
!     ATMOSPHERIC FLUXES.

DO JLEV=KLEV-1,KTDIA,-1
  DO JLON=KIDIA,KFDIA
    ZN3(JLON,JLEV)=ZN3(JLON,JLEV)+ZSUB2(JLON,JLEV)*ZN3(JLON,JLEV+1)
    ZN4(JLON,JLEV)=ZN4(JLON,JLEV)+ZSUB2(JLON,JLEV)*ZN4(JLON,JLEV+1)

    PDIFTQ(JLON,JLEV)=ZKSRV(JLON,JLEV)*(ZN4(JLON,JLEV)&
     & -ZN4(JLON,JLEV+1))+ZFCQNG(JLON,JLEV)  
    PDIFTS(JLON,JLEV)=ZKSRV(JLON,JLEV)*(ZN3(JLON,JLEV)-ZN3(JLON,JLEV+1))
  ENDDO
ENDDO

!     UPPER BOUNDARY CONDITION.

DO JLON=KIDIA,KFDIA
  PDIFTQ(JLON,KTDIA-1)=0.0_JPRB
  PDIFTS(JLON,KTDIA-1)=0.0_JPRB
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACDIFSP',1,ZHOOK_HANDLE)
END SUBROUTINE ACDIFSP
