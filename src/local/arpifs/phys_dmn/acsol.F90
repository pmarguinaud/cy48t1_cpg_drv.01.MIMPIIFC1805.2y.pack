!OPTIONS XOPT(NOEVAL)
SUBROUTINE ACSOL ( YDPHY,YDPHY1,KIDIA,KFDIA,KLON,&
 !-----------------------------------------------------------------------
 ! - INPUT  1D
 & PARG,PD2,PGZ0F,PGZ0HF,PGZ0RLF,PLSM,PIVEG,PLAI,PALBNS,&
 & PRHONS,PSAB,PSNS,PTS,PVEG0,PWP,PWPI,PWS,PWSI,&
 ! - INPUT  LOGIQUE
 & LDHMT,&
 ! - OUTPUT 1D .
 & PC1,PC2,PC3,PCG,PCN,PCT,&
 & PNEIJG,PNEIJV,&
 & PWFC,PWPMX,PWSEQ,PWSMX,PWWILT)

!**** *ACSOL  * - DETERMINATION DES CARACTERISTIQUES DE SURFACE.

!     Sujet.
!     ------

!     - ROUTINE DE CALCUL ACTIF .
!       DETERMINATION DES CARACTERISTIQUES DE SURFACE CONCERNANT LES
!       COEFFICIENTS THERMO-HYDRIQUES, LES TENEURS EN EAU MAXIMALES
!       DES RESERVOIRS PROFOND ET SUPERFICIEL AINSI QUE LA TENEUR EN
!       EAU A L'EQUILIBRE DU RESERVOIR SUPERFICIEL.
!       CAS DE L EAU: ON POSE C1=C2=1.

!**   Interface.
!     ----------
!        *CALL* *ACSOL*

!-----------------------------------------------------------------------
! WARNING: THE ENGLISH VERSION OF VARIABLES' NAMES IS TO BE READ IN THE
!          "APLPAR" CODE.
!-----------------------------------------------------------------------

! -   ARGUMENTS D'ENTREE.
!     -------------------

! - NOM DES PARAMETRES DE DIMENSIONNEMENT DE LA PHYSIQUE.

! KIDIA      : INDICE DE DEPART DES BOUCLES VECTORISEES SUR L'HORIZONT..
! KFDIA      : INDICE DE FIN DES BOUCLES VECTORISEES SUR L'HORIZONTALE.
! KLON       : DIMENSION HORIZONTALE DES TABLEAUX.

! - NOM DES CLES LOGIQUES

! LDHMT      : ETAPE DANS LE CALCUL DES CHAMPS AUX HAUTEURS METEO

! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
!   CATEGORIE).

! - 1D (PROGNOSTIQUE) .

! PRHONS     : DENSITE DE LA NEIGE.
! PSNS       : MASSE DE NEIGE PAR UNITE DE SURFACE.
! PTS        : TEMPERATURE DE SURFACE.
! PWP        : CONTENU EN EAU DU RESERVOIR PROFOND.
! PWPI       : CONTENU EN EAU GELEE DU RESERVOIR PROFOND.
! PWS        : CONTENU EN EAU DU RESERVOIR DE SURFACE.
! PWSI       : CONTENU EN EAU GELEE DU RESERVOIR DE SURFACE.

! - 1D (GEOGRAPHIQUE) .

! PARG       : POURCENTAGE D'ARGILE DANS LA MAILLE.
! PD2        : EPAISSEUR DU RESERVOIR PROFOND.
! PGZ0F      : GRAVITE * LONGUEUR DE RUGOSITE.
! PGZ0HF     : GRAVITE * LONGUEUR DE RUGOSITE THERMIQUE.
! PGZ0RLF    : GRAVITE * LONGUEUR DE RUGOSITE DUE AU RELIEF.
! PLSM       : INDICE TERRE/MER.
! PIVEG      : TYPE DE SURFACE (MER, BANQUISE OU GLACIER, TERRE).
! PLAI       : INDICE FOLIAIRE.
! PALBNS     : ALBEDO DE LA NEIGE.
! PSAB       : POURCENTAGE DE SABLE DANS LA MAILLE.
! PVEG0      : PROPORTION DE SOL COUVERT PAR LA VEGETATION.

!-----------------------------------------------------------------------

! -   ARGUMENTS DE SORTIE.
!     --------------------

! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
!   CATEGORIE).

! - 1D (DIAGNOSTIQUE) .

! PC1        : COEFFICIENT HYDRIQUE REPRESENTANT L'INTENSITE AVEC LAQUEL
!              LES FLUX DE SURFACE PARTICIPENT A L'EVOLUTION DE WS.
! PC2        : COEFFICIENT HYDRIQUE TRADUISANT LA RAPIDITE DES TRANSFERT
!              D'EAU ENTRE LES RESERVOIRS PROFOND ET DE SURFACE.
! PC3        : COEFFICIENT UTILE POUR LE CALCUL DU DRAINAGE
! PCG        : COEFFICIENT THERMIQUE DU SOL NU.
! PCN        : COEFFICIENT THERMIQUE DE LA NEIGE.
! PCT        : COEFFICIENT THERMIQUE DU MILIEU SOL-VEGETATION.
! PNEIJG     : FRACTION DE NEIGE RECOUVRANT LE SOL.
! PNEIJV     : FRACTION DE NEIGE RECOUVRANT LA VEGETATION.
! PWFC       : TENEUR EN EAU A LA CAPACITE AUX CHAMPS.
! PWPMX      : TENEUR EN EAU MAXIMALE DU RESERVOIR PROFOND.
! PWSEQ      : TENEUR EN EAU A L'EQUILIBRE POUR LE RESERVOIR DE SURFACE.
!              (EQUILIBRE ENTRE FORCES DE GRAVITE ET CAPILLARITE).
! PWSMX      : TENEUR EN EAU MAXIMALE DU RESERVOIR DE SURFACE.
! PWWILT     : TENEUR EN EAU CORRESPONDANT AU POINT DE FLETRISSEMENT.

!-----------------------------------------------------------------------

! -   ARGUMENTS IMPLICITES.
!     ---------------------

!-----------------------------------------------------------------------

!     Externes :  ACSOLW
!     ---------

!     Methode.
!     --------

!     Auteur.
!     -------
!      91-12, J. Noilhan.

!     Modifications.
!     --------------
!      R. El Khatib : 01-08-07 Pruning options
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!      2011-06: M. Jerczynski - some cleaning to meet norms
!      2012-12: E. Bazile and E. Brun - change of soil inertia for ice Cap.
!      2018-09, J. Masek: Calculation of snow fractions for LVGSN=T moved
!        here from APLPAR. Coding of ALARO-1 fixes for LZ0HSREL=T.
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMPHY1  , ONLY : TPHY1
USE YOMPHY   , ONLY : TPHY
USE YOMCST   , ONLY : RG       ,RTT
USE YOMCLI   , ONLY : STHER

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TPHY)        ,INTENT(IN)    :: YDPHY
TYPE(TPHY1)       ,INTENT(IN)    :: YDPHY1
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PARG(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PD2(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGZ0F(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGZ0HF(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGZ0RLF(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIVEG(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAI(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBNS(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHONS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSAB(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSNS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVEG0(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWP(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWPI(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSI(KLON) 
LOGICAL           ,INTENT(IN)    :: LDHMT 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PC1(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PC2(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PC3(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCG(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCN(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCT(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNEIJG(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNEIJV(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWFC(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWPMX(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWSEQ(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWSMX(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWWILT(KLON) 
REAL(KIND=JPRB) :: ZWSAT(KLON)

INTEGER(KIND=JPIM) :: JLON

REAL(KIND=JPRB) :: Z1, Z2, ZA, ZARG, ZB, ZBWA, ZBWB, ZC1MAX,&
 & ZC1SAT, ZC1X2, ZC1Y2, ZC2REF, ZC3REF, ZCGSAT, &
 & ZCRIN, ZD1, ZD2, ZDELTA, ZEPS, ZEPS1, ZEPW, &
 & ZLYMY1, ZNEIJC, ZNEIJV, ZP, ZSAB, ZTS, ZTSVAR, &
 & ZWP, ZWPI, ZWPSWS, ZWS, ZWSSWS, ZZ0V, ZZA, &
 & ZZB, ZZC  
REAL(KIND=JPRB) :: ZCOEF,ZCFN,ZCK,Z0CR,ZUZ0CN
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "acsolw.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACSOL',0,ZHOOK_HANDLE)
ASSOCIATE(NTVGLA=>YDPHY1%NTVGLA, GC2REF=>YDPHY1%GC2REF, G2CGSAT=>YDPHY1%G2CGSAT, &
 & NTVMER=>YDPHY1%NTVMER, G1C1SAT=>YDPHY1%G1C1SAT, WCRIN=>YDPHY1%WCRIN, &
 & XCRINR=>YDPHY1%XCRINR, EC2REF=>YDPHY1%EC2REF, RCTVEG=>YDPHY1%RCTVEG, &
 & RD1=>YDPHY1%RD1, GC1S3=>YDPHY1%GC1S3, GC1S2=>YDPHY1%GC1S2, &
 & GC1S1=>YDPHY1%GC1S1, G1B=>YDPHY1%G1B, GSNC2=>YDPHY1%GSNC2, &
 & GC1S4=>YDPHY1%GC1S4, GC1=>YDPHY1%GC1, GSNC1=>YDPHY1%GSNC1, GC3=>YDPHY1%GC3, &
 & GC2=>YDPHY1%GC2, GC1Y1=>YDPHY1%GC1Y1, G3CGSAT=>YDPHY1%G3CGSAT, &
 & RCGMAX=>YDPHY1%RCGMAX, WCRINC=>YDPHY1%WCRINC, G2P=>YDPHY1%G2P, &
 & WCRING=>YDPHY1%WCRING, G1P=>YDPHY1%G1P, EA=>YDPHY1%EA, LIMW=>YDPHY1%LIMW, &
 & GA=>YDPHY1%GA, LIMC=>YDPHY1%LIMC, GTSVAP=>YDPHY1%GTSVAP, &
 & G1CGSAT=>YDPHY1%G1CGSAT, RC1MAX=>YDPHY1%RC1MAX, GC32=>YDPHY1%GC32, &
 & GC31=>YDPHY1%GC31, XCRINV=>YDPHY1%XCRINV, G2C1SAT=>YDPHY1%G2C1SAT, &
 & GCONV=>YDPHY1%GCONV, G2B=>YDPHY1%G2B, LC1VAP=>YDPHY1%LC1VAP, &
 & RCTGLA=>YDPHY1%RCTGLA, ALB1=>YDPHY1%ALB1, ALB2=>YDPHY1%ALB2, &
 & ALRCN1=>YDPHY1%ALRCN1, ALRCN2=>YDPHY1%ALRCN2, &
 & RLAIMX=>YDPHY1%RLAIMX, RLAI=>YDPHY1%RLAI, &
 & LSNV=>YDPHY%LSNV, LGLACIERS=>YDPHY%LGLACIERS, LFGEL=>YDPHY%LFGEL, &
 & LVGSN=>YDPHY%LVGSN, LZ0HSREL=>YDPHY%LZ0HSREL, LCOEFKSURF=>YDPHY%LCOEFKSURF)
!-----------------------------------------------------------------------
ZEPS=1.E-1_JPRB
ZEPS1=1.E-5_JPRB
ZEPW=1.E-2_JPRB
ZTSVAR=MAX(0.0_JPRB,SIGN(1.0_JPRB,GTSVAP))
ZD1=RD1*GCONV

! GOING TO GEOPOTENTIAL FOR CONSTANTS RELATED TO ROUGHNESS LENGTH OF SNOW
Z0CR=RG*ALRCN1
ZUZ0CN=1.0_JPRB/(RG*ALRCN2)

!     ------------------------------------------------------------------
!     I - CALCULS DES TENEURS EN EAU DU SOL CARACTERISTIQUES
!         --------------------------------------------------
!         COMPUTING CHARACTERISTIC SOIL MOISTURES
!         ---------------------------------------

CALL ACSOLW(YDPHY1,KIDIA,KFDIA,KLON,&
 & PARG,PD2,PLSM,PIVEG,PSAB,&
 & LDHMT,&
 & PWFC,PWPMX,ZWSAT,PWSMX,PWWILT)

IF ( .NOT. LDHMT ) THEN
  
!     ------------------------------------------------------------------
!     II - CALCUL DES AUTRES CHAMPS
!          ------------------------
!          COMPUTING OTHER FIELDS
!          ----------------------
  
  DO JLON=KIDIA,KFDIA
  
! ** INTERMEDIAIRES NE DEPENDANT QUE DE LA TEXTURE **
  
    ZSAB = MAX(ZEPS,PSAB(JLON))
    ZARG = MAX(ZEPS,PARG(JLON))
  
    ZA = GA*(ZARG**EA)
    ZB = G1B*ZARG+G2B
    ZP = G1P*ZARG+G2P
    ZC1SAT = (G1C1SAT*ZARG+G2C1SAT)*GC2
    ZCGSAT = G1CGSAT*ZSAB + G2CGSAT*ZARG + G3CGSAT
    ZC2REF = GC2REF*(ZARG**EC2REF)
    ZC3REF = GC31*(ZARG**GC32)
  
! ** WSEQ **
  
    ZD2=PD2(JLON)*GCONV
  
    ZWP = PWP(JLON)/ZD2
    IF (LFGEL) THEN
      ZWPSWS = MAX(ZEPW,ZWP/MAX(ZEPW,ZWSAT(JLON)-PWPI(JLON)/ZD2))
    ELSE
      ZWPSWS = MAX(ZEPW,ZWP/ZWSAT(JLON))
    ENDIF
    IF ( NINT(PIVEG(JLON)) == NTVGLA ) ZWPSWS = 1.0_JPRB
    PWSEQ(JLON) = (ZWPSWS-ZA*(ZWPSWS**ZP)*(1.0_JPRB-(ZWPSWS**(GC3*ZP))))&
     & * PWSMX(JLON)  
  
! ** CAS DE LA MER - SEA **
  
    IF ( NINT(PIVEG(JLON)) == NTVMER ) THEN
      PSNS(JLON) = 0.0_JPRB
      PNEIJG(JLON) = 0.0_JPRB
      PNEIJV(JLON) = 0.0_JPRB
      PCN(JLON) = RCTGLA
      PCG(JLON) = ZCGSAT
      PCT(JLON) = PCG(JLON)
      PC1(JLON) = 100._JPRB*RD1
      PC2(JLON) = 1.0_JPRB
      PC3(JLON) = 1.0_JPRB
  
    ELSE
  
! ** EXTENSION DE LA COUVERTURE NEIGEUSE - SNOW COVER **
  
      IF (LSNV) THEN
        PSNS(JLON) = MAX(0.0_JPRB,PSNS(JLON))
        PNEIJG(JLON) = PSNS(JLON)&
         & /(PSNS(JLON)+WCRIN*(1.0_JPRB+PGZ0RLF(JLON)*XCRINR))  
        ZZ0V = SQRT(MAX(0.0_JPRB,PGZ0F(JLON)**2-PGZ0RLF(JLON)**2))/RG
        ZCRIN = MAX(WCRING,PRHONS(JLON)*XCRINV*ZZ0V)
        PNEIJV(JLON) = PNEIJG(JLON)*(PSNS(JLON)+WCRING)/(PSNS(JLON)+ZCRIN)
      ELSEIF (LVGSN) THEN
        IF (LZ0HSREL.AND.LCOEFKSURF) THEN
          PNEIJG(JLON)=PLSM(JLON)*PSNS(JLON)/(PSNS(JLON)+WCRIN*         &
           & (1.0_JPRB+ZUZ0CN*PGZ0HF(JLON)/STHER))
        ELSE
          PNEIJG(JLON)=PLSM(JLON)*PSNS(JLON)/(PSNS(JLON)+WCRIN)
        ENDIF
        ZCFN =(1.0_JPRB-MAX(0.0_JPRB,MIN(1.0_JPRB,PLAI(JLON)/RLAIMX)))* &
         & MAX(0.0_JPRB,SIGN(1.0_JPRB,PLAI(JLON)-RLAI))+                &
         & (1.0_JPRB-MAX(0.0_JPRB,SIGN(1.0_JPRB,PLAI(JLON)-RLAI)))
        ZCK  =MAX(0.0_JPRB,(ALB1-MAX(ALB2,PALBNS(JLON))))/(ALB1-ALB2)
        ZCOEF=ZCFN*ZCK+(1.0_JPRB-ZCK)
        PNEIJV(JLON)=PNEIJG(JLON)*ZCOEF
      ENDIF
  
      IF ( NINT(PIVEG(JLON)) /= NTVGLA ) THEN
  
! ** COEFFICIENTS THERMIQUES
!    DU SOL NU CG, DE LA NEIGE CN ET DU MILIEU SOL-NEIGE-VEGETATION CT **
! ** THERMICAL COEFFICIENTS
!    OF BARE SOIL CG, OF SNOW CN AND OF SOIL-SNOW-VEGETATION SYSTEM CT **
  
        ZBWB = GC1*ZB/LOG(GC2)
        IF (LIMW) ZWPSWS = MAX(ZEPW,PWWILT(JLON)/ZWSAT(JLON),ZWPSWS)
        PCG(JLON) = ZCGSAT/(ZWPSWS**ZBWB)
        IF (LFGEL) THEN
          ZWPI=PWPI(JLON)/MAX(PWP(JLON)+PWPI(JLON),ZEPS1)
          PCG(JLON)=1.0_JPRB/ ( (1-ZWPI)/PCG(JLON) + ZWPI/RCTGLA )
        ENDIF
        IF (LIMC)  PCG(JLON) = MIN(RCGMAX,PCG(JLON))
        IF (LSNV) THEN
          ZNEIJC = MIN(1.0_JPRB,PSNS(JLON)/WCRINC)
          ZNEIJV = MIN(ZNEIJC,PNEIJV(JLON))
          PCN(JLON) = 2.0_JPRB*SQRT( GSNC1/(PRHONS(JLON)**GSNC2) )
          PCT(JLON) = 1.0_JPRB/&
           & ( (1.0_JPRB-PVEG0(JLON))*(1.0_JPRB-ZNEIJC) /PCG(JLON)&
           & +(1.0_JPRB-PVEG0(JLON))*   ZNEIJC   /PCN(JLON)&
           & +  PVEG0(JLON)   *   ZNEIJV   /PCN(JLON)&
           & +  PVEG0(JLON)   *(1.0_JPRB-ZNEIJV) /&
           & RCTVEG(NINT(PIVEG(JLON))) )  
        ELSE
          PCT(JLON) = 1.0_JPRB/((1.0_JPRB-PVEG0(JLON)) /PCG(JLON)&
           & +   PVEG0(JLON)/RCTVEG(NINT(PIVEG(JLON))))  
        ENDIF
  
! ** COEFFICIENT HYDRAULIQUE C1 - HYDRAULIC COEFFICIENT C1 **
  
        ZBWA = GC1*ZB+1.0_JPRB
        IF (LFGEL) THEN
          ZWS=PWS(JLON)+PWSI(JLON)
        ELSE
          ZWS=PWS(JLON)
        ENDIF
        ZWSSWS = MAX(ZEPW,ZWS/PWSMX(JLON))
        IF (LIMW.AND.(.NOT.LC1VAP))&
         & ZWSSWS = MAX(ZEPW,PWWILT(JLON)/ZWSAT(JLON),ZWSSWS)  
        PC1(JLON) = ZC1SAT*RD1/(ZWSSWS**ZBWA)
        IF (LIMC.AND.(.NOT.LC1VAP))PC1(JLON) = MIN(RC1MAX*RD1,PC1(JLON))
  
        IF (LC1VAP.AND.((ZWS/ZD1) < PWWILT(JLON))) THEN
          ZTS=ZTSVAR * MAX(RTT,PTS(JLON)) + (1.0_JPRB-ZTSVAR) * (RTT-GTSVAP)
          ZC1MAX=(GC1S1*PWWILT(JLON)+GC1S2)*ZTS+ (GC1S3*PWWILT(JLON)+GC1S4)
          ZC1X2=PWWILT(JLON)*ZD1
          ZC1Y2=ZC1SAT/((PWWILT(JLON)/ZWSAT(JLON))**ZBWA)
          ZC1MAX=MAX(ZC1Y2,ZC1MAX)
          ZZA=LOG(ZC1Y2/GC1Y1)
          ZLYMY1=LOG(ZC1MAX/GC1Y1)
          ZZB = -2.0_JPRB * ZC1X2 * ZLYMY1
          ZZC = ZC1X2**2 * ZLYMY1
          ZDELTA = ZZB**2 - 4._JPRB * ZZA * ZZC
          ZDELTA=MAX(0.0_JPRB,ZDELTA)
          Z1=(-ZZB-SQRT(ZDELTA))/(2.0_JPRB*ZZA)
          Z2= Z1**2 / ZLYMY1
          PC1(JLON)= ZC1MAX * RD1 * EXP(-(ZWS-Z1)**2/Z2)
        ENDIF
  
! ** COEFFICIENTS HYDRAULIQUES C2, C3 - HYDRAULIC COEFFICIENTS C2, C3 **
        IF (LFGEL) THEN
          ZWP=ZWP+PWPI(JLON)/ZD2
        ENDIF
        IF (LIMW) ZWP=MAX(PWWILT(JLON),ZWP)
        PC2(JLON)=MAX(0.0_JPRB,ZC2REF*ZWP/(ZWSAT(JLON)-ZWP+ZEPW))
        PC3(JLON)=ZC3REF/PD2(JLON)
  
! ** CAS DE LA GLACE - ICE CAP **
  
      ELSE
        PCG(JLON) = RCTGLA
        IF (LGLACIERS) THEN
           ZNEIJC = MIN(1.0_JPRB,PSNS(JLON)/WCRINC)
           PCN(JLON) = 2.0_JPRB*SQRT( GSNC1/(PRHONS(JLON)**GSNC2) )
           PCT(JLON) = 1.0_JPRB/&
             & ( (1.0_JPRB-ZNEIJC) /PCG(JLON)&
             & + ZNEIJC   /PCN(JLON))
        ELSE
           PCN(JLON) = RCTGLA
           PCT(JLON) = RCTGLA
        ENDIF
        PC1(JLON) = 0.0_JPRB
        PC2(JLON) = 0.0_JPRB
        PC3(JLON) = 0.0_JPRB
      ENDIF
    ENDIF
  
  ENDDO
  
ENDIF 

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACSOL',1,ZHOOK_HANDLE)

END SUBROUTINE ACSOL
