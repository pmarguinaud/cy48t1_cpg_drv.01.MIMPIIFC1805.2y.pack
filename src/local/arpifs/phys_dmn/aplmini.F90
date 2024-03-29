!OPTIONS XOPT(NOEVAL)
SUBROUTINE APLMINI ( YDML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
 !-----------------------------------------------------------------------
 ! - INPUT 2D .
 & PAPRS,PAPRSF,PCP,PQIMP,PQLMP,PR,PTMP,&
 & PIPOI,PLHS,PLHV,PNEBM,PPOID,PTCORR,&
 ! - INPUT 1D .
 & PDECRD,&
 ! - OUTPUT 2D .
 & PMELNET,PMELGET)

!**** *APLMINI * - MINI MICROPHYSIC: CALCUL DU GEL/FONTE


!     Sujet.
!     ------
!     - ROUTINE DE CALCUL ACTIF .
!       CALCUL DU GEL/FONTE.
!     - COMPUTATION OF FREEZING/MELTING.

!**   Interface.
!     ----------
!        *CALL* *APLMINI*

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
! KTDIA      : INDICE DE DEPART DES BOUCLES VERTICALES (1 EN GENERAL).
! KLEV       : DIMENSION VERTICALE DES TABLEAUX "FULL LEVEL".

! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
!   CATEGORIE).

! - 2D (0:KLEV) .

! PAPRS      : PRESSION AUX DEMI-NIVEAUX.

! - 2D (1:KLEV) .

! PAPRSF     : PRESSION AUX NIVEAUX PLEIN.
! PCP        : CHALEUR MASSIQUE A PRESSION CONSTANTE.
! PQIMP      : RATIO OF SUSPENDED ICE
! PQLMP      : RATIO OF SUSPENDED LIQUID WATER
! PR         : CONSTANTE DES GAZ POUR L'AIR.
! PTMP       : TEMPERATURE.
! PIPOI      : FACTEUR DE CONVERSION DE D_P FLUX VERS D_T Q_X.
! PLHS       : CHALEUR LATENTE DE SUBLIMATION.
! PLHV       : CHALEUR LATENTE DE VAPORISATION.
! PNEBM      : NEBULOSITE POUR LA MICROPHYSIQUE.
! PPOID      : FACTEUR DE CONVERSION DE D_T Q_X VERS D_P FLUX.
! PTCORR     : CORRECTION DE TEMPERATURE POUR LE NUAGE CONVECTIF.

! - 1D .

! PDECRD     : DECORRELATION DEPTH FOR CLOUD OVERLAPS [Pa].

!-----------------------------------------------------------------------

! -   ARGUMENTS DE SORTIE.
!     --------------------

! - 2D (KLEV) .

! PMELNET    : NET MELTING RATE.
!-----------------------------------------------------------------------

! -   ARGUMENTS IMPLICITES.
!     ---------------------

! COMMON/YOMPHY /
! COMMON/YOMCST /
! COMMON/YOMPHY0/
! COMMON/YOMPHY2/

!-----------------------------------------------------------------------

!     Externes.
!     ---------
!        *ACACON*
!        *ACEVMEL*

!     Methode.
!     --------

!     Auteur.
!     -------
!        07-12, R. Brozkova.

!     Modifications.
!     --------------
!        08-10, Passing q vapour as argument to ACEVMEL - B. Catry
!        K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!        J. Van den Bergh: adapting for prognostic graupel 
!        J. Masek  (Apr 2016): Exponential-random cloud overlap with variable
!                              decorrelation depth.
!        L. Gerard (Sep 2019): New temperature dependency for autoconversion.
!-----------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST    ,ONLY : RTT


!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQIMP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQLMP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PR(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIPOI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLHS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLHV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNEBM(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPOID(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTCORR(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDECRD(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMELNET(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMELGET(KLON,KLEV)

!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZFPLSL(KLON,0:KLEV), ZFPLSN(KLON,0:KLEV),ZFPLSG(KLON,0:KLEV)

REAL(KIND=JPRB) :: ZIPB(KLON),ZIPBR(KLON),ZIPH(KLON),ZIPHR(KLON),ZDQ(KLON)&
 & ,ZTST(KLON),ZEXPA(KLON),ZEXPN(KLON),ZLSCP(KLON),ZRPCEFF(KLON)&
 & ,ZACORL(KLON),ZACONI(KLON),ZACOGI(KLON),ZACONL(KLON),ZACOGL(KLON)&
 & ,ZEVAR(KLON),ZEVAN(KLON),ZEVAG(KLON),ZFONTN(KLON),ZFONTG(KLON)&
 & ,ZZACORL(KLON),ZZACOGL(KLON),ZZACONL(KLON)&
 & ,ZZACONI(KLON),ZZACOGI(KLON)&
 & ,ZSTAL3(KLON),ZHPLSL(KLON)&
 & ,ZSTAN3(KLON),ZSTAG3(KLON),ZHPLSN(KLON),ZHPLSG(KLON)&
 & ,ZLPLSL(KLON),ZLPLSN(KLON),ZLPLSG(KLON)&
 & ,ZQRST(KLON),ZQNST(KLON),ZQGST(KLON),ZQLST(KLON),ZQIST(KLON)&
 & ,ZIPOI(KLON),ZPOID(KLON),ZPART(KLON),ZTEST(KLON),ZQVST(KLON)&
 & ,ZZRHO(KLON),ZZPRSF(KLON),ZLSM(KLON),ZNEIJ(KLON),ZTS(KLON)

REAL(KIND=JPRB) :: ZIPLSLO(KLON),ZIPLSLE(KLON),ZIPLSNO(KLON),ZIPLSNE(KLON)&
 & ,ZPRPLO(KLON),ZPRPLE(KLON),ZNEBLOC(KLON)&
 & ,ZOPLSLO(KLON),ZOPLSNO(KLON),ZOPLSLE(KLON),ZOPLSNE(KLON)&
 & ,ZIPLSGO(KLON),ZIPLSGE(KLON),ZOPLSGO(KLON),ZOPLSGE(KLON)

INTEGER(KIND=JPIM) :: JLEV, JLON

REAL(KIND=JPRB) :: ZEPS1, ZEXTMP, ZCLOV,&
 & ZMUL, ZZNEBLOC, ZW, ZOPRSL, ZOPRSN,ZOPRSG, ZBETABR1, ZBETABR2,&
 & ZZPRPLE, ZZTEST, ZIPLSL, ZIPLSN,ZIPLSG, ZNEBLI, ZNEBLOCO,&
 & ZUMNLOC
REAL(KIND=JPRB), PARAMETER :: ZALPHA=-0.1572_JPRB, ZBETA=-4.9632_JPRB

LOGICAL :: LLSIMP

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "acacon.intfb.h"
#include "acevmel.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('APLMINI',0,ZHOOK_HANDLE)
ASSOCIATE(REVAPN=>YDML_PHY_MF%YRPHY0%REVAPN, LAB12=>YDML_PHY_MF%YRPHY%LAB12, &
 & LGRAPRO=>YDML_PHY_MF%YRPHY%LGRAPRO)
!-----------------------------------------------------------------------

!*
!     ------------------------------------------------------------------
!     I - CONDITION A LA LIMITE SUPERIEURE, SUBSTITUTION P-Z  ET PARAMETRES.

!         UPPER BOUNDARY CONDITION P-Z EXCHANGE AND PARAMETERS.

! - TEMPORAIRE(S) 1D .

! ZIPH      : INVERSE DE LA PRESSION AU SOMMET DE LA COUCHE.
!           : INVERSE PRESSURE AT THE TOP OF THE LAYER.

LLSIMP=.TRUE.

ZEPS1=1.E-10_JPRB

ZEXTMP=0.0231_JPRB

DO JLON=KIDIA,KFDIA
  ZFPLSL(JLON,KTDIA-1)=0.0_JPRB
  ZFPLSN(JLON,KTDIA-1)=0.0_JPRB
  ZIPH(JLON)=1.0_JPRB/MAX(ZEPS1,PAPRS(JLON,KTDIA-1))
  ZIPHR(JLON)=ZIPH(JLON)

  IF (LAB12) THEN
    ZIPHR(JLON)=ZIPH(JLON)*SQRT(SQRT(ZIPH(JLON)))
  ENDIF
ENDDO

DO JLON=KIDIA,KFDIA
  ZNEBLOC(JLON)=PNEBM(JLON,KTDIA)
  ZPRPLO(JLON)=0.0_JPRB
  ZPRPLE(JLON)=0.0_JPRB
  ZIPLSLO(JLON)=0.0_JPRB
  ZIPLSLE(JLON)=0.0_JPRB
  ZIPLSNO(JLON)=0.0_JPRB
  ZIPLSNE(JLON)=0.0_JPRB
ENDDO

IF (LGRAPRO) THEN
  DO JLON=KIDIA,KFDIA
    ZFPLSG(JLON,KTDIA-1)=0.0_JPRB
    ZIPLSGO(JLON)=0.0_JPRB
    ZIPLSGE(JLON)=0.0_JPRB
  ENDDO
ENDIF
!*
!     ------------------------------------------------------------------
!     II - CALCULS PROPREMENTS DITS DANS UNE BOUCLE VERTICALE OU
!     L'INFORMATION SE TRANSMET DE COUCHE A COUCHE.

!          EFFECTIVE CALCULATIONS IN A VERTICAL LOOP WHERE THE
!     INFORMATION IS PASSED FROM LAYER TO LAYER.

DO JLEV=KTDIA,KLEV

!     CALCULS PRELIMINAIRES ET PREPARATION DE LA DEPENDANCE EN TEMPERATURE
!     POUR LES CONSTANTES NEIGEUSES (LES DEPENDANCES, CHOISIES D'APRES
!     LOPEZ (2002), SONT REGROUPEES).

!     PRELIMINARY COMPUTATIONS AND PREPARATIOB OF THE TEMPERATURE DEPENDENCY
!     FOR THE SNOW-RELATED CONSTANTS (ALL DEPENDENCIES, TAKEN FROM
!     LOPEZ (2002), ARE MERGED).

! - TEMPORAIRE(S) 1D.

! ZQIST     : VALEUR COURRANTE DE QI.
!           : RUNNING VALUE OF QI.
! ZQLST     : VALEUR COURRANTE DE QL.
!           : RUNNING VALUE OF QL.
! ZTST      : VALEUR COURRANTE DE T.
!           : RUNNING VALUE OF T.
! ZDQ       : COPIE LOCALE DE PDQ.
!           : LOCAL COPY OF PDQ.
! ZIPOI     : COPIE LOCALE DE PIPOI.
!           : LOCAL COPY OF PIPOI.
! ZPOID     : COPIE LOCALE DE PPOID.
!           : LOCAL COPY OF PPOID.
! ZZPRSF    : COPIE LOCALE DE PAPRSF.
!           : LOCAL COPY OF PAPRSF.

! ZEXPN     : FACTEUR COMMUN DE DEPENDANCE EN TEMPERATURE (CROISSANT AVEC T).
!           : COMMON TEMPERATURE DEPENDENCY FACTOR (INCREASING WITH T).
! ZIPB      : INVERSE DE LA PRESSION A LA BASE DE LA COUCHE.
!           : INVERSE PRESSURE AT THE BOTTOM OF THE LAYER.
! ZLSCP     : CHALEUR LATENTE DE FUSION DIVISEE PAR CP.
!           : MELTING LATENT HEAT DIVIDED BY CP.


  DO JLON=KIDIA,KFDIA

!     COPIES LOCALES.
!     LOCAL COPIES.

    ZMUL=1.0_JPRB/MAX(ZEPS1,ZNEBLOC(JLON))
    ZQLST(JLON)=PQLMP(JLON,JLEV)*ZMUL
    ZQIST(JLON)=PQIMP(JLON,JLEV)*ZMUL
    ZTST(JLON)=PTMP(JLON,JLEV)
    ZDQ(JLON)=0.0_JPRB
    ZIPOI(JLON)=PIPOI(JLON,JLEV)
    ZPOID(JLON)=PPOID(JLON,JLEV)
    ZZPRSF(JLON)=PAPRSF(JLON,JLEV)
    ZZRHO(JLON)=PAPRSF(JLON,JLEV)/(PR(JLON,JLEV)*PTMP(JLON,JLEV))
    ZRPCEFF(JLON)=REVAPN

    IF (LAB12) THEN
      ZRPCEFF(JLON)=0.9694_JPRB*REVAPN
    ENDIF

!     DEPENDANCES EN TEMPERATURE (POUR LA NEIGE ET LA GLACE) & AUTRES CHAMPS 1D.
!     TEMPERATURE DEPENDENCIES (FOR SNOW AND ICE) AND OTHER 1D ARRAYS.

    ZEXPA(JLON)=(1.0_JPRB-TANH(ZALPHA*(PTMP(JLON,JLEV)-RTT)+ZBETA))*0.5_JPRB
    ZEXPN(JLON)=MIN(1.0_JPRB,EXP(ZEXTMP*(PTMP(JLON,JLEV)-RTT)))
    ZIPB(JLON)=1.0_JPRB/PAPRS(JLON,JLEV)
    ZIPBR(JLON)=ZIPB(JLON)

    IF (LAB12) THEN
      ZIPBR(JLON)=ZIPB(JLON)*SQRT(SQRT(ZIPB(JLON)))
    ENDIF

    ZLSCP(JLON)=(PLHS(JLON,JLEV)-PLHV(JLON,JLEV))/PCP(JLON,JLEV)

  ENDDO

!     AUTOCONVERSION (INCLUANT LE PROCESSUS DE WEGENER-BERGERON-FINDEISEN).
!     AUTOCONVERSION (INCLUDING THE WEGENER-BERGERON-FINDEISEN PROCESS).

! - TEMPORAIRE(S) 1D.

! ZPART     : COPIE LOCALE DE LA NEBULOSITE.
!           : LOCAL COPY OF CLOUDINESS.
! ZQNST     : VALEUR COURRANTE DE QN.
!           : RUNNING VALUE OF QN.
! ZQRST     : VALEUR COURRANTE DE QR.
!           : RUNNING VALUE OF QR.

! ZACONI    : INCREMENT D'AUTOCONVERSION GLACE-NEIGE.
!           : ICE-SNOW AUTOCONVERSION INCREMENT.
! ZACORL    : INCREMENT D'AUTOCONVERSION EAU-PLUIE.
!           : WATER-RAIN AUTOCONVERSION INCREMENT.
! ZACONL    : INCREMENT D'AUTOCONVERSION EAU-NEIGE (PROCESSUS W-B-F).
!           : WATER-SNOW AUTOCONVERSION INCREMENT (W-B-F PROCESS).

  DO JLON=KIDIA,KFDIA
    ZQRST(JLON)=0.0_JPRB
    ZQNST(JLON)=0.0_JPRB
    ZPART(JLON)=ZNEBLOC(JLON)
    ZLSM(JLON)=0.0_JPRB
    ZNEIJ(JLON)=0.0_JPRB
    ZTS(JLON)=PTMP(JLON,KLEV)
  ENDDO

  CALL ACACON(YDML_PHY_MF,KIDIA,KFDIA,KLON,&
!   INPUT
    &           LLSIMP,&
    &           ZEXPA,ZEXPN,ZLSCP,ZPART,ZLSM,ZNEIJ,ZTS,&
!   INPUT/OUTPUT
    &           ZQIST,ZQLST,ZQNST,ZQRST,ZQGST,ZTST,&
!   OUTPUT
    &           ZACONI,ZACONL,ZACORL,ZACOGI,ZACOGL)

!     SEDIMENTATION STATISTIQUE.
!     STATISTICAL SEDIMENTATION.

  DO JLON=KIDIA,KFDIA
    ZZACONI(JLON)=ZACONI(JLON)
    ZZACORL(JLON)=ZACORL(JLON)
    ZZACONL(JLON)=ZACONL(JLON)

    ZACONI(JLON)=ZACONI(JLON)*ZPART(JLON)
    ZACORL(JLON)=ZACORL(JLON)*ZPART(JLON)
    ZACONL(JLON)=ZACONL(JLON)*ZPART(JLON)

    ZFPLSL(JLON,JLEV)=ZFPLSL(JLON,JLEV-1)+PPOID(JLON,JLEV)*ZACORL(JLON)
    ZFPLSN(JLON,JLEV)=ZFPLSN(JLON,JLEV-1)+PPOID(JLON,JLEV)&
     &*(ZACONL(JLON)+ZACONI(JLON))

    IF (LGRAPRO) THEN
      ZZACOGI(JLON)=ZACOGI(JLON)
      ZZACOGL(JLON)=ZACOGL(JLON)
      ZACOGI(JLON)=ZACOGI(JLON)*ZPART(JLON)
      ZACOGL(JLON)=ZACOGL(JLON)*ZPART(JLON)
      ZFPLSG(JLON,JLEV)=ZFPLSG(JLON,JLEV-1)+PPOID(JLON,JLEV)&
       &*(ZACOGL(JLON)+ZACOGI(JLON))
    ENDIF
  ENDDO

!     CALCULS D'EVAPORATION PUIS DE FONTE DES FLUX DE PRECIPITATIONS.

!     EVAPORATION AND MELTING COMPUTATIONS FOR THE PRECIPITATION FLUXES.

! - TEMPORAIRE(S) 1D.

!           : RAIN EVAPORATION INCREMENT.
! ZFONTN     : INCREMENT DE FONTE DE LA NEIGE (-GEL DE LA PLUIE).
!           : SNOW MELT (-RAIN FREEZING) INCREMENT.
! ZTEST     : TEST POUR LA PERPETUATION DU FLUX DE PRECIPITATIONS.
!           : TEST FOR THE PERPETUATION OF THE PRECIPITATION FLUX.

! PASS (CLOUDY PART AND SEEDED - FOR THE MELTING)

  DO JLON=KIDIA,KFDIA
    ZPART(JLON)=ZNEBLOC(JLON)+(1.0_JPRB-ZNEBLOC(JLON))*ZPRPLE(JLON)
    ZW=ZNEBLOC(JLON)/MAX(ZEPS1,ZPART(JLON))
    ZHPLSL(JLON)=ZIPLSLO(JLON)*ZPRPLO(JLON)*ZW+ZIPLSLE(JLON)*&
     &(1.0_JPRB-ZW) 
    ZHPLSN(JLON)=ZIPLSNO(JLON)*ZPRPLO(JLON)*ZW+ZIPLSNE(JLON)*&
     &(1.0_JPRB-ZW) 
    ZLPLSL(JLON)=ZHPLSL(JLON)+PPOID(JLON,JLEV)*ZW*ZZACORL(JLON)
    ZLPLSN(JLON)=ZHPLSN(JLON)+PPOID(JLON,JLEV)*ZW&
     &*(ZZACONL(JLON)+ZZACONI(JLON))

    IF (LGRAPRO) THEN
      ZHPLSG(JLON)=ZIPLSGO(JLON)*ZPRPLO(JLON)*ZW+ZIPLSGE(JLON)*&
       &(1.0_JPRB-ZW)
      ZLPLSG(JLON)=ZHPLSG(JLON)+PPOID(JLON,JLEV)*ZW&
       &*(ZZACOGL(JLON)+ZZACOGI(JLON))
    ENDIF

    ZDQ(JLON)=0.0_JPRB
    ZQRST(JLON)=0.0_JPRB
    ZQNST(JLON)=0.0_JPRB
    ZQGST(JLON)=0.0_JPRB
    ZQVST(JLON)=0.0_JPRB
    ZTST(JLON)=ZTST(JLON)+PTCORR(JLON,JLEV)
    ZSTAL3(JLON)=1.0_JPRB
    ZSTAN3(JLON)=1.0_JPRB
    IF (LGRAPRO) THEN
      ZSTAG3(JLON) = 1.0_JPRB
    ENDIF
  ENDDO

  CALL ACEVMEL (YDML_PHY_MF,KIDIA,KFDIA,KLON,&
!   INPUT
    &            LLSIMP,&
    &            ZDQ,ZEXPN,ZHPLSL,ZHPLSN,ZHPLSG,ZIPB,ZIPBR,ZIPH,ZIPHR,ZIPOI,&
    &            ZLPLSL,ZLPLSN,ZLPLSG,ZPOID,ZQNST,ZQRST,ZQGST,ZQVST,ZSTAL3,ZSTAN3,ZSTAG3,ZTST,&
    &            ZZPRSF,ZLSCP,ZZRHO,ZRPCEFF,&
!   OUTPUT
    &            ZEVAN,ZEVAR,ZEVAG,ZFONTN,ZFONTG,ZTEST)

  DO JLON=KIDIA,KFDIA
    IF (LGRAPRO) THEN
      ZOPRSL=MAX(0.0_JPRB,ZLPLSL(JLON)+ZPOID(JLON)*(ZFONTN(JLON)+ZFONTG(JLON)))/&
       &      MAX(ZEPS1,ZLPLSL(JLON))
      ZOPRSG=MAX(0.0_JPRB,ZLPLSG(JLON)-ZPOID(JLON)*ZFONTG(JLON))/&
       &      MAX(ZEPS1,ZLPLSG(JLON))

      ZOPLSGE(JLON)=ZOPRSG*ZIPLSGE(JLON)

      ZOPLSGO(JLON)=ZOPRSG*(ZIPLSGO(JLON)*ZPRPLO(JLON)+PPOID(JLON,JLEV)&
       &      *(ZZACOGL(JLON)+ZZACOGI(JLON)))

      ZFONTG(JLON)=ZFONTG(JLON)*ZPART(JLON)
    ELSE
      ZOPRSL=MAX(0.0_JPRB,ZLPLSL(JLON)+ZPOID(JLON)*ZFONTN(JLON))/&
       &      MAX(ZEPS1,ZLPLSL(JLON))
    ENDIF
    ZOPRSN=MAX(0.0_JPRB,ZLPLSN(JLON)-ZPOID(JLON)*ZFONTN(JLON))/&
     &      MAX(ZEPS1,ZLPLSN(JLON))


    ZOPLSLE(JLON)=ZOPRSL*ZIPLSLE(JLON)
    ZOPLSNE(JLON)=ZOPRSN*ZIPLSNE(JLON)

    ZOPLSLO(JLON)=ZOPRSL*(ZIPLSLO(JLON)*ZPRPLO(JLON)+PPOID(JLON,JLEV)&
     &      *ZZACORL(JLON))
    ZOPLSNO(JLON)=ZOPRSN*(ZIPLSNO(JLON)*ZPRPLO(JLON)+PPOID(JLON,JLEV)&
     &      *(ZZACONL(JLON)+ZZACONI(JLON)))

    ZFONTN(JLON)=ZFONTN(JLON)*ZPART(JLON)
  ENDDO


!     VALEURS DES PSEUDO-FLUX, FIN DU CALCUL DE SEDIMENTATION, VALEURS DES
!     FLUX ET ROTATION A LA FIN DU CALCUL DE COUCHE.

!     VALUES OF PSEUDO-FLUXES, END OF THE SEDIMENTATION COMPUTATION, VALUES
!     OF THE FLUXES AND SWAP AT THE END OF THE LAYER'S COMPUTATIONS.

  DO JLON=KIDIA,KFDIA
    PMELNET(JLON,JLEV)=ZFONTN(JLON)-ZACONL(JLON)
    IF (LGRAPRO) THEN
      PMELNET(JLON,JLEV)=PMELNET(JLON,JLEV)+ZFONTG(JLON)-ZACOGL(JLON)
      ZFPLSL(JLON,JLEV)=MAX(0.0_JPRB,ZFPLSL(JLON,JLEV)+(ZFONTN(JLON)+ZFONTG(JLON))&
       &*PPOID(JLON,JLEV))
      ZFPLSG(JLON,JLEV)=MAX(0.0_JPRB,ZFPLSG(JLON,JLEV)-ZFONTG(JLON)&
       &*PPOID(JLON,JLEV))
    ELSE
      ZFPLSL(JLON,JLEV)=MAX(0.0_JPRB,ZFPLSL(JLON,JLEV)+ZFONTN(JLON)&
       &*PPOID(JLON,JLEV))
    ENDIF
    ZFPLSN(JLON,JLEV)=MAX(0.0_JPRB,ZFPLSN(JLON,JLEV)-ZFONTN(JLON)&
     &*PPOID(JLON,JLEV))
    ZIPH(JLON)=ZIPB(JLON)
    ZIPHR(JLON)=ZIPBR(JLON)
  ENDDO


!     RECOUVREMENT DES ZONES NUAGEUSES ET PLUVIEUSES.
!     OVERLAP OF CLOUDY AND RAINY AREAS.

  IF(JLEV /= KLEV) THEN

    IF (LGRAPRO) THEN
      DO JLON=KIDIA,KFDIA
        ZCLOV=EXP(-(PAPRSF(JLON,JLEV+1)-PAPRSF(JLON,JLEV))/PDECRD(JLON))

        ZZPRPLE=ZPRPLE(JLON)
        ZZNEBLOC=ZNEBLOC(JLON)
        ZNEBLOC(JLON)=PNEBM(JLON,JLEV+1)
        ZNEBLOCO=ZCLOV*ZNEBLOC(JLON)
        ZUMNLOC=1.0_JPRB-ZNEBLOC(JLON)

      ! UPDATING OF SEEDED FRACTIONS
      ! ALPHA 
        ZMUL=1.0_JPRB/MAX(ZEPS1,(1.0_JPRB-ZNEBLOCO))
        ZPRPLE(JLON)=((ZNEBLOCO+ZZNEBLOC-MIN(ZNEBLOCO,ZZNEBLOC))&
         &*(1.0_JPRB-ZZPRPLE)-ZNEBLOCO+ZZPRPLE)*ZMUL
        ZPRPLE(JLON)=MAX(0.0_JPRB,MIN(1.0_JPRB,ZPRPLE(JLON)))
      ! BETA
        ZBETABR1=(1.0_JPRB-ZCLOV)*(ZZNEBLOC+ZZPRPLE*(ZNEBLOC(JLON)&
         &-ZZNEBLOC))
        ZBETABR2=MIN(ZNEBLOCO,ZZNEBLOC)*(1.0_JPRB-ZZPRPLE)+ZZPRPLE&
         &*ZNEBLOC(JLON)
        ZNEBLI=1.0_JPRB/MAX(ZEPS1,ZNEBLOC(JLON))
        ZPRPLO(JLON)=ZBETABR1*ZMUL+ZBETABR2*(ZNEBLI-(1.0_JPRB-ZCLOV)*ZMUL)
        ZPRPLO(JLON)=MAX(0.0_JPRB,MIN(1.0_JPRB,ZPRPLO(JLON)))

      ! UPDATING OF PARTIAL FLUXES IN SEEDED AREAS:
        ZIPLSL=MAX(0.0_JPRB,ZFPLSL(JLON,JLEV))
        ZIPLSN=MAX(0.0_JPRB,ZFPLSN(JLON,JLEV))
        ZIPLSG=MAX(0.0_JPRB,ZFPLSG(JLON,JLEV))

        ZMUL=1.0_JPRB/MAX(ZEPS1,ZNEBLOC(JLON)*ZPRPLO(JLON))
        ZIPLSLO(JLON)=MAX(0.0_JPRB,(ZIPLSL-ZOPLSLE(JLON)*ZPRPLE(JLON)&
         &*ZUMNLOC))*ZMUL
        ZIPLSNO(JLON)=MAX(0.0_JPRB,(ZIPLSN-ZOPLSNE(JLON)*ZPRPLE(JLON)&
         &*ZUMNLOC))*ZMUL
        ZIPLSGO(JLON)=MAX(0.0_JPRB,(ZIPLSG-ZOPLSGE(JLON)*ZPRPLE(JLON)&
          &*ZUMNLOC))*ZMUL

      ! TWO CASES, DISTINGUISHED BY TEST ON SIGN(ZCLOV*ZNEBLOC-ZZNEBLOC)
        ZZTEST=MAX(0.0_JPRB,SIGN(1.0_JPRB,ZNEBLOCO-ZZNEBLOC))

        ZMUL=1.0_JPRB/MAX(ZEPS1,ZBETABR1+ZUMNLOC*ZBETABR2*ZNEBLI)
        ZIPLSLO(JLON)=ZZTEST*ZIPLSLO(JLON)+(1.0_JPRB-ZZTEST)*ZMUL&
         &*((ZCLOV*ZUMNLOC+(1.0_JPRB-ZCLOV)*ZZNEBLOC)&
         &*ZOPLSLO(JLON)&
         &+(ZZPRPLE*(1.0_JPRB-ZCLOV)*(1.0_JPRB-ZZNEBLOC))&
         &*ZOPLSLE(JLON))
        ZIPLSNO(JLON)=ZZTEST*ZIPLSNO(JLON)+(1.0_JPRB-ZZTEST)*ZMUL&
         &*((ZCLOV*ZUMNLOC+(1.0_JPRB-ZCLOV)*ZZNEBLOC)&
         &*ZOPLSNO(JLON)&
         &+(ZZPRPLE*(1.0_JPRB-ZCLOV)*(1.0_JPRB-ZZNEBLOC))&
         &*ZOPLSNE(JLON))
        ZIPLSGO(JLON)=ZZTEST*ZIPLSGO(JLON)+(1.0_JPRB-ZZTEST)*ZMUL&
         &*((ZCLOV*ZUMNLOC+(1.0_JPRB-ZCLOV)*ZZNEBLOC)&
         &*ZOPLSLO(JLON)&
         &+(ZZPRPLE*(1.0_JPRB-ZCLOV)*(1.0_JPRB-ZZNEBLOC))& 
         &*ZOPLSGE(JLON))

        ZMUL=1.0_JPRB/MAX(ZEPS1,ZUMNLOC*ZPRPLE(JLON))
        ZIPLSLE(JLON)=MAX(0.0_JPRB,(ZIPLSL-ZIPLSLO(JLON)*ZPRPLO(JLON)&
         &*ZNEBLOC(JLON)))*ZMUL
        ZIPLSNE(JLON)=MAX(0.0_JPRB,(ZIPLSN-ZIPLSNO(JLON)*ZPRPLO(JLON)&
         &*ZNEBLOC(JLON)))*ZMUL
        ZIPLSGE(JLON)=MAX(0.0_JPRB,(ZIPLSG-ZOPLSGO(JLON)*ZPRPLO(JLON)&
          &*ZNEBLOC(JLON)))*ZMUL

        ZIPLSLE(JLON)=ZZTEST*ZOPLSLE(JLON)+(1.0_JPRB-ZZTEST)*ZIPLSLE(JLON)
        ZIPLSNE(JLON)=ZZTEST*ZOPLSNE(JLON)+(1.0_JPRB-ZZTEST)*ZIPLSNE(JLON)
        ZIPLSGE(JLON)=ZZTEST*ZOPLSGE(JLON)+(1.0_JPRB-ZZTEST)*ZIPLSGE(JLON)
      ENDDO
    ELSE
      DO JLON=KIDIA,KFDIA

        ! INTRODUCING THE DECORRELATION DEPTH

        ZCLOV=EXP(-(PAPRSF(JLON,JLEV+1)-PAPRSF(JLON,JLEV))/PDECRD(JLON))

        ZZPRPLE=ZPRPLE(JLON)
        ZZNEBLOC=ZNEBLOC(JLON)
        ZNEBLOC(JLON)=PNEBM(JLON,JLEV+1)
        ZNEBLOCO=ZCLOV*ZNEBLOC(JLON)
        ZUMNLOC=1.0_JPRB-ZNEBLOC(JLON)

        ! UPDATING OF SEEDED FRACTIONS
        ! ALPHA 
        ZMUL=1.0_JPRB/MAX(ZEPS1,(1.0_JPRB-ZNEBLOCO))
        ZPRPLE(JLON)=((ZNEBLOCO+ZZNEBLOC-MIN(ZNEBLOCO,ZZNEBLOC))&
         &*(1.0_JPRB-ZZPRPLE)-ZNEBLOCO+ZZPRPLE)*ZMUL
        ZPRPLE(JLON)=MAX(0.0_JPRB,MIN(1.0_JPRB,ZPRPLE(JLON)))
        ! BETA
        ZBETABR1=(1.0_JPRB-ZCLOV)*(ZZNEBLOC+ZZPRPLE*(ZNEBLOC(JLON)&
         &-ZZNEBLOC))
        ZBETABR2=MIN(ZNEBLOCO,ZZNEBLOC)*(1.0_JPRB-ZZPRPLE)+ZZPRPLE&
         &*ZNEBLOC(JLON)
        ZNEBLI=1.0_JPRB/MAX(ZEPS1,ZNEBLOC(JLON))
        ZPRPLO(JLON)=ZBETABR1*ZMUL+ZBETABR2*(ZNEBLI-(1.0_JPRB-ZCLOV)*ZMUL)
        ZPRPLO(JLON)=MAX(0.0_JPRB,MIN(1.0_JPRB,ZPRPLO(JLON)))

        ! UPDATING OF PARTIAL FLUXES IN SEEDED AREAS:
        ZIPLSL=MAX(0.0_JPRB,ZFPLSL(JLON,JLEV))
        ZIPLSN=MAX(0.0_JPRB,ZFPLSN(JLON,JLEV))

        ZMUL=1.0_JPRB/MAX(ZEPS1,ZNEBLOC(JLON)*ZPRPLO(JLON))
        ZIPLSLO(JLON)=MAX(0.0_JPRB,(ZIPLSL-ZOPLSLE(JLON)*ZPRPLE(JLON)&
         &*ZUMNLOC))*ZMUL
        ZIPLSNO(JLON)=MAX(0.0_JPRB,(ZIPLSN-ZOPLSNE(JLON)*ZPRPLE(JLON)&
         &*ZUMNLOC))*ZMUL

        ! TWO CASES, DISTINGUISHED BY TEST ON SIGN(ZCLOV*ZNEBLOC-ZZNEBLOC)
        ZZTEST=MAX(0.0_JPRB,SIGN(1.0_JPRB,ZNEBLOCO-ZZNEBLOC))

        ZMUL=1.0_JPRB/MAX(ZEPS1,ZBETABR1+ZUMNLOC*ZBETABR2*ZNEBLI)
        ZIPLSLO(JLON)=ZZTEST*ZIPLSLO(JLON)+(1.0_JPRB-ZZTEST)*ZMUL&
         &*((ZCLOV*ZUMNLOC+(1.0_JPRB-ZCLOV)*ZZNEBLOC)&
         &*ZOPLSLO(JLON)&
         &+(ZZPRPLE*(1.0_JPRB-ZCLOV)*(1.0_JPRB-ZZNEBLOC))&
         &*ZOPLSLE(JLON))
        ZIPLSNO(JLON)=ZZTEST*ZIPLSNO(JLON)+(1.0_JPRB-ZZTEST)*ZMUL&
         &*((ZCLOV*ZUMNLOC+(1.0_JPRB-ZCLOV)*ZZNEBLOC)&
         &*ZOPLSNO(JLON)&
         &+(ZZPRPLE*(1.0_JPRB-ZCLOV)*(1.0_JPRB-ZZNEBLOC))&
         &*ZOPLSNE(JLON))

        ZMUL=1.0_JPRB/MAX(ZEPS1,ZUMNLOC*ZPRPLE(JLON))
        ZIPLSLE(JLON)=MAX(0.0_JPRB,(ZIPLSL-ZIPLSLO(JLON)*ZPRPLO(JLON)&
         &*ZNEBLOC(JLON)))*ZMUL
        ZIPLSNE(JLON)=MAX(0.0_JPRB,(ZIPLSN-ZIPLSNO(JLON)*ZPRPLO(JLON)&
         &*ZNEBLOC(JLON)))*ZMUL

        ZIPLSLE(JLON)=ZZTEST*ZOPLSLE(JLON)+(1.0_JPRB-ZZTEST)*ZIPLSLE(JLON)
        ZIPLSNE(JLON)=ZZTEST*ZOPLSNE(JLON)+(1.0_JPRB-ZZTEST)*ZIPLSNE(JLON)

      ENDDO
    ENDIF
  ENDIF
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('APLMINI',1,ZHOOK_HANDLE)
END SUBROUTINE APLMINI 


