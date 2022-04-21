!OPTIONS XOPT(NOEVAL)
SUBROUTINE ACNEBCOND_NEW ( YDCST, YDRIP,YDML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,LDREDPR,&
 !-----------------------------------------------------------------------
 ! - INPUT  1D VERTICAL.
 & PHUC,PVETAF,&
 !-----------------------------------------------------------------------
 ! - INPUT  2D .
 & PAPHI,PAPHIF,PAPRSF,PCP,PR,PDELP,PRH,PBLH,PQ,PQI,PQL,PQW,PT,&
 ! - INPUT-OUTPUT  2D .
 & PNCV,&
 ! - INPUT  1D .
 & PGM,PTS,&
 ! - OUTPUT 2D .
 & PQCS,PNEBCOND,PHCRICS,PRHOUT,PQSATS,PRMF,PQCS0,PNEBS0)

!**** *ACNEBCOND * - CALCUL DE NEBULOSITE ET HUMIDITE CRITIQUE POUR LA
!                    CONDENSATION RESOLUE

!     Sujet.
!     ------
!     - ROUTINE DE CALCUL ACTIF .
!       CALCUL DE NEBULOSITE ET HUMIDITE CRITIQUE POUR LA CONDENSATION RESOLUE.
!     - COMPUTATION OF CLOUDINESS AND CRITICAL HUMIDITY FOR THE RESOLVED
!              CONDENSATION .

!**   Interface.
!     ----------
!        *CALL* *ACNEBCOND*

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

! LOGICAL KEY
! LDREDPR    : T FOR REDUCED PROTECTION, F for FULL PROTECTION OF PNCV 

! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
!   CATEGORIE).

! - 1D (1:KLEV) .

! PHUC       : PROFIL DE BASE POUR L'HUMIDITE RELATIVE CRITIQUE.
! PVETAF     : LA COORDONNEE VERTICALE AUX COUCHES.

! - 2D (0:KLEV) .

! PAPHI      : HALF LEVEL GEOPOTENTIAL.

! - 2D (1:KLEV) .

! PAPHIF     : FULL LEVEL GEOPOTENTIAL.
! PAPRSF     : PRESSION AUX NIVEAUX PLEIN.
! PCP        : CHALEUR MASSIQUE A PRESSION CONSTANTE.
! PR         : CONSTANTE DE L'AIR HUMIDE.
! PDELP      : EPAISSEUR EN PRESSION DE LA COUCHE.
! PQ         : HUMIDITE SPECIFIQUE DE LA VAPEUR D'EAU.
! PQW        : HUMIDITE SPECIFIQUE DU THERMOMETRE MOUILLE.
! PT         : TEMPERATURE.
! PNCV       : CONVECTIVE CLOUDINESS (HISTORIC, INPUT)
!               -> EFFECTIVELY PROTECTED FRACTION (OUTPUT) 


! ADDITIONAL PROGNOSTIC VARIABLES

! PQI        : RATIO OF SUSPENDED ICE.
! PQL        : RATIO OF SUSPENDED LIQUID WATER.

! - 1D .

! PGM        : FACTEUR D'ECHELLE.
! PTS        : TEMPERATURE DE SURFACE.

!-----------------------------------------------------------------------

! -   ARGUMENTS DE SORTIE.
!     --------------------

! - 2D (KLEV) .

! PQCS       : CONTENU STRATIFORME EN CONDENSAT NUAGEUX (option SMITH).
! PNEBCOND   : NEBULOSITE POUR LA CONDENSATION RESOLUE.
! PHCRICS    : HUMIDITE CRITIQUE RELATIVE POUR LA CONDENSATION RESOLUE.
! PQSATS     : HUMIDITE SPECIFIQUE DE SATURATION.
! PRMF       : PROPORTION DE LA GLACE.
! PQCS0      : CONTENU STRATIFORME EN CONDENSAT NUAGEUX POUR RAYONNEMENT.
! PNEBS0     : NEBULOSITE PARTIELLE STRATIFORME POUR RAYONNEMENT.

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

!     Methode.
!     --------

!     Auteur.
!     -------
!        07-02, R. Brozkova.

!     Modifications.
!     --------------
!        09-10, L. Bengtsson: Rasch Kristjansson (RK) scheme
!        K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!        2011-06: M. Jerczynski - some cleaning to meet norms
!        K-I Ivarsson 2011-05: Some modifications to RK-scheme
!        2012-02: R. Brozkova: Smith scheme and output for radiation.
!        2012-06: R. Brozkova: XR scheme - better RHcrit dependency on dx.
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!        R. Brozkova (Oct 2014) parameters QXRAL_ADJ and ADJTAU for XR option.
!        L. Gerard (Apr 2016) XR: reduced protection LDREDPR still allowing 
!                                 condensation; simplified PQCS calculation 
!                                 (remove unnecessary phase separation).
!                             SMG: Fix uninitialized ZLV, ZLS.
!        P. Marguinaud (Oct 2016) : Port to single precision
!        R. Brozkova (Sep 2018): Fixes in thermodynamic adjustment - deep
!                                convective condensates protection.
!-------------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE PARKIND1             , ONLY : JPIM     ,JPRB      ,JPRD
USE YOMHOOK              , ONLY : LHOOK,   DR_HOOK

USE YOMCST               , ONLY : TCST
USE YOMRIP               , ONLY : TRIP
!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
TYPE (TCST), INTENT (IN) :: YDCST
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA
LOGICAL, INTENT(IN)              :: LDREDPR
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHUC(KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVETAF(KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHI(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PR(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQL(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQW(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNCV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQCS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNEBCOND(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PHCRICS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQSATS(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRMF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBLH(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRHOUT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQCS0(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PNEBS0(KLON,KLEV)







INTEGER(KIND=JPIM) :: JLEV, JLON, IITER


REAL(KIND=JPRB) :: ZEPS1, ZEPS2, ZEPS3, ZEPS4, ZEPS5, &
 & ZQXRAL, &
 & ZLESEFR, ZLESEFS, &
 & &
 & &
 & &
 & &
 & &
 & ZWEIGHT, ZSQRT6, ZRDSRV

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "fctdoi.ycst.h"
#include "fcttrm.ycst.h"

!-----------------------------------------------------------------------

#include "acnebsm.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACNEBCOND',0,ZHOOK_HANDLE)
ASSOCIATE(REFLRHC=>YDML_PHY_MF%YRPHY0%REFLRHC, RDPHIC=>YDML_PHY_MF%YRPHY0%RDPHIC, &
 & QXRAL=>YDML_PHY_MF%YRPHY0%QXRAL, QXRCDIL=>YDML_PHY_MF%YRPHY0%QXRCDIL, &
 & ADJTAU=>YDML_PHY_MF%YRPHY0%ADJTAU, &
 & NSMTPB=>YDML_PHY_MF%YRPHY0%NSMTPB, RHCEXPDX=>YDML_PHY_MF%YRPHY0%RHCEXPDX, &
 & RDTFAC=>YDML_PHY_MF%YRPHY0%RDTFAC, NSMTPA=>YDML_PHY_MF%YRPHY0%NSMTPA, RSMDNEBX=>YDML_PHY_MF%YRPHY0%RSMDNEBX, &
 & RETAMIN=>YDML_PHY_MF%YRPHY0%RETAMIN, GRHCMOD=>YDML_PHY_MF%YRPHY0%GRHCMOD, RHCRIT1=>YDML_PHY_MF%YRPHY0%RHCRIT1, &
 & TEQH=>YDML_PHY_MF%YRPHY0%TEQH, RHCRIT2=>YDML_PHY_MF%YRPHY0%RHCRIT2, QXRAL_ADJ=>YDML_PHY_MF%YRPHY0%QXRAL_ADJ, &
 & SCLESPS=>YDML_PHY_MF%YRPHY0%SCLESPS, SCLESPR=>YDML_PHY_MF%YRPHY0%SCLESPR, HUCRED=>YDML_PHY_MF%YRPHY0%HUCRED, &
 & TSPHY=>YDML_PHY_MF%YRPHY2%TSPHY, &
 & LSMGCDEV=>YDML_PHY_MF%YRPHY%LSMGCDEV, NSMDNEB=>YDML_PHY_MF%YRPHY%NSMDNEB, L3MT=>YDML_PHY_MF%YRPHY%L3MT, &
 & LSMTPS=>YDML_PHY_MF%YRPHY%LSMTPS, NSMTBOT=>YDML_PHY_MF%YRPHY%NSMTBOT, LXRCDEV=>YDML_PHY_MF%YRPHY%LXRCDEV, &
 & LRKCDEV=>YDML_PHY_MF%YRPHY%LRKCDEV, LSMITH_CDEV=>YDML_PHY_MF%YRPHY%LSMITH_CDEV, &
 & RSTATI=>YDRIP%RSTATI, YDPHY0=>YDML_PHY_MF%YRPHY0)
!-----------------------------------------------------------------------

!*
!     ------------------------------------------------------------------
!     I - CALCUL DE PARAMETRES DERIVES ET CONSTANTES DE SECURITE.

!         COMPUTATION OF DERIVED PARAMETERS AND SECURITY CONSTANTS.

IF (JPRB == JPRD) THEN
  ZEPS1=1.E-14_JPRB
  ZEPS5=1.E-12_JPRB
ELSE
  ZEPS1=1.E-06_JPRB
  ZEPS5=1.E-06_JPRB
ENDIF
ZEPS2=1.E-02_JPRB
ZEPS3=1.E-20_JPRB
ZEPS4=1.E-10_JPRB
ZSQRT6 = SQRT(6._JPRB)
ZRDSRV = YDCST%RD/YDCST%RV
ZWEIGHT=1.0_JPRB-EXP(-TSPHY/ADJTAU)

ZQXRAL=QXRAL_ADJ
ZLESEFR=1._JPRB/SCLESPR**RHCEXPDX
ZLESEFS=1._JPRB/SCLESPS**RHCEXPDX
IITER=3

!*
!     ------------------------------------------------------------------
!     II - CALCULS PRELIMINAIRES DE LA PROPORTION DE LA GLACE.

!          PRELIMINARY COMPUTATIONS OF ICE PROPORTION.


DO JLEV=KTDIA,KLEV
!DEC$ IVDEP
  DO JLON=KIDIA,KFDIA
    PRMF(JLON,JLEV)=FONICE(PT(JLON,JLEV),YDPHY0%RDTFAC)
  ENDDO
ENDDO

!*
!     ------------------------------------------------------------------
!     III - CALCULS DE CONDENSATION-EVAPORATION.

!           CONDENSATION-EVAPORATION COMPUTATIONS.

 ! LXRCDEV

 ! LSMGCDEV

 ! LRKCDEV

 ! LXRCDEV 

 ! LSMGCDEV
 ! LRKCDEV 


!  ------------------------------------------------------------------
!           IVd - CALCUL DE LA NEBULOSITE (FORMULE SMITH).

!           CLOUDINESS DIAGNOSED BY SMITH SCHEME.
! ---------------------------------------------------------------------

!    1. Computation following Smith QJRMS 1990 paper: routine ACNEBSM

  CALL ACNEBSM (YDCST, YDML_PHY_MF%YRPHY0, KIDIA, KFDIA, KLON, KTDIA, KLEV,&
   & PT, PQ, PQL, PQI, &
   & PAPHI, PAPRSF, PCP, PR, &
   & PGM, PVETAF, &
   & PQCS, PNEBCOND, PHCRICS, PRMF ,PQSATS )

  

    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA
        PNEBCOND(JLON,JLEV)=MIN(1.0_JPRB-ZEPS5,MAX(PNEBCOND(JLON,JLEV),ZEPS5))
        PNEBS0(JLON,JLEV)=PNEBCOND(JLON,JLEV)  
        PQCS0(JLON,JLEV)=PQCS(JLON,JLEV)
        PRHOUT(JLON,JLEV)=PQ(JLON,JLEV)/PQSATS(JLON,JLEV)
      ENDDO
    ENDDO 

   ! L3MT
 
 ! LSMITH_CDEV
!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACNEBCOND',1,ZHOOK_HANDLE)
END SUBROUTINE ACNEBCOND_NEW

