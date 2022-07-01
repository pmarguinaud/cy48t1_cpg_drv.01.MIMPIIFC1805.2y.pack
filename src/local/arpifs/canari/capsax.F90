SUBROUTINE CAPSAX(ROBHDR,ROBODY,YDDIMV,YDSURF,KPROMA,KNBPT,KTASK,PSP_CI,PTT,PQT,PSPS,PSLN,&
 & PESIG,PCAGUE,PGELAT,PGELAM,PMORO,PMLSM)
!****-------------------------------------------------------------------
!****   CAPSAX : ROUTINE GERANT L'EXECUTION DE L'ANALYSE DE PRESSION DE SURFACE.
!****-------------------------------------------------------------------
!***  ROUTINES APPELANTES: - CAPOTX -
!***  -------------------
!***  ARGUMENTS:  IN       KPROMA : surdimensionnement horizontal
!***  ---------   IN       KNBPT  : nombre reel de points traites
!***              IN       KTASK  : numero de la tache (multitasking)
!***              OUT      PSP_CI :
!***              IN       PTT, PQT :
!***              IN/OUT   PSPS, PSLN :
!***              IN       PESIG, PCAGUE, PGELAT, PGELAM : 

!     COMMONS: - QADOCK - CONTIENT LA TABLE DE CONTINGENCE
!     -------             PREDICTEURS/PREDICTANTS ET DES CONTRAINTES
!                         POUR LA SELECTION DES OBSERVATIONS.
!              - QAPAVU - DIMENSIONS ET INDICE RELATIFS A *QAVARA*.
!              - QADIAG - CHAMPS POINTS DE GRILLE SPECIAUX DE L'ANALYSE.
!              - QACOST - PARAMETRES STATISTIQUES DE L'ANALYSE.
!              - YOMDIM - DIMENSION DES TABLEAUX DE TRAVAIL.
!              - PARDIMO- PARAMETERS POUR LES TABLEAUX D'OBSERVATIONS.

!     EXTERNES: - CANEVA - ROUTINE EFFECTUANT L'ANALYSE.
!     --------  - GPRCP  - ROUTINE CALCULANT R, CP ET KAPPA.
!****--------------------------------------------------------------------

USE YOMDIMV            , ONLY : TDIMV
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE PARDIMO  , ONLY : JPNOTP
USE QAPAVU   , ONLY : JPPXZ
USE QAPREX   , ONLY : NPREA    ,NBPREA   ,NAPS
USE QADOCK   , ONLY : QDSTVA   ,MINMA    ,QCORMIN  ,QDELPI
USE QADIAG   , ONLY : RCATSL   ,RCARES   ,RCASIG
USE QACOST   , ONLY : TYPE_QACOST
USE IFS_DBASE_VIEW_MOD, ONLY: IFS_DBASE_VIEW

!**---------------------------------------------------------------------

IMPLICIT NONE

TYPE(IFS_DBASE_VIEW), INTENT(INOUT) :: ROBHDR,ROBODY
TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KNBPT
INTEGER(KIND=JPIM),INTENT(IN)    :: KTASK
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSP_CI(KPROMA,YDSURF%YSP_CID%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQT(KPROMA)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPS(KPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSLN(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESIG(2*YDDIMV%NFLEVG+5,KNBPT)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCAGUE(2,KNBPT)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAT(KNBPT)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KNBPT)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMORO(KNBPT)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMLSM(KNBPT)
!**---------------------------------------------------------------------
!     INM    : NOMBRE DE PREDICTEURS PAR PREDICTANT.
!     INBPR  : NOMBRE DE PREDICTEURS PAR POINTS D'ANALYSE.
INTEGER(KIND=JPIM) :: INM(1), INBPR(KNBPT)
INTEGER(KIND=JPIM) :: INLXI, INXI(1)

!     ZSIGPS : EN GENERAL, ECART TYPES D'ERREUR DE L'EBAUCHE
!     ZRESIA : RESIDUS ANALYSES.
!     ZSIGA  : ECART TYPES D'ERREUR D'ANALYSE.
!     ZQI    : EAU GLACE.
!     ZQL    : EAU LIQUIDE.
!     ZQR    : PLUIE.
!     ZQS    : NEIGE.
!     ZQG    : GRAUPEL.
!     ZCP    : CP DEPENDANT DES DIFFERENTES PHASES DE L'EAU.
!     ZR     : R DEPENDANT DES DIFFERENTES PHASES DE L'EAU.
!     ZKAP   : R/CP DEPENDANT DES DIFFERENTES PHASES DE L'EAU.

REAL(KIND=JPRB) :: ZDSTVA(JPNOTP), ZSIGPS(KNBPT), ZRESIA(KNBPT), ZSIGA(KNBPT)
REAL(KIND=JPRB) :: ZCP(KNBPT), ZR(KNBPT), ZKAP(KNBPT)
REAL(KIND=JPRB) :: ZQI(KNBPT), ZQL(KNBPT),ZQR(KNBPT), ZQS(KNBPT),ZQG(KNBPT)

INTEGER(KIND=JPIM) :: IANTY, IQUOI, IXX(JPNOTP), JROF

TYPE(TYPE_QACOST) :: YLCOST(KNBPT)

REAL(KIND=JPRB) :: ZCORMIN, ZDELPI, ZXXX(JPNOTP)
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!**---------------------------------------------------------------------

#include "caneva.intfb.h"
#include "gprcp_qlirsg.intfb.h"

!**---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CAPSAX',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & YSP_CI=>YDSURF%YSP_CI, YSP_CID=>YDSURF%YSP_CID)

!**----------------------------------------------------------------------
!**       1. Initialisations.
!**          ---------------

RCATSL(:,1,KTASK) = 0.0_JPRB
RCARES(:,1,KTASK) = 0.0_JPRB
RCASIG(:,1,KTASK) = 1.0_JPRB

ZSIGPS(:) = PESIG(NFLEVG,:)

IQUOI = -1
IANTY = 0

INLXI = 1
INXI(1) = JPPXZ

ZDELPI = QDELPI(NAPS)
ZCORMIN = QCORMIN(NAPS)

ZDSTVA = 0.0_JPRB

ZDSTVA(1) = QDSTVA(1)
ZDSTVA(4) = QDSTVA(4)
ZDSTVA(5) = QDSTVA(5)
ZDSTVA(6) = QDSTVA(6)

INM(1) = MINMA(NAPS)

!**---------------------------------------------------------------------
!**       2. Lancement de l'analyse par interpolation optimale.
!**          -------------------------------------------------

CALL CANEVA(ROBHDR, ROBODY, IQUOI, IANTY, KTASK, ZXXX, ZDSTVA(1), ZDELPI, IXX, IXX,&
 & ZCORMIN, PCAGUE, PGELAT, PGELAM, PMORO, PMLSM, PSLN, KNBPT, 1,&
 & INXI(1), INLXI, NPREA(1,NAPS), NBPREA(NAPS), INM(1),&
 & ZSIGPS(1), ZRESIA(1), ZSIGA(1), INBPR(1), YLCOST)

!**---------------------------------------------------------------------
!**       3. Calcul du champ de pression analyse.
!**          -----------------------------------

CALL GPRCP_QLIRSG(KNBPT,1,KNBPT,1,PQ=PQT,PQI=ZQI,PQL=ZQL,PQR=ZQR,PQS=ZQS,PQG=ZQG, &
 & PCP=ZCP,PR=ZR,PKAP=ZKAP)

DO JROF = 1 , KNBPT
  RCATSL(JROF,1,KTASK) = REAL(INBPR(JROF),JPRB)
  RCARES(JROF,1,KTASK) = ZRESIA(JROF)
  RCASIG(JROF,1,KTASK) = ZSIGA(JROF)/ZSIGPS(JROF)
  PSP_CI(JROF,YSP_CI%YCI(1)%MP0) = ZSIGPS(JROF)
  PSP_CI(JROF,YSP_CI%YCI(2)%MP0) = ZSIGA(JROF)
  PSLN(JROF) = PSLN(JROF) + ZRESIA(JROF)/(ZR(JROF)*PTT(JROF))
  PSPS(JROF) = EXP(PSLN(JROF))
ENDDO

!**---------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CAPSAX',1,ZHOOK_HANDLE)
END SUBROUTINE CAPSAX
