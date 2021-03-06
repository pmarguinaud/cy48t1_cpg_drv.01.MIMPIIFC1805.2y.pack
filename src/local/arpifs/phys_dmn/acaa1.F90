!OPTIONS XOPT(NOEVAL)
SUBROUTINE ACAA1 ( YDML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
 !-----------------------------------------------------------------------
 ! - INPUT  2D .
 & PAPRS,PCOEFN,PCP,PQL,PQI,PT,&
 ! -- OUTPUT 2-D for acdifv2.f90 
 & PALPHA1,PCOEFA,PLVT,PQICE)

!!     Sujet.
!     ------
!     - ROUTINE DE CALCUL ACTIF .
!       CALCULS DE ALPHA1 ET DE A*F2(Q1) POUR EXPRIMER LE FLUX
!       TURBULENT DU CONDENSAT EN FONCTION DES FLUX DE THETAL ET QT

!**   Interface.
!     ----------
!        *CALL* *ACAA1*

! -   ARGUMENTS D'ENTREE.
!     -------------------

! - NOM DES PARAMETRES DE DIMENSIONNEMENT DE LA PHYSIQUE.

! KIDIA      : INDICE DE DEPART DES BOUCLES VECTORISEES SUR L'HORIZONT..
! KFDIA      : INDICE DE FIN DES BOUCLES VECTORISEES SUR L'HORIZONTALE.
! KLON       : DIMENSION HORIZONTALE DES TABLEAUX.
! KTDIA      : INDICE DE DEPART DES BOUCLES VERTICALES (1 EN GENERAL).
! KLEV       : DIMENSION VERTICALE DES TABLEAUX "FULL LEVEL".

! - 2D (1:KLEV) .

! PAPRS      : PRESSION AUX DEMI-NIVEAUX.
! PCOEFN     : COEFFICIENT STATISTIQUE POUR LES FLUX D'EAUX CONDENSEES.
! PCP        : CHALEUR MASSIQUE A PRESSION CONSTANTE DE L'AIR.
! PQL        : HUMIDITE SPECIFIQUE DE L'EAU LIQUIDE.
! PQI        : HUMIDITE SPECIFIQUE DE L'EAU GLACE.
! PT         : TEMPERATURE.

!-----------------------------------------------------------------------

! -   ARGUMENTS DE SORTIE.
!     --------------------

! - 2D (0:KLEV) .

! PALPHA1    : ALPHA1
! PCOEFA     : A*F2(Q1)
! PLVT       : CHALEUR LATENTE
! PQICE      : PROPORTION D'EAU SOLIDE
!-----------------------------------------------------------------------

! -   ARGUMENTS IMPLICITES.
!     ---------------------

! COMMON/YOMCST /

!-----------------------------------------------------------------------

!     Auteur.
!     -------
!        2013-01, J.F. Gueremy

!     Modified.
!     ---------
!
!-----------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY :  RV       ,RCPV     ,RETV     ,RCW      ,&
 & RCS      ,RLVTT    ,RLSTT    ,RTT      ,RALPW    ,RBETW    ,&
 & RGAMW    ,RALPS    ,RBETS    ,RGAMS    ,RALPD    ,RBETD    ,&
 & RGAMD    ,RDT

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOEFN(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQL(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQI(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)

! - OUTPUT  2D for acdifv2.f90
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PALPHA1(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PCOEFA(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PLVT(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PQICE(KLON,0:KLEV)

INTEGER(KIND=JPIM) :: JLEV, JLON

REAL(KIND=JPRB) :: ZTH, ZLV, ZLS, ZCPH, ZTLI, ZICE, ZLIQ, &
 & ZEISP, ZELSP, ZQSATI, ZQSATL, ZSRVCPT2, ZDQSATI, ZDQSATL

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "fcttrm.func.h"
#include "fctdoi.func.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACAA1',0,ZHOOK_HANDLE)
ASSOCIATE(RDTFAC=>YDML_PHY_MF%YRPHY0%RDTFAC, YDPHY0=>YDML_PHY_MF%YRPHY0)
!-----------------------------------------------------------------------


DO JLEV=KTDIA,KLEV-1
  DO JLON=KIDIA,KFDIA
    ZTH=(PT(JLON,JLEV)+PT(JLON,JLEV+1))*0.5_JPRB
    ZLV = FOLH(ZTH,0.0_JPRB)
    ZLS = FOLH(ZTH,1.0_JPRB)
    PQICE(JLON,JLEV) = FONICE(ZTH,YDPHY0%RDTFAC)
    ZCPH=(PCP(JLON,JLEV)+PCP(JLON,JLEV+1))*0.5_JPRB
    ZTLI=ZTH-0.5_JPRB*(ZLV*(PQL(JLON,JLEV)+PQL(JLON,JLEV+1))&
     & +ZLS*(PQI(JLON,JLEV)+PQI(JLON,JLEV+1)))/ZCPH
    ZICE=PQICE(JLON,JLEV)
    ZLIQ=1.0_JPRB-ZICE
    PLVT(JLON,JLEV)=ZLV*ZLIQ+ZLS*ZICE
    ZEISP = FOEW(ZTLI,1.0_JPRB)/PAPRS(JLON,JLEV) 
    ZELSP = FOEW(ZTLI,0.0_JPRB)/PAPRS(JLON,JLEV) 
    ZQSATI = FOQS(ZEISP)
    ZQSATL = FOQS(ZELSP)
    ZSRVCPT2=1.0_JPRB/(RV*ZTLI*ZTLI*ZCPH)
    ZDQSATI=ZQSATI*FOLH(ZTLI,1.0_JPRB)*ZSRVCPT2
    ZDQSATL=ZQSATL*FOLH(ZTLI,0.0_JPRB)*ZSRVCPT2
    PALPHA1(JLON,JLEV)=(ZLIQ*ZDQSATL+ZICE*ZDQSATI)
    PCOEFA(JLON,JLEV)=(PCOEFN(JLON,JLEV)+PCOEFN(JLON,JLEV+1))*0.5_JPRB&
     & /((1.0_JPRB+ZDQSATL*ZLV)*ZLIQ&
     &  +(1.0_JPRB+ZDQSATI*ZLS)*ZICE)
  ENDDO
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACAA1',1,ZHOOK_HANDLE)
END SUBROUTINE ACAA1

