SUBROUTINE ACAJUCV( YDPHY0,KIDIA,  KFDIA, KLON,  KTDIA, KLEV, KLAJU,&
!-----------------------------------------------------------------------
! - INPUT  2D .
& PAPRS, PALPH, PDELP, PLNPR, &
! - INPUT/OUTPUT  2D .
& PT)


!**** *ACAJUCV - AJUSTEMENT CONVECTIF SEC.

!     Sujet.
!     ------
!
!     - ROUTINE DE CALCUL ACTIF .
!       CORRECTION DU PROFIL DE TEMPERATURE PAR AJUSTEMENT
!       CONVECTIF SEC. IL Y A UN NIVEAU "KLAJU" QUI DEFINIT UN
!       NIVEAU AU DESSUS DUQUEL ON A UN AJUSTEMENT "CLASSIQUE"
!       (POUR JLEV<KLAJU) : AJUSTEMENT DES QUE LE PROFIL EST
!       SURADIABATIQUE ; CORRECTION AU NIVEAU DE "AJ1MEPS=1-epsilon".
!       AU DESSOUS DE CE NIVEAU (POUR JLEV>KLAJU) ON PERMET DE PLUS
!       EN PLUS LA PRESENCE DE SURADIABATISME : TEST DE DECLENCHEMENT
!       VARIANT LINEAIREMENT DE "1" A "AJ1PEPS" ENTRE KLAJU ET KLEV ;
!       NIVEAU DE CORRECTION VARIANT DE "AJ1MEPS" A "1" ENTRE KLAJU
!       ET KLEV.
!       > SI LE PROCESSUS ITERATIF NE CONVERGE PAS (NAJITER), ON SORT
!       > AVEC LE PROFIL D'ENTREE NON CORRIGE.

!     - CORRECTION OF THE VERTICAL PROFILE OF TEMPERATURE TO MAKE
!       A DRY CONVECTIVE ADJUSTMENT. FOR JLEV<KLAJU THE ADJUSTMENT
!       IS ACTIVE AS SOON AS THE VERTICAL GRADIENT IS SUPER-
!       ADIABATIC AND AN ITERATIVE CORRECTION IS MADE UP TO
!       "AJ1MEPS=1-epsilon" LEVEL. FOR THE LOWER PART OF THE
!       ATMOSPHERE (FOR THE CLP : JLEV>KLAJU) SUPER-ADIABATIC
!       GRADIENT ARE ALOWED : THE CRITICAL GRADIENT GOES LINEARLY
!       FROM "1" TO "AJ1PEPS"  BETWEEN KLAJU AND KLEV ; THE LEVEL
!       OF THE CORRECTION VARIES LINEARLY FROM "AJ1MEPS" TO "1"
!       BETWEEN KLAJU AND KLEV.
!       > IF THE MAXIMUM ACCOUNT OF ITERATION IS REACHED (NAJITER),
!       > THE INPUT TEMPERATURE ARRAY "PT" IS PUT INTO THE OUTPUT
!       > "PT" ARRAY.

!**   Interface.
!     ----------
!        *CALL* *ACAJUCV*

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
! KLAJU      : INDICE DU NIVEAU POUR LA LIMITATION INF. DE L'AJUSTEMENT
!              (ON PERMET DES SURADIABATISMES DANS LA CLP : JLEV>KLAJU)

! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
!   CATEGORIE).

! - 2D (0:KLEV) .

! PAPRS      : PRESSION AUX DEMI-NIVEAUX.

! - 2D (1:KLEV) .INTEGER

! PALPH      : LOG(PAPRS(JLEV)/PAPRSF(JLEV)).
! PDELP      : EPAISSEUR DE LA COUCHE : PAPRS(JLEV)-PAPRS(JLEV-1)
! PLNPR      : LOG(PAPRS(JLEV)/PAPRS(JLEV-1)).
! PT         : TEMPERATURE.

!-----------------------------------------------------------------------

! -   ARGUMENTS DE SORTIE.
!     --------------------

! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
!   CATEGORIE).

! - 2D (1:KLEV) .

! PT      : LA TEMPERATURE APRES AJUSTEMENT CONVECTIF SEC.
!-----------------------------------------------------------------------

!     Externes.
!     ---------

!     Methode.
!     --------

!     Auteur.
!     -------
!        JF Gueremy, 2009, after
!        98-12, P. Marquet (from the sub-routine AJUCON of  )

!     Modifications.
!     --------------
!        09-05-07, - JF Gueremy. phasage cy32.
!        11-11-21, - JF Gueremy. phasage cy37.
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : RD       ,RCPD
USE YOMLUN   , ONLY : NULOUT
USE YOMPHY0  , ONLY : TPHY0

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TPHY0)       ,INTENT(IN):: YDPHY0
INTEGER(KIND=JPIM),INTENT(IN):: KFDIA
INTEGER(KIND=JPIM),INTENT(IN):: KIDIA
INTEGER(KIND=JPIM),INTENT(IN):: KLAJU
INTEGER(KIND=JPIM),INTENT(IN):: KLEV
INTEGER(KIND=JPIM),INTENT(IN):: KLON
INTEGER(KIND=JPIM),INTENT(IN):: KTDIA

REAL(KIND=JPRB),INTENT(IN)    :: PAPRS(KLON,0:KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PALPH(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PDELP(KLON,KLEV)
REAL(KIND=JPRB),INTENT(IN)    :: PLNPR(KLON,KLEV)
REAL(KIND=JPRB),INTENT(INOUT) :: PT   (KLON,KLEV)

!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZALPHA(KLON,KLEV), ZBETA (KLON,KLEV), ZDELTA(KLON,KLEV)
REAL(KIND=JPRB) :: ZAJCV1(KLON,KLEV), ZAJCV2(KLON,KLEV), ZAJCV3(KLON,KLEV)
REAL(KIND=JPRB) :: ZTAJU(KLON,KLEV),  ZMODIF(KLON)

INTEGER(KIND=JPIM) :: IMODIF, ITDIA, JITER, JLEV, JLON
REAL(KIND=JPRB) :: ZDEN23, ZFACT, ZPROF1, ZPROF2, ZY, ZZAJ1, ZZAJ23, ZZIND, ZZIND4
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACAJUCV',0,ZHOOK_HANDLE)
ASSOCIATE(AJ1PEPS=>YDPHY0%AJ1PEPS, NAJITER=>YDPHY0%NAJITER, &
 & AJ1MEPS=>YDPHY0%AJ1MEPS)
!-----------------------------------------------------------------------

!     ------------------------------------------------------------------
!     I - CALCULS PRELIMINAIRES
!     ------------------------------------------------------------------

!     Tableaux de travail pour l'ajustement convectif sec :
!     (ZALPHA,ZDELTA,ZBETA).

IF (KTDIA == 1) THEN
  DO JLON=KIDIA,KFDIA
    ZALPHA(JLON,KTDIA)=1.0_JPRB
    ZDELTA(JLON,KTDIA)=2.0_JPRB
    ZBETA (JLON,KTDIA)=ZDELTA(JLON,KTDIA)-ZALPHA(JLON,KTDIA)
  ENDDO ! JLON
ENDIF

IF (KTDIA == 1) THEN
  ITDIA=KTDIA+1
ELSE
  ITDIA=KTDIA
ENDIF

DO JLEV=ITDIA,KLEV
  DO JLON=KIDIA,KFDIA
    ZALPHA(JLON,JLEV)=PALPH(JLON,JLEV)
    ZDELTA(JLON,JLEV)=PLNPR(JLON,JLEV)
    ZBETA (JLON,JLEV)=ZDELTA(JLON,JLEV)-ZALPHA(JLON,JLEV)
  ENDDO ! JLON
ENDDO   ! JLEV

! Tableaux de travail pour l'ajustement convectif sec :

!      ZAJCV1 : pour le test de suradiabatisme
!      ZAJCV2 : dans la correction pour JLEV
!      ZAJCV3 : dans la correction pour JLEV+1
!      ZZIND  : constant a  0. pour p < PAPRS(KLAJU)
!               entre 0. et 1. pour PAPRS(KLAJU) < p < PAPRS(KLEV)
!               egal     a  1. pour p = PAPRS(KLEV)

!      Et on prend ZZIND4=ZZIND**4 pour avoir un profil continu
!      et derivable en PAPRS(KLAJU), avec des gradients fortement
!      suradiabatiques preserves uniquement tres pres de la
!      surface p=PAPRS(KLEV). La formule suivante pour ZPROF1
!      "(AJ1PEPS-1.)*ZZIND4+1." varie entre 1. et AJ1PEPS entre
!      PAPRS(KLAJU) et PAPRS(KLEV) en suivant ZZIND4. On autorise
!      ainsi des gradients ZPROF1*RD/CPD (en valeur absolue) qui
!      sont tres suradiabatiques la ou ZPROF1 est tres superieur
!      a 1 (pour AJ1PEPS eleve et dans les basses couches).

ZFACT=RD/RCPD
DO JLEV=KTDIA,KLEV-1
  DO JLON=KIDIA,KFDIA
    ZZIND=MIN(1.0_JPRB,&
     &MAX( (PAPRS(JLON,JLEV)-PAPRS(JLON,KLAJU))/&
     &     (PAPRS(JLON,KLEV)-PAPRS(JLON,KLAJU)), 0.0_JPRB) )
    ZZIND4=ZZIND*ZZIND*ZZIND*ZZIND
    ZPROF1=(AJ1PEPS-1.0_JPRB)*ZZIND4+1.0_JPRB
    ZPROF2=(1.0_JPRB-AJ1MEPS)*ZZIND4+AJ1MEPS
    ZZAJ1 =(1.0_JPRB-ZBETA (JLON,JLEV+1)*ZFACT*ZPROF1)/&
     &     (1.0_JPRB+ZALPHA(JLON,JLEV)  *ZFACT*ZPROF1)
    ZZAJ23=(1.0_JPRB-ZBETA (JLON,JLEV+1)*ZFACT*ZPROF2)/&
     &     (1.0_JPRB+ZALPHA(JLON,JLEV)  *ZFACT*ZPROF2)
    ZDEN23=PDELP(JLON,JLEV+1)+PDELP(JLON,JLEV)*ZZAJ23
    ZAJCV1(JLON,JLEV)=ZZAJ1
    ZAJCV2(JLON,JLEV)=1.0_JPRB/ZDEN23
    ZAJCV3(JLON,JLEV)=ZZAJ23/ZDEN23
  ENDDO ! JLON
ENDDO   ! JLEV

! On remplit le tableau ZTAJU par le profil d'entree PT :
DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    ZTAJU(JLON,JLEV)=PT(JLON,JLEV)
  ENDDO ! JLON
ENDDO   ! JLEV

!     ------------------------------------------------------------------
!     II - CORRECTION DU PROFIL DE TEMPERATURE PAR AJUSTEMENT
!          CONVECTIF SEC. ON REMPLI LE TABLEAU "ZTAJU" QUI
!          SERA GARDE EN CAS D'AJUSTEMENT (ZMODIF=0), OU
!          QUI SERA REMPLACE PAR LE PROFIL D'ENTREE "PT" SINON.
!     ------------------------------------------------------------------

DO JITER=1,NAJITER
  IMODIF=0
  DO JLON=KIDIA,KFDIA
    ZMODIF(JLON)=0.0_JPRB
  ENDDO  ! JLON
  DO JLEV=KLEV-1,KTDIA,-1
    DO JLON=KIDIA,KFDIA
      IF(ZTAJU(JLON,JLEV)  < ZTAJU(JLON,JLEV+1)*ZAJCV1(JLON,JLEV)) THEN
        IMODIF=1
        ZMODIF(JLON)=1.0_JPRB
        ZY=ZTAJU(JLON,JLEV)  *ZAJCV2(JLON,JLEV)-&
         & ZTAJU(JLON,JLEV+1)*ZAJCV3(JLON,JLEV)
        ZTAJU(JLON,JLEV)  =ZTAJU(JLON,JLEV)  -ZY*PDELP(JLON,JLEV+1)
        ZTAJU(JLON,JLEV+1)=ZTAJU(JLON,JLEV+1)+ZY*PDELP(JLON,JLEV)
      ENDIF
    ENDDO ! JLON
  ENDDO   ! JLEV
  IF (IMODIF == 0) GOTO 111
ENDDO     ! JITER

111 CONTINUE

!     ------------------------------------------------------------------
!     III - ON REMPLI LE TABLEAU "PT" EN CAS D'AJUSTEMENT
!           (ZMODIF=0.), OU ON GARDE "PT" SINON (UNIQUEMENT
!           POUR LES COLONNES NON AJUSTEES : ZMODIF=1).
!           ON IMPRIME UN DIAGNOSTIC DE PROBLEME EN CAS DE
!           NON CONVERGENCE.
!     ------------------------------------------------------------------

DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    PT(JLON,JLEV) =   PT (JLON,JLEV) *        ZMODIF(JLON)&
     &               + ZTAJU(JLON,JLEV) * (1.0_JPRB-ZMODIF(JLON))
  ENDDO ! JLON
ENDDO   ! JLEV

DO JLON=KIDIA,KFDIA
  IF(ZMODIF(JLON) /= 0.0_JPRB) THEN
    WRITE(UNIT=NULOUT,FMT='('' !! ACAJUCV : PB CONV. : '',I5)') JLON
  ENDIF
ENDDO ! JLON

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACAJUCV',1,ZHOOK_HANDLE)
END SUBROUTINE ACAJUCV
