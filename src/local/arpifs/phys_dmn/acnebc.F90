!OPTIONS XOPT(NOEVAL)
SUBROUTINE ACNEBC ( YDPHY0,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
 !-----------------------------------------------------------------------
 ! - INPUT  2D .
 & PFPLCL,PFPLCN,KNLAB,&
 ! - INPUT  1D .
 & KNND,&
 ! - OUTPUT 1D .
 & PTCCH,PSCCH,PBCCH)

!**** *ACNEBC * - CALCUL DE LA NEBULOSITE CONVECTIVE.

!     Sujet.
!     ------
!     - ROUTINE DE CALCUL ACTIF .
!       CALCUL DE LA NEBULOSITE CONVECTIVE A PARTIR DU FLUX DE
!       PRECIPITATION ET DE L'INDICE DE STABILITE VERTICALE .
!     - COMPUTATION OF THE CONVECTIVE CLOUDINESS AS A FUNCTION OF THE
!       PRECIPITATION FLUX AND OF THE VERTICAL STABILITY INDEX .

!**   Interface.
!     ----------
!        *CALL* *ACNEBC*

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

! PFPLCL     : PRECIPITATIONS CONVECTIVES SOUS FORME LIQUIDE.
! PFPLCN     : PRECIPITATIONS CONVECTIVES SOUS FORME NEIGE.

! - 2D (1:KLEV) .

! KNLAB      : INDICE D'INSTABILITE CONVECT. (1 DANS LE NUAGE, 0 SINON).

! - 1D (DIAGNOSTIQUE) .

! KNND       : OCCURRENCE DE PRECIPITATIONS CONVECTIVES.

!-----------------------------------------------------------------------

! -   ARGUMENTS DE SORTIE.
!     --------------------

! - 1D (PROGNOSTIQUE) .

! PTCCH      : PSEUDO-HISTORICAL ARRAY FOR TOTAL CONVECTIVE CLOUDINESS (1D).
! PSCCH      : PSEUDO-HISTORICAL ARRAY FOR CONVECTIVE CLOUD SUMMIT (1D).
! PBCCH      : PSEUDO-HISTORICAL ARRAY FOR CONVECTIVE CLOUD BASE (1D).

!-----------------------------------------------------------------------

! -   ARGUMENTS IMPLICITES.
!     ---------------------

! COMMON/YOMPHY /
! COMMON/YOMCST /
! COMMON/YOMPHY0/

!-----------------------------------------------------------------------

!     Externes.
!     ---------

!     Methode.
!     --------

!     Auteur.
!     -------
!        91-06, M. Deque.

!     Modifications.
!     --------------
!        97-08, PZx variables > Zx, introduce CDLOCK - J.M. Piriou.
!   2003-04-20, M. Bellus: new pseudo-historical arrays PTCCH, PSCCH
!                  and PBCCH (previously contained in PNEBH)
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMPHY0  , ONLY : TPHY0

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TPHY0)       ,INTENT(IN)    :: YDPHY0
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLCL(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLCN(KLON,0:KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNLAB(KLON,KLEV) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNND(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTCCH(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSCCH(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PBCCH(KLON) 

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JLON
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACNEBC',0,ZHOOK_HANDLE)
ASSOCIATE(SPNBCO=>YDPHY0%SPNBCO, SNNBCO=>YDPHY0%SNNBCO, SXNBCO=>YDPHY0%SXNBCO)
!-----------------------------------------------------------------------

!*
!     ------------------------------------------------------------------
!     I - CALCUL DE LA NEBULOSITE TOTALE.

!         COMPUTATION OF THE TOTAL CLOUD COVER.

DO JLON=KIDIA,KFDIA
  PTCCH(JLON)=KNND(JLON)*MIN(SXNBCO,SNNBCO+SPNBCO&
   & *(PFPLCL(JLON,KLEV)+PFPLCN(JLON,KLEV)))  
  PSCCH(JLON)=0.0_JPRB
  PBCCH(JLON)=0.0_JPRB
ENDDO
DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    PSCCH(JLON)=MAX(PSCCH(JLON),(1.0_JPRB-PSCCH(JLON))&
     & *REAL(JLEV*KNLAB(JLON,JLEV),JPRB))  
    PBCCH(JLON)=MAX(PBCCH(JLON),REAL(JLEV*KNLAB(JLON,JLEV),JPRB))
  ENDDO
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACNEBC',1,ZHOOK_HANDLE)
END SUBROUTINE ACNEBC
