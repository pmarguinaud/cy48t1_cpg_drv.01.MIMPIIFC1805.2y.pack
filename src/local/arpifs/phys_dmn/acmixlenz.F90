!OPTIONS XOPT(NOEVAL)
SUBROUTINE ACMIXLENZ ( YDPHY,YDPHY0,KIDIA,KFDIA,KLON,KTDIA,KLEV,LDEML,&
 !-----------------------------------------------------------------------
 ! - INPUT  2D .
 & PAPHI,PAPHIF,&
 ! - INPUT 1D .
 & PBLH,PGZ0,PGZ0H,&
 ! - OUTPUT 2D .
 & PLMU,PLMT)

!**** *ACMIXLENZ * - CALCUL DE LA LONGUEUR DE MELANGE .
!                 COMPUTE THE MIXING LENGTHS.

! -   ARGUMENTS D'ENTREE.
!     -------------------

! - NOM DES PARAMETRES DE DIMENSIONNEMENT DE LA PHYSIQUE.

! KIDIA      : INDICE DE DEPART DES BOUCLES VECTORISEES SUR L'HORIZONT..
! KFDIA      : INDICE DE FIN DES BOUCLES VECTORISEES SUR L'HORIZONTALE.
! KLON       : DIMENSION HORIZONTALE DES TABLEAUX.
! KTDIA     : INDICE DE DEPART DES BOUCLES VERTICALES (1 EN GENERAL)
!              POUR LES CALCULS DE TURBULENCE.
! KTDIA     : START OF THE VERTICAL LOOP FOR THE TURBUL COMPUTATIONS.
! KLEV       : DIMENSION VERTICALE DES TABLEAUX "FULL LEVEL".
! LDEML     : KEY FOR e-MIXING LENGTH PROFILE.

! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
!   CATEGORIE).

! - 2D (0:KLEV) .

! PAPHI      : GEOPOTENTIEL AUX DEMI-NIVEAUX.

! - 2D (1:KLEV) .

! PAPHIF     : GEOPOTENTIEL AUX NIVEAUX DES COUCHES.

! - 1D (PROGNOSTIQUE) .

! PBLH       : HAUTEUR DE LA COUCHE LIMITE.
! PGZ0       : G FOIS LA LONGUEUR DE RUGOSITE COURANTE.
! PGZ0H      : G*LONGUEUR DE RUGOSITE THERMIQUE COURANTE (SI KVCLIV >= 8)

!-----------------------------------------------------------------------

! -   ARGUMENTS DE SORTIE.
!     --------------------

! - 2D (0:KLEV) .
! PLMU        : LONGUEUR DE MELANGE POUR LE VENT.
! PLMT        : LONGUEUR DE MELANGE POUR T ET Q.

!-----------------------------------------------------------------------

! -   ARGUMENTS IMPLICITES.
!     ---------------------

! COMMON/YOMPHY0/
! COMMON/YOMCST /

!-----------------------------------------------------------------------

!     Externes.
!     ---------

!     Methode.
!     --------

!     Auteur.
!     -------
!        2003-11, J.M. Piriou: move lines of source codes from ACCOEFK to the present routine.

!     Modifications.
!     --------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        E. BAZILE     30-06-2004 Correction for USURID=0. and possible use of interactive PBL heigth.
!        J. CEDILNIK   03-03-2006 Reunification of computation of mixig lengths
!                                 for momentum and heat.
!        K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!        F. Vana + I. Bastak Duran 07-Oct-2009 C3 computation in case of TKE
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : RG       
USE YOMPHY   , ONLY : TPHY
USE YOMPHY0  , ONLY : TPHY0

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TPHY)        ,INTENT(IN)    :: YDPHY
TYPE(TPHY0)       ,INTENT(IN)    :: YDPHY0
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA 
LOGICAL           ,INTENT(IN)    :: LDEML
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHI(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBLH(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGZ0(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGZ0H(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLMU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLMT(KLON,0:KLEV) 

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JLON

REAL(KIND=JPRB) :: ZRLTLU,ZGLU,ZGLT,ZZ,ZPFD,ZZT,ZEDIFV,ZPFH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZDIFV(KLON),ZDIFT(KLON)

REAL(KIND=JPRB) :: ZEXPU,ZEXPT

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACMIXLENZ',0,ZHOOK_HANDLE)
ASSOCIATE(VKARMN=>YDPHY0%VKARMN, EDD=>YDPHY0%EDD, C3TKEFREE=>YDPHY0%C3TKEFREE, &
 & ALMAV=>YDPHY0%ALMAV, BEDIFV=>YDPHY0%BEDIFV, A0ML_AU=>YDPHY0%A0ML_AU, &
 & A0ML_AT=>YDPHY0%A0ML_AT, USURID=>YDPHY0%USURID, A0ML_BT=>YDPHY0%A0ML_BT, &
 & A0ML_BU=>YDPHY0%A0ML_BU, &
 & LPRGML=>YDPHY%LPRGML)
!-----------------------------------------------------------------------

! *
! ------------------------------------------------------------------
! I - CALCUL DES PARAMETRES DERIVES ET CONSTANTE DE SECURITE (POUR
! LE CARRE DU CISAILLEMENT DE VENT).

! COMPUTATION OF DERIVED PARAMETERS AND SECURITY CONSTANTS (FOR
! THE SQUARE OF THE WIND SHEAR).

PLMU(:,:)=0.0_JPRB
PLMT(:,:)=0.0_JPRB

     IF (LDEML) THEN 
!      vertical profile of PLMU/PLMT goes from 1 at surf. to C3TKEFREE above PBL 
       ZEDIFV=BEDIFV
       ZRLTLU=C3TKEFREE
     ELSE
!      This is the old version of Mixing length computation
       ZEDIFV=BEDIFV/SQRT(3._JPRB)
       ZRLTLU=1.5_JPRB*EDD
     ENDIF

DO JLON=KIDIA,KFDIA
  ZDIFV(JLON)=1.0_JPRB/(RG*VKARMN*PBLH(JLON))
  ZDIFT(JLON)=ZDIFV(JLON)/SQRT(3._JPRB)
ENDDO

IF(USURID == 0.0_JPRB) THEN
  ZEDIFV=BEDIFV
  ZDIFT(KIDIA:KFDIA)=ZDIFV(KIDIA:KFDIA)
ENDIF

! *
! ------------------------------------------------------------------
! II - PASSAGE EN GEOPOTENTIEL POUR LES LONGUEURS DE MELANGE
! ASYMPTOTIQUES.

! GOING TO GEOPOTENTIAL FOR ASYMPTOTIC MIXING LENGTHS.

ZGLU=RG*ALMAV
ZGLT=ZGLU*ZRLTLU

! *
! ------------------------------------------------------------------
! III - BOUCLE PASSIVE SUR LES NIVEAUX VERTICAUX.

! PASSIVE LOOP ON VERTICAL LEVELS.

DO JLEV=KTDIA,KLEV-1
  !
  ! *
  ! ------------------------------------------------------------------
  ! IV - CALCULS PROPREMENT DITS.
  !
  ! EFFECTIVE CALCULATIONS.
  !
  !
!DEC$ IVDEP
  DO JLON=KIDIA,KFDIA
    !
    ! CALCULS GEOMETRIQUES, LE COEFFICIENT MULTIPLICATEUR ZPFD
    ! (EVENTUELLEMENT DIFFERENT DE 1) POUR LE PROFIL VERTICAL DE LA
    ! DIFFUSION DEPEND DU CHOIX DE UHDIFV ET DE BEDIFV.
    ! GEOMETRIC COMPUTATIONS, THE MULTIPLYING FACTOR ZPFD (POTENTIALLY
    ! DIFFERENT OF 1) FOR THE DIFFUSION'S VERTICAL PROFILE DEPENDS UPON
    ! THE CHOICE OF UHDIFV AND BEDIFV.
    !
    ZZ=VKARMN*(0.5_JPRB*(PAPHIF(JLON,JLEV)+PAPHIF(JLON,JLEV+1))&
     & -PAPHI(JLON,KLEV)+PGZ0(JLON))  
    ZZT=VKARMN*(0.5_JPRB*(PAPHIF(JLON,JLEV)+PAPHIF(JLON,JLEV+1))&
     & -PAPHI(JLON,KLEV)+PGZ0H(JLON))  

     IF (.NOT.LPRGML) THEN
!      This is the old version of Mixing length computation
       ZPFD=MAX(0.0_JPRB,BEDIFV+(1.0_JPRB-BEDIFV)/(1.0_JPRB+(ZDIFV(JLON)*ZZ)**2))
       ZPFH=MAX(0.0_JPRB,ZEDIFV+(1.0_JPRB-ZEDIFV)/(1.0_JPRB+(ZDIFT(JLON)*ZZT)**2))
       PLMU(JLON,JLEV)=ZZ*ZGLU*ZPFD/(ZGLU+ZZ)/RG
       PLMT(JLON,JLEV)=ZZT*ZGLT*ZPFH/(ZGLT+ZZT)/RG

     ELSE
!      new version of mix length computation
       ZEXPU=EXP(-A0ML_AU*SQRT((ZZ/RG)/(VKARMN*PBLH(JLON)))+A0ML_BU)
       ZEXPT=EXP(-A0ML_AT*SQRT((ZZT/RG)/(VKARMN*PBLH(JLON)))+A0ML_BT) 
       PLMU(JLON,JLEV)=(ZZ/(1.0_JPRB+(ZZ/ZGLU)*((1.0_JPRB+ZEXPU)/&
        & (BEDIFV+ZEXPU))))/RG
       PLMT(JLON,JLEV)=(ZZT/(1.0_JPRB+(ZZT/ZGLT)*((1.0_JPRB+ZEXPT)/&
        & (ZEDIFV+ZEXPT))))/RG
     ENDIF
  ENDDO
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACMIXLENZ',1,ZHOOK_HANDLE)
END SUBROUTINE ACMIXLENZ
