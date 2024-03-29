!OPTIONS XOPT(NOEVAL)
SUBROUTINE ACDRAG ( YDML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
 !-----------------------------------------------------------------------
 ! - INPUT  2D .
 & PAPRS,PAPRSF,PDELP,PNBVNO,PRDELP,PU,PV,&
 ! - INPUT  1D .
 & PRCORI,PGETRL,PGWDCS,PVRLAN,PVRLDI,&
 ! - OUTPUT 2D .
 & PSTRDU,PSTRDV,PRAPTRAJ)

!**** *ACDRAG * - EFFET DES ONDES DE GRAVITE OROGRAPHIQUES.

!     Sujet.
!     ------
!     - ROUTINE DE CALCUL ACTIF .
!       CALCUL DES EFFETS SUR LA QUANTITE DE MOUVEMENT DU DEFERLEMENT
!       DES ONDES DE GRAVITE D'ORIGINE OROGRAPHIQUE .
!     - COMPUTATION OF OROGRAPHIC GRAVITY-WAVES' BREAKING EFFECTS ON
!       MOMENTUM .

!**   Interface.
!     ----------
!        *CALL* *ACDRAG*

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
! PNBVNO     : CARRE DE FR. BRUNT-VAISALA DIVISEE PAR G FOIS LA DENSITE.

! - 2D (1:KLEV) .

! PAPRSF     : PRESSION AUX NIVEAUX DES COUCHES.
! PDELP      : EPAISSEUR EN PRESSION DE LA COUCHE.
! PRDELP     : INVERSE DE L'EPAISSEUR EN PRESSION DE LA COUCHE.
! PU         : COMPOSANTE EN X DU VENT.
! PV         : COMPOSANTE EN Y DU VENT.

! - 1D (GEOGRAPHIQUE) .

! PRCORI     : FACTEUR DE CORIOLIS.
! PGETRL     : ECART TYPE DE L'OROGRAPHIE (EN J/KG).
! PVRLAN     : ANISOTROPIE DU RELIEF SOUS MAILLE (1 = ISOTROPE).
! PVRLDI     : ANGLE DU GRADIENT MAXIMUM DU RELIEF AVEC L'AXE DES X.

! - 1D (DIAGNOSTIQUE) .

! PGWDCS     : DENSITE EN SURFACE POUR LE DRAG OROGRAPHIQUE.

!-----------------------------------------------------------------------

! -   ARGUMENTS DE SORTIE.
!     --------------------

! - 2D (0:KLEV) .

! PSTRDU     : FLUX "GRAVITY WAVE DRAG" "U".
! PSTRDV     : FLUX "GRAVITY WAVE DRAG" "V".

!-----------------------------------------------------------------------

! -   ARGUMENTS IMPLICITES.
!     ---------------------

!-----------------------------------------------------------------------

!     Externes.
!     ---------

!     Methode.
!     --------

!     Auteur.
!     -------
!      89-12, J.F. Geleyn.

!     Modifications.
!     --------------
!      04-03, New options LGLT, LNEWD, GWDPROF, GWDVALI - J.F.Geleyn F.Bouyssel
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      2006-03-03 R.Brozkova - improved form drag along the slopes when LNEWD=.TRUE.
!      K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!      O.Riviere 0ct 09 - Save ZRAPP for simpl phys under LGWDSPNL
!      R. El Khatib 20-Oct-2014 Fix a broken vectorization + workaround against Intel bug
!-----------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE PARKIND1  ,ONLY : JPIM     ,JPRB    ,JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : RPI      ,RG

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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNBVNO(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCORI(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGETRL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWDCS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRLAN(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRLDI(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRDU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTRDV(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAPTRAJ(KLON,0:KLEV) 


!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZIPOI(KLON,KLEV),ZNFNO(KLON,0:KLEV),ZPOID(KLON,KLEV)&
 & ,ZRAPP(KLON,0:KLEV),ZU(KLON,KLEV),ZV(KLON,KLEV)&
 & ,ZWGHT(KLON,0:KLEV),ZDS1(KLON),ZDS2(KLON),ZFACS(KLON)&
 & ,ZOBST(KLON),ZPR1(KLON),ZPR2(KLON),ZREF(KLON),ZSTM(KLON)&
 & ,ZSTS(KLON),ZSTUS(KLON),ZSTVS(KLON),ZSUM(KLON)&
 & ,ZSUMD(KLON),ZSUMF(KLON),ZSUMU(KLON),ZSUMV(KLON)&
 & ,ZSUSR(KLON),ZTST1(KLON),ZTST2(KLON),ZUSR(KLON)&
 & ,ZUSUR(KLON),ZVSUR(KLON),ZWGHTK(KLON),ZZAA(KLON),ZZBB(KLON)&
 & ,ZZCR(KLON),ZZCU(KLON),ZZCV(KLON)  
 
INTEGER(KIND=JPIM) :: ITOP, JLEV, JLON

REAL(KIND=JPRB) :: ZA, ZALP1, ZALP2, ZALPHA, ZALTI, ZARGLI, ZAUSR,&
 & ZBETA, ZD, ZDEL1, ZDEL2, ZDIFSR, ZDIFSU, ZDIFSV,&
 & ZDU, ZDV, ZDX, ZEPS1, ZEPS2, ZEPS3, ZEPS4,&
 & ZEPS5, ZEPS6, ZFORCA, ZFORCB, ZFORCL, ZGDT,&
 & ZGDTI, ZLIFT, ZLIFTG, ZPBL, ZPBLF, ZPRRAP, ZRZ,&
 & ZTEMPF(KLON), ZTEST, ZTUNE, ZUSTAR, ZVSTAR, ZX, ZZPR, ZRZR
REAL(KIND=JPRB) :: ZCOS(KLON)
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACDRAG',0,ZHOOK_HANDLE)
ASSOCIATE(HOBST=>YDML_PHY_MF%YRPHY0%HOBST, GWDBC=>YDML_PHY_MF%YRPHY0%GWDBC, GWDCD=>YDML_PHY_MF%YRPHY0%GWDCD, &
 & GWDLT=>YDML_PHY_MF%YRPHY0%GWDLT, GWDVALI=>YDML_PHY_MF%YRPHY0%GWDVALI, GWDPROF=>YDML_PHY_MF%YRPHY0%GWDPROF, &
 & GWDAMP=>YDML_PHY_MF%YRPHY0%GWDAMP, GWDSE=>YDML_PHY_MF%YRPHY0%GWDSE, &
 & TSPHY=>YDML_PHY_MF%YRPHY2%TSPHY, &
 & NPHYREP=>YDML_PHY_MF%YRPHY%NPHYREP, LNEWD=>YDML_PHY_MF%YRPHY%LNEWD, LGLT=>YDML_PHY_MF%YRPHY%LGLT)
!-----------------------------------------------------------------------


PRAPTRAJ(:,:)=0.0
!*
!     ------------------------------------------------------------------
!     I - CALCUL DES PARAMETRES DERIVES, CONSTANTES DE SECURITE (POUR LE
!     FACTEUR D'ANISOTROPIE, LE CARRE DU MODULE DU VENT, LA STABILITE
!     STATIQUE NORMALISEE, UN SEUIL DE PRECISION SUR DES PRESSIONS
!     RECALCULEES, LE MODULE DE VENT ET UN NOMBRE SANS DIMENSION) AINSI
!     QU'UN ARGUMENT LIMITE POUR LA FONCTION SINUS.

!         COMPUTATION OF DERIVED PARAMETERS, SECURITY CONSTANTS (FOR THE
!     ANISOTROPY FACTOR, THE SQUARE OF THE WIND SPEED, THE NORMALISED
!     STATIC STABILITY A PRECISION THRESHOLD FOR RECOMPUTED PRESSURES,
!     THE WIND SPEED AND AN ADIMENSIONALISED PARAMETER) AS WELL AS A LIMIT
!     ARGUMENT FOR THE SINE FUNCTION.

ZGDT=RG*TSPHY
ZGDTI=1.0_JPRB/ZGDT

!     LES CONSTANTES 'ALP' ET 'DEL' POUR LE CALCUL ANISOTROPE RESULTENT
!     D'APPROXIMATIONS ANALYTIQUES D'INTEGRALES ELLIPTIQUES.
!     THE 'ALP' AND 'DEL' VALUES ARE RESULTING FROM ANALYTICAL FITS OF
!     ELLIPTIC INTEGRALS FOR THE ANISOTROPIC CASE.

IF (GWDPROF == 0.0_JPRB) THEN
  ZZPR=0.5_JPRB
ELSE
  ZZPR=((1.0_JPRB+GWDPROF)*LOG(1.0_JPRB+GWDPROF)-GWDPROF)/GWDPROF**2
ENDIF

ZLIFT=GWDLT*TSPHY/(HOBST*ZZPR)
ZLIFTG=GWDLT/(HOBST*ZZPR)

IF (.NOT.LGLT) THEN
  ZLIFTG=0.0_JPRB
ENDIF

ZALP1=0.75_JPRB+0.25_JPRB*LOG(16._JPRB)
ZALP2=3.5_JPRB*RPI-LOG(16._JPRB)-8._JPRB
ZDEL1=2.75_JPRB-0.75_JPRB*LOG(16._JPRB)
ZDEL2=-0.25_JPRB*RPI-LOG(16._JPRB)+4._JPRB

ZEPS1=1.E-06_JPRB
ZEPS2=1.E-04_JPRB
ZEPS3=1.E-09_JPRB
ZEPS4=1.E-05_JPRB
ZEPS5=1.E-02_JPRB

IF (JPRB == JPRD) THEN
  ZEPS6=1.E-12_JPRB
ELSE
  ZEPS6=1.E-6_JPRB
ENDIF

ZARGLI=1.E+06_JPRB

ZBETA=MIN(GWDVALI,1.0_JPRB-ZEPS6)

!*
!     ------------------------------------------------------------------
!     II - CALCULS PRELIMINAIRES DE L'EPAISSEUR EN PRESSION DES COUCHES
!     DIVISEE PAR G*DT ET DE SON INVERSE POUR CHAQUE NIVEAU AINSI QUE
!     MISE A ZERO DES POIDS POUR L'INTEGRATION SUR L'EPAISSEUR DE
!     L'OBSTACLE OROGRAPHIQUE.

!          PRELIMINARY COMPUTATIONS FOR THE LAYER'S PRESSURE THICKNESS
!     DIVIDED BY G*DT AND ITS INVERSE FOR ALL LEVELS AS WELL AS SETTING
!     TO ZERO THE WEIGHTS FOR THE INTEGRAL OVER THE DEPTH OF THE
!     OROGRAPHIC OBSTACLE.

! - TEMPORAIRE(S) 2D (1:KLEV) .

! ZPOID     : DP/(RG*DT) POUR UNE COUCHE ET UN PAS DE TEMPS DONNES.
!            : DP/(RG*DT) FOR A GIVEN LEVEL AND A GIVEN TIME STEP.
! ZIPOI     : INVERSE DE ZPOID.
!            : INVERSE OF ZPOID.
! ZWGHT     : POIDS POUR L'INTEGRALE SUR LA HAUTEUR DE L'OBSTACLE.
!            : WEIGHTS FOR THE INTEGRAL OVER THE DEPTH OF THE OBSTACLE.

DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    ZPOID(JLON,JLEV)=PDELP(JLON,JLEV)*ZGDTI
    ZIPOI(JLON,JLEV)=PRDELP(JLON,JLEV)*ZGDT
    ZWGHT(JLON,JLEV)=0.0_JPRB
  ENDDO
ENDDO

!*
!     ------------------------------------------------------------------
!     III - CALCULS EN SURFACE ET INITIALISATION POUR LA BOUCLE
!     ASCENDANTE.

!           SURFACE COMPUTATIONS AND INITIALIZATION FOR THE ASCENDING
!     LOOP.

! - TEMPORAIRE(S) 2D (0:KLEV) .

! ZNFNO     : FREQUENCE DE BRUNT-VAISALA EFFECTIVE ET NORMEE (N/RHO/G).
!            : NORMALISED EFECTIVE BRUNT-VAISALA FREQUENCY (N/RHO/G).
! ZRAPP     : RAPPORT DU FLUX EN ALTITUDE A CELUI AU SOL.
!            : RATIO OF THE UPPER AIR FLUX TO THE SURFACE ONE.
! ZU        : COMPOSANTE EN X DU VENT EFFECTIF.
!            : X COMPONENT OF THE EFFECTIVE WIND.
! ZV        : COMPOSANTE EN Y DU VENT EFFECTIF.
!            : Y COMPONENT OF THE EFFECTIVE WIND.

!     MOYENNES VERTICALES SUR L'EPAISSEUR DE L'OBSTACLE OROGRAPHIQUE ET
!     CALCUL DES CARACTERISTIQUES (VENT ET STABILITE) EFFECTIVES PAR
!     COMBINAISON LINEAIRE ENTRE CES MOYENNES (CHOISIES EN SURFACE) ET
!     LES CARACTERISTIQUES INITIALES (RETROUVEES A LA HAUTEUR EFFECTIVE
!     DE L'OBSTACLE).
!     VERTICAL AVERAGING OVER THE THICKNESS OF THE OROGRAPHIC OBSTACLE
!     AND CALCULATION OF THE EFFECTIVE VALUES (WIND AND STABILITY)
!     THROUGH LINEAR COMBINATIONS BETWEEN THE AVERAGES (CHOSEN AT THE
!     SURFACE) AND THE ORIGINAL VALUES (OBTAINED AGAIN AT THE EFFECTIVE
!     HEIGHT OF THE OBSTACLE).

! - TEMPORAIRE(S) 1D .

! ZOBST     : HAUTEUR EFFECTIVE DE L'OBSTACLE EN COORDONNEE PRESSION.
!            : EFFECTIVE HEIGHT OF THE OBSTACLE IN PRESSURE COORDINATE.
! ZSUMF     : MOYENNE PONDEREE POUR "N/RHO/G" AU CARRE (<0 ADMIS).
!            : AVERAGED VALUE FOR A SQUARED "N/RHO/G" (<0 ACCEPTED).
! ZSUMU     : MOYENNE PONDEREE POUR LE VENT "U".
!            : AVERAGED VALUE FOR THE "U" WIND.
! ZSUMV     : MOYENNE PONDEREE POUR LE VENT "V".
!            : AVERAGED VALUE FOR THE "V" WIND.

!     SOMMATION PONDEREE.
!     WEIGHTED SUMMATION.

DO JLON=KIDIA,KFDIA
  ZOBST(JLON)=MIN(PAPRS(JLON,KLEV)&
   & -PAPRSF(JLON,KTDIA),MAX(PAPRS(JLON,KLEV)&
   & -PAPRSF(JLON,KLEV),HOBST*PGETRL(JLON)*PGWDCS(JLON)))  
  ZWGHT(JLON,KLEV)=(PAPRS(JLON,KLEV)-PAPRSF(JLON,KLEV))/ZOBST(JLON)
  ZWGHT(JLON,KLEV)=ZWGHT(JLON,KLEV)*(1.0_JPRB-ZBETA)
  ZWGHTK(JLON)=ZWGHT(JLON,KLEV)
  ZSUMU(JLON)=ZWGHT(JLON,KLEV)*PU(JLON,KLEV)
  ZSUMV(JLON)=ZWGHT(JLON,KLEV)*PV(JLON,KLEV)
  ZSUMF(JLON)=ZWGHT(JLON,KLEV)*PNBVNO(JLON,KLEV)
ENDDO

!     ZTEST ET ITOP POUR EVITER DES CALCULS INUTILES.
!     ZTEST AND ITOP TO AVOID UNNECESSARY COMPUTATIONS.

ZTEST=1.0_JPRB
ITOP=KLEV

DO JLEV=KLEV-1,KTDIA,-1

  IF (ZTEST /= 0.0_JPRB) THEN
    ITOP=JLEV

    DO JLON=KIDIA,KFDIA
      ZWGHT(JLON,JLEV)=MAX(0.0_JPRB,(PAPRSF(JLON,JLEV+1)&
       & -MAX(PAPRSF(JLON,JLEV),PAPRS(JLON,KLEV)&
       & -ZOBST(JLON)))/ZOBST(JLON))  
      ZALTI=MIN(1.0_JPRB,(PAPRS(JLON,KLEV)-PAPRS(JLON,JLEV))/ZOBST(JLON))
      ZWGHT(JLON,JLEV)=ZWGHT(JLON,JLEV)&
       & *(1.0_JPRB-ZBETA*((1.0_JPRB-ZALTI)/(1.0_JPRB+GWDPROF*ZALTI)))  
      ZWGHTK(JLON)=ZWGHTK(JLON)+ZWGHT(JLON,JLEV)
      ZSUMU(JLON)=ZSUMU(JLON)+0.5_JPRB*ZWGHT(JLON,JLEV)*(PU(JLON,JLEV)&
       & +PU(JLON,JLEV+1))  
      ZSUMV(JLON)=ZSUMV(JLON)+0.5_JPRB*ZWGHT(JLON,JLEV)*(PV(JLON,JLEV)&
       & +PV(JLON,JLEV+1))  
      ZSUMF(JLON)=ZSUMF(JLON)+ZWGHT(JLON,JLEV)*PNBVNO(JLON,JLEV)
    ENDDO
  ENDIF

  ZTEST=0.0_JPRB
  DO JLON=KIDIA,KFDIA
    ZTEST=ZTEST+ZWGHT(JLON,JLEV)
  ENDDO

ENDDO

DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    ZWGHT(JLON,JLEV)=ZWGHT(JLON,JLEV)*(1.0_JPRB/ZWGHTK(JLON))
  ENDDO
ENDDO

DO JLON=KIDIA,KFDIA
  ZSUMU(JLON)=ZSUMU(JLON)*(1.0_JPRB/ZWGHTK(JLON))
  ZSUMV(JLON)=ZSUMV(JLON)*(1.0_JPRB/ZWGHTK(JLON))
  ZSUMF(JLON)=ZSUMF(JLON)*(1.0_JPRB/ZWGHTK(JLON))
ENDDO

IF(NPHYREP /= 0 .AND. NPHYREP /= -2) ITOP=KTDIA

!     COMBINAISON LINEAIRE.
!     LINEAR COMBINATION.

DO JLON=KIDIA,KFDIA
  ZNFNO(JLON,KLEV)=SQRT(MAX(0.0_JPRB,ZSUMF(JLON)))
  ZPBLF=1.0_JPRB-MIN(1.0_JPRB,(PAPRS(JLON,KLEV)-PAPRSF(JLON,KLEV))/ZOBST(JLON))
  ZU(JLON,KLEV)=PU(JLON,KLEV)+ZPBLF*(ZSUMU(JLON)-PU(JLON,KLEV))
  ZV(JLON,KLEV)=PV(JLON,KLEV)+ZPBLF*(ZSUMV(JLON)-PV(JLON,KLEV))
ENDDO

DO JLEV=KLEV-1,KTDIA,-1
  DO JLON=KIDIA,KFDIA
    ZPBL=1.0_JPRB-MIN(1.0_JPRB,(PAPRS(JLON,KLEV)-PAPRS(JLON,JLEV))/ZOBST(JLON))
    ZNFNO(JLON,JLEV)=SQRT(MAX(0.0_JPRB,PNBVNO(JLON,JLEV)+ZPBL&
     & *(ZSUMF(JLON)-PNBVNO(JLON,JLEV))))  
    ZPBLF=1.0_JPRB-MIN(1.0_JPRB,(PAPRS(JLON,KLEV)-PAPRSF(JLON,JLEV))/ZOBST(JLON))
    ZU(JLON,JLEV)=PU(JLON,JLEV)+ZPBLF*(ZSUMU(JLON)-PU(JLON,JLEV))
    ZV(JLON,JLEV)=PV(JLON,JLEV)+ZPBLF*(ZSUMV(JLON)-PV(JLON,JLEV))
  ENDDO
ENDDO

DO JLON=KIDIA,KFDIA

!     CALCUL D'UN VENT FICTIF CORRESPONDANT DANS LE CAS ISOTROPE AU
!     DRAG ANISOTROPE.
!     COMPUTATION OF A FICTITIOUS WIND THAT WOULD CORRESPOND, IN THE
!     ISOTROPIC CASE, TO THE ANISOTROPIC DRAG.

! - TEMPORAIRE(S) 1D .

! ZUSUR     : COMPOSANTE X DU VENT FICTIF EN SURFACE.
!            : X COMPONENT OF THE FICTITIOUS SURFACE WIND.
! ZVSUR     : COMPOSANTE Y DU VENT FICTIF EN SURFACE.
!            : Y COMPONENT OF THE FICTITIOUS SURFACE WIND.

! Compiler bug
#ifdef __INTEL_COMPILER
  ZCOS(JLON)=COS(2.0_JPRB*PVRLDI(JLON))
ENDDO
DO JLON=KIDIA,KFDIA
  ZUSTAR=ZSUMU(JLON)*ZCOS(JLON)+ZSUMV(JLON)*SIN(2.0_JPRB*PVRLDI(JLON))
  ZVSTAR=ZSUMU(JLON)*SIN(2.0_JPRB*PVRLDI(JLON))-ZSUMV(JLON)*ZCOS(JLON)
#else
  ZUSTAR=ZSUMU(JLON)*COS(2.0_JPRB*PVRLDI(JLON))+ZSUMV(JLON)*SIN(2.0_JPRB*PVRLDI(JLON))
  ZVSTAR=ZSUMU(JLON)*SIN(2.0_JPRB*PVRLDI(JLON))-ZSUMV(JLON)*COS(2.0_JPRB*PVRLDI(JLON))
#endif
  ZA=PVRLAN(JLON)**2+0.5_JPRB*(4._JPRB*(1.0_JPRB-PVRLAN(JLON))*(1.0_JPRB+ZALP1&
   & *PVRLAN(JLON))+PVRLAN(JLON)*LOG(1.0_JPRB/MAX(ZEPS1,PVRLAN(JLON)))&
   & *(1.0_JPRB+ZALP2*PVRLAN(JLON)))/RPI  
  ZD=0.5_JPRB*(4._JPRB*(1.0_JPRB-PVRLAN(JLON))*(1.0_JPRB+ZDEL1*PVRLAN(JLON))-3._JPRB&
   & *PVRLAN(JLON)*LOG(1.0_JPRB/MAX(ZEPS1,PVRLAN(JLON)))*(1.0_JPRB+ZDEL2&
   & *PVRLAN(JLON)))/RPI  
  ZUSUR(JLON)=ZA*ZSUMU(JLON)+ZD*ZUSTAR
  ZVSUR(JLON)=ZA*ZSUMV(JLON)+ZD*ZVSTAR

!     PARTIE "LINEAIRE" DU CALCUL.
!     "LINEAR" PART OF THE COMPUTATION.

! - TEMPORAIRE(S) 1D .

! ZUSR      : MODULE DU VENT FICTIF EN SURFACE.
!            : SPEED OF FICTITIOUS SURFACE WIND.
! ZFACS     : FACTEUR AU SOL POUR LE CALCUL DE ZRAPP.
!            : SURFACE FACTOR FOR THE COMPUTATION OF ZRAPP.
! ZSTUS     : VALEUR "LINEAIRE" DU FLUX EN "U" AU SOL.
!            : "LINEAR" VALUE OF THE "U" SURFACE FLUX.
! ZSTVS     : VALEUR "LINEAIRE" DU FLUX EN "V" AU SOL.
!            : "LINEAR" VALUE OF THE "V" SURFACE FLUX.
! ZSUSR     : PRODUIT SCALAIRE VENT EFFECTIF PAR VENT FICTIF EN SURF..
!            : SURF. SCALAR PRODUCT OF EFFECTIVE AND FICITIOUS WINDS.

  ZSUSR(JLON)=MAX(ZEPS2,ZSUMU(JLON)*ZUSUR(JLON)+ZSUMV(JLON)*ZVSUR(JLON))
  ZUSR(JLON)=SQRT(ZUSUR(JLON)**2+ZVSUR(JLON)**2)
  ZFACS(JLON)=ZNFNO(JLON,KLEV)/ZSUSR(JLON)**3
  ZRAPP(JLON,KLEV)=1.0_JPRB
  ZSTUS(JLON)=GWDSE*PGWDCS(JLON)**2*ZNFNO(JLON,KLEV)*PGETRL(JLON)*ZUSUR(JLON)
  ZSTVS(JLON)=GWDSE*PGWDCS(JLON)**2*ZNFNO(JLON,KLEV)*PGETRL(JLON)*ZVSUR(JLON)


!     PARTIE "AMPLIFICATION RESONANTE" DU CALCUL.
!     "RESONNANT AMPLIFICATION" PART OF THE COMPUTATION.

! - TEMPORAIRE(S) 1D .

! ZTST1     : VALEUR BINAIRE DE TEST POUR LA PREMIERE DEPOSITION.
!            : BINARY TEST VALUE FOR THE FIRST DEPOSITION.
! ZPR1      : PRESSION DU NIVEAU DE PREMIERE DEPOSITION.
!            : PRESSURE AT THE LEVEL OF FIRST DEPOSITION.
! ZSUM      : INTEGRALE POUR LE CALCUL DE LA PHASE DE L'ONDE.
!            : INTEGRAL FOR THE CALCULATION OF THE WAVE'S PHASE.

  ZTST1(JLON)=1.0_JPRB
  ZPR1(JLON)=PAPRSF(JLON,KLEV)
  ZSUM(JLON)=(PAPRS(JLON,KLEV)-PAPRSF(JLON,KLEV))*ZNFNO(JLON,KLEV)&
   & *ZUSR(JLON)/ZSUSR(JLON)  

!     PARTIE "REFLEXION ET DISSIPATION AVAL" DU CALCUL.
!     "REFLEXION AND DOWNSTREAM DISSIPATION" PART OF THE COMPUTATION.

! - TEMPORAIRE(S) 1D .

! ZTST2     : VALEUR BINAIRE DE TEST POUR LA PREMIERE REFLEXION.
!            : BINARY TEST VALUE FOR THE FIRST REFLEXION.
! ZPR2      : PRESSION DU NIVEAU DE PREMIERE REFLEXION.
!            : PRESSURE AT THE LEVEL OF FIRST REFLEXION.
! ZREF      : VALEUR RELATIVE DES FLUX A LA PREMIERE REFLEXION.
!            : RELATIVE FLUX VALUE AT FIRST DEPOSITION.

  ZTST2(JLON)=1.0_JPRB
  ZPR2(JLON)=PAPRSF(JLON,KLEV)
  ZREF(JLON)=1.0_JPRB
  ZTUNE=GWDBC*(ZUSR(JLON)**2/ZSUSR(JLON))**2
  ZTEMPF(JLON)=ZUSR(JLON)/MAX(ZEPS5,ZTUNE*HOBST*PGETRL(JLON)&
   & *PGWDCS(JLON)*ZNFNO(JLON,KLEV))  
  ZZBB(JLON)=MAX(0.0_JPRB,1.0_JPRB-MAX(ZTEMPF(JLON),ZEPS6))
  ZZAA(JLON)=GWDCD*ZZBB(JLON)*ZTUNE*(1.0_JPRB-ZZBB(JLON))

  IF (LNEWD) THEN
    ZSTUS(JLON)=ZSTUS(JLON)*(1.0_JPRB-ZZBB(JLON))/MAX(1.0_JPRB,ZTEMPF(JLON))
    ZSTVS(JLON)=ZSTVS(JLON)*(1.0_JPRB-ZZBB(JLON))/MAX(1.0_JPRB,ZTEMPF(JLON))
    ZZAA(JLON)=GWDCD*ZZBB(JLON)
  ENDIF

ENDDO

DO JLEV=0,KLEV
!DEC$ IVDEP
  DO JLON=KIDIA,KFDIA
    PRAPTRAJ(JLON,JLEV)=GWDSE*PGWDCS(JLON)**2*ZNFNO(JLON,KLEV)*PGETRL(JLON)
  ENDDO
  IF (LNEWD) THEN
    DO JLON=KIDIA,KFDIA
      PRAPTRAJ(JLON,JLEV)=PRAPTRAJ(JLON,JLEV)*(1.0_JPRB-ZZBB(JLON))/MAX(1.0_JPRB,ZTEMPF(JLON))
    ENDDO
  ENDIF
ENDDO
  
!*
!     ------------------------------------------------------------------
!     IV - CALCULS ASCENDANTS DE L'ESTIMATION LINEAIRE DU FLUX RELATIF
!     ET DU REPERAGE DES PREMIERS NIVEAUX DE DEPOSITION ET DE REFLEXION
!     (SI CE DERNIER EXISTE).

!          ASCENDING COMPUTATIONS OF THE LINEAR ESTIMATE OF THE RELATIVE
!     FLUX AND CHECKING OF THE FIRST LEVELS OF DEPOSITION AND OF
!     REFLEXION (IF THE LATTER EXISTS).

!     CALCULS DANS LA BOUCLE VERTICALE.
!     COMPUTATIONS IN THE VERTICAL LOOP.

DO JLEV=KLEV-1,KTDIA,-1
  DO JLON=KIDIA,KFDIA

!     PARTIE "LINEAIRE" DU CALCUL.
!     "LINEAR" PART OF THE COMPUTATION.

    ZAUSR=0.5_JPRB*((ZU(JLON,JLEV)+ZU(JLON,JLEV+1))*ZUSUR(JLON)&
     & +(ZV(JLON,JLEV)+ZV(JLON,JLEV+1))*ZVSUR(JLON))  
    ZPRRAP=ZFACS(JLON)*ZAUSR**3/MAX(ZEPS3,ZNFNO(JLON,JLEV))
    ZRAPP(JLON,JLEV)=MAX(0.0_JPRB,MIN(ZRAPP(JLON,JLEV+1),ZPRRAP))

!     PARTIE "AMPLIFICATION RESONANTE" DU CALCUL.
!     "RESONNANT AMPLIFICATION" PART OF THE COMPUTATION.

    ZTST1(JLON)=ZTST1(JLON)*MAX(0.0_JPRB,SIGN(1.0_JPRB,ZPRRAP-1.0_JPRB))
    ZPR1(JLON)=ZPR1(JLON)+ZTST1(JLON)*(PAPRSF(JLON,JLEV)-ZPR1(JLON))
    ZSUM(JLON)=ZSUM(JLON)+ZTST1(JLON)*(PAPRSF(JLON,JLEV+1)&
     & -PAPRSF(JLON,JLEV))*ZNFNO(JLON,JLEV)*ZUSR(JLON)&
     & /MAX(ZEPS2,ZAUSR)  

!     PARTIE "REFLEXION ET DISSIPATION AVAL" DU CALCUL.
!     "REFLEXION AND DOWNSTREAM DISSIPATION" PART OF THE COMPUTATION.

    ZTST2(JLON)=ZTST2(JLON)*(1.0_JPRB-MAX(0.0_JPRB,SIGN(1.0_JPRB,0.0_JPRB&
     & -ZNFNO(JLON,JLEV))))  
    ZPR2(JLON)=ZPR2(JLON)+ZTST2(JLON)*(PAPRSF(JLON,JLEV)-ZPR2(JLON))
    ZREF(JLON)=ZREF(JLON)+ZTST2(JLON)*(ZRAPP(JLON,JLEV)-ZREF(JLON))
  ENDDO
ENDDO

!     CONDITION A LA LIMITE SUPERIEURE ET PREPARATION DE LA BOUCLE
!     DESCENDANTE.
!     UPPER BOUNDARY CONDITION AND PREPARATION OF THE DESCENDING LOOP.

!DEC$ IVDEP
DO JLON=KIDIA,KFDIA

!     PARTIE "LINEAIRE" DU CALCUL.
!     "LINEAR" PART OF THE COMPUTATION.

  PSTRDU(JLON,KTDIA-1)=0.0_JPRB
  PSTRDV(JLON,KTDIA-1)=0.0_JPRB
  ZRAPP(JLON,KTDIA-1)=0.0_JPRB
  PRAPTRAJ(JLON,KTDIA-1)=0.0_JPRB

!     PARTIE "AMPLIFICATION RESONANTE" DU CALCUL.
!     "RESONNANT AMPLIFICATION" PART OF THE COMPUTATION.

! - TEMPORAIRE(S) 1D .

! ZSTM      : VALEUR MAXIMALE RELATIVE DU FLUX.
!            : MAXIMUM RELATIVE VALUE OF THE FLUX.
! ZDS1      : DERIVEE EN PRESSION POUR LA PARTIE "PIEGEE".
!            : PRESSURE DERIVATIVE FOR THE "TRAPPED" PART.

  ZSTM(JLON)=1.0_JPRB/SQRT(1.0_JPRB+MAX(0.0_JPRB,SIGN(1.0_JPRB,ZPR1(JLON)-ZPR2(JLON)&
   & -ZEPS4))*GWDAMP*(2.0_JPRB*SIN(MIN(ZARGLI,ZSUM(JLON)))&
   & +GWDAMP))  

  IF (LNEWD) THEN
    ZSTM(JLON)=ZSTM(JLON)/(1.0_JPRB-ZZBB(JLON))
    ZREF(JLON)=ZREF(JLON)/(1.0_JPRB-ZZBB(JLON))
  ENDIF

  ZDS1(JLON)=MAX(0.0_JPRB,ZSTM(JLON)-1.0_JPRB)/(PAPRS(JLON,KLEV)-ZPR1(JLON))

!     PARTIE "REFLEXION ET DISSIPATION AVAL" DU CALCUL.
!     "REFLEXION AND DOWNSTREAM DISSIPATION" PART OF THE COMPUTATION.

! - TEMPORAIRE(S) 1D .

! ZSTS      : VALEUR RELATIVE DU FLUX A LA REFLEXION (SI ELLE A LIEU).
!            : RELATIVE VALUE OF THE FLUX AT REFLEXION (IF ANY).
! ZDS2      : DERIVEE EN PRESSION POUR LA PARTIE "AVAL".
!            : PRESSURE DERIVATIVE FOR THE "DOWNSTREAM" PART.

  ZSTS(JLON)=MIN(ZREF(JLON)*(1.0_JPRB-ZTST2(JLON)),ZSTM(JLON))
  ZDS2(JLON)=ZSTS(JLON)/(PAPRS(JLON,KLEV)-ZPR2(JLON))
ENDDO

!*
!     ------------------------------------------------------------------
!     V - BOUCLE DESCENDANTE POUR LA COMBINAISON DES EFFETS ET LES
!     SECURITES NUMERIQUES CONDUISANT AU CALCUL DEFINITIF DES FLUX.

!         DESCENDING LOOP FOR THE COMBINATION OF EFFECTS AND THE
!     NUMERICAL SECURITIES LEADING TO THE FINAL FLUX CALCULATION.

!     FLUX PROVISOIRES PRENANT EN COMPTE LES TROIS EFFETS.
!     PROVISIONAL FLUXES TAKING INTO ACOUNT ALL THREE EFFECTS.

DO JLON=KIDIA,KFDIA
  ZZCR(JLON)=ZRAPP(JLON,KLEV)*ZSTM(JLON)
ENDDO

DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    ZRAPP(JLON,JLEV)=ZDS2(JLON)*MAX(0.0_JPRB,PAPRS(JLON,JLEV)&
     & -ZPR2(JLON))&
     & +MAX(0.0_JPRB,MIN(ZSTM(JLON),ZRAPP(JLON,JLEV)&
     & +ZDS1(JLON)*MAX(0.0_JPRB,PAPRS(JLON,JLEV)&
     & -ZPR1(JLON)))-ZSTS(JLON))  
    ZRZ=(PAPRS(JLON,KLEV)-PAPRS(JLON,JLEV))&
     & /MAX(ZEPS4,ZZBB(JLON)*HOBST*PGETRL(JLON)&
     & *PGWDCS(JLON))

    IF (.NOT.LNEWD) THEN
      ZZCR(JLON)=ZRAPP(JLON,JLEV)
      ZRZR=1.0_JPRB
      ZRZ=MIN(1.0_JPRB,ZRZ)
    ELSE
      ZRZR=SQRT(ZZBB(JLON)/(1.0_JPRB+ZRZ**2*ZZBB(JLON)*(1.0_JPRB-ZZBB(JLON))))
    ENDIF

    ZRAPP(JLON,JLEV)=ZRAPP(JLON,JLEV)+ZZCR(JLON)*(ZZAA(JLON)&
     & *SQRT(MAX(0.0_JPRB,(1.0_JPRB-ZRZ*ZRZR)**3/(1.0_JPRB+ZRZ*GWDPROF&
     & *ZZBB(JLON)))))  
  ENDDO
ENDDO

DO JLON=KIDIA,KFDIA
  ZZCR(JLON)=ZRAPP(JLON,KTDIA-1)
ENDDO

DO JLEV=KTDIA,KLEV
!DEC$ IVDEP
  DO JLON=KIDIA,KFDIA
    ZDIFSR=ZRAPP(JLON,JLEV)-ZZCR(JLON)
    ZZCR(JLON)=ZRAPP(JLON,JLEV)
    ZALTI=MIN(1.0_JPRB,(PAPRS(JLON,KLEV)-PAPRSF(JLON,JLEV))/ZOBST(JLON))
    ZRAPP(JLON,JLEV)=ZRAPP(JLON,JLEV-1)&
     & +ZDIFSR*(1.0_JPRB+ZLIFTG*(1.0_JPRB-ZALTI)/(1.0_JPRB+ZALTI*GWDPROF))&
     & *(1.0_JPRB-ZBETA*((1.0_JPRB-ZALTI)/(1.0_JPRB+GWDPROF*ZALTI)))  
    PSTRDU(JLON,JLEV)=ZRAPP(JLON,JLEV)*ZSTUS(JLON)
    PSTRDV(JLON,JLEV)=ZRAPP(JLON,JLEV)*ZSTVS(JLON)
    PRAPTRAJ(JLON,JLEV)=PRAPTRAJ(JLON,JLEV)*ZRAPP(JLON,JLEV)
  ENDDO
ENDDO

!     MISE EN PSEUDO-IMPLICITE PAR RAPPORT AU VENT DE SURFACE EFFECTIF
!     PROJETE SUR LA DIRECTION DU STRESS.
!     PSEUDO-IMPLICIT ALGORITHM WITH RESPECT TO THE EFFECTIVE SURFACE
!     WIND PROJECTED ON THE STRESS DIRECTION.

! - TEMPORAIRE(S) 1D .

! ZSUMD     : MOYENNE DES TENDANCES DU VENT PROJETE SUR LE STRESS.
!            : AVERAGED WIND TENDENCIES PROJECTED ON THE STRESS' DIREC..

!     SOMMATION POUR UNE MOYENNE DE TENDANCES PROVISOIREMENT NORMEE.
!     SUMMATION FOR A TEMPORARY NORMALISED AVERAGE OF TENDENCIES.

DO JLON=KIDIA,KFDIA
  ZSUMD(JLON)=0.5_JPRB*ZWGHT(JLON,ITOP)*ZIPOI(JLON,ITOP)&
   & *(ZRAPP(JLON,ITOP)-ZRAPP(JLON,ITOP-1))  
ENDDO

DO JLEV=ITOP+1,KLEV
  DO JLON=KIDIA,KFDIA
    ZSUMD(JLON)=ZSUMD(JLON)+0.5_JPRB*(ZWGHT(JLON,JLEV)&
     & +ZWGHT(JLON,JLEV-1))*ZIPOI(JLON,JLEV)&
     & *(ZRAPP(JLON,JLEV)-ZRAPP(JLON,JLEV-1))  
  ENDDO
ENDDO

DO JLON=KIDIA,KFDIA
  ZSUMD(JLON)=ZSUMD(JLON)+0.5_JPRB*ZWGHT(JLON,KLEV)*ZIPOI(JLON,KLEV)&
   & *(ZRAPP(JLON,KLEV)-ZRAPP(JLON,KLEV-1))  
ENDDO

!     VALEUR FINALE DE LA MOYENNE, FACTEUR DE REDUCTION LINEAIRE ET SON
!     APPLICATION.
!     FINAL AVERAGED VALUE, LINEAR REDUCTION FACTOR AND ITS USE.

DO JLON=KIDIA,KFDIA
  ZSUMD(JLON)=ZSUMD(JLON)*SQRT(ZSTUS(JLON)**2+ZSTVS(JLON)**2)
  ZRAPP(JLON,KLEV)=1.0_JPRB/(1.0_JPRB+ZSUMD(JLON)*ZUSR(JLON)/ZSUSR(JLON))
ENDDO

DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    PSTRDU(JLON,JLEV)=PSTRDU(JLON,JLEV)*ZRAPP(JLON,KLEV)
    PSTRDV(JLON,JLEV)=PSTRDV(JLON,JLEV)*ZRAPP(JLON,KLEV)
    PRAPTRAJ(JLON,JLEV)=PRAPTRAJ(JLON,JLEV)*ZRAPP(JLON,KLEV)
  ENDDO
ENDDO

!     SECURITE NUMERIQUE (UTILISATION D'UN ALGORITHME INSPIRE DE CELUI
!     DES ROUTINES DE CONVECTION PROFONDE A FLUX DE MASSE).
!     NUMERICAL SECURITY (USE OF AN ALGORITHM INSPIRED BY THE ONE OF THE
!     MASS-FLUX TYPE DEEP CONVECTION ROUTINES).

DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    ZDU=ZIPOI(JLON,JLEV)*(PSTRDU(JLON,JLEV-1)-PSTRDU(JLON,JLEV))
    ZDV=ZIPOI(JLON,JLEV)*(PSTRDV(JLON,JLEV-1)-PSTRDV(JLON,JLEV))
    ZDX=ZIPOI(JLON,JLEV)*(ZUSUR(JLON)*PSTRDU(JLON,JLEV)&
     & +ZVSUR(JLON)*PSTRDV(JLON,JLEV))  
    ZX=ZUSUR(JLON)*PU(JLON,JLEV)+ZVSUR(JLON)*PV(JLON,JLEV)
    ZALPHA=1.0_JPRB/(1.0_JPRB+ZDX/MAX(ZEPS2,ZX))
    PSTRDU(JLON,JLEV)=PSTRDU(JLON,JLEV-1)-ZPOID(JLON,JLEV)*ZALPHA*ZDU
    PSTRDV(JLON,JLEV)=PSTRDV(JLON,JLEV-1)-ZPOID(JLON,JLEV)*ZALPHA*ZDV
  ENDDO
ENDDO

DO JLON=KIDIA,KFDIA
  ZZCU(JLON)=PSTRDU(JLON,KTDIA-1)
  ZZCV(JLON)=PSTRDV(JLON,KTDIA-1)
ENDDO

DO JLEV=KTDIA,KLEV
!DEC$ IVDEP
  DO JLON=KIDIA,KFDIA
    ZDIFSU=PSTRDU(JLON,JLEV)-ZZCU(JLON)
    ZDIFSV=PSTRDV(JLON,JLEV)-ZZCV(JLON)
    ZZCU(JLON)=PSTRDU(JLON,JLEV)
    ZZCV(JLON)=PSTRDV(JLON,JLEV)
    ZALTI=MIN(1.0_JPRB,(PAPRS(JLON,KLEV)-PAPRSF(JLON,JLEV))/ZOBST(JLON))
    ZFORCL=ZLIFT*PRCORI(JLON)*(1.0_JPRB-ZALTI)/(1.0_JPRB+ZALTI*GWDPROF)
    ZFORCA=ZFORCL/(1.0_JPRB+0.25_JPRB*ZFORCL**2)
    ZFORCB=0.5_JPRB*ZFORCL*ZFORCA
    PSTRDU(JLON,JLEV)=PSTRDU(JLON,JLEV-1)+ZDIFSU+ZPOID(JLON,JLEV)&
     & *(ZFORCB*PU(JLON,JLEV)-ZFORCA*PV(JLON,JLEV))  
    PSTRDV(JLON,JLEV)=PSTRDV(JLON,JLEV-1)+ZDIFSV+ZPOID(JLON,JLEV)&
     & *(ZFORCB*PV(JLON,JLEV)+ZFORCA*PU(JLON,JLEV))  
  ENDDO
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACDRAG',1,ZHOOK_HANDLE)
END SUBROUTINE ACDRAG
