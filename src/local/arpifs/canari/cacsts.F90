!option! -O nomove
SUBROUTINE CACSTS(YDSURF,YDRIP,YDPHY,YDPHY1,KPROMA,KNBPT,KTASK,PT2INC,PH2INC,&
 & PSP_SB,PSP_SG,PSP_RR,PSP_CI,PSP_X2,PSD_VF,PSD_VV,PSD_VX,&
 & PGELAT,PGELAM,PGEMU)

!****---------------------------------------------------------------------------
!****   CACSTS : INITIALISE LES CHAMPS DE SURFACE
!****   ------
!****---------------------------------------------------------------------------

!**  BUT : INITIALISE LES CHAMPS DE SURFACE PROGNOSTIQUES
!**  ---
!**  SEQUENCE D'APPEL :
!**  ----------------
!**        CALL CACSTS(....)

!**  ARGUMENTS D'ENTREE :
!**  ------------------
!**        - IMPLICITE - COMMONS  YOMDPHY YOMCST YOMPHY YOMPHY1
!**                               QACLIM  QACTEX YOMDAG QACVEG
!**        - EXPLICITE - KPROMA : surdimensionnement horizontal
!**                      KNBPT  : nombre reel de points traites
!**                      KTASK  : numero de la tache (multitasking)
!**                      PT2INC : increment d'analyse de T2m
!**                      PH2INC : increment d'analyse de Hu2m
!**                      PSP_SB,PSP_SG,PSP_RR,PSD_VF,PSD_VV,PSD_VX,PSP_CI,PSP_X2    : 
!**                      buffer des champs pdg de l'analyse
!**                      PGELAM, PGELAT, PGEMU : coordonnees geographiques

!**  ARGUMENTS DE SORTIE : 
!**  -------------------

!**  COMMONS UTILISES : - 
!**  ----------------   - 
!**                     - YOMDPHY- DIMENSIONS DES TABLEAUX DE LA PHYSIQUE.
!**                     - YOMCST - CONSTANTES FONDAMENTALES D'ARPEGE.
!**                     - YOMPHY - CLES DE LA PHYSIQUE.
!**                     - YOMDAG - DATE, ECHEANCE ET TYPE DE L'EBAUCHE.
!**                     - YOMPHY1- CONSTANTES PHYSIQUES DU SOL.
!**                     - QACLIM - CLEF D'EXISTENCE DES CHAMPS CLIM.
!**                     - QACTEX - CLES DE L'ANALYSE.
!**                     - QACVEG - CONTROLE DE L'ANALYSE DE SURFACE AVEC ISBA.

!**  EXTERNES : FCTVEG - ACSOLW - TSL
!**  --------

!**  ALGORITHME : - INITIALISE LA TEMPERATURE DE SURFACE.
!**  ----------   - INITIALISE LA TEMPERATURE PROFONDE.
!**               - INITIALISE LE RESERVOIR DE SURFACE.
!**               - INITIALISE LE RESERVOIR PROFOND.
!**               - CORRIGE LA QUANTITE DE NEIGE.

!****  Auteurs :   CB 01/91, BU 05/92, VC 05/93, DG 03/94, PA 09/95, DG 05/96

!****  Modifications:
!       F. Taillefer 09/02 : mise a jour constantes surface selon SST
!       F. Bouyssel 02/04 : Seuil utilisant l'angle zenithal solaire
!       E. Bazile 01/2007 : Parametre pour la correction PSNS et WPI
!       M.Hamrud      01-Jul-2006  Revised surface fields
!       A.Trojakova   27-Jun-2007 bugfixing ZV10M (surface pointers)
!       K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!       F. Taillefer 04/16 : code cleaning/reordering
!***-----------------------------------------------------------------

USE YOMRIP             , ONLY : TRIP
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : RPI      ,RTT      ,RG
USE YOMPHY   , ONLY : TPHY
USE YOMPHY1  , ONLY : TPHY1
USE YOMDAG   , ONLY : NECHGU
USE YOMCLI   , ONLY : YRCLI
USE QACLIM   , ONLY : LCLIM
USE QACTEX   , ONLY : LAEICS   ,RCLIMCA  ,RCLISST  ,LECSST   ,&
 & RSNSA, RSNSB, RWPIA, RWPIB
USE QACVEG   , ONLY : VGAT1    ,VGAT2    ,VGAT3    ,VGAH1    ,&
 & VGAH2    ,VGAH3    ,VGBT1    ,VGBT2    ,VGBT3    ,&
 & VGBH1    ,VGBH2    ,VGBH3    ,VGCT1    ,VGCT2    ,&
 & VGCH1    ,VGCH2    ,SIGT2MP  ,SIGHP1   ,SIGHP2   ,&
 & SIGT2MR  ,SIGH2MR  ,                              &
 & ADWR     ,NLISSEW  ,LSGOBS   ,LIMVEG   ,MINDJ    ,&
 & V10MX    ,SPRECIP  ,SWFC     ,LHUMID   ,LISSEW   ,&
 & SIGT2MO  ,SIGH2MO  ,ANEBUL   ,NNEBUL   ,SNEIGT   ,&
 & NNEIGT   ,SNEIGW   ,NNEIGW   ,SCOEFT   ,SCOEFH   ,&
 & RCLIMTS  ,RCLIMTP  ,RCLIMWS  ,RCLIMWP  ,RCLIMSN  ,&
 & RCLIMN   ,RCLIMV   ,SEVAP    ,SICE     ,SMU0
USE IOGRIDA_MOD, ONLY : RUNDEF

!--------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(TPHY)  ,INTENT(INOUT) :: YDPHY
TYPE(TPHY1) ,INTENT(INOUT) :: YDPHY1
TYPE(TRIP)  ,INTENT(INOUT) :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNBPT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTASK 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT2INC(KNBPT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PH2INC(KNBPT) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSP_SB(KPROMA,YDSURF%YSP_SBD%NLEVS,YDSURF%YSP_SBD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSP_SG(KPROMA,YDSURF%YSP_SGD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSP_RR(KPROMA,YDSURF%YSP_RRD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSP_CI(KPROMA,YDSURF%YSP_CID%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSP_X2(KPROMA,YDSURF%YSP_X2D%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_VF(KPROMA,YDSURF%YSD_VFD%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSD_VV(KPROMA,YDSURF%YSD_VVD%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSD_VX(KPROMA,YDSURF%YSD_VXD%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAT(KNBPT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KNBPT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KNBPT)
!--------------------------------------------------------------------
REAL(KIND=JPRB) :: ZITS(KNBPT), ZITP(KNBPT), ZIWS(KNBPT), ZIWP(KNBPT)
REAL(KIND=JPRB) :: ZIVEG(KPROMA), ZWFC(KPROMA), ZWWILT(KPROMA)
REAL(KIND=JPRB) :: ZWSMX(KPROMA), ZWPMX(KPROMA), ZWPI(KPROMA), ZWSAT(KPROMA)

INTEGER(KIND=JPIM) :: IDJ, IH, JROF

LOGICAL :: LLDHMT

REAL(KIND=JPRB) :: ZCLIM, ZCLIMCA, ZCOEF, ZCWPH, ZCWPT, ZCWSH,&
 & ZCWST, ZDW, ZECHGU, ZEPS, ZEPW, ZEVAP, ZGEL, &
 & ZH2D, ZHEFF, ZHSA, ZHSP, ZLAISRS, ZMSN, ZNEIG, &
 & ZPDN, ZPDS, ZPDT, ZPRECIP, ZSNA, ZSNC, ZT2D, &
 & ZTEFF, ZTINER, ZTPC, ZTSC, ZV10M, ZVEG, ZWP, &
 & ZWPA, ZWPC, ZWPD, ZWPD1, ZWPD2, ZWPDX, ZWPMIN, &
 & ZWPR, ZWSA, ZWSC, ZWSD, ZWSD1, ZWSD2, ZWSMSPI, &
 & ZWSR, ZZH, ZZT, ZDACW, ZMU0, ZPDM, ZPDV
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!--------------------------------------------------------------------

#include "acsolw.intfb.h"
#include "cacdgu.intfb.h"
#include "tsl.intfb.h"

#include "fctveg.func.h"

!--------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CACSTS',0,ZHOOK_HANDLE)
ASSOCIATE(WCRIN=>YDPHY1%WCRIN, NTVGLA=>YDPHY1%NTVGLA, RD1=>YDPHY1%RD1, &
 & SODELX=>YDPHY1%SODELX, GCONV=>YDPHY1%GCONV, WPMX=>YDPHY1%WPMX, &
 & TMERGL=>YDPHY1%TMERGL, RTINER=>YDPHY1%RTINER, WSMX=>YDPHY1%WSMX, &
 & RZHZ0G=>YDPHY1%RZHZ0G, &
 & YSP_SG=>YDSURF%YSP_SG, YSD_VXD=>YDSURF%YSD_VXD, YSD_VV=>YDSURF%YSD_VV, &
 & YSP_RRD=>YDSURF%YSP_RRD, YSP_SGD=>YDSURF%YSP_SGD, YSD_VVD=>YDSURF%YSD_VVD, &
 & YSD_VX=>YDSURF%YSD_VX, YSP_X2D=>YDSURF%YSP_X2D, YSP_RR=>YDSURF%YSP_RR, &
 & YSP_CID=>YDSURF%YSP_CID, YSP_SBD=>YDSURF%YSP_SBD, YSD_VF=>YDSURF%YSD_VF, &
 & YSD_VFD=>YDSURF%YSD_VFD, YSP_SB=>YDSURF%YSP_SB, YSP_CI=>YDSURF%YSP_CI, &
 & YSP_X2=>YDSURF%YSP_X2, &
 & LFGEL=>YDPHY%LFGEL, LSOLV=>YDPHY%LSOLV)
!--------------------------------------------------------------------

!**-----------------------------------------------------------------------
!**  - 1 - Initialisations, diagnostics.
!**        ----------------------------

!*   1.1  Constantes

ZEPS   = 1.E-13_JPRB
ZEPW   = 1.E-3_JPRB
ZECHGU = REAL(NECHGU,JPRB) * 3600._JPRB
!*    Seuil d'evaporation min. pour analyse de W (SEVAP en mm/jour)
IF (LAEICS) ZEVAP  = -SEVAP / 86400._JPRB

!*   1.2  Initialisation des variables intermediaires pour ISBA

IF ( LSOLV ) THEN
  DO JROF = 1 , KNBPT
    ZIVEG(JROF) = ANINT(PSD_VV(JROF,YSD_VV%YIVEG%MP))
    IF (LFGEL) THEN
      ZWPI(JROF) = PSP_SB(JROF,1,YSP_SB%YTL%MP0)
    ELSE
      ZWPI (JROF) = 0.0_JPRB
    ENDIF
  ENDDO
  LLDHMT = .FALSE.
  CALL ACSOLW(YDPHY1,1,KNBPT,KPROMA,PSD_VV(1,YSD_VV%YARG%MP),PSD_VV(1,YSD_VV%YD2%MP),&
   & PSD_VF(1,YSD_VF%YLSM%MP),ZIVEG,PSD_VV(1,YSD_VV%YSAB%MP),LLDHMT,&
   & ZWFC,ZWPMX,ZWSAT,ZWSMX,ZWWILT)
ENDIF

!*   1.3  Diagnostics sur les champs du guess

IF ( LAEICS ) THEN
  CALL CACDGU (KTASK, 'TS GUESS        ',38, PSP_RR(1,YSP_RR%YT%MP0),&
   & KPROMA, KNBPT, 1, 1, PGELAM, PGELAT)
  CALL CACDGU (KTASK, 'TP GUESS        ',41, PSP_SB(1,1,YSP_SB%YT%MP0),&
   & KPROMA, KNBPT, 1, 1, PGELAM, PGELAT)
  CALL CACDGU (KTASK, 'WS GUESS        ',44, PSP_RR(1,YSP_RR%YW%MP0),&
   & KPROMA, KNBPT, 1, 1, PGELAM, PGELAT)
  CALL CACDGU (KTASK, 'WP GUESS        ',47, PSP_SB(1,1,YSP_SB%YQ%MP0),&
   & KPROMA, KNBPT, 1, 1, PGELAM, PGELAT)
ENDIF

!**---------------------------------------------------------------------
!**  - 2 - Calcul des champs analyses.
!**        --------------------------

!ocl nopreex
DO JROF = 1 , KNBPT

! stockage du guess
  ZITS(JROF) = PSP_RR(JROF,YSP_RR%YT%MP0)
  IF ( LSOLV ) THEN
    ZITP(JROF) = PSP_SB(JROF,1,YSP_SB%YT%MP0)
    ZIWS(JROF) = PSP_RR(JROF,YSP_RR%YW%MP0)
    ZIWP(JROF) = PSP_SB(JROF,1,YSP_SB%YQ%MP0)
  ENDIF

  IH = 0
  IDJ= 0
  ZPRECIP= 0.0_JPRB
  ZV10M  = 0.0_JPRB

! analyse de surface ou rappel vers la climatologie, sur terre
  IF (PSD_VF(JROF,YSD_VF%YLSM%MP) > 0.5_JPRB.AND.(LAEICS.OR.RCLIMCA > 0.0_JPRB)) THEN

    ZVEG  = PSD_VF(JROF,YSD_VF%YVEG%MP)
    ZNEIG = MAX(0.0_JPRB,PSP_SG(JROF,YSP_SG%YF%MP0)/(PSP_SG(JROF,YSP_SG%YF%MP0)+WCRIN))

! analyse de surface
    IF ( LAEICS ) THEN

!*   2.O  Initialisations pour l'analyse de surface

      IH = 0
      IDJ= 0
      ZPRECIP= 0.0_JPRB
      ZV10M  = 0.0_JPRB

! conditions locales d'analyse effective des champs de surface
! calcul du temps solaire local utile
      IF ( LSOLV ) THEN
        ZV10M=SQRT(PSP_CI(JROF,YSP_CI%YCI(6)%MP0)**2+PSP_CI(JROF,YSP_CI%YCI(7)%MP0)**2)
        ZPRECIP = PSP_X2(JROF,YSP_X2%YX2(1)%MP0)+ PSP_X2(JROF,YSP_X2%YX2(2)%MP0)&
         & + PSP_X2(JROF,YSP_X2%YX2(3)%MP0)+ PSP_X2(JROF,YSP_X2%YX2(4)%MP0)  
        CALL TSL(YDRIP,IH,IDJ,ZMU0,PGELAT(JROF),PGELAM(JROF),PGEMU(JROF))

        ! Temporary fix
        IF (IH == 0) IH = 24

        ZDACW = MIN(1.0_JPRB,MAX(0.0_JPRB,ABS(REAL(NINT(ZIVEG(JROF))-NTVGLA,JPRB))))&
         & * MIN(1.0_JPRB,MAX(0.0_JPRB,REAL(IH,JPRB)))&
         & * MIN(1.0_JPRB,MAX(0.0_JPRB,REAL(IDJ,JPRB)/REAL(MINDJ,JPRB)))&
         & * MIN(1.0_JPRB,MAX(0.0_JPRB,1.0_JPRB-ZV10M/(V10MX+ZEPS)))&
         & * MIN(1.0_JPRB,MAX(0.0_JPRB,1.0_JPRB-ZPRECIP/(SPRECIP+ZEPS)))&
         & * MIN(1.0_JPRB,MAX(0.0_JPRB,1.0_JPRB-ZWPI(JROF)/SICE))  

! coefficients : dependance par rapport a l'angle zenithal solaire
        IF ( SMU0 > ZEPS ) THEN
          ZPDM=0.5_JPRB*(1.0_JPRB+TANH(SMU0*(ZMU0-0.5_JPRB)))
        ELSE
          ZPDM=1.0_JPRB
        ENDIF
        ZDACW = ZDACW * ZPDM

! coefficients : dependance par rapport a l'evaporation de surface
        IF ( SEVAP > ZEPS ) THEN
          ZPDV=MIN(1.0_JPRB,MAX(0.0_JPRB,PSP_X2(JROF,YSP_X2%YX2(6)%MP0)/ZEVAP))
        ELSE
          ZPDV=1.0_JPRB
        ENDIF
        ZDACW = ZDACW * ZPDV

! coefficients : dependance par rapport a la nebulosite
        IF ( ANEBUL > ZEPS ) THEN
          ZPDN=1.0_JPRB-ANEBUL*(PSP_X2(JROF,YSP_X2%YX2(5)%MP0)/ZECHGU)**NNEBUL
        ELSE
          ZPDN=1.0_JPRB
        ENDIF
        ZDACW = ZDACW * ZPDN
      ELSE
        ZDACW = 0.0_JPRB
      ENDIF

! increments de temperature et d'humidite relative a 2m utiles
      ZT2D = PT2INC(JROF)
      ZH2D = PH2INC(JROF)

!*   2.1  Analyse de temperature

! report de l'increment de temperature a 2m sur Ts et Tp avec amortissement
      IF ( LSOLV ) THEN
        ZTINER = SODELX(1)/SODELX(0)
        IF (NNEIGT <= 0.OR. ZNEIG < ZEPS) THEN
          ZPDT= 1.0_JPRB
        ELSEIF (SNEIGT < ZEPS) THEN
          ZPDT= 0.0_JPRB
        ELSE
          ZPDT= (1.0_JPRB- MIN(ZNEIG,SNEIGT)/SNEIGT)**NNEIGT
        ENDIF
      ELSE
        ZTINER = RTINER
        ZPDT = 1.0_JPRB
      ENDIF
      PSP_RR(JROF,YSP_RR%YT%MP0) = PSP_RR(JROF,YSP_RR%YT%MP0) + ZT2D*ZPDT
      PSP_SB(JROF,1,YSP_SB%YT%MP0) = PSP_SB(JROF,1,YSP_SB%YT%MP0) + ZT2D*ZPDT/ZTINER

!*   2.2  Analyse d'humidite par interpolation optimale pour ISBA

      IF ( LSOLV ) THEN

! coefficients : dependance principale par rapport a la vegetation
        ZLAISRS = PSD_VV(JROF,YSD_VV%YLAI%MP)/MAX(1.0_JPRB,PSD_VV(JROF,YSD_VV%YRSMIN%MP))
        ZCWST = VGST(IH,ZVEG)
        ZCWSH = VGSH(IH,ZVEG)
        ZCWPT = VGPT1(IH,ZVEG) + ZLAISRS*VGPT2(IH,ZVEG)
        ZCWPH = VGPH1(IH,ZVEG) + ZLAISRS*VGPH2(IH,ZVEG)

! coefficients : dependance par rapport a la texture
        ZDW = (ZWFC(JROF)-ZWWILT(JROF))/ADWR

! coefficients : dependance par rapport aux erreurs d'observation
        IF ( LSGOBS ) THEN
          ZZT = G(SIGH2MO/SIGH2MP(IH,ZVEG),SIGT2MO/SIGT2MP(IH))&
           & / G(SIGH2MR/SIGH2MP(IH,ZVEG),SIGT2MR/SIGT2MP(IH))  
          ZZH = G(SIGT2MO/SIGT2MP(IH),SIGH2MO/SIGH2MP(IH,ZVEG))&
           & / G(SIGT2MR/SIGT2MP(IH),SIGH2MR/SIGH2MP(IH,ZVEG))  
        ELSE
          ZZT=1.0_JPRB
          ZZH=1.0_JPRB
        ENDIF

! coefficients : dependance par rapport a la couverture neigeuse
        IF (NNEIGW <= 0.OR. ZNEIG < ZEPS) THEN
          ZPDS= 1.0_JPRB
        ELSEIF (SNEIGW < ZEPS) THEN
          ZPDS= 0.0_JPRB
        ELSE
          ZPDS= (1.0_JPRB- MIN(ZNEIG,SNEIGW)/SNEIGW)**NNEIGW
        ENDIF
        ZDACW = ZDACW * ZPDS

! calcul des increments bruts pour ws=Ws/ds/ro, wp=Wp/dp/ro
! coefficients finaux
        ZCWST = ZCWST * ZDW * ZZT * ZDACW
        ZCWSH = ZCWSH * ZDW * ZZH * ZDACW
        ZCWPT = ZCWPT * ZDW * ZZT * ZDACW
        ZCWPH = ZCWPH * ZDW * ZZH * ZDACW
! limitation eventuelle des increments de T2m et H2m
! limitation de la valeur absolue des increments
        IF (SIGT2MO < 0.0_JPRB) ZT2D=MAX(SIGT2MO,MIN(-SIGT2MO,ZT2D))
        IF (SIGH2MO < 0.0_JPRB) ZH2D=MAX(SIGH2MO,MIN(-SIGH2MO,ZH2D))
! retrait du biais moyen
! soustraction du biais moyen si SCOEF(T/H) <> 1
        PSP_CI(JROF,YSP_CI%YCI(12)%MP0)=&
         & PSP_CI(JROF,YSP_CI%YCI(12)%MP0)*(1.0_JPRB-SCOEFT)+ZT2D*SCOEFT  
        PSP_CI(JROF,YSP_CI%YCI(13)%MP0)=&
         & PSP_CI(JROF,YSP_CI%YCI(13)%MP0)*(1.0_JPRB-SCOEFH)+ZH2D*SCOEFH  
        ZTEFF = ZT2D - PSP_CI(JROF,YSP_CI%YCI(12)%MP0)
        ZHEFF = ZH2D - PSP_CI(JROF,YSP_CI%YCI(13)%MP0)
! si le biais courant est inferieur au biais moyen on le met a zero
!                IF (ABS(ZT2D).LT.ABS(PSP_CI(JROF,YSP_CI%YCI(12)%MP0)) ZTEFF = 0.
!                IF (ABS(ZH2D).LT.ABS(PSP_CI(JROF,YSP_CI%YCI(13)%MP0)) ZHEFF = 0.
! si le biais courant est inferieur au biais effectif on le garde
        IF ( (SCOEFT /= 0.0_JPRB) .OR. (SCOEFH /= 0.0_JPRB) ) THEN
        IF (ABS(ZT2D) < ABS(ZTEFF)) ZTEFF = ZT2D
        IF (ABS(ZH2D) < ABS(ZHEFF)) ZHEFF = ZH2D
        ZT2D = ZTEFF
        ZH2D = ZHEFF
        ENDIF

! increments bruts
        ZWSD = ZCWST*ZT2D + ZCWSH*ZH2D
        ZWPD = ZCWPT*ZT2D + ZCWPH*ZH2D

! limitations sur les corrections
! pas d'analyse de ws si pas d'evaporation sur sol nu
        IF (PSP_X2(JROF,YSP_X2%YX2(6)%MP0)-PSP_X2(JROF,YSP_X2%YX2(7)%MP0) >= 0.0_JPRB)THEN
          ZWSD = 0.0_JPRB
        ENDIF
! analyse de wp limitee pour assurer veg*wwilt <= wp <= SWFC*wfc
        IF ( LIMVEG ) THEN
          ZWPR = ZIWP(JROF)/(PSD_VV(JROF,YSD_VV%YD2%MP)*GCONV)
          IF ( ZWPR > ZWFC(JROF)*SWFC ) THEN
            IF ( LHUMID ) THEN
              ZWPD = MIN(0.0_JPRB,ZWPD)
            ELSE
              ZWPD = 0.0_JPRB
            ENDIF
          ELSEIF ( ZWPR < ZWWILT(JROF)*ZVEG ) THEN
            IF ( LHUMID ) THEN
              ZWPD = MAX(0.0_JPRB,ZWPD)
            ELSE
              ZWPD = 0.0_JPRB
            ENDIF
          ELSE
            ZWPD1 = ZWWILT(JROF)*ZVEG -ZWPR
            ZWPD2 = ZWFC(JROF)*SWFC   -ZWPR
            ZWPD = MAX(ZWPD1,MIN(ZWPD2,ZWPD))
          ENDIF
! analyse de ws limitee pour assurer veg*wwilt <= ws <= SWFC*wfc
          ZWSR = ZIWS(JROF)/(RD1*GCONV)
          IF ( ZWSR > ZWFC(JROF)*SWFC ) THEN
            IF ( LHUMID ) THEN
              ZWSD = MIN(0.0_JPRB,ZWSD)
            ELSE
              ZWSD = 0.0_JPRB
            ENDIF
          ELSEIF ( ZWSR < ZWWILT(JROF)*ZVEG) THEN
            IF (LHUMID) THEN
              ZWSD = MAX(0.0_JPRB,ZWSD)
            ELSE
              ZWSD = 0.0_JPRB
            ENDIF
          ELSE
            ZWSD1 = ZWWILT(JROF)*ZVEG -ZWSR
            ZWSD2 = ZWFC(JROF)*SWFC   -ZWSR
            ZWSD = MAX(ZWSD1,MIN(ZWSD2,ZWSD))
          ENDIF
        ENDIF

! lissage des increments d'analyse de wp
        IF ( LISSEW ) THEN
          ZWPDX = ZWPD
          IF ( NLISSEW >= 3 ) THEN
            ZWPD =.25_JPRB*(PSP_CI(JROF,YSP_CI%YCI(11)%MP0)+PSP_CI(JROF,YSP_CI%YCI(10)%MP0)+&
             & PSP_CI(JROF,YSP_CI%YCI(9)%MP0)+ZWPDX)  
          ELSE
            ZWPD = 0.0_JPRB
          ENDIF
          IF ( NLISSEW >= 2 ) THEN
            PSP_CI(JROF,YSP_CI%YCI(11)%MP0)=PSP_CI(JROF,YSP_CI%YCI(10)%MP0)
          ENDIF
          IF ( NLISSEW >= 1 ) THEN
            PSP_CI(JROF,YSP_CI%YCI(10)%MP0)=PSP_CI(JROF,YSP_CI%YCI(9)%MP0)
          ENDIF
          PSP_CI(JROF,YSP_CI%YCI(9)%MP0)=ZWPDX
        ENDIF

! report des increments sur Ws, Wp
        ZWSA = PSP_RR(JROF,YSP_RR%YW%MP0) + ZWSD*RD1*GCONV
        ZWPA = PSP_SB(JROF,1,YSP_SB%YQ%MP0) + ZWPD*PSD_VV(JROF,YSD_VV%YD2%MP)*GCONV
        PSP_RR(JROF,YSP_RR%YW%MP0) = MAX(0.0_JPRB,MIN(ZWSMX(JROF),ZWSA))
! contenu en eau total minimum
        ZWPMIN = MAX(PSP_RR(JROF,YSP_RR%YW%MP0),ZEPW*PSD_VV(JROF,YSD_VV%YD2%MP)*GCONV)
        PSP_SB(JROF,1,YSP_SB%YQ%MP0) = MAX(ZWPMIN,MIN(ZWPMX(JROF),ZWPA))

!*   2.3  Analyse d'humidite elementaire

      ELSE

        ZWSMSPI = WSMX/RPI

! report de l'increment d'humidite a 2m sur l'humidite de surface Hs
        ZHSP = 0.5_JPRB * (1.0_JPRB - COS(PSP_RR(JROF,YSP_RR%YW%MP0)/ZWSMSPI))
        ZHSA = MAX(0.0_JPRB,MIN(1.0_JPRB,ZHSP+ZH2D))
        PSP_RR(JROF,YSP_RR%YW%MP0) = ZWSMSPI * ACOS(1.0_JPRB-2.0_JPRB*ZHSA)

! report de l'increment de Ws sur Wp
        ZWSD = PSP_RR(JROF,YSP_RR%YW%MP0)-ZIWS(JROF)
        PSP_SB(JROF,1,YSP_SB%YQ%MP0) = MAX(0.0_JPRB,MIN(WPMX,PSP_SB(JROF,1,YSP_SB%YQ%MP0)+ZWSD))
      ENDIF

    ENDIF

!*   2.4  Rappel vers la climatologie

! mise a jour des champs climatologiques
    IF ( .NOT. LCLIM ) THEN
      ZTSC = ZITS(JROF)
      ZTPC = ZITP(JROF)
      ZWSC = ZIWS(JROF)
      ZWPC = ZIWP(JROF)
      ZSNC = PSP_SG(JROF,YSP_SG%YF%MP0)
    ELSEIF ( LSOLV ) THEN
      ZTSC = PSD_VX(JROF,YSD_VX%YTSC%MP)
      ZTPC = PSD_VX(JROF,YSD_VX%YTPC%MP)
      ZWSC = PSD_VX(JROF,YSD_VX%YPWS%MP) * ZWSMX(JROF)
      ZWPC = PSD_VX(JROF,YSD_VX%YPWP%MP) * ZWPMX(JROF)
      ZSNC = PSD_VX(JROF,YSD_VX%YSNO%MP)
    ELSE
      ZTSC = PSD_VX(JROF,YSD_VX%YTSC%MP)
      ZTPC = PSD_VX(JROF,YSD_VX%YTPC%MP)
      ZWSC = PSD_VX(JROF,YSD_VX%YPWS%MP) * WSMX
      ZWPC = PSD_VX(JROF,YSD_VX%YPWP%MP) * WPMX
      ZSNC = PSD_VX(JROF,YSD_VX%YSNO%MP)
    ENDIF

    ZCLIM = RCLIMCA /(1.0_JPRB+RCLIMN*ZNEIG)
! Rappel de Ts
    ZCLIMCA = RCLIMTS * ZCLIM
    PSP_RR(JROF,YSP_RR%YT%MP0) = (1.0_JPRB-ZCLIMCA)*PSP_RR(JROF,YSP_RR%YT%MP0) + ZCLIMCA*ZTSC
! Rappel de Tp
    ZCLIMCA = RCLIMTP * RCLIMCA
    PSP_SB(JROF,1,YSP_SB%YT%MP0) = (1.0_JPRB-ZCLIMCA)*PSP_SB(JROF,1,YSP_SB%YT%MP0)+ ZCLIMCA*ZTPC
! Rappel de Ws
    ZCLIMCA = RCLIMWS * ZCLIM
    ZCLIMCA = ZCLIMCA* ZVEG + MIN(1.0_JPRB,RCLIMV*ZCLIMCA)* (1.0_JPRB-ZVEG)
    PSP_RR(JROF,YSP_RR%YW%MP0) = (1.0_JPRB-ZCLIMCA)*PSP_RR(JROF,YSP_RR%YW%MP0)+ ZCLIMCA*ZWSC
! Rappel de Wp
    ZCLIMCA = RCLIMWP * ZCLIM
    ZCLIMCA = ZCLIMCA* ZVEG + MIN(1.0_JPRB,RCLIMV*ZCLIMCA)* (1.0_JPRB-ZVEG)
    IF ( LSOLV .AND. LFGEL ) THEN
      ZWP = PSP_SB(JROF,1,YSP_SB%YQ%MP0)
      ZGEL = ZWPI(JROF) / MAX(ZWP+ZWPI(JROF),ZEPS)
      ZWPC = ZWPC * (1.0_JPRB - MAX(0.0_JPRB,MIN(1.0_JPRB,ZGEL)))
      ZWPC=MAX(ZWPMIN,ZWPC)
    ENDIF
    PSP_SB(JROF,1,YSP_SB%YQ%MP0) = (1.0_JPRB-ZCLIMCA)*PSP_SB(JROF,1,YSP_SB%YQ%MP0)+ ZCLIMCA*ZWPC

! rappel de Sn avec correction eventuelle pour fonte
    ZCLIMCA = RCLIMSN * RCLIMCA
    ZSNA = (1.0_JPRB-ZCLIMCA)*PSP_SG(JROF,YSP_SG%YF%MP0)+ZCLIMCA*ZSNC
    ZCOEF= RSNSA/21600._JPRB * ZECHGU
    ZMSN = MAX (0.0_JPRB, ZCOEF*(PSP_CI(JROF,YSP_CI%YCI(4)%MP0)-RTT))**RSNSB
    PSP_SG(JROF,YSP_SG%YF%MP0) = MAX (ZSNA-ZMSN ,0.0_JPRB)
    IF (LFGEL) THEN
      ZCOEF= RWPIA/21600._JPRB * ZECHGU
      ZMSN = MAX (0.0_JPRB, ZCOEF*(PSP_CI(JROF,YSP_CI%YCI(4)%MP0)-RTT))**RWPIB
      PSP_SB(JROF,1,YSP_SB%YTL%MP0)=MAX (ZWPI(JROF)-ZMSN ,0.0_JPRB)
      PSP_SB(JROF,1,YSP_SB%YQ%MP0)=PSP_SB(JROF,1,YSP_SB%YQ%MP0)-PSP_SB(JROF,1,YSP_SB%YTL%MP0)+ZWPI(JROF)
    ENDIF

!*   2.5  Calculations over sea

  ELSEIF ( PSD_VF(JROF,YSD_VF%YLSM%MP) <= 0.5_JPRB  .AND. PSP_CI(JROF,YSP_CI%YCI(8)%MP0) /= RUNDEF ) THEN

!*  Climatological SST over sea
!     (Hirlam : replaced by ECMWF SST / combine Ts and SST only with old_surface (LAEICS=T))
    IF ( (RCLISST /= 0.0_JPRB .AND. .NOT.LECSST) .OR.&
      &  (RCLISST /= 0.0_JPRB .AND. LAEICS .AND. LECSST) ) THEN
      PSP_RR(JROF,YSP_RR%YT%MP0) = (1.0_JPRB-RCLISST)*PSP_RR(JROF,YSP_RR%YT%MP0)&
       &                          + RCLISST*PSP_CI(JROF,YSP_CI%YCI(8)%MP0)  
      IF (LSOLV) PSP_SB(JROF,1,YSP_SB%YT%MP0) = PSP_RR(JROF,YSP_RR%YT%MP0)
    ENDIF

!*  Mise a jour des constantes de surface, en fonction de la banquise

    IF ( PSP_RR(JROF,YSP_RR%YT%MP0) <= TMERGL ) THEN
      PSD_VF(JROF,YSD_VF%YALBF%MP)  = YRCLI%SALBB
      PSD_VF(JROF,YSD_VF%YEMISF%MP) = YRCLI%SEMIB
      PSD_VF(JROF,YSD_VF%YZ0F%MP )  = YRCLI%SZZ0B*RG
      PSD_VV(JROF,YSD_VV%YZ0H%MP )  = RZHZ0G*YRCLI%SZZ0B*RG
    ELSE
      PSD_VF(JROF,YSD_VF%YALBF%MP)  = YRCLI%SALBM
      PSD_VF(JROF,YSD_VF%YEMISF%MP) = YRCLI%SEMIM
    ENDIF

  ENDIF

!*   2.7  Calcul des increments d'analyse

  ZITS(JROF) = PSP_RR(JROF,YSP_RR%YT%MP0) - ZITS(JROF)
  IF ( LSOLV ) THEN
    ZITP(JROF) = PSP_SB(JROF,1,YSP_SB%YT%MP0) - ZITP(JROF)
    ZIWS(JROF) = PSP_RR(JROF,YSP_RR%YW%MP0) - ZIWS(JROF)
    ZIWP(JROF) = PSP_SB(JROF,1,YSP_SB%YQ%MP0) - ZIWP(JROF)
  ENDIF

ENDDO

!**-----------------------------------------------------------------------
!**  - 3 - Impression statistiques champs analyses.
!**        ---------------------------------------

IF ( LAEICS ) THEN

!*   3.1  Champs bruts

  CALL CACDGU (KTASK, 'TS ANALYSE      ',39, PSP_RR(1,YSP_RR%YT%MP0),&
   & KPROMA, KNBPT, 1, 1, PGELAM, PGELAT)
  CALL CACDGU (KTASK, 'TP ANALYSE      ',42, PSP_SB(1,1,YSP_SB%YT%MP0),&
   & KPROMA, KNBPT, 1, 1, PGELAM, PGELAT)
  CALL CACDGU (KTASK, 'WS ANALYSE      ',45, PSP_RR(1,YSP_RR%YW%MP0),&
   & KPROMA, KNBPT, 1, 1, PGELAM, PGELAT)
  CALL CACDGU (KTASK, 'WP ANALYSE      ',48, PSP_SB(1,1,YSP_SB%YQ%MP0),&
   & KPROMA, KNBPT, 1, 1, PGELAM, PGELAT)

!*   3.2  Increments

  CALL CACDGU (KTASK, 'TS ANA-GUE      ',40, ZITS,&
   & KPROMA, KNBPT, 1, 1, PGELAM, PGELAT)
  CALL CACDGU (KTASK, 'TP ANA-GUE      ',43, ZITP,&
   & KPROMA, KNBPT, 1, 1, PGELAM, PGELAT)
  CALL CACDGU (KTASK, 'WS ANA-GUE      ',46, ZIWS,&
   & KPROMA, KNBPT, 1, 1, PGELAM, PGELAT)
  CALL CACDGU (KTASK, 'WP ANA-GUE      ',49, ZIWP,&
   & KPROMA, KNBPT, 1, 1, PGELAM, PGELAT)

ENDIF

!**---------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CACSTS',1,ZHOOK_HANDLE)
END SUBROUTINE CACSTS
