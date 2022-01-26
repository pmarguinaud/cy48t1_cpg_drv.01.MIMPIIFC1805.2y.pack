!OPTIONS XOPT(NOEVAL)
SUBROUTINE ACMIXELEN(YGFL,YDPHY,YDPHY0,&
 & KIDIA,KFDIA,KLON,KTDIA,KLEV,KSTEP,KMLTYPE,&
 !-----------------------------------------------------------------------
 ! - INPUT  2D .
 & PAPHI,PAPHIF,PAPRS,PAPRSF,PT,PQV,PQLI,PQICE,PQRAIN,PQSNOW,PTKE,&
 & PR,PU,PV,PALPH,PLNPR,&
 & PMN2PP,PFMGST,PFPLSH,PFPLCH, &
 ! - INPUT 1D .
 & PGZ0,PTS,PCDN,PBLH,&
 ! - INPUT/OUTPUT 2D .
 & PLMU,PLMT,PLMLF,&
 ! - OUTPUT
 & PLML,PLMLTILD,PRRCOR,&
 ! LOGICAL
 & LDMAF)

!**** *ACMIXELEN * - COMPUTATION OF MIXING LENGTHS FOR TKE SCHEMES.

! -   ARGUMENTS D'ENTREE.
!     -------------------

! - NOM DES PARAMETRES DE DIMENSIONNEMENT DE LA PHYSIQUE.

! KIDIA     : INDICE DE DEPART DES BOUCLES VECTORISEES SUR L'HORIZONT..
! KFDIA     : INDICE DE FIN DES BOUCLES VECTORISEES SUR L'HORIZONTALE.
! KLON      : DIMENSION HORIZONTALE DES TABLEAUX.
! KTDIA     : INDICE DE DEPART DES BOUCLES VERTICALES (1 EN GENERAL)
!             POUR LES CALCULS DE TURBULENCE.
! KTDIA     : START OF THE VERTICAL LOOP FOR THE TURBUL COMPUTATIONS.
! KLEV      : DIMENSION VERTICALE DES TABLEAUX "FULL LEVEL".
! KSTEP     : NUMERO DU PAS DE TEMPS.
! KMLTYPE   : TYPE OF MIXING LENGTH:
!               0 - Geleyn-Cedilnik
!               1 - BL89 (Bougeault-Lacarrere 1989)
!               2 - BL89 in unstable case, GB08 (Grisogono and Belusic 2008)
!                   in stable case, the latter is generalized Deardorff (1980)

! - NOM DES VARIABLES DE LA PHYSIQUE (PAR ORDRE ALPHABETIQUE DANS CHAQUE
!   CATEGORIE).

! - 2D (0:KLEV) .

! PAPHI      : GEOPOTENTIEL AUX DEMI-NIVEAUX.
! PAPRS      : PRESSION AUX DEMI-NIVEAUX.
! PFMGST     : STABILITY FUNCTION F_m FOR MOIST GUSTINESS CORRECTION.
! PFPLSH     : HISTORICAL FIELD FOR STRATIFORM PRECIPITATION.
! PFPLCH     : HISTORICAL FIELD FOR CONVECTIVE PRECIPITATION.

! - 2D (1:KLEV) .

! PAPHIF     : GEOPOTENTIEL AUX NIVEAUX DES COUCHES.
! PAPRSF     : PRESSION AUX NIVEAUX DES COUCHES.
! PT         : TEMPERATURE.
! PQV        : PROGNOSTIC VATER VAPOUR.
! PQLI       : PROGNOSTIC SUSPENDED CLOUD WATER.
! PQICE      : PROGNOSTIC SUSPENDED CLOUD ICE.
! PQRAIN     : PROGNOSTIC LIQUID PRECIPITATIONS.
! PQSNOW     : PROGNOSTIC FROZEN PRECIPITATIONS.
! PTKE       : PROGNOSTIC FULL LEVEL TKE.
! PR         : CONSTANTE DES GAZ POUR L'AIR.
! PU         : COMPOSANTE EN X DU VENT.
! PV         : COMPOSANTE EN Y DU VENT.
! PALPH      : LOG(PAPRS(JLEV)/PAPRSF(JLEV)) (FOR HYDROSTATICS).
! PLNPR      : LOG(PAPRS(JLEV)/PAPRS(JLEV-1)) (FOR HYDROSTATICS).


! - 1D (PROGNOSTIQUE) .

! PGZ0       : G FOIS LA LONGUEUR DE RUGOSITE COURANTE.
! PTS        : TEMPERATURE DE SURFACE.
! PCDN       : NEUTRAL SURFACE DRAG COEFFICIENT.
! PBLH       : PLANETARY BOUNDARY LAYER HEIGHT.

!-----------------------------------------------------------------------

! -   ARGUMENTS D'ENTREE/SORTIE.
!     --------------------------

! - 2D (0:KLEV) .
! PLMU        : LONGUEUR DE MELANGE POUR LE VENT.
! PLMT        : LONGUEUR DE MELANGE POUR T ET Q.

! PMN2PP      : MOIST BRUNT-VAISALA FREQUENCY

! - 2D (1:KLEV) .
! PLMLF      : PROGNOSTIC MIX. LENGTH ON FULL LEVELS

! -   ARGUMENTS D'SORTIE.
!     -------------------

! PLML       : TKE type mixing length (on half levels)
! PLMLTILD   : 'STATIONARY' TKE type mixing length (on half levels) 
! PRRCOR     : MOIST GUSTINESS CORRECTION COEFFICIENT

! LOGICAL
! LDMAF      : MIXING LENGTH FOR MOIST AF
!-----------------------------------------------------------------------

! -   ARGUMENTS IMPLICITES.
!     ---------------------

! COMMON/YOMPHY0/
! COMMON/YOMCST /

!-----------------------------------------------------------------------

!     Externals.
!     ---------

!     Method.
!     --------
!     Bougeault and Lacarrere (1989), Deadorff (1974),
!     Cheng, Canuto and Howard (2002), ...
!     and the implementation note:
!     F. Vana: Turbulent mixing length closure for the ALARO physics.

!     Author.
!     -------
!        F. Vana (LACE/CHMI) 

!     Modifications.
!     --------------
!        Original version 20-Sep-2007
!        F. Vana 26-Mar-2008: consistency with eTKE + more options
!        F. Vana 21-Jan-2009: phased with CY35T1
!        F. Vana 17-Jun-2009: extended by X-term computation
!        F. Vana 26-Jun-2009: reworked EL2, EL3, EL4, EL5
!        F. Vana 20-Aug-2009: L-type mixing length as direct output
!                             plus special 1st timestep treatment
!        I. Bastak Duran 12-Oct-2009: updated modules for eTKE
!        F. Vana 03-Jun-2010: All computations rewritten to 'moist'
!                               Brunt-Vaissala freq.
!        2011-06: M. Jerczynski - some cleaning to meet norms
!        I. Bastak-Duran, F. Vana 27-Jan-2012: Various fixes
!        I. Bastak-Duran, R. Brozkova 19-Jan-2013: Multiple corrections.
!        I. Bastak-Duran, R. Brozkova 23-Oct-2014: BV computation correction. 
!        M. Hrastinski, J. Masek Oct-2020:
!          Moved moist gustiness correction from ACMRIP, fixed EL1 option:
!          corrected computation of BL89 integrals with added shear term,
!          geometric averaging of L_up and L_down, correct conversion between
!          L_TKE and lm, smooth transition to kappa.(z+z0) in surface layer,
!          no parcel crossing. EL2 option recoded as combination of EL1 in
!          unstable case and Grisogono and Belusic (2008) extension of
!          Deardorff (1980) mixing length in stable case. Removed
!          experimental options EL3-EL5.
! End modifications
!-------------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK

USE YOMPHY   , ONLY : TPHY
USE YOMCST   , ONLY : RG       ,RATM     ,RKAPPA  ,RETV
USE YOMPHY0  , ONLY : TPHY0
USE YOM_YGFL , ONLY : TYPE_GFLD

IMPLICIT NONE

TYPE(TPHY)        ,INTENT(IN)    :: YDPHY
TYPE(TPHY0)       ,INTENT(IN)    :: YDPHY0
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLTYPE
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHI(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHIF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQLI(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQICE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQRAIN(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSNOW(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTKE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PR(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLNPR(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMN2PP(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFMGST(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLSH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLCH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGZ0(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCDN(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBLH(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLMU(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLMT(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLMLF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLML(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLMLTILD(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRRCOR(KLON,0:KLEV)
LOGICAL           ,INTENT(IN)    :: LDMAF 

INTEGER(KIND=JPIM) :: JLEV, JLON, JJLEV

REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZINTE(KLON),ZLWORK(KLON)

REAL(KIND=JPRB) :: ZEPS,ZEPS1,ZEPSLML,ZX
REAL(KIND=JPRB) :: ZQV,ZQC,ZTESTM,ZTEST0,ZTEST1,ZKONV,ZRKONV
REAL(KIND=JPRB) :: ZDPHI,ZGLU,ZZ
REAL(KIND=JPRB) :: ZCD,ZKUROV,ZFP,ZLMGB08
REAL(KIND=JPRB) :: ZF0,ZF1,ZFMIN,ZINTE_NEG,ZPOTE,ZSIG,ZTKES,ZWEIGHT

REAL(KIND=JPRB) :: ZVPT(KLON,KLEV+1)
REAL(KIND=JPRB) :: ZLM(KLON,KLEV)
REAL(KIND=JPRB) :: ZLUP(KLON,KLEV)
REAL(KIND=JPRB) :: ZLDOWN(KLON,KLEV)
REAL(KIND=JPRB) :: ZMOIST(KLON,KLEV) 
REAL(KIND=JPRB) :: ZC3(KLON,0:KLEV) 
REAL(KIND=JPRB) :: ZLMU(KLON,0:KLEV) 
REAL(KIND=JPRB) :: ZDZZ (KLON,0:KLEV)
REAL(KIND=JPRB) :: ZDZZ2(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZH(KLON,KLEV)
REAL(KIND=JPRB) :: ZSHEAR2(KLON,0:KLEV) 
REAL(KIND=JPRB) :: ZU(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZCISA(KLON,0:KLEV),ZDPHIA(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZRTIA(KLON,0:KLEV)

!---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ACMIXELEN',0,ZHOOK_HANDLE)
ASSOCIATE(GCISMIN=>YDPHY0%GCISMIN, VKARMN=>YDPHY0%VKARMN, &
 & C3TKEFREE=>YDPHY0%C3TKEFREE, TKEMULT=>YDPHY0%TKEMULT, &
 & ALMAV=>YDPHY0%ALMAV, ALMAVE=>YDPHY0%ALMAVE, &
 & NUPTKE=>YDPHY0%NUPTKE, C_EPSILON=>YDPHY0%C_EPSILON, &
 & RRSCALE=>YDPHY0%RRSCALE,RRGAMMA=>YDPHY0%RRGAMMA,UTILGUST=>YDPHY0%UTILGUST, &
 & ETKE_C0SHEAR=>YDPHY0%ETKE_C0SHEAR, ETKE_R1SIM=>YDPHY0%ETKE_R1SIM, &
 & ETKE_R2SIM=>YDPHY0%ETKE_R2SIM, ETKE_GB08A=>YDPHY0%ETKE_GB08A, &
 & ETKE_GB08B=>YDPHY0%ETKE_GB08B, YTKE_NL=>YGFL%YTKE_NL, &
 & LCOEFK_ML=>YDPHY%LCOEFK_ML, LCOEFKTKE=>YDPHY%LCOEFKTKE, &
 & LCOEFK_PL=>YDPHY%LCOEFK_PL, LRRGUST=>YDPHY%LRRGUST)
!---------------------------------------------------------------------

! 0. Initial settings

!   security constants
ZEPS     = EPSILON(1.0_JPRB)*100.0_JPRB
ZEPS1    = SQRT(TINY(1.0_JPRB))*100.0_JPRB
ZFMIN    = 1.0E-10_JPRB
ZINTE_NEG=-1.0_JPRB
ZEPSLML  = 1.0E-08_JPRB

! 0.1 Setup of constants

! conversion factor L_TKE -> l_m
ZKONV=NUPTKE*NUPTKE*NUPTKE/C_EPSILON
! reverse value of ZKONV
ZRKONV=1._JPRB/ZKONV

ZGLU=RG*ALMAV

! 0.2 Setup of fields

DO JLEV=KTDIA,KLEV-1
  DO JLON=KIDIA,KFDIA
    ! calculation of full level spacings and half level heights
    ZDZZ(JLON,JLEV)=(PAPHIF(JLON,JLEV)-PAPHIF(JLON,JLEV+1))/RG
    ZH  (JLON,JLEV)=(PAPHI (JLON,JLEV)-PAPHI (JLON,KLEV  ))/RG
  ENDDO
ENDDO
ZDZZ(KIDIA:KFDIA,KLEV)=(PAPHIF(KIDIA:KFDIA,KLEV)-PAPHI(KIDIA:KFDIA,KLEV))/RG
ZH  (KIDIA:KFDIA,KLEV)=0.0_JPRB

! calculation of ZDZZ**2
ZDZZ2(KIDIA:KFDIA,KTDIA:KLEV)=ZDZZ(KIDIA:KFDIA,KTDIA:KLEV)**2

! calculation of surface wind speed (needed for moist gustiness correction)
DO JLON=KIDIA,KFDIA
  ZDPHI=PAPHIF(JLON,KLEV)-PAPHI(JLON,KLEV)
  ZDPHIA(JLON,KLEV)=ZDPHI
  ZCISA(JLON,KLEV)=PU(JLON,KLEV)**2+PV(JLON,KLEV)**2+&
   &(GCISMIN*ZDPHI/RG)**2
  ZU(JLON,KLEV)=SQRT(ZCISA(JLON,KLEV))
ENDDO

DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    ! calculation of moist contribution to theta
    ZQV=PQV(JLON,JLEV)
    ZQC=PQLI(JLON,JLEV)+PQICE(JLON,JLEV)+PQRAIN(JLON,JLEV)+PQSNOW(JLON,JLEV)
    ZMOIST(JLON,JLEV)=(1._JPRB+RETV*ZQV-ZQC)
  ENDDO
ENDDO

! Richardson number and Brunt Vaissala (moist) frequencey
DO JLEV=KTDIA,KLEV-1
  DO JLON=KIDIA,KFDIA
    ! shear computation
    ZZ=VKARMN*(0.5_JPRB*(PAPHIF(JLON,JLEV)+PAPHIF(JLON,JLEV+1))&
           & -PAPHI(JLON,KLEV)+PGZ0(JLON))
    ZDPHI=ZDZZ(JLON,JLEV)*RG
    ZDPHIA(JLON,JLEV)=PAPHIF(JLON,JLEV)-PAPHIF(JLON,JLEV+1)
    ZCISA(JLON,JLEV)=(PU(JLON,JLEV)-PU(JLON,JLEV+1))**2 &
     &  +(PV(JLON,JLEV)-PV(JLON,JLEV+1))**2 &
     &  +(GCISMIN*ZDPHIA(JLON,JLEV)/(RG*(1.0_JPRB+ZZ/ZGLU)))**2
    ! above surface wind speed
    ZU(JLON,JLEV)=SQRT(ZCISA(JLON,JLEV))
  ENDDO
ENDDO

! Computation of shear for generalized BL89 method (QR2017)
DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    ZSHEAR2(JLON,JLEV)=ZCISA(JLON,JLEV)/ZDZZ2(JLON,JLEV)
  ENDDO
ENDDO

IF (KMLTYPE == 1.OR.KMLTYPE == 2) THEN  ! options where BL89 is included

  ! virtual potential temperature at the surface
  DO JLON=KIDIA,KFDIA
    ZVPT(JLON,KLEV+1)=PTS(JLON)*ZMOIST(JLON,KLEV)* &
     & (RATM/PAPRS(JLON,KLEV))**RKAPPA
  ENDDO

  ! virtual potential temperature at full levels
  DO JLEV=KTDIA,KLEV
    DO JLON=KIDIA,KFDIA
      ZVPT(JLON,JLEV)=PT(JLON,JLEV)*ZMOIST(JLON,JLEV)* &
       & (RATM/PAPRSF(JLON,JLEV))**RKAPPA
    ENDDO
  ENDDO

ENDIF

! C3 computation
ZC3(KIDIA:KFDIA,KTDIA:KLEV-1)=PLMT(KIDIA:KFDIA,KTDIA:KLEV-1)/ &
 &                            PLMU(KIDIA:KFDIA,KTDIA:KLEV-1)
ZC3(KIDIA:KFDIA,KLEV)=1._JPRB

!---------------------------------------------------------------------
! 1. Mixing length computation

IF (LCOEFK_PL) THEN
  PLMLF(KIDIA:KFDIA,KTDIA:KLEV)=MAX(PLMLF(KIDIA:KFDIA,KTDIA:KLEV),ZEPS)
ENDIF

IF ((KSTEP == 0).AND.(YTKE_NL%NREQIN == 0)) THEN

  ! special treatment for the first timestep
  DO JLEV=KTDIA,KLEV-1
    DO JLON=KIDIA,KFDIA
      ZLMU(JLON,JLEV)=PLMU(JLON,JLEV)
    ENDDO
  ENDDO

  PLMLTILD(KIDIA:KFDIA,KTDIA:KLEV-1)=ZLMU(KIDIA:KFDIA,KTDIA:KLEV-1)*ZRKONV
  PLML(KIDIA:KFDIA,KTDIA:KLEV-1)=PLMLTILD(KIDIA:KFDIA,KTDIA:KLEV-1)

  IF (LCOEFK_PL) THEN
    DO JLEV=KTDIA+1,KLEV-1
      DO JLON=KIDIA,KFDIA
        PLMLF(JLON,JLEV)=(PALPH(JLON,JLEV)/PLNPR(JLON,JLEV)*PLMLTILD(JLON,JLEV-1) +&
     &          (1.0_JPRB-PALPH(JLON,JLEV)/PLNPR(JLON,JLEV))*PLMLTILD(JLON,JLEV))
      ENDDO
    ENDDO
  ENDIF

ELSE

  ! any other obvious timestep or when TKE is initialized

  !---------------------------------------------------------------------
  ! 1.1. computation of the BL89 mixing length

  IF (KMLTYPE == 1.OR.KMLTYPE == 2) THEN

    ! 1.1.0. main loop over model full levels

    DO JLEV=KTDIA,KLEV

      ! 1.1.1. mixing length for a downwards displacement

      ZINTE (KIDIA:KFDIA)=TKEMULT*PTKE(KIDIA:KFDIA,JLEV)
      ZLWORK(KIDIA:KFDIA)=0._JPRB
      ZTESTM             =1._JPRB
      DO JJLEV=JLEV,KLEV
        IF (ZTESTM <= 0._JPRB) EXIT
        ZTESTM=0._JPRB
        DO JLON=KIDIA,KFDIA
          ZTEST0=0.5_JPRB+SIGN(0.5_JPRB,ZINTE(JLON))
          ZTESTM=ZTESTM+ZTEST0
          ZF0   =(RG*(ZVPT(JLON,JLEV)/ZVPT(JLON,JJLEV  )-1.0_JPRB)+       &
           & ETKE_C0SHEAR*SQRT(PTKE(JLON,JJLEV)*ZSHEAR2(JLON,JJLEV)))* &
           & ZDZZ(JLON,JJLEV)
          ZF1   =(RG*(ZVPT(JLON,JLEV)/ZVPT(JLON,JJLEV+1)-1.0_JPRB)+       &
           & ETKE_C0SHEAR*SQRT(PTKE(JLON,JJLEV)*ZSHEAR2(JLON,JJLEV)))* &
           & ZDZZ(JLON,JJLEV)
          ZPOTE =0.5_JPRB*(ZF0+ZF1)
          ZTEST1=0.5_JPRB+SIGN(0.5_JPRB,MAX(2.0_JPRB*ZINTE(JLON),ZF0)*    &
           & (2.0_JPRB*ZINTE(JLON)-ZF0)-2.0_JPRB*ZINTE(JLON)*ZF1)
          ZX=(SQRT(ABS(ZF0*ZF0+2.0_JPRB*ZINTE(JLON)*(ZF1-ZF0)))-ZF0+      &
           & ZEPS1*ZINTE(JLON))/                                          &
           & SIGN(MAX(ABS(ZF1-ZF0),ZEPS1*(ABS(ZF0)+ZFMIN)),ZF1-ZF0)
          ZLWORK(JLON)=ZLWORK(JLON)+ZTEST0*(ZTEST1+(1.0_JPRB-ZTEST1)*ZX)* &
           & ZDZZ(JLON,JJLEV)
          ZINTE(JLON) =ZTEST0*ZTEST1*(ZINTE(JLON)-ZPOTE)+                 &
           & (1.0_JPRB-ZTEST0*ZTEST1)*ZINTE_NEG
        ENDDO
      ENDDO

      ! 1.1.2. intermediate storage (and security) of L_down

      ZLDOWN(KIDIA:KFDIA,JLEV)=ZLWORK(KIDIA:KFDIA)

      ! 1.1.3. mixing length for an upwards displacement
      !        (Note: The uppermost level computation is skipped,
      !        so results in L_up=0 there.)

      ZINTE (KIDIA:KFDIA)=TKEMULT*PTKE(KIDIA:KFDIA,JLEV)
      ZLWORK(KIDIA:KFDIA)=0._JPRB
      ZTESTM             =1._JPRB
      DO JJLEV=JLEV,KTDIA+1,-1
        IF (ZTESTM <= 0._JPRB) EXIT
        ZTESTM=0._JPRB
        DO JLON=KIDIA,KFDIA
          ZTEST0=0.5_JPRB+SIGN(0.5_JPRB,ZINTE(JLON))
          ZTESTM=ZTESTM+ZTEST0
          ZF0   =(RG*(1.0_JPRB-ZVPT(JLON,JLEV)/ZVPT(JLON,JJLEV  ))+       &
           & ETKE_C0SHEAR*SQRT(PTKE(JLON,JJLEV)*ZSHEAR2(JLON,JJLEV-1)))* &
           & ZDZZ(JLON,JJLEV-1)
          ZF1   =(RG*(1.0_JPRB-ZVPT(JLON,JLEV)/ZVPT(JLON,JJLEV-1))+       &
           & ETKE_C0SHEAR*SQRT(PTKE(JLON,JJLEV)*ZSHEAR2(JLON,JJLEV-1)))* &
           & ZDZZ(JLON,JJLEV-1)
          ZPOTE =0.5_JPRB*(ZF0+ZF1)
          ZTEST1=0.5_JPRB+SIGN(0.5_JPRB,MAX(2.0_JPRB*ZINTE(JLON),ZF0)*    &
           & (2.0_JPRB*ZINTE(JLON)-ZF0)-2.0_JPRB*ZINTE(JLON)*ZF1)
          ZX=(SQRT(ABS(ZF0*ZF0+2.0_JPRB*ZINTE(JLON)*(ZF1-ZF0)))-ZF0+      &
           & ZEPS1*ZINTE(JLON))/ &
           & SIGN(MAX(ABS(ZF1-ZF0),ZEPS1*(ABS(ZF0)+ZFMIN)),ZF1-ZF0)
          ZLWORK(JLON)=ZLWORK(JLON)+ZTEST0*(ZTEST1+(1.0_JPRB-ZTEST1)*ZX)* &
           & ZDZZ(JLON,JJLEV-1)
          ZINTE(JLON) =ZTEST0*ZTEST1*(ZINTE(JLON)-ZPOTE)+                 &
           & (1.0_JPRB-ZTEST0*ZTEST1)*ZINTE_NEG
        ENDDO
      ENDDO

      ! 1.1.4. intermediate storage of L_up

      ZLUP(KIDIA:KFDIA,JLEV)=ZLWORK(KIDIA:KFDIA)

    ! end of the loop on the vertical levels
    ENDDO

    ! 1.1.5. treat crossing parcels

    ! L_down
    DO JLEV=KTDIA+1,KLEV
      DO JLON=KIDIA,KFDIA
        ZLDOWN(JLON,JLEV)=MAX(ZLDOWN(JLON,JLEV),ZLDOWN(JLON,JLEV-1)+ &
         & (PAPHIF(JLON,JLEV)-PAPHIF(JLON,JLEV-1))/RG,ZEPS)
      ENDDO
    ENDDO

    ! L_up
    DO JLEV=KLEV-1,KTDIA,-1
      DO JLON=KIDIA,KFDIA
        ZLUP(JLON,JLEV)=MAX(ZLUP(JLON,JLEV),ZLUP(JLON,JLEV+1)+ &
         & (PAPHIF(JLON,JLEV+1)-PAPHIF(JLON,JLEV))/RG,ZEPS)
      ENDDO
    ENDDO

    ! 1.1.6. averaging of L_up and L_down, application of minimum value
    !        in the free atmosphere
    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA
        ZLM(JLON,JLEV)=MAX(SQRT(ZLUP(JLON,JLEV)*ZLDOWN(JLON,JLEV)),ALMAVE)
      ENDDO
    ENDDO
   
    ! 1.1.7 conversion of full level L-type mixing length to half level l-type
    DO JLEV=KTDIA,KLEV-1
      DO JLON=KIDIA,KFDIA
        ZLMU(JLON,JLEV)=0.5_JPRB*(ZLM(JLON,JLEV)+ZLM(JLON,JLEV+1))*ZKONV
      ENDDO
    ENDDO

  ENDIF  ! KMLTYPE == 1.OR.KMLTYPE == 2


  !---------------------------------------------------------------------
  ! 2. Combinations of the different types of mixing lengths into final one

  SELECT CASE (KMLTYPE)

    ! 2.0. unmodified Geleyn-Cedilnik (l_gc) solution
    CASE(0)
      ZLMU(KIDIA:KFDIA,KTDIA:KLEV-1)=PLMU(KIDIA:KFDIA,KTDIA:KLEV-1)

    ! 2.1. BL89 - nothing to do
    CASE(1)

    ! 2.2. BL89 in unstable case, GB08 in stable case
    CASE(2)
      DO JLEV=KTDIA,KLEV-1
        DO JLON=KIDIA,KFDIA

          ! stability flag:
          !   ZSIG=0 for PMN2PP <= 0
          !   ZSIG=1 for PMN2PP >  0
          ZSIG=0.5_JPRB-SIGN(0.5_JPRB,-PMN2PP(JLON,JLEV))

          ! half level TKE
          ZTKES=0.5_JPRB*(PTKE(JLON,JLEV)+PTKE(JLON,JLEV+1))

          ! ZLMGB89 is always computed, but used only for PMN2PP > 0
          ZLMGB08=SQRT(ZTKES*MIN( &
           & ETKE_GB08A*ETKE_GB08A/MAX(ZEPS1,ABS(PMN2PP(JLON,JLEV))), &
           & ETKE_GB08B*ETKE_GB08B/ZSHEAR2(JLON,JLEV)))

          ! minimum value and conversion from L-type mixing length to l-type
          ZLMGB08=MAX(ZLMGB08,ALMAVE)*ZKONV

          ! final combination of BL89 and GB08 (half level)
          ZLMU(JLON,JLEV)=(1._JPRB-ZSIG)*ZLMU(JLON,JLEV)+ZSIG*ZLMGB08
        ENDDO
      ENDDO

  ENDSELECT
 
  ! push TKE based mixing lengths to kappa.(z+z0) in the surface layer
  IF (KMLTYPE == 1.OR.KMLTYPE == 2) THEN
    DO JLEV=KTDIA,KLEV-1
      DO JLON=KIDIA,KFDIA
        ZWEIGHT=MAX(0._JPRB,MIN(1._JPRB, &
         & (ETKE_R2SIM-ZH(JLON,JLEV)/PBLH(JLON))/(ETKE_R2SIM-ETKE_R1SIM)))
        ZWEIGHT=ZWEIGHT*ZWEIGHT*(3._JPRB-2._JPRB*ZWEIGHT)
        ZLMU(JLON,JLEV)=ZLMU(JLON,JLEV)+ &
         & ZWEIGHT*(VKARMN*(ZH(JLON,JLEV)+PGZ0(JLEV)/RG)-ZLMU(JLON,JLEV))
      ENDDO
    ENDDO
  ENDIF

  PLMLTILD(KIDIA:KFDIA,KTDIA:KLEV-1)=ZLMU(KIDIA:KFDIA,KTDIA:KLEV-1)*ZRKONV
ENDIF 

IF (LCOEFK_PL.AND.(.NOT.LDMAF)) THEN
  IF (KSTEP == 0) THEN
    DO JLEV=KTDIA+1,KLEV
      DO JLON=KIDIA,KFDIA
        PLMLF(JLON,JLEV)=(PALPH(JLON,JLEV)/PLNPR(JLON,JLEV)*PLMLTILD(JLON,JLEV-1) +&
     &          (1.0_JPRB-PALPH(JLON,JLEV)/PLNPR(JLON,JLEV))*PLMLTILD(JLON,JLEV))
      ENDDO
    ENDDO
  ELSE
    ! conversion of prognostic L on full level to half level l_m
    DO JLEV=KTDIA,KLEV-1
      DO JLON=KIDIA,KFDIA
        PLML(JLON,JLEV)=(PLMLF(JLON,JLEV)+PLMLF(JLON,JLEV+1))*0.5_JPRB
        ZLMU(JLON,JLEV)=PLML(JLON,JLEV)*ZKONV
      ENDDO
    ENDDO
  ENDIF
ELSE
  PLML(KIDIA:KFDIA,KTDIA:KLEV-1)=PLMLTILD(KIDIA:KFDIA,KTDIA:KLEV-1)
ENDIF


!---------------------------------------------------------------------
! 3. Computation and application of moist gustiness correction 

IF (LRRGUST) THEN

  ! computation part

  ! surface
  DO JLON=KIDIA,KFDIA
    ZCD=PFMGST(JLON,KLEV)*PCDN(JLON)
    ZFP=MAX(0.0_JPRB,PFPLSH(JLON,KLEV)+PFPLCH(JLON,KLEV))
    PRRCOR(JLON,KLEV)=SQRT(SQRT(1.0_JPRB+((((ZFP/(ZFP+RRSCALE))**RRGAMMA)* &
     & UTILGUST)**2)/(ZCD*ZU(JLON,KLEV)**2)))
  ENDDO

  ! remaining levels
  DO JLEV=KTDIA,KLEV-1
    DO JLON=KIDIA,KFDIA
      ZDPHI=PAPHIF(JLON,JLEV)-PAPHIF(JLON,JLEV+1)
      ZRTIA(JLON,JLEV)=2.0_JPRB/(PR(JLON,JLEV)*PT(JLON,JLEV)&
      &+PR(JLON,JLEV+1)*PT(JLON,JLEV+1))
      ZKUROV=MAX(PFMGST(JLON,JLEV)*ZU(JLON,JLEV)* (PLMU(JLON,JLEV)*RG)**2 &
       & *PAPRS(JLON,JLEV)*ZRTIA(JLON,JLEV)/ZDPHI**2,ZEPS)
      ZFP=MAX(0.0_JPRB,PFPLSH(JLON,JLEV)+PFPLCH(JLON,JLEV))
      PRRCOR(JLON,JLEV)=SQRT(SQRT(1.0_JPRB+((((ZFP/(ZFP+RRSCALE))**RRGAMMA) &
       & *UTILGUST)**2)*PAPRS(JLON,JLEV)*ZRTIA(JLON,JLEV)/(ZU(JLON,JLEV)*ZKUROV)))
    ENDDO
  ENDDO

  ! application part
  DO JLEV=KTDIA,KLEV-1
    DO JLON=KIDIA,KFDIA
      ZLMU(JLON,JLEV)=ZLMU(JLON,JLEV)*PRRCOR(JLON,JLEV)
      PLML(JLON,JLEV)=PLML(JLON,JLEV)*PRRCOR(JLON,JLEV)
    ENDDO
  ENDDO

ENDIF

!---------------------------------------------------------------------
! 4. Treatment for upper and lower bounds and final export
!    to output arrays

! lower bound (note, that TOMS solver requires full level value at the bb)
ZLMU(KIDIA:KFDIA,KLEV)= PGZ0(KIDIA:KFDIA)*VKARMN/RG
PLML(KIDIA:KFDIA,KLEV)=ZRKONV*ZLMU(KIDIA:KFDIA,KLEV)
PLMLTILD(KIDIA:KFDIA,KLEV)=PLML(KIDIA:KFDIA,KLEV)

IF (KSTEP == 0 .AND. LCOEFK_PL) THEN
  PLMLF(KIDIA:KFDIA,KLEV)=(PALPH(KIDIA:KFDIA,KLEV)/PLNPR(KIDIA:KFDIA,KLEV)*&
   & PLML(KIDIA:KFDIA,KLEV-1) +&
   & (1.0_JPRB-PALPH(KIDIA:KFDIA,KLEV)/PLNPR(KIDIA:KFDIA,KLEV))*&
   & PLML(KIDIA:KFDIA,KLEV))
  PLMLF(KIDIA:KFDIA,KTDIA)=PLMLF(KIDIA:KFDIA,KTDIA+1)
  PLMLF(KIDIA:KFDIA,KTDIA:KLEV)=MAX(ZEPSLML,PLMLF(KIDIA:KFDIA,KTDIA:KLEV))
ENDIF

! uppermost level is probably not needed but anyway...
ZLMU(KIDIA:KFDIA,KTDIA-1)=ZLMU(KIDIA:KFDIA,KTDIA)
! is needed for 3D turbulence
PLML(KIDIA:KFDIA,KTDIA-1)=PLML(KIDIA:KFDIA,KTDIA)

PLMLTILD(KIDIA:KFDIA,KTDIA-1)=PLMLTILD(KIDIA:KFDIA,KTDIA)

! Finaly overwrite the original mixing length by the new one.
PLMU(KIDIA:KFDIA,KTDIA-1:KLEV)=ZLMU(KIDIA:KFDIA,KTDIA-1:KLEV)

! mixing length for heat
IF (LCOEFKTKE) THEN
  ! no vertical profile of C3
  PLMT(KIDIA:KFDIA,KTDIA:KLEV)=C3TKEFREE*PLMU(KIDIA:KFDIA,KTDIA:KLEV)
ELSE
  ! vertical profile of C3
  PLMT(KIDIA:KFDIA,KTDIA:KLEV)=ZC3(KIDIA:KFDIA,KTDIA:KLEV)*&
  & PLMU(KIDIA:KFDIA,KTDIA:KLEV)
ENDIF

! ----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACMIXELEN',1,ZHOOK_HANDLE)

END SUBROUTINE ACMIXELEN
