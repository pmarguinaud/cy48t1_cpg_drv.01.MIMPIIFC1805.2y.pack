SUBROUTINE CACLSST(YDGEOMETRY,YDSURF,YDPHY1)
!****--------------------------------------------------------------------------
!****   CACLSST : CONSTITUTION D'UN CHAMP DE RAPPEL CLIMATOLOGIQUE 
!****   -------    POUR LA TEMPERATURE DE SURFACE DE LA MER
!****
!****   Auteur :       F. TAILLEFER  07/99
!****--------------------------------------------------------------------------
!  BUT : LECTURE DANS LE FICHIER "ICMSHxxxxSST" D'UN CHAMP REFERENCE DE SST
!  ---   (L'ANALYSE NESDIS, RECUPEREE DANS LA BDAP ET INTERPOLEE SUR LA GRILLE
!         DU MODELE PAR UNE CONFIGURATION 931, EXECUTEE AU PREALABLE).
!        MISE A JOUR AVEC NOTRE CLIM SUR LA BANQUISE.
!        LE CHAMP AINSI ELABORE SERVIRA DE RAPPEL AU CHAMP DE SST ANALYSE ICI
!        DANS CANARI.
!        CE CHAMP EST STOCKE DANS LE BUFFER DES CHAMPS POINT DE GRILLE DU
!        MODELE.
!        SI LE SWITCH QUI COMMANDE L'ELABORATION DE CE CHAMP (NSSTLIS) N'EST
!        PAS ACTIVE, ON PREND LA CLIMATOLOGIE (CONF. 923) PARTOUT.

!  ROUTINE APPELANTE :    *CANARI*
!  -----------------

!  SOUS-PROGRAMMES APPELES :  routines d'utilisation des fichiers LFI
!  -----------------------    INCGPF DISGRID
!                             FAITOU FADIES FACILE FAIRME
!                             CHIEN CCHIEN UPDCAL

!  ARGUMENT D'ENTREE :
!  -----------------

!  ARGUMENT DE SORTIE :
!  ------------------

!  MODIFICATIONS :
!  -------------
!    M.Hamrud      01-Oct-2003 CY28 Cleaning
!    M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!    M.Hamrud      01-Jul-2006 Revised surface fields
!    K. Yessad (Jan 2010): externalisation of group EGGX in XRD/IFSAUX
!    R. El Khatib : 23-Apr-2010 use disgrid_mod instead of disgrid
!    T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!    K. Yessad (July 2014): Move some variables.
!      ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE YOMVERT  , ONLY : VP00
USE YOMLUN   , ONLY : NULOUT   ,NULCL1   ,NULCL2
USE YOMMP0   , ONLY : MYPROC   ,NPROC
USE YOMCT0   , ONLY : CNMEXP   ,NQUAD    ,LELAM
USE YOMPHY1  , ONLY : TPHY1
USE YOMRIP0  , ONLY : NINDAT
USE YOMTAG   , ONLY : MTAGGSUM
USE QACTEX   , ONLY : NSSTLIS  ,RCLISST  ,NSEAICE  ,LECSST
USE QACLIM   , ONLY : LCLIM
USE DISGRID_MOD, ONLY : DISGRID_SEND, DISGRID_RECV
USE IOGRIDA_MOD, ONLY : RUNDEF
USE MPL_MODULE

!      ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY),INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)   ,INTENT(INOUT) :: YDSURF
TYPE(TPHY1)   ,INTENT(INOUT) :: YDPHY1
REAL(KIND=JPRB) :: ZICE(YDGEOMETRY%YRGEM%NGPTOT), ZICE2(YDGEOMETRY%YRGEM%NGPTOT), ZSST(YDGEOMETRY%YRGEM%NGPTOT), ZCLIM(YDGEOMETRY%YRGEM%NGPTOT),&
 & ZGEO(YDGEOMETRY%YRGEM%NGPTOT),&
 & ZFIELD(YDGEOMETRY%YRGEM%NGPTOT), ZHIS(YDGEOMETRY%YRGEM%NGPTOTG), ZGAUSS(YDGEOMETRY%YRDIM%NDGSAG:YDGEOMETRY%YRDIM%NDGENG),&
 & ZDATO(3),&
 & ZMASK(YDGEOMETRY%YRGEM%NGPTOT)  

INTEGER(KIND=JPIM) :: IDATEF(11), ILMO(12)
CHARACTER (LEN = 22) ::  CLNOMF
LOGICAL :: LLEX

INTEGER(KIND=JPIM) :: IEND, IHIS, INDEX, IA2, IM2, IJ2, ILDAT, ITAG,&
 & IPROLEC, IREP, IREP2, IST, JL, JSTGLO, IADAT, IFDAT, IA1, IM1, IJ1, IINF, IBLK,&
 & INECSST

REAL(KIND=JPRB) :: ZSL1, ZSL2, ZSL3, ZLICE, ZCAL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZFMIN, ZFMAX, ZFMEAN
LOGICAL         :: LLUNDEF,LDGELATO
REAL(KIND=JPRB) :: ZLUNDEF

!      ------------------------------------------------------------------

#include "fcttim.func.h"
#include "chien.h"
#include "cchien.intfb.h"
#include "updcal.intfb.h"

!      ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CACLSST',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDVAB=>YDGEOMETRY%YRVAB&
& )
ASSOCIATE(TMERGL=>YDPHY1%TMERGL,   NDGENG=>YDDIM%NDGENG, NDGLG=>YDDIM%NDGLG, NDGSAG=>YDDIM%NDGSAG,                 &
& NDLON=>YDDIM%NDLON, NPROMA=>YDDIM%NPROMA, NSMAX=>YDDIM%NSMAX,   NFLEVG=>YDDIMV%NFLEVG,   NGPTOT=>YDGEM%NGPTOT,   &
& NHTYP=>YDGEM%NHTYP,   NLOENG=>YDGEM%NLOENG, NMENG=>YDGEM%NMENG, NSTTYP=>YDGEM%NSTTYP,   RLOCEN=>YDGEM%RLOCEN,    &
& RMUCEN=>YDGEM%RMUCEN, RSTRET=>YDGEM%RSTRET,   SD_VX=>YDSURF%SD_VX, SP_CI=>YDSURF%SP_CI, YSD_VX=>YDSURF%YSD_VX,   &
& YSP_CI=>YDSURF%YSP_CI)
!      ------------------------------------------------------------------

!**----------------------------------------------------------------------------
!**   1.    Recuperation du champ climatologique de SST.
!**         -------------------------------------------

WRITE(NULOUT,'(/,''entree subroutine CACLSST '')')

IF (RCLISST == 0.0_JPRB ) THEN
  WRITE(NULOUT,*) '   on ne veut pas faire de rappel climatologique '
  WRITE(NULOUT,*) '     --> inutile d''elaborer le champ ! '
  WRITE(NULOUT,*) '  '
  IF (LHOOK) CALL DR_HOOK('CACLSST',1,ZHOOK_HANDLE)
  RETURN
ENDIF

IF (.NOT.LCLIM) THEN
  WRITE(NULOUT,*) '   la lecture des fichiers climatologiques a ete impossible'
  WRITE(NULOUT,*) '     --> aucun rappel ne sera applique !'
  WRITE(NULOUT,*) '  '
  IF (LHOOK) CALL DR_HOOK('CACLSST',1,ZHOOK_HANDLE)
  RETURN
ENDIF

IST=1
DO JSTGLO = 1,NGPTOT,NPROMA
  IEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
  IBLK=(JSTGLO-1)/NPROMA+1
  DO JL = IST,IEND
    ZCLIM(JSTGLO+JL-1)=SD_VX(JL,YSD_VX%YTSC%MP,IBLK)
    ZGEO(JSTGLO+JL-1) =SD_VX(JL,YSD_VX%YORO%MP,IBLK)
  ENDDO
ENDDO

IPROLEC=1

!**----------------------------------------------------------------------------
!**   2.    Traitement de la partie banquise.
!**         --------------------------------

WRITE(NULOUT,*) ' Traitement de la banquise '

IF (NSEAICE == 0) THEN
  WRITE(NULOUT,*) '   on ne veut pas lire les champs de glace de mer '
  ZICE (:)=-999.9_JPRB
  ZICE2(:)=-999.9_JPRB
ELSE

  !**   2.2 CONCENTRATION GLACE DE MER
  !ouverture fichier de concentration de glace de mer et controle validite
  CLNOMF='ICMSH'//CNMEXP(1:4)//'SEAICE_CONC'
  INQUIRE(FILE=CLNOMF,EXIST=LLEX)
  IF (.NOT.LLEX) THEN
    WRITE(NULOUT,*)'   fichier de concentration de glace de mer absent'
    WRITE(NULOUT,*)'   pas d''actualisation de la concentration de banquise '
    NSEAICE=0
    ZICE (:)=-999.9_JPRB
  ELSE
    INDEX=10
    IF (MYPROC == IPROLEC) THEN
      WRITE(NULOUT,*)'   controle fichier fourni '
      CALL FAITOU(IREP,NULCL1,.TRUE.,CLNOMF,'UNKNOWN',.TRUE.,.FALSE.,0,&
       & INDEX,IHIS,'CADRE SEAICE')  
      IINF=-1
      IF (.NOT.LELAM) THEN
        CALL CHIEN('CADRE SEAICE    ',NSTTYP,RMUCEN,RLOCEN,RSTRET,NSMAX,NDGLG,&
         & NDLON,NLOENG,NMENG,NHTYP,NFLEVG,VP00,YDVAB%VALH,YDVAB%VBH,NQUAD,&
         & IINF,NDGSAG,NDGENG,1.0E-7_JPRB,LLEX,NULOUT)  
      ELSE
        CALL CCHIEN(YDGEOMETRY,'CADRE SEAICE',NULCL1,IINF)
        LLEX=.FALSE.
      ENDIF
      CALL FADIES(IREP,NULCL1,IDATEF)
      IADAT=IDATEF(1)*10000+IDATEF(2)*100+IDATEF(3)
      CALL FALAIS(IREP,NULCL1,'DATA.DATE',ZDATO,3)
      IFDAT=NINT(ZDATO(1))*10000+NINT(ZDATO(2))*100+NINT(ZDATO(3))
      IA1=NCCAA(NINDAT)
      IM1=NMM(NINDAT)
      IJ1=NDD(NINDAT)
      CALL UPDCAL(IJ1,IM1,IA1,-NSEAICE,IJ2,IM2,IA2,ILMO,99)
      ILDAT=IA2*10000+IM2*100+IJ2
      WRITE(NULOUT,'(''    date ANALYSE : '',I8,/,''    date fichier glace : '',&
       & I8,/,'' date donnees glace : '',I8)') NINDAT,IADAT,IFDAT  
      IF ( IADAT < ILDAT .OR. IADAT > NINDAT ) THEN
        WRITE(NULOUT,'(''  retard autorise : '',I2,'' jours donc '',&
         & ''MAUVAISE DATE fichier'',/,''    --> banquise non '',&
         & ''actualisee !'')') NSEAICE  
        NSEAICE=0
        OPEN(UNIT=63,FORM='FORMATTED')
        WRITE(63,'('' PB DATE FICHIER GLACE '')')
        CLOSE(63)
      ELSE
        WRITE(NULOUT,'(''  retard autorise : '',I2,'' jours donc OK '')') NSEAICE
      ENDIF
      !**   2.3  lecture du champ de concentration de glace de mer
      IF (NSEAICE > 0) THEN
        CALL FACILE(IREP,NULCL1,'SURF',0,'SEAICE.CONC',ZHIS,.FALSE.)
        CALL DISGRID_SEND(YDGEOMETRY,1,ZHIS,1,ZICE)
      ENDIF
      CALL FAIRME(IREP,NULCL1,'UNKNOWN')
    ENDIF
    IF (NPROC > 1) THEN
      ITAG=MTAGGSUM+147
      CALL MPL_BROADCAST(NSEAICE,KROOT=IPROLEC,KTAG=ITAG,CDSTRING='CACLSST:')
      IF (MYPROC /= IPROLEC .AND. NSEAICE > 0) CALL DISGRID_RECV(YDGEOMETRY,IPROLEC,1,ZICE,1)
    ENDIF
  ENDIF

  !**   2.4 HAUTEUR GLACE DE MER
  !ouverture fichier de hauteur de glace de mer
  CLNOMF='ICMSH'//CNMEXP(1:4)//'SEAICE_THICK'
  INQUIRE(FILE=CLNOMF,EXIST=LLEX)
  IF (.NOT.LLEX) THEN
    WRITE(NULOUT,*)'   fichier de hauteur de glace de mer absent'
    ZICE2(:)=-999.9_JPRB
  ELSE
    INDEX=10
    IF (MYPROC == IPROLEC) THEN
      WRITE(NULOUT,*)'   controle fichier fourni '
      CALL FAITOU(IREP2,NULCL1,.TRUE.,CLNOMF,'UNKNOWN',.TRUE.,.FALSE.,0,&
       & INDEX,IHIS,'CADRE SEAICE')  
      !**   2.5  lecture du champ de hauteur de glace de mer
      CALL FACILE(IREP2,NULCL1,'SURF',0,'SEAICE.THICK',ZHIS,.FALSE.)
      CALL DISGRID_SEND(YDGEOMETRY,1,ZHIS,1,ZICE2)
      CALL FAIRME(IREP2,NULCL1,'UNKNOWN')
    ENDIF
    IF (NPROC > 1) THEN
      ITAG=MTAGGSUM+147
      CALL MPL_BROADCAST(NSEAICE,KROOT=IPROLEC,KTAG=ITAG,CDSTRING='CACLSST:')
      IF (MYPROC /= IPROLEC .AND. NSEAICE > 0) CALL DISGRID_RECV(YDGEOMETRY,IPROLEC,1,ZICE2,1)
    ENDIF
  ENDIF

ENDIF      ! fin NSEAICE > 0

!**----------------------------------------------------------------------------
!**   3.    Traitement de la partie SST NESDIS.
!**         ----------------------------------

WRITE(NULOUT,*) ' Traitement SST '
ZMASK(:)=0.0_JPRB

!**   3.1  on ne veut pas utiliser de SST NESDIS

IF (NSSTLIS == 0) THEN
  WRITE(NULOUT,*) '   on ne veut pas utiliser de produit de SST  '
ELSE

!**   3.2  ouverture fichier SST et controles validite

  CLNOMF='ICMSH'//CNMEXP(1:4)//'SST'
  INQUIRE(FILE=CLNOMF,EXIST=LLEX)
  IF (.NOT.LLEX) THEN
    WRITE(NULOUT,*)'   fichier du produit de SST absent'
    WRITE(NULOUT,*)'     --> pas de rappel vers cette SST !'
    NSSTLIS=0
  ELSE
    INDEX=10
    IF (MYPROC == IPROLEC) THEN
      WRITE(NULOUT,*)'   controle fichier fourni '
      CALL FAITOU(IREP,NULCL2,.TRUE.,CLNOMF,'UNKNOWN',.TRUE.,.FALSE.,0,&
       & INDEX,IHIS,'CADRE SST')  
      IINF=-1
      IF (.NOT.LELAM) THEN
        CALL CHIEN('CADRE SST       ',NSTTYP,RMUCEN,RLOCEN,RSTRET,NSMAX,NDGLG,NDLON,&
         & NLOENG,NMENG,NHTYP,NFLEVG,VP00,YDVAB%VALH,YDVAB%VBH,NQUAD,IINF,&
         & NDGSAG,NDGENG,1.0E-7_JPRB,LLEX,NULOUT)  
      ELSE
        CALL CCHIEN(YDGEOMETRY,'CADRE SST',NULCL2,IINF)
        LLEX=.FALSE.
      ENDIF
      CALL FADIES(IREP,NULCL2,IDATEF)
      IFDAT=IDATEF(1)*10000+IDATEF(2)*100+IDATEF(3)
      IA1=NCCAA(NINDAT)
      IM1=NMM(NINDAT)
      IJ1=NDD(NINDAT)
      CALL UPDCAL(IJ1,IM1,IA1,-NSSTLIS,IJ2,IM2,IA2,ILMO,99)
      ILDAT=IA2*10000+IM2*100+IJ2
      WRITE(NULOUT,'(''    date ANALYSE : '',I8,''    date SST NESDIS : '',&
       & I8)') NINDAT,IFDAT  
      IF ( IFDAT < ILDAT .OR. IFDAT > NINDAT ) THEN
        WRITE(NULOUT,'(''  retard autorise : '',I2,'' jours donc '',&
         & ''MAUVAISE DATE fichier'',/,''    --> application '',&
         & ''rappel climatologique classique !'')') NSSTLIS  
        NSSTLIS=0
        OPEN(UNIT=64,FORM='FORMATTED')
        WRITE(64,'('' PB DATE FICHIER SST NESDIS '')')
        CLOSE(64)
      ELSE
        WRITE(NULOUT,'(''  retard autorise : '',I2,'' jours donc OK '')') NSSTLIS
      ENDIF

!**   3.3  lecture champ SST de rappel

      IF (NSSTLIS > 0) THEN
        LLUNDEF = .TRUE.
        ZLUNDEF = RUNDEF
        CALL FACILO (IREP,NULCL2,'SURF',0,'SST.CLIM.',ZHIS,.FALSE.,&
                     & LLUNDEF,ZLUNDEF)
        ZSST(1:NGPTOT)=ZHIS(1:NGPTOT)
        CALL DISGRID_SEND(YDGEOMETRY,1,ZHIS,1,ZSST)
      ENDIF
      CALL FAIRME(IREP,NULCL2,'UNKNOWN')
    ENDIF

    IF (NPROC > 1) THEN
      ITAG=MTAGGSUM+148
      CALL MPL_BROADCAST(NSSTLIS,KROOT=IPROLEC,KTAG=ITAG,CDSTRING='CACLSST:')
      IF (MYPROC /= IPROLEC .AND. NSSTLIS > 0) CALL DISGRID_RECV(YDGEOMETRY,IPROLEC,1,ZSST,1)
    ENDIF

  ENDIF             ! fin fichier present
ENDIF               ! fin NSSTLIS /= 0

!**----------------------------------------------------------------------------
!**   4.    Elaboration du champ de rappel definitif.
!**         ----------------------------------------

ZSL1=271.4_JPRB         ! limite validite SST NESDIS interpolee  (~ -1.75 degC)
ZCAL=0.05_JPRB          ! enveloppe estimee erreurs interpolations
ZSL2=TMERGL-ZCAL        ! transformation mer en banquise  (~ -1.97 degC)
ZSL3=TMERGL+ZCAL        ! transformation mer en eau       (~ -1.87 degC)
ZLICE=60._JPRB          ! seuil pourcentage occurence glace --> banquise

! on ne fait pas de distinction terre/mer car les valeurs sur terre ne
! sont pas importantes (elles ne seront jamais utilisees !)

!**   4.1  mise au niveau modele de la SST NESDIS (conservation nature)

!first check if gelato sea ice model is activated or not
LDGELATO = .FALSE.
CALL TEST_GELATO(LDGELATO)

IF (NSSTLIS > 0) THEN
  DO JL = 1,NGPTOT
    IF (ZSST(JL) > ZSL1) ZMASK(JL)=1.0_JPRB
    ZHIS(JL)=MAX(ZSL3,ZSST(JL))
  ENDDO
ENDIF

!**   4.2  cas sans banquise actualisee

IF (NSEAICE == 0) THEN
  ZFIELD(:)=ZCLIM(:)
  DO JL = 1,NGPTOT
    IF (ZCLIM(JL) > TMERGL .AND. ZMASK(JL) > 0.5_JPRB) ZFIELD(JL)=ZHIS(JL)
  ENDDO
ELSEIF (LDGELATO) THEN
  DO JL = 1,NGPTOT
    ZFIELD(JL)=ZHIS(JL)
  ENDDO
ELSE

!**   4.3  cas avec actualisation de la banquise

  DO JL = 1,NGPTOT
    IF (ZICE(JL) < ZLICE) THEN
      IF (ZMASK(JL) > 0.5_JPRB) THEN
        ZFIELD(JL)=ZHIS(JL)
      ELSE
        ZFIELD(JL)=MAX(ZCLIM(JL),ZSL3)
      ENDIF
    ELSE
      ZFIELD(JL)=MIN(ZCLIM(JL),ZSL2)
    ENDIF
  ENDDO

ENDIF

!**----------------------------------------------------------------------------
!**   5.    Read ECMWF SST and SIC and replace climatological SST field with 
!**         ECMWF SST over open sea, NOT over lakes and sea ice
!**         -------------------------------------------

IF (LECSST) THEN
  INECSST = 1
ELSE
  INECSST = 0
ENDIF

IF (INECSST == 1) THEN

  WRITE(NULOUT,*)' Try to read ECMWF SST and SIC'
  CLNOMF='ICMSHANALESST'
  WRITE(NULOUT,*)'      CLNOMF:',CLNOMF
  IPROLEC=1
  INQUIRE(FILE=CLNOMF,EXIST=LLEX)
  IF (.NOT.LLEX) THEN
    WRITE(NULOUT,*)' No ECMWF SST and SIC'
    LECSST = .FALSE.
  ELSE
    IHIS = 0
    IF (MYPROC == IPROLEC) THEN
      CALL FAITOU(IREP,NULCL1,.TRUE.,CLNOMF,'UNKNOWN',.TRUE.,&
       & .FALSE.,1,1,IHIS,'CADRE SST')  
      CALL FADIES(IREP,NULCL1,IDATEF)
      IFDAT=IDATEF(1)*10000+IDATEF(2)*100+IDATEF(3)
      IA1=NCCAA(NINDAT)
      IM1=NMM(NINDAT)
      IJ1=NDD(NINDAT)
!**    CALL UPDCAL(IA1,IM1,IJ1,-NSSTLIS,IJ2,IM2,IA2,ILMO,99)
!**    ILDAT=IA2*10000+IM2*100+IJ2
      WRITE(NULOUT,'(''    date ANALYSE : '',I8,''    date SST ECMWF : '',&
       & I8)') NINDAT,IFDAT  
!**    IF ( IFDAT < ILDAT .OR. IFDAT > NINDAT ) THEN
!**        WRITE(NULOUT,'(''  retard autorise : '',I2,'' jours donc '', &
!**         & ''MAUVAISE DATE fichier'',/,''    --> application '', &
!**         & ''rappel climatologique classique !'')') NSSTLIS  
!**        NSSTLIS=0
!**      ELSE
!**        WRITE(NULOUT,'(''  retard autorise : '',I2,'' jours donc OK '')') NSSTLIS
!**      ENDIF

!**   Read ECMWF SST and SIC
      WRITE(NULOUT,*) 'Read SEA.TEMPERA'
      CALL FACILE(IREP,NULCL1,'SURF',0,'SEA.TEMPERA',ZHIS,.FALSE.)
      CALL DISGRID_SEND(YDGEOMETRY,1,ZHIS,1,ZSST)

      ZFMIN = MINVAL(ZSST)
      ZFMAX = MAXVAL(ZSST)
      ZFMEAN = SUM(ZSST)/FLOAT(NGPTOT)
      WRITE(NULOUT,*)'  ZSST - MIN  =', ZFMIN 
      WRITE(NULOUT,*)'  ZSST - MAX  =', ZFMAX 
      WRITE(NULOUT,*)'  ZSST - MEAN =', ZFMEAN 

      WRITE(NULOUT,*) 'Read SEA.ICECONC'
      CALL FACILE(IREP,NULCL1,'SURF',0,'SEA.ICECONC',ZHIS,.FALSE.)
      CALL DISGRID_SEND(YDGEOMETRY,1,ZHIS,2,ZICE)
 
      ZFMIN = MINVAL(ZICE)
      ZFMAX = MAXVAL(ZICE)
      ZFMEAN = SUM(ZICE)/FLOAT(NGPTOT)
      WRITE(NULOUT,*)'  ZICE - MIN  =', ZFMIN 
      WRITE(NULOUT,*)'  ZICE - MAX  =', ZFMAX 
      WRITE(NULOUT,*)'  ZICE - MEAN =', ZFMEAN 
      CALL FAIRME(IREP,NULCL1,'UNKNOWN')
    ENDIF

    IF (NPROC > 1) THEN
      ITAG=MTAGGSUM+146
      CALL MPL_BROADCAST(INECSST,KROOT=IPROLEC,KTAG=ITAG,CDSTRING='CACLSST:')
      IF (MYPROC /= IPROLEC .AND. INECSST>0) THEN
        CALL DISGRID_RECV(YDGEOMETRY,IPROLEC,1,ZSST,1)
        CALL DISGRID_RECV(YDGEOMETRY,IPROLEC,1,ZICE,2)
      ENDIF
    ENDIF

  ENDIF

!** Replace climatological SST with ECMWF SST,
!** where this field NOT is missing (= -9999.0)
!** 
  IF (LECSST) THEN     
    DO JL = 1,NGPTOT
      IF (ZSST(JL) > 0.0_JPRB) ZFIELD(JL)=ZSST(JL)
    ENDDO
  ENDIF

ENDIF

!**----------------------------------------------------------------------------
!**   6.    Ecriture du champ dans le buffer.
!**         --------------------------------

! position du champ de rappel SST dans le buffer
DO JSTGLO = 1,NGPTOT,NPROMA
  IEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
  IBLK=(JSTGLO-1)/NPROMA+1
  DO JL = IST,IEND
    SP_CI(JL,YSP_CI%YCI(8)%MP0,IBLK)=ZFIELD(JSTGLO+JL-1)    !SST
    SP_CI(JL,YSP_CI%YCI(14)%MP0,IBLK)=  ZICE  (JSTGLO+JL-1) !SIC 
    SP_CI(JL,YSP_CI%YCI(15)%MP0,IBLK)=  ZICE2 (JSTGLO+JL-1) !SIT
  ENDDO
ENDDO

WRITE(NULOUT,'('' Ecriture champ de rappel SST terminee'',/)')

!**----------------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CACLSST',1,ZHOOK_HANDLE)
END SUBROUTINE CACLSST

