SUBROUTINE PPOBSAC(YDCST, YDGP5,KDCOUNT,CDPAR,PINTH,KOBSTYPE,KCOUNT,&
 & PALTI,KDPP,PXPP,LDDIRCLSMOD)
!****-------------------------------------------------------------------
!****  PPOBSAC : CANARI SURFACE POST-PROCESSING
!****            FIND OUT AN EQUIVALENT MODEL PARAMETER FOR A GIVEN OBS
!****  called from SURFACEO
!****-------------------------------------------------------------------

!     WHERE YDGP5   : MODEL VARIABLE AT OBS LOCATION
!           KDCOUNT : SECOND DIMENSION (Data)
!           CDPAR   : PARAMETER FOR POST-PROCESSING
!           PINTH   : INTERPOLATION HEIGHT
!           KOBSTYPE: OBS. TYPE
!           KCOUNT  : COUNTER OF DATA PER OBS
!           PALTI   : OBS. OROGRAPHY
!           KDPP    : Number of fields, dimension of PXPP.
!           PXPP    : POST-PROCESSED VALUES
!

!**   EXTERNALS : ACHMT  ACSOLW
!     ---------

!     MODIFIED
!     --------
!     R. El Khatib 01-08-07 : Pruning options
!     J.M. Piriou / J.F Geleyn  01-08-23 : New definition and use of the thermal mixing length.
!     M. Bellus  2003-04-24 : Change of arguments in ACHMT calling
!                             + removing LSOLV and some cleaning
!     M. Tudor   2003-10-21 : Change args in ACHMT, reintroduce LSOLV
!     M.Hamrud      01-Jul-2006 Revised surface fields
!     A. Alias   2007-10-11 : Change args in ACHMT, PRTI is introduced
!     C. Payan   2008-11    : Neutral Wind
!     K. Yessad (Jul 2009)  : remove CDLOCK + some cleanings
!     J.Hague    21-Mar-11  : YDGOM Derived Type added
!     C. Soci    27-02-2012 : LCORGRADT option for T2m
!     A. Geer    26-Mar-12  : GOM modernisation
!     A. Geer    27-Jul-12  : Remove LSM rejection hack that hijacked VF_EMISF
!     T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     A. Geer     12-Aug-2015 Pre-OOPS: - profile variables intent(inout) - Meteo-France to fix.
!                                       - Remove direct acces to GOM arrays
!     T. Montmerle  Aug 2018: Better use of objects in arguments
!****-------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMANCS  , ONLY : RMDI
USE YOMCST   , ONLY : TCST
USE YOMSTA   , ONLY : RDTDZ1
USE YOMOBS   , ONLY : LSCATT_NEUTRAL ,LCORGRADT
USE YOMPHY1  , ONLY : YRPHY1
USE YOMPHY2  , ONLY : YRPHY2
USE YOMPHY   , ONLY : YRPHY
USE YOMCOCTP , ONLY : NSYNOP, NSCATT
USE YOMDB
USE GOM_PLUS , ONLY : TYPE_GOM_PLUS, IH

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(TYPE_GOM_PLUS), INTENT(IN)  :: YDGP5
INTEGER(KIND=JPIM),INTENT(IN)    :: KDCOUNT
CHARACTER(LEN=*)  ,INTENT(IN)    :: CDPAR 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PINTH(YDGP5%NDLEN) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KOBSTYPE 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCOUNT(YDGP5%NDLEN) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALTI(YDGP5%NDLEN) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDPP
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXPP(YDGP5%NDLEN,KDCOUNT,KDPP) 
LOGICAL           ,INTENT(IN)    :: LDDIRCLSMOD
!-----------------------------------------------------------------------
REAL(KIND=JPRB) :: ZDPHIT(YDGP5%NDLEN),                ZDPHIV(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZARG(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZWFC(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZGZ0R(YDGP5%NDLEN),                 ZNEIJG(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZNEIJV(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZNBVNO(YDGP5%NDLEN,0:YDGP5%NFLEVG)       
REAL(KIND=JPRB) :: ZMRIPP(YDGP5%NDLEN,0:YDGP5%NFLEVG)
REAL(KIND=JPRB) :: ZFPLSH(YDGP5%NDLEN,0:YDGP5%NFLEVG),       ZFPLCH(YDGP5%NDLEN,0:YDGP5%NFLEVG)
REAL(KIND=JPRB) :: ZCD(YDGP5%NDLEN),                   ZCDN(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZCDROV(YDGP5%NDLEN),                ZCH(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZCHROV(YDGP5%NDLEN),                ZCPS(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZDQSTS(YDGP5%NDLEN),                ZGWDCS(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZGZ0(YDGP5%NDLEN),                  ZHQ(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZGZ0H(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZHU(YDGP5%NDLEN),                   ZNEIJ(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZQCLS(YDGP5%NDLEN),                 ZQS(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZQSATS(YDGP5%NDLEN),                ZRHCLS(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZRS(YDGP5%NDLEN), ZRTI(YDGP5%NDLEN),      ZSTAB(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZTCLS(YDGP5%NDLEN),                 ZUCLS(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZVCLS(YDGP5%NDLEN),                 ZVEG(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZNUCLS(YDGP5%NDLEN),                ZNVCLS(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZXDROV(YDGP5%NDLEN),                ZXHROV(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZURAF(YDGP5%NDLEN),                 ZVRAF(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZPCLS(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZLSM(YDGP5%NDLEN), ZSNS(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZGZ0F(YDGP5%NDLEN), ZGZ0HF(YDGP5%NDLEN)
REAL(KIND=JPRB) :: ZRHMAX,ZRHMIN,ZINCR(YDGP5%NDLEN),ZZARG,ZEPS
REAL(KIND=JPRB) :: ZSPEED5
INTEGER(KIND=JPIM) :: JROF,IMXCNT,JC,ILEN
LOGICAL :: LLDCLS,LLDHMT,LLZ0H
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "abor1.intfb.h"
#include "achmt.intfb.h"
#include "actkehmt.intfb.h"
#include "chksurf.intfb.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PPOBSAC',0,ZHOOK_HANDLE)
ASSOCIATE(EWFC=>YRPHY1%EWFC, GWFC=>YRPHY1%GWFC, &
 & HVCLS=>YRPHY2%HVCLS, HTCLS=>YRPHY2%HTCLS, &
 & LFGEL=>YRPHY%LFGEL, LSOLV=>YRPHY%LSOLV, LCOEFKSURF=>YRPHY%LCOEFKSURF)
!-----------------------------------------------------------------------

!**---------------------------------------------------------------------
!**  1.  PREPARATION DU CALCUL DES PARAMETRES AUX HAUTEURS METEOS.
!**      --------------------------------------------------------
ILEN=YDGP5%NDLEN
IMXCNT = MAXVAL(KCOUNT(1:ILEN))
IF(IMXCNT == 0 .AND. LHOOK) CALL DR_HOOK('PPOBSAC',1,ZHOOK_HANDLE)
IF(IMXCNT == 0) RETURN
ZNEIJV(:)=0.0_JPRB
ZNEIJG(:)=0.0_JPRB
ZGZ0R(:)=0.0_JPRB
ZFPLSH(:,:)=0.0_JPRB
ZFPLCH(:,:)=0.0_JPRB
ZLSM(:)=0.0_JPRB
ZSNS(:)=0.0_JPRB
ZGZ0F(:)=0.0_JPRB
ZGZ0HF(:)=0.0_JPRB

IF (CDPAR == 'U10' .OR. CDPAR == 'FFS' .OR. CDPAR == 'T2' .OR. CDPAR == 'H2' .OR. CDPAR == 'TS') THEN
  ! THESE ARE ACCEPTABLE PARAMETERS
ELSE
  WRITE(NULOUT,'('' *** REQUESTED PARAMETER='',A)') CDPAR
  CALL ABOR1('ERROR PPOBSAC: REQUESTED PARAMETER')
ENDIF

!**---------------------------------------------------------------------
!**  1.  RECUPERATION DES CHAMPS COUCHE LIMITE (SI SURFEX)
!**      -------------------------------------------------

IF (LDDIRCLSMOD) THEN
  ZTCLS(1:ILEN) = YDGP5%TCLS(1:ILEN,IH)
  ZRHCLS(1:ILEN)= YDGP5%HUCLS(1:ILEN,IH)
  ZUCLS(1:ILEN) = YDGP5%UCLS(1:ILEN,IH)
  ZVCLS(1:ILEN) = YDGP5%VCLS(1:ILEN,IH)
  ZNUCLS(1:ILEN)= YDGP5%NUCLS(1:ILEN,IH)
  ZNVCLS(1:ILEN)= YDGP5%NVCLS(1:ILEN,IH)
ELSE

!**---------------------------------------------------------------------
!**  2.  CALCUL DES PARAMETRES AUX HAUTEURS METEOS.
!**      -----------------------------------------

  ZDPHIT(1:ILEN) = YDCST%RG*HTCLS
  ZDPHIV(1:ILEN) = YDCST%RG*HVCLS

  IF(KOBSTYPE == NSYNOP)THEN
    DO JROF = 1 , ILEN
      ZDPHIV(JROF)=PINTH(JROF)
    ENDDO
  ENDIF 

  ZEPS=1.E-1_JPRB

  LLDCLS=.FALSE.
  LLDHMT=.TRUE.

! Calculs preliminaires pour le schema sol-vegetation (WFC,HV)
  IF ( LSOLV ) THEN
    ZARG=YDGP5%ARG(:,IH)
    DO JROF = 1 , ILEN
      LLZ0H = YDGP5%Z0H(JROF,IH) /= YDGP5%MISSING_VALUE
      IF ( YDGP5%LS(JROF,IH) <= 0.5_JPRB ) THEN
        ZWFC(JROF) = 1.0_JPRB
      ELSE
        ZZARG = MAX(ZEPS,ZARG(JROF))
        ZWFC(JROF) = GWFC*(ZZARG**EWFC)
      ENDIF
    ENDDO
  ENDIF

!**---------------------------------------------------------------------
!**  - 2 - Appel pour le calcul des parametres aux hauteurs meteo.
!**        ------------------------------------------------------

  IF (LCOEFKSURF) THEN
    CALL ACTKEHMT(1,ILEN,YDGP5%NDLEN,YDGP5%NFLEVG,LLZ0H,&
    & YDGP5%GEOPH(:,:,IH),YDGP5%GEOPF(:,:,IH),YDGP5%PRESH(:,:,IH),YDGP5%PRESF(:,:,IH),&
    & YDGP5%CP(:,:,IH),YDGP5%QF(:,:,IH),YDGP5%RF(:,:,IH),YDGP5%TF(:,:,IH),&
    & YDGP5%UF(:,:,IH),YDGP5%VF(:,:,IH),ZFPLSH,ZFPLCH,ZDPHIT,ZDPHIV,YDGP5%Z0(:,IH),&
    & YDGP5%Z0H(:,IH),ZGZ0R,YDGP5%HV(:,IH),YDGP5%LS(:,IH),ZNEIJG,ZNEIJV,YDGP5%SN(:,IH),&
    & YDGP5%TS(:,IH),YDGP5%VEG(:,IH),ZWFC,YDGP5%WS(:,IH),YDGP5%WSI(:,IH),&
    & LLDCLS,LLDHMT,&
    & ZNBVNO,ZMRIPP,ZCD,ZCDN,ZCDROV,ZCH,ZCHROV,ZCPS,ZDQSTS,ZGWDCS,&
    & ZGZ0,ZGZ0H,ZHQ,ZHU,ZNEIJ,ZQCLS,ZQS,ZQSATS,ZRHCLS,ZRS,ZRTI,&
    & ZSTAB,ZTCLS,ZUCLS,ZVCLS,ZNUCLS,ZNVCLS,ZPCLS,ZVEG,ZXDROV,ZXHROV,ZURAF,ZVRAF)
  ELSE
    CALL CHKSURF(1,ILEN,YDGP5,ZLSM,ZSNS,ZGZ0F,ZGZ0HF)
    CALL ACHMT(YDCST,1,ILEN,YDGP5%NDLEN,YDGP5%NFLEVG,LLZ0H,&
    & YDGP5%GEOPH(:,:,IH),YDGP5%GEOPF(:,:,IH),YDGP5%PRESH(:,:,IH),YDGP5%PRESF(:,:,IH),&
    & YDGP5%CP(:,:,IH),YDGP5%QF(:,:,IH),YDGP5%RF(:,:,IH),YDGP5%TF(:,:,IH),&
    & YDGP5%UF(:,:,IH),YDGP5%VF(:,:,IH),ZFPLSH,ZFPLCH,ZDPHIT,ZDPHIV,ZGZ0F,&
    & ZGZ0HF,ZGZ0R,YDGP5%HV(:,IH),ZLSM,ZNEIJG,ZNEIJV,ZSNS,YDGP5%TS(:,IH),&
    & YDGP5%VEG(:,IH),ZWFC,YDGP5%WS(:,IH),YDGP5%WSI(:,IH),&
    & LLDCLS,LLDHMT,&
    & ZNBVNO,ZMRIPP,ZCD,ZCDN,ZCDROV,ZCH,ZCHROV,ZCPS,ZDQSTS,ZGWDCS,&
    & ZGZ0,ZGZ0H,ZHQ,ZHU,ZNEIJ,ZQCLS,ZQS,ZQSATS,ZRHCLS,ZRS,ZRTI,&
    & ZSTAB,ZTCLS,ZUCLS,ZVCLS,ZNUCLS,ZNVCLS,ZPCLS,ZVEG,ZXDROV,ZXHROV,ZURAF,ZVRAF)
  ENDIF

ENDIF

!**---------------------------------------------------------------------
!**  - 3 - Stockage/initialisation des champs diagnostics
!**        necessaires a l'analyse.
!**        ----------------------------------------------

ZRHMAX=1.0_JPRB
ZRHMIN=0.0_JPRB

IF (CDPAR == 'H2') THEN
  PXPP(1:ILEN,1,1) = MAX(ZRHMIN,MIN(ZRHMAX,ZRHCLS(1:ILEN)))
  DO JC = 2,IMXCNT
    WHERE(KCOUNT(1:ILEN) >= JC) PXPP(1:ILEN,JC,1)=PXPP(1:ILEN,1,1)
  ENDDO
ELSEIF (CDPAR == 'U10' .OR. CDPAR == 'FFS')     THEN
  IF(KOBSTYPE == NSCATT .AND. LSCATT_NEUTRAL) THEN
    PXPP(1:ILEN,1,1) = ZNUCLS(1:ILEN)
    PXPP(1:ILEN,1,2) = ZNVCLS(1:ILEN)
  ELSE
    PXPP(1:ILEN,1,1) = ZUCLS(1:ILEN)
    PXPP(1:ILEN,1,2) = ZVCLS(1:ILEN)
  ENDIF
  DO JC = 2,IMXCNT
    WHERE(KCOUNT(1:ILEN) >= JC)
      PXPP(1:ILEN,JC,1)=PXPP(1:ILEN,1,1)
      PXPP(1:ILEN,JC,2)=PXPP(1:ILEN,1,2)
    ENDWHERE
  ENDDO
  IF (CDPAR == 'FFS') THEN
    DO JC=1,IMXCNT
      DO JROF=1,ILEN
        IF (JC > KCOUNT(JROF)) EXIT
        ZSPEED5 = SQRT(PXPP(JROF,JC,1)**2+PXPP(JROF,JC,2)**2)
        PXPP(JROF,JC,1) = MAX(ZSPEED5,1.E-10_JPRB)
        PXPP(JROF,JC,2) = RMDI
      ENDDO
    ENDDO
  ENDIF
ELSEIF (CDPAR == 'T2') THEN
  IF (LCORGRADT) THEN
    ZINCR(1:ILEN) = RDTDZ1*(YDGP5%OROG(1:ILEN,IH)-PALTI(1:ILEN))/YDCST%RG
    WHERE(ABS(PALTI(1:ILEN)) > 90000._JPRB) ZINCR(1:ILEN)=0.0_JPRB
    PXPP(1:ILEN,1,1) = ZTCLS(1:ILEN)-ZINCR(1:ILEN)
  ELSE
    PXPP(1:ILEN,1,1) = ZTCLS(1:ILEN)
  ENDIF
  DO JC = 2,IMXCNT
    WHERE(KCOUNT(1:ILEN) >= JC) PXPP(1:ILEN,JC,1)=PXPP(1:ILEN,1,1)
  ENDDO
ELSEIF (CDPAR == 'TS') THEN
!  cas des SST : l'altitude de l'obs est consideree a 0 metre.
  ZINCR(1:ILEN)=RDTDZ1*YDGP5%OROG(1:ILEN,IH)/YDCST%RG
  PXPP(1:ILEN,1,1) = YDGP5%SST(1:ILEN,IH)-ZINCR(1:ILEN)
  DO JC = 2,IMXCNT
    WHERE(KCOUNT(1:ILEN) >= JC) PXPP(1:ILEN,JC,1)=PXPP(1:ILEN,1,1)
  ENDDO
ELSE
  WRITE(NULOUT,'('' *** ERROR IN PPOBSAC'')')
  WRITE(NULOUT,'('' *** REQUESTED PARAMETER='',A)') CDPAR
ENDIF

!**---------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('PPOBSAC',1,ZHOOK_HANDLE)

END SUBROUTINE PPOBSAC
