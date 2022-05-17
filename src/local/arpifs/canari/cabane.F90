SUBROUTINE CABANE(YDGEOMETRY)
!****---------------------------------------------------------------------------
!****   CABANE : INITIALISATION DES PARAMETRES DU MODELE STATISTIQUE
!****   ------
!****---------------------------------------------------------------------------
!  BUT : DEUXIEME ETAPE POUR L'INITIALISATION DES PARAMETRES DU MODELE
!  ---   STATISTIQUE, MISE A JOUR DE TOUS LES COEFFICIENTS.

!  ROUTINE APPELANTE :    *CALIFE*
!  -----------------
!  SOUS-PROGRAMMES APPELES :   GPPRE - CABINE -
!  -----------------------
!  ARGUMENTS D'ENTREE :   common  QAREF
!  ------------------
!  ARGUMENTS DE SORTIE :   common  QACOSS
!  -------------------

!  MODIFICATIONS :
!  -------------
!   M.Hamrud      01-Oct-2003 CY28 Cleaning
!   R. El Khatib 01-04-2011 Bugfix for B-level distribution
!   G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!   K. Yessad (Nov 2012): call GPHPRE.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!-------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
USE YOMANCS  , ONLY : RADIANS  ,RDEGREES
USE YOMCT0   , ONLY : LELAM
USE YOMMP0   , ONLY : MY_REGION_EW
USE QALORI   , ONLY : QCAGUE   ,RCALPH
USE QAPASS   , ONLY : JPNVS    ,JPBLS
USE QAREF    , ONLY : REF_STAT ,REF_A_SST,REF_A_T2 ,REF_A_H2 ,REF_A_SN ,&
 & REF_A_VOR1,REF_A_DIV1,REF_NU_BL,REF_KP_BL,&
 & REF_S_SN ,REF_S_T2 ,REF_S_H2 ,REF_S_V1 ,&
 & REF_MU1  ,REF_MU2  ,REF_MU3  ,REF_NU1  ,&
 & REF_NU2  ,REF_NU3  ,REF_PHUD ,REF_PHHU ,&
 & REF_PVH  ,REF_COEFN,REF_COEFT,REF_COEFS,REF_AP_SN  
USE QACOSS   , ONLY : XPPRS    ,XPLRS    ,XPOHZ    ,&
 & XPOHUD   ,XPOSST   ,XPOHHU   ,XRMU     ,XRNU     ,&
 & XPOVZ    ,XPOVH    ,XPOHT2   ,XPOHH2   ,&
 & XPOVORBL ,XPODIVBL ,XNUBL    ,XKPBL    ,&
 & XPOHSN   ,XPOHPSN  

!-------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN)   :: YDGEOMETRY
INTEGER(KIND=JPIM) , ALLOCATABLE :: ISTAGP(:)
REAL(KIND=JPRB) :: ZPHILO(JPNVS,-90:90), ZHULOC(JPNVS,-90:90)
REAL(KIND=JPRB) :: ZPRESH(0:YDGEOMETRY%YRDIMV%NFLEVG), ZPRESF(YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IENDLAT, ILAT, IND, IROF, ISTLAT,&
 & JB, JC, JGL, JL, JLIMNS, JLIMTR, JLON  

REAL(KIND=JPRB) :: ZCOEFN, ZCOEFS, ZCOEFT, ZFPHHU, ZFPHUD,&
 & ZFPHZ, ZFPVH, ZFPVZ, ZGM, ZXL, ZCOEF  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------------

#include "cabine.intfb.h"
#include "gphpre.intfb.h"

!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CABANE',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,    &
& YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDVAB=>YDGEOMETRY%YRVAB)
ASSOCIATE(NDGENL=>YDDIM%NDGENL,   NFLEVG=>YDDIMV%NFLEVG,   NSTAGP=>YDGEM%NSTAGP,   NONL=>YDMP%NONL, NPTRFLOFF=>YDMP%NPTRFLOFF&
& )

!**----------------------------------------------------------------------------
!**  - 1 - Initialisation de coefficients.
!**        ------------------------------

ZFPHZ  = 1.0_JPRB
ZFPHUD = REF_PHUD
ZFPHHU = REF_PHHU

ZFPVZ = 1.0_JPRB
ZFPVH = REF_PVH

ZCOEFN = REF_COEFN
ZCOEFT = REF_COEFT
ZCOEFS = REF_COEFS

!**----------------------------------------------------------------------------
!**  - 2 - Initialisation des parametres du modele statistique.
!**        ---------------------------------------------------

DO JB = 1 , JPNVS
  XPPRS(JB) = REF_STAT(JPNVS+1-JB,1)
  XPLRS(JB) = LOG(REF_STAT(JPNVS+1-JB,1))
ENDDO
DO JB = 1 , JPNVS
  DO JL = -JPBLS , JPBLS+1
    XPOHZ(JL,JB)  = ZFPHZ  * REF_STAT(JPNVS+1-JB,7)
    XPOHUD(JL,JB) = ZFPHUD * REF_STAT(JPNVS+1-JB,7)
    XPOHHU(JL,JB) = ZFPHHU * REF_STAT(JPNVS+1-JB,7)
    ZXL  = REAL(JL,JPRB) - 0.5_JPRB
! coefficient d'ageostrophisme par latitude
    XRMU(JL,JB) = (REF_MU1-MAX(0.0_JPRB,REF_MU2*(XPPRS(JB)-REF_MU3)&
     & / (XPPRS(JPNVS)-REF_MU3)))&
     & * SIN(RADIANS*ZXL*10._JPRB)  
! rapport E.T. vent divergent / E.T. vent rotationnel par latitude
    XRNU(JL,JB) = (REF_NU1+MAX(0.0_JPRB,REF_NU2*(XPPRS(JB)-REF_NU3)&
     & / (XPPRS(JPNVS)-REF_NU3)))&
     & * MAX(COS(RADIANS*ZXL*10._JPRB),SQRT(3._JPRB)*0.5_JPRB)  
    XPOVZ(JL,JB) = ZFPVZ * REF_STAT(JPNVS+1-JB,6)
    XPOVH(JL,JB) = ZFPVH * REF_STAT(JPNVS+1-JB,6)
  ENDDO
ENDDO

IF (LELAM) THEN
!pour neutraliser la multiplication en fin de programme et dans CABINE
  ZGM   = 3.5_JPRB
  ZCOEF = EXP(RCALPH*(ZGM-1.0_JPRB/ZGM))
ELSE
  ZCOEF = 1.0_JPRB
ENDIF

XPOSST   = REF_A_SST * ZCOEF
XPOHT2   = REF_A_T2  * ZCOEF
XPOHH2   = REF_A_H2  * ZCOEF
XPOHSN   = REF_A_SN  * ZCOEF
XPOHPSN  = REF_AP_SN
XPOVORBL = REF_A_VOR1* ZCOEF
XPODIVBL = REF_A_DIV1* ZCOEF
XNUBL    = REF_NU_BL
XKPBL    = REF_KP_BL

REF_S_SN = REF_S_SN * ZCOEF**2
REF_S_T2 = REF_S_T2 * ZCOEF**2
REF_S_H2 = REF_S_H2 * ZCOEF**2
REF_S_V1 = REF_S_V1 * ZCOEF**2

JLIMNS = 25
JLIMTR = 15

DO JC = 1 , JPNVS
  DO JL = -90 , -JLIMNS
    ZPHILO(JC,JL) = REF_STAT(JPNVS+1-JC,2) * ZCOEFS
    ZHULOC(JC,JL) = REF_STAT(JPNVS+1-JC,5)
  ENDDO
  DO JL = -(JLIMNS-1) , -(JLIMTR+1)
    ZPHILO(JC,JL) = REF_STAT(JPNVS+1-JC,2) * ( ZCOEFS +&
     & (ZCOEFT-ZCOEFS) * REAL(JL+JLIMNS,JPRB) / REAL(JLIMNS-JLIMTR,JPRB) )  
    ZHULOC(JC,JL) = REF_STAT(JPNVS+1-JC,5)
  ENDDO
  DO JL = -JLIMTR , JLIMTR
    ZPHILO(JC,JL) = REF_STAT(JPNVS+1-JC,2) * ZCOEFT
    ZHULOC(JC,JL) = REF_STAT(JPNVS+1-JC,5)
  ENDDO
  DO JL = JLIMTR+1 , JLIMNS-1
    ZPHILO(JC,JL) = REF_STAT(JPNVS+1-JC,2) * ( ZCOEFT +&
     & (ZCOEFN-ZCOEFT) * REAL(JL-JLIMTR,JPRB) / REAL(JLIMNS-JLIMTR,JPRB) )  
    ZHULOC(JC,JL) = REF_STAT(JPNVS+1-JC,5)
  ENDDO
  DO JL = JLIMNS , 90
    ZPHILO(JC,JL) = REF_STAT(JPNVS+1-JC,2) * ZCOEFN
    ZHULOC(JC,JL) = REF_STAT(JPNVS+1-JC,5)
  ENDDO
ENDDO

!**----------------------------------------------------------------------------
!**  - 3 - Interpolations verticales.
!**        -------------------------

!  ILAT est l'indice de la latitude la plus proche dans ZPHILO et ZHULOC.
!  Calcul par GPPRE des niveaux pressions de la verticale de reference
!  ou calculer les ecart-types, et interpolations par CABINE.

! NB : cela fait que les data sont en niveaux sigma variable avec YRVAB%VAH et YRVAB%VBH !
ZPRESH(NFLEVG)=101325._JPRB
CALL GPHPRE(1,NFLEVG,1,1,YDVAB,ZPRESH,PRESF=ZPRESF)

ISTLAT=1
IENDLAT=NDGENL

IF (LELAM) THEN
  ALLOCATE(ISTAGP(ISTLAT:IENDLAT))
  IROF = 1
  DO JGL=ISTLAT,IENDLAT
    ISTAGP(JGL)=IROF
    IROF=IROF+NONL(NPTRFLOFF+JGL,MY_REGION_EW)
  ENDDO
ENDIF
DO JGL = ISTLAT , IENDLAT
  IROF=NSTAGP(JGL)-1
  IF (LELAM) IROF=ISTAGP(JGL)-1
  DO JLON = 1,NONL(NPTRFLOFF+JGL,MY_REGION_EW)
    IND = IROF+JLON
    ZGM = YDGSGEOM_NB%GM(IND)
    IF (LELAM) ZGM = 3.5_JPRB
    QCAGUE(2,IND) = EXP( RCALPH * (1.0_JPRB/ZGM - ZGM) )
    QCAGUE(1,IND) = QCAGUE(2,IND)**2
    ILAT = NINT(RDEGREES*YDGSGEOM_NB%GELAT(IND))
    CALL CABINE(YDGEOMETRY%YRDIMV,ZPHILO(1,ILAT),ZHULOC(1,ILAT),ZPRESF,IND)
  ENDDO
ENDDO
IF (ALLOCATED(ISTAGP)) DEALLOCATE (ISTAGP)

!**----------------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CABANE',1,ZHOOK_HANDLE)
END SUBROUTINE CABANE
