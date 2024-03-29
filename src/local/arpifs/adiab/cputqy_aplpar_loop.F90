SUBROUTINE CPUTQY_APLPAR_LOOP(YDDYN, YDDIMV,YDGMV, &
 & YGFL,YDPTRSLB1,YDPHY,KPROMA, KST, KPROF, KFLEV, PDT, &
 !  Variables 2D Input
 & PTENDGFL,&
 !  Variables 2D In/Out
 & PB1, PGMVT1, PGFLT1)

!     -----------------------------------------------------------------
!**** *CPUTQY_APLPAR_LOOP*  EVOLUTION DE T, VENT , Q , QI , QL, QR, QS
!               ET  TKE.
!     -----------------------------------------------------------------

!     ARGUMENTS D'ENTREE
!     ------------------
!       KPROMA : DIMENSION HORIZONTALE
!       KSTART : DEBUT DE LA BOUCLE HORIZONTALE
!       KPROF : BORNE HORIZONTALE DES CALCULS
!       KFLEV : DIMENSION ET BORNE VERTICALE
!       PDT : PAS DE TEMPS EFFECTIF

!       PTENDH (KPROMA,KFLEV) : TENDANCE DE L'ENTHALPIE
!       PTENDQ ( IDEM ) : TENDANCE DE L'EAU VAPEUR
!       PTENDQI ( IDEM ) : TENDANCE DE L'EAU GLACE
!       PTENDQL ( IDEM ) : TENDANCE DE L'EAU LIQUIDE
!       PTENDQR ( IDEM ) : TENDANCE DE L'EAU PRECIPITANTE LIQUIDE
!       PTENDQS ( IDEM ) : TENDANCE DE L'EAU PRECIPITANTE SOLIDE
!       PTENDU, PTENDV (KPROMA,KFLEV) : TENDANCE DU VENT
!       PTENDU_ZDEC, PTENDV_ZDEC (KPROMA,KFLEV) : TENDANCE DU VENT POUR CALCULER ZDEC
!       PTENDTKE : TENDANCE DE L'ENERGIE KINETIQUE TURB.
!       PTENDT   : TENDANCE DE LA TEMPERATURE SI ON N'UTILISE PAS
!                  LES FLUXES, MAIS LES TENDANCES (HIRLAM OPTION)
!       PCP9 ( IDEM ) : CHALEUR SPECIFIQUE A T=9
!       PDELP9 ( IDEM ) : EPAISSEUR DE LA COUCHE A T=9
!       PTT9, PUT9, PVT9 (IDEM) : ETAT THERMODYNAMIQUE A T=9

!     ARGUMENTS IMPLICITES
!     --------------------
!       CONSTANTS UNIVERSELLES  = COMMON /YOMCST/:  RG,RCPD,RCPV,RCS,RCW

!     SORTIES
!     -------
!       PTT1, PQT1 (KPROMA,KFLEV): VARIABLES A FAIRE EVOLUER
!       PQIT1,PQLT1 (IDEM) : EAU SOL. ET LIQ. A FAIRE EVOLUER
!       PQST1,PQRT1 (IDEM) : EAU SOL. ET LIQ. PRECIPITANTE A FAIRE EVOLUER
!       PUT1, PVT1 (IDEM)    : VENT A FAIRE EVOLUER
!       PFDIS(KPROMA,0:KFLEV): FLUX ENTHALPIE DU A LA DISSIPATION EC
!       PTKET1               : ENERGIE KINETIQUE TURB. A FAIRE EVOLUER

! --- AUTHOR.
!     -------
!      E.Bazile .

! --- MODIFICATIONS.
!     --------------
!      07/10/04 Modification par F. Bouyssel (Lopez Ql/Qi/Qr prog. scheme)
!      25/01/06 Modification par F. Bouyssel (PTENDQS, PQST1)
!      01/06/06 Modification par B. Sass (PTENDTKE,PTKET1,PTENDT) 
!      30/10/06 Modification par F. Bouyssel (Bug correction found by M.Bellus)
!      02/10/06 Modification par M. Bellus (PTENDPTKE,LPTKE,L3MT,LSTRAPRO);
!               + corrected bug in ZTDCP comp. (missing snow, rain parts)
!      01/02/07 Modification par M.Janousek (ALARO pTKE is not *Dt, adding
!               pTKE tendency to T and enthalpy diss.flux)
!      26/04/10 Modification par Y.Bouteloup (Only one call to cputqy in mf_phys)
!      2013-11, D. Degrauwe: Generality in terms of GFL; added tendency 
!                            of D (NH).
!-----------------------------------------------------------------------

USE YOMDIMV  , ONLY : TDIMV
USE YOMGMV   , ONLY : TGMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCT0   , ONLY : LSLAG, LTWOTL, LNHDYN
USE YOMCST   , ONLY : RG, RCPD
USE YOMPHY   , ONLY : TPHY
USE YOM_YGFL , ONLY : TYPE_GFLD
USE PTRSLB1  , ONLY : TPTRSLB1
USE YOMDYN   , ONLY : TDYN
 
!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TGMV)        ,INTENT(IN)    :: YDGMV
TYPE(TPHY)        ,INTENT(IN)    :: YDPHY
TYPE(TPTRSLB1)    ,INTENT(IN)    :: YDPTRSLB1
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDGFL(KPROMA,KFLEV,YGFL%NUMFLDS)

REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1(KPROMA,YDPTRSLB1%NFLDSLB1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1(KPROMA,KFLEV,YDGMV%YT1%NDIM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLT1(KPROMA,KFLEV,YGFL%NDIM1)

#include "cp_ptrslb1.intfb.h"

!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JLEV, JROF, JPGFL, JGFL
INTEGER(KIND=JPIM) :: IMP1EXPLICIT (14)
INTEGER(KIND=JPIM) :: ISLB1U9  ,ISLB1V9  ,ISLB1T9  ,ISLB1GFL9, ISLB1VD9
INTEGER(KIND=JPIM) :: IPGFL(YGFL%NUMFLDS)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPUTQY_APLPAR_LOOP',0,ZHOOK_HANDLE)

ASSOCIATE(NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, NFLSA=>YDDIMV%NFLSA, NFLSUL=>YDDIMV%NFLSUL)

! * Calculation of IPGFL, since the old pointers
!   MSLB1[X]9 (=MSLB1GFL9+IP[X]) do not exist any longer in PTRSLB1.

! usefull pointer for new version of cputqy

DO JGFL=1,NUMFLDS
  IF ((YCOMP(JGFL)%MP1 > 0) .AND. (YCOMP(JGFL)%MP_SL1 > 0)) THEN
     IPGFL(YCOMP(JGFL)%MP1) = (YCOMP(JGFL)%MP_SL1-1)*(KFLEV+2*NFLSUL)
  ENDIF   
ENDDO  

IF (LSLAG) CALL CP_PTRSLB1(YDDYN,YDPTRSLB1,ISLB1U9,ISLB1V9,ISLB1T9,ISLB1VD9,ISLB1GFL9)

IMP1EXPLICIT = [YGFL%YL%MP1, YGFL%YI%MP1, YGFL%YS%MP1, YGFL%YR%MP1, YGFL%YG%MP1, &
& YGFL%YTKE%MP1, YGFL%YEFB1%MP1, YGFL%YEFB2%MP1, YGFL%YEFB3%MP1, YGFL%YLCONV%MP1, &
& YGFL%YICONV%MP1, YGFL%YRCONV%MP1, YGFL%YSCONV%MP1, YGFL%YQ%MP1]

DO JGFL=1,NUMFLDS
  IF (ALL (IMP1EXPLICIT /= YCOMP(JGFL)%MP1)) THEN
    IF (YCOMP(JGFL)%LT1) THEN
      IF (LSLAG .AND. YCOMP(JGFL)%LADV) THEN
        JPGFL=IPGFL(YCOMP(JGFL)%MP1)
        DO JLEV=1,KFLEV
          DO JROF=KST,KPROF
            PB1(JROF,ISLB1GFL9+JPGFL+JLEV-NFLSA)=&
             & PB1(JROF,ISLB1GFL9+JPGFL+JLEV-NFLSA)&
             & + PDT*PTENDGFL(JROF,JLEV,YCOMP(JGFL)%MP1)
          ENDDO
        ENDDO
      ELSE
        DO JLEV=1,KFLEV
          DO JROF=KST,KPROF
            PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1)=&
             & PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1)&
             & + PDT*PTENDGFL(JROF,JLEV,YCOMP(JGFL)%MP1)
          ENDDO
        ENDDO
      ENDIF
    ENDIF 
  ENDIF
ENDDO

!-----------------------------------------------------------------------

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('CPUTQY_APLPAR_LOOP',1,ZHOOK_HANDLE)

END SUBROUTINE CPUTQY_APLPAR_LOOP
