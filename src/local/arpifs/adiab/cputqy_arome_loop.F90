SUBROUTINE CPUTQY_AROME_LOOP(YDMODEL, YDDIMV, YDGMV, YDCPG_BNDS, YDCPG_OPTS, PDT, KPGFL, KPTR, &
& PTENDGFLR, PB1, PGMVT1, PGFLT1)

!     --------------------------------------------
!**** *CPUTQY_AROME_LOOP*  EVOLUTION DE T,  et des  QI et de la TKE
!     --------------------------------------------

!     ARGUMENTS D'ENTREE
!     ------------------
!       KPROMA : DIMENSION HORIZONTALE
!       KSTART : DEBUT DE LA BOUCLE HORIZONTALE
!       KPROF : BORNE HORIZONTALE DES CALCULS
!       KFLEV : DIMENSION ET BORNE VERTICALE
!       KFLDN, KFLDX : SECONDE DIMENSION de PEXT1 (/= selon SL ou pas)
!       PDT : PAS DE TEMPS EFFECTIF

!       PTTS (KPROMA,KFLEV) : TENDANCE DE TEMPERATURE 
!       PRS ( KPROMA,KFLEV,NRR ) : TENDANCE DES Qi

!       KPGFL : pointer for PB1
!       KPTR : pointer for PTENDGFLR

!     ARGUMENTS IMPLICITES
!     --------------------

!     SORTIES
!     -------
!       PTT1(NPROMA,KFLEV): VARIABLES A FAIRE EVOLUER
!       PNEBH (NVCLIN_AROME): tableau contenant les Qi a faire evoluer
!       PQV1 : Qv au temps 1
!       PQL1 : Ql au temps 1
!       PQI1 : Qi au temps 1

!     Modifications
!     -------------
!   03/03/03 Creation par Y. Seity et S. Malardel
!   04/07/27 Modif par G. Hello (to run Hydrostatic case)
!   05/09/28 Modif par Y. Seity (Add PTENDEXT)
!   25-Jan-08 K. Yessad: bug correction for EXT.
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      07/12/18 Modif par Y. Seity (Add Hail PQH1)
!   26-Apr-2010 Modif par Y.Bouteloup (Only one call to cputqy_arome in mf_phys)
!   29-Apr-2019 R. El Khatib Use PTENDGFLR + pointer adresses KPTR to avoid a big array copy
!   18-Feb-2020 Y. Seity bf for non advected prognostic GFL (IF ordering)
! ----------------------------------------------------------------------------

USE YOMDIMV   , ONLY : TDIMV
USE YOMGMV    , ONLY : TGMV
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK
USE YOMCT0    , ONLY : LNHDYN, LSLAG
USE TYPE_MODEL, ONLY : MODEL
USE CPG_OPTS_TYPE_MOD, ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE

! ----------------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL)       ,INTENT(IN)    :: YDMODEL
TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TGMV)        ,INTENT(IN)    :: YDGMV
TYPE(CPG_BNDS_TYPE),INTENT(IN)   :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE),INTENT(IN)   :: YDCPG_OPTS
INTEGER(KIND=JPIM),INTENT(IN)    :: KPGFL(YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
INTEGER(KIND=JPIM),INTENT(IN)    :: KPTR(YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDGFLR(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
REAL(KIND=JPRB)   ,TARGET, INTENT(INOUT) :: PB1(YDCPG_OPTS%KLON,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1) 
REAL(KIND=JPRB)   ,TARGET, INTENT(INOUT) :: PGMVT1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDGMV%YT1%NDIM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLT1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1)

#include "cp_ptrslb1.intfb.h"

INTEGER(KIND=JPIM) :: IMP1EXPLICIT (8)

INTEGER(KIND=JPIM) :: ISLB1U9  ,ISLB1V9  ,ISLB1T9  ,ISLB1GFL9, ISLB1VD9

INTEGER(KIND=JPIM) :: JLEV, JROF, JGFL, IPGFL
REAL(KIND=JPRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPUTQY_AROME_LOOP', 0, ZHOOK_HANDLE)
ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL, NDIM1=>YGFL%NDIM1, NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & NFLSA=>YDDIMV%NFLSA, YT1=>YDGMV%YT1, NFLDSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1)
! ----------------------------------------------------------------------------

IF (LSLAG) CALL CP_PTRSLB1(YDMODEL%YRML_DYN%YRDYN, YDMODEL%YRML_DYN%YRPTRSLB1, ISLB1U9, ISLB1V9, &
           & ISLB1T9, ISLB1VD9, ISLB1GFL9)

IMP1EXPLICIT = [YGFL%YL%MP1, YGFL%YI%MP1, YGFL%YS%MP1, YGFL%YR%MP1, &
              & YGFL%YG%MP1, YGFL%YH%MP1, YGFL%YTKE%MP1, YGFL%YQ%MP1]

! Increment PGFLT1 for non-advected GFL, PB1 for advected GFL.
DO JGFL=1,NUMFLDS
  IF (ALL (IMP1EXPLICIT /= YCOMP(JGFL)%MP1)) THEN
    IF (YCOMP(JGFL)%LT1 .AND. KPTR(YCOMP(JGFL)%MP1) > 0) THEN
      IF (LSLAG .AND. YCOMP(JGFL)%LADV) THEN
        IPGFL=KPGFL(YCOMP(JGFL)%MP1)
        DO JLEV=1,YDCPG_OPTS%KFLEVG
          DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            PB1(JROF,ISLB1GFL9+IPGFL+JLEV-NFLSA)=&
             & PB1(JROF,ISLB1GFL9+IPGFL+JLEV-NFLSA)&
             & + PDT*PTENDGFLR(JROF,JLEV,KPTR(YCOMP(JGFL)%MP1))
          ENDDO
        ENDDO
      ELSE
        DO JLEV=1,YDCPG_OPTS%KFLEVG
          DO JROF=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
            PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1)=&
             & PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1)&
             & + PDT*PTENDGFLR(JROF,JLEV,KPTR(YCOMP(JGFL)%MP1))
          ENDDO
        ENDDO
      ENDIF ! End of non-zero diabatic tendency for this GFL
    ENDIF
  ENDIF
ENDDO

! ----------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPUTQY_AROME_LOOP', 1, ZHOOK_HANDLE)
END SUBROUTINE CPUTQY_AROME_LOOP
