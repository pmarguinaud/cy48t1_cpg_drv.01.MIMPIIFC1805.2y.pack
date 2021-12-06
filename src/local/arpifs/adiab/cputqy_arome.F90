SUBROUTINE CPUTQY_AROME(YDMODEL, YDDIMV,YDGMV, KPROMA, KSTART, KPROF, KFLEV, PDT, KPGFL, KPTR, &
 !  Variables 2D Input
 & PTENDT, PTENDGFLR, PTENDU, PTENDV, PTENDD,&
 !  Variables 2D In/Out
 & PB1, PGMVT1, PGFLT1)

!     --------------------------------------------
!**** *CPUTQY_AROME*  EVOLUTION DE T,  et des  QI et de la TKE
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

! ----------------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPGFL(YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
INTEGER(KIND=JPIM),INTENT(IN)    :: KPTR(YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDT(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDGFLR(KPROMA,KFLEV,YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDU(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDV(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTENDD(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,TARGET, INTENT(INOUT) :: PB1(KPROMA,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1) 
REAL(KIND=JPRB)   ,TARGET, INTENT(INOUT) :: PGMVT1(KPROMA,KFLEV,YDGMV%YT1%NDIM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLT1(KPROMA,KFLEV,YDMODEL%YRML_GCONF%YGFL%NDIM1)

#include "cp_ptrslb1.intfb.h"

! ----------------------------------------------------------------------------
!  Local pointers LSLAG choice
REAL(KIND=JPRB),DIMENSION(:,:),    POINTER :: ZTT1
REAL(KIND=JPRB),DIMENSION(:,:),    POINTER :: ZUT1
REAL(KIND=JPRB),DIMENSION(:,:),    POINTER :: ZVT1
REAL(KIND=JPRB),DIMENSION(:,:),    POINTER :: ZDT1

INTEGER(KIND=JPIM) :: ISLB1U9  ,ISLB1V9  ,ISLB1T9  ,ISLB1GFL9, ISLB1VD9

INTEGER(KIND=JPIM) :: JLEV, JROF, JGFL, IPGFL
REAL(KIND=JPRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPUTQY_AROME',0,ZHOOK_HANDLE)
ASSOCIATE(NDIM1=>YDMODEL%YRML_GCONF%YGFL%NDIM1, NUMFLDS=>YDMODEL%YRML_GCONF%YGFL%NUMFLDS, YCOMP=>YDMODEL%YRML_GCONF%YGFL%YCOMP, &
 & NFLSA=>YDDIMV%NFLSA, &
 & YT1=>YDGMV%YT1, &
 & NFLDSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1)
! ----------------------------------------------------------------------------

IF (LSLAG) CALL CP_PTRSLB1(YDMODEL%YRML_DYN%YRDYN,YDMODEL%YRML_DYN%YRPTRSLB1,ISLB1U9,ISLB1V9,ISLB1T9,ISLB1VD9,ISLB1GFL9)

!  ------------------------------------------------------------
!   1. CALCUL DE L'EVOLUTION DES GMV (TEMPERATURE,VENT, et w)
!  ------------------------------------------------------------

! -------------------------------------------------------------------
! Define local pointers LSLAG choice
!  GMV variables
IF (LSLAG) THEN
  ZTT1 => PB1(:,ISLB1T9+1-NFLSA:)
  ZUT1 => PB1(:,ISLB1U9+1-NFLSA:)
  ZVT1 => PB1(:,ISLB1V9+1-NFLSA:)
  ZDT1 => PB1(:,ISLB1VD9+1-NFLSA:)
ELSE
  ZTT1 => PGMVT1(:,:,YT1%MT)
  ZUT1 => PGMVT1(:,:,YT1%MU)
  ZVT1 => PGMVT1(:,:,YT1%MV)
  ZDT1 => PGMVT1(:,:,YT1%MSVD)
ENDIF  

DO JLEV=1,KFLEV
  DO JROF=KSTART,KPROF
    ZTT1(JROF,JLEV)=ZTT1(JROF,JLEV) + PDT*PTENDT(JROF,JLEV)
    ZUT1(JROF,JLEV)=ZUT1(JROF,JLEV) + PDT*PTENDU(JROF,JLEV)
    ZVT1(JROF,JLEV)=ZVT1(JROF,JLEV) + PDT*PTENDV(JROF,JLEV)
    IF(LNHDYN) THEN
       ZDT1(JROF,JLEV)=ZDT1(JROF,JLEV) + PTENDD(JROF,JLEV)*PDT  
    ENDIF     
  ENDDO
ENDDO

!  ------------------------------------------------------------
!   2. CALCUL DE L'EVOLUTION DES GFL ( HYDROMETEORES et EXTRA-GFL ) 
!  ------------------------------------------------------------

! should be identical to what is in CPUTQY.

IF (LSLAG) THEN
  ! Increment PGFLT1 for non-advected GFL, PB1 for advected GFL.
  DO JGFL=1,NUMFLDS
    IF (YCOMP(JGFL)%LT1) THEN
      IF (KPTR(YCOMP(JGFL)%MP1) > 0) THEN
        IF (YCOMP(JGFL)%LADV) THEN
          IPGFL=KPGFL(YCOMP(JGFL)%MP1)
          DO JLEV=1,KFLEV
            DO JROF=KSTART,KPROF
              PB1(JROF,ISLB1GFL9+IPGFL+JLEV-NFLSA)=&
               & PB1(JROF,ISLB1GFL9+IPGFL+JLEV-NFLSA)&
               & + PDT*PTENDGFLR(JROF,JLEV,KPTR(YCOMP(JGFL)%MP1))
            ENDDO
          ENDDO
        ELSE
          DO JLEV=1,KFLEV
            DO JROF=KSTART,KPROF
              PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1)=&
               & PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1)&
               & + PDT*PTENDGFLR(JROF,JLEV,KPTR(YCOMP(JGFL)%MP1))
            ENDDO
          ENDDO
        ENDIF ! End of non-zero diabatic tendency for this GFL
      ENDIF
    ENDIF
  ENDDO
ELSE
  ! Increment PGFLT1.
  DO JGFL=1,NUMFLDS
    IF (YCOMP(JGFL)%LT1) THEN
      IF (KPTR(YCOMP(JGFL)%MP1) > 0) THEN
        DO JLEV=1,KFLEV
          DO JROF=KSTART,KPROF
            PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1)=&
             & PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1)&
             & + PDT*PTENDGFLR(JROF,JLEV,KPTR(YCOMP(JGFL)%MP1))
          ENDDO
        ENDDO
      ENDIF  ! End of non-zero diabatic tendency for this GFL
    ENDIF
  ENDDO
ENDIF

! ----------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPUTQY_AROME',1,ZHOOK_HANDLE)
END SUBROUTINE CPUTQY_AROME
