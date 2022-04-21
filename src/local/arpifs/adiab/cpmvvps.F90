!OCL  NOEVAL
SUBROUTINE CPMVVPS(YDVAB,KPROMA,KSTART,KPROF,KFLEV,PDT,&
 !  Variables 2D Input
 & PFP,&
 !  Variables 1D Input
 & PRES, PFEVL, PFEVN,&
 !  Variables 2D In/Out
 & PEVEL, PSDVBC, PSPT1)

!     ------------------------------------------------------------------
!     MODIFICATION DE LA VITESSE VERTICALE ET DE LA PRESSION DE SURFACE
!     DANS LE CAS NDPSFI = 1.
!     ------------------------------------------------------------------

!      VOIR DOCUMENTATION, INTERFACE PHYSICO-DYNAMIQUE
!                            -----------------

!     ARGUMENTS D ENTREE.
!     ------------------.
!       KPROMA : DIMENSION HORIZONTALE.
!       KSTART : DEBUT DE LA BOUCLE HORIZONTALE.
!       KPROF : BORNE HORIZONTALE DES CALCULS.
!       KFLEV : DIMENSION ET BORNE VERTICALE.
!       PDT   : Delta t (SL2TL or first timestep) or 2 Delta t.

! --- INPUT 2D.
!     --------.
!       PFP : FLUX TOTAL DE PRECIPITATIONS LIQUIDES ET NEIGEUSES.

! --- INPUT 1D.
!     --------.
!       PRES : PRESSION DE SURFACE A L'INSTANT OU EST CALCULEE LA PHYSIQUE.
!              Surface pressure at the same instant as non lagged physics.
!       PFEVL, PFEVN ( IDEM) : FLUX D'EVAPORATION.

!     ARGUMENTS IMPLICITES
!     --------------------
!       CONSTANTES UNIVERSELLES = COMMON /YOMCST/: RG
!       DECOUPAGE VERTICAL      = YOMGEM: YRVAB%YRVAB%VBH

!     SORTIES
!     -------
!       PEVEL (KPROMA,0:KFLEV) : VITESSE VERTICALE GENERALISEE.
!       PSDVBC (IDEM) : Integral of divergence term, including
!                       the "lrubc" and "delta m=1" contributions
!                       of (etadot DP/DETA), but not the
!                       "delta m=1" physics.
!       PSPT1 (KPROMA) : surface pressure or log of pressure buffer.

!     AUTEUR : E.BAZILE  JUIN 93.
!     ------

!     INSPIRE DE CPATY ECRIT EN SON TEMPS PAR A.JOLY

!     Modifications:
!     --------------
!      K. Yessad (Dec 2008): remove dummy CDLOCK and useless dummy arg
!      K. Yessad (Jan 2011): remove useless calculations.
!      K. Yessad (Feb 2018): remove deep-layer formulations.
!    ------------------------------------------------------------------

USE YOMVERT  , ONLY : TVAB
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCST   , ONLY : RG

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVAB)        ,INTENT(IN)    :: YDVAB
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFP(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRES(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFEVL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFEVN(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVEL(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSDVBC(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPT1(KPROMA) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZFE(KPROMA,0:KFLEV)
REAL(KIND=JPRB) :: ZEVELS(KPROMA)
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPMVVPS',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!     1. - MODIFY :
!          * The divergence vertical integral term.
!          * etadot (d prehyd/d eta).
!     ------------------------------------------------------------------

! * ZFE : FLUX D EVAPORATION TOTAL
DO JROF = KSTART, KPROF
  ZFE(JROF,KFLEV)=PFEVL(JROF)+PFEVN(JROF)
  ZEVELS(JROF)=ZFE(JROF,KFLEV)+PFP(JROF,KFLEV)
ENDDO
ZFE(KSTART:KPROF,0:KFLEV-1)=0.0_JPRB

! * 1.1 Surface/bottom values:

DO JROF = KSTART, KPROF
  PEVEL(JROF,KFLEV) = PEVEL(JROF,KFLEV) + RG*ZEVELS(JROF)
  PSDVBC(JROF,KFLEV) = PSDVBC(JROF,KFLEV) + RG*ZEVELS(JROF)
ENDDO

! * 1.2 Other levels:

DO JLEV = 1, KFLEV-1
  DO JROF = KSTART, KPROF
    PEVEL(JROF,JLEV) = PEVEL(JROF,JLEV) + RG*( YDVAB%VBH(JLEV)* ZEVELS(JROF) )  
    PSDVBC(JROF,JLEV)= PSDVBC(JROF,JLEV) + RG*( YDVAB%VBH(JLEV)* ZEVELS(JROF) )
  ENDDO
ENDDO

!    ------------------------------------------------------------------
!     2. - ADD PHYSICS TO PSPT1.
!    ------------------------------------------------------------------

DO JROF = KSTART, KPROF
  PSPT1(JROF)=PSPT1(JROF)-PDT*RG*(PFP(JROF,KFLEV)+ZFE(JROF,KFLEV))/PRES(JROF)  
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPMVVPS',1,ZHOOK_HANDLE)
END SUBROUTINE CPMVVPS
