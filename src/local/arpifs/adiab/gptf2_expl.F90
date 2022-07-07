SUBROUTINE GPTF2_EXPL(YDGEOMETRY,&
 ! --- INPUT ---------------------------------------------------------
 & YDML_GCONF,YDDYN,KST,KEN,LDFSTEP,&
 ! --- INPUT-OUTPUT --------------------------------------------------
 & PGFL, &
 & P0SP  , P0SPL , P0SPM , P9SP  , P9SPL , P9SPM, &
 & P0DIV , P0NHX , P0SPD , P0SPDL, P0SPDM, P0SVD , P0SVDL, &
 & P0SVDM, P0T   , P0TL  , P0TM  , P0U   , P0V   , P9DIV,  &
 & P9NHX , P9SPD , P9SPDL, P9SPDM, P9SVD , P9SVDL, P9SVDM, &
 & P9T   , P9TL  , P9TM  , P9U   , P9V)

!**** *GPTF2_EXPL* - Timefilter part 2

!     Purpose.
!     --------
!           Performs part 2 of the time-filtering.
!           - leap frog + ldfstep=true : pxt9. = pxt0.
!           - leap frog + ldfstep=false: pxt9. = pxt9. + eps2*pxt0.
!           - sl2tl     + ldfstep=true : put9=put0 and pvt9=pvt0 only.
!           - sl2tl     + ldfstep=false: nothing.

!**   Interface.
!     ----------
!        *CALL* *GPTF2_EXPL(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!         KST          : start of work
!         KEN          : depth of work
!         LDLSTEP      : check on the last time step.

!        INPUT/OUTPUT:
!         PGMV         : "t" and "t-dt" upper air GMV variables.
!         PGMVS        : "t" and "t-dt" surface GMV variables.
!         PGFL         : "t" and "t-dt" GFL variables.

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud  *ECMWF*
!      Original : 92-02-01

! Modifications.
! --------------
!   Modified 02-09-30 by P. Smolikova (variable d4 in NH)
!   Modified 13-11-02 K. YESSAD : cleanings + improve vectorization.
!   Modified 2003-07-17 C. Fischer - psvdauxt0* come from pgmv
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!   08-Jun-2004 J. Masek   NH cleaning (LVSLWBC)
!   Modified Nov 2007 N. Wedi: bug correction
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2008): remove dummy CDLOCK
! End Modifications
!------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCT0       , ONLY : LTWOTL, LNHDYN
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : NVDVAR

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEN 
LOGICAL           ,INTENT(IN)    :: LDFSTEP 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB),OPTIONAL, INTENT(INOUT),DIMENSION(YDGEOMETRY%YRDIM%NPROMA)   :: P0SP  , P0SPL , P0SPM , P9SP  , P9SPL , P9SPM 
REAL(KIND=JPRB),OPTIONAL, INTENT(INOUT),DIMENSION(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) :: P0DIV , P0NHX , P0SPD , P0SPDL, P0SPDM, P0SVD , P0SVDL
REAL(KIND=JPRB),OPTIONAL, INTENT(INOUT),DIMENSION(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) :: P0SVDM, P0T   , P0TL  , P0TM  , P0U   , P0V   , P9DIV 
REAL(KIND=JPRB),OPTIONAL, INTENT(INOUT),DIMENSION(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) :: P9NHX , P9SPD , P9SPDL, P9SPDM, P9SVD , P9SVDL, P9SVDM
REAL(KIND=JPRB),OPTIONAL, INTENT(INOUT),DIMENSION(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) :: P9T   , P9TL  , P9TM  , P9U   , P9V   

#include "gptf2_expl_2tl.intfb.h"

!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JL,JGFL,JK
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPTF2_EXPL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YGFL=>YDML_GCONF%YGFL,YDDIMF=>YDML_GCONF%YRDIMF)
ASSOCIATE(NDIM=>YGFL%NDIM, NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFTHER=>YDDIMF%NFTHER, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & REPS2=>YDDYN%REPS2, REPSM2=>YDDYN%REPSM2)
!     ------------------------------------------------------------------


!*       1. PERFORM TIME FILTER (PART 2)           
!           ----------------------------           
           
!*        1.1   TWO TIME LEVEL           
           
IF (LTWOTL) THEN

  CALL GPTF2_EXPL_2TL (YDGEOMETRY, KST, KEN, LDFSTEP, P0U, P0V, P9U, P9V)

ELSE

!*        1.2   THREE TIME LEVEL

  IF(LDFSTEP) THEN
    DO JK=1,NFLEVG
      DO JL=KST,KEN
        P9U(JL,JK)   = P0U(JL,JK)
        P9V(JL,JK)   = P0V(JL,JK)
        P9DIV(JL,JK) = P0DIV(JL,JK)
      ENDDO
      IF(NFTHER >= 1) THEN
        DO JL=KST,KEN
          P9T(JL,JK)  = P0T(JL,JK)
          P9TL(JL,JK) = P0TL(JL,JK)
          P9TM(JL,JK) = P0TM(JL,JK)
        ENDDO
      ENDIF
    ENDDO

    IF(LNHDYN) THEN
      DO JK=1,NFLEVG
        DO JL=KST,KEN
          P9SPD(JL,JK) = P0SPD(JL,JK)
          P9SVD(JL,JK) = P0SVD(JL,JK)
          P9SPDL(JL,JK) = P0SPDL(JL,JK)
          P9SVDL(JL,JK) = P0SVDL(JL,JK)
          P9SPDM(JL,JK) = P0SPDM(JL,JK)
          P9SVDM(JL,JK) = P0SVDM(JL,JK)
        ENDDO
      ENDDO

      IF( NVDVAR==4 .OR. NVDVAR==5 ) THEN
        DO JK=1,NFLEVG
          DO JL=KST,KEN
            P9NHX(JL,JK) = P0NHX(JL,JK)
          ENDDO
        ENDDO
      ENDIF

    ENDIF

    DO JL=KST,KEN
      P9SP(JL)  = P0SP(JL)
      P9SPL(JL) = P0SPL(JL)
      P9SPM(JL) = P0SPM(JL)
    ENDDO

  ELSE

    DO JK=1,NFLEVG
      DO JL=KST,KEN
        P9U(JL,JK)   = P9U(JL,JK)  +REPS2*P0U(JL,JK)
        P9V(JL,JK)   = P9V(JL,JK)  +REPS2*P0V(JL,JK)
        P9DIV(JL,JK) = P9DIV(JL,JK)+REPS2*P0DIV(JL,JK)
      ENDDO
      IF(NFTHER >= 1) THEN
        DO JL=KST,KEN
          P9T(JL,JK)  = P9T(JL,JK)  +REPS2*P0T(JL,JK)
          P9TL(JL,JK) = P9TL(JL,JK) +REPS2*P0TL(JL,JK)
          P9TM(JL,JK) = P9TM(JL,JK) +REPS2*P0TM(JL,JK)
        ENDDO
      ENDIF
    ENDDO

    IF(LNHDYN) THEN
      DO JK=1,NFLEVG
        DO JL=KST,KEN
          P9SPD(JL,JK)  =&
           & P9SPD(JL,JK)+REPS2*P0SPD(JL,JK)
          P9SVD(JL,JK)  =&
           & P9SVD(JL,JK)+REPS2*P0SVD(JL,JK)
          P9SPDL(JL,JK) =&
           & P9SPDL(JL,JK)+REPS2*P0SPDL(JL,JK)
          P9SVDL(JL,JK) =&
           & P9SVDL(JL,JK)+REPS2*P0SVDL(JL,JK)
          P9SPDM(JL,JK) =&
           & P9SPDM(JL,JK)+REPS2*P0SPDM(JL,JK)
          P9SVDM(JL,JK) =&
           & P9SVDM(JL,JK)+REPS2*P0SVDM(JL,JK)
        ENDDO
      ENDDO
    ENDIF

    DO JL=KST,KEN
      P9SP(JL)  = P9SP(JL) +REPS2*P0SP(JL)
      P9SPL(JL) = P9SPL(JL)+REPS2*P0SPL(JL)
      P9SPM(JL) = P9SPM(JL)+REPS2*P0SPM(JL)
    ENDDO

  ENDIF

  IF (LDFSTEP) THEN

    DO JGFL=1,NUMFLDS
      IF(YCOMP(JGFL)%LT9) THEN
        DO JK=1,NFLEVG
          DO JL=KST,KEN
            PGFL(JL,JK,YCOMP(JGFL)%MP9) = PGFL(JL,JK,YCOMP(JGFL)%MP)
          ENDDO
        ENDDO
      ENDIF
    ENDDO

  ELSE

    DO JGFL=1,NUMFLDS
      IF(YCOMP(JGFL)%LT9) THEN
        DO JK=1,NFLEVG
          DO JL=KST,KEN
            PGFL(JL,JK,YCOMP(JGFL)%MP9) = PGFL(JL,JK,YCOMP(JGFL)%MP9)+&
             & REPSM2*PGFL(JL,JK,YCOMP(JGFL)%MP)  
          ENDDO
        ENDDO
      ENDIF
    ENDDO

  ENDIF
  

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPTF2_EXPL',1,ZHOOK_HANDLE)
END SUBROUTINE GPTF2_EXPL
