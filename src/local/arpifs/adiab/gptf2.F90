SUBROUTINE GPTF2(YDGEOMETRY,YDGMV,&
 ! --- INPUT ---------------------------------------------------------
 & YDML_GCONF,YDDYN,KST,KEN,LDFSTEP,&
 ! --- INPUT-OUTPUT --------------------------------------------------
 & PGMV,PGMVS,PGFL, &
 & P0SP  , P0SPL , P0SPM , P9SP  , P9SPL , P9SPM, &
 & P0DIV , P0NHX , P0SPD , P0SPDL, P0SPDM, P0SVD , P0SVDL, &
 & P0SVDM, P0T   , P0TL  , P0TM  , P0U   , P0V   , P9DIV,  &
 & P9NHX , P9SPD , P9SPDL, P9SPDM, P9SVD , P9SVDL, P9SVDM, &
 & P9T   , P9TL  , P9TM  , P9U   , P9V)

!**** *GPTF2* - Timefilter part 2

!     Purpose.
!     --------
!           Performs part 2 of the time-filtering.
!           - leap frog + ldfstep=true : pxt9. = pxt0.
!           - leap frog + ldfstep=false: pxt9. = pxt9. + eps2*pxt0.
!           - sl2tl     + ldfstep=true : put9=put0 and pvt9=pvt0 only.
!           - sl2tl     + ldfstep=false: nothing.

!**   Interface.
!     ----------
!        *CALL* *GPTF2(...)

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
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCT0       , ONLY : LTWOTL, LNHDYN
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : NVDVAR

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEN 
LOGICAL           ,INTENT(IN)    :: LDFSTEP 
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV) 
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB),POINTER,OPTIONAL,INTENT(INOUT),CONTIGUOUS,DIMENSION(:)   :: P0SP  , P0SPL , P0SPM , P9SP  , P9SPL , P9SPM 
REAL(KIND=JPRB),POINTER,OPTIONAL,INTENT(INOUT),CONTIGUOUS,DIMENSION(:,:) :: P0DIV , P0NHX , P0SPD , P0SPDL, P0SPDM, P0SVD , P0SVDL
REAL(KIND=JPRB),POINTER,OPTIONAL,INTENT(INOUT),CONTIGUOUS,DIMENSION(:,:) :: P0SVDM, P0T   , P0TL  , P0TM  , P0U   , P0V   , P9DIV 
REAL(KIND=JPRB),POINTER,OPTIONAL,INTENT(INOUT),CONTIGUOUS,DIMENSION(:,:) :: P9NHX , P9SPD , P9SPDL, P9SPDM, P9SVD , P9SVDL, P9SVDM
REAL(KIND=JPRB),POINTER,OPTIONAL,INTENT(INOUT),CONTIGUOUS,DIMENSION(:,:) :: P9T   , P9TL  , P9TM  , P9U   , P9V   

!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JL,JGFL,JK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:)   :: Z0SP  , Z0SPL , Z0SPM , Z9SP  , Z9SPL , Z9SPM 
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z0DIV , Z0NHX , Z0SPD , Z0SPDL, Z0SPDM, Z0SVD , Z0SVDL
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z0SVDM, Z0T   , Z0TL  , Z0TM  , Z0U   , Z0V   , Z9DIV 
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z9NHX , Z9SPD , Z9SPDL, Z9SPDM, Z9SVD , Z9SVDL, Z9SVDM
REAL(KIND=JPRB),POINTER,CONTIGUOUS,DIMENSION(:,:) :: Z9T   , Z9TL  , Z9TM  , Z9U   , Z9V   

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPTF2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YGFL=>YDML_GCONF%YGFL,YDDIMF=>YDML_GCONF%YRDIMF)
ASSOCIATE(NDIM=>YGFL%NDIM, NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFTHER=>YDDIMF%NFTHER, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & REPS2=>YDDYN%REPS2, REPSM2=>YDDYN%REPSM2, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0, &
 & YT9=>YDGMV%YT9)
!     ------------------------------------------------------------------

Z0SP    => NULL ()
Z0SPL   => NULL ()
Z0SPM   => NULL ()
Z9SP    => NULL ()
Z9SPL   => NULL ()
Z9SPM   => NULL ()

Z0DIV   => NULL ()
Z0NHX   => NULL ()
Z0SPD   => NULL ()
Z0SPDL  => NULL ()
Z0SPDM  => NULL ()
Z0SVD   => NULL ()
Z0SVDL  => NULL ()
Z0SVDM  => NULL ()
Z0T     => NULL ()
Z0TL    => NULL ()
Z0TM    => NULL ()
Z0U     => NULL ()
Z0V     => NULL ()
Z9DIV   => NULL ()
Z9NHX   => NULL ()
Z9SPD   => NULL ()
Z9SPDL  => NULL ()
Z9SPDM  => NULL ()
Z9SVD   => NULL ()
Z9SVDL  => NULL ()
Z9SVDM  => NULL ()
Z9T     => NULL ()
Z9TL    => NULL ()
Z9TM    => NULL ()
Z9U     => NULL ()
Z9V     => NULL ()
           
IF (PRESENT (PGMVS)) THEN           
  Z0SP  => PGMVS(:,YT0%MSP)                       
  Z0SPL => PGMVS(:,YT0%MSPL)                       
  Z0SPM => PGMVS(:,YT0%MSPM)                       
  Z9SP  => PGMVS(:,YT9%MSP)                       
  Z9SPL => PGMVS(:,YT9%MSPL)                       
  Z9SPM => PGMVS(:,YT9%MSPM)                       
ELSEIF (PRESENT (P0SP) .AND.  PRESENT (P0SPL) .AND. PRESENT (P0SPM) .AND. &
      & PRESENT (P9SP ) .AND.  PRESENT (P9SPL) .AND.  PRESENT (P9SPM)) THEN
  Z0SP  => P0SP  
  Z0SPL => P0SPL 
  Z0SPM => P0SPM 
  Z9SP  => P9SP  
  Z9SPL => P9SPL 
  Z9SPM => P9SPM 
ELSE
  CALL ABOR1 ('GPTF2: PGMVS OR P0SP/P0SPL/P0SPM/P9SP/P9SPL/P9SPM REQUIRED')
ENDIF

IF (PRESENT (PGMV)) THEN
  Z0DIV  => PGMV(:,:,YT0%MDIV)                       
  Z0NHX  => PGMV(:,:,YT0%MNHX)                       
  Z0SPD  => PGMV(:,:,YT0%MSPD)                       
  Z0SPDL => PGMV(:,:,YT0%MSPDL)                       
  Z0SPDM => PGMV(:,:,YT0%MSPDM)                       
  Z0SVD  => PGMV(:,:,YT0%MSVD)                       
  Z0SVDL => PGMV(:,:,YT0%MSVDL)                       
  Z0SVDM => PGMV(:,:,YT0%MSVDM)                       
  Z0T    => PGMV(:,:,YT0%MT)                       
  Z0TL   => PGMV(:,:,YT0%MTL)                       
  Z0TM   => PGMV(:,:,YT0%MTM)                       
  Z0U    => PGMV(:,:,YT0%MU)                       
  Z0V    => PGMV(:,:,YT0%MV)                       
  Z9DIV  => PGMV(:,:,YT9%MDIV)                       
  Z9NHX  => PGMV(:,:,YT9%MNHX)                       
  Z9SPD  => PGMV(:,:,YT9%MSPD)                       
  Z9SPDL => PGMV(:,:,YT9%MSPDL)                       
  Z9SPDM => PGMV(:,:,YT9%MSPDM)                       
  Z9SVD  => PGMV(:,:,YT9%MSVD)                       
  Z9SVDL => PGMV(:,:,YT9%MSVDL)                       
  Z9SVDM => PGMV(:,:,YT9%MSVDM)                       
  Z9T    => PGMV(:,:,YT9%MT)                       
  Z9TL   => PGMV(:,:,YT9%MTL)                       
  Z9TM   => PGMV(:,:,YT9%MTM)                       
  Z9U    => PGMV(:,:,YT9%MU)                       
  Z9V    => PGMV(:,:,YT9%MV)                       
ELSEIF (PRESENT (P0DIV ) .AND.  PRESENT (P0NHX ) .AND.  PRESENT (P0SPD ) .AND.  PRESENT (P0SPDL) .AND. &
      & PRESENT (P0SPDM) .AND.  PRESENT (P0SVD ) .AND.  PRESENT (P0SVDL) .AND.  PRESENT (P0SVDM) .AND. &
      & PRESENT (P0T   ) .AND.  PRESENT (P0TL  ) .AND.  PRESENT (P0TM  ) .AND.  PRESENT (P0U   ) .AND. &
      & PRESENT (P0V   ) .AND.  PRESENT (P9DIV ) .AND.  PRESENT (P9NHX ) .AND.  PRESENT (P9SPD ) .AND. &
      & PRESENT (P9SPDL) .AND.  PRESENT (P9SPDM) .AND.  PRESENT (P9SVD ) .AND.  PRESENT (P9SVDL) .AND. &
      & PRESENT (P9SVDM) .AND.  PRESENT (P9T   ) .AND.  PRESENT (P9TL  ) .AND.  PRESENT (P9TM  ) .AND. &
      & PRESENT (P9U   ) .AND.  PRESENT (P9V   )) THEN
  Z0DIV  => P0DIV  
  Z0NHX  => P0NHX  
  Z0SPD  => P0SPD  
  Z0SPDL => P0SPDL 
  Z0SPDM => P0SPDM 
  Z0SVD  => P0SVD  
  Z0SVDL => P0SVDL 
  Z0SVDM => P0SVDM 
  Z0T    => P0T    
  Z0TL   => P0TL   
  Z0TM   => P0TM   
  Z0U    => P0U    
  Z0V    => P0V    
  Z9DIV  => P9DIV  
  Z9NHX  => P9NHX  
  Z9SPD  => P9SPD  
  Z9SPDL => P9SPDL 
  Z9SPDM => P9SPDM 
  Z9SVD  => P9SVD  
  Z9SVDL => P9SVDL 
  Z9SVDM => P9SVDM 
  Z9T    => P9T    
  Z9TL   => P9TL   
  Z9TM   => P9TM   
  Z9U    => P9U    
  Z9V    => P9V    
ELSE
  WRITE (0, *) " P0DIV  = ", PRESENT (P0DIV ) 
  WRITE (0, *) " P0NHX  = ", PRESENT (P0NHX ) 
  WRITE (0, *) " P0SPD  = ", PRESENT (P0SPD ) 
  WRITE (0, *) " P0SPDL = ", PRESENT (P0SPDL) 
  WRITE (0, *) " P0SPDM = ", PRESENT (P0SPDM) 
  WRITE (0, *) " P0SVD  = ", PRESENT (P0SVD ) 
  WRITE (0, *) " P0SVDL = ", PRESENT (P0SVDL) 
  WRITE (0, *) " P0SVDM = ", PRESENT (P0SVDM) 
  WRITE (0, *) " P0T    = ", PRESENT (P0T   ) 
  WRITE (0, *) " P0TL   = ", PRESENT (P0TL  ) 
  WRITE (0, *) " P0TM   = ", PRESENT (P0TM  ) 
  WRITE (0, *) " P0U    = ", PRESENT (P0U   ) 
  WRITE (0, *) " P0V    = ", PRESENT (P0V   ) 
  WRITE (0, *) " P9DIV  = ", PRESENT (P9DIV ) 
  WRITE (0, *) " P9NHX  = ", PRESENT (P9NHX ) 
  WRITE (0, *) " P9SPD  = ", PRESENT (P9SPD ) 
  WRITE (0, *) " P9SPDL = ", PRESENT (P9SPDL) 
  WRITE (0, *) " P9SPDM = ", PRESENT (P9SPDM) 
  WRITE (0, *) " P9SVD  = ", PRESENT (P9SVD ) 
  WRITE (0, *) " P9SVDL = ", PRESENT (P9SVDL) 
  WRITE (0, *) " P9SVDM = ", PRESENT (P9SVDM) 
  WRITE (0, *) " P9T    = ", PRESENT (P9T   ) 
  WRITE (0, *) " P9TL   = ", PRESENT (P9TL  ) 
  WRITE (0, *) " P9TM   = ", PRESENT (P9TM  ) 
  WRITE (0, *) " P9U    = ", PRESENT (P9U   ) 
  WRITE (0, *) " P9V    = ", PRESENT (P9V   )
  CALL FLUSH (0)
  CALL ABOR1 ('GPTF2: PGMV OR P0DIV/P0NHX/P0SPD/P0SPDL/P0SPDM/P0SVD/P0SVDL/P0SVDM/P0T/P0TL/P0TM/P0U/P0V/'//&
            & 'P9DIV/P9NHX/P9SPD/P9SPDL/P9SPDM/P9SVD/P9SVDL/P9SVDM/P9T/P9TL/P9TM/P9U/P9V REQUIRED')
ENDIF

           
           
           
!*       1. PERFORM TIME FILTER (PART 2)           
!           ----------------------------           
           
!*        1.1   TWO TIME LEVEL           
           
IF (LTWOTL) THEN

  IF(LDFSTEP) THEN
    DO JK=1,NFLEVG
      DO JL=KST,KEN
        Z9U(JL,JK) = Z0U(JL,JK)
        Z9V(JL,JK) = Z0V(JL,JK)
      ENDDO
    ENDDO
  ENDIF

ELSE

!*        1.2   THREE TIME LEVEL

  IF(LDFSTEP) THEN
    DO JK=1,NFLEVG
      DO JL=KST,KEN
        Z9U(JL,JK)   = Z0U(JL,JK)
        Z9V(JL,JK)   = Z0V(JL,JK)
        Z9DIV(JL,JK) = Z0DIV(JL,JK)
      ENDDO
      IF(NFTHER >= 1) THEN
        DO JL=KST,KEN
          Z9T(JL,JK)  = Z0T(JL,JK)
          Z9TL(JL,JK) = Z0TL(JL,JK)
          Z9TM(JL,JK) = Z0TM(JL,JK)
        ENDDO
      ENDIF
    ENDDO

    DO JGFL=1,NUMFLDS
      IF(YCOMP(JGFL)%LT9) THEN
        DO JK=1,NFLEVG
          DO JL=KST,KEN
            PGFL(JL,JK,YCOMP(JGFL)%MP9) = PGFL(JL,JK,YCOMP(JGFL)%MP)
          ENDDO
        ENDDO
      ENDIF
    ENDDO

    IF(LNHDYN) THEN
      DO JK=1,NFLEVG
        DO JL=KST,KEN
          Z9SPD(JL,JK) = Z0SPD(JL,JK)
          Z9SVD(JL,JK) = Z0SVD(JL,JK)
          Z9SPDL(JL,JK) = Z0SPDL(JL,JK)
          Z9SVDL(JL,JK) = Z0SVDL(JL,JK)
          Z9SPDM(JL,JK) = Z0SPDM(JL,JK)
          Z9SVDM(JL,JK) = Z0SVDM(JL,JK)
        ENDDO
      ENDDO

      IF( NVDVAR==4 .OR. NVDVAR==5 ) THEN
        DO JK=1,NFLEVG
          DO JL=KST,KEN
            Z9NHX(JL,JK) = Z0NHX(JL,JK)
          ENDDO
        ENDDO
      ENDIF

    ENDIF

    DO JL=KST,KEN
      Z9SP(JL)  = Z0SP(JL)
      Z9SPL(JL) = Z0SPL(JL)
      Z9SPM(JL) = Z0SPM(JL)
    ENDDO

  ELSE

    DO JK=1,NFLEVG
      DO JL=KST,KEN
        Z9U(JL,JK)   = Z9U(JL,JK)  +REPS2*Z0U(JL,JK)
        Z9V(JL,JK)   = Z9V(JL,JK)  +REPS2*Z0V(JL,JK)
        Z9DIV(JL,JK) = Z9DIV(JL,JK)+REPS2*Z0DIV(JL,JK)
      ENDDO
      IF(NFTHER >= 1) THEN
        DO JL=KST,KEN
          Z9T(JL,JK)  = Z9T(JL,JK)  +REPS2*Z0T(JL,JK)
          Z9TL(JL,JK) = Z9TL(JL,JK) +REPS2*Z0TL(JL,JK)
          Z9TM(JL,JK) = Z9TM(JL,JK) +REPS2*Z0TM(JL,JK)
        ENDDO
      ENDIF
    ENDDO

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

    IF(LNHDYN) THEN
      DO JK=1,NFLEVG
        DO JL=KST,KEN
          Z9SPD(JL,JK)  =&
           & Z9SPD(JL,JK)+REPS2*Z0SPD(JL,JK)
          Z9SVD(JL,JK)  =&
           & Z9SVD(JL,JK)+REPS2*Z0SVD(JL,JK)
          Z9SPDL(JL,JK) =&
           & Z9SPDL(JL,JK)+REPS2*Z0SPDL(JL,JK)
          Z9SVDL(JL,JK) =&
           & Z9SVDL(JL,JK)+REPS2*Z0SVDL(JL,JK)
          Z9SPDM(JL,JK) =&
           & Z9SPDM(JL,JK)+REPS2*Z0SPDM(JL,JK)
          Z9SVDM(JL,JK) =&
           & Z9SVDM(JL,JK)+REPS2*Z0SVDM(JL,JK)
        ENDDO
      ENDDO
    ENDIF

    DO JL=KST,KEN
      Z9SP(JL)  = Z9SP(JL) +REPS2*Z0SP(JL)
      Z9SPL(JL) = Z9SPL(JL)+REPS2*Z0SPL(JL)
      Z9SPM(JL) = Z9SPM(JL)+REPS2*Z0SPM(JL)
    ENDDO

  ENDIF

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPTF2',1,ZHOOK_HANDLE)
END SUBROUTINE GPTF2
