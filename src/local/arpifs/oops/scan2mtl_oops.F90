SUBROUTINE SCAN2MTL_OOPS(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,CDCONF,LD_TST_GPGFL,LD_DFISTEP,PTRAJEC_OOPS,YDPHYSMWAVE5)

!**** *SCAN2MTL* - Interface to computations in grid-point space

!     Purpose. Interface to computations in grid-point space
!     --------

!     Interface.
!     ----------
!        *CALL* *SCAN2MTL(....)

!        Explicit arguments :
!        --------------------  
!        CDCONF       : configuration of work (see doc.)
!        LD_TST_GPGFL : if T, do timestepping on grid-point GFL.

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   See includes below.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud   *ECMWF*
!      Original 91-04-01 Mats Hamrud

!     Modifications.
!     --------------
!      O. Marsden June 2016 : same as SCAN2MTL with GOM removed, could be remerged
!                             with SCAN2MTL when that has been cleaned up
!    -------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE FIELDS_MOD         , ONLY : FIELDS
USE MTRAJ_MOD          , ONLY : MTRAJ
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE YOMCT0             , ONLY : LELAM, LECMWF
USE YOMCT3             , ONLY : NSTEP
USE YOMVRTL            , ONLY : LIDMODEL
USE YOMLCZ             , ONLY : LFORCEWR, GPFORCEU, GPFORCEV, GPFORCET, GPFORCEQ, GPFORCESP
USE YOMTRAJ_OOPS       , ONLY : TRAJ_TYPE_OOPS
USE RANDOM_NUMBERS_MIX , ONLY : GAUSSIAN_DISTRIBUTION
USE YOE_PHYS_MWAVE     , ONLY : TEPHYSMWAVE

!    -------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)       , INTENT(IN)    :: YDGEOMETRY
TYPE(FIELDS)         , INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ)          , INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL)           ,INTENT(INOUT) :: YDMODEL
CHARACTER(LEN=9)     , INTENT(IN)    :: CDCONF
LOGICAL              , INTENT(IN)    :: LD_TST_GPGFL
LOGICAL              , INTENT(IN)    :: LD_DFISTEP
TYPE(TRAJ_TYPE_OOPS) , INTENT(INOUT) :: PTRAJEC_OOPS
TYPE(TEPHYSMWAVE)    , INTENT(INOUT) :: YDPHYSMWAVE5
!    -------------------------------------------------------------------
REAL(KIND=JPRB), ALLOCATABLE :: ZPERTTSG(:), ZPERTTS(:)
INTEGER(KIND=JPIM) :: IBL, ICEND, ICSTA, IOFF, JKGLO, JLEV, JL,&
                    & JGFL, IGPCOMP

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!    -------------------------------------------------------------------

#include "gp_model_tl.intfb.h"
#include "reset_spert.intfb.h"
#include "lcnorggtl.intfb.h"

!    -------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SCAN2MTL_OOPS',0,ZHOOK_HANDLE)
ASSOCIATE(YDGFL5=>YDMTRAJ%YRGFL5,YDGMV5=>YDMTRAJ%YRGMV5, YDGFL=>YDFIELDS%YRGFL,YDGMV=>YDFIELDS%YRGMV,                               &
& YDSURF=>YDFIELDS%YRSURF, YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM,                             &
& YDGSGEOM=>YDGEOMETRY%YRGSGEOM, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, LEGBRAD=>YDMODEL%YRML_PHY_EC%YREPHY%LEGBRAD,                  &
& SD_VN=>YDFIELDS%YRSURF%SD_VN, YSD_VN=>YDFIELDS%YRSURF%YSD_VN,   YGFL=>YDMODEL%YRML_GCONF%YGFL, YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,&
& YDRIP=>YDMODEL%YRML_GCONF%YRRIP)

!    -------------------------------------------------------------------

!    -------------------------------------------------------------------
!*       1.    INITIAL SETUP, MAIN LOOP CONTROL.
!    -------------------------------------------------------------------

! Total number of grid points - ndglg*ndlon in global / ndguxg*ndlon in LAM
! used for model computations
IF (LELAM) THEN
  IGPCOMP=MIN(YDGEM%NGPTOT,YDGEM%NGPTOT_CAP)
ELSE
  IGPCOMP=YDGEM%NGPTOT
ENDIF


!    -------------------------------------------------------------------
!*       3.    INITIAL GRIDPOINT COMPUTATIONS
!    -------------------------------------------------------------------


!*       3.x IDENTITY MODEL
IF (CDCONF(3:3) /= '0')THEN
  IF (LIDMODEL .AND. CDCONF(6:6)=='V') THEN
    YDGMV5%GMV_DEPART(:,:,1:YDGMV%YT0%NDIM,:) = YDGMV%GMV(:,:,1:YDGMV%YT0%NDIM,:)
    YDGMV5%GMVS_DEPART(:,1:YDGMV%YT0%NDIMS,:) = YDGMV%GMVS(:,1:YDGMV%YT0%NDIMS,:)
    YDGFL5%GFL_DEPART(:,:,1:YGFL%NDIM0,:) = YDGFL%GFL(:,:,1:YGFL%NDIM0,:)
  ENDIF
ENDIF

CALL GSTATS(34,0)

!     ------------------------------------------------------------------
!*       4.0   INITIALIZE TANGENT-LINEAR SURFACE FIELDS 
!     ------------------------------------------------------------------

IF (NSTEP == 0 .AND. YDEPHY%LEPHYS) CALL RESET_SPERT(YDGEOMETRY,YDSURF,YDMODEL%YRML_DYN%YRDYN,NSTEP)

!*       4.4   GRID PONT CALCULATIONS FOR SIGMA_B OR SIGMA_A

IF (CDCONF(6:6) == 'A' .OR. CDCONF(6:6) == 'B') THEN

  ALLOCATE(ZPERTTSG(YDDIM%NDGLG*YDDIM%NDLON))
  ALLOCATE(ZPERTTS(YDGEM%NGPTOT))
  CALL GSTATS(1860,0)
  CALL GAUSSIAN_DISTRIBUTION (ZPERTTSG,YDMODEL%YRML_PHY_STOCH%YR_RANDOM_STREAMS%SCAN2MTL)
  CALL GSTATS(1860,1)
  ZPERTTS(1:YDGEM%NGPTOT)=ZPERTTSG(YDGSGEOM_NB%NUNIQUEGP(1:YDGEM%NGPTOT))
  DEALLOCATE(ZPERTTSG)
  CALL GSTATS(1046,0)

  IF(.NOT.LECMWF)THEN
  ! V. Guidard Temporary fix, should solved a better way
  YGFL%YO3%MP=6
  YGFL%YO3%MP5=6
  ! End VG 
  ENDIF
  CALL GSTATS(1046,1)
  DEALLOCATE(ZPERTTS)
ENDIF

!     ------------------------------------------------------------------
!*       5.    GRIDPOINT MODEL (TANGENT LINEAR).
!     ------------------------------------------------------------------

IF(LEGBRAD)SD_VN(:,YSD_VN%YACCPR5%MP,:)=SD_VN(:,YSD_VN%YACCPR%MP,:)

IF(LECMWF.AND.YDEPHY%LEMWAVE)THEN
  !!update YDPHYSMWAVE status as appropriate, to be used below COBSALL
  YDPHYSMWAVE5%LPHYS_MWAVE_FILLED_IN = .TRUE.
ENDIF

IF (CDCONF(4:4) /= '0') THEN
  CALL GP_MODEL_TL(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,CDCONF,LD_DFISTEP, 0._JPRB,&
   & PTRAJEC_OOPS%SURFACE%Y_SURF,PTRAJEC_OOPS=PTRAJEC_OOPS, &
   & YDPHYSMWAVE5=YDPHYSMWAVE5)
ENDIF


!     ------------------------------------------------------------------
!*       6.    RESET SURFACE-FIELD PERTURBATIONS
!     ------------------------------------------------------------------

IF (NSTEP == YDRIP%NSTOP .AND. YDEPHY%LEPHYS) CALL RESET_SPERT(YDGEOMETRY,YDSURF,YDMODEL%YRML_DYN%YRDYN,NSTEP)


!    -------------------------------------------------------------------
!*    7. NORM LOCALIZATION
!    -------------------------------------------------------------------

IF(CDCONF(5:5) == 'L')THEN

  IF(.NOT.ALLOCATED(YDFIELDS%YRGMV%GMVT1))  ALLOCATE(YDFIELDS%YRGMV%GMVT1(YDDIM%NPROMA,YDDIMV%NFLEVG,YDGMV%YT1%NDIM,YDDIM%NGPBLKS))
  IF(.NOT.ALLOCATED(YDFIELDS%YRGMV%GMVT1S)) ALLOCATE(YDFIELDS%YRGMV%GMVT1S(YDDIM%NPROMA,YDGMV%YT1%NDIMS,YDDIM%NGPBLKS))
  IF(.NOT.ALLOCATED(YDFIELDS%YRGFL%GFLT1))  ALLOCATE(YDFIELDS%YRGFL%GFLT1(YDDIM%NPROMA,YDDIMV%NFLEVG,YGFL%NDIM1,YDDIM%NGPBLKS))

  ICSTA=1

  CALL GSTATS(1896,0)
  DO JKGLO=1,YDGEM%NGPTOT,YDDIM%NPROMA
    ICEND=MIN(YDDIM%NPROMA,YDGEM%NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/YDDIM%NPROMA+1
    IOFF=JKGLO

    IF(LFORCEWR) THEN
      DO JLEV=1,YDDIMV%NFLEVG
        DO JL=ICSTA,ICEND
          GPFORCEU(JL,JLEV,IBL,0)=YDGMV%GMV(JL,JLEV,YDGMV%YT0%MU,IBL)
          GPFORCEV(JL,JLEV,IBL,0)=YDGMV%GMV(JL,JLEV,YDGMV%YT0%MV,IBL)
          GPFORCET(JL,JLEV,IBL,0)=YDGMV%GMV(JL,JLEV,YDGMV%YT0%MT,IBL)
          GPFORCEQ(JL,JLEV,IBL,0)=YDGFL%GFL(JL,JLEV,YGFL%YQ%MP,IBL)
          GPFORCESP(JL,IBL,0)    =YDGMV%GMVS(JL,YDGMV%YT0%MSP,IBL)
        ENDDO
      ENDDO
    ENDIF

    CALL LCNORGGTL(YDGEOMETRY,ICSTA,ICEND,&
     & YDGSGEOM(IBL)%GELAM,YDGSGEOM(IBL)%GELAT,YDGSGEOM(IBL)%GM,&
     & YDGMV%GMV(:,:,YDGMV%YT0%MU,IBL) ,YDGMV%GMV(:,:,YDGMV%YT0%MV,IBL) ,YDGMV%GMV(:,:,YDGMV%YT0%MT,IBL),&
     & YDGFL%GFL(:,:,YGFL%YQ%MP,IBL) ,YDGMV%GMVS(:,YDGMV%YT0%MSP,IBL),&
     & YDGMV%GMVT1(:,:,YDGMV%YT1%MU,IBL) ,YDGMV%GMVT1(:,:,YDGMV%YT1%MV,IBL) ,YDGMV%GMVT1(:,:,YDGMV%YT1%MT,IBL),&
     & YDGFL%GFLT1(:,:,YGFL%YQ%MP1,IBL) ,YDGMV%GMVT1S(:,YDGMV%YT1%MSP,IBL))

  ENDDO
  CALL GSTATS(1896,1)
ENDIF
CALL GSTATS(34,1)

!     ------------------------------------------------------------------
!*       8.    "TIME-STEPPING" FOR GRID-POINT GFL FIELDS
!     ------------------------------------------------------------------

IF (LD_TST_GPGFL) THEN
  DO JGFL=1,YGFL%NUMFLDS
    IF(YGFL%YCOMP(JGFL)%LGP .AND. YGFL%YCOMP(JGFL)%LT1) THEN
      IF (YGFL%YCOMP(JGFL)%NCOUPLING==0) THEN
!$OMP WORKSHARE
        YDGFL%GFL(:,:,YGFL%YCOMP(JGFL)%MP,:) = YDGFL%GFLT1(:,:,YGFL%YCOMP(JGFL)%MP1,:)
!$OMP END WORKSHARE
      ENDIF
    ENDIF
  ENDDO
ENDIF

!    -------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SCAN2MTL_OOPS',1,ZHOOK_HANDLE)
END SUBROUTINE SCAN2MTL_OOPS
