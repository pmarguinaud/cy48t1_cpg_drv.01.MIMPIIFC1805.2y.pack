SUBROUTINE SCAN2MAD(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDGOM5,YDGOM,CDCONF,LD_TST_GPGFL,PTRAJEC,PII0,YDACV)

!**** *SCAN2MAD* - Grid-point space computations (adjoint)

!     Purpose. Computations in grid-point space
!     --------

!     Interface.
!     ----------
!        *CALL* *SCAN2MAD(....)

!        Explicit arguments :  CDCONF - configuration of work (see doc.)
!        --------------------

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   see includes below.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud   *ECMWF*
!      Original : 91-04-01

!     Modifications.
!     --------------
!      Modified 02-09-30 P. Smolikova : Interface to SC2CGAP for d4 in NH
!      Modified 03-06-02 C. Soci : nullify GPP when reallocated
!      C. Fischer : 03-06-03 use igpcomp/ngptot_cap in call to cobs for LAM/Jo
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib 04-05-24 Allocate MASK_SL1 anytime because it is argument of subroutine
!      C Soci & C Fischer 04-06-18 GMVs set to zero up to NGPTOT for aladin/ifs compliency
!      Y.Tremolet    01-Aug-2005 Number sections consistently with SCAN2M
!      D Salmond 21-09-05 Fix for LSPRT=.true.
!      E. Holm   06-07-11 Move NSTEP=0 to CVARGPAD and remove unused option CDCONF(6)=W.
!      K. Yessad Dec 2008: merge the different (E)SLEXTPOL.. -> (E)SLEXTPOL.
!      K. Yessad (Jan 2010): revised code (first draft)
!      H. Hersbach Dec 2009: reset surface perts (moved from ec_phys_ad)
!      O.Riviere Oct 2010 : move reset_spert under key lephys
!      G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!      J.Hague Mar 2011: Replace part with COBSALLAD for Observation Interpolation
!      A. Geer Jan 2012: COBSALL before and after the GP model for diagnostic physics
!      F. vana 19-03-2012 : all time-steps treated consistently in terms Tv -> T conversion
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      F. Vana 28-Nov-2013 : Redesigned trajectory handling
!      K. Yessad (July 2014): Move some variables.
!      A. Geer Aug 2013/Mar 2016: Single call to cobsall now; use supergom
!      A. Geer 05 Apr 2016 Use LTESTADJM to prevent use of GOMs in scan2mtl/ad
!      S. Massart 19-Feb-2019 : Solar constant optimisation
!    -------------------------------------------------------------------

USE TYPE_MODEL    , ONLY : MODEL
USE GEOMETRY_MOD  , ONLY : GEOMETRY
USE FIELDS_MOD    , ONLY : FIELDS
USE MTRAJ_MOD     , ONLY : MTRAJ
USE PARKIND1      , ONLY : JPIM, JPRB
USE YOMHOOK       , ONLY : LHOOK, DR_HOOK
USE YOMCT0        , ONLY : LSPRT, LECMWF
USE YOMCT3        , ONLY : NSTEP
USE YOMDYNA       , ONLY : LGRADSP
USE YOMECTAB      , ONLY : YECVAR
USE YOMVRTL       , ONLY : LIDMODEL
USE YOMVRTLX      , ONLY : L801TL
USE TESTVAR_MIX   , ONLY : LTESTADJM
USE TRAJECTORY_MOD, ONLY : LTRAJGP, GET_TRAJ_GRID
USE YOMTRAJ       , ONLY : MSTEPTRAJR, TRAJ_TYPE, LPRTTRAJ
USE YOMLUN        , ONLY : NULOUT
USE SUPERGOM_CLASS, ONLY : CLASS_SUPERGOM
USE TYPE_ACV      , ONLY : ACV_CONTAINER

!    -------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY  !! INOUT needed for GP_MODEL_AD
TYPE(FIELDS)        ,INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ)         ,INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL)         ,INTENT(INOUT) :: YDMODEL
TYPE(CLASS_SUPERGOM),INTENT(IN)    :: YDGOM5
TYPE(CLASS_SUPERGOM),INTENT(INOUT) :: YDGOM
CHARACTER(LEN=9)    ,INTENT(IN)    :: CDCONF
LOGICAL             ,INTENT(IN)    :: LD_TST_GPGFL
TYPE(TRAJ_TYPE)     ,INTENT(IN)    :: PTRAJEC
REAL(KIND=JPRB)    , INTENT(INOUT) :: PII0
TYPE(ACV_CONTAINER) ,INTENT(INOUT), OPTIONAL :: YDACV
!    -------------------------------------------------------------------
INTEGER(KIND=JPIM) :: IBL, ICEND, IEND,ICSTA,IOFF,JKGLO,IST, ISTEPR
INTEGER(KIND=JPIM) :: JGFL, JROF, JLEV

REAL(KIND=JPRB), ALLOCATABLE :: ZSD_VD(:,:,:) 

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!    -------------------------------------------------------------------

#include "abor1.intfb.h"
#include "ctvtot5.intfb.h"
#include "ctvtotad.intfb.h"
#include "gp_model_ad.intfb.h"
#include "reset_spert.intfb.h"
#include "lcnorggad.intfb.h"

!    -------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SCAN2MAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDGFL5=>YDMTRAJ%YRGFL5, YDGMV5=>YDMTRAJ%YRGMV5, YDGFL=>YDFIELDS%YRGFL, &
 & YDGMV=>YDFIELDS%YRGMV, YDSURF=>YDFIELDS%YRSURF, &
 & YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM, YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YGFL=>YDMODEL%YRML_GCONF%YGFL, &
 & YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY)

ASSOCIATE(NDIM0=>YGFL%NDIM0, NDIM5=>YGFL%NDIM5, NUMFLDS=>YGFL%NUMFLDS, &
 & YCOMP=>YGFL%YCOMP, YQ=>YGFL%YQ, &
 & NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LEPHYS=>YDEPHY%LEPHYS, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & GFL=>YDGFL%GFL, GFLT1=>YDGFL%GFLT1, &
 & GFL_DEPART=>YDGFL5%GFL_DEPART, &
 & GMV=>YDGMV%GMV, GMVS=>YDGMV%GMVS, GMVT1=>YDGMV%GMVT1, GMVT1S=>YDGMV%GMVT1S, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, &
 & YT0=>YDGMV%YT0, YT1=>YDGMV%YT1, &
 & GMVS_DEPART=>YDGMV5%GMVS_DEPART, GMV_DEPART=>YDGMV5%GMV_DEPART, &
 & YT5=>YDGMV5%YT5, &
 & NSTOP=>YDRIP%NSTOP, &
 & LEGBRAD=>YDEPHY%LEGBRAD, &
 & SD_VD=>YDSURF%SD_VD, SD_VN=>YDSURF%SD_VN, &
 & YSD_VN=>YDSURF%YSD_VN, YSD_VDD=>YDSURF%YSD_VDD, YSD_VND=>YDSURF%YSD_VND)

!    -------------------------------------------------------------------

!    -------------------------------------------------------------------
!*       1.    INITIAL SETUP, MAIN LOOP CONTROL.
!    -------------------------------------------------------------------
ISTEPR=MSTEPTRAJR(NSTEP)

!    -------------------------------------------------------------------
!*       2.   TRAJECTORY HANDLING
!    -------------------------------------------------------------------

IF((CDCONF(4:4) /= '0'.OR. CDCONF(6:6) == 'G' .OR. CDCONF(6:6) == 'X').OR.LSPRT) THEN
  IF( LTRAJGP )THEN
    ! retrieve grid-point trajectory
    IF (.NOT.ALLOCATED(YDGMV5%GMV5))  ALLOCATE(YDGMV5%GMV5(NPROMA,NFLEVG,YT5%NDIM,NGPBLKS))
    IF (.NOT.ALLOCATED(YDGMV5%GMV5S)) ALLOCATE(YDGMV5%GMV5S(NPROMA,YT5%NDIMS,NGPBLKS))
    IF (.NOT.ALLOCATED(YDGFL5%GFL5))  ALLOCATE(YDGFL5%GFL5(NPROMA,NFLEVG,NDIM5,NGPBLKS))
    CALL GSTATS(16,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,ICEND,IBL)
    DO JKGLO=1,NGPTOT,NPROMA
      ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      IBL=(JKGLO-1)/NPROMA+1
      IF(.NOT.ASSOCIATED(PTRAJEC%MAIN(IBL)%GMV))&
         & CALL ABOR1('READING TRAJEC%MAIN: TRAJEC%MAIN%GMV NOT ALLOCATED')  
      YDGMV5%GMV5(1:ICEND,:,1:YT5%NDIM,IBL)=PTRAJEC%MAIN(IBL)%GMV(1:ICEND,:,1:YT5%NDIM)
      YDGMV5%GMV5S(1:ICEND,1:YT5%NDIMS,IBL)=PTRAJEC%MAIN(IBL)%GMVS(1:ICEND,1:YT5%NDIMS)
      YDGFL5%GFL5(1:ICEND,:,1:NDIM5,IBL)=PTRAJEC%MAIN(IBL)%GFL(1:ICEND,:,1:NDIM5)
    ENDDO
!$OMP END PARALLEL DO
    CALL GSTATS(16,1)
    IF (LPRTTRAJ)  WRITE(NULOUT,*)'GREPTRAJ GET_TRAJ_GRID GRID istep=',ISTEPR
  ELSE
    IF(.NOT.ALLOCATED(YDGMV5%GMV5)) CALL ABOR1('SCAN2MAD:GMV5 NOT ALLOCATED')
  ENDIF
ENDIF

! If LSPRT=.T. trajectory temperature contains Tv, so it is converted to T here

IF (LSPRT) THEN
  CALL GSTATS(1044,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IBL,IST,IEND)
  DO JKGLO=1,NGPTOT,NPROMA
    IBL=(JKGLO-1)/NPROMA+1
    IST =1
    IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    CALL CTVTOT5(YDGEOMETRY,1,IEND,YDGMV5%GMV5(:,:,YT5%MT,IBL),YDGFL5%GFL5(:,:,YQ%MP5,IBL))
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1044,1)
ENDIF

!     ------------------------------------------------------------------
!*       8.    "TIME-STEPPING" FOR GRID-POINT GFL FIELDS
!     ------------------------------------------------------------------

IF (LD_TST_GPGFL) THEN
  CALL GSTATS(1044,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IBL,JGFL)
  DO JKGLO=1,NGPTOT,NPROMA
    IBL=(JKGLO-1)/NPROMA+1
    DO JGFL=1,NUMFLDS
      IF(YCOMP(JGFL)%LGP.AND.YCOMP(JGFL)%LT1) THEN
        IF (YCOMP(JGFL)%NCOUPLING==0) THEN
          GFLT1(:,:,YCOMP(JGFL)%MP1,IBL) = GFL(:,:,YCOMP(JGFL)%MP,IBL)
          GFL(:,:,YCOMP(JGFL)%MP,IBL)= 0._JPRB
        ENDIF
      ENDIF
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1044,1)
ENDIF

! Following loop moved from gp_model_ad
CALL GSTATS(1044,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,IBL,JGFL)
DO JKGLO=1,NGPTOT,NPROMA
  IBL=(JKGLO-1)/NPROMA+1
  IF( .NOT. LGRADSP ) THEN
    GMV(:,:,1:YT0%NDIM,IBL)=0.0_JPRB
  ELSE
    GMV(:,:,1:YT0%NDIM-2,IBL)=0.0_JPRB
    IF (NSTEP == NSTOP-1) GMV(:,:,YT0%NDIM-1:YT0%NDIM,IBL)=0.0_JPRB
  ENDIF
  GMVS(:,1:YT0%NDIMS,IBL)=0.0_JPRB
  DO JGFL=1,NUMFLDS
    IF(YCOMP(JGFL)%LT1) THEN
      IF (LTESTADJM .AND. YCOMP(JGFL)%LGP .AND. .NOT.LD_TST_GPGFL) THEN
        ! GR: in test of adjoint we start from a non-zero state so
        ! we must avoid zeroing if there was no swapping to GFLT1 before...      
        ! to be cleaned...
      ELSE
         GFL(:,:,YCOMP(JGFL)%MP,IBL) = 0.0_JPRB
      ENDIF
      IF(YCOMP(JGFL)%LCDERS) THEN
        GFL(:,:,YCOMP(JGFL)%MPL,IBL) = 0.0_JPRB
        GFL(:,:,YCOMP(JGFL)%MPM,IBL) = 0.0_JPRB
      ENDIF
    ENDIF
    IF (YCOMP(JGFL)%LDIAG .AND. YCOMP(JGFL)%NCOUPLING==0) THEN
      ! Initialise adjoint to zero for diagnostic fields
      GFL(:,:,YCOMP(JGFL)%MP,IBL)= 0._JPRB
    ENDIF
  ENDDO
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1044,1)

!    -------------------------------------------------------------------
!*    7. NORM LOCALIZATION
!    -------------------------------------------------------------------

IF(CDCONF(5:5) == 'L')THEN
  ICSTA=1
  CALL GSTATS(1045,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,ICEND,IBL,IOFF)
  DO JKGLO=1,NGPTOT,NPROMA
    ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1
    IOFF=JKGLO
    CALL LCNORGGAD(YDGEOMETRY,ICSTA,ICEND,&
     & YDGSGEOM(IBL)%GELAM,YDGSGEOM(IBL)%GELAT,YDGSGEOM(IBL)%GM,&
     & GMV(:,:,YT0%MU,IBL) ,GMV(:,:,YT0%MV,IBL) ,GMV(:,:,YT0%MT,IBL),&
     & GFL(:,:,YQ%MP,IBL) ,GMVS(:,YT0%MSP,IBL),&
     & GMVT1(:,:,YT1%MU,IBL) ,GMVT1(:,:,YT1%MV,IBL) ,GMVT1(:,:,YT1%MT,IBL),&
     & GFLT1(:,:,YQ%MP1,IBL) ,GMVT1S(:,YT1%MSP,IBL))
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1045,1)
ENDIF

!     ------------------------------------------------------------------
!*       6.    INITIALIZE ADJOINT SURFACE FIELDS
!     ------------------------------------------------------------------

IF (NSTEP == NSTOP .AND. LEPHYS) CALL RESET_SPERT(YDGEOMETRY,YDSURF,YDMODEL%YRML_DYN%YRDYN,NSTEP)

!     ------------------------------------------------------------------
!*       4.    COMPARISON WITH OBSERVATIONS
!     ------------------------------------------------------------------

IF(LECMWF.AND.YDEPHY%LEMWAVE)YDFIELDS%YEC_PHYS_FIELDS%YRPHYSMWAVE%PHYS_MWAVE=0.0_JPRB
IF ((CDCONF(6:6) == 'V' .OR. CDCONF(6:6) == 'F') .AND. .NOT.LTESTADJM) THEN
  ALLOCATE (ZSD_VD(NPROMA,YSD_VDD%NDIM,NGPBLKS))
  ZSD_VD = 0.0_JPRB
  IF (LEGBRAD) SD_VN(:,YSD_VN%YACCPR5%MP,:)=0.0_JPRB
  CALL YDGOM%MODEL_IN_AD(YDGEOMETRY,YDGMV,YDSURF,YDMODEL%YRML_GCONF%YRDIMACV,YECVAR,&
 &                       YGFL,YDFIELDS%YEC_PHYS_FIELDS%YRPHYSMWAVE,GFL,GMV,GMVS,&
 &                       ZSD_VD,SD_VN,YDACV=YDACV)
ENDIF

!     ------------------------------------------------------------------
!*       5.    GRIDPOINT MODEL (ADJOINT).
!     ------------------------------------------------------------------

CALL GSTATS(35,0)

IF (CDCONF(4:4) /= '0') THEN
  CALL GP_MODEL_AD(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,CDCONF,PII0,YDFIELDS%YRSURF,&
    & PTRAJEC=PTRAJEC,YDACV=YDACV)
ENDIF

CALL GSTATS(35,1)

IF ((CDCONF(6:6) == 'V' .OR. CDCONF(6:6) == 'F') .AND. .NOT.LTESTADJM) THEN
  ! Cobsall uses the before-timestep version of all fields except those that
  ! are diagnostic from the model. For some fields the two time-levels are
  ! not kept separate, so they require special treatment. 
IF(LEGBRAD)THEN
  SD_VN(:,YSD_VN%YACCPR%MP,:)=SD_VN(:,YSD_VN%YACCPR%MP,:)+SD_VN(:,YSD_VN%YACCPR5%MP,:)
  SD_VN(:,YSD_VN%YACCPR5%MP,:)=0.0_JPRB
ENDIF
  SD_VD = SD_VD + ZSD_VD
  DEALLOCATE(ZSD_VD)
ENDIF

!    -------------------------------------------------------------------
!*       4.0 RESET ADJOINT SURFACE FIELDS 
!    -------------------------------------------------------------------

IF (NSTEP == 0 .AND. LEPHYS) CALL RESET_SPERT(YDGEOMETRY,YDSURF,YDMODEL%YRML_DYN%YRDYN,NSTEP)

!    -------------------------------------------------------------------
!*       3.    INITIAL GRIDPOINT COMPUTATIONS
!    -------------------------------------------------------------------

!*       3.2 CONVERSION OF TTx ARRAY FROM 'TV' TO T

CALL GSTATS(35,0)

IF (LSPRT) THEN
  CALL GSTATS(1047,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IBL,IST,IEND)
  DO JKGLO=1,NGPTOT,NPROMA
    IBL=(JKGLO-1)/NPROMA+1
    IST =1
    IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    CALL CTVTOTAD(YDGEOMETRY,YDMODEL%YRML_DYN%YRDYN,1,IEND,&
     & GMV(:,:,YT0%MT,IBL),GFL(:,:,YQ%MP,IBL),&
     & YDGMV5%GMV5(:,:,YT5%MT,IBL),YDGFL5%GFL5(:,:,YQ%MP5,IBL))
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1047,1)
ENDIF

CALL GSTATS(35,1)

!*       3.3 SETTING SENSITIVITY TO Q TO ZERO WHEN Q IS VERY SMALL
!             (to avoid spurious sensitivities to Q)

IF (L801TL) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IBL,IST,IEND)
  DO JKGLO=1,NGPTOT,NPROMA
    IBL=(JKGLO-1)/NPROMA+1
    IST =1
    IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    DO JLEV=1,NFLEVG
      DO JROF=IST,IEND
        IF (YDGFL5%GFL5(JROF,JLEV,YQ%MP5,IBL) < 0.0001_JPRB) THEN
          GFL(JROF,JLEV,YQ%MP,IBL) = 0.0_JPRB
        ENDIF
      ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDIF

!*       3.x IDENTITY MODEL
IF(CDCONF(3:3) /= '0')THEN
  IF (LIDMODEL .AND. CDCONF(6:6) == 'V') THEN
!$OMP WORKSHARE
    GMV(:,:,1:YT0%NDIM,:) = GMV_DEPART(:,:,1:YT0%NDIM,:)
    GMVS(:,1:YT0%NDIMS,:) = GMVS_DEPART(:,1:YT0%NDIMS,:)
    GFL(:,:,1:NDIM0,:) = GFL_DEPART(:,:,1:NDIM0,:)
!$OMP END WORKSHARE
  ENDIF
ENDIF

!    -------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SCAN2MAD',1,ZHOOK_HANDLE)
END SUBROUTINE SCAN2MAD
