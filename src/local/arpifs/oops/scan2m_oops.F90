#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE SCAN2M_OOPS(YDGEOMETRY,YDFIELDS,YDGMV5,YDMODEL,CDCONF,LD_TST_GPGFL,LD_DFISTEP,YDSPEC,PTRAJEC_OOPS)

!****-------------------------------------------------------------------
!**** *SCAN2M* - Grid-point space computations
!****-------------------------------------------------------------------
!     Purpose.   Computations in grid-point space
!     --------   

!**   Interface.
!     ----------
!        *CALL* *SCAN2M (..)
!        Explicit arguments :  
!        --------------------  

!     INPUT:
!     ------
!        CDCONF       : configuration of work (see doc.)
!        LD_TST_GPGFL : if T, do timestepping on grid-point GFL.

!     INPUT/OUTPUT:
!     -------------


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
!        Mats Hamrud   *ECMWF*

! Modifications
! -------------
!   J. Vivoda  03-2002     PC schemes for NH dynamics (LPC_XXXX keys)
!   P. Smolikova 02-09-30 : variable d4 in NH (GPPAUX etc.)
!   G. Desroziers 02-09-16 : also treat Aladin mean wind trajectory
!   R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!   G.Mozdzynski: 02-10-01 support for radiation on-demand comms
!   C. Fischer   03-07-16: gppaux becomes part of gmv
!   O.Spaniel    03-04-15: cleaning-a same named entity from modules
!   R. El Khatib : 03-04-17 Fullpos improvemnts
!   C. Fischer : 03-06-03 - use igpcomp/ngptot_cap in call to cobs for LAM/Jo
!   R. El Khatib : 03-08-08 fullpos & gfl
!   G. Radnoti: 03-10-01 : update calls to eslextpol1
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!   05-Jan-2004 Y. Seity   time stepping only for not coupled GFL(bug LGPQ=T)
!   M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!   27-Feb-2004 M. Tudor   time stepping only for during predictor     
!   10-Jun-2004 J. Masek   NH cleaning (LPC_NOTR, LFULLIMP)
!   01-Jul-2004 K. Yessad  Make clearer the tests about PC schemes.
!   R. El Khatib : 04-10-21 Fullpos fix for extra gfl
!   Y. Seity : 04-11-16: new arguments to CTVTOT for R and Cp calc.
!   14-Mar-2005 R. Engelen   Fullpos of greenhouse gases
!   14-Mar-2005 A. Untch   Fullpos of aerosols
!   Y.Tremolet    01-Aug-2005 Number sections consistently with TL and AD
!   Y.Bouteloup : 05-01-19 : Add QL QI in the call to capotx
!   Y.Bouteloup : 05-01-31 : Add YR in the call to VPOS
!   G. Hello: 05-01-26 Fullpos fix for arome gfl
!   M.Hamrud  15-Jan-2005  Revised GPRCP
!   G. Desroziers and K. Yessad (sept 2005):
!    - split option LTRAJHR into LTRAJHR_ALTI and LTRAJHR_SURF.
!    - adapt option LTRAJHR to METEO-FRANCE configurations.
!   R. El Ouaraini & R. El Khatib : 05-07-28 Allocate CUF fields and CALL fpinvtrcuf.
!   R. El Khatib  20-May-2005 NFPWIDE moved to YOMWFPDS
!   Y.Bouteloup : 05-12-25 : Add QR QS in the call to capotx
!   O.Spaniel :   22-02-2006 : Phasing CY31
!   M. Bellus :   29-Sep-2006 : Add prognostic convection variables (ALARO-0)
!    - YDAL, YDOM, YUAL, YUOM, YUEN, YUNEBH in the call to VPOS
!   E. Holm      05-10-20: Move jbvcoord_interpolate to cvargptl
!   08-08-2006  S. Serrar diagnostic tracers (TRAC) added to GFL fields
!   M.Hamrud      01-Jul-2006  Revised surface fields
!   K. Yessad (oct 2006): bug correction for LTRAJHR in part 2.
!   J.Haseler 27-Feb-2007  Generalise read of gridpoint trajectory
!   S. Serrar 22-03-2007: changes for the post-processing of ERA40 GFL fields
!   S. Serrar 17-07-2007 post-processing methane related GFL fields
!   S. Serrar 02-05-2008 test on LERA40 removed
!   Y. Seity  11-01-2008: bf for non-coupled GFL not initialized in E zone
!   R. El Khatib  17-Oct-2008 bf for ALADIN fullpos with equal distribution
!                             of C+I+E (RDISTR_E=1.)
!   K. Yessad Dec 2008: merge SLCOMM+SLCOMM1 -> SLCOMM.
!   K. Yessad Dec 2008: merge the different (E)SLEXTPOL.. -> (E)SLEXTPOL.
!   K. Yessad (Aug 2009): rewrite vertical interpolator in post-processing.
!   K. Yessad (Jan 2010): revised code (first draft)
!   G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!   G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM and TOROG
!   J.Hague (Oct 2011): CALL COBSALL moved here from SCAN2MTL
!   J.Hague (Oct 2011): VSURF (Surface derived type), GFL, GMV, GMVS passed down
!   J.Hague (Oct 2011):-L_OK set to true after first iteration (moved from GP_MODEL_TL) 
!   A. Geer  (Jan 2012): COBSALL before and after the GP model for diagnostic physics
!   A. Geer  (Mar 2012): GOM rewrite
!   G. Mozdzynski (May 2012): further cleaning
!   K. Yessad (Dec 2011)   : various contributions.
!   C. Soci  27-02-2012    : add Model orog and lsm for Mesan corr function
!   U. Andrae (Apr 2012)   : add graupels (ZGT) in the call to capotx
!   R. El Khatib 18-Jul-2012 Fullpos move away
!   K. Yessad (Nov 2012): simplify testings.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   F. Vana  28-Nov-2013   : Redesigned trajectory handling 
!   R. El Khatib 02-Oct-2014 gp_model allocating large arrays on heap or stack
!   K. Yessad (July 2014): Move some variables.
!   R. El Khatib 12-Aug-2016 cleaning
!   F. Suzat     08-Apr-2018 broke this routine to compile because oopsification of capotx
!   H Petithomme (Dec 2020): introduce gp_model_drv and gpopertraj, allocate t1 arrays here now
!-----------------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE FIELDS_MOD         , ONLY : FIELDS
USE SURFACE_FIELDS_MIX , ONLY : GPOPER,TSURF
USE YOMGMV             , ONLY : TGMV
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE YOMLUN             , ONLY : NULOUT 
USE YOMCT0             , ONLY : LELAM,LCANARI
USE YOMCT3             , ONLY : NSTEP
USE QAGPSF             , ONLY : RSGPHI, RSGHUM
USE QALORI             , ONLY : QCAGUE
USE YOMLCZ             , ONLY : LFORCEWR, GPFORCEU, GPFORCEV, GPFORCET, GPFORCEQ, GPFORCESP  
USE TESTVAR_MIX        , ONLY : LTESTINC
USE QASSET             , ONLY : ESIG
USE YOMCST             , ONLY : RG
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)
USE YOMTRAJ_OOPS       , ONLY : TRAJ_TYPE_OOPS
USE YOMTRAJ            , ONLY : NTRAJ_CST, NGP5,MKINDTRAJ
USE YOMDYN,ONLY: TDYN
USE YOMDIM,ONLY: TDIM

!     ------------------------------------------------------------------

IMPLICIT NONE
TYPE(GEOMETRY)       , INTENT(IN)    :: YDGEOMETRY
TYPE(FIELDS)         , INTENT(INOUT) :: YDFIELDS
TYPE(TGMV)           , INTENT(INOUT) :: YDGMV5
TYPE(MODEL)           ,INTENT(INOUT) :: YDMODEL
CHARACTER(LEN=9)     , INTENT(IN)    :: CDCONF
LOGICAL              , INTENT(IN)    :: LD_TST_GPGFL
LOGICAL              , INTENT(IN)    :: LD_DFISTEP
TYPE(SPECTRAL_FIELD) , INTENT(INOUT) :: YDSPEC
TYPE(TRAJ_TYPE_OOPS) , INTENT(INOUT) :: PTRAJEC_OOPS
!     ------------------------------------------------------------------

!    GRIBEX cannot handle more than one type of REAL and thus cannot accept MKINDTRAJ arrays

REAL(KIND=JPRB),ALLOCATABLE :: ZUT(:,:),ZVT(:,:),ZTT(:,:),ZQT(:,:),ZLT(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZIT(:,:),ZRT(:,:),ZST(:,:),ZGT(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZSP_SB(:,:,:),ZSP_SG(:,:,:),ZSP_RR(:,:),ZSP_CI(:,:),ZSP_X2(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZSD_VF(:,:),ZSD_VV(:,:),ZSD_VX(:,:)
REAL(KIND=JPRB),DIMENSION(:),ALLOCATABLE :: ZPS,ZRCORI,ZGEMU,ZMORO,ZMLSM
REAL(KIND=JPRB),DIMENSION(:),ALLOCATABLE :: ZGM,ZGELAT,ZGELAM,ZGNORDL,ZGNORDM
REAL(KIND=JPRB),DIMENSION(:,:),ALLOCATABLE :: ZCAGUE,ZESIG,ZSGPHI,ZSGHUM

INTEGER(KIND=JPIM) :: IBL, ICEND, IEND, INBPT, INDEX,&
 & IOFF, IST, ITASK, JKGLO, JROF, JLEV, JL,&
 & IGPCOMP, JGFL

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "user_clock.h"

#include "abor1.intfb.h"
#include "gpnorm3.intfb.h"
#include "gp_model_drv.intfb.h"
#include "gprcp_pgfl.intfb.h"
#include "sujbvcoord.intfb.h"
#include "gp_derivatives.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SCAN2M_OOPS',0,ZHOOK_HANDLE)
ASSOCIATE(YDGFL=>YDFIELDS%YRGFL,YDGMV=>YDFIELDS%YRGMV, YDSURF=>YDFIELDS%YRSURF,   LEGBRAD=>YDMODEL%YRML_PHY_EC%YREPHY%LEGBRAD,   &
& YDDIM=>YDGEOMETRY%YRDIM,  YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, YDRIP=>YDMODEL%YRML_GCONF%YRRIP,                 &
& YGFL=>YDMODEL%YRML_GCONF%YGFL,YDDYN=>YDMODEL%YRML_DYN%YRDYN)

ASSOCIATE(NDIM=>YGFL%NDIM, NUMFLDS=>YGFL%NUMFLDS,   YQ=>YGFL%YQ,   YRSPEC=>YGFL%YRSPEC,  NGPBLKS=>YDDIM%NGPBLKS, &
& NPROMA=>YDDIM%NPROMA, NFLEVG=>YDDIMV%NFLEVG,   NGPTOT=>YDGEM%NGPTOT, NGPTOT_CAP=>YDGEM%NGPTOT_CAP,             &
& GFL=>YDGFL%GFL,   GMV=>YDGMV%GMV, GMVS=>YDGMV%GMVS, YT0=>YDGMV%YT0,   NGP5_OOPS=>YDSURF%NGP5_OOPS,             &
& SD_VN=>YDSURF%SD_VN,   YSD_VN=>YDSURF%YSD_VN)
!     ------------------------------------------------------------------
!*       1.    INITIAL SETUP, MAIN LOOP CONTROL.
!     ------------------------------------------------------------------

!*       1.1   INITIAL SETUP

! Total number of grid points - ndglg*ndlon in global / ndguxg*ndlon in LAM
! used for model computations
IF (LELAM) THEN
  IGPCOMP=MIN(NGPTOT,NGPTOT_CAP)
ELSE
  IGPCOMP=NGPTOT
ENDIF

!     ------------------------------------------------------------------
!*       2.    INITIAL GRIDPOINT COMPUTATIONS
!     ------------------------------------------------------------------

IF (.NOT. LFORCEWR) CALL GSTATS(8,0)

IF (CDCONF(3:3) == 'C' .OR. CDCONF(3:3) == 'D') THEN
  CALL ABOR1('SCAN2M: CDCONF(3:3)= C or D not supported')
ENDIF



!*       2.1  FORCE

IF (LFORCEWR) THEN
  DO JKGLO=1,NGPTOT,NPROMA
    ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1
    DO JLEV=1,NFLEVG
      DO JL=1,ICEND
        GPFORCEU(JL,JLEV,IBL,0)=GMV(JL,JLEV,YT0%MU,IBL)
        GPFORCEV(JL,JLEV,IBL,0)=GMV(JL,JLEV,YT0%MV,IBL)
        GPFORCET(JL,JLEV,IBL,0)=GMV(JL,JLEV,YT0%MT,IBL)
        GPFORCEQ(JL,JLEV,IBL,0)=GFL(JL,JLEV,YQ%MP,IBL)
      ENDDO
    ENDDO
    DO JL=1,ICEND
      GPFORCESP(JL,IBL,0)    =GMVS(JL,YT0%MSP,IBL)
    ENDDO
  ENDDO
  IF (LHOOK) CALL DR_HOOK('SCAN2M_OOPS',1,ZHOOK_HANDLE)
  RETURN
ENDIF

IF(CDCONF(3:3) /= '0' .AND. YRSPEC%LGP) THEN
  DO JKGLO=1,NGPTOT,NPROMA
    IST =1
    IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1
    CALL GPRCP_PGFL(NPROMA,IST,IEND,NFLEVG,PGFL=GFL(1,1,1,IBL),PR=GFL(1,1,YRSPEC%MP,IBL))
  ENDDO
  CALL GPNORM3(YDGEOMETRY,GFL,NDIM,YRSPEC%MP,YRSPEC%CNAME)
ENDIF
IF (CDCONF(3:3) /= '0') THEN
  CALL GP_DERIVATIVES(YDGEOMETRY,YDFIELDS,YGFL)
ENDIF


!     ------------------------------------------------------------------
!*       3.    TRAJECTORY HANDLING
!     ------------------------------------------------------------------

! Store surface and grid-point trajectory
! ---------------------------------------

IF(CDCONF(1:1)=='T')THEN
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,IBL,IEND)
    DO JKGLO=1,NGPTOT,NPROMA
      IBL=(JKGLO-1)/NPROMA+1
      IEND=MIN(NPROMA,NGPTOT-JKGLO+1)

      CALL GPOPERTRAJ(YDDIM,YDDYN,"TRAJSTORECST",YDSURF,IBL,IEND,&
       & PTRAJEC_OOPS%CST(IBL)%MIX_CST,NTRAJ_CST)
      IF (PTRAJEC_OOPS%CST(IBL)%LASTCHUNK)&
        WRITE(NULOUT,*)'GREPTRAJ TRAJEC%CST: constants saved in SCAN2M '
    ENDDO
!$OMP END PARALLEL DO
ENDIF


IF(CDCONF(1:1)=='T')THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IBL,IEND)
      DO JKGLO=1,NGPTOT,NPROMA
        IBL=(JKGLO-1)/NPROMA+1
        IEND=MIN(NPROMA,NGPTOT-JKGLO+1)

        CALL GPOPERTRAJ(YDDIM,YDDYN,"TRAJSTORE",YDSURF,IBL,IEND,&
         & PTRAJEC_OOPS%SRFC(IBL)%MIX_SRFC,NGP5_OOPS)
      ENDDO
!$OMP END PARALLEL DO
ENDIF

!     ------------------------------------------------------------------
!*       5.    GRIDPOINT MODEL.
!     ------------------------------------------------------------------
! GBRAD: Save Surface field for GOM
IF(LEGBRAD)SD_VN(:,YSD_VN%YACCPR5%MP,:)=SD_VN(:,YSD_VN%YACCPR%MP,:)

IF (CDCONF(4:4) /= '0') THEN
  CALL GSTATS(1228,0)
  IF(.NOT.ALLOCATED(YDGMV%GMVT1))  THEN
    ALLOCATE(YDGMV%GMVT1(NPROMA,NFLEVG,YDGMV%YT1%NDIM,NGPBLKS))
    ALLOCATE(YDGMV%GMVT1S(NPROMA,YDGMV%YT1%NDIMS,NGPBLKS))
    ALLOCATE(YDGFL%GFLT1(NPROMA,NFLEVG,YGFL%NDIM1,NGPBLKS))
  ENDIF

  !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(IBL)
  DO IBL=1,NGPBLKS
    YDGMV%GMVT1(:,:,:,IBL)=0.0_JPRB
    YDGMV%GMVT1S(:,:,IBL)=0.0_JPRB
    YDGFL%GFLT1(:,:,:,IBL)=0.0_JPRB
  ENDDO
  !$OMP END PARALLEL DO

  CALL GSTATS(1228,1)

  IF(.NOT.LCANARI)&
   & CALL GP_MODEL_DRV(YDGEOMETRY,YDFIELDS,YDMODEL,CDCONF,LD_DFISTEP,PTRAJEC_OOPS=PTRAJEC_OOPS)
ENDIF

!     ------------------------------------------------------------------
!*       7.    GRID POINT CALCULATIONS FOR ANALYSIS.
!     ------------------------------------------------------------------

!*       7.1  GRID PONT CALCULATIONS FOR JG

IF (CDCONF(6:6) == 'H') THEN
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,ICEND,IBL,IOFF)
  DO JKGLO=1,NGPTOT,NPROMA
    ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1

    IOFF=JKGLO
    CALL SUJBVCOORD(YDGEOMETRY,ICEND,JKGLO,&
     & GMV(1,1,YT0%MT,IBL),GFL(1,1,YQ%MP,IBL),GMVS(1,YT0%MSP,IBL))  
  ENDDO
!$OMP END PARALLEL DO
ENDIF

!*       7.2  GRID POINT CALCULATIONS FOR CANARI

IF (CDCONF(6:6) == '1') THEN
! F. SUZAT should never append...
      CALL ABOR1('F. SUZAT commented to compile, should never happend')

ENDIF

CALL GSTATS(8,1)

!     ------------------------------------------------------------------
!*       8.    "TIME-STEPPING" FOR GRID-POINT GFL FIELDS
!     ------------------------------------------------------------------

IF (LD_TST_GPGFL) THEN
  CALL GSTATS(1032,0)
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,ICEND,IBL,JGFL,JLEV,JL)
  DO JKGLO=1,IGPCOMP,NPROMA
    IBL=(JKGLO-1)/NPROMA+1
    ICEND=MIN(NPROMA,IGPCOMP-JKGLO+1)
    DO JGFL=1,NUMFLDS
      IF(YGFL%YCOMP(JGFL)%LGP.AND.YGFL%YCOMP(JGFL)%LT1) THEN
        IF (YGFL%YCOMP(JGFL)%NCOUPLING==0) THEN
          YDFIELDS%YRGFL%GFL(1:ICEND,1:NFLEVG,YGFL%YCOMP(JGFL)%MP,IBL) =&
           & YDFIELDS%YRGFL%GFLT1(1:ICEND,1:NFLEVG,YGFL%YCOMP(JGFL)%MP1,IBL)
        ENDIF
      ENDIF
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  CALL GSTATS(1032,1)
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SCAN2M_OOPS',1,ZHOOK_HANDLE)
CONTAINS
  ! version of GPOPERTRAJ simplified compared to SCAN2M
  ! PFIELD needs a copy from ZTEMP because of its kind MKINDTRAJ /= JPRB
  SUBROUTINE GPOPERTRAJ(YDDIM,YDDYN,CDNAME,YDSURF,IBL,KF,PFIELD,NF)
    TYPE(TDIM),INTENT(IN) :: YDDIM
    TYPE(TDYN),INTENT(IN) :: YDDYN
    TYPE(TSURF),INTENT(INOUT) :: YDSURF
    INTEGER(KIND=JPIM),INTENT(IN) :: IBL,KF,NF
    CHARACTER(LEN=*),INTENT(IN) :: CDNAME
    REAL(KIND=MKINDTRAJ),INTENT(OUT) :: PFIELD(:,:)

    REAL(KIND=JPRB) :: ZTEMP(YDDIM%NPROMA,NF)

    CALL GPOPER(YDDIM,YDDYN,CDNAME,YDSURF,KBL=IBL,PFIELD=ZTEMP)

    PFIELD(1:KF,:) = ZTEMP(1:KF,:)
  END SUBROUTINE
END SUBROUTINE SCAN2M_OOPS
