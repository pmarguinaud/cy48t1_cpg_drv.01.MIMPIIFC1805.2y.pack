SUBROUTINE OBSOP_CONV(ROBHDR,ROBODY,SATHDR,SATBODY,ROBSU,SCATTBODY,YDGP5,VARNOS_TO_PROCESS,LDSCREEN,&
  & LDFINAL,YDPHY,YDSET,PHOFX,YDGP_TL,YDGP_AD)

! Observation operator for most "conventional" data (roughly, anything observing a simple
! geophysical equivalent of the model state for which PPNEW can be used). So this
! does include a few satellite types like SCATT and SATOB.

!**   MODIFICATIONS
!     -------------

! A. Geer          01 Oct 2015  FIRST VERSION: code moved out of old "hop"
! T. Montmerle     04 Nov 2015 : add calls to TL/AD
! C. Payan         03 Dec 2015 : TL/AD updates (calls to ppnew)
! P. Lean          19 Sep 2015 Make use of new odb interface features
! F.Suzat          8  Apr 2018 report 43t2bf to make canari running 
! M. Rennie        20 Jun 2018 Avoid filling in Aeolus L2C analysis info with OOPS
! T.Montmerle         Aug 2018 : simplify calls to PPOBSAC*
! M. Rennie        06 Dec 2018 Ensure AD zeroed for Aeolus
! F. Duruisseau    27-dec-2018  Allow to use conv with SATEM (BAYRAD)
! B. Ingleby       16 Jan 2019 Add gom_plus route for surface T, q and wind
! S. Quesada-Ruiz / R. Dragani  14 Mar 2019 Add AK to O3
! B. Ingleby        9 Oct 2019 Add lapse rate adjustment for T2m
  
USE PARKIND1           , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE YOMCT0             , ONLY : LECMWF, LCANARI, NCONF, L_OOPS
USE YOMCST             , ONLY : RG, TCST, TCST_INIT
USE ALGORITHM_STATE_MOD, ONLY : GET_NUPTRA, L_OBS_IN_FC
USE YOMOBS             , ONLY : LCAPACH, LCACHMT, LSCATT_NEUTRAL, &
 &                              L_SCREEN_LEVEL_OBSOP, R4DLAPSET2
USE VARNO_MODULE       , ONLY : NVNUMAX, VARNO
USE YOMCOCTP           , ONLY : NAIREP, NSATOB, NTEMP, NPILOT, NSCATT,&
 &                              NSYNOP, NDRIBU, NPAOB, NISAC, NMDEHS
USE YOMDB
USE YOMANCS            , ONLY : RMDI
USE YOMSATS            , ONLY : SATOBGRP_TABLE
USE ENKF_MIX           , ONLY : LENKF
USE GOM_PLUS           , ONLY : TYPE_GOM_PLUS, IH
USE IFS_DBASE_VIEW_MOD , ONLY : IFS_DBASE_VIEW
USE OBSOP_SETS         , ONLY : TYPE_SET_INFO
USE YOMPHY             , ONLY : TPHY
USE YOMDPHY            , ONLY : YRDPHY
USE PARCMA             , ONLY : JPMX_AK
USE YOMLUN   , ONLY : NULOUT
USE YOMCMBDY , ONLY : NCMDPABP ,NCMDPAOC

!     -------------------------------------------------------------------------

IMPLICIT NONE

TYPE(IFS_DBASE_VIEW) , INTENT(INOUT) :: ROBHDR,ROBODY,SATHDR,SATBODY,ROBSU,SCATTBODY
TYPE(TYPE_GOM_PLUS)  , INTENT(IN)    :: YDGP5     ! Model variables at observation locations
INTEGER(KIND=JPIM)   , INTENT(INOUT) :: VARNOS_TO_PROCESS(:,:)
LOGICAL              , INTENT(IN)    :: LDSCREEN  ! Indicates that one-off screening and other preparations
LOGICAL              , INTENT(IN)    :: LDFINAL   ! Indicates a final analysis traj
TYPE(TPHY)           , INTENT(IN)    :: YDPHY
TYPE(TYPE_SET_INFO)  , INTENT(IN)    :: YDSET
REAL(KIND=JPRB)      , INTENT(INOUT) :: PHOFX(YDGP5%NDLEN,YDSET%MXBDY)
TYPE(TYPE_GOM_PLUS)  , INTENT(IN)   , OPTIONAL :: YDGP_TL   ! model variables - TL
TYPE(TYPE_GOM_PLUS)  , INTENT(INOUT), OPTIONAL :: YDGP_AD   ! model variables - adjoint

!   ZINTH   :   Interpolation height
!   ZINTHUV :   Interpolation height, for surface wind data
!   ZOROGOB :   Elevation*RG of observation
REAL(KIND=JPRB), DIMENSION(YDGP5%NDLEN) :: ZINTHUV,ZOROGOB,ZCMALT,ZT2MADJ
REAL(KIND=JPRB), ALLOCATABLE :: ZINTH(:,:)

!   ZVERTP  : Vertical position of obs datum (Pressure/Height/Channel)
!   ZVERTP_L: Vertical position of obs datum (Lower level)
REAL(KIND=JPRD), DIMENSION(YDGP5%NDLEN,YDSET%MXBDY) :: ZVERTP,ZVERTP_L

!   IPNLV   : Pointer to locations in PHOFX
INTEGER(KIND=JPIM), DIMENSION(YDGP5%NDLEN,YDSET%MXBDY) :: IPNLV

!   ICMBDY  :  Number of data entries in each report in this KSET
!   ICOUNT  :  Number of data per obs, for call to PP-routines
INTEGER(KIND=JPIM), DIMENSION(YDGP5%NDLEN) :: ICMBDY,ICOUNT,ICOUNTA

REAL(KIND=JPRB), ALLOCATABLE :: ZXPP(:,:,:),ZVPOBS(:,:),ZVPOBS_L(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZXPP5(:,:,:)

!   ZUV  : U and V on model or RTTOV pressure levels
!   ZW   : Weights used for mean layer winds
!   ZTAU : Transmittances on RTTOV pressure levels
REAL(KIND=JPRB), ALLOCATABLE :: ZUV(:,:,:), ZW(:,:), ZTAU(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZUV5(:,:,:), ZW5(:,:)

!   ZPRES_AK : Pressure levels for averaging kernel calculations
!   ZWEIGHT_AK: Weights for averaging kernel calculations
!   ZPRIOR_AK: Prior profile for averaging kernel calculations
!   ZPAK    : List of pressures for averaging kernel calculations
!   ZWAK    : List of weights for averaging kernel calculations
!   ZAPAK   : List of prior values for averaging kernel calculations
REAL(KIND=JPRB), ALLOCATABLE :: ZPRES_AK(:,:,:),ZWEIGHT_AK(:,:,:),ZPRIOR_AK(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZPAK(:,:,:),ZWAK(:,:,:),ZAPAK(:,:,:)

!   IAK    : number of elements for averaging kernel calculations
INTEGER(KIND=JPIM), ALLOCATABLE :: INAK(:,:), IAK(:,:)

! Indices of first element of ODB array columns
INTEGER(KIND=JPIM) :: IDX_PAK_1
INTEGER(KIND=JPIM) :: IDX_WAK_1
INTEGER(KIND=JPIM) :: IDX_APAK_1

CHARACTER (LEN = 10) :: CLV
INTEGER(KIND=JPIM) :: IVARL,IVARNO_INDEX,IBODY,IDSTA
INTEGER(KIND=JPIM) :: IMXBDY,IMXCOUNT,IBDY,IDF,IUPTRA

INTEGER(KIND=JPIM) :: JOBS,JBODY,JVNM,JNLV
REAL(KIND=JPRB) :: Z500HPA,ZLAYL(2)
REAL(KIND=JPRB) :: ZHLS,ZAZIM,ZSIN_AZIM,ZCOS_AZIM
LOGICAL :: LLAPACH,LL_LAYER,LL_WIND, LL_HLS, &
 & LL_SFC_OP,LL_SPEED,LL_MEAN_UV,&
 & LL_AIREP, LL_NONLINEAR,LLREO3_AK, LL_AK, LL_AK_MR

TYPE (TCST) :: YLCST

LOGICAL :: LLFIRST_AIREP
LOGICAL :: LLTL, LLAD, LLDIRECT

REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "append_num.h"

#include "abor1.intfb.h"
#include "calver.intfb.h"
#include "hinth.intfb.h"
#include "meanuv_average.intfb.h"
#include "meanuv_averagetl.intfb.h"
#include "meanuv_averagead.intfb.h"
#include "meanuv_weights.intfb.h"
#include "meanuv_weightstl.intfb.h"
#include "meanuv_weightsad.intfb.h"
#include "ppobsac.intfb.h"
#include "ppobsacad.intfb.h"
#include "ppobsactl.intfb.h"
#include "ppobsap.intfb.h"
#include "ppnew.intfb.h"
#include "ppobsn.intfb.h"
#include "map_varno_to_nvar.intfb.h"
#include "obsop_varno_subset.intfb.h"
#include "select_closest_scatt.intfb.h"
#include "hretr_conv.intfb.h"
#include "reo3_ak_op.intfb.h"
#include "reo3_ak_tl.intfb.h"
#include "reo3_ak_ad.intfb.h"
#include "reo3_ak_mr_op.intfb.h"
#include "reo3_ak_mr_tl.intfb.h"
#include "reo3_ak_mr_ad.intfb.h"

!     -------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OBSOP_CONV',0,ZHOOK_HANDLE)

YLCST = TCST_INIT ()

LLTL=PRESENT(YDGP_TL)
LLAD=PRESENT(YDGP_AD)
LLDIRECT=.NOT.(LLTL.OR.LLAD)
IUPTRA = GET_NUPTRA()
LL_AK = .FALSE.
LLREO3_AK = .FALSE.
LL_AK_MR = .TRUE. ! Averaging Kernels are related to Mixing Ratios (Set it to .FALSE. if AK are related to Partial Columns)

! Once-only preparations and screening
IF(LDSCREEN.OR.LCANARI) CALL HRETR_CONV(ROBHDR,ROBODY,SATHDR,SATBODY,ROBSU,YDSET,YDGP5)

Z500HPA=50000._JPRB

! Flags
LLAPACH=LCAPACH .AND. (YDSET%OBSTYPE == NTEMP .OR. YDSET%OBSTYPE == NPILOT)
LL_AIREP = YDSET%OBSTYPE == NAIREP
LL_SFC_OP=ANY(YDSET%OBSTYPE == (/NSYNOP,NDRIBU,NTEMP,NPILOT,NPAOB,NSCATT/))

IF(YDSET%OBSTYPE == NSATOB) THEN
  LL_MEAN_UV  = TRIM(SATOBGRP_TABLE(YDSET%GROUP)%OBS_OPER) == 'Gauss average'&
    & .OR.TRIM(SATOBGRP_TABLE(YDSET%GROUP)%OBS_OPER) == 'Boxcar average'
ELSE
  LL_MEAN_UV = .FALSE.
ENDIF

! Orography
DO JOBS = 1,YDGP5%NDLEN
  ICMBDY(JOBS)=ROBODY%REPORT_END_ROW(JOBS)+1-ROBODY%REPORT_START_ROW(JOBS)
ENDDO

IMXBDY=MAXVAL(ICMBDY)

ZCMALT  = ROBHDR%COPY_COLUMN_REAL("stalt@hdr")
ZOROGOB = ZCMALT * RG
WHERE (ZCMALT == RMDI .OR. ABS(ZCMALT) >  10000._JPRB) ZOROGOB = RMDI

! Extract vertical positions i.e. pressure/height/channel
! (ZVERTP) These will
! constitude the 'request' for PP-calculations, later.
IPNLV(:,1:IMXBDY)=0
ZVERTP(:,1:IMXBDY)=RMDI
ZVERTP_L(:,1:IMXBDY)=RMDI

CALL ROBODY%FILL_2D(ZVERTP,"vertco_reference_1@body",RMDI=REAL(RMDI,JPRD))
CALL ROBODY%FILL_2D(ZVERTP_L,"vertco_reference_2@body",RMDI=REAL(RMDI,JPRD))

!* Some problematic observations are ignored by setting VARNOS_TO_PROCESS=0

IF (LLDIRECT) THEN

  DO JBODY=1,IMXBDY
    DO JOBS = 1,YDGP5%NDLEN

      IF(JBODY > ICMBDY(JOBS)) CYCLE

      IF (.NOT.LCANARI.AND.VARNOS_TO_PROCESS(JOBS,JBODY)==VARNO%SFALL) VARNOS_TO_PROCESS(JOBS,JBODY)=0

      ! Missing vertical position
      IF(ZVERTP(JOBS,JBODY) == RMDI) VARNOS_TO_PROCESS(JOBS,JBODY)=0
      IF(ZVERTP(JOBS,JBODY) <= 0. .AND. VARNOS_TO_PROCESS(JOBS,JBODY) /= VARNO%PS)&
        & VARNOS_TO_PROCESS(JOBS,JBODY)=0

      ! Missing lower vertical position
      IF(YDSET%OBSTYPE /= NISAC) THEN
        IF(VARNOS_TO_PROCESS(JOBS,JBODY) == VARNO%DZ .OR.&
          & VARNOS_TO_PROCESS(JOBS,JBODY) == VARNO%RHLAY .OR.&
          & VARNOS_TO_PROCESS(JOBS,JBODY) == VARNO%O3LAY) THEN
          IF(ZVERTP_L(JOBS,JBODY) == RMDI .OR. ZVERTP_L(JOBS,JBODY) <= 0.)&
            & VARNOS_TO_PROCESS(JOBS,JBODY)=0
        ENDIF
      ENDIF

      ! Discard thickness data above model's top
      IF(VARNOS_TO_PROCESS(JOBS,JBODY) == VARNO%DZ .AND.&
        & ZVERTP(JOBS,JBODY) < YDGP5%PRESF(JOBS,1,IH)) VARNOS_TO_PROCESS(JOBS,JBODY)=0
    ENDDO

    ! SCAT requires at least two ambiguous winds (if processed)
    IF( .NOT. LENKF) THEN
      IF(YDSET%OBSTYPE == NSCATT)&
        & WHERE(VARNOS_TO_PROCESS(:,JBODY) == VARNO%SCATU &
        & .AND. ICMBDY < 4) VARNOS_TO_PROCESS(:,JBODY)=0
    ENDIF

    ! To reproduce current operations we must not compute
    ! R/S T2m departures, not to affect 'Combined flagging'
    IF (LECMWF) THEN
      IF(YDSET%OBSTYPE == NTEMP)&
        & WHERE(VARNOS_TO_PROCESS(:,JBODY) == VARNO%T2M) VARNOS_TO_PROCESS(:,JBODY)=0
    ELSEIF(NCONF == 131) THEN
      IF(YDSET%OBSTYPE == NTEMP .OR. YDSET%OBSTYPE == NPILOT)&
        & WHERE(VARNOS_TO_PROCESS(:,JBODY) == VARNO%RH2M .OR.&
        & VARNOS_TO_PROCESS(:,JBODY) == VARNO%T2M) VARNOS_TO_PROCESS(:,JBODY)=0
    ENDIF

    ! Discard ambiguous winds in CANARI
    IF (LCANARI) THEN
      WHERE(&
        & VARNOS_TO_PROCESS(:,JBODY) == VARNO%SCATV .OR.&
        & VARNOS_TO_PROCESS(:,JBODY) == VARNO%SCATU) VARNOS_TO_PROCESS(:,JBODY)=0
    ENDIF

  ENDDO

ELSEIF (LLTL) THEN

  DO JBODY=1,IMXBDY
    DO JOBS = 1,YDGP5%NDLEN
      IF (VARNOS_TO_PROCESS(JOBS,JBODY)==VARNO%SDEPTH .OR. &
        & VARNOS_TO_PROCESS(JOBS,JBODY)==VARNO%SFALL .OR. &
        & VARNOS_TO_PROCESS(JOBS,JBODY)==VARNO%RR) VARNOS_TO_PROCESS(JOBS,JBODY)=0
    ENDDO
  ENDDO

ENDIF


!* Surface wind data : Find interpolation heights

IF(LL_SFC_OP) CALL HINTH(ROBHDR,ROBSU,&
  & YDGP5%NDLEN,YDGP5%NDLEN,YDSET%MXBDY,IMXBDY,ICMBDY,YDSET%OBSTYPE,&
  & YDGP5%APHI(:,YDGP5%NFLEVG,IH),ZINTHUV,VARNOS_TO_PROCESS)

IF(YDSET%OBSTYPE == NSCATT) THEN
!* Quadratic SCATT cost function requires that one of the two
!   ambiguous winds is chosen a priori.
!   Use the one closest to the high resolution trajectory.
!   For the trajectory itself (apart from last one) we need to compute H(x) for all ambiguous winds

  IF( .NOT.L_OBS_IN_FC() ) THEN
  !IF(NCONF == 131) THEN !CP would be more adequate?

    DO JBODY=1,IMXBDY
      DO JOBS = 1,YDGP5%NDLEN
        IF(JBODY > ICMBDY(JOBS)) CYCLE
        IBODY = ROBODY%REPORT_START_ROW(JOBS)+JBODY-1
        IF(.NOT.ANY(ROBODY%BODY%VARNO(IBODY) == &
         & (/VARNO%SCATU,VARNO%SCATV,VARNO%U10M,VARNO%V10M/))) THEN ! Only use ambiguous winds or de-aliased solution  
          VARNOS_TO_PROCESS(JOBS,JBODY) = 0
        ELSE
          IF (.NOT.ANY(ROBODY%BODY%VARNO(IBODY) ==  (/VARNO%SCATU,VARNO%SCATV/))) CYCLE
          ! Check bit corresponding to current outer loop number in ambig_select, bit set in SELECT_CLOSEST_SCATT in outer loop
          IF(IBITS(NINT(SCATTBODY%SCATT_BODY%AMBIG_SELECT(IBODY)),IUPTRA,1) /= 1)  THEN ! Only use closest ambiguous wind
            VARNOS_TO_PROCESS(JOBS,JBODY) = 0
          ENDIF
        ENDIF
      ENDDO
    ENDDO

  ENDIF
ENDIF


!* loop on VARNOs
LLFIRST_AIREP=.TRUE.
DO JVNM=1,NVNUMAX

  CALL MAP_VARNO_TO_NVAR(JVNM,CLV,IVARL,IVARNO_INDEX)
  CALL OBSOP_VARNO_SUBSET(VARNOS_TO_PROCESS,YDGP5%NDLEN,IMXBDY,JVNM,ICOUNT,IMXCOUNT,IPNLV)
  IF(IMXCOUNT == 0) CYCLE

  IF(JVNM==VARNO%V .OR. JVNM==VARNO%V10M .OR. JVNM==VARNO%SCATV) CYCLE
!F.S add the following line because otherwise we erase the work done by obsop_radar
!in hop.F90 (RMDI is put into PHOFX)
  IF(JVNM==VARNO%DOPP .OR. JVNM==VARNO%REFL .OR. JVNM==VARNO%RAWBT) CYCLE !radar and bayrad)

  LL_WIND = JVNM==VARNO%U   .OR. JVNM==VARNO%U10M .OR. JVNM==VARNO%SCATU
  LL_SPEED= JVNM==VARNO%FF  .OR. JVNM==VARNO%SCATWS
!CP ll_speed definition (and name) misleading as only surface speed is treated in practice
!CP  or something that I do not understand if altitude
  LL_NONLINEAR = JVNM==VARNO%RH2M .OR. LL_SPEED
  LL_HLS = JVNM==VARNO%LOS

  ! Identify and deal with layer quantities
  LL_LAYER = ( JVNM == VARNO%DZ    .OR. JVNM == VARNO%PWC .OR.&
             & JVNM == VARNO%RHLAY .OR. JVNM == VARNO%O3LAY )
  IF(LL_LAYER) THEN
    IF (JVNM == VARNO%RHLAY) THEN
      ZLAYL(1:2)=(/0.5_JPRB,0.5_JPRB/)
    ELSEIF(JVNM == VARNO%DZ) THEN
      ZLAYL(1:2)=(/1.0_JPRB,-1.0_JPRB/)
    ELSE
      ZLAYL(1:2)=(/-1.0_JPRB,1.0_JPRB/)
    ENDIF
        
    !********* SQR,RD ********
      
    IF (JVNM == VARNO%O3LAY) THEN
      ALLOCATE(ZPRES_AK(YDGP5%NDLEN,YDSET%MXBDY,JPMX_AK),&
       & ZWEIGHT_AK(YDGP5%NDLEN,YDSET%MXBDY,JPMX_AK),&
       & ZPRIOR_AK(YDGP5%NDLEN,YDSET%MXBDY,JPMX_AK),&
       & INAK(YDGP5%NDLEN,YDSET%MXBDY))
      INAK(:,1:YDSET%MXBDY)=0
      ZPRES_AK(:,1:YDSET%MXBDY,1:JPMX_AK)=RMDI
      ZWEIGHT_AK(:,1:YDSET%MXBDY,1:JPMX_AK)=0.0_JPRB
      ZPRIOR_AK(:,1:YDSET%MXBDY,1:JPMX_AK)=0.0_JPRB


      ! Get index to first element of each required ODB array column
      IDX_PAK_1  = SATBODY%GET_COLUMN_INDEX("pak_1@resat_averaging_kernel")
      IDX_WAK_1  = SATBODY%GET_COLUMN_INDEX("wak_1@resat_averaging_kernel")
      IDX_APAK_1 = SATBODY%GET_COLUMN_INDEX("apak_1@resat_averaging_kernel")
   
      DO JBODY=1,ROBODY%MAX_DATUM_PER_REPORT
        DO JOBS = 1,YDGP5%NDLEN
          IF(JBODY > ROBODY%REPORT_END_ROW(JOBS)+1-ROBODY%REPORT_START_ROW(JOBS)) CYCLE 
          IBODY = ROBODY%REPORT_START_ROW(JOBS)+(JBODY-1)
    
          LLREO3_AK = (SATHDR%HDR%RETRTYPE(JOBS) == 1)
          IF (LLREO3_AK) THEN
            INAK(JOBS,JBODY)  = SATBODY%RESAT_AVERAGING_KERNEL%NAK(IBODY)
	      
            !IF(INAK(JOBS,JBODY) > 0 .AND. ABS(INAK(JOBS,JBODY)) < 1000) THEN    ! As written in obsop_composition
            IF(INAK(JOBS,JBODY) > 0 .AND. ABS(INAK(JOBS,JBODY)) <= JPMX_AK) THEN           
              IF (.NOT.LL_AK) LL_AK=.TRUE.
              ZPRES_AK(JOBS,JBODY,1:INAK(JOBS,JBODY))  =&
              & SATBODY%DATA(IBODY,IDX_PAK_1:IDX_PAK_1+INAK(JOBS,JBODY)-1)
              ZWEIGHT_AK(JOBS,JBODY,1:INAK(JOBS,JBODY))  =&
              & SATBODY%DATA(IBODY,IDX_WAK_1:IDX_WAK_1+INAK(JOBS,JBODY)-1)
              ZPRIOR_AK(JOBS,JBODY,1:INAK(JOBS,JBODY))  =&
              & SATBODY%DATA(IBODY,IDX_APAK_1:IDX_APAK_1+INAK(JOBS,JBODY)-1)
            ENDIF
          ENDIF
   
        ENDDO 
      ENDDO 

      ! One-per varno allocations and fills
      ALLOCATE(IAK(YDGP5%NDLEN,IMXCOUNT))
      IAK(:,:)=0.0_JPRB
      ALLOCATE(ZPAK(YDGP5%NDLEN,IMXCOUNT,JPMX_AK))
      ZPAK(:,:,:)=50000.0_JPRB
      ALLOCATE(ZWAK(YDGP5%NDLEN,IMXCOUNT,JPMX_AK))
      ZWAK(:,:,:)=0.0_JPRB
      ALLOCATE(ZAPAK(YDGP5%NDLEN,IMXCOUNT,JPMX_AK))
      ZAPAK(:,:,:)=0.0_JPRB

      DO JNLV=1,IMXCOUNT
        DO JOBS=1,YDGP5%NDLEN
          IF(JNLV <= ICOUNT(JOBS)) THEN
            IBDY=IPNLV(JOBS,JNLV)
 
            IAK(JOBS,JNLV)=INAK(JOBS,IBDY)
            !IF(IAK(JOBS,JNLV) > 0) THEN       ! As written in obsop_composition
            IF(IAK(JOBS,JNLV) > 0  .AND. IAK(JOBS,JNLV) <= JPMX_AK) THEN
              ZPAK(JOBS,JNLV,1:IAK(JOBS,JNLV)) =&
                  & ZPRES_AK(JOBS,IBDY,1:INAK(JOBS,IBDY))
              ZWAK(JOBS,JNLV,1:IAK(JOBS,JNLV)) =&
                  & ZWEIGHT_AK(JOBS,IBDY,1:INAK(JOBS,IBDY))
              ZAPAK(JOBS,JNLV,1:IAK(JOBS,JNLV)) =&
                  & ZPRIOR_AK(JOBS,IBDY,1:INAK(JOBS,IBDY))
            ENDIF
          ENDIF
        ENDDO
      ENDDO
   
    ENDIF   ! JVNM == VARNO%O3LAY    
    
        
  ENDIF

  ! Arrays to hold PP-values (ZXPP) and contracted lists of
  ! vertical positions (ZVPOBS, ZVPOBS_L)
  IDF=2
  ALLOCATE(ZXPP(YDGP5%NDLEN,IMXCOUNT,IDF))
  IF(.NOT.LLAD)THEN
    ZXPP=RMDI
  ELSE
    ZXPP=0.0_JPRB
  ENDIF
  ALLOCATE(ZINTH(YDGP5%NDLEN,IMXCOUNT))
  ALLOCATE(ZVPOBS(YDGP5%NDLEN,IMXCOUNT))
  ALLOCATE(ZVPOBS_L(YDGP5%NDLEN,IMXCOUNT))
  ZVPOBS(:,:)=Z500HPA
  ZVPOBS_L(:,:)=Z500HPA

  IF (.NOT.LLDIRECT) THEN
    ALLOCATE(ZXPP5(YDGP5%NDLEN,IMXCOUNT,IDF))
    ZXPP5=RMDI
  ENDIF

  DO JNLV=1,IMXCOUNT
    DO JOBS=1,YDGP5%NDLEN
      IF(JNLV <= ICOUNT(JOBS)) THEN
        IBDY=IPNLV(JOBS,JNLV)
        ZVPOBS(JOBS,JNLV)=ZVERTP(JOBS,IBDY)
        ZVPOBS_L(JOBS,JNLV)=ZVERTP_L(JOBS,IBDY)
      ENDIF
    ENDDO
  ENDDO

  !* AD preliminary computations
  ICOUNTA = ICOUNT
  IF (LLAD) THEN

    !  Some non-linear AD operators require PP5
    IF(LL_NONLINEAR) THEN
      IF(LL_WIND .OR. LL_SPEED) THEN
        DO JNLV=1,IMXCOUNT
          ZINTH(:,JNLV) = ZINTHUV(:)
        ENDDO
      ELSE
        ZINTH(:,:) =2.0_JPRB*RG
      ENDIF
      IF(LCACHMT) THEN
        CALL PPOBSAC (YLCST, YDGP5,IMXCOUNT,CLV,ZINTH(:,1),YDSET%OBSTYPE,&
          & ICOUNT,ZOROGOB,IDF,ZXPP5,.FALSE.)
      ENDIF
    ENDIF

    ! Transfer PHOFX to their PP-values location
    ! Design note: "if" tests are knowingly kept inside this loop to preserve readability/maintainability
    DO JNLV=1,IMXCOUNT
      DO JOBS=1,YDGP5%NDLEN
        IF(JNLV > ICOUNT(JOBS)) CYCLE
        IBDY=IPNLV(JOBS,JNLV)
        IBODY = ROBODY%REPORT_START_ROW(JOBS)+(IBDY-1)

        IF(LL_WIND.OR.LL_SPEED) THEN
          IF (ROBHDR%HDR%CODETYPE(JOBS) == NMDEHS ) THEN
            ! Do not compute TL/AD for modes passive data 
            ! Get status of data
            IF (JNLV == 1) ICOUNTA(JOBS)=0
            IDSTA = IBITS(INT(ROBODY%BODY%DATUM_STATUS(IBODY)),NCMDPABP, NCMDPAOC)
            IF (IDSTA /= 1) THEN
              ICOUNTA(JOBS) = ICOUNTA(JOBS) + 1
              ! Send adjoint TB into the adjoint obsop
              IF(LLAD) THEN
                ZXPP(JOBS,ICOUNTA(JOBS),1) = PHOFX(JOBS,IBDY)
                ZXPP(JOBS,ICOUNTA(JOBS),2) = PHOFX(JOBS,IBDY+1)
              ENDIF
            ENDIF
           ELSE
           ZXPP(JOBS,JNLV,1) = PHOFX(JOBS,IBDY)
           ZXPP(JOBS,JNLV,2) = PHOFX(JOBS,IBDY+1)
           ENDIF
        ELSEIF(LL_HLS) THEN
          ! Horizontal line-of-sight wind operator adjoint
          ! Used for Aeolus Doppler Wind Lidar
          ZAZIM=SATHDR%SAT%AZIMUTH(JOBS)
          ZSIN_AZIM=SIN(ZAZIM)
          ZCOS_AZIM=COS(ZAZIM)
          ZXPP(JOBS,JNLV,1) = -PHOFX(JOBS,IBDY)*ZSIN_AZIM
          ZXPP(JOBS,JNLV,2) = -PHOFX(JOBS,IBDY)*ZCOS_AZIM
          PHOFX(JOBS,IBDY)=0.0_JPRB
        ELSEIF(LL_LAYER .AND. (.NOT. LL_AK)) THEN
          ZXPP(JOBS,JNLV,1) = ZLAYL(1)*PHOFX(JOBS,IBDY)
          ZXPP(JOBS,JNLV,2) = ZLAYL(2)*PHOFX(JOBS,IBDY)
          PHOFX(JOBS,IBDY)=0.0_JPRB
        ELSE
          ZXPP (JOBS,JNLV,1) = PHOFX(JOBS,IBDY)
        ENDIF

      ENDDO
    ENDDO

  ENDIF !end LLAD (preliminary computations)

  !*  Upper air operators

  IF( (JVNM==VARNO%U .AND. .NOT.LL_MEAN_UV)  .OR.&
    & JVNM==VARNO%LOS.OR.&
    & JVNM==VARNO%T  .OR. JVNM==VARNO%RH  .OR. (&
    & JVNM==VARNO%FF .AND. YDSET%OBSTYPE /= NSCATT) .OR.&
    & JVNM==VARNO%RHLAY .OR. JVNM==VARNO%PWC.OR. JVNM==VARNO%Q  .OR.&
    & JVNM==VARNO%Z  .OR. JVNM==VARNO%DZ .OR. & 
    & (JVNM==VARNO%O3LAY .AND. (.NOT. LL_AK)) .OR.& ! SQR/RD: split O3 with AKs from the rest
    & JVNM==VARNO%SATCL ) THEN

    IF (.NOT.(LLDIRECT.AND.LLAPACH)) THEN
      IF (LLTL .OR. LLAD) THEN
!      WRITE(NULOUT,*) 'ADJOINT ICOUNT = ',ICOUNT(:),' ICOUNTA ',ICOUNTA(:) 
       CALL PPNEW(YDGP5,ZVPOBS,ICOUNTA,ZXPP,CDPAR=CLV,YDGP_TL=YDGP_TL,YDGP_AD=YDGP_AD)
      ELSE
!      WRITE(NULOUT,*) 'DIRECT ICOUNT = ',ICOUNT(:),' ICOUNTA ',ICOUNTA(:)       
       CALL PPNEW(YDGP5,ZVPOBS,ICOUNT,ZXPP,CDPAR=CLV,YDGP_TL=YDGP_TL,YDGP_AD=YDGP_AD)
      ENDIF
      IF(LL_LAYER) CALL PPNEW(YDGP5,ZVPOBS_L,ICOUNT,ZXPP(:,:,2:2),CDPAR=CLV,YDGP_TL=YDGP_TL,YDGP_AD=YDGP_AD)
    ELSE
      CALL PPOBSAP(YDGP5,ROBODY,YDSET%MXBDY,IMXCOUNT,IPNLV,CLV,ICOUNT,ZVPOBS,ZOROGOB,IDF,ZXPP,YDPHY)
    ENDIF




  !* SQR/RD: Ozone with AKs
  ELSEIF(JVNM==VARNO%O3LAY .AND. LL_AK) THEN


    !* THIS IS A PATCH TO AVOID HAVING OBSERVATION PRESSURES BELOW THE SURFACE MODEL PRESSURE
    !* START PATCH
    DO JNLV=1,IMXCOUNT
      DO JOBS=1,YDGP5%NDLEN
        IF(IAK(JOBS,JNLV) > 0  .AND. IAK(JOBS,JNLV) <= JPMX_AK) THEN
          ZPAK(JOBS,JNLV,IAK(JOBS,JNLV)) = MIN( ZPAK(JOBS,JNLV,IAK(JOBS,JNLV)), YDGP5%PRESH(JOBS,YDGP5%NFLEVG,IH)  ) 
        ENDIF
      ENDDO
    ENDDO
    !* END PATCH
        
    IF (LLDIRECT) THEN

      IF (LL_AK_MR) THEN 
        !* REO3_AK_MR_OP considers Averaging Kernel related to Mixing Ratios (idem for TL/AD)
        CALL REO3_AK_MR_OP(YDGP5%NFLEVG,YDGP5%NPPM,YDGP5%NDLEN,YDGP5%NDLEN,IMXCOUNT,ICOUNT,0,&
               & ZVPOBS, ZVPOBS_L,&
               & IAK,ZWAK,ZAPAK,ZPAK,YDGP5%PRESH(:,:,IH),YDGP5%PRESF(:,:,IH),YDGP5%PXP(:,:,:,IH),&
               & YDGP5%PXPD(:,:,:,IH),YDGP5,IDF,ZXPP)
      ELSE
        !* REO3_AK_OP considers Averaging Kernel related to Partial Columns (idem for TL/AD)
        CALL REO3_AK_OP(YDGP5%NFLEVG,YDGP5%NPPM,YDGP5%NDLEN,YDGP5%NDLEN,IMXCOUNT,ICOUNT,0,&
               & IAK,ZWAK,ZAPAK,ZPAK,YDGP5%PRESH(:,:,IH),YDGP5%PRESF(:,:,IH),YDGP5%PXP(:,:,:,IH),&
               & YDGP5%PXPD(:,:,:,IH),YDGP5,IDF,ZXPP)
      ENDIF

    ELSEIF (LLTL) THEN

      IF (LL_AK_MR) THEN 
        CALL REO3_AK_MR_TL(YDGP5%NFLEVG,YDGP5%NPPM,YDGP5%NDLEN,YDGP5%NDLEN,IMXCOUNT,ICOUNT,0,&
               & ZVPOBS, ZVPOBS_L,&
               & IAK,ZWAK,ZPAK,YDGP_TL%PXP(:,:,:,IH),YDGP_TL%PXPD(:,:,:,IH),YDGP_TL,&
               & IDF,ZXPP,&
               & YDGP5%PRESH(:,:,IH),YDGP5%PRESF(:,:,IH),YDGP5%PXP(:,:,:,IH),YDGP5%PXPD(:,:,:,IH),YDGP5)
      ELSE
        CALL REO3_AK_TL(YDGP5%NFLEVG,YDGP5%NPPM,YDGP5%NDLEN,YDGP5%NDLEN,IMXCOUNT,ICOUNT,0,&
               & IAK,ZWAK,ZPAK,YDGP_TL%PXP(:,:,:,IH),YDGP_TL%PXPD(:,:,:,IH),YDGP_TL,&
               & IDF,ZXPP,&
               & YDGP5%PRESH(:,:,IH),YDGP5%PRESF(:,:,IH),YDGP5%PXP(:,:,:,IH),YDGP5%PXPD(:,:,:,IH),YDGP5)
      ENDIF

    ELSEIF (LLAD) THEN

      IF (LL_AK_MR) THEN 
        CALL REO3_AK_MR_AD(YDGP5%NFLEVG,YDGP5%NPPM,YDGP5%NDLEN,YDGP5%NDLEN,IMXCOUNT,ICOUNT,0,&
               & ZVPOBS, ZVPOBS_L,&
               & IAK,ZWAK,ZPAK,YDGP_AD%PXP(:,:,:,IH),YDGP_AD%PXPD(:,:,:,IH),YDGP_AD,&
               & IDF,ZXPP,&                
               & YDGP5%PRESH(:,:,IH),YDGP5%PRESF(:,:,IH),YDGP5%PXP(:,:,:,IH),YDGP5%PXPD(:,:,:,IH),YDGP5)
      ELSE
        CALL REO3_AK_AD(YDGP5%NFLEVG,YDGP5%NPPM,YDGP5%NDLEN,YDGP5%NDLEN,IMXCOUNT,ICOUNT,0,&
               & IAK,ZWAK,ZPAK,YDGP_AD%PXP(:,:,:,IH),YDGP_AD%PXPD(:,:,:,IH),YDGP_AD,&
               & IDF,ZXPP,&                
               & YDGP5%PRESH(:,:,IH),YDGP5%PRESF(:,:,IH),YDGP5%PXP(:,:,:,IH),YDGP5%PXPD(:,:,:,IH),YDGP5)
      ENDIF

    ENDIF




  !*  Surface operators

  ELSEIF((JVNM==VARNO%U10M .OR. JVNM==VARNO%SCATU) .OR. JVNM==VARNO%T2M .OR.&
    & JVNM==VARNO%RH2M .OR. JVNM==VARNO%Q2M .OR.&
    & JVNM==VARNO%SCATWS .OR. JVNM==VARNO%TS .OR. (&
    & JVNM==VARNO%FF  .AND. YDSET%OBSTYPE == NSCATT)) THEN

    ! Hack to get scatterometer windspeeds treated right
    IF(JVNM==VARNO%FF) CLV='FFS'

    IF(LL_WIND .OR. LL_SPEED) THEN
      DO JNLV=1,IMXCOUNT
        ZINTH(:,JNLV) = ZINTHUV(:)
      ENDDO
    ELSE
      ZINTH(:,:) =2.0_JPRB*RG
    ENDIF

    ! AJGDB ...maybe:
    IF(.NOT.LL_SFC_OP) CALL ABOR1('OBSOP_CONV requires LL_SFC_OP=.T. to run PPNEW here')
    
    IF(.NOT.LCACHMT) THEN

      IF (L_SCREEN_LEVEL_OBSOP .AND. YDSET%OBSTYPE == NSYNOP .AND.&
         & (JVNM==VARNO%T2M .OR. JVNM==VARNO%Q2M .OR. JVNM==VARNO%U10M)) THEN
        IF (LLDIRECT) THEN
          ZT2MADJ = 0.0
          WHERE (ZCMALT /= RMDI) ZT2MADJ = (ZCMALT-YDGP5%OROG(:,IH)/RG)*R4DLAPSET2
          IF (JVNM==VARNO%T2M) ZXPP(:,1,1) = YDGP5%T2M(:,IH) - ZT2MADJ(:)
          IF (JVNM==VARNO%Q2M) ZXPP(:,1,1) = YDGP5%Q2M(:,IH)
          IF (JVNM==VARNO%U10M) ZXPP(:,1,1) = YDGP5%U10M(:,IH)
          IF (JVNM==VARNO%U10M) ZXPP(:,1,2) = YDGP5%V10M(:,IH)
        ELSEIF (LLTL) THEN
          IF (JVNM==VARNO%T2M) ZXPP(:,1,1) = YDGP_TL%T2M(:,IH)
          IF (JVNM==VARNO%Q2M) ZXPP(:,1,1) = YDGP_TL%Q2M(:,IH)
          IF (JVNM==VARNO%U10M) ZXPP(:,1,1) = YDGP_TL%U10M(:,IH)
          IF (JVNM==VARNO%U10M) ZXPP(:,1,2) = YDGP_TL%V10M(:,IH)
        ELSE  ! Adjoint case
          IF (JVNM==VARNO%T2M) YDGP_AD%T2M(:,IH) = YDGP_AD%T2M(:,IH) +  ZXPP(:,1,1)
          IF (JVNM==VARNO%Q2M) YDGP_AD%Q2M(:,IH) = YDGP_AD%Q2M(:,IH) +  ZXPP(:,1,1)
          IF (JVNM==VARNO%U10M) YDGP_AD%U10M(:,IH) = YDGP_AD%U10M(:,IH) + ZXPP(:,1,1)
          IF (JVNM==VARNO%U10M) YDGP_AD%V10M(:,IH) = YDGP_AD%V10M(:,IH) + ZXPP(:,1,2)
        ENDIF
      ELSE

        CALL PPNEW(YDGP5,ZINTH,ICOUNT,ZXPP,CDPAR=CLV,&
          & LDNEUTRAL=(YDSET%OBSTYPE == NSCATT .AND. LSCATT_NEUTRAL), &
          & LDRELCURR=(YDSET%OBSTYPE == NSCATT),&
          & YDGP_TL=YDGP_TL,YDGP_AD=YDGP_AD)
      ENDIF

    ELSE
      IF (LLDIRECT) THEN

        CALL PPOBSAC (YLCST, YDGP5,IMXCOUNT,CLV,ZINTH(:,1),YDSET%OBSTYPE,&
          & ICOUNT,ZOROGOB,IDF,ZXPP,YRDPHY%LDIRCLSMOD)

      ELSEIF (LLTL) THEN

        IF(LL_NONLINEAR) CALL PPOBSAC (YLCST, YDGP5,IMXCOUNT,CLV,ZINTH(:,1),YDSET%OBSTYPE,&
          & ICOUNT,ZOROGOB,IDF,ZXPP5,.FALSE.)

       CALL PPOBSACTL(ROBHDR,YDGP5,IMXCOUNT,CLV,ZINTH(:,1),&
          & YDSET%OBSTYPE,ICOUNT,YDGP_TL,IDF,ZXPP,ZXPP5)

      ELSEIF (LLAD) THEN

        CALL PPOBSACAD(ROBHDR,YDSET,IMXCOUNT,CLV,ZINTH(:,1),&
          & ICOUNT,YDGP5,YDGP_AD,IDF,ZXPP,ZXPP5)

      ENDIF

    ENDIF

  ELSEIF(JVNM==VARNO%SFALL ) THEN

    ! Surface operators CANARI
    IF (LLDIRECT) CALL PPOBSN(YDGP5,YDGP5%NDLEN,IMXCOUNT,CLV,ICOUNT,IDF,ZXPP)

  ELSEIF(JVNM==VARNO%PS) THEN

    ! Surface pressure, with height as vertical coordinate ZVPOBS is height
    CALL PPNEW(YDGP5,ZVPOBS,ICOUNT,ZXPP,CDPAR=CLV,LDHEIGHT=.TRUE.,KSLCT=1,YDGP_TL=YDGP_TL,YDGP_AD=YDGP_AD)

  ELSEIF(JVNM==VARNO%U .AND. LL_MEAN_UV) THEN

    ! Mean layer wind (for SATOB/AMVs)
    SELECT CASE (TRIM(SATOBGRP_TABLE(YDSET%GROUP)%OBS_OPER))
      CASE ('Gauss average','Boxcar average')
        ! Function averaging
        ALLOCATE(ZUV(YDGP5%NDLEN,YDGP5%NFLEVG,2))
        ALLOCATE(ZW(YDGP5%NDLEN,YDGP5%NFLEVG))
        ALLOCATE(ZTAU(YDGP5%NDLEN,IMXCOUNT,YDGP5%NFLEVG))

        IF (LLDIRECT) THEN

          CALL MEANUV_WEIGHTS(YDGP5%NDLEN,IMXCOUNT,YDGP5%NFLEVG,&
            & SATOBGRP_TABLE(YDSET%GROUP)%OBS_OPER,ZTAU,ZVPOBS,&
            & YDGP5%PRESF(:,:,IH),SATOBGRP_TABLE(YDSET%GROUP)%SIGMA,ZW)

          ZUV(:,1:YDGP5%NFLEVG,1) = YDGP5%UF(:,1:YDGP5%NFLEVG,IH)
          ZUV(:,1:YDGP5%NFLEVG,2) = YDGP5%VF(:,1:YDGP5%NFLEVG,IH)

          CALL MEANUV_AVERAGE(YDGP5%NDLEN,IMXCOUNT,YDGP5%NFLEVG,ZUV,ZW,ZXPP)

        ELSE
          ALLOCATE(ZUV5(YDGP5%NDLEN,YDGP5%NFLEVG,2))
          ALLOCATE(ZW5(YDGP5%NDLEN,YDGP5%NFLEVG))

          IF (LLTL) THEN
            CALL MEANUV_WEIGHTSTL(YDGP5%NDLEN,YDGP5%NDLEN,IMXCOUNT,YDGP5%NFLEVG,&
              & SATOBGRP_TABLE(YDSET%GROUP)%OBS_OPER,ZTAU,ZVPOBS,&
              & YDGP5%PRESF(:,:,IH),YDGP_TL%PRESF(:,:,IH),&
              & SATOBGRP_TABLE(YDSET%GROUP)%SIGMA,ZW5,ZW)

            ZUV5(:,1:YDGP5%NFLEVG,1) = YDGP5%UF(:,1:YDGP5%NFLEVG,IH)
            ZUV5(:,1:YDGP5%NFLEVG,2) = YDGP5%VF(:,1:YDGP5%NFLEVG,IH)
            ZUV(:,1:YDGP5%NFLEVG,1) = YDGP_TL%UF(:,1:YDGP5%NFLEVG,IH)
            ZUV(:,1:YDGP5%NFLEVG,2) = YDGP_TL%VF(:,1:YDGP5%NFLEVG,IH)

            CALL MEANUV_AVERAGETL(YDGP5%NDLEN,YDGP5%NDLEN,IMXCOUNT,YDGP5%NFLEVG,&
              & ZUV5,ZUV,ZW5,ZW,ZXPP)

          ELSEIF (LLAD) THEN
            ZUV = 0.0_JPRB
            ZW = 0.0_JPRB

            CALL MEANUV_WEIGHTS(YDGP5%NDLEN,IMXCOUNT,YDGP5%NFLEVG,&
              & SATOBGRP_TABLE(YDSET%GROUP)%OBS_OPER,ZTAU,ZVPOBS,&
              & YDGP5%PRESF(:,:,IH),SATOBGRP_TABLE(YDSET%GROUP)%SIGMA,ZW5)

            ZUV5(:,1:YDGP5%NFLEVG,1) = YDGP5%UF(:,1:YDGP5%NFLEVG,IH)
            ZUV5(:,1:YDGP5%NFLEVG,2) = YDGP5%VF(:,1:YDGP5%NFLEVG,IH)

            CALL MEANUV_AVERAGEAD(YDGP5%NDLEN,YDGP5%NDLEN,IMXCOUNT,YDGP5%NFLEVG,&
              & ZUV5,ZUV,ZW5,ZW,ZXPP)

            YDGP_AD%UF(:,1:YDGP5%NFLEVG,IH) = ZUV(:,1:YDGP5%NFLEVG,1)
            YDGP_AD%VF(:,1:YDGP5%NFLEVG,IH) = ZUV(:,1:YDGP5%NFLEVG,2)

            CALL MEANUV_WEIGHTSAD(YDGP5%NDLEN,YDGP5%NDLEN,IMXCOUNT,YDGP5%NFLEVG,&
              & SATOBGRP_TABLE(YDSET%GROUP)%OBS_OPER,ZTAU,ZVPOBS,&
              & YDGP5%PRESF(:,:,IH),YDGP_AD%PRESF(:,:,IH),&
              & SATOBGRP_TABLE(YDSET%GROUP)%SIGMA,ZW5,ZW)
          ENDIF
          DEALLOCATE(ZUV5,ZW5)
        ENDIF
        DEALLOCATE(ZUV,ZW,ZTAU)

      CASE DEFAULT
        CALL ABOR1("Unimplemented observation operator choice for AMVs")
    END SELECT

  ELSEIF(JVNM==VARNO%SOILM) THEN

    IF (LLDIRECT) THEN
      ! ASCAT Soil Moisture
      DO JNLV=1,IMXCOUNT
        DO JOBS=1,YDGP5%NDLEN
          IF(JNLV > ICOUNT(JOBS)) CYCLE
          ! nullify soil moisture if Temperature is below 0 or snow is present
          ! under these conditions soil moisture is not correct
          ! AJGDB this looks like a bug as the near-surface level is at YDGP5%NFLEVG.
          ! This code is checking TOA temperature...
          IF ((YDGP5%TF(JOBS,1,IH) > 273.15).AND.(YDGP5%SN(JOBS,IH) < 0.1)) THEN
            ZXPP(JOBS,JNLV,1)=MAX(1.E-10_JPRB,YDGP5%WS(JOBS,IH))
          ENDIF
        ENDDO
      ENDDO
    ENDIF

  ENDIF

  !*           CANARI CALLS CALVER

  IF(LCANARI) THEN
    IF(JVNM==VARNO%U  .OR. JVNM==VARNO%T  .OR. JVNM==VARNO%RH  .OR.(&
      & JVNM==VARNO%U10M .OR. JVNM==VARNO%SCATU).OR. JVNM==VARNO%T2M .OR. JVNM==VARNO%RH2M .OR.&
      & JVNM==VARNO%Z  .OR. JVNM==VARNO%TS .OR. JVNM==VARNO%SFALL .OR.&
      & JVNM==VARNO%RR .OR. JVNM==VARNO%Q  .OR. JVNM==VARNO%PWC ) THEN

       CALL CALVER(ROBHDR,ROBODY,YDGP5,YDGP5%NDLEN,YDSET%MXBDY,YDSET%ID,IMXCOUNT,&
                 & IPNLV,ZVPOBS,ICOUNT,YDGP5%PRESF(:,:,IH),0)
    ENDIF
  ENDIF

  ! Transfer PP-values (i.e. the model equivalents of the observations)
  ! to their location in PHOFX
  ! Design note: "if" tests are knowingly kept inside this loop to preserve readability/maintainability
  IF (.NOT.LLAD) THEN

    DO JNLV=1,IMXCOUNT
      DO JOBS=1,YDGP5%NDLEN
        IF (JNLV > ICOUNT(JOBS)) CYCLE
        IBDY=IPNLV(JOBS,JNLV)

        IF(LL_WIND) THEN
          PHOFX(JOBS,IBDY  )=ZXPP(JOBS,JNLV,1)
          PHOFX(JOBS,IBDY+1)=ZXPP(JOBS,JNLV,2)
        ELSEIF(JVNM==VARNO%LOS) THEN
          ! Horizontal line-of-sight wind operator
          ! Used for Aeolus Doppler Wind Lidar
          ZAZIM=SATHDR%SAT%AZIMUTH(JOBS)
          ZSIN_AZIM=SIN(ZAZIM)
          ZCOS_AZIM=COS(ZAZIM)
          ZHLS=-ZXPP(JOBS,JNLV,1)*ZSIN_AZIM-ZXPP(JOBS,JNLV,2)*ZCOS_AZIM
          PHOFX(JOBS,IBDY)=ZHLS

          ! Fill in the Aeolus L2C product information
          IF (LDSCREEN.AND..NOT.LLTL) THEN  !first trajectory
            SATHDR%AEOLUS_L2C%U_FG(JOBS)=ZXPP(JOBS,JNLV,1) !u component
            SATHDR%AEOLUS_L2C%V_FG(JOBS)=ZXPP(JOBS,JNLV,2) !v component
            SATHDR%AEOLUS_L2C%HLOS_FG(JOBS)=ZHLS !HLOS
          ENDIF

          IF (.NOT.L_OOPS.AND.LDFINAL.AND..NOT.LLTL) THEN  !final trajectory
            ! OOPS runs a "final" trajectory with only CCMA, hence causing
            ! problems for filling in Aeolus L2C data
            SATHDR%AEOLUS_L2C%U_AN(JOBS)=ZXPP(JOBS,JNLV,1) !u component
            SATHDR%AEOLUS_L2C%V_AN(JOBS)=ZXPP(JOBS,JNLV,2) !v component
            SATHDR%AEOLUS_L2C%HLOS_AN(JOBS)=ZHLS !HLOS
          ENDIF

        ELSEIF(LL_LAYER .AND. (.NOT. LL_AK)) THEN
          PHOFX(JOBS,IBDY) = ZLAYL(1)*ZXPP(JOBS,JNLV,1) &
                          & +ZLAYL(2)*ZXPP(JOBS,JNLV,2)
        ELSE
          PHOFX(JOBS,IBDY)=ZXPP(JOBS,JNLV,1)
        ENDIF

      ENDDO
    ENDDO
  ENDIF

  ! Deallocate
  DEALLOCATE(ZXPP,ZINTH,ZVPOBS,ZVPOBS_L)
  IF (.NOT.LLDIRECT) DEALLOCATE(ZXPP5)

ENDDO

IF(YDSET%OBSTYPE == NSCATT) THEN
  IF(L_OBS_IN_FC()) THEN
    CALL SELECT_CLOSEST_SCATT(ROBODY,SCATTBODY,YDGP5,YDSET,PHOFX)
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('OBSOP_CONV',1,ZHOOK_HANDLE)

END SUBROUTINE OBSOP_CONV
