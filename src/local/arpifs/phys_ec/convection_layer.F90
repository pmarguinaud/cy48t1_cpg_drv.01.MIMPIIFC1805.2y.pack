SUBROUTINE CONVECTION_LAYER(YDSURF, &
 ! Input quantities
  & YDMODEL,KDIM, LDSLPHY, STATE, TENDENCY_CML, TENDENCY_TMP, PAUX, PPERT, &
 ! Input/Output quantities
  & LLKEYS, PDIAG, AUXL, PERTL, FLUX, PSURF, GEMSL, &
 ! Output tendencies
  & TENDENCY_LOC)

!**** *CONVECTION_LAYER* - Layer routine calling full convection scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! LDSLPHY  : Key for SL phys
! state    : Derived variable for model state
! tendency_cml : D. V. for stored model tendencies from processes before 
! tendency_tmp : D. V. for model tendencies (entering convection) from processes before 
! PAUX     : Derived variables for general auxiliary quantities
! PPERT    : Derived variable for incoming perturbations etc... 

!     ==== Input/output ====
! LLKEYS       : Derived variable with keys
! PDIAG        : Derived variable for diagnostic quantities
! AUXL         : Derived variables for local auxiliary quantities
! PERTL        : Derived variables for local perturbartions
! FLUX         : Derived variable for fluxes
! PSURF        : Derived variables for general surface quantities
! GEMSL        : Derived variables for local GEMS quantities

!    ==== Output tendencies from convection ====
! tendency_loc :  Output process tendencies


!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!      Original : 2012-11-22  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     M. Leutbecher & S.-J. Lock (Jan 2016) Introduced SPP scheme (LSPP)
!     --------------

!-----------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
   & AUX_DIAG_TYPE, AUX_DIAG_LOCAL_TYPE, SURF_AND_MORE_TYPE, &
   & KEYS_LOCAL_TYPE, PERTURB_LOCAL_TYPE, FLUX_TYPE, GEMS_LOCAL_TYPE, &
   & PERTURB_TYPE
USE SPP_MOD  , ONLY : YSPP
USE YOMDYNCORE, ONLY : RPLDARE, RPLRG

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(MODEL) ,INTENT(INOUT) :: YDMODEL
TYPE (DIMENSION_TYPE)          , INTENT (IN)   :: KDIM
LOGICAL                        , INTENT (IN)   :: LDSLPHY
TYPE (STATE_TYPE)              , INTENT (IN)   :: STATE
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY_CML
TYPE (STATE_TYPE)              , INTENT (IN)   :: TENDENCY_TMP
TYPE (AUX_TYPE)                , INTENT (IN)   :: PAUX
TYPE (PERTURB_TYPE)            , INTENT (IN)   :: PPERT
TYPE (KEYS_LOCAL_TYPE)         , INTENT(INOUT) :: LLKEYS
TYPE (AUX_DIAG_TYPE)           , INTENT(INOUT) :: PDIAG
TYPE (AUX_DIAG_LOCAL_TYPE)     , INTENT(INOUT) :: AUXL
TYPE (PERTURB_LOCAL_TYPE)      , INTENT(INOUT) :: PERTL
TYPE (FLUX_TYPE)               , INTENT(INOUT) :: FLUX
TYPE (SURF_AND_MORE_TYPE)      , INTENT(INOUT) :: PSURF
TYPE (GEMS_LOCAL_TYPE)         , INTENT(INOUT) :: GEMSL
TYPE (STATE_TYPE)              , INTENT(INOUT) :: TENDENCY_LOC
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JARP, JL
REAL(KIND=JPRB)    :: ZGP2DSPP(KDIM%KLON, YSPP%N2D)          !SPP pattern

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "cucalln.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CONVECTION_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2)
ASSOCIATE(TSPHY=>YDPHY2%TSPHY, &
 & YSD_VN=>YDSURF%YSD_VN)
!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL CUCALLN

DO JARP=1, YSPP%N2D
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    ZGP2DSPP(JL,JARP)=PPERT%PGP2DSPP(JL,1,JARP)
  ENDDO
ENDDO

CALL CUCALLN &
  & (RPLDARE, RPLRG, YDMODEL%YRML_PHY_EC%YRTHF,YDMODEL%YRCST,YDMODEL%YRML_PHY_RAD%YRERAD,YDMODEL%YRML_PHY_SLIN,YDMODEL%YRML_PHY_EC,YDMODEL%YRML_GCONF%YGFL, &
  & KDIM%KIDIA  , KDIM%KFDIA , KDIM%KLON  , KDIM%KSMAX  , KDIM%KLEV, YSPP%N2D, &
  & LLKEYS%LLLAND, LDSLPHY,&
  & TSPHY,YDMODEL%YRML_PHY_G%YRVDF%RVDIFTS,&
  & STATE%T    , STATE%Q   , STATE%U   , STATE%V,   AUXL%ZLISUM,&
  & PAUX%PVERVEL, FLUX%PDIFTQ, FLUX%PDIFTS,  PAUX%PAPRS,&
  & PAUX%PRSF1  , PAUX%PRS1  , PAUX%PGEOM1, PAUX%PGEOMH, PAUX%PGAW,&
  & PERTL%ZCUCONVCA, ZGP2DSPP, &
  & TENDENCY_CML, TENDENCY_TMP, TENDENCY_LOC, & 
  & PSURF%GSD_VN%PACPR(:),&
  & AUXL%ITOPC  , AUXL%IBASC , PDIAG%ITYPE, &
  & PDIAG%ICBOT , PDIAG%ICTOP, AUXL%IBOTSC, LLKEYS%LLCUM , LLKEYS%LLSC,&
  & LLKEYS%LLSHCV,&
  & AUXL%ZLCRIT_AER,&
  & PDIAG%ZLU   , PDIAG%ZLUDE, PDIAG%ZLUDELI, PDIAG%ZSNDE, PDIAG%PMFU , PDIAG%PMFD,&
  & FLUX%PDIFCQ , FLUX%PDIFCS, FLUX%PFHPCL, FLUX%PFHPCN,&
  & FLUX%PFPLCL , FLUX%PFPLCN, PDIAG%ZLRAIN,PDIAG%PRSUD,& 
  & FLUX%PSTRCU, FLUX%PSTRCV, FLUX%PFCCQL, FLUX%PFCCQN,&
  & PDIAG%PMFUDE_RATE ,    PDIAG%PMFDDE_RATE ,   PDIAG%PCAPE,  PERTL%ZWMEAN,  PERTL%ZDISSCU,&
  & GEMSL%ITRAC  , GEMSL%ZCEN  , GEMSL%ZTENC, GEMSL%ZSCAV )


!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CONVECTION_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECTION_LAYER
