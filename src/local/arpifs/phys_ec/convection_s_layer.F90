SUBROUTINE CONVECTION_S_LAYER(YDTHF, YDCST, YDSURF, &
 ! Input quantities
  & YDERAD,YDML_PHY_SLIN,YDML_PHY_EC,YDPHY2,KDIM, LDSLPHY, LDRAIN1D, STATE, TENDENCY_CML, TENDENCY_TMP, PAUX, PVDIFTS,&
 ! Input/Output quantities
  & LLKEYS, PDIAG, AUXL, FLUX, PSURF, GEMSL, &
 ! Output tendencies
  & TENDENCY_LOC)

!**** *CONVECTION_LAYER* - Layer routine calling simplified convection scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! LDSLPHY  : Key for SL phys
! LDRAIN1D : Key for 1D rain assim
! state    : Derived variable for  model state
! tendency_cml : D. V. for model resulting tendencies
! tendency_tmp : D. V. for model tendencies entering convection
! PAUX     : Derived variables for general auxiliary quantities

!     ==== Input/output ====
! LLKEYS       : Derived variable with keys
! PDIAG        : Derived variable for diagnostic quantities
! AUXL         : Derived variables for local auxiliary quantities
! FLUX         : Derived variable for fluxes
! PSURF        : Derived variables for general surface quantities
! GEMSL        : Derived variables for local GEMS quantities

!    ==== Output tendencies from convection ====
! tendency_loc :  Derived variables with process tendencies


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
!     --------------
!     F. Vana  Oct-2013  Bugfix

!-----------------------------------------------------------------------

USE MODEL_PHYSICS_ECMWF_MOD      , ONLY : MODEL_PHYSICS_ECMWF_TYPE
USE MODEL_PHYSICS_SIMPLINEAR_MOD , ONLY : MODEL_PHYSICS_SIMPLINEAR_TYPE
USE YOERAD                       , ONLY : TERAD
USE SURFACE_FIELDS_MIX           , ONLY : TSURF
USE YOMCST                       , ONLY : TCST
USE YOETHF                       , ONLY : TTHF
USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
   &                  AUX_DIAG_TYPE, AUX_DIAG_LOCAL_TYPE, SURF_AND_MORE_TYPE, &
   &                  KEYS_LOCAL_TYPE, FLUX_TYPE, GEMS_LOCAL_TYPE
USE YOMPHY2  , ONLY : TPHY2

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TTHF)                    ,INTENT(IN)    :: YDTHF
TYPE(TCST)                    ,INTENT(IN)    :: YDCST
TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(TERAD), INTENT(INOUT) :: YDERAD
TYPE(MODEL_PHYSICS_ECMWF_TYPE),INTENT(INOUT)     :: YDML_PHY_EC
TYPE(MODEL_PHYSICS_SIMPLINEAR_TYPE),INTENT(INOUT):: YDML_PHY_SLIN
TYPE(TPHY2)               , INTENT(INOUT) :: YDPHY2
TYPE (DIMENSION_TYPE)     , INTENT (IN)   :: KDIM
LOGICAL                   , INTENT (IN)   :: LDSLPHY
LOGICAL                   , INTENT (IN)   :: LDRAIN1D
TYPE (STATE_TYPE)         , INTENT (IN)   :: STATE
TYPE (STATE_TYPE)         , INTENT (IN)   :: TENDENCY_CML
TYPE (STATE_TYPE)         , INTENT (IN)   :: TENDENCY_TMP
TYPE (AUX_TYPE)           , INTENT (IN)   :: PAUX
REAL(KIND=JPRB)           , INTENT (IN)   :: PVDIFTS
TYPE (KEYS_LOCAL_TYPE)    , INTENT(INOUT) :: LLKEYS
TYPE (AUX_DIAG_TYPE)      , INTENT(INOUT) :: PDIAG
TYPE (AUX_DIAG_LOCAL_TYPE), INTENT(INOUT) :: AUXL
TYPE (FLUX_TYPE)          , INTENT(INOUT) :: FLUX
TYPE (SURF_AND_MORE_TYPE) , INTENT(INOUT) :: PSURF
TYPE (GEMS_LOCAL_TYPE)    , INTENT(INOUT) :: GEMSL
TYPE (STATE_TYPE)         , INTENT(INOUT) :: TENDENCY_LOC
!-----------------------------------------------------------------------
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "cucalln2.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CONVECTION_S_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(TSPHY=>YDPHY2%TSPHY, &
 & YSD_VN=>YDSURF%YSD_VN)
!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL CUCALLN2

CALL CUCALLN2 &
  & (YDTHF, YDCST, YDERAD,YDML_PHY_SLIN,YDML_PHY_EC, &
  & KDIM%KIDIA  , KDIM%KFDIA , KDIM%KLON  , KDIM%KLEV,&
  & LLKEYS%LLLAND, LDSLPHY, LDRAIN1D, &
  & TSPHY,PVDIFTS,&
  & STATE%T     , STATE%Q  , STATE%U    , STATE%V,&
  & PAUX%PVERVEL, FLUX%PDIFTQ, FLUX%PDIFTS, PAUX%PAPRS,&
  & PAUX%PRSF1  , PAUX%PRS1  , PAUX%PGEOM1, PAUX%PGEOMH, PAUX%PGAW,&
  & TENDENCY_LOC%T, TENDENCY_CML%T, TENDENCY_TMP%T,&
  & TENDENCY_LOC%Q, TENDENCY_CML%Q, TENDENCY_TMP%Q,&
  & TENDENCY_LOC%U, TENDENCY_CML%U, TENDENCY_TMP%U,&
  & TENDENCY_LOC%V ,TENDENCY_CML%V, TENDENCY_TMP%V,&
  & PSURF%GSD_VN%PACPR(:),&
  & AUXL%ITOPC  , AUXL%IBASC , PDIAG%ITYPE,&
  & PDIAG%ICBOT , PDIAG%ICTOP, AUXL%IBOTSC, LLKEYS%LLCUM , LLKEYS%LLSC,&
  & PDIAG%ZLU   , PDIAG%ZLUDE, PDIAG%PMFU , PDIAG%PMFD,&
  & FLUX%PDIFCQ , FLUX%PDIFCS, FLUX%PFHPCL, FLUX%PFHPCN,&
  & FLUX%PFPLCL , FLUX%PFPLCN, FLUX%PSTRCU, FLUX%PSTRCV, FLUX%PFCCQL, FLUX%PFCCQN,&
  & PDIAG%PMFUDE_RATE ,    PDIAG%PMFDDE_RATE ,   PDIAG%PCAPE,&
  & GEMSL%ITRAC  , GEMSL%ZCEN  , GEMSL%ZTENC,  GEMSL%ZSCAV )


!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CONVECTION_S_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECTION_S_LAYER
