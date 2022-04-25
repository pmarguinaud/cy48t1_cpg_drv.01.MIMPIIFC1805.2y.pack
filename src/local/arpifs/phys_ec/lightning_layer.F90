SUBROUTINE LIGHTNING_LAYER(&
 ! Input quantities
 &  YDTHF, YDCST, YDCHEM,YDEPHY,YGFL,KDIM, PAUX, LLKEYS, STATE, PDIAG, FLUX,&
 ! Output quantity
 &  VNOEMI)

!**** *LIGHTNING_LAYER* - Layer routine calling lightning scheme

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variables for general auxiliary quantities
! LLKEYS   : Derived variable with keys
! state    : Derived variable for model state
! tendency : D. V. for model tendencies (entering convection) from processes before 
! PDIAG    : Derived variable for diagnostic quantities
! FLUX     : Derived variable for fluxes

!    ==== Output tendencies from convection ====
! PNOEMI   :  3D NO Emissions in            kg/m*2/s


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
!     P. Lopez 08/10/2015 : Added separate routine for lightning NOx computations.

!-----------------------------------------------------------------------

USE YOEPHY   , ONLY : TEPHY
USE YOM_YGFL , ONLY : TYPE_GFLD
USE PARKIND1  ,ONLY : JPIM ,   JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMPHYDER ,ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
                    & AUX_DIAG_TYPE, KEYS_LOCAL_TYPE, FLUX_TYPE
USE YOMCHEM  , ONLY : TCHEM
USE VARIABLE_MODULE, ONLY: VARIABLE_3D
USE YOMCST    , ONLY : TCST
USE YOETHF    , ONLY : TTHF  
USE YOMDYNCORE, ONLY : RPLDARE, RPLRG

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TTHF)                      ,INTENT(IN)     :: YDTHF
TYPE(TCST)                      ,INTENT(IN)     :: YDCST
TYPE(TCHEM)                     ,INTENT(INOUT)  :: YDCHEM
TYPE(TEPHY)                     ,INTENT(INOUT)  :: YDEPHY
TYPE(TYPE_GFLD)                 ,INTENT(INOUT)  :: YGFL
TYPE (DIMENSION_TYPE)          , INTENT (IN)    :: KDIM
TYPE (AUX_TYPE)                , INTENT (IN)    :: PAUX
TYPE (KEYS_LOCAL_TYPE)         , INTENT (IN)    :: LLKEYS
TYPE (STATE_TYPE)              , INTENT (IN)    :: STATE
TYPE (AUX_DIAG_TYPE)           , INTENT (INOUT) :: PDIAG
TYPE (FLUX_TYPE)               , INTENT (IN)    :: FLUX
TYPE (VARIABLE_3D)             , INTENT (INOUT) :: VNOEMI

!-----------------------------------------------------------------------

REAL(KIND=JPRB) :: ZLIGH_EMI(KDIM%KLON), ZWMFU(KDIM%KLON)
REAL(KIND=JPRB) :: ZCTOPH(KDIM%KLON), ZPRECMX(KDIM%KLON), ZICETOT(KDIM%KLON), ZCDEPTH(KDIM%KLON)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL :: LLLINOX(KDIM%KLON)

!-----------------------------------------------------------------------

#include "culight.intfb.h"
#include "culinox.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LIGHTNING_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(LCHEM_LIGHT=>YDCHEM%LCHEM_LIGHT)
!     ------------------------------------------------------------------

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL CULIGHT

CALL CULIGHT (RPLDARE, RPLRG, YDTHF, YDCST,  YDEPHY,YGFL, KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV, &
  &    PAUX%PGAW,       PAUX%PGELAT,&
  &    PAUX%PRSF1,      PAUX%PRS1,      PAUX%PAPHI,   PAUX%PAPHIF,  LLKEYS%LLLAND,&
  &    STATE%T,         PDIAG%ZLU,      PDIAG%PMFU,   PDIAG%PCAPE,&
  &    FLUX%PFPLCL,     FLUX%PFPLCN,&
  &    LLKEYS%LLCUM,    PDIAG%ICBOT,    PDIAG%ICTOP,&
  ! Outputs
  &    LLLINOX,   PDIAG%PLIGH_TOT, PDIAG%PLIGH_CTG, ZCTOPH,&
  &    ZPRECMX,   ZICETOT,   ZCDEPTH,  ZWMFU)


!*         2.     CALL CULINOX (LIGHTNING NOx COMPUTATIONS)

IF (LCHEM_LIGHT) THEN 

  CALL CULINOX ( YDCHEM, KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV, &
    &    PAUX%PGELAT,     PAUX%PRS1,    PAUX%PAPHI,   PAUX%PAPHIF, &
    &    LLKEYS%LLLAND,   LLLINOX, &
    &    STATE%T,         PDIAG%ICTOP,&
    &    PDIAG%PLIGH_TOT, PDIAG%PLIGH_CTG, &
    ! Outputs
    &    VNOEMI%P, ZLIGH_EMI)

ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LIGHTNING_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE LIGHTNING_LAYER
