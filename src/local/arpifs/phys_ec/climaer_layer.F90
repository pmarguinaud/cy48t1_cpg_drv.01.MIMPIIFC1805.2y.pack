SUBROUTINE CLIMAER_LAYER(YDSURF,YDML_PHY_RAD,YDPHY2,KDIM,LDPHYLIN,PAUX,STATE,PSURF,GEMSL)

!**** *CLIMAER_LAYER* - Layer routine calling aerosol climatology

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDIM     : Derived variable for dimensions
! PAUX     : Derived variable for general auxiliary quantities
! state    : Derived variable for model state

!     ==== Input/output ====
! PSURF    : Derived variable for general surface quantities
! GEMSL    : Derived variables for local GEMS quantities


!        IMPLICIT ARGUMENTS :   NONE
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
!      Original : 2012-12-04  F. VANA (c) ECMWF

!     MODIFICATIONS.
!     --------------
!     R Hogan  3 Oct 2016: Visibility for MACC aerosol climatology
!     R Hogan 28 Mar 2017: RAD_CONFIG now in YRADIATION object
!-----------------------------------------------------------------------

USE MODEL_PHYSICS_RADIATION_MOD , ONLY : MODEL_PHYSICS_RADIATION_TYPE
USE YOECLDP                     , ONLY : NCLDQR, NCLDQS, NCLDQI, NCLDQL
USE SURFACE_FIELDS_MIX          , ONLY : TSURF
USE PARKIND1                    , ONLY : JPIM , JPRB
USE YOMHOOK                     , ONLY : LHOOK, DR_HOOK
USE YOMPHYDER                   , ONLY : DIMENSION_TYPE, STATE_TYPE, AUX_TYPE, &
 &                                       SURF_AND_MORE_TYPE, GEMS_LOCAL_TYPE
USE YOMPHY2                     , ONLY : TPHY2
USE YOE_AERODIAG                , ONLY : JPAERO_WVL_AOD
USE YOMCST                      , ONLY : RG, YRCST
USE YOETHF                      , ONLY : YRTHF
USE RADIATION_AEROSOL_OPTICS    , ONLY : AEROSOL_SW_EXTINCTION

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF)                       , INTENT (INOUT) :: YDSURF
TYPE(MODEL_PHYSICS_RADIATION_TYPE), INTENT (INOUT) :: YDML_PHY_RAD
TYPE(TPHY2)                       , INTENT (INOUT) :: YDPHY2
LOGICAL                           , INTENT (IN)    :: LDPHYLIN
TYPE (DIMENSION_TYPE)             , INTENT (IN)    :: KDIM
TYPE (AUX_TYPE)                   , INTENT (IN)    :: PAUX
TYPE (STATE_TYPE)                 , INTENT (IN)    :: STATE
TYPE (SURF_AND_MORE_TYPE)         , INTENT (INOUT) :: PSURF
TYPE (GEMS_LOCAL_TYPE)            , INTENT (INOUT) :: GEMSL
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JL, JK, JAER
REAL(KIND=JPRB) :: ZDUM(KDIM%KLON,KDIM%KLEV)
REAL(KIND=JPRB) :: ZMACCAER(KDIM%KLON,KDIM%KLEV,YDML_PHY_RAD%YRERAD%NMCVAR) ! MACC aerosol layer mass (kg/m2)
REAL(KIND=JPRB) :: ZINV_LAYER_MASS(KDIM%KLON) ! Inverse of mass of air in a layer (m2/kg)
REAL(KIND=JPRB) :: ZMACCAER_MMR(KDIM%KLON,YDML_PHY_RAD%YRADIATION%RAD_CONFIG%N_AEROSOL_TYPES) ! MACC aerosol mixing ratio (g/kg)

REAL(KIND=JPRB) :: ZTH(KDIM%KLON,KDIM%KLEV+1) ! Half level temperature (K)
REAL(KIND=JPRB) :: ZQAER(KDIM%KLON,6,KDIM%KLEV), ZTAER(KDIM%KLON,6)
REAL(KIND=JPRB) :: ZQSAT(KDIM%KLON,KDIM%KLEV) ! Saturation specific humidity (kg/kg)
REAL(KIND=JPRB) :: ZRH(KDIM%KLON)             ! Relative humidity in lowest model level

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

#include "radaca.intfb.h"
#include "claervis.intfb.h"
#include "satur.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CLIMAER_LAYER',0,ZHOOK_HANDLE)
ASSOCIATE(LAERCLIM=>YDML_PHY_RAD%YRERAD%LAERCLIM, LECSRAD=>YDML_PHY_RAD%YRERAD%LECSRAD, &
 & NCSRADF=>YDML_PHY_RAD%YRERAD%NCSRADF, YDRADIATION=>YDML_PHY_RAD%YRADIATION,&
 & YSD_VD=>YDSURF%YSD_VD, YSP_RR=>YDSURF%YSP_RR, &
 & TSPHY=>YDPHY2%TSPHY)
!     ------------------------------------------------------------------

!*         0.     Preliminary computation

DO JK=2,KDIM%KLEV
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    ZTH(JL,JK)=(STATE%T(JL,JK-1)*PAUX%PAPRSF(JL,JK-1)&
      & *(PAUX%PAPRSF(JL,JK)-PAUX%PAPRS(JL,JK))&
      & +STATE%T(JL,JK)*PAUX%PAPRSF(JL,JK)*(PAUX%PAPRS(JL,JK)-PAUX%PAPRSF(JL,JK-1)))&
      & *(1.0_JPRB/(PAUX%PAPRS(JL,JK)*(PAUX%PAPRSF(JL,JK)-PAUX%PAPRSF(JL,JK-1))))  
  ENDDO
ENDDO
DO JL=KDIM%KIDIA,KDIM%KFDIA
  ZTH(JL,1)=STATE%T(JL,1)-PAUX%PAPRSF(JL,1)*(STATE%T(JL,1)-ZTH(JL,2))&
    & /(PAUX%PAPRSF(JL,1)-PAUX%PAPRS(JL,2))  
  ZTH(JL,KDIM%KLEV+1)=PSURF%GSP_RR%PT_T9(JL)
ENDDO

!*         1.     UNROLL THE DERIVED STRUCTURES AND CALL RADACA & CLAERVIS

! RADACA does too much for its use here: we only really need aerosol
! from one scheme at the surface yet it computes aerosols from all
! schemes at all heights. Something to improve in future.
CALL  RADACA ( YDML_PHY_RAD%YREAERD,YDML_PHY_RAD%YRERAD, KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV, &
   & PAUX%PAPRS, PAUX%PGELAM, PAUX%PGEMU, &
   & PAUX%PCLON, PAUX%PSLON , ZTH  , &
   & ZQAER, ZMACCAER, ZDUM , GEMSL%ZCLAERS , &
   & PSURF%GSD_VD%PAERO_WVL_DIAG(1,JPAERO_WVL_AOD)%P(:), PSURF%GSD_VD%PODSS(:), PSURF%GSD_VD%PODDU(:),&
   & PSURF%GSD_VD%PODOM(:), PSURF%GSD_VD%PODBC(:), PSURF%GSD_VD%PODSU(:))

IF ((YDML_PHY_RAD%YREAERATM%LAERCCN .OR. YDML_PHY_RAD%YREAERATM%LAERRRTM .OR. YDML_PHY_RAD%YRERAD%NAERMACC == 1) &
     &  .AND. .NOT. YDML_PHY_RAD%YRERAD%LUSEPRE2017RAD) THEN
  ! We are using the MACC/CAMS aerosol climatology AND the newer
  ! radiation scheme: overwrite the extinction coefficient in
  ! GEMSL%ZCLAERS with the value for the MACC/CAMS aerosol climatology
  CALL SATUR(YRTHF, YRCST, KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV, KDIM%KLEV, LDPHYLIN, &
       &  PAUX%PAPRSF, STATE%T, ZQSAT, 2)
  DO JL=KDIM%KIDIA,KDIM%KFDIA
    ZRH(JL) = STATE%Q(JL,KDIM%KLEV)/ZQSAT(JL,KDIM%KLEV)
    ZINV_LAYER_MASS(JL) = RG / (PAUX%PAPRSF(JL,KDIM%KLEV)-PAUX%PAPRSF(JL,KDIM%KLEV-1))
  ENDDO

  ! Convert aerosol layer mass in kg/m2 to mass mixing ratio in g/kg
  !
  ! NOTE: This whole mechanism appears to be predicated on the assumption that the set of
  !   aerosol types used here (unconditionally from the climatology?) is the same as that
  !   used in the radiation scheme (from prognostics if LAERRRTM=true).
  !
  !   If radiation is using prognostic aerosol but visibility is using the climatology,
  !   don't we then really need separate instances of RAD_CONFIG?
  !
  !   I've patched it up to not actually crash at least when N_AEROSOL_TYPES /= NMCVAR, but
  !   this will only behave correctly if the types in the climatology match the first
  !   NMCVAR types in the prognostic configuration...
  !
  !   Zak 2018-11-12
  DO JAER = 1,MIN(YDML_PHY_RAD%YRERAD%NMCVAR,YDRADIATION%RAD_CONFIG%N_AEROSOL_TYPES)
    DO JL = KDIM%KIDIA,KDIM%KFDIA
      ZMACCAER_MMR(JL,JAER) = 1000.0_JPRB*ZMACCAER(JL,KDIM%KLEV,JAER)*ZINV_LAYER_MASS(JL)
    ENDDO
  ENDDO
  IF (YDRADIATION%RAD_CONFIG%N_AEROSOL_TYPES > YDML_PHY_RAD%YRERAD%NMCVAR) THEN
    ZMACCAER_MMR(KDIM%KIDIA:KDIM%KFDIA,YDML_PHY_RAD%YRERAD%NMCVAR+1:YDRADIATION%RAD_CONFIG%N_AEROSOL_TYPES) = 0._JPRB
  ENDIF


  ! Use new radiation scheme to compute extinction coefficient in the
  ! 10th shortwave band, which includes 550 nm; by providing mass
  ! mixing ratio in g/kg we get output in km-1, which is what CLAERVIS
  ! expects. This non-SI stuff wastes everyone's time and should be
  ! fixed!
  CALL AEROSOL_SW_EXTINCTION(KDIM%KLON, KDIM%KIDIA, KDIM%KFDIA, &
       &  YDRADIATION%RAD_CONFIG, 10, ZMACCAER_MMR, ZRH, GEMSL%ZCLAERS)

ENDIF

CALL  CLAERVIS ( KDIM%KIDIA, KDIM%KFDIA, KDIM%KLON, KDIM%KLEV, &
     & GEMSL%ZCLAERS, STATE%A, STATE%CLD(:,:,NCLDQL), STATE%CLD(:,:,NCLDQI),&
     & STATE%CLD(:,:,NCLDQR),  STATE%CLD(:,:,NCLDQS), PAUX%PAPRSF, STATE%T, &
     & GEMSL%ZVISICL )     

!-- for diagnostics only: store the 2D distribution of climatological aerosols in EXTRA fields    
IF (LECSRAD .AND. NCSRADF == 1 .AND. LAERCLIM) THEN
  DO JAER=1,6
    ZTAER(KDIM%KIDIA:KDIM%KFDIA,JAER)=0._JPRB        
    DO JK=1,KDIM%KLEV
      DO JL=KDIM%KIDIA,KDIM%KFDIA    
        ZTAER(JL,JAER)=ZTAER(JL,JAER)+ZQAER(JL,JAER,JK)
      ENDDO
    ENDDO
    DO JL=KDIM%KIDIA,KDIM%KFDIA
      PSURF%GSD_XA%PGROUP(JL,JAER+8,14)=PSURF%GSD_XA%PGROUP(JL,JAER+8,14)+ZTAER(JL,JAER)*TSPHY
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CLIMAER_LAYER',1,ZHOOK_HANDLE)
END SUBROUTINE CLIMAER_LAYER
