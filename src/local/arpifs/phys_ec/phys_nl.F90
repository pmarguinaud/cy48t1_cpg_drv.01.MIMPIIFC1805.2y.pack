SUBROUTINE PHYS_NL (YDTHF,YDCST,YDERAD,YDML_PHY_SLIN,YDML_PHY_EC,YDRIP,KLEV, KLON, PVDIFTS, &        
                 &  PT                , &
                 &  PQ                , &
                 &  PP                , &
                 &  PPH               , &
                 &  PTENT             , &
                 &  PTENQ             , & 
                 &  PDIFTS            , &
                 &  PDIFTQ            , &
                 &  PSTRTU            , &
                 &  PSTRTV            , &
                 &  PRAIN3D           , &
                 &  PSNOW3D           , &
                 &  PDQL              , &
                 &  PDQI              , &
                 &  PDCC              , &
                 &  PFLS              , &
                 &  PFLC)
   

!* Physics part of the observation operator: Mass-flux convection & large scale precip: 

!     AUTHOR.
!     -------
!*     Original (06/11/97) : J.-F. MAHFOUF

!     MODIFICATIONS.
!     --------------
!*     Modified (31/03/04) : New Interface w/ 1D-Var Obs. operator, clean-up. P.BAUER
!*     Modified (01/10/08) : LLRAIN1D name change to reflect usage. AGEER
!      K. Yessad (July 2014): Move some variables.
!------------------------------------------------------------------------

USE MODEL_PHYSICS_ECMWF_MOD      , ONLY : MODEL_PHYSICS_ECMWF_TYPE
USE MODEL_PHYSICS_SIMPLINEAR_MOD , ONLY : MODEL_PHYSICS_SIMPLINEAR_TYPE
USE YOERAD                       , ONLY : TERAD
USE PARKIND1                     , ONLY : JPIM, JPRB
USE YOMHOOK                      , ONLY : LHOOK,   DR_HOOK
USE YOMRIP                       , ONLY : TRIP
USE YOMCST                       , ONLY : TCST
USE YOETHF                       , ONLY : TTHF

!------------------------------------------------------------------------

IMPLICIT NONE

TYPE(TTHF)          ,INTENT(IN)                              :: YDTHF
TYPE(TCST)          ,INTENT(IN)                              :: YDCST
TYPE(TERAD)                        ,INTENT(INOUT)            :: YDERAD
TYPE(MODEL_PHYSICS_ECMWF_TYPE)     ,INTENT(INOUT)            :: YDML_PHY_EC
TYPE(MODEL_PHYSICS_SIMPLINEAR_TYPE),INTENT(INOUT)            :: YDML_PHY_SLIN
TYPE(TRIP)                         ,INTENT(INOUT)            :: YDRIP
INTEGER (KIND=JPIM), INTENT ( IN)                            :: KLEV, KLON
REAL(KIND=JPRB)    ,INTENT  (IN)                             :: PVDIFTS
REAL    (KIND=JPRB), INTENT ( IN),   DIMENSION (KLON,KLEV)   :: PP, PT, PQ
REAL    (KIND=JPRB), INTENT ( IN),   DIMENSION (KLON,KLEV+1) :: PPH
REAL    (KIND=JPRB), INTENT (INOUT), DIMENSION (KLON,KLEV)   :: PTENQ, PTENT
REAL    (KIND=JPRB), INTENT ( IN),   DIMENSION (KLON)        :: PDIFTS, PDIFTQ, PSTRTU, PSTRTV
REAL    (KIND=JPRB), INTENT (OUT),   DIMENSION (KLON,KLEV)   :: PRAIN3D, PSNOW3D 
REAL    (KIND=JPRB), INTENT (OUT),   DIMENSION (KLON,KLEV)   :: PDQL, PDQI, PDCC
REAL    (KIND=JPRB), INTENT (OUT),   DIMENSION (KLON)        :: PFLS, PFLC

!------------------------------------------------------------------------

REAL (KIND=JPRB), DIMENSION (KLON,KLEV)   :: ZTPROC, ZQPROC, ZUPROC, ZVPROC
REAL (KIND=JPRB), DIMENSION (KLON,KLEV)   :: ZU, ZV, ZW, ZTENU, ZTENV, ZQSAT, ZCAPE
REAL (KIND=JPRB), DIMENSION (KLON,KLEV)   :: ZTTENU, ZTTENV, ZTTENT, ZTTENQ
REAL (KIND=JPRB), DIMENSION (KLON,KLEV)   :: ZGEOM1, ZMFUDE_RATE, ZMFDDE_RATE
REAL (KIND=JPRB), DIMENSION (KLON,KLEV)   :: ZLU, ZLUDE, ZMFU, ZMFD 
REAL (KIND=JPRB), DIMENSION (KLON,KLEV)   :: ZCOVPTOT
REAL (KIND=JPRB), DIMENSION (KLON,KLEV+1) :: ZGEOMH, ZFPLCL, ZFPLCN, ZFHPCL, ZFHPCN 
REAL (KIND=JPRB), DIMENSION (KLON,KLEV+1) :: ZFHPSL, ZFHPSN, ZFPLSL, ZFPLSN
REAL (KIND=JPRB), DIMENSION (KLON,KLEV+1) :: ZDIFTQ, ZDIFTS, ZDIFCQ, ZDIFCS
REAL (KIND=JPRB), DIMENSION (KLON,KLEV+1) :: ZSTRCU, ZSTRCV, ZFCCQL, ZFCCQN

REAL (KIND=JPRB), DIMENSION (KLON)        :: ZARPRC, ZGAW
! Allocated for subroutine calls
REAL (KIND=JPRB), DIMENSION (1,1,1) :: ZCEN5, ZTENC5
REAL (KIND=JPRB), DIMENSION (1)     :: ZSCAV5

REAL (KIND=JPRB) :: ZTSPHY

INTEGER (KIND=JPIM), DIMENSION (KLON) :: ITOPC, IBASC, JTOPC, JBASEC, JBOTSC, JTYPE
INTEGER (KIND=JPIM) :: ISTEP, ISTART, IIDIA, IFDIA, ITDIA
INTEGER (KIND=JPIM) :: JL, JK, ITRAC
     
LOGICAL, DIMENSION (KLON) :: LLLAND, LLCUM, LLSC 
LOGICAL :: LLSLPHY, LLRAIN1D
      
REAL (KIND=JPRB) :: ZHOOK_HANDLE

!------------------------------------------------------------------------

#include "satur.intfb.h"
#include "cucalln2.intfb.h"
#include "cloudst.intfb.h"

!------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PHYS_NL',0,ZHOOK_HANDLE)
ASSOCIATE(TSTEP=>YDRIP%TSTEP)
!------------------------------------------------------------------------

IIDIA  =    1
IFDIA  = KLON
ITDIA  =    1
ISTEP  =    1
ISTART =    1
ITOPC  = KLEV
IBASC  =    1

LLRAIN1D = .TRUE.
LLSLPHY = .FALSE.
LLLAND  = .FALSE. ! Force the estimation of precipitation over oceans

!   TIME STEP
ZTSPHY = TSTEP

ZMFUDE_RATE (:,:) = 0.0_JPRB
ZMFDDE_RATE (:,:) = 0.0_JPRB
         
ZU    (:,:) = 0.0_JPRB
ZV    (:,:) = 0.0_JPRB 
ZW    (:,:) = 0.0_JPRB 
ZTENU (:,:) = 0.0_JPRB
ZTENV (:,:) = 0.0_JPRB

!* geopotential
ZGEOMH (:,KLEV+1) = 0.0_JPRB
DO JK = KLEV,1,-1
   ZGEOMH (:,JK) = ZGEOMH (:,JK+1) + &
  & 287.0_JPRB*PT(:,JK)*(1.0_JPRB+0.608_JPRB*PQ(:,JK))/PP(:,JK)*(PPH(:,JK+1)-PPH(:,JK))
ENDDO               
DO JK=1,KLEV
   ZGEOM1 (:,JK) = 0.5_JPRB * (ZGEOMH (:,JK+1) + ZGEOMH (:,JK))
ENDDO

DO JL = 1, KLON
   ZDIFTS (JL,KLEV+1) = PDIFTS (JL)
   ZDIFTQ (JL,KLEV+1) = PDIFTQ (JL)
ENDDO

!* Specific humidity at saturation
CALL SATUR (YDTHF, YDCST, IIDIA, IFDIA, KLON, ITDIA, KLEV, YDML_PHY_SLIN%YREPHLI%LPHYLIN, PP, PT, ZQSAT, 2)

!* Simplified convection scheme     
!!! this is set to a constant appropriate for T799, needs fixing !!!
ZGAW(:)=1.2E-6_JPRB
!!! this is not correct, needs fixing !!!
ITRAC=0
ZTTENT(:,:)=PTENT(:,:)
ZTTENQ(:,:)=PTENQ(:,:)
ZTTENU(:,:)=ZTENU(:,:)
ZTTENV(:,:)=ZTENV(:,:)
CALL CUCALLN2 (YDTHF,YDCST,YDERAD,YDML_PHY_SLIN,YDML_PHY_EC   , &
            &  IIDIA  , IFDIA , KLON  , KLEV                  , &
            &  LLLAND, LLSLPHY, LLRAIN1D                      , &
            &  ZTSPHY  ,  PVDIFTS                             , &
            &  PT     , PQ    , ZU    , ZV                    , &
            &  ZW     , ZDIFTQ, ZDIFTS, PPH                   , &
            &  PP     , PPH   , ZGEOM1, ZGEOMH                ,  ZGAW, &
            &  ZTPROC , ZTTENT, PTENT , ZQPROC, ZTTENQ, PTENQ , &
            &  ZUPROC , ZTTENU, ZTENU , ZVPROC, ZTTENV, ZTENV , &
            &  ZARPRC                                         , &
            &  ITOPC  , IBASC , JTYPE                         , &
            &  JBASEC , JTOPC , JBOTSC, LLCUM , LLSC          , & 
            &  ZLU    , ZLUDE , ZMFU  , ZMFD                  , &
            &  ZDIFCQ , ZDIFCS, ZFHPCL, ZFHPCN                , &
            &  ZFPLCL , ZFPLCN, ZSTRCU, ZSTRCV, ZFCCQL, ZFCCQN, &
            &  ZMFUDE_RATE , ZMFDDE_RATE , ZCAPE              , &
            &  ITRAC   , ZCEN5  , ZTENC5 , ZSCAV5)
PTENT(:,:)=PTENT(:,:)+ZTPROC(:,:)
PTENQ(:,:)=PTENQ(:,:)+ZQPROC(:,:)
ZTENU(:,:)=ZTENU(:,:)+ZUPROC(:,:)
ZTENV(:,:)=ZTENV(:,:)+ZVPROC(:,:)

!* Simplified clouds scheme
CALL CLOUDST (YDML_PHY_SLIN%YREPHLI,YDML_PHY_SLIN%YRPHNC,YDML_PHY_EC%YRECLD,YDML_PHY_EC%YRECLDP, &
            & IIDIA,  IFDIA,  KLON,    ITDIA, KLEV, &
            & LLRAIN1D, ZTSPHY                    , &
            & PPH    , PP  ,    PQ,    ZQSAT, PT  , &
            & ZLUDE  , ZLU                        , &
            & ZTPROC , PTENT  , ZQPROC  , PTENQ   , &
            & PDCC   , PDQI,   PDQL               , &
            & ZFPLSL , ZFPLSN, ZFHPSL, ZFHPSN     , &
            & ZCOVPTOT)
PTENT(:,:)=PTENT(:,:)+ZTPROC(:,:)
PTENQ(:,:)=PTENQ(:,:)+ZQPROC(:,:)

!* 2D/3D and strat./conv. precipitation flux 
DO JL = 1, KLON
!    DO JK = 1, KLEV + 1
!       ZFPLCL (JL,JK) = MAX (0.0_JPRB, ZFPLCL (JL,JK))
!       ZFPLCN (JL,JK) = MAX (0.0_JPRB, ZFPLCN (JL,JK))
!       ZFPLSL (JL,JK) = MAX (0.0_JPRB, ZFPLSL (JL,JK))
!       ZFPLSN (JL,JK) = MAX (0.0_JPRB, ZFPLSN (JL,JK))
!    ENDDO
    DO JK = 1, KLEV
       PRAIN3D (JL,JK) = (ZFPLSL (JL,JK  ) + ZFPLCL (JL,JK  ) &
                     & +  ZFPLSL (JL,JK+1) + ZFPLCL (JL,JK+1))*0.5_JPRB
       PSNOW3D (JL,JK) = (ZFPLSN (JL,JK  ) + ZFPLCN (JL,JK  ) &
                     & +  ZFPLSN (JL,JK+1) + ZFPLCN (JL,JK+1))*0.5_JPRB
    ENDDO
    PFLC (JL) = ZFPLCL (JL,KLEV+1) + ZFPLCN (JL,KLEV+1)
    PFLS (JL) = ZFPLSL (JL,KLEV+1) + ZFPLSN (JL,KLEV+1)
ENDDO

!------------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('PHYS_NL',1,ZHOOK_HANDLE)
END SUBROUTINE PHYS_NL      
