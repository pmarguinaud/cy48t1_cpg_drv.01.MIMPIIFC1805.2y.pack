SUBROUTINE BRI2ACCONV(YDML_PHY_MF,YDEGEO, &
 & KIDIA,KFDIA,&
 & KLON,KLEV,&
 & PGM,PAPRSF,PZZF,&
 & PT,PRV,PRC,PRI,PRHOREF,&
 & PU,PV,PW,&
 & PMF,PTTEN,PRVTEN,PRCTEN,PRITEN,&
 & PPRTEN,PPRSTEN)  

! Purpose:
! --------
!*** BRI2ACCON - Bridge to MNH convection call (deep and shallow)
!***           - This routine allow to go towards MNH world.
!***           - Interface routine on ARP/ALD side.
 
! Interface:
! ----------
 
!-  Dimensions
 
!KIDIA,KFDIA : Start/end of horizontal loop
!KLON : Horizontal dimension (NPROMA in CPG)
!KLEV : Vertical dimension
 
!- 2D (1:KLEV)
 
!PAPRSF : Pressure on full levels 
!PZZF : Altitude of full levels
!PT : Temperature
!PRV : Water vapor (rv)
!PRC : Cloud contain (rc)
!PRI : Ice (ri)
!PRHOREF : rho
!PU : X-Component of wind
!PV : Y-component of wind
!PW : Vertical velocity
 
! Output arguments
! ----------------
 
!-2D (1:KLEV)
 
!PMF   : Convective Mass flux for subgrid condensation
!PTTEN : Convective T tendency (K/s)
!PRVTEN : Convective r_v tendency (1/s)
!PRCTEN : Convective r_c tendency (1/s)
!PRITEN : Convective r_i tendency (1/s)

!-1D
!PPRTEN : Total surf precipitation tendency (m/s)
!PPRSTEN : Solid surf precipitation tendency (m/s)

!  Method
!  ------
!  None
 
!  Externals
!  ---------
!  AC_CONV_MNH (interface array on MNH side)
 
! Author
! -----
! G. Hello : 04-02-01 
 
! Modifications
! -------------
!  04-05-04 T.Kovacic: added PPRTEN and  PPRSTEN to argument list 
!  M.Hamrud      01-Oct-2003 CY28 Cleaning
!  T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
! End modifications
!----------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE YOMCT0  , ONLY : LELAM
USE YEMGEO  , ONLY : TEGEO

IMPLICIT NONE

!Input variables

TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
TYPE(TEGEO), INTENT(IN)          :: YDEGEO
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGM(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZZF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRI(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHOREF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PW(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTTEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRVTEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRCTEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRITEN(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPRTEN(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PPRSTEN(KLON) 
!Local variables

REAL(KIND=JPRB) :: ZDTCONV
REAL(KIND=JPRB) :: ZDXDY(KLON) ! grid area (m2)
REAL(KIND=JPRB) :: ZCAPE(KLON)
REAL(KIND=JPRB) :: ZT(KLON,KLEV), ZPRLFLX(KLON,KLEV),&
 & ZPRSFLX(KLON,KLEV), ZUMF(KLON,KLEV), ZDMF(KLON,KLEV)  

INTEGER(KIND=JPIM) :: ICLTOP(KLON), ICLBAS(KLON)

INTEGER(KIND=JPIM) :: JLON, JLEV
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "aro_conv_mnh.h"


IF (LHOOK) CALL DR_HOOK('BRI2ACCONV',0,ZHOOK_HANDLE)
ASSOCIATE(TEQK=>YDML_PHY_MF%YRPHY0%TEQK, &
 & TSPHY=>YDML_PHY_MF%YRPHY2%TSPHY, &
 & LDOWN=>YDML_PHY_MF%YRCVMNH%LDOWN, NIICE=>YDML_PHY_MF%YRCVMNH%NIICE, LDEEP=>YDML_PHY_MF%YRCVMNH%LDEEP, &
 & NSETENS=>YDML_PHY_MF%YRCVMNH%NSETENS, OTADJD=>YDML_PHY_MF%YRCVMNH%OTADJD, LSETTADJ=>YDML_PHY_MF%YRCVMNH%LSETTADJ, &
 & LREFRESH_ALL=>YDML_PHY_MF%YRCVMNH%LREFRESH_ALL, LDIAGCONV=>YDML_PHY_MF%YRCVMNH%LDIAGCONV, &
 & LSHALLOW=>YDML_PHY_MF%YRCVMNH%LSHALLOW, OTADJS=>YDML_PHY_MF%YRCVMNH%OTADJS, &
 & EDELY=>YDEGEO%EDELY, EDELX=>YDEGEO%EDELX)
! Convective time step
! --------------------

ZDTCONV=TSPHY

! Grid size
! ---------

IF(LELAM) THEN
  DO JLON=KIDIA,KFDIA
    ZDXDY(JLON)=EDELX*EDELY/PGM(JLON)**2
  ENDDO
ELSE
  DO JLON=KIDIA,KFDIA
    ZDXDY(JLON)=1.0_JPRB/(TEQK*PGM(JLON))**2
  ENDDO
ENDIF

! Temperature on reverse order for MNH
! ------------------------------------

DO JLON=KIDIA,KFDIA
  DO JLEV=1,KLEV
    ZT(JLON,JLEV)    = PT(JLON,JLEV)
  ENDDO
ENDDO

! Cross towards MNH world for KFB scheme
! --------------------------------------

CALL ARO_CONV_MNH(KIDIA,KFDIA,KLON,KLEV,LDEEP,&
 & LSHALLOW,LDIAGCONV,LSETTADJ,OTADJD,OTADJS,ZDTCONV,&
 & NSETENS,NIICE,LREFRESH_ALL,LDOWN,ZDXDY,&
 & PAPRSF,PZZF,&
 & ZT,PRV,PRC,&
 & PRI,PRHOREF,&
 & PU,PV,PW,&
 & PTTEN,PRVTEN,&
 & PRCTEN,&
 & PRITEN,&
 & PMF,&
 & PPRTEN,PPRSTEN,&
 & ZPRLFLX,ZPRSFLX,&
 & ZUMF,ZDMF,&
 & ZCAPE,ICLTOP,ICLBAS)  

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('BRI2ACCONV',1,ZHOOK_HANDLE)
END SUBROUTINE BRI2ACCONV
