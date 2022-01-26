SUBROUTINE SURF_IDEAL_FLUX(YDRIP,YDPHY0,YDPHYDS,LDAROME,KIDIA, KFDIA, KLON, PZ, PRHOA, PSFORC, &
          &       PTN,PTS, PLSM, PQS, PU, PV, PTHETAS, PSFTH, PSFTQ, PSFU, PSFV)

!     ------------------------------------------------------------------
!           Computes the surface fluxes for the temperature, 
!             vapor and  horizontal components of the wind  
!                         in ideal cases
!           Computes also a skin temperature consistent with the surface flux 
!     ------------------------------------------------------------------

!       Original : Yves Bouteloup
!             
!              modified by E. Bazile 18/07/2013 Ustar and Ts forcing
!                          E. Bazile 16/03/2016 Phasing and cleaning for Ts.


USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK
USE YOMRIP   , ONLY : TRIP
USE YOMLSFORC, ONLY : NT_SH_ADV_TIME,NSH_FORC_DEB,NSH_FORC_NUM, &
                    & NT_LH_ADV_TIME,NLH_FORC_DEB,NLH_FORC_NUM, &
                    & RZ0_FORC, &
                    & NT_US_ADV_TIME,NUS_FORC_DEB, NUS_FORC_NUM, &
                    & NT_TS_ADV_TIME,NTS_FORC_DEB, NTS_FORC_NUM, &
                    & LMUSCLFA,NMUSCLFA
USE YOMPHYDS , ONLY : TPHYDS

USE YOMCST   , ONLY : RG, RD, RV, RPI, RCPD, RLVTT, RALPD, RALPS, RALPW, &
 &          RBETD, RBETS, RBETW, RCPV, RCS, RCW, RETV, RGAMD, RGAMS, RGAMW, &
 &          RLSTT, RTT
USE YOMPHY0   ,ONLY : TPHY0

!     ------------------------------------------------------------------

IMPLICIT NONE
TYPE(TPHY0)       ,INTENT(IN)    :: YDPHY0
TYPE(TPHYDS)      ,INTENT(IN)    :: YDPHYDS
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSFORC(KLON,YDPHYDS%NSFORC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHOA (KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQS(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLSM(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTS(KLON)   ! Surface temperature 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTN(KLON)   ! PT(KLEV) in APLPAR
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTHETAS(KLON)  ! ThetaV en surface
LOGICAL           ,INTENT(IN)    :: LDAROME

REAL(KIND=JPRB)   ,INTENT(OUT)    :: PSFU(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PSFV(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PSFTH(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PSFTQ(KLON)

!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JLON ,JT, IINDA, IINDB, JITER
REAL(KIND=JPRB) :: ZA, ZB, ZEPS, ZINT, ZWIND(KLON), ZZ(KLON), ZZ0(KLON)
REAL(KIND=JPRB) :: ZUSTAR(KLON),ZLMO(KLON),ZRS(KLON),ZQS(KLON),ZQ0(KLON),ZE0(KLON),ZLH(KLON)

REAL(KIND=JPRB) :: ZSFBUO(KLON)   ! Buoyancy flux 
REAL(KIND=JPRB) :: ZVP1,ZVP2,ZVP3,ZVPT0,ZVS,ZMAV  ! Constants for Stevens case (to convert 
                                         ! buoyancy flux into latent and sensible heat flux
                                         
REAL(KIND=JPRB) :: ZE1,ZEP2,ZSFC,ZSFC_AIR,ZHFX,ZQFX
REAL(KIND=JPRB), SAVE :: ZSKT=288._JPRB

REAL(KIND=JPRB) :: ZHOOK_HANDLE,ZDELTA,ZUNDEF

!     ------------------------------------------------------------------
#include "fcttrm.func.h"
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SURF_IDEAL_FLUX',0,ZHOOK_HANDLE)
ASSOCIATE(VKARMN=>YDPHY0%VKARMN, &
 & NSFORC=>YDPHYDS%NSFORC, RSTATI=>YDRIP%RSTATI)
!     ------------------------------------------------------------------

! Lecture de la Ts prescrite utile pour calculer sigma*Ts**4
    IF (NTS_FORC_NUM > 0 ) THEN
      DO JLON=KIDIA,KFDIA
         PTS(JLON)=PSFORC(JLON,NTS_FORC_DEB+NTS_FORC_NUM-1)
      ENDDO
      IINDA = 1
      IINDB = 1
      DO JT=2,NTS_FORC_NUM
         !!! Compute tendency for RSTATI (current time) in between 
         !!! NT_TS_ADV_TIME(JT-1) and NT_TS_ADV_TIME(JT)
         IF (  NT_TS_ADV_TIME(JT)>=INT(RSTATI)) THEN
            IINDA = NTS_FORC_DEB + JT-2
            IINDB = NTS_FORC_DEB + JT-1
            ZINT=REAL(NT_TS_ADV_TIME(JT),JPRB)-REAL(NT_TS_ADV_TIME(JT-1),JPRB)
!DEC$ IVDEP
            DO JLON=KIDIA,KFDIA
              ZA = ( PSFORC(JLON,IINDB)-PSFORC(JLON,IINDA) )/ZINT
              ZB = PSFORC(JLON,IINDA)
              PTS(JLON) = ZA*(RSTATI-REAL(NT_TS_ADV_TIME(JT-1),JPRB))+ZB
            ENDDO
            EXIT
         ENDIF
      ENDDO
      DO JLON=KIDIA,KFDIA
         ZDELTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-PTS(JLON)))
         ZLH(JLON)   = FOLH (PTS(JLON),ZDELTA)
      ENDDO
    ELSE
      DO JLON=KIDIA,KFDIA
!         ZDELTA=MAX(0.0_JPRB,SIGN(1.0_JPRB,RTT-PTN(JLON)))
         ZDELTA=0. ! pour retrouver exactement les memes resultats mais pas propre !
         ZLH(JLON)   = FOLH (PTN(JLON),ZDELTA)
      ENDDO
    ENDIF ! fin de ts 

IF (NLH_FORC_NUM<0) THEN  ! Buoyancy forcing only one surface flux to be split into SH and LH

! Constants initialization
   
    ZVP1=0.6112_JPRB
    ZVP2=17.67_JPRB
    ZVP3=29.65_JPRB
    ZVPT0=273.15_JPRB
    ZVS=0.01_JPRB
    ZMAV=0.9_JPRB
    ZEP2 = RD/RV
    
! Surface forcing via buoyancy flux, to convert to sensible and latent heat flux
!  For Stevens case 
    ZDELTA=0
    DO JLON=KIDIA,KFDIA
       ZSFBUO(JLON) = PSFORC(JLON,NSH_FORC_DEB+NSH_FORC_NUM-1)
    ENDDO
! Presently no time modulation of surface buoyancy flux (only Stevens case)
    
! Computation of vapor saturation pressure

    ZE1  = ZVP1*EXP(ZVP2*(ZSKT-ZVPT0)/(ZSKT-ZVP3)) 
    ZSFC = ZEP2*ZE1/(100._JPRB-ZE1)
    ZSFC_AIR = ZMAV*ZSFC
    
! Computation of flux 

    ZQFX = ZVS*(ZSFC_AIR-PQS(1))
    ZHFX = ZSFBUO(1)*PTN(1)/RG - 0.608_JPRB*PTN(1)*ZQFX
    
    ZHFX = PRHOA(1)*RCPD*ZHFX

    PSFTH(:) = ZHFX

    ZQFX = ZQFX * PRHOA(1) * ZLH(1)
    PSFTQ(:) = ZQFX
    
ELSE  ! NLH_FORC_NUM>=0
!  Sensible heat  ==> enthalpy flux
    ZDELTA=0
    DO JLON=KIDIA,KFDIA
       PSFTH(JLON) = PSFORC(JLON,NSH_FORC_DEB+NSH_FORC_NUM-1)
       ZLH(JLON)   = FOLH (PTN(JLON),ZDELTA)
    ENDDO

    IINDA = 1
    IINDB = 1
    DO JT=2,NSH_FORC_NUM
!!! Compute tendency for RSTATI (current time) in between 
!!! NT_SH_ADV_TIME(JT-1) and NT_SH_ADV_TIME(JT)

      IF (  NT_SH_ADV_TIME(JT)>=INT(RSTATI)) THEN
      IINDA = NSH_FORC_DEB + JT-2
      IINDB = NSH_FORC_DEB + JT-1
      ZINT=REAL(NT_SH_ADV_TIME(JT),JPRB)-REAL(NT_SH_ADV_TIME(JT-1),JPRB)
!DEC$ IVDEP
       DO JLON=KIDIA,KFDIA
       
          ZA = ( PSFORC(JLON,IINDB)-PSFORC(JLON,IINDA) )/ZINT
          ZB = PSFORC(JLON,IINDA)
          PSFTH(JLON) = ZA*(RSTATI-REAL(NT_SH_ADV_TIME(JT-1),JPRB))+ZB

       ENDDO
      EXIT
      ENDIF
    ENDDO

! Latent heat  ==> q flux
    DO JLON=KIDIA,KFDIA
       PSFTQ(JLON) = PSFORC(JLON,NLH_FORC_DEB+NLH_FORC_NUM-1)
    ENDDO

    IINDA = 1
    IINDB = 1
    DO JT=2,NLH_FORC_NUM
!!! Compute tendency for RSTATI (current time) in between 
!!! NT_LH_ADV_TIME(JT-1) and NT_LH_ADV_TIME(JT)

      IF (  NT_LH_ADV_TIME(JT)>=INT(RSTATI)) THEN
      IINDA = NLH_FORC_DEB + JT-2
      IINDB = NLH_FORC_DEB + JT-1
      ZINT=REAL(NT_LH_ADV_TIME(JT),JPRB)-REAL(NT_LH_ADV_TIME(JT-1),JPRB)
!DEC$ IVDEP
       DO JLON=KIDIA,KFDIA
          ZA = ( PSFORC(JLON,IINDB)-PSFORC(JLON,IINDA) )/ZINT
          ZB = PSFORC(JLON,IINDA)
          PSFTQ(JLON) = ZA*(RSTATI-REAL(NT_LH_ADV_TIME(JT-1),JPRB))+ZB
       ENDDO
      EXIT
      ENDIF
    ENDDO

ENDIF  ! NLH_FORC_NUM<0

!  Wind Forcing

ZEPS = (RV-RD)/RD
DO JLON=KIDIA,KFDIA
  ZZ(JLON)    = PZ(JLON)/RG
  ZZ0(JLON)   = RZ0_FORC
  ZWIND(JLON) = SQRT(PU(JLON)*PU(JLON)+PV(JLON)*PV(JLON))
ENDDO  

!* water mixing ratio
!
ZRS(:) = 0._JPRB
ZQS(:)  = PQS(:) 

WHERE (ZQS(:)/=0._JPRB) ZRS(:) = 1._JPRB/(1._JPRB/ZQS(:) - 1._JPRB)

!* cinematic surface fluxes
ZQ0(:) = PSFTH(:) / RCPD / PRHOA(:)
ZE0(:) = PSFTQ(:)        / PRHOA(:) / ZLH(:)

ZUSTAR(:) = 0._JPRB
ZUNDEF = 1.E+20_JPRB
ZLMO(:) = ZUNDEF

! Ustar forcing
IF (NUS_FORC_NUM > 0 ) THEN
   DO JLON=KIDIA,KFDIA
      ZUSTAR(JLON)=PSFORC(JLON,NUS_FORC_DEB+NUS_FORC_NUM-1)
   ENDDO
   IINDA = 1
   IINDB = 1
   DO JT=2,NUS_FORC_NUM
      !!! Compute tendency for RSTATI (current time) in between 
      !!! NT_TS_ADV_TIME(JT-1) and NT_TS_ADV_TIME(JT)
      IF (  NT_US_ADV_TIME(JT)>=INT(RSTATI)) THEN
         IINDA = NUS_FORC_DEB + JT-2
         IINDB = NUS_FORC_DEB + JT-1
         ZINT=REAL(NT_US_ADV_TIME(JT),JPRB)-REAL(NT_US_ADV_TIME(JT-1),JPRB)
         DO JLON=KIDIA,KFDIA
            ZA = ( PSFORC(JLON,IINDB)-PSFORC(JLON,IINDA) )/ZINT
            ZB = PSFORC(JLON,IINDA)
            ZUSTAR(JLON) = ZA*(RSTATI-REAL(NT_US_ADV_TIME(JT-1),JPRB))+ZB
         ENDDO
         EXIT
      ENDIF
   ENDDO
ELSE
   ZUSTAR(:) = USTAR(YDPHY0,ZWIND(:),ZZ(:),ZZ0(:),ZLMO(:))
   ZUSTAR(:) = MAX ( ZUSTAR(:), 0.01_JPRB )

   DO JITER=1,10
      ZLMO  (:) = LMO  (YDPHY0,ZUSTAR(:),PTHETAS(:),ZRS(:),ZQ0(:),ZE0(:))
      ZUSTAR(:) = USTAR(YDPHY0,ZWIND(:),ZZ(:),ZZ0(:),ZLMO(:))
   ENDDO
   ZUSTAR(:) = MAX ( ZUSTAR(:), 0.01_JPRB )
ENDIF !! Fin forcage ustar

! Computation of skin temperature consistent with surface flux  

ZHFX = PSFTH(1)/PRHOA(1)/RCPD
ZSKT = PTN(1) + ZHFX/ZUSTAR(1)

IF(LMUSCLFA) THEN 
   CALL LFAECRR(NMUSCLFA,'SKT',ZSKT,1)
   CALL LFAECRR(NMUSCLFA,'LHSF',PSFTQ(1),1)
   CALL LFAECRR(NMUSCLFA,'SHSF',PSFTH(1),1)
ENDIF   

! Overwrite APLPAR surface temperature only over land !
IF (NTS_FORC_NUM==0) THEN
   DO JLON=KIDIA,KFDIA
      PTS(JLON)=ZSKT*PLSM(JLON)+(1.-PLSM(JLON))*PTS(JLON)
   ENDDO
ENDIF

! Some conversion ...

PSFU(:) = 0._JPRB
PSFV(:) = 0._JPRB
IF (LDAROME) THEN
    PSFTH(:) = PSFTH(:) / RCPD / PRHOA(:)
    PSFTQ(:) = PSFTQ(:) / ZLH(:)/ PRHOA(:)
   WHERE (ZWIND>0.)
     PSFU(:) = - PRHOA(:) * ZUSTAR(:)**2 * PU(:) / ZWIND(:)
     PSFV(:) = - PRHOA(:) * ZUSTAR(:)**2 * PV(:) / ZWIND(:)
   ENDWHERE
ELSE
   PSFTH(:) = -PSFTH(:) 
   PSFTQ(:) = -PSFTQ(:) / ZLH(:)
   WHERE (ZWIND>0.)
     PSFU(:) = + PRHOA(:) * ZUSTAR(:)**2 * PU(:) / ZWIND(:)
     PSFV(:) = + PRHOA(:) * ZUSTAR(:)**2 * PV(:) / ZWIND(:)
   ENDWHERE
ENDIF


END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURF_IDEAL_FLUX',1,ZHOOK_HANDLE)

CONTAINS 
!-------------------------------------------------------------------------------

FUNCTION PAULSON_PSIM(PZ_O_LMO)

  USE PARKIND1  ,ONLY : JPIM     ,JPRB
  USE YOMCST   , ONLY : RPI
  USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK
  
  REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)    :: PZ_O_LMO
  REAL(KIND=JPRB), DIMENSION(SIZE(PZ_O_LMO,1)) :: PAULSON_PSIM

  REAL(KIND=JPRB), DIMENSION(SIZE(PZ_O_LMO,1)) :: ZX
  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('SURF_IDEAL_FLUX:PAULSON_PSIM',0,ZHOOK_HANDLE)
  ZX=1.
  WHERE ( PZ_O_LMO(:) < 0. )
    ZX=(1.-15.*PZ_O_LMO)**(0.25)
    PAULSON_PSIM(:) = LOG( (1.+ZX**2)*(1+ZX)**2/8. ) - 2.*ATAN(ZX) + RPI/2.
  ELSEWHERE
    PAULSON_PSIM(:) = - 4.7 * PZ_O_LMO
  ENDWHERE
  IF (LHOOK) CALL DR_HOOK('SURF_IDEAL_FLUX:PAULSON_PSIM',1,ZHOOK_HANDLE)
END FUNCTION PAULSON_PSIM

!-------------------------------------------------------------------------------

FUNCTION LMO(YDPHY0,PUSTAR,PTHETA,PRV,PSFTH,PSFRV)

  USE PARKIND1  ,ONLY : JPIM     ,JPRB
  USE YOMCST   , ONLY : RG, RD, RV
  USE YOMPHY0   ,ONLY : TPHY0
  USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK
  
  TYPE(TPHY0)                  , INTENT(IN)   :: YDPHY0
  REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)   :: PUSTAR
  REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)   :: PTHETA
  REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)   :: PRV
  REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)   :: PSFTH
  REAL(KIND=JPRB), DIMENSION(:), INTENT(IN)   :: PSFRV
  REAL(KIND=JPRB), DIMENSION(SIZE(PUSTAR))    :: LMO

  REAL(KIND=JPRB), DIMENSION(SIZE(PUSTAR))   :: ZTHETAV
  REAL(KIND=JPRB)                            :: ZEPS,ZUNDEF
  REAL(KIND=JPRB) :: ZHOOK_HANDLE


  IF (LHOOK) CALL DR_HOOK('SURF_IDEAL_FLUX:LMO',0,ZHOOK_HANDLE)
  ASSOCIATE(VKARMN=>YDPHY0%VKARMN)
  ZEPS   = (RV-RD)/RD
  ZUNDEF = 1.E+20_JPRB

  ZTHETAV(:) = PTHETA(:) * ( 1._JPRB +ZEPS * PRV(:))

  LMO(:) = ZUNDEF
  WHERE ( PSFTH(:)/ZTHETAV(:)+ZEPS*PSFRV(:)/=0._JPRB )&
    & LMO(:) = - MAX(PUSTAR(:),1.E-6_JPRB)**3&
    & / ( VKARMN * (  RG / ZTHETAV(:)    * PSFTH(:)&
    &                         + RG * ZEPS * PSFRV(:) )  ) 

  WHERE(ABS(LMO)>10000._JPRB) LMO=ZUNDEF
  END ASSOCIATE
  IF (LHOOK) CALL DR_HOOK('SURF_IDEAL_FLUX:LMO',1,ZHOOK_HANDLE)

END FUNCTION LMO

!-------------------------------------------------------------------------------

FUNCTION USTAR(YDPHY0,PWIND,PZ,PZ0,PLMO)

  USE PARKIND1  ,ONLY : JPIM     ,JPRB
  USE YOMPHY0   ,ONLY : TPHY0
  USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK
  
  TYPE(TPHY0),INTENT(IN)                    :: YDPHY0
  REAL(KIND=JPRB), DIMENSION(:), INTENT(IN) :: PWIND
  REAL(KIND=JPRB), DIMENSION(:), INTENT(IN) :: PZ
  REAL(KIND=JPRB), DIMENSION(:), INTENT(IN) :: PZ0
  REAL(KIND=JPRB), DIMENSION(:), INTENT(IN) :: PLMO
  REAL(KIND=JPRB), DIMENSION(SIZE(PZ))      :: USTAR

  REAL(KIND=JPRB), DIMENSION(SIZE(PZ))  :: ZZ_O_LMO
  REAL(KIND=JPRB), DIMENSION(SIZE(PZ))  :: ZZ0_O_LMO
  REAL(KIND=JPRB)                       :: ZUNDEF
  REAL(KIND=JPRB) :: ZHOOK_HANDLE

!* purely unstable case
  IF (LHOOK) CALL DR_HOOK('SURF_IDEAL_FLUX:USTAR',0,ZHOOK_HANDLE)
  ASSOCIATE(VKARMN=>YDPHY0%VKARMN)
  ZUNDEF = 1.E+20_JPRB
  USTAR (:) = 0.
  ZZ_O_LMO (:) = ZUNDEF
  ZZ0_O_LMO(:) = ZUNDEF
  
!* general case
  WHERE(ABS(PLMO) > 1.E-20 .AND. PLMO/=ZUNDEF)
    ZZ_O_LMO  = PZ(:)  / PLMO(:)
    ZZ0_O_LMO = PZ0(:) / PLMO(:)
    USTAR(:) = PWIND&
            &     * VKARMN / ( LOG(PZ(:)/PZ0(:))&
            &     - PAULSON_PSIM(ZZ_O_LMO(:))&
            &                  + PAULSON_PSIM(ZZ0_O_LMO(:)) ) 
  ENDWHERE

!* purely neutral case
  WHERE(PLMO==ZUNDEF)
    ZZ_O_LMO = 0.
    USTAR(:) = PWIND&
           &        * VKARMN / LOG(PZ(:)/PZ0(:)) 
  ENDWHERE
  
  END ASSOCIATE
  IF (LHOOK) CALL DR_HOOK('SURF_IDEAL_FLUX:USTAR',1,ZHOOK_HANDLE)

END FUNCTION USTAR

END SUBROUTINE SURF_IDEAL_FLUX
