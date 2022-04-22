!OPTIONS XOPT(HSFUN)
#ifdef RS6K
@PROCESS HOT(NOVECTOR) NOSTRICT
#endif
SUBROUTINE CLOUD_SATADJ(YDECLDP, YDEPHLI, &
 & KIDIA,    KFDIA,    KLON,    KLEV, &
 & PTSPHY, PAP,  PAPH, &
 & PT, PQ, PA, &
 & PL, PI, & 
 & TENDENCY_T, TENDENCY_Q, TENDENCY_A, &
 & TENDENCY_L, TENDENCY_I,&
 & PHRSW,    PHRLW, &
 & PMFU,     PMFD, &
 & PVERVEL, &
 & TENDENCY_LOC_T, TENDENCY_LOC_Q, TENDENCY_LOC_A, &
 & TENDENCY_LOC_L, TENDENCY_LOC_I &
 & )
 
!===============================================================================
!**** *CLOUD_SATADJ* -  ROUTINE FOR PARAMATERIZATION OF CLOUD PROCESSES
!                  FOR PROGNOSTIC CLOUD SCHEME
!!
!     R.Forbes     (E.C.M.W.F.)
!!
!     PURPOSE
!     -------
!     Evaporation/condensation of cloud water in connection
!     with heating/cooling such as by subsidence/ascent
!!
!     INTERFACE.
!     ----------
!     *CLOUD_SATADJ* IS CALLED FROM *CALLPAR*
!     THE ROUTINE TAKES ITS INPUT FROM THE LONG-TERM STORAGE:
!     T,Q,L,PHI AND DETRAINMENT OF CLOUD WATER FROM THE
!     CONVECTIVE CLOUDS (MASSFLUX CONVECTION SCHEME), BOUNDARY
!     LAYER TURBULENT FLUXES OF HEAT AND MOISTURE, RADIATIVE FLUXES,
!     OMEGA.
!     IT MODIFIES TENDENCIES OF MODEL VARIABLES T AND Q
!        AS WELL AS CLOUD LIQUID, ICE, CLOUD FRACTION
!!
!     EXTERNALS.
!     ----------
!          NONE
!!
!     MODIFICATIONS.
!     -------------
!        01-10-2016 : R. Forbes  Duplicate of cond/evap from cloudsc.F90
!
!
!     REFERENCES.
!     ----------
!     Tietdke MWR 1993 - original description of the cloud parametrization
!     Jakob PhD 2000
!     Gregory et al. (2000) QJRMS
!     Tompkins el al. (2007) QJRMS - ice supersaturation parametrization
!!
!===============================================================================

USE YOECLDP  , ONLY : TECLDP
USE YOEPHLI  , ONLY : TEPHLI
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCST   , ONLY : RG, RD, RCPD, RETV, RTT, RLVTT, RLSTT, RLMLT, RV, YRCST
USE YOETHF   , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 & R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RALFDCP, RTWAT, RTICE, RTICECU, &
 & RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2, YRTHF

IMPLICIT NONE

!-------------------------------------------------------------------------------
!                 Declare input/output arguments
!-------------------------------------------------------------------------------
 
TYPE(TECLDP)      ,INTENT(IN)    :: YDECLDP
TYPE(TEPHLI)      ,INTENT(IN)    :: YDEPHLI
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON             ! Number of grid points
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV             ! Number of levels
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY           ! Physics timestep
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV)   ! Pressure on full levels
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1)! Pressure on half levels
! Initial state
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)    ! Temperature
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV)    ! Humidity
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV)    ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(IN)    :: PL(KLON,KLEV)    ! Cloud liquid
REAL(KIND=JPRB)   ,INTENT(IN)    :: PI(KLON,KLEV)    ! Cloud ice
! Input tendency to add to initial state
REAL(KIND=JPRB)   ,INTENT(IN)    :: TENDENCY_T(KLON,KLEV) ! Temperature
REAL(KIND=JPRB)   ,INTENT(IN)    :: TENDENCY_Q(KLON,KLEV) ! Humidity
REAL(KIND=JPRB)   ,INTENT(IN)    :: TENDENCY_A(KLON,KLEV) ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(IN)    :: TENDENCY_L(KLON,KLEV) ! Cloud liquid
REAL(KIND=JPRB)   ,INTENT(IN)    :: TENDENCY_I(KLON,KLEV) ! Cloud ice
! Input tendencies from other processes
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRSW(KLON,KLEV)   ! Short-wave heating rate
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHRLW(KLON,KLEV)   ! Long-wave heating rate
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFU(KLON,KLEV)    ! Conv. mass flux up
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMFD(KLON,KLEV)    ! Conv. mass flux down
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV) ! Vertical velocity
! Local output tendencies from cloud scheme
REAL(KIND=JPRB)   ,INTENT(OUT)   :: TENDENCY_LOC_T(KLON,KLEV) ! Temperature
REAL(KIND=JPRB)   ,INTENT(OUT)   :: TENDENCY_LOC_Q(KLON,KLEV) ! Humidity
REAL(KIND=JPRB)   ,INTENT(OUT)   :: TENDENCY_LOC_A(KLON,KLEV) ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(OUT)   :: TENDENCY_LOC_L(KLON,KLEV) ! Cloud liquid
REAL(KIND=JPRB)   ,INTENT(OUT)   :: TENDENCY_LOC_I(KLON,KLEV) ! Cloud ice

!-------------------------------------------------------------------------------
!                       Declare local variables
!-------------------------------------------------------------------------------

REAL(KIND=JPRB) :: ZTP1(KLON,KLEV)   
REAL(KIND=JPRB) :: ZQ(KLON,KLEV)
REAL(KIND=JPRB) :: ZA(KLON,KLEV)
REAL(KIND=JPRB) :: ZL(KLON,KLEV)
REAL(KIND=JPRB) :: ZI(KLON,KLEV)
REAL(KIND=JPRB) :: ZLI(KLON,KLEV)
REAL(KIND=JPRB) :: ZLIQFRAC(KLON,KLEV) ! cloud liquid water fraction: ql/(ql+qi)
REAL(KIND=JPRB) :: ZICEFRAC(KLON,KLEV) ! cloud ice water fraction: qi/(ql+qi)
REAL(KIND=JPRB) :: ZQSMIX(KLON,KLEV)   ! diagnostic mixed phase saturation 
REAL(KIND=JPRB) :: ZQSLIQ(KLON,KLEV)   ! liquid water saturation
REAL(KIND=JPRB) :: ZQSICE(KLON,KLEV)   ! ice water saturation
REAL(KIND=JPRB) :: ZFOEEWMT(KLON,KLEV)
REAL(KIND=JPRB) :: ZFOEEW(KLON,KLEV)
REAL(KIND=JPRB) :: ZFOEELIQT(KLON,KLEV)
REAL(KIND=JPRB) :: ZFOEALFA(KLON,KLEV+1)

REAL(KIND=JPRB) :: ZACOND(KLON)
REAL(KIND=JPRB) :: ZLCOND1(KLON) 
REAL(KIND=JPRB) :: ZLCOND2(KLON)
REAL(KIND=JPRB) :: ZLCOND1L(KLON) 
REAL(KIND=JPRB) :: ZLCOND2L(KLON)
REAL(KIND=JPRB) :: ZLCOND1I(KLON) 
REAL(KIND=JPRB) :: ZLCOND2I(KLON)
REAL(KIND=JPRB) :: ZLEVAPL(KLON)
REAL(KIND=JPRB) :: ZLEVAPI(KLON)
REAL(KIND=JPRB) :: ZFOKOOP(KLON)
REAL(KIND=JPRB) :: ZLICLD(KLON)
REAL(KIND=JPRB) :: ZLIQCLD(KLON)
REAL(KIND=JPRB) :: ZICECLD(KLON)
REAL(KIND=JPRB) :: ZDQS(KLON)
REAL(KIND=JPRB) :: ZTOLD(KLON)
REAL(KIND=JPRB) :: ZQOLD(KLON)  
REAL(KIND=JPRB) :: ZDTGDP(KLON) 
REAL(KIND=JPRB) :: ZRDTGDP(KLON)  
REAL(KIND=JPRB) :: ZLDEFR(KLON)
REAL(KIND=JPRB) :: ZRHO(KLON)
REAL(KIND=JPRB) :: ZGDP(KLON)
REAL(KIND=JPRB) :: ZDA(KLON)
REAL(KIND=JPRB) :: ZDP(KLON)
REAL(KIND=JPRB) :: ZSUPSAT(KLON)
REAL(KIND=JPRB) :: ZSUPSATL(KLON)
REAL(KIND=JPRB) :: ZSUPSATI(KLON)
REAL(KIND=JPRB) :: ZSUPSATA(KLON)
REAL(KIND=JPRB) :: ZDQSLIQDT(KLON)
REAL(KIND=JPRB) :: ZDQSICEDT(KLON)
REAL(KIND=JPRB) :: ZDQSMIXDT(KLON)
REAL(KIND=JPRB) :: ZCORQSLIQ(KLON)
REAL(KIND=JPRB) :: ZCORQSICE(KLON) 
REAL(KIND=JPRB) :: ZCORQSMIX(KLON)
REAL(KIND=JPRB) :: ZEVAPLIMLIQ(KLON)
REAL(KIND=JPRB) :: ZEVAPLIMICE(KLON)
REAL(KIND=JPRB) :: ZEVAPLIMMIX(KLON)

REAL(KIND=JPRB) :: ZANEW
REAL(KIND=JPRB) :: ZLEVAP  
REAL(KIND=JPRB) :: ZDTFORC
REAL(KIND=JPRB) :: ZDTDIAB
REAL(KIND=JPRB) :: ZMFDN
REAL(KIND=JPRB) :: ZALFA
REAL(KIND=JPRB) :: ZALFAW
REAL(KIND=JPRB) :: ZCOR
REAL(KIND=JPRB) :: ZCDMAX
REAL(KIND=JPRB) :: ZLCONDLIM
REAL(KIND=JPRB) :: ZDPMXDT
REAL(KIND=JPRB) :: ZDTDP
REAL(KIND=JPRB) :: ZEPSEC
REAL(KIND=JPRB) :: ZFAC
REAL(KIND=JPRB) :: ZFACI
REAL(KIND=JPRB) :: ZFACW
REAL(KIND=JPRB) :: ZQE
REAL(KIND=JPRB) :: ZQTMST
REAL(KIND=JPRB) :: ZRDCP
REAL(KIND=JPRB) :: ZRHC
REAL(KIND=JPRB) :: ZSIGK
REAL(KIND=JPRB) :: ZWTOT
REAL(KIND=JPRB) :: ZZDL
REAL(KIND=JPRB) :: ZZZDT
REAL(KIND=JPRB) :: ZQP1ENV
REAL(KIND=JPRB) :: ZTMPA
REAL(KIND=JPRB) :: ZEPSILON
REAL(KIND=JPRB) :: ZHOOK_HANDLE

LOGICAL :: LLFLAG(KLON)

INTEGER(KIND=JPIM) :: ICALL, IK, JK, JL

#include "cuadjtq.intfb.h"

!DIR$ VFUNCTION EXPHF
#include "fcttre.func.h"
#include "fccld.func.h"

!===============================================================================
IF (LHOOK) CALL DR_HOOK('CLOUD_SATADJ',0,ZHOOK_HANDLE)

ASSOCIATE(NCLDTOP=>YDECLDP%NCLDTOP, &
 & NSSOPT=>YDECLDP%NSSOPT, & 
 & RAMID=>YDECLDP%RAMID,   & 
 & RAMIN=>YDECLDP%RAMIN,   &
 & RLMIN=>YDECLDP%RLMIN,   &
 & RTHOMO=>YDECLDP%RTHOMO, &
 & RKOOPTAU=>YDECLDP%RKOOPTAU )

!===============================================================================


!######################################################################
!
!             1.  *** INITIAL VALUES FOR VARIABLES ***
!
!######################################################################

! Define a small number
ZEPSILON=100._JPRB*EPSILON(ZEPSILON)

! ---------------------
! Some simple constants
! ---------------------
ZQTMST  = 1.0_JPRB/PTSPHY
ZRDCP   = RD/RCPD
ZEPSEC  = 1.E-14_JPRB

! ------------------------------
! Updated state initialization 
! ------------------------------
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    ZTP1(JL,JK)= PT(JL,JK) + PTSPHY*TENDENCY_T(JL,JK)
    ZQ(JL,JK)  = PQ(JL,JK) + PTSPHY*TENDENCY_Q(JL,JK) 
    ZA(JL,JK)  = PA(JL,JK) + PTSPHY*TENDENCY_A(JL,JK)
    ZL(JL,JK)  = PL(JL,JK) + PTSPHY*TENDENCY_L(JL,JK)
    ZI(JL,JK)  = PI(JL,JK) + PTSPHY*TENDENCY_I(JL,JK)
  ENDDO
ENDDO

! ------------------------------------------------------------------
! Zero tendency arrays at model top levels where no cloud processes 
! ------------------------------------------------------------------
DO JK=1,NCLDTOP-1
  DO JL=KIDIA,KFDIA
    TENDENCY_LOC_T(JL,JK) = 0.0_JPRB
    TENDENCY_LOC_Q(JL,JK) = 0.0_JPRB 
    TENDENCY_LOC_A(JL,JK) = 0.0_JPRB
    TENDENCY_LOC_L(JL,JK) = 0.0_JPRB
    TENDENCY_LOC_I(JL,JK) = 0.0_JPRB
  ENDDO
ENDDO

! ------------------------------
! Define saturation values
! ------------------------------
DO JK=1,KLEV

    DO JL=KIDIA,KFDIA
      
      ! Note: ZQSMIX must be calculated on all levels for CUADJTQ
      
      !----------------------------------------
      ! old *diagnostic* mixed phase saturation
      !---------------------------------------- 
      ZFOEALFA(JL,JK) = FOEALFA(ZTP1(JL,JK))
      ZFOEEWMT(JL,JK) = MIN(FOEEWM(ZTP1(JL,JK))/PAP(JL,JK),0.5_JPRB)
      ZQSMIX(JL,JK)   = ZFOEEWMT(JL,JK)
      ZQSMIX(JL,JK)   = ZQSMIX(JL,JK)/(1.0_JPRB-RETV*ZQSMIX(JL,JK))

      !---------------------------------------------
      ! ice saturation T<273K
      ! liquid water saturation for T>273K 
      !---------------------------------------------
      ZALFA           = FOEDELTA(ZTP1(JL,JK))
      ZFOEEW(JL,JK)   = MIN((ZALFA*FOEELIQ(ZTP1(JL,JK))+ &
                        & (1.0_JPRB-ZALFA)*FOEEICE(ZTP1(JL,JK))) &
                        & /PAP(JL,JK),0.5_JPRB)
      ZQSICE(JL,JK)   = ZFOEEW(JL,JK)/(1.0_JPRB-RETV*ZFOEEW(JL,JK))

      !----------------------------------
      ! liquid water saturation
      !---------------------------------- 
      ZFOEELIQT(JL,JK)= MIN(FOEELIQ(ZTP1(JL,JK))/PAP(JL,JK),0.5_JPRB)
      ZQSLIQ(JL,JK)   = ZFOEELIQT(JL,JK)
      ZQSLIQ(JL,JK)   = ZQSLIQ(JL,JK)/(1.0_JPRB-RETV*ZQSLIQ(JL,JK))
    ENDDO

ENDDO ! on JK


!----------------------------------------------------------------------
!
!                   START OF VERTICAL LOOP OVER LEVELS
!
!----------------------------------------------------------------------

DO JK=NCLDTOP,KLEV

  !----------------------------------------------------------------------
  !
  ! 1. INITIALIZE VARIABLES
  !
  !----------------------------------------------------------------------

  
  DO JL=KIDIA,KFDIA

    ZLDEFR(JL)   = 0.0_JPRB                                
    
    !-------------------------
    ! derived variables needed
    !-------------------------

    ZDP(JL)     = PAPH(JL,JK+1)-PAPH(JL,JK)     ! dp
    ZGDP(JL)    = RG/ZDP(JL)                    ! g/dp
    ZRHO(JL)    = PAP(JL,JK)/(RD*ZTP1(JL,JK))   ! p/RT air density
    ZDTGDP(JL)  = PTSPHY*ZGDP(JL)               ! dt g/dp
    ZRDTGDP(JL) = ZDP(JL)*(1.0_JPRB/(PTSPHY*RG))  ! 1/(dt g/dp)

    !------------------------------------------
    ! Ensure cloud fraction is between 0 and 1
    !------------------------------------------
    ZA(JL,JK)=MAX(0.0_JPRB,MIN(1.0_JPRB,ZA(JL,JK)))

    !-------------------------------------------------------------------
    ! Calculate liq/ice fractions (no longer a diagnostic relationship)
    !-------------------------------------------------------------------
    ZLI(JL,JK)=ZL(JL,JK)+ZI(JL,JK)
    IF (ZLI(JL,JK)>RLMIN) THEN
      ZLIQFRAC(JL,JK)=ZL(JL,JK)/ZLI(JL,JK)
      ZICEFRAC(JL,JK)=1.0_JPRB-ZLIQFRAC(JL,JK)
    ELSE
      ZLIQFRAC(JL,JK)=0.0_JPRB
      ZICEFRAC(JL,JK)=0.0_JPRB
    ENDIF

    !------------------------------------
    ! Calculate dqs/dT correction factor
    !------------------------------------
    ! Reminder: RETV=RV/RD-1
    
    ! liquid
    ZFACW         = R5LES/((ZTP1(JL,JK)-R4LES)**2)
    ZCOR          = 1.0_JPRB/(1.0_JPRB-RETV*ZFOEELIQT(JL,JK))
    ZDQSLIQDT(JL) = ZFACW*ZCOR*ZQSLIQ(JL,JK)
    ZCORQSLIQ(JL) = 1.0_JPRB+RALVDCP*ZDQSLIQDT(JL)

    ! ice
    ZFACI         = R5IES/((ZTP1(JL,JK)-R4IES)**2)
    ZCOR          = 1.0_JPRB/(1.0_JPRB-RETV*ZFOEEW(JL,JK))
    ZDQSICEDT(JL) = ZFACI*ZCOR*ZQSICE(JL,JK)
    ZCORQSICE(JL) = 1.0_JPRB+RALSDCP*ZDQSICEDT(JL)

    ! diagnostic mixed
    ZALFAW        = ZFOEALFA(JL,JK)
    ZFAC          = ZALFAW*ZFACW+(1.0_JPRB-ZALFAW)*ZFACI
    ZCOR          = 1.0_JPRB/(1.0_JPRB-RETV*ZFOEEWMT(JL,JK))
    ZDQSMIXDT(JL) = ZFAC*ZCOR*ZQSMIX(JL,JK)
    ZCORQSMIX(JL) = 1.0_JPRB+FOELDCPM(ZTP1(JL,JK))*ZDQSMIXDT(JL)

    ! evaporation/sublimation limits
    ZEVAPLIMMIX(JL) = MAX((ZQSMIX(JL,JK)-ZQ(JL,JK))/ZCORQSMIX(JL),0.0_JPRB)
    ZEVAPLIMLIQ(JL) = MAX((ZQSLIQ(JL,JK)-ZQ(JL,JK))/ZCORQSLIQ(JL),0.0_JPRB)
    ZEVAPLIMICE(JL) = MAX((ZQSICE(JL,JK)-ZQ(JL,JK))/ZCORQSICE(JL),0.0_JPRB)

    !--------------------------------
    ! in-cloud consensate amount
    !--------------------------------
    ZTMPA = 1.0_JPRB/MAX(ZA(JL,JK),ZEPSEC)
    ZLIQCLD(JL) = ZL(JL,JK)*ZTMPA
    ZICECLD(JL) = ZI(JL,JK)*ZTMPA
    ZLICLD(JL)  = ZLIQCLD(JL)+ZICECLD(JL)

  ENDDO
    

  !======================================================================
  !
  !
  !  2.  SUPERSATURATION ADJUSTMENT DUE TO CHANGES IN WATER VAPOUR
  !
  !
  !======================================================================
  ! Note that the supersaturation adjustment is made with respect to 
  ! liquid saturation:  when T>0C 
  ! ice saturation:     when T<0C
  !                     with an adjustment made to allow for ice 
  !                     supersaturation in the clear sky
  ! Note also that the KOOP factor automatically clips the supersaturation
  ! to a maximum set by the liquid water saturation mixing ratio
  ! important for temperatures near to but below 0C
  !----------------------------------------------------------------------- 

!DIR$ NOFUSION
  !------------------------------------------------------------------------
  ! 3.1.1 Calculate Koop supersaturation limit function
  !  FOEELIQ(PTARE) = R2ES*EXP(R3LES*(PTARE-RTT)/(PTARE-R4LES))
  !  FOEEICE(PTARE) = R2ES*EXP(R3IES*(PTARE-RTT)/(PTARE-R4IES))
  !  Used in function: 
  !  FOKOOP(PTARE) = MIN(RKOOP1-RKOOP2*PTARE,FOEELIQ(PTARE)/FOEEICE(PTARE))
  ! Needs to be set for all temperatures
  !------------------------------------------------------------------------
  DO JL=KIDIA,KFDIA
    ZFOKOOP(JL)=FOKOOP(ZTP1(JL,JK))
  ENDDO

  DO JL=KIDIA,KFDIA

    IF (ZTP1(JL,JK)>=RTT .OR. NSSOPT==0) THEN
      ZFAC  = 1.0_JPRB
      ZFACI = 1.0_JPRB
    ELSE
      ZFAC  = ZA(JL,JK)+ZFOKOOP(JL)*(1.0_JPRB-ZA(JL,JK))
      ZFACI = PTSPHY/RKOOPTAU
    ENDIF

    !-------------------------------------------------------------------
    ! 3.1.2 Calculate supersaturation wrt Koop including dqs/dT 
    !       correction factor
    !-------------------------------------------------------------------

    ! Calculate supersaturation to add to cloud
    IF (ZA(JL,JK) > 1.0_JPRB-RAMIN) THEN
      ZSUPSAT(JL) = MAX((ZQ(JL,JK)-ZFAC*ZQSICE(JL,JK))/ZCORQSICE(JL)&
     &                  ,0.0_JPRB)
    ELSE
      ! Calculate environmental humidity supersaturation
      ZQP1ENV = (ZQ(JL,JK) - ZA(JL,JK)*ZQSICE(JL,JK))/ &
     & SIGN(MAX(ABS(1.0_JPRB-ZA(JL,JK)),ZEPSILON),1.0_JPRB-ZA(JL,JK))
      ZSUPSAT(JL) = MAX((1.0_JPRB-ZA(JL,JK))*(ZQP1ENV-ZFAC*ZQSICE(JL,JK))&
     &                  /ZCORQSICE(JL),0.0_JPRB)
    ENDIF 
    
    !-------------------------------------------------------------------
    ! Here the supersaturation is turned into liquid water
    ! However, if the temperature is below the threshold for homogeneous
    ! freezing then the supersaturation is turned instantly to ice.
    !--------------------------------------------------------------------

    IF (ZSUPSAT(JL) > ZEPSEC) THEN

      IF (ZTP1(JL,JK) > RTHOMO) THEN
        ! Turn supersaturation into liquid water        
        ZSUPSATL(JL) = ZSUPSAT(JL)
        ZSUPSATI(JL) = 0.0_JPRB
      ELSE
        ! Turn supersaturation into ice water        
        ZSUPSATI(JL) = ZSUPSAT(JL)
        ZSUPSATL(JL) = 0.0_JPRB
      ENDIF

      ! Increase cloud amount using RKOOPTAU timescale
      ZSUPSATA(JL) = (1.0_JPRB-ZA(JL,JK))*ZFACI
 
    ELSE

      ZSUPSATL(JL) = 0.0_JPRB
      ZSUPSATI(JL) = 0.0_JPRB
      ZSUPSATA(JL) = 0.0_JPRB

    ENDIF

  ENDDO ! on JL


  !======================================================================
  !
  !
  !  3.  CONDENSATION/EVAPORATION DUE TO DQSAT/DT
  !
  !
  !======================================================================
  !  calculate dqs/dt
  !  Note: For the separate prognostic Qi and Ql, one would ideally use
  !  Qsat/DT wrt liquid/Koop here, since the physics is that new clouds
  !  forms by liquid droplets [liq] or when aqueous aerosols [Koop] form.
  !  These would then instantaneous freeze if T<-38C or lead to ice growth 
  !  by deposition in warmer mixed phase clouds.  However, since we do 
  !  not have a separate prognostic equation for in-cloud humidity or a 
  !  statistical scheme approach in place, the depositional growth of ice 
  !  in the mixed phase can not be modelled and we resort to supersaturation  
  !  wrt ice instanteously converting to ice over one timestep 
  !  (see Tompkins et al. QJRMS 2007 for details)
  !  Thus for the initial implementation the diagnostic mixed phase is 
  !  retained for the moment, and the level of approximation noted.  
  !----------------------------------------------------------------------

  DO JL=KIDIA,KFDIA
    ZDTDP   = ZRDCP*ZTP1(JL,JK)/PAP(JL,JK)
    ZDPMXDT = ZDP(JL)*ZQTMST
    ZMFDN   = 0.0_JPRB
    IF(JK < KLEV) ZMFDN=PMFU(JL,JK+1)+PMFD(JL,JK+1)
    ZWTOT   = PVERVEL(JL,JK)+0.5_JPRB*RG*(PMFU(JL,JK)+PMFD(JL,JK)+ZMFDN)
    ZWTOT   = MIN(ZDPMXDT,MAX(-ZDPMXDT,ZWTOT))
    ZZZDT   = PHRSW(JL,JK)+PHRLW(JL,JK)
    ZDTDIAB = MIN(ZDPMXDT*ZDTDP,MAX(-ZDPMXDT*ZDTDP,ZZZDT))&
                    & *PTSPHY+RALFDCP*ZLDEFR(JL)  
! Note: ZLDEFR should be set to the difference between the mixed phase functions
! in the convection and cloud scheme, but this is not calculated, so is zero and
! the functions must be the same
    ZDTFORC = ZDTDP*ZWTOT*PTSPHY+ZDTDIAB
    ZQOLD(JL)   = ZQSMIX(JL,JK)
    ZTOLD(JL)   = ZTP1(JL,JK)
    ZTP1(JL,JK) = ZTP1(JL,JK)+ZDTFORC
    ZTP1(JL,JK) = MAX(ZTP1(JL,JK),160.0_JPRB)
    LLFLAG(JL)  = .TRUE.
  ENDDO

  IK=JK
  ICALL=5
  CALL CUADJTQ &
   & (YRTHF, YRCST, YDEPHLI,KIDIA,    KFDIA,   KLON,     KLEV,    IK,&
   & PAP(KIDIA,JK), ZTP1,    ZQSMIX,   LLFLAG,  ICALL)  

  DO JL=KIDIA,KFDIA
    ZDQS(JL)      = ZQSMIX(JL,JK)-ZQOLD(JL)
    ZQSMIX(JL,JK) = ZQOLD(JL)
    ZTP1(JL,JK)   = ZTOLD(JL)
  ENDDO

  !----------------------------------------------------------------------
  ! 3a  ZDQS(JL) > 0:  EVAPORATION OF CLOUDS
  ! ----------------------------------------------------------------------
  ! Assumes a delta function of cloud condensate, no change to cloud cover

  DO JL=KIDIA,KFDIA

    IF (ZDQS(JL) > 0.0_JPRB) THEN

      ZLEVAP = ZA(JL,JK)*MIN(ZDQS(JL),ZLICLD(JL))
      ZLEVAP = MIN(ZLEVAP,ZEVAPLIMMIX(JL))
      ZLEVAP = MIN(ZLEVAP,MAX(ZQSMIX(JL,JK)-ZQ(JL,JK),0.0_JPRB))

      ! For first guess call
      ZLEVAPL(JL) = ZLIQFRAC(JL,JK)*ZLEVAP
      ZLEVAPI(JL) = ZICEFRAC(JL,JK)*ZLEVAP

    ELSE
    
      ZLEVAPL(JL) = 0.0_JPRB
      ZLEVAPI(JL) = 0.0_JPRB
    
    ENDIF

  ENDDO

  !----------------------------------------------------------------------
  ! 3b ZDQS(JL) < 0: FORMATION OF CLOUDS
  !----------------------------------------------------------------------

  ! (1) Increase of cloud water in existing clouds

  DO JL=KIDIA,KFDIA
    
    ! Initialise output arrays
    ZLCOND1L(JL) = 0.0_JPRB
    ZLCOND1I(JL) = 0.0_JPRB

    IF(ZA(JL,JK) > ZEPSEC.AND.ZDQS(JL) <= -RLMIN) THEN

      ZLCOND1(JL)=MAX(-ZDQS(JL),0.0_JPRB) !new limiter

      !old limiter (significantly improves upper tropospheric humidity rms)
      IF(ZA(JL,JK) > 0.99_JPRB) THEN
        ZCOR=1.0_JPRB/(1.0_JPRB-RETV*ZQSMIX(JL,JK))
        ZCDMAX=(ZQ(JL,JK)-ZQSMIX(JL,JK))/&
         & (1.0_JPRB+ZCOR*ZQSMIX(JL,JK)*FOEDEM(ZTP1(JL,JK)))  
      ELSE
        ZCDMAX=(ZQ(JL,JK)-ZA(JL,JK)*ZQSMIX(JL,JK))/ZA(JL,JK)
      ENDIF
      ZLCOND1(JL)=MAX(MIN(ZLCOND1(JL),ZCDMAX),0.0_JPRB)
      ! end old limiter
      
      ZLCOND1(JL)=ZA(JL,JK)*ZLCOND1(JL)
      IF(ZLCOND1(JL) < RLMIN) ZLCOND1(JL)=0.0_JPRB
      
      !-------------------------------------------------------------------------
      ! All increase goes into liquid unless so cold cloud homogeneously freezes
      ! Include new liquid formation in first guess value, otherwise liquid 
      ! remains at cold temperatures until next timestep.
      !-------------------------------------------------------------------------
      IF (ZTP1(JL,JK)>RTHOMO) THEN
        ZLCOND1L(JL) = ZLCOND1(JL)
      ELSE
        ZLCOND1I(JL) = ZLCOND1(JL)
      ENDIF
    ENDIF
  ENDDO

  ! (2) Generation of new clouds (da/dt>0)
  
  DO JL=KIDIA,KFDIA

    ZLCOND2L(JL) = 0.0_JPRB
    ZLCOND2I(JL) = 0.0_JPRB
    ZACOND(JL)   = 0.0_JPRB
    
    IF(ZDQS(JL) <= -RLMIN .AND. ZA(JL,JK)<1.0_JPRB-ZEPSEC) THEN

      !---------------------------
      ! Critical relative humidity
      !---------------------------
      ZRHC=RAMID
      ZSIGK=PAP(JL,JK)/PAPH(JL,KLEV+1)
      ! Increase RHcrit to 1.0 towards the surface (eta>0.8)
      IF(ZSIGK > 0.8_JPRB) THEN
        ZRHC=RAMID+(1.0_JPRB-RAMID)*((ZSIGK-0.8_JPRB)/0.2_JPRB)**2
      ENDIF

      !---------------------------
      ! Supersaturation options
      !---------------------------      
      IF (NSSOPT==0) THEN 
        ! No scheme
        ZQE=(ZQ(JL,JK)-ZA(JL,JK)*ZQSICE(JL,JK))/&
            & MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))  
        ZQE=MAX(0.0_JPRB,ZQE)
      ELSEIF (NSSOPT==1) THEN 
        ! Tompkins 
        ZQE=(ZQ(JL,JK)-ZA(JL,JK)*ZQSICE(JL,JK))/&
            & MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))  
        ZQE=MAX(0.0_JPRB,ZQE)
      ELSEIF (NSSOPT==2) THEN 
        ! Lohmann and Karcher
        ZQE=ZQ(JL,JK)  
      ELSEIF (NSSOPT==3) THEN 
        ! Gierens
        ZQE=ZQ(JL,JK)+ZLI(JL,JK)
      ENDIF

      IF (ZTP1(JL,JK)>=RTT .OR. NSSOPT==0) THEN 
        ! No ice supersaturation allowed
        ZFAC=1.0_JPRB        
      ELSE
        ! Ice supersaturation
        ZFAC=ZFOKOOP(JL)
      ENDIF

      IF(ZQE >= ZRHC*ZQSICE(JL,JK)*ZFAC.AND.ZQE<ZQSICE(JL,JK)*ZFAC) THEN
        ! note: not **2 on 1-a term if ZQE is used. 
        ! Added correction term ZFAC to numerator 15/03/2010
        ZACOND(JL)=-(1.0_JPRB-ZA(JL,JK))*ZFAC*ZDQS(JL)/&
         &MAX(2.0_JPRB*(ZFAC*ZQSICE(JL,JK)-ZQE),ZEPSEC)

        ZACOND(JL)=MIN(ZACOND(JL),1.0_JPRB-ZA(JL,JK))  !PUT THE LIMITER BACK

        ! Linear term:
        ! Added correction term ZFAC 15/03/2010
        ZLCOND2(JL)=-ZFAC*ZDQS(JL)*0.5_JPRB*ZACOND(JL) !mine linear

        ! new limiter formulation
        ZZDL=2.0_JPRB*(ZFAC*ZQSICE(JL,JK)-ZQE)/MAX(ZEPSEC,1.0_JPRB-ZA(JL,JK))
        ! Added correction term ZFAC 15/03/2010
        IF (ZFAC*ZDQS(JL)<-ZZDL) THEN
          ! ZLCONDLIM=(ZA(JL,JK)-1.0_JPRB)*ZDQS(JL)-ZQSICE(JL,JK)+ZQ(JL,JK)
          ZLCONDLIM=(ZA(JL,JK)-1.0_JPRB)*ZFAC*ZDQS(JL)- &
     &               ZFAC*ZQSICE(JL,JK)+ZQ(JL,JK)
          ZLCOND2(JL)=MIN(ZLCOND2(JL),ZLCONDLIM)
        ENDIF
        ZLCOND2(JL)=MAX(ZLCOND2(JL),0.0_JPRB)

        IF(ZLCOND2(JL) < RLMIN .OR. (1.0_JPRB-ZA(JL,JK))<ZEPSEC ) THEN
          ZLCOND2(JL) = 0.0_JPRB
          ZACOND(JL) = 0.0_JPRB
        ENDIF
        IF(ZLCOND2(JL) == 0.0_JPRB) ZACOND(JL)=0.0_JPRB

        !------------------------------------------------------------------------
        ! All increase goes into liquid unless so cold cloud homogeneously freezes
        ! Include new liquid formation in first guess value, otherwise liquid 
        ! remains at cold temperatures until next timestep.
        !------------------------------------------------------------------------
        IF (ZTP1(JL,JK)>RTHOMO) THEN
          ZLCOND2L(JL) = ZLCOND2(JL)
        ELSE
          ZLCOND2I(JL) = ZLCOND2(JL)
        ENDIF

      ENDIF
    ENDIF
  ENDDO

  !------------------------------------------------------------------
  ! Add tendencies for first guess call of cloud scheme only
  !------------------------------------------------------------------

  DO JL=KIDIA,KFDIA
  
    ! Temperature tendency
    TENDENCY_LOC_T(JL,JK) = &
        &   RALVDCP*(ZLCOND1L(JL)+ZLCOND2L(JL)+ZSUPSATL(JL)-ZLEVAPL(JL))*ZQTMST &
        & + RALSDCP*(ZLCOND1I(JL)+ZLCOND2I(JL)+ZSUPSATI(JL)-ZLEVAPI(JL))*ZQTMST
  
    ! Humidity tendency
    TENDENCY_LOC_Q(JL,JK) = &
        &  (ZLEVAPL(JL)-ZLCOND1L(JL)-ZLCOND2L(JL)-ZSUPSATL(JL) &
        & + ZLEVAPI(JL)-ZLCOND1I(JL)-ZLCOND2I(JL)-ZSUPSATI(JL))*ZQTMST

    ! Cloud cover (ensure between bounds of 0 and 1)
    ZANEW = ZA(JL,JK)+ZACOND(JL)+ZSUPSATA(JL)
    ZANEW = MIN(ZANEW,1.0_JPRB)
    IF (ZANEW<RAMIN) ZANEW=0.0_JPRB
    ZDA(JL)=ZANEW-ZA(JL,JK)
    TENDENCY_LOC_A(JL,JK) = ZDA(JL)*ZQTMST

    ! Cloud liquid and ice phase
    TENDENCY_LOC_L(JL,JK) = (ZLCOND1L(JL)+ZLCOND2L(JL)+ZSUPSATL(JL) &
        & - ZLEVAPL(JL))*ZQTMST
    TENDENCY_LOC_I(JL,JK) = (ZLCOND1I(JL)+ZLCOND2I(JL)+ZSUPSATI(JL) &
        & - ZLEVAPI(JL))*ZQTMST

  ENDDO

ENDDO ! on vertical level JK

!===============================================================================
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CLOUD_SATADJ',1,ZHOOK_HANDLE)
END SUBROUTINE CLOUD_SATADJ
