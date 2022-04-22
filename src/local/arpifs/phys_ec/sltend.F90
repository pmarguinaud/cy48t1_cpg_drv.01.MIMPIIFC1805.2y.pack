#ifdef RS6K
@PROCESS HOT(NOVECTOR) NOSTRICT
#endif
SUBROUTINE SLTEND( YDMODEL,KIDIA, KFDIA, KLON, KLEV, KTEND,&
 & PT, PQ, PA, PRSF1, PRS1, &
 & PDYNU, PDYNV, PDYNT, PDYNQ, PDYNO3, PDYNA, &
 & PCENU, PCENV, PCENT, PCENQ, PCENO3, PCENA, &
 & PTENU, PTENV, PTENT, PTENQ, PTENO3, PTENA, &
 & PPHTENU, PPHTENV, PPHTENT, PPHTENGFL, &
 & PVFU, PVFV, PVFT, PVFQ, PVFA, &
 & PSAVTEND, PGFLSLP, PEXTRA, KFLDX)  


!**** *SLTEND * - compute sl tendencies for next time step and the 
!                 updates of the current profile. 
!               - remove supersaturation on final profiles and 
!                 save it for the next time step
!     PURPOSE.
!     --------
!               second order physics, computing the total dynamical and physical 
!               tendencies and the contribution of te physical tendencies
!               to be used in the next time step as part of the first guess
!               and as part of the final profile.

!**   Interface.
!     ----------
!        *CALL* *SLTEND*

! INPUT:

! KIDIA      : START OF HORIZONTAL LOOP
! KFDIA      : END   OF HORIZONTAL LOOP
! KLON       : HORIZONTAL DIMENSION
! KLEV       : END OF VERTICAL LOOP AND VERTICAL DIMENSION
! KTEND      : NUMBER OF PARAMETERS IN PSAVTEND
! KTRAC      : NUMBER OF ACTIVE TRACERS IN PVFC: NGHG+NAER+NGFL_EXT

! PT         : TEMPERATURE AT TIME t
! PQ         : HUMIDITY AT TIME t
! PA         : CLOUD FRACTION AT TIME t

! PRSF1      : PROVISIONAL T+DT PRESSURE ON FULL LEVELS
! PRS1       : PROVISIONAL T+DT PRESSURE ON HALF LEVELS

! time level (t+1)-(t)
! PDYNU      : EUL. DYNAMICAL TENDENCY OF U-COMP. OF WIND.
! PDYNV      : EUL. DYNAMICAL TENDENCY OF V-COMP. OF WIND.
! PDYNT      : EUL. DYNAMICAL TENDENCY OF TEMPERATURE.
! PDYNQ      : EUL. DYNAMICAL TENDENCY OF HUMIDITY
! PDYNO3     : EUL. DYNAMICAL TENDENCY OF OZONE MIXING RATIO (ECMWF PROG. OZONE)
! PDYNA      : EUL. DYNAMICAL TENDENCY OF CLOUD FRACTION

! (note: the updated total cloud liquid water and ice are computed in CLOUDSC !)

! time level (t) at departure point
! PPHTENU    : 0.5*PHYSICAL TENDENCY OF U-COMP. OF WIND
! PPHTENV    : 0.5*PHYSICAL TENDENCY OF V-COMP. OF WIND
! PPHTENT    : 0.5*PHYSICAL TENDENCY OF TEMPERATURE
! PPHTENGFL    : 0.5*PHYSICAL TENDENCY OF GFL FIELDS with LPHY=T

! arriv. point only at time level (t+1)
! PVFU       : ARRIV. POINT TENDENCY FOR VDIF + GWDRAG OF U-COMP. OF WIND
! PVFV       : ARRIV. POINT TENDENCY FOR VDIF + GWDRAG OF V-COMP. OF WIND
! PVFT       : ARRIV. POINT TENDENCY FOR VDIF OF TEMPERATURE
! PVFQ       : ARRIV. POINT TENDENCY FOR VDIF OF HUMIDITY
! PVFA       : ARRIV. POINT TENDENCY FOR VDIF OF CLOUD FRACTION

! UPDATED:

! input : arriv. point only at time level (t+1),
! PCENU      : TENDENCY OF U-COMP. OF WIND.
! PCENV      : TENDENCY OF V-COMP. OF WIND.
! PCENT      : TENDENCY OF TEMPERATURE.
! PCENQ      : TENDENCY OF HUMIDITY
! PCENO3     : TENDENCY OF OZONE MIXING RATIO (ECMWF PROG. OZONE)
! PCENL      : TENDENCY OF CLOUD WATER
! output: complete arriv. + departure point tendencies including the dynamical tendencies
! PTENU      : TENDENCY OF U-COMP. OF WIND.
! PTENV      : TENDENCY OF V-COMP. OF WIND.
! PTENT      : TENDENCY OF TEMPERATURE.
! PTENQ      : TENDENCY OF HUMIDITY
! PTENO3     : TENDENCY OF OZONE MIXING RATIO (ECMWF PROG. OZONE)
! PTENL      : TENDENCY OF CLOUD WATER

! PSAVTEND   : ARRAY OF TENDENCIES + AUX. SUPERSATURATION TO BE SAVED FOR NEXT TIME STEP

!-----------------------------------------------------------------------
!     Method. See documentation.
!   
!       This subroutine returns the following difference of slow physics tendencies for a field X:
!              PTENX=0.5*[-(dX^(+)/dt)_slowPhys+(dX^(-)/dt)_slowphys|D] where D denotes interpolation
!                                                                       to the departure point.
!       The above quantity is added in total (cummulative) tendency PCENX when update_state() which is 
!       called immediately after to return the required tendency when LSLPHY=T:
!              PTENX+PCENX=(dX^(+)/dt)_fastPhys + 0.5*[(dX^(+)/dt)_slowPhys + (dX^(-)/dt)_slowphys|D]   
!
!     -------
!     Modifications.
!     --------------
!     ORIGINAL 2001-06-25, Nils Wedi
!     M.Hamrud    : 03-08-01 GFL introduction
!     M.Hamrud    : 01-Oct-2003 CY28 Cleaning
!     A.Untch     : March-2004 Introduced EXTRA GFL fields (GFL_EXT)
!     M.Ko"hler   : 03-12-2004 VDF cloud tendency added for moist
!                              advection-diffusion PBL
!     A.Untch     : 12-03-2005 Aerosols as named GFL fields
!     J.Flemming  : 11-04-2005 Aerosols replaced with reactive gases
!     A.Tompkins  : CY31R2     Changes for supersaturation code
!     P.Bechtold  : 11-12-2005 Reorganize GFL/Tracer part
!     M.Janiskova : 21-12-2005 Modified condition for using cloud tendency
!     D.Salmond     22-Nov-2005 Mods for coarser/finer physics
!     S. Serrar     07-Sep-2006 tracers added for diagnostics (GEMS)
!     P.Bechtold    04-10-2006 Remove GFL tracers as SLTEND does not give positive definit results
!     A.Tompkins/R.Forbes : 15-03-2010  New CLV cloud variables added
!     R.Forbes    : 01-03-2011 Changed supersat liquid temperature threshold
!     R.Forbes    : 01-10-2011 Limited supersaturation to avoid excessive values
!     K. Yessad (July 2014): Move some variables.
!     M.Diamantakis/F. Vana  : 18-10-2013 More freedom to use GFL attributes and LSLPHY independently
!     J.Hague : 21-10-2014 Vector Optimisation for Cray
!     F. Vana  05-Mar-2015  Support for single precision
!     S. Malardel: Nov 2017 Fix pointers for PPHTENGFL to allow LSLPHY if ICI
!     M. Diamantakis (June 2018): Add LPHY control for q
!     R. Forbes 01-Nov-2018  Pass in arrays for cloud budget update 
!-----------------------------------------------------------------------

USE TYPE_MODEL, ONLY : MODEL
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK
USE YOMSLPHY  , ONLY : RSLWX
USE YOMCST    , ONLY : RG, RLVTT, RLSTT, RTT
USE YOETHF    , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 &                     R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RTWAT, RTICE, RTICECU, &
 &                     RTWAT_RTICE_R, RTWAT_RTICECU_R, RKOOP1, RKOOP2  
USE YOMCT3    , ONLY : NSTEP

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSF1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRS1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNO3(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDYNA(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCENU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCENV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCENT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCENQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCENO3(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCENA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENO3(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPHTENU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPHTENV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPHTENT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PPHTENGFL(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NUMFLDS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFU(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFV(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFT(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVFA(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSAVTEND(KLON,KLEV,KTEND) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGFLSLP(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NDIMSLP) 
! Extra fields for diagnostics
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEXTRA(KLON,KLEV,KFLDX) ! extra fields
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX ! Number of extra fields

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JK, JL, IK, ICALL, IU, IV, IT

LOGICAL :: LLFLAG(KLON)

! required to record supersaturation adjustments
REAL(KIND=JPRB) :: ZSUBQ(KLON)
REAL(KIND=JPRB) :: ZSUBT(KLON)
REAL(KIND=JPRB) :: ZSUBA(KLON)
REAL(KIND=JPRB) :: ZADJQ(KLON,KLEV)
REAL(KIND=JPRB) :: ZADJT(KLON,KLEV)
REAL(KIND=JPRB) :: ZADJA(KLON,KLEV)
REAL(KIND=JPRB) :: ZADJL(KLON,KLEV)
REAL(KIND=JPRB) :: ZADJI(KLON,KLEV)

! required for vertical integral cloud budget
REAL(KIND=JPRB) :: ZDP(KLON)
REAL(KIND=JPRB) :: ZRHO(KLON)
REAL(KIND=JPRB) :: ZDZ(KLON)

! adjustment variables
REAL(KIND=JPRB) :: ZQOLD(KLON),ZPP(KLON),ZQAD(KLON)
REAL(KIND=JPRB) :: ZCONS,ZCONS2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZFAC, ZFACCC, ZEPSEC, ZCCICETAU, &
                 & ZRG_R, ZQP1ENV

REAL(KIND=JPRB) :: ZZZFAC(KLON), ZZFAC(KLON), ZPT(KLON), ZZSUBA(KLON)

REAL(KIND=JPRB) :: ZEPSILON

INTEGER(KIND=JPIM) :: JX(KLON)
INTEGER(KIND=JPIM) :: JJ, JJJ

!-----------------------------------------------------------------------

#include "cuadjtq.intfb.h"
#include "fcttre.func.h"
#include "fccld.func.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SLTEND',0,ZHOOK_HANDLE)
ASSOCIATE(YGFL=>YDMODEL%YRML_GCONF%YGFL,YDECLDP=>YDMODEL%YRML_PHY_EC%YRECLDP,YDRIP=>YDMODEL%YRML_GCONF%YRRIP, &
 & YDSLPHY=>YDMODEL%YRML_PHY_G%YRSLPHY,YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2)
ASSOCIATE(NDIMSLP=>YGFL%NDIMSLP, NUMFLDS=>YGFL%NUMFLDS, YA=>YGFL%YA, &
 & YO3=>YGFL%YO3, YQ=>YGFL%YQ, &
 & NCLDTOP=>YDECLDP%NCLDTOP, NSSOPT=>YDECLDP%NSSOPT, RAMIN=>YDECLDP%RAMIN, &
 & RKOOPTAU=>YDECLDP%RKOOPTAU, RTHOMO=>YDECLDP%RTHOMO, &
 & LCLDBUD_VERTINT=>YDECLDP%LCLDBUD_VERTINT, &
 & NSTART=>YDRIP%NSTART, &
 & MSAT_SAVTEND=>YDSLPHY%MSAT_SAVTEND, TSPHY=>YDPHY2%TSPHY)
!     ------------------------------------------------------------------

ZEPSILON=100._JPRB*EPSILON(ZEPSILON)
ZEPSEC=1.E-14_JPRB
ZCONS=1.0_JPRB/TSPHY
ZCCICETAU=(1.0_JPRB-EXP(-TSPHY/RKOOPTAU)) ! factor for CC supersat adjustment
ZRG_R=1.0_JPRB/RG

ZCONS2=1.0_JPRB
IF (NSTEP==NSTART) THEN
  ZCONS2=0.0_JPRB
ENDIF

IU=1
IV=2
IT=3

LLFLAG(:)=.TRUE.


DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PSAVTEND(JL,JK,MSAT_SAVTEND)=0.0_JPRB
    ZADJT(JL,JK)=0.0_JPRB
    ZADJQ(JL,JK)=0.0_JPRB
    ZADJA(JL,JK)=0.0_JPRB
    ZADJL(JL,JK)=0.0_JPRB
    ZADJI(JL,JK)=0.0_JPRB
  ENDDO
ENDDO

!-----------------------------------------------------------------------
! note that the supersaturation adjustment is made with respect to 
! liquid saturation:  when T>0C 
! ice saturation:     when T<0C
!                     with an adjustment made to allow for ice 
!                     supersaturation in the clear sky
!
! Note also that the KOOP factor automatically clips the supersaturation
! to a maximum set by the liquid water saturation mixing ratio
! important for temperatures near to but below 0C
!
! This treatment mimics that at the start of CLOUDSC
!----------------------------------------------------------------------- 
DO JK=NCLDTOP,KLEV
  DO JL=KIDIA,KFDIA
    IF (NSTEP==NSTART) THEN
      ZSUBT(JL)=PT(JL,JK)+TSPHY*PCENT(JL,JK)
      ZSUBQ(JL)=PQ(JL,JK)+TSPHY*PCENQ(JL,JK)
      ZSUBA(JL)=PA(JL,JK)+TSPHY*PCENA(JL,JK)
    ELSE
      ZSUBT(JL)=PT(JL,JK)+RSLWX*TSPHY*PCENT(JL,JK)+&
       & TSPHY*(1.0_JPRB-RSLWX)*(PVFT(JL,JK)+PDYNT(JL,JK))+PPHTENT(JL,JK)
      IF (YQ%LPHY) THEN
        ZSUBQ(JL)=PQ(JL,JK)+RSLWX*TSPHY*PCENQ(JL,JK)+&
          & TSPHY*(1.0_JPRB-RSLWX)*(PVFQ(JL,JK)+PDYNQ(JL,JK))+PPHTENGFL(JL,JK,YQ%MP)
      ELSE
        ZSUBQ(JL)=PQ(JL,JK)+TSPHY*PCENQ(JL,JK)
      ENDIF
! MD When Yxx%LPHY false no memory is allocated for that part of the array and is not
!    needed by gpaddslphy()
      IF (YA%LPHY) THEN
        ZSUBA(JL)=PA(JL,JK)+RSLWX*TSPHY*PCENA(JL,JK)+&
         & TSPHY*(1.0_JPRB-RSLWX)*(PVFA(JL,JK)+PDYNA(JL,JK))+PPHTENGFL(JL,JK,YA%MP)  
      ELSE
        ZSUBA(JL)=PA(JL,JK)+TSPHY*PCENA(JL,JK)
      ENDIF
    ENDIF

    ZQOLD(JL)=ZSUBQ(JL)
    ZSUBT(JL)=MAX(ZSUBT(JL),160.0_JPRB) !safety
    ZPP(JL)=PRSF1(JL,JK)
  ENDDO

  IK=1
  ICALL=5
  CALL CUADJTQ &
   & (YDMODEL%YRML_PHY_EC%YRTHF, YDMODEL%YRCST, YDMODEL%YRML_PHY_SLIN%YREPHLI, KIDIA,    KFDIA,    KLON,     1,        IK,&
   &  ZPP,      ZSUBT,    ZSUBQ,    LLFLAG,   ICALL)  

!J---Vectorised code start----
! Gather Loop
  JJJ=0
  DO JL=KIDIA,KFDIA
    IF (PT(JL,JK) < RTT .AND. NSSOPT /= 0) THEN
      JJJ=JJJ+1
      JX(JJJ)=JL
      ZPT(JJJ)=PT(JL,JK)
      ZZSUBA(JJJ)=ZSUBA(JL)
    ENDIF
  ENDDO
! Vector Loop
  DO JJ=1,JJJ
    ZZFAC(JJ)=ZZSUBA(JJ)+FOKOOP(ZPT(JJ))*(1.0_JPRB-ZZSUBA(JJ))
  ENDDO
! Scatter Loop
  ZZZFAC(:)=1.0_JPRB
  DO JJ=1,JJJ
    JL=JX(JJ)
    ZZZFAC(JL)=ZZFAC(JJ)
  ENDDO
!J---Vectorised code end ----

  
!DEC$ IVDEP
  DO JL=KIDIA,KFDIA
    !----------------------------
    ! supersaturation adjustments
    !----------------------------
    IF (PT(JL,JK)>=RTT .OR. NSSOPT==0) THEN
! Moved above
!J    ZFAC=1.0_JPRB
      ZFACCC=1.0_JPRB
    ELSE
! Moved above
!J    ZFAC=ZSUBA(JL)+FOKOOP(PT(JL,JK))*(1.0_JPRB-ZSUBA(JL))
      ZFACCC=ZCCICETAU
    ENDIF
    ZFAC=ZZZFAC(JL)

    ! Calculate supersaturation to add to cloud
    IF (ZSUBA(JL) > 1.0_JPRB-RAMIN) THEN
      ZQAD(JL)=ZCONS*MAX(ZQOLD(JL)-ZFAC*ZSUBQ(JL),0.0_JPRB)
    ELSE
      ! Calculate environmental humidity supersaturation
      ZQP1ENV = (ZQOLD(JL) - ZSUBA(JL)*ZSUBQ(JL))/ &
       & SIGN(MAX(ABS(1.0_JPRB-ZSUBA(JL)),ZEPSILON),1.0_JPRB-ZSUBA(JL))
      ZQAD(JL)=ZCONS*MAX((1.0_JPRB-ZSUBA(JL))*(ZQP1ENV-ZFAC*ZSUBQ(JL)),0.0_JPRB)
    ENDIF 

    !-------------------------------------------------------------------
    ! Here the supersaturation is turned into liquid water
    ! However, if the temperature is below the threshold for homogeneous
    ! freezing then the supersaturation is turned instantly to ice.
    !--------------------------------------------------------------------
!    IF (ZQAD(JL)>ZEPSEC) THEN 
      IF (PT(JL,JK) > RTHOMO) THEN
        ZADJT(JL,JK)=RALVDCP*ZQAD(JL)
      ELSE
        ZADJT(JL,JK)=RALSDCP*ZQAD(JL)
      ENDIF
!    ENDIF
    ZADJQ(JL,JK)=ZQAD(JL)
    ! any surplus is added to the cloud at the next time-step
    ! to ensure consistent treatment from cloud scheme
    PSAVTEND(JL,JK,MSAT_SAVTEND)=ZQAD(JL)*TSPHY
  ENDDO

ENDDO

!---save physics tendencies for next time step(t+1-->t)

DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PSAVTEND(JL,JK,IU)=(1.0_JPRB-RSLWX)*(PCENU(JL,JK)-PVFU(JL,JK)-PDYNU(JL,JK))
    PSAVTEND(JL,JK,IV)=(1.0_JPRB-RSLWX)*(PCENV(JL,JK)-PVFV(JL,JK)-PDYNV(JL,JK))
    PSAVTEND(JL,JK,IT)=(1.0_JPRB-RSLWX)*(PCENT(JL,JK)-PVFT(JL,JK)-PDYNT(JL,JK))
    IF (YO3%LPHY) PGFLSLP(JL,JK,YO3%MPSLP)=(1.0_JPRB-RSLWX)*(PCENO3(JL,JK)-PDYNO3(JL,JK))
    IF (YQ%LPHY)PGFLSLP(JL,JK,YQ%MPSLP)=(1.0_JPRB-RSLWX)*(PCENQ(JL,JK)-PVFQ(JL,JK)-PDYNQ(JL,JK))
    IF (YA%LPHY)PGFLSLP(JL,JK,YA%MPSLP)=(1.0_JPRB-RSLWX)*(PCENA(JL,JK)-PVFA(JL,JK)-PDYNA(JL,JK))
  ENDDO
ENDDO


!---create updates for current time-step(t+1)
DO JK=1,KLEV
  DO JL=KIDIA,KFDIA
    PTENU(JL,JK)=-ZCONS2*PSAVTEND(JL,JK,IU)+ZCONS*PPHTENU(JL,JK)
    PTENV(JL,JK)=-ZCONS2*PSAVTEND(JL,JK,IV)+ZCONS*PPHTENV(JL,JK)
    ! add supersaturation adjustement to t+1 terms only
    PTENT(JL,JK)=-ZCONS2*PSAVTEND(JL,JK,IT)+ZADJT(JL,JK)&
     & +ZCONS*PPHTENT(JL,JK)  
! Need to reset cloud fraction to avoid incorrect accumulation 
! when PTENA updated by update_state()
    IF(YQ%LPHY) THEN
      PTENQ(JL,JK)=-ZCONS2*PGFLSLP(JL,JK,YQ%MPSLP)-ZADJQ(JL,JK)&
        & +ZCONS*PPHTENGFL(JL,JK,YQ%MP)        
    ELSE
      PTENQ(JL,JK)=0.0_JPRB
    ENDIF
    IF(YA%LPHY) THEN
      PTENA(JL,JK)=-ZCONS2*PGFLSLP(JL,JK,YA%MPSLP)+ZCONS*PPHTENGFL(JL,JK,YA%MP)
    ELSE
      PTENA(JL,JK)=0.0_JPRB
    ENDIF
    IF(YO3%LPHY) THEN
      PTENO3(JL,JK)=-ZCONS2*PGFLSLP(JL,JK,YO3%MPSLP)+ZCONS*PPHTENGFL(JL,JK,YO3%MP)
    ELSE
      PTENO3(JL,JK)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

!-----------------------------------------------------------------
! Vertical integral of all cloud process terms in one 3D field 
!-----------------------------------------------------------------
! The following is needed when the sequential saturation adjustment is included
!IF (LLCLDBUD_VERTINT) THEN
!  DO JK=1,KLEV
!    DO JL=KIDIA,KFDIA
! 
!      ZDP(JL)     = PRS1(JL,JK+1)-PRS1(JL,JK)   ! dp
!      ZRHO(JL)    = PRSF1(JL,JK)/(RD*PT(JL,JK)) ! p/RT air density
!      ZDZ(JL)     = ZDP(JL)/(ZRHO(JL)*RG)       ! Layer depth (m)

!      PEXTRA(JL,5,1)  = PEXTRA(JL,5,1)  + ZADJA(JL,JK)*ZDZ(JL)  ! + Supersat clipping for cloud fraction
!      PEXTRA(JL,21,1) = PEXTRA(JL,21,1) + ZADJL(JL,JK)*ZDZ(JL)  ! + Supersat clipping for liquid
!      PEXTRA(JL,46,1) = PEXTRA(JL,46,1) + ZADJI(JL,JK)*ZDZ(JL)  ! + Supersat clipping for ice
!    ENDDO
!  ENDDO
!ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SLTEND',1,ZHOOK_HANDLE)
END SUBROUTINE SLTEND
