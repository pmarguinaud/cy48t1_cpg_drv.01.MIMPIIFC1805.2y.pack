SUBROUTINE ACNORGWD(YDNORGWD, KIDIA, KFDIA, KLON, KLEV, PDTIME, &
                    & PP, PGEMU, PTT, PUU, PVV, PVOVO, PPREC, &
                    & PTEND_U, PTEND_V)


!**** *ACNORGWD * - NON-OROGRAPHIC GRAVITY WAVE DRAG

! SUJET.
!     ------

! Parametrization of the momentum flux deposition due to a discrete number of gravity waves


!**   INTERFACE.
!     ----------
!        *CALL* *ACNORGWD*

!-----------------------------------------------------------------------

! -   INPUT ARGUMENTS.
!     ----------------

! KIDIA                 : INDICE DE DEPART DES BOUCLES VECTORISEES SUR L'HORIZONT..
! KFDIA                 : INDICE DE FIN DES BOUCLES VECTORISEES SUR L'HORIZONTALE.
! KLON                  : DIMENSION HORIZONTALE DES TABLEAUX.
! KTDIA                 : INDICE DE DEPART DES BOUCLES VERTICALES (1 EN GENERAL).
! KLEV                  : DIMENSION VERTICALE DES TABLEAUX "FULL LEVEL".
! PDTIME                : PAS DE TEMPS DE LA PHYSIQUE (TSPHY)
! PP(KLON, KLEV)        : FULL-LEVEL PRESSURE (convention : 1=SURFACE, KLEV=TOA)
! PGEMU(KLON)           : SINE OF GEOGRAPHICAL LATITUDE  
! PVOVO(KLON, KLEV)     : FULL-LEVEL RELATIVE VORTICITY
! PTT(KLON, KLEV)       : FULL-LEVEL TEMPERATURE
! PUU(KLON, KLEV)       : FULL-LEVEL ZONAL WIND zonal wind at full levels
! PVV(KLON, KLEV)       : FULL-LEVEL MERID WIND
! PPREC(KLON)           : PRECIPITATION

!-----------------------------------------------------------------------

! -   OUTPUT ARGUMENTS
!     ----------------


!-----------------------------------------------------------------------

! -   INPUT/OUTPUT ARGUMENTS
!     ----------------------

! PTEND_U(KLON, KLEV)   : FULL-LEVEL TENDENCY ON ZONAL WIND (m/s2)
! PTEND_V(KLON, KLEV)   : FULL-LEVEL TENDENCY ON MERID WIND (m/s2)

!-----------------------------------------------------------------------

! METHOD.
!     --------

! Based on :
! Lott et al. (2012) A stochastic parameterization of non-orographic gravity waves
! Lott et Guez (2013) A stochastic parameterization of non-orographic gravity waves due to convection
! De la Camara et Lott (2015) A parameterization of gravity waves emitted by fronts and jets


! AUTHOR.
!     -------
! 2012-09, F. Lott


! MODIFICATIONS.
!     --------------
! 2016-01, D. Saint-Martin : adaptation for ARPEGE (norm doctor)
!                            inclusion of extratropical sources
! 2016-08                  : sources from frontogenesis only extratropical
! 2017-02, J.F. Gueremy    : numerical security in Log argument for LL_USE_NORMAL_C
! 2017-11, D. Saint-Martin : numerical optimization
! 2018-01, D. Saint-Martin : diagnostics XIOS
! 2018-09, D. Saint-Martin : adpated for CY46
!-----------------------------------------------------------------------
    
  USE PARKIND1,  ONLY : JPIM, JPRB
  USE YOMHOOK,   ONLY : LHOOK, DR_HOOK
  USE YOMCST,    ONLY : RG, RD, RCPD, RPI, RLVTT, ROMEGA
  USE YOMNORGWD, ONLY : TNORGWD

!#ifdef WXIOS
!  USE MODD_IO_SURF_ARO, ONLY : LXIOS
!#endif

  IMPLICIT NONE

  TYPE(TNORGWD)     , INTENT(IN)    :: YDNORGWD
  INTEGER(KIND=JPIM), INTENT(IN)    :: KIDIA
  INTEGER(KIND=JPIM), INTENT(IN)    :: KFDIA
  INTEGER(KIND=JPIM), INTENT(IN)    :: KLON
  INTEGER(KIND=JPIM), INTENT(IN)    :: KLEV
  REAL(KIND=JPRB),    INTENT(IN)    :: PDTIME
  REAL(KIND=JPRB),    INTENT(IN)    :: PP(KLON, KLEV)
  REAL(KIND=JPRB),    INTENT(IN)    :: PGEMU(KLON)
  REAL(KIND=JPRB),    INTENT(IN)    :: PVOVO(KLON, KLEV)
  REAL(KIND=JPRB),    INTENT(IN)    :: PTT(KLON, KLEV)
  REAL(KIND=JPRB),    INTENT(IN)    :: PUU(KLON, KLEV)
  REAL(KIND=JPRB),    INTENT(IN)    :: PVV(KLON, KLEV)
  REAL(KIND=JPRB),    INTENT(IN)    :: PPREC(KLON)
  REAL(KIND=JPRB),    INTENT(INOUT) :: PTEND_U(KLON, KLEV)
  REAL(KIND=JPRB),    INTENT(INOUT) :: PTEND_V(KLON, KLEV)

!-----------------------------------------------------------------------
  
  INTEGER(KIND=JPIM) :: JLON, JLEV

  INTEGER(KIND=JPIM), PARAMETER :: INK = 2, INP = 2, INO = 2, INW = INK * INP * INO
  INTEGER(KIND=JPIM) :: JK, JP, JO, JW
  
  REAL(KIND=JPRB) :: ZCPHA ! Intrinsic/absolute phase velocity
  REAL(KIND=JPRB) :: ZK(INW, KLON) ! Horizontal wavenumber amplitude
  REAL(KIND=JPRB) :: ZP(INW, KLON) ! Horizontal wavenumber angle
  REAL(KIND=JPRB) :: ZPCOS(INW, KLON) ! Cosinus of horizontal wavenumber angle
  REAL(KIND=JPRB) :: ZPSIN(INW, KLON) ! Sinus of horizontal wavenumber angle
  REAL(KIND=JPRB) :: ZO(INW, KLON) ! Absolute frequency

  REAL(KIND=JPRB) :: ZOM(INW, KLON) ! Waves intr. freq. at the half-lev surrounding the full level
  REAL(KIND=JPRB) :: ZOP(INW, KLON)
  
  REAL(KIND=JPRB) :: ZWWM(INW, KLON) ! Waves vertical vel. at the 2 half-lev surrounding the full level
  REAL(KIND=JPRB) :: ZWWP(INW, KLON)
  
  INTEGER(KIND=JPIM) :: ILAUNCH ! Launching altitude index
  INTEGER(KIND=JPIM) :: ILTROP ! tropo altitude

  REAL(KIND=JPRB) :: ZRUWP(INW, KLON) ! temporary x-fluxes for each wave at half-levels
  REAL(KIND=JPRB) :: ZRVWP(INW, KLON) ! temporary y-fluxes for each wave at half-levels
  REAL(KIND=JPRB) :: ZRUW(KLON, KLEV+1) ! total x-fluxes at half-levels
  REAL(KIND=JPRB) :: ZRVW(KLON, KLEV+1) ! total y-fluxes at half-levels

  REAL(KIND=JPRB) :: ZABSRUW(KLON, KLEV+1)  ! averaged absolute x-fluxes at half-levels
  REAL(KIND=JPRB) :: ZD_RUW(KLON, KLEV)     ! net x-fluxes at full levels (ARPEGE convention)
  REAL(KIND=JPRB) :: ZZD_RUW(KLON, KLEV)    ! net x-fluxes at full levels
  REAL(KIND=JPRB) :: ZD_ABSRUW(KLON, KLEV)  ! averaged absolute x-fluxes at full levels (ARPEGE convention)
  REAL(KIND=JPRB) :: ZZD_ABSRUW(KLON, KLEV) ! averaged absolute x-fluxes at full levels

  REAL(KIND=JPRB) :: ZH0 ! Characteristic height of the atmosphere
  REAL(KIND=JPRB) :: ZPR ! Reference pressure
  REAL(KIND=JPRB) :: ZTR ! Reference temperature

  REAL(KIND=JPRB) :: ZH(KLON, KLEV+1)  ! Log-pressure altitude
  REAL(KIND=JPRB) :: ZUH(KLON, KLEV+1) ! Zonal winds at half-levels 
  REAL(KIND=JPRB) :: ZVH(KLON, KLEV+1) ! Merid winds at half-levels
  REAL(KIND=JPRB) :: ZPH(KLON, KLEV+1) ! Pressure at half-levels
  REAL(KIND=JPRB) :: ZBV(KLON, KLEV+1) ! Brunt-Vaisala freq. at half-levels
  REAL(KIND=JPRB) :: ZBVLOW(KLON)      ! Brunt-Vaisala freq. averaged between launch and ltrop altitude
    
  REAL(KIND=JPRB), PARAMETER :: ZPSEC = 1.E-6 ! Security to avoid division by 0 pressure
  REAL(KIND=JPRB), PARAMETER :: ZBVSEC = 1.E-5 ! Security to avoid negative BVF
  REAL(KIND=JPRB), PARAMETER :: ZOISEC = 1.E-6 ! Security to avoid 0 intrinsic frequency
  REAL(KIND=JPRB), PARAMETER :: ZSHSEC = 1.E-12 ! Security to avoid 0 wind vertical shear
  REAL(KIND=JPRB), PARAMETER :: ZRICMX = 1.E3 ! Max Richardson number
  REAL(KIND=JPRB), PARAMETER :: ZCOSEC = 1.E-8 ! Security to avoid 0 Coriolis parameter

  REAL(KIND=JPRB) :: ZSHEAR(KLON, KLEV) ! zonal wind vertical shear 
  REAL(KIND=JPRB) :: ZRIC(KLON, KLEV)   ! Richardson number
  REAL(KIND=JPRB) :: ZFL(KLON) ! Launched EP flux from front and jets sources 
  REAL(KIND=JPRB) :: ZFCOR(KLON) ! Coriolis parameter 
  REAL(KIND=JPRB) :: ZX

  REAL(KIND=JPRB) :: ZRUW0_RAND ! Launched EP flux from random tropical sources
  REAL(KIND=JPRB) :: ZRUW0_PREC ! Launched EP flux from convective sources
  REAL(KIND=JPRB) :: ZRUW0_BACK ! Launched EP flux from 'background' sources
  REAL(KIND=JPRB) :: ZRUW0_FRON ! Launched EP flux from 'front and jets' sources

  REAL(KIND=JPRB) :: ZZRUW0_RAND(KLON) ! Total launched EP flux from random tropical sources
  REAL(KIND=JPRB) :: ZZRUW0_PREC(KLON) ! Total launched EP flux from convective sources
  REAL(KIND=JPRB) :: ZZRUW0_BACK(KLON) ! Total launched EP flux from 'background' sources
  REAL(KIND=JPRB) :: ZZRUW0_FRON(KLON) ! Total launched EP flux from 'front and jets' sources

  REAL(KIND=JPRB) :: ZALPH_RAND ! Fraction of launched EP flux coming from random tropical sources
  REAL(KIND=JPRB) :: ZALPH_PREC ! Fraction of launched EP flux coming from convective sources
  REAL(KIND=JPRB) :: ZALPH_BACK ! Fraction of launched EP flux coming from 'background' sources
  REAL(KIND=JPRB) :: ZALPH_FRON ! Fraction of launched EP flux coming from 'front and jets' sources

  REAL(KIND=JPRB) :: ZBV3
  REAL(KIND=JPRB) :: ZOM4
  REAL(KIND=JPRB) :: ZOP3
  REAL(KIND=JPRB) :: ZTMP1
  REAL(KIND=JPRB) :: ZTMP2

  ! Temporary : to save old 'random' settings
  LOGICAL :: LL_SMOOTH_BV
  LOGICAL :: LL_SECU_BV
  LOGICAL :: LL_USE_INTR_FREQ
  LOGICAL :: LL_USE_NEW_ANGLE
  LOGICAL :: LL_USE_NORMAL_C

  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  
!-----------------------------------------------------------------------
 
#include "abor1.intfb.h"
! #include "mse_xios_send.h"

!-----------------------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('ACNORGWD',0,ZHOOK_HANDLE)

  ASSOCIATE(NORGWD_SCHEME => YDNORGWD%NORGWD_SCHEME, &
    & NORGWD_RUWMAX => YDNORGWD%NORGWD_RUWMAX, &
    & NORGWD_SAT => YDNORGWD%NORGWD_SAT, &
    & NORGWD_RDISS => YDNORGWD%NORGWD_RDISS, &
    & NORGWD_DELTAT => YDNORGWD%NORGWD_DELTAT, &
    & NORGWD_KMIN => YDNORGWD%NORGWD_KMIN, &
    & NORGWD_KMAX => YDNORGWD%NORGWD_KMAX, &
    & NORGWD_CMIN => YDNORGWD%NORGWD_CMIN, & 
    & NORGWD_CMAX => YDNORGWD%NORGWD_CMAX, &
    & NORGWD_NLAUNCH => YDNORGWD%NORGWD_NLAUNCH, &
    & NORGWD_PRMAX => YDNORGWD%NORGWD_PRMAX, &
    & NORGWD_DZ => YDNORGWD%NORGWD_DZ, &
    & NORGWD_NTROPO => YDNORGWD%NORGWD_NTROPO, &
    & NORGWD_GB => YDNORGWD%NORGWD_GB, & 
    & NORGWD_GFRON => YDNORGWD%NORGWD_GFRON, &
    & NORGWD_DZFRON => YDNORGWD%NORGWD_DZFRON)

  !-----------------------------------------------------------------
  !  1. INITIALISATIONS
  !-----------------------------------------------------------------

    IF (NORGWD_DELTAT < PDTIME) THEN
      CALL ABOR1('ACNORGWD : NORGWD_DELTAT < PDTIME !')
    ENDIF

    IF (KLEV < INW) THEN
      CALL ABOR1('ACNORGWD : PROBLEM WITH RANDOM NUMBERS !')
    ENDIF

    IF (NORGWD_SCHEME == 'RAND') THEN
      LL_SMOOTH_BV     = .FALSE.
      LL_SECU_BV       = .FALSE.
      LL_USE_INTR_FREQ = .FALSE.
      LL_USE_NEW_ANGLE = .FALSE.
      LL_USE_NORMAL_C = .FALSE.
      ZALPH_RAND     = 1.0
      ZALPH_PREC     = 0.0
      ZALPH_BACK     = 0.0
      ZALPH_FRON     = 0.0
    ELSEIF (NORGWD_SCHEME == 'PREC') THEN
      LL_SMOOTH_BV     = .TRUE.
      LL_SECU_BV       = .TRUE.
      LL_USE_INTR_FREQ = .TRUE.
      LL_USE_NEW_ANGLE = .TRUE.
      LL_USE_NORMAL_C  = .TRUE.
      ZALPH_RAND     = 0.0
      ZALPH_PREC     = 1.0
      ZALPH_BACK     = 0.0
      ZALPH_FRON     = 0.0
    ELSEIF (NORGWD_SCHEME == 'ALLS') THEN
      LL_SMOOTH_BV     = .TRUE.
      LL_SECU_BV       = .TRUE.
      LL_USE_INTR_FREQ = .TRUE.
      LL_USE_NEW_ANGLE = .TRUE.
      LL_USE_NORMAL_C  = .TRUE.      
      ZALPH_RAND     = 0.0
      ZALPH_PREC     = 1.0
      ZALPH_BACK     = 1.0
      ZALPH_FRON     = 1.0
    ELSE
      CALL ABOR1('ACNORGWD : NORGWD_SCHEME NOT SPECIFIED')
    ENDIF

    ! NORGWD_RDISS : Dissipation coefficient
    ! NORGWD_SAT : Saturation parameter S_c
    ! NORGWD_DZ : characteristic depth of the source
    ! NORGWD_PRMAX : maximum of rain for which the theory applies (in kg/m^2/s)
    ! NORGWD_KMIN : min horizontal wavenumbers
    ! NORGWD_KMAX : max horizontal wavenumbers
    ! NORGWD_CMIN : min absolute/intrinsic phase vel.
    ! NORGWD_CMAX : max absolute/intrinsic phase vel.
    ! NORGWD_DELTAT : time scale of the life cycle of the waves parameterized
    ! NORGWD_RUWMAX : max EP-Flux at launch altitude (use for random or convective GW waves)
    ! NORGWD_DZFRON : characteristic depth of the source (front and jets source)
    ! NORGWD_GFRON : parameter G_0 (~1) that controls the amplitude of the EP flux emitted by fronts and jets
    ! NORGWD_GB : parameter that controls the amplitude of the EP flux emitted by 'background' sources

    ILAUNCH = NORGWD_NLAUNCH
    ILTROP  = NORGWD_NTROPO

    ZTR = 240.0_JPRB
    ZPR = 101300.0_JPRB
    ZH0 = RD * ZTR / RG

    !-------------------------------------------------------------    
    ! 2. EVALUATION OF THE BACKGROUND FLOW AT HALF-LEVELS
    !-------------------------------------------------------------

    ! Pressure and inverse of pressure at half-levels
    DO JLON = KIDIA, KFDIA
      DO JLEV = 2, KLEV
        ZPH(JLON, JLEV) = EXP((LOG(PP(JLON, JLEV)) + LOG(PP(JLON, JLEV - 1))) / 2.)
      ENDDO
      ZPH(JLON, KLEV + 1) = 0.
      ZPH(JLON, 1) = 2. * PP(JLON, 1) - ZPH(JLON, 2)
    ENDDO

    ! Log pressure vert. coordinate
    DO JLON = KIDIA, KFDIA
      DO JLEV = 1, KLEV + 1 
        ZH(JLON, JLEV) = ZH0 * LOG(ZPR / (ZPH(JLON, JLEV) + ZPSEC))
      ENDDO
    ENDDO

    ! BV frequency (ZBV) at half-levels
    ! ZUH USED IS AS A TEMPORARY ARRAY DOWN TO WINDS
    DO JLON = KIDIA, KFDIA 
      DO JLEV = 2, KLEV
        ZUH(JLON, JLEV) = 0.5 * (PTT(JLON, JLEV) + PTT(JLON, JLEV - 1)) &
         & * RD**2 / RCPD / ZH0**2 + (PTT(JLON, JLEV) &
         & - PTT(JLON, JLEV - 1)) / (ZH(JLON, JLEV) - ZH(JLON, JLEV - 1)) * RD / ZH0
      ENDDO
      ZBVLOW(JLON) = 0.5 * (PTT(JLON, ILTROP) + PTT(JLON, ILAUNCH)) &
       & * RD**2 / RCPD / ZH0**2 + (PTT(JLON, ILTROP) &
       & - PTT(JLON, ILAUNCH))/(ZH(JLON, ILTROP)- ZH(JLON, ILAUNCH)) * RD / ZH0
      ZUH(JLON, 1) = ZUH(JLON, 2)
      ZUH(JLON, KLEV + 1) = ZUH(JLON, KLEV)
      ZBV(JLON, 1) = ZUH(JLON, 2)
      ZBV(JLON, KLEV + 1) = ZUH(JLON, KLEV)
      DO JLEV = 2, KLEV
        IF (LL_SMOOTH_BV) THEN
          ZBV(JLON, JLEV) = ( ZUH(JLON, JLEV+1) + 2.*ZUH(JLON, JLEV) + ZUH(JLON, JLEV-1) ) / 4.
        ELSE
          ZBV(JLON, JLEV) = ZUH(JLON, JLEV)
        ENDIF
      ENDDO
    ENDDO  
    IF (LL_SECU_BV) THEN
      ZBV = MAX(SQRT(MAX(ZBV, 0.)), ZBVSEC)
      ZBVLOW = MAX(SQRT(MAX(ZBVLOW, 0.)), ZBVSEC)
    ELSE
      ZBV = SQRT(MAX(ZBV, ZBVSEC))
      ZBVLOW = SQRT(MAX(ZBVLOW, ZBVSEC))
    ENDIF

    ! WINDS (ZUH and ZVH) at half-levels
    DO JLON = KIDIA, KFDIA 
      DO JLEV = 2, KLEV
        ZUH(JLON, JLEV) = 0.5 * (PUU(JLON, JLEV) + PUU(JLON, JLEV - 1))
        ZVH(JLON, JLEV) = 0.5 * (PVV(JLON, JLEV) + PVV(JLON, JLEV - 1))
      ENDDO
      ZUH(JLON, 1) = 0.
      ZVH(JLON, 1) = 0.
      ZUH(JLON, KLEV + 1) = PUU(JLON, KLEV)
      ZVH(JLON, KLEV + 1) = PVV(JLON, KLEV)
    ENDDO

    !-------------------------------------------------------------
    ! 3. WAVES CHARACTERISTICS CHOSEN RANDOMLY
    !-------------------------------------------------------------

    ! The mod function of weird arguments are used to produce the waves characteristics
    ! in an almost stochastic way
    JW = 0
    DO JP = 1, INP
      DO JK = 1, INK
        DO JO = 1, INO
          JW = JW + 1
          DO JLON = KIDIA, KFDIA
            ! Angle (0 or PI so far)
            IF (LL_USE_NEW_ANGLE) THEN
              ZP(JW, JLON) = (SIGN(1.0_JPRB, 0.5_JPRB - MOD(PTT(JLON, JW) * 10._JPRB, 1.0_JPRB)) + 1.0_JPRB) &
               &   * RPI / 2.0_JPRB
            ELSE
              ZP(JW, JLON) = 2.0_JPRB * RPI * REAL(JP - 1.0_JPRB) / REAL(INP)
            ENDIF 
            ZPCOS(JW, JLON) = COS(ZP(JW, JLON))
            ZPSIN(JW, JLON) = SIN(ZP(JW, JLON))
            ! Horizontal wavenumber amplitude
            ZK(JW, JLON) = NORGWD_KMIN + (NORGWD_KMAX - NORGWD_KMIN) * MOD(PTT(JLON, JW) * 100.0_JPRB, 1.0_JPRB)
            ! Horizontal phase speed
            IF (LL_USE_NORMAL_C) THEN
              ! Horizontal phase speed : normal sampling with cmax as standard dev
              ZTMP1 = MAX(ZSHSEC,MIN(1.0_JPRB-ZSHSEC,MOD(PTT(JLON, JW)**2, 1.0_JPRB)))
              ZTMP2 = MOD(PTT(JLON, JW)**3, 1.0_JPRB)
              ZCPHA = NORGWD_CMAX * SQRT(-2.0_JPRB * LOG(ZTMP1))*COS(2.0_JPRB*RPI*ZTMP2)
            ELSE
              ZCPHA = NORGWD_CMIN + (NORGWD_CMAX - NORGWD_CMIN) * MOD(PTT(JLON, JW)**2, 1.0_JPRB)
            ENDIF
            ! Absolute frequency is imposed
            ZO(JW, JLON) = ZCPHA * ZK(JW, JLON)
            ! Intrinsic frequency is imposed
            IF (LL_USE_INTR_FREQ) THEN 
              ZO(JW, JLON) = ZO(JW, JLON) &
               &    + ZK(JW, JLON) * ZPCOS(JW, JLON) * ZUH(JLON, ILAUNCH) &
               &    + ZK(JW, JLON) * ZPSIN(JW, JLON) * ZVH(JLON, ILAUNCH)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    !-------------------------------------------------------------
    ! 4. COMPUTE THE FLUXES
    !-------------------------------------------------------------

    ! 4.1  EP flux at launching altitude

    ! Richardson number (ZRIC)
    ! only zonal wind shear is used
    DO JLON = KIDIA, KFDIA
      DO JLEV = 1, KLEV
        ZSHEAR(JLON, JLEV) = ((ZUH(JLON, JLEV + 1) - ZUH(JLON, JLEV)) / (ZH(JLON, JLEV + 1) - ZH(JLON, JLEV)))**2
        ZSHEAR(JLON, JLEV) = MAX(ZSHSEC, ZSHEAR(JLON, JLEV))
      ENDDO
      DO JLEV = 1, KLEV
        ZRIC(JLON, JLEV) = MIN(ZRICMX, (ZBV(JLON, JLEV)*ZBV(JLON, JLEV)) / ZSHEAR(JLON, JLEV))
      ENDDO
    ENDDO

    ! Launch EP-flux emitted by fronts and jets
    ! formula (4) in de la Camara and Lott (GRL, 2015) 
    DO JLON = KIDIA, KFDIA 
      ZFL(JLON) = 0.0
      ZFCOR(JLON) = 2.0_JPRB * ROMEGA * PGEMU(JLON)
      ZFCOR(JLON) = MAX(ABS(ZFCOR(JLON)), ZCOSEC)
      DO JLEV = 1, KLEV
        ZX = MIN(PVOVO(JLON, JLEV)**2 / ZFCOR(JLON), ZFCOR(JLON))
        ZFL(JLON) = ZFL(JLON) + (ZH(JLON, JLEV+1) - ZH(JLON, JLEV)) &
         & *EXP(- ZH(JLON, JLEV) / ZH0)*ZBV(JLON, JLEV)*ZX*EXP(-RPI*SQRT(ZRIC(JLON, JLEV)))
      ENDDO
      ZFL(JLON) = 0.25*NORGWD_GFRON*NORGWD_DZFRON*ZFL(JLON)
    ENDDO

    DO JLON = KIDIA, KFDIA 
      
      ZZRUW0_RAND(JLON) = 0.0_JPRB
      ZZRUW0_PREC(JLON) = 0.0_JPRB
      ZZRUW0_BACK(JLON) = 0.0_JPRB
      ZZRUW0_FRON(JLON) = 0.0_JPRB
 
      DO JW = 1, INW
        ! Evaluate intrinsic frequency at launching altitude:
        ZOP(JW, JLON) = ZO(JW, JLON) &
         & - ZK(JW, JLON) * ZPCOS(JW, JLON) * ZUH(JLON, ILAUNCH) &
         & - ZK(JW, JLON) * ZPSIN(JW, JLON) * ZVH(JLON, ILAUNCH)

        ! RANDOM TROPICAL SOURCES (Lott et al., GRL, 2012)
        ZRUW0_RAND = NORGWD_RUWMAX &
         & * MOD(100. * (PUU(JLON, JW)**2 + PVV(JLON, JW)**2), 1.0_JPRB) &
         & * COS(ASIN(PGEMU(JLON)))**8
        
        ! CONVECTIVE SOURCES (Lott and Guez, JGR, 2013)
        ZRUW0_PREC = NORGWD_RUWMAX &
         & * (RD / RCPD / ZH0 * RLVTT * NORGWD_PRMAX * TANH(PPREC(JLON) / NORGWD_PRMAX))**2 &
         & * ZK(JW, JLON)**3 / NORGWD_KMIN / ZBVLOW(JLON) / MAX(ABS(ZOP(JW, JLON)), ZOISEC)**3 &
         & * EXP(- ZBVLOW(JLON)**2 / MAX(ABS(ZOP(JW, JLON)), ZOISEC)**2 * ZK(JW, JLON)**2 &
         & * NORGWD_DZ**2)
          
        ! RANDOM 'BACKGROUND' SOURCES (de la Camara et al, JGR, 2014)
        ZRUW0_BACK = NORGWD_GB * MOD(100. * (PUU(JLON, JW)**2 + PVV(JLON, JW)**2), 1.0)

        ! FRONTS AND JETS SOURCES (de la Camara and Lott, GRL, 2015)
        ZRUW0_FRON = (1 - COS(ASIN(PGEMU(JLON)))**12)*ZFL(JLON)
          
        ! SUM ALL SOURCES
        ZWWP(JW, JLON) = ZALPH_RAND*ZRUW0_RAND + ZALPH_PREC*ZRUW0_PREC + ZALPH_BACK*ZRUW0_BACK + ZALPH_FRON*ZRUW0_FRON

        ! DIAGNOSE ABSOLUTE LAUNCHED MOMENTUM FLUXES (Pa)
        ZZRUW0_RAND(JLON) = ZZRUW0_RAND(JLON) + ZALPH_RAND*ZRUW0_RAND
        ZZRUW0_PREC(JLON) = ZZRUW0_PREC(JLON) + ZALPH_PREC*ZRUW0_PREC
        ZZRUW0_BACK(JLON) = ZZRUW0_BACK(JLON) + ZALPH_BACK*ZRUW0_BACK
        ZZRUW0_FRON(JLON) = ZZRUW0_FRON(JLON) + ZALPH_FRON*ZRUW0_FRON

        ZRUWP(JW, JLON) = SIGN(1.0_JPRB, ZOP(JW, JLON))*ZPCOS(JW, JLON) * ZWWP(JW, JLON)
        ZRVWP(JW, JLON) = SIGN(1.0_JPRB, ZOP(JW, JLON))*ZPSIN(JW, JLON) * ZWWP(JW, JLON)
      ENDDO
      
      ! DIAGNOSE AVERAGED ABSOLUTE LAUNCHED MOMENTUM FLUXES (Pa)
      ZZRUW0_RAND(JLON) = ZZRUW0_RAND(JLON) / INW
      ZZRUW0_PREC(JLON) = ZZRUW0_PREC(JLON) / INW
      ZZRUW0_BACK(JLON) = ZZRUW0_BACK(JLON) / INW
      ZZRUW0_FRON(JLON) = ZZRUW0_FRON(JLON) / INW

    ENDDO

    !  4.2 Uniform values below the launching altitude

    DO JLON = KIDIA, KFDIA 
      DO JLEV = 1, ILAUNCH
        ZRUW(JLON, JLEV) = 0.0
        ZRVW(JLON, JLEV) = 0.0
        DO JW = 1, INW
          ZRUW(JLON, JLEV) = ZRUW(JLON, JLEV) + ZRUWP(JW, JLON)
          ZRVW(JLON, JLEV) = ZRVW(JLON, JLEV) + ZRVWP(JW, JLON)
        ENDDO
      ENDDO
    ENDDO

    !  4.3 Loop over altitudes, with passage from one level to the
    !      next done by i) conserving the EP flux, ii) dissipating
    !      a little, iii) testing critical levels, and vi) testing
    !      the breaking.
    !      W(KB)ARNING: AJLEV THE PHYSICS IS HERE

    DO JLON = KIDIA, KFDIA
      DO JLEV = ILAUNCH, KLEV - 1

        ZBV3 = ((ZBV(JLON, JLEV + 1) + ZBV(JLON, JLEV)) / 2.)   &
              & *((ZBV(JLON, JLEV + 1) + ZBV(JLON, JLEV)) / 2.) &
              & *((ZBV(JLON, JLEV + 1) + ZBV(JLON, JLEV)) / 2.)

        DO JW = 1, INW
          
          ZOM(JW, JLON) = ZOP(JW, JLON)
          ZWWM(JW, JLON) = ZWWP(JW, JLON)
        
          ! Intrinsic Frequency
          ZOP(JW, JLON) = ZO(JW, JLON) - ZK(JW, JLON) * ZPCOS(JW, JLON) * ZUH(JLON, JLEV + 1) &
           & - ZK(JW, JLON) * ZPSIN(JW, JLON) * ZVH(JLON, JLEV + 1) 

          ! No breaking
          ! Dissipation
          ! Critical levels (forced to zero if intrinsic frequency changes sign)
          ! Saturation         

          ZOM4 = MAX(ABS(ZOP(JW, JLON) + ZOM(JW, JLON)) / 2.0_JPRB, ZOISEC)   &
                & *MAX(ABS(ZOP(JW, JLON) + ZOM(JW, JLON)) / 2.0_JPRB, ZOISEC) &
                & *MAX(ABS(ZOP(JW, JLON) + ZOM(JW, JLON)) / 2.0_JPRB, ZOISEC) &
                & *MAX(ABS(ZOP(JW, JLON) + ZOM(JW, JLON)) / 2.0_JPRB, ZOISEC)

          ZOP3 = ABS(ZOP(JW, JLON)) &
               &*ABS(ZOP(JW, JLON)) &
               &*ABS(ZOP(JW, JLON))

          ZWWP(JW, JLON) = MIN(ZWWM(JW, JLON)                                     &
           & * EXP(- 2.0_JPRB * NORGWD_RDISS * ZPR / (ZPH(JLON, JLEV + 1) + ZPH(JLON, JLEV))       &
           & * ZBV3                      &
           & / ZOM4          &
           &  * (ZK(JW, JLON)*ZK(JW, JLON)*ZK(JW, JLON)) * (ZH(JLON, JLEV + 1) - ZH(JLON, JLEV))),           &
           & MAX(0.0_JPRB, SIGN(1.0_JPRB, ZOP(JW, JLON) * ZOM(JW, JLON)))                  &
           & * ZOP3 / ZBV(JLON, JLEV + 1)                         &
           & * EXP(- ZH(JLON, JLEV + 1) / ZH0) * NORGWD_SAT*NORGWD_SAT*NORGWD_KMIN*NORGWD_KMIN / &
           & ( ZK(JW, JLON)*ZK(JW, JLON)*ZK(JW, JLON)*ZK(JW, JLON) ))

        ENDDO
       
        DO JW = 1, INW
          ZRUWP(JW, JLON) = SIGN(1.0_JPRB, ZOP(JW, JLON))*ZPCOS(JW, JLON) * ZWWP(JW, JLON)
          ZRVWP(JW, JLON) = SIGN(1.0_JPRB, ZOP(JW, JLON))*ZPSIN(JW, JLON) * ZWWP(JW, JLON)
        ENDDO
     
        ZRUW(JLON, JLEV + 1) = 0.
        ZRVW(JLON, JLEV + 1) = 0.
        ZABSRUW(JLON, JLEV + 1) = 0.

        DO JW = 1, INW
          ZRUW(JLON, JLEV + 1) = ZRUW(JLON, JLEV + 1) + ZRUWP(JW, JLON) 
          ZRVW(JLON, JLEV + 1) = ZRVW(JLON, JLEV + 1) + ZRVWP(JW, JLON) 
          ZABSRUW(JLON, JLEV + 1) = ZABSRUW(JLON, JLEV + 1) + ZWWP(JW, JLON)
        ENDDO

        ZABSRUW(JLON, JLEV + 1) = ZABSRUW(JLON, JLEV + 1) / INW

      ENDDO
    ENDDO

    !-------------------------------------------------------------
    ! 5 CALCUL DES TENDANCES
    !-------------------------------------------------------------

    ! 5.1 Rectification des fluxs au sommet et dans les basses couches:

    DO JLON = KIDIA, KFDIA
      ZRUW(JLON, KLEV + 1) = 0.
      ZRVW(JLON, KLEV + 1) = 0.
      ZABSRUW(JLON, KLEV + 1) = 0.
      ZRUW(JLON, 1) = ZRUW(JLON, ILAUNCH)
      ZRVW(JLON, 1) = ZRVW(JLON, ILAUNCH)
      ZABSRUW(JLON, 1) = ZABSRUW(JLON, ILAUNCH)
    ENDDO
    
    DO JLON = KIDIA, KFDIA
      DO JLEV = 1, ILAUNCH
        ZRUW(JLON, JLEV) = ZRUW(JLON, ILAUNCH+1)
        ZRVW(JLON, JLEV) = ZRVW(JLON, ILAUNCH+1)
        ZABSRUW(JLON, JLEV) = ZABSRUW(JLON, ILAUNCH+1)
      ENDDO
    ENDDO

    DO JLON = KIDIA, KFDIA
      DO JLEV = 1, KLEV
        ZZD_RUW(JLON, JLEV) = 0.5*(ZRUW(JLON, JLEV) + ZRUW(JLON, JLEV+1))
        ZZD_ABSRUW(JLON, JLEV) = 0.5*(ZABSRUW(JLON, JLEV) + ZABSRUW(JLON, JLEV+1))
      ENDDO
      DO JLEV = 1, KLEV
        ZD_RUW(JLON, JLEV)    = ZZD_RUW(JLON, KLEV-JLEV+1)
        ZD_ABSRUW(JLON, JLEV) = ZZD_ABSRUW(JLON, KLEV-JLEV+1)
      ENDDO

    ENDDO

    ! AR-1 RECURSIVE FORMULA (13)
    DO JLON = KIDIA, KFDIA
       DO JLEV = 1, KLEV
          PTEND_U(JLON, JLEV) = (1.-PDTIME/NORGWD_DELTAT) * PTEND_U(JLON, JLEV) + PDTIME/NORGWD_DELTAT/REAL(INW) * &
         &      RG * (ZRUW(JLON, JLEV + 1) - ZRUW(JLON, JLEV)) &
         &      / (ZPH(JLON, JLEV + 1) - ZPH(JLON, JLEV)) * PDTIME
          PTEND_V(JLON, JLEV) = (1.-PDTIME/NORGWD_DELTAT) * PTEND_V(JLON, JLEV) + PDTIME/NORGWD_DELTAT/REAL(INW) * &
         &      RG * (ZRVW(JLON, JLEV + 1) - ZRVW(JLON, JLEV)) &
         &      / (ZPH(JLON, JLEV + 1) - ZPH(JLON, JLEV)) * PDTIME
       ENDDO
    ENDDO

!#ifdef WXIOS
!    IF (LXIOS) THEN
!      CALL MSE_XIOS_SEND('nogwduflx_rand', ZZRUW0_RAND(KIDIA:KFDIA))
!      CALL MSE_XIOS_SEND('nogwduflx_prec', ZZRUW0_PREC(KIDIA:KFDIA))
!      CALL MSE_XIOS_SEND('nogwduflx_back', ZZRUW0_BACK(KIDIA:KFDIA))
!      CALL MSE_XIOS_SEND('nogwduflx_fron', ZZRUW0_FRON(KIDIA:KFDIA))
!      CALL MSE_XIOS_SEND('nogwduflx_net', PFIELD2=ZD_RUW(KIDIA:KFDIA,:))
!      CALL MSE_XIOS_SEND('nogwduflx_abs', PFIELD2=ZD_ABSRUW(KIDIA:KFDIA,:))
!    ENDIF
!#endif

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACNORGWD',1,ZHOOK_HANDLE)

END SUBROUTINE ACNORGWD
