SUBROUTINE AROINI_SURFC (KSV, KSWB, PCO2, PMU0, PDIR_ALB, PSCA_ALB, PEMIS, &
                       & PTSRAD, PSW_BANDS, OMCC03, PITM, YDSURF_ATM_TURB, LDOUT)
USE PARKIND1, ONLY : JPRB, JPIM
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     #######################################################
!
!!****  *INI_SURF* - routine to initialize externalised surface things for AROME
!!                   SURFEX is initialized here
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      S.Malardel   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/2004
!!      A.Alias     06/2007 - setup for ARPEGE (Gaussian Grid)
!!      A.Alias     09/2007 - arguments added to AROINI_SURF (KULOUT/PLOCEN/PSTRET)
!!      Y.Seity     01/2008 - add consistency checks for date
!!      Y.Seity     09/2008 - modifications for DUSTS
!!      A.Alias     04/2008 - init surface file is from now in FA format for a ARPEGE run
!!      Y.Seity     10/2009 - Add missed deallocations
!!     R. El Khatib 06/2010 - CDSURFEX_FNAME
!!     P.Marguinaud 07/2010 - Read only for NBLOCK==1 and distribute on other procs
!!     P.Marguinaud 08/2010 - Cleaning
!!     P.Marguinaud 09/2010 - Broadcast small parameters
!!     A.Voldoire   12/2010 - test on PLOCEN fixed and CALL FLUSH added
!!     A.Alias      12/2010 - Modification to deal with FA file for ALADIN
!!     B.Decharme   12/2010 - Put Earth System Model key into SURFEX
!!     A.Alias      06/2011 - Bugfix : no open of file .des when LLFA
!!     P.Marguinaud 07/2011 - Use YSURFEX_CACHE_IN
!!     P.Marguinaud 01/2012 - Create from aroini_surf.F90
!!     P.Marguinaud 04/2012 - Disable call to SURFEX cache dump
!!     P.Marguinaud 09/2012 - More initialization arguments + make them optional
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!
USE MODD_SURFEX_ARO, ONLY : YSURFEX_ARO_ALL, YSURFEX_ARO_CUR
USE MODD_SURF_PAR,   ONLY   : XUNDEF, NUNDEF
!
USE MODD_IO_SURF_ARO, ONLY : NGPTOT, NPROMA, NBLOCKTOT, NGPTOT_CAP, &
                           & SURFEX_FIELD_BUF_DEALLOC,              &
                           & YSURFEX_CACHE_IN,                      &
                           & NDATE_IN, XTIME_IN,          MYPROC
USE MODD_FROMMPA,  ONLY   : CSV_MSE, NSV_CHEMBEG_MSE, NSV_CHEMEND_MSE, &
                            NSV_AERBEG_MSE, NSV_AEREND_MSE,            &
                            NSV_SLTBEG_MSE, NSV_SLTEND_MSE,            &
                            NSV_DSTBEG_MSE, NSV_DSTEND_MSE,            &
                            NSV_DSTDEPBEG_MSE, NSV_DSTDEPEND_MSE,      &
                            NSV_CO2_MSE
#ifdef USE_SODA
USE MODD_ASSIM,      ONLY : LASSIM,CASSIM_ISBA,XF,YF_PATCH,&
                          & INCV,NVAR,NOBSTYPE,XAT2M_ISBA,&
                          & XAHU2M_ISBA,XVAR,XOBS
#endif
USE MODD_SURF_ATM_TURB_n, ONLY : SURF_ATM_TURB_t
USE MODD_SURF_ATM, ONLY : LDRAG_COEF_ARP
USE MODI_GET_FRAC_N
USE MODI_READ_ALL_NAMELISTS
USE MODI_INI_CSTS
USE MODI_INIT_SURF_ATM_N
USE MODD_IO_SURF_FA,  ONLY : NUNIT_FA
!
USE MODI_INI_SUN
USE MODI_INI_SW_SETUP
USE MODI_ABOR1_SFX
USE MODD_PREP_SNOW, ONLY : NIMPUR
USE MODD_TYPE_DATE_SURF, ONLY : DATE

!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER (KIND=JPIM),    INTENT (IN),  OPTIONAL         :: KSV                               ! number of passiv scalars for chemical scheme
INTEGER (KIND=JPIM),    INTENT (IN),  OPTIONAL         :: KSWB                              ! number of SW bands
REAL (KIND=JPRB),       INTENT (IN),  OPTIONAL         :: PCO2                              ! CO2 concentration (kg/kg)
REAL (KIND=JPRB),       INTENT (IN),  OPTIONAL, TARGET :: PMU0      (:)   !   (NGPTOT)      ! cos of zenithal angle at t  +dt
REAL (KIND=JPRB),       INTENT (OUT), OPTIONAL         :: PDIR_ALB  (:,:) !   (NGPTOT,KSWB) ! direct albedo for each band
REAL (KIND=JPRB),       INTENT (OUT), OPTIONAL         :: PSCA_ALB  (:,:) !   (NGPTOT,KSWB) ! diffuse albedo for each band
REAL (KIND=JPRB),       INTENT (OUT), OPTIONAL         :: PEMIS     (:)   !   (NGPTOT)      ! emissivity
REAL (KIND=JPRB),       INTENT (OUT), OPTIONAL         :: PTSRAD    (:)   !   (NGPTOT)      ! radiative temperature
REAL (KIND=JPRB),       INTENT (IN),  OPTIONAL, TARGET :: PSW_BANDS (:)   !   (KSWB)        ! centers of spectral bands
REAL (KIND=JPRB),       INTENT (OUT), OPTIONAL         :: PITM      (:)   !   (NGPTOT)      ! land-sea mask
LOGICAL,                INTENT (IN),  OPTIONAL         :: OMCC03                            ! Oceanic variables are given by the coupler
TYPE (SURF_ATM_TURB_t), INTENT (IN),  OPTIONAL         :: YDSURF_ATM_TURB                   ! Atmospheric turbulence parameters
LOGICAL,                INTENT (IN),  OPTIONAL         :: LDOUT
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER            :: IGPCOMP        ! DIM of array excepted upper E aladin zone

REAL, DIMENSION(NGPTOT) :: ZCO2
REAL, DIMENSION(NGPTOT) :: ZRHODREF
REAL, DIMENSION(NGPTOT) :: ZZENITH
REAL, DIMENSION(NGPTOT) :: ZAZIM
REAL, DIMENSION(NGPTOT) :: ZSEA, ZTOWN, ZWATER, ZNATURE
CHARACTER(LEN=6), DIMENSION(:) , ALLOCATABLE  :: YSV
INTEGER :: JLAYER
LOGICAL :: LLOUT,LLAND_USE
INTEGER :: IREP
!
INTEGER       :: IBLOCK
INTEGER       :: ISV
INTEGER       :: ISWB
REAL, POINTER :: ZSW_BANDS (:)
REAL, TARGET  :: ZSW_BANDS0 (1)
REAL, POINTER :: ZMU0 (:)
REAL, TARGET  :: ZMU00 (NGPTOT)
REAL          :: ZCO20
INTEGER,SAVE  :: ZIVAR=0 
INTEGER       :: IRESP
TYPE(SURF_ATM_TURB_t) :: YLSURF_ATM_TURB
REAL, DIMENSION(NGPTOT,NIMPUR)  :: ZIMPWET       
REAL, DIMENSION(NGPTOT,NIMPUR)  :: ZIMPDRY    
TYPE(DATE) :: TDATE_END
!
!-------------------------------------------------------------------------------

!
! ** IO initialisation
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AROINI_SURFC',0,ZHOOK_HANDLE)

LLOUT = .TRUE.
IF (PRESENT (LDOUT)) LLOUT = LDOUT

ISV = 0
IF (PRESENT (KSV)) ISV = KSV

IF (PRESENT (KSWB)) THEN
  ISWB = KSWB
  ZSW_BANDS => PSW_BANDS
ELSE
  ZSW_BANDS0 (1) = 999.
  ISWB = 1
  ZSW_BANDS => ZSW_BANDS0
ENDIF

IF (PRESENT (PMU0)) THEN
  ZMU0 => PMU0
ELSE
  ZMU00 = 1.
  ZMU0 => ZMU00
ENDIF

IF (PRESENT (PCO2)) THEN
  ZCO20 = PCO2
ELSE
  ZCO20 = 0.
ENDIF

IGPCOMP=MIN(NGPTOT,NGPTOT_CAP)

ALLOCATE(YSV(ISV))

!
! ZRHODREF used only for chemistry, not used in AROME prototype
ZRHODREF(:)=1.16
!
ZCO2(:) = ZCO20
!
!
! Astronomic parameter initialisation
!
ZZENITH (:) = ACOS( ZMU0 (:) )
!
ZAZIM   (:) = 0. !Azimuthal angle is not used in SURFEX
!
!YLSURF_ATM_TURB%LPTKE=.FALSE.      ! default: don't use Pseudo-TKE
YLSURF_ATM_TURB%LCOEFKTKE=.FALSE.  ! default: don't use TOUCANS turbulence
IF (PRESENT(YDSURF_ATM_TURB)) YLSURF_ATM_TURB=YDSURF_ATM_TURB  ! atmospheric turbulence parameters

! Externalised Surface Initialisation
!
! fixed options : SIZE(CSV_MSE,1) passive scalar
!                 6 short-wave spectral bands
!
IF (ISV > 0) THEN
  YSV(:) = "     "
  DO JLAYER=1,SIZE(CSV_MSE,1)
    YSV(JLAYER)(1:6) = CSV_MSE(JLAYER)(1:6)
  END DO
  DO JLAYER=NSV_DSTBEG_MSE,NSV_DSTEND_MSE
    YSV(JLAYER) = TRIM(CSV_MSE(JLAYER))
  END DO
  DO JLAYER=NSV_DSTDEPBEG_MSE,NSV_DSTDEPEND_MSE
    YSV(JLAYER) = TRIM(CSV_MSE(JLAYER))
  END DO
  DO JLAYER=NSV_CHEMBEG_MSE,NSV_CHEMEND_MSE
    YSV(JLAYER) = '#'//TRIM(CSV_MSE(JLAYER))
  END DO
  DO JLAYER=NSV_AERBEG_MSE,NSV_AEREND_MSE
    YSV(JLAYER) = '@'//TRIM(CSV_MSE(JLAYER))
  END DO
  IF (NSV_CO2_MSE > 0) YSV(NSV_CO2_MSE) = TRIM(CSV_MSE(NSV_CO2_MSE))
ENDIF

!
! Model number initialisation (or later NPROMA bloc number)
!
IF (LLOUT) THEN
  IF (PRESENT (PDIR_ALB)) PDIR_ALB = XUNDEF
  IF (PRESENT (PSCA_ALB)) PSCA_ALB = XUNDEF
  IF (PRESENT (PEMIS   )) PEMIS    = XUNDEF
  IF (PRESENT (PTSRAD  )) PTSRAD   = XUNDEF
  IF (PRESENT (PITM    )) PITM     = 0.
ENDIF
!init LLAND_USE
LLAND_USE=.FALSE.
!
ZSEA    = 0.
ZNATURE = 0.
ZWATER  = 0.
ZTOWN   = 0.

IF (NBLOCKTOT >= 1) THEN
  CALL PROCESS_SURFEX_BLOCK (1)
ENDIF


DO IBLOCK = 2, NBLOCKTOT
  CALL PROCESS_SURFEX_BLOCK (IBLOCK)
END DO


IF (LLOUT .AND. PRESENT (PITM)) &
& PITM(:) = ZNATURE(:) + ZTOWN(:)

DEALLOCATE(YSV)

IF (NBLOCKTOT >= 1) THEN
  YSURFEX_ARO_CUR => YSURFEX_ARO_ALL (1)
ENDIF

IF (LHOOK) CALL DR_HOOK('AROINI_SURFC',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE PROCESS_SURFEX_BLOCK(IBLOCK)
USE MODD_IO_SURF_ARO, ONLY :  NINDX1, NINDX2, NKPROMA, NBLOCK
USE MODE_MODELN_SURFEX_HANDLER, ONLY : ICURRENT_MODEL
INTEGER, INTENT(IN) :: IBLOCK
!
REAL, DIMENSION(NPROMA,ISWB) :: ZW1, ZW2
REAL, DIMENSION(NPROMA)      :: ZW3, ZW4, ZW5
INTEGER                           :: IIVAR
INTEGER                           :: IOBS

  ZW1 = 0._JPRB
  ZW2 = 0._JPRB
  ZW3 = 0._JPRB
  ZW4 = 0._JPRB
  ZW5 = 0._JPRB

  ICURRENT_MODEL = NBLOCK

  NBLOCK   = IBLOCK
  NINDX1   = 1+(NBLOCK-1)*NPROMA
  NINDX2   = MIN(NBLOCK*NPROMA,IGPCOMP)
  NKPROMA  = NINDX2-NINDX1+1

  YSURFEX_ARO_CUR => YSURFEX_ARO_ALL (NBLOCK)
  TDATE_END%YEAR = NDATE_IN(1)
  TDATE_END%MONTH = NDATE_IN(2)
  TDATE_END%DAY = NDATE_IN(3)


  CALL INIT_SURF_ATM_N (YSURFEX_ARO_CUR, 'AROME ', 'ALL', LLAND_USE, NKPROMA, ISV, ISWB,            &
              & YSV, ZCO2(NINDX1:NINDX2), ZIMPWET(NINDX1:NINDX2,:),&
              &  ZIMPDRY(NINDX1:NINDX2,:), ZRHODREF(NINDX1:NINDX2), &
              & ZZENITH(NINDX1:NINDX2), ZAZIM(NINDX1:NINDX2), ZSW_BANDS,&
              & ZW1(1:NKPROMA,1:ISWB), ZW2(1:NKPROMA,1:ISWB),&
              & ZW3(1:NKPROMA), ZW4(1:NKPROMA), ZW5(1:NKPROMA),&
              & NDATE_IN(1), NDATE_IN(2), NDATE_IN(3), XTIME_IN,TDATE_END,   &
              & YLSURF_ATM_TURB, &
              & '                            ','      ',  'OK')

  IF (LLOUT) THEN
    DO JLAYER =1, ISWB
      IF (PRESENT (PDIR_ALB)) PDIR_ALB(NINDX1:NINDX2,JLAYER) = ZW1(1:NKPROMA,JLAYER)
      IF (PRESENT (PSCA_ALB)) PSCA_ALB(NINDX1:NINDX2,JLAYER) = ZW2(1:NKPROMA,JLAYER)
      IF (PRESENT (PEMIS   )) PEMIS (NINDX1:NINDX2) = ZW3(1:NKPROMA)
      IF (PRESENT (PTSRAD  )) PTSRAD(NINDX1:NINDX2) = ZW4(1:NKPROMA)
    END DO

  ENDIF

  CALL GET_FRAC_N(YSURFEX_ARO_CUR%U, 'AROME ',NKPROMA,ZSEA   (NINDX1:NINDX2),ZWATER(NINDX1:NINDX2), &
                  ZNATURE(NINDX1:NINDX2),ZTOWN (NINDX1:NINDX2)  )

#ifdef USE_SODA
  ! Assign the control/perturbed state vectors if running EKF
  IF ( LASSIM .AND. CASSIM_ISBA == "EKF" ) THEN
    IF (SUM(INCV) /= NVAR) THEN
      WRITE(*,*) 'INCONSISTENCY in set-up of CONTROL VARIABLES',SUM(INCV),NVAR
      CALL ABOR1_SFX('INCONSISTENCY in set-up of CONTROL VARIABLES')
    ENDIF

    IF ( ZIVAR == 0 ) THEN   
      ALLOCATE(XF(YSURFEX_ARO_CUR%U%NSIZE_NATURE,YSURFEX_ARO_CUR%IM%I%NPATCH,NVAR+1,NVAR))
      ALLOCATE(YF_PATCH(YSURFEX_ARO_CUR%U%NSIZE_NATURE,YSURFEX_ARO_CUR%IM%I%NPATCH,NVAR+1,NOBSTYPE))
    ENDIF
    IF (( ZIVAR > 0 ) .AND. (ZIVAR <= (NVAR + 1) ))  THEN
      ! Set the global state values for this control value
      ! The control run has index ZIVAR=1
      DO IOBS = 1,NOBSTYPE
        SELECT CASE (TRIM(XOBS(IOBS)))
          CASE("T2M")
            YF_PATCH(:,:,ZIVAR,IOBS) = XAT2M_ISBA(:,:)
          CASE("HU2M")
            YF_PATCH(:,:,ZIVAR,IOBS) = XAHU2M_ISBA(:,:)
          CASE("WG1")
            YF_PATCH(:,:,ZIVAR,IOBS) = YSURFEX_ARO_CUR%IM%I%XWG(:,1,:)
          CASE DEFAULT
            CALL ABOR1_SFX("Mapping of "//XOBS(IOBS)//" is not defined in AROINI_SURFC!")
        END SELECT
      ENDDO

      ! Prognostic fields for assimilation (Control vector)
      DO IIVAR = 1,NVAR
        SELECT CASE (TRIM(XVAR(IIVAR)))
          CASE("TG1")
            XF(:,:,ZIVAR,IIVAR) = YSURFEX_ARO_CUR%IM%I%XTG(:,1,:)
          CASE("TG2")
            XF(:,:,ZIVAR,IIVAR) = YSURFEX_ARO_CUR%IM%I%XTG(:,2,:)
          CASE("WG1")
            XF(:,:,ZIVAR,IIVAR) = YSURFEX_ARO_CUR%IM%I%XWG(:,1,:)
          CASE("WG2")
            XF(:,:,ZIVAR,IIVAR) = YSURFEX_ARO_CUR%IM%I%XWG(:,2,:)
          CASE DEFAULT
            CALL ABOR1_SFX("Mapping of "//TRIM(XVAR(IIVAR))//" is not defined in AROINI_SURFC!")
        END SELECT
      ENDDO
    ENDIF
    ZIVAR=ZIVAR+1
  ENDIF  
#endif

END SUBROUTINE PROCESS_SURFEX_BLOCK

END SUBROUTINE AROINI_SURFC

