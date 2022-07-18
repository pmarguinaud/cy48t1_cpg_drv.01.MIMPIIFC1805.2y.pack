MODULE TURBHALO_MOD

! 3D TURBULENCE SCHEME SUPPORT : COMPUTATION OF HORIZONTAL GRADIENTS.

! GPTURB3D     : Fill a KPROMA-sized core array with fields to differentiate
! EGPTURB3DLAG : Compute gradients from the NGPTOT-sized core+halo array in the LAM model.
!                 which has a regular grid
! GEOTURB :  Compute the geopotential and horinzontal wind for altitude correction 
! 
!     Author.
!     -------
!        Rachel Honnert & Ryad El Khatib *Meteo-France*  28-Nov-2019

! Modifications
! -------------

!-----------------------------------------------------------------------------

USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK

IMPLICIT NONE

SAVE

CONTAINS

SUBROUTINE GPTURB3D(KST,KEND,KPROMA,KFLDCORE,KFLEVG,KDIMGMV,YDT0,PGMV,PRTL,PRTM,PCORE,PRPARITY,PPHIH,PPHIFL,PPHIFM,PUVH,PRT)

USE TYPE_GMVS, ONLY : TYPE_T0
USE INTDYN_MOD,ONLY : YYTHW0
USE YOMDYNA  , ONLY : YRDYNA
!USE YOMCT0       , ONLY : LSPRT
!USE YOMCST       , ONLY : RD, RV

IMPLICIT NONE

! YDT0 : fields pointers in PGMV
! PRPARITY : parity array of output fields (for the global model only)

INTEGER(KIND=JPIM), INTENT(IN)  :: KST
INTEGER(KIND=JPIM), INTENT(IN)  :: KEND
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM), INTENT(IN)  :: KFLDCORE
INTEGER(KIND=JPIM), INTENT(IN)  :: KFLEVG
INTEGER(KIND=JPIM), INTENT(IN)  :: KDIMGMV
TYPE(TYPE_T0),      INTENT(IN)  :: YDT0
REAL(KIND=JPRB),    INTENT(IN)  :: PGMV(KPROMA,KFLEVG,KDIMGMV)
REAL(KIND=JPRB),    INTENT(IN)  :: PRTL(KPROMA,KFLEVG) 
REAL(KIND=JPRB),    INTENT(IN)  :: PRTM(KPROMA,KFLEVG) 
REAL(KIND=JPRB),    INTENT(OUT) :: PRPARITY(KFLDCORE*KFLEVG)
REAL(KIND=JPRB),    INTENT(OUT) :: PCORE(KPROMA,KFLDCORE*KFLEVG)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)  :: PPHIH(KPROMA,0:KFLEVG) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)  :: PPHIFL(KPROMA,KFLEVG) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)  :: PPHIFM(KPROMA,KFLEVG) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)  :: PUVH(KPROMA,0:KFLEVG,YYTHW0%NDIM) 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)  :: PRT(KPROMA,KFLEVG) 
REAL(KIND=JPRB) :: ZUL, ZVL, ZVOR, ZDIV
REAL(KIND=JPRB) :: ZTEMP_DUDX, ZTEMP_DVDY, ZTEMP_DTDX, ZTEMP_DTDY, ZTEMP_DRDX, ZTEMP_DRDY ! horinzontal gradients 
INTEGER(KIND=JPIM) :: JLEV, JROF
REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZDELTA_PHIF
REAL(KIND=JPRB) :: ZDELTA_MUF
REAL(KIND=JPRB) :: ZDELTA_MVF
REAL(KIND=JPRB) :: ZDELTA_MTF
REAL(KIND=JPRB) :: ZDELTA_MRF
REAL(KIND=JPRB) :: ZEPS

!-----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TURBHALO_MOD:GPTURB3D',0,ZHOOK_HANDLE)

! For now:
PRPARITY(:)=HUGE(1._JPRB)
!
IF(PRESENT(PPHIH))THEN
 ZEPS=1.E-12
 DO JLEV=1,KFLEVG
  DO JROF=KST,KEND
    ZDELTA_PHIF=PPHIH(JROF,JLEV)-PPHIH(JROF,JLEV-1)
    ZDELTA_MUF=PUVH(JROF,JLEV,YYTHW0%M_UH)-PUVH(JROF,JLEV-1,YYTHW0%M_UH)
    ZDELTA_MVF=PUVH(JROF,JLEV,YYTHW0%M_VH)-PUVH(JROF,JLEV-1,YYTHW0%M_VH)
    !rek ZDELTA_MTF=PGMV(JROF,JLEV,YDT0%MT)-PGMV(JROF,JLEV-1,YDT0%MT)
    !rek ZDELTA_MRF=PRT(JROF,JLEV)-PRT(JROF,JLEV-1)
    ZDELTA_MTF=PGMV(JROF,JLEV,YDT0%MT)-PGMV(JROF,MAX(1,JLEV-1),YDT0%MT)
    ZDELTA_MRF=PRT(JROF,JLEV)-PRT(JROF,MAX(1,JLEV-1))
    IF (ZDELTA_PHIF > ZEPS ) CALL ABOR1("TURBHALO.F90 : ABORT d(gz)/dz > 0")
    ZUL =PGMV(JROF,JLEV,YDT0%MUL)
    ZVL =PGMV(JROF,JLEV,YDT0%MVL)
    ZDIV=PGMV(JROF,JLEV,YDT0%MDIV)
    ZVOR=PGMV(JROF,JLEV,YDT0%MVOR)
    ! du/dx_z
    ZTEMP_DUDX=ZUL       -PPHIFL(JROF,JLEV)*ZDELTA_MUF/ZDELTA_PHIF
    ! du/dy_z
    !ZTEMP_DUDY=ZVL-ZVOR -PPHIFM(JROF,JLEV)*ZDELTA_MUF/ZDELTA_PHIF
    ! dv/dx_z
    !ZTEMP_DVDX=ZVL      -PPHIFL(JROF,JLEV)*ZDELTA_MVF/ZDELTA_PHIF
    ! dv/dy_z
    ZTEMP_DVDY=ZDIV-ZUL  -PPHIFM(JROF,JLEV)*ZDELTA_MVF/ZDELTA_PHIF
    ZTEMP_DTDX=PGMV(JROF,JLEV,YDT0%MTL) -PPHIFM(JROF,JLEV)*ZDELTA_MTF/ZDELTA_PHIF
    ZTEMP_DTDY=PGMV(JROF,JLEV,YDT0%MTM) -PPHIFM(JROF,JLEV)*ZDELTA_MTF/ZDELTA_PHIF
    ZTEMP_DRDX=PRTL(JROF,JLEV) -PPHIFM(JROF,JLEV)*ZDELTA_MRF/ZDELTA_PHIF
    ZTEMP_DRDY=PRTM(JROF,JLEV) -PPHIFM(JROF,JLEV)*ZDELTA_MRF/ZDELTA_PHIF
    ! *         Fill PCORE
    PCORE(JROF,JLEV)=ZTEMP_DUDX
    PCORE(JROF,JLEV+KFLEVG)=ZTEMP_DVDY
    PCORE(JROF,JLEV+2*KFLEVG)=ZTEMP_DTDX
    PCORE(JROF,JLEV+3*KFLEVG)=ZTEMP_DTDY
    PCORE(JROF,JLEV+4*KFLEVG)=ZTEMP_DRDX
    PCORE(JROF,JLEV+5*KFLEVG)=ZTEMP_DRDY
  ENDDO
 ENDDO
ELSE
 DO JLEV=1,KFLEVG
  DO JROF=KST,KEND
    ZUL =PGMV(JROF,JLEV,YDT0%MUL)
    ZVL =PGMV(JROF,JLEV,YDT0%MVL)
    ZDIV=PGMV(JROF,JLEV,YDT0%MDIV)
    ZVOR=PGMV(JROF,JLEV,YDT0%MVOR)
    ! du/dx_z
    ZTEMP_DUDX=ZUL       
    ! du/dy_z
    !ZTEMP_DUDY=ZVL-ZVOR
    ! dv/dx_z
    !ZTEMP_DVDX=ZVL     
    ! dv/dy_z
    ZTEMP_DVDY=ZDIV-ZUL 
    ZTEMP_DTDX=PGMV(JROF,JLEV,YDT0%MTL)
    ZTEMP_DTDY=PGMV(JROF,JLEV,YDT0%MTM)
    ZTEMP_DRDX=PRTL(JROF,JLEV)
    ZTEMP_DRDY=PRTM(JROF,JLEV)
    ! *         Fill PCORE
    PCORE(JROF,JLEV)=ZTEMP_DUDX
    PCORE(JROF,JLEV+KFLEVG)=ZTEMP_DVDY
    PCORE(JROF,JLEV+2*KFLEVG)=ZTEMP_DTDX
    PCORE(JROF,JLEV+3*KFLEVG)=ZTEMP_DTDY
    PCORE(JROF,JLEV+4*KFLEVG)=ZTEMP_DRDX
    PCORE(JROF,JLEV+5*KFLEVG)=ZTEMP_DRDY
  ENDDO
 ENDDO
ENDIF

IF (LHOOK) CALL DR_HOOK('TURBHALO_MOD:GPTURB3D',1,ZHOOK_HANDLE)

END SUBROUTINE GPTURB3D

SUBROUTINE EGPTURB3DLAG(KST,KEND,KPROMA,KFLDCORE,KFLEVG,KASLB1,KL0,PHALO,YDEGEO,KGRADIENTS,PGRADIENT)

! KL0 : adresses of the 4 western points "13", "7", "9", "15" in the halo,
!       while point "1" is the originating model gridpoint :

!                     13       5       6      14

!                      7       1       2       8
!                                 
!                      9       3       4      10

!                     15      11      12      16

USE YEMGEO   , ONLY : TEGEO

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)  :: KST
INTEGER(KIND=JPIM), INTENT(IN)  :: KEND
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROMA
INTEGER(KIND=JPIM), INTENT(IN)  :: KFLDCORE
INTEGER(KIND=JPIM), INTENT(IN)  :: KFLEVG
INTEGER(KIND=JPIM), INTENT(IN)  :: KASLB1
INTEGER(KIND=JPIM), INTENT(IN)  :: KL0(KPROMA,4)
REAL(KIND=JPRB),    INTENT(IN)  :: PHALO(KASLB1,KFLDCORE*KFLEVG)
TYPE(TEGEO),        INTENT(IN)  :: YDEGEO
INTEGER(KIND=JPIM), INTENT(IN)  :: KGRADIENTS
REAL(KIND=JPRB),    INTENT(OUT) :: PGRADIENT(KPROMA,KGRADIENTS,KFLEVG)

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TURBHALO_MOD:EGPTURB3DLAG',0,ZHOOK_HANDLE)

ASSOCIATE(EDELY=>YDEGEO%EDELY,EDELX=>YDEGEO%EDELX)

DO JLEV=1,KFLEVG
  DO JROF=KST,KEND
   ! du/dx_z
   PGRADIENT(JROF,1,JLEV)= PHALO(KL0(JROF,2)+1,JLEV)  
   ! dv/dy_z
   PGRADIENT(JROF,2,JLEV)= PHALO(KL0(JROF,2)+1,JLEV+KFLEVG)
   PGRADIENT(JROF,3,JLEV)= PHALO(KL0(JROF,2)+1,JLEV+2*KFLEVG) ! dt/dx
   PGRADIENT(JROF,4,JLEV)= PHALO(KL0(JROF,2)+1,JLEV+3*KFLEVG) ! dt/dy
   PGRADIENT(JROF,5,JLEV)= PHALO(KL0(JROF,2)+1,JLEV+4*KFLEVG) ! drv/dx
   PGRADIENT(JROF,6,JLEV)= PHALO(KL0(JROF,2)+1,JLEV+5*KFLEVG) ! drv/dy
   ! d²u/dx²_z
   PGRADIENT(JROF,7,JLEV)=&
         & (PHALO(KL0(JROF,2)+2,JLEV)- &
         &  PHALO(KL0(JROF,2),JLEV))/(2._JPRB*EDELX)
   ! d²v/dy²_z
   PGRADIENT(JROF,8,JLEV)=&
         & (PHALO(KL0(JROF,1)+1,KFLEVG+JLEV)-&
         &  PHALO(KL0(JROF,3)+1,KFLEVG+JLEV))/(2._JPRB*EDELY)
  ENDDO
ENDDO

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('TURBHALO_MOD:EGPTURB3DLAG',1,ZHOOK_HANDLE)

END SUBROUTINE EGPTURB3DLAG

!============================================================================

SUBROUTINE GEOTURB(&
 ! ----- INPUT ---------------------------------------------------------------
 & YDGEOMETRY,KPROMA,KSTART,KEND,POROG,POROGL,POROGM,&
 & PSP,PSPL,PSPM,PT,PTL,PTM,&
 & PQ,PQL,PQM,PL,PLL,PLM,PI,PIL,PIM,PR,PS,PG,PUF,PVF,&
 & PSPD,PSPDL,PSPDM,&
 ! ----- OUTPUT --------------------------------------------------------------
 & PPHIH,PPHIFL,PPHIFM,PUVH,PRT,PRTL,PRTM)

! -----
! REMARKS:
!  - Variables PLL, PLM, PIL, PIM are not used currently, but they
!    might be needed in corrected version of GPRT.
!  - The structure of this routine must remain consistent with the one of CPG_GP.
!  - In the NHQE model, prognostic T is a modified temperature.
! -----

! GEOTURB - Diagnose NHX-term

! Purpose
! -------
!   Diagnose NHX-term

! Interface
! ---------
! * INPUT:
!   YDGEOMETRY   : structure containing all geometry.
!   KPROMA       : length of work
!   KSTART       : start of work
!   KEND         : end of work
!   POROG        : surface geopotential
!   POROGL       : zonal gradient of surface geopotential
!   POROGM       : meridional gradient of surface geopotential
!   PSP          : ln(prehyds)
!   PSPL         : zonal gradient of ln(prehyds)
!   PSPM         : meridional gradient of ln(prehyds)
!   PT           : temperature T
!   PTL          : zonal gradient of temperature
!   PTM          : meridional gradient of temperature
!   PQ           : specific humidity q
!   PQL          : zonal gradient of q
!   PQM          : meridional gradient of q
!   PL           : specific mass of liquid ql
!   PLL          : zonal gradient of ql
!   PLM          : meridional gradient of ql
!   PI           : specific mass of ice qi
!   PIL          : zonal gradient of qi
!   PIM          : meridional gradient of qi
!   PR           : specific mass of rain qr
!   PS           : specific mass of snow qs
!   PG           : specific mass of graupel qg
!   PUF          : U-wind at full levels
!   PVF          : V-wind at full levels
!   PSPD         : 
!   PSPDL        : 
!   PSPDM        : 

! * OUTPUT:

!  PPHIH          : geopotantial at half layer (required for NHEE only)
!  PPHIFL         : zonal gradient of geopotantial (required for NHEE only)
!  PPHIFM         : meridional gradient of geopotantial (required for NHEE only)
!  PUVH           : Wind at half layer

! Externals
! ---------

! Method
! ------

! Reference
! ---------

! Author
! ------
!  02-Dec-2019 Rachel Honnert & Fabrice Voitus 

! Modifications
! -------------


! End Modifications
!---------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMDIMV      , ONLY : TDIMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : RD, RV
USE YOMCT0       , ONLY : LSPRT, LNHEE
!USE YOMCVER      , ONLY : LVERTFE, LVFE_X_TERM
USE INTDYN_MOD   , ONLY :  YYTHW0, YYTXYBDER, YYTXYB
USE YOMDYNA      , ONLY : YRDYNA

! -----------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY),    INTENT(IN)       :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)       :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)       :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)       :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)       :: POROG (KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: POROGL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: POROGM(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PSP   (KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PSPL  (KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PSPM  (KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PT    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PTL   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PTM   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PQ    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PQL   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PQM   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PL    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PLL   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)       :: PLM   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)       :: PI    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PIL   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)       :: PIM   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)       :: PR    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PS    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PG    (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PUF   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PVF   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PSPD  (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PSPDL (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)       :: PSPDM (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)      :: PPHIH (KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)      :: PPHIFL(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)      :: PPHIFM(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)      :: PUVH  (KPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTHW0%NDIM) 
REAL(KIND=JPRB)   ,INTENT(OUT)      :: PRT   (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)  ! RT
REAL(KIND=JPRB)   ,INTENT(OUT)      :: PRTL  (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)  ! zonal gradient of RT
REAL(KIND=JPRB)   ,INTENT(OUT)      :: PRTM  (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)  ! meridional gradient of RT

! -----------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,JLEV

REAL(KIND=JPRB) :: ZHOOK_HANDLE

REAL(KIND=JPRB) :: ZPIS  (KPROMA)          ! surface pressure prehyds
REAL(KIND=JPRB) :: ZPISL (KPROMA)          ! zonal derivative of prehyds 
REAL(KIND=JPRB) :: ZPISM (KPROMA)          ! meridional derivative of prehyds

REAL(KIND=JPRB) :: ZXYBDER(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYBDER%NDIM) ! cf. PXYBDER in GPGRXYB
REAL(KIND=JPRB) :: ZXYB  (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)  ! contains "delta", "alpha"
REAL(KIND=JPRB) :: ZPIH  (KPROMA, 0:YDGEOMETRY%YRDIMV%NFLEVG)! prehyd at half levels
REAL(KIND=JPRB) :: ZPIF  (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! prehyd at full levels
REAL(KIND=JPRB) :: ZCP   (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! specific heat of air cp
REAL(KIND=JPRB) :: ZR    (KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! gas "constant" of air R
REAL(KIND=JPRB) :: ZKAPPA(KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)  ! kappa = R/cp

REAL(KIND=JPRB) :: ZPHIF(KPROMA, YDGEOMETRY%YRDIMV%NFLEVG)     ! merid comp of geopotential (full lay)
REAL(KIND=JPRB) :: ZPHIHL(KPROMA, 0:YDGEOMETRY%YRDIMV%NFLEVG)    ! zonal comp of grad(gz) (half lay)
REAL(KIND=JPRB) :: ZPHIHM(KPROMA, 0:YDGEOMETRY%YRDIMV%NFLEVG)    ! merid comp of grad(gz) (half lay)

REAL(KIND=JPRB) :: ZNHPREF(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! "pre" at full levels.
REAL(KIND=JPRB) :: ZNHPPI (KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! "pre/prehyd" at full levels.
REAL(KIND=JPRB) :: ZRNHPPI(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! "prehyd/pre" at full levels.
REAL(KIND=JPRB) :: ZNHPREL(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! zon comp of "grad pre" (full lev)
REAL(KIND=JPRB) :: ZNHPREM(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! mer comp of "grad pre" (full lev)
REAL(KIND=JPRB) :: ZLNNHPREFL(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) ! zon comp grad(ZNHPREF)/ZNHPREF
REAL(KIND=JPRB) :: ZLNNHPREFM(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG) ! mer comp grad(ZNHPREF)/ZNHPREF
REAL(KIND=JPRB) :: ZQCHAL(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)     ! zonal comp grad(log(pre/prehyd))
REAL(KIND=JPRB) :: ZQCHAM(KPROMA,YDGEOMETRY%YRDIMV%NFLEVG)     ! merid comp grad(log(pre/prehyd))
REAL(KIND=JPRB) :: ZRRED0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
! -----------------------------------------------------------------------------

#include "gphpre.intfb.h"
#include "gpgrxyb.intfb.h"
#include "gprcp_qlirsg.intfb.h"
#include "gprt.intfb.h"
#include "gphlwi.intfb.h"
#include "gphluv.intfb.h"
#include "gnhpre.intfb.h"
!#include "abor1.intfb.h"
#include "gnhgrpre.intfb.h"
#include "gpgeo.intfb.h"
#include "gpgrgeo.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TURBHALO_MOD:GEOTURB',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,YDDIMV=>YDGEOMETRY%YRDIMV)

! -----------------------------------------------------------------------------
! -----
! computation of pressures
! -----
! convert ln(prehyds) to prehyds (same for derivatives)
DO JROF=KSTART,KEND
  ZPIS (JROF)=EXP(PSP(JROF))
  ZPISL(JROF)=ZPIS(JROF)*PSPL(JROF)
  ZPISM(JROF)=ZPIS(JROF)*PSPM(JROF)
  ! half level pressure
  ZPIH(JROF,NFLEVG)=ZPIS(JROF)
ENDDO

CALL GPHPRE(KPROMA,NFLEVG,KSTART,KEND,YDVAB,ZPIH,PXYB=ZXYB,PRESF=ZPIF)

! additional auxiliary quantities (grad(delta) and grad(alpha) for ex.)
CALL GPGRXYB(KPROMA,KSTART,KEND,NFLEVG,.FALSE.,YDVAB,ZPISL,ZPISM,ZXYB,ZXYBDER)

! -----
! computation of half level wind 
! -----

! compute interpolation weights
CALL GPHLWI(YDDIMV,KPROMA,KSTART,KEND,ZXYB(:,:,YYTXYB%M_LNPR),ZXYB(:,:,YYTXYB%M_ALPH),PUVH(:,:,YYTHW0%M_WWI))
! interpolate wind into half levels
CALL GPHLUV(YDDIMV,KPROMA,KSTART,KEND,PUF,PVF,PUVH)

! -----
! computation of R, cp, kappa
! -----

CALL GPRCP_QLIRSG(KPROMA,KSTART,KEND,NFLEVG,PQ=PQ,PQI=PI,PQL=PL,PQR=PR,PQS=PS,PQG=PG,&
 & PCP=ZCP,PR=ZR,PKAP=ZKAPPA)  

! -----
! computation of RT, grad(RT)
! -----

CALL GPRT(LSPRT,KPROMA,KSTART,KEND,NFLEVG,RD,RV,ZR,PT,PTL,&
 & PTM,PQL,PQM,PRT,PRTL,PRTM)

! -----
! computation of pre/prehyd and some other "pressure departure" quantities for NHEE model
! -----

CALL GNHPRE(YRDYNA%NPDVAR,YDGEOMETRY,KPROMA,NFLEVG,KSTART,KEND,PSPD,ZPIF,PNHPREF=ZNHPREF,PNHPPI=ZNHPPI,PRNHPPI=ZRNHPPI)

DO JLEV=1,NFLEVG
 DO JROF=KSTART,KEND
  ZRRED0(JROF,JLEV)=ZR(JROF,JLEV)/ZNHPPI(JROF,JLEV)
 ENDDO
ENDDO

CALL GNHGRPRE(YRDYNA%NPDVAR,YDGEOMETRY,KPROMA,NFLEVG,KSTART,KEND,ZXYB(:,:,YYTXYB%M_RTGR),ZPISL,ZPISM,&
   & ZNHPREF,PSPDL,PSPDM,ZNHPREL,ZNHPREM,ZLNNHPREFL,ZLNNHPREFM,PQCHAL=ZQCHAL,PQCHAM=ZQCHAM)

! -----
! computation of phi=gz (geopotential)
! -----

PPHIH(KSTART:KEND,NFLEVG)=POROG(KSTART:KEND)
CALL GPGEO(KPROMA,KSTART,KEND,NFLEVG,PPHIH,ZPHIF,PT,&
 & ZRRED0,ZXYB(:,:,YYTXYB%M_LNPR),ZXYB(:,:,YYTXYB%M_ALPH),&
 & YDGEOMETRY%YRVERT_GEOM )

! -----
! computation of grad(gz)
! -----

CALL GPGRGEO(YDGEOMETRY,KPROMA,KSTART,KEND,NFLEVG,&
   & PRT,PRTL,PRTM,&
   & ZXYB(:,:,YYTXYB%M_LNPR),ZXYB(:,:,YYTXYB%M_ALPH),ZXYBDER,&
   & POROGL,POROGM,&
   & PPHIFL,PPHIFM,ZPHIHL,ZPHIHM,&
   & LDNHEE=LNHEE,PRNHPPI=ZRNHPPI,PQCHAL=ZQCHAL,PQCHAM=ZQCHAM)
!
!
!! -----------------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TURBHALO_MOD:GEOTURB',1,ZHOOK_HANDLE)

END SUBROUTINE GEOTURB


END MODULE TURBHALO_MOD
