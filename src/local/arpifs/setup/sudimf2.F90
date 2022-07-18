SUBROUTINE SUDIMF2(YDML_GCONF)

!**** *SUDIMF2* - Set up number of fields dimensions part 2.

!     Purpose.
!     --------
!           Initialization of YOMDIMF, and some prints
!**   Interface.
!     ----------
!        *CALL* *SUDIMF2*

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Tomas Wilhelmsson *ECMWF* 
!      Original based on SUDIM2 : 2013-08-02 

! Modifications
! -------------
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMLUN   , ONLY : NULOUT   
USE YOMCT0   , ONLY : NCONF, LNHDYN, LR3D, LR2D
USE YOMDYNA  , ONLY : YRDYNA
USE PARDIM   , ONLY : JPNPPM

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
INTEGER(KIND=JPIM) :: IFD2DD, IFD2DU, IWIND
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUDIMF2',0,ZHOOK_HANDLE)
ASSOCIATE(NUMSPFLDS=>YDML_GCONF%YGFL%NUMSPFLDS, &
 & LADER=>YDML_GCONF%YRDIMF%LADER, LSPT=>YDML_GCONF%YRDIMF%LSPT, LUVDER=>YDML_GCONF%YRDIMF%LUVDER, &
 & LVOR=>YDML_GCONF%YRDIMF%LVOR, NF3D=>YDML_GCONF%YRDIMF%NF3D, NFC2D=>YDML_GCONF%YRDIMF%NFC2D, &
 & NFD2D=>YDML_GCONF%YRDIMF%NFD2D, NFTHER=>YDML_GCONF%YRDIMF%NFTHER, &
 & NS1D=>YDML_GCONF%YRDIMF%NS1D, NS2D=>YDML_GCONF%YRDIMF%NS2D, NS3D=>YDML_GCONF%YRDIMF%NS3D)
!     ------------------------------------------------------------------
!*       1. DEFINE SOME SWITCHES.
!        --------------------------

!     IWIND=0 allows to suppress the wind from the model working arrays
IWIND =1
IF (ANY(SPREAD(NCONF,1,4) == (/923,931,932,933/))) THEN
  IWIND =0
ENDIF

IF(LR3D .OR. ANY(SPREAD(NCONF,1,2) == (/701,901/))) THEN  
  IFD2DD=1
  IF(YRDYNA%LRUBC) THEN
    IFD2DU=1
  ELSE
    IFD2DU=0
  ENDIF
ELSEIF(LR2D) THEN  
  IFD2DD=1
  IFD2DU=0
ELSEIF (ANY(SPREAD(NCONF,1,4) == (/923,931,932,933/))) THEN
  IFD2DD=0
  IFD2DU=0
ENDIF

NFD2D=IFD2DD+IFD2DU

!     ------------------------------------------------------------------
!*       3. Initialize the number of fields in model
!           ----------------------------------------

NFTHER=0
IF(LSPT)  NFTHER=NFTHER+1
IF(LNHDYN) THEN
  NFTHER=NFTHER+2
  ! Add X term among NFTHER
  IF(YRDYNA%LNHX) NFTHER=NFTHER+1
ENDIF

NF3D  =2*IWIND+NFTHER+NUMSPFLDS

IF(.NOT.LADER)THEN
  LUVDER=.FALSE.
ENDIF

! ------------------------------------------------------------------
!*       4. Some printings.
!           ---------------

WRITE(NULOUT,'('' PRINTINGS IN SUDIMF2: NUMBER OF FIELDS '')')
WRITE(UNIT=NULOUT,FMT='('' NFTHER='',I6,'' NF3D  ='',I6,'' NFD2D ='',I6,&
 & '' NFC2D ='',I6,'' NPPM = '',I6)') NFTHER,NF3D,NFD2D,NFC2D,JPNPPM
WRITE(UNIT=NULOUT,FMT='('' LVOR  ='',L2,'' LADER ='',L2,'' LUVDER='',L2)') LVOR,LADER,LUVDER  

!     ------------------------------------------------------------------
!*       6. Initialize Laplace space dimensioning.
!           --------------------------------------

NS3D  =NF3D
NS2D  =NFD2D+NFC2D
NS1D  =2

WRITE(NULOUT,'('' PRINTINGS IN SUDIMF2: LAPLACE SPACE DIMENSIONING '')')
WRITE(NULOUT,'('' NS3D  ='',I6,'' NS2D  ='',I6,'' NS1D  ='',I6)')&
 & NS3D,NS2D,NS1D

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUDIMF2',1,ZHOOK_HANDLE)
END SUBROUTINE SUDIMF2
