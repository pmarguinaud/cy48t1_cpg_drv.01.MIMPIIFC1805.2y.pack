SUBROUTINE INCLI0(YDGEOMETRY,YDSURF,YDGFL,YDEPHY,YDML_PHY_MF,YDMCC)

!**** *INCLI0*

!     PURPOSE.
!     --------

!     This routine prepares the calculation of model climatological constants.

!**   INTERFACE.
!     ----------

!     CALL INCLI0

!     METHOD.
!     -------

!     The parameters are read on *NULNAM* .Then the memory is allocated.
!     Finally, the active routine is called.

!     EXTERNALS.
!     ----------
!      see below

!     AUTHORS.
!     --------
!      M. DEQUE   1 FEB 91.

!     MODIFICATIONS.
!     --------------
!      Y. Bouteloup : 02-03-29 INCLI8 => add A,B and C for climatological ozone profile
!      D. Giard     : 02-11-25 new formulation for extension zone
!      D. Giard     : 02-12-16 LNEWORO3 -> LNEWENV
!      D. Giard     : 02-12-16 spectral smoothing for orography
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      D. Giard     : 03-01-13 spectral smoothing for orography (tunings)
!      D. Giard     : 03-06-16 aqua-planet : (E)INCLI10
!      D. Giard     : 04-09-15 update (parts 2-10)
!      D. Giard     : 05-04-04 semi-envelope and old ARPEGE cost functions
!                              removed, polynomial cost function added, simple 
!                              or no spectral fit for ARPEGE added 
!      M. Hortal    : 05-09-07 spectral smoothing of orography by a diffusion operator 5del16
!      F. Taillefer : 09-06-02 add LZ0THER
!      K. Essaouini : 10-01-25 add LIPGD in namelist
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!      K. Yessad (Jan 2012): remove old 923
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE YOEPHY       , ONLY : TEPHY
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMLUN   , ONLY : NULOUT, NULNAM
USE YOMCT0   , ONLY : LELAM
USE YOMMP0   , ONLY : NPROC
USE YOMMCC   , ONLY : TMCC
USE YOMCST   , ONLY : RPI
USE YOMCLA   , ONLY : NLISSZ, NLISSR, FACZ0, FENVN, FENVS, FACE, QMAX, QMIN,&
 & HDIM, HMIN, QPOWER, QCONST, XINCOC, SCEXT, LNORO, LNLSM, LKEYF,&
 & LNEWORO, LNEWORO2, LIPGD, NLISSP, FLISA, FLISB, LSPSMORO
USE YOMCLI   , ONLY : YRCLI
USE YOM_YGFL , ONLY : TYPE_GFLD

IMPLICIT NONE

TYPE(GEOMETRY),INTENT(INOUT)    :: YDGEOMETRY
TYPE(TSURF)   ,INTENT(INOUT) :: YDSURF
TYPE(TYPE_GFLD)    ,INTENT(INOUT):: YDGFL
TYPE(TEPHY)   ,INTENT(INOUT)    :: YDEPHY
TYPE(TMCC)    ,INTENT(INOUT)    :: YDMCC
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
INTEGER(KIND=JPIM) :: IADL(YDGEOMETRY%YRDIM%NDGLG)
INTEGER(KIND=JPIM) :: IDINT, IPINT, IX, J

REAL(KIND=JPRB) :: ZCV, ZEPS, ZLATN, ZLATS, ZLONE, ZLONW, ZX
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
#include "sualclia.intfb.h"
#include "dealclia.intfb.h"
#include "eincli1.intfb.h"
#include "eincli2.intfb.h"
#include "eincli3.intfb.h"
#include "eincli4.intfb.h"
#include "eincli5.intfb.h"
#include "eincli6.intfb.h"
#include "eincli7.intfb.h"
#include "eincli8.intfb.h"
#include "eincli9.intfb.h"
#include "eincli10.intfb.h"
#include "einter0.intfb.h"
#include "incli1.intfb.h"
#include "incli2.intfb.h"
#include "incli3.intfb.h"
#include "incli4.intfb.h"
#include "incli5.intfb.h"
#include "incli6.intfb.h"
#include "incli7.intfb.h"
#include "incli8.intfb.h"
#include "incli9.intfb.h"
#include "incli10.intfb.h"
#include "inter0.intfb.h"
#include "posnam.intfb.h"

#include "namcla.nam.h"
#include "namcli.nam.h"

LOGICAL :: LIEEE,LGLOBE,LZ0THER
INTEGER (KIND=JPIM) :: NPINT,NDATX,NDATY,NAEROF,NSLICE
REAL (KIND=JPRB) :: ELONSW,ELATSW,ELONNE,ELATNE,SVEG,SFCZ0,RSTR,RSWR

!     ------------------------------------------------------------------
!*
!     1. SET INITIAL VALUES.
!        -------------------

!        1.0 Environment

IF (LHOOK) CALL DR_HOOK('INCLI0',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM, &
  & YDMP=>YDGEOMETRY%YRMP, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDPHY=>YDML_PHY_MF%YRPHY)

ASSOCIATE(NDGLG=>YDDIM%NDGLG, NDGUNG=>YDDIM%NDGUNG, NDGUXG=>YDDIM%NDGUXG, &
 & NDLON=>YDDIM%NDLON, NDLUNG=>YDDIM%NDLUNG, NDLUXG=>YDDIM%NDLUXG, &
 & NMSMAX=>YDDIM%NMSMAX, &
 & NSTAGP=>YDGEM%NSTAGP, NTSTAGP=>YDGEM%NTSTAGP, RSTRET=>YDGEM%RSTRET, &
 & N923=>YDMCC%N923,LVGSN=>YDPHY%LVGSN, LSOLV=>YDPHY%LSOLV)

IF ( NPROC > 1 ) CALL ABOR1(' ONLY ONE PROC. FOR 923 !')
IF (LELAM) THEN
  IADL(1:NDGLG)=NTSTAGP(1:NDGLG)-1
ELSE
  IADL(1:NDGLG)=NSTAGP(1:NDGLG)-1
ENDIF  

ZEPS=1.E-12_JPRB

!        1.1 Set implicit default values

!*** Set values for *COMMON* *YOMCLA*
!  Parameters for orography
!  FENVN  : Factor for the envelope orography on pole of interest.
!           CV=0. OM=SQRT(2.)
!  FENVS  : Factor for the envelope orography on the antipode of the pole of
!           interest. CV=0. OM=0.
!  LNORO  : Key for reading a new orography on a separate file
!  LNLSM  : Key for reading a new land-sea mask on a separate file
!  LIPGD  : Key to take a new GP orographie from a FA PGD file
!  LKEYF  : Key for spectral fit on the orography, only for ALADIN
!  NLISSR : Number of calls to ELISLAP to smooth the orography (ALADIN only)
FENVN = 0._JPRB
FENVS = 0._JPRB
LNORO = .FALSE.
LNLSM = .FALSE.
LIPGD = .FALSE.
LKEYF = .TRUE.
NLISSR= 0
!  Parameters for fitted orography
!  LNEWORO : Switch to the new cost function (Bouteloup)
!  LNEWORO2: Switch to the new cost function (Jerczynski)
!  LSPSMORO: Switch for spectral smoothing of orography by a diffusion operator (linear grid)
!  QMAX    : Maximum value of weight
!  QMIN    : Min--------------------
!  HMIN    : Reference height in the formulation of weight
!  HDIM    : Scaling factor for the orography deviation
!  XINCOC  : Part of the fraction of sea in the weight
!  SCEXT   : Scaling factor for the weight of the extension zone
!  Defaults are ARPEGE ones for the constants
LNEWORO = .FALSE.
LNEWORO2= .FALSE.
LSPSMORO= .FALSE.
QMAX    = 4._JPRB
QMIN    = 2._JPRB
HMIN    = 150._JPRB
HDIM    = 1._JPRB
XINCOC  = 0.0_JPRB
SCEXT   = 0.0_JPRB
!  Additional parameters for the Jerczynski cost function
!  QPOWER  : Exponent in the additional term
!  QCONST  : Scaling factor in the additional term
!  FACE    : 1. for initial formulation , 0. for polynomial one
!  Defaults are ALADIN/LACE ones.
QPOWER = 3.5_JPRB
QCONST = 0.4_JPRB
FACE   = 1.0_JPRB
!  Parameters for the spectral smoothing of orography
!  NLISSP : type of smoothing
!   0 : no, 1 : importing a lower resolution one , 2 : in optimization
!  FLISA  : tuning parameter if NLISSP=2
!  FLISB  : threshold in spectral space if NLISSP=2
NLISSP = 0
FLISB  = 0.0_JPRB
FLISA  = 0.0_JPRB
!  Parameters for roughness length
!  NLISSZ : Number of calls to LISLAP to smooth the roughness length
!  FACZ0  : Scaling factor for the orographic part of Z0 ;CV=1. OM=0.42
NLISSZ = 3
FACZ0  = 1.0_JPRB

!*** Set values for *COMMON* *YOMCLI* 
!  NPINT  : size of the interpolation box ; >= minimum computed by (E)INTER0
!  LIEEE  : if ieee format is used
!  LGLOBE : if global dataset
!  NDATX  : x-size of the dataset (longitude) 
!  NDATY  : y-size of the dataset (latitude) 
!  NAEROF : index of aerosol input file 0=tegen 1=camsaod 2=camsmmr 
!  ELONSW , ELATSW , ELONNE , ELATNE : for a local dataset (INCLI5)
!           latitudes and longitudes of the SW and NE corners (in degrees) 
!  EDLON , EDLAT  : resolution of the dataset (in degrees)       - local
!  NGLOBX, NGLOBY : corresponding numbers of points on the globe - local
!  SVEG   : threshold for significant vegetation cover 
!  LZ0THER: .FALSE. if no orographic part in the thermic roughness length
!  SFCZ0  : scaling factor for the secondary part of z0 (urban., veget.)
!  RSTR   : threshold for temperature-RTT in presence of snow ((E)INCLI6)
!           fixed value of SST-RTT for aqua-planet ((E)INCLI10)
!  RSWR   : threshold for mean moisture   in presence of snow ((E)INCLI6)
!  NSLICE : number of packet used to slice the domain (aladin only)
!  Variables in YOMCLI but not in NAMCLI are set in VAL923
YRCLI%LIEEE  = .TRUE.
YRCLI%LGLOBE = .TRUE.
YRCLI%NPINT = 1
YRCLI%NDATX = 1
YRCLI%NDATY = 1
YRCLI%NAEROF= 0
YRCLI%ELONSW =   0._JPRB
YRCLI%ELATSW = -90._JPRB
YRCLI%ELONNE = 360._JPRB
YRCLI%ELATNE =  90._JPRB
YRCLI%EDLON = 360._JPRB
YRCLI%EDLAT = 180._JPRB
YRCLI%NGLOBX = 1
YRCLI%NGLOBY = 1
YRCLI%SVEG = 0.02_JPRB
YRCLI%LZ0THER = .TRUE.
YRCLI%SFCZ0 = 1._JPRB
YRCLI%RSTR = 5.00_JPRB
YRCLI%RSWR = 0.02_JPRB
YRCLI%NSLICE = 1
IF (.NOT.LSOLV) YRCLI%SVEG = 0._JPRB

!        1.2 Modify default values for YOMCLI according to YOMMCC

IF (N923  ==  1) THEN
  YRCLI%LIEEE =.TRUE.
  YRCLI%NDATX = 8640
  YRCLI%NDATY = 4320
ELSEIF (N923  ==  2) THEN
  IF (LSOLV) THEN
    YRCLI%LIEEE =.TRUE.
    YRCLI%NDATX = 360
    YRCLI%NDATY = 180
  ELSE
    YRCLI%LIEEE =.FALSE.
    YRCLI%NDATX = 432
    YRCLI%NDATY = 216
  ENDIF
ELSEIF (N923  ==  3) THEN
  YRCLI%LIEEE =.TRUE.
  YRCLI%NDATX = 432
  YRCLI%NDATY = 216
ELSEIF (N923  ==  4) THEN
  YRCLI%LIEEE =.TRUE.
  YRCLI%NDATX = 360
  YRCLI%NDATY = 180
ELSEIF (N923  ==  6) THEN
  YRCLI%LIEEE =.TRUE.
  YRCLI%NDATX = 360
  YRCLI%NDATY = 180
ELSEIF (N923  ==  8) THEN    
  YRCLI%LIEEE =.FALSE.
  YRCLI%NDATX = 144
  YRCLI%NDATY = 73                 
ELSEIF (N923  ==  9) THEN    
  YRCLI%LIEEE =.FALSE.
  YRCLI%NDATX = 72
  YRCLI%NDATY = 45                 
ENDIF
YRCLI%NGLOBX = YRCLI%NDATX
YRCLI%NGLOBY = YRCLI%NDATY
YRCLI%EDLON = 360._JPRB/REAL(YRCLI%NDATX,JPRB)
YRCLI%EDLAT = 180._JPRB/REAL(YRCLI%NDATY,JPRB)
IF (LELAM) THEN
  CALL EINTER0(YDGEOMETRY,YRCLI%NGLOBX,YRCLI%NGLOBY,YRCLI%NPINT)
ELSE
  CALL INTER0(NDGLG,NDLON,YRCLI%NGLOBY,YRCLI%NGLOBX,RSTRET,YRCLI%NPINT)
ENDIF
IF (N923  ==  5) THEN
  YRCLI%LGLOBE=.FALSE.
  YRCLI%LIEEE =.TRUE.
  YRCLI%NPINT = 1
  YRCLI%NDATX = 860
  YRCLI%NDATY = 420
  YRCLI%ELONSW = -25.0_JPRB
  YRCLI%ELATSW =  30.0_JPRB
  YRCLI%EDLON = 0.1_JPRB
  YRCLI%EDLAT = 0.1_JPRB
  YRCLI%ELONNE = YRCLI%ELONSW + YRCLI%EDLON*YRCLI%NDATX
  YRCLI%ELATNE = YRCLI%ELATSW + YRCLI%EDLAT*YRCLI%NDATY
  YRCLI%NGLOBX = NINT(360._JPRB/YRCLI%EDLON)
  YRCLI%NGLOBY = NINT(180._JPRB/YRCLI%EDLAT)
ENDIF

!        1.3 Read namelist NAMCLI and modify YOMCLI accordingly

LIEEE     = YRCLI%LIEEE
LGLOBE    = YRCLI%LGLOBE
LZ0THER   = YRCLI%LZ0THER
NPINT     = YRCLI%NPINT
NDATX     = YRCLI%NDATX
NDATY     = YRCLI%NDATY
NAEROF    = YRCLI%NAEROF
NSLICE    = YRCLI%NSLICE
ELONSW    = YRCLI%ELONSW 
ELATSW    = YRCLI%ELATSW 
ELONNE    = YRCLI%ELONNE 
ELATNE    = YRCLI%ELATNE 
SVEG      = YRCLI%SVEG 
SFCZ0     = YRCLI%SFCZ0 
RSTR      = YRCLI%RSTR 
RSWR      = YRCLI%RSWR

CALL POSNAM(NULNAM,'NAMCLI')
READ(NULNAM,NAMCLI)

YRCLI%LIEEE     = LIEEE
YRCLI%LGLOBE    = LGLOBE
YRCLI%LZ0THER   = LZ0THER
YRCLI%NPINT     = NPINT
YRCLI%NDATX     = NDATX
YRCLI%NDATY     = NDATY
YRCLI%NAEROF    = NAEROF
YRCLI%NSLICE    = NSLICE
YRCLI%ELONSW    = ELONSW 
YRCLI%ELATSW    = ELATSW 
YRCLI%ELONNE    = ELONNE 
YRCLI%ELATNE    = ELATNE 
YRCLI%SVEG      = SVEG 
YRCLI%SFCZ0     = SFCZ0 
YRCLI%RSTR      = RSTR 
YRCLI%RSWR      = RSWR

!*** Miscellaneous
YRCLI%RSTR=MAX(0.0_JPRB,YRCLI%RSTR)
YRCLI%RSWR=MAX(0.0_JPRB,MIN(1.0_JPRB,YRCLI%RSWR))

!*** Resolution
YRCLI%NDATX= MAX(YRCLI%NDATX,1)
YRCLI%NDATY= MAX(YRCLI%NDATY,1)
!  Global dataset (described by NDAT-X/Y)
IF (YRCLI%LGLOBE) THEN
  YRCLI%NGLOBX = YRCLI%NDATX
  YRCLI%NGLOBY = YRCLI%NDATY
  YRCLI%EDLON = 360._JPRB/REAL(YRCLI%NDATX,JPRB)
  YRCLI%EDLAT = 180._JPRB/REAL(YRCLI%NDATY,JPRB)
  YRCLI%ELONSW =   0._JPRB
  YRCLI%ELATSW = -90._JPRB
  YRCLI%ELONNE = 360._JPRB
  YRCLI%ELATNE =  90._JPRB
!  Local dataset (described by NDAT-X/Y, ELON-SW/NE, ELAT-SW/NE)
ELSE
  YRCLI%EDLON = (YRCLI%ELONNE-YRCLI%ELONSW)/REAL(YRCLI%NDATX,JPRB)
  YRCLI%EDLAT = (YRCLI%ELATNE-YRCLI%ELATSW)/REAL(YRCLI%NDATY,JPRB)
  IF (YRCLI%EDLON <= ZEPS .OR. YRCLI%EDLAT <= ZEPS)&
   & CALL ABOR1(' PROBLEM WITH LOCAL DATASET CORNERS')  
  IF (YRCLI%EDLON >= 360._JPRB .AND. YRCLI%EDLAT >= 180._JPRB)&
   & CALL ABOR1(' GLOBAL DATASET DECLARED AS LOCAL ?')  
  YRCLI%NGLOBX = NINT(360._JPRB/YRCLI%EDLON)
  YRCLI%NGLOBY = NINT(180._JPRB/YRCLI%EDLAT)
ENDIF

!*** Interpolation box or operator
IF (N923 == 5 .OR. N923 == 7) THEN
! Control in parts 5 or 7
  YRCLI%NPINT = MIN(2,MAX(1,YRCLI%NPINT))
ELSE
! Size for parts 1, 2, 3, 4, 6, 8, 9, 10 and old options B and R
!  Minimum size, as a function of model geometry and data resolution
  IPINT = 1
  IF (LELAM) THEN
    CALL EINTER0(YDGEOMETRY,YRCLI%NGLOBX,YRCLI%NGLOBY,IPINT)
  ELSE
    CALL INTER0(NDGLG,NDLON,YRCLI%NGLOBY,YRCLI%NGLOBX,RSTRET,IPINT)
  ENDIF
!  Control (NPINT must be odd and at least IPINT)
  YRCLI%NPINT = (YRCLI%NPINT/2)*2+1
  YRCLI%NPINT = MAX(IPINT,YRCLI%NPINT)
!  Temporary warning for interpolators
  IF (MOD(YRCLI%NDATY,2) == 1) WRITE(UNIT=NULOUT,&
   & FMT='(''CAUTION : INTERPOLATORS MAY NOT SUIT ODD NDATY'')')
ENDIF

!*** Control of boundaries
!  ALADIN
IF (LELAM .AND. .NOT.YRCLI%LGLOBE) THEN
  IDINT= YRCLI%NPINT/2
  ZCV = 180._JPRB/RPI
!  Latitude boundary
  ZLATN = 0.0_JPRB
  ZLATS = RPI
  DO J=NDLUNG,NDLUXG
    ZLATN = MAX(ZLATN, YDGSGEOM_NB%GELAT(J+IADL(NDGUXG)) )
    ZLATS = MIN(ZLATS, YDGSGEOM_NB%GELAT(J+IADL(NDGUNG)) )
  ENDDO
  ZLATN = ZLATN*ZCV + (IDINT+1)*YRCLI%EDLAT
  ZLATS = ZLATS*ZCV - (IDINT+1)*YRCLI%EDLAT
  IF (N923 == 5 .OR. N923 == 7) THEN
    IF (YRCLI%ELATSW  >=  ZLATN) THEN
      WRITE(UNIT=NULOUT,FMT='(''BAD POSITION OF DATASET :'',&
       & '' ELATSW ='',F8.2,'' >'',F8.2,'' (dg)'')') YRCLI%ELATSW,ZLATN  
      IF (LHOOK) CALL DR_HOOK('INCLI0',1,ZHOOK_HANDLE)
      RETURN
    ENDIF
    IF (YRCLI%ELATNE  <=  ZLATS) THEN
      WRITE(UNIT=NULOUT,FMT='(''BAD POSITION OF DATASET :'',&
       & '' ELATNE ='',F8.2,'' <'',F8.2,'' (dg)'')') YRCLI%ELATNE,ZLATS  
      IF (LHOOK) CALL DR_HOOK('INCLI0',1,ZHOOK_HANDLE)
      RETURN
    ENDIF
  ELSE
    IF (YRCLI%ELATSW  >  ZLATS) THEN
      WRITE(UNIT=NULOUT,FMT='(''BAD POSITION OF DATASET :'',&
       & '' ELATSW ='',F8.2,'' >'',F8.2,'' (dg)'')') YRCLI%ELATSW,ZLATS  
      CALL ABOR1('INCLI0')
    ENDIF
    IF (YRCLI%ELATNE  <  ZLATN) THEN
      WRITE(UNIT=NULOUT,FMT='(''BAD POSITION OF DATASET :'',&
       & '' ELATNE ='',F8.2,'' <'',F8.2,'' (dg)'')') YRCLI%ELATNE,ZLATN  
      CALL ABOR1('INCLI0')
    ENDIF
  ENDIF
!  Longitude boundary
  ZX = 1.0_JPRB/MAX(ZEPS, SIN(RPI-MAX(ABS(ZLATN),ABS(ZLATS))/ZCV))
  IX = MIN(YRCLI%NGLOBX/4,INT(ZX))
  ZLONW = YDGSGEOM_NB%GELAM(     1+IADL(NDGUXG))*ZCV
  ZLONE = YDGSGEOM_NB%GELAM(NDLUXG+IADL(NDGUXG))*ZCV
  ZLONW = ZLONW-(1.0_JPRB-SIGN(1.0_JPRB,180._JPRB-ZLONW))*180._JPRB
  ZLONE = ZLONE-(1.0_JPRB-SIGN(1.0_JPRB,180._JPRB-ZLONE))*180._JPRB
  ZLONW = ZLONW - (IX+IDINT)*YRCLI%EDLON
  ZLONE = ZLONE + (IX+IDINT)*YRCLI%EDLON
  IF (N923 == 5 .OR. N923 == 7) THEN
    IF (YRCLI%ELONSW  >=  ZLONE) THEN
      WRITE(UNIT=NULOUT,FMT='(''BAD POSITION OF DATASET :'',&
       & '' ELONSW ='',F8.2,'' >'',F8.2,'' (dg)'')') YRCLI%ELONSW,ZLONE  
      IF (LHOOK) CALL DR_HOOK('INCLI0',1,ZHOOK_HANDLE)
      RETURN
    ENDIF
    IF (YRCLI%ELONNE  <=  ZLONW) THEN
      WRITE(UNIT=NULOUT,FMT='(''BAD POSITION OF DATASET :'',&
       & '' ELONNE ='',F8.2,'' <'',F8.2,'' (dg)'')') YRCLI%ELONNE,ZLONW  
      IF (LHOOK) CALL DR_HOOK('INCLI0',1,ZHOOK_HANDLE)
      RETURN
    ENDIF
  ELSE
    IF (YRCLI%ELONSW  >  ZLONW) THEN
      WRITE(UNIT=NULOUT,FMT='(''BAD POSITION OF DATASET :'',&
       & '' ELONSW ='',F8.2,'' >'',F8.2,'' (dg)'')') YRCLI%ELONSW,ZLONW  
      CALL ABOR1('INCLI0')
    ENDIF
    IF (YRCLI%ELONNE  <  ZLONE) THEN
      WRITE(UNIT=NULOUT,FMT='(''BAD POSITION OF DATASET :'',&
       & '' ELONNE ='',F8.2,'' <'',F8.2,'' (dg)'')') YRCLI%ELONNE,ZLONE  
      CALL ABOR1('INCLI0')
    ENDIF
  ENDIF
ENDIF
!  ARPEGE
IF (.NOT.LELAM .AND. .NOT.YRCLI%LGLOBE .AND. N923 /= 5 .AND. N923 /= 7) THEN
  CALL ABOR1(' GLOBAL DATASETS REQUIRED')
ENDIF

!*** Temporary control for the albedos of vegetation and bare-ground
IF (LELAM .AND. (N923 == 2 .OR. N923 == 4 .OR. N923 ==5) .AND. .NOT.LVGSN)&
 & CALL ABOR1(' SET LVGSN=.T. IN NAMPHY FOR SAFETY')  

!        1.4 Print dataset description

WRITE(UNIT=NULOUT,FMT='('' COMMON YOMCLI''/&
 & '' LIEEE='',L2,'' LGLOBE='',L2,'' NDATX='',I4,'' NDATY='',I4,'' NAEROF='',I4,&
 & '' NPINT='',I3,'' EDLON (dg)='',F5.2,'' EDLAT (dg)='',F5.2,/,&
 & '' SVEG = '',F5.2,'' LZ0THER = '',L2,'' SFCZ0 = '',F5.2,&
 & '' RSTR = '',F3.0,''K  RSWR = '',F4.2,'' NSLICE = '',I4)')&
 & YRCLI%LIEEE,YRCLI%LGLOBE,YRCLI%NDATX,YRCLI%NDATY,YRCLI%NAEROF,YRCLI%NPINT,YRCLI%EDLON,YRCLI%EDLAT,&
 & YRCLI%SVEG,YRCLI%LZ0THER,YRCLI%SFCZ0,YRCLI%RSTR,YRCLI%RSWR,YRCLI%NSLICE  
IF (.NOT.YRCLI%LGLOBE) THEN
  WRITE(UNIT=NULOUT,FMT='(&
   & '' ELATSW (dg)='',F8.2,'' ELONSW (dg)='',F8.2,&
   & '' ELATNE (dg)='',F8.2,'' ELONNE (dg)='',F8.2,&
   & '' NGLOBX='',I5,'' NGLOBY='',I5)')&
   & YRCLI%ELATSW,YRCLI%ELONSW,YRCLI%ELATNE,YRCLI%ELONNE,YRCLI%NGLOBX,YRCLI%NGLOBY  
ENDIF

!        1.5  Read and check topography description (if required)

IF ( N923 == 1 ) THEN

  CALL POSNAM(NULNAM,'NAMCLA')
  READ(NULNAM,NAMCLA)

  SCEXT  = MAX(0.0_JPRB,SCEXT)
  NLISSP = MAX(0,MIN(2,NLISSP))

  ZEPS=1.E-4_JPRB

  IF (NLISSP == 1 .AND. .NOT. LNORO) THEN
    CALL ABOR1(' IMPORTING SPECTRAL OROGRAPHY REQUIRES LNORO=.TRUE. !')
  ENDIF
  IF ((LNEWORO.OR.LNEWORO2.OR.LNORO) .AND. .NOT.LKEYF .AND. .NOT.LIPGD) THEN
    CALL ABOR1(' IMPORTING OR OPTIMIZING ONLY SPECTRAL OROGRAPHY !')
  ENDIF
  IF (NLISSP == 2 .AND. .NOT.(LNEWORO.OR.LNEWORO2)) THEN
    CALL ABOR1(' SPECTRAL COST FUNCTION ALONE FAILS !')
  ENDIF
  IF (LIPGD .AND. .NOT.LELAM) THEN
    CALL ABOR1(' IMPORTING FROM A PGD FILE ALLOWED ONLY FOR ALADIN ')
  ENDIF

  IF (LNORO) THEN
    IF (.NOT.LELAM) THEN
      LNEWORO  = .FALSE.
      LNEWORO2 = .FALSE.
    ENDIF
  ELSE
    LNLSM = .FALSE.
  ENDIF

  WRITE(UNIT=NULOUT,FMT='('' COMMON YOMCLA''/&
  & '' LKEYF '',L2,'' LNEWORO '',L2,'' LNEWORO2 '',L2,'' LSPSMORO '',L2,&
  & '' XINCOC '',F8.3,'' SCEXT '',F8.3,/,&
  & '' QMAX '',F6.1,'' QMIN '',F6.1,'' HDIM '',F6.2,'' HMIN '',F6.0,&
  & '' QPOWER '',F6.3,'' QCONST '',F6.3,'' FACE '',F6.3,/,&
  & '' FENVN '',F6.3,'' FENVS '',F6.3,/,&
  & '' LNORO '',L2,'' LNLSM '',L2,'' LIPGD '',L2,/,&
  & '' NLISSZ '',I2,'' FACZ0'',F6.3,/,'' NLISSR '',I2,&
  & '' NLISSP '',I2,'' FLISA '',F6.3,'' FLISB '',F6.1,//)')&
  & LKEYF,LNEWORO,LNEWORO2,LSPSMORO,XINCOC,SCEXT,QMAX,QMIN,HDIM,HMIN,QPOWER,QCONST,&
  & FACE,FENVN,FENVS,LNORO,LNLSM,LIPGD,NLISSZ,FACZ0,NLISSR,NLISSP,FLISA,FLISB

  IF ((MOD(NDLON,NMSMAX) > 2) .AND. LSPSMORO) THEN
    CALL ABOR1(' THE USE OF LSPSMORO REQUIRES A LINEAR GRID !')
  ENDIF
  IF (LSPSMORO .AND. (LNEWORO.OR.LNEWORO2)) THEN
    CALL ABOR1('ONLY ONE TYPE OF SPECTRAL SMOOTHING AT THE TIME !')
  ENDIF
  IF (LNEWORO .AND. LNEWORO2) THEN
    CALL ABOR1('ONLY ONE TYPE OF SPECTRAL SMOOTHING AT THE TIME !')
  ENDIF

ENDIF

!     ------------------------------------------------------------------
!*
!     2  ACTIVE CALCULATION.
!        -------------------

!*       2.1   MEMORY ALLOCATION

IF ( N923 == 1 ) THEN
  CALL SUALCLIA(YDGEOMETRY)
ENDIF

!*       2.2   NEW C923

!  TOPOGRAPHY
IF (N923 == 1) THEN
  IF (LELAM) THEN
    CALL EINCLI1(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF)
  ELSE
    CALL INCLI1(YDGEOMETRY)
  ENDIF
ENDIF
!  FIXED FIELDS
IF (N923 == 2) THEN
  IF (LELAM) THEN
    CALL EINCLI2(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF)
  ELSE
    CALL INCLI2(YDGEOMETRY,YDPHY)
  ENDIF
ENDIF
!  OLD CLIMATOLOGY FOR SOIL TEMPERATURE AND MOISTURE
IF (N923 == 3) THEN
  IF (LELAM) THEN
    CALL EINCLI3(YDGEOMETRY,YDGFL,YDSURF,YDEPHY,YDML_PHY_MF)
  ELSE
    CALL INCLI3(YDGEOMETRY,YDSURF,YDPHY,YDML_PHY_MF%YRPHY1)
  ENDIF
ENDIF
!  OTHER MONTHLY VARYING FIELDS - VEGETATION
IF (N923 == 4) THEN
  IF (LELAM) THEN
    CALL EINCLI4(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF)
  ELSE
    CALL INCLI4(YDGEOMETRY)
  ENDIF
ENDIF
!  LOCAL MODIFICATION OF FIELDS ON LAND
IF (N923 == 5) THEN
  IF (LELAM) THEN
    CALL EINCLI5(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF)
  ELSE
    CALL INCLI5(YDGEOMETRY)
  ENDIF
ENDIF
!  NEW CLIMATOLOGY FOR SOIL TEMPERATURE AND MOISTURE
IF (N923 == 6) THEN
  IF (LELAM) THEN
    CALL EINCLI6(YDGEOMETRY,YDGFL,YDSURF,YDEPHY,YDML_PHY_MF)
  ELSE
    CALL INCLI6(YDGEOMETRY,YDSURF,YDPHY,YDML_PHY_MF%YRPHY1)
  ENDIF
ENDIF
!  LOCAL MODIFICATION OF FIELDS ON SEA/LAKES
IF (N923 == 7) THEN
  IF (LELAM) THEN
    CALL EINCLI7(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF)
  ELSE
    CALL INCLI7(YDGEOMETRY,YDML_PHY_MF%YRPHY1)
  ENDIF
ENDIF
!  A, B, C COEFFICIENTS FOR OZONE PROFILES
IF (N923 == 8) THEN
  IF (LELAM) THEN
    CALL EINCLI8(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF)
  ELSE
    CALL INCLI8(YDGEOMETRY)
  ENDIF
ENDIF
!  AEROSOLS (NEW)
IF (N923 == 9) THEN
  IF (LELAM) THEN
    CALL EINCLI9(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF)
  ELSE
    CALL INCLI9(YDGEOMETRY)
  ENDIF
ENDIF
!  AQUA-PLANET
IF (N923 == 10) THEN
  IF (LELAM) THEN
    CALL EINCLI10(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF)
  ELSE
    CALL INCLI10(YDGEOMETRY,YDML_PHY_MF%YRPHY1)
  ENDIF
ENDIF
!  SPECTRAL EMISSIVITY
IF (N923 == 11) THEN
  CALL ABOR1('INCLI0 : PART 11 NOT YET READY')
ENDIF

!*       2.4   MEMORY DEALLOCATION

IF ( N923 == 1 ) THEN
  CALL DEALCLIA
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INCLI0',1,ZHOOK_HANDLE)
END SUBROUTINE INCLI0
