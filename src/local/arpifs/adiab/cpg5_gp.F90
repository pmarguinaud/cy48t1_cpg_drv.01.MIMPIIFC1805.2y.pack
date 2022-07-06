SUBROUTINE CPG5_GP(YDGEOMETRY,YDGMV5,YDSURF,&
 !---------------------------------------------------------------------
 ! - INPUT .
 & YDML_GCONF,YDPHLC,YDSIMPHL,YDML_DYN,KST,KEND,KSTGLO,&
 & KIBL,&
 & PTRAJ_PHYS,PTRAJ_SRFC,PTRAJ_CST,&
 !---------------------------------------------------------------------
 ! - INPUT/OUTPUT .
 & PGFL5,PGMV5,PGMV5S,&
 !---------------------------------------------------------------------
 ! - OUTPUT .
 & POROGL,POROGM,&
 & PGMVT95,PGFLT95,&
 & PRE5,PRE5L,PRE5M,PRE5F,PXYB5,PRPREF5,&
 & PUVH5,PHI5,PHIF5,PRCP5,PRRED5,PCTY5,PRT5,PRT5L,PRT5M,&
 & PRE95,PRE95F,PXYB95,PRRED95,PHI95,PHIF95,PRCP95,&
 & PGPPC5,PQS,PQS1,PTS,PTS1,PSNS,PSD_VF5,PQCHA5L,PQCHA5M,&
 & PXYBDER5,PRNHPPI5,PHI5FL,PHI5FM,PNH15L,PNH15M)

!**** *CPG5_GP* - Grid point calculations:
!                initial part of not lagged grid-point calculations
!                for the trajectory in TL and AD codes.

!     Purpose.
!     --------
!           Grid point calculations:
!           initial part of not lagged grid-point calculations,
!           for the trajectory in TL and AD codes.
!           - get data in buffers.
!           - multiply p-order horizontal derivatives by M**p.
!           - calls some GP... routines to initialise some auxiliary variables.

!           This routine has the same structure as CPG_GP, but it works on
!           trajectory arrays (the name of which ending by 5 or 95)
!           and does some calls to RDPHTR.. routines. Only a subset of CPG_GP
!           is done under CPG5_GP. The call of CPG5_GP is due to the fact that,
!           when computing the trajectory in the direct code, all the variables
!           computed by CPG_GP are not saved on buffers (too expensive in memory):
!           some of them need to be recomputed by CPG5_GP when calling the TL
!           or AD code.

!           Abbreviation "vwv" stands for "vertical wind variable".

!**   Interface.
!     ----------
!        *CALL* *CPG5_GP(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        KST       : first element of work.
!        KEND      : last element of work.
!        KSTGLO    : global offset.
!        KIBL      : index into YDGSGEOM instance in YDGEOMETRY
!        PTRAJ_PHYS : structure for NL trajectory used in physics relevant to t
!        PTRAJ_SRFC : structure for NL trajectory with surface fields
!        PTRAJ_CST  : structure for NL trajectory with climate fields

!     INPUT/OUTPUT:
!     -------------
!        PGFL5     : unified_treatment grid-point fields at t
!        PGMV5     : upper air GMV variables at time t and t-dt.
!        PGMV5S    : surface GMV variables at time t and t-dt.

!     OUTPUT:
!     -------
!        POROGL    : zonal component of "grad(surf orography)"
!        POROGM    : meridian component of "grad(surf orography)"
!        PGMVT95   : GMV at time t-dt at full levels.
!        PGFLT95   : GFL at time t-dt at full levels.
!        PRE5      : hydrostatic pressure "prehyd" at half levels at time t.
!        PRE5L     : zonal component of "grad prehyds" at t.
!        PRE5M     : meridian component of "grad prehyds" at t.
!        PRE5F     : hydrostatic pressure "prehyd" at full levels at time t.
!        PXYB5     : contains pressure depth, "delta", "alpha" at t.
!        PRPREF5   : 1/PRE5F.
!        PUVH5     : horizontal wind at time t at half levels.
!        PHI5      : geopotential height "gz" at t at half levels.
!        PHIF5     : geopotential height "gz" at t at full levels.
!        PRCP5     : contains "cp", "R" and "Kap=R/Cp" at t.
!        PRRED5    : reduced gas constant by (1.+Prond) in NH at t.
!        PCTY5     : contains vertical velocities, vertical integral of divergence at t.
!        PRT5      : RT at full levels at t.
!        PRT5L     : zonal component of "grad RT" at full levels at t.
!        PRT5M     : meridian component of "grad RT" at full levels at t.
!        PRE95     : hydrostatic pressure "prehyd" at half levels at t-dt.
!        PRE95F    : hydrostatic pressure "prehyd" at full levels at time t-dt.
!        PXYB95    : contains pressure depth, "delta", "alpha" at t-dt.
!        PRRED95   : reduced gas constant by (1.+Prond) in NH at t-dt.
!        PHI95     : geopotential height "gz" at t-dt at half levels.
!        PHIF95    : geopotential height "gz" at t-dt at full levels.
!        PRCP95    : contains "cp", "R" and "Kap=R/Cp" at t-dt.
!        PGPPC5    : buffer containing some fields used in the PC schemes.
!        PQS       : surface specific humidity at t or t-dt accord. to ltwotl.
!        PQS1      : surface specific humidity at t+dt.
!        PTS       : surface temperature at t or t-dt accord. to ltwotl.
!        PTS1      : surface temperature at t+dt.
!        PSNS      : mass of snow per unit surf at t or t-dt accord. to ltwotl.
!        PSD_VF5   : surface fields(recovered by GET_TRAJ_SFC).
!        PQCHA5L   : zonal comp grad(log(pre/prehyd)).
!        PQCHA5M   : merid comp grad(log(pre/prehyd)).
!        PXYBDER5  : cf. PXYBDER in GPGRXYB (time t).
!        PRNHPPI5  : "prehyd/pre" at full levels (time t).
!        PHI5FL    : zonal comp grad(gz) at full levels
!        PHI5FM    : merid comp grad(gz) at full levels
!        PNH15L    : zonal comp of RT grad(log(prehyd/pre)) at full layers
!        PNH15M    : merid comp of RT grad(log(prehyd/pre)) at full layers

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE documentation vol 2 ch 1 and vol 3 ch 6

!     Author.
!     -------
!        K. YESSAD, after part 3 of CPGTL and CPGAD.
!        Original : 08 Jul 2004

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K.Yessad (Aug 2009): remove LPC_OLD in TL and AD codes.
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures + cleanings.
!   O.Riviere (Feb 11): add QL/QI in TL/AD
!   K. Yessad (Nov 2011): various contributions.
!   K. Yessad (Nov 2012): simplify testings.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   F. Vana  28-Nov-2013 : Redesigned trajectory handling.
!   K. Yessad (July 2014): Move some variables.
!   K. Yessad (June 2017): Introduce NHQE model.
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   H Petithomme (Dec 2020): optimisation on gphpre 3TL
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE SURFACE_FIELDS_OPER_TRAJ , ONLY : GPOPER
USE YOMGMV   , ONLY : TGMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCT0   , ONLY : LNHDYN, LTWOTL, LSPRT
USE YOMDYNA  , ONLY : LRUBC, LPC_FULL
!!! USE YOMCT3   , ONLY : NSTEP
USE YOMCST   , ONLY : RD, RV
USE YOMCVER  , ONLY : LVERTFE
USE YOMSIMPHL, ONLY : TSIMPHL
USE YOPHLC   , ONLY : TPHLC
USE INTDYN_MOD,ONLY : YYTHW5, YYTCTY5, YYTXYBDER5,&
 & YYTRCP5, YYTRCP95, YYTXYB5, YYTXYB95, YYTGMVT95, YYTGFLT95
USE YOMTRAJ  , ONLY : TRAJ_PHYS_TYPE, TRAJ_CST_TYPE, TRAJ_SRFC_TYPE,&
 & LPRTTRAJ
USE YOMLUN   , ONLY : NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)          ,INTENT(INOUT) :: YDGMV5
TYPE(TSURF)         ,INTENT(IN)    :: YDSURF
TYPE(MODEL_DYNAMICS_TYPE),INTENT(INOUT):: YDML_DYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TPHLC)         ,INTENT(INOUT) :: YDPHLC
TYPE(TSIMPHL)       ,INTENT(INOUT) :: YDSIMPHL
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KST 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KSTGLO 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KIBL
TYPE(TRAJ_PHYS_TYPE),INTENT(IN)    :: PTRAJ_PHYS
TYPE(TRAJ_SRFC_TYPE),INTENT(IN)    :: PTRAJ_SRFC
TYPE(TRAJ_CST_TYPE) ,INTENT(IN)    :: PTRAJ_CST
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PGFL5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM5) 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PGMV5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV5%YT5%NDIM) 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PGMV5S(YDGEOMETRY%YRDIM%NPROMA,YDGMV5%YT5%NDIMS) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: POROGL(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: POROGM(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PGMVT95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTGMVT95%NDIM)
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PGFLT95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTGFLT95%NDIM)
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PRE5(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PRE5L(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PRE5M(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PRE5F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PXYB5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB5%NDIM) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PRPREF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PUVH5(YDGEOMETRY%YRDIM%NPROMNH,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTHW5%NDIM)
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PHI5(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PHIF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PRCP5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP5%NDIM)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PRRED5(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PCTY5(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTCTY5%NDIM)
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PRT5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PRT5L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PRT5M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PRE95(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PRE95F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PXYB95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB95%NDIM) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PRRED95(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PHI95(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PHIF95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PRCP95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP95%NDIM)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PGPPC5(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRGPPC%NFGPPC+1) ! curr. not used
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PQS(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PQS1(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PTS(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PTS1(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PSNS(YDGEOMETRY%YRDIM%NPROMM) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PSD_VF5(YDGEOMETRY%YRDIM%NPROMA,YDSURF%YSD_VFD%NDIM)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PQCHA5L(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PQCHA5M(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PXYBDER5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYBDER5%NDIM)
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PRNHPPI5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PHI5FL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PHI5FM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PNH15L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PNH15M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: IFLAG
!!! INTEGER(KIND=JPIM) :: INHFIELDS
INTEGER(KIND=JPIM) :: JLEV, JROF, KNGP5, KNTRAJ_CST

LOGICAL :: LLSTR

REAL(KIND=JPRB) :: ZGPH5L(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZGPH5M(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZEVT0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZPSPT95(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZNHPPI5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)    ! "pre/prehyd" at full levels (time t).
REAL(KIND=JPRB) :: ZOROGLL(YDGEOMETRY%YRDIM%NPROMNH)
REAL(KIND=JPRB) :: ZOROGMM(YDGEOMETRY%YRDIM%NPROMNH)
REAL(KIND=JPRB) :: ZOROGLM(YDGEOMETRY%YRDIM%NPROMNH)

REAL(KIND=JPRB), ALLOCATABLE :: ZTEMP(:,:),ZTEMPC(:,:)

REAL(KIND=JPRB) :: ZDUM3(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,3)
REAL(KIND=JPRB) :: ZDUM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDUMSPT(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZEPS,ZEPS2
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "gpcty.intfb.h"
#include "gpgeo.intfb.h"
#include "gpgrgeo.intfb.h"
#include "gpgrxyb.intfb.h"
#include "gphluv.intfb.h"
#include "gphlwi.intfb.h"
#include "gpmpfc5.intfb.h"
#include "gphpre.intfb.h"
#include "gprcp_qlirsg.intfb.h"
#include "gprt.intfb.h"
#include "gp_spv.intfb.h"
#include "rdphtrajm.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPG5_GP',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDOROG=>YDGEOMETRY%YROROG(KIBL), YDPTRGPPC=>YDML_DYN%YRPTRGPPC, &
 & YDDYN=>YDML_DYN%YRDYN,YGFL=>YDML_GCONF%YGFL)

ASSOCIATE(NDIM5=>YGFL%NDIM5, YI=>YGFL%YI, YL=>YGFL%YL, YQ=>YGFL%YQ, &
 & NPROMA=>YDDIM%NPROMA, NPROMM=>YDDIM%NPROMM, NPROMNH=>YDDIM%NPROMNH, &
 & RSTRET=>YDGEM%RSTRET, &
 & YSD_VFD=>YDSURF%YSD_VFD, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & YT5=>YDGMV5%YT5, &
 & LTRAJPS=>YDSIMPHL%LTRAJPS, LPROCLDTL=>YDSIMPHL%LPROCLDTL, &
 & LSIMPH=>YDSIMPHL%LSIMPH, &
 & LSPHLC=>YDPHLC%LSPHLC, &
 & NFGPPC=>YDPTRGPPC%NFGPPC, &
 & SITLAF=>YDDYN%SITLAF)
!     ------------------------------------------------------------------

!*       0.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

ZEPS=100.0_JPRB*TINY(1.0_JPRB)
ZEPS2=1.E-12_JPRB
LLSTR=(ABS(RSTRET-1.0_JPRB)>ZEPS)
ZDUM3(:,:,:)=0.0_JPRB
ZDUM(:,:)=0.0_JPRB
ZDUMSPT(:)=0.0_JPRB

! Protection against negative condensed water species after advection
IF (LPROCLDTL) THEN
  IF (LTWOTL) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        IF (PGFL5(JROF,JLEV,YL%MP5) < ZEPS2) PGFL5(JROF,JLEV,YL%MP5)=0.0_JPRB
        IF (PGFL5(JROF,JLEV,YI%MP5) < ZEPS2) PGFL5(JROF,JLEV,YI%MP5)=0.0_JPRB
      ENDDO
    ENDDO
  ELSE
    ! ky: do we have to modify PGFL5(:,:,YL%MP5) and PGFL5(:,:,YI%MP5) too?

    ! ky: the following piece of code found in CY38 is not correct, work must
    !     be done on PGFLT95 after calling RDPHTRAJM.
    !!bugged DO JLEV=1,NFLEVG
    !!bugged   DO JROF=KST,KEND
    !!bugged     IF (PGFL5(JROF,JLEV,YL%MP9) < ZEPS2) PGFL5(JROF,JLEV,YL%MP9)=0.0_JPRB
    !!bugged     IF (PGFL5(JROF,JLEV,YI%MP9) < ZEPS2) PGFL5(JROF,JLEV,YI%MP9)=0.0_JPRB
    !!bugged   ENDDO
    !!bugged ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       1.    INTERFACE TO GLOBAL ARRAYS/WORKFILES, AND READ TRAJECTORY.
!              ----------------------------------------------------------

!*       1.1    Geographical derivatives of orography
!               -------------------------------------

!DIR$ IVDEP
!CDIR NODEP
DO JROF=KST,KEND
  POROGL(JROF)=YDOROG%OROGL(JROF)*YDGSGEOM%GM(JROF)
  POROGM(JROF)=YDOROG%OROGM(JROF)*YDGSGEOM%GM(JROF)
ENDDO

IF(LNHDYN) THEN
  ! ky: curr. not used; may be necessary later in NH
  DO JROF=KST,KEND
    ZOROGLL(JROF)=YDOROG%OROGLL(JROF)*YDGSGEOM%GM(JROF)**2
    ZOROGMM(JROF)=YDOROG%OROGMM(JROF)*YDGSGEOM%GM(JROF)**2
    ZOROGLM(JROF)=YDOROG%OROGLM(JROF)*YDGSGEOM%GM(JROF)**2
  ENDDO
ENDIF

!*       1.2    Initialization of NHS traj array (iteration)
!               --------------------------------------------

!*       1.2.1  Read stored model trajectory at t-dt.

! simplified physics
IF (LTRAJPS) THEN
  CALL RDPHTRAJM(YDGEOMETRY,YDSIMPHL,KST,KEND,PTRAJ_PHYS,&
   & PGMVT95(1,1,YYTGMVT95%M_U),PGMVT95(1,1,YYTGMVT95%M_V),PGMVT95(1,1,YYTGMVT95%M_T),&
   & PGFLT95(1,1,YYTGFLT95%M_Q),PGFLT95(1,1,YYTGFLT95%M_L),PGFLT95(1,1,YYTGFLT95%M_I),ZPSPT95)
  IF (LPROCLDTL .AND. .NOT.LTWOTL) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        IF (PGFLT95(JROF,JLEV,YYTGFLT95%M_L) < ZEPS2) PGFLT95(JROF,JLEV,YYTGFLT95%M_L)=0.0_JPRB
        IF (PGFLT95(JROF,JLEV,YYTGFLT95%M_I) < ZEPS2) PGFLT95(JROF,JLEV,YYTGFLT95%M_I)=0.0_JPRB
      ENDDO
    ENDDO
  ENDIF
  PQS(KST:KEND) =PTRAJ_PHYS%PQSSMF(KST:KEND)
  PTS(KST:KEND) =PTRAJ_PHYS%PTSMF(KST:KEND)
  PSNS(KST:KEND)=PTRAJ_PHYS%PSNSMF(KST:KEND)
  PQS1(KST:KEND)=PTRAJ_PHYS%PQSS1MF(KST:KEND)
  PTS1(KST:KEND)=PTRAJ_PHYS%PTS1MF(KST:KEND)
  IF (LPRTTRAJ.AND.PTRAJ_PHYS%LASTCHUNK) WRITE(NULOUT,*)'GREPTRAJ READ TRAJ_PHYS in CPG5_GP'
ENDIF

! NHS predictive step
!!! Obsolete way to manage NH trajectory, to recode properly if developpement of some TL/AD NH code.
!!! IF (LNHDYN.AND.LTRAJNH) THEN
!!!   INHFIELDS=NG3NH95*NFLEVG+1
!!!   CALL RDNHTRAJM(YDGEOMETRY,YDML_DYN%YRTNH,NSTEP,KST,KEND,KSTGLO,INHFIELDS,&
!!!    & PGMVT95(1,1,YYTGMVT95%M_SPD),PGMVT95(1,1,YYTGMVT95%M_SVD),PGMVT95(1,1,YYTGMVT95%M_DIV),&
!!!    & PGMVT95(1,1,YYTGMVT95%M_U),PGMVT95(1,1,YYTGMVT95%M_V),PGMVT95(1,1,YYTGMVT95%M_T),ZPSPT95)
!!! ENDIF

!*       1.3    Get surface variables from incore traj
!               --------------------------------------

IF (LSIMPH.OR.LSPHLC) THEN
  ! Working arrays required for different KIND conversion
  KNGP5=SIZE(PTRAJ_SRFC%MIX_SRFC,2)
  KNTRAJ_CST=SIZE(PTRAJ_CST%MIX_CST,2)
  ALLOCATE (ZTEMP(NPROMA,KNGP5))
  ALLOCATE (ZTEMPC(NPROMA,KNTRAJ_CST))
  ZTEMP(:,:)=0.0_JPRB
  ZTEMPC(:,:)=0.0_JPRB
  ZTEMP(KST:KEND,:)=PTRAJ_SRFC%MIX_SRFC(KST:KEND,:)
  ZTEMPC(KST:KEND,:)=PTRAJ_CST%MIX_CST(KST:KEND,:)
  CALL  GPOPER(YDGEOMETRY%YRDIM,YDDYN,'GETTRAJ',YDSURF,PSD_VF=PSD_VF5,PFIELD=ZTEMP,PFIELD2=ZTEMPC)
  IF (LPRTTRAJ.AND.PTRAJ_SRFC%LASTCHUNK) WRITE(NULOUT,*)'GREPTRAJ READ TRAJ_CST and TRAJ_SRFC in CPG5_GP'
  DEALLOCATE(ZTEMP,ZTEMPC)
ENDIF

!*       1.4    Map factor
!               ----------

IF(LLSTR) THEN
  IFLAG=0
  CALL GPMPFC5(YDGMV5,YDML_GCONF,NPROMA,NFLEVG,KST,KEND,IFLAG,YDGSGEOM%GM,PGMV5,PGMV5S,PGFL5)
ENDIF

!*       1.5    surface pressure variable and surface variables.
!               ------------------------------------------------

CALL GP_SPV(LNHDYN,LPC_FULL,LTWOTL,YDGEOMETRY,YDDYN,YDSIMPHL,.TRUE.,KST,KEND,&
 & PGMV5S(1,YT5%MSP),PGMV5S(1,YT5%MSPL),PGMV5S(1,YT5%MSPM),&
 & ZPSPT95,ZDUMSPT,ZDUMSPT,&
 & PRE5,PRE5L,PRE5M,PRE95,ZDUMSPT,ZDUMSPT)

!     ------------------------------------------------------------------

!*       3.    INITIALISE AUXILIARY VARIABLES.
!              -------------------------------

!*       3.1   TIME t0 CALCULATIONS.
!              ---------------------

!*       3.1.1  COMPUTE PRE5, PXYB5, PRE5F.

CALL GPHPRE(NPROMA,NFLEVG,KST,KEND,YDVAB,PRE5,PXYB=PXYB5,PRESF=PRE5F)

!*       3.1.2 COMPUTE PRPREF5 AND PXYBDER5.

IF(LVERTFE) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      PRPREF5(JROF,JLEV)=1.0_JPRB/PRE5F(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

CALL GPGRXYB(NPROMA,KST,KEND,NFLEVG,.TRUE.,YDVAB,PRE5L,PRE5M,PXYB5,PXYBDER5)

!*       3.1.3  Compute R, Cp and Kap=R/Cp

CALL GPRCP_QLIRSG(NPROMA,KST,KEND,NFLEVG,PQ=PGFL5(1,1,YQ%MP5),&
 & PCP=PRCP5(1,1,YYTRCP5%M_CP),PR=PRCP5(1,1,YYTRCP5%M_R),PKAP=PRCP5(1,1,YYTRCP5%M_KAP))

!*       3.1.4  Computations of RT and its derivatives

CALL GPRT(LSPRT,NPROMA,KST,KEND,NFLEVG,RD,RV,PRCP5(1,1,YYTRCP5%M_R),&
 & PGMV5(1,1,YT5%MT),PGMV5(1,1,YT5%MTL),&
 & PGMV5(1,1,YT5%MTM),PGFL5(1,1,YQ%MP5L),PGFL5(1,1,YQ%MP5M),&
 & PRT5,PRT5L,PRT5M)

!*       3.1.5  Compute quantities for NH model.

IF (LNHDYN) THEN
  CALL ABOR1(' CPG5_GP: NH model not yet coded for AD and TL ')
  ! not yet coded (see CPG_GP)
  ! call GNHPRE
  ! calculation of PRRED5,ZRTR5
  ! call GNHPREH
  ! call GNHGRPRE
ELSE
  ZNHPPI5(1:NPROMA,1:NFLEVG)=1.0_JPRB
  PRNHPPI5(1:NPROMA,1:NFLEVG)=1.0_JPRB
ENDIF

!*       3.1.6  Solve continuity equation
!               and compute some vertical velocities.

CALL GPCTY(YDVFE,NPROMA,KST,KEND,NFLEVG,LRUBC,YDVAB,YDVETA,&
 & PGMV5(1,1,YT5%MU),PGMV5(1,1,YT5%MV),PGMV5(1,1,YT5%MDIV),ZEVT0,&
 & PXYB5,PRE5L,PRE5M,PRPREF5,PCTY5)

!*       3.1.7  Integrate Hydrostatics
!               Computes the geopotential height "gz" and its gradient.

! * "gz" at full levels and half levels.

PHI5(KST:KEND,NFLEVG)=YDOROG%OROG(KST:KEND)

IF(LNHDYN) THEN
  CALL GPGEO(NPROMA,KST,KEND,NFLEVG,PHI5,PHIF5,&
   & PGMV5(1,1,YT5%MT),PRRED5,PXYB5(1,1,YYTXYB5%M_LNPR),PXYB5(1,1,YYTXYB5%M_ALPH),&
   & YDGEOMETRY%YRVERT_GEOM)
ELSE
  CALL GPGEO(NPROMA,KST,KEND,NFLEVG,PHI5,PHIF5,&
   & PGMV5(1,1,YT5%MT),PRCP5(1,1,YYTRCP5%M_R),PXYB5(1,1,YYTXYB5%M_LNPR),PXYB5(1,1,YYTXYB5%M_ALPH),&
   & YDGEOMETRY%YRVERT_GEOM)
ENDIF

! * "grad gz" at full levels and half levels.
                                                                                
CALL GPGRGEO(YDGEOMETRY,NPROMA,KST,KEND,NFLEVG,&
 & PRT5,PRT5L,PRT5M,&
 & PXYB5(1,1,YYTXYB5%M_LNPR),PXYB5(1,1,YYTXYB5%M_ALPH),PXYBDER5,&
 & POROGL,POROGM,&
 & PHI5FL,PHI5FM,ZGPH5L,ZGPH5M,&
 & LDNHEE=LNHDYN,PRNHPPI=PRNHPPI5,PQCHAL=PQCHA5L,PQCHAM=PQCHA5M,&
 & PNH1L=PNH15L,PNH1M=PNH15M)

!*       3.1.8  Compute half-level winds for non-hydrostatic.

IF(LNHDYN) THEN
  ! ky: if VFE, condition of call to GPHLWI+GPHLUV has to be updated later
  !     when implementing NH+VFE TL and AD codes.
  CALL GPHLWI(YDGEOMETRY%YRDIMV,NPROMA,KST,KEND,PXYB5(1,1,YYTXYB5%M_LNPR),PXYB5(1,1,YYTXYB5%M_ALPH),PUVH5(1,1,YYTHW5%M_WWI))
  CALL GPHLUV(YDGEOMETRY%YRDIMV,NPROMA,KST,KEND,PGMV5(1,1,YT5%MU),PGMV5(1,1,YT5%MV),PUVH5)
ENDIF

!*       3.2   TIME t9 CALCULATIONS.
!              ---------------------

IF (.NOT.LTWOTL) THEN

!*       3.2.1  COMPUTE PRE95, PXYB95, PRE95F.

  IF (LSIMPH.OR.LNHDYN) THEN
    ! note: upstream code (cpgtl, mf_phystl) only uses LNPR, DELP and ALPH from PXYB95
    CALL GPHPRE(NPROMA,NFLEVG,KST,KEND,YDVAB,PRE95,PXYB=PXYB95,PRESF=PRE95F,LRTGR=.FALSE.,LRPP=.FALSE.)
  ENDIF

!*       3.2.3  Compute R, Cp and Kap=Cp/R

  IF (LSIMPH) THEN
    CALL GPRCP_QLIRSG(NPROMA,KST,KEND,NFLEVG,PQ=PGFLT95(1,1,YYTGFLT95%M_Q),&
    & PCP=PRCP95(1,1,YYTRCP95%M_CP),PR=PRCP95(1,1,YYTRCP95%M_R))
  ENDIF

!*       3.2.5  Compute quantities for NH model.

  IF(LNHDYN.AND.LSIMPH) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        PRRED95(JROF,JLEV)=PRCP95(JROF,JLEV,YYTRCP95%M_R)/(1.0_JPRB+PGMVT95(JROF,JLEV,YYTGMVT95%M_SPD)&
         & *SITLAF(JLEV)/PRE95F(JROF,JLEV))
      ENDDO
    ENDDO
  ENDIF

!*       3.2.7  Integrate Hydrostatics
!               Computes the geopotential height "gz" and its gradient.

  IF (LSIMPH) THEN
    PHI95(KST:KEND,NFLEVG)=YDOROG%OROG(KST:KEND)
  ENDIF

  IF(LNHDYN) THEN
    IF(LSIMPH) THEN
      CALL GPGEO(NPROMA,KST,KEND,NFLEVG,&
       & PHI95,PHIF95,PGMVT95(1,1,YYTGMVT95%M_T),PRRED95,PXYB95(1,1,YYTXYB95%M_LNPR),PXYB95(1,1,YYTXYB95%M_ALPH),&
       & YDGEOMETRY%YRVERT_GEOM)
    ENDIF
  ELSE
    IF (LSIMPH) THEN
      CALL GPGEO(NPROMA,KST,KEND,NFLEVG,&
       & PHI95,PHIF95,PGMVT95(1,1,YYTGMVT95%M_T),PRCP95(1,1,YYTRCP95%M_R),&
       & PXYB95(1,1,YYTXYB95%M_LNPR),PXYB95(1,1,YYTXYB95%M_ALPH),&
       & YDGEOMETRY%YRVERT_GEOM)
    ENDIF
  ENDIF

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG5_GP',1,ZHOOK_HANDLE)
END SUBROUTINE CPG5_GP
