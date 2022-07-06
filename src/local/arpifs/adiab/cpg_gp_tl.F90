#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE CPG_GP_TL(YDGEOMETRY,YDGMV,YDGMV5,&
 !---------------------------------------------------------------------
 ! - INPUT .
 & YDEPHY,YDML_GCONF,YDML_DYN,YDSIMPHL,KST,KEND,LDLFSTEP,PTE,KIBL,&
 & PGFL,&
 !---------------------------------------------------------------------
 ! - INPUT/OUTPUT .
 & PGMV,PGMVS,&
 !---------------------------------------------------------------------
 ! - OUTPUT .
 & PRE0,PRE0L,PRE0M,PRE0F,PXYB0,&
 & PHI0,PHIF0,PRCP0,PCTY0,PKENE0,PRT0L,PRT0M,&
 & PRE9,PRE9F,PXYB9,PHI9,PHIF9,PRCP9,&
 & PGPPC,PB2,PGMVT1,PGMVT1S,PGFLT1,PATND,&
 !---------------------------------------------------------------------
 ! - TRAJECTORY (INPUT) .
 & PGFL5,PGMV5,PGMV5S,&
 & PRE5,PRE5L,PRE5M,PRE5F,PXYB5,PRCP5,&
 & PCTY5,PRT5,PRT5L,PRT5M,&
 & PRPREF5,&
 & PXYBDER5,PRNHPPI5,PHI5FL,PHI5FM,&
 & PGMVT95,PRE95,PXYB95,PRCP95&
 & )

!**** *CPG_GP_TL* - Tan. lin. grid point calculations.
!                   initial part of not lagged grid-point calculations.

!     Purpose.
!     --------
!           Tan. lin. grid point calculations.
!           initial part of not lagged grid-point calculations.
!           - get data in buffers.
!           - multiply p-order horizontal derivatives by M**p.
!           - second part of the temporal filter.
!           - calls some GP... routines to initialise some auxiliary variables.
!           - sets-up and PB2.

!           Remark: the way of coding is the most consistent as possible
!           with the way of coding of CPG_GP. Some items present in CPG_GP
!           and not coded in CPG_GP_TL are mentioned by comments.

!**   Interface.
!     ----------
!        *CALL* *CPG_GP_TL(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        KST       : first element of work.
!        KEND      : last element of work.
!        LDLFSTEP  : .T.: first time-step (T if NSTEP <= NSTART)
!        PTE       : 1. or 0. according to different configurations.
!        KIBL      : index into YRGSGEOM instance in YDGEOMETRY
!        PGFL      : unified_treatment grid-point fields at t

!     INPUT/OUTPUT:
!     -------------
!        PGMV      : upper air GMV variables at time t and t-dt.
!        PGMVS     : surface GMV variables at time t and t-dt.

!     OUTPUT:
!     -------
!        PRE0      : hydrostatic pressure "prehyd" at half levels at time t.
!        PRE0L     : zonal component of "grad prehyds" at t.
!        PRE0M     : meridian component of "grad prehyds" at t.
!        PRE0F     : hydrostatic pressure "prehyd" at full levels at time t.
!        PXYB0     : contains pressure depth, "delta", "alpha" at t.
!        PUVH0     : horizontal wind at time t at half levels.
!        PHI0      : geopotential height "gz" at t at half levels.
!        PHIF0     : geopotential height "gz" at t at full levels.
!        PRCP0     : contains "cp", "R" and "Kap=R/Cp" at t.
!        PCTY0     : contains vertical velocities, vertical integral of divergence at t.
!        PKENE0    : kinetic energy at t.
!        PRT0L     : zonal component of "grad RT" at full levels.
!        PRT0M     : meridian component of "grad RT" at full levels.
!        PRE9      : hydrostatic pressure "prehyd" at half levels at t-dt.
!        PRE9F     : hydrostatic pressure "prehyd" at full levels at time t-dt.
!        PXYB9     : contains pressure depth, "delta", "alpha" at t-dt.
!        PHI9      : geopotential height "gz" at t-dt at half levels.
!        PHIF9     : geopotential height "gz" at t-dt at full levels.
!        PRCP9     : contains "cp", "R" and "Kap=R/Cp" at t-dt.
!        PGPPC     : buffer containing some fields used in the PC schemes.
!        PB2       : "SLB2" buffer.
!        PGMVT1    : upper air GMV variables at t+dt.
!        PGMVT1S   : surface GMV variables at t+dt.
!        PGFLT1    : GFL variables at t+dt.
!        PATND     : adiabatic Lagrangian tendencies.

!     TRAJECTORY (INPUT):
!     -------------------
!        PGFL5     : unified_treatment grid-point fields at t.
!        PGMV5     : upper air GMV variables at time t and t-dt.
!        PGMV5S    : surface GMV variables at time t and t-dt.
!        PRE5      : hydrostatic pressure "prehyd" at half levels at time t.
!        PRE5L     : zonal component of "grad prehyds" at t.
!        PRE5M     : meridian component of "grad prehyds" at t.
!        PRE5F     : hydrostatic pressure "prehyd" at full levels at time t.
!        PXYB5     : contains pressure depth, "delta", "alpha" at t.
!        PRCP5     : contains "cp", "R" and "Kap=R/Cp" at t.
!        PCTY5     : contains vertical velocities, vertical integral of divergence at t.
!        PRT5      : RT at full levels at t.
!        PRT5L     : zonal component of "grad RT" at full levels at t.
!        PRT5M     : meridian component of "grad RT" at full levels at t.
!        PRPREF5   : 1/PRE5F.
!        PXYBDER5  : cf. PXYBDER in GPGRXYB (time t).
!        PRNHPPI5  : "prehyd/pre" at full levels (time t).
!        PHI5FL    : zonal comp grad(gz) at full levels
!        PHI5FM    : merid comp grad(gz) at full levels
!        PGMVT95   : GMV at time t-dt at full levels.
!        PRE95     : hydrostatic pressure "prehyd" at half levels at t-dt.
!        PXYB95    : contains pressure depth, "delta", "alpha" at t-dt.
!        PRCP95    : contains "cp", "R" and "Kap=R/Cp" at t-dt.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        K. Yessad, after CPGTL (old part 4 of CPGTL)
!        Original     : 12 Jul 2004.

! Modifications
! -------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2008): remove dummy CDLOCK + cleanings
!   K. Yessad (Aug 2009): remove LPC_OLD in TL and AD codes.
!   K. Yessad (Nov 2009): cleanings, DT/Dt now pre-computed in CPG_GP_TL.
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   K. Yessad (Nov 2011): various contributions.
!   K. Yessad (Nov 2012): simplify testings.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   K. Yessad (Feb 2018): remove "NYC-NH" commented lines.
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
! End Modifications
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOEPHY       , ONLY : TEPHY
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE SC2PRG_MOD   , ONLY : SC2PRG
USE YOMCT0       , ONLY : LTWOTL, LSPRT
USE YOMCST       , ONLY : RD, RV
USE YOMDYNA      , ONLY : LGRADSP, LNHX, LGWADV
USE YOMCT0       , ONLY : LSLAG, LTWOTL, LNHDYN
USE YOMCVER      , ONLY : LVERTFE, LVFE_GW
USE YOMSIMPHL    , ONLY : TSIMPHL
USE INTDYN_MOD   , ONLY : YYTTND, YYTCTY5, YYTCTY0, YYTXYBDER0, YYTXYBDER5,&
 & YYTRCP5, YYTRCP95, YYTRCP0, YYTRCP9, YYTXYB5, YYTXYB95, YYTXYB0, YYTXYB9,&
 & YYTGMVT95

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)            ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)                ,INTENT(INOUT) :: YDGMV
TYPE(TGMV)                ,INTENT(INOUT) :: YDGMV5
TYPE(TEPHY)               ,INTENT(INOUT) :: YDEPHY
TYPE(MODEL_DYNAMICS_TYPE) ,INTENT(INOUT) :: YDML_DYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TSIMPHL)             ,INTENT(INOUT) :: YDSIMPHL
INTEGER(KIND=JPIM)        ,INTENT(IN)    :: KST 
INTEGER(KIND=JPIM)        ,INTENT(IN)    :: KEND 
LOGICAL                   ,INTENT(IN)    :: LDLFSTEP 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PTE 
INTEGER(KIND=JPIM)        ,INTENT(IN)    :: KIBL
REAL(KIND=JPRB)           ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB)           ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV) 
REAL(KIND=JPRB)           ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PRE0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PRE0L(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PRE0M(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PRE0F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PXYB0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB0%NDIM) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PHI0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PHIF0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PRCP0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP0%NDIM)
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PCTY0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTCTY0%NDIM)
REAL(KIND=JPRB)           ,INTENT(INOUT) :: PKENE0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PRT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PRT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PRE9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PRE9F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PXYB9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB9%NDIM) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PHI9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PHIF9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PRCP9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP9%NDIM)
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PGPPC(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRGPPC%NFGPPC+1)
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PGMVT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%YT1%NDIM) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PGMVT1S(YDGEOMETRY%YRDIM%NPROMA,YDGMV%YT1%NDIMS) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PGFLT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM1) 
REAL(KIND=JPRB)           ,INTENT(OUT)   :: PATND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTTND%NDIM)
REAL(KIND=JPRB)           ,INTENT(IN)    :: PGFL5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM5)
REAL(KIND=JPRB)           ,INTENT(IN)    :: PGMV5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV5%YT5%NDIM)
REAL(KIND=JPRB)           ,INTENT(IN)    :: PGMV5S(YDGEOMETRY%YRDIM%NPROMA,YDGMV5%YT5%NDIMS)
REAL(KIND=JPRB)           ,INTENT(IN)    :: PRE5(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PRE5L(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PRE5M(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PRE5F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PXYB5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB5%NDIM) 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PRCP5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP5%NDIM)
REAL(KIND=JPRB)           ,INTENT(IN)    :: PCTY5(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTCTY5%NDIM)
REAL(KIND=JPRB)           ,INTENT(IN)    :: PRT5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PRT5L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PRT5M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PRPREF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)           ,INTENT(IN)    :: PXYBDER5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYBDER5%NDIM)
REAL(KIND=JPRB)           ,INTENT(IN)    :: PRNHPPI5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)           ,INTENT(IN)    :: PHI5FL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)           ,INTENT(IN)    :: PHI5FM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)           ,INTENT(IN)    :: PGMVT95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTGMVT95%NDIM)
REAL(KIND=JPRB)           ,INTENT(IN)    :: PRE95(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PXYB95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB95%NDIM) 
REAL(KIND=JPRB)           ,INTENT(IN)    :: PRCP95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP95%NDIM)
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZGPHL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZGPHM(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZRT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZRPRE0F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDUM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZXYBDER0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYBDER0%NDIM)
REAL(KIND=JPRB) :: ZPSGRTL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZSGRTL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPSGRTM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZSGRTM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: Z_PHI0FL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_PHI0FM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: JLEV, JROF, ILVVEL, IFLAG

LOGICAL :: LLSTR
LOGICAL :: LLNHDYN
REAL(KIND=JPRB) :: ZEPS

REAL(KIND=JPRB), POINTER :: ZT0DIV(:,:), ZT0SGRTL(:,:), ZT0SGRTM(:,:), ZT0SP(:)
REAL(KIND=JPRB), POINTER :: ZT0SPL(:), ZT0SPM(:), ZT0T(:,:), ZT0TL(:,:)
REAL(KIND=JPRB), POINTER :: ZT0TM(:,:), ZT0U(:,:), ZT0V(:,:), ZT1SP(:)
REAL(KIND=JPRB), POINTER :: ZT1SPD(:,:), ZT1SVD(:,:), ZT1T(:,:), ZT1U(:,:)
REAL(KIND=JPRB), POINTER :: ZT1V(:,:), ZT9SGRTL(:,:), ZT9SGRTM(:,:), ZT9SP(:)
REAL(KIND=JPRB), POINTER :: ZT9SPD(:,:), ZT9SVD(:,:), ZT9T(:,:), ZT9U(:,:)
REAL(KIND=JPRB), POINTER :: ZT9V(:,:)
REAL(KIND=JPRB), POINTER :: ZPQ(:,:), ZPQ9(:,:), ZPQL(:,:), ZPQM(:,:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gpctytl.intfb.h"
#include "gpgeotl.intfb.h"
#include "gpgrgeotl.intfb.h"
#include "gpgrptl.intfb.h"
#include "gpgrxybtl.intfb.h"
#include "gpinislb.intfb.h"
#include "gpmpfc.intfb.h"
#include "gphpretl.intfb.h"
#include "gprcptl.intfb.h"
#include "gprttl.intfb.h"
#include "gp_spvtl.intfb.h"
#include "gptf2.intfb.h"
#include "gp_tndlagadiab_uv_tl.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPG_GP_TL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE, YDSTA=>YDGEOMETRY%YRSTA, &
 & YDLAP=>YDGEOMETRY%YRLAP, YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDVSPLIP=>YDGEOMETRY%YRVSPLIP,  &
 & YDVSLETA=>YDGEOMETRY%YRVSLETA, YDHSLMER=>YDGEOMETRY%YRHSLMER, YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL),  &
 & YDSPGEOM=>YDGEOMETRY%YSPGEOM, YDDYN=>YDML_DYN%YRDYN,YDPTRSLB2=>YDML_DYN%YRPTRSLB2, &
 & YGFL=>YDML_GCONF%YGFL,YDPTRGPPC=>YDML_DYN%YRPTRGPPC)

ASSOCIATE(NDIM=>YGFL%NDIM, NDIM1=>YGFL%NDIM1, NDIM5=>YGFL%NDIM5, YQ=>YGFL%YQ, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, LSIMPH=>YDSIMPHL%LSIMPH, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
 & RSTRET=>YDGEM%RSTRET, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0, &
 & YT1=>YDGMV%YT1, YT9=>YDGMV%YT9, &
 & YT5=>YDGMV5%YT5, &
 & NFGPPC=>YDPTRGPPC%NFGPPC, &
 & MSLB2VVEL=>YDPTRSLB2%MSLB2VVEL, NFLDSLB2=>YDPTRSLB2%NFLDSLB2)

CALL SC2PRG(YQ%MP     ,PGFL     ,ZPQ)
CALL SC2PRG(YQ%MP9    ,PGFL     ,ZPQ9)
CALL SC2PRG(YQ%MPL    ,PGFL     ,ZPQL)
CALL SC2PRG(YQ%MPM    ,PGFL     ,ZPQM)
CALL SC2PRG(YT0%MDIV  ,PGMV    ,ZT0DIV)
CALL SC2PRG(YT0%MSGRTL,PGMV    ,ZT0SGRTL)
CALL SC2PRG(YT0%MSGRTM,PGMV    ,ZT0SGRTM)
CALL SC2PRG(YT0%MSP   ,PGMVS   ,ZT0SP)
CALL SC2PRG(YT0%MSPL  ,PGMVS   ,ZT0SPL)
CALL SC2PRG(YT0%MSPM  ,PGMVS   ,ZT0SPM)
CALL SC2PRG(YT0%MT    ,PGMV    ,ZT0T)
CALL SC2PRG(YT0%MTL   ,PGMV    ,ZT0TL)
CALL SC2PRG(YT0%MTM   ,PGMV    ,ZT0TM)
CALL SC2PRG(YT0%MU    ,PGMV    ,ZT0U)
CALL SC2PRG(YT0%MV    ,PGMV    ,ZT0V)
CALL SC2PRG(YT1%MSP   ,PGMVT1S ,ZT1SP)
CALL SC2PRG(YT1%MSPD  ,PGMVT1  ,ZT1SPD)
CALL SC2PRG(YT1%MSVD  ,PGMVT1  ,ZT1SVD)
CALL SC2PRG(YT1%MT    ,PGMVT1  ,ZT1T)
CALL SC2PRG(YT1%MU    ,PGMVT1  ,ZT1U)
CALL SC2PRG(YT1%MV    ,PGMVT1  ,ZT1V)
CALL SC2PRG(YT9%MSGRTL,PGMV    ,ZT9SGRTL)
CALL SC2PRG(YT9%MSGRTM,PGMV    ,ZT9SGRTM)
CALL SC2PRG(YT9%MSP   ,PGMVS   ,ZT9SP)
CALL SC2PRG(YT9%MSPD  ,PGMV    ,ZT9SPD)
CALL SC2PRG(YT9%MSVD  ,PGMV    ,ZT9SVD)
CALL SC2PRG(YT9%MT    ,PGMV    ,ZT9T)
CALL SC2PRG(YT9%MU    ,PGMV    ,ZT9U)
CALL SC2PRG(YT9%MV    ,PGMV    ,ZT9V)
!     ------------------------------------------------------------------

!*       0.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

ZEPS=100.0_JPRB*TINY(1.0_JPRB)
LLSTR=(ABS(RSTRET-1.0_JPRB)>ZEPS)
LLNHDYN=.FALSE. ! TL of NH not coded

!     ------------------------------------------------------------------

!*       1.    INTERFACE TO GLOBAL ARRAYS/WORK FILES.
!              --------------------------------------

PGPPC(:,:)=0.0_JPRB

!*       1.1  TENDENCY COUPLING FOR SURFACE PRESSURE.

! Not coded in the TL code (CALL GPTENCTL).

!*       1.4  PART 2 OF TIME FILTER.

IF(NCURRENT_ITER == 0) THEN
  CALL GPTF2(YDGEOMETRY,YDGMV,&
   ! --- INPUT ---------------------------------------------------------
   & YDML_GCONF,YDDYN,KST,KEND,LDLFSTEP,&
   ! --- INPUT for P.T0.; INPUT/OUTPUT for P.T9. -----------------------
   & PGMV,PGMVS,PGFL)
ENDIF

!*       1.5  MAP FACTOR AND SURFACE PRESSURE VARIABLES.

! * Geographical derivatives of orography:
!   calculation has already been done in CPG5_GP, so no need to do it here (contrary to
!   what is done in CPG_GP).

IF(LLSTR) THEN
  IFLAG=0
  CALL GPMPFC(YDGMV,YDML_GCONF,YDDYN,NPROMA,NFLEVG,KST,KEND,IFLAG,YDGSGEOM%GM,PGMV,PGMVS,PGFL)
ENDIF

CALL GP_SPVTL(YDGEOMETRY,YDDYN,YDSIMPHL,KST,KEND,&
 & ZT0SP,ZT0SPL,ZT0SPM,ZT9SP,&
 & PRE0,PRE0L,PRE0M,PRE9,&
 & PGMV5S(1,YT5%MSPL),PGMV5S(1,YT5%MSPM),PRE5,PRE95)

!*       1.6  GET SURFACE VARIABLES.

! Not coded in the TL code 

!*       1.7  SEMI-LAGRANGIAN MASS FIX.

! Not coded in the TL code.

!*       1.8  OZONE PHYSICO-CHEMICAL PROPERTIES AND SUBGRID SURFACE TEMPERATURE

! Not coded in the TL code (CALL GPINOZSTDMTL).

!*       1.9  NUDGING.

! Not coded in the TL code (CALL UPDSSTTL).

!     ------------------------------------------------------------------

!*       3.    INITIALISE AUXILIARY VARIABLES.
!              -------------------------------

!*       3.1   TIME t0 CALCULATIONS.

!*       3.1.1 COMPUTE PRE0, PXYB0, PRE0F.

CALL GPHPRETL(NPROMA,NFLEVG,KST,KEND,YDVAB,PRE0,PRE5,PXYB=PXYB0,PXYB5=PXYB5,PRESF=PRE0F)

!*       3.1.2 COMPUTE ZRPRE0F AND ZXYBDER0.

! * Case of finite elements:
IF(LVERTFE) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      ZRPRE0F(JROF,JLEV)=-PRE0F(JROF,JLEV)*PRPREF5(JROF,JLEV)**2
    ENDDO
  ENDDO
ENDIF

CALL GPGRXYBTL(NPROMA,KST,KEND,NFLEVG,YDVAB,PRE0L,PRE0M,PXYB0,ZXYBDER0,&
 & PRE5L,PRE5M,PXYB5,PXYBDER5)

!*       3.1.3 COMPUTE "R", "Cp" AND "Kap=R/Cp".

CALL GPRCPTL(NPROMA,KST,KEND,NFLEVG,&
 & PQ=ZPQ,PCP=PRCP0(1,1,YYTRCP0%M_CP),PR=PRCP0(1,1,YYTRCP0%M_R),PKAP=PRCP0(1,1,YYTRCP0%M_KAP),&
 & PCP5=PRCP5(1,1,YYTRCP5%M_CP),PR5=PRCP5(1,1,YYTRCP5%M_R))

!*       3.1.4 COMPUTES "RT" AND ITS GRADIENT.

CALL GPRTTL(LSPRT,NPROMA,KST,KEND,NFLEVG,RD,RV,&
 & PRCP5(1,1,YYTRCP5%M_R),PGMV5(1,1,YT5%MT),PGMV5(1,1,YT5%MTL),PGMV5(1,1,YT5%MTM),&
 & PGFL5(1,1,YQ%MP5L),PGFL5(1,1,YQ%MP5M),&
 & PRCP0(1,1,YYTRCP0%M_R),ZT0T,ZT0TL,ZT0TM,&
 & ZPQL,ZPQM,ZRT0,PRT0L,PRT0M)

!*       3.1.7 COMPUTATION OF SOME VERTICAL VELOCITIES.
!*             Solve continuity equation
!              Computation of the vertical velocities
!              "etapt (d prehyd / d eta)" and "omega/prehyd".

CALL GPCTYTL(YDVFE,NPROMA,KST,KEND,NFLEVG,YDVAB,YDVETA,&
 & ZT0U,ZT0V,ZT0DIV,&
 & PXYB0,PRE0L,PRE0M,ZRPRE0F,PCTY0,&
 & PGMV5(1,1,YT5%MU),PGMV5(1,1,YT5%MV),PGMV5(1,1,YT5%MDIV),&
 & PXYB5,PRE5L,PRE5M,PCTY5(1,0,YYTCTY5%M_PSDIV),PRPREF5)

!*       3.1.8 COMPUTES THE GEOPOTENTIAL HEIGHT "gz" AND ITS GRADIENT.

! * "gz" at full levels and half levels.

PHI0(KST:KEND,NFLEVG)=0.0_JPRB
CALL GPGEOTL(NPROMA,KST,KEND,NFLEVG,&
 & PHI0,PHIF0,ZT0T,PRCP0(1,1,YYTRCP0%M_R),PXYB0(1,1,YYTXYB0%M_LNPR),PXYB0(1,1,YYTXYB0%M_ALPH),&
 & PGMV5(1,1,YT5%MT),PRCP5(1,1,YYTRCP5%M_R),PXYB5(1,1,YYTXYB5%M_LNPR),PXYB5(1,1,YYTXYB5%M_ALPH),&
 & YDGEOMETRY%YRVERT_GEOM)  

! * "grad gz" at full levels and half levels

CALL GPGRGEOTL(YDGEOMETRY,LLNHDYN,&
 & NPROMA,KST,KEND,NFLEVG,&
 & ZT0T,PRCP0(1,1,YYTRCP0%M_R),PRT0L,PRT0M,&
 & PXYB0(1,1,YYTXYB0%M_LNPR),PXYB0(1,1,YYTXYB0%M_ALPH),ZXYBDER0,&
 & Z_PHI0FL,Z_PHI0FM,ZGPHL,ZGPHM,&
 & PRNHPPI5,PGMV5(1,1,YT5%MT),PRCP5(1,1,YYTRCP5%M_R),PRT5L,PRT5M,&
 & PXYB5(1,1,YYTXYB5%M_LNPR),PXYB5(1,1,YYTXYB5%M_ALPH),PXYBDER5 )

!*       3.1.9 COMPUTES KINETIC ENERGY.

DO JLEV=1,NFLEVG
  DO JROF=KST,KEND
    PKENE0(JROF,JLEV)=(ZT0U(JROF,JLEV)*PGMV5(JROF,JLEV,YT5%MU)&
     & +ZT0V(JROF,JLEV)*PGMV5(JROF,JLEV,YT5%MV))  
  ENDDO
ENDDO

!*       3.1.12 COMPUTES THE PRESSURE FORCE FOR THE RHS OF MOMENTUM EQN.

CALL GPGRPTL(&
 ! --- INPUT ----------------------------------------------------------------
 & NPROMA,&
 & KST,KEND,NFLEVG,&
 & ZRT0,PRT0L,PRT0M,PRE0L,PRE0M,PXYB0,ZXYBDER0,&
 & ZGPHL,ZGPHM,Z_PHI0FL,Z_PHI0FM,&
 ! --- OUTPUT ---------------------------------------------------------------
 & ZPSGRTL,ZPSGRTM,&
 ! --- INPUT-TRAJECTORY -----------------------------------------------------
 & PRT5,PRT5L,PRT5M,PRE5L,PRE5M,PXYB5,PXYBDER5,&
 & PHI5FL,PHI5FM&
 & )

! store for next time-step
IF( LGRADSP ) THEN

  ! adjust field with filter result
  ZSGRTL(KST:KEND,1:NFLEVG)=ZPSGRTL(KST:KEND,1:NFLEVG) - (&
   & ZT9SGRTL(KST:KEND,1:NFLEVG) - ZT0SGRTL(KST:KEND,1:NFLEVG))
  ZSGRTM(KST:KEND,1:NFLEVG)=ZPSGRTM(KST:KEND,1:NFLEVG) - (&
   & ZT9SGRTM(KST:KEND,1:NFLEVG) - ZT0SGRTM(KST:KEND,1:NFLEVG))

  ! store for next time-step unfiltered fields
  ZT0SGRTL(KST:KEND,1:NFLEVG)=ZPSGRTL(KST:KEND,1:NFLEVG)
  ZT0SGRTM(KST:KEND,1:NFLEVG)=ZPSGRTM(KST:KEND,1:NFLEVG)

  ! set new values for current timestep
  ZPSGRTL(KST:KEND,1:NFLEVG)=ZSGRTL(KST:KEND,1:NFLEVG)
  ZPSGRTM(KST:KEND,1:NFLEVG)=ZSGRTM(KST:KEND,1:NFLEVG)

ENDIF

!*     3.1.20a COMPUTES [D V/Dt].

! Pressure gradient term + explicit Coriolis + Rayleigh friction
CALL GP_TNDLAGADIAB_UV_TL(YDGEOMETRY,YDGMV,YDEPHY,YDDYN,KST,KEND,YDGSGEOM%RCORI,YDGSGEOM%GNORDL,YDGSGEOM%GNORDM,&
 & ZPSGRTL,ZPSGRTM,PGMV,&
 & PATND(1,1,YYTTND%M_TNDU),PATND(1,1,YYTTND%M_TNDV),&
 & PATND(1,1,YYTTND%M_TNDU_NOC),PATND(1,1,YYTTND%M_TNDV_NOC))

!*     3.1.24 COMPUTES [DT/Dt].

DO JLEV=1,NFLEVG
  DO JROF=KST,KEND
    PATND(JROF,JLEV,YYTTND%M_TNDT)=&
     & PRCP0(JROF,JLEV,YYTRCP0%M_KAP)*PGMV5(JROF,JLEV,YT5%MT)*PCTY5(JROF,JLEV,YYTCTY5%M_VVEL)&
     & +PRCP5(JROF,JLEV,YYTRCP5%M_KAP)*ZT0T(JROF,JLEV)*PCTY5(JROF,JLEV,YYTCTY5%M_VVEL)&
     & +PRCP5(JROF,JLEV,YYTRCP5%M_KAP)*PGMV5(JROF,JLEV,YT5%MT)*PCTY0(JROF,JLEV,YYTCTY0%M_VVEL)
  ENDDO
ENDDO

!*       3.2   TIME t9 CALCULATIONS.

IF(.NOT.LTWOTL .AND. NCURRENT_ITER == 0) THEN

!*       3.2.1 COMPUTE PRE9, PXYB9, PRE9F.

  IF (LSIMPH) THEN
    CALL GPHPRETL(NPROMA,NFLEVG,KST,KEND,YDVAB,PRE9,PRE95,PXYB=PXYB9,PXYB5=PXYB95,PRESF=PRE9F)
  ENDIF

  ! * Case of finite elements: not done in the TL code (ZRPRE9F is useless).

!*       3.2.3 COMPUTE "R", "Cp" AND "Kap=R/Cp".

  IF(LSIMPH) THEN
    CALL GPRCPTL(NPROMA,KST,KEND,NFLEVG,PQ=ZPQ9,PCP=PRCP9(1,1,YYTRCP9%M_CP),PR=PRCP9(1,1,YYTRCP9%M_R))
  ENDIF

!*       3.2.8 COMPUTES THE GEOPOTENTIAL HEIGHT "gz" AND ITS GRADIENT.

  ! * "gz" at full levels and half levels.

  IF (LSIMPH) THEN
    PHI9(KST:KEND,NFLEVG)=0.0_JPRB
    CALL GPGEOTL(NPROMA,KST,KEND,NFLEVG,PHI9,PHIF9,&
     & ZT9T,PRCP9(1,1,YYTRCP9%M_R),PXYB9(1,1,YYTXYB9%M_LNPR),PXYB9(1,1,YYTXYB9%M_ALPH),&
     & PGMVT95(1,1,YYTGMVT95%M_T),PRCP95(1,1,YYTRCP95%M_R),PXYB95(1,1,YYTXYB95%M_LNPR),PXYB95(1,1,YYTXYB95%M_ALPH),&
     & YDGEOMETRY%YRVERT_GEOM)  
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       5.    TIME DIMENSION.
!              ---------------

CALL GPINISLB(&
 ! --- INPUT ---------------------------------------------------------
 & LGWADV, LNHDYN, LNHX, LSLAG, LTWOTL, LVERTFE, LVFE_GW, &
 & YDGEOMETRY,YGFL,YDDYN,YDPTRSLB2,KST,KEND,PTE,.TRUE.,&
 & PGFL,&
 & ZT9U,ZT9V,ZT9T,&
 & ZT9SPD,ZT9SVD,ZDUM,&
 & ZT9SP,&
 & ZDUM,ZDUM,&
 ! --- OUTPUT --------------------------------------------------------
 & ZT1U,ZT1V,ZT1T,PGFLT1,&
 & ZT1SPD,ZT1SVD,ZDUM,&
 & ZT1SP,&
 & PB2,ZDUM,ZDUM,ZDUM,ZDUM)

! * Compute perturbation of vertical velocity:
!   (not coded in GPINISLB because non linear, allows to avoid to code
!   GPINISLBTL).
DO JLEV=1,NFLEVG
  ILVVEL=MSLB2VVEL+JLEV-1
  DO JROF=KST,KEND
    PB2(JROF,ILVVEL)=PCTY5(JROF,JLEV,YYTCTY5%M_VVEL)*PRE0F(JROF,JLEV)&
     & +PRE5F(JROF,JLEV)*PCTY0(JROF,JLEV,YYTCTY0%M_VVEL)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG_GP_TL',1,ZHOOK_HANDLE)
END SUBROUTINE CPG_GP_TL
