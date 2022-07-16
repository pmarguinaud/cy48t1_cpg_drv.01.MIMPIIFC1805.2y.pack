SUBROUTINE LACDYNAD(YDGEOMETRY,YDGMV,YDGMV5,YDEPHY,YDML_GCONF,YDML_DYN,YDPHY,KST,KPROF,PBETADT,PDT,PSLHDA,PSLHDD0,&
 & KIBL,POROGL,POROGM,&
 & PGFL,PRE0,PRDELP0,PCTY0,PATND,KSETTLOFF,PWRL95,&
 & PGMV,PGMVS,PB1,PB2,PGMVT1,PGMVT1S,&
 & PGFL5,PRE5,PRDELP5,PCTY5,PGMV5,PGMV5S)  

!**** *LACDYNAD* Semi-Lagrangian scheme.  (adjoint version)
!                Computation of the t and t-dt useful quantities
!                 at grid-points.

!     Purpose.
!     --------

!          Dynamic non-linear computations in grid-point space
!          for hydrostatic and NH primitive equations and SL scheme.
!          (adjoint code of NH currently not yet coded).

!          Additional remarks:
!          - notation "prehyds" is for hydrostatic surface pressure.
!          - for input and output upper air variables, values are at full levels
!            if no other information is provided.

!          This subroutine fills the semi-Lagrangian buffers to be
!          interpolated.

!**   Interface.
!     ----------
!        *CALL* *LACDYNAD(..)

!        Explicit arguments :
!        --------------------

!        INPUT:      
!          KST     - first element of work.
!          KPROF   - depth of work.
!          PBETADT - BETADT or 0 according to configuration.
!          PDT     - For a leap-frog scheme (three time level scheme):
!                     'dt' at the first time-step, '2 dt' otherwise.
!                    For a 2TL SL scheme: timestep 'dt'.
!          PSLHDA  - Scaling factor of the deformation in f(d) function
!                    (including the model resolution correction)
!          PSLHDD0 - Treshold for deformation tensor enhancement
!          KIBL    - index into YDGSGEOM instance in YDGEOMETRY
!          POROGL  - zonal component of the orography gradient.
!          POROGM  - meridian component of the orography gradient.
!          PWRL95  - trajectory vertical velocity from t-dt


!        INPUT in TL:   ('5' for corresponding trajectory variables)
!          PGFL    - GFL variables at full levels.
!          PRE0    - hydrostatic pressure "prehyd" at half levels at t.
!          PRDELP0 - 1/(pressure depth of layers) at t.
!          PCTY0   - contains vertical velocities, vertical integral of divergence at t.
!          PATND   - adiabatic Lagrangian tendencies.

!        INPUT-OUTPUT in TL:   ('5' for corresponding trajectory variables)
!          PGMV    - GMV variables at t-dt and t.
!          PGMVS   - GMVS variables at t-dt and t.
!          PB1     - "SLBUF1" buffer for interpolations.
!          PB2     - "SLBUF2" buffer.
!          PGMVT1  - GMV variables at t+dt.
!          PGMVT1S - GMVS variables at t+dt.
!          PGFL5 to PGMV5S: trajectory for PGFL to PGMV5S.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           Called by CPG_DYN_AD.

!     Reference.
!     ----------
!             Arpege documentation about semi-lagrangian scheme.

!     Author.
!     -------
!      C. Temperton (ECMWF)
!      Original : 99/10/13.

!     Modifications.
!     --------------
!   Modified 02-01-22 by C. Temperton: case LSPRT=TRUE.
!   O.Spaniel    : 03-04-15 cleaning-a same named entity from modules
!   M.Hamrud      01-Oct-2003 CY28 Cleaning
!   K.Yessad      19-Aug-2004 Correct some intent attributes (PDT,PGFL).
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!   F. Vana  13-Jan-2009: computation of KAPPA for SLHD
!   K. Yessad (Nov 2009): cleanings, DT/Dt now pre-computed in CPG_GP_AD.
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   K. Yessad (Dec 2011): various contributions.
!   M. Diamantakis (Feb 2014): add code for LSETTLSVF option
!   F. Vana  13-Feb-2014: Two more arguments to LATTE_KAPPAAD
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   F. Vana    21-Nov-2017: Option LSLDP_CURV
!  ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD     , ONLY : MODEL_DYNAMICS_TYPE
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOEPHY                 , ONLY : TEPHY
USE YOMPHY                 , ONLY : TPHY
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE YOMGMV                 , ONLY : TGMV
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK
USE YOMDYNA                , ONLY : LSLHD, LSLHD_STATIC
USE INTDYN_MOD             , ONLY : YYTTND, YYTCTY5, YYTCTY0

!  ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV5
TYPE(TEPHY)       ,INTENT(IN)    :: YDEPHY
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
TYPE(TPHY)        ,INTENT(IN)    :: YDPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KSETTLOFF(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBETADT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDA(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLHDD0(YDGEOMETRY%YRDIM%NPROMA)
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWRL95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRE0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRDELP0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCTY0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTCTY0%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PATND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTTND%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%YT1%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1S(YDGEOMETRY%YRDIM%NPROMA,YDGMV%YT1%NDIMS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM5)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE5(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCTY5(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTCTY5%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV5%YT5%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV5S(YDGEOMETRY%YRDIM%NPROMA,YDGMV5%YT5%NDIMS)
!     ------------------------------------------------------------------
! - computed in LASURE:
REAL(KIND=JPRB) :: ZREDIV(YDGEOMETRY%YRDIM%NPROMA)
! - computed in LASSIE or LANHSI:
REAL(KIND=JPRB) :: ZTOD0  (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZGAGT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZGAGT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZBDT(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZSDIV0(YDGEOMETRY%YRDIM%NPROMA)

LOGICAL :: LL2TLFF1

REAL(KIND=JPRB) :: ZBT, ZDTS2, ZESGM, ZESGP
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!  ------------------------------------------------------------------

#include "lassiead.intfb.h"
#include "lasure.intfb.h"
#include "lattesad.intfb.h"
#include "latte_kappaad.intfb.h"
#include "lattexad.intfb.h"
#include "laventad.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LACDYNAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDPTRSLB1=>YDML_DYN%YRPTRSLB1,YDPTRSLB2=>YDML_DYN%YRPTRSLB2, &
 & YGFL=>YDML_GCONF%YGFL)

ASSOCIATE(NDIM=>YGFL%NDIM, NDIM5=>YGFL%NDIM5, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT1=>YDGMV%YT1, &
 & YT5=>YDGMV5%YT5, &
 & NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS:
!              ----------------------------

! Remark: T/L of LASURE is not needed

CALL LASURE(YDGEOMETRY,YDEPHY,YDML_DYN%YRDYN,YDML_DYN%YREDYN,YDPHY,KST,KPROF,PBETADT,PDT,KIBL,&
 & ZDTS2,ZBT,LL2TLFF1,ZBDT,ZREDIV,ZESGP,ZESGM)

ZTOD0(:,:)=0.0_JPRB
ZGAGT0L(:,:)=0.0_JPRB
ZGAGT0M(:,:)=0.0_JPRB
ZSDIV0(:)=0.0_JPRB

!     ------------------------------------------------------------------

!*       6.    COMPUTATION OF "KAPPA" 
!              ----------------------

IF (LSLHD.AND.(.NOT.LSLHD_STATIC)) THEN
  CALL LATTE_KAPPAAD(YDGEOMETRY,YDGMV,YDGMV5,YDML_DYN%YRDYN,YDPTRSLB2,KST,KPROF,ZDTS2,PSLHDA,PSLHDD0,PGMV,PGMVS,PGMV5,PGMV5S,PB2)
ENDIF

!     ------------------------------------------------------------------

!*       5.    COMPUTATION OF THE 2D-EQUATIONS RIGHT-HAND SIDE TERMS.
!              ------------------------------------------------------

CALL LATTESAD(YDGEOMETRY,YDGMV,YDGMV5,YDML_GCONF%YRRIP,YDML_DYN,KST,KPROF,ZDTS2,ZBDT,ZESGP,ZESGM,&
 & POROGL,POROGM,ZSDIV0,PCTY0(1,0,YYTCTY0%M_PSDVBC),PRE0,PGMVS,&
 & PGMV,PGMVT1S,PB1,PB2,&
 & PCTY5(1,0,YYTCTY5%M_PSDVBC),PRE5,PGMV5,PGMV5S)

!     ------------------------------------------------------------------

!*       4.    COMPUTATION OF THE 3D-EQUATIONS RIGHT-HAND SIDE TERMS.
!              ------------------------------------------------------

CALL LATTEXAD(YDGEOMETRY,YDGMV,YDML_GCONF,YDML_DYN,KST,KPROF,ZDTS2,ZBT,ZBDT,ZESGP,ZESGM,&
 & KIBL,POROGL,POROGM,&
 & ZGAGT0L,ZGAGT0M,ZTOD0,PRDELP0,PCTY0(1,0,YYTCTY0%M_EVEL),&
 & PGFL,PATND,PGMV,PGMVT1,PB1,PB2,&
 & PRDELP5,PCTY5(1,0,YYTCTY5%M_EVEL))

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF THE WIND COMPONENTS NECESSARY FOR SL TRAJECTORY.
!              ---------------------------------------------------------------

CALL LAVENTAD(YDGEOMETRY,YDGMV,YDML_GCONF%YRRIP,YDML_DYN,KST,KPROF,KIBL,LL2TLFF1,ZDTS2,&
 & PRDELP0,PCTY0(1,0,YYTCTY0%M_EVEL),PATND,KSETTLOFF,&
 & PWRL95,PGMV,PB1,PB2,PRDELP5,PCTY5(1,0,YYTCTY5%M_EVEL))

!     ------------------------------------------------------------------

!*       2.    COMPUTATION OF THE LINEAR TERMS FOR SEMI-IMPLICIT SCHEME.
!              ---------------------------------------------------------

CALL LASSIEAD(YDGEOMETRY,YDGMV,YDGMV5,YGFL,YDML_DYN%YRDYN,KST,KPROF,YDGSGEOM%RCORI,PGMV,PGMVS,PGFL,&
 & ZSDIV0,ZTOD0,ZGAGT0L,ZGAGT0M,PGMV5,PGFL5)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LACDYNAD',1,ZHOOK_HANDLE)
END SUBROUTINE LACDYNAD
