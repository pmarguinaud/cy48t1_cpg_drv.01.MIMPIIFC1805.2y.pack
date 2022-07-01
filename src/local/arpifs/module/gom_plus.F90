MODULE GOM_PLUS
!
! Alan Geer   1-June-2015     ECMWF
!
! Purpose.
! --------
!
! Maintains and fills a data structure containing the model fields interpolated to
! observation locations, arranged in KSET order, and with additional re-generated
! model fields such as geopotential, full and half pressure. This makes 
! extensive use of forecast model code like GPGEO - now we need to work out how
! to OOPSify this!
!
! Replaces code previously known as "the preint routines".
!
! AJGDB could be made more automatic in the future (using the "model fields" concept
! of a few years ago) but for now it is quite hand-coded.
!
! AJGDB2 we might aspire to change the variable names to something more meaningful
! and/or consistent with the names in the GOMS, but for now these names are, to
! prevent confusion, identical to the old variables used around hop and elsewhere
!
! (GOM**6) Copy GOM variables from/to hop*.F90 "packets" 
! 1) First follow the instructions at the top of gom_mod.F90 to add the new variable in the GOMs.
! 2) Add a pointer of the appropriate name, shape and type in the definition of type_gom_plus
! 3) In gom_plus_alloc, point the new variable into a storage array of the appropriate shape, at the next available index
! 4) Also in gom_plus_alloc, increase the size of the storage area as appropriate (e.g. insfc, inua1 etc.)
! 5) Nullify the new pointer in gom_plus_destroy
! 6) Extract the variable from the GOMs in gom_plus_fill
! 7) Put back to the GOMs in gom_plus_fill_ad
!
! Possibly out of date but it might help: (see gom_mod for up-to-date accurate info)
!
!        WHERE KDLEN   = FIRST DIMENSION (OBSERVATIONS)
!              KSET    = NUMBER (NUMERO) OF THE OBS SET TO BE TREATED
!              POROG   = OROGRAPHY AT OBS. PONTS
!              PLS5    = Land/Sea Mask
!              PCI5    = Sea-ice fraction
!              PAL5    = Albedo
!              PTS5    = Surface skin temperature
!              PWS5    = Soil wetness
!              PSN5    = Snow depth
!              PU10N5  = U-COMP OF 10M NEUTRAL WIND
!              PV10N5  = V-COMP OF 10M NEUTRAL WIND
!              PSST5   = SURFACE TEMPERATURE OVER SEA (duplicated obs)
!              PU05    = U-COMP OF OCEAN CURRENT
!              PV05    = V-COMP OF OCEAN CURRENT
!              PCHN5   = WAVE CHARNOCK
!              PRF5    = R
!              PVEG5   = VEGETATION
!              QN      = Negative values of humidity set to zero (AJGDB - commonality with hop call to qneglim?)
!              PWL5    = SKIN RESERVOIR
!              PZ05    = ROUGHNESS LENGTH
!              PCPTGZ5 = STATIC ENERGY (LOWER LEVEL)
!              PWET5   = WETNESS
!              PQS5    = SATURATION HUMIDITY
!              PCPTS5  = STATIC ENERGY AT THE SURFACE
!              PRI5    = RICHARDSON NO.
!              PBN5    = STABILITY FUNCTION (NEUTRAL)
!              PBM5    = DRAG COEFFICIENT FOR MOMENTUM
!              PBH5    = DRAG COEFFICIENT FOR HEAT
!              PLNPR5  = LOG. OF THE RATIO OF PRESSURE
!              PALPH5  = COEFFICIENTS OF THE HYDROSTATICS
!              PTCLS5    = CLS TEMPERATURE
!              PRESH5  = HALF LEVEL PRESS. AT OBS. POINTS
!              PRESF5  = FULL LEVEL PRESS. AT OBS. POINTS
!              PGEOPH5 = HALF LEVEL GEOP. AT OBS. POINTS
!              PGEOPF5 = FULL LEVEL GEOP. AT OBS. POINTS
!              PAPHI   = SOME SORT OF MUTANT GEOPOTENTIAL USED BY SOME SURFACE OPERATORS
!              PXP5    = INTERMEDIATE ARRAY (FOR INTERPOLATION IN P)
!              PXPD5   = ...........................................
!              PXZ5    = INTERMEDIATE ARRAY (FOR INTERPOLATION IN Z)
!              PXZD5   = ...........................................
!              PCP5    = Cp
!              PRF5    = R
!              PKAP5   = Kappa whatever that is
!              PUF5    = U VALUES (MODEL) AT OB. POINS
!              PVF5    = V VALUES (MODEL) AT OB. POINS
!              PTF5    = T VALUES (MODEL) AT OB. POINS
!              PQF5    = Q VALUES (MODEL) AT OB. POINS
!              PO3F5   = O3 VALUES (MODEL) AT OB. POINS
!              PCLWF5  = CLW VALUES (MODEL) AT OB. POINS
!              PCIWF5  = CIW VALUES (MODEL) AT OB. POINS
!              PCCF5   = CC VALUES (MODEL) AT OB. POINS
!              PCSWF5    = S (snow) VALUES (MODEL) AT OB. POINS
!              PCRWF5    = R (rain) VALUES (MODEL) AT OB. POINS
!              PCGWF5    = G (graupel) VALUES (MODEL) AT OB. POINS
!              PCLWD5  = Diagnostic CLW from callpar
!              PCIWD5  =      "     CIW      "
!              PCCD5   =      "     CC       "
!              PRFL5   = Diagnostic stratiform rain flux from callpar
!              PSFL5   =      "     stratiform snow flux      "
!              PCRFL5  =      "     convective rain flux      "
!              PCSFL5  =      "     convective snow flux      "
!              PPFRC5  =      "     large-scale precipitation fraction
!              PTSTAR5 = T*
!              PT05    = T0
!              PXYB5   = Output from GPHPRE
!              PAL     = Albedo AT OBS. PONTS
!              VGEOM   = vertical model geometry profiles

! Modifications.
! --------------
!     A.Geer      17-Apr-2015   First version
!     A.Geer      31-Dec-2015   Finally removed the spurious vertical dimension 0:NFLEVG
!     N.Bormann    1-Feb-2016   Slant-path and add lat/lon/hgt
!     A.Geer      06-Apr-2016   Copy, dot-product for test harness/hop_driver
!     A.Geer      01-Jul-2016   Last of 3 supergom/slant-path/vgeom merge bugs
!     P.Lean      17-Nov-2017   Optimisation for slant-path ; change loop ordering
!     P.Lopez     24-May-2017   Moved accumulated precipitation from model to GOM.
!     A.Geer      10-Apr-2018   Rationalise treatment of -ve Q using QNEGLIM which preserve TL/AD sensitivity
!     F.Duruisseau 24-Sep-2018  Add hydrometeors sum (resol + conv) for BAYRAD
!     B.Ingleby   09-Jan-2019   Add 2m humidity and 10m wind to GOM.
!     Y.Hirahara  01-Apr-2019   Add GOM variables for CMEM (Soil, Snow, Lake) 
!     A.Geer      21-May-2019   NetCDF format gom_plus_dump to create colocated model datasets
!     R. El Khatib -10-Apr-2020 Remove gfortran workaround now that this version of the code requires a higher fixed version
!     F. Suzat    07-06-21      Correct back memory leak inserted by error 

use parkind1  , only : jpim, jpib, jprb
use yomhook   , only : lhook, dr_hook
use gom_mod   , only : type_gom, gom_get, gom_set, gid, gom_is_interpolated, &
  &                    gom_missing_value, gom_subvar_number, gom_nflevg, gom_subvar_codes_get, &
  &                    gom_nhoriz, gom_nlev_soil
use yomsats   , only: rdxmax_slant
use yomobs    , only : lcachmt, lvdftraj, lvdfmin
use yomancs   , only : rmdi
use intdyn_mod, only : yytxyb
use yomsta    , only : nlextrap
use yomct0    , only : lecmwf, l_oops, lelam
!use yoephy    , only : yrephy
!use yomgem    , only : tgem
use yomcst    , only : rg, rpi
use yomcoctp  , only : nscatt
!use yomphy    , only : yrphy
use yomvert   , only : tvertical_geom,alloc_copy_vertical_geom,dealloc_vertical_geom
use pardim    , only : jpnppm

implicit none

private

public type_gom_plus, gom_plus_create, gom_plus_create_ad, gom_plus_dim_alloc, gom_plus_bgobs, &
  & gom_plus_destroy, gom_plus_dump, gems_type_ghg, gems_type_aero, gems_type_chem, gom_plus_read_dump, &
  & gom_plus_dump_open, gom_plus_dump_close, type_gom_plus_dump, &
  & gom_plus_dotprod, ih, gom_plus_copy, gom_plus_zero, gom_plus_make_slant_path_tl, &
  & gom_plus_make_slant_path_ad, gom_plus_make_slant_path

type type_gom_plus
  integer(kind=jpim) :: ndlen, nflevg, nhoriz, nppm, nxyb, ngems, nlev_soil
  logical            :: lsrg, lhydro_sum, lphys, lcanari
  real(kind=jprb) :: missing_value
  real(kind=jprb), pointer, contiguous :: lat   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: lon   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: timestep(:,:)  => null()
  real(kind=jprb), pointer, contiguous :: orog  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: al    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: tstar (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: t0    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: ts    (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: ls    (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: ci    (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: ws    (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: sn    (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: u10n  (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: v10n  (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: u0    (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: v0    (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: chn   (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: tcls  (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: hucls (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: ucls  (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: vcls  (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: nucls (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: nvcls (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: z0    (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: wl    (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: sst   (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: veg   (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: cptgz (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: wet   (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: qs    (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: cpts  (:,:)    => null()  
  real(kind=jprb), pointer, contiguous :: ri    (:,:)    => null()  
  real(kind=jprb), pointer, contiguous :: bn    (:,:)    => null()  
  real(kind=jprb), pointer, contiguous :: bm    (:,:)    => null()  
  real(kind=jprb), pointer, contiguous :: bh    (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: upd   (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: arg   (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: hv    (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: z0h   (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: wsi   (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: rrr   (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: rsn   (:,:)    => null() 
  real(kind=jprb), pointer, contiguous :: es    (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: prec_accum (:,:) => null()
  real(kind=jprb), pointer, contiguous :: sdfor (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: t2m   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: d2m   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: q2m   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: u10m  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: v10m  (:,:)    => null()
  !LLL
  real(kind=jprb), pointer, contiguous :: soty  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: lail  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: laih  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: sndt  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: sndp  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: cvl   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: cvh   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: tvl   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: tvh   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: tlk   (:,:)    => null()
!canari
  real(kind=jprb), pointer, contiguous :: cori  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: eh    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: ez    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: est   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: esn   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: et2   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: eh2   (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: ev1   (:,:)    => null()
!canari
  real(kind=jprb), pointer, contiguous :: presh (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: presf (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: geoph (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: geopf (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: aphi  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: cp    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: rf    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: kap   (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: uf    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: vf    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: tf    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: qf    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: qn    (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: o3f   (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: clwf  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: ciwf  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: ccf   (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: cswf  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: crwf  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: cgwf  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: clwd  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: ciwd  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: clw_sum  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: ciw_sum  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: crw_sum  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: csw_sum  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: ccd   (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: phys_type(:,:,:)  => null() !KaLo
  real(kind=jprb), pointer, contiguous :: rfl   (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: sfl   (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: crfl  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: csfl  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: pfrc  (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: hgt   (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: tsoil (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: qsoil (:,:,:)  => null()
  real(kind=jprb), pointer, contiguous :: gems  (:,:,:,:) => null()
  integer(kind=jpim), pointer, contiguous :: gems_igrib (:)  => null()
  integer(kind=jpim), pointer, contiguous :: gems_type (:)   => null()
  real(kind=jprb), pointer, contiguous :: pxp   (:,:,:,:) => null()
  real(kind=jprb), pointer, contiguous :: pxpd  (:,:,:,:) => null()
  real(kind=jprb), pointer, contiguous :: pxz   (:,:,:,:) => null()
  real(kind=jprb), pointer, contiguous :: pxzd  (:,:,:,:) => null()
  real(kind=jprb), pointer, contiguous :: pxyb  (:,:,:,:) => null()
! Zonal + meridonal component of geographical north direction (angle):
  real(kind=jprb), pointer, contiguous :: gnordl  (:,:)    => null()
  real(kind=jprb), pointer, contiguous :: gnordm  (:,:)    => null()
  real(kind=jprb), private, pointer, contiguous :: store_sfc(:,:,:) => null()
  real(kind=jprb), private, pointer, contiguous :: store_ua0(:,:,:,:) => null()
  real(kind=jprb), private, pointer, contiguous :: store_ua1(:,:,:,:) => null()
  real(kind=jprb), private, pointer, contiguous :: store_soil(:,:,:,:) => null()
  real(kind=jprb), private, pointer, contiguous :: store_ppm(:,:,:,:,:) => null()
  type(tvertical_geom)    , pointer :: vgeom => null()
end type

integer(kind=jpim), parameter :: ndump_filename_len=50

type type_gom_plus_dump
  integer(kind=jpim) :: file_id
  character(len=ndump_filename_len) :: clfilename
  logical            :: lnetcdf
  integer(kind=jpim) :: nsets
  integer(kind=jpim) :: counter
  integer(kind=jpim), pointer :: ipos_set(:) => null()
  integer(kind=jpim), pointer :: dimid(:) => null()
end type

save

! In lieu of anything more complicated (a future enhancement...) this is 
! how we reduce the 2D-permitting arrays to 1D for the majority of obsops:
integer(kind=jpim), parameter :: ih=1

integer(kind=jpim), parameter :: gems_type_ghg = 1
integer(kind=jpim), parameter :: gems_type_aero = 2
integer(kind=jpim), parameter :: gems_type_chem = 3

#include "abor1.intfb.h"
#include "ctstar.intfb.h"
#include "gpgeo.intfb.h"
#include "gprcp_qlirsg.intfb.h"
#include "gphpre.intfb.h"
#include "ppinit.intfb.h"
#include "ppinitz.intfb.h"
#include "ctstartl.intfb.h"
#include "calc_geom_height.intfb.h"
#include "calc_geom_height_tl.intfb.h"
#include "calc_geom_height_ad.intfb.h"
#include "gpgeotl.intfb.h"
#include "gprcptl.intfb.h"
#include "gphpretl.intfb.h"
#include "ppinittl.intfb.h"
#include "ppinitztl.intfb.h"
#include "surf_inq.h"
#include "exchco.intfb.h"
#include "qneglim.intfb.h"
#include "qneglimtl.intfb.h"
#include "rousea.intfb.h"
#include "z0sea.intfb.h"
#include "exchco_vdf.intfb.h"
#include "surbound.intfb.h"
#include "surboundtl.intfb.h"
#include "z0seatl.intfb.h"
#include "exchco_vdftl.intfb.h"
#include "exchcotl.intfb.h"
#include "rouseatl.intfb.h"
#include "gprcpad.intfb.h"
#include "gphpread.intfb.h"
#include "z0seaad.intfb.h"
#include "exchco_vdfad.intfb.h"
#include "ppinitad.intfb.h"
#include "exchcoad.intfb.h"
#include "gpgeoad.intfb.h"
#include "surboundad.intfb.h"
#include "ppinitza.intfb.h"
#include "ctstarad.intfb.h"
#include "rouseaad.intfb.h"
#include "netcdf.inc"

contains

! -------------------------------------------------------------
! Create and fill the structure - direct or TL version
!
! kdlen  - length of obs dimension in arrays
! kposm  - array of GOM IDs of these observations
! kobtype- observation type, used to determine if a GOM variable is available
!
! -------------------------------------------------------------
subroutine gom_plus_create(ydephy,ydml_phy_mf,gp,g,kdlen,kposm,kobtype,ydvert_geom,gp5)

use model_physics_mf_mod , only : model_physics_mf_type
use yoephy               , only : tephy
type(tephy)                ,intent(inout) :: ydephy
type(model_physics_mf_type),intent(inout):: ydml_phy_mf
type(type_gom_plus), intent(inout) :: gp
type(type_gom),      intent(in)    :: g
integer(kind=jpim),  intent(in)    :: kdlen
integer(kind=jpim),  intent(in)    :: kposm(kdlen)
integer(kind=jpim),  intent(in)    :: kobtype
type(tvertical_geom),intent(in)    :: ydvert_geom 
type(type_gom_plus), optional, intent(in) :: gp5

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_create',0,zhook_handle)

call gom_plus_dim_alloc(gp,g,kdlen,kposm,kobtype)
call alloc_copy_vertical_geom( gp%vgeom, ydvert_geom%yrvab,  ydvert_geom%yrveta,ydvert_geom%yrvfe) 
call gom_plus_zero(gp) ! AJGDB why is this necessary? Came in with 43r1_v3 it seems. Niels? 
call gom_plus_fill(ydephy,ydml_phy_mf%yrphy,gp,g,kposm,kobtype,gp5)
if (present(gp5)) then
  call gom_plus_derived_tl(ydephy,ydml_phy_mf%yrphy0,ydml_phy_mf%yrphy1,gp,gp5,kobtype)
else
  call gom_plus_derived(ydephy,ydml_phy_mf%yrphy0,ydml_phy_mf%yrphy1,gp,kobtype)
endif

if (lhook) call dr_hook('gom_plus_create',1,zhook_handle)
end subroutine gom_plus_create

! -------------------------------------------------------------
! Create and fill the structure
!
! kdlen  - length of obs dimension in arrays
! kposm  - array of GOM IDs of these observations
! kobtype- observation type, used to determine if a GOM variable is available
!
! Adjoint is done in 2 calls: 1) with ld_alloc_only=.true. 2) with ld_alloc_only=.false.
!
! -------------------------------------------------------------
subroutine gom_plus_create_ad(ydephy,ydphy0,ydphy1,gp,g,kdlen,kposm,kobtype,&
  & gp5,ydvert_geom,ld_alloc_only)

use yoephy  , only : tephy
use yomphy0 , only : tphy0
use yomphy1 , only : tphy1
type(tephy)         ,intent(inout) :: ydephy
type(tphy0)         ,intent(inout) :: ydphy0
type(tphy1)         ,intent(inout) :: ydphy1
type(type_gom_plus), intent(inout) :: gp
type(type_gom),      intent(inout) :: g
integer(kind=jpim),  intent(in)    :: kdlen
integer(kind=jpim),  intent(in)    :: kposm(kdlen)
integer(kind=jpim),  intent(in)    :: kobtype
type(type_gom_plus), intent(in)    :: gp5
type(tvertical_geom),intent(in)    :: ydvert_geom
logical,             intent(in)    :: ld_alloc_only

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_create_ad',0,zhook_handle)

if (ld_alloc_only) then
  ! First call - set up and zero arrays to receive adjoints
  call gom_plus_dim_alloc(gp,g,kdlen,kposm,kobtype)
  call alloc_copy_vertical_geom( gp%vgeom, ydvert_geom%yrvab,  ydvert_geom%yrveta,ydvert_geom%yrvfe) 
  call gom_plus_zero(gp)
else
  ! Second call - adjoint of gom plus; adjoints returned to goms.
  call gom_plus_derived_ad(ydephy,ydphy0,ydphy1,gp,gp5,kobtype)
  call gom_plus_fill_ad(gp,g,kposm,kobtype,gp5)
 ! call gom_plus_zero(gp) ! Nice but expensive /MH
endif

if (lhook) call dr_hook('gom_plus_create_ad',1,zhook_handle)
end subroutine gom_plus_create_ad

! -------------------------------------------------------------
! (Create and fill the gom plus structure) - special behaviour for bgobs
! which fills the gom plus arrays itself but still needs
! to generate the "derived" stuff
! -------------------------------------------------------------
subroutine gom_plus_bgobs(ydephy,ydphy0,ydphy1,gp,kobtype,ydvert_geom,gp5)

use yoephy  , only : tephy
use yomphy0 , only : tphy0
use yomphy1 , only : tphy1
type(tephy)         ,intent(inout) :: ydephy
type(tphy0)         ,intent(inout) :: ydphy0
type(tphy1)         ,intent(inout) :: ydphy1
type(type_gom_plus) ,intent(inout) :: gp
integer(kind=jpim)  ,intent(in)    :: kobtype
type(tvertical_geom),intent(in)    :: ydvert_geom 
type(type_gom_plus) ,optional, intent(in) :: gp5

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_bgobs',0,zhook_handle)

if (present(gp5)) then
  call gom_plus_derived_tl(ydephy,ydphy0,ydphy1,gp,gp5,kobtype)
else
  call gom_plus_derived(ydephy,ydphy0,ydphy1,gp,kobtype)
endif

if (lhook) call dr_hook('gom_plus_bgobs',1,zhook_handle)
end subroutine gom_plus_bgobs

! -------------------------------------------------------------
! Dimension the gom plus and call lower level routine to allocate.
! -------------------------------------------------------------
subroutine gom_plus_dim_alloc(gp,g,kdlen,kposm,kobtype)

type(type_gom_plus) ,intent(inout) :: gp
type(type_gom)      ,intent(in)    :: g
integer(kind=jpim),  intent(in)    :: kdlen
integer(kind=jpim),  intent(in)    :: kposm(kdlen)
integer(kind=jpim),  intent(in)    :: kobtype

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_create',0,zhook_handle)

! Determine dimensions from the GOM itself, not from global modules that should be private
gp%nppm      = jpnppm ! (!!!should be private) used by ppinint and pp routines
gp%ndlen     = kdlen
gp%nflevg    = gom_nflevg(g)     ! number of vertical levels
gp%nlev_soil = gom_nlev_soil(g)  ! number of vertical levels (soil)
if(kdlen > 0) then
  if(kposm(1)/=0) then
    ! number of horiz points in 2D goms, obtained from the GOMs
    ! Assume that nhoriz is constant within one set of obs (enforced by ecset.F90)
    gp%nhoriz = gom_nhoriz(g,kposm(1))
  else
    ! AJGDB hack to allow bgobs to call us with a fake kposm (all set to zero)
    gp%nhoriz=1
  endif
else
  gp%nhoriz=0
endif
gp%nxyb=yytxyb%ndim !! again a sign of non-modular code

gp%missing_value = gom_missing_value

! Currently, the following variables are only allocated when they are present in the GOMS. It is intended
! to make this fully automatic so that ALL variables coming from GOMS are only allocated if the
! GOM itself is available.
gp%ngems = 0
if(gom_is_interpolated(g, gid%ghg,  kobtype))       gp%ngems=gp%ngems+gom_subvar_number(g, gid%ghg)
if(gom_is_interpolated(g, gid%aero, kobtype))       gp%ngems=gp%ngems+gom_subvar_number(g, gid%aero)
if(gom_is_interpolated(g, gid%chem_assim, kobtype)) gp%ngems=gp%ngems+gom_subvar_number(g, gid%chem_assim)

gp%lsrg = ( gom_is_interpolated(g, gid%s, kobtype) &
     & .or. gom_is_interpolated(g, gid%r, kobtype) &
     & .or. gom_is_interpolated(g, gid%g, kobtype) )

gp%lhydro_sum = ( gom_is_interpolated(g, gid%l_sum, kobtype) &
           & .or. gom_is_interpolated(g, gid%i_sum, kobtype) &
           & .or. gom_is_interpolated(g, gid%r_sum, kobtype) &
           & .or. gom_is_interpolated(g, gid%s_sum, kobtype))

gp%lphys = gom_is_interpolated(g, gid%phys, kobtype)

gp%lcanari = ( gom_is_interpolated(g, gid%eh, kobtype) .and. & 
           &   gom_is_interpolated(g, gid%ez, kobtype) )

call gom_plus_alloc(gp, kdlen, gp%nflevg, gp%nhoriz, gp%nppm, gp%nxyb, gp%lsrg, gp%lhydro_sum, & 
           & gp%lphys, gp%lcanari, gp%ngems, gp%nlev_soil)

if (lhook) call dr_hook('gom_plus_dim_alloc',1,zhook_handle)
end subroutine gom_plus_dim_alloc

! --------------------------------------------------------------
! Fill gom_plus from gom. This is intended to be called
! in both direct and TL cases, and to be ignorant of whether
! it is called for TL or direct usage.
! --------------------------------------------------------------
subroutine gom_plus_fill(ydephy,ydphy,gp,g,kposm,kobtype,gp5)

use yoephy , only : tephy
use yomphy , only : tphy
type(tephy)          ,intent(inout) :: ydephy
type(tphy)           ,intent(inout) :: ydphy
type(type_gom_plus)  ,intent(inout) :: gp
type(type_gom)       ,intent(in)    :: g 
integer(kind=jpim)   ,intent(in)    :: kposm(gp%ndlen)
integer(kind=jpim)   ,intent(in)    :: kobtype
type(type_gom_plus),optional,intent(in)    :: gp5

real(kind=jprb), dimension(gp%ndlen,gp%nhoriz) :: zcvh, ztvh, zcvl, ztvl
real(kind=jprb) :: zrtfreezsice
real(kind=jprb),allocatable :: zrvcov(:)
integer(kind=jpim) :: jobs, jhoriz, igems, igems_start, isubvar,ivtypes
logical :: llgot, llgot1, llgot2, llgot3, ll_tl, llgot_u

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_fill',0,zhook_handle)

zrtfreezsice = 271.46_jprb
!CP to check if "if..." is really required. Seems ok side MF in 43t2
!if ((ydphy%lmphys .or. ydephy%lephys).and.&
! & (ydphy%lrayfm15 .or. .not. ydphy%lray))&
call surf_inq(ydephy%ysurf,knvtypes=ivtypes)
allocate(zrvcov(0:ivtypes))
call surf_inq(ydephy%ysurf,prtfreezsice=zrtfreezsice,prvcov=zrvcov)

do jobs = 1,gp%ndlen

  ! Latitude/longitude
  llgot = gom_get(g, kposm(jobs), gid%lat, gp%lat(jobs,:))
  if ((.not.llgot)) then
    call abor1('gom_plus: latitude undefined')
  endif
  llgot = gom_get(g, kposm(jobs), gid%lon, gp%lon(jobs,:))
  if ((.not.llgot)) then
    call abor1('gom_plus: longitude undefined')
  endif

  ! Time(step)
  llgot = gom_get(g, kposm(jobs), gid%timestep, gp%timestep(jobs,:))

  ! orography
  llgot = gom_get(g, kposm(jobs), gid%orog, gp%orog(jobs,:))

  ! compas
  if (lelam) then
    llgot = gom_get(g, kposm(jobs), gid%gnordl, gp%gnordl(jobs,:))
    llgot = gom_get(g, kposm(jobs), gid%gnordm, gp%gnordm(jobs,:))
  endif

  ! albedo 
  llgot = gom_get(g, kposm(jobs), gid%vf_albf, gp%al(jobs,:))

  ! AJGDB U10N, V10N, WS, TS, Z0, SN in preint/s used to have "set_missing=0.0" only in the TL
  ! AJGDB ultimately no-one should rely on this as in the future, where missing, the relevant gom_plus pointer 
  ! will be NULL
  llgot = gom_get(g, kposm(jobs), gid%rr_t, gp%ts(jobs,:), pset_missing=0.0_jprb)

  llgot = gom_get(g, kposm(jobs), gid%sb_q, gp%ws(jobs,:), pset_missing=0.0_jprb)
  if(.not.llgot) then
    llgot = gom_get(g, kposm(jobs), gid%rr_w, gp%ws(jobs,:), pset_missing=0.0_jprb)
  endif

  ! soil (layers)
  llgot = gom_get(g, kposm(jobs), gid%sb_tly, gp%tsoil(jobs,1:gp%nlev_soil,:), pset_missing=0.0_jprb)
  llgot = gom_get(g, kposm(jobs), gid%sb_qly, gp%qsoil(jobs,1:gp%nlev_soil,:), pset_missing=0.0_jprb)

  ! lake 
  llgot = gom_get(g, kposm(jobs), gid%sl_lmlt, gp%tlk(jobs,:), pset_missing=0.0_jprb)

  llgot = gom_get(g, kposm(jobs), gid%sg_f, gp%sn(jobs,:),   pset_missing=0.0_jprb)
  llgot = gom_get(g, kposm(jobs), gid%sg_r, gp%sndt(jobs,:), pset_missing=0.0_jprb)

  llgot = gom_get(g, kposm(jobs), gid%vd_10nu, gp%u10n(jobs,:), pset_missing=0.0_jprb)
  llgot = gom_get(g, kposm(jobs), gid%vd_10nv, gp%v10n(jobs,:), pset_missing=0.0_jprb)
  llgot = gom_get(g, kposm(jobs), gid%vd_upd,  gp%upd(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%vd_2t,   gp%t2m(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%vd_2d,   gp%d2m(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%vd_2q,   gp%q2m(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%vd_10u,  gp%u10m(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%vd_10v,  gp%v10m(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%vf_z0f,  gp%z0(jobs,:), pset_missing=0.0_jprb)
  if (lecmwf) then
    ! AJGDB please, can't we agree on units????
    gp%z0(jobs,:)=gp%z0(jobs,:)/rg
  endif

  ! Soil type
  llgot = gom_get(g, kposm(jobs), gid%vf_soty,  gp%soty(jobs,:))

  ! LAI
  llgot = gom_get(g, kposm(jobs), gid%vf_lail,  gp%lail(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%vf_laih,  gp%laih(jobs,:))

  ! lsm and sea ice
  !CP shifted after ts initialization as it is required for gp%ci (MF recipe)
  !CP I agree, MF should implement the EC's way (but the seaice field is still missing) 
  llgot = gom_get(g, kposm(jobs), gid%vf_lsm,  gp%ls(jobs,:))
  if (lecmwf) then
    llgot = gom_get(g, kposm(jobs), gid%vf_ci, gp%ci(jobs,:))
  else
    gp%ci(jobs,:)=0.0_jprb
    !CP Old .not. ld_surf case preferred: 
    where(gp%ts(jobs,:) < zrtfreezsice) gp%ci(jobs,:) = 1._jprb - gp%ls(jobs,:)
    !where(gp%ts(jobs,:) < zrtfreezsice) gp%ci(jobs,:) = 1._jprb
  endif

  gp%wl(jobs,:) = 1.0_jprb ! gom modernisation: this used to be done in the old cobs.f90 but seemed rather pointless

  llgot = gom_get(g, kposm(jobs), gid%vf_ucur,  gp%u0(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%vf_vcur,  gp%v0(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%vf_sdfor, gp%sdfor(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%ws_char,  gp%chn(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%cl_tcls,  gp%tcls(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%cl_hucls, gp%hucls(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%cl_ucls, gp%ucls(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%cl_vcls, gp%vcls(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%cl_nucls, gp%nucls(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%cl_nvcls, gp%nvcls(jobs,:))

  llgot = gom_get(g, kposm(jobs),   gid%rr_t_ocean, gp%sst(jobs,:)) 
  if (.not.llgot) then
    ! surface temperature needs to be stored in psst5 whether rr_t_ocean
    ! (surf. temp. interpolated with lsm) is used or not
    gp%sst(jobs,:)=gp%ts(jobs,:)
  endif

  if (.not. l_oops .and. ydephy%lephys) then
    if (lecmwf) then
      ! Fudge veg, because it is only used for surbound, which does not reflect
      ! the ecmwf surface parameterization, anyhow. if it is needed for some
      ! more serious purpose, the arrays tvl,tvh,cvl,cvh should be ouput 
      ! AJGDB basically ECMWF should just send VEG through the GOMS rather
      ! than trying to fudge it here. At least fudge it in the model, not here!
      llgot  = gom_get(g, kposm(jobs), gid%vf_cvh, zcvh(jobs,:)) 
      llgot1 = gom_get(g, kposm(jobs), gid%vf_tvh, ztvh(jobs,:))
      llgot2 = gom_get(g, kposm(jobs), gid%vf_cvl, zcvl(jobs,:))
      llgot3 = gom_get(g, kposm(jobs), gid%vf_tvl, ztvl(jobs,:))
       if(llgot.and.llgot1.and.llgot2.and.llgot3) then
         gp%veg(jobs,:)=&
           & zcvh(jobs,:)*zrvcov(nint(ztvh(jobs,:)))+&
           & zcvl(jobs,:)*zrvcov(nint(ztvl(jobs,:)))  
       else
         gp%veg(jobs,:)=0.0_jprb
       endif
    else
      llgot = gom_get(g, kposm(jobs), gid%vf_veg, gp%veg(jobs,:))
    endif
  else
    ! Fudge Veg for OOPS and IFS without Physics 
    gp%veg(jobs,:)=0.0_jprb
  endif

  llgot = gom_get(g, kposm(jobs), gid%vf_cvl, gp%cvl(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%vf_cvh, gp%cvh(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%vf_tvl, gp%tvl(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%vf_tvh, gp%tvh(jobs,:))

  ! Canari land-surface variables from PPOBSAC (AJGDB) 
  llgot = gom_get(g, kposm(jobs), gid%vv_arg, gp%arg(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%vv_hv,  gp%hv (jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%vv_z0h, gp%z0h(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%vv_sab, gp%wsi(jobs,:))

  ! Canari snow stuff (?) from PPOBSN (AJGDB) Note that "SN" is already available
  llgot = gom_get(g, kposm(jobs), gid%x2_prwa, gp%rrr(jobs,:))
  llgot = gom_get(g, kposm(jobs), gid%x2_prsn, gp%rsn(jobs,:))

  ! surface pressure - presh processing is completed in "gom_plus_derived"
  llgot = gom_get(g, kposm(jobs), gid%sp, gp%presh(jobs,gp%nflevg,:))
  if (.not.llgot) then
    call abor1('gom_plus: surface pressure undefined')
  endif

  llgot = gom_get(g, kposm(jobs), gid%u, gp%uf(jobs,:,:))
  llgot_u = llgot
  if(.not. llgot) then
    if(gom_get(g, kposm(jobs), gid%ul, gp%uf(jobs,gp%nflevg,:))) then
      ! This situation occurs for observations where interpolating
      ! the whole vertical profile would be a waste.
      do jhoriz=1,gp%nhoriz
        gp%uf(jobs,1:gp%nflevg-1,jhoriz) = gp%uf(jobs,gp%nflevg,jhoriz)
      enddo
    else
      call abor1('gom_plus: gp%uf undefined')
    endif
  endif

  llgot = gom_get(g, kposm(jobs), gid%v, gp%vf(jobs,1:gp%nflevg,:))
  if(.not. llgot) then
    if(gom_get(g, kposm(jobs), gid%vl, gp%vf(jobs,gp%nflevg,:))) then
      do jhoriz=1,gp%nhoriz
        gp%vf(jobs,1:gp%nflevg-1,jhoriz) = gp%vf(jobs,gp%nflevg,jhoriz)
      enddo
    else
      call abor1('gom_plus: gp%vf undefined')
    endif
  endif

  llgot = gom_get(g, kposm(jobs), gid%t, gp%tf(jobs,1:gp%nflevg,:))
  if(.not. llgot) then
    if(gom_get(g, kposm(jobs), gid%tl, gp%tf(jobs,gp%nflevg,:))) then
      do jhoriz=1,gp%nhoriz
        gp%tf(jobs,1:gp%nflevg-1,jhoriz) = gp%tf(jobs,gp%nflevg,jhoriz)
      enddo
    else
      call abor1('gom_plus: gp%tf undefined')
    endif
  endif

  llgot = gom_get(g, kposm(jobs), gid%q, gp%qf(jobs,1:gp%nflevg,:))
  if(.not. llgot) then
    if(gom_get(g, kposm(jobs), gid%ql, gp%qf(jobs,gp%nflevg,:))) then
      do jhoriz=1,gp%nhoriz
        gp%qf(jobs,1:gp%nflevg-1,jhoriz) = gp%qf(jobs,gp%nflevg,jhoriz)
      enddo
    else
      call abor1('gom_plus: gp%qf undefined')
    endif
  endif

  ! initialisation of ozone for rttov
  ! AJGDB gp%o3f(jobs,1:gp%nflevg,:) = 1.e-07 ! AJGDB this probably doesn't work as it will be overwritten...
  ! AJGDB let's see if this works - we need to tell later on (e.g. hradp_ml) if we have a valid O3 field or not:
  ! llgot = gom_get(g, kposm(jobs), gid%o3, gp%o3f(jobs,1:nflevg,:), pset_missing=0.0_jprb)
  llgot = gom_get(g, kposm(jobs), gid%o3, gp%o3f(jobs,1:gp%nflevg,:))

  ! AJGDB L,I,A were always set to 0 when missing in direct and TL case - this could be misleading and we 
  ! could try to remove it in the future. 
  if (.not.LECMWF) then
    llgot = gom_get(g, kposm(jobs), gid%l, gp%clwf(jobs,1:gp%nflevg,:))!CPtoclean, pset_missing=0.0_jprb) !CP/AJGBD GOM_MISSING_VALUE used instead of
    llgot = gom_get(g, kposm(jobs), gid%i, gp%ciwf(jobs,1:gp%nflevg,:))!CPtoclean, pset_missing=0.0_jprb)
    llgot = gom_get(g, kposm(jobs), gid%a, gp%ccf(jobs,1:gp%nflevg,:))!CPtoclean, pset_missing=0.0_jprb)
  else
    ! AJGDB - the slant path code (used only at ECMWF) is incompatible with GOM_MISSING_VALUE
    llgot = gom_get(g, kposm(jobs), gid%l, gp%clwf(jobs,1:gp%nflevg,:), pset_missing=0.0_jprb)
    llgot = gom_get(g, kposm(jobs), gid%i, gp%ciwf(jobs,1:gp%nflevg,:), pset_missing=0.0_jprb)
    llgot = gom_get(g, kposm(jobs), gid%a, gp%ccf(jobs,1:gp%nflevg,:), pset_missing=0.0_jprb)
  endif

  if (gom_is_interpolated(g, gid%s, kobtype))&
    & llgot = gom_get(g, kposm(jobs), gid%s, gp%cswf(jobs,1:gp%nflevg,:), pset_missing=0.0_jprb)

  if (gom_is_interpolated(g, gid%r, kobtype))&
    & llgot = gom_get(g, kposm(jobs), gid%r, gp%crwf(jobs,1:gp%nflevg,:), pset_missing=0.0_jprb)

  if (gom_is_interpolated(g, gid%g, kobtype))&
    & llgot = gom_get(g, kposm(jobs), gid%g, gp%cgwf(jobs,1:gp%nflevg,:), pset_missing=0.0_jprb)

  ! conv and resol for bayrad

  if (gom_is_interpolated(g, gid%l_sum, kobtype)) &
    & llgot = gom_get(g, kposm(jobs), gid%l_sum, gp%clw_sum(jobs,1:gp%nflevg,:), pset_missing=0.0_jprb)
  
  if (gom_is_interpolated(g, gid%i_sum, kobtype)) &
    & llgot = gom_get(g, kposm(jobs), gid%i_sum, gp%ciw_sum(jobs,1:gp%nflevg,:), pset_missing=0.0_jprb)

  if (gom_is_interpolated(g, gid%r_sum, kobtype)) &
    & llgot = gom_get(g, kposm(jobs), gid%r_sum, gp%crw_sum(jobs,1:gp%nflevg,:), pset_missing=0.0_jprb)

  if (gom_is_interpolated(g, gid%s_sum, kobtype)) &
    & llgot = gom_get(g, kposm(jobs), gid%s_sum, gp%csw_sum(jobs,1:gp%nflevg,:), pset_missing=0.0_jprb)


  ! AJGDB we should consider unifying some of these all-sky microwave variables with the above
  if (gom_is_interpolated(g, gid%phys, kobtype)) then
    llgot = gom_get(g, kposm(jobs), gid%phys, 1, gp%clwd(jobs,1:gp%nflevg,:))
    llgot = gom_get(g, kposm(jobs), gid%phys, 2, gp%ciwd(jobs,1:gp%nflevg,:))
    llgot = gom_get(g, kposm(jobs), gid%phys, 3, gp%ccd (jobs,1:gp%nflevg,:))
    llgot = gom_get(g, kposm(jobs), gid%phys, 4, gp%rfl (jobs,1:gp%nflevg,:))
    llgot = gom_get(g, kposm(jobs), gid%phys, 5, gp%sfl (jobs,1:gp%nflevg,:))
    llgot = gom_get(g, kposm(jobs), gid%phys, 6, gp%crfl(jobs,1:gp%nflevg,:))
    llgot = gom_get(g, kposm(jobs), gid%phys, 7, gp%csfl(jobs,1:gp%nflevg,:))
    llgot = gom_get(g, kposm(jobs), gid%phys, 8, gp%pfrc(jobs,1:gp%nflevg,:))
    llgot = gom_get(g, kposm(jobs), gid%phys, 9, gp%phys_type(jobs,1:gp%nflevg,:)) !KaLo
  endif

  if (gom_is_interpolated(g, gid%prec_accum, kobtype)) then
    llgot = gom_get(g, kposm(jobs), gid%prec_accum, gp%prec_accum(jobs,:))
  endif

  igems = 0
  if(gom_is_interpolated(g, gid%ghg, kobtype)) then
    igems_start = igems+1
    do isubvar = 1,gom_subvar_number(g, gid%ghg)
      igems = igems+1
      llgot = gom_get(g, kposm(jobs), gid%ghg, isubvar, gp%gems(jobs,1:gp%nflevg,:,igems))  
      gp%gems_type(igems) = gems_type_ghg   
    enddo
    call gom_subvar_codes_get(g, gid%ghg, gp%gems_igrib(igems_start:igems)) 
  endif

  if(gom_is_interpolated(g, gid%aero, kobtype)) then
    igems_start = igems+1
    do isubvar = 1,gom_subvar_number(g, gid%aero)
      igems = igems+1
      llgot = gom_get(g, kposm(jobs), gid%aero, isubvar, gp%gems(jobs,1:gp%nflevg,:,igems))  
      gp%gems_type(igems) = gems_type_aero   
    enddo
    call gom_subvar_codes_get(g, gid%aero, gp%gems_igrib(igems_start:igems)) 
  endif

  if(gom_is_interpolated(g, gid%chem_assim, kobtype)) then
    igems_start = igems+1
    do isubvar = 1,gom_subvar_number(g, gid%chem_assim)
      igems = igems+1
      llgot = gom_get(g, kposm(jobs), gid%chem_assim, isubvar, gp%gems(jobs,1:gp%nflevg,:,igems))  
      gp%gems_type(igems) = gems_type_chem  
    enddo
    call gom_subvar_codes_get(g, gid%chem_assim, gp%gems_igrib(igems_start:igems)) 
  endif

  if (gp%lcanari) then
    llgot = gom_get(g, kposm(jobs), gid%cori, gp%cori(jobs,:))
    llgot = gom_get(g, kposm(jobs), gid%eh, gp%eh(jobs,1:gp%nflevg,:))
    llgot = gom_get(g, kposm(jobs), gid%ez, gp%ez(jobs,1:gp%nflevg,:))
    llgot = gom_get(g, kposm(jobs), gid%est, gp%est(jobs,:))
    llgot = gom_get(g, kposm(jobs), gid%esn, gp%esn(jobs,:))
    llgot = gom_get(g, kposm(jobs), gid%et2, gp%et2(jobs,:))
    llgot = gom_get(g, kposm(jobs), gid%eh2, gp%eh2(jobs,:))
    llgot = gom_get(g, kposm(jobs), gid%ev1, gp%ev1(jobs,:))
  endif

  ! Rotation of winds on grids that don't align with N and E axes
  ll_tl = present(gp5)
  if (lelam.and.ll_tl) then
    call rotate_uv_tlad(jobs, gp, gp5, llgot_u, ll_tl)
  endif
enddo

if(allocated(zrvcov))deallocate(zrvcov)


if (lhook) call dr_hook('gom_plus_fill',1,zhook_handle)
end subroutine gom_plus_fill

! --------------------------------------------------------------
! Fill gom_plus from gom 
! (in adjoint sense it is of course, fill gom from gom_plus)
!
! Not all variables available in the direct and TL have yet
! been added here; they will need to be if used in the future.
! --------------------------------------------------------------
subroutine gom_plus_fill_ad(gp,g,kposm,kobtype,gp5)

type(type_gom_plus),  intent(inout) :: gp
type(type_gom),       intent(inout) :: g 
integer(kind=jpim),   intent(in)    :: kposm(gp%ndlen)
integer(kind=jpim),   intent(in)    :: kobtype
type(type_gom_plus),  intent(in)    :: gp5

integer(kind=jpim) :: jobs
integer(kind=jpim) :: igems, igems_start, isubvar
logical :: llset, ll_tl, ll_u
real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_fill_ad',0,zhook_handle)

do  jobs = 1,gp%ndlen

  ! Rotation of winds on grids that don't align with N and E axes
  if ( lelam ) then
    ll_tl=.false.
    ll_u = gom_is_interpolated(g,gid%u,kobtype)
    if ((.not.ll_u).and.(.not.gom_is_interpolated(g,gid%ul,kobtype))) &
     & call abor1("gom_plus_fill_ad: this should never happen.")
    call rotate_uv_tlad(jobs, gp, gp5, ll_u, ll_tl)
  endif

  ! AJGDB this is all horrible: 
  if(lecmwf) gp%z0(jobs,:) = gp%z0(jobs,:)/rg

  llset = gom_set(g, kposm(jobs), gid%rr_t,    gp%ts(jobs,:),   ldsum=.true.)
  llset = gom_set(g, kposm(jobs), gid%vf_z0f,  gp%z0(jobs,:),   ldsum=.true.)
  llset = gom_set(g, kposm(jobs), gid%vd_10nu, gp%u10n(jobs,:), ldsum=.true.)
  llset = gom_set(g, kposm(jobs), gid%vd_10nv, gp%v10n(jobs,:), ldsum=.true.)
  llset = gom_set(g, kposm(jobs), gid%vd_2t,   gp%t2m(jobs,:),  ldsum=.true.)
  llset = gom_set(g, kposm(jobs), gid%vd_2d,   gp%d2m(jobs,:),  ldsum=.true.)
  llset = gom_set(g, kposm(jobs), gid%vd_2q,   gp%q2m(jobs,:),  ldsum=.true.)
  llset = gom_set(g, kposm(jobs), gid%vd_10u,  gp%u10m(jobs,:), ldsum=.true.)
  llset = gom_set(g, kposm(jobs), gid%vd_10v,  gp%v10m(jobs,:), ldsum=.true.)
  llset = gom_set(g, kposm(jobs), gid%sb_q,    gp%ws(jobs,:),   ldsum=.true.)
  if(.not.llset) then
    llset = gom_set(g, kposm(jobs), gid%rr_w,  gp%ws(jobs,:),   ldsum=.true.)
  endif
  llset = gom_set(g, kposm(jobs), gid%sg_f,    gp%sn(jobs,:),   ldsum=.true.)

  llset = gom_set(g, kposm(jobs), gid%u, gp%uf(jobs,1:gp%nflevg,:), ldsum=.true.)
  if(.not.llset) then
    llset = gom_set(g, kposm(jobs), gid%ul, sum(gp%uf(jobs,1:gp%nflevg,:),dim=1), ldsum=.true.)
  endif

  llset = gom_set(g, kposm(jobs), gid%v, gp%vf(jobs,1:gp%nflevg,:), ldsum=.true.)
  if(.not.llset) then
    llset = gom_set(g, kposm(jobs), gid%vl, sum(gp%vf(jobs,1:gp%nflevg,:),dim=1), ldsum=.true.)
  endif

  llset = gom_set(g, kposm(jobs), gid%t, gp%tf(jobs,1:gp%nflevg,:), ldsum=.true.)
  if(.not.llset) then
    llset = gom_set(g, kposm(jobs), gid%tl, sum(gp%tf(jobs,1:gp%nflevg,:),dim=1), ldsum=.true.)
  endif

  llset = gom_set(g, kposm(jobs), gid%q, gp%qf(jobs,1:gp%nflevg,:), ldsum=.true.)
  if(.not.llset) then
    llset = gom_set(g, kposm(jobs), gid%ql, sum(gp%qf(jobs,1:gp%nflevg,:),dim=1), ldsum=.true.)
  endif

  llset = gom_set(g, kposm(jobs), gid%o3, gp%o3f (jobs,1:gp%nflevg,:), ldsum=.true.)
  llset = gom_set(g, kposm(jobs), gid%l,  gp%clwf(jobs,1:gp%nflevg,:), ldsum=.true.)
  llset = gom_set(g, kposm(jobs), gid%i,  gp%ciwf(jobs,1:gp%nflevg,:), ldsum=.true.)
  llset = gom_set(g, kposm(jobs), gid%a,  gp%ccf (jobs,1:gp%nflevg,:), ldsum=.true.)

  igems = 0
  if(gom_is_interpolated(g, gid%ghg, kobtype)) then
    igems_start = igems+1
    do isubvar = 1,gom_subvar_number(g, gid%ghg)
      igems = igems+1
      llset = gom_set(g, kposm(jobs), gid%ghg, isubvar, gp%gems(jobs,1:gp%nflevg,:,igems), ldsum=.true.)  
    enddo
  endif

  if(gom_is_interpolated(g, gid%aero, kobtype)) then
    igems_start = igems+1
    do isubvar = 1,gom_subvar_number(g, gid%aero)
      igems = igems+1
      llset = gom_set(g, kposm(jobs), gid%aero, isubvar, gp%gems(jobs,1:gp%nflevg,:,igems), ldsum=.true.)  
    enddo
  endif

  if(gom_is_interpolated(g, gid%chem_assim, kobtype)) then
    igems_start = igems+1
    do isubvar = 1,gom_subvar_number(g, gid%chem_assim)
      igems = igems+1
      llset = gom_set(g, kposm(jobs), gid%chem_assim, isubvar, gp%gems(jobs,1:gp%nflevg,:,igems), ldsum=.true.)  
    enddo
  endif

  if (gom_is_interpolated(g, gid%phys, kobtype)) then
    llset = gom_set(g, kposm(jobs), gid%phys, 1, gp%clwd(jobs,1:gp%nflevg,:), ldsum=.true.)
    llset = gom_set(g, kposm(jobs), gid%phys, 2, gp%ciwd(jobs,1:gp%nflevg,:), ldsum=.true.)
    llset = gom_set(g, kposm(jobs), gid%phys, 3, gp%ccd (jobs,1:gp%nflevg,:), ldsum=.true.)
    llset = gom_set(g, kposm(jobs), gid%phys, 4, gp%rfl (jobs,1:gp%nflevg,:), ldsum=.true.)
    llset = gom_set(g, kposm(jobs), gid%phys, 5, gp%sfl (jobs,1:gp%nflevg,:), ldsum=.true.)
    llset = gom_set(g, kposm(jobs), gid%phys, 6, gp%crfl(jobs,1:gp%nflevg,:), ldsum=.true.)
    llset = gom_set(g, kposm(jobs), gid%phys, 7, gp%csfl(jobs,1:gp%nflevg,:), ldsum=.true.)
    llset = gom_set(g, kposm(jobs), gid%phys, 8, gp%pfrc(jobs,1:gp%nflevg,:), ldsum=.true.)
    llset = gom_set(g, kposm(jobs), gid%phys, 9, gp%phys_type(jobs,1:gp%nflevg,:), ldsum=.true.) !KaLo
  endif

!LLL
  if (gom_is_interpolated(g, gid%prec_accum, kobtype)) then
    llset = gom_set(g, kposm(jobs), gid%prec_accum, gp%prec_accum(jobs,:), ldsum=.true.)
  endif
!LLL

  llset = gom_set(g, kposm(jobs), gid%sp, gp%presh(jobs,gp%nflevg,:), ldsum=.true.)

enddo


if (lhook) call dr_hook('gom_plus_fill_ad',1,zhook_handle)
end subroutine gom_plus_fill_ad

! --------------------------------------------------------------
! Compute quantities derived from GOM variables
! AJGDB we have two options here: (a) keep the model isolated from the observation
! world and interpolate this stuff across in the GOMS and possible do vertical
! interpolation in the GOMS (at extra cost)  or (b) accept
! that we need to use model functions to manipulate model variables in a way
! consistent with the model, which is a key part of many osbervation operators.
! --------------------------------------------------------------
subroutine gom_plus_derived(ydephy,ydphy0,ydphy1,gp,kobtype)

use yoephy  , only : tephy
use yomphy0 , only : tphy0
use yomphy1 , only : tphy1
type(tephy)          ,intent(inout) :: ydephy
type(tphy0)          ,intent(inout) :: ydphy0
type(tphy1)          ,intent(inout) :: ydphy1
type(type_gom_plus),  intent(inout) :: gp
integer(kind=jpim),   intent(in)    :: kobtype

logical :: llvdf5
integer(kind=jpim) :: jhoriz, jlev
real(kind=jprb) :: zaphih (gp%ndlen,0:gp%nflevg,gp%nhoriz) ! geopotential on model half levels for surface fields
real(kind=jprb) :: zuf(gp%ndlen),zvf(gp%ndlen),zh10(gp%ndlen)
real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_derived',0,zhook_handle)

do jhoriz = 1,gp%nhoriz

  ! These lines are hopefully not needed at some stage.
  gp%aphi(:,:,jhoriz) = 0.0_jprb         ! Uninitialised otherwise
  gp%qn(:,:,jhoriz) = 0.0_jprb           ! Uninitialised otherwise

  ! AJGDB to bit-replicate existing code, no TL of this exists, though that is obviously wrong!
  ! (either L,I,A is never used in assimilation, not thoroughly tested, or there never are negative
  ! values in the direct)
  if (lecmwf) then
    where ( gp%clwf(:,1:gp%nflevg,jhoriz) < 0.0_jprb) gp%clwf(:,1:gp%nflevg,jhoriz) = 0.0_jprb !CPtoclean after testing (set to gom_missing_value)
    where ( gp%ciwf(:,1:gp%nflevg,jhoriz) < 0.0_jprb) gp%ciwf(:,1:gp%nflevg,jhoriz) = 0.0_jprb
    where ( gp%ccf(:,1:gp%nflevg,jhoriz) < 0.0_jprb) gp%ccf(:,1:gp%nflevg,jhoriz) = 0.0_jprb
  endif

  gp%presh(:,gp%nflevg,jhoriz) = exp(gp%presh(:,gp%nflevg,jhoriz))

  call gphpre(gp%ndlen,gp%nflevg,1,gp%ndlen,gp%vgeom%yrvab,gp%presh(:,:,jhoriz),&
    & pxyb=gp%pxyb(:,:,:,jhoriz), presf=gp%presf(:,:,jhoriz))

  call gprcp_qlirsg(gp%ndlen,1,gp%ndlen,gp%nflevg, pq=gp%qf(:,1:gp%nflevg,jhoriz),&
    & pcp=gp%cp(:,:,jhoriz), pr=gp%rf(:,:,jhoriz), pkap=gp%kap(:,:,jhoriz))

  ! snow depth calculation
  gp%sndp(:,jhoriz)=gp%sn(:,jhoriz)*0.001_JPRB

  ! geopotential calculation
  do jlev=0,gp%nflevg-1
    gp%geoph(:,jlev,jhoriz)=0.0_jprb
  enddo
  gp%geoph(:,gp%nflevg,jhoriz)=gp%orog(:,jhoriz)

  call gpgeo(gp%ndlen,1,gp%ndlen,gp%nflevg,gp%geoph(:,:,jhoriz),gp%geopf(:,:,jhoriz),&
    & gp%tf(:,1:gp%nflevg,jhoriz), gp%rf(:,:,jhoriz),&
    & gp%pxyb(:,:,yytxyb%m_lnpr,jhoriz), gp%pxyb(:,:,yytxyb%m_alph,jhoriz),&
    & gp%vgeom)

  call calc_geom_height(gp%ndlen, gp%nflevg, gp%lat(:,jhoriz),&
    & gp%geopf(:,:,jhoriz),gp%presf(:,:,jhoriz),gp%tf(:,1:gp%nflevg,jhoriz),&
    & gp%qf(:,1:gp%nflevg,jhoriz),&
    & gp%hgt(:,1:gp%nflevg,jhoriz))

  ! compution of aux. arrays for vertical interpolations
  ! AJGDB some of tthis stuff may belong more logically in a pp-preparation routine

  ! vertical interpolations in p
  call ppinit(gp%ndlen,1,gp%ndlen,gp%nflevg,gp%nppm,gp%presh(:,:,jhoriz),&
    & gp%presf(:,:,jhoriz),gp%pxp(:,:,:,jhoriz),gp%pxpd(:,:,:,jhoriz))

  ! vertical interpolations in z
  call ppinitz(gp%ndlen,1,gp%ndlen,gp%nflevg,gp%nppm,gp%geoph(:,:,jhoriz),&
    & gp%geopf(:,:,jhoriz),gp%pxz(:,:,:,jhoriz),gp%pxzd(:,:,:,jhoriz))

  ! compute t* and t0
  call ctstar(gp%ndlen,1,gp%ndlen,&
    & gp%tf(:,nlextrap,jhoriz),gp%presh(:,gp%nflevg,jhoriz),&
    & gp%presf(:,nlextrap,jhoriz),gp%orog(:,jhoriz),&
    & gp%tstar(:,jhoriz),gp%t0(:,jhoriz))
   
  ! Create Q with negative humidities removed but treated gracefully (preserving some TL/AD sensitivity)
  gp%qn(:,:,jhoriz)=gp%qf(:,:,jhoriz)
  call qneglim(gp%qn(:,:,jhoriz),1,gp%ndlen,1,gp%nflevg,1,gp%ndlen,1,gp%nflevg)

enddo

do jhoriz = 1,gp%nhoriz

  ! AJGDB I don't understand why the surface operators need their own mutant
  ! geopotential, but OK, we might as well make it:
  zaphih(:,gp%nflevg,jhoriz)=0.0_jprb
  call gpgeo(gp%ndlen,1,gp%ndlen,gp%nflevg,zaphih(:,:,jhoriz),gp%aphi(:,:,jhoriz),&
    & gp%tf(:,1:gp%nflevg,jhoriz),gp%rf(:,:,jhoriz),&
    & gp%pxyb(:,:,yytxyb%m_lnpr,jhoriz), gp%pxyb(:,:,yytxyb%m_alph,jhoriz),&
    & gp%vgeom)

enddo

! AJGDB more exceptionalism - we shouldn't need to check all these logicals...
if (lecmwf .or. (.not. lcachmt)) then

  ! AJGDB (next two lines) we need to try to remove this kind of module-driven exceptionalism. 
  ! If this stuff is strictly limited to scatterometer data it may be best part of a self-contained 
  ! scatterometer observation operator. If it is potentially useful for other obstypes, why 
  ! not just make it universally available?
  llvdf5=.false.
  llvdf5=(lvdftraj.and.(kobtype == nscatt)) ! currently, limited usage to scatt data

  do jhoriz = 1,gp%nhoriz

    ! correct for ocean current
    if (ydephy%lecurr) then
      zuf=gp%uf(:,gp%nflevg,jhoriz)-gp%u0(:,jhoriz)
      zvf=gp%vf(:,gp%nflevg,jhoriz)-gp%v0(:,jhoriz)
    else
      zuf=gp%uf(:,gp%nflevg,jhoriz)
      zvf=gp%vf(:,gp%nflevg,jhoriz)
    endif

    ! compute the roughness length over sea
    if (.not.lecmwf) then
      call rousea(gp%ndlen,1,gp%ndlen,gp%aphi(:,gp%nflevg,jhoriz),&
        & gp%uf(:,gp%nflevg,jhoriz),gp%vf(:,gp%nflevg,jhoriz),gp%ls(:,jhoriz),gp%z0(:,jhoriz))  
    endif
    if (llvdf5) then
      zh10(1:gp%ndlen)=10.0_jprb*rg
      call z0sea(gp%ndlen,1,gp%ndlen,zh10(:),&
        & gp%u10n(:,jhoriz),gp%v10n(:,jhoriz),gp%ls(:,jhoriz),gp%chn(:,jhoriz),gp%z0(:,jhoriz))  
    endif

    ! new thermodynamic variable and boundary conditions.
    call surbound(&
      & ydephy,ydphy1,gp%ndlen,1,gp%ndlen,&
      & gp%aphi(:,gp%nflevg,jhoriz),gp%qn(:,gp%nflevg,jhoriz),gp%tf(:,gp%nflevg,jhoriz),&
      & gp%ts(:,jhoriz),gp%presh(:,gp%nflevg,jhoriz),&
      & gp%ws(:,jhoriz),gp%veg(:,jhoriz),gp%wl(:,jhoriz),&
      & gp%sn(:,jhoriz),gp%ls(:,jhoriz),&
      & gp%es(:,jhoriz),gp%cptgz(:,jhoriz),gp%wet(:,jhoriz),gp%qs(:,jhoriz),gp%cpts(:,jhoriz))

    ! computation of the modified roughness lengths.
    if (llvdf5) then
      call exchco_vdf(gp%ndlen,1,gp%ndlen,&
        & zuf(:),zvf(:),gp%u10n(:,jhoriz),gp%v10n(:,jhoriz),&
        & gp%aphi(:,gp%nflevg,jhoriz),gp%z0(:,jhoriz),&
        & gp%ri(:,jhoriz),gp%bn(:,jhoriz),gp%bm(:,jhoriz))
    else
      call exchco(&
        & ydphy0,gp%ndlen,1,gp%ndlen,&
        & gp%qn(:,gp%nflevg,jhoriz),gp%tf(:,gp%nflevg,jhoriz),&
        & gp%aphi(:,gp%nflevg,jhoriz),zaphih(:,gp%nflevg,jhoriz),gp%rf(:,gp%nflevg,jhoriz),&
        & gp%ts(:,jhoriz),zuf(:),zvf(:),gp%z0(:,jhoriz),gp%wet(:,jhoriz),gp%qs(:,jhoriz),&
        & gp%ri(:,jhoriz),gp%bn(:,jhoriz),gp%bm(:,jhoriz),gp%bh(:,jhoriz))
    endif

  enddo
endif

if (lhook) call dr_hook('gom_plus_derived',1,zhook_handle)
end subroutine gom_plus_derived

! --------------------------------------------------------------
! Compute quantities derived from GOM variables - TL
! --------------------------------------------------------------
subroutine gom_plus_derived_tl(ydephy,ydphy0,ydphy1,gp,gp5,kobtype)
 
use yoephy , only : tephy
use yomphy0 , only : tphy0
use yomphy1 , only : tphy1
type(tephy)         ,intent(inout) :: ydephy
type(tphy0)         ,intent(inout) :: ydphy0
type(tphy1)         ,intent(inout) :: ydphy1
type(type_gom_plus), intent(inout) :: gp
type(type_gom_plus), intent(in)    :: gp5
integer(kind=jpim),  intent(in)    :: kobtype

logical :: llvdf5
logical :: llun5(gp%ndlen) ! Allow perturbations in UN10
integer(kind=jpim) :: jhoriz, jlev
real(kind=jprb), dimension(gp%ndlen,0:gp%nflevg) :: zaphih, zaphih5 ! geopotential for surface fields
real(kind=jprb), dimension(gp%ndlen,gp%nflevg) :: zaphif5  
real(kind=jprb), dimension(gp%ndlen) :: zuf,zvf,zh10,zuf5,zvf5,zh105

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_derived_tl',0,zhook_handle)

do jhoriz=1,gp%nhoriz

  ! AJGDB to bit-replicate prior code, no TL of this direct code exists, though that is obviously wrong!
  ! (either L,I,A is never used in assimilation, not thoroughly tested, or there never are negative 
  ! values in the direct)
  !where ( gp%clwf(:,1:gp%nflevg,jhoriz) < 0.0_jprb) gp%clwf(:,1:gp%nflevg,jhoriz) = 0.0_jprb
  !where ( gp%ciwf(:,1:gp%nflevg,jhoriz) < 0.0_jprb) gp%ciwf(:,1:gp%nflevg,jhoriz) = 0.0_jprb
  !where ( gp%ccf(:,1:gp%nflevg,jhoriz) < 0.0_jprb) gp%ccf(:,1:gp%nflevg,jhoriz) = 0.0_jprb

  gp%presh(:,gp%nflevg,jhoriz) = gp5%presh(:,gp%nflevg,jhoriz)*gp%presh(:,gp%nflevg,jhoriz)

  call gphpretl(gp%ndlen,gp%nflevg,1,gp%ndlen,gp%vgeom%yrvab,gp%presh(:,:,jhoriz),gp5%presh(:,:,jhoriz),&
    & pxyb=gp%pxyb(:,:,:,jhoriz),pxyb5=gp5%pxyb(:,:,:,jhoriz),presf=gp%presf(:,:,jhoriz))

  call gprcptl (gp%ndlen,1,gp%ndlen,gp%nflevg,pq=gp%qf(:,1:gp%nflevg,jhoriz),&
              &pcp=gp%cp(:,:,jhoriz),pr=gp%rf(:,:,jhoriz),pkap=gp%kap(:,:,jhoriz),&
              &pcp5=gp5%cp(:,:,jhoriz),pr5=gp5%rf(:,:,jhoriz))

  ! geop. calculation
  do jlev=0,gp%nflevg
    gp%geoph(1:gp%ndlen,jlev,jhoriz)=0.0_jprb
  enddo

  call gpgeotl (gp%ndlen,1,gp%ndlen,gp%nflevg,&
    & gp%geoph(:,:,jhoriz),gp%geopf(:,:,jhoriz),&
    & gp%tf(:,1:gp%nflevg,jhoriz),gp%rf(:,:,jhoriz),&
    & gp%pxyb(:,:,yytxyb%m_lnpr,jhoriz),gp%pxyb(:,:,yytxyb%m_alph,jhoriz),&
    & gp5%tf(:,1:gp%nflevg,jhoriz),gp5%rf(:,:,jhoriz),&
    & gp5%pxyb(:,:,yytxyb%m_lnpr,jhoriz),gp5%pxyb(:,:,yytxyb%m_alph,jhoriz),gp%vgeom)

  call calc_geom_height_tl(gp%ndlen, gp%nflevg, gp5%lat(:,jhoriz),&
    & gp%geopf(:,:,jhoriz),gp%presf(:,:,jhoriz),gp%tf(:,1:gp%nflevg,jhoriz),&
    & gp%qf(:,1:gp%nflevg,jhoriz),&
    & gp%hgt(:,1:gp%nflevg,jhoriz),&
    & gp5%geopf(:,:,jhoriz),gp5%presf(:,:,jhoriz),gp5%tf(:,1:gp%nflevg,jhoriz),&
    & gp5%qf(:,1:gp%nflevg,jhoriz))

  ! compution of aux. arrays for vertical interpolations

  ! vertical interpolations in p
  call ppinittl (gp%ndlen,1,gp%ndlen,gp%nflevg,gp%nppm,gp%presh(:,:,jhoriz),gp%presf(:,:,jhoriz),&
    & gp%pxp(:,:,:,jhoriz),gp%pxpd(:,:,:,jhoriz),gp5%presh(:,:,jhoriz),gp5%presf(:,:,jhoriz),&
    & gp5%pxp(:,:,:,jhoriz))

  ! vertical interpolations in z
  call ppinitztl (gp%ndlen,1,gp%ndlen,gp%nflevg,gp%nppm,&
    & gp%geoph(:,:,jhoriz),gp%geopf(:,:,jhoriz),&
    & gp%pxz(:,:,:,jhoriz),gp%pxzd(:,:,:,jhoriz),&
    & gp5%pxz(:,:,:,jhoriz))  

  ! compute t* and t0
  call ctstartl (gp%ndlen,1,gp%ndlen,&
    & gp%tf(:,nlextrap,jhoriz),gp%presh(:,gp%nflevg,jhoriz),gp%presf(:,nlextrap,jhoriz),&
    & gp%tstar(:,jhoriz),gp%t0(:,jhoriz),&
    & gp5%tf(:,nlextrap,jhoriz),gp5%presh(:,gp%nflevg,jhoriz),&
    & gp5%presf(:,nlextrap,jhoriz))  

  ! Create Q with negative humidities removed but treated gracefully (preserving some TL/AD sensitivity)    
  gp%qn(:,:,jhoriz)=gp%qf(:,:,jhoriz)
  call qneglimtl(gp%qn(:,:,jhoriz),gp5%qf(:,:,jhoriz),1,gp%ndlen,1,gp%nflevg,1,gp%ndlen,1,gp%nflevg)
  
enddo

! AJGDB more exceptionalism - we shouldn't need to check all these logicals...
if(lecmwf .or. (.not. lcachmt)) then

  ! AJGDB (next two lines) we need to try to remove this kind of module-driven exceptionalism. 
  ! If this stuff is strictly limited to scatterometer data it may be best part of a self-contained 
  ! scatterometer observation operator. If it is potentially useful for other obstypes, why 
  ! not just make it universally available?
  llvdf5=.false.
  llvdf5=(lvdftraj.and.(kobtype == nscatt)) ! currently, limited usage to scatt data

  do jhoriz=1,gp%nhoriz

    ! AJGDB that extra surface geop. calculation again. Is this really necessary given what we do above?
    zaphih5(:,gp%nflevg)=0.0_jprb
    zaphih (:,gp%nflevg)=0.0_jprb
    zaphif5 (:,gp%nflevg)=0.0_jprb ! zaphih5 needs to be recomputed; gp5 is IN only
    call gpgeo(gp%ndlen,1,gp%ndlen,gp%nflevg,zaphih5(:,:),zaphif5(:,:),&
      & gp5%tf(:,1:gp%nflevg,jhoriz),gp5%rf(:,:,jhoriz),gp5%pxyb(:,:,yytxyb%m_lnpr,jhoriz),&
      & gp5%pxyb(:,:,yytxyb%m_alph,jhoriz),&
      & gp%vgeom)
    call gpgeotl(gp%ndlen,1,gp%ndlen,gp%nflevg,&
      & zaphih(:,:),gp%aphi(:,:,jhoriz),&
      & gp%tf(:,1:gp%nflevg,jhoriz),gp%rf(:,:,jhoriz),gp%pxyb(:,:,yytxyb%m_lnpr,jhoriz),&
      & gp%pxyb(:,:,yytxyb%m_alph,jhoriz),gp5%tf(:,1:gp%nflevg,jhoriz),gp5%rf(:,:,jhoriz),&
      & gp5%pxyb(:,:,yytxyb%m_lnpr,jhoriz),gp5%pxyb(:,:,yytxyb%m_alph,jhoriz),gp%vgeom)  
    ! correct for ocean current
    if (ydephy%lecurr) then
      zuf5=gp5%uf(:,gp%nflevg,jhoriz)-gp5%u0(:,jhoriz)
      zvf5=gp5%vf(:,gp%nflevg,jhoriz)-gp5%v0(:,jhoriz)
    else
      zuf5=gp5%uf(:,gp%nflevg,jhoriz)
      zvf5=gp5%vf(:,gp%nflevg,jhoriz)
    endif
    zuf=gp%uf(:,gp%nflevg,jhoriz)
    zvf=gp%vf(:,gp%nflevg,jhoriz)

    ! compute the roughness length over sea
    if(.not.lecmwf) then
      call rouseatl(gp%ndlen,1,gp%ndlen,&
        & gp5%aphi(:,gp%nflevg,jhoriz),gp5%uf(:,gp%nflevg,jhoriz),gp5%vf(:,gp%nflevg,jhoriz),gp5%ls(:,jhoriz),&
        & gp%aphi (:,gp%nflevg,jhoriz),gp%uf (:,gp%nflevg,jhoriz),gp%vf (:,gp%nflevg,jhoriz),&
        & gp%z0(:,jhoriz))  
    endif
    if (llvdf5) then
      zh105=10.0_jprb*rg
      zh10 = 0.0_jprb
      call z0seatl(gp%ndlen,1,gp%ndlen,&
        & zh105,gp5%u10n(:,jhoriz),gp5%v10n(:,jhoriz),gp5%ls(:,jhoriz),gp5%chn(:,jhoriz),&
        & zh10 ,gp%u10n(:,jhoriz) ,gp%v10n(:,jhoriz),&
        & gp%z0(:,jhoriz))
    endif

    ! new thermodynamic variable and boundary conditions.
    call surboundtl(ydephy,ydphy1,gp%ndlen,1,gp%ndlen,&
      & gp5%qn(:,gp%nflevg,jhoriz),gp5%tf(:,gp%nflevg,jhoriz),gp5%ts(:,jhoriz),gp5%presh(:,gp%nflevg,jhoriz),&
      & gp%qn(:,gp%nflevg,jhoriz),gp%tf(:,gp%nflevg,jhoriz),gp%ts(:,jhoriz),gp%presh(:,gp%nflevg,jhoriz),&
      & gp%aphi(:,gp%nflevg,jhoriz),&
      & gp5%ws(:,jhoriz),gp5%veg(:,jhoriz),gp5%wl(:,jhoriz),gp5%sn(:,jhoriz),gp5%ls(:,jhoriz),&
      & gp%ws(:,jhoriz),gp%veg(:,jhoriz),gp%wl(:,jhoriz),gp%sn(:,jhoriz),&
      & gp5%es(:,jhoriz),gp5%wet(:,jhoriz),gp5%qs(:,jhoriz),&
      & gp%wet(:,jhoriz),gp%qs(:,jhoriz),gp%cptgz(:,jhoriz),gp%cpts(:,jhoriz))

    ! computation of the modified roughness lengths.
    if(llvdf5) then
      if(any(gp5%upd(:,jhoriz)==gom_missing_value)) then
        llun5=.false.
      else
        llun5=lvdfmin.and.(gp5%upd(:,jhoriz) > 0.5_jprb)
      endif
      call exchco_vdftl(gp%ndlen,1,gp%ndlen,llun5,&
        & zuf5,zvf5,gp5%u10n(:,jhoriz),gp5%v10n(:,jhoriz),&
        & gp5%aphi(:,gp%nflevg,jhoriz),gp5%z0(:,jhoriz),&
        & zuf,zvf,gp%u10n(:,jhoriz) ,gp%v10n(:,jhoriz) ,&
        & gp%aphi (:,gp%nflevg,jhoriz),gp%z0(:,jhoriz) ,&
        & gp%bn(:,jhoriz) ,gp%bm(:,jhoriz))
    else
      call exchcotl(ydphy0,gp%ndlen,1,gp%ndlen,&
        & gp5%qn(:,gp%nflevg,jhoriz),gp5%tf(:,gp%nflevg,jhoriz),&
        & gp5%aphi(:,gp%nflevg,jhoriz),zaphih5(:,gp%nflevg),gp5%rf(:,gp%nflevg,jhoriz),&
        & gp5%ts(:,jhoriz),gp%qn(:,gp%nflevg,jhoriz),gp%tf(:,gp%nflevg,jhoriz),&
        & gp%aphi(:,gp%nflevg,jhoriz),zaphih(:,gp%nflevg),gp%rf(:,gp%nflevg,jhoriz),gp%ts(:,jhoriz),&
        & zuf5,zvf5,gp5%z0(:,jhoriz),&
        & zuf,zvf,gp%z0(:,jhoriz),&
        & gp5%wet(:,jhoriz),gp5%qs(:,jhoriz),&
        & gp%wet(:,jhoriz),gp%qs(:,jhoriz),&
        & gp5%ri(:,jhoriz),&
        & gp%bn(:,jhoriz),gp%bm(:,jhoriz),gp%bh(:,jhoriz))
    endif

  enddo
endif

if (lhook) call dr_hook('gom_plus_derived_tl',1,zhook_handle)
end subroutine gom_plus_derived_tl

! --------------------------------------------------------------
! Compute quantities derived from GOM variables - Adjoint
! --------------------------------------------------------------
subroutine gom_plus_derived_ad(ydephy,ydphy0,ydphy1,gp,gp5,kobtype)

use yoephy , only : tephy
use yomphy0 , only : tphy0
use yomphy1 , only : tphy1
type(tephy)         ,intent(inout) :: ydephy
type(tphy0)         ,intent(inout) :: ydphy0
type(tphy1)         ,intent(inout) :: ydphy1
type(type_gom_plus), intent(inout) :: gp
type(type_gom_plus), intent(in)    :: gp5
integer(kind=jpim),  intent(in)    :: kobtype

logical :: llvdf5
logical :: llun5(gp%ndlen) ! Allow perturbations in UN10
integer(kind=jpim) :: jhoriz, jlev
real(kind=jprb), dimension(gp%ndlen,0:gp%nflevg) :: zaphih, zaphih5 ! geopotential on model levels for surface fields
real(kind=jprb), dimension(gp%ndlen,gp%nflevg) :: zaphif5  
real(kind=jprb), dimension(gp%ndlen) :: zuf,zvf,zh10,zuf5,zvf5,zh105,zqs,zwet

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_derived_ad',0,zhook_handle)

if(lecmwf .or. (.not. lcachmt)) then

  ! AJGDB (next two lines) we need to try to remove this kind of module-driven exceptionalism. 
  ! If this stuff is strictly limited to scatterometer data it may be best part of a self-contained 
  ! scatterometer observation operator. If it is potentially useful for other obstypes, why 
  ! not just make it universally available?
  llvdf5=.false.
  llvdf5=(lvdftraj.and.(kobtype == nscatt)) ! currently, limited usage to scatt data

  do jhoriz=1,gp%nhoriz

    zqs (1:gp%ndlen)=0.0_jprb
    zwet(1:gp%ndlen)=0.0_jprb
    zuf (1:gp%ndlen)=0.0_jprb
    zvf (1:gp%ndlen)=0.0_jprb
    zh10(1:gp%ndlen)=0.0_jprb
    zaphih5(1:gp%ndlen,gp%nflevg)=0.0_jprb
    zaphih(1:gp%ndlen,0:gp%nflevg)=0.0_jprb
    zaphif5(1:gp%ndlen,gp%nflevg)=0.0_jprb

    ! Correct for ocean current
    if (ydephy%lecurr) then
      zuf5(1:gp%ndlen)=gp5%uf(1:gp%ndlen,gp%nflevg,jhoriz)-gp5%u0(1:gp%ndlen,jhoriz)
      zvf5(1:gp%ndlen)=gp5%vf(1:gp%ndlen,gp%nflevg,jhoriz)-gp5%v0(1:gp%ndlen,jhoriz)
    else
      zuf5(1:gp%ndlen)=gp5%uf(1:gp%ndlen,gp%nflevg,jhoriz)
      zvf5(1:gp%ndlen)=gp5%vf(1:gp%ndlen,gp%nflevg,jhoriz)
    endif

    call gpgeo(gp%ndlen,1,gp%ndlen,gp%nflevg,zaphih5(:,:),zaphif5(:,:),&
      & gp5%tf(:,1:gp%nflevg,jhoriz),gp5%rf(:,:,jhoriz),gp5%pxyb(:,:,yytxyb%m_lnpr,jhoriz),&
      & gp5%pxyb(:,:,yytxyb%m_alph,jhoriz),&
      &gp%vgeom)

    ! Thermodynamic variables and drag coefficients
    if(llvdf5) then
      if(any(gp5%upd(:,jhoriz)==gom_missing_value)) then
        llun5=.false.
      else
        llun5=lvdfmin.and.(gp5%upd(:,jhoriz) > 0.5_jprb)
      endif
      call exchco_vdfad(gp%ndlen,1,gp%ndlen,llun5,&
        & zuf5,zvf5,gp5%u10n(:,jhoriz),gp5%v10n(:,jhoriz),&
        & zaphif5(:,gp%nflevg),gp5%z0(:,jhoriz),&
        & zuf,zvf,gp%u10n(:,jhoriz) ,gp%v10n(:,jhoriz) ,&
        & gp%aphi (:,gp%nflevg,jhoriz),gp%z0(:,jhoriz) ,&
        & gp%bn(:,jhoriz) ,gp%bm(:,jhoriz))
    else
      call exchcoad(ydphy0,gp%ndlen,1,gp%ndlen,&
        & gp5%qn(:,gp%nflevg,jhoriz),gp5%tf(:,gp%nflevg,jhoriz),&
        & zaphif5(:,gp%nflevg),zaphih5(:,gp%nflevg),gp5%rf(:,gp%nflevg,jhoriz),&
        & gp5%ts(:,jhoriz),gp%qn(:,gp%nflevg,jhoriz),gp%tf(:,gp%nflevg,jhoriz),&
        & gp%aphi(:,gp%nflevg,jhoriz),zaphih(:,gp%nflevg),gp%rf(:,gp%nflevg,jhoriz),gp%ts(:,jhoriz),&
        & zuf5(:),zvf5(:),gp5%z0(:,jhoriz),&
        & gp%uf(:,gp%nflevg,jhoriz),gp%vf(:,gp%nflevg,jhoriz),gp%z0(:,jhoriz),&
        & gp5%wet(:,jhoriz),gp5%qs(:,jhoriz),&
        & gp%wet(:,jhoriz),gp%qs(:,jhoriz),&
        & gp5%ri(:,jhoriz),&
        & gp%bn(:,jhoriz),gp%bm(:,jhoriz),gp%bh(:,jhoriz))
    endif

    call surboundad(ydephy,ydphy1,gp%ndlen,1,gp%ndlen,&
      & gp5%qn(:,gp%nflevg,jhoriz),gp5%tf(:,gp%nflevg,jhoriz),gp5%ts(:,jhoriz),gp5%presh(:,gp%nflevg,jhoriz),&
      & gp%qn(:,gp%nflevg,jhoriz),gp%tf(:,gp%nflevg,jhoriz),gp%ts(:,jhoriz),gp%presh(:,gp%nflevg,jhoriz),&
      & gp%aphi(:,gp%nflevg,jhoriz),&
      & gp5%ws(:,jhoriz),gp5%veg(:,jhoriz),gp5%wl(:,jhoriz),gp5%sn(:,jhoriz),gp5%ls(:,jhoriz),&
      & gp%ws(:,jhoriz),gp%veg(:,jhoriz),gp%wl(:,jhoriz),gp%sn(:,jhoriz),&
      & gp5%es(:,jhoriz),gp5%wet(:,jhoriz),gp5%qs(:,jhoriz),&
      & gp%wet(:,jhoriz),gp%qs(:,jhoriz),&
      & gp%cptgz(:,jhoriz),gp%cpts(:,jhoriz))

    ! Adjoint of the roughness length over sea
    if(.not.lecmwf) then
      call rouseaad(gp%ndlen,1,gp%ndlen,&
        & zaphif5(:,gp%nflevg),gp5%uf(:,gp%nflevg,jhoriz),gp5%vf(:,gp%nflevg,jhoriz),gp5%ls(:,jhoriz),&
        & gp%aphi(:,gp%nflevg,jhoriz),gp%uf(:,gp%nflevg,jhoriz),gp%vf(:,gp%nflevg,jhoriz),&
        & gp%z0(:,jhoriz))  
    endif
    if (llvdf5) then
      zh105(1:gp%ndlen)=10.0_jprb*rg
      call z0seaad(gp%ndlen,1,gp%ndlen,&
        & zh105,gp5%u10n(:,jhoriz),gp5%v10n(:,jhoriz),gp5%ls(:,jhoriz),gp5%chn(:,jhoriz),&
        & zh10,gp%u10n(:,jhoriz),gp%v10n(:,jhoriz) ,&
        & gp%z0(:,jhoriz))  
      zh10 (1:gp%ndlen)= 0.0_jprb
    endif

    ! Correct for ocean current
    ! AJGDB I am suspicious this (copied from old preint routines) is not TL/AD correct,
    ! or is at least confusing - see the TL equivalent is keyed on lecurr.
    if (llvdf5) then
      gp%uf(1:gp%ndlen,gp%nflevg,jhoriz)=gp%uf(1:gp%ndlen,gp%nflevg,jhoriz)+zuf(1:gp%ndlen)
      gp%vf(1:gp%ndlen,gp%nflevg,jhoriz)=gp%vf(1:gp%ndlen,gp%nflevg,jhoriz)+zvf(1:gp%ndlen)
    endif

    ! Adjoint of of geopotential on model levels
    call gpgeoad(gp%ndlen,1,gp%ndlen,gp%nflevg,zaphih,gp%aphi(:,:,jhoriz),&
      & gp%tf(:,1:gp%nflevg,jhoriz),gp%rf(:,:,jhoriz),&
      & gp%pxyb(:,:,yytxyb%m_lnpr,jhoriz),gp%pxyb(:,:,yytxyb%m_alph,jhoriz),gp5%tf(:,1:gp%nflevg,jhoriz),&
      & gp5%rf(:,:,jhoriz),gp5%pxyb(:,:,yytxyb%m_lnpr,jhoriz),gp5%pxyb(:,:,yytxyb%m_alph,jhoriz),gp%vgeom)  

  enddo
endif

do jhoriz=1,gp%nhoriz

  ! Create Q with negative humidities removed but treated gracefully (preserving some TL/AD sensitivity)    
  call qneglimtl(gp%qn(:,:,jhoriz),gp5%qf(:,:,jhoriz),1,gp%ndlen,1,gp%nflevg,1,gp%ndlen,1,gp%nflevg)
  gp%qf(:,:,jhoriz) = gp%qf(:,:,jhoriz) + gp%qn(:,:,jhoriz)
  
  call ctstarad(gp%ndlen,1,gp%ndlen,&
    & gp%tf(:,nlextrap,jhoriz),gp%presh(:,gp%nflevg,jhoriz),gp%presf(:,nlextrap,jhoriz),&
    & gp%tstar(:,jhoriz),gp%t0(:,jhoriz),&
    & gp5%tf(:,nlextrap,jhoriz),gp5%presh(:,gp%nflevg,jhoriz),&
    & gp5%presf(:,nlextrap,jhoriz))  

  call ppinitza(gp%ndlen,1,gp%ndlen,gp%nflevg,gp%nppm,&
    & gp%geoph(:,:,jhoriz),gp%geopf(:,:,jhoriz),&
    & gp%pxz(:,:,:,jhoriz),gp%pxzd(:,:,:,jhoriz),&
    & gp5%pxz(:,:,:,jhoriz))  

  call ppinitad(gp%ndlen,1,gp%ndlen,gp%nflevg,gp%nppm,gp%presh(:,:,jhoriz),gp%presf(:,:,jhoriz),&
    & gp%pxp(:,:,:,jhoriz),gp%pxpd(:,:,:,jhoriz),gp5%presh(:,:,jhoriz),gp5%presf(:,:,jhoriz),&
    & gp5%pxp(:,:,:,jhoriz))

  call calc_geom_height_ad(gp%ndlen, gp%nflevg, gp5%lat(:,jhoriz),&
    & gp%geopf(:,:,jhoriz),gp%presf(:,:,jhoriz),gp%tf(:,1:gp%nflevg,jhoriz),&
    & gp%qf(:,1:gp%nflevg,jhoriz),&
    & gp%hgt(:,1:gp%nflevg,jhoriz),&
    & gp5%geopf(:,:,jhoriz),gp5%presf(:,:,jhoriz),gp5%tf(:,1:gp%nflevg,jhoriz),&
    & gp5%qf(:,1:gp%nflevg,jhoriz))

  call gpgeoad(gp%ndlen,1,gp%ndlen,gp%nflevg,&
    & gp%geoph(:,:,jhoriz),gp%geopf(:,:,jhoriz),&
    & gp%tf(:,1:gp%nflevg,jhoriz),gp%rf(:,:,jhoriz),&
    & gp%pxyb(:,:,yytxyb%m_lnpr,jhoriz),gp%pxyb(:,:,yytxyb%m_alph,jhoriz),&
    & gp5%tf(:,1:gp%nflevg,jhoriz),gp5%rf(:,:,jhoriz),&
    & gp5%pxyb(:,:,yytxyb%m_lnpr,jhoriz),gp5%pxyb(:,:,yytxyb%m_alph,jhoriz), gp%vgeom)

  call gprcpad (gp%ndlen,1,gp%ndlen,gp%nflevg,pq=gp%qf(:,1:gp%nflevg,jhoriz),&
    & pcp=gp%cp(:,:,jhoriz),pr=gp%rf(:,:,jhoriz),&
    & pkap=gp%kap(:,:,jhoriz),pcp5=gp5%cp(:,:,jhoriz),pr5=gp5%rf(:,:,jhoriz))

  call gphpread(gp%ndlen,gp%nflevg,1,gp%ndlen,gp%vgeom%yrvab,gp%presh(:,:,jhoriz),gp5%presh(:,:,jhoriz),&
    & pxyb=gp%pxyb(:,:,:,jhoriz),pxyb5=gp5%pxyb(:,:,:,jhoriz),presf=gp%presf(:,:,jhoriz))

  gp%presh(:,gp%nflevg,jhoriz)=gp5%presh(:,gp%nflevg,jhoriz)*gp%presh(:,gp%nflevg,jhoriz)

enddo

if (lhook) call dr_hook('gom_plus_derived_ad',1,zhook_handle)
end subroutine gom_plus_derived_ad

! -------------------------------------------------------------
! Allocate memory

! AJGDB allocate everything, just like in hop where these were automatic variables
! for the future we might want only to allocate what we need, looking at
! the GOM first.

! -------------------------------------------------------------
subroutine gom_plus_alloc(gp, kdlen, kflevg, khoriz, kppm, kxyb, &
 & ldsrg, ldhydro_sum, ldphys, ldcanari, kgems, klev_soil)

type(type_gom_plus), intent(inout) :: gp
integer(kind=jpim),  intent(in)    :: kdlen, kflevg, khoriz, kppm, kxyb
logical, intent(in)                :: ldsrg, ldhydro_sum, ldphys, ldcanari
integer(kind=jpim), intent(in)     :: kgems 
integer(kind=jpim), intent(in)     :: klev_soil

real(kind=jprb)    :: zhook_handle
integer(kind=jpim) :: insfc, inua0, inua1, inppvar, insoil

if (lhook) call dr_hook('gom_plus_alloc',0,zhook_handle)

! Always attempt to deallocate first, as the gom_plus can be re-used multiple times
! in the kset loop. AJGDB: future performance saving: _reuse_ storage if at all 
! possible
call gom_plus_destroy(gp)

! Number of components must be managed manually by people adding new
! fields to the gom_plus: (AJGDB - all this could be done more elegantly 
! and safely in future using "fields" code to manage these pointers)
insfc = 66
if (lelam) insfc = insfc + 2
allocate(gp%store_sfc(kdlen, khoriz, insfc))

gp%orog       => gp%store_sfc(:,:,1)
gp%al         => gp%store_sfc(:,:,2)
gp%tstar      => gp%store_sfc(:,:,3)
gp%t0         => gp%store_sfc(:,:,4)
gp%ts         => gp%store_sfc(:,:,5)
gp%ls         => gp%store_sfc(:,:,6)
gp%ci         => gp%store_sfc(:,:,7)
gp%ws         => gp%store_sfc(:,:,8)
gp%sn         => gp%store_sfc(:,:,9)
gp%u10n       => gp%store_sfc(:,:,10)
gp%v10n       => gp%store_sfc(:,:,11)
gp%u0         => gp%store_sfc(:,:,12)
gp%v0         => gp%store_sfc(:,:,13)
gp%chn        => gp%store_sfc(:,:,14)
gp%tcls       => gp%store_sfc(:,:,15)
gp%hucls      => gp%store_sfc(:,:,16)
gp%z0         => gp%store_sfc(:,:,17)
gp%wl         => gp%store_sfc(:,:,18)
gp%sst        => gp%store_sfc(:,:,19)
gp%veg        => gp%store_sfc(:,:,20)
gp%cptgz      => gp%store_sfc(:,:,21)
gp%wet        => gp%store_sfc(:,:,22)
gp%qs         => gp%store_sfc(:,:,23)
gp%cpts       => gp%store_sfc(:,:,24)  
gp%ri         => gp%store_sfc(:,:,25)  
gp%bn         => gp%store_sfc(:,:,26)  
gp%bm         => gp%store_sfc(:,:,27)  
gp%bh         => gp%store_sfc(:,:,28) 
gp%upd        => gp%store_sfc(:,:,29) 
gp%arg        => gp%store_sfc(:,:,30) 
gp%hv         => gp%store_sfc(:,:,31) 
gp%z0h        => gp%store_sfc(:,:,32) 
gp%wsi        => gp%store_sfc(:,:,33) 
gp%rrr        => gp%store_sfc(:,:,34) 
gp%rsn        => gp%store_sfc(:,:,35) 
gp%es         => gp%store_sfc(:,:,36) 
gp%lat        => gp%store_sfc(:,:,37) 
gp%lon        => gp%store_sfc(:,:,38) 
gp%ucls       => gp%store_sfc(:,:,39) 
gp%vcls       => gp%store_sfc(:,:,40) 
gp%nucls      => gp%store_sfc(:,:,41) 
gp%nvcls      => gp%store_sfc(:,:,42) 
gp%prec_accum => gp%store_sfc(:,:,43) 
gp%sdfor      => gp%store_sfc(:,:,44) 
gp%t2m        => gp%store_sfc(:,:,45) 
gp%d2m        => gp%store_sfc(:,:,46) 
gp%q2m        => gp%store_sfc(:,:,47) 
gp%u10m       => gp%store_sfc(:,:,48) 
gp%v10m       => gp%store_sfc(:,:,49) 
gp%soty       => gp%store_sfc(:,:,50)
gp%lail       => gp%store_sfc(:,:,51)
gp%laih       => gp%store_sfc(:,:,52)
gp%sndt       => gp%store_sfc(:,:,53)
gp%sndp       => gp%store_sfc(:,:,54)
gp%cvl        => gp%store_sfc(:,:,55)
gp%cvh        => gp%store_sfc(:,:,56)
gp%tvl        => gp%store_sfc(:,:,57)
gp%tvh        => gp%store_sfc(:,:,58)
gp%tlk        => gp%store_sfc(:,:,59)
gp%timestep   => gp%store_sfc(:,:,60)
!LLL
!canari
gp%cori       => gp%store_sfc(:,:,61) 
gp%est        => gp%store_sfc(:,:,62) 
gp%esn        => gp%store_sfc(:,:,63) 
gp%et2        => gp%store_sfc(:,:,64) 
gp%eh2        => gp%store_sfc(:,:,65) 
gp%ev1        => gp%store_sfc(:,:,66) 
!canari
!compas
if (lelam) then
  gp%gnordl     => gp%store_sfc(:,:,67)
  gp%gnordm     => gp%store_sfc(:,:,68)
endif

inua0 = 2
allocate(gp%store_ua0(kdlen,0:kflevg, khoriz, inua0))

inua1 = 16
if( count( (/ldphys,ldsrg,ldhydro_sum,ldcanari/) ) >1 ) &
 &  call abor1('gom_plus: gom has got mutually incompatible settings')
if(ldsrg) inua1=inua1+3
if(ldhydro_sum) inua1=inua1+4
if(ldphys) inua1=inua1+9 !8 KaLo
if(ldcanari) inua1=inua1+2
allocate(gp%store_ua1(kdlen, kflevg, khoriz, inua1))

gp%presf => gp%store_ua1(:,:,:,1)
gp%geopf => gp%store_ua1(:,:,:,2)
gp%aphi  => gp%store_ua1(:,:,:,3)
gp%cp    => gp%store_ua1(:,:,:,4)
gp%rf    => gp%store_ua1(:,:,:,5)
gp%kap   => gp%store_ua1(:,:,:,6)

gp%uf  => gp%store_ua1(:,:,:,7)
gp%vf  => gp%store_ua1(:,:,:,8)
gp%tf  => gp%store_ua1(:,:,:,9)
gp%qf  => gp%store_ua1(:,:,:,10)
gp%qn  => gp%store_ua1(:,:,:,11)
gp%o3f => gp%store_ua1(:,:,:,12)
gp%clwf=> gp%store_ua1(:,:,:,13)
gp%ciwf=> gp%store_ua1(:,:,:,14)
gp%ccf => gp%store_ua1(:,:,:,15)
gp%hgt => gp%store_ua1(:,:,:,16)
if (ldsrg) then
  gp%cswf=> gp%store_ua1(:,:,:,17)
  gp%crwf=> gp%store_ua1(:,:,:,18)
  gp%cgwf=> gp%store_ua1(:,:,:,19)
endif
if (ldhydro_sum) then
  gp%clw_sum=> gp%store_ua1(:,:,:,17)
  gp%ciw_sum=> gp%store_ua1(:,:,:,18)
  gp%crw_sum=> gp%store_ua1(:,:,:,19)
  gp%csw_sum=> gp%store_ua1(:,:,:,20)
endif
if (ldphys) then
  gp%clwd=> gp%store_ua1(:,:,:,17)
  gp%ciwd=> gp%store_ua1(:,:,:,18)
  gp%ccd => gp%store_ua1(:,:,:,19)
  gp%rfl => gp%store_ua1(:,:,:,20)
  gp%sfl => gp%store_ua1(:,:,:,21)
  gp%crfl=> gp%store_ua1(:,:,:,22)
  gp%csfl=> gp%store_ua1(:,:,:,23)
  gp%pfrc=> gp%store_ua1(:,:,:,24)
  gp%phys_type=> gp%store_ua1(:,:,:,25) !KaLo
endif
if (ldcanari) then
  gp%eh  => gp%store_ua1(:,:,:,17) 
  gp%ez  => gp%store_ua1(:,:,:,18) 
endif

! AJGDB p296 of Metcalf et al. suggests that pointer assignment should preserve the
! array bounds of the target. In practice it does not and it has to be specified
! explicitly with (1:,0:,1:). Is this a compiler bug?
! AJGDB2 Could be perhaps be done (neater?) with real, dimension(1:,0:,1:)::presh??
! AJGDB3 also most of the 0:kflevg dimensions are for the convenience of the
! vertical interpolation, not the user!

gp%presh(1:,0:,1:)=> gp%store_ua0(:,:,:,1)
gp%geoph(1:,0:,1:)=> gp%store_ua0(:,:,:,2)

! store_soil
insoil = 2
allocate(gp%store_soil(kdlen, klev_soil, khoriz, insoil))
gp%tsoil=> gp%store_soil(:,:,:,1)
gp%qsoil=> gp%store_soil(:,:,:,2)


inppvar = 4
allocate(gp%store_ppm(kdlen,0:kflevg, kppm, khoriz, inppvar))
gp%pxp (1:,0:,1:,1:)=> gp%store_ppm(:,:,:,:,1)
gp%pxpd(1:,0:,1:,1:)=> gp%store_ppm(:,:,:,:,2)
gp%pxz (1:,0:,1:,1:)=> gp%store_ppm(:,:,:,:,3)
gp%pxzd(1:,0:,1:,1:)=> gp%store_ppm(:,:,:,:,4)

allocate(gp%pxyb(kdlen, kflevg, kxyb, khoriz))

allocate(gp%gems(kdlen,kflevg,khoriz,kgems))
allocate(gp%gems_igrib(kgems))
allocate(gp%gems_type(kgems))


!Allocate the vertical geometry type; allocations inside type happen in the _fill(_ad) routine
allocate(gp%vgeom)


if (lhook) call dr_hook('gom_plus_alloc',1,zhook_handle)
end subroutine gom_plus_alloc

! -------------------------------------------------------------
! Deallocate memory and nullify dangling pointers, in case
! the gom plus is re-used.
! -------------------------------------------------------------
subroutine gom_plus_destroy(gp)

type(type_gom_plus), intent(inout) :: gp

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_destroy',0,zhook_handle)

if(associated(gp%store_sfc)) deallocate(gp%store_sfc)
if(associated(gp%store_ua0)) deallocate(gp%store_ua0)
if(associated(gp%store_ua1)) deallocate(gp%store_ua1)
if(associated(gp%store_soil))deallocate(gp%store_soil)
if(associated(gp%store_ppm)) deallocate(gp%store_ppm)
if(associated(gp%gems))      deallocate(gp%gems)
if(associated(gp%gems_igrib))deallocate(gp%gems_igrib)
if(associated(gp%gems_type)) deallocate(gp%gems_type)
if(associated(gp%pxyb))      deallocate(gp%pxyb)

if (associated(gp%vgeom) ) then
  call dealloc_vertical_geom(gp%vgeom)
  deallocate(gp%vgeom)
  nullify(gp%vgeom)
end if


nullify(gp%lat )
nullify(gp%lon )
nullify(gp%timestep)
nullify(gp%orog )
nullify(gp%al   )
nullify(gp%tstar)
nullify(gp%t0   )
nullify(gp%ts   ) 
nullify(gp%ls   ) 
nullify(gp%ci   ) 
nullify(gp%ws   ) 
nullify(gp%sn   ) 
nullify(gp%u10n ) 
nullify(gp%v10n ) 
nullify(gp%u0   ) 
nullify(gp%v0   ) 
nullify(gp%chn  ) 
nullify(gp%tcls ) 
nullify(gp%hucls) 
nullify(gp%ucls) 
nullify(gp%vcls) 
nullify(gp%nucls) 
nullify(gp%nvcls) 
nullify(gp%z0   ) 
nullify(gp%wl   ) 
nullify(gp%sst  ) 
nullify(gp%veg  ) 
nullify(gp%cptgz) 
nullify(gp%wet  ) 
nullify(gp%qs   ) 
nullify(gp%cpts )  
nullify(gp%ri   )  
nullify(gp%bn   )  
nullify(gp%bm   )  
nullify(gp%bh   ) 
nullify(gp%upd  ) 
nullify(gp%arg  ) 
nullify(gp%hv   ) 
nullify(gp%z0h  ) 
nullify(gp%wsi  ) 
nullify(gp%rrr  ) 
nullify(gp%rsn  ) 
nullify(gp%es   )
nullify(gp%prec_accum)
nullify(gp%sdfor)
nullify(gp%t2m  )
nullify(gp%d2m  )
nullify(gp%q2m  )
nullify(gp%u10m )
nullify(gp%v10m )
nullify(gp%presh)
nullify(gp%presf)
nullify(gp%geoph)
nullify(gp%geopf)
nullify(gp%aphi )
nullify(gp%cp   )
nullify(gp%rf   )
nullify(gp%kap  )
nullify(gp%uf   )
nullify(gp%vf   )
nullify(gp%tf   )
nullify(gp%qf   )
nullify(gp%qn   )
nullify(gp%o3f  )
nullify(gp%clwf )
nullify(gp%ciwf )
nullify(gp%ccf  )
nullify(gp%hgt  )
nullify(gp%cswf )
nullify(gp%crwf )
nullify(gp%cgwf )
nullify(gp%clw_sum )
nullify(gp%ciw_sum )
nullify(gp%crw_sum )
nullify(gp%csw_sum )
nullify(gp%clwd )
nullify(gp%ciwd )
nullify(gp%ccd  )
nullify(gp%phys_type ) !KaLo
nullify(gp%rfl  )
nullify(gp%sfl  )
nullify(gp%crfl )
nullify(gp%csfl )
nullify(gp%pfrc )
nullify(gp%pxp  )
nullify(gp%pxpd )
nullify(gp%pxz  )
nullify(gp%pxzd )
nullify(gp%prec_accum)
nullify(gp%sdfor)
nullify(gp%t2m  )
nullify(gp%soty )
nullify(gp%lail )
nullify(gp%laih )
nullify(gp%sndt )
nullify(gp%sndp )
nullify(gp%cvl  )
nullify(gp%cvh  )
nullify(gp%tvl  )
nullify(gp%tvh  )
nullify(gp%tsoil)
nullify(gp%qsoil)
nullify(gp%tlk  )
! canari
nullify(gp%cori  )
nullify(gp%eh    )
nullify(gp%ez    )
nullify(gp%est   )
nullify(gp%esn   )
nullify(gp%et2   )
nullify(gp%eh2   )
nullify(gp%ev1   )
nullify(gp%gnordl)
nullify(gp%gnordm)

if (lhook) call dr_hook('gom_plus_destroy',1,zhook_handle)
end subroutine gom_plus_destroy

! --------------------------------------------------------------
! Dump GOM plus to a file: open and create file
! --------------------------------------------------------------
subroutine gom_plus_dump_open(gpd, cdtag, kproc, knsets, kflevg, ld_netcdf_format)

type(type_gom_plus_dump), intent(inout) :: gpd
character(len=*)        , intent(in)    :: cdtag  ! to distinguish e.g. TL and traj dumps
integer(kind=jpim)      , intent(in)    :: kproc  ! MPI processor ID
integer(kind=jpim)      , intent(in)    :: knsets ! Number of sets
integer(kind=jpim)      , intent(in)    :: kflevg ! Number of vertical levels
logical                 , intent(in)    :: ld_netcdf_format ! T=NetCDF ; F=F77 binary

integer(kind=jpib)  :: ipos_index
integer(kind=jpim) :: ioerr, idimids(2), idimids0(2), idimid(1)
integer(kind=jpim) :: id_pressure, id_pressure0, id_obs

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_dump_open',0,zhook_handle)

gpd%lnetcdf = ld_netcdf_format
gpd%counter = 0

if (.not.gpd%lnetcdf) then

  ! Binary F77 format dump = complete dump of everything for offline debugging

  write(gpd%clfilename,'("gom_plus_dump_",A,"_pe_",i4.4)') cdtag, kproc
  open(newunit=gpd%file_id,file=trim(gpd%clfilename),status='replace',action='write',&
    & form='unformatted',access='stream')

  write(gpd%file_id) knsets

  ! Create an empty index
  inquire(gpd%file_id, pos=ipos_index) ! not used yet - needs F2008
  gpd%nsets=knsets
  allocate(gpd%ipos_set(knsets))
  gpd%ipos_set = -1_JPIM
  write(gpd%file_id) gpd%ipos_set

else

  ! NetCDF dump format is partial, aimed at creating offline profile datasets
  ! at observation locations, using the IFS as an interpolation tool.

  write(gpd%clfilename,'("gom_plus_dump_",A,"_pe_",i4.4,".nc")') cdtag, kproc
  ioerr = nf_create(trim(gpd%clfilename),nf_clobber,gpd%file_id)
  ioerr = nf_def_dim(gpd%file_id, "pressure", kflevg, id_pressure)
  ioerr = nf_def_dim(gpd%file_id, "half_pressure", kflevg+1, id_pressure0)
  ioerr = nf_def_dim(gpd%file_id, "observation", nf_unlimited, id_obs)
  idimid   = (/ id_obs /)
  idimids  = (/ id_pressure,  id_obs /)
  idimids0 = (/ id_pressure0, id_obs /)

  allocate(gpd%dimid(26))

  ioerr = nf_def_var(gpd%file_id, "latitude",             nf_DOUBLE, 1, idimid,  gpd%dimid(1))
  ioerr = nf_def_var(gpd%file_id, "longitude",            nf_DOUBLE, 1, idimid,  gpd%dimid(2))
  ioerr = nf_def_var(gpd%file_id, "timestep",             nf_DOUBLE, 1, idimid,  gpd%dimid(3))
  ioerr = nf_def_var(gpd%file_id, "orography",            nf_DOUBLE, 1, idimid,  gpd%dimid(4))
  ioerr = nf_def_var(gpd%file_id, "tsfc",                 nf_DOUBLE, 1, idimid,  gpd%dimid(5))
  ioerr = nf_def_var(gpd%file_id, "lsm",                  nf_DOUBLE, 1, idimid,  gpd%dimid(6))
  ioerr = nf_def_var(gpd%file_id, "seaice",               nf_DOUBLE, 1, idimid,  gpd%dimid(7))
  ioerr = nf_def_var(gpd%file_id, "snow_water_equivalent",nf_DOUBLE, 1, idimid,  gpd%dimid(8))
  ioerr = nf_def_var(gpd%file_id, "t2m",                  nf_DOUBLE, 1, idimid,  gpd%dimid(9))
  ioerr = nf_def_var(gpd%file_id, "u10",                  nf_DOUBLE, 1, idimid,  gpd%dimid(10))
  ioerr = nf_def_var(gpd%file_id, "v10",                  nf_DOUBLE, 1, idimid,  gpd%dimid(11))

  ioerr = nf_def_var(gpd%file_id, "pressure_full",        nf_DOUBLE, 2, idimids,  gpd%dimid(12))
  ioerr = nf_def_var(gpd%file_id, "pressure_half",        nf_DOUBLE, 2, idimids0, gpd%dimid(13))
  ioerr = nf_def_var(gpd%file_id, "geopotential_full",    nf_DOUBLE, 2, idimids,  gpd%dimid(14))
  ioerr = nf_def_var(gpd%file_id, "geopotential_half",    nf_DOUBLE, 2, idimids0, gpd%dimid(15))

  ioerr = nf_def_var(gpd%file_id, "t",                    nf_DOUBLE, 2, idimids,  gpd%dimid(16))
  ioerr = nf_def_var(gpd%file_id, "q",                    nf_DOUBLE, 2, idimids,  gpd%dimid(17))
  ioerr = nf_def_var(gpd%file_id, "o3",                   nf_DOUBLE, 2, idimids,  gpd%dimid(18))
  ioerr = nf_def_var(gpd%file_id, "clw",                  nf_DOUBLE, 2, idimids,  gpd%dimid(19))
  ioerr = nf_def_var(gpd%file_id, "ciw",                  nf_DOUBLE, 2, idimids,  gpd%dimid(20))
  ioerr = nf_def_var(gpd%file_id, "cc",                   nf_DOUBLE, 2, idimids,  gpd%dimid(21))
  ioerr = nf_def_var(gpd%file_id, "rain_ls",              nf_DOUBLE, 2, idimids,  gpd%dimid(22))
  ioerr = nf_def_var(gpd%file_id, "snow_ls",              nf_DOUBLE, 2, idimids,  gpd%dimid(23))
  ioerr = nf_def_var(gpd%file_id, "precip_fraction_ls",   nf_DOUBLE, 2, idimids,  gpd%dimid(24))
  ioerr = nf_def_var(gpd%file_id, "rain_cv",              nf_DOUBLE, 2, idimids,  gpd%dimid(25))
  ioerr = nf_def_var(gpd%file_id, "snow_cv",              nf_DOUBLE, 2, idimids,  gpd%dimid(26))

  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(1), 'units', 3, 'rad')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(2), 'units', 3, 'rad')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(3), 'units', 1, ' ')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(4), 'units', 7, 'm^2 s-2')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(5), 'units', 1, 'K')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(6), 'units', 3, '0-1')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(7), 'units', 3, '0-1')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(8), 'units', 1, 'm')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(9), 'units', 1, 'K')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(10), 'units', 5, 'm s-1')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(11), 'units', 5, 'm s-1')

  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(12), 'units', 2, 'Pa')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(13), 'units', 2, 'Pa')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(14), 'units', 7, 'm^2 s-2')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(15), 'units', 7, 'm^2 s-2')

  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(16), 'units', 1, 'K')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(17), 'units', 7, 'kg kg-1')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(18), 'units', 7, 'kg kg-1')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(19), 'units', 7, 'kg kg-1')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(20), 'units', 7, 'kg kg-1')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(21), 'units', 3, '0-1')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(22), 'units', 10, 'kg m-2 s-1')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(23), 'units', 10, 'kg m-2 s-1')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(24), 'units', 3,  '0-1')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(25), 'units', 10, 'kg m-2 s-1')
  ioerr = nf_put_att_text(gpd%file_id, gpd%dimid(26), 'units', 10, 'kg m-2 s-1')

  ioerr = nf_enddef(gpd%file_id)

endif

if (lhook) call dr_hook('gom_plus_dump_open',1,zhook_handle)
end subroutine gom_plus_dump_open

! --------------------------------------------------------------
! Dump GOM plus to a file: used for debugging and test harnesses
! --------------------------------------------------------------
subroutine gom_plus_dump(gp, ktslot, gpd)

type(type_gom_plus),      intent(in)    :: gp
integer(kind=jpim),       intent(in)    :: ktslot
type(type_gom_plus_dump), intent(inout) :: gpd

integer(kind=jpim) :: ierr, istart(2), icount(2), icount0(2), jobs

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_dump',0,zhook_handle)

! Write one set to file
if (.not.gpd%lnetcdf) then

  ! F77 binary format

  ! Counter increases once per set
  gpd%counter = gpd%counter + 1
  inquire(gpd%file_id, pos=gpd%ipos_set(gpd%counter))
  write(gpd%file_id) gpd%counter ! Acts as a kind of record header / cross-check

  write(gpd%file_id) gp%ndlen, gp%nflevg, gp%nhoriz, gp%nppm, gp%nxyb, &
    & gp%ngems, associated(gp%cswf), associated(gp%clwd), &
    & gp%missing_value, gp%nlev_soil
  write(gpd%file_id) gp%store_sfc
  write(gpd%file_id) gp%store_ua0
  write(gpd%file_id) gp%store_ua1
  write(gpd%file_id) gp%store_soil
  write(gpd%file_id) gp%store_ppm
  write(gpd%file_id) gp%gems
  write(gpd%file_id) gp%gems_igrib
  write(gpd%file_id) gp%gems_type
  write(gpd%file_id) gp%pxyb

else

  ! NetCDF format

  istart(1) = 1
  icount(1) = gp%nflevg
  icount(2) = 1
  icount0(1) = gp%nflevg+1
  icount0(2) = 1

  do jobs = 1,gp%ndlen

    ! Counter increases once per obs
    gpd%counter = gpd%counter + 1
    istart(2) = gpd%counter

    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(1),  gpd%counter, 1, gp%lat(jobs,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(2),  gpd%counter, 1, gp%lon(jobs,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(3),  gpd%counter, 1, gp%timestep(jobs,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(4),  gpd%counter, 1, gp%orog(jobs,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(5),  gpd%counter, 1, gp%ts(jobs,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(6),  gpd%counter, 1, gp%ls(jobs,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(7),  gpd%counter, 1, gp%ci(jobs,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(8),  gpd%counter, 1, gp%sn(jobs,1))

    ! AJGDB these are as used in RTTOV/RTTOV_SCATT noting than model T2M and neutral winds are now available
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(9),  gpd%counter, 1, gp%tf(jobs,gp%nflevg,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(10), gpd%counter, 1, gp%uf(jobs,gp%nflevg,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(11), gpd%counter, 1, gp%vf(jobs,gp%nflevg,1))

    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(12), istart, icount,  gp%presf(jobs,1:gp%nflevg,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(13), istart, icount0, gp%presh(jobs,0:gp%nflevg,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(14), istart, icount,  gp%geopf(jobs,1:gp%nflevg,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(15), istart, icount0, gp%geoph(jobs,0:gp%nflevg,1))

    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(16), istart, icount,  gp%tf(jobs,1:gp%nflevg,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(17), istart, icount,  gp%qf(jobs,1:gp%nflevg,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(18), istart, icount,  gp%o3f(jobs,1:gp%nflevg,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(19), istart, icount,  gp%clwd(jobs,1:gp%nflevg,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(20), istart, icount,  gp%ciwd(jobs,1:gp%nflevg,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(21), istart, icount,  gp%ccd(jobs,1:gp%nflevg,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(22), istart, icount,  gp%rfl(jobs,1:gp%nflevg,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(23), istart, icount,  gp%sfl(jobs,1:gp%nflevg,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(24), istart, icount,  gp%pfrc(jobs,1:gp%nflevg,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(25), istart, icount,  gp%crfl(jobs,1:gp%nflevg,1))
    ierr = nf_put_vara_double(gpd%file_id, gpd%dimid(26), istart, icount,  gp%csfl(jobs,1:gp%nflevg,1))

  enddo

endif

if (lhook) call dr_hook('gom_plus_dump',1,zhook_handle)
end subroutine gom_plus_dump

! --------------------------------------------------------------
! Dump GOM plus to a file: close file
! --------------------------------------------------------------
subroutine gom_plus_dump_close(gpd)

type(type_gom_plus_dump), intent(inout) :: gpd

integer(kind=jpim) :: ioerr

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_dump_open',0,zhook_handle)

if (.not.gpd%lnetcdf) then

  ! Go back and fill the index with the starting point of each set record
  ! AJGDB Under Cray, this is the only working method I have found so far.
  ! I would like to be able to use "write(iu, pos=ipos_index)"
  close(gpd%file_id)
  open(newunit=gpd%file_id,file=gpd%clfilename,status='old',action='write',&
    & form='unformatted',access='stream')

  write(gpd%file_id) gpd%nsets
  write(gpd%file_id) gpd%ipos_set
  deallocate(gpd%ipos_set)

  close(gpd%file_id)

else

  deallocate(gpd%dimid)
  ioerr = nf_close(gpd%file_id)

endif

if (lhook) call dr_hook('gom_plus_dump_close',1,zhook_handle)
end subroutine gom_plus_dump_close

! -------------------------------------------------------------------
! Read GOM plus from dump file: used for debugging and test harnesses
! Note: only works with binary dump format
! -------------------------------------------------------------------
subroutine gom_plus_read_dump(gp, ku)

type(type_gom_plus), intent(inout) :: gp
integer(kind=jpim), intent(in)     :: ku

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_read_dump',0,zhook_handle)

read(ku) gp%ndlen, gp%nflevg, gp%nhoriz, gp%nppm, gp%nxyb, gp%ngems, &
  & gp%lsrg, gp%lhydro_sum, gp%lphys, gp%lcanari, gp%missing_value, gp%nlev_soil

call gom_plus_alloc(gp, gp%ndlen, gp%nflevg, gp%nhoriz, gp%nppm, gp%nxyb, gp%lsrg, gp%lhydro_sum, &
  & gp%lphys, gp%lcanari, gp%ngems, gp%nlev_soil)

read(ku) gp%store_sfc
read(ku) gp%store_ua0
read(ku) gp%store_ua1
read(ku) gp%store_soil
read(ku) gp%store_ppm
read(ku) gp%gems
read(ku) gp%gems_igrib
read(ku) gp%gems_type
read(ku) gp%pxyb

if (lhook) call dr_hook('gom_plus_read_dump',1,zhook_handle)
end subroutine gom_plus_read_dump

! --------------------------------------------------------------
! Zero all fields in a gom_plus
! --------------------------------------------------------------
subroutine gom_plus_zero(gp)

type(type_gom_plus), intent(inout) :: gp

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_zero',0,zhook_handle)

gp%store_sfc=0.0_jprb
gp%store_ua0=0.0_jprb
gp%store_ua1=0.0_jprb
gp%store_soil=0.0_jprb
gp%store_ppm=0.0_jprb
gp%gems=0.0_jprb
gp%pxyb=0.0_jprb

if (lhook) call dr_hook('gom_plus_zero',1,zhook_handle)
end subroutine gom_plus_zero

! --------------------------------------------------------------
! Copy a gom_plus (pointers would be incorrect with Fortran
! default assignment)
!
! khoriz_reduce is optional and just for the slant path to stop
! that duplicating almost all this code and generating bugs (3 
! so far)
!
! --------------------------------------------------------------
subroutine gom_plus_copy(gp,gp2,khoriz_reduce)

type(type_gom_plus), intent(inout) :: gp
type(type_gom_plus), intent(in)    :: gp2
integer(kind=jpim), optional, intent(in) :: khoriz_reduce

integer(kind=jpim) :: ihoriz
real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_copy',0,zhook_handle)

ihoriz = gp2%nhoriz
if(present(khoriz_reduce)) ihoriz=khoriz_reduce

call gom_plus_alloc(gp, gp2%ndlen, gp2%nflevg, ihoriz, gp2%nppm, &
  & gp2%nxyb, associated(gp2%cswf), associated(gp2%clw_sum), &
  & associated(gp2%clwd), associated(gp2%eh), gp2%ngems, gp2%nlev_soil)

! dimensions and other settings
gp%ndlen        =gp2%ndlen
gp%nflevg       =gp2%nflevg
gp%nhoriz       =ihoriz
gp%nppm         =gp2%nppm
gp%nxyb         =gp2%nxyb
gp%ngems        =gp2%ngems
gp%lsrg         =gp2%lsrg
gp%lhydro_sum   =gp2%lhydro_sum
gp%lphys        =gp2%lphys
gp%lcanari      =gp2%lcanari
gp%missing_value=gp2%missing_value
gp%nlev_soil=gp2%nlev_soil

if(.not.present(khoriz_reduce)) then

  ! Main storage arrays for atmospheric data
  gp%store_sfc =gp2%store_sfc
  gp%store_ua0 =gp2%store_ua0
  gp%store_ua1 =gp2%store_ua1
  gp%store_soil=gp2%store_soil
  gp%store_ppm =gp2%store_ppm
  gp%gems      =gp2%gems
  gp%pxyb      =gp2%pxyb

endif

! Other info that's passed around with the gom_plus
gp%gems_igrib=gp2%gems_igrib
gp%gems_type =gp2%gems_type
call alloc_copy_vertical_geom( gp%vgeom, gp2%vgeom%yrvab, gp2%vgeom%yrveta, gp2%vgeom%yrvfe )

if (lhook) call dr_hook('gom_plus_copy',1,zhook_handle)
end subroutine gom_plus_copy

! --------------------------------------------------------------
! Dot product
! --------------------------------------------------------------
subroutine gom_plus_dotprod(gp,gp2,prod)

type(type_gom_plus), intent(in)  :: gp
type(type_gom_plus), intent(in)  :: gp2
real(kind=jprb),     intent(out) :: prod

real(kind=jprb) :: zhook_handle

if (lhook) call dr_hook('gom_plus_dotprod',0,zhook_handle)

if(gp%ndlen     /= gp2%ndlen)     call abor1('GOM_PLUS_DOTPROD: Incompatible gom_plus - ndlen') 
if(gp%nflevg    /= gp2%nflevg)    call abor1('GOM_PLUS_DOTPROD: Incompatible gom_plus - nflevg') 
if(gp%nhoriz    /= gp2%nhoriz)    call abor1('GOM_PLUS_DOTPROD: Incompatible gom_plus - nhoriz') 
if(gp%nppm      /= gp2%nppm)      call abor1('GOM_PLUS_DOTPROD: Incompatible gom_plus - nppm') 
if(gp%nxyb      /= gp2%nxyb)      call abor1('GOM_PLUS_DOTPROD: Incompatible gom_plus - nxyb') 
if(gp%ngems     /= gp2%ngems)     call abor1('GOM_PLUS_DOTPROD: Incompatible gom_plus - ngems') 
if(gp%nlev_soil /= gp2%nlev_soil) call abor1('GOM_PLUS_DOTPROD: Incompatible gom_plus - nlev_soil')

prod =        sum(gp%store_sfc*gp2%store_sfc)
prod = prod + sum(gp%store_ua0*gp2%store_ua0)
prod = prod + sum(gp%store_ua1*gp2%store_ua1)
prod = prod + sum(gp%store_soil*gp2%store_soil)
prod = prod + sum(gp%store_ppm*gp2%store_ppm)
prod = prod + sum(gp%gems*gp2%gems)
prod = prod + sum(gp%pxyb*gp2%pxyb)

if (lhook) call dr_hook('gom_plus_dotprod',1,zhook_handle)
end subroutine gom_plus_dotprod


! Subroutine to convert a 2d gom_plus to a 1d slant-path profile.
!     ------------------------------------------------------------------
subroutine gom_plus_make_slant_path(gpfull, gpslant, pzenith)

type(type_gom_plus),intent(in)    :: gpfull   ! Input 2d-gom_plus
type(type_gom_plus),intent(inout) :: gpslant  ! Output 1d slant-path profile
real(kind=jprb),intent(in)   :: pzenith(gpfull%ndlen) ! Zenith angle in [deg]

real(kind=jprb)    :: zdx ! Spacing between 2d-gom profiles [m]
integer(kind=jpim) :: iobs, ilev,invars,jvar
integer(kind=jpim) :: ix1s(gpfull%ndlen,gpfull%nflevg), ix2s(gpfull%ndlen,gpfull%nflevg)
integer(kind=jpim) :: ix1, ix2
real(kind=jprb)    :: zw1s(gpfull%ndlen,gpfull%nflevg), zw2s(gpfull%ndlen,gpfull%nflevg)
real(kind=jprb)    :: zw1, zw2
real(kind=jprb)    :: zxloc, z_hgt(gpfull%ndlen)

real(kind=jprb) :: zhook_handle

!=============================================================

if (lhook) call dr_hook('gom_plus_make_slant_path',0,zhook_handle)

! Allocate the output gom_plus structure
zdx = rdxmax_slant / real(gpfull%nhoriz-1)
if (gpfull%nhoriz > 0) then
  call gom_plus_copy(gpslant,gpfull,khoriz_reduce=1)
else
  call abor1('make_slant_path: input structure has horizontal dimension 0 - this should not happen!')
endif

! Loop through observations and derive the slant path
do iobs = 1, gpfull%ndlen
  if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then 
    ! Calculate  displacement in plane-parallel 2d atmosphere and interpolate.
    ! Limit z_hgt to be above zero, so that zxloc stays within the covered plane.
    z_hgt(iobs) = max(0.0_jprb,gpfull%hgt(iobs,gpfull%nflevg,1))
  endif
enddo


do ilev = gpfull%nflevg, 1, -1          ! Loop from bottom to top.
  do iobs = 1, gpfull%ndlen
    if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then 

      zxloc =  z_hgt(iobs) * tan(pzenith(iobs) * rpi / 180.0_jprb)

      ix1s(iobs,ilev) = min(1 + int(zxloc/zdx), gpfull%nhoriz - 1)
      ix2s(iobs,ilev) = ix1s(iobs,ilev) + 1

      zw2 = min(1.0_jprb,(zxloc - (zdx * (ix1s(iobs,ilev)-1)))/zdx)
      zw1 = 1.0_jprb - zw2

      zw2s(iobs,ilev) = zw2
      zw1s(iobs,ilev) = zw1
      

      ! Interpolate height of next level above for (approximate) calculation of displacement at next level
      if (ilev > 1) then
        z_hgt(iobs) = zw1 * gpfull%hgt(iobs,ilev-1,ix1s(iobs,ilev)) + zw2 * gpfull%hgt(iobs,ilev-1,ix2s(iobs,ilev))
      endif

    endif
  enddo
enddo

      ! Linearly interpolate all variables.
      ! Note we might want to distinguish between different interpolation methods
      ! here in the future, e.g., nearest-point for cloud/rain fields, but that's
      ! not done yet. 

invars = size(gpslant%store_sfc(1,1,:))
do jvar=1,invars
  do iobs = 1, gpfull%ndlen
    if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then 
      
    ! Neglecting surface displacement here - maybe deal with it later
      gpslant%store_sfc(iobs,ih,jvar) = gpfull%store_sfc(iobs,1,jvar)
    else
      gpslant%store_sfc(iobs,ih,jvar)   = gpfull%store_sfc(iobs,1,jvar)
    endif
  enddo
enddo

invars = size(gpslant%store_soil(1,1,1,:))
do jvar=1,invars
  do ilev = 1, gpfull%nlev_soil
    do iobs = 1, gpfull%ndlen
      if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
      ! Neglecting soil displacement here - maybe deal with it later
        gpslant%store_soil(iobs,ilev,ih,jvar) = gpfull%store_soil(iobs,ilev,1,jvar)
      else
        gpslant%store_soil(iobs,ilev,ih,jvar) = gpfull%store_soil(iobs,ilev,1,jvar)
      endif
    enddo
  enddo
enddo

invars = size(gpslant%store_ua1(1,1,1,:))
do jvar=1,invars
  do ilev = 1, gpfull%nflevg
    do iobs = 1, gpfull%ndlen
      if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
        ix1 = ix1s(iobs,ilev)
        ix2 = ix2s(iobs,ilev)

        zw2 = zw2s(iobs,ilev)
        zw1 = zw1s(iobs,ilev)
        gpslant%store_ua1(iobs,ilev,ih,jvar) = zw1 * gpfull%store_ua1(iobs,ilev,ix1,jvar) + zw2 * gpfull%store_ua1(iobs,ilev,ix2,jvar)
      else
        ! Just use first profile if zenith angle is missing or really nadir
        gpslant%store_ua1(iobs,ilev,ih,jvar) = gpfull%store_ua1(iobs,ilev,1,jvar)
      endif
    enddo
  enddo
enddo

invars = size(gpslant%pxp(1,1,:,1))
do jvar=1,invars
  do ilev = 1, gpfull%nflevg
    do iobs = 1, gpfull%ndlen
      if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
        ix1 = ix1s(iobs,ilev)
        ix2 = ix2s(iobs,ilev)

        zw2 = zw2s(iobs,ilev)
        zw1 = zw1s(iobs,ilev)
        gpslant%pxp(iobs,ilev,jvar,ih)       = zw1 * gpfull%pxp(iobs,ilev,jvar,ix1)       + zw2 * gpfull%pxp(iobs,ilev,jvar,ix2)
        gpslant%pxz(iobs,ilev,jvar,ih)       = zw1 * gpfull%pxz(iobs,ilev,jvar,ix1)       + zw2 * gpfull%pxz(iobs,ilev,jvar,ix2)
        gpslant%pxpd(iobs,ilev-1,jvar,ih)    = zw1 * gpfull%pxpd(iobs,ilev-1,jvar,ix1)    + zw2 * gpfull%pxpd(iobs,ilev-1,jvar,ix2)
        gpslant%pxzd(iobs,ilev-1,jvar,ih)    = zw1 * gpfull%pxzd(iobs,ilev-1,jvar,ix1)    + zw2 * gpfull%pxzd(iobs,ilev-1,jvar,ix2)
      else
        ! Just use first profile if zenith angle is missing or really nadir
        gpslant%pxp(iobs,ilev,jvar,ih)     = gpfull%pxp(iobs,ilev,jvar,1)    
        gpslant%pxz(iobs,ilev,jvar,ih)     = gpfull%pxz(iobs,ilev,jvar,1)   
        gpslant%pxpd(iobs,ilev,jvar,ih)    = gpfull%pxpd(iobs,ilev,jvar,1) 
        gpslant%pxzd(iobs,ilev,jvar,ih)    = gpfull%pxzd(iobs,ilev,jvar,1)
        if (ilev == 1) then
          gpslant%pxp(iobs,0,jvar,ih)     = gpfull%pxp(iobs,0,jvar,1)    
          gpslant%pxz(iobs,0,jvar,ih)     = gpfull%pxz(iobs,0,jvar,1)   
          gpslant%pxpd(iobs,0,jvar,ih)    = gpfull%pxpd(iobs,0,jvar,1) 
          gpslant%pxzd(iobs,0,jvar,ih)    = gpfull%pxzd(iobs,0,jvar,1)
        endif
      endif
    enddo
  enddo
enddo

invars = size(gpslant%store_ua0(1,1,1,:))
do jvar=1,invars
  do ilev = 1, gpfull%nflevg
    do iobs = 1, gpfull%ndlen
      if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
        ix1 = ix1s(iobs,ilev)
        ix2 = ix2s(iobs,ilev)

        zw2 = zw2s(iobs,ilev)
        zw1 = zw1s(iobs,ilev)
        gpslant%store_ua0(iobs,ilev,ih,jvar) = zw1 * gpfull%store_ua0(iobs,ilev,ix1,jvar) + zw2 * gpfull%store_ua0(iobs,ilev,ix2,jvar)
        if (ilev == 1) then
          gpslant%store_ua0(iobs,0,ih,jvar) = zw1 * gpfull%store_ua0(iobs,0,ix1,jvar) + zw2 * gpfull%store_ua0(iobs,0,ix2,jvar)
        endif
      else
        ! Just use first profile if zenith angle is missing or really nadir
        gpslant%store_ua0(iobs,ilev,ih,jvar) = gpfull%store_ua0(iobs,ilev,1,jvar)
        if (ilev == 1) then
          gpslant%store_ua0(iobs,0,ih,jvar) = gpfull%store_ua0(iobs,0,1,jvar)
        endif
      endif
    enddo
  enddo
enddo

invars = size(gpslant%pxyb(1,1,:,1))
do jvar=1,invars
  do ilev = 1, gpfull%nflevg
    do iobs = 1, gpfull%ndlen
      if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
        ix1 = ix1s(iobs,ilev)
        ix2 = ix2s(iobs,ilev)

        zw2 = zw2s(iobs,ilev)
        zw1 = zw1s(iobs,ilev)
        gpslant%pxyb(iobs,ilev,jvar,ih)      = zw1 * gpfull%pxyb(iobs,ilev,jvar,ix1)      + zw2 * gpfull%pxyb(iobs,ilev,jvar,ix2)
      else
        ! Just use first profile if zenith angle is missing or really nadir
        gpslant%pxyb(iobs,ilev,jvar,ih)      = gpfull%pxyb(iobs,ilev,jvar,1)
      endif
    enddo
  enddo
enddo

if (gpfull%ngems > 0) then
  invars = size(gpslant%gems(1,1,1,:))
  do jvar=1,invars
    do ilev = 1, gpfull%nflevg
      do iobs = 1, gpfull%ndlen
        if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
          ix1 = ix1s(iobs,ilev)
          ix2 = ix2s(iobs,ilev)
          
          zw2 = zw2s(iobs,ilev)
          zw1 = zw1s(iobs,ilev)
          gpslant%gems(iobs,ilev,ih,jvar)    = zw1 * gpfull%gems(iobs,ilev,ix1,jvar)      + zw2 * gpfull%gems(iobs,ilev,ix2,jvar)
        else
          ! Just use first profile if zenith angle is missing or really nadir
          gpslant%gems(iobs,ilev,ih,jvar) = gpfull%gems(iobs,ilev,1,jvar) 
        endif
      enddo
    enddo
  enddo
endif

if (lhook) call dr_hook('gom_plus_make_slant_path',1,ZHOOK_HANDLE)

end subroutine gom_plus_make_slant_path


! TL of Subroutine to convert a 2d gom_plus to a 1d slant-path profile.
!     ------------------------------------------------------------------
subroutine gom_plus_make_slant_path_tl(gpfull, gpslant, gpfull5, pzenith)

type(type_gom_plus),intent(in)    :: gpfull   ! Input 2d-gom_plus
type(type_gom_plus),intent(inout) :: gpslant  ! Output 1d slant-path profile
type(type_gom_plus),intent(in)    :: gpfull5   ! Input 2d-gom_plus
real(kind=jprb),intent(in)   :: pzenith(gpfull%ndlen) ! Zenith angle in [deg]

real(kind=jprb)    :: zdx ! Spacing between 2d-gom profiles [m]
integer(kind=jpim) :: iobs, ilev
integer(kind=jpim) :: ix1s(gpfull%ndlen,gpfull5%nflevg), ix2s(gpfull%ndlen,gpfull5%nflevg)
integer(kind=jpim) :: ix1, ix2,invars,jvar
real(kind=jprb)    :: zw1(gpfull%ndlen,gpfull5%nflevg), zw2(gpfull%ndlen,gpfull5%nflevg)
real(kind=jprb)    :: zw1_5(gpfull%ndlen,gpfull5%nflevg), zw2_5(gpfull%ndlen,gpfull5%nflevg)
real(kind=jprb)    :: zxloc, z_hgt
real(kind=jprb)    :: zxloc5, z_hgt5

real(kind=jprb) :: zhook_handle

!=============================================================

if (lhook) call dr_hook('gom_plus_make_slant_path_tl',0,zhook_handle)

! Allocate the output gom_plus structure
zdx = rdxmax_slant / real(gpfull%nhoriz-1)
if (gpfull5%nhoriz > 0 .and. gpfull%nhoriz > 0) then
  call gom_plus_copy(gpslant,gpfull,khoriz_reduce=1)
else
  call abor1('make_slant_path_tl: input structure has horizontal dimension 0 - this should not happen!')
endif
if (gpfull5%nhoriz /= gpfull%nhoriz .or. gpfull5%nflevg /= gpfull%nflevg )then
  call abor1('make_slant_path_tl: input TL and direct GOM-plus inconsistent!')
endif

! Loop through observations and derive the slant path
do iobs = 1, gpfull5%ndlen
  if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then 

    ! Neglecting surface displacement here - maybe deal with it later
    gpslant%store_sfc(iobs,ih,:)  = gpfull%store_sfc(iobs,1,:)

    ! Neglecting soil displacement here - maybe deal with it later
    gpslant%store_soil(iobs,:,ih,:)  = gpfull%store_soil(iobs,:,1,:)

    ! Calculate  displacement in plane-parallel 2d atmosphere and interpolate.
    if(gpfull5%hgt(iobs,gpfull5%nflevg,1) >= 0.0_jprb)then
      z_hgt5 = gpfull5%hgt(iobs,gpfull5%nflevg,1)
      z_hgt  = gpfull%hgt(iobs,gpfull5%nflevg,1)
    else
      z_hgt5 = 0.0_jprb
      z_hgt  = 0.0_jprb
    endif

    do ilev = gpfull5%nflevg, 1, -1          ! Loop from bottom to top.

      zxloc5 =  z_hgt5 * tan(pzenith(iobs) * rpi / 180.0_jprb)
      zxloc  =  z_hgt  * tan(pzenith(iobs) * rpi / 180.0_jprb)

      ix1s(iobs,ilev) = min(1 + int(zxloc5/zdx), gpfull5%nhoriz - 1)
      ix2s(iobs,ilev) = ix1s(iobs,ilev) + 1

      zw2_5(iobs,ilev) = (zxloc5-(zdx*(ix1s(iobs,ilev)-1)))/zdx
      zw2(iobs,ilev)   = zxloc/zdx
      if(zw2_5(iobs,ilev) > 1.0_jprb)then
        zw2_5(iobs,ilev) = 1.0_jprb
        zw2(iobs,ilev)   = 0.0_jprb
      endif

      zw1_5(iobs,ilev) = 1.0_jprb - zw2_5(iobs,ilev)
      zw1(iobs,ilev)   = - zw2(iobs,ilev)

      ! Interpolate height of next level above for (approximate) calculation of displacement at next level
      if (ilev > 1) then
        z_hgt5 = zw1_5(iobs,ilev) * gpfull5%hgt(iobs,ilev-1,ix1s(iobs,ilev)) + zw2_5(iobs,ilev) * gpfull5%hgt(iobs,ilev-1,ix2s(iobs,ilev))
        z_hgt  = zw1_5(iobs,ilev) * gpfull%hgt(iobs,ilev-1,ix1s(iobs,ilev))  + zw2_5(iobs,ilev) * gpfull%hgt(iobs,ilev-1,ix2s(iobs,ilev))&
         &     + zw1(iobs,ilev)   * gpfull5%hgt(iobs,ilev-1,ix1s(iobs,ilev)) + zw2(iobs,ilev)   * gpfull5%hgt(iobs,ilev-1,ix2s(iobs,ilev))
      endif

    enddo
  endif
enddo

invars = size(gpslant%store_ua1(1,1,1,:))
do jvar=1,invars
  do ilev = 1, gpfull5%nflevg
    do iobs = 1, gpfull5%ndlen
      if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
        ix1 = ix1s(iobs,ilev)
        ix2 = ix2s(iobs,ilev)

        gpslant%store_ua1(iobs,ilev,ih,jvar) = zw1_5(iobs,ilev) * gpfull%store_ua1(iobs,ilev,ix1,jvar)&
         &                             + zw2_5(iobs,ilev) * gpfull%store_ua1(iobs,ilev,ix2,jvar)&
         &                             + zw1(iobs,ilev) * gpfull5%store_ua1(iobs,ilev,ix1,jvar) &
         &                             + zw2(iobs,ilev) * gpfull5%store_ua1(iobs,ilev,ix2,jvar)

      else
        ! Just use first profile if zenith angle is missing or really nadir
        gpslant%store_ua1(iobs,ilev,ih,jvar) = gpfull%store_ua1(iobs,ilev,1,jvar)
      endif
    enddo
  enddo
enddo

invars = size(gpslant%pxp(1,1,:,1))
do jvar=1,invars
  do ilev = 1, gpfull5%nflevg
    do iobs = 1, gpfull5%ndlen
      if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
        ix1 = ix1s(iobs,ilev)
        ix2 = ix2s(iobs,ilev)

        gpslant%pxp(iobs,ilev,jvar,ih)       = zw1_5(iobs,ilev) * gpfull%pxp(iobs,ilev,jvar,ix1)      &
         &                             + zw2_5(iobs,ilev) * gpfull%pxp(iobs,ilev,jvar,ix2)      &
         &                             + zw1(iobs,ilev) * gpfull5%pxp(iobs,ilev,jvar,ix1)       &
         &                             + zw2(iobs,ilev) * gpfull5%pxp(iobs,ilev,jvar,ix2)
        gpslant%pxz(iobs,ilev,jvar,ih)       = zw1_5(iobs,ilev) * gpfull%pxz(iobs,ilev,jvar,ix1)      &
         &                             + zw2_5(iobs,ilev) * gpfull%pxz(iobs,ilev,jvar,ix2)      &
         &                             + zw1(iobs,ilev) * gpfull5%pxz(iobs,ilev,jvar,ix1)       &
         &                             + zw2(iobs,ilev) * gpfull5%pxz(iobs,ilev,jvar,ix2)
        gpslant%pxpd(iobs,ilev-1,jvar,ih)     = zw1_5(iobs,ilev) * gpfull%pxpd(iobs,ilev-1,jvar,ix1)   &
         &                             + zw2_5(iobs,ilev) * gpfull%pxpd(iobs,ilev-1,jvar,ix2)   &
         &                             + zw1(iobs,ilev) * gpfull5%pxpd(iobs,ilev-1,jvar,ix1)    &
         &                             + zw2(iobs,ilev) * gpfull5%pxpd(iobs,ilev-1,jvar,ix2)
        gpslant%pxzd(iobs,ilev-1,jvar,ih)    = zw1_5(iobs,ilev) * gpfull%pxzd(iobs,ilev-1,jvar,ix1)   &
         &                             + zw2_5(iobs,ilev) * gpfull%pxzd(iobs,ilev-1,jvar,ix2)   &
         &                             + zw1(iobs,ilev) * gpfull5%pxzd(iobs,ilev-1,jvar,ix1)    &
         &                             + zw2(iobs,ilev) * gpfull5%pxzd(iobs,ilev-1,jvar,ix2)
       else
        ! Just use first profile if zenith angle is missing or really nadir
         gpslant%pxp(iobs,ilev,jvar,ih)     = gpfull%pxp(iobs,ilev,jvar,1)    
         gpslant%pxz(iobs,ilev,jvar,ih)     = gpfull%pxz(iobs,ilev,jvar,1)   
         gpslant%pxpd(iobs,ilev,jvar,ih)    = gpfull%pxpd(iobs,ilev,jvar,1) 
         gpslant%pxzd(iobs,ilev,jvar,ih)    = gpfull%pxzd(iobs,ilev,jvar,1)
         if(ilev==1)then
           gpslant%pxp(iobs,0,jvar,ih)     = gpfull%pxp(iobs,0,jvar,1)    
           gpslant%pxz(iobs,0,jvar,ih)     = gpfull%pxz(iobs,0,jvar,1)   
           gpslant%pxpd(iobs,0,jvar,ih)    = gpfull%pxpd(iobs,0,jvar,1) 
           gpslant%pxzd(iobs,0,jvar,ih)    = gpfull%pxzd(iobs,0,jvar,1)
         endif
      endif
    enddo
  enddo
enddo

invars = size(gpslant%store_ua0(1,1,1,:))
do jvar=1,invars
  do ilev = 1, gpfull5%nflevg
    do iobs = 1, gpfull5%ndlen
      if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
        ix1 = ix1s(iobs,ilev)
        ix2 = ix2s(iobs,ilev)

        gpslant%store_ua0(iobs,ilev,ih,jvar) = zw1_5(iobs,ilev) * gpfull%store_ua0(iobs,ilev,ix1,jvar)&
         &                             + zw2_5(iobs,ilev) * gpfull%store_ua0(iobs,ilev,ix2,jvar)&
         &                             + zw1(iobs,ilev) * gpfull5%store_ua0(iobs,ilev,ix1,jvar) &
         &                             + zw2(iobs,ilev) * gpfull5%store_ua0(iobs,ilev,ix2,jvar)

        if (ilev == 1) then
          gpslant%store_ua0(iobs,0,ih,jvar) = zw1_5(iobs,ilev) * gpfull%store_ua0(iobs,0,ix1,jvar)&
           &                             + zw2_5(iobs,ilev) * gpfull%store_ua0(iobs,0,ix2,jvar)&
           &                             + zw1(iobs,ilev) * gpfull5%store_ua0(iobs,0,ix1,jvar) &
           &                             + zw2(iobs,ilev) * gpfull5%store_ua0(iobs,0,ix2,jvar)
        endif
      else
        ! Just use first profile if zenith angle is missing or really nadir
        gpslant%store_ua0(iobs,ilev,ih,jvar) = gpfull%store_ua0(iobs,ilev,1,jvar)
        if (ilev == 1) then
          gpslant%store_ua0(iobs,0,ih,jvar) = gpfull%store_ua0(iobs,0,1,jvar)
        endif
      endif
    enddo
  enddo
enddo

invars = size(gpslant%pxyb(1,1,:,1))
do jvar=1,invars
  do ilev = 1, gpfull5%nflevg
    do iobs = 1, gpfull5%ndlen
      if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
        ix1 = ix1s(iobs,ilev)
        ix2 = ix2s(iobs,ilev)

        gpslant%pxyb(iobs,ilev,jvar,ih)      = zw1_5(iobs,ilev) * gpfull%pxyb(iobs,ilev,jvar,ix1)     &
         &                             + zw2_5(iobs,ilev) * gpfull%pxyb(iobs,ilev,jvar,ix2)     &
         &                             + zw1(iobs,ilev) * gpfull5%pxyb(iobs,ilev,jvar,ix1)      &
         &                             + zw2(iobs,ilev) * gpfull5%pxyb(iobs,ilev,jvar,ix2)
      else
        ! Just use first profile if zenith angle is missing or really nadir
        gpslant%pxyb(iobs,ilev,jvar,ih)      = gpfull%pxyb(iobs,ilev,jvar,1)
      endif
    enddo
  enddo
enddo

if (gpfull%ngems > 0) then
  invars = size(gpslant%gems(1,1,1,:))
  do jvar=1,invars
    do ilev = 1, gpfull5%nflevg
      do iobs = 1, gpfull5%ndlen
        if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
          ix1 = ix1s(iobs,ilev)
          ix2 = ix2s(iobs,ilev)
          
          gpslant%gems(iobs,ilev,ih,jvar)    = zw1_5(iobs,ilev) * gpfull%gems(iobs,ilev,ix1,jvar)     &
           &                           + zw2_5(iobs,ilev) * gpfull%gems(iobs,ilev,ix2,jvar)     &
           &                           + zw1(iobs,ilev) * gpfull5%gems(iobs,ilev,ix1,jvar)      &
           &                           + zw2(iobs,ilev) * gpfull5%gems(iobs,ilev,ix2,jvar)
        else
          ! Just use first profile if zenith angle is missing or really nadir
          gpslant%gems(iobs,ilev,ih,jvar) = gpfull%gems(iobs,ilev,1,jvar) 
        endif
      enddo
    enddo
  enddo
endif



if (lhook) call dr_hook('gom_plus_make_slant_path_tl',1,zhook_handle)

end subroutine gom_plus_make_slant_path_tl

! AD of Subroutine to convert a 2d gom_plus to a 1d slant-path profile.
!     ------------------------------------------------------------------
subroutine gom_plus_make_slant_path_ad(gpfull, gpslant, gpfull5, pzenith, ld_alloc)

type(type_gom_plus),intent(inout) :: gpfull   ! Gradient of 2d-gom_plus
type(type_gom_plus),intent(inout) :: gpslant  ! Input gradients of 1d slant-path profile
type(type_gom_plus),intent(in)    :: gpfull5  ! Input 2d-gom_plus (direct)
real(kind=jprb),intent(in)        :: pzenith(gpfull5%ndlen) ! Zenith angle in [deg]
logical, intent(in)               :: ld_alloc

real(kind=jprb)    :: zdx ! Spacing between 2d-gom profiles [m]
integer(kind=jpim) :: iobs, ilev, invars, jvar
integer(kind=jpim) :: ix1, ix2
integer(kind=jpim) :: ix1s(gpfull%ndlen,gpfull5%nflevg), ix2s(gpfull%ndlen,gpfull5%nflevg)
real(kind=jprb)    :: zw1(gpfull%ndlen,gpfull5%nflevg), zw2(gpfull%ndlen,gpfull5%nflevg)
real(kind=jprb)    :: zw1_5(gpfull%ndlen,gpfull5%nflevg), zw2_5(gpfull%ndlen,gpfull5%nflevg)
real(kind=jprb)    :: zw2_5a(gpfull%ndlen,gpfull5%nflevg)
real(kind=jprb)    :: zxloc, z_hgt
real(kind=jprb)    :: zxloc5, z_hgt5

real(kind=jprb) :: zhook_handle

!=============================================================

if (lhook) call dr_hook('gom_plus_make_slant_path_ad',0,zhook_handle)

! allocate the output gom_plus structure
zdx = rdxmax_slant / real(gpfull%nhoriz-1)
if (gpfull5%nhoriz <= 0 .or. gpfull%nhoriz <= 0) then
  call abor1('make_slant_path_ad: input structure has horizontal dimension 0 - this should not happen!')
endif
if (gpfull5%nhoriz /= gpfull%nhoriz .or. gpfull5%nflevg /= gpfull%nflevg )then
  call abor1('make_slant_path_ad: input AD and direct GOM-plus inconsistent!')
endif

if(ld_alloc)then
  call gom_plus_copy(gpslant,gpfull,khoriz_reduce=1)
  call gom_plus_zero(gpslant)
else

  ! Loop through observations and derive the slant path

  do iobs = 1, gpfull%ndlen
    if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then 
 
      ! Calculate  displacement in plane-parallel 2d atmosphere and interpolate.
      z_hgt5 = max(0.0_jprb,gpfull5%hgt(iobs,gpfull5%nflevg,1))

      do ilev = gpfull5%nflevg, 1, -1          ! Loop from bottom to top.

        zxloc5 =  z_hgt5 * tan(pzenith(iobs) * rpi / 180.0_jprb)

        ix1s(iobs,ilev) = min(1 + int(zxloc5/zdx), gpfull5%nhoriz - 1)
        ix2s(iobs,ilev) = ix1s(iobs,ilev) + 1

        zw2_5a(iobs,ilev) = (zxloc5 - (zdx * (ix1s(iobs,ilev)-1)))/zdx
        zw2_5(iobs,ilev) = min(1.0_jprb,zw2_5a(iobs,ilev))
        zw1_5(iobs,ilev) = 1.0_jprb - zw2_5(iobs,ilev)

        ! Interpolate height of next level above for (approximate) calculation of displacement at next level
        if (ilev > 1) then
          z_hgt5 = zw1_5(iobs,ilev) * gpfull5%hgt(iobs,ilev-1,ix1s(iobs,ilev)) + zw2_5(iobs,ilev) * gpfull5%hgt(iobs,ilev-1,ix2s(iobs,ilev))
        endif

      enddo
    endif
  enddo

  zw1(:,:) = 0.0_jprb
  zw2(:,:) = 0.0_jprb

  if (gpfull%ngems > 0) then
    invars = size(gpslant%gems(1,1,1,:))
    do jvar=1,invars
      do ilev = 1, gpfull5%nflevg
        do iobs = 1, gpfull5%ndlen
          if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
            ix1 = ix1s(iobs,ilev)
            ix2 = ix2s(iobs,ilev)
            gpfull%gems(iobs,ilev,ix1,jvar) = gpfull%gems(iobs,ilev,ix1,jvar) + zw1_5(iobs,ilev) *gpslant%gems(iobs,ilev,ih,jvar)
            gpfull%gems(iobs,ilev,ix2,jvar) = gpfull%gems(iobs,ilev,ix2,jvar) + zw2_5(iobs,ilev) *gpslant%gems(iobs,ilev,ih,jvar)
            zw1(iobs,ilev) = zw1(iobs,ilev) + gpfull5%gems(iobs,ilev,ix1,jvar) * gpslant%gems(iobs,ilev,ih,jvar)
            zw2(iobs,ilev) = zw2(iobs,ilev) + gpfull5%gems(iobs,ilev,ix2,jvar) * gpslant%gems(iobs,ilev,ih,jvar)  
          else
          ! Just use first profile if zenith angle is missing or really nadir
            gpfull%gems(iobs,ilev,1,jvar) = gpfull%gems(iobs,ilev,1,jvar)+gpslant%gems(iobs,ilev,ih,jvar)
          endif
          gpslant%gems(iobs,ilev,ih,jvar) = 0.0_jprb
        enddo
      enddo
    enddo
  endif

  invars = size(gpslant%pxyb(1,1,:,1))
  do jvar=1,invars
    do ilev = 1, gpfull5%nflevg
      do iobs = 1, gpfull5%ndlen
        if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
          ix1 = ix1s(iobs,ilev)
          ix2 = ix2s(iobs,ilev)
        
          gpfull%pxyb(iobs,ilev,jvar,ix1) = gpfull%pxyb(iobs,ilev,jvar,ix1) + zw1_5(iobs,ilev)*gpslant%pxyb(iobs,ilev,jvar,ih)
          gpfull%pxyb(iobs,ilev,jvar,ix2) = gpfull%pxyb(iobs,ilev,jvar,ix2) + zw2_5(iobs,ilev)*gpslant%pxyb(iobs,ilev,jvar,ih)
          zw1(iobs,ilev) = zw1(iobs,ilev) + gpfull5%pxyb(iobs,ilev,jvar,ix1)*gpslant%pxyb(iobs,ilev,jvar,ih)
          zw2(iobs,ilev) = zw2(iobs,ilev) + gpfull5%pxyb(iobs,ilev,jvar,ix2)*gpslant%pxyb(iobs,ilev,jvar,ih)

        else
        ! Just use first profile if zenith angle is missing or really nadir
          gpfull%pxyb(iobs,ilev,jvar,1) = gpfull%pxyb(iobs,ilev,jvar,1)+gpslant%pxyb(iobs,ilev,jvar,ih)
        endif
        gpslant%pxyb(iobs,ilev,jvar,ih) = 0.0_jprb
      enddo
    enddo
  enddo

  invars = size(gpslant%store_ua0(1,1,1,:))
  do jvar=1,invars
    do ilev = 1, gpfull5%nflevg
      do iobs = 1, gpfull5%ndlen
        if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
          ix1 = ix1s(iobs,ilev)
          ix2 = ix2s(iobs,ilev)

          gpfull%store_ua0(iobs,ilev,ix1,jvar) = gpfull%store_ua0(iobs,ilev,ix1,jvar) + zw1_5(iobs,ilev)*gpslant%store_ua0(iobs,ilev,ih,jvar)
          gpfull%store_ua0(iobs,ilev,ix2,jvar) = gpfull%store_ua0(iobs,ilev,ix2,jvar) + zw2_5(iobs,ilev)*gpslant%store_ua0(iobs,ilev,ih,jvar)
          zw1(iobs,ilev) = zw1(iobs,ilev) + gpfull5%store_ua0(iobs,ilev,ix1,jvar)*gpslant%store_ua0(iobs,ilev,ih,jvar)
          zw2(iobs,ilev) = zw2(iobs,ilev) + gpfull5%store_ua0(iobs,ilev,ix2,jvar)*gpslant%store_ua0(iobs,ilev,ih,jvar)

          if (ilev == 1) then
            gpfull%store_ua0(iobs,0,ix1,jvar) = gpfull%store_ua0(iobs,0,ix1,jvar) + zw1_5(iobs,ilev)*gpslant%store_ua0(iobs,0,ih,jvar)
            gpfull%store_ua0(iobs,0,ix2,jvar) = gpfull%store_ua0(iobs,0,ix2,jvar) + zw2_5(iobs,ilev)*gpslant%store_ua0(iobs,0,ih,jvar)
            zw1(iobs,ilev) = zw1(iobs,ilev) + gpfull5%store_ua0(iobs,0,ix1,jvar)*gpslant%store_ua0(iobs,0,ih,jvar)
            zw2(iobs,ilev) = zw2(iobs,ilev) + gpfull5%store_ua0(iobs,0,ix2,jvar)*gpslant%store_ua0(iobs,0,ih,jvar)
          endif
        else
        ! Just use first profile if zenith angle is missing or really nadir
          gpfull%store_ua0(iobs,ilev,1,jvar) = gpfull%store_ua0(iobs,ilev,1,jvar)+gpslant%store_ua0(iobs,ilev,ih,jvar)
          if (ilev == 1) then
            gpfull%store_ua0(iobs,0,1,jvar) = gpfull%store_ua0(iobs,0,1,jvar)+gpslant%store_ua0(iobs,0,ih,jvar)
          endif
        endif
        gpslant%store_ua0(iobs,0,ih,jvar) = 0.0_jprb
      enddo
    enddo
  enddo

  invars = size(gpslant%pxp(1,1,:,1))
  do jvar=1,invars
    do ilev = 1, gpfull5%nflevg
      do iobs = 1, gpfull5%ndlen
        if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
          ix1 = ix1s(iobs,ilev)
          ix2 = ix2s(iobs,ilev)

          gpfull%pxp(iobs,ilev,jvar,ix1) = gpfull%pxp(iobs,ilev,jvar,ix1) + zw1_5(iobs,ilev)*gpslant%pxp(iobs,ilev,jvar,ih)
          gpfull%pxp(iobs,ilev,jvar,ix2) = gpfull%pxp(iobs,ilev,jvar,ix2) + zw2_5(iobs,ilev)*gpslant%pxp(iobs,ilev,jvar,ih)
          zw1(iobs,ilev) = zw1(iobs,ilev) + gpfull5%pxp(iobs,ilev,jvar,ix1)*gpslant%pxp(iobs,ilev,jvar,ih)
          zw2(iobs,ilev) = zw2(iobs,ilev) + gpfull5%pxp(iobs,ilev,jvar,ix2)*gpslant%pxp(iobs,ilev,jvar,ih)
          gpslant%pxp(iobs,ilev,jvar,ih) = 0.0_jprb

          gpfull%pxz(iobs,ilev,jvar,ix1) = gpfull%pxz(iobs,ilev,jvar,ix1) + zw1_5(iobs,ilev)*gpslant%pxz(iobs,ilev,jvar,ih)
          gpfull%pxz(iobs,ilev,jvar,ix2) = gpfull%pxz(iobs,ilev,jvar,ix2) + zw2_5(iobs,ilev)*gpslant%pxz(iobs,ilev,jvar,ih)
          zw1(iobs,ilev) = zw1(iobs,ilev) + gpfull5%pxz(iobs,ilev,jvar,ix1)*gpslant%pxz(iobs,ilev,jvar,ih)
          zw2(iobs,ilev) = zw2(iobs,ilev) + gpfull5%pxz(iobs,ilev,jvar,ix2)*gpslant%pxz(iobs,ilev,jvar,ih)
          gpslant%pxz(iobs,ilev,jvar,ih) = 0.0_jprb

          gpfull%pxpd(iobs,ilev-1,jvar,ix1) = gpfull%pxpd(iobs,ilev-1,jvar,ix1) + zw1_5(iobs,ilev)*gpslant%pxpd(iobs,ilev-1,jvar,ih)
          gpfull%pxpd(iobs,ilev-1,jvar,ix2) = gpfull%pxpd(iobs,ilev-1,jvar,ix2) + zw2_5(iobs,ilev)*gpslant%pxpd(iobs,ilev-1,jvar,ih)
          zw1(iobs,ilev) = zw1(iobs,ilev) + gpfull5%pxpd(iobs,ilev-1,jvar,ix1)*gpslant%pxpd(iobs,ilev-1,jvar,ih)
          zw2(iobs,ilev) = zw2(iobs,ilev) + gpfull5%pxpd(iobs,ilev-1,jvar,ix2)*gpslant%pxpd(iobs,ilev-1,jvar,ih)
          gpslant%pxpd(iobs,ilev-1,jvar,ih) = 0.0_jprb

          gpfull%pxzd(iobs,ilev-1,jvar,ix1) = gpfull%pxzd(iobs,ilev-1,jvar,ix1) + zw1_5(iobs,ilev)*gpslant%pxzd(iobs,ilev-1,jvar,ih)
          gpfull%pxzd(iobs,ilev-1,jvar,ix2) = gpfull%pxzd(iobs,ilev-1,jvar,ix2) + zw2_5(iobs,ilev)*gpslant%pxzd(iobs,ilev-1,jvar,ih)
          zw1(iobs,ilev) = zw1(iobs,ilev) + gpfull5%pxzd(iobs,ilev-1,jvar,ix1)*gpslant%pxzd(iobs,ilev-1,jvar,ih)
          zw2(iobs,ilev) = zw2(iobs,ilev) + gpfull5%pxzd(iobs,ilev-1,jvar,ix2)*gpslant%pxzd(iobs,ilev-1,jvar,ih)
          gpslant%pxzd(iobs,ilev-1,jvar,ih) = 0.0_jprb
        else
        ! Just use first profile if zenith angle is missing or really nadir
          gpfull%pxp(iobs,ilev,jvar,1) = gpfull%pxp(iobs,ilev,jvar,1) + gpslant%pxp(iobs,ilev,jvar,ih)
          gpslant%pxp(iobs,ilev,jvar,ih) = 0.0_jprb
          gpfull%pxz(iobs,ilev,jvar,1) = gpfull%pxz(iobs,ilev,jvar,1) + gpslant%pxz(iobs,ilev,jvar,ih)
          gpslant%pxz(iobs,ilev,jvar,ih) = 0.0_jprb
          gpfull%pxpd(iobs,ilev,jvar,1) = gpfull%pxpd(iobs,ilev,jvar,1) + gpslant%pxpd(iobs,ilev,jvar,ih)
          gpslant%pxpd(iobs,ilev,jvar,ih) = 0.0_jprb
          gpfull%pxzd(iobs,ilev,jvar,1) = gpfull%pxzd(iobs,ilev,jvar,1) + gpslant%pxzd(iobs,ilev,jvar,ih)
          gpslant%pxzd(iobs,ilev,jvar,ih) = 0.0_jprb
          if(ilev==1)then
            gpfull%pxp(iobs,0,jvar,1) = gpfull%pxp(iobs,0,jvar,1) + gpslant%pxp(iobs,0,jvar,ih)
            gpslant%pxp(iobs,0,jvar,ih) = 0.0_jprb
            gpfull%pxz(iobs,0,jvar,1) = gpfull%pxz(iobs,0,jvar,1) + gpslant%pxz(iobs,0,jvar,ih)
            gpslant%pxz(iobs,0,jvar,ih) = 0.0_jprb
            gpfull%pxpd(iobs,0,jvar,1) = gpfull%pxpd(iobs,0,jvar,1) + gpslant%pxpd(iobs,0,jvar,ih)
            gpslant%pxpd(iobs,0,jvar,ih) = 0.0_jprb
            gpfull%pxzd(iobs,0,jvar,1) = gpfull%pxzd(iobs,0,jvar,1) + gpslant%pxzd(iobs,0,jvar,ih)
            gpslant%pxzd(iobs,0,jvar,ih) = 0.0_jprb
          endif
        endif
      enddo
    enddo
  enddo

  invars = size(gpslant%store_ua1(1,1,1,:))
  do jvar=1,invars
    do ilev = 1, gpfull5%nflevg
      do iobs = 1, gpfull5%ndlen
        if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then
          ix1 = ix1s(iobs,ilev)
          ix2 = ix2s(iobs,ilev)
        
          gpfull%store_ua1(iobs,ilev,ix1,jvar) = gpfull%store_ua1(iobs,ilev,ix1,jvar) + zw1_5(iobs,ilev)*gpslant%store_ua1(iobs,ilev,ih,jvar)
          gpfull%store_ua1(iobs,ilev,ix2,jvar) = gpfull%store_ua1(iobs,ilev,ix2,jvar) + zw2_5(iobs,ilev)*gpslant%store_ua1(iobs,ilev,ih,jvar)
          zw1(iobs,ilev) = zw1(iobs,ilev) + gpfull5%store_ua1(iobs,ilev,ix1,jvar)*gpslant%store_ua1(iobs,ilev,ih,jvar)
          zw2(iobs,ilev) = zw2(iobs,ilev) + gpfull5%store_ua1(iobs,ilev,ix2,jvar)*gpslant%store_ua1(iobs,ilev,ih,jvar)
        else
        ! Just use first profile if zenith angle is missing or really nadir
          gpfull%store_ua1(iobs,ilev,1,jvar) = gpfull%store_ua1(iobs,ilev,1,jvar) + gpslant%store_ua1(iobs,ilev,ih,jvar)
        endif
        gpslant%store_ua1(iobs,ilev,ih,jvar) = 0.0_jprb
      enddo
    enddo
  enddo

    ! Neglecting soil displacement here - maybe deal with it later
  do iobs = 1, gpfull5%ndlen
    gpfull%store_soil(iobs,:,1,:) = gpfull%store_soil(iobs,:,1,:)+gpslant%store_soil(iobs,:,ih,:)
    gpslant%store_soil(iobs,:,ih,:)  = 0.0_jprb
  enddo


    ! Neglecting surface displacement here - maybe deal with it later
  do iobs = 1, gpfull5%ndlen
    gpfull%store_sfc(iobs,1,:) = gpfull%store_sfc(iobs,1,:)+gpslant%store_sfc(iobs,ih,:)
    gpslant%store_sfc(iobs,ih,:)  = 0.0_jprb
  enddo

  do iobs = 1, gpfull5%ndlen
    if(pzenith(iobs) /= rmdi .and. pzenith(iobs) /= 0.0_jprb) then 
      zxloc = 0.0_jprb
      z_hgt = 0.0_jprb
      do ilev = 1,gpfull5%nflevg
        ix1 = ix1s(iobs,ilev)
        ix2 = ix2s(iobs,ilev)
        if (ilev > 1) then       
          gpfull%hgt(iobs,ilev-1,ix1) = gpfull%hgt(iobs,ilev-1,ix1) + zw1_5(iobs,ilev)*z_hgt
          gpfull%hgt(iobs,ilev-1,ix2) = gpfull%hgt(iobs,ilev-1,ix2) + zw2_5(iobs,ilev)*z_hgt
          zw1(iobs,ilev) = zw1(iobs,ilev) + gpfull5%hgt(iobs,ilev-1,ix1)*z_hgt
          zw2(iobs,ilev) = zw2(iobs,ilev) + gpfull5%hgt(iobs,ilev-1,ix2)*z_hgt
          z_hgt = 0.0_jprb
        endif
        zw2(iobs,ilev) = zw2(iobs,ilev) - zw1(iobs,ilev)
        zw1(iobs,ilev) = 0.0_jprb
        if(zw2_5a(iobs,ilev) > 1.0_jprb)then
          zw2(iobs,ilev)   = 0.0_jprb
        endif
        zxloc = zw2(iobs,ilev)/zdx
        zw2(iobs,ilev)  = 0.0_jprb 
        z_hgt = z_hgt +  tan(pzenith(iobs) * rpi / 180.0_jprb)*zxloc
      enddo
      if(gpfull5%hgt(iobs,gpfull5%nflevg,1) >= 0.0_jprb)then
        gpfull%hgt(iobs,gpfull5%nflevg,1) = gpfull%hgt(iobs,gpfull5%nflevg,1)+z_hgt
      endif
      
    endif
  enddo
          

endif


if (lhook) call dr_hook('gom_plus_make_slant_path_ad',1,zhook_handle)

end subroutine gom_plus_make_slant_path_ad

! -------------------------------------------------------------
! Apply wind rotation - TL and AD version
!
! -------------------------------------------------------------
subroutine rotate_uv_tlad(kobs,gp,gp5,ld_u,ld_tl)
integer(kind=jpim)   ,intent(in)    :: kobs
type(type_gom_plus)  ,intent(inout) :: gp
type(type_gom_plus)  ,intent(in)    :: gp5
logical              ,intent(in)    :: ld_u
logical              ,intent(in)    :: ld_tl

integer(kind=jpim) :: ilev, ilevbegin, ilevend
real(kind=jprb)    :: zhook_handle

!  ------------------------------------------------------------------
if (lhook) call dr_hook('gom_plus:rotate_uv_tlad',0,zhook_handle)

if (ld_u) then
   ilevbegin=1
   ilevend=gp%nflevg
else
   ilevbegin=gp%nflevg-1
   ilevend=gp%nflevg
endif

do ilev = ilevbegin, ilevend
  call xfactor_rotate_uv_tlad(gp%nhoriz, gp%uf(kobs,ilev,:), gp%vf(kobs,ilev,:), &
   & gp5%gnordl(kobs,:), gp5%gnordm(kobs,:), ld_tl)
enddo

if (lhook) call dr_hook('gom_plus:rotate_uv_tlad',1,zhook_handle)

end subroutine rotate_uv_tlad

! -------------------------------------------------------------
! Calculate wind rotation - tl and ad version
!
! -------------------------------------------------------------
subroutine xfactor_rotate_uv_tlad(khoriz,pu,pv,pnordl,pnordm,ld_tl)

integer(kind=jpim),intent(in)  :: khoriz
real(kind=jprb), intent(inout) :: pu(khoriz)
real(kind=jprb), intent(inout) :: pv(khoriz)
real(kind=jprb), intent(in)    :: pnordl(khoriz)
real(kind=jprb), intent(in)    :: pnordm(khoriz)
logical,         intent(in)    :: ld_tl

real(kind=jprb) ,dimension(khoriz) :: zdl2, zdm2, zul, zvl

if(ld_tl) then
  zdl2(:) = pu(:) * pnordm(:) - pv(:) * pnordl(:)
  zdm2(:) = pu(:) * pnordl(:) + pv(:) * pnordm(:)
else
  zdl2(:) = pu(:) * pnordm(:) + pv(:) * pnordl(:)
  zdm2(:) = - pu(:) * pnordl(:) + pv(:) * pnordm(:)
endif
pu(:)=zdl2(:)
pv(:)=zdm2(:)

end subroutine xfactor_rotate_uv_tlad

end module gom_plus


