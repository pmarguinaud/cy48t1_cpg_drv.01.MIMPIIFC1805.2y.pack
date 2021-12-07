#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE MF_PHYS_MPA(YDGEOMETRY, YDCPG_DIM, YDCPG_MISC, YDCPG_PHY0, YDCPG_PHY9, YDMF_PHYS,  YDMF_PHYS_TMP, &
& YDCPG_DYN0, YDCPG_DYN9, YDMF_PHYS_SURF, YDVARS, YDGMV, YDSURF, YDCFU, YDXFU, YDMODEL, LDCONFX, PDTPHY, &
& PGFL, PKOZO, PGP2DSDT, PB1, PB2, PGMVT1, PGFLT1, PGPAR, PTRAJ_PHYS, YDDDH, PFTCNS)

!**** *MF_PHYS_MPA* METEO-FRANCE PHYSICS.

!     Purpose.
!     --------
!         Call METEO-FRANCE physics and physical tendencies.

!**   Interface.
!     ----------
!        *CALL* *MF_PHYS_MPA(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!        KBL       : NPROMA-packets number
!        KGPCOMP   : total number of grid points in the domain
!        KST       : first element of work.
!        KEND      : last element of work.
!        KSTGLO    : global offset.
!        LDCONFX   : (see in CPG)
!        PDTPHY    : timestep used in the physics.
!        KIBL      : index into YRCSGEOM/YRGSGEOM types in YDGEOMETRY
!        POROGL,POROGM: components of grad(orography).
!        PGFL      : GFL at time t and t-dt.
!        PKOZO     : fields for photochemistery of ozon.
!        PGP2DSDT  : stochastic physics random pattern.

!     INPUT/OUTPUT:
!     -------------
!        PB1       : "SLB1"-buffer, used for interpolations in the SL scheme.
!        PB2       : "SLB2"-buffer.
!        PGFLT1    : GFL t+dt
!        PGPAR     : surface fields for AROME.
!        PGMU0     : COSINE OF SOLAR ZENITH ANGLE, APPROXIMATE ACTUAL VALUE
!                    linear T_e correction
!                    linear T_e correction

!     OUTPUT:
!     -------
!        PDHSF     : distribution of horizontal mean weights used for
!                    simplified radiation scheme.
!        ---------------------- output of aplpar ------------------------------
!        PFCQNG    : pseudo-flux of water to correct for Q<0.
!        PDIFCQLC to PFCNEGQSC:
!        ---------------------- end of output of aplpar -----------------------
!        PTENDU    : "U"-wind tendency due to physics.
!        PTENDV    : "V"-wind tendency due to physics.
!        PDIAGH    : Add Hail diagnostic PDIAGH (AROME)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!       2000-12-04: F. Bouyssel & J.M. Piriou

!     Modifications.
!     --------------
!       04-Mar-2009 A.Alias : call CPTEND/INITAPLPAR modified to add
!                         Humidity Mesopheric flux (ZFRMQ).
!                     and IVCLIA removed and call to CPNUDG modified as
!                         Nuding mask is now in SD_VF group
!                         call HL_APLPAR modified to add PFCQNG for acdifus
!                         call APL_AROME modified to add Sulfate/Volcano aerosols for radaer
!       K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!       2009-10-15 Y. Bouteloup : Store radiative cloud water and ice in GFL (YIRAD and YLRAD)
!       F. Vana   15-Oct-2009 : NSPLTHOI option
!       K.Yessad (Feb 2010): use YM_RADTC and RFORADTC
!       2010-03-26 Y. Bouteloup : Store radiative cloud water and ice in GFL (AROME case)
!       2010-04-26 Y. Bouteloup : Only one call to cputqy, cputqys and cputqy_arome
!            This need the use of ZTENDGFL as argument of cptend, cptend_new and apl_arome.
!       2010-05-11 F. Bouyssel : Use of PINDX, PINDY
!       2010-05-28 C. Geijo    : Fix error in IPTR array element referencing 
!       2010-06-21 O.Riviere/F. Bouyssel : Fix to have Ts evolving in Fa files with Surfex
!       Dec 2010 A.Alias   : ZMU0N added to call CPPHINP/APLPAR/APL_AROME/HL_APLPAR
!                            CALL to CPNUDG with or with LMSE (A.Voldoire)
!       K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!       L. Bengtsson-Sedlar & F. Vana 18-Feb-2011 : CA scheme for convection
!       F. Vana   22-Feb-2011 : 3D turbulence
!       2011-02-01 M. Mokhtari: Add LMDUST and PEXTT9 and PEXTT0 IN APLPAR
!                             (treatment of the desert aerosols) 
!       2011-03 A.Alias  : new argument to  for sunshine hours YSD_VD%YSUND 
!                      CPNUDG if LMSE=.T. or LMSE=.F. (bugfix)
!                      debug ozone GFL (IPO3) (D. St-Martin)
!                      Humidity Mesopheric flux (ZFRMQ) added in CPTEND_NEW
!       F.Bouyssel (26-03-2011): Fix to have Snow in hist file with surfex
!       2011-06: M. Jerczynski - some cleaning to meet norms
!       E. Bazile 2011-08-26 : Output for MUSC 1D with LFA files with WRITEPHYSIO
!         used previously for extracting profiles from 3D (now also available for AROME).
!       K. Yessad (Dec 2011): use YDOROG, YDGSGEOM and YDCSGEOM.
!       2011-11-21 JF Gueremy : dry convective adjustment (LAJUCV)
!       F. Vana  26-Jan-2012 : historic Qs for TOM's BBC.
!       F.Bouttier Jul 2012: stochastic physics for AROME
!       Z. SASSI  : 07-Mar-2013   INITIALIZING THE WEIGHT VECTORS PDHSF(NPROMA)
!       [DISTRIBUTION OF HORIZONTAL MEANS WEIGHTS]
!       F. Vana  28-Nov-2013 : Redesigned trajectory handling
!       T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!       2013-11, D. Degrauwe: Flexible interface CPTEND_FLEX.
!       2013-11, J. Masek: Passing intermittency arrays for ACRANEB2.
!       K. Yessad (July 2014): Move some variables.
!       2016-04, J. Masek: Passing sunshine duration to APL_AROME.
!       2016-09, M. Mokhtari & A. Ambar: replacement of ZEXT and ZEZDIAG by PGFL
!                                        in aplpar.F90 argument.
!       2016-10, P. Marguinaud : Port to single precision
!       K. Yessad (Dec 2016): Prune obsolete options.
!       K. Yessad (June 2017): Introduce NHQE model.
!       2017-09, J. Masek: Shifted dimensioning of PGMU0.
!       K. Yessad (Feb 2018): remove deep-layer formulations.
!       K. Yessad (Apr 2018): introduce key L_RDRY_VD (ensure consistent definition of "dver" everywhere).
!       2018-09, F. Duruisseau: add rconv and sconv in gfl for bayrad
!       2018-09, R. Brozkova: Passing of diagnostic hail, global normal
!         irradiance and mean radiant temperature from APLPAR.
!       2018-09, D. St-Martin : add NOGWD inputs in aplpar
!       2018-09, M. Michou : add ARPEGE-Climat chemistry call in aplpar  
!   R. El Khatib 27-02-2019 Use pointer function SC2PRG to avoid bounds violation
!       2019-05, I. Etchevers : add visibilities and precipitation type
!   R. El Khatib 27-02-2019 memory bandwidth savings.
!   R. El Khatib 30-Oct-2018 IMAXDRAFT
!       2019-09, M. Hrastinski: Dataflow for TKE and TTE terms in ALARO DDH (PFTCNS).
!       2019-09, J. Masek: Modified call to APL_AROME (added argument NFRRC).
!       2019-12, Y. Bouteloup: Introduction of ZTENDU and ZTENDV for computation of ZDEC in cputqy
!                diferent from PTENDU and PTENDV in the case of the use of Tiedtke scheme to avoid double counting
!       2020-12, U. Andrae : Introduce SPP for HARMONIE-AROME
!       2021-01, R. Brozkova: ALARO graupel fix.
! End Modifications
!-------------------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE MF_PHYS_TYPE_MOD   , ONLY : MF_PHYS_TYPE, MF_PHYS_TMP_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_DYN_TYPE, CPG_PHY_TYPE, &
                              & CPG_MISC_TYPE
USE CPG_DIM_TYPE_MOD   , ONLY : CPG_DIM_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD,ONLY : MF_PHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMGMV             , ONLY : TGMV
USE YOMCFU             , ONLY : TCFU
USE YOMXFU             , ONLY : TXFU
USE TYPE_MODEL         , ONLY : MODEL
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK
USE SC2PRG_MOD         , ONLY : SC2PRG

USE YOMCT0             , ONLY : LSLAG, LTWOTL, LNHDYN, LAROME, LSFORCS, LNHQE
USE YOMCVER            , ONLY : LVERTFE  ,LVFE_GWMPA 
USE YOMDYNA            , ONLY : LGWADV, L_RDRY_VD
USE YOMSCM             , ONLY : LGSCM
USE YOMCST             , ONLY : RG, RD
USE YOMLSFORC          , ONLY : LMUSCLFA
USE YOMSPSDT           , ONLY : YSPPT
USE YOMTRAJ            , ONLY : TRAJ_PHYS_TYPE

USE DDH_MIX            , ONLY : TYP_DDH
USE INTFLEX_MOD        , ONLY : LINTFLEX, TYPE_INTPROCSET, NEWINTPROCSET, CLEANINTPROCSET
USE SPP_MOD , ONLY : YSPP_CONFIG,YSPP
USE MF_PHYS_STATE_TYPE_MOD &
                       , ONLY : MF_PHYS_STATE_TYPE
!     -------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_DIM_TYPE),INTENT(IN)    :: YDCPG_DIM
TYPE(CPG_MISC_TYPE),INTENT(INOUT):: YDCPG_MISC
TYPE(CPG_PHY_TYPE),INTENT(INOUT), TARGET :: YDCPG_PHY0
TYPE(CPG_PHY_TYPE),INTENT(INOUT), TARGET :: YDCPG_PHY9
TYPE(MF_PHYS_TYPE),INTENT(INOUT) :: YDMF_PHYS
TYPE(MF_PHYS_TMP_TYPE),INTENT(INOUT) :: YDMF_PHYS_TMP
TYPE(CPG_DYN_TYPE),INTENT(INOUT), TARGET :: YDCPG_DYN0
TYPE(CPG_DYN_TYPE),INTENT(INOUT), TARGET :: YDCPG_DYN9
TYPE(MF_PHYS_SURF_TYPE),INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
TYPE(TCFU)        ,INTENT(INOUT) :: YDCFU
TYPE(TXFU)        ,INTENT(INOUT) :: YDXFU
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
LOGICAL           ,INTENT(IN)    :: LDCONFX
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTPHY 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKOZO(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG*YDMODEL%YRML_PHY_G%YRDPHY%NVCLIS+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP2DSDT(YDCPG_DIM%KLON,YSPPT%YGPSDT(1)%NG2D)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1(YDCPG_DIM%KLON,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB2(YDCPG_DIM%KLON,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDGMV%YT1%NDIM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLT1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGPAR(YDCPG_DIM%KLON,YDMODEL%YRML_PHY_MF%YRPARAR%NGPAR+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFTCNS(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,6)
TYPE (TRAJ_PHYS_TYPE), INTENT(INOUT) :: PTRAJ_PHYS
TYPE(TYP_DDH)     ,INTENT(INOUT) :: YDDDH

!     ------------------------------------------------------------------
LOGICAL :: LLDIAB
LOGICAL :: LL_SAVE_PHSURF
LOGICAL :: LLXFUMSE

INTEGER(KIND=JPIM) :: IPTREXT,IEFB1,IEFB2,IEFB3
INTEGER(KIND=JPIM) :: IPTR(YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
INTEGER(KIND=JPIM) :: IPTRLIMA
INTEGER(KIND=JPIM) :: IRR ! pointer of 1st hydrometeors in ZTENDGFLR
INTEGER(KIND=JPIM) :: IPTRTKE ! pointer of TKE in ZTENDGFLR

INTEGER(KIND=JPIM) :: IPGFL(YDMODEL%YRML_GCONF%YGFL%NUMFLDS)

INTEGER(KIND=JPIM) :: INSTEP_DEB,INSTEP_FIN
INTEGER(KIND=JPIM) :: JLEV, JGFL
INTEGER(KIND=JPIM) :: JROF
INTEGER(KIND=JPIM) :: ISLB1U9  ,ISLB1V9  ,ISLB1T9  ,ISLB1GFL9, ISLB1VD9

!     --- UPPER AIR PHYSICAL TENDENCIES.
REAL(KIND=JPRB) :: ZTENDH(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)     ! Enthalpy tendency.
REAL(KIND=JPRB) :: ZTENDPTKE(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)  ! Pseudo progn. TKE

! GFL tendencies for APL_AROME (assumes YDMODEL%YRML_GCONF%YGFL%NUMFLDS>=YDMODEL%YRML_PHY_MF%YRPARAR%NRR)
! for now, use Jovi's trick :
REAL(KIND=JPRB), TARGET :: ZTENDGFLR(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,0:YDMODEL%YRML_GCONF%YGFL%NUMFLDS)
REAL(KIND=JPRB) :: ZTENDGFL(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NUMFLDS)   ! GFL tendencies

!     --- UPPER AIR PHYSICAL TENDENCIES FOR AROME.
!       (the previous one are not used in AROME)
REAL(KIND=JPRB) :: ZTENDT (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)        ! temperature tendency
REAL(KIND=JPRB) :: ZTENDW (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)        ! W  tendency
REAL(KIND=JPRB) :: ZTENDD (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)        ! d  tendency
REAL(KIND=JPRB) :: ZDIFEXT(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)     ! Extra-GFL fluxes.

REAL(KIND=JPRB) :: ZTENDU (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! U tendency without deep convection contribution
REAL(KIND=JPRB) :: ZTENDV (YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)    ! V tendency without deep convection contribution


!     ---FOR AROME PHYSICS  ---
REAL(KIND=JPRB) :: ZGWT1(YDCPG_DIM%KLON,0:YDCPG_DIM%KFLEVG)     ! vertical velocity calculated by cputqy_arome before convertion in d
REAL(KIND=JPRB) :: ZTT1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)       ! Temperature at t1

! ZRTT1: appropriate version of R*T at t1 for gnhgw2svd
!  Version of R must be consistent with definition of vertical divergence.
REAL(KIND=JPRB) :: ZRTT1(YDCPG_DIM%KLON,YDCPG_DIM%KFLEVG)

! required for INTFLEX
TYPE(TYPE_INTPROCSET) :: YLPROCSET

! SPP
REAL(KIND=JPRB) :: ZGP2DSPP(YDCPG_DIM%KLON,YSPP%N2D)

INTEGER(KIND=JPIM) :: IMAXDRAFT

REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDR    (:,:,:) 
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDLIMA (:,:,:) 
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDTKE  (:,:) 
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDEFB1 (:,:) 
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDEFB2 (:,:)
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDEFB3 (:,:)
REAL(KIND=JPRB), POINTER, CONTIGUOUS :: ZTENDEXT  (:,:,:)

REAL(KIND=JPRB), POINTER :: ZP1CHEM0(:,:,:), ZP1CHEM9(:,:,:)
REAL(KIND=JPRB), POINTER :: ZP1EXT0(:,:,:), ZP1EXT9(:,:,:)
REAL(KIND=JPRB), POINTER :: ZP1LIMA0(:,:,:), ZP1LIMA9(:,:,:)
REAL(KIND=JPRB), POINTER :: ZP1EZDIAG(:,:)
REAL(KIND=JPRB), POINTER :: ZP1NOGW0(:,:), ZP1NOGW9(:,:), ZP2NOGW0(:,:), ZP2NOGW9(:,:)

TYPE (MF_PHYS_STATE_TYPE) :: YLMF_PHYS_STATE

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "apl_arome.intfb.h"
#include "cpphinp.intfb.h"
#include "cptend_flex.intfb.h"
#include "cputqy_arome.intfb.h"
#include "cputqy.intfb.h"
#include "cp_ptrslb1.intfb.h"
#include "gnhgw2svdarome.intfb.h"
#include "mf_phys_init.intfb.h"
#include "writephysio.intfb.h"
#include "mf_phys_nhqe_part1.intfb.h"
#include "mf_phys_nhqe_part2.intfb.h"
#include "mf_phys_save_phsurf_part1.intfb.h"
#include "mf_phys_save_phsurf_part2.intfb.h"
#include "mf_phys_transfer.intfb.h"
#include "mf_phys_precips.intfb.h"
#include "apl_arome_calc_iptr.intfb.h"
#include "apl_arome_calc_ipgfl.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MF_PHYS_MPA',0,ZHOOK_HANDLE)
ASSOCIATE(YDCSGEOM=> YDGEOMETRY%YRCSGEOM(YDCPG_DIM%KBL), YDGSGEOM=>YDGEOMETRY%YRGSGEOM(YDCPG_DIM%KBL), &
& YDOROG=>YDGEOMETRY%YROROG(YDCPG_DIM%KBL), YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,  &
& YDPTRSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1, YDRIP=>YDMODEL%YRML_GCONF%YRRIP, YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH,  &
& YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY, YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY,     &
& YDLDDH=>YDMODEL%  YRML_DIAG%YRLDDH, YGFL=>YDMODEL%YRML_GCONF%YGFL,           &
& YDEPHY=> YDMODEL%YRML_PHY_EC%YREPHY, YDPARAR=>YDMODEL%YRML_PHY_MF%YRPARAR)

ASSOCIATE(MVTS=>YDPARAR%MVTS, CMF_UPDRAFT=>YDPARAR%CMF_UPDRAFT, NTSSG=>YDDPHY%NTSSG, LMSE=>YDARPHY%LMSE,    &
& LMFSHAL=>YDARPHY%LMFSHAL, YEZDIAG=>YGFL%YEZDIAG, YEXT=>YGFL%YEXT, YNOGW=>YGFL%YNOGW, YCHEM=>YGFL%YCHEM, &
& YSP_SBD=>YDSURF%YSP_SBD, LEDR=>YDPHY%LEDR, LAGPHY=>YDEPHY%LAGPHY, YLIMA=>YGFL%YLIMA)

CALL SC2PRG(1,YCHEM(:)%MP ,YDMODEL%YRML_GCONF%YGFL%NCHEM,   PGFL,ZP1CHEM0)    !  YDVARS%CHEM(1)%T0
CALL SC2PRG(1,YCHEM(:)%MP9,YDMODEL%YRML_GCONF%YGFL%NCHEM,   PGFL,ZP1CHEM9)    !  YDVARS%CHEM(1)%T9
CALL SC2PRG(1,YLIMA(:)%MP ,YDMODEL%YRML_GCONF%YGFL%NLIMA,   PGFL,ZP1LIMA0)    !  YDVARS%LIMA(1)%T0
CALL SC2PRG(1,YLIMA(:)%MP9,YDMODEL%YRML_GCONF%YGFL%NLIMA,   PGFL,ZP1LIMA9)    !  YDVARS%LIMA(1)%T9
CALL SC2PRG(1,YEXT(:)%MP  ,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT,PGFL,ZP1EXT0 )    !  YDVARS%EXT(1)%T0
CALL SC2PRG(1,YEXT(:)%MP9 ,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT,PGFL,ZP1EXT9 )    !  YDVARS%EXT(1)%T9
CALL SC2PRG(1,YEZDIAG(:)%MP,PGFL,ZP1EZDIAG)
CALL SC2PRG(1,YNOGW(:)%MP ,PGFL ,ZP1NOGW0)    !  YDVARS%NOGW(1)%T0
CALL SC2PRG(1,YNOGW(:)%MP9,PGFL ,ZP1NOGW9)    !  YDVARS%NOGW(1)%T9 
CALL SC2PRG(2,YNOGW(:)%MP ,PGFL ,ZP2NOGW0)    !  YDVARS%NOGW(2)%T0
CALL SC2PRG(2,YNOGW(:)%MP9,PGFL ,ZP2NOGW9)    !  YDVARS%NOGW(2)%T9 

CALL YLMF_PHYS_STATE%INIT (LTWOTL, YDCPG_DYN0, YDCPG_DYN9, YDCPG_PHY0, YDCPG_PHY9, YDVARS, YDMF_PHYS_SURF, &
                  & P1EXT0=ZP1EXT0, P1EXT9=ZP1EXT9, P1CHEM0=ZP1CHEM0, P1CHEM9=ZP1CHEM9, &
                  & P1NOGW0=ZP1NOGW0, P1NOGW9=ZP1NOGW9, P2NOGW0=ZP2NOGW0, P2NOGW9=ZP2NOGW9, &
                  & P1LIMA0=ZP1LIMA0, P1LIMA9=ZP1LIMA9)

!     ------------------------------------------------------------------

!        0.    constructor for procset
IF (LINTFLEX) YLPROCSET=NEWINTPROCSET()

!        1.    Preliminary calculations necessary
!              for all types of physics.
!              ------------------------------------


INSTEP_DEB=1
INSTEP_FIN=1

! initialisation for surfex if XFU
LLXFUMSE=.FALSE.
IF (LDCONFX) THEN
  LLXFUMSE=.TRUE.
ENDIF

! SPP 
IF ( YSPP_CONFIG%LSPP ) THEN
 DO JROF=1,YSPP%N2D
   ZGP2DSPP(:,JROF) = YSPP%GP_ARP(JROF)%GP2D(:,1,YDCPG_DIM%KBL)
 ENDDO
ENDIF

! Complete physics is called.
LLDIAB=(.NOT.LAGPHY)

! In the NHQE model, MF_PHYS_MPA enters with Tt and grad(Tt), where Tt = T * exp(-(R/cp) log(pre/prehyd)).
! But calculations of MF_PHYS_MPA must use T and grad(T).
! So we do a conversion Tt -> T.
IF (LNHQE) THEN
  CALL MF_PHYS_NHQE_PART1 (YDGEOMETRY, YDCPG_DIM, YDMF_PHYS_TMP, YDVARS, YDMODEL, PGFL)
ENDIF

IF (LLDIAB) THEN
  CALL CPPHINP(YDGEOMETRY, YDMODEL, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDGSGEOM%GEMU, YDGSGEOM%GELAM, YDVARS%U%T0, YDVARS%V% &
  & T0, YDVARS%Q%T0, YDVARS%Q%DL, YDVARS%Q%DM, YDVARS%CVGQ%DL, YDVARS%CVGQ%DM, YDCPG_PHY0%XYB%RDELP, &
  & YDCPG_DYN0%CTY%EVEL, YDVARS%CVGQ%T0, YDMF_PHYS_TMP%RDG%MU0, YDMF_PHYS_TMP%RDG%MU0LU, YDMF_PHYS_TMP%RDG%MU0M, YDMF_PHYS_TMP%RDG%MU0N, YDMF_PHYS_TMP%RDG%CVGQ)
  YDMF_PHYS_TMP%RDG%LCVQ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=YDMF_PHYS_TMP%RDG%CVGQ(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
ENDIF

! * In some cases, some pseudo-historic surface buffers (like z0) should
!   not be modified between the entrance and the output of MF_PHYS_MPA
!   (this is the case for example if LDCONFX=T).
!   For the time being, we must save:
!   - HV (group VV) : resistance to evapotranspiration
!   - Z0F (group VD): gravity * surface roughness length
!   - Z0H (group VV): gravity * roughness length for heat
!   - PBLH (group VH): PBL height
!   - SPSH (group VH):
!   - QSH (group VH):

LL_SAVE_PHSURF = .FALSE.

IF (LLDIAB) THEN
  LL_SAVE_PHSURF=LDCONFX
  IF (LL_SAVE_PHSURF) THEN
    CALL MF_PHYS_SAVE_PHSURF_PART1 (YDCPG_DIM, YDMF_PHYS_TMP, YDMF_PHYS_SURF, YDVARS, YDSURF, YDMODEL)
  ENDIF
ENDIF


CALL MF_PHYS_INIT ( YGFL, YDARPHY, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA, YDCPG_DIM%KLON, YDCPG_DIM%KFLEVG, NTSSG, YSP_SBD%NLEVS,&
  & YDMF_PHYS%OUT%DIFCQ , YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS , ZDIFEXT, YDMF_PHYS%OUT%DIFTQ , YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL,&
  & YDMF_PHYS%OUT%DIFTS , YDMF_PHYS%OUT%FCCQL , YDMF_PHYS%OUT%FCCQN , YDMF_PHYS%OUT%FCSQL , YDMF_PHYS%OUT%FCSQN , YDMF_PHYS%OUT%FCQNG , YDMF_PHYS%OUT%FCQNNG,&
  & YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG,YDMF_PHYS%OUT%FCQGNG,&
  & YDMF_PHYS%OUT%FPLCL , YDMF_PHYS%OUT%FPLCN , YDMF_PHYS%OUT%FPLCG , YDMF_PHYS%OUT%FPLCH, YDMF_PHYS%OUT%FPLSL , YDMF_PHYS%OUT%FPLSN ,YDMF_PHYS%OUT%FPLSG ,&
  & YDMF_PHYS%OUT%FPLSH,  YDMF_PHYS%OUT%FRSO  , YDMF_PHYS%OUT%FRSOC ,&
  & YDMF_PHYS%OUT%FRTH  , YDMF_PHYS%OUT%FRTHC , YDMF_PHYS%OUT%STRCU , YDMF_PHYS%OUT%STRCV , YDMF_PHYS%OUT%STRDU , YDMF_PHYS%OUT%STRDV , YDMF_PHYS%OUT%STRTU ,&
  & YDMF_PHYS%OUT%STRTV , YDMF_PHYS%OUT%STRMU , YDMF_PHYS%OUT%STRMV , YDMF_PHYS%OUT%FRMH  , &
  & YDMF_PHYS%OUT%DIFCQLC,YDMF_PHYS%OUT%DIFCQIC,YDMF_PHYS%OUT%FIMCC,&
  & YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC,YDMF_PHYS%OUT%FCNEGQIC,YDMF_PHYS%OUT%FCNEGQRC,YDMF_PHYS%OUT%FCNEGQSC,&
  & YDMF_PHYS%OUT%FCHOZ , &
  & YDCPG_MISC%NEB   , YDCPG_MISC%QICE  , YDCPG_MISC%QLI   , &
  & YDCPG_MISC%RH    , YDMF_PHYS%OUT%ALB , &
  & YDMF_PHYS%OUT%CT    , YDMF_PHYS%OUT%FCHSP , YDMF_PHYS%OUT%FCLL  , YDMF_PHYS%OUT%FCLN  ,&
  & YDMF_PHYS%OUT%FCS   , YDMF_PHYS%OUT%FEVL  , YDMF_PHYS%OUT%FEVN  , YDMF_PHYS%OUT%FEVV  , YDMF_PHYS%OUT%FLASH , YDMF_PHYS%OUT%FTR   , YDMF_PHYS%OUT%FLWSP ,&
  & YDMF_PHYS%OUT%FONTE , YDMF_PHYS%OUT%FGEL  , YDMF_PHYS%OUT%FGELS ,&
  & YDMF_PHYS%OUT%FRSGNI, YDMF_PHYS%OUT%FRSDNI, YDMF_PHYS%OUT%FRSODS, YDMF_PHYS%OUT%FRSOPS, YDMF_PHYS%OUT%FRSOPT, YDMF_PHYS%OUT%FRSOLU, YDMF_PHYS%OUT%FRTHDS,&
  & YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN,YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN,YDMF_PHYS%OUT%FPEVPSL,YDMF_PHYS%OUT%FPEVPSN,YDMF_PHYS%OUT%FPEVPSG,YDMF_PHYS%OUT%FPEVPCL,&
  & YDMF_PHYS%OUT%FPEVPCN,YDMF_PHYS%OUT%FPEVPCG, &
  & YDMF_PHYS%OUT%GZ0   , YDMF_PHYS%OUT%GZ0H  , YDCPG_MISC%QS    , YDMF_PHYS%OUT%RUISL ,&
  & YDMF_PHYS%OUT%RUISP , YDMF_PHYS%OUT%RUISS , YDMF_PHYS%OUT%UCLS  , YDMF_PHYS%OUT%VCLS  , YDMF_PHYS%OUT%NUCLS , YDMF_PHYS%OUT%NVCLS , YDMF_PHYS%OUT%TCLS  , YDMF_PHYS%OUT%MRT,&
  & YDMF_PHYS%OUT%QCLS  , YDMF_PHYS%OUT%RHCLS , YDCPG_MISC%CLCT  , YDMF_PHYS%OUT%CLCH  , YDMF_PHYS%OUT%CLCM  , YDMF_PHYS%OUT%CLCL  , YDMF_PHYS%OUT%CLCC  ,&
  & YDMF_PHYS%OUT%CAPE  , YDMF_PHYS%OUT%CTOP  , YDMF_PHYS%OUT%CLPH  , YDMF_PHYS%OUT%VEIN  , YDMF_PHYS%OUT%UGST  , YDMF_PHYS%OUT%VGST  ,&
  & YDMF_PHYS%OUT%DIAGH , YDMF_PHYS%OUT%EDR, YDMF_PHYS%OUT%VISICLD, YDMF_PHYS%OUT%VISIHYD, YDMF_PHYS%OUT%MXCLWC, YDMF_PHYS%OUT%MOCON)

CALL APL_AROME_CALC_IPGFL (YDGEOMETRY, YDCPG_DIM, YDMODEL, IPGFL)

CALL MF_PHYS_TRANSFER (YDCPG_DIM, YDVARS, YDMODEL)

CALL APL_AROME_CALC_IPTR (YDMODEL, IEFB1, IEFB2, IEFB3, IPTR, IPTREXT, IPTRLIMA, IPTRTKE, IRR)

! If an incorrect address is used, then the initialization below will detect it :
ZTENDGFLR(:,:,0)=HUGE(1._JPRB)

IF (LMFSHAL .AND. CMF_UPDRAFT=='DUAL') THEN
  IMAXDRAFT=3
ELSE
  IMAXDRAFT=0
ENDIF

ZTENDR    => ZTENDGFLR (:, :, IRR:IRR+YDMODEL%YRML_PHY_MF%YRPARAR%NRR-1)
ZTENDLIMA => ZTENDGFLR (:, :, IPTRLIMA:IPTRLIMA+YDMODEL%YRML_GCONF%YGFL%NLIMA-1)
ZTENDTKE  => ZTENDGFLR (:, :, IPTRTKE)
ZTENDEFB1 => ZTENDGFLR (:, :, IEFB1)
ZTENDEFB2 => ZTENDGFLR (:, :, IEFB2)
ZTENDEFB3 => ZTENDGFLR (:, :, IEFB3)
ZTENDEXT  => ZTENDGFLR (:, :, IPTREXT:IPTREXT+YDMODEL%YRML_GCONF%YGFL%NGFL_EXT-1)

CALL APL_AROME(YDGEOMETRY, YDCPG_DIM, YLMF_PHYS_STATE, YDCPG_MISC, YDMF_PHYS, YDMF_PHYS_TMP, &
 & YDMF_PHYS_SURF, YDVARS, YDSURF, YDCFU, YDXFU, YDMODEL, IMAXDRAFT, NTSSG,      &
 & YDCFU%NFRRC, PDTPHY, LLXFUMSE, YDCSGEOM%RINDX, YDCSGEOM%RINDY, YDGSGEOM%GEMU, YDGSGEOM%GELAM,    &
 & YDOROG%OROG, YDGSGEOM%GM, YDGSGEOM%GECLO, YDGSGEOM%GESLO, PGP2DSDT, ZGP2DSPP, PGPAR,       &
 & ZTENDT, ZTENDR, ZTENDW, ZTENDLIMA, ZTENDTKE, &
 & ZTENDEFB1, ZTENDEFB2, ZTENDEFB3, ZTENDEXT, ZP1EZDIAG, YLPROCSET, YDDDH )

!Save surface temperature
IF (LMSE.OR.LSFORCS) THEN
  IF (LLXFUMSE) THEN
    DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS_SURF%GSP_RR%PT_T0(JROF)=PGPAR(JROF,MVTS)
    ENDDO
  ELSE
    DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      YDMF_PHYS_SURF%GSP_RR%PT_T1(JROF)=PGPAR(JROF,MVTS)
    ENDDO
  ENDIF 
ENDIF  
!      4.2  COMPUTE THE PHYS. TENDENCY FOR "T" AND "w"
!           ------------------------------------------

IF (LVERTFE.AND.LVFE_GWMPA) THEN
  ! * case LVFE_GWMPA not yet coded.
  !   (in this case ZGWT1 must be computed at full levels and
  !   not at half levels)
  CALL ABOR1(' MF_PHYS_MPA: case LVFE_GWMPA not yet coded if LMPA=T!')
ENDIF

! * compute ZTT1:
IF (LSLAG.AND.LTWOTL) THEN
  DO JLEV=1,YDCPG_DIM%KFLEVG
    DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZTT1(JROF,JLEV)=YDVARS%T%T0(JROF,JLEV)+PDTPHY*ZTENDT(JROF,JLEV)
    ENDDO
  ENDDO
ELSE
  DO JLEV=1,YDCPG_DIM%KFLEVG
    DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZTT1(JROF,JLEV)=YDVARS%T%T9(JROF,JLEV)+PDTPHY*ZTENDT(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

! * compute ZGWT1 = tendency of gw:
IF (LNHDYN) THEN
  ! Valid for LVFE_GWMPA=F only; ZGWT1 assumed to be half level values.
  DO JLEV=1,YDCPG_DIM%KFLEVG-1
    DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
      ZGWT1(JROF,JLEV)=0.5_JPRB*RG*(ZTENDW(JROF,JLEV)+ZTENDW(JROF,JLEV+1))
    ENDDO
  ENDDO
  DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
    ZGWT1(JROF,YDCPG_DIM%KFLEVG)=0.0_JPRB
    ZGWT1(JROF,0)=0.0_JPRB
  ENDDO
ENDIF

! * convert gw tendency in d tendency:
IF(LNHDYN) THEN

  IF (LGWADV) THEN
    ZTENDD(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=ZGWT1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
  ELSE

    ! * Provide the appropriate version of (RT) at t+dt for GNHGW2SVDAROME:
    IF (L_RDRY_VD) THEN
      ! Use Rd because "dver" is currently defined with Rd.
      ZRTT1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)=RD*ZTT1(YDCPG_DIM%KIDIA:YDCPG_DIM%KFDIA,1:YDCPG_DIM%KFLEVG)
    ELSE
      ! Use "moist R" because "dver" is defined with "moist R".
      ! Unfortunately, R(t+dt) is not yet available there, use R(t) instead.
      ! "Moist R" tendency is neglected in the below call to GNHGW2SVDAROME.
      DO JLEV=1,YDCPG_DIM%KFLEVG
        DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
          ZRTT1(JROF,JLEV)=YDCPG_DYN0%RCP%R(JROF,JLEV)*ZTT1(JROF,JLEV)
        ENDDO
      ENDDO
    ENDIF

    ! * Do conversion:
    IF (LSLAG.AND.LTWOTL) THEN
      CALL GNHGW2SVDAROME(YDGEOMETRY,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_PHY0%PREHYDF,YDCPG_PHY0%XYB%LNPR,ZRTT1,YDCPG_PHY0%PREF,ZGWT1,&
       & ZTENDD)  
    ELSE
      CALL GNHGW2SVDAROME(YDGEOMETRY,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,YDCPG_PHY9%PREHYDF,YDCPG_PHY9%XYB%LNPR,ZRTT1,YDCPG_PHY9%PREF,ZGWT1,&
       & ZTENDD)  
    ENDIF

  ENDIF
ELSE
  ZTENDD=0.0_JPRB
ENDIF

!      4.3  PUT THE TENDENCIES IN PB1/GFLT1/GMVT1.
!           --------------------------------------


IF (LSLAG) CALL CP_PTRSLB1(YDMODEL%YRML_DYN%YRDYN,YDPTRSLB1,ISLB1U9,ISLB1V9,ISLB1T9,ISLB1VD9,ISLB1GFL9)
IF ( LINTFLEX ) THEN

  ! Set GFL tendencies to 0
  ZTENDGFL(:,:,:) = 0.0_JPRB

  CALL CPTEND_FLEX( YDLDDH, YDMDDH, YGFL, YDPHY, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, YDCPG_DIM%KFDIA,&
  & YDCPG_DIM%KFLEVG, YDGSGEOM%GNORDL, &
  & YDGSGEOM%GNORDM, YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_STATE%YCPG_PHY%XYB%RDELP,       &
  & YLMF_PHYS_STATE%YCPG_DYN%RCP%CP, YLMF_PHYS_STATE%U, YLMF_PHYS_STATE%V, YLMF_PHYS_STATE%T,     &
  & YLMF_PHYS_STATE%YGSP_RR%T, PGFL, YLPROCSET, YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDH, &
  & ZTENDGFL, YDMF_PHYS%OUT%FHSCL, YDMF_PHYS%OUT%FHSCN, YDMF_PHYS%OUT%FHSSL, YDMF_PHYS%OUT%FHSSN, &
  & YDMF_PHYS%OUT%FHPCL, YDMF_PHYS%OUT%FHPCN, YDMF_PHYS%OUT%FHPSL, YDMF_PHYS%OUT%FHPSN, &
  & PFEPFP =YDMF_PHYS%OUT%FEPFP, PFCMPCQ=YDMF_PHYS%OUT%FCMPCQ, PFCMPSN=YDMF_PHYS%OUT%FCMPSN, &
  & PFCMPSL=YDMF_PHYS%OUT%FCMPSL, YDDDH=YDDDH )
  
  CALL CPUTQY(YDGEOMETRY%YRDIMV, YDGMV, YGFL, YDPTRSLB1, YDPHY, YDCPG_DIM%KLON, YDCPG_DIM%KIDIA, &
  & YDCPG_DIM%KFDIA, YDCPG_DIM%KFLEVG, PDTPHY, &
  & IPGFL, ISLB1T9, ISLB1U9, ISLB1V9, ISLB1VD9, ISLB1GFL9, ZTENDH, ZTENDT, YDMF_PHYS%OUT%TENDU,    &
  & YDMF_PHYS%OUT%TENDV, ZTENDU, ZTENDV, ZTENDD, ZTENDGFL, YLMF_PHYS_STATE%YCPG_DYN%RCP%CP,        &
  & YLMF_PHYS_STATE%YCPG_PHY%XYB%DELP, YLMF_PHYS_STATE%T, YLMF_PHYS_STATE%U, YLMF_PHYS_STATE%V,    &
  & PB1, PGMVT1, PGFLT1, YDMF_PHYS%OUT%FDIS)    
   
ELSE

  ! start ZTENDGFLR at 1 because it is dimensionned (:,:,0:n)
  CALL CPUTQY_AROME(YDMODEL,YDGEOMETRY%YRDIMV,YDGMV,YDCPG_DIM%KLON,YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA,&
     & YDCPG_DIM%KFLEVG,PDTPHY,IPGFL,IPTR,&
     & ZTENDT, ZTENDGFLR(:,:,1:), YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDD ,&
     & PB1, PGMVT1, PGFLT1)
ENDIF
  

!     ------------------------------------------------------------------ 
!*       5.    Final calculations.
!              -------------------

! * Restore the initial value of some pseudo-historical surface buffers
!   if relevant.
IF (LLDIAB) THEN
  IF (LL_SAVE_PHSURF) THEN
    CALL MF_PHYS_SAVE_PHSURF_PART2 (YDCPG_DIM, YDMF_PHYS_TMP, YDMF_PHYS_SURF, YDVARS, YDSURF, YDMODEL)
  ENDIF
ENDIF

!-------------------------------------------------
! Extract Single Column Model profiles from 3D run or 
! write LFA file for MUSC (1D model)
!-------------------------------------------------
IF(LGSCM.OR.LMUSCLFA) THEN
  IF (LAROME) THEN
    DO JLEV=1,YDCPG_DIM%KFLEVG
      DO JROF=YDCPG_DIM%KIDIA,YDCPG_DIM%KFDIA
        YDCPG_MISC%NEB(JROF,JLEV)=YDVARS%A%T1(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
  CALL WRITEPHYSIO(YDGEOMETRY, YDCPG_MISC, YDCPG_PHY0, YDCPG_PHY9, YDMF_PHYS, YDCPG_DYN0, YDCPG_DYN9, &
  & YDMF_PHYS_SURF, YDVARS, YDSURF, YDDPHY, YDRIP, YDMODEL%YRML_PHY_MF, YDCPG_DIM%KFDIA, YDCPG_DIM%KIDIA, YDCPG_DIM%KGL1, YDCPG_DIM%KGL2,        &
  & YDCPG_DIM%KSTGLO, YDCPG_DIM%KSTEP, NTSSG, YSP_SBD%      NLEVS, YDGSGEOM%GELAM, YDGSGEOM%GEMU, YDGSGEOM%GM, YDOROG%    &
  & OROG, YDGSGEOM%RCORI, YDCSGEOM%RATATH, YDCSGEOM%RATATX, YDGSGEOM%       GECLO, YDGSGEOM%GESLO, YDMF_PHYS_TMP  )
ENDIF

IF (LEDR) THEN
  YDMF_PHYS_SURF%GSD_DI%PXEDR(:,:)=YDMF_PHYS%OUT%EDR(:,:)
ENDIF

CALL MF_PHYS_PRECIPS (YDCPG_DIM, YDMF_PHYS_TMP, YDMF_PHYS_SURF, YDMODEL)

! Restore Tt and grad(Tt) for NHQE model.
IF (LNHQE) THEN
  CALL MF_PHYS_NHQE_PART2 (YDGEOMETRY, YDCPG_DIM, YDMF_PHYS_TMP, YDVARS)
ENDIF

!     ------------------------------------------------------------------

!       6. destructor for procset
IF (LINTFLEX) CALL CLEANINTPROCSET(YLPROCSET)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('MF_PHYS_MPA',1,ZHOOK_HANDLE)
END SUBROUTINE MF_PHYS_MPA
