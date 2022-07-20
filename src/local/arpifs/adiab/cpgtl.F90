#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE CPGTL(YDGEOMETRY,YDGMV,YDGMV5,YDSURF5,&
 !---------------------------------------------------------------------
 ! - INPUT .
 & YDMODEL,LDFSTEP,LDUSEPB1,LDSTR,&
 & KGPCOMP,KSTGLO,KIBL,&
 & PDT,PDTPHY,PTE,PBETADT,YDSL,&
 & PFORCEU,PFORCEV,PFORCET,PFORCEQ,&
 !- --------------------------------------------------------------------
 ! - INPUT/OUTPUT .
 & PGMV,PGMVS,PGFL,&
 !---------------------------------------------------------------------
 ! - OUTPUT .
 & KSETTLOFF,PB1,PB2,PGMVT1,PGMVT1S,PGFLT1,&
 !---------------------------------------------------------------------
 ! - TRAJECTORY .
 & PB15,PGMV5,PGMV5S,PGFL5,PTRAJ_PHYS,PTRAJ_PHYS_TLAD,&
 & PTRAJ_SRFC,PTRAJ_CST,PTRAJ_SLAG)

!**** *CPGTL* - Tan. lin. grid point calculations.

!     Purpose.
!     --------
!           Tan. lin. grid point calculations.

!**   Interface.
!     ----------
!        *CALL* *CPGTL(...)*

!        Explicit arguments :
!        --------------------

!      INPUT:
!      ------
!       LDFSTEP        : .T. if first step.
!       LDUSEPB1       : .T. if updating PB1
!       LDSTR          : .T. if stretched grid
!       KGPCOMP        : total number of grid points:
!                        ndglg*ndlon in global / ndguxg*ndlon in LAM
!                        used for model computations.
!       KSTGLO         : global offset.
!       KIBL           : index into YRGSGEOM/YRCSGEOM instances in YDGEOMETRY
!       PDT            : For a leap-frog scheme (three time level scheme):
!                        'dt' at the first time-step, '2 dt' otherwise.
!                        For a 2TL SL scheme: timestep 'dt'.
!       PDTPHY         : timestep used in the physics.
!       PTE            : 1. or 0. according to different configurations.
!       PBETADT        : BETADT or 0. according to different configurations.
!       YDSL           : SL_STRUCT definition

!      INPUT/OUTPUT:
!      -------------
!       PGMV           : upper air GMV variables at t and t-dt.
!       PGMVS          : surface GMV variables at t and t-dt.
!       PGFL           : unified_treatment grid-point fields at t.

!      OUTPUT:
!      -------
!       PB1            : "SLBUF1" buffer for interpolations in SL scheme.
!       PB2            : "SLBUF2" buffer.
!       PGMVT1         : upper air GMV variables buffer.
!       PGMVT1S        : surface GMV variables buffer.
!       PGFLT1         : GFL variables buffer.

!      TRAJECTORY:
!      -----------
!       PB15           : "SLBUF15" buffer.
!       PGMV5          : upper air GMV variables at t and t-dt.
!       PGMV5S         : surface GMV variables at t and t-dt.
!       PGFL5          : unified_treatment grid-point fields at t.

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
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 88-02-04

! Modifications
! -------------
!   17-Apr-2007 S. Ivatek-S Over dimensioning of ZGPNH and ZGPNH5 to NFGPNH+1,
!                          boundary checking problem if NFGPNH=0 bf
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2008): remove dummy CDLOCK + cleanings
!   F. Vana  13-Jan-2009: SLHDA and SLHDD0 for SLHD KAPPA
!   K. Yessad (Aug 2009): remove NTRSLTYPE/=2 cases
!   K. Yessad (Aug 2009): remove LPC_OLD in TL and AD codes.
!   K. Yessad (Nov 2009): cleanings, DT/Dt now pre-computed in CPG_GP_TL.
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures + cleanings.
!   O.Riviere (Feb 2011): add QL/QI in TL/AD.
!   G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!   T.Wilhelmsson (Apr 2012): Move OpenMP loop to cpg_drv_tl
!   F. Vana  28-Nov-2013 : Redesigned trajectory handling
!   M. Diamantakis (Feb 2014): add code for LSETTLSVF option
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (Dec 2016): Prune obsolete options.
!   E. Arbogast (Sep 2018): SURF comes from trajectory and is intent IN
! End Modifications
!     ------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMGMV   , ONLY : TGMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE INTDYN_MOD,ONLY : YYTTND, YYTHW5, YYTCTY5, YYTCTY0,&
 & YYTXYBDER5, YYTRCP5, YYTRCP95, YYTRCP0, YYTRCP9,&
 & YYTXYB0, YYTXYB5, YYTXYB9, YYTXYB95, YYTGMVT95, YYTGFLT95
USE EINT_MOD , ONLY : SL_STRUCT
USE YOMTRAJ  , ONLY : TRAJ_PHYS_TYPE, TRAJ_PHYS_TLAD_TYPE,&
 & TRAJ_SRFC_TYPE, TRAJ_CST_TYPE, TRAJ_SLAG_TYPE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)           ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)               ,INTENT(INOUT) :: YDGMV
TYPE(TGMV)               ,INTENT(INOUT) :: YDGMV5
TYPE(TSURF)              ,INTENT(IN)    :: YDSURF5 ! Meta data for trajectory surface fields
TYPE(MODEL)              ,INTENT(INOUT) :: YDMODEL
LOGICAL                  ,INTENT(IN)    :: LDFSTEP
LOGICAL                  ,INTENT(IN)    :: LDUSEPB1
LOGICAL                  ,INTENT(IN)    :: LDSTR
INTEGER(KIND=JPIM)       ,INTENT(IN)    :: KGPCOMP
INTEGER(KIND=JPIM)       ,INTENT(IN)    :: KSTGLO
INTEGER(KIND=JPIM)       ,INTENT(OUT)   :: KSETTLOFF(YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM)       ,INTENT(IN)    :: KIBL
REAL(KIND=JPRB)          ,INTENT(IN)    :: PDT
REAL(KIND=JPRB)          ,INTENT(IN)    :: PDTPHY
REAL(KIND=JPRB)          ,INTENT(IN)    :: PTE
REAL(KIND=JPRB)          ,INTENT(IN)    :: PBETADT
REAL(KIND=JPRB)          ,INTENT(IN)    :: PFORCEU(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)          ,INTENT(IN)    :: PFORCEV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)          ,INTENT(IN)    :: PFORCET(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)          ,INTENT(IN)    :: PFORCEQ(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
TYPE(SL_STRUCT)          ,INTENT(INOUT) :: YDSL
REAL(KIND=JPRB)          ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)          ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS)
REAL(KIND=JPRB)          ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB)          ,INTENT(OUT)   :: PB1(YDSL%NASLB1,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)          ,INTENT(OUT)   :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)          ,INTENT(OUT)   :: PGMVT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%YT1%NDIM)
REAL(KIND=JPRB)          ,INTENT(OUT)   :: PGMVT1S(YDGEOMETRY%YRDIM%NPROMA,YDGMV%YT1%NDIMS)
REAL(KIND=JPRB)          ,INTENT(OUT)   :: PGFLT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM1)
REAL(KIND=JPRB)          ,INTENT(OUT)   :: PB15(YDSL%NASLB1,YDMODEL%YRML_DYN%YRPTRSLB15%NFLDSLB15)
REAL(KIND=JPRB)          ,INTENT(INOUT) :: PGMV5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV5%YT5%NDIM)
REAL(KIND=JPRB)          ,INTENT(INOUT) :: PGMV5S(YDGEOMETRY%YRDIM%NPROMA,YDGMV5%YT5%NDIMS)
REAL(KIND=JPRB)          ,INTENT(INOUT) :: PGFL5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NDIM5)
TYPE(TRAJ_PHYS_TYPE)     ,INTENT(IN)    :: PTRAJ_PHYS
TYPE(TRAJ_PHYS_TLAD_TYPE),INTENT(INOUT) :: PTRAJ_PHYS_TLAD
TYPE(TRAJ_SRFC_TYPE)     ,INTENT(IN)    :: PTRAJ_SRFC
TYPE(TRAJ_CST_TYPE)      ,INTENT(IN)    :: PTRAJ_CST
TYPE(TRAJ_SLAG_TYPE)     ,INTENT(IN)    :: PTRAJ_SLAG
!     ------------------------------------------------------------------
! - derivatives of orography.
REAL(KIND=JPRB) :: ZOROGL(YDGEOMETRY%YRDIM%NPROMA)                          ! cf. POROGL in CPG5_GP.
REAL(KIND=JPRB) :: ZOROGM(YDGEOMETRY%YRDIM%NPROMA)                          ! cf. POROGM in CPG5_GP.
! - trajectory quantities at time t.
REAL(KIND=JPRB) :: ZPRE5(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)                  ! cf. PRE5 in CPG5_GP.
REAL(KIND=JPRB) :: ZPRE5F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                   ! cf. PRE5F in CPG5_GP.
REAL(KIND=JPRB) :: ZPRE5L(YDGEOMETRY%YRDIM%NPROMA)                          ! cf. PRE5L in CPG5_GP.
REAL(KIND=JPRB) :: ZPRE5M(YDGEOMETRY%YRDIM%NPROMA)                          ! cf. PRE5M in CPG5_GP.
REAL(KIND=JPRB) :: ZRPREF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                  ! cf. PRPREF5 in CPG5_GP.
REAL(KIND=JPRB) :: ZRT5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                     ! cf. PRT5 in CPG5_GP.
REAL(KIND=JPRB) :: ZRT5L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                    ! cf. PRT5L in CPG5_GP.
REAL(KIND=JPRB) :: ZRT5M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                    ! cf. PRT5M in CPG5_GP.
REAL(KIND=JPRB) :: ZXYB5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB5%NDIM)       ! cf. PXYB5 in CPG5_GP.
REAL(KIND=JPRB) :: ZRRED5(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)                  ! cf. PRRED5 in CPG5_GP.
REAL(KIND=JPRB) :: ZUVH5(YDGEOMETRY%YRDIM%NPROMNH,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTHW5%NDIM)     ! cf. PUVH5 in CPG5_GP.
REAL(KIND=JPRB) :: ZRCP5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP5%NDIM)       ! cf. PRCP5 in CPG5_GP.
REAL(KIND=JPRB) :: ZPHI5(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)                  ! cf. PHI5 in CPG5_GP.
REAL(KIND=JPRB) :: ZPHIF5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                   ! cf. PHIF5 in CPG5_GP.
REAL(KIND=JPRB) :: ZCTY5(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTCTY5%NDIM)     ! cf. PCTY5 in CPG5_GP.
REAL(KIND=JPRB) :: ZSD_VF5(YDGEOMETRY%YRDIM%NPROMA,YDSURF5%YSD_VFD%NDIM)            ! Surface traj (VARSF)
REAL(KIND=JPRB) :: ZQCHA5L(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)                 ! cf. PQCHA5L in CPG5_GP.
REAL(KIND=JPRB) :: ZQCHA5M(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)                 ! cf. PQCHA5M in CPG5_GP.
REAL(KIND=JPRB) :: ZXYBDER5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYBDER5%NDIM) ! cf. PXYBDER5 in CPG5_GP.
REAL(KIND=JPRB) :: ZRNHPPI5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                 ! cf. PRNHPPI5 in CPG5_GP.
REAL(KIND=JPRB) :: ZPHI5FL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                  ! cf. PHI5FL in CPG5_GP.
REAL(KIND=JPRB) :: ZPHI5FM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                  ! cf. PHI5FM in CPG5_GP.
REAL(KIND=JPRB) :: ZNH15L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                   ! cf. PNH15L in CPG5_GP.
REAL(KIND=JPRB) :: ZNH15M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                   ! cf. PNH15M in CPG5_GP.
! - trajectory quantities at time t-dt.
REAL(KIND=JPRB) :: ZPRE95(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)                 ! cf. PRE95 in CPG5_GP.
REAL(KIND=JPRB) :: ZPRE95F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                  ! cf. PRE95F in CPG5_GP.
REAL(KIND=JPRB) :: ZXYB95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB95%NDIM)     ! cf. PXYB95 in CPG5_GP.
REAL(KIND=JPRB) :: ZGFLT95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTGFLT95%NDIM)   ! cf. PGFLT95 in CPG5_GP.
REAL(KIND=JPRB) :: ZGMVT95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTGMVT95%NDIM)   ! cf. PGMVT95 in CPG5_GP.
REAL(KIND=JPRB) :: ZRCP95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP95%NDIM)     ! cf. PRCP95 in CPG5_GP.
REAL(KIND=JPRB) :: ZPHI95(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)                 ! cf. PHI95 in CPG5_GP.
REAL(KIND=JPRB) :: ZPHIF95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                  ! cf. PHIF95 in CPG5_GP.
REAL(KIND=JPRB) :: ZRRED95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                  ! cf. PRRED95 in CPG5_GP.
! - auxiliary quantities at time t.
REAL(KIND=JPRB) :: ZPRE0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)                  ! cf. PRE0 in CPG_GP_TL.
REAL(KIND=JPRB) :: ZPRE0F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                   ! cf. PRE0F in CPG_GP_TL.
REAL(KIND=JPRB) :: ZPRE0L(YDGEOMETRY%YRDIM%NPROMA)                          ! cf. PRE0L in CPG_GP_TL.
REAL(KIND=JPRB) :: ZPRE0M(YDGEOMETRY%YRDIM%NPROMA)                          ! cf. PRE0M in CPG_GP_TL.
REAL(KIND=JPRB) :: ZRT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                    ! cf. PRT0L in CPG_GP_TL.
REAL(KIND=JPRB) :: ZRT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                    ! cf. PRT0M in CPG_GP_TL.
REAL(KIND=JPRB) :: ZXYB0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB0%NDIM)       ! cf. PXYB0 in CPG_GP_TL.
!NYC-NH REAL(KIND=JPRB) :: ZUVH0(NPROMNH,0:NFLEVG,YYTHW0%NDIM) ! cf. PUVH0 in CPG_GP_TL.
REAL(KIND=JPRB) :: ZRCP0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP0%NDIM)       ! cf. PRCP0 in CPG_GP_TL.
REAL(KIND=JPRB) :: ZPHI0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)                  ! cf. PHI0 in CPG_GP_TL.
REAL(KIND=JPRB) :: ZPHIF0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                   ! cf. PHIF0 in CPG_GP_TL.
REAL(KIND=JPRB) :: ZCTY0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG,YYTCTY0%NDIM)     ! cf. PCTY0 in CPG_GP_TL.
REAL(KIND=JPRB) :: ZKENE0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                   ! cf. PKENE0 in CPG_GP_TL.
! - adiabatic Lagrangian tendencies.
REAL(KIND=JPRB) :: ZATND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTTND%NDIM)
! - quantities for physics.
REAL(KIND=JPRB) :: ZQS(YDGEOMETRY%YRDIM%NPROMM)                             ! cf. PQS in MF_PHYS_TL.
REAL(KIND=JPRB) :: ZQS1(YDGEOMETRY%YRDIM%NPROMM)                            ! cf. PQS1 in MF_PHYS_TL.
REAL(KIND=JPRB) :: ZTS(YDGEOMETRY%YRDIM%NPROMM)                             ! cf. PTS in MF_PHYS_TL.
REAL(KIND=JPRB) :: ZTS1(YDGEOMETRY%YRDIM%NPROMM)                            ! cf. PTS1 in MF_PHYS_TL.
REAL(KIND=JPRB) :: ZSNS(YDGEOMETRY%YRDIM%NPROMM)                            ! cf. PSNS in MF_PHYS_TL.
! - auxiliary quantities at time t-dt.
REAL(KIND=JPRB) :: ZPRE9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)                  ! cf. PRE9 in CPG_GP_TL.
REAL(KIND=JPRB) :: ZPRE9F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                   ! cf. PRE9F in CPG_GP_TL.
REAL(KIND=JPRB) :: ZXYB9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB9%NDIM)       ! cf. PXYB9 in CPG_GP_TL.
REAL(KIND=JPRB) :: ZRCP9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTRCP9%NDIM)       ! cf. PCP9 in CPG_GP_TL.
REAL(KIND=JPRB) :: ZPHI9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)                  ! cf. PHI9 in CPG_GP_TL.
REAL(KIND=JPRB) :: ZPHIF9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)                   ! cf. PHIF9 in CPG_GP_TL.
! - other quantities.
REAL(KIND=JPRB) :: ZGPPC(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRGPPC%NFGPPC+1)                  ! cf. PGPPC in CPG_GP.
REAL(KIND=JPRB) :: ZGPPC5(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRGPPC%NFGPPC+1)                 ! cf. PGPPC5 in CPG5_GP.
REAL(KIND=JPRB) :: ZSLBUF1AU(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB1%NFLDSLB1)              ! local version of SLBUF1 
REAL(KIND=JPRB) :: ZSLBUF1AU5(YDGEOMETRY%YRDIM%NPROMA,YDMODEL%YRML_DYN%YRPTRSLB15%NFLDSLB15)            ! local version of SLBUF15 

REAL(KIND=JPRB) :: ZDUM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDUMSPT(YDGEOMETRY%YRDIM%NPROMA)

INTEGER(KIND=JPIM) :: JROF,JFLD
INTEGER(KIND=JPIM) :: IEND, IENDC, IST, ISTC
INTEGER(KIND=JPIM) :: IBL
INTEGER(KIND=JPIM) :: ICEND, ICENDE
INTEGER(KIND=JPIM) :: IFLAG
INTEGER(KIND=JPIM) :: IOFF

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "cpg5_gp.intfb.h"
#include "cpg_dyn_tl.intfb.h"
#include "cpg_end_tl.intfb.h"
#include "cpg_gp_tl.intfb.h"
#include "gpmpfc5.intfb.h"
#include "mf_phystl.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CPGTL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 &  YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDDYN=>YDMODEL%YRML_DYN%YRDYN,YDPTRSLB1=>YDMODEL%YRML_DYN%YRPTRSLB1, &
 & YDPTRSLB2=>YDMODEL%YRML_DYN%YRPTRSLB2, &
 & YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL,YDPTRSLB15=>YDMODEL%YRML_DYN%YRPTRSLB15,YGFL=>YDMODEL%YRML_GCONF%YGFL, &
 & YDPTRGPPC=>YDMODEL%YRML_DYN%YRPTRGPPC)

ASSOCIATE(NDIM=>YGFL%NDIM, NDIM1=>YGFL%NDIM1, NDIM5=>YGFL%NDIM5, &
 & NPROMA=>YDDIM%NPROMA, NPROMM=>YDDIM%NPROMM, NPROMNH=>YDDIM%NPROMNH, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, SLHDA=>YDDYN%SLHDA, SLHDD0=>YDDYN%SLHDD0, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & GM=>YDGSGEOM%GM, NGPLAT=>YDGSGEOM%NGPLAT, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0, &
 & YT1=>YDGMV%YT1, YT9=>YDGMV%YT9, &
 & YT5=>YDGMV5%YT5, LSIMPH=>YDSIMPHL%LSIMPH, &
 & NFGPPC=>YDPTRGPPC%NFGPPC, &
 & NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & NFLDSLB15=>YDPTRSLB15%NFLDSLB15, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2)
!     ------------------------------------------------------------------

!*       1.    Preliminary calculations.
!              -------------------------

ZDUM(:,:)=0.0_JPRB
ZDUMSPT(:)=0.0_JPRB
ZGPPC(:,:)=0.0_JPRB
ZGPPC5(:,:)=0.0_JPRB

IST=1
ISTC=1

IOFF=KSTGLO
IBL=(KSTGLO-1)/NPROMA+1
ICEND=MIN(NPROMA,KGPCOMP-KSTGLO+1)
ICENDE=MIN(NPROMA,NGPTOT-KSTGLO+1)
IEND=ICEND
IENDC=ICENDE

!     ------------------------------------------------------------------

!*       3.    Gridpoint computations for the trajectory
!              -----------------------------------------

CALL CPG5_GP(YDGEOMETRY,YDGMV5,YDSURF5,YDMODEL%YRML_GCONF,YDMODEL%YRML_PHY_SLIN%YRPHLC,YDSIMPHL,YDMODEL%YRML_DYN,&
 !---------------------------------------------------------------------
 ! - INPUT .
 & IST,IEND,KSTGLO,&
 & KIBL,&
 & PTRAJ_PHYS,PTRAJ_SRFC,PTRAJ_CST,&
 !---------------------------------------------------------------------
 ! - INPUT/OUTPUT .
 & PGFL5,PGMV5,PGMV5S,&
 !---------------------------------------------------------------------
 ! - OUTPUT .
 & ZOROGL,ZOROGM,&
 & ZGMVT95,ZGFLT95,&
 & ZPRE5,ZPRE5L,ZPRE5M,ZPRE5F,ZXYB5,ZRPREF5,&
 & ZUVH5,ZPHI5,ZPHIF5,ZRCP5,ZRRED5,ZCTY5,ZRT5,ZRT5L,ZRT5M,&
 & ZPRE95,ZPRE95F,ZXYB95,ZRRED95,ZPHI95,ZPHIF95,ZRCP95,&
 & ZGPPC5,ZQS,ZQS1,ZTS,ZTS1,ZSNS,ZSD_VF5,ZQCHA5L,ZQCHA5M,&
 & ZXYBDER5,ZRNHPPI5,ZPHI5FL,ZPHI5FM,ZNH15L,ZNH15M)

!     ------------------------------------------------------------------

!*       4.    Gridpoint computations before physics for perturbations
!              at t-dt and t
!              -------------

CALL CPG_GP_TL(YDGEOMETRY,YDGMV,YDGMV5,&
 !---------------------------------------------------------------------
 ! - INPUT .
 & YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_GCONF,YDMODEL%YRML_DYN,YDSIMPHL,IST,IEND,LDFSTEP,PTE,KIBL,&
 & PGFL,&
 !---------------------------------------------------------------------
 ! - INPUT/OUTPUT .
 & PGMV,PGMVS,&
 !---------------------------------------------------------------------
 ! - OUTPUT .
 & ZPRE0,ZPRE0L,ZPRE0M,ZPRE0F,ZXYB0,&
 !NYC-NH & ZUVH0, &
 & ZPHI0,ZPHIF0,ZRCP0,ZCTY0,ZKENE0,ZRT0L,ZRT0M,&
 & ZPRE9,ZPRE9F,ZXYB9,ZPHI9,ZPHIF9,ZRCP9,&
 & ZGPPC,PB2,&
 & PGMVT1,PGMVT1S,PGFLT1,&
 & ZATND,&
 !---------------------------------------------------------------------
 ! - TRAJECTORY (INPUT) .
 & PGFL5,PGMV5,PGMV5S,&
 & ZPRE5,ZPRE5L,ZPRE5M,ZPRE5F,ZXYB5,ZRCP5,&
 !NYC-NH & ZRRED5, &
 & ZCTY5,ZRT5,ZRT5L,ZRT5M,&
 !NYC-NH & ZUVH5, &
 & ZRPREF5,&
 !NYC-NH & ZGPPC5,&
 !NYC-NH & ZQCHA5L,ZQCHA5M, &
 & ZXYBDER5,ZRNHPPI5,ZPHI5FL,ZPHI5FM,&
 !NYC-NH & ZNH15L,ZNH15M, &
 & ZGMVT95,ZPRE95,ZXYB95,ZRCP95&
 !NYC-NH & ZRRED95 &
 & )

!     ------------------------------------------------------------------

!*       5.    Call physics.
!              -------------

! * Set ZSLBUF1AU and ZSLBUF1AU5 to zero.
IF (LDUSEPB1) THEN
  DO JFLD=1,NFLDSLB1
    DO JROF=KSTGLO,MIN(KSTGLO-1+NPROMA,KGPCOMP)
      ZSLBUF1AU(JROF-KSTGLO+1,JFLD)=0.0_JPRB
    ENDDO
  ENDDO
  DO JFLD=1,NFLDSLB15
    DO JROF=KSTGLO,MIN(KSTGLO-1+NPROMA,KGPCOMP)
      ZSLBUF1AU5(JROF-KSTGLO+1,JFLD)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF

! * Unlagged ECMWF physics (CALL EC_PHYS_TL) or
!   split ECMWF physics (CALL EC_PHYS_LSLPHY_TL) not yet coded.

! * Unlagged MF physics.
!   Remark! the management of physics for LPC_FULL PC scheme
!   (CALL CPG_PT_TL after CALL MF_PHYSTL) is not yet coded.

IF (LSIMPH.AND.(NCURRENT_ITER == 0)) THEN

  CALL MF_PHYSTL(YDGEOMETRY,YDGMV,YDGMV5,YDSURF5,&
   & YDMODEL,IST,IEND,KSTGLO,PDTPHY,&
   & KIBL,&
   & PGMV,PGFL,&
   & ZRCP0,ZXYB0,ZPHI0,ZPHIF0,ZPRE0,ZPRE0F,ZCTY0,&
   & ZRCP9,ZXYB9,ZPHI9,ZPHIF9,ZPRE9,ZPRE9F,&
   & ZQS,ZQS1,ZSNS,ZTS,ZTS1,&
   & PFORCEU,PFORCEV,PFORCET,PFORCEQ, &
   & ZSLBUF1AU,PGMVT1,PGFLT1,&
   & PGMV5,PGFL5,&
   & ZRCP5,ZXYB5,ZPHI5,ZPHIF5,ZPRE5,ZPRE5F,ZCTY5,&
   & ZGMVT95,ZGFLT95,ZRCP95,ZXYB95,ZPHI95,ZPHIF95,ZPRE95,ZPRE95F,ZSD_VF5,&
   & PTRAJ_PHYS,PTRAJ_PHYS_TLAD)

ENDIF

!     ------------------------------------------------------------------

!*       6.    Dynamics.
!              ---------

CALL CPG_DYN_TL(YDGEOMETRY,YDGMV,YDGMV5,YDSURF5,&
 ! --- INPUT ----------------------------------------------------------
 & YDMODEL,IST,IEND,KSTGLO,PBETADT,PDT,&
 & SLHDA(:,IBL),SLHDD0(:,IBL),&
 & KIBL,ZOROGL,ZOROGM,&
 & ZPRE0,ZXYB0(1,1,YYTXYB0%M_RDELP),ZPHI0,ZPHIF0,ZCTY0,ZKENE0,ZATND,&
 & PGFL,&
 ! --- INPUT/OUTPUT ---------------------------------------------------
 & KSETTLOFF,PGMV,PGMVS,&
 & ZSLBUF1AU,PB2,PGMVT1,PGMVT1S,PGFLT1,&
 ! --- TRAJECTORY-INPUT -----------------------------------------------
 & ZPRE5,ZXYB5(1,1,YYTXYB5%M_RDELP),ZPHI5,ZCTY5,ZSD_VF5,&
 ! --- TRAJECTORY-INOUT -----------------------------------------------
 & ZSLBUF1AU5,PGMV5,PGMV5S,PGFL5,PTRAJ_SLAG)

! * Update PB1 and PB15.
IF (LDUSEPB1) THEN
  DO JFLD=1,NFLDSLB1
    DO JROF=KSTGLO,MIN(KSTGLO-1+NPROMA,KGPCOMP)
      PB1(YDSL%NSLCORE(JROF),JFLD)=ZSLBUF1AU(JROF-KSTGLO+1,JFLD)
    ENDDO
  ENDDO
  DO JFLD=1,NFLDSLB15
    DO JROF=KSTGLO,MIN(KSTGLO-1+NPROMA,KGPCOMP)
      PB15(YDSL%NSLCORE(JROF),JFLD)=ZSLBUF1AU5(JROF-KSTGLO+1,JFLD)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!*       7.    Final part of vector gridpoint operations.
!              ------------------------------------------

CALL CPG_END_TL(YDGEOMETRY,YDGMV,&
 !---------------------------------------------------------------------
 ! - INPUT .
 & YDMODEL%YRML_GCONF,YDDYN,IST,IEND,GM,&
 !---------------------------------------------------------------------
 ! - INPUT/OUTPUT .
 & PGFL,PGMV,PGMVS)

!     ------------------------------------------------------------------

!*       8.    Division by M or M**2 for derivatives: trajectory
!              -------------------------------------------------

IF(LDSTR) THEN
  IFLAG=1
  CALL GPMPFC5(YDGMV5,YDMODEL%YRML_GCONF,NPROMA,NFLEVG,IST,IEND,IFLAG,GM,PGMV5,PGMV5S,PGFL5)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPGTL',1,ZHOOK_HANDLE)
END SUBROUTINE CPGTL
