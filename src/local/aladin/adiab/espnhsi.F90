SUBROUTINE ESPNHSI(&
 ! --- INPUT -----------------------------------------------------------------
 & YDCST, YDGEOMETRY,YDLDDH,YDRIP,YDDYN,YDEDYN,KM,KMLOC,KSTA,KEND,LDONEM,&
 ! --- INOUT -----------------------------------------------------------------
 & PSPVORG,PSPDIVG,PSPTG,PSPSPG,PSPSPDG,PSPSVDG,&
 & PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,PSPTNDSI_SPDG,PSPTNDSI_SVDG)

!**** *ESPNHSI* - SPECTRAL SPACE COMPUTATIONS FOR ALADIN

! !! The structure of ESPNHSI must remain very close to SPNHSI !!
!    The only differences between ESPNHSI and SPNHSI are the following ones:
!    * use of YDLEP%RLEPDIM instead of RLAPDI.
!    * use of YDLEP%RLEPINM instead of RLAPIN.
!    * for geometries with significant variations of the mapping factor M,
!      use the ALADIN designed option LESIDG instead of LSIDG.
!    * LIMPF not coded in ALADIN for the time being.

!     Purpose.
!     --------

!     Semi-implicit scheme in the NH model with (Qcha, dver or d4) as
!     NH prognostic variables, in ALADIN.

!     When the constraint C1 is matched, I_NITERHELM=0; if not matched
!     I_NITERHELM>0.
!     The predictor has dver(t+dt) (or d4(t+dt)), D'(t+dt), zeta'(t+dt),
!     Qcha(t+dt), T(t+dt), log(prehyds)(t+dt) as unknown.
!     The corrector steps (JITER>0) work in an incremental manner,
!     the unknowns are (delta X) = X(jiter) - X(iter=0), where
!     X=dver(t+dt), D'(t+dt), etc...
!     The incremental method makes the calculation of the RHS simpler
!     in the corrector steps; this is equivalent to replace:
!      * zeta_star by 0.
!      * Dcha_star_star by 0.
!      * Dprim_star_star by
!        beta**2 (Delta t)**2 vnabla'**2 * C**2 * COR * Mbar**2 D'(t+dt,jiter-1)

!**   Interface.
!     ----------
!        *CALL* *ESPNHSI(...)

!        Explicit arguments :
!        --------------------
!        * INPUT:
!          KM          - Zonal wavenumber
!          KMLOC       - Zonal wavenumber (DM-local numbering)
!          KSTA        - First column processed
!          KEND        - Last column processed
!          LDONEM      - T if only one m if processed
!        * INPUT/OUTPUT:
!          PSPVORG     - Vorticity columns
!          PSPDIVG     - Divergence columns
!          PSPTG       - Temperature columns
!          PSPSPG      - Surface Pressure
!          PSPSPDG     - Pressure departure variable columns
!          PSPSVDG     - Vertical divergence variable columns

!     Method.
!     -------
!       There is a great description in ESPCM

!     Externals.
!     ----------
!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Sylvie Malardel (original, March 1992)
!        Karim YESSAD (Feb 2005, moving the NH SI scheme part in ESPNHSI)

!     Modifications.
!     --------------
!        K. Yessad Aug-2007: update comments and LIMPF markers according
!            to what is done in SPNHSI.
!        K. Yessad Aug-2007: iterative algorithm to relax the C1 constraint
!        A. Bogatchev Jun-2008: phasing cy34
!        I. Santos, I. Martinez Mar2010: variable map factor for RTM.
!        R. El Khatib and F. Voitus (Mar 2011): optimisations.
!        K. Yessad (Feb 2012): tests on CDCONF and LL3D in the caller, simplifications.
!        A.Bogatchev 11-04-2013 phasing cy40, coherence with modified modules
!                               and renamed namelists and functions
!        B. Bochenek (Apr 2015): Phasing: move some variables.
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!        K. Yessad (Dec 2016): Prune obsolete options.
!        K. Yessad (June 2017): Vertical-dependent SITRA.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST   , ONLY : TCST
USE YOMMP0   , ONLY : MYSETV
USE YOMDYNA  , ONLY : YRDYNA
USE YOMDYN   , ONLY : TDYN
USE YEMDYN   , ONLY : TEDYN
USE YOMLDDH  , ONLY : TLDDH
USE YOMRIP   , ONLY : TRIP
USE YOMCVER      , ONLY : LVERTFE, LVFE_LAPL_BC

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY), INTENT(IN)       :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT)    :: YDDYN
TYPE(TEDYN)    ,INTENT(INOUT)    :: YDEDYN
TYPE(TLDDH)    ,INTENT(INOUT)    :: YDLDDH
TYPE(TRIP)     ,INTENT(INOUT)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KM
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOC
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
LOGICAL           ,INTENT(IN)    :: LDONEM 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVORG(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIVG(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTG(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPG(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPDG(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSVDG(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_VORG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_DIVG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_TG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_SPDG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_SVDG(:,:)

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZZSPVORG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZZSPDIVG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZZSPTG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZZSPSPG(KSTA:KEND)
REAL(KIND=JPRB) :: ZZSPSPDG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZZSPSVDG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)

! ZZSPGDIVG to ZZSPSVD0G are useful only if NITERHELM>0
REAL(KIND=JPRB) :: ZZSPGDIVG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZZSPDIV0G(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZZSPGDIV0G(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZZSPVOR0G(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZZSPSVD0G(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)

REAL(KIND=JPRB) :: ZSDIV(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZST(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSP(KSTA:KEND)
REAL(KIND=JPRB) :: ZSDIVP(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSPDIVP(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSVED(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZR4D(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZR3D(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZR2D(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZR1D(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSRHS(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZWORK(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSRHS2(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSNHP(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)

REAL(KIND=JPRB), ALLOCATABLE :: ZSDIVPL(:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPDIVPL(:,:,:)

REAL(KIND=JPRB) :: ZPD(0:YDGEOMETRY%YRDIM%NSMAX), ZPE(0:YDGEOMETRY%YRDIM%NSMAX), ZPF(0:YDGEOMETRY%YRDIM%NSMAX)

INTEGER(KIND=JPIM) :: IM, IN, IOFF, ISP, ISPCOL, JLEV, JSP,&
 & IS0 ,IS02 ,ISE ,JN ,II
INTEGER(KIND=JPIM) :: I_NITERHELM,JITER

REAL(KIND=JPRB) :: ZBDT, ZBDT2, ZRRT, ZCC, ZBTCM, ZBTCMH
                                                                                
LOGICAL :: LLNH_NOC1

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "mxmaop.h"

#include "mxptma.h"
#include "mxture.h"
#include "mxturs.h"

#include "abor1.intfb.h"
#include "siseve.intfb.h"
#include "si_cccor.intfb.h"
#include "sidd.intfb.h"
#include "sigam.intfb.h"
#include "siptp.intfb.h"
#include "sitnu.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ESPNHSI',0,ZHOOK_HANDLE)
ASSOCIATE(YDLAP=>YDGEOMETRY%YRLAP, YDLEP=>YDGEOMETRY%YRELAP)
ASSOCIATE( &
 & SIMO=>YDDYN%SIMO, &
 & SIMI=>YDDYN%SIMI, &
 & SIFACI=>YDDYN%SIFACI, &
 & SITR=>YDDYN%SITR, &
 & SIVP=>YDDYN%SIVP, &
 & TDT=>YDRIP%TDT, &
 & RBTS2=>YDDYN%RBTS2, &
 & NITERHELM=>YDDYN%NITERHELM, &
 & LIMPF=>YDDYN%LIMPF, &
 & SIHEG=>YDDYN%SIHEG, &
 & SIHEG2=>YDDYN%SIHEG2, &
 & SIHEGB=>YDDYN%SIHEGB, &
 & SIHEGB2=>YDDYN%SIHEGB2, &
 & NISNAX=>YDGEOMETRY%YREDIM%NISNAX, &
 & NSMAX=>YDGEOMETRY%YRDIM%NSMAX, &
 & NPTRSV=>YDGEOMETRY%YRMP%NPTRSV, NPTRSVF=>YDGEOMETRY%YRMP%NPTRSVF, &
 & RSTRET=>YDGEOMETRY%YRGEM%RSTRET, &
 & LESIDG=>YDEDYN%LESIDG, &
 & NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & LRSIDDH=>YDLDDH%LRSIDDH, &
 & ESCGMAP=>YDGEOMETRY%YSPGEOM%ESCGMAP)
!     ------------------------------------------------------------------

! * Case where the C1 constraint is not matched.
LLNH_NOC1=.NOT.(YRDYNA%NDLNPR == 1)
IF (LLNH_NOC1) THEN
  I_NITERHELM=NITERHELM
ELSE
  I_NITERHELM=0
ENDIF

! LESIDG: like LSIDG but for ALADIN in the case where:
!  - large domain where the variations of the mapping factor M are large.
!  - type of projections where M is quasi constant on a latitude and
!    M can be written as a low order polynomial of mu=sin(latitude).

IF (LIMPF) CALL ABOR1(' ESPNHSI: LIMPF=T not yet available in ALADIN')

!     ------------------------------------------------------------------

!*       1.    MEMORY TRANSFER.
!              ----------------

IF (LIMPF) ZZSPVORG(1:NFLEVG,KSTA:KEND)=PSPVORG(1:NFLEVG,KSTA:KEND)
ZZSPDIVG(1:NFLEVG,KSTA:KEND)=PSPDIVG(1:NFLEVG,KSTA:KEND)
ZZSPTG  (1:NFLEVG,KSTA:KEND)=PSPTG  (1:NFLEVG,KSTA:KEND)
ZZSPSPG (         KSTA:KEND)=PSPSPG (         KSTA:KEND)
ZZSPSPDG(1:NFLEVG,KSTA:KEND)=PSPSPDG(1:NFLEVG,KSTA:KEND)
ZZSPSVDG(1:NFLEVG,KSTA:KEND)=PSPSVDG(1:NFLEVG,KSTA:KEND)

IF (LRSIDDH) THEN
  ! DDH memory transfer.
  IF (LIMPF) PSPTNDSI_VORG(1:NFLEVG,KSTA:KEND)=-PSPVORG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_DIVG(1:NFLEVG,KSTA:KEND)=-PSPDIVG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_TG  (1:NFLEVG,KSTA:KEND)=-PSPTG  (1:NFLEVG,KSTA:KEND)
  PSPTNDSI_SPDG(1:NFLEVG,KSTA:KEND)=-PSPSPDG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_SVDG(1:NFLEVG,KSTA:KEND)=-PSPSVDG(1:NFLEVG,KSTA:KEND)
  !the case of surface pressure has not been treated yet
ENDIF

!     ------------------------------------------------------------------

!*       2.    SEMI-IMPLICIT SPECTRAL COMPUTATIONS.
!              ------------------------------------

!*        2.1  Preliminary initialisations.

IF (LDONEM) THEN
  IOFF=NPTRSVF(MYSETV)-1
ELSE
  IOFF=NPTRSV(MYSETV)-1
ENDIF
ISPCOL=KEND-KSTA+1

ZBDT=RBTS2*TDT
ZBDT2=(ZBDT*RSTRET)**2

ZRRT =  1.0_JPRB/(YDCST%RD*SITR)
ZCC  =  YDCST%RCPD/YDCST%RCVD
ZBTCM=  ZBDT*ZBDT*YDCST%RD*SITR*ZCC*RSTRET*RSTRET
ZBTCMH= ZBTCM*YDCST%RG*YDCST%RG*ZRRT*ZRRT

IF (LESIDG) THEN
  ! ky: incorrect use of NSE0L; use the ALADIN counterpart NESE0L instead.
  !     (its set-up and allocation remain to be coded).
  IS0=YDLAP%NSE0L(KMLOC)
  IS02=0

  IF (KM == 0) THEN
    II = 2
  ELSE
    II = 4
  ENDIF

  ! Allocations
  ALLOCATE(ZSDIVPL(NFLEVG,0:NISNAX(KM),4))
  ALLOCATE(ZSPDIVPL(NFLEVG,0:NISNAX(KM),4))

  DO JN=0,NSMAX
    ZPD(JN)=ESCGMAP(1)
    ZPE(JN)=ESCGMAP(2)
    ZPF(JN)=ESCGMAP(3)
  ENDDO
ENDIF

! ky: if LESIDG, missing there the norther shift to have C+I and C+I+E
!     centered on the same reference latitude.

!*        2.2  Set up helper arrays for implicit Coriolis case.

! IF (LIMPF) THEN
!   LIMPF not coded for ALADIN-NH.
! ENDIF

DO JITER=0,I_NITERHELM

  !*        2.3  Computes right-hand side of Helmholtz equation.

  IF (JITER == 0) THEN

    ! * Provides:
    !   a/ - [ gamma [SITR Qcha_rhs - T_rhs] - (Rd SITR)/(SIPR) logprehyd_rhs ]
    !      in array ZSDIV.
    !   b/ (g TRSI)/(H TARSI) LLstar Qcha_rhs in array ZSVED.
    CALL SIDD(LVERTFE, YRDYNA, YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZSDIV,ZSVED,ZZSPSPDG,ZZSPTG,ZZSPSPG,ISPCOL)

    ! * Provides Dprim_star_star in ZSDIV.
    DO JLEV=1,NFLEVG
!cdir nodep
      DO JSP=KSTA,KEND
        IN=YDLAP%NVALUE(JSP+IOFF)
        IM=YDLEP%MVALUE(JSP+IOFF)
        ISP=YDLEP%NPME(IM)+IN
        ZSDIV(JLEV,JSP)=ZZSPDIVG(JLEV,JSP)-ZBDT*YDLEP%RLEPDIM(ISP)*ZSDIV(JLEV,JSP)
      ENDDO
    ENDDO

  ELSE

    ! * Compute C**2 * COR * Mbar**2 D'(t+dt,jiter-1)
    CALL SI_CCCOR(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ISPCOL,ZZSPGDIVG,ZSDIV)

    ! * Multiply by beta**2 (Delta t)**2 vnabla'**2
    DO JLEV=1,NFLEVG
!cdir nodep
      DO JSP=KSTA,KEND
        IN=YDLAP%NVALUE(JSP+IOFF)
        IM=YDLEP%MVALUE(JSP+IOFF)
        ISP=YDLEP%NPME(IM)+IN
        ZSDIV(JLEV,JSP)=ZBDT*ZBDT*YDLEP%RLEPDIM(ISP)*ZSDIV(JLEV,JSP)
      ENDDO
    ENDDO

  ENDIF

  ! * Provides Dprim_star_star_star (in ZSDIV) from Dprim_star_star (in ZSDIV)
  ! IF (LIMPF) THEN
  !   LIMPF not coded for ALADIN-NH.
  ! ENDIF

  ! * Provides Dprim_star_star in ZR1D.
  !   and Dcha_star_star in ZR2D.
  IF (JITER == 0) THEN
    DO JLEV=1,NFLEVG
!cdir nodep
      DO JSP=KSTA,KEND
        ZR2D(JLEV,JSP)=ZZSPSVDG(JLEV,JSP)-ZBDT*ZSVED(JLEV,JSP)
      ENDDO
    ENDDO
  ELSE
    ! Dcha_star_star is replaced by zero.
    ZR2D(1:NFLEVG,KSTA:KEND)=0.0_JPRB
  ENDIF
  ZR1D(1:NFLEVG,KSTA:KEND)=ZSDIV(1:NFLEVG,KSTA:KEND)

  IF (LESIDG) THEN
    ! Multiply ZR1D by M**2 (CALL MXPTMA using SCGMAP, see SPCSI),
    !  then divide it by RSTRET**2.
    ZSDIVPL(:,:,:)=0.0_JPRB
    ZSPDIVPL(:,:,:)=0.0_JPRB
    DO JN=0,NISNAX(KM)
      ISE=KSTA+4*JN
      ZSDIVPL(:,JN,1:4)=ZR1D(:,ISE:ISE+3)
    ENDDO
    CALL MXPTMA(NISNAX(KM)+1,NFLEVG,NFLEVG,II,&
     & ZPD,ZPE,ZPF,ZPE,ZPF,&
     & ZSDIVPL,ZSPDIVPL)

    DO JN=0,NISNAX(KM)
      ISE=KSTA+4*JN
      ZR1D(:,ISE:ISE+3)=ZSPDIVPL(:,JN,1:4)/(RSTRET*RSTRET)
    ENDDO
  ENDIF

  ! * Provides "LLstar Dprim_star_star" (mult by a constant coefficient)
  !   in ZSRHS.
  CALL SISEVE(LVERTFE, LVFE_LAPL_BC, YRDYNA,YDGEOMETRY,YDDYN,1,NFLEVG,ZR1D,ZSRHS,ISPCOL)

  ! * Provides "Tau Dprim_star_star" (mult by a constant coefficient) in ZST.
  CALL SITNU(LVERTFE, YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZR1D,ZST,ZSP,ISPCOL)

  ! * Provides "LLstar Tau Dprim_star_star" (mult by a constant coefficient)
  !   in ZWORK.
  CALL SISEVE(LVERTFE, LVFE_LAPL_BC, YRDYNA,YDGEOMETRY,YDDYN,1,NFLEVG,ZST,ZWORK,ISPCOL)
  DO JLEV=1,NFLEVG
!cdir nodep
    DO JSP=KSTA,KEND
      ZWORK(JLEV,JSP)=ZWORK(JLEV,JSP)*YDCST%RCVD*ZRRT
    ENDDO
  ENDDO

  ! * Provides
  !   "[I - beta**2 (Delta t)**2 C**2 Mbar**2 vnabla'**2] Dcha_star_star"
  !   in ZSVED.
  IF (JITER == 0) THEN
    IF (LESIDG) THEN
      ! First compute ZR3D=ZR2D*YDLEP%RLEPDIM.
      ! Multiply ZR3D by M**2 (CALL MXPTMA using SCGMAP, see SPCSI),
      !  put the result in a separate array ZR4D.
      ! Compute ZSVED using ZR2D and ZR4D:
      !  ZSVED=ZR2D-(ZBTCM/(RSTRET*RSTRET))*ZR4D.
      DO JLEV=1,NFLEVG
        DO JSP=KSTA,KEND
          IN=YDLAP%NVALUE(JSP+IOFF)
          IM=YDLEP%MVALUE(JSP+IOFF)
          ISP=YDLEP%NPME(IM)+IN
          ZR3D(JLEV,JSP)=YDLEP%RLEPDIM(ISP)*ZR2D(JLEV,JSP)
        ENDDO
      ENDDO
      ZSDIVPL(:,:,:)=0.0_JPRB
      ZSPDIVPL(:,:,:)=0.0_JPRB
      DO JN=0,NISNAX(KM)
        ISE=KSTA+4*JN
        ZSDIVPL(:,JN,1:4)=ZR3D(:,ISE:ISE+3)
      ENDDO
      CALL MXPTMA(NISNAX(KM)+1,NFLEVG,NFLEVG,II,&
       & ZPD,ZPE,ZPF,ZPE,ZPF,&
       & ZSDIVPL,ZSPDIVPL)
      DO JN=0,NISNAX(KM)
        ISE=KSTA+4*JN
        ZR4D(:,ISE:ISE+3)=ZSPDIVPL(:,JN,1:4)
      ENDDO
      DO JLEV=1,NFLEVG
        DO JSP=KSTA,KEND
          ZSVED(JLEV,JSP)=ZR2D(JLEV,JSP)-(ZBTCM/(RSTRET*RSTRET))*ZR4D(JLEV,JSP)
        ENDDO
      ENDDO
    ELSE
      DO JLEV=1,NFLEVG
!cdir nodep
        DO JSP=KSTA,KEND
          IN=YDLAP%NVALUE(JSP+IOFF)
          IM=YDLEP%MVALUE(JSP+IOFF)
          ISP=YDLEP%NPME(IM)+IN
          ZSVED(JLEV,JSP)=ZR2D(JLEV,JSP)*(1.0_JPRB-ZBTCM*YDLEP%RLEPDIM(ISP))
        ENDDO
      ENDDO
    ENDIF
  ELSE
    ! ZSVED is simply 0 in this case.
    ZSVED(1:NFLEVG,KSTA:KEND)=0.0_JPRB
  ENDIF

  ! * Provides "SIFAC * RHS of the Helmholtz eqn" in ZSRHS2.
  DO JLEV=1,NFLEVG
!cdir nodep
    DO JSP=KSTA,KEND
      ZSRHS2(JLEV,JSP)=ZBTCMH*&
       & (ZSRHS(JLEV,JSP)-ZWORK(JLEV,JSP))+&
       & ZSVED(JLEV,JSP)  
    ENDDO
  ENDDO

  ! * Case LIMPF: add terms containing Dcha_star_star(m,n),
  !   Dcha_star_star(m,n-2), Dcha_star_star(m,n+2) in the RHS.
  ! IF (LIMPF .AND. (JITER == 0)) THEN
  !   LIMPF not coded for ALADIN-NH.
  ! ENDIF

  ! * Multiply by array SIFACI to obtain the RHS of the Helmholtz eqn,
  !   stored in ZSRHS.
  CALL MXMAOP(SIFACI,1,NFLEVG,ZSRHS2,1,NFLEVG,ZSRHS,1,NFLEVG,&
   & NFLEVG,NFLEVG,ISPCOL)

  !*        2.4  Solve Helmholtz equation (compute dver(t+dt))

  ! Current space --> vertical eigenmodes space
  ! (multiply by matrix Q stored in SIMI).

  CALL MXMAOP(SIMI,1,NFLEVG,ZSRHS,1,NFLEVG,ZSDIVP,1,NFLEVG,&
   & NFLEVG,NFLEVG,ISPCOL)

  ! Inversion of Helmholtz equation itself.

  IF (LESIDG) THEN
    ! - memory transfer zsdivp -> zsdivpl
    ! - km such as "n always >0": call to MXTURS, result in zspdivpl
    ! - km such as "n can be 0": two calls to MXTURE, result in zspdivpl
    ! - memory transfer zspdivpl -> zzspsvdg
    ZSDIVPL(:,:,:)=0.0_JPRB
    ZSPDIVPL(:,:,:)=0.0_JPRB
    DO JN=0,NISNAX(KM)
      ISE=KSTA+4*JN
      ZSDIVPL(:,JN,1:4)=ZSDIVP(:,ISE:ISE+3)
    ENDDO
    IF (KM > 0) THEN
      ! Inversion of a symmetric penta-diagonal matrix,
      !  provides grad'**2*Q*dver(t+dt).
      CALL MXTURS(NISNAX(KM)+1,NFLEVG,NFLEVG,II,&
       & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
       & ZSDIVPL,ZSPDIVPL)
      ! Multiply by grad'**(-2) and put the result in ZSPDIVP.
      DO JN=0,NISNAX(KM)
        ISE=KSTA+4*JN
        IN = YDLAP%NVALUE(ISE+IOFF)
        IM = YDLEP%MVALUE(ISE+IOFF)
        ISP = YDLEP%NPME(IM)+IN
        ZSPDIVP(:,ISE:ISE+3)=ZSPDIVPL(:,JN,1:4)*YDLEP%RLEPINM(ISP)
      ENDDO
    ELSE
      ! Inversion of a non-symmetric penta-diagonal matrix,
      !  provides Q*dver(t+dt).
      CALL MXTURE(NISNAX(KM)+1,NFLEVG,NFLEVG,II,-2,.TRUE.,&
       & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
       & ZSDIVPL,ZSPDIVPL)
      CALL MXTURE(NISNAX(KM)+1,NFLEVG,NFLEVG,II,3,.FALSE.,&
       & SIHEG(1,IS0+1,1),SIHEG2(1,IS02+1,2),&
       & SIHEG2(1,IS02+1,3),ZSDIVPL,ZSPDIVPL)
      ! Put the result in ZSPDIVP.
      DO JN=0,NISNAX(KM)
        ISE=KSTA+4*JN
        ZSPDIVP(:,ISE:ISE+3)=ZSPDIVPL(:,JN,1:4)
      ENDDO
    ENDIF
  ELSE

    ! Case designed when no stretching:

    IF (LIMPF) THEN
      ! LIMPF not coded for ALADIN-NH.
    ELSE

      ! Inversion of a diagonal matrix, provides Q*dver(t+dt).
      DO JLEV=1,NFLEVG
!cdir nodep
        DO JSP=KSTA,KEND
          IN=YDLAP%NVALUE(JSP+IOFF)
          IM=YDLEP%MVALUE(JSP+IOFF)
          ISP=YDLEP%NPME(IM)+IN
          ZSPDIVP(JLEV,JSP)=ZSDIVP(JLEV,JSP)&
           & /(1.0_JPRB-ZBDT2*SIVP(JLEV)*YDLEP%RLEPDIM(ISP))
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! Vertical eigenmodes space --> current space.
  ! (multiply by matrix Q**(-1) stored in SIMO).
  ! Provides dver(t+dt).
   CALL MXMAOP(SIMO,1,NFLEVG,ZSPDIVP,1,NFLEVG,ZZSPSVDG,1,&
   & NFLEVG,NFLEVG,NFLEVG,ISPCOL)

  !*        2.5  Recover the other prognostic variables
  !              (successively D, T, "spd" and log(prehyds)).

  DO JSP=KSTA,KEND
    ZSP(JSP) = 0.0_JPRB
  ENDDO

  ! * Provides "Gamma d(t+dt)" in array ZWORK.
  CALL SIGAM(LVERTFE, YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZWORK,ZZSPSVDG,ZSP,ISPCOL,NFLEVG)

  ! * Provides
  !   " beta**2 (Delta t)**2 vnabla'**2 (- SITR Gamma + C**2) d(t+dt)
  !   + Dprim_star_star "
  !   in array ZSRHS.
  DO JLEV=1,NFLEVG
!cdir nodep
    DO JSP=KSTA,KEND
      IN=YDLAP%NVALUE(JSP+IOFF)
      IM=YDLEP%MVALUE(JSP+IOFF)
      ISP=YDLEP%NPME(IM)+IN
      ZSRHS(JLEV,JSP)=(ZZSPSVDG(JLEV,JSP)*YDCST%RD*ZCC&
       & -ZWORK(JLEV,JSP))*SITR*ZBDT*ZBDT*&
       & YDLEP%RLEPDIM(ISP)+ZR1D(JLEV,JSP)  
    ENDDO
  ENDDO

  ! * Provides D'(t+dt).
  IF (LESIDG) THEN
    ! Divide by the pentadiagonal operator
    ! [I - beta**2 (Delta t)**2 vnabla'**2 C**2 M**2]
    ! (LU decomposition to be done in SUNHHEG,
    ! calls to MXTURS or MXTURE according to KM to be done).
    ! For the KM where "YDLEP%RLEPDIM" is never 0
    ! it is desirable to work with the symmetric
    ! operator [vnabla'**(-2) - beta**2 (Delta t)**2 C**2 M**2]
    ! after having multiplied the RHS "ZSRHS" by "YDLEP%RLEPINM".
    ZSDIVPL(:,:,:)=0.0_JPRB
    ZSPDIVPL(:,:,:)=0.0_JPRB
    IF (KM > 0) THEN
      ! Multiply the RHS by grad'**(-2).
      DO JN=0,NISNAX(KM)
        ISE=KSTA+4*JN
        IN = YDLAP%NVALUE(ISE+IOFF)
        IM = YDLEP%MVALUE(ISE+IOFF)
        ISP = YDLEP%NPME(IM)+IN
        ZSDIVPL(:,JN,1:4)=ZSRHS(:,ISE:ISE+3)*YDLEP%RLEPINM(ISP)
      ENDDO
      ! Inversion of a symmetric penta-diagonal matrix.
      CALL MXTURS(NISNAX(KM)+1,NFLEVG,NFLEVG,II,&
       & SIHEGB(1,IS0+1,1),SIHEGB(1,IS0+1,2),SIHEGB(1,IS0+1,3),&
       & ZSDIVPL,ZSPDIVPL)
    ELSE
      DO JN=0,NISNAX(KM)
        ISE=KSTA+4*JN
        ZSDIVPL(:,JN,1:4)=ZSRHS(:,ISE:ISE+3)
      ENDDO
      ! Inversion of a non-symmetric penta-diagonal matrix.
      CALL MXTURE(NISNAX(KM)+1,NFLEVG,NFLEVG,II,-2,.TRUE.,&
       & SIHEGB(1,IS0+1,1),SIHEGB(1,IS0+1,2),SIHEGB(1,IS0+1,3),&
       & ZSDIVPL,ZSPDIVPL)
      CALL MXTURE(NISNAX(KM)+1,NFLEVG,NFLEVG,II,3,.FALSE.,&
       & SIHEGB(1,IS0+1,1),SIHEGB2(1,IS02+1,2),&
       & SIHEGB2(1,IS02+1,3),ZSDIVPL,ZSPDIVPL)
    ENDIF
    ! Put the result in ZZSPDIVG.
    DO JN=0,NISNAX(KM)
      ISE=KSTA+4*JN
      ZZSPDIVG(:,ISE:ISE+3)=ZSPDIVPL(:,JN,1:4)
    ENDDO
  ELSE

    ! * Retrieve D'(t+dt).

    IF (LIMPF) THEN
      ! LIMPF not coded for ALADIN-NH.
    ELSE
      ! * Division by the diagonal operator
      !   [I - beta**2 (Delta t)**2 vnabla'**2 C**2 RSTRET**2].
      DO JLEV=1,NFLEVG
!cdir nodep
        DO JSP=KSTA,KEND
          IN=YDLAP%NVALUE(JSP+IOFF)
          IM=YDLEP%MVALUE(JSP+IOFF)
          ISP=YDLEP%NPME(IM)+IN
          ZZSPDIVG(JLEV,JSP)=ZSRHS(JLEV,JSP)/(1.0_JPRB-ZBTCM*YDLEP%RLEPDIM(ISP))
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! * Provides Mbar**2 D'(t+dt).
  IF (LESIDG) THEN
    ! Apply a spectral multiplication by M**2
    ! (CALL MXPTMA using SCGMAP, see SPCSI).
    ZSDIVPL(:,:,:)=0.0_JPRB
    ZSPDIVPL(:,:,:)=0.0_JPRB
    DO JN=0,NISNAX(KM)
      ISE=KSTA+4*JN
      ZSDIVPL(:,JN,1:4)=ZZSPDIVG(:,ISE:ISE+3)
    ENDDO
    CALL MXPTMA(NISNAX(KM)+1,NFLEVG,NFLEVG,II,&
     & ZPD,ZPE,ZPF,ZPE,ZPF,&
     & ZSDIVPL,ZSPDIVPL)
    DO JN=0,NISNAX(KM)
      ISE=KSTA+4*JN
      ZWORK(:,ISE:ISE+3)=ZSPDIVPL(:,JN,1:4)
    ENDDO
  ELSE
    DO JLEV=1,NFLEVG
!cdir nodep
      DO JSP=KSTA,KEND
        ZWORK(JLEV,JSP)=ZZSPDIVG(JLEV,JSP)*RSTRET*RSTRET
      ENDDO
    ENDDO
  ENDIF

  !*       2.6  Increment vorticity for case LIMPF=T

  ! IF (LIMPF) THEN
  !   LIMPF not coded for ALADIN-NH.
  ! ENDIF

  !*       2.7  Savings.

  IF (I_NITERHELM > 0) THEN
    IF (JITER == 0) THEN
      ! Save D'(t+dt,iter=0), Mbar**2 D'(t+dt,iter=0),
      !  zeta'(t+dt,iter=0), dver(t+dt,iter=0).
      ZZSPDIV0G(1:NFLEVG,KSTA:KEND)=ZZSPDIVG(1:NFLEVG,KSTA:KEND)
      ZZSPGDIV0G(1:NFLEVG,KSTA:KEND)=ZWORK(1:NFLEVG,KSTA:KEND)
      IF (LIMPF) ZZSPVOR0G(1:NFLEVG,KSTA:KEND)=ZZSPVORG(1:NFLEVG,KSTA:KEND)
      ZZSPSVD0G(1:NFLEVG,KSTA:KEND)=ZZSPSVDG(1:NFLEVG,KSTA:KEND)
    ENDIF

    ! Save Mbar**2 D'(t+dt,jiter) in ZZSPGDIVG
    IF (JITER == 0) THEN
      ZZSPGDIVG(1:NFLEVG,KSTA:KEND)=ZWORK(1:NFLEVG,KSTA:KEND)
    ELSE
      DO JLEV=1,NFLEVG
        DO JSP=KSTA,KEND
          ZZSPGDIVG(JLEV,JSP)=ZWORK(JLEV,JSP)+ZZSPGDIV0G(JLEV,JSP)
        ENDDO
      ENDDO
    ENDIF
  ENDIF

ENDDO ! End loop on JITER

!*       2.8  Add the increment and recover X(t+dt)
!             for X=D',Mbar**2 D',zeta',dver, and deallocate.

IF (I_NITERHELM > 0) THEN

  DO JLEV=1,NFLEVG
!cdir nodep
    DO JSP=KSTA,KEND
      ZZSPDIVG(JLEV,JSP)=ZZSPDIVG(JLEV,JSP)+ZZSPDIV0G(JLEV,JSP)
      ZWORK(JLEV,JSP)=ZWORK(JLEV,JSP)+ZZSPGDIV0G(JLEV,JSP)
      ZZSPSVDG(JLEV,JSP)=ZZSPSVDG(JLEV,JSP)+ZZSPSVD0G(JLEV,JSP)
    ENDDO
  ENDDO
  IF (LIMPF) THEN
    DO JLEV=1,NFLEVG
!cdir nodep
      DO JSP=KSTA,KEND
        ZZSPVORG(JLEV,JSP)=ZZSPVORG(JLEV,JSP)+ZZSPVOR0G(JLEV,JSP)
      ENDDO
    ENDDO
  ENDIF

ENDIF

!*       2.9  Increment T, log(prehyds) and spd.

! * Provides some intermediate quantities allowing to compute
!   T(t+dt) (in ZST), log(prehyds)(t+dt) (in ZSP), spd(t+dt) (in ZSNHP).
CALL SIPTP(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZWORK,ZZSPSVDG,ZSNHP,ZST,ZSP,ISPCOL)

! * Provides T(t+dt) (in ZZSPTG) and spd(t+dt) (in ZZSPSPDG).
DO JLEV=1,NFLEVG
!cdir nodep
  DO JSP=KSTA,KEND
    ZZSPTG(JLEV,JSP)=ZZSPTG(JLEV,JSP)-ZBDT*ZST(JLEV,JSP)
    ZZSPSPDG(JLEV,JSP)=ZZSPSPDG(JLEV,JSP)-ZBDT*ZSNHP(JLEV,JSP)
  ENDDO
ENDDO

! * Provides log(prehyds)(t+dt) (in ZZSPSPG).
!cdir nodep
DO JSP=KSTA,KEND
  ZZSPSPG(JSP)=ZZSPSPG(JSP)-ZBDT*ZSP(JSP)
ENDDO

IF (LESIDG) THEN
  DEALLOCATE(ZSDIVPL)
  DEALLOCATE(ZSPDIVPL)
ENDIF

! ky: if LESIDG, missing there the reverse shift.

!     ------------------------------------------------------------------

!*       3.    MEMORY TRANSFER AND DDH-SI UPDATE.
!              ----------------------------------

IF (LIMPF) PSPVORG(1:NFLEVG,KSTA:KEND)=ZZSPVORG(1:NFLEVG,KSTA:KEND)
PSPDIVG(1:NFLEVG,KSTA:KEND)=ZZSPDIVG(1:NFLEVG,KSTA:KEND)
PSPTG  (1:NFLEVG,KSTA:KEND)=ZZSPTG  (1:NFLEVG,KSTA:KEND)
PSPSPG (         KSTA:KEND)=ZZSPSPG (         KSTA:KEND)
PSPSPDG(1:NFLEVG,KSTA:KEND)=ZZSPSPDG(1:NFLEVG,KSTA:KEND)
PSPSVDG(1:NFLEVG,KSTA:KEND)=ZZSPSVDG(1:NFLEVG,KSTA:KEND)

IF (LRSIDDH) THEN
  IF (LIMPF) PSPTNDSI_VORG(1:NFLEVG,KSTA:KEND)=&
   & PSPTNDSI_VORG(1:NFLEVG,KSTA:KEND)+PSPVORG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_DIVG(1:NFLEVG,KSTA:KEND)=PSPTNDSI_DIVG(1:NFLEVG,KSTA:KEND)&
   & + PSPDIVG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_TG(1:NFLEVG,KSTA:KEND)=PSPTNDSI_TG(1:NFLEVG,KSTA:KEND)&
   & + PSPTG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_SPDG(1:NFLEVG,KSTA:KEND)=PSPTNDSI_SPDG(1:NFLEVG,KSTA:KEND)&
   & + PSPSPDG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_SVDG(1:NFLEVG,KSTA:KEND)=PSPTNDSI_SVDG(1:NFLEVG,KSTA:KEND)&
   & + PSPSVDG(1:NFLEVG,KSTA:KEND)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ESPNHSI',1,ZHOOK_HANDLE)
END SUBROUTINE ESPNHSI
