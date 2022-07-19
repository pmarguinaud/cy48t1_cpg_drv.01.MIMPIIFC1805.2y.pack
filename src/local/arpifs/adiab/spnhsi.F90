SUBROUTINE SPNHSI(&
 ! --- INPUT -----------------------------------------------------------------
 & YDCST, YDGEOMETRY,YDLDDH,YDRIP,YDDYN,KM,KMLOC,KSTA,KEND,LDONEM,&
 ! --- INOUT -----------------------------------------------------------------
 & PSPVORG,PSPDIVG,PSPTG,PSPSPG,PSPSPDG,PSPSVDG,&
 & PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,PSPTNDSI_SPDG,PSPTNDSI_SVDG)

!**** *SPNHSI* - SPECTRAL SPACE SEMI-IMPLICIT COMPUTATIONS FOR NH MODEL.

!     Purpose.
!     --------

!     Semi-implicit scheme in the NH model with (Qcha, dver or d4) as
!     NH prognostic variables.

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
!        *CALL* *SPNHSI(..)

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

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        K. Yessad (after SPCSI and ESPC).
!        Original : 09-Dec-2004

!     Modifications.
!     --------------
!   K. Yessad Aug-2007: iterative algorithm to relax the C1 constraint
!   K. Yessad Aug-2007: LIMPF in NH model.
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Aug 2009): remove LSITRIC option
!   K. Yessad (Feb 2012): tests on LL3D, LLDOSI in the caller, simplifications.
!   G. Mozdzynski (Oct 2012): improve NH performance at scale
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (July 2014): Move some variables.
!   O. Marsden (May 2016): Remove redundant geometry argument
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Vertical-dependent SITRA.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : TCST
USE YOMMP0       , ONLY : MYSETV
USE YOMDYN       , ONLY : TDYN
USE YOMDYNA      , ONLY : YRDYNA
USE YOMLDDH      , ONLY : TLDDH
USE YOMRIP       , ONLY : TRIP
USE YOMCVER      , ONLY : LVERTFE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TLDDH)       ,INTENT(IN)    :: YDLDDH
TYPE(TRIP)        ,INTENT(IN)    :: YDRIP
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

REAL(KIND=JPRB),ALLOCATABLE :: ZZSPGDIVG(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZZSPDIV0G(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZZSPGDIV0G(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZZSPVOR0G(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZZSPSVD0G(:,:)

REAL(KIND=JPRB) :: ZSDIV(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZST(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSP(KSTA:KEND)
REAL(KIND=JPRB) :: ZSDIVP(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSPDIVP(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSVED(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZR1D(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZR1DPRIM(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZR2D(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZR3D(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZR4D(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSRHS(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZWORK(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSRHS2(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSRHS2_INC(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSNHP(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSDIVPL(YDGEOMETRY%YRDIMV%NFLEVG,KM:YDGEOMETRY%YRDIM%NSMAX,2)
REAL(KIND=JPRB) :: ZSPDIVPL(YDGEOMETRY%YRDIMV%NFLEVG,KM:YDGEOMETRY%YRDIM%NSMAX,2)

REAL(KIND=JPRB) :: ZALPHA (KM:YDGEOMETRY%YRDIM%NSMAX+1)
REAL(KIND=JPRB) :: ZDENIM (KM:YDGEOMETRY%YRDIM%NSMAX+1)
REAL(KIND=JPRB) :: ZEPSI  (KM:YDGEOMETRY%YRDIM%NSMAX)
REAL(KIND=JPRB) :: ZFPLUS (KM:YDGEOMETRY%YRDIM%NSMAX+1)
REAL(KIND=JPRB) :: ZFMINUS(KM:YDGEOMETRY%YRDIM%NSMAX+1)

REAL(KIND=JPRB) :: ZSIVP2(YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: II, IN, IOFF, IS0, IS02, ISE, ISPCOL, JLEV, JSP, JN
INTEGER(KIND=JPIM) :: ILO, JL, IRSP, IMSP
INTEGER(KIND=JPIM) :: I_NITERHELM, JITER

REAL(KIND=JPRB) :: ZBDT, ZBDT2, ZRRT, ZCC, ZBTCM, ZBTCMH
REAL(KIND=JPRB) :: ZEM, ZAL, ZEN, ZF, ZTEMP

LOGICAL :: LLNH_NOC1

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "mxmaop.h"

#include "mxptma.h"
#include "mxture.h"
#include "mxturs.h"
#include "siseve.intfb.h"
#include "si_cccor.intfb.h"
#include "sidd.intfb.h"
#include "sigam.intfb.h"
#include "simplico.h"
#include "si_mxptco.h"
#include "siptp.intfb.h"
#include "sitnu.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPNHSI',0,ZHOOK_HANDLE)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
 & YDDIMV=>YDGEOMETRY%YRDIMV, &
 & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,  &
 & YDLAP=>YDGEOMETRY%YRLAP, YDSPGEOM=>YDGEOMETRY%YSPGEOM)

ASSOCIATE(NSMAX=>YDDIM%NSMAX, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LIMPF=>YDDYN%LIMPF, LSIDG=>YDDYN%LSIDG, NITERHELM=>YDDYN%NITERHELM, &
 & RBTS2=>YDDYN%RBTS2, SIFACI=>YDDYN%SIFACI, SIHEG=>YDDYN%SIHEG, &
 & SIHEG2=>YDDYN%SIHEG2, SIHEGB=>YDDYN%SIHEGB, SIHEGB2=>YDDYN%SIHEGB2, &
 & SIMI=>YDDYN%SIMI, SIMO=>YDDYN%SIMO, SITR=>YDDYN%SITR, &
 & SIVP=>YDDYN%SIVP, &
 & RSTRET=>YDGEM%RSTRET, &
 & LRSIDDH=>YDLDDH%LRSIDDH, &
 & NPTRSV=>YDMP%NPTRSV, NPTRSVF=>YDMP%NPTRSVF, &
 & TDT=>YDRIP%TDT, &
 & SCGMAP=>YDSPGEOM%SCGMAP)

!     ------------------------------------------------------------------

! * Case where the C1 constraint is not matched.
LLNH_NOC1=.NOT.(YRDYNA%NDLNPR == 1)
IF (LLNH_NOC1) THEN
  I_NITERHELM=NITERHELM
ELSE
  I_NITERHELM=0
ENDIF

!     ------------------------------------------------------------------

!*       1.    MEMORY TRANSFER.
!              ----------------

  IF (LIMPF) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
        ZZSPVORG(JLEV,JSP)=PSPVORG(JLEV,JSP)  
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDIF
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      ZZSPDIVG(JLEV,JSP)=PSPDIVG(JLEV,JSP)
      ZZSPTG  (JLEV,JSP)=PSPTG  (JLEV,JSP)
      ZZSPSPDG(JLEV,JSP)=PSPSPDG(JLEV,JSP)
      ZZSPSVDG(JLEV,JSP)=PSPSVDG(JLEV,JSP)
    ENDDO
    ZZSPSPG (JSP)=PSPSPG (JSP)
  ENDDO
!$OMP END PARALLEL DO
IF (LRSIDDH) THEN
  ! DDH memory transfer
  IF (LIMPF) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
        PSPTNDSI_VORG(JLEV,JSP)=-PSPVORG(JLEV,JSP)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDIF
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      PSPTNDSI_DIVG(JLEV,JSP)=-PSPDIVG(JLEV,JSP)
      PSPTNDSI_TG  (JLEV,JSP)=-PSPTG  (JLEV,JSP)
      PSPTNDSI_SPDG(JLEV,JSP)=-PSPSPDG(JLEV,JSP)
      PSPTNDSI_SVDG(JLEV,JSP)=-PSPSVDG(JLEV,JSP)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
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
ZSIVP2(1:NFLEVG)=YDCST%RD*SITR*(YDCST%RCPD/YDCST%RCVD)

IF (LSIDG) THEN
  II=MIN(KM,1)+1
  IS0=YDLAP%NSE0L(KMLOC)
  IS02=0
ENDIF

IF (I_NITERHELM > 0) THEN
  ALLOCATE(ZZSPGDIVG(1:NFLEVG,KSTA:KEND))
  ALLOCATE(ZZSPDIV0G(1:NFLEVG,KSTA:KEND))
  ALLOCATE(ZZSPGDIV0G(1:NFLEVG,KSTA:KEND))
  ALLOCATE(ZZSPVOR0G(1:NFLEVG,KSTA:KEND))
  ALLOCATE(ZZSPSVD0G(1:NFLEVG,KSTA:KEND))
ENDIF

!*        2.2  Set up helper arrays for implicit Coriolis case.
!              (identical to the part 2.2 of SPCSI).

IF (LIMPF) THEN
  ZEM=REAL(KM,JPRB)
  ZAL=2.0_JPRB*ZBDT*YDCST%ROMEGA*ZEM
  ILO=KM
  IF (KM == 0) THEN
    ZALPHA(0)=0.0_JPRB
    ZDENIM(0)=0.0_JPRB
    ZEPSI(0)=0.0_JPRB
    ILO=1
  ENDIF
  DO JL=ILO,NSMAX
    ZEN=REAL(JL,JPRB)
    ZALPHA(JL)=ZAL/(ZEN*(ZEN+1.0_JPRB))
    ZDENIM(JL)=1.0_JPRB/(1.0_JPRB+ZALPHA(JL)**2)
    ZEPSI(JL)=SQRT((ZEN*ZEN-ZEM*ZEM)/(4.0_JPRB*ZEN*ZEN-1.0_JPRB))
  ENDDO
  ZALPHA(NSMAX+1)=0.0_JPRB
  ZDENIM(NSMAX+1)=0.0_JPRB

  IF (KM == 0) THEN
    ZFPLUS(0)=0.0_JPRB
    ZFMINUS(0)=0.0_JPRB
  ENDIF
  ZF=2.0_JPRB*ZBDT*YDCST%ROMEGA
  DO JL=ILO,NSMAX-1
    ZEN=REAL(JL,JPRB)
    ZFPLUS(JL)=ZF*ZEN*ZEPSI(JL+1)/(ZEN+1.0_JPRB)
    ZFMINUS(JL)=ZF*(ZEN+1.0_JPRB)*ZEPSI(JL)/ZEN
  ENDDO
  ZEN=REAL(NSMAX,JPRB)
  ZFPLUS(NSMAX)=0.0_JPRB
  ZFMINUS(NSMAX)=ZF*(ZEN+1.0_JPRB)*ZEPSI(NSMAX)/ZEN
  ZFPLUS(NSMAX+1)=0.0_JPRB
  ZFMINUS(NSMAX+1)=0.0_JPRB
ENDIF

DO JITER=0,I_NITERHELM

  !*        2.3  Computes right-hand side of Helmholtz equation.

  IF (JITER == 0) THEN

    ! * Provides:
    !   a/ - [ gamma [SITR Qcha_rhs - T_rhs] - (Rd SITR)/(SIPR) logprehyd_rhs ]
    !      in array ZSDIV.
    !   b/ (g TRSI)/(H TARSI) LLstar Qcha_rhs in array ZSVED.
    CALL SIDD(YRDYNA, YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZSDIV,ZSVED,ZZSPSPDG,ZZSPTG,ZZSPSPG,ISPCOL)

    ! * Provides Dprim_star_star in ZSDIV.
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV,IN)
      DO JSP=KSTA,KEND
        DO JLEV=1,NFLEVG
        IN=YDLAP%NVALUE(JSP+IOFF)
        ZSDIV(JLEV,JSP)=ZZSPDIVG(JLEV,JSP)-ZBDT*YDLAP%RLAPDI(IN)*ZSDIV(JLEV,JSP)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO

  ELSE

    ! * Compute C**2 * COR * Mbar**2 D'(t+dt,jiter-1)
    CALL SI_CCCOR(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ISPCOL,ZZSPGDIVG,ZSDIV)

    ! * Multiply by beta**2 (Delta t)**2 vnabla'**2
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV,IN)
      DO JSP=KSTA,KEND
        DO JLEV=1,NFLEVG
        IN=YDLAP%NVALUE(JSP+IOFF)
        ZSDIV(JLEV,JSP)=ZBDT*ZBDT*YDLAP%RLAPDI(IN)*ZSDIV(JLEV,JSP)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO

  ENDIF

  ! * Provides Dprim_star_star_star (in ZSDIV) from Dprim_star_star (in ZSDIV)
  !   (the quantity added is the same one as in the end of part 2.3 of SPCSI).
  IF (LIMPF) THEN
    IF (JITER == 0) THEN

      IF (KM > 0) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,JLEV,IRSP,IMSP,ZTEMP)
          DO JN=KM,NSMAX
            DO JLEV=1,NFLEVG
            IRSP=KSTA+(JN-KM)*2
            IMSP=IRSP+1
            ZTEMP=ZDENIM(JN)*(ZZSPVORG(JLEV,IMSP)&
             & +ZALPHA(JN)*ZZSPVORG(JLEV,IRSP))
            ZZSPVORG(JLEV,IRSP)=ZDENIM(JN)*(ZZSPVORG(JLEV,IRSP)&
             & -ZALPHA(JN)*ZZSPVORG(JLEV,IMSP))
            ZZSPVORG(JLEV,IMSP)=ZTEMP
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !        Add [F] * result to rhs of Helmholtz equation

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,JLEV,IRSP,IMSP)
        DO JN=KM+1,NSMAX
          DO JLEV=1,NFLEVG
          IRSP=KSTA+(JN-KM)*2
          IMSP=IRSP+1
          ZSDIV(JLEV,IRSP)=ZSDIV(JLEV,IRSP)&
           & + ZFMINUS(JN)*ZZSPVORG(JLEV,IRSP-2)
          ZSDIV(JLEV,IMSP)=ZSDIV(JLEV,IMSP)&
           & + ZFMINUS(JN)*ZZSPVORG(JLEV,IMSP-2)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,JLEV,IRSP,IMSP)
        DO JN=KM,NSMAX-1
          DO JLEV=1,NFLEVG
          IRSP=KSTA+(JN-KM)*2
          IMSP=IRSP+1
          ZSDIV(JLEV,IRSP)=ZSDIV(JLEV,IRSP)+ ZFPLUS(JN)*ZZSPVORG(JLEV,IRSP+2)
          ZSDIV(JLEV,IMSP)=ZSDIV(JLEV,IMSP)+ ZFPLUS(JN)*ZZSPVORG(JLEV,IMSP+2)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO

    ELSE

      ! zeta_star is replaced by zero.
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
        DO JSP=KSTA,KEND
          DO JLEV=1,NFLEVG
            ZZSPVORG(JLEV,JSP)=0.0_JPRB
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

    ENDIF
  ENDIF

  ! * Provides Dprim_star_star (resp. Dprim_star_star_star) in ZR1D and
  !   ZR1DPRIM if LIMPF=F (resp. LIMPF=T); and Dcha_star_star in ZR2D.
  IF (JITER == 0) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
      DO JSP=KSTA,KEND
        DO JLEV=1,NFLEVG
        ZR2D(JLEV,JSP)=ZZSPSVDG(JLEV,JSP)-ZBDT*ZSVED(JLEV,JSP)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ELSE
    ! Dcha_star_star is replaced by zero.
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
      DO JSP=KSTA,KEND
        DO JLEV=1,NFLEVG
          ZR2D(JLEV,JSP)=0.0_JPRB
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
  ENDIF
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
        ZR1D(JLEV,JSP)=ZSDIV(JLEV,JSP)
        ZR1DPRIM(JLEV,JSP)=ZSDIV(JLEV,JSP)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO

  IF (LSIDG) THEN
    ! Multiply ZR1D by M**2 (CALL MXPTMA using SCGMAP),
    !  then divide it by RSTRET**2.
    !  Use ZSDIVPL and ZSPDIVPL as intermediate work arrays.
!$OMP WORKSHARE
    ZSDIVPL(:,:,:)=0.0_JPRB
    ZSPDIVPL(:,:,:)=0.0_JPRB
!$OMP END WORKSHARE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,ISE)
    DO JN=KM,NSMAX
      ISE=KSTA+2*(JN-KM)
      ZSDIVPL(:,JN,1:2)=ZR1D(:,ISE:ISE+1)
    ENDDO
!$OMP END PARALLEL DO
    CALL MXPTMA(NSMAX+1-KM,NFLEVG,NFLEVG,II,SCGMAP(IS0+1,1),&
     & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
     & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
     & ZSDIVPL,ZSPDIVPL)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,ISE)
    DO JN=KM,NSMAX
      ISE=KSTA+2*(JN-KM)
      ZR1D(:,ISE:ISE+1)=ZSPDIVPL(:,JN,1:2)/(RSTRET*RSTRET)
    ENDDO
!$OMP END PARALLEL DO
  ENDIF

  ! * Provides "LLstar Dprim_star_star" if LIMPF=F, or
  !   "LLstar Dprim_star_star_star" if LIMPF=T
  !   (mult by a constant coefficient) in ZSRHS.
  CALL SISEVE(YRDYNA,YDGEOMETRY,YDDYN,1,NFLEVG,ZR1D,ZSRHS,ISPCOL)

  ! * Provides "Tau Dprim_star_star" if LIMPF=F, or
  !   "Tau Dprim_star_star_star" if LIMPF=T (mult by a constant coefficient)
  !   in ZST.
  CALL SITNU(LVERTFE, YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZR1D,ZST,ZSP,ISPCOL)

  ! * Provides "LLstar Tau Dprim_star_star" if LIMPF=F, or
  !   "LLstar Tau Dprim_star_star_star" if LIMPF=T (mult by a constant
  !   coefficient) in ZWORK.
  CALL SISEVE(YRDYNA,YDGEOMETRY,YDDYN,1,NFLEVG,ZST,ZWORK,ISPCOL)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV,IN)
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
      ZWORK(JLEV,JSP)=ZWORK(JLEV,JSP)*YDCST%RCVD*ZRRT
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

  ! * Provides
  !   "[I - beta**2 (Delta t)**2 C**2 Mbar**2 vnabla'**2] Dcha_star_star"
  !   in ZSVED.
  IF (JITER == 0) THEN
    IF (LSIDG) THEN
      ! First compute ZR3D=ZR2D*RLAPDI.
      ! Multiply ZR3D by M**2 (CALL MXPTMA using SCGMAP),
      !  put the result in a separate array ZR4D,
      !  use ZSDIVPL and ZSPDIVPL as intermediate work arrays.
      ! Compute ZSVED using ZR2D and ZR4D:
      !  ZSVED=ZR2D-(ZBTCM/(RSTRET*RSTRET))*ZR4D.
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV,IN)
        DO JSP=KSTA,KEND
          DO JLEV=1,NFLEVG
          IN=YDLAP%NVALUE(JSP+IOFF)
          ZR3D(JLEV,JSP)=YDLAP%RLAPDI(IN)*ZR2D(JLEV,JSP)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
!$OMP WORKSHARE
      ZSDIVPL(:,:,:)=0.0_JPRB
      ZSPDIVPL(:,:,:)=0.0_JPRB
!$OMP END WORKSHARE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,ISE)
      DO JN=KM,NSMAX
        ISE=KSTA+2*(JN-KM)
        ZSDIVPL(:,JN,1:2)=ZR3D(:,ISE:ISE+1)
      ENDDO
!$OMP END PARALLEL DO
      CALL MXPTMA(NSMAX+1-KM,NFLEVG,NFLEVG,II,SCGMAP(IS0+1,1),&
       & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
       & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
       & ZSDIVPL,ZSPDIVPL)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,ISE)
      DO JN=KM,NSMAX
        ISE=KSTA+2*(JN-KM)
        ZR4D(:,ISE:ISE+1)=ZSPDIVPL(:,JN,1:2)
      ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
        DO JSP=KSTA,KEND
          DO JLEV=1,NFLEVG
          ZSVED(JLEV,JSP)=ZR2D(JLEV,JSP)-(ZBTCM/(RSTRET*RSTRET))*ZR4D(JLEV,JSP)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV,IN)
        DO JSP=KSTA,KEND
          DO JLEV=1,NFLEVG
          IN=YDLAP%NVALUE(JSP+IOFF)
          ZSVED(JLEV,JSP)=ZR2D(JLEV,JSP)-ZBTCM*ZR2D(JLEV,JSP)*YDLAP%RLAPDI(IN)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    ENDIF
  ELSE
    ! ZSVED is simply 0 in this case.
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
      DO JSP=KSTA,KEND
        DO JLEV=1,NFLEVG
          ZSVED(JLEV,JSP)=0.0_JPRB
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
  ENDIF

  ! * Provides "SIFAC * RHS of the Helmholtz eqn" in ZSRHS2.
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
      ZSRHS2(JLEV,JSP)=ZBTCMH*&
       & (ZSRHS(JLEV,JSP)-ZWORK(JLEV,JSP))+&
       & ZSVED(JLEV,JSP)  
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

  ! * Case LIMPF: add terms containing Dcha_star_star(m,n),
  !   Dcha_star_star(m,n-2), Dcha_star_star(m,n+2) in the RHS.

  IF (LIMPF .AND. (JITER == 0)) THEN
    CALL SI_MXPTCO(KM,NSMAX,NFLEVG,NFLEVG,ZF,ZALPHA(KM),&
     & ZDENIM(KM),ZEPSI(KM),ZR2D,ZSRHS2_INC)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
      DO JSP=KSTA,KEND
        DO JLEV=1,NFLEVG
          ZSRHS2(JLEV,JSP)=ZSRHS2(JLEV,JSP)&
           & +ZSRHS2_INC(JLEV,JSP)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
  ENDIF

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

  IF (LSIDG) THEN

    ! Case designed for stretching :
    ! Use ZSDIVPL and ZSPDIVPL as intermediate work arrays.

!$OMP WORKSHARE
    ZSDIVPL(:,:,:)=0.0_JPRB
    ZSPDIVPL(:,:,:)=0.0_JPRB
!$OMP END WORKSHARE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,ISE)
    DO JN=KM,NSMAX
      ISE=KSTA+2*(JN-KM)
      ZSDIVPL(:,JN,1:2)=ZSDIVP(:,ISE:ISE+1)
    ENDDO
!$OMP END PARALLEL DO
    IF (KM > 0) THEN
      ! Inversion of a symmetric penta-diagonal matrix,
      !  provides grad'**2*Q*dver(t+dt).
      CALL MXTURS(NSMAX+1-KM,NFLEVG,NFLEVG,II,&
       & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
       & ZSDIVPL,ZSPDIVPL)  
      ! Multiply by grad'**(-2) and put the result in ZSPDIVP.
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,ISE)
      DO JN=KM,NSMAX
        ISE=KSTA+2*(JN-KM)
        ZSPDIVP(:,ISE:ISE+1)=ZSPDIVPL(:,JN,1:2)*YDLAP%RLAPIN(JN)
      ENDDO
!$OMP END PARALLEL DO
    ELSE
      ! Inversion of a non-symmetric penta-diagonal matrix,
      !  provides Q*dver(t+dt).
      CALL MXTURE(NSMAX+1-KM,NFLEVG,NFLEVG,II,-2,.TRUE.,&
       & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
       & ZSDIVPL,ZSPDIVPL)  
      CALL MXTURE(NSMAX+1-KM,NFLEVG,NFLEVG,II,3,.FALSE.,&
       & SIHEG(1,IS0+1,1),SIHEG2(1,IS02+1,2),&
       & SIHEG2(1,IS02+1,3),ZSDIVPL,ZSPDIVPL)  
      ! Put the result in ZSPDIVP.
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,ISE)
      DO JN=KM,NSMAX
        ISE=KSTA+2*(JN-KM)
        ZSPDIVP(:,ISE:ISE+1)=ZSPDIVPL(:,JN,1:2)
      ENDDO
!$OMP END PARALLEL DO
    ENDIF

  ELSE

    ! Case designed when no stretching :

    IF (LIMPF) THEN
      ! Solve complex pentadiagonal system
      CALL SIMPLICO(KM,NSMAX,NFLEVG,NFLEVG,ZALPHA(KM),&
       & ZDENIM(KM),ZFPLUS(KM),ZFMINUS(KM),SIVP,YDLAP%RLAPDI(0:NSMAX),&
       & ZBDT2,ZSDIVP,ZSPDIVP)
    ELSE
      ! Inversion of a diagonal matrix, provides Q*dver(t+dt).
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
        DO JSP=KSTA,KEND
          DO JLEV=1,NFLEVG
          ZSPDIVP(JLEV,JSP)=ZSDIVP(JLEV,JSP)&
           & /(1.0_JPRB-ZBDT2*SIVP(JLEV)*YDLAP%RLAPDI(YDLAP%NVALUE(JSP+IOFF)))
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
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
  !   + Dprim_star_star " if LIMPF=F, or
  !   " beta**2 (Delta t)**2 vnabla'**2 (- SITR Gamma + C**2) d(t+dt)
  !   + Dprim_star_star_star " if LIMPF=T
  !   in array ZSRHS.
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV,IN)
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
      IN=YDLAP%NVALUE(JSP+IOFF)
      ZSRHS(JLEV,JSP)=(ZZSPSVDG(JLEV,JSP)*YDCST%RD*ZCC&
       & -ZWORK(JLEV,JSP))*SITR*ZBDT*ZBDT*&
       & YDLAP%RLAPDI(IN)+ZR1DPRIM(JLEV,JSP)  
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

  ! * Provides D'(t+dt).
  IF (LSIDG) THEN

    ! * Divide by the pentadiagonal operator
    !   [I - beta**2 (Delta t)**2 vnabla'**2 C**2 M**2]
    !   For KM>0 one works with the symmetric
    !   operator [vnabla'**(-2) - beta**2 (Delta t)**2 C**2 M**2]
    !   after having multiplied the RHS "ZSRHS" by "RLAPIN".
    !   Use ZSDIVPL and ZSPDIVPL as intermediate work arrays.

!$OMP WORKSHARE
    ZSDIVPL(:,:,:)=0.0_JPRB
    ZSPDIVPL(:,:,:)=0.0_JPRB
!$OMP END WORKSHARE
    IF (KM > 0) THEN
      ! Multiply the RHS by grad'**(-2).
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,ISE)
      DO JN=KM,NSMAX
        ISE=KSTA+2*(JN-KM)
        ZSDIVPL(:,JN,1:2)=ZSRHS(:,ISE:ISE+1)*YDLAP%RLAPIN(JN)
      ENDDO
!$OMP END PARALLEL DO
      ! Inversion of a symmetric penta-diagonal matrix.
      CALL MXTURS(NSMAX+1-KM,NFLEVG,NFLEVG,II,&
       & SIHEGB(1,IS0+1,1),SIHEGB(1,IS0+1,2),SIHEGB(1,IS0+1,3),&
       & ZSDIVPL,ZSPDIVPL)  
    ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,ISE)
      DO JN=KM,NSMAX
        ISE=KSTA+2*(JN-KM)
        ZSDIVPL(:,JN,1:2)=ZSRHS(:,ISE:ISE+1)
      ENDDO
!$OMP END PARALLEL DO
      ! Inversion of a non-symmetric penta-diagonal matrix.
      CALL MXTURE(NSMAX+1-KM,NFLEVG,NFLEVG,II,-2,.TRUE.,&
       & SIHEGB(1,IS0+1,1),SIHEGB(1,IS0+1,2),SIHEGB(1,IS0+1,3),&
       & ZSDIVPL,ZSPDIVPL)  
      CALL MXTURE(NSMAX+1-KM,NFLEVG,NFLEVG,II,3,.FALSE.,&
       & SIHEGB(1,IS0+1,1),SIHEGB2(1,IS02+1,2),&
       & SIHEGB2(1,IS02+1,3),ZSDIVPL,ZSPDIVPL)  
    ENDIF
    ! Put the result in ZZSPDIVG.
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,ISE)
    DO JN=KM,NSMAX
      ISE=KSTA+2*(JN-KM)
      ZZSPDIVG(:,ISE:ISE+1)=ZSPDIVPL(:,JN,1:2)
    ENDDO
!$OMP END PARALLEL DO

  ELSE

    ! * Retrieve D'(t+dt).
    IF (LIMPF) THEN
      ! * Solve complex pentadiagonal system.
      CALL SIMPLICO(KM,NSMAX,NFLEVG,NFLEVG,ZALPHA(KM),&
       & ZDENIM(KM),ZFPLUS(KM),ZFMINUS(KM),ZSIVP2,YDLAP%RLAPDI(0:NSMAX),&
       & ZBDT2,ZSRHS,ZZSPDIVG)
    ELSE
      ! * Division by the diagonal operator
      !   [I - beta**2 (Delta t)**2 vnabla'**2 C**2 RSTRET**2].
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV,IN)
      DO JLEV=1,NFLEVG
        DO JSP=KSTA,KEND
          IN=YDLAP%NVALUE(JSP+IOFF)
          ZZSPDIVG(JLEV,JSP)=ZSRHS(JLEV,JSP)/(1.0_JPRB-ZBTCM*YDLAP%RLAPDI(IN))
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    ENDIF

  ENDIF

  ! * Provides Mbar**2 D'(t+dt).
  IF (LSIDG) THEN
    ! Multiply ZZSPDIVG by M**2 (penta-diagonal operator).
    ! Use ZSDIVPL and ZSPDIVPL as intermediate work arrays.
!$OMP WORKSHARE
    ZSDIVPL(:,:,:)=0.0_JPRB
    ZSPDIVPL(:,:,:)=0.0_JPRB
!$OMP END WORKSHARE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,ISE)
    DO JN=KM,NSMAX
      ISE=KSTA+2*(JN-KM)
      ZSDIVPL(:,JN,1:2)=ZZSPDIVG(:,ISE:ISE+1)
    ENDDO
!$OMP END PARALLEL DO
    CALL MXPTMA(NSMAX+1-KM,NFLEVG,NFLEVG,II,SCGMAP(IS0+1,1),&
     & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
     & SCGMAP(IS0+1,2),SCGMAP(IS0+1,3),&
     & ZSDIVPL,ZSPDIVPL)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,ISE)
    DO JN=KM,NSMAX
      ISE=KSTA+2*(JN-KM)
      ZWORK(:,ISE:ISE+1)=ZSPDIVPL(:,JN,1:2)
    ENDDO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
      DO JSP=KSTA,KEND
        DO JLEV=1,NFLEVG
        ZWORK(JLEV,JSP)=ZZSPDIVG(JLEV,JSP)*RSTRET*RSTRET
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDIF

  !*       2.6  Increment vorticity for case LIMPF=T
  !             (identical to part 2.6 of SPCSI).

    IF (LIMPF) THEN
      IF (KM == 0) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,JLEV,IRSP)
        DO JN=2,NSMAX
          IRSP=KSTA+(JN-KM)*2
          DO JLEV=1,NFLEVG
            ZZSPVORG(JLEV,IRSP)=ZZSPVORG(JLEV,IRSP)&
             & -ZDENIM(JN)*ZFMINUS(JN)*ZZSPDIVG(JLEV,IRSP-2)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,JLEV,IRSP)
        DO JN=1,NSMAX-1
          IRSP=KSTA+(JN-KM)*2
          DO JLEV=1,NFLEVG
            ZZSPVORG(JLEV,IRSP)=ZZSPVORG(JLEV,IRSP)&
             & -ZDENIM(JN)*ZFPLUS(JN)*ZZSPDIVG(JLEV,IRSP+2)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,JLEV,IRSP,IMSP)
        DO JN=KM+1,NSMAX
          IRSP=KSTA+(JN-KM)*2
          IMSP=IRSP+1
          DO JLEV=1,NFLEVG
          ZZSPVORG(JLEV,IRSP)=ZZSPVORG(JLEV,IRSP)&
           & -ZDENIM(JN)*ZFMINUS(JN)*(ZZSPDIVG(JLEV,IRSP-2)&
           & -ZALPHA(JN)*ZZSPDIVG(JLEV,IMSP-2))
          ZZSPVORG(JLEV,IMSP)=ZZSPVORG(JLEV,IMSP)&
           & -ZDENIM(JN)*ZFMINUS(JN)*(ZZSPDIVG(JLEV,IMSP-2)&
           & +ZALPHA(JN)*ZZSPDIVG(JLEV,IRSP-2))
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JN,JLEV,IRSP,IMSP)
        DO JN=KM,NSMAX-1
          IRSP=KSTA+(JN-KM)*2
          IMSP=IRSP+1
          DO JLEV=1,NFLEVG
          ZZSPVORG(JLEV,IRSP)=ZZSPVORG(JLEV,IRSP)&
           & -ZDENIM(JN)*ZFPLUS(JN)*(ZZSPDIVG(JLEV,IRSP+2)&
           & -ZALPHA(JN)*ZZSPDIVG(JLEV,IMSP+2))
          ZZSPVORG(JLEV,IMSP)=ZZSPVORG(JLEV,IMSP)&
           & -ZDENIM(JN)*ZFPLUS(JN)*(ZZSPDIVG(JLEV,IMSP+2)&
           & +ZALPHA(JN)*ZZSPDIVG(JLEV,IRSP+2))
        ENDDO
       ENDDO
!$OMP END PARALLEL DO
      ENDIF
    ENDIF

  !*       2.7  Savings.

  IF (I_NITERHELM > 0) THEN
    IF (JITER == 0) THEN
      ! Save D'(t+dt,iter=0), Mbar**2 D'(t+dt,iter=0),
      !  zeta'(t+dt,iter=0), dver(t+dt,iter=0).
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JSP)
      DO JLEV=1,NFLEVG
        DO JSP=KSTA,KEND
        ZZSPDIV0G(JLEV,JSP)=ZZSPDIVG(JLEV,JSP)
        ZZSPGDIV0G(JLEV,JSP)=ZWORK(JLEV,JSP)
        ZZSPSVD0G(JLEV,JSP)=ZZSPSVDG(JLEV,JSP)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      IF (LIMPF) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JSP)
        DO JLEV=1,NFLEVG
          DO JSP=KSTA,KEND
            ZZSPVOR0G(JLEV,JSP)=ZZSPVORG(JLEV,JSP)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDIF
    ENDIF

    ! Save Mbar**2 D'(t+dt,jiter) in ZZSPGDIVG
    IF (JITER == 0) THEN
!$OMP WORKSHARE
      ZZSPGDIVG(1:NFLEVG,KSTA:KEND)=ZWORK(1:NFLEVG,KSTA:KEND)
!$OMP END WORKSHARE
    ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JSP)
      DO JLEV=1,NFLEVG
        DO JSP=KSTA,KEND
          ZZSPGDIVG(JLEV,JSP)=ZWORK(JLEV,JSP)+ZZSPGDIV0G(JLEV,JSP)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
    ENDIF
  ENDIF

ENDDO ! End loop on JITER

!*       2.8  Add the increment and recover X(t+dt)
!             for X=D',Mbar**2 D',zeta',dver, and deallocate.

IF (I_NITERHELM > 0) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JSP)
  DO JLEV=1,NFLEVG
    DO JSP=KSTA,KEND
      ZZSPDIVG(JLEV,JSP)=ZZSPDIVG(JLEV,JSP)+ZZSPDIV0G(JLEV,JSP)
      ZWORK(JLEV,JSP)=ZWORK(JLEV,JSP)+ZZSPGDIV0G(JLEV,JSP)
      ZZSPSVDG(JLEV,JSP)=ZZSPSVDG(JLEV,JSP)+ZZSPSVD0G(JLEV,JSP)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  IF (LIMPF) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JSP)
    DO JLEV=1,NFLEVG
      DO JSP=KSTA,KEND
        ZZSPVORG(JLEV,JSP)=ZZSPVORG(JLEV,JSP)+ZZSPVOR0G(JLEV,JSP)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDIF

  DEALLOCATE(ZZSPGDIVG)
  DEALLOCATE(ZZSPDIV0G)
  DEALLOCATE(ZZSPGDIV0G)
  DEALLOCATE(ZZSPVOR0G)
  DEALLOCATE(ZZSPSVD0G)

ENDIF

!*       2.9  Increment T, log(prehyds) and spd.

! * Provides some intermediate quantities allowing to compute
!   T(t+dt) (in ZST), log(prehyds)(t+dt) (in ZSP), spd(t+dt) (in ZSNHP).
CALL SIPTP(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZWORK,ZZSPSVDG,ZSNHP,ZST,ZSP,ISPCOL)

! * Provides T(t+dt) (in ZZSPTG) and spd(t+dt) (in ZZSPSPDG).
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
    ZZSPTG(JLEV,JSP)=ZZSPTG(JLEV,JSP)-ZBDT*ZST(JLEV,JSP)
    ZZSPSPDG(JLEV,JSP)=ZZSPSPDG(JLEV,JSP)-ZBDT*ZSNHP(JLEV,JSP)
  ENDDO
ENDDO
!$OMP END PARALLEL DO

! * Provides log(prehyds)(t+dt) (in ZZSPSPG).
DO JSP=KSTA,KEND
  ZZSPSPG(JSP)=ZZSPSPG(JSP)-ZBDT*ZSP(JSP)
ENDDO

!     ------------------------------------------------------------------

!*       3.    MEMORY TRANSFER AND DDH-SI UPDATE.
!              ----------------------------------
IF (LIMPF) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      PSPVORG(JLEV,JSP)=ZZSPVORG(JLEV,JSP)  
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDIF
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
DO JSP=KSTA,KEND
  DO JLEV=1,NFLEVG
    PSPDIVG(JLEV,JSP)=ZZSPDIVG(JLEV,JSP)
    PSPTG  (JLEV,JSP)=ZZSPTG  (JLEV,JSP)
    PSPSPDG(JLEV,JSP)=ZZSPSPDG(JLEV,JSP)
    PSPSVDG(JLEV,JSP)=ZZSPSVDG(JLEV,JSP)
  ENDDO
  PSPSPG (JSP)=ZZSPSPG (JSP)
ENDDO
!$OMP END PARALLEL DO

IF (LRSIDDH) THEN
  IF (LIMPF) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
        PSPTNDSI_VORG(JLEV,JSP)=PSPTNDSI_VORG(JLEV,JSP) + PSPVORG(JLEV,JSP)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDIF
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEV)
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      PSPTNDSI_DIVG(JLEV,JSP)=PSPTNDSI_DIVG(JLEV,JSP) + PSPDIVG(JLEV,JSP)
      PSPTNDSI_TG  (JLEV,JSP)=PSPTNDSI_TG  (JLEV,JSP) + PSPTG  (JLEV,JSP)
      PSPTNDSI_SPDG(JLEV,JSP)=PSPTNDSI_SPDG(JLEV,JSP) + PSPSPDG(JLEV,JSP)
      PSPTNDSI_SVDG(JLEV,JSP)=PSPTNDSI_SVDG(JLEV,JSP) + PSPSVDG(JLEV,JSP)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('SPNHSI',1,ZHOOK_HANDLE)
END SUBROUTINE SPNHSI
