SUBROUTINE SUNHSI_TESTCONV(YDCST, YDGEOMETRY,YDRIP,YDDYN)

!**** *SUNHSI_TESTCONV*   - Convergence test for the iterative Helmholtz
!                  solver in case of unsatisfied condition (C1).

!     Purpose.    Calculates eigenvalues of an iteration matrix
!                 in case the iterative solution of the Helmholz equation 
!                 is applied in the semi-implicit (for ex. VFE NHEE).
!
!       Calculates the convergence matrix ZMCONV as:
!          ZMCONV = ZMA^{-1} (Id + ZMB ZMH^{-1} ZMF) ZMC, where
!          ZMA = 1-beta delta(t)^2 C^2 nabla^2 m^2,
!          ZMB = beta delta(t)^2 nabla^2 (-RT* G* + C^2),
!          ZMC = - beta delta(t)^2 nabla^2 m^2 COR, where
!             COR = - RT* (G*S* - S* - G* + N*) calculated in SI_CCCOR,
!          ZMF = beta delta(t)^2 m^2 (T*/Te*) (1/H^2) Lv* (-RT* S* + C^2),
!          ZMH = 1-beta delta(t)^2 m^2 SIB nabla^2.
!          ZMHINV = ZMH^{-1} is calculated using factorization as 
!             SIMO (1-beta delta(t)^2 m^2 SIVP nabla^2)^(-1) SIMI SIFACI.
!       Determines eigenvalues of matrix ZMCONV and tests if 
!       |lambda(MCONV)|<1 for all eigenvalues lambda of ZMCONV.
!       If YES => the iterative process converges,
!       if NO  => the iterative process diverges.
! ---------------------------------------------------------

!**   Interface.
!     ----------
!        *CALL* *SUNHSI_TESTCONV

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE/ALADIN documentation

!     Author.
!     -------
!      P. Smolikova and J. Vivoda (2013)

!     Modifications.
!     --------------
!      K. Yessad (July 2014): Move some variables.
!      K. Yessad (June 2017): Vertical-dependent SITRA.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCT0   , ONLY : LELAM    ,LTWOTL
USE YOMDYN   , ONLY : TDYN
USE YOMRIP   , ONLY : TRIP
USE YOMLUN   , ONLY : NULOUT   ,NULERR
USE YOMCST   , ONLY : TCST
USE YOMDYNA  , ONLY : YRDYNA
USE YOMCVER      , ONLY : LVERTFE, LVFE_LAPL_BC

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT) :: YDDYN
TYPE(TRIP)     ,INTENT(INOUT) :: YDRIP
REAL(KIND=JPRB) :: ZID  (YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZCOR (YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZTDT, ZBT2, ZC2, ZM2, ZN2, ZRH2
REAL(KIND=JPRB) :: ZINV, ZAUX, ZMAXFN
REAL(KIND=JPRB) :: ZHOOK_HANDLE

REAL(KIND=JPRB) :: ZGAM(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTAU(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZNU(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSP(YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZMB(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZMC(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZMF(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZMBG(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZMFG(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZMHINV(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZMAUX1(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZMAUX2(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZMAUX3(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZMAUX4(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZMAUX5(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZMCONV(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)

! Arrays for eigenvalues calc.
REAL(KIND=JPRB) :: ZFR(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZFI(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZFN(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZMO(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZWO(YDGEOMETRY%YRDIMV%NFLEVG+1)
INTEGER(KIND=JPIM) :: IWO(YDGEOMETRY%YRDIMV%NFLEVG+1)
INTEGER(KIND=JPIM) :: IER

INTEGER(KIND=JPIM) :: IM, IN, JN, JL1, JL2, JMLOC

!      ----------------------------------------------------------------

#include "eigsol.h"
#include "siseve.intfb.h"
#include "sigam.intfb.h"
#include "sitnu.intfb.h"
#include "si_cccor.intfb.h"

!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUNHSI_TESTCONV',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDLAP=>YDGEOMETRY%YRLAP, &
& YDELAP=>YDGEOMETRY%YRELAP)
ASSOCIATE(NSMAX=>YDDIM%NSMAX, NUMP=>YDDIM%NUMP,   NFLEVG=>YDDIMV%NFLEVG,   RBTS2=>YDDYN%RBTS2, SIFACI=>YDDYN%SIFACI, &
& SIMI=>YDDYN%SIMI, SIMO=>YDDYN%SIMO,   SITR=>YDDYN%SITR, SIVP=>YDDYN%SIVP,   RSTRET=>YDGEM%RSTRET,                  &
& MYMS=>YDLAP%MYMS, RLAPDI=>YDLAP%RLAPDI,   TSTEP=>YDRIP%TSTEP)
!      ----------------------------------------------------------------

!*     1.    PRELIMINARY CALCULATIONS.
!      -------------------------------

IF (LTWOTL) THEN
  ZTDT=TSTEP
ELSE
  ZTDT=2.0_JPRB*TSTEP
ENDIF

WRITE(NULOUT,'('' SUNHSI_TESTCONV: COMPUTES EIGENVALUES OF THE ITERATION MATRIX MCONV'')')
ZBT2 = (RBTS2*ZTDT)**2
ZM2  = RSTRET*RSTRET
ZC2  = YDCST%RD*SITR*YDCST%RCPD/YDCST%RCVD
ZN2  = YDCST%RG*YDCST%RG/(YDCST%RCPD*SITR)
ZRH2 = YDCST%RG*YDCST%RG/(YDCST%RD*YDCST%RD*SITR*SITR)

!      ----------------------------------------------------------------

!*     2.    AUXILIARY MATRICES CALCULATION
!      -------------------------------------

! identity matrix computation
DO JL1=1,NFLEVG
  DO JL2=1,NFLEVG
    ZID(JL2,JL1)=1.0_JPRB-MIN(ABS(REAL(JL2-JL1,JPRB)),1.0_JPRB)
  ENDDO
ENDDO

ZSP(1:NFLEVG)=0.0_JPRB
CALL SIGAM(LVERTFE, YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZGAM,ZID,ZSP,NFLEVG,NFLEVG)
CALL SITNU(LVERTFE, YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZID,ZTAU,ZNU,NFLEVG)

! matrices independent of spectral wavenumbers
DO JL1=1,NFLEVG
  DO JL2=1,NFLEVG
    ZMBG(JL1,JL2) = ZC2*ZID(JL1,JL2)-SITR*ZGAM(JL1,JL2)
    ZMFG(JL1,JL2) = ZBT2*ZRH2*ZM2*&
     & (ZC2*ZID(JL1,JL2)-YDCST%RCPD*ZTAU(JL1,JL2))
  ENDDO
ENDDO
CALL SISEVE(LVERTFE, LVFE_LAPL_BC, YRDYNA,YDGEOMETRY,YDDYN,1,NFLEVG,ZMFG,ZMF,NFLEVG)
CALL SI_CCCOR(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,NFLEVG,ZID,ZCOR)
CALL MXMAOP(SIMI,1,NFLEVG,SIFACI,1,NFLEVG,ZMAUX1,1,NFLEVG,NFLEVG,NFLEVG,NFLEVG)

!      ----------------------------------------------------------------

!*     3.    DEPENDENCY ON SPECTRAL WAVENUMBERS
!      ----------------------------------------

ZMAXFN=0.0_JPRB

IF (LELAM) THEN

  DO JMLOC=1,NUMP
    IM=MYMS(JMLOC)
    DO JN=0,YDELAP%NCPLM(IM)-1
      IN=YDELAP%NPME(IM)+JN

      ZAUX=ZBT2*YDELAP%RLEPDIM(IN)

      DO JL1=1,NFLEVG
        DO JL2=1,NFLEVG
          ZMB(JL1,JL2) = ZAUX*ZMBG(JL1,JL2)
          ZMC(JL1,JL2) = -ZAUX*ZM2*ZCOR(JL1,JL2)
        ENDDO
      ENDDO

      ! 3.1 Define ZMCONV
      ! -----------------

      DO JL1=1,NFLEVG
        DO JL2=1,NFLEVG
          ZMAUX2(JL1,JL2)=ZMAUX1(JL1,JL2)/(1.0_JPRB-ZAUX*ZM2*SIVP(JL1))
        ENDDO
      ENDDO

      CALL MXMAOP(SIMO,1,NFLEVG,ZMAUX2,1,NFLEVG,ZMHINV,1,NFLEVG,NFLEVG,NFLEVG,NFLEVG)
      CALL MXMAOP(ZMHINV,1,NFLEVG,ZMF,1,NFLEVG,ZMAUX3,1,NFLEVG,NFLEVG,NFLEVG,NFLEVG)
      CALL MXMAOP(ZMB,1,NFLEVG,ZMAUX3,1,NFLEVG,ZMAUX4,1,NFLEVG,NFLEVG,NFLEVG,NFLEVG)

      ZINV=1.0_JPRB/(1.0_JPRB-ZC2*ZM2*ZAUX)
      DO JL1=1,NFLEVG
        DO JL2=1,NFLEVG
          ZMAUX5(JL1,JL2)=ZINV*(ZMAUX4(JL1,JL2)+ZID(JL1,JL2))
        ENDDO
      ENDDO

      CALL MXMAOP(ZMAUX5,1,NFLEVG,ZMC,1,NFLEVG,ZMCONV,1,NFLEVG,NFLEVG,NFLEVG,NFLEVG)

      ! 3.2 Determine eigenvalues of ZMCONV and test them
      ! -------------------------------------------------
      CALL EIGSOL(NFLEVG,NFLEVG,ZMCONV,ZFR,ZFI,0,ZMO,IWO,ZWO,IER)

      DO JL1=1,NFLEVG
        ZFN(JL1)=SQRT(ZFR(JL1)*ZFR(JL1)+ZFI(JL1)*ZFI(JL1))
        IF (ZFN(JL1)>1.0_JPRB) THEN
          WRITE(NULOUT,*) 'SUNHSI_TESTCONV: Iterative Helmholtz solver alert!!!'
          WRITE(NULOUT,*) 'JL1,ZFR,ZFI,ZFN=',JL1,ZFR(JL1),ZFI(JL1),ZFN(JL1)
          WRITE(NULOUT,*) 'Convergence not garanteed |eigval(M)|>1'
          WRITE(NULERR,*) 'SUNHSI_TESTCONV: Iterative Helmholtz solver alert!!!'
          WRITE(NULERR,*) 'JL1,ZFR,ZFI,ZFN=',JL1,ZFR(JL1),ZFI(JL1),ZFN(JL1)
          WRITE(NULERR,*) 'Convergence not garanteed |eigval(M)|>1'
        ENDIF
      ENDDO

      ZMAXFN=MAX(ZMAXFN,MAXVAL(ZFN))

    ENDDO
  ENDDO

ELSE

  DO JMLOC=1,NUMP
    IM=MYMS(JMLOC)
    DO JN=IM,NSMAX
      IN=JN

      ZAUX=ZBT2*RLAPDI(IN)

      DO JL1=1,NFLEVG
        DO JL2=1,NFLEVG
          ZMB(JL1,JL2) = ZAUX*ZMBG(JL1,JL2)
          ZMC(JL1,JL2) = -ZAUX*ZM2*ZCOR(JL1,JL2)
        ENDDO
      ENDDO

      ! 3.1 Define ZMCONV
      ! -----------------

      DO JL1=1,NFLEVG
        DO JL2=1,NFLEVG
          ZMAUX2(JL1,JL2)=ZMAUX1(JL1,JL2)/(1.0_JPRB-ZAUX*ZM2*SIVP(JL1))
        ENDDO
      ENDDO

      CALL MXMAOP(SIMO,1,NFLEVG,ZMAUX2,1,NFLEVG,ZMHINV,1,NFLEVG,NFLEVG,NFLEVG,NFLEVG)
      CALL MXMAOP(ZMHINV,1,NFLEVG,ZMF,1,NFLEVG,ZMAUX3,1,NFLEVG,NFLEVG,NFLEVG,NFLEVG)
      CALL MXMAOP(ZMB,1,NFLEVG,ZMAUX3,1,NFLEVG,ZMAUX4,1,NFLEVG,NFLEVG,NFLEVG,NFLEVG)

      ZINV=1.0_JPRB/(1.0_JPRB-ZC2*ZM2*ZAUX)
      DO JL1=1,NFLEVG
        DO JL2=1,NFLEVG
          ZMAUX5(JL1,JL2)=ZINV*(ZMAUX4(JL1,JL2)+ZID(JL1,JL2))
        ENDDO
      ENDDO

      CALL MXMAOP(ZMAUX5,1,NFLEVG,ZMC,1,NFLEVG,ZMCONV,1,NFLEVG,NFLEVG,NFLEVG,NFLEVG)

      ! 3.2 Determine eigenvalues of ZMCONV and test them
      ! -------------------------------------------------
      CALL EIGSOL(NFLEVG,NFLEVG,ZMCONV,ZFR,ZFI,0,ZMO,IWO,ZWO,IER)

      DO JL1=1,NFLEVG
        ZFN(JL1)=SQRT(ZFR(JL1)*ZFR(JL1)+ZFI(JL1)*ZFI(JL1))
        IF (ZFN(JL1)>1.0_JPRB) THEN
          WRITE(NULOUT,*) 'SUNHSI_TESTCONV: Iterative Helmholtz solver alert!!!'
          WRITE(NULOUT,*) 'JL1,ZFR,ZFI,ZFN=',JL1,ZFR(JL1),ZFI(JL1),ZFN(JL1)
          WRITE(NULOUT,*) 'Convergence not garanteed |eigval(M)|>1'
          WRITE(NULERR,*) 'SUNHSI_TESTCONV: Iterative Helmholtz solver alert!!!'
          WRITE(NULERR,*) 'JL1,ZFR,ZFI,ZFN=',JL1,ZFR(JL1),ZFI(JL1),ZFN(JL1)
          WRITE(NULERR,*) 'Convergence not garanteed |eigval(M)|>1'
        ENDIF
      ENDDO

      ZMAXFN=MAX(ZMAXFN,MAXVAL(ZFN))

    ENDDO
  ENDDO

ENDIF

WRITE(NULOUT,*) 'SUNHSI_TESTCONV: Highest eigenvalue of MCONV for current processor:', ZMAXFN

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUNHSI_TESTCONV',1,ZHOOK_HANDLE)
END SUBROUTINE SUNHSI_TESTCONV
