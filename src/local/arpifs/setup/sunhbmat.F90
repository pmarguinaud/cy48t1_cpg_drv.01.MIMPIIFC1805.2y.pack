SUBROUTINE SUNHBMAT(YDCST, YDGEOMETRY,YDDYN)

!**** *SUNHBMAT*  Initialize vertical structure matrix for NHEE model.

!     Purpose.
!     --------
!      Initialize vertical structure matrix BB used for semi-implicit
!      scheme in the NHEE model.
!      This matrix is stored in the array SIB.

!      BB = C**2
!       * [ I - beta**2 (Delta t)**2 C**2 (SITR/SITRAM) (LLstar/H**2) ]**(-1)
!       * [ I + beta**2 (Delta t)**2 N**2 (SITR/SITRAM) TTstar ]
!       = C**2
!       * [ I - beta**2 (Delta t)**2 C**2 (LLstarstar/H**2) ]**(-1)
!       * [ I + beta**2 (Delta t)**2 N**2 TTstarstar ]
!      where:
!       C = sqrt( Rd * SITR * ( cpdry/cvdry ) )
!       H = Rd * SITR / g
!       N = g / sqrt( cpdry * SITR )
!       LLstar is a tridiagonal "Laplacian" operator.
!       LLstarstar=(SITR/SITRAM) LLstar.
!       When constraint C2 is matched,
!        TTstar is a tridiagonal operator equal to (I + LLstar * QQstar)
!         (QQstar is the diagonal matrix "deltastar - 2*alphastar")
!        TTstarstar=(SITR/SITRAM) TTstar writes (SITR/SITRAM + LLstarstar * QQstar)
!       When constraint C2 is not matched, TTstarstar (computed by SI_CCCOR),
!         has a more tricky expression involving LLstarstar, gamma, tau.
!         TTstarstar is not tridiagonal in this case.

!      The intermediate factor:
!       [ I - beta**2 (Delta t)**2 C**2 (LLstarstar/H**2) ]
!       is stored in SIFAC.

!      The intermediate factor:
!       [ I - beta**2 (Delta t)**2 C**2 (LLstarstar/H**2) ]**(-1)
!       is stored in SIFACI. It is required in the RHS of the Helmholtz equation.

!      The LHS of the Helmholtz equation writes:
!       (I - beta**2 (Delta t)**2 BB Mbabar2 vnablaprim**2) svd
!      where
!       Mbabar2 = c**2 if LSIDG=F
!       Mbabar2 = M**2 if LSIDG=T (in this case a pentadiagonal operator
!        must be computed in the routine SUNHHEG).
!       "svd" is the vertical divergence variable.

!      Replace "beta" by "(1+epsilon_uncentering) * beta"
!       in the previous formulae
!       when there is uncentering (in the SL scheme).

!     Remark: the upper boundary condition is well specified at least
!      for models with "prehyd_top = 0" but it is not sure that it
!      is well specied if "prehyd_top > 0" (more generally it is
!      not currently desirable to run the NH model with "prehyd_top > 0").
!      At least for prehyd_top = 0, QQstar(1,1) must be set to zero.

!**   Interface.
!     ----------
!        *CALL* *SUNHBMAT

!     Explicit arguments :  None
!     --------------------

!     Implicit arguments :
!     --------------------
!        Matrix SIB, SIFAC, SIFACI in YOMDYN

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
!        K. YESSAD (after routines SUBMAT and SUEHEL).

!     Modifications.
!     --------------
!   Original : 09-12-2004
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Feb 2012): use RBTS2.
!   P. Smolikova (Nov 2013): bugfix to compute TTstar when constraint C2 is not matched.
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Vertical-dependent SITRA.
!   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMDYNA  , ONLY : YRDYNA
USE YOMDYN   , ONLY : TDYN
USE YOMCST   , ONLY : TCST
USE YOMCVER  , ONLY : LVERTFE, NVFE_ORDER, LVFE_LAPL_BC
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : LOUTPUT, NPRINTLEV
USE YOMVERT  , ONLY : VP00,TOPPRES

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT):: YDDYN

INTEGER(KIND=JPIM) :: JL1, JL2, IFIRST

REAL(KIND=JPRB) :: ZBT2
REAL(KIND=JPRB) :: ZCOEF
REAL(KIND=JPRB) :: ZC2, ZH, ZH2, ZN2
REAL(KIND=JPRB) :: ZD(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSIB(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSIBIN(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSIBOUT(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZCOR(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZX, ZA, ZMAE, ZRMSE
REAL(KIND=JPRB) :: ZFUNC(1:YDGEOMETRY%YRDIMV%NFLEVG),ZANA(1:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZVFE(1:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZRATIO_SITR(YDGEOMETRY%YRDIMV%NFLEVG)
LOGICAL :: LLNH_NOC1C2

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "mxmaop.h"
#include "minv_caller.h"
#include "si_cccor.intfb.h"
#include "siseve.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUNHBMAT',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,   RBTS2=>YDDYN%RBTS2,   SIALPH=>YDDYN%SIALPH, SIB=>YDDYN%SIB,     &
& SIFAC=>YDDYN%SIFAC, SIFACI=>YDDYN%SIFACI,   SILNPR=>YDDYN%SILNPR, SITIME=>YDDYN%SITIME, SITR=>YDDYN%SITR,   &
& SITRAM=>YDDYN%SITRAM)
!     ------------------------------------------------------------------

!*       1.    INITIALIZE.
!              -----------

! * ZD contains the identity matrix "I":
DO JL1=1,NFLEVG
  DO JL2=1,NFLEVG
    ZD(JL2,JL1)=1.0_JPRB-MIN(ABS(REAL(JL2-JL1,JPRB)),1.0_JPRB)
  ENDDO
ENDDO

! * Compute the constants C**2, H**2 and N**2:
ZC2=YDCST%RD*SITR*YDCST%RCPD/YDCST%RCVD
ZH=YDCST%RD*SITR/YDCST%RG
ZH2=ZH*ZH
ZN2=(YDCST%RG*YDCST%RG)/(YDCST%RCPD*SITR)

! * Compute ZBT2 = (1+epsilon_uncentering)**2 * BETADT**2 * TSTEP**2:
ZBT2=(RBTS2*SITIME)**2

! * In the case where prehyd_top=0, a particular treatment
!   is applied to some quantities at the first full level;
!   for ex QQstar(1,1)=0
!   instead of
!   QQstar(1,1)=deltastar(1)-2*alphastar(1)
IF(YDVAB%VAH(0)+YDVAB%VBH(0)*VP00 <= TOPPRES)THEN
  IFIRST=2
ELSE
  IFIRST=1
ENDIF

!      ---------------------------------------------------------------

!*       2.    COMPUTES THE MATRIX SIFAC.
!              --------------------------

! * Apply the tridiagonal "Laplacian" operator LLstarstar to ZD.
CALL SISEVE(LVERTFE, LVFE_LAPL_BC, YRDYNA,YDGEOMETRY,YDDYN,1,NFLEVG,ZD,SIFAC,NFLEVG)

! * Multiply by
! "-(1+epsilon_uncentering)**2 beta**2 (Delta t)**2 C**2 (1/H**2)"
! then add the identity matrix:
ZCOEF=ZBT2*ZC2/ZH2
DO JL1=1,NFLEVG
  DO JL2=1,NFLEVG
    SIFAC(JL2,JL1)=ZD(JL2,JL1)-ZCOEF*SIFAC(JL2,JL1)
  ENDDO
ENDDO

!      ---------------------------------------------------------------

!*       3.    INVERTS THE MATRIX SIFAC.
!              --------------------------

CALL MINV_CALLER(.TRUE.,NFLEVG,SIFAC,SIFACI)

!      ---------------------------------------------------------------

!*       4.    COMPUTES THE MATRIX SIB.
!              ------------------------

LLNH_NOC1C2=.NOT.(YRDYNA%NDLNPR == 1)

! * Calculation of TTstar:
IF (LLNH_NOC1C2) THEN

  ! Constraint C2 is not matched, general formula for TTstar (= (g**2/(N**2 C**2)) LLstar* A2star)
  CALL SI_CCCOR(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,NFLEVG,ZD,ZCOR,POU2=ZSIBOUT)

ELSE

  ! Constraint C2 is matched, tridiagonal TTstar (=I + LLstar * QQstar)

  ! * Apply the diagonal operator QQstar:
  ZSIBIN(1:NFLEVG,1:NFLEVG)=0.0_JPRB
  DO JL1=IFIRST,NFLEVG
    ZSIBIN(JL1,JL1)=SILNPR(JL1)-2.0_JPRB*SIALPH(JL1)
  ENDDO

  ! * Apply the tridiagonal operator LLstar:
  CALL SISEVE(LVERTFE, LVFE_LAPL_BC, YRDYNA,YDGEOMETRY,YDDYN,1,NFLEVG,ZSIBIN,ZSIBOUT,NFLEVG)

  ! * Add (SITR/SITRAM)*Id:
  DO JL1=1,NFLEVG
    ZRATIO_SITR(JL1)=SITR/SITRAM(JL1)
    ZSIBOUT(JL1,JL1)=ZRATIO_SITR(JL1)+ZSIBOUT(JL1,JL1)
  ENDDO

ENDIF

! * Multiply by
! "(1+epsilon_uncentering)**2 beta**2 (Delta t)**2 N**2",
! then adds the identity matrix,
! then multiply by C**2:
ZCOEF=ZBT2*ZN2
DO JL1=1,NFLEVG
  DO JL2=1,NFLEVG
    ZSIBOUT(JL2,JL1)=ZC2*(ZD(JL2,JL1)+ZCOEF*ZSIBOUT(JL2,JL1))
  ENDDO
ENDDO

! * Multiply by SIFACI:
CALL MXMAOP(SIFACI,1,NFLEVG,ZSIBOUT,1,NFLEVG,ZSIB,1,NFLEVG,&
 & NFLEVG,NFLEVG,NFLEVG)

!      ----------------------------------------------------------------

!*       5.    STORAGE OF MATRIX.
!              ------------------

DO JL1=1,NFLEVG
  DO JL2=1,NFLEVG
    SIB(JL2,JL1)=ZSIB(JL2,JL1)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*       6.   BASIC TEST OF VERTICAL LAPLACIAN
!             --------------------------------

IF (LVERTFE .AND. LOUTPUT .AND. NPRINTLEV>=1) THEN
  ! acuracy of laplacian operator against analytical function
  DO JL1=1,NFLEVG
    ZX = YDVETA%VFE_ETAF(JL1)
    ZA = YDCST%RPI*ZX

    ! function
    ZFUNC(JL1) = (SIN(ZA)**3)*COS(ZA)

    ! analytical solution
    ZANA(JL1)  = ZX*(6.0_JPRB*YDCST%RPI*COS(ZA)**2*SIN(ZA)**2 &
     & - 2.0_JPRB*YDCST%RPI*SIN(ZA)**4 + ZX*(6.0_JPRB*YDCST%RPI**2*COS(ZA)**3*SIN(ZA) &
     & - 10.0_JPRB*YDCST%RPI**2*COS(ZA)*SIN(ZA)**3))
  ENDDO

  CALL SISEVE(LVERTFE, LVFE_LAPL_BC, YRDYNA,YDGEOMETRY,YDDYN,1,NFLEVG,ZFUNC,ZVFE,1)

  DO JL1=1,NFLEVG
    WRITE(NULOUT,'("(LAPL ANALYSIS) :",I5,1X,4F15.12)') JL1, ZFUNC(JL1), &
     & ZVFE(JL1), ZANA(JL1), ZVFE(JL1) - ZANA(JL1)
  ENDDO

  ZMAE  = 0.0_JPRB
  ZRMSE = 0.0_JPRB
  DO JL1=1,NFLEVG
    ZMAE  = ZMAE  + ABS(ZVFE(JL1) - ZANA(JL1))
    ZRMSE = ZRMSE + (ZVFE(JL1) - ZANA(JL1))**2
  ENDDO
  ZMAE  = ZMAE/REAL(NFLEVG,JPRB)
  ZRMSE = SQRT(ZRMSE/REAL(NFLEVG,JPRB))

  WRITE(NULOUT,'("LAPLACIAN ", I4, I4, " MAE: ",F15.12," RMSE: ",F15.12)') &
   & NFLEVG, NVFE_ORDER, ZMAE, ZRMSE
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUNHBMAT',1,ZHOOK_HANDLE)
END SUBROUTINE SUNHBMAT
