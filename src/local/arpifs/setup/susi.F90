#ifdef VPP
!OCL SCALAR
#endif
SUBROUTINE SUSI(YDCST, YDGEOMETRY,YDRIP,YDDYN,YDEDYN,KULOUT)

!**** *SUSI*   - Initialize constants and control for the dynamics:
!                semi-implicit scheme operators.

!     Purpose.
!     --------
!        Computes the operator B=mu*nu+gamma*tau which is used in
!        the semi-implicit scheme in order to work on the Hough space.
!        The inverse 1/B and also the eigenvalues and the eigenvectors
!        of B are computed.
!        * B is stored in array SIB.
!        * The eigenvalues of B are stored in array SIVP.
!        * The eigenvectors of B are stored in array SIMO.
!        * The inverse of SIMO is stored in array SIMI.
!        * The inverse of B is stored in array SIBI.
!        * The call to SUHEG fills the arrays SIHEG and SIHEG2.
!        In the hydrostatic model, only the content of SIHEG and SIHEG2
!         depend on "Delta t".

!**   Interface.
!     ----------
!        *CALL* *SUSI(KULOUT)

!     Explicit arguments :
!     --------------------
!        KULOUT : Logical unit for the output

!     Implicit arguments :
!     --------------------

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
!        K. YESSAD (after old part 7 of SUDYN).

! Modifications
! -------------
!   Original : 08-Dec-2004
!   I. Santos, I. Martinez (Mar 2010): Variable map factor for RTM.
!   K. Yessad (Dec 2011): Use GPHPRE.
!   M. Fisher   7-March-2012 Move SIBI out of Jb (to allow late Jb setup)
!   J. Vivoda and P. Smolikova (Sep 2020): Add analysis of stability.
!   H Petithomme (Dec 2020): optimisation on gphpre
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCVER,ONLY: LVERTFE
USE YOMDYN   , ONLY : TDYN
USE YEMDYN   , ONLY : TEDYN
USE YOMCST   , ONLY : TCST       
USE YOMRIP   , ONLY : TRIP
USE INTDYN_MOD,ONLY : YYTXYB

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY), INTENT(INOUT)    :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT)    :: YDDYN
TYPE(TEDYN)    ,INTENT(INOUT)    :: YDEDYN
TYPE(TRIP)     ,INTENT(INOUT)    :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JPLL
PARAMETER(JPLL=1)

REAL(KIND=JPRB) :: ZPRES(JPLL,0:YDGEOMETRY%YRDIMV%NFLEVG),ZPRESF(JPLL,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZXYB(JPLL,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)
REAL(KIND=JPRB) :: ZWO(YDGEOMETRY%YRDIMV%NFLEVG+1),ZA(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZFR(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZFI(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZMO(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IWO(YDGEOMETRY%YRDIMV%NFLEVG+1)
INTEGER(KIND=JPIM) :: IER, ISIGN, ISTOP, JL1, JL2, JLEV, JLON, J

REAL(KIND=JPRB) :: ZSUM, ZTES

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "eigsol.h"
#include "minv_caller.h"
#include "gphpre.intfb.h"
#include "scordo.intfb.h"
#include "submat.intfb.h"
#include "sudyn_stability.intfb.h"
#include "sueheg.intfb.h"
#include "suheg.intfb.h"

!     ------------------------------------------------------------------


IF (LHOOK) CALL DR_HOOK('SUSI',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIMV=>YDGEOMETRY%YRDIMV,YDVAB=>YDGEOMETRY%YRVAB)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, LSIDG=>YDDYN%LSIDG, LDYN_STABAN=>YDDYN%LDYN_STABAN,   SIALPH=>YDDYN%SIALPH, &
& SIB=>YDDYN%SIB, SIBI=>YDDYN%SIBI,   SIDELP=>YDDYN%SIDELP, SIDPHI=>YDDYN%SIDPHI, SILNPR=>YDDYN%SILNPR,      &
& SIMI=>YDDYN%SIMI, SIMO=>YDDYN%SIMO, SIPR=>YDDYN%SIPR, SIRDEL=>YDDYN%SIRDEL,   SITIME=>YDDYN%SITIME,        &
& SITLAF=>YDDYN%SITLAF, SITLAH=>YDDYN%SITLAH,   SITR=>YDDYN%SITR, SIVP=>YDDYN%SIVP, VNORM=>YDDYN%VNORM,      &
& LESIDG=>YDEDYN%LESIDG)
!     ------------------------------------------------------------------

!*       1.    INITIALIZE SEMI-IMPLICIT AND VERTICAL NORMAL MODES
!              FOR THE 3D PRIMITIVE EQUATION MODEL.
!              --------------------------------------------------

!*       1.1   INITIALIZE SEMI-IMPLICIT COMPUTATION.

IER=0
ZTES=0.0_JPRB
ISTOP=0

DO JLON=1,JPLL
  ZPRES(JLON,NFLEVG)=SIPR
ENDDO

! alpha in SI only used in NH
CALL GPHPRE(JPLL,NFLEVG,JPLL,JPLL,YDVAB,ZPRES,PXYB=ZXYB,PRESF=ZPRESF,LALPHA=.NOT.LVERTFE,&
 & LRTGR=.FALSE.,LRPP=.FALSE.)

DO JLEV=1,NFLEVG
  SIDELP(JLEV)=ZXYB(1,JLEV,YYTXYB%M_DELP)
  SIRDEL(JLEV)=ZXYB(1,JLEV,YYTXYB%M_RDELP)
  SILNPR(JLEV)=ZXYB(1,JLEV,YYTXYB%M_LNPR)
  SIALPH(JLEV)=ZXYB(1,JLEV,YYTXYB%M_ALPH)
  SITLAH(JLEV)=ZPRES(1,JLEV)
ENDDO
DO JLEV=1,NFLEVG
  SITLAF(JLEV)=ZXYB(1,JLEV,YYTXYB%M_DELP)/ZXYB(1,JLEV,YYTXYB%M_LNPR)
  SIDPHI(JLEV)=YDCST%RD*SITR*ZXYB(1,JLEV,YYTXYB%M_LNPR)
ENDDO

DO JLEV=1,NFLEVG
  WRITE(UNIT=KULOUT,FMT='('' JLEV = '',I3.3,&
   & '' SITLAF = '',E14.8,'' SIDPHI = '',E14.8)')&
   & JLEV,SITLAF(JLEV),SIDPHI(JLEV)  
ENDDO

SITLAH(0)=ZPRES(1,0)

!*       1.2   INITIALIZE MATRIX B (SIB)

CALL SUBMAT(YDCST, YDGEOMETRY,YDDYN)
  
!*       1.3   SOLVE EIGEN PROBLEM

DO JL1=1,NFLEVG
  DO JL2=1,NFLEVG
    ZA(JL1,JL2)=SIB(JL1,JL2)
  ENDDO
ENDDO

CALL EIGSOL(NFLEVG,NFLEVG,ZA,ZFR,ZFI,1,ZMO,IWO,ZWO,IER)
IF(IER /= 0)THEN
  CALL ABOR1('SUSI: ABORT IN SUSI AFTER CALL TO RG')
ENDIF
CALL SCORDO(NFLEVG,ZFR,ZFI,ZMO,ZA)

DO JL1=1,NFLEVG
  DO JL2=1,NFLEVG
    SIMO(JL1,NFLEVG+1-JL2)=ZMO(JL1,JL2)
  ENDDO
ENDDO

DO JL1=1,NFLEVG
  SIVP(NFLEVG+1-JL1)=ZFR(JL1)
ENDDO

ZTES=MAXVAL(ZFI)

!        1.4   NORMALIZE THE VERTICAL MODES
!                               --
! NEW SCALING     , > X(J)**2 *DP(J)/PSOL = VNORM
!                               --
!  (PRESERVING THE AMPLITUDE OF THE EIGENMODE AS THE VERTICAL RESOLUTION
!   -NFLEVG- CHANGES , C.TEMPERTON  09/91)

DO JL1=1,NFLEVG
  ZSUM = 0._JPRB
  DO JL2=1,NFLEVG
    ZSUM=ZSUM+SIMO(JL2,JL1)*SIMO(JL2,JL1)*SIDELP(JL2)
  ENDDO
  ZSUM = ZSUM/SIPR/VNORM
  ZSUM = SQRT(ZSUM)
  DO JL2=1,NFLEVG
    SIMO(JL2,JL1) = SIMO(JL2,JL1)/ZSUM
  ENDDO
ENDDO

!        1.5   CHECK SIGNS OF VERTICAL MODES

DO JL1=1,NFLEVG
  ISIGN = +1
  IF (SIMO(NFLEVG,JL1) < 0.0_JPRB) ISIGN = -ISIGN
  IF (ISIGN == -1) THEN
    DO JL2=1,NFLEVG
      SIMO(JL2,JL1) = -SIMO(JL2,JL1)
    ENDDO
  ENDIF
ENDDO

!*       1.6   INVERSE MODE MATRIX

CALL MINV_CALLER(.TRUE.,NFLEVG,SIMO,SIMI)

!*       1.7   PRINT OUT INFORMATION ON VERTICAL MODES

WRITE(KULOUT,481) SITR,SIPR,(J,YDVAB%VAH(J-1),YDVAB%VBH(J-1),SIVP(J),J=1,NFLEVG)
481 FORMAT('1VERTICAL MODES FOR SITR=',F6.1,' K;   SIPR=',F8.0,&
 & ' PA'///,' VERTICAL RESOLUTION         A(K)            B(K)',&
 & 26X,'EIGENVALUES '/&
 & (20X,I4,G13.5,4X,G13.5,18X,G13.5))  
WRITE(KULOUT,482) ZTES
482 FORMAT(///'0LARGEST IMAGINARY PART OF EIGENVALUES:',G13.5)
IF (ZTES /= 0.0_JPRB) THEN
  ! * Print the imaginary part of eigenvalues.
  WRITE(UNIT=KULOUT,FMT='('' IMAGINARY PARTS OF EIGENVALUES'')')
  WRITE(KULOUT,FMT='(14(1X,E8.2))')(ZFI(NFLEVG+1-J),J=1,NFLEVG)
  ! * For the non-zero imaginary parts, print the eigenvectors.
  DO JL2=1,NFLEVG
    IF (ABS(ZFI(NFLEVG+1-JL2)) > 0.0_JPRB) THEN
      WRITE(UNIT=KULOUT,FMT='(1X)')
      WRITE(KULOUT,FMT='(A,I4)') ' SIMO for VP number ',JL2
      WRITE(KULOUT,FMT='(14(1X,E8.2))')(SIMO(JL1,JL2),JL1=1,NFLEVG)
    ENDIF
  ENDDO
  WRITE(UNIT=KULOUT,FMT='(1X)')
  ISTOP=1
ENDIF

!*       1.8   SOLVE INVERSE HELMOLTZ EQUATION IF LSIDG=.T.

IF (LSIDG) THEN
  SITIME=-1.0_JPRB
  CALL SUHEG(YDGEOMETRY,YDRIP,YDDYN)
ENDIF
IF (LESIDG) THEN
  SITIME=-1.0_JPRB
  CALL SUEHEG(YDGEOMETRY,YDDYN,YDEDYN,YDRIP)
ENDIF

!*       1.9   CALCULATE INVERSE OF THE B MATRIX.

CALL MINV_CALLER(.TRUE.,NFLEVG,SIB,SIBI)

  !*       1.10  ISSUE THE VERDICT ON B MATRIX.

IF (ISTOP == 1) THEN
  CALL ABOR1("SUSI: Eigenvalues of B are not all real > 0")
ELSE
  WRITE(KULOUT,*) ' SUSI: B has all its eigenvalues real and >0; that is OK. ' 
ENDIF

!*       2.    CALCULATE EIGENVALUES FOR THE ANALYSIS OF STABILITY 
!              FOR THE LINEAR MODEL.
!              --------------------------------------------------

IF (LDYN_STABAN) THEN
  CALL SUDYN_STABILITY(YDCST, YDGEOMETRY,YDDYN,KULOUT)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSI',1,ZHOOK_HANDLE)
END SUBROUTINE SUSI
