#ifdef VPP
!OCL SCALAR
#endif
SUBROUTINE SUNHSI(YDCST, YDGEOMETRY,YDRIP,YDDYN,YDEDYN,KULOUT,LDEB)

!**** *SUNHSI*   - Initialize constants and control for the dynamics:
!                  Semi-implicit scheme operators.
!                  NHEE model.

!     Purpose.
!     --------
!        Computes the operator B=mu*nu+gamma*tau which is used in
!        the semi-implicit scheme
!        in order to work on the Hough space.
!        The inverse 1/B and also the eigenvalues and the eigenvectors
!        of B are computed.
!        * B is stored in array SIB.
!        * The eigenvalues of B are stored in array SIVP.
!        * The eigenvectors of B are stored in array SIMO.
!        * The inverse of SIMO is stored in array SIMI.
!        * The inverse of B is stored in array SIBI.
!        * The call to SUNHHEG fills the arrays SIHEG and SIHEG2.
!        In the NHEE model, all these quantities depend on Delta t.

!**   Interface.
!     ----------
!        *CALL* *SUNHSI(KULOUT)

!     Explicit arguments :
!     --------------------
!        KULOUT : Logical unit for the output
!        LDEB   : .T.: Computations of Part 1.1 are done.

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
!        Original : 08-Dec-2004

! Modifications
! -------------
!   Modified : 09-Oct-2007 by K. YESSAD: possibility to have a specific
!              value of (LVERTFE,NDLNPR) in the SI linear model.
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   I. Santos, I. Martinez (Mar 2010): Variable map factor for RTM.
!   K. Yessad (Dec 2011): Use GPHPRE.
!   M. Fisher   7-March-2012 Move SIBI out of Jb (to allow late Jb setup)
!   P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!   K. Yessad (July 2014): Move some variables.
!   K. Yessad (Sept 2014): Add printing of ZCOR2.
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Introduce NHQE model.
!   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!   J. Vivoda and P. Smolikova (Sep 2020): Add analysis of stability.
!   H Petithomme (Dec 2020): optimization on gphpre
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
USE YOMMP0   , ONLY : NPRINTLEV
USE YOMDYNA  , ONLY : YRDYNA
USE YOMDYN   , ONLY : TDYN
USE YEMDYN   , ONLY : TEDYN
USE YOMCST   , ONLY : TCST
USE INTDYN_MOD,ONLY : YYTXYB
USE YOMRIP   , ONLY : TRIP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)  :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(INOUT):: YDDYN
TYPE(TEDYN)       ,INTENT(INOUT):: YDEDYN
TYPE(TRIP)        ,INTENT(INOUT):: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)  :: KULOUT
LOGICAL           ,INTENT(IN)  :: LDEB

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JPLL
PARAMETER(JPLL=1)

REAL(KIND=JPRB) :: ZPRES(JPLL,0:YDGEOMETRY%YRDIMV%NFLEVG),ZPRESF(JPLL,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZXYB(JPLL,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)
REAL(KIND=JPRB) :: ZWO(YDGEOMETRY%YRDIMV%NFLEVG+1),ZA(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZFR(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZFI(YDGEOMETRY%YRDIMV%NFLEVG),ZMO(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZFN(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZID(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZCOR(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG),ZCOR2(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IWO(YDGEOMETRY%YRDIMV%NFLEVG+1)
INTEGER(KIND=JPIM) :: IER, ISIGN, ISTOP, JL1, JL2, JLEV, JLON, J

REAL(KIND=JPRB) :: ZSUM, ZTES, ZTER

LOGICAL :: LLNH_NOC1C2

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "eigsol.h"
#include "minv_caller.h"
#include "gphpre.intfb.h"
#include "scordo.intfb.h"
#include "si_cccor.intfb.h"
#include "sudyn_stability.intfb.h"
#include "suenhheg.intfb.h" 
#include "sunhbmat.intfb.h"
#include "sunhheg.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUNHSI',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
  & YDMP=>YDGEOMETRY%YRMP, YDVAB=>YDGEOMETRY%YRVAB)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, LSIDG=>YDDYN%LSIDG, LDYN_STABAN=>YDDYN%LDYN_STABAN, &
 & SIALPH=>YDDYN%SIALPH, SIB=>YDDYN%SIB, SIBI=>YDDYN%SIBI, &
 & SIDELP=>YDDYN%SIDELP, SIDPHI=>YDDYN%SIDPHI, SILNPR=>YDDYN%SILNPR, &
 & SIMI=>YDDYN%SIMI, SIMO=>YDDYN%SIMO, SIPR=>YDDYN%SIPR, SIRDEL=>YDDYN%SIRDEL, &
 & SITIME=>YDDYN%SITIME, SITLAF=>YDDYN%SITLAF, SITLAH=>YDDYN%SITLAH, &
 & SITR=>YDDYN%SITR, SITRAM=>YDDYN%SITRAM, SIVP=>YDDYN%SIVP, VNORM=>YDDYN%VNORM, &
 & LESIDG=>YDEDYN%LESIDG, &
 & TDT=>YDRIP%TDT)
!     ------------------------------------------------------------------

!*       1.    INITIALIZE SEMI-IMPLICIT AND VERTICAL NORMAL MODES
!              FOR THE 3D PRIMITIVE EQUATION MODEL.
!              --------------------------------------------------

!*       1.1   INITIALIZE SEMI-IMPLICIT COMPUTATION.

ISTOP = 0

IF (LDEB) THEN

! * This part does not depend on Delta t.

  DO JLON=1,JPLL
    ZPRES(JLON,NFLEVG)=SIPR
  ENDDO

  CALL GPHPRE(JPLL,NFLEVG,JPLL,JPLL,YDVAB,ZPRES,PXYB=ZXYB,PRESF=ZPRESF,LRTGR=.FALSE.,&
    LRPP=.FALSE.)

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

  SITIME=-1.0_JPRB

  ! * Case where the C1 or C2 constraints are not matched.
  LLNH_NOC1C2=.NOT.(YRDYNA%NDLNPR == 1)
  IF (LLNH_NOC1C2) THEN

    ! * Compute COR * I = -RT* A1
    !   Additionally, compute COR2 * I = (g**2/(N**2 C**2))(1/r) Lv* A2 = (Cv/rR) Lv* A2
    DO JLEV=1,NFLEVG
      DO JLON=1,NFLEVG
        ZID(JLON,JLEV)=1.0_JPRB-MIN(ABS(REAL(JLON-JLEV,JPRB)),1.0_JPRB)
      ENDDO
    ENDDO
    CALL SI_CCCOR(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,NFLEVG,ZID,ZCOR,POU2=ZCOR2)

    IF (NPRINTLEV > 0) THEN
      WRITE(KULOUT,'(A)') ''
      WRITE(KULOUT,'(A)') ' SUNHSI: Content of COR:'
      DO JLEV=1,NFLEVG
        DO JLON=1,NFLEVG
          WRITE(KULOUT,'(A,I3,A,I3,A,E10.4)')&
           & ' COR(',JLON,',',JLEV,') = ',ZCOR(JLON,JLEV)
        ENDDO
        WRITE(KULOUT,'(A)') ''
      ENDDO
      WRITE(KULOUT,'(A)') ''
      WRITE(KULOUT,'(A)') ' SUNHSI: Content of COR2:'
      DO JLEV=1,NFLEVG
        DO JLON=1,NFLEVG
          WRITE(KULOUT,'(A,I3,A,I3,A,E10.4)')&
           & ' COR2(',JLON,',',JLEV,') = ',ZCOR2(JLON,JLEV)
        ENDDO
        WRITE(KULOUT,'(A)') ''
      ENDDO
    ENDIF

    ! * Compute the eigenvalues of COR.
    CALL EIGSOL(NFLEVG,NFLEVG,ZCOR,ZFR,ZFI,1,ZMO,IWO,ZWO,IER)

    DO JLEV=1,NFLEVG
      ZFN(JLEV)=SQRT(ZFR(JLEV)*ZFR(JLEV)+ZFI(JLEV)*ZFI(JLEV))
    ENDDO

    WRITE(KULOUT,'(A)') ''
    WRITE(KULOUT,'(A,E10.4)') ' SUNHSI: Highest eigenvalue module of COR:', MAXVAL(ZFN)

    ! * Compute the eigenvalues of COR2-Id.
    DO JLEV=1,NFLEVG
      DO JLON=1,NFLEVG
        ZID(JLON,JLEV)=ZID(JLON,JLEV)*SITR/SITRAM(JLEV)
      ENDDO
    ENDDO
    CALL EIGSOL(NFLEVG,NFLEVG,ZCOR2-ZID,ZFR,ZFI,1,ZMO,IWO,ZWO,IER)

    DO JLEV=1,NFLEVG
      ZFN(JLEV)=SQRT(ZFR(JLEV)*ZFR(JLEV)+ZFI(JLEV)*ZFI(JLEV))
    ENDDO

    WRITE(KULOUT,'(A)') ''
    WRITE(KULOUT,'(A,E10.4)') ' SUNHSI: Highest eigenvalue module of COR2 - I:', MAXVAL(ZFN)

  ENDIF

ENDIF

IF(SITIME /= TDT)THEN

  SITIME=TDT

  !*       1.2   INITIALIZE MATRIX B (SIB)

  CALL SUNHBMAT(YDCST, YDGEOMETRY,YDDYN)

  !*       1.3   SOLVE EIGEN PROBLEM

  IER=0
  ZTES=0.0_JPRB

  DO JL1=1,NFLEVG
    DO JL2=1,NFLEVG
      ZA(JL1,JL2)=SIB(JL1,JL2)
    ENDDO
  ENDDO

  CALL EIGSOL(NFLEVG,NFLEVG,ZA,ZFR,ZFI,1,ZMO,IWO,ZWO,IER)
  IF(IER /= 0)THEN
    CALL ABOR1('SUNHSI: ABORT 1.3 IN SUNHSI AFTER CALL TO RG')
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
  ZTER=MINVAL(ZFR)

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

  DO JL1=1,NFLEVG
    DO JL2=1,NFLEVG
      ZMO(JL1,JL2) = SIMO(JL1,JL2)
    ENDDO
  ENDDO

  !*       1.6   INVERSE MODE MATRIX

  CALL MINV_CALLER(.TRUE.,NFLEVG,SIMO,SIMI)

  !*       1.7   PRINT OUT INFORMATION ON VERTICAL MODES

  WRITE(KULOUT,*) ''
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
        WRITE(KULOUT,'(''IMAGINARY PART OF EIGENVALUE :'',E20.14)') ZFI(NFLEVG+1-JL2)
        WRITE(UNIT=KULOUT,FMT='(1X)')
        WRITE(KULOUT,FMT='(A,I4)') ' SIMO for VP number ',JL2
        WRITE(KULOUT,FMT='(14(1X,E8.2))')(SIMO(JL1,JL2),JL1=1,NFLEVG)
      ENDIF
    ENDDO
    ISTOP=2
  ELSE
    IF (ZTER <= 0.0_JPRB) ISTOP=1
  ENDIF

  !*       1.8   SOLVE INVERSE HELMOLTZ EQUATION IF LSIDG=.T.

  IF (LSIDG) THEN
    CALL SUNHHEG(YDGEOMETRY,YDDYN)
  ENDIF

  IF (LESIDG) THEN
    CALL SUENHHEG(YDGEOMETRY,YDDYN,YDEDYN)
  ENDIF

  !*       1.9   CALCULATE INVERSE OF THE B MATRIX.

  CALL MINV_CALLER(.TRUE.,NFLEVG,SIB,SIBI)

ENDIF

  !*       1.10  ISSUE THE VERDICT ON B MATRIX.

WRITE(UNIT=KULOUT,FMT='(1X)')
IF (ISTOP == 1) THEN
  WRITE(KULOUT,*) ' WARNING SUNHSI: B has some real eigenvalues < or = 0! &
   & Linear system may be unstable! '
ELSEIF (ISTOP == 2) THEN
  WRITE(KULOUT,*) " WARNING SUNHSI: Eigenvalues of B are not all real > 0"
ENDIF


!*       2.    CALCULATE EIGENVALUES FOR THE ANALYSIS OF STABILITY 
!              FOR THE LINEAR MODEL.
!              --------------------------------------------------

IF(LDYN_STABAN)THEN
  CALL SUDYN_STABILITY(YDCST, YDGEOMETRY,YDDYN,KULOUT)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUNHSI',1,ZHOOK_HANDLE)
END SUBROUTINE SUNHSI
