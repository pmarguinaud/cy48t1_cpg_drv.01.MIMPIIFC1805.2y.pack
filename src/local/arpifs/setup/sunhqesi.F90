#ifdef VPP
!OCL SCALAR
#endif
SUBROUTINE SUNHQESI(YDCST, YDGEOMETRY,YDRIP,YDDYN,KULOUT,LDEB)

!**** *SUNHQESI*   - Initialize constants and control for the dynamics:
!                    Semi-implicit scheme operators.
!                    NHQE model.

!     Purpose.
!     --------
!        Computes the operator B for NHQE model, which is used in
!        the semi-implicit scheme in order to work on the Hough space.
!        The inverse 1/B and also the eigenvalues and the eigenvectors
!        of B are computed.
!        * B is stored in array SIB.
!        * The eigenvalues of B are stored in array SIVP.
!        * The eigenvectors of B are stored in array SIMO.
!        * The inverse of SIMO is stored in array SIMI.
!        * The inverse of B is stored in array SIBI.
!        * The call to SUHEG fills the arrays SIHEG and SIHEG2.
!        In the NHQE model, all these quantities depend on Delta t.

!**   Interface.
!     ----------
!        *CALL* *SUNHQESI(...)

!     Explicit arguments :
!     --------------------
!        YDGEOMETRY   : structure for geometry                   (IN)
!        YDRIP        : structure for timestep                   (INOUT)
!        YDDYN        : structure for dynamics                   (INOUT)
!        KULOUT       : Logical unit for the output              (IN)
!        LDEB         : .T.: Computations of Part 1.1 are done.  (IN)

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
!        Documentation (IDSI)

!     Author.
!     -------
!        K. YESSAD (after SUSI and SUNHSI).
!        Original : march 2017

! Modifications
! -------------
!   H Petithomme (Dec 2020): optimisation on gphpre
!     ------------------------------------------------------------------

USE YOMRIP       , ONLY : TRIP
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
USE YOMDYNA  , ONLY : LNHQE_C2
USE YOMDYN   , ONLY : TDYN
USE YOMCST   , ONLY : TCST
USE YOMCVER  , ONLY : LVERTFE
USE INTDYN_MOD,ONLY : YYTXYB

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)     :: YDGEOMETRY
TYPE(TRIP)        ,INTENT(INOUT)  :: YDRIP
TYPE(TDYN)        ,INTENT(INOUT)  :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)     :: KULOUT
LOGICAL           ,INTENT(IN)     :: LDEB

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JPLL
PARAMETER(JPLL=1)

REAL(KIND=JPRB) :: ZPRES(JPLL,0:YDGEOMETRY%YRDIMV%NFLEVG),ZPRESF(JPLL,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZXYB(JPLL,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)
REAL(KIND=JPRB) :: ZWO(YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZA(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZFR(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZFI(YDGEOMETRY%YRDIMV%NFLEVG),ZMO(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZIDI(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZIDO(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSKAP(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZLAP(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZLAPKAP(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZLAPKAP_GE(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZLAPKAP_C2(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IWO(YDGEOMETRY%YRDIMV%NFLEVG+1)
INTEGER(KIND=JPIM) :: IER, ISIGN, JL1, JL2, JLEV, JLON, J

REAL(KIND=JPRB) :: ZDET, ZSUM, ZTES, ZTER

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "eigsol.h"
#include "minv_caller.h"
#include "gphpre.intfb.h"
#include "scordo.intfb.h"
#include "sunhqebmat.intfb.h"
#include "suheg.intfb.h"
#include "siskap.intfb.h"
#include "siskapi.intfb.h"
#include "silkap.intfb.h"
#include "silkapi.intfb.h"
#include "sinhqe_seve.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUNHQESI',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIMV=>YDGEOMETRY%YRDIMV,YDVAB=>YDGEOMETRY%YRVAB)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & LSIDG=>YDDYN%LSIDG, SIALPH=>YDDYN%SIALPH, SIB=>YDDYN%SIB, SIBI=>YDDYN%SIBI, &
 & SIDELP=>YDDYN%SIDELP, SIDPHI=>YDDYN%SIDPHI, SILNPR=>YDDYN%SILNPR, &
 & SIMI=>YDDYN%SIMI, SIMO=>YDDYN%SIMO, SIPR=>YDDYN%SIPR, SIRDEL=>YDDYN%SIRDEL, &
 & SITIME=>YDDYN%SITIME, SITLAF=>YDDYN%SITLAF, SITLAH=>YDDYN%SITLAH, &
 & SITR=>YDDYN%SITR, SIVP=>YDDYN%SIVP, VNORM=>YDDYN%VNORM, &
 & TDT=>YDRIP%TDT)

!     ------------------------------------------------------------------

!*       1.    INITIALIZE SEMI-IMPLICIT AND VERTICAL NORMAL MODES
!              FOR THE 3D PRIMITIVE EQUATION MODEL.
!              --------------------------------------------------

!*       1.1   INITIALIZE SEMI-IMPLICIT COMPUTATION.

IF (LDEB) THEN

! * This part does not depend on Delta t.

  DO JLON=1,JPLL
    ZPRES(JLON,NFLEVG)=SIPR
  ENDDO

  CALL GPHPRE(JPLL,NFLEVG,JPLL,JPLL,YDVAB,ZPRES,PXYB=ZXYB,PRESF=ZPRESF,LRTGR=.FALSE.,&
   & LRPP=.FALSE.)

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

ENDIF

IF(SITIME /= TDT)THEN

  SITIME=TDT
  ZDET=0.0_JPRB

  !*       1.2   INITIALIZE MATRIX B (SIB)

  CALL SUNHQEBMAT(YDCST, YDGEOMETRY,YDDYN)

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
    CALL ABOR1('SUNHQESI: ABORT 1.3 IN SUNHQESI AFTER CALL TO RG')
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

  !*       1.6   INVERSE MODE MATRIX

  CALL MINV_CALLER(.TRUE.,NFLEVG,SIMO,SIMI)

  !*       1.7   PRINT OUT INFORMATION ON VERTICAL MODES

  ! * print SIB itself and additional checkings for "small" NFLEVG.

  IF (NFLEVG <= 15) THEN

    ! print SIB:
    WRITE(KULOUT,*) ''
    WRITE(KULOUT,*) ' SUNHQESI: coefficients of B '
    DO JL2=1,NFLEVG
      WRITE(KULOUT,'(15(1X,E9.3))') (SIB(J,JL2),J=1,NFLEVG)
    ENDDO
    WRITE(KULOUT,*) ''

    ! ZIDI contains the identity matrix "I":
    DO JL1=1,NFLEVG
      DO JL2=1,NFLEVG
        ZIDI(JL1,JL2)=1.0_JPRB-MIN(ABS(REAL(JL2-JL1,JPRB)),1.0_JPRB)
      ENDDO
    ENDDO

    ! check Sstar_kap**-1 Sstar_kap = Id
    IF (.NOT.LVERTFE) THEN
      CALL SISKAP(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,NFLEVG,ZIDI,ZSKAP)
      CALL SISKAPI(YDGEOMETRY,YDDYN,1,NFLEVG,NFLEVG,ZSKAP,ZIDO)
      WRITE(KULOUT,*) ''
      WRITE(KULOUT,*) ' SUNHQESI: coefficients of Sstar_kap**-1 Sstar_kap '
      DO JL2=1,NFLEVG
        WRITE(KULOUT,'(15(1X,E9.3))') (ZIDO(J,JL2),J=1,NFLEVG)
      ENDDO
      WRITE(KULOUT,*) ''
    ENDIF

    ! check - [ (Cpd/(Rd*Rd*Tstar)) (gamma tau + Rd Tstar nu) + QQstar ] LLstar = Id
    IF (LNHQE_C2) THEN
      WRITE(KULOUT,*) ''
      CALL SINHQE_SEVE(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZIDI,ZLAP,NFLEVG,LD_LSTAR=.TRUE.)
      CALL SILKAPI(YDCST, YDGEOMETRY,YDDYN,0,1,NFLEVG,NFLEVG,ZLAP,ZIDO)
      WRITE(KULOUT,*) ' SUNHQESI: coefficients of LLstar**-1 LLstar '
      WRITE(KULOUT,*) ' Check that they are close to Id (otherwise constraint C2 is not matched) '
      DO JL2=1,NFLEVG
        WRITE(KULOUT,'(15(1X,E9.3))') (ZIDO(J,JL2),J=1,NFLEVG)
      ENDDO
      WRITE(KULOUT,*) ''
    ENDIF

    ! check (Tstar/Tastar) [ LLstar + (Rd/Cpd - 1) d_star - Sstar_kap LLstar ] = 0 or very close to 0
    IF (LNHQE_C2) THEN
      CALL SILKAP(YDCST, YDGEOMETRY,YDDYN,.TRUE.,1,NFLEVG,NFLEVG,ZIDI,ZLAPKAP_C2)
      CALL SILKAP(YDCST, YDGEOMETRY,YDDYN,.FALSE.,1,NFLEVG,NFLEVG,ZIDI,ZLAPKAP_GE)
      WRITE(KULOUT,*) ''
      WRITE(KULOUT,*) ' SUNHQESI: coefficients of (Tstar/Tastar) [ LLstar + (Rd/Cpd - 1) d_star - Sstar_kap LLstar ] '
      WRITE(KULOUT,*) ' Check that they are close to 0 (otherwise constraint C2 is not matched) '
      DO JL2=1,NFLEVG
        WRITE(KULOUT,'(15(1X,E9.3))') (ZLAPKAP_GE(J,JL2)-ZLAPKAP_C2(J,JL2),J=1,NFLEVG)
      ENDDO
      WRITE(KULOUT,*) ''
    ENDIF

    ! check LLsstar_kap**-1 LLsstar_kap = Id
    CALL SILKAP(YDCST, YDGEOMETRY,YDDYN,LNHQE_C2,1,NFLEVG,NFLEVG,ZIDI,ZLAPKAP)
    CALL SILKAPI(YDCST, YDGEOMETRY,YDDYN,1,1,NFLEVG,NFLEVG,ZLAPKAP,ZIDO)
    WRITE(KULOUT,*) ''
    WRITE(KULOUT,*) ' SUNHQESI: coefficients of LLsstar_kap**-1 LLsstar_kap '
    DO JL2=1,NFLEVG
      WRITE(KULOUT,'(15(1X,E9.3))') (ZIDO(J,JL2),J=1,NFLEVG)
    ENDDO
    WRITE(KULOUT,*) ''

    ! check [LLsstar_kap**-1 Sstar_kap] [Sstar_kap**-1 LLsstar_kap] = Id
    IF (.NOT.LVERTFE) THEN
      CALL SILKAP(YDCST, YDGEOMETRY,YDDYN,LNHQE_C2,1,NFLEVG,NFLEVG,ZIDI,ZLAPKAP)
      CALL SISKAPI(YDGEOMETRY,YDDYN,1,NFLEVG,NFLEVG,ZLAPKAP,ZSKAP)
      CALL SILKAPI(YDCST, YDGEOMETRY,YDDYN,2,1,NFLEVG,NFLEVG,ZSKAP,ZIDO)
      WRITE(KULOUT,*) ''
      WRITE(KULOUT,*) ' SUNHQESI: coefficients of [LLsstar_kap**-1 Sstar_kap] [Sstar_kap**-1 LLsstar_kap] '
      DO JL2=1,NFLEVG
        WRITE(KULOUT,'(15(1X,E9.3))') (ZIDO(J,JL2),J=1,NFLEVG)
      ENDDO
      WRITE(KULOUT,*) ''
    ENDIF

  ENDIF

  ! * print information on vertical modes

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
    WRITE(UNIT=KULOUT,FMT='(1X)')
    WRITE(KULOUT,*) " WARNING: SUNHQESI: Eigenvalues of B are not all real > 0"
  ELSE
    IF (ZTER <= 0.0_JPRB) THEN
      WRITE(KULOUT,*) ' SUNHQESI: B has some real eigenvalues < or = 0! Linear system may be unstable! '
    ELSE
      WRITE(KULOUT,*) ' SUNHQESI: B has all its eigenvalues real and >0; that is OK. '
    ENDIF
  ENDIF

  !*       1.8   SOLVE INVERSE HELMOLTZ EQUATION IF LSIDG=.T.

  IF (LSIDG) THEN
    CALL SUHEG(YDGEOMETRY,YDRIP,YDDYN)
  ENDIF

  !*       1.9   CALCULATE INVERSE OF THE B MATRIX.

  CALL MINV_CALLER(.TRUE.,NFLEVG,SIB,SIBI)

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUNHQESI',1,ZHOOK_HANDLE)
END SUBROUTINE SUNHQESI
