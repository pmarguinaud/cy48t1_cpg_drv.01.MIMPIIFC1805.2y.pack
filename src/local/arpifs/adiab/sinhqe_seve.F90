SUBROUTINE SINHQE_SEVE(YDCST, YDGEOMETRY,YDDYN,KLEV,KLON,PV1,PV2,KNLON,LD_LSTAR)

!**** *SINHQE_SEVE* - SEcond VErtical derivative operator in Semi - Implicit for NHQE model

!     Purpose.
!     --------
!           Operator L_star or L_star_star to compute vertical Laplacian in NHQE non-hydrostatic
!           dynamics: formulation for semi - implicit

!           The discretisation of the vertical Laplacian must remain
!           consistent with the discretisation of the equivalent quantity
!           in the non linear part (pieces of codes spread in GNHQE_PREH,
!           GNHQE_TNDLAGADIAB_GW and GNHQE_TNDLAGADIAB_SVD).

!           This discretisation may differ from the NHEE one (in particuliar at levels 1, 2 and NFLEVG),
!           this is why the code is in a separate routine.

!           This vertical Laplacian is applied to a quantity "X" (in PV1)
!           which matches the following conditions:
!           X_top = 0
!           X_surf = X(l=nflevg)
!           (see for example the calculation of the top and bottom
!           pressure departure variable in GNHQE_PREH).

!**   Interface.
!     ----------
!        *CALL* *SINHQE_SEVE(...)

!        Explicit arguments :
!        --------------------
!        KLEV   : DISTANCE IN MEMORY BETWEEN VALUES OF THE VARIABLES      (IN)
!                 PV1 OR PV2 AT THE VERTICAL
!        KLON   : DISTANCE IN MEMORY BETWEEN VALUES OF THE VARIABLES      (IN)
!                 PV1 OR PV2 AT THE SAME LEVEL
!        KNLON  : NUMBER OF VERTICAL COLUMNS TREATED                      (IN)

!           TYPICAL VALUES ARE  (KLEV,KLON)=(NDLSUR,1)  FOR GRID POINT ARRAY
!                               (KLEV,KLON)=(1,NFLSUR)  FOR SPECTRAL ARRAY

!        PV1   : INPUT VARIABLE                                           (IN)
!        PV2   : DERIVED VARIABLE BY VERTICAL LAPLACIAN                   (OUT)
!        LD_LSTAR : T: operator L_star                                    (OPT IN)
!                   F: operator L_star_star

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Arpege/Aladin documentation

!     Author.
!     -------
!      Radmila Bubnova GMAP/COMPAS, stage MRE
!      Original : 93-03-21

!     Modifications.
!     --------------
!      P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!      K. Yessad (Dec 2016): Prune obsolete options.
!      K. Yessad (June 2017): Vertical-dependent SITRA; introduce NHQE model.
!      J. Vivoda and P. Smolikova (Sep 2020): VFE pruning.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : TCST
USE YOMDYN       , ONLY : TDYN
USE YOMCVER      , ONLY : LVFE_LAPL_BC

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNLON 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV1(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PV2(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
LOGICAL,OPTIONAL  ,INTENT(IN)    :: LD_LSTAR

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZA(YDGEOMETRY%YRDIMV%NFLEVG), ZC(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPRTOP

INTEGER(KIND=JPIM) :: IBOT1, IBOT2, ILEV1, ILEV2, ILEV3, ITOP1, ITOP2, ITOP3
INTEGER(KIND=JPIM) :: JLEV, JLON, IDT

REAL(KIND=JPRB) :: ZRATIO_SITR(YDGEOMETRY%YRDIMV%NFLEVG)
LOGICAL         :: LL_LSTAR

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SINHQE_SEVE',0,ZHOOK_HANDLE)

ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & SILNPR=>YDDYN%SILNPR,SIALPH=>YDDYN%SIALPH, &
 & SIPR=>YDDYN%SIPR, SITLAF=>YDDYN%SITLAF, &
 & SITLAH=>YDDYN%SITLAH,SITR=>YDDYN%SITR,SITRAM=>YDDYN%SITRAM)
!     ------------------------------------------------------------------

IF (LVFE_LAPL_BC) CALL ABOR1(' SINHQE_SEVE: VFE Laplacian not yet implemented')


!*       1.    BOUNDARY CONDITIONS FOR FD
!              --------------------------

! --- Finite differences for top and bottom b.c. of L_star ---

! ** Top (full level number 1):

ZPRTOP=YDVAB%VAH(0)+YDVAB%VBH(0)*SIPR
ZA(1)=ZPRTOP/(SILNPR(1)*(SITLAF(1)-ZPRTOP))
ZC(1)=SITLAF(2)/(SITLAF(2) - SITLAF(1)*SILNPR(1)*(SIALPH(2)-SILNPR(2)))
DO JLON = 1, KNLON
  ITOP1 = 1 + (JLON-1)*KLON
  ITOP2 = ITOP1 + KLEV
  PV2(ITOP1) = -ZA(1)*PV1(ITOP1)+ZC(1)*(PV1(ITOP2)-PV1(ITOP1))
ENDDO

! ** Full level number 2:
ZA(2)=SITLAF(1)/(SILNPR(2)*(SITLAF(2)/SILNPR(1)-SITLAF(1)*(SIALPH(2)-SILNPR(2))))
ZC(2)=SITLAF(3)/(SILNPR(2)*(SITLAF(3)-SITLAF(2)))
DO JLON = 1, KNLON
  ITOP1 = 1 + (JLON-1)*KLON
  ITOP2 = ITOP1 + KLEV
  ITOP3 = ITOP2 + KLEV
  PV2(ITOP2) = ZA(2)*(PV1(ITOP1)-PV1(ITOP2))+ZC(2)*(PV1(ITOP3)-PV1(ITOP2))
ENDDO

! ** Bottom (full level number NFLEVG):
ZA(NFLEVG)=SITLAF(NFLEVG-1)/(SILNPR(NFLEVG)*(SITLAF(NFLEVG)-SITLAF(NFLEVG-1)))
ZC(NFLEVG)=YDCST%RKAPPA*SITLAH(NFLEVG)/(SILNPR(NFLEVG)*(YDCST%RKAPPA*SITLAH(NFLEVG)+(1.0_JPRB-YDCST%RKAPPA)*SITLAF(NFLEVG)))
DO JLON = 1, KNLON
  IBOT1 = 1 + (NFLEVG-2)*KLEV + (JLON-1)*KLON
  IBOT2 = IBOT1 + KLEV
  PV2(IBOT2) = ZA(NFLEVG)*(PV1(IBOT1)-PV1(IBOT2)) - ZC(NFLEVG)*PV1(IBOT2)
ENDDO

!     ------------------------------------------------------------------

!*       2.    INNER LAYERS (AND BOUNDARY CONDITIONS FOR VFE)
!              ----------------------------------------------

! * Compute the coefficients A and C at full levels 2 to nflevg-1:
DO JLEV = 3, NFLEVG-1
  ZA(JLEV) = SITLAF(JLEV-1)/(SILNPR(JLEV)*(SITLAF(JLEV)-SITLAF(JLEV-1)))
  ZC(JLEV) = SITLAF(JLEV+1)/(SILNPR(JLEV)*(SITLAF(JLEV+1)-SITLAF(JLEV)))
ENDDO

! * Other Levels:
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JLON,ILEV1,ILEV2,ILEV3)
DO JLEV = 3, NFLEVG-1
  DO JLON = 1, KNLON
    ILEV1 = 1 + (JLEV-2)*KLEV + (JLON-1)*KLON
    ILEV2 = ILEV1 + KLEV
    ILEV3 = ILEV2 + KLEV
    PV2(ILEV2) = ZA(JLEV)*(PV1(ILEV1)-PV1(ILEV2))+ZC(JLEV)*(PV1(ILEV3)-PV1(ILEV2))
  ENDDO
ENDDO
!$OMP END PARALLEL DO

!     ------------------------------------------------------------------

!*       3.    MULTIPLY BY SITR/SITRAM
!              ----------------------

IF(PRESENT(LD_LSTAR)) THEN
  LL_LSTAR=LD_LSTAR
ELSE
  LL_LSTAR=.FALSE.
ENDIF

IF(.NOT.LL_LSTAR) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JLON,IDT)
  DO JLEV=1,NFLEVG
    ZRATIO_SITR(JLEV)=SITR/SITRAM(JLEV)
    DO JLON=1,KNLON
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      PV2(IDT)=ZRATIO_SITR(JLEV)*PV2(IDT)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SINHQE_SEVE',1,ZHOOK_HANDLE)
END SUBROUTINE SINHQE_SEVE
