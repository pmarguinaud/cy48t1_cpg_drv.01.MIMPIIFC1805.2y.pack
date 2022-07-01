SUBROUTINE GNHQE_CONV_TEMPE(&
 ! ----- INPUT ---------------------------------------------------------------
 & YDGEOMETRY,LD_TT_TO_T,KDIM_GFL,KST,KEND,PQCHA,PSP, &
 ! ----- INPUT-OUTPUT --------------------------------------------------------
 & PT, &
 ! ----- OPTIONAL INPUT ------------------------------------------------------
 & KGFLTYP,PGFL,KDDER,PQCHAL,PQCHAM, &
 ! ----- OPTIONAL INPUT-OUTPUT -----------------------------------------------
 & PTL,PTM)

!* GNHQE_CONV_TEMPE - Conversion between Tt and T for NHQE model.

! Purpose
! -------
!   Computes Tt from T, or T from Tt.

! Interface
! ---------
!  * INPUT:
!     YDGEOMETRY      : structure containing geometry
!     LD_TT_TO_T      : T: conversion Tt -> T; F: conversion T -> Tt.
!     KDIM_GFL        : last dimension of PGFL
!     KST             : start of work.
!     KEND            : end of work.
!     PQCHA           : log(pre/prehyd)
!     PSP             : log(prehyds)

!  * INPUT-OUTPUT:
!     PT              : T or Tt

!  * OPTIONAL INPUT:
!     KGFLTYP         : cf. KGFLTYP in GPRCP.
!     PGFL            : unified_treatment grid-point fields at t
!     KDDER           : KDDER present and >0: conversion applied to grad(T) too,
!                       with some approximations depending on KDDER.
!     PQCHAL          : zonal comp of grad(log(pre/prehyd))
!     PQCHAM          : merid comp of grad(log(pre/prehyd))

!  * OPTIONAL INPUT-OUTPUT:
!     PTL             : zonal comp of grad(temperature variable)
!     PTM             : merid comp of grad(temperature variable)

! Externals
! ---------

! Method
! ------
!   See documentation

! Reference
! ---------

! Author
! ------
!   K. Yessad (METEO-FRANCE/CNRM/GMAP)
!   Date: June 2017

! Modifications
! -------------
!   H Petithomme (Dec 2020): delete useless ZXYB
!------------------------------------------------------------------------------

USE YOM_YGFL     , ONLY : TYPE_GFLD
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE INTDYN_MOD   , ONLY : YYTXYB
USE YOMCST       , ONLY : RD, RCPD

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY) ,INTENT(IN)             :: YDGEOMETRY
LOGICAL        ,INTENT(IN)             :: LD_TT_TO_T
INTEGER(KIND=JPIM),INTENT(IN)          :: KDIM_GFL
INTEGER(KIND=JPIM),INTENT(IN)          :: KST
INTEGER(KIND=JPIM),INTENT(IN)          :: KEND
REAL(KIND=JPRB),INTENT(IN)             :: PQCHA(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(IN)             :: PSP(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB),INTENT(INOUT)          :: PT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KGFLTYP
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,KDIM_GFL)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KDDER
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PQCHAL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PQCHAM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PTL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),OPTIONAL,INTENT(INOUT) :: PTM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!------------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV
LOGICAL :: LLDER1,LLDER2,LLMOIST_RCP

REAL(KIND=JPRB) :: ZPREH  (YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPREF  (YDGEOMETRY%YRDIM%NPROMA,  YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZKAP   (YDGEOMETRY%YRDIM%NPROMA,  YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZEQCHAF(YDGEOMETRY%YRDIM%NPROMA,  YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------

#include "gnhpre.intfb.h"
#include "gphpre.intfb.h"
#include "gprcp_pgfl.intfb.h"

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHQE_CONV_TEMPE',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA,NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,YDVAB=>YDGEOMETRY%YRVAB)

!------------------------------------------------------------------------------

!* 1. GRID POINT COMPUTATIONS
!     -----------------------

IF (PRESENT(KDDER)) THEN
  LLDER1=(PRESENT(KDDER).AND.(KDDER == 1).AND.PRESENT(PTL).AND.PRESENT(PTM))
  LLDER2=(PRESENT(KDDER).AND.(KDDER == 2).AND.PRESENT(PQCHAL).AND.PRESENT(PQCHAM).AND.PRESENT(PTL).AND.PRESENT(PTM))
ELSE
  LLDER1=.FALSE.
  LLDER2=.FALSE.
ENDIF
LLMOIST_RCP=(PRESENT(KGFLTYP).AND.PRESENT(PGFL))

! * Compute pressures:
ZPREH(KST:KEND,NFLEVG)=EXP(PSP(KST:KEND))
CALL GPHPRE(NPROMA,NFLEVG,KST,KEND,YDVAB,ZPREH,PRESF=ZPREF)

! * Compute Kap=R/Cp:
IF (LLMOIST_RCP) THEN
  CALL GPRCP_PGFL(NPROMA,KST,KEND,NFLEVG,PGFL=PGFL,PKAP=ZKAP,KGFLTYP=KGFLTYP)
ELSE
  ZKAP(KST:KEND,1:NFLEVG)=RD/RCPD
ENDIF

! * Compute exp(Kap*Qcha):
CALL GNHPRE(YDGEOMETRY,NPROMA,NFLEVG,KST,KEND,PQCHA,ZPREF,PKAP=ZKAP,PEQCHAF=ZEQCHAF)

IF ( LD_TT_TO_T ) THEN

  ! * Transform T: compute T from Tt:
  DO JLEV=1,NFLEVG
    PT(KST:KEND,JLEV)=PT(KST:KEND,JLEV)*ZEQCHAF(KST:KEND,JLEV)
  ENDDO

  ! * Transform grad(T): compute grad(T) from grad(Tt) with some approximations
  IF (LLDER1) THEN
    DO JLEV=1,NFLEVG
      PTL(KST:KEND,JLEV)=PTL(KST:KEND,JLEV)*ZEQCHAF(KST:KEND,JLEV) 
      PTM(KST:KEND,JLEV)=PTM(KST:KEND,JLEV)*ZEQCHAF(KST:KEND,JLEV) 
    ENDDO
  ELSEIF (LLDER2) THEN
    DO JLEV=1,NFLEVG
      PTL(KST:KEND,JLEV)=PTL(KST:KEND,JLEV)*ZEQCHAF(KST:KEND,JLEV) &
       & +ZKAP(KST:KEND,JLEV)*PT(KST:KEND,JLEV)*PQCHAL(KST:KEND,JLEV)
      PTM(KST:KEND,JLEV)=PTM(KST:KEND,JLEV)*ZEQCHAF(KST:KEND,JLEV) &
       & +ZKAP(KST:KEND,JLEV)*PT(KST:KEND,JLEV)*PQCHAM(KST:KEND,JLEV)
    ENDDO
  ENDIF

ELSE

  ! * Transform T: compute Tt from T:
  DO JLEV=1,NFLEVG
    PT(KST:KEND,JLEV)=PT(KST:KEND,JLEV)/ZEQCHAF(KST:KEND,JLEV)
  ENDDO

  ! * Transform grad(T): compute grad(Tt) from grad(T) with some approximations
  IF (LLDER1) THEN
    DO JLEV=1,NFLEVG
      PTL(KST:KEND,JLEV)=PTL(KST:KEND,JLEV)/ZEQCHAF(KST:KEND,JLEV)
      PTM(KST:KEND,JLEV)=PTM(KST:KEND,JLEV)/ZEQCHAF(KST:KEND,JLEV)
    ENDDO
  ELSEIF (LLDER2) THEN
    DO JLEV=1,NFLEVG
      PTL(KST:KEND,JLEV)=PTL(KST:KEND,JLEV)/ZEQCHAF(KST:KEND,JLEV) &
       & -ZKAP(KST:KEND,JLEV)*PT(KST:KEND,JLEV)*PQCHAL(KST:KEND,JLEV)
      PTM(KST:KEND,JLEV)=PTM(KST:KEND,JLEV)/ZEQCHAF(KST:KEND,JLEV) &
       & -ZKAP(KST:KEND,JLEV)*PT(KST:KEND,JLEV)*PQCHAM(KST:KEND,JLEV)
    ENDDO
  ENDIF

ENDIF  ! LD_TT_TO_T

!------------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHQE_CONV_TEMPE',1,ZHOOK_HANDLE)

END SUBROUTINE GNHQE_CONV_TEMPE
