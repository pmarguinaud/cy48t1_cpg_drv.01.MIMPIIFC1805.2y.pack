SUBROUTINE GNHQE_CONV_NHVAR(&
 ! ----- INPUT ---------------------------------------------------------------
 & YDGEOMETRY,YDGFL,LDMODEL_TO_FILE,&
 ! ----- INPUT-OUTPUT --------------------------------------------------------
 & PU,PV,PDIV,PT,PTL,PTM,PSPD,PSPDL,PSPDM,PSVD,PSP,PSPL,PSPM,&
 & PQ,PQL,PQM,PL,PLL,PLM,PI,PIL,PIM,PR,PS,PG,&
 ! ----- OUTPUT --------------------------------------------------------------
 & PSNHX&
 &)

!* GNHQE_CONV_NHVAR - Conversion of NH variables (model to file and vice-versa)
!                     for NHQE model.
!                     Grid-point part.

! Purpose
! -------
!   Converts prognostic variables between model and file for NHQE model.
!   Files contain T, "pre-prehyd" and "-G.dw". 
!   Model prognostic variables are Tt, log(pre/prehyd), d4=dver+NHX
!   Modifies spectral buffers:
!     SPT
!     SPSPD
!     SPSVD
!     SPNHX (only for d4 + file ---> model conversion)

! Interface
! ---------
!  * INPUT:
!     YDGEOMETRY      : structure containing geometry
!     YDGFL           : structure containing GFL
!     LDMODEL_TO_FILE : switch for model to file / file to model conversion

!  * INPUT-OUTPUT:
!     PU              : U-wind
!     PV              : V-wind
!     PDIV            : horizontal divergence
!     PT              : temperature equation variable
!     PTL             : zonal comp of grad(temperature equation variable)
!     PTM             : merid comp of grad(temperature equation variable)
!     PSPD            : NH pressure departure variable "spd"
!     PSPDL           : zonal comp of grad(spd)
!     PSPDM           : merid comp of grad(spd)
!     PSVD            : NH vertical divergence variable "svd"
!     PSP             : log(prehyds)
!     PSPL            : zonal comp of grad(log(prehyds))
!     PSPM            : merid comp of grad(log(prehyds))
!     PQ              : moisture
!     PQL             : zonal comp of grad(moisture)
!     PQM             : merid comp of grad(moisture)
!     PL              : liquid water
!     PLL             : zonal comp of grad(liquid water)
!     PLM             : merid comp of grad(liquid water)
!     PI              : ice
!     PIL             : zonal comp of grad(ice)
!     PIM             : merid comp of grad(ice)
!     PR              : rain
!     PS              : snow
!     PG              : graupels

!  * OUTPUT:
!     PSNHX           : "NHX" divergence variable for case nvdvar=4.

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
!   Date: March 2017

! Modifications
! -------------
!   H Petithomme (Dec 2020): optimisation on gphpre and gnhqe_xxd
!------------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMDYNA      , ONLY : NVDVAR, NPDVAR, L_RDRY_VD
USE YOM_YGFL     , ONLY : TYPE_GFLD
USE INTDYN_MOD   , ONLY : YYTXYB
use yomcver,only: lvertfe

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY) ,INTENT(IN)    :: YDGEOMETRY
TYPE(TYPE_GFLD),INTENT(IN)    :: YDGFL
LOGICAL        ,INTENT(IN)    :: LDMODEL_TO_FILE 
REAL(KIND=JPRB),INTENT(INOUT) :: PU   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PV   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PDIV (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PT   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PTL  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PTM  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PSPD (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PSPDL(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PSPDM(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PSVD (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PSP  (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB),INTENT(INOUT) :: PSPL (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB),INTENT(INOUT) :: PSPM (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB),INTENT(INOUT) :: PQ   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PQL  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PQM  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PL   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PLL  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PLM  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PI   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PIL  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PIM  (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PR   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PS   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(INOUT) :: PG   (YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB),INTENT(OUT)   :: PSNHX(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)

!------------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JIST,JLEV
INTEGER(KIND=JPIM) :: IST,IEND,IPROMA,IBL
LOGICAL :: LLNHQE_BALANCE

REAL(KIND=JPRB) :: ZOROGL (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZOROGM (YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZPDEP  (YDGEOMETRY%YRDIM%NPROMA,  YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZXYB   (YDGEOMETRY%YRDIM%NPROMA,  YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB%NDIM)
REAL(KIND=JPRB) :: ZPREH  (YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPREF  (YDGEOMETRY%YRDIM%NPROMA,  YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZNHX   (YDGEOMETRY%YRDIM%NPROMA,  YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZR     (YDGEOMETRY%YRDIM%NPROMA,  YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZKAP   (YDGEOMETRY%YRDIM%NPROMA,  YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZRT    (YDGEOMETRY%YRDIM%NPROMA,  YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZEQCHAF(YDGEOMETRY%YRDIM%NPROMA,  YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------

#include "abor1.intfb.h"
#include "gnhqe_nhx.intfb.h"
#include "gnhpre.intfb.h"
#include "gphpre.intfb.h"
#include "gprcp_qlirsg.intfb.h"
#include "gpskap.intfb.h"

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHQE_CONV_NHVAR',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA,NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & NGPTOT=>YDGEOMETRY%YRGEM%NGPTOT, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM,YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, &
 & YDOROG=>YDGEOMETRY%YROROG,YDVAB=>YDGEOMETRY%YRVAB)

!------------------------------------------------------------------------------

IF (.NOT.(NPDVAR == 2 .AND. NVDVAR == 4)) THEN
  CALL ABOR1(' GNHQE_CONV_NHVAR: wrong values for NPDVAR and NVDVAR')
ENDIF
IF (L_RDRY_VD) THEN
  ! Definition of "dver" in NHQE model uses R_moist.
  CALL ABOR1(' GNHQE_CONV_NHVAR: L_RDRY_VD should be F')
ENDIF

!* 1. GRID POINT COMPUTATIONS
!     -----------------------

!* 1.0 MULTIPLICATION BY MAP FACTOR

! orography
DO JIST=1,NGPTOT,NPROMA
  IST    = JIST                      ! first element
  IEND   = MIN(IST+NPROMA-1,NGPTOT)  ! last element
  IPROMA = IEND-IST+1                ! actual number of elements
  IBL=(JIST-1)/NPROMA+1
  ZOROGL(IST:IEND)=YDOROG(IBL)%OROGL(1:IPROMA)*YDGSGEOM(IBL)%GM(1:IPROMA)
  ZOROGM(IST:IEND)=YDOROG(IBL)%OROGM(1:IPROMA)*YDGSGEOM(IBL)%GM(1:IPROMA)
ENDDO

! horizontal wind:
DO JLEV=1,NFLEVG
  PU(1:NGPTOT,JLEV)=PU(1:NGPTOT,JLEV)*YDGSGEOM_NB%GM(1:NGPTOT)
  PV(1:NGPTOT,JLEV)=PV(1:NGPTOT,JLEV)*YDGSGEOM_NB%GM(1:NGPTOT)
ENDDO

! divergence:
DO JLEV=1,NFLEVG
  PDIV(1:NGPTOT,JLEV)=PDIV(1:NGPTOT,JLEV)*YDGSGEOM_NB%GM(1:NGPTOT)*YDGSGEOM_NB%GM(1:NGPTOT)
ENDDO

! surface pressure (GMVS)
PSPL(1:NGPTOT)=PSPL(1:NGPTOT)*YDGSGEOM_NB%GM(1:NGPTOT)
PSPM(1:NGPTOT)=PSPM(1:NGPTOT)*YDGSGEOM_NB%GM(1:NGPTOT)

! GMV:
DO JLEV = 1,NFLEVG
  PTL  (1:NGPTOT,JLEV)=PTL  (1:NGPTOT,JLEV)*YDGSGEOM_NB%GM(1:NGPTOT)
  PTM  (1:NGPTOT,JLEV)=PTM  (1:NGPTOT,JLEV)*YDGSGEOM_NB%GM(1:NGPTOT)
  PSPDL(1:NGPTOT,JLEV)=PSPDL(1:NGPTOT,JLEV)*YDGSGEOM_NB%GM(1:NGPTOT)
  PSPDM(1:NGPTOT,JLEV)=PSPDM(1:NGPTOT,JLEV)*YDGSGEOM_NB%GM(1:NGPTOT)
ENDDO ! JLEV

! GFL:
DO JLEV = 1,NFLEVG
  IF (YDGFL%YQ%LSP) THEN
    PQM(1:NGPTOT,JLEV)=PQM(1:NGPTOT,JLEV)*YDGSGEOM_NB%GM(1:NGPTOT)
    PQL(1:NGPTOT,JLEV)=PQL(1:NGPTOT,JLEV)*YDGSGEOM_NB%GM(1:NGPTOT)
  ENDIF
  IF (YDGFL%YL%LSP) THEN
    PLM(1:NGPTOT,JLEV)=PLM(1:NGPTOT,JLEV)*YDGSGEOM_NB%GM(1:NGPTOT)
    PLL(1:NGPTOT,JLEV)=PLL(1:NGPTOT,JLEV)*YDGSGEOM_NB%GM(1:NGPTOT)
  ENDIF
  IF (YDGFL%YI%LSP) THEN
    PIM(1:NGPTOT,JLEV)=PIM(1:NGPTOT,JLEV)*YDGSGEOM_NB%GM(1:NGPTOT)
    PIL(1:NGPTOT,JLEV)=PIL(1:NGPTOT,JLEV)*YDGSGEOM_NB%GM(1:NGPTOT)
  ENDIF
ENDDO ! JLEV

LLNHQE_BALANCE=.TRUE.

IF ( LDMODEL_TO_FILE ) THEN

  !* 1.1 TRANSFORM FROM MODEL TO FILE

  ! loop through sections of length NPROMA
  DO JIST=1,NGPTOT,NPROMA

    IST    = JIST                      ! first element
    IEND   = MIN(IST+NPROMA-1,NGPTOT)  ! last element
    IPROMA = IEND-IST+1                ! actual number of elements
    IBL=(JIST-1)/NPROMA+1

    !* 1.1.1 COMPUTE PRESSURES

    ZPREH(1:IPROMA,NFLEVG)=EXP(PSP(IST:IEND))
    CALL GPHPRE(IPROMA,NFLEVG,1,IPROMA,YDVAB,ZPREH(1:IPROMA,:),PXYB=ZXYB(1:IPROMA,:,:),&
      PRESF=ZPREF(1:IPROMA,:),LDELP=.FALSE.,LALPHA=.FALSE.,LRTGR=.FALSE.,LRPP=.FALSE.)

    !* 1.1.2 COMPUTE "R" AND "kap=R/Cp".

    CALL GPRCP_QLIRSG(IPROMA,1,IPROMA,NFLEVG,&
     & PQ=PQ(IST:IEND,:),PQI=PI(IST:IEND,:),PQL=PL(IST:IEND,:),&
     & PQR=PR(IST:IEND,:),PQS=PS(IST:IEND,:),PQG=PG(IST:IEND,:),&
     & PR=ZR(1:IPROMA,:),PKAP=ZKAP(1:IPROMA,:))

    !* 1.1.3 COMPUTE "R Tt".

    DO JLEV=1,NFLEVG
      ZRT(1:IPROMA,JLEV)=ZR(1:IPROMA,JLEV)*PT(IST:IEND,JLEV)
    ENDDO

    !* 1.1.4 COMPUTE "exp((R/cp)Qcha)".
    CALL GNHPRE(NPDVAR,YDGEOMETRY,IPROMA,NFLEVG,1,IPROMA,PSPD(IST:IEND,:),ZPREF(1:IPROMA,:), &
     & PKAP=ZKAP(1:IPROMA,:),PEQCHAF=ZEQCHAF(1:IPROMA,:))

    !* 1.1.5 COMPUTE TERM NHX

    CALL GNHQE_NHX(YDGEOMETRY,IPROMA,1,IPROMA,YDOROG(IBL)%OROG,&
     & ZOROGL(IST:IEND),ZOROGM(IST:IEND),&
     & PSP(IST:IEND),PSPL(IST:IEND),PSPM(IST:IEND),&
     & PT(IST:IEND,:),PTL(IST:IEND,:),PTM(IST:IEND,:),&
     & PQ(IST:IEND,:),PQL(IST:IEND,:),PQM(IST:IEND,:),&
     & PL(IST:IEND,:),PLL(IST:IEND,:),PLM(IST:IEND,:),&
     & PI(IST:IEND,:),PIL(IST:IEND,:),PIM(IST:IEND,:),&
     & PR(IST:IEND,:),PS(IST:IEND,:),PG(IST:IEND,:),&
     & PU(IST:IEND,:),PV(IST:IEND,:),PDIV(IST:IEND,:),&
     & ZNHX(1:IPROMA,:),PPIH=ZPREH(1:IPROMA,:))

    !* 1.1.6 TRANSFORM PD: compute "pre-prehyd"

    DO JLEV=1,NFLEVG
      PSPD(IST:IEND,JLEV)=( EXP(PSPD(IST:IEND,JLEV))-1.0_JPRB )*ZPREF(1:IPROMA,JLEV)
    ENDDO

    !* 1.1.7 TRANSFORM VD: compute -G.dw

    DO JLEV=1,NFLEVG
      PSVD(IST:IEND,JLEV) = (PSVD(IST:IEND,JLEV)-ZNHX(1:IPROMA,JLEV))*&
       & ZRT(1:IPROMA,JLEV)*ZXYB(1:IPROMA,JLEV,YYTXYB%M_LNPR)
    ENDDO

    !* 1.1.8 TRANSFORM T: compute T from Tt

    DO JLEV=1,NFLEVG
      PT(IST:IEND,JLEV)=PT(IST:IEND,JLEV)*ZEQCHAF(1:IPROMA,JLEV)
    ENDDO

  ENDDO  ! JIST = 1,NGPTOT,NPROMA

ELSEIF (.NOT.LLNHQE_BALANCE) THEN

  !* 1.2 TRANSFORM FROM FILE TO MODEL

  ! for this option file 'pre-prehyd' and '-G dw' is taken into account,
  ! but we do not guarantee that S_kap D + d4 = 0.

  ! loop through sections of length NPROMA
  DO JIST=1,NGPTOT,NPROMA

    IST    = JIST                      ! first element
    IEND   = MIN(IST+NPROMA-1,NGPTOT)  ! last element
    IPROMA = IEND-IST+1                ! actual number of elements
    IBL=(JIST-1)/NPROMA+1

    !* 1.2.1 COMPUTE PRESSURES

    ZPREH(1:IPROMA,NFLEVG)=EXP(PSP(IST:IEND))
    CALL GPHPRE(IPROMA,NFLEVG,1,IPROMA,YDVAB,ZPREH(1:IPROMA,:),PXYB=ZXYB(1:IPROMA,:,:),&
      PRESF=ZPREF(1:IPROMA,:),LDELP=.FALSE.,LALPHA=.FALSE.,LRPP=.FALSE.)

    !* 1.2.2 Transform PD

    DO JLEV=1,NFLEVG

      ! save "pre-prehyd" in ZPDEP 
      ZPDEP(1:IPROMA,JLEV)=PSPD(IST:IEND,JLEV)
      PSPD(IST:IEND,JLEV) = LOG( PSPD(IST:IEND,JLEV)/ZPREF(1:IPROMA,JLEV)+1.0_JPRB )  

      ! For the horizontal derivatives of "spd" the formula used is:
      !  dQcha = [d(pre-prehyd) - (pre-prehyd) (dprehyd/prehyd)]/[(pre-prehyd)+prehyd]
      PSPDL(IST:IEND,JLEV)=&
       & ( PSPDL(IST:IEND,JLEV)&
       & - ZPDEP(1:IPROMA,JLEV)*ZXYB(1:IPROMA,JLEV,YYTXYB%M_RTGR)*PSPL(IST:IEND) )&
       & /(ZPDEP(1:IPROMA,JLEV)+ZPREF(1:IPROMA,JLEV))
      PSPDM(IST:IEND,JLEV)=&
       & ( PSPDM(IST:IEND,JLEV)&
       & - ZPDEP(1:IPROMA,JLEV)*ZXYB(1:IPROMA,JLEV,YYTXYB%M_RTGR)*PSPM(IST:IEND) )&
       & /(ZPDEP(1:IPROMA,JLEV)+ZPREF(1:IPROMA,JLEV))

    ENDDO ! JLEV

    !* 1.2.3 COMPUTE "R" AND "kap=R/Cp".

    CALL GPRCP_QLIRSG(IPROMA,1,IPROMA,NFLEVG,&
     & PQ=PQ(IST:IEND,:),PQI=PI(IST:IEND,:),PQL=PL(IST:IEND,:),&
     & PQR=PR(IST:IEND,:),PQS=PS(IST:IEND,:),PQG=PG(IST:IEND,:),&
     & PR=ZR(1:IPROMA,:),PKAP=ZKAP(1:IPROMA,:))

    !* 1.2.4 COMPUTE "exp((R/cp)Qcha)".
    CALL GNHPRE(NPDVAR,YDGEOMETRY,IPROMA,NFLEVG,1,IPROMA,PSPD(IST:IEND,:),ZPREF(1:IPROMA,:), &
     & PKAP=ZKAP(1:IPROMA,:),PEQCHAF=ZEQCHAF(1:IPROMA,:))

    !* 1.2.5 TRANSFORM T: compute Tt from T

    DO JLEV=1,NFLEVG
      PT(IST:IEND,JLEV)=PT(IST:IEND,JLEV)/ZEQCHAF(1:IPROMA,JLEV)
    ENDDO

    !* 1.2.6 TRANSFORM grad(T): compute grad(Tt) from grad(T) with some approximations

    DO JLEV=1,NFLEVG
      PTL(IST:IEND,JLEV)=PTL(IST:IEND,JLEV)/ZEQCHAF(1:IPROMA,JLEV) &
       & -ZKAP(1:IPROMA,JLEV)*PT(IST:IEND,JLEV)*PSPDL(IST:IEND,JLEV)
      PTM(IST:IEND,JLEV)=PTM(IST:IEND,JLEV)/ZEQCHAF(1:IPROMA,JLEV) &
       & -ZKAP(1:IPROMA,JLEV)*PT(IST:IEND,JLEV)*PSPDM(IST:IEND,JLEV)
    ENDDO

    !* 1.2.7 COMPUTE "R Tt".

    DO JLEV=1,NFLEVG
      ZRT(1:IPROMA,JLEV)=ZR(1:IPROMA,JLEV)*PT(IST:IEND,JLEV)
    ENDDO

    !* 1.2.8 Transform VD: compute "dver"

    DO JLEV=1,NFLEVG
      PSVD(IST:IEND,JLEV)=PSVD(IST:IEND,JLEV)/( ZRT(1:IPROMA,JLEV)*ZXYB(1:IPROMA,JLEV,YYTXYB%M_LNPR) )
    ENDDO

    !* 1.2.9 COMPUTE TERM NHX THEN ADD IT TO dver

    CALL GNHQE_NHX(YDGEOMETRY,IPROMA,1,IPROMA,YDOROG(IBL)%OROG,&
     & ZOROGL(IST:IEND),ZOROGM(IST:IEND),&
     & PSP(IST:IEND),PSPL(IST:IEND),PSPM(IST:IEND),&
     & PT(IST:IEND,:),PTL(IST:IEND,:),PTM(IST:IEND,:),&
     & PQ(IST:IEND,:),PQL(IST:IEND,:),PQM(IST:IEND,:),&
     & PL(IST:IEND,:),PLL(IST:IEND,:),PLM(IST:IEND,:),&
     & PI(IST:IEND,:),PIL(IST:IEND,:),PIM(IST:IEND,:),&
     & PR(IST:IEND,:),PS(IST:IEND,:),PG(IST:IEND,:),&
     & PU(IST:IEND,:),PV(IST:IEND,:),PDIV(IST:IEND,:),&
     & ZNHX(1:IPROMA,:),PPIH=ZPREH(1:IPROMA,:))

    DO JLEV=1,NFLEVG
      PSVD(IST:IEND,JLEV)=PSVD(IST:IEND,JLEV)+ZNHX(1:IPROMA,JLEV)
      PSNHX(IST:IEND,JLEV) = ZNHX(1:IPROMA,JLEV)
    ENDDO ! JLEV

  ENDDO ! JIST

ELSEIF (LLNHQE_BALANCE) THEN

  !* 1.3 TRANSFORM FROM FILE TO MODEL

  ! For this option file 'pre-prehyd' and '-G dw' are ignored,
  ! we assume d4 = - S_kap D, pre-prehyd=0, Tt=T.
  ! This is the only case where PDIV is needed.

  ! loop through sections of length NPROMA
  DO JIST=1,NGPTOT,NPROMA

    IST    = JIST                      ! first element
    IEND   = MIN(IST+NPROMA-1,NGPTOT)  ! last element
    IPROMA = IEND-IST+1                ! actual number of elements
    IBL=(JIST-1)/NPROMA+1

    !* 1.3.1 COMPUTE PRESSURES

    ZPREH(1:IPROMA,NFLEVG)=EXP(PSP(IST:IEND))
    CALL GPHPRE(IPROMA,NFLEVG,1,IPROMA,YDVAB,ZPREH(1:IPROMA,:),PXYB=ZXYB(1:IPROMA,:,:),&
      PRESF=ZPREF(1:IPROMA,:),LALPHA=.NOT.LVERTFE,LRTGR=.FALSE.,LRPP=.FALSE.)

    !* 1.3.2 Transform PD (assume pre-prehyd=0).

    DO JLEV=1,NFLEVG
      PSPD(IST:IEND,JLEV) = 0.0_JPRB
      PSPDL(IST:IEND,JLEV) = 0.0_JPRB
      PSPDM(IST:IEND,JLEV) = 0.0_JPRB
    ENDDO ! JLEV

    !* 1.3.8 Transform VD: compute "dver" using relationship dver=-S_kap D

    ! * compute d4 = - S_kap DIV.
    CALL GPSKAP(YDGEOMETRY,.TRUE.,IPROMA,1,IPROMA,ZXYB(1:IPROMA,:,:),PDIV(IST:IEND,1:NFLEVG),PSVD(IST:IEND,1:NFLEVG))
    PSVD(IST:IEND,1:NFLEVG)=-PSVD(IST:IEND,1:NFLEVG)

    !* 1.3.9 COMPUTE TERM NHX

    CALL GNHQE_NHX(YDGEOMETRY,IPROMA,1,IPROMA,YDOROG(IBL)%OROG,&
     & ZOROGL(IST:IEND),ZOROGM(IST:IEND),&
     & PSP(IST:IEND),PSPL(IST:IEND),PSPM(IST:IEND),&
     & PT(IST:IEND,:),PTL(IST:IEND,:),PTM(IST:IEND,:),&
     & PQ(IST:IEND,:),PQL(IST:IEND,:),PQM(IST:IEND,:),&
     & PL(IST:IEND,:),PLL(IST:IEND,:),PLM(IST:IEND,:),&
     & PI(IST:IEND,:),PIL(IST:IEND,:),PIM(IST:IEND,:),&
     & PR(IST:IEND,:),PS(IST:IEND,:),PG(IST:IEND,:),&
     & PU(IST:IEND,:),PV(IST:IEND,:),PDIV(IST:IEND,:),&
     & ZNHX(1:IPROMA,:),PPIH=ZPREH(1:IPROMA,:))

    DO JLEV=1,NFLEVG
      PSNHX(IST:IEND,JLEV) = ZNHX(1:IPROMA,JLEV)
    ENDDO ! JLEV

  ENDDO ! JIST

ENDIF  ! LDMODEL_TO_FILE

!------------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHQE_CONV_NHVAR',1,ZHOOK_HANDLE)

END SUBROUTINE GNHQE_CONV_NHVAR
