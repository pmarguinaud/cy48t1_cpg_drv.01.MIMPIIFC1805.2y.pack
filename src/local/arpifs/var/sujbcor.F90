SUBROUTINE SUJBCOR(YDGEOMETRY,CDJBTYPE,KPRT,PCOVTPS,PCOVP,PCOVU,PCOVD,PCOVQ,&
 & PCOVO3,KLEV,KNUMB,YD_JB_STRUCT)  

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMJG    , ONLY : TYPE_JB_STRUCT
USE YOM_GRIB_CODES  , ONLY : NGRBVO, NGRBD, NGRBT, NGRBQ, NGRBO3, NGRBLNSP
USE YOMCST   , ONLY : RPI, RA
USE YOMLUN   , ONLY : NULOUT, NULERR, RESERVE_LUN, FREE_LUN
USE YOMMP0   , ONLY : MYPROC, MYSETW

!**** *SUJBCOR* - Prepare spectral correlations for background constraint Jb

!     Purpose. Prepare the error covariance model.
!     --------

!     Interface. CALL SUJBCOR(...)
!     ----------

!     Explicit arguments : In : CDJBTYPE = config string for Jb 
!     -------------------- In : KPRT = print level
!                          In : PCOVx   = covariances for parameter x

!     Implicit arguments : common YOMJG and related files.
!     --------------------

!     Method.
!     -------

!     Externals: see below.
!     ----------

!     Author : Francois Bouttier *ECMWF*  96-02-07
!     --------

!     Modifications :
!     ---------------
!      E.Holm     : 01-04-03  Move compact support SUJBCOSU to SUJBDAT
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!      M.Fisher   : 05-08-16  FGSCNM is no longer a structure.
!      K.Yessad   : Nov 2006  Check that the matrix entering EIGENMD is positive definite.
!      L.Berre    : Nov 2010 Fix for the bench
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!      P.Marguinaud (Nov 2012): Fix unallocated array arguments
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KNUMB 
CHARACTER(LEN=10)   ,INTENT(IN)    :: CDJBTYPE 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KPRT 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PCOVTPS(KLEV+1,KLEV+1,0:KNUMB-1) 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PCOVP(KLEV,KLEV,0:KNUMB-1) 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PCOVU(KLEV,KLEV,0:KNUMB-1) 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PCOVD(KLEV,KLEV,0:KNUMB-1) 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PCOVQ(KLEV,KLEV,0:KNUMB-1) 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PCOVO3(KLEV,KLEV,0:KNUMB-1) 
TYPE(TYPE_JB_STRUCT),INTENT(INOUT) :: YD_JB_STRUCT

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IPTX(10)

REAL(KIND=JPRB) :: ZCORTPS(YDGEOMETRY%YRDIMV%NFLEVG+1,YDGEOMETRY%YRDIMV%NFLEVG+1),ZCORO3(KLEV,KLEV),&
 & ZCORP(KLEV,KLEV),ZCORU(KLEV,KLEV),&
 & ZCORD(KLEV,KLEV),ZCORQ(KLEV,KLEV)  
REAL(KIND=JPRB) :: ZEIGVAL(YDGEOMETRY%YRDIMV%NFLEVG),ZEIGVAL2(YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZTCOVTPS(YDGEOMETRY%YRDIMV%NFLEVG+1,YDGEOMETRY%YRDIMV%NFLEVG+1),&
 & ZTCOVO3(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZTCOVP(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG),ZTCOVU(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZTCOVD(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG),ZTCOVQ(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)  
REAL(KIND=JPRB) :: ZSPECM(YDGEOMETRY%YRDIMV%NFLEVG,0:KNUMB-1),ZSPECD (YDGEOMETRY%YRDIMV%NFLEVG,0:KNUMB-1),&
 & ZSPECU(YDGEOMETRY%YRDIMV%NFLEVG,0:KNUMB-1),ZSPECQ (YDGEOMETRY%YRDIMV%NFLEVG,0:KNUMB-1),&
 & ZSPECT(YDGEOMETRY%YRDIMV%NFLEVG,0:KNUMB-1),ZSPECO3(YDGEOMETRY%YRDIMV%NFLEVG,0:KNUMB-1),&
 & ZSPECPS(0:KNUMB-1)  
REAL(KIND=JPRB) :: ZSPEC(YDGEOMETRY%YRDIMV%NFLEVG,0:YDGEOMETRY%YRDIM%NSMAX),&
 & ZSPECS(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIM%NSPEC2G),ZSPLSC(YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZREEL(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)  
REAL(KIND=JPRB) :: ZTCORTPS(YDGEOMETRY%YRDIMV%NFLEVG+1,YDGEOMETRY%YRDIMV%NFLEVG+1),&
 & ZTCORO3(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZTCORP(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG),ZTCORU(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZTCORQ(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG),ZTCORD(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)  
REAL(KIND=JPRB), ALLOCATABLE :: ZFGEGVNS (:,:), ZFGEDVNS (:,:,:),&
 & ZFGEGVND (:,:), ZFGEDVND (:,:,:),&
 & ZFGEGVNQ (:,:), ZFGEDVNQ (:,:,:),&
 & ZFGEGVNU (:,:), ZFGEDVNU (:,:,:),&
 & ZFGEGVNO3(:,:), ZFGEDVNO3(:,:,:)  

INTEGER(KIND=JPIM) :: IFIL, IFLEV, INLOC,&
 & IVARMAX, IPTJB, JFIELD,&
 & J1, J2, JFIL, JGL, JJ, JLEV, JN  
INTEGER(KIND=JPIM) :: IULTMP

LOGICAL :: LLBAL96TOT,LLOZJB, LL2D
LOGICAL :: LLEVI,LLEVN,LLEVZ

REAL(KIND=JPRB) :: ZSUM, ZX, ZZZ
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "eigenmd.intfb.h"
#include "gathereigmd.intfb.h"
#include "speree.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUJBCOR',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,    &
& YDLAP=>YDGEOMETRY%YRLAP, YDCSGLEG=>YDGEOMETRY%YRCSGLEG)
ASSOCIATE(NDGLG=>YDDIM%NDGLG, NSMAX=>YDDIM%NSMAX, NUMP=>YDDIM%NUMP,   NFLEVG=>YDDIMV%NFLEVG,   NSTAGP=>YDGEM%NSTAGP,   &
& RLAPDI=>YDLAP%RLAPDI,   NPROCM=>YDMP%NPROCM)

!     -----------------------------------------------------------------

!        0. Preliminary processing

LLBAL96TOT=(CDJBTYPE=='STABAL96').OR.(CDJBTYPE=='TOTENRGY')
LLOZJB=ANY(YD_JB_STRUCT%SPJB_VARS_INFO(:)%IGRIBCODE==NGRBO3)

WRITE(NULOUT,*) 'SUJBCOR: Initializing Jb correlations, ','config=',CDJBTYPE
IF (LLBAL96TOT) THEN
  IF((KLEV /= NFLEVG).OR.(KNUMB-1 /= NSMAX))&
   & CALL ABOR1(' SUJBCOR: BAD INTERFACE TO SUJBCOR')  
ELSE
  CALL ABOR1(' SUJBCOR: Unknown Jb correlation config')
ENDIF

!    ------------------------------------------------------------------

CALL GSTATS(1928,0)
!        1. Compute the vertical correlation matrices

ALLOCATE (ZFGEGVNS(NFLEVG,NUMP))
ALLOCATE (ZFGEDVNS(NFLEVG,NFLEVG,NUMP))
ALLOCATE (ZFGEGVND(NFLEVG,NUMP))
ALLOCATE (ZFGEDVND(NFLEVG,NFLEVG,NUMP))
ALLOCATE (ZFGEGVNQ(NFLEVG,NUMP))
ALLOCATE (ZFGEDVNQ(NFLEVG,NFLEVG,NUMP))
ALLOCATE (ZFGEGVNU(NFLEVG+1,NUMP))
ALLOCATE (ZFGEDVNU(NFLEVG+1,NFLEVG+1,NUMP))

IF (LLOZJB) THEN
  ALLOCATE (ZFGEGVNO3(NFLEVG,NUMP))
  ALLOCATE (ZFGEDVNO3(NFLEVG,NFLEVG,NUMP))
ENDIF

INLOC=0
DO JN=0,KNUMB-1

! Work is distributed among PE's in the same way as for 
! spectral wavenumbers

  IF (MYSETW/=NPROCM(JN)) CYCLE
  INLOC=INLOC+1

!         1.1 Normalize covariances

  DO J1=1,NFLEVG
#ifdef VPP
! Is the next directive really relevent ?? REK
!OCL SCALAR
#endif
     DO J2=1,NFLEVG
       ZCORP(J1,J2) = PCOVP (J1,J2,JN)/SQRT(PCOVP(J1,J1,JN)*PCOVP(J2,J2,JN))
       ZCORD(J1,J2) = PCOVD (J1,J2,JN)/SQRT(PCOVD(J1,J1,JN)*PCOVD(J2,J2,JN))
       ZCORU(J1,J2) = PCOVU(J1,J2,JN)/SQRT(PCOVU(J1,J1,JN)*PCOVU(J2,J2,JN))
       ZCORQ(J1,J2) = PCOVQ (J1,J2,JN)/SQRT(PCOVQ(J1,J1,JN)*PCOVQ(J2,J2,JN))
     ENDDO

  ENDDO
  DO J1=1,NFLEVG+1
#ifdef VPP
! Is the next directive really relevent ?? REK
!OCL SCALAR
#endif
    DO J2=1,NFLEVG+1
      ZCORTPS(J1,J2) = PCOVTPS(J1,J2,JN)&
       & /SQRT(PCOVTPS(J1,J1,JN)*PCOVTPS(J2,J2,JN))  
    ENDDO
  ENDDO
  IF(LLOZJB)THEN
    DO J1=1,NFLEVG
#ifdef VPP
! Is the next directive really relevent ?? REK
!OCL SCALAR
#endif
      DO J2=1,NFLEVG
        ZCORO3(J1,J2) = PCOVO3(J1,J2,JN)&
         & /SQRT(PCOVO3(J1,J1,JN)*PCOVO3(J2,J2,JN))  
      ENDDO
    ENDDO
  ENDIF

!          1.2  Load Eigenvalues and eigenvectors of inverse correl matrices
!               into Jb work arrays

  CALL EIGENMD(ZCORP,ZFGEDVNS(1,1,INLOC),ZEIGVAL,NFLEVG,LLEVI,LLEVN,LLEVZ)
  IF (LLEVI) THEN
    WRITE(NULERR,*) ' SUJBCOR: ZCORP has at least one not real eigenvalue.'
    WRITE(NULERR,*) ' for total wavenumber =',JN
  ENDIF
  IF (LLEVN) THEN
    WRITE(NULERR,*) ' SUJBCOR: ZCORP has at least one real eigenvalue < 0.'
    WRITE(NULERR,*) ' for total wavenumber =',JN
  ENDIF
  IF (LLEVZ) THEN
    WRITE(NULERR,*) ' SUJBCOR: ZCORP has at least one real eigenvalue = 0.'
    WRITE(NULERR,*) ' for total wavenumber =',JN
  ENDIF
  IF (LLEVI.OR.LLEVN.OR.LLEVZ) CALL FLUSH(NULERR)
  IF (LLEVI.OR.LLEVN.OR.LLEVZ) CALL ABOR1(' SUJBCOR: ABOR1 nr 1.2a')
  ZFGEGVNS(:,INLOC)=1.0_JPRB/SQRT(ZEIGVAL(:))

  CALL EIGENMD(ZCORD,ZFGEDVND(1,1,INLOC),ZEIGVAL,NFLEVG,LLEVI,LLEVN,LLEVZ)
  IF (LLEVI) THEN
    WRITE(NULERR,*) ' SUJBCOR: ZCORD has at least one not real eigenvalue.'
    WRITE(NULERR,*) ' for total wavenumber =',JN
  ENDIF
  IF (LLEVN) THEN
    WRITE(NULERR,*) ' SUJBCOR: ZCORD has at least one real eigenvalue < 0.'
    WRITE(NULERR,*) ' for total wavenumber =',JN
  ENDIF
  IF (LLEVZ) THEN
    WRITE(NULERR,*) ' SUJBCOR: ZCORD has at least one real eigenvalue = 0.'
    WRITE(NULERR,*) ' for total wavenumber =',JN
  ENDIF
  IF (LLEVI.OR.LLEVN.OR.LLEVZ) CALL FLUSH(NULERR)
  IF (LLEVI.OR.LLEVN.OR.LLEVZ) CALL ABOR1(' SUJBCOR: ABOR1 nr 1.2b')
  ZFGEGVND(:,INLOC)=1.0_JPRB/SQRT(ZEIGVAL(:))

  CALL EIGENMD(ZCORQ,ZFGEDVNQ(1,1,INLOC),ZEIGVAL,NFLEVG,LLEVI,LLEVN,LLEVZ)
  IF (LLEVI) THEN
    WRITE(NULERR,*) ' SUJBCOR: ZCORQ has at least one not real eigenvalue.'
    WRITE(NULERR,*) ' for total wavenumber =',JN
  ENDIF
  IF (LLEVN) THEN
    WRITE(NULERR,*) ' SUJBCOR: ZCORQ has at least one real eigenvalue < 0.'
    WRITE(NULERR,*) ' for total wavenumber =',JN
  ENDIF
  IF (LLEVZ) THEN
    WRITE(NULERR,*) ' SUJBCOR: ZCORQ has at least one real eigenvalue = 0.'
    WRITE(NULERR,*) ' for total wavenumber =',JN
  ENDIF
  IF (LLEVI.OR.LLEVN.OR.LLEVZ) CALL FLUSH(NULERR)
  IF (LLEVI.OR.LLEVN.OR.LLEVZ) CALL ABOR1(' SUJBCOR: ABOR1 nr 1.2c')
  ZFGEGVNQ(:,INLOC)=1.0_JPRB/SQRT(ZEIGVAL(:))

  CALL EIGENMD(ZCORTPS,ZFGEDVNU(1,1,INLOC),ZEIGVAL2,NFLEVG+1,LLEVI,LLEVN,LLEVZ)
  IF (LLEVI) THEN
    WRITE(NULERR,*) ' SUJBCOR: ZCORTPS has at least one not real eigenvalue.'
    WRITE(NULERR,*) ' for total wavenumber =',JN
  ENDIF
  IF (LLEVN) THEN
    WRITE(NULERR,*) ' SUJBCOR: ZCORTPS has at least one real eigenvalue < 0.'
    WRITE(NULERR,*) ' for total wavenumber =',JN
  ENDIF
  IF (LLEVZ) THEN
    WRITE(NULERR,*) ' SUJBCOR: ZCORTPS has at least one real eigenvalue = 0.'
    WRITE(NULERR,*) ' for total wavenumber =',JN
  ENDIF
  IF (LLEVI.OR.LLEVN.OR.LLEVZ) CALL FLUSH(NULERR)
  IF (LLEVI.OR.LLEVN.OR.LLEVZ) CALL ABOR1(' SUJBCOR: ABOR1 nr 1.2d')
  ZFGEGVNU(:,INLOC)=1.0_JPRB/SQRT(ZEIGVAL2(:))

  IF(LLOZJB)THEN
    CALL EIGENMD(ZCORO3,ZFGEDVNO3(1,1,INLOC),ZEIGVAL,NFLEVG,LLEVI,LLEVN,LLEVZ)
    IF (LLEVI) THEN
      WRITE(NULERR,*) ' SUJBCOR: ZCORO3 has at least one not real eigenvalue.'
      WRITE(NULERR,*) ' for total wavenumber =',JN
    ENDIF
    IF (LLEVN) THEN
      WRITE(NULERR,*) ' SUJBCOR: ZCORO3 has at least one real eigenvalue < 0.'
      WRITE(NULERR,*) ' for total wavenumber =',JN
    ENDIF
    IF (LLEVZ) THEN
      WRITE(NULERR,*) ' SUJBCOR: ZCORO3 has at least one real eigenvalue = 0.'
      WRITE(NULERR,*) ' for total wavenumber =',JN
    ENDIF
    IF (LLEVI.OR.LLEVN.OR.LLEVZ) CALL FLUSH(NULERR)
    IF (LLEVI.OR.LLEVN.OR.LLEVZ) CALL ABOR1(' SUJBCOR: ABOR1 nr 1.2e')
    ZFGEGVNO3(:,INLOC)=1.0_JPRB/SQRT(ZEIGVAL(:))
  ENDIF

ENDDO
CALL GSTATS(1928,1)

! Gather global FGExxxx arrays on each PE. This routine is only coded
! to work if KNUMB-1==NSMAX, like the rest of the sujbcor routine.

IF (LLOZJB) THEN
  CALL GATHEREIGMD(YDGEOMETRY,YD_JB_STRUCT,ZFGEGVNS,ZFGEDVNS,ZFGEGVND,ZFGEDVND,ZFGEGVNQ,&
   & ZFGEDVNQ,ZFGEGVNU,ZFGEDVNU,LLOZJB,ZFGEGVNO3,ZFGEDVNO3)
ELSE
  CALL GATHEREIGMD(YDGEOMETRY,YD_JB_STRUCT,ZFGEGVNS,ZFGEDVNS,ZFGEGVND,ZFGEDVND,ZFGEGVNQ,&
   & ZFGEDVNQ,ZFGEGVNU,ZFGEDVNU,LLOZJB)
ENDIF

DEALLOCATE (ZFGEGVNS, ZFGEDVNS, ZFGEGVND, ZFGEDVND)
DEALLOCATE (ZFGEGVNQ, ZFGEDVNQ, ZFGEGVNU, ZFGEDVNU)
IF (LLOZJB) THEN
  DEALLOCATE (ZFGEGVNO3, ZFGEDVNO3)
ENDIF

!    ------------------------------------------------------------------

!        2. Compute the horizontal correlation matrices

!         2.1 Extract the horizontal variance spectra, level by level

CALL GSTATS(1929,0)
DO J1=1,NFLEVG+1
  DO J2=1,NFLEVG+1
    ZTCOVTPS(J1,J2)=SUM(PCOVTPS(J1,J2,:))
  ENDDO
ENDDO
DO J1=1,NFLEVG
  DO J2=1,NFLEVG
    ZTCOVP(J1,J2) = SUM(PCOVP (J1,J2,:))
    ZTCOVD(J1,J2) = SUM(PCOVD (J1,J2,:))
    ZTCOVU(J1,J2) = SUM(PCOVU (J1,J2,:))
    ZTCOVQ(J1,J2) = SUM(PCOVQ (J1,J2,:))
  ENDDO
ENDDO
IF (LLOZJB) THEN
  DO J1=1,NFLEVG
    DO J2=1,NFLEVG
      ZTCOVO3(J1,J2) = SUM(PCOVO3(J1,J2,:))
    ENDDO
  ENDDO
ENDIF

DO JJ=1,NFLEVG
  DO JN=0,KNUMB-1
    ZSPECM(JJ,JN) = PCOVP  (JJ,JJ,JN)/ZTCOVP(JJ,JJ)
    ZSPECD(JJ,JN) = PCOVD  (JJ,JJ,JN)/ZTCOVD(JJ,JJ)
    ZSPECU(JJ,JN) = PCOVU  (JJ,JJ,JN)/ZTCOVU(JJ,JJ)
    ZSPECQ(JJ,JN) = PCOVQ  (JJ,JJ,JN)/ZTCOVQ(JJ,JJ)
    ZSPECT(JJ,JN) = PCOVTPS(JJ,JJ,JN)/ZTCOVTPS(JJ,JJ)
  ENDDO
ENDDO
IF (LLOZJB) THEN
  DO JJ=1,NFLEVG
    DO JN=0,KNUMB-1
      ZSPECO3(JJ,JN) = PCOVO3 (JJ,JJ,JN)/ZTCOVO3(JJ,JJ)
    ENDDO
  ENDDO
ENDIF
ZSPECPS(:)=PCOVTPS(NFLEVG+1,NFLEVG+1,:)/ZTCOVTPS(NFLEVG+1,NFLEVG+1)

!     Loop on variables :
WRITE(NULOUT,*) ' Variable index definition in Jb :'
WRITE(NULOUT,*) '     1 : vorticity'
WRITE(NULOUT,*) '     2 : unbalanced divergence'
WRITE(NULOUT,*) '     4 : specific humidity'
WRITE(NULOUT,*) '     5 : unbalanced temperature'
WRITE(NULOUT,*) '     6 : surface pressure '
IVARMAX=6
IF(LLOZJB)THEN
  WRITE(NULOUT,*) '     7 : unbalanced ozone '
  IVARMAX=7
ENDIF
!     Variable pointers (note variable 6 has only one level)
DO JFIL=1,6
  IPTX(JFIL)=(JFIL-1)*NFLEVG
ENDDO
IPTX(7)=5*NFLEVG+1

HORIZ_LOOP : DO JFIL=1,IVARMAX

!         2.2 Prepare horizontal correlation spectra for all variables

  ZSPEC(:,:)=0.0_JPRB
  DO JN=0,NSMAX
    IF (JFIL==1) ZSPEC(:,JN)=ZSPECM (:,JN)
    IF (JFIL==2) ZSPEC(:,JN)=ZSPECD (:,JN)
    IF (JFIL==3) ZSPEC(:,JN)=ZSPECT (:,JN)
    IF (JFIL==6) ZSPEC(1,JN)=ZSPECPS(JN)
    IF (JFIL==7) ZSPEC(:,JN)=ZSPECO3(:,JN)
    IF (JFIL==4) ZSPEC(:,JN)=ZSPECQ(:,JN)
    IF (JFIL==5) ZSPEC(:,JN)=ZSPECT(:,JN)
  ENDDO

!          2.3 Load inverse square root of the horizontal correlation
!              spectra into Jb work arrays - remember the 2n+1 factor !

  LL2D=.FALSE.
  IPTJB=HUGE(IPTJB)
  DO JFIELD=1,YD_JB_STRUCT%CONFIG%N_SPJB_VARS
    IF (JFIL==1 .AND. YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBVO) THEN
      IPTJB=YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (JFIL==2 .AND. YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBD) THEN
      IPTJB=YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (JFIL==3 .AND. YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBT) THEN
      IPTJB=YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (JFIL==4 .AND. YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBQ) THEN
      IPTJB=YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (JFIL==5 .AND. YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBT) THEN
      IPTJB=YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (JFIL==6 .AND. YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBLNSP) THEN
      LL2D=.TRUE.
      IPTJB=YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (JFIL==7 .AND. YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBO3) THEN
      IPTJB=YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ENDIF
  ENDDO

  IF (IPTJB==HUGE(IPTJB)) THEN
    WRITE (NULOUT,*) 'FAILED TO MATCH JFIL=',JFIL,' WITH IPTJB'
    CALL ABOR1 (' SUJBCOR: FAILED TO MATCH JFIL WITH IPTJB')
  ENDIF

  DO JN=0,NSMAX
    IF (LL2D) THEN
      YD_JB_STRUCT%JB_DATA%FGSCNM2D(JN,IPTJB)=&
          & 1.0_JPRB / SQRT( ZSPEC(1,JN)/(2.0_JPRB*JN+1.0_JPRB) )
    ELSE
      YD_JB_STRUCT%JB_DATA%FGSCNM3D(:,JN,IPTJB)=&
          & 1.0_JPRB / SQRT( ZSPEC(:,JN)/(2.0_JPRB*JN+1.0_JPRB) )
    ENDIF
  ENDDO

ENDDO HORIZ_LOOP

!         2.3 R.I.P. Old separable horizontal correlations (ex-SUGERCO)

!     Separable correls setup is commented out because it can no longer
!     run owing to the suppression of SUVPECM (because ECMWF OI
!     stat-files are no longer supported). However, it may be useful
!     later if another formulation of vertical correlations is designed.
!     The parameters should be defined in SUJB with namelist NAMJB.
!               Summary of old code for GAUSSIAN formulation :
!CCC      Definition of exponential length-scales a in e(-r2/2a2)
!CC       RPORVO=RPORDI=RPORPP=RPORTT=RPORPS=800000. RPORQQ=300000.
!CCC      Conversion into spectrum using polar autocorrelation field
!CC       ZR=RA*(ACOS(YRCSGLEG%RMU(JGL)))
!CC       JGL latitude loop : ZC(JGL)=EXP(-(ZR**2)/(2*ZLENSC**2))
!CC       JGL,JLEV,JLON loop :  ZREEL(JLON+NSTAGP(JGL),JLEV) =ZC(JGL)
!CC       CALL REESPE(NFLEVG,NFLEVG,ZSPECS,ZREEL)
!CC       JLEV,JN loop : INM=1+JN*2
!CC       ZSPECS(JLEV,INM)=ZSPECS(JLEV,INM)*SQRT(2.*JN+1.)
!CCC      Extrapolate spectrum up to NCMAX in case of stretched model
!CC       for vort and div : ZSPEC(JLEV,INM)=-ZSPEC(JLEV,INM)*YRLAP%RLAPDI(JN)
!                Summary of old code for PARAMETRIC SPECTRUM formulation :
!CCC      Definition of spectrum parameters for all levels
!CC       XNCO0(:)=15 XNCO1(:)=2. XCOSP0(:)=2. XCOSP1(:)=3. VCO0(:)=.1
!CC       ZSPEC=(1./(1.+(REAL(JN)/XNCO0(JLEV))**(XCOSP0(JLEV)+XCOSP1(JLEV)) ))
!CC  S    *(VCO0(JLEV)+(REAL(JN)/XNCO1(JLEV))**XCOSP1(JLEV))
!CC       ZSPEC(JLEV,1)=0.
!CC       for mass : ZSPEC(JLEV,INM)=-ZSPEC(JLEV,INM)*RCAPIN(JN)
!CC                  ZSPEC(JLEV,1)=ZSPEC(JLEV,1+2)

!     ------------------------------------------------------------------

!         8. Diagnose implied vertical correlations

DO J2=1,NFLEVG
  DO J1=1,NFLEVG
    ZTCORP(J1,J2)=ZTCOVP(J1,J2)/SQRT(ZTCOVP(J1,J1)*ZTCOVP(J2,J2))
    ZTCORD(J1,J2)=ZTCOVD(J1,J2)/SQRT(ZTCOVD(J1,J1)*ZTCOVD(J2,J2))
    ZTCORU(J1,J2)=ZTCOVU(J1,J2)/SQRT(ZTCOVU(J1,J1)*ZTCOVU(J2,J2))
    ZTCORQ(J1,J2)=ZTCOVQ(J1,J2)/SQRT(ZTCOVQ(J1,J1)*ZTCOVQ(J2,J2))
  ENDDO
ENDDO
IF(LLOZJB)THEN
  DO J2=1,NFLEVG
    DO J1=1,NFLEVG
      ZTCORO3(J1,J2)=ZTCOVO3(J1,J2)/SQRT(ZTCOVO3(J1,J1)*ZTCOVO3(J2,J2))
    ENDDO
  ENDDO
ENDIF
DO J2=1,NFLEVG+1
  DO J1=1,NFLEVG+1
    ZTCORTPS(J1,J2)=ZTCOVTPS(J1,J2)/SQRT(ZTCOVTPS(J1,J1)*ZTCOVTPS(J2,J2))
  ENDDO
ENDDO

!            8.3 Print implied vertical correlation matrices

WRITE(NULOUT,*) ' Diagnostics of the average ','vertical correlations :'
IF (KPRT >= 1) THEN
  831 FORMAT(1X,100I3)
  832 FORMAT(/1X,100I3)
  WRITE(NULOUT,'(/35X,'' Vorticity correlations (*100) '')')
  WRITE(NULOUT,832) (JLEV, JLEV=1,NFLEVG)
  DO J1=1,NFLEVG
    WRITE(NULOUT,831) (NINT(100._JPRB*ZTCORP(J2,J1)),J2=1,NFLEVG)
  ENDDO

  WRITE(NULOUT,'(/35X,'' Unbal Div correlations (*100) '')')
  WRITE(NULOUT,832) (JLEV, JLEV=1,NFLEVG)
  DO J1=1,NFLEVG
    WRITE(NULOUT,831) (NINT(100._JPRB*ZTCORD(J2,J1)),J2=1,NFLEVG)
  ENDDO
  WRITE(NULOUT,'(/35X,'' unbalanced'')')
  WRITE(NULOUT,'(/35X,'' T,Ps field  correlations (*100) '')')
  WRITE(NULOUT,832) (JLEV, JLEV=1,NFLEVG+1)
  DO J1=1,NFLEVG+1
    WRITE(NULOUT,831) (NINT(100._JPRB*ZTCORTPS(J2,J1)),J2=1,NFLEVG+1)
  ENDDO
  IF(LLOZJB) THEN
    WRITE(NULOUT,'(/35X,'' Unbal ozone correlations (*100) '')')
    WRITE(NULOUT,832) (JLEV, JLEV=1,NFLEVG)
    DO J1=1,NFLEVG
      WRITE(NULOUT,831) (NINT(100._JPRB*ZTCORO3(J2,J1)),J2=1,NFLEVG)
    ENDDO
  ENDIF

  WRITE(NULOUT,'(/35X,'' Q humidity field correlations(*100)'')')
  WRITE(NULOUT,832) (JLEV, JLEV=1,NFLEVG)
  DO J1=1,NFLEVG
    WRITE(NULOUT,831) (NINT(100._JPRB*ZTCORQ(J2,J1)),J2=1,NFLEVG)
  ENDDO
ENDIF
CALL GSTATS(1929,1)

!            9. Print implied horizontal correlation structures

!     Warning: this will make sense only if FGSCNM contains Ch^-1/2, not Ch

WRITE(NULOUT,*) ' Diagnostics of the horizontal correlations :'
VARIABLE_LOOP : DO JFIL=1,IVARMAX
  IFLEV=NFLEVG
  IF(JFIL == 6) IFLEV=1

  LL2D=.FALSE.
  IPTJB=HUGE(IPTJB)
  DO JFIELD=1,YD_JB_STRUCT%CONFIG%N_SPJB_VARS
    IF (JFIL==1 .AND. YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBVO) THEN
      IPTJB=YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (JFIL==2 .AND. YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBD) THEN
      IPTJB=YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (JFIL==3 .AND. YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBT) THEN
      IPTJB=YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (JFIL==4 .AND. YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBQ) THEN
      IPTJB=YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (JFIL==5 .AND. YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBT) THEN
      IPTJB=YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (JFIL==6 .AND. YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBLNSP) THEN
      LL2D=.TRUE.
      IPTJB=YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (JFIL==7 .AND. YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBO3) THEN
      IPTJB=YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ENDIF
  ENDDO

  IF (IPTJB==HUGE(IPTJB)) THEN
    WRITE (NULOUT,*) 'FAILED TO MATCH JFIL=',JFIL,' WITH IPTJB'
    CALL ABOR1 (' SUJBCOR: FAILED TO MATCH JFIL WITH IPTJB')
  ENDIF


  IF (KPRT>=1) THEN
!         Print actual length scale of spectrally fitted
!         correlation function. L**2=-2*F(0)/Laplace(F)(0)
!         (The '2' should be there in the 2D-case (Daley p.110)
    DO JLEV=1,IFLEV
      IFIL=JLEV+IPTX(JFIL)
      ZSUM=0._JPRB
      DO JN=0,NSMAX
        IF (LL2D) THEN
          ZSUM=ZSUM+RLAPDI(JN)*REAL(2*JN+1,JPRB)/&
                               & (YD_JB_STRUCT%JB_DATA%FGSCNM2D(JN,IPTJB)**2)  
        ELSE
          ZSUM=ZSUM+RLAPDI(JN)*REAL(2*JN+1,JPRB)/&
                               & (YD_JB_STRUCT%JB_DATA%FGSCNM3D(JLEV,JN,IPTJB)**2)  
        ENDIF
      ENDDO
      ZSPLSC(JLEV)=SQRT( -2.0_JPRB/ ZSUM )
    ENDDO
    WRITE(NULOUT,*) 'VARIABLE ',JFIL,' CORREL LENGTH SCALES :'
    WRITE(NULOUT,FMT='(1X,7F10.1)') ZSPLSC(1:IFLEV)
  ENDIF

  IF (KPRT>=2) THEN
!         Print spectra - big output !
    WRITE(NULOUT,*) 'Printing correlation spectra for each ',&
     & ' variable & level on file jbdiag.lst (big !)'  
    IULTMP = RESERVE_LUN()
    OPEN(IULTMP,FILE='jbdiag.lst',POSITION='APPEND')
    DO JLEV=1,IFLEV
      WRITE(IULTMP,*) 'CORRELATION SPECTRUM FOR VARIABLE ',&
       & JFIL,' AT LEVEL ',JLEV,' :'  
      WRITE(IULTMP,*) ' log10(n) log10(var)  n        var'
      IFIL=JLEV+IPTX(JFIL)
      DO JN=0,NSMAX
        IF (LL2D) THEN
          ZX=1.0_JPRB/(YD_JB_STRUCT%JB_DATA%FGSCNM2D(JN,IPTJB)**2)
        ELSE
          ZX=1.0_JPRB/(YD_JB_STRUCT%JB_DATA%FGSCNM3D(JLEV,JN,IPTJB)**2)
        ENDIF
        WRITE(IULTMP,FMT='(E11.4,1X,E11.4,1X,I3,1X,E11.4)')&
         & LOG10(MAX(TINY(ZZZ),REAL(JN,JPRB))),LOG10(ZX),JN,ZX  
      ENDDO
    ENDDO
    CLOSE(IULTMP)
    CALL FREE_LUN(IULTMP)
  ENDIF

  IF(KPRT >= 2.AND.MYPROC == 1)THEN
!         Print gridpoint structure for polar autocorrelation - big output !
!         Grid-point values are displayed at truncation NSMAX because
!         SPEREE cannot work at NCMAX in stretched model.
    WRITE(NULOUT,*) 'Printing gridpoint correlations for each ',&
     & ' variable & level on file jbdiag.lst (big !)'  
    IULTMP = RESERVE_LUN()
    OPEN(IULTMP,FILE='jbdiag.lst',POSITION='APPEND')
    ZSPECS(:,:)=0.0_JPRB
    DO JLEV=1,IFLEV
      IFIL=JLEV+IPTX(JFIL)
      DO JN=0,NSMAX
        IF (LL2D) THEN
          ZSPECS(JLEV,1+JN*2)=  SQRT( REAL(2*JN+1,JPRB))&
                             & /(YD_JB_STRUCT%JB_DATA%FGSCNM2D(JN,IPTJB)**2)  
        ELSE
          ZSPECS(JLEV,1+JN*2)=  SQRT( REAL(2*JN+1,JPRB))&
                             & /(YD_JB_STRUCT%JB_DATA%FGSCNM3D(JLEV,JN,IPTJB)**2)  
        ENDIF
      ENDDO
    ENDDO
    CALL SPEREE(YDGEOMETRY,NFLEVG,IFLEV,ZSPECS,ZREEL)
    DO JLEV=1,IFLEV
      WRITE(IULTMP,*) 'GRIDPOINT CORRELATION FOR VARIABLE ',&
       & JFIL,' AT LEVEL ',JLEV,' :'  
      WRITE(IULTMP,*) '  dist(km)    correl       mu       colat'
      DO JGL=1,NDGLG
        WRITE(IULTMP,FMT='(4(1X,E11.4))')&
         & RA*ACOS(YDCSGLEG%RMU(JGL))/1000._JPRB,&
         & ZREEL(1+NSTAGP(JGL),JLEV),&
         & YDCSGLEG%RMU(JGL),ACOS(YDCSGLEG%RMU(JGL))*180._JPRB/RPI  
      ENDDO
    ENDDO
    CLOSE(IULTMP)
    CALL FREE_LUN(IULTMP)
  ENDIF
ENDDO VARIABLE_LOOP

!---------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUJBCOR',1,ZHOOK_HANDLE)
END SUBROUTINE SUJBCOR



