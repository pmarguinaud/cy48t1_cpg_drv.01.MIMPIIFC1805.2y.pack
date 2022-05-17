SUBROUTINE SUJBDAT(YDGEOMETRY,YDDIMF,YDDYN,CDJBTYPE,KPRT,PCOVTPS,PCOVP,PCOVU,PCOVD,PCOVQ,&
                 & PCOVO3,KLEV,KNUMB,YD_JB_STRUCT)

!**** *SUJBDAT* - Input data for the background constraint Jb

!     Purpose. Read and prepare data for the error covariance model.
!     --------
 
!     Interface. CALL SUJBDAT
!     ----------
 
!     Implicit arguments : common YOMJG and related files.
!     --------------------

!     Explicit arguments : In : CDJBTYPE = type of setups
!     -------------------- In : KPRT = print level
!                          Out : PCOVx   = covariances for parameter x

!     Method.
!     -------
!     0. Initialize Jb type and general-purpose data.
!     1. Read data from Jbtype-dependent files.
!     2. Interpolate covariances to vertical model geometry.
!     3. Perform a set of Jbtype-dependent monkey business like :
!           Map covariances onto model spectral representation.
!           Correct (T,Ps) covariances for hydrostatic balance.
!     4. Convert (T,Ps) into P covariances, simultaneously computing
!        the P->(T,Ps) operator.
!     5. Make covariance matrices compactly supported

!     Externals: see below.
!     ----------

!     Author : Francois Bouttier *ECMWF*  96-02-07
!     -------- using plenty of Erik Andersson's code from SUNSFCE

!     Modifications :
!     ---------------
!      01-04-02 E.Holm  : Move call to SUJBCOSU here from SUJBCOR
!                         (SUJBCOSU now does 3D compact support).
!      01-08-07 R. El Khatib : Pruning options
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE YOMDIMF      , ONLY : TDIMF
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMJG    , ONLY : TYPE_JB_STRUCT
USE YOM_GRIB_CODES  , ONLY : NGRBO3
USE YOMLUN   , ONLY : NULOUT, RESERVE_LUN, FREE_LUN
USE YOMDYN   , ONLY : TDYN
USE YOMMP0   , ONLY : NPROC, MYPROC
USE YOMCST   , ONLY : RCPD, RLVTT
USE YOMTAG   , ONLY : MTAGFCE
USE MPL_MODULE,ONLY : MPL_BROADCAST

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDIMF)         ,INTENT(INOUT) :: YDDIMF
TYPE(TDYN)          ,INTENT(INOUT) :: YDDYN
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KNUMB 
CHARACTER(LEN=10)   ,INTENT(IN)    :: CDJBTYPE 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KPRT 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PCOVTPS(KLEV+1,KLEV+1,0:KNUMB-1) 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PCOVP(KLEV,KLEV,0:KNUMB-1) 
REAL(KIND=JPRB)     ,INTENT(OUT)   :: PCOVU(KLEV,KLEV,0:KNUMB-1) 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PCOVD(KLEV,KLEV,0:KNUMB-1) 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PCOVQ(KLEV,KLEV,0:KNUMB-1) 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PCOVO3(KLEV,KLEV,0:KNUMB-1) 
TYPE(TYPE_JB_STRUCT),INTENT(IN)    :: YD_JB_STRUCT

!        Input covariance arrays
REAL(KIND=JPRB),ALLOCATABLE :: ZFCOVTPS(:,:,:),ZFCOVD(:,:,:),&
  & ZFCOVQ(:,:,:)  ,ZFCOVP(:,:,:),ZFCOVO3(:,:,:)  
 INTEGER(KIND=JPIM) :: IDIMF(2)
 
!        GSA file-specific data
CHARACTER :: CLID*10,CLCOM*70,CLFILE*40
REAL(KIND=JPRB),ALLOCATABLE :: ZPDAT(:)
 
!       Specific data for (T,Ps) adjustment and q resetting
REAL(KIND=JPRB) :: ZFAC
INTEGER(KIND=JPIM) :: ILV100
INTEGER(KIND=JPIM), PARAMETER :: JPKZ = KIND(ZFAC)

!       Specific data for (T,Ps)->(P) conversion and vertical interpolation
REAL(KIND=JPRB) :: ZSIG(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZVI(YDGEOMETRY%YRDIMV%NFLEVG+1,YDGEOMETRY%YRDIMV%NFLEVG+1)

!       Specific data for total energy weights
REAL(KIND=JPRB),ALLOCATABLE :: ZTEPCOVP(:,:,:),ZTEPCOVD(:,:,:),&
 & ZTEPCOVQ(:,:,:),ZTEPCOVTPS(:,:,:)  

INTEGER(KIND=JPIM) :: ICHKWD, ICOR1, ICOR2, IDATE, IDIM1, IDIM2,&
 & ILENDEF, IMN, INBMAT, INBSET, INLEVF, INMAXF,&
 & IOMASTER, IORIG, IPAR1, IPAR2, ISETDIST,&
 & ITIME, ITYPDI1, ITYPDI2, ITYPMAT, IWEIGHT, J1, JJ, JK, JLEV, JN  
INTEGER(KIND=JPIM) :: IULTMP

LOGICAL :: LLO3JB
REAL(KIND=JPRB) :: ZDUMMY
REAL(KIND=JPRB) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif

 
#include "abor1.intfb.h"
#include "commjbdat.intfb.h"
#include "sujbcosu.intfb.h"
#include "sujbcovsignal.intfb.h"
#include "sujbcovnoise.intfb.h"
#include "fltbgcalc.intfb.h"
! #include "suprecov.intfb.h"

!     -----------------------------------------------------------------

!        0. Describe initialization method

IF (LHOOK) CALL DR_HOOK('SUJBDAT',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDLAP=>YDGEOMETRY%YRLAP)
ASSOCIATE(NSMAX=>YDDIM%NSMAX,   NFLEVG=>YDDIMV%NFLEVG,   SIDELP=>YDDYN%SIDELP, SIRPRG=>YDDYN%SIRPRG, &
& SIRPRN=>YDDYN%SIRPRN,   SITR=>YDDYN%SITR,   RLAPIN=>YDLAP%RLAPIN)
 LLO3JB=ANY(YD_JB_STRUCT%SPJB_VARS_INFO(:)%IGRIBCODE==NGRBO3)
 
 IF((KLEV /= NFLEVG).OR.(KNUMB-1 /= NSMAX))&
  & CALL ABOR1('BAD INTERFACE TO SUJBDAT')  
 
 WRITE(NULOUT,*) 'SUJBDAT: Jb data is initialized as follows:'
 IF (CDJBTYPE=='STABAL96') THEN
   WRITE(NULOUT,*) '   covs implied by statistical balance.'
 ELSEIF (CDJBTYPE=='TOTENRGY') THEN
   WRITE(NULOUT,*) '   covs implied by total energy weights.'
 ELSE
   CALL ABOR1('Unknown Jb data config')
 ENDIF
  
 !     -----------------------------------------------------------------
 
 !        1. Read raw data
 
 !           1.0.0 Read GSA file preamble to determine data resolution
 
 IOMASTER=1
 IF (CDJBTYPE=='STABAL96') CLFILE='stabal96.cv'
 IF (MYPROC==IOMASTER) THEN
   CALL GSTATS(1705,0)
   WRITE(NULOUT,*) '  Opening input file ',CLFILE
   IULTMP = RESERVE_LUN()

#ifdef LITTLE_ENDIAN
   OPEN(IULTMP,FILE=CLFILE,FORM='unformatted',CONVERT='big_endian')
#else
   OPEN(IULTMP,FILE=CLFILE,FORM='unformatted')
#endif

   READ(IULTMP) CLID
   WRITE(NULOUT,*) 'GSA ID=',CLID
   IF (CDJBTYPE/=CLID) CALL ABOR1('Bad Jb ID in GSA file')
   READ(IULTMP) CLCOM
   WRITE(NULOUT,*) ' Description : ',CLCOM
   READ(IULTMP) IORIG,IDATE,ITIME,INBSET
   WRITE(NULOUT,*) ' Center=',IORIG,' Date=',IDATE,' Time=',ITIME
 !       GSA set 1 header : vertical level definition
   READ(IULTMP) INBMAT,IWEIGHT,ITYPMAT,ISETDIST,ILENDEF
   READ(IULTMP) IDIM1,IDIM2,IPAR1,IPAR2,ICOR1,ICOR2
   INMAXF=INBMAT-1                  ! matrices go from n=0 to itruncf
   INLEVF=IDIM1
 !       matrix dimension is inlevf except if the param is (T,lnps) :
   IF ((IPAR1==2).OR.(IPAR1==13).OR.(IPAR1==14)) INLEVF=IDIM1-1
   WRITE(NULOUT,*) ' Cov file is at nlev=',INLEVF,' nmax=',INMAXF
   IF (INLEVF/=NFLEVG)CALL ABOR1('Error - Jb covs must have NFLEVG levels')
 !         We still broadcast inlevf because one day we may like to
 !         reimplement some vertical interpolation of correlations.
 !         But for the time being it does not work. /fb
   IDIMF(1)=INMAXF
   IDIMF(2)=INLEVF
   CLOSE(IULTMP)
   CALL FREE_LUN(IULTMP)
   CALL GSTATS(1705,1)
 ENDIF
 
 IF (NPROC>1) THEN
   CALL MPL_BROADCAST(IDIMF,KROOT=IOMASTER,KTAG=MTAGFCE+146,&
     & CDSTRING='SUJBDAT')
   INMAXF=IDIMF(1)
   INLEVF=IDIMF(2)
 ENDIF
 
 !           1.0.1 Allocate raw data arrays
 
 9990 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)
 ALLOCATE(ZPDAT(NFLEVG+1))
 ALLOCATE(ZFCOVTPS(NFLEVG+1,NFLEVG+1,0:INMAXF))
 ALLOCATE(ZFCOVQ(NFLEVG,NFLEVG,0:INMAXF))
 WRITE(NULOUT,9990) 'ZPDAT     ',SIZE(ZPDAT),SHAPE(ZPDAT)
 WRITE(NULOUT,9990) 'ZFCOVTPS  ',SIZE(ZFCOVTPS),SHAPE(ZFCOVTPS)
 WRITE(NULOUT,9990) 'ZFCOVQ    ',SIZE(ZFCOVQ),SHAPE(ZFCOVQ)
 IF(CDJBTYPE=='STABAL96')THEN
   ALLOCATE(ZFCOVD (NFLEVG,NFLEVG,0:INMAXF))
   ALLOCATE(ZFCOVP (NFLEVG,NFLEVG,0:INMAXF))
   WRITE(NULOUT,9990) 'ZFCOVD    ',SIZE(ZFCOVD),SHAPE(ZFCOVD)
   WRITE(NULOUT,9990) 'ZFCOVP    ',SIZE(ZFCOVP),SHAPE(ZFCOVP)
   IF (LLO3JB) THEN
     ALLOCATE(ZFCOVO3 (NFLEVG,NFLEVG,0:INMAXF))
     WRITE(NULOUT,9990) 'ZFCOVO3   ',SIZE(ZFCOVO3),SHAPE(ZFCOVO3)
   ENDIF
 ENDIF
 
 !           1.1 Read GSA covariances on master PE
 
 CALL GSTATS(1705,0)
 IF (MYPROC==IOMASTER) THEN
 
 !         1.1.1. Open & read GSA file header
   IF (CDJBTYPE=='STABAL96') THEN
     IULTMP = RESERVE_LUN()
#ifdef LITTLE_ENDIAN
     OPEN(IULTMP,FILE=CLFILE,FORM='unformatted',CONVERT='big_endian')
#else
     OPEN(IULTMP,FILE=CLFILE,FORM='unformatted')
#endif
     READ(IULTMP) CLID
     WRITE(NULOUT,*) 'GSA ID=',CLID
     IF(CLID/=CDJBTYPE) CALL ABOR1('Bad ID in GSA file')
     READ(IULTMP) CLCOM
     WRITE(NULOUT,*) ' Description : ',CLCOM
     READ(IULTMP) IORIG,IDATE,ITIME,INBSET
     WRITE(NULOUT,*) ' Orig=',IORIG,' Date=',IDATE,' Time=',ITIME
 
 !         1.1.2. check GSA set 1 header : vertical level definition
     READ(IULTMP) INBMAT,IWEIGHT,ITYPMAT,ISETDIST,ILENDEF
     READ(IULTMP) IDIM1,IDIM2,IPAR1,IPAR2,ITYPDI1,ITYPDI2
     IF ((IDIM1/=IDIM2).OR.(IPAR1/=IPAR2)&
      & .OR.(ITYPDI1/=ITYPDI2)) CALL ABOR1('Nonsymmetric matrix')  
     WRITE(NULOUT,*) ' Matrix dimension=',IDIM1
     IF(IDIM1<=0) CALL ABOR1(' Bad matrix dimensions')
     IF(IDIM1/=NFLEVG) CALL ABOR1('code/data dim mismatch')
     READ(IULTMP)
     READ(IULTMP) (ZPDAT(JJ),JJ=1,IDIM1)
     IF(KPRT>=1) WRITE(NULOUT,*) ' Matrix input levels=',ZPDAT(1:IDIM1)
     IF(ITYPDI1/=1) CALL ABOR1('Matrix not on pressure levels')
     IF(IPAR1/=4) CALL ABOR1('Not vorticity in GSA set 1')
   ENDIF
 
 !        1.1.4 read STABAL 96 GSA covariances (except 1st header)
   IF (CDJBTYPE=='STABAL96') THEN
     WRITE(NULOUT,*) '1.Reading STABAL96 GSA covariances'
 !         GSA set 1, vort covs
     DO JN=0,INMAXF
       READ(IULTMP) ZDUMMY
       READ(IULTMP) ((ZFCOVP(JJ,JK,JN),JJ=1,NFLEVG),JK=1,NFLEVG),ICHKWD
       IF(ICHKWD /= 3141592) CALL ABOR1('Bad GSA control word 1')
     ENDDO
 !         GSA set 2, unbal div covs
     READ(IULTMP)
     READ(IULTMP) IDIM1,IDIM2,IPAR1,IPAR2,ITYPDI1,ITYPDI2
     IF(IPAR1/=12) CALL ABOR1('Not unbal div in GSA set 2')
     READ(IULTMP)
     READ(IULTMP)
     DO JN=0,INMAXF
       READ(IULTMP) ZDUMMY
       READ(IULTMP) ((ZFCOVD(JJ,JK,JN),JJ=1,NFLEVG),JK=1,NFLEVG),ICHKWD
       IF(ICHKWD /= 3141592) CALL ABOR1('Bad GSA control word 2')
     ENDDO
 !         GSA set 3, unbal (T,ps) covs
     READ(IULTMP)
     READ(IULTMP) IDIM1,IDIM2,IPAR1,IPAR2,ITYPDI1,ITYPDI2
     IF(IPAR1/=14) CALL ABOR1('Not unbal T,lnps in GSA set 3')
     IF(IDIM1/=NFLEVG+1) CALL ABOR1('code/data dim mismatch')
     READ(IULTMP)
     READ(IULTMP)
     DO JN=0,INMAXF
       READ(IULTMP) ZDUMMY
       READ(IULTMP) ((ZFCOVTPS(JJ,JK,JN),JJ=1,NFLEVG+1),JK=1,NFLEVG+1),ICHKWD
       IF(ICHKWD /= 3141592) CALL ABOR1('Bad GSA control word 3')
     ENDDO
 !         GSA set 4, q covs
     READ(IULTMP)
     READ(IULTMP) IDIM1,IDIM2,IPAR1,IPAR2,ITYPDI1,ITYPDI2
     IF(IPAR1/=3) CALL ABOR1('Not unbal q in GSA set 4')
     READ(IULTMP)
     READ(IULTMP)
     DO JN=0,INMAXF
       READ(IULTMP) ZDUMMY
       READ(IULTMP) ((ZFCOVQ(JJ,JK,JN),JJ=1,NFLEVG),JK=1,NFLEVG),ICHKWD
       IF(ICHKWD /= 3141592) CALL ABOR1('Bad GSA control word 4')
     ENDDO
 !         GSA set 5, unbal ozone covs
     IF (LLO3JB) THEN
       READ(IULTMP)
       READ(IULTMP) IDIM1,IDIM2,IPAR1,IPAR2,ITYPDI1,ITYPDI2
       IF(IPAR1/=18) CALL ABOR1('Not unbal ozone in GSA set 5')
       READ(IULTMP)
       READ(IULTMP)
       DO JN=0,INMAXF
         READ(IULTMP) ZDUMMY
         READ(IULTMP) ((ZFCOVO3(JJ,JK,JN),JJ=1,NFLEVG),JK=1,NFLEVG),ICHKWD
         IF(ICHKWD /= 3141592) CALL ABOR1('Bad GSA control word 5')
       ENDDO
     ENDIF
     CLOSE(IULTMP)
     CALL FREE_LUN(IULTMP)
   ENDIF
 
   IF(CDJBTYPE=='STABAL96') THEN
     WRITE(NULOUT,*) 'STABAL96 GSA data read ok.'
   ELSEIF(CDJBTYPE=='TOTENRGY') THEN
     WRITE(NULOUT,*) 'Initializing with total energy weights'
   ELSE
     CALL ABOR1('Unknown CDJBTYPE for data input')
   ENDIF
 ENDIF
 CALL GSTATS(1705,1)
  
 !        1.2  Broadcast covariances on all PEs
 
 IF (NPROC>1) THEN
   CALL COMMJBDAT(CDJBTYPE,NFLEVG,INMAXF,ZPDAT,&
    & ZFCOVTPS,ZFCOVQ,ZFCOVD,ZFCOVP,ZFCOVO3,YD_JB_STRUCT)
 ENDIF
 
 !        1.9 Determine model vertical geometry
 
 IF (CDJBTYPE=='STABAL96') THEN
   ! This call does not accomplish anything - commented out /MH
!!$   CALL SUPRECOV(YDGEOMETRY,YDDIMF,YDDYN,KPRT,.FALSE.,ZSIG,.FALSE.,ZVI,ZPDAT,NFLEVG+1,ILV100)
!!$   WRITE(NULOUT,*) ' 100hPa is just below level no.',ILV100
   ILV100=0    ! temporary patch for historical reasons
 ENDIF
 
 !---------------------------------------------------------------------------
 
 !       2. Load covariance arrays
 
 IF (CDJBTYPE=='STABAL96') THEN
   WRITE(NULOUT,*) ' No vertical interpolation of covariances.'
   IF (NFLEVG/=NFLEVG) CALL ABOR1(&
    & 'Vertical interpolation of Jb-STABAL is not implemented.')  
 
 !         P is vor, D is unbal div, TPS is unbal (T,Lnps), Q is q, O3 is unbal ozone.
   IMN=MIN(NSMAX,INMAXF)
   PCOVP(:,:,0:IMN)  =ZFCOVP(:,:,0:IMN)
   PCOVD(:,:,0:IMN)  =ZFCOVD(:,:,0:IMN)
   PCOVTPS(:,:,0:IMN)=ZFCOVTPS(:,:,0:IMN)
   PCOVQ(:,:,0:IMN)  =ZFCOVQ(:,:,0:IMN)
   IF (LLO3JB) PCOVO3(:,:,0:IMN)  =ZFCOVO3(:,:,0:IMN)
 ENDIF
 
 !         2.5  Prepare total energy weights if applicable
 
 IF (YD_JB_STRUCT%CONFIG%LJBENER.OR.(YD_JB_STRUCT%CONFIG%RAMENER>0.0_JPRB)) THEN
 !       fill arrays with inverse of SCALP (total energy weights)
 !       times 2n+1 times arbitrary scaling
   ALLOCATE( ZTEPCOVP(KLEV,KLEV,0:NSMAX),&
    & ZTEPCOVD(KLEV,KLEV,0:NSMAX),ZTEPCOVQ(KLEV,KLEV,0:NSMAX),&
    & ZTEPCOVTPS(KLEV+1,KLEV+1,0:NSMAX) )  
   ZTEPCOVP  (:,:,:)=0.0_JPRB
   ZTEPCOVD  (:,:,:)=0.0_JPRB
   ZTEPCOVQ  (:,:,:)=0.0_JPRB
   ZTEPCOVTPS(:,:,:)=0.0_JPRB
   DO JLEV=1,NFLEVG
     DO JN=1,NSMAX
       ZFAC = REAL(2*JN+1,JPKZ)/&
            & (YD_JB_STRUCT%CONFIG%RSCENER*(YD_JB_STRUCT%CONFIG%REDNMC**2))
       ZTEPCOVP  (JLEV,JLEV,JN)=-ZFAC/(RLAPIN(JN)*SIDELP(JLEV)/2.0_JPRB)
       ZTEPCOVD  (JLEV,JLEV,JN)=-ZFAC/(RLAPIN(JN)*SIDELP(JLEV)/2.0_JPRB)
       ZTEPCOVQ  (JLEV,JLEV,JN)= ZFAC
       ZTEPCOVTPS(JLEV,JLEV,JN)= ZFAC/(RCPD*SIDELP(JLEV)/SITR/2.0_JPRB)
     ENDDO
     IF (YD_JB_STRUCT%CONFIG%COEQTERMJB/=0) ZTEPCOVQ(JLEV,JLEV,:)=ZTEPCOVQ(JLEV,JLEV,:)/&
      & (YD_JB_STRUCT%CONFIG%COEQTERMJB*( (RLVTT**2)/(RCPD*SITR) )*SIDELP(JLEV)/2.0_JPRB)  
   ENDDO
   DO JN=0,NSMAX
     ZFAC = REAL(2*JN+1,JPKZ)/&
          & (YD_JB_STRUCT%CONFIG%RSCENER*(YD_JB_STRUCT%CONFIG%REDNMC**2))
     ZTEPCOVTPS(NFLEVG+1,NFLEVG+1,:)=ZFAC/(SIRPRG/SIRPRN/2.0_JPRB)
   ENDDO
 ENDIF
 !        2.6   Load energy coeffs or interpolate with static Jb ones
  
 IF ((CDJBTYPE=='STABAL96').AND.(YD_JB_STRUCT%CONFIG%RAMENER>0.0_JPRB)) THEN
   WRITE(NULOUT,*) '**WARNING** mixing stabal96 with energy coeffs'
   PCOVP  (:,:,:)=(1.0_JPRB-YD_JB_STRUCT%CONFIG%RAMENER)*&
                & PCOVP(:,:,:)+YD_JB_STRUCT%CONFIG%RAMENER*ZTEPCOVP(:,:,:)
   PCOVD  (:,:,:)=(1.0_JPRB-YD_JB_STRUCT%CONFIG%RAMENER)*&
                & PCOVD(:,:,:)+YD_JB_STRUCT%CONFIG%RAMENER*ZTEPCOVD(:,:,:)
   PCOVQ  (:,:,:)=(1.0_JPRB-YD_JB_STRUCT%CONFIG%RAMENER)*&
                & PCOVQ(:,:,:)+YD_JB_STRUCT%CONFIG%RAMENER*ZTEPCOVQ(:,:,:)
   PCOVTPS(:,:,:)=(1.0_JPRB-YD_JB_STRUCT%CONFIG%RAMENER)*&
                & PCOVTPS(:,:,:)+YD_JB_STRUCT%CONFIG%RAMENER*ZTEPCOVTPS(:,:,:)
 ELSEIF (CDJBTYPE=='TOTENRGY') THEN
   PCOVP  (:,:,:)=ZTEPCOVP  (:,:,:)
   PCOVD  (:,:,:)=ZTEPCOVD  (:,:,:)
   PCOVQ  (:,:,:)=ZTEPCOVQ  (:,:,:)
   PCOVTPS(:,:,:)=ZTEPCOVTPS(:,:,:)
 ENDIF
 !------------------------------------------------------------------------
 
 !       3. Miscellaneous covariance preprocessing
 
 WRITE(NULOUT,*) '3.Covariance preprocessing :'
 
 !          3.2 Wavenumber-wise extrapolations
 
 IF((CDJBTYPE=='STABAL96').OR.(CDJBTYPE=='TOTENRGY'))THEN
   IF(NSMAX>INMAXF)THEN
     WRITE(NULOUT,*) '  Extrapolate statistics beyond n=',INMAXF
     DO JN=INMAXF+1,NSMAX
       PCOVTPS(:,:,JN)=PCOVTPS(:,:,INMAXF)
        PCOVQ  (:,:,JN)=PCOVQ  (:,:,INMAXF)
     ENDDO
     IF((CDJBTYPE=='STABAL96').OR.(CDJBTYPE=='TOTENRGY'))THEN
       DO JN=INMAXF+1,NSMAX
         PCOVP  (:,:,JN)=PCOVP  (:,:,INMAXF)
         PCOVD  (:,:,JN)=PCOVD  (:,:,INMAXF)
       ENDDO
       IF (LLO3JB) THEN
         DO JN=INMAXF+1,NSMAX
           PCOVO3 (:,:,JN)=PCOVO3 (:,:,INMAXF)
         ENDDO
       ENDIF
     ENDIF
   ENDIF
 ENDIF
 
 !          3.4 Reset q correlations to zero in the stratosphere (above 100hPa)
 
 !     ILV100 has been precomputed in SUPRECOV
 DO J1=1,ILV100
   PCOVQ(J1,J1+1:NFLEVG,:)=0.0_JPRB
   PCOVQ(J1+1:NFLEVG,J1,:)=0.0_JPRB
 ENDDO
 !------------------------------------------------------------------------
 
 !       5.  Any other covariance preprocessing
 
 IF((CDJBTYPE=='STABAL96').OR.(CDJBTYPE=='TOTENRGY')) THEN
 
 !mf     Security : force covariances for vor/div at n=0 to n=1 values
 
   PCOVP(:,:,0)=PCOVP(:,:,1)
   PCOVD(:,:,0)=PCOVD(:,:,1)
 !       Dummy covariances for old mass
   PCOVU(:,:,:)=0.0_JPRB
   DO JJ=1,NFLEVG
     PCOVU(JJ,JJ,:)=1.0_JPRB
   ENDDO
 ENDIF
 
 !       5.1  Make covariance matrices compactly supported
 
 IF (YD_JB_STRUCT%CONFIG%LCORCOSU) THEN
   IF((CDJBTYPE=='STABAL96').OR.(CDJBTYPE=='TOTENRGY')) THEN
     CALL SUJBCOSU(YDGEOMETRY,YD_JB_STRUCT%CONFIG%NJBVERB,PCOVTPS,PCOVP,PCOVD,PCOVQ,PCOVO3,&
                 & SIZE(PCOVP,1),SIZE(PCOVP,3),YD_JB_STRUCT)
   ENDIF
 ENDIF

 !       5.2  Calculate an objective filter for ensemble variances
 
 IF (YD_JB_STRUCT%CONFIG%LFLTBGCALC) THEN
     CALL SUJBCOVNOISE(YDGEOMETRY,PCOVTPS,PCOVP,PCOVD,PCOVQ,&
      & SIZE(PCOVP,1),SIZE(PCOVP,3))  
     CALL SUJBCOVSIGNAL(YDGEOMETRY,YDDIMF)
     CALL FLTBGCALC(YDGEOMETRY%YRDIMV,YDDIMF)
 ENDIF

 IF (ALLOCATED(ZFCOVTPS)) DEALLOCATE(ZFCOVTPS)
 IF (ALLOCATED(ZFCOVD)) DEALLOCATE(ZFCOVD)
 IF (ALLOCATED(ZFCOVQ)) DEALLOCATE(ZFCOVQ)
 IF (ALLOCATED(ZFCOVP)) DEALLOCATE(ZFCOVP)
 IF (ALLOCATED(ZFCOVO3)) DEALLOCATE(ZFCOVO3)
 IF (ALLOCATED(ZPDAT)) DEALLOCATE(ZPDAT)
 IF (ALLOCATED(ZTEPCOVP)) DEALLOCATE(ZTEPCOVP)
 IF (ALLOCATED(ZTEPCOVD)) DEALLOCATE(ZTEPCOVD)
 IF (ALLOCATED(ZTEPCOVQ)) DEALLOCATE(ZTEPCOVQ)
 IF (ALLOCATED(ZTEPCOVTPS)) DEALLOCATE(ZTEPCOVTPS)
 
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUJBDAT',1,ZHOOK_HANDLE)
END SUBROUTINE SUJBDAT
