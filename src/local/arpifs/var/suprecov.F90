SUBROUTINE SUPRECOV(YDCST,YDGEOMETRY,YDDIMF,YDDYN,KPRT,LDSIG,PSIG,LDINTP,PVI,PRESS,KNBP,KLV100)

!**** *SUPRECOV* - Prepare covariance preprocessing operators

!     Purpose. Prepare data for the error covariance model.
!     --------

!     Interface. CALL SUPRECOV(...)
!     ----------

!     Implicit arguments : model vertical geometry in commons
!     --------------------

!     Explicit arguments : In : KPRT   = diagnostics setups
!     -------------------- In : LDSIG  = to compute SIGAM hydrostatic bal matrix
!                         Out : PSIG(NFLEVG,NFLEVG+1) = SIGAM (T,ps)->P operator
!                          In : LDINTP = to compute internal postprocessing
!                                        vertical interpolation matrix
!                         Out : PVI(NFLEVG+1,KNBP) = vertical interpol matrix
!                          In : PRESS(KNBP) = pressures to interpolate from
!                          In : KNBP  = number of p levels to interpolate from
!                         Out : KLV100 = index of lowest level above 100hPa

!     Method.
!     -------
!     1. Fill (T,Ps)->P matrix using SIGAM operator if requested
!     2. Fill (old levels)->(new levels) interpolation matrix using
!               IFS internal postprocessing operator

!     Externals: see below.
!     ----------

!     Author : Francois Bouttier *ECMWF*  96-03-19
!     -------- using plenty of Erik Andersson's code from SUNSFCE

!     Modifications :
!     ---------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad (Nov 2012): call GPHPRE.
!     ------------------------------------------------------------------

USE YOMDYN       , ONLY : TDYN
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMDIMF      , ONLY : TDIMF
USE YOMLUN       , ONLY : NULOUT
USE YOMVERT      , ONLY : VP00
USE PARDIM       , ONLY : JPNPPM
USE YOMCST       , ONLY : TCST
USE YOMCVER      , ONLY : LVERTFE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCST),        INTENT(IN)    :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDIMF)       ,INTENT(INOUT) :: YDDIMF
TYPE(TDYN)        ,INTENT(INOUT) :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KNBP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPRT 
LOGICAL           ,INTENT(IN)    :: LDSIG 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSIG(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG+1) 
LOGICAL           ,INTENT(IN)    :: LDINTP 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVI(YDGEOMETRY%YRDIMV%NFLEVG+1,KNBP) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRESS(KNBP) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLV100 

!     ------------------------------------------------------------------

!       Specific data for (T,Ps)->(P) conversion
REAL(KIND=JPRB) :: ZSP1(YDGEOMETRY%YRDIMV%NFLEVG),ZSP2(YDGEOMETRY%YRDIMV%NFLEVG),ZSPS1(1)

!       Specific data to use internal post-processing for interpolation
REAL(KIND=JPRB) :: ZPPP(KNBP,YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZPRSHI(KNBP,0:KNBP), ZPRSFI(KNBP,KNBP)
REAL(KIND=JPRB) :: ZPXP(KNBP,0:KNBP,JPNPPM), ZPXPD(KNBP,0:KNBP,JPNPPM)
REAL(KIND=JPRB) :: ZVAL(KNBP,0:KNBP), ZVALO(KNBP,YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZPRESF(YDGEOMETRY%YRDIMV%NFLEVG+1),ZPRESH(0:YDGEOMETRY%YRDIMV%NFLEVG+1)
INTEGER(KIND=JPIM) :: ILEVB(KNBP,YDGEOMETRY%YRDIMV%NFLEVG+1,JPNPPM)
LOGICAL :: LLBELO(KNBP,YDGEOMETRY%YRDIMV%NFLEVG+1),   LLBLOW(YDGEOMETRY%YRDIMV%NFLEVG+1),&
 & LLBELS(KNBP,YDGEOMETRY%YRDIMV%NFLEVG+1),   LLBLES(YDGEOMETRY%YRDIMV%NFLEVG+1)  

INTEGER(KIND=JPIM) :: J1, J2, JLEV
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gphpre.intfb.h"
#include "ppflev.intfb.h"
#include "ppinit.intfb.h"
#include "ppintp.intfb.h"
#include "sigam.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUPRECOV',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB)
ASSOCIATE( NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)
!     ------------------------------------------------------------------

!              1. Prepare the matrix of operator SIGAM

IF(LDSIG)THEN
  WRITE(NULOUT,*) '   Initializing SIGAM operator matrix'
  DO J2=1,NFLEVG
    ZSP1(:)=0.0_JPRB
    ZSP1(J2)=1.0_JPRB
    ZSPS1(1)=0._JPRB
    CALL SIGAM(LVERTFE, YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZSP2,ZSP1,ZSPS1,1,NFLEVG)
    PSIG(:,J2)=ZSP2(:)
  ENDDO
  ZSP1(:)=0.0_JPRB
  ZSPS1(1)=1._JPRB
  CALL SIGAM(LVERTFE, YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZSP2,ZSP1,ZSPS1,1,NFLEVG)
  PSIG(:,NFLEVG+1)=ZSP2(:)
  IF(KPRT >= 2)THEN
    WRITE(NULOUT,*) ' SUPRECOV: SIGAM MATRIX '
    DO J1=1,NFLEVG
      WRITE(NULOUT,*) J1
      WRITE(NULOUT,'(17(1X,F8.0))') (PSIG(J1,J2),J2=1,NFLEVG+1)
    ENDDO
  ENDIF
ENDIF

!              2. Prepare the matrix of operator PPINTP

IF(LDINTP)THEN
  WRITE(NULOUT,*) '   Initializing PPINTP operator matrix'
!       Define full and half level pressures
  DO JLEV=1,KNBP
    ZPRSFI(:,JLEV) = PRESS(JLEV)
  ENDDO
  DO JLEV=1,KNBP-1
    ZPRSHI(:,JLEV) = 0.5_JPRB*(ZPRSFI(:,JLEV)+ZPRSFI(:,JLEV+1))
  ENDDO
  ZPRSHI(:,0) = ZPRSFI(:,1) + 0.5_JPRB*(ZPRSFI(:,1)-ZPRSFI(:,2))
  ZPRSHI(:,KNBP) = ZPRSFI(:,KNBP) +0.5_JPRB*(ZPRSFI(:,KNBP)-ZPRSFI(:,KNBP-1))
!       Prepare post-processing on standard levels
  CALL PPINIT(KNBP,1,KNBP,KNBP,JPNPPM,ZPRSHI,ZPRSFI,ZPXP,ZPXPD)
  ZPRESH(NFLEVG)=VP00
  CALL GPHPRE(1,NFLEVG,1,1,YDVAB,ZPRESH,PRESF=ZPRESF)
  ZPRESF(NFLEVG+1)=VP00
  DO JLEV=1,NFLEVG+1
    ZPPP(:,JLEV) = ZPRESF(JLEV)
  ENDDO
  CALL PPFLEV(KNBP,1,KNBP,KNBP,NFLEVG+1,JPNPPM,ZPPP,ZPRSHI,ZPRSFI,&
   & ILEVB,LLBELO,LLBELS,LLBLOW,LLBLES)  
!       Set LLBLOW=.FALSE. so that extrapolation takes place
  LLBLOW(:)=.FALSE.

!       Fill the vertical interpolation matrix
  ZVAL(:,:) = 0.0_JPRB
  DO J1=1,KNBP
    ZVAL(J1,J1) = 1.0_JPRB
  ENDDO
  CALL PPINTP(KNBP,1,KNBP,KNBP,NFLEVG+1,1,JPNPPM,ILEVB,2,&
   & LLBELO,LLBLOW,ZPPP,ZPXP,ZPXPD,ZVAL,ZVALO)  
  IF(KPRT >= 2)THEN
    WRITE(NULOUT,*) ' SUPRECOV: VERTICAL INTERPOL MATRIX*100'
    DO J1=1,KNBP
      WRITE(NULOUT,'(1X,33I3)')&
       & J1,(NINT(ZVALO(J1,J2)*100._JPRB),J2=1,NFLEVG+1)  
    ENDDO
  ENDIF
  DO JLEV=1,NFLEVG+1
    PVI(JLEV,:) = ZVALO(:,JLEV)
  ENDDO
ENDIF

!              3. Find lowest level above 100hPa

ZPRESH(NFLEVG)=VP00
CALL GPHPRE(1,NFLEVG,1,1,YDVAB,ZPRESH,PRESF=ZPRESF)
KLV100=0
DO JLEV=1,NFLEVG
  IF(ZPRESF(JLEV) < 100.E02_JPRB) KLV100=JLEV
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUPRECOV',1,ZHOOK_HANDLE)
END SUBROUTINE SUPRECOV
