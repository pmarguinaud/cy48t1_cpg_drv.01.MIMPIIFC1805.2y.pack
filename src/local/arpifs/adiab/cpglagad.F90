SUBROUTINE CPGLAGAD(YDGEOMETRY,YDGMV,YDML_GCONF,YDDYN,YDPTRSLB2,KGPCOMP,&
 & PB2,PGMV,PGMVS,PGFL,PGMVT1,PGMVT1S,PGFLT1)  

!**** *CPGLAGAD* - Grid point calculations lagged part.adjoint.

!     Purpose.
!     --------
!           Grid point calculations lagged part.AD

!**   Interface.
!     ----------
!        *CALL* *CPGLAGAD(...)*

!        Explicit arguments :
!        --------------------

!        KGPCOMP   : number of elements of arrays for which computations
!                    are performed (used in message passing version)
!        PB2       : "SLB2" buffer.
!        PGMV      : "t-dt" and "t" GMV upper air variables.
!        PGMVS     : "t-dt" and "t" GMV surface variables.
!        PGFL      : "t-dt" and "t" GFL variables.
!        PGMVT1    : "t+dt" GMV upper air variables.
!        PGMVT1S   : "t+dt" GMV surface variables.
!        PGFLT1    : "t+dt" GFL variables.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Called by GP_MODEL_AD.

!     Reference.
!     ----------
!        ARPEGE documentation vol 2 ch 1 and vol 3 ch 6

!     Author.
!     -------
!      Jean-Noel Thepaut  *ECMWF*
!      Original : 92-04-08

! Modifications
! -------------
!   C.Temperton: 01-01-26 add the case LSPRT=.TRUE.
!   R. El Khatib : 01-08-07 Pruning options
!   03-2002 J.Vivoda PC schemes for NH dynamics (LPC_XXXX keys)
!   01-Oct-2003 M. Hamrud  CY28 Cleaning
!   09-Jun-2004 J. Masek   NH cleaning (LPC_NOTR, LFULLIMP)
!   01-Jul-2004 K. Yessad  Make clearer the tests for PC scheme.
!   09-Apr-2006 M. Jidane  R dep. on q variables under key
!   K. Yessad (Dec 2008): cleanings
!   G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!   F. Vana    19-03-2012: fix of LSPRT option (correct increments)
!   K. Yessad (Nov 2012): simplify testings.
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
! End Modifications
!------------------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV   , ONLY : TGMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE PTRSLB2  , ONLY : TPTRSLB2
USE YOMCST   , ONLY : RD, RV
USE YOMCT0   , ONLY : LNHDYN, LTWOTL, LSPRT
USE YOMCT3   , ONLY : NSTEP
USE YOMDYN   , ONLY : TDYN
USE YOMDYNA   ,ONLY : LGRADSP, LDRY_ECMWF

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)               ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)                   ,INTENT(INOUT) :: YDGMV
TYPE(TDYN)                   ,INTENT(INOUT) :: YDDYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT) :: YDML_GCONF
TYPE(TPTRSLB2)               ,INTENT(INOUT) :: YDPTRSLB2
INTEGER(KIND=JPIM)           ,INTENT(IN)    :: KGPCOMP
REAL(KIND=JPRB)              ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDPTRSLB2%NFLDSLB2,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)              ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV, &
 &                                                  YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)              ,INTENT(INOUT) :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)              ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM, &
 &                                                  YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)              ,INTENT(INOUT) :: PGMVT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%YT1%NDIM,&
 &                                                    YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)              ,INTENT(IN)    :: PGMVT1S(YDGEOMETRY%YRDIM%NPROMA,YDGMV%YT1%NDIMS,YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)              ,INTENT(INOUT) :: PGFLT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM1, &
 &                                                    YDGEOMETRY%YRDIM%NGPBLKS)
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: IEND, ILEV, IST, JLEV, JROF, JGFL, JKGLO, IBL

LOGICAL :: LLLSTEP, LLSTR
REAL(KIND=JPRB) :: ZEPS
REAL(KIND=JPRB) :: ZR5   ,ZR

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gptf1ad.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CPGLAGAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
 & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDGSGEOM=>YDGEOMETRY%YRGSGEOM, &
 & YGFL=>YDML_GCONF%YGFL)

ASSOCIATE(NDIM=>YGFL%NDIM, NDIM1=>YGFL%NDIM1, NUMFLDS=>YGFL%NUMFLDS, &
 & YCOMP=>YGFL%YCOMP, YQ=>YGFL%YQ, &
 & NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
 & RSTRET=>YDGEM%RSTRET, &
 & NDIMGMV=>YDGMV%NDIMGMV, NDIMGMVS=>YDGMV%NDIMGMVS, YT0=>YDGMV%YT0, &
 & YT1=>YDGMV%YT1, YT9=>YDGMV%YT9, &
 & MSLB2PDSI=>YDPTRSLB2%MSLB2PDSI, MSLB2Q15=>YDPTRSLB2%MSLB2Q15, &
 & MSLB2SPSI=>YDPTRSLB2%MSLB2SPSI, MSLB2T15=>YDPTRSLB2%MSLB2T15, &
 & MSLB2TSI=>YDPTRSLB2%MSLB2TSI, MSLB2USI=>YDPTRSLB2%MSLB2USI, &
 & MSLB2VDSI=>YDPTRSLB2%MSLB2VDSI, MSLB2VSI=>YDPTRSLB2%MSLB2VSI, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2)
!     ------------------------------------------------------------------

!*       1.    INITIAL COMPUTATIONS
!              --------------------

ZEPS=100.0_JPRB*TINY(1.0_JPRB)
LLSTR=(ABS(RSTRET-1.0_JPRB)>ZEPS)

IF(NSTEP < YDML_GCONF%YRRIP%NSTOP-1) THEN
  LLLSTEP=.FALSE.
ELSE
  LLLSTEP=.TRUE.
ENDIF

IST=1

!     ------------------------------------------------------------------

!*       2.    LOOP ON JKGLO
!              -------------

CALL GSTATS(1012,0)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,IEND,IBL,JLEV,JROF,ILEV,ZR,ZR5,JGFL)
DO JKGLO=1,KGPCOMP,NPROMA

  IEND=MIN(NPROMA,KGPCOMP-JKGLO+1)
  IBL=(JKGLO-1)/NPROMA+1

  !*    2.6  DIVIDE U AND V BY MAP FACTOR
  !          ----------------------------

  IF(LLSTR) THEN
    DO JLEV=1,NFLEVG
      DO JROF=IST,IEND
        PGMVT1(JROF,JLEV,YT1%MU,IBL) = PGMVT1(JROF,JLEV,YT1%MU,IBL)/YDGSGEOM(IBL)%GM(JROF)
        PGMVT1(JROF,JLEV,YT1%MV,IBL) = PGMVT1(JROF,JLEV,YT1%MV,IBL)/YDGSGEOM(IBL)%GM(JROF)
      ENDDO
    ENDDO
  ENDIF

  !*    2.5  TRANSFER FROM BUFFER TO T+DT ARRAYS
  !          -----------------------------------

  IF (LSPRT) THEN

    DO JLEV=1,NFLEVG
      ILEV=JLEV-1
      DO JROF=IST,IEND
        ZR5=RD+(RV-RD)*PB2(JROF,MSLB2Q15+ILEV,IBL)
        IF (.NOT. LDRY_ECMWF) THEN
          ZR=PB2(JROF,MSLB2T15+ILEV,IBL)*PGMVT1(JROF,JLEV,YT1%MT,IBL)/RD
          PGFLT1(JROF,JLEV,YQ%MP1,IBL)=PGFLT1(JROF,JLEV,YQ%MP1,IBL)+(RV-RD)*ZR
        ENDIF
        PGMVT1(JROF,JLEV,YT1%MT,IBL)=ZR5*PGMVT1(JROF,JLEV,YT1%MT,IBL)/RD
      ENDDO
    ENDDO

  ENDIF

  DO JLEV=1,NFLEVG
    ILEV = JLEV-1
    DO JROF=IST,IEND
      PB2(JROF,MSLB2USI+ILEV,IBL)=PB2(JROF,MSLB2USI+ILEV,IBL)+PGMVT1(JROF,JLEV,YT1%MU,IBL)
      PB2(JROF,MSLB2VSI+ILEV,IBL)=PB2(JROF,MSLB2VSI+ILEV,IBL)+PGMVT1(JROF,JLEV,YT1%MV,IBL)
      PB2(JROF,MSLB2TSI+ILEV,IBL)=PB2(JROF,MSLB2TSI+ILEV,IBL)+PGMVT1(JROF,JLEV,YT1%MT,IBL)
    ENDDO
    IF(LNHDYN) THEN
      DO JROF=IST,IEND
        PB2(JROF,MSLB2PDSI+ILEV,IBL)= PB2(JROF,MSLB2PDSI+ILEV,IBL)+&
         & PGMVT1(JROF,JLEV,YT1%MSPD,IBL)
        PB2(JROF,MSLB2VDSI+ILEV,IBL)= PB2(JROF,MSLB2VDSI+ILEV,IBL)+&
         & PGMVT1(JROF,JLEV,YT1%MSVD,IBL)
      ENDDO
    ENDIF
  ENDDO
  DO JROF=IST,IEND
    PB2(JROF,MSLB2SPSI,IBL) = PB2(JROF,MSLB2SPSI,IBL)+PGMVT1S(JROF,YT1%MSP,IBL)
  ENDDO

  !*    2.3   TIME FILTER (PART 1).
  !           ---------------------

  IF(NCURRENT_ITER == 0)THEN
    IF(LTWOTL) THEN
      IF(LLLSTEP) THEN
        PGMV(:,:,YT0%NDIM+1:,IBL) = 0.0_JPRB
        PGMVS(:,YT0%NDIMS+1:,IBL) = 0.0_JPRB
        DO JGFL=1,NUMFLDS
          IF(YCOMP(JGFL)%LT9) THEN
            PGFL(:,:,YCOMP(JGFL)%MP9,IBL) = 0.0_JPRB
          ENDIF
        ENDDO
      ENDIF
      DO JLEV=1,NFLEVG
        DO JROF=IST,IEND
          PGMV(JROF,JLEV,YT0%MU,IBL)=PGMV(JROF,JLEV,YT0%MU,IBL)+PGMV(JROF,JLEV,YT9%MU,IBL)
          PGMV(JROF,JLEV,YT9%MU,IBL)=0.0_JPRB
          PGMV(JROF,JLEV,YT0%MV,IBL)=PGMV(JROF,JLEV,YT0%MV,IBL)+PGMV(JROF,JLEV,YT9%MV,IBL)
          PGMV(JROF,JLEV,YT9%MV,IBL)=0.0_JPRB
        ENDDO
      ENDDO
    ELSE
      CALL GPTF1AD(YDGEOMETRY,YDGMV,YDML_GCONF,YDDYN,IST,IEND,LLLSTEP,PGMV(1,1,1,IBL),PGMVS(1,1,IBL),PGFL(1,1,1,IBL))
    ENDIF
  ENDIF

  !*    2.2   UPDATE U AND V FOR LGPSTRESS OPTION
  !           -----------------------------------

  IF((LTWOTL.AND.NCURRENT_ITER == 0).OR.(NCURRENT_ITER>0))THEN
    IF(LGRADSP)THEN
      DO JLEV=1,NFLEVG
        DO JROF=IST,IEND
          PGMV(JROF,JLEV,YT0%MSGRTL,IBL)=PGMV(JROF,JLEV,YT0%MSGRTL,IBL)+PGMV(JROF,JLEV,YT9%MSGRTL,IBL)
          PGMV(JROF,JLEV,YT9%MSGRTL,IBL)=0.0_JPRB
          PGMV(JROF,JLEV,YT0%MSGRTM,IBL)=PGMV(JROF,JLEV,YT0%MSGRTM,IBL)+PGMV(JROF,JLEV,YT9%MSGRTM,IBL)
          PGMV(JROF,JLEV,YT9%MSGRTM,IBL)=0.0_JPRB
        ENDDO
      ENDDO
    ENDIF
  ENDIF

ENDDO ! JKGLO
!$OMP END PARALLEL DO
CALL GSTATS(1012,1)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPGLAGAD',1,ZHOOK_HANDLE)
END SUBROUTINE CPGLAGAD
