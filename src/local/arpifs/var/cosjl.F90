SUBROUTINE COSJL(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL)

!**** *COSJL*  - LINEAR COST FUNCTION IN SPECTRAL SPACE

!     Purpose.
!     --------
!     EVALUATE THE COST-FUNCTION AND THE GRADIENT

!**   Interface.
!     ----------
!        *CALL* *COSJL

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        None

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
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15 under the name COSTRA

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib : 01-08-07 Pruning options
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      A. Geer      27-Jul-2015  VarBC is now an object passed by argument, for OOPS
!      O. Marsden   August 2016  Removed use of SPA3

!     ------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE FIELDS_MOD , ONLY : FIELDS
USE MTRAJ_MOD  , ONLY : MTRAJ
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMCOSJR , ONLY : FSPALT ,FSPSUR
USE YOMSP5   , ONLY : SPA5
USE YOMLCZ   , ONLY : LOCNORM
USE YOMVERT  , ONLY : VP00
USE YOMMP0   , ONLY : MYSETW
USE MPL_MODULE
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)
USE VARBC_CLASS, ONLY : CLASS_VARBC

IMPLICIT NONE

!     ----------------------------------------

TYPE(GEOMETRY)   ,INTENT(INOUT) :: YDGEOMETRY  !! INOUT needed for STEPOAD call
TYPE(FIELDS)     ,INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ)      ,INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL)      ,INTENT(INOUT) :: YDMODEL

CHARACTER (LEN = 9) :: CLCONF

REAL(KIND=JPRB), DIMENSION (YDGEOMETRY%YRGEM%NGPTOT) :: ZGDSP,ZGDSP2
REAL(KIND=JPRB), DIMENSION (6) :: ZZ

INTEGER(KIND=JPIM) :: ILEV, JLEV, JMLOC, JVAL, JVAR

REAL(KIND=JPRB) :: ZO3
REAL(KIND=JPRB) :: ZSGLOB,ZSDOM,ZVERT,ZRAP
REAL(KIND=JPRB) :: ZVOR,ZDIV,ZT,ZQ,ZSP
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "cossmq.intfb.h"
#include "reespe.intfb.h"
#include "speree.intfb.h"
#include "spnorm.intfb.h"
#include "stepo.intfb.h"
#include "stepoad.intfb.h"
#include "stepotl.intfb.h"

!     ----------------------------------------

!     1. Starts COSJL
!     ---------------

IF (LHOOK) CALL DR_HOOK('COSJL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,   YDDIMV=>YDGEOMETRY%YRDIMV, YDGEM=>YDGEOMETRY%YRGEM,  YDMP=>YDGEOMETRY%YRMP,  &
& YDLAP=>YDGEOMETRY%YRLAP,    YDSTA=>YDGEOMETRY%YRSTA, YDDYN=>YDMODEL%YRML_DYN%YRDYN,YGFL=>YDMODEL%YRML_GCONF%YGFL&
& )

ASSOCIATE(YO3=>YGFL%YO3,   NSMAX=>YDDIM%NSMAX, NSPEC2=>YDDIM%NSPEC2, NUMP=>YDDIM%NUMP,   NFLEVL=>YDDIMV%NFLEVL,   &
& SIDELP=>YDDYN%SIDELP, SIPR=>YDDYN%SIPR,   NGPTOT=>YDGEM%NGPTOT,   NASM0=>YDLAP%NASM0,   MYLEVS=>YDMP%MYLEVS,    &
& NPROCM=>YDMP%NPROCM, NPSP=>YDMP%NPSP,   STTEM=>YDSTA%STTEM)
ZSGLOB=0.0_JPRB
ZSDOM=0.0_JPRB
ZVERT=0.0_JPRB
ZRAP=0.0_JPRB

ZGDSP(:)=0.0_JPRB
ZGDSP2(:)=0.0_JPRB

CLCONF(1:9)='0GB0L0AA0'

CALL COSSMQ(YDGEOMETRY,YDDYN,ZSGLOB,ZSDOM,ZVERT)
ZRAP=ZSGLOB*SIPR/(ZSDOM*ZVERT)
!     --------------------------------------------------------------

!     2. COMPUTATION OF THE COST FUNCTION AND ITS GRADIENT.
!     -----------------------------------------------------

!     2.1 Put trajectory
!     ------------------

YDFIELDS%YRSPEC%VOR(:,:) =SPA5%VOR(:,:)  !! was SPA3
YDFIELDS%YRSPEC%DIV(:,:) =SPA5%DIV(:,:)
YDFIELDS%YRSPEC%T(:,:)   =SPA5%T(:,:)
YDFIELDS%YRSPEC%Q(:,:)   =SPA5%Q(:,:)
IF (YO3%LSP) THEN
  YDFIELDS%YRSPEC%O3(:,:)  =SPA5%O3(:,:)
ENDIF
IF (NPSP == 1) THEN
  YDFIELDS%YRSPEC%SP(:)=SPA5%SP(:)
ENDIF

CALL SPNORM(YDGEOMETRY,YDMODEL%YRML_GCONF,YDFIELDS%YRSPEC)

!     2.2 Subtract reference values (standard atmosphere) to traj
!     -----------------------------------------------------------

IF (MYSETW == NPROCM(0)) THEN
  YDFIELDS%YRSPEC%T(:,NASM0(0))=YDFIELDS%YRSPEC%T(:,NASM0(0))-STTEM(:) !! was SPA3
ENDIF

! Log of surface pressure
CALL SPEREE(YDGEOMETRY,1,1,YDFIELDS%YRSPEC%SP,ZGDSP)  !! was SPA3
DO JVAL=1,NGPTOT
  ZGDSP(JVAL)=EXP(ZGDSP(JVAL))-VP00
ENDDO
CALL REESPE(YDGEOMETRY,1,1,YDFIELDS%YRSPEC%SP,ZGDSP)  !! was SPA3

!     2.3 Mask the difference
!     -----------------------

IF(LOCNORM) THEN
  CALL STEPO(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,'0K0000000')
  CALL STEPOTL(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,CLCONF)
ENDIF

CALL SPNORM(YDGEOMETRY,YDMODEL%YRML_GCONF,YDFIELDS%YRSPEC)

!     2.4 Cost function & Gradient
!     ----------------------------

DO JLEV=1,NFLEVL
  ILEV=MYLEVS(JLEV)

!      2.4.1 3D cost-function & Gradients
!      ----------------------------------
  ZVOR=SIDELP(ILEV)/VP00/1.0_JPRB
  ZDIV=SIDELP(ILEV)/VP00/1.0_JPRB
  ZT=SIDELP(ILEV)/VP00/1.0_JPRB
  ZQ=SIDELP(ILEV)/VP00/1.E-4_JPRB
  ZO3=0.0_JPRB

  IF(LOCNORM) THEN
    ZVOR=ZVOR*ZRAP
    ZDIV=ZDIV*ZRAP
    ZT=ZT*ZRAP
    ZQ=ZQ*ZRAP
    ZO3=ZO3*ZRAP
  ENDIF

  IF (MYSETW == NPROCM(0)) THEN
    FSPALT(ILEV,1,1)=YDFIELDS%YRSPEC%VOR(JLEV,NASM0(0))*ZVOR  !! was SPA3
    FSPALT(ILEV,2,1)=YDFIELDS%YRSPEC%DIV(JLEV,NASM0(0))*ZDIV
    FSPALT(ILEV,3,1)=YDFIELDS%YRSPEC%T (JLEV,NASM0(0))*ZT
    FSPALT(ILEV,4,1)=YDFIELDS%YRSPEC%Q (JLEV,NASM0(0))*ZQ
    IF (YO3%LSP) THEN
      FSPALT(ILEV,5,1)=YDFIELDS%YRSPEC%O3(JLEV,NASM0(0))*ZO3
    ELSE
      FSPALT(ILEV,5,1)=0.0_JPRB
    ENDIF
  ELSE
    DO JMLOC=1,NUMP
      FSPALT(ILEV,1:5,JMLOC)=0.0_JPRB
    ENDDO
  ENDIF
!*3D Gradients*!        
  YDFIELDS%YRSPEC%VOR(JLEV,:)=0.0_JPRB  !! was SPA3
  YDFIELDS%YRSPEC%DIV(JLEV,:)=0.0_JPRB
  YDFIELDS%YRSPEC%T(JLEV,:)=0.0_JPRB
  YDFIELDS%YRSPEC%Q(JLEV,:)=0.0_JPRB
  YDFIELDS%YRSPEC%O3(JLEV,:)=0.0_JPRB
  IF(MYSETW == NPROCM(0)) THEN
    YDFIELDS%YRSPEC%VOR(JLEV,NASM0(0))=ZVOR
    YDFIELDS%YRSPEC%DIV(JLEV,NASM0(0))=ZDIV
    YDFIELDS%YRSPEC%T(JLEV,NASM0(0))=ZT
    YDFIELDS%YRSPEC%Q(JLEV,NASM0(0))=ZQ
    IF (YO3%LSP) THEN
      YDFIELDS%YRSPEC%O3(JLEV,NASM0(0))=ZO3
    ENDIF
  ENDIF

ENDDO

!      2.4.2 Cost and gradient for Surface pressure
!      --------------------------------------------

IF (NPSP == 1) THEN
  ZSP=1/100._JPRB
  IF(LOCNORM) THEN
    ZSP=ZSP*ZRAP
  ENDIF

!* 2DCost function *!

  IF (MYSETW == NPROCM(0)) THEN
    FSPSUR(1,1)=YDFIELDS%YRSPEC%SP(NASM0(0))*ZSP  !! was SPA3
  ELSE
    DO JMLOC=1,NUMP
      FSPSUR(1,JMLOC)=0.0_JPRB
    ENDDO
  ENDIF
!* 2DGradients *!
  YDFIELDS%YRSPEC%SP(:)=0.0_JPRB  !! was SPA3
  ZGDSP(:)=0.0_JPRB
  ZGDSP2(:)=0.0_JPRB
!* Gradient w.r.t. Log(Ps) directly in gridpoint space *!
  CALL SPEREE(YDGEOMETRY,1,1,SPA5%SP,ZGDSP)
  DO JVAL=1,NGPTOT
    ZGDSP2(JVAL)=ZSP*EXP(ZGDSP(JVAL))
  ENDDO
  CALL REESPE(YDGEOMETRY,1,1,YDFIELDS%YRSPEC%SP,ZGDSP2)  !! was SPA3

  IF(MYSETW == NPROCM(0)) THEN
    YDFIELDS%YRSPEC%SP(2*(NSMAX+1):NSPEC2)=YDFIELDS%YRSPEC%SP(2*(NSMAX+1):NSPEC2)*2.0_JPRB  !! was SPA3
  ELSE
    YDFIELDS%YRSPEC%SP(:)=YDFIELDS%YRSPEC%SP(:)*2.0_JPRB
  ENDIF

ENDIF
!     2.5 Adjoint Mask
!     ----------------

IF(LOCNORM) THEN
  WRITE(NULOUT,*) 'GRADIENT NORMS BEFORE LOCAL MASKING:'
  CALL SPNORM(YDGEOMETRY,YDMODEL%YRML_GCONF,YDFIELDS%YRSPEC)
  CALL STEPOAD(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,CLCONF)
  WRITE(NULOUT,*) 'GRADIENT NORMS AFTER LOCAL MASKING:'
  CALL SPNORM(YDGEOMETRY,YDMODEL%YRML_GCONF,YDFIELDS%YRSPEC)
ELSE
  WRITE(NULOUT,*) 'NORMS OF THE ADDED GRADIENT:'
ENDIF

!     ------------------------------------------------------------------

!     3. Writes diagnostics
!     ---------------------

!    Values are communicated here to achieve global values

ZZ(:) = 0.0_JPRB
WRITE(NULOUT,*) 'Contributions to linear cost-function :'
DO JLEV=1,NFLEVL
  ZZ(1:5) = 0.0_JPRB
  DO JMLOC=1,NUMP
    DO JVAR=1,4
      ZZ(JVAR)=ZZ(JVAR)+FSPALT(JLEV,JVAR,JMLOC)
    ENDDO
    IF (YO3%LSP) THEN
      ZZ(5)=ZZ(5)+FSPALT(JLEV,5,JMLOC)
    ENDIF
  ENDDO
  WRITE(NULOUT,311) JLEV,ZZ(1:5)
  311 FORMAT('LEV=',I3,' Cost U=',G13.3,' V=',G13.3,' T=',G13.3,&
   & ' Q=',G13.3,'O3=',G13.3)  
ENDDO

IF (NPSP == 1) THEN
  ZZ(6)=0.0_JPRB
  DO JMLOC=1,NUMP
    ZZ(6)=ZZ(6)+FSPSUR(1,JMLOC)
  ENDDO
  WRITE(NULOUT,312) ZZ(6)
  312 FORMAT(' Cost Ps=',G14.3)
ENDIF

CALL MPL_ALLREDUCE(ZZ,'SUM',LDREPROD=.FALSE.,CDSTRING='COSJL:')

WRITE(NULOUT,*) 'Total  U=',ZZ(1),' V=',ZZ(2),' T=',&
 & ZZ(3),' Q=',ZZ(4),' O3=',ZZ(5),' Ps=',ZZ(6)  

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('COSJL',1,ZHOOK_HANDLE)
END SUBROUTINE COSJL
