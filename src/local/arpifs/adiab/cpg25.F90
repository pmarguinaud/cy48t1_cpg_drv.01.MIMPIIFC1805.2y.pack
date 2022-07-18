SUBROUTINE CPG25(YDGEOMETRY,YDRIP,YDML_DYN,CDCONF,KNUMB,KSTGLO,PB15,PB2,&
 & KIBL,&
 & PUT0,PVT0,PUT0L,PVT0L,PDIVT0,PVORT0,&
 & PSPT0,PSPT0L,PSPT0M,PTRAJ_SLAG)  

!**** *CPG25* - Grid point calculations - 2-D model.

!     Purpose.
!     --------
!           Grid point calculations in dynamics - 2-D MODEL.
!     ****  Either repeats some trajectory calculations for
!     ****  the adjoint of the semi-Lagrangian scheme - or
!     ****  reads them back from additional trajectory storage.

!**   Interface.
!     ----------
!        *CALL* *CPG25(...)

!        Explicit arguments :
!        --------------------

!        CDCONF    : configuration of work                       (I)
!        KNUMB     : number of elements of arrays for which          
!                    computations are performed (MP version)     (I)
!        KSTGLO    : global offset (MP version)                  (I)
!        PB15       : SL buffer for quantities to be interpolated (O)
!        PB2       : Buffer for t+dt quantities                  (O)
!        KIBL      : index into YRGSGEOM instance in YDGEOMETRY  (I)
!        PUT0      : zonal wind time t                           (I/O)
!        PVT0      : meridian wind time t                        (I/O)
!        PUT0L     : zonal derivative of zonal wind time t       (I/O)
!        PVT0L     : zonal derivative of meridian wind time t    (I/O)
!        PDIVT0    : divergence  time t                          (I/O)
!        PVORT0    : vorticity  time t                           (I/O)
!        PSPT0     : equivalent height time t                    (I)
!        PSPT0L    : zonal derivative of PSPT0                   (I/O)
!        PSPT0M    : meridian derivative of PSPT0                (I/O)
!        PTRAJ_SLAG : stored trajectory from NL model            (I)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Calls RDSLTRAJ to get stored trajectory information.
!        Calls LACDYNSHW if SL scheme.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Clive Temperton *ECMWF* (after CPG2)
!        Original : 98-12-22

!     Modifications.
!     --------------
!        Modified 02-11-14 by K. YESSAD  : cleanings.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!        F. Vana  28-Nov-2013 : Redesigned trajectory handling.
!        K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD , ONLY : MODEL_DYNAMICS_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCT0   , ONLY : LTWOTL
USE YOMCT3   , ONLY : NSTEP
USE YOMRIP   , ONLY : TRIP

USE YOMVWRK  , ONLY : NTRSLTYPE
USE YOMTRAJ  , ONLY : TRAJ_SLAG_TYPE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL_DYNAMICS_TYPE),INTENT(INOUT):: YDML_DYN
TYPE(TRIP)          ,INTENT(INOUT) :: YDRIP
CHARACTER(LEN=1)    ,INTENT(IN)    :: CDCONF 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KNUMB 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KSTGLO 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PB15(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB15%NFLDSLB15) 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2) 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KIBL
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PUT0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PVT0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PUT0L(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PVT0L(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PDIVT0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PVORT0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PSPT0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PSPT0L(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)     ,INTENT(INOUT) :: PSPT0M(YDGEOMETRY%YRDIM%NPROMA) 
TYPE(TRAJ_SLAG_TYPE),INTENT(IN)    :: PTRAJ_SLAG

!     ------------------------------------------------------------------

! * time 9:
REAL(KIND=JPRB) :: ZDIVT9  (YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZUT9    (YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZVT9    (YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZSPT9   (YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZSPT9L  (YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZSPT9M  (YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZSPNLT9 (YDGEOMETRY%YRDIM%NPROMA)
! * other:
REAL(KIND=JPRB) :: ZDUM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZOROGL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZOROGM(YDGEOMETRY%YRDIM%NPROMA)

INTEGER(KIND=JPIM) :: IEND, IST

LOGICAL :: LLREAD

REAL(KIND=JPRB) :: ZDT
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "lacdynshw.intfb.h"
#include "rdsltraj2.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CPG25',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, &
 & YDGEM=>YDGEOMETRY%YRGEM, YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDOROG=>YDGEOMETRY%YROROG(KIBL), &
 & YDPTRSLB2=>YDML_DYN%YRPTRSLB2,YDPTRSLB15=>YDML_DYN%YRPTRSLB15)

ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & MSLB1SP05=>YDPTRSLB15%MSLB1SP05, MSLB1SP95=>YDPTRSLB15%MSLB1SP95, &
 & MSLB1U05=>YDPTRSLB15%MSLB1U05, MSLB1U95=>YDPTRSLB15%MSLB1U95, &
 & MSLB1UR05=>YDPTRSLB15%MSLB1UR05, MSLB1V05=>YDPTRSLB15%MSLB1V05, &
 & MSLB1V95=>YDPTRSLB15%MSLB1V95, MSLB1VR05=>YDPTRSLB15%MSLB1VR05, &
 & NFLDSLB15=>YDPTRSLB15%NFLDSLB15, &
 & MSLB2U15=>YDPTRSLB2%MSLB2U15, MSLB2URL5=>YDPTRSLB2%MSLB2URL5, &
 & MSLB2USI5=>YDPTRSLB2%MSLB2USI5, MSLB2V15=>YDPTRSLB2%MSLB2V15, &
 & MSLB2VRL5=>YDPTRSLB2%MSLB2VRL5, MSLB2VSI5=>YDPTRSLB2%MSLB2VSI5, &
 & NFLDSLB2=>YDPTRSLB2%NFLDSLB2, &
 & TDT=>YDRIP%TDT)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS.
!              ----------------------------

IF    (CDCONF == 'A')THEN
  ZDT=TDT
ELSEIF(CDCONF == 'B')THEN
  ZDT=0.0_JPRB
ELSE
  CALL ABOR1(' CPG25 - UNKNOWN CONFIGURATION ')
ENDIF

IST =1
IEND=KNUMB

!     ------------------------------------------------------------------

!*       2.    INTERFACE TO GLOBAL ARRAYS.
!              ---------------------------

!       * Preliminary fill buffers with zero.

PB15(IST:IEND,:)=0.0_JPRB
PB2(IST:IEND,:)=0.0_JPRB

IF (NTRSLTYPE <= 1) THEN

!*       2.1   CASE NTRSLTYPE <=1.

!*         2.1.1   GET t-dt VARIABLES.

  LLREAD= (NSTEP > 0)

  IF(LLREAD) THEN
    CALL RDSLTRAJ2(YDDIM,IST,IEND,PTRAJ_SLAG,ZUT9,ZVT9,ZSPNLT9,&
     & ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,ZDUM,&
     & ZDUM,ZDUM)  
  ENDIF

!*         2.1.3   MAP FACTOR

  PUT0  (IST:IEND)=PUT0  (IST:IEND)* YDGSGEOM%GM(IST:IEND)
  PVT0  (IST:IEND)=PVT0  (IST:IEND)* YDGSGEOM%GM(IST:IEND)
  PUT0L (IST:IEND)=PUT0L (IST:IEND)*(YDGSGEOM%GM(IST:IEND)**2)
  PVT0L (IST:IEND)=PVT0L (IST:IEND)*(YDGSGEOM%GM(IST:IEND)**2)
  PDIVT0(IST:IEND)=PDIVT0(IST:IEND)*(YDGSGEOM%GM(IST:IEND)**2)
  PVORT0(IST:IEND)=PVORT0(IST:IEND)*(YDGSGEOM%GM(IST:IEND)**2)
  PSPT0L(IST:IEND)=PSPT0L(IST:IEND)* YDGSGEOM%GM(IST:IEND)
  PSPT0M(IST:IEND)=PSPT0M(IST:IEND)* YDGSGEOM%GM(IST:IEND)
  IF (LLREAD.AND..NOT.YDML_DYN%YRDYNA%LELTRA) THEN
    ZUT9  (IST:IEND)=ZUT9  (IST:IEND)* YDGSGEOM%GM(IST:IEND)
    ZVT9  (IST:IEND)=ZVT9  (IST:IEND)* YDGSGEOM%GM(IST:IEND)
  ENDIF
  IF (.NOT.LTWOTL) THEN
    ZDIVT9(IST:IEND)=ZDIVT9(IST:IEND)*(YDGSGEOM%GM(IST:IEND)**2)
    ZSPT9L(IST:IEND)=ZSPT9L(IST:IEND)* YDGSGEOM%GM(IST:IEND)
    ZSPT9M(IST:IEND)=ZSPT9M(IST:IEND)* YDGSGEOM%GM(IST:IEND)
  ENDIF
  ZOROGL(IST:IEND)=YDOROG%OROGL(IST:IEND)* YDGSGEOM%GM(IST:IEND)
  ZOROGM(IST:IEND)=YDOROG%OROGM(IST:IEND)* YDGSGEOM%GM(IST:IEND)

!*         2.1.4   SHALLOW WATER MODEL: SEMI-LAGRANGIAN SCHEME.

!     * Fill buffers.

  ZDUM(1:NPROMA)=0.0_JPRB
  CALL LACDYNSHW(YDGEM,YDML_DYN%YRDYN,CDCONF,NPROMA,IST,IEND,ZDT,&
   & ZUT9,ZVT9,ZSPT9,ZDIVT9,ZSPT9L,ZSPT9M,&
   & PUT0,PVT0,PSPT0,PDIVT0,PSPT0L,PSPT0M,&
   & YDGSGEOM,YDOROG%OROG,ZOROGL,ZOROGM,&
   & ZSPNLT9,&
   & PB15(1,MSLB1U95),PB15(1,MSLB1V95),&
   & PB15(1,MSLB1SP95),PB15(1,MSLB1U05),&
   & PB15(1,MSLB1V05),PB15(1,MSLB1SP05),&
   & PB2(1,MSLB2U15),PB2(1,MSLB2V15),ZDUM,&
   & PB15(1,MSLB1UR05),PB15(1,MSLB1VR05),&
   & PB2(1,MSLB2USI5),PB2(1,MSLB2VSI5),&
   & PB2(1,MSLB2URL5),PB2(1,MSLB2VRL5))

ELSEIF (NTRSLTYPE == 2) THEN

!*       2.2   CASE NTRSLTYPE =2.

  CALL RDSLTRAJ2(YDDIM,IST,IEND,PTRAJ_SLAG,ZDUM,ZDUM,ZDUM,&
   & PB15(1,MSLB1U95),PB15(1,MSLB1V95),&
   & PB15(1,MSLB1SP95),PB15(1,MSLB1U05),&
   & PB15(1,MSLB1V05),PB15(1,MSLB1SP05),&
   & PB15(1,MSLB1UR05),PB15(1,MSLB1VR05),&
   & PB2(1,MSLB2URL5),PB2(1,MSLB2VRL5),&
   & PB2(1,MSLB2USI5),PB2(1,MSLB2VSI5),&
   & PB2(1,MSLB2U15),PB2(1,MSLB2V15) )  

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPG25',1,ZHOOK_HANDLE)
END SUBROUTINE CPG25
