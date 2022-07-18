SUBROUTINE CTVTOTTL(YDGEOMETRY,YDDYN,KSTART,KPROF,PT,PQ,PT5,PQ5)

!**** *CTVTOTTL* - Computes T from 'TV' (TV defined as RT/Rd)

!     Purpose.
!     --------
!           COMPUTES T FROM 'TV'    (T/L version)
!            (and T5 from TV5)

!**   Interface.
!     ----------
!        *CALL* *CTVTOTTL(...)

!        Explicit arguments :
!        --------------------
!                              PT(NPROMA,NFLEVG)    - output: T
!                                                  - input: TV
!                              PQ(NPROMA,NFLEVG)    - Q

!        Implicit arguments :
!        --------------------

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
!        Clive Temperton  *ECMWF*

!     Modifications.
!     --------------
!        Original : 01-07-24
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Jidane  08-04-2006 : Reintro of R dep. on q variables under key
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMCT0 , ONLY : L_OOPS
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCST   , ONLY : RD, RV
USE YOMDYN   , ONLY : TDYN
USE YOMDYNA  , ONLY : YRDYNA

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(INOUT) :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PT5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZR   ,ZR5

INTEGER(KIND=JPIM) :: IEND, IST, JLEV, JROF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CTVTOTTL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, NFLEVG=>YDDIMV%NFLEVG)
!     ------------------------------------------------------------------

!*       1.    COMPUTES T FROM 'TV'
!              --------------------

IST =KSTART
IEND=KPROF

DO JLEV = 1, NFLEVG
  DO JROF = IST, IEND
    ZR5=RD+(RV-RD)*PQ5(JROF,JLEV)
    IF(.NOT.L_OOPS)PT5(JROF,JLEV)=RD*PT5(JROF,JLEV)/ZR5
    PT(JROF,JLEV)=RD*PT(JROF,JLEV)/ZR5
    IF (.NOT. YRDYNA%LDRY_ECMWF) THEN
      ZR=(RV-RD)*PQ(JROF,JLEV)
      PT(JROF,JLEV)=PT(JROF,JLEV)-ZR*PT5(JROF,JLEV)/ZR5
    ENDIF
  ENDDO
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CTVTOTTL',1,ZHOOK_HANDLE)
END SUBROUTINE CTVTOTTL
