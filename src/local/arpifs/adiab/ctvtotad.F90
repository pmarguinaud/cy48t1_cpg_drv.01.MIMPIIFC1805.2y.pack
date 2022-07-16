SUBROUTINE CTVTOTAD(YDGEOMETRY,YDDYN,KSTART,KPROF,PT,PQ,PT5,PQ5)

!**** *CTVTOTAD* - Computes T from 'TV' (TV defined as RT/Rd)

!     Purpose.
!     --------
!           COMPUTES T FROM 'TV'  (adjoint version)

!**   Interface.
!     ----------
!        *CALL* *CTVTOTAD(...)

!        Explicit arguments :
!        --------------------
!                              PT(NPROMA,NFLEVG)    - input: T
!                                                  - output: TV
!                              PQ(NPROMA,NFLEVG)    - Q
!                              PT5, PQ5 - same but for trajectory
!        (N.B. assumes PT5 has ** ALREADY ** been converted to T)

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
!        Original : 01-01-29
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Jidane 09-04-2006 : R with q variables under key
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCST   , ONLY : RD, RV
USE YOMDYN   , ONLY : TDYN
USE YOMDYNA  , ONLY : LDRY_ECMWF

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(INOUT) :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQ(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZR5, ZR0

INTEGER(KIND=JPIM) :: IEND, IST, JLEV, JROF
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CTVTOTAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NPROMA=>YDDIM%NPROMA,NFLEVG=>YDDIMV%NFLEVG)
!     ------------------------------------------------------------------

!*       1.    COMPUTES T FROM 'TV'
!              --------------------

IST =KSTART
IEND=KPROF

DO JLEV = 1, NFLEVG
  DO JROF = IST, IEND
    ZR5=RD+(RV-RD)*PQ5(JROF,JLEV)
    IF (.NOT. LDRY_ECMWF) THEN
      ZR0=-PT5(JROF,JLEV)*PT(JROF,JLEV)/ZR5
      PQ(JROF,JLEV)=PQ(JROF,JLEV)+(RV-RD)*ZR0
    ENDIF
    PT(JROF,JLEV)=RD*PT(JROF,JLEV)/ZR5
  ENDDO
ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CTVTOTAD',1,ZHOOK_HANDLE)
END SUBROUTINE CTVTOTAD

