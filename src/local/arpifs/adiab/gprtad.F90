SUBROUTINE GPRTAD(LDSPRT,KPROMA,KSTART,KEND,KLEV,PRD,PRV,&
 & PR5,PT5,PT5L,PT5M,PQ5L,PQ5M,&
 & PR,PT,PTL,PTM,PQL,PQM,PRT,PRTL,PRTM)  

!*    *GPRTAD* 

!     Purpose
!     -------  To calculate RT and its derivates (adjoint code)

!     Interface
!     ---------

!     Explicit arguments
!     ------------------
!     Input:
!    -------
!              LDSPRT   : .TRUE. if PTL and PTM already contain
!                         the derivatives of 'TV')
!              KPROMA   : Horizontal dimension
!              KSTART   : Start index
!              KEND     : End index
!              KLEV     : number of levels
!              PRD      : Rd
!              PRV      : Rv
!              PR       : R
!              PT       : T
!              PTL ,PTM : Horizontal derivatives of T
!              PQL, PQM : Horizontal derivatives of q
!     Trajectory:
!    ------------
!              PR5       : R
!              PT5       : T
!              PT5L ,PT5M : Horizontal derivatives of T
!              PQ5L, PQ5M : Horizontal derivatives of q
!     Output:
!    --------
!              PRT       : RT
!              PRTL,PRTM : Horizontal derivatives of RT

!     Author
!     ------
!           C. Fischer *METEO-FRANCE*
!     Modifications
!     -------------
!           Original: 02/06/28
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M. Jidane 09-04-2006  R with q variables under key
!     K. Yessad (Dec 2008): remove dummy CDLOCK
!----------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
USE YOMDYNA  , ONLY : LDRY_ECMWF

!----------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
LOGICAL           ,INTENT(IN)    :: LDSPRT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRD 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PR5(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT5(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT5L(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT5M(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ5L(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ5M(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PR(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PT(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTL(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTM(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQL(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQM(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRT(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRTL(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRTM(KPROMA,KLEV) 

!----------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JLH
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPRTAD',0,ZHOOK_HANDLE)

!----------------------------------------------------------
!---------------------------------------------------------------
!*       1. Compute RT and its derivatives

DO JLEV=1,KLEV
  DO JLH=KSTART,KEND
    IF (LDSPRT) THEN
      PTL(JLH,JLEV) =PTL(JLH,JLEV) +PRD*PRTL(JLH,JLEV)
      PTM(JLH,JLEV) =PTM(JLH,JLEV) +PRD*PRTM(JLH,JLEV)
      PRTL(JLH,JLEV)=0.0_JPRB
      PRTM(JLH,JLEV)=0.0_JPRB
    ELSE
      PT(JLH,JLEV)  =PT(JLH,JLEV) +&
       & (PRV-PRD)*PQ5M(JLH,JLEV)*PRTM(JLH,JLEV)  
      PTM(JLH,JLEV) =PTM(JLH,JLEV) +&
       & PR5(JLH,JLEV)*PRTM(JLH,JLEV)  
      IF (.NOT. LDRY_ECMWF) THEN
        PQM(JLH,JLEV) =PQM(JLH,JLEV) +&
         & (PRV-PRD)*PT5(JLH,JLEV)*PRTM(JLH,JLEV)  
        PR(JLH,JLEV)  =PR(JLH,JLEV) +&
         & PT5M(JLH,JLEV)*PRTM(JLH,JLEV)  
      ENDIF
      PRTM(JLH,JLEV)=0.0_JPRB
      PT(JLH,JLEV)  =PT(JLH,JLEV) +&
       & (PRV-PRD)*PQ5L(JLH,JLEV)*PRTL(JLH,JLEV)  
      PTL(JLH,JLEV) =PTL(JLH,JLEV) +&
       & PR5(JLH,JLEV)*PRTL(JLH,JLEV)  
      IF (.NOT. LDRY_ECMWF) THEN
        PQL(JLH,JLEV) =PQL(JLH,JLEV) +&
         & (PRV-PRD)*PT5(JLH,JLEV)*PRTL(JLH,JLEV)  
        PR(JLH,JLEV)  =PR(JLH,JLEV) +&
         & PT5L(JLH,JLEV)*PRTL(JLH,JLEV)  
      ENDIF
      PRTL(JLH,JLEV)=0.0_JPRB
    ENDIF
    PT(JLH,JLEV) =PT(JLH,JLEV) +PR5(JLH,JLEV)*PRT(JLH,JLEV)
    IF (.NOT. LDRY_ECMWF) PR(JLH,JLEV)=PR(JLH,JLEV)+PT5(JLH,JLEV)*PRT(JLH,JLEV)
    PRT(JLH,JLEV)=0.0_JPRB
  ENDDO
ENDDO

!----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPRTAD',1,ZHOOK_HANDLE)
END SUBROUTINE GPRTAD
