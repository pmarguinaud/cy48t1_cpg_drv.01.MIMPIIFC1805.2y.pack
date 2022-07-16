SUBROUTINE GPRTTL(LDSPRT,KPROMA,KSTART,KEND,KLEV,PRD,PRV,&
 & PR5,PT5,PT5L,PT5M,PQ5L,PQ5M,&
 & PR,PT,PTL,PTM,PQL,PQM,PRT,PRTL,PRTM)  

!*    *GPRTTL* 

!     Purpose
!     -------  To calculate RT and its derivates (TL code)

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
!        M.Jidane  08-04-2006 : Reintro of R dep. on q variables under key
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PR(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTL(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQL(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRT(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRTL(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRTM(KPROMA,KLEV) 

!----------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JLH
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!----------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPRTTL',0,ZHOOK_HANDLE)
!----------------------------------------------------------
!---------------------------------------------------------------
!*       1. Compute RT and its derivatives

DO JLEV=1,KLEV
  DO JLH=KSTART,KEND
    PRT(JLH,JLEV)= PR5(JLH,JLEV)*PT(JLH,JLEV)
    IF (.NOT. LDRY_ECMWF) PRT(JLH,JLEV)=PRT(JLH,JLEV)+PR(JLH,JLEV)*PT5(JLH,JLEV)
    IF (LDSPRT) THEN
      PRTL(JLH,JLEV)=PRD*PTL(JLH,JLEV)
      PRTM(JLH,JLEV)=PRD*PTM(JLH,JLEV)
    ELSE
      PRTL(JLH,JLEV)= (PRV-PRD)*PT(JLH,JLEV)*PQ5L(JLH,JLEV)&
       & + PR5(JLH,JLEV)*PTL(JLH,JLEV)
      IF (.NOT.LDRY_ECMWF) PRTL(JLH,JLEV)=PRTL(JLH,JLEV)+(PRV-PRD)*PT5(JLH,JLEV)*PQL(JLH,JLEV)&
       & + PR(JLH,JLEV)*PT5L(JLH,JLEV)
      PRTM(JLH,JLEV)= (PRV-PRD)*PT(JLH,JLEV)*PQ5M(JLH,JLEV)&
       & + PR5(JLH,JLEV)*PTM(JLH,JLEV)
      IF (.NOT.LDRY_ECMWF) PRTM(JLH,JLEV)=PRTM(JLH,JLEV)+(PRV-PRD)*PT5(JLH,JLEV)*PQM(JLH,JLEV)&
       & + PR(JLH,JLEV)*PT5M(JLH,JLEV)
    ENDIF
  ENDDO
ENDDO

!----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPRTTL',1,ZHOOK_HANDLE)
END SUBROUTINE GPRTTL
