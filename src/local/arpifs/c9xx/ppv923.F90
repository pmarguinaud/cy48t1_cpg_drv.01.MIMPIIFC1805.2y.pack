SUBROUTINE PPV923(YDGEOMETRY,KJ,KLOND,KLATD,KINDMIN,PMD,PLON,PLAT)

!**** *PPV923*

!     PURPOSE.
!     --------

!       Searching for suitable data when there is a conflict
!       between initial and final masks.

!**   INTERFACE.
!     ----------

!      CALL PPV923(KJ,KLOND,KLATD,KINDMIN,PMD,PLON,PLAT)

!        KJ      = index of the point to treat for the final grid
!        KLOND   = number of longitudes of the initial grid
!        KLATD   = number of latitudes of the initial grid
!        KINDMIN = position of the nearest point (output)
!        PMD     = initial mask for missing data
!        PLON    = longitudes of the initial grid
!        PLAT    = latitudes of the initial grid

!     METHOD.
!     -------

!     Calculation of the distance between the ARPEGE/ALADIN gridpoint
!     and all the points of the dataset. Then research of the nearest one.

!     angular distance = arcos(sin(latA).sin(latB)+cos(latA).cos(latB).cos(lonB-lonA))

!     EXTERNALS.
!     ----------

!     AUTHORS.
!     --------

!       F. TAILLEFER    02/12/2005

!     MODIFICATIONS.
!     --------------
!     G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!     T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCT0   , ONLY : LELAM
USE YOMCLI   , ONLY : YRCLI

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KLOND,KLATD,KJ

INTEGER(KIND=JPIM),INTENT(OUT) :: KINDMIN(2)

REAL(KIND=JPRB),INTENT(IN) :: PMD(KLOND,KLATD)
REAL(KIND=JPRB),INTENT(IN) :: PLON(KLOND,KLATD),PLAT(KLOND,KLATD)

INTEGER(KIND=JPIM) :: JX,JY,I1,IXFING,IX,IY

REAL(KIND=JPRB) :: ZD(KLOND,KLATD)
REAL(KIND=JPRB) :: Z1(KLOND,KLATD),Z2(KLOND,KLATD),Z3(KLOND,KLATD)
REAL(KIND=JPRB) :: Z4(KLOND,KLATD),Z5(KLOND,KLATD)
REAL(KIND=JPRB) :: ZLONF,ZLATF,ZS1,ZC1,ZEPS
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPV923',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB)
ASSOCIATE(NDLUNG=>YDDIM%NDLUNG, NDLUXG=>YDDIM%NDLUXG, &
 & NTSTAGP=>YDGEM%NTSTAGP)

IF (LELAM) THEN
  ZEPS=1.E-8_JPRB
  IXFING=NDLUXG-NDLUNG+1
  I1=INT((REAL(KJ-1,JPRB)+ZEPS)/REAL(IXFING,JPRB))
  IY=I1+1
  IX=KJ-(I1*IXFING)
  I1=IX+NTSTAGP(IY)-1
ELSE
  I1=KJ
ENDIF
ZLONF=YDGSGEOM_NB%GELAM(I1)
ZLATF=YDGSGEOM_NB%GELAT(I1)

ZS1=MAX(-1.0_JPRB,MIN(1.0_JPRB,SIN(ZLATF)))
ZC1=MAX(-1.0_JPRB,MIN(1.0_JPRB,COS(ZLATF)))

ZD(:,:)=1.E10_JPRB
DO JY = 1,KLATD
  DO JX = 1,KLOND
    IF (PMD(JX,JY) < YRCLI%SMASK) THEN
      Z1(JX,JY) = MAX(-1.0_JPRB,MIN(1.0_JPRB,SIN(PLAT(JX,JY))))
      Z2(JX,JY) = MAX(-1.0_JPRB,MIN(1.0_JPRB,ZS1*Z1(JX,JY)))
      Z3(JX,JY) = MAX(-1.0_JPRB,MIN(1.0_JPRB,COS(PLAT(JX,JY))))
      Z4(JX,JY) = MAX(-1.0_JPRB,MIN(1.0_JPRB,COS(PLON(JX,JY)-ZLONF)))
      Z5(JX,JY) = MAX(-1.0_JPRB,MIN(1.0_JPRB,ZC1*Z3(JX,JY)*Z4(JX,JY)))
      ZD(JX,JY) = ACOS(Z2(JX,JY)+Z5(JX,JY))
    ENDIF
  ENDDO
ENDDO

KINDMIN(:) = MINLOC(ZD)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('PPV923',1,ZHOOK_HANDLE)
END SUBROUTINE PPV923
