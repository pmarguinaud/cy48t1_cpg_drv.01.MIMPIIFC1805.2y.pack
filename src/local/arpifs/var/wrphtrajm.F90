SUBROUTINE WRPHTRAJM(YDGEOMETRY,YDSIMPHL,KSTA,KEND,PTRAJ_PHYS,&
 & PUT95,PVT95,PTT95,PQT95,PQLT95,PQIT95,PSP95)  

!**** *WRPHTRAJM* - Write out the trajectory at t-dt to work file

!     Purpose.  Write out the trajectory at t-dt to work file

!**   Interface.
!     ----------
!        *CALL* *WRPHTRAJM(...)*

!        Explicit arguments :  
!        --------------------  KSTA    - start adress         (input)
!                              KEND    - end adress           (input)
!                              PTRAJ_PHYS - trajectory        (output)
!                              PUT95 to PSP95: trajectory to be written
!                                        in a traject. field  (input)

!        Implicit arguments :  None
!        --------------------

!     Method.
!     -------

!     Externals. FILLB    - interface to packing routines
!     ---------- 

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Jean-Francois Mahfouf  *ECMWF*  95-10-18

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!      O.Riviere  Feb 11 Storage of PQL5 and PQI5
!      R. El Khatib  21-Apr-2011 Use fillb3 for 3D fields
!      F. Vana  28-Nov-2013 : Redesigned trajectory handling
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMSIMPHL, ONLY : TSIMPHL
USE YOMTRAJ  , ONLY : TRAJ_PHYS_TYPE
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)       ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSIMPHL)        ,INTENT(IN)    :: YDSIMPHL
INTEGER(KIND=JPIM)   ,INTENT(IN)    :: KSTA 
INTEGER(KIND=JPIM)   ,INTENT(IN)    :: KEND 
TYPE (TRAJ_PHYS_TYPE),INTENT(INOUT) :: PTRAJ_PHYS
REAL(KIND=JPRB)      ,INTENT(IN)    :: PUT95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)      ,INTENT(IN)    :: PVT95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)      ,INTENT(IN)    :: PTT95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)      ,INTENT(IN)    :: PQT95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)      ,INTENT(IN)    :: PQLT95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)      ,INTENT(IN)    :: PQIT95(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)      ,INTENT(IN)    :: PSP95(YDGEOMETRY%YRDIM%NPROMA) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WRPHTRAJM',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & NPROMA=>YDDIM%NPROMA, &
 & LPROCLDTL=>YDSIMPHL%LPROCLDTL)
!     ------------------------------------------------------------------

!*       1.    WRITE OUT GRIDPOINT ARRAY.
!              --------------------------

!                 UPPER AIR FIELDS
PTRAJ_PHYS%PUT9MF5(KSTA:KEND,1:NFLEVG)   =PUT95(KSTA:KEND,1:NFLEVG)
PTRAJ_PHYS%PVT9MF5(KSTA:KEND,1:NFLEVG)   =PVT95(KSTA:KEND,1:NFLEVG)
PTRAJ_PHYS%PTT9MF5(KSTA:KEND,1:NFLEVG)   =PTT95(KSTA:KEND,1:NFLEVG)
PTRAJ_PHYS%PQT9MF5(KSTA:KEND,1:NFLEVG)   =PQT95(KSTA:KEND,1:NFLEVG)
IF(LPROCLDTL) THEN
  PTRAJ_PHYS%PQLT9MF5(KSTA:KEND,1:NFLEVG)=PQLT95(KSTA:KEND,1:NFLEVG)
  PTRAJ_PHYS%PQIT9MF5(KSTA:KEND,1:NFLEVG)=PQIT95(KSTA:KEND,1:NFLEVG)
ENDIF

!                 SURFACE 2D FIELDS
PTRAJ_PHYS%PSP9MF5(KSTA:KEND) = PSP95(KSTA:KEND)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('WRPHTRAJM',1,ZHOOK_HANDLE)
END SUBROUTINE WRPHTRAJM

