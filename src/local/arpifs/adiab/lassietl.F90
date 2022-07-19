SUBROUTINE LASSIETL(YDGEOMETRY,YDGMV,YDGMV5,YGFL,YDDYN,KST,KPROF,PRCORI,PGMV,PGMVS,PGFL,&
 & PSDIV0,PTOD0,PGAGT0L,PGAGT0M,PGMV5,PGFL5)

!**** *LASSIETL* Semi-Lagrangian scheme.  (tangent-linear version)
!                Computation of linear terms used in the semi-implicit scheme. 

!     Purpose.
!     --------
!        Computation of linear terms used in the semi-implicit scheme:
!        Nabla(Gamma*T+Mu*Pi), (Tau*D) and (Nu*D).

!**   Interface.
!     ----------
!        *CALL* *LASSIETL(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST      - first element of work.
!          KPROF    - depth of work.
!          PRCORI   - Coriolis parameter.
!          PGMV     - GMV variables at t-dt and t.
!          PGMVS    - GMVS variables at t-dt and t.
!          PGFL     - unified_treatment grid-point (GFL) fields.

!        OUTPUT:
!          PSDIV0   - semi-implicit term at time t for continuity equation
!                     (Nu*D).
!          PTOD0    - semi-implicit term at time t for temperature equation
!                     (Tau*D).
!          PGAGT0L  - semi-implicit term at time t for U-wind equation
!                     (zonal component of Nabla(Gamma*T+Mu*Pi)).
!          PGAGT0M  - semi-implicit term at time t for V-wind equation
!                     (meridian component of Nabla(Gamma*T+Mu*Pi)).

!        TRAJECTORY INPUT:
!          PGMV5    - trajectory for PGMV.
!          PGFL5    - trajectory for PGFL.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           Called by LACDYNTL.

!     Reference.
!     ----------
!             Arpege documentation about semi-lagrangian scheme.

!     Author.
!     -------
!        C. Temperton (ECMWF)

!     Modifications.
!     --------------
!     Original : 99/07/20
!     Modified 01-01-23 by C. Temperton: case LSPRT=.T. corrected.
!     Modified by A. Untch 2000-12 : Vertical finite element scheme
!     M.Hamrud      01-Oct-2003 CY28 Cleaning
!     M.Jidane  08-04-2006 : Reintro of R dep. on q variables under key
!     K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!     K. Yessad (Nov 2009): cleanings, remove useless calculations.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : RD, RV, YRCST
USE YOMCT0       , ONLY : LTWOTL, LSPRT
USE YOMDYN       , ONLY : TDYN
USE YOM_YGFL     , ONLY : TYPE_GFLD
USE YOMDYNA      , ONLY : YRDYNA
USE YOMCVER      , ONLY : LVERTFE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV5
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCORI(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMVS(YDGEOMETRY%YRDIM%NPROMA,YDGMV%NDIMGMVS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YGFL%NDIM)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSDIV0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTOD0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAGT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGAGT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMV5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV5%YT5%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YGFL%NDIM5)
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZR0  (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZR5  (YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSDIV5(YDGEOMETRY%YRDIM%NPROMA),ZTOD5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: IPROFS, JLEV, JROF

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "sigam.intfb.h"
#include "sitnu.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LASSIETL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV)
ASSOCIATE(YQ=>YGFL%YQ,   NPROMA=>YDDIM%NPROMA,   NFLEVG=>YDDIMV%NFLEVG,   LIMPF=>YDDYN%LIMPF,   YT0=>YDGMV%YT0,   &
& YT5=>YDGMV5%YT5)
!     ------------------------------------------------------------------

IF (.NOT.LTWOTL) THEN
  CALL ABOR1(' LASSIETL: CALLED WITH LTWOTL=.F.')
ENDIF

IPROFS=KPROF-KST+1

!     ------------------------------------------------------------------

!      1.   TRAJECTORY.

! * Trajectory
!  - Computation of Nu*D (SI term for continuity equation)
!     and Tau*D (SI term for temperature equation).
IF (LSPRT) THEN
  CALL SITNU(LVERTFE, YRCST,YDGEOMETRY,YDDYN,NPROMA,1,PGMV5(KST,1,YT5%MDIV),ZTOD5(KST,1),ZSDIV5(KST),IPROFS)
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      ZR5(JROF,JLEV)=RD+(RV-RD)*PGFL5(JROF,JLEV,YQ%MP5)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!      2.   TL CODE.

!     * Computation of Nu*D, Tau*D, Nabla(Gamma*T+Mu*Pi).

!   - Computation of Nu*D (SI term for continuity equation)
!     and Tau*D (SI term for temperature equation).
CALL SITNU(LVERTFE, YRCST,YDGEOMETRY,YDDYN,NPROMA,1,PGMV(KST,1,YT0%MDIV),PTOD0(KST,1),PSDIV0(KST),IPROFS)
!   - Computation of Nabla(Gamma*T+Mu*Pi) (SI term for momentum equation).
CALL SIGAM(YRCST,YDGEOMETRY,YDDYN,NPROMA,1,PGAGT0L(KST,1),PGMV(KST,1,YT0%MTL),PGMVS(KST,YT0%MSPL),&
 & IPROFS,NFLEVG)
CALL SIGAM(YRCST,YDGEOMETRY,YDDYN,NPROMA,1,PGAGT0M(KST,1),PGMV(KST,1,YT0%MTM),PGMVS(KST,YT0%MSPM),&
 & IPROFS,NFLEVG)

!     * Modify (Tau*D) if (LSPRT=.T.) in temperature equation to compensate
!       for later multiplication by R/Rd:

IF (LSPRT) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      PTOD0(JROF,JLEV)=RD*PTOD0(JROF,JLEV)/ZR5(JROF,JLEV)  
    ENDDO
    IF (.NOT. YRDYNA%LDRY_ECMWF) THEN
      DO JROF=KST,KPROF
        ZR0(JROF,JLEV)=(RV-RD)*PGFL(JROF,JLEV,YQ%MP)
        PTOD0(JROF,JLEV)=PTOD0(JROF,JLEV)-RD*ZTOD5(JROF,JLEV)*ZR0(JROF,JLEV)&
          & /(ZR5(JROF,JLEV)*ZR5(JROF,JLEV)) 
      ENDDO
    ENDIF
  ENDDO
ENDIF

!     * Add semi-implicit Coriolis terms to Nabla(Gamma*T+Mu*Pi)
!       if required (LIMPF=.T.).

IF (LIMPF) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      PGAGT0L(JROF,JLEV)=PGAGT0L(JROF,JLEV)-PRCORI(JROF)*PGMV(JROF,JLEV,YT0%MV)
      PGAGT0M(JROF,JLEV)=PGAGT0M(JROF,JLEV)+PRCORI(JROF)*PGMV(JROF,JLEV,YT0%MU)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LASSIETL',1,ZHOOK_HANDLE)
END SUBROUTINE LASSIETL

