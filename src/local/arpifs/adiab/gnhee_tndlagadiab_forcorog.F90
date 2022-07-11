!OCL  NOEVAL
SUBROUTINE GNHEE_TNDLAGADIAB_FORCOROG(&
 ! --- INPUT -----------------------------------------------------------------
 & LDVEREGINT,YDCST,YDGEOMETRY,KSTART,KPROF,&
 & POROGL,POROGM,POROGLM,POROGLL,POROGMM,&
 & PUH,PVH,PTNDUF,PTNDVF,&
 & PEVEL,PRDELP,PLNPR,PALPH,&
 ! --- OUTPUT ----------------------------------------------------------------
 & PATND13)  

!**** *GNHEE_TNDLAGADIAB_FORCOROG* - Computation of the adiabatic contribution of
!                                    [1/g**2] [D (S V grad(Phi_s}) / Dt] in the NHEE model

!     Purpose.
!     --------

!     Computation of the adiabatic contribution of [1/g**2] [D (S V grad(Phi_s}) / Dt] in the NHEE model,
!     at full levels.
!     [D (S V grad(Phi_s}) / Dt] is the sum of three contributions:
!      * [etadot DS/Deta V grad(Phi_s)]: contains the vertical advection of S.
!      * [S DV/Dt grad(Phi_s)]: contains the adiabatic RHS of momentum equation.
!      * [S V.V.grad(grad(Phi_s))]: upper air Jacobian term, uses second order derivatives of "surf orography"

!**   Interface.
!     ----------
!        *CALL* *GNHEE_TNDLAGADIAB_FORCOROG(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!           YDGEOMETRY : structure containing geometry.
!           KSTART     : start of work.
!           KPROF      : working length.
!           POROGL     : zonal comp. of grad(surf orography).
!           POROGM     : merid comp. of grad(surf orography).
!           POROGLM    : -I
!           POROGLL    :  I- second order derivatives of "surf orography".
!           POROGMM    : -I
!           PUF        : full level U wind.
!           PVF        : full level V wind.
!           PTNDUF     : full level [DU/Dt].
!           PTNDVF     : full level [DV/Dt].
!           PEVEL      : [etadot dprehyd/deta] at half levels.
!           PRDELP     : 1/[Delta prehyd] at full levels.
!           PLNPR      : "delta" at full levels.
!           PALPH      : "alpha" at full levels.

!         * OUTPUT:
!           PATND13    : [1/g**2] [D (S V grad(Phi_s}) / Dt]

!        Implicit arguments :   None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    None.
!     ----------

!     Reference.
!     ----------
!        See documentations Arpege about the NH model.

!     Author.
!     -------
!        K. YESSAD.
!        Original : August 2018

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : TCST
!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL            ,INTENT(IN)   :: LDVEREGINT
TYPE(TCST)         ,INTENT(IN)   :: YDCST
TYPE(GEOMETRY)     ,INTENT(IN)   :: YDGEOMETRY
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KSTART 
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KPROF 
REAL(KIND=JPRB)    ,INTENT(IN)   :: POROGL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)   :: POROGM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)   :: POROGLM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)   :: POROGLL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)   :: POROGMM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PUH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PVH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PTNDUF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PTNDVF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PEVEL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PRDELP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PLNPR(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PALPH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    ,INTENT(OUT)  :: PATND13(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,JLEV

REAL(KIND=JPRB) :: ZRG2,ZJACOB
REAL(KIND=JPRB) :: ZTNDWWI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTNDUH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTNDVH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZZ1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZZ2(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZZ3(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

#include "gphluv_expl.intfb.h"
#include "gphlwi.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHEE_TNDLAGADIAB_FORCOROG',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA,NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,&
& YDVAB=>YDGEOMETRY%YRVAB,YDDIMV=>YDGEOMETRY%YRDIMV)

!     ------------------------------------------------------------------

ZRG2=YDCST%RG*YDCST%RG

!    ------------------------------------------------------------------

!*       1.    PRELIMINARY CALCULATIONS D[U,V]/DT AT HALF-LEVELS
!              ---------------------------------------------

! compute interpolation weights
CALL GPHLWI(YDDIMV,NPROMA,KSTART,KPROF,PLNPR,&
 & PALPH,ZTNDWWI,LDVERINT=LDVEREGINT)
! interpolate wind tendencies at half levels
CALL GPHLUV_EXPL(YDDIMV,NPROMA,KSTART,KPROF,PTNDUF,PTNDVF,ZTNDWWI,ZTNDUH,ZTNDVH)

!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF [etadot DS/Deta V grad(Phi_s)]
!              ---------------------------------------------

DO JLEV=1,NFLEVG-1
  DO JROF=KSTART,KPROF
    ZZ1(JROF,JLEV)=0.5_JPRB*PEVEL(JROF,JLEV) &
     & *((YDVAB%VRATH(JLEV+1)-YDVAB%VRATH(JLEV))*PRDELP(JROF,JLEV+1)&
     & +(YDVAB%VRATH(JLEV)-YDVAB%VRATH(JLEV-1))*PRDELP(JROF,JLEV))  &
     & *(PUH(JROF,JLEV)*POROGL(JROF)+PVH(JROF,JLEV)*POROGM(JROF))
  ENDDO
ENDDO
ZZ1(KSTART:KPROF,NFLEVG)=0.0_JPRB

!     ------------------------------------------------------------------

!*       2.    COMPUTATION OF [S DV/Dt grad(Phi_s)]
!              ------------------------------------

DO JLEV=1,NFLEVG
  DO JROF=KSTART,KPROF
    ZZ2(JROF,JLEV)=YDVAB%VRATH(JLEV)*(ZTNDUH(JROF,JLEV)*POROGL(JROF) &
                  & + ZTNDVH(JROF,JLEV)*POROGM(JROF))
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF [S V.V.grad(grad(Phi_s))]
!              ----------------------------------------

DO JLEV=1,NFLEVG
  DO JROF=KSTART,KPROF
    ZJACOB=POROGLL(JROF)*PUH(JROF,JLEV)*PUH(JROF,JLEV)&
     & +POROGMM(JROF)*PVH(JROF,JLEV)*PVH(JROF,JLEV)&
     & +2.0_JPRB*POROGLM(JROF)*PUH(JROF,JLEV)*PVH(JROF,JLEV)
    ZZ3(JROF,JLEV)=YDVAB%VRATH(JLEV)*ZJACOB
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*       4.    FILL PATND13
!              ------------

DO JLEV=1,NFLEVG
  DO JROF=KSTART,KPROF
    PATND13(JROF,JLEV)=(ZZ1(JROF,JLEV)+ZZ2(JROF,JLEV)+ZZ3(JROF,JLEV))/ZRG2
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHEE_TNDLAGADIAB_FORCOROG',1,ZHOOK_HANDLE)
END SUBROUTINE GNHEE_TNDLAGADIAB_FORCOROG
