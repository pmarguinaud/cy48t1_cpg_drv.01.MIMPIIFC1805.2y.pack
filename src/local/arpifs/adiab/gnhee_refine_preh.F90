!OCL  NOEVAL
SUBROUTINE GNHEE_REFINE_PREH(&
 ! --- INPUT -----------------------------------------------------------------
 & YDCST, YDGEOMETRY,KSTART,KPROF,&
 & POROGL,POROGM,POROGLM,POROGLL,POROGMM,&
 & PUS,PVS,&
 & PREF,PREH,PDEP,PWH2F,&
 ! --- INPUT-OUTPUT ----------------------------------------------------------
 & PSGRTL,PSGRTM,PSGRTSL,PSGRTSM, &
 & PTNDU,PTNDV,PTNDU_NOC,PTNDV_NOC,&
 & PTNDUS,PTNDVS,PNHPREH)  

!**** *GNHEE_REFINE_PREH* - Using an iterative algorithm:
!                           * refine calculation of half level "pre" and full level [Delta pre]
!                           * refine calculation of pressure gradient term in RHS of DV/Dt.
!                           in NHEE model.

!     Purpose.
!     --------

!     Using an iterative algorithm:
!      * refine calculation of half level "pre" and full level [Delta pre]
!      * refine calculation of pressure gradient term in RHS of DV/Dt, and increment DV/Dt.
!     in NHEE model.

!     In particular, final value of [Delta pre]_surf/[Delta prehyd]_surf must match:
!      [Delta pre]_surf/[Delta prehyd]_surf = 1 + (1/g**2) D(Vs grad(Phi_s))/Dt

!**   Interface.
!     ----------
!        *CALL* *GNHEE_REFINE_PREH(...)

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
!           PUS        : surface U wind.
!           PVS        : surface V wind.
!           PREF       : "prehyd" at full levels.
!           PREH       : "prehyd" at half levels.
!           PDEP       : (pre-prehyd) at full levels.
!           PWH2F      : (pre-prehyd) at full levels.

!         * INPUT-OUTPUT:
!           PSGRTL     : zonal component of the pressure force grad.
!           PSGRTM     : meridian component of the pressure force grad.
!           PSGRTSL    : zonal component of the surface pressure force grad.
!           PSGRTSM    : meridian component of the surface pressure force grad.
!           PTNDU      : [DU/Dt].
!           PTNDV      : [DV/Dt].
!           PTNDU_NOC  : [DU/Dt] without Coriolis term.
!           PTNDV_NOC  : [DV/Dt] without Coriolis term.
!           PTNDUS     : [DU/Dt]_surf.
!           PTNDVS     : [DV/Dt]_surf.
!           PNHPREH    : total pressure "pre" at half levels.



!         * OUTPUT:

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
!        Original : Sept 2018

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : TCST
USE YOMCVER      , ONLY : LVERTFE, LVFE_GW

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGLM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGLL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGMM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDEP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWH2F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSGRTL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSGRTM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSGRTSL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSGRTSM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTNDU(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTNDV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTNDU_NOC(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTNDV_NOC(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTNDUS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTNDVS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNHPREH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)


!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF
REAL(KIND=JPRB) :: ZRG2,ZJACOB,ZACUGF,ZACVGF,ZWDENO,ZINS
REAL(KIND=JPRB) :: ZDELPSDPREHYD_NEW(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZDELPSDPREHYD_OLD(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZINCRL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZINCRM(YDGEOMETRY%YRDIM%NPROMA)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

#include "abor1.intfb.h"    

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHEE_REFINE_PREH',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA,NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

ZRG2=YDCST%RG*YDCST%RG

!     ------------------------------------------------------------------

IF (LVERTFE.AND.LVFE_GW) CALL ABOR1(' GNHEE_REFINE_PREH 3.1.20c: option not yet coded!')

! * Update [Delta (pre-prehyd)]_surf/[Delta prehyd]_surf:

DO JROF=KSTART,KPROF
  ZJACOB=POROGLL(JROF)*PUS(JROF)*PUS(JROF) &
   & + POROGMM(JROF)*PVS(JROF)*PVS(JROF) &
   & + 2.0_JPRB*POROGLM(JROF)*PUS(JROF)*PVS(JROF)
  ZACUGF=PTNDUS(JROF)*POROGL(JROF)
  ZACVGF=PTNDVS(JROF)*POROGM(JROF)
  ZWDENO=1.0_JPRB+PWH2F(JROF,NFLEVG,2)*(((POROGL(JROF)**2)+(POROGM(JROF)**2))/ZRG2)
  ZINS  =(ZACUGF+ZACVGF+ZJACOB)/ZRG2 
  ZDELPSDPREHYD_OLD(JROF)=PDEP(JROF,NFLEVG)*((PREH(JROF,NFLEVG)/PREF(JROF,NFLEVG))-1.0_JPRB) &
             & /(PREH(JROF,NFLEVG)-PREF(JROF,NFLEVG))
  ZDELPSDPREHYD_NEW(JROF)=ZDELPSDPREHYD_OLD(JROF) &
             & +((ZINS-ZDELPSDPREHYD_OLD(JROF))/ZWDENO)
ENDDO

! * Update surface total pressure using bottom condition.

DO JROF=KSTART,KPROF
  PNHPREH(JROF,NFLEVG)=(PDEP(JROF,NFLEVG)+PREF(JROF,NFLEVG)) &
        &+(1.0_JPRB+ZDELPSDPREHYD_NEW(JROF))*(PREH(JROF,NFLEVG)-PREF(JROF,NFLEVG))
ENDDO

! * Update horizontal pressure gradient term and DV/Dt at surface

DO JROF=KSTART,KPROF
   ZINCRL(JROF)=PWH2F(JROF,NFLEVG,2)*(ZDELPSDPREHYD_NEW(JROF) &
               & -ZDELPSDPREHYD_OLD(JROF))*POROGL(JROF)
   ZINCRM(JROF)=PWH2F(JROF,NFLEVG,2)*(ZDELPSDPREHYD_NEW(JROF) & 
               & -ZDELPSDPREHYD_OLD(JROF))*POROGM(JROF)

   PSGRTL(JROF,NFLEVG)=PSGRTL(JROF,NFLEVG)+ZINCRL(JROF)
   PSGRTM(JROF,NFLEVG)=PSGRTM(JROF,NFLEVG)+ZINCRM(JROF)
   PTNDU(JROF,NFLEVG)=PTNDU(JROF,NFLEVG)-ZINCRL(JROF)
   PTNDU_NOC(JROF,NFLEVG)=PTNDU_NOC(JROF,NFLEVG)-ZINCRL(JROF)
   PTNDV(JROF,NFLEVG)=PTNDV(JROF,NFLEVG)-ZINCRM(JROF)
   PTNDV_NOC(JROF,NFLEVG)=PTNDV_NOC(JROF,NFLEVG)-ZINCRM(JROF)
   
   PSGRTSL(JROF)=PSGRTL(JROF,NFLEVG)
   PSGRTSM(JROF)=PSGRTM(JROF,NFLEVG)
   PTNDUS(JROF)=PTNDUS(JROF)-ZINCRL(JROF)
   PTNDVS(JROF)=PTNDVS(JROF)-ZINCRM(JROF)
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHEE_REFINE_PREH',1,ZHOOK_HANDLE)
END SUBROUTINE GNHEE_REFINE_PREH
