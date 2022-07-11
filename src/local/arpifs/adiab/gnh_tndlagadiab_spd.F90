!OCL  NOEVAL
SUBROUTINE GNH_TNDLAGADIAB_SPD(   KPDVAR, YDGEOMETRY, KFLEV, KPROMA, KSTART, KPROF, PKAP, P3DIVG, &
& PVVEL, PTNDPD  )  

!**** *GNH_TNDLAGADIAB_SPD* - Computation of the adiabatic contribution of
!                             [D (spd) / Dt] in the NH model.
!                             "spd" is the pressure departure variable.

!     Purpose.
!     --------

!     Computation of the adiabatic contribution of [D (spd) / Dt] in the NH model
!     at full levels.

!     For spd=Qcha, expression of this tendency is:
!      [D (Qcha) / Dt]_adiab = - (cp/cv) D3 - omega/prehyd

!**   Interface.
!     ----------
!        *CALL* *GNH_TNDLAGADIAB_SPD(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!           KFLEV     : number of levels.
!           KPROMA    : horizontal dimension.
!           KSTART    : start of work.
!           KPROF     : working length.
!           PKAP      : Kappa=R/cp at full levels.
!           P3DIVG    : total divergence D3 at full levels.
!           PVVEL     : omega/prehyd at full levels.

!         * OUTPUT:
!           PTNDPD    : tendency [D (spd) / Dt]_adiab at full levels.

!         Remark: input data PKAP,P3DIVG,PVVEL are at time t.

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
!        K. YESSAD (after routine GNHPDVD).
!        Original : 08-Dec-2004

!     Modifications.
!     --------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK


!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM) ,INTENT(IN)   :: KPDVAR
TYPE(GEOMETRY)     ,INTENT(IN)   :: YDGEOMETRY
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KFLEV
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KPROMA 
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KSTART 
INTEGER(KIND=JPIM) ,INTENT(IN)   :: KPROF 
REAL(KIND=JPRB)    ,INTENT(IN)   :: PKAP(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)   :: P3DIVG(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(IN)   :: PVVEL(KPROMA,KFLEV)
REAL(KIND=JPRB)    ,INTENT(OUT)  :: PTNDPD(KPROMA,KFLEV)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,JLEV
REAL(KIND=JPRB) :: ZCPCV

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNH_TNDLAGADIAB_SPD', 0, ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF [D (spd) / Dt]_adiab
!              -----------------------------------

IF (KPDVAR==2) THEN
  DO JLEV=1,KFLEV
    DO JROF=KSTART,KPROF
      ZCPCV=1.0_JPRB/(1.0_JPRB-PKAP(JROF,JLEV))
      PTNDPD(JROF,JLEV)=-(ZCPCV*P3DIVG(JROF,JLEV)+PVVEL(JROF,JLEV))
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNH_TNDLAGADIAB_SPD', 1, ZHOOK_HANDLE)
END SUBROUTINE GNH_TNDLAGADIAB_SPD
