SUBROUTINE LACDYNSHWAD(YDGEM,YDDYN,CDCONF,KPROMA,KSTART,KPROF,PDT,&
 & PUT9,PVT9,&
 & PUT0,PVT0,PSPT0,PDIVT0,PSPT0L,PSPT0M,&
 & PSPT5,PDIVT5,&
 & YDGSGEOM,POROG,POROGL,POROGM,&
 & PSPNLT9,&
 & PUL9,PVL9,PSPL9,PUL0,PVL0,PSPL0,PSPT1,&
 & PURL0,PVRL0,PUSI,PVSI,PURL,PVRL)  

USE YOMGEM   , ONLY : TGEM
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCT0   , ONLY : LTWOTL
USE YOMCT3   , ONLY : NSTEP
USE YOMDYNA  , ONLY : LELTRA, LSETTLS, LSETTLST
USE YOMDYN   , ONLY : TDYN
USE YOMGSGEOM, ONLY : TGSGEOM

!**** *LACDYNSHWAD* - Semi-Lagrangian scheme. - 2-D shallow-water model.
!                     Computation of the t and t-dt useful quantities
!                     at grid-points.  (Adjoint version)

!     Purpose.
!     --------
!          This subroutine computes the equation quantities at t-dt and t
!          at each grid-point of the colocation grid (Gauss grid).

!**   Interface.
!     ----------
!        *CALL* *LACDYNSHWAD(...)

!        Explicit arguments :
!        --------------------

!        INPUT OR INPUT/OUTPUT MATCHING WITH INPUT IN THE TL CODE:
!          CDCONF  - configuration.
!          KPROMA  - horizontal dimension.
!          KSTART  - first element of work.
!          KPROF   - depth of work.
!          PDT     - For a leap-frog scheme (three time level scheme):
!                     'dt' at the first time-step, '2 dt' otherwise.
!                    For a 2TL SL scheme: timestep 'dt'.
!          PUT9    - U-component of the wind at t-dt.
!          PVT9    - V-component of the wind at t-dt.
!          PUT0    - U-component of the wind at t.
!          PVT0    - V-component of the wind at t.
!          PSPT0   - equivalent height "Phi" at t.
!          PDIVT0  - divergence at t.
!          PSPT0L  - zonal derivative of equivalent height "Phi" at t.
!          PSPT0M  - meridian derivative of equivalent height "Phi" at t.
!     ---------------------- trajectory variables -----------------
!          PSPT5   - equivalent height "Phi" at t.
!          PDIVT5  - divergence at t.
!     -------------------------------------------------------------
!          YDGSGEOM- Grid point geometry.
!          POROG   - grid-point surface orography
!          POROGL  - zonal component of the orography gradient
!          POROGM  - merid component of the orography gradient

!        INPUT/OUTPUT MATCHING WITH INPUT/OUTPUT OR OUTPUT IN THE TL CODE:
!          PSPNLT9 - nonlinear term in continuity equation at t-dt.
!          PUL9    - SL quantity to be interpolated at O for U-wind eqn.
!          PVL9    - SL quantity to be interpolated at O for V-wind eqn.
!          PSPL9   - SL quantity to be interpolated at O for continuity eqn.
!          PUL0    - second SL quantity to be
!                    interpolated (linearly) at O if KWLAG=3,
!                    for the U-wind equation.
!          PVL0    - second SL quantity to be
!                    interpolated (linearly) at O if KWLAG=3,
!                    for the V-wind equation.
!          PSPL0   - second SL quantity to be
!                    interpolated (linearly) at O if NVLAG=3,
!                    for the continuity equation.
!          PSPT1   - t+dt term and other final point terms for continuity eqn.
!     ---------------------- trajectory variables -----------------
!          PSP095  - as PSPL0 but for trajectory.
!     -------------------------------------------------------------
!          PURL0   - U-component of the wind used for the
!                    research of the medium and origin points.
!                    (to be interpolated)
!          PVRL0   - V-component of the wind used for the
!                    research of the medium and origin points.
!                    (to be interpolated)
!          PUSI    - Part of the t semi-implicit term evaluated at the
!                    final point in the U-wind equation.
!          PVSI    - Part of the t semi-implicit term evaluated at the
!                    final point in the V-wind equation.
!          PURL    - U-component of the wind used for the
!                    research of the medium and origin points.
!                    (to be used at the arrival point)
!          PVRL    - V-component of the wind used for the
!                    research of the medium and origin points.
!                    (to be used at the arrival point)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!      None
!      Called by CPG2AD.

!     Reference.
!     ----------
!        Arpege documentation about semi-lagrangian scheme.

!     Author.
!     -------
!      C. Temperton after LACDYNSHW.
!      Original : 98-12-30

!     Modifications.
!     --------------
!      Modified 01-08-30 by K. YESSAD: pruning and some other cleanings.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TGEM)        ,INTENT(IN)    :: YDGEM
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
CHARACTER(LEN=1)  ,INTENT(IN)    :: CDCONF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUT9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVT9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUT0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVT0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPT0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDIVT0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPT0L(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPT0M(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIVT5(KPROMA) 
TYPE(TGSGEOM)     ,INTENT(IN)    :: YDGSGEOM
REAL(KIND=JPRB)   , INTENT(IN)   :: POROG(KPROMA)
REAL(KIND=JPRB)   , INTENT(IN)   :: POROGL(KPROMA)
REAL(KIND=JPRB)   , INTENT(IN)   :: POROGM(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPNLT9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUL9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVL9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPL9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUL0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVL0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPL0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT1(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PURL0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVRL0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PURL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVRL(KPROMA) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZBDT(KPROMA)
REAL(KIND=JPRB) :: ZGAGT0L(KPROMA),ZGAGT0M(KPROMA)
REAL(KIND=JPRB) :: ZMOY1SP(KPROMA)
REAL(KIND=JPRB) :: ZMOYSPNL(KPROMA)
REAL(KIND=JPRB) :: ZSPNLT0(KPROMA)
REAL(KIND=JPRB) :: ZSPSI(KPROMA)

REAL(KIND=JPRB) :: ZBETA, ZBT, ZCMSLP, ZDTS2, ZESGM, ZESGP
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LACDYNSHWAD',0,ZHOOK_HANDLE)
ASSOCIATE(BETADT=>YDDYN%BETADT, LIMPF=>YDDYN%LIMPF, LSIDG=>YDDYN%LSIDG, &
 & NVLAG=>YDDYN%NVLAG, NWLAG=>YDDYN%NWLAG, RCMSLP0=>YDDYN%RCMSLP0, &
 & SIVP=>YDDYN%SIVP, VESL=>YDDYN%VESL, &
 & RSTRET=>YDGEM%RSTRET)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS:
!              ----------------------------

ZDTS2 = 0.5_JPRB*PDT
ZBETA = 0.0_JPRB
IF (CDCONF == 'A'.OR.CDCONF == 'B') THEN
  IF (LTWOTL) THEN
    ZBETA = 1.0_JPRB
  ELSE
    ZBETA = BETADT
  ENDIF
ELSE
  CALL ABOR1(' LACDYNSHWAD: Not existing CDCONF ')
ENDIF
ZBT = ZDTS2*ZBETA

IF (LSIDG) THEN
  ZBDT(KSTART:KPROF)=ZBT*SIVP(1)
ELSE
  ZBDT(KSTART:KPROF)=ZBT*SIVP(1)*RSTRET*RSTRET&
   & /(YDGSGEOM%GM(KSTART:KPROF)*YDGSGEOM%GM(KSTART:KPROF))  
ENDIF

ZESGM=1.0_JPRB-VESL
ZESGP=1.0_JPRB+VESL

ZCMSLP=RCMSLP0

!     ------------------------------------------------------------------

!*       5.    COMPUTATION OF THE CONTINUITY EQUATION RIGHT-HAND SIDE TERMS.
!              -------------------------------------------------------------

IF (LTWOTL) THEN

!       * Two-time-level scheme:
!         Is coded only for conventional formulation of continuity eqn (NVLAG>0)

  IF (NVLAG == 2) THEN
    PSPT0(KSTART:KPROF)=PSPT0(KSTART:KPROF)+PSPL9(KSTART:KPROF)
    ZSPSI(KSTART:KPROF)=-ZESGM*PSPL9(KSTART:KPROF)
    ZMOYSPNL(KSTART:KPROF)=&
     & ZESGM*PSPL9(KSTART:KPROF)+ZESGP*PSPT1(KSTART:KPROF)  
    ZSPNLT0(KSTART:KPROF)=0.0_JPRB
  ELSEIF(NVLAG == 3) THEN
    PSPT0(KSTART:KPROF)=PSPT0(KSTART:KPROF)+PSPL9(KSTART:KPROF)
    IF (LSETTLS) THEN
      ZMOYSPNL(KSTART:KPROF)=PSPL0(KSTART:KPROF)
      ZSPSI(KSTART:KPROF)=-ZESGM*PSPL0(KSTART:KPROF)
      ZSPNLT0(KSTART:KPROF)=ZESGP*PSPT1(KSTART:KPROF)
    ELSE
      ZSPSI(KSTART:KPROF)=-ZESGM*PSPL0(KSTART:KPROF)
      ZMOYSPNL(KSTART:KPROF)=&
       & ZESGM*PSPL0(KSTART:KPROF)+ZESGP*PSPT1(KSTART:KPROF)  
      ZSPNLT0(KSTART:KPROF)=0.0_JPRB
    ENDIF
  ENDIF
  IF (LELTRA) THEN
    ZMOY1SP(KSTART:KPROF)=0.0_JPRB
  ELSE
    ZMOY1SP(KSTART:KPROF)=PSPNLT9(KSTART:KPROF)
    ZSPSI(KSTART:KPROF)=ZSPSI(KSTART:KPROF)+PSPNLT9(KSTART:KPROF)
  ENDIF
  PSPNLT9(KSTART:KPROF)=0.0_JPRB
  IF ((NSTEP <= 0).OR.LELTRA) THEN
    ZMOY1SP(KSTART:KPROF)=&
     & ZMOY1SP(KSTART:KPROF)+ZMOYSPNL(KSTART:KPROF)&
     & +ZSPNLT0(KSTART:KPROF)  
    ZSPSI(KSTART:KPROF)=&
     & ZSPSI(KSTART:KPROF)+ZMOYSPNL(KSTART:KPROF)&
     & +ZSPNLT0(KSTART:KPROF)  
  ELSE
    IF (LSETTLS) THEN
      ZMOY1SP(KSTART:KPROF)=&
       & ZMOY1SP(KSTART:KPROF)+(1.0_JPRB+ZESGM)*ZMOYSPNL(KSTART:KPROF)&
       & +ZSPNLT0(KSTART:KPROF)  
      ZSPSI(KSTART:KPROF)=&
       & ZSPSI(KSTART:KPROF)+(1.0_JPRB+ZESGM)*ZMOYSPNL(KSTART:KPROF)&
       & +ZSPNLT0(KSTART:KPROF)  
      PSPNLT9(KSTART:KPROF)=PSPNLT9(KSTART:KPROF)-ZMOYSPNL(KSTART:KPROF)
    ELSE
      ZMOY1SP(KSTART:KPROF)=&
       & ZMOY1SP(KSTART:KPROF)+1.5_JPRB*ZMOYSPNL(KSTART:KPROF)  
      ZSPSI(KSTART:KPROF)=ZSPSI(KSTART:KPROF)+1.5_JPRB*ZMOYSPNL(KSTART:KPROF)
      PSPNLT9(KSTART:KPROF)=&
       & PSPNLT9(KSTART:KPROF)-0.5_JPRB*ZMOYSPNL(KSTART:KPROF)  
    ENDIF
  ENDIF
  PDIVT0(KSTART:KPROF)=&
   & PDIVT0(KSTART:KPROF)+ZBDT(KSTART:KPROF)*ZSPSI(KSTART:KPROF)&
   & -ZDTS2*(PSPT5(KSTART:KPROF)-POROG(KSTART:KPROF))&
   & *ZMOY1SP(KSTART:KPROF)  
  PSPT0(KSTART:KPROF)=PSPT0(KSTART:KPROF)&
   & -ZDTS2*PDIVT5(KSTART:KPROF)*ZMOY1SP(KSTART:KPROF)  
  PUT0(KSTART:KPROF)=PUT0(KSTART:KPROF)&
   & +ZDTS2*ZCMSLP*POROGL(KSTART:KPROF)*ZMOY1SP(KSTART:KPROF)  
  PVT0(KSTART:KPROF)=PVT0(KSTART:KPROF)&
   & +ZDTS2*ZCMSLP*POROGM(KSTART:KPROF)*ZMOY1SP(KSTART:KPROF)  

ELSE    !!! AD not done yet !!!

ENDIF

!     ------------------------------------------------------------------

!*       4.    COMPUTATION OF THE MOMENTUM EQUATION RIGHT-HAND SIDE TERMS.
!              -----------------------------------------------------------

IF (LTWOTL) THEN

!       * Two-time-level scheme:
!         The non linear terms ZMOY1U+PUSI, ZMOY1V+PVSI are assumed to be
!         zero in this case (BETADT=1, LADVF=T or LIMPF=T), so simplifications 
!         are already made. 

  IF (NWLAG == 2) THEN
    PUSI(KSTART:KPROF)=ZESGP*PUSI(KSTART:KPROF)
    PVSI(KSTART:KPROF)=ZESGP*PVSI(KSTART:KPROF)
    PUT0(KSTART:KPROF)=PUT0(KSTART:KPROF)+PUL9(KSTART:KPROF)
    PVT0(KSTART:KPROF)=PVT0(KSTART:KPROF)+PVL9(KSTART:KPROF)
    PUSI(KSTART:KPROF)=PUSI(KSTART:KPROF)-ZESGM*PUL9(KSTART:KPROF)
    PVSI(KSTART:KPROF)=PVSI(KSTART:KPROF)-ZESGM*PVL9(KSTART:KPROF)
  ELSEIF(NWLAG == 3) THEN
    PUSI(KSTART:KPROF)=ZESGP*PUSI(KSTART:KPROF)
    PVSI(KSTART:KPROF)=ZESGP*PVSI(KSTART:KPROF)
    PUT0(KSTART:KPROF)=PUT0(KSTART:KPROF)+PUL9(KSTART:KPROF)
    PVT0(KSTART:KPROF)=PVT0(KSTART:KPROF)+PVL9(KSTART:KPROF)
    PUSI(KSTART:KPROF)=PUSI(KSTART:KPROF)-ZESGM*PUL0(KSTART:KPROF)
    PVSI(KSTART:KPROF)=PVSI(KSTART:KPROF)-ZESGM*PVL0(KSTART:KPROF)
  ENDIF
  ZGAGT0L(KSTART:KPROF)=ZBT*PUSI(KSTART:KPROF)
  ZGAGT0M(KSTART:KPROF)=ZBT*PVSI(KSTART:KPROF)

ELSE        !!! AD not done yet !!!

ENDIF

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF THE WIND COMPONENTS NECESSARY FOR SL TRAJECTORY.
!              ---------------------------------------------------------------

IF (LTWOTL) THEN

! * 2TL SL scheme.
  IF (LSETTLST) THEN
    PUT0(KSTART:KPROF)=PUT0(KSTART:KPROF)+PURL(KSTART:KPROF)
    PVT0(KSTART:KPROF)=PVT0(KSTART:KPROF)+PVRL(KSTART:KPROF)
  ELSE
    PURL0(KSTART:KPROF)=PURL0(KSTART:KPROF)+PURL(KSTART:KPROF)
    PVRL0(KSTART:KPROF)=PVRL0(KSTART:KPROF)+PVRL(KSTART:KPROF)
  ENDIF
  IF (LELTRA) THEN
    PUT0(KSTART:KPROF)=PUT0(KSTART:KPROF)+PURL0(KSTART:KPROF)&
     & -0.5_JPRB*PDT*YDGSGEOM%RCORI(KSTART:KPROF)*PVRL0(KSTART:KPROF)  
    PVT0(KSTART:KPROF)=PVT0(KSTART:KPROF)+PVRL0(KSTART:KPROF)&
     & +0.5_JPRB*PDT*YDGSGEOM%RCORI(KSTART:KPROF)*PURL0(KSTART:KPROF)  
    PSPT0L(KSTART:KPROF)=PSPT0L(KSTART:KPROF)-0.5_JPRB*PDT*PURL0(KSTART:KPROF)
    PSPT0M(KSTART:KPROF)=PSPT0M(KSTART:KPROF)-0.5_JPRB*PDT*PVRL0(KSTART:KPROF)
  ELSE
    IF (NSTEP <= 0) THEN
      PUT0(KSTART:KPROF)=PUT0(KSTART:KPROF)+PURL0(KSTART:KPROF)
      PVT0(KSTART:KPROF)=PVT0(KSTART:KPROF)+PVRL0(KSTART:KPROF)
    ELSE
      IF (LSETTLST) THEN
        PUT0(KSTART:KPROF)=PUT0(KSTART:KPROF)+2.0_JPRB*PURL0(KSTART:KPROF)
        PVT0(KSTART:KPROF)=PVT0(KSTART:KPROF)+2.0_JPRB*PVRL0(KSTART:KPROF)
        PUT9(KSTART:KPROF)=PUT9(KSTART:KPROF)-PURL0(KSTART:KPROF)
        PVT9(KSTART:KPROF)=PVT9(KSTART:KPROF)-PVRL0(KSTART:KPROF)
      ELSE
        PUT0(KSTART:KPROF)=PUT0(KSTART:KPROF)+1.5_JPRB*PURL0(KSTART:KPROF)
        PVT0(KSTART:KPROF)=PVT0(KSTART:KPROF)+1.5_JPRB*PVRL0(KSTART:KPROF)
        PUT9(KSTART:KPROF)=PUT9(KSTART:KPROF)-0.5_JPRB*PURL0(KSTART:KPROF)
        PVT9(KSTART:KPROF)=PVT9(KSTART:KPROF)-0.5_JPRB*PVRL0(KSTART:KPROF)
      ENDIF
    ENDIF
  ENDIF

ELSE       !!! AD not done yet !!!

ENDIF

!     ------------------------------------------------------------------

!*       2.    COMPUTATION OF THE LINEAR TERMS FOR SEMI-IMPLICIT SCHEME.
!              ---------------------------------------------------------

!     * Momentum equation.
!       Add semi-implicit Coriolis terms if required (LIMPF=.T.).

IF (LIMPF) THEN
  PUT0(KSTART:KPROF)=PUT0(KSTART:KPROF)&
   & +YDGSGEOM%RCORI(KSTART:KPROF)*ZGAGT0M(KSTART:KPROF)  
  PVT0(KSTART:KPROF)=PVT0(KSTART:KPROF)&
   & -YDGSGEOM%RCORI(KSTART:KPROF)*ZGAGT0L(KSTART:KPROF)  
ENDIF
PSPT0L(KSTART:KPROF)=PSPT0L(KSTART:KPROF)+ZGAGT0L(KSTART:KPROF)
PSPT0M(KSTART:KPROF)=PSPT0M(KSTART:KPROF)+ZGAGT0M(KSTART:KPROF)

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LACDYNSHWAD',1,ZHOOK_HANDLE)
END SUBROUTINE LACDYNSHWAD
