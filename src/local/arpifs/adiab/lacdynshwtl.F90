SUBROUTINE LACDYNSHWTL(YDGEM,YDDYN,CDCONF,KPROMA,KSTART,KPROF,PDT,&
 & PUT9,PVT9,PSPT9L,PSPT9M,&
 & PUT0,PVT0,PSPT0,PDIVT0,PSPT0L,PSPT0M,&
 & PUT5,PVT5,PSPT5,PDIVT5,PSPT5L,PSPT5M,&
 & YDGSGEOM,POROG,POROGL,POROGM,&
 & PSPNLT9,PUT95,PVT95,PSPNLT95,&
 & PUL9,PVL9,PSPL9,PUL0,PVL0,PSPL0,PSPT1,&
 & PURL0,PVRL0,PUSI,PVSI,PURL,PVRL,&
 & PUL95,PVL95,PSPL95,PUL05,PVL05,PSPL05,&
 & PUT15,PVT15,PURL05,PVRL05,PUSI5,PVSI5,PURL5,PVRL5)  

USE YOMGEM   , ONLY : TGEM
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCT0   , ONLY : LTWOTL
USE YOMCT3   , ONLY : NSTEP
USE YOMDYNA  , ONLY : YRDYNA
USE YOMDYN   , ONLY : TDYN
USE YOMGSGEOM, ONLY : TGSGEOM
USE YOMVWRK  , ONLY : NTRSLTYPE

!**** *LACDYNSHWTL* - Semi-Lagrangian scheme. - 2-D shallow-water model.
!                     Computation of the t and t-dt useful quantities
!                     at grid-points.  (Tangent-linear version)

!     Purpose.
!     --------
!          This subroutine computes the equation quantities at t-dt and t
!          at each grid-point of the colocation grid (Gauss grid).

!**   Interface.
!     ----------
!        *CALL* *LACDYNSHWTL(...)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          CDCONF  - configuration.
!          KPROMA  - horizontal dimension.
!          KSTART  - first element of work.
!          KPROF   - depth of work.
!          PDT     - For a leap-frog scheme (three time level scheme):
!                     'dt' at the first time-step, '2 dt' otherwise.
!                    For a 2TL SL scheme: timestep 'dt'.
!          PUT9    - U-component of the wind at t-dt.
!          PVT9    - V-component of the wind at t-dt.
!          PSPT9L  - zonal derivative of equivalent height "Phi" at t-dt.
!          PSPT9M  - meridian derivative of equivalent height "Phi" at t-dt.
!          PUT0    - U-component of the wind at t.
!          PVT0    - V-component of the wind at t.
!          PSPT0   - equivalent height "Phi" at t.
!          PDIVT0  - divergence at t.
!          PSPT0L  - zonal derivative of equivalent height "Phi" at t.
!          PSPT0M  - meridian derivative of equivalent height "Phi" at t.
!     ---------------------- trajectory variables -----------------
!          PUT5    - U-component of the wind at t.
!          PVT5    - V-component of the wind at t.
!          PSPT5   - equivalent height "Phi" at t.
!          PDIVT5  - divergence at t.
!          PSPT5L  - zonal derivative of equivalent height "Phi" at t.
!          PSPT5M  - meridian derivative of equivalent height "Phi" at t.
!          PUT95   - U-component of the wind at t-dt.
!          PVT95   - V-component of the wind at t-dt.
!          PSPNLT95- nonlinear term in continuity equation at t-dt.
!     -------------------------------------------------------------
!          YDGSGEOM  : Grid point geometry.
!          POROG   - grid-point surface orography
!          POROGL  - zonal component of the orography gradient
!          POROGM  - merid component of the orography gradient

!        INPUT/OUTPUT:
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
!          PUL95   - as PUL9 but for trajectory.
!          PVL95   - as PVL9 but for trajectory.
!          PSPL95  - as PSPL9 but for trajectory.
!          PUL05   - as PUL0 but for trajectory.
!          PVL05   - as PVL0 but for trajectory.
!          PSP095  - as PSPL0 but for trajectory.
!     -------------------------------------------------------------

!        OUTPUT:
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
!     ---------------------- trajectory variables -----------------
!          PUT15   - t+dt term and other final point terms for U-wind eqn.
!          PVT15   - t+dt term and other final point terms for V-wind eqn.
!          PURL05  - as PURL0 but for trajectory.
!          PVRL05  - as PVRL0 but for trajectory.
!          PUSI5   - as PUSI but for trajectory.
!          PVSI5   - as PVSI but for trajectory.
!          PURL5   - as PURL but for trajectory.
!          PVRL5   - as PVRL but for trajectory.
!     -------------------------------------------------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!      None
!      Called by CPG2TL.

!     Reference.
!     ----------
!        Arpege documentation about semi-lagrangian scheme.

!     Author.
!     -------
!      C. Temperton after LACDYNSHW.
!      Original : 98-12-14

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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVT9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT9L(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT9M(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVT0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIVT0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0L(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT0M(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVT5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIVT5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT5L(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPT5M(KPROMA) 
TYPE(TGSGEOM)     ,INTENT(IN)    :: YDGSGEOM
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPNLT9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT95(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVT95(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPNLT95(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUL9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVL9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPL9(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUL0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVL0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPL0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPT1(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PURL0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVRL0(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSI(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PURL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVRL(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUL95(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVL95(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPL95(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUL05(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVL05(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPL05(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUT15(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVT15(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PURL05(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVRL05(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUSI5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVSI5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PURL5(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PVRL5(KPROMA) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZBDT(KPROMA)
REAL(KIND=JPRB) :: ZGAGT9L(KPROMA),ZGAGT9M(KPROMA)
REAL(KIND=JPRB) :: ZGAGT0L(KPROMA),ZGAGT0M(KPROMA)
REAL(KIND=JPRB) :: ZGAGT5L(KPROMA),ZGAGT5M(KPROMA)
REAL(KIND=JPRB) :: ZMOY1SP(KPROMA)
REAL(KIND=JPRB) :: ZMOYSPNL(KPROMA)
REAL(KIND=JPRB) :: ZSPNLT0(KPROMA)
REAL(KIND=JPRB) :: ZMOY1SP5(KPROMA)
REAL(KIND=JPRB) :: ZMOYSPNL5(KPROMA)
REAL(KIND=JPRB) :: ZSPSI(KPROMA)
REAL(KIND=JPRB) :: ZSPSI5(KPROMA)

REAL(KIND=JPRB) :: ZBETA, ZBT, ZCMSLP, ZDTS2, ZESGM, ZESGP
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LACDYNSHWTL',0,ZHOOK_HANDLE)
ASSOCIATE(BETADT=>YDDYN%BETADT, LADVF=>YDDYN%LADVF, LIMPF=>YDDYN%LIMPF, &
 & LSIDG=>YDDYN%LSIDG, NVLAG=>YDDYN%NVLAG, NWLAG=>YDDYN%NWLAG, &
 & RCMSLP0=>YDDYN%RCMSLP0, SIVP=>YDDYN%SIVP, VESL=>YDDYN%VESL, &
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
  CALL ABOR1(' LACDYNSHWTL: Not existing CDCONF ')
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

!*       2.    COMPUTATION OF THE LINEAR TERMS FOR SEMI-IMPLICIT SCHEME.
!              ---------------------------------------------------------

!     * Momentum equation.
!       Add semi-implicit Coriolis terms if required (LIMPF=.T.).

IF (NTRSLTYPE <= 1) THEN
  ZGAGT5L(KSTART:KPROF)=PSPT5L(KSTART:KPROF)
  ZGAGT5M(KSTART:KPROF)=PSPT5M(KSTART:KPROF)
  IF (LIMPF) THEN
    ZGAGT5L(KSTART:KPROF)=ZGAGT5L(KSTART:KPROF)&
     & -YDGSGEOM%RCORI(KSTART:KPROF)*PVT5(KSTART:KPROF)  
    ZGAGT5M(KSTART:KPROF)=ZGAGT5M(KSTART:KPROF)&
     & +YDGSGEOM%RCORI(KSTART:KPROF)*PUT5(KSTART:KPROF)  
  ENDIF
ENDIF

ZGAGT0L(KSTART:KPROF)=PSPT0L(KSTART:KPROF)
ZGAGT0M(KSTART:KPROF)=PSPT0M(KSTART:KPROF)
IF (.NOT.LTWOTL) THEN
  ZGAGT9L(KSTART:KPROF)=PSPT9L(KSTART:KPROF)
  ZGAGT9M(KSTART:KPROF)=PSPT9M(KSTART:KPROF)
ENDIF
IF (LIMPF) THEN
  ZGAGT0L(KSTART:KPROF)=ZGAGT0L(KSTART:KPROF)&
   & -YDGSGEOM%RCORI(KSTART:KPROF)*PVT0(KSTART:KPROF)  
  ZGAGT0M(KSTART:KPROF)=ZGAGT0M(KSTART:KPROF)&
   & +YDGSGEOM%RCORI(KSTART:KPROF)*PUT0(KSTART:KPROF)  
  IF (.NOT.LTWOTL) THEN
    ZGAGT9L(KSTART:KPROF)=ZGAGT9L(KSTART:KPROF)&
     & -YDGSGEOM%RCORI(KSTART:KPROF)*PVT9(KSTART:KPROF)  
    ZGAGT9M(KSTART:KPROF)=ZGAGT9M(KSTART:KPROF)&
     & +YDGSGEOM%RCORI(KSTART:KPROF)*PUT9(KSTART:KPROF)  
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF THE WIND COMPONENTS NECESSARY FOR SL TRAJECTORY.
!              ---------------------------------------------------------------

IF (LTWOTL) THEN

! * 2TL SL scheme.
  IF (YRDYNA%LELTRA) THEN
    PURL0(KSTART:KPROF) = PUT0(KSTART:KPROF)&
     & +0.5_JPRB*PDT*(YDGSGEOM%RCORI(KSTART:KPROF)*PVT0(KSTART:KPROF)&
     & -PSPT0L(KSTART:KPROF))  
    PVRL0(KSTART:KPROF) = PVT0(KSTART:KPROF)&
     & -0.5_JPRB*PDT*(YDGSGEOM%RCORI(KSTART:KPROF)*PUT0(KSTART:KPROF)&
     & +PSPT0M(KSTART:KPROF))  
  ELSE
    IF (NSTEP <= 0) THEN
      PURL0(KSTART:KPROF)=PUT0(KSTART:KPROF)
      PVRL0(KSTART:KPROF)=PVT0(KSTART:KPROF)
    ELSE
      IF (YRDYNA%LSETTLST) THEN
        PURL0(KSTART:KPROF)=2.0_JPRB*PUT0(KSTART:KPROF)-PUT9(KSTART:KPROF)
        PVRL0(KSTART:KPROF)=2.0_JPRB*PVT0(KSTART:KPROF)-PVT9(KSTART:KPROF)
      ELSE
        PURL0(KSTART:KPROF)=1.5_JPRB*PUT0(KSTART:KPROF)&
         & -0.5_JPRB*PUT9(KSTART:KPROF)  
        PVRL0(KSTART:KPROF)=1.5_JPRB*PVT0(KSTART:KPROF)&
         & -0.5_JPRB*PVT9(KSTART:KPROF)  
      ENDIF
    ENDIF
  ENDIF
  IF (YRDYNA%LSETTLST) THEN
    PURL(KSTART:KPROF)=PUT0(KSTART:KPROF)
    PVRL(KSTART:KPROF)=PVT0(KSTART:KPROF)
  ELSE
    PURL(KSTART:KPROF)=PURL0(KSTART:KPROF)
    PVRL(KSTART:KPROF)=PVRL0(KSTART:KPROF)
  ENDIF

  IF (NTRSLTYPE <= 1) THEN
    IF (YRDYNA%LELTRA) THEN
      PURL05(KSTART:KPROF) = PUT5(KSTART:KPROF)&
       & +0.5_JPRB*PDT*(YDGSGEOM%RCORI(KSTART:KPROF)*PVT5(KSTART:KPROF)&
       & -PSPT5L(KSTART:KPROF))  
      PVRL05(KSTART:KPROF) = PVT5(KSTART:KPROF)&
       & -0.5_JPRB*PDT*(YDGSGEOM%RCORI(KSTART:KPROF)*PUT5(KSTART:KPROF)&
       & +PSPT5M(KSTART:KPROF))  
    ELSE
      IF (NSTEP <= 0) THEN
        PURL05(KSTART:KPROF)=PUT5(KSTART:KPROF)
        PVRL05(KSTART:KPROF)=PVT5(KSTART:KPROF)
      ELSE
        IF (YRDYNA%LSETTLST) THEN
          PURL05(KSTART:KPROF)=2.0_JPRB*PUT5(KSTART:KPROF)-PUT95(KSTART:KPROF)
          PVRL05(KSTART:KPROF)=2.0_JPRB*PVT5(KSTART:KPROF)-PVT95(KSTART:KPROF)
        ELSE
          PURL05(KSTART:KPROF)=1.5_JPRB*PUT5(KSTART:KPROF)&
           & -0.5_JPRB*PUT95(KSTART:KPROF)  
          PVRL05(KSTART:KPROF)=1.5_JPRB*PVT5(KSTART:KPROF)&
           & -0.5_JPRB*PVT95(KSTART:KPROF)  
        ENDIF
      ENDIF
    ENDIF
    IF (YRDYNA%LSETTLST) THEN
      PURL5(KSTART:KPROF)=PUT5(KSTART:KPROF)
      PVRL5(KSTART:KPROF)=PVT5(KSTART:KPROF)
    ELSE
      PURL5(KSTART:KPROF)=PURL05(KSTART:KPROF)
      PVRL5(KSTART:KPROF)=PVRL05(KSTART:KPROF)
    ENDIF
  ENDIF

ELSE       !!! TL not done yet !!!

ENDIF

!     ------------------------------------------------------------------

!*       4.    COMPUTATION OF THE MOMENTUM EQUATION RIGHT-HAND SIDE TERMS.
!              -----------------------------------------------------------

IF (LTWOTL) THEN

!       * Two-time-level scheme:
!         The non linear terms ZMOY1U+PUSI, ZMOY1V+PVSI are assumed to be
!         zero in this case (BETADT=1, LADVF=T or LIMPF=T), so simplifications 
!         are already made. 

  PUSI(KSTART:KPROF)=ZBT*ZGAGT0L(KSTART:KPROF)
  PVSI(KSTART:KPROF)=ZBT*ZGAGT0M(KSTART:KPROF)
  IF (NWLAG == 2) THEN
    PUL9(KSTART:KPROF)=PUL9(KSTART:KPROF)+PUT0(KSTART:KPROF)&
     & -ZESGM*PUSI(KSTART:KPROF)  
    PVL9(KSTART:KPROF)=PVL9(KSTART:KPROF)+PVT0(KSTART:KPROF)&
     & -ZESGM*PVSI(KSTART:KPROF)  
    PUSI(KSTART:KPROF)=ZESGP*PUSI(KSTART:KPROF)
    PVSI(KSTART:KPROF)=ZESGP*PVSI(KSTART:KPROF)
  ELSEIF(NWLAG == 3) THEN
    PUL9(KSTART:KPROF)=PUL9(KSTART:KPROF)+PUT0(KSTART:KPROF)
    PVL9(KSTART:KPROF)=PVL9(KSTART:KPROF)+PVT0(KSTART:KPROF)
    PUL0(KSTART:KPROF)=PUL0(KSTART:KPROF)-ZESGM*PUSI(KSTART:KPROF)
    PVL0(KSTART:KPROF)=PVL0(KSTART:KPROF)-ZESGM*PVSI(KSTART:KPROF)
    PUSI(KSTART:KPROF)=ZESGP*PUSI(KSTART:KPROF)
    PVSI(KSTART:KPROF)=ZESGP*PVSI(KSTART:KPROF)
  ENDIF

  IF (NTRSLTYPE <= 1) THEN
    PUSI5(KSTART:KPROF)=ZBT*ZGAGT5L(KSTART:KPROF)
    PVSI5(KSTART:KPROF)=ZBT*ZGAGT5M(KSTART:KPROF)
    IF (NWLAG == 2) THEN
      PUL95(KSTART:KPROF)=PUL95(KSTART:KPROF)+PUT5(KSTART:KPROF)&
       & -ZESGM*PUSI5(KSTART:KPROF)  
      PVL95(KSTART:KPROF)=PVL95(KSTART:KPROF)+PVT5(KSTART:KPROF)&
       & -ZESGM*PVSI5(KSTART:KPROF)  
    ELSEIF(NWLAG == 3) THEN
      PUL95(KSTART:KPROF)=PUL95(KSTART:KPROF)+PUT5(KSTART:KPROF)
      PVL95(KSTART:KPROF)=PVL95(KSTART:KPROF)+PVT5(KSTART:KPROF)
      PUL05(KSTART:KPROF)=PUL05(KSTART:KPROF)-ZESGM*PUSI5(KSTART:KPROF)
      PVL05(KSTART:KPROF)=PVL05(KSTART:KPROF)-ZESGM*PVSI5(KSTART:KPROF)
    ENDIF
  ENDIF

ELSE        !!! TL not done yet !!!

ENDIF

IF (NTRSLTYPE <= 1) THEN
  IF(LADVF) THEN
    PUT15(KSTART:KPROF)=PUT15(KSTART:KPROF)-YDGSGEOM%GOMVRL(KSTART:KPROF)
    PVT15(KSTART:KPROF)=PVT15(KSTART:KPROF)-YDGSGEOM%GOMVRM(KSTART:KPROF)
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       5.    COMPUTATION OF THE CONTINUITY EQUATION RIGHT-HAND SIDE TERMS.
!              -------------------------------------------------------------

IF (LTWOTL) THEN

!       * Two-time-level scheme:
!         Is coded only for conventional formulation of continuity eqn (NVLAG>0)

  ZMOY1SP(KSTART:KPROF)=-ZDTS2*(PDIVT0(KSTART:KPROF)&
   & *(PSPT5(KSTART:KPROF)-POROG(KSTART:KPROF))&
   & +PDIVT5(KSTART:KPROF)*PSPT0(KSTART:KPROF))&
   & +ZDTS2*ZCMSLP*(PUT0(KSTART:KPROF)*POROGL(KSTART:KPROF)&
   & +PVT0(KSTART:KPROF)*POROGM(KSTART:KPROF))  
  ZSPSI(KSTART:KPROF)=ZBDT(KSTART:KPROF)*PDIVT0(KSTART:KPROF)
  IF ((NSTEP <= 0).OR.YRDYNA%LELTRA) THEN
    ZMOYSPNL(KSTART:KPROF)=ZMOY1SP(KSTART:KPROF)+ZSPSI(KSTART:KPROF)
    ZSPNLT0(KSTART:KPROF)=ZMOY1SP(KSTART:KPROF)+ZSPSI(KSTART:KPROF)
  ELSE
    IF (YRDYNA%LSETTLS) THEN
      ZMOYSPNL(KSTART:KPROF)=&
       & (1.0_JPRB+ZESGM)*(ZMOY1SP(KSTART:KPROF)+ZSPSI(KSTART:KPROF))&
       & -PSPNLT9(KSTART:KPROF)  
      ZSPNLT0(KSTART:KPROF)=ZMOY1SP(KSTART:KPROF)+ZSPSI(KSTART:KPROF)
    ELSE
      ZMOYSPNL(KSTART:KPROF)=&
       & 1.5_JPRB*(ZMOY1SP(KSTART:KPROF)+ZSPSI(KSTART:KPROF))&
       & -0.5_JPRB*PSPNLT9(KSTART:KPROF)  
    ENDIF
  ENDIF
  PSPNLT9(KSTART:KPROF)=ZMOY1SP(KSTART:KPROF)+ZSPSI(KSTART:KPROF)
  IF (NVLAG == 2) THEN
    PSPT1(KSTART:KPROF)=PSPT1(KSTART:KPROF)+ZESGP*ZMOYSPNL(KSTART:KPROF)
    PSPL9(KSTART:KPROF)=PSPL9(KSTART:KPROF)+PSPT0(KSTART:KPROF)&
     & -ZESGM*ZSPSI(KSTART:KPROF)+ZESGM*ZMOYSPNL(KSTART:KPROF)  
  ELSEIF(NVLAG == 3) THEN
    IF (YRDYNA%LSETTLS) THEN
      PSPL0(KSTART:KPROF)=PSPL0(KSTART:KPROF)&
       & +ZMOYSPNL(KSTART:KPROF)-ZESGM*ZSPSI(KSTART:KPROF)  
      PSPT1(KSTART:KPROF)=PSPT1(KSTART:KPROF)+ZESGP*ZSPNLT0(KSTART:KPROF)
    ELSE
      PSPL0(KSTART:KPROF)=PSPL0(KSTART:KPROF)&
       & +ZESGM*ZMOYSPNL(KSTART:KPROF)-ZESGM*ZSPSI(KSTART:KPROF)  
      PSPT1(KSTART:KPROF)=PSPT1(KSTART:KPROF)+ZESGP*ZMOYSPNL(KSTART:KPROF)
    ENDIF
    PSPL9(KSTART:KPROF)=PSPL9(KSTART:KPROF)+PSPT0(KSTART:KPROF)
  ENDIF

  IF (NTRSLTYPE <= 1) THEN
    ZMOY1SP5(KSTART:KPROF)=-ZDTS2*PDIVT5(KSTART:KPROF)&
     & *(PSPT5(KSTART:KPROF)-POROG(KSTART:KPROF))&
     & +ZDTS2*ZCMSLP*(PUT5(KSTART:KPROF)*POROGL(KSTART:KPROF)&
     & +PVT5(KSTART:KPROF)*POROGM(KSTART:KPROF))  
    ZSPSI5(KSTART:KPROF)=ZBDT(KSTART:KPROF)*PDIVT5(KSTART:KPROF)
    IF ((NSTEP <= 0).OR.YRDYNA%LELTRA) THEN
      ZMOYSPNL5(KSTART:KPROF)=ZMOY1SP5(KSTART:KPROF)+ZSPSI5(KSTART:KPROF)
    ELSE
      IF (YRDYNA%LSETTLS) THEN
        ZMOYSPNL5(KSTART:KPROF)=&
         & (1.0_JPRB+ZESGM)*(ZMOY1SP5(KSTART:KPROF)+ZSPSI5(KSTART:KPROF))&
         & -PSPNLT95(KSTART:KPROF)  
      ELSE
        ZMOYSPNL5(KSTART:KPROF)=&
         & 1.5_JPRB*(ZMOY1SP5(KSTART:KPROF)+ZSPSI5(KSTART:KPROF))&
         & -0.5_JPRB*PSPNLT95(KSTART:KPROF)  
      ENDIF
    ENDIF
    IF (NVLAG == 2) THEN
      PSPL95(KSTART:KPROF)=PSPL95(KSTART:KPROF)&
       & +(PSPT5(KSTART:KPROF)-(1.0_JPRB-ZCMSLP)*POROG(KSTART:KPROF))&
       & -ZESGM*ZSPSI5(KSTART:KPROF)+ZESGM*ZMOYSPNL5(KSTART:KPROF)  
    ELSEIF(NVLAG == 3) THEN
      IF (YRDYNA%LSETTLS) THEN
        PSPL05(KSTART:KPROF)=PSPL05(KSTART:KPROF)&
         & +ZMOYSPNL5(KSTART:KPROF)-ZESGM*ZSPSI5(KSTART:KPROF)  
      ELSE
        PSPL05(KSTART:KPROF)=PSPL05(KSTART:KPROF)&
         & +ZESGM*ZMOYSPNL5(KSTART:KPROF)-ZESGM*ZSPSI5(KSTART:KPROF)  
      ENDIF
      PSPL95(KSTART:KPROF)=PSPL95(KSTART:KPROF)&
       & +(PSPT5(KSTART:KPROF)-(1.0_JPRB-ZCMSLP)*POROG(KSTART:KPROF))  
    ENDIF
  ENDIF

ELSE    !!! TL not done yet !!!

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LACDYNSHWTL',1,ZHOOK_HANDLE)
END SUBROUTINE LACDYNSHWTL
