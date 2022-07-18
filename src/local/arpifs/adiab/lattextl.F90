SUBROUTINE LATTEXTL(YDGEOMETRY,YDGMV,YDLDDH,YDML_GCONF,YDML_DYN,KST,KPROF,PDTS2,PBT,PBDT,PESGP,PESGM,&
 & KIBL,POROGL,POROGM,&
 & PGAGT0L,PGAGT0M,PTOD0,PRDELP,PEVEL,&
 & PGFL,PATND,PGMV,PGMVT1,PB1,PB2,&
 & PRDELP5,PEVEL5)

!**** *LATTEXTL* Semi-Lagrangian scheme.  (tangent-linear version)
!                Computation of the t and t-dt useful quantities
!                 at grid-points. Equations for tri-dimensional
!                 variables.

!     Purpose.
!     --------
!       TL of LATTEX: see LATTEX for more information.
!       Remarks:
!        - TL of SL3TL not yet coded.
!        - TL of SL2TL: only K[X]LAG = 3 is coded.
!        - LATTEX_DNT is its own tangent linear code.

!**   Interface.
!     ----------
!        *CALL* *LATTEXTL(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST     - first element of work.
!          KPROF   - depth of work.
!          PDTS2   - 0.5*time step for the first time-integration step of
!                    a leap-frog scheme or all time-integration steps of
!                    a two-time level scheme; time step for the following
!                    time-integration steps of a leap-frog scheme.
!          PBT     - PDTS2*BETADT (BETADT is in YOMDYN) or zero 
!                    according to the value of configuration.
!          PBDT    - PBT if semi-implicit scheme with unreduced
!                    divergence, PBT*(c**2/GM**2) if semi-implicit
!                    scheme with reduced divergence.
!          PESGP   - (1 + uncentering factor).
!          PESGM   - (1 - uncentering factor).
!          KIBL    - index into YRGSGEOM instance in YDGEOMETRY
!          POROGL  - zonal component of the surface orography gradient.
!          POROGM  - meridian component of the surface orography gradient.
!          PGAGT0L - semi-implicit term at time t for U-wind equation.
!          PGAGT0M - semi-implicit term at time t for V-wind equation.
!          PTOD0   - semi-implicit term at time t for temperature equation.
!          PRDELP  - "1/(pressure depth of layers)" at t.
!          PEVEL   - "etadot (d prehyd/d eta)" at half levels at t.
!          PGFL    - GFL at full levels at t.
!          PATND   - adiabatic Lagrangian tendencies.

!        INPUT/OUTPUT:
!          PGMV    - GMV variables at t-dt and t.
!          PGMVT1  - GMV variables at t+dt.
!          PB1     - "SLBUF1" buffer for interpolations.
!          PB2     - "SLBUF2" buffer.

!        TRAJECTORY INPUT:
!          PRDELP5 - "1/(pressure depth of layers)" at t.
!          PEVEL5  - "etadot (d prehyd/d eta)" at half levels at t.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           none
!           Called by LACDYNTL.

!     Reference.
!     ----------
!             Arpege documentation about semi-Lagrangian scheme.

!     Author.
!     -------
!      C. Temperton (ECMWF) 
!      Original : 99/08/10.

!     Modifications.
!     --------------
!      A.Untch  June-2005: TL of Rayleigh friction
!      K.Yessad  Mar-2006: minor reorganisation to match better with LATTEX.
!      K.Yessad  Mar-2007: make the code more modular, cf. the direct code.
!      K. Yessad Aug 2008: rationalisation of dummy argument interfaces
!      K. Yessad Nov 2008: rationalisation of dummy argument interfaces
!      P. Bechtold+A. Untch 26-10-2008: add LEGWWMS switch for non-orogr. GWD
!      K.Yessad (Aug 2009): remove NTRSLTYPE/=2 cases
!      F. Vana  15-Oct-2009: update of dummy args with respect to NSPLTHOI
!      K. Yessad (Nov 2009): cleanings, DT/Dt now pre-computed in CPG_GP_TL.
!      K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!      K. Yessad (Dec 2011): various contributions.
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE MODEL_DYNAMICS_MOD     , ONLY : MODEL_DYNAMICS_TYPE
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMLDDH                , ONLY : TLDDH
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE YOMGMV                 , ONLY : TGMV
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK
USE YOMCST                 , ONLY : RD
USE YOMSTA                 , ONLY : RTSUR
USE YOMCT0                 , ONLY : LTWOTL
USE YOMCVER                , ONLY : LVERTFE

USE INTDYN_MOD             , ONLY : YYTTND

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TLDDH)       ,INTENT(IN)    :: YDLDDH
TYPE(MODEL_DYNAMICS_TYPE),INTENT(IN):: YDML_DYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDTS2 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBDT(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESGP 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PESGM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAGT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGAGT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTOD0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVEL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NDIM) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PATND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTTND%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMVT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%YT1%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB1%NFLDSLB1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB2(YDGEOMETRY%YRDIM%NPROMA,YDML_DYN%YRPTRSLB2%NFLDSLB2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP5(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVEL5(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZMOY1U(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZMOY1V(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZMOY1T(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZWT0(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZWT5(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZXT9_DUM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZCXNLT9_DUM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZUSELESS(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
REAL(KIND=JPRB) :: ZDUM_SIDDH(1,1) ! not used in this call to LATTEX_DNT
REAL(KIND=JPRB) :: ZMIXNL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IPX
INTEGER(KIND=JPIM) :: JLEV,JGFL,JROF
LOGICAL :: LLSETTLSW

REAL(KIND=JPRB) :: ZCMSLP
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "lattex_dnt.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LATTEXTL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDSTA=>YDGEOMETRY%YRSTA, &
 & YDGSGEOM=>YDGEOMETRY%YRGSGEOM(KIBL), YDOROG=>YDGEOMETRY%YROROG(KIBL), YDDYN=>YDML_DYN%YRDYN, &
 & YDPTRSLB1=>YDML_DYN%YRPTRSLB1, &
 & YDPTRSLB2=>YDML_DYN%YRPTRSLB2,YGFL=>YDML_GCONF%YGFL)

ASSOCIATE(NDIM=>YGFL%NDIM, NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFLEN=>YDDIMV%NFLEN, NFLEVG=>YDDIMV%NFLEVG, NFLSA=>YDDIMV%NFLSA, &
 & LADVF=>YDDYN%LADVF, NTLAG=>YDDYN%NTLAG, NWLAG=>YDDYN%NWLAG, &
 & RCMSLP0=>YDDYN%RCMSLP0, RCORDIF=>YDDYN%RCORDIF, RCORDIH=>YDDYN%RCORDIH, &
 & NDIMGMV=>YDGMV%NDIMGMV, YT0=>YDGMV%YT0, YT1=>YDGMV%YT1, YT9=>YDGMV%YT9, &
 & MSLB1GFL9=>YDPTRSLB1%MSLB1GFL9, MSLB1T0=>YDPTRSLB1%MSLB1T0, &
 & MSLB1T9=>YDPTRSLB1%MSLB1T9, MSLB1U0=>YDPTRSLB1%MSLB1U0, &
 & MSLB1U9=>YDPTRSLB1%MSLB1U9, MSLB1V0=>YDPTRSLB1%MSLB1V0, &
 & MSLB1V9=>YDPTRSLB1%MSLB1V9, NFLDSLB1=>YDPTRSLB1%NFLDSLB1, &
 & MSLB2TSI=>YDPTRSLB2%MSLB2TSI, MSLB2USI=>YDPTRSLB2%MSLB2USI, &
 & MSLB2VSI=>YDPTRSLB2%MSLB2VSI, NFLDSLB2=>YDPTRSLB2%NFLDSLB2 )
!     ------------------------------------------------------------------

!*      1.  PRELIMINARY INITIALISATIONS.
!       --------------------------------

ZCMSLP=RCMSLP0/(RD*RTSUR)
ZMIXNL=1.0_JPRB

!     ------------------------------------------------------------------

!*      2.  TREATMENT OF GMV VARIABLES.
!       -------------------------------

!*       2.1   Momentum equation.

! * LSETTLS is replaced by LLSETTLSW=.FALSE. for wind-eqn if VESL>0 because
!   stable extrapolation deteriorates scores without improving stability.
LLSETTLSW=YDML_DYN%YRDYNA%LSETTLS.AND.(.NOT.(PESGP > PESGM))

DO JLEV=1,NFLEVG

  DO JROF=KST,KPROF
    PB2(JROF,MSLB2USI+JLEV-1)=PBT*PGAGT0L(JROF,JLEV)
    PB2(JROF,MSLB2VSI+JLEV-1)=PBT*PGAGT0M(JROF,JLEV)
  ENDDO

  ! * Pressure gradient term and Rayleigh friction.
  DO JROF=KST,KPROF
    ZMOY1U(JROF,JLEV)=PDTS2*PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)
    ZMOY1V(JROF,JLEV)=PDTS2*PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)
  ENDDO

  ! * Explicit coriolis term.
  IF (.NOT.LADVF) THEN
    DO JROF=KST,KPROF
      ZMOY1U(JROF,JLEV)=ZMOY1U(JROF,JLEV)+PDTS2*YDGSGEOM%RCORI(JROF)*PGMV(JROF,JLEV,YT0%MV)
      ZMOY1V(JROF,JLEV)=ZMOY1V(JROF,JLEV)-PDTS2*YDGSGEOM%RCORI(JROF)*PGMV(JROF,JLEV,YT0%MU)
    ENDDO
  ENDIF

ENDDO

IF (LTWOTL) THEN

  CALL LATTEX_DNT(YDGEOMETRY,YDLDDH,YDML_GCONF%YRRIP,YDDYN,KST,KPROF,LLSETTLSW,NWLAG,PESGP,PESGM,&
   & PGMV(1,1,YT0%MU),ZXT9_DUM,ZMOY1U,ZMIXNL,&
   & PB2(1,MSLB2USI),PGMV(1,1,YT9%MUNL),PGMVT1(1,1,YT1%MU),&
   & PB1(1,MSLB1U0),PB1(1,MSLB1U9),ZUSELESS,ZCXNLT9_DUM,&
   & ZDUM_SIDDH,ZDUM_SIDDH,ZDUM_SIDDH,ZDUM_SIDDH)

  CALL LATTEX_DNT(YDGEOMETRY,YDLDDH,YDML_GCONF%YRRIP,YDDYN,KST,KPROF,LLSETTLSW,NWLAG,PESGP,PESGM,&
   & PGMV(1,1,YT0%MV),ZXT9_DUM,ZMOY1V,ZMIXNL,&
   & PB2(1,MSLB2VSI),PGMV(1,1,YT9%MVNL),PGMVT1(1,1,YT1%MV),&
   & PB1(1,MSLB1V0),PB1(1,MSLB1V9),ZUSELESS,ZCXNLT9_DUM,&
   & ZDUM_SIDDH,ZDUM_SIDDH,ZDUM_SIDDH,ZDUM_SIDDH)

ELSE         !!! T/L OF THREE-TIME-LEVEL SCHEME NOT DONE !!!

  ! not coded (call to LATTEX_TNT)

ENDIF

!*       2.2   Temperature equation.

DO JLEV=1,NFLEVG
  DO JROF=KST,KPROF
    PB2(JROF,MSLB2TSI+JLEV-1)=PBDT(JROF)*PTOD0(JROF,JLEV)
  ENDDO
  IF(LVERTFE) THEN
    DO JROF=KST,KPROF
      ZWT5(JROF)=PEVEL5(JROF,JLEV)
      ZWT0(JROF)=PEVEL(JROF,JLEV)
    ENDDO
  ELSE
    DO JROF=KST,KPROF
      ZWT5(JROF)=0.5_JPRB*(PEVEL5(JROF,JLEV)+PEVEL5(JROF,JLEV-1))  
      ZWT0(JROF)=0.5_JPRB*(PEVEL(JROF,JLEV)+PEVEL(JROF,JLEV-1))  
    ENDDO
  ENDIF
  DO JROF=KST,KPROF
    ZMOY1T(JROF,JLEV)=&
     & PDTS2*PATND(JROF,JLEV,YYTTND%M_TNDT)&
     & +PDTS2*ZCMSLP*RCORDIF(JLEV)*PGMV(JROF,JLEV,YT0%MU)*POROGL(JROF)&
     & +PDTS2*ZCMSLP*RCORDIF(JLEV)*PGMV(JROF,JLEV,YT0%MV)*POROGM(JROF)&
     & +PDTS2*ZCMSLP*(RCORDIH(JLEV)-RCORDIH(JLEV-1))&
     & *( ZWT0(JROF)*PRDELP5(JROF,JLEV)+ZWT5(JROF)*PRDELP (JROF,JLEV) )&
     & *YDOROG%OROG(JROF)  
  ENDDO
ENDDO

IF (LTWOTL) THEN

  CALL LATTEX_DNT(YDGEOMETRY,YDLDDH,YDML_GCONF%YRRIP,YDDYN,KST,KPROF,YDML_DYN%YRDYNA%LSETTLS,NTLAG,PESGP,PESGM,&
   & PGMV(1,1,YT0%MT),ZXT9_DUM,ZMOY1T,ZMIXNL,&
   & PB2(1,MSLB2TSI),PGMV(1,1,YT9%MTNL),PGMVT1(1,1,YT1%MT),&
   & PB1(1,MSLB1T0),PB1(1,MSLB1T9),ZUSELESS,ZCXNLT9_DUM,&
   & ZDUM_SIDDH,ZDUM_SIDDH,ZDUM_SIDDH,ZDUM_SIDDH)

ELSE         !!! T/L OF THREE-TIME-LEVEL SCHEME NOT DONE !!!

  ! not coded (call to LATTEX_TNT)

ENDIF

!*       2.3   Anhydrostatic variables equations: "pressure departure".

!NYC-NH IF (LNHDYN) THEN
!NYC-NH   TL not yet coded.
!NYC-NH ENDIF

!*       2.4a  Anhydrostatic variables equations: "vertical divergence".

!NYC-NH IF (LNHDYN.AND.(.NOT.LGWADV)) THEN
!NYC-NH   TL not yet coded.
!NYC-NH ENDIF

!*       2.4b  Anhydrostatic variables equations: "gw".

!NYC-NH IF (LNHDYN.AND.LGWADV) THEN
!NYC-NH   TL not yet coded.
!NYC-NH ENDIF

!     ------------------------------------------------------------------

!*      3.  TREATMENT OF GFL VARIABLES.
!       --------------------------------

IF (LTWOTL) THEN

  DO JGFL=1,NUMFLDS
    IF(YCOMP(JGFL)%LADV) THEN
      IPX=(YCOMP(JGFL)%MP_SL1-1)*(NFLEN-NFLSA+1)
      DO JLEV=1,NFLEVG
        DO JROF=KST,KPROF
          PB1(JROF,MSLB1GFL9+IPX+JLEV-NFLSA) =&
           & PB1(JROF,MSLB1GFL9+IPX+JLEV-NFLSA)&
           & +PGFL(JROF,JLEV,YCOMP(JGFL)%MP)
        ENDDO
      ENDDO
    ENDIF
  ENDDO

ELSE    !!! T/L OF THREE-TIME-LEVEL SCHEME NOT DONE !!!

ENDIF

!     ------------------------------------------------------------------

!*      4.  DEALLOCATIONS.
!       ------------------

! Empty for the time being.

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LATTEXTL',1,ZHOOK_HANDLE)
END SUBROUTINE LATTEXTL

