SUBROUTINE MF_PHYS_PREP(YDGEOMETRY,YDCPG_DYN0,YDCPG_DYN9,YDGMV,&
 ! --- INPUT -----------------------------------------------------------------
 & YDARPHY,KST,KEND,KIBL,POROGL,POROGM,PGMV,PR0,PR9,&
 & PGWFT0,PGWFT9,PGWFL,PGWFM,PNHPRE0F,PNHPRE9F,PNHPRE0H,PNHPRE9H,&
 & PRE0F,PRE9F,PRE0,PRE9,PXYB0,PXYB9,&
 ! --- OUTPUT ----------------------------------------------------------------
 & PW0,PW9,PW0L,PW0M,PRE0F_PHY,PRE9F_PHY,PRE0_PHY,PRE9_PHY,&
 & PREHYD0F_PHY,PREHYD9F_PHY,PREHYD0_PHY,PREHYD9_PHY,&
 & PXYB0_PHY,PXYB9_PHY)

!**** *MF_PHYS_PREP* .

!     Purpose.
!     --------
!         Does preparations for MF physics (some memory transfers for example)

!**   Interface.
!     ----------
!        *CALL* *MF_PHYS_PREP(...)*

!        Explicit arguments :
!        --------------------

!     INPUT:
!     ------
!      KST              : first element of work.
!      KEND             : last element of work.
!      KIBL             : index into YROROG types in YDGEOMETRY
!      POROGL,POROGM    : components of grad(orography).
!      PGMV             : upper air GMV variables at time t and t-dt.
!      PR0              : air constant "R" at t.
!      PR9              : air constant "R" at t-dt.
!      PGWFT0           : [gw] at full levels at t.
!      PGWFT9           : [gw] at full levels at t-dt.
!      PGWFL            : zonal comp grad(gw) at full levels at time t.
!      PGWFM            : merid comp grad(gw) at full levels at time t.
!      PNHPRE0F         : "pre" at full levels (time t).
!      PNHPRE9F         : "pre" at full levels (time t-dt).
!      PNHPRE0H         : "pre" at half levels (time t).
!      PNHPRE9H         : "pre" at half levels (time t-dt).
!      PRE0F            : "prehyd" at full levels at time t.
!      PRE9F            : "prehyd" at full levels at time t-dt.
!      PRE0             : "prehyd" at half levels at time t.
!      PRE9             : "prehyd" at half levels at time t-dt.
!      PXYB0            : contains pressure depth, "delta", "alpha" at time t.
!      PXYB9            : contains pressure depth, "delta", "alpha" at time t-dt.

!     OUTPUT:
!     -------
!      PW0              : "w" at full levels at time t.
!      PW9              : "w" at full levels at time t-dt.
!      PW0L             : zonal comp grad(w) at full levels at t.
!      PW0M             : merid comp grad(w) at full levels at t.
!      PRE0F_PHY        : input "pre" for AROME at full levels at time t.
!      PRE9F_PHY        : input "pre" for AROME at full levels at time t-dt.
!      PRE0_PHY         : input "pre" for AROME at half levels at time t.
!      PRE9_PHY         : input "pre" for AROME at half levels at time t-dt.
!      PREHYD0F_PHY     : input "prehyd" for phys. at full levels at t.
!      PREHYD9F_PHY     : input "prehyd" for phys. at full levels at t-dt.
!      PREHYD0_PHY      : input "prehyd" for phys. at half levels at t.
!      PREHYD9_PHY      : input "prehyd" for phys. at half levels at t-dt.
!      PXYB0_PHY        : contains pressure depth, "delta", "alpha" for physics input at t.
!      PXYB9_PHY        : contains pressure depth, "delta", "alpha" for physics input at t-dt.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        K. Yessad (March 2008), after some CPG pieces of code.

!     Modifications.
!     --------------
!   K. Yessad (March 2009): correct false comments for LRWSDLG=T
!   K. Yessad (Jan 2011): introduce INTDYN_MOD structures.
!   F. Vana  22-Feb_2011:  Horizontal derivatives of w for 3D TKE shear
!   R. El Khatib 16-Mar-2012 Cleaning
!   K. Yessad (June 2017): Introduce NHQE model.
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!     -------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE CPG_DYN_TYPE_MOD,ONLY : CPG_DYN_TYPE
USE YOMGMV   , ONLY : TGMV
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCT0   , ONLY : LTWOTL, LNHDYN
USE YOMCST   , ONLY : RG
USE YOMARPHY , ONLY : TARPHY
USE YOMDYNA  , ONLY : NVDVAR
USE INTDYN_MOD,ONLY : YYTXYB0_PHY,YYTXYB9_PHY,YYTXYB0  ,YYTXYB9

!     -------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_DYN_TYPE),INTENT(INOUT) :: YDCPG_DYN0
TYPE(CPG_DYN_TYPE),INTENT(INOUT) :: YDCPG_DYN9
TYPE(TGMV)        ,INTENT(INOUT) :: YDGMV
TYPE(TARPHY)      ,INTENT(INOUT) :: YDARPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIBL
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGMV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDGMV%NDIMGMV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PR0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PR9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWFT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWFT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWFL(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWFM(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHPRE0F(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHPRE9F(YDGEOMETRY%YRDIM%NPROMNH,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHPRE0H(YDGEOMETRY%YRDIM%NPROMNH,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHPRE9H(YDGEOMETRY%YRDIM%NPROMNH,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRE0F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRE9F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRE0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRE9(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXYB0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB0%NDIM) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXYB9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB9%NDIM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PW0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PW9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PW0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PW0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0F_PHY(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE9F_PHY(YDGEOMETRY%YRDIM%NPROMM9,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE0_PHY(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRE9_PHY(YDGEOMETRY%YRDIM%NPROMM9,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PREHYD0F_PHY(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PREHYD9F_PHY(YDGEOMETRY%YRDIM%NPROMM9,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PREHYD0_PHY(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PREHYD9_PHY(YDGEOMETRY%YRDIM%NPROMM9,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXYB0_PHY(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB0_PHY%NDIM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXYB9_PHY(YDGEOMETRY%YRDIM%NPROMM9,YDGEOMETRY%YRDIMV%NFLEVG,YYTXYB9_PHY%NDIM) 
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZUSG
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: JLEV, JROF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MF_PHYS_PREP',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,&
 & YDOROG=>YDGEOMETRY%YROROG(KIBL),YDVAB=>YDGEOMETRY%YRVAB)
ASSOCIATE(YT0=>YDGMV%YT0, YT9=>YDGMV%YT9, &
 & LMPA=>YDARPHY%LMPA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NPROMNH=>YDDIM%NPROMNH, NPROMA=>YDDIM%NPROMA, &
 & NPROMM9=>YDDIM%NPROMM9, NDIMGMV=>YDGMV%NDIMGMV, &
 & POROG=>YDOROG%OROG)
!     ------------------------------------------------------------------

!        1.    Input for AROME physics: transfer "w" and "pre"
!              -----------------------------------------------

! * Transfer "w" and "grad(w)"

ZUSG=1.0_JPRB/RG

PW0(KST:KEND,1:NFLEVG)=PGWFT0(KST:KEND,1:NFLEVG)*ZUSG
IF (LNHDYN .AND. NVDVAR==5) THEN
  ! compute w from W:
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      PW0(JROF,JLEV)=PW0(JROF,JLEV)+ZUSG*YDVAB%VRATF(JLEV)* &
       & (PGMV(JROF,JLEV,YDGMV%YT0%MU)*POROGL(JROF) &
       & +PGMV(JROF,JLEV,YDGMV%YT0%MV)*POROGM(JROF))
    ENDDO
  ENDDO
ENDIF

IF (LMPA.AND.(.NOT.LTWOTL)) THEN
  PW9(KST:KEND,1:NFLEVG)=PGWFT9(KST:KEND,1:NFLEVG)*ZUSG
  IF (LNHDYN .AND. (NVDVAR==5)) THEN
    ! compute w from W:
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        PW9(JROF,JLEV)=PW9(JROF,JLEV)+ZUSG*YDVAB%VRATF(JLEV)* &
         &(PGMV(JROF,JLEV,YDGMV%YT9%MU)*POROGL(JROF) &
         & +PGMV(JROF,JLEV,YDGMV%YT9%MV)*POROGM(JROF))
      ENDDO
    ENDDO
  ENDIF
ELSE
  PW9(:,:)=0.0_JPRB
ENDIF

IF (LNHDYN) THEN
  PW0L(KST:KEND,1:NFLEVG)=PGWFL(KST:KEND,1:NFLEVG)*ZUSG
  PW0M(KST:KEND,1:NFLEVG)=PGWFM(KST:KEND,1:NFLEVG)*ZUSG
  ! for NVDVAR==5, a transformation grad(gW) towards grad(gw) would be
  ! required (ignored for the time being)
ELSE
  PW0L(:,:)=0.0_JPRB
  PW0M(:,:)=0.0_JPRB
ENDIF

! * Transfer "pre"

PRE0F_PHY(:,:)=0.0_JPRB
PRE0_PHY(:,:)=0.0_JPRB
PRE9F_PHY(:,:)=0.0_JPRB
PRE9_PHY(:,:)=0.0_JPRB

IF (LMPA) THEN
  IF (LNHDYN.AND.LTWOTL) THEN
    PRE0F_PHY(KST:KEND,1:NFLEVG)=PNHPRE0F(KST:KEND,1:NFLEVG)
    PRE0_PHY(KST:KEND,0:NFLEVG)=PNHPRE0H(KST:KEND,0:NFLEVG)
  ELSEIF (LNHDYN.AND.(.NOT.LTWOTL)) THEN
    PRE9F_PHY(KST:KEND,1:NFLEVG)=PNHPRE9F(KST:KEND,1:NFLEVG)
    PRE9_PHY(KST:KEND,0:NFLEVG)=PNHPRE9H(KST:KEND,0:NFLEVG)
  ELSEIF ((.NOT.LNHDYN).AND.LTWOTL) THEN
    PRE0F_PHY(KST:KEND,1:NFLEVG)=PRE0F(KST:KEND,1:NFLEVG)
    PRE0_PHY(KST:KEND,0:NFLEVG)=PRE0(KST:KEND,0:NFLEVG)
  ELSEIF ((.NOT.LNHDYN).AND.(.NOT.LTWOTL)) THEN
    PRE9F_PHY(KST:KEND,1:NFLEVG)=PRE9F(KST:KEND,1:NFLEVG)
    PRE9_PHY(KST:KEND,0:NFLEVG)=PRE9(KST:KEND,0:NFLEVG)
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!        2.    store variables based on "prehyd" 
!              ---------------------------------

! "t" variables.
PREHYD0F_PHY(KST:KEND,1:NFLEVG)=PRE0F(KST:KEND,1:NFLEVG)
PREHYD0_PHY(KST:KEND,0:NFLEVG)=PRE0(KST:KEND,0:NFLEVG)
PXYB0_PHY(KST:KEND,1:NFLEVG,YYTXYB0_PHY%M_DELP)=PXYB0(KST:KEND,1:NFLEVG,YYTXYB0%M_DELP)
PXYB0_PHY(KST:KEND,1:NFLEVG,YYTXYB0_PHY%M_RDELP)=PXYB0(KST:KEND,1:NFLEVG,YYTXYB0%M_RDELP)
PXYB0_PHY(KST:KEND,1:NFLEVG,YYTXYB0_PHY%M_LNPR)=PXYB0(KST:KEND,1:NFLEVG,YYTXYB0%M_LNPR)
PXYB0_PHY(KST:KEND,1:NFLEVG,YYTXYB0_PHY%M_ALPH)=PXYB0(KST:KEND,1:NFLEVG,YYTXYB0%M_ALPH)

! "t-dt" variables.
IF (.NOT.LTWOTL) THEN
  PREHYD9F_PHY(KST:KEND,1:NFLEVG)=PRE9F(KST:KEND,1:NFLEVG)
  PREHYD9_PHY(KST:KEND,0:NFLEVG)=PRE9(KST:KEND,0:NFLEVG)
  PXYB9_PHY(KST:KEND,1:NFLEVG,YYTXYB9_PHY%M_DELP)=PXYB9(KST:KEND,1:NFLEVG,YYTXYB9%M_DELP)
  PXYB9_PHY(KST:KEND,1:NFLEVG,YYTXYB9_PHY%M_RDELP)=PXYB9(KST:KEND,1:NFLEVG,YYTXYB9%M_RDELP)
  PXYB9_PHY(KST:KEND,1:NFLEVG,YYTXYB9_PHY%M_LNPR)=PXYB9(KST:KEND,1:NFLEVG,YYTXYB9%M_LNPR)
  PXYB9_PHY(KST:KEND,1:NFLEVG,YYTXYB9_PHY%M_ALPH)=PXYB9(KST:KEND,1:NFLEVG,YYTXYB9%M_ALPH)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('MF_PHYS_PREP',1,ZHOOK_HANDLE)
END SUBROUTINE MF_PHYS_PREP