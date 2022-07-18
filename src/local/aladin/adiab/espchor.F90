SUBROUTINE ESPCHOR(YDGEOMETRY,YDMODEL,KSTA,KEND,KOFF,&
 & PSPVORG,PSPDIVG,PSPTG,PSPSPDG,PSPSVDG,PSPSNHXG,PSPGFL,PSPSPG,&
 & PSPTNDHD_VOR,PSPTNDHD_DIV,PSPTNDHD_T,&
 & PSPTNDHD_SPD,PSPTNDHD_SVD,PSPTNDHD_SNHX,PSPTNDHD_GFL)

!**** *ESPCHOR* - SPECTRAL SPACE COMPUTATIONS FOR ALADIN
!                 HORIZONTAL DIFFUSION

!     Purpose.
!     --------
!       Compute the celebrated spectral part of the equations:
!        horizontal diffusion

!**   Interface.
!     ----------
!        *CALL* *ESPCHOR(...)

!        Explicit arguments :
!        --------------------
!         INPUT:
!          KSTA         : First column processed
!          KEND         : Last columns processed
!          KOFF         : Offset for arrays YDLEP%MVALUE and NVALUE

!         INOUT:
!          PSPVORG      : Vorticity columns
!          PSPDIVG      : Divergence columns
!          PSPTG        : Temperature columns
!          PSPSPDG      : NH pressure departure var.
!          PSPSVDG      : NH vertical divergence var.
!          PSPSNHXG     : NH "X" part divergence var.
!          PSPGFL       : GFL
!          PSPSPG       : Surface Pressure
!          PSPTNDHD_[X] : HD tendencies of variable [X] (for DDH)

!     Method.
!     -------
!       There is a great description in ESPCM

!     Externals.
!     ----------
!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Sylvie Malardel (original, March 1992)
!        Karim YESSAD (Feb 2005, moving the horiz. diffusion part in ESPCHOR)

!     Modifications.
!     --------------
!        Modified : 07-10-03 R. Radu optional spectral nudging 
!        Modified : 07-12-18 M. Deque isotropic spectral nudging
!        K. Yessad (Sep 2008): prune enhanced diffusion (lfrein).
!        F. Vana  13-Jan-2009: better treatment for LSLHD_STATIC
!        F. Vana  13-Oct-2009: Optimization (improved vect., OpenMP)
!        Y. Seity 6-Fev-2010: Optimisations of spectral coupling part
!        K. Yessad (Jan 2011): new architecture for LBC modules and set-up.
!        R. El Khatib and F. Voitus (Mar 2011): optimisations.
!        K. Yessad (Feb 2012): new sponge; tests on CDCONF in the caller.
!        A.bogatchev 11-04-2013 - Phasing cy40, coherence with modified modules
!                                 and renamed namelista and functions
!        B. Bochenek (Apr 2015): Phasing: move some variables.
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     ------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCT0   , ONLY : LNHDYN


!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTA
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
INTEGER(KIND=JPIM),INTENT(IN)    :: KOFF
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVORG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIVG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPDG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSVDG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSNHXG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPGFL(:,:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPG(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDHD_VOR(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDHD_DIV(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDHD_T(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDHD_SPD(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDHD_SVD(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDHD_SNHX(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDHD_GFL(:,:,:)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IM, IN, JLEV, JSP, JGFL
INTEGER(KIND=JPIM) :: ISP_VEC(KSTA:KEND)

REAL(KIND=JPRB) :: ZCOEF(KSTA:KEND)
REAL(KIND=JPRB) :: ZDT
LOGICAL :: LLDO

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ESPCHOR',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIMV=>YDGEOMETRY%YRDIMV, YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH, YDMDDH=>YDMODEL%YRML_DIAG%YRMDDH, &
 & YDML_DYN=>YDMODEL%YRML_DYN, YDML_GCONF=>YDMODEL%YRML_GCONF, &
 & YDLAP=>YDGEOMETRY%YRLAP, YDLEP=>YDGEOMETRY%YRELAP)
ASSOCIATE( &
 & LSPT=>YDML_GCONF%YRDIMF%LSPT, RCORDIT=>YDML_DYN%YRDYN%RCORDIT, TDT=>YDML_GCONF%YRRIP%TDT, &
 & YCOMP=>YDML_GCONF%YGFL%YCOMP, YQ=>YDML_GCONF%YGFL%YQ, YL=>YDML_GCONF%YGFL%YL, YI=>YDML_GCONF%YGFL%YI, &
 & YA=>YDML_GCONF%YGFL%YA, YO3=>YDML_GCONF%YGFL%YO3, &
 & RDIVORE=>YDML_DYN%YREDYN%RDIVORE, RDIDIVE=>YDML_DYN%YREDYN%RDIDIVE, RDITE=>YDML_DYN%YREDYN%RDITE, &
 & RDIGFLE=>YDML_DYN%YREDYN%RDIGFLE, RDIPDE=>YDML_DYN%YREDYN%RDIPDE, RDIVDE=>YDML_DYN%YREDYN%RDIVDE, &
 & RDISPE=>YDML_DYN%YREDYN%RDISPE, RDSVORE=>YDML_DYN%YREDYN%RDSVORE, RDSDIVE=>YDML_DYN%YREDYN%RDSDIVE, &
 & NFLEVG=>YDDIMV%NFLEVG, RDSVDE=>YDML_DYN%YREDYN%RDSVDE, LRHDDDH=>YDLDDH%LRHDDDH, &
 & MHDDDH_Q=>YDMDDH%MHDDDH_Q, LNSPONGE=>YDML_DYN%YRSPNG%LNSPONGE, RSPONGF=>YDML_DYN%YRSPNG%RSPONGF )
!     ------------------------------------------------------------------

!*       1.    Horizontal diffusion
!              --------------------

DO JSP=KSTA,KEND
  IN=YDLAP%NVALUE(JSP+KOFF)
  IM=YDLEP%MVALUE(JSP+KOFF)
  ISP_VEC(JSP)=YDLEP%NPME(IM)+IN
  ! ZCOEF=0. for (m,n)=(0,0), 1. otherwise.
  ZCOEF(JSP)=MAX(MIN(REAL(IN,JPRB)+REAL(IM,JPRB),1._JPRB),0._JPRB)
ENDDO

ZDT=ABS(TDT)
CALL GSTATS(1031,0)

! Memory transfert for spectral HD tendencies
IF (LRHDDDH) THEN
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      PSPTNDHD_VOR(JLEV,JSP)=-PSPVORG(JLEV,JSP)
      PSPTNDHD_DIV(JLEV,JSP)=-PSPDIVG(JLEV,JSP)
      PSPTNDHD_T(JLEV,JSP)=-PSPTG(JLEV,JSP)
      IF (LNHDYN) THEN
        PSPTNDHD_SPD(JLEV,JSP)=-PSPSPDG(JLEV,JSP)
        PSPTNDHD_SVD(JLEV,JSP)=-PSPSVDG(JLEV,JSP)
        IF (YDMODEL%YRML_DYN%YRDYNA%LNHX) PSPTNDHD_SNHX(JLEV,JSP)=-PSPSNHXG(JLEV,JSP)
      ENDIF
    ENDDO
  ENDDO
  IF (YQ%LSP) THEN
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
        PSPTNDHD_GFL(JLEV,JSP,MHDDDH_Q)=-PSPGFL(JLEV,JSP,YQ%MPSP)
      ENDDO
    ENDDO
  ENDIF
ENDIF

!           1.1  Horizontal diffusion.

! * VOR,DIV:
DO JSP=KSTA,KEND
  DO JLEV=1,NFLEVG
!cdir nodep
    PSPVORG(JLEV,JSP)=PSPVORG(JLEV,JSP)&
     & /(1.0_JPRB+ZDT*RDIVORE(JLEV,ISP_VEC(JSP)))
    PSPDIVG(JLEV,JSP)=PSPDIVG(JLEV,JSP)&
     & /(1.0_JPRB+ZDT*RDIDIVE(JLEV,ISP_VEC(JSP)))
  ENDDO
ENDDO
! * Additional diffusion
IF (YDMODEL%YRML_DYN%YRDYNA%LSLHD_W.AND.(.NOT.YDMODEL%YRML_DYN%YRDYNA%LSLHD_STATIC)) THEN
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
!cdir nodep
      PSPVORG(JLEV,JSP)=PSPVORG(JLEV,JSP)&
       & /(1.0_JPRB+ZDT*RDSVORE(JLEV,ISP_VEC(JSP)))
      PSPDIVG(JLEV,JSP)=PSPDIVG(JLEV,JSP)&
       & /(1.0_JPRB+ZDT*RDSDIVE(JLEV,ISP_VEC(JSP)))
    ENDDO
  ENDDO
ENDIF
! * Temperature:
IF (LSPT) THEN
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
!cdir nodep
      PSPTG  (JLEV,JSP)=(PSPTG  (JLEV,JSP)&
       & +RCORDIT(JLEV)*PSPSPG(JSP)*ZDT*RDITE(JLEV,ISP_VEC(JSP)))&
       & /(1.0_JPRB+ZDT*RDITE(JLEV,ISP_VEC(JSP)))
    ENDDO
  ENDDO
ENDIF
! * NH variables:
IF(LNHDYN) THEN
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
!cdir nodep
      PSPSPDG(JLEV,JSP)=PSPSPDG(JLEV,JSP)&
       & /(1.0_JPRB+ZDT*RDIPDE(JLEV,ISP_VEC(JSP)))
      PSPSVDG(JLEV,JSP)=PSPSVDG(JLEV,JSP)&
       & /(1.0_JPRB+ZDT*RDIVDE(JLEV,ISP_VEC(JSP)))
    ENDDO
  ENDDO
  IF(YDMODEL%YRML_DYN%YRDYNA%LNHX)THEN
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
!cdir nodep
        PSPSNHXG(JLEV,JSP)= PSPSNHXG(JLEV,JSP)&
         & /(1.0_JPRB+ZDT*RDIVDE(JLEV,ISP_VEC(JSP)))
      ENDDO
    ENDDO
  ENDIF
  ! * Additional diffusion of VD done in the direct integration only.
  IF (YDMODEL%YRML_DYN%YRDYNA%LSLHD_SVD.AND.(.NOT.YDMODEL%YRML_DYN%YRDYNA%LSLHD_STATIC)) THEN
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
!cdir nodep
        PSPSVDG(JLEV,JSP)=PSPSVDG(JLEV,JSP)&
         & /(1.0_JPRB+ZDT*RDSVDE(JLEV,ISP_VEC(JSP)))
      ENDDO
    ENDDO
  ENDIF
  ! * Additional diffusion of aux VD done in the direct integration only.
  IF (YDMODEL%YRML_DYN%YRDYNA%LNHX.AND.YDMODEL%YRML_DYN%YRDYNA%LSLHD_SVD.AND.(.NOT.YDMODEL%YRML_DYN%YRDYNA%LSLHD_STATIC)) THEN
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
!cdir nodep
        PSPSNHXG(JLEV,JSP)= PSPSNHXG(JLEV,JSP)&
         & /(1.0_JPRB+ZDT*RDSVDE(JLEV,ISP_VEC(JSP)))
      ENDDO
    ENDDO
  ENDIF
ENDIF

! * GFL variables:
DO JGFL=1,YDML_GCONF%YGFL%NUMFLDS
  IF(YCOMP(JGFL)%LSP) THEN
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
        PSPGFL(JLEV,JSP,YCOMP(JGFL)%MPSP) =&
         & PSPGFL(JLEV,JSP,YCOMP(JGFL)%MPSP)/&
         & (1.0_JPRB+ZDT*RDIGFLE(JLEV,ISP_VEC(JSP),YCOMP(JGFL)%MPSP))
      ENDDO
    ENDDO
  ENDIF
ENDDO

! * Log(prehyds):
!cdir nodep
DO JSP=KSTA,KEND
  PSPSPG (JSP)=PSPSPG (JSP)/(1.0_JPRB+ZDT*RDISPE(ISP_VEC(JSP)))
ENDDO

CALL GSTATS(1031,1)

!           1.2  Sponge (= additional diffusion near the model top).

IF (LNSPONGE) THEN

  ! * 1.2.1: Sponge for GMV.
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
!cdir nodep
      PSPVORG(JLEV,JSP)=((1.0_JPRB-ZCOEF(JSP))+ZCOEF(JSP)*RSPONGF(JLEV,1))*PSPVORG(JLEV,JSP)
      PSPDIVG(JLEV,JSP)=((1.0_JPRB-ZCOEF(JSP))+ZCOEF(JSP)*RSPONGF(JLEV,1))*PSPDIVG(JLEV,JSP)
    ENDDO
    IF (LSPT) THEN
      DO JLEV=1,NFLEVG
!cdir nodep
        PSPTG(JLEV,JSP)=((1.0_JPRB-ZCOEF(JSP))+ZCOEF(JSP)*RSPONGF(JLEV,1))*PSPTG(JLEV,JSP)
      ENDDO
    ENDIF
    IF (LNHDYN) THEN
      DO JLEV=1,NFLEVG
!cdir nodep
        PSPSPDG(JLEV,JSP)=((1.0_JPRB-ZCOEF(JSP))+ZCOEF(JSP)*RSPONGF(JLEV,1))*PSPSPDG(JLEV,JSP)
        PSPSVDG(JLEV,JSP)=((1.0_JPRB-ZCOEF(JSP))+ZCOEF(JSP)*RSPONGF(JLEV,1))*PSPSVDG(JLEV,JSP)
      ENDDO
    ENDIF
    IF (LNHDYN.AND.YDMODEL%YRML_DYN%YRDYNA%LNHX) THEN
      DO JLEV=1,NFLEVG
!cdir nodep
        PSPSNHXG(JLEV,JSP)=((1.0_JPRB-ZCOEF(JSP))+ZCOEF(JSP)*RSPONGF(JLEV,1))*PSPSNHXG(JLEV,JSP)
      ENDDO
    ENDIF
  ENDDO

  ! * 1.2.2: Sponge for GFL.
  DO JGFL=1,YDML_GCONF%YGFL%NUMFLDS
    IF(YCOMP(JGFL)%LSP) THEN
      ! sponge is applied to q,ql,qi,qa,O3 only.
      LLDO=( (YCOMP(JGFL)%MPSP == YQ%MPSP) .OR.&
       & (YCOMP(JGFL)%MPSP == YL%MPSP) .OR. (YCOMP(JGFL)%MPSP == YI%MPSP) .OR.&
       & (YCOMP(JGFL)%MPSP == YA%MPSP) .OR. (YCOMP(JGFL)%MPSP == YO3%MPSP) )
      IF (LLDO) THEN
        DO JSP=KSTA,KEND
          DO JLEV=1,NFLEVG
!cdir nodep
            PSPGFL(JLEV,JSP,YCOMP(JGFL)%MPSP) =&
             & ((1.0_JPRB-ZCOEF(JSP))+ZCOEF(JSP)*RSPONGF(JLEV,2))&
             & * PSPGFL(JLEV,JSP,YCOMP(JGFL)%MPSP)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO

ENDIF ! LNSPONGE

! Memory transfert for spectral HD tendencies
IF (LRHDDDH) THEN
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      PSPTNDHD_VOR(JLEV,JSP)=PSPTNDHD_VOR(JLEV,JSP)+PSPVORG(JLEV,JSP)
      PSPTNDHD_DIV(JLEV,JSP)=PSPTNDHD_DIV(JLEV,JSP)+PSPDIVG(JLEV,JSP)
      PSPTNDHD_T(JLEV,JSP)=PSPTNDHD_T(JLEV,JSP)+PSPTG(JLEV,JSP)
      IF (LNHDYN) THEN
        PSPTNDHD_SPD(JLEV,JSP)=PSPTNDHD_SPD(JLEV,JSP)+PSPSPDG(JLEV,JSP)
        PSPTNDHD_SVD(JLEV,JSP)=PSPTNDHD_SVD(JLEV,JSP)+PSPSVDG(JLEV,JSP)
        IF (YDMODEL%YRML_DYN%YRDYNA%LNHX) PSPTNDHD_SNHX(JLEV,JSP)=PSPTNDHD_SNHX(JLEV,JSP)+PSPSNHXG(JLEV,JSP)
      ENDIF
    ENDDO
  ENDDO
  IF (YQ%LSP) THEN
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
        PSPTNDHD_GFL(JLEV,JSP,MHDDDH_Q)=PSPTNDHD_GFL(JLEV,JSP,1)&
         & + PSPGFL(JLEV,JSP,YQ%MPSP)
      ENDDO
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ESPCHOR',1,ZHOOK_HANDLE)
END SUBROUTINE ESPCHOR
