SUBROUTINE ESPCSI(&
 ! --- INPUT -----------------------------------------------------------------
 & YDCST,YDGEOMETRY,YDLDDH,YDRIP,YDDYN,YDEDYN,KM,KMLOC,KSTA,KEND,LDONEM,&
 ! --- INOUT -----------------------------------------------------------------
 & PSPVORG,PSPDIVG,PSPTG,PSPSPG,&
 & PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG)  

!**** *ESPCSI* - SPECTRAL SPACE COMPUTATIONS FOR ALADIN
!                 SEMI-IMPLICIT SCHEME (HYDROSTATIC MODEL)

! !! The structure of ESPCSI must remain very close to SPCSI !!
!    The only differences between ESPCSI and SPCSI are the following ones:
!    * use of YDLEP%RLEPDIM instead of RLAPDI.
!    * use of YDLEP%RLEPINM instead of RLAPIN.
!    * for geometries with significant variations of the mapping factor M,
!      use the ALADIN designed option LESIDG instead of LSIDG.
!    * LIMPF currently not coded; this option has probably a sense for LAM
!      domains where the zonal coordinate exactly follows geographic latitudes,
!      and only in this case.

!     Purpose.
!     --------
!       Compute the celebrated spectral part of the equations:
!        semi-implicit correction for hydrostatic model

!**   Interface.
!     ----------
!        *CALL* *ESPCSI(...)

!        Explicit arguments :
!        --------------------
!         INPUT:
!          KM        - Zonal wavenumber
!          KMLOC   - Zonal wavenumber (DM-local numbering)
!          KSTA      - first column processed
!          KEND      - last column processed
!          LDONEM    - T if only one m if processed

!         INOUT:
!          PSPVORG   - Vorticity columns
!          PSPDIVG   - divergence columns
!          PSPTG     - temperature columns
!          PSPSPG    - surface pressure
!          PSPTNDSI_VORG -vorticity   SI tendencies
!          PSPTNDSI_DIVG -Divergence  SI tendencies
!          PSPTNDSI_TG -Temperature SI tendencies

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
!        Karim YESSAD (Feb 2005, moving the hydrostatic SI part in ESPCSI)

!     Modifications.
!     --------------
!        K. YESSAD: (Aug 2005). SI scheme in ALADIN done like in ARPEGE.
!        A. Bogatchev: (Jun 2008). Phasing cy34
!        F. Vana and NEC: (September 2008) - optimization 
!        K. Yessad (Aug 2009): remove LSITRIC option.
!        I. Santos, I. Martinez (Mar 2010): Variable map factor for RTM.
!        R. El Khatib and F. Voitus (Mar 2011): optimisations.
!        K. Yessad (Feb 2012): tests on CDCONF and LL3D in the caller, simplifications.
!        B. Bochenek 15-04-2013 - Phasing cy40, coherence with modified modules
!        B. Bochenek (Apr 2015): Phasing: move some variables.
!        O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMMP0   , ONLY : MYSETV
USE YOMDYN   , ONLY : TDYN
USE YEMDYN   , ONLY : TEDYN
USE YOMLDDH  , ONLY : TLDDH
USE YOMRIP   , ONLY : TRIP
USE YOMCST   , ONLY : TCST
USE YOMCVER      , ONLY : LVERTFE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCST),     INTENT(IN)    :: YDCST
TYPE(GEOMETRY), INTENT(IN)    :: YDGEOMETRY
TYPE(TLDDH)    ,INTENT(IN)    :: YDLDDH
TYPE(TRIP)     ,INTENT(IN)    :: YDRIP
TYPE(TDYN)     ,INTENT(IN)    :: YDDYN
TYPE(TEDYN)    ,INTENT(IN)    :: YDEDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KM
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOC
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
LOGICAL           ,INTENT(IN)    :: LDONEM 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVORG(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIVG(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTG(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPG(:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_VORG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_DIVG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_TG(:,:)

!     ------------------------------------------------------------------

REAL(KIND=JPRB), ALLOCATABLE :: ZSDIVP (:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPDIVP(:,:)

REAL(KIND=JPRB) :: ZSDIV(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZHELP(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZST(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSP(KSTA:KEND)

REAL(KIND=JPRB), ALLOCATABLE :: ZSDIVPL (:,:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPDIVPL(:,:,:)

REAL(KIND=JPRB) :: ZPD(0:YDGEOMETRY%YRDIM%NSMAX), ZPE(0:YDGEOMETRY%YRDIM%NSMAX), ZPF(0:YDGEOMETRY%YRDIM%NSMAX)

INTEGER(KIND=JPIM) :: II, IN, IM, ISP, IOFF, IS0, IS02, ISE, ISPCOL, JLEV, JN, JSP

REAL(KIND=JPRB) :: ZBDT, ZBDT2
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "mxmaop.h"

#include "mxptma.h"
#include "mxture.h"
#include "mxturs.h"

#include "abor1.intfb.h"
#include "sigam.intfb.h"
#include "sitnu.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ESPCSI',0,ZHOOK_HANDLE)
ASSOCIATE(YDLAP=>YDGEOMETRY%YRLAP, YDLEP=>YDGEOMETRY%YRELAP)
ASSOCIATE( &
 & SIMO=>YDDYN%SIMO, &
 & SIMI=>YDDYN%SIMI, &
 & SIVP=>YDDYN%SIVP, &
 & SIHEG=>YDDYN%SIHEG, &
 & SIHEG2=>YDDYN%SIHEG2, &
 & TDT=>YDRIP%TDT, &
 & RBTS2=>YDDYN%RBTS2, &
 & LIMPF=>YDDYN%LIMPF, &
 & NISNAX=>YDGEOMETRY%YREDIM%NISNAX, &
 & NSMAX=>YDGEOMETRY%YRDIM%NSMAX, &
 & NPTRSV=>YDGEOMETRY%YRMP%NPTRSV, NPTRSVF=>YDGEOMETRY%YRMP%NPTRSVF, NSPEC2V=>YDGEOMETRY%YRMP%NSPEC2V, &
 & NSPEC2VF=>YDGEOMETRY%YRMP%NSPEC2VF, &
 & RSTRET=>YDGEOMETRY%YRGEM%RSTRET, &
 & LESIDG=>YDEDYN%LESIDG, &
 & NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & LRSIDDH=>YDLDDH%LRSIDDH, &
 & ESCGMAP=>YDGEOMETRY%YSPGEOM%ESCGMAP)
!     ------------------------------------------------------------------

! LESIDG: like LSIDG but for ALADIN in the case where:
!  - large domain where the variations of the mapping factor M are large.
!  - type of projections where M is quasi constant on a latitude and
!    M can be written as a low order polynomial of mu=sin(latitude).

IF (LIMPF) CALL ABOR1(' ESPCSI: LIMPF=T not yet available in ALADIN')

!     ------------------------------------------------------------------

!*       1.    MEMORY TRANSFER.
!              ----------------

IF (LRSIDDH) THEN
  ! DDH memory transfer
  IF (LIMPF) PSPTNDSI_VORG(1:NFLEVG,KSTA:KEND)=-PSPVORG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_DIVG(1:NFLEVG,KSTA:KEND)=-PSPDIVG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_TG  (1:NFLEVG,KSTA:KEND)=-PSPTG  (1:NFLEVG,KSTA:KEND)
  !the case of surface pressure has not been treated yet
ENDIF

!     ------------------------------------------------------------------

!*       2.    SEMI-IMPLICIT SPECTRAL COMPUTATIONS.
!              ------------------------------------

ALLOCATE(ZSDIVP (NFLEVG,MAX(NSPEC2V,NSPEC2VF)))
ALLOCATE(ZSPDIVP(NFLEVG,MAX(NSPEC2V,NSPEC2VF)))
 
!*        2.1  Preliminary initialisations.
 
IF (LDONEM) THEN
  IOFF=NPTRSVF(MYSETV)-1
ELSE
  IOFF=NPTRSV(MYSETV)-1
ENDIF
ISPCOL=KEND-KSTA+1

ZBDT=RBTS2*TDT
ZBDT2=(ZBDT*RSTRET)**2

IF (LESIDG) THEN
  ALLOCATE(ZSDIVPL (NFLEVG,0:NISNAX(KM),4))
  ALLOCATE(ZSPDIVPL(NFLEVG,0:NISNAX(KM),4))
  DO JN=0,NSMAX
    ZPD(JN)=ESCGMAP(1) 
    ZPE(JN)=ESCGMAP(2)
    ZPF(JN)=ESCGMAP(3)
  ENDDO
ENDIF

!*        2.3  Computes right-hand side of Helmholtz equation.

IF (.NOT.LDONEM) CALL GSTATS(1655,0)
CALL SIGAM(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZSDIV,PSPTG(:,KSTA:KEND),PSPSPG(KSTA:KEND),ISPCOL,NFLEVG)
IF (.NOT.LDONEM) CALL GSTATS(1655,1)
IF (.NOT.LDONEM) CALL GSTATS(1656,0)

IF (LESIDG) THEN
  !  km such as "n always >0": zsdiv=rlepinm*zzspdivg-zbdt*zsdiv
  !  km such as "n can be 0": zsdiv=zzspdivg-zbdt*rlepdim*zsdiv
  IF (KM > 0) THEN
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
        IN=YDLAP%NVALUE(JSP+IOFF)
        IM=YDLEP%MVALUE(JSP+IOFF)
        ISP=YDLEP%NPME(IM)+IN
        ZSDIV(JLEV,JSP)=YDLEP%RLEPINM(ISP)*PSPDIVG(JLEV,JSP)&
         & -ZBDT*ZSDIV(JLEV,JSP)
      ENDDO
    ENDDO 
  ELSE
    DO JSP=KSTA,KEND
      DO JLEV=1,NFLEVG
        IN=YDLAP%NVALUE(JSP+IOFF)
        IM=YDLEP%MVALUE(JSP+IOFF)
        ISP=YDLEP%NPME(IM)+IN
        ZSDIV(JLEV,JSP)=PSPDIVG(JLEV,JSP)&
         & -ZBDT*YDLEP%RLEPDIM(ISP)*ZSDIV(JLEV,JSP)
      ENDDO
    ENDDO
  ENDIF

ELSE

!cdir noloopchg
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      IN=YDLAP%NVALUE(JSP+IOFF)
      IM=YDLEP%MVALUE(JSP+IOFF)
      ISP=YDLEP%NPME(IM)+IN
      ZSDIV(JLEV,JSP)=PSPDIVG(JLEV,JSP)-ZBDT*YDLEP%RLEPDIM(ISP)*ZSDIV(JLEV,JSP)
    ENDDO
  ENDDO

ENDIF

! some LIMPF=T code here: LIMPF=T not yet coded in ALADIN.

!*        2.4  Solve Helmholtz equation

!        Current space --> vertical eigenmodes space.

IF (.NOT.LDONEM) CALL GSTATS(2020,0)
CALL MXMAOP(SIMI,1,NFLEVG,ZSDIV,1,NFLEVG,ZSDIVP(:,KSTA:KEND),1,NFLEVG,&
 & NFLEVG,NFLEVG,ISPCOL)
IF (.NOT.LDONEM) CALL GSTATS(2020,1)

IF (LESIDG) THEN
  ! - memory transfer zsdivp -> zsdivpl
  ! - km such as "n always >0": call to MXTURS, result in zspdivpl
  ! - km such as "n can be 0": two calls to MXTURE, result in zspdivpl
  ! - memory transfer zspdivpl -> zspdivp
  IS0=YDLAP%NSE0L(KMLOC)
  IS02=0
  ZSDIVPL(:,:,:)=0.0_JPRB
  ZSPDIVPL(:,:,:)=0.0_JPRB

  DO JN=0,NISNAX(KM)
    ISE=KSTA+4*JN
    ZSDIVPL(:,JN,1:4)=ZSDIVP(:,ISE:ISE+3)
  ENDDO

  IF (KM > 0) THEN

    !               Inversion of a symmetric matrix.
    II = 4
    CALL MXTURS(NISNAX(KM)+1,NFLEVG,NFLEVG,II,&
     & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
     & ZSDIVPL,ZSPDIVPL)
  ELSE

    !               Inversion of a non-symmetric matrix.
    II = 2
    CALL MXTURE(NISNAX(KM)+1,NFLEVG,NFLEVG,II,-2,.TRUE.,&
     & SIHEG(1,IS0+1,1),SIHEG(1,IS0+1,2),SIHEG(1,IS0+1,3),&
     & ZSDIVPL,ZSPDIVPL)
    CALL MXTURE(NISNAX(KM)+1,NFLEVG,NFLEVG,II,3,.FALSE.,&
     & SIHEG(1,IS0+1,1),SIHEG2(1,IS02+1,2),&
     & SIHEG2(1,IS02+1,3),ZSDIVPL,ZSPDIVPL)
  ENDIF

  DO JN=0,NISNAX(KM)
    ISE=KSTA+4*JN
    ZSPDIVP(:,ISE:ISE+3)=ZSPDIVPL(:,JN,1:4)
  ENDDO

ELSE

  !             Case with NO Stretching :

  IF (LIMPF) THEN
    ! LIMPF not yet coded in ALADIN.
  ELSE

    !                 Inversion of a diagonal system (Helmholtz equation)
    !                 --> (SIMI*DIVprim(t+dt)).

!cdir noloopchg
    DO JSP=KSTA,KEND
!cdir nodep
      DO JLEV=1,NFLEVG
        IN=YDLAP%NVALUE(JSP+IOFF)
        IM=YDLEP%MVALUE(JSP+IOFF)
        ISP=YDLEP%NPME(IM)+IN
        ZSPDIVP(JLEV,JSP)=ZSDIVP(JLEV,JSP)&
         & /(1.0_JPRB-ZBDT2*SIVP(JLEV)*YDLEP%RLEPDIM(ISP))
      ENDDO
    ENDDO
  ENDIF
ENDIF

!           Vertical eigenmodes space --> current space.

IF (.NOT.LDONEM) CALL GSTATS(2020,0)
CALL MXMAOP(SIMO,1,NFLEVG,ZSPDIVP(:,KSTA:KEND),1,NFLEVG,PSPDIVG(:,KSTA:KEND),1,&
 & NFLEVG,NFLEVG,NFLEVG,ISPCOL)
IF (.NOT.LDONEM) CALL GSTATS(2020,1)

IF (LESIDG) THEN
  ! - set zsdivpl and zspdivpl to zero.
  ! - memory transfer zzspdivg -> zsdivpl
  ! - call to MXPTMA: zspdivpl=(M**2)*zsdivpl
  ! - zhelp=zsdivpl
  ZSDIVPL(:,:,:)=0.0_JPRB
  ZSPDIVPL(:,:,:)=0.0_JPRB

  !           Reorganisation of ZSDIVP (Back to the USSR)

  DO JN=0,NISNAX(KM)
    ISE=KSTA+4*JN
    ZSDIVPL(:,JN,1:4)=PSPDIVG(:,ISE:ISE+3)
  ENDDO

  !        ZSPDIV=(DIVprim(t+dt)) --> ZPSPDIVG=(GMBAR**2 * DIVprim(t+dt)).
  CALL MXPTMA(NISNAX(KM)+1,NFLEVG,NFLEVG,II,&
   & ZPD,ZPE,ZPF,ZPE,ZPF,&
   & ZSDIVPL,ZSPDIVPL)

  !           Reorganisation of ZSPDIVPL

  DO JN=0,NISNAX(KM)
    ISE=KSTA+4*JN
    ZHELP(:,ISE:ISE+3)=ZSPDIVPL(:,JN,1:4)
  ENDDO

ELSE

  !       ZZSPDIVG=(DIVprim(t+dt)) --> ZHELP=(GMBAR**2 * DIVprim(t+dt)) .
  !       where GMBAR=max(M)

  IF (.NOT.LDONEM) CALL GSTATS(1656,0)
!cdir noloopchg
  DO JSP=KSTA,KEND
!cdir nodep
    DO JLEV=1,NFLEVG
      ZHELP(JLEV,JSP)=PSPDIVG(JLEV,JSP)*RSTRET*RSTRET
    ENDDO
  ENDDO
  IF (.NOT.LDONEM) CALL GSTATS(1656,1)

ENDIF

!       If LESIDG:
!         (GM**2 * DIVprim(t+dt)) --> [ tau * (GM**2 * DIVprim(t+dt)) ]
!                                 and [  nu * (GM**2 * DIVprim(t+dt)) ]
!       or if not LESIDG:
!         (GMBAR**2 * DIVprim(t+dt)) --> [ tau * (GMBAR**2 * DIVprim(t+dt)) ]
!                                    and [  nu * (GMBAR**2 * DIVprim(t+dt)) ]
!         where GMBAR=max(M)
 
IF (.NOT.LDONEM) CALL GSTATS(1657,0)
CALL SITNU(LVERTFE, YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ZHELP,ZST,ZSP,ISPCOL)
IF (.NOT.LDONEM) CALL GSTATS(1657,1)

!*       2.5  Increment Temperature and surface pressure
    
!cdir noloopchg
IF (.NOT.LDONEM) CALL GSTATS(1656,0)
DO JSP=KSTA,KEND
!cdir nodep
  DO JLEV=1,NFLEVG
    PSPTG(JLEV,JSP)=PSPTG(JLEV,JSP)-ZBDT*ZST(JLEV,JSP)
  ENDDO
  PSPSPG(JSP)=PSPSPG(JSP)-ZBDT*ZSP(JSP)
ENDDO
IF (.NOT.LDONEM) CALL GSTATS(1656,1)
!cdir nodep

DEALLOCATE(ZSDIVP)
DEALLOCATE(ZSPDIVP)

!*       2.6  Increment vorticity

! LIMPF=T not yet coded in ALADIN.

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF SI TERM AT t+dt FOR DDH.
!              ---------------------------------------

IF (LRSIDDH) THEN
  IF (LIMPF) PSPTNDSI_VORG(1:NFLEVG,KSTA:KEND)=&
   & PSPTNDSI_VORG(1:NFLEVG,KSTA:KEND) + PSPVORG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_DIVG(1:NFLEVG,KSTA:KEND)=PSPTNDSI_DIVG(1:NFLEVG,KSTA:KEND)&
   & + PSPDIVG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_TG(1:NFLEVG,KSTA:KEND)=PSPTNDSI_TG(1:NFLEVG,KSTA:KEND)&
   & + PSPTG(1:NFLEVG,KSTA:KEND)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ESPCSI',1,ZHOOK_HANDLE)
END SUBROUTINE ESPCSI
