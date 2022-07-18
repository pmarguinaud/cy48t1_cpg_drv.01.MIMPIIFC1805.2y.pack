SUBROUTINE SPCHORAD(YDGEOMETRY,YDML_GCONF,YDDYN,KSTA,KEND,PSPVOR,PSPDIV,PSPT,PSPGFL,PSPSP,PSPTALLG)

!     ------------------------------------------------------------------

!**** *SPCHORAD* - HORIZONTAL SPECTRAL SPACE COMPUTATIONS - ADJOINT

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SPCHORAD(..)

!        Explicit arguments :
!        --------------------  
!                              KSTA    - First column processed
!                              KEND    - Last columns processed
!                              PSPVOR - Vorticity columns
!                              PSPDIV - Divergence columns
!                              PSPT   - Temperature columns
!                              PSPGFL - GFL columns
!                              PSPSP  - Surface Pressure
!                              PSPTALLG - Temperature for all levels

!     Method.
!     -------
!      SPCHORAD  works in the "normal" spectral domain where a
!                given processor owns NFLEVL levels and NSPEC2
!                spectral coefficients. NSPEC2 depends on MYSETW only,
!                and NFLEVL depends on MYSETV only

!     Externals.
!     ----------
!        MXTURHD

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Lars Isaksen *ECMWF*
!      Original   : 96-12-12  Based on spcad.F and spcdm.F plus 
!                             modifications by J.Latour

!     Modifications.
!     --------------
!      K. YESSAD 03-04-16: bug correction in call MXTURHD.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K.Yessad (Dec 2003): cleaning in horizontal diffusion
!      K. Yessad 15-May-2006: memory optimisations for stretched geometry
!      K. Yessad (Sep 2008): prune enhanced diffusion (lfrein).
!      F. Vana 13-Jan_2009 : supported diffusion when LSLHD (not yet for NH)
!      K. Yessad (Feb 2012): tests on CDCONF and NCONF in the caller.
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK
USE YOMMP0                 , ONLY : MYSETV, NPRTRV
USE YOMDYNA                , ONLY : YRDYNA
USE YOMDYN                 , ONLY : TDYN
USE YOMNUD                 , ONLY : LNUDG

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVOR(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIV(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPT(:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPGFL(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSP(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPTALLG(:,:) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZSPX(YDGEOMETRY%YRDIMV%NFLEVL*YDML_GCONF%YRDIMF%NS3D+1,0:YDGEOMETRY%YRDIM%NSMAX,2)
REAL(KIND=JPRB) :: ZSPY(YDGEOMETRY%YRDIMV%NFLEVL*YDML_GCONF%YRDIMF%NS3D+1,0:YDGEOMETRY%YRDIM%NSMAX,2)
REAL(KIND=JPRB),ALLOCATABLE :: ZSPXTG(:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZSPX_HDS(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZSPY_HDS(:,:,:)

INTEGER(KIND=JPIM) :: IJSE, ILEV0, ILL, IM, IN, IS0, INUM_GFL0, INUM_THE0,&
 & JLEV, JLEVG, JMLOC, JN, JSP, JGFL, IVTHS

REAL(KIND=JPRB) :: ZLAPDI, ZLAPIN
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "gathert.intfb.h"
#include "mxturhd.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPCHORAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
  & YDMP=>YDGEOMETRY%YRMP, YDLAP=>YDGEOMETRY%YRLAP, YDRIP=>YDML_GCONF%YRRIP,YGFL=>YDML_GCONF%YGFL, &
  & YDDIMF=>YDML_GCONF%YRDIMF)
ASSOCIATE(NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & NSMAX=>YDDIM%NSMAX, NUMP=>YDDIM%NUMP, &
 & LSPT=>YDDIMF%LSPT, NFTHER=>YDDIMF%NFTHER, NS3D=>YDDIMF%NS3D, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLEVL=>YDDIMV%NFLEVL, &
 & LSTRHD=>YDDYN%LSTRHD, RCORDIT=>YDDYN%RCORDIT, RDHI=>YDDYN%RDHI, &
 & RDHS=>YDDYN%RDHS, RDIDIV=>YDDYN%RDIDIV, RDIGFL=>YDDYN%RDIGFL, &
 & RDISP=>YDDYN%RDISP, RDITG=>YDDYN%RDITG, RDIVOR=>YDDYN%RDIVOR, &
 & RDSDIV=>YDDYN%RDSDIV, RDSVOR=>YDDYN%RDSVOR, &
 & NPTRLL=>YDMP%NPTRLL, &
 & TDT=>YDRIP%TDT)
!     ------------------------------------------------------------------

!  ILEV0 is the ordinal of the first level for this processor in the interval 1:NFLEVG
ILEV0=NPTRLL(MYSETV) - 1

!     ------------------------------------------------------------------

!*       3     NUDGING.
!              --------

IF (LNUDG) THEN
  CALL ABOR1('SPCHORAD: NO ADJOINT OF NUDGING' )
ENDIF

!     ------------------------------------------------------------------

!*       2.    MAIN HORIZONTAL DIFFUSION.
!              --------------------------

! no sponge in the adjoint model

IF (LSTRHD) THEN

  !           2.1  Diffusion when stretching.

  IF (YRDYNA%LSLHD_W.AND.(.NOT.YRDYNA%LSLHD_STATIC)) THEN

    ! * Additional diffusion
    !   The order in ZSPX_HDS, ZSPY_HDS and also in RDHS is assumed to be
    !   the following one: VOR,DIV

    ! * First compute the number of prognostic variables on which the HDS
    !   diffusion is applied (currently IVTHS=(ILL/NFLEVL) is between 0 and 2).
    !   This piece of code must remain consistent with the equivalent one done
    !   in SUHDU and SUALDYNB.
    IVTHS=0
    IF (YRDYNA%LSLHD_W) IVTHS=IVTHS+2
    ILL=IVTHS*NFLEVL

    ! * Allocate ZSPX_HDS,ZSPY_HDS, then fill them with zeros.
    ALLOCATE(ZSPX_HDS(ILL,0:NSMAX,2))
    ALLOCATE(ZSPY_HDS(ILL,0:NSMAX,2))
    ZSPX_HDS(:,:,:)=0.0_JPRB
    ZSPY_HDS(:,:,:)=0.0_JPRB

    DO JMLOC=1,NUMP

      IM=YDLAP%MYMS(JMLOC)

      DO JN=IM,NSMAX
        IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
        ! * Vor,Div.
        ZLAPDI=YDLAP%RLAPDI(JN)
        DO JLEV=1,NFLEVL
          ZSPY_HDS(JLEV       ,JN,1)=PSPVOR(JLEV,IJSE  )*ZLAPDI
          ZSPY_HDS(JLEV       ,JN,2)=PSPVOR(JLEV,IJSE+1)*ZLAPDI
          ZSPY_HDS(JLEV+NFLEVL,JN,1)=PSPDIV(JLEV,IJSE  )*ZLAPDI
          ZSPY_HDS(JLEV+NFLEVL,JN,2)=PSPDIV(JLEV,IJSE+1)*ZLAPDI
        ENDDO
      ENDDO

      IS0=YDLAP%NSE0L(JMLOC)
      CALL MXTURHD(NSMAX+1-IM,ILL,ILL,-3,.TRUE.,&
       & RDHS(1,IS0+1,1),RDHS(1,IS0+1,3),&
       & ZSPY_HDS(1,IM,1),ZSPX_HDS(1,IM,1))  
      CALL MXTURHD(NSMAX+1-IM,ILL,ILL,2,.FALSE.,&
       & RDHS(1,IS0+1,1),RDHS(1,IS0+1,2),&
       & ZSPY_HDS(1,IM,1),ZSPX_HDS(1,IM,1))  

      IF (IM > 0) THEN
        CALL MXTURHD(NSMAX+1-IM,ILL,ILL,-3,.TRUE.,&
         & RDHS(1,IS0+1,1),RDHS(1,IS0+1,3),&
         & ZSPY_HDS(1,IM,2),ZSPX_HDS(1,IM,2))  
        CALL MXTURHD(NSMAX+1-IM,ILL,ILL,2,.FALSE.,&
         & RDHS(1,IS0+1,1),RDHS(1,IS0+1,2),&
         & ZSPY_HDS(1,IM,2),ZSPX_HDS(1,IM,2))  
      ENDIF

      DO JN=IM,NSMAX
        IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
        ! * Vor,Div.
        ZLAPIN=YDLAP%RLAPIN(JN)
        DO JLEV=1,NFLEVL
          PSPVOR(JLEV,IJSE  )=ZSPX_HDS(JLEV         ,JN,1)*ZLAPIN
          PSPVOR(JLEV,IJSE+1)=ZSPX_HDS(JLEV         ,JN,2)*ZLAPIN
          PSPDIV(JLEV,IJSE  )=ZSPX_HDS(JLEV+  NFLEVL,JN,1)*ZLAPIN
          PSPDIV(JLEV,IJSE+1)=ZSPX_HDS(JLEV+  NFLEVL,JN,2)*ZLAPIN
        ENDDO
      ENDDO
            
    ENDDO

    DEALLOCATE(ZSPX_HDS)
    DEALLOCATE(ZSPY_HDS)

  ENDIF

  ! * The order in ZSPX, ZSPY and also in RHDI is assumed to be the
  !   following one for hydrostatic model:
  !   - VOR,DIV.
  !   - Thermodynamic variables.
  !   - GFL variables.
  !   - GMVS variables.
  INUM_THE0=2
  INUM_GFL0=2+NFTHER

  IF (LSPT.AND.NPRTRV > 1) ALLOCATE (ZSPXTG(NFLEVG,0:NSMAX,2,NUMP))

  ILL=NFLEVL*NS3D+1

  ZSPX(:,:,:)=0.0_JPRB
  ZSPY(:,:,:)=0.0_JPRB

  DO JMLOC=1,NUMP

    IM=YDLAP%MYMS(JMLOC)

    DO JN=IM,NSMAX
      IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
      ZSPY(ILL,JN,1)=PSPSP(IJSE  )
      ZSPY(ILL,JN,2)=PSPSP(IJSE+1)
    ENDDO

    DO JN=IM,NSMAX
      IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
      ZLAPDI=YDLAP%RLAPDI(YDLAP%NVALUE(IJSE))
      DO JLEV=1,NFLEVL
        ZSPY(JLEV          ,JN,1)=PSPVOR(JLEV,IJSE  )*ZLAPDI
        ZSPY(JLEV          ,JN,2)=PSPVOR(JLEV,IJSE+1)*ZLAPDI
        ZSPY(JLEV +  NFLEVL,JN,1)=PSPDIV(JLEV,IJSE  )*ZLAPDI
        ZSPY(JLEV +  NFLEVL,JN,2)=PSPDIV(JLEV,IJSE+1)*ZLAPDI
      ENDDO
      IF (LSPT) THEN
        DO JLEV=1,NFLEVL
          ZSPY(JLEV +INUM_THE0*NFLEVL,JN,1)=PSPT(JLEV,IJSE  )
          ZSPY(JLEV +INUM_THE0*NFLEVL,JN,2)=PSPT(JLEV,IJSE+1)
        ENDDO
      ENDIF
      DO JGFL=1,NUMFLDS
        IF (YCOMP(JGFL)%LSP) THEN
          DO JLEV=1,NFLEVL
            ZSPY(JLEV+(INUM_GFL0+YCOMP(JGFL)%MPSP-1)*NFLEVL,JN,1)=&
             & PSPGFL(JLEV,IJSE  ,YCOMP(JGFL)%MPSP)  
            ZSPY(JLEV+(INUM_GFL0+YCOMP(JGFL)%MPSP-1)*NFLEVL,JN,2)=&
             & PSPGFL(JLEV,IJSE+1,YCOMP(JGFL)%MPSP)  
          ENDDO
        ENDIF
      ENDDO

      ! Also in the case NPRTRV.GT.1 we have to use the temperatures
      ! from all levels for all coefficients (PSPTALLG)

      IF (LSPT) THEN
        DO JLEVG=1,NFLEVG
          IF (NPRTRV > 1) THEN
            ZSPY(ILL,JN,1)=ZSPY(ILL,JN,1)+RCORDIT(JLEVG)*PSPTALLG(JLEVG,IJSE  )  
            ZSPY(ILL,JN,2)=ZSPY(ILL,JN,2)+RCORDIT(JLEVG)*PSPTALLG(JLEVG,IJSE+1)  
          ELSE
            ZSPY(ILL,JN,1)=ZSPY(ILL,JN,1)+RCORDIT(JLEVG)*PSPT(JLEVG,IJSE  )  
            ZSPY(ILL,JN,2)=ZSPY(ILL,JN,2)+RCORDIT(JLEVG)*PSPT(JLEVG,IJSE+1)  
          ENDIF
        ENDDO
      ENDIF

    ENDDO

    IS0=YDLAP%NSE0L(JMLOC)
    CALL MXTURHD(NSMAX+1-IM,ILL,ILL,-3,.TRUE.,&
     & RDHI(1,IS0+1,1),RDHI(1,IS0+1,3),ZSPY(1,IM,1),ZSPX(1,IM,1))  
    CALL MXTURHD(NSMAX+1-IM,ILL,ILL,2,.FALSE.,&
     & RDHI(1,IS0+1,1),RDHI(1,IS0+1,2),ZSPY(1,IM,1),ZSPX(1,IM,1))  

    IF (IM > 0) THEN
      CALL MXTURHD(NSMAX+1-IM,ILL,ILL,-3,.TRUE.,&
       & RDHI(1,IS0+1,1),RDHI(1,IS0+1,3),ZSPY(1,IM,2),ZSPX(1,IM,2))  
      CALL MXTURHD(NSMAX+1-IM,ILL,ILL,2,.FALSE.,&
       & RDHI(1,IS0+1,1),RDHI(1,IS0+1,2),ZSPY(1,IM,2),ZSPX(1,IM,2))  
    ENDIF

    IF (LSPT.AND.NPRTRV > 1) THEN
      ! Store local diffused temperature values. After the 
      ! loop they are communicated to build up global version
      ZSPXTG(ILEV0+1:ILEV0+NFLEVL,:,:,JMLOC)=&
       & ZSPX(INUM_THE0*NFLEVL+1:(INUM_THE0+1)*NFLEVL,:,:)  
    ENDIF

    DO JN=IM,NSMAX
      IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
      PSPSP(IJSE  )=ZSPX(ILL,JN,1)
      PSPSP(IJSE+1)=ZSPX(ILL,JN,2)
    ENDDO

    DO JN=IM,NSMAX
      IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
      ZLAPIN=YDLAP%RLAPIN(JN)
      IF (LSPT.AND.NPRTRV == 1) THEN
        DO JLEV=1,NFLEVL
          PSPSP(IJSE  )=PSPSP(IJSE  )-RCORDIT(JLEV)*ZSPX(JLEV+INUM_THE0*NFLEVL,JN,1)  
          PSPSP(IJSE+1)=PSPSP(IJSE+1)-RCORDIT(JLEV)*ZSPX(JLEV+INUM_THE0*NFLEVL,JN,2)  
        ENDDO
      ENDIF
      DO JLEV=1,NFLEVL
        PSPVOR(JLEV,IJSE  )=ZSPX(JLEV         ,JN,1)*ZLAPIN
        PSPVOR(JLEV,IJSE+1)=ZSPX(JLEV         ,JN,2)*ZLAPIN
        PSPDIV(JLEV,IJSE  )=ZSPX(JLEV+  NFLEVL,JN,1)*ZLAPIN
        PSPDIV(JLEV,IJSE+1)=ZSPX(JLEV+  NFLEVL,JN,2)*ZLAPIN
      ENDDO
      IF (LSPT) THEN
        DO JLEV=1,NFLEVL
          PSPT  (JLEV,IJSE  )=ZSPX(JLEV+INUM_THE0*NFLEVL,JN,1)
          PSPT  (JLEV,IJSE+1)=ZSPX(JLEV+INUM_THE0*NFLEVL,JN,2)
        ENDDO
      ENDIF
      DO JGFL=1,NUMFLDS
        IF (YCOMP(JGFL)%LSP) THEN
          DO JLEV=1,NFLEVL
            PSPGFL(JLEV,IJSE  ,YCOMP(JGFL)%MPSP)=&
             & ZSPX(JLEV+(INUM_GFL0+YCOMP(JGFL)%MPSP-1)*NFLEVL,JN,1)  
            PSPGFL(JLEV,IJSE+1,YCOMP(JGFL)%MPSP)=&
             & ZSPX(JLEV+(INUM_GFL0+YCOMP(JGFL)%MPSP-1)*NFLEVL,JN,2)  
          ENDDO
        ENDIF
      ENDDO
    ENDDO

  ENDDO

  IF (LSPT.AND.NPRTRV > 1) THEN

    ! Surface Pressure temperature update has been delayed
    ! in order to reduce communication. First we gather all 
    ! diffused temperature levels on all PE's within an A-set

    CALL GATHERT(YDGEOMETRY,ZSPXTG)
    DO JMLOC=1,NUMP
      IM=YDLAP%MYMS(JMLOC)
      DO JN=IM,NSMAX
        IJSE=YDLAP%NASM0(IM)+2*(JN-IM)
        ZLAPIN=YDLAP%RLAPIN(JN)
        DO JLEVG=1,NFLEVG
          PSPSP(IJSE  )=PSPSP(IJSE  )-RCORDIT(JLEVG)*ZSPXTG(JLEVG,JN,1,JMLOC)  
          PSPSP(IJSE+1)=PSPSP(IJSE+1)-RCORDIT(JLEVG)*ZSPXTG(JLEVG,JN,2,JMLOC)  
        ENDDO
      ENDDO
    ENDDO
    DEALLOCATE (ZSPXTG)
  ENDIF

ELSE

!           2.2  Diffusion when no stretching.

  CALL GSTATS(1035,0)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSP,JLEVG,IN,JGFL)
  DO JSP=KSTA,KEND
    IN=YDLAP%NVALUE(JSP)
    PSPSP(JSP)=PSPSP(JSP)/(1.0_JPRB+TDT*RDISP(IN))
    IF (LSPT) THEN
      IF (NPRTRV > 1) THEN
        DO JLEVG=1,NFLEVG
          PSPSP(JSP)=PSPSP(JSP)&
           & +RCORDIT(JLEVG)*PSPTALLG(JLEVG,JSP)*TDT*&
           & RDITG(JLEVG,IN)/(1.0_JPRB+TDT*RDITG(JLEVG,IN))  
        ENDDO
      ELSE
        DO JLEVG=1,NFLEVG
          PSPSP(JSP)=PSPSP(JSP)&
           & +RCORDIT(JLEVG)*PSPT(JLEVG,JSP)*TDT*&
           & RDITG(JLEVG,IN)/(1.0_JPRB+TDT*RDITG(JLEVG,IN))  
        ENDDO
      ENDIF
    ENDIF
    IF (YRDYNA%LSLHD_W.AND.(.NOT.YRDYNA%LSLHD_STATIC)) THEN
      ! * Additional diffusion 
!CDIR NODEP
      DO JLEV=1,NFLEVL
        PSPVOR(JLEV,JSP)=PSPVOR(JLEV,JSP)/(1.0_JPRB+TDT*RDSVOR(JLEV,IN))
        PSPDIV(JLEV,JSP)=PSPDIV(JLEV,JSP)/(1.0_JPRB+TDT*RDSDIV(JLEV,IN))
      ENDDO
    ENDIF
    DO JLEV=1,NFLEVL
      PSPVOR(JLEV,JSP)=PSPVOR(JLEV,JSP)/(1.0_JPRB+TDT*RDIVOR(JLEV,IN))
      PSPDIV(JLEV,JSP)=PSPDIV(JLEV,JSP)/(1.0_JPRB+TDT*RDIDIV(JLEV,IN))
    ENDDO
    IF (LSPT) THEN
      DO JLEV=1,NFLEVL
        PSPT(JLEV,JSP)=PSPT(JLEV,JSP)/(1.0_JPRB+TDT*RDITG(ILEV0+JLEV,IN))
      ENDDO
    ENDIF
    DO JGFL=1,NUMFLDS
      IF(YCOMP(JGFL)%LSP) THEN
        DO JLEV=1,NFLEVL
          PSPGFL(JLEV,JSP,YCOMP(JGFL)%MPSP) =&
           & PSPGFL(JLEV,JSP,YCOMP(JGFL)%MPSP)/&
           & (1.0_JPRB+TDT*RDIGFL(JLEV,IN,YCOMP(JGFL)%MPSP))  
        ENDDO
      ENDIF
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

  CALL GSTATS(1035,1)

ENDIF ! LSTRHD

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPCHORAD',1,ZHOOK_HANDLE)
END SUBROUTINE SPCHORAD
