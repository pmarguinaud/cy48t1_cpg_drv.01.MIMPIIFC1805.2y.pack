SUBROUTINE ESPCHORAD(YDGEOMETRY,YDMODEL,KSTA,KEND,LDONEM,&
 & PSPVORG,PSPDIVG,PSPTG,PSPGFL,PSPSPG,&
 & PSPSPDG,PSPSVDG)  

!**** *ESPCHORAD* - SPECTRAL SPACE COMPUTATIONS FOR ALADIN - ADJOINT
!                   HORIZONTAL DIFFUSION 

!     Purpose.
!     --------
!       Compute the adjoint of the celebrated spectral part of the equations:
!        horizontal diffusion

!**   Interface.
!     ----------
!        *CALL* *ESPCHORAD(...)

!        Explicit arguments :
!        --------------------  KSTA   - first column processed
!                              KEND   - last column processed
!                              LDONEM - T if only one m if processed
!                              PSPVORG- sensitivity to vorticity - columns
!                              PSPDIVG- sens. to divergence - columns
!                              PSPTG  - sens. to temperature - columns
!                              PSPGFL - GFL
!                              PSPSPG - sens. to surface pressure
!                              PSPSPDG- sens. to nhs pressure departure - col.
!                              PSPSVDG- sens. to nhs vertical divergence - col.

!     Method.
!     -------
!       There is a great description in ESPCM

!     Externals.
!     ----------
!     Reference.
!     ----------
!        ARPEGE/ALADIN documentation

!     Author.
!     -------
!        Andras Horanyi (original, July 1995)
!        Karim YESSAD (Feb 2005, moving the horiz. diffusion part in ESPCHORAD)

!     Modifications.
!     --------------
!   K. Yessad (Sep 2008): prune enhanced diffusion (lfrein).
!   F. Vana 13-Jan-2009: supported diffusion if LSLHD
!   F. Vana  13-Oct-2009: Optimization (improved vect., OpenMP)
!   K. Yessad (Feb 2012): tests on CDCONF in the caller.
!   A.Bogatchev 11-04-2013 phasinc cy40, coherence with modified modules
!                          and renamed nsmelists and functions
!   B. Bochenek (Apr 2015): Phasing: move some variables.
!   O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     ------------------------------------------------------------------

USE TYPE_MODEL , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMMP0   , ONLY : MYSETV
USE YOMCT0   , ONLY : LNHDYN
USE YOMDYNA  , ONLY : LSLHD_STATIC ,LSLHD_W
USE OML_MOD  , ONLY : OML_GET_NUM_THREADS

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT)   :: YDGEOMETRY
TYPE(MODEL),    INTENT(INOUT)   :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
LOGICAL           ,INTENT(IN)    :: LDONEM 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVORG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIVG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPGFL(:,:,:) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPG(KSTA:KEND) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPDG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSVDG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IM, IN, IOFF, JLEV, JSP, JGFL
INTEGER(KIND=JPIM) :: ITHRDS,IBLKS,JI_KSTA,I_KEND
INTEGER(KIND=JPIM) :: ISP_VEC(KSTA:KEND)

REAL(KIND=JPRB) :: ZDT
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

! The use of nvalue must match the definition in suelap.F

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ESPCHORAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDYN=>YDMODEL%YRML_DYN%YRDYN, YDEDYN=>YDMODEL%YRML_DYN%YREDYN, YDML_GCONF=>YDMODEL%YRML_GCONF, &
 & YDLAP=>YDGEOMETRY%YRLAP, YDLEP=>YDGEOMETRY%YRELAP)
ASSOCIATE( &
 & LSPT=>YDML_GCONF%YRDIMF%LSPT, RCORDIT=>YDDYN%RCORDIT, TDT=>YDML_GCONF%YRRIP%TDT, &
 & YCOMP=>YDML_GCONF%YGFL%YCOMP, &
 & RDIVORE=>YDEDYN%RDIVORE, RDIDIVE=>YDEDYN%RDIDIVE, RDSVORE=>YDEDYN%RDSVORE, &
 & RDSDIVE=>YDEDYN%RDSDIVE, RDITE=>YDEDYN%RDITE, RDIGFLE=>YDEDYN%RDIGFLE, &
 & RDIPDE=>YDEDYN%RDIPDE, RDIVDE=>YDEDYN%RDIVDE, RDISPE=>YDEDYN%RDISPE, &
 & NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, NPTRSV=>YDGEOMETRY%YRMP%NPTRSV, NPTRSVF=>YDGEOMETRY%YRMP%NPTRSVF )
!     ------------------------------------------------------------------

IF (LDONEM) THEN
  IOFF=NPTRSVF(MYSETV)-1
ELSE
  IOFF=NPTRSV(MYSETV)-1
ENDIF

!     ------------------------------------------------------------------

!*       2.    Horizontal diffusion
!              --------------------

ZDT=ABS(TDT)
CALL GSTATS(1035,0)

DO JSP = KSTA, KEND
  IN=YDLAP%NVALUE(JSP+IOFF)
  IM=YDLEP%MVALUE(JSP+IOFF)
  ISP_VEC(JSP)=YDLEP%NPME(IM)+IN
ENDDO

!$OMP PARALLEL&
!$OMP& PRIVATE(JSP,JLEV,IN,IM,JGFL,ITHRDS,IBLKS,JI_KSTA,I_KEND)
ITHRDS=OML_GET_NUM_THREADS()
IBLKS=(KEND-KSTA+ITHRDS)/ITHRDS

! * Log(prehyds):
!cdir nodep
!$OMP DO SCHEDULE(STATIC)
DO JSP=KSTA,KEND
  PSPSPG (JSP)=PSPSPG (JSP)/(1.0_JPRB+ZDT*RDISPE(ISP_VEC(JSP)))
ENDDO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO JI_KSTA =KSTA,KEND,IBLKS
  I_KEND=MIN(KEND,JI_KSTA+IBLKS-1)

  ! * VOR,DIV:
  IF (LSLHD_W.AND.(.NOT.LSLHD_STATIC)) THEN
    ! * Additional diffusion
    DO JSP=JI_KSTA,I_KEND
      DO JLEV=1,NFLEVG
!cdir nodep
        PSPVORG(JLEV,JSP)=PSPVORG(JLEV,JSP)&
         & /(1.0_JPRB+ZDT*RDSVORE(JLEV,ISP_VEC(JSP)))
        PSPDIVG(JLEV,JSP)=PSPDIVG(JLEV,JSP)&
         & /(1.0_JPRB+ZDT*RDSDIVE(JLEV,ISP_VEC(JSP)))
      ENDDO
    ENDDO
  ENDIF
  DO JSP=JI_KSTA,I_KEND
    DO JLEV=1,NFLEVG
!cdir nodep
      PSPVORG(JLEV,JSP)=PSPVORG(JLEV,JSP)&
       & /(1.0_JPRB+ZDT*RDIVORE(JLEV,ISP_VEC(JSP)))
      PSPDIVG(JLEV,JSP)=PSPDIVG(JLEV,JSP)&
       & /(1.0_JPRB+ZDT*RDIDIVE(JLEV,ISP_VEC(JSP)))
    ENDDO
  ENDDO

  ! * Temperature:
  IF (LSPT) THEN
    DO JSP=JI_KSTA,I_KEND
      DO JLEV=1,NFLEVG
!cdir nodep
        PSPSPG (JSP)=PSPSPG (JSP) +&
         & RCORDIT(JLEV)*PSPTG(JLEV,JSP)*ZDT*RDITE(JLEV,ISP_VEC(JSP))&
         & /(1.0_JPRB+ZDT*RDITE(JLEV,ISP_VEC(JSP)))  
        PSPTG  (JLEV,JSP)=PSPTG(JLEV,JSP)&
         & /(1.0_JPRB+ZDT*RDITE(JLEV,ISP_VEC(JSP)))
      ENDDO
    ENDDO
  ENDIF

  ! * NH variables:
  IF (LNHDYN) THEN
    DO JSP=JI_KSTA,I_KEND
      DO JLEV=1,NFLEVG
!cdir nodep
        PSPSPDG(JLEV,JSP)=PSPSPDG(JLEV,JSP)&
         & /(1.0_JPRB+ZDT*RDIPDE(JLEV,ISP_VEC(JSP)))
        PSPSVDG(JLEV,JSP)=PSPSVDG(JLEV,JSP)&
         & /(1.0_JPRB+ZDT*RDIVDE(JLEV,ISP_VEC(JSP)))
      ENDDO
    ENDDO
  ENDIF

  ! * GFL variables:
  DO JGFL=1,YDML_GCONF%YGFL%NUMFLDS
    IF(YCOMP(JGFL)%LSP) THEN
      DO JSP=JI_KSTA,I_KEND
        DO JLEV=1,NFLEVG
          PSPGFL(JLEV,JSP,YCOMP(JGFL)%MPSP) =&
           & PSPGFL(JLEV,JSP,YCOMP(JGFL)%MPSP)/&
           & (1.0_JPRB+ZDT*RDIGFLE(JLEV,ISP_VEC(JSP),YCOMP(JGFL)%MPSP))  
        ENDDO
      ENDDO
    ENDIF
  ENDDO

ENDDO
!$OMP END DO

!$OMP END PARALLEL

CALL GSTATS(1035,1)

!     ------------------------------------------------------------------


END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ESPCHORAD',1,ZHOOK_HANDLE)
END SUBROUTINE ESPCHORAD
