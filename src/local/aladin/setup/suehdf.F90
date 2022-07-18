!option! -O nomove
SUBROUTINE SUEHDF(YDGEOMETRY,YDDYN,YDEDYN,YDML_GCONF)

!**** *SUEHDF*   - Initialize horizontal diffusion for the dynamics
!                  in the Aladin geometry.
!     Purpose.
!     --------


!**   Interface.
!     ----------
!        *CALL* *SUEHDF

!         Some remarks:
!          - The code of suehdf.F90 is valid for LAM models.
!          - Several options are available for the vertical profile (NPROFILEHD=1, 2 or 3).
!          - One option is currently available for the horizontal function:
!            * LRDISPE_EC=F: uses BDISPE and is designed for MF and all partners using LAM model.
!            * LRDISPE_EC=T: is currently not implemented in LAM models.
!          - For GFL variables, the old passive scalars are now considered
!            as obsolete features, and there are not any longer diffused;
!            GFL variables other than "q" and "O3" are not diffused.

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE/ALADIN DOCUMENTATION

!     Author.
!     -------
!      I. Ihasz CNRM/GMAP MICECO
!      Original     : 92-11-13

!     Modifications.
!     --------------
!      F. Vana      : 03-08-25 new vertical profile of diffusion
!      G. Radnoti    : 03-09-23 GFL introduction
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      F. Vana and K. YESSAD: diffusion for SLHD + cleanings.
!      K. Yessad: (jul 2004) bound the increase of hor diffusion
!       in the upper stratosphere; ZPKDV renamed into Z_PDILEV for
!       consistency with suhdf.F90.
!      F. Vana      : 29-Oct-2004 cleaning + new setup of HD
!      F. Vana      : 08-Feb-2005 splitting of the global LSLHD key.
!      F. Vana      : 13-Jan-2009 removing argument + LSLHD_STATIC
!      F. Vana + NEC: 28-Apr-2009 OpenMP
!      K. Yessad (Jan 2012): shallow-water code removed, more flexibility in vert profiles.
!      A. Bogatchev 12-04-2013 phasing cy40, coherence with modified modules
!                              and renamed namelists and functions
!                              O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMLUN   , ONLY : NULOUT
USE YOMCT0   , ONLY : LNHDYN
USE YOMMP0   , ONLY : NPROC
USE YOMDYNA  , ONLY : YRDYNA
USE YOMDYN   , ONLY : TDYN
USE YEMDYN   , ONLY : TEDYN

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT) :: YDDYN
TYPE(TEDYN)    ,INTENT(INOUT) :: YDEDYN
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
REAL(KIND=JPRB) :: Z_PDILEV(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_PDILEV_SLD(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_PDILEVS(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZX, ZEFOLV, ZEFOLD
REAL(KIND=JPRB) :: ZFACX, ZFACY, ZFAC

INTEGER(KIND=JPIM) :: INSHDI(YDGEOMETRY%YRDIMV%NFLEVG),INSHDS(YDGEOMETRY%YRDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: INMSHDI(YDGEOMETRY%YRDIMV%NFLEVG),INMSHDS(YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IM, IN, IL, JLEV, JMLOC, JN, ISTDV2
INTEGER(KIND=JPIM) :: IGGX, IGGY, ISHIFTX, ISHIFTY

LOGICAL :: LLSLHD, LLSLHD_W ,LLSLHD_T ,LLSLHD_Q ,LLSLHD_O3 ,LLSLHD_SVD ,LLSLHD_SPD
LOGICAL :: LLMESO

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "suhdvp.intfb.h"
#include "suehdvpn.intfb.h"

!     ------------------------------------------------------------------

! * FUNCTIONS AND FUNCTION ARGUMENTS:
REAL(KIND=JPRB) :: BDISPE
INTEGER(KIND=JPIM) :: KNSMAX, KNMSMAX, KN, KM
REAL(KIND=JPRB) :: PEXPDH

! * HORIZONTAL DIFFUSION SHAPE FUNCTIONS:
BDISPE(KM,KN,KNMSMAX,KNSMAX,PEXPDH)=MAX(0.0_JPRB,( ( SQRT(&
 & REAL(KM,JPRB)**2/MAX(0.5_JPRB,REAL(KNMSMAX,JPRB))**2&
 & + REAL(KN,JPRB)**2/MAX(0.5_JPRB,REAL(KNSMAX,JPRB))**2)&
 &- YDDYN%FRANDH )/( 1.0_JPRB-YDDYN%FRANDH ) ) )**PEXPDH

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUEHDF',0,ZHOOK_HANDLE)

ASSOCIATE( &
 & LSPT=>YDML_GCONF%YRDIMF%LSPT, &
 & REXPDH=>YDDYN%REXPDH, &
 & HRDIRVOR=>YDDYN%HRDIRVOR, &
 & HRDIRDIV=>YDDYN%HRDIRDIV, &
 & HRDIRT=>YDDYN%HRDIRT, &
 & HRDIRQ=>YDDYN%HRDIRQ, &
 & HRDIRO3=>YDDYN%HRDIRO3, &
 & HRDIRPD=>YDDYN%HRDIRPD, &
 & HRDIRVD=>YDDYN%HRDIRVD, &
 & HRDIRSP=>YDDYN%HRDIRSP, &
 & FRANDH=>YDDYN%FRANDH, &
 & REXPDHS=>YDDYN%REXPDHS, &
 & HRDSRVOR=>YDDYN%HRDSRVOR, &
 & HRDSRDIV=>YDDYN%HRDSRDIV, &
 & HRDSRVD=>YDDYN%HRDSRVD, &
 & NPROFILEHD=>YDDYN%NPROFILEHD, &
 & SLEVDH1=>YDDYN%SLEVDH1, &
 & SLEVDHS1=>YDDYN%SLEVDHS1, &
 & SLEVDH2=>YDDYN%SLEVDH2, &
 & SLEVDHS2=>YDDYN%SLEVDHS2, &
 & LRDISPE_EC=>YDDYN%LRDISPE_EC, &
 & YQ=>YDML_GCONF%YGFL%YQ, YO3=>YDML_GCONF%YGFL%YO3, YQ_NL=>YDML_GCONF%YGFL%YQ_NL, YO3_NL=>YDML_GCONF%YGFL%YO3_NL, &
 & RDIVORE=>YDEDYN%RDIVORE, RDIDIVE=>YDEDYN%RDIDIVE, RDITE=>YDEDYN%RDITE, &
 & RDIGFLE=>YDEDYN%RDIGFLE, RDIPDE=>YDEDYN%RDIPDE, RDIVDE=>YDEDYN%RDIVDE, &
 & RDISPE=>YDEDYN%RDISPE, RDSVORE=>YDEDYN%RDSVORE, RDSDIVE=>YDEDYN%RDSDIVE, &
 & NSECPLG=>YDGEOMETRY%YREDIM%NSECPLG, &
 & NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,RDSVDE=>YDEDYN%RDSVDE,  &
 & YRLAP=>YDGEOMETRY%YRLAP, YRLEP=>YDGEOMETRY%YRELAP, &
 & NSMAX=>YDGEOMETRY%YRDIM%NSMAX, NMSMAX=>YDGEOMETRY%YRDIM%NMSMAX, NUMP=>YDGEOMETRY%YRDIM%NUMP)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY INITIALISATIONS AND VERTICAL PROFILES.
!              --------------------------------------------------

! * 1.1: Initialization of SLHD keys

LLSLHD    = (.NOT.YRDYNA%LSLHD_STATIC) .AND. YRDYNA%LSLHD
LLSLHD_W  = (.NOT.YRDYNA%LSLHD_STATIC) .AND. YRDYNA%LSLHD_W
LLSLHD_T  = (.NOT.YRDYNA%LSLHD_STATIC) .AND. YRDYNA%LSLHD_T
LLSLHD_Q  = (.NOT.YRDYNA%LSLHD_STATIC) .AND. YQ_NL%LSLHD
LLSLHD_O3 = (.NOT.YRDYNA%LSLHD_STATIC) .AND. YO3_NL%LSLHD
LLSLHD_SVD= (.NOT.YRDYNA%LSLHD_STATIC) .AND. YRDYNA%LSLHD_SVD
LLSLHD_SPD= (.NOT.YRDYNA%LSLHD_STATIC) .AND. YRDYNA%LSLHD_SPD

! * 1.2: Vertical profiles Z_PDILEV..

CALL SUHDVP(YDGEOMETRY,YDDYN,Z_PDILEV,Z_PDILEV_SLD,Z_PDILEVS,LLMESO,ISTDV2)

! * 1.3: Vertical profiles INSHD. and INMSHD..

WRITE(NULOUT,'(A)')
WRITE(NULOUT,'(A)') ' --- SUEHDVPN:'
WRITE(NULOUT,'(A)') ' VERTICAL PROFILES N(l) AND Nm(l) FOR HORIZONTAL DIFFUSION: HDI'
CALL SUEHDVPN(YDGEOMETRY,YDDYN,REXPDH,SLEVDH1,SLEVDH2,INSHDI,INMSHDI)
IF (LLSLHD) THEN
  WRITE(NULOUT,'(A)') ' VERTICAL PROFILES N(l) AND Nm(l) FOR HORIZONTAL DIFFUSION: HDS'
  CALL SUEHDVPN(YDGEOMETRY,YDDYN,REXPDHS,SLEVDHS1,SLEVDHS2,INSHDS,INMSHDS)
ELSE
  ! not used in this case.
  INSHDS(:)=NSMAX
  INMSHDS(:)=NMSMAX
ENDIF

!     ------------------------------------------------------------------

!*       2.    MAIN HORIZONTAL DIFFUSION.
!              --------------------------

IF(ALLOCATED(YDEDYN%RDIGFLE)) RDIGFLE(:,:,:) = 0.0_JPRB

!*       2.1   Computation of RDI[X]E at ECMWF:
!nyc IF (LRDISPE_EC) THEN
!nyc   not yet implemented in LAM models.
!nyc ENDIF

!*       2.2   Computation of RDI[X]E for METEO-FRANCE type set-up:

IF (.NOT.LRDISPE_EC) THEN

  IF (NPROFILEHD == 1 .OR. NPROFILEHD == 2 .OR. NPROFILEHD == 3) THEN

    !*       2.2.1 3D variables:
    IGGX=NINT(REAL(2*NMSMAX,JPRB)/3._JPRB)-1
    IGGY=NINT(REAL(2*NSMAX,JPRB)/3._JPRB)-1
    ISHIFTX=NINT(0.83_JPRB*REAL(NMSMAX,JPRB))-IGGX
    ISHIFTY=NINT(0.83_JPRB*REAL(NSMAX,JPRB))-IGGY
!$OMP PARALLEL DO &
!$OMP& PRIVATE(JLEV,JMLOC,JN,IM,IN,ZFACX,ZFACY,ZFAC)
    DO JLEV=1,NFLEVG

      DO JMLOC=1,NUMP
        IM=YRLAP%MYMS(JMLOC)

#ifdef NECSX
!cdir NOVECTOR
#else
!ocl noeval, nopreex
#endif
        DO JN=0,YRLEP%NCPLM(IM)-1

          IN=YRLEP%NPME(IM)+JN
          IF(YRDYNA%LGRADSP) THEN
            IF(IM + ISHIFTX > IGGX) THEN
              ZFACX = 32._JPRB*EXP(-0.5_JPRB*REAL(IM+ISHIFTX-NMSMAX,JPRB)**2/&
               &REAL(IM+ISHIFTX-IGGX,JPRB)**2)
            ELSE
              ZFACX = 0._JPRB
            ENDIF
            IF(JN+ISHIFTY > IGGY) THEN
              ZFACY = 32._JPRB*EXP(-0.5_JPRB*REAL(JN+ISHIFTY-NSMAX,JPRB)**2/&
                &REAL(JN+ISHIFTY-IGGY,JPRB)**2)
            ELSE
              ZFACY = 0._JPRB
            ENDIF
            ZFAC=(ZFACX*ZFACX+ZFACY*ZFACY)
            YDEDYN%REFILV(JLEV,IN)=1._JPRB/(1._JPRB+ZFAC)
            YDEDYN%REFILD(JLEV,IN)=1._JPRB/(1._JPRB+ZFAC/5000._JPRB)
          ENDIF

          IF(YRDYNA%LGRADSP) THEN
            IF(IM + ISHIFTX > IGGX) THEN
              ZFACX = 32._JPRB*EXP(-0.5_JPRB*REAL(IM+ISHIFTX-NMSMAX,JPRB)**2/REAL(IM+ISHIFTX-IGGX,JPRB)**2)
            ELSE
              ZFACX = 0._JPRB
            ENDIF
            IF(JN+ISHIFTY > IGGY) THEN
              ZFACY = 32._JPRB*EXP(-0.5_JPRB*REAL(JN+ISHIFTY-NSMAX,JPRB)**2/REAL(JN+ISHIFTY-IGGY,JPRB)**2)
            ELSE
              ZFACY = 0._JPRB
            ENDIF
            ZFAC=(ZFACX*ZFACX+ZFACY*ZFACY)
            YDEDYN%REFILV(JLEV,IN)=1._JPRB/(1._JPRB+ZFAC)
            YDEDYN%REFILD(JLEV,IN)=1._JPRB/(1._JPRB+ZFAC/5000._JPRB)
          ENDIF

          ! * vorticity
          IF (LLSLHD_W) THEN
            RDIVORE(JLEV,IN)=Z_PDILEV_SLD(JLEV)*BDISPE(IM,JN,INMSHDI(JLEV),INSHDI(JLEV),REXPDH)*HRDIRVOR
          ELSE
            RDIVORE(JLEV,IN)=Z_PDILEV(JLEV)*BDISPE(IM,JN,INMSHDI(JLEV),INSHDI(JLEV),REXPDH)*HRDIRVOR
          ENDIF

          ! * divergence:
          IF (LLSLHD_W) THEN
            RDIDIVE(JLEV,IN)=Z_PDILEV_SLD(JLEV)*BDISPE(IM,JN,INMSHDI(JLEV),INSHDI(JLEV),REXPDH)*HRDIRDIV
          ELSE
            RDIDIVE(JLEV,IN)=Z_PDILEV(JLEV)*BDISPE(IM,JN,INMSHDI(JLEV),INSHDI(JLEV),REXPDH)*HRDIRDIV
          ENDIF

          ! * temperature:
          IF (LSPT) THEN
            IF (LLSLHD_T) THEN
              RDITE(JLEV,IN)=Z_PDILEV_SLD(JLEV)*BDISPE(IM,JN,INMSHDI(JLEV),INSHDI(JLEV),REXPDH)*HRDIRT
            ELSE
              RDITE(JLEV,IN)=Z_PDILEV(JLEV)*BDISPE(IM,JN,INMSHDI(JLEV),INSHDI(JLEV),REXPDH)*HRDIRT
            ENDIF
          ENDIF

          ! * NH pressure departure variable:
          IF(LNHDYN) THEN
            IF (LLSLHD_SPD) THEN
              RDIPDE(JLEV,IN)=Z_PDILEV_SLD(JLEV)*BDISPE(IM,JN,INMSHDI(JLEV),INSHDI(JLEV),REXPDH)*HRDIRPD
            ELSE
              RDIPDE(JLEV,IN)=Z_PDILEV(JLEV)*BDISPE(IM,JN,INMSHDI(JLEV),INSHDI(JLEV),REXPDH)*HRDIRPD
            ENDIF
          ENDIF

          ! * NH vertical divergence variable:
          IF(LNHDYN) THEN
            IF (LLSLHD_SVD) THEN
              RDIVDE(JLEV,IN)=Z_PDILEV_SLD(JLEV)*BDISPE(IM,JN,INMSHDI(JLEV),INSHDI(JLEV),REXPDH)*HRDIRVD
            ELSE
              RDIVDE(JLEV,IN)=Z_PDILEV(JLEV)*BDISPE(IM,JN,INMSHDI(JLEV),INSHDI(JLEV),REXPDH)*HRDIRVD
            ENDIF
          ENDIF

          ! -- GFL variables:

          ! * specific humidity:
          IF (YQ%LSP) THEN
            IF (LLSLHD_Q) THEN
              RDIGFLE(JLEV,IN,YQ%MPSP)=Z_PDILEV_SLD(JLEV)*BDISPE(IM,JN,INMSHDI(JLEV),INSHDI(JLEV),REXPDH)*HRDIRQ
            ELSE
              RDIGFLE(JLEV,IN,YQ%MPSP)=Z_PDILEV(JLEV)*BDISPE(IM,JN,INMSHDI(JLEV),INSHDI(JLEV),REXPDH)*HRDIRQ
            ENDIF
          ENDIF

          ! * ozone:
          IF (YO3%LSP) THEN
            IF (LLSLHD_O3) THEN
              RDIGFLE(JLEV,IN,YO3%MPSP)=Z_PDILEV_SLD(JLEV)*BDISPE(IM,JN,INMSHDI(JLEV),INSHDI(JLEV),REXPDH)*HRDIRO3
            ELSE
              RDIGFLE(JLEV,IN,YO3%MPSP)=Z_PDILEV(JLEV)*BDISPE(IM,JN,INMSHDI(JLEV),INSHDI(JLEV),REXPDH)*HRDIRO3
            ENDIF
          ENDIF

          ! * other GFL variables: they are currently not diffused.

          ! -- end of GFL variables:

        ENDDO ! JN
      ENDDO ! JMLOC
    ENDDO ! JLEV
!$OMP END PARALLEL DO

    !*       2.2.2 2D variables:

    ! Remark: no diffusion is generally done on this variable,
    ! so the SLHD scheme has not been implemented.

    DO JMLOC=1,NUMP
      IM=YRLAP%MYMS(JMLOC)
#ifdef NECSX
!cdir NOVECTOR
#else
!ocl noeval, nopreex
#endif
      DO JN=0,YRLEP%NCPLM(IM)-1
        IN=YRLEP%NPME(IM)+JN
        RDISPE(IN)=BDISPE(IM,JN,NMSMAX,NSMAX,REXPDH)*HRDIRSP
      ENDDO ! JN
    ENDDO ! JMLOC

  ELSEIF (NPROFILEHD == 4) THEN

    ! not implemented in LAM models.

  ENDIF ! NPROFILEHD

ENDIF ! .NOT.LRDISPE_EC

!*       2.3   Computation of RDS[X]E for VOR,DIV,VD for METEO-FRANCE type set-up:

IF ((.NOT.LRDISPE_EC).AND.(LLSLHD_W.OR.LLSLHD_SVD)) THEN

  IF (NPROFILEHD == 1 .OR. NPROFILEHD == 2 .OR. NPROFILEHD == 3) THEN

!$OMP PARALLEL DO &
!$OMP& PRIVATE(JLEV,JMLOC,JN,IM,IN)
    DO JLEV=1,NFLEVG
    
      DO JMLOC=1,NUMP
        IM=YRLAP%MYMS(JMLOC)

#ifdef NECSX
!cdir NOVECTOR
#else 
!ocl noeval, nopreex
#endif
        DO JN=0,YRLEP%NCPLM(IM)-1

          IN=YRLEP%NPME(IM)+JN

          ! * vorticity and divergence:
          IF (LLSLHD_W) THEN
            RDSVORE(JLEV,IN)=Z_PDILEVS(JLEV)*BDISPE(IM,JN,INMSHDS(JLEV),INSHDS(JLEV),REXPDHS)*HRDSRVOR
            RDSDIVE(JLEV,IN)=Z_PDILEVS(JLEV)*BDISPE(IM,JN,INMSHDS(JLEV),INSHDS(JLEV),REXPDHS)*HRDSRDIV
          ENDIF

          ! * NH vertical divergence variable:
          IF(LNHDYN.AND.LLSLHD_SVD) THEN
            RDSVDE(JLEV,IN)=Z_PDILEVS(JLEV)*BDISPE(IM,JN,INMSHDS(JLEV),INSHDS(JLEV),REXPDHS)*HRDSRVD
          ENDIF

        ENDDO ! JN
      ENDDO ! JMLOC
    ENDDO ! JLEV
!$OMP END PARALLEL DO

  ELSEIF (NPROFILEHD == 4) THEN

    ! not implemented in LAM models but obsolescent.

  ENDIF ! NPROFILEHD

ENDIF ! .NOT.LRDISPE_EC

!*       2.4   Diagnostic of e-folding time:

! Printings are done for (m,n) giving the maximum for RDI[X].
! Printings are done for NPROC=1 only because RDI[X] is DM-local (no DM-global RDI[X]G available).
IF (NPROC == 1) THEN
  WRITE(NULOUT,*) ' '
  DO JLEV=NFLEVG,1,-1
    IL=JLEV
    IF (IL <= ISTDV2) THEN
      ZX=MAXVAL(RDIVORE(JLEV,1:NSECPLG))
      IF(ZX > 1E-7_JPRB) THEN
        ZEFOLV=1.0_JPRB/(ZX*3600._JPRB)
      ELSE
        ZEFOLV=0._JPRB
      ENDIF
      ZX=MAXVAL(RDIDIVE(JLEV,1:NSECPLG))
      IF(ZX > 1E-7_JPRB) THEN
        ZEFOLD=1.0_JPRB/(ZX*3600._JPRB)
      ELSE
        ZEFOLD=0._JPRB
      ENDIF
      WRITE(NULOUT,*) '  E-FOLDING TIME FOR VORTICITY=',ZEFOLV,' AND FOR DIVERGENCE=',ZEFOLD,' HOURS AT LEVEL ',IL
    ENDIF
  ENDDO
  WRITE(NULOUT,*) ' '
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUEHDF',1,ZHOOK_HANDLE)
END SUBROUTINE SUEHDF

