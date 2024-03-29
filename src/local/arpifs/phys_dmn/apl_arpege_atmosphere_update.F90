SUBROUTINE APL_ARPEGE_ATMOSPHERE_UPDATE (YDMF_PHYS_BASE_STATE, YDMF_PHYS_NEXT_STATE, YDGEOMETRY,      &
& YDCPG_BNDS, YDCPG_OPTS, YDCPG_MISC, YDMF_PHYS, YDCPG_DYN0, YDMF_PHYS_SURF, YDVARS, YDSURF, YDMODEL, &
& LDCONFX, PDTPHY, YDDDH, PDIFEXT, PDSA_CPS, PMSC_FRMQ, PPFL_FEFB1, PPFL_FEFB2, PPFL_FEFB3,           &
& PPFL_FTKE, PPFL_FTKEI)

USE PARKIND1, ONLY : JPIM, JPRB

USE MF_PHYS_BASE_STATE_TYPE_MOD &
                       , ONLY : MF_PHYS_BASE_STATE_TYPE
USE MF_PHYS_NEXT_STATE_TYPE_MOD &
                       , ONLY : MF_PHYS_NEXT_STATE_TYPE
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE CPG_OPTS_TYPE_MOD   , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_MISC_TYPE, CPG_DYN_TYPE, &
                              & CPG_SL2_TYPE, CPG_GPAR_TYPE
USE MF_PHYS_TYPE_MOD   , ONLY : MF_PHYS_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_MISC_TYPE, CPG_DYN_TYPE, &
                              & CPG_SL2_TYPE, CPG_GPAR_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD,ONLY : MF_PHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE TYPE_MODEL         , ONLY : MODEL
USE DDH_MIX            , ONLY : TYP_DDH
USE YOMCT0             , ONLY : LTWOTL


IMPLICIT NONE

TYPE (MF_PHYS_BASE_STATE_TYPE), INTENT(IN)    :: YDMF_PHYS_BASE_STATE
TYPE (MF_PHYS_NEXT_STATE_TYPE), INTENT(INOUT) :: YDMF_PHYS_NEXT_STATE
TYPE(GEOMETRY),                 INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_BNDS_TYPE),            INTENT(IN)    :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE),            INTENT(IN)    :: YDCPG_OPTS
TYPE(CPG_MISC_TYPE),            INTENT(INOUT) :: YDCPG_MISC
TYPE(MF_PHYS_TYPE),             INTENT(INOUT) :: YDMF_PHYS
TYPE(CPG_DYN_TYPE),             INTENT(IN)    :: YDCPG_DYN0
TYPE(MF_PHYS_SURF_TYPE),        INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES),          INTENT(INOUT) :: YDVARS
TYPE(TSURF),                    INTENT(IN)    :: YDSURF
TYPE(MODEL),                    INTENT(IN)    :: YDMODEL
LOGICAL,                        INTENT(IN)    :: LDCONFX
REAL(KIND=JPRB),                INTENT(IN)    :: PDTPHY
TYPE(TYP_DDH),                  INTENT(INOUT) :: YDDDH
REAL(KIND=JPRB),                INTENT(IN)    :: PDIFEXT(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB),                INTENT(IN)    :: PDSA_CPS (YDCPG_OPTS%KLON)
REAL(KIND=JPRB),                INTENT(IN)    :: PMSC_FRMQ (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(IN)    :: PPFL_FEFB1 (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(IN)    :: PPFL_FEFB2 (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(IN)    :: PPFL_FEFB3 (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(IN)    :: PPFL_FTKE (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(IN)    :: PPFL_FTKEI (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)

#include "cpmvvps.intfb.h"
#include "cptend_new.intfb.h"
#include "cputqy_aplpar_expl.intfb.h"

REAL(KIND=JPRB) :: ZMSC_FHP (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZPFL_FP (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDD (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDEFB1  (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDEFB2  (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDEFB3  (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDEXT(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB) :: ZTENDG     (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDH(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDI     (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDICONV (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDL     (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDLCONV (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDQ     (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDR     (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDRCONV (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDS     (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDSCONV (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDTKE   (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDU (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTENDV (YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)

INTEGER (KIND=JPIM) :: JLEV
INTEGER (KIND=JPIM) :: JROF


!        2.7  Computation of tendencies T,u,v and Q.
!             --------------------------------------

! Set GFL tendencies to 0

ZTENDEFB1   = 0._JPRB
ZTENDEFB2   = 0._JPRB
ZTENDEFB3   = 0._JPRB
ZTENDG      = 0._JPRB
ZTENDI      = 0._JPRB
ZTENDQ      = 0._JPRB
ZTENDR      = 0._JPRB
ZTENDS      = 0._JPRB
ZTENDL      = 0._JPRB
ZTENDICONV  = 0._JPRB
ZTENDLCONV  = 0._JPRB
ZTENDRCONV  = 0._JPRB
ZTENDSCONV  = 0._JPRB
ZTENDTKE    = 0._JPRB
! ky: non-zero option not yet coded for the time being.
ZTENDD      = 0.0_JPRB


! * CPTEND+CPUTQY = Old( CPATY + CPDUP + CPDTHP )
! Calcul des tendances de T , U et de Q et modifications
! eventuelles de W et de OMEGA/P


CALL CPTEND_NEW( YDMODEL, YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, YDVARS%GEOMETRY%GNORDL%T0, &
& YDVARS%GEOMETRY%GNORDM%T0, YDMF_PHYS%OUT%DIFCQ, YDMF_PHYS%OUT%DIFCQN, YDMF_PHYS%OUT%DIFCQL, YDMF_PHYS%OUT%DIFCS,           &
& PDIFEXT, YDMF_PHYS%OUT%DIFTQ, YDMF_PHYS%OUT%DIFTQN, YDMF_PHYS%OUT%DIFTQL, YDMF_PHYS%OUT%DIFTS, YDMF_PHYS%OUT%FCCQL,        &
& YDMF_PHYS%OUT%FCCQN, YDMF_PHYS%OUT%FCSQL, YDMF_PHYS%OUT%FCSQN, YDMF_PHYS%OUT%FPLSL, YDMF_PHYS%OUT%FPLSN,                   &
& YDMF_PHYS%OUT%FPLSG, YDMF_PHYS%OUT%FPLCL, YDMF_PHYS%OUT%FPLCN, YDMF_PHYS%OUT%FPLCG, YDMF_PHYS%OUT%FPEVPSL,                 &
& YDMF_PHYS%OUT%FPEVPSN, YDMF_PHYS%OUT%FPEVPSG, YDMF_PHYS%OUT%FPEVPCL, YDMF_PHYS%OUT%FPEVPCN, YDMF_PHYS%OUT%FPEVPCG,         &
& YDMF_PHYS%OUT%FPFPSL, YDMF_PHYS%OUT%FPFPSN, YDMF_PHYS%OUT%FPFPSG, YDMF_PHYS%OUT%FPFPCL, YDMF_PHYS%OUT%FPFPCN,              &
& YDMF_PHYS%OUT%FCQLNG, YDMF_PHYS%OUT%FCQNNG, YDMF_PHYS%OUT%FCQRNG, YDMF_PHYS%OUT%FCQSNG, YDMF_PHYS%OUT%FCQGNG,              &
& YDMF_PHYS%OUT%FCQNG, YDMF_PHYS%OUT%FRMH, PMSC_FRMQ, YDMF_PHYS%OUT%FRSO, YDMF_PHYS%OUT%FRTH, YDMF_PHYS%OUT%STRCU,           &
& YDMF_PHYS%OUT%STRCV, YDMF_PHYS%OUT%STRDU, YDMF_PHYS%OUT%STRDV, YDMF_PHYS%OUT%STRTU, YDMF_PHYS%OUT%STRTV,                   &
& YDMF_PHYS%OUT%STRMU, YDMF_PHYS%OUT%STRMV, YDMF_PHYS%OUT%DIFCQLC, YDMF_PHYS%OUT%DIFCQIC, YDMF_PHYS%OUT%FIMCC,               &
& YDMF_PHYS%OUT%FEDQLC, YDMF_PHYS%OUT%FEDQIC, YDMF_PHYS%OUT%FEDQRC, YDMF_PHYS%OUT%FEDQSC, YDMF_PHYS%OUT%FCNEGQLC,            &
& YDMF_PHYS%OUT%FCNEGQIC, YDMF_PHYS%OUT%FCNEGQRC, YDMF_PHYS%OUT%FCNEGQSC, PPFL_FTKE, PPFL_FTKEI,                             &
& PPFL_FEFB1, PPFL_FEFB2, PPFL_FEFB3, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%RDELP,       &
& YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V,  &
& YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%Q, YDMF_PHYS_BASE_STATE%I, YDMF_PHYS_BASE_STATE%L, YDVARS%LCONV%T0,           &
& YDVARS%ICONV%T0, YDVARS%RCONV%T0, YDVARS%SCONV%T0, YDMF_PHYS_BASE_STATE%R, YDMF_PHYS_BASE_STATE%S,                         &
& YDMF_PHYS_BASE_STATE%G, PDSA_CPS, YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDMF_PHYS%OUT%FHSCL, YDMF_PHYS%OUT%FHSCN,                &
& YDMF_PHYS%OUT%FHSSL, YDMF_PHYS%OUT%FHSSN, YDMF_PHYS%OUT%FHSSG, YDMF_PHYS%OUT%FHPCL, YDMF_PHYS%OUT%FHPCN,                   &
& YDMF_PHYS%OUT%FHPCG, YDMF_PHYS%OUT%FHPSL, YDMF_PHYS%OUT%FHPSN, YDMF_PHYS%OUT%FHPSG, ZMSC_FHP,                              &
& ZPFL_FP, YDMF_PHYS%OUT%FEPFP, YDMF_PHYS%OUT%FCMPCQ, YDMF_PHYS%OUT%FCMPSN, YDMF_PHYS%OUT%FCMPSL, YDMF_PHYS%OUT%TENDU,       &
& YDMF_PHYS%OUT%TENDV, ZTENDU, ZTENDV, ZTENDH, ZTENDQ, ZTENDI, ZTENDL, ZTENDLCONV, ZTENDICONV,                               &
& ZTENDRCONV, ZTENDSCONV, ZTENDR, ZTENDS, ZTENDG, ZTENDTKE, ZTENDEFB1, ZTENDEFB2, ZTENDEFB3,                                 &
& ZTENDEXT, YDDDH)


IF (LTWOTL) THEN
ELSE    
  IF ( (YDMODEL%YRML_PHY_MF%YRPHY%NDPSFI==1)) THEN
!     PFEPFP was ZFEPFP in CPTEND_NEW, before, ZFEPFP still in CPFHPFS
    DO JLEV= 0, YDCPG_OPTS%KFLEVG 
      DO JROF = 1, YDCPG_OPTS%KLON
        YDMF_PHYS%OUT%FEPFP (JROF,JLEV) = 0.0_JPRB
        YDMF_PHYS%OUT%FCMPCQ(JROF,JLEV) = 0.0_JPRB
        YDMF_PHYS%OUT%FCMPSN(JROF,JLEV) = 0.0_JPRB
        YDMF_PHYS%OUT%FCMPSL(JROF,JLEV) = 0.0_JPRB
      ENDDO
    ENDDO
  ENDIF
ENDIF 

!        2.8  Modification of vertical velocities
!             by some physics output when required.
!             -------------------------------------


! * MODIFICATION DE LA VITESSE VERTICALE ET DE LA TENDANCE DE
! PRESSION DE SURFACE SI NDPSFI=1 ( MASSE VARIABLE ).
! Ajout de la physique dans l'equation de continuite/Add physics
! in continuity equation.

IF (YDMODEL%YRML_PHY_MF%YRPHY%NDPSFI == 1) THEN
  CALL CPMVVPS(YDGEOMETRY%YRVAB, YDCPG_OPTS%KLON, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KFLEVG, &
  & PDTPHY, ZPFL_FP, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD(:, YDCPG_OPTS%KFLEVG), YDMF_PHYS%OUT%FEVL,     &
  & YDMF_PHYS%OUT%FEVN, YDCPG_DYN0%CTY%EVEL, YDCPG_DYN0%CTY%PSDVBC, YDMF_PHYS_NEXT_STATE%SP)
ENDIF

!        2.9  Computation of evolution of T, u, v and Q.
!             ------------------------------------------

! Calcul de T , Q et du Vent a l'instant 1

CALL CPUTQY_APLPAR_EXPL(YDCPG_BNDS, YDCPG_OPTS, YDMF_PHYS_NEXT_STATE, YDMF_PHYS_BASE_STATE, YDVARS, YDMODEL%YRML_PHY_MF%YRPHY, &
& PDTPHY, ZTENDH, YDMF_PHYS%OUT%TENDU, YDMF_PHYS%OUT%TENDV, ZTENDU, ZTENDV, ZTENDD, ZTENDEFB1,                                 &
& ZTENDEFB2, ZTENDEFB3, ZTENDG, ZTENDICONV, ZTENDI, ZTENDLCONV, ZTENDL, ZTENDQ, ZTENDRCONV, ZTENDR,                            &
& ZTENDSCONV, ZTENDS, ZTENDTKE, YDMF_PHYS%OUT%FDIS)



END SUBROUTINE
