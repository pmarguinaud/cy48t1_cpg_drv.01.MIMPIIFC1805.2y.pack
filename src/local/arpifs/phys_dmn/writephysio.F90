SUBROUTINE WRITEPHYSIO(YDGEOMETRY,&
 & YDCPG_MISC,YDCPG_PHY0,YDCPG_PHY9,YDMF_PHYS,YDCPG_DYN0,YDCPG_DYN9,YDMF_PHYS_SURF,YDVARS,&
 & YDSURF,YDDPHY,YDRIP,YDML_PHY_MF,KEND, KST, KGL1, KGL2, KSTGLO,&
 & KSTEP  , KSGST  , KCSS   , PGELAM, PGEMU,&
 & PGM, POROG, PRCORI, PRATATH, PRATATX,&
 & PCLON  , PSLON, YDMF_PHYS_TMP, YDAPLPAR_TMP )

! --------------------------------------------------------------
! **** *writephysio* Write SCM input/output of physics on file.
! --------------------------------------------------------------
! Subject:
!     Write on files all inputs and some outputs from MF_PHYS
!     for profiles matching a given criteria.

! Method:
! Authors:  1999-03, J.M. Piriou / F. Bouyssel.
! Modifications:
!   2000-10, J.M. Piriou: use of a namelist for SCM extractions.
!   2002-03, J.M. Piriou: same arguments for writephysio as for MF_PHYS.
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!   2004-03: T. Kovacic : extracted subroutine WRITEPROFILE
!   2004-10, F. Bouyssel : cleaning
!   2007-05, F. Bouyssel : add PFCLL, PFCLN
!   2011-08, E. Bazile : add fields for WRITEMUSC (old recopie).
!   2013-11, Y. Bouteloup : EFB fields
!   K. Yessad (July 2014): Move some variables.
! --------------------------------------------------------------

USE MF_PHYS_TYPE_MOD   , ONLY : MF_PHYS_TYPE, APLPAR_TMP_TYPE, MF_PHYS_TMP_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_DYN_TYPE, CPG_PHY_TYPE, &
                              & CPG_MISC_TYPE
USE MF_PHYS_SURFACE_TYPE_MOD,ONLY : MF_PHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1  ,ONLY : JPIM    ,JPRB
USE YOMHOOK   ,ONLY : LHOOK   ,DR_HOOK
USE YOMLUN    ,ONLY : NULOUT
USE YOMRIP    ,ONLY : TRIP
USE YOMSCM    ,ONLY : NFRSCM, NSCMTS, NSCM_SPACE_S, GSCM_LON1, GSCM_LON2&
 & , GSCM_LAT1, GSCM_LAT2, GSCM_RADIUS, NSCM_ADD_SAMPL , LGSCM 
USE YOMDPHY   ,ONLY : TDPHY
USE YOMCST    ,ONLY : RA  
USE YOMCT3    ,ONLY : NSTEP
USE YOMLSFORC, ONLY : LMUSCLFA

! --------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(CPG_MISC_TYPE),INTENT(INOUT):: YDCPG_MISC
TYPE(CPG_PHY_TYPE),INTENT(INOUT) :: YDCPG_PHY0
TYPE(CPG_PHY_TYPE),INTENT(INOUT) :: YDCPG_PHY9
TYPE(MF_PHYS_TYPE),INTENT(INOUT) :: YDMF_PHYS
TYPE(CPG_DYN_TYPE),INTENT(INOUT) :: YDCPG_DYN0
TYPE(CPG_DYN_TYPE),INTENT(INOUT) :: YDCPG_DYN9
TYPE(MF_PHYS_SURF_TYPE),INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
TYPE(TDPHY)       ,INTENT(INOUT) :: YDDPHY
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
TYPE(TRIP)        ,INTENT(INOUT) :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGL1
INTEGER(KIND=JPIM),INTENT(IN)    :: KGL2
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTGLO
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP,  KSGST, KCSS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGM(YDGEOMETRY%YRDIM%NPROMA)

REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCORI(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRATATH(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRATATX(YDGEOMETRY%YRDIM%NPROMA)


REAL(KIND=JPRB),INTENT(IN)    :: PCLON(YDGEOMETRY%YRDIM%NPROMM)
REAL(KIND=JPRB),INTENT(IN)    :: PSLON(YDGEOMETRY%YRDIM%NPROMM)
TYPE(MF_PHYS_TMP_TYPE),INTENT(IN) :: YDMF_PHYS_TMP
TYPE(APLPAR_TMP_TYPE), INTENT(IN),OPTIONAL :: YDAPLPAR_TMP

! --------------------------------------------------------------
INTEGER(KIND=JPIM) :: JSTEP,JLON
LOGICAL :: LLOK, LLTIME_OK
REAL(KIND=JPRB) :: ZLAT_MODELE,ZLON_MODELE,ZD,ZHOUR,ZRR,ZPI,ZREQUIRED
REAL(KIND=JPRB) :: ZHOOK_HANDLE

! --------------------------------------------------------------

#include "abor1.intfb.h"
#include "writeprofile.intfb.h"
#include "writemusc.intfb.h"

! --------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WRITEPHYSIO', 0, ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
  & YDMP=>YDGEOMETRY%YRMP, YDVAB=>YDGEOMETRY%YRVAB, YDVETA=>YDGEOMETRY%YRVETA, YDVFE=>YDGEOMETRY%YRVFE)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, NPROMM=>YDDIM%NPROMM, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NTSSG=>YDDPHY%NTSSG, &
 & RSTATI=>YDRIP%RSTATI, TSTEP=>YDRIP%TSTEP)

! --------------------------------------------------------------

IF (LGSCM) THEN
! -------------------------------------------------
! Define criteria for choosing profiles.
! -------------------------------------------------

!-------------------------------------------------
! 1. Time sampling.
!-------------------------------------------------
ZHOUR=RSTATI/3600._JPRB ! time since the beginning of the prediction in hours.
IF(NSCMTS(0) == 0) THEN
  !
  ! -------------------------------------------------
  ! One time step over nfrscm.
  ! -------------------------------------------------
  !
  LLTIME_OK=MOD(NSTEP,NFRSCM) == 0
  WRITE(NULOUT,*) 'writephysio: nscmts(0) == 0: lltime_ok=',LLTIME_OK
ELSEIF(NSCMTS(0) > 0) THEN
  !
  ! -------------------------------------------------
  ! Explicit list of time steps.
  ! -------------------------------------------------
  !
  LLTIME_OK=.FALSE.
  DO JSTEP=1,NSCMTS(0)
    IF(NFRSCM*NSCMTS(JSTEP) == NSTEP) LLTIME_OK=.TRUE.
  ENDDO
  WRITE(NULOUT,*) 'writephysio: nscmts(0) > 0: NSTEP=',NSTEP,' lltime_ok=',LLTIME_OK
ELSE
  !
  ! -------------------------------------------------
  ! Explicit list of hours.
  ! -------------------------------------------------
  !
  LLTIME_OK=.FALSE.
  DO JSTEP=1,-NSCMTS(0)
    ZREQUIRED=-REAL(NFRSCM,JPRB)*NSCMTS(JSTEP) ! required time in hours.
    IF(ABS(ZHOUR-ZREQUIRED) <= TSTEP) LLTIME_OK=.TRUE.
  ENDDO
  WRITE(NULOUT,*) 'writephysio: expl list hours: lltime_ok=',LLTIME_OK
ENDIF

!-------------------------------------------------
! 2. Space sampling.
!-------------------------------------------------

! RAPPEL: en t042: 8000 pdg, en T199: 180000 pdg.

IF (LLTIME_OK) THEN
  DO JLON=KST,KEND
    !
    ! -------------------------------------------------
    ! Compute model grid-point lat/lon coordinates.
    ! -------------------------------------------------
    !
    ZLAT_MODELE=ASIN(PGEMU(JLON))
    ZPI=4._JPRB*ATAN(1.0_JPRB)
    ZLON_MODELE=MODULO(PGELAM(JLON)+ZPI,2.0_JPRB*ZPI)-ZPI
    !
    ! -------------------------------------------------
    ! Is this profile required by the user?
    ! -------------------------------------------------
    !
    LLOK=.FALSE.
    IF(NSCM_SPACE_S == 0) THEN
      !
      ! -------------------------------------------------
      ! Choose profiles close to a given location.
      ! -------------------------------------------------
      !
      ! Compute othodromic distance (in meters) between searched point and model grid-point.
      !
      ZD=COS(ZLAT_MODELE)*COS(GSCM_LAT1)*COS(ZLON_MODELE-GSCM_LON1)&
       & + SIN(ZLAT_MODELE)*SIN(GSCM_LAT1)  
      ZD=ACOS(ZD)*RA
      LLOK=ZD < GSCM_RADIUS

    ELSEIF(NSCM_SPACE_S == 1) THEN
      !
      ! -------------------------------------------------
      ! Choose profiles inside a given box (in (lon, lat)),
      ! with additional regular sampling inside that box.
      ! -------------------------------------------------
      !
      LLOK=ZLON_MODELE > GSCM_LON1&
       & .AND. ZLON_MODELE < GSCM_LON2&
       & .AND. ZLAT_MODELE > GSCM_LAT1&
       & .AND. ZLAT_MODELE < GSCM_LAT2&
       & .AND.MOD(JLON,NSCM_ADD_SAMPL) == 0  
    ELSEIF(NSCM_SPACE_S == 10) THEN
      !
      ! -------------------------------------------------
      ! Choose profiles depending on physical output.
      ! -------------------------------------------------
      !
      ZRR=YDMF_PHYS%OUT%FPLCL(JLON,NFLEVG)+YDMF_PHYS%OUT%FPLCN(JLON,NFLEVG)
      LLOK=ZRR > 0.0_JPRB
    ELSEIF(NSCM_SPACE_S < 0) THEN
      !
      ! -------------------------------------------------
      ! Choose one profile over -NSCM_SPACE_S.
      ! -------------------------------------------------
      !
      LLOK=MOD(JLON,-NSCM_SPACE_S) == 0
    ELSE
      WRITE (NULOUT,*) 'writephysio/ERROR: NSCM_SPACE_S unexpected!...'
      WRITE (NULOUT,*) NSCM_SPACE_S
      CALL ABOR1('writephysio/ERROR: NSCM_SPACE_S unexpected!')
    ENDIF
    IF (LLOK) THEN
      !
      ! -------------------------------------------------
      ! The grid-point is a desired profile.
      ! -------------------------------------------------
      !
      CALL WRITEPROFILE(YDGEOMETRY, &
      & YDCPG_MISC,YDMF_PHYS,YDMF_PHYS_TMP,YDCPG_DYN0,YDMF_PHYS_SURF,YDVARS,&
      & YDSURF, YDDPHY, YDRIP, YDML_PHY_MF, 'Profile.', JLON,               &
      & ZLAT_MODELE, ZLON_MODELE, ZD, PGELAM, PGEMU, PGM, POROG, &
      & PRCORI, PRATATH, PRATATX, &
      & YDAPLPAR_TMP)
    ENDIF
  ENDDO
ENDIF

ENDIF ! (LGSCM)

IF (LMUSCLFA) THEN
   CALL WRITEMUSC(&
   & YDCPG_MISC,YDCPG_PHY0,YDCPG_PHY9,YDMF_PHYS,YDCPG_DYN0,YDCPG_DYN9,YDMF_PHYS_SURF,YDVARS,   &
   & YDVETA, YDRIP, YDML_PHY_MF, KST, KEND, NPROMA, 1, NFLEVG, KGL1, KGL2, KSTEP, KSGST,       &
   & KCSS, PCLON, PSLON, PGELAM, PGEMU, PGM, POROG, 'WRITEMUSC', YDMF_PHYS_TMP, YDAPLPAR_TMP)

ENDIF

! --------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('WRITEPHYSIO', 1, ZHOOK_HANDLE)
END SUBROUTINE WRITEPHYSIO
