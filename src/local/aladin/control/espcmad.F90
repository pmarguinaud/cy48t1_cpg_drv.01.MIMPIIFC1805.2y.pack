SUBROUTINE ESPCMAD(YDGEOMETRY,YDMODEL,CDCONF,YDSP)

!**** *ESPCMAD* - Interface to spectral calculations in LAM models (AD code)

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *ESPCMAD

!        Explicit arguments :
!        --------------------
!          CDCONF      - configuration of work (see doc.)
!                        CDCONF='0': nothing done
!                        CDCONF='A': SI scheme and horizontal diffusion done.
!                        CDCONF='I': SI scheme done; no horizontal diffusion.

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------  

!     Reference.
!     ----------
!        ARPEGE/ALADIN documentation

!     Author.
!     -------
!      Andras Horanyi CNRM/GMAP/EXT
!      Original : 95-07-27

!     Modifications.
!     --------------
!      P. Smolikova 02-09-30 : interface to ETRMTOS,ETRSTOM for d4 in NH
!      A. Bogatchev 03-09-24 : introducing GFL structure
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      A. Bogatchev 04-10-07  removed interafce for aladin routines
!      K.Yessad     : Feb 2005: split ESPCAD into ESPCHORAD,ESPCSIAD,ESPNHSIAD
!      M.Janousek   : Sep 2005: replace ETRSTOM/ETRMTOS by TRSTOM/TRMTOS
!      K. YESSAD: (Aug 2005). SI scheme in ALADIN done like in ARPEGE.
!      K. Yessad (Feb 2012): simplifications for tests and CDCONF.
!      B. Bochenek 11-04-2013 - Phasing cy40, coherence with modified modules
!      O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     ------------------------------------------------------------------

USE TYPE_MODEL , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCT0   , ONLY : LR3D     ,LNHDYN
USE YOMLUN   , ONLY : NULOUT
USE YOMDYNA  , ONLY : NVDVAR   ,LGWADV
USE SPECTRAL_FIELDS_DATA, ONLY: SPECTRAL_FIELD


!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT)    :: YDGEOMETRY
TYPE(MODEL)     ,INTENT(INOUT)   :: YDMODEL
CHARACTER(LEN=1)  ,INTENT(IN) :: CDCONF 
TYPE(SPECTRAL_FIELD), INTENT(INOUT) :: YDSP

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: INH , INHX
INTEGER(KIND=JPIM) :: IM

LOGICAL :: LLESIDG
LOGICAL :: LLONEM
LOGICAL :: LLDOSI, LLDOHD

REAL(KIND=JPRB) ::  ZSPVORG(YDGEOMETRY%YRDIMV%NFLEVG,MAX(YDGEOMETRY%YRMP%NSPEC2V,YDGEOMETRY%YRMP%NSPEC2VF))
REAL(KIND=JPRB) ::  ZSPDIVG(YDGEOMETRY%YRDIMV%NFLEVG,MAX(YDGEOMETRY%YRMP%NSPEC2V,YDGEOMETRY%YRMP%NSPEC2VF))
REAL(KIND=JPRB) ::  ZSPTG  (YDGEOMETRY%YRDIMV%NFLEVG,MAX(YDGEOMETRY%YRMP%NSPEC2V,YDGEOMETRY%YRMP%NSPEC2VF))
REAL(KIND=JPRB) ::  ZSPSPG (MAX(YDGEOMETRY%YRMP%NSPEC2V,YDGEOMETRY%YRMP%NSPEC2VF))
REAL(KIND=JPRB) ::  ZSPGFLG(YDGEOMETRY%YRDIMV%NFLEVG,MAX(YDGEOMETRY%YRMP%NSPEC2V,YDGEOMETRY%YRMP%NSPEC2VF), &
                          & MAX(1,YDMODEL%YRML_GCONF%YGFL%NUMSPFLDS))

REAL(KIND=JPRB), ALLOCATABLE :: ZSPSPDG(:,:), ZSPSVDG(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZSPNHXG(:,:)

REAL(KIND=JPRB),ALLOCATABLE :: ZZPS(:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "espchorad.intfb.h"
#include "espcsiad.intfb.h"
#include "trmtos.intfb.h"
#include "trstom.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ESPCMAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDLAP=>YDGEOMETRY%YRLAP, YDELAP=>YDGEOMETRY%YRELAP, &
 & YDDYN=>YDMODEL%YRML_DYN%YRDYN, YDEDYN=>YDMODEL%YRML_DYN%YREDYN, YDML_GCONF=>YDMODEL%YRML_GCONF)
ASSOCIATE( &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NSPEC2V=>YDMP%NSPEC2V, NSPEC2VF=>YDMP%NSPEC2VF, &
 & LIMPF=>YDDYN%LIMPF, NSPEC2=>YDDIM%NSPEC2)
!     ------------------------------------------------------------------

!*       0.    CONTROL ALLOCATIONS
!              -------------------

ALLOCATE(ZZPS(NSPEC2))
IF(YDMP%NPSP==1)ZZPS=YDSP%SP

IF (LNHDYN) THEN
  INH=1
  IF (NVDVAR == 4 .AND. .NOT.LGWADV) THEN
    INHX=1
  ELSE
    INHX=0
  ENDIF
ELSE
  INH=0
  INHX=0
ENDIF

!     ------------------------------------------------------------------

!*       1.    LOOP OVER ZONAL WAVENUMBER M.
!              -----------------------------

! future LESIDG: like LSIDG but for ALADIN in the case where:
!  - large domain where the variations of the mapping factor M are large.
!  - type of projections where M is quasi constant on a latitude and
!    M can be written as a low order polynomial of mu=sin(latitude).
LLESIDG=.FALSE.
                                                                                
IF (LLESIDG.OR.LIMPF) THEN
  LLONEM=.TRUE.
ELSE
  LLONEM=.FALSE.
ENDIF

IF (CDCONF == 'P') THEN

  !       * FILTERING OF FULL-POS SPECTRAL FIELDS.

  WRITE (NULOUT,*) 'NO ADJOINT FULL-POS AVAILABLE'
  CALL ABOR1('ESPCMAD: NO ADJOINT POST-PROCESSING')

ELSEIF (CDCONF == 'A'.OR.CDCONF == 'I') THEN

 IF (LR3D) THEN

  ALLOCATE(ZSPSPDG(NFLEVG,MAX(NSPEC2V,NSPEC2VF)*INH))
  ALLOCATE(ZSPSVDG(NFLEVG,MAX(NSPEC2V,NSPEC2VF)*INH))
  ALLOCATE(ZSPNHXG(NFLEVG,MAX(NSPEC2V,NSPEC2VF)*INHX))

  IF (LLONEM) THEN
    ! llonem not implemented so no false illusions are permitted
    WRITE (NULOUT,*) 'M PER M TREATMENT NOT IMPLEMENTED FOR ',&
     & 'NPROC > 1  IN ALADIN'
    CALL ABOR1('ESPCMAD: LLONEM=.T. IMPOSSIBLE')
  ELSE
    ! dummy argument KM not used in ESPCSIAD.
    IM=-999
  ENDIF

  CALL GSTATS(86,0)

  ! ky: this call to TRMTOS has to be moved later between ESPCSIAD
  !     and ESPCHORAD to be consistent with what is done in SPCMAD.
  IF (YDML_GCONF%YGFL%NUMSPFLDS > 0) THEN
    CALL TRMTOS(YDGEOMETRY,PSPVOR=YDSP%VOR,PSPDIV=YDSP%DIV,PSPT=YDSP%T,&
     & PSPSPD=YDSP%SPD,PSPSVD=YDSP%SVD,PSPSNHX=YDSP%NHX,&
     & PSPGFL=YDSP%GFL,PSPSP=ZZPS,&
     & PSPVORG=ZSPVORG,PSPDIVG=ZSPDIVG,PSPTG=ZSPTG,&
     & PSPSPDG=ZSPSPDG,PSPSVDG=ZSPSVDG,PSPSNHXG=ZSPNHXG,&
     & PSPGFLG=ZSPGFLG,PSPSPG=ZSPSPG)
  ELSE
    CALL TRMTOS(YDGEOMETRY,PSPVOR=YDSP%VOR,PSPDIV=YDSP%DIV,PSPT=YDSP%T,&
     & PSPSPD=YDSP%SPD,PSPSVD=YDSP%SVD,PSPSNHX=YDSP%NHX,&
     & PSPSP=ZZPS,&
     & PSPVORG=ZSPVORG,PSPDIVG=ZSPDIVG,PSPTG=ZSPTG,&
     & PSPSPDG=ZSPSPDG,PSPSVDG=ZSPSVDG,PSPSNHXG=ZSPNHXG,&
     & PSPSPG=ZSPSPG)
  ENDIF

  CALL GSTATS(86,2)

  CALL GSTATS(36,0)

  ! * HORIZONTAL DIFFUSION.
  LLDOHD=(CDCONF == 'A')

  !   All spectral columns are done in one go
  IF (LLDOHD .AND. NSPEC2V > 0) THEN
    CALL ESPCHORAD(YDGEOMETRY,YDMODEL,1,NSPEC2V,LLONEM,&
     & ZSPVORG,ZSPDIVG,ZSPTG,ZSPGFLG,ZSPSPG,ZSPSPDG,ZSPSVDG)
  ENDIF

  ! * SPECTRAL NUDGING.
  !   no spectral nudging in the adjoint code.

  ! * SEMI-IMPLICIT SCMEME.
  LLDOSI=(CDCONF == 'A'.OR.CDCONF == 'I')

  !   All spectral columns are done in one go
  IF (LLDOSI .AND. NSPEC2V > 0) THEN
    IF (LNHDYN) THEN
      CALL ABOR1('ESPCMAD: Adjoint of NH model not coded')
    ELSE 
      CALL ESPCSIAD(YDGEOMETRY,YDDYN,YDML_GCONF%YRRIP,IM,1,NSPEC2V,LLONEM,ZSPVORG,ZSPDIVG,ZSPTG,ZSPSPG)
    ENDIF
  ENDIF

  CALL GSTATS(36,2)

  CALL GSTATS(86,3)

  IF (YDML_GCONF%YGFL%NUMSPFLDS > 0) THEN
    CALL TRSTOM(YDGEOMETRY,PSPVORG=ZSPVORG,PSPDIVG=ZSPDIVG,PSPTG=ZSPTG,&
     & PSPSPDG=ZSPSPDG,PSPSVDG=ZSPSVDG,PSPSNHXG=ZSPNHXG,&
     & PSPGFLG=ZSPGFLG,PSPSPG=ZSPSPG,&
     & PSPVOR=YDSP%VOR,PSPDIV=YDSP%DIV,PSPT=YDSP%T,&
     & PSPSPD=YDSP%SPD,PSPSVD=YDSP%SVD,PSPSNHX=YDSP%NHX,&
     & PSPGFL=YDSP%GFL,PSPSP=ZZPS)
  ELSE
    CALL TRSTOM(YDGEOMETRY,PSPVORG=ZSPVORG,PSPDIVG=ZSPDIVG,PSPTG=ZSPTG,&
     & PSPSPDG=ZSPSPDG,PSPSVDG=ZSPSVDG,PSPSNHXG=ZSPNHXG,&
     & PSPSPG=ZSPSPG,&
     & PSPVOR=YDSP%VOR,PSPDIV=YDSP%DIV,PSPT=YDSP%T,&
     & PSPSPD=YDSP%SPD,PSPSVD=YDSP%SVD,PSPSNHX=YDSP%NHX,&
     & PSPSP=ZZPS)
  ENDIF

  CALL GSTATS(86,1)
  CALL GSTATS(36,3)

  CALL GSTATS(36,1)

  IF (ALLOCATED(ZSPSPDG)) DEALLOCATE(ZSPSPDG)
  IF (ALLOCATED(ZSPSVDG)) DEALLOCATE(ZSPSVDG)
  IF (ALLOCATED(ZSPNHXG)) DEALLOCATE(ZSPNHXG)

 ELSE

  CALL ABOR1('ESPCMAD MUST NOT BE CALLED FOR THIS VALUE OF NCONF')

 ENDIF ! LR3D

ELSE

  IF (CDCONF /= '0') THEN
    WRITE (NULOUT,*) 'CDCONF(9:9) = ', CDCONF
    CALL ABOR1('ESPCMAD: CONFIGURATION NOT ALLOWED: CHECK CDCONF')
  ENDIF

ENDIF ! CDCONF

IF(YDMP%NPSP==1)YDSP%SP=ZZPS
DEALLOCATE(ZZPS)
!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ESPCMAD',1,ZHOOK_HANDLE)
END SUBROUTINE ESPCMAD
