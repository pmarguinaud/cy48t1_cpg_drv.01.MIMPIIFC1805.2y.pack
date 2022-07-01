SUBROUTINE PPNEW(YDGP5,PRESO,KCOUNT,PXPP,CDPAR,KGRIBCODE,KHORIZ,KSLCT,&
  & LDEXTLI,LDHEIGHT,LDNEUTRAL,LDRELCURR,LDGEMS_PWC,PES5OUT,YDGP_TL,YDGP_AD)

! Vertical interpolation and conversion of model profiles into simple geophysical variables.
! This is a replacement for the old "PPOBSA" or "PPOBSAZ". Done as part of OOPS observation 
! cleaning. But note the lower-level routines are only partly refactored.

!     AUTHOR.
!     -------
!      A. GEER          ECMWF        30-SEP-2015    refactor PPOBSA and other originals

!     MODIFICATIONS.
!     --------------
!      C. Payan    01-DEC-2015 add calls to tl/ad
!      B. Ingleby  18-Feb-2016 Use Sonntag option for forward RH calculation
!      A. Geer     28-Jun-2016 Protect against missing CDPAR and IGRIBCODE not found
!      B. Ingleby  09-Oct-2019 Print CDPAR if not recognised, return early if Q2

!     ------------------------------------------------------------------

USE PARKIND1    , ONLY : JPIM, JPRB
USE YOMHOOK     , ONLY : LHOOK, DR_HOOK
USE GOM_PLUS    , ONLY : TYPE_GOM_PLUS, IH
USE INTDYN_MOD  , ONLY : YYTXYB
USE YOMLUN      , ONLY : NULERR
USE YOMANCS     , ONLY : RMDI

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TYPE_GOM_PLUS),INTENT(IN)                      :: YDGP5           ! Model profiles at obs locations
REAL(KIND=JPRB)    ,INTENT(IN)                      :: PRESO(:,:)      ! Output pressure/height levels (IDLEN,IDCOUNT) [Pa?/m?]
INTEGER(KIND=JPIM) ,INTENT(IN)                      :: KCOUNT(:)       ! Maximum count per obs (IDLEN)
REAL(KIND=JPRB)    ,INTENT(INOUT), TARGET           :: PXPP(:,:,:)   ! Processed ouput (IDLEN,IDCOUNT,IDPP)
INTEGER(KIND=JPIM) ,INTENT(IN)   , OPTIONAL         :: KGRIBCODE ! Specify output variable: GRIB code - atmos composition only
CHARACTER(LEN=*)   ,INTENT(IN)   , OPTIONAL         :: CDPAR  ! Specify output variable: old "CLV" character code that will be replaced
INTEGER(KIND=JPIM) ,INTENT(IN)   , OPTIONAL         :: KHORIZ ! For 2D gom_plus, specify the location to extract (default is KHORIZ=IH=1)
INTEGER(KIND=JPIM) ,INTENT(IN)   , OPTIONAL         :: KSLCT 
LOGICAL            ,INTENT(IN)   , OPTIONAL         :: LDEXTLI 
LOGICAL            ,INTENT(IN)   , OPTIONAL         :: LDHEIGHT ! Vertical interpolation in height rather than pressure, which is default
LOGICAL            ,INTENT(IN)   , OPTIONAL         :: LDNEUTRAL, LDRELCURR ! Scatterometer surface wind options
LOGICAL            ,INTENT(IN)   , OPTIONAL         :: LDGEMS_PWC ! Composition variables: total column (true) or vertical profile (false or not present)
REAL(KIND=JPRB)    ,INTENT(OUT)  , OPTIONAL, TARGET :: PES5OUT(:) 
TYPE(TYPE_GOM_PLUS),INTENT(IN)   , OPTIONAL         :: YDGP_TL   ! model variables - TL
TYPE(TYPE_GOM_PLUS),INTENT(INOUT), OPTIONAL         :: YDGP_AD   ! model variables - adjoint

INTEGER(KIND=JPIM) :: IDLEN, IDCOUNT, IDPP, IVAR, INVAR, IGEMS
INTEGER(KIND=JPIM) :: IMXCNT, JLEVP, JC, JOBS, KH
! Extra 10 is a zone to allow us to more safely identify inadequate length in dimension IDPP
!CP is this really safe? At least dim(IVARLIST)=2*size(PXPP,3), u/v case
INTEGER(KIND=JPIM) :: IVARLIST(SIZE(PXPP,3)+10)
REAL(KIND=JPRB)    :: ZLNPRESO(YDGP5%NDLEN,SIZE(PRESO,2))
INTEGER(KIND=JPIM) :: ILEVB(YDGP5%NDLEN,SIZE(PRESO,2),YDGP5%NPPM)
REAL(KIND=JPRB), DIMENSION(YDGP5%NDLEN,YDGP5%NFLEVG) :: ZR05,ZR0
LOGICAL :: LLBELO(YDGP5%NDLEN,SIZE(PRESO,2))
LOGICAL :: LLBELS(YDGP5%NDLEN,SIZE(PRESO,2))
LOGICAL :: LLBLOW(SIZE(PRESO,2))
LOGICAL :: LLBLES(SIZE(PRESO,2))
LOGICAL :: LLPRES, LLEXTLI, LL_SPEED, LL_NONLINEAR
LOGICAL :: LLNEUTRAL,LLRELCURR
LOGICAL :: LLGEMS_PWC
LOGICAL :: LLTL, LLAD, LLDIRECT
LOGICAL :: LLGRIBFOUND, LLPRESO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL(KIND=JPRB), PARAMETER :: ZRHMAX=1.2_JPRB, ZRHMIN=-0.2_JPRB
REAL(KIND=JPRB), POINTER :: ZFOR_INTERP5(:,:), ZFOR_INTERP(:,:)
REAL(KIND=JPRB), TARGET  :: ZES5(YDGP5%NDLEN)
REAL(KIND=JPRB), POINTER :: ZESOUT5(:)
REAL(KIND=JPRB), POINTER :: ZXPP5(:,:,:)
REAL(KIND=JPRB) :: ZSPEED5

TYPE TYPE_PPVAR
  INTEGER(KIND=JPIM) :: U=1
  INTEGER(KIND=JPIM) :: V=2
  INTEGER(KIND=JPIM) :: T=3
  INTEGER(KIND=JPIM) :: H=4
  INTEGER(KIND=JPIM) :: Q=5
  INTEGER(KIND=JPIM) :: PWC=6
  INTEGER(KIND=JPIM) :: Z=7
  INTEGER(KIND=JPIM) :: O3=8
  INTEGER(KIND=JPIM) :: TO3=9
  INTEGER(KIND=JPIM) :: CC=10
  INTEGER(KIND=JPIM) :: CLW=11
  INTEGER(KIND=JPIM) :: CIW=12
  INTEGER(KIND=JPIM) :: TCC=13
  INTEGER(KIND=JPIM) :: U10=14
  INTEGER(KIND=JPIM) :: V10=15
  INTEGER(KIND=JPIM) :: T2=16
  INTEGER(KIND=JPIM) :: H2=17
  INTEGER(KIND=JPIM) :: PS=18
  INTEGER(KIND=JPIM) :: COMPOSITION=19
END TYPE
TYPE(TYPE_PPVAR) :: PPVAR

#include "abor1.intfb.h"
!from ppobsa (altitude):
#include "gprcp_qlirsg.intfb.h"
#include "gprcptl.intfb.h"
#include "gprcpad.intfb.h"
#include "ppcc.intfb.h"
#include "ppcctl.intfb.h"
#include "ppccad.intfb.h"
#include "ppclw.intfb.h"
#include "ppclwtl.intfb.h"
#include "ppclwad.intfb.h"
#include "ppflev.intfb.h"
#include "ppgeop.intfb.h"
#include "ppgeoptl.intfb.h"
#include "ppgeopad.intfb.h"
#include "pppwc.intfb.h"
#include "pppwctl.intfb.h"
#include "pppwcad.intfb.h"
#include "ppq.intfb.h"
#include "ppqtl.intfb.h"
#include "ppqad.intfb.h"
#include "pprh.intfb.h"
#include "pprhtl.intfb.h"
#include "pprhad.intfb.h"
#include "ppt.intfb.h"
#include "ppttl.intfb.h"
#include "pptad.intfb.h"
#include "pptcc.intfb.h"
#include "pptcctl.intfb.h"
!#include "pptccad.intfb.h"
#include "ppuv.intfb.h"
#include "ppuvtl.intfb.h"
#include "ppuvad.intfb.h"
!from ppobsaz (pressure):
#include "ppps.intfb.h"
#include "pppstl.intfb.h"
#include "pppsad.intfb.h"
#include "ppzhlev.intfb.h"
!from ppobsas (surface):
#include "pprh2m.intfb.h"
#include "pprh2mtl.intfb.h"
#include "pprh2mad.intfb.h"
#include "ppt2m.intfb.h"
#include "ppt2mtl.intfb.h"
#include "ppt2mad.intfb.h"
#include "ppuv10m.intfb.h"
#include "ppuv10mtl.intfb.h"
#include "ppuv10mad.intfb.h"

IF (LHOOK) CALL DR_HOOK('PPNEW',0,ZHOOK_HANDLE)

LLTL=PRESENT(YDGP_TL)
LLAD=PRESENT(YDGP_AD)
LLDIRECT=.NOT.(LLTL.OR.LLAD)

IDLEN   = YDGP5%NDLEN
IDCOUNT = SIZE(PRESO,2)
IDPP    = SIZE(PXPP,3) 

IF(ANY(SHAPE(PRESO)  /= (/IDLEN,IDCOUNT/))) CALL ABOR1('PPNEW: PRESO WRONG SHAPE')
IF(ANY(SHAPE(KCOUNT) /= (/IDLEN/))) CALL ABOR1('PPNEW: KCOUNT WRONG SHAPE')
IF(ANY(SHAPE(PXPP)   /= (/IDLEN,IDCOUNT,IDPP/))) CALL ABOR1('PPNEW: PXPP WRONG SHAPE')

IMXCNT = MAXVAL(KCOUNT)
IF(IMXCNT == 0) THEN 
  IF (LHOOK) CALL DR_HOOK('PPNEW',1,ZHOOK_HANDLE)
  RETURN
ENDIF

IF(PRESENT(KHORIZ))THEN
  ! Location in a 2D gom
  KH=KHORIZ
ELSE
  ! Normal 1D goms (IH=1)
  KH=IH
ENDIF

LLGEMS_PWC=.FALSE.
IF(PRESENT(LDGEMS_PWC)) LLGEMS_PWC=LDGEMS_PWC
LLEXTLI=.FALSE.
IF(PRESENT(LDEXTLI)) LLEXTLI=LDEXTLI
LLPRES=.TRUE.
IF(PRESENT(LDHEIGHT)) LLPRES=.NOT.LDHEIGHT
IF(.NOT.LLPRES) THEN
  IF(.NOT.PRESENT(KSLCT)) CALL ABOR1('PPNEW IN HEIGHT MODE: SELECT INTERPOLATION TYPE')
ENDIF
LLNEUTRAL=.FALSE.
IF(PRESENT(LDNEUTRAL)) LLNEUTRAL=LDNEUTRAL
LLRELCURR=.FALSE.
IF(PRESENT(LDRELCURR)) LLRELCURR=LDRELCURR
LL_SPEED=.FALSE.
LL_NONLINEAR=.FALSE.
LL_SPEED= CDPAR=='FFS' .OR. CDPAR=='FF'
LL_NONLINEAR = LL_SPEED
IF(PRESENT(CDPAR)) LL_NONLINEAR = LL_NONLINEAR .OR. (CDPAR=='H2')

! Create list of non-constituent variables required. 
! AJGDB it would be good to eliminate these character codes and just use varnos.
! AJGDB to create the longer lists of variables, the calling routine should be
! responsible (e.g. BG1, TVC, etc.)
IVARLIST = -1
IF(PRESENT(CDPAR)) THEN

  IF(CDPAR=='Q2') THEN  ! Q2m needs L_SCREEN_LEVEL_OBSOP=true, newer obs operator
    IF (LHOOK) CALL DR_HOOK('PPNEW',1,ZHOOK_HANDLE)
    RETURN
  ENDIF

  IF(CDPAR=='PWC') THEN
    IVARLIST(1:1)=(/PPVAR%PWC/)
  ELSEIF(CDPAR=='PS') THEN  
    IVARLIST(1:1)=(/PPVAR%PS/)
  ELSEIF(CDPAR=='T') THEN  
    IVARLIST(1:1)=(/PPVAR%T/)
  ELSEIF(CDPAR=='Q') THEN  
    IVARLIST(1:1)=(/PPVAR%Q/)
  ELSEIF(CDPAR=='TO3') THEN  
    IVARLIST(1:1)=(/PPVAR%TO3/)
  ELSEIF(CDPAR=='CC') THEN  
    IVARLIST(1:1)=(/PPVAR%CC/)
  ELSEIF(CDPAR=='TCC') THEN  
    IVARLIST(1:1)=(/PPVAR%TCC/)
  ELSEIF(CDPAR=='CLW') THEN  
    IVARLIST(1:1)=(/PPVAR%CLW/)
  ELSEIF(CDPAR=='CIW') THEN  
    IVARLIST(1:1)=(/PPVAR%CIW/)
  ELSEIF(CDPAR=='T2') THEN  
    IVARLIST(1:1)=(/PPVAR%T2/)
  ELSEIF(CDPAR=='H2') THEN  
    IVARLIST(1:1)=(/PPVAR%H2/)
  ELSEIF(CDPAR=='TV1') THEN  
    IVARLIST(1:3)=(/PPVAR%T,PPVAR%Q,PPVAR%O3/)
  ELSEIF(CDPAR=='RTL') THEN
    IVARLIST(1:4)=(/PPVAR%T,PPVAR%Q,PPVAR%O3,PPVAR%Z/)
  ELSEIF(CDPAR=='TVC') THEN  
    IVARLIST(1:6)=(/PPVAR%T,PPVAR%Q,PPVAR%O3,PPVAR%CC,PPVAR%CLW,PPVAR%CIW/)
  ELSEIF(CDPAR=='TQL') THEN  
    IVARLIST(1:4)=(/PPVAR%T,PPVAR%Q,PPVAR%O3,PPVAR%CLW/)
  ELSEIF(CDPAR=='BG1') THEN  
    IVARLIST(1:9)=(/PPVAR%U,PPVAR%V,PPVAR%T,PPVAR%H,PPVAR%Q,PPVAR%PWC,PPVAR%Z,PPVAR%O3,PPVAR%TO3/)
  ELSEIF(CDPAR=='H' .OR. CDPAR=='LH') THEN  
    IVARLIST(1:1)=(/PPVAR%H/)
  ELSEIF(CDPAR=='Z' .OR. CDPAR=='DZ') THEN  
    IVARLIST(1:1)=(/PPVAR%Z/)
  ELSEIF(CDPAR=='U10' .OR. CDPAR=='FFS') THEN  
    IVARLIST(1:2)=(/PPVAR%U10,PPVAR%V10/)
  ELSEIF(CDPAR=='U' .OR. CDPAR=='FF' .OR. CDPAR=='HLS') THEN  
    IVARLIST(1:2)=(/PPVAR%U,PPVAR%V/)
  ENDIF

ELSEIF(PRESENT(KGRIBCODE)) THEN

  IF(KGRIBCODE==-1) THEN
    ! Extract all composition variables
    IVARLIST(1:YDGP5%NGEMS)=PPVAR%COMPOSITION
  ELSE
    ! Extract a single composition variable by its GRIB code
    IVARLIST(1:1)=(/PPVAR%COMPOSITION/)
  ENDIF

ELSE
  CALL ABOR1('SPECIFY EITHER CDPAR OR KGRIBCODE')
ENDIF

INVAR = COUNT(IVARLIST/=-1)
IF(INVAR > IDPP) CALL ABOR1('PPNEW: PP VARIABLE DIMENSION INADEQUATE')
IF(INVAR == 0) CALL ABOR1('PPNEW: NO RECOGNISED VARIABLE SPECIFIED, CDPAR='//CDPAR)

IF(LLPRES) THEN 

  ! Find pressure level below specified (observed) pressures
  CALL PPFLEV(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,YDGP5%NPPM,&
    & PRESO,YDGP5%PRESH(:,:,KH),YDGP5%PRESF(:,:,KH),&
    & ILEVB,LLBELO,LLBELS,LLBLOW,LLBLES)  

  ! ln(p)
  LLPRESO=.FALSE.
  IF(PRESENT(CDPAR)) THEN
    IF(CDPAR == 'TV1' .OR. CDPAR == 'TVC' .OR. CDPAR == 'TQL') THEN
      LLPRESO=.TRUE.
    ENDIF
  ENDIF

  IF(LLPRESO) THEN
    DO JLEVP=1,IMXCNT
      ZLNPRESO(:,JLEVP) = LOG(PRESO(1,JLEVP))
    ENDDO
  ELSE
    DO JLEVP=1,IMXCNT
      ZLNPRESO(:,JLEVP) = LOG(PRESO(:,JLEVP))
    ENDDO
  ENDIF

ELSE

  ! Find height level below specified (observed) heights
  CALL PPZHLEV(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,YDGP5%NPPM,&
  & PRESO,YDGP5%GEOPH(:,:,KH),YDGP5%GEOPF(:,:,KH),&
  & ILEVB,LLBELO,LLBELS,LLBLOW,LLBLES)  
  
ENDIF


IF ((LLTL.OR.LLAD).AND.LL_NONLINEAR) THEN
  ! To allow non-linear operators to generate a trajectory observation equivalent
  ALLOCATE(ZXPP5(IDLEN,IDCOUNT,IDPP))
  ZXPP5=RMDI !CP : or 0.0_JPRB as ZXPP?
ELSE
  ! For linear case in operators that are sometimes non-linear
  ZXPP5=>PXPP
ENDIF

IF( ANY(IVARLIST == PPVAR%T) .OR. ANY(IVARLIST == PPVAR%Z) ) THEN
  ! Note: only the model knows "YGFL%YQ%LTHERMACT" so here we always assume it is true.
  ! To do better, we would need to bring that information here through GOM and GOM_PLUS
  CALL GPRCP_QLIRSG(IDLEN,1,IDLEN,YDGP5%NFLEVG,PQ=YDGP5%QF(:,1:YDGP5%NFLEVG,KH),PR=ZR05,LDTHERMACT=.TRUE.)

  IF (LLTL) THEN
    CALL GPRCPTL(IDLEN,1,IDLEN,YDGP5%NFLEVG,PQ=YDGP_TL%QF(:,1:YDGP5%NFLEVG,KH),PR=ZR0,LDTHERMACT=.TRUE.)
  ENDIF
  IF (LLAD) ZR0=0.0_JPRB
ENDIF

! AJGDB most of these lower-level routines could be rationalised into just a single one
! with options. Those old routines are another example of "copy and paste" coding,
! and can be refactored when time allows.
DO IVAR = 1,INVAR

  !CP some times initialization of one array only is required, or not required (surface), function of ivar value
  !CP   but is it really useful to distinct these different cases? Or to repeat this initialization
  !CP   in each if block (as it is/was done in ppobsaad)?

  IF(IVARLIST(IVAR) == PPVAR%U .OR. IVARLIST(IVAR)==PPVAR%U10) THEN
    IF(IVARLIST(IVAR) == PPVAR%U) THEN

      ! Note that "V" is also covered here and ignored when it comes up next
      IF (LLDIRECT) THEN
        CALL PPUV(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
          & LLBELO,LLBLOW,ZLNPRESO,YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),&
          & YDGP5%UF(:,:,KH),YDGP5%VF(:,:,KH),PXPP(:,:,IVAR),PXPP(:,:,IVAR+1))

      ELSEIF (LLTL) THEN
        CALL PPUVTL(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
          & LLBELO, LLBLOW,ZLNPRESO,YDGP_TL%PXP(:,:,:,KH),YDGP_TL%PXPD(:,:,:,KH),&
          & YDGP_TL%UF(:,:,KH),YDGP_TL%VF(:,:,KH),PXPP(:,:,IVAR),PXPP(:,:,IVAR+1),&
          & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%UF(:,:,KH),YDGP5%VF(:,:,KH))

      ELSEIF (LLAD) THEN
        YDGP_AD%UF(:,:,KH)=0.0_JPRB
        YDGP_AD%VF(:,:,KH)=0.0_JPRB

        CALL PPUVAD(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
          & LLBELO,LLBLOW,ZLNPRESO,YDGP_AD%PXP(:,:,:,KH),YDGP_AD%PXPD(:,:,:,KH),&
          & YDGP_AD%UF(:,:,KH),YDGP_AD%VF(:,:,KH),PXPP(:,:,IVAR),PXPP(:,:,IVAR+1),&
          & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%UF(:,:,KH),YDGP5%VF(:,:,KH))

      ENDIF

    ELSEIF(IVARLIST(IVAR)==PPVAR%U10) THEN

      ! Note that "V10" is also covered here and ignored when it comes up next
  
      IF(LLDIRECT .OR. LL_SPEED)THEN
        CALL PPUV10M(IDLEN,1,IDLEN,LLNEUTRAL,LLRELCURR,PRESO(:,1),&
        & YDGP5%UF(:,YDGP5%NFLEVG,KH),YDGP5%VF(:,YDGP5%NFLEVG,KH),YDGP5%APHI(:,YDGP5%NFLEVG,KH),&
        & YDGP5%RI(:,KH),YDGP5%BN(:,KH),YDGP5%BM(:,KH),YDGP5%U0(:,KH),YDGP5%V0(:,KH),&
        & ZXPP5(:,:,IVAR),ZXPP5(:,:,IVAR+1))
        DO JC=2,IMXCNT
          WHERE(KCOUNT(:) >= JC)
            ZXPP5(:,JC,IVAR)=ZXPP5(:,1,IVAR)
            ZXPP5(:,JC,IVAR+1)=ZXPP5(:,1,IVAR+1)
          ENDWHERE
        ENDDO
      ENDIF

      IF (LLTL) THEN
        CALL PPUV10MTL(IDLEN,1,IDLEN,LLNEUTRAL,PRESO(:,1),&
          & YDGP5%UF(:,YDGP5%NFLEVG,KH),YDGP5%VF(:,YDGP5%NFLEVG,KH),YDGP5%APHI(:,YDGP5%NFLEVG,KH),&
          & YDGP_TL%UF(:,YDGP5%NFLEVG,KH),YDGP_TL%VF(:,YDGP5%NFLEVG,KH),YDGP_TL%APHI(:,YDGP5%NFLEVG,KH),&
          & YDGP5%RI(:,KH),YDGP5%BN(:,KH),YDGP5%BM(:,KH),&
          & YDGP_TL%BN(:,KH),YDGP_TL%BM(:,KH),YDGP5%U0(:,KH),YDGP5%V0(:,KH),&
          & PXPP(:,:,IVAR),PXPP(:,:,IVAR+1))

      ELSEIF (LLAD) THEN
        IF (LL_SPEED) THEN
          DO JC=1,IMXCNT
            DO JOBS=1,IDLEN
              IF(JC > KCOUNT(JOBS)) CYCLE
              ZSPEED5 = SQRT(ZXPP5(JOBS,JC,IVAR)**2+ZXPP5(JOBS,JC,IVAR+1)**2)
              IF (ZSPEED5 >= 1.E-10_JPRB) THEN
                PXPP(JOBS,JC,IVAR) = ZXPP5(JOBS,JC,1)/ZSPEED5*PXPP(JOBS,JC,IVAR)
                PXPP(JOBS,JC,IVAR+1) = ZXPP5(JOBS,JC,2)/ZSPEED5*PXPP(JOBS,JC,IVAR+1)
              ENDIF
            ENDDO
          ENDDO
        ENDIF

        CALL PPUV10MAD(IDLEN,1,IDLEN,LLNEUTRAL,PRESO(:,1),&
          & YDGP5%UF(:,YDGP5%NFLEVG,KH),YDGP5%VF(:,YDGP5%NFLEVG,KH),YDGP5%APHI(:,YDGP5%NFLEVG,KH),&
          & YDGP_AD%UF(:,YDGP5%NFLEVG,KH),YDGP_AD%VF(:,YDGP5%NFLEVG,KH),YDGP_AD%APHI(:,YDGP5%NFLEVG,KH),&
          & YDGP5%RI(:,KH),YDGP5%BN(:,KH),YDGP5%BM(:,KH),&
          & YDGP_AD%BN(:,KH),YDGP_AD%BM(:,KH),YDGP5%U0(:,KH),YDGP5%V0(:,KH),&
          & PXPP(:,:,IVAR),PXPP(:,:,IVAR+1))

      ENDIF

    ENDIF

    IF (LL_SPEED.AND..NOT.LLAD) THEN
      IF (LLDIRECT) THEN
        DO JC=1,IMXCNT
          DO JOBS=1,IDLEN
            IF (JC > KCOUNT(JOBS)) CYCLE
            ZSPEED5 = SQRT(PXPP(JOBS,JC,IVAR)**2+PXPP(JOBS,JC,IVAR+1)**2)
            PXPP(JOBS,JC,IVAR) = MAX(ZSPEED5,1.E-10_JPRB)
            PXPP(JOBS,JC,IVAR+1) = RMDI
          ENDDO
        ENDDO
      ELSEIF (LLTL) THEN
        DO JC=1,IMXCNT
          DO JOBS=1,IDLEN
            IF (JC > KCOUNT(JOBS)) CYCLE
            ZSPEED5 = SQRT(ZXPP5(JOBS,JC,IVAR)**2+ZXPP5(JOBS,JC,IVAR+1)**2)
            IF (ZSPEED5 < 1.E-10_JPRB) THEN
              PXPP(JOBS,JC,IVAR) = 0.0_JPRB
            ELSE
                PXPP(JOBS,JC,IVAR) = (&
                & ZXPP5(JOBS,JC,IVAR)*PXPP(JOBS,JC,IVAR) +&
                & ZXPP5(JOBS,JC,IVAR+1)*PXPP(JOBS,JC,IVAR+1)   )/ZSPEED5
            ENDIF
            PXPP(JOBS,JC,IVAR+1) = RMDI
          ENDDO
        ENDDO
      ENDIF
    ENDIF

  ELSEIF(IVARLIST(IVAR) == PPVAR%T) THEN

    IF (LLDIRECT) THEN
      CALL PPT(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & PRESO,ZLNPRESO,LLBELO,LLBELS,LLBLOW,LLBLES,.TRUE.,YDGP5%OROG(:,KH),&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%TSTAR(:,KH),YDGP5%T0(:,KH),ZR05,YDGP5%TF(:,:,KH),&
        & PXPP(:,:,IVAR))

    ELSEIF (LLTL) THEN
      CALL PPTTL(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & PRESO,ZLNPRESO,LLBELO,LLBELS,LLBLOW,LLBLES,.TRUE.,YDGP5%OROG(:,KH),&
        & YDGP_TL%PXP(:,:,:,KH),YDGP_TL%PXPD(:,:,:,KH),YDGP_TL%TSTAR(:,KH),YDGP_TL%T0(:,KH),ZR0,YDGP_TL%TF(:,:,KH),&
        & PXPP(:,:,IVAR),&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%TSTAR(:,KH),YDGP5%T0(:,KH),ZR05,YDGP5%TF(:,:,KH))

    ELSEIF (LLAD) THEN
      CALL  PPTAD(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & PRESO,ZLNPRESO,LLBELO,LLBELS,LLBLOW,LLBLES,.TRUE.,YDGP5%OROG(:,KH),&
        & YDGP_AD%PXP(:,:,:,KH),YDGP_AD%PXPD(:,:,:,KH),YDGP_AD%TSTAR(:,KH),YDGP_AD%T0(:,KH),ZR0,YDGP_AD%TF(:,:,KH),&
        & PXPP(:,:,IVAR),&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%TSTAR(:,KH),YDGP5%T0(:,KH),ZR05,YDGP5%TF(:,:,KH))

    ENDIF

  ELSEIF(IVARLIST(IVAR) == PPVAR%H) THEN

    IF (LLDIRECT) THEN
      CALL PPRH(.TRUE.,IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & PRESO,LLBELO,LLBLOW,.TRUE.,ZRHMAX,ZRHMIN,&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%QF(:,:,KH),YDGP5%TF(:,:,KH),&
        & PXPP(:,:,IVAR))
      ! LDSONNTAG=.TRUE. for use in calculation of background RH for Sonde and 
      ! Aircraft, no TL/AD provided/used since q is analysis variable

    ELSEIF (LLTL) THEN
      CALL PPRHTL(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & PRESO,LLBELO,LLBLOW,ZRHMAX,ZRHMIN,&
        & YDGP_TL%PXP(:,:,:,KH),YDGP_TL%PXPD(:,:,:,KH),YDGP_TL%QF(:,:,KH),YDGP_TL%TF(:,:,KH),&
        & PXPP(:,:,IVAR),&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%QF(:,:,KH),YDGP5%TF(:,:,KH))

    ELSEIF (LLAD) THEN
      CALL PPRHAD(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & PRESO,LLBELO,LLBLOW,ZRHMAX,ZRHMIN,&
        & YDGP_AD%PXP(:,:,:,KH),YDGP_AD%PXPD(:,:,:,KH),YDGP_AD%QF(:,:,KH),YDGP_AD%TF(:,:,KH), &
        & PXPP(:,:,IVAR),&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%QF(:,:,KH),YDGP5%TF(:,:,KH))

    ENDIF

  ELSEIF(IVARLIST(IVAR) == PPVAR%Q .OR. IVARLIST(IVAR) == PPVAR%O3) THEN

    ! Generic tracers all use "PPQ" routine
    IF (IVARLIST(IVAR) == PPVAR%Q) THEN
      ZFOR_INTERP5=>YDGP5%QF(:,:,KH)
      IF (LLTL) THEN
        ZFOR_INTERP=>YDGP_TL%QF(:,:,KH)
      ELSEIF (LLAD) THEN
        ZFOR_INTERP=>YDGP_AD%QF(:,:,KH)
      ENDIF
    ELSE !PPVAR%O3
      ZFOR_INTERP5=>YDGP5%O3F(:,:,KH)
      IF (LLTL) THEN
        ZFOR_INTERP=>YDGP_TL%O3F(:,:,KH)
      ELSEIF (LLAD) THEN
        ZFOR_INTERP=>YDGP_AD%O3F(:,:,KH)
      ENDIF
    ENDIF

    IF (LLDIRECT) THEN
      CALL PPQ(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & PRESO,LLBELO,LLBLOW,&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),ZFOR_INTERP5,&
        & PXPP(:,:,IVAR))

    ELSEIF (LLTL) THEN
      CALL PPQTL(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & PRESO,LLBELO,LLBLOW,&
        & YDGP_TL%PXP(:,:,:,KH),YDGP_TL%PXPD(:,:,:,KH),ZFOR_INTERP,&
        & PXPP(:,:,IVAR),&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),ZFOR_INTERP5)

    ELSEIF (LLAD) THEN
      CALL PPQAD(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & PRESO,LLBELO,LLBLOW,&
        & YDGP_AD%PXP(:,:,:,KH),YDGP_AD%PXPD(:,:,:,KH),ZFOR_INTERP,&
        & PXPP(:,:,IVAR),&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),ZFOR_INTERP5)

    ENDIF

  ELSEIF(IVARLIST(IVAR) == PPVAR%PWC .OR. IVARLIST(IVAR) == PPVAR%TO3) THEN

    ! Total column from p=preso to p=0.

    IF (IVARLIST(IVAR) == PPVAR%PWC) THEN
      ZFOR_INTERP5=>YDGP5%QF(:,:,KH)
      IF (LLTL) THEN
        ZFOR_INTERP=>YDGP_TL%QF(:,:,KH)
      ELSEIF (LLAD) THEN
        ZFOR_INTERP=>YDGP_AD%QF(:,:,KH)
      ENDIF
    ELSE !PPVAR%TO3
      ZFOR_INTERP5=>YDGP5%O3F(:,:,KH)
      IF (LLTL) THEN
        ZFOR_INTERP=>YDGP_TL%O3F(:,:,KH)
      ELSEIF (LLAD) THEN
        ZFOR_INTERP=>YDGP_AD%O3F(:,:,KH)
      ENDIF
    ENDIF

    IF (LLDIRECT) THEN
      CALL PPPWC(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & PRESO,ZLNPRESO,LLBELS,LLBLES,&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),ZFOR_INTERP5,&
        & PXPP(:,:,IVAR))

    ELSEIF (LLTL) THEN
!     ZFOR_INTERP(:,0)=0.0_JPRB
      CALL PPPWCTL(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & PRESO,ZLNPRESO,LLBELS,LLBLES,&
        & YDGP_TL%PXP(:,:,:,KH),YDGP_TL%PXPD(:,:,:,KH),ZFOR_INTERP,&
        & PXPP(:,:,IVAR),&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),ZFOR_INTERP5)

    ELSEIF (LLAD) THEN
      CALL PPPWCAD(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & PRESO,ZLNPRESO,LLBELS,LLBLES,&
        & YDGP_AD%PXP(:,:,:,KH),YDGP_AD%PXPD(:,:,:,KH),ZFOR_INTERP,&
        & PXPP(:,:,IVAR),&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),ZFOR_INTERP5)

    ENDIF

  ELSEIF(IVARLIST(IVAR) == PPVAR%Z) THEN

    IF (LLDIRECT) THEN
      CALL PPGEOP(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & PRESO,ZLNPRESO,LLBELO,LLBELS,LLBLOW,LLBLES,.TRUE.,YDGP5%OROG(:,KH),YDGP5%VGEOM,&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%TSTAR(:,KH),ZR05,&
        & YDGP5%PXYB(:,:,YYTXYB%M_LNPR,KH),YDGP5%PXYB(:,:,YYTXYB%M_ALPH,KH),YDGP5%TF(:,:,KH),&
        & PXPP(:,:,IVAR))

    ELSEIF (LLTL) THEN
      CALL PPGEOPTL(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & ZLNPRESO,LLBELO,LLBELS,LLBLOW,LLBLES,.TRUE.,YDGP5%OROG(:,KH),YDGP5%VGEOM,&
        & YDGP_TL%PXP(:,:,:,KH),YDGP_TL%PXPD(:,:,:,KH),YDGP_TL%TSTAR(:,KH),ZR0,&
        & YDGP_TL%PXYB(:,:,YYTXYB%M_LNPR,KH),YDGP_TL%PXYB(:,:,YYTXYB%M_ALPH,KH),YDGP_TL%TF(:,:,KH),&
        & PXPP(:,:,IVAR),&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%TSTAR(:,KH),ZR05,&
        & YDGP5%PXYB(:,:,YYTXYB%M_LNPR,KH),YDGP5%PXYB(:,:,YYTXYB%M_ALPH,KH),YDGP5%TF(:,:,KH))

    ELSEIF (LLAD) THEN
      CALL PPGEOPAD(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & ZLNPRESO,LLBELO,LLBELS,LLBLOW,LLBLES,.TRUE.,YDGP5%OROG(:,KH),YDGP5%VGEOM,&
        & YDGP_AD%PXP(:,:,:,KH),YDGP_AD%PXPD(:,:,:,KH),YDGP_AD%TSTAR(:,KH),ZR0,&
        & YDGP_AD%PXYB(:,:,YYTXYB%M_LNPR,KH),YDGP_AD%PXYB(:,:,YYTXYB%M_ALPH,KH),YDGP_AD%TF(:,:,KH),&
        & PXPP(:,:,IVAR),&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%TSTAR(:,KH),ZR05,&
        & YDGP5%PXYB(:,:,YYTXYB%M_LNPR,KH),YDGP5%PXYB(:,:,YYTXYB%M_ALPH,KH),YDGP5%TF(:,:,KH))

    ENDIF
 
  ELSEIF(IVARLIST(IVAR) == PPVAR%CC) THEN

    IF (LLDIRECT) THEN
      CALL PPCC(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & YDGP5%CCF(:,:,KH),&
        & PXPP(:,:,IVAR))

    ELSEIF (LLTL) THEN
      CALL PPCCTL(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & PRESO,LLBELO,LLBELS,LLBLOW,&
        & YDGP_TL%PXP(:,:,:,KH),YDGP_TL%PXPD(:,:,:,KH),YDGP_TL%CCF(:,:,KH),&
        & PXPP(:,:,IVAR),&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%CCF(:,:,KH))

    ELSEIF (LLAD) THEN
      CALL PPCCAD(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & YDGP_AD%CCF(:,:,KH),&
        & PXPP(:,:,IVAR),&
        & YDGP5%CCF(:,:,KH))

    ENDIF

  ELSEIF(IVARLIST(IVAR) == PPVAR%CLW .OR. IVARLIST(IVAR) == PPVAR%CIW) THEN

!CP just for checking
    !IF(IVARLIST(IVAR) == PPVAR%CLW) ZFOR_INTERP=>YDGP5%CLWF(:,:,KH)
    !IF(IVARLIST(IVAR) == PPVAR%CIW) ZFOR_INTERP=>YDGP5%CIWF(:,:,KH)
!CP end checking => to remove

    IF (IVARLIST(IVAR) == PPVAR%CLW) THEN
      ZFOR_INTERP5=>YDGP5%CLWF(:,:,KH)
      IF (LLTL) THEN
        ZFOR_INTERP=>YDGP_TL%CLWF(:,:,KH)
      ELSEIF (LLAD) THEN
        ZFOR_INTERP=>YDGP_AD%CLWF(:,:,KH)
      ENDIF
    ELSE !PPVAR%CIWF
      ZFOR_INTERP5=>YDGP5%CIWF(:,:,KH)
      IF (LLTL) THEN
        ZFOR_INTERP=>YDGP_TL%CIWF(:,:,KH)
      ELSEIF (LLAD) THEN
        ZFOR_INTERP=>YDGP_AD%CIWF(:,:,KH)
      ENDIF
    ENDIF

    IF (LLDIRECT) THEN
      CALL PPCLW(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & ZFOR_INTERP5,&
        & PXPP(:,:,IVAR))

    ELSEIF (LLTL) THEN
      CALL PPCLWTL(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & PRESO,LLBELO,LLBELS,LLBLOW,&
        & YDGP_TL%PXP(:,:,:,KH),YDGP_TL%PXPD(:,:,:,KH),ZFOR_INTERP,&
        & PXPP(:,:,IVAR),&
        & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),ZFOR_INTERP5)

    ELSEIF (LLAD) THEN
      CALL PPCLWAD(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
        & ZFOR_INTERP,&
        & PXPP(:,:,IVAR))

    ENDIF

  ELSEIF(IVARLIST(IVAR) == PPVAR%TCC) THEN

    IF (LLDIRECT) THEN
      CALL PPTCC(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,&
        & YDGP5%CCF(:,:,KH),&
        & PXPP(:,:,IVAR))

    ELSEIF (LLTL) THEN
      CALL PPTCCTL(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,&
        & YDGP_TL%CCF(:,:,KH),&
        & PXPP(:,:,IVAR),&
        & YDGP5%CCF(:,:,KH))

    ELSEIF (LLAD) THEN
      WRITE(NULERR,'('' *** ERROR IN PPOBSAAD: TOTAL CLOUD'',&
        & '' COVER OPERATOR NOT DONE YET '')')
      CALL ABOR1('ERROR PPOBSAAD: PPTCCAD DOES NOT EXIST')
!      CALL PPTCCAD(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,&
!        & YDGP_AD%CCF(:,1:YDGP5%NFLEVG,KH),&
!        & PXPP(:,:,IVAR),&
!        & YDGP5%CCF(:,:,KH))

    ENDIF

  ELSEIF(IVARLIST(IVAR)==PPVAR%PS) THEN

    IF (LLDIRECT) THEN
      CALL PPPS(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,KSLCT,&
        & LLEXTLI,PRESO,LLBELO,LLBELS,LLBLOW,LLBLES,.TRUE.,YDGP5%OROG(:,KH),&
        & YDGP5%GEOPH(:,:,KH),YDGP5%GEOPF(:,:,KH),YDGP5%PXZ(:,:,:,KH),&
        & YDGP5%PXZD(:,:,:,KH),YDGP5%PRESH(:,:,KH),YDGP5%PRESF(:,:,KH),&
        & YDGP5%TSTAR(:,KH),&
        & PXPP(:,:,IVAR))

    ELSEIF (LLTL) THEN
      CALL PPPSTL(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,KSLCT,&
        & LLEXTLI,PRESO,LLBELO,LLBELS,LLBLOW,LLBLES,.TRUE.,YDGP5%OROG(:,KH),&
        & YDGP_TL%GEOPH(:,:,KH),YDGP_TL%GEOPF(:,:,KH),YDGP_TL%PXZ(:,:,:,KH),&
        & YDGP_TL%PXZD(:,:,:,KH),YDGP_TL%PRESH(:,:,KH),YDGP_TL%PRESF(:,:,KH),&
        & YDGP_TL%TSTAR(:,KH),&
        & PXPP(:,:,IVAR),&
        & YDGP5%GEOPH(:,:,KH),YDGP5%GEOPF(:,:,KH),YDGP5%PXZ(:,:,:,KH),&
        & YDGP5%PXZD(:,:,:,KH),YDGP5%PRESH(:,:,KH),YDGP5%PRESF(:,:,KH),&
        & YDGP5%TSTAR(:,KH))

    ELSEIF (LLAD) THEN
      CALL PPPSAD(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,KSLCT,&
        & LLEXTLI,PRESO,LLBELO,LLBELS,LLBLOW,LLBLES,.TRUE.,YDGP5%OROG(:,KH),&
        & YDGP_AD%GEOPH(:,:,KH),YDGP_AD%GEOPF(:,:,KH),YDGP_AD%PXZ(:,:,:,KH),&
        & YDGP_AD%PXZD(:,:,:,KH),YDGP_AD%PRESH(:,:,KH),YDGP_AD%PRESF(:,:,KH),&
        & YDGP_AD%TSTAR(:,KH),&
        & PXPP(:,:,IVAR),&
        & YDGP5%GEOPH(:,:,KH),YDGP5%GEOPF(:,:,KH),YDGP5%PXZ(:,:,:,KH),&
        & YDGP5%PXZD(:,:,:,KH),YDGP5%PRESH(:,:,KH),YDGP5%PRESF(:,:,KH),&
        & YDGP5%TSTAR(:,KH))

    ENDIF

  ELSEIF(IVARLIST(IVAR)==PPVAR%T2) THEN

    IF (LLDIRECT) THEN
      CALL PPT2M(IDLEN,1,IDLEN,PRESO(:,1),&
        & YDGP5%QF(:,YDGP5%NFLEVG,KH),YDGP5%APHI(:,YDGP5%NFLEVG,KH),&
        & YDGP5%CPTGZ(:,KH),YDGP5%CPTS(:,KH),&
        & YDGP5%RI(:,KH),YDGP5%BN(:,KH),YDGP5%BH(:,KH),&
        & PXPP(:,:,IVAR))
      DO JC=2,IMXCNT
        WHERE(KCOUNT(:)>=JC) PXPP(:,JC,IVAR)=PXPP(:,1,IVAR)
      ENDDO

    ELSEIF (LLTL) THEN
      CALL PPT2MTL(IDLEN,1,IDLEN,PRESO(:,1),&
        & YDGP5%QF(:,YDGP5%NFLEVG,KH),YDGP5%APHI(:,YDGP5%NFLEVG,KH),&
        & YDGP_TL%QF(:,YDGP5%NFLEVG,KH),YDGP_TL%APHI(:,YDGP5%NFLEVG,KH),&
        & YDGP5%CPTGZ(:,KH),YDGP5%CPTS(:,KH),&
        & YDGP_TL%CPTGZ(:,KH),YDGP_TL%CPTS(:,KH),&
        & YDGP5%RI(:,KH),YDGP5%BN(:,KH),YDGP5%BH(:,KH),&
        & YDGP_TL%BN(:,KH),YDGP_TL%BH(:,KH),&
        & PXPP(:,:,IVAR))

    ELSEIF (LLAD) THEN
      CALL PPT2MAD(IDLEN,1,IDLEN,PRESO(:,1),&
        & YDGP5%QF(:,YDGP5%NFLEVG,KH),YDGP5%APHI(:,YDGP5%NFLEVG,KH),&
        & YDGP_AD%QF(:,YDGP5%NFLEVG,KH),YDGP_AD%APHI(:,YDGP5%NFLEVG,KH),&
        & YDGP5%CPTGZ(:,KH),YDGP5%CPTS(:,KH),&
        & YDGP_AD%CPTGZ(:,KH),YDGP_AD%CPTS(:,KH),&
        & YDGP5%RI(:,KH),YDGP5%BN(:,KH),YDGP5%BH(:,KH),&
        & YDGP_AD%BN(:,KH),YDGP_AD%BH(:,KH),&
        & PXPP(:,:,IVAR))

    ENDIF

  ELSEIF(IVARLIST(IVAR)==PPVAR%H2) THEN

!CP it is not clear how to manage this block: inside lldirect statement or outside as here
!CP   it seems if (direct) then PES5OUT present else tl or ad, it is ZES5
    IF(PRESENT(PES5OUT))THEN
      ZESOUT5=>PES5OUT
    ELSE
      ZESOUT5=>ZES5
    ENDIF

    ! H2 Always NONLINEAR so always call

    CALL PPRH2M(IDLEN,1,IDLEN,ZRHMAX,ZRHMIN,&
      & YDGP5%PRESF(:,YDGP5%NFLEVG,KH),YDGP5%QF(:,YDGP5%NFLEVG,KH),YDGP5%TF(:,YDGP5%NFLEVG,KH),&
      & ZESOUT5(:),&
      & ZXPP5(:,:,IVAR))
    DO JC=2,IMXCNT
      WHERE(KCOUNT(:)>=JC) ZXPP5(:,JC,IVAR)=ZXPP5(:,1,IVAR)
    ENDDO

    IF (LLTL) THEN
      CALL PPRH2MTL(IDLEN,1,IDLEN,ZRHMAX,ZRHMIN,&
        & YDGP5%PRESF(:,YDGP5%NFLEVG,KH),YDGP5%QF(:,YDGP5%NFLEVG,KH),YDGP5%TF(:,YDGP5%NFLEVG,KH),&
        & YDGP_TL%PRESF(:,YDGP5%NFLEVG,KH),YDGP_TL%QF(:,YDGP5%NFLEVG,KH),YDGP_TL%TF(:,YDGP5%NFLEVG,KH),&
        & ZXPP5,ZESOUT5(:),&
        & PXPP(:,:,IVAR))

    ELSEIF (LLAD) THEN
      CALL PPRH2MAD(IDLEN,1,IDLEN,ZRHMAX,ZRHMIN,&
        & YDGP5%PRESF(:,YDGP5%NFLEVG,KH),YDGP5%QF(:,YDGP5%NFLEVG,KH),YDGP5%TF(:,YDGP5%NFLEVG,KH),&
        & YDGP_AD%PRESF(:,YDGP5%NFLEVG,KH),YDGP_AD%QF(:,YDGP5%NFLEVG,KH),YDGP_AD%TF(:,YDGP5%NFLEVG,KH),&
        & ZXPP5,ZESOUT5(:),&
        & PXPP(:,:,IVAR))

    ENDIF

  ELSEIF(IVARLIST(IVAR)==PPVAR%COMPOSITION) THEN

    ! Composition variables: either extract a single variable specified by a single GRIB
    ! code, or all composition variables (GRIB code set to -1)
    ! Use LDGEMS_PWC to chose between total column output or vertical interpolation.
    LLGRIBFOUND=.FALSE.
    DO IGEMS=1,YDGP5%NGEMS
      IF((YDGP5%GEMS_IGRIB(IGEMS) == KGRIBCODE) .OR. (KGRIBCODE==-1 .AND. IGEMS==IVAR)) THEN
        LLGRIBFOUND=.TRUE.
        IF(LLGEMS_PWC) THEN
          IF (LLDIRECT) THEN
            CALL PPPWC(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
              & PRESO,ZLNPRESO,LLBELS,LLBLES,&
              & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%GEMS(:,:,KH,IGEMS),&
              & PXPP(:,:,IVAR))
          ELSEIF (LLTL) THEN
            CALL PPPWCTL(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
              & PRESO,ZLNPRESO,LLBELS,LLBLES,&
              & YDGP_TL%PXP(:,:,:,KH),YDGP_TL%PXPD(:,:,:,KH),YDGP_TL%GEMS(:,:,KH,IGEMS),&
              & PXPP(:,:,IVAR),&
              & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%GEMS(:,:,KH,IGEMS))     
          ELSEIF (LLAD) THEN
            CALL PPPWCAD(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
              & PRESO,ZLNPRESO,LLBELS,LLBLES,&
              & YDGP_AD%PXP(:,:,:,KH),YDGP_AD%PXPD(:,:,:,KH),YDGP_AD%GEMS(:,:,KH,IGEMS),&
              & PXPP(:,:,IVAR),&
              & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%GEMS(:,:,KH,IGEMS))
          ENDIF
        ELSE
          IF (LLDIRECT) THEN
            CALL PPQ(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
              & PRESO,LLBELO,LLBLOW,&
              & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%GEMS(:,:,KH,IGEMS),&
              & PXPP(:,:,IVAR))
          ELSEIF (LLTL) THEN
            CALL PPQTL(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
              & PRESO,LLBELO,LLBLOW,&
              & YDGP_TL%PXP(:,:,:,KH),YDGP_TL%PXPD(:,:,:,KH),YDGP_TL%GEMS(:,:,KH,IGEMS),&
              & PXPP(:,:,IVAR),&
              & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%GEMS(:,:,KH,IGEMS))     
          ELSEIF (LLAD) THEN
            CALL PPQAD(IDLEN,1,IDLEN,YDGP5%NFLEVG,IMXCNT,1,YDGP5%NPPM,ILEVB,&
              & PRESO,LLBELO,LLBLOW,&
              & YDGP_AD%PXP(:,:,:,KH),YDGP_AD%PXPD(:,:,:,KH),YDGP_AD%GEMS(:,:,KH,IGEMS),&
              & PXPP(:,:,IVAR),&
              & YDGP5%PXP(:,:,:,KH),YDGP5%PXPD(:,:,:,KH),YDGP5%GEMS(:,:,KH,IGEMS))
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    IF(.NOT.LLGRIBFOUND) CALL ABOR1('PPNEW: REQUESTED GRIB CODE NOT AVAILABLE')

  ENDIF

  IF (LLAD .AND. ANY(IVARLIST(IVAR)==(/PPVAR%T,PPVAR%Z/))) THEN
    CALL GPRCPAD(IDLEN,1,IDLEN,YDGP5%NFLEVG,PQ=YDGP_AD%QF(:,:,KH),PR=ZR0,LDTHERMACT=.TRUE.)
  ENDIF


ENDDO

! Deallocate
IF ((LLTL.OR.LLAD).AND.LL_NONLINEAR) THEN
  DEALLOCATE(ZXPP5)
ENDIF

IF (LHOOK) CALL DR_HOOK('PPNEW',1,ZHOOK_HANDLE)
END SUBROUTINE PPNEW
