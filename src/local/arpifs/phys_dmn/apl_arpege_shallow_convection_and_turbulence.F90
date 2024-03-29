SUBROUTINE APL_ARPEGE_SHALLOW_CONVECTION_AND_TURBULENCE (YDMF_PHYS_BASE_STATE, YDCPG_BNDS,           &
& YDCPG_OPTS, YDCPG_MISC, YDMF_PHYS, YDCPG_DYN0, YDMODEL, YDDDH, KNLAB_CVPP, PCDROV, PCHROV, PCOEFN, &
& PCONDCVPPI, PCONDCVPPL, PDIFCVPPQ, PDIFCVPPS, PEPS, PFLU_CD, PFLU_CH, PKQLROV, PKQROV, PKTROV,     &
& PKUROV, PMSC_LSCPE, PNBVNO, PNEBS, PNEBS0, PNEB_CVPP, PPFL_FPLCH, PPFL_FTKE, PPFL_FTKEI,           &
& PPRODTH_CVPP, PQI, PQIC, PQL, PQLC, PQLIS, PQLIS0, PQLI_CVPP, PQV, PTKE1, PXTROV, PXUROV)

USE PARKIND1, ONLY : JPIM, JPRB

USE MF_PHYS_BASE_STATE_TYPE_MOD &
                       , ONLY : MF_PHYS_BASE_STATE_TYPE
USE CPG_OPTS_TYPE_MOD   , ONLY : CPG_BNDS_TYPE, CPG_OPTS_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_MISC_TYPE, CPG_DYN_TYPE, &
                              & CPG_SL2_TYPE, CPG_GPAR_TYPE
USE MF_PHYS_TYPE_MOD   , ONLY : MF_PHYS_TYPE
USE CPG_TYPE_MOD       , ONLY : CPG_MISC_TYPE, CPG_DYN_TYPE, &
                              & CPG_SL2_TYPE, CPG_GPAR_TYPE
USE TYPE_MODEL         , ONLY : MODEL
USE DDH_MIX            , ONLY : TYP_DDH


IMPLICIT NONE

TYPE (MF_PHYS_BASE_STATE_TYPE), INTENT(IN)    :: YDMF_PHYS_BASE_STATE
TYPE(CPG_BNDS_TYPE),            INTENT(IN)    :: YDCPG_BNDS
TYPE(CPG_OPTS_TYPE),            INTENT(IN)    :: YDCPG_OPTS
TYPE(CPG_MISC_TYPE),            INTENT(INOUT) :: YDCPG_MISC
TYPE(MF_PHYS_TYPE),             INTENT(INOUT) :: YDMF_PHYS
TYPE(CPG_DYN_TYPE),             INTENT(IN)    :: YDCPG_DYN0
TYPE(MODEL),                    INTENT(IN)    :: YDMODEL
TYPE(TYP_DDH),                  INTENT(INOUT) :: YDDDH
INTEGER(KIND=JPIM),             INTENT(OUT)   :: KNLAB_CVPP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(IN)    :: PCDROV(YDCPG_OPTS%KLON)
REAL(KIND=JPRB),                INTENT(IN)    :: PCHROV(YDCPG_OPTS%KLON)
REAL(KIND=JPRB),                INTENT(OUT)   :: PCOEFN(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(INOUT) :: PCONDCVPPI(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(INOUT) :: PCONDCVPPL(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(INOUT) :: PDIFCVPPQ(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(INOUT) :: PDIFCVPPS(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(IN)    :: PEPS
REAL(KIND=JPRB),                INTENT(IN)    :: PFLU_CD (YDCPG_OPTS%KLON)
REAL(KIND=JPRB),                INTENT(IN)    :: PFLU_CH (YDCPG_OPTS%KLON)
REAL(KIND=JPRB),                INTENT(OUT)   :: PKQLROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(OUT)   :: PKQROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(OUT)   :: PKTROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(OUT)   :: PKUROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(IN)    :: PMSC_LSCPE (YDCPG_OPTS%KLON, 1:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(OUT)   :: PNBVNO(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(OUT)   :: PNEBS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(OUT)   :: PNEBS0(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(INOUT) :: PNEB_CVPP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(IN)    :: PPFL_FPLCH (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(OUT)   :: PPFL_FTKE (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(OUT)   :: PPFL_FTKEI (YDCPG_OPTS%KLON, 0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(INOUT) :: PPRODTH_CVPP(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(IN)    :: PQI(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(IN)    :: PQIC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(IN)    :: PQL(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(IN)    :: PQLC(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(OUT)   :: PQLIS(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(OUT)   :: PQLIS0(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(INOUT) :: PQLI_CVPP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(IN)    :: PQV(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(OUT)   :: PTKE1(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(OUT)   :: PXTROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB),                INTENT(OUT)   :: PXUROV(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)

#include "actke.intfb.h"
#include "acvppkf.intfb.h"

REAL(KIND=JPRB) :: ZNLAB(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZNLABCVP(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZTPRDY(YDCPG_OPTS%KLON,YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZCDROV_SAVE(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZCHROV_SAVE(YDCPG_OPTS%KLON)
REAL(KIND=JPRB) :: ZKTROV_SAVE(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
REAL(KIND=JPRB) :: ZKUROV_SAVE(YDCPG_OPTS%KLON,0:YDCPG_OPTS%KFLEVG)
INTEGER(KIND=JPIM) :: INND(YDCPG_OPTS%KLON)

INTEGER (KIND=JPIM) :: JLEV
INTEGER (KIND=JPIM) :: JLON


! SHALLOW CONVECTION AND TURBULENCE

! shallow convection
CALL ACVPPKF( YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA, YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDMODEL%YRML_PHY_MF%YRTOPH%NTCVIM,              &
& YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, &
& YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R, YDMF_PHYS_BASE_STATE%T, PQV, PQL, PQI, YDMF_PHYS_BASE_STATE%U,                                   &
& YDMF_PHYS_BASE_STATE%V, YDCPG_DYN0%CTY%VVEL(:, 1:), YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%CP, YDMF_PHYS_BASE_STATE%TKE,                   &
& PDIFCVPPQ, PDIFCVPPS, PCONDCVPPL, PCONDCVPPI, PPRODTH_CVPP, KNLAB_CVPP, PQLI_CVPP, PNEB_CVPP, INND                                    &
&                                                                           )
  

YDCPG_MISC%QICE(:,:)= PQI(:,:)
YDCPG_MISC%QLI(:,:) = PQL(:,:)

! Computation of the 2 znlab used in acbl89 (link between deep/shallow convection and turbulence)
    
ZNLABCVP(:,:) = 1.0_JPRB
    
IF(YDMODEL%YRML_PHY_MF%YRPHY%LNEBN.OR.YDMODEL%YRML_PHY_MF%YRPHY%LRRGUST) THEN
  DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
    DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
      ZNLABCVP(JLON,JLEV) = ZNLABCVP(JLON,JLEV)&
      & *MAX(0.0_JPRB,SIGN(1.0_JPRB,PPFL_FPLCH(JLON,JLEV)-PPFL_FPLCH(JLON,JLEV-1)-PEPS))
      ZNLAB(JLON,JLEV) = REAL(KNLAB_CVPP(JLON,JLEV),JPRB)
    ENDDO
  ENDDO
ENDIF
DO JLEV=YDCPG_OPTS%KTDIA,YDCPG_OPTS%KFLEVG
  DO JLON=YDCPG_BNDS%KIDIA,YDCPG_BNDS%KFDIA
    ZNLAB(JLON,JLEV) = REAL(KNLAB_CVPP(JLON,JLEV),JPRB)
  ENDDO
ENDDO

! turbulence
CALL ACTKE ( YDMODEL%YRML_DIAG%YRLDDH, YDMODEL%YRML_DIAG%YRMDDH, YDMODEL%YRML_PHY_MF, YDCPG_BNDS%KIDIA,                           &
& YDCPG_BNDS%KFDIA, YDCPG_OPTS%KLON, YDMODEL%YRML_PHY_MF%YRTOPH%NTCOEF, YDMODEL%YRML_PHY_MF%YRTOPH%NTCOET,                        &
& YDCPG_OPTS%KFLEVG, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHI, YDMF_PHYS_BASE_STATE%YCPG_DYN%PHIF, YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYD, &
& YDMF_PHYS_BASE_STATE%YCPG_PHY%PREHYDF, YDMF_PHYS_BASE_STATE%YCPG_PHY%XYB%DELP, YDMF_PHYS_BASE_STATE%YCPG_DYN%RCP%R,             &
& YDMF_PHYS_BASE_STATE%T, YDMF_PHYS_BASE_STATE%U, YDMF_PHYS_BASE_STATE%V, PQV, PQIC, PQLC,                                        &
& PMSC_LSCPE, PFLU_CD, PFLU_CH, YDMF_PHYS%OUT%GZ0, YDMF_PHYS_BASE_STATE%YGSP_RR%T, YDCPG_MISC%QS, YDCPG_MISC%QICE,                &
& YDCPG_MISC%QLI, YDMF_PHYS_BASE_STATE%TKE, PPRODTH_CVPP, ZNLAB, ZNLABCVP, PKTROV, PKQROV, PKQLROV,                               &
& PKUROV, PXTROV, PXUROV, PNBVNO, PNEBS, PQLIS, PNEBS0, PQLIS0, PCOEFN, PPFL_FTKE, PPFL_FTKEI, PTKE1,                             &
& ZTPRDY, YDMF_PHYS%OUT%EDR, YDDDH)

YDMF_PHYS%OUT%CLPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)=MIN(YDMODEL%YRML_PHY_MF%YRPHY0%XMAXLM,MAX(YDMODEL%YRML_PHY_MF%YRPHY0%XMINLM,YDMF_PHYS%OUT%CLPH(YDCPG_BNDS%KIDIA:YDCPG_BNDS%KFDIA)))
  

     !-------------------------------------------------
     ! Store diffusion coefficients in trajectory in temporary variables
     ! before final writing.
     !-------------------------------------------------
ZKTROV_SAVE(:,:)=PKTROV(:,:)
ZKUROV_SAVE(:,:)=PKUROV(:,:)
ZCDROV_SAVE(:)=PCDROV(:)
ZCHROV_SAVE(:)=PCHROV(:)



END SUBROUTINE
