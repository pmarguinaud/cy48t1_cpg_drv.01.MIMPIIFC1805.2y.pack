SUBROUTINE ESPNHQESI(&
 ! --- INPUT -----------------------------------------------------------------
 & YDCST, YDGEOMETRY,YDLDDH,YDRIP,YDDYN,KM,KMLOC,KSTA,KEND,LDONEM,LDDO_PARTSI,&
 ! --- INOUT -----------------------------------------------------------------
 & PSPVORG,PSPDIVG,PSPTG,PSPSPG,PSPSPDG,PSPSVDG,PSPBUFG,&
 & PSPTNDSI_VORG,PSPTNDSI_DIVG,PSPTNDSI_TG,PSPTNDSI_SPDG,PSPTNDSI_SVDG)

!**** *ESPNHQESI* - SPECTRAL SPACE COMPUTATIONS FOR LAM MODEL
!                   SEMI-IMPLICIT SCHEME (NHQE MODEL)

! !! The structure of ESPNHQESI must remain very close to SPNHQESI !!
!    Differences between ESPNHQESI and SPNHQESI are for example:
!    * use of YDGEOMETRY%YRELAP instead of YDGEOMETRY%YRLAP for Laplacian and some indexations.
!    * no use of option LSIDG.
!    * no use of option LIMPF.

!     Purpose.
!     --------
!       Semi-implicit scheme (inversion of linear system) in the NHQE model.
!       Calculation of a diagnostic value of Qcha.
!       Algorithm requires several iterations.

!**   Interface.
!     ----------
!        *CALL* *ESPNHQESI(...)

!        Explicit arguments :
!        --------------------
!         INPUT:
!          YDGEOMETRY      : Structure containing geometry
!          YDLDDH          : Structure containing DDH logical variables
!          YDRIP           : Structure containing timestep
!          YDDYN           : Structure containing dynamics
!          KM              : Zonal wavenumber "m"
!          KMLOC           : Zonal wavenumber "m" (DM-local numbering)
!          KSTA            : First column processed
!          KEND            : Last column processed
!          LDONEM          : T if only one "m" is processed
!          LDDO_PARTSI     : To select active parts according to iteration.

!         INOUT:
!          PSPVORG         : Vorticity columns
!          PSPDIVG         : Divergence columns
!          PSPTG           : Temperature columns
!          PSPSPG          : Surface pressure
!          PSPSPDG         : Pressure departure variable (Qcha) columns
!          PSPSVDG         : Vertical divergence variable (d4) columns
!          PSPBUFG         : Buffer to communicate data between the different iterations
!                            As input:
!                             PSPBUFG(.,.,1) contains DIV(iter=i-1)
!                             PSPBUFG(.,.,2) contains "Z16" of predictor (iter=0), if "iter>0"
!                             PSPBUFG(.,.,3) contains S_kap DIV(iter=i-1)
!                            As output:
!                             PSPBUFG(.,.,1) stores DIV(iter=i)
!                             PSPBUFG(.,.,2) stores "Z16" of predictor, if "iter=0"
!          PSPTNDSI_VORG   : Vorticity SI tendencies for DDH
!          PSPTNDSI_DIVG   : Divergence SI tendencies for DDH
!          PSPTNDSI_TG     : Temperature SI tendencies for DDH
!          PSPTNDSI_SPDG   : Pressure departure variable (Qcha) SI tendencies for DDH
!          PSPTNDSI_SVDG   : Vertical divergence variable (d4) SI tendencies for DDH

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation (IDSI)

!     Author.
!     -------
!        Karim YESSAD (METEO-FRANCE/CNRM/GMAP)
!        Date: march 2017

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM     ,JPRB
USE YOMHOOK      , ONLY : LHOOK,   DR_HOOK
USE YOMCST       , ONLY : TCST
USE YOMMP0       , ONLY : MYSETV
USE YOMDYN       , ONLY : TDYN
USE YOMLDDH      , ONLY : TLDDH
USE YOMRIP       , ONLY : TRIP
USE YOMCVER      , ONLY : LVERTFE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TLDDH)       ,INTENT(INOUT) :: YDLDDH
TYPE(TRIP)        ,INTENT(INOUT) :: YDRIP
TYPE(TDYN)        ,INTENT(INOUT) :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KM
INTEGER(KIND=JPIM),INTENT(IN)    :: KMLOC
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTA
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND
LOGICAL           ,INTENT(IN)    :: LDONEM
LOGICAL           ,INTENT(IN)    :: LDDO_PARTSI(3)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPVORG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPDIVG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPG(:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSPDG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPSVDG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPBUFG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND,3)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_VORG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_DIVG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_TG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_SPDG(:,:)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPTNDSI_SVDG(:,:)

!     ------------------------------------------------------------------

!! REAL(KIND=JPRB) :: ZZSPVORG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZZSPDIVG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZZSPTG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZZSPSPG(KSTA:KEND)
REAL(KIND=JPRB) :: ZZSPSPDG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZZSPSVDG(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)

REAL(KIND=JPRB) :: Z12(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: Z13(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: Z14(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: Z15(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: Z16(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: Z20(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: Z21(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: Z22(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: Z23(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: Z24(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: Z25(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: Z30(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: Z31(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: Z32(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: Z33(KSTA:KEND)
REAL(KIND=JPRB) :: Z34(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)

REAL(KIND=JPRB) :: ZSRHS(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZSRHS_VES(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZZSPDIVG_VES(YDGEOMETRY%YRDIMV%NFLEVG,KSTA:KEND)
REAL(KIND=JPRB) :: ZLAPDI(KSTA:KEND)

INTEGER(KIND=JPIM) :: IN, IM, ISP, IOFF, ISPCOL, JLEV, JSP

INTEGER(KIND=JPIM) :: IPART_2A,IPART_2BC,IPART_2D

REAL(KIND=JPRB) :: ZBDT, ZBDT2, ZH2, ZMBABAR
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "mxmaop.h"

#include "abor1.intfb.h"
#include "sigam.intfb.h"
#include "sitnu.intfb.h"
#include "siskap.intfb.h"
!!! #include "silkapi.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ESPNHQESI',0,ZHOOK_HANDLE)

ASSOCIATE( &
 & NPTRSV=>YDGEOMETRY%YRMP%NPTRSV, NPTRSVF=>YDGEOMETRY%YRMP%NPTRSVF, &
 & RSTRET=>YDGEOMETRY%YRGEM%RSTRET, &
 & NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & SIMO=>YDDYN%SIMO,SIMI=>YDDYN%SIMI,SIVP=>YDDYN%SIVP, &
 & SITR=>YDDYN%SITR, &
 & SI_ILAPKSSI=>YDDYN%SI_ILAPKSSI, &
 & RBTS2=>YDDYN%RBTS2, &
 & LIMPF=>YDDYN%LIMPF, &
 & LSIDG=>YDDYN%LSIDG, &
 & TDT=>YDRIP%TDT, &
 & LRSIDDH=>YDLDDH%LRSIDDH)
 
!     ------------------------------------------------------------------

!*       0.    PRELIMINARY CALCULATIONS.
!              -------------------------

IF (LIMPF) CALL ABOR1(' ESPNHQESI: LIMPF=T not available in LAM models')
IF (LSIDG) CALL ABOR1(' ESPNHQESI: LSIDG=T not available in LAM models')

IPART_2A=1
IPART_2BC=2
IPART_2D=3

IF (LDONEM) THEN
  IOFF=NPTRSVF(MYSETV)-1
ELSE
  IOFF=NPTRSV(MYSETV)-1
ENDIF
ISPCOL=KEND-KSTA+1

ZBDT=RBTS2*TDT
ZBDT2=(ZBDT*RSTRET)**2

ZH2=((YDCST%RD*SITR)/YDCST%RG)**2
ZMBABAR=RSTRET*RSTRET

!cdir noloopchg
DO JSP=KSTA,KEND
  IN=YDGEOMETRY%YRLAP%NVALUE(JSP+IOFF)
  IM=YDGEOMETRY%YRELAP%MVALUE(JSP+IOFF)
  ISP=YDGEOMETRY%YRELAP%NPME(IM)+IN
  ZLAPDI(JSP)=YDGEOMETRY%YRELAP%RLEPDIM(ISP)
ENDDO

!     ------------------------------------------------------------------

!*       1.    MEMORY TRANSFER.
!              ----------------

!! IF (LIMPF) ZZSPVORG(1:NFLEVG,KSTA:KEND)=PSPVORG(1:NFLEVG,KSTA:KEND)
ZZSPDIVG(1:NFLEVG,KSTA:KEND)=PSPDIVG(1:NFLEVG,KSTA:KEND)
ZZSPTG  (1:NFLEVG,KSTA:KEND)=PSPTG  (1:NFLEVG,KSTA:KEND)
ZZSPSPG (         KSTA:KEND)=PSPSPG (         KSTA:KEND)
ZZSPSVDG(1:NFLEVG,KSTA:KEND)=PSPSVDG(1:NFLEVG,KSTA:KEND)

IF (LRSIDDH) THEN
  ! DDH memory transfer.
  !! IF (LIMPF) PSPTNDSI_VORG(1:NFLEVG,KSTA:KEND)=-PSPVORG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_DIVG(1:NFLEVG,KSTA:KEND)=-PSPDIVG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_TG  (1:NFLEVG,KSTA:KEND)=-PSPTG  (1:NFLEVG,KSTA:KEND)
  PSPTNDSI_SVDG(1:NFLEVG,KSTA:KEND)=-PSPSVDG(1:NFLEVG,KSTA:KEND)
  !the case of surface pressure has not been treated yet
ENDIF

!     ------------------------------------------------------------------

!*       2.    SEMI-IMPLICIT SPECTRAL COMPUTATIONS.
!              ------------------------------------

!*        2.A  Preliminary initialisations, and non iterative part of RHS of Helmholtz equation.

IF (LDDO_PARTSI(IPART_2A)) THEN

  CALL SIGAM(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,Z12,ZZSPTG,ZZSPSPG,ISPCOL,NFLEVG)

  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      Z13(JLEV,JSP)=ZLAPDI(JSP)*Z12(JLEV,JSP)
      Z14(JLEV,JSP)=ZLAPDI(JSP)*ZZSPSVDG(JLEV,JSP)
    ENDDO
  ENDDO

  !!! CALL SILKAPI(YDGEOMETRY,YDDYN,KOPT=1,1,NFLEVG,ISPCOL,Z14,Z15)
  CALL MXMAOP(SI_ILAPKSSI(1,1,1),1,NFLEVG,Z14,1,NFLEVG,Z15,1,NFLEVG,NFLEVG,NFLEVG,ISPCOL)

  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      Z16(JLEV,JSP)=ZZSPDIVG(JLEV,JSP)-ZBDT*Z13(JLEV,JSP)-ZH2*Z15(JLEV,JSP)
      PSPBUFG(JLEV,JSP,2)=Z16(JLEV,JSP)
    ENDDO
  ENDDO

ENDIF

!*        2.B  Iterative part of RHS of Helmholtz equation, and complete RHS.

IF (LDDO_PARTSI(IPART_2BC)) THEN

  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      Z20(JLEV,JSP)=PSPBUFG(JLEV,JSP,1)
    ENDDO
  ENDDO

  CALL SISKAP(YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,ISPCOL,Z20,Z21)

  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      Z22(JLEV,JSP)=Z21(JLEV,JSP)-PSPBUFG(JLEV,JSP,3)
    ENDDO
  ENDDO

  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      Z23(JLEV,JSP)=ZMBABAR*Z22(JLEV,JSP)
    ENDDO
  ENDDO

  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      Z24(JLEV,JSP)=ZLAPDI(JSP)*Z23(JLEV,JSP)
    ENDDO
  ENDDO

  !!! CALL SILKAPI(YDGEOMETRY,YDDYN,KOPT=1,1,NFLEVG,ISPCOL,Z24,Z25)
  CALL MXMAOP(SI_ILAPKSSI(1,1,1),1,NFLEVG,Z24,1,NFLEVG,Z25,1,NFLEVG,NFLEVG,NFLEVG,ISPCOL)

  ! ZSRHS contains complete RHS of Helmholtz equation
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      ZSRHS(JLEV,JSP)=PSPBUFG(JLEV,JSP,2)+ZH2*Z25(JLEV,JSP)
    ENDDO
  ENDDO

ENDIF

!*        2.C  Invert Helmholtz equation through vertical eigenmodes space.

IF (LDDO_PARTSI(IPART_2BC)) THEN

  ! * go into vertical eigenmodes space.
  CALL MXMAOP(SIMI,1,NFLEVG,ZSRHS,1,NFLEVG,ZSRHS_VES,1,NFLEVG,NFLEVG,NFLEVG,ISPCOL)

  ! * solve Helmholtz equation.
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      ZZSPDIVG_VES(JLEV,JSP)=ZSRHS_VES(JLEV,JSP)/(1.0_JPRB-ZBDT2*SIVP(JLEV)*ZLAPDI(JSP))
    ENDDO
  ENDDO

  ! * return into current space.
  CALL MXMAOP(SIMO,1,NFLEVG,ZZSPDIVG_VES,1,NFLEVG,ZZSPDIVG,1,NFLEVG,NFLEVG,NFLEVG,ISPCOL)

  ! * save DIV(current iter)
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      PSPBUFG(JLEV,JSP,1)=ZZSPDIVG(JLEV,JSP)
    ENDDO
  ENDDO

ENDIF

!*        2.D  Solve the other equations.

IF (LDDO_PARTSI(IPART_2D)) THEN

  ! * provides final value of DIV.
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      ZZSPDIVG(JLEV,JSP)=PSPBUFG(JLEV,JSP,1)
    ENDDO
  ENDDO

  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      Z30(JLEV,JSP)=ZMBABAR*PSPBUFG(JLEV,JSP,3)
      Z31(JLEV,JSP)=Z30(JLEV,JSP)+ZZSPSVDG(JLEV,JSP)
    ENDDO
  ENDDO

  ! * provides final value of Qcha.
  !!! CALL SILKAPI(YDGEOMETRY,YDDYN,KOPT=1,1,NFLEVG,ISPCOL,Z31,ZZSPSPDG)
  CALL MXMAOP(SI_ILAPKSSI(1,1,1),1,NFLEVG,Z31,1,NFLEVG,ZZSPSPDG,1,NFLEVG,NFLEVG,NFLEVG,ISPCOL)

  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      ZZSPSPDG(JLEV,JSP)=(ZH2/(ZBDT*YDCST%RD*SITR))*ZZSPSPDG(JLEV,JSP)
    ENDDO
  ENDDO

  ! * provides final value of d4.
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      ZZSPSVDG(JLEV,JSP)=-Z30(JLEV,JSP)
    ENDDO
  ENDDO

  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      Z32(JLEV,JSP)=ZMBABAR*ZZSPDIVG(JLEV,JSP)
    ENDDO
  ENDDO

  CALL SITNU(LVERTFE, YDCST, YDGEOMETRY,YDDYN,1,NFLEVG,Z32,Z34,Z33,ISPCOL)

  ! * provides final value of log(prehyds).
  DO JSP=KSTA,KEND
    ZZSPSPG(JSP)=ZZSPSPG(JSP)-ZBDT*Z33(JSP)
  ENDDO

  ! * provides final value of T.
  DO JSP=KSTA,KEND
    DO JLEV=1,NFLEVG
      ZZSPTG(JLEV,JSP)=ZZSPTG(JLEV,JSP)-ZBDT*Z34(JLEV,JSP)
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

!*       3.    COMPUTATION OF SI TERM AT t+dt FOR DDH.
!              ---------------------------------------

!! IF (LIMPF) PSPVORG(1:NFLEVG,KSTA:KEND)=ZZSPVORG(1:NFLEVG,KSTA:KEND)
PSPDIVG(1:NFLEVG,KSTA:KEND)=ZZSPDIVG(1:NFLEVG,KSTA:KEND)
PSPTG  (1:NFLEVG,KSTA:KEND)=ZZSPTG  (1:NFLEVG,KSTA:KEND)
PSPSPG (         KSTA:KEND)=ZZSPSPG (         KSTA:KEND)
PSPSPDG(1:NFLEVG,KSTA:KEND)=ZZSPSPDG(1:NFLEVG,KSTA:KEND)
PSPSVDG(1:NFLEVG,KSTA:KEND)=ZZSPSVDG(1:NFLEVG,KSTA:KEND)

IF (LRSIDDH) THEN
  !! IF (LIMPF) PSPTNDSI_VORG(1:NFLEVG,KSTA:KEND)=PSPTNDSI_VORG(1:NFLEVG,KSTA:KEND)+PSPVORG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_DIVG(1:NFLEVG,KSTA:KEND)=PSPTNDSI_DIVG(1:NFLEVG,KSTA:KEND)+PSPDIVG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_TG(1:NFLEVG,KSTA:KEND)=PSPTNDSI_TG(1:NFLEVG,KSTA:KEND)+PSPTG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_SVDG(1:NFLEVG,KSTA:KEND)=PSPTNDSI_SVDG(1:NFLEVG,KSTA:KEND)+PSPSVDG(1:NFLEVG,KSTA:KEND)
  PSPTNDSI_SPDG(1:NFLEVG,KSTA:KEND)=0._JPRB
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ESPNHQESI',1,ZHOOK_HANDLE)
END SUBROUTINE ESPNHQESI
