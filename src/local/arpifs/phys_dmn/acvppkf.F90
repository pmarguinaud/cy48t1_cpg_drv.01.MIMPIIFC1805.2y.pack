!OPTIONS XOPT(NOEVAL)
!-----------------------------------------------------------------
SUBROUTINE ACVPPKF( YDCST, YDML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV, &
 !-----------------------------------------------------------------
 ! - INPUT  2D .
 & PAPRSF, PAPHIF, PDELP, PR, PT, PQ, &
 & PQL, PQI, PU, PV, PVERVEL, PCP, PTKE, &
 ! - OUTPUT 2D .
 & PDIFCQ, PDIFCS, PFCCQL, PFCCQN, PPRODTH, &
 & KNLAB, PQCPP, PNEBPP,&
 ! - OUTPUT 1D .
 & KNND)

!-----------------------------------------------------------------

!  Authors  : E. Bazile and P. Bechtold  (CNRM/GMAP et L.A.)

!-----------------------------------------------------------------

!  Modified : 
!  05/2002    phased with CONVECTION call for IFS/ECMWF 
!            (routine cucalln.F90 calling both Tiedtke convection scheme
!             and present scheme)
!             ouput of present scheme (updraft QL and QV) provides also
!             necessary parameters for Tiedtke prognostic cloud scheme
!  03/2002  P. Marquet.  new  ZFHMLTS, ZFHEVPP in CPFHPRS (for Lopez)
!  03/2002  P. Marquet.  new  LKFDEEP, LKFSHAL
!  03/2002  P. Marquet.  "call deep_convection" 
!                      > "call convection" (deep + shallow)
!  09/2006  E. Bazile : Appel de la routine de shallow convection d'AROME
!                       uniquement
!  04/2008  E. Bazile : calcul du terme de production thermique PPROTH
!  10/2008  Y. Bouteloup & F. Bouyssel : Correction of bugs in initialization
!  07/2009  E. Bazile : TKE en entree de KFB et W fct de W_conv
!  K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!  10/2009  F. Bouyssel : Limitation on maximal TKE value
!  02/2010  E. Bazile : Correction for W without TKE scheme.
!  04/2010  F. Bouyssel : Bug correction on KNLAB computation
!  09/2010  O. Spaniel : Bug correction in expression SQRT(MIN)
!  04/2011  F. Bouyssel : Correction of a jlon loop (kidia,kfdia)
!  12/2012  E. Bazile   : Modif of W_turb and qc and cc fct of mass flux.

!  Peter.Bechtold@ecmwf.int

! Sequence de  routines :
! aplpar > acvppkf > convection_shal

! iv)  Momentum transport:
!      Option LLUVTRANS: c'est possible d'utiliser maintenant
!           mais pas encore bien teste. Donc par defaut mettre
!           LLUVTRANS=.FALSE.

! vi)   Traceurs passifs - chimie: 
!      Cette partie est utilisee uniquement dans MOCAGE et dans MESONH
!      Si on ne veut pas de transport de traceurs (ex. Ozone,CO) dans 
!      ARPEGE/ECMWF IFS, mettre tout siplement OCHTRANS=FALSE et KCH1=0 
!      (nombre de traceurs). PCH1 (traceur) et PCH1TEN (tendance 
!      convective du traceur) ont alors les dimensions
!      (KLON,KLEV,KCH1=0) qui ne prennent pas de place.

!-----------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCST  , ONLY :  TCST
USE YOMLSFORC, ONLY : LMUSCLFA,NMUSCLFA
!-----------------------------------------------------------------

IMPLICIT NONE

TYPE (TCST), INTENT (IN) :: YDCST
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KLON
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KTDIA 
INTEGER(KIND=JPIM) ,INTENT(IN)    :: KLEV
REAL(KIND=JPRB)    ,INTENT(IN)    :: PAPRSF (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PAPHIF (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PDELP  (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PR     (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PT     (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PQ     (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PQL    (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PQI    (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PU     (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PV     (KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PVERVEL(KLON,KLEV) 
REAL(KIND=JPRB)    ,INTENT(IN)    :: PCP    (KLON,KLEV)
REAL(KIND=JPRB)    ,INTENT(IN)    :: PTKE   (KLON,KLEV)
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PDIFCQ (KLON,0:KLEV) 
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PDIFCS (KLON,0:KLEV) 
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PFCCQL (KLON,0:KLEV) 
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PFCCQN (KLON,0:KLEV) 
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PPRODTH(KLON,0:KLEV)
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PQCPP  (KLON,KLEV)
REAL(KIND=JPRB)    ,INTENT(INOUT) :: PNEBPP (KLON,KLEV)
INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KNLAB  (KLON,KLEV) 
INTEGER(KIND=JPIM) ,INTENT(OUT)   :: KNND   (KLON) 

!-----------------------------------------------------------------

LOGICAL :: LLREFRESH_ALL, LLDOWN, LLUVTRANS, LLOCHTRANS, LLCONDWT

REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL(KIND=JPRB) :: ZDTCONV, ZVMD, ZWMD, ZSMD, ZTDCP, ZEPS, ZDQCDT, ZDTLDT

INTEGER(KIND=JPIM) :: JLON, JLEV, I_KBDIA, I_KCH1, IKICE

INTEGER(KIND=JPIM) :: I_KCOUNT(KLON)
INTEGER(KIND=JPIM) :: I_KCLTOP(KLON)
INTEGER(KIND=JPIM) :: I_KCLBAS(KLON)

REAL(KIND=JPRB) :: ZCAPE(KLON)

REAL(KIND=JPRB) :: ZW     (KLON,KLEV)
REAL(KIND=JPRB) :: ZDTDT  (KLON,KLEV)
REAL(KIND=JPRB) :: ZDQVDT (KLON,KLEV)
REAL(KIND=JPRB) :: ZDQLDT (KLON,KLEV)
REAL(KIND=JPRB) :: ZDQIDT (KLON,KLEV)
REAL(KIND=JPRB) :: ZDUDT  (KLON,KLEV)
REAL(KIND=JPRB) :: ZDVDT  (KLON,KLEV)
REAL(KIND=JPRB) :: ZUMF   (KLON,KLEV)
REAL(KIND=JPRB) :: ZUQV   (KLON,KLEV)
REAL(KIND=JPRB) :: ZUQL   (KLON,KLEV)
REAL(KIND=JPRB) :: ZDPSG  (KLON,KLEV)
REAL(KIND=JPRB) :: ZLV    (KLON,KLEV)
REAL(KIND=JPRB) :: ZLS    (KLON,KLEV)
REAL(KIND=JPRB) :: ZQC    (KLON,KLEV)
REAL(KIND=JPRB) :: ZBETA  (KLON,KLEV)
REAL(KIND=JPRB) :: ZAPHIF (KLON,KLEV)

REAL(KIND=JPRB) :: ZTHETA (KLON,0:KLEV+1)
REAL(KIND=JPRB) :: ZRHO   (KLON,0:KLEV+1)

REAL(KIND=JPRB) :: ZCH1   (KLON,KLEV,0)
REAL(KIND=JPRB) :: ZCH1TEN(KLON,KLEV,0)
REAL(KIND=JPRB) :: ZTKECLS(KLON)
REAL(KIND=JPRB) :: ZUMFMAX(KLON)

!-----------------------------------------------------------------

#include "convection_shal.h"
#include "fcttrm.ycst.h"
#include "wrscmr.intfb.h"

!-----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACVPPKF',0,ZHOOK_HANDLE)
ASSOCIATE(RKFBNBX=>YDML_PHY_MF%YRPHY0%RKFBNBX, RKFBTAU=>YDML_PHY_MF%YRPHY0%RKFBTAU, RQLCR=>YDML_PHY_MF%YRPHY0%RQLCR, &
 & RPRTH=>YDML_PHY_MF%YRPHY0%RPRTH, ECTMIN=>YDML_PHY_MF%YRPHY0%ECTMIN, AECLS4=>YDML_PHY_MF%YRPHY0%AECLS4, &
 & TSPHY=>YDML_PHY_MF%YRPHY2%TSPHY, &
 & XNHGAM=>YDML_PHY_MF%YRCVMNH%XNHGAM, XA25=>YDML_PHY_MF%YRCVMNH%XA25, XTFRZ1=>YDML_PHY_MF%YRCVMNH%XTFRZ1, &
 & XTFRZ2=>YDML_PHY_MF%YRCVMNH%XTFRZ2, XENTR=>YDML_PHY_MF%YRCVMNH%XENTR, LSMOOTH=>YDML_PHY_MF%YRCVMNH%LSMOOTH, &
 & XZLCL=>YDML_PHY_MF%YRCVMNH%XZLCL, OTADJS=>YDML_PHY_MF%YRCVMNH%OTADJS, XBW=>YDML_PHY_MF%YRCVMNH%XBW, &
 & XDTPERT=>YDML_PHY_MF%YRCVMNH%XDTPERT, XCRAD=>YDML_PHY_MF%YRCVMNH%XCRAD, XBTPERT=>YDML_PHY_MF%YRCVMNH%XBTPERT, &
 & XAW=>YDML_PHY_MF%YRCVMNH%XAW, XCDEPTH=>YDML_PHY_MF%YRCVMNH%XCDEPTH, XZPBL=>YDML_PHY_MF%YRCVMNH%XZPBL, &
 & XATPERT=>YDML_PHY_MF%YRCVMNH%XATPERT, XSTABT=>YDML_PHY_MF%YRCVMNH%XSTABT, LSETTADJ=>YDML_PHY_MF%YRCVMNH%LSETTADJ, &
 & XCDEPTH_D=>YDML_PHY_MF%YRCVMNH%XCDEPTH_D, XWTRIG=>YDML_PHY_MF%YRCVMNH%XWTRIG, XSTABC=>YDML_PHY_MF%YRCVMNH%XSTABC, &
 & LECT=>YDML_PHY_MF%YRPHY%LECT, LCVDD=>YDML_PHY_MF%YRPHY%LCVDD)
!-----------------------------------------------------------------

ZVMD=YDCST%RCPV-YDCST%RCPD
ZWMD=YDCST%RCW-YDCST%RCPD
ZSMD=YDCST%RCS-YDCST%RCPD

ZUMFMAX(:)=1.E-12_JPRB
ZW(:,:)=0.0_JPRB
ZTKECLS(:)=0.0_JPRB
IF (LECT) THEN
  DO JLEV=1,KLEV
    DO JLON=KIDIA,KFDIA
      ZW(JLON,JLEV) = SQRT(MIN(3.0_JPRB,MAX(ECTMIN,PTKE(JLON,JLEV)))/AECLS4)
    ENDDO
  ENDDO
  DO JLON=KIDIA,KFDIA
     ZTKECLS(JLON)=PTKE(JLON,KLEV) 
  ENDDO
ENDIF
DO JLEV=1,KLEV
  DO JLON=KIDIA,KFDIA
    ZDPSG (JLON,JLEV) = PDELP(JLON,JLEV)/YDCST%RG
    ZAPHIF(JLON,JLEV) = PAPHIF(JLON,JLEV)/YDCST%RG
  ENDDO
ENDDO

DO JLEV=1,KLEV
  DO JLON=KIDIA,KFDIA
    ZDTDT (JLON,JLEV) = 0.0_JPRB
    ZDQVDT(JLON,JLEV) = 0.0_JPRB
    ZDQLDT(JLON,JLEV) = 0.0_JPRB
    ZDQIDT(JLON,JLEV) = 0.0_JPRB
    ZDUDT (JLON,JLEV) = 0.0_JPRB
    ZDVDT (JLON,JLEV) = 0.0_JPRB
    ZUMF  (JLON,JLEV) = 0.0_JPRB
    ZUQV  (JLON,JLEV) = 0.0_JPRB
    ZUQL  (JLON,JLEV) = 0.0_JPRB
    ZLV   (JLON,JLEV) = FOLH(PT(JLON,JLEV),0.0_JPRB)
    ZLS   (JLON,JLEV) = FOLH(PT(JLON,JLEV),1.0_JPRB)
    ZQC   (JLON,JLEV) = PQL(JLON,JLEV)+PQI(JLON,JLEV)
  ENDDO
ENDDO

DO JLON=1,KLON
  I_KCOUNT(JLON) = 0
  ZCAPE   (JLON) = 0.0_JPRB
ENDDO

I_KBDIA=1
IKICE=1
I_KCH1=0
LLREFRESH_ALL=.TRUE.
LLDOWN=LCVDD
LLUVTRANS=.FALSE. ! not yet well tested but possible to use
LLOCHTRANS=.FALSE.
ZDTCONV=TSPHY

! - - - - - - - - - - - - - -
! Arguments de : CONVECTION
! - - - - - - - - - - - - - -

! KLON      ! horizontal dimension
! KLEV      ! vertical dimension
! KIDIA     ! value of the first point in x
! KFDIA     ! value of the last point in x
! KBDIA     ! vertical  computations start at KBDIA (that is at least 1)
! KTDIA     ! vertical computations can belimited to KLEV+1-KTDIA (default=1)

! ZDTCONV   ! Interval of time between two calls of the deep convection scheme
! LLREFRESH_ALL ! refresh or not all tendencies  at every call
! LLODOWN   ! take or not convectivedowndrafts into account
! IKICE     ! flag for ice ( 1 = yes , 0 = no ice )

! PQ     ! grid scale water vapor (kg/kg)
! PQL    ! grid scale r_c  (kg/kg)
! PQI    ! grid scale r_i  (kg/kg)
! PU     ! grid scale wind in x (m/s)
! PV     ! grid scale wind in y (m/s)
! ZW     ! grid scale vertical velocity (m/s)

! KCOUNT ! convective counter (recompute tendency or keep it)
! ZDTDT  ! convective temperat. tendency (K/s)
! ZDQVDT ! convective r_v tendency (1/s)
! ZDQLDT ! convective r_c tendency (1/s)

! Diagnostic variables:

! ZUMF   ! updraft mass flux   (kg/s m2)
! ZCAPE  ! CAPE (J/kg)
! KCLTOP ! cloud top level  (number of model level)
! KCLBAS ! cloud base level (number of model level)
!        ! they are given a value of 0 if no convection
! Attention default value of KCLTOP,KCLBAS = 1 if no convection

! Momentum transport:

! LLUVTRANS ! compute convective tendencies for horizontal wind
! ZDUDT  ! convective tendency for u (m/s^2)
! ZDVDT  ! convective tendency for v (m/s^2)

! Chemical Tracers:

! LLOCHTRANS ! flag to compute convective transport for chemical tracer
! KCH1       ! number of species
! ZCH1       ! grid scale chemical species
! ZCH1TEN    ! chemical convective tendency (1/s)

CALL  CONVECTION_SHAL( KLON, KLEV, KIDIA, KFDIA, I_KBDIA, KTDIA,&
 & ZDTCONV, LLREFRESH_ALL, LLDOWN, IKICE,&
 & LSETTADJ, OTADJS, LSMOOTH,&
 & XA25, XCRAD, XCDEPTH, XCDEPTH_D, XDTPERT, XATPERT, XBTPERT,XENTR,&
 & XZLCL, XZPBL, XWTRIG, XNHGAM, XTFRZ1, XTFRZ2, XSTABT, XSTABC, XAW, XBW,&
 & PAPRSF, ZAPHIF, ZTKECLS,&
 & PT, PQ, PQL, PQI, PU, PV, ZW,&
 & I_KCOUNT, ZDTDT, ZDQVDT, ZDQLDT, ZDQIDT,&
 & ZUMF, ZCAPE, I_KCLTOP, I_KCLBAS,&
 & ZUQV, ZUQL,&
 & LLUVTRANS, ZDUDT, ZDVDT,&
 & LLOCHTRANS, I_KCH1, ZCH1, ZCH1TEN )

IF(LMUSCLFA) CALL WRSCMR(NMUSCLFA,'ZMF_shal',ZUMF,KLON,KLEV)

! Calcul de la production thermique pour la TKE
PPRODTH(:,:)=0.0_JPRB
IF (RPRTH > 0._JPRB) THEN
  DO JLEV=1,KLEV
    DO JLON=KIDIA,KFDIA
      ZBETA (JLON,JLEV) = (YDCST%RATM/PAPRSF(JLON,JLEV))**(YDCST%RKAPPA)
      ZTHETA(JLON,JLEV) = PT(JLON,JLEV)*ZBETA(JLON,JLEV)
      ZRHO  (JLON,JLEV) = PAPRSF(JLON,JLEV)/(PR(JLON,JLEV)*PT(JLON,JLEV))
    ENDDO
  ENDDO
  DO JLON=KIDIA,KFDIA
    ZTHETA(JLON,0)      = ZTHETA(JLON,1)
    ZRHO  (JLON,0)      = ZRHO  (JLON,1)
    ZTHETA(JLON,KLEV+1) = ZTHETA(JLON,KLEV)
    ZRHO  (JLON,KLEV+1) = ZRHO  (JLON,KLEV)
  ENDDO

  DO JLEV=1,KLEV
    DO JLON=KIDIA,KFDIA
      ZTDCP=ZVMD*ZDQVDT(JLON,JLEV)+ZWMD*ZDQLDT(JLON,JLEV)+ZSMD*ZDQIDT(JLON,JLEV)
      ZDQCDT=ZDQLDT(JLON,JLEV)+ZDQIDT(JLON,JLEV)
      ZDTLDT=ZBETA (JLON,JLEV)&
        &  * ( ZDTDT(JLON,JLEV) + ZLV(JLON,JLEV)/PCP(JLON,JLEV)&
        &    * ( ZQC(JLON,JLEV)*ZTDCP/PCP(JLON,JLEV)-ZDQCDT ) )
      PPRODTH(JLON,JLEV)=PPRODTH(JLON,JLEV-1)-ZDPSG(JLON,JLEV)*ZDTLDT
    ENDDO
  ENDDO

  DO JLEV=0,KLEV
    DO JLON=KIDIA,KFDIA
      PPRODTH(JLON,JLEV)=PPRODTH(JLON,JLEV)*YDCST%RG*4._JPRB&
        & / ( ZRHO  (JLON,JLEV) + ZRHO  (JLON,JLEV+1) )&
        & / ( ZTHETA(JLON,JLEV) + ZTHETA(JLON,JLEV+1) )
    ENDDO
  ENDDO
ENDIF ! Fin du calcul de la production thermique pour la TKE

DO JLON=KIDIA,KFDIA
  DO JLEV=1,KLEV
    KNLAB(JLON,JLEV)=1-MAX(0,MIN(1,(I_KCLTOP(JLON)-JLEV)*(I_KCLBAS(JLON)-JLEV)))
  ENDDO
  KNND(JLON)=MIN(1,I_KCLTOP(JLON)-1)
ENDDO

! Calcul de la nebulosite et de l'eau condensee
ZEPS=1.E-12_JPRB
IF (RKFBTAU > 0._JPRB) THEN
  DO JLEV=1,KLEV
!DEC$ IVDEP
    DO JLON=KIDIA,KFDIA
      ZDQCDT=MAX(0.0_JPRB,ZDQLDT(JLON,JLEV)+ZDQIDT(JLON,JLEV))
      PQCPP (JLON,JLEV)=RKFBTAU*ZDQCDT*FLOAT(KNLAB(JLON,JLEV))
      PNEBPP(JLON,JLEV)=MAX(ZEPS,MIN(RKFBNBX,PQCPP(JLON,JLEV)/RQLCR))
      ZUMFMAX(JLON)=MAX(ZUMFMAX(JLON),ZUMF(JLON,JLEV))
    ENDDO
  ENDDO
  IF (.NOT.LSMOOTH) THEN
    DO JLEV=1,KLEV
!DEC$ IVDEP
      DO JLON=KIDIA,KFDIA
        PQCPP (JLON,JLEV)=PQCPP(JLON,JLEV)&
         & *MIN(1.0_JPRB,ZUMF(JLON,JLEV)/ZUMFMAX(JLON))
        PNEBPP(JLON,JLEV)=MAX(ZEPS,MIN(RKFBNBX,PQCPP(JLON,JLEV)/RQLCR))
      ENDDO
    ENDDO
  ENDIF
ENDIF ! Fin du calcul de la nebulosite et de l'eau condensee

LLCONDWT=.FALSE.
IF (.NOT.LLCONDWT) THEN
  DO JLEV=1,KLEV
    DO JLON=KIDIA,KFDIA
      ZDQVDT(JLON,JLEV)=ZDQVDT(JLON,JLEV)+ZDQLDT(JLON,JLEV)+ZDQIDT(JLON,JLEV)
      ZDTDT(JLON,JLEV)=ZDTDT(JLON,JLEV) - (ZLV(JLON,JLEV)*ZDQLDT(JLON,JLEV)&
        & + ZLS(JLON,JLEV)*ZDQIDT(JLON,JLEV)) / PCP(JLON,JLEV)
      ZDQLDT(JLON,JLEV)=0.0_JPRB
      ZDQIDT(JLON,JLEV)=0.0_JPRB
    ENDDO
  ENDDO
ELSE
  DO JLEV=1,KLEV
    DO JLON=KIDIA,KFDIA
      PFCCQL(JLON,JLEV)=PFCCQL(JLON,JLEV-1)+ZDPSG(JLON,JLEV)*ZDQLDT(JLON,JLEV)
      PFCCQN(JLON,JLEV)=PFCCQN(JLON,JLEV-1)+ZDPSG(JLON,JLEV)*ZDQIDT(JLON,JLEV)
    ENDDO
  ENDDO
ENDIF

DO JLEV=1,KLEV
!DEC$ IVDEP
  DO JLON=KIDIA,KFDIA
    PDIFCQ(JLON,JLEV)=PDIFCQ(JLON,JLEV-1)-ZDPSG(JLON,JLEV)*ZDQVDT(JLON,JLEV)
    ZTDCP=ZVMD*ZDQVDT(JLON,JLEV)+ZWMD*ZDQLDT(JLON,JLEV)+ZSMD*ZDQIDT(JLON,JLEV)
    PDIFCS(JLON,JLEV)=PDIFCS(JLON,JLEV-1)-ZDPSG(JLON,JLEV)&
     & * (ZDTDT(JLON,JLEV)*(PCP(JLON,JLEV)+TSPHY*ZTDCP)+PT(JLON,JLEV)*ZTDCP)
  ENDDO
ENDDO

DO JLEV=1,KLEV
  DO JLON=KIDIA,KFDIA
    PDIFCQ(JLON,JLEV)=PDIFCQ(JLON,JLEV)-PFCCQL(JLON,JLEV)-PFCCQN(JLON,JLEV)
    PDIFCS(JLON,JLEV)=PDIFCS(JLON,JLEV)&
     & + YDCST%RLVZER*PFCCQL(JLON,JLEV) + YDCST%RLSZER*PFCCQN(JLON,JLEV)
  ENDDO
ENDDO

!-----------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACVPPKF',1,ZHOOK_HANDLE)
END SUBROUTINE ACVPPKF
