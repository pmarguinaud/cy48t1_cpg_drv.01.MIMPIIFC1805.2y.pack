!OPTIONS XOPT(HSFUN)
SUBROUTINE VDFHGHTN (YDEPHLI , YDEPHY  , KIDIA   , KFDIA   , KLON    , KLEV   , KHPBL    , PTMST,&
                   & PTM1    , PQM1    , PLM1    , PIM1   , PAM1,&
                   & PAPHM1  , PAPM1   , PGEOM1  , PGEOH,&
                   & PKMFL   , PKHFL   , PKQFL   , PKHVFL , PMFLX,&
                   & PSLGUH  , PQTUH   , PZINV   , PWUAVG , PZCLDBASE,&
                   & PBIR    , LDNODECP, LDRUNDRY, KPBLTYPE)  
!     ------------------------------------------------------------------

!**   *VDFHGHTN* - DETERMINES THE PBL-HEIGHT AND STRONG UPDRAFT FIELDS
!                  USING A ENTRAINING PARCEL ASCENT METHOD.

!     A.P. SIEBESMA    30/06/99   Original (dry)
!     M. Ko"hler        3/12/2004 Moist Version
!     P. Lopez         02/06/2005 Removed useless option LPHYLIN
!     N.Semane+P.Becht 08/06/2012 Scaling for small planet
!     N. Semane+P.Bechtold     04-10-2012 Add RPLRG/RPLDARE factors for small planet
!     R. El Khatib 06-Oct-2014 R. El Khatib disable a potentially risky vectorization with Intel 
!     F. Vana         17-Dec-2015  Support for single precision 
!     P. Bechtold+I. Sandu 26/03/2018 Pass buoyancy flux, remove mass flux scaling

!     PURPOSE
!     -------

!     DETERMINE PBL HEIGHT AND UPDRAFT FIELDS

!     INTERFACE
!     ---------

!     *VDFHGHTN* IS CALLED BY *VDFMAIN*

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        START POINT
!     *KFDIA*        END POINT
!     *KLEV*         NUMBER OF LEVELS
!     *KLON*         NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):

!     *PTMST*        DOUBLE TIME STEP (SINGLE AT 1TH STEP)        S
!     *PTM1*         TEMPERATURE AT T-1                           K
!     *PQM1*         SPECIFIC HUMUDITY AT T-1                     KG/KG
!     *PLM1*         SPECIFIC CLOUD LIQUID WATER AT T-1           KG/KG
!     *PIM1*         SPECIFIC CLOUD ICE AT T-1                    KG/KG
!     *PAM1*         CLOUD FRACTION AT T-1                        KG/KG
!     *PAPHM1*       PRESSURE AT HALF LEVEL AT T-1                PA
!     *PAPM1*        PRESSURE AT FULL LEVEL AT T-1                PA
!     *PGEOM1*       GEOPOTENTIAL AT T-1                          M2/S2
!     *PGEOH*        GEOPOTENTIAL AT HALF LEVEL                   M2/S2
!     *PKMFL*        SURFACE KINEMATIC MOMENTUM FLUX              M2/S2  
!     *PKHFL*        SURFACE KINEMATIC HEAT FLUX                  K*M/S
!     *PKQFL*        SURFACE KINEMATIC MOISTURE FLUX              M/S
!     *PKHVFL*       SURFACE BUOYANCY FLUX                        K*M/S
!     *PBIR*         BUOYANCY-FLUX INTEGRAL RATIO (-N/P)
!                    USED FOR DECOUPLING CRITERIA

!     INPUT PARAMETERS (LOGICAL):

!     *LDNODECP*     TRUE:  NEVER DECOUPLE
!                    FALSE: MAYBE DECOUPLE
!     *LDRUNDRY*     TRUE:  RUN PARCEL WITHOUT CONDENSATION
!                    FALSE: RUN PARCEL WITH CONDENSATION

!     OUTPUT PARAMETERS (REAL):

!     *PSLGUH*       UPDRAFT GENERALIZED LIQUID STATIC ENERGY (SLG)
!                    AT HALF LEVEL                                M2/S2
!     *PQTUH*        UPDRAFT SPECIFIC TOTAL WATER AT HALF LEVEL   KG/KG
!     *PZINV*        PBL HEIGHT                                   M
!     *PMFLX*        PBL MASS FLUX                                KGM-2S-2


!     OUTPUT PARAMETERS (INTEGER):

!     *KHPBL*        HIGHEST HALF LEVEL BELOW PBL HEIGHT, AND
!                    PBL TOP FULL LEVEL (PZINV IS WITHIN THAT LAYER)
!     *KPBLTYPE*     0: stable PBL
!                    1: dry convective PBL (no cloud below parcel top)
!                    2: cloudy PBL ("stratocumulus")
!                    3: dry PBL below convection ("cumulus")

!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------

USE YOEPHLI   , ONLY : TEPHLI
USE PARKIND1  , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK
USE YOMCST    , ONLY : RG, RD, RCPD, RETV, RLVTT, RLSTT, RTT, RV      , YRCST
USE PARPHY    , ONLY : RKAP
USE YOETHF    , ONLY : R2ES, R3LES, R3IES, R4LES, R4IES, R5LES, R5IES, &
 &                     RVTMP2, R5ALVCP, R5ALSCP, RALVDCP, RALSDCP, RTWAT, RTICE, RTICECU, &
 &                     RTWAT_RTICE_R, RTWAT_RTICECU_R, YRTHF
USE YOMDYNCORE, ONLY : RPLRG, RPLDARE
USE YOEPHY    , ONLY : TEPHY

IMPLICIT NONE

!*         0.1    GLOBAL VARIABLES

TYPE(TEPHLI)      ,INTENT(INOUT) :: YDEPHLI
TYPE(TEPHY)       ,INTENT(INOUT) :: YDEPHY
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KHPBL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOM1(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKMFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKHFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKQFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKHVFL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMFLX(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSLGUH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQTUH(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZINV(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWUAVG(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZCLDBASE(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBIR(KLON) 
LOGICAL           ,INTENT(IN)    :: LDNODECP(KLON) 
LOGICAL           ,INTENT(IN)    :: LDRUNDRY(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KPBLTYPE(KLON) 

!*         0.2    LOCAL VARIABLES

LOGICAL ::            LLDONE(KLON)       , LLDOIT(KLON)       , LLCLOUD(KLON)
REAL(KIND=JPRD) ::    ZUST  (KLON)       , ZWS   (KLON)       , &
                    & ZSLGENH(KLON,0:KLEV),ZQLENH(KLON,0:KLEV), ZQIENH(KLON,0:KLEV), &
                    & ZQTENH(KLON,0:KLEV), ZTVEN(KLON,KLEV)   , ZQTM1 (KLON,KLEV)  , &
                    & ZSLGM1(KLON,KLEV)

REAL(KIND=JPRD) ::    ZWU2H (KLON,0:KLEV), ZWUH  (KLON,0:KLEV), ZQCUH (KLON,0:KLEV), &
                    & ZQUH  (KLON,0:KLEV), ZTUH  (KLON,0:KLEV), &
                    & ZEPS  (KLON,0:KLEV), ZDETR (KLON,0:KLEV)
REAL(KIND=JPRB) ::    ZPH   (KLON)       , ZTTEMP(KLON,KLEV)  , ZQTEMP(KLON,KLEV)

INTEGER(KIND=JPIM) :: IS, JK, JL

REAL(KIND=JPRD) ::    ZQEXC   , ZTEXC   , ZTVUF   , ZQUF    , ZQCUF   , ZQLUF   , &
                    & ZQIUF   , ZQTUF   , ZTUF    , ZSLGUF  , ZMIX(KLON,0:KLEV) , &
                    & ZDZ     , ZBUOF(KLON,KLEV)  , &
                    & ZZ      , ZCPM    , ZTOTW2(KLON)      , ZTOTP(KLON)       , &
                    & ZCONS10 , ZTVMEAN , ZTEMP   , ZRTVLIM , ZRG     , &
                    & ZMFMAX  ,ZMFS(KLON)

!          Model parameters

REAL(KIND=JPRD) ::    ZCM     , ZFACEXC , ZCLDDEPTH,ZTVLIM  , &
                    & ZENTSTPC1         , ZENTSTPC2 ,&
                    & ZW2THRESH         , ZSTABTHRESH       , ZEISTHRESH        , &
                    & ZBIRTHRESH        , ZDETRC

!          VARIABLES FOR CLOUD BASE ESTIMATION

REAL(KIND=JPRD) ::    ZQS(KLON,0:KLEV)  , ZALFAW  , ZFACW   , ZFACI   , ZFAC    , &
                    & ZESDP   , ZCOR    , ZDQSDTEMP(KLON)   , ZSTABILITY(KLON)  , &
                    & ZDZCLOUD(KLON)    , ZEIS(KLON), ZT850 , ZGAMMA850

REAL(KIND=JPRB) ::    ZREPUST 

INTEGER(KIND=JPIM) :: ICLDBASE(KLON), I700(KLON), I850(KLON)
REAL(KIND=JPRB) ::  ZEPSILON
REAL(KIND=JPRB) ::    ZHOOK_HANDLE

#include "surf_inq.h"

#include "cuadjtq.intfb.h"

!DIR$ VFUNCTION EXPHF
#include "fcttre.func.h"


!     ------------------------------------------------------------------

!*         1.     INITIALIZE VARIABLES
!                 --------------------

IF (LHOOK) CALL DR_HOOK('VDFHGHTN',0,ZHOOK_HANDLE)

ZEPSILON    = 100._JPRB*EPSILON(1._JPRB) ! numerical security constant
ZENTSTPC1   = 0.8_JPRB     ! prefactor of the entrainment function
ZENTSTPC2   = 2.E-4_JPRB*RPLRG   ! constant entrainment function
ZDETRC      = 0.0003_JPRB*RPLRG  ! detrainment (1/z)
ZCM         = 0.1_JPRB     ! prefactor of the mass flux initialization
                           ! (consistency with vdfexcu!)
ZFACEXC     = 2.0_JPRB     ! prefactor of T, q and w,up excess values at L60
!ZW2THRESH  = -1._JPRB     ! threshold parcel vertical velocity squared [m2/s2]
ZW2THRESH   = 0._JPRB      ! threshold parcel vertical velocity squared [m2/s2]
ZCLDDEPTH   = 2000._JPRB/RPLRG ! threshold cloud thickness for stcu/cu transition [m]
ZSTABTHRESH = 20._JPRB     ! threshold stability (Klein & Hartmann criteria) [K]
ZEISTHRESH  = 7.0_JPRB     ! threshold stability (Wood & Bretherton) [K]
ZBIRTHRESH  = 0.1_JPRB     ! threshold BIR (TKE decoupling criteria) [1]
ZTVLIM      = 0.1_JPRB     ! cloud fraction limit in Tv,env calculation

CALL SURF_INQ(YDEPHY%YSURF,PREPUST=ZREPUST)

! optimization
ZRG         = 1.0_JPRB/RG
ZRTVLIM     = 1.0_JPRB/ZTVLIM

!DIR$ LOOP_INFO EST_TRIPS(16)
DO JL=KIDIA,KFDIA
  PZCLDBASE(JL)  = -100._JPRB  ! default value: no PBL cloud
  PWUAVG(JL)     = 0.0_JPRB
  KHPBL(JL)      = 0
  KPBLTYPE(JL)   = 0
  ICLDBASE(JL)   = 0           ! default value: no PBL cloud
  LLCLOUD(JL)    = .FALSE.
ENDDO
DO JK=0,KLEV
!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA
    PSLGUH(JL,JK)= 0.0_JPRB
    PQTUH(JL,JK) = 0.0_JPRB
    PMFLX(JL,JK) = 0.0_JPRB
    ZTUH(JL,JK)  = 0.0_JPRB
    ZQUH(JL,JK)  = 0.0_JPRB
    ZQCUH(JL,JK) = 0.0_JPRB
    ZEPS(JL,JK)  = 0.0_JPRB
    ZDETR(JL,JK) = 0.0_JPRB
    ZWU2H(JL,JK) = 0.0_JPRB
  ENDDO
ENDDO

!     -----------------------------------------------------------------

!*         2.     PREPARE FIELDS ON HALF LEVELS BY LINEAR INTERPOLATION
!*                OF CONSERVED VARIABLES
!                 -----------------------------------------------------

!*         2.1  full level cpm, slg, qt and Tv

  DO JK=1,KLEV
!DIR$ LOOP_INFO EST_TRIPS(16)
    DO JL=KIDIA,KFDIA
      ZCPM          = RCPD * ( 1.0_JPRB + RVTMP2 * PQM1(JL,JK) )
      ZSLGM1(JL,JK) = ZCPM * PTM1(JL,JK) + PGEOM1(JL,JK) &
                  & - RLVTT * PLM1(JL,JK) - RLSTT * PIM1(JL,JK)  
      ZQTM1 (JL,JK) = PQM1(JL,JK) + PLM1(JL,JK) + PIM1(JL,JK)

!          parcel goes through cloud portion of environment
!          (added ql loading; ql,cld=ql,mean/fc; qv = qsat) 
!          safety: fc>0.1; linear interpolation between overcast 
!                  and cloudy portion for 0<fc<0.1
!                  guaranteed to be < tv from mean conditions

!          grid box mean virtual effect
      ZTVMEAN       = PTM1(JL,JK) * ( 1.0_JPRB + RETV * PQM1(JL,JK) &
                  & - PLM1(JL,JK) - PIM1(JL,JK) )       !qli loading  
      ZTVEN(JL,JK) = ZTVMEAN
    ENDDO
  ENDDO


!*         2.2  half-level environment interpolation (qt, ql, qi, slg)
!*              attention:  not good to interpolate everything independently
!*              better:     interpolate conserved variables and derive rest!!!

  JK=KLEV-1
!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA
    ZQTENH(JL,JK) = ( ZQTM1(JL,JK+1) *(PGEOH(JL,JK-1)-PGEOH(JL,JK  )) &
                & +   ZQTM1(JL,JK)   *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(PGEOH(JL,JK-1)-PGEOH(JL,JK+1))
    ZQLENH(JL,JK) = ( PLM1(JL,JK+1)  *(PGEOH(JL,JK-1)-PGEOH(JL,JK  )) &
                & +   PLM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(PGEOH(JL,JK-1)-PGEOH(JL,JK+1))
    ZQIENH(JL,JK) = ( PIM1(JL,JK+1)  *(PGEOH(JL,JK-1)-PGEOH(JL,JK  )) &
                & +   PIM1(JL,JK)    *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(PGEOH(JL,JK-1)-PGEOH(JL,JK+1))
    ZSLGENH(JL,JK)= ( ZSLGM1(JL,JK+1)*(PGEOH(JL,JK-1)-PGEOH(JL,JK  )) &
                & +   ZSLGM1(JL,JK)  *(PGEOH(JL,JK  )-PGEOH(JL,JK+1)) &
                &   )                /(PGEOH(JL,JK-1)-PGEOH(JL,JK+1))
  ENDDO


!     -----------------------------------------------------------------

!*         3.     CALCULATE EXCESS VALUES AT LOWEST HALF MODEL LEVEL
!                 --------------------------------------------------

!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA
 

!*         3.0  surface initializaion

    ZUST  (JL)  = MAX( SQRT(PKMFL(JL)), ZREPUST )               !u* (repust=10e-4)
    PKHVFL(JL)  = (PKHFL(JL) + RETV * PTM1(JL,KLEV) * PKQFL(JL))/(RPLRG*RPLDARE) !w'theta,v'
    KHPBL (JL)  = KLEV

    IF ( PKHVFL(JL) >= 0.0_JPRB ) THEN

      LLDONE(JL)       = .TRUE.
      PZINV (JL)       = 0.0_JPRB   ! stable, therefore 0 PBL depth, no updraft

    ELSE

      LLDONE(JL)       = .FALSE.
      ZQUH  (JL,KLEV)  = 0.0_JPRB
      ZQCUH (JL,KLEV)  = 0.0_JPRB
      PQTUH (JL,KLEV)  = 0.0_JPRB
      PSLGUH(JL,KLEV)  = 0.0_JPRB
      ZTUH  (JL,KLEV)  = 0.0_JPRB
      ZBUOF (JL,KLEV)  = 0.0_JPRB


!*         3.1  sigma-w-L60 (ignore 1-z/zi term)

      ZWS(JL)          = 1.2_JPRB &
       & * ( ZUST(JL)**3 &
       & - 1.5_JPRB * RKAP * PKHVFL(JL) * PGEOH(JL,KLEV-1) / PTM1(JL,KLEV-1) &
       & ) ** ( 1.0_JPRB/3._JPRB )                         ! Kolmogorov 1/3-power


!*         3.2  excess values
      
      ZWU2H(JL,KLEV-1) = ( ZFACEXC * ZWS(JL) )**2         ! parcel wup initialization
      ZTEXC            = (- ZFACEXC * PKHFL(JL)) / (ZWS(JL)*RPLRG*RPLDARE) 
      ZQEXC            = (- ZFACEXC * PKQFL(JL)) / (ZWS(JL)*RPLRG*RPLDARE) 
      ZTEXC            = MAX(ZTEXC, 0.0_JPRB)
      ZQEXC            = MAX(ZQEXC, 0.0_JPRB)
      ZWU2H(JL,KLEV-1) = MIN(ZWU2H(JL,KLEV-1),100.0_JPRB) !10 m/s  limit (pure safety)
      ZTEXC            = MIN(ZTEXC           , 10.0_JPRB) !10 K    limit
      ZQEXC            = MIN(ZQEXC           , 0.01_JPRB) !10 g/kg limit
      PQTUH(JL,KLEV-1) = ZQTENH(JL,KLEV-1) + ZQEXC
      ZQCUH(JL,KLEV-1) = ZQLENH(JL,KLEV-1) + ZQIENH(JL,KLEV-1)
      ZQUH (JL,KLEV-1) = PQTUH(JL,KLEV-1)  - ZQCUH(JL,KLEV-1)
      ZCPM             = RCPD * ( 1.0_JPRB + RVTMP2 * ZQUH(JL,KLEV-1) )
      PSLGUH(JL,KLEV-1)= ZSLGENH(JL,KLEV-1) + ZCPM * ZTEXC
      ZTUH (JL,KLEV-1) = ( PSLGUH (JL,KLEV-1) - PGEOH(JL,KLEV-1) &
                       & + RLVTT*ZQLENH(JL,KLEV-1) + RLSTT*ZQIENH(JL,KLEV-1) &
                       & ) / ZCPM
    ENDIF

  ENDDO


!     -----------------------------------------------------------------

!*         4.     VERTICAL ASCENT UNTIL VELOCITY BECOMES NEGATIVE
!                 -----------------------------------------------

  DO JK=KLEV-2,1,-1
    IS=0
! Disable vectorization with Intel compiler as it may conflict with fp speculation,
! leading to non-reproducible results.
!DEC$ NOVECTOR
    DO JL=KIDIA,KFDIA
      IF (.NOT. LLDONE(JL)) THEN
        IS            = IS+1


!*         4.1  parcel entrainment

        ZWUH(JL,JK+1) = SQRT( MAX( ZWU2H(JL,JK+1), 0.01_JPRB) ) ! w,up > 0.1 m/s (safety)
        ZEPS(JL,JK+1) = ZENTSTPC2 &
                    & + ZENTSTPC1 / ( PGEOM1(JL,JK+1)*ZRG )
        ZDZ           = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG


!*         4.2  ascent of slg and qt (exact)
 
        ZMIX(JL,JK+1) = EXP( - ZDZ * ZEPS(JL,JK+1) )
        PQTUH(JL,JK)  = ( PQTUH (JL,JK+1) - ZQTM1 (JL,JK+1) ) * ZMIX(JL,JK+1) &
                    & + ZQTM1 (JL,JK+1)
        PSLGUH(JL,JK) = ( PSLGUH(JL,JK+1) - ZSLGM1(JL,JK+1) ) * ZMIX(JL,JK+1) &
                    & + ZSLGM1(JL,JK+1)


!*         4.3  condensation - diagnose T, qv, ql
!*         (Taylor, becomes more inaccurate for states far from q*
!*          -> double condensation step)

!          cuadjtq initialization (assume qv=qt)

        ZQUH(JL,JK)   = PQTUH(JL,JK)
        ZQCUH(JL,JK)  = 0.0_JPRB

!          cuadjtq initialization (assume qv=qv(jk+1) for speed)

!       IF ( ZQUH(JL,JK+1) < PQTUH(JL,JK) ) THEN
!         ZQUH(JL,JK) = ZQUH(JL,JK+1)
!         ZQCUH(JL,JK)= PQTUH(JL,JK) - ZQUH(JL,JK+1)
!       ENDIF

        ZCPM          = RCPD * ( 1.0_JPRB + RVTMP2 * ZQUH(JL,JK) )
        ZTUH(JL,JK)   = ( PSLGUH(JL,JK) - PGEOH(JL,JK) + RLVTT*ZQCUH(JL,JK) ) &
                    & / ZCPM      ! assume liquid phase!
        ZPH(JL)       = PAPHM1(JL,JK)
        ZQTEMP(JL,JK) = ZQUH(JL,JK)
        ZTTEMP(JL,JK) = ZTUH(JL,JK)
      ENDIF
      LLDOIT(JL)      = .NOT. LLDONE(JL)
      IF ( LDRUNDRY(JL) ) THEN    ! condensation not done for dry updraft
        LLDOIT(JL)    = .FALSE.
      ENDIF
    ENDDO

    CALL CUADJTQ &
     & ( YRTHF, YRCST, YDEPHLI,  KIDIA,    KFDIA,    KLON,     KLEV,    JK,&
     &   ZPH,      ZTTEMP,   ZQTEMP,   LLDOIT,  4)  

!DIR$ LOOP_INFO EST_TRIPS(16)
    DO JL=KIDIA,KFDIA
      IF ( LLDOIT(JL) ) THEN
        IF ( ZQTEMP(JL,JK) < PQTUH(JL,JK) ) THEN !allow evaporation up to qt
          ZQUH(JL,JK) = ZQTEMP(JL,JK)
          ZQCUH(JL,JK)= PQTUH(JL,JK) - ZQUH(JL,JK)
          ZTUH(JL,JK) = ZTTEMP(JL,JK)
        ELSE                          !case where qv(initial)<qt but qv(final)>qt
          ZQUH(JL,JK) = PQTUH(JL,JK)  !(unusual!)
          ZQCUH(JL,JK)= 0.0_JPRB
          ZCPM        = RCPD * ( 1.0_JPRB + RVTMP2 * ZQUH(JL,JK) )
          ZTUH(JL,JK) = ( PSLGUH(JL,JK) - PGEOH(JL,JK) + RLVTT*ZQCUH(JL,JK) ) &
                    & / ZCPM
        ENDIF
      ENDIF
    ENDDO


!*         4.4  ESTIMATION OF PARCEL CLOUD BASE (INTERPOLATED)

!DIR$ LOOP_INFO EST_TRIPS(16)
    DO JL=KIDIA,KFDIA
      IF ( ZQCUH(JL,JK) > 0.0_JPRB  .AND.  .NOT. LLCLOUD(JL) ) THEN
        LLCLOUD(JL)   = .TRUE.
        ICLDBASE(JL)  = JK

!... calculate qsat (mixed phase)
        ZQS(JL,JK)    = FOEEWM(REAL(ZTUH(JL,JK),JPRB))/PAPHM1(JL,JK)
        ZQS(JL,JK)    = MIN(0.5_JPRB,ZQS(JL,JK))
        ZQS(JL,JK)    = ZQS(JL,JK)/(1.0_JPRB-RETV*ZQS(JL,JK))

!... calculate dqs/dT correction factor (mixed phase)
        ZALFAW        = FOEALFA(REAL(ZTUH(JL,JK),JPRB))
        ZFACW         = R5LES/((ZTUH(JL,JK)-R4LES)**2)
        ZFACI         = R5IES/((ZTUH(JL,JK)-R4IES)**2)
        ZFAC          = ZALFAW*ZFACW+(1.0_JPRB-ZALFAW)*ZFACI
        ZESDP         = FOEEWM(REAL(ZTUH(JL,JK),JPRB))/PAPHM1(JL,JK)
        ZCOR          = 1.0_JPRB/(1.0_JPRB-RETV*ZESDP)
        ZDQSDTEMP(JL) = ZFAC*ZCOR*ZQS(JL,JK)

        PZCLDBASE(JL) = PGEOH(JL,ICLDBASE(JL))*ZRG &
                    & - ZQCUH(JL,JK) * ( RCPD*ZRG / ZDQSDTEMP(JL) + RLVTT*ZRG )
        IF ( ICLDBASE(JL) < KLEV ) THEN
          PZCLDBASE(JL) = MAX( PZCLDBASE(JL), PGEOH(JL,ICLDBASE(JL)+1)*ZRG )
        ELSE
          PZCLDBASE(JL) = MAX( PZCLDBASE(JL), 0.0_JPRB )
        ENDIF
      ENDIF


!*         4.5  buoyancy
!*         (at full level k+1 from interpolation of slg, q, qc)

      IF ( .NOT. LLDONE(JL) ) THEN

        ZSLGUF        = 0.5_JPRB * ( PSLGUH(JL,JK) + PSLGUH(JL,JK+1) )
        ZQTUF         = 0.5_JPRB * ( PQTUH (JL,JK) + PQTUH (JL,JK+1) )
        ZQCUF         = 0.5_JPRB * ( ZQCUH (JL,JK) + ZQCUH (JL,JK+1) )
        IF ( JK == ICLDBASE(JL) ) THEN
          ZQCUF = ZQCUF * ( PGEOH(JL,JK)*ZRG - PZCLDBASE(JL)      ) &
                &       / ( PGEOH(JL,JK)*ZRG - PGEOH(JL,JK+1)*ZRG )
        ENDIF
        ZQUF          = ZQTUF - ZQCUF
        ZCPM          =  RCPD  * ( 1.0_JPRB + RVTMP2 * ZQUF )
        ZTUF          = ( ZSLGUF - PGEOM1(JL,JK+1) & ! preliminary estimate:
                    & + RLVTT * ZQCUF ) / ZCPM       ! all liquid 
        ZALFAW        = FOEALFA(REAL( ZTUF,JPRB) )
        ZQLUF         = ZALFAW            * ZQCUF
        ZQIUF         = (1.0_JPRB-ZALFAW) * ZQCUF
        ZTUF          = ( ZSLGUF - PGEOM1(JL,JK+1) &
                    & + RLVTT * ZQLUF + RLSTT * ZQIUF ) / ZCPM  
        ZTVUF         = ZTUF * ( 1.0_JPRB + RETV * ZQUF - ZQCUF )  
        ZBUOF(JL,JK+1)= RG * ( ZTVUF - ZTVEN(JL,JK+1) ) / ZTVEN(JL,JK+1)


!*         4.6  kinetic energy equation (exact)

        ZDZ           = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG
        ZTEMP         = ZBUOF(JL,JK+1)/ZEPS(JL,JK+1)
        ZWU2H(JL,JK)  = ( ZWU2H(JL,JK+1) - ZTEMP ) * ZMIX(JL,JK+1)**2 + ZTEMP


!*         4.7  inversion height at w=0  (lin. interpolation in w^2)
        
        IF ( ZWU2H(JL,JK) < 0.0_JPRB  .AND.  ZWU2H(JL,JK+1) > 0.0_JPRB ) THEN
          KHPBL(JL)   = JK+1
          PZINV(JL)   = PGEOH(JL,JK+1) * ZRG &
                    & + ZDZ * ZWU2H(JL,JK+1) / ( ZWU2H(JL,JK+1) - ZWU2H(JL,JK) )
        ENDIF

        IF ( ZWU2H(JL,JK) < ZW2THRESH ) THEN   !allow parcel to overcome layers 
          LLDONE(JL)  = .TRUE.                 !with small negative kin. energy
        ENDIF                                  !but remember last w(z)=0 (z=PZINV)

      ENDIF
    ENDDO
    IF (IS == 0) EXIT
  ENDDO


!     -----------------------------------------------------------------

!*         5.     PBL TYPE (0, 1, 2)
!                 ------------------

!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA
    IF ( PZCLDBASE(JL) > PZINV(JL) .OR. ICLDBASE(JL) == 0 ) THEN
      PZCLDBASE(JL) = -100._JPRB          !default value (no cloud)
      KPBLTYPE(JL)  = 1                   !dry convective PBL
      ZDZCLOUD(JL)  = 0.0_JPRB            !cloud thickness
    ELSE
      KPBLTYPE(JL)  = 2                   !cloudy PBL
      ZDZCLOUD(JL)  = PZINV(JL) - PZCLDBASE(JL) !cloud thickness
    ENDIF
    IF ( PZINV(JL) == 0.0_JPRB ) THEN
      KPBLTYPE(JL)  = 0                   !dry stable PBL
    ENDIF
  ENDDO


!     -----------------------------------------------------------------

!*         6.     PARCEL MASS FLUX
!                 ----------------

!*         6.1  mass flux initialization
!               (updraft fraction * scale factor * sigma-w at L60 * rho)
!               (ignore u* term following Anton's suggestion
!               - consistency with similarity theory close to the surface)

!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA
    PMFLX(JL,:)        = 0.0_JPRB
    IF ( PZINV(JL) > ZEPSILON ) THEN
      ZZ               = PGEOH(JL,KLEV-1)*ZRG
      PMFLX(JL,KLEV-1) = ZCM * 1.2_JPRB &
       & * ( - 1.5_JPRB * RKAP * PKHVFL(JL) / PTM1(JL,KLEV) * RG *ZZ ) &
       & ** (1.0_JPRB/3._JPRB) 
      PMFLX(JL,KLEV-1) = PMFLX(JL,KLEV-1) &
       & * SQRT( 1.0_JPRB - ZZ / PZINV(JL))
      PMFLX(JL,KLEV-1) = PMFLX(JL,KLEV-1) &
       & * PAPHM1(JL,KLEV-1)/RD/PTM1(JL,KLEV)          !rho  
    ENDIF
  ENDDO


!*         6.2  mass flux entrainment (forward, implicit or exact;
!               forward and exact is OK, implicit is instable for eps*dz>1)

  DO JK=KLEV-2,1,-1
!DIR$ LOOP_INFO EST_TRIPS(16)
    DO JL=KIDIA,KFDIA
      IF ( PZINV(JL) > ZEPSILON ) THEN
        IF ( PGEOH(JL,JK-1)*ZRG < PZINV(JL) ) THEN !exclude entrainment layer
          IF ( PGEOM1(JL,JK+1)*ZRG > PZCLDBASE(JL) ) THEN
            ZDETR(JL,JK+1) = ZDETRC                !detrainment above cloud base
          ELSE
            ZDETR(JL,JK+1) = 0.0_JPRB
          ENDIF
          ZDZ          = (PGEOH(JL,JK) - PGEOH(JL,JK+1))*ZRG
          PMFLX(JL,JK) = PMFLX(JL,JK+1) *      EXP( ZDZ * ( ZEPS(JL,JK+1) - ZDETR(JL,JK+1) ) ) !exact


!*         6.3  limit mass flux covering 50% area (M<rho*w,up*0.5)
!               (detrainment is initiated if strong w,up slowdown)

          PMFLX(JL,JK) = MIN( PMFLX(JL,JK), ZWUH(JL,JK) * 0.5_JPRB &
                           & *PAPHM1(JL,JK)/RD/PTM1(JL,JK) )         !rho
        ENDIF
      ENDIF
    ENDDO
  ENDDO


!*         6.4  mass flux limit according to CFL criteria
!               (reduce M profile uniformly by maximum excess)

  ZCONS10 = 1.0_JPRB/(RG*PTMST)
!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA
    ZMFS(JL) = 1.0_JPRB  ! mass flux scaling value (reduction)
  ENDDO
  DO JK=1,KLEV-1
!DIR$ LOOP_INFO EST_TRIPS(16)
    DO JL=KIDIA,KFDIA
      IF ( PGEOH(JL,JK-1)*ZRG < PZINV(JL) ) THEN
        ZMFMAX = (PAPM1(JL,JK+1)-PAPM1(JL,JK)) * ZCONS10
        ZMFMAX = 2.0_JPRB * ZMFMAX

        IF ( PMFLX(JL,JK) > ZMFMAX ) THEN
          ZMFS(JL) = MIN(ZMFS(JL),ZMFMAX/PMFLX(JL,JK))
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  DO JK=1,KLEV
!DIR$ LOOP_INFO EST_TRIPS(16)
    DO JL=KIDIA,KFDIA
      IF ( PGEOH(JL,JK-1)*ZRG < PZINV(JL) ) THEN
        PMFLX(JL,JK) = PMFLX(JL,JK)*ZMFS(JL)
      ENDIF
    ENDDO
  ENDDO

!     -----------------------------------------------------------------

!*         7.     STRATOCUMULUS - SHALLOW CUMULUS CRITERIA: 
!*                * CLOUD THICKNESS = 1000M
!*                * STABILITY = 15K
!*                * TKE DECOUPLING, BIR = 0.1
!                 -----------------------------------------

!*         7.1  stability criteria == theta(700hPa) - theta(sfc)

!          find index I700 of pressure closest to 700hPa

!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA
    I700(JL) = 0
    I850(JL) = 0
  ENDDO
  DO JK=1,KLEV
!DIR$ LOOP_INFO EST_TRIPS(16)
    DO JL=KIDIA,KFDIA
      IF ( I700(JL) == 0  .AND.  PAPM1(JL,JK) > 70000.0_JPRB ) THEN
        I700(JL) = JK
      ENDIF
    ENDDO
!DIR$ LOOP_INFO EST_TRIPS(16)
    DO JL=KIDIA,KFDIA
      IF ( I850(JL) == 0  .AND.  PAPM1(JL,JK) > 85000.0_JPRB ) THEN
        I850(JL) = JK
      ENDIF
    ENDDO
  ENDDO
!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA
    IF ( I700(JL) > 1 ) THEN
       IF ( ABS(PAPM1(JL,I700(JL)-1)-70000.0_JPRB) <  &
          & ABS(PAPM1(JL,I700(JL))  -70000.0_JPRB)   ) THEN
          I700(JL) = I700(JL)-1
       ENDIF
    ENDIF
    IF ( I850(JL) > 1 ) THEN
       IF( ABS(PAPM1(JL,I850(JL)-1)-85000.0_JPRB) <  &
         & ABS(PAPM1(JL,I850(JL))  -85000.0_JPRB)   ) THEN 
          I850(JL) = I850(JL)-1
       ENDIF
    ENDIF
    IF ( I700(JL) > 0 ) THEN
      ZSTABILITY(JL) = PTM1(JL,I700(JL)) * ( 1.0E5_JPRB/PAPM1(JL,I700(JL)) ) ** (RD/RCPD) &
                   & - PTM1(JL,KLEV)     * ( 1.0E5_JPRB/PAPM1(JL,KLEV) )     ** (RD/RCPD)  
    ELSE
      ZSTABILITY(JL) = 0.0_JPRB
    ENDIF
  ENDDO


!*         7.2  Estimated Inversion Strength (EIS) criteria from Wood & Bretherton (2006)

!          find index I700 of pressure closest to 700hPa

!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA
    IF ( I700(JL) > 0  .AND.  I850(JL) > 0  .AND.  KPBLTYPE(JL) == 2 ) THEN
      ZT850      = ( PTM1(JL,KLEV) + PTM1(JL,I700(JL)) ) / 2 
      JK         = I850(JL)
!          qsat (full level)
      ZQS(JL,JK) = FOEEWM(REAL(ZT850,JPRB))/PAPM1(JL,JK)
      ZQS(JL,JK) = MIN(0.5_JPRB,ZQS(JL,JK))
      ZQS(JL,JK) = ZQS(JL,JK)/(1.0_JPRB-RETV*ZQS(JL,JK))

      ZGAMMA850  = RG/RCPD* (1 - ( 1 + RLVTT   *ZQS(JL,JK) / (     RD*ZT850   ) ) &
                             & / ( 1 + RLVTT**2*ZQS(JL,JK) / (RCPD*RV*ZT850**2) ) )
      ZEIS(JL)   = ZSTABILITY(JL) - ZGAMMA850 * ( PGEOM1(JL,I700(JL))/RG - PZINV(JL) )
    ELSE
      ZEIS(JL)   = 0.0_JPRB
    ENDIF
  ENDDO


!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA

!*         7.3  decoupling criteria (options)

    IF ( .NOT. LDNODECP(JL)  .AND.  KPBLTYPE(JL) == 2 ) THEN

!..........stability criteria (Klein & Hartmann 1993)
!     IF ( ZSTABILITY(JL) < ZSTABTHRESH ) THEN 

!..........stability criteria (Wood & Bretherton 2006)
      IF ( ZEIS(JL) < ZEISTHRESH ) THEN 

!..........cloud thickness criteria
!     IF ( ZDZCLOUD(JL) > ZCLDDEPTH ) THEN

!..........TKE decoupling (or cloud thickness criteria for safety)
!     IF ( PBIR(JL) > ZBIRTHRESH .OR. ZDZCLOUD(JL) > ZCLDDEPTH ) THEN

!..........always decouple
!     IF ( .TRUE. ) THEN

        KPBLTYPE(JL) = 3                  !decouple - PBL type 3

      ENDIF
    ENDIF

    IF ( KPBLTYPE(JL) == 3 ) THEN

!     KHPBL(JL)     = ICLDBASE(JL) + 1    !level below updraft cloud
!     PZINV(JL)     = PZCLDBASE(JL)       !cloud base
      KHPBL(JL)     = ICLDBASE(JL)        !lowest level of updraft cloud
      PZINV(JL)     = PGEOH(JL,ICLDBASE(JL))*ZRG !lowest cloud level
      PZCLDBASE(JL) = -100._JPRB          !default value (no cloud)
    ENDIF

  ENDDO


!     -----------------------------------------------------------------

!*         8.     MEAN PARCEL W (USED IN CLOUD VARIANCE DISSIPATION)
!                 --------------------------------------------------

!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA
    ZTOTW2(JL) = 0.0_JPRB
    ZTOTP(JL)  = 0.0_JPRB
  ENDDO
  DO JK=1,KLEV-1
!DIR$ LOOP_INFO EST_TRIPS(16)
    DO JL=KIDIA,KFDIA
      IF ( JK >= KHPBL(JL) ) THEN   ! only up to PBL top (to cloud base if decoupled)
        ZTOTW2(JL) = ZTOTW2(JL) + ( PAPM1(JL,JK+1) - PAPM1(JL,JK) ) * MAX(ZWU2H(JL,JK),0.0_JPRB)
        ZTOTP(JL)  = ZTOTP(JL)  + ( PAPM1(JL,JK+1) - PAPM1(JL,JK) )
      ENDIF
    ENDDO
  ENDDO
!DIR$ LOOP_INFO EST_TRIPS(16)
  DO JL=KIDIA,KFDIA
    IF (ZTOTP(JL) > 0.0_JPRB) THEN  ! protect against single layer PBL
      PWUAVG(JL) = SQRT ( ZTOTW2(JL) / ZTOTP(JL)  )
    ELSE
      PWUAVG(JL) = SQRT ( MAX(ZWU2H(JL,KLEV-1),0.0_JPRB) )
    ENDIF
  ENDDO


!     -----------------------------------------------------------------

!*         9.     CLEAN-UP (no impact!)
!                 --------

  DO JK=0,KLEV
!DIR$ LOOP_INFO EST_TRIPS(16)
    DO JL=KIDIA,KFDIA
      IF ( JK <= KHPBL(JL) ) THEN  ! entrainment layer: M=0 (K does entrainment)
        PMFLX(JL,JK)  = 0.0_JPRB   ! important for shallow convection criteria
        PQTUH(JL,JK)  = 0.0_JPRB
        PSLGUH(JL,JK) = 0.0_JPRB
      ENDIF
    ENDDO
  ENDDO

IF (LHOOK) CALL DR_HOOK('VDFHGHTN',1,ZHOOK_HANDLE)
END SUBROUTINE VDFHGHTN
