SUBROUTINE ACTRGRC( YDECUMF,YDPHY0,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
!--------------------------------------------------------------------------
! - INPUT 2D.
                   & PAPHI,PAPHIF,PAPRS,PAPRSF,PCP,PDELP, &
                   & PQ,PQC,PQSAT,PQW, PT,PTVE, PTW, &
                   & PIFHCS, &
! - INPUT 1D.
                   & PIB, &
! - INPUT/OUTPUT 2D
                   & KNASC, PQU, PSU, PTU, PTVU, &
! - OUTPUT 1D.
                   & KNND,KNNS, KUSLT, KLFCT, PPLCL, PQUSL, PTLCL, PTUSL)
!**** *ACTRGRC * USL TRIGGERING CALCULATION FOR ACCSU
!                TO BE CALLED BY ACCSU
!--------------------------------------------------------------------------
! -   INPUT ARGUMENTS.
!     ----------------

! - PHYSICS DIMENSIONNING PARAMETERS

! KIDIA      : INDICE DE DEPART DES BOUCLES VECTORISEES SUR L'HORIZONT..
! KFDIA      : INDICE DE FIN DES BOUCLES VECTORISEES SUR L'HORIZONTALE.
! KLON       : DIMENSION HORIZONTALE DES TABLEAUX.
! KTDIA      : INDICE DE DEPART DES BOUCLES VERTICALES (1 EN GENERAL).
! KLEV       : DIMENSION VERTICALE DES TABLEAUX "FULL LEVEL".

! - PHYSICS VARIABLES (ALPHABETICALLY IN EACH CATEGORY)

! - 2D (0:KLEV) .
! PAPHI      : HALF LEVEL GEOPOTENTIAL
! PAPRS      : HALF LEVEL PRESSURE

! - 2D (1:KLEV) .
! PAPHIF     : FULL LEVEL GEOPOTENTIAL
! PAPRSF     : FULL LEVEL PRESSURE
! PCP        : CHALEUR MASSIQUE A PRESSION CONSTANTE DE L'AIR.
! PDELP      : PRESSURE DIFFERENCE OVER THE LAYER
! PQ         : MICROPHYSICAL WATER VAPOUR SPECIFIC HUMIDITY
! PQC        : MICROPHYSICAL RESOLVED CLOUD CONDENSATE (LIQUID+ICE)
! PQSAT      : SATURATION SPECIFIC MOISTURE
! PQW        : WET BULB SPECIFIC MOISTURE
! PT         : TEMPERATURE.
! PTVE       : ENVIRONMENT VIRTUAL TEMPERATURE
! PTW        : WET BULB TEMPERATURE.
! PIFHCS     : RESOLVED CONDENSATION D(Fcs)*L/c_p as COMPUTED IN ACCSU
! - 1D
! PIB        : INTEGRATED BUOYANCY (to skip columns with no CAPE)
!
! -   INPUT/OUTPUT ARGUMENTS
!     ----------------------
!--------------------------------------------------------------------------
! -   OUTPUT ARGUMENTS
!     ----------------
! - 2D (1:KLEV) (ONLY AT CERTAIN LEVELS)
! KNASC  : ASCENT MEMBERSHIP INDICATOR
! PSU    : DRY STATIC HEAT PERTURBATION OF THE CLOUDY ASCENT 
! PTU    : TEMPERATURE PERTURBATION OF THE CLOUDY ASCENT.
! PQU    : WATER VAPOUR SPECIFIC HUMIDITY PERTURBATION OF THE CLOUDY ASCENT.
!
!--------------------------------------------------------------------------
! -   DIAGNOSTICS
!     -----------
! 
!--------------------------------------------------------------------------
!      BASE VERSION:  ACUSLTRIG MAY 10 2011, L. GERARD
!                     ACUSLTRIG60 FEBRUARY 2014, L. GERARD
!                     ACUSLTRIG50 JUNE 2014, L. GERARD
!                     ACTRGRC MARCH 2016, L. GERARD
!      MODIFICATIONS
!      -------------
!--------------------------------------------------------------------------
USE PARKIND1, ONLY: JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE YOMCST   , ONLY : RATM, RG  ,RD       ,RV       ,RCPD     ,RDT     ,&
            &RCPV     ,RETV     ,RCW      ,RCS      ,RKAPPA   ,RLVTT   ,&
            &RLSTT    ,RTT      ,RALPW    ,RBETW    ,RGAMW    ,&
            &RALPS    ,RBETS    ,RGAMS    ,RALPD    ,RBETD    ,&
            &RGAMD    ,RLVZER   ,RLSZER
USE YOMPHY0  , ONLY : TPHY0
USE YOECUMF  , ONLY : TECUMF

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
TYPE(TECUMF)      ,INTENT(IN) :: YDECUMF
TYPE(TPHY0)       ,INTENT(IN) :: YDPHY0
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KTDIA

REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHI(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHIF(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQC(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSAT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQW(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTVE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTW(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIFHCS(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PIB(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PPLCL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PQUSL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PTUSL(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PTLCL(KLON) 

REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PQU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PSU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PTU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PTVU(KLON,KLEV)

INTEGER(KIND=JPIM),INTENT(INOUT)   :: KNASC(KLON,KLEV) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KNND(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KNNS(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KUSLT(KLON) 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLFCT(KLON) 



! 2D REAL LOCALS

REAL(KIND=JPRB) :: ZA12(KLON,KLEV), ZA13(KLON,KLEV), &
    &ZQM(KLON,KLEV),ZQCM(KLON,KLEV),ZTHVE(KLON,KLEV),&
    &ZTM(KLON,KLEV),ZTHESE(KLON,KLEV)

! 1D REAL LOCALS
REAL(KIND=JPRB) :: ZCONTINUE(KLON),  ZS14(KLON),  &
    &ZS6(KLON), ZDTVLCL(KLON),ZDTVFC(KLON),ZKICKRC(KLON),&
    &ZSFCS(KLON),ZSUM(KLON), ZSUM1(KLON), ZSDP(KLON), &
    &ZSQ(KLON),ZSQC(KLON),ZST(KLON), ZTHEQSL(KLON), ZTHVSL(KLON),&
    &ZTCIN(KLON),ZTVKICK(KLON), &
    &ZKICKSC(KLON),ZDPHISC(KLON)

! REAL SCALARS
REAL(KIND=JPRB) :: ZB, ZBAS, ZBLUE,  ZCPM,&
    & ZDELTA, ZDEEPER, ZDP, ZDP0, ZDP1,ZDP2, ZDPHI,&
    & ZESL,&
    & ZEW, ZEXP1,ZEXP2, ZHPHI, ZIALBU, ZINAC9, ZINCIN,&
    & ZLIQ,ZLM, ZNO,&
    & ZOK, ZOK0, ZOK1, ZOK2, ZOK3, ZOK4, ZOK5,&
    & ZPLCL,ZQSL, ZR, ZTSL,ZTV,&
    & ZZ,ZZ2 !, Z2BSA
! real Local parameters
REAL(KIND=JPRB) :: ZCPVMD, ZDPMIX2,ZDPMX2I,&
             & ZKAPPAI, ZRDSRG, ZRDSRV, ZTHRLCL
REAL (KIND=JPRB), PARAMETER :: ZEPS1=1.E-10_JPRB
REAL (KIND=JPRB), PARAMETER :: ZEPSQV=1.E-10_JPRB
REAL (KIND=JPRB), PARAMETER :: ZQSMAX=0.2_JPRB
REAL (KIND=JPRB), PARAMETER :: Z3I=1._JPRB/3._JPRB
REAL (KIND=JPRB), PARAMETER :: ZDPHISCMN=981._JPRB ! g*dz for shallow
REAL (KIND=JPRB), PARAMETER :: ZTRDTX=4._JPRB ! Max T kick
REAL (KIND=JPRB), PARAMETER :: ZTRBDP=300.E2_JPRB ! dp above LCL

! 1D INTEGER LOCALS
INTEGER(KIND=JPIM) :: JHIBA(KLON), JLCL(KLON), JLCLT(KLON),&
        & JLFC(KLON), JLOBA(KLON), JM(KLON), JNOCUT(KLON),&
        & JUSLSC(KLON),JLCLSC(KLON),JLFCSC(KLON)

! INTEGER SCALARS
INTEGER(KIND=JPIM) :: IASC, IBAS, JLON,JLEV,JL0,JL1, IOK, INO
INTEGER(KIND=JPIM) ::  ILEV, ILEV1, JLFCMIN, JUSLMIN, JUSLMAX
! HARD-CODED TEST PARAMETER
REAL(KIND=JPRB) :: ZTRGIBMN
INTEGER(KIND=JPIM) ::  JKEEP9
!     LOCAL LOGICAL SCALARS
LOGICAL :: LLDO, LLFC, LLTRIG, LLHPHICU, LLUSLMIN

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!---------------------------------------------------------------------------
#include "fcttrm.func.h"
#include "fctdoi.func.h"
!---------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACTRGRC',0,ZHOOK_HANDLE)
ASSOCIATE(GCVALBU=>YDPHY0%GCVALBU, GTRBRC=>YDPHY0%GTRBRC,&
 & GTRGDPMIX=>YDPHY0%GTRGDPMIX, GTRGDPHIMN=>YDPHY0%GTRGDPHIMN, &
 & GTRGPUSLMN=>YDPHY0%GTRGPUSLMN, GTRGAIN=>YDPHY0%GTRGAIN, &
 & GTRTHRS=>YDPHY0%GTRTHRS, GTRTHCK=>YDPHY0%GTRTHCK,&
 & RDTFAC=>YDPHY0%RDTFAC, NJKT1=>YDECUMF%NJKT1)
!------------------------------------------------------------------------
! CONSTANTS AND DERIVED PARAMETERS
!------------------------------------------------------------------------
JKEEP9=1 ! FORCE TRIGGERING PREVIOUSLY ACTIVE
ZTRGIBMN=1._JPRB ! MIN VALUE OF Integrated buoyancy
LLHPHICU=GTRGDPHIMN> 981._JPRB ! at least 100m or ignore.
ZIALBU=1._JPRB/GCVALBU
!---------------------------------------------------------------------------
! INITIALIZE OUTPUT VALUES
KNND(KIDIA:KFDIA)=0   ! triggering indicator
JNOCUT(KIDIA:KFDIA)=0 ! used for earlier triggered ascents (knasc on input)
KUSLT(KIDIA:KFDIA)=1 ! final triggering USL level
JLCLT(KIDIA:KFDIA)=1 ! final LCL level
KLFCT(KIDIA:KFDIA)=1 ! final LFC level
! Initialize also for final diagnostics:
! TEST IF THERE CAN BE TRIGGERING IN THE CURRENT SET:
IF (MAXVAL(PIB(KIDIA:KFDIA)) > ZTRGIBMN) THEN
!---------------------------------------------------------------------------
! DERIVED PARAMETERS
ZCPVMD=RCPV-RCPD
ZRDSRV=RD/RV
ZDPMIX2=GTRGDPMIX*0.5_JPRB
ZDPMX2I=1._JPRB/ZDPMIX2
ZKAPPAI=RCPD/RD
ZRDSRG=RD/RG

JUSLMIN=NJKT1 !=std 300hPa (KLEV/3): highest lev to search top of USL: 
              ! anyway LOOPSLUP will adapt it directly, using GTRGPUSLMN
JLFCMIN=NJKT1 !=std 300hPa (KLEV/4) ! Highest level to search for LFC 

! ZTHRLCL: THRESHOLD FOR RES COND JUST ABOVE LCL
ZTHRLCL=GTRTHRS*GTRBRC
!---------------------------------------------------------------------------
! 0) Auxiliary arrays
!    ZTHESE: Saturated Equivalent potential temperature
!    THIS FORMULA BETTER NOT USED BELOW 200hPA : MOREOVER AT VERY 
!    LOW PRESSURE WE OFTEN GET ZR>1 => PREVENTED BY FORBIDDING QSAT>0.5
DO JLEV=KTDIA,KLEV
   DO JLON=KIDIA,KFDIA
      ZZ=MIN(PQSAT(JLON,JLEV),ZQSMAX)
      ZR=ZZ/(1._JPRB-ZZ)
      ZEXP1=0.2854_JPRB*(1._JPRB-0.28_JPRB*ZR)
      ZEXP2=(3376._JPRB/PT(JLON,JLEV)-2.54_JPRB)*ZR*(1._JPRB+0.81_JPRB*ZR)
      ZTHESE(JLON,JLEV)= PT(JLON,JLEV)*(RATM/PAPRSF(JLON,JLEV))**ZEXP1&
                       & * EXP(ZEXP2)
   ENDDO
ENDDO
!---------------------------------------------------------------------------
! 1) Compute properties of parcels mixed over GTRGDPMIX, 
!    starting from highest usl top
! JUSLMIN: highest level to consider for usl top.
! should be related to GTRGpuslmin, dp from surface or given pressure level
ZSUM(KIDIA:KFDIA)=PAPRS(KIDIA:KFDIA,KLEV)-GTRGPUSLMN
JLOBA(:)=JUSLMIN
JHIBA(:)=KLEV
ZSQ(:)=0._JPRB
ZSQC(:)=0._JPRB
ZST(:)=0._JPRB
ZSDP(:)=0._JPRB
ZQM(:,:)=0._JPRB
ZTM(:,:)=0._JPRB
ZQCM(:,:)=0._JPRB
DO JLEV=KLEV,KTDIA,-1
  DO JLON=KIDIA,KFDIA
    ZZ=RATM/PAPRSF(JLON,JLEV)
    ZTHVE(JLON,JLEV)=PTVE(JLON,JLEV)*(ZZ**RKAPPA)
  ENDDO
ENDDO

! 1a) Accumulate upwards half updraught-source-layers from klev
! Initialize everything to values in KLEV layer:
DO JLON=KIDIA,KFDIA
 ZDP=PDELP(JLON,KLEV)*0.5_JPRB
 ZSDP(JLON)=ZDP
 ZST(JLON)=PCP(JLON,KLEV)*PT(JLON,KLEV)*ZDP
 ZSQ(JLON)=PQ(JLON,KLEV)*ZDP
 ZSQC(JLON)=PQC(JLON,KLEV)*ZDP
 JM(JLON)=KLEV
ENDDO
! upwards loop
LOOPSLUP: DO JLEV=KLEV-1,JUSLMIN,-1
         LLDO=.FALSE.
         LLUSLMIN=.TRUE. ! will remain .T. only when all ZB=0
         DO JLON=KIDIA,KFDIA
! Do not search usl level higher than pressure level in ZSUM
! ZB is 0 when current level is too high
            ZB=MAX(0._JPRB,SIGN(1._JPRB,PAPRSF(JLON,JM(JLON))-ZSUM(JLON)))
            LLUSLMIN=LLUSLMIN.AND.(ZB<0.5_JPRB)
! Weight of coming layer JLEV
! A negative zdp should induce an abort: inadequate discretisation
! but with the DO WHILE , this should never occur
            ZDP=MAX(0._JPRB,MIN(ZDPMIX2-ZSDP(JLON),PDELP(JLON,JLEV)))*ZB
! Accumulate layer JLEV
            ZST(JLON)=ZST(JLON)+ZDP*PCP(JLON,JLEV)*PT(JLON,JLEV)
            ZSQ(JLON)=ZSQ(JLON)+ZDP*PQ(JLON,JLEV)
            ZSQC(JLON)=ZSQC(JLON)+ZDP*PQC(JLON,JLEV)
            ZSDP(JLON)=ZSDP(JLON)+ZDP
! ZOK passes to 1 if accumulation over ZDPMIX2 complete
! with fine vertical discretization, zeps1 should be small for accuracy
            ZOK=MAX(0._JPRB,SIGN(1._JPRB,ZSDP(JLON)-ZDPMIX2+ZEPS1))
! MEMORIES:
            ZS14(JLON)=ZDP
            ZSUM1(JLON)=ZOK
            LLDO=LLDO.OR.(ZOK>0.5_JPRB)
         ENDDO
! Update the averages at level JM when ZOK is 1
! Pass to next usl when ZOK is 1
! Iteration required as long as (one of the ZOK)=1, i.e. LLDO
         DO WHILE (LLDO)
          LLDO=.FALSE.
          DO JLON=KIDIA,KFDIA
            ZDP=ZS14(JLON)
            ZOK=ZSUM1(JLON)
            ZZ=ZOK*ZDPMX2I
            ZQM(JLON,JM(JLON))=ZSQ(JLON)*ZZ
            ZCPM=RCPD+ZCPVMD*ZQM(JLON,JM(JLON))
            ZTM(JLON,JM(JLON))=ZST(JLON)*ZZ/ZCPM
            ZQCM(JLON,JM(JLON))=ZSQC(JLON)*ZZ
            ILEV=JM(JLON)
            ILEV1=ILEV-1
            ZDP0=ZOK*PDELP(JLON,ILEV)*0.5_JPRB
            ZDP1=ZOK*PDELP(JLON,ILEV1)*0.5_JPRB
            ZDP2=ZOK*(PDELP(JLON,JLEV)-ZDP) ! remainder of current JLEV
            ZZ=ZDP1+ZDP0
            ZDP2=MAX(0._JPRB,MIN(ZDP2, ZZ))
            ZDP=ZDP+ZDP2
            ZSDP(JLON)=ZSDP(JLON)-ZZ+ZDP2
            ZSQ(JLON)=ZSQ(JLON)-ZDP0*PQ(JLON,ILEV)-ZDP1*PQ(JLON,ILEV1)&
                              &+ZDP2*PQ(JLON,JLEV)
            ZST(JLON)=ZST(JLON)-ZDP0*PCP(JLON,ILEV)*PT(JLON,ILEV)&
                              &-ZDP1*PCP(JLON,ILEV1)*PT(JLON,ILEV1)&
                              &+ZDP2*PCP(JLON,JLEV)*PT(JLON,JLEV)
            ZSQC(JLON)=ZSQC(JLON)-ZDP0*PQC(JLON,ILEV)-ZDP1*PQC(JLON,ILEV1)&
                              &+ZDP2*PQC(JLON,JLEV)
! Security
             ZSQC(JLON)=MAX(ZSQC(JLON),0._JPRB)
! Update the averages at level JM when ZOK is 1
            JM(JLON)=JM(JLON)-NINT(ZOK)
! Recompute ZDP2: iterate as long as >0
            ZDP2=PDELP(JLON,JLEV)-ZDP
            ZOK=ZOK*MAX(0._JPRB,SIGN(1._JPRB,ZDP2-ZEPS1))
            ZS14(JLON)=ZDP
            ZSUM1(JLON)=ZOK
            LLDO=LLDO.OR.(ZOK>0.5_JPRB)
          ENDDO ! JLON
         ENDDO ! WHILE LLDO
         IF (LLUSLMIN) EXIT LOOPSLUP
ENDDO LOOPSLUP ! JLEV
DO JLON=KIDIA,KFDIA
! Highest USL mid-level to consider
        JM(JLON)=JM(JLON)+1
! Memorize it in JHIBA
        JHIBA(JLON)=JM(JLON)
! Initialize the further accumulations
        ZDP=PDELP(JLON,JM(JLON))*0.5_JPRB
        ZSQ(JLON)=PQ(JLON,JM(JLON))*ZDP
        ZST(JLON)=PCP(JLON,JM(JLON))*PT(JLON,JM(JLON))*ZDP
        ZSQC(JLON)=PQC(JLON,JM(JLON))*ZDP
        ZSDP(JLON)=ZDP
ENDDO
JUSLMIN=MAX(JUSLMIN,MINVAL(JHIBA(KIDIA:KFDIA)))
! 1b) Accumulate downwards the lower half of each usl
LOOPSLDN: DO JLEV=JUSLMIN,KLEV
         LLDO=.FALSE.
         DO JLON=KIDIA,KFDIA
! ZZ is 0 as long as JLEV is above OR AT current JM 
            ZZ=MAX(0._JPRB,SIGN(1._JPRB,JLEV-JM(JLON)-0.5_JPRB))
! Weight of coming layer
            ZDP=MAX(0._JPRB,MIN(ZDPMIX2-ZSDP(JLON),PDELP(JLON,JLEV)))*ZZ
! Accumulate layer JLEV
            ZST(JLON)=ZST(JLON)+ZDP*PCP(JLON,JLEV)*PT(JLON,JLEV)
            ZSQ(JLON)=ZSQ(JLON)+ZDP*PQ(JLON,JLEV)
            ZSQC(JLON)=ZSQC(JLON)+ZDP*PQC(JLON,JLEV)
            ZSDP(JLON)=ZSDP(JLON)+ZDP
! ZOK passes to 1 if accumulation over ZDPMIX2 complete
            ZOK=MAX(0._JPRB,SIGN(1._JPRB,ZSDP(JLON)-ZDPMIX2+ZEPS1))
! Memories
            ZS14(JLON)=ZDP
            ZSUM1(JLON)=ZOK
            LLDO=LLDO.OR.(ZOK>0.5_JPRB)
         ENDDO
         DO WHILE (LLDO)
             LLDO=.FALSE.
             DO JLON=KIDIA,KFDIA
                 ZOK=ZSUM1(JLON)
                 ZDP=ZS14(JLON)
! Update the averages when ZOK is 1
                 ZZ=ZOK*ZDPMX2I*0.5_JPRB
                 ZZ2=1._JPRB-0.5_JPRB*ZOK
                 ZQM(JLON,JM(JLON))=ZQM(JLON,JM(JLON))*ZZ2+ZSQ(JLON)*ZZ
                 ZCPM=RCPD+ZCPVMD*ZQM(JLON,JM(JLON))
                 ZTM(JLON,JM(JLON))=ZTM(JLON,JM(JLON))*ZZ2+ZST(JLON)*ZZ/ZCPM
                 ZQCM(JLON,JM(JLON))=ZQCM(JLON,JM(JLON))*ZZ2+ZSQC(JLON)*ZZ
! Pass to next usl when ZOK is 1
                 ILEV=JM(JLON)
                 ILEV1=ILEV+1
                 ZDP0=ZOK*PDELP(JLON,ILEV)*0.5_JPRB
                 ZDP1=ZOK*PDELP(JLON,ILEV1)*0.5_JPRB
                 ZDP2=ZOK*(PDELP(JLON,JLEV)-ZDP) ! remainder of current JLEV
                 ZZ=ZDP1+ZDP0
                 ZDP2=MAX(0._JPRB,MIN(ZDP2, ZZ))
                 ZDP=ZDP+ZDP2
                 ZSDP(JLON)=ZSDP(JLON)-ZZ+ZDP2
                 ZSQ(JLON)=ZSQ(JLON)-ZDP0*PQ(JLON,ILEV)-ZDP1*PQ(JLON,ILEV1)&
                                   &+ZDP2*PQ(JLON,JLEV)
                 ZST(JLON)=ZST(JLON)-ZDP0*PCP(JLON,ILEV)*PT(JLON,ILEV)&
                                   &-ZDP1*PCP(JLON,ILEV1)*PT(JLON,ILEV1)&
                                   &+ZDP2*PCP(JLON,JLEV)*PT(JLON,JLEV)
                 ZSQC(JLON)=ZSQC(JLON)-ZDP0*PQC(JLON,ILEV)-ZDP1*&
                                   &PQC(JLON,ILEV1)+ZDP2*PQC(JLON,JLEV)
! Security
             ZSQC(JLON)=MAX(ZSQC(JLON),0._JPRB)
! JM should not be able to reach more than KLEV
                 JM(JLON)=JM(JLON)+NINT(ZOK)
! Recompute ZDP2: iterate as long as >0
                 ZDP2=PDELP(JLON,JLEV)-ZDP
                 ZOK=ZOK*MAX(0._JPRB,SIGN(1._JPRB,ZDP2-ZEPS1))
                 ZS14(JLON)=ZDP
                 ZSUM1(JLON)=ZOK
                 LLDO=LLDO.OR.(ZOK>0.5_JPRB)
             ENDDO ! JLON
         ENDDO ! WHILE LLDO
ENDDO LOOPSLDN ! JLEV
! Memorize lowest base in JLOBA
DO JLON=KIDIA,KFDIA
   JLOBA(JLON)=JM(JLON)-1
ENDDO
JUSLMAX=MAXVAL(JLOBA(KIDIA:KFDIA))

! Now ZS*M(JM) has the mean properties of parcel around level JM
! Mixing GTRGDPMIX/2 below and above it.
! Hence no more need to consider usltop but we can concentrate 
! the parcel at this level.
! Valid JM values at each vertical are between jhiba and jloba.
! To ease block calculation, copy values at jloba below,
! and values at jhiba above, to fill juslmax to juslmin.
DO JLEV=JUSLMIN,JUSLMAX
   DO JLON=KIDIA,KFDIA
      ILEV=MAX(JHIBA(JLON),MIN(JLOBA(JLON),JLEV))
      ZQM(JLON,JLEV)=ZQM(JLON,ILEV)
      ZQCM(JLON,JLEV)=ZQCM(JLON,ILEV)
      ZTM(JLON,JLEV)=ZTM(JLON,ILEV)
   ENDDO
ENDDO
!----------------------------------------
! 2) Compute LCL temperature and pressure for each usl particle
ZA12(:,:)=0._JPRB
ZA13(:,:)=1._JPRB ! Initialize p_LCL to 1 Pa to prevent division by zero
DO JLEV=JUSLMAX,JUSLMIN,-1
   DO JLON=KIDIA,KFDIA
! Re-evaporate any condensate:
      ZDELTA=FONICE(ZTM(JLON,JLEV),YDPHY0%RDTFAC)
      ZLM=FOLH(ZTM(JLON,JLEV),ZDELTA)
! FORBID COMPLETELY DRY USL (CRASH IN LOG)
      ZQSL=MAX(ZEPSQV,ZQM(JLON,JLEV)+ZQCM(JLON,JLEV))
      ZCPM=RCPD+ZCPVMD*ZQSL
      ZTSL=ZTM(JLON,JLEV)-ZLM/ZCPM*ZQCM(JLON,JLEV)
      ZESL=ZQSL*PAPRSF(JLON,JLEV)/(ZQSL+ZRDSRV*(1._JPRB-ZQSL))
      ZEW=FOEW(ZTSL,FONICE(ZTSL,YDPHY0%RDTFAC))
! Check saturation: in this case take blue point of this level and safi
! (ever seen a cumulonimbus starting at the surface ?)
! ZBLUE is 1 if saturation.
      ZBLUE=MAX(0._JPRB, SIGN(1._JPRB, ZESL-ZEW))
      ZQSL=ZQSL*(1._JPRB-ZBLUE) + PQW(JLON,JLEV)*ZBLUE
      ZTSL=ZTSL*(1._JPRB-ZBLUE) + PTW(JLON,JLEV)*ZBLUE
! IF zblue=1, zesl no longer used.
! If not saturated, compute LCL
! STORE TLCL in ZA12, PLCL in ZA13
      ZA12(JLON,JLEV)=(1._JPRB-ZBLUE)*(2840._JPRB/(3.5_JPRB*LOG(ZTSL)&
                   & -LOG(ZESL*0.01_JPRB)-4.805_JPRB) + 55._JPRB)&
                   & +ZBLUE*ZTSL
      ZA13(JLON,JLEV)=PAPRSF(JLON,JLEV)*(ZA12(JLON,JLEV)/ZTSL)**ZKAPPAI
      ZTM(JLON,JLEV)=ZTSL
      ZQM(JLON,JLEV)=ZQSL
   ENDDO
ENDDO

!-----------------------------------------------------
! 3. ACCUMULATE DRY AND MOIST CIN STARTING AT EACH USL
!-----------------------------------------------------
ZTVKICK=0._JPRB ! TO REDO FOR EACH USL ? NOT IF TRIGGERED
! INITIALIZE SHALLOW CLOUD VALUES
ZKICKSC=0._JPRB
ZDPHISC=0._JPRB
JUSLSC=1
JLFCSC=1
JLCLSC=1
! PRECOMPUTE REFERENCE ACCUMULATIONS
ZSFCS(KIDIA:KFDIA)=PIFHCS(KIDIA:KFDIA,KLEV)
ZSDP(KIDIA:KFDIA)=PDELP(KIDIA:KFDIA,KLEV)
DO JLEV=KLEV-1,KTDIA,-1
  DO JLON=KIDIA,KFDIA
! Limit rc accumulation to a dpx=gtrthck above surface 
     ZZ=MAX(0._JPRB,SIGN(1._JPRB,&
              & PAPRSF(JLON,JLEV)-PAPRS(JLON,KLEV)+GTRTHCK))
     ZSFCS(JLON)=ZSFCS(JLON) + PIFHCS(JLON,JLEV)*ZZ
     ZSDP(JLON)=ZSDP(JLON)+PDELP(JLON,JLEV)*ZZ
  ENDDO
ENDDO
ZSFCS(KIDIA:KFDIA)=ZSFCS(KIDIA:KFDIA)/ZSDP(KIDIA:KFDIA)
! Compute Thresholds
DO JLON=KIDIA,KFDIA
! Reference kick, to reduce afterwards following WLCL
   ZKICKRC(JLON)=GTRGAIN*(MAX(0._JPRB,ZSFCS(JLON))-GTRTHRS)
ENDDO

LOOPUSL:DO JL0=JUSLMAX, JUSLMIN,-1     ! vertical loop for each usl
!*INDICATORS:
! ZCONTINUE(JLON): AS SOON AS ZERO, THE CURRENT USL IS CANCELLED
!                  OR TRIGGERED: CAN EXIT INNER LOOP LOOPLFC (ON JL1)
! KNND(JLON): IF 1 IT MEANS CURRENT COLUMN HAS BEEN TRIGGERED 
!                  CAN EXIT OUTER LOOP LOOPUSL (on JL0)
! JNOCUT(KLON)=ZINAC9: IF 1 MEANS CURRENT COLUMN IS PROTECTED BECAUSE
!            WE PASSED THE LCL AND UPDRAUGHT WAS ALREADY ACTIVE
!            => FORCES ALL FURTHER CONTROLS TO OK.
! LLTRIG: T WHEN ALL COLUMNS HAVE NEEN TRIGGERED 
! ZOK (CURRENT LEVEL): REMAINS 1 WHEN TRIGGERING HAS OCCURRED
!
!*MEMORIES:
! JLCL(KLON)
! JLFC(KLON)
! ZDTVLCL(KLON): BUOYANCY DEFICIT AT LCL (if <0)
! ZDTVFC(KLON): RESOLVED CONDENSATION BUOYANCY KICK AT LCL
! ZTCIN(KLON): REQUIRED BUOYANCY TO PASS MOIST CIN FROM LCL TO LFC
 LLFC=.FALSE. ! .T. when all columns from JL0 reach LFC OR stopped 
 ZS6=0._JPRB
 ZSDP=0._JPRB
 ZST=0._JPRB
 ZSQC=0._JPRB
 ZS14=0._JPRB
 ZSQ=0._JPRB
 ZDTVFC=0._JPRB
 ZDTVLCL=0._JPRB
 ZTCIN=0._JPRB
 ZSDP=0._JPRB ! reset to zero, case LLTKICK only
 ZSUM=0._JPRB  ! passes to 1 when reached the LFC
 DO JLON=KIDIA,KFDIA
    ILEV=MIN(JL0,JLOBA(JLON)) ! effective starting level at jlon
    JM(JLON)=ILEV
    JLCL(JLON)=0 ! LCL level for usl=jl0
    JLFC(JLON)=0 ! LFC level for usl=jl0
    ZTV=(1._JPRB+RETV*ZQM(JLON,ILEV))*ZTM(JLON,ILEV)
    ZTHVSL(JLON)=ZTV*(RATM/PAPRSF(JLON,ILEV))**RKAPPA
    ZR=ZQM(JLON,ILEV)/(1._JPRB-ZQM(JLON,ILEV))
    ZEXP1=0.2854_JPRB*(1._JPRB-0.28*ZR)
    ZEXP2=(3376._JPRB/ZTM(JLON,ILEV)-2.54_JPRB)*ZR&
     & *(1._JPRB+0.81_JPRB*ZR)
    ZTHEQSL(JLON)=ZTM(JLON,ILEV)*(RATM/PAPRSF(JLON,ILEV))**ZEXP1&
     & * EXP(ZEXP2)
    ZCONTINUE(JLON)=MAX(0._JPRB, SIGN(1._JPRB, PIB(JLON)-ZTRGIBMN))
 ENDDO
 LLTRIG=.TRUE. ! will remain .T. when all colums have triggered 
 LOOPLFC:DO JL1=JL0, JLFCMIN,-1       ! vertical accumulation for JL0
   DO JLON=KIDIA,KFDIA               ! Vectorizable Horizontal loop
! ZOK0=0 if JL0  is below JM or JL1 at current JM, 1 if above JM
! IF this column already triggered, maintain ZOK0=0
       ZOK0=MAX(0._JPRB,SIGN(1._JPRB,&
        & REAL(JM(JLON)-JL0,JPRB)+0.5_JPRB))*(1-KNND(JLON))
!        ZOK0=MAX(0._JPRB,SIGN(1._JPRB,&
!            & REAL((JM(JLON)-JL1)*(1-KNND(JLON)),JPRB)-0.5_JPRB))
       ZPLCL=ZA13(JLON,JL0)
!               ===============================
!                   DRY ASCENT BELOW LCL
!               ===============================
! ZOK1 is 1 while layer JL1 is in ascent and entirely below LCL
       ZOK1=ZOK0*MAX(0._JPRB,SIGN(1._JPRB,PAPRS(JLON,JL1-1)-ZPLCL))
       ZB=ZTHVSL(JLON)/ZTHVE(JLON,JL1)-1._JPRB
       ZDPHI=(PAPHI(JLON,JL1-1)-PAPHI(JLON,JL1))
! ZOK2=1 if part or all of JL1 layer is above LCL: 
!        then interpolate and continue with Moist CIN
       ZOK2=ZOK0*(1._JPRB-ZOK1)
       JLCL(JLON)=NINT(ZOK2)*MAX(JLCL(JLON),JL1)
! ZOK3=1 If JL1 layer totally above LCL (then reset ZOK2=0)
       ZOK3=ZOK0*MAX(0._JPRB,SIGN(1._JPRB,ZPLCL-PAPRS(JLON,JL1)))
       ZOK2=ZOK2*(1._JPRB-ZOK3)
 
 
! If LCL level (zok2=1) was already active, ignore the dry cin criterion
! IF active from previous dt: knasc=1, knnd still 0, zok2=1
! then set kuslt directly as well as the LCL 
! IBAS is 1 if jl1 is the LCL and is previously active and we have not
!       triggered earlier. Then JNOCUT passes to 1 and remains 1
! OTHERWISE, IBAS REMAINS ZERO AND JNOCUT TOO
!---------------------------------------------------------------------
! KNASC will contain on entry INAC9<<<
       IBAS=(1-KNND(JLON))*KNASC(JLON,JL1)*JKEEP9*NINT(ZOK2)
       KUSLT(JLON)=MAX(KUSLT(JLON),JL0*IBAS)
       JLCLT(JLON)=MAX(JLCLT(JLON),JLCL(JLON)*IBAS)
       JNOCUT(JLON)=MAX(JNOCUT(JLON),IBAS)
! ABOVE LCL (ZOK3=1): if zcontinue still 1 or knasc=1 (zinac9), 
!       it implies we have passed dry cin and LCL buoyancy criteria
       ZINAC9=REAL(JNOCUT(JLON),JPRB)
       ZOK=MAX(ZINAC9,ZCONTINUE(JLON))*ZOK3
!------------------------------------------------------
 
!                     ==============
!                       LCL KICKING 
!                     ==============
!                JL1 IS THE LCL WHEN ZOK2==1
! COMPUTE MISSING BUOYANCY OF LCL (KEEP IN MEMORY)
       ZLIQ=PQC(JLON,JL1)*ZOK2 ! just to compensate mean qc at lcl
       ZDTVLCL(JLON)=ZOK2*(ZA12(JLON,JL0)*&
        & (1._JPRB-ZLIQ+RETV*(ZQM(JLON,JL0)-ZLIQ))&
        & -PTVE(JLON,JL1))+(1._JPRB-ZOK2)*ZDTVLCL(JLON)
! KICK FROM INTEGRATED RESOLVED CONDENSATION FROM Surface to psurf-dpx
       ZDTVFC(JLON)=(1._JPRB-ZOK2)*ZDTVFC(JLON)+ZOK2*ZKICKRC(JLON)
!               ====================================
!                MOIST ASCENT above LCL (set ZOK3=1)
!               ====================================
 
       ZOK3=MAX(ZOK2,ZOK3)
       ZB=(ZTHEQSL(JLON)/ZTHESE(JLON,JL1)-1._JPRB)*ZDPHI
! ZZ is 1 above the LFC 
       ZZ=SIGN(1._JPRB,ZB)*ZOK3
       ZSUM(JLON)= MAX(ZSUM(JLON), ZZ)
! ACCUMULATE moist cin (<0) in zs6 and dphi in ZSDP 
! and compute mean T and condensate
! ZSUM passes to 1 once we have passed the CIN 
! KEEP ZINCIN=1 IF LFC=LCL
       ZINCIN=MAX(ZOK2,ZOK3*(1._JPRB-ZSUM(JLON)))
       ZS6(JLON)=ZS6(JLON)+ZINCIN*MIN(0._JPRB,ZB)
       ZS14(JLON)=ZS14(JLON)+ZINCIN*ZDPHI
       ZST(JLON)=ZOK2*PT(JLON,JL1) +(1._JPRB-ZOK2)*ZST(JLON)
! TAKE INCREMENT OF RES COND 
       ZSQC(JLON)=ZSQC(JLON)+ZOK3*PIFHCS(JLON,JL1)
       ZSDP(JLON)=ZSDP(JLON)+ZOK3*PDELP(JLON,JL1)
! ZB>0 if we have reached the LFC: then finished for this usl
! Update ZSUM and JLFC
       JLFC(JLON)=MAX(JLFC(JLON),NINT(ZZ)*JL1)
! STORE ZTCIN UP TO  THE LFC: OLD ZSUM WAS ZERO,PASSES TO 1
       ZTCIN(JLON)=ZINCIN*MAX(0._JPRB, -ZS6(JLON)/MAX(1._JPRB,ZS14(JLON)))&
          &      *ZST(JLON)*ZIALBU&
          &      +(1._JPRB-ZINCIN)*ZTCIN(JLON)
! UP TO LFC COMPARE ZDTVFC +ZDTVLCL WITH ZTCIN 
       ZOK4=MAX(ZINAC9,SIGN(1._JPRB,ZDTVLCL(JLON)+ZDTVFC(JLON)-ZTCIN(JLON)))
! if ZOK4=0 between LCL and LFC (ZINCIN=1) => cancel
! UNLESS SUFFICIENT UPWARDS RESOLVED VELOCITY
       ZZ=MAX(0._JPRB,1._JPRB-ZINCIN*(1._JPRB-ZOK4))
       ZCONTINUE(JLON)=ZCONTINUE(JLON)*ZZ
! As long as ZSUM=0 we must continue and triggering not confirmed
       ZOK=ZOK*ZZ*ZSUM(JLON)
       IF (LLHPHICU) THEN ! REQUIRE MIN ELEVATION ABOVE LCL
! elevation above LCL:
            ZHPHI=PAPHIF(JLON,JL1)-PAPHIF(JLON,MAX(1,JLCL(JLON)))
! ZZ is 0 between (LFC: zsum=1) and GTRGDPHIMN (1 outside): if there ZB<0, 
! and zsum=1, then cancel unless already active 
           ZZ=1._JPRB-MAX(0._JPRB,SIGN(1._JPRB,GTRGDPHIMN-ZHPHI))*ZSUM(JLON)
           ZZ2=MAX(0._JPRB, SIGN(1._JPRB,ZB))
! ZCONTINUE goes to zero as soon as we must cancel current usl:
!       both ZZ2 and ZZ are 0
! AT TIME OF CANCELLING, 
! CHECK SHALLOW HEIGHT LARGER THAN EARLIER BEFORE CANCELING
            ZDEEPER=MAX(0._JPRB,ZCONTINUE(JLON)-MAX(ZZ,ZZ2,ZINAC9))&
                   & *MAX(0._JPRB,SIGN(1._JPRB,ZHPHI-ZDPHISC(JLON))) 
            ZNO=1._JPRB-ZDEEPER
!-------------------------------------------
! UPDATE THE DEEPEST SHALLOW CLOUD VALUES
! Shallow cloud depth
            ZDPHISC(JLON)=ZDPHISC(JLON)*ZNO + ZHPHI*ZDEEPER
! Shallow cloud kick
            ZKICKSC(JLON)=ZKICKSC(JLON)*ZNO+ZDEEPER*&
            & MIN(ZTCIN(JLON),ZDTVLCL(JLON)+MAX(0._JPRB,ZDTVFC(JLON)))
            IOK=NINT(ZDEEPER)
            INO=1-IOK
! Shallow cloud usl
            JUSLSC(JLON)=JUSLSC(JLON)*INO+JL0*IOK
! Shallow cloud lcl
            JLCLSC(JLON)=JLCLSC(JLON)*INO+JLCL(JLON)*IOK
! shallow cloud lfc
            JLFCSC(JLON)=JLFCSC(JLON)*INO+JLFC(JLON)*IOK
!-------------------------------------------
! NOW CANCEL TOO SHALLOW CLOUD
            ZCONTINUE(JLON)=ZCONTINUE(JLON)* MAX(ZZ,ZZ2,ZINAC9)
! Either KNND=1 or ZZ=1=ZZ2=ZSUM(JLON) - criterion does not play if jnocut
            ZOK=ZOK*MAX(ZINAC9,REAL(KNND(JLON),JPRB),ZZ*ZZ2*ZSUM(JLON))
       ENDIF !LLHPHICU
! Test that there is resolved condensation below a min thickness above LCL,
! otherwise this was a fake triggering
! This is directly given by the contents of PQC 
! ZZ goes to 1 when we are enough above LCL
       ZZ=MAX(0._JPRB,SIGN(1._JPRB,ZSDP(JLON)-ZTRBDP))
! ZZ2 is 1 as long as ZSQC too small : ZTHRLCL 
       ZZ2=MAX(0._JPRB,SIGN(1._JPRB,ZTHRLCL*ZSDP(JLON)-ZSQC(JLON)))
! Memorize humidity criterion for the case zinac9=1
! Triggered if zz=1 and zz2=0: zsq=1
       ZSQ(JLON)=ZZ*(1._JPRB-ZZ2)
       ZOK5=MAX(ZINAC9,ZSQ(JLON))*ZZ
! Cancel if both zz and zz2 are 1
       ZCONTINUE(JLON)=ZCONTINUE(JLON)*(1._JPRB-ZZ*ZZ2)
       ZOK=ZOK*ZOK5
! Do not continue if triggered (then ZOK has remained 1 !)
       ZCONTINUE(JLON)=ZCONTINUE(JLON)*(1._JPRB-ZOK)
! Condition to pass to next usl
       LLFC=LLFC .AND. (ZCONTINUE(JLON)<0.5_JPRB)
       IOK=NINT(ZOK)
! IF ZOK=1 and knnd still 0 we have triggered NOW: the base is JLCL(JLON)
! if knnd=1 and zok2=1, we are at the lcl and triggered by history
       IBAS=IOK*(1-KNND(JLON))
! When IBAS=1, it means current usl is the triggering one
       KUSLT(JLON)=IBAS*JL0+(1-IBAS)*KUSLT(JLON)
       JLCLT(JLON)=IBAS*JLCL(JLON)+(1-IBAS)*JLCLT(JLON)
       KLFCT(JLON)=IBAS*JLFC(JLON)+(1-IBAS)*KLFCT(JLON)
       KNND(JLON)=MAX(KNND(JLON),IOK)
! IF protected (jnocut=1) but negative kick: what ? limit it or not ?
! RESET THE POSITIVE KICKS IF MOISTURE CRITERION NOT FULFILLED (case jnocut)
! ZSQ is 1 or zero
! REDUCE THE KICKS (+ and -) WHEN W0 AT LCL LARGE
! Limit the kick to (a multiple of) ZTCIN The one necessary to pass the CIN
       ZTV=MIN(ZDTVLCL(JLON)+ZDTVFC(JLON),ZTCIN(JLON))
       ZZ=SIGN(1._JPRB, ZTV)
       ZTV=MIN(ZTRDTX,ABS(ZTV))*ZZ
       ZTVKICK(JLON)=IBAS*MIN(10._JPRB*ZSQ(JLON), ZTV)&
                & +(1-IBAS)*ZTVKICK(JLON)
       LLTRIG=LLTRIG .AND. KNND(JLON)==1
   ENDDO !JLON
   IF (LLFC) EXIT LOOPLFC ! then pass to next usl unless all triggered
 ENDDO LOOPLFC
 IF (LLTRIG) EXIT LOOPUSL
ENDDO LOOPUSL

! CHECK SELECTION OF SHALLOW CLOUD: deeper than 100m and no deep found
DO JLON=KIDIA,KFDIA
   KNNS(JLON)=(1-KNND(JLON))*NINT(MAX(0._JPRB,&
              &         SIGN(1._JPRB, ZDPHISC(JLON)-ZDPHISCMN)))
   ZTVKICK(JLON)=ZTVKICK(JLON)*KNND(JLON)+ZKICKSC(JLON)*KNNS(JLON) 
   JLCLT(JLON)=JLCLT(JLON)*KNND(JLON)+JLCLSC(JLON)*KNNS(JLON)
   KUSLT(JLON)=MAX(1,KUSLT(JLON)*KNND(JLON)+JUSLSC(JLON)*KNNS(JLON))
   KLFCT(JLON)=KLFCT(JLON)*KNND(JLON)+JLFCSC(JLON)*KNNS(JLON)
ENDDO
! Now update KNASC and output arrays
KNASC(KIDIA:KFDIA,KLEV)=0
DO JLEV=KLEV-1,KTDIA+1,-1
   DO JLON=KIDIA,KFDIA
      IASC=MAX(0,SIGN(1,JLCLT(JLON)-JLEV)) ! *KNND(JLON)
      KNASC(JLON,JLEV)=IASC
      ZBAS=REAL(KNASC(JLON,JLEV)*(1-KNASC(JLON,JLEV+1)),JPRB)
      ILEV=KUSLT(JLON)
      PQU(JLON,JLEV)=ZBAS*(ZQM(JLON,ILEV)-PQ(JLON,JLEV))
      PTU(JLON,JLEV)=ZBAS*(ZA12(JLON,ILEV)-PT(JLON,JLEV))
      PSU(JLON,JLEV)=ZBAS*(RCPD*PTU(JLON,JLEV)+ZCPVMD*(&
         & PTU(JLON,JLEV)*(PQU(JLON,JLEV)+PQ(JLON,JLEV))&
         & +PT(JLON,JLEV)*PQU(JLON,JLEV)))
!--------------------------------------------------
! if ZTVKICK=0 it means we did not find the base 
      PTVU(JLON,JLEV)=ZBAS*ZTVKICK(JLON)
   ENDDO
ENDDO
DO JLON=KIDIA,KFDIA
    ILEV=KUSLT(JLON)
! RETURN VALUES 
    PPLCL(JLON)=PAPRSF(JLON,ILEV)
    PTLCL(JLON)=ZA12(JLON,ILEV)
    PTUSL(JLON)=ZTM(JLON,ILEV)
    PQUSL(JLON)=ZQM(JLON,ILEV)
ENDDO
ENDIF ! MAXVAL(PIB) TEST
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACTRGRC',1,ZHOOK_HANDLE)
 END SUBROUTINE ACTRGRC
