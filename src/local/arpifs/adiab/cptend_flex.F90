SUBROUTINE CPTEND_FLEX ( YDLDDH,YDMDDH,YGFL,YDPHY,KPROMA, KSTART, KPROF, KFLEV,PGNORDL,PGNORDM,&
 !  2D input variables
 & PDELP,PRDELP , PCP, & 
 ! state variables
 & PU, PV, PT, PTS, PGFL, & 
 ! set of processes
 & YDPROCSET ,& 
 ! output tendencies
 & PTENDU , PTENDV , PTENDH, PTENDGFL , & 
 ! 2D diagnostic output
 & PFHSCL , PFHSCN , PFHSSL , PFHSSN,&
 & PFHPCL , PFHPCN , PFHPSL , PFHPSN,&
 & PFHP   , PFP    , PFEPFP , PFCMPCQ, PFCMPSN, PFCMPSL, &
 & YDDDH)

!     ------------------------------------------------------------------
!     CALCULATION OF U, V, ENTHALPY AND GFL TENDENCIES FROM
!     FLUXES AND PARTIAL TENDENCIES
!     ------------------------------------------------------------------

!     SEE PHYSICS-DYNAMICS INTERFACE DOCUMENTATION (DEGRAUWE, GELEYN AND BOUYSSEL, 2010)
!                            -----------------

!     IMPORTANT NOTE: THE THERMODYNAMIC VARIABLE IS ENTHALPY
!                     (MORE PRECISELY THE MOIST STATIC ENERGY)
! 
!     Author
!     ------
!       2013-11, D. Degrauwe
!-----------------------------------------------------------------------

USE YOMMDDH  , ONLY : TMDDH
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK
USE YOMCST   , ONLY : RG       ,RV       ,&
 &                    RCPV     ,RETV     ,RCW      ,RCS      ,RLVTT    ,&
 &                    RLSTT    ,RTT      ,RALPW    ,RBETW    ,RGAMW    ,&
 &                    RALPS    ,RBETS    ,RGAMS    ,RALPD    ,RBETD    ,&
 &                    RGAMD    ,RCPD     
USE YOMPHY     , ONLY : TPHY
USE YOMLDDH    , ONLY : TLDDH
USE DDH_MIX    , ONLY : ADD_FIELD_3D, NEW_ADD_FIELD_3D, TYP_DDH
USE INTFLEX_MOD, ONLY : TYPE_INTPROCSET, TYPE_INTPROC
USE YOM_YGFL   , ONLY : TYPE_GFLD
USE YOMLSFORC  , ONLY : LMUSCLFA,NMUSCLFA
USE YOMCT0     , ONLY : LTWOTL
USE YOMLUN     , ONLY : NULOUT

!-----------------------------------------------------------------------


IMPLICIT NONE

! input: dimensions
TYPE(TLDDH)       ,INTENT(IN)    :: YDLDDH
TYPE(TMDDH)       ,INTENT(IN)    :: YDMDDH
TYPE(TPHY)        ,INTENT(IN)    :: YDPHY
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGNORDL (KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGNORDM (KPROMA)
! input: state variables
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELP   (KPROMA,0:KFLEV+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP  (KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCP     (KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU      (KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV      (KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT      (KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS     (KPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL    (KPROMA,KFLEV,YGFL%NUMFLDS)
! input: physics contributions (set of processes)
TYPE(TYPE_INTPROCSET),INTENT(IN) :: YDPROCSET
! output: tendencies of prognostic variables
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDU  (KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDV  (KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDH  (KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDGFL(KPROMA,KFLEV,YGFL%NDIM1)
! output: diagnostics
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHSCL  (KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHSCN  (KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHSSL  (KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHSSN  (KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPCL  (KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPCN  (KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPSL  (KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFHPSN  (KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT), OPTIONAL   :: PFHP    (KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT), OPTIONAL   :: PFP     (KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFEPFP  (KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCMPCQ(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCMPSN(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFCMPSL(KPROMA,0:KFLEV)
TYPE(TYP_DDH), INTENT(INOUT)     :: YDDDH

!-----------------------------------------------------------------------

! auxiliary variables
LOGICAL :: LLDPSFI

REAL(KIND=JPRB) :: ZGSDP (KPROMA,KFLEV)                 ! g/dp
REAL(KIND=JPRB) :: ZTI   (KPROMA,0:KFLEV)               ! temperature at half levels
REAL(KIND=JPRB) :: ZLIFT (KPROMA,0:KFLEV)               ! total absolute precipitation == lifting term in barycentric system
REAL(KIND=JPRB) :: ZQXH  (KPROMA,0:KFLEV,YGFL%NDIM)     ! mass fractions at half levels

REAL(KIND=JPRB) :: ZFLUXU  (KPROMA,0:KFLEV)             ! total U momentum flux
REAL(KIND=JPRB) :: ZFLUXV  (KPROMA,0:KFLEV)             ! total V momentum flux
REAL(KIND=JPRB) :: ZFLUXH  (KPROMA,0:KFLEV)             ! total enthalpy flux
REAL(KIND=JPRB) :: ZFLUXGFL(KPROMA,0:KFLEV,YGFL%NDIM1)  ! total GFL flux

REAL(KIND=JPRB) :: ZCP  (KPROMA,0:KFLEV)                ! kind of total cp (including delta_m influence)
REAL(KIND=JPRB) :: ZCPP (KPROMA,0:KFLEV)                ! contribution of precipitating species to total cp (needed for diagnostics only)
REAL(KIND=JPRB) :: ZCPD (KPROMA,0:KFLEV)                ! contribution of diffusive species to total cp (needed for diagnostics only)
REAL(KIND=JPRB) :: ZJTMP(KPROMA,0:KFLEV)                ! partial enthalpy flux

INTEGER(KIND=JPIM) :: JLEV, JROF
INTEGER(KIND=JPIM) :: JFIELD, JPROCESS, JFIELD2
INTEGER(KIND=JPIM) :: JGFL, JGFLTARGET
REAL(KIND=JPRB) :: ZVAR(KPROMA,1:KFLEV)                ! auxiliary variable
REAL(KIND=JPRB) :: ZVAR2(KPROMA,1:KFLEV)               ! auxiliary variable 2
REAL(KIND=JPRB) :: ZFLX(KPROMA,0:KFLEV)                ! auxiliary flux
REAL(KIND=JPRB) :: ZFLX2(KPROMA,0:KFLEV)               ! auxiliary flux 2
REAL(KIND=JPRB) :: ZUZGEO(KPROMA,1:KFLEV)              ! u w.r.t. geographic north
REAL(KIND=JPRB) :: ZVMGEO(KPROMA,1:KFLEV)              ! v w.r.t. geographic north
REAL(KIND=JPRB) :: ZQ1(KPROMA,KFLEV), ZQ2(KPROMA,KFLEV)

TYPE(TYPE_INTPROC), POINTER :: YLPROC
REAL(KIND=JPRB),    POINTER :: ZVAL(:,:), ZVAL2(:,:)

CHARACTER(LEN=2) :: CLFLEXDIA(YGFL%NUMFLDS)

INTEGER(KIND=JPIM) :: IPGFL(YGFL%NUMFLDS)              ! pointer to t0 or t9, depending on LTWOTL

REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "wrscmr.intfb.h"
#include "fcttrm.func.h"

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CPTEND_FLEX',0,ZHOOK_HANDLE)

ASSOCIATE(NDIM=>YGFL%NDIM, NDIM1=>YGFL%NDIM1, NUMFLDS=>YGFL%NUMFLDS, &
 & YCOMP=>YGFL%YCOMP, YG=>YGFL%YG, YH=>YGFL%YH, YI=>YGFL%YI, YL=>YGFL%YL, &
 & YQ=>YGFL%YQ, YR=>YGFL%YR, YS=>YGFL%YS, &
 & LFLEXDIA=>YDLDDH%LFLEXDIA, LSDDH=>YDLDDH%LSDDH, LDDH_OMP=>YDLDDH%LDDH_OMP, &
 & LPHSPSH=>YDPHY%LPHSPSH, NDPSFI=>YDPHY%NDPSFI)
!-----------------------------------------------------------------------

! 
! 0. Preparations
! 

! 0.1. GFL pointers at previous timestep
DO JGFL=1,NUMFLDS
  IF (YCOMP(JGFL)%LACTIVE) THEN
    IF (LTWOTL) THEN
      IPGFL(JGFL)=YCOMP(JGFL)%MP
    ELSE
      IPGFL(JGFL)=YCOMP(JGFL)%MP9
    ENDIF
  ENDIF
ENDDO

! 0.2. pass variables to ddh
IF (LFLEXDIA) THEN
  ! GFL variables
  !     get names and values for flexible diagnostics; this is quite ugly,
  !     but it's because ddh uses its own names, instead of using gfl names
  
  ! set names to unknown
  CLFLEXDIA(:)='**'
  
  ! set known names
  IF (YQ%LT1) CLFLEXDIA(YQ%MP1)='QV'
  IF (YL%LT1) CLFLEXDIA(YL%MP1)='QL'
  IF (YI%LT1) CLFLEXDIA(YI%MP1)='QI'
  IF (YR%LT1) CLFLEXDIA(YR%MP1)='QR'
  IF (YS%LT1) CLFLEXDIA(YS%MP1)='QS'
  IF (YG%LT1) CLFLEXDIA(YG%MP1)='QG'
  IF (YH%LT1) CLFLEXDIA(YH%MP1)='QH'
  
  ! loop over all GFL variables
  DO JGFL=1,NUMFLDS
    IF (YCOMP(JGFL)%LT1) THEN
      ! if name is unknown, set it to beginning of CNAME attribute
      IF ( CLFLEXDIA(YCOMP(JGFL)%MP1)=='**' ) THEN
        CLFLEXDIA(YCOMP(JGFL)%MP1)=YCOMP(JGFL)%CNAME(1:2)
      ENDIF
      ! pass values to ddh
      ZVAR(KSTART:KPROF,1:KFLEV)= PDELP(KSTART:KPROF,1:KFLEV)*PGFL(KSTART:KPROF,1:KFLEV,IPGFL(JGFL))
      IF (LDDH_OMP) THEN
        CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAR,'V'//CLFLEXDIA(YCOMP(JGFL)%MP1),YDDDH)
      ELSE
        CALL ADD_FIELD_3D(YDLDDH,ZVAR,'V'//CLFLEXDIA(YCOMP(JGFL)%MP1),'V','ARP',.TRUE.,.TRUE.)
      ENDIF
    ENDIF
  ENDDO

  ! Other variables
  ! 
  
  ! enthalpy value
  ZVAR(KSTART:KPROF,1:KFLEV)=&
   & PDELP(KSTART:KPROF,1:KFLEV)*PCP(KSTART:KPROF,1:KFLEV)*PT(KSTART:KPROF,1:KFLEV)
  IF (LDDH_OMP) THEN
    CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAR,'VCT',YDDDH)
  ELSE
    CALL ADD_FIELD_3D(YDLDDH,ZVAR,'VCT','V','ARP',.TRUE.,.TRUE.)
  ENDIF
  
  ! rotate wind to geographic north
  DO JLEV = 1, KFLEV
    DO JROF=KSTART,KPROF
      ZUZGEO(JROF,JLEV) =PGNORDM(JROF)*PU(JROF,JLEV)-PGNORDL(JROF)*PV(JROF,JLEV)
      ZVMGEO(JROF,JLEV) =PGNORDL(JROF)*PU(JROF,JLEV)+PGNORDM(JROF)*PV(JROF,JLEV)
    ENDDO
  ENDDO
  
! pass wind values to diagnostics
  IF (LDDH_OMP) THEN
    ZVAR(KSTART:KPROF,1:KFLEV)=&
     & PDELP(KSTART:KPROF,1:KFLEV)*ZUZGEO(KSTART:KPROF,1:KFLEV)
    CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAR,'VUU',YDDDH)
    ZVAR(KSTART:KPROF,1:KFLEV)=&
     & PDELP(KSTART:KPROF,1:KFLEV)*ZVMGEO(KSTART:KPROF,1:KFLEV)
    CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAR,'VVV',YDDDH)

    ! kinetic energy value
    ZVAR(KSTART:KPROF,1:KFLEV)=&
     & PDELP(KSTART:KPROF,1:KFLEV)*0.5_JPRB*(&
     & PU(KSTART:KPROF,1:KFLEV)**2+PV(KSTART:KPROF,1:KFLEV)**2)
    CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAR,'VKK',YDDDH)
  ELSE
    ZVAR(KSTART:KPROF,1:KFLEV)=&
     & PDELP(KSTART:KPROF,1:KFLEV)*ZUZGEO(KSTART:KPROF,1:KFLEV)
    CALL ADD_FIELD_3D(YDLDDH,ZVAR,'VUU','V','ARP',.TRUE.,.TRUE.)
    ZVAR(KSTART:KPROF,1:KFLEV)=&
     & PDELP(KSTART:KPROF,1:KFLEV)*ZVMGEO(KSTART:KPROF,1:KFLEV)
    CALL ADD_FIELD_3D(YDLDDH,ZVAR,'VVV','V','ARP',.TRUE.,.TRUE.)

    ! kinetic energy value
    ZVAR(KSTART:KPROF,1:KFLEV)=&
     & PDELP(KSTART:KPROF,1:KFLEV)*0.5_JPRB*(&
     & PU(KSTART:KPROF,1:KFLEV)**2+PV(KSTART:KPROF,1:KFLEV)**2)
    CALL ADD_FIELD_3D(YDLDDH,ZVAR,'VKK','V','ARP',.TRUE.,.TRUE.)
  ENDIF 
  
ENDIF

! set logical for delta_m
LLDPSFI=(NDPSFI==1)

! 0.3. initialize cumulative fields
!
ZLIFT=0.0_JPRB
ZFLUXU=0.0_JPRB
ZFLUXV=0.0_JPRB
ZFLUXH=0.0_JPRB
PTENDU=0.0_JPRB
PTENDV=0.0_JPRB
PTENDH=0.0_JPRB
DO JGFL=1,NUMFLDS
  IF (YCOMP(JGFL)%LT1) THEN
    ZFLUXGFL(:,:,YCOMP(JGFL)%MP1)=0.0_JPRB
    PTENDGFL(:,:,YCOMP(JGFL)%MP1)=0.0_JPRB
  ENDIF
ENDDO

! 0.4. calculation of temperatures on half levels
IF (LPHSPSH) THEN
  DO JROF = KSTART, KPROF
    ZTI(JROF,0)=PT(JROF,1)
    ZTI(JROF,KFLEV)=PT(JROF,KFLEV)
  ENDDO
ELSE
  DO JROF = KSTART, KPROF
    ZTI(JROF,0)=PT(JROF,1)
    ZTI(JROF,KFLEV)=PTS(JROF)
  ENDDO
ENDIF
DO JLEV= 1, KFLEV-1
  DO JROF = KSTART, KPROF
    ZTI(JROF,JLEV)=(PT(JROF,JLEV)+PT(JROF,JLEV+1))*0.5_JPRB
  ENDDO
ENDDO

! 0.5. calculate conversions factors between fluxes and tendencies
DO JLEV = 1, KFLEV
  DO JROF = KSTART, KPROF
    ZGSDP(JROF,JLEV)=RG*PRDELP(JROF,JLEV)
  ENDDO
ENDDO

! 0.6. mass fractions at half levels; only needed if delta_m=1
IF ( LLDPSFI ) THEN
  DO JGFL=1,NUMFLDS
    IF ( YCOMP(JGFL)%LWATER .AND. YCOMP(JGFL)%LT1 ) THEN
      DO JROF = KSTART, KPROF
        ZQXH (JROF,0,YCOMP(JGFL)%MP) = 0.0_JPRB
        ZQXH (JROF,KFLEV,YCOMP(JGFL)%MP) = PGFL (JROF,KFLEV,IPGFL(JGFL))
      ENDDO
      DO JLEV= 1, KFLEV-1
        DO JROF = KSTART, KPROF
          ZQXH(JROF,JLEV, YCOMP(JGFL)%MP) = (PGFL (JROF,JLEV,IPGFL(JGFL))+&
           & PGFL (JROF,JLEV+1,IPGFL(JGFL)))*0.5_JPRB
        ENDDO
      ENDDO
    ENDIF
  ENDDO
ENDIF

! 0.7. calculate kind-of total cp (including delta_m influence)
ZCP=RCPD
IF ( LLDPSFI ) THEN
  DO JGFL=1,NUMFLDS
    IF (YCOMP(JGFL)%LWATER .AND. YCOMP(JGFL)%LT1 ) THEN
      DO JLEV= 0, KFLEV
        DO JROF = KSTART, KPROF
          ZCP(JROF,JLEV)=ZCP(JROF,JLEV)+(YCOMP(JGFL)%RCP-RCPD)*ZQXH(JROF,JLEV,YCOMP(JGFL)%MP)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
ENDIF

! 
! 1. main loop over processes and fields
! 
DO JPROCESS=1,YDPROCSET%NPROCESS

  ! pointer to process
  YLPROC=>YDPROCSET%RPROCESS(JPROCESS)

  ! loop over fields belonging to this process
  DO JFIELD=1,YLPROC%NFIELD
  
    ! pointer to numerical values
    ZVAL=>YLPROC%RFIELD(JFIELD)%RVAL
    
    ! switch on type of variable
    SELECT CASE (YLPROC%RFIELD(JFIELD)%CVAR)
    
!
! 1.1. momentum u
!
      CASE ('U') ! momentum flux/tendency
      
        IF (YLPROC%RFIELD(JFIELD)%CTYPE == 'F' ) THEN ! flux
        
          ZFLUXU(KSTART:KPROF,:)=ZFLUXU(KSTART:KPROF,:)+ZVAL(KSTART:KPROF,:)
          
          IF (LFLEXDIA) THEN
            ! find corresponding V field
            DO JFIELD2=1,YLPROC%NFIELD
              IF (YLPROC%RFIELD(JFIELD2)%CVAR=='V'.AND.&
               & YLPROC%RFIELD(JFIELD)%CNAME==YLPROC%RFIELD(JFIELD2)%CNAME ) EXIT
            ENDDO
            IF (JFIELD2<=YLPROC%NFIELD) THEN
              ! transform to geographic north
              ZVAL2=>YLPROC%RFIELD(JFIELD2)%RVAL
              DO JLEV=0,KFLEV
                DO JROF=KSTART,KPROF
                  ZFLX(JROF,JLEV) = PGNORDM(JROF)*ZVAL(JROF,JLEV)-&
                   & PGNORDL(JROF)*ZVAL2(JROF,JLEV)
                  ZFLX2(JROF,JLEV) = PGNORDL(JROF)*ZVAL(JROF,JLEV)+&
                   & PGNORDM(JROF)*ZVAL2(JROF,JLEV)
                ENDDO
              ENDDO
              ! U and V fluxes
              IF (LDDH_OMP) THEN
                CALL NEW_ADD_FIELD_3D(YDMDDH,ZFLX,'FUU'//YLPROC%RFIELD(JFIELD)%CNAME,YDDDH)
                CALL NEW_ADD_FIELD_3D(YDMDDH,ZFLX2,'FVV'//YLPROC%RFIELD(JFIELD)%CNAME,YDDDH)
              ELSE
                CALL ADD_FIELD_3D(YDLDDH,ZFLX,'FUU'//YLPROC%RFIELD(JFIELD)%CNAME,'F','ARP',.TRUE.,.TRUE.)
                CALL ADD_FIELD_3D(YDLDDH,ZFLX2,'FVV'//YLPROC%RFIELD(JFIELD)%CNAME,'F','ARP',.TRUE.,.TRUE.)
              ENDIF
              ! kinetic energy flux
              DO JLEV=1,KFLEV
                DO JROF=KSTART,KPROF
                  ZVAR(JROF,JLEV)=-( ZUZGEO(JROF,JLEV)*&
                   & (ZFLX(JROF,JLEV)-ZFLX(JROF,JLEV-1)) + ZVMGEO(JROF,JLEV)&
                   & *(ZFLX2(JROF,JLEV)-ZFLX2(JROF,JLEV-1)) )
                ENDDO
              ENDDO
              IF (LDDH_OMP) THEN
                CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAR,'TKK'//YLPROC%RFIELD(JFIELD)%CNAME,YDDDH)
              ELSE
                CALL ADD_FIELD_3D(YDLDDH,ZVAR,'TKK'//YLPROC%RFIELD(JFIELD)%CNAME,'T','ARP',.TRUE.,.TRUE.)
              ENDIF
            ELSE
              WRITE (NULOUT,*) 'CPTEND_FLEX: did not find corresponding v-field'
            ENDIF
          ENDIF
          
        ELSE ! tendency
        
          PTENDU(KSTART:KPROF,:)=PTENDU(KSTART:KPROF,:)+ZVAL(KSTART:KPROF,:)

          IF (LFLEXDIA) THEN
            ! find corresponding V field
            DO JFIELD2=1,YLPROC%NFIELD
              IF (YLPROC%RFIELD(JFIELD2)%CVAR=='V'.AND.&
               & YLPROC%RFIELD(JFIELD)%CNAME==YLPROC%RFIELD(JFIELD2)%CNAME ) EXIT
            ENDDO
            IF (JFIELD2<=YLPROC%NFIELD) THEN
              ! transform to geographic north
              ZVAL2=>YLPROC%RFIELD(JFIELD2)%RVAL
              DO JLEV=1,KFLEV
                DO JROF=KSTART,KPROF
                  ZVAR(JROF,JLEV) = PGNORDM(JROF)*ZVAL(JROF,JLEV)-&
                   & PGNORDL(JROF)*ZVAL2(JROF,JLEV)
                  ZVAR2(JROF,JLEV) = PGNORDL(JROF)*ZVAL(JROF,JLEV)+&
                   & PGNORDM(JROF)*ZVAL2(JROF,JLEV)
                ENDDO
              ENDDO
              ! U and V tendencies
              IF (LDDH_OMP) THEN
                CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAR,'TUU'//YLPROC%RFIELD(JFIELD)%CNAME,YDDDH)
                CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAR2,'TVV'//YLPROC%RFIELD(JFIELD)%CNAME,YDDDH)
              ELSE
                CALL ADD_FIELD_3D(YDLDDH,ZVAR,'TUU'//YLPROC%RFIELD(JFIELD)%CNAME,'F','ARP',.TRUE.,.TRUE.)
                CALL ADD_FIELD_3D(YDLDDH,ZVAR2,'TVV'//YLPROC%RFIELD(JFIELD)%CNAME,'F','ARP',.TRUE.,.TRUE.)
              ENDIF
              ! kinetic energy tendency
              DO JLEV=1,KFLEV
                DO JROF=KSTART,KPROF
                  ZVAR(JROF,JLEV)=-( ZUZGEO(JROF,JLEV)*ZVAR(JROF,JLEV) +&
                   & ZVMGEO(JROF,JLEV)*ZVAR2(JROF,JLEV) )
                 ENDDO
              ENDDO
              IF (LDDH_OMP) THEN
                CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAR,'TKK'//YLPROC%RFIELD(JFIELD)%CNAME,YDDDH)
              ELSE
                CALL ADD_FIELD_3D(YDLDDH,ZVAR,'TKK'//YLPROC%RFIELD(JFIELD)%CNAME,'T','ARP',.TRUE.,.TRUE.)
              ENDIF
            ENDIF
          ENDIF
          
        ENDIF !

!
! 1.2. momentum v
!
      CASE ('V') ! momentum flux/tendency
      
        IF (YLPROC%RFIELD(JFIELD)%CTYPE == 'F' ) THEN
          ZFLUXV(KSTART:KPROF,:)=ZFLUXV(KSTART:KPROF,:)+ZVAL(KSTART:KPROF,:)
        ELSE
          PTENDV(KSTART:KPROF,:)=PTENDV(KSTART:KPROF,:)+ZVAL(KSTART:KPROF,:)
        ENDIF
        ! diagnostics are treated under U

!
! 1.3. enthalpy
!
      CASE ('H') ! enthalpy flux/tendency
      
        IF (YLPROC%RFIELD(JFIELD)%CTYPE == 'F' ) THEN
          DO JLEV= 0, KFLEV
            DO JROF = KSTART, KPROF
              ZFLUXH(JROF,JLEV)=ZFLUXH(JROF,JLEV)+ZVAL(JROF,JLEV)
            ENDDO
          ENDDO
          IF (LFLEXDIA) THEN
            IF (LDDH_OMP) THEN
              CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAL,'FCT'//YLPROC%RFIELD(JFIELD)%CNAME,YDDDH)
            ELSE
              CALL ADD_FIELD_3D(YDLDDH,ZVAL,'FCT'//YLPROC%RFIELD(JFIELD)%CNAME,'F','ARP',.TRUE.,.TRUE.)
            ENDIF
          ENDIF 
        ELSE
          DO JLEV= 1, KFLEV
            DO JROF = KSTART, KPROF
              PTENDH(JROF,JLEV)=PTENDH(JROF,JLEV)+ZVAL(JROF,JLEV)
            ENDDO
          ENDDO
          IF (LFLEXDIA) THEN
            IF (LDDH_OMP) THEN
              CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAL,'TCT'//YLPROC%RFIELD(JFIELD)%CNAME,YDDDH)
            ELSE
              CALL ADD_FIELD_3D(YDLDDH,ZVAL,'TCT'//YLPROC%RFIELD(JFIELD)%CNAME,'T','ARP',.TRUE.,.TRUE.)
            ENDIF
          ENDIF 
        ENDIF

!
! 1.4. GFL variable
!        
      CASE ('G')
        ! GFL variables: water and others
        JGFL=YLPROC%RFIELD(JFIELD)%JGFL

        SELECT CASE (YLPROC%RFIELD(JFIELD)%CTYPE)
        
          CASE ('P')  ! precipitation flux
          
            ! influence on enthalpy and lifting term
            IF ( LLDPSFI ) THEN
              ZFLX(KSTART:KPROF,:)=(YCOMP(JGFL)%RCP-ZCP(KSTART:KPROF,:))*&
               & ZVAL(KSTART:KPROF,:)*ZTI(KSTART:KPROF,:)
              ZLIFT(KSTART:KPROF,:)=ZLIFT(KSTART:KPROF,:)+ZVAL(KSTART:KPROF,:)
            ELSE
              ZFLX(KSTART:KPROF,:)=(YCOMP(JGFL)%RCP-RCPD)*&
               & ZVAL(KSTART:KPROF,:)*ZTI(KSTART:KPROF,:)
            ENDIF
            ZFLUXH(KSTART:KPROF,:)=ZFLUXH(KSTART:KPROF,:)+ZFLX(KSTART:KPROF,:)

            IF (LFLEXDIA) THEN
              IF (LDDH_OMP) THEN
                CALL NEW_ADD_FIELD_3D(YDMDDH,ZFLX,'FCT'//YLPROC%RFIELD(JFIELD)%CNAME,YDDDH)
              ELSE
                CALL ADD_FIELD_3D(YDLDDH,ZFLX,'FCT'//YLPROC%RFIELD(JFIELD)%CNAME,'F','ARP',.TRUE.,.TRUE.)
              ENDIF
            ENDIF
            
            ! influence on species tendency
            IF ( YCOMP(JGFL)%LT1 ) THEN
              ZFLUXGFL(KSTART:KPROF,:,YCOMP(JGFL)%MP1)=ZFLUXGFL(KSTART:KPROF,:,YCOMP(JGFL)%MP1)+ZVAL(KSTART:KPROF,:)
              IF (LFLEXDIA) THEN
                IF(LDDH_OMP) THEN
                  CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAL,'F'//CLFLEXDIA(YCOMP(JGFL)%MP1)//YLPROC%RFIELD(JFIELD)%CNAME,&
                   & YDDDH)
                ELSE
                  CALL ADD_FIELD_3D(YDLDDH,ZVAL,'F'//CLFLEXDIA(YCOMP(JGFL)%MP1)//YLPROC%RFIELD(JFIELD)%CNAME,&
                   & 'F','ARP',.TRUE.,.TRUE.)
                ENDIF
              ENDIF
            ENDIF

          CASE ('D','C')  ! diffusive flux or corrective flux
          
            IF ( YCOMP(JGFL)%LT1 ) THEN
              ZFLUXGFL(KSTART:KPROF,:,YCOMP(JGFL)%MP1)=ZFLUXGFL(KSTART:KPROF,:,YCOMP(JGFL)%MP1)+ZVAL(KSTART:KPROF,:)
              IF (LFLEXDIA) THEN
                IF(LDDH_OMP) THEN
                  CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAL,'F'//CLFLEXDIA(YCOMP(JGFL)%MP1)//YLPROC%RFIELD(JFIELD)%CNAME,&
                 & YDDDH)
                ELSE
                  CALL ADD_FIELD_3D(YDLDDH,ZVAL,'F'//CLFLEXDIA(YCOMP(JGFL)%MP1)//YLPROC%RFIELD(JFIELD)%CNAME,&
                 & 'F','ARP',.TRUE.,.TRUE.)
                ENDIF
              ENDIF
            ENDIF

          CASE ('S')  ! pseudoflux
          
            JGFLTARGET=YLPROC%RFIELD(JFIELD)%JGFLTARGET

            ! source species
            IF ( YCOMP(JGFL)%LT1 ) THEN
              ZFLUXGFL(KSTART:KPROF,:,YCOMP(JGFL)%MP1)=ZFLUXGFL(KSTART:KPROF,:,YCOMP(JGFL)%MP1)+ZVAL(KSTART:KPROF,:)
              IF (LFLEXDIA) THEN
                IF(LDDH_OMP) THEN
                  CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAL,'F'//CLFLEXDIA(YCOMP(JGFL)%MP1)//YLPROC%RFIELD(JFIELD)%CNAME,&
                 & YDDDH)
                ELSE
                  CALL ADD_FIELD_3D(YDLDDH,ZVAL,'F'//CLFLEXDIA(YCOMP(JGFL)%MP1)//YLPROC%RFIELD(JFIELD)%CNAME,&
                 & 'F','ARP',.TRUE.,.TRUE.)
                ENDIF
              ENDIF
            ENDIF
            
            ! target species
            IF ( YCOMP(JGFLTARGET)%LT1 ) THEN
              ZFLUXGFL(KSTART:KPROF,:,YCOMP(JGFLTARGET)%MP1)=ZFLUXGFL(KSTART:KPROF,:,YCOMP(JGFLTARGET)%MP1)&
               & - ZVAL(KSTART:KPROF,:)
              IF (LFLEXDIA) THEN
                IF(LDDH_OMP) THEN
                  ZFLX(KSTART:KPROF,:)=-ZVAL(KSTART:KPROF,:)
                  CALL NEW_ADD_FIELD_3D(YDMDDH,ZFLX,'F'//CLFLEXDIA(YCOMP(JGFLTARGET)%MP1)//YLPROC%RFIELD(JFIELD)%CNAME,&
                   & YDDDH)
                 ELSE
                  ZFLX(KSTART:KPROF,:)=-ZVAL(KSTART:KPROF,:)
                  CALL ADD_FIELD_3D(YDLDDH,ZFLX,'F'//CLFLEXDIA(YCOMP(JGFLTARGET)%MP1)//YLPROC%RFIELD(JFIELD)%CNAME,&
                   & 'F','ARP',.TRUE.,.TRUE.)
                ENDIF
              ENDIF
            ENDIF
            
            ! enthalpy
            ZFLX(KSTART:KPROF,:)=(YCOMP(JGFL)%RLZER-YCOMP(JGFLTARGET)%RLZER)*ZVAL(KSTART:KPROF,:)
            ZFLUXH(KSTART:KPROF,:)=ZFLUXH(KSTART:KPROF,:)+ZFLX(KSTART:KPROF,:)
            IF (LFLEXDIA) THEN
              IF(LDDH_OMP) THEN
               CALL NEW_ADD_FIELD_3D(YDMDDH,ZFLX,'FCT'//YLPROC%RFIELD(JFIELD)%CNAME,YDDDH)
              ELSE
                CALL ADD_FIELD_3D(YDLDDH,ZFLX,'FCT'//YLPROC%RFIELD(JFIELD)%CNAME,'F','ARP',.TRUE.,.TRUE.)
              ENDIF
            ENDIF 

          CASE ('T')  ! general GFL tendency
          
            IF (YCOMP(JGFL)%LT1) THEN
              ! species contribution
              PTENDGFL(KSTART:KPROF,:,YCOMP(JGFL)%MP1)=PTENDGFL(KSTART:KPROF,:,YCOMP(JGFL)%MP1) +&
               & ZVAL(KSTART:KPROF,:)
              ! diagnostic
              IF (LFLEXDIA) THEN
                IF(LDDH_OMP) THEN
                  CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAL,'T'//CLFLEXDIA(YCOMP(JGFL)%MP1)//YLPROC%RFIELD(JFIELD)%CNAME,&
                 & YDDDH)
                ELSE
                  CALL ADD_FIELD_3D(YDLDDH,ZVAL,'T'//CLFLEXDIA(YCOMP(JGFL)%MP1)//YLPROC%RFIELD(JFIELD)%CNAME,&
                 & 'T','ARP',.TRUE.,.TRUE.)
                ENDIF
              ENDIF
  
              ! enthalpy contribution
              IF (YCOMP(JGFL)%LWATER) THEN
                ZVAR(KSTART:KPROF,:)=YCOMP(JGFL)%RLZER*ZVAL(KSTART:KPROF,:)
                PTENDH(KSTART:KPROF,:)=PTENDH(KSTART:KPROF,:)+ZVAR(KSTART:KPROF,:)
                IF (LFLEXDIA) THEN
                  IF(LDDH_OMP) THEN
                    CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAR,'T'//CLFLEXDIA(YCOMP(JGFL)%MP1)//&
                    &YLPROC%RFIELD(JFIELD)%CNAME,YDDDH)
                  ELSE
                    CALL ADD_FIELD_3D(YDLDDH,ZVAR,'T'//CLFLEXDIA(YCOMP(JGFL)%MP1)//&
                    &YLPROC%RFIELD(JFIELD)%CNAME,'T','ARP',.TRUE.,.TRUE.)
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
            
          CASE ('F') ! general GFL flux
            IF ( YCOMP(JGFL)%LT1 ) THEN
              ZFLUXGFL(KSTART:KPROF,:,YCOMP(JGFL)%MP1)=ZFLUXGFL(KSTART:KPROF,:,YCOMP(JGFL)%MP1) +&
                   & ZVAL(KSTART:KPROF,:)
              IF (LFLEXDIA) THEN
                IF(LDDH_OMP) THEN
                  CALL NEW_ADD_FIELD_3D(YDMDDH,ZVAR,'F'//CLFLEXDIA(YCOMP(JGFL)%MP1)//&
                   &YLPROC%RFIELD(JFIELD)%CNAME,YDDDH)
                 ELSE
                  CALL ADD_FIELD_3D(YDLDDH,ZVAR,'F'//CLFLEXDIA(YCOMP(JGFL)%MP1)//&
                   &YLPROC%RFIELD(JFIELD)%CNAME,'F','ARP',.TRUE.,.TRUE.)
                ENDIF
              ENDIF 
            ENDIF
            
          CASE DEFAULT
            CALL ABOR1('CPTEND_FLEX : no such field type')
            
        END SELECT
        
      CASE DEFAULT
        CALL ABOR1('CPTEND_FLEX : no such variable type')
        
    END SELECT
  ENDDO ! loop over fields
ENDDO   ! loop over processes

! 
! 2. treatment of ZLIFT
! 

! contribution of compensating flux (zlift) to humidities
IF ( LLDPSFI ) THEN
  DO JGFL=1,NUMFLDS
    IF ( YCOMP(JGFL)%LWATER .AND. YCOMP(JGFL)%LT1 ) THEN
      ZFLX(KSTART:KPROF,:)=-ZQXH(KSTART:KPROF,:,YCOMP(JGFL)%MP)*ZLIFT(KSTART:KPROF,:)
      ZFLUXGFL(KSTART:KPROF,:,YCOMP(JGFL)%MP1)=ZFLUXGFL(KSTART:KPROF,:,YCOMP(JGFL)%MP1)+ZFLX(KSTART:KPROF,:)
      IF (LFLEXDIA) THEN
        IF(LDDH_OMP) THEN
          CALL NEW_ADD_FIELD_3D(YDMDDH,ZFLX,'F'//CLFLEXDIA(YCOMP(JGFL)%MP1)//&
         &'LIFT',YDDDH)
        ELSE
          CALL ADD_FIELD_3D(YDLDDH,ZFLX,'F'//CLFLEXDIA(YCOMP(JGFL)%MP1)//&
         &'LIFT','F','ARP',.TRUE.,.TRUE.)
        ENDIF
      ENDIF
    ENDIF
  ENDDO
ENDIF

! 
! 3. transform fluxes into tendencies
! 

! u, v, T
DO JLEV = 1, KFLEV
  DO JROF = KSTART, KPROF
    PTENDU(JROF,JLEV) = PTENDU(JROF,JLEV) + ZGSDP(JROF,JLEV)*(&
     & (ZFLUXU (JROF,JLEV-1) - ZFLUXU (JROF,JLEV)))
    PTENDV(JROF,JLEV) = PTENDV(JROF,JLEV) + ZGSDP(JROF,JLEV)*(&
     & (ZFLUXV (JROF,JLEV-1) - ZFLUXV (JROF,JLEV)))
    PTENDH(JROF,JLEV) = PTENDH(JROF,JLEV) + ZGSDP(JROF,JLEV)*(&
     &  ZFLUXH  (JROF,JLEV-1) - ZFLUXH  (JROF,JLEV))
  ENDDO
ENDDO

! GFL
DO JGFL=1,NUMFLDS
  IF ( YCOMP(JGFL)%LT1 ) THEN
    DO JLEV= 1, KFLEV
      DO JROF = KSTART, KPROF
        PTENDGFL(JROF,JLEV,YCOMP(JGFL)%MP1)=PTENDGFL(JROF,JLEV,YCOMP(JGFL)%MP1)+ZGSDP(JROF,JLEV)*(&
         & ZFLUXGFL(JROF,JLEV-1,YCOMP(JGFL)%MP1)-ZFLUXGFL(JROF,JLEV,YCOMP(JGFL)%MP1) )
      ENDDO
    ENDDO
  ENDIF
ENDDO

! calculation of pfp (total relative precip. flux), required further in mfphys ...
IF ( LLDPSFI .AND. PRESENT (PFP)) THEN
  PFP=1.0_JPRB
  DO JGFL=1,NUMFLDS
    ! daand: not sure if LPRECIP key should be used here (or anywhere else!?)
    !IF (YCOMP(JGFL)%LWATER .AND. YCOMP(JGFL)%LT1 .AND. YCOMP(JGFL)%LPRECIP ) THEN
    IF (YCOMP(JGFL)%LWATER .AND. YCOMP(JGFL)%LT1 ) THEN
      DO JLEV= 0, KFLEV
        DO JROF = KSTART, KPROF
          PFP(JROF,JLEV)=PFP(JROF,JLEV)-ZQXH(JROF,JLEV,YCOMP(JGFL)%MP)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  DO JLEV= 0, KFLEV
    DO JROF = KSTART, KPROF
      PFP(JROF,JLEV)=PFP(JROF,JLEV)*ZLIFT(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

! 
! 4. diagnostics
! 

IF (LSDDH) THEN
  ! diagnostics, these should be replaced by add_field_3d calls in the main loop over processes above
  PFHSCL(:,:) = 0.0_JPRB
  PFHSCN(:,:) = 0.0_JPRB
  PFHSSL(:,:) = 0.0_JPRB
  PFHSSN(:,:) = 0.0_JPRB
  PFHPCL(:,:) = 0.0_JPRB
  PFHPCN(:,:) = 0.0_JPRB
  PFHPSL(:,:) = 0.0_JPRB
  PFHPSN(:,:) = 0.0_JPRB
  PFEPFP(:,:) = 0.0_JPRB
  IF (PRESENT (PFHP)) PFHP(:,:)   = 0.0_JPRB
  PFCMPCQ(:,:)= 0.0_JPRB
  PFCMPSN(:,:)= 0.0_JPRB
  PFCMPSL(:,:)= 0.0_JPRB

  ! contributions of precipitating and diffusive species to total cp
  ZCPP=0.0_JPRB
  ZCPD=0.0_JPRB
  IF ( LLDPSFI ) THEN
    DO JGFL=1,NUMFLDS
      IF ( YCOMP(JGFL)%LWATER .AND. YCOMP(JGFL)%LT1 ) THEN
        IF ( YCOMP(JGFL)%LPRECIP )  THEN
          DO JLEV= 0, KFLEV
            DO JROF = KSTART, KPROF
              ZCPP(JROF,JLEV)=ZCPP(JROF,JLEV)+(YCOMP(JGFL)%RCP-RCPD)*ZQXH(JROF,JLEV,YCOMP(JGFL)%MP)
            ENDDO
          ENDDO
        ELSE
          DO JLEV= 0, KFLEV
            DO JROF = KSTART, KPROF
              ZCPD(JROF,JLEV)=ZCPD(JROF,JLEV)+(YCOMP(JGFL)%RCP-RCPD)*ZQXH(JROF,JLEV,YCOMP(JGFL)%MP)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
  
  !reloop over processes
  DO JPROCESS=1,YDPROCSET%NPROCESS
    YLPROC=>YDPROCSET%RPROCESS(JPROCESS)
    DO JFIELD=1,YLPROC%NFIELD
      ZVAL=>YLPROC%RFIELD(JFIELD)%RVAL
      IF (YLPROC%RFIELD(JFIELD)%CVAR=='G') THEN
        JGFL=YLPROC%RFIELD(JFIELD)%JGFL
        IF (YLPROC%RFIELD(JFIELD)%CTYPE=='P') THEN
          ! precipitation flux
          DO JLEV= 0, KFLEV
            DO JROF = KSTART, KPROF
              ZJTMP(JROF,JLEV)=(YCOMP(JGFL)%RCP-RCPD-ZCPP(JROF,JLEV))*&
               & ZVAL(JROF,JLEV)*ZTI(JROF,JLEV)
            ENDDO
          ENDDO
          IF (LSDDH.AND.YLPROC%RFIELD(JFIELD)%CNAME=='PFPLCL') PFHSCL=PFHSCL+ZJTMP
          IF (YLPROC%RFIELD(JFIELD)%CNAME=='PFPLCN') PFHSCN=PFHSCN+ZJTMP
          IF (YLPROC%RFIELD(JFIELD)%CNAME=='PFPLSL') PFHSSL=PFHSSL+ZJTMP
          IF (YLPROC%RFIELD(JFIELD)%CNAME=='PFPLSN') PFHSSN=PFHSSN+ZJTMP
          IF (PRESENT (PFHP)) PFHP=PFHP+ZJTMP
        ELSEIF (YLPROC%RFIELD(JFIELD)%CTYPE=='S') THEN
          ! pseudoflux
          JGFLTARGET=YLPROC%RFIELD(JFIELD)%JGFLTARGET
          DO JLEV= 0, KFLEV
            DO JROF = KSTART, KPROF
              ZJTMP(JROF,JLEV)=(YCOMP(JGFL)%RLZER-YCOMP(JGFLTARGET)%RLZER)*&
               & ZVAL(JROF,JLEV)
            ENDDO
          ENDDO  
          IF (YLPROC%RFIELD(JFIELD)%CNAME=='PFCCQL') PFHPCL=PFHPCL+ZJTMP
          IF (YLPROC%RFIELD(JFIELD)%CNAME=='PFECL')  PFHPCL=PFHPCL+ZJTMP
          IF (YLPROC%RFIELD(JFIELD)%CNAME=='PFCSQL') PFHPSL=PFHPSL+ZJTMP
          IF (YLPROC%RFIELD(JFIELD)%CNAME=='PFESL')  PFHPSL=PFHPSL+ZJTMP
          IF (YLPROC%RFIELD(JFIELD)%CNAME=='PFCCQN') PFHPCN=PFHPCN+ZJTMP
          IF (YLPROC%RFIELD(JFIELD)%CNAME=='PFECN')  PFHPCN=PFHPCN+ZJTMP
          IF (YLPROC%RFIELD(JFIELD)%CNAME=='PFCSQN') PFHPSN=PFHPSN+ZJTMP
          IF (YLPROC%RFIELD(JFIELD)%CNAME=='PFESN')  PFHPSN=PFHPSN+ZJTMP
          IF (PRESENT (PFHP)) PFHP=PFHP+ZJTMP
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  ! contribution of zlift to enthalpy
  IF (LLDPSFI) THEN
    DO JLEV= 0, KFLEV
      DO JROF = KSTART, KPROF
        PFEPFP(JROF,JLEV)=-ZCPD(JROF,JLEV)*ZLIFT(JROF,JLEV)*ZTI(JROF,JLEV)
        IF (PRESENT (PFHP)) PFHP(JROF,JLEV)=PFHP(JROF,JLEV)+PFEPFP(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF

  ! compensating fluxes of diff. species
  IF ( LLDPSFI ) THEN
    DO JLEV= 0, KFLEV
      DO JROF = KSTART, KPROF
        PFCMPCQ(JROF,JLEV)=ZQXH(JROF,JLEV,YQ%MP)*ZLIFT(JROF,JLEV)
      ENDDO
    ENDDO
    DO JLEV= 0, KFLEV
      DO JROF = KSTART, KPROF
        PFCMPSN(JROF,JLEV)=ZQXH(JROF,JLEV,YI%MP)*ZLIFT(JROF,JLEV)
      ENDDO
    ENDDO
    DO JLEV= 0, KFLEV
      DO JROF = KSTART, KPROF
        PFCMPSL(JROF,JLEV)=ZQXH(JROF,JLEV,YL%MP)*ZLIFT(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF
  
ENDIF !LSDDH

! 
! 5. other business: should this be done in cptend_flex ???
! 
IF(LMUSCLFA) THEN

  !-------------------------------------------------
  ! Single column model MUSC diagnostics.
  ! Q1 and Q2.
  !-------------------------------------------------

  DO JLEV=1,KFLEV
    DO JROF=KSTART,KPROF
      ZQ1(JROF,JLEV)=(PTENDH(JROF,JLEV)-(RCPV-RCPD)*PTENDGFL(JROF,JLEV,YQ%MP1)*PT(JROF,JLEV))/PCP(JROF,JLEV)
      ZQ2(JROF,JLEV)=-FOLH(PT(JROF,JLEV),MAX(0._JPRB,SIGN(1._JPRB,RTT-PT(JROF,JLEV))))&
        & *PTENDGFL(JROF,JLEV,YQ%MP1)/PCP(JROF,JLEV)
    ENDDO
  ENDDO
  CALL WRSCMR(NMUSCLFA,'Q1',ZQ1,KPROF,KFLEV)
  CALL WRSCMR(NMUSCLFA,'Q2',ZQ2,KPROF,KFLEV)

ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPTEND_FLEX',1,ZHOOK_HANDLE)
END SUBROUTINE CPTEND_FLEX
