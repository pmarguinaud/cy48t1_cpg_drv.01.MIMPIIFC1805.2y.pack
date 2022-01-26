SUBROUTINE DIAGFLASH(YDCFU,KST,KEND,KPROMA,KLEV,KSTEP,&
                     & PQL,PQI,PQR,PQS,PQG,PQH,PDELP,PT,PVW,PLSM,PFLASH)


!******** FPDIAGFLASH  ************
!
!      PURPOSE:
!      --------
!      Compute lightning densities
!
!      INTERFACE:
!      ----------     
!      *CALL DIAGFLASH*
!      
!
!        EXPLICIT ARGUMENTS:
!        -------------------
!          INPUT:
!        KST     : start of work
!        KEND    : end of work
!        KPROMA  : dimension of work
!        KLEV    : number of levels
!        KSTEP   : actual step of model
!        PQL     : cloud water (g/kg)
!        PQI     : cloud ice (g/kg)
!        PQR     : rain (g/kg)
!        PQS     : snow (g/kg) 
!        PQG     : graupel (and hail) (g/kg)
!        PQH     : hail (g/kg)
!        PDELP   : layer thickness (Pa) 
!        PT      : temperature (K)
!        PVW     : vertical velocity (m/s)
!        PLSM    : land sea mask
!
!          OUTPUT:
!        PFLASH  : flash density for given column
!
!        IMPLICIT ARGUMENTS:
!        -------------------
!           NONE
!
!      METHOD:
!      -------
!        different methods to compute flash rate can be chosen, dependent on  
!        NMTFLASH:
!        
!         = 1  McCaul et. al  (2009) part 1: flash rate is based on upward flux of 
!              precipitation ice hydrometeors (graupel [and hail if available])
!               in the mixed phase region
!	  = 2  McCaul et. al  (2009) part 2: flash rate is based on vertically integrated amount of ice hydrometeors              
!         = 3  McCaul et. al  (2009) part 3: blending version 1 and 2
!
!      REFERENCE:
!      ----------
!
!      NMTFLASH=1-3: McCaul Jr, Eugene W., et al. "Forecasting lightning threat using cloud-resolving model simulations." 
!                    Weather and Forecasting 24.3 (2009): 709-729.McCaul et al.    
!
!      AUTHOR:
!      -------
!        Christoph Wittmann *ZAMG*
!
!      MODIFICATIONS:
!      --------------
!
!        original version: June 2013
!


USE PARKIND1, ONLY: JPIM, JPRB
USE YOMHOOK,  ONLY: LHOOK, DR_HOOK
USE YOMCST,   ONLY: RTT, RG 
USE YOMCFU   , ONLY : TCFU        


IMPLICIT NONE

TYPE(TCFU),INTENT(IN) :: YDCFU
INTEGER(KIND=JPIM),INTENT(IN) :: KST,KEND,KPROMA,KLEV,KSTEP
REAL(KIND=JPRB),INTENT(IN) :: PQI(KPROMA,KLEV),PQL(KPROMA,KLEV),PQS(KPROMA,KLEV),PQR(KPROMA,KLEV)  
REAL(KIND=JPRB),INTENT(IN) :: PQG(KPROMA,KLEV),PQH(KPROMA,KLEV), PDELP(KPROMA,KLEV)        
REAL(KIND=JPRB),INTENT(IN) :: PVW(KPROMA,KLEV),PLSM(KPROMA),PT(KPROMA,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PFLASH(KPROMA)

!local declarations
!INTEGER (KIND=JPIM),POINTER ::NMTFLASH
INTEGER (KIND=JPIM) :: JLON, JLEV
REAL(KIND=JPRB) :: ZHOOK_HANDLE
REAL (KIND=JPRB) :: ZINT(KPROMA), ZFLASH1(KPROMA), ZFLASH2(KPROMA)
REAL (KIND=JPRB) :: ZTEST, ZTEST1, ZTEST2, ZTEST3 
REAL (KIND=JPRB) :: ZQILMIN, ZTMIX, ZVW 

!additional declarations/thresholds -> should be moved to namelist later
REAL (KIND=JPRB) :: ZRFLASHQIL,ZRFLASHTMXPL,ZRFLASHTCRIT  !!!!!ZRFLASHCAL


IF (LHOOK) CALL DR_HOOK('DIAGFLASH',0,ZHOOK_HANDLE)
! ------------------------------------------------

ASSOCIATE(NMTFLASH=>YDCFU%NMTFLASH,CALFLASH1=>YDCFU%CALFLASH1,CALFLASH2=>YDCFU%CALFLASH2)

ZFLASH1(:)=0.0_JPRB
ZFLASH2(:)=0.0_JPRB

ZRFLASHQIL=1.E-9_JPRB
ZRFLASHTMXPL=15.0_JPRB
ZRFLASHTCRIT=-5.0_JPRB
ZQILMIN=ZRFLASHQIL      ! min value for qi+ql to detect "inside cloud" 

! -> to get flashes / ( km**2 * model_time_step) we have multiplication by TSTEP in  cpcfu.F90 !!!!!!!!

ZTMIX=RTT-ZRFLASHTMXPL   

IF (NMTFLASH == 1 .OR. NMTFLASH == 3) THEN
  ZINT(:)=0.0_JPRB   

  ! first: find level closest to -15 degree isotherm inside a cloud
  DO JLEV=KST+2,KLEV,1
    DO JLON=KST,KEND
      
      ! test 1: inside a cloud? , test 2: temperture below -15 , test 3: temperatur above -15
      ! test= 1+2+3 should garantue that we've just passed the -15 degree isoterm inside a cloud
      
      ZTEST1=MAX(0.0_JPRB,SIGN(1.0_JPRB,PQI(JLON,JLEV)+PQL(JLON,JLEV)-ZQILMIN))         
      ZTEST2=MAX(0.0_JPRB,SIGN(1.0_JPRB,ZTMIX-PT(JLON,JLEV-1)))   
      ZTEST3=MAX(0.0_JPRB,SIGN(1.0_JPRB,PT(JLON,JLEV)-ZTMIX)) 
      ZTEST=ZTEST1*ZTEST2*ZTEST3
  
      ! if all test are successful we save the level (just below -15) inside a cloud we need 
      ! vertical velocity and graupe+hail content
      ! therefore we take avaerage from level below and above isotherm
      ! and multiply with (positiv) vertical velocity (average from level below and above isotherm)
      
      ZVW=MAX(0.0_JPRB,0.5_JPRB*(PVW(JLON,JLEV-1)+PVW(JLON,JLEV)))
 
      !the following line should be active max 1 time during vertical loop
     
      ZINT(JLON)=ZINT(JLON)+ZTEST*0.5_JPRB*(PQG(JLON,JLEV-1)+PQG(JLON,JLEV)+PQH(JLON,JLEV-1)+PQH(JLON,JLEV))*ZVW
     
    ENDDO
  ENDDO

  DO JLON=KST,KEND

    ZFLASH1(JLON)=CALFLASH1*ZINT(JLON)
  ENDDO
ENDIF  

IF (NMTFLASH == 2 .OR. NMTFLASH == 3) THEN
  ZINT(:)=0.0_JPRB
  ! we need to integrate solid hydrometeors inside a cloud
  DO JLEV=KLEV,1,-1
    DO JLON=KST, KEND
      ! if sum of cloud condensate exceeds a certain threshold we assume that we're inside a cloud
      ZTEST1=MAX(0.0_JPRB,SIGN(1.0_JPRB,PQI(JLON,JLEV)+PQL(JLON,JLEV)-ZQILMIN))
      ZINT(JLON)=ZINT(JLON)+ZTEST1*PDELP(JLON,JLEV)/RG*&
                 &(PQI(JLON,JLEV)+PQS(JLON,JLEV)+PQG(JLON,JLEV)+PQH(JLON,JLEV))
    ENDDO
  ENDDO
  
  DO JLON=KST,KEND
    ZFLASH2(JLON)=CALFLASH2*ZINT(JLON)
  ENDDO

ENDIF  ! end of NTMFLASH

! final calculation:

IF (NMTFLASH == 3 ) THEN
        ! McCaul version 3 is simply a weighted combination of version 1 and 2 (-> McCaul 2009)
  DO JLON=KST,KEND
    PFLASH(JLON)=0.95_JPRB*ZFLASH1(JLON)+0.05_JPRB*ZFLASH2(JLON)
  ENDDO
 
ELSE

  DO JLON=KST,KEND
    PFLASH(JLON)=ZFLASH1(JLON)+ZFLASH2(JLON) 
  ENDDO
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DIAGFLASH',1,ZHOOK_HANDLE)


END SUBROUTINE DIAGFLASH 
