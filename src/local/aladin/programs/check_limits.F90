PROGRAM CHECK_LIMITS 

!      -----------------------------------------------------------------
! program for checking limits in aladin files
! before running this program copy fnam1 onto FNAME3!
! On output FNAME3 will contain the checking values of FNAME1
! for some surface variables
!      -----------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB, JPRD
USE YOMLUN   , ONLY : NULOUT,  NULNAM
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK

USE YOMRIP0  , ONLY : RTIMST 

!      -----------------------------------------------------------------

IMPLICIT NONE

CHARACTER (LEN = 14) ::  CL_FNAME1,CL_FNAME3
INTEGER(KIND=JPIM), PARAMETER      ::   JPMXLEV=200
    
INTEGER(KIND=JPIM) :: I_SSABN
INTEGER(KIND=JPIM) :: INXLON,I
INTEGER(KIND=JPIM) :: IDATEF1(11),IDATEF3(11),ILMO(12)
INTEGER(KIND=JPIM) :: INLOPA(322)
INTEGER(KIND=JPIM) :: ITYPTR, ITRONC
INTEGER(KIND=JPIM) :: INLATI, INIVER
INTEGER(KIND=JPIM) :: ISSSS, IUSSSS, ISHOUR, ISDAY, INC, IUDATE
INTEGER(KIND=JPIM) :: IA, IM, ID, INRGRPOI, I_NDGL, I_NDLON, I_NTVMER,I_NTVGLA
INTEGER(KIND=JPIM) :: IREP,INMAX
         
REAL(KIND=JPRB)    ::   RDAY
REAL(KIND=JPRB)    :: Z_RINC, ZVT1,ZVT3,ZREF
REAL(KIND=JPRB)    :: ZSINLA(322),ZVALH(0:100),ZVBH(0:100)     
REAL(KIND=JPRB)    :: ZSUR1(:),ZSUR3(:),Z_WSAT(:),ZD2(:),Z_SAB(:),Z_VEG(:)
REAL(KIND=JPRB)    :: Z_LAI(:),ZLSM(:),ZIVEG(:),Z_WPMX(:),Z_WSMX(:),Z_WP(:),Z_WS(:) 
REAL(KIND=JPRB)    :: Z_GCONV,Z_GWPIMX,Z_RD1,Z_WPI
REAL(KIND=JPRB)    :: Z_WSI
REAL(KIND=JPRB)    :: Z_PWFC,Z_G1WSAT,Z_G2WSAT,ZSAB
REAL(KIND=JPRB)    :: Z_SZZ0N,Z_SDEPN,ZEPS,ZPEPS
         
ALLOCATABLE :: ZSUR1,ZSUR3,Z_WSAT,ZLSM,ZIVEG,ZD2,Z_SAB,Z_VEG,Z_LAI,Z_WPMX,Z_WSMX,Z_WP,Z_WS
          
LOGICAL :: LLERR,LLMAP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
         
!      -----------------------------------------------------------------
NAMELIST/NAMCHECK_LIMITS/ CL_FNAME1,CL_FNAME3 
!      -----------------------------------------------------------------

#include "echien.h"
#include "sucst_ifsaux.h"

#include "sudyncore.intfb.h"
#include "updcal.intfb.h"

#include "fcttim.func.h"

!      -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CHECK_LIMITS',0,ZHOOK_HANDLE)
!      -----------------------------------------------------------------

ISHOUR=3600
ISDAY=3600*24

ZPEPS=1.0E-7_JPRB
ZEPS=1.0E-1_JPRB

I_NTVGLA=2
I_NTVMER=1
Z_PWFC=1._JPRB
Z_G1WSAT=-.1080E-02_JPRB
Z_G2WSAT=0.4943E+00_JPRB 

!      constants from YOMCLI 

Z_SZZ0N=0.001_JPRB 
Z_SDEPN=0.10_JPRB
I_SSABN=6._JPRB

!      constants from YOMPHY1
 
Z_GCONV=1.E+3_JPRB
Z_GWPIMX=150._JPRB
Z_RD1=1.E-2_JPRB

!  reading namelist (file names and zsign to have adding or subtracting) 

READ (4,NAMCHECK_LIMITS)

!  opening fa files

INMAX=1
       
CALL FAITOU(IREP,11,.TRUE.,CL_FNAME1,'OLD',.TRUE.,.FALSE.,0,INMAX,&
 & INMAX,'CADRE LECTURE   ')  
       
IF (IREP/=0) THEN
       
  WRITE(*,*) 'pb with opening ',CL_FNAME1,' error code:',IREP
  STOP
       
ENDIF 
       
INMAX=1 
       
CALL FAITOU(IREP,13,.TRUE.,CL_FNAME3,'OLD',.TRUE.,.FALSE.,0,INMAX,&
 & INMAX,'CADRE LECTURE1  ')  
 
IF (IREP/=0) THEN
 
  WRITE(*,*) 'pb with opening ',CL_FNAME3,' error code:',IREP
  STOP
 
ENDIF
 
WRITE(*,*) 'Pass faitou'
          
!      checking validity time of files

CALL FADIES(IREP,11,IDATEF1)
 
IF (IREP/=0) THEN
 
  WRITE(*,*) 'pb with fadies ',CL_FNAME1,' error code:',IREP
  STOP
 
ENDIF
 
ISSSS=IDATEF1(7)*ISHOUR+IDATEF1(4)*ISHOUR+IDATEF1(5)*60
IUSSSS=MOD(ISSSS,ISDAY)
Z_RINC=REAL(ISSSS,JPRB)/REAL(ISDAY,JPRB)
INC=INT(Z_RINC)
 
IF (INC /= 0) THEN
 
  CALL UPDCAL(IDATEF1(3),IDATEF1(2),IDATEF1(1),INC,&
   & IDATEF1(3),IDATEF1(2),IDATEF1(1),ILMO,6)  
 
ENDIF
 
IUDATE=IDATEF1(1)*10000+IDATEF1(2)*100+IDATEF1(3)
 
LLERR=.TRUE.
OPEN (NULNAM,FILE='fort.4',STATUS='OLD',FORM='FORMATTED',ERR=904)
LLERR=.FALSE.
904  CONTINUE
IF (LLERR) THEN
  WRITE(NULOUT,'(''check_limits : namelist - problem to open fort.4'')')
  CALL ABOR1('CHECK_LIMITS : namelist - problem to open fort.4')
ENDIF

CALL SUDYNCORE          ! yomdyncore
RDAY=86400._JPRB
CLOSE(NULNAM)

CALL SUCST_IFSAUX
ZVT1 = RTIMST
CALL FADIES(IREP,13,IDATEF3)
 
IF (IREP/=0) THEN
 
  WRITE(*,*) 'pb with fadies ',CL_FNAME3,' error code:',IREP
  STOP
 
ENDIF
 
ISSSS=IDATEF3(7)*ISHOUR+IDATEF3(4)*ISHOUR+IDATEF3(5)*60
IUSSSS=MOD(ISSSS,ISDAY)
Z_RINC=REAL(ISSSS,JPRB)/REAL(ISDAY,JPRB)
INC=INT(Z_RINC)
 
IF (INC /= 0) THEN
 
  CALL UPDCAL(IDATEF3(3),IDATEF3(2),IDATEF3(1),INC,&
   & IDATEF3(3),IDATEF3(2),IDATEF3(1),ILMO,6)  
 
ENDIF
 
IUDATE=IDATEF3(1)*10000+IDATEF3(2)*100+IDATEF3(3)
ID=NDD(IUDATE)
IM=NMM(IUDATE)
IA=NCCAA(IUDATE)
 
ZVT3 = RTIME(IA,IM,ID,IUSSSS,RDAY)
 
WRITE(*,*) 'validity times:'
WRITE(*,*) ZVT1
WRITE(*,*) ZVT3
 
IF (ZVT1/=ZVT3) THEN
 
  WRITE(*,*) 'validity times are not the same'
  STOP
 
ENDIF

!      checking cadre
INIVER=JPMXLEV
       
CALL ECHIEN('CADRE LECTURE   ',ITYPTR,LLMAP,ITRONC,INLATI,&
 & INXLON,INLOPA,ZSINLA,INIVER,ZREF,ZVALH,ZVBH,1,ZPEPS,6)  
      
I_NDGL=INLATI
I_NDLON=INXLON
           
WRITE(*,*) 'Pass echien CADRE LECTURE '
       
CALL ECHIEN('CADRE LECTURE1  ',ITYPTR,LLMAP,ITRONC,INLATI,&
 & INXLON,INLOPA,ZSINLA,INIVER,ZREF,ZVALH,ZVBH,0,ZPEPS,6)  

WRITE(*,*) 'Pass echien CADRE LECTURE1 '

!      Number of grid points

INRGRPOI=I_NDGL*I_NDLON
         
ALLOCATE(ZSUR1(INRGRPOI),ZSUR3(INRGRPOI),Z_WP(INRGRPOI),Z_WS(INRGRPOI),Z_WPMX(INRGRPOI),ZLSM(INRGRPOI))
ALLOCATE(Z_WSMX(INRGRPOI),Z_SAB(INRGRPOI),Z_VEG(INRGRPOI),Z_WSAT(INRGRPOI),ZD2(INRGRPOI),ZIVEG(INRGRPOI))
ALLOCATE(Z_LAI(INRGRPOI))

!      0.1 reading ZLSM land-sea mask from 'SURFIND.TERREMER'

ZLSM=0._JPRB

CALL FACILE(IREP,11,'SURF',1,'IND.TERREMER',ZLSM,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in reading surf IND.TERREMER in file1',IREP
  STOP
 
ENDIF

!      0.2 reading VEG vegetation fraction from 'SURFPROP.VEGETAT'

Z_VEG=0._JPRB

CALL FACILE(IREP,11,'SURF',1,'PROP.VEGETAT',Z_VEG,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in reading surf PROP.VEGETAT in file1',IREP
  STOP

ENDIF

!      0.2.1 reading ZIVEG real va from 'SURFIND.VEG.DOMI'

ZIVEG=0._JPRB

CALL FACILE(IREP,11,'SURF',1,'IND.VEG.DOMI',ZIVEG,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in reading SURF IND.VEG.DOMI in file1',IREP
  STOP

ENDIF

!      0.2.4 reading SAB percentage of sand  SURFPROP.SABLE

Z_SAB=I_SSABN

CALL FACILE(IREP,11,'SURF',1,'PROP.SABLE  ',Z_SAB,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in reading SURF PROP.SABLE   in file1',IREP
  STOP

ENDIF

!      0.2.5 reading ZD2 useful soil depth SURFEPAIS.SOL      

ZD2=Z_SDEPN

CALL FACILE(IREP,11,'SURF',1,'EPAIS.SOL   ',ZD2,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in reading SURF EPAIS.SOL    in file1',IREP
  STOP

ENDIF

!      0.2.6 reading LAI leaf-area index SURFIND.FOLIAIRE 

Z_LAI=0._JPRB

CALL FACILE(IREP,11,'SURF',1,'IND.FOLIAIRE',Z_LAI,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in reading SURF IND.FOLIAIRE in file1',IREP
  STOP

ENDIF

!       calculate  wsat(SAB,ZLSM,ZIVEG) in this program wsat(i)       

DO I=1,INRGRPOI

  IF ((ZLSM(I)<0.5_JPRB).OR.(ZIVEG(I)==REAL(I_NTVGLA,JPRB))) THEN

    Z_WSAT(I)=Z_PWFC

  ELSE

    ZSAB=MAX(ZEPS,Z_SAB(I))
    Z_WSAT(I)=Z_G1WSAT*ZSAB+Z_G2WSAT 

  ENDIF

ENDDO
      
ZSUR1=0._JPRB
ZSUR3=0._JPRB
       
CALL FACILE(IREP,11,'PROF',1,'RESERV.EAU  ',ZSUR1,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in reading prof RESERV.EAU in file1',IREP
  STOP

ENDIF

DO I=1,INRGRPOI

  Z_WP(I)=ZSUR1(I)
             
  Z_WPMX(I)=Z_GCONV*ZD2(I)*Z_WSAT(I)

  Z_WP(I)=MAX(0._JPRB,MIN(Z_WP(I),Z_WPMX(I)))
                 
  IF (NINT(ZIVEG(I))==I_NTVMER) THEN
                     
    Z_WP(I)=Z_WPMX(I)
                    
  ELSEIF (NINT(ZIVEG(I))==I_NTVGLA) THEN
                 
    Z_WP(I)=0._JPRB
                    
  ENDIF

  ZSUR3(I)=Z_WP(I)  

ENDDO 

CALL FAIENC(IREP,13,'PROF',1,'RESERV.EAU  ',ZSUR3,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in writing prof RESERV.EAU in file3',IREP
  STOP

ENDIF

!      end of reading and writing of second prof field

!      3    reading and writing of third prof field

ZSUR1=0._JPRB
ZSUR3=0._JPRB

CALL FACILE(IREP,11,'PROF',1,'RESERV.GLACE',ZSUR1,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in reading prof RESERV.GLACE in file1',IREP
  STOP

ENDIF

DO I=1,INRGRPOI

  Z_WPI=ZSUR1(I)

  Z_WPI=MAX(0._JPRB,MIN(Z_WPI,Z_GWPIMX,Z_WPMX(I))-Z_WP(I))
                
  IF (NINT(ZIVEG(I))==I_NTVMER) THEN
                     
    Z_WPI=0._JPRB
                    
  ELSEIF (NINT(ZIVEG(I))==I_NTVGLA) THEN
                 
    Z_WPI=MIN(Z_GWPIMX,Z_WPMX(I))
                    
  ENDIF

  ZSUR3(I)=Z_WPI 

ENDDO

CALL FAIENC(IREP,13,'PROF',1,'RESERV.GLACE',ZSUR3,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in writing prof RESERV.GLACE in file3',IREP
  STOP

ENDIF

!      end of reading and writing of third prof field
!      4 1   reading and writing of first surface field

ZSUR1=0._JPRB
ZSUR3=0._JPRB

CALL FACILE(IREP,11,'SURF',1,'RESERV.NEIGE',ZSUR1,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in reading surf RESERV.NEIGE in file1',IREP
  STOP

ENDIF

DO I=1,INRGRPOI

  ZSUR3(I)=MAX(0._JPRB,ZSUR1(I)) 

ENDDO

CALL FAIENC(IREP,13,'SURF',1,'RESERV.NEIGE',ZSUR3,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in reading surf RESERV.NEIGE in file3',IREP
  STOP

ENDIF

ZSUR1=0._JPRB
ZSUR3=0._JPRB

CALL FACILE(IREP,11,'SURF',1,'RESERV.EAU  ',ZSUR1,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in reading surf RESERV.EAU in file1',IREP
  STOP

ENDIF

DO I=1,INRGRPOI

  Z_WS(I)=ZSUR1(I)
  Z_WSMX(I)=Z_GCONV*Z_RD1*Z_WSAT(I)

  IF (NINT(ZIVEG(I))==I_NTVMER) THEN
                     
    Z_WS(I)=Z_WSMX(I)
                    
  ELSEIF (NINT(ZIVEG(I))==I_NTVGLA) THEN
                 
    Z_WS(I)=0._JPRB
                    
  ENDIF

  ZSUR3(I)=Z_WS(I)

ENDDO

CALL FAIENC(IREP,13,'SURF',1,'RESERV.EAU  ',ZSUR3,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in writing surf RESERV.EAU in file3',IREP
  STOP

ENDIF

!      end of reading and writing of third surface field
!      7 4   reading and writing of fourth surface field

ZSUR1=0._JPRB
ZSUR3=0._JPRB

CALL FACILE(IREP,11,'SURF',1,'RESERV.INTER',ZSUR1,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in reading surf RESERV.INTER in file1',IREP
  STOP

ENDIF

!       Wl=0 up to now set in 927            

CALL FAIENC(IREP,13,'SURF',1,'RESERV.INTER',ZSUR3,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in reading surf RESERV.INTER in file3',IREP
  STOP

ENDIF

!      end of reading and writing of fourth surface field
!      8 5   reading and writing of fifth surface field

ZSUR1=0._JPRB
ZSUR3=0._JPRB

CALL FACILE(IREP,11,'SURF',1,'RESERV.GLACE',ZSUR1,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in reading surf RESERV.GLACE in file1',IREP
  STOP

ENDIF

DO I=1,INRGRPOI

  Z_WSI=ZSUR1(I)

  Z_WSI=MAX(0._JPRB,MIN(Z_WSI,Z_WSMX(I))-Z_WS(I))
                
  IF (NINT(ZIVEG(I))==I_NTVMER) THEN
                     
    Z_WSI=0._JPRB
                    
  ELSEIF (NINT(ZIVEG(I))==I_NTVGLA) THEN
                 
    Z_WSI=Z_WPMX(I)
                    
  ENDIF

  ZSUR3(I)=Z_WSI

ENDDO

CALL FAIENC(IREP,13,'SURF',1,'RESERV.GLACE',ZSUR3,.FALSE.)

IF (IREP/=0) THEN

  WRITE(*,*) 'pb in writing surf RESERV.GLACE in file3',IREP
  STOP
ENDIF

!      end of reading and writing of fifth surface field

!      close file

CALL FAIRME(IREP,11,'KEEP')
CALL FAIRME(IREP,13,'KEEP')

WRITE(*,*) 'all files are checked'
STOP

!      -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CHECK_LIMITS',1,ZHOOK_HANDLE)
END PROGRAM CHECK_LIMITS
