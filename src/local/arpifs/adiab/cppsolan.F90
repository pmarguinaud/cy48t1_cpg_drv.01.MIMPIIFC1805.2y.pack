SUBROUTINE CPPSOLAN(YDCST,KLON,KST,KEND,PGEMU,PGELAM,PMU0)

!**** *CPPSOLAN* - Computation of mean solar angle

!     Purpose. Computation of mean solar angle for obtaining 
!     -------- the radiation coeeficients for simpl.rad.scheme

!**   Interface.
!     ----------
!       *CALL* *CPPSOLAN(...)*

!     Explicit arguments: KST  - start address   (input)
!     ------------------- KEND - end address     (input)

!     Implicit arguments: None
!     ------------------

!     Method:
!     -------

!     Author:
!     -------
!      M. Janiskova - 96-05

!     Modifications:
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      P. Bechtold   14/05/2012 replace 86400 by RDAY
!      K. Yessad (July 2014): Move some variables.
!   ----------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMDIM   , ONLY : TDIM
USE YOMCST   , ONLY : TCST

USE YOMRIP0  , ONLY : NSSSSS, RTIMST
USE YOMCST   , ONLY : RPI, RDAY, REA, REPSM
USE YOMDYNCORE,ONLY : LAPE

!   ----------------------------------------------------------------

IMPLICIT NONE

TYPE(TCST)        ,INTENT(IN)    :: YDCST
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMU0(KLON) 

!   ----------------------------------------------------------------

INTEGER(KIND=JPIM) :: IPR, ITIME, IZT, JROF, JSTEP, I_NFIN

REAL(KIND=JPRB) :: ZMU0, ZMXX, ZNUM, ZRCODEC, ZRCOVSR, ZRDECLI,&
 & ZREQTIM, ZRHGMT, ZRSIDEC, ZRSIVSR, ZRSOVR,&
 & ZRSTATI, ZRTIMTR, ZRWSOVR, ZTETA, ZTSTEP  
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!   ----------------------------------------------------------------

#include "fctast.func.h"

!   ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CPPSOLAN',0,ZHOOK_HANDLE)

!   ----------------------------------------------------------------

ZTSTEP=900._JPRB
I_NFIN=NINT(RDAY/ZTSTEP)

DO JROF=KST,KEND
  ZMU0=0.0_JPRB
  ZNUM=0.0_JPRB

  DO JSTEP=1,I_NFIN
    ITIME=NINT(ZTSTEP)
    IZT=ITIME*JSTEP
    ZRSTATI=REAL(IZT,JPRB)
    ZRTIMTR=RTIMST+ZRSTATI
    IPR=0
    ZRHGMT=REAL(MOD(NINT(ZRSTATI)+NSSSSS,NINT(RDAY)),JPRB)

    ZTETA=RTETA(ZRTIMTR)
    IF( LAPE ) THEN
      ZRDECLI=RDSAQUA(ZTETA)
      ZREQTIM=RETAQUA(ZTETA)
    ELSE
      ZRDECLI=RDS(ZTETA)
      ZREQTIM=RET(ZTETA)
    ENDIF
    ZRSOVR =ZREQTIM+ZRHGMT
    ZRWSOVR=ZRSOVR*2.0_JPRB*RPI/RDAY

    ZRCODEC=COS(ZRDECLI)
    ZRSIDEC=SIN(ZRDECLI)

    ZRCOVSR=COS(ZRWSOVR)
    ZRSIVSR=SIN(ZRWSOVR)

!*       1.1   Astronomy.

    ZMXX=MAX( ZRSIDEC*PGEMU(JROF)&
     & -ZRCODEC*ZRCOVSR*SQRT(1.0_JPRB-PGEMU(JROF)**2)&
     & *COS(PGELAM(JROF))&
     & +ZRCODEC*ZRSIVSR*SQRT(1.0_JPRB-PGEMU(JROF)**2)&
     & *SIN(PGELAM(JROF))&
     & ,0.0_JPRB)  
    ZMU0=ZMU0+ZMXX
    IF (ZMXX /= 0.0_JPRB) THEN
      ZNUM=ZNUM+1
    ENDIF

  ENDDO
  PMU0(JROF)=ZMU0/MAX(1.0_JPRB,ZNUM)
ENDDO

!   ----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPPSOLAN',1,ZHOOK_HANDLE)
END SUBROUTINE CPPSOLAN
