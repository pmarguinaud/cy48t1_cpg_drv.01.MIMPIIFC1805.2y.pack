SUBROUTINE SUJBSTD(YDGEOMETRY,YDFIELDS,YD_TRAJ,CDJBTYPE,KPRT,PCOVTPS,PCOVP,PCOVD,PCOVQ,&
 & PCOVO3,KLEV,KNUMB,YD_JB_STRUCT)  

!**** *SUJBSTD* - Prepare gdpt standard errors for background constraint Jb

!     Purpose. Prepare the error covariance model.
!     --------

!     Interface. CALL SUJBSTD
!     ----------

!     Explicit arguments : In : CDJBTYPE = config string for Jb
!     -------------------- In : KPRT = print level
!                          In : PCOVx   = covariances for parameter x
!                               (only used by CDJBTYPE=NONSEP93)

!     Implicit arguments : common YOMJG and related files.
!     --------------------

!     Method.
!     -------
!         Compute total variance vertical profiles (sometimes using
!          an f-plane balance assumption).
!         Empirically reduce NMC-method variances to agree with 6-h
!          forecast std devs.
!         Compute length-scales as a parameter in the correlation model:
!          f(r)=A.exp(-r^2/(2a^2)) --> a=\sqrt{ \frac{-2 f(0)}{f''(0)} }
!          There is a 2 because we are in 2-D.

!     Externals: see below.
!     ----------

!     Author : Francois Bouttier *ECMWF*  96-02-07
!     --------

!     Modifications :
!     ---------------
!      01-08-07 R. El Khatib : Pruning options 
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib 22-Mar-2012 Fix uninitialized variables
!      K. Yessad (Nov 2012): call GPHPRE.
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.

!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE FIELDS_MOD , ONLY : FIELDS
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMMP0   , ONLY : NPROC
USE YOMJG    , ONLY : TYPE_JB_STRUCT
USE YOM_GRIB_CODES  , ONLY : NGRBVO, NGRBD, NGRBT, NGRBQ, NGRBO3, NGRBLNSP
USE YOMVERT  , ONLY : VP00
USE YOMLUN   , ONLY : NULOUT
USE MPL_MODULE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(FIELDS)        ,INTENT(INOUT) :: YDFIELDS
TYPE(FIELDS)        ,INTENT(INOUT) :: YD_TRAJ
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KNUMB 
CHARACTER(LEN=10)   ,INTENT(IN)    :: CDJBTYPE 
INTEGER(KIND=JPIM)  ,INTENT(IN)    :: KPRT 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PCOVTPS(KLEV+1,KLEV+1,0:KNUMB-1) 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PCOVP(KLEV,KLEV,0:KNUMB-1) 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PCOVD(KLEV,KLEV,0:KNUMB-1) 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PCOVQ(KLEV,KLEV,0:KNUMB-1) 
REAL(KIND=JPRB)     ,INTENT(IN)    :: PCOVO3(KLEV,KLEV,0:KNUMB-1) 
TYPE(TYPE_JB_STRUCT),INTENT(INOUT) :: YD_JB_STRUCT
!     ------------------------------------------------------------------
REAL(KIND=JPRB) :: ZVARP(YDGEOMETRY%YRDIMV%NFLEVG),ZVARD(YDGEOMETRY%YRDIMV%NFLEVG),ZVARQ(YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZVARO3(YDGEOMETRY%YRDIMV%NFLEVG),&
 & ZVART(YDGEOMETRY%YRDIMV%NFLEVG),ZVARTU(YDGEOMETRY%YRDIMV%NFLEVG),ZVARSP,ZVARSPU  
REAL(KIND=JPRB) :: ZVARVO(YDGEOMETRY%YRDIMV%NFLEVG),ZVARU(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPRESF(YDGEOMETRY%YRDIMV%NFLEVG+1),ZPRESH(0:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZFUDGE, ZFCEMN_VO, ZFCEMN_D, ZFCEMN_TU, ZFCEMN_O3, ZFCEMN_Q
INTEGER(KIND=JPIM), PARAMETER :: JPKZ = KIND(ZFUDGE)

INTEGER(KIND=JPIM) :: IL, ILV950, J1, JJ, JLEV, IPTJBVO, IPTJBD,&
                    & IPTJBTU, IPTJBQ, IPTJBO3, IPTJBPSU,&
                    & JFIELD
LOGICAL :: LLO3JB
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "gphpre.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUJBSTD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, YDLAP=>YDGEOMETRY%YRLAP, YDVAB=>YDGEOMETRY%YRVAB&
& )

ASSOCIATE(NSMAX=>YDDIM%NSMAX,   NFLEVG=>YDDIMV%NFLEVG,   RLAPIN=>YDLAP%RLAPIN)
!     ------------------------------------------------------------------

!        0. Preliminary checks

LLO3JB=ANY(YD_JB_STRUCT%SPJB_VARS_INFO(:)%IGRIBCODE==NGRBO3)

WRITE(NULOUT,*) 'SUJBSTD: Jb std errors, config=',CDJBTYPE
IF ((CDJBTYPE=='STABAL96').OR.(&
   & CDJBTYPE=='TOTENRGY') ) THEN  
  IF((KLEV /= NFLEVG).OR.(KNUMB-1 /= NSMAX))&
   & CALL ABOR1('BAD INTERFACE TO SUJBCOR')  
ELSE
  CALL ABOR1('SUJBCOR: Unknown Jb std errors config')
ENDIF

!-----------------------------------------------------------------------

!        2. Standard errors, STABAL96-style (John Derber's balance)

IF ((CDJBTYPE=='STABAL96').OR.(CDJBTYPE=='TOTENRGY')) THEN
  WRITE(NULOUT,*) ' Setup vertical profiles from covariances'

!         2.1 Total variances implied by the 24/48 NMC method

!        Note that total T,ps,P and D are NOT Jb variables.
!        However, total D can be safely identified with  D.
  DO J1=1,NFLEVG
    ZVARVO(J1) = SUM(PCOVP(J1,J1,:))
    ZVARD(J1)  = SUM(PCOVD(J1,J1,:))
    ZVARTU(J1) = SUM(PCOVTPS(J1,J1,:))
    ZVARQ(J1)  = SUM(PCOVQ(J1,J1,:))
  ENDDO
  IF(LLO3JB)THEN
    DO J1=1,NFLEVG
      ZVARO3(J1) = SUM(PCOVO3(J1,J1,:))
    ENDDO
  ENDIF
  ZVARSPU=SUM(PCOVTPS(NFLEVG+1,NFLEVG+1,:))

!         2.2 Estimate variances for missing fields

!        The following massaging of the sigma-b profiles has no impact
!       on the Jb term, but it will affect the threshold for screening
!       and VarQC in the next analysis cycle. Weaknesses in the currently
!       assumed sigma-b's are compensated by the tuning of rejection
!       thresholds.  In the longer term, genuine 6-h sigma-b profiles shall 
!       be read from separate files, produced with the help of 
!       Hollingsworth/Lonnberg diagnostics. In the meantime, do not change
!       these values without retuning the screening/QC thresholds
!       accordingly. /FB

  WRITE(NULOUT,*) ' Total (T,ps) profile is given by',' scaling (T,ps)u'
  ZVART(:)=ZVARTU(:)*4._JPRB
  ZVARSP=ZVARSPU*4._JPRB

  WRITE(NULOUT,*) ' Wind profile is given by ',&
   & ' inverse Laplacian of vorticity'  
!       The factor 0.5 is necessary in 2-D.
  DO J1=1,NFLEVG
    ZVARU(J1)  = SUM(-0.5_JPRB*RLAPIN(0:NSMAX)*PCOVP(J1,J1,0:NSMAX))
  ENDDO

  WRITE(NULOUT,*) ' Total mass profile is given by ',&
   & ' scaled streamfunction above 950hPa'  
!       Find lowest level above 950hPa in std atm (alt 630m)
  ZPRESH(NFLEVG)=VP00
  CALL GPHPRE(1,NFLEVG,1,1,YDVAB,ZPRESH,PRESF=ZPRESF)
  ILV950=-HUGE(ILV950)
  DO JLEV=1,NFLEVG
    IF(ZPRESF(JLEV)<95000._JPRB) ILV950=JLEV
  ENDDO
  WRITE(NULOUT,*) '   profile is extrapolated below level ',ILV950
  DO J1=1,NFLEVG
    IL=MIN(J1,ILV950)
!         Do not try to understand this formula
    ZFUDGE=(1.6_JPRB+REAL(IL,JPKZ)/REAL(NFLEVG,JPKZ)*(1.3_JPRB-1.6_JPRB))*&
     & 0.00005_JPRB  
    ZVARP(J1)  = SUM( (RLAPIN(0:NSMAX)**2)*PCOVP(IL,IL,0:NSMAX) )&
     & * (ZFUDGE**2)  
  ENDDO

!         2.3 Load vertical profile arrays:
!             REDNMC is an empirical reduction to get the NMC-method
!             variances to agree with typical 6-hour forecast errors.
!             The f-plane approximation is ok for z except near the ground.

  IPTJBVO  = HUGE(IPTJBVO)
  IPTJBD   = HUGE(IPTJBD)
  IPTJBTU  = HUGE(IPTJBTU)
  IPTJBQ   = HUGE(IPTJBQ)
  IPTJBO3  = HUGE(IPTJBO3)
  IPTJBPSU = HUGE(IPTJBPSU)
  DO JFIELD=1,YD_JB_STRUCT%CONFIG%N_SPJB_VARS
    IF (YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBVO) THEN
      IPTJBVO = YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBD) THEN
      IPTJBD  = YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBT) THEN
      IPTJBTU = YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBQ) THEN
      IPTJBQ  = YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBO3) THEN
      IPTJBO3 = YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSEIF (YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE==NGRBLNSP) THEN
      IPTJBPSU= YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IPTJB
    ELSE
      WRITE (NULOUT,*) 'DON''T KNOW WHAT TO DO WITH GRIB CODE ',&
                      & YD_JB_STRUCT%SPJB_VARS_INFO(JFIELD)%IGRIBCODE
    ENDIF
  ENDDO

  DO JJ=1,NFLEVG
    IF (IPTJBVO/=HUGE(IPTJBVO)) THEN
      YD_JB_STRUCT%JB_DATA%FCEMN3D(JJ,IPTJBVO )=&
                     & YD_JB_STRUCT%CONFIG%REDNMC*SQRT(ZVARVO(JJ))
    ENDIF
    IF (IPTJBD/=HUGE(IPTJBD)) THEN
      YD_JB_STRUCT%JB_DATA%FCEMN3D(JJ,IPTJBD  )=&
                     & YD_JB_STRUCT%CONFIG%REDNMC*SQRT(ZVARD(JJ))
    ENDIF
    YD_JB_STRUCT%FCEMN%U(JJ) =YD_JB_STRUCT%CONFIG%REDNMC*SQRT(ZVARU(JJ))
    YD_JB_STRUCT%FCEMN%Z(JJ) =YD_JB_STRUCT%CONFIG%REDNMC*SQRT(ZVARP(JJ))
    YD_JB_STRUCT%FCEMN%T(JJ) =YD_JB_STRUCT%CONFIG%REDNMC*SQRT(ZVART(JJ))
    IF (IPTJBTU/=HUGE(IPTJBTU)) THEN
      YD_JB_STRUCT%JB_DATA%FCEMN3D(JJ,IPTJBTU )=&
                     & YD_JB_STRUCT%CONFIG%REDNMC*SQRT(ZVARTU(JJ))
    ENDIF
    IF (IPTJBQ/=HUGE(IPTJBQ)) THEN
      YD_JB_STRUCT%JB_DATA%FCEMN3D(JJ,IPTJBQ  )= SQRT(ZVARQ(JJ)) ! Q not rescaled
    ENDIF
    IF (IPTJBO3/=HUGE(IPTJBO3)) THEN
      YD_JB_STRUCT%JB_DATA%FCEMN3D(JJ,IPTJBO3)=&
                     & YD_JB_STRUCT%CONFIG%REDNMC*SQRT(ZVARO3(JJ))
    ENDIF
  ENDDO
  IF (IPTJBPSU/=HUGE(IPTJBPSU)) THEN
    YD_JB_STRUCT%JB_DATA%FCEMN2D(IPTJBPSU)= YD_JB_STRUCT%CONFIG%REDNMC*SQRT(ZVARSPU)
  ENDIF
  YD_JB_STRUCT%FCEMN%PS = YD_JB_STRUCT%CONFIG%REDNMC*SQRT(ZVARSP)

  IF(KPRT>=1)THEN
    WRITE(NULOUT,*) 'Jb preliminary vertical stdev profiles:'
    WRITE(NULOUT,*) '(-999.0 indicates variable is not in Jb)'
    IF (IPTJBPSU/=HUGE(IPTJBPSU)) THEN
      WRITE(NULOUT,*) ' unbal ps : ',&
             & YD_JB_STRUCT%JB_DATA%FCEMN2D(IPTJBPSU)*VP00
    ELSE
      WRITE(NULOUT,*) ' unbal ps : ',-999.0_JPRB
    ENDIF
    WRITE(NULOUT,*) ' total ps : ',YD_JB_STRUCT%FCEMN%PS *VP00
    WRITE(NULOUT,*) '       Vor      unbal Div  '//&
     & '    Wind        mass    '//&
     & '     T       unbal T    '//&
     & '     Q       unbal O3'  
    DO JJ=1,NFLEVG
      IF (IPTJBVO/=HUGE(IPTJBVO)) THEN
        ZFCEMN_VO = YD_JB_STRUCT%JB_DATA%FCEMN3D(JJ,IPTJBVO)
      ELSE
        ZFCEMN_VO = -999.0_JPRB
      ENDIF
      IF (IPTJBD/=HUGE(IPTJBD)) THEN
        ZFCEMN_D = YD_JB_STRUCT%JB_DATA%FCEMN3D(JJ,IPTJBD)
      ELSE
        ZFCEMN_D = -999.0_JPRB
      ENDIF
      IF (IPTJBTU/=HUGE(IPTJBTU)) THEN
        ZFCEMN_TU = YD_JB_STRUCT%JB_DATA%FCEMN3D(JJ,IPTJBTU)
      ELSE
        ZFCEMN_TU = -999.0_JPRB
      ENDIF
      IF (IPTJBQ/=HUGE(IPTJBQ)) THEN
        ZFCEMN_Q = YD_JB_STRUCT%JB_DATA%FCEMN3D(JJ,IPTJBQ)
      ELSE
        ZFCEMN_Q = -999.0_JPRB
      ENDIF
      IF (IPTJBO3/=HUGE(IPTJBO3)) THEN
        ZFCEMN_O3 = YD_JB_STRUCT%JB_DATA%FCEMN3D(JJ,IPTJBO3)
      ELSE
        ZFCEMN_O3 = -999.0_JPRB
      ENDIF
      WRITE(NULOUT,'(I3,8(1X,G11.4))') JJ,ZFCEMN_VO,&
       & ZFCEMN_D,YD_JB_STRUCT%FCEMN%U(JJ) ,YD_JB_STRUCT%FCEMN%Z(JJ),&
       & YD_JB_STRUCT%FCEMN%T(JJ),ZFCEMN_TU,ZFCEMN_Q,&
       & ZFCEMN_O3
    ENDDO
  ENDIF
ENDIF

!-----------------------------------------------------------------------

!        9. Initialize gridpoint sigma-b's from external file

WRITE(NULOUT,*) '2.Read sigma-b fields from GRIB gdpt file'
CALL FLUSH(NULOUT)
IF (NPROC>1) THEN
  CALL GSTATS(721,0)
  CALL MPL_BARRIER(CDSTRING='SUJBSTD:')
  CALL GSTATS(721,1)
ENDIF

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUJBSTD',1,ZHOOK_HANDLE)
END SUBROUTINE SUJBSTD
