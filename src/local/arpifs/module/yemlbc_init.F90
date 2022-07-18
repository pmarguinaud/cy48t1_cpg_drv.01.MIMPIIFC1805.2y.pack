MODULE YEMLBC_INIT

! Purpose :
! -------
!    Forcing a LAM model by another model: part 0A
!    - forcing by lateral boundary conditions
!    - pressure tendency coupling
!    - spectral nudging

! Interface :
! ---------
!    Empty.

! External :
! --------
!    None.

! Method :
! ------
!    See Documentation.

! Reference :
! ---------

! Author :
! ------
!    K. YESSAD (CNRM/GMAP) after YOMCT0, YEMCT0, YEMSPCPL, SUCT0, SUECT0.

! Modifications :
! -------------
! Original : December 2010
! H. Dhouioui (Sep 2017) renamed from yemgeolbc.F90 
! P. Smolikova (Sep 2020): NBZONU for LUNBC option.
!-----------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMLUN   , ONLY : NULNAM   ,NULOUT
USE YOMCT0   , ONLY : NCONF    ,LNHDYN   ,LELAM
USE YOMDYNA  , ONLY : YRDYNA


IMPLICIT NONE
SAVE

!=============================================================================

!      1.    TYPE DEFINITION
!            ---------------

! none for the time being

!=============================================================================

!      2.    DECLARATIONS
!            ------------

!      2.1   Coupling of surface pressure tendency

! LTENC   : TRUE if tendency coupling of surface pressure is switched on
! LALLTC  : used together with LTENC when LTENC=.T.
!           - no meaning for quadratic tendency coupling, where just t1 coupling
!             is applied at every NEFRCL time step
!           - for lin. tendency coupling:
!             TRUE if tendency coupling of surf. pres. at every step
!             FALSE if tend. coupl., except at every NEFRCL time steps
!             when just t1 coupling

LOGICAL :: LTENC
LOGICAL :: LALLTC
LOGICAL :: LRFIRST   ! Force reading of first coupling file (usually, it is the
                     ! same as the initial conditions file)

!      2.2   Lateral forcing

!   NBICOU  : controls coupling of wind components (GMV)
!   NBICOT  : controls coupling of temperature (GMV)
!   NBICPD  : controls coupling of pressure departure variable (GMV)
!   NBICVD  : controls coupling of vertical divergence variable (GMV)
!   NBICNHX : controls coupling of "NHX" term (GMV)
!   NBICOP  : controls coupling of surface pressure (GMVS)
!   Possible value for the NBIC[X] variables:
!   * 0: no coupling
!   * 1: default coupling
!   * 2: specific coupling function

!   NECRIPL : controls timelevel of coupling
!   * 0: coupling at t
!   * 1: coupling at t+dt
!   * 2: coupling at t and t+dt
!   LQCPL   : if T (resp. F), quadratic (resp. linear) temporal interpolation
!   LCCPL   : if T (resp. F), cubic (resp. linear) temporal interpolation
!   NECOTL  : Controls the coupling in the tangent linear model:
!             0   - no coupling
!             < 0 - coupling with the Davies relaxation
!             -1  - the linear time interpolation is switched on
!             > 0 - coupling other than Davies relaxation (to be implemented)
!   NECOAD  : cf. NECOTL but for the adjoint model.
!   LE0COTA : TRUE if the relaxation in the I+E zone is towards
!             nullified boundary conditions (assumed exact);
!             FALSE if the relaxation is towards a predefined forcing.
!             LE0COTA has the same meaning in the TL and AD models
!             (respectively, relaxation to 0 perturbation or sensitivity)
!   LEREADINI: TRUE if initial historical file has to be read
!              (used for E501, E801 - not a namelist parameter)
!   N1LSG   : if =1, the gradient with respect to the large scale
!             coupling data has to be written into a file.
!   NFRLSG  : frequency of writing large scale gradients (time or steps)
!   NLSGTS  : array containing large scale gradients timesteps
!             * if NLSGTS(0)=0 action if MOD(JSTEP,NFRLSG)=0
!             * if NLSGTS(0)>0 NLSGTS(0) significant numbers in NLSGTS
!               are then considered and action for JSTEP=NLSGTS(.)*NFRLSG
!             * if NLSGTS(0)<0 action for JSTEP=(NLSGTS(.)/DELTAT)*NFRLSG
!   LRDLSG  : switch for using boundary data perturbation in conf 1
!             for sensitivity forecast run.
!   JPALFNM : Dimension for reading alpha function parameters

INTEGER(KIND=JPIM) :: NBICOU
INTEGER(KIND=JPIM) :: NBICOT
INTEGER(KIND=JPIM) :: NBICPD
INTEGER(KIND=JPIM) :: NBICVD
INTEGER(KIND=JPIM) :: NBICNHX
INTEGER(KIND=JPIM) :: NBICOP
INTEGER(KIND=JPIM) :: NECRIPL
LOGICAL :: LQCPL
LOGICAL :: LCCPL
INTEGER(KIND=JPIM) :: NECOTL
INTEGER(KIND=JPIM) :: NECOAD
LOGICAL :: LE0COTA
LOGICAL :: LEREADINI
INTEGER(KIND=JPIM), PARAMETER :: JPLSGT=20
INTEGER(KIND=JPIM) :: N1LSG
INTEGER(KIND=JPIM) :: NFRLSG
INTEGER(KIND=JPIM) :: NLSGTS(0:JPLSGT)
LOGICAL :: LRDLSG
INTEGER(KIND=JPIM), PARAMETER :: JPALFNM=31

!      2.3   Spectral nudging

!   LESPCPL  : control of spectral nudging

LOGICAL :: LESPCPL
LOGICAL :: LSPTENC

!      2.4   Upper nesting boundary conditions

LOGICAL :: LUNBC               ! controls upper nesting boundary conditions

!=============================================================================

CONTAINS

!      3.    SET-UP

SUBROUTINE SUELBC_INIT

!--------------------------------------------------------------------------
! Sets-up part 0A of forcing a LAM model by another model
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!--------------------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"

#include "nemelbc0a.nam.h"

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YEMLBC_INIT:SUELBC_INIT',0,ZHOOK_HANDLE)
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------

!      Part A: default values.

! * Tendency coupling
LTENC=.FALSE.
LALLTC=.FALSE.
LRFIRST=.TRUE.

! * Lateral forcing
IF (LELAM) THEN
  NBICOU=1
  NBICOT=1
  IF (LNHDYN) THEN
    NBICPD=1
    NBICVD=1
    IF(YRDYNA%LNHX) THEN
      NBICNHX=1
    ELSE
      NBICNHX=0
    ENDIF
  ELSE
    NBICPD=0
    NBICVD=0
    NBICNHX=0
  ENDIF
  NBICOP=1
  NECRIPL=1
  LQCPL=.FALSE.
  LCCPL=.FALSE.
  NECOAD=0
  NECOTL=0
  LE0COTA=.FALSE.
  LEREADINI=.TRUE.
  N1LSG=0
  NFRLSG=1
  DO J=0,JPLSGT
    NLSGTS(J)=0
  ENDDO
  LRDLSG=.FALSE.
ELSE
  NBICOU=0
  NBICOT=0
  NBICOP=0
  NBICPD=0
  NBICVD=0
  NBICNHX=0
  NECRIPL=1
  LQCPL=.FALSE.
  LCCPL=.FALSE.
  NECOAD=0
  NECOTL=0
  LE0COTA=.FALSE.
  LEREADINI=.FALSE.
  N1LSG=0
  NFRLSG=1
  DO J=0,JPLSGT
    NLSGTS(J)=0
  ENDDO
  LRDLSG=.FALSE.
ENDIF

! * Spectral nudging
LESPCPL=.FALSE.
LSPTENC=.FALSE.

! * Upper nesting boundary conditions
LUNBC=.FALSE.

!--------------------------------------------------------------------------

!      Part B: namelist reading.

IF (LELAM) THEN
  CALL POSNAM(NULNAM,'NEMELBC0A')
  READ(NULNAM,NEMELBC0A)
ENDIF

!--------------------------------------------------------------------------

!      Part C: checkings.

! * checkings on tendency coupling
IF (.NOT.LTENC.AND.LALLTC) THEN
  CALL ABOR1('SUELBC_INIT: ABOR1 CALLED: .NOT.LTENC.AND.LALLTC')
ENDIF
IF (LTENC.AND.LRDLSG) THEN
  CALL ABOR1('SUELBC_INIT: ABOR1 CALLED: LTENC.AND.LRDLSG')
ENDIF

! * checkings on lateral boundary forcing
IF (NECRIPL /= 0.AND.NECRIPL /= 1.AND.NECRIPL /= 2) THEN
  CALL ABOR1('SUELBC_INIT: IMPROPER VALUE FOR NECRIPL')
ENDIF
IF (LQCPL.AND.NCONF == 701) THEN
  CALL ABOR1('SUELBC_INIT: LQCPL=.T. NOT PREPARED FOR NCONF=701')
ENDIF
IF (NECRIPL /= 0.AND.( NBICOU /= NBICOT.OR.NBICOU /= NBICOP ) ) THEN
  WRITE(NULOUT,*) 'T1 COUPLING DOESNT LET COUPLING FUNCTIONS DIFFER'
  CALL ABOR1 ('SUELBC_INIT: NECRIPL /= 0 and improper values of some NBIC..')
ENDIF

!--------------------------------------------------------------------------

!      Part D: printings.

IF (LELAM) THEN
  WRITE(NULOUT,*) ' --- PRINTINGS IN SUELBC_INIT --- '
  WRITE(NULOUT,*) ' LTENC = ',LTENC,' LALLTC = ',LALLTC
  WRITE(NULOUT,*) ' LSPTENC = ',LSPTENC
  WRITE(NULOUT,*) ' NBICOU = ',NBICOU,' NBICOT = ',NBICOT,' NBICOP = ',NBICOP
  WRITE(NULOUT,*) ' NBICPD = ',NBICPD,' NBICVD = ',NBICVD,&
   & ' NBICNHX = ',NBICNHX,' LQCPL = ',LQCPL,' LCCPL = ',LCCPL
  WRITE(UNIT=NULOUT,FMT='('' NECRIPL = '',I2)') NECRIPL
  WRITE(NULOUT,*) ' NECOAD = ',NECOAD,' NECOTL = ',NECOTL,' LE0COTA = ',LE0COTA
  WRITE(NULOUT,FMT='('' NFRLSG = '',I2,'' N1LSG = '',I2)') NFRLSG,N1LSG
  WRITE(NULOUT,*) ' NLSGTS =  ',NLSGTS(0),(NLSGTS(J),J=1,ABS(NLSGTS(0)))
  WRITE(NULOUT,FMT='('' LRDLSG = '',L2)') LRDLSG
  WRITE(NULOUT,*) 'LESPCPL = ',LESPCPL
  WRITE(NULOUT,*) 'LUNBC = ',LUNBC
  WRITE(NULOUT,*) ' '
ENDIF

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YEMLBC_INIT:SUELBC_INIT',1,ZHOOK_HANDLE)
END SUBROUTINE SUELBC_INIT

!=============================================================================

END MODULE YEMLBC_INIT
