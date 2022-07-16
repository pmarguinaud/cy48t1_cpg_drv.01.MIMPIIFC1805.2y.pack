SUBROUTINE SUVWRK(KULOUT)

!**** *SURWRK*   - Routine to initialize path names for trajectory
!                -                              and for cma files
!     Purpose.
!     --------
!           Initialize path names
!       update 2013: initialization of one parameter related to SL traj. handling
!**   Interface.
!     ----------
!        *CALL* *SUVWRK(...) (from SU0YOMA)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMVWRK

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Jean-Noel Thepaut *ECMWF*  92-07-13

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      F. Vana       23-Oct-2013 removing useless code, adding security
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMLUN   , ONLY : NULNAM
USE YOMCT0   , ONLY : LR2D
USE YOMVWRK  , ONLY : NTRSLTYPE
USE YOMDYNA  , ONLY : LELTRA

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KULOUT 

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

#include "posnam.intfb.h"

#include "namvwrk.nam.h"

!      ----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVWRK',0,ZHOOK_HANDLE)
!      ----------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!              -------------------

!        1.1 Set implicit default value

NTRSLTYPE=2

!      ----------------------------------------------------------------
!*       2.    Modifies default values.
!              ------------------------

CALL POSNAM(NULNAM,'NAMVWRK')
READ(NULNAM,NAMVWRK)
!      ----------------------------------------------------------------
!*       3.    Check the namelist value
!              ------------------------

IF ((NTRSLTYPE==0).AND.((.NOT.LR2D).OR.(.NOT.LELTRA))) &
 & CALL ABOR1('NTRSLTYPE=0 could only be used in 2D model with LELTRA.')
 
IF ((NTRSLTYPE==1).AND.(.NOT.LR2D)) &
 & CALL ABOR1('NTRSLTYPE=1 could only be used in 2D model.')

IF ((NTRSLTYPE < 0).OR.(NTRSLTYPE > 2)) &
 & CALL ABOR1('NTRSLTYPE value is only allowed to be 0, 1 or 2.')

!      ----------------------------------------------------------------
!*       4.    Print final value.
!              -------------------

WRITE(UNIT=KULOUT,FMT='('' COMMON YOMVWRK '')')
WRITE(UNIT=KULOUT,FMT='('' NTRSLTYPE ='',I2)') NTRSLTYPE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVWRK',1,ZHOOK_HANDLE)
END SUBROUTINE SUVWRK
