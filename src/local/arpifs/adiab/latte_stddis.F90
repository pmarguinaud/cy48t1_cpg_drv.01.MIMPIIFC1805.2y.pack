SUBROUTINE LATTE_STDDIS(YDGEOMETRY,YDVARS,YDCPG_SL2,KST,KPROF,PDT,LDCOMADH,LDCOMADV)

!**** *LATTE_STDDIS* input for corrected SL scheme COMAD.
!                    Computation of the stretching/shrinking functions
!                    in the fluid (based on the derivative of the velocity,
!                    cf definition of the 1D divergence/deformation)
!                    Functions are stored in PB2 to be used later under CALL_SL.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *LATTE_STDDIS(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KST      - first element of work.
!          KPROF    - depth of work.
!          PDT      - time step 
!          LDCOMADH - COMAD for horizontal SL interpolations
!          LDCOMADV - COMAD for vertical   SL interpolations
!          PGMV     - variables at time t, only wind derivatives will be used in PGMV

!        INPUT/OUTPUT:
!          PB2      - "SLBUF2" buffer.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Arpege documentation about semi-Lagrangian scheme.

!     Author.
!     -------
!        Original (S.Malardel): november 2013

!     Modifications.
!     --------------

!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE FIELD_VARIABLES_MOD,ONLY : FIELD_VARIABLES
USE CPG_TYPE_MOD, ONLY : CPG_SL2_TYPE

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(FIELD_VARIABLES) ,INTENT(INOUT)  :: YDVARS
TYPE(CPG_SL2_TYPE)    ,INTENT(INOUT)  :: YDCPG_SL2
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT 
LOGICAL           ,INTENT(IN)    :: LDCOMADH,LDCOMADV
!     -------------------------------------------------------
INTEGER(KIND=JPIM) :: JLEV

! stretching/shrinking deformation along zonal direction
REAL(KIND=JPRB)    :: ZSTDDISU(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
! stretching/shrinking deformation along meridional direction
REAL(KIND=JPRB)    :: ZSTDDISV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
! stretching/shrinking deformation along vertical direction
REAL(KIND=JPRB)    :: ZSTDDISW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
 
REAL(KIND=JPRB)    :: ZHOOK_HANDLE

!     -------------------------------------------------------

#include "gp_stddis.intfb.h"

!     -------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('LATTE_STDDIS',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV)
ASSOCIATE(NPROMA=>YDDIM%NPROMA,   NFLEVG=>YDDIMV%NFLEVG)
!     -------------------------------------------------------

!*       1.   Stretching functions computation

CALL GP_STDDIS(NPROMA,KST,KPROF,NFLEVG,PDT,&
 & LDCOMADH,LDCOMADV,&
 & YDVARS%U%DL,YDVARS%DIV%T0,&
 & ZSTDDISU,ZSTDDISV,ZSTDDISW)

!*       2.   Fill PB2 with functions

DO JLEV=1,NFLEVG
  YDCPG_SL2%STDDISU(KST:KPROF,JLEV)=ZSTDDISU(KST:KPROF,JLEV)
  YDCPG_SL2%STDDISV(KST:KPROF,JLEV)=ZSTDDISV(KST:KPROF,JLEV)
  YDCPG_SL2%STDDISW(KST:KPROF,JLEV)=ZSTDDISW(KST:KPROF,JLEV)
ENDDO

!    --------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('LATTE_STDDIS',1,ZHOOK_HANDLE)
END SUBROUTINE LATTE_STDDIS
