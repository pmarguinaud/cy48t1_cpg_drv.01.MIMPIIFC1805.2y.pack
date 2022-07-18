SUBROUTINE GPHVINT(KPROMA,KFLEV,KSTART,KPROF,PRESH,PRESF,PXYB,PWEIG)

!**** *GPHVINT* - Computes half and full level vertical 
!                interpolating weights

!     Purpose.
!     --------
!           Computes inerpolating weights at half and full model levels.

!**   Interface.
!     ----------
!        *CALL* *GPHINT(...)

!        Explicit arguments :
!        --------------------

!          KPROMA    : horizontal dimensioning                                (in)
!          KFLEV     : vertical dimensioning                                  (in)
!          KSTART    : start of work                                          (in)
!          KPROF     : depth of work                                          (in)
!          PRESH     : half level pressure                                    (in)
!          PRESF     : full level pressure                                    (in)
!          PXYB      : contains pressure depth, "delta", "alpha"              (in)
!          PWWI      : vertical interpolating weights                         (out)

!        Implicit arguments :  NONE.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      F. VOITUS (Feb 2020) after GPPRE, GPPREH, GPXYB and GPPREF.

!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK
USE YOMDYNA   , ONLY : YRDYNA
USE INTDYN_MOD, ONLY : YYTXYB

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM)         ,INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM)         ,INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM)         ,INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM)         ,INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)            ,INTENT(IN)    :: PRESH(KPROMA,0:KFLEV)
REAL(KIND=JPRB)            ,INTENT(IN)    :: PRESF(KPROMA,KFLEV)
REAL(KIND=JPRB)            ,INTENT(IN)    :: PXYB(KPROMA,KFLEV,YYTXYB%NDIM)
REAL(KIND=JPRB)            ,INTENT(OUT)   :: PWEIG(KPROMA,KFLEV,3)


!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JLON

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPHVINT',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!*    COMPUTES VERTICAL INTERPOLATING WEIGHTS
!     ---------------------------------------


IF(.NOT.YRDYNA%LVEREGINT) THEN

  DO JLEV=1,KFLEV
    DO JLON=KSTART,KPROF
      PWEIG(JLON,JLEV,1)=PXYB(JLON,JLEV,YYTXYB%M_ALPH)/PXYB(JLON,JLEV,YYTXYB%M_LNPR)
      PWEIG(JLON,JLEV,2)=1._JPRB-(PXYB(JLON,JLEV,YYTXYB%M_ALPH)/PXYB(JLON,JLEV,YYTXYB%M_LNPR))
    ENDDO
  ENDDO
        
  DO JLEV=1,KFLEV-1
    DO JLON=KSTART,KPROF
      PWEIG(JLON,JLEV,3)=(PXYB(JLON,JLEV+1,YYTXYB%M_LNPR)-PXYB(JLON,JLEV+1,YYTXYB%M_ALPH))&
       & /(PXYB(JLON,JLEV+1,YYTXYB%M_LNPR)-PXYB(JLON,JLEV+1,YYTXYB%M_ALPH)+PXYB(JLON,JLEV,YYTXYB%M_ALPH))
    ENDDO
  ENDDO
  DO JLON=KSTART,KPROF
    !* bottom half level:
    PWEIG(JLON,KFLEV,3)=0.0_JPRB
  ENDDO

ELSE

  DO JLEV=2,KFLEV-1
    DO JLON=KSTART,KPROF
      PWEIG(JLON,JLEV,1)=0.5_JPRB*(PRESF(JLON,JLEV)-PRESF(JLON,JLEV-1))/PXYB(JLON,JLEV,YYTXYB%M_DELP)
      PWEIG(JLON,JLEV,2)=0.5_JPRB*(PRESF(JLON,JLEV+1)-PRESF(JLON,JLEV))/PXYB(JLON,JLEV,YYTXYB%M_DELP)  
    ENDDO
  ENDDO
  DO JLEV=1,KFLEV-1
    DO JLON=KSTART,KPROF
      PWEIG(JLON,JLEV,3)=0.5_JPRB
    ENDDO
  ENDDO 
  DO JLON=KSTART,KPROF
    !* uppermost full level: 
    PWEIG(JLON,1,1)=(PRESF(JLON,1)-PRESH(JLON,0))/PXYB(JLON,1,YYTXYB%M_DELP)
    PWEIG(JLON,1,2)=0.5_JPRB*(PRESF(JLON,2)-PRESF(JLON,1))/PXYB(JLON,1,YYTXYB%M_DELP)
    !* lowest full level: 
    PWEIG(JLON,KFLEV,1)=0.5_JPRB*(PRESF(JLON,KFLEV)-PRESF(JLON,KFLEV-1))/PXYB(JLON,KFLEV,YYTXYB%M_DELP)
    PWEIG(JLON,KFLEV,2)=(PRESH(JLON,KFLEV)-PRESF(JLON,KFLEV))/PXYB(JLON,KFLEV,YYTXYB%M_DELP)
    !* bottom half level:
    PWEIG(JLON,KFLEV,3)=0.0_JPRB 
  ENDDO

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPHVINT',1,ZHOOK_HANDLE)
END SUBROUTINE GPHVINT
