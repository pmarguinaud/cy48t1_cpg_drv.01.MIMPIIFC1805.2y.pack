SUBROUTINE GPHPRETL(KPROMA,KFLEV,KSTART,KPROF,YDVAB,PRESH,PRESH5,PXYB,PXYB5,PRESF)

!**** *GPHPRETL* - Computes half and full level pressure (TL code)
!                  Modern version of former GPPREHTL+GPXYBTL+GPPREFTL

!     Purpose.
!     --------
!           Computes pressures at half and full model levels.

!**   Interface.
!     ----------
!        *CALL* *GPHPRETL(...)

!        Explicit arguments :
!        --------------------

!          KPROMA    : horizontal dimensioning                                (in)
!          KFLEV     : vertical dimensioning                                  (in)
!          KSTART    : start of work                                          (in)
!          KPROF     : depth of work                                          (in)
!          YDVAB     : contains information about hybrid vertical coordinate  (in)
!          PRESH     : half level pressure                                    (inout)
!          PRESH5    : half level pressure (trajectory)                       (in)
!          PXYB      : contains pressure depth, "delta", "alpha"              (opt out)
!          PXYB5     : contains pressure depth, "delta", "alpha" (trajectory) (opt in)
!          PRESF     : full level pressure                                    (opt out)

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
!      K. YESSAD (Sep 2011) after GPPREHTL, GPXYBTL and GPPREFTL.

!     Modifications.
!     --------------
!   K. Yessad (Dec 2016): Prune obsolete options.
!   H Petithomme (Dec 2020): optimisation with pointer
!     ------------------------------------------------------------------

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK
USE YOMDYNA   , ONLY : YRDYNA
USE YOMCVER   , ONLY : LVERTFE
USE YOMCST    , ONLY : RD, RCVD
USE YOMVERT   , ONLY : TVAB, TOPPRES
USE INTDYN_MOD, ONLY : YYTXYB

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM)         ,INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM)         ,INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM)         ,INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM)         ,INTENT(IN)    :: KPROF 
TYPE(TVAB)                 ,INTENT(IN)    :: YDVAB
REAL(KIND=JPRB)            ,INTENT(INOUT) :: PRESH(KPROMA,0:KFLEV)
REAL(KIND=JPRB)            ,INTENT(IN)    :: PRESH5(KPROMA,0:KFLEV)
REAL(KIND=JPRB),OPTIONAL,TARGET,INTENT(OUT) :: PXYB(KPROMA,KFLEV,YYTXYB%NDIM)
REAL(KIND=JPRB),OPTIONAL   ,INTENT(IN)    :: PXYB5(KPROMA,KFLEV,YYTXYB%NDIM)
REAL(KIND=JPRB),OPTIONAL   ,INTENT(OUT)   :: PRESF(KPROMA,KFLEV)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IFIRST, JLEV, JLON
REAL(KIND=JPRB) :: ZPRESF,ZPRESFD,ZPRESF5,ZPRESF5D,ZMUL
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZXYB(:,:,:)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPHPRETL',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    PRELIMINARY CALCULATIONS
!              ------------------------

IF ((PRESENT(PXYB).AND.PRESENT(PXYB5)).OR.(PRESENT(PRESF).AND.PRESENT(PXYB5))) THEN

  ! This is introduced to get rid of the implicit
  ! assumption that the top level input for pressure is 0 hPa.
  ! The first block if is for economy (no do loop start up) and the second for safety.
  IF(PRESH5(KSTART,0) <= TOPPRES)THEN
    IFIRST=2
  ELSE
    IFIRST=1
    DO JLON=KSTART,KPROF
      IF(PRESH5(JLON,0) <= TOPPRES)THEN
        IFIRST=2
        EXIT
      ENDIF
    ENDDO
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       2.    COMPUTES HALF LEVEL PRESSURES
!              -----------------------------

! * compute upper-air PRESH from surface PRESH.
DO JLEV=0,KFLEV-1
  DO JLON=KSTART,KPROF
    PRESH(JLON,JLEV)=YDVAB%VBH(JLEV)*PRESH(JLON,KFLEV)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*       3.    COMPUTES PXYB
!              -------------

IF (PRESENT(PXYB)) THEN
  ZXYB => PXYB(:,:,:)
ELSE
  ALLOCATE(ZXYB(KPROMA,KFLEV,YYTXYB%NDIM))
ENDIF

IF ((PRESENT(PXYB).AND.PRESENT(PXYB5)).OR.(PRESENT(PRESF).AND.PRESENT(PXYB5))) THEN

  ZXYB(:,:,:)=0.0_JPRB

  IF(LVERTFE) THEN

    DO JLEV=1,KFLEV
      !DIR$ IVDEP
      !CDIR NODEP
      DO JLON=KSTART,KPROF
        ! Trajectory:
        ZPRESF5=YDVAB%VAF(JLEV)+YDVAB%VBF(JLEV)*PRESH5(JLON,KFLEV)
        ZPRESF5D=1.0_JPRB/ZPRESF5
        ! TL:
        ZXYB(JLON,JLEV,YYTXYB%M_DELP)=YDVAB%VDELB(JLEV)*PRESH(JLON,KFLEV)
        ZXYB(JLON,JLEV,YYTXYB%M_RDELP)=-ZXYB(JLON,JLEV,YYTXYB%M_DELP)/PXYB5(JLON,JLEV,YYTXYB%M_DELP)**2
        ZPRESF=YDVAB%VBF(JLEV)*PRESH(JLON,KFLEV)
        ZPRESFD=-(ZPRESF5D*ZPRESF5D)*ZPRESF
        ZXYB(JLON,JLEV,YYTXYB%M_LNPR)=ZXYB(JLON,JLEV,YYTXYB%M_DELP)*ZPRESF5D &
         & -PXYB5(JLON,JLEV,YYTXYB%M_DELP)*ZPRESF*ZPRESF5D**2
        ZXYB(JLON,JLEV,YYTXYB%M_RTGR)=YDVAB%VBF(JLEV)*ZPRESFD
        ZXYB(JLON,JLEV,YYTXYB%M_ALPH)=(PRESH5(JLON,JLEV)-ZPRESF5)*ZPRESFD &
         & +(PRESH(JLON,JLEV)-ZPRESF)*ZPRESF5D
      ENDDO
    ENDDO

  ELSE

    IF(YRDYNA%NDLNPR == 0) THEN

      DO JLEV=IFIRST,KFLEV
        !DIR$ IVDEP
        !CDIR NODEP
        DO JLON=KSTART,KPROF
          ZXYB(JLON,JLEV,YYTXYB%M_DELP)=PRESH(JLON,JLEV)-PRESH(JLON,JLEV-1)
          ZXYB(JLON,JLEV,YYTXYB%M_RDELP)=(-1.0_JPRB/PXYB5(JLON,JLEV,YYTXYB%M_DELP)**2)*ZXYB(JLON,JLEV,YYTXYB%M_DELP)
          ZXYB(JLON,JLEV,YYTXYB%M_LNPR)=PRESH(JLON,JLEV)/PRESH5(JLON,JLEV)&
           & -PRESH(JLON,JLEV-1)/PRESH5(JLON,JLEV-1)
          ZXYB(JLON,JLEV,YYTXYB%M_RPRE)=-PRESH(JLON,JLEV)*PXYB5(JLON,JLEV,YYTXYB%M_RPRE)**2
          ZXYB(JLON,JLEV,YYTXYB%M_ALPH)= &
           & -PRESH5(JLON,JLEV-1)*PXYB5(JLON,JLEV,YYTXYB%M_RDELP)*ZXYB(JLON,JLEV,YYTXYB%M_LNPR) &
           & -PRESH5(JLON,JLEV-1)*ZXYB(JLON,JLEV,YYTXYB%M_RDELP)*PXYB5(JLON,JLEV,YYTXYB%M_LNPR) &
           & -PRESH(JLON,JLEV-1)*PXYB5(JLON,JLEV,YYTXYB%M_RDELP)*PXYB5(JLON,JLEV,YYTXYB%M_LNPR)
          ZXYB(JLON,JLEV,YYTXYB%M_RPP)= -PXYB5(JLON,JLEV,YYTXYB%M_RPP)* &
           & ( PRESH(JLON,JLEV)*PXYB5(JLON,JLEV,YYTXYB%M_RPRE) &
           & + PRESH(JLON,JLEV-1)*PXYB5(JLON,JLEV-1,YYTXYB%M_RPRE) )
          ZXYB(JLON,JLEV,YYTXYB%M_RTGR)=ZXYB(JLON,JLEV,YYTXYB%M_RDELP)&
           & *(YDVAB%VDELB(JLEV)+YDVAB%VC(JLEV)*PXYB5(JLON,JLEV,YYTXYB%M_LNPR)*PXYB5(JLON,JLEV,YYTXYB%M_RDELP)) &
           & +PXYB5(JLON,JLEV,YYTXYB%M_RDELP) &
           & *YDVAB%VC(JLEV)*PXYB5(JLON,JLEV,YYTXYB%M_LNPR)*ZXYB(JLON,JLEV,YYTXYB%M_RDELP) &
           & +PXYB5(JLON,JLEV,YYTXYB%M_RDELP) &
           & *YDVAB%VC(JLEV)*ZXYB(JLON,JLEV,YYTXYB%M_LNPR)*PXYB5(JLON,JLEV,YYTXYB%M_RDELP)
        ENDDO
      ENDDO

      DO JLEV=1,IFIRST-1
        !DIR$ IVDEP
        !CDIR NODEP
        DO JLON=KSTART,KPROF
          ZXYB(JLON,JLEV,YYTXYB%M_DELP)=PRESH(JLON,JLEV)-PRESH(JLON,JLEV-1)
          ZXYB(JLON,JLEV,YYTXYB%M_RDELP)=(-1.0_JPRB/PXYB5(JLON,JLEV,YYTXYB%M_DELP)**2)*ZXYB(JLON,JLEV,YYTXYB%M_DELP)
          ZXYB(JLON,JLEV,YYTXYB%M_LNPR)=1.0_JPRB/PRESH5(JLON,1)*PRESH(JLON,1)
          ZXYB(JLON,JLEV,YYTXYB%M_RPRE)=-PRESH(JLON,JLEV)*PXYB5(JLON,JLEV,YYTXYB%M_RPRE)**2
          ZXYB(JLON,JLEV,YYTXYB%M_ALPH)=0.0_JPRB
          ZXYB(JLON,JLEV,YYTXYB%M_RPP)=-PRESH(JLON,JLEV)*PXYB5(JLON,JLEV,YYTXYB%M_RPRE)**2/TOPPRES
          ZXYB(JLON,JLEV,YYTXYB%M_RTGR)=ZXYB(JLON,JLEV,YYTXYB%M_RDELP)*YDVAB%VDELB(JLEV)
        ENDDO
      ENDDO

    ELSEIF(YRDYNA%NDLNPR == 1) THEN

      DO JLEV=IFIRST,KFLEV
        !DIR$ IVDEP
        !CDIR NODEP
        DO JLON=KSTART,KPROF
          ZXYB(JLON,JLEV,YYTXYB%M_DELP)=PRESH(JLON,JLEV)-PRESH(JLON,JLEV-1)
          ZXYB(JLON,JLEV,YYTXYB%M_RDELP)=-1.0_JPRB/(PXYB5(JLON,JLEV,YYTXYB%M_DELP)**2) *ZXYB(JLON,JLEV,YYTXYB%M_DELP)
          ZXYB(JLON,JLEV,YYTXYB%M_RPP)=-1.0_JPRB/(PRESH5(JLON,JLEV)*PRESH5(JLON,JLEV-1)) *&
           & (PRESH(JLON,JLEV)*PRESH5(JLON,JLEV-1)+PRESH5(JLON,JLEV)*PRESH(JLON,JLEV-1))
          ZXYB(JLON,JLEV,YYTXYB%M_LNPR)=(0.5_JPRB*PXYB5(JLON,JLEV,YYTXYB%M_DELP)/PXYB5(JLON,JLEV,YYTXYB%M_RPP) &
           & *ZXYB(JLON,JLEV,YYTXYB%M_RPP)+ZXYB(JLON,JLEV,YYTXYB%M_DELP))*SQRT(PXYB5(JLON,JLEV,YYTXYB%M_RPP))
          ZXYB(JLON,JLEV,YYTXYB%M_ALPH)= &
           & -PRESH5(JLON,JLEV-1)*PXYB5(JLON,JLEV,YYTXYB%M_RDELP)*ZXYB(JLON,JLEV,YYTXYB%M_LNPR) &
           & -PRESH5(JLON,JLEV-1)*ZXYB(JLON,JLEV,YYTXYB%M_RDELP)*PXYB5(JLON,JLEV,YYTXYB%M_LNPR) &
           & -PRESH(JLON,JLEV-1)*PXYB5(JLON,JLEV,YYTXYB%M_RDELP)*PXYB5(JLON,JLEV,YYTXYB%M_LNPR)
          ZXYB(JLON,JLEV,YYTXYB%M_RTGR)=ZXYB(JLON,JLEV,YYTXYB%M_RDELP)&
           & *(YDVAB%VDELB(JLEV)+YDVAB%VC(JLEV)*PXYB5(JLON,JLEV,YYTXYB%M_LNPR)*PXYB5(JLON,JLEV,YYTXYB%M_RDELP)) &
           & + PXYB5(JLON,JLEV,YYTXYB%M_RDELP)*( &
           & YDVAB%VC(JLEV)*ZXYB(JLON,JLEV,YYTXYB%M_LNPR)*PXYB5(JLON,JLEV,YYTXYB%M_RDELP) + &
           & YDVAB%VC(JLEV)*PXYB5(JLON,JLEV,YYTXYB%M_LNPR)*ZXYB(JLON,JLEV,YYTXYB%M_RDELP) )
          ZXYB(JLON,JLEV,YYTXYB%M_RPRE)=-1.0_JPRB/(PRESH5(JLON,JLEV))**2*PRESH(JLON,JLEV)
        ENDDO
      ENDDO

      DO JLEV=1,IFIRST-1
        !DIR$ IVDEP
        !CDIR NODEP
        DO JLON=KSTART,KPROF
          ZXYB(JLON,JLEV,YYTXYB%M_DELP)=PRESH(JLON,JLEV)
          ZXYB(JLON,JLEV,YYTXYB%M_RDELP)=-1.0_JPRB/(PXYB5(JLON,JLEV,YYTXYB%M_DELP))**2 * ZXYB(JLON,JLEV,YYTXYB%M_DELP)
          ZXYB(JLON,JLEV,YYTXYB%M_LNPR)=0.0_JPRB
          ZXYB(JLON,JLEV,YYTXYB%M_ALPH)=0.0_JPRB
          ZXYB(JLON,JLEV,YYTXYB%M_RTGR)=ZXYB(JLON,JLEV,YYTXYB%M_RDELP)*YDVAB%VDELB(JLEV)
          ZXYB(JLON,JLEV,YYTXYB%M_RPRE)=-1.0_JPRB/(PRESH5(JLON,1))**2 * PRESH(JLON,1)
          ZXYB(JLON,JLEV,YYTXYB%M_RPP)=2.0_JPRB*PXYB5(JLON,JLEV,YYTXYB%M_LNPR)*PXYB5(JLON,JLEV,YYTXYB%M_RDELP)* &
           & (PXYB5(JLON,JLEV,YYTXYB%M_LNPR)*ZXYB(JLON,JLEV,YYTXYB%M_RDELP) +&
           & ZXYB(JLON,JLEV,YYTXYB%M_LNPR)*PXYB5(JLON,JLEV,YYTXYB%M_RDELP) )
        ENDDO
      ENDDO

    ENDIF ! NDLNPR

  ENDIF ! LVERTFE

ENDIF

!     ------------------------------------------------------------------

!*       4.    COMPUTES FULL LEVEL PRESSURES
!              -----------------------------

IF (PRESENT(PRESF).AND.PRESENT(PXYB5)) THEN

  IF (LVERTFE) THEN
    DO JLEV=1,KFLEV
      PRESF(KSTART:KPROF,JLEV)=YDVAB%VBF(JLEV)*PRESH(KSTART:KPROF,KFLEV)
    ENDDO
  ELSE
    IF (YRDYNA%NDLNPR == 0) THEN
      IF (YRDYNA%LAPRXPK) THEN
        DO JLEV=1,KFLEV
          DO JLON=KSTART,KPROF
            PRESF(JLON,JLEV)=(PRESH(JLON,JLEV-1)+PRESH(JLON,JLEV))*0.5_JPRB
          ENDDO
        ENDDO
      ELSE
        DO JLEV=1,KFLEV
          DO JLON=KSTART,KPROF
            PRESF(JLON,JLEV)=EXP(-PXYB5(JLON,JLEV,YYTXYB%M_ALPH))*PRESH(JLON,JLEV)&
             & -EXP(-PXYB5(JLON,JLEV,YYTXYB%M_ALPH))*PRESH5(JLON,JLEV)*ZXYB(JLON,JLEV,YYTXYB%M_ALPH)
          ENDDO
        ENDDO
      ENDIF
    ELSEIF (YRDYNA%NDLNPR == 1) THEN
      DO JLEV=IFIRST,KFLEV
        DO JLON=KSTART,KPROF
          PRESF(JLON,JLEV)=(1.0_JPRB-PXYB5(JLON,JLEV,YYTXYB%M_ALPH))*PRESH(JLON,JLEV) &
           & -PRESH5(JLON,JLEV)*ZXYB(JLON,JLEV,YYTXYB%M_ALPH)
        ENDDO
      ENDDO
      ZMUL=1.0_JPRB/(2.0_JPRB+RCVD/RD)
      DO JLEV=1,IFIRST-1
        DO JLON=KSTART,KPROF
          PRESF(JLON,JLEV)=PRESH(JLON,JLEV)*ZMUL
        ENDDO
      ENDDO
    ENDIF ! NDLNPR
  ENDIF ! LVERTFE

ENDIF

IF (.NOT.PRESENT(PXYB)) DEALLOCATE(ZXYB)

IF (LHOOK) CALL DR_HOOK('GPHPRETL',1,ZHOOK_HANDLE)
END SUBROUTINE GPHPRETL
