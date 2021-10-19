SUBROUTINE GPHPRE(KPROMA,KFLEV,K1,K2,VAB,PRESH,PXYB,PRESF,LHSET,LDELP,LALPHA,LRTGR,LRPP)

!**** *GPHPRE* - Computes half and full level pressure
!                Modern version of former GPPRE.
!                Modern version of former GPPREH+GPXYB+GPPREF

!     Purpose.
!     --------
!           Computes pressures at half and full model levels.

!**   Interface.
!     ----------
!        *CALL* *GPHPRE(...)

!        Explicit arguments :
!        --------------------

!          KPROMA    : horizontal dimensioning                                (in)
!          KFLEV     : vertical dimensioning                                  (in)
!          KSTART    : start of work                                          (in)
!          KPROF     : depth of work                                          (in)
!          YDVAB     : contains information about hybrid vertical coordinate  (in)
!          PRESH     : half level pressure                                    (inout)
!          PXYB      : contains pressure depth, "delta", "alpha"              (opt out)
!          PRESF     : full level pressure                                    (opt out)
!          LDELP,LALPHA,... : activation keys for partial computations        (opt in)

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
!      K. YESSAD (Sep 2011) after GPPRE, GPPREH, GPXYB and GPPREF.

!     Modifications.
!     --------------
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (Mar 2017): Introduce NDLNPR=2 for NHQE model.
!   H Petithomme (Dec 2020): add options, use of pointers, group VFE tests
!     ------------------------------------------------------------------

USE PARKIND1,ONLY: JPIM,JPRB
USE YOMHOOK,ONLY: LHOOK,DR_HOOK
USE YOMDYNA,ONLY: LAPRXPK,NDLNPR,RHYDR0
USE YOMCVER,ONLY: LVERTFE
USE YOMCST,ONLY: RD,RCVD
USE YOMVERT,ONLY: TVAB,TOPPRES
USE INTDYN_MOD,ONLY: YYTXYB
USE CODETOOLS,ONLY: SETDEFAULTL

IMPLICIT NONE

TYPE(TVAB),INTENT(IN) :: VAB
INTEGER(KIND=JPIM),INTENT(IN) :: KPROMA,KFLEV,K1,K2
LOGICAL,OPTIONAL,INTENT(IN) :: LHSET,LDELP,LALPHA,LRTGR,LRPP
REAL(KIND=JPRB),INTENT(INOUT) :: PRESH(KPROMA,0:KFLEV)
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(OUT) :: PXYB(KPROMA,KFLEV,YYTXYB%NDIM)
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(OUT) :: PRESF(KPROMA,KFLEV)

INTEGER(KIND=JPIM) :: IFIRST,JL,JP
LOGICAL :: LLHSET,LLDELP,LLALPHA,LLRTGR,LLRPP,LTEST
REAL(KIND=JPRB) :: ZPRE(KPROMA),ZHOOK_HANDLE
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZXYB(:,:,:),ZPRESF(:,:)
REAL(KIND=JPRB),ALLOCATABLE,TARGET :: ZZXYB(:,:,:),ZZPRESF(:,:)
REAL(KIND=JPRB), POINTER :: ZDELP (:,:), ZLNPR (:,:), ZRDELP (:,:), ZALPH (:,:), ZRTGR (:,:), ZRPRE (:,:), ZRPP (:,:)  

IF (LHOOK) CALL DR_HOOK('GPHPRE',0,ZHOOK_HANDLE)

LLHSET = SETDEFAULTL(.FALSE.,LHSET)

IF (.NOT.LLHSET) THEN
  DO JL=0,KFLEV-1
    PRESH(K1:K2,JL) = VAB%VAH(JL)+VAB%VBH(JL)*PRESH(K1:K2,KFLEV)
  ENDDO
ENDIF

ASSOCIATE(X=>YYTXYB)

IF (LVERTFE) THEN
  IF (PRESENT(PXYB)) THEN
    LLDELP = SETDEFAULTL(.TRUE.,LDELP)
    LLALPHA = SETDEFAULTL(.TRUE.,LALPHA)
    LLRTGR = SETDEFAULTL(.TRUE.,LRTGR)

    ZDELP  => PXYB (:,:,X%M_DELP )
    ZLNPR  => PXYB (:,:,X%M_LNPR )
    ZRDELP => PXYB (:,:,X%M_RDELP)
    ZALPH  => PXYB (:,:,X%M_ALPH )
    ZRTGR  => PXYB (:,:,X%M_RTGR )
    ZRPRE  => PXYB (:,:,X%M_RPRE )
    ZRPP   => PXYB (:,:,X%M_RPP  )

    IF (PRESENT(PRESF)) THEN
      ZPRESF => PRESF(:,:)
    ELSE
      ALLOCATE(ZZPRESF(KPROMA,KFLEV))
      ZPRESF => ZZPRESF
    ENDIF

    ! useless:
    !PXYB(:,:,:) = 0

    DO JL=1,KFLEV
      !DIR$ IVDEP
      !CDIR NODEP
      DO JP=K1,K2
        PXYB(JP,JL,X%M_DELP) = VAB%VDELA(JL)+VAB%VDELB(JL)*PRESH(JP,KFLEV)
        ZPRESF(JP,JL) = VAB%VAF(JL)+VAB%VBF(JL)*PRESH(JP,KFLEV)
        ZPRE(JP) = 1._JPRB/ZPRESF(JP,JL)
        PXYB(JP,JL,X%M_LNPR) = PXYB(JP,JL,X%M_DELP)*ZPRE(JP)
      ENDDO

      IF (LLDELP) THEN
        !DIR$ IVDEP
        !CDIR NODEP
        DO JP=K1,K2
          PXYB(JP,JL,X%M_RDELP) = 1._JPRB/PXYB(JP,JL,X%M_DELP)
        ENDDO
      ENDIF

      IF (LLALPHA) THEN
        !DIR$ IVDEP
        !CDIR NODEP
        DO JP=K1,K2
          PXYB(JP,JL,X%M_ALPH) = (PRESH(JP,JL)-ZPRESF(JP,JL))*ZPRE(JP)
        ENDDO
      ENDIF

      IF (LLRTGR) THEN
        !DIR$ IVDEP
        !CDIR NODEP
        DO JP=K1,K2
          PXYB(JP,JL,X%M_RTGR) = VAB%VBF(JL)*ZPRE(JP)
        ENDDO
      ENDIF
    ENDDO
  ELSE IF (PRESENT(PRESF)) THEN
    DO JL=1,KFLEV
      !DIR$ IVDEP
      !CDIR NODEP
      DO JP=K1,K2
        PRESF(JP,JL) = VAB%VAF(JL)+VAB%VBF(JL)*PRESH(JP,KFLEV)
      ENDDO
    ENDDO
  ENDIF
ELSE
  IF (PRESENT(PXYB)) THEN
    LLDELP = SETDEFAULTL(.TRUE.,LDELP)
    LLALPHA = SETDEFAULTL(.TRUE.,LALPHA)
    LLRTGR = SETDEFAULTL(.TRUE.,LRTGR)
    LLRPP = SETDEFAULTL(.TRUE.,LRPP)

    ! broader condition for computing alpha: 
    IF (PRESENT(PRESF)) LLALPHA = LLALPHA.OR.NDLNPR == 1.OR.NDLNPR == 2.OR..NOT.LAPRXPK

    ZXYB => PXYB(:,:,:)
  ELSE
    LLDELP = .FALSE.
    LLRTGR = .FALSE.
    LLRPP = .FALSE.

    ! reduced condition for computing alpha: 
    LLALPHA = PRESENT(PRESF).AND.(NDLNPR == 1.OR.NDLNPR == 2.OR..NOT.LAPRXPK)

    IF (LLALPHA) THEN
      ALLOCATE(ZZXYB(KPROMA,KFLEV,X%NDIM))
      ZXYB => ZZXYB
    ENDIF
  ENDIF

  ZDELP  => ZXYB (:,:,X%M_DELP )
  ZLNPR  => ZXYB (:,:,X%M_LNPR )
  ZRDELP => ZXYB (:,:,X%M_RDELP)
  ZALPH  => ZXYB (:,:,X%M_ALPH )
  ZRTGR  => ZXYB (:,:,X%M_RTGR )
  ZRPRE  => ZXYB (:,:,X%M_RPRE )
  ZRPP   => ZXYB (:,:,X%M_RPP  )

  IF (PRESENT(PXYB).OR.LLALPHA) THEN
    ! pressure at top
    ZPRE(K1:K2) = VAB%VAH(0)+VAB%VBH(0)*PRESH(K1:K2,KFLEV)

    DO JP=K1,K2
      IF (ZPRE(JP) <= TOPPRES) EXIT
    END DO

    IF (JP > K2) THEN
      IFIRST = 1
    ELSE
      IFIRST = 2
    ENDIF

    ! useless:
    !ZXYB(:,:,:) = 0

    IF (NDLNPR == 0) THEN
      IF (IFIRST == 2) THEN
        ZXYB(K1:K2,1,X%M_LNPR) = LOG(PRESH(K1:K2,1)/TOPPRES)

        IF (LLDELP.OR.LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZXYB(JP,1,X%M_DELP) = PRESH(JP,1)-PRESH(JP,0)
            ZXYB(JP,1,X%M_RDELP) = 1._JPRB/ZXYB(JP,1,X%M_DELP)
          ENDDO
        ENDIF

        IF (LLALPHA) ZXYB(K1:K2,1,X%M_ALPH) = RHYDR0

        IF (LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZXYB(JP,1,X%M_RTGR) = ZXYB(JP,1,X%M_RDELP)*VAB%VDELB(1)
          ENDDO
        ENDIF

        IF (LLRPP) THEN
          DO JP=K1,K2
            ZXYB(JP,1,X%M_RPRE) = 1._JPRB/PRESH(JP,1)
            ZXYB(JP,1,X%M_RPP) = 1._JPRB/(PRESH(JP,1)*TOPPRES)
          ENDDO
        ENDIF
      ENDIF

      ZPRE(K1:K2) = 1._JPRB/PRESH(K1:K2,IFIRST-1)

      DO JL=IFIRST,KFLEV
        ZXYB(K1:K2,JL,X%M_LNPR) = LOG(PRESH(K1:K2,JL)*ZPRE(K1:K2))

        IF (LLDELP.OR.LLALPHA.OR.LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZXYB(JP,JL,X%M_DELP) = PRESH(JP,JL)-PRESH(JP,JL-1)
            ZXYB(JP,JL,X%M_RDELP) = 1._JPRB/ZXYB(JP,JL,X%M_DELP)
          ENDDO
        ENDIF

        IF (LLALPHA) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZXYB(JP,JL,X%M_ALPH) = 1._JPRB-PRESH(JP,JL-1)*ZXYB(JP,JL,X%M_RDELP)*&
              ZXYB(JP,JL,X%M_LNPR)
          ENDDO
        ENDIF

        IF (LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZXYB(JP,JL,X%M_RTGR) = ZXYB(JP,JL,X%M_RDELP)*(VAB%VDELB(JL)+&
              VAB%VC(JL)*ZXYB(JP,JL,X%M_LNPR)*ZXYB(JP,JL,X%M_RDELP))
          ENDDO
        ENDIF

        IF (LLRPP) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZXYB(JP,JL,X%M_RPRE) = 1._JPRB/PRESH(JP,JL)
            ZXYB(JP,JL,X%M_RPP) = ZXYB(JP,JL,X%M_RPRE)*ZPRE(JP)
            ZPRE(JP) = ZXYB(JP,JL,X%M_RPRE)
          ENDDO
        ELSE
          ZPRE(K1:K2) = 1._JPRB/PRESH(K1:K2,JL)
        ENDIF
      ENDDO
    ELSE IF (NDLNPR == 1.OR.NDLNPR == 2) THEN
      LTEST = LLDELP.OR.LLALPHA.OR.LLRTGR

      DO JL=IFIRST,KFLEV
        ! optim: keep lnpr separate but consistent
        IF (LTEST) THEN
          DO JP=K1,K2
            ZXYB(JP,JL,X%M_DELP) = PRESH(JP,JL)-PRESH(JP,JL-1)
            ZXYB(JP,JL,X%M_RDELP) = 1._JPRB/ZXYB(JP,JL,X%M_DELP)
            ZXYB(JP,JL,X%M_RPP) = 1._JPRB/(PRESH(JP,JL)*PRESH(JP,JL-1))
            ZXYB(JP,JL,X%M_LNPR) = ZXYB(JP,JL,X%M_DELP)*SQRT(ZXYB(JP,JL,X%M_RPP))
          ENDDO
        ELSE
          DO JP=K1,K2
            ZXYB(JP,JL,X%M_LNPR) = (PRESH(JP,JL)-PRESH(JP,JL-1))/&
              SQRT(PRESH(JP,JL)*PRESH(JP,JL-1))
          ENDDO
        ENDIF

        IF (LLALPHA) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZXYB(JP,JL,X%M_ALPH) = 1._JPRB-PRESH(JP,JL-1)*ZXYB(JP,JL,X%M_RDELP)*&
              ZXYB(JP,JL,X%M_LNPR)
          ENDDO
        ENDIF

        IF (LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZXYB(JP,JL,X%M_RTGR) = ZXYB(JP,JL,X%M_RDELP)*(VAB%VDELB(JL)+&
              VAB%VC(JL)*ZXYB(JP,JL,X%M_LNPR)*ZXYB(JP,JL,X%M_RDELP))
          ENDDO
        ENDIF

        IF (LLRPP) ZXYB(K1:K2,JL,X%M_RPRE) = 1._JPRB/PRESH(K1:K2,JL)
      ENDDO

      IF (IFIRST == 2) THEN
        !DIR$ IVDEP
        !CDIR NODEP
        DO JP=K1,K2
          ZXYB(JP,1,X%M_DELP) = PRESH(JP,1)
          ZXYB(JP,1,X%M_RDELP) = 1._JPRB/ZXYB(JP,1,X%M_DELP)
        ENDDO

        IF (NDLNPR == 1) THEN
          ZXYB(K1:K2,1,X%M_LNPR) = 2._JPRB+RCVD/RD
        ELSE
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZXYB(JP,1,X%M_LNPR) = 1._JPRB+ZXYB(JP,2,X%M_LNPR)*&
              (ZXYB(JP,1,X%M_DELP)/ZXYB(JP,2,X%M_DELP))*SQRT(PRESH(JP,2)/PRESH(JP,1))
          ENDDO
        ENDIF

        IF (LLALPHA) ZXYB(K1:K2,1,X%M_ALPH) = 1._JPRB

        IF (LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZXYB(JP,1,X%M_RTGR) = ZXYB(JP,1,X%M_RDELP)*VAB%VDELB(1)
          ENDDO
        ENDIF

        IF (LLRPP) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZXYB(JP,1,X%M_RPRE) = 1._JPRB/PRESH(JP,1)
            ZXYB(JP,1,X%M_RPP) = (ZXYB(JP,1,X%M_LNPR)*ZXYB(JP,1,X%M_RDELP))**2
          ENDDO
        ENDIF
      ENDIF
    ENDIF
  ENDIF

  IF (PRESENT(PRESF)) THEN
    IF (NDLNPR == 0) THEN
      IF (LAPRXPK) THEN
        DO JL=1,KFLEV
          PRESF(K1:K2,JL) = (PRESH(K1:K2,JL-1)+PRESH(K1:K2,JL))*0.5_JPRB
        ENDDO
      ELSE
        DO JL=1,KFLEV
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            PRESF(JP,JL) = EXP(-ZXYB(JP,JL,X%M_ALPH))*PRESH(JP,JL)
          ENDDO
        ENDDO
      ENDIF
    ELSE IF (NDLNPR == 1.OR.NDLNPR == 2) THEN
      DO JL=IFIRST,KFLEV
        !DIR$ IVDEP
        !CDIR NODEP
        DO JP=K1,K2
          PRESF(JP,JL) = (1._JPRB-ZXYB(JP,JL,X%M_ALPH))*PRESH(JP,JL)
        ENDDO
      ENDDO

      IF (IFIRST == 2) THEN
        !DIR$ IVDEP
        !CDIR NODEP
        DO JP=K1,K2
          PRESF(JP,1) = PRESH(JP,1)/ZXYB(JP,1,X%M_LNPR)
        ENDDO
      ENDIF
    ENDIF
  ENDIF

ENDIF
END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('GPHPRE',1,ZHOOK_HANDLE)
END SUBROUTINE

