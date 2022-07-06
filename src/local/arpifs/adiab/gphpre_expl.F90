SUBROUTINE GPHPRE_EXPL(LDAPRXPK, LDVERTFE, KDLNPR, RHYDR0, TOPPRES, YDCST,KPROMA,KFLEV,K1,K2,VAB,PRESH,PRESF,LHSET,LDELP,LALPHA,LRTGR,LRPP,&
                     & PDELP, PLNPR, PRDELP, PALPH, PRTGR, PRPRE, PRPP)

!**** *GPHPRE_EXPL* - Computes half and full level pressure
!                Modern version of former GPPRE.
!                Modern version of former GPPREH+GPXYB+GPPREF

!     Purpose.
!     --------
!           Computes pressures at half and full model levels.

!**   Interface.
!     ----------
!        *CALL* *GPHPRE_EXPL(...)

!        Explicit arguments :
!        --------------------

!          KPROMA    : horizontal dimensioning                                (in)
!          KFLEV     : vertical dimensioning                                  (in)
!          KSTART    : start of work                                          (in)
!          KPROF     : depth of work                                          (in)
!          YDVAB     : contains information about hybrid vertical coordinate  (in)
!          PRESH     : half level pressure                                    (inout)
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

USE PARKIND1,ONLY: JPIM, JPRB
USE YOMHOOK, ONLY: LHOOK, DR_HOOK

USE YOMCST,    ONLY: TCST
USE YOMVERT,   ONLY: TVAB

USE CODETOOLS, ONLY: SETDEFAULTL





IMPLICIT NONE

LOGICAL, INTENT (IN) :: LDAPRXPK
LOGICAL, INTENT (IN) :: LDVERTFE
INTEGER (KIND=JPIM), INTENT (IN) :: KDLNPR
REAL (KIND=JPRB), INTENT (IN) :: RHYDR0
REAL (KIND=JPRB), INTENT (IN) :: TOPPRES
TYPE(TCST),                     INTENT(IN)    :: YDCST
TYPE(TVAB),                     INTENT(IN)    :: VAB
INTEGER(KIND=JPIM),             INTENT(IN)    :: KPROMA,KFLEV,K1,K2
LOGICAL,OPTIONAL,               INTENT(IN)    :: LHSET,LDELP,LALPHA,LRTGR,LRPP
REAL(KIND=JPRB),                INTENT(INOUT) :: PRESH(KPROMA,0:KFLEV)
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(OUT)   :: PRESF(KPROMA,KFLEV)
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(OUT)   :: PDELP(KPROMA,KFLEV), PLNPR(KPROMA,KFLEV), PRDELP(KPROMA,KFLEV), PALPH(KPROMA,KFLEV)
REAL(KIND=JPRB),TARGET,OPTIONAL,INTENT(OUT)   :: PRTGR(KPROMA,KFLEV), PRPRE(KPROMA,KFLEV), PRPP(KPROMA,KFLEV)

INTEGER(KIND=JPIM) :: IFIRST,JL,JP
LOGICAL :: LLHSET,LLDELP,LLALPHA,LLRTGR,LLRPP,LTEST,LLXYB
REAL(KIND=JPRB) :: ZPRE(KPROMA)
REAL(KIND=JPRB),CONTIGUOUS,POINTER :: ZPRESF(:,:)
REAL(KIND=JPRB),TARGET :: ZZPRESF(KPROMA,KFLEV)
REAL(KIND=JPRB),TARGET :: ZZDELP(KPROMA,KFLEV), ZZLNPR(KPROMA,KFLEV), ZZRDELP(KPROMA,KFLEV), ZZALPH(KPROMA,KFLEV)
REAL(KIND=JPRB),TARGET :: ZZRTGR(KPROMA,KFLEV), ZZRPRE(KPROMA,KFLEV), ZZRPP(KPROMA,KFLEV)
REAL(KIND=JPRB), POINTER :: ZDELP (:,:), ZLNPR (:,:), ZRDELP (:,:), ZALPH (:,:), ZRTGR (:,:), ZRPRE (:,:), ZRPP (:,:)  

REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('GPHPRE_EXPL',0,ZHOOK_HANDLE)

LLHSET = SETDEFAULTL(.FALSE.,LHSET)
LLXYB  = &
& (PRESENT (PDELP) .AND.  PRESENT (PLNPR) .AND.  PRESENT (PRDELP) .AND.  PRESENT (PALPH) .AND. &
&  PRESENT (PRTGR) .AND.  PRESENT (PRPRE) .AND.  PRESENT (PRPP))

ZPRESF => NULL ()

IF (.NOT.LLHSET) THEN
  DO JL=0,KFLEV-1
    PRESH(K1:K2,JL) = VAB%VAH(JL)+VAB%VBH(JL)*PRESH(K1:K2,KFLEV)
  ENDDO
ENDIF

IF (LDVERTFE) THEN
  IF (LLXYB) THEN
    LLDELP = SETDEFAULTL(.TRUE.,LDELP)
    LLALPHA = SETDEFAULTL(.TRUE.,LALPHA)
    LLRTGR = SETDEFAULTL(.TRUE.,LRTGR)

    ZDELP  => PDELP 
    ZLNPR  => PLNPR 
    ZRDELP => PRDELP
    ZALPH  => PALPH 
    ZRTGR  => PRTGR 
    ZRPRE  => PRPRE 
    ZRPP   => PRPP  

    IF (PRESENT(PRESF)) THEN
      ZPRESF => PRESF(:,:)
    ELSE
      ZPRESF => ZZPRESF
    ENDIF

    DO JL=1,KFLEV
      !DIR$ IVDEP
      !CDIR NODEP
      DO JP=K1,K2
        ZDELP(JP,JL) = VAB%VDELA(JL)+VAB%VDELB(JL)*PRESH(JP,KFLEV)
        ZPRESF(JP,JL) = VAB%VAF(JL)+VAB%VBF(JL)*PRESH(JP,KFLEV)
        ZPRE(JP) = 1._JPRB/ZPRESF(JP,JL)
        ZLNPR(JP,JL) = ZDELP(JP,JL)*ZPRE(JP)
      ENDDO

      IF (LLDELP) THEN
        !DIR$ IVDEP
        !CDIR NODEP
        DO JP=K1,K2
          ZRDELP(JP,JL) = 1._JPRB/ZDELP(JP,JL)
        ENDDO
      ENDIF

      IF (LLALPHA) THEN
        !DIR$ IVDEP
        !CDIR NODEP
        DO JP=K1,K2
          ZALPH(JP,JL) = (PRESH(JP,JL)-ZPRESF(JP,JL))*ZPRE(JP)
        ENDDO
      ENDIF

      IF (LLRTGR) THEN
        !DIR$ IVDEP
        !CDIR NODEP
        DO JP=K1,K2
          ZRTGR(JP,JL) = VAB%VBF(JL)*ZPRE(JP)
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
  IF (LLXYB) THEN
    LLDELP = SETDEFAULTL(.TRUE.,LDELP)
    LLALPHA = SETDEFAULTL(.TRUE.,LALPHA)
    LLRTGR = SETDEFAULTL(.TRUE.,LRTGR)
    LLRPP = SETDEFAULTL(.TRUE.,LRPP)

    ! broader condition for computing alpha: 
    IF (PRESENT(PRESF)) LLALPHA = LLALPHA.OR.KDLNPR == 1.OR.KDLNPR == 2.OR..NOT.LDAPRXPK

    ZDELP  => PDELP  
    ZLNPR  => PLNPR  
    ZRDELP => PRDELP 
    ZALPH  => PALPH  
    ZRTGR  => PRTGR  
    ZRPRE  => PRPRE  
    ZRPP   => PRPP   

  ELSE
    LLDELP = .FALSE.
    LLRTGR = .FALSE.
    LLRPP = .FALSE.

    ! reduced condition for computing alpha: 
    LLALPHA = PRESENT(PRESF).AND.(KDLNPR == 1.OR.KDLNPR == 2.OR..NOT.LDAPRXPK)

    IF (LLALPHA) THEN
      ZDELP  => ZZDELP  
      ZLNPR  => ZZLNPR  
      ZRDELP => ZZRDELP 
      ZALPH  => ZZALPH  
      ZRTGR  => ZZRTGR  
      ZRPRE  => ZZRPRE  
      ZRPP   => ZZRPP   
    ENDIF
  ENDIF

  IF (LLXYB.OR.LLALPHA) THEN
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

    IF (KDLNPR == 0) THEN
      IF (IFIRST == 2) THEN
        ZLNPR(K1:K2,1) = LOG(PRESH(K1:K2,1)/TOPPRES)

        IF (LLDELP.OR.LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZDELP(JP,1) = PRESH(JP,1)-PRESH(JP,0)
            ZRDELP(JP,1) = 1._JPRB/ZDELP(JP,1)
          ENDDO
        ENDIF

        IF (LLALPHA) ZALPH(K1:K2,1) = RHYDR0

        IF (LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZRTGR(JP,1) = ZRDELP(JP,1)*VAB%VDELB(1)
          ENDDO
        ENDIF

        IF (LLRPP) THEN
          DO JP=K1,K2
            ZRPRE(JP,1) = 1._JPRB/PRESH(JP,1)
            ZRPP(JP,1) = 1._JPRB/(PRESH(JP,1)*TOPPRES)
          ENDDO
        ENDIF
      ENDIF

      ZPRE(K1:K2) = 1._JPRB/PRESH(K1:K2,IFIRST-1)

      DO JL=IFIRST,KFLEV
        ZLNPR(K1:K2,JL) = LOG(PRESH(K1:K2,JL)*ZPRE(K1:K2))

        IF (LLDELP.OR.LLALPHA.OR.LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZDELP(JP,JL) = PRESH(JP,JL)-PRESH(JP,JL-1)
            ZRDELP(JP,JL) = 1._JPRB/ZDELP(JP,JL)
          ENDDO
        ENDIF

        IF (LLALPHA) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZALPH(JP,JL) = 1._JPRB-PRESH(JP,JL-1)*ZRDELP(JP,JL)*&
              ZLNPR(JP,JL)
          ENDDO
        ENDIF

        IF (LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZRTGR(JP,JL) = ZRDELP(JP,JL)*(VAB%VDELB(JL)+&
              VAB%VC(JL)*ZLNPR(JP,JL)*ZRDELP(JP,JL))
          ENDDO
        ENDIF

        IF (LLRPP) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZRPRE(JP,JL) = 1._JPRB/PRESH(JP,JL)
            ZRPP(JP,JL) = ZRPRE(JP,JL)*ZPRE(JP)
            ZPRE(JP) = ZRPRE(JP,JL)
          ENDDO
        ELSE
          ZPRE(K1:K2) = 1._JPRB/PRESH(K1:K2,JL)
        ENDIF
      ENDDO
    ELSE IF (KDLNPR == 1.OR.KDLNPR == 2) THEN
      LTEST = LLDELP.OR.LLALPHA.OR.LLRTGR

      DO JL=IFIRST,KFLEV
        ! optim: keep lnpr separate but consistent
        IF (LTEST) THEN
          DO JP=K1,K2
            ZDELP(JP,JL) = PRESH(JP,JL)-PRESH(JP,JL-1)
            ZRDELP(JP,JL) = 1._JPRB/ZDELP(JP,JL)
            ZRPP(JP,JL) = 1._JPRB/(PRESH(JP,JL)*PRESH(JP,JL-1))
            ZLNPR(JP,JL) = ZDELP(JP,JL)*SQRT(ZRPP(JP,JL))
          ENDDO
        ELSE
          DO JP=K1,K2
            ZLNPR(JP,JL) = (PRESH(JP,JL)-PRESH(JP,JL-1))/&
              SQRT(PRESH(JP,JL)*PRESH(JP,JL-1))
          ENDDO
        ENDIF

        IF (LLALPHA) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZALPH(JP,JL) = 1._JPRB-PRESH(JP,JL-1)*ZRDELP(JP,JL)*&
              ZLNPR(JP,JL)
          ENDDO
        ENDIF

        IF (LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZRTGR(JP,JL) = ZRDELP(JP,JL)*(VAB%VDELB(JL)+&
              VAB%VC(JL)*ZLNPR(JP,JL)*ZRDELP(JP,JL))
          ENDDO
        ENDIF

        IF (LLRPP) ZRPRE(K1:K2,JL) = 1._JPRB/PRESH(K1:K2,JL)
      ENDDO

      IF (IFIRST == 2) THEN
        !DIR$ IVDEP
        !CDIR NODEP
        DO JP=K1,K2
          ZDELP(JP,1) = PRESH(JP,1)
          ZRDELP(JP,1) = 1._JPRB/ZDELP(JP,1)
        ENDDO

        IF (KDLNPR == 1) THEN
          ZLNPR(K1:K2,1) = 2._JPRB+YDCST%RCVD/YDCST%RD
        ELSE
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZLNPR(JP,1) = 1._JPRB+ZLNPR(JP,2)*&
              (ZDELP(JP,1)/ZDELP(JP,2))*SQRT(PRESH(JP,2)/PRESH(JP,1))
          ENDDO
        ENDIF

        IF (LLALPHA) ZALPH(K1:K2,1) = 1._JPRB

        IF (LLRTGR) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZRTGR(JP,1) = ZRDELP(JP,1)*VAB%VDELB(1)
          ENDDO
        ENDIF

        IF (LLRPP) THEN
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            ZRPRE(JP,1) = 1._JPRB/PRESH(JP,1)
            ZRPP(JP,1) = (ZLNPR(JP,1)*ZRDELP(JP,1))**2
          ENDDO
        ENDIF
      ENDIF
    ENDIF
  ENDIF

  IF (PRESENT(PRESF)) THEN
    IF (KDLNPR == 0) THEN
      IF (LDAPRXPK) THEN
        DO JL=1,KFLEV
          PRESF(K1:K2,JL) = (PRESH(K1:K2,JL-1)+PRESH(K1:K2,JL))*0.5_JPRB
        ENDDO
      ELSE
        DO JL=1,KFLEV
          !DIR$ IVDEP
          !CDIR NODEP
          DO JP=K1,K2
            PRESF(JP,JL) = EXP(-ZALPH(JP,JL))*PRESH(JP,JL)
          ENDDO
        ENDDO
      ENDIF
    ELSE IF (KDLNPR == 1.OR.KDLNPR == 2) THEN
      DO JL=IFIRST,KFLEV
        !DIR$ IVDEP
        !CDIR NODEP
        DO JP=K1,K2
          PRESF(JP,JL) = (1._JPRB-ZALPH(JP,JL))*PRESH(JP,JL)
        ENDDO
      ENDDO

      IF (IFIRST == 2) THEN
        !DIR$ IVDEP
        !CDIR NODEP
        DO JP=K1,K2
          PRESF(JP,1) = PRESH(JP,1)/ZLNPR(JP,1)
        ENDDO
      ENDIF
    ENDIF
  ENDIF

ENDIF

IF (LHOOK) CALL DR_HOOK('GPHPRE_EXPL',1,ZHOOK_HANDLE)

END SUBROUTINE

