SUBROUTINE CPPHINP(LDVERTFE,YDGEOMETRY,YDMODEL,KST,KEND,&
 & PGEMU,PGELAM,&
 & PUT0,PVT0,PQT0,PQT0L,PQT0M,PQSLT0L,PQSLT0M,&
 & PRDELP0,PEVEL0,PCVGQSL,&
 & PMU0,PMU0LU,PMU0M,PMU0N,PCVGQ)  

!**** *CPPHINP*  - ComPute PHysical INPut.

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *CPPHINP*

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------

!     Modifications.
!     --------------
!      2002-05, K. YESSAD: consistency with finite element scheme.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      2004-11, Y. Seity: Do not compute Convergence of Humidity for Arome
!      2005-10  Y. Bouteloup : Modification for computation of CVGQ
!        2011-03  Y. Seity: add PMU0N from ARPEGE-climat
!        Y.Bouteloup 21-Jan-2011 : Compute pmu0m as average of pmu0 instead of 
!                                  pmu0(mean time)
!      2011-05-10 E. Bazile: PMU0M=PMU0 for NRADFR=1
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY

!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL           ,INTENT(IN)    :: LDVERTFE
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL)       ,INTENT(IN)    :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSLT0L(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSLT0M(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVEL0(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVGQSL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMU0(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMU0LU(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMU0M(YDGEOMETRY%YRDIM%NPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCVGQ(YDGEOMETRY%YRDIM%NPROMM,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMU0N(YDGEOMETRY%YRDIM%NPROMA)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZED
REAL(KIND=JPRB) :: ZEPS,ZABSTSPHY,ZUSTSPHY
REAL(KIND=JPRB) :: Z1MU0M(YDGEOMETRY%YRDIM%NPROMA),Z2MU0M(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZCOS (YDGEOMETRY%YRDIM%NPROMA), ZSIN (YDGEOMETRY%YRDIM%NPROMA)

LOGICAL :: LLCOMPUTE_CVGQ

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPPHINP',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,&
 &  YDDYN=>YDMODEL%YRML_DYN%YRDYN,YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,YDSIMPHL=>YDMODEL%YRML_PHY_MF%YRSIMPHL, &
 & YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YDARPHY=>YDMODEL%YRML_PHY_MF%YRARPHY, &
  & YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD,YDERIP=>YDMODEL%YRML_PHY_RAD%YRERIP,YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, &
  & YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY)
ASSOCIATE(&
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NCOMP_CVGQ=>YDDYN%NCOMP_CVGQ, &
 & LEPHYS=>YDEPHY%LEPHYS, &
 & NRADFR=>YDERAD%NRADFR, &
 & RCODECM=>YDERIP%RCODECM, RSIDECM=>YDERIP%RSIDECM, RCOVSRM=>YDERIP%RCOVSRM, &
 & RSIVSRM=>YDERIP%RSIVSRM, &
 & RCODEC=>YDRIP%RCODEC, RCODECF=>YDRIP%RCODECF, RCODECLU=>YDRIP%RCODECLU, &
 & RCODECN=>YDRIP%RCODECN, RCOVSR=>YDRIP%RCOVSR, RCOVSRF=>YDRIP%RCOVSRF, &
 & RCOVSRLU=>YDRIP%RCOVSRLU, RCOVSRN=>YDRIP%RCOVSRN, RSIDEC=>YDRIP%RSIDEC, &
 & RSIDECF=>YDRIP%RSIDECF, RSIDECLU=>YDRIP%RSIDECLU, RSIDECN=>YDRIP%RSIDECN, &
 & RSIVSR=>YDRIP%RSIVSR, RSIVSRF=>YDRIP%RSIVSRF, RSIVSRLU=>YDRIP%RSIVSRLU, &
 & RSIVSRN=>YDRIP%RSIVSRN, &
 & LMSE=>YDARPHY%LMSE, LMPA=>YDARPHY%LMPA, &
 & TSPHY=>YDPHY2%TSPHY, &
 & LSIMPH=>YDSIMPHL%LSIMPH, &
 & LRMU0M=>YDPHY%LRMU0M, LRAYLU=>YDPHY%LRAYLU, LRAYFM=>YDPHY%LRAYFM, &
 & LMPHYS=>YDPHY%LMPHYS, LRAYFM15=>YDPHY%LRAYFM15)
!     ------------------------------------------------------------------

!*       1.1   Astronomy.

DO JROF=KST,KEND
  ZCOS (JROF) = COS(PGELAM(JROF))
ENDDO
DO JROF=KST,KEND
  ZSIN (JROF) = SIN(PGELAM(JROF)) 
ENDDO

!DEC$ IVDEP
DO JROF=KST,KEND
  PMU0(JROF)=MAX( RSIDEC*PGEMU(JROF)&
   & -RCODEC*RCOVSR*SQRT(1.0_JPRB-PGEMU(JROF)**2)*ZCOS (JROF)&
   & +RCODEC*RSIVSR*SQRT(1.0_JPRB-PGEMU(JROF)**2)*ZSIN (JROF)&
   & ,0.0_JPRB)  
ENDDO

IF(LMSE)THEN
!DEC$ IVDEP
  DO JROF=KST,KEND
    PMU0N(JROF)=MAX( RSIDECN*PGEMU(JROF)&
     & -RCODECN*RCOVSRN*SQRT(1.0_JPRB-PGEMU(JROF)**2)*ZCOS (JROF)&
     & +RCODECN*RSIVSRN*SQRT(1.0_JPRB-PGEMU(JROF)**2)*ZSIN (JROF)&
     & ,0.0_JPRB)
  ENDDO
ELSE
  DO JROF=KST,KEND
    PMU0N(JROF)=0.0
  ENDDO
ENDIF


!*       Lunar astronomy.

IF(LRAYLU) THEN
!DEC$ IVDEP
  DO JROF=KST,KEND
    PMU0LU(JROF)=MAX( RSIDECLU*PGEMU(JROF)&
     & -RCODECLU*RCOVSRLU*SQRT(1.0_JPRB-PGEMU(JROF)**2)*ZCOS (JROF)&
     & +RCODECLU*RSIVSRLU*SQRT(1.0_JPRB-PGEMU(JROF)**2)*ZSIN (JROF)&
     & ,0.0_JPRB)  
  ENDDO
ENDIF

IF (LEPHYS.OR.(LMPHYS.AND.LRAYFM.AND.(.NOT.LRMU0M))) THEN
!DEC$ IVDEP
  DO JROF=KST,KEND
    PMU0M(JROF)=MAX( RSIDECM*PGEMU(JROF)&
     & -RCODECM*RCOVSRM*SQRT(1.0_JPRB-PGEMU(JROF)**2)*ZCOS (JROF)&
     & +RCODECM*RSIVSRM*SQRT(1.0_JPRB-PGEMU(JROF)**2)*ZSIN (JROF)&
     & ,0.0_JPRB)  
  ENDDO
  IF (NRADFR == 1 ) THEN
     DO JROF=KST,KEND
        PMU0M(JROF)=PMU0(JROF)
     ENDDO
  ENDIF
ELSEIF (LMPHYS.AND.LRAYFM.AND.LRMU0M) THEN
!DEC$ IVDEP
  DO JROF=KST,KEND
    Z1MU0M(JROF) = MAX( RSIDECM*PGEMU(JROF)&
     & -RCODECM*RCOVSRM*SQRT(1.0_JPRB-PGEMU(JROF)**2)*ZCOS (JROF)&
     & +RCODECM*RSIVSRM*SQRT(1.0_JPRB-PGEMU(JROF)**2)*ZSIN (JROF)&
     & ,0.0_JPRB)  
    Z2MU0M(JROF)=0.5_JPRB*(MAX( RSIDECF*PGEMU(JROF)&
     & -RCODECF*RCOVSRF*SQRT(1.0_JPRB-PGEMU(JROF)**2)*ZCOS (JROF)&
     & +RCODECF*RSIVSRF*SQRT(1.0_JPRB-PGEMU(JROF)**2)*ZSIN (JROF)&
     & ,0.0_JPRB)  + PMU0(JROF))
     PMU0M(JROF) = MAX(Z1MU0M(JROF),Z2MU0M(JROF))
  ENDDO

ELSEIF (LMPHYS.AND.LRAYFM15) THEN
  CALL ABOR1('RADIATION SCHEME FROM CYCLE 15 NO LONGER AVAILABLE')
ELSE
  DO JROF=KST,KEND
    PMU0M(JROF)=0.0_JPRB
  ENDDO
ENDIF

!*    1.2   Convergence of humidity and physical wind components.

! ky: variable to be put later in a module, saying if the CVGQ
! calculation is required or not (its definition must be the same
! everywhere in the code).
! For the time being this calculation is needed if MF physics
! (other than AROME physics) is activated.
LLCOMPUTE_CVGQ=(LMPHYS.OR.LSIMPH).AND.(.NOT.LMPA)

IF (LLCOMPUTE_CVGQ) THEN

  IF ((NCOMP_CVGQ == 0) .OR. (NCOMP_CVGQ == 1)) THEN
  
    ! * Eulerian convergence:

    ! a/ Horizontal part of the moisture convergence:
    IF (NCOMP_CVGQ == 0) THEN
      DO JLEV=1,NFLEVG
        DO JROF=KST,KEND
          PCVGQ(JROF,JLEV)=&
           & -PQT0L(JROF,JLEV)*PUT0(JROF,JLEV)&
           & -PQT0M(JROF,JLEV)*PVT0(JROF,JLEV)  
        ENDDO
      ENDDO
    ELSE
      DO JLEV=1,NFLEVG
        DO JROF=KST,KEND
          PCVGQ(JROF,JLEV)=&
           & -PQSLT0L(JROF,JLEV)*PUT0(JROF,JLEV)&
           & -PQSLT0M(JROF,JLEV)*PVT0(JROF,JLEV)  
        ENDDO
      ENDDO
    ENDIF

    ! b/ Vertical part of the moisture convergence:
    IF (LDVERTFE) THEN

      ! * "pevel0" contains full-layer "etadot (d prehyd / d eta)"
      !   The current discretisation of the vertical advection of a variable X
      !   which is proposed is:
      !   [etadot d X/d eta](l) = 0.5 * [1 / Delta eta](l)
      !    * [etadot d prehyd/d eta](l) * (X(l+1)-X(l-1))
      !   Remark K.Y.: this discretisation is not fully consistent with the
      !    discretisation of a vertical advection provided by ECMWF between
      !    CY29 and CY29R2 (first coded in CPDYN, then transferred to CPEULDYN),
      !    and it is desirable (but not urgent) to update this discretisation
      !    in the future to make it consistent with CPEULDYN (in this case
      !    one must use routine VERDER to compute the vertical derivative
      !    [d X/d eta](l)).

      ! * Layer 1:

      DO JROF=KST,KEND
        PCVGQ(JROF,1  )=PCVGQ(JROF,1  )+&
         & (PQT0(JROF,1)-PQT0  (JROF,2))&
         & *PEVEL0(JROF,1)*PRDELP0(JROF,1  )  
      ENDDO

      ! * Layers 2 to L-1:

      DO JLEV=2,NFLEVG-1
        DO JROF=KST,KEND
          ZED=0.5_JPRB*PEVEL0(JROF,JLEV)
 
          PCVGQ(JROF,JLEV  )=PCVGQ(JROF,JLEV  )+&
           & (PQT0(JROF,JLEV-1)-PQT0  (JROF,JLEV+1))&
           & *ZED*PRDELP0(JROF,JLEV  )  
        ENDDO
      ENDDO

      ! * Layer L:

      DO JROF=KST,KEND
        PCVGQ(JROF,NFLEVG  )=PCVGQ(JROF,NFLEVG  )+&
         & (PQT0(JROF,NFLEVG-1)-PQT0  (JROF,NFLEVG))&
         & *PEVEL0(JROF,NFLEVG)*PRDELP0(JROF,NFLEVG  )  
      ENDDO

    ELSE

      ! * "pevel0" contains half-layer "etadot (d prehyd / d eta)":
      !   The discretisation of the vertical advection of a variable X writes:
      !   [etadot d X/d eta](l) =
      !    0.5 * [1 / Delta eta](l) *
      !    { [etadot d prehyd/d eta](lbar) * (X(l+1)-X(l))
      !    + [etadot d prehyd/d eta](lbar-1) * (X(l)-X(l-1)) }

      DO JLEV=1,NFLEVG-1
        DO JROF=KST,KEND
          ZED=0.5_JPRB*PEVEL0(JROF,JLEV)
  
          PCVGQ(JROF,JLEV  )=PCVGQ(JROF,JLEV  )+&
           & (PQT0(JROF,JLEV  )-PQT0  (JROF,JLEV+1))&
           & *ZED*PRDELP0(JROF,JLEV  )  

          PCVGQ(JROF,JLEV+1)=PCVGQ(JROF,JLEV+1)+&
            & (PQT0(JROF,JLEV  )-PQT0  (JROF,JLEV+1))&
            & *ZED*PRDELP0(JROF,JLEV+1)  
        ENDDO
      ENDDO
    ENDIF
  
  ELSE
  
    ! * case NCOMP_CVGQ=2: Semi-Lagrangian convergence:

    ZEPS=1.E-12
    ZABSTSPHY=MAX(ZEPS,ABS(TSPHY))
    IF (TSPHY >=0) THEN
      ZUSTSPHY=1.0_JPRB/ZABSTSPHY
    ELSE
      ZUSTSPHY=-1.0_JPRB/ZABSTSPHY
    ENDIF

    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        PCVGQ(JROF,JLEV)=PCVGQSL(JROF,JLEV)*ZUSTSPHY
      ENDDO
    ENDDO
  ENDIF     

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPPHINP',1,ZHOOK_HANDLE)
END SUBROUTINE CPPHINP
