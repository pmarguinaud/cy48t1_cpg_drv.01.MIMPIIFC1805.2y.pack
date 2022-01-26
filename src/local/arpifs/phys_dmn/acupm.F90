!OPTIONS XOPT(NOEVAL)
SUBROUTINE ACUPM( YDML_PHY_MF,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
! ----------------------------------------------------------------------
! - INPUT  2D .
                &PCP, PLHS, PLHV,&
                &PEVEL0, PIPOI, PPOID, PUDAL, PUDOM,&
                &PFPEVPSL,PFPEVPSN,PFPEVPSG,PFPFPSL,PFPFPSN,PFPFPSG,&
                &PFPLSL,PFPLSN,PFPLSG,&
                &PZFCQL,PZFCQI,PFCCQL,PFCCQN,&
! - OUTPUT 2D
                &POME,&
                &PFHP, &
! - INPUT/OUTPUT  2D .
                &PZQV,PZQL,PZQI,PZQR,PZQS,PZQG,PZT,&
                &PFCSQL,PFCSQN,&
                &PFCQVNG,PFCQRNG,PFCQSNG,PFCQGNG)
! **** *ACUPM * - UPDATE AFTER MICROPHYSICS - 3MT

! **   Interface.
!      ----------
! *CALL* *ACUPM*
! ----------------------------------------------------------------------
! - INPUT ARGUMENTS
! -------------------

! - DIMENSIONING PARAMETERS FOR PHYSICS

! KIDIA      : INDICE DE DEPART DES BOUCLES VECTORISEES SUR L'HORIZONT..
! KFDIA      : INDICE DE FIN DES BOUCLES VECTORISEES SUR L'HORIZONTALE.
! KLON       : DIMENSION HORIZONTALE DES TABLEAUX.
! KTDIA      : INDICE DE DEPART DES BOUCLES VERTICALES (NIVEL DU SOMMET
!            : DU NUAGE).
! KLEV       : DIMENSION VERTICALE DES TABLEAUX "FULL LEVEL".


! - PHYSICAL VARIABLES

! - 2D (1:KLEV) .
! PCP
! PLHS
! PLHV
! PEVEL0    : ETA_dot*d(Pi)/d(ETA)at full levels
! PIPOI     : g*dt/delta p
! PPOID     : delta p/g*dt
! PUDAL     : updraught mesh fraction
! PUDOM     : relative updraught vertical velocity

! - 2D FLUXES (0:KLEV) .
! PFPEVPSL  : resolved evaporation of rain.
! PFPEVPSN  : resolved evaporation of snow.
! PFPFPSL   : resolved autoconversion to rain.
! PFPFPSN   : resolved autoconversion to snow.
! PFPLSL    : flux of rain.
! PFPLSN    : flux of snow.
! PZFCQL    : total liquid condensation flux.
! PZFCQI    : total solid condensation flux.
! PFCCQL    : convective liquid condensation flux.
! PFCCQN    : convective solid condensation flux.

! ----------------------------------------------------------------------
! -   OUTPUT ARGUMENTS 
!     ----------------
! - 1D
! - 2D
! POME   : Vertical velocity *dt in updraught environment (Pa)
! PFHP   : Latent heat flux associated to precipitation
!           (limited "downdraught closure" edition)
! ----------------------------------------------------------------------
! -   INPUT/OUTPUT ARGUMENTS 
!     ----------------------
! - 2D (1:KLEV)
! PZQI :  cloud ice of the cascade.
! PZQL :  cloud water of the cascade.
! PZQV :  water vapour of the cascade.
! PZQR :  rain of the cascade.
! PZQS :  snow of the cascade.
! PZT  :  temperature of the cascade.

! - 2D FLUXES (0:KLEV)
! PFCSQL  : resolved liquid condensation flux.
! PFCSQN  : resolved solid condensation flux.
! PFCQVNG : water vapour negative values
! PFCQRNG : rain negative values
! PFCQSNG : snow negative values
! PFCQGNG : graupel negative values.

! ----------------------------------------------------------------------
!     Author.
!     -------
!     01-2007, L. Gerard 

!     Modifications.
!     --------------
!     04-2007, L. Gerard 
!     A.Bogatchev : 26-Oct-2007 removing of severe breaks of coding rools
!     K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!     R.Brozkova  : 29-Sep-2009 cascade update and cleaning LUDEN
! ----------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE PARKIND1 , ONLY: JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(IN):: YDML_PHY_MF
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KTDIA
REAL(KIND=JPRB) ,INTENT(IN) :: PCP(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PLHS(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PLHV(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PEVEL0(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PIPOI(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PPOID(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PUDAL(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PUDOM(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PFPEVPSL(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PFPEVPSN(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PFPEVPSG(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PFPFPSL(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PFPFPSN(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PFPFPSG(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PFPLSL(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PFPLSN(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PFPLSG(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PZFCQL(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PZFCQI(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PFCCQL(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PFCCQN(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: POME(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PFHP(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PZQI(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PZQL(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PZQV(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PZQR(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PZQS(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PZQG(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PZT(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PFCSQL(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PFCSQN(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PFCQVNG(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PFCQRNG(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PFCQSNG(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PFCQGNG(KLON,0:KLEV)

! ----------------------------------------------------------------------

INTEGER (KIND=JPIM) :: JLEV, JLON
REAL(KIND=JPRB)  :: ZDDEVC, ZDEVPL, ZDEVPN, ZUDAL, ZQX1, ZDQR, ZDQS, ZDQG, ZDQC,&
                   &ZDFCQL, ZDFCQI, ZQV0, ZDQV
REAL(KIND=JPRB)  :: ZFCQRNG(KLON,0:KLEV),ZFCQSNG(KLON,0:KLEV),ZFCQGNG(KLON,0:KLEV),&
                   &ZFCQVNG(KLON,0:KLEV),ZFCQIDM(KLON,0:KLEV),&
                   &ZFCQLDM(KLON,0:KLEV)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACUPM',0,ZHOOK_HANDLE)
ASSOCIATE(GDDEVF=>YDML_PHY_MF%YRPHY0%GDDEVF, GCVALMX=>YDML_PHY_MF%YRPHY0%GCVALMX, &
 & TSPHY=>YDML_PHY_MF%YRPHY2%TSPHY, &
 & LCDDEVPRO=>YDML_PHY_MF%YRPHY%LCDDEVPRO, LNSDO=>YDML_PHY_MF%YRPHY%LNSDO, &
 & LGRAPRO=>YDML_PHY_MF%YRPHY%LGRAPRO)
! ----------------------------------------------------------------------

! *
! ------------------------------------------------------------------
!     I - INITIALISATIONS


! EFFECT OF GDDEVF on INITIAL DD T-PROFILE:
IF (LCDDEVPRO) THEN
   ZDDEVC=1.0_JPRB-GDDEVF
ELSEIF (LNSDO) THEN
   ZDDEVC=1._JPRB
ELSE
   ZDDEVC=0.0_JPRB
ENDIF

ZFCQRNG(:,:)=0.0_JPRB
ZFCQSNG(:,:)=0.0_JPRB
IF (LGRAPRO) THEN
 ZFCQGNG(:,:)=0.0_JPRB
ENDIF
ZFCQVNG(:,:)=0.0_JPRB
ZFCQLDM(:,:)=0.0_JPRB
ZFCQIDM(:,:)=0.0_JPRB
  ! ----------------------------------------------
  !  UPDATE VARIABLES BY MICROPHYSICS CONTRIBUTION
  ! ----------------------------------------------

!cdir unroll=8
    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA
        ZQX1=PZQR(JLON,JLEV)-PIPOI(JLON,JLEV)*(0.0_JPRB&
          & -(PFPFPSL(JLON,JLEV) - PFPFPSL(JLON,JLEV-1))&
          & +(PFPLSL(JLON,JLEV) - PFPLSL(JLON,JLEV-1))&
          & +(PFPEVPSL(JLON,JLEV)-PFPEVPSL(JLON,JLEV-1)) )
        PZQR(JLON,JLEV)=MAX(0.0_JPRB,ZQX1)
        ZDQR=MAX(0.0_JPRB,ZQX1)-ZQX1
        ZFCQRNG(JLON,JLEV)=ZFCQRNG(JLON,JLEV-1)-ZDQR*PPOID(JLON,JLEV)
        PFCQRNG(JLON,JLEV)=PFCQRNG(JLON,JLEV)+ZFCQRNG(JLON,JLEV)
        ZQX1=PZQS(JLON,JLEV)-PIPOI(JLON,JLEV)*(0.0_JPRB&
          & -(PFPFPSN(JLON,JLEV) - PFPFPSN(JLON,JLEV-1))&
          & +(PFPLSN(JLON,JLEV) - PFPLSN(JLON,JLEV-1))&
          & +(PFPEVPSN(JLON,JLEV)-PFPEVPSN(JLON,JLEV-1)) )
        PZQS(JLON,JLEV)=MAX(0.0_JPRB,ZQX1)
        ZDQS=MAX(0.0_JPRB,ZQX1)-ZQX1
        ZFCQSNG(JLON,JLEV)=ZFCQSNG(JLON,JLEV-1)-ZDQS*PPOID(JLON,JLEV)
        PFCQSNG(JLON,JLEV)=PFCQSNG(JLON,JLEV)+ZFCQSNG(JLON,JLEV)
        IF (LGRAPRO) THEN
          ZQX1=PZQG(JLON,JLEV)-PIPOI(JLON,JLEV)*(0.0_JPRB  &
            & -(PFPFPSG(JLON,JLEV) - PFPFPSG(JLON,JLEV-1)) &
            & +(PFPLSG(JLON,JLEV) - PFPLSG(JLON,JLEV-1))   &
            & +(PFPEVPSG(JLON,JLEV) - PFPEVPSG(JLON,JLEV-1)) )
          PZQG(JLON,JLEV)=MAX(0.0_JPRB,ZQX1)
          ZDQG=MAX(0.0_JPRB,ZQX1)-ZQX1
          ZFCQGNG(JLON,JLEV)=ZFCQGNG(JLON,JLEV-1)-ZDQG*PPOID(JLON,JLEV)
          PFCQGNG(JLON,JLEV)=PFCQGNG(JLON,JLEV)+ZFCQGNG(JLON,JLEV)
          ZDQC=ZDQR+ZDQS+ZDQG
        ELSE
          ZDQC=ZDQR+ZDQS
        ENDIF

        ZDFCQL=PZFCQL(JLON,JLEV)-PZFCQL(JLON,JLEV-1)&
          & - PFCSQL(JLON,JLEV)+PFCSQL(JLON,JLEV-1)&
          & - PFCCQL(JLON,JLEV)+PFCCQL(JLON,JLEV-1)
        ZFCQLDM(JLON,JLEV)=ZFCQLDM(JLON,JLEV-1)+ZDFCQL
        ZDFCQI=PZFCQI(JLON,JLEV)-PZFCQI(JLON,JLEV-1)&
          & - PFCSQN(JLON,JLEV)+PFCSQN(JLON,JLEV-1)&
          & - PFCCQN(JLON,JLEV)+PFCCQN(JLON,JLEV-1)
        ZFCQIDM(JLON,JLEV)=ZFCQIDM(JLON,JLEV-1)+ZDFCQI
        ZQX1=PZQV(JLON,JLEV)-PIPOI(JLON,JLEV)*(0.0_JPRB&
          & -(PFPEVPSL(JLON,JLEV)-PFPEVPSL(JLON,JLEV-1))&
          & -(PFPEVPSN(JLON,JLEV)-PFPEVPSN(JLON,JLEV-1))&
          & +(ZDFCQL+ZDFCQI) )
        IF (LGRAPRO) THEN
          ZQX1=ZQX1-PIPOI(JLON,JLEV)*(0.0_JPRB  &
            & -(PFPEVPSG(JLON,JLEV)-PFPEVPSG(JLON,JLEV-1)))
        ENDIF
        ZQV0=ZQX1-PIPOI(JLON,JLEV)*(0.0_JPRB&
          & -ZFCQVNG(JLON,JLEV-1)-ZFCQRNG(JLON,JLEV-1)&
          & -ZFCQSNG(JLON,JLEV-1))
        IF (LGRAPRO) THEN
          ZQV0=ZQV0-PIPOI(JLON,JLEV)*(0.0_JPRB - ZFCQGNG(JLON,JLEV-1))
        ENDIF
        PZQV(JLON,JLEV)=MAX(0.0_JPRB,ZQV0-ZDQC)
        ZDQV=MAX(0.0_JPRB,ZQV0-ZDQC)-ZQX1
        ZFCQVNG(JLON,JLEV)=ZFCQVNG(JLON,JLEV-1)-ZDQV*PPOID(JLON,JLEV)
        PFCQVNG(JLON,JLEV)= PFCQVNG(JLON,JLEV)+ZFCQVNG(JLON,JLEV)
        PZQL(JLON,JLEV)=PZQL(JLON,JLEV)-PIPOI(JLON,JLEV)*(&
          & (PFPFPSL(JLON,JLEV)-PFPFPSL(JLON,JLEV-1))&
          & -ZDFCQL )
        PZQI(JLON,JLEV)=PZQI(JLON,JLEV)-PIPOI(JLON,JLEV)*(&
          & (PFPFPSN(JLON,JLEV)-PFPFPSN(JLON,JLEV-1))&
          & -ZDFCQI )

        PZT(JLON,JLEV)=PZT(JLON,JLEV)-PIPOI(JLON,JLEV)/PCP(JLON,JLEV)&
         & *(0.0_JPRB-PLHV(JLON,JLEV)*ZDFCQL-PLHS(JLON,JLEV)*ZDFCQI)

      ENDDO
    ENDDO

    !   STRATIFORM CORRECTION OF THE TOTAL CONDENSATION FLUX
    !   ------------------------------------

    DO JLEV=KTDIA,KLEV
      DO JLON=KIDIA,KFDIA
        PFCSQL(JLON,JLEV)=PFCSQL(JLON,JLEV)+ZFCQLDM(JLON,JLEV)
        PFCSQN(JLON,JLEV)=PFCSQN(JLON,JLEV)+ZFCQIDM(JLON,JLEV)
      ENDDO
    ENDDO
! ------------------------------------------------------------------
!     II - CALCULATIONS

PFHP(KIDIA:KFDIA,KTDIA-1)=0._JPRB
DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    !
    !  ------------------------------------------------------------
    !   COMPUTE THE TEMPERATURE SINK SOURCE FOR DOWNDRAUGHT CLOSURE
    !  ------------------------------------------------------------
    !
    ZDEVPL=PFPEVPSL(JLON,JLEV)-PFPEVPSL(JLON,JLEV-1)
    ZDEVPN=PFPEVPSN(JLON,JLEV)-PFPEVPSN(JLON,JLEV-1)
    PFHP(JLON,JLEV)= PFHP(JLON,JLEV-1) + (&
      &  PLHV(JLON,JLEV)*ZDEVPL +PLHS(JLON,JLEV)*ZDEVPN )
    IF (LGRAPRO) THEN
      PFHP(JLON,JLEV)=PFHP(JLON,JLEV) + ( &
        & PLHS(JLON,JLEV)*(PFPEVPSG(JLON,JLEV)-PFPEVPSG(JLON,JLEV-1)) )
    ENDIF
! Effect of 1-gddevf on profile
! if PFHP is bigger at level l than l-1, there is residual cooling of layer l
    PZT(JLON,JLEV)= PZT(JLON,JLEV)-ZDDEVC*(PIPOI(JLON,JLEV)/PCP(JLON,JLEV))*&
      &(PFHP(JLON,JLEV)-PFHP(JLON,JLEV-1))
    IF (LNSDO)&
      & PZQV(JLON,JLEV)=PZQV(JLON,JLEV)+ZDDEVC*PIPOI(JLON,JLEV)*(ZDEVPL+ZDEVPN)
    !
    ! Adapt values to updraught environment for ACMODO
    !
    ZUDAL=MIN(GCVALMX,PUDAL(JLON,JLEV))
! Mean Updraught environment vertical velocity:
    POME(JLON,JLEV)=(PEVEL0(JLON,JLEV)&
         &   -ZUDAL*PUDOM(JLON,JLEV))*TSPHY
  ENDDO
ENDDO

! ----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACUPM',1,ZHOOK_HANDLE)
END SUBROUTINE ACUPM
