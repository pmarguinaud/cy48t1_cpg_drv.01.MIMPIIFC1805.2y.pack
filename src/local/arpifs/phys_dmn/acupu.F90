!OPTIONS XOPT(NOEVAL)
SUBROUTINE ACUPU( YDPHY,YDPHY0,KIDIA,KFDIA,KLON,KTDIA,KLEV,&
! ----------------------------------------------------------------------
! - INPUT  2D .
                &PUDAL,PUNEBH,PDETFI,&
                &PDPOID,PIPOI,PLHV,PLHS,PCP,&
                &PDIFCQ,PDIFCQI,PDIFCQL,PDIFCS,&
                &PFCCQL,PFCCQN,PFCSQL,PFCSQN, &
! - OUTPUT  1D .
                &PSIGP, PSIGPC,&
! - INPUT/OUTPUT  2D .
                &PNEBE,PZQI,PZQL,PZQV,PZT,&
                &PFCQVNG,PFCQING,PFCQLNG)
! **** *ACUPU * - UPDATE INTERNAL STATE AFTER UPDRAUGHT
! **   Interface.
!      ----------
! *CALL* *ACUPU*
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
! PUDAL:  prognostic updraught mesh fraction
! PDPOID: DP/(RG*DT) weight.
! PIPOI:  Inverse of DP/(RG*DT).
! PLHV :  latent heat evaporation.
! PLHS :  latent heat sublimation.
! PCP  :  specific heat at constant pressure.
 
! - 2D FLUXES (0:KLEV) .
! PDIFCQ  : water vapour transport by updraft.
! PDIFCQI : cloud ice transport by updraft.
! PDIFCQL : cloud water transport by updraft.
! PDIFCS  : enthalpy transport by updraft.
! PFCCQL
! PFCCQN
! PFCSQL
! PFCSQN

! ----------------------------------------------------------------------
! -   OUTPUT ARGUMENTS (I/O because initialized to zero in aplpar)
!     ----------------
! - 1D
! PSIGP     : PRECIPITATION MESH FRACTION
! PSIGPC    : CONVECTIVE FRACTION OF THE PRECIPITATION FLUX
! - 2D
! ----------------------------------------------------------------------
! -   INPUT/OUTPUT ARGUMENTS (I/O because update here)
!     ----------------------
! - 2D (1:KLEV)
! PUNEBH : pseudo-historic detrainment mesh fraction
! PDETFI : instantaneous detrainment mesh fraction
! PNEBE  : input=strat cloud fraction, output=equiv cloud fraction.
! PZQI :  cloud ice of the cascade.
! PZQL :  cloud water of the cascade.
! PZQV :  water vapour of the cascade.
! PZT  :  temperature of the cascade.
 
! - 2D FLUXES (0:KLEV)
! PFCQING : cloud ice negative values
! PFCQLNG : cloud water negative values
! PFCQVNG : water vapour negative values
! ----------------------------------------------------------------------
!     Author.
!     -------
!     01-2007, L. Gerard 

!     Modifications.
!     --------------
!     04-2007, L. Gerard 
!     K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!     R. Brozkova (Jan 2013): moving LNEBINS to module.
!     L. Gerard (Apr 2016):  New formulation of equiv cloud fraction for LCVCSD;
!       LLFIX: more physical expression of  mutually exclusive conv and strat 
!       cloud fractions; reorganize vertical loop.
! ----------------------------------------------------------------------

USE PARKIND1 , ONLY: JPIM, JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE YOMPHY   , ONLY : TPHY
USE YOMPHY0  , ONLY : TPHY0

! ----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TPHY)         ,INTENT(IN) :: YDPHY
TYPE(TPHY0)        ,INTENT(IN) :: YDPHY0
INTEGER(KIND=JPIM), INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
INTEGER(KIND=JPIM), INTENT(IN) :: KLON
INTEGER(KIND=JPIM), INTENT(IN) :: KTDIA
REAL(KIND=JPRB) ,INTENT(IN) :: PUDAL(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PDPOID(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PIPOI(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PLHV(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PLHS(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PCP(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PDIFCQ(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PDIFCQI(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PDIFCQL(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PDIFCS(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PFCCQL(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PFCCQN(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PFCSQL(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(IN) :: PFCSQN(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(OUT) :: PSIGP(KLON)
REAL(KIND=JPRB) ,INTENT(OUT) :: PSIGPC(KLON)
REAL(KIND=JPRB) ,INTENT(IN)  :: PUNEBH(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(IN)  :: PDETFI(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PNEBE(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PZQI(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PZQL(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PZQV(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PZT(KLON,KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PFCQING(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PFCQLNG(KLON,0:KLEV)
REAL(KIND=JPRB) ,INTENT(INOUT) :: PFCQVNG(KLON,0:KLEV)

! ----------------------------------------------------------------------

INTEGER (KIND=JPIM) :: JLEV, JLON
REAL (KIND=JPRB) :: ZFCQING(KLON,0:KLEV)
REAL (KIND=JPRB) :: ZFCQLNG(KLON,0:KLEV)
REAL (KIND=JPRB) :: ZFCQVNG(KLON,0:KLEV)
REAL (KIND=JPRB) :: ZPX(KLON), ZCC1(KLON), ZCS1(KLON)
REAL (KIND=JPRB) :: ZNEBC(KLON),ZFRCO(KLON)
REAL(KIND=JPRB)  :: ZEPS0, ZNEI,ZNEQ,&
                 &  ZZ,ZZ1,ZZ2, ZCC,ZCS,ZNC,ZNT,ZDQIT,ZDQLT,&
                 &  ZQX0, ZQX1, ZDQI, ZDQL, ZDQV, ZDQC, ZQV0
LOGICAL :: LLFIX
REAL(KIND=JPRB) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ACUPU',0,ZHOOK_HANDLE)
ASSOCIATE(GRRMINA=>YDPHY0%GRRMINA, &
 & LCVCSD=>YDPHY%LCVCSD, LNEBINS=>YDPHY%LNEBINS)
! ----------------------------------------------------------------------

! *
! ------------------------------------------------------------------
!     I - INITIALISATIONS

LLFIX=LCVCSD
ZEPS0=1.E-12_JPRB

ZFCQING=0._JPRB
ZFCQLNG=0._JPRB
ZFCQVNG=0._JPRB

    !  -------------------------------------------
    !  UPDATE VARIABLES BY UPDRAUGHT CONTRIBUTION
    !  -------------------------------------------
!cdir unroll=8
DO JLEV=KTDIA,KLEV
  DO JLON=KIDIA,KFDIA
    ZDQIT=PIPOI(JLON,JLEV)*(0.0_JPRB&
      & -(PFCCQN(JLON,JLEV)-PFCCQN(JLON,JLEV-1)) )
    ZQX0=PZQI(JLON,JLEV)
    ZQX1=ZQX0-ZDQIT-PIPOI(JLON,JLEV)*&
      &  (PDIFCQI(JLON,JLEV)-PDIFCQI(JLON,JLEV-1))
    PZQI(JLON,JLEV)=MAX(0.0_JPRB,ZQX1)
    ZDQI=MAX(0.0_JPRB,ZQX1)-ZQX1
    ZFCQING(JLON,JLEV)=ZFCQING(JLON,JLEV-1)-ZDQI*PDPOID(JLON,JLEV)
    PFCQING(JLON,JLEV)=PFCQING(JLON,JLEV)+ZFCQING(JLON,JLEV)

    ZDQLT=PIPOI(JLON,JLEV)*(0.0_JPRB&
      & -(PFCCQL(JLON,JLEV)-PFCCQL(JLON,JLEV-1)) )
    ZQX0=PZQL(JLON,JLEV)
    ZQX1=ZQX0-ZDQLT-PIPOI(JLON,JLEV)*&
      &  (PDIFCQL(JLON,JLEV)-PDIFCQL(JLON,JLEV-1))
    PZQL(JLON,JLEV)=MAX(0.0_JPRB,ZQX1)
    ZDQL=MAX(0.0_JPRB,ZQX1)-ZQX1
    ZFCQLNG(JLON,JLEV)=ZFCQLNG(JLON,JLEV-1)-ZDQL*PDPOID(JLON,JLEV)
    PFCQLNG(JLON,JLEV)=PFCQLNG(JLON,JLEV)+ZFCQLNG(JLON,JLEV)
    ZDQC=ZDQI+ZDQL

    ZQX0=PZQV(JLON,JLEV)
    ZQX1=ZQX0-PIPOI(JLON,JLEV)*(&
      &  (PDIFCQ(JLON,JLEV)-PDIFCQ(JLON,JLEV-1))&
      &  +(PFCCQN(JLON,JLEV)-PFCCQN(JLON,JLEV-1))&
      & +(PFCCQL(JLON,JLEV)-PFCCQL(JLON,JLEV-1)) )
    ZQV0=ZQX1-PIPOI(JLON,JLEV)*(0.0_JPRB&
      & -ZFCQVNG(JLON,JLEV-1)-ZFCQING(JLON,JLEV-1)&
      & -ZFCQLNG(JLON,JLEV-1))
    PZQV(JLON,JLEV)=MAX(0.0_JPRB,ZQV0-ZDQC)
    ZDQV=MAX(0.0_JPRB,ZQV0-ZDQC)-ZQX1
    ZFCQVNG(JLON,JLEV)=ZFCQVNG(JLON,JLEV-1)-ZDQV*PDPOID(JLON,JLEV)
    PFCQVNG(JLON,JLEV)=PFCQVNG(JLON,JLEV)+ZFCQVNG(JLON,JLEV)

    PZT(JLON,JLEV)=PZT(JLON,JLEV)-(PIPOI(JLON,JLEV)*&
      &  (PDIFCS(JLON,JLEV)-PDIFCS(JLON,JLEV-1))&
      &  +PLHS(JLON,JLEV)*ZDQIT+PLHV(JLON,JLEV)*ZDQLT)/PCP(JLON,JLEV)
  ENDDO
ENDDO

! ------------------------------------------------------------------
!    II - CALCULATIONS
!
! REM: AT ENTRY PUDAL is guaranteed to be in [0, GCVALMX]

ZPX=0.0_JPRB
PSIGP(KIDIA:KFDIA)=1.0_JPRB
ZCC1=0.0_JPRB
ZCS1=0.0_JPRB
DO JLEV=KTDIA,KLEV
  IF (LNEBINS) THEN
     DO JLON=KIDIA,KFDIA
        ZNEBC(JLON)=PDETFI(JLON,JLEV)+PUDAL(JLON,JLEV)
     ENDDO
  ELSE
     DO JLON=KIDIA,KFDIA
        ZNEBC(JLON)=PUNEBH(JLON,JLEV)+PUDAL(JLON,JLEV)
     ENDDO
  ENDIF
  DO JLON=KIDIA,KFDIA
! ZFRCO = Convective fraction of the condensate at each level
!         based on the ratio of *local* CONDENSATION
     ZCC= PFCCQL(JLON,JLEV)+PFCCQN(JLON,JLEV)
     ZCS= PFCSQL(JLON,JLEV)+PFCSQN(JLON,JLEV)
     ZNC=MAX(0.0_JPRB,ZCC-ZCC1(JLON))
     ZNT=ZNC+MAX(0.0_JPRB,ZCS-ZCS1(JLON))
     ZCC1(JLON)=ZCC
     ZCS1(JLON)=ZCS
! ZFRCO is reset to zero if total cond or ZNEBC is too small
     ZZ2=MAX(0.0_JPRB,SIGN(1.0_JPRB,ZNT-ZEPS0))&
        & *MAX(0.0_JPRB,SIGN(1.0_JPRB,ZNEBC(JLON)-GRRMINA))
     ZFRCO(JLON)=ZZ2*ZNC/(ZNT+(1.0_JPRB-ZZ2))
  ENDDO
  IF(LLFIX) THEN
   DO JLON=KIDIA,KFDIA
! STRATIFORM FRACTION: (1-n_c)*PNEBE
! THIS ENSURES THAT THEIR SUM IS THE (UPDATED) TOTAL CLOUDINESS <=1
     PNEBE(JLON,JLEV)=PNEBE(JLON,JLEV)*(1._JPRB-ZNEBC(JLON))
     ZPX(JLON)=MAX(ZPX(JLON),&
                  & PNEBE(JLON,JLEV)+ZNEBC(JLON)-PUDAL(JLON,JLEV))
   ENDDO
  ELSE
!----------
! Old way to re-scale (nc,ns)->(nc',ns') ensuring nc'+ns'<=1
! (conceptually wrong, kept for compatibility)
    DO JLON=KIDIA,KFDIA
        ZZ=MAX(0.0_JPRB, 1.0_JPRB-PUDAL(JLON,JLEV))
        PSIGP(JLON)=MIN(PSIGP(JLON),ZZ)
        ZZ1=ZNEBC(JLON)+PNEBE(JLON,JLEV)
        ZZ2=ZZ1-ZNEBC(JLON)*PNEBE(JLON,JLEV)
        ZPX(JLON)=MAX(ZPX(JLON),ZZ2)
        ZZ=MAX(0.0_JPRB,SIGN(1.0_JPRB,ZZ1-GRRMINA))
        ZNEBC(JLON)=ZZ*ZNEBC(JLON)*ZZ2/(ZZ1+(1.0_JPRB-ZZ))
        PNEBE(JLON,JLEV)=ZZ*PNEBE(JLON,JLEV)*ZZ2/(ZZ1+(1.0_JPRB-ZZ))
   ENDDO
  ENDIF
! COMPUTE EQUIVALENT CLOUD FRACTION FOR MICROPHYSICS
! ne: equivalent cloud fraction for microphysics
!     nc, ns are now disjunct (nc+ns=nt<=1)
  IF (LCVCSD) THEN
!-------------
! CSD FORMULATION, COMPATIBLE WITH SMALL ZFRCO DESPITE LARGE NC
! ne=(1-ZFRCO)*nt+nc or ne=nt (if nt too small or ZFRCO<nc/nt)
!-------------
     DO JLON=KIDIA,KFDIA
        ZNT=PNEBE(JLON,JLEV)+ZNEBC(JLON)
        ZZ=MAX(0._JPRB,SIGN(1._JPRB, ZNT-GRRMINA))
        PNEBE(JLON,JLEV)=ZNT*(1._JPRB-ZZ*MAX(0._JPRB,&
                     & ZFRCO(JLON)-ZNEBC(JLON)/(ZNT+1._JPRB-ZZ)))
     ENDDO
  ELSE
!-------------
! OLD FORMULATION:
! 1/ne= (1-ZFRCO)^2/ns + ZFRCO^2/nc
! if ns or nc too small, use the bigger of them
!--------------
     DO JLON=KIDIA,KFDIA
        ZZ1=MAX(0.0_JPRB,SIGN(1.0_JPRB,PNEBE(JLON,JLEV)-GRRMINA))
        ZZ2=MAX(0.0_JPRB,SIGN(1.0_JPRB,ZNEBC(JLON)-GRRMINA))
! case when ns' or nc' or both are too small: the smaller one or zero
        ZNEQ=ZZ1*(1.0_JPRB-ZZ2)*PNEBE(JLON,JLEV)&
             &+ZZ2*(1.0_JPRB-ZZ1)*ZNEBC(JLON)
! case when both nc' and ns' are big enough: combine them
        ZZ=(1.0_JPRB-ZFRCO(JLON))
! 1/ne is then (if zz1=0=zz2)
        ZNEI= ZZ*ZZ/(PNEBE(JLON,JLEV)+(1.0_JPRB-ZZ1))&
             &+ZFRCO(JLON)*ZFRCO(JLON)/(ZNEBC(JLON)+(1.0_JPRB-ZZ2)) 
!   \--> this must be greater than 1 !
! Put the result ne (equiv cloud fraction) in PNEBE
        PNEBE(JLON,JLEV)=ZNEQ + ZZ1*ZZ2/MAX(0.0_JPRB,ZNEI)
     ENDDO
  ENDIF ! NOT LCVCSD
ENDDO ! JLEV
! OUTPUT CONVECTIVE CONDENSATION RATIO and PRECIPITATION FRACTION
DO JLON=KIDIA, KFDIA
   ZZ= MAX(0.0_JPRB, PFCCQL(JLON,KLEV)+PFCCQN(JLON,KLEV))
   ZZ1=ZZ+MAX(0.0_JPRB, PFCSQL(JLON,KLEV)+PFCSQN(JLON,KLEV))
   ZZ2=MAX(0.0_JPRB,SIGN(1.0_JPRB,ZZ1-ZEPS0))
   PSIGPC(JLON)=ZZ2*ZZ/(MAX(ZZ1,ZEPS0))
! Precipitation area available for dd:
   PSIGP(JLON)=MIN(PSIGP(JLON), ZPX(JLON))
ENDDO

! ----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACUPU',1,ZHOOK_HANDLE)
END SUBROUTINE ACUPU
