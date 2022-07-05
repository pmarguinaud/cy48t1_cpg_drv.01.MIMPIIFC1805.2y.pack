!OCL  NOEVAL
SUBROUTINE GNH_TNDLAGADIAB_GW(&
 ! --- INPUT -----------------------------------------------------------------
 & YDCST, YDGEOMETRY, KFLEV,KPROMA,KSTART,KPROF,&
 ! --- OPTIONAL INPUT --------------------------------------------------------
 & POROGL,POROGM,POROGLM,POROGLL,POROGMM,&
 & PUS,PVS,PTNDUS,PTNDVS,&
 & PDEP,PDEPS,PREF,PREH,PDELP,PRDELP,PDELNHPRE,PEVEL,&
 ! --- OPTIONAL OUTPUT -------------------------------------------------------
 & PTNDGWH_LAP,PTNDGWH_OTH,&
 & PTNDGWF_LAP,PTNDGWF_OTH,PDBBC&
 & )  

!**** *GNH_TNDLAGADIAB_GW* - Computation of the adiabatic contribution of
!                            [D (Gw) / Dt] in the NH model.

!     Purpose.
!     --------

!     Computation of the adiabatic contribution of [D (Gw) / Dt] in the NH model.
!     Upper air and surface [D (Gw) / Dt] are computed.

!     Upper air [D (Gw) / Dt] for thin layer equations:
!      [D (Gw) / Dt]_adiab = g**2 ( (d pre/ d prehyd) - 1 )

!     Surface [D (Gw) / Dt] for thin layer equations:
!      a temporal derivation of the lower boundary condition
!      [Gw]_surf = V_surf grad[Phi_surf]
!      is done; that yields:
!      [D [Gw]_surf / Dt]_adiab = [D V_surf / Dt] grad[Phi_surf]
!                                 + V_surf . [V_surf . grad[grad[Phi_surf]]]
!     Calculations at the surface are useless for LGWADV=T.

!     For the finite difference vertical discretisation, the upper air
!      [D (Gw) / Dt] is computed at half levels.

!     For the finite element vertical discretisation, the upper air
!      [D (Gw) / Dt] is computed at full levels.

!     There is a separate treatment for the Laplacian part and the non-Laplacian part
!     (the non-Laplacian part is generally zero in the thin layer system of eqns).

!**   Interface.
!     ----------
!        *CALL* *GNH_TNDLAGADIAB_GW(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!           KFLEV     : number of levels.
!           KPROMA    : horizontal dimension.
!           KSTART    : start of work.
!           KPROF     : working length.

!         * OPTIONAL INPUT:
!           POROGL    : zonal comp. of grad(surf orography).
!           POROGM    : merid comp. of grad(surf orography).
!           POROGLM   : -I
!           POROGLL   :  I- second order derivatives of "surf orography".
!           POROGMM   : -I
!           PUS       : surface U wind.
!           PVS       : surface V wind.
!           PTNDUS    : [DU/Dt]_surf.
!           PTNDVS    : [DV/Dt]_surf.
!           PDEP      : (pre - prehyd) at full levels.
!           PDEPS     : surface (pre - prehyd).
!           PREF      : "prehyd" at full levels.
!           PREH      : "prehyd" at half levels.
!           PDELP     : "Delta prehyd" at full levels.
!           PRDELP    : 1/[Delta prehyd] at full levels.
!           PDELNHPRE : "Delta pre" at full levels.
!           PEVEL     : "etadot dprehyd/deta".

!         * OPTIONAL OUTPUT:
!           Partial contributions are kept separate to allow individual VFE/FD
!           treatment for each term (when possible).
!           PTNDGWH_LAP: [D (Gw) / Dt]_adiab at half levels (Laplacian part).
!           PTNDGWH_OTH: [D (Gw) / Dt]_adiab at half levels (other contribs).
!           PTNDGWF_LAP: [D (Gw) / Dt]_adiab at full levels (Laplacian part).
!           PTNDGWF_OTH: [D (Gw) / Dt]_adiab at full levels (other contribs).
!           PDBBC      : [D (Gw)_surf / Dt]_adiab.

!        Implicit arguments :   None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    None.
!     ----------

!     Reference.
!     ----------
!        See documentations Arpege about the NH model.

!     Author.
!     -------
!        K. YESSAD (after routine GNHPDVD).
!        Original : 08-Dec-2004

!     Modifications.
!     --------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Dec 2009): LRWSDLW,LRWSDLR,LRWSDLG=T,T,T in NH model for LGWADV=F.
!   K. Yessad (Jan 2011): bug correction for (LRWSDLR,LGWADV)=(T,F).
!   P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!   K. Yessad (Dec 2016): Prune obsolete options.
!   K. Yessad (June 2017): Introduce NHQE model.
!   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   J. Vivoda and P. Smolikova (Sep 2020): VFE pruning.
!   H Petithomme (Dec 2020): test cleaning, simplification
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : TCST
USE YOMDYNA      , ONLY : LGWADV
USE YOMCVER      , ONLY : LVERTFE, LVFE_GW, LVFE_LAPL_BC, LVFE_ECMWF

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: POROGL(KPROMA)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: POROGM(KPROMA)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: POROGLM(KPROMA)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: POROGLL(KPROMA)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: POROGMM(KPROMA)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PUS(KPROMA)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PVS(KPROMA)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PTNDUS(KPROMA)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PTNDVS(KPROMA)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PDEP(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PDEPS(KPROMA)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PREF(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PREH(KPROMA,0:KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PDELP(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PRDELP(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PDELNHPRE(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(IN)    :: PEVEL(KPROMA,0:KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)   :: PTNDGWH_LAP(KPROMA,0:KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)   :: PTNDGWH_OTH(KPROMA,0:KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)   :: PTNDGWF_LAP(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)   :: PTNDGWF_OTH(KPROMA,KFLEV)
REAL(KIND=JPRB),OPTIONAL,INTENT(OUT)   :: PDBBC(KPROMA)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,JLEV
REAL(KIND=JPRB) :: ZRG2,ZJACOB,ZACUGF,ZACVGF
REAL(KIND=JPRB) :: ZIN (KPROMA,0:KFLEV+1)
REAL(KIND=JPRB) :: ZOUT(KPROMA,0:KFLEV  )

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

#include "abor1.intfb.h"
#include "verdisint.intfb.h"

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNH_TNDLAGADIAB_GW',0,ZHOOK_HANDLE)

ZRG2=YDCST%RG*YDCST%RG

!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF [D [Gw]_surf / Dt]_adiab
!              ---------------------------------------

IF (.NOT.LGWADV) THEN
  IF (.NOT.(PRESENT(PDBBC))) &
   & CALL ABOR1(' GNH_TNDLAGADIAB_GW 3a: missing optional output arguments!')

  IF (.NOT.(PRESENT(POROGL).AND.PRESENT(POROGM).AND.PRESENT(POROGLL) &
   & .AND.PRESENT(POROGMM).AND.PRESENT(POROGLM).AND.PRESENT(PUS) &
   & .AND.PRESENT(PVS).AND.PRESENT(PTNDUS).AND.PRESENT(PTNDVS))) &
   & CALL ABOR1(' GNH_TNDLAGADIAB_GW 1a: missing optional input arguments!')

  DO JROF=KSTART,KPROF
    ZJACOB=POROGLL(JROF)*PUS(JROF)*PUS(JROF) &
     & + POROGMM(JROF)*PVS(JROF)*PVS(JROF) &
     & + 2.0_JPRB*POROGLM(JROF)*PUS(JROF)*PVS(JROF)
    ZACUGF=PTNDUS(JROF)*POROGL(JROF)
    ZACVGF=PTNDVS(JROF)*POROGM(JROF)
    PDBBC(JROF)=-ZACUGF-ZACVGF-ZJACOB
  ENDDO

!     ------------------------------------------------------------------

!*       2.    COMPUTATION OF UPPER AIR [D (Gw) / Dt]_adiab
!              --------------------------------------------

  ! ** Treatment of the Laplacian term.

  IF (LVERTFE.AND..NOT.LVFE_LAPL_BC) THEN
    ! * PTNDGWH_LAP is not computed without LVFE_LAPL_BC

    IF (.NOT.(PRESENT(PDEP).AND.PRESENT(PREF).AND.PRESENT(PREH))) &
     & CALL ABOR1(' GNH_TNDLAGADIAB_GW 2Ab: missing optional input arguments!')
    IF (.NOT.(PRESENT(PTNDGWH_LAP))) &
     & CALL ABOR1(' GNH_TNDLAGADIAB_GW 2Ac: missing optional output arguments!')

    ! * PTNDGWH_LAP is not computed (excepted at half levels 0 and 1);

    ! * Half levels 1 and KFLEV-1.
    DO JLEV=1,KFLEV-1,KFLEV-2
      DO JROF=KSTART,KPROF
        PTNDGWH_LAP(JROF,JLEV)=ZRG2*(PDEP(JROF,JLEV+1)-PDEP(JROF,JLEV)) &
         & /(PREF(JROF,JLEV+1)-PREF(JROF,JLEV))
      ENDDO
    ENDDO

    ! * Top of the model.
    !   [d (pre-prehyd) / d prehyd]_[l=1] is computed by
    !   [ (pre-prehyd)_[l=1] ]/[prehyd_[l=1] - prehyd_[top]]
    !   and (pre-prehyd)_[top] is assumed to be equal to zero.
    DO JROF=KSTART,KPROF
      PTNDGWH_LAP(JROF,0)=ZRG2*PDEP(JROF,1)/(PREF(JROF,1)-PREH(JROF,0))
    ENDDO

    ! * Surface.
    DO JROF=KSTART,KPROF
      PTNDGWH_LAP(JROF,KFLEV)=-PDBBC(JROF)
    ENDDO

  ELSEIF (.NOT.LVERTFE) THEN

    IF (.NOT.(PRESENT(PDEP).AND.PRESENT(PREF).AND.PRESENT(PREH))) &
     & CALL ABOR1(' GNH_TNDLAGADIAB_GW 2Ad: missing optional input arguments!')
    IF (.NOT.(PRESENT(PTNDGWH_LAP))) &
     & CALL ABOR1(' GNH_TNDLAGADIAB_GW 2Ae: missing optional output arguments!')

    ! * Half levels 1 to nflevg-1.
    DO JLEV=1,KFLEV-1
      DO JROF=KSTART,KPROF
        PTNDGWH_LAP(JROF,JLEV)=ZRG2*(PDEP(JROF,JLEV+1)-PDEP(JROF,JLEV)) &
         & /(PREF(JROF,JLEV+1)-PREF(JROF,JLEV))
      ENDDO
    ENDDO

    ! * Top of the model.
    !   [d (pre-prehyd) / d prehyd]_[l=1] is computed by
    !   [ (pre-prehyd)_[l=1] ]/[prehyd_[l=1] - prehyd_[top]]
    !   and (pre-prehyd)_[top] is assumed to be equal to zero.
    DO JROF=KSTART,KPROF
      PTNDGWH_LAP(JROF,0)=ZRG2*PDEP(JROF,1)/(PREF(JROF,1)-PREH(JROF,0))
    ENDDO

    ! * Surface.
    DO JROF=KSTART,KPROF
      PTNDGWH_LAP(JROF,KFLEV)=-PDBBC(JROF)
    ENDDO

  ENDIF

  ! ** Treatment of the other terms: zero.

  IF (LVFE_GW) THEN
    IF (.NOT.PRESENT(PTNDGWF_OTH)) &
     & CALL ABOR1(' GNH_TNDLAGADIAB_GW 2Af: missing optional output arguments!')
    PTNDGWF_OTH(KSTART:KPROF,1:KFLEV)=0.0_JPRB
  ELSE
    IF (.NOT.PRESENT(PTNDGWH_OTH)) &
     & CALL ABOR1(' GNH_TNDLAGADIAB_GW 2Ag: missing optional output arguments!')
    PTNDGWH_OTH(KSTART:KPROF,0:KFLEV)=0.0_JPRB
  ENDIF

ELSE

  ! ** Treatment of the Laplacian term, and fill PTNDGW[F,H]_OTH with zeros.

  IF (LVFE_GW) THEN

    IF (.NOT.(PRESENT(PDELNHPRE).AND.PRESENT(PDELP))) &
     & CALL ABOR1(' GNH_TNDLAGADIAB_GW 2Ba: missing optional input arguments!')
    IF (.NOT.(PRESENT(PTNDGWF_LAP).AND.PRESENT(PTNDGWF_OTH))) &
     & CALL ABOR1(' GNH_TNDLAGADIAB_GW 2Bb: missing optional output arguments!')

    DO JLEV=1,KFLEV
      DO JROF=KSTART,KPROF
        PTNDGWF_LAP(JROF,JLEV)=ZRG2*(PDELNHPRE(JROF,JLEV)/PDELP(JROF,JLEV) - 1.0_JPRB)
      ENDDO
    ENDDO
    PTNDGWF_OTH(KSTART:KPROF,1:KFLEV)=0.0_JPRB

  ELSE

    IF (.NOT.(PRESENT(PDEP).AND.PRESENT(PREF).AND.PRESENT(PREH))) &
     & CALL ABOR1(' GNH_TNDLAGADIAB_GW 2Bd: missing optional input arguments!')
    IF (.NOT.(PRESENT(PTNDGWH_LAP).AND.PRESENT(PTNDGWH_OTH))) &
     & CALL ABOR1(' GNH_TNDLAGADIAB_GW 2Be: missing optional output arguments!')

    ! * Half levels 1 to nflevg-1.

    IF(LVERTFE .AND. .NOT.LVFE_ECMWF)THEN

      ! pressure departure is derivated using VFE operator from full levels to half levels
      ! but still through transformation from d4->gw and back is FD style
      ZIN(KSTART:KPROF,0)=0.0_JPRB
      ZIN(KSTART:KPROF,1:KFLEV)=PDEP(KSTART:KPROF,1:KFLEV)
      ZIN(KSTART:KPROF,KFLEV+1)=0.0_JPRB
      CALL VERDISINT(YDGEOMETRY%YRVFE,'HDER','01',KPROMA,KSTART,KPROF,KFLEV,ZIN,ZOUT)

      DO JLEV=1,KFLEV-1
        DO JROF=KSTART,KPROF
          PTNDGWH_LAP(JROF,JLEV)=ZRG2*ZOUT(JROF,JLEV)*(YDGEOMETRY%YRVETA%VFE_ETAF(JLEV+1) &
           & - YDGEOMETRY%YRVETA%VFE_ETAF(JLEV)) / (PREF(JROF,JLEV+1)-PREF(JROF,JLEV))
        ENDDO
      ENDDO

    ELSE
      DO JLEV=1,KFLEV-1
        DO JROF=KSTART,KPROF
          PTNDGWH_LAP(JROF,JLEV)=ZRG2*(PDEP(JROF,JLEV+1)-PDEP(JROF,JLEV)) &
           & /(PREF(JROF,JLEV+1)-PREF(JROF,JLEV))
        ENDDO
      ENDDO
    ENDIF

    PTNDGWH_OTH(KSTART:KPROF,1:KFLEV-1)=0.0_JPRB

    ! * Top of the model.
    !   [d (pre-prehyd) / d prehyd]_[top] is computed by
    !   [ (pre-prehyd)_[l=1] ]/[prehyd_[l=1] - prehyd_[top]]
    !   and (pre-prehyd)_[top] is assumed to be equal to zero.
    DO JROF=KSTART,KPROF
      PTNDGWH_LAP(JROF,0)=ZRG2*PDEP(JROF,1)/(PREF(JROF,1)-PREH(JROF,0))
    ENDDO
    PTNDGWH_OTH(KSTART:KPROF,0)=0.0_JPRB
  ENDIF
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNH_TNDLAGADIAB_GW',1,ZHOOK_HANDLE)
END SUBROUTINE GNH_TNDLAGADIAB_GW
