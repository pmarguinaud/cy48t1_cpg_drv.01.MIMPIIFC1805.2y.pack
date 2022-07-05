!OCL  NOEVAL
SUBROUTINE GNH_TNDLAGADIAB_SVD(&
 ! --- INPUT -----------------------------------------------------------------
 & YDCST, YDGEOMETRY, KFLEV,KPROMA,KST,KPROF,&
 & PUF,PVF,PUH,PVH,PRDPHI,PGWH,PGWF,PGWHL,PGWHM,PGWFL,PGWFM,PDVER,PDIV,P3DIVG,&
 & PTNDGWH_LAP,PTNDGWH_OTH,PTNDGWF_OTH,&
 & PDEP,PREH,PREF,PDELP,PLNPR,&
 ! --- OUTPUT ----------------------------------------------------------------
 & PTNDVD&
 & )  

!**** *GNH_TNDLAGADIAB_SVD* - Computation of the adiabatic contribution of
!                             [D (svd) / Dt] in the NHEE model.
!                             "svd" is the vertical divergence variable.

!     Purpose.
!     --------

!     Computation of the adiabatic contribution of [D (svd) / Dt] in the NHEE model
!     at full levels.
!     For svd=d4 (NVDVAR=4) the lagrangian tendency currently computed
!     is [D (dver) / Dt]_adiab (and not [D (dqua) / Dt]_adiab).

!      For svd=dver, expression of this tendency is:
!       [D (dver) / Dt]_adiab =
!        dver (D-D3)
!        - [pre/(R T)] [ d [ [D (Gw) / Dt]_adiab ] / d prehyd ]
!        + [pre/(R T)] grad(Gw).[ d V / d prehyd ]
!      which can be rewritten:
!       [D (dver) / Dt]_adiab =
!        dver (D-D3)
!        - [pre/(R T prehyd)] [ d [ [D (Gw) / Dt]_adiab ] / d log(prehyd) ]
!        + [pre/(R T prehyd)] grad(Gw).[ d V / d log(prehyd) ]
!      "R" is the version of R (may be Rdry or Rmoist) used in the definition of vertical divergence "dver".

!      Calculation of the Laplacien term with VFE:
!      the Laplacian term:
!       d [ d (pre - prehyd) / d prehyd ]/ d eta
!      is rewritten as the sum of two terms:
!       (1/prehyd) * d[(pre - prehyd)/prehyd]/d eta * d[(d eta/d prehyd) prehyd**2]/d eta
!       + d eta/d(log(prehyd)) (d**2 [(pre - prehyd)/prehyd]/d eta**2)

!**   Interface.
!     ----------
!        *CALL* *GNH_TNDLAGADIAB_SVD(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!           KFLEV     : number of levels.
!           KPROMA    : horizontal dimension.
!           KST       : start of work.
!           KPROF     : working length.
!           PUF       : U-wind at full levels.
!           PVF       : V-wind at full levels.
!           PUH       : U-wind at half levels.
!           PVH       : V-wind at half levels.
!           PRDPHI    : term "pre / (R T prehyd delta)" at full levels.
!                       "R" is the version of R (may be Rdry or Rmoist) used in the definition of vertical divergence "dver".
!           PGWH      : Gw at half levels.
!           PGWF      : Gw at full levels.
!           PGWHL     : zonal comp. of grad(Gw) at half levels.
!           PGWHM     : merid comp. of grad(Gw) at half levels.
!           PGWFL     : zonal comp. of grad(Gw) at full levels.
!           PGWFM     : merid comp. of grad(Gw) at full levels.
!           PDVER     : vertical divergence "dver" at full levels.
!           PDIV      : horizontal divergence at full levels.
!           P3DIVG    : total divergence D3 at full levels.
!           PTNDGWH_LAP: tendency [D (Gw) / Dt]_adiab at half levels (Lapl term)
!           PTNDGWH_OTH: tendency [D (Gw) / Dt]_adiab at half levels (other)
!           PTNDGWF_OTH: tendency [D (Gw) / Dt]_adiab at full levels (other)
!           PDEP      : (pre-prehyd) at full levels.
!           PREH      : prehyd at half levels.
!           PREF      : prehyd at full levels.
!           PDELP     : [Delta prehyd] at full levels.
!           PLNPR     : quantity "delta" at full levels.

!         * OUTPUT:
!           PTNDVD    : tendency [D (svd) / Dt]_adiab at full levels.

!         Remark: input data PUF to PTNDGWH are at time t.

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
!   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!   K. Yessad (Feb 2018): remove deep-layer formulations.
!   J. Vivoda and P. Smolikova (Sep 2020): VFE pruning.
!   H Petithomme (Dec 2020): delete ZRDPHI
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCVER      , ONLY : LVFE_LAPL_BC,LVFE_GW
USE YOMCST       , ONLY : TCST
USE YOMDYNA      , ONLY : NVDVAR

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUF(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVF(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUH(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVH(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDPHI(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWH(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWF(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWHL(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWHM(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWFL(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWFM(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDVER(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIV(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P3DIVG(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTNDGWH_LAP(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTNDGWH_OTH(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTNDGWF_OTH(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDEP(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREH(KPROMA,0:KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREF(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELP(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLNPR(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTNDVD(KPROMA,KFLEV)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,JLEV
REAL(KIND=JPRB) :: ZRG2, ZDETA, Z2SAG
REAL(KIND=JPRB) :: ZLVGW(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZOTH(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZLAPL(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZDX(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZDDX(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZPI2MDETA(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZIN(KPROMA,0:KFLEV+1)
REAL(KIND=JPRB) :: ZDUF(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZDVF(KPROMA,KFLEV)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNH_TNDLAGADIAB_SVD',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       0.    PRELIMINARY CALCULATIONS
!              ------------------------

ZRG2=YDCST%RG*YDCST%RG
Z2SAG=2.0_JPRB/(YDCST%RG*YDCST%RA)

!     ------------------------------------------------------------------

!*       1.    LAPLACIAN CONTRIBUTION OF (D Gw/D t)
!              ------------------------------------

ZLAPL=0.0_JPRB

IF (LVFE_LAPL_BC) THEN

  DO JLEV=1,KFLEV
    DO JROF=KST,KPROF
      ZIN(JROF,JLEV)=PDEP(JROF,JLEV)/PREF(JROF,JLEV)
    ENDDO
  ENDDO
  ZIN(KST:KPROF,0)=0.0_JPRB
  ZIN(KST:KPROF,KFLEV+1)=ZIN(KST:KPROF,KFLEV)
  CALL VERDISINT(YDGEOMETRY%YRVFE,'FDER','00',KPROMA,KST,KPROF,KFLEV,ZIN,ZDX)
  CALL VERDISINT(YDGEOMETRY%YRVFE,'DDER','01',KPROMA,KST,KPROF,KFLEV,ZIN,ZDDX)

  ! d/deta (pi^2/[d pi/d eta])
  DO JLEV=1,KFLEV
    DO JROF=KST,KPROF
      ZIN(JROF,JLEV)=(PREF(JROF,JLEV)**2)/(PDELP(JROF,JLEV)*YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV))
    ENDDO
  ENDDO
  DO JROF=KST,KPROF
    ZIN(JROF,0)=PREH(JROF,0)*PREH(JROF,0)/(PDELP(JROF,1)*YDGEOMETRY%YRVETA%VFE_RDETAH(1))
    ZIN(JROF,KFLEV+1)=PREH(JROF,KFLEV)*PREH(JROF,KFLEV) &
     & /(PDELP(JROF,KFLEV)*YDGEOMETRY%YRVETA%VFE_RDETAH(KFLEV))
  ENDDO
  CALL VERDISINT(YDGEOMETRY%YRVFE,'FDER','00',KPROMA,KST,KPROF,KFLEV,ZIN,ZPI2MDETA)

  IF(LVFE_LAPL_BC)THEN

    DO JLEV=1,KFLEV
      ZDETA=1.0_JPRB/YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)
      DO JROF=KST,KPROF
        ZLAPL(JROF,JLEV) = -ZRG2*PRDPHI(JROF,JLEV)* ( &
         & ZPI2MDETA(JROF,JLEV)*ZDX(JROF,JLEV)*ZDETA/PREF(JROF,JLEV) &
         & + ZDDX(JROF,JLEV)*ZDETA*ZDETA/PLNPR(JROF,JLEV))
      ENDDO
    ENDDO

  ELSE

    DO JLEV=2,KFLEV-1
      ZDETA=1.0_JPRB/YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)
      DO JROF=KST,KPROF
        ZLAPL(JROF,JLEV) = -ZRG2*PRDPHI(JROF,JLEV)* ( &
         & ZPI2MDETA(JROF,JLEV)*ZDX(JROF,JLEV)*ZDETA/PREF(JROF,JLEV) &
         & + ZDDX(JROF,JLEV)*ZDETA*ZDETA/PLNPR(JROF,JLEV))
      ENDDO
    ENDDO

    DO JROF=KST,KPROF
      ZLAPL(JROF,1)=-PRDPHI(JROF,1)* &
       & (PTNDGWH_LAP(JROF,1)-PTNDGWH_LAP(JROF,0))
      ZLAPL(JROF,KFLEV)=-PRDPHI(JROF,KFLEV)* &
       & (PTNDGWH_LAP(JROF,KFLEV)-PTNDGWH_LAP(JROF,KFLEV-1))
    ENDDO

  ENDIF

ELSE

  DO JLEV=1,KFLEV
    DO JROF=KST,KPROF
      ZLAPL(JROF,JLEV)=-PRDPHI(JROF,JLEV)* &
       & (PTNDGWH_LAP(JROF,JLEV)-PTNDGWH_LAP(JROF,JLEV-1))
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

!*       2.    OTHER CONTRIBUTIONS OF (D Gw/D t)
!              ---------------------------------

! remark ky: zero for the time being.
! alternate formulations to be tested in the future may use PTNDGWH_OTH or PTNDGWF_OTH.
ZOTH(KST:KPROF,1:KFLEV)=0.0_JPRB

!     ------------------------------------------------------------------

!*       3.    Z-TERM
!              ------

IF (NVDVAR == 3 .OR. NVDVAR == 4) THEN

  ! Compute Z-term = (p/(RT)) grad(Gw).(dV/dprehyd)
  ! warning: formerly LVFE_Z_TERM
  IF (LVFE_GW) THEN
    ! vertical derivatives of wind
    ZIN(KST:KPROF,0)       = 0.0_JPRB
    ZIN(KST:KPROF,KFLEV+1) = 0.0_JPRB
    ZIN(KST:KPROF,1:KFLEV) = PUF(KST:KPROF,1:KFLEV)
    CALL VERDISINT(YDGEOMETRY%YRVFE,'FDER','11',KPROMA,KST,KPROF,KFLEV,ZIN,ZDUF)
    ZIN(KST:KPROF,1:KFLEV) = PVF(KST:KPROF,1:KFLEV)
    CALL VERDISINT(YDGEOMETRY%YRVFE,'FDER','11',KPROMA,KST,KPROF,KFLEV,ZIN,ZDVF)

    DO JLEV=1,KFLEV
      DO JROF=KST,KPROF
        ZLVGW(JROF,JLEV)=PRDPHI(JROF,JLEV)/YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)*&
         & (PGWFL(JROF,JLEV)*ZDUF(JROF,JLEV)+PGWFM(JROF,JLEV)*ZDVF(JROF,JLEV))
      ENDDO
    ENDDO
  ELSE
    DO JLEV=1,KFLEV
      DO JROF=KST,KPROF
        ZLVGW(JROF,JLEV)=((PUH(JROF,JLEV)-PUF(JROF,JLEV))*PGWHL(JROF,JLEV)&
         & +(PVH(JROF,JLEV)-PVF(JROF,JLEV))*PGWHM(JROF,JLEV)&
         & +(PUF(JROF,JLEV)-PUH(JROF,JLEV-1))*PGWHL(JROF,JLEV-1)&
         & +(PVF(JROF,JLEV)-PVH(JROF,JLEV-1))*PGWHM(JROF,JLEV-1)&
         & )*PRDPHI(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF

ENDIF

!     ------------------------------------------------------------------

!*       4.    COMPUTATION OF [D (svd) / Dt]_adiab: SUM OF ALL CONTRIBS
!              --------------------------------------------------------

IF (NVDVAR == 3 .OR. NVDVAR == 4) THEN
  DO JLEV=1,KFLEV
    DO JROF=KST,KPROF
      PTNDVD(JROF,JLEV)=PDVER(JROF,JLEV)*(PDIV(JROF,JLEV)-P3DIVG(JROF,JLEV)) &
        & +ZLAPL(JROF,JLEV)+ZOTH(JROF,JLEV)+ZLVGW(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNH_TNDLAGADIAB_SVD',1,ZHOOK_HANDLE)
END SUBROUTINE GNH_TNDLAGADIAB_SVD
