!OCL  NOEVAL
SUBROUTINE GNHEE_TNDLAGADIAB_SVD(&
 ! --- INPUT -----------------------------------------------------------------
 & LDVFE_LAPL_BC, LDVERTFE, KVDVAR, YDCST, YDGEOMETRY,KST,KPROF,&
 & POROGL,POROGM,POROGLM,POROGLL,POROGMM,&
 & PEVEL,PRDELP,PLNPR,PALPH,PREF,PREH,PDEP,&
 & PUF,PVF,PUH,PVH,PTNDUF,PTNDVF,PTNDUS,PTNDVS,&
 & PRDPHI,PGWHL,PGWHM,PGWFL,PGWFM,&
 & PDVER,PDVERW,PDIV,P3DIVG,&
 ! --- OUTPUT ----------------------------------------------------------------
 & PTNDVD,PDBBC&
 & )  

!**** *GNHEE_TNDLAGADIAB_SVD* - Computation of the adiabatic contribution of
!                               [D (svd) / Dt] in the NHEE model.
!                               "svd" is the vertical divergence variable.

!     This alternate version of GNH_TNDLAGADIAB_SVD directly computes the Laplacian term
!     and does not require input value of upper air [ d [ [D (Gw) / Dt]_adiab ].
!     Structure matches GNHQE_TNDLAGADIAB_SVD.
!     Called only if LGWADV=F

!     VFE code not yet provided for the time being.

!     Purpose.
!     --------

!     Computation of the adiabatic contribution of [D (svd) / Dt] in the NHEE model
!     at full levels.
!     For svd=d4 (NVDVAR=4) the lagrangian tendency currently computed
!     is [D (dver) / Dt]_adiab (and not [D (dqua) / Dt]_adiab).
!     For svd=d5 (NVDVAR=5) the lagrangian tendency currently computed
!     is [D (d13) / Dt]_adiab (and not [D (dqua) / Dt]_adiab).

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

!      Calculation of the Laplacien term with VFE: not yet coded in GNHEE_TNDLAGADIAB_SVD.

!**   Interface.
!     ----------
!        *CALL* *GNHEE_TNDLAGADIAB_SVD(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!           YDGEOMETRY : structure containing geometry.
!           KST        : start of work.
!           KPROF      : working length.
!           POROGL     : zonal comp. of grad(surf orography).
!           POROGM     : merid comp. of grad(surf orography).
!           POROGLM    : -I
!           POROGLL    :  I- second order derivatives of "surf orography".
!           POROGMM    : -I
!           PEVEL      : [etadot dprehyd/deta] at half levels.
!           PRDELP     : [1/Delta prehyd] at full levels.
!           PLNPR      : "delta" at full levels.
!           PALPH      : "alpha" at full levels.
!           PREF       : "prehyd" at full levels.
!           PREH       : "prehyd" at half levels.
!           PDEP       : (pre-prehyd) at full levels.
!           PUF        : U-wind at full levels.
!           PVF        : V-wind at full levels.
!           PUH        : U-wind at half levels.
!           PVH        : V-wind at half levels.
!           PTNDUF     : [DU/Dt] at full levels (used for NVDVAR=5).
!           PTNDVF     : [DV/Dt] at full levels (used for NVDVAR=5).
!           PTNDUS     : [DU/Dt]_surf.
!           PTNDVS     : [DV/Dt]_surf.
!           PRDPHI     : term "pre / (R T prehyd delta)" at full levels.
!                        "R" is the version of R (may be Rdry or Rmoist) used in the definition of vertical divergence "dver".
!           PGWHL      : zonal comp. of grad(Gw) at half levels.
!           PGWHM      : merid comp. of grad(Gw) at half levels.
!           PGWFL      : zonal comp. of grad(Gw) at full levels.
!           PGWFM      : merid comp. of grad(Gw) at full levels.
!           PDVER      : vertical divergence "dver" at full levels.
!           PDVERW     : modified vertical divergence "d13" at full levels (for NVDVAR=5).
!           PDIV       : horizontal divergence at full levels.
!           P3DIVG     : total divergence D3 at full levels.

!         * OUTPUT:
!           PTNDVD     : tendency [D (svd) / Dt]_adiab at full levels.
!           PDBBC      : [D (Gw)_surf / Dt]_adiab.


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
!        See documentations Arpege about the NHEE model.

!     Author.
!     -------
!        K. YESSAD and F. VOITUS (after routine GNH_TNDLAGADIAB_SVD and GNHQE_TNDLAGADIAB_SVD).
!        Original : July 2018

!     Modifications.
!     --------------
!       K. Yessad (Aug 2018): NVDVAR=13 and 14.
!       J. Vivoda and P. Smolikova (Sep 2020): VFE pruning.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK

USE YOMCST       , ONLY : TCST


!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL, INTENT (IN) :: LDVFE_LAPL_BC
LOGICAL, INTENT (IN) :: LDVERTFE
INTEGER (KIND=JPIM), INTENT (IN) :: KVDVAR
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGLM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGLL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGMM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVEL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLNPR(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDEP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTNDUF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTNDVF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTNDUS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTNDVS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDPHI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWHL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWHM(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWFM(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDVER(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDVERW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P3DIVG(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTNDVD(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDBBC(YDGEOMETRY%YRDIM%NPROMA)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,JLEV
REAL(KIND=JPRB) :: ZRG2,ZJACOB,ZACUGF,ZACVGF
REAL(KIND=JPRB) :: ZLVGW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZLAPL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZIN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZOUT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZINS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB) :: ZIN_Z(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZDUF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDVF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPART_RESCALED_PDEPF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZATND13(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

! other contributions, currently unused
!REAL(KIND=JPRB) :: ZOTH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "gnhee_lapl.intfb.h"
#include "gnhee_pdep_orog.intfb.h"
#include "gnhee_tndlagadiab_forcorog.intfb.h"
#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHEE_TNDLAGADIAB_SVD',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA,NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

!*       1.    LAPLACIAN CONTRIBUTION OF (D Gw/D t)
!              ------------------------------------

ZRG2=YDCST%RG*YDCST%RG

IF (KVDVAR==3 .OR. KVDVAR==4 .OR. KVDVAR==5) THEN
  ! coded in GNH_TNDLAGADIAB_SVD but not easy to transpose there; to be studied later.
  IF (LDVERTFE) CALL ABOR1(' GNHEE_TNDLAGADIAB_SVD 1.1: option not yet coded')

  ! * Compute [exp(Qcha) - 1]=(pre-prehyd)/prehyd at full levels:
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      ZIN(JROF,JLEV)=PDEP(JROF,JLEV)/PREF(JROF,JLEV)
    ENDDO
  ENDDO

  IF (KVDVAR == 3 .OR. KVDVAR == 4) THEN
    ! * Compute [(1/g) D w_surf/Dt]
    DO JROF=KST,KPROF
      ZJACOB=POROGLL(JROF)*PUH(JROF,NFLEVG)*PUH(JROF,NFLEVG) &
       & + POROGMM(JROF)*PVH(JROF,NFLEVG)*PVH(JROF,NFLEVG) &
       & + 2.0_JPRB*POROGLM(JROF)*PUH(JROF,NFLEVG)*PVH(JROF,NFLEVG)
      ZACUGF=PTNDUS(JROF)*POROGL(JROF)
      ZACVGF=PTNDVS(JROF)*POROGM(JROF)
      ZINS(JROF)=(ZACUGF+ZACVGF+ZJACOB)/ZRG2
      PDBBC(JROF)=-ZACUGF-ZACVGF-ZJACOB
    ENDDO
  ELSEIF (KVDVAR==5) THEN
    ! * Bottom boundary condition is 0
    DO JROF=KST,KPROF
      ZINS(JROF)=0.0_JPRB
      PDBBC(JROF)=0.0_JPRB
    ENDDO

    ! * Compute half level (1/g**2) * D(S V grad(Phi_s))/Dt:
    CALL GNHEE_TNDLAGADIAB_FORCOROG(YDGEOMETRY,KST,KPROF,&
     & POROGL,POROGM,POROGLM,POROGLL,POROGMM,&
     & PUH,PVH,PTNDUF,PTNDVF,PEVEL,PRDELP,PLNPR,PALPH,&
     & ZATND13)
  
    ! * Compute the equivalent '(p-Pi)/Pi' generated by term 'S V grad(Phi_s)'.
    CALL GNHEE_PDEP_OROG(YDGEOMETRY,KST,KPROF,ZATND13,PREF,ZPART_RESCALED_PDEPF)

    ! * Correct ZIN:
    DO JLEV=1,NFLEVG 
      DO JROF=KST,KPROF
        ZIN(JROF,JLEV)=ZIN(JROF,JLEV)-ZPART_RESCALED_PDEPF(JROF,JLEV)
      ENDDO
    ENDDO
  ENDIF

  ! * Apply Laplacian operator:
  !   ZOUT contains ( [1/G**2] * [ d [ [D (Gw) / Dt]_adiab ] / d log(prehyd) ] )
  CALL GNHEE_LAPL(LDVFE_LAPL_BC,YDGEOMETRY,KST,KPROF,ZIN,ZINS,PLNPR,PALPH,PREF,PREH,ZOUT)

  ! * Compute - [pre/(prehyd RT)] [ d [ [D (Gw) / Dt]_adiab ] / d log(prehyd) ]
  !   i.e. - ([G**2] [PRDPHI] [delta]) * ZOUT
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      ZLAPL(JROF,JLEV)=-(ZRG2*PRDPHI(JROF,JLEV)*PLNPR(JROF,JLEV))*ZOUT(JROF,JLEV)
    ENDDO
  ENDDO

!     ------------------------------------------------------------------

!*       2.    OTHER CONTRIBUTIONS OF (D Gw/D t)
!              ---------------------------------

! remark ky: zero for the time being.
!ZOTH(KST:KPROF,1:NFLEVG)=0.0_JPRB

!     ------------------------------------------------------------------

!*       3.    Z-TERM
!              ------

  ! Compute (p/RT) grad(Gw).(dV/dprehyd), or (p/RT) grad(GW).(dV/dprehyd) if NVDVAR=5
  ! warning: formerly LVFE_Z_TERM (dead case LVERTFE)
  IF (LDVERTFE) THEN
    ! vertical derivatives of wind
    ZIN_Z(KST:KPROF,0)       = 0.0_JPRB
    ZIN_Z(KST:KPROF,NFLEVG+1) = 0.0_JPRB
    ZIN_Z(KST:KPROF,1:NFLEVG) = PUF(KST:KPROF,1:NFLEVG)
    CALL VERDISINT(YDGEOMETRY%YRVFE,'FDER','11',NPROMA,KST,KPROF,NFLEVG,ZIN_Z,ZDUF)
    ZIN_Z(KST:KPROF,1:NFLEVG) = PVF(KST:KPROF,1:NFLEVG)
    CALL VERDISINT(YDGEOMETRY%YRVFE,'FDER','11',NPROMA,KST,KPROF,NFLEVG,ZIN_Z,ZDVF)

    DO JLEV=1,NFLEVG
      DO JROF=KST,KPROF
        ZLVGW(JROF,JLEV)=PRDPHI(JROF,JLEV)/YDGEOMETRY%YRVETA%VFE_RDETAH(JLEV)*&
         & (PGWFL(JROF,JLEV)*ZDUF(JROF,JLEV)+PGWFM(JROF,JLEV)*ZDVF(JROF,JLEV))
      ENDDO
    ENDDO
  ELSE
    DO JLEV=1,NFLEVG
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

IF (KVDVAR == 3 .OR. KVDVAR == 4) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      PTNDVD(JROF,JLEV)=PDVER(JROF,JLEV)*(PDIV(JROF,JLEV)-P3DIVG(JROF,JLEV)) &
        & +ZLAPL(JROF,JLEV)+ZLVGW(JROF,JLEV) !+ZOTH(JROF,JLEV)
    ENDDO
  ENDDO
ELSEIF (KVDVAR==5) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KPROF
      PTNDVD(JROF,JLEV)=PDVERW(JROF,JLEV)*(PDIV(JROF,JLEV)-P3DIVG(JROF,JLEV)) &
        & +ZLAPL(JROF,JLEV)+ZLVGW(JROF,JLEV) !+ZOTH(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHEE_TNDLAGADIAB_SVD',1,ZHOOK_HANDLE)
END SUBROUTINE GNHEE_TNDLAGADIAB_SVD
