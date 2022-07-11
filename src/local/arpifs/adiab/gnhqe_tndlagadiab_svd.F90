!OCL  NOEVAL
SUBROUTINE GNHQE_TNDLAGADIAB_SVD(&
 ! --- INPUT -----------------------------------------------------------------
 & LDVFE_LAPL_BC, LDVERTFE, YDCST, YDGEOMETRY,KSTART,KPROF,&
 & POROGL,POROGM,POROGLM,POROGLL,POROGMM,&
 & PLNPR,PALPH,PREF,PREH,&
 & PUF,PVF,PUH,PVH,PTNDUS,PTNDVS,&
 & PEQCHAF,PKAPF,&
 & PRDPHI,PGWHL,PGWHM,PDVER,PNHXS,&
 ! --- OUTPUT ----------------------------------------------------------------
 & PTNDVD,PDBBC)

!**** *GNHQE_TNDLAGADIAB_SVD* - Computation of the adiabatic contribution of
!                               [D (svd) / Dt] in the NHQE model.
!                               "svd" is the vertical divergence variable.

!     Purpose.
!     --------

!     Computation of the adiabatic contribution of [D (svd) / Dt] in the NHQE model at full levels.
!     For svd=d4 (NVDVAR=4) the lagrangian tendency currently computed
!     is [D (dver) / Dt]_adiab (and not [D (dqua) / Dt]_adiab).

!     For svd=dver, expression of this tendency is:
!       [D (dver) / Dt]_adiab
!        = - dver (dver+NHXS) - [1/(R Tt)] [ d [ [D (Gw) / Dt]_adiab ] / d log(prehyd) ] + [1/(R Tt)] grad(Gw).[ d V / d log(prehyd) ]
!        = - dver ( d4 -NHXD) - [1/(R Tt)] [ d [ [D (Gw) / Dt]_adiab ] / d log(prehyd) ] + [1/(R Tt)] grad(Gw).[ d V / d log(prehyd) ]

!     Calculation of the Laplacien term with VFE: to be studied later, not yet coded.

!**   Interface.
!     ----------
!        *CALL* *GNHQE_TNDLAGADIAB_SVD(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!           YDGEOMETRY : structure containing geometry.
!           KSTART     : start of work.
!           KPROF      : working length.
!           POROGL     : zonal comp. of grad(surf orography).
!           POROGM     : merid comp. of grad(surf orography).
!           POROGLM    : -I
!           POROGLL    :  I- second order derivatives of "surf orography".
!           POROGMM    : -I
!           PLNPR      : "delta" at full levels.
!           PALPH      : "alpha" at full levels.
!           PREF       : "prehyd" at full levels.
!           PREH       : "prehyd" at half levels.
!           PUF        : U-wind at full levels.
!           PVF        : V-wind at full levels.
!           PUH        : U-wind at half levels.
!           PVH        : V-wind at half levels.
!           PTNDUS     : [DU/Dt]_surf.
!           PTNDVS     : [DV/Dt]_surf.
!           PEQCHAF    : exp(Kap Qcha) at full levels.
!           PKAPF      : "Kap=R/Cp" at full levels.
!           PRDPHI     : term "1 / (R Tt delta)" at full levels.
!           PGWHL      : zonal comp. of grad(Gw) at half levels.
!           PGWHM      : merid comp. of grad(Gw) at half levels.
!           PDVER      : vertical divergence "dver" at full levels.
!           PNHXS      : NHXS-term at full levels.

!         * OUTPUT:
!           PTNDVD     : tendency [D (svd) / Dt]_adiab at full levels.
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
!        K. YESSAD and F. VOITUS (after routine GNH_TNDLAGADIAB_SVD).
!        Original : March 2017

!     Modifications.
!     --------------
!     K. Yessad (March 2018): deep revision
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK

USE YOMCST       , ONLY : TCST

!     ------------------------------------------------------------------

IMPLICIT NONE

LOGICAL, INTENT (IN) :: LDVFE_LAPL_BC
LOGICAL, INTENT (IN) :: LDVERTFE
TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGLM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGLL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGMM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLNPR(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTNDUS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTNDVS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEQCHAF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDPHI(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWHL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGWHM(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDVER(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNHXS(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTNDVD(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDBBC(YDGEOMETRY%YRDIM%NPROMA)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,JLEV
REAL(KIND=JPRB) :: ZRG2,ZJACOB,ZACUGF,ZACVGF

REAL(KIND=JPRB) :: ZLVGW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZOTH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZLAPL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZIN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZOUT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZINS(YDGEOMETRY%YRDIM%NPROMA)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "gnhqe_lkap.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHQE_TNDLAGADIAB_SVD',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA,NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

!*       1.    LAPLACIAN CONTRIBUTION OF (D Gw/D t)
!              ------------------------------------

ZRG2=YDCST%RG*YDCST%RG

! warning: formerly LVFE_LAPL
IF (LDVERTFE) CALL ABOR1(' GNHQE_TNDLAGADIAB_SVD 1.1: option not yet coded')

  ! ky: in these calculations, for the time being, a term containing [D (1/Kap)/ D log(prehyd)] has been neglected.

! * Compute [(1/g) D w_surf/Dt]
DO JROF=KSTART,KPROF
  ZJACOB=POROGLL(JROF)*PUH(JROF,NFLEVG)*PUH(JROF,NFLEVG) &
   & + POROGMM(JROF)*PVH(JROF,NFLEVG)*PVH(JROF,NFLEVG) &
   & + 2.0_JPRB*POROGLM(JROF)*PUH(JROF,NFLEVG)*PVH(JROF,NFLEVG)
  ZACUGF=PTNDUS(JROF)*POROGL(JROF)
  ZACVGF=PTNDVS(JROF)*POROGM(JROF)
  ZINS(JROF)=(ZACUGF+ZACVGF+ZJACOB)/ZRG2
  PDBBC(JROF)=-ZACUGF-ZACVGF-ZJACOB
ENDDO

! * Compute [exp(Kap Qcha) - 1] at full levels:
DO JLEV=1,NFLEVG
  DO JROF=KSTART,KPROF
    ZIN(JROF,JLEV)=PEQCHAF(JROF,JLEV)-1.0_JPRB
  ENDDO
ENDDO

! * Apply Laplacian operator:
!   ZOUT contains ( [Kap/G**2] * [ d [ [D (Gw) / Dt]_adiab ] / d log(prehyd) ] )
CALL GNHQE_LKAP(LDVFE_LAPL_BC,YDGEOMETRY,KSTART,KPROF,ZIN,ZINS,PLNPR,PALPH,PREF,PREH,PKAPF,ZOUT)

! * Compute - [1/(R Tt)] [ d [ [D (Gw) / Dt]_adiab ] / d log(prehyd) ]
!   i.e. - ([G**2/Kap] [PRDPHI] [delta]) * ZOUT
DO JLEV=1,NFLEVG
  DO JROF=KSTART,KPROF
    ZLAPL(JROF,JLEV)=-(ZRG2*PRDPHI(JROF,JLEV)*PLNPR(JROF,JLEV)/PKAPF(JROF,JLEV))*ZOUT(JROF,JLEV)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*       2.    OTHER CONTRIBUTIONS OF (D Gw/D t)
!              ---------------------------------

! ky: the currently neglected term containing [D (1/Kap)/ D log(prehyd)] could go there, or could be added to ZLAPL.

ZOTH(KSTART:KPROF,1:NFLEVG)=0.0_JPRB

!     ------------------------------------------------------------------

!*       3.    Z-TERM
!              ------

! Compute (1/(R Tt)) grad(Gw).(dV/delta)
! warning: formerly LVFE_Z_TERM (dead case LVERTFE)
IF (LDVERTFE) CALL ABOR1(' GNHQE_TNDLAGADIAB_SVD 3.1: option not yet coded')

DO JLEV=1,NFLEVG
  DO JROF=KSTART,KPROF
    ZLVGW(JROF,JLEV)=((PUH(JROF,JLEV)-PUF(JROF,JLEV))*PGWHL(JROF,JLEV)&
     & +(PVH(JROF,JLEV)-PVF(JROF,JLEV))*PGWHM(JROF,JLEV)&
     & +(PUF(JROF,JLEV)-PUH(JROF,JLEV-1))*PGWHL(JROF,JLEV-1)&
     & +(PVF(JROF,JLEV)-PVH(JROF,JLEV-1))*PGWHM(JROF,JLEV-1)&
     & )*PRDPHI(JROF,JLEV)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*       4.    COMPUTATION OF [D (svd) / Dt]_adiab: SUM OF ALL CONTRIBS
!              --------------------------------------------------------

DO JLEV=1,NFLEVG
  DO JROF=KSTART,KPROF
    PTNDVD(JROF,JLEV)=-PDVER(JROF,JLEV)*(PDVER(JROF,JLEV)+PNHXS(JROF,JLEV)) &
      & +ZLAPL(JROF,JLEV)+ZOTH(JROF,JLEV)+ZLVGW(JROF,JLEV)
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHQE_TNDLAGADIAB_SVD',1,ZHOOK_HANDLE)
END SUBROUTINE GNHQE_TNDLAGADIAB_SVD
