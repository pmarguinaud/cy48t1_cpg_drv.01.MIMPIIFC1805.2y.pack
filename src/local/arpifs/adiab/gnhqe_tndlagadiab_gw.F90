!OCL  NOEVAL
SUBROUTINE GNHQE_TNDLAGADIAB_GW(&
 ! --- INPUT -----------------------------------------------------------------
 & YDCST, YDGEOMETRY,KSTART,KPROF,&
 & POROGL,POROGM,POROGLM,POROGLL,POROGMM,&
 & PUS,PVS,PTNDUS,PTNDVS,&
 & PLNPR,PALPH,PREF,PREH,&
 & PEQCHAF,PKAPF,&
 ! --- OUTPUT ----------------------------------------------------------------
 & PTNDGW)  

!**** *GNHQE_TNDLAGADIAB_GW* - Computation of the adiabatic contribution of
!                              [D (Gw) / Dt] in the NHQE model.

!     Purpose.
!     --------

!     Computation of the adiabatic contribution of [D (Gw) / Dt] in the NHQE model.
!     Upper air and surface [D (Gw) / Dt] are computed.

!     Upper air [D (Gw) / Dt] for thin layer equations:
!      [D (Gw) / Dt]_adiab = G**2/Kap ( d (prehyd (EQ-1))/ d prehyd) + (Kap-1)*(EQ-1) )
!     where Kap=R/cp and EQ=exp(Kap*Qcha).

!     Surface [D (Gw) / Dt] for thin layer equations:
!      a temporal derivation of the lower boundary condition
!      [Gw]_surf = V_surf grad[Phi_surf]
!      is done; that yields:
!      [D [Gw]_surf / Dt]_adiab = [D V_surf / Dt] grad[Phi_surf]
!                                 + V_surf . [V_surf . grad[grad[Phi_surf]]]

!     In the NHQE model the first-order derivative to compute [D (Gw) / Dt] from (EQ-1)
!     is not convenient to discretise, it is computed as a product between a "G-type" vertical integral
!     and a "L_kap"-type double derivative.

!**   Interface.
!     ----------
!        *CALL* *GNHQE_TNDLAGADIAB_GW(...)

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
!           PUS        : surface U wind.
!           PVS        : surface V wind.
!           PTNDUS     : [DU/Dt]_surf.
!           PTNDVS     : [DV/Dt]_surf.
!           PLNPR      : "delta" at full levels.
!           PALPH      : "alpha" at full levels.
!           PREF       : "prehyd" at full levels.
!           PREH       : "prehyd" at half levels.
!           PEQCHAF    : exp(Kap Qcha) at full levels.
!           PKAPF      : "Kap=R/Cp" at full levels.

!         * OUTPUT:
!           PTNDGW     : LVFE_GW=F: [D (Gw) / Dt]_adiab at half levels 0 to NFLEVG-1.
!                        LVFE_GW=T: [D (Gw) / Dt]_adiab at full levels 1 to NFLEVG.

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
!        K. YESSAD and F. VOITUS (after routine GNH_TNDLAGADIAB_GW).
!        Original : March 2017

!     Modifications.
!     --------------
!     K. Yessad (March 2018): deep revision
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMCST       , ONLY : TCST
USE YOMCVER      , ONLY : LVFE_LAPL_BC, LVFE_GW

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCST)        ,INTENT(IN)    :: YDCST
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGLM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGLL(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROGMM(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTNDUS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTNDVS(YDGEOMETRY%YRDIM%NPROMA)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLNPR(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREH(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEQCHAF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKAPF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTNDGW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF,JLEV
REAL(KIND=JPRB) :: ZRG2,ZJACOB,ZACUGF,ZACVGF
REAL(KIND=JPRB) :: ZIN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZOUT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZINS(YDGEOMETRY%YRDIM%NPROMA)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

! -----------------------------------------------------------------------------

#include "abor1.intfb.h"
#include "gnhqe_lkap.intfb.h"    

! -----------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GNHQE_TNDLAGADIAB_GW',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA,NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

ZRG2=YDCST%RG*YDCST%RG

!     ------------------------------------------------------------------

!*       1.    COMPUTATION OF [D [Gw]_surf / Dt]_adiab.
!              ----------------------------------------

! * Compute [(1/g) D w_surf/Dt]
DO JROF=KSTART,KPROF
  ZJACOB=POROGLL(JROF)*PUS(JROF)*PUS(JROF) &
   & + POROGMM(JROF)*PVS(JROF)*PVS(JROF) &
   & + 2.0_JPRB*POROGLM(JROF)*PUS(JROF)*PVS(JROF)
  ZACUGF=PTNDUS(JROF)*POROGL(JROF)
  ZACVGF=PTNDVS(JROF)*POROGM(JROF)
  ZINS(JROF)=(ZACUGF+ZACVGF+ZJACOB)/ZRG2
ENDDO

!     ------------------------------------------------------------------

!*       2.   COMPUTATION OF UPPER AIR [D (Gw) / Dt]_adiab.
!             ---------------------------------------------

! ky: a term containing [d (1/Kap)/ d log(prehyd)] has provisionally be neglected there.

! * Compute [exp(Kap Qcha) - 1] at full levels:
DO JLEV=1,NFLEVG
  DO JROF=KSTART,KPROF
    ZIN(JROF,JLEV)=PEQCHAF(JROF,JLEV)-1.0_JPRB
  ENDDO
ENDDO

! * Apply Laplacian (L_kap) operator:  
!   ZOUT contains ( [Kap/G**2] * [ d [ [D (Gw) / Dt]_adiab ] / d log(prehyd) ] )
! warning: formerly LVFE_LAPL
IF (LVFE_LAPL_BC) THEN
  ! Only VFD calculation of Laplacian is currently available.
  CALL ABOR1(' GNHQE_TNDLAGADIAB_GW 2a: option not yet coded!')
ELSE
  CALL GNHQE_LKAP(YDGEOMETRY,KSTART,KPROF,ZIN,ZINS,PLNPR,PALPH,PREF,PREH,PKAPF,ZOUT)
ENDIF

! * Multiply ZOUT by [G**2/Kap]:
DO JLEV=1,NFLEVG
  DO JROF=KSTART,KPROF
    ZOUT(JROF,JLEV)=ZRG2*ZOUT(JROF,JLEV)/PKAPF(JROF,JLEV)
  ENDDO
ENDDO

! * Integrate, to compute upper-air [D (Gw) / Dt]_adiab:
!   This is a "G-type" vertical integral (cf. GPGEO, GPGW).
IF (LVFE_GW) THEN
  ! In this case, [D (Gw) / Dt]_adiab is at full levels (in PTNDGW).
  ! VFE integration to compute PTNDGW is not yet coded.
  CALL ABOR1(' GNHQE_TNDLAGADIAB_GW 2b: option not yet coded!')
ELSE
  ! In this case, [D (Gw) / Dt]_adiab is at half levels (in PTNDGW).
  ! Caution: top value is in PTNDGW(.,1), lbar value is in PTNDGW(.,jlev=lbar+1).
  ! It is not necessary to save [D (Gw) / Dt]_surf

  ! Half level lbar=L-1:
  DO JROF=KSTART,KPROF
    PTNDGW(JROF,NFLEVG)=ZRG2*ZINS(JROF)-PLNPR(JROF,NFLEVG)*ZOUT(JROF,NFLEVG)
  ENDDO

  ! Half level lbar=L-2 to 0:
  DO JLEV=NFLEVG-1,1,-1
    DO JROF=KSTART,KPROF
      PTNDGW(JROF,JLEV)=PTNDGW(JROF,JLEV+1)-PLNPR(JROF,JLEV)*ZOUT(JROF,JLEV)
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GNHQE_TNDLAGADIAB_GW',1,ZHOOK_HANDLE)
END SUBROUTINE GNHQE_TNDLAGADIAB_GW
