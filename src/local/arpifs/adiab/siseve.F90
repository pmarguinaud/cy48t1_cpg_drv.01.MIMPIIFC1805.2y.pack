SUBROUTINE SISEVE(YDDYNA, YDGEOMETRY,YDDYN,KLEV,KLON,PV1,PV2,KNLON,LD_LSTAR)

!**** *SISEVE* - SEcond VErtical derivative operator in Semi - Implicit

!     Purpose.
!     --------
!           Operator L_star or L_star_star to compute vertical Laplacian in non-hydrostatic
!           dynamics: formulation for semi - implicit

!           The discretisation of the vertical Laplacian must remain
!           consistent with the discretisation of the equivalent quantity
!           in the non linear part (pieces of codes spread in GNHPREH,
!           GNH_TNDLAGADIAB_GW and GNH_TNDLAGADIAB_SVD).

!           This vertical Laplacian is applied to a quantity "X" (in PV1)
!           which matches the following conditions:
!           X_top = 0
!           X_surf = X(l=nflevg)
!           (see for example the calculation of the top and bottom
!           pressure departure variable in GNHPREH).

!**   Interface.
!     ----------
!        *CALL* *SISEVE(...)

!        Explicit arguments :
!        --------------------
!        KLEV   : DISTANCE IN MEMORY BETWEEN VALUES OF THE VARIABLES      (IN)
!                 PV1 OR PV2 AT THE VERTICAL
!        KLON   : DISTANCE IN MEMORY BETWEEN VALUES OF THE VARIABLES      (IN)
!                 PV1 OR PV2 AT THE SAME LEVEL
!        KNLON  : NUMBER OF VERTICAL COLUMNS TREATED                      (IN)

!           TYPICAL VALUES ARE  (KLEV,KLON)=(NDLSUR,1)  FOR GRID POINT ARRAY
!                               (KLEV,KLON)=(1,NFLSUR)  FOR SPECTRAL ARRAY

!        PV1   : INPUT VARIABLE                                           (IN)
!        PV2   : DERIVED VARIABLE BY VERTICAL LAPLACIAN                   (OUT)
!        LD_LSTAR : T: operator L_star                                    (OPT IN)
!                   F: operator L_star_star

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Arpege/Aladin documentation

!     Author.
!     -------
!      Radmila Bubnova GMAP/COMPAS, stage MRE
!      Original : 93-03-21

!     Modifications.
!     --------------
!      P. Smolikova and J. Vivoda (Oct 2013): new options for VFE-NH
!      K. Yessad (Dec 2016): Prune obsolete options.
!      K. Yessad (June 2017): Vertical-dependent SITRA.
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!      F. Voitus (Feb 2020): new refinement option for VFD-NH
!      J. Vivoda and P. Smolikova (Sep 2020): VFE pruning.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK
USE YOMDYN       , ONLY : TDYN
USE YOMCVER      , ONLY : LVERTFE, LVFE_LAPL_BC
USE YOMDYNA      , ONLY : TDYNA

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TDYNA), INTENT (IN) :: YDDYNA
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TDYN)        ,INTENT(IN)    :: YDDYN
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNLON 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV1(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PV2(1+(YDGEOMETRY%YRDIMV%NFLEVG-1)*KLEV+(KNLON-1)*KLON) 
LOGICAL,OPTIONAL  ,INTENT(IN)    :: LD_LSTAR

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZA(YDGEOMETRY%YRDIMV%NFLEVG),  ZC(YDGEOMETRY%YRDIMV%NFLEVG), ZRM(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZAS(YDGEOMETRY%YRDIMV%NFLEVG), ZBS(YDGEOMETRY%YRDIMV%NFLEVG),ZCS(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZDS(YDGEOMETRY%YRDIMV%NFLEVG), ZES(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPRTOP,ZDPRE,ZDETA, ZQBBC
REAL(KIND=JPRB) :: Z0,Z1,Z2,Z3,Z4,Z5

INTEGER(KIND=JPIM) :: ILEV0, ILEV1, ILEV2, ILEV3, ILEV4
INTEGER(KIND=JPIM) :: IBOT0, IBOT1, IBOT2, IBOT3
INTEGER(KIND=JPIM) :: ITOP0, ITOP1, ITOP2, ITOP3
INTEGER(KIND=JPIM) :: JLEV, JLON, IDT

REAL(KIND=JPRB) :: ZIN1(KNLON,0:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB) :: ZINDETA(KNLON,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZINDETA2(KNLON,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZP2MDETA(KNLON,YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZRATIO_SITR(YDGEOMETRY%YRDIMV%NFLEVG)
LOGICAL         :: LL_LSTAR

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "verdisint.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SISEVE',0,ZHOOK_HANDLE)

ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE)

ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,   SILNPR=>YDDYN%SILNPR,   SIPR=>YDDYN%SIPR, SIRDEL=>YDDYN%SIRDEL, &
& SITLAF=>YDDYN%SITLAF,   SITR=>YDDYN%SITR,SITRAM=>YDDYN%SITRAM,     SIWEIG=>YDDYN%SIWEIG,SIRSLP=>YDDYN%SIRSLP&
& )
!     ------------------------------------------------------------------

!*       1.    BOUNDARY CONDITIONS FOR FD
!              --------------------------

IF( .NOT.(LVERTFE.AND.LVFE_LAPL_BC) )THEN

  ! --- Finite differences ---
  IF (YDDYNA%LNHEE_REFINE_SILAPL) THEN
  ! * Compute the coefficients "A" and "C" at full level 1:
  !   "A" must be zero if the top pressure of the model is zero.
     ZPRTOP = YDVAB%VAH(0)+YDVAB%VBH(0)*SIPR
     ZA(1)  = ZPRTOP/(SILNPR(1)*(SITLAF(1)-ZPRTOP))
     ZC(1)  = SITLAF(2)/(SILNPR(1)*(SITLAF(2)-SITLAF(1)))

     Z0     =  SIRSLP*YDVAB%VRATH(0)*(YDVAB%VRATH(0)*SIWEIG(1,1)-YDVAB%VRATH(1)*SIWEIG(1,1)*SIWEIG(1,3))
     Z1     =  SIRSLP*YDVAB%VRATH(1)*(YDVAB%VRATH(1)*(SIWEIG(2,1)*(1.0_JPRB-SIWEIG(1,3)) &
            & +SIWEIG(1,3)*SIWEIG(1,2))-YDVAB%VRATH(0)*SIWEIG(1,2))  
     Z2     =  SIRSLP*YDVAB%VRATH(2)*YDVAB%VRATH(1)*SIWEIG(2,2)*(1.0_JPRB-SIWEIG(1,3))

     ZAS(1) =  0.0_JPRB
     ZBS(1) = -Z1*SITLAF(1)/(SILNPR(1)*(SITLAF(2)-SITLAF(1))) &
            & -Z0/SILNPR(1)
     ZCS(1) =  Z1*SITLAF(2)/(SILNPR(1)*(SITLAF(2)-SITLAF(1))) &
            & -Z2*SITLAF(2)/(SILNPR(1)*(SITLAF(3)-SITLAF(2)))
     ZDS(1) =  Z2*SITLAF(3)/(SILNPR(1)*(SITLAF(3)-SITLAF(2)))
     ZES(1) =  0.0_JPRB      

     ! * Top:
     DO JLON = 1, KNLON
       ITOP1 = 1 + (JLON-1)*KLON
       ITOP2 = ITOP1 + KLEV
       ITOP3 = ITOP2 + KLEV
       PV2(ITOP1) = -ZA(1)*PV1(ITOP1)+ZC(1)*(PV1(ITOP2)-PV1(ITOP1)) &
                  & + ZBS(1)*PV1(ITOP1) + ZCS(1)*PV1(ITOP2) + ZDS(1)*PV1(ITOP3)
     ENDDO

  ! * Compute the coefficients "A" and "C" at full level 2:
  !   "A" must be zero if the top pressure of the model is zero.
     
     ZA(2)  = SITLAF(1)/(SILNPR(2)*(SITLAF(2)-SITLAF(1)))
     ZC(2)  = SITLAF(3)/(SILNPR(2)*(SITLAF(3)-SITLAF(2)))

     Z0     =  SIRSLP*YDVAB%VRATH(0)*YDVAB%VRATH(1)*SIWEIG(1,1)*SIWEIG(1,3)
     Z1     =  SIRSLP*YDVAB%VRATH(1)*(YDVAB%VRATH(1)*(SIWEIG(2,1)*(1.0_JPRB-SIWEIG(1,3))+SIWEIG(1,3)*SIWEIG(1,2)) &
            & -YDVAB%VRATH(2)*SIWEIG(2,1)*SIWEIG(2,3))
     Z2     =  SIRSLP*YDVAB%VRATH(2)*(YDVAB%VRATH(2)*(SIWEIG(3,1)*(1.0_JPRB-SIWEIG(2,3))+SIWEIG(2,3)*SIWEIG(2,2)) &
            & -YDVAB%VRATH(1)*SIWEIG(2,2)*(1.0_JPRB-SIWEIG(1,3)))
     Z3     =  SIRSLP*YDVAB%VRATH(3)*YDVAB%VRATH(2)*SIWEIG(3,2)*(1.0_JPRB-SIWEIG(2,3))     

     ZAS(2) =  Z1*SITLAF(1)/(SILNPR(2)*(SITLAF(2)-SITLAF(1)))-(Z0/SILNPR(2))
     ZBS(2) = -Z2*SITLAF(2)/(SILNPR(2)*(SITLAF(3)-SITLAF(2))) &
            & -Z1*SITLAF(2)/(SILNPR(2)*(SITLAF(2)-SITLAF(1)))
     ZCS(2) =  Z2*SITLAF(3)/(SILNPR(2)*(SITLAF(3)-SITLAF(2))) &
            & -Z3*SITLAF(3)/(SILNPR(2)*(SITLAF(4)-SITLAF(3)))
     ZDS(2) =  Z3*SITLAF(4)/(SILNPR(2)*(SITLAF(4)-SITLAF(3))) 
     ZES(2) =  0.0_JPRB

     DO JLON = 1, KNLON
       ITOP0 = 1 + (JLON-1)*KLON
       ITOP1 = ITOP0 + KLEV
       ITOP2 = ITOP1 + KLEV
       ITOP3 = ITOP2 + KLEV
       PV2(ITOP1) = ZA(2)*(PV1(ITOP0)-PV1(ITOP1))             &
                  & + ZC(2)*(PV1(ITOP2)-PV1(ITOP1))           &
                  & + ZAS(2)*PV1(ITOP0) + ZBS(2)*PV1(ITOP1)   &
                  & + ZCS(2)*PV1(ITOP2) + ZDS(2)*PV1(ITOP3)
     ENDDO


     ZQBBC = 0.0_JPRB
     IF (YDDYNA%LNHEE_REFINE_PREH_BBC) ZQBBC=1.0_JPRB

   ! * Compute the coefficients "A" and "C" at full level nflevg-1:
     ZA(NFLEVG-1) = SITLAF(NFLEVG-2)/(SILNPR(NFLEVG-1)*(SITLAF(NFLEVG-1)-SITLAF(NFLEVG-2)))
     ZC(NFLEVG-1) = SITLAF(NFLEVG)/(SILNPR(NFLEVG-1)*(SITLAF(NFLEVG)-SITLAF(NFLEVG-1))) 

     Z0  = SIRSLP*YDVAB%VRATH(NFLEVG)*YDVAB%VRATH(NFLEVG-1)*SIWEIG(NFLEVG,2)*(1.0_JPRB-SIWEIG(NFLEVG-1,3))
     Z1  = SIRSLP*YDVAB%VRATH(NFLEVG-1)*(YDVAB%VRATH(NFLEVG-1)*(SIWEIG(NFLEVG,1)*(1.0_JPRB-SIWEIG(NFLEVG-1,3)) &
         & +SIWEIG(NFLEVG-1,3)*SIWEIG(NFLEVG-1,2))-YDVAB%VRATH(NFLEVG-2)*SIWEIG(NFLEVG-1,2)*(1.0_JPRB-SIWEIG(NFLEVG-2,3)))
     Z2  = SIRSLP*YDVAB%VRATH(NFLEVG-2)*(YDVAB%VRATH(NFLEVG-2)*(SIWEIG(NFLEVG-1,1)*(1.0_JPRB-SIWEIG(NFLEVG-2,3))&
         & +SIWEIG(NFLEVG-2,3)*SIWEIG(NFLEVG-2,2))-YDVAB%VRATH(NFLEVG-1)*SIWEIG(NFLEVG-1,1)*SIWEIG(NFLEVG-1,3))
     Z3  = SIRSLP*YDVAB%VRATH(NFLEVG-2)*YDVAB%VRATH(NFLEVG-3)*SIWEIG(NFLEVG-2,1)*SIWEIG(NFLEVG-2,3)
     Z4  =-ZQBBC*(SIRSLP*YDVAB%VRATH(NFLEVG-1)*YDVAB%VRATH(NFLEVG)*SIWEIG(NFLEVG,1)/ &
         &(1.0_JPRB+SIRSLP*YDVAB%VRATH(NFLEVG)*YDVAB%VRATH(NFLEVG)*SIWEIG(NFLEVG,2)))
     Z5 = Z4 + (1.0_JPRB-ZQBBC)*(1.0_JPRB-(SITLAF(NFLEVG-1)/SITLAF(NFLEVG)))

     ZES(NFLEVG-1)= Z3*SITLAF(NFLEVG-3)/(SILNPR(NFLEVG-1)*(SITLAF(NFLEVG-2)-SITLAF(NFLEVG-3)))
     ZAS(NFLEVG-1)= Z2*SITLAF(NFLEVG-2)/(SILNPR(NFLEVG-1)*(SITLAF(NFLEVG-1)-SITLAF(NFLEVG-2)))  &
                  &-Z3* SITLAF(NFLEVG-2)/(SILNPR(NFLEVG-1)*(SITLAF(NFLEVG-2)-SITLAF(NFLEVG-3))) 
     ZBS(NFLEVG-1)=-(Z1+Z4*Z0)*SITLAF(NFLEVG-1)/(SILNPR(NFLEVG-1)*(SITLAF(NFLEVG)-SITLAF(NFLEVG-1))) &
                  &-Z2*SITLAF(NFLEVG-1)/(SILNPR(NFLEVG-1)*(SITLAF(NFLEVG-1)-SITLAF(NFLEVG-2)))
     ZCS(NFLEVG-1)=(Z1+Z5*Z0)*SITLAF(NFLEVG)/(SILNPR(NFLEVG-1)*(SITLAF(NFLEVG)-SITLAF(NFLEVG-1))) 
     ZDS(NFLEVG-1)= 0.0_JPRB    
 
     DO JLON = 1, KNLON
       IBOT1 = 1 + (NFLEVG-3)*KLEV + (JLON-1)*KLON
       IBOT2 = IBOT1 + KLEV
       IBOT3 = IBOT2 + KLEV
       IBOT0 = IBOT1 - KLEV
       PV2(IBOT2) = ZA(NFLEVG-1)*(PV1(IBOT1)-PV1(IBOT2))   &
                  & + ZC(NFLEVG-1)*(PV1(IBOT3)-PV1(IBOT2)) &
                  & + ZES(NFLEVG-1)*PV1(IBOT0) + ZAS(NFLEVG-1)*PV1(IBOT1) &
                  & + ZBS(NFLEVG-1)*PV1(IBOT2) + ZCS(NFLEVG-1)*PV1(IBOT3) 
     ENDDO
     
   ! * Compute the coefficients "A" and "C" at full level nflevg:
   !   "C" must be zero.
     ZA(NFLEVG)=SITLAF(NFLEVG-1)/(SILNPR(NFLEVG)*(SITLAF(NFLEVG)-SITLAF(NFLEVG-1)))
     ZC(NFLEVG)=0.0_JPRB

     ZCS(NFLEVG)= 0.0_JPRB
     ZDS(NFLEVG)= 0.0_JPRB

     Z0  =-SIRSLP*YDVAB%VRATH(NFLEVG)*YDVAB%VRATH(NFLEVG-1)*SIWEIG(NFLEVG,2)*(1.0_JPRB-SIWEIG(NFLEVG-1,3))
     Z1  = SIRSLP*YDVAB%VRATH(NFLEVG-1)*YDVAB%VRATH(NFLEVG-1)*(SIWEIG(NFLEVG,1)*(1.0_JPRB-SIWEIG(NFLEVG-1,3)) &
         & +SIWEIG(NFLEVG-1,3)*SIWEIG(NFLEVG-1,2))
     Z2  = SIRSLP*YDVAB%VRATH(NFLEVG-1)*YDVAB%VRATH(NFLEVG-2)*SIWEIG(NFLEVG-1,1)*SIWEIG(NFLEVG-1,3)
     
     ZAS(NFLEVG)=(Z1-Z4*Z0)*SITLAF(NFLEVG-1)/(SILNPR(NFLEVG)*(SITLAF(NFLEVG)-SITLAF(NFLEVG-1))) &
                &-Z2*SITLAF(NFLEVG-1)/(SILNPR(NFLEVG)*(SITLAF(NFLEVG-1)-SITLAF(NFLEVG-2)))
     ZBS(NFLEVG)=(Z5*Z0-Z1)*SITLAF(NFLEVG)/(SILNPR(NFLEVG)*(SITLAF(NFLEVG)-SITLAF(NFLEVG-1)))
     ZES(NFLEVG)=Z2*SITLAF(NFLEVG-2)/(SILNPR(NFLEVG)*(SITLAF(NFLEVG-1)-SITLAF(NFLEVG-2))) 

     ! * Bottom:
     DO JLON = 1, KNLON
       IBOT1 = 1 + (NFLEVG-2)*KLEV + (JLON-1)*KLON
       IBOT2 = IBOT1 + KLEV
       IBOT0 = IBOT1 - KLEV
       PV2(IBOT2) = ZA(NFLEVG)*(PV1(IBOT1)-PV1(IBOT2))&
        & -ZA(NFLEVG)*(SITLAF(NFLEVG)/SITLAF(NFLEVG-1) - 1.0_JPRB)*PV1(IBOT2) &
        & +ZES(NFLEVG)*PV1(IBOT0)+ZAS(NFLEVG)*PV1(IBOT1)+ZBS(NFLEVG)*PV1(IBOT2)
     ENDDO

  ELSE ! NOT LNHEE_REFINE_SILAPL

  ! * Compute the coefficients "A" and "C" at full level 1:
  !   "A" must be zero if the top pressure of the model is zero.
     ZPRTOP=YDVAB%VAH(0)+YDVAB%VBH(0)*SIPR
     ZA(1)=ZPRTOP/(SILNPR(1)*(SITLAF(1)-ZPRTOP))
     ZC(1)=SITLAF(2)/(SILNPR(1)*(SITLAF(2)-SITLAF(1)))
   
  ! * Compute the coefficients "A" and "C" at full level nflevg:
  !   "C" must be zero.
     ZA(NFLEVG)=SITLAF(NFLEVG-1)/(SILNPR(NFLEVG)*(SITLAF(NFLEVG)-SITLAF(NFLEVG-1)))
     ZC(NFLEVG)=0.0_JPRB
  
  !* Top:
     DO JLON = 1, KNLON
       ITOP1 = 1 + (JLON-1)*KLON
       ITOP2 = ITOP1 + KLEV
       PV2(ITOP1) = -ZA(1)*PV1(ITOP1)+ZC(1)*(PV1(ITOP2)-PV1(ITOP1))
     ENDDO

  ! * Bottom:
     DO JLON = 1, KNLON
       IBOT1 = 1 + (NFLEVG-2)*KLEV + (JLON-1)*KLON
       IBOT2 = IBOT1 + KLEV
       PV2(IBOT2) = ZA(NFLEVG)*(PV1(IBOT1)-PV1(IBOT2))&
        & - ZA(NFLEVG)*(SITLAF(NFLEVG)/SITLAF(NFLEVG-1) - 1.0_JPRB)*PV1(IBOT2)
     ENDDO  
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       2.    INNER LAYERS (AND BOUNDARY CONDITIONS FOR VFE)
!              ----------------------------------------------

IF(LVERTFE)THEN

  ! * precompute pi(l)/[d pi/d eta](l):
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV)
  DO JLEV=1,NFLEVG
    ZRM(JLEV) = SIRDEL(JLEV)/YDVETA%VFE_RDETAH(JLEV)
    ZA(JLEV)  = (SITLAF(JLEV)**2)*ZRM(JLEV)
    ZC(JLEV)  = ZA(JLEV)*ZRM(JLEV)
  ENDDO
!$OMP END PARALLEL DO

  ! * transform PV1 into form appropriate for VERDISINT:
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JLON,IDT)
  DO JLEV=1,NFLEVG
    DO JLON=1,KNLON
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      ZIN1(JLON,JLEV)=PV1(IDT) ! array to be derivated
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

  ! * term: 1/[d pi/d eta] * d/deta ( pi^2/[d pi/d eta] ) * ( d X/deta )
  ! and
  ! * term: (pi/[d pi/d eta])^2  d^2 X/deta^2
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JLON,IDT)
  DO JLON=1,KNLON
    ZIN1(JLON,0       )=0.0_JPRB
    ZIN1(JLON,NFLEVG+1)=0.0_JPRB
  ENDDO
!$OMP END PARALLEL DO
  CALL VERDISINT(YDVFE,'FDER','11',KNLON,1,KNLON,NFLEVG,ZIN1,ZINDETA)
  CALL VERDISINT(YDVFE,'DDER','01',KNLON,1,KNLON,NFLEVG,ZIN1,ZINDETA2)

  ! * d/deta (p^2/[d pi/d eta])
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JLON)
  DO JLEV=1,NFLEVG
    DO JLON=1,KNLON
      ZIN1(JLON,JLEV)=ZA(JLEV)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO
  ZIN1(1:KNLON,0)=0.0_JPRB
  ZIN1(1:KNLON,NFLEVG+1)=2.0_JPRB*SIPR
  CALL VERDISINT(YDVFE,'FDER','01',KNLON,1,KNLON,NFLEVG,ZIN1,ZP2MDETA)

  ! * compute PV2 = 2 * pi_l/[d pi/d eta]_l * dX/deta 
  !   + (pi_l/[d pi/d eta]_l)**2 * d^2 X/deta^2
  IF(LVFE_LAPL_BC)THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JLON,IDT)
    DO JLEV=1,NFLEVG
      DO JLON=1,KNLON
        IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
        PV2(IDT)= ZRM(JLEV)*ZP2MDETA(JLON,JLEV)*ZINDETA(JLON,JLEV) +&
         & ZC(JLEV)*ZINDETA2(JLON,JLEV)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ELSE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JLON,IDT)
    DO JLEV=2,NFLEVG-1
      DO JLON=1,KNLON
        IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
        PV2(IDT)= ZRM(JLEV)*ZP2MDETA(JLON,JLEV)*ZINDETA(JLON,JLEV) +&
         & ZC(JLEV)*ZINDETA2(JLON,JLEV)
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDIF

ELSE

  ! --- Finite differences ---
  IF (YDDYNA%LNHEE_REFINE_SILAPL) THEN
    ! * Compute the coefficients A and C at full levels 2 to nflevg-1:
    DO JLEV = 3, NFLEVG-2
       ZA(JLEV) = SITLAF(JLEV-1)/(SILNPR(JLEV)*(SITLAF(JLEV)-SITLAF(JLEV-1)))
       ZC(JLEV) = SITLAF(JLEV+1)/(SILNPR(JLEV)*(SITLAF(JLEV+1)-SITLAF(JLEV)))
    ENDDO

    DO JLEV=3,NFLEVG-2
       Z1 = SIRSLP*YDVAB%VRATH(JLEV+1)*YDVAB%VRATH(JLEV)*SIWEIG(JLEV+1,2)*(1.0_JPRB-SIWEIG(JLEV,3))
       Z2 = SIRSLP*YDVAB%VRATH(JLEV)*(YDVAB%VRATH(JLEV)*(SIWEIG(JLEV+1,1)*(1.0_JPRB-SIWEIG(JLEV,3)) &
        & +SIWEIG(JLEV,3)*SIWEIG(JLEV,2))-YDVAB%VRATH(JLEV-1)*SIWEIG(JLEV,2)*(1.0_JPRB-SIWEIG(JLEV-1,3))) 
       Z3 = SIRSLP*YDVAB%VRATH(JLEV-1)*(YDVAB%VRATH(JLEV-1)*(SIWEIG(JLEV,1)*(1.0_JPRB-SIWEIG(JLEV-1,3)) &
        & +SIWEIG(JLEV-1,3)*SIWEIG(JLEV-1,2))-YDVAB%VRATH(JLEV)*SIWEIG(JLEV,1)*SIWEIG(JLEV,3))
       Z4 = SIRSLP*YDVAB%VRATH(JLEV-1)*YDVAB%VRATH(JLEV-2)*SIWEIG(JLEV-1,1)*SIWEIG(JLEV-1,3)
      
       ZES(JLEV) =  Z4*SITLAF(JLEV-2)/(SILNPR(JLEV)*(SITLAF(JLEV-1)-SITLAF(JLEV-2)))
       ZAS(JLEV) =  Z3*SITLAF(JLEV-1)/(SILNPR(JLEV)*(SITLAF(JLEV)-SITLAF(JLEV-1)))  &
                 & -Z4*SITLAF(JLEV-1)/(SILNPR(JLEV)*(SITLAF(JLEV-1)-SITLAF(JLEV-2)))
       ZBS(JLEV) = -Z2*SITLAF(JLEV)/(SILNPR(JLEV)*(SITLAF(JLEV+1)-SITLAF(JLEV)))    &
                 & -Z3*SITLAF(JLEV)/(SILNPR(JLEV)*(SITLAF(JLEV)-SITLAF(JLEV-1)))  
       ZCS(JLEV) =  Z2*SITLAF(JLEV+1)/(SILNPR(JLEV)*(SITLAF(JLEV+1)-SITLAF(JLEV)))  &
                 & -Z1*SITLAF(JLEV+1)/(SILNPR(JLEV)*(SITLAF(JLEV+2)-SITLAF(JLEV+1)))
       ZDS(JLEV) =  Z1*SITLAF(JLEV+2)/(SILNPR(JLEV)*(SITLAF(JLEV+2)-SITLAF(JLEV+1)))
    ENDDO

  ! * Other Levels:
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JLON,ILEV0,ILEV1,ILEV2,ILEV3,ILEV4)
    DO JLEV = 3, NFLEVG-2
      DO JLON = 1, KNLON
        ILEV0 = 1 + (JLEV-3)*KLEV + (JLON-1)*KLON
        ILEV1 = ILEV0 + KLEV
        ILEV2 = ILEV1 + KLEV
        ILEV3 = ILEV2 + KLEV
        ILEV4 = ILEV3 + KLEV
        PV2(ILEV2) = ZA(JLEV)*(PV1(ILEV1)-PV1(ILEV2))    &
         & + ZC(JLEV)*(PV1(ILEV3)-PV1(ILEV2))            &
         & + ZAS(JLEV)*PV1(ILEV1) + ZBS(JLEV)*PV1(ILEV2) &
         & + ZCS(JLEV)*PV1(ILEV3) + ZDS(JLEV)*PV1(ILEV4) &
         & + ZES(JLEV)*PV1(ILEV0)  
      ENDDO
    ENDDO    
!$OMP END PARALLEL DO

  ELSE !NOT LNHEE_REFINE_SILAPL

  ! * Compute the coefficients A and C at full levels 2 to nflevg-1:
    DO JLEV = 2, NFLEVG-1
       ZA(JLEV) = SITLAF(JLEV-1)/(SILNPR(JLEV)*(SITLAF(JLEV)-SITLAF(JLEV-1)))
       ZC(JLEV) = SITLAF(JLEV+1)/(SILNPR(JLEV)*(SITLAF(JLEV+1)-SITLAF(JLEV)))
    ENDDO

  ! * Other Levels:
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JLON,ILEV1,ILEV2,ILEV3)
    DO JLEV = 2, NFLEVG-1
      DO JLON = 1, KNLON
        ILEV1 = 1 + (JLEV-2)*KLEV + (JLON-1)*KLON
        ILEV2 = ILEV1 + KLEV
        ILEV3 = ILEV2 + KLEV
        PV2(ILEV2) = ZA(JLEV)*(PV1(ILEV1)-PV1(ILEV2))&
         & + ZC(JLEV)*(PV1(ILEV3)-PV1(ILEV2))
      ENDDO
    ENDDO
!$OMP END PARALLEL DO
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       3.    MULTIPLY BY SITR/SITRAM
!              ----------------------

IF(PRESENT(LD_LSTAR)) THEN
  LL_LSTAR=LD_LSTAR
ELSE
  LL_LSTAR=.FALSE.
ENDIF

IF(.NOT.LL_LSTAR) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JLEV,JLON,IDT)
  DO JLEV=1,NFLEVG
    ZRATIO_SITR(JLEV)=SITR/SITRAM(JLEV)
    DO JLON=1,KNLON
      IDT=1+(JLEV-1)*KLEV+(JLON-1)*KLON
      PV2(IDT)=ZRATIO_SITR(JLEV)*PV2(IDT)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SISEVE',1,ZHOOK_HANDLE)
END SUBROUTINE SISEVE
