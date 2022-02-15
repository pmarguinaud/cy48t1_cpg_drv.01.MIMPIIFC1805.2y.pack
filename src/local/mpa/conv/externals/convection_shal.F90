!-------------------------------------------------------------------------------
!   ############################################################################
SUBROUTINE CONVECTION_SHAL( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,&
 & PDTCONV, LD_OREFRESH_ALL, LD_ODOWN, KICE,&
 & LD_OSETTADJ, PTADJS, LSMOOTH, &
 & PA25, PCRAD, PCDEPTH, PCDEPTH_D, PDTPERT, PATPERT, PBTPERT, PENTR, &
 & PZLCL, PZPBL, PWTRIG, PNHGAM, PTFRZ1, PTFRZ2, PSTABT, PSTABC, PAW, PBW, &
 & PPABS, PZZ, PTKECLS,&
 & PT, PRV, PRC, PRI, PU, PV, PW,&
 & KCOUNT, PTTEN, PRVTEN, PRCTEN, PRITEN,&
 & PUMF, PCAPE, KCLTOP, KCLBAS,&
 & PURV, PURCI,&
 & LD_OUVTRANS, PUTEN, PVTEN,&
 & LD_OCHTRANS, KCH1, PCH1, PCH1TEN )
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!   ############################################################################

!!**** Interface routine to the fast Meso-NH shallow-convection code developed for ECMWF/ARPEGE
!!     having a structure typical for operational routines
!!     
!!     Transformations necessary to call shallow code
!!     - skip input vertical arrays/levels : bottom=1, top=KLEV
!!     - transform specific humidities in mixing ratio
!!
!!
!!    PURPOSE
!!    -------
!!      The routine interfaces the MNH convection code as developed for operational
!!      forecast models like ECMWF/ARPEGE or HIRLAM with the typical Meso-NH array structure
!!      Calls the deep and/or shallow convection routine
!!
!!
!!**  METHOD
!!    ------
!!     Returns one tendency for shallow+deep convection but each part can
!!     be activated/desactivated separately
!!     For deep convection one can enable up to 3 additional ensemble members
!!     - this substantially improves the smoothness of the scheme and reduces
!!       allows for runs with different cloud radii (entrainment rates) and
!!       reduces the arbitrariness inherent to convective trigger condition
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!    CONVECT_SHALLOW
!!    SU_CONVPAR, SU_CONVPAR1
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/12/98
!!      modified    20/03/2002 by P. Marquet : transformed for ARPEGE/Climat
!!                             (tsmbkind.h, REAL_B, INTEGER_M, , "& &",
!!                              _ZERO_, _ONE_, _HALF_)
!!      modified    11/04/O2 allow for ensemble of deep updrafts/downdrafts
!!      modified    20/07/2009 : E. Bazile Parameter for Wlcl XAW, XBW, for
!!                         the temp. pert. and LSMOOTH. Default value are:
!!                         XAW=0. XBW=1. XATPERT=0. XBTPERT=1. LSMOOTH=T
!!      modified    30/01/2012 : F. Bouyssel. Remove Meso-Nh initialization
!!                         inducing non-reproducibility under Open-MP
!!
!!    REFERENCE
!!    ---------
!!    Bechtold et al., 2001, Quart. J. Roy. Meteor. Soc., Vol 127, pp 869-886: 
!!           A mass flux convection scheme for regional and global models.
!!
!-------------------------------------------------------------------------------

USE MODD_CONVPAR_SHAL, ONLY : LLSMOOTH, XA25, XATPERT, XAW, XBTPERT, XBW, XCDEPTH, XCDEPTH_D, XCRAD, XDTPERT, XENTR, &
& XNHGAM, XSTABC, XSTABT, XTFRZ1, XTFRZ2, XWTRIG, XZLCL, XZPBL

IMPLICIT NONE


!*       0.1   Declarations of dummy arguments :

INTEGER, INTENT(IN)    :: KLON ! horizontal dimension
INTEGER, INTENT(IN)    :: KLEV ! vertical dimension
INTEGER, INTENT(IN)    :: KIDIA ! value of the first point in x
INTEGER, INTENT(IN)    :: KFDIA ! value of the last point in x
INTEGER, INTENT(IN)    :: KBDIA ! vertical  computations start at
INTEGER, INTENT(IN)    :: KTDIA ! vertical computations can be
REAL   ,INTENT(IN)    :: PDTCONV ! Interval of time between two
LOGICAL           ,INTENT(IN)    :: LD_OREFRESH_ALL ! refresh or not all 
LOGICAL           ,INTENT(IN)    :: LD_ODOWN ! take or not convective
INTEGER, INTENT(IN)    :: KICE ! flag for ice ( 1 = yes, 
LOGICAL           ,INTENT(IN)    :: LD_OSETTADJ ! logical to set convective
LOGICAL, INTENT(IN)    :: LSMOOTH ! Supposed to be necessary ...
REAL    ,INTENT(IN)    :: PTADJS ! user defined shal. adjustment time (s)
REAL    ,INTENT(IN)    :: PA25, PCRAD, PCDEPTH, PCDEPTH_D, PDTPERT, &
         & PATPERT, PBTPERT, &
         & PENTR, PZLCL, PZPBL, PWTRIG, PNHGAM, PTFRZ1, PTFRZ2, &
         & PSTABT, PSTABC, PAW, PBW
REAL    ,INTENT(IN)    :: PPABS(KLON,KLEV) ! grid scale pressure (Pa)
REAL    ,INTENT(IN)    :: PZZ(KLON,KLEV) ! geopotential (m2/s2) 
REAL    ,INTENT(IN)    :: PTKECLS(KLON) ! TKE in the CLS (m2/s2) 
REAL    ,INTENT(IN)    :: PT(KLON,KLEV) ! grid scale T at time t  (K)
REAL    ,INTENT(IN)    :: PRV(KLON,KLEV) ! grid scale water vapor  (kg/kg)
REAL    ,INTENT(IN)    :: PRC(KLON,KLEV) ! grid scale r_c (kg/kg)
REAL    ,INTENT(IN)    :: PRI(KLON,KLEV) ! grid scale r_i (kg/kg)
REAL    ,INTENT(IN)    :: PU(KLON,KLEV) ! grid scale horiz. wind u (m/s) 
REAL    ,INTENT(IN)    :: PV(KLON,KLEV) ! grid scale horiz. wind v (m/s)
REAL    ,INTENT(IN)    :: PW(KLON,KLEV) ! grid scale vertical velocity (m/s)
INTEGER, INTENT(INOUT) :: KCOUNT(KLON) ! convective counter(recompute
REAL    ,INTENT(INOUT) :: PTTEN(KLON,KLEV) ! convective temperat. tendency (K/s)
REAL    ,INTENT(INOUT) :: PRVTEN(KLON,KLEV) ! convective r_v tendency (1/s)
REAL    ,INTENT(INOUT) :: PRCTEN(KLON,KLEV) ! convective r_c tendency (1/s)
REAL    ,INTENT(INOUT) :: PRITEN(KLON,KLEV) ! convective r_i tendency (1/s)
REAL    ,INTENT(INOUT) :: PUMF(KLON,KLEV) ! updraft mass flux   (kg/s m2)
REAL    ,INTENT(OUT)   :: PCAPE(KLON) ! CAPE (J/kg)
INTEGER, INTENT(INOUT) :: KCLTOP(KLON) ! cloud top level (number of model level)
INTEGER, INTENT(INOUT) :: KCLBAS(KLON) ! cloud base level(number of model level)
REAL    ,INTENT(INOUT) :: PURV(KLON,KLEV) ! water vapor in updraft (kg/kg)
REAL    ,INTENT(INOUT) :: PURCI(KLON,KLEV) ! total condensate in updraft (kg/kg)
LOGICAL           ,INTENT(IN)    :: LD_OUVTRANS ! flag to compute convective
REAL    ,INTENT(INOUT) :: PUTEN(KLON,KLEV) ! convecctive u tendency (m/s^2)
REAL    ,INTENT(INOUT) :: PVTEN(KLON,KLEV) ! convecctive v tendency (m/s^2)
LOGICAL           ,INTENT(IN)    :: LD_OCHTRANS ! flag to compute convective
INTEGER, INTENT(IN)    :: KCH1 ! number of species
REAL    ,INTENT(IN)    :: PCH1(KLON,KLEV,KCH1) ! grid scale chemical species
REAL    ,INTENT(INOUT) :: PCH1TEN(KLON,KLEV,KCH1) ! chemical convective tendency
!*       0.2   Declarations of local variables :

INTEGER  :: JI, JK, JKP, JN  ! loop index

! Local arrays (upside/down) necessary for change of ECMWF arrays to convection arrays
REAL , DIMENSION(KLON,KLEV) :: ZT     ! grid scale T at time t  (K)
REAL , DIMENSION(KLON,KLEV) :: ZRV    ! grid scale water vapor  (kg/kg)
REAL , DIMENSION(KLON,KLEV) :: ZRC    ! grid scale r_c mixing ratio (kg/kg)
REAL , DIMENSION(KLON,KLEV) :: ZRI    ! grid scale r_i mixing ratio (kg/kg)
REAL , DIMENSION(KLON,KLEV) :: ZU     ! grid scale horiz. wind u (m/s) 
REAL , DIMENSION(KLON,KLEV) :: ZV     ! grid scale horiz. wind v (m/s)
REAL , DIMENSION(KLON,KLEV) :: ZW     ! grid scale vertical velocity (m/s)
REAL , DIMENSION(KLON,KLEV) :: ZPABS  ! grid scale pressure (Pa)
REAL , DIMENSION(KLON,KLEV) :: ZZZ    ! height of model layer (m) 

REAL , DIMENSION(KLON,KLEV) :: ZTTEN  ! convective temperat. tendency (K/s)
REAL , DIMENSION(KLON,KLEV) :: ZRVTEN ! convective r_v tendency (1/s)
REAL , DIMENSION(KLON,KLEV) :: ZRCTEN ! convective r_c tendency (1/s)
REAL , DIMENSION(KLON,KLEV) :: ZRITEN ! convective r_i tendency (1/s)
REAL , DIMENSION(KLON,KLEV) :: ZUTEN  ! convective u tendency (m/s^2)
REAL , DIMENSION(KLON,KLEV) :: ZVTEN  ! convective m tendency (m/s^2)
REAL , DIMENSION(KLON,KLEV) :: ZUMF   ! updraft mass flux   (kg/s m2)
REAL , DIMENSION(KLON,KLEV) :: ZURV   ! water vapor in updrafts (kg/kg)
REAL , DIMENSION(KLON,KLEV) :: ZURCI  ! total condensate in updrafts (kg/kg)
INTEGER,  DIMENSION(KLON)   :: ICLTOP ! cloud top level (number of model level)
INTEGER,  DIMENSION(KLON)   :: ICLBAS ! cloud base level(number of model level)
REAL , DIMENSION(KLON,KLEV,KCH1):: ZCH1     ! grid scale chemical species
REAL , DIMENSION(KLON,KLEV,KCH1):: ZCH1TEN  ! chemical convective tendency

! special for shallow convection
REAL , DIMENSION(:,:), ALLOCATABLE   :: ZTTENS, ZRVTENS, ZRCTENS, ZRITENS, &
 & ZUMFS
REAL , DIMENSION(:,:,:), ALLOCATABLE :: ZCH1TENS
INTEGER,  DIMENSION(:), ALLOCATABLE  :: ICLBASS, ICLTOPS


!*       1.   Allocate 2D (horizontal, vertical) arrays and additional ensemble arrays
!             ------------------------------------------------------------------------

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CONVECTION_SHAL',0,ZHOOK_HANDLE)
ALLOCATE( ZTTENS(KLON,KLEV) ) 
ALLOCATE( ZRVTENS(KLON,KLEV) ) 
ALLOCATE( ZRCTENS(KLON,KLEV) )
ALLOCATE( ZRITENS(KLON,KLEV) ) 
ALLOCATE( ZCH1TENS(KLON,KLEV,KCH1) ) 
ALLOCATE( ZUMFS(KLON,KLEV) )
ALLOCATE( ICLBASS(KLON) )
ALLOCATE( ICLTOPS(KLON) )


KCLTOP(:)  = 1 ! set default value when no convection
KCLBAS(:)  = 1 ! can be changed  depending on user
ICLTOP(:)  = 1 
ICLBAS(:)  = 1 
ICLTOPS(:) = 1 
ICLBASS(:) = 1 


!*       2.   Flip arrays upside-down as  first vertical level in convection is 1
!             --------------------------------------------------------------------

DO JK = 1, KLEV
  JKP = KLEV - JK + 1
  DO JI = KIDIA, KFDIA
    ZPABS(JI,JKP) = PPABS(JI,JK)
    ZZZ(JI,JKP)   = PZZ(JI,JK)
    ZT(JI,JKP)    = PT(JI,JK)
    ZRV(JI,JKP)   = PRV(JI,JK) / ( 1.0 - PRV(JI,JK) ) ! transform specific humidity
    ZRC(JI,JKP)   = PRC(JI,JK) / ( 1.0 - PRC(JI,JK) ) ! in mixing ratio
    ZRI(JI,JKP)   = PRI(JI,JK) / ( 1.0 - PRI(JI,JK) ) 
    ZU(JI,JKP)    = PU(JI,JK)
    ZV(JI,JKP)    = PV(JI,JK)
    ZW(JI,JKP)    = PW(JI,JK) 
  ENDDO
ENDDO
IF ( LD_OCHTRANS ) THEN
  DO JK = 1, KLEV
    JKP = KLEV - JK + 1
    DO JN = 1, KCH1
      DO JI = KIDIA, KFDIA
        ZCH1(JI,JKP,JN) = PCH1(JI,JK,JN)
      ENDDO
    ENDDO
  ENDDO
ENDIF

  KCOUNT(:)     =0
  ZTTEN(:,:)    =0.0
  ZRVTEN(:,:)   =0.0
  ZRCTEN(:,:)   =0.0
  ZRITEN(:,:)   =0.0
  ZUTEN(:,:)    =0.0
  ZVTEN(:,:)    =0.0
  ZUMF(:,:)     =0.0
  ZURV(:,:)     =0.0
  ZURCI(:,:)    =0.0
  ZCH1TEN(:,:,:)=0.0
  PCAPE(:)      =0.0


!*       4.b  Call shallow convection routine
!             -------------------------------

  CALL INI_CONVPAR 
!  CALL INI_CONVPAR_SHAL
  XA25=PA25
  XCRAD=PCRAD
  XCDEPTH=PCDEPTH
  XCDEPTH_D=PCDEPTH_D
  XDTPERT=PDTPERT
  XATPERT=PATPERT
  XBTPERT=PBTPERT
  XENTR=PENTR
  XZLCL=PZLCL
  XZPBL=PZPBL
  XWTRIG=PWTRIG
  XNHGAM=PNHGAM
  XTFRZ1=PTFRZ1
  XTFRZ2=PTFRZ2
  XSTABT=PSTABT
  XSTABC=PSTABC
  XAW=PAW
  XBW=PBW
  LLSMOOTH=LSMOOTH

  CALL SHALLOW_CONVECTION( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,        &
   & PDTCONV, KICE, LD_OSETTADJ, PTADJS,         &
   & ZPABS, ZZZ,PTKECLS,                         &
   & ZT, ZRV, ZRC, ZRI, ZW,                      &
   & ZTTENS, ZRVTENS, ZRCTENS, ZRITENS,          &
   & ICLTOPS, ICLBASS, ZUMFS, &
   & LD_OCHTRANS, KCH1, ZCH1, ZCH1TENS )

DO JK = 1, KLEV
  DO JI = KIDIA, KFDIA
    ZTTEN(JI,JK) = ZTTENS(JI,JK)
    ZRVTEN(JI,JK) = ZRVTENS(JI,JK)
    ZRCTEN(JI,JK) = ZRCTENS(JI,JK)
    ZRITEN(JI,JK) = ZRITENS(JI,JK)
    ZUMF(JI,JK)   = ZUMFS(JI,JK)
  ENDDO
ENDDO
DO JI = KIDIA, KFDIA
  ICLTOP(JI)   = MAX(ICLTOP(JI), ICLTOPS(JI))
  ICLBAS(JI)   = MAX(ICLBAS(JI), ICLBASS(JI))
ENDDO

!*       6.  Reflip arrays to ECMWF/ARPEGE vertical structure
!            change mixing ratios to sepcific humidity

DO JK = 1, KLEV
  JKP = KLEV - JK + 1
  DO JI = KIDIA, KFDIA
    PTTEN(JI,JK)  = ZTTEN(JI,JKP)
   ! don't transform back to specific hum, does not conserve integrals
    PRVTEN(JI,JK) = ZRVTEN(JI,JKP) ! / ( 1.0 + ZRV(JI,JKP) ) ** 2
    PRCTEN(JI,JK) = ZRCTEN(JI,JKP) ! / ( 1.0 + ZRC(JI,JKP) ) ** 2
    PRITEN(JI,JK) = ZRITEN(JI,JKP) ! / ( 1.0 + ZRI(JI,JKP) ) ** 2
    PUTEN(JI,JK)  = ZUTEN(JI,JKP)
    PVTEN(JI,JK)  = ZVTEN(JI,JKP)
    PUMF(JI,JK)   = ZUMF(JI,JKP)
    PURV(JI,JK)   = ZURV(JI,JKP) / ( 1.0 + ZURV(JI,JKP) )
    PURCI(JI,JK)  = ZURCI(JI,JKP)/ ( 1.0 + ZURCI(JI,JKP) )
  ENDDO
ENDDO

DO JI = KIDIA, KFDIA
  JK = ICLTOP(JI)
  KCLTOP(JI) = KLEV - JK + 1
  JK = ICLBAS(JI)
  KCLBAS(JI) = KLEV - JK + 1
  IF ( ICLTOP(JI) == 1 ) KCLTOP(JI) = 1
  IF ( ICLBAS(JI) == 1 ) KCLBAS(JI) = 1
ENDDO

IF ( LD_OCHTRANS ) THEN
  DO JK = 1, KLEV
    JKP = KLEV - JK + 1
    DO JN = 1, KCH1
      DO JI = KIDIA, KFDIA
        PCH1TEN(JI,JK,JN) = ZCH1TEN(JI,JKP,JN)
      ENDDO
    ENDDO
  ENDDO
ENDIF
   
!*       7.  Deallocate local arrays

DEALLOCATE( ICLBASS )
DEALLOCATE( ICLTOPS )
DEALLOCATE( ZUMFS )
DEALLOCATE( ZCH1TENS ) 
DEALLOCATE( ZRCTENS )
DEALLOCATE( ZRITENS ) 
DEALLOCATE( ZTTENS )
DEALLOCATE( ZRVTENS ) 

IF (LHOOK) CALL DR_HOOK('CONVECTION_SHAL',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECTION_SHAL

!----------------------------------------------------------------------------
