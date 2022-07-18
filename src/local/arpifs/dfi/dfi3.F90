SUBROUTINE DFI3(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,YDVARBC,LDIAB,KDIGLST,KDFI,LDIFH,PTSTEP,LDMPHYS,LDEPHYS,&
 & KNINDAT,KNSSSSS,PSPA3,PSPA2,PSPA1,YDTCV,PGFLGPQ,YDGOM5,YDODB,YDJOT)

!**** *DFI3*  - Controls integration for digital filter initialization
!                   on lowest level

!     Purpose.
!     --------
!     CONTROLS  DIGITAL FILTER INITIALIZATION ON LEVEL 3 (LOWEST LEVEL)

!**   Interface.
!     ----------
!        *CALL* *DFI3(...)*

!        Explicit arguments :
!        --------------------
!        LDIAB: controls if filter is diabatic
!        KDIGLST: number of integration steps in one direction
!        KDFI   : indicates type of integration session (backw.,
!                   forw./1 , forw./2 , forw./1+2 )
!        LDIFH: controls horizontal diffusion backwards
!        PSTEP  : time-step
!        LDMPHYS,LDEPHYS: original values of physics switches
!        KNINDAT,KNSSSSS: original starting model time (date,second)
!        PSPA3...: accumulation arrays of spectral state variables

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!     ARPEGE/ALADIN DOCUMENTATION

!     Author.
!     -------
!      Gabor Radnoti GMAP/MICECO
!      Original : 92-12-24

! Modifications
! -------------
!   F. Vana  13-Jan-2009 : Dummy update for SUHDU
!   A.Bogatchev&C.Fischer : 09-06-2009 restore accumulation in DF
!   K. Essaouini & R. El Khatib 23-Jul-2009 Fix restored accumulation in DF
!   K. Yessad (Dec 2009): CLCONF(4:4)='A' for corrector step too.
!   B. Bochenek (Oct 2013): CLCONF(7:7)='D' for DFI
!   K. Yessad (July 2014): Move some variables.
!   F. Vana   05-Mar-2015  Support for single precision
!   A. Geer   27-Jul-2015  VarBC is now an object passed by argument, for OOPS
!   O. Marsden   Aug 2016  Removed use of SPA3
!   K. Yessad (June 2017): Introduce NHQE model.
!   K. Yessad (June 2018): Alternate NHEE SI scheme elimination.
!   F. Suzat Add YDGOM optional arguments to STEPO 
! End Modifications
!------------------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE FIELDS_MOD         , ONLY : FIELDS
USE MTRAJ_MOD          , ONLY : MTRAJ
USE PARKIND1           , ONLY : JPRD, JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK
USE JO_TABLE_MOD       , ONLY : JO_TABLE
USE YOMCT0             , ONLY : LTWOTL, LELAM, LOBSC1, LNHEE, LNHQE
USE YOMCT2             , ONLY : NSTAR2, NSTOP2
USE YOMCT3             , ONLY : NSTEP
USE YOMINI             , ONLY : LSCRINI
USE YOMLUN             , ONLY : NULOUT, NULERR
USE YOMECTAB           , ONLY : MTSLOTNO
USE YOMJCDFI           , ONLY : RACCSPA3, RACCSPA2, RACCSPA1
USE YOMMP0             , ONLY : MYPROC, MYSETV
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)
USE VARBC_CLASS        , ONLY : CLASS_VARBC
USE TOVSCV_MOD         , ONLY : TOVSCV
USE SUPERGOM_CLASS     , ONLY : CLASS_SUPERGOM
USE DBASE_MOD          , ONLY : DBASE
USE YOMDYNA            , ONLY : YRDYNA

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(INOUT) :: YDGEOMETRY  !! INOUT needed for call to STEPO
TYPE(FIELDS)      ,INTENT(INOUT) :: YDFIELDS
TYPE(MTRAJ)       ,INTENT(INOUT) :: YDMTRAJ
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
TYPE(CLASS_VARBC) ,INTENT(INOUT) :: YDVARBC
LOGICAL           ,INTENT(IN)    :: LDIAB
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIGLST
INTEGER(KIND=JPIM),INTENT(INOUT) :: KDFI
LOGICAL           ,INTENT(IN)    :: LDIFH 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSTEP 
LOGICAL           ,INTENT(IN)    :: LDMPHYS 
LOGICAL           ,INTENT(IN)    :: LDEPHYS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNINDAT
INTEGER(KIND=JPIM),INTENT(IN)    :: KNSSSSS
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPA3(YDGEOMETRY%YRDIMV%NFLSUR,YDGEOMETRY%YRDIM%NSPEC2,YDMODEL%YRML_GCONF%YRDIMF%NS3D)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPA2(YDGEOMETRY%YRDIM%NSPEC2,YDMODEL%YRML_GCONF%YRDIMF%NFD2D)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPA1(YDGEOMETRY%YRDIMV%NFLEVL,2)
TYPE(TOVSCV)      ,INTENT(IN),    OPTIONAL :: YDTCV
REAL(KIND=JPRB)   ,INTENT(INOUT), OPTIONAL :: PGFLGPQ(:,:) 
TYPE(CLASS_SUPERGOM),INTENT(INOUT), OPTIONAL :: YDGOM5
CLASS(DBASE),       INTENT(INOUT),OPTIONAL :: YDODB
TYPE(JO_TABLE),     INTENT(INOUT),OPTIONAL :: YDJOT

!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: IACTIM, IDEC, ISIGN, ISTART, ISTEP, ISTOP
INTEGER(KIND=JPIM) :: JSITER, JST, JSTEP, IGPP, ISUB

LOGICAL :: LLWRITE, LLOBSC1, LLSLOT

REAL(KIND=JPRD) :: ZT1, ZT2
REAL(KIND=JPRB), ALLOCATABLE :: ZGFLGP(:,:)    !  Buffer to store midle value of grib point GFL
CHARACTER :: CLCONF*9, CLDAYF*39, CLTIMEOD*10, CL1*1, CL4*1, CL6*1, CL9*1
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "user_clock.h"

#include "abor1.intfb.h"
#include "chkobtim.intfb.h"
#include "digfil.intfb.h"
#include "digp.intfb.h"
#include "elsrw.intfb.h"
#include "reast.intfb.h"
#include "obsv.intfb.h"
#include "stepo.intfb.h"
#include "sueheg.intfb.h"
#include "suhdu.intfb.h"
#include "suheg.intfb.h"
#include "sunhsi.intfb.h"
#include "sunheesi.intfb.h"
#include "sunhqesi.intfb.h"
#include "updobs.intfb.h"
#include "updtim.intfb.h"
#include "sualldfi.intfb.h"
#include "dealldfi.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('DFI3',0,ZHOOK_HANDLE)
ASSOCIATE(YDGFL5=>YDMTRAJ%YRGFL5,YDGMV5=>YDMTRAJ%YRGMV5, YDGFL=>YDFIELDS%YRGFL,YDGMV=>YDFIELDS%YRGMV, &
 &  YDSURF=>YDFIELDS%YRSURF, YDDIM=>YDGEOMETRY%YRDIM, YDDIMV=>YDGEOMETRY%YRDIMV, &
 & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDDYN=>YDMODEL%YRML_DYN%YRDYN,YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY, &
 & YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YDEDYN=>YDMODEL%YRML_DYN%YREDYN, &
 & YGFL=>YDMODEL%YRML_GCONF%YGFL,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY,YDDIMF=>YDMODEL%YRML_GCONF%YRDIMF)

ASSOCIATE(NUMGPFLDS=>YGFL%NUMGPFLDS, &
 & NSPEC2=>YDDIM%NSPEC2, &
 & NFD2D=>YDDIMF%NFD2D, NS3D=>YDDIMF%NS3D, &
 & NFLEVG=>YDDIMV%NFLEVG, NFLEVL=>YDDIMV%NFLEVL, NFLSUR=>YDDIMV%NFLSUR, &
 & LRHDI_LASTITERPC=>YDDYN%LRHDI_LASTITERPC, LSIDG=>YDDYN%LSIDG, &
 & LSTRHD=>YDDYN%LSTRHD, NCURRENT_ITER=>YDDYN%NCURRENT_ITER, &
 & NSITER=>YDDYN%NSITER, &
 & LESIDG=>YDEDYN%LESIDG, &
 & LEPHYS=>YDEPHY%LEPHYS, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & NBSETSP=>YDMP%NBSETSP, &
 & NSTOP=>YDRIP%NSTOP, TDT=>YDRIP%TDT, TSTEP=>YDRIP%TSTEP, &
 & LMPHYS=>YDPHY%LMPHYS)
!     ------------------------------------------------------------------

!* 0.  Allocate ZGFLGP array
!      ---------------------------------

IGPP = NUMGPFLDS*NFLEVG
ALLOCATE(ZGFLGP(NGPTOT,IGPP))

!* 1.  Set up for different DFI sessions
!      ---------------------------------

CLCONF(1:1)='0'
CLCONF(2:2)='A'
CLCONF(3:3)='A'
CLCONF(4:4)='A'
CLCONF(5:5)='0'
CLCONF(6:6)='0'
CLCONF(7:7)='D'
CLCONF(8:8)='A'
! CLCONF(9:9)='A' if horizontal diffusion required, 'I' else

LLWRITE=.FALSE.
LLOBSC1=.FALSE.
IACTIM=1

!                              *********************************
!        1.1   DFI session:    *        BACKWARD               *
!                              *********************************

IF (KDFI == 1) THEN
  ISTART= 0
  TSTEP=-PTSTEP
  LMPHYS=.FALSE.
  LEPHYS=.FALSE.
  CLCONF(9:9)='I'
  IF (LDIAB) THEN
    ISTOP= KDIGLST-1
  ELSE
    ISTOP= KDIGLST
  ENDIF
!                              *********************************
!        1.2   DFI session:    *          FORWARD/1:           *
!                              *********************************
ELSEIF (KDFI == 2) THEN
  ISTART= 0
  TSTEP=PTSTEP
  LMPHYS=LDMPHYS
  LEPHYS=LDEPHYS
  CLCONF(9:9)='A'
  IF (LDIAB) THEN
    ISTOP= KDIGLST-1
  ELSE
    ISTOP= KDIGLST
  ENDIF
!                              *********************************
!        1.3   DFI session:    *        FORWARD/2:             *
!                              *********************************
ELSEIF (KDFI == 3) THEN
  ISTART= 0
  ISTOP= KDIGLST
  TSTEP=PTSTEP
  IF (LDIAB) THEN
    LMPHYS=LDMPHYS
    LEPHYS=LDEPHYS
    CLCONF(9:9)='A'
  ELSE
    LMPHYS=.FALSE.
    LEPHYS=.FALSE.
    CLCONF(9:9)='I'
  ENDIF
!                              *********************************
!        1.4   DFI session:    *        FORWARD/1 & 2:         *
!                              *********************************
ELSEIF (KDFI == 4) THEN
  ISTART= 0
  ISTOP= 2*KDIGLST
  TSTEP=PTSTEP
  IF (LDIAB) THEN
    LMPHYS=LDMPHYS
    LEPHYS=LDEPHYS
    CLCONF(9:9)='A'
  ELSE
    LMPHYS=.FALSE.
    LEPHYS=.FALSE.
    CLCONF(9:9)='I'
  ENDIF
  LLWRITE=.TRUE.
  LLOBSC1=LOBSC1 .AND. LSCRINI
!                              *********************************
!        1.5   DFI session:    *   ONE STEP FORWARD FOR MDDFI  *
!                              *********************************
ELSEIF (KDFI == 5) THEN
  ISTART= 0
  ISTOP= 0
  TSTEP=PTSTEP
  LMPHYS=.FALSE.
  LEPHYS=.FALSE.
  CLCONF(9:9)='I'
!                              *********************************
!        1.6   DFI session:    *   BACKWARD AND FORWARD        *
!                              *********************************
ELSEIF (KDFI == 6) THEN
  ISTART=0
  ISTOP=2*KDIGLST
  IF (LDIAB) THEN
    TSTEP=PTSTEP
    LMPHYS=LDMPHYS
    LEPHYS=LDEPHYS
    ISIGN=1
    CLCONF(9:9)='A'
  ELSE
    TSTEP=-PTSTEP
    LMPHYS=.FALSE.
    LEPHYS=.FALSE.
    ISIGN=-1
    CLCONF(9:9)='I'
  ENDIF

ELSE
  CALL ABOR1(' UNKNOWN OPTION FOR DFI3')
ENDIF
IF (LDIFH) CLCONF(9:9)='A'

CL1=CLCONF(1:1)
CL4=CLCONF(4:4)
CL6=CLCONF(6:6)
CL9=CLCONF(9:9)

CALL REAST(YDRIP,KDFI,KDIGLST,KNINDAT,KNSSSSS,IDEC)

!-----------------------------------------------------------------

!* 2.  Temporal loop
!      -------------

NSTAR2=ISTART
NSTOP2=ISTOP
IF (LLOBSC1) NSTOP =ISTOP

DO JSTEP=ISTART,ISTOP

  CALL USER_CLOCK(PTOTAL_CP=ZT1)

!        2.1  Filter accumulation if necessary

  JST = JSTEP
  ISUB=1
  IF (KDFI == 1)                                 JST=-JSTEP
  IF ((KDFI == 2.AND.LDIAB).OR.(KDFI == 4).OR.(KDFI == 6))JST=JSTEP-KDIGLST
  IF ((KDFI == 2.AND..NOT.LDIAB).OR.(KDFI == 5)) JST=JSTEP-KDIGLST-1
  IF ( (KDFI == 1.AND..NOT.LDIAB) .OR.(&
   & KDFI == 2.AND.(LDIAB.OR.JSTEP > 0)) .OR.(&
   & KDFI == 3.AND.(LDIAB.OR.JSTEP > 0)) .OR.(&
   & KDFI == 4) .OR. (KDFI == 6) ) THEN
    CALL SUALLDFI(YDGEOMETRY,YDDIMF,NULOUT)
  IF (PRESENT(PGFLGPQ)) THEN
      CALL DIGFIL(YDGEOMETRY,YDFIELDS%YRGFL,YDMODEL%YRML_GCONF,JST,KDIGLST,ISUB,&
           &      YDFIELDS%YRSPEC%SP3D,YDFIELDS%YRSPEC%SP2D,YDFIELDS%YRSPEC%SP1D,LDIAB,PGFLGPQ)
  ELSE
      CALL DIGFIL(YDGEOMETRY,YDFIELDS%YRGFL,YDMODEL%YRML_GCONF,JST,KDIGLST,ISUB,&
           &      YDFIELDS%YRSPEC%SP3D,YDFIELDS%YRSPEC%SP2D,YDFIELDS%YRSPEC%SP1D,LDIAB)
  ENDIF
    PSPA3(1:NFLEVL,:,:)=PSPA3(1:NFLEVL,:,:) + RACCSPA3(1:NFLEVL,:,:,1)
    IF (MYSETV  ==  NBSETSP) THEN
     PSPA2(:,:)=PSPA2(:,:) + RACCSPA2(:,:,1)
    ENDIF
    PSPA1(1:NFLEVL,:)=PSPA1(1:NFLEVL,:) + RACCSPA1(1:NFLEVL,:)
    CALL DEALLDFI(NULOUT)
  ENDIF

!    Store GFLGP values at middle time 
   IF (JSTEP == ISTOP/2) THEN
      CALL DIGP(YDGEOMETRY,YDFIELDS%YRGFL,YGFL,1,IGPP,ZGFLGP)
   ENDIF   

!        2.2  Coupling

  IF (LELAM) THEN
    IF (KDFI == 4) THEN
      NSTEP=JSTEP
    ELSEIF (KDFI == 6) THEN
      NSTEP=ISIGN*JSTEP+MAX(0,ISIGN)*IDEC
    ELSE
      NSTEP=JST
    ENDIF
    NSTAR2=0
    CALL ELSRW(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,.TRUE.)
    NSTAR2=ISTART
  ENDIF

!        2.3  Current value of the time step length

  IF (JSTEP == 0.OR. LTWOTL) THEN
    TDT=TSTEP
  ELSE
    TDT=2.0_JPRB*TSTEP
  ENDIF
  CALL UPDTIM(YDGEOMETRY,YDFIELDS%YRSURF,YDMODEL,JSTEP,TDT,TSTEP,.FALSE.)

!        2.4  Reset semi-implicit solver in the multilevel model
!             Reset horizontal diffusion

  IF (LELAM) THEN
    ! * SI scheme
    IF (LNHEE.AND.(.NOT.YRDYNA%LSI_NHEE)) THEN
      CALL SUNHSI(YDMODEL%YRCST,YDGEOMETRY,YDRIP,YDDYN,YDEDYN,NULOUT,.FALSE.)
    ELSEIF (LNHEE.AND.YRDYNA%LSI_NHEE) THEN
      CALL SUNHEESI(YDMODEL%YRCST,YDGEOMETRY,YDRIP,YDDYN,NULOUT,.FALSE.)
    ELSEIF (LNHQE.AND.(.NOT.YRDYNA%LNHQE_SIHYD)) THEN
      CALL SUNHQESI(YDMODEL%YRCST,YDGEOMETRY,YDRIP,YDDYN,NULOUT,.FALSE.)
    ELSE
      IF (LESIDG) CALL SUEHEG(YDGEOMETRY,YDDYN,YDEDYN,YDRIP)
    ENDIF
  ELSE
    ! * SI scheme
    IF (LNHEE.AND.(.NOT.YRDYNA%LSI_NHEE)) THEN
      CALL SUNHSI(YDMODEL%YRCST,YDGEOMETRY,YDRIP,YDDYN,YDEDYN,NULOUT,.FALSE.)
    ELSEIF (LNHEE.AND.YRDYNA%LSI_NHEE) THEN
      CALL SUNHEESI(YDMODEL%YRCST,YDGEOMETRY,YDRIP,YDDYN,NULOUT,.FALSE.)
    ELSEIF (LNHQE.AND.(.NOT.YRDYNA%LNHQE_SIHYD)) THEN
      CALL SUNHQESI(YDMODEL%YRCST,YDGEOMETRY,YDRIP,YDDYN,NULOUT,.FALSE.)
    ELSE
      IF (LSIDG) CALL SUHEG(YDGEOMETRY,YDRIP,YDDYN)
    ENDIF
    ! * Horizontal diffusion scheme
    IF (CL9 == 'A') THEN
      IF (LSTRHD) CALL SUHDU(YDGEOMETRY,YDMODEL%YRML_GCONF,YDDYN)
    ENDIF
  ENDIF

!        2.5  Time step

!       Reset YOMCT3
  NSTEP=JSTEP

! Write fields at half the time-span for finalisation
  IF (LLWRITE.AND.NSTEP == IDEC) THEN
    CLCONF(1:1)='A'
  ELSE
    CLCONF(1:1)=CL1
  ENDIF

! Prepare screening
  IF (LLOBSC1) THEN
    LLSLOT=.FALSE.
    CALL CHKOBTIM(YDGEOMETRY,YDRIP,LLSLOT,IACTIM,YDODB)
    IF(LLSLOT) THEN
      CLCONF(6:6)='V'
    ELSE
      CLCONF(6:6)=CL6
    ENDIF
  ENDIF

! Time-step
  IF (NSITER > 0) THEN

    ! Non-hydrostatic or RUBC or predictor-corrector time-step

    NCURRENT_ITER=0
    IF (LRHDI_LASTITERPC) THEN
      CLCONF(9:9)='I'
    ELSE
      CLCONF(9:9)=CL9
    ENDIF
    CALL STEPO(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,CLCONF,YDJOT=YDJOT)
    CLCONF(1:1)=CL1
    CLCONF(6:6)=CL6

    DO JSITER=1,NSITER-1
      NCURRENT_ITER=JSITER
      IF (LRHDI_LASTITERPC) THEN
        CLCONF(9:9)='I'
      ELSE
        CLCONF(9:9)=CL9
      ENDIF
      CALL STEPO(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,CLCONF,YDJOT=YDJOT)
    ENDDO

    NCURRENT_ITER=NSITER
    CLCONF(4:4)='A'
    CLCONF(9:9)=CL9
    CALL STEPO(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,CLCONF,YDJOT=YDJOT)
    CLCONF(4:4)=CL4

  ELSE
    ! Ordinary time step
    NCURRENT_ITER=0
    CALL STEPO(YDGEOMETRY,YDFIELDS,YDMTRAJ,YDMODEL,CLCONF,YDJOT=YDJOT,YDVARBC=YDVARBC,&
& YDGOM5=YDGOM5,YDODB=YDODB)
  ENDIF

!        2.6  Printings

  IF (KDFI == 1) THEN
    ISTEP=-NSTEP
  ELSEIF (KDFI == 2) THEN
    ISTEP=NSTEP
  ELSEIF (KDFI == 3) THEN
    ISTEP=NSTEP+KDIGLST
  ELSEIF (KDFI == 4) THEN
    ISTEP=NSTEP
  ELSEIF (KDFI == 5) THEN
    ISTEP=-1
  ELSEIF (KDFI == 6) THEN
    ISTEP=ISIGN*NSTEP
  ENDIF

  CALL USER_CLOCK(PTOTAL_CP=ZT2)
  CALL DATE_AND_TIME(TIME=CLTIMEOD)
  IF ((KDFI == 1).AND.(ISTEP == 0)) THEN
    WRITE (UNIT=CLDAYF,FMT=&
     & '(1X,A,'':'',A,'':'',A,1X,''DFI STEP   -'',I1,&
     & ''/'',I3,'' +CPU='',F6.3)')&
     & CLTIMEOD(1:2),CLTIMEOD(3:4),CLTIMEOD(5:6),&
     & ISTEP,2*KDIGLST,ZT2-ZT1  
  ELSE
    WRITE (UNIT=CLDAYF,FMT=&
     & '(1X,A,'':'',A,'':'',A,1X,''DFI STEP '',SP,I4,SS,&
     & ''/'',I3,'' +CPU='',F6.3)')&
     & CLTIMEOD(1:2),CLTIMEOD(3:4),CLTIMEOD(5:6),&
     & ISTEP,2*KDIGLST,ZT2-ZT1  
  ENDIF
  IF(MYPROC==1) WRITE(NULERR,*) CLDAYF
  WRITE (UNIT=NULOUT, FMT='(A39)') CLDAYF

ENDDO

!-----------------------------------------------------------------

!* 3.  Reset YOMCT3
!      ------------

! Last call for screening
IF (LLOBSC1) THEN
  CLCONF='0AA00V000'
  WRITE(NULOUT,'('' NSTEP ='',I6,'' OBSV    '',A9)') NSTEP,CLCONF
  CALL OBSV(YDEPHY,YDMODEL%YRML_PHY_MF,YDJOT,YDVARBC,YDTCV,YDGOM5,YDODB,MTSLOTNO,'DI')
  CALL UPDOBS(YDRIP)
ENDIF

!-----------------------------------------------------------------

!* 4.  Reset GFLGP to middle window values
!      -----------------------------------

CALL DIGP(YDGEOMETRY,YDFIELDS%YRGFL,YGFL,-1,IGPP,ZGFLGP)
DEALLOCATE(ZGFLGP)
NSTEP=0

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('DFI3',1,ZHOOK_HANDLE)
END SUBROUTINE DFI3
