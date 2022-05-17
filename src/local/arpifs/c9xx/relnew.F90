SUBROUTINE RELNEW(YDGEOMETRY,LDLSM,LDFIT,KPDG,PZR,PZM)

!**** *RELNEW*

!     PURPOSE.
!     --------

!     Reads a new orography ,checks the *CADRE* , and makes a spectral
!     truncation. May also read a new land-sea mask, according to LDLSM.

!**   INTERFACE.
!     ----------

!     CALL RELNEW(...)

!         LDLSM : Option for reading a new land-sea mask
!         LDFIT : Option for making spectral fit of orography
!         KPDG  : Size of the field
!         PZR   : Orography
!         PZM   : Land-sea mask
!     Output in PTRSPOR : XREV, XRNM (in m)  

!     METHOD.
!     -------

!     The fields are read on unit 'Neworog'.
!     The grid point values are returned in PZR/XREV and PZM.
!     For ARPEGE, the spectral coefficients for orography are stored in XRNM.

!     EXTERNALS.
!     ----------

!     FAITOU , FACILE , FACIES , FANION , FAIRME, FALIMU
!     CHIEN  , CCHIEN 
!     REESPE , SPEREE , EREESPE , SPEREEE , SPREORD
!     LFINFO

!     AUTHORS.
!     --------
!      M. DEQUE  92-02-15

!     MODIFICATIONS.
!     --------------
!      D. Giard    : 02-12-16 Importing an orography with a reduced 
!                             spectral resolution (from R. EL Khatib)
!      M.Hamrud   01-Oct-2003 CY28 Cleaning
!      D. Paradis & R. El Khatib : 04-07-22 GRIBEX
!      D. Giard    : 05-04-04 Update and cleaning
!      F.Taillefer : 09-10-23 Read pdg orog from PGD file
!      S. Riette   : 09-11-02 Build land/sea mask from PGD file
!      K. Yessad (Jan 2010): externalisation of group EGGX in XRD/IFSAUX
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!      F.Taillefer (sept 2016): adapt to new PGD fields name (surfex V8)
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMCT0   , ONLY : NQUAD, LELAM
USE YOMCST   , ONLY : RG
USE YOMVERT  , ONLY : VP00
USE YOMCLA   , ONLY : LIPGD, NLISSP
USE PTRSPOR  , ONLY : XRNM, XREV

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KPDG
LOGICAL           ,INTENT(IN)    :: LDLSM
LOGICAL           ,INTENT(INOUT) :: LDFIT
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PZR(KPDG)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZM(KPDG)

!     ------------------------------------------------------------------

REAL(KIND=JPRB),ALLOCATABLE :: ZSINLA(:),ZAHYBR(:),ZBHYBR(:)
REAL(KIND=JPRB) :: ZRNM(YDGEOMETRY%YRDIM%NSPEC2),&
 & ZW(YDGEOMETRY%YRDIM%NDLON*(YDGEOMETRY%YRDIM%NDGLG+2)),&
 & ZW2(YDGEOMETRY%YRDIM%NDLON*(YDGEOMETRY%YRDIM%NDGLG+2))
REAL(KIND=JPRB) :: ZEPS,ZSLAPO,ZCLOPO,ZSLOPO,ZCODIL,ZPREF
REAL(KIND=JPRB) :: ZRG
REAL(KIND=JPRB) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: IMAXLEV, IMAXTRUNC, IMAXGL, IMAXLON
INTEGER(KIND=JPIM),ALLOCATABLE :: INLOPA(:),INOZPA(:),IMENG(:)
INTEGER(KIND=JPIM) :: IADL(YDGEOMETRY%YRDIM%NDGLG)
INTEGER(KIND=JPIM) :: IADR, IARI, IARP, IDEB, IINF, ILONG, IMES,&
 & IMSMAX, INBITS, INDEX1, INDEX2, INDIC, INGRIB, INIVER, INLATI,&
 & INUM, INXLON, IPOSEX, IPUILA, IREP, IREPS, ISMAX, ISTRON, ITRONC,&
 & ITYPTR, J, JJ, JX, JY, IREP1, IREP2, IPOSEX1, IPOSEX2, ILONG1, ILONG2

CHARACTER (LEN = 16) :: CLNOMC
CHARACTER (LEN = 4)  :: CLPREF
CHARACTER (LEN = 12) :: CLFIELD

LOGICAL :: LLFICP, LLNOTF, LLKEEP, LLEXIST, LLCOSP, LL_EZONE_NOTIN_PGD

!     ------------------------------------------------------------------

#include "chien.h"

#include "abor1.intfb.h"
#include "cchien.intfb.h"
#include "ereespe.intfb.h"
#include "esperee.intfb.h"
#include "reespe.intfb.h"
#include "speree.intfb.h"
#include "spreord.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RELNEW',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM,   YDVAB=>YDGEOMETRY%YRVAB&
& )
ASSOCIATE(NDGENG=>YDDIM%NDGENG, NDGLG=>YDDIM%NDGLG, NDGNH=>YDDIM%NDGNH,   NDGSAG=>YDDIM%NDGSAG, NDGUNG=>YDDIM%NDGUNG, &
& NDGUXG=>YDDIM%NDGUXG,   NDLON=>YDDIM%NDLON, NDLUNG=>YDDIM%NDLUNG, NDLUXG=>YDDIM%NDLUXG,   NMSMAX=>YDDIM%NMSMAX,     &
& NSMAX=>YDDIM%NSMAX,   NFLEVG=>YDDIMV%NFLEVG,   NHTYP=>YDGEM%NHTYP, NLOENG=>YDGEM%NLOENG, NMENG=>YDGEM%NMENG,        &
& NSTAGP=>YDGEM%NSTAGP, NSTTYP=>YDGEM%NSTTYP, NTSTAGP=>YDGEM%NTSTAGP,   RLOCEN=>YDGEM%RLOCEN, RMUCEN=>YDGEM%RMUCEN,   &
& RSTRET=>YDGEM%RSTRET)
!     ------------------------------------------------------------------

IF (LELAM) THEN
  IADL(1:NDGLG)=NTSTAGP(1:NDGLG)-1
ELSE
  IADL(1:NDGLG)=NSTAGP(1:NDGLG)-1
ENDIF  

IREP=0
INUM=7
IMES=1
IARP=2
IARI=0
IDEB=0
CLNOMC='Cadre.orog      '
ZEPS=1.E-6_JPRB
LLFICP=.FALSE.
IREP1=0
IREP2=0

!     1. OPEN THE FILE AND CHECK GEOMETRY.
!        ---------------------------------

! Open the file with the new orography
CALL FAITOU(IREP,INUM,.TRUE.,'Neworog','OLD',.TRUE.,.TRUE.,IMES,&
 & IARP,IARI,CLNOMC)  

! Check the spectral resolution if required and limit controls
! or
! Check if Ezone is in PGD
IF (NLISSP == 1 .OR. (LELAM .AND. LIPGD)) THEN
  CALL FALIMU(IMAXLEV,IMAXTRUNC,IMAXGL,IMAXLON)
  IF (LELAM) THEN
    ALLOCATE ( INLOPA (IMAXGL) )
    ALLOCATE ( INOZPA (IMAXGL) )
    ALLOCATE ( ZSINLA (IMAXGL) )
  ELSE
    ALLOCATE ( INLOPA (NDGNH) )
    ALLOCATE ( INOZPA (NDGNH) )
    ALLOCATE ( IMENG  (NDGNH) )
    ALLOCATE ( ZSINLA (NDGNH) )
  ENDIF
  ALLOCATE ( ZAHYBR (0:IMAXLEV) )
  ALLOCATE ( ZBHYBR (0:IMAXLEV) )
  CALL FACIES(CLNOMC,ITYPTR,ZSLAPO,ZCLOPO,ZSLOPO,ZCODIL,ITRONC,&
   & INLATI,INXLON,INLOPA,INOZPA,ZSINLA,INIVER,ZPREF,&
   & ZAHYBR,ZBHYBR,LLKEEP)

  ! Is E-zone in PGD ?
  IF (LELAM .AND. LIPGD) THEN
    LL_EZONE_NOTIN_PGD = INLATI == INLOPA(6) .OR. INXLON == INLOPA(4)  ! E-zone size == 0 // no E-zone
  ENDIF
  ! smoother orography
  IF (NLISSP == 1) THEN
    WRITE(NULOUT,&
         &'('' IMPORTING OROGRAPHY WITH REDUCED SPECTRAL RESOLUTION '')')
    IF (LELAM) THEN
      ISMAX = NSMAX
      IMSMAX= NMSMAX
      NSMAX = ITRONC 
      NMSMAX= -ITYPTR
      WRITE(NULOUT,'('' NSMAX NMSMAX = '',2I6,'' / '',2I6)')&
           & NSMAX,NMSMAX,ISMAX,IMSMAX
    ELSE
      ISMAX = NSMAX
      NSMAX = ITRONC 
      IMENG(:)=NMENG(1:NDGNH)
      NMENG(1:NDGNH)=INOZPA(:)    
      WRITE(NULOUT,'('' NSMAX = '',I6,'' / '',I6)') NSMAX,ISMAX
    ENDIF
  ENDIF

  DEALLOCATE ( INLOPA )
  DEALLOCATE ( INOZPA )
  DEALLOCATE ( ZSINLA )
  DEALLOCATE ( ZAHYBR )
  DEALLOCATE ( ZBHYBR )
ENDIF

! Check geometry
IF (LELAM) THEN
! ALADIN
  IINF=-2
  IF (.NOT.LIPGD) CALL CCHIEN(YDGEOMETRY,CLNOMC,INUM,IINF)
  IF (NLISSP == 1 .AND. IINF == -3) THEN
    CALL ABOR1('RELNEW : EXTENSION ZONES CANNOT DIFFER HERE ')
  ENDIF
ELSE
! ARPEGE
  IINF=-1
  CALL CHIEN(CLNOMC,NSTTYP,RMUCEN,RLOCEN,RSTRET,NSMAX,NDGLG,NDLON,&
   & NLOENG,NMENG,NHTYP,NFLEVG,VP00,YDVAB%VALH,YDVAB%VBH,&
   & NQUAD,IINF,NDGSAG,NDGENG,ZEPS,LLFICP,NULOUT)  
  IF (LLFICP) IDEB= NLOENG(1)
ENDIF

! Restore spectral resolution if required
IF (NLISSP == 1) THEN
  IF (LELAM) THEN
    NSMAX = ISMAX
    NMSMAX= IMSMAX
  ELSE
    NSMAX = ISMAX
    NMENG(1:NDGNH)=IMENG(:)
    DEALLOCATE ( IMENG )
  ENDIF
ENDIF

!     2. READ AND TRANSFORM OROGRAPHY.
!        -----------------------------

IF (LELAM) THEN
! ALADIN

! Check gridpoint orography
  CALL LFINFO(IREP,INUM,'SURFGEOPOTENTIEL',ILONG,IPOSEX)
  IF (LIPGD) THEN
    CALL LFINFO(IREP1,INUM,'S1D_ZS',ILONG1,IPOSEX1)
    CALL LFINFO(IREP2,INUM,'SFX.ZS',ILONG2,IPOSEX2)
    IREP=IREP1+IREP2
    ILONG=ILONG1+ILONG2
    IPOSEX=IPOSEX1+IPOSEX2
    IF (ILONG1>0) CLPREF='S1D_'
    IF (ILONG2>0) CLPREF='SFX.'
  ENDIF
  LLNOTF=IREP == 0.AND.ILONG == 0.AND.IPOSEX == 0
  IF (LLNOTF) THEN
!   Gridpoint orography not available : use spectral one if possible
    IF (NLISSP == 1) THEN
      CALL ABOR1('RELNEW: GRIDPOINT OROGRAPHY REQUIRED IF SMOOTHING ')
    ENDIF
!   Check spectral orography and GRIB coding
    CALL FANION(IREPS,INUM,'SPECSURF',1,'GEOPOTENTIEL',LLEXIST,LLCOSP,&
     & INGRIB,INBITS,ISTRON,IPUILA)
    IF (LLEXIST) THEN
      LDFIT=.FALSE.
!     Old spectral ordering in file
      IF (.NOT.(INGRIB == -1 .OR. INGRIB == 3)) THEN
        CALL FACILE(IREP,INUM,'SPECSURF',1,'GEOPOTENTIEL',XRNM,.TRUE.)
        XRNM(:)= XRNM(:)/RG
        CALL SPREORD(YDGEOMETRY%YRDIM,YDGEOMETRY%YREDIM,YDGEOMETRY%YRELAP,1,XRNM,ZRNM,.TRUE.)
!     New spectral ordering in file
      ELSE
        CALL FACILE(IREP,INUM,'SPECSURF',1,'GEOPOTENTIEL',ZRNM,.TRUE.)
        ZRNM(:)=ZRNM(:)/RG
        CALL SPREORD(YDGEOMETRY%YRDIM,YDGEOMETRY%YREDIM,YDGEOMETRY%YRELAP,1,XRNM,ZRNM,.FALSE.)
      ENDIF
!     Compute gridpoint orography
      CALL ESPEREE(YDGEOMETRY,1,1,ZRNM,XREV)
      INDIC= 0
      DO JY=NDGUNG,NDGUXG
        DO JX=NDLUNG,NDLUXG
          INDIC= INDIC+1
          IADR= JX+IADL(JY)
          PZR(INDIC)= XREV(IADR)
        ENDDO
      ENDDO
    ELSE
      CALL ABOR1('RELNEW: NO OROGRAPHY AVAILABLE ')
    ENDIF
  ELSE
!   Gridpoint orography available
    IF (LIPGD) THEN
      CALL FACILE(IREP,INUM,CLPREF,1,'ZS',ZW2,.FALSE.)
!    A PGD file may not contain extension zone, we may have to shift data
      ZW(:)=0._JPRB
      IF (LL_EZONE_NOTIN_PGD) THEN
        DO JY=NDGUXG,NDGUNG,-1
          DO JX=NDLUXG,NDLUNG,-1
            ZW((JY-1)*NDLON+JX)=ZW2((JY-NDGUNG)*(NDLUXG-NDLUNG+1)+JX-NDLUNG+1)
          ENDDO
        ENDDO
      ELSE  ! E-zone is present
        ZW(:) = ZW2(:)
      ENDIF
      ZRG=1._JPRB
    ELSE
      CALL FACILE(IREP,INUM,'SURF',1,'GEOPOTENTIEL',ZW,.FALSE.)
      ZRG=RG
    ENDIF
!   Reading a smoother orography : no further spectral optimization  
    IF (NLISSP == 1) THEN
      DO JY=1,NDGLG
      DO JX=1,NDLON
        IADR=JX+IADL(JY)
        XREV(IADR)=ZW((JY-1)*NDLON+JX)/ZRG
      ENDDO
      ENDDO
      LDFIT=.FALSE.
      CALL EREESPE(YDGEOMETRY,1,1,ZRNM,XREV)
      CALL SPREORD(YDGEOMETRY%YRDIM,YDGEOMETRY%YREDIM,YDGEOMETRY%YRELAP,1,XRNM,ZRNM,.FALSE.)
      CALL ESPEREE(YDGEOMETRY,1,1,ZRNM,XREV)
!   Standard case : spectral fit and computation of XREV outside
    ELSE
      LDFIT=.TRUE.
    ENDIF 
    IF (IINF == -3) THEN
      PZR(1:KPDG)= ZW(1:KPDG)/ZRG
    ELSE
      INDEX1= 0
      INDEX2= (NDGUNG-1)*NDLON+NDLUNG-1
      DO JY=1,NDGUXG-NDGUNG+1
        DO JX=1,NDLUXG-NDLUNG+1
          PZR(JX+INDEX1)= ZW(JX+INDEX2)/ZRG
        ENDDO
        INDEX1= INDEX1+NDLUXG-NDLUNG+1
        INDEX2= INDEX2+NDLON
      ENDDO
    ENDIF
  ENDIF

ELSE
! ARPEGE

! Check gridpoint orography
  CALL LFINFO(IREP,INUM,'SURFGEOPOTENTIEL',ILONG,IPOSEX)
  LLNOTF=IREP == 0.AND.ILONG == 0.AND.IPOSEX == 0
  IF (LLNOTF) THEN
!   Gridpoint orography not available : use spectral one if possible
    IF (NLISSP == 1) THEN
      CALL ABOR1('RELNEW: GRIDPOINT OROGRAPHY REQUIRED IF SMOOTHING ')
    ENDIF
!   Check spectral orography and GRIB coding
    CALL FANION(IREP,INUM,'SPECSURF',1,'GEOPOTENTIEL',LLEXIST,LLCOSP,&
     & INGRIB,INBITS,ISTRON,IPUILA)
    IF (LLEXIST) THEN
!     Old spectral ordering in file
      IF (.NOT.(INGRIB == -1 .OR. INGRIB == 3)) THEN
        CALL FACILE(IREP,INUM,'SPECSURF',1,'GEOPOTENTIEL',XRNM,.TRUE.)
        XRNM(:)= XRNM(:)/RG
        CALL SPREORD(YDGEOMETRY%YRDIM,YDGEOMETRY%YREDIM,YDGEOMETRY%YRELAP,1,XRNM,ZRNM,.TRUE.)
!     New spectral ordering in file
      ELSE
        CALL FACILE(IREP,INUM,'SPECSURF',1,'GEOPOTENTIEL',ZRNM,.TRUE.)
        ZRNM(:)=ZRNM(:)/RG
        CALL SPREORD(YDGEOMETRY%YRDIM,YDGEOMETRY%YREDIM,YDGEOMETRY%YRELAP,1,XRNM,ZRNM,.FALSE.)
      ENDIF
!     Compute gridpoint orography
      CALL SPEREE(YDGEOMETRY,1,1,ZRNM,XREV)
      INDIC=0
      DO JJ=1,NDGLG
        DO J=1,NLOENG(JJ)
          INDIC=INDIC+1
          IADR= J+IADL(JJ)
          PZR(INDIC)=XREV(IADR)
        ENDDO
      ENDDO
    ELSE
      CALL ABOR1('RELNEW: NO OROGRAPHY AVAILABLE ')
    ENDIF
  ELSE
!   Gridpoint orography available
    CALL FACILE(IREP,INUM,'SURF',1,'GEOPOTENTIEL',ZW,.FALSE.)
    PZR(1:KPDG)= ZW(1+IDEB:KPDG+IDEB)/RG
    INDIC= 0
    DO JJ=1,NDGLG
      DO J=1,NLOENG(JJ)
        INDIC= INDIC+1
        IADR= J+IADL(JJ)
        XREV(IADR)= PZR(INDIC)
      ENDDO
    ENDDO
    CALL REESPE(YDGEOMETRY,1,1,ZRNM,XREV)
    CALL SPREORD(YDGEOMETRY%YRDIM,YDGEOMETRY%YREDIM,YDGEOMETRY%YRELAP,1,XRNM,ZRNM,.FALSE.)
    CALL SPEREE(YDGEOMETRY,1,1,ZRNM,XREV)
    INDIC=0
    DO JJ=1,NDGLG
      DO J=1,NLOENG(JJ)
        INDIC= INDIC+1
        IADR= J+IADL(JJ)
        PZR(INDIC)= XREV(IADR)
      ENDDO
    ENDDO
  ENDIF
! No further spectral fit for ARPEGE
  LDFIT=.FALSE.

ENDIF

!     3. READ AND TRANSFORM LAND-SEA MASK.
!        ---------------------------------

IF (LDLSM) THEN

! ALADIN
  IF (LELAM) THEN
    IF (LIPGD) THEN
!    Field is not present when there's no sea point in domain
      CLFIELD=CLPREF//'COVER001'
      CALL LFINFO(IREP,INUM,CLFIELD,ILONG,IPOSEX)
      IF(IREP == 0.AND.ILONG == 0.AND.IPOSEX == 0) THEN
        ZW(:)=0._JPRB
      ELSE
        CALL FACILE(IREP,INUM,CLPREF,1,'COVER001',ZW,.FALSE.)
      ENDIF
!    Field is not present when there's no lake point in domain
      CLFIELD=CLPREF//'COVER002'
      CALL LFINFO(IREP,INUM,CLFIELD,ILONG,IPOSEX)
      IF(IREP == 0.AND.ILONG == 0.AND.IPOSEX == 0) THEN
        ZW2(:)=0._JPRB
      ELSE
        CALL FACILE(IREP,INUM,CLPREF,1,'COVER002',ZW2,.FALSE.)
      ENDIF
!    A PGD file may not contain extension zone, we may have to shift data
      IF (LL_EZONE_NOTIN_PGD) THEN
        DO JY=NDGUXG,NDGUNG,-1
          DO JX=NDLUXG,NDLUNG,-1
            ZW((JY-1)*NDLON+JX)=&
             & ANINT(1._JPRB-(ZW((JY-NDGUNG)*(NDLUXG-NDLUNG+1)+JX-NDLUNG+1)&
             & +ZW2((JY-NDGUNG)*(NDLUXG-NDLUNG+1)+JX-NDLUNG+1)))
!!! this formula presumes that SMASK=0.5
          ENDDO
        ENDDO
      ELSE  ! E-zone is present
        ZW(:) = ANINT(1._JPRB - (ZW(:) + ZW2(:)))
      ENDIF
    ELSE
      CALL FACILE(IREP,INUM,'SURF',1,'IND.TERREMER',ZW,.FALSE.)
    ENDIF
    IF (IINF == -3) THEN
      PZM(1:KPDG)= 1.0_JPRB-ZW(1:KPDG)
    ELSE
      INDEX1= 0
      INDEX2= (NDGUNG-1)*NDLON+NDLUNG-1
      DO JY=1,NDGUXG-NDGUNG+1
        DO JX=1,NDLUXG-NDLUNG+1
          PZM(JX+INDEX1)= 1.0_JPRB-ZW(JX+INDEX2)
        ENDDO
        INDEX1= INDEX1+NDLUXG-NDLUNG+1
        INDEX2= INDEX2+NDLON
      ENDDO
    ENDIF
! ARPEGE
  ELSE
    CALL FACILE(IREP,INUM,'SURF',1,'IND.TERREMER',ZW,.FALSE.)
    PZM(1:KPDG)= 1.0_JPRB-ZW(1+IDEB:KPDG+IDEB)
  ENDIF

ENDIF

CALL FAIRME(IREP,INUM,'KEEP')

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RELNEW',1,ZHOOK_HANDLE)
END SUBROUTINE RELNEW


