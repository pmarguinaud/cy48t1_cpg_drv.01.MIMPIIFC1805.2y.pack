SUBROUTINE EINCLI10(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF)

!**** *EINCLI10*

!     PURPOSE.
!     --------

!     This routine calculates climatological fields for an aqua-planet.
!      . monthly SST is interpolated from input data, or imported, or fixed 
!        by namelist (constant)  
!      . roughness lengths, albedo and emissivity are set according to 
!        the sea/sea-ice limit.
!      . other fields are constants.

!**   INTERFACE.
!     ----------

!     CALL EINCLI10

!     METHOD.
!     -------

!     The input dataset is read on unit 11 (file sst_GL).
!     If not available, the FA file Newsst is read or 
!     a constant value is imposed from namelist (RTT+RSTR, RSTR in NAMCLI).
!     The current month is read in NAMRIP.
!     The case of several deep layers (NCSS > 1) is not handled.
!     The values on the subgrid is obtained by a 12-point 
!     interpolation operator applied on the source grid.

!     EXTERNALS.
!     ----------

!     AUTHORS.
!     --------

!      D. Giard   00-03-07   

!     MODIFICATIONS.
!     --------------
!      D. Giard   03-06-16 update  
!      D. Giard   04-09-15 adding ozone and aerosols, phasing
!      D. Giard   04-11-16 bugfix
!      D. Giard   05-04-04 GRIBEX
!      D. Giard   05-04-07 no packing for roughness lengths, new EBICLI
!     JD. Gril    05-05-12 : add LMRT=true => NROTEQ/ZSINLA(1) = -2
!     JD. Gril    06-02-01 : Test for LMRT=.T. and NCADFORM=0 => abort
!      A. Bogatchev 13-04-11 phasing cy40, coherence with modified modules
!                            and renamed modules and functions
!      B. Bochenek (Apr 2015): Phasing: move some variables.
!      O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE YOEPHY       , ONLY : TEPHY
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB     ,JPRD 
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMLUN   , ONLY : NULOUT
USE YOMVERT  , ONLY : VP00
USE YOMCST   , ONLY : RG       ,RTT
USE YOMRIP0  , ONLY : NINDAT
USE YOMCLI   , ONLY : YRCLI
USE YOMOPH0  , ONLY : NCADFORM
USE YOM_YGFL , ONLY : TYPE_GFLD

!     ------------------------------------------------------------------

IMPLICIT NONE

! JPFIX : number of fixed fields
! JPVAR : number of monthly fields
! JPFLD : number of fields
! JPBX/JPBY : number of extra longitudes/latitudes on each side of the data
!          grid required by interpolation (4 or 12 points -> 2/2)

TYPE(GEOMETRY), INTENT(INOUT)   :: YDGEOMETRY
TYPE(TYPE_GFLD)    ,INTENT(INOUT):: YDGFL
TYPE(TEPHY)    ,INTENT(INOUT)   :: YDEPHY
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
INTEGER(KIND=JPIM) :: JPBX
INTEGER(KIND=JPIM) :: JPBY
INTEGER(KIND=JPIM) :: JPFIX
INTEGER(KIND=JPIM) :: JPFLD
INTEGER(KIND=JPIM) :: JPVAR

PARAMETER(JPFIX=32,JPVAR=7,JPFLD=JPFIX+JPVAR,JPBX=2,JPBY=2)

REAL(KIND=JPRB) :: ZS(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON,JPFLD),ZSINLA(YDGEOMETRY%YRDIM%NDGNH+1),&
 & ZINT((YDGEOMETRY%YRDIM%NDGUXG-YDGEOMETRY%YRDIM%NDGUNG+1)*(YDGEOMETRY%YRDIM%NDLUXG-YDGEOMETRY%YRDIM%NDLUNG+1)),&
 & ZIMP(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON)
REAL(KIND=JPRB),ALLOCATABLE :: ZFLD(:,:),ZRNM(:)
REAL(KIND=JPRD),ALLOCATABLE :: ZAUX(:)
REAL(KIND=JPRB) :: ZTR, ZCODIL, ZOZA, ZOZB, ZOZC, ZAES, ZAEL, ZAEU, ZAED
REAL(KIND=JPRB) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: INIV0,INIV(JPFLD),INLOPA(YDGEOMETRY%YRDIM%NDGLG+2),INOZPA(YDGEOMETRY%YRDIM%NDGNH+1)
INTEGER(KIND=JPIM) :: IDATF(11),IDATI(11)
INTEGER(KIND=JPIM) :: IARI, IARP, IBCSP, IBPDG, ICO, IDATX, IDATY,&
 & IFLD, IGRIBC, IGRIBO, IINF, ILAP, IMES, IMOD,&
 & INDEX1, INDEX2, INUI, INUM, IOS, IREP, ITFING,&
 & ITRON, ITRONC, ITYPTR, IXFING, IYFING, IZ,&
 & J, JJ, JX, JY

CHARACTER :: CLPRE0*8,CLPRE(JPFLD)*8,CLSUF0*12,CLSUF(JPFLD)*12
CHARACTER :: CLNOMC*16,CLNOMI*16,CLNOMF*12,CLFORM*12

LOGICAL :: LLBIP(JPFLD),LLWRI(JPFLD),LLPAC(JPFLD)
LOGICAL :: LLCOSP, LLIMST, LLIMP, LLKEEP, LLSST, LLOPN

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "cchien.intfb.h"
#include "ebicli.intfb.h"
#include "echk923.intfb.h"
#include "egeo923.intfb.h"
#include "einter2.intfb.h"
!#include "einter8.intfb.h"

#include "fcttim.func.h"

!     ------------------------------------------------------------------
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EINCLI10',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
 & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YREGEO=>YDGEOMETRY%YREGEO, &
 & YDVAB=>YDGEOMETRY%YRVAB, &
 & YDPHY1=>YDML_PHY_MF%YRPHY1)

ASSOCIATE(TMERGL=>YDPHY1%TMERGL, &
 & NBZONL=>YDGEOMETRY%YREDIM%NBZONL, NBZONG=>YDGEOMETRY%YREDIM%NBZONG, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NSEFRE=>YDDIM%NSEFRE, NDGLG=>YDDIM%NDGLG, NSMAX=>YDDIM%NSMAX, &
 & NDLUXG=>YDDIM%NDLUXG, NDGUNG=>YDDIM%NDGUNG, NMSMAX=>YDDIM%NMSMAX, &
 & NDGUXG=>YDDIM%NDGUXG, NSPEC2G=>YDDIM%NSPEC2G, NDLUNG=>YDDIM%NDLUNG, &
 & NDGNH=>YDDIM%NDGNH, NDLON=>YDDIM%NDLON, &
 & ELON1=>YREGEO%ELON1, ELON0=>YREGEO%ELON0, ELATC=>YREGEO%ELATC, &
 & ELON2=>YREGEO%ELON2, LMAP=>YREGEO%LMAP, EXWN=>YREGEO%EXWN, LMRT=>YREGEO%LMRT, &
 & ELX=>YREGEO%ELX, ELAT1=>YREGEO%ELAT1, ELAT0=>YREGEO%ELAT0, &
 & ELONC=>YREGEO%ELONC, ELAT2=>YREGEO%ELAT2, EYWN=>YREGEO%EYWN, &
 & ERPK=>YREGEO%ERPK, EDELY=>YREGEO%EDELY, EDELX=>YREGEO%EDELX, &
 & ELY=>YREGEO%ELY )
!     1. SET INITIAL VALUES AND DEFINE THE SUBGRID.
!        ------------------------------------------

!     1.1 Reference values

ZTR=RTT+YRCLI%RSTR
ZOZA=0.066_JPRB
ZOZB=3166._JPRB
ZOZC=3.0_JPRB
ZAES=6.5E-3_JPRB
ZAEL=0.043_JPRB
ZAEU=2.8E-3_JPRB
ZAED=0.026_JPRB


!     1.2 Dimensions for data and interpolation.

ICO=YRCLI%NPINT*YRCLI%NPINT
IXFING=NDLUXG-NDLUNG+1
IYFING=NDGUXG-NDGUNG+1
ITFING=IXFING*IYFING
IDATX=YRCLI%NDATX+2*JPBX
IDATY=YRCLI%NDATY+2*JPBY


!     1.3 Type of input data files.

IF (YRCLI%LIEEE) THEN
  CLFORM='UNFORMATTED'
ELSE
  CLFORM='FORMATTED'
ENDIF


!     1.4 Model grid.

CALL EGEO923(YDGEOMETRY%YRGEM)


!     1.5 ALADIN files features

! Cadre
CLNOMC='Const.Clim.Surfa'
LLKEEP=.TRUE.
ITYPTR=-NMSMAX
ITRONC=NSMAX
ZCODIL=0.0_JPRB
IF (.NOT.LMAP) ZCODIL=-1.0_JPRB

DO J=1,NDGNH+1
  INOZPA(J)=0
  ZSINLA(J)=0.0_JPRB
ENDDO
DO J=1,NDGLG+2
  INLOPA(J) = 0
ENDDO

IF(NCADFORM == 0) THEN
  IF (LMRT) THEN
    CALL ABOR1('SUFRAME : Mercator Rot/Til and OLD CADRE not allowed !')
  ENDIF
  ZSINLA(1)  = 0.0_JPRB              ! ex-NROTEQ
  ZSINLA(2)  = 0.0_JPRB              ! ex-ELONR
  ZSINLA(3)  = 0.0_JPRB              ! ex-ELATR
  ZSINLA(4)  = ELON1
  ZSINLA(5)  = ELAT1
  ZSINLA(6)  = ELON2
  ZSINLA(7)  = ELAT2
  ZSINLA(8)  = ELON0
  ZSINLA(9)  = ELAT0
  ZSINLA(10) = ERPK
  ZSINLA(11) = 0.0_JPRB              ! ex-NSOTRP
  ZSINLA(12) = 0.0_JPRB              ! ex-NGIV0
  ZSINLA(13) = ELX
  ZSINLA(14) = ELY
  ZSINLA(15) = EDELX
  ZSINLA(16) = EDELY
  ZSINLA(17) = EXWN
  ZSINLA(18) = EYWN
ELSE
  ZSINLA(1)  = -1.0_JPRB
  IF (LMRT) ZSINLA(1)  = ZSINLA(1) - 1.0_JPRB ! Mercator Rot/Tilted = -2
  ZSINLA(2)  = ERPK
  ZSINLA(3)  = ELON0
  ZSINLA(4)  = ELAT0
  ZSINLA(5)  = ELONC
  ZSINLA(6)  = ELATC
  ZSINLA(7)  = EDELX
  ZSINLA(8)  = EDELY
  ZSINLA(9)  = ELX
  ZSINLA(10) = ELY
  ZSINLA(11) = EXWN
  ZSINLA(12) = EYWN
  ZSINLA(13) = ELON1
  ZSINLA(14) = ELAT1
  ZSINLA(15) = ELON2
  ZSINLA(16) = ELAT2
  ZSINLA(17) = 0.0_JPRB              ! free
  ZSINLA(18) = 0.0_JPRB              ! free
ENDIF

INLOPA(1) = 10
INLOPA(2) = 1
INLOPA(3) = NDLUNG
INLOPA(4) = NDLUXG
INLOPA(5) = NDGUNG
INLOPA(6) = NDGUXG
INLOPA(7) = NBZONL
INLOPA(8) = NBZONG

! Date
DO J=1,11
  IDATF(J)=0
ENDDO
IDATF(1)=1
IDATF(2)=NMM(NINDAT)
IDATF(3)=15
IDATF(6)=1

! Output file
INUM=3
IMES=1
IARP=33
IREP=0
IARI=0
CLNOMF='Const.Clim'
LLIMST=.TRUE.

! Default for GRIB coding
IGRIBO=0
IGRIBC=0
IBPDG=0
IBCSP=0
ITRON=0
ILAP=0
IMOD=0


!     1.6 Names of ARPEGE fields

CLPRE0='SPECSURF'
DO J=1,JPFLD
  CLPRE(J)='SURF'
ENDDO
CLPRE(JPFIX-1)='PROF'
CLPRE(JPFIX  )='RELA'
CLPRE(JPFLD-1)='PROF'
CLPRE(JPFLD  )='RELA'

! Constant fields over sea
! computed by EINCLI1
CLSUF0='GEOPOTEN    '
CLSUF( 1)='GEOPOTENTIEL'
CLSUF( 2)='IND.TERREMER'
CLSUF( 3)='ET.GEOPOTENT'
CLSUF( 4)='VAR.GEOP.ANI'
CLSUF( 5)='VAR.GEOP.DIR'
CLSUF( 6)='PROP.TERRE  '
CLSUF( 7)='PROP.URBANIS'
CLSUF( 8)='Z0REL.FOIS.G'
! computed by EINCLI2
CLSUF( 9)='IND.VEG.DOMI'
CLSUF(10)='ALBEDO.SOLNU'
CLSUF(11)='ALBEDO.COMPL'
CLSUF(12)='PROP.ARGILE '
CLSUF(13)='PROP.SABLE  '
CLSUF(14)='EPAI.SOL.MAX'
CLSUF(15)='EPAIS.SOL   '
CLSUF(16)='PROP.VEG.MAX'
! computed by EINCLI4
CLSUF(17)='PROP.VEGETAT'
CLSUF(18)='IND.FOLIAIRE'
CLSUF(19)='RESI.STO.MIN'
CLSUF(20)='Z0VEG.FOIS.G'
CLSUF(21)='ALBEDO.VEG'
! computed by EINCLI8
CLSUF(22)='A.OF.OZONE  '
CLSUF(23)='B.OF.OZONE  '
CLSUF(24)='C.OF.OZONE  '
! computed by EINCLI9
CLSUF(25)='AEROS.SEA   '
CLSUF(26)='AEROS.LAND  '
CLSUF(27)='AEROS.SOOT  '
CLSUF(28)='AEROS.DESERT'
! computed by EINCLI3
CLSUF(29)='RESERV.NEIGE'
CLSUF(30)='PROP.RMAX.EA'
CLSUF(31)='PROP.RMAX.EA'
CLSUF(32)='PROP.RMAX.EA'

! Varying fields over sea
CLSUF(JPFIX+1)='ALBEDO      '
CLSUF(JPFIX+2)='EMISSIVITE  '
CLSUF(JPFIX+3)='Z0.FOIS.G   '
CLSUF(JPFIX+4)='GZ0.THERM   '
CLSUF(JPFIX+5)='TEMPERATURE '
CLSUF(JPFIX+6)='TEMPERATURE '
CLSUF(JPFIX+7)='TEMPERATURE '

INIV0=0
DO J=1,JPFLD
  INIV(J)=0
ENDDO
INIV(JPFIX-3)=1
INIV(JPFIX-2)=1
INIV(JPFIX-1)=1
INIV(JPFIX  )=1
INIV(JPFLD-2)=1
INIV(JPFLD-1)=1
INIV(JPFLD  )=1

DO J=1,JPFLD
  LLWRI(J)=.TRUE.
  LLPAC(J)=.TRUE.
ENDDO
LLPAC( 8)=.FALSE.
LLPAC(20)=.FALSE.
LLPAC(JPFIX+3)=.FALSE.
LLPAC(JPFIX+4)=.FALSE.

DO J=1,JPFIX
  LLBIP(J)=.FALSE.
ENDDO
DO J=JPFIX+1,JPFLD
  LLBIP(J)=.TRUE.
ENDDO


!     ------------------------------------------------------------------

!     2. FIRST STEPS.
!        ------------

!     2.1 Open output file

! Define the CADRE
CALL FACADE(CLNOMC,ITYPTR,0.0_JPRB,0.0_JPRB,0.0_JPRB,ZCODIL,&
 & ITRONC,NDGLG,NDLON,INLOPA,INOZPA,ZSINLA,&
 & NFLEVG,VP00,YDVAB%VALH,YDVAB%VBH,LLKEEP)

! Open the file
CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'NEW',.TRUE.,LLIMST,&
 & IMES,IARP,IARI,CLNOMC)

! Write the date
CALL FANDAR(IREP,INUM,IDATF)

!  Get GRIB options
CALL FAVEUR(IREP,INUM,IGRIBC,IBPDG,IBCSP,ITRON,ILAP,IMOD)

! New spectral ordering only if packing type is -1 or 3
IF (IGRIBC == -1 .OR. IGRIBC == 3) THEN
  IGRIBO=-1
ELSE
  IGRIBO=0
ENDIF


!     2.2 Final values for fixed fields

DO J=1,NDGLG*NDLON
! computed by EINCLI1
  ZS(J, 1)= 0.0_JPRB
  ZS(J, 2)= 0.0_JPRB
  ZS(J, 3)= 0.0_JPRB
  ZS(J, 4)= 1.0_JPRB
  ZS(J, 5)= 0.0_JPRB
  ZS(J, 6)= 0.0_JPRB
  ZS(J, 7)= 0.0_JPRB
  ZS(J, 8)= RG*YRCLI%SZZ0M
! computed by EINCLI2
  ZS(J, 9)= REAL(YRCLI%NTPMER,JPRB)
  ZS(J,10)= YRCLI%SALBM 
  ZS(J,11)= YRCLI%SALBM 
  ZS(J,12)= YRCLI%SARGN 
  ZS(J,13)= YRCLI%SSABN 
  ZS(J,14)= YRCLI%SDEPX 
  ZS(J,15)= YRCLI%SDEPX
  ZS(J,16)= 0.0_JPRB
! computed by EINCLI4
  ZS(J,17)= 0.0_JPRB
  ZS(J,18)= 0.0_JPRB
  ZS(J,19)= YRCLI%SRSMX
  ZS(J,20)= RG*YRCLI%SZZ0N
  ZS(J,21)= 0.0_JPRB
! computed by EINCLI8
  ZS(J,22)= ZOZA
  ZS(J,23)= ZOZB
  ZS(J,24)= ZOZC
! computed by EINCLI9
  ZS(J,25)= ZAES
  ZS(J,26)= ZAEL
  ZS(J,27)= ZAEU
  ZS(J,28)= ZAED
! computed by EINCLI3
  ZS(J,29)= 0.0_JPRB
  ZS(J,30)= 1.0_JPRB
  ZS(J,31)= 1.0_JPRB
  ZS(J,32)= 1.0_JPRB
ENDDO


!     2.3 Default values for monthly fields

DO J=1,NDGLG*NDLON
  ZS(J,JPFIX+1)= YRCLI%SALBM
  ZS(J,JPFIX+2)= YRCLI%SEMIM
  ZS(J,JPFIX+3)= RG*YRCLI%SZZ0M
  ZS(J,JPFIX+4)= RG*YRCLI%SZZ0M
  ZS(J,JPFIX+5)= ZTR
  ZS(J,JPFIX+6)= ZTR
  ZS(J,JPFIX+7)= ZTR
ENDDO

!     ------------------------------------------------------------------

!     3. READING DATA.
!        -------------

!     3.1 Check input datasets

!  Sea Surface Temperature (interpolated)
LLSST=.FALSE.
LLOPN=.FALSE.
IOS= 0
INQUIRE(FILE='sst_GL',IOSTAT=IOS,EXIST=LLSST,OPENED=LLOPN)
LLSST= LLSST .AND. (IOS == 0) .AND. .NOT.LLOPN

!  Sea Surface Temperature (imported)
LLIMP=.FALSE.
LLOPN=.FALSE.
IOS= 0
INQUIRE(FILE='Newsst',IOSTAT=IOS,EXIST=LLIMP,OPENED=LLOPN)
LLIMP= LLIMP .AND. (IOS == 0) .AND. .NOT.LLOPN


!     3.2 Read and interpolate regular input data

IF (LLSST) THEN
  WRITE(NULOUT,'('' EINCLI10 : SST INTERPOLATED FROM FILE sst_GL'')')

! Read input data
  ALLOCATE ( ZFLD(IDATX,IDATY) )
  OPEN(UNIT=11,FILE='sst_GL',FORM=CLFORM)
  IF (YRCLI%LIEEE) THEN
    ALLOCATE ( ZAUX(YRCLI%NDATX*YRCLI%NDATY) )
    READ(11) ZAUX
    DO JJ=YRCLI%NDATY+JPBY,1+JPBY,-1
      DO J=1+JPBX,YRCLI%NDATX+JPBX
        IZ=J-JPBX+(YRCLI%NDATY+JPBY-JJ)*YRCLI%NDATX
        ZFLD(J,JJ)=ZAUX(IZ)
      ENDDO
    ENDDO
    DEALLOCATE ( ZAUX )
  ELSE
    READ(11,*)((ZFLD(J,JJ),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=YRCLI%NDATY+JPBY,1+JPBY,-1)
  ENDIF
  CLOSE(11)

! Interpolate
! 12 points
  CALL EINTER2(YDGEOMETRY%YREGEO,ZFLD,IDATY,IDATX,1,ZINT,&
   & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)
! 4 points
!  CALL EINTER8(ZFLD,IDATY,IDATX,1,ZINT,&
!  & IYFING,IXFING,ITFING,EDELX,EDELY,NPINT,ICO,NULOUT)
  DEALLOCATE ( ZFLD )


!     3.3 Import SST from a FA file

ELSEIF (LLIMP) THEN
  WRITE(NULOUT,'('' EINCLI10 : SST IMPORTED FROM FILE Newsst'')')

! File characteristics
  IREP= 0
  INUI= 10
  IMES= 1
  IARP= 1
  IARI= 0
  CLNOMI='Cadre.orog      '
  LLIMST=.TRUE.

! Open the file
  CALL FAITOU(IREP,INUI,.TRUE.,'Newsst','OLD',.TRUE.,LLIMST,&
   & IMES,IARP,IARI,CLNOMI)

! Check geometry
  IINF=-2
  CALL CCHIEN(YDGEOMETRY,CLNOMI,INUI,IINF)

! Check date
  CALL FADIES(IREP,INUI,IDATI)
  WRITE(NULOUT,'('' EINCLI10 : MONTHS OF MODEL AND FILE '',2I3.2)')&
  & IDATF(2),IDATI(2)

! Read SST
  LLCOSP=.FALSE.
  CALL FACILE(IREP,INUI,'SURF',1,'TEMPERATURE ',ZIMP,LLCOSP)
  IF (IINF == -3) THEN
    ZINT(1:ITFING)= ZIMP(1:ITFING)
  ELSE
    INDEX1= 0
    INDEX2= (NDGUNG-1)*NDLON+NDLUNG-1
    DO JY=1,NDGUXG-NDGUNG+1
      DO JX=1,NDLUXG-NDLUNG+1
        ZINT(JX+INDEX1)= ZIMP(JX+INDEX2)
      ENDDO
      INDEX1= INDEX1+NDLUXG-NDLUNG+1
      INDEX2= INDEX2+NDLON
    ENDDO
  ENDIF

! Close file
  CALL FAIRME(IREP,INUI,'KEEP')


!     3.4 Constant SST

ELSE
  WRITE(NULOUT,'('' EINCLI10 : CONSTANT SST '',F7.3)') ZTR
  DO J=JPFIX+1,JPFLD
    LLBIP(J)=.FALSE.
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!     4. CORRECT FIELDS ACCORDING TO SEA-ICE.
!        ------------------------------------

IF (LLSST .OR. LLIMP) THEN
  DO J=1,ITFING
    ZS(J,JPFIX+5)= ZINT(J)
    ZS(J,JPFIX+6)= ZINT(J)
    ZS(J,JPFIX+7)= ZINT(J)
    IF (ZINT(J) <= TMERGL) THEN
      ZS(J,JPFIX+1)= YRCLI%SALBB
      ZS(J,JPFIX+2)= YRCLI%SEMIB
      ZS(J,JPFIX+3)= RG*YRCLI%SZZ0B
      ZS(J,JPFIX+4)= RG*YRCLI%SZZ0B
    ENDIF
  ENDDO
ENDIF

!     ------------------------------------------------------------------

!     5. WRITE.
!        ------

!     5.1 Spectral orography

! Orography is not packed at all, in order to keep the 2 representations 
! of this field in spectral agreement.
! However both are 0 here, the problem is just the size of spectral arrays

! Adapt array size to GRIB coding
SELECT CASE (IGRIBO)
CASE (0)
  ALLOCATE ( ZRNM (NSEFRE) ) 
CASE (-1)
  ALLOCATE ( ZRNM (NSPEC2G) )
CASE DEFAULT
  CALL ABOR1('INCLI10: WRONG VALUE IGRIBO')
END SELECT

! Change GRIB level of coding to "none"
CALL FAGOTE(IREP,INUM,IGRIBO,IBPDG,IBCSP,ITRON,ILAP,IMOD)

! Initialize and write spectral orography
ZRNM(:)= 0.0_JPRB
LLCOSP=.TRUE.
CALL FAIENC(IREP,INUM,CLPRE0,INIV0,CLSUF0,ZRNM,LLCOSP)

!  Reset previous GRIB options
CALL FAGOTE(IREP,INUM,IGRIBC,IBPDG,IBCSP,ITRON,ILAP,IMOD)


!     5.2 Gridpoint fields

IFLD=JPFLD
CALL EBICLI(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF,IFLD,INIV,CLPRE,CLSUF,INUM,ZS,LLBIP,LLWRI,LLPAC)


!     5.3 Controls

CALL ECHK923(YDGEOMETRY%YRDIM,INUM) 


!     5.4 Close file

CALL LFILAF(IREP,INUM,.TRUE.)
CALL FAIRME(IREP,INUM,'KEEP')

!     ------------------------------------------------------------------


END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EINCLI10',1,ZHOOK_HANDLE)
END SUBROUTINE EINCLI10
