SUBROUTINE INCLI10(YDGEOMETRY,YDPHY1)

!**** *INCLI10*

!     PURPOSE.
!     --------

!     This routine calculates climatological fields for an aqua-planet.
!      . monthly SST is interpolated from input data, or imported, or fixed 
!        by namelist (constant)  
!      . roughness lengths, albedo and emissivity are set according to 
!        the sea/sea-ice limit.
!      . other fields are constants.
!     The target grid is a Gaussian grid with rotation of the pole, and 
!     "Schmidt" stretching as well as -if asked- a reduction of the 
!     number of grid points towards the poles of the final representation.

!**   INTERFACE.
!     ----------

!     CALL INCLI10

!     METHOD.
!     -------

!     The input dataset is read on unit 11 (file sst_GL).
!     If not available, the FA file Newsst is read or 
!     a constant value is imposed from namelist (RTT+RSTR, RSTR in NAMCLI).
!     The current month is read in NAMRIP.
!     The case of several deep layers (NCSS > 1) is not handled.
!     The values on the subgrid is obtained by a 12-point 
!     interpolation operator applied on the source grid.
!     The values on the Gaussian grid are the averages of the values 
!     in the boxes (a box is the part of the subgrid corresponding 
!     to a point of the Gaussian grid).

!     EXTERNALS.
!     ----------

!      GEO923
!      INTER2 / INTER8
!      CHIEN
!      CHK923
!      FA-LFI package

!     AUTHORS.
!     --------

!      D. Giard   00-03-07 from INCLI*   

!     MODIFICATIONS.
!     --------------
!      D. Giard   03-06-16 update  
!      D. Giard   04-09-15 adding ozone and aerosols, phasing
!      D. Giard   04-11-16 bugfix
!      D. Giard   05-04-04 GRIBEX
!      K. Yessad (Jan 2010): externalisation of group EGGX in XRD/IFSAUX
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMDIL   , ONLY : SLAPO, SLOPO, GLOPO, FACDI
USE YOMCT0   , ONLY : NQUAD
USE YOMVERT  , ONLY : VP00
USE YOMCST   , ONLY : RG, RTT
USE YOMPHY1  , ONLY : TPHY1
USE YOMRIP0  , ONLY : NINDAT
USE YOMCLI   , ONLY : YRCLI

!     ------------------------------------------------------------------

IMPLICIT NONE

! JPFIX : number of fixed fields
! JPVAR : number of monthly fields
! JPFLD : number of fields
! JPBX/JPBY : number of extra longitudes/latitudes on each side of the data
!          grid required by interpolation (4 or 12 points -> 2/2)

TYPE(GEOMETRY), INTENT(IN)   :: YDGEOMETRY
TYPE(TPHY1)    ,INTENT(INOUT):: YDPHY1
INTEGER(KIND=JPIM) :: JPBX
INTEGER(KIND=JPIM) :: JPBY
INTEGER(KIND=JPIM) :: JPFIX
INTEGER(KIND=JPIM) :: JPFLD
INTEGER(KIND=JPIM) :: JPVAR

PARAMETER(JPFIX=32,JPVAR=7,JPFLD=JPFIX+JPVAR,JPBX=2,JPBY=2)

REAL(KIND=JPRB) :: ZS(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON,JPFLD)
REAL(KIND=JPRB) :: ZINT(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON),ZIMP((YDGEOMETRY%YRDIM%NDGLG+2)*YDGEOMETRY%YRDIM%NDLON)
REAL(KIND=JPRB) :: ZBO(YDGEOMETRY%YRDIM%NDGLG+1),ZLONG(YDGEOMETRY%YRDIM%NDGLG),ZMU(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON),&
 & ZSLA(YDGEOMETRY%YRDIM%NDGLG*YRCLI%NPINT),ZSLO(YDGEOMETRY%YRDIM%NDLON*YRCLI%NPINT,YDGEOMETRY%YRDIM%NDGLG),&
 & ZCLO(YDGEOMETRY%YRDIM%NDLON*YRCLI%NPINT,YDGEOMETRY%YRDIM%NDGLG),&
 & ZSINLA(YDGEOMETRY%YRDIM%NDGNH)
REAL(KIND=JPRB), ALLOCATABLE :: ZFLD(:,:),ZAUX(:),ZRNM(:)
REAL(KIND=JPRB) :: ZEPS, ZTR, ZOZA, ZOZB, ZOZC, ZAES, ZAEL, ZAEU, ZAED
REAL(KIND=JPRB) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: INIV(0:JPFLD),INLOPA(YDGEOMETRY%YRDIM%NDGNH),INOZPA(YDGEOMETRY%YRDIM%NDGNH)
INTEGER(KIND=JPIM) :: IDATF(11),IDATI(11)
INTEGER(KIND=JPIM) :: IARI, IARP, IBCSP, IBPDG, ICO, IDATX, IDATY, IDEB,&
 & IGRIBC, IGRIBO, IINF, ILAP, ILENE, ILENT, IMES, IMOD,&
 & INIQ, INJQ, INUI, INUM, IOS, IREP, ITRON, IZ, J, JJ

CHARACTER :: CLPRE(0:JPFLD)*8,CLSUF(0:JPFLD)*12
CHARACTER :: CLNOMC*16,CLNOMI*16,CLNOMF*12,CLFORM*12

LOGICAL :: LLPACK(0:JPFLD)
LOGICAL :: LLCOSP, LLIMST, LLIMP, LLKEEP, LLSST, LLOPN, LLPOLE

!     ------------------------------------------------------------------

#include "chien.h"

#include "chk923.intfb.h"
#include "geo923.intfb.h"
#include "inter2.intfb.h"
!#include "inter8.intfb.h"
#include "abor1.intfb.h"

#include "fcttim.func.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INCLI10',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
& YDVAB=>YDGEOMETRY%YRVAB)
ASSOCIATE(TMERGL=>YDPHY1%TMERGL,   NDGENG=>YDDIM%NDGENG, NDGLG=>YDDIM%NDGLG, NDGNH=>YDDIM%NDGNH,   NDGSAG=>YDDIM%NDGSAG, &
& NDLON=>YDDIM%NDLON, NSEFRE=>YDDIM%NSEFRE,   NSMAX=>YDDIM%NSMAX, NSPEC2G=>YDDIM%NSPEC2G,   NFLEVG=>YDDIMV%NFLEVG,       &
& NHTYP=>YDGEM%NHTYP, NLOENG=>YDGEM%NLOENG, NMENG=>YDGEM%NMENG,   NSTTYP=>YDGEM%NSTTYP, RLOCEN=>YDGEM%RLOCEN,            &
& RMUCEN=>YDGEM%RMUCEN,   RSTRET=>YDGEM%RSTRET)
!     ------------------------------------------------------------------

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
ZEPS=1.E-6_JPRB

!     1.2 Dimensions for data and interpolation.

ICO=YRCLI%NPINT*YRCLI%NPINT
INJQ=NDGLG*YRCLI%NPINT
INIQ=NDLON*YRCLI%NPINT
ILENT=NDGLG*NDLON
IDATX=YRCLI%NDATX+2*JPBX
IDATY=YRCLI%NDATY+2*JPBY

!     1.3 Type of input data files.

IF (YRCLI%LIEEE) THEN
  CLFORM='UNFORMATTED'
ELSE
  CLFORM='FORMATTED'
ENDIF

!     1.4 Model grid.

ILENE=0
CALL GEO923(YDGEOMETRY,YRCLI%NPINT,ILENE,ZBO,ZLONG,ZMU,ZSLA,ZSLO,ZCLO)

!     1.5 ARPEGE files features

! Cadre
CLNOMC='Const.Clim.Surfa'
LLKEEP=.TRUE.
DO J=1,NDGNH
  INLOPA(J)=NLOENG(J)
  INOZPA(J)=NMENG(J)
  ZSINLA(J)=YDCSGLEG%RMU(J)
ENDDO

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

CLPRE( 0)='SPECSURF'
DO J=1,JPFLD
  CLPRE(J)='SURF'
ENDDO
CLPRE(JPFIX-1)='PROF'
CLPRE(JPFIX  )='RELA'
CLPRE(JPFLD-1)='PROF'
CLPRE(JPFLD  )='RELA'

! Constant fields over sea
! computed by INCLI1
CLSUF( 0)='GEOPOTEN    '
CLSUF( 1)='GEOPOTENTIEL'
CLSUF( 2)='IND.TERREMER'
CLSUF( 3)='ET.GEOPOTENT'
CLSUF( 4)='VAR.GEOP.ANI'
CLSUF( 5)='VAR.GEOP.DIR'
CLSUF( 6)='PROP.TERRE  '
CLSUF( 7)='PROP.URBANIS'
CLSUF( 8)='Z0REL.FOIS.G'
! computed by INCLI2
CLSUF( 9)='IND.VEG.DOMI'
CLSUF(10)='ALBEDO.SOLNU'
CLSUF(11)='ALBEDO.COMPL'
CLSUF(12)='PROP.ARGILE '
CLSUF(13)='PROP.SABLE  '
CLSUF(14)='EPAI.SOL.MAX'
CLSUF(15)='EPAIS.SOL   '
CLSUF(16)='PROP.VEG.MAX'
! computed by INCLI4
CLSUF(17)='PROP.VEGETAT'
CLSUF(18)='IND.FOLIAIRE'
CLSUF(19)='RESI.STO.MIN'
CLSUF(20)='Z0VEG.FOIS.G'
CLSUF(21)='ALBEDO.VEG'
! computed by INCLI8
CLSUF(22)='A.OF.OZONE  '
CLSUF(23)='B.OF.OZONE  '
CLSUF(24)='C.OF.OZONE  '
! computed by INCLI9
CLSUF(25)='AEROS.SEA   '
CLSUF(26)='AEROS.LAND  '
CLSUF(27)='AEROS.SOOT  '
CLSUF(28)='AEROS.DESERT'
! computed by INCLI3
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

DO J=0,JPFLD
  INIV(J)=0
ENDDO
INIV(JPFIX-3)=1
INIV(JPFIX-2)=1
INIV(JPFIX-1)=1
INIV(JPFIX  )=1
INIV(JPFLD-2)=1
INIV(JPFLD-1)=1
INIV(JPFLD  )=1

! No packing for orography and roughness lengths
LLPACK(:)=.TRUE.
LLPACK( 0)=.FALSE.
LLPACK( 1)=.FALSE.
LLPACK( 8)=.FALSE.
LLPACK(20)=.FALSE.
LLPACK(JPFIX+3)=.FALSE.
LLPACK(JPFIX+4)=.FALSE.

!     ------------------------------------------------------------------

!     2. FIRST STEPS.
!        ------------

!     2.1 Open output file

! Define the CADRE
CALL FACADE(CLNOMC,NSTTYP,SLAPO,GLOPO,SLOPO,FACDI,&
 & NSMAX,NDGLG,NDLON,INLOPA,INOZPA,ZSINLA,&
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

DO J=1,ILENE
! computed by INCLI1
  ZS(J, 1)= 0.0_JPRB
  ZS(J, 2)= 0.0_JPRB
  ZS(J, 3)= 0.0_JPRB
  ZS(J, 4)= 1.0_JPRB
  ZS(J, 5)= 0.0_JPRB
  ZS(J, 6)= 0.0_JPRB
  ZS(J, 7)= 0.0_JPRB
  ZS(J, 8)= RG*YRCLI%SZZ0M
! computed by INCLI2
  ZS(J, 9)= REAL(YRCLI%NTPMER,JPRB)
  ZS(J,10)= YRCLI%SALBM 
  ZS(J,11)= YRCLI%SALBM 
  ZS(J,12)= YRCLI%SARGN 
  ZS(J,13)= YRCLI%SSABN 
  ZS(J,14)= YRCLI%SDEPX 
  ZS(J,15)= YRCLI%SDEPX
  ZS(J,16)= 0.0_JPRB
! computed by INCLI4
  ZS(J,17)= 0.0_JPRB
  ZS(J,18)= 0.0_JPRB
  ZS(J,19)= YRCLI%SRSMX
  ZS(J,20)= RG*YRCLI%SZZ0N
  ZS(J,21)= 0.0_JPRB
! computed by INCLI8
  ZS(J,22)= ZOZA
  ZS(J,23)= ZOZB
  ZS(J,24)= ZOZC
! computed by INCLI9
  ZS(J,25)= ZAES
  ZS(J,26)= ZAEL
  ZS(J,27)= ZAEU
  ZS(J,28)= ZAED
! computed by INCLI3
  ZS(J,29)= 0.0_JPRB
  ZS(J,30)= 1.0_JPRB
  ZS(J,31)= 1.0_JPRB
  ZS(J,32)= 1.0_JPRB
ENDDO

!     2.3 Default values for monthly fields

DO J=1,ILENE
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
  WRITE(NULOUT,'('' INCLI10 : SST INTERPOLATED FROM FILE sst_GL'')')

! Read input data
  ALLOCATE ( ZFLD(IDATX,IDATY) )
  OPEN(UNIT=11,FILE='sst_GL',FORM=CLFORM)
  IF (YRCLI%LIEEE) THEN
    ALLOCATE ( ZAUX(YRCLI%NDATX*YRCLI%NDATY) )
    READ(11) ZAUX
    DO JJ=1+JPBY,YRCLI%NDATY+JPBY
      DO J=1+JPBX,YRCLI%NDATX+JPBX
        IZ=J-JPBX+(JJ-1-JPBY)*YRCLI%NDATX
        ZFLD(J,JJ)=ZAUX(IZ)
      ENDDO
    ENDDO
    DEALLOCATE ( ZAUX )
  ELSE
    READ(11,*)((ZFLD(J,JJ),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=1+JPBY,YRCLI%NDATY+JPBY)
  ENDIF
  CLOSE(11)

! Interpolate
! 12 points
  CALL INTER2(NDGLG,NDLON,IDATY,IDATX,1,YRCLI%NPINT,ILENE,ICO,INJQ,INIQ,&
   & NLOENG(1:),ILENT,ZINT,ZFLD,ZSLA,ZSLO,ZCLO)
! 4 points
!  CALL INTER8(NDGLG,NDLON,IDATY,IDATX,1,NPINT,ILENE,ICO,INJQ,INIQ,&
!  & NLOENG(1:),ILENT,ZINT,ZFLD,ZSLA,ZSLO,ZCLO)
  DEALLOCATE ( ZFLD )

!     3.3 Import SST from FA file

ELSEIF (LLIMP) THEN
  WRITE(NULOUT,'('' INCLI10 : SST IMPORTED FROM FILE Newsst'')')

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
  IINF=-1
  LLPOLE=.FALSE.
  CALL CHIEN(CLNOMI,NSTTYP,RMUCEN,RLOCEN,RSTRET,NSMAX,NDGLG,NDLON,&
   & NLOENG,NMENG,NHTYP,NFLEVG,VP00,YDVAB%VALH,YDVAB%VBH,&
   & NQUAD,IINF,NDGSAG,NDGENG,ZEPS,LLPOLE,NULOUT)
  IF (LLPOLE) IDEB= NLOENG(1)

! Check date
  CALL FADIES(IREP,INUI,IDATI)
  WRITE(NULOUT,'('' INCLI10 : MONTHS OF MODEL AND FILE :'',2I3.2)')&
  & IDATF(2),IDATI(2)

! Read SST
  LLCOSP=.FALSE.
  CALL FACILE(IREP,INUI,'SURF',1,'TEMPERATURE ',ZIMP,LLCOSP)
  ZINT(1:ILENE)=ZIMP(1+IDEB:ILENE+IDEB)

! Close file
  CALL FAIRME(IREP,INUI,'KEEP')

!     3.4 Constant SST

ELSE
  WRITE(NULOUT,'('' INCLI10 : CONSTANT SST '',F7.3)') ZTR
ENDIF

!     ------------------------------------------------------------------

!     4. CORRECT FIELDS ACCORDING TO SEA-ICE.
!        ------------------------------------

IF (LLSST .OR. LLIMP) THEN
  DO J=1,ILENE
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

LLCOSP=.TRUE.

! Orography is not packed at all, in order to keep the 2 representations 
! of this field in spectral agreement. Change GRIB level of coding set "none"
CALL FAGOTE(IREP,INUM,IGRIBO,IBPDG,IBCSP,ITRON,ILAP,IMOD)

! Adapt array size to GRIB coding
SELECT CASE (IGRIBO)
CASE (0)
  ALLOCATE ( ZRNM (NSEFRE) ) 
CASE (-1)
  ALLOCATE ( ZRNM (NSPEC2G) )
CASE DEFAULT
  CALL ABOR1('INCLI10: WRONG VALUE IGRIBO')
END SELECT

! Initialize and write
ZRNM(:)=0.0_JPRB
CALL FAIENC(IREP,INUM,CLPRE(0),INIV,CLSUF(0),ZRNM,LLCOSP)

! Restore defaults
DEALLOCATE ( ZRNM )
LLCOSP=.FALSE.
CALL FAGOTE(IREP,INUM,IGRIBC,IBPDG,IBCSP,ITRON,ILAP,IMOD)

!     5.2 Gridpoint fields

LLCOSP=.FALSE.
DO J=1,JPFLD
  IF (.NOT.LLPACK(J)) THEN
!   Do not pack roughness lengths
    CALL FAGOTE(IREP,INUM,IGRIBO,IBPDG,IBCSP,ITRON,ILAP,IMOD)
  ENDIF

  CALL FAIENC(IREP,INUM,CLPRE(J),INIV(J),CLSUF(J),ZS(1,J),LLCOSP)

  IF (.NOT.LLPACK(J)) THEN
!   Reset packing after the treatment of roughness lengths
    CALL FAGOTE(IREP,INUM,IGRIBC,IBPDG,IBCSP,ITRON,ILAP,IMOD)
  ENDIF
ENDDO

!     5.3 Controls

CALL CHK923(YDGEOMETRY,INUM) 

!     5.4 Close file

CALL LFILAF(IREP,INUM,.TRUE.)
CALL FAIRME(IREP,INUM,'KEEP')

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INCLI10',1,ZHOOK_HANDLE)
END SUBROUTINE INCLI10
