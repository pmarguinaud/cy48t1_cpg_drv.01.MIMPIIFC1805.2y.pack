SUBROUTINE EINCLI3(YDGEOMETRY,YDGFL,YDSURF,YDEPHY,YDML_PHY_MF)

!**** *EINCLI3*

!     PURPOSE.
!     --------

!     This routine calculates for each of the 12 months five fields :
!    -climatological surface temperature
!    -climatological deep soil temperature -the "depth" being a model parameter
!                   to be accounted for in the parameterization of the routine-
!    -climatological surface relative moisture content
!    -climatological deep soil relative moisture content -same remark as above-
!    -an empirically estimated equivalent water depth snow climatology
!     It also calculates the relaxation values for deep temperature and soil
!     moisture in the same way as the deep soil values.
!     Albedo, emissivity and roughness over sea are modified according to the 
!     extension of sea-ice, and snow cover is modified according to land use 
!     type (ice cap or not), provided that the corresponding fields are present 
!     in the initial clim file.

!**   INTERFACE.
!     ----------

!     CALL EINCLI3

!     METHOD.
!     -------

!     The  source dataset is read on unit 2.
!     Input data are coming from :
!    -an old NCAR climatology for surface temperature and moisture
!    -AMIP1 mean data for sea surface temperature and sea-ice extent
!     ALL input fields are organized from 0 to 360(E) deg in longitude
!                                    from North to South in latitude

!     The target grid is a regular grid (geographical or in plane projection).
!     A subgrid of the final grid is introduced for interpolation
!     The values on the subgrid is obtained by a 12-point interpolation operator
!     applied on the source grid.
!     The values on the final grid are the averages of the values in the
!     boxes (a box is the part of the subgrid corresponding to a point of the
!     final grid).
!     Empirical formulae calculate the soil moisture and the snow amount from
!     temperature, moisture, and ice cap fraction.
!     The deep soil values are obtained by filtering the surface values with a
!     factor EXP(-(1+i)H/SQRT(T)) in Fourier space, where T is the period (in
!     days) of the wave and H the depth between the mid-surface layer and the
!     mid-nth layer scaled by the e-folding depth of the diurnal cycle (T=1 and
!     H=1 yield EXP(-(1+i)) )

!     The possible improvements to be developed are:
!     *The use of a database for monthly snow cover (or, best, snow amount)
!     *The use of a database for monthly deep soil moisture
!     *The use of a better database for soil temperature
!     *The surface temperature interpolation scheme with a 3-bin mask
!       (land / sea / sea-ice)

!     EXTERNALS.
!     ----------

!     CCHIEN
!     EGEO923
!     EINTER6
!     EINTER8
!     INCLAG
!     EBICLI (previously EBIEN)
!     ELECI
!     FA-LFI package (FAITOU,LFILAF,FAIRME)

!     AUTHORS.
!     --------
!      L. Gerard 30/05/97 from INCLI3, EINCLIA

!     MODIFICATIONS.
!     --------------
!     R. El Khatib : 01-12-06 Cleaning sm variables
!     S. Alexandru & N. Pristov 00-06-01 : use of the allocatable arrays
!     D. Giard 00-03-10 : roughness length of sea-ice;  preparing the case when 
!                         the initial albedo does not fit the minimum sea-ice 
!                         extension; small corrections for fields headers
!     M.Hamrud 03-10-01 : CY28 Cleaning
!     D. Giard 04-09-15 : cosmetic changes after automatic ones 
!     D. Giard 05-04-07 : no packing for roughness lengths, new EBICLI
!     A.Bogatchev 05-10-13 : replace YOMDPHY,NCSS, with SURFACE_FIELDS_MIX, YSP_SBD 
!     A.Bogatchev 13-04-11 : phasing cy40, coherence with modified modules 
!                            and renamed namelists and functions
!     O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!------------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE YOEPHY       , ONLY : TEPHY
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1  ,ONLY : JPIM, JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCLI   , ONLY : YRCLI
USE YOMCST   , ONLY : RPI      ,RG       ,RV       ,RCPV     ,&
 & RETV     ,RCW      ,RCS      ,RLVTT    ,RLSTT    ,&
 & RTT      ,RALPW    ,RBETW    ,RGAMW    ,RALPS    ,&
 & RBETS    ,RGAMS    ,RALPD    ,RBETD    ,RGAMD
USE YOMLUN   , ONLY : NULOUT
USE YOMSTA   , ONLY : RDTDZ1
USE YOM_YGFL , ONLY : TYPE_GFLD

!------------------------------------------------------------------------

IMPLICIT NONE

! JPNFIX : number of fixed fields in dataset
! JPNVAR : number of monthly fields in dataset
! JPNFLD : number of fields in dataset
! JPBX/JPBY : number of extra longitudes/latitudes on each side of the data
!          grid required by interpolation (4 or 12 points -> 2/2)

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TYPE_GFLD)    ,INTENT(INOUT):: YDGFL
TYPE(TSURF), INTENT(INOUT) :: YDSURF
TYPE(TEPHY)    ,INTENT(INOUT) :: YDEPHY
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
INTEGER(KIND=JPIM) :: JPBX
INTEGER(KIND=JPIM) :: JPBY
INTEGER(KIND=JPIM) :: JPNC
INTEGER(KIND=JPIM) :: JPNFIX
INTEGER(KIND=JPIM) :: JPNFLD
INTEGER(KIND=JPIM) :: JPNVAR

REAL(KIND=JPRB) :: PPRELA
REAL(KIND=JPRB) :: PPTSHT
REAL(KIND=JPRB) :: PPWSHT

PARAMETER(JPNFIX=3,JPNVAR=3,JPNFLD=JPNFIX+12*JPNVAR,JPBX=2,JPBY=2)
PARAMETER(JPNC=9,PPTSHT=0.0786_JPRB,PPWSHT=0.1969_JPRB,PPRELA=4._JPRB)

! Input grid geometry:
REAL(KIND=JPRB) ,ALLOCATABLE :: ZFLD(:,:,:)
REAL(KIND=JPRD) ,ALLOCATABLE :: ZRBUF(:,:)
! Final grid geometry:
REAL(KIND=JPRB) ,ALLOCATABLE :: ZALB(:), ZEMI(:), ZGZ0(:),&
 & ZGEO(:), ZLSM(:), ZITP(:)  
REAL(KIND=JPRB) ,ALLOCATABLE :: ZCMP(:,:), ZRES(:,:)
! Writing
REAL(KIND=JPRB) :: ZEXT(YDGEOMETRY%YRDIM%NDLON*YDGEOMETRY%YRDIM%NDGLG,0:YDSURF%YSP_SBD%NLEVS+JPNC)
REAL(KIND=JPRB) :: ZTSHT(YDSURF%YSP_SBD%NLEVS+1)

REAL(KIND=JPRB) :: FOHU
REAL(KIND=JPRB) :: FOSN
REAL(KIND=JPRB) :: FOWR

REAL(KIND=JPRB) :: PHUARG, PTDARG, PXARG, PYARG, ZEPS, ZGST, ZHU,&
 & ZNEI1, ZNEI2, ZNEI3, ZNEI4, ZNEI5, ZTD, ZTS, ZWR, ZWSHT, ZXSN, ZYSN
REAL(KIND=JPRB) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: IDATEF(11),INIVL(0:YDSURF%YSP_SBD%NLEVS+JPNC)
INTEGER(KIND=JPIM) :: IARI, IARP, IBIT, ICO, ICS, IDATX, IDATY,&
 & IFL, IFLD, IGRB, IINF, ILAP, ILEN2, IMD, IMES, IMF, INFL,&
 & INIV, INUM, IPNF, IREP, ISN, ITFING, ITRN, IWR, IXFING, IYFING,&
 & J, JC, JCS, JF, JJ, JM  

CHARACTER :: CLPREF(0:YDSURF%YSP_SBD%NLEVS+JPNC)*8,&
           & CLSUFF(0:YDSURF%YSP_SBD%NLEVS+JPNC)*12
CHARACTER :: CLNOMC*16,CLNOMF*16,CLFORM*12

LOGICAL :: LLBIP(0:YDSURF%YSP_SBD%NLEVS+JPNC),LLWRI(0:YDSURF%YSP_SBD%NLEVS+JPNC),&
         & LLPAC(0:YDSURF%YSP_SBD%NLEVS+JPNC)
LOGICAL :: LLCOSP, LLIMST, LLNEW, LLTYP

!------------------------------------------------------------------------

#include "fcttrm.func.h"

! Relative humidity function (from Td and T).
FOHU(PTDARG,PTARG)=ES(PTDARG)/ES(PTARG)

! Soil relative water content function (from surface relative humidity).
FOWR(PHUARG)=ACOS(1.0_JPRB-2.0_JPRB*PHUARG)/RPI

! Empirical snow climatology function.
FOSN(PXARG,PYARG)=MAX(0.0_JPRB,(ZNEI2*PXARG)/(ZNEI3-PXARG)&
 & +ZNEI4*(1.0_JPRB-COS(RPI*PYARG))/2.0_JPRB)  

!------------------------------------------------------------------------

#include "cchien.intfb.h"
#include "ebicli.intfb.h"
#include "egeo923.intfb.h"
#include "einter6.intfb.h"
#include "einter8.intfb.h"
#include "eleci.intfb.h"
#include "inclag.intfb.h"

!------------------------------------------------------------------------
!------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EINCLI3',0,ZHOOK_HANDLE)
ASSOCIATE(SODELX=>YDML_PHY_MF%YRPHY1%SODELX, TMERGL=>YDML_PHY_MF%YRPHY1%TMERGL, &
 & EDELY=>YDGEOMETRY%YREGEO%EDELY, EDELX=>YDGEOMETRY%YREGEO%EDELX, &
 & YSP_SBD=>YDSURF%YSP_SBD, &
 & NDGLG=>YDGEOMETRY%YRDIM%NDGLG, NDLUXG=>YDGEOMETRY%YRDIM%NDLUXG, NDGUNG=>YDGEOMETRY%YRDIM%NDGUNG, &
 & NDGUXG=>YDGEOMETRY%YRDIM%NDGUXG, NDLUNG=>YDGEOMETRY%YRDIM%NDLUNG, NDLON=>YDGEOMETRY%YRDIM%NDLON, &
 & LSOLV=>YDML_PHY_MF%YRPHY%LSOLV)
!     1. SET INITIAL VALUES.
!        -------------------

!     1.1 Number and position of fields

!  IPNF= Number of ARPEGE fields to be derived
!  ICS = Number of soil layers for temperature
!  IWR = Position of moisture fields
!  ISN = Position of snow fields
LLNEW= YSP_SBD%NLEVS > 1.OR. LSOLV
IPNF=JPNC+YSP_SBD%NLEVS
IF (LLNEW) THEN
!  ISBA or several layers in soil
!  Relaxation fields are usually useless -> copy of the first layer in soil
!  Standard version of ISBA : YSP_SBD%NLEVS=1
!   1 -> ICS     : surface and deep soil temperature
!   ICS+1, ICS+2 : surface and deep soil moisture
!   ICS+3        : snow
!   ICS+4, ICS+5 : relaxation temperature and moisture (usually unused)
! ICS=YSP_SBD%NLEVS+1 IWR=ICS+1 ISN=IWR+2
  ICS=YSP_SBD%NLEVS+1
  IWR=ICS+1
  ISN=ICS+3
  DO J=1,ICS-1
    ZTSHT(J)=(SODELX(J-1)+SODELX(J))/(2.0_JPRB*SQRT(365._JPRB))
  ENDDO
ELSE
!  Old surface scheme : YSP_SBD%NLEVS=1
!   1 -> 3 : surface, deep, relaxation soil temperature
!   4 -> 7 : surface, deep, relaxation soil moisture, snow
!  ICS=2  IWR=4=ICS+2  ISN=7=IWR+3
  ICS=3
  IWR=4
  ISN=7
  ZTSHT(1)=PPTSHT
  ZTSHT(2)=PPRELA*PPTSHT
ENDIF

!     1.2 Miscellaneous

!  Standard vertical gradient of temperature.
ZGST=RDTDZ1/RG
!  1./ZNEI1 indicates the temperature deficit time soil relative wetness
!  to saturate the quantity of non-permanent snow.
ZNEI1=0.04_JPRB
!  Relative snow equivalent for a small partial permanent ice cover.
ZNEI2=10._JPRB
!  ZNEI2/(ZNEI3-1) gives the maximum snow equivalent for total permanent
!  ice cover.
ZNEI3=1.0001_JPRB
!  Maximum non-permanent snow value.
ZNEI4=1000._JPRB
!  Permanent snow value on ice cap.
ZNEI5=10000._JPRB

!     1.3 Interpolation

ICO=YRCLI%NPINT*YRCLI%NPINT
!  Initial grid:
IDATY=YRCLI%NDATY+2*JPBY
IDATX=YRCLI%NDATX+2*JPBX
!  Final grid:
IXFING=NDLUXG-NDLUNG+1
IYFING=NDGUXG-NDGUNG+1
ITFING=IXFING*IYFING

!     1.4 Initial ALADIN clim and data files

INUM=3
IMES=1
IARP=20
IREP=0
IARI=0
LLIMST=.TRUE.
CLNOMC='Const.Clim.Surfa'

ZEPS=1.E-10_JPRB
IINF=-1
INIV=0
LLCOSP=.FALSE.
IGRB=0
IBIT=0
ITRN=0
ILAP=0

DO J=1,11
  IDATEF(J)=0
ENDDO
IDATEF(1)=1
IDATEF(2)=1
IDATEF(3)=15
IDATEF(6)=1

IF (YRCLI%LIEEE) THEN
  CLFORM='UNFORMATTED'
ELSE
  CLFORM='FORMATTED'
ENDIF

!     1.5 Geometry : initialize YEMCLI and YOMDIL

CALL EGEO923(YDGEOMETRY%YRGEM)

!     1.6 Final biperiodization and writing

DO J=0,YSP_SBD%NLEVS+JPNC
  LLBIP(J)=.TRUE.
  LLWRI(J)=.TRUE.
  LLPAC(J)=.TRUE.
  INIVL(J)=1
ENDDO
LLBIP(0)=.FALSE.
LLWRI(0)=.FALSE.
INIVL(0)=0

DO JC=0,YSP_SBD%NLEVS+JPNC
  DO J=1,NDGLG*NDLON
    ZEXT(J,JC)=0.0_JPRB
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!     2. READING  DATA.
!        --------------

!     2.1 Source data

!  Fixed fields:
!   Fraction of land , Orography*g , Fraction of the land covered by ice caps
!  Monthly fields:
!   Temperature and Td at mean orography level (K)

!  Allocate array for initial data
ALLOCATE ( ZFLD(IDATX,IDATY,JPNFLD) )
!  Allocate array for reading input file
ALLOCATE ( ZRBUF(IDATX,IDATY) )

OPEN(UNIT=10,FILE='N108_GL',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  DO JF=1,JPNFLD
    DO JJ=YRCLI%NDATY+JPBY,1+JPBY,-1
      READ(10) (ZRBUF(J,JJ),J=1+JPBX,YRCLI%NDATX+JPBX)
    ENDDO
    ZFLD(1+JPBX:YRCLI%NDATX+JPBX,1+JPBY:YRCLI%NDATY+JPBY,JF) = &
         & ZRBUF(1+JPBX:YRCLI%NDATX+JPBX,1+JPBY:YRCLI%NDATY+JPBY)
  ENDDO
ELSE
  DO JF=1,JPNFLD
    READ(10,*)((ZFLD(J,JJ,JF),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=YRCLI%NDATY+JPBY,1+JPBY,-1)
  ENDDO
ENDIF
CLOSE(10)

!  T and Td are reduced to sea level before interpolation
DO JF=JPNFIX+1,JPNFLD
  DO JJ=JPBY+1,JPBY+YRCLI%NDATY
    DO J=JPBX+1,JPBX+YRCLI%NDATX
      ZFLD(J,JJ,JF)= ZFLD(J,JJ,JF) - ZGST*ZFLD(J,JJ,2)
    ENDDO
  ENDDO
ENDDO

!     2.2 Initial ALADIN clim file

!  Open file and check CADRE compatibility
CLNOMF='Const.Clim'
CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'OLD',.TRUE.,LLIMST,&
 & IMES,IARP,IARI,CLNOMC)  
CALL CCHIEN(YDGEOMETRY,CLNOMC,INUM,IINF)

!  Allocate arrays for initial fields
ALLOCATE ( ZALB(ITFING) )
ALLOCATE ( ZEMI(ITFING) )
ALLOCATE ( ZGZ0(ITFING) )
ALLOCATE ( ZGEO(ITFING) )
ALLOCATE ( ZLSM(ITFING) )
ALLOCATE ( ZITP(ITFING) )

!  Read land-sea mask and height
INIV=0
CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'SURF',INIV,'IND.TERREMER',ZLSM,NULOUT)
CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'SURF',INIV,'GEOPOTENTIEL',ZGEO,NULOUT)

!  Put land-sea mask on C+I domain for EBICLI
ZEXT(1:ITFING,0)=ZLSM(1:ITFING)

!  Read land-use type, albedo, roughness and emissivity if available
LLTYP=.TRUE.
CALL FANION(IREP,INUM,'SURF',INIV,'IND.VEG.DOMI',LLTYP,LLCOSP,&
 & IGRB,IBIT,ITRN,ILAP)  
IF (LLTYP) THEN
  CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'SURF',INIV,'IND.VEG.DOMI',ZITP,NULOUT)
  CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'SURF',INIV,'ALBEDO      ',ZALB,NULOUT)
  CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'SURF',INIV,'EMISSIVITE  ',ZEMI,NULOUT)
  CALL ELECI(YDGEOMETRY%YRDIM,IREP,INUM,'SURF',INIV,'Z0.FOIS.G   ',ZGZ0,NULOUT)
  WRITE(NULOUT,'('' EINCLI3 :'',&
   & '' MODIFICATION OF SNOW COVER ACCORDING TO ICE CAP,'',&
   & '' OF ALBEDO, EMISSIVITY, ROUGHNESS ACCORDING TO SEA-ICE'')')  
ELSE
  WRITE(NULOUT,'('' EINCLI3 : LAND USE TYPE NOT AVAILABLE'')')
ENDIF

!  Close file
CALL FAIRME(IREP,INUM,'KEEP')

!     ------------------------------------------------------------------

!     3. INTERPOLATION.
!        --------------

!  ZCMP(J,1) interpolated land sea mask from the source grid
!    (may be different from the ARPEGE land/sea mask)
!  ZCMP(J,2) unused, and not computed (could be the interpolated orography,
!    but this is unnecessary)
!  ZCMP(J,3) interpolated ice cap fraction
!  ZCMP(J,4/6/8/...) interpolated monthly surface temperature
!  ZCMP(J,5/7/9/...) interpolated monthly surface Td

!  Allocate array for interpolated fields
ALLOCATE ( ZCMP(ITFING,JPNFLD) )

!  Interpolation of mask ( no mask on sea)
CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,1),IDATY,IDATX,1,ZCMP(1,1),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  Interpolation of glacier fraction (4 points - no mask on sea)
CALL EINTER8(YDGEOMETRY%YREGEO,ZFLD(1,1,3),IDATY,IDATX,1,ZCMP(1,3),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!  Interpolation of T and Td (12/4 points - mask on sea)
INFL=JPNFLD-JPNFIX
CALL EINTER6(YDGEOMETRY%YREGEO,ZFLD(1,1,JPNFIX+1),IDATY,IDATX,INFL,ZCMP(1,JPNFIX+1),&
 & IYFING,IXFING,ITFING,EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT,&
 & ZFLD(1,1,1),ZLSM(1),YRCLI%SMASK)  

!  Deallocate array for initial data
DEALLOCATE ( ZFLD )
!  Deallocate array for reading input
DEALLOCATE ( ZRBUF )

!     ------------------------------------------------------------------

!     4. FINAL CALCULATIONS.
!        -------------------

!  Monthly dependent values.
!   IMF=(JMOIS-1)*IPNF address of month in ARGEGE fields
!   IMD=(JMOIS-1)*JPNVAR + JPNFIX address of month in data
!    +(1-ICS) -> Surface Temperature, followed by deep values
!    +1+ICS   -> Soil moisture
!    +2+ICS   -> Deep soil moisture
!    +3+ICS   -> Snow depth
!    +4+ICS   -> Relaxation temperature
!    +5+ICS   -> Relaxation soil moisture

!     4.1 Allocate array for final fields (C+I zone)

ALLOCATE ( ZRES(ITFING,(YSP_SBD%NLEVS+JPNC)*12) )

!     4.2 Surface fields

DO JM=1,12
  IMF=(JM-1)*IPNF
  IMD=JPNFIX+(JM-1)*JPNVAR
  DO J=1,ITFING
    IF (ZLSM(J) >= YRCLI%SMASK) THEN
      ZTS=ZCMP(J,IMD+1)+ZGST*ZGEO(J)
      ZRES(J,IMF+1)=ZTS
      ZTD=ZCMP(J,IMD+2)+ZGST*ZGEO(J)
      ZHU=FOHU(MIN(ZTD,ZTS),ZTS)
      ZWR=FOWR(MAX(0.0_JPRB,MIN(1.0_JPRB,ZHU)))
      ZXSN=MAX(0.0_JPRB,MIN(1.0_JPRB,ZCMP(J,3)/MAX(ZEPS,ZCMP(J,1))))
      ZYSN=MAX(0.0_JPRB,MIN(1.0_JPRB,ZNEI1*ZWR*(RTT-ZTS)))
      ZRES(J,IMF+IWR)=ZWR
      ZRES(J,IMF+ISN)=FOSN(ZXSN,ZYSN)
    ELSE
      ZRES(J,IMF+1)=ZCMP(J,IMD+3)+ZGST*ZGEO(J)
      ZRES(J,IMF+IWR)=1.0_JPRB
      ZRES(J,IMF+ISN)=0.0_JPRB
    ENDIF
  ENDDO
ENDDO

!     4.3 Time lag and attenuation of deep soil values.

ILEN2=IPNF*11+2

!  Deep temperature ; relaxation temperature if old scheme
DO JCS=1,ICS-1
  CALL INCLAG(ZRES(1,JCS),ZTSHT(JCS),ITFING,ITFING,ILEN2,IPNF)
ENDDO

!  Deep moisture ; relaxation moisture if old scheme
ZWSHT=PPWSHT
CALL INCLAG(ZRES(1,IWR),ZWSHT,ITFING,ITFING,ILEN2,IPNF)
IF (.NOT.LLNEW) THEN
  ZWSHT=PPRELA*PPWSHT
  CALL INCLAG(ZRES(1,IWR+1),ZWSHT,ITFING,ITFING,ILEN2,IPNF)
ENDIF

!  No lag over sea
DO JM=1,12
  IMF=(JM-1)*IPNF
  DO JCS=2,ICS
    DO J=1,ITFING
      IF (ZLSM(J) < YRCLI%SMASK) ZRES(J,IMF+JCS)=ZRES(J,IMF+1)
    ENDDO
  ENDDO
  DO J=1,ITFING
    IF (ZLSM(J) < YRCLI%SMASK) ZRES(J,IMF+IWR+1)=ZRES(J,IMF+IWR)
  ENDDO
  IF (.NOT.LLNEW) THEN
    DO J=1,ITFING
      IF (ZLSM(J) < YRCLI%SMASK) ZRES(J,IMF+IWR+2)=ZRES(J,IMF+IWR)
    ENDDO
  ENDIF
ENDDO

!  Relaxation temperature and moisture (new schames)
IF (LLNEW) THEN
  DO JM=1,12
    IMF=(JM-1)*IPNF
    DO J=1,ITFING
      ZRES(J,IMF+ICS+4)=ZRES(J,IMF+2    )
      ZRES(J,IMF+ICS+5)=ZRES(J,IMF+IWR+1)
    ENDDO
  ENDDO
ENDIF

!     4.4 Corrections according to land use type

IF (LLTYP) THEN
  DO JM=1,12
    IMF=(JM-1)*IPNF
    DO J=1,ITFING
!  On land : snow fraction on ice-cap
      IF (NINT(ZITP(J)) == YRCLI%NTPGLA) THEN
        ZRES(J,IMF+ISN)= ZNEI5
      ELSE
        ZRES(J,IMF+ISN)= MIN(ZNEI4,ZRES(J,IMF+ISN))
      ENDIF
!  On water : albedo, emissivity and roughness length
      IF (NINT(ZITP(J)) == YRCLI%NTPMER .OR. NINT(ZITP(J)) == YRCLI%NTPLAC) THEN
        IF (ZRES(J,IMF+1) <= TMERGL) THEN
          ZRES(J,IMF+IPNF-2)= YRCLI%SALBB
          ZRES(J,IMF+IPNF-1)= YRCLI%SEMIB
          ZRES(J,IMF+IPNF  )= YRCLI%SZZ0B*RG
        ELSE
          ZRES(J,IMF+IPNF-2)= YRCLI%SALBM
          ZRES(J,IMF+IPNF-1)= YRCLI%SEMIM
          ZRES(J,IMF+IPNF  )= YRCLI%SZZ0M*RG
        ENDIF
      ELSE
        ZRES(J,IMF+IPNF-2)= ZALB(J)
        ZRES(J,IMF+IPNF-1)= ZEMI(J)
        ZRES(J,IMF+IPNF  )= ZGZ0(J)
      ENDIF
    ENDDO
  ENDDO
ENDIF

!     4.5 Deallocate temporary arrays

DEALLOCATE ( ZCMP )
DEALLOCATE ( ZITP )
DEALLOCATE ( ZLSM )
DEALLOCATE ( ZGEO )
DEALLOCATE ( ZGZ0 )
DEALLOCATE ( ZEMI )
DEALLOCATE ( ZALB )

!     ------------------------------------------------------------------

!     5. WRITING.
!        --------

!     5.1 Name of records

CLPREF(0    )='SURF    '
CLSUFF(0    )='IND.TERREMER'
CLPREF(1    )='SURF    '
CLPREF(2    )='PROF    '
CLSUFF(1    )='TEMPERATURE '
CLSUFF(2    )='TEMPERATURE '
CLPREF(IWR  )='SURF    '
CLPREF(IWR+1)='PROF    '
CLSUFF(IWR  )='PROP.RMAX.EA'
CLSUFF(IWR+1)='PROP.RMAX.EA'
CLPREF(ISN  )='SURF    '
CLSUFF(ISN  )='RESERV.NEIGE'
IF (LLNEW) THEN
  CLPREF(ICS+4)='RELA    '
  CLPREF(ICS+5)='RELA    '
  CLSUFF(ICS+4)='TEMPERATURE '
  CLSUFF(ICS+5)='PROP.RMAX.EA'
  DO JCS=3,ICS
    WRITE(CLPREF(JCS),'(''LEV'',I1)') JCS
    CLSUFF(JCS  )='TEMPERATURE'
  ENDDO
ELSE
  CLPREF(3    )='RELA    '
  CLPREF(IWR+2)='RELA    '
  CLSUFF(3    )='TEMPERATURE '
  CLSUFF(IWR+2)='PROP.RMAX.EA'
ENDIF

CLSUFF(IPNF-2)='ALBEDO      '
CLSUFF(IPNF-1)='EMISSIVITE  '
CLSUFF(IPNF  )='Z0.FOIS.G   '

DO JCS=IPNF-2,IPNF
  CLPREF(JCS)='SURF    '
  INIVL(JCS)=0
ENDDO

LLPAC(IPNF)=.FALSE.

!     5.2 Final biperiodization and writing on CLIM files as results.

!  Effective number of fields
IF (LLTYP) THEN
  IFLD=IPNF
ELSE
  IFLD=IPNF-3
ENDIF

IFL=0
DO JM=1,12

  WRITE(UNIT=CLNOMF,FMT='(''Const.Clim.'',I2.2)') JM
  IDATEF(2)=JM
  CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'OLD',.TRUE.,LLIMST,&
   & IMES,IARP,IARI,CLNOMC)  
  CALL FANDAR(IREP,INUM,IDATEF)

  DO JC=1,IFLD
    IFL=IFL+1
    DO J=1,ITFING
      ZEXT(J,JC)=ZRES(J,IFL)
    ENDDO
  ENDDO
  IFL=IFL+IPNF-IFLD

  CALL EBICLI(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF,IFLD+1,INIVL(0),CLPREF(0),CLSUFF(0),INUM,ZEXT(1,0),&
   & LLBIP(0),LLWRI(0),LLPAC(0))

  IF (JM == 1) THEN
    CALL LFILAF(IREP,INUM,.TRUE.)
    LLIMST=.FALSE.
  ENDIF
  CALL FAIRME(IREP,INUM,'KEEP')

ENDDO

DEALLOCATE ( ZRES )

!------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EINCLI3',1,ZHOOK_HANDLE)
END SUBROUTINE EINCLI3

