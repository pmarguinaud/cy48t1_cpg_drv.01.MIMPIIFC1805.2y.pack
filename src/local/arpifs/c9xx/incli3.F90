SUBROUTINE INCLI3(YDGEOMETRY,YDSURF,YDPHY,YDPHY1)

!**** *INCLI3*

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
!     Albedo and emissivity over sea are modified according to the extension
!     of sea-ice, and snow cover is modified according to land use type (ice
!     cap or not), provided that these fields are present in the initial file.

!**   INTERFACE.
!     ----------

!     CALL INCLI3

!     METHOD.
!     -------

!     The  source dataset is read on unit 2.
!     Input data are coming from :
!    -an old NCAR climatology for surface temperature and moisture
!       and also temperature of the sea-ice
!    -AMIP1 mean data for sea surface temperature and sea-ice extent
!       and also temperature of the lake (reduced to sea level)
!     The target grid is a Gaussian grid with rotation of the pole and "Schmidt"
!     stretching as well as -if asked- reduction of the number of grid points
!     towards the poles of the final representation.
!     A subgrid of the Gaussian grid is introduced for interpolation
!     The values on the subgrid is obtained by a 12-point interpolation operator
!     applied on the source grid.
!     The interpolation is done separately for sea/sea-ice and land temperature.
!     The values on the Gaussian grid are the averages of the values in the
!     boxes (a box is the part of the subgrid corresponding to a point of the
!     Gaussian grid).
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
!     *The sea surface temperature interpolation scheme with a 3-bin mask
!       (land / sea / sea-ice)

!     EXTERNALS.
!     ----------

!     GEO923
!     INTER6
!     INTER8
!     INCLAG
!     FA-LFI package (FAITOU,FACILE,FAIENC,LFILAF,FAIRME)
!     CHIEN

!     AUTHORS.
!     --------
!      M. DEQUE      MARCH 1997

!     MODIFICATIONS.
!     --------------
!      D. Giard  14/11/2001 : roughness length of sea-ice;  preparing the case 
!        when the initial albedo does not fit the minimum sea-ice extension; 
!        small corrections for fields headers; more allocatable arrays
!      R. El Khatib 18/08/2003 : Roughness lengths not packed
!      M.Hamrud  01/10/2003 : CY28 Cleaning
!      D. Giard  15/09/2004 : cosmetic changes after automatic ones 
!      K. Yessad (Jan 2010): externalisation of group EGGX in XRD/IFSAUX
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMPHY   , ONLY : TPHY
USE YOMCT0   , ONLY : NQUAD
USE YOMCST   , ONLY : RPI, RG, RV, RCPV, RETV, RCW, RCS, RLVTT, RLSTT,&
 & RTT, RALPW, RBETW, RGAMW, RALPS, RBETS, RGAMS, RALPD, RBETD, RGAMD  
USE YOMVERT  , ONLY : VP00
USE YOMPHY1  , ONLY : TPHY1
USE YOMSTA   , ONLY : RDTDZ1
USE YOMCLI   , ONLY : YRCLI  

IMPLICIT NONE

! JPNFIX : number of fixed fields in dataset
! JPNVAR : number of monthly fields in dataset
! JPNFLD : number of fields in dataset
! JPBX/JPBY : number of extra longitudes/latitudes on each side of the data
!          grid required by interpolation (4 or 12 points -> 2/2)

TYPE(GEOMETRY),INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)   ,INTENT(INOUT) :: YDSURF
TYPE(TPHY)    ,INTENT(INOUT) :: YDPHY
TYPE(TPHY1)   ,INTENT(INOUT) :: YDPHY1
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

REAL(KIND=JPRB) ,ALLOCATABLE :: ZFLD(:,:,:)
REAL(KIND=JPRB) ,ALLOCATABLE :: ZALB(:), ZEMI(:), ZGZ0(:),&
 & ZGEO(:), ZLSM(:), ZITP(:)  
REAL(KIND=JPRB) ,ALLOCATABLE :: ZCMP(:,:), ZRES(:,:)
REAL(KIND=JPRB) :: ZTSHT(YDSURF%YSP_SBD%NLEVS+1),ZBO(YDGEOMETRY%YRDIM%NDGLG+1),ZLONG(YDGEOMETRY%YRDIM%NDGLG),&
 & ZMU(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON)&
 & ,ZSLA(YDGEOMETRY%YRDIM%NDGLG*YRCLI%NPINT),ZSLO(YDGEOMETRY%YRDIM%NDLON*YRCLI%NPINT,YDGEOMETRY%YRDIM%NDGLG),&
 & ZCLO(YDGEOMETRY%YRDIM%NDLON*YRCLI%NPINT,YDGEOMETRY%YRDIM%NDGLG)  
REAL(KIND=JPRB) :: PHUARG, PTDARG, PXARG, PYARG, ZEPS, ZGST,&
 & ZHU, ZNEI1, ZNEI2, ZNEI3, ZNEI4, ZNEI5, ZTD, ZTS, ZWR, ZWSHT, ZXSN, ZYSN  

REAL(KIND=JPRB) :: FOHU
REAL(KIND=JPRB) :: FOSN
REAL(KIND=JPRB) :: FOWR
REAL(KIND=JPRB) :: ZHOOK_HANDLE

CHARACTER :: CLPREF(YDSURF%YSP_SBD%NLEVS+JPNC)*8,CLSUFF(YDSURF%YSP_SBD%NLEVS+JPNC)*12&
 & ,CLNOMC*16,CLNOMF*16,CLFORM*12  

INTEGER(KIND=JPIM) :: IDATEF(11),INIVL(0:YDSURF%YSP_SBD%NLEVS+JPNC)
INTEGER(KIND=JPIM) :: IARI, IARP, IBIT, ICO, ICS, IDATX, IDATY,&
 & IFL, IFLD, IGRB, IINF, ILAP, ILEN2, ILENE,&
 & ILENT, IMD, IMES, IMF, INFL, INIQ, INIV,&
 & INJQ, INUM, IPNF, IPOSI, IREP, ISN, ITP,&
 & ITRN, IWR, J, JCS, JF, JJ, JMOIS  
INTEGER(KIND=JPIM) :: INGRIG,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL

LOGICAL :: LLCOSP, LLIMST, LLNEW, LLPOLE, LLTYP
LOGICAL :: LLPACK(YDSURF%YSP_SBD%NLEVS+JPNC)

!     ------------------------------------------------------------------

#include "fcttrm.func.h"

! Relative humidity function (from Td and T).
FOHU(PTDARG,PTARG)=ES(PTDARG)/ES(PTARG)

! Soil relative water content function (from surface relative humidity).
FOWR(PHUARG)=ACOS(1.0_JPRB-2.0_JPRB*PHUARG)/RPI

! Empirical snow climatology function.
FOSN(PXARG,PYARG)=MAX(0.0_JPRB,(ZNEI2*PXARG)/(ZNEI3-PXARG)&
 & +ZNEI4*(1.0_JPRB-COS(RPI*PYARG))/2.0_JPRB)  

!     ------------------------------------------------------------------

#include "chien.h"

#include "abor1.intfb.h"
#include "geo923.intfb.h"
#include "inclag.intfb.h"
#include "inter6.intfb.h"
#include "inter8.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INCLI3',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB)
ASSOCIATE(SODELX=>YDPHY1%SODELX, TMERGL=>YDPHY1%TMERGL,   NDGENG=>YDGEOMETRY%YRDIM%NDGENG, NDGLG=>YDGEOMETRY%YRDIM%NDGLG, &
& NDGSAG=>YDGEOMETRY%YRDIM%NDGSAG,   NDLON=>YDGEOMETRY%YRDIM%NDLON, NSMAX=>YDGEOMETRY%YRDIM%NSMAX,                        &
& NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, LSOLV=>YDPHY%LSOLV,   NHTYP=>YDGEOMETRY%YRGEM%NHTYP, NLOENG=>YDGEOMETRY%YRGEM%NLOENG, &
& NMENG=>YDGEOMETRY%YRGEM%NMENG,   NSTTYP=>YDGEOMETRY%YRGEM%NSTTYP, RLOCEN=>YDGEOMETRY%YRGEM%RLOCEN,                      &
& RMUCEN=>YDGEOMETRY%YRGEM%RMUCEN,   RSTRET=>YDGEOMETRY%YRGEM%RSTRET, YSP_SBD=>YDSURF%YSP_SBD)
!     ------------------------------------------------------------------

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
!  ZNEI2/(ZNEI3-1.) gives the maximum snow equivalent for total permanent
!  ice cover.
ZNEI3=1.0001_JPRB
!  Maximum non-permanent snow value.
ZNEI4=1000._JPRB
!  Permanent snow value on ice cap.
ZNEI5=10000._JPRB

!     1.3 Interpolation

ICO=YRCLI%NPINT*YRCLI%NPINT
INJQ=NDGLG*YRCLI%NPINT
INIQ=NDLON*YRCLI%NPINT
ILENT=NDGLG*NDLON
IDATY=YRCLI%NDATY+2*JPBY
IDATX=YRCLI%NDATX+2*JPBX

!     1.4 ARPEGE and data files

INUM=3
IMES=1
IARP=20
IREP=0
IARI=0
LLIMST=.TRUE.
CLNOMC='Const.Clim.Surfa'

ZEPS=1.E-10_JPRB
IINF=-1
LLPOLE=.TRUE.
INIV=0
LLCOSP=.FALSE.
IGRB=0
IBIT=0
ITRN=0
ILAP=0
LLPACK(:)=.TRUE.

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

!     1.5 Geometry

ILENE=0
CALL GEO923(YDGEOMETRY,YRCLI%NPINT,ILENE,ZBO,ZLONG,ZMU,ZSLA,ZSLO,ZCLO)

!     ------------------------------------------------------------------

!     2. READING  DATA.
!        --------------

!     2.1 Source data

!  Fixed fields:
!   Fraction of land , Orography*g , Fraction of the land covered by ice caps
!  Monthly fields:
!   Temperature and Td at mean orography level, SST (K)

!  Allocate array for initial data
ALLOCATE ( ZFLD(IDATX,IDATY,JPNFLD) )

OPEN(UNIT=10,FILE='N108_GL',FORM=CLFORM)
IF (YRCLI%LIEEE) THEN
  DO JF=1,JPNFLD
    DO JJ=1+JPBY,YRCLI%NDATY+JPBY
      READ(10) (ZFLD(J,JJ,JF),J=1+JPBX,YRCLI%NDATX+JPBX)
    ENDDO
  ENDDO
ELSE
  DO JF=1,JPNFLD
    READ(10,*)((ZFLD(J,JJ,JF),J=1+JPBX,YRCLI%NDATX+JPBX),JJ=1+JPBY,YRCLI%NDATY+JPBY)
  ENDDO
ENDIF
CLOSE(10)

!  T and Td are reduced to sea level before interpolation
DO JMOIS=1,12
  IPOSI=JPNFIX+(JMOIS-1)*JPNVAR
  DO JJ=JPBY+1,JPBY+YRCLI%NDATY
    DO J=JPBX+1,JPBX+YRCLI%NDATX
      ZFLD(J,JJ,IPOSI+1)= ZFLD(J,JJ,IPOSI+1) - ZGST*ZFLD(J,JJ,2)
      ZFLD(J,JJ,IPOSI+2)= ZFLD(J,JJ,IPOSI+2) - ZGST*ZFLD(J,JJ,2)
    ENDDO
  ENDDO
ENDDO

!     2.2 Initial ARPEGE clim file

!  Open file and check geometry
CLNOMF='Const.Clim'
CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'OLD',.TRUE.,LLIMST,&
 & IMES,IARP,IARI,CLNOMC)  
CALL CHIEN(CLNOMC,NSTTYP,RMUCEN,RLOCEN,RSTRET,NSMAX,&
 & NDGLG,NDLON,NLOENG,NMENG,NHTYP,NFLEVG,VP00,YDVAB%VALH,YDVAB%VBH,&
 & NQUAD,IINF,NDGSAG,NDGENG,ZEPS,LLPOLE,NULOUT)  
IF (LLPOLE) CALL ABOR1(' CLIM. FILES MUST NOT HAVE POLES !')

!  Allocate arrays for initial fields
ALLOCATE ( ZALB(ILENT) )
ALLOCATE ( ZEMI(ILENT) )
ALLOCATE ( ZGZ0(ILENT) )
ALLOCATE ( ZGEO(ILENT) )
ALLOCATE ( ZLSM(ILENT) )
ALLOCATE ( ZITP(ILENT) )

!  Read land-sea mask and height
INIV=0
CALL FACILE(IREP,INUM,'SURF',INIV,'IND.TERREMER',ZLSM,LLCOSP)
CALL FACILE(IREP,INUM,'SURF',INIV,'GEOPOTENTIEL',ZGEO,LLCOSP)

!  Read land use type, albedo roughness and emissivity, if already computed
LLTYP=.TRUE.
CALL FANION(IREP,INUM,'SURF',INIV,'IND.VEG.DOMI',LLTYP,LLCOSP,&
 & IGRB,IBIT,ITRN,ILAP)  
IF (LLTYP) THEN
  CALL FACILE(IREP,INUM,'SURF',INIV,'IND.VEG.DOMI',ZITP,LLCOSP)
  CALL FACILE(IREP,INUM,'SURF',INIV,'ALBEDO      ',ZALB,LLCOSP)
  CALL FACILE(IREP,INUM,'SURF',INIV,'EMISSIVITE  ',ZEMI,LLCOSP)
  CALL FACILE(IREP,INUM,'SURF',INIV,'Z0.FOIS.G   ',ZGZ0,LLCOSP)
  WRITE(NULOUT,'('' INCLI3 :'',&
   & '' MODIFICATION OF SNOW COVER ACCORDING TO ICE CAP,'',&
   & '' OF ALBEDO, EMISSIVITY, ROUGHNESS ACCORDING TO SEA-ICE'')')  
ELSE
  WRITE(NULOUT,'('' INCLI3 : LAND USE TYPE NOT AVAILABLE'')')
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
!  ZCMP(J,4/7/10/...) interpolated monthly surface temperature
!  ZCMP(J,5/8/11/...) interpolated monthly surface Td
!  ZCMP(J,6/9/12/...) interpolated monthly SST

!  Allocate array for interpolated fields
ALLOCATE ( ZCMP(ILENT,JPNFLD) )

!  Interpolation of mask (4 points -  no mask on sea)
CALL INTER8(NDGLG,NDLON,IDATY,IDATX,1,YRCLI%NPINT,ILENE,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZCMP(1,1),ZFLD(1,1,1),&
 & ZSLA,ZSLO,ZCLO)  

!  Interpolation of glacier fraction (4 points - no mask on sea)
CALL INTER8(NDGLG,NDLON,IDATY,IDATX,1,YRCLI%NPINT,ILENE,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,ZCMP(1,3),ZFLD(1,1,3),&
 & ZSLA,ZSLO,ZCLO)  

!  Interpolation of T, Td, and SST (12/4 points - mask on sea)
INFL=JPNFLD-JPNFIX
CALL INTER6(NDGLG,NDLON,IDATY,IDATX,INFL,YRCLI%NPINT,ILENE,&
 & ICO,INJQ,INIQ,NLOENG(1:),ILENT,&
 & ZCMP(1,JPNFIX+1),ZFLD(1,1,JPNFIX+1),ZSLA,ZSLO,ZCLO,&
 & ZFLD(1,1,1),ZLSM(1),YRCLI%SMASK)  

!  Deallocate array for initial data
DEALLOCATE ( ZFLD )

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

ALLOCATE ( ZRES(ILENT,(YSP_SBD%NLEVS+JPNC)*12) )

!     4.2 Surface fields

DO JMOIS=1,12
  IMF=(JMOIS-1)*IPNF
  IMD=JPNFIX+(JMOIS-1)*JPNVAR
  DO J=1,ILENE
    IF (ZLSM(J) >= YRCLI%SMASK) THEN
      ZTS=ZCMP(J,IMD+1)+ZGST*ZGEO(J)
      ZTD=ZCMP(J,IMD+2)+ZGST*ZGEO(J)
      ZRES(J,IMF+1)=ZTS
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
  CALL INCLAG(ZRES(1,JCS),ZTSHT(JCS),ILENE,ILENT,ILEN2,IPNF)
ENDDO

!  Deep moisture ; relaxation moisture if old scheme
ZWSHT=PPWSHT
CALL INCLAG(ZRES(1,IWR),ZWSHT,ILENE,ILENT,ILEN2,IPNF)
IF (.NOT.LLNEW) THEN
  ZWSHT=PPRELA*PPWSHT
  CALL INCLAG(ZRES(1,IWR+1),ZWSHT,ILENE,ILENT,ILEN2,IPNF)
ENDIF

!  No lag over sea
DO JMOIS=1,12
  IMF=(JMOIS-1)*IPNF
  DO JCS=2,ICS
    DO J=1,ILENE
      IF (ZLSM(J) < YRCLI%SMASK) ZRES(J,IMF+JCS)=ZRES(J,IMF+1)
    ENDDO
  ENDDO
  DO J=1,ILENE
    IF (ZLSM(J) < YRCLI%SMASK) ZRES(J,IMF+IWR+1)=ZRES(J,IMF+IWR)
  ENDDO
  IF (.NOT.LLNEW) THEN
    DO J=1,ILENE
      IF (ZLSM(J) < YRCLI%SMASK) ZRES(J,IMF+IWR+2)=ZRES(J,IMF+IWR)
    ENDDO
  ENDIF
ENDDO

!  Relaxation temperature and moisture (new schemes)
IF (LLNEW) THEN
  DO JMOIS=1,12
    IMF=(JMOIS-1)*IPNF
    DO J=1,ILENE
      ZRES(J,IMF+ICS+4)=ZRES(J,IMF+2    )
      ZRES(J,IMF+ICS+5)=ZRES(J,IMF+IWR+1)
    ENDDO
  ENDDO
ENDIF

!     4.4 Corrections according to land use type

IF (LLTYP) THEN
  DO JMOIS=1,12
    IMF=(JMOIS-1)*IPNF
    DO J=1,ILENE
      ITP= NINT(ZITP(J))
!  Increased snow fraction on ice-cap
      IF (ITP == YRCLI%NTPGLA) THEN
        ZRES(J,IMF+ISN)= ZNEI5
      ELSE
        ZRES(J,IMF+ISN)= MIN(ZNEI4,ZRES(J,IMF+ISN))
      ENDIF
!  On water : albedo, emissivity and roughness length
      IF (ITP == YRCLI%NTPMER .OR. ITP == YRCLI%NTPLAC) THEN
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

DO JCS=1,IPNF-3
  INIVL(JCS)=1
ENDDO
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
LLPACK(IPNF)=.FALSE.
DO JCS=IPNF-2,IPNF
  CLPREF(JCS)='SURF    '
  INIVL(JCS)=0
ENDDO

!     5.2 Final writing on ARPEGE files as results.

!  Effective number of fields
IF (LLTYP) THEN
  IFLD=IPNF
ELSE
  IFLD=IPNF-3
ENDIF

IFL=0
DO JMOIS=1,12
  WRITE(UNIT=CLNOMF,FMT='(''Const.Clim.'',I2.2)') JMOIS
  IDATEF(2)=JMOIS
  CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'OLD',.TRUE.,LLIMST,IMES,&
   & IARP,IARI,CLNOMC)  
  CALL FANDAR(IREP,INUM,IDATEF)
  DO JF=1,IFLD
    IFL=IFL+1
    IF (.NOT.LLPACK(JF)) THEN
!     Do not pack roughness lengths
      CALL FAVEUR(IREP,INUM,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
      INGRIG=0
      CALL FAGOTE(IREP,INUM,INGRIG,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
    ENDIF
    CALL FAIENC(IREP,INUM,CLPREF(JF),INIVL(JF),CLSUFF(JF),ZRES(1,IFL),&
     & LLCOSP)  
    IF (.NOT.LLPACK(JF)) THEN
!     Reset packing after the treatment of roughness lengths
      CALL FAGOTE(IREP,INUM,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
    ENDIF
  ENDDO
  IFL=IFL+IPNF-IFLD
  IF (JMOIS == 1) THEN
    CALL LFILAF(IREP,INUM,.TRUE.)
    LLIMST=.FALSE.
  ENDIF
  CALL FAIRME(IREP,INUM,'KEEP')
ENDDO

DEALLOCATE ( ZRES )

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INCLI3',1,ZHOOK_HANDLE)
END SUBROUTINE INCLI3

