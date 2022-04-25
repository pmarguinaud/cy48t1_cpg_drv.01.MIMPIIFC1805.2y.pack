SUBROUTINE EINCLI9(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF)

!**** *EINCLI9*

!     PURPOSE.
!     --------

!      Interpolation or computation of aerosol fields

!**   INTERFACE.
!     ----------

!     CALL EINCLI9

!     METHOD.
!     -------

!      The data are read on unit 11, if available, then interpolated.
!      Finally they are written on ALADIN clim file (unit 3).

!     EXTERNALS.
!     ----------

!     ABOR1
!     EGEO923
!     EINTER2
!     EBICLI
!     FA-LFI package (FAITOU,LFILAF,FAIRME)
!     CCHIEN

!     AUTHORS.
!     --------
!      D. GIARD SEPT. 04 (from EINCLIR modified by F. BOUYSSEL)

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 01-12-06 Cleaning sm variables
!      Modified 03-10-01 by M.Hamrud : CY28 Cleaning
!      Modified 04-06-15 by F. Bouyssel : change of format and scaling
!                      for input climatological data, AEROS.URBAN
!                      by AEROS.SOOT accordingly
!      Modified 04-09-15 by D. Giard : new version renamed EINCLI9 and 
!                      phased EINCLIR kept for older versions
!      Modified 04-11-30 Cleaning pre29
!      Modified 05-04-07 by D. Giard : new call to EBICLI
!      O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!      Modified 13-03-18 by L. Rontu : CAMS AOD import in addition to Tegen
!     ------------------------------------------------------------------

USE MODEL_PHYSICS_MF_MOD , ONLY : MODEL_PHYSICS_MF_TYPE
USE YOEPHY       , ONLY : TEPHY
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMCLI   , ONLY : YRCLI
USE YOMLUN   , ONLY : NULOUT
USE YOM_YGFL , ONLY : TYPE_GFLD

!     ------------------------------------------------------------------

IMPLICIT NONE

!     JPBY : Extra-latitudes
!          2 when INTER2 or INTER3 are used (even NDATY)
!          1 when INTER4 or INTER5 are used (odd  NDATY)
!     JPNCH: Number of fields to be processed.
!     NAEROF: 0 = Tegen AOD550, 1 = CAMS AOD550, 2 = CAMS MMR (2 not yet active)
TYPE(GEOMETRY), INTENT(INOUT)   :: YDGEOMETRY
TYPE(TYPE_GFLD)    ,INTENT(INOUT):: YDGFL
TYPE(TEPHY)    ,INTENT(INOUT)   :: YDEPHY
TYPE(MODEL_PHYSICS_MF_TYPE),INTENT(INOUT):: YDML_PHY_MF
INTEGER(KIND=JPIM) :: JPBX
INTEGER(KIND=JPIM) :: JPBY
INTEGER(KIND=JPIM) :: JPNCH
!PARAMETER (JPBX=2,JPBY=2,JPNCH=4)
PARAMETER (JPBX=2,JPBY=2,JPNCH=11)

REAL(KIND=JPRB) :: ZS(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON,JPNCH)
REAL(KIND=JPRB),ALLOCATABLE :: ZFLD(:,:), ZRES(:), ZZLA(:), ZZLO(:)

REAL(KIND=JPRB) :: ZSCALE
REAL(KIND=JPRB) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: INIVL(JPNCH)
INTEGER(KIND=JPIM) :: IARI, IARP, ICO, IDATX, IDATY, IFLD, IINF, IMES,&
 & INUM, IOS, IREP, ITFING, IXFING, IYFING, J, JF, JJ, NNDATY, NNDATX

CHARACTER :: CLPRE(JPNCH)*8,CLSUF(JPNCH)*12
CHARACTER :: CLNOMC*16, CLNOMF*10, CLFORM*12
CHARACTER :: CID*4

LOGICAL :: LLBIP(JPNCH),LLWRI(JPNCH),LLPAC(JPNCH)
LOGICAL :: LLIMST, LLDATA, LLOPEN

INTEGER :: JJPNCH
!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "cchien.intfb.h"
#include "ebicli.intfb.h"
#include "egeo923.intfb.h"
#include "einter2.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EINCLI9',0,ZHOOK_HANDLE)

ASSOCIATE( &
 & NDGUNG=>YDGEOMETRY%YRDIM%NDGUNG, NDGUXG=>YDGEOMETRY%YRDIM%NDGUXG, NDLUNG=>YDGEOMETRY%YRDIM%NDLUNG, &
 & NDLUXG=>YDGEOMETRY%YRDIM%NDLUXG, NDGLG=>YDGEOMETRY%YRDIM%NDGLG, NDLON=>YDGEOMETRY%YRDIM%NDLON, &
 & EDELX=>YDGEOMETRY%YREGEO%EDELX, EDELY=>YDGEOMETRY%YREGEO%EDELY)
!     ------------------------------------------------------------------
!*
!     1. SET INITIAL VALUES.
!        -------------------

!     Constants

ZSCALE=1._JPRB

!     Grids.

ICO=YRCLI%NPINT*YRCLI%NPINT
IDATX=YRCLI%NDATX+2*JPBX
IDATY=YRCLI%NDATY+2*JPBY
IXFING=NDLUXG-NDLUNG+1
IYFING=NDGUXG-NDGUNG+1
ITFING=IXFING*IYFING

CALL EGEO923(YDGEOMETRY%YRGEM)

!     ALADIN fields

JJPNCH=JPNCH

IF(YRCLI%NAEROF==1.OR.YRCLI%NAEROF==0) THEN ! AOD
   JJPNCH=4
ENDIF

IF(YRCLI%NAEROF==1.OR.YRCLI%NAEROF==2) THEN ! CAMS
   NNDATX=YRCLI%NDATX
   NNDATY=YRCLI%NDATY+1
ENDIF

DO J=1,JJPNCH
  INIVL(J) =0
  LLBIP(J)=.TRUE.
  LLWRI(J)=.TRUE.
  LLPAC(J)=.TRUE.
  CLPRE(J)='SURF'
ENDDO
IF (YRCLI%NAEROF==0.OR.YRCLI%NAEROF==1) THEN
CLSUF(1)='AEROS.SEA   '
CLSUF(2)='AEROS.LAND  '
CLSUF(3)='AEROS.SOOT  '
CLSUF(4)='AEROS.DESERT'
!ELSE IF (NAEROF==2) THEN  !not yet possible
!SS1,SS2,SS3,DD1,DD2,DD3,OM1,OM2,BC1,BC2,SU
!   CLSUF(1)='AEROMMR.SS1 '
!   CLSUF(2)='AEROMMR.SS2 '
!   CLSUF(3)='AEROMMR.SS3 '
!   CLSUF(4)='AEROMMR.DD1 '
!   CLSUF(5)='AEROMMR.DD2 '
!   CLSUF(6)='AEROMMR.DD3 '
!   CLSUF(7)='AEROMMR.OM1 '
!   CLSUF(8)='AEROMMR.OM2 '
!   CLSUF(9)='AEROMMR.BC1 '
!   CLSUF(10)='AEROMMR.BC2 '
!   CLSUF(11)='AEROMMR.SU  '
ELSE
  CALL ABOR1('EINCLI9 : AEROSOL INPUT FILE NOT DEFINED!')
ENDIF
ZS(:,:)=0.0_JPRB

!     ALADIN file identifiers

INUM=3
IREP=0
IARI=0
IARP=1
CLNOMF='Const.Clim'
CLNOMC='Const.Clim.Surfa'
LLIMST=.TRUE.
IMES=1

IINF=0

!     Data type

IF (YRCLI%LIEEE) THEN
  CLFORM='UNFORMATTED'
  CALL ABOR1('EINCLI9 : STANDARD INPUT FILES NOT YET READY !')
ELSE
  CLFORM='FORMATTED'
ENDIF

!     ------------------------------------------------------------------
!*
!     2. OPEN FILES.
!        -----------

!     ALADIN clim file

CALL FAITOU(IREP,INUM,.TRUE.,CLNOMF,'OLD',.TRUE.,LLIMST,&
 & IMES,IARP,IARI,CLNOMC)  
CALL CCHIEN(YDGEOMETRY,CLNOMC,INUM,IINF)

!     External input data

LLDATA=.FALSE.
LLOPEN=.FALSE.
IOS= 0
INQUIRE(FILE='aero_GL',IOSTAT=IOS,EXIST=LLDATA,OPENED=LLOPEN)
LLDATA= LLDATA .AND. (IOS == 0) .AND. .NOT.LLOPEN

IF (LLDATA) THEN
  OPEN(UNIT=11,FILE='aero_GL',FORM=CLFORM)
ELSE
  CALL ABOR1('EINCLI9 : INPUT DATA NOT AVAILABLE !')
ENDIF

!     ------------------------------------------------------------------
!*
!     3. READ AND INTERPOLATE DATA.
!        --------------------------

!     Allocate temporary space

ALLOCATE ( ZFLD(IDATX,IDATY) )
!ALLOCATE ( ZZLA(NDATY) )
ALLOCATE ( ZZLO(YRCLI%NDATX) )
ALLOCATE ( ZRES(ITFING) )

!     Read file header
!     but the lats and lons are not used beyond this point

IF (YRCLI%NAEROF.eq.1.or.YRCLI%NAEROF.eq.2) THEN
!CAMS
   ALLOCATE ( ZZLA(NNDATY) )
   READ(11,*) CID, NNDATY
   WRITE(NULOUT,*) "eincli9 cams: NNDATY NDATY ",CID,NNDATY,YRCLI%NDATY
   READ(11,*) ZZLA
   READ(11,*) CID, NNDATX
   WRITE(NULOUT,*) "eincli9 cams: NNDATX NDATX ",CID,NNDATX,YRCLI%NDATX
   READ(11,*) ZZLO
ELSEIF (YRCLI%NAEROF.eq.0) THEN
!TEGEN
   ALLOCATE ( ZZLA(YRCLI%NDATY) )
READ(11,9901) ZZLA
READ(11,9901) ZZLO
9901 FORMAT(5F15.4)
ENDIF

!     Loop on fields

DO JF=1,JJPNCH

  ZFLD(:,:)=0.0_JPRB
  ZRES(:)=0.0_JPRB

!     Read data
IF (YRCLI%NAEROF.eq.1.or.YRCLI%NAEROF.eq.2) THEN
!CAMS AOD
  WRITE(NULOUT,*) 'eincli9 read cams aerosol: IY,IX,AERO'
  DO JJ=JPBY+YRCLI%NDATY,JPBY+1,-1
     DO J=JPBX+1,JPBX+YRCLI%NDATX
        READ (11,*) ZFLD(J,JJ)
        WRITE(NULOUT,*) JJ,J,ZFLD(J,JJ)
     ENDDO
  ENDDO
!in the CAMS files, there is data for 61 latitudes, not 60 so
!read the lat=-90 values to the last row (JPBY+2) already filled:
     DO J=JPBX+1,JPBX+YRCLI%NDATX
        READ (11,*) ZFLD(J,JPBY+1)
        WRITE(NULOUT,*) JPBY+1,J,ZFLD(J,JPBY+1)
     ENDDO
ELSEIF (YRCLI%NAEROF.eq.0) THEN
!TEGEN AOD
  DO JJ=JPBY+YRCLI%NDATY,JPBY+1,-1
    READ (11,9902) (ZFLD(J,JJ),J=JPBX+1,JPBX+YRCLI%NDATX)
  ENDDO
ENDIF

!     Horizontal interpolation (12 points)
  IFLD=1
  CALL EINTER2(YDGEOMETRY%YREGEO,ZFLD,IDATY,IDATX,IFLD,ZRES,IYFING,IXFING,ITFING,&
   & EDELX,EDELY,YRCLI%NPINT,ICO,NULOUT)  

!     Copy and scale interpolated data
  DO J=1,ITFING
    ZS(J,JF)=MAX(0.0_JPRB,ZRES(J)/ZSCALE)
  ENDDO

ENDDO
9902 FORMAT(12F10.6)

!     Deallocate temporary space

DEALLOCATE ( ZRES )
DEALLOCATE ( ZFLD )
DEALLOCATE ( ZZLA )
DEALLOCATE ( ZZLO )

!     Close input file

CLOSE(11)

!     ------------------------------------------------------------------
!*
!     4. WRITE FIELDS.
!        -------------

!     Biperiodize and write fields

IFLD=JJPNCH
CALL EBICLI(YDGEOMETRY,YDGFL,YDEPHY,YDML_PHY_MF,IFLD,INIVL,CLPRE,CLSUF,INUM,ZS,LLBIP,LLWRI,LLPAC)

!     Close ALADIN clim file

CALL LFILAF (IREP,INUM,.TRUE.)
CALL FAIRME (IREP,INUM,'KEEP')

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EINCLI9',1,ZHOOK_HANDLE)
END SUBROUTINE EINCLI9
