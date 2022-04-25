SUBROUTINE CHK923(YDGEOMETRY,KNUM)

!**** *CHK923*

!     PURPOSE
!     -------

!     This routine checks climatological fields on unit KNUM.

!**   INTERFACE
!     ---------

!     CALL CHK923(...) (called by INCLIn, n > 3)

!     METHOD
!     ------

!     EXTERNALS
!     ---------

!     FA-LFI package (FAITOU,FACILE,FAIENC,LFILAF,FAIRME)
!     CHIEN
!     ABOR1

!     AUTHORS
!     -------
!     D. Giard 97-03-13

!     MODIFICATIONS
!     -------------
!     M.Hamrud 03-10-01 CY28 Cleaning
!     K. Stadlbacher 03-01-07 bug correction
!     D. Giard 04-09-15 phasing and translating
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMLUN   , ONLY : NULOUT
USE YOMCLI   , ONLY : YRCLI

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KNUM

! JPTYVE : number of land-use types, at least 3 :
!   1 -> sea , 2 -> ice , 3 -> desert or low vegetation
!    (forests are 4, lakes are 1 or 5)
! JPNFIX : number of fixed climatological fields to check
! JPNCLI : number of monthly climatological fields (relaxation values)
! JPNVAR : number of monthly climatological fields (surface characteristics)


INTEGER(KIND=JPIM) :: JPNCLI
INTEGER(KIND=JPIM) :: JPNFIX
INTEGER(KIND=JPIM) :: JPNTOT
INTEGER(KIND=JPIM) :: JPNVAR
INTEGER(KIND=JPIM) :: JPTYVE

PARAMETER(JPTYVE=5,JPNFIX=9,JPNCLI=7,JPNVAR=8)
PARAMETER(JPNTOT=JPNFIX+JPNCLI+JPNVAR)


REAL(KIND=JPRB) :: ZRES(YDGEOMETRY%YRDIM%NDGLG*YDGEOMETRY%YRDIM%NDLON,JPNTOT),ZS(JPNTOT,JPTYVE,4)
REAL(KIND=JPRB) :: ZEPS, ZLND
REAL(KIND=JPRB) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: IPTS(JPTYVE)
INTEGER(KIND=JPIM) :: I, INIV, IPT, IREP, IV, IVG, IVX, J, JF, JJ, JS, JT

CHARACTER :: CLPREF(JPNTOT)*4,CLSUFF(JPNTOT)*12

LOGICAL :: LLCOSP

#include "abor1.intfb.h"

!     CHARACTER :: CLNOMC*16,CLNOMF*10

!     ------------------------------------------------------------------

!     1. SET INITIAL VALUES.
!        -------------------

IF (LHOOK) CALL DR_HOOK('CHK923',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NDGLG=>YDDIM%NDGLG, NDLON=>YDDIM%NDLON, &
 & NLOENG=>YDGEM%NLOENG)

DO JT=1,JPTYVE
  IPTS(JT)=0
  DO J=1,JPNTOT
    ZS(J,JT,1)=+1.E+10_JPRB
    ZS(J,JT,2)=-1.E+10_JPRB
    ZS(J,JT,3)=0.0_JPRB
    ZS(J,JT,4)=0.0_JPRB
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!     2. READ ARPEGE FIELDS.
!        -------------------

!     IARI=0
!     IARP=8
!     IMES=1
!     CLNOMF='Const.Clim'
!     CLNOMC='Const.Clim.Surfa'
!     LLIMST=.TRUE.

!     IINF=-1
!     ZEPS=1.E-10
!     LLPOLE=.TRUE. 

IREP=0
INIV=0
LLCOSP=.FALSE.

!     CALL FAITOU(IREP,KNUM,.TRUE.,CLNOMF,'OLD',.TRUE.,LLIMST,&
!      & IMES,IARP,IARI,CLNOMC)
!     CALL CHIEN(CLNOMC,NSTTYP,RMUCEN,RLOCEN,RSTRET,NSMAX,&
!      & NDGLG,NDLON,NLOENG,NMENG,RMU,NHTYP,NFLEVG,VP00,VALH,VBH,&
!      & NQUAD,IINF,NDGSAG,NDGENG,ZEPS,LLPOLE,NULOUT)

DO J=1,JPNTOT
  CLPREF(J)='SURF'
ENDDO
DO J=2,4,2
  CLPREF(JPNFIX+J  )='PROF'
  CLPREF(JPNFIX+J+1)='RELA'
ENDDO
! computed by INCLI1
CLSUFF(1)='PROP.TERRE  '
! computed by INCLI2
CLSUFF(2)='IND.VEG.DOMI'
CLSUFF(3)='ALBEDO.SOLNU'
CLSUFF(4)='EMISSIVITE  '
CLSUFF(5)='EPAI.SOL.MAX'
CLSUFF(6)='PROP.ARGILE '
CLSUFF(7)='PROP.SABLE  '
CLSUFF(8)='PROP.VEG.MAX'
CLSUFF(9)='EPAIS.SOL   '
! computed by INCLI3
CLSUFF(JPNFIX+1)='TEMPERATURE '
CLSUFF(JPNFIX+2)='TEMPERATURE '
CLSUFF(JPNFIX+3)='TEMPERATURE '
CLSUFF(JPNFIX+4)='PROP.RMAX.EA'
CLSUFF(JPNFIX+5)='PROP.RMAX.EA'
CLSUFF(JPNFIX+6)='PROP.RMAX.EA'
CLSUFF(JPNFIX+7)='RESERV.NEIGE'
! computed by INCLI4
CLSUFF(JPNFIX+JPNCLI+1)='PROP.VEGETAT'
CLSUFF(JPNFIX+JPNCLI+2)='Z0.FOIS.G   '
CLSUFF(JPNFIX+JPNCLI+3)='ALBEDO      '
CLSUFF(JPNFIX+JPNCLI+4)='IND.FOLIAIRE'
CLSUFF(JPNFIX+JPNCLI+5)='RESI.STO.MIN'
CLSUFF(JPNFIX+JPNCLI+6)='GZ0.THERM   '
CLSUFF(JPNFIX+JPNCLI+7)='Z0VEG.FOIS.G'
CLSUFF(JPNFIX+JPNCLI+8)='ALBEDO.VEG  '

DO J=1,JPNTOT
  CALL FACILE(IREP,KNUM,CLPREF(J),INIV,CLSUFF(J),ZRES(1,J),LLCOSP)
ENDDO

!     ------------------------------------------------------------------

!     3. CHECK FIELDS.
!     ----------------

IVX= 8
IVG= JPNFIX+JPNCLI+1
ZEPS= 1.E-10_JPRB

I=0
DO JJ=1,NDGLG
  DO J=1,NLOENG(JJ)
    I=I+1
    IV=NINT(ZRES(I,2))
!  fraction of vegetation over land
    IF (IV > 2) THEN
      ZLND=MAX(ZEPS,ZRES(I,1))
      ZRES(I,IVX)=MIN(1.0_JPRB,ZRES(I,IVX)/ZLND)
      ZRES(I,IVG)=MIN(1.0_JPRB,ZRES(I,IVG)/ZLND)
    ENDIF
!  loop over land use type
    DO JT=1,JPTYVE
      IF (IV == JT) THEN
        IPTS(JT)=IPTS(JT)+1
!  loop over fields
        DO JF=3,JPNTOT
          ZS(JF,JT,1)= MIN(ZS(JF,JT,1),ZRES(I,JF))
          ZS(JF,JT,2)= MAX(ZS(JF,JT,2),ZRES(I,JF))
          ZS(JF,JT,3)= ZS(JF,JT,3) + ZRES(I,JF)
          ZS(JF,JT,4)= ZS(JF,JT,4) + ZRES(I,JF)**2
        ENDDO
!  fraction of desert
        IF (ZRES(I,JPNFIX+JPNCLI+1) < YRCLI%SVEG)ZS(2,JT,3)=ZS(2,JT,3)+1.0_JPRB
      ENDIF
    ENDDO
  ENDDO
ENDDO

IPT=0
DO JT=1,JPTYVE
  IPT=IPT+IPTS(JT)
  ZS(2,JT,3)= ZS(2,JT,3)/MAX(1,IPTS(JT))
  DO JF=3,JPNTOT
    ZS(JF,JT,3)= ZS(JF,JT,3)/MAX(1,IPTS(JT))
    ZS(JF,JT,4)=&
     & SQRT(MAX(0.0_JPRB,ZS(JF,JT,4)/MAX(1,IPTS(JT)) - ZS(JF,JT,3)**2 ))  
  ENDDO
ENDDO
IPT=MAX(1,IPT)
IF (IPT /= I) CALL ABOR1(' CHK923 : UNKNOWN LAND USE TYPES !' )

!     ------------------------------------------------------------------

!     4. WRITE RESULTS.
!     -----------------

WRITE(NULOUT,FMT=411)
DO JT=1,JPTYVE
  WRITE(NULOUT,FMT=412) JT,IPTS(JT),IPT,ZS(2,JT,3)
  DO J=3,JPNTOT
    WRITE(NULOUT,FMT=413) CLPREF(J),CLSUFF(J),(ZS(J,JT,JS),JS=1,4)
  ENDDO
ENDDO

411 FORMAT(/,' CONFIGURATION 923',/,' STATISTICS ON SURFACE FIELDS')  
412 FORMAT(//,' LAND USE TYPE :',I2,/,&
 & ' number of points :',I6,'/',I6,' proportion of desert :',F5.2,/,&
 & '    field ',10X,'min ',6X,'max ',6X,'mean',6X,'std ',3X)  
413 FORMAT(1X,A4,A12,4G10.3)

!      CALL FAIRME(IREP,KNUM,'KEEP')

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CHK923',1,ZHOOK_HANDLE)
END SUBROUTINE CHK923
