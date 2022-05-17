SUBROUTINE CASMSWI(YDGEOMETRY,YDSURF,YDPHY1)

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK
USE YOMCT0   , ONLY : LELAM
USE YOMPHY1  , ONLY : TPHY1
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : MYPROC   ,NPROC

USE QACVEG   , ONLY : NR_SM_WP ,RA_SM_WP

USE DISGRID_MOD, ONLY : DISGRID_SEND, DISGRID_RECV
USE DIWRGRID_MOD, ONLY : DIWRGRID_SEND, DIWRGRID_RECV

!**** **  - Spatial smoothing of SWI (Soil Wetness Index) for changing Wp

!     Purpose.
!     --------
!           Spatial smoothing of SWI (Soil Wetness Index) then changing of Wp 
!           (Total soil water content) in Canari OI 

!**   Interface.
!     ----------
!        *CALL* *CASMSWI

!        Explicit arguments :  
!        ------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        ACSOLW
!        DIWRGRID
!        LISLAP &(OR) ELISLAP
!        DISGRID

!     Reference.
!     ----------
!        HPOS, WRGRIDUA, WRGP2FA & INCLI1 CY24T1

!     Author.
!     -------
!        *Meteo-France*

!     Modifications.
!     --------------
!        Original : 27-09-02 S. Ivatek-Sahdan
!        Modified :
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      01-Jul-2006 Revised surface fields
!        A. Trojakova : 30/05/2006 new argument for ELISLAP
!               due to Bf for map factor addressing on multiproc
!        Apr 2008  K. Yessad: use DISGRID instead of DISGRID_C + cleanings
!        K. Yessad (Jul 2009): remove CDLOCK + some cleanings
!        R. El Khatib : 23-Apr-2010 use disgrid_mod & diwrgrid_mod
!        G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM and TCSGLEG
!        T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!        E.Dutra/G.Arduini (Jan 2018) Change in SP_SG dimensions
!     ------------------------------------------------------------------

IMPLICIT NONE
TYPE(GEOMETRY),INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)   ,INTENT(INOUT) :: YDSURF
TYPE(TPHY1)   ,INTENT(INOUT) :: YDPHY1
REAL(KIND=JPRB) :: ZBO(YDGEOMETRY%YRDIM%NDGLG)
REAL(KIND=JPRB) :: ZD2(YDGEOMETRY%YRGEM%NGPTOT), ZDSW(YDGEOMETRY%YRGEM%NGPTOT), ZDSWN(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZOCE(YDGEOMETRY%YRGEM%NGPTOT), ZSSW(YDGEOMETRY%YRGEM%NGPTOT), ZSWI(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZVEG(YDGEOMETRY%YRGEM%NGPTOT), ZWFC(YDGEOMETRY%YRGEM%NGPTOT), ZWPMX(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZWSAT(YDGEOMETRY%YRGEM%NGPTOT), ZWSMX(YDGEOMETRY%YRGEM%NGPTOT), ZWWILT(YDGEOMETRY%YRGEM%NGPTOT)
REAL(KIND=JPRB) :: ZEPW, ZFDSW, ZLSM, ZSD, ZMASK
REAL(KIND=JPRB) :: ZGOCE(YDGEOMETRY%YRGEM%NGPTOTG), ZGSWI(YDGEOMETRY%YRGEM%NGPTOTG), ZGM(YDGEOMETRY%YRGEM%NGPTOTG)

INTEGER(KIND=JPIM) :: IEND, IIEND, IIST, INDICE, INSIGN, IBLK, IST
INTEGER(KIND=JPIM) :: JFLD, JJ, JKGLO, JROF

LOGICAL :: LLHMT

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "acsolw.intfb.h"
#include "elislap.intfb.h"
#include "lislap.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CASMSWI',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDGEM=>YDGEOMETRY%YRGEM, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDCSGLEG=>YDGEOMETRY%YRCSGLEG&
& )
ASSOCIATE(GCONV=>YDPHY1%GCONV,   NGPTOTG=>YDGEM%NGPTOTG, NLOENG=>YDGEM%NLOENG, NGPTOT=>YDGEM%NGPTOT,               &
& RSTRET=>YDGEM%RSTRET,   NDGLG=>YDDIM%NDGLG, NDLON=>YDDIM%NDLON, NDGNH=>YDDIM%NDGNH,   NPROMA=>YDDIM%NPROMA,      &
& YSP_SG=>YDSURF%YSP_SG, YSD_VV=>YDSURF%YSD_VV, YSP_SB=>YDSURF%YSP_SB,   SD_VF=>YDSURF%SD_VF, SP_RR=>YDSURF%SP_RR, &
& YSP_RR=>YDSURF%YSP_RR,   YSD_VF=>YDSURF%YSD_VF, SD_VV=>YDSURF%SD_VV, SP_SB=>YDSURF%SP_SB,   SP_SG=>YDSURF%SP_SG  &
& )
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!*       1.    READ BUFFER AND COMPUTATION OF SWI AND MASK FOR SMOOTHING
!              ---------------------------------------------------------

WRITE (NULOUT,*) ' #######  SUBROUTINE CASMSWI  ####### '
ZMASK=0.5_JPRB
ZEPW = 1.E-3_JPRB
ZOCE(:)=1.0_JPRB
LLHMT=.TRUE.
IST=1
DO JKGLO=1,NGPTOT,NPROMA
  IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
  IBLK=(JKGLO-1)/NPROMA+1
  DO JROF=IST,IEND
!   Land/sea mask
    ZLSM=SD_VF(JROF,YSD_VF%YLSM%MP  ,IBLK)
!   Snow water content
    ZSD=SP_SG(JROF,1,YSP_SG%YF%MP0 ,IBLK) ! Assumes only 1 snow layer active
!   Total soil frozen water content
    ZFDSW=SP_SB(JROF,1,YSP_SB%YTL%MP0,IBLK)
!   Fraction of vegetation
    ZVEG(JROF+JKGLO-1)=SD_VF(JROF,YSD_VF%YVEG%MP,IBLK)
!   Soil depth
    ZD2(JROF+JKGLO-1)=SD_VV(JROF,YSD_VV%YD2 %MP,IBLK)
!   Superficial soil liquid water content
    ZSSW(JROF+JKGLO-1)=SP_RR(JROF,YSP_RR%YW%MP0 ,IBLK)
!   Total soil liquid water content
    ZDSW(JROF+JKGLO-1)=SP_SB(JROF,1,YSP_SB%YQ%MP0,IBLK)
!   New mask, if it is land, without snow and no frozen soil
!   water content, input for (e)lislap is set to _ZERO_
    IF ((ZLSM == 1.0_JPRB) .AND. (ZSD == 0.0_JPRB) .AND.&
     & (ZFDSW == 0.0_JPRB)) ZOCE(JROF+JKGLO-1)=0.0_JPRB
  ENDDO

! Calculation of soil water contents at wilting point and field capacity
  IIST=IST+JKGLO-1
  IIEND=IEND+JKGLO-1
  CALL ACSOLW(YDPHY1,IST,IEND,IEND,SD_VV(IST:IEND,YSD_VV%YARG%MP,IBLK),&
   & SD_VV(IST:IEND,YSD_VV%YD2%MP,IBLK),SD_VF(IST:IEND,YSD_VF%YLSM%MP,IBLK),&
   & SD_VV(IST:IEND,YSD_VV%YIVEG%MP,IBLK),SD_VV(IST:IEND,YSD_VV%YSAB%MP,IBLK),&
   & LLHMT,ZWFC(IIST:IIEND),ZWPMX(IIST:IIEND),&
   & ZWSAT(IIST:IIEND),ZWSMX(IIST:IIEND),ZWWILT(IIST:IIEND))
  ZWFC(IIST:IIEND)=ZWFC(IIST:IIEND)*GCONV*SD_VV(IST:IEND,YSD_VV%YD2%MP,IBLK)
  ZWWILT(IIST:IIEND)=ZWWILT(IIST:IIEND)*GCONV*SD_VV(IST:IEND,YSD_VV%YD2%MP,IBLK)
ENDDO

! Calculation of soil wetness index (SWI)
ZSWI(1:NGPTOT)=(ZDSW(1:NGPTOT)-ZWWILT(1:NGPTOT))/(ZWFC(1:NGPTOT)&
 & -ZWWILT(1:NGPTOT))

!     ------------------------------------------------------------------
!*       2.    DISTRIBUTE AND COLLECT GRID-POINT FIELDS
!              ----------------------------------------

IF (NPROC > 1) THEN
  IF (MYPROC /= 1) THEN
    JFLD=1
    CALL DIWRGRID_SEND(YDGEOMETRY%YRGEM,1,1,ZSWI(1:NGPTOT),JFLD)
    JFLD=JFLD+1
    CALL DIWRGRID_SEND(YDGEOMETRY%YRGEM,1,1,ZOCE(1:NGPTOT),JFLD)
    JFLD=JFLD+1
    CALL DIWRGRID_SEND(YDGEOMETRY%YRGEM,1,1,YDGSGEOM_NB%GM(1:NGPTOT),JFLD)
  ELSE
    JFLD=1
    CALL DIWRGRID_RECV(YDGEOMETRY,1,ZSWI(1:NGPTOT),JFLD,ZGSWI)
    JFLD=JFLD+1
    CALL DIWRGRID_RECV(YDGEOMETRY,1,ZOCE(1:NGPTOT),JFLD,ZGOCE)
    JFLD=JFLD+1
    CALL DIWRGRID_RECV(YDGEOMETRY,1,YDGSGEOM_NB%GM(1:NGPTOT),JFLD,ZGM)
  ENDIF
ELSE
  ZGSWI(1:NGPTOT)=ZSWI(1:NGPTOT)
  ZGOCE(1:NGPTOT)=ZOCE(1:NGPTOT)
  ZGM(1:NGPTOT)=YDGSGEOM_NB%GM(1:NGPTOT)
ENDIF

!     ------------------------------------------------------------------
!*       3.    SMOOTHING ON PE 1
!              -----------------

IF (MYPROC == 1) THEN
  IF (.NOT.LELAM) THEN
    DO JJ=1,NDGLG
      INDICE=MIN(JJ,NDGLG+1-JJ)
      INSIGN=SIGN(1,NDGNH-JJ)
      ZBO(JJ)=REAL(INSIGN,JPRB)*ASIN(YDCSGLEG%RMU(INDICE))
    ENDDO
    DO JJ=1,NR_SM_WP
      CALL LISLAP(ZGSWI(1:NGPTOTG),NLOENG(1:NDGLG),NDGLG,NDLON,NGPTOTG,&
       & YDCSGLEG%RSQM2(1:NDGLG),ZBO(1:NDGLG),RSTRET,ZGOCE(1:NGPTOTG),ZMASK,RA_SM_WP)
    ENDDO
  ELSE
    DO JJ=1,NR_SM_WP
      CALL ELISLAP(YDGEOMETRY%YREGEO,ZGSWI(1:NGPTOTG),NDGLG,NDLON,NGPTOTG,&
       & ZGOCE(1:NGPTOTG),ZGM(1:NGPTOTG),ZMASK,RA_SM_WP)
    ENDDO
  ENDIF
ENDIF

!     ------------------------------------------------------------------
!*       4.    DISTRIBUTE AND COLLECT GRID-POINT FIELD
!              ---------------------------------------

IF (NPROC > 1) THEN
  JFLD=1
  IF (MYPROC == 1) THEN
    CALL DISGRID_SEND(YDGEOMETRY,1,ZGSWI,JFLD,ZSWI)
  ELSE
    CALL DISGRID_RECV(YDGEOMETRY,1,1,ZSWI,JFLD)
  ENDIF
ELSE
  ZSWI(1:NGPTOT)=ZGSWI(1:NGPTOT)
ENDIF

!     ------------------------------------------------------------------
!*       5.    MODIFICATION OF TOTAL SOIL LIQUID WATER CONTENT
!              -----------------------------------------------

DO JJ=1,NGPTOT
  IF (ZOCE(JJ) == 0.0_JPRB) THEN
!   changes for points over the land, no frozen deep soil wetness
!   and no snow
    ZDSWN(JJ)=ZSWI(JJ)*ZWFC(JJ)+(1.0_JPRB-ZSWI(JJ))*ZWWILT(JJ)
!   Test that new field is inside of limits like it is in ISBA
    IF (ZDSWN(JJ) <= ZVEG(JJ)*ZWWILT(JJ)) THEN
      ZDSW(JJ)=MAX(ZDSW(JJ),ZDSWN(JJ))
    ELSEIF (ZDSWN(JJ) >= ZWFC(JJ)) THEN
      ZDSW(JJ)=MIN(ZDSW(JJ),ZDSWN(JJ))
    ELSE
      ZDSW(JJ)=ZDSWN(JJ)
    ENDIF
    ZDSW(JJ)=MAX(ZDSW(JJ),MAX(ZSSW(JJ),ZEPW*ZD2(JJ)*GCONV))
  ENDIF
ENDDO

! Put the changed field ZDSW in GPPBUF on each PE
IST=1
DO JKGLO=1,NGPTOT,NPROMA
  IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
  IBLK=(JKGLO-1)/NPROMA+1
  DO JROF=IST,IEND
    SP_SB(JROF,1,YSP_SB%YQ%MP0,IBLK)=ZDSW(JROF+JKGLO-1)
  ENDDO
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CASMSWI',1,ZHOOK_HANDLE)
END SUBROUTINE CASMSWI

