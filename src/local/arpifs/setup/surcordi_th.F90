!OCL NOEVAL
SUBROUTINE SURCORDI_TH(YDGEOMETRY,YDDYN)

!**** *SURCORDI_TH*   - Computation of RCORDIT, RCORDIH, RCORDIF

!     Purpose.
!     --------
!      The main difference from the original alpha formulation is
!      that this one makes T-alpha*ln(ps) behaving more like potential
!      temperature. 

!**   Interface.
!     ----------
!        *CALL* *SURCORDI_TH

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE/ALADIN DOCUMENTATION

!     Author.
!     -------
!        K. Vana , written after original SURCORDI

!     Modifications.
!     --------------
!        Original : 5-Jun-2013
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMMP0   , ONLY : NPRINTLEV
USE YOMLUN   , ONLY : NULOUT
USE YOMDYN   , ONLY : TDYN
USE YOMCST   , ONLY : RD       ,RKAPPA
USE YOMSTA   , ONLY : RTSUR
USE YOMVERT  , ONLY : VP00

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT) :: YDDYN
INTEGER(KIND=JPIM) :: JLEV

REAL(KIND=JPRB) :: ZTBH
REAL(KIND=JPRB) :: ZPPH(0:YDGEOMETRY%YRDIMV%NFLEVG), ZPP(YDGEOMETRY%YRDIMV%NFLEVG), ZPHH(0:YDGEOMETRY%YRDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURCORDI_TH',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDSTA=>YDGEOMETRY%YRSTA)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG,   RCORDIF=>YDDYN%RCORDIF,   RCORDIH=>YDDYN%RCORDIH, RCORDIT=>YDDYN%RCORDIT,    &
& STPHI=>YDSTA%STPHI, STTEM=>YDSTA%STTEM)
!     ------------------------------------------------------------------

!*       1.    Computation of RCORDIT, RCORDIH, RCORDIF:

! * Remark: this piece of code is fully consistent with
!   lvertfe=F, ndlnpr=0, laprxpk=T
!   where prehyd(l)=0.5*prehyd(lbar-1)+0.5*prehyd(lbar)
!   and not completely consistent with the other cases
!   for lvertfe,ndlnpr,laprxpk.

IF (NFLEVG > 1) THEN
  RCORDIH(0)=0._JPRB
  DO JLEV=NFLEVG,1,-1
    RCORDIT(JLEV)=RKAPPA*(YDVAB%VBH(JLEV-1)+YDVAB%VBH(JLEV))*0.5_JPRB*STTEM(JLEV)*EXP(STPHI(JLEV)/(RD*RTSUR))
    RCORDIF(JLEV)=RCORDIT(JLEV)
    ZPP (JLEV)=0.5_JPRB*(YDVAB%VAH(JLEV)+YDVAB%VAH(JLEV-1)+(YDVAB%VBH(JLEV)+YDVAB%VBH(JLEV-1))*VP00)
    ZPPH(JLEV)=          YDVAB%VAH(JLEV)                  + YDVAB%VBH(JLEV)*VP00
    IF(JLEV==NFLEVG) THEN
      ZPHH(JLEV)=0._JPRB
      ZTBH=RTSUR
    ELSE
      ZPHH(JLEV)=ZPHH(JLEV+1)-RD*STTEM(JLEV+1)/ZPP(JLEV+1)*(ZPPH(JLEV)-ZPPH(JLEV+1))
      ZTBH=STTEM(JLEV+1)+(STTEM(JLEV)-STTEM(JLEV+1))/(STPHI(JLEV)-STPHI(JLEV+1))*(ZPHH(JLEV)-STPHI(JLEV+1))
    ENDIF
    RCORDIH(JLEV)=RKAPPA*YDVAB%VBH(JLEV)*ZTBH*EXP(ZPHH(JLEV)/(RD*RTSUR))
  ENDDO
ELSE
  RCORDIT(1)=1._JPRB
ENDIF

! * Printings:
                                                                                
IF (NPRINTLEV > 0) THEN
  WRITE(NULOUT,'(A)')
  WRITE(NULOUT,'(A)') ' --- SURCORDI_TH:'
  WRITE(NULOUT,'(A)') ' RCORDIT:'
  WRITE(NULOUT,'(4(1X,E19.13))') (RCORDIT(JLEV), JLEV=1,NFLEVG)
  IF (NFLEVG > 1) THEN
    WRITE(NULOUT,'(A)') ' RCORDIH:'
    WRITE(NULOUT,'(4(1X,E19.13))') (RCORDIH(JLEV), JLEV=0,NFLEVG)
    WRITE(NULOUT,'(A)') ' RCORDIF:'
    WRITE(NULOUT,'(4(1X,E19.13))') (RCORDIF(JLEV), JLEV=1,NFLEVG)
  ENDIF
  CALL FLUSH(NULOUT)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURCORDI_TH',1,ZHOOK_HANDLE)
END SUBROUTINE SURCORDI_TH
