SUBROUTINE SUNH_VERTFE3D(YDGEOMETRY)

!**** *SUNH_VERTFE3D*  - Setup VERTical Finite Element Derivative based on
!                        cubic (3) elements: first order derivative.
!                        Version designed for the NH model.

!     Purpose.
!     --------
!           Compute matrix for the vertical differentiation based on
!           cubic B-spline finite elements

!**   Interface.
!     ----------

!     *CALL* SUNH_VERTFE3D

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!        MINV

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Jozef Vivoda (after SUVERTFE3D) SHMI/LACE
!        Original : 2006-06

!     Modifications.
!     --------------
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   K. Yessad (Sep 2008): use VFE_ETAH and VFE_ETAF
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMMP0   , ONLY : NPRINTLEV, MYPROC
USE YOMLUN   , ONLY : NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
INTEGER(KIND=JPIM) :: IDIM, J, JK, IKC, IJD, JR, IKR, J1, J2, J3, JC, II
INTEGER(KIND=JPIM) :: IROW
REAL(KIND=JPRB)    :: ZDELTA(-3:YDGEOMETRY%YRDIMV%NFLEVG+3)
REAL(KIND=JPRB)    :: ZFAA(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZFBB(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZFCC(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZFDD(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZDAA(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZDBB(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZDCC(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZDDD(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZWAA(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZWBB(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZWCC(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZWDD(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZA(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZAINV(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZWORK(2*(YDGEOMETRY%YRDIMV%NFLEVG+3))
REAL(KIND=JPRB)    :: ZPROD(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZB(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZDERI(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZSF(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZSFINV(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZSD(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZEPS, ZALPHA(4,4), ZBETA(2,2), ZHH(4,4)
REAL(KIND=JPRB)    :: ZRHS(2)
REAL(KIND=JPRB)    :: ZBMAX, ZD1, ZD2, ZD3, ZD4, ZDETB, ZDEL, ZDEL2, ZDEL3
REAL(KIND=JPRB)    :: ZDDD1, ZDDD2, ZDDD3, ZDDD4, ZDET
REAL(KIND=JPRB)    :: ZAAA2, ZBBB2, ZCCC2, ZAAA3, ZBBB3, ZCCC3
REAL(KIND=JPRB)    :: FHH, FDD, FPOLY, FDERI
REAL(KIND=JPRB)    :: ZARG, ZW2T, ZW1T, ZSUM
! REAL(KIND=JPRB)    :: ZCONST
REAL(KIND=JPRB)    :: ZW2B, ZW1B
REAL(KIND=JPRB)    :: XA,XB,XC,XD,YA,YB,YC,YD, XX

LOGICAL   :: LL_CHANGE_BASE_V1       ! change of boundary conditions
LOGICAL   :: LLGP_AT01               ! grid points at eta=0 and eta=1

REAL(KIND=JPRB)    :: ZHOOK_HANDLE

#include "minv.h"

!     ------------------------------------------------------------------

! indefiniteINTEGRAL(EJ*EI)dx
FHH(XX,XA,XB,XC,XD,YA,YB,YC,YD)=&
 &   XX*(  XA*YA&
 &   + XX*( (XA*YB+XB*YA)/2._JPRB&
 &   + XX*( (XA*YC+XB*YB+XC*YA)/3._JPRB&
 &   + XX*( (XA*YD+XB*YC+XC*YB+XD*YA)/4._JPRB&
 &   + XX*( (      XB*YD+XC*YC+XD*YB)/5._JPRB&
 &   + XX*( (            XC*YD+XD*YC)/6._JPRB&
 & + XX*                     XD*YD /7._JPRB))))))  

! indefiniteINTEGRAL(EJ*EI')dx                '
FDD(XX,XA,XB,XC,XD,YA,YB,YC,YD)=&
 &   XX*(     XA*YB&
 &   + XX*( (2._JPRB*XA*YC+XB*YB         )/2._JPRB&
 &   + XX*( (3._JPRB*XA*YD+2._JPRB*XB*YC+XC*YB)/3._JPRB&
 &   + XX*( (3._JPRB*XB*YD+2._JPRB*XC*YC+XD*YB)/4._JPRB&
 &   + XX*( (3._JPRB*XC*YD+2._JPRB*XD*YC      )/5._JPRB&
 & + XX*( (3._JPRB*XD*YD               )/6._JPRB))))))

FPOLY(XX,XA,XB,XC,XD)=XA+XX*(XB+XX*(XC+XD*XX))
FDERI(XX,XA,XB,XC,XD)=XB+XX*(2._JPRB*XC+3._JPRB*XD*XX)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUNH_VERTFE3D',0,ZHOOK_HANDLE)
ASSOCIATE(YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)!     ------------------------------------------------------------------

IDIM=NFLEVG+4
ZEPS=100.0_JPRB*TINY(1.0_JPRB)

LL_CHANGE_BASE_V1 = .FALSE.
LLGP_AT01 = .TRUE.

! compute distance between nodes (full eta levels)

IF( LLGP_AT01 )THEN
  ZDELTA(0)=(YDVETA%VFE_ETAF(1)-YDVETA%VFE_ETAH(0))
ELSE
  ZDELTA(0)=2._JPRB*(YDVETA%VFE_ETAF(1)-YDVETA%VFE_ETAH(0))
ENDIF

ZDELTA(-3)=ZDELTA(0)
ZDELTA(-2)=ZDELTA(0)
ZDELTA(-1)=ZDELTA(0)
DO JK=1,NFLEVG-1
  ZDELTA(JK)=YDVETA%VFE_ETAF(JK+1)-YDVETA%VFE_ETAF(JK)
ENDDO
IF( LLGP_AT01 )THEN
  ZDELTA(NFLEVG)=(YDVETA%VFE_ETAH(NFLEVG) -YDVETA%VFE_ETAF(NFLEVG))
ELSE
  ZDELTA(NFLEVG)=2._JPRB*(YDVETA%VFE_ETAH(NFLEVG) -YDVETA%VFE_ETAF(NFLEVG))
ENDIF
ZDELTA(NFLEVG+1)=ZDELTA(NFLEVG)
ZDELTA(NFLEVG+2)=ZDELTA(NFLEVG)
ZDELTA(NFLEVG+3)=ZDELTA(NFLEVG)

! compute basis functions    

ZBMAX=4._JPRB/6._JPRB 
  
DO JK=-1,NFLEVG+2
  ZD1=ZDELTA(JK-2)
  ZD2=ZDELTA(JK-1)
  ZD3=ZDELTA(JK  )
  ZD4=ZDELTA(JK+1)

  ! matrix alpha
  ZALPHA(1,1)=0.0_JPRB     
  ZALPHA(1,2)=0.0_JPRB      
  ZALPHA(1,3)=ZD3**3      
  ZALPHA(1,4)=ZD4**3+3._JPRB*ZD3*ZD4**2+3._JPRB*ZD3**2*ZD4  
  ZALPHA(2,1)=ZD1**2+2._JPRB*ZD1*ZD2
  ZALPHA(2,2)=ZD2**2
  ZALPHA(2,3)=ZD3**2
  ZALPHA(2,4)=ZD4**2+2._JPRB*ZD3*ZD4
  ZALPHA(3,1)= ZD1
  ZALPHA(3,2)= ZD2
  ZALPHA(3,3)=-ZD3
  ZALPHA(3,4)=-ZD4
  ZALPHA(4,1)=ZD1**3+3._JPRB*ZD1**2*ZD2+3._JPRB*ZD1*ZD2**2
  ZALPHA(4,2)=ZD2**3
  ZALPHA(4,3)=0.0_JPRB
  ZALPHA(4,4)=0.0_JPRB

  ! matrix beta
  ZBETA(1,1)=ZALPHA(2,2)-ZALPHA(2,1)*ZALPHA(4,2)/ZALPHA(4,1)
  ZBETA(1,2)=ZALPHA(2,3)-ZALPHA(2,4)*ZALPHA(1,3)/ZALPHA(1,4)
  ZBETA(2,1)=ZALPHA(3,2)-ZALPHA(3,1)*ZALPHA(4,2)/ZALPHA(4,1)
  ZBETA(2,2)=ZALPHA(3,3)-ZALPHA(3,4)*ZALPHA(1,3)/ZALPHA(1,4)

  ! right-hand side vector
  ZRHS(1)=-ZBMAX*(ZALPHA(2,1)/ZALPHA(4,1)+ZALPHA(2,4)/ZALPHA(1,4))
  ZRHS(2)=-ZBMAX*(ZALPHA(3,1)/ZALPHA(4,1)+ZALPHA(3,4)/ZALPHA(1,4))

  ! determinant of beta 
  ZDETB=ZBETA(1,1)*ZBETA(2,2)-ZBETA(1,2)*ZBETA(2,1)

  ! coefficients 
  ZDDD2=(ZRHS(1)*ZBETA(2,2)-ZRHS(2)*ZBETA(1,2))/ZDETB
  ZDDD3=(ZRHS(2)*ZBETA(1,1)-ZRHS(1)*ZBETA(2,1))/ZDETB
  ZDDD1=(ZBMAX - ZALPHA(4,2)*ZDDD2)/ZALPHA(4,1)    
  ZDDD4=(ZBMAX - ZALPHA(1,3)*ZDDD3)/ZALPHA(1,4) 

  ZAAA2=ZD1**3*ZDDD1
  ZBBB2=3._JPRB*ZD1**2*ZDDD1
  ZCCC2=3._JPRB*ZD1   *ZDDD1

  ZAAA3=ZD4**3*ZDDD4
  ZBBB3=3._JPRB*ZD4**2*ZDDD4
  ZCCC3=3._JPRB*ZD4   *ZDDD4

  ! coefficients of explicit cubic polynomials in each interval
  ! (x=0 at the left of the corresponding interval)
  
  ZFAA(1,JK)=0.0_JPRB
  ZFBB(1,JK)=0.0_JPRB
  ZFCC(1,JK)=0.0_JPRB
  ZFDD(1,JK)=ZDDD1
 
  ZFAA(2,JK)=ZAAA2
  ZFBB(2,JK)=ZBBB2
  ZFCC(2,JK)=ZCCC2
  ZFDD(2,JK)=ZDDD2  
  
  ZDEL=ZDELTA(JK)
  ZDEL2=ZDEL*ZDEL
  ZDEL3=ZDEL2*ZDEL
  ZFAA(3,JK)=ZAAA3+ZBBB3*ZDEL+ZCCC3*ZDEL2+ZDDD3*ZDEL3
  ZFBB(3,JK)=-ZBBB3-2._JPRB*ZCCC3*ZDEL-3._JPRB*ZDDD3*ZDEL2
  ZFCC(3,JK)=ZCCC3+3._JPRB*ZDDD3*ZDEL
  ZFDD(3,JK)=-ZDDD3

  ZDEL=ZDELTA(JK+1)
  ZDEL2=ZDEL*ZDEL
  ZDEL3=ZDEL2*ZDEL
  ZFAA(4,JK)=ZDDD4*ZDEL3
  ZFBB(4,JK)=-3._JPRB*ZDDD4*ZDEL2
  ZFCC(4,JK)= 3._JPRB*ZDDD4*ZDEL
  ZFDD(4,JK)=-ZDDD4 
ENDDO

! derivative basis
ZDAA(:,:)=ZFAA(:,:)
ZDBB(:,:)=ZFBB(:,:)
ZDCC(:,:)=ZFCC(:,:)
ZDDD(:,:)=ZFDD(:,:)

! weighting functions
ZWAA(:,:)=ZFAA(:,:)
ZWBB(:,:)=ZFBB(:,:)
ZWCC(:,:)=ZFCC(:,:)
ZWDD(:,:)=ZFDD(:,:)

IF( LL_CHANGE_BASE_V1 )THEN

  !----------------------------
  ! change of weighting function
  !----------------------------
  DO J=1,4
    ZWAA(J,-1)=0.0_JPRB
    ZWBB(J,-1)=0.0_JPRB
    ZWCC(J,-1)=0.0_JPRB
    ZWDD(J,-1)=0.0_JPRB
  ENDDO
  DO J=1,2
    ZWAA(J,0)=0.0_JPRB
    ZWBB(J,0)=0.0_JPRB
    ZWCC(J,0)=0.0_JPRB
    ZWDD(J,0)=0.0_JPRB
  ENDDO

  ZWAA(1,1)=0.0_JPRB
  ZWBB(1,1)=0.0_JPRB
  ZWCC(1,1)=0.0_JPRB
  ZWDD(1,1)=0.0_JPRB

  DO J=2,4
    ZWAA(J,NFLEVG+2)=0.0_JPRB
    ZWBB(J,NFLEVG+2)=0.0_JPRB
    ZWCC(J,NFLEVG+2)=0.0_JPRB
    ZWDD(J,NFLEVG+2)=0.0_JPRB
  ENDDO
  DO J=3,4
    ZWAA(J,NFLEVG+1)=0.0_JPRB
    ZWBB(J,NFLEVG+1)=0.0_JPRB
    ZWCC(J,NFLEVG+1)=0.0_JPRB
    ZWDD(J,NFLEVG+1)=0.0_JPRB
  ENDDO
  ZWAA(4,NFLEVG)=0.0_JPRB
  ZWBB(4,NFLEVG)=0.0_JPRB
  ZWCC(4,NFLEVG)=0.0_JPRB
  ZWDD(4,NFLEVG)=0.0_JPRB

ENDIF

! COMPUTE MATRIX SF ( SF(I,J)=BASISFUNC_J(X_I) )
!  (from FE space of the function), to be differenciated into physical space
! ---------------------------------------------

ZSF(:,:)=0.0_JPRB
DO JR=1,IDIM
  DO JC=1,IDIM
    IKC=JC-2
    IF(JC == JR-1) ZSF(JR,JC)=ZFAA(4,IKC)
    IF(JC == JR  ) ZSF(JR,JC)=ZFAA(3,IKC)
    IF(JC == JR+1) ZSF(JR,JC)=ZFAA(2,IKC)
  ENDDO
ENDDO

! Instead of the value at node -1 and NFLEVG+2 we should use the derivative
!   at node 1 and at NFLEVG:

!---------
! TBC - derivative at the model top prescribed
!---------
! index where the value is prescribed (if it is in nodal point of course)
! IROW - the index in which derivative is in input vector
! II   - the nodal point index on which derivative is valid
!      - if II=0 and ARG=ZDELTA(0)/2 then the values are evaluated at model top
!        where eta=0

IROW = 1
II   = 0
IF( LLGP_AT01 )THEN
  ZARG = 0.0_JPRB
ELSE
  ZARG = 0.5_JPRB*ZDELTA(0)
ENDIF
ZSF(IROW,:)=0.0_JPRB
ZSF(IROW,1)=FDERI(ZARG,ZFAA(4, II-1),ZFBB(4, II-1),ZFCC(4, II-1),ZFDD(4, II-1))
ZSF(IROW,2)=FDERI(ZARG,ZFAA(3, II  ),ZFBB(3, II  ),ZFCC(3, II  ),ZFDD(3, II  ))
ZSF(IROW,3)=FDERI(ZARG,ZFAA(2, II+1),ZFBB(2, II+1),ZFCC(2, II+1),ZFDD(2, II+1))
ZSF(IROW,4)=FDERI(ZARG,ZFAA(1, II+2),ZFBB(1, II+2),ZFCC(1, II+2),ZFDD(1, II+2))

!---------
! TBC - value at the model top prescribed
!---------
IROW = 2
II   = 0
IF( LLGP_AT01 )THEN
  ZARG = 0.0_JPRB
ELSE
  ZARG = 0.5_JPRB*ZDELTA(0)
ENDIF
ZSF(IROW,:)=0.0_JPRB
ZSF(IROW,1)=FPOLY(ZARG,ZFAA(4, II-1),ZFBB(4, II-1),ZFCC(4, II-1),ZFDD(4, II-1))
ZSF(IROW,2)=FPOLY(ZARG,ZFAA(3, II  ),ZFBB(3, II  ),ZFCC(3, II  ),ZFDD(3, II  ))
ZSF(IROW,3)=FPOLY(ZARG,ZFAA(2, II+1),ZFBB(2, II+1),ZFCC(2, II+1),ZFDD(2, II+1))
ZSF(IROW,4)=FPOLY(ZARG,ZFAA(1, II+2),ZFBB(1, II+2),ZFCC(1, II+2),ZFDD(1, II+2))

!---------
! BBC - value at the model bottom prescribed
!---------
IROW  = IDIM-1
II    = NFLEVG
IF( LLGP_AT01 )THEN
  ZARG = ZDELTA(NFLEVG)
ELSE
  ZARG = 0.5_JPRB*ZDELTA(NFLEVG)
ENDIF
ZSF(IROW,:)=0.0_JPRB
ZSF(IROW,IDIM-3)=FPOLY(ZARG,ZFAA(4,II-1),ZFBB(4,II-1),ZFCC(4,II-1),ZFDD(4,II-1))
ZSF(IROW,IDIM-2)=FPOLY(ZARG,ZFAA(3,II  ),ZFBB(3,II  ),ZFCC(3,II  ),ZFDD(3,II  ))
ZSF(IROW,IDIM-1)=FPOLY(ZARG,ZFAA(2,II+1),ZFBB(2,II+1),ZFCC(2,II+1),ZFDD(2,II+1))
ZSF(IROW,IDIM  )=FPOLY(ZARG,ZFAA(1,II+2),ZFBB(1,II+2),ZFCC(1,II+2),ZFDD(1,II+2))

!---------
! BBC - derivative at the model bottom prescribed
!---------
IROW  = IDIM
II=NFLEVG
IF( LLGP_AT01 )THEN
  ZARG = ZDELTA(NFLEVG)
ELSE
  ZARG = 0.5_JPRB*ZDELTA(NFLEVG)
ENDIF
ZSF(IROW,:)=0.0_JPRB
ZSF(IROW,IDIM-3)=FDERI(ZARG,ZFAA(4,II-1),ZFBB(4,II-1),ZFCC(4,II-1),ZFDD(4,II-1))
ZSF(IROW,IDIM-2)=FDERI(ZARG,ZFAA(3,II  ),ZFBB(3,II  ),ZFCC(3,II  ),ZFDD(3,II  ))
ZSF(IROW,IDIM-1)=FDERI(ZARG,ZFAA(2,II+1),ZFBB(2,II+1),ZFCC(2,II+1),ZFDD(2,II+1))
ZSF(IROW,IDIM  )=FDERI(ZARG,ZFAA(1,II+2),ZFBB(1,II+2),ZFCC(1,II+2),ZFDD(1,II+2))

! invert SF:       (ZSFINV projects the function to be integrated
!                  from physical space onto FE space)
ZDET=0.0_JPRB
ZSFINV(:,:)=ZSF(:,:)*3._JPRB/2._JPRB
CALL MINV(ZSFINV,IDIM,IDIM,ZWORK,ZDET,ZEPS,0,1)
ZSFINV(:,:)=ZSFINV(:,:)*3._JPRB/2._JPRB

! COMPUTE MATRIX SD ( SD(I,J)=BASISDER_J(X_I) )
!  To project the derivative from FE space to physical space 
! -------------------------------------------

ZSD(:,:)=0.0_JPRB
DO JR=1,IDIM
  DO JC=1,IDIM
    IKC=JC-2
    IF(JC == JR-1) ZSD(JR,JC)=ZDAA(4,IKC)
    IF(JC == JR  ) ZSD(JR,JC)=ZDAA(3,IKC)
    IF(JC == JR+1) ZSD(JR,JC)=ZDAA(2,IKC)
  ENDDO
ENDDO

!----------------------------------
! in output I want the derivatives:
! i) at the model bottom (item IDIM)
! i) at the model top    (item 1)
!----------------------------------
!---------
! BBC - value at the model bottom (eta=1)
!---------
II = NFLEVG
IF( LLGP_AT01 )THEN
  ZARG = ZDELTA(NFLEVG)
ELSE
  ZARG = 0.5_JPRB*ZDELTA(NFLEVG)
ENDIF
ZSD(IDIM,:)=0.0_JPRB
ZSD(IDIM,IDIM-3)=FPOLY(ZARG,ZDAA(4,II-1),ZDBB(4,II-1),ZDCC(4,II-1),ZDDD(4,II-1))
ZSD(IDIM,IDIM-2)=FPOLY(ZARG,ZDAA(3,II  ),ZDBB(3,II  ),ZDCC(3,II  ),ZDDD(3,II  ))
ZSD(IDIM,IDIM-1)=FPOLY(ZARG,ZDAA(2,II+1),ZDBB(2,II+1),ZDCC(2,II+1),ZDDD(2,II+1))
ZSD(IDIM,IDIM  )=FPOLY(ZARG,ZDAA(1,II+2),ZDBB(1,II+2),ZDCC(1,II+2),ZDDD(1,II+2))

!---------
! TBC - value at the model top (eta=0)
!---------
II = 0
IF( LLGP_AT01 )THEN
  ZARG = 0.0_JPRB
ELSE
  ZARG = 0.5_JPRB*ZDELTA(0)
ENDIF
ZSD(1,:)=0.0_JPRB
ZSD(1,1)=FPOLY(ZARG,ZDAA(4, II-1),ZDBB(4, II-1),ZDCC(4, II-1),ZDDD(4, II-1))
ZSD(1,2)=FPOLY(ZARG,ZDAA(3, II  ),ZDBB(3, II  ),ZDCC(3, II  ),ZDDD(3, II  ))
ZSD(1,3)=FPOLY(ZARG,ZDAA(2, II+1),ZDBB(2, II+1),ZDCC(2, II+1),ZDDD(2, II+1))
ZSD(1,4)=FPOLY(ZARG,ZDAA(1, II+2),ZDBB(1, II+2),ZDCC(1, II+2),ZDDD(1, II+2))

! COMPUTE MATRIX A
! ----------------

ZA(:,:)=0.0_JPRB

! JR is the row number of matrix A
DO JR=1,IDIM
  ! IKR is the index of the test function
  IKR=JR-2
  ! ZHH is an auxiliary 4x4 matrix
  ZHH(:,:)=0.0_JPRB
  ! J is the index of the interval belonging to function IKR
  DO J=1,4
    ! IJD is the absolute number of interval J
    IJD=IKR-3+J
    ! Integration is from node 0 to node NFLEVG+1
    IF(IJD >= 0 .AND. IJD <= NFLEVG) THEN
      DO JK=1,4
        ! IKC is the JKth basis number covering interval J
        IKC=IKR+J-JK
        IF(IKC >= -1 .AND. IKC <= NFLEVG+2) THEN
          ZHH(J,JK)=FHH(ZDELTA(IJD) &
           & ,ZWAA(J,IKR),ZWBB(J,IKR),ZWCC(J,IKR),ZWDD(J,IKR) &
           & ,ZFAA(JK,IKC),ZFBB(JK,IKC),ZFCC(JK,IKC),ZFDD(JK,IKC))  
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  ZA(JR,JR  )=ZHH(1,1)+ZHH(2,2)+ZHH(3,3)+ZHH(4,4)
  IF(JR+1 <= IDIM) ZA(JR,JR+1)=ZHH(2,1)+ZHH(3,2)+ZHH(4,3)  
  IF(JR+2 <= IDIM) ZA(JR,JR+2)=ZHH(3,1)+ZHH(4,2)  
  IF(JR+3 <= IDIM) ZA(JR,JR+3)=ZHH(4,1)
  IF(JR-1 >= 1)    ZA(JR,JR-1)=ZHH(1,2)+ZHH(2,3)+ZHH(3,4)
  IF(JR-2 >= 1)    ZA(JR,JR-2)=ZHH(1,3)+ZHH(2,4)
  IF(JR-3 >= 1)    ZA(JR,JR-3)=ZHH(1,4)     
ENDDO

! invert A: 
ZDET=0.0_JPRB
ZAINV(:,:)=ZA(:,:)*REAL(NFLEVG,JPRB)*7._JPRB/2._JPRB
CALL MINV(ZAINV,IDIM,IDIM,ZWORK,ZDET,ZEPS,0,1)
ZAINV(:,:)=ZAINV(:,:)*REAL(NFLEVG,JPRB)*7._JPRB/2._JPRB

DO J1=1,IDIM
  DO J2=1,IDIM
    ZPROD(J1,J2)=0.0_JPRB
    DO J3=1,IDIM
      ZPROD(J1,J2)=ZPROD(J1,J2)+ZAINV(J1,J3)*ZA(J3,J2)
    ENDDO
  ENDDO
ENDDO

! COMPUTE MATRIX B
! ----------------

ZB(:,:)=0.0_JPRB

! JR is the row number of matrix B
DO JR=1,IDIM
  ! IKR is the index of the test function
  IKR=JR-2
  ! ZHH is an auxiliary 4x4 matrix
  ZHH(:,:)=0.0_JPRB
  ! J is the index of the interval belonging to function IKR
  DO J=1,4
    ! IJD is the absolute number of interval J
    IJD=IKR-3+J
    ! Integration is from node 0 to node NFLEVG+1
    IF(IJD >= 0 .AND. IJD <= NFLEVG) THEN
      DO JK=1,4
        ! IKC is the JKth basis number covering interval J
        IKC=IKR+J-JK
        IF(IKC >= -1 .AND. IKC <= NFLEVG+2) THEN
          ZHH(J,JK)=FDD(ZDELTA(IJD) &
           & ,ZWAA(J,IKR),ZWBB(J,IKR),ZWCC(J,IKR),ZWDD(J,IKR)&
           & ,ZFAA(JK,IKC),ZFBB(JK,IKC),ZFCC(JK,IKC),ZFDD(JK,IKC))  
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  ZB(JR,JR  )=ZHH(1,1)+ZHH(2,2)+ZHH(3,3)+ZHH(4,4)
  IF(JR+1 <= IDIM) ZB(JR,JR+1)=ZHH(2,1)+ZHH(3,2)+ZHH(4,3)  
  IF(JR+2 <= IDIM) ZB(JR,JR+2)=ZHH(3,1)+ZHH(4,2)  
  IF(JR+3 <= IDIM) ZB(JR,JR+3)=ZHH(4,1)
  IF(JR-1 >= 1)    ZB(JR,JR-1)=ZHH(1,2)+ZHH(2,3)+ZHH(3,4)
  IF(JR-2 >= 1)    ZB(JR,JR-2)=ZHH(1,3)+ZHH(2,4)
  IF(JR-3 >= 1)    ZB(JR,JR-3)=ZHH(1,4)     

ENDDO

! COMPUTE PRODUCT SD*A^-1*B*SF^-1=ZDERI
! ------------------------------------

ZDERI(:,:)=0.0_JPRB
DO JC=1,IDIM
  DO JR=1,IDIM
    DO JK=1,IDIM
      ZDERI(JR,JC)=ZDERI(JR,JC)+ZB(JR,JK)*ZSFINV(JK,JC)
    ENDDO
  ENDDO
ENDDO

ZPROD(:,:)=0.0_JPRB
DO JC=1,IDIM
  DO JR=1,IDIM
    DO JK=1,IDIM
      ZPROD(JR,JC)=ZPROD(JR,JC)+ZAINV(JR,JK)*ZDERI(JK,JC)
    ENDDO
  ENDDO
ENDDO

ZDERI(:,:)=0.0_JPRB
DO JC=1,IDIM
  DO JR=1,IDIM
    DO JK=1,IDIM
      ZDERI(JR,JC)=ZDERI(JR,JC)+ZSD(JR,JK)*ZPROD(JK,JC)
    ENDDO
  ENDDO
ENDDO

! Boundary values chosen: linear extrapolation to nodes 0 and NFLEVG+1
!                         derivative at nodes 1 and NFLEVG
 
! linear extrapolation into nodal points at model boundaries

!--------------
! f_top  = 3/2 f_1 - 1/2 f_2
! f_surf = 3/2 f_L - 1/2 f_L-1
!--------------
ZW1T = +(ZDELTA(1)+ZDELTA(0))/ZDELTA(1)
ZW2T = - ZDELTA(0)/ZDELTA(1)

II=NFLEVG
ZW1B = +(ZDELTA(II-1)+ZDELTA(II))/ZDELTA(II-1)
ZW2B = - ZDELTA(II)/ZDELTA(II-1)

!--------------
! f_top  = f_1 
! f_surf = f_L 
!--------------
! ZW1 =  1.0_JPRB
! ZW2 =  0.0_JPRB

II=NFLEVG
DO JR=1,NFLEVG
  DO JC=3,NFLEVG-2
    YDVFE%RDERI(JR,JC)=ZDERI(JR+2,JC+2)
  ENDDO
  YDVFE%RDERI(JR,1)=ZDERI(JR+2,3)+ZW1T*ZDERI(JR+2,2)&
   & -ZDERI(JR+2,1)/(YDVETA%VFE_ETAF(2)-YDVETA%VFE_ETAF(1))
  YDVFE%RDERI(JR,2)=ZDERI(JR+2,4)+ZW2T*ZDERI(JR+2,2)&
   & +ZDERI(JR+2,1)/(YDVETA%VFE_ETAF(2)-YDVETA%VFE_ETAF(1))
  YDVFE%RDERI(JR,II-1)=ZDERI(JR+2,II+1)+ZW2B*ZDERI(JR+2,II+3)&
   & +ZDERI(JR+2,II+4)/(YDVETA%VFE_ETAF(II-1)-YDVETA%VFE_ETAF(II))
  YDVFE%RDERI(JR,II  )=ZDERI(JR+2,II+2)+ZW1B*ZDERI(JR+2,II+3)&
   & -ZDERI(JR+2,II+4)/(YDVETA%VFE_ETAF(II-1)-YDVETA%VFE_ETAF(II))
ENDDO

IF(MYPROC == 1.AND.(NPRINTLEV >= 1)) THEN

  WRITE(NULOUT,*) ""
  WRITE(NULOUT,*) " PRINTINGS IN SUNH_VERTFE3D "

  ! check the sum of RDERI rows
  WRITE(NULOUT,*) "(I) SUM OF RDERI ROWS (THE SAME AS f(x)=1)"
  DO JR=1,NFLEVG
    ZSUM=0.0_JPRB
    DO JC=1,NFLEVG
      ZSUM = ZSUM + YDVFE%RDERI(JR,JC)
    ENDDO
    WRITE(NULOUT,'(I2.2,1X,E20.14)') JR,ZSUM
  ENDDO

  ! print RDERI
  WRITE(NULOUT,*) ' MATRIX RDERI FOR NH MODEL:'
  DO JR=1,NFLEVG
    WRITE(NULOUT,'(8(1X,E9.3))') (YDVFE%RDERI(JR,JC),JC=1,NFLEVG)
  ENDDO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUNH_VERTFE3D',1,ZHOOK_HANDLE)
END SUBROUTINE SUNH_VERTFE3D
