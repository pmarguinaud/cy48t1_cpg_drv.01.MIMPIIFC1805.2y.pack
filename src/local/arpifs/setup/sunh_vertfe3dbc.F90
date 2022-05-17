SUBROUTINE SUNH_VERTFE3DBC(YDGEOMETRY)

!**** *SUNH_VERTFE3DBC*  - Setup VERTical Finite Element Derivative based on
!                          cubic (3) elements: first order derivative.
!                          Version designed for the NH model.
!                          Takes account of the upper and lower boundary cond.

!     Purpose.
!     --------
!           Compute matrix for the vertical differentiation based on
!           cubic B-spline finite elements

!**   Interface.
!     ----------

!     *CALL* SUNH_VERTFE3DBC

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
!        Original : 2008-01

!     Modifications.
!     --------------
!   K. Yessad (Sep 2008): use VFE_ETAH and VFE_ETAF
!   T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!   K. Yessad (July 2014): Move some variables.
!   P. Smolikova (Sep 2020): VFE pruning.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMMP0   , ONLY : NPRINTLEV, MYPROC
USE YOMLUN   , ONLY : NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
INTEGER(KIND=JPIM) :: IDIM, J, JK, IKC, IJD, JR, IKR, J1, J2, J3, JC
INTEGER(KIND=JPIM) :: II, IROW, INM(2)
REAL(KIND=JPRB)    :: ZDELTA(-3:YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZFAA(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZFBB(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZFCC(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZFDD(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZDAA(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZDBB(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZDCC(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZDDD(4,-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZA(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZAINV(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZWORK(2*(YDGEOMETRY%YRDIMV%NFLEVG+3))
REAL(KIND=JPRB)    :: ZPROD(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZB(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZDERI(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZSF(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZSFINV(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZSD(YDGEOMETRY%YRDIMV%NFLEVG+4,YDGEOMETRY%YRDIMV%NFLEVG+4)
REAL(KIND=JPRB)    :: ZIN(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZOU1(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    :: ZOU2(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZOUD(YDGEOMETRY%YRDIMV%NFLEVG,YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZEPS, ZALPHA(4,4), ZBETA(2,2), ZHH(4,4)
REAL(KIND=JPRB)    :: ZRHS(2)
REAL(KIND=JPRB)    :: ZD0, ZD0H
REAL(KIND=JPRB)    :: ZBMAX, ZD1, ZD2, ZD3, ZD4, ZDETB, ZDEL, ZDEL2, ZDEL3
REAL(KIND=JPRB)    :: ZDDD1, ZDDD2, ZDDD3, ZDDD4, ZDET
REAL(KIND=JPRB)    :: ZAAA2, ZBBB2, ZCCC2, ZAAA3, ZBBB3, ZCCC3
REAL(KIND=JPRB)    :: FHH, FDD, FPOLY, FDERI
REAL(KIND=JPRB)    :: XA,XB,XC,XD,YA,YB,YC,YD,XX
REAL(KIND=JPRB)    :: ZARG, ZSUM

LOGICAL            :: LLGP_AT01               ! grid points at eta=0 and eta=1

REAL(KIND=JPRB)    :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "verdisint.intfb.h"

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

IF (LHOOK) CALL DR_HOOK('SUNH_VERTFE3DBC',0,ZHOOK_HANDLE)
ASSOCIATE(YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)!     ------------------------------------------------------------------

IDIM=NFLEVG+4
ZEPS=100.0_JPRB*TINY(1.0_JPRB)

LLGP_AT01 = .TRUE.

! compute distance between nodes (full eta levels)

IF( LLGP_AT01 )THEN
  ZDELTA(0)=(YDVETA%VFE_ETAF(1)-YDVETA%VFE_ETAH(0))
ELSE
  ZDELTA(0)=2._JPRB*(YDVETA%VFE_ETAF(1)-YDVETA%VFE_ETAH(0))
ENDIF

ZD0=ZDELTA(0)
ZD0H=0.5_JPRB*ZD0
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
ZDELTA(NFLEVG+4)=ZDELTA(NFLEVG)

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

! COMPUTE MATRIX SF ( SF(I,J)=BASISFUNC_J(X_I) )
!  from FE space of the function to be differenciated into physical space
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

!-------------------------
! the boundary conditions
!-------------------------

! CONDITION 1 - derivative at eta=0 prescribed

! in the case of LLGP_AT01=T it is a nodal point 0
! in the case of LLGP_AT01=F it is a middle of interval (0,1)
! derivative at model top diminish to 0
! IROW - the index in which derivative is in input vector
! II   - the nodal point index on which derivative is valid
!      - if II=0 and ARG=ZDELTA(0)/2 then the values are evaluated at model top
!        where eta=0 in the case of LLGP_AT01=F

IROW = 1
II   = 0
IF( LLGP_AT01 )THEN
  ZARG = 0.0_JPRB
ELSE
  ZARG = 0.5_JPRB*ZDELTA(II)
ENDIF
ZSF(IROW,:)=0.0_JPRB
ZSF(IROW,1)=FDERI(ZARG,ZFAA(4, II-1),ZFBB(4, II-1),ZFCC(4, II-1),ZFDD(4, II-1))
ZSF(IROW,2)=FDERI(ZARG,ZFAA(3, II  ),ZFBB(3, II  ),ZFCC(3, II  ),ZFDD(3, II  ))
ZSF(IROW,3)=FDERI(ZARG,ZFAA(2, II+1),ZFBB(2, II+1),ZFCC(2, II+1),ZFDD(2, II+1))
ZSF(IROW,4)=FDERI(ZARG,ZFAA(1, II+2),ZFBB(1, II+2),ZFCC(1, II+2),ZFDD(1, II+2))

! CONDITION 2 - functional value at eta=0 is prescribed

IROW = 2
II   = 0
IF( LLGP_AT01 )THEN
  ZARG = 0.0_JPRB
ELSE
  ZARG = 0.5_JPRB*ZDELTA(II)
ENDIF
ZSF(IROW,:)=0.0_JPRB
ZSF(IROW,1)=FPOLY(ZARG,ZFAA(4, II-1),ZFBB(4, II-1),ZFCC(4, II-1),ZFDD(4, II-1))
ZSF(IROW,2)=FPOLY(ZARG,ZFAA(3, II  ),ZFBB(3, II  ),ZFCC(3, II  ),ZFDD(3, II  ))
ZSF(IROW,3)=FPOLY(ZARG,ZFAA(2, II+1),ZFBB(2, II+1),ZFCC(2, II+1),ZFDD(2, II+1))
ZSF(IROW,4)=FPOLY(ZARG,ZFAA(1, II+2),ZFBB(1, II+2),ZFCC(1, II+2),ZFDD(1, II+2))

! CONDITION 3 - functional value at eta=1 is prescribed

IROW  = IDIM-1
II    = NFLEVG
IF( LLGP_AT01 )THEN
  ZARG = ZDELTA(II)
ELSE
  ZARG = 0.5_JPRB*ZDELTA(II)
ENDIF
ZSF(IROW,:)=0.0_JPRB
ZSF(IROW,IDIM-3)=FPOLY(ZARG,ZFAA(4,II-1),ZFBB(4,II-1),ZFCC(4,II-1),ZFDD(4,II-1))
ZSF(IROW,IDIM-2)=FPOLY(ZARG,ZFAA(3,II  ),ZFBB(3,II  ),ZFCC(3,II  ),ZFDD(3,II  ))
ZSF(IROW,IDIM-1)=FPOLY(ZARG,ZFAA(2,II+1),ZFBB(2,II+1),ZFCC(2,II+1),ZFDD(2,II+1))
ZSF(IROW,IDIM  )=FPOLY(ZARG,ZFAA(1,II+2),ZFBB(1,II+2),ZFCC(1,II+2),ZFDD(1,II+2))

! CONDITION 4 - derivative value at eta=1 is prescribed

IROW  = IDIM
II    = NFLEVG
IF( LLGP_AT01 )THEN
  ZARG = ZDELTA(II)
ELSE
  ZARG = 0.5_JPRB*ZDELTA(II)
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
           & ,ZFAA(J,IKR),ZFBB(J,IKR),ZFCC(J,IKR),ZFDD(J,IKR) &
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
           & ,ZFAA(J,IKR),ZFBB(J,IKR),ZFCC(J,IKR),ZFDD(J,IKR)&
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

! check dimensions of RDERB 
INM = SHAPE( YDVFE%RDERB )
IF( INM(1) /= NFLEVG )THEN
  CALL ABOR1('(E) in SUNH_VERTFE3DBC: FIRST DIMENSION OF RDERB MUST BE NFLEVG ')
ENDIF
IF( INM(2) /= NFLEVG+2 )THEN
  CALL ABOR1('(E) in SUNH_VERTFE3DBC: 2ND DIMENSION OF RDERB MUST BE NFLEVG+2 ')
ENDIF

!------------------------------------------------ 
! Packing :

! ZDERI ( L+4 , L+4 )   -->     RDERB( L , L+2 )
!------------------------------------------------
! RDERB has L rows and L+2 columns

! RDERB acts on L+2 input vector f( f0, f1, f2, f3, ... , fL, fL+1)
! and provides L derivatives on full levels ( d1, d2, d3, ... , dL )

! WARNING:
! in this packing we assume that derivatives d(eta=0) and d(eta=1) are 0.
  
DO JR=1,NFLEVG
  DO JC=1,NFLEVG+2
    YDVFE%RDERB(JR,JC)=ZDERI(JR+2,JC+1)
  ENDDO
ENDDO

IF(MYPROC == 1.AND.(NPRINTLEV >= 1)) THEN
  
  WRITE(NULOUT,*) ""
  WRITE(NULOUT,*) " PRINTINGS IN SUNH_VERTFE3DBC "
  
  ! check the sum of RDERB rows (check if RDERB sum of row values is 0,
  !  equivalent to apply operator on f(eta) = constant).
  WRITE(NULOUT,*) ""
  WRITE(NULOUT,*) "(I) SUM OF RDERB ROWS (THE SAME AS f(x)=1)"
  DO JR=1,NFLEVG
    ZSUM=0.0_JPRB
    DO JC=1,NFLEVG+2
      ZSUM = ZSUM + YDVFE%RDERB(JR,JC)
    ENDDO
    WRITE(NULOUT,'(I2.2,1X,E20.14)') JR,ZSUM
  ENDDO

  ! check the sum of RDERB columns.
  WRITE(NULOUT,*) ""
  WRITE(NULOUT,*) "(I) SUM OF RDERB COLUMNS WEIGHTED BY [Delta eta]"
  DO JC=1,NFLEVG+2
    ZSUM=0.0_JPRB
    DO JR=1,NFLEVG
      ZSUM = ZSUM + YDVFE%RDERB(JR,JC)*(YDVETA%VFE_ETAH(JR)-YDVETA%VFE_ETAH(JR-1))
    ENDDO
    WRITE(NULOUT,'(I2.2,1X,E20.14)') JC,ZSUM
  ENDDO

  ! print RDERB
  WRITE(NULOUT,*) ""
  WRITE(NULOUT,*) ' MATRIX RDERB FOR NH MODEL:'
  DO JR=1,NFLEVG
    WRITE(NULOUT,'(8(1X,E9.3))') (YDVFE%RDERB(JR,JC),JC=1,NFLEVG+2)
  ENDDO

  ! check that, when we apply RDERB then RINTE, we obtain something not
  ! too far from identity, and that the management of boundaries is OK
  ! (for example, not too big drift at the bottom when prescribing the top value).
  ZIN(1:NFLEVG,1:NFLEVG+2)=0.0_JPRB
  DO JR=1,NFLEVG
    ZIN(JR,JR+1)=1.0_JPRB
  ENDDO
  CALL VERDISINT(YDVFE,'FDER','XX',NFLEVG,1,NFLEVG,NFLEVG,ZIN,ZOU1)
  CALL VERDISINT(YDVFE,'ITOP','XX',NFLEVG,1,NFLEVG,NFLEVG,ZOU1,ZOU2(1:NFLEVG,2:NFLEVG+2))
  ZOU2(1:NFLEVG,1)=0.0_JPRB
  ZOUD(1:NFLEVG,1:NFLEVG+2)=ZOU2(1:NFLEVG,1:NFLEVG+2)-ZIN(1:NFLEVG,1:NFLEVG+2)
  ! print ZOUD
  WRITE(NULOUT,*) ""
  WRITE(NULOUT,*) ' MATRIX RINTE(from top)*RDERB-IDENTITY FOR NH MODEL:'
  DO JR=1,NFLEVG
    WRITE(NULOUT,'(8(1X,E9.3))') (ZOUD(JR,JC),JC=1,NFLEVG+2)
  ENDDO

  WRITE(NULOUT,*) ""

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUNH_VERTFE3DBC',1,ZHOOK_HANDLE)
END SUBROUTINE SUNH_VERTFE3DBC
