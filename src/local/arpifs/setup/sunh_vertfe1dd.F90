SUBROUTINE SUNH_VERTFE1DD(YDGEOMETRY)

!**** *SUNH_VERTFE1DD*  - Setup VERTical Finite Element  Second Derivative based
!                         on linear (1) elements
!                         Version designed for the NH model.

!      weak form is used !!

!     Purpose.
!     --------
!           Compute matrix for the vertical second order differentiation based
!           on linear elements

!**   Interface.
!     ----------

!     *CALL* SUNH_VERTFE1DD

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

!     Modifications.
!     --------------
!        Original : 2006-06
!        K. Yessad (Sep 2008): use VFE_ETAH and VFE_ETAF.
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!        K. Yessad (July 2014): Move some variables.
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
REAL(KIND=JPRB)    :: ZDELTA(-1:YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZFAA(2,0:YDGEOMETRY%YRDIMV%NFLEVG+1), ZFBB(2,0:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB)    :: ZDAA(2,0:YDGEOMETRY%YRDIMV%NFLEVG+1), ZDBB(2,0:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB)    :: ZWAA(2,0:YDGEOMETRY%YRDIMV%NFLEVG+1), ZWBB(2,0:YDGEOMETRY%YRDIMV%NFLEVG+1)
REAL(KIND=JPRB)    :: ZA(YDGEOMETRY%YRDIMV%NFLEVG+2,YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZAINV(YDGEOMETRY%YRDIMV%NFLEVG+2,YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZWORK(2*(YDGEOMETRY%YRDIMV%NFLEVG+3))
REAL(KIND=JPRB)    :: ZPROD(YDGEOMETRY%YRDIMV%NFLEVG+2,YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZB(YDGEOMETRY%YRDIMV%NFLEVG+2,YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZDERI(YDGEOMETRY%YRDIMV%NFLEVG+2,YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZSF(YDGEOMETRY%YRDIMV%NFLEVG+2,YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZSFINV(YDGEOMETRY%YRDIMV%NFLEVG+2,YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZSD(YDGEOMETRY%YRDIMV%NFLEVG+2,YDGEOMETRY%YRDIMV%NFLEVG+2)
REAL(KIND=JPRB)    :: ZEPS, ZHH(2,2)
REAL(KIND=JPRB)    :: ZBMAX, ZD2, ZD3
REAL(KIND=JPRB)    :: ZDET
REAL(KIND=JPRB)    :: FHH, FDD, FPOLY
REAL(KIND=JPRB)    :: ZW2, ZW1, ZSUM
REAL(KIND=JPRB)    :: XA,XB,YA,YB,XX
REAL(KIND=JPRB)    :: ZHOOK_HANDLE

#include "minv.h"

!     ------------------------------------------------------------------

! indefiniteINTEGRAL(EJ*EI)dx
FHH(XX,XA,XB,YA,YB)=&
 & XX*( XA*YA +&
 & XX*((XA*YB+XB*YA)/2._JPRB +&
 & XX*( XB*YB/3._JPRB )))           

! indefiniteINTEGRAL(EJ'*EI')dx
FDD(XX,XA,XB,YA,YB)= XX*( XB*YB )

FPOLY(XX,XA,XB)=XA+XX*XB
!FDERI(XX,XA,XB)=XB

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUNH_VERTFE1DD',0,ZHOOK_HANDLE)
ASSOCIATE(YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG)!     ------------------------------------------------------------------

IDIM=NFLEVG+2
ZEPS=100.0_JPRB*TINY(1.0_JPRB)

! compute distance between nodes (full eta levels)
ZDELTA( 0)=YDVETA%VFE_ETAF(1)-YDVETA%VFE_ETAH(0)
ZDELTA(-1)=ZDELTA(0)
ZDELTA(-1)=ZDELTA(0)
DO JK=1,NFLEVG-1
  ZDELTA(JK)=YDVETA%VFE_ETAF(JK+1)-YDVETA%VFE_ETAF(JK)
ENDDO
ZDELTA(NFLEVG)=YDVETA%VFE_ETAH(NFLEVG) -YDVETA%VFE_ETAF(NFLEVG)
ZDELTA(NFLEVG+1)=ZDELTA(NFLEVG)
ZDELTA(NFLEVG+2)=ZDELTA(NFLEVG)

! compute basis functions    

ZBMAX=1.0_JPRB
  
DO JK=0,NFLEVG+1
  ZD2=ZDELTA(JK-1)
  ZD3=ZDELTA(JK  )

  ! rising linear function
  ZFAA(1,JK)=0.0_JPRB
  ZFBB(1,JK)=ZBMAX/ZD2
  
  ! decreasing linear function
  ZFAA(2,JK)= ZBMAX
  ZFBB(2,JK)=-ZBMAX/ZD3

ENDDO

! derivative basis
ZDAA(:,:)=ZFAA(:,:)
ZDBB(:,:)=ZFBB(:,:)

! weighting functions
ZWAA(:,:)=ZFAA(:,:)
ZWBB(:,:)=ZFBB(:,:)

! COMPUTE MATRIX SF ( SF(I,J)=BASISFUNC_J(X_I) )  (from FE space of the function
!                                                 to be differenciated into physical space)
! ---------------------------------------------

ZSF(:,:)=0.0_JPRB
DO JR=1,IDIM
  IKC = JR-1
  ZSF(JR,JR)=FPOLY(0.0_JPRB,ZFAA(2,IKC),ZFBB(2,IKC))
ENDDO

! invert SF:       (ZSFINV projects the function to be integrated
!                  from physical space onto FE space)
ZDET=0.0_JPRB
ZSFINV(:,:)=ZSF(:,:)*3._JPRB/2._JPRB
CALL MINV(ZSFINV,IDIM,IDIM,ZWORK,ZDET,ZEPS,0,1)
ZSFINV(:,:)=ZSFINV(:,:)*3._JPRB/2._JPRB

! COMPUTE MATRIX SD ( SD(I,J)=BASISDER_J(X_I) )  (To project the derivative from FE space
!                                                to physical space 
! -------------------------------------------

ZSD(:,:)=0.0_JPRB
DO JR=1,IDIM
  IKC = JR-1
  ZSD(JR,JR)=FPOLY(0.0_JPRB,ZDAA(2,IKC),ZDBB(2,IKC))
ENDDO

! COMPUTE MATRIX A
! ----------------

ZA(:,:)=0.0_JPRB

! JR is the row number of matrix A
DO JR=1,IDIM
  ! IKR is the index of the test function
  IKR=JR-1
  ! ZHH is an auxiliary 2x2 matrix
  ZHH(:,:)=0.0_JPRB
  ! J is the index of the interval belonging to function IKR
  DO J=1,2
    ! IJD is the absolute number of interval J
    IJD=IKR-2+J
    ! Integration is from node 0 to node NFLEVG+1
    IF(IJD >= 0 .AND. IJD <= NFLEVG) THEN
      DO JK=1,2
        ! IKC is the JKth basis number covering interval J
        IKC=IKR+J-JK
        IF(IKC >= 0.AND. IKC <= NFLEVG+1) THEN
          ZHH(J,JK)=FHH(ZDELTA(IJD),ZWAA(J,IKR),ZWBB(J,IKR) &
           & ,ZFAA(JK,IKC),ZFBB(JK,IKC))  
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  ZA(JR,JR)=ZHH(1,1)+ZHH(2,2) 
  IF(JR+1 <= IDIM) ZA(JR,JR+1)=ZHH(2,1)          
  IF(JR-1 >= 1)    ZA(JR,JR-1)=ZHH(1,2)          
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

! JR is the row number of matrix A
DO JR=1,IDIM
  ! IKR is the index of the test function
  IKR=JR-1
  ! ZHH is an auxiliary 2x2 matrix
  ZHH(:,:)=0.0_JPRB

  ! J is the index of the interval belonging to function IKR
  DO J=1,2
    ! IJD is the absolute number of interval J
    IJD=IKR-2+J
    ! Integration is from node 0 to node NFLEVG+1
    IF(IJD >= 0 .AND. IJD <= NFLEVG) THEN
      DO JK=1,2
        ! IKC is the JKth basis number covering interval J
        IKC=IKR+J-JK
        IF(IKC >= 0.AND. IKC <= NFLEVG+1) THEN
          ZHH(J,JK)=FDD(ZDELTA(IJD),ZWAA(J,IKR),ZWBB(J,IKR)&
           & ,ZFAA(JK,IKC),ZFBB(JK,IKC))  
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  ZB(JR,JR  )=-(ZHH(1,1)+ZHH(2,2))
  IF(JR+1 <= IDIM) ZB(JR,JR+1)=-ZHH(2,1)          
  IF(JR-1 >= 1)    ZB(JR,JR-1)=-ZHH(1,2)          
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
ZW1 =  1.5_JPRB
ZW2 = -0.5_JPRB

II  = NFLEVG
DO JR=1,NFLEVG
  DO JC=2,NFLEVG-1
    YDVFE%RDDERI(JR,JC)=ZDERI(JR+1,JC+1)
  ENDDO
  YDVFE%RDDERI(JR,1)=ZDERI(JR+1,2)+ZW1*ZDERI(JR+1,1)
  YDVFE%RDDERI(JR,2)=ZDERI(JR+1,3)+ZW2*ZDERI(JR+1,1)
  YDVFE%RDDERI(JR,II-1)=ZDERI(JR+1,II  )+ZW2*ZDERI(JR+1,II+2)
  YDVFE%RDDERI(JR,II  )=ZDERI(JR+1,II+1)+ZW1*ZDERI(JR+1,II+2)
ENDDO

! check the sum of RDDERI rows
WRITE(NULOUT,*) "(I) SUM OF RDDERI ROWS (THE SAME AS DERIVATIVE OF f(x)=1)"
DO JR=1,NFLEVG
  ZSUM=0.0_JPRB
  DO JC=1,NFLEVG
    ZSUM = ZSUM + YDVFE%RDDERI(JR,JC)
  ENDDO
  WRITE(NULOUT,'(I2.2,1X,F12.8)') JR,ZSUM
ENDDO

IF(MYPROC == 1.AND.(NPRINTLEV >= 1)) THEN

  WRITE(NULOUT,*) ""
  WRITE(NULOUT,*) " PRINTINGS IN SUNH_VERTFE1DD "

  ! check the sum of RDDERI rows
  WRITE(NULOUT,*) "(I) SUM OF RDDERI ROWS (THE SAME AS DERIVATIVE OF f(x)=1)"
  DO JR=1,NFLEVG
    ZSUM=0.0_JPRB
    DO JC=1,NFLEVG
      ZSUM = ZSUM + YDVFE%RDDERI(JR,JC)
    ENDDO
    WRITE(NULOUT,'(I2.2,1X,E20.14)') JR,ZSUM
  ENDDO

  ! print RDDERI
  WRITE(NULOUT,*) ' MATRIX RDDERI FOR NH MODEL:'
  DO JR=1,NFLEVG
    WRITE(NULOUT,'(8(1X,E9.3))') (YDVFE%RDDERI(JR,JC),JC=1,NFLEVG)
  ENDDO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUNH_VERTFE1DD',1,ZHOOK_HANDLE)
END SUBROUTINE SUNH_VERTFE1DD
