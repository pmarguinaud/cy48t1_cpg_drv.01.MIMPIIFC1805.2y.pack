SUBROUTINE EWRLSGRAD(YDGEOMETRY,YDML_GCONF,YDML_LBC,YDELBC_FIELDS)

!**** *EWRLSGRAD* Writes the large scale coupling data from gt3buf
!                into an FA-file

!     Purpose.
!     --------
!     Save the adjoint with respect to the large scale gridpoint coupling
!     into an FA-file (== gradient to LBC)
!     The routine is shaped in such a way that it should work both in
!     shared and distributed memory environment (but DM not yet tested...)

!**   Interface.
!     ----------
!        *CALL* *EWRLSGRAD*

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        Too many

!     Method.
!     -------
!     The "all in one" technique is favoured: everything is done inside
!     ewrlsgrad:
!     a.- copy the local adjoint coupling data from the buffer to a local
!         array (contains E+I zone only) via esc2r
!     b.- open the gradient FA-file
!     c.- if more than 1 proc are used, those which are not iomaster send
!         their data immediatly
!     d.- the iomaster waits for receiving all the messages (if any),
!         and copies its own data to the final gridpoint array zreelg
!     e.- the fields are written in the file, field by field and level
!         by level (horizontal gridpoint fields)

!     For sure, the technical choices are far from optimum as far as
!     distributed input/ouput are concerned, because only one and sole
!     processor is writing after recollecting every thing from the others.
!     But it should be kept in mind that this routine only is called by
!     configuration 801 and is not devoted to operations.

!     Externals.
!     ----------
!       MPL_SEND,MPL_RECV,MPL_PROBE,FAITOU,FANDAR,FAISAN,FAIRME
!       ESC2R

!     Reference.
!     ----------
!       None

!     Author.
!     -------
!      Claude Fischer    * CNRM/GMAP *
!      Original : 99-03-15

!     Modifications.
!     --------------
!      C. Fischer     01-03-15  new MPL message passing
!      P. Smolikova   02-09-30  interface to SC2CGAP for d4 in NH
!      A. Bogatchev   03-10-06  remove LSPQ,LSPL,LSPI, introducing GFL structure instead
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R.ZAABOUL     08-Dec-2003 Introduction of ZGT3GMV, ZGT3GMVS and ZGT3GFL
!      C. Fischer    04-02-25 Remove ng3d3 and ng2d3 but the routine remains
!                             gmv and gfl uncompliant (sorry ...)
!      K. Yessad (Jan 2011): new architecture for LBC modules and set-up.
!      K. Yessad (May 2012): externalisation of coupling.
!      B. Bochenek 11-04-2013 - Phasing cy40, coherence with modified modules
!      R. El Khatib 28-Jan-2015 Optimizations
!      B. Bochenek (Apr 2015): Phasing: move some variables.
!      O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMLUN   ,ONLY : NULOUT   ,NULUSR3
USE YOMCT0   ,ONLY : CNMEXP
USE YOMDFI   ,ONLY : NSTDFI   ,NSTDFIA
USE YOMCT3   ,ONLY : NSTEP
USE YOMDYNA  ,ONLY : YRDYNA
USE YOMINI   ,ONLY : LBIAS
USE YOMRIP0  ,ONLY : NINDAT   ,NSSSSS
USE YOMMP0   ,ONLY : LSPLIT   ,MYPROC    ,NPROC    ,NPRGPNS  ,NPRGPEW, LOUTPUT ,NPRINTLEV
USE YOMTAG   ,ONLY : MTAGDISTGP
USE MPL_MODULE
USE YEMLBC_INIT,ONLY : LTENC    ,LALLTC   ,LQCPL    ,LCCPL    ,LRDLSG   ,NECRIPL
USE YEMLBC_MODEL,ONLY : TELBC_MODEL
USE YEMLBC_FIELDS,ONLY : TELBC_FIELDS
USE YOMOPH0  , ONLY : CNMCA    ,LINC

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT)   :: YDGEOMETRY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
TYPE(TELBC_MODEL),INTENT(INOUT):: YDML_LBC
TYPE(TELBC_FIELDS),INTENT(INOUT):: YDELBC_FIELDS

REAL(KIND=JPRB) :: ZGT3GMV(YDGEOMETRY%YRDIM%NDLON,YDGEOMETRY%YRDIMV%NFLEVG,&
 & YDELBC_FIELDS%YYTGMVCPL%NDIM)
REAL(KIND=JPRB) :: ZGT3GMVS(YDGEOMETRY%YRDIM%NDLON,YDELBC_FIELDS%YYTGMVSCPL%NDIM)
REAL(KIND=JPRB) :: ZGT3GFL(YDGEOMETRY%YRDIM%NDLON,YDGEOMETRY%YRDIMV%NFLEVG,YDELBC_FIELDS%NDIMCPL)
REAL(KIND=JPRB),ALLOCATABLE:: ZBUFT3(:)
REAL(KIND=JPRB),ALLOCATABLE:: ZBUFR(:)
REAL(KIND=JPRB),ALLOCATABLE:: ZREELG(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: ZREELV(:)

INTEGER(KIND=JPIM) :: IDATEF(11)

CHARACTER :: CLFNAM*16,CLF1*3,CLF2*4,CLF3*5,CLINC*4
CHARACTER :: CLLEV*4
CHARACTER (LEN = 12) :: CLFIELD(YDELBC_FIELDS%YYTGMVCPL%NDIM+YDELBC_FIELDS%YYTGMVSCPL%NDIM+&
 & YDELBC_FIELDS%NDIMCPL)

INTEGER(KIND=JPIM) :: IDGLL, IDIM2D, IDIM3D, IDIM_BUFR, IDIM_BUFT3,&
 & IFIELDOFF, IFIELDS, IFLD,&
 & IFSTLAT, IFSTLATOFF, IGGP, IGP, IGPTOTH,&
 & ILENMES, ILENR,&
 & INBARI, INBARP, INIMES, INUML,&
 & IOFFSET, IOMASTER,&
 & IREP, ISENDER,&
 & ITAG, ITIMLEV, ITRUE_FIELDS,&
 & JFLD, JGL, JLEV, JLON, JPROC, I_KINC, IINC  , ISWAP

LOGICAL :: LLSPEC
LOGICAL :: LLTENC,LLALLTC
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "esc2r.h"

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EWRLSGRAD',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YGFL=>YDML_GCONF%YGFL,YDRIP=>YDML_GCONF%YRRIP)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & NDGLG=>YDDIM%NDGLG, NDGLL=>YDDIM%NDGLL, NDGSAL=>YDDIM%NDGSAL, &
 & NDGENL=>YDDIM%NDGENL, NDLON=>YDDIM%NDLON, &
 & NFRSTLAT=>YDMP%NFRSTLAT, &
 & NGFL_EXT=>YGFL%NGFL_EXT, YQ=>YGFL%YQ, YL=>YGFL%YL, YI=>YGFL%YI, &
 & TSTEP=>YDRIP%TSTEP, NSTOP=>YDRIP%NSTOP )
!     ------------------------------------------------------------------

!*   1. START UP
!    1.1 Initialize and test message passing options

WRITE(NULOUT,*) ' *** EWRLSGRAD *** '

CALL ABOR1('EWRLSGRAD: GFL NOT YET FULLY CODED')

IOMASTER=1
INUML=0
IF (NPRGPEW > 1) THEN
  WRITE(NULOUT,*) ' EWRLSGRAD: ERROR THIS ROUTINE DOES NOT ',&
   & 'ALLOW NPRGPEW STRICTLY BIGGER THAN 1'  
  CALL ABOR1('EWRLSGRAD: ABORT CALLED')
ENDIF
IF (LSPLIT) THEN
  WRITE(NULOUT,*) ' EWRLSGRAD: ERROR THIS ROUTINE DOES NOT ',&
   & 'ALLOW LSPLIT TO BE TRUE (FIXED NDLON IS REQUIRED)'  
  CALL ABOR1('EWRLSGRAD: ABORT CALLED')
ENDIF

!    1.2 Initialize and test grid point options

IF (YRDYNA%LRUBC.OR.LQCPL) THEN
  WRITE(NULOUT,*) ' EWRLSGRAD: LRUBC AND LQCPL ARE FORBIDDEN'
  CALL ABOR1('EWRLSGRAD: ABORT CALLED')
ENDIF
IF (NECRIPL /= 1) THEN
  WRITE(NULOUT,*) ' EWRLSGRAD: NECRIPL MUST BE EQUAL TO 1'
  CALL ABOR1('EWRLSGRAD: ABORT CALLED')
ENDIF
IF (YL%LSP.OR.YI%LSP) THEN
  WRITE(NULOUT,*) ' EWRLSGRAD: NO LIQUID OR ICE WATER WILL BE ',&
   & 'WRITTEN IN COUPLING GRADIENT FILE'  
ENDIF
IF (NGFL_EXT > 0) THEN
  WRITE(NULOUT,*) ' EWRLSGRAD: NO EXTRA-GFL WILL BE ',&
   & 'WRITTEN IN COUPLING GRADIENT FILE'  
ENDIF

IFIELDS=(YDELBC_FIELDS%YYTGMVCPL%NDIM+YDELBC_FIELDS%NDIMCPL)*NFLEVG+&
  & YDELBC_FIELDS%YYTGMVSCPL%NDIM
IOFFSET=0
IDIM3D=0
IDIM2D=0

ZGT3GMV=0.0_JPRB
ZGT3GMVS=0.0_JPRB
ZGT3GFL=0.0_JPRB

CLFIELD(:)='LSCGRAD_UNKN'

! CLFIELD for GMV:
DO JFLD=1,YDELBC_FIELDS%YYTGMVCPL%NDIM
  IOFFSET=IOFFSET+1
  CLFIELD(IOFFSET)=YDELBC_FIELDS%CCFIELD_GMV(JFLD)
ENDDO

! CLFIELD for GFL (only specific humidity for the time being):
IF ( YQ%LSP ) THEN
  IOFFSET=IOFFSET+1
  CLFIELD(IOFFSET)=YDELBC_FIELDS%CCFIELD_GFL(1)
ENDIF
IDIM3D=IOFFSET

! CLFIELD for GMVS:
DO JFLD=1,YDELBC_FIELDS%YYTGMVSCPL%NDIM
  IOFFSET=IOFFSET+1
  CLFIELD(IOFFSET)=YDELBC_FIELDS%CCFIELD_GMVS(JFLD)
ENDDO
IDIM2D=IOFFSET-IDIM3D

! ITRUE_FIELDS = number of coupled fields.
ITRUE_FIELDS=IDIM3D*NFLEVG+IDIM2D
IOFFSET=0

IF (NPRINTLEV >= 1) THEN
  WRITE(NULOUT,*) ' Total number of fields in GMVCPL+GMVSCPL+GFLCPL : ',&
   & IFIELDS,' Number of coupled fields : ',ITRUE_FIELDS  
  WRITE(NULOUT,*) ' Split into 3D fields (total)     : ',&
   & YDELBC_FIELDS%YYTGMVCPL%NDIM+YDELBC_FIELDS%NDIMCPL  ,' (really coupled)         : ',IDIM3D  
  WRITE(NULOUT,*) ' Split into 2D fields (total)     : ',&
   & YDELBC_FIELDS%YYTGMVSCPL%NDIM  ,' (really coupled)         : ',IDIM2D  
ENDIF

!     ------------------------------------------------------------------

!*   2. COPY THE PROCESSOR LOCAL LBC-COUPLING BUFFER

!    2.2 Read in gt3buf all the local rows of latitude
!         by the way, the buffer-type of organisation of arpege/aladin
!         is replaced by a horizontal staggering

IDIM_BUFT3=ITRUE_FIELDS*NDGLL*NDLON
ALLOCATE(ZBUFT3(IDIM_BUFT3))
ZBUFT3(:)=0.0_JPRB
ITIMLEV=1
DO JGL=NDGSAL,NDGENL
  IOFFSET=(JGL-NDGSAL)*NDLON

  IF ((.NOT.LQCPL).AND.LRDLSG) THEN
    CALL ABOR1('EWRLSGRAD: case where LQCPL=F with LRDLSG=T not implemented')
  ENDIF

  LLTENC=LTENC   ! ky: LLTENC=LTENC or LLTENC=.FALSE. in this call to ESC2R?
  LLALLTC=LALLTC ! ky: LLALLTC=LALLTC or LLALLTC=.FALSE. in this call to ESC2R?

  ! * externalisable part of coupling (temporal interpolation).
  !   This call to ESC2R does not treat horizontal derivatives (useless here)
  !   contrary to the one called by ECOUPL1.
  ISWAP=SIZE(YDELBC_FIELDS%GMVSCPL,DIM=3)
  CALL ESC2R(NDLON,NFLEVG,ISWAP,YDELBC_FIELDS%YYTGMVCPL%NDIM,YDELBC_FIELDS%YYTGMVSCPL%NDIM,&
   & YDELBC_FIELDS%NDIMCPL,ITIMLEV,&
   & NSTOP,NSTDFI,NSTDFIA,NSTEP,YDML_LBC%NETLS1,YDML_LBC%NEFRCL,LBIAS,LQCPL,LCCPL,LLTENC,LLALLTC,&
   & TSTEP,YDML_LBC%EWB,YDML_LBC%EWBDFIFW,YDML_LBC%EWBDFIBW,LDDFISTEP=.FALSE.,&
   & PGMVCPL=YDELBC_FIELDS%GMVCPL,PGMVSCPL=YDELBC_FIELDS%GMVSCPL,PGFLCPL=YDELBC_FIELDS%GFLCPL,&
   & PGT3GMV=ZGT3GMV,PGT3GMVS=ZGT3GMVS,PGT3GFL=ZGT3GFL)

  ! * GMV:
  DO JFLD=1,YDELBC_FIELDS%YYTGMVCPL%NDIM
    DO JLEV=1,NFLEVG
      DO JLON=1,NDLON
        ZBUFT3(JLON+IOFFSET)=ZGT3GMV(JLON,JLEV,JFLD)
      ENDDO
      IOFFSET=IOFFSET+NDGLL*NDLON
    ENDDO
  ENDDO

  ! * GFL (only specific humidity for the time being):
  IF (YQ%LSP) THEN
    DO JLEV=1,NFLEVG
      DO JLON=1,NDLON
        ZBUFT3(JLON+IOFFSET)=ZGT3GFL(JLON,JLEV,YQ%MPSP)
      ENDDO
      IOFFSET=IOFFSET+NDGLL*NDLON
    ENDDO
  ENDIF

  ! * GMVS:
  DO JLON=1,NDLON
    ZBUFT3(JLON+IOFFSET)=ZGT3GMVS(JLON,YDELBC_FIELDS%YYTGMVSCPL%MSP)
  ENDDO
ENDDO
IOFFSET=0

!     ------------------------------------------------------------------

!*   3. COMMUNICATE THE COUPLING DATA BETWEEN PROCS AND WRITE TO FILE
!    3.1 Open gradient file

IF (MYPROC == IOMASTER) THEN
! Open the LS gradient file
  INIMES=1
  INBARP=IFIELDS
  INBARI=0
  CLF1='LSG'
  CLF2=CNMEXP(1:4)
  CLF3='GT3B+'
  IF(LINC) THEN
    I_KINC=NINT(REAL(NSTEP,JPRB)*TSTEP/3600._JPRB)
  ELSE
    I_KINC=NSTEP
    IF(NSTEP==1) I_KINC=0
  ENDIF
  WRITE(CLINC,'(I4.4)') I_KINC
  WRITE(CLFNAM,'(A3,A4,A5,A4)') CLF1,CLF2,CLF3,CLINC
  CALL FAITOU(IREP,NULUSR3,.TRUE.,CLFNAM,'UNKNOWN',&
   & .TRUE.,.TRUE.,INIMES,INBARP,INBARI,CNMCA)  
! Setup date
  IDATEF(1) = NINDAT/10000
  IDATEF(2) = (NINDAT-10000*IDATEF(1))/100
  IDATEF(3) = NINDAT-10000*IDATEF(1)-100*IDATEF(2)
  IDATEF(4) = NSSSSS/3600
  IDATEF(5) = 0
  IDATEF(6) = 1
  IINC=NINT(REAL(NSTEP,JPRB)*TSTEP/3600._JPRB)
  IDATEF(7) = IINC
  IDATEF(8) = 0
  IDATEF(9) = 10
  IF(NSTEP==1) IDATEF(9) = 1
  IDATEF(10)= 0
  IDATEF(11)= 0
  CALL FANDAR(IREP,NULUSR3,IDATEF)
ENDIF

!    3.2 Send data if more than 1 proc

INUML=INUML+1
ITAG=MTAGDISTGP+INUML*IOMASTER*NPROC+MYPROC
IF ((NPROC > 1).AND.(MYPROC /= IOMASTER)) THEN

  CALL MPL_SEND(ZBUFT3(1:IDIM_BUFT3),KDEST=IOMASTER,KTAG=ITAG,&
   & CDSTRING='EWRLSGRAD:')  

ENDIF

!    3.3 Receive data and copy to the gridpoint array of iomaster

IF (MYPROC == IOMASTER) THEN
! begin of test and work only if iomaster

  ALLOCATE(ZREELG(NDGLG*NDLON,ITRUE_FIELDS))
  ALLOCATE(ZREELV(NDGLG*NDLON))
  ZREELG(:,:)=0.0_JPRB
  ZREELV(:)=0.0_JPRB

  IF (NPROC > 1) THEN

    IDIM_BUFR=NDGLG*NDLON*ITRUE_FIELDS
    ALLOCATE(ZBUFR(IDIM_BUFR))
    ZBUFR(:)=0.0_JPRB

    DO JPROC=1,NPRGPNS
      WRITE(NULOUT,*)'JPROC= ',JPROC
      IF (JPROC /= IOMASTER) THEN
        ITAG=MTAGDISTGP+INUML*IOMASTER*NPROC+JPROC
        ILENMES=IDIM_BUFR
        CALL MPL_RECV(ZBUFR(1:ILENMES),KSOURCE=JPROC,KTAG=ITAG,&
         & KOUNT=ILENR,KFROM=ISENDER,CDSTRING='EWRLSGRAD:')  
        IDGLL=ILENR/(NDLON*ITRUE_FIELDS)

        IFSTLAT=NFRSTLAT(JPROC)
        IFSTLATOFF=(IFSTLAT-1)*NDLON

        DO JFLD=1,ITRUE_FIELDS
          DO JGL=1,IDGLL
            DO JLON=1,NDLON
              IGP=JLON+(JGL-1)*NDLON+IFSTLATOFF
              IGGP=JLON+(JGL-1)*NDLON+(JFLD-1)*NDLON*IDGLL
              ZREELG(IGP,JFLD)=ZBUFR(IGGP)
            ENDDO
          ENDDO
        ENDDO

      ENDIF
    ENDDO

    DEALLOCATE(ZBUFR)

  ENDIF

  DO JFLD=1,ITRUE_FIELDS
    DO JGL=1,NDGLL
      DO JLON=1,NDLON
        IGP=JLON+(JGL-1)*NDLON
        IGGP=JLON+(JGL-1)*NDLON+(JFLD-1)*NDLON*NDGLL
        ZREELG(IGP,JFLD)=ZBUFT3(IGGP)
      ENDDO
    ENDDO
  ENDDO

!    3.4 Write the LS data into FA-file

  IGPTOTH=NDGLG*NDLON
  IFIELDOFF=0
  LLSPEC=.FALSE.

  DO JFLD=1,IDIM3D
    DO JLEV=1,NFLEVG
      WRITE(CLLEV,'(A1,I3.3)') 'S',JLEV
      IFLD=JLEV+(JFLD-1)*NFLEVG
      ZREELV(1:IGPTOTH)=ZREELG(1:IGPTOTH,IFLD)
      CALL FAIENC(IREP,NULUSR3,CLLEV,JLEV,CLFIELD(JFLD),ZREELV,LLSPEC)
      IF (LOUTPUT)&
       & WRITE(NULOUT,'(2X,A,I3.3,'' WRITTEN TO ALADIN FILE '')')&
       & CLFIELD(JFLD),JLEV  
    ENDDO
  ENDDO
  IFIELDOFF=IFIELDOFF+NFLEVG*IDIM3D
  DO JFLD=1,IDIM2D
    CLLEV='SURF'
    ZREELV(1:IGPTOTH)=ZREELG(1:IGPTOTH,JFLD+IFIELDOFF)
    CALL FAIENC(IREP,NULUSR3,CLLEV,1,CLFIELD(JFLD+IDIM3D),ZREELV,LLSPEC)
    IF (LOUTPUT)&
     & WRITE(NULOUT,'(2X,A,'' WRITTEN TO ALADIN FILE '')')&
     & CLFIELD(JFLD+IDIM3D)  
  ENDDO

!    3.5 Close the gradient file

  CALL FAIRME(IREP,NULUSR3,'UNKNOWN')

  DEALLOCATE(ZREELG)
  DEALLOCATE(ZREELV)

! end of test on iomaster
ENDIF

DEALLOCATE(ZBUFT3)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('EWRLSGRAD',1,ZHOOK_HANDLE)
END SUBROUTINE EWRLSGRAD
