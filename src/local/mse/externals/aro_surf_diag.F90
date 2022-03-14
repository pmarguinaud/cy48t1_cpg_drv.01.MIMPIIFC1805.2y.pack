SUBROUTINE ARO_SURF_DIAG (LDMPA, LDPGDFWR, LDHISFWR, KSTOP, KSTEP, &
                        & KFREQFS, KYPROC, KLON, PINRAIN, PINSNOW, &
                        & CDFILE)
USE PARKIND1, ONLY : JPRB, JPIM
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###################################
!
!!****  *ARO_SURF_DIAG* - Interface pour preparer l'appel des output de la surface externalisee
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    MODIFICATIONS
!!    -------------
!!      July 2007 : A.Alias ajout du traitement de la grille de GAUSS
!!      June 2008 : Y.Seity codage sur 4 chiffres de l'echeance du AROMOUT
!!      Dec  2008 : A.Alias modifications to write surface files in FA format
!!      Juin 2009 : B. Decharme Bug calcul NKPROMA at NBLOCK=0
!!      Aug  2009 : A.Alias modifications in writing surface files in FA format at NSTOP
!!                          adding arguments in aro_surf_diag
!!      Oct  2009 : Y.Seity Remove writting of INPR* and ACPR* fields in surfex file
!!      Jul  2010 : P. Marguinaud & R. El Khatib : I/O savings + V7 phasing
!!      Aug  2010 : P.Marguinaud : More arguments to goto_surfex
!!      Dec  2010 : A.ALias : Modifications to deal with SURFEX file in FA for ALADIN
!!                  SURFEX restart file added
!!       Jul 2011 : P.Marguinaud : use YSURFEX_CACHE_OUT
!!       Sep 2012 : P.Marguinaud : more initialization arguments + make them optional
!!       Dec 2013 : F.Taillefer  : enable selection of output fields to write
!!       Fev 2021 : A.Napoly  : add init of LFIRST_WRITE and NCPT_WRITE for V9 compatibility
!!       Fev 2021 : R. El Khatib : bugfix to allow mpi tasks fully in extension zone
!-------------------------------------------------------------------------------

USE MODD_SURFEX_ARO, ONLY : YSURFEX_ARO_CUR, YSURFEX_ARO_ALL
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_IO_SURF_ARO, ONLY :                                         &
   & NGPTOTG,                                                        &
   &                         NPROMA, NBLOCK,         NGPTOT,         &
   & NINDX1, NINDX2, NKPROMA, NBLOCKTOT,                             &
   &                                  NGPTOT_CAP,                    &
   &             SURFEX_FIELD_BUF_DEALLOC,  SURFEX_FIELD_BUF_DUMP,   &
   & SURFEX_FIELD_BUF_SET_RECORD, SURFEX_FIELD_BUF_PREALLOC,         &
   & SURFEX_FIELD_BUF_DEALLOC, SURFEX_FIELD_BUF_ADD,                 &
   & YSURFEX_CACHE_OUT, MYPROC

USE MODD_IO_SURF_FA, ONLY : CFILEOUT_FA, NUNIT_FA, &
                          & CDNOMC, IVERBFA,       &
                          & CFILEIN_FA, LFANOCOMPACT 
USE MODI_FLAG_UPDATE
USE MODI_FMWRIT
USE MODI_DIAG_SURF_ATM_N
USE MODI_WRITE_DIAG_SURF_ATM_N
USE MODI_WRITE_PGD_SURF_ATM_N
USE MODI_WRITE_SURF_ATM_N
USE MODD_WRITE_SURF_ATM, ONLY : LFIRST_WRITE,NCPT_WRITE
!
IMPLICIT NONE
!
LOGICAL,             INTENT(IN), OPTIONAL :: LDMPA                 ! value of LMPA
LOGICAL,             INTENT(IN), OPTIONAL :: LDPGDFWR              ! LPGDFWR from yommse; write PGD fields
LOGICAL,             INTENT(IN), OPTIONAL :: LDHISFWR              ! LHISFWR from yommse; write historic & diagnostic data
INTEGER(KIND=JPIM),  INTENT(IN), OPTIONAL :: KSTOP
INTEGER(KIND=JPIM),  INTENT(IN), OPTIONAL :: KSTEP
INTEGER(KIND=JPIM),  INTENT(IN), OPTIONAL :: KFREQFS
INTEGER(KIND=JPIM),  INTENT(IN), OPTIONAL :: KYPROC
INTEGER(KIND=JPIM),  INTENT(IN), OPTIONAL :: KLON
REAL(KIND=JPRB),     INTENT(IN), OPTIONAL :: PINRAIN (:) ! (KLON)  ! liquid precipitation (kg/m2/s)
REAL(KIND=JPRB),     INTENT(IN), OPTIONAL :: PINSNOW (:) ! (KLON)  ! solid precipitation (kg/m2/s)
CHARACTER (LEN=*),   INTENT(IN), OPTIONAL :: CDFILE

#include "abor1.intfb.h"

CHARACTER (LEN=28) :: YFMFILE   ! name of the OUTPUT FM-file
CHARACTER (LEN=5)  :: YNUMBER   ! character string for the OUTPUT FM-file number
INTEGER            :: ILEN
INTEGER            :: IGRID          ! IGRID : grid indicator
INTEGER            :: ILENCH         ! ILENCH : length of comment string 
!
CHARACTER(LEN=16)  :: YRECFM         ! Name of the article to be written
CHARACTER(LEN=100) :: YCOMMENT       ! Comment string
CHARACTER(LEN=10)  :: YBIBUSER       ! user library
CHARACTER(LEN=28)  :: YDADNAME       ! dad name
INTEGER            :: IRESP, ININAR, IOUT, INB
INTEGER            :: IWORK
! Variables for FA files
INTEGER                            :: ITYPTR, ITRONC, INLATI, INXLON
INTEGER, DIMENSION(:), ALLOCATABLE :: INLOPA, INOZPA
REAL                               :: ZSLAPO, ZCLOPO, ZSLOPO, ZCODIL, ZREFER
REAL,    DIMENSION(:), ALLOCATABLE :: ZSINLA, ZAHYBR, ZBHYBR
REAL,    DIMENSION(:), ALLOCATABLE :: ZRAIN, ZSNOW
! 
LOGICAL                            :: LLPGDFWR   
LOGICAL                            :: LLHISFWR    
LOGICAL                            :: LLMPA
INTEGER (KIND=JPIM)                :: IYPROC
!
!-------------------------------------------------------------------------------

! Initialisation of writting counter
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('ARO_SURF_DIAG',0,ZHOOK_HANDLE)

LLMPA = .FALSE.
IF (PRESENT (LDMPA)) LLMPA = LDMPA

LLPGDFWR = .FALSE.
IF (PRESENT (LDPGDFWR)) LLPGDFWR = LDPGDFWR

LLHISFWR = .FALSE.
IF (PRESENT (LDHISFWR)) LLHISFWR = LDHISFWR

IF (PRESENT (KYPROC)) THEN
  IYPROC = KYPROC
ELSE
  IYPROC = MYPROC
ENDIF

DO NBLOCK = 1, NBLOCKTOT
  YSURFEX_ARO_CUR => YSURFEX_ARO_ALL (NBLOCK)
  CALL FLAG_UPDATE(YSURFEX_ARO_CUR%DLO,              &
                 & YSURFEX_ARO_CUR%DUO,                 &
                 & .FALSE.,.TRUE.,                      &
                 & YSURFEX_ARO_CUR%DUO%LPROVAR_TO_DIAG, &
                 & YSURFEX_ARO_CUR%DUO%LSELECT)
ENDDO

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IF KYPROC is the writting PROC (KYPROC=1), then
! open the MNH file
! write some basic REC
! close the MNH file
!
IF (IYPROC == 1) THEN
!
  IF (PRESENT (KFREQFS)) THEN
    IOUT=KFREQFS
    WRITE (YNUMBER,FMT="('.',I4.4)") IOUT
    YFMFILE=ADJUSTL('AROMOUT_'//YNUMBER)
  ELSEIF (PRESENT (CDFILE)) THEN
    YFMFILE = CDFILE
  ENDIF

  YRECFM='STORAGE_TYPE'
  YCOMMENT=' '
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  CALL SURFEX_FIELD_BUF_ADD (YSURFEX_CACHE_OUT, 'SU', YRECFM)

  YRECFM='MASDEV'
  YCOMMENT=' '
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  CALL SURFEX_FIELD_BUF_ADD (YSURFEX_CACHE_OUT, 47, YRECFM)
!
  YRECFM='BUGFIX'
  YCOMMENT=' '
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  CALL SURFEX_FIELD_BUF_ADD (YSURFEX_CACHE_OUT, 0, YRECFM)
!
  YRECFM='BIBUSER'
  YCOMMENT=' '
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  YBIBUSER = ' '
  CALL SURFEX_FIELD_BUF_ADD (YSURFEX_CACHE_OUT, YBIBUSER, YRECFM)
!
  YRECFM='PROGRAM'
  YCOMMENT=' '
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  YBIBUSER = ' '
  CALL SURFEX_FIELD_BUF_ADD (YSURFEX_CACHE_OUT, 'AROME', YRECFM)
!
  YRECFM='SURF'
  YCOMMENT=' '
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  CALL SURFEX_FIELD_BUF_ADD (YSURFEX_CACHE_OUT, 'EXTE', YRECFM)
!
  YRECFM='MY_NAME'
  YCOMMENT=' '
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  CALL SURFEX_FIELD_BUF_ADD (YSURFEX_CACHE_OUT, YFMFILE, YRECFM)
!
  YRECFM='DAD_NAME'
  YCOMMENT=' '
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  YDADNAME = ' '
  CALL SURFEX_FIELD_BUF_ADD (YSURFEX_CACHE_OUT, YDADNAME, YRECFM)
!
  YRECFM='L1D'
  YCOMMENT=' '
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  YDADNAME = ' '
  CALL SURFEX_FIELD_BUF_ADD (YSURFEX_CACHE_OUT, .FALSE., YRECFM)
!
  YRECFM='L2D'
  YCOMMENT=' '
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  YDADNAME = ' '
  CALL SURFEX_FIELD_BUF_ADD (YSURFEX_CACHE_OUT, .FALSE., YRECFM)
!
  YRECFM='PACK'
  YCOMMENT=' '
  IGRID=0
  ILENCH=LEN(YCOMMENT)
  YDADNAME = ' '
  CALL SURFEX_FIELD_BUF_ADD (YSURFEX_CACHE_OUT, .FALSE., YRECFM)
!
!
ENDIF
!
LFANOCOMPACT=.FALSE.
!

!* computes surface diagnostics
!
IF (LLMPA) THEN
  DO NBLOCK = 1,NBLOCKTOT
    YSURFEX_ARO_CUR => YSURFEX_ARO_ALL (NBLOCK)
    CALL DIAG_SURF_ATM_N(YSURFEX_ARO_CUR, &
                       & 'AROME ')
  ENDDO
ENDIF
!
!* writes the surface fields
!
CALL SURFEX_FIELD_BUF_SET_RECORD (YSURFEX_CACHE_OUT, .FALSE.)

! if NBLOCKTOT==0 the the task is fully in the extension zone. REK
IF (NBLOCKTOT > 0) THEN
DO NBLOCK = 0,NBLOCKTOT

  IWORK=MAX(1,NBLOCK)
  NINDX1=1+(IWORK-1)*NPROMA
  NINDX2=MIN(IWORK*NPROMA,MIN(NGPTOT,NGPTOT_CAP))
  NKPROMA=NINDX2-NINDX1+1

  IF (NBLOCK==0) THEN
    YSURFEX_ARO_CUR => YSURFEX_ARO_ALL (1)
  ELSE
    YSURFEX_ARO_CUR => YSURFEX_ARO_ALL (NBLOCK)
  ENDIF

  IF (LLHISFWR) &
    CALL WRITE_SURF_ATM_N(YSURFEX_ARO_CUR, 'AROME ','ALL',.FALSE.)

  IF (LLPGDFWR) &
    CALL WRITE_PGD_SURF_ATM_N(YSURFEX_ARO_CUR, 'AROME ')

  IF (LLHISFWR) &
    CALL WRITE_DIAG_SURF_ATM_N(YSURFEX_ARO_CUR, 'AROME ','ALL')

  IF (NBLOCK == 0) THEN
    CALL SURFEX_FIELD_BUF_PREALLOC (YSURFEX_CACHE_OUT)
    CALL SURFEX_FIELD_BUF_SET_RECORD (YSURFEX_CACHE_OUT, .TRUE.)
  ENDIF

ENDDO
ENDIF


LFIRST_WRITE = .FALSE.  !Set to FALSE after first writing step
NCPT_WRITE=0            !Set to 0 after each writing step

IF (LHOOK) CALL DR_HOOK('ARO_SURF_DIAG',1,ZHOOK_HANDLE)

END SUBROUTINE ARO_SURF_DIAG
