! Rank and shape definitions for simple templating

MODULE VARIABLE_MODULE
  ! Base class definition of VARIABLE types that manages configuration
  ! and metadata for individual variables and associates them with the
  ! respective FIELDS objects that store the data.

USE PARKIND1, ONLY: JPIM, JPRB
USE FIELD_MODULE, ONLY: FIELD_2D, FIELD_3D, FIELD_4D

IMPLICIT NONE

TYPE, ABSTRACT :: VARIABLE_BASE
  ! Description and definition of a scientific variable that stores
  ! its associated data in one or more underlying fields, eg. for timestepping.

  ! Generic metadata like names and IDs
  CHARACTER(LEN=16)  :: NAME                  ! Primary name used for indexing
  CHARACTER(LEN=16)  :: CNAME     = ''        ! ARPEGE field name
  INTEGER(KIND=JPIM) :: IGRBCODE  = -999      ! GRIB code

  ! Flags that define the behaviour of the field variable
  LOGICAL            :: LACTIVE   = .FALSE.   ! Field in use
  LOGICAL            :: LT1       = .FALSE.   ! Field in t+dt GFL
  LOGICAL            :: LT9       = .FALSE.   ! Field in t-dt GFL
  LOGICAL            :: LPH9      = .FALSE.   ! Field in t-dt physics
  LOGICAL            :: LDL       = .FALSE.   ! Field has zontal derivative
  LOGICAL            :: LDM       = .FALSE.   ! Field has meridional derivative
  LOGICAL            :: LDL9      = .FALSE.   ! Field has zontal derivative at t-dt
  LOGICAL            :: LDM9      = .FALSE.   ! Field has meridional derivative at t-dt
  LOGICAL            :: LADV      = .FALSE.   ! Field advected or not
  LOGICAL            :: LWATER    = .FALSE.   
  REAL (KIND=JPRB)   :: RCP       = 0._JPRB

  ! TODO: Storage backend for Atlas (guard by #ifdef)

CONTAINS
  PROCEDURE(VARIABLE_BASE_FINAL), DEFERRED :: FINAL
END TYPE VARIABLE_BASE

ABSTRACT INTERFACE
  SUBROUTINE VARIABLE_BASE_FINAL(SELF)
    IMPORT :: VARIABLE_BASE
    CLASS(VARIABLE_BASE) :: SELF
  END SUBROUTINE VARIABLE_BASE_FINAL
END INTERFACE

TYPE, EXTENDS(VARIABLE_BASE) :: VARIABLE_2D
  ! TODO: Allocation-specific metadata, like shapes and dimensions
  ! Note that storing things like NLEV would break templating

  ! Array view pointers, to be set up from associated fields
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: P(:)  => NULL()  ! Basic field at t
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: T0(:) => NULL()  ! Basic field at t (alias of P)
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: T1(:) => NULL()  ! Basic field at t+dt
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: T9(:) => NULL()  ! Basic field at t-dt
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PH9(:)=> NULL()  ! Basic field for physics
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DL(:) => NULL()  ! Zonal derivative field
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DM(:) => NULL()  ! Meridional derivative field
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DL9(:) => NULL()  ! Zonal derivative field at t-dt
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DM9(:) => NULL()  ! Meridional derivative field at t-dt

  ! Pointers to associated FIELD objects
  TYPE(FIELD_2D), POINTER :: FT0 => NULL()  ! Basic field at t
  TYPE(FIELD_2D), POINTER :: FT1 => NULL()  ! Basic field at t+dt
  TYPE(FIELD_2D), POINTER :: FT9 => NULL()  ! Basic field at t-dt
  TYPE(FIELD_2D), POINTER :: FPH9 => NULL() ! Basic field for physics
  TYPE(FIELD_2D), POINTER :: FDL => NULL()  ! Zonal derivative field
  TYPE(FIELD_2D), POINTER :: FDM => NULL()  ! Meridional derivative field
  TYPE(FIELD_2D), POINTER :: FDL9 => NULL() ! Zonal derivative field at t-dt
  TYPE(FIELD_2D), POINTER :: FDM9 => NULL() ! Meridional derivative field at t-dt

CONTAINS
  PROCEDURE :: UPDATE_VIEW => VARIABLE_2D_UPDATE_VIEW
  PROCEDURE :: CLONE => VARIABLE_2D_CLONE
  PROCEDURE :: FINAL => VARIABLE_2D_FINAL
END TYPE VARIABLE_2D
TYPE, EXTENDS(VARIABLE_BASE) :: VARIABLE_3D
  ! TODO: Allocation-specific metadata, like shapes and dimensions
  ! Note that storing things like NLEV would break templating

  ! Array view pointers, to be set up from associated fields
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: P(:,:)  => NULL()  ! Basic field at t
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: T0(:,:) => NULL()  ! Basic field at t (alias of P)
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: T1(:,:) => NULL()  ! Basic field at t+dt
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: T9(:,:) => NULL()  ! Basic field at t-dt
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PH9(:,:)=> NULL()  ! Basic field for physics
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DL(:,:) => NULL()  ! Zonal derivative field
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DM(:,:) => NULL()  ! Meridional derivative field
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DL9(:,:) => NULL()  ! Zonal derivative field at t-dt
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DM9(:,:) => NULL()  ! Meridional derivative field at t-dt

  ! Pointers to associated FIELD objects
  TYPE(FIELD_3D), POINTER :: FT0 => NULL()  ! Basic field at t
  TYPE(FIELD_3D), POINTER :: FT1 => NULL()  ! Basic field at t+dt
  TYPE(FIELD_3D), POINTER :: FT9 => NULL()  ! Basic field at t-dt
  TYPE(FIELD_3D), POINTER :: FPH9 => NULL() ! Basic field for physics
  TYPE(FIELD_3D), POINTER :: FDL => NULL()  ! Zonal derivative field
  TYPE(FIELD_3D), POINTER :: FDM => NULL()  ! Meridional derivative field
  TYPE(FIELD_3D), POINTER :: FDL9 => NULL() ! Zonal derivative field at t-dt
  TYPE(FIELD_3D), POINTER :: FDM9 => NULL() ! Meridional derivative field at t-dt

CONTAINS
  PROCEDURE :: UPDATE_VIEW => VARIABLE_3D_UPDATE_VIEW
  PROCEDURE :: CLONE => VARIABLE_3D_CLONE
  PROCEDURE :: FINAL => VARIABLE_3D_FINAL
END TYPE VARIABLE_3D
TYPE, EXTENDS(VARIABLE_BASE) :: VARIABLE_4D
  ! TODO: Allocation-specific metadata, like shapes and dimensions
  ! Note that storing things like NLEV would break templating

  ! Array view pointers, to be set up from associated fields
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: P(:,:,:)  => NULL()  ! Basic field at t
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: T0(:,:,:) => NULL()  ! Basic field at t (alias of P)
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: T1(:,:,:) => NULL()  ! Basic field at t+dt
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: T9(:,:,:) => NULL()  ! Basic field at t-dt
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: PH9(:,:,:)=> NULL()  ! Basic field for physics
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DL(:,:,:) => NULL()  ! Zonal derivative field
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DM(:,:,:) => NULL()  ! Meridional derivative field
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DL9(:,:,:) => NULL()  ! Zonal derivative field at t-dt
  REAL(KIND=JPRB), POINTER, CONTIGUOUS :: DM9(:,:,:) => NULL()  ! Meridional derivative field at t-dt

  ! Pointers to associated FIELD objects
  TYPE(FIELD_4D), POINTER :: FT0 => NULL()  ! Basic field at t
  TYPE(FIELD_4D), POINTER :: FT1 => NULL()  ! Basic field at t+dt
  TYPE(FIELD_4D), POINTER :: FT9 => NULL()  ! Basic field at t-dt
  TYPE(FIELD_4D), POINTER :: FPH9 => NULL() ! Basic field for physics
  TYPE(FIELD_4D), POINTER :: FDL => NULL()  ! Zonal derivative field
  TYPE(FIELD_4D), POINTER :: FDM => NULL()  ! Meridional derivative field
  TYPE(FIELD_4D), POINTER :: FDL9 => NULL() ! Zonal derivative field at t-dt
  TYPE(FIELD_4D), POINTER :: FDM9 => NULL() ! Meridional derivative field at t-dt

CONTAINS
  PROCEDURE :: UPDATE_VIEW => VARIABLE_4D_UPDATE_VIEW
  PROCEDURE :: CLONE => VARIABLE_4D_CLONE
  PROCEDURE :: FINAL => VARIABLE_4D_FINAL
END TYPE VARIABLE_4D

INTERFACE VARIABLE_2D
  MODULE PROCEDURE :: VARIABLE_2D_INIT
  ! MODULE PROCEDURE :: VARIABLE_FROM_NAMELIST
END INTERFACE VARIABLE_2D
INTERFACE VARIABLE_3D
  MODULE PROCEDURE :: VARIABLE_3D_INIT
  ! MODULE PROCEDURE :: VARIABLE_FROM_NAMELIST
END INTERFACE VARIABLE_3D
INTERFACE VARIABLE_4D
  MODULE PROCEDURE :: VARIABLE_4D_INIT
  ! MODULE PROCEDURE :: VARIABLE_FROM_NAMELIST
END INTERFACE VARIABLE_4D

INTERFACE ARGUMENT_VALUE
  ! Helper interface to resolve values of optional arguments with defaults
  MODULE PROCEDURE ARGUMENT_VALUE_REAL
  MODULE PROCEDURE ARGUMENT_VALUE_INTEGER
  MODULE PROCEDURE ARGUMENT_VALUE_STRING
  MODULE PROCEDURE ARGUMENT_VALUE_LOGICAL
END INTERFACE ARGUMENT_VALUE

CONTAINS

  FUNCTION ARGUMENT_VALUE_REAL(ARG, DEFAULT) RESULT(VAL)
    ! Helper function to resolve value for optional argument
    REAL(KIND=JPRB), OPTIONAL, INTENT(IN) :: ARG
    REAL(KIND=JPRB), INTENT(IN) :: DEFAULT
    REAL(KIND=JPRB) :: VAL

    IF (PRESENT(ARG)) THEN
      VAL = ARG
    ELSE
      VAL = DEFAULT
    END IF
  END FUNCTION ARGUMENT_VALUE_REAL

  FUNCTION ARGUMENT_VALUE_INTEGER(ARG, DEFAULT) RESULT(VAL)
    ! Helper function to resolve value for optional argument
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: ARG
    INTEGER(KIND=JPIM), INTENT(IN) :: DEFAULT
    INTEGER(KIND=JPIM) :: VAL

    IF (PRESENT(ARG)) THEN
      VAL = ARG
    ELSE
      VAL = DEFAULT
    END IF
  END FUNCTION ARGUMENT_VALUE_INTEGER

  FUNCTION ARGUMENT_VALUE_STRING(ARG, DEFAULT) RESULT(VAL)
    ! Helper function to resolve value for optional argument
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: ARG
    CHARACTER(LEN=*), INTENT(IN) :: DEFAULT
    CHARACTER(:), ALLOCATABLE :: VAL

    IF (PRESENT(ARG)) THEN
      ALLOCATE(VAL, SOURCE=ARG)
    ELSE
      ALLOCATE(VAL, SOURCE=DEFAULT)
    END IF
  END FUNCTION ARGUMENT_VALUE_STRING

  FUNCTION ARGUMENT_VALUE_LOGICAL(ARG, DEFAULT) RESULT(VAL)
    ! Helper function to resolve value for optional argument
    LOGICAL, OPTIONAL, INTENT(IN) :: ARG
    LOGICAL, INTENT(IN) :: DEFAULT
    LOGICAL :: VAL

    IF (PRESENT(ARG)) THEN
      VAL = ARG
    ELSE
      VAL = DEFAULT
    END IF
  END FUNCTION ARGUMENT_VALUE_LOGICAL


  FUNCTION VARIABLE_2D_INIT(NAME, CNAME, IGRBCODE, LACTIVE, LADV, LT1, LT9, LPH9, LDL, LDM, LDL9, LDM9, LWATER, RCP) RESULT(SELF)
    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_vALUE, IEEE_QUIET_NAN
    ! Templated constructor that creates new instances
    TYPE(VARIABLE_2D) :: SELF
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: CNAME
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: IGRBCODE
    LOGICAL, OPTIONAL, INTENT(IN) :: LACTIVE, LADV, LT1, LT9, LPH9, LDL, LDM, LDL9, LDM9, LWATER
    REAL(KIND=JPRB), OPTIONAL, INTENT (IN) :: RCP

    SELF%NAME = NAME
    SELF%CNAME = ARGUMENT_VALUE(ARG=CNAME, DEFAULT=NAME)
    SELF%IGRBCODE = ARGUMENT_VALUE(ARG=IGRBCODE, DEFAULT=-999)
    SELF%LACTIVE = ARGUMENT_VALUE(ARG=LACTIVE, DEFAULT=.FALSE.)
    SELF%LADV = ARGUMENT_VALUE(ARG=LADV, DEFAULT=.FALSE.)
    SELF%LT1 = ARGUMENT_VALUE(ARG=LT1, DEFAULT=.FALSE.)
    SELF%LT9 = ARGUMENT_VALUE(ARG=LT9, DEFAULT=.FALSE.)
    SELF%LPH9 = ARGUMENT_VALUE(ARG=LPH9, DEFAULT=.FALSE.)
    SELF%LDL = ARGUMENT_VALUE(ARG=LDL, DEFAULT=.FALSE.)
    SELF%LDM = ARGUMENT_VALUE(ARG=LDM, DEFAULT=.FALSE.)
    SELF%LDL9 = ARGUMENT_VALUE(ARG=LDL9, DEFAULT=.FALSE.)
    SELF%LDM9 = ARGUMENT_VALUE(ARG=LDM9, DEFAULT=.FALSE.)
    SELF%LWATER = ARGUMENT_VALUE(ARG=LWATER, DEFAULT=.FALSE.)
    SELF%RCP = ARGUMENT_VALUE(ARG=RCP, DEFAULT=IEEE_VALUE (RCP, IEEE_QUIET_NAN))
  END FUNCTION VARIABLE_2D_INIT
  FUNCTION VARIABLE_3D_INIT(NAME, CNAME, IGRBCODE, LACTIVE, LADV, LT1, LT9, LPH9, LDL, LDM, LDL9, LDM9, LWATER, RCP) RESULT(SELF)
    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_vALUE, IEEE_QUIET_NAN
    ! Templated constructor that creates new instances
    TYPE(VARIABLE_3D) :: SELF
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: CNAME
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: IGRBCODE
    LOGICAL, OPTIONAL, INTENT(IN) :: LACTIVE, LADV, LT1, LT9, LPH9, LDL, LDM, LDL9, LDM9, LWATER
    REAL(KIND=JPRB), OPTIONAL, INTENT (IN) :: RCP

    SELF%NAME = NAME
    SELF%CNAME = ARGUMENT_VALUE(ARG=CNAME, DEFAULT=NAME)
    SELF%IGRBCODE = ARGUMENT_VALUE(ARG=IGRBCODE, DEFAULT=-999)
    SELF%LACTIVE = ARGUMENT_VALUE(ARG=LACTIVE, DEFAULT=.FALSE.)
    SELF%LADV = ARGUMENT_VALUE(ARG=LADV, DEFAULT=.FALSE.)
    SELF%LT1 = ARGUMENT_VALUE(ARG=LT1, DEFAULT=.FALSE.)
    SELF%LT9 = ARGUMENT_VALUE(ARG=LT9, DEFAULT=.FALSE.)
    SELF%LPH9 = ARGUMENT_VALUE(ARG=LPH9, DEFAULT=.FALSE.)
    SELF%LDL = ARGUMENT_VALUE(ARG=LDL, DEFAULT=.FALSE.)
    SELF%LDM = ARGUMENT_VALUE(ARG=LDM, DEFAULT=.FALSE.)
    SELF%LDL9 = ARGUMENT_VALUE(ARG=LDL9, DEFAULT=.FALSE.)
    SELF%LDM9 = ARGUMENT_VALUE(ARG=LDM9, DEFAULT=.FALSE.)
    SELF%LWATER = ARGUMENT_VALUE(ARG=LWATER, DEFAULT=.FALSE.)
    SELF%RCP = ARGUMENT_VALUE(ARG=RCP, DEFAULT=IEEE_VALUE (RCP, IEEE_QUIET_NAN))
  END FUNCTION VARIABLE_3D_INIT
  FUNCTION VARIABLE_4D_INIT(NAME, CNAME, IGRBCODE, LACTIVE, LADV, LT1, LT9, LPH9, LDL, LDM, LDL9, LDM9, LWATER, RCP) RESULT(SELF)
    USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_vALUE, IEEE_QUIET_NAN
    ! Templated constructor that creates new instances
    TYPE(VARIABLE_4D) :: SELF
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: CNAME
    INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: IGRBCODE
    LOGICAL, OPTIONAL, INTENT(IN) :: LACTIVE, LADV, LT1, LT9, LPH9, LDL, LDM, LDL9, LDM9, LWATER
    REAL(KIND=JPRB), OPTIONAL, INTENT (IN) :: RCP

    SELF%NAME = NAME
    SELF%CNAME = ARGUMENT_VALUE(ARG=CNAME, DEFAULT=NAME)
    SELF%IGRBCODE = ARGUMENT_VALUE(ARG=IGRBCODE, DEFAULT=-999)
    SELF%LACTIVE = ARGUMENT_VALUE(ARG=LACTIVE, DEFAULT=.FALSE.)
    SELF%LADV = ARGUMENT_VALUE(ARG=LADV, DEFAULT=.FALSE.)
    SELF%LT1 = ARGUMENT_VALUE(ARG=LT1, DEFAULT=.FALSE.)
    SELF%LT9 = ARGUMENT_VALUE(ARG=LT9, DEFAULT=.FALSE.)
    SELF%LPH9 = ARGUMENT_VALUE(ARG=LPH9, DEFAULT=.FALSE.)
    SELF%LDL = ARGUMENT_VALUE(ARG=LDL, DEFAULT=.FALSE.)
    SELF%LDM = ARGUMENT_VALUE(ARG=LDM, DEFAULT=.FALSE.)
    SELF%LDL9 = ARGUMENT_VALUE(ARG=LDL9, DEFAULT=.FALSE.)
    SELF%LDM9 = ARGUMENT_VALUE(ARG=LDM9, DEFAULT=.FALSE.)
    SELF%LWATER = ARGUMENT_VALUE(ARG=LWATER, DEFAULT=.FALSE.)
    SELF%RCP = ARGUMENT_VALUE(ARG=RCP, DEFAULT=IEEE_VALUE (RCP, IEEE_QUIET_NAN))
  END FUNCTION VARIABLE_4D_INIT

  FUNCTION VARIABLE_2D_CLONE(SELF) RESULT(NEWOBJ)
    ! Clone (deep-copy) this VARIABLE object to replicate associated FIELD objects
    !
    ! This is required create per-thread replication of the data view pointers
    ! under the fields associated with this VARIABLE.
    CLASS(VARIABLE_2D) :: SELF
    TYPE(VARIABLE_2D) :: NEWOBJ

    NEWOBJ = SELF
    NEWOBJ%NAME = SELF%NAME
    NEWOBJ%CNAME = SELF%CNAME
    IF (ASSOCIATED(SELF%FT0))  NEWOBJ%FT0 => SELF%FT0%CLONE()
    IF (ASSOCIATED(SELF%FT1))  NEWOBJ%FT1 => SELF%FT1%CLONE()
    IF (ASSOCIATED(SELF%FT9))  NEWOBJ%FT9 => SELF%FT9%CLONE()
    IF (ASSOCIATED(SELF%FPH9)) NEWOBJ%FPH9=> SELF%FPH9%CLONE()
    IF (ASSOCIATED(SELF%FDL))  NEWOBJ%FDL => SELF%FDL%CLONE()
    IF (ASSOCIATED(SELF%FDM))  NEWOBJ%FDM => SELF%FDM%CLONE()
    IF (ASSOCIATED(SELF%FDL9)) NEWOBJ%FDL9=> SELF%FDL9%CLONE()
    IF (ASSOCIATED(SELF%FDM9)) NEWOBJ%FDM9=> SELF%FDM9%CLONE()
  END FUNCTION VARIABLE_2D_CLONE
  FUNCTION VARIABLE_3D_CLONE(SELF) RESULT(NEWOBJ)
    ! Clone (deep-copy) this VARIABLE object to replicate associated FIELD objects
    !
    ! This is required create per-thread replication of the data view pointers
    ! under the fields associated with this VARIABLE.
    CLASS(VARIABLE_3D) :: SELF
    TYPE(VARIABLE_3D) :: NEWOBJ

    NEWOBJ = SELF
    NEWOBJ%NAME = SELF%NAME
    NEWOBJ%CNAME = SELF%CNAME
    IF (ASSOCIATED(SELF%FT0))  NEWOBJ%FT0 => SELF%FT0%CLONE()
    IF (ASSOCIATED(SELF%FT1))  NEWOBJ%FT1 => SELF%FT1%CLONE()
    IF (ASSOCIATED(SELF%FT9))  NEWOBJ%FT9 => SELF%FT9%CLONE()
    IF (ASSOCIATED(SELF%FPH9)) NEWOBJ%FPH9=> SELF%FPH9%CLONE()
    IF (ASSOCIATED(SELF%FDL))  NEWOBJ%FDL => SELF%FDL%CLONE()
    IF (ASSOCIATED(SELF%FDM))  NEWOBJ%FDM => SELF%FDM%CLONE()
    IF (ASSOCIATED(SELF%FDL9)) NEWOBJ%FDL9=> SELF%FDL9%CLONE()
    IF (ASSOCIATED(SELF%FDM9)) NEWOBJ%FDM9=> SELF%FDM9%CLONE()
  END FUNCTION VARIABLE_3D_CLONE
  FUNCTION VARIABLE_4D_CLONE(SELF) RESULT(NEWOBJ)
    ! Clone (deep-copy) this VARIABLE object to replicate associated FIELD objects
    !
    ! This is required create per-thread replication of the data view pointers
    ! under the fields associated with this VARIABLE.
    CLASS(VARIABLE_4D) :: SELF
    TYPE(VARIABLE_4D) :: NEWOBJ

    NEWOBJ = SELF
    NEWOBJ%NAME = SELF%NAME
    NEWOBJ%CNAME = SELF%CNAME
    IF (ASSOCIATED(SELF%FT0))  NEWOBJ%FT0 => SELF%FT0%CLONE()
    IF (ASSOCIATED(SELF%FT1))  NEWOBJ%FT1 => SELF%FT1%CLONE()
    IF (ASSOCIATED(SELF%FT9))  NEWOBJ%FT9 => SELF%FT9%CLONE()
    IF (ASSOCIATED(SELF%FPH9)) NEWOBJ%FPH9=> SELF%FPH9%CLONE()
    IF (ASSOCIATED(SELF%FDL))  NEWOBJ%FDL => SELF%FDL%CLONE()
    IF (ASSOCIATED(SELF%FDM))  NEWOBJ%FDM => SELF%FDM%CLONE()
    IF (ASSOCIATED(SELF%FDL9)) NEWOBJ%FDL9=> SELF%FDL9%CLONE()
    IF (ASSOCIATED(SELF%FDM9)) NEWOBJ%FDM9=> SELF%FDM9%CLONE()
  END FUNCTION VARIABLE_4D_CLONE

  SUBROUTINE VARIABLE_2D_UPDATE_VIEW(SELF, BLOCK_INDEX)
    ! Update the internal data view pointers of all associated fields
    CLASS(VARIABLE_2D) :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    REAL(KIND=JPRB), TARGET, SAVE :: ZDUM(1)

    ! Set on-object data view pointers from storage FIELDS
    IF (ASSOCIATED(SELF%FT0))  THEN
      SELF%T0 => SELF%FT0%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%T0 => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FT1))  THEN
      SELF%T1 => SELF%FT1%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%T1 => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FT9))  THEN
      SELF%T9 => SELF%FT9%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%T9 => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FPH9)) THEN
      SELF%PH9=> SELF%FPH9%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%PH9 => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FDL))  THEN
      SELF%DL => SELF%FDL%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%DL => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FDM))  THEN
      SELF%DM => SELF%FDM%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%DM => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FDL9)) THEN
      SELF%DL9=> SELF%FDL9%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%DL9 => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FDM9)) THEN
      SELF%DM9=> SELF%FDM9%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%DM9 => ZDUM
    ENDIF
    SELF%P => SELF%T0  ! Alias T0 pointer
  END SUBROUTINE VARIABLE_2D_UPDATE_VIEW
  SUBROUTINE VARIABLE_3D_UPDATE_VIEW(SELF, BLOCK_INDEX)
    ! Update the internal data view pointers of all associated fields
    CLASS(VARIABLE_3D) :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    REAL(KIND=JPRB), TARGET, SAVE :: ZDUM(1, 1)

    ! Set on-object data view pointers from storage FIELDS
    IF (ASSOCIATED(SELF%FT0))  THEN
      SELF%T0 => SELF%FT0%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%T0 => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FT1))  THEN
      SELF%T1 => SELF%FT1%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%T1 => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FT9))  THEN
      SELF%T9 => SELF%FT9%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%T9 => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FPH9)) THEN
      SELF%PH9=> SELF%FPH9%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%PH9 => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FDL))  THEN
      SELF%DL => SELF%FDL%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%DL => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FDM))  THEN
      SELF%DM => SELF%FDM%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%DM => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FDL9)) THEN
      SELF%DL9=> SELF%FDL9%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%DL9 => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FDM9)) THEN
      SELF%DM9=> SELF%FDM9%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%DM9 => ZDUM
    ENDIF
    SELF%P => SELF%T0  ! Alias T0 pointer
  END SUBROUTINE VARIABLE_3D_UPDATE_VIEW
  SUBROUTINE VARIABLE_4D_UPDATE_VIEW(SELF, BLOCK_INDEX)
    ! Update the internal data view pointers of all associated fields
    CLASS(VARIABLE_4D) :: SELF
    INTEGER(KIND=JPIM), INTENT(IN) :: BLOCK_INDEX
    REAL(KIND=JPRB), TARGET, SAVE :: ZDUM(1, 1, 1)

    ! Set on-object data view pointers from storage FIELDS
    IF (ASSOCIATED(SELF%FT0))  THEN
      SELF%T0 => SELF%FT0%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%T0 => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FT1))  THEN
      SELF%T1 => SELF%FT1%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%T1 => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FT9))  THEN
      SELF%T9 => SELF%FT9%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%T9 => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FPH9)) THEN
      SELF%PH9=> SELF%FPH9%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%PH9 => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FDL))  THEN
      SELF%DL => SELF%FDL%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%DL => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FDM))  THEN
      SELF%DM => SELF%FDM%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%DM => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FDL9)) THEN
      SELF%DL9=> SELF%FDL9%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%DL9 => ZDUM
    ENDIF
    IF (ASSOCIATED(SELF%FDM9)) THEN
      SELF%DM9=> SELF%FDM9%GET_VIEW(BLOCK_INDEX)
    ELSE
      SELF%DM9 => ZDUM
    ENDIF
    SELF%P => SELF%T0  ! Alias T0 pointer
  END SUBROUTINE VARIABLE_4D_UPDATE_VIEW

  SUBROUTINE VARIABLE_2D_FINAL(SELF)
    ! Templated destructor that cleans up an object instance
    CLASS(VARIABLE_2D) :: SELF

    IF (ASSOCIATED(SELF%FT0)) THEN
      CALL SELF%FT0%FINAL()
      DEALLOCATE(SELF%FT0)
    END IF
    IF (ASSOCIATED(SELF%FT1)) THEN
      CALL SELF%FT1%FINAL()
      DEALLOCATE(SELF%FT1)
    END IF
    IF (ASSOCIATED(SELF%FT9)) THEN
      CALL SELF%FT9%FINAL()
      DEALLOCATE(SELF%FT9)
    END IF
    IF (ASSOCIATED(SELF%FPH9)) THEN
      CALL SELF%FPH9%FINAL()
      DEALLOCATE(SELF%FPH9)
    END IF
    IF (ASSOCIATED(SELF%FDL)) THEN
      CALL SELF%FDL%FINAL()
      DEALLOCATE(SELF%FDL)
    END IF
    IF (ASSOCIATED(SELF%FDM)) THEN
      CALL SELF%FDM%FINAL()
      DEALLOCATE(SELF%FDM)
    END IF
    IF (ASSOCIATED(SELF%FDL9)) THEN
      CALL SELF%FDL9%FINAL()
      DEALLOCATE(SELF%FDL9)
    END IF
    IF (ASSOCIATED(SELF%FDM9)) THEN
      CALL SELF%FDM9%FINAL()
      DEALLOCATE(SELF%FDM9)
    END IF
  END SUBROUTINE VARIABLE_2D_FINAL
  SUBROUTINE VARIABLE_3D_FINAL(SELF)
    ! Templated destructor that cleans up an object instance
    CLASS(VARIABLE_3D) :: SELF

    IF (ASSOCIATED(SELF%FT0)) THEN
      CALL SELF%FT0%FINAL()
      DEALLOCATE(SELF%FT0)
    END IF
    IF (ASSOCIATED(SELF%FT1)) THEN
      CALL SELF%FT1%FINAL()
      DEALLOCATE(SELF%FT1)
    END IF
    IF (ASSOCIATED(SELF%FT9)) THEN
      CALL SELF%FT9%FINAL()
      DEALLOCATE(SELF%FT9)
    END IF
    IF (ASSOCIATED(SELF%FPH9)) THEN
      CALL SELF%FPH9%FINAL()
      DEALLOCATE(SELF%FPH9)
    END IF
    IF (ASSOCIATED(SELF%FDL)) THEN
      CALL SELF%FDL%FINAL()
      DEALLOCATE(SELF%FDL)
    END IF
    IF (ASSOCIATED(SELF%FDM)) THEN
      CALL SELF%FDM%FINAL()
      DEALLOCATE(SELF%FDM)
    END IF
    IF (ASSOCIATED(SELF%FDL9)) THEN
      CALL SELF%FDL9%FINAL()
      DEALLOCATE(SELF%FDL9)
    END IF
    IF (ASSOCIATED(SELF%FDM9)) THEN
      CALL SELF%FDM9%FINAL()
      DEALLOCATE(SELF%FDM9)
    END IF
  END SUBROUTINE VARIABLE_3D_FINAL
  SUBROUTINE VARIABLE_4D_FINAL(SELF)
    ! Templated destructor that cleans up an object instance
    CLASS(VARIABLE_4D) :: SELF

    IF (ASSOCIATED(SELF%FT0)) THEN
      CALL SELF%FT0%FINAL()
      DEALLOCATE(SELF%FT0)
    END IF
    IF (ASSOCIATED(SELF%FT1)) THEN
      CALL SELF%FT1%FINAL()
      DEALLOCATE(SELF%FT1)
    END IF
    IF (ASSOCIATED(SELF%FT9)) THEN
      CALL SELF%FT9%FINAL()
      DEALLOCATE(SELF%FT9)
    END IF
    IF (ASSOCIATED(SELF%FPH9)) THEN
      CALL SELF%FPH9%FINAL()
      DEALLOCATE(SELF%FPH9)
    END IF
    IF (ASSOCIATED(SELF%FDL)) THEN
      CALL SELF%FDL%FINAL()
      DEALLOCATE(SELF%FDL)
    END IF
    IF (ASSOCIATED(SELF%FDM)) THEN
      CALL SELF%FDM%FINAL()
      DEALLOCATE(SELF%FDM)
    END IF
    IF (ASSOCIATED(SELF%FDL9)) THEN
      CALL SELF%FDL9%FINAL()
      DEALLOCATE(SELF%FDL9)
    END IF
    IF (ASSOCIATED(SELF%FDM9)) THEN
      CALL SELF%FDM9%FINAL()
      DEALLOCATE(SELF%FDM9)
    END IF
  END SUBROUTINE VARIABLE_4D_FINAL

END MODULE VARIABLE_MODULE
