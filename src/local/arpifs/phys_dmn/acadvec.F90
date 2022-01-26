SUBROUTINE ACADVEC(YDPHY0,YDPHY2,KIDIA,KFDIA,KLON,KTDIA,KLEV,PDELP,PAPHI,PVAR,PVARW,PFVAR)

! ========================================================

!   Compute vertical flux of a variable with vertical speed PVARW

! ========================================================

!   Auteur: Yves Bouteloup, CNRM/GMAP

!   Date: 2010-12   

!     Modifications.
!     --------------
!         2012-06-22, J.M. Piriou: if LADVLIM=.TRUE. minimise p2 to the value such that 
!                     the transport scheme cannot transport beyond a distance of w*deltat, in a time-step.

! ========================================================

! ---------------
! INPUT VARIABLES
! ---------------

! KIDIA, : START OF HORIZONTAL LOOP 
! KFDIA  : END OF HORIZONTAL LOOP     
! KLON   : HORIZONTAL DIMENSION
! KTDIA : INDICE DE DEPART DES BOUCLES VERTICALES (1 EN GENERAL).
! KLEV  : VERTICAL DIMENSION
! PDELP  : LAYER THICKNESS IN PRESSURE UNITS.
! PAPHI  : GEOPOTENTIAL ON HALF-LEVELS.
! PVAR   : INPUT VARIABLE ON FULL LEVEL
! PVARW  : VERTICAL SPEED OF VARIABLE PVAR ON FULL LEVEL (positive from top to bottom)


! ---------------
! OUTPUT VARIABLE
! ---------------

! PFVAR  : FLUX OF VARIABLE PVAR ON HALF LEVEL (Positive from top to bottom)

! ========================================================

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE YOMPHY0   , ONLY : TPHY0
USE YOMPHY2   , ONLY : TPHY2
USE YOMCST    , ONLY : RG


IMPLICIT NONE

TYPE(TPHY0)       ,INTENT(IN)     :: YDPHY0
TYPE(TPHY2)       ,INTENT(IN)     :: YDPHY2
INTEGER(KIND=JPIM),INTENT(IN)     :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)     :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)     :: KLON 
INTEGER(KIND=JPIM), INTENT(IN)    :: KTDIA
INTEGER(KIND=JPIM),INTENT(IN)     :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)     :: PVAR (KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)     :: PVARW(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)     :: PDELP(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)     :: PAPHI(KLON,0:KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)    :: PFVAR(KLON,0:KLEV) 

REAL(KIND=JPRB) :: ZP1, ZP2(KLON), ZDZV(KLON), ZDZL(KLON), ZEPS

REAL(KIND=JPRB) :: ZFUP(KLON,0:KLEV), ZFDN(KLON,0:KLEV), ZALTIH(KLON,0:KLEV)       
REAL(KIND=JPRB) :: ZDPSG(KLON,KLEV)


INTEGER(KIND=JPIM) :: JLEV, JLON
REAL(KIND=JPRB) :: ZHOOK_HANDLE,ZLN_NEGLIG,ZUSCFL

! --------------------------------------------------------

!     CHECK RELIABILITY OF INPUT ARGUMENTS.

IF (LHOOK) CALL DR_HOOK('ACADVEC',0,ZHOOK_HANDLE)
ASSOCIATE(LADVLIM=>YDPHY0%LADVLIM, &
 & TSPHY=>YDPHY2%TSPHY)
! --------------------------------------------------------

DO JLEV = 0, KLEV
  ZFDN(:,JLEV) = 0.0_JPRB
  ZFUP(:,JLEV) = 0.0_JPRB
ENDDO

!- - - - - - - - - - - - - - -
IF (TSPHY > 0.0_JPRB) THEN
!- - - - - - - - - - - - - - -

! Some intialisation
  DO JLEV = KTDIA, KLEV
    DO JLON = KIDIA, KFDIA
      ZDPSG(JLON,JLEV) = PDELP(JLON,JLEV) / RG
    ENDDO        
  ENDDO  
  DO JLEV = KTDIA-1, KLEV
    DO JLON = KIDIA, KFDIA
      ZALTIH(JLON,JLEV) = PAPHI(JLON,JLEV) / RG 
    ENDDO
  ENDDO
  ZEPS = 1.E-20_JPRB
  ZLN_NEGLIG=LOG(1.0E-2_JPRB)


! First loop : From top to bottom

  DO JLEV = KTDIA, KLEV-1
    DO JLON = KIDIA, KFDIA
      ZDZV(JLON) = MAX(0._JPRB,-PVARW(JLON,JLEV))*TSPHY
      ZDZL(JLON) = ZALTIH(JLON,JLEV-1) - ZALTIH(JLON,JLEV)
    ENDDO
    IF(LADVLIM) THEN
      DO JLON = KIDIA, KFDIA
        ZUSCFL=ZDZL(JLON)/MAX(ZEPS,ZDZV(JLON))
        ZP2(JLON)  = MAX(0._JPRB,MIN(EXP(ZLN_NEGLIG*ZUSCFL),1._JPRB - ZUSCFL))
      ENDDO
    ELSE
      DO JLON = KIDIA, KFDIA
        ZP2(JLON)  = MAX(0._JPRB,1._JPRB - ZDZL(JLON)/MAX(ZEPS,ZDZV(JLON)))
      ENDDO
    ENDIF
    DO JLON = KIDIA, KFDIA
      ZP1  = MIN(1._JPRB , ZDZV(JLON)/ZDZL(JLON))
      ZFDN(JLON,JLEV) = ZP1*ZDPSG(JLON,JLEV)*PVAR(JLON,JLEV)/TSPHY&
        &+ ZP2(JLON)*ZFDN(JLON,JLEV-1)
    ENDDO
  ENDDO    

! Second loop : From bottom to top

  DO JLEV = KLEV, KTDIA+1, -1
    DO JLON = KIDIA, KFDIA
      ZDZV(JLON) = MAX(0._JPRB,PVARW(JLON,JLEV))*TSPHY
      ZDZL(JLON) = ZALTIH(JLON,JLEV-1) - ZALTIH(JLON,JLEV)
    ENDDO
    IF(LADVLIM) THEN
      DO JLON = KIDIA, KFDIA
        ZUSCFL=ZDZL(JLON)/MAX(ZEPS,ZDZV(JLON))
        ZP2(JLON)  = MAX(0._JPRB,MIN(EXP(ZLN_NEGLIG*ZUSCFL),1._JPRB - ZUSCFL))
      ENDDO
    ELSE
      DO JLON = KIDIA, KFDIA
        ZP2(JLON)  = MAX(0._JPRB,1._JPRB - ZDZL(JLON)/MAX(ZEPS,ZDZV(JLON)))
      ENDDO
    ENDIF
    DO JLON = KIDIA, KFDIA
      ZP1 = MIN(1._JPRB , ZDZV(JLON)/ZDZL(JLON))
      ZFUP(JLON,JLEV-1) = ZP1*ZDPSG(JLON,JLEV)*PVAR(JLON,JLEV)/TSPHY&
        &+ ZP2(JLON)*ZFUP(JLON,JLEV)
    ENDDO
  ENDDO    


!- - - - - - - - - - - - - - - - - - - - - - -
ENDIF  ! End of test on TSPHY > 0.0_JPRB
!- - - - - - - - - - - - - - - - - - - - - - -

!  Final loop construction of net flux

DO JLEV = KTDIA, KLEV-1
  DO JLON = KIDIA, KFDIA
    PFVAR(JLON,JLEV) = ZFDN(JLON,JLEV) - ZFUP(JLON,JLEV) 
  ENDDO
ENDDO   

! --------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ACADVEC',1,ZHOOK_HANDLE)
END SUBROUTINE ACADVEC
