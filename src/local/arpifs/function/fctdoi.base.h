!*   ------------------------------------------------------------------
!
!     FONCTION DE LA PARTITION EAU/GLACE
!     ICE/LIQUID WATER PARTITION FUNCTION
!
REAL(KIND=JPRB) :: PTARG1,P_RDTFAC, FONICE
!
FONICE ( PTARG1, P_RDTFAC ) = 1.0_JPRB - EXP ( - (_RTT_-MIN(_RTT_,PTARG1))**2 &
  & * (1.0_JPRB/(2.0_JPRB*(_RDT_*P_RDTFAC)**2)) )
!*
!     -----------------------------------------------------------------
