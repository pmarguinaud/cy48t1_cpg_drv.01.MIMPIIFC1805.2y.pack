!*
!     ------------------------------------------------------------------

!     This COMDECK includes the Thermodynamical functions for the cy39
!       ECMWF Physics package.
!       Consistent with YOMCST Basic physics constants, assuming the
!       partial pressure of water vapour is given by a first order
!       Taylor expansion of Qs(T) w.r.t. to Temperature, using constants
!       in YOETHF
!       Two sets of functions are available. In the first set only the
!       cases water or ice are distinguished by temperature.  This set 
!       consists of the functions FOEDELTA,FOEEW,FOEDE and FOELH.
!       The second set considers, besides the two cases water and ice 
!       also a mix of both for the temperature range _RTICE_ < T < _RTWAT_.
!       This set contains FOEALFA,FOEEWM,FOEDEM,FOELDCPM and FOELHM.
!       FKOOP modifies the ice saturation mixing ratio for homogeneous 
!       nucleation. FOE_DEWM_DT provides an approximate first derivative
!       of FOEEWM.

!       Depending on the consideration of mixed phases either the first 
!       set (e.g. surface, post-processing) or the second set 
!       (e.g. clouds, condensation, convection) should be used.

!     ------------------------------------------------------------------
!     *****************************************************************

!                NO CONSIDERATION OF MIXED PHASES

!     *****************************************************************
REAL(KIND=JPRB) :: FOEDELTA
REAL(KIND=JPRB) :: PTARE
FOEDELTA (PTARE) = MAX (0.0_JPRB,SIGN(1.0_JPRB,PTARE-_RTT_))

!                  FOEDELTA = 1    water
!                  FOEDELTA = 0    ice

!     THERMODYNAMICAL FUNCTIONS .

!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOEEW,FOEDE,FOEDESU,FOELH,FOELDCP
FOEEW ( PTARE ) = _R2ES_*EXP (&
  &(_R3LES_*FOEDELTA(PTARE)+_R3IES_*(1.0_JPRB-FOEDELTA(PTARE)))*(PTARE-_RTT_)&
&/ (PTARE-(_R4LES_*FOEDELTA(PTARE)+_R4IES_*(1.0_JPRB-FOEDELTA(PTARE)))))

FOEDE ( PTARE ) = &
  &(FOEDELTA(PTARE)*_R5ALVCP_+(1.0_JPRB-FOEDELTA(PTARE))*_R5ALSCP_)&
&/ (PTARE-(_R4LES_*FOEDELTA(PTARE)+_R4IES_*(1.0_JPRB-FOEDELTA(PTARE))))**2

FOEDESU ( PTARE ) = &
  &(FOEDELTA(PTARE)*_R5LES_+(1.0_JPRB-FOEDELTA(PTARE))*_R5IES_)&
&/ (PTARE-(_R4LES_*FOEDELTA(PTARE)+_R4IES_*(1.0_JPRB-FOEDELTA(PTARE))))**2

FOELH ( PTARE ) =&
         &FOEDELTA(PTARE)*_RLVTT_ + (1.0_JPRB-FOEDELTA(PTARE))*_RLSTT_

FOELDCP ( PTARE ) = &
         &FOEDELTA(PTARE)*_RALVDCP_ + (1.0_JPRB-FOEDELTA(PTARE))*_RALSDCP_

!     *****************************************************************

!           CONSIDERATION OF MIXED PHASES

!     *****************************************************************

!     FOEALFA is calculated to distinguish the three cases:

!                       FOEALFA=1            water phase
!                       FOEALFA=0            ice phase
!                       0 < FOEALFA < 1      mixed phase

!               INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOEALFA
FOEALFA (PTARE) = MIN(1.0_JPRB,((MAX(_RTICE_,MIN(_RTWAT_,PTARE))-_RTICE_)&
 &*_RTWAT_RTICE_R_)**2) 


!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOEEWM,FOEDEM,FOELDCPM,FOELHM,FOE_DEWM_DT
FOEEWM ( PTARE ) = _R2ES_ *&
     &(FOEALFA(PTARE)*EXP(_R3LES_*(PTARE-_RTT_)/(PTARE-_R4LES_))+&
  &(1.0_JPRB-FOEALFA(PTARE))*EXP(_R3IES_*(PTARE-_RTT_)/(PTARE-_R4IES_)))

FOE_DEWM_DT( PTARE ) = _R2ES_ * ( &
     & _R3LES_*FOEALFA(PTARE)*EXP(_R3LES_*(PTARE-_RTT_)/(PTARE-_R4LES_)) &
     &    *(_RTT_-_R4LES_)/(PTARE-_R4LES_)**2 + &
     & _R3IES_*(1.0-FOEALFA(PTARE))*EXP(_R3IES_*(PTARE-_RTT_)/(PTARE-_R4IES_)) &
     &    *(_RTT_-_R4IES_)/(PTARE-_R4IES_)**2)

FOEDEM ( PTARE ) = FOEALFA(PTARE)*_R5ALVCP_*(1.0_JPRB/(PTARE-_R4LES_)**2)+&
             &(1.0_JPRB-FOEALFA(PTARE))*_R5ALSCP_*(1.0_JPRB/(PTARE-_R4IES_)**2)

FOELDCPM ( PTARE ) = FOEALFA(PTARE)*_RALVDCP_+&
            &(1.0_JPRB-FOEALFA(PTARE))*_RALSDCP_

FOELHM ( PTARE ) =&
         &FOEALFA(PTARE)*_RLVTT_+(1.0_JPRB-FOEALFA(PTARE))*_RLSTT_


!     Temperature normalization for humidity background change of variable
!        INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOETB
FOETB ( PTARE )=FOEALFA(PTARE)*_R3LES_*(_RTT_-_R4LES_)*(1.0_JPRB/(PTARE-_R4LES_)**2)+&
             &(1.0_JPRB-FOEALFA(PTARE))*_R3IES_*(_RTT_-_R4IES_)*(1.0_JPRB/(PTARE-_R4IES_)**2)

!     ------------------------------------------------------------------
!     *****************************************************************

!           CONSIDERATION OF DIFFERENT MIXED PHASE FOR CONV

!     *****************************************************************

!     FOEALFCU is calculated to distinguish the three cases:

!                       FOEALFCU=1            water phase
!                       FOEALFCU=0            ice phase
!                       0 < FOEALFCU < 1      mixed phase

!               INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOEALFCU 
FOEALFCU (PTARE) = MIN(1.0_JPRB,((MAX(_RTICECU_,MIN(_RTWAT_,PTARE))&
&-_RTICECU_)*_RTWAT_RTICECU_R_)**2) 


!     Pressure of water vapour at saturation
!        INPUT : PTARE = TEMPERATURE
REAL(KIND=JPRB) :: FOEEWMCU,FOEDEMCU,FOELDCPMCU,FOELHMCU
FOEEWMCU ( PTARE ) = _R2ES_ *&
     &(FOEALFCU(PTARE)*EXP(_R3LES_*(PTARE-_RTT_)/(PTARE-_R4LES_))+&
  &(1.0_JPRB-FOEALFCU(PTARE))*EXP(_R3IES_*(PTARE-_RTT_)/(PTARE-_R4IES_)))

FOEDEMCU ( PTARE )=FOEALFCU(PTARE)*_R5ALVCP_*(1.0_JPRB/(PTARE-_R4LES_)**2)+&
             &(1.0_JPRB-FOEALFCU(PTARE))*_R5ALSCP_*(1.0_JPRB/(PTARE-_R4IES_)**2)

FOELDCPMCU ( PTARE ) = FOEALFCU(PTARE)*_RALVDCP_+&
            &(1.0_JPRB-FOEALFCU(PTARE))*_RALSDCP_

FOELHMCU ( PTARE ) =&
         &FOEALFCU(PTARE)*_RLVTT_+(1.0_JPRB-FOEALFCU(PTARE))*_RLSTT_
!     ------------------------------------------------------------------

!     Pressure of water vapour at saturation
!     This one is for the WMO definition of saturation, i.e. always
!     with respect to water.
!     
!     Duplicate to FOEELIQ and FOEEICE for separate ice variable
!     FOEELIQ always respect to water 
!     FOEEICE always respect to ice 
!     (could use FOEEW and FOEEWMO, but naming convention unclear)
!     FOELSON returns e wrt liquid water using D Sonntag (1994, Met. Zeit.)
!      - now recommended for use with radiosonde data (WMO CIMO guide, 2014)
!      unlike the FOEE functions does not include 1/(_RETV_+1.0_JPRB) factor

REAL(KIND=JPRB) :: FOEEWMO, FOEELIQ, FOEEICE, FOELSON 
FOEEWMO( PTARE ) = _R2ES_*EXP(_R3LES_*(PTARE-_RTT_)/(PTARE-_R4LES_))
FOEELIQ( PTARE ) = _R2ES_*EXP(_R3LES_*(PTARE-_RTT_)/(PTARE-_R4LES_))
FOEEICE( PTARE ) = _R2ES_*EXP(_R3IES_*(PTARE-_RTT_)/(PTARE-_R4IES_))
FOELSON( PTARE ) = EXP( -6096.9385_JPRB/PTARE + 21.2409642_JPRB &
	             - 2.711193E-2_JPRB * PTARE    &
                     + 1.673952E-5_JPRB * PTARE**2 &
		     + 2.433502_JPRB * LOG(PTARE))

REAL(KIND=JPRB) :: FOEEWM_V,FOEEWMCU_V,FOELES_V,FOEIES_V
REAL(KIND=JPRB) :: EXP1,EXP2
      FOELES_V(PTARE)=_R3LES_*(PTARE-_RTT_)/(PTARE-_R4LES_)
      FOEIES_V(PTARE)=_R3IES_*(PTARE-_RTT_)/(PTARE-_R4IES_)
      FOEEWM_V( PTARE,EXP1,EXP2 )=_R2ES_*(FOEALFA(PTARE)*EXP1+ &
          & (1.0_JPRB-FOEALFA(PTARE))*EXP2)
      FOEEWMCU_V ( PTARE,EXP1,EXP2 ) = _R2ES_*(FOEALFCU(PTARE)*EXP1+&
          &(1.0_JPRB-FOEALFCU(PTARE))*EXP2)

