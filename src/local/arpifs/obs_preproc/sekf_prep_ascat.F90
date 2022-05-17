SUBROUTINE SEKF_PREP_ASCAT(YDGEOMETRY,YDRIP)

!**** *SEKF_PREP_ASCAT*  - match ASCAT soil moisture 'observations' with the corresponding model field

!     Purpose.
!     --------
!     (1) extract ASCAT soil moisture observations from ODB and collocate with model time steps / grid boxes
!     (2) calculate first guess departures and store them in ODB


!**   Interface.
!     ----------
!        *CALL* *SEKF_PREP_ASCAT

!     Explicit arguments :           None
!     --------------------

!     Implicit arguments :      None
!     --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        K. Sciapl  *ECMWF*

!     Modifications.
!     --------------
!     2008-08-03 - K.Scipal - new sql query
!     G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM and TCSGLEG
!     T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     K. Yessad (July 2014): Move some variables.
!     P. de Rosnay (November 2019): introduce ASCAT-C for monitoring only
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK
USE YOMCT0   , ONLY : LSCREEN
USE YOMRIP0  , ONLY : NSSSSS, NINDAT
USE YOMRIP   , ONLY : TRIP
USE YOMLUN   , ONLY : NULOUT
USE YOMDB
USE YOMANCS  , ONLY : RMDI, NMDI, RADIANS, RDEGREES
USE YOMCST   , ONLY : RPI
USE YOMCOCTP , ONLY : NSCAT3
USE YOMSEKF  , ONLY : TYPE_ASCAT_OBS, YLASCATOBS, VOL_SM_FG, NOBS_ASCAT, NOBS_ASCAT_TOT
USE ECSORT_MIX,ONLY : KEYSORT
USE YOMMP0   , ONLY : MYPROC

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN)   :: YDGEOMETRY
TYPE(TRIP)     ,INTENT(INOUT):: YDRIP
INTEGER(KIND=JPIM)              :: JLENG, JBLOCK, JL, JKGLO
INTEGER(KIND=JPIM)              :: JTST, IMODTIME
INTEGER(KIND=JPIM)              :: INFO(1), ILEN, IKINFOLEN
INTEGER(KIND=JPIM)              :: IASCAT, IOBS, IOBS_DAT, IOBS_ETM, IRET
INTEGER(KIND=JPIM)              :: JALL, IND, II, JJ, IGP
INTEGER(KIND=JPIM)              :: IOBST, IOBST_HH, IOBST_MN, IOBST_DD, IOBST_MM, IOBST_YY, IDIFF
INTEGER(KIND=JPIM)              :: IMODT_DD, IMODT_YY, IMODT_MM, IDOFF
INTEGER(KIND=JPIM)              :: IGP_N, ISUNIQ, IS, IGUNIQ, IG

INTEGER(KIND=JPIM), ALLOCATABLE :: IKINFO(:), INDGP(:,:), IKEY(:,:)

REAL(KIND=JPRB)                 :: ZHOOK_HANDLE
REAL(KIND=JPRB)                 :: ZLAT0, ZLATI0, ZLATE
REAL(KIND=JPRB)                 :: ZOBSLAT, ZOBSLON, ZGPLON, ZGPLAT
REAL(KIND=JPRB)                 :: ZFGDEP


TYPE(TYPE_ASCAT_OBS), ALLOCATABLE :: YLASCATOBSTMP(:), YLASCATOBSTMP2(:)

!     ------------------------------------------------------------------

#include "getdb.intfb.h"
#include "putdb.intfb.h"

#include "openmp_obs.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SEKF_PREP_ASCAT',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDGEM=>YDGEOMETRY%YRGEM, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDCSGLEG=>YDGEOMETRY%YRCSGLEG&
& )
ASSOCIATE(NDGLG=>YDDIM%NDGLG, NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA,   NGPTOT=>YDGEM%NGPTOT, NGPTOTG=>YDGEM%NGPTOTG, &
& NLOENG=>YDGEM%NLOENG,   NSTART=>YDRIP%NSTART, NSTOP=>YDRIP%NSTOP, TSTEP=>YDRIP%TSTEP)
!     ------------------------------------------------------------------

WRITE(NULOUT,*) 'START SEKF_PREP_ASCAT'


!*     0.0 Set up
!*

ALLOCATE(NOBS_ASCAT(NPROMA,NGPBLKS))
NOBS_ASCAT(:,:) = 0


!*     1.0 Set up lookuptable for observation/model collocation
!*
!*        A global lookup tabel is generated and for each grid
!*        point processed on MYPROC, the PROMA/BLOCK index is stored.
!*        The rest is filled with missing value indicators


ALLOCATE(INDGP(NGPTOTG,2))
INDGP(:,:) = NMDI

ZLATE=ASIN(YDCSGLEG%RMU(1))
ZLAT0=ZLATE*RDEGREES*1._JPRB
ZLATI0=-2._JPRB*ZLAT0/REAL(NDGLG-1)

!WRITE(NULOUT,*) 'CONS: ', RDEGREES, RADIANS, RPI
!WRITE(NULOUT,*) 'NLOENG: ', NLOENG
!WRITE(NULOUT,*) 'ZLATE: ', ZLATE, ZLAT0, ZLATI0, NLOENG(1)


DO JKGLO=1,NGPTOT,NPROMA

   JBLOCK=(JKGLO-1)/NPROMA + 1
   JLENG=MIN(NPROMA,NGPTOT-JKGLO+1)

   DO JL=1,JLENG

      ZGPLON=YDGSGEOM_NB%GELAM(JKGLO+JL-1)*RDEGREES
      ZGPLAT=YDGSGEOM_NB%GELAT(JKGLO+JL-1)*RDEGREES
      IF (ZGPLON<0.) ZGPLON=ZGPLON+360.

      II=NINT((ZGPLAT-ZLAT0)/ZLATI0+1._JPRB)
      JJ=NINT((ZGPLON*RADIANS*FLOAT(NLOENG(II)))/(RPI*2._JPRB)+1._JPRB)
      IGP=SUM(NLOENG(1:II-1))+JJ

      INDGP(IGP,1)=JL
      INDGP(IGP,2)=JBLOCK

   ENDDO
ENDDO


!*     2.0 Read ASCAT soil moisture observations from ODB
!*

WRITE(NULOUT,*) 'START SEKF_PREP_ASCAT ODB ACCESS'


!* set variables required for ODB access
!* LSCREEN temporarily has to be set to true, otherwise data can not be read
!* will be set to false again, after closing ODB

LSCREEN=.TRUE.
ILEN = 0
IKINFOLEN=2
ALLOCATE(IKINFO(IKINFOLEN))


CALL GETDB('ASCATSM',1,ILEN,INFO,0, -1, -1, -1, -1, NSCAT3, -1)


IASCAT=0

IF (ILEN>0) THEN

   ALLOCATE(YLASCATOBSTMP(ILEN))

   DO IOBS = 1, ILEN

      !* 1. set the link header - body (containing the observations)
      !* 2. we only need transformed soil moisture so look for variablenumber 180
      !* 3. avoid missing data (i.e. not over land, bad data, snow/freezing etc.)


      !* initialise observation entry with missing value
      !* later we will dump missing observations
      !* observations are missing if they aer
      !* 1. invalid
      !* 2. out of the time window
      !* 3. processed on a different processor

      YLASCATOBSTMP(IOBS)%ZOBS = RMDI



!         WRITE(NULOUT,*) ' SEKF_PREP_ASCAT REPORTYPE',ROBHDR(IOBS,MDB_REPORTYPE_AT_HDR)
      DO JALL = MLNKH2B(IOBS),MLNKH2B(IOBS+1)-1
         IF (ROBODY(JALL,MDBVNM) == 180 .AND.  ROBHDR(IOBS,MDB_REPORTYPE_AT_HDR) /= 9011 ) THEN
            IF (ROBODY(JALL,MDBVAR) /= RMDI)  THEN

               !* get the location for each observation into the red. gaussian grid
               !* for the time being nearest neighbour is used

               ZOBSLAT   = ROBHDR(IOBS,MDBLAT) * RDEGREES
               ZOBSLON   = ROBHDR(IOBS,MDBLON) * RDEGREES
               IF (ZOBSLON < 0.) ZOBSLON=ZOBSLON+360.

               II=NINT((ZOBSLAT-ZLAT0)/ZLATI0+1._JPRB)
               JJ=NINT((ZOBSLON*RADIANS*FLOAT(NLOENG(II)))/(RPI*2._JPRB)+1._JPRB)
               IGP=SUM(NLOENG(1:II-1))+JJ

               !* if the location is processed on this processors then continue
               !* (this is done by checking the lookup table generated under step
               !* 1.0 for a valid entry)

               IF (INDGP(IGP,1) /= NMDI) THEN

                  !* get the nearest model time step

                  IOBS_DAT   = ROBHDR(IOBS,MDBDAT)
                  IOBS_ETM   = ROBHDR(IOBS,MDBETM)

                  IOBST_YY = INT(IOBS_DAT / 10000)
                  IOBST_MM = INT(IOBS_DAT / 100)-IOBST_YY*100
                  IOBST_DD = INT(IOBS_DAT)-IOBST_YY*10000-IOBST_MM*100

                  IMODT_YY = INT(NINDAT / 10000)
                  IMODT_MM = INT(NINDAT / 100)-IMODT_YY*100
                  IMODT_DD = INT(NINDAT)-IMODT_YY*10000-IMODT_MM*100

                  IRET=0
                  CALL DAYDIFF(IOBST_YY,IOBST_MM,IOBST_DD,IMODT_YY,IMODT_MM,IMODT_DD,IDOFF,IRET)

                  IOBST_HH = INT(IOBS_ETM / 10000)
                  IOBST_MN = INT((IOBS_ETM-IOBST_HH*10000)/100)

                  IOBST = IOBST_HH*60*60 + IOBST_MN*60 + IDOFF*24*60*60
                  IMODTIME = NSSSSS
                  IDIFF = IOBST - NSSSSS
                  JTST=IDIFF/TSTEP

                  IF ((JTST >= NSTART).AND.( JTST <= NSTOP)) THEN

                     IASCAT=IASCAT+1

                     ! WRITE(NULOUT,*) 'ASCAT JTST  : ', IOBS_DAT, IOBS_ETM, IDIFF, JTST
                     ! WRITE(NULOUT,*) 'ASCAT OBS: ', IASCAT, IOBS_DAT, IOBS_ETM, &
                     !     & ZOBSLAT, ZOBSLON, ROBODY(JALL,MDBVAR)

                     !* write data into temporary array.


                     YLASCATOBSTMP(IOBS)%ZOBS    = ROBODY(JALL,MDBVAR)
                     YLASCATOBSTMP(IOBS)%IBODY   = JALL
                     YLASCATOBSTMP(IOBS)%IHDR    = IOBS
                     YLASCATOBSTMP(IOBS)%ISTEP   = JTST
                     YLASCATOBSTMP(IOBS)%IPROMA  = INDGP(IGP,1)
                     YLASCATOBSTMP(IOBS)%IGPBLCK = INDGP(IGP,2)
                     YLASCATOBSTMP(IOBS)%IGP     = IGP

                     ! WRITE(NULOUT,*) 'ASCAT STR   : ', YLASCATOBSTMP(IOBS)

                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDIF


!*     3.0 Calculate First Guess Departures & Discard all bad observations


IF (IASCAT>0) THEN

   ALLOCATE(YLASCATOBSTMP2(IASCAT))

   IND=0
   DO IOBS = 1, ILEN
      IF (YLASCATOBSTMP(IOBS)%ZOBS /= RMDI) THEN

         IND=IND+1

         YLASCATOBSTMP2(IND)=YLASCATOBSTMP(IOBS)
         ZFGDEP = YLASCATOBSTMP2(IND)%ZOBS-&
          & VOL_SM_FG(YLASCATOBSTMP2(IND)%IPROMA,1,YLASCATOBSTMP2(IND)%IGPBLCK,YLASCATOBSTMP2(IND)%ISTEP,1)

         ROBODY(YLASCATOBSTMP2(IND)%IBODY,MDBOMF)=ZFGDEP

       ! WRITE(NULOUT,*) 'ASCAT OBS    :', YLASCATOBSTMP2(IND)
       ! WRITE(NULOUT,*) 'ASCAT FG DEP : ', &
       !   & VOL_SM_FG(YLASCATOBSTMP2(IND)%IPROMA,1,YLASCATOBSTMP2(IND)%IGPBLCK,YLASCATOBSTMP2(IND)%ISTEP,1), &
       !   & YLASCATOBSTMP2(IND)%ZOBS, ZFGDEP

      ENDIF
   ENDDO

ENDIF


!*     4.  Average, if more than one obs per grid box and timestep is found

NOBS_ASCAT_TOT = 0

IF (IASCAT>1) THEN


   !* 4.1 Sort Observations
   !*       Sort key is composed of global index and tsetp

   ALLOCATE(IKEY(IASCAT,2))

   DO IOBS = 1, IASCAT
      IKEY(IOBS,1)=IOBS
      IKEY(IOBS,2)=YLASCATOBSTMP2(IOBS)%IGP*1000+YLASCATOBSTMP2(IOBS)%ISTEP
   ENDDO

   CALL KEYSORT(IRET,IKEY,IASCAT,KEY=2,TRANSPOSED=.FALSE.)


   !* 4.2  Find multiple observations

   !         Loop through observations and look if consecutive entries have the same key
   !         i.e. they belong to the same gridpoint and timestep
   !         if they do then add the observations up, increase the obs counter and set
   !         the respective obs to missing so we can dump it later
   !         if they don't then take the average and start with the next grid point


   ISUNIQ = 1
   IG=IKEY(1,1)
   IGUNIQ=IG
   IGP_N = 1

   DO IS = 2, IASCAT

      IG=IKEY(IS,1)

      IF (IKEY(ISUNIQ,2) /= IKEY(IS,2)) THEN

        YLASCATOBSTMP2(IGUNIQ)%ZOBS = YLASCATOBSTMP2(IGUNIQ)%ZOBS /  IGP_N
!        WRITE(NULOUT,*) 'ASCAT GP UNIQ : ', YLASCATOBSTMP2(IGUNIQ)

        IGP_N = 1
        IGUNIQ=IG
        ISUNIQ=IS
        NOBS_ASCAT_TOT = NOBS_ASCAT_TOT+1


      ELSE

        YLASCATOBSTMP2(IGUNIQ)%ZOBS = YLASCATOBSTMP2(IGUNIQ)%ZOBS + YLASCATOBSTMP2(IG)%ZOBS
        IGP_N = IGP_N + 1
!       WRITE(NULOUT,*) 'ASCAT GP  : ',IGP_N, YLASCATOBSTMP2(IG)
        YLASCATOBSTMP2(IG)%ZOBS = RMDI

     ENDIF

  ENDDO

  !* The last entry has to be dealt with seperatley

  YLASCATOBSTMP2(IGUNIQ)%ZOBS = YLASCATOBSTMP2(IGUNIQ)%ZOBS /  IGP_N
  NOBS_ASCAT_TOT = NOBS_ASCAT_TOT+1
!  WRITE(NULOUT,*) 'ASCAT GP UNIQ : ', YLASCATOBSTMP2(IGUNIQ)


  !* 4.3  Dump multiple observations
  !       & store number of observations per gridbox (can be more than one, i.e. at different tsteps)

  ALLOCATE(YLASCATOBS(NOBS_ASCAT_TOT))

  IND=0
  DO IOBS = 1, IASCAT
      IF (YLASCATOBSTMP2(IOBS)%ZOBS /= RMDI) THEN

         IND=IND+1
         YLASCATOBS(IND)=YLASCATOBSTMP2(IOBS)

         NOBS_ASCAT(YLASCATOBS(IND)%IPROMA,YLASCATOBS(IND)%IGPBLCK) =&
          & NOBS_ASCAT(YLASCATOBS(IND)%IPROMA,YLASCATOBS(IND)%IGPBLCK)+1

         ! WRITE(NULOUT,*) 'ASCAT GP OBS   :', YLASCATOBS(IND)

      ENDIF
   ENDDO

ENDIF

!*     4.  In the rare case that only one obs is found just copy the array

IF (IASCAT==1) THEN

   IF (YLASCATOBSTMP2(1)%ZOBS /= RMDI) THEN

      NOBS_ASCAT_TOT = 1
      ALLOCATE(YLASCATOBS(NOBS_ASCAT_TOT))

      YLASCATOBS(1)=YLASCATOBSTMP2(1)

      NOBS_ASCAT(YLASCATOBS(1)%IPROMA,YLASCATOBS(1)%IGPBLCK) =&
          & NOBS_ASCAT(YLASCATOBS(1)%IPROMA,YLASCATOBS(1)%IGPBLCK)+1

      ! WRITE(NULOUT,*) 'ASCAT GP OBS :', YLASCATOBS(IND)

   ENDIF

ENDIF

WRITE(NULOUT,*) 'ASCAT OBS PROCESSED (PROC/TOTOBS/VALID/TOTOBS after GPAVG): ', MYPROC, ILEN, IASCAT, NOBS_ASCAT_TOT


!*     5.0  Clean up


!*    5.1 PUT THE DB AWAY

IRET=0
CALL PUTDB('ASCATSM',1,IRET,INFO,0)
LSCREEN=.FALSE.


!*    5.2 Deallocate arrays

IF (ALLOCATED(YLASCATOBSTMP)) DEALLOCATE(YLASCATOBSTMP)
IF (ALLOCATED(YLASCATOBSTMP2)) DEALLOCATE(YLASCATOBSTMP2)
IF (ALLOCATED(IKEY)) DEALLOCATE(IKEY)
IF (ALLOCATED(INDGP)) DEALLOCATE(INDGP)
IF (ALLOCATED(IKINFO)) DEALLOCATE(IKINFO)


!*   5.3 Flush output

CALL FLUSH(NULOUT)


!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SEKF_PREP_ASCAT',1,ZHOOK_HANDLE)

#include "openmp_obs_undef.h"

END SUBROUTINE SEKF_PREP_ASCAT
