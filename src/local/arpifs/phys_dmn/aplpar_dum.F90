!OPTIONS XOPT(NOEVAL)
SUBROUTINE APLPAR_DUM(&
   & YDGEOMETRY,         & 
   & YDMF_PHYS,          & 
   & YDCPG_DYN,          & 
   & YDMF_PHYS_SURF,     & 
   & YDVARS,             & 
   & YDSURF,             & 
   & YDXFU,              & 
   & YDCFU,              & 
   & YDMODEL,            & 
   & KIDIA,              & 
   & KFDIA,              & 
   & KLON,               & 
   & KTDIA,              & 
   & KLEV,               & 
   & KSTGLO,             & 
   & KVCLIS,             & 
   & KVCLIV,             & 
   & KSTEP,              & 
   & KSGST,              & 
   & KCSS,               & 
   & KBL,                & 
   & KGPCOMP,            & 
   & KNFRRC,             & 
#ifdef UNDEF
   & PDT,                & 
   & PINDX,              & 
   & PINDY,              & 
   & LDXFUMSE,           & 
   & PAPHI,              & 
   & PAPRS,              & 
   & PAPHIF,             & 
   & PAPRSF,             & 
   & PALPH,              & 
   & PDELP,              & 
   & PLNPR,              & 
   & PRDELP,             & 
   & PRCORI,             & 
   & PSV,                & 
   & PU,                 & 
   & PV,                 & 
   & PT,                 & 
   & PQ,                 & 
   & PQI,                & 
   & PQL,                & 
   & PQSRES,             & 
   & PQRRES,             & 
   & PG,                 & 
   & PTKE,               & 
   & PEFB1,              & 
   & PEFB2,              & 
   & PEFB3,              & 
   & PCVV,               & 
   & PO3,                & 
   & PCHEM,              & 
   & PNOGWDU,            & 
   & PNOGWDV,            & 
   & PGFL,               & 
   & PVORT0,             & 
   & PCP,                & 
   & PCVGQ,              & 
   & PR,                 & 
   & PKOZO,              & 
   & PFPLCH,             & 
   & PFPLSH,             & 
   & PEVEL0,             & 
   & PGPAR,              & 
   & PSNS,               & 
   & PALBNS,             & 
   & PRHONS,             & 
   & PTP,                & 
   & PTS,                & 
   & PWL,                & 
   & PWP,                & 
   & PWPI,               & 
   & PWS,                & 
#endif
   & PWSI,               & 
   & PVERVEL,            & 
   & PEMTD,              & 
   & PEMTU,              & 
   & PTRSO,              & 
   & PVRMOON,            & 
   & PMU0,               & 
   & PMU0LU,             & 
   & PMU0M,              & 
   & PMU0N,              & 
   & PGELAM,             & 
   & PGEMU,              & 
   & PGM,                & 
   & POMPAC,             & 
   & PAC,                & 
   & PCOR,               & 
   & PMMU0,              & 
   & PDHSF,              & 
   & PRAB3C,             & 
   & PRAB3N,             & 
   & PRAB4C,             & 
   & PRAB4N,             & 
   & PRAB6C,             & 
   & PRAB6N,             & 
   & PRAT1C,             & 
   & PRAT1N,             & 
   & PRAT2C,             & 
   & PRAT2N,             & 
   & PRAT3C,             & 
   & PRAT3N,             & 
   & PRAT4C,             & 
   & PRAT4N,             & 
   & PRAT5C,             & 
   & PRAT5N,             & 
   & POROG,              & 
   & PWW,                & 
   & PDIV,               & 
   & PUL,                & 
   & PVL,                & 
   & PWWL,               & 
   & PWWM,               & 
   & PGDEOSI,            & 
   & PGUEOSI,            & 
   & PGMU0,              & 
   & PGMU0_MIN,          & 
   & PGMU0_MAX,          & 
   & PGDEOTI,            & 
   & PGDEOTI2,           & 
   & PGUEOTI,            & 
   & PGUEOTI2,           & 
   & PGEOLT,             & 
   & PGEOXT,             & 
   & PGRPROX,            & 
   & PGMIXP,             & 
   & PGFLUXC,            & 
   & PGRSURF,            & 
   & PDIFSV,             & 
   & PFRMQ,              & 
   & PCPS,               & 
   & PLHS,               & 
   & PRS,                & 
   & PLH,                & 
   & PLSCPE,             & 
   & PNEB,               & 
   & PQICE,              & 
   & PQLI,               & 
   & PQRCONV1,           & 
   & PQSCONV1,           & 
   & PQSAT,              & 
   & PQW,                & 
   & PRH,                & 
   & PTW,                & 
   & PCD,                & 
   & PCDN,               & 
   & PCH,                & 
   & PC1,                & 
   & PC2,                & 
   & PEMIS,              & 
   & PFEVI,              & 
   & PFTKE,              & 
   & PFTKEI,             & 
   & PFEFB1,             & 
   & PFEFB2,             & 
   & PFEFB3,             & 
   & PNEIJ,              & 
   & PVEG,               & 
   & PQS,                & 
   & PQSATS,             & 
   & PCLCT,              & 
   & KCLPH,              & 
   & PDPRECIPS,          & 
   & PDPRECIPS2,         & 
   & PTENDPTKE,          & 
   & PKUROV_H,           & 
   & PKTROV_H,           & 
   & PTENDEXT_DEP,       & 
   & PTRAJ_PHYS,         & 
   & PTDISS,             & 
   & YDDDH,              & 
   & PFTCNS,             & 
   & PGP2DSPP)           

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE MF_PHYS_TYPE_MOD   , ONLY : MF_PHYS_TYPE
USE CPG_TYPE_MOD   , ONLY : CPG_DYN_TYPE
USE MFPHYS_SURFACE_TYPE_MOD,ONLY : MFPHYS_SURF_TYPE
USE FIELD_VARIABLES_MOD, ONLY : FIELD_VARIABLES
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMXFU             , ONLY : TXFU
USE TYPE_MODEL         , ONLY : MODEL
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK    ,DR_HOOK
USE YOMMP0             , ONLY : NPRINTLEV
USE YOMVERT            , ONLY : VP00
USE YOMCST             , ONLY : RG       ,RSIGMA   ,RV       ,RD       ,&
                              & RCPV     ,RETV     ,RCW      ,RCS      ,RLVTT    ,&
                              & RLSTT    ,RTT      ,RALPW    ,RBETW    ,RGAMW    ,&
                              & RALPS    ,RBETS    ,RGAMS    ,RALPD    ,RBETD    ,&
                              & RGAMD    ,RCPD     ,RATM     ,RKAPPA   ,RLMLT
USE YOMCT0             , ONLY : LCALLSFX ,LSFORCS, LELAM
USE YOMDYNA            , ONLY : L3DTURB
USE YOMRIP0            , ONLY : NINDAT
USE DDH_MIX            , ONLY : ADD_FIELD_3D ,ADD_FIELD_2D ,NTOTSVAR , & 
                              & NTOTSURF ,NTOTSVFS, NEW_ADD_FIELD_3D, NEW_ADD_FIELD_2D, &
                              & TYP_DDH
USE YOMLUN             , ONLY : NULOUT
USE YOMLSFORC          , ONLY : LMUSCLFA,NMUSCLFA
USE YOMTRAJ            , ONLY : TRAJ_PHYS_TYPE
USE YOMCFU             , ONLY : TCFU 
USE SPP_MOD  , ONLY : YSPP, YSPP_CONFIG

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MF_PHYS_TYPE),INTENT(INOUT) :: YDMF_PHYS
TYPE(CPG_DYN_TYPE),INTENT(INOUT) :: YDCPG_DYN
TYPE(MFPHYS_SURF_TYPE),INTENT(INOUT) :: YDMF_PHYS_SURF
TYPE(FIELD_VARIABLES),INTENT(INOUT) :: YDVARS
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
TYPE(TXFU)        ,INTENT(INOUT) :: YDXFU
TYPE(TCFU)        ,INTENT(INOUT) :: YDCFU
TYPE(MODEL)       ,INTENT(INOUT),TARGET :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
INTEGER(KIND=JPIM),INTENT(IN)    :: KBL
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPCOMP
INTEGER(KIND=JPIM),INTENT(IN)    :: KNFRRC
INTEGER(KIND=JPIM),INTENT(IN)    :: KVCLIS
INTEGER(KIND=JPIM),INTENT(IN)    :: KSGST
INTEGER(KIND=JPIM),INTENT(IN)    :: KCSS
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KTDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTGLO
INTEGER(KIND=JPIM),INTENT(IN)    :: KVCLIV
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTEP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT
LOGICAL           ,INTENT(IN)    :: LDXFUMSE
REAL(KIND=JPRB)   ,INTENT(IN)    :: PINDX(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PINDY(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHI(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHIF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPH(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDELP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLNPR(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRDELP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRCORI(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSV(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQSRES(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQRRES(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PG(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTKE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEFB1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEFB2(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEFB3(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCVV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQI(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQL(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PO3(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCHEM(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NCHEM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNOGWDU(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNOGWDV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVORT0(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCP(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCVGQ(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PR(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEVEL0(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PKOZO(KLON,KLEV,KVCLIS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLCH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFPLSH(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSNS(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGPAR(KLON,YDMODEL%YRML_PHY_MF%YRPARAR%NGPAR+1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBNS(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHONS(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTP(KLON,KCSS)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTS(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWL(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWP(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWPI(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWS(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWSI(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVERVEL(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEMTD(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEMTU(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTRSO(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVRMOON(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0LU(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0M(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0N(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGELAM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEMU(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: POROG(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWW(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIV(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUL(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVL(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWWL(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWWM(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGDEOSI(KLON,0:KLEV,2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGUEOSI(KLON,0:KLEV,2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMU0(KLON,0:YDMODEL%YRML_PHY_MF%YRPHY%NSORAYFR-1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMU0_MIN(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMU0_MAX(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGDEOTI(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGDEOTI2(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGUEOTI(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGUEOTI2(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGEOLT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGEOXT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGRPROX(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGMIXP(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLUXC(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGRSURF(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: POMPAC(KLON,0:KLEV,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PAC(KLEV+1,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOR(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMMU0(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDHSF(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAB3C(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAB3N(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAB4C(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAB4N(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAB6C(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAB6N(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAT1C(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAT1N(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAT2C(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAT2N(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAT3C(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAT3N(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAT4C(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAT4N(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAT5C(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRAT5N(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDIFSV(KLON,0:KLEV,1:YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRMQ(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFTKE(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFTKEI(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFEFB1(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFEFB2(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFEFB3(KLON,0:KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCPS(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLHS(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRS(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLH(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PLSCPE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNEB(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQICE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQLI(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQRCONV1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQSCONV1(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQSAT(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQW(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRH(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTW(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCD(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCDN(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PCH(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PC1(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PC2(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEMIS(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFEVI(KLON,KSGST+1)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PNEIJ(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVEG(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQS(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQSATS(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCLCT(KLON)
INTEGER(KIND=JPIM),INTENT(OUT)   :: KCLPH(KLON)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PFTCNS(KLON,0:KLEV,6)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTENDPTKE(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKUROV_H(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PKTROV_H(KLON,KLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENDEXT_DEP(KLON,KLEV,YDMODEL%YRML_GCONF%YGFL%NGFL_EXT)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDPRECIPS(KLON,YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS%NDTPREC)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDPRECIPS2(KLON,YDMODEL%YRML_PHY_MF%YRPHY%YRDPRECIPS%NDTPREC2)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTDISS(KLON,KLEV)
TYPE (TRAJ_PHYS_TYPE), INTENT(INOUT) :: PTRAJ_PHYS
TYPE(TYP_DDH), INTENT(INOUT)     :: YDDDH
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGP2DSPP(KLON,YSPP%N2D)

END SUBROUTINE APLPAR_DUM
