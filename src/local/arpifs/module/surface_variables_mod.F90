
MODULE SURFACE_VARIABLES_MOD
  ! The SURFACE_VARIABLES type provides namespaced access to individual
  ! groups of surface VARIABLE objects via group-specific container
  ! types.
  !
  ! The generated group types contain the set of VARIABLE objects
  ! configured for each surface group. The VARIABLES in turn that hold
  ! metadata for surface fields and provide access to the underlying
  ! storage FIELD objects and respective data pointers for
  ! thread-parallel regions.

USE PARKIND1, ONLY: JPIM, JPRB
USE VARIABLE_MODULE, ONLY: VARIABLE_2D, VARIABLE_3D
USE FIELD_MODULE, ONLY: FIELD_3D, FIELD_4D

IMPLICIT NONE

! Prognostic variable group types
TYPE SURFACE_VARIABLE_GROUP_SOILB
  ! Prognostic surface variable group
  TYPE(VARIABLE_3D) :: VT   ! temperature
  TYPE(VARIABLE_3D) :: VQ   ! liquid water content
  TYPE(VARIABLE_3D) :: VTL   ! ice water content (for MF)

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_4D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_SOILB

TYPE SURFACE_VARIABLE_GROUP_SNOWG
  ! Prognostic surface variable group
  TYPE(VARIABLE_3D) :: VF   ! content of surface snow
  TYPE(VARIABLE_3D) :: VA   ! snow albedo
  TYPE(VARIABLE_3D) :: VR   ! snow density
  TYPE(VARIABLE_3D) :: VT   ! total albedo (diagnostic for MF for LVGSN)
  TYPE(VARIABLE_3D) :: VW   ! Liquid water content

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_4D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_SNOWG

TYPE SURFACE_VARIABLE_GROUP_LAKEB
  ! Prognostic surface variable group
  TYPE(VARIABLE_2D) :: VLICT   ! lake ice temperature
  TYPE(VARIABLE_2D) :: VLMLT   ! lake mixed-layer temperature
  TYPE(VARIABLE_2D) :: VLTLT   ! lake total layer temperature
  TYPE(VARIABLE_2D) :: VLBLT   ! lake bottom layer temperature
  TYPE(VARIABLE_2D) :: VLSHF   ! lake shape factor
  TYPE(VARIABLE_2D) :: VLICD   ! lake ice depth
  TYPE(VARIABLE_2D) :: VLMLD   ! lake mixed-layer depth

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_LAKEB

TYPE SURFACE_VARIABLE_GROUP_RESVR
  ! Prognostic surface variable group
  TYPE(VARIABLE_2D) :: VT   ! skin temperature (Ts)
  TYPE(VARIABLE_2D) :: VW   ! skin water content (Wskin) at ECMWF superficial reservoir water content (Ws) at MF
  TYPE(VARIABLE_2D) :: VFC   ! skin water content (Wl) at MF
  TYPE(VARIABLE_2D) :: VIC   ! superficial reservoir ice
  TYPE(VARIABLE_2D) :: VFP1   ! interpolated Ts for 2nd part of 927-FULLPOS

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_RESVR

TYPE SURFACE_VARIABLE_GROUP_CLS
  ! Prognostic surface variable group
  TYPE(VARIABLE_2D) :: VTCLS   ! 2m temperature
  TYPE(VARIABLE_2D) :: VHUCLS   ! 2m humidity
  TYPE(VARIABLE_2D) :: VUCLS   ! 10m U-wind
  TYPE(VARIABLE_2D) :: VVCLS   ! 10m V-wind
  TYPE(VARIABLE_2D) :: VNUCLS   ! 10m neutral U-wind
  TYPE(VARIABLE_2D) :: VNVCLS   ! 10m neutral V-wind

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_CLS

TYPE SURFACE_VARIABLE_GROUP_OML
  ! Prognostic surface variable group
  TYPE(VARIABLE_3D) :: VTO   ! temperature
  TYPE(VARIABLE_3D) :: VSO   ! salinity
  TYPE(VARIABLE_3D) :: VUO   ! U velocity
  TYPE(VARIABLE_3D) :: VVO   ! V velocity

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_4D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_OML

TYPE SURFACE_VARIABLE_GROUP_EXTRP
  ! Prognostic surface variable group

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_4D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_EXTRP

TYPE SURFACE_VARIABLE_GROUP_XTRP2
  ! Prognostic surface variable group

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_XTRP2

TYPE SURFACE_VARIABLE_GROUP_CANRI
  ! Prognostic surface variable group

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_CANRI


! Diagnostic variable group types
TYPE SURFACE_VARIABLE_GROUP_VARSF
  ! Diagnostic surface variable group
  TYPE(VARIABLE_2D) :: VZ0F   ! gravity * surface roughness length
  TYPE(VARIABLE_2D) :: VALBF   ! surface shortwave albedo
  TYPE(VARIABLE_2D) :: VEMISF   ! surface longwave emissivity
  TYPE(VARIABLE_2D) :: VGETRL   ! standard deviation of orography
  TYPE(VARIABLE_2D) :: VLSM   ! land-sea mask
  TYPE(VARIABLE_2D) :: VVEG   ! vegetation cover
  TYPE(VARIABLE_2D) :: VVRLAN   ! anisotropy of the sub-grid scale orography
  TYPE(VARIABLE_2D) :: VVRLDI   ! angle of the direction of orography with the x axis
  TYPE(VARIABLE_2D) :: VSIG   ! characteristic orographic slope
  TYPE(VARIABLE_2D) :: VALBSF   ! soil shortwave albedo
  TYPE(VARIABLE_2D) :: VLAN   ! fraction of land
  TYPE(VARIABLE_2D) :: VSST   ! (open) sea surface temperature
  TYPE(VARIABLE_2D) :: VSSS   ! sea surface salinity
  TYPE(VARIABLE_2D) :: VLZ0H   ! logarithm of roughness length for heat
  TYPE(VARIABLE_2D) :: VCVL   ! low vegetation cover
  TYPE(VARIABLE_2D) :: VCVH   ! high vegetation cover
  TYPE(VARIABLE_2D) :: VTVL   ! low vegetation type
  TYPE(VARIABLE_2D) :: VTVH   ! high vegetation type
  TYPE(VARIABLE_2D) :: VLAIL   ! low vegetation LAI
  TYPE(VARIABLE_2D) :: VLAIH   ! high vegetation LAI
  TYPE(VARIABLE_2D) :: VSOTY   ! soil type
  TYPE(VARIABLE_2D) :: VCLK   ! lake cover
  TYPE(VARIABLE_2D) :: VDL   ! lake depth
  TYPE(VARIABLE_2D) :: VCI   ! sea ice fraction
  TYPE(VARIABLE_2D) :: VUCUR   ! U-component of the ocean current
  TYPE(VARIABLE_2D) :: VVCUR   ! V-component of the ocean current
  TYPE(VARIABLE_2D) :: VZ0RLF   ! gravity * vegetation roughness length
  TYPE(VARIABLE_2D) :: VCO2O   ! oceanic CO2 flux
  TYPE(VARIABLE_2D) :: VCO2B   ! biosphere CO2 flux
  TYPE(VARIABLE_2D) :: VCO2A   ! anthropogenic CO2 flux
  TYPE(VARIABLE_2D) :: VCO2F   ! CO2 fire emissions
  TYPE(VARIABLE_2D) :: VCGPP   ! GPP bias correction factor
  TYPE(VARIABLE_2D) :: VCREC   ! REC bias correction factor
  TYPE(VARIABLE_2D) :: VCH4AG   ! CH4 surface fluxes - aggregated field
  TYPE(VARIABLE_2D) :: VCH4F   ! CH4 fire emissions
  TYPE(VARIABLE_2D) :: VSDFOR   ! SD filtered orography
  TYPE(VARIABLE_2D) :: VALUVP   ! MODIS-derived parallel albedo for shortwave radiation
  TYPE(VARIABLE_2D) :: VALUVD   ! MODIS-derived diffuse albedo for shortwave radiation
  TYPE(VARIABLE_2D) :: VALNIP   ! MODIS-derived parallel albedo for longwave radiation
  TYPE(VARIABLE_2D) :: VALNID   ! MODIS-derived diffuse albedo for longwave radiation
  TYPE(VARIABLE_2D) :: VFP1   ! surface orography in the 2nd part of FULLPOS-927
  TYPE(VARIABLE_2D) :: VBCBF   ! black carbon biogenic
  TYPE(VARIABLE_2D) :: VBCFF   ! black carbon fossil fuel
  TYPE(VARIABLE_2D) :: VBCGF   ! black carbon GFED
  TYPE(VARIABLE_2D) :: VOMBF   ! organic matter biogenic
  TYPE(VARIABLE_2D) :: VOMFF   ! organic matter fossil fuel
  TYPE(VARIABLE_2D) :: VOMGF   ! organic matter GFED
  TYPE(VARIABLE_2D) :: VINJF   ! height of maximum injection for biomass burning emissions
  TYPE(VARIABLE_2D) :: VSO2L   ! sulphate low-level
  TYPE(VARIABLE_2D) :: VSO2H   ! sulphate higher-level
  TYPE(VARIABLE_2D) :: VSO2DD   ! sulphate dry dep velocity
  TYPE(VARIABLE_2D) :: VSOGF   ! sulphate GFED
  TYPE(VARIABLE_2D) :: VSOA   ! secondary organic
  TYPE(VARIABLE_2D) :: VVOLC   ! volcanic continuous
  TYPE(VARIABLE_2D) :: VVOLE   ! volcanic explosive
  TYPE(VARIABLE_2D) :: VDMSO   ! oceanic DMS
  TYPE(VARIABLE_2D) :: VSOACO   ! SOA from CO
  TYPE(VARIABLE_2D) :: VURBF   ! Urban fraction
  TYPE(VARIABLE_2D) :: VVOLCALTI   ! Altitude of volcanoes
  TYPE(VARIABLE_2D) :: VFCA1   ! Fraction of calcite over dust 1st bin
  TYPE(VARIABLE_2D) :: VFCA2   ! Fraction of calcite over dust 2nd bin
  TYPE(VARIABLE_2D) :: VAERDEP   ! dust emission potential
  TYPE(VARIABLE_2D) :: VAERLTS   ! dust lifting threshold speed
  TYPE(VARIABLE_2D) :: VAERSCC   ! dust soil clay content
  TYPE(VARIABLE_2D) :: VDSF   ! dust source function
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:) :: VCHEMFLX   ! chemistry emissions input
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:) :: VCHEMFLXO   ! total chemistry flux (emissions + deposition)
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:) :: VCHEMDV   ! chemistry deposition velocity
  TYPE(VARIABLE_2D) :: VNUDM   ! nudging mask

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VARSF

TYPE SURFACE_VARIABLE_GROUP_VCLIH
  ! Diagnostic surface variable group
  TYPE(VARIABLE_2D) :: VTCCH   ! total convective cloudiness
  TYPE(VARIABLE_2D) :: VSCCH   ! convective cloud summit
  TYPE(VARIABLE_2D) :: VBCCH   ! convective cloud base
  TYPE(VARIABLE_2D) :: VPBLH   ! PBL height
  TYPE(VARIABLE_2D) :: VSPSH   ! variable for prognostic convection scheme (ALARO)
  TYPE(VARIABLE_2D) :: VQSH   ! surface moisture historic variable (used by TOUCANS)
  TYPE(VARIABLE_2D) :: VPCL   ! 
  TYPE(VARIABLE_2D) :: VPSL   ! 
  TYPE(VARIABLE_2D) :: VPCN   ! 
  TYPE(VARIABLE_2D) :: VPSN   ! 
  TYPE(VARIABLE_2D) :: VEVA   ! 

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VCLIH

TYPE SURFACE_VARIABLE_GROUP_VCLIK
  ! Diagnostic surface variable group
  TYPE(VARIABLE_2D) :: VUDGRO   ! ud top position (accsu)

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VCLIK

TYPE SURFACE_VARIABLE_GROUP_VCLIP
  ! Diagnostic surface variable group
  TYPE(VARIABLE_2D) :: VTPC   ! climatological deep layer temperature
  TYPE(VARIABLE_2D) :: VWPC   ! climatological deep layer moisture

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VCLIP

TYPE SURFACE_VARIABLE_GROUP_VCLIV
  ! Diagnostic surface variable group
  TYPE(VARIABLE_2D) :: VARG   ! silt percentage within soil
  TYPE(VARIABLE_2D) :: VSAB   ! percentage of sand within the soil
  TYPE(VARIABLE_2D) :: VD2   ! soil depth
  TYPE(VARIABLE_2D) :: VIVEG   ! type of vegetation
  TYPE(VARIABLE_2D) :: VRSMIN   ! stomatal minimum resistance
  TYPE(VARIABLE_2D) :: VLAI   ! leaf area index
  TYPE(VARIABLE_2D) :: VHV   ! resistance to evapotranspiration
  TYPE(VARIABLE_2D) :: VZ0H   ! gravity * roughness length for heat
  TYPE(VARIABLE_2D) :: VALS   ! albedo of bare ground
  TYPE(VARIABLE_2D) :: VALV   ! albedo of vegetation

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VCLIV

TYPE SURFACE_VARIABLE_GROUP_VCLIA
  ! Diagnostic surface variable group
  TYPE(VARIABLE_2D) :: VSEA   ! aerosol sea
  TYPE(VARIABLE_2D) :: VLAN   ! aerosol land
  TYPE(VARIABLE_2D) :: VSOO   ! aerosol soot
  TYPE(VARIABLE_2D) :: VDES   ! aerosol desert
  TYPE(VARIABLE_2D) :: VSUL   ! aerosol sulfate
  TYPE(VARIABLE_2D) :: VVOL   ! aerosol volcano

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VCLIA

TYPE SURFACE_VARIABLE_GROUP_VCLIN
  ! Diagnostic surface variable group
  TYPE(VARIABLE_2D) :: VTOP   ! index of convective cloud top
  TYPE(VARIABLE_2D) :: VBAS   ! index of convective cloud base
  TYPE(VARIABLE_2D) :: VACPR   ! averaged convective precipitaion rate
  TYPE(VARIABLE_2D) :: VACCPR   ! accumulated total precipitaion for assimilation
  TYPE(VARIABLE_2D) :: VACCPR5   ! accumulated total precipitaion for assimilation (trajectory)

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VCLIN

TYPE SURFACE_VARIABLE_GROUP_VDIAGO2
  ! Diagnostic surface variable group
  TYPE(VARIABLE_2D) :: VOCDEP   ! bottom layer depth
  TYPE(VARIABLE_2D) :: VUSTRC   ! taux clim.
  TYPE(VARIABLE_2D) :: VVSTRC   ! tauy clim.

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VDIAGO2

TYPE SURFACE_VARIABLE_GROUP_VDIAGO3
  ! Diagnostic surface variable group
  TYPE(VARIABLE_3D) :: VDIFM   ! viscosity
  TYPE(VARIABLE_3D) :: VDIFT   ! diff. coef. of temp
  TYPE(VARIABLE_3D) :: VDIFS   ! diff. coef. of salinity
  TYPE(VARIABLE_3D) :: VADVT   ! correction term for temp.
  TYPE(VARIABLE_3D) :: VADVS   ! correction term for sal.
  TYPE(VARIABLE_3D) :: VTRI0   ! coef. for solving matrix.
  TYPE(VARIABLE_3D) :: VTRI1   ! coef. for solving matrix.
  TYPE(VARIABLE_3D) :: VSWDK   ! radiation term
  TYPE(VARIABLE_3D) :: VZO   ! depth of layer
  TYPE(VARIABLE_3D) :: VHO   ! depth of interface layer
  TYPE(VARIABLE_3D) :: VDO   ! layer thickness
  TYPE(VARIABLE_3D) :: VHO_INV   ! 1 / YHO
  TYPE(VARIABLE_3D) :: VUOC   ! U velocity clim.
  TYPE(VARIABLE_3D) :: VVOC   ! V velocity clim.

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_4D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VDIAGO3

TYPE SURFACE_VARIABLE_GROUP_VDIAG
  ! Diagnostic surface variable group
  TYPE(VARIABLE_2D) :: VLSP   ! Large scale precipitation
  TYPE(VARIABLE_2D) :: VCP   ! Convective precipitation
  TYPE(VARIABLE_2D) :: VSF   ! Snowfall
  TYPE(VARIABLE_2D) :: VFZRA   ! Freezing rain
  TYPE(VARIABLE_2D) :: VBLD   ! Boundary layer dissipation
  TYPE(VARIABLE_2D) :: VSSHF   ! Surface sensible heat flux
  TYPE(VARIABLE_2D) :: VSLHF   ! Surface latent heat flux
  TYPE(VARIABLE_2D) :: VNEE   ! Surface net ecosystem exchange of CO2
  TYPE(VARIABLE_2D) :: VGPP   ! Surface gross primary production of CO2
  TYPE(VARIABLE_2D) :: VREC   ! Surface ecosystem respiration of CO2
  TYPE(VARIABLE_2D) :: VMSL   ! Mean sea level pressure
  TYPE(VARIABLE_2D) :: VSP   ! Surface pressure
  TYPE(VARIABLE_2D) :: VTCC   ! Total cloud cover
  TYPE(VARIABLE_2D) :: V10U   ! U-wind at 10 m
  TYPE(VARIABLE_2D) :: V10V   ! V-wind at 10 m
  TYPE(VARIABLE_2D) :: V2T   ! Temperature at 2 m
  TYPE(VARIABLE_2D) :: V2D   ! Dewpoint temperature at 2 m
  TYPE(VARIABLE_2D) :: V2Q   ! Specific humidity at 2 m
  TYPE(VARIABLE_2D) :: VSSR   ! Surface solar radiation
  TYPE(VARIABLE_2D) :: VSTR   ! Surface thermal radiation
  TYPE(VARIABLE_2D) :: VTSR   ! Top solar radiation
  TYPE(VARIABLE_2D) :: VTTR   ! Top thermal radiation
  TYPE(VARIABLE_2D) :: VEWSS   ! Instantaneous surface U-wind stress
  TYPE(VARIABLE_2D) :: VNSSS   ! Instantaneous surface V-wind stress
  TYPE(VARIABLE_2D) :: VE   ! Water evaporation
  TYPE(VARIABLE_2D) :: VPEV   ! Potential evaporation
  TYPE(VARIABLE_2D) :: VCCC   ! Convective cloud cover
  TYPE(VARIABLE_2D) :: VLCC   ! Low cloud cover
  TYPE(VARIABLE_2D) :: VMCC   ! Medium cloud cover
  TYPE(VARIABLE_2D) :: VHCC   ! High cloud cover
  TYPE(VARIABLE_2D) :: VLGWS   ! Zonal gravity wave stress
  TYPE(VARIABLE_2D) :: VMGWS   ! Meridian gravity wave stress
  TYPE(VARIABLE_2D) :: VGWD   ! Gravity wave dissipation
  TYPE(VARIABLE_2D) :: VMX2T   ! Maximum temperature at 2 m
  TYPE(VARIABLE_2D) :: VMN2T   ! Minimum temperature at 2 m
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:) :: VMX2T6   ! Bins for maximum temperature at 2 m since last 6 hours
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:) :: VMN2T6   ! Bins for minimum temperature at 2 m since last 6 hours
  TYPE(VARIABLE_2D) :: VRO   ! Runoff (total)
  TYPE(VARIABLE_2D) :: VSRO   ! Runoff surface
  TYPE(VARIABLE_2D) :: VSSRO   ! Runoff sub-surface
  TYPE(VARIABLE_2D) :: VALB   ! (surface shortwave) albedo
  TYPE(VARIABLE_2D) :: VIEWSS   ! Instantaneous surface zonal component of stress
  TYPE(VARIABLE_2D) :: VINSSS   ! Instantaneous surface meridian component of stress
  TYPE(VARIABLE_2D) :: VISSHF   ! Instantaneous surface heat flux
  TYPE(VARIABLE_2D) :: VIE   ! Instantaneous surface moisture flux
  TYPE(VARIABLE_2D) :: VINEE   ! Instantaneous net ecosystem exchange of CO2
  TYPE(VARIABLE_2D) :: VIGPP   ! Instantaneous gross primary production of CO2
  TYPE(VARIABLE_2D) :: VIREC   ! Instantaneous ecosystem respiration of CO2
  TYPE(VARIABLE_2D) :: VCSF   ! Convective snow fall
  TYPE(VARIABLE_2D) :: VLSSF   ! Large scale snowfall
  TYPE(VARIABLE_2D) :: VMXTPR   ! Max precip rate since last post-processing
  TYPE(VARIABLE_2D) :: VMNTPR   ! Min precip rate since last post-processing
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:) :: VMXTPR6   ! Max precip rate in last 6 hours
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:) :: VMNTPR6   ! Min precip rate in last 6 hours
  TYPE(VARIABLE_2D) :: VTPR   ! Total precipitation rate
  TYPE(VARIABLE_2D) :: VLSRR   ! Large scale rain rate
  TYPE(VARIABLE_2D) :: VCRR   ! Convective rain rate
  TYPE(VARIABLE_2D) :: VLSSFR   ! Large scale snowfall rate
  TYPE(VARIABLE_2D) :: VCSFR   ! Convective snowfall rate
  TYPE(VARIABLE_2D) :: VPTYPE   ! Precipitation type
  TYPE(VARIABLE_2D) :: VILSPF   ! Large-scale precipitation fraction (inst.)
  TYPE(VARIABLE_2D) :: VZ0F   ! Gravity * surface roughness length
  TYPE(VARIABLE_2D) :: VLZ0H   ! Logarithm of z0 times heat flux
  TYPE(VARIABLE_2D) :: VVIWVE   ! Vertical integral of eastward water vapour flux
  TYPE(VARIABLE_2D) :: VVIWVN   ! Vertical integral of northward water vapour flux
  TYPE(VARIABLE_2D) :: VTCW   ! Total water content in a vertical column
  TYPE(VARIABLE_2D) :: VTCWV   ! Total water vapor content in a vertical column
  TYPE(VARIABLE_2D) :: VTCLW   ! Total liquid water content in a vertical column
  TYPE(VARIABLE_2D) :: VTCIW   ! Total ice water content in a vertical column
  TYPE(VARIABLE_2D) :: VTCRW   ! Total rain water content in a vertical column
  TYPE(VARIABLE_2D) :: VTCSW   ! Total snow water content in a vertical column
  TYPE(VARIABLE_2D) :: VTCSLW   ! Total supercooled liquid water content in a vertical column
  TYPE(VARIABLE_2D) :: VSSRD   ! Downward surface solar radiation
  TYPE(VARIABLE_2D) :: VSTRD   ! Downward surface thermic radiation
  TYPE(VARIABLE_2D) :: VSSRDC   ! Clear-sky downward surface solar radiation
  TYPE(VARIABLE_2D) :: VSTRDC   ! Claer-sky downward surface thermal radiation
  TYPE(VARIABLE_2D) :: VBLH   ! Height of boundary layer
  TYPE(VARIABLE_2D) :: VSUND   ! Sunshine duration
  TYPE(VARIABLE_2D) :: VSPAR   ! Surface downward PARadiation
  TYPE(VARIABLE_2D) :: VSUVB   ! Surface downward UV-B radiation
  TYPE(VARIABLE_2D) :: VSFDIR   ! Surface total sky direct downward SW radiation
  TYPE(VARIABLE_2D) :: VSCDIR   ! Surface clear-sky direct downward SW radiation
  TYPE(VARIABLE_2D) :: VSDSRP   ! Surface total-sky direct beam downward SW radiation
  TYPE(VARIABLE_2D) :: VCAPE   ! Conv.avail.potential energy (CAPE)
  TYPE(VARIABLE_2D) :: VCAPES   ! CAPE-Shear
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:) :: VMXCAP6   ! Bins for maximum CAPE in last 6 hours
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:) :: VMXCAPS6   ! Bins for maximum CAPE-Shear in last 6 hours
  TYPE(VARIABLE_2D) :: VTSRC   ! Top solar radiation clear sky
  TYPE(VARIABLE_2D) :: VTTRC   ! Top thermal radiation clear sky
  TYPE(VARIABLE_2D) :: VSSRC   ! Surface solar radiation clear sky
  TYPE(VARIABLE_2D) :: VSTRC   ! Surface thermal radiation clear sky
  TYPE(VARIABLE_2D) :: VES   ! Evaporation of snow
  TYPE(VARIABLE_2D) :: VSMLT   ! Snow melt
  TYPE(VARIABLE_2D) :: V10FG   ! Wind gust at 10 m (max since previous pp)
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:) :: V10FG6   ! Bins for wind gust at 10 m (max since last 6 hours)
  TYPE(VARIABLE_2D) :: V10FGCV   ! convective wind gust at 10m for current time level (m/s)
  TYPE(VARIABLE_2D) :: VI10FG   ! Wind gust at 10 m ("instantaneous")
  TYPE(VARIABLE_2D) :: VLSPF   ! Large scale precipitation fraction
  TYPE(VARIABLE_2D) :: VTCO3   ! Total ozone content in a vertical column
  TYPE(VARIABLE_2D) :: VVIMD   ! Vertically integrated mass divergence
  TYPE(VARIABLE_2D) :: VSPARC   ! Surface clear-sky parallel radiation
  TYPE(VARIABLE_2D) :: VSTINC   ! Top of atmosphere incident solar radiation
  TYPE(VARIABLE_2D) :: VCBASE   ! Cloud base level
  TYPE(VARIABLE_2D) :: V0DEGL   ! Zero deg. level
  TYPE(VARIABLE_2D) :: VVISIH   ! Horizontal visibility
  TYPE(VARIABLE_2D) :: VCIN   ! CIN
  TYPE(VARIABLE_2D) :: VKINDEX   ! Convective K-Index
  TYPE(VARIABLE_2D) :: VTTINDEX   ! Convective TT-Index
  TYPE(VARIABLE_2D) :: VCBASEA   ! Cloud base aviation
  TYPE(VARIABLE_2D) :: VCTOPC   ! Cloud top convective
  TYPE(VARIABLE_2D) :: VZTWETB0   ! Height of 0 deg wet bulb temperature
  TYPE(VARIABLE_2D) :: VZTWETB1   ! Height of 1 deg wet bulb temperature
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:) :: VTCGHG   ! Total column greenhouse gases
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:) :: VTCCHEM   ! Total column chemistry
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:,:) :: VAERODIAG   ! Per-aerosol-type diagnostics
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:,:) :: VAERO_WVL_DIAG   ! Per-wavelength aerosol optical diagnostics
  TYPE(VARIABLE_2D) :: V100U   ! 100m zonal wind
  TYPE(VARIABLE_2D) :: V100V   ! 100m meridional wind
  TYPE(VARIABLE_2D) :: V200U   ! 200m zonal wind
  TYPE(VARIABLE_2D) :: V200V   ! 200m meridional wind
  TYPE(VARIABLE_2D) :: VZUST   ! Friction velocity
  TYPE(VARIABLE_2D) :: V10NU   ! 10m zonal neutral wind
  TYPE(VARIABLE_2D) :: V10NV   ! 10m meridional neutral wind
  TYPE(VARIABLE_2D) :: VDNDZN   ! Minimum vertical refractivity gradient
  TYPE(VARIABLE_2D) :: VDNDZA   ! Mean vertical refractivity gradient
  TYPE(VARIABLE_2D) :: VDCTB   ! Duct base height
  TYPE(VARIABLE_2D) :: VTPLB   ! Trapping layer base height
  TYPE(VARIABLE_2D) :: VTPLT   ! Trapping layer top height
  TYPE(VARIABLE_2D) :: VODSS   ! optical depth sea salt aerosols
  TYPE(VARIABLE_2D) :: VODDU   ! optical depth dust aerosols
  TYPE(VARIABLE_2D) :: VODOM   ! optical depth organic m. aerosols
  TYPE(VARIABLE_2D) :: VODBC   ! optical depth black C aerosols
  TYPE(VARIABLE_2D) :: VODSU   ! optical depth sulphate aerosols
  TYPE(VARIABLE_2D) :: VODNI   ! optical depth nitrate aerosols
  TYPE(VARIABLE_2D) :: VODAM   ! optical depth ammonium aerosols
  TYPE(VARIABLE_2D) :: VODSOA   ! optical depth secondary organic aerosols
  TYPE(VARIABLE_2D) :: VODVFA   ! optical depth volcanic flying ash
  TYPE(VARIABLE_2D) :: VODVSU   ! optical depth volcanic sulphate aerosols
  TYPE(VARIABLE_2D) :: VODTOACC   ! optical depth total aerosol accumulated
  TYPE(VARIABLE_2D) :: VAEPM1   ! particulate matter le 1 um
  TYPE(VARIABLE_2D) :: VAEPM25   ! particulate matter le 2.5um
  TYPE(VARIABLE_2D) :: VAEPM10   ! particulate matter le 10 um
  TYPE(VARIABLE_2D) :: VUVBED   ! UV biologically effective dose
  TYPE(VARIABLE_2D) :: VUVBEDCS   ! UV biologically effective dose clear sky
  TYPE(VARIABLE_2D) :: VLITOTI   ! instantaneous total lightning flash density
  TYPE(VARIABLE_2D) :: VLICGI   ! instantaneous cloud-to-ground lightning flash density
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:) :: VLITOTA6   ! Bins for averaged total lightning over last 6 hours
  TYPE(VARIABLE_2D), ALLOCATABLE, DIMENSION(:) :: VLICGA6   ! Bins for averaged cloud-to-ground lightning over last 6 hours

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VDIAG

TYPE SURFACE_VARIABLE_GROUP_SATSIM
  ! Diagnostic surface variable group
  TYPE(VARIABLE_3D) :: VCLBT   ! Cloudy brightness temperature
  TYPE(VARIABLE_3D) :: VCSBT   ! Clear-sky brightness temperature

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_4D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_SATSIM

TYPE SURFACE_VARIABLE_GROUP_WAVES
  ! Diagnostic surface variable group
  TYPE(VARIABLE_2D) :: VCHAR   ! Charnock parameter as modified by the wave model.
  TYPE(VARIABLE_2D) :: VUSTOKES   ! U-component of the surface Stokes drift.
  TYPE(VARIABLE_2D) :: VVSTOKES   ! V-component of the surface Stokes drift.
  TYPE(VARIABLE_2D) :: VPHIOC   ! Energy flux to ocean.
  TYPE(VARIABLE_2D) :: VPHIAW   ! Energy flux to ocean waves.
  TYPE(VARIABLE_2D) :: VTAUOC   ! Momentum flux to ocean.
  TYPE(VARIABLE_2D) :: VEMEAN   ! Wave variance.
  TYPE(VARIABLE_2D) :: VFMEAN   ! Wave mean frequency.

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_WAVES

TYPE SURFACE_VARIABLE_GROUP_WAM
  ! Diagnostic surface variable group
  TYPE(VARIABLE_2D) :: VU10N   ! 10m neutral wind U-component passed to the wave model (WAM).
  TYPE(VARIABLE_2D) :: VV10N   ! 10m neutral wind V-component passed to the wave model (WAM).
  TYPE(VARIABLE_2D) :: VRHO   ! surface density passed to the wave model (WAM).
  TYPE(VARIABLE_2D) :: VZIL   ! ZI/L passed to the wave model (used for gustiness in WAM).
  TYPE(VARIABLE_2D) :: VCIF   ! Sea ice fraction passed to the wave model (WAM).
  TYPE(VARIABLE_2D) :: VCLK   ! Lake cover passed to the wave model (WAM).
  TYPE(VARIABLE_2D) :: VUCURW   ! Ocean current    U-component passed to the wave model (WAM).
  TYPE(VARIABLE_2D) :: VVCURW   ! Ocean current    V-component passed to the wave model (WAM).

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_WAM

TYPE SURFACE_VARIABLE_GROUP_VEXTRA
  ! Diagnostic surface variable group

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_4D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VEXTRA

TYPE SURFACE_VARIABLE_GROUP_VEXTRDI
  ! Diagnostic surface variable group
  TYPE(VARIABLE_3D) :: VXEDR   ! Eddy diffusivity rate

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_4D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VEXTRDI

TYPE SURFACE_VARIABLE_GROUP_VPRECIP
  ! Diagnostic surface variable group
  TYPE(VARIABLE_3D) :: VPRECIP   ! Diagnostic of precipitations type

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_4D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VPRECIP

TYPE SURFACE_VARIABLE_GROUP_VPRECIP2
  ! Diagnostic surface variable group
  TYPE(VARIABLE_3D) :: VPRECIP2   ! Diagnostic of precipitations type

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_4D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VPRECIP2

TYPE SURFACE_VARIABLE_GROUP_VEXTR2
  ! Diagnostic surface variable group

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VEXTR2

TYPE SURFACE_VARIABLE_GROUP_SFORC
  ! Diagnostic surface variable group

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_SFORC

TYPE SURFACE_VARIABLE_GROUP_SFLUX
  ! Diagnostic surface variable group

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_SFLUX

TYPE SURFACE_VARIABLE_GROUP_VO3ABC
  ! Diagnostic surface variable group
  TYPE(VARIABLE_2D) :: VA   ! A climatological ozone profile
  TYPE(VARIABLE_2D) :: VB   ! B climatological ozone profile
  TYPE(VARIABLE_2D) :: VC   ! C climatological ozone profile

  ! Store a field encapsualting the entire group array for backward compatibility
  TYPE(FIELD_3D) :: F_GROUP
END TYPE SURFACE_VARIABLE_GROUP_VO3ABC


TYPE SURFACE_VARIABLES
  ! Global variable and field storage for surface variables

  ! Prognostic variable groups
  TYPE(SURFACE_VARIABLE_GROUP_SOILB) :: GSP_SB
  TYPE(SURFACE_VARIABLE_GROUP_SNOWG) :: GSP_SG
  TYPE(SURFACE_VARIABLE_GROUP_LAKEB) :: GSP_SL
  TYPE(SURFACE_VARIABLE_GROUP_RESVR) :: GSP_RR
  TYPE(SURFACE_VARIABLE_GROUP_CLS) :: GSP_CL
  TYPE(SURFACE_VARIABLE_GROUP_OML) :: GSP_OM
  TYPE(SURFACE_VARIABLE_GROUP_EXTRP) :: GSP_EP
  TYPE(SURFACE_VARIABLE_GROUP_XTRP2) :: GSP_X2
  TYPE(SURFACE_VARIABLE_GROUP_CANRI) :: GSP_CI

  ! Diagnostic variable groups
  TYPE(SURFACE_VARIABLE_GROUP_VARSF) :: GSD_VF
  TYPE(SURFACE_VARIABLE_GROUP_VCLIH) :: GSD_VH
  TYPE(SURFACE_VARIABLE_GROUP_VCLIK) :: GSD_VK
  TYPE(SURFACE_VARIABLE_GROUP_VCLIP) :: GSD_VP
  TYPE(SURFACE_VARIABLE_GROUP_VCLIV) :: GSD_VV
  TYPE(SURFACE_VARIABLE_GROUP_VCLIA) :: GSD_VA
  TYPE(SURFACE_VARIABLE_GROUP_VCLIN) :: GSD_VN
  TYPE(SURFACE_VARIABLE_GROUP_VDIAGO2) :: GSD_V2
  TYPE(SURFACE_VARIABLE_GROUP_VDIAGO3) :: GSD_V3
  TYPE(SURFACE_VARIABLE_GROUP_VDIAG) :: GSD_VD
  TYPE(SURFACE_VARIABLE_GROUP_SATSIM) :: GSD_SM
  TYPE(SURFACE_VARIABLE_GROUP_WAVES) :: GSD_WS
  TYPE(SURFACE_VARIABLE_GROUP_WAM) :: GSD_WW
  TYPE(SURFACE_VARIABLE_GROUP_VEXTRA) :: GSD_XA
  TYPE(SURFACE_VARIABLE_GROUP_VEXTRDI) :: GSD_DI
  TYPE(SURFACE_VARIABLE_GROUP_VPRECIP) :: GSD_XP
  TYPE(SURFACE_VARIABLE_GROUP_VPRECIP2) :: GSD_XP2
  TYPE(SURFACE_VARIABLE_GROUP_VEXTR2) :: GSD_X2
  TYPE(SURFACE_VARIABLE_GROUP_SFORC) :: GSD_SFO
  TYPE(SURFACE_VARIABLE_GROUP_SFLUX) :: GSD_SFL
  TYPE(SURFACE_VARIABLE_GROUP_VO3ABC) :: GSD_VC
END TYPE SURFACE_VARIABLES

END MODULE SURFACE_VARIABLES_MOD
