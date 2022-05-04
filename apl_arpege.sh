#!/bin/bash

set -x

for f in \
  src/local/arpifs/phys_dmn/apl_arpege_init.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_init_surfex.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_oceanic_fluxes.F90 \
  src/local/arpifs/phys_dmn/apl_wind_gust.F90 \
  src/local/arpifs/phys_dmn/mf_phys_mocon.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_shallow_convection_and_turbulence.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_albedo_computation.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_aerosols_for_radiation.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_cloudiness.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_radiation.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_soil_hydro.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_deep_convection.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_surface.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_precipitation.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_hydro_budget.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_dprecips.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_atmosphere_update.F90 \
  src/local/arpifs/adiab/cputqy_aplpar_expl.F90 \
  src/local/arpifs/adiab/acctnd0.F90 \
  src/local/arpifs/adiab/cputqy0.F90 \
  src/local/arpifs/phys_dmn/mf_phys_transfer.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_surface_update.F90 \
  src/local/arpifs/phys_dmn/apl_arpege.F90 
do
./scripts/apl_arpege.pl $f 
done

grep _parallel src/local/arpifs/phys_dmn/apl_arpege_parallel.F90
