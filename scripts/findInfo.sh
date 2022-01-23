#!/bin/bash

rm -f h.pl h.xml

./scripts/findVarsInfo.pl src/local/arpifs/module/field_variables_mod.F90
./scripts/findFieldInfo.pl \
  src/local/arpifs/module/cpg_type_mod.F90  src/local/arpifs/module/mf_phys_type_mod.F90  \
  src/local/arpifs/module/mf_phys_base_state_type_mod.F90 src/local/arpifs/module/mf_phys_next_state_type_mod.F90 \
  src/local/arpifs/module/mf_phys_surface_type_mod.F90 src/local/arpifs/module/surface_variables_mod.F90 \
  src/local/arpifs/module/surface_views_module.F90 
