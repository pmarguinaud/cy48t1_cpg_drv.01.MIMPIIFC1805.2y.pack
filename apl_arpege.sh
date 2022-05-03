#!/bin/bash

set -x

for f in \
  src/local/arpifs/phys_dmn/apl_arpege.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_init.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_init_surfex.F90 \
  src/local/arpifs/phys_dmn/apl_arpege_oceanic_fluxes.F90
do
./scripts/apl_arpege.pl $f
done
