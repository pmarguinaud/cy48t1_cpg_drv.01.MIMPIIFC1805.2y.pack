#!/bin/bash

for f in acbl89.F90 acevolet.F90 aclender.F90 actke.F90 acturb.F90 fl2hl.F90 hl2fl.F90
do
  ./scripts/ppsd.pl src/local/arpifs/phys_dmn/$f
done


